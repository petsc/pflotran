
Module pflow_step

#include "include/finclude/petscsnes.h"
  use petscsnes
  use pflow_gridtype_module
  
  implicit none
  
  private
#include "definitions.h"

  public pflowgrid_update_dt
  public pflowGrid_step
  public pflowGrid_update



Contains
!#include "pflowgrid_update_dt.F90"

  subroutine pflowgrid_update_dt(grid, timestep,its)
  implicit none

  type(pflowGrid), intent(inout) :: grid
   type(time_stepping_context), intent(inout) :: timestep
  
  integer, intent(in) :: its

  real*8 :: fac,dtt,up,utmp,uc,ut,uus
  
  if (timestep%iaccel == 0) return
  
    fac = 0.5d0
    if (its >= timestep%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = timestep%dpmxe/(grid%dpmax+0.1)
      utmp = 0.D0 !grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = 0.D0 !grid%dcmxe/(grid%dcmax+1.d-6)
      uus= timestep%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,uus)
    endif

    dtt = fac * grid%dt * (1.d0 + ut)
  

  
  if (dtt > 2.d0 * grid%dt) dtt = 2.d0 * grid%dt 
  if (dtt > timestep%dt_max) dtt = timestep%dt_max
  if(dtt>.25d0*grid%t .and. grid%t>1.d-2) dtt=.25d0*grid%t
  grid%dt = dtt

  end subroutine pflowgrid_update_dt

!======================================================================

!#include "pflowgrid_step.F90"

  subroutine pflowGrid_step(grid,pflowsolv,timestep ,ntstep,kplt,iplot,iflgcut,ihalcnt,its)
  
  use IMS_module, only : pflow_ims_step_maxchange
  use pflow_output_module
  use ims_module
  use pflow_solv_module
  
  implicit none

  type(pflowGrid), intent(inout) :: grid
  type(pflow_solver_context), intent(inout) :: pflowsolv
   type(time_stepping_context), intent(inout) :: timestep  

  integer :: ierr, ihalcnt, ns !i,m,n,ix,jy,kz,j
  integer :: its,kplt,iplot,ntstep !,idpmax,idtmpmax,idcmax
  integer :: icut, iflgcut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  integer update_reason
  real*8 :: tsrc
  type(pflow_localpatch_info), pointer :: locpat
    
  locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr
    
! real*8, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  its = 0
  icut = 0

!  print *,'Stepping...'
  ! Perform some global-to-local scatters to update the ghosted vectors.
  ! We have to do this so that the routines for calculating the residual
  ! and the Jacobian will have the ghost points they need.
  ! Note that we don't do the global-to-local scatter for the ppressure 
  ! vector, as that needs to be done within the residual calculation routine
  ! because that routine may get called several times during one Newton step
  ! if a method such as line search is being used.
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                            grid%porosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                          grid%porosity_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                            grid%tor_loc, ierr)

  call DAGlobalToLocalEnd(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                          grid%tor_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                            grid%icap_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                          grid%icap_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%ithrm, INSERT_VALUES, &
                            grid%ithrm_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%ithrm, INSERT_VALUES, &
                          grid%ithrm_loc, ierr)

  !print *,'Step :: End trans distribute'
  grid%t = grid%t + grid%dt
  grid%flowsteps = grid%flowsteps + 1

! Adjust time step to plot time
  if (grid%t + 0.2*grid%dt >= timestep%tplot(kplt)) then
    grid%t = grid%t - grid%dt
    grid%dt = timestep%tplot(kplt) - grid%t
    if (grid%dt > timestep%dt_max) then
      grid%dt = timestep%dt_max
      grid%t = grid%t + grid%dt
    else
      grid%t = timestep%tplot(kplt)
      iplot = 1
    endif
  else if (grid%flowsteps == grid%stepmax) then
    iplot = 1
  endif

! source/sink time step control
  if (grid%nblksrc > 0) then
    ns = 1
    tsrc = grid%timesrc(grid%isrc1,ns)
    if (grid%t >= tsrc ) then
      if(grid%t > tsrc +1D2) then
        grid%t = grid%t - grid%dt
        grid%dt = tsrc - grid%t
        grid%t = tsrc
      endif
      grid%isrc1 = grid%isrc1 + 1
    endif
  endif
  
  if (iflgcut == 1) then
    ihalcnt = ihalcnt + 1
    if(ihalcnt > grid%ndtcmx) then
      iflgcut = 0
      ihalcnt = 0
    endif
  endif
  
  !set maximum time step
  if (ntstep > timestep%nstpmax) then
  ! do nothing - keep current dt_max
  else if (grid%t > timestep%tstep(ntstep)) then
    ntstep = ntstep + 1
    if (ntstep <= timestep%nstpmax) then
      timestep%dt_max = timestep%dtstep(ntstep)
    else
      timestep%dt_max = timestep%dtstep(timestep%nstpmax)
    endif
  endif
  
  do
! call VecView(grid%r,PETSC_VIEWER_STDOUT_WORLD,ierr) 
! call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr) 
 !print *, 'Begin solving'
	grid%iphch=0
      if(grid%use_ksp == PETSC_TRUE)then
         call pflow_solve(grid,pflowsolv, its,snes_reason,ierr)
 	   else 
         call SNESSolve(pflowsolv%snes, PETSC_NULL, grid%xx, ierr)
      endif

   !print *,'pflow_step, finish SNESSolve'
   call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   if(grid%use_ksp /= PETSC_TRUE) call SNESGetIterationNumber(pflowsolv%snes, its, ierr)
!   call KSPGetIterationNumber(grid%ksp, kspits, ierr)
  
!   printout residual and jacobian to screen for testing purposes
!   call VecView(grid%r,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   call MatView(grid%J,PETSC_VIEWER_STDOUT_WORLD,ierr)

!   call SNESDefaultMonitor(grid%snes, its, r2norm, ierr)

    if(grid%use_ksp /= PETSC_TRUE)call SNESGetConvergedReason(pflowsolv%snes, snes_reason, ierr)

!   parameter (SNES_CONVERGED_ITERATING         =  0)
!   parameter (SNES_CONVERGED_FNORM_ABS         =  2)
!   parameter (SNES_CONVERGED_FNORM_RELATIVE    =  3)
!   parameter (SNES_CONVERGED_PNORM_RELATIVE    =  4)
!   parameter (SNES_CONVERGED_GNORM_ABS         =  5)
!   parameter (SNES_CONVERGED_TR_REDUCTION      =  6)
!   parameter (SNES_CONVERGED_TR_DELTA          =  7)

!   parameter (SNES_DIVERGED_FUNCTION_COUNT     = -2)
!   parameter (SNES_DIVERGED_FNORM_NAN          = -4)
!   parameter (SNES_DIVERGED_MAX_IT             = -5)
!   parameter (SNES_DIVERGED_LS_FAILURE         = -6)
!   parameter (SNES_DIVERGED_TR_REDUCTION       = -7)
!   parameter (SNES_DIVERGED_LOCAL_MIN          = -8)


!******************************************************************
! since calculation on saturation has special requirments, need special treatment
! update_reason
    !call PETScBarrier(PETSC_NULL_OBJECT, ierr)
    update_reason = 1

    call IMS_Update_Reason(update_reason, grid, locpat)
      if (grid%myrank==0) print *,'update_reason: ',update_reason, its
    

!******************************************************************
    
    if(snes_reason < 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      iflgcut = 1

      if(icut > grid%icut_max .or. grid%dt<1d-9 ) then
!       call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)
        if(grid%myrank == 0) then
!         t = pflowgrid_get_t(grid)
          print *,"icut_max exceeded. icut= ",icut, "t= ",grid%t/grid%tconv,grid%dt, &
            ". Stopping execution."
          print *, "Final state:"
        endif
        iplot = 1
 !       call pflow_output(grid%ppressure,grid%ttemp,grid%conc,grid%phis,grid%porosity, &
 !       grid%perm_xx, grid%perm_yy, grid%perm_zz, &
 !       grid%porosity, grid%sat, grid%vl, &
 !       grid%c_nat,grid%vl_nat,grid%p_nat,grid%t_nat,grid%s_nat,grid%phis_nat,grid%por_nat, &
 !       grid%ibrkface, grid%jh2o, grid%nphase, grid%nmax, &
 !       grid%snes, &
 !       grid%t, grid%dt, grid%tconv, grid%flowsteps, grid%rk, &
 !       grid%k1brk,grid%k2brk,grid%j1brk,grid%j2brk,grid%i1brk,grid%i2brk, &
 !       grid%nx,grid%ny,grid%nz,grid%nxy,grid%dx0,grid%dy0,grid%dz0,grid%x,grid%y,grid%z, &
 !       grid%da_nphase_dof,grid%da_1_dof,grid%da_3np_dof,grid%da_ndof,grid%ndof, &
 !       kplt,iplot,grid%iprint,grid%ibrkcrv, &
 !       grid%itecplot,grid%write_init,grid%myrank)
		
	!	    call pflow_output(grid,kplt,iplot)
        ! The above line won't work when scope is restricted properly!
        ! Replace this with a different function call!
     !   call pflowgrid_destroy(grid)
        call PetscFinalize()
        stop
      endif

      grid%t = grid%t - grid%dt
      grid%dt = 0.5d0 * grid%dt
      grid%t = grid%t + grid%dt
    
      if(grid%myrank == 0) write(*,'('' -> Cut time step: snes='',i3, &
  &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '',1pe12.4,i2)') &
      snes_reason,icut,grid%icutcum,grid%t/grid%tconv,grid%dt/grid%tconv,iflgcut

      call pflow_ims_timecut(grid, locpat)
		
	
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  grid%newtcum = grid%newtcum + its
  grid%icutcum = grid%icutcum + icut

! print screen output
  if (grid%myrank == 0) then
    if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
      write(*, '(/," FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
 &    " snes: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4,"]")') &
      grid%flowsteps,grid%t/grid%tconv,grid%dt/grid%tconv,timestep%tunit, &
      snes_reason,its,grid%newtcum,icut,grid%icutcum
    
      write(IUNIT2, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
 &    " snes: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4,"]")') &
      grid%flowsteps,grid%t/grid%tconv,grid%dt/grid%tconv,timestep%tunit, &
      snes_reason,its,grid%newtcum,icut,grid%icutcum
    endif
  endif
  
  ! calculate maxium changes in fields over a time step

     call pflow_ims_step_maxchange(grid)
	 
			
	 if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
 &      " dsmx= ",1pe12.4)') &
        grid%dpmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
 &      " dsmx= ",1pe12.4)') &
        grid%dpmax,grid%dsmax
      endif
    endif
   
 
  if (grid%myrank == 0 .and. mod(grid%flowsteps,grid%imod) == 0) then
    print *, ""
  endif

  end subroutine pflowGrid_step

!==========================================================================
  
subroutine pflowGrid_update (grid)
  
  use ims_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), pointer :: locpat
  
  locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr
  

  ! update solution vector and physical properties (VecCopy(x,y): y=x)
  
  
  ! call VecView(grid%ppressure,PETSC_VIEWER_STDOUT_WORLD,ierr)
  call pflow_update_ims(grid, locpat)
  
  
end subroutine pflowGrid_update

!======================================================================

end module pflow_step
