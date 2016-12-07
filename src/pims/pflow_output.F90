
module pflow_output_module

  public
  
  contains

 !======================================================================

  subroutine porperm_out(t, dt, tconv, kplt, nx, ny, nz, nmax, &
                         x, y, z, flowsteps, scat_1dof, da_1_dof, &
                         porosity, por_nat, perm, perm_nat, myrank)
#include "include/finclude/petscsnes.h"
  use petscsnes
  implicit none

#include "definitions.h"
  
  integer, intent(in) :: kplt, nx, ny, nz, nmax, flowsteps, myrank
  real*8, intent(in) :: t, dt, tconv, x(:), y(:), z(:)
  Vec, intent(inout) :: porosity, por_nat, perm, perm_nat
  Vec :: por_all, perm_all
  
  DA :: da_1_dof
  
  VecScatter :: scat_1dof
  
  integer :: ierr, n
  
  real*8 :: fac, tyr
  
  real*8, pointer :: por_p(:), perm_p(:)
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab

#include "definitions.h"

  tab = char(9)
  q = '","'
  tyr =t/tconv

! write out porosity, permeability

! porosity field
  call DAGlobalToNaturalBegin(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)
print *, 'Scatter:: por'
! permeability field
! call DAGlobalToNaturalBegin(da_3, perm, INSERT_VALUES, perm_nat, ierr)
! call DAGlobalToNaturalEnd(da_3, perm, INSERT_VALUES, perm_nat, ierr)
  call DAGlobalToNaturalBegin(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)

! call VecScatterCreate(scatter, ierr)
#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(por_nat, por_all, ierr); CHKERRQ(ierr)
! call VecScatterCreateToZero(por_nat, scatter, por_all, ierr)
  call VecScatterBegin(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
  call VecScatterEnd(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
! call VecConvertMPIToMPIZero(perm_nat, perm_all, ierr); CHKERRQ(ierr)
  call VecScatterBegin(perm_nat, perm_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
  call VecScatterEnd(perm_nat, perm_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
#else
! call VecConvertMPIToSeqAll(por_nat, por_all, ierr); CHKERRQ(ierr)
! call VecScatterCreateToAll(por_nat, scatter, por_all, ierr)

! call VecScatterBegin(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
!      scat_1dof, ierr)
! call VecScatterEnd(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
!      scat_1dof, ierr)

! call VecConvertMPIToSeqAll(perm_nat, perm_all, ierr); CHKERRQ(ierr)

  call VecScatterCreateToAll(por_nat, scat_1dof, por_all, ierr)
  call VecScatterBegin(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
  call VecScatterEnd(por_nat, por_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
! call VecScatterDestroy(scat_1dof, ierr)

  call VecScatterCreateToAll(perm_nat, scat_1dof, perm_all, ierr)
  call VecScatterBegin(perm_nat, perm_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
  call VecScatterEnd(perm_nat, perm_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_1dof, ierr)
  call VecScatterDestroy(scat_1dof, ierr)
#endif

  if(myrank == 0) then
    
    if (kplt < 10) then
      write(fname,'(a10,i1,a4)') 'pflow_por0', kplt, '.dat'
    else
      write(fname,'(a9,i2,a4)') 'pflow_por', kplt, '.dat'
    endif
    
    write(*,*) '--> write output file: ',fname
    
    open(unit=IUNIT3,file=fname,action="write")
      
    call VecGetArrayF90(por_all, por_p, ierr)
    call VecGetArrayF90(perm_all, perm_p, ierr)
      
    fac = 1.d12
      
    if (nx+ny+nz .gt. nx*ny*nz) then ! 1D
        
      write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," y", &
 &    " dt =",1pe12.4," sec")') &
      flowsteps,t,tyr,dt
      write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
 &    "  -por- ","  -perm- ")')
!&    "  -por-      -permx-      -permy-      -permz-&
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (nz == 1) then ! 2D x-y

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,ny
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (ny == 1) then ! 2D x-z

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else ! 3D

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,nx,ny,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
    endif
    call VecRestoreArrayF90(por_all, por_p, ierr)
    call VecRestoreArrayF90(perm_all, perm_p, ierr)
    close (IUNIT3)
  endif
  call VecDestroy(por_all, ierr)
  call VecDestroy(perm_all, ierr)

  end subroutine porperm_out




 subroutine pflow_var_output(grid,timestep,kplt,iplot)
  
#include "include/finclude/petscda.h"
  use petscda
  use pflow_gridtype_module
  implicit none

   type(time_stepping_context), intent(inout) :: timestep

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  
  integer, intent(inout) :: iplot, kplt
   
  Vec :: vl_nat, vl_all 
  real*8, pointer :: p_p(:), t_p(:), xx_p(:), phis_p(:), vl_p(:), var_p(:)
  real*8, pointer :: x_p(:), iphase_p(:)
  real*8, pointer :: fldflx(:), fldvol(:)
  real*8 xxx(1:grid%ndof)
  real*8, allocatable:: vvar(:), vavel(:, :)
  integer :: nxm1,nym1,nzm1, ndex, na, nv,ng
  integer :: ierr,mm, nvar_out, status(MPI_STATUS_SIZE)
  integer :: i,ip,j,jn,k,ibrk,n,nc,nn, ipha, jj, ifound
  
  real*8 :: tyr, sum1, sum2, sum1v, sum2v, area, vol, vel, vf, iipha
  integer, save :: icall, icall_brk, icall_plt
   VecScatter :: scat_3npdof
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab
  
  character*6, allocatable :: var_name(:)!2+grid%nphase*(2+grid%nspec))

  type(pflow_localpatch_info), pointer :: locpat
  
  data icall/1/, icall_brk/1/, icall_plt/1/
  
  nvar_out= grid%nphase*2+1
  allocate(var_name(nvar_out))
  
  var_name(1) = 'T'
  var_name(2) = 'p'
  do i = 1,grid%nphase
    write(var_name(i+2),'(a2,i3,a1)')'p(',i,')'
  !  print *, var_name(i+2)
   enddo
  tab = char(9)
  q = '","'
  tyr = grid%t/grid%tconv

  allocate(vvar(1:grid%size_var_use))
  
  !save icall, icall_brk, icall_plt
  
! ibrkcrv = 0, no time-history
!         = 1, time-history output
  
! iplot  = 0, only time-history
!        = 1, spatial plot

! iprint = -1, none
!        =  0, p, t, s, c
!        =  1, vl
!        =  2, porosity, (permeability)
  !return
  
  if ((grid%ibrkcrv == 0 .and. iplot == 0) .or. grid%iprint == -1) then
    if (grid%iprint==-1 .and. iplot==1) then
      kplt = kplt + 1
      iplot = 0
    endif
    return
 endif
  
     
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1 .and. icall == 1) then
    icall = 0
    if (grid%myrank == 0) then
      write(fname,'(a9,a4)') 'pflow_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' iplot= ',iplot
      open(unit=IUNIT4,file=fname,action="write")
      write(IUNIT4,'("%#t            dt",100i12)') (i,i=1,grid%ibrkcrv), &
      (i,i=1,grid%ibrkcrv)
    endif
    if (iplot == 0) return
 endif
  
  if (grid%ndof == 1 .and. iplot == 0) return ! no time-history plot for 1 dof

! Create Natural Vec for output: use VecDuplicate here?
  
  ! time-history plots



! plot spatial data (iplot > 0)
  
    call VecGetArrayF90(grid%xx,  xx_p, ierr)
    locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr

    var_p => locpat%var
	
	if (kplt < 10) then
      write(fname,'(a7,i1,a4)') 'pflow00', kplt, '.dat'
    else if (kplt < 100) then
      write(fname,'(a6,i2,a4)') 'pflow0', kplt, '.dat'
    else
      write(fname,'(a5,i3,a4)') 'pflow', kplt, '.dat'
    endif
    if(grid%myrank==0) write(*,*) '--> write output file: ',fname
    
    if(grid%myrank==0) open(unit=IUNIT3,file=fname,action="write")
      
 
   if(grid%myrank==0)then
    write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,timestep%tunit
!    write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))') "x",q,"y",q,"z",q, var_name(1), ((q,var_name(i)),i=2,nvar_out),'"'
    write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, '' , J='',i4,'' K='',i4)') tyr,grid%nx,grid%ny,grid%nz
    endif 
	
	   		 
      	ndex=1
				do na=0, grid%nmax-1
				      ifound=0
					  if(grid%myrank==0) then
				       ifound=0
					   do n=ndex,locpat%nlmax
					   	 if(locpat%nL2A(n) == na) then
					      ng=locpat%nL2G(n)
						  ifound=1
						  jn= 1 + (n-1)*grid%ndof 
					      nv = 1 + (ng-1) * grid%size_var_node
					     
					      xxx= xx_p(jn:jn+grid%ndof-1)
						  vvar= var_p(nv:nv + grid%size_var_use-1)
					      ndex=n
						  exit	 
					    endif
					   enddo
					   if(ifound==0)then  
					   
					   	call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, na, PETSC_COMM_WORLD, status,ierr) 					
					    call MPI_Recv(vvar, grid%size_var_use ,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr) 
                       endif
					    
					  else
					    do n=ndex,locpat%nlmax
					      if(locpat%nL2A(n) == na) then
					        jn= 1 + (n-1)*grid%ndof
							 ng=locpat%nL2G(n)
					        nv = 1 + (ng-1) * grid%size_var_node
		        
						   xxx= xx_p(jn:jn+grid%ndof-1)
						   call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
				           call MPI_Send(var_p(nv:nv + grid%size_var_use-1), grid%size_var_use ,MPI_DOUBLE_PRECISION, 0,&
                                na+na, PETSC_COMM_WORLD, ierr)   
					       ndex=n
						   exit
					      endif
					   enddo  
				     endif     	     
                
				  call MPI_Barrier(PETSC_COMM_WORLD, ierr)
				
				  if(grid%myrank ==0 )then
                       write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), &
                          vvar(1:2+grid%nphase)! , & ! Saturations
                          !vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
                          !vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
                  
                   endif
  				enddo
			   	
	
		
  !   write out initial conditions

    if (grid%write_init == 1) then
				
		 if(grid%myrank ==0)then
			  fname = 'pflow_init0.dat'
			  write(*,*) '--> write output file: ',fname
			  close (IUNIT3)
			  open(unit=IUNIT3,file=fname,action="write")
!			  write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q, var_name(1), ((q,var_name(i)),i=2,nvar_out),'"'
		 endif
		
			   		 
  	ndex=1
				do na=0, grid%nmax-1
				      ifound=0
					  if(grid%myrank==0) then
				       ifound=0
					   do n=ndex,locpat%nlmax
					   	 if(locpat%nL2A(n) == na) then
					      ifound=1
						  jn= 1 + (n-1)*grid%ndof 
						  ng=locpat%nL2G(n)
					      nv = 1 + (ng-1) * grid%size_var_node
					      xxx= xx_p(jn:jn+grid%ndof-1)
						  vvar= var_p(nv:nv + grid%size_var_use-1)
					      ndex=n
						  exit	 
					    endif
					   enddo
					   if(ifound==0)then  
					   	call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, na, PETSC_COMM_WORLD, status,ierr) 					
					    call MPI_Recv(vvar, grid%size_var_use ,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr) 
                       endif
					    
					  else
					    do n=ndex,locpat%nlmax
					      if(locpat%nL2A(n) == na) then
					        jn= 1 + (n-1)*grid%ndof
							ng=locpat%nL2G(n) 
					        nv = 1 + (ng-1) * grid%size_var_node
						    xxx= xx_p(jn:jn+grid%ndof-1)
						   call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
				           call MPI_Send(var_p(nv:nv + grid%size_var_use-1), grid%size_var_use ,MPI_DOUBLE_PRECISION,&
                            0, na+na, PETSC_COMM_WORLD, ierr)   
					       ndex=n
						   exit
					      endif
					   enddo  
				     endif     	     
                
				  call MPI_Barrier(PETSC_COMM_WORLD, ierr)
				  
				  if(grid%myrank ==0 )then
                        write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), vvar(1:2+grid%nphase)
                          
                  
                    endif
   				enddo
			  
	endif
    
      
    
 
  ! if (iprint == 0) call VecScatterDestroy(scat_1dof, ierr)
 ! call VecScatterDestroy(scat_nph, ierr)

  close (IUNIT3)

    call VecRestoreArrayF90(grid%xx,  xx_p, ierr)
  !  call VecRestoreArrayF90(grid%var, var_p, ierr)
    nullify(var_p)
   
  if (grid%iprint >= 1) then

!   if (ibrkcrv == 0) then

!   velocity fields

 call DACreateNaturalVector(grid%da_3np_dof, vl_nat, ierr)
	  
		   

    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
    vl_nat,ierr)
    call DAGlobalToNaturalEnd  (grid%da_3np_dof,grid%vl,INSERT_VALUES, &
    vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must 
! be followed by a call to VecDestroy(x_all)

 







#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(vl_nat, vl_all, ierr)
  call VecScatterCreateToZero(vl_nat, scat_3npdof, vl_all, ierr)
  call VecScatterBegin(vl_nat, vl_all, INSERT_VALUES, &
       SCATTER_FORWARD, scat_3npdof, ierr)
  call VecScatterEnd(vl_nat, vl_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_3npdof, ierr)
#else
! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr)
  call VecScatterCreateToAll(vl_nat, scat_3npdof, vl_all, ierr)
  call VecScatterBegin(vl_nat, vl_all, INSERT_VALUES, &
       SCATTER_FORWARD, scat_3npdof, ierr)
  call VecScatterEnd(vl_nat, vl_all, INSERT_VALUES, SCATTER_FORWARD, &
       scat_3npdof, ierr)
 call VecScatterDestroy(scat_3npdof, ierr)
#endif




  ! endif
    
    if(grid%myrank == 0 .and. kplt > 0) then
       
      call VecGetArrayF90(vl_all, vl_p, ierr)
      allocate(vavel(3, grid%nphase))
	  vavel=0.D0       
      ! write x,y,z-velocities
      ! note ordering of velocities: v(p,d,n) = v(p+(d-1)*Np+(n-1)*3*Np)
      if (kplt < 10) then
        write(fname,'(a10,i1,a4)') 'pflow_vel0', kplt, '.dat'
      else
        write(fname,'(a9,i2,a4)') 'pflow_vel', kplt, '.dat'
      endif
      write(*,*) '--> write output file: ',fname
      open(unit=IUNIT3,file=fname,action="write")
      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,timestep%tunit
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlx',q,'vgx',q,'vly',q,'vgy',q,'vlz',q, &
          'vgz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
                       !x-connections
            nc = n
            nxm1 = nc-1
            ip = 0
            do ipha=1, grid%nphase
			  if (i>1 .and. i<grid%nx) then
               vavel(ip+1, ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
               + vl_p(ipha+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              else if (i==1) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              else if (i==grid%nx) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              endif
            enddo 
			
            !y-connections
            nc = n
            nym1 = nc-grid%nx
            ip = 1
            do ipha=1, grid%nphase
			if (j>1 .and. j<grid%ny) then
              vavel(ip+1,ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(ipha+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            else if (j==1) then
              vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (j==grid%ny) then
              vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            endif
            enddo
			
            !z-connections
            nc = n
            nzm1 = nc-grid%nxy
            ip = 2
			do ipha=1, grid%nphase
            if (k>1 .and. k<grid%nz) then
              vavel(ip+1,ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(ipha+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            else if (k==1) then
              vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (k==grid%nz) then
              vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            endif
            enddo
			
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
            ((vavel(ip,jj),jj=1,grid%nphase),ip=1,3)
            
!           write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
!           (vl_p(1+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
!           vl_p(2+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,ip=0,2)
          enddo
        enddo
      enddo
      close(IUNIT3)
     call VecRestoreArrayF90(vl_all, vl_p, ierr)
    endif
   ! call VecRestoreArrayF90(vl_all, vl_p, ierr)
	 call VecDestroy(vl_all, ierr)
     call VecDestroy(vl_nat, ierr)
  endif
  
! if (ibrkcrv >= 0 .or. iprint >= 1) then
!   call VecDestroy(vl_all, ierr)
 !   call VecScatterDestroy(scat_3dof, ierr)
!   call VecDestroy(vl_nat, ierr)
! endif
 
  ! if (grid%iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)

  kplt = kplt + 1
  iplot = 0

! call VecDestroy(grid%c_nat, ierr)
! call VecDestroy(grid%phis_nat, ierr)
! call VecDestroy(grid%t_nat, ierr)
! call VecDestroy(grid%por_nat, ierr)
! call VecDestroy(grid%p_nat, ierr)
! call VecDestroy(grid%s_nat, ierr)
! call VecDestroy(grid%vl_nat, ierr)
! call VecDestroy(grid%x_nat, ierr)

  
  
 ! if (iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)
!call VecScatterView(scat_1dof,PETSC_VIEWER_STDOUT_SELF)
  end subroutine pflow_var_output






  end module pflow_output_module
  
