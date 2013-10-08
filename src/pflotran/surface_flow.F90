#ifdef SURFACE_FLOW

module Surface_Flow_module

  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
!#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceFlowSetup, &
         SurfaceFlowKinematic, &
         SurfaceFlowDiffusion, &
         SurfaceFlowUpdateSolution, &
         SurfaceFlowRHSFunction, &
         SurfaceFlowComputeMaxDt, &
         SurfaceFlowGetTecplotHeader, &
         SurfaceFlowCreateSurfSubsurfVec, &
         SurfaceFlowCreateSurfSubsurfVecNew, &
         SurfaceFlowUpdateSubsurfSS, &
         SurfaceFlowUpdateSurfBC, &
         SurfaceFlowSurf2SubsurfFlux, &
         SurfaceFlowGetSubsurfProp, &
         SurfaceFlowUpdateAuxVars, &
         SurfaceFlowUpdateSurfState, &
         SurfaceFlowUpdateSubsurfBC, &
         SurfaceFlowUpdateSurfStateNew

contains

! ************************************************************************** !
!> This routine sets up surfaceflow type
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowSetup(surf_realization)

  use Surface_Realization_class
  
  type(surface_realization_type) :: surf_realization

  call SurfaceFlowSetPlotVariables(surf_realization)
  
end subroutine SurfaceFlowSetup

! ************************************************************************** !
!> This routine adds variables to be printed to list
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/30/12
! ************************************************************************** !
subroutine SurfaceFlowSetPlotVariables(surf_realization)
  
  use Surface_Realization_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => surf_realization%output_option%output_variable_list
  
  name = 'H'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_HEAD)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceFlowSetPlotVariables

! ************************************************************************** !
!> This routine computes the internal flux term for the residual under
!! kinematic-wave assumption.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowKinematic(hw_up, &
                                mannings_up, &
                                hw_dn, &
                                mannings_dn, &
                                slope, &
                                length, &
                                option, &
                                vel, &
                                Res)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: slope
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s
  
  PetscReal :: flux       ! units: m^2/s
  
  ! initialize
  flux = 0.d0
  vel  = 0.d0
  
  if (slope<0.d0) then
    vel =  sqrt(dabs(slope))/mannings_up*((hw_up)**(2.d0/3.d0))
    flux=  hw_up*vel
  else
    vel = -sqrt(dabs(slope))/mannings_dn*((hw_dn)**(2.d0/3.d0))
    flux=  hw_dn*vel
  endif

  Res(1) = flux*length

end subroutine SurfaceFlowKinematic

! ************************************************************************** !
!> This routine computes the internal flux term for the residual under
!! diffusion-wave assumption.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceFlowDiffusion(hw_up, &
                                zc_up, &
                                mannings_up, &
                                hw_dn, &
                                zc_dn, &
                                mannings_dn, &
                                dist, &
                                length, &
                                option, &
                                vel, &
                                Res)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: head_up, head_dn
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s

  PetscReal :: flux       ! units: m^2/s
  PetscReal :: Cd
  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: dhead

  ! initialize
  flux = 0.d0
  Cd = 1.0d0

  head_up = hw_up + zc_up
  head_dn = hw_dn + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (hw_up>0.d0) then
      hw_half = hw_up
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    if (hw_dn>0.d0) then
      hw_half = hw_dn
    else
      hw_half = 0.d0
    endif
  endif
  
  dhead=head_up-head_dn
  if(abs(dhead)<eps) then
    dhead=0.d0
    vel = 0.d0
  else
    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)
  endif

  flux = hw_half*vel
  Res(1) = flux*length

end subroutine SurfaceFlowDiffusion

! ************************************************************************** !
!> This routine updates data in module after a successful time step
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/22/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateSolution(surf_realization)

  use Surface_Realization_class
  use Surface_Field_module

  implicit none

  type(surface_realization_type)   :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode                   :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr)

end subroutine SurfaceFlowUpdateSolution

! ************************************************************************** !
!> This routine provides the function evaluation for PETSc TSSolve()
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceFlowRHSFunction(ts,t,xx,ff,surf_realization,ierr)

  use Surface_Realization_class
  use Surface_Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Connection_module
  use Grid_module
  use Coupler_module
  use Surface_Field_module
  use Debug_module
  use Surface_Global_Aux_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  PetscViewer :: viewer

  type(discretization_type), pointer        :: discretization
  type(surface_field_type), pointer         :: surf_field
  type(option_type), pointer                :: option
  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(coupler_type), pointer               :: boundary_condition
  type(coupler_type), pointer               :: source_sink
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)

  PetscInt  :: local_id_up, local_id_dn, local_id
  PetscInt  :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscInt  :: iconn
  PetscInt  :: sum_connection
  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: qsrc, qsrc_flow

  character(len=MAXSTRINGLENGTH)       :: string,string2

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc

  surf_field              => surf_realization%surf_field
  discretization          => surf_realization%discretization
  option                  => surf_realization%option
  patch                   => surf_realization%patch
  grid                    => patch%grid
  surf_global_aux_vars    => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc
  
  surf_realization%iter_count = surf_realization%iter_count+1
  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

  call DiscretizationGlobalToLocal(discretization,xx,surf_field%flow_xx_loc,NFLOWDOF)
  ! Then, update the aux vars
  call SurfaceFlowUpdateAuxVars(surf_realization)

  call VecGetArrayF90(ff,ff_p, ierr)
  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p, ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  ff_p = 0.d0
  Res  = 0.d0

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      if(surf_global_aux_vars(ghosted_id_up)%head(1)<0.d0 .or. &
         surf_global_aux_vars(ghosted_id_dn)%head(1)<0.d0) then
        write(*,*),'In SurfaceFlowFlux: ', surf_global_aux_vars(ghosted_id_up)%head(1), &
          surf_global_aux_vars(ghosted_id_dn)%head(1),ghosted_id_up,ghosted_id_dn
          option%io_buffer='stopping: -ve head values '
          call printErrMsg(option)
      endif

      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          option%io_buffer='Explicit Surface flow not implemented for ' // &
            'Kinematic wave'
          call printErrMsg(option)
        case (DIFFUSION_WAVE)
#if 1
        call SurfaceFlowFlux(surf_global_aux_vars(ghosted_id_up), &
                             zc(ghosted_id_up), &
                             mannings_loc_p(ghosted_id_up), &
                             surf_global_aux_vars(ghosted_id_dn), &
                             zc(ghosted_id_dn), &
                             mannings_loc_p(ghosted_id_dn), &
                             dist, cur_connection_set%area(iconn), &
                             option,vel,Res)
#endif
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(RICHARDS_PRESSURE_DOF,sum_connection) = Res(1)

      vel = patch%internal_velocities(1,sum_connection)
      Res(1) = patch%surf_internal_fluxes(RICHARDS_PRESSURE_DOF,sum_connection)

      if (local_id_up>0) then
        ff_p(local_id_up) = ff_p(local_id_up) - Res(1)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        ff_p(local_id_dn) = ff_p(local_id_dn) + Res(1)/area_p(local_id_dn)
      endif


    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

#if 1
      call SurfaceFlowBCFlux(boundary_condition%flow_condition%itype, &
                         surf_global_aux_vars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)
#endif

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(RICHARDS_PRESSURE_DOF,sum_connection) = Res(1)
      vel = patch%boundary_velocities(1,sum_connection)
      Res(1) = patch%surf_boundary_fluxes(RICHARDS_PRESSURE_DOF,sum_connection)
      
      ff_p(local_id) = ff_p(local_id) + Res(1)/area_p(local_id)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    qsrc_flow = 0.d0
    if(source_sink%flow_condition%rate%itype/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%dataset%rarray(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc = m^3/sec
          qsrc = qsrc_flow*area_p(local_id)
        case(HET_VOL_RATE_SS)
          ! qsrc = m^3/sec
          qsrc = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)*area_p(local_id)
        case default
          option%io_buffer = 'Source/Sink flow condition type not recognized'
          call printErrMsg(option)
      end select
      
      ff_p(local_id) = ff_p(local_id) + qsrc/area_p(local_id)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(ff,ff_p, ierr)
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)

  if (surf_realization%debug%vecview_solution) then
    string = 'Surf_xx_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    string = 'Surf_ff_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(ff,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  endif

end subroutine SurfaceFlowRHSFunction

! ************************************************************************** !
!> This routine maximum allowable 'dt' for surface flow model.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceFlowComputeMaxDt(surf_realization,max_allowable_dt)

  use Water_EOS_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_Global_Aux_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: max_allowable_dt
  PetscReal :: dt

  PetscReal, pointer :: mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc

  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p, ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
        case (DIFFUSION_WAVE)
          call SurfaceFlowFlux(surf_global_aux_vars(ghosted_id_up), &
                               zc(ghosted_id_up), &
                               mannings_loc_p(ghosted_id_up), &
                               surf_global_aux_vars(ghosted_id_dn), &
                               zc(ghosted_id_dn), &
                               mannings_loc_p(ghosted_id_dn), &
                               dist, cur_connection_set%area(iconn), &
                               option,vel,Res)
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(RICHARDS_PRESSURE_DOF,sum_connection) = Res(1)
      if(abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt,dt)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope_dn = dz/dist

      call SurfaceFlowBCFlux(boundary_condition%flow_condition%itype, &
                         surf_global_aux_vars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)
      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(RICHARDS_PRESSURE_DOF,sum_connection) = Res(1)

      if(abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt,dt)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)

end subroutine SurfaceFlowComputeMaxDt

! ************************************************************************** !
!> This routine computes the internal flux term for under
!! diffusion-wave assumption.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceFlowFlux(surf_global_aux_var_up, &
                         zc_up, &
                         mannings_up, &
                         surf_global_aux_var_dn, &
                         zc_dn, &
                         mannings_dn, &
                         dist, &
                         length, &
                         option, &
                         vel, &
                         Res)

  use Surface_Global_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  type(surface_global_auxvar_type) :: surf_global_aux_var_up
  type(surface_global_auxvar_type) :: surf_global_aux_var_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn


  PetscReal :: head_up, head_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s

  PetscReal :: flux       ! units: m^2/s
  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: dhead

  ! initialize
  flux = 0.d0

  head_up = surf_global_aux_var_up%head(1) + zc_up
  head_dn = surf_global_aux_var_dn%head(1) + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (surf_global_aux_var_up%head(1)>eps) then
      hw_half = surf_global_aux_var_up%head(1)
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    if (surf_global_aux_var_dn%head(1)>eps) then
      hw_half = surf_global_aux_var_dn%head(1)
    else
      hw_half = 0.d0
    endif
  endif
  
  dhead=head_up-head_dn
  if(abs(dhead)<eps) then
    dhead=0.d0
    vel = 0.d0
  else
    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)
  endif

  flux = hw_half*vel
  Res(TH_PRESSURE_DOF) = flux*length

end subroutine SurfaceFlowFlux

! ************************************************************************** !
!> This routine computes the boundary term surface water equation
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceFlowBCFlux(ibndtype, &
                           surf_global_aux_var, &
                           slope, &
                           mannings, &
                           length, &
                           option, &
                           vel, &
                           Res)

  use Option_module
  use Surface_Global_Aux_module
  
  implicit none

  type(option_type) :: option
  type(surface_global_auxvar_type) :: surf_global_aux_var
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt  :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: Res(1:option%nflowdof) 

  PetscInt :: pressure_bc_type
  PetscReal :: head

  flux = 0.d0
  vel = 0.d0
  
  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_aux_var%head(1)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
      else
        vel = -sqrt(dabs(slope))/mannings*((head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = head*vel
  Res(RICHARDS_PRESSURE_DOF) = flux*length

end subroutine SurfaceFlowBCFlux

! ************************************************************************** !
!> This routine updates auxiliary variables
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceFlowUpdateAuxVars(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module
  use Surface_Global_Aux_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_ss(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(1)
  PetscReal :: xxss(1)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc
  surf_global_aux_vars_ss => patch%surf_aux%SurfaceGlobal%aux_vars_ss
  
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    surf_global_aux_vars(ghosted_id)%head(1) = xx_loc_p(ghosted_id)
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)
   
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)
  ! Boundary aux vars
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

        select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,HET_DIRICHLET)
            xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(1) = xx_loc_p(ghosted_id)
        end select
      
      surf_global_aux_vars_bc(sum_connection)%head(1) = xxbc(1)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/Sink aux vars
  ! source/sinks
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      xxss = xx_loc_p(ghosted_id)
      surf_global_aux_vars_ss(sum_connection)%head(1) = xxss(1)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

end subroutine SurfaceFlowUpdateAuxVars

! ************************************************************************** !
!> This routine surface flow tecplot file header
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/29/12
! ************************************************************************** !
function SurfaceFlowGetTecplotHeader(surf_realization,icolumn)

  use Surface_Realization_class
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: SurfaceFlowGetTecplotHeader
  type(surface_realization_type) :: surf_realization
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: string, string2

  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [m]"'')')
  endif
  string = trim(string) // trim(string2)
#ifdef GLENN_NEW_IO
  !call OutputOptionAddPlotVariable(realization%output_option,PRESSURE, &
  !                           ZERO_INTEGER,ZERO_INTEGER)
#endif

  SurfaceFlowGetTecplotHeader = string

end function SurfaceFlowGetTecplotHeader

! ************************************************************************** !
!> This routine get soil properties of the top-most soil layer from the 
!! subsurface domain.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceFlowGetSubsurfProp(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: qsrc_p(:),vec_p(:)
  PetscInt :: local_id,iconn,sum_connection,ghosted_id
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: yy_p(:)
  PetscReal, pointer :: zz_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: dist_p(:)
  PetscReal :: dist_x,dist_y,dist_z,dist
  PetscInt :: size
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
  
  ! Update the surface BC
  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  sum_connection = 0
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      cur_connection_set => coupler%connection_set
      
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then

        ! perm_x
        call VecGetArrayF90(field%perm_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_xx_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                             surf_field%subsurf_temp_vec_1dof, &
                             surf_field%perm_xx, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_xx, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! perm_y
        call VecGetArrayF90(field%perm_yy_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_yy_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%perm_yy, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_yy, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! perm_z
        call VecGetArrayF90(field%perm_zz_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_zz_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%perm_zz, &
                            INSERT_VALUES, SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_zz, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! por
        call VecGetArrayF90(field%porosity_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%porosity_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%por, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%por, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! icap: ID of saturation function
        call VecGetArrayF90(field%icap_loc,icap_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          !vec_p(iconn)=patch%sat_func_id(ghosted_id)
          vec_p(iconn)=icap_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%icap_loc,icap_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%icap_loc, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%icap_loc, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
        ! x
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%x(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_xx, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_xx, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! y
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%y(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_yy, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_yy, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! z
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%z(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_zz, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_zz, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    coupler => coupler%next
  enddo

  call VecGetArrayF90(surf_field%subsurf_xx,xx_p,ierr)
  call VecGetArrayF90(surf_field%subsurf_yy,yy_p,ierr)
  call VecGetArrayF90(surf_field%subsurf_zz,zz_p,ierr)
  call VecGetArrayF90(surf_field%perm_xx,perm_xx_p,ierr)
  call VecGetArrayF90(surf_field%perm_yy,perm_yy_p,ierr)
  call VecGetArrayF90(surf_field%perm_zz,perm_zz_p,ierr)
  call VecGetArrayF90(surf_field%Dq,Dq_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist_gravity,dist_p,ierr)

  do local_id=1,surf_grid%nlmax
    dist_x = -(xx_p(local_id) - surf_grid%x(local_id))
    dist_y = -(yy_p(local_id) - surf_grid%y(local_id))
    dist_z = -(zz_p(local_id) - surf_grid%z(local_id))
      
    dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z)

    Dq_p(local_id) = (perm_xx_p(local_id)*abs(dist_x)/dist + &
                      perm_yy_p(local_id)*abs(dist_y)/dist + &
                      perm_zz_p(local_id)*abs(dist_z)/dist)/dist
    dist_p(local_id) = (dist_x*option%gravity(1)+ &
                        dist_y*option%gravity(2)+ &
                        dist_z*option%gravity(3))
  enddo

  call VecRestoreArrayF90(surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_xx,xx_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_yy,yy_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_zz,zz_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_xx,perm_xx_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_yy,perm_yy_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_zz,perm_zz_p,ierr)
  call VecRestoreArrayF90(surf_field%Dq,Dq_p,ierr)

  surf_realization%first_time=PETSC_FALSE

end subroutine SurfaceFlowGetSubsurfProp

! ************************************************************************** !
!> This routine prescribes the strength of source/sink between surface and
!! subsurface volume.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 06/06/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateSubsurfSS(realization,surf_realization,dt)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_loc_p(:),vec_p(:)
  PetscReal :: dt
  PetscInt :: local_id
  PetscInt :: iconn
  PetscReal :: den
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)

  call density(option%reference_temperature,option%reference_pressure,den)

  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then
      
        coupler_found = PETSC_TRUE
        
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                             surf_field%exchange_subsurf_2_surf, &
                             surf_field%subsurf_temp_vec_1dof, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%exchange_subsurf_2_surf, &
                           surf_field%subsurf_temp_vec_1dof, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,coupler%connection_set%num_connections
          coupler%flow_aux_real_var(ONE_INTEGER,iconn)=-vec_p(iconn)/dt*den
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)

        call VecSet(surf_field%exchange_subsurf_2_surf,0.d0,ierr)
      endif

    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_ss.'
    call printErrMsg(option)
  endif
  
end subroutine SurfaceFlowUpdateSubsurfSS

! ************************************************************************** !
!> This routine scatters pressure values for first soil layer from subsurface
!! model to surface model. When this routine is called for the first time,
!! following subsurface properties are send from subsurface model to surface
!! model:
!!  - Permeabilities in X,Y,Z directions.
!!  - Distance of soil control volume to center of surface cell (DX,DY,DZ).
!!  - Porosity.
!!  - Saturation function ID of soil control volume.
!!
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 06/06/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateSurfBC(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: qsrc_p(:),vec_p(:)
  PetscInt :: local_id,iconn,sum_connection,ghosted_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: yy_p(:)
  PetscReal, pointer :: zz_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: dist_p(:)
  PetscReal :: dist_x,dist_y,dist_z,dist
  PetscInt :: size
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
  
  ! Update the surface BC
  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  sum_connection = 0
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      cur_connection_set => coupler%connection_set
      
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then

        call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)

        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)

        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    coupler => coupler%next
  enddo

end subroutine SurfaceFlowUpdateSurfBC

! ************************************************************************** !
!> This routine computes flux between surface and subsurface model.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 06/06/12
! ************************************************************************** !
subroutine SurfaceFlowSurf2SubsurfFlux(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Water_EOS_module
  use Saturation_Function_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal :: den          ! density      [kg/m^3]
  PetscInt :: local_id,iconn
  
  PetscReal, pointer :: hw_p(:)   ! head [m]
  PetscReal, pointer :: press_sub_p(:) ! Pressure [Pa]
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: mass_p(:)
  PetscReal, pointer :: area_p(:)
  PetscReal, pointer :: dist_p(:)
  
  PetscReal :: press_surf
  PetscReal :: dphi
  PetscReal :: press
  PetscReal :: sat
  PetscReal :: kr
  PetscReal :: ds_dp
  PetscReal :: dkr_dp
  PetscReal :: por
  PetscReal :: perm
  PetscBool :: saturated
  PetscReal :: sat_pressure
  PetscReal :: pw
  PetscReal :: visl
  PetscReal :: dvis_dp
  PetscReal :: dvis_dt
  PetscReal :: v_darcy
  PetscReal :: v_darcy_max
  PetscReal :: gravity
  PetscReal :: press_up, press_dn
    
  PetscBool :: coupler_found = PETSC_FALSE
  PetscBool :: v_darcy_limit

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  call density(option%reference_temperature,option%reference_pressure,den)

  call VecGetArrayF90(surf_field%press_subsurf,press_sub_p,ierr)
  call VecGetArrayF90(surf_field%flow_xx,hw_p,ierr)
  call VecGetArrayF90(surf_field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(surf_field%Dq,Dq_p,ierr)
  call VecGetArrayF90(surf_field%exchange_subsurf_2_surf,mass_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  ! Update the surface BC
  coupler_list => surf_patch%source_sinks
  coupler => coupler_list%first
  v_darcy_max=0.d0
  v_darcy_limit=PETSC_FALSE
  
  !
  !            SURFACE
  !              (DN)
  !   ---------------------------------
  !           SUBSURFACE         ///\\\
  !              (UP)
  !

  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if(StringCompare(coupler%name,'from_subsurface_ss')) then

      if (.not.associated(coupler%flow_aux_real_var)) then
        option%io_buffer='SURF_SOURCE_SINK: from_subsurface_ss was does not ' //&
          ' a heterogeneous flow condition associated with it'
        call printErrMsg(option)
      endif

      if(coupler%flow_condition%rate%itype/=HET_VOL_RATE_SS) then
        option%io_buffer = 'Flow condition from_subsurface_ss should be ' // &
          'heterogeneous_volumetric_rate'
        call printErrMsg(option)
      endif
      
      do local_id=1,surf_grid%nlmax
        press_surf=hw_p(local_id)*(abs(option%gravity(3)))*den+option%reference_pressure

        press_up = press_sub_p(local_id)
        press_dn = press_surf
        gravity = dist_p(local_id)*den
        
        dphi = press_up - press_dn + gravity
        
        ! check if there is standing water on the surface
        if (dphi<0.d0 .and. press_dn - option%reference_pressure<eps) then
          dphi= 0.d0
        endif
        
        ! exfiltration will only occur if subsurface pressure is greater
        ! than reference pressure
        if (dphi>=0 .and. press_up-option%reference_pressure<eps ) then
          dphi = 0.d0
        endif
        
        if(dphi>0.d0) then
          press = press_up
        else
          press = press_dn
        endif
        
        call SaturationFunctionCompute( &
          press,sat,kr,ds_dp,dkr_dp,&
          patch%saturation_function_array(int(icap_loc_p(local_id)))%ptr, &
          0.d0,0.d0,saturated,option)
        
        if(saturated) then
          pw = press
        else
          pw = option%reference_pressure
        endif
                                           
        call psat(option%reference_temperature,sat_pressure,ierr)
        call VISW(option%reference_temperature,pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr)

        v_darcy = Dq_p(local_id)*kr/visl*dphi
        if (v_darcy<=0.d0) then
          ! Flow is happening from surface to subsurface
          if ( abs(v_darcy) > hw_p(local_id)/option%surf_flow_dt ) then
            v_darcy = -hw_p(local_id)/option%surf_flow_dt
            v_darcy_limit=PETSC_TRUE
          endif
        else
          ! Exfiltration is occuring
          !v_darcy=0.d0
        endif
        
        mass_p(local_id)=mass_p(local_id)+v_darcy*area_p(local_id)*option%surf_flow_dt
        !coupler%flow_aux_real_var(ONE_INTEGER,local_id)=v_darcy
        coupler%flow_aux_real_var(ONE_INTEGER,local_id)=0.d0
        hw_p(local_id) = hw_p(local_id) + v_darcy*option%surf_flow_dt
        if(hw_p(local_id)<1.d-15) hw_p(local_id) = 0.d0
        if(abs(v_darcy)>v_darcy_max) v_darcy_max=v_darcy
      enddo

    endif
    
    coupler => coupler%next
  enddo
  
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)
  call VecRestoreArrayF90(surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call VecRestoreArrayF90(surf_field%exchange_subsurf_2_surf,mass_p,ierr)
  call VecRestoreArrayF90(surf_field%Dq,Dq_p,ierr)  
  call VecRestoreArrayF90(surf_field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%flow_xx,hw_p,ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf,press_sub_p,ierr)

end subroutine SurfaceFlowSurf2SubsurfFlux

! ************************************************************************** !
!> This routine creates a MPI vector to keep accumulated volume of water
!! exchanged between surface and subsurface model.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 10/05/12
! ************************************************************************** !
subroutine SurfaceFlowCreateSurfSubsurfVec(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_loc_p(:),vec_p(:)
  PetscReal :: den          ! density      [kg/m^3]
  PetscInt :: local_id,i
  PetscInt :: iconn
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then
        if(coupler%flow_condition%rate%itype/=HET_MASS_RATE_SS) then
          option%io_buffer = 'Flow condition from_surface_ss should be ' // &
            'heterogeneous_mass_rate'
          call printErrMsg(option)
        endif
        coupler_found = PETSC_TRUE
        if(surf_realization%first_time) then
          call VecCreate(option%mycomm,surf_field%subsurf_temp_vec_1dof,ierr)
          call VecSetSizes(surf_field%subsurf_temp_vec_1dof, &
                           coupler%connection_set%num_connections,PETSC_DECIDE,ierr)
          call VecSetFromOptions(surf_field%subsurf_temp_vec_1dof,ierr)
        endif
      endif
    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_ss.'
    call printErrMsg(option)
  endif

end subroutine SurfaceFlowCreateSurfSubsurfVec

! ************************************************************************** !
!> This routine creates a MPI vector to exchanged data from subsurface model
!! to surface model.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 07/29/13
! ************************************************************************** !
subroutine SurfaceFlowCreateSurfSubsurfVecNew(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_loc_p(:),vec_p(:)
  PetscReal :: den          ! density      [kg/m^3]
  PetscInt :: local_id,i
  PetscInt :: iconn
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  coupler_list => patch%boundary_conditions
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      ! Find the Source/Sink from the list of Source/Sinks
      if(StringCompare(coupler%name,'from_surface_bc')) then
        if(coupler%flow_condition%pressure%itype/=HET_SURF_SEEPAGE_BC) then
          option%io_buffer = 'Flow condition from_surface_bc should be ' // &
            'heterogeneous_surface_seepage'
          call printErrMsg(option)
        endif
        coupler_found = PETSC_TRUE
        if(surf_realization%first_time) then
          call VecCreate(option%mycomm,surf_field%subsurf_temp_vec_1dof,ierr)
          call VecSetSizes(surf_field%subsurf_temp_vec_1dof, &
                           coupler%connection_set%num_connections,PETSC_DECIDE,ierr)
          call VecSetFromOptions(surf_field%subsurf_temp_vec_1dof,ierr)
          call VecSet(surf_field%subsurf_temp_vec_1dof,0.d0,ierr)

          call VecDuplicate(surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_avg_vdarcy, ierr)
        endif
      endif
    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_bc.'
    call printErrMsg(option)
  endif

end subroutine SurfaceFlowCreateSurfSubsurfVecNew

! ************************************************************************** !
!> This routine gets updated values of standing water at the end of 
!! subsurface-flow model timestep.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 07/30/13
! ************************************************************************** !
subroutine SurfaceFlowUpdateSurfState(realization, surf_realization, dt)

  use Connection_module
  use Coupler_module
  use Discretization_module
  use DM_Kludge_module
  use Grid_module
  use Option_module
  use Patch_module
  use Realization_class
  use Realization_Base_class
  use String_module
  use Surface_Field_module
  use Surface_Realization_class
  use Water_EOS_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization
  PetscReal                      :: dt

  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(connection_set_type), pointer  :: cur_connection_set
  type(dm_ptr_type), pointer          :: dm_ptr
  type(grid_type),pointer             :: grid,surf_grid
  type(option_type), pointer          :: option
  type(patch_type),pointer            :: patch,surf_patch
  type(surface_field_type),pointer    :: surf_field

  PetscInt                :: count
  PetscInt                :: ghosted_id
  PetscInt                :: iconn
  PetscInt                :: local_id
  PetscInt                :: sum_connection

  PetscReal               :: den
  PetscReal, pointer      :: avg_vdarcy_p(:)   ! avg darcy velocity [m/s]
  PetscReal, pointer      :: hw_p(:)           ! head [m]
  PetscReal, pointer      :: surfpress_p(:)
  PetscErrorCode          :: ierr

  PetscBool :: coupler_found = PETSC_FALSE

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization, ONEDOF)
  option     => realization%option
  patch      => realization%patch
  surf_field => surf_realization%surf_field
  surf_grid  => surf_realization%discretization%grid
  
  sum_connection = 0

  coupler_list => patch%boundary_conditions
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      cur_connection_set => coupler%connection_set

      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_bc')) then
      
        coupler_found = PETSC_TRUE

        ! Save the updated surface-pressure BC from subsurface model into Vec
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,surfpress_p,ierr)
        do iconn = 1,coupler%connection_set%num_connections
          sum_connection = sum_connection + 1
          surfpress_p(iconn) = coupler%flow_aux_real_var(ONE_INTEGER,sum_connection)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,surfpress_p,ierr)

      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_bc.'
    call printErrMsg(option)
  endif

  call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                       surf_field%subsurf_temp_vec_1dof, &
                       surf_field%work, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                     surf_field%subsurf_temp_vec_1dof, &
                     surf_field%work, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)

  call density(option%reference_temperature,option%reference_pressure,den)

  call VecGetArrayF90(surf_field%flow_xx, hw_p, ierr)
  call VecGetArrayF90(surf_field%work, surfpress_p, ierr)
  count = 0
  do ghosted_id = 1,surf_grid%ngmax

    local_id = surf_grid%nL2G(ghosted_id)
    if(local_id <= 0) cycle

    count = count + 1
    hw_p(ghosted_id) = (surfpress_p(count)-option%reference_pressure)/ &
                        (abs(option%gravity(3)))/den
    if(hw_p(ghosted_id)<1.d-15) hw_p(ghosted_id) = 0.d0

  enddo
  call VecRestoreArrayF90(surf_field%flow_xx, hw_p, ierr)
  call VecRestoreArrayF90(surf_field%work, surfpress_p, ierr)

end subroutine SurfaceFlowUpdateSurfState

! ************************************************************************** !
!> This routine gets updated values of standing water at the end of 
!! subsurface-flow model timestep.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 07/30/13
! ************************************************************************** !
subroutine SurfaceFlowUpdateSurfStateNew(surf_realization)

  use Connection_module
  use Coupler_module
  use Discretization_module
  use DM_Kludge_module
  use Grid_module
  use Option_module
  use Patch_module
  use Realization_class
  use Realization_Base_class
  use String_module
  use Surface_Field_module
  use Surface_Realization_class
  use Water_EOS_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(surface_realization_type) :: surf_realization

  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(connection_set_type), pointer  :: cur_connection_set
  type(dm_ptr_type), pointer          :: dm_ptr
  type(grid_type),pointer             :: grid,surf_grid
  type(option_type), pointer          :: option
  type(patch_type),pointer            :: patch,surf_patch
  type(surface_field_type),pointer    :: surf_field

  PetscInt                :: count
  PetscInt                :: ghosted_id
  PetscInt                :: iconn
  PetscInt                :: local_id
  PetscInt                :: sum_connection

  PetscReal               :: den
  PetscReal, pointer      :: avg_vdarcy_p(:)   ! avg darcy velocity [m/s]
  PetscReal, pointer      :: hw_p(:)           ! head [m]
  PetscReal, pointer      :: surfpress_p(:)
  PetscErrorCode          :: ierr

  PetscBool :: coupler_found = PETSC_FALSE

  option     => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_grid  => surf_realization%discretization%grid
  
  call density(option%reference_temperature,option%reference_pressure,den)

  call VecGetArrayF90(surf_field%flow_xx, hw_p, ierr)
  call VecGetArrayF90(surf_field%press_subsurf, surfpress_p, ierr)
  count = 0
  do ghosted_id = 1,surf_grid%ngmax

    local_id = surf_grid%nG2L(ghosted_id)
    if(local_id <= 0) cycle

    count = count + 1
    hw_p(ghosted_id) = (surfpress_p(count)-option%reference_pressure)/ &
                        (abs(option%gravity(3)))/den
    if(hw_p(ghosted_id)<1.d-15) hw_p(ghosted_id) = 0.d0

  enddo
  call VecRestoreArrayF90(surf_field%flow_xx, hw_p, ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf, surfpress_p, ierr)

  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%flow_xx, &
                                   surf_field%flow_xx_loc, &
                                   NFLOWDOF)
  call SurfaceFlowUpdateAuxVars(surf_realization)

end subroutine SurfaceFlowUpdateSurfStateNew

! ************************************************************************** !
!> This routine updates the pressure BC for subsurface model at the end of
!!  surface-flow model timestep.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 07/30/13
! ************************************************************************** !
subroutine SurfaceFlowUpdateSubsurfBC(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: hw_p(:)
  PetscReal, pointer :: surfpress_p(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: count
  PetscInt :: iconn
  PetscReal :: den
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)

  call density(option%reference_temperature,option%reference_pressure,den)

  call VecGetArrayF90(surf_field%flow_xx, hw_p, ierr)
  call VecGetArrayF90(surf_field%work, surfpress_p, ierr)
  count = 0
  do ghosted_id = 1,surf_grid%ngmax
    local_id = surf_grid%nL2G(ghosted_id)
    if(local_id <= 0) cycle
    count = count + 1
    surfpress_p(count) = hw_p(ghosted_id)*(abs(option%gravity(3)))*den + &
                        option%reference_pressure
  enddo
  call VecRestoreArrayF90(surf_field%work, surfpress_p, ierr)
  call VecRestoreArrayF90(surf_field%flow_xx, hw_p, ierr)


  coupler_list => patch%boundary_conditions
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_bc')) then

        coupler_found = PETSC_TRUE

        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                             surf_field%work, &
                             surf_field%subsurf_temp_vec_1dof, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%work, &
                           surf_field%subsurf_temp_vec_1dof, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,surfpress_p,ierr)
        do iconn = 1,coupler%connection_set%num_connections
          coupler%flow_aux_real_var(ONE_INTEGER,iconn) = surfpress_p(iconn)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,surfpress_p,ierr)

        call VecSet(surf_field%exchange_subsurf_2_surf,0.d0,ierr)
      endif

    endif

    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_bc.'
    call printErrMsg(option)
  endif

end subroutine SurfaceFlowUpdateSubsurfBC

end module Surface_Flow_module

#endif
