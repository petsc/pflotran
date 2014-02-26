#ifdef SURFACE_FLOW

module Surface_TH_module

  use Surface_Global_Aux_module
  use Surface_TH_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceTHSetup, &
         SurfaceTHRHSFunction, &
         SurfaceTHComputeMaxDt, &
         SurfaceTHUpdateAuxVars, &
         SurfaceTHUpdateSolution, &
         SurfaceTHUpdateTemperature, &
         SurfaceTHUpdateSurfState

contains

! ************************************************************************** !

subroutine SurfaceTHSetup(surf_realization)
  ! 
  ! This routine sets up surface_TH_type
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
 
  implicit none
  
  type(surface_realization_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(coupler_type), pointer :: initial_condition
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase
  
  
  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
    
  patch%surf_aux%SurfaceTH => SurfaceTHAuxCreate(option)

  ! allocate auxvar data structures for all grid cells
  allocate(Surf_TH_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SurfaceTHAuxVarInit(Surf_TH_auxvars(ghosted_id),option)
  enddo

  patch%surf_aux%SurfaceTH%auxvars => Surf_TH_auxvars
  patch%surf_aux%SurfaceTH%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_conditions%first

  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then 
    allocate(Surf_TH_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_auxvars_bc(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%auxvars_bc => Surf_TH_auxvars_bc
  endif
  patch%surf_aux%SurfaceTH%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(Surf_TH_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_auxvars_ss(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%auxvars_ss => Surf_TH_auxvars_ss
  endif
  patch%surf_aux%SurfaceTH%num_aux_ss = sum_connection

  call SurfaceTHSetPlotVariables(surf_realization)

end subroutine SurfaceTHSetup

! ************************************************************************** !

subroutine SurfaceTHSetPlotVariables(surf_realization)
  ! 
  ! This routine adds variables to be printed to list
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 
  
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

  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_TEMPERATURE)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceTHSetPlotVariables

! ************************************************************************** !

subroutine SurfaceTHRHSFunction(ts,t,xx,ff,surf_realization,ierr)
  ! 
  ! This routine provides the function evaluation for PETSc TSSolve()
  ! Author: Gautam Bisht, LBNL
  ! 

  use EOS_Water_module
  use Connection_module
  use Surface_Realization_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(coupler_type), pointer               :: source_sink
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set

  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt :: istart, iend

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: qsrc, qsrc_flow
  PetscReal :: esrc
  PetscReal :: den

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH)       :: string,string2

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss

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

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Then, update the aux vars
  ! RTM: This includes calculation of the accumulation terms, correct?
  call SurfaceTHUpdateAuxVars(surf_realization)
  ! override flags since they will soon be out of date  
  patch%surf_aux%SurfaceTH%auxvars_up_to_date = PETSC_FALSE

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
      
      call SurfaceTHFlux(surf_auxvars(ghosted_id_up), &
                         surf_global_auxvars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(TH_PRESSURE_DOF)
      patch%surf_internal_fluxes(TH_TEMPERATURE_DOF,sum_connection) = Res(TH_TEMPERATURE_DOF)

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) - Res(:)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
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
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         surf_auxvars_bc(sum_connection), &
                         surf_global_auxvars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(TH_PRESSURE_DOF)
      patch%surf_boundary_fluxes(TH_TEMPERATURE_DOF,sum_connection) = Res(TH_TEMPERATURE_DOF)
      
      iend = local_id_dn*option%nflowdof
      istart = iend-option%nflowdof+1
      ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    if(source_sink%flow_condition%rate%itype/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%dataset%rarray(1)
      
    if (source_sink%flow_condition%rate%itype == ENERGY_RATE_SS) &
      esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)

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
      
      esrc = 0.d0
      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (ENERGY_RATE_SS)
          esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (HET_ENERGY_RATE_SS)
          esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ff_p(istart) = ff_p(istart) + qsrc/area_p(local_id)
      ! RTM: TODO: What should the density term and specific heat capactiy be
      ! in the freezing case?
      ! I think using the weighted average of liquid and ice densities and Cwi 
      ! is correct here, but I should check.
      ff_p(iend) = ff_p(iend) + esrc + &
                    surf_global_auxvars_ss(local_id)%den_kg(1)* &
                    (surf_global_auxvars_ss(local_id)%temp(1) + 273.15d0)* &
                    surf_auxvars(local_id)%Cwi* &
                    qsrc/area_p(local_id)
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

end subroutine SurfaceTHRHSFunction

! ************************************************************************** !

subroutine SurfaceTHComputeMaxDt(surf_realization,max_allowable_dt)
  ! 
  ! This routine maximum allowable 'dt' for explicit time scheme.
  ! Author: Gautam Bisht, LBNL
  ! 

  use EOS_Water_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
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

  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)

  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
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

  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc

  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p, ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10
  vel = 0.d0

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
      
      call SurfaceTHFlux(surf_auxvars(ghosted_id_up), &
                         surf_global_auxvars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(:,sum_connection) = Res(:)
      if(abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt, dt)
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
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         surf_auxvars_bc(sum_connection), &
                         surf_global_auxvars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(:,sum_connection) = Res(:)

      if(abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt, dt)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)
  
end subroutine SurfaceTHComputeMaxDt

! ************************************************************************** !

subroutine SurfaceTHFlux(surf_auxvar_up, &
                         surf_global_auxvar_up, &
                         zc_up, &
                         mannings_up, &
                         surf_auxvar_dn, &
                         surf_global_auxvar_dn, &
                         zc_dn, &
                         mannings_dn, &
                         dist, &
                         length, &
                         option, &
                         vel, &
                         Res)
  ! 
  ! This routine computes the internal flux term for under
  ! diffusion-wave assumption.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 08/03/12
  ! 

  use Surface_TH_Aux_module
  use Surface_Global_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_auxvar_up
  type(Surface_TH_auxvar_type) :: surf_auxvar_dn
  type(surface_global_auxvar_type) :: surf_global_auxvar_up
  type(surface_global_auxvar_type) :: surf_global_auxvar_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn

  PetscReal :: head_up, head_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s
  
  PetscReal :: flux       ! units: m^2/s
  PetscReal :: hw_half
  PetscReal :: hw_liq_half
  PetscReal :: mannings_half
  PetscReal :: unfrozen_fraction_half
  PetscReal :: dhead
  PetscReal :: den_aveg
  PetscReal :: temp_half
  PetscReal :: dtemp
  PetscReal :: Cw
  PetscReal :: k_therm

  ! initialize
  flux = 0.d0

  ! Flow equation
  head_up = surf_global_auxvar_up%head(1) + zc_up
  head_dn = surf_global_auxvar_dn%head(1) + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    temp_half = surf_global_auxvar_up%temp(1) + 273.15d0
    unfrozen_fraction_half = surf_auxvar_up%unfrozen_fraction
    if (surf_global_auxvar_up%head(1)>eps) then
      hw_half = surf_global_auxvar_up%head(1)
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    temp_half = surf_global_auxvar_dn%temp(1) + 273.15d0
    unfrozen_fraction_half = surf_auxvar_dn%unfrozen_fraction
    if (surf_global_auxvar_dn%head(1)>eps) then
      hw_half = surf_global_auxvar_dn%head(1)
    else
      hw_half = 0.d0
    endif
  endif

  ! Find pressure head at interface for LIQUID fraction
  hw_liq_half = unfrozen_fraction_half*hw_half

  ! Compute pressure difference
  dhead = head_up - head_dn

  if (abs(dhead) < eps) then
    dhead = 0.d0
    vel = 0.d0
  else
    ! RTM: We modify the term raised to the power 2/3 (the "hydraulic radius") 
    ! by the (upwinded) unfrozen fraction.  For a wide rectangular channel, 
    ! hydraulic radius (which is a measure of the "efficiency" of the channel) 
    ! is often taken to be the flow depth, so I believe this makes sense. (?)
    ! The actual total head term ('hw_half' here) is NOT modified by the 
    ! unfrozen fraction: though the ice is immobile, its weight does 
    ! contribute to the pressure head.
    vel = (hw_liq_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)

     !RTM: Original code for when freezing is not considered is
!    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
!          dhead/(abs(dhead)**(1.d0/2.d0))* &
!          1.d0/(dist**0.5d0)
  endif

  flux = hw_liq_half*vel
  Res(TH_PRESSURE_DOF) = flux*length
  
  ! Temperature equation
  ! RTM: k_therm is the weighted average of the liquid and ice thermal 
  ! conductivities.  For the density and specific heat capacity in the 
  ! advection term, we want these for liquid water ONLY, as the ice portion 
  ! is immobile and thus should not make up part of the advection term. We 
  ! also multiply the ponded water depth (hw_half) by the unfrozen fraction 
  ! in the advection term but NOT the conduction term.
  ! We do the same in SurfaceTHBCFlux().

  ! Average density
  ! Here we only consider the LIQUID fraction.
  den_aveg = (surf_auxvar_up%den_water_kg + &
              surf_auxvar_dn%den_water_kg)/2.d0
  den_aveg = (surf_global_auxvar_up%den_kg(1) + &
              surf_global_auxvar_dn%den_kg(1))/2.d0
  ! Temperature difference
  dtemp = surf_global_auxvar_up%temp(1) - surf_global_auxvar_dn%temp(1)

  ! Note, Cw and k_therm are same for up and downwind
  Cw = surf_auxvar_up%Cw
  k_therm = surf_auxvar_up%k_therm
  
  ! Unfrozen fraction multiplies hw_half in advection term, but does NOT affect the 
  ! conduction therm.  
  ! RTM: Brookfield et al. 2009 also has dispersion term, which we are not using.
  Res(TH_TEMPERATURE_DOF) = (den_aveg*vel*temp_half*Cw*hw_liq_half + &
                             k_therm*dtemp/dist*hw_half)*length

end subroutine SurfaceTHFlux

! ************************************************************************** !

subroutine SurfaceTHBCFlux(ibndtype, &
                           surf_auxvar, &
                           surf_global_auxvar, &
                           slope, &
                           mannings, &
                           length, &
                           option, &
                           vel, &
                           Res)
  ! 
  ! This routine computes flux for boundary cells.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Option_module
  
  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_auxvar
  type(surface_global_auxvar_type) :: surf_global_auxvar
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt  :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: den_aveg

  PetscInt :: pressure_bc_type
  PetscReal :: head
  PetscReal :: head_liq

  flux = 0.d0
  vel = 0.d0
  
  ! RTM: I've multiplied the head (ponded water depth, actually) by the 
  ! unfrozen fraction.  I believe this makes sense, but I should think a bit 
  ! more about what a "zero gradient" condition means in the case of freezing
  ! surface water.

  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_auxvar%head(1)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
        head_liq = 0.d0
      else
        head_liq = surf_auxvar%unfrozen_fraction * head
        vel = -sqrt(dabs(slope))/mannings*(head_liq**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = head_liq*vel
  Res(TH_PRESSURE_DOF) = flux*length

  ! Temperature
  ! RTM: See note about in SufaceTHFlux() about how frozen/unfrozen are handled here.
  Res(TH_TEMPERATURE_DOF) = surf_global_auxvar%den_kg(1)* &
                            (surf_global_auxvar%temp(1) + 273.15d0)* &
                            surf_auxvar%Cwi* &
                            vel*head_liq*length

end subroutine SurfaceTHBCFlux

! ************************************************************************** !

subroutine SurfaceTHUpdateAuxVars(surf_realization)
  ! 
  ! This routine updates auxiliary variables
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr
  PetscReal :: den

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_th_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_th_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_th_auxvars_ss => patch%surf_aux%SurfaceTH%auxvars_ss
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss
  
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call SurfaceTHAuxVarCompute(xx_loc_p(istart:iend), &
                                surf_th_auxvars(ghosted_id), &
                                surf_global_auxvars(ghosted_id), &
                                option)
    ! [rho*h*T*Cwi]
    xx_loc_p(istart+1) = surf_global_auxvars(ghosted_id)%den_kg(1)* &
                         xx_loc_p(istart)* &
                         (surf_global_auxvars(ghosted_id)%temp(1) + 273.15d0)* &
                         surf_th_auxvars(ghosted_id)%Cwi
  enddo
   
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

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,HET_DIRICHLET)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      surf_global_auxvars_bc(sum_connection)%temp(1) = xxbc(2)
      call SurfaceTHAuxVarCompute(xxbc, &
                                  surf_th_auxvars_bc(sum_connection), &
                                  surf_global_auxvars_bc(sum_connection), &
                                  option)

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

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1

      if (associated(source_sink%flow_condition%temperature)) then
        if(source_sink%flow_condition%temperature%itype/=HET_DIRICHLET) then
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        else
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        endif
      else
        tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+1)
        tsrc1 = surf_global_auxvars(ghosted_id)%temp(1)
      endif

      xxss = xx_loc_p(istart:iend)
      xxss(2) = tsrc1

      surf_global_auxvars_ss(sum_connection)%temp(1) = tsrc1
      call SurfaceTHAuxVarCompute(xxss, &
                                  surf_th_auxvars_ss(sum_connection), &
                                  surf_global_auxvars_ss(sum_connection), &
                                  option)
    enddo
    source_sink => source_sink%next
  enddo

  patch%surf_aux%SurfaceTH%auxvars_up_to_date = PETSC_TRUE

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

end subroutine SurfaceTHUpdateAuxVars

! ************************************************************************** !

subroutine SurfaceTHUpdateTemperature(surf_realization)
  ! 
  ! This routine updates the temperature after TSSolve.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/25/13
  ! 

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module
  use EOS_Water_module

  implicit none

  type(surface_realization_type) :: surf_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: temp
  PetscInt :: iter
  PetscInt :: niter
  PetscReal :: den
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss
  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc

  ! Number of iterations to solve for T^{t+1}
  ! T^{t+1,m} = (rho Cwi hw T)^{t+1} / rho^{t+1,m-1} Cw)^{t} (hw)^{t+1}
  ! niter = max(m)
  niter = 20

  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr)

  ! Update internal aux vars
  do ghosted_id = 1,grid%ngmax
    local_id = grid%nG2L(ghosted_id)
    if(local_id>0) then
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      if (xx_loc_p(istart) < 1.d-15) then
        temp = option%reference_temperature
      else
        ! T^{t+1,m} = (rho Cwi hw T)^{t+1} / rho^{t+1,m-1} Cw)^{t} (hw)^{t+1}
        do iter = 1,niter
          temp = xx_loc_p(iend)/xx_loc_p(istart)/ &
                  surf_global_auxvars(local_id)%den_kg(1)/ &
                  surf_auxvars(local_id)%Cwi - 273.15d0
          call EOSWaterdensity(temp,option%reference_pressure,den)
          surf_global_auxvars(local_id)%den_kg(1) = den
        enddo
      endif
      surf_global_auxvars(ghosted_id)%temp(1) = temp
    endif
  enddo

  ! Update boundary aux vars
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      surf_global_auxvars_bc(sum_connection)%temp(1) = &
        surf_global_auxvars(ghosted_id)%temp(1)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Update source/sink aux vars
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      surf_global_auxvars_ss(sum_connection)%temp(1) = &
        surf_global_auxvars(ghosted_id)%temp(1)

    enddo
    source_sink => source_sink%next
  enddo

  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr)

end subroutine SurfaceTHUpdateTemperature

! ************************************************************************** !

subroutine SurfaceTHUpdateSurfState(surf_realization)
  ! 
  ! This routine updates the states for surface-model at the end of
  ! subsurface-model timestep.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/25/13
  ! 

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
  use EOS_Water_module

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
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)

  PetscInt                :: count
  PetscInt                :: ghosted_id
  PetscInt                :: local_id
  PetscInt                :: ibeg
  PetscInt                :: iend
  PetscInt                :: iconn
  PetscInt                :: sum_connection

  PetscReal               :: den
  PetscReal, pointer      :: avg_vdarcy_p(:)   ! avg darcy velocity [m/s]
  PetscReal, pointer      :: xx_p(:)           ! head [m]
  PetscReal, pointer      :: surfpress_p(:)
  PetscReal, pointer      :: surftemp_p(:)
  PetscReal               :: Cwi
  PetscReal               :: temp_K
  PetscErrorCode          :: ierr

  PetscBool :: coupler_found = PETSC_FALSE

  patch      => surf_realization%patch
  option     => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_grid  => surf_realization%discretization%grid
  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars

  call VecGetArrayF90(surf_field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(surf_field%press_subsurf, surfpress_p, ierr)
  call VecGetArrayF90(surf_field%temp_subsurf, surftemp_p, ierr)

  count = 0
  do ghosted_id = 1,surf_grid%ngmax

    local_id = surf_grid%nG2L(ghosted_id)
    if(local_id <= 0) cycle

    iend = ghosted_id*option%nflowdof
    ibeg = iend - 1

    ! Compute density
    count = count + 1
    call EOSWaterdensity(surftemp_p(count),option%reference_pressure,den)
    xx_p(ibeg) = (surfpress_p(count)-option%reference_pressure)/ &
                        (abs(option%gravity(3)))/den
    if(xx_p(ibeg)<1.d-15) then
      xx_p(ibeg) = 0.d0
      xx_p(iend) = option%reference_temperature
    else
      Cwi = surf_auxvars(ghosted_id)%Cwi
      temp_K = surftemp_p(count) + 273.15d0
      xx_p(iend) = den*Cwi*temp_K*xx_p(ibeg)
    endif

  enddo
  call VecRestoreArrayF90(surf_field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf, surfpress_p, ierr)
  call VecRestoreArrayF90(surf_field%temp_subsurf, surftemp_p, ierr)

  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%flow_xx, &
                                   surf_field%flow_xx_loc, &
                                   NFLOWDOF)
  call SurfaceTHUpdateAuxVars(surf_realization)

end subroutine SurfaceTHUpdateSurfState

! ************************************************************************** !

subroutine SurfaceTHUpdateSolution(surf_realization)
  ! 
  ! This routine updates solution after a successful time step
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Surface_Realization_class
  use Surface_Field_module

  implicit none

  type(surface_realization_type)   :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode                   :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr)

end subroutine SurfaceTHUpdateSolution


end module Surface_TH_module

#endif
