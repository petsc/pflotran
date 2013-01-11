module Realization_module

  use Option_module
  use Output_Aux_module
  use Input_module
  use Region_module
  use Condition_module
  use Constraint_module
  use Material_module
  use Saturation_Function_module
  use Dataset_Aux_module
  use Fluid_module
  use Discretization_module
  use Field_module
  use Debug_module
  use Velocity_module
  use Waypoint_module
  use Output_Aux_module
  
  use Reaction_Aux_module
  
  use Level_module
  use Patch_module
  
  implicit none

private

#include "definitions.h"
  type, public :: realization_type

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    type(level_list_type), pointer :: level_list
    type(patch_type), pointer :: patch
      ! This is simply used to point to the patch associated with the current 
      ! level.  In other words: Never trust this!

    type(option_type), pointer :: option

    type(input_type), pointer :: input
    type(field_type), pointer :: field
    type(flow_debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option

    type(region_list_type), pointer :: regions
    type(condition_list_type), pointer :: flow_conditions
    type(tran_condition_list_type), pointer :: transport_conditions
    type(tran_constraint_list_type), pointer :: transport_constraints
    
    type(reaction_type), pointer :: reaction
    
    type(material_property_type), pointer :: material_properties
    type(material_property_ptr_type), pointer :: material_property_array(:)
    type(fluid_property_type), pointer :: fluid_properties
    type(fluid_property_type), pointer :: fluid_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    type(dataset_type), pointer :: datasets
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    
    type(velocity_dataset_type), pointer :: velocity_dataset
    
    type(waypoint_list_type), pointer :: waypoints
    
  end type realization_type

  interface RealizationCreate
    module procedure RealizationCreate1
    module procedure RealizationCreate2
  end interface
  
  public :: RealizationCreate, &
            RealizationDestroy, &
            RealizationProcessCouplers, &
            RealizationInitAllCouplerAuxVars, &
            RealizationProcessConditions, &
            RealizationUpdate, &
            RealizationAddWaypointsToList, &
            RealizationCreateDiscretization, &
            RealizationLocalizeRegions, &
            RealizationAddCoupler, &
            RealizationAddStrata, &
            RealizationAddObservation, &
            RealizUpdateUniformVelocity, &
            RealizationRevertFlowParameters, &
            RealizationGetDataset, &
            RealizGetDatasetValueAtCell, &
            RealizationSetDataset, &
            RealizationPrintCouplers, &
            RealizationInitConstraints, &
            RealProcessMatPropAndSatFunc, &
            RealProcessFluidProperties, &
            RealizationUpdateProperties, &
            RealizationCountCells, &
            RealizationPrintGridStatistics, &
            RealizationSetUpBC4Faces, &
            RealizatonPassPtrsToPatches, &
            RealLocalToLocalWithArray, &
            RealizationCalculateCFL1Timestep
 
contains
  
! ************************************************************************** !
!
! RealizationCreate1: Allocates and initializes a new Realization object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function RealizationCreate1()

  implicit none
  
  type(realization_type), pointer :: RealizationCreate1
  
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  nullify(option)
  RealizationCreate1 => RealizationCreate2(option)
  
end function RealizationCreate1  

! ************************************************************************** !
!
! RealizationCreate2: Allocates and initializes a new Realization object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function RealizationCreate2(option)

  implicit none
  
  type(option_type), pointer :: option
  
  type(realization_type), pointer :: RealizationCreate2
  
  type(realization_type), pointer :: realization
  
  allocate(realization)
  realization%discretization => DiscretizationCreate()
  if (associated(option)) then
    realization%option => option
  else
    realization%option => OptionCreate()
  endif
  nullify(realization%input)
  realization%field => FieldCreate()
  realization%debug => DebugCreateFlow()
  realization%output_option => OutputOptionCreate()

  realization%level_list => LevelCreateList()

  allocate(realization%regions)
  call RegionInitList(realization%regions)

  allocate(realization%flow_conditions)
  call FlowConditionInitList(realization%flow_conditions)
  allocate(realization%transport_conditions)
  call TranConditionInitList(realization%transport_conditions)
  allocate(realization%transport_constraints)
  call TranConstraintInitList(realization%transport_constraints)

  nullify(realization%material_properties)
  nullify(realization%material_property_array)
  nullify(realization%fluid_properties)
  nullify(realization%fluid_property_array)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  nullify(realization%datasets)
  nullify(realization%velocity_dataset)
  
  nullify(realization%reaction)

  nullify(realization%patch)
  nullify(realization%level_list)
  nullify(realization%waypoints)

  RealizationCreate2 => realization
  
end function RealizationCreate2 

! ************************************************************************** !
!
! RealizationCreateDiscretization: Creates grid
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationCreateDiscretization(realization)

  use Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Grid_module, only : UGridEnsureRightHandRule
  use Structured_Grid_module, only : StructGridCreateTVDGhosts
  use MFD_module
  use Coupler_module
  use Discretization_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  PetscErrorCode :: ierr
  PetscInt, allocatable :: int_tmp(:)
  PetscInt :: test,j, num_LP_dof
  PetscOffset :: i_da
  PetscReal, pointer :: real_tmp(:)
  type(dm_ptr_type), pointer :: dm_ptr
  Vec :: is_bnd_vec
  PetscInt :: ivar

  option => realization%option
  field => realization%field
 
  discretization => realization%discretization
  
  call DiscretizationCreateDMs(discretization,option)

  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity0, &
                                  GLOBAL,option)
 
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%tortuosity0)
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%volume)

  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%work)
  
  ! 1 degree of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity_loc, &
                                  LOCAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%tortuosity_loc)

  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%work_loc)
  
  if (option%nflowdof > 0) then

    ! 1-dof global  
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_xx)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_yy)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_zz)
    if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
      call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                         field%perm0_xz)
      call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                         field%perm0_xy)
      call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                         field%perm0_yz)
    endif
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm_pow)

    ! 1-dof local
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%ithrm_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%icap_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%iphas_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%iphas_old_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_xx_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_yy_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_zz_loc)
    if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
      call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                         field%perm_xz_loc)
      call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                         field%perm_xy_loc)
      call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                         field%perm_yz_loc)
    endif

    ! ndof degrees of freedom, global
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx, &
                                    GLOBAL,option)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_yy)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_dxx)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_r)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_accum)

    ! ndof degrees of freedom, local
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx_loc, &
                                    LOCAL,option)
  endif

  if (option%ntrandof > 0) then
    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
      ! ndof degrees of freedom, global
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_r)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_accum)

      ! ndof degrees of freedom, local
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                      LOCAL,option)
                                      
      if (realization%reaction%use_log_formulation) then
        call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                           field%tran_log_xx)
        call DiscretizationDuplicateVector(discretization,field%tran_xx_loc, &
                                           field%tran_work_loc)
      endif
 
    else ! operator splitting
      ! ndof degrees of freedom, global
      ! create the 1 dof vector for solving the individual linear systems
      call DiscretizationCreateVector(discretization,ONEDOF,field%tran_rhs_coef, &
                                      GLOBAL,option)
      ! create the ntran dof vector for storage of the solution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_rhs)

      ! ndof degrees of freedom, local
      ! again, just for storage of the current colution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                      LOCAL,option)

    endif
    
  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      grid => discretization%grid
      ! set up nG2L, nL2G, etc.
      call GridMapIndices(grid, discretization%dm_1dof%sgdm, &
                          discretization%stencil_type,&
                          discretization%lsm_flux_method, &
                          option)
      if (option%itranmode == EXPLICIT_ADVECTION) then
        call StructGridCreateTVDGhosts(grid%structured_grid, &
                                       realization%reaction%naqcomp, &
                                       field%tran_xx, &
                                       discretization%dm_1dof%sgdm, &
                                       field%tvd_ghosts, &
                                       discretization%tvd_ghost_scatter, &
                                       option)
      endif
      call GridComputeSpacing(grid,option)
      call GridComputeCoordinates(grid,discretization%origin,option)
      call GridComputeVolumes(grid,field%volume,option)
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option)
      if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
          call GridComputeCell2FaceConnectivity(grid, discretization%MFD, option)
      end if
      if (discretization%lsm_flux_method) then
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                          is_bnd_vec)
        call GridComputeNeighbors(grid,field%work_loc,option)
        call DiscretizationLocalToLocal(discretization,field%work_loc,is_bnd_vec,ONEDOF)
        call GridSaveBoundaryCellInfo(discretization%grid,is_bnd_vec,option)
        call VecDestroy(is_bnd_vec,ierr)
      endif
    case(UNSTRUCTURED_GRID)
      grid => discretization%grid
      ! set up nG2L, NL2G, etc.
      call UGridMapIndices(grid%unstructured_grid, &
                           discretization%dm_1dof%ugdm, &
                           grid%nG2L,grid%nL2G,grid%nG2A)
      call GridComputeCoordinates(grid,discretization%origin,option, & 
                                    discretization%dm_1dof%ugdm) 
      if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
        call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                      grid%y,grid%z,grid%nG2A,grid%nL2G, &
                                      option)
      endif
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option, &
                                      discretization%dm_1dof%ugdm) 
      call GridComputeVolumes(grid,field%volume,option)
  end select 
 
  ! Vectors with face degrees of freedom
#ifdef DASVYAT
   if (discretization%itype == STRUCTURED_GRID_MIMETIC) then

     if (option%nflowdof > 0) then

       num_LP_dof = (grid%nlmax_faces + grid%nlmax)*option%nflowdof
   
       call VecCreateMPI(option%mycomm, num_LP_dof, &
                    PETSC_DETERMINE,field%flow_xx_faces,ierr)



       call VecSetBlockSize(field%flow_xx_faces,option%nflowdof,ierr)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_r_faces)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_dxx_faces)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_yy_faces)


       call VecCreateSeq(PETSC_COMM_SELF, (grid%ngmax_faces + grid%ngmax)*option%nflowdof, field%flow_xx_loc_faces, ierr)
       call VecSetBlockSize(field%flow_xx_loc_faces,option%nflowdof,ierr)

!       call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*option%nflowdof, field%flow_r_loc_faces, ierr)
!       call VecSetBlockSize(field%flow_r_loc_faces,NFLOWDOF,ierr)

!       call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*option%nflowdof, field%flow_bc_loc_faces, ierr)
!       call VecSetBlockSize(field%flow_bc_loc_faces,NFLOWDOF,ierr)

        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%flow_r_loc_faces) 

        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%flow_bc_loc_faces)
   
        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%work_loc_faces)

!       call VecGetArrayF90(field%volume, real_tmp, ierr)
!       call VecRestoreArrayF90(field%volume, real_tmp, ierr)
 


     end if

     call RealizationCreatenG2LP(realization)


     dm_ptr => DiscretizationGetDMPtrFromIndex(discretization, NFLOWDOF)
     call GridComputeGlobalCell2FaceConnectivity(grid, discretization%MFD, dm_ptr%sgdm, NFLOWDOF, option)
  
   end if
#endif
 
  ! initialize to -999.d0 for check later that verifies all values 
  ! have been set
  call VecSet(field%porosity0,-999.d0,ierr)

  ! Allocate vectors to hold temporally average output quantites
  if(realization%output_option%aveg_output_variable_list%nvars>0) then

    field%nvars = realization%output_option%aveg_output_variable_list%nvars
    allocate(field%avg_vars_vec(field%nvars))

    do ivar=1,field%nvars
      call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                         field%avg_vars_vec(ivar))
      call VecSet(field%avg_vars_vec(ivar),0.d0,ierr)
    enddo
  endif
       

end subroutine RealizationCreateDiscretization

! ************************************************************************** !
!
! RealizationLocalizeRegions: Localizes regions within each patch
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationLocalizeRegions(realization)

  use Option_module
  use String_module
  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type (region_type), pointer :: cur_region, cur_region2
  type(option_type), pointer :: option

  option => realization%option

  ! check to ensure that region names are not duplicated
  cur_region => realization%regions%first
  do
    if (.not.associated(cur_region)) exit
    cur_region2 => cur_region%next
    do
      if (.not.associated(cur_region2)) exit
      if (StringCompare(cur_region%name,cur_region2%name,MAXWORDLENGTH)) then
        option%io_buffer = 'Duplicate region names: ' // trim(cur_region%name)
        call printErrMsg(option)
      endif
      cur_region2 => cur_region2%next
    enddo
    cur_region => cur_region%next
  enddo

  call PatchLocalizeRegions(realization%patch,realization%regions, &
                            realization%option)

end subroutine RealizationLocalizeRegions

! ************************************************************************** !
!
! RealizatonPassPtrsToPatches: Sets patch%field => realization%field
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine RealizatonPassPtrsToPatches(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  realization%patch%field => realization%field
  realization%patch%datasets => realization%datasets
  realization%patch%reaction => realization%reaction
  
end subroutine RealizatonPassPtrsToPatches

! ************************************************************************** !
!
! RealizationAddCoupler: Adds a copy of a coupler to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddCoupler(realization,coupler)

  use Coupler_module

  implicit none
  
  type(realization_type) :: realization
  type(coupler_type), pointer :: coupler
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: patch
  
  type(coupler_type), pointer :: new_coupler
  
  patch => realization%patch
  
  ! only add to flow list for now, since they will be split out later
  new_coupler => CouplerCreate(coupler)
  select case(coupler%itype)
    case(BOUNDARY_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%boundary_conditions)
    case(INITIAL_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%initial_conditions)
    case(SRC_SINK_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%source_sinks)
  end select
  nullify(new_coupler)

  call CouplerDestroy(coupler)
 
end subroutine RealizationAddCoupler

subroutine RealizationCreatenG2LP(realization)

    use Grid_module

    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"


#include "definitions.h"

#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"

   
    type(realization_type) :: realization  


    type(option_type), pointer :: option
    type(discretization_type), pointer :: discretization 
    type(grid_type), pointer :: grid
    PetscInt :: DOF, global_offset, icell, num_ghosted
    PetscErrorCode :: ierr

    Vec :: vec_LP_cell_id
    Vec :: vec_LP_cell_id_loc

    IS :: is_ghosted, is_global
    VecScatter :: VC_global2ghosted 

    PetscScalar, pointer :: lp_cell_ids(:), lp_cell_ids_loc(:) 
    PetscInt, pointer :: int_tmp_gh(:), int_tmp_gl(:)

    option => realization%option
    discretization => realization%discretization
    grid => discretization%grid
  
    global_offset = 0
    grid%global_faces_offset = 0
    grid%global_cell_offset = 0

    allocate(grid%nG2LP(grid%ngmax))

    call MPI_Exscan(grid%nlmax_faces, grid%global_faces_offset, &
                      ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)

    call MPI_Exscan(grid%nlmax, grid%global_cell_offset, &
                      ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)


    global_offset = grid%global_faces_offset + grid%global_cell_offset

    call DiscretizationCreateVector(discretization,ONEDOF,vec_LP_cell_id, &
                                      GLOBAL,option)

    call DiscretizationCreateVector(discretization,ONEDOF,vec_LP_cell_id_loc, &
                                      LOCAL,option)

    call VecGetArrayF90(vec_LP_cell_id, lp_cell_ids, ierr)     

     do icell = 1, grid%nlmax
        grid%nG2LP(grid%nL2G(icell)) = global_offset + grid%nlmax_faces + icell - 1
        lp_cell_ids(icell) = global_offset + grid%nlmax_faces + icell
     end do
          
    call VecRestoreArrayF90(vec_LP_cell_id, lp_cell_ids, ierr)     

    allocate(int_tmp_gh(grid%ngmax - grid%nlmax))
    allocate(int_tmp_gl(grid%ngmax - grid%nlmax))

    num_ghosted = 1
    do icell = 1, grid%ngmax
       if (grid%nG2L(icell) < 1) then
          int_tmp_gh(num_ghosted) = icell - 1
          int_tmp_gl(num_ghosted) = grid%nG2P(icell)
          num_ghosted = num_ghosted + 1
       end if
    end do


    call ISCreateBlock(option%mycomm, ONEDOF, grid%ngmax - grid%nlmax, &
                     int_tmp_gh, PETSC_COPY_VALUES, is_ghosted, ierr)
    call ISCreateBlock(option%mycomm, ONEDOF, grid%ngmax - grid%nlmax, &
                     int_tmp_gl, PETSC_COPY_VALUES, is_global, ierr)



    call VecScatterCreate(vec_LP_cell_id, is_global, vec_LP_cell_id_loc, is_ghosted, VC_global2ghosted, ierr) 

    deallocate(int_tmp_gh)
    deallocate(int_tmp_gl)

    call VecScatterBegin(VC_global2ghosted, vec_LP_cell_id, vec_LP_cell_id_loc, &
                              INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(VC_global2ghosted, vec_LP_cell_id, vec_LP_cell_id_loc, &
                              INSERT_VALUES,SCATTER_FORWARD,ierr)
 

!    call DiscretizationGlobalToLocal(discretization, vec_LP_cell_id, vec_LP_cell_id_loc, ONEDOF)

    call VecGetArrayF90(vec_LP_cell_id_loc, lp_cell_ids_loc, ierr)


     do icell = 1, grid%ngmax
      if (grid%nG2L(icell) < 1) grid%nG2LP(icell) = lp_cell_ids_loc(icell) - 1
    end do

    call VecRestoreArrayF90(vec_LP_cell_id_loc, lp_cell_ids_loc, ierr)


    call VecDestroy(vec_LP_cell_id, ierr)
    call VecDestroy(vec_LP_cell_id_loc, ierr)
    call VecScatterDestroy(VC_global2ghosted , ierr)
    call ISDestroy(is_ghosted, ierr)
    call ISDestroy(is_global, ierr)
!    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

end subroutine


! ************************************************************************** !
!
! RealizationAddStrata: Adds a copy of a strata to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddStrata(realization,strata)

  use Strata_module

  implicit none
  
  type(realization_type) :: realization
  type(strata_type), pointer :: strata
  
  type(strata_type), pointer :: new_strata
  
  new_strata => StrataCreate(strata)
  call StrataAddToList(new_strata,realization%patch%strata)
  nullify(new_strata)
  
  call StrataDestroy(strata)
 
end subroutine RealizationAddStrata

! ************************************************************************** !
!
! RealizationAddObservation: Adds a copy of a observation object to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddObservation(realization,observation)

  use Observation_module

  implicit none
  
  type(realization_type) :: realization
  type(observation_type), pointer :: observation
  
  type(observation_type), pointer :: new_observation
  
  new_observation => ObservationCreate(observation)
  call ObservationAddToList(new_observation, &
                            realization%patch%observation)
  nullify(new_observation)

  call ObservationDestroy(observation)
 
end subroutine RealizationAddObservation

! ************************************************************************** !
!
! RealizationProcessCouplers: Sets connectivity and pointers for couplers
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationProcessCouplers(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  call PatchProcessCouplers( realization%patch,realization%flow_conditions, &
                             realization%transport_conditions, &
                             realization%option)
  
end subroutine RealizationProcessCouplers

! ************************************************************************** !
!
! RealizationProcessConditions: Sets up auxiliary data associated with 
!                               conditions
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine RealizationProcessConditions(realization)

  use Dataset_module
  
  implicit none
  
  type(realization_type) :: realization

  call DatasetProcessDatasets(realization%datasets,realization%option)
  
  if (realization%option%nflowdof > 0) then
    call RealProcessFlowConditions(realization)
  endif
  if (realization%option%ntrandof > 0) then
    call RealProcessTranConditions(realization)
  endif
 
end subroutine RealizationProcessConditions

! ************************************************************************** !
!
! RealProcessMatPropAndSatFunc: Sets up linkeage between material properties
!                               and saturation function, auxiliary arrays
!                               and datasets
! author: Glenn Hammond
! date: 01/21/09, 01/12/11
!
! ************************************************************************** !
subroutine RealProcessMatPropAndSatFunc(realization)

  use String_module
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  type(patch_type), pointer :: patch
  character(len=MAXSTRINGLENGTH) :: string

  option => realization%option
  patch => realization%patch
  
  ! organize lists
  call MaterialPropConvertListToArray(realization%material_properties, &
                                      realization%material_property_array, &
                                      option)
  call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                      realization%saturation_function_array, &
                                      option)

  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  patch%material_properties => realization%material_properties
  call MaterialPropConvertListToArray(patch%material_properties, &
                                      patch%material_property_array, &
                                      option)
  patch%saturation_functions => realization%saturation_functions
  call SaturatFuncConvertListToArray(patch%saturation_functions, &
                                      patch%saturation_function_array, &
                                      option)
    
  cur_material_property => realization%material_properties                            
  do                                      
    if (.not.associated(cur_material_property)) exit

    ! obtain saturation function id
    if (option%iflowmode /= NULL_MODE) then
      cur_material_property%saturation_function_id = &
        SaturationFunctionGetID(realization%saturation_functions, &
                                cur_material_property%saturation_function_name, &
                                cur_material_property%name,option)
    endif
    
    ! if named, link dataset to property
    if (.not.StringNull(cur_material_property%porosity_dataset_name)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),POROSITY'
      cur_material_property%porosity_dataset => &
        DatasetGetPointer(realization%datasets, &
                          cur_material_property%porosity_dataset_name, &
                          string,option)
    endif
    if (.not.StringNull(cur_material_property%permeability_dataset_name)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY'
      cur_material_property%permeability_dataset => &
        DatasetGetPointer(realization%datasets, &
                          cur_material_property%permeability_dataset_name, &
                          string,option)
    endif
    
    cur_material_property => cur_material_property%next
  enddo
  
  
end subroutine RealProcessMatPropAndSatFunc

! ************************************************************************** !
!
! RealProcessFluidProperties: Sets up linkeage with fluid properties
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
subroutine RealProcessFluidProperties(realization)
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscBool :: found
  type(option_type), pointer :: option
  type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  
  found = PETSC_FALSE
  cur_fluid_property => realization%fluid_properties                            
  do                                      
    if (.not.associated(cur_fluid_property)) exit
    found = PETSC_TRUE
    select case(trim(cur_fluid_property%phase_name))
      case('LIQUID_PHASE')
        cur_fluid_property%phase_id = LIQUID_PHASE
      case('GAS_PHASE')
        cur_fluid_property%phase_id = GAS_PHASE
      case default
        cur_fluid_property%phase_id = LIQUID_PHASE
    end select
    cur_fluid_property => cur_fluid_property%next
  enddo
  
  if (option%ntrandof > 0 .and. .not.found) then
    option%io_buffer = 'A fluid property must be present in input file' // &
                       ' for solute transport'
  endif
  
end subroutine RealProcessFluidProperties

! ************************************************************************** !
!
! RealProcessFlowConditions: Sets linkage of flow conditions to dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine RealProcessFlowConditions(realization)

  use Dataset_module

  implicit none

  type(realization_type) :: realization
  
  type(flow_condition_type), pointer :: cur_flow_condition
  type(flow_sub_condition_type), pointer :: cur_flow_sub_condition
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: i
  type(dataset_type), pointer :: dataset
  
  option => realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_flow_condition => realization%flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    !TODO(geh): could destroy the time_series here if dataset allocated
    select case(option%iflowmode)
      case(G_MODE)
      case(RICHARDS_MODE,MIS_MODE)
        do i = 1, size(cur_flow_condition%sub_condition_ptr)
          ! check for dataset in flow_dataset
          if (associated(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          flow_dataset%dataset)) then
            dataset_name = cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          flow_dataset%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                                flow_dataset%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_flow_condition%name)
            dataset => &
              DatasetGetPointer(realization%datasets,dataset_name,string,option)
            cur_flow_condition%sub_condition_ptr(i)%ptr%flow_dataset%dataset => &
              dataset
            nullify(dataset)
          endif
          if (associated(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          datum%dataset)) then
            dataset_name = cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          datum%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                                datum%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_flow_condition%name)
            dataset => &
              DatasetGetPointer(realization%datasets,dataset_name,string,option)
            cur_flow_condition%sub_condition_ptr(i)%ptr%datum%dataset => &
              dataset
            nullify(dataset)
          endif
          if (associated(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          gradient%dataset)) then
            dataset_name = cur_flow_condition%sub_condition_ptr(i)%ptr% &
                          gradient%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_flow_condition%sub_condition_ptr(i)%ptr% &
                                gradient%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_flow_condition%name)
            dataset => &
              DatasetGetPointer(realization%datasets,dataset_name,string,option)
            cur_flow_condition%sub_condition_ptr(i)%ptr%gradient%dataset => &
              dataset
            nullify(dataset)
          endif
        enddo
    end select
    cur_flow_condition => cur_flow_condition%next
  enddo

end subroutine RealProcessFlowConditions

! ************************************************************************** !
!
! RealProcessTranConditions: Sets up auxiliary data associated with 
!                            transport conditions
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine RealProcessTranConditions(realization)

  use String_module
  use Reaction_module
  use Constraint_module
  
  implicit none
  
  type(realization_type) :: realization
  
  
  PetscBool :: found
  type(option_type), pointer :: option
  type(tran_condition_type), pointer :: cur_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(tran_constraint_type), pointer :: cur_constraint, another_constraint
  
  option => realization%option
  
  ! check for duplicate constraint names
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
      another_constraint => cur_constraint%next
      ! now compare names
      found = PETSC_FALSE
      do
        if (.not.associated(another_constraint)) exit
        if (StringCompare(cur_constraint%name,another_constraint%name, &
            MAXWORDLENGTH)) then
          found = PETSC_TRUE
        endif
        another_constraint => another_constraint%next
      enddo
      if (found) then
        option%io_buffer = 'Duplicate transport constraints named "' // &
                 trim(cur_constraint%name) // '"'
        call printErrMsg(realization%option)
      endif
    cur_constraint => cur_constraint%next
  enddo
  
  ! initialize constraints
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
    call ReactionProcessConstraint(realization%reaction, &
                                   cur_constraint%name, &
                                   cur_constraint%aqueous_species, &
                                   cur_constraint%minerals, &
                                   cur_constraint%surface_complexes, &
                                   cur_constraint%colloids, &
                                   realization%option)
    cur_constraint => cur_constraint%next
  enddo
  
  ! tie constraints to couplers, if not already associated
  cur_condition => realization%transport_conditions%first
  do

    if (.not.associated(cur_condition)) exit
    cur_constraint_coupler => cur_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      if (.not.associated(cur_constraint_coupler%aqueous_species)) then
        cur_constraint => realization%transport_constraints%first
        do
          if (.not.associated(cur_constraint)) exit
          if (StringCompare(cur_constraint%name, &
                             cur_constraint_coupler%constraint_name, &
                             MAXWORDLENGTH)) then
            cur_constraint_coupler%aqueous_species => cur_constraint%aqueous_species
            cur_constraint_coupler%minerals => cur_constraint%minerals
            cur_constraint_coupler%surface_complexes => cur_constraint%surface_complexes
            cur_constraint_coupler%colloids => cur_constraint%colloids
            exit
          endif
          cur_constraint => cur_constraint%next
        enddo
        if (.not.associated(cur_constraint_coupler%aqueous_species)) then
          option%io_buffer = 'Transport constraint "' // &
                   trim(cur_constraint_coupler%constraint_name) // &
                   '" not found in input file constraints.'
          call printErrMsg(realization%option)
        endif
      endif
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    if (associated(cur_condition%constraint_coupler_list%next)) then ! more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    cur_condition => cur_condition%next
  enddo
 
  ! final details for setup
  cur_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_condition)) exit
    ! is the condition transient?
    if (associated(cur_condition%constraint_coupler_list%next)) then ! more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    ! set pointer to first constraint coupler
    cur_condition%cur_constraint_coupler => cur_condition%constraint_coupler_list
    
    cur_condition => cur_condition%next
  enddo

end subroutine RealProcessTranConditions

! ************************************************************************** !
!
! RealizationInitConstraints: Initializes constraint concentrations
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine RealizationInitConstraints(realization)

  implicit none

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchInitConstraints(cur_patch,realization%reaction, &
                                realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo            
 
end subroutine RealizationInitConstraints

! ************************************************************************** !
!
! RealizationPrintCouplers: Print boundary and initial condition data
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine RealizationPrintCouplers(realization)

  use Coupler_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(coupler_type), pointer :: cur_coupler
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
 
  option => realization%option
  reaction => realization%reaction
 
  if (.not.OptionPrintToFile(option)) return
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      cur_coupler => cur_patch%initial_conditions%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo
     
      cur_coupler => cur_patch%boundary_conditions%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo
     
      cur_coupler => cur_patch%source_sinks%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo            
 
end subroutine RealizationPrintCouplers

! ************************************************************************** !
!
! RealizationPrintCoupler: Prints boundary and initial condition coupler 
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine RealizationPrintCoupler(coupler,reaction,option)

  use Coupler_module
  use Reaction_module
  
  implicit none
  
  type(coupler_type) :: coupler
  type(option_type) :: option
  type(reaction_type), pointer :: reaction
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(region_type), pointer :: region
  type(tran_constraint_coupler_type), pointer :: constraint_coupler
   
98 format(40('=+'))
99 format(80('-'))
100 format(a)
  
  flow_condition => coupler%flow_condition
  tran_condition => coupler%tran_condition
  region => coupler%region

  write(option%fid_out,*)
  write(option%fid_out,98)


  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      string = 'Initial Condition'
    case(BOUNDARY_COUPLER_TYPE)
      string = 'Boundary Condition'
    case(SRC_SINK_COUPLER_TYPE)
      string = 'Source Sink'
  end select
  write(option%fid_out,'(/,2x,a,/)') trim(string)

  write(option%fid_out,99)
101 format(5x,'     Flow Condition: ',2x,a)
  if (associated(flow_condition)) write(option%fid_out,101) trim(flow_condition%name)
102 format(5x,'Transport Condition: ',2x,a)
  if (associated(tran_condition)) write(option%fid_out,102) trim(tran_condition%name)
103 format(5x,'             Region: ',2x,a)
  if (associated(region)) write(option%fid_out,103) trim(region%name)
  write(option%fid_out,99)
  
  if (associated(flow_condition)) then
    call FlowConditionPrint(flow_condition,option)
  endif
  if (associated(tran_condition)) then
    constraint_coupler => tran_condition%cur_constraint_coupler
    write(option%fid_out,'(/,2x,''Transport Condition: '',a)') &
      trim(tran_condition%name)
    if (associated(reaction)) then
      call ReactionPrintConstraint(constraint_coupler,reaction,option)
      write(option%fid_out,'(/)')
      write(option%fid_out,99)
    endif
  endif
 
end subroutine RealizationPrintCoupler

! ************************************************************************** !
!
! RealizationInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationInitAllCouplerAuxVars(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  !geh: Must update conditions prior to initializing the aux vars.  
  !     Otherwise, datasets will not have been read for routines such as
  !     hydrostatic and auxvars will be initialized to garbage.
  call FlowConditionUpdate(realization%flow_conditions,realization%option, &
                           realization%option%time)
  call TranConditionUpdate(realization%transport_conditions, &
                           realization%option, &
                           realization%option%time)  
  call PatchInitAllCouplerAuxVars(realization%patch,realization%option)
   
end subroutine RealizationInitAllCouplerAuxVars

! ************************************************************************** !
!
! RealizUpdateAllCouplerAuxVars: Updates auxiliary variables associated 
!                                  with couplers in lis
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizUpdateAllCouplerAuxVars(realization,force_update_flag)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  PetscBool :: force_update_flag

  call PatchUpdateAllCouplerAuxVars(realization%patch,force_update_flag, &
                                    realization%option)

end subroutine RealizUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! RealizationUpdate: Update parameters in realization (e.g. conditions, bcs, srcs)
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine RealizationUpdate(realization)

  implicit none
  
  type(realization_type) :: realization
  
  PetscBool :: force_update_flag = PETSC_FALSE
  
  ! must update conditions first
  call FlowConditionUpdate(realization%flow_conditions,realization%option, &
                           realization%option%time)
  call TranConditionUpdate(realization%transport_conditions, &
                           realization%option, &
                           realization%option%time)
  call RealizUpdateAllCouplerAuxVars(realization,force_update_flag)
  if (associated(realization%velocity_dataset)) then
    call VelocityDatasetUpdate(realization%option,realization%option%time, &
                               realization%velocity_dataset)
    call RealizUpdateUniformVelocity(realization)
  endif
! currently don't use aux_vars, just condition for src/sinks
!  call RealizationUpdateSrcSinks(realization)

end subroutine RealizationUpdate

! ************************************************************************** !
!
! RealizationRevertFlowParameters: Assigns initial porosity/perms to vecs
! author: Glenn Hammond
! date: 05/09/08
!
! ************************************************************************** !
subroutine RealizationRevertFlowParameters(realization)

  use Option_module
  use Field_module
  use Discretization_module

  implicit none
  
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                           field%perm_xx_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                           field%perm_yy_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                           field%perm_zz_loc,ONEDOF)   
  endif   
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%porosity_loc,ONEDOF)
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                   field%tortuosity_loc,ONEDOF)
                           
end subroutine RealizationRevertFlowParameters

! ************************************************************************** !
!
! RealizUpdateUniformVelocity: Assigns uniform velocity for transport
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizUpdateUniformVelocity(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  call PatchUpdateUniformVelocity(realization%patch, &
                                  realization%velocity_dataset%cur_value, &
                                      realization%option)
 
end subroutine RealizUpdateUniformVelocity

! ************************************************************************** !
!
! RealizationAddWaypointsToList: Creates waypoints associated with source/sinks
!                             boundary conditions, etc. and add to list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RealizationAddWaypointsToList(realization)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(realization_type) :: realization
  
  type(waypoint_list_type), pointer :: waypoint_list
  type(flow_condition_type), pointer :: cur_flow_condition
  type(tran_condition_type), pointer :: cur_tran_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => realization%option
  waypoint_list => realization%waypoints
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_output = realization%output_option%print_final
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in RealizationAddWaypointsToList'
    call printErrMsg(option)
  endif

  ! add update of flow conditions
  cur_flow_condition => realization%flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    if (cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        !TODO(geh): check if this updated more than simply the flow_dataset (i.e. datum and gradient)
        !geh: followup - no, datum/gradient are not considered.  Should they be considered?
        call FlowConditionDatasetGetTimes(option,sub_condition,final_time, &
                                          times)
        if (size(times) > 1000) then
          option%io_buffer = 'For flow condition "' // &
            trim(cur_flow_condition%name) // &
            '" dataset "' // trim(sub_condition%name) // &
            '", the number of times is excessive for synchronization ' // &
            'with waypoints.'
          call printErrMsg(option)
        endif
        do itime = 1, size(times)
          waypoint => WaypointCreate()
          waypoint%time = times(itime)
          waypoint%update_conditions = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        enddo
        deallocate(times)
        nullify(times)
      enddo
    endif
    cur_flow_condition => cur_flow_condition%next
  enddo
      
  ! add update of transport conditions
  cur_tran_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_tran_condition)) exit
    if (cur_tran_condition%is_transient) then
      cur_constraint_coupler => cur_tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        if (cur_constraint_coupler%time > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = cur_constraint_coupler%time
          waypoint%update_conditions = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
    cur_tran_condition => cur_tran_condition%next
  enddo

  ! add update of velocity fields
  if (associated(realization%velocity_dataset)) then
    if (realization%velocity_dataset%times(1) > 1.d-40 .or. &
        size(realization%velocity_dataset%times) > 1) then
      do itime = 1, size(realization%velocity_dataset%times)
        waypoint => WaypointCreate()
        waypoint%time = realization%velocity_dataset%times(itime)
        waypoint%update_conditions = PETSC_TRUE
        call WaypointInsertInList(waypoint,waypoint_list)
      enddo
    endif
  endif

  ! add waypoints for periodic output
  if (realization%output_option%periodic_output_time_incr > 0.d0 .or. &
      realization%output_option%periodic_tr_output_time_incr > 0.d0) then

    if (realization%output_option%periodic_output_time_incr > 0.d0) then
      ! standard output
      temp_real = 0.d0
      do
        temp_real = temp_real + realization%output_option%periodic_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,realization%waypoints)
      enddo
    endif
    
    if (realization%output_option%periodic_tr_output_time_incr > 0.d0) then
      ! transient observation output
      temp_real = 0.d0
      do
        temp_real = temp_real + realization%output_option%periodic_tr_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_tr_output = PETSC_TRUE 
        call WaypointInsertInList(waypoint,realization%waypoints)
      enddo
    endif

  endif

end subroutine RealizationAddWaypointsToList

! ************************************************************************** !
!
! RealizationGetDataset: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationGetDataset(realization,vec,ivar,isubvar,isubvar1)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  call PatchGetDataset(realization%patch,realization%field, &
                       realization%reaction,realization%option, &
                       realization%output_option,vec,ivar,isubvar,isubvar1)

end subroutine RealizationGetDataset

! ************************************************************************** !
!
! RealizGetDatasetValueAtCell: Extracts variables indexed by ivar and isubvar
!                              from a realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id, &
                                     isubvar1)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetDatasetValueAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: ghosted_id
  
  PetscReal :: value
  
  value = PatchGetDatasetValueAtCell(realization%patch,realization%field, &
                                     realization%reaction, &
                                     realization%option, &
                                     realization%output_option, &
                                     ivar,isubvar,ghosted_id,isubvar1)
  RealizGetDatasetValueAtCell = value

end function RealizGetDatasetValueAtCell

! ************************************************************************** !
!
! RealizationSetDataset: Sets variables indexed by ivar and isubvar in a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationSetDataset(realization,vec,vec_format,ivar,isubvar)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetDataset(realization%patch,realization%field, &
                       realization%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetDataset

! ************************************************************************** !
!
! RealizationUpdateProperties: Updates coupled properties at each grid cell
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RealizationUpdateProperties(realization)

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal :: min_value  
  PetscInt :: ivalue
  PetscErrorCode :: ierr
  
  option => realization%option
    
  call RealizationUpdatePropertiesPatch(realization)
  
  ! perform check to ensure that porosity is bounded between 0 and 1
  ! since it is calculated as 1.d-sum_volfrac, it cannot be > 1
  call VecMin(realization%field%porosity_loc,ivalue,min_value,ierr)
  if (min_value < 0.d0) then
    write(option%io_buffer,*) 'Sum of mineral volume fractions has ' // &
      'exceeded 1.d0 at cell (note PETSc numbering): ', ivalue
    call printErrMsg(option)
  endif
   
end subroutine RealizationUpdateProperties

! ************************************************************************** !
!
! RealizationUpdatePropertiesPatch: Updates coupled properties at each grid cell 
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RealizationUpdatePropertiesPatch(realization)

  use Grid_module
  use Reactive_Transport_Aux_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:) 
  type(discretization_type), pointer :: discretization

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl
  PetscReal :: sum_volfrac
  PetscReal :: scale, porosity_scale, volfrac_scale
  PetscBool :: porosity_updated
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: porosity_loc_p(:), porosity0_p(:)
  PetscReal, pointer :: tortuosity_loc_p(:), tortuosity0_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  reaction => realization%reaction
  grid => patch%grid
  material_property_array => realization%material_property_array
  rt_auxvars => patch%aux%RT%aux_vars

  if (.not.associated(patch%imat)) then
    option%io_buffer = 'Materials IDs not present in run.  Material ' // &
      ' properties cannot be updated without material ids ask Glenn'
    call printErrMsg(option)
  endif

  porosity_updated = PETSC_FALSE
  if (reaction%update_porosity) then
    porosity_updated = PETSC_TRUE
  
    call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

    if (reaction%mineral%nkinmnrl > 0) then
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)

        ! Go ahead and compute for inactive cells since their porosity does
        ! not matter (avoid check on active/inactive)
        sum_volfrac = 0.d0
        do imnrl = 1, reaction%mineral%nkinmnrl
          sum_volfrac = sum_volfrac + &
                        rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)
        enddo 
        porosity_loc_p(ghosted_id) = max(1.d0-sum_volfrac, &
                                         reaction%minimum_porosity)
      enddo
    endif

    call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    
  endif
  
  if ((porosity_updated .and. &
       (reaction%update_tortuosity .or. &
        reaction%update_permeability)) .or. &
      ! if porosity ratio is used in mineral surface area update, we must
      ! recalculate it every time.
      (reaction%update_mineral_surface_area .and. &
       reaction%update_mnrl_surf_with_porosity)) then
    call GridVecGetArrayF90(grid,field%porosity0,porosity0_p,ierr)
    call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = porosity_loc_p(ghosted_id)/porosity0_p(local_id)
    enddo
    call GridVecRestoreArrayF90(grid,field%porosity0,porosity0_p,ierr)
    call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
  endif      

  if (reaction%update_mineral_surface_area) then
    porosity_scale = 1.d0
!   if (option%update_mnrl_surf_with_porosity) then
    if (reaction%update_mnrl_surf_with_porosity) then
      ! placing the get/restore array calls within the condition will
      ! avoid improper access.
      call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    endif
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (reaction%update_mnrl_surf_with_porosity) then
        porosity_scale = vec_p(local_id)** &
             reaction%mineral%kinmnrl_surf_area_porosity_pwr(imnrl)
!geh: srf_area_vol_frac_pwr must be defined on a per mineral basis, not
!     solely material type.
!          material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_porosity_pwr
      endif
      do imnrl = 1, reaction%mineral%nkinmnrl
        if (rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) > 0.d0) then
          volfrac_scale = (rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)/ &
                         rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl))** &
             reaction%mineral%kinmnrl_surf_area_vol_frac_pwr(imnrl)
!geh: srf_area_vol_frac_pwr must be defined on a per mineral basis, not
!     solely material type.
!            material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_volfrac_pwr
          rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
            rt_auxvars(ghosted_id)%mnrl_area0(imnrl)*porosity_scale*volfrac_scale
        else
          rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
            rt_auxvars(ghosted_id)%mnrl_area0(imnrl)
        endif
      enddo
    enddo
    if (reaction%update_mnrl_surf_with_porosity) then
      call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
    endif

    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                     field%tortuosity_loc,ONEDOF)
  endif
      
  if (reaction%update_tortuosity) then
    call GridVecGetArrayF90(grid,field%tortuosity_loc,tortuosity_loc_p,ierr)  
    call GridVecGetArrayF90(grid,field%tortuosity0,tortuosity0_p,ierr)  
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = vec_p(local_id)** &
        material_property_array(patch%imat(ghosted_id))%ptr%tortuosity_pwr
      tortuosity_loc_p(ghosted_id) = tortuosity0_p(local_id)*scale
    enddo
    call GridVecRestoreArrayF90(grid,field%tortuosity_loc,tortuosity_loc_p,ierr)  
    call GridVecRestoreArrayF90(grid,field%tortuosity0,tortuosity0_p,ierr)  
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)

    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                     field%tortuosity_loc,ONEDOF)
  endif
      
  if (reaction%update_permeability) then
    call GridVecGetArrayF90(grid,field%perm0_xx,perm0_xx_p,ierr)
    call GridVecGetArrayF90(grid,field%perm0_zz,perm0_zz_p,ierr)
    call GridVecGetArrayF90(grid,field%perm0_yy,perm0_yy_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_zz_loc,perm_zz_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_yy_loc,perm_yy_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = vec_p(local_id)** &
              material_property_array(patch%imat(ghosted_id))%ptr%permeability_pwr
      perm_xx_loc_p(ghosted_id) = perm0_xx_p(local_id)*scale
      perm_yy_loc_p(ghosted_id) = perm0_yy_p(local_id)*scale
      perm_zz_loc_p(ghosted_id) = perm0_zz_p(local_id)*scale
    enddo
    call GridVecRestoreArrayF90(grid,field%perm0_xx,perm0_xx_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm0_zz,perm0_zz_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm0_yy,perm0_yy_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_zz_loc,perm_zz_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_yy_loc,perm_yy_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)

    call DiscretizationLocalToLocal(discretization,field%perm_xx_loc, &
                                    field%perm_xx_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%perm_yy_loc, &
                                    field%perm_yy_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%perm_zz_loc, &
                                    field%perm_zz_loc,ONEDOF)
  endif  
  
end subroutine RealizationUpdatePropertiesPatch

! ************************************************************************** !
!
! RealLocalToLocalWithArray: Takes an F90 array that is ghosted
!                            and updates the ghosted values
! author: Glenn Hammond
! date: 06/09/11
!
! ************************************************************************** !
subroutine RealLocalToLocalWithArray(realization,array_id)

  use Grid_module

  implicit none

  type(realization_type) :: realization
  PetscInt :: array_id
  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field

  field => realization%field
  patch => realization%patch

  grid => patch%grid
  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%imat,field%work_loc, &
                                     grid%ngmax)
    case(SATURATION_FUNCTION_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%sat_func_id, &
                                     field%work_loc, grid%ngmax)
  end select

  call DiscretizationLocalToLocal(realization%discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)

  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%imat,field%work_loc, &
                                      grid%ngmax)
    case(SATURATION_FUNCTION_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%sat_func_id, &
                                      field%work_loc, grid%ngmax)
  end select

end subroutine RealLocalToLocalWithArray

! ************************************************************************** !
!
! RealizationCountCells: Counts # of active and inactive grid cells 
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine RealizationCountCells(realization,global_total_count, &
                                 global_active_count,total_count,active_count)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  PetscInt :: global_total_count
  PetscInt :: global_active_count
  PetscInt :: total_count
  PetscInt :: active_count
  
  PetscInt :: patch_total_count
  PetscInt :: patch_active_count
  PetscInt :: temp_int_in(2), temp_int_out(2)
  PetscErrorCode :: ierr
  
  type(patch_type), pointer :: patch
  
  total_count = 0
  active_count = 0
    
  patch => realization%patch
  call PatchCountCells(patch,patch_total_count,patch_active_count)
  total_count = total_count + patch_total_count
  active_count = active_count + patch_active_count
  
  temp_int_in(1) = total_count
  temp_int_in(2) = active_count
  call MPI_Allreduce(temp_int_in,temp_int_out,TWO_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,realization%option%mycomm,ierr)
  global_total_count = temp_int_out(1)
  global_active_count = temp_int_out(2)

end subroutine RealizationCountCells

subroutine RealizationSetUpBC4Faces(realization)




  use Connection_module
  use Coupler_module
  use Patch_module
  use Grid_module
  use Field_module
  use MFD_Aux_module
  

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization

#ifdef DASVYAT

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  type(coupler_type), pointer ::  boundary_condition

  PetscReal, pointer :: bc_faces_p(:), xx_faces_p(:)
  PetscInt :: iconn, sum_connection, bc_type, bound_id
  PetscInt :: local_id, ghosted_id, ghost_face_id, j, jface, local_face_id
  PetscErrorCode :: ierr


  patch => realization%patch
  grid => patch%grid
  field => realization%field



  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecGetArrayF90(field%flow_xx_faces, xx_faces_p, ierr)

  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    bc_type = boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF)

    do iconn = 1, boundary_condition%numfaces_set
      sum_connection = sum_connection + 1

      local_id = boundary_condition%region%cell_ids(iconn)
      ghosted_id = grid%nL2G(local_id)

      aux_var => grid%MFD%aux_vars(local_id)
      do j = 1, aux_var%numfaces
        ghost_face_id = aux_var%face_id_gh(j)
        local_face_id = grid%fG2L(ghost_face_id)
        conn => grid%faces(ghost_face_id)%conn_set_ptr
        jface = grid%faces(ghost_face_id)%id
        if (boundary_condition%faces_set(iconn) == ghost_face_id) then
           if ((bc_type == DIRICHLET_BC).or.(bc_type == HYDROSTATIC_BC)  &
              .or.(bc_type == SEEPAGE_BC).or.(bc_type == CONDUCTANCE_BC) ) then
                    bc_faces_p(ghost_face_id) = boundary_condition%flow_aux_real_var(1,iconn)*conn%area(jface)
                    xx_faces_p(local_face_id) = boundary_condition%flow_aux_real_var(1,iconn)
           else if ((bc_type == NEUMANN_BC)) then
                    bc_faces_p(ghost_face_id) = boundary_condition%flow_aux_real_var(1,iconn)*conn%area(jface)
                    bound_id = grid%fL2B(local_face_id)
                    if (bound_id>0) then
                        patch%boundary_velocities(realization%option%nphase, bound_id) = &
                                             boundary_condition%flow_aux_real_var(1,iconn)
                    end if
           end if 
  !            bc_faces_p(ghost_face_id) = conn%cntr(3,jface)*conn%area(jface) 
        end if
      end do
    end do
    boundary_condition => boundary_condition%next
  end do


  call VecRestoreArrayF90(field%flow_xx_faces, xx_faces_p, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)



#ifdef DASVYAT
!   write(*,*) "Boundary faces"
!   do iconn = 1, grid%boundary_connection_set_list%first%num_connections 
!      write(*,*) "bound_flux", iconn, patch%boundary_velocities(realization%option%nphase, iconn),&
!               grid%boundary_connection_set_list%first%cntr(1,iconn), &
!               grid%boundary_connection_set_list%first%cntr(2,iconn), &
!               grid%boundary_connection_set_list%first%cntr(3,iconn)
!   end do  
!   read(*,*)
#endif

#endif

end subroutine RealizationSetUpBC4Faces

! ************************************************************************** !
!
! RealizationPrintGridStatistics: Prints statistics regarding the numerical
!                                 discretization 
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine RealizationPrintGridStatistics(realization)

  use Grid_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid

  PetscInt :: i1, i2, i3
  PetscReal :: r1, r2, r3
  PetscInt :: global_total_count, global_active_count
  PetscInt :: total_count, active_count
  PetscReal :: total_min, total_max, total_mean, total_variance
  PetscReal :: active_min, active_max, active_mean, active_variance
  PetscInt :: inactive_histogram(12), temp_int_out(12)
  PetscReal :: inactive_percentages(12)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%patch%grid

  ! print # of active and inactive grid cells
  call RealizationCountCells(realization,global_total_count, &
                             global_active_count,total_count,active_count)
  r1 = dble(total_count)
  call OptionMaxMinMeanVariance(r1,total_max, &
                                total_min,total_mean, &
                                total_variance,PETSC_TRUE,option)
  r1 = dble(active_count)
  call OptionMaxMinMeanVariance(r1,active_max, &
                                active_min,active_mean, &
                                active_variance,PETSC_TRUE,option)
                  
  r1 = dble(active_count) / dble(total_count)    
  inactive_histogram = 0                          
  if (r1 >= (1.d0-1.d-8)) then
    inactive_histogram(12) = 1
  else if (r1 >= .9d0 .and. r1 < (1.d0-1.d-8)) then
    inactive_histogram(11) = 1
  else if (r1 >= .8d0 .and. r1 < .9d0) then
    inactive_histogram(10) = 1
  else if (r1 >= .7d0 .and. r1 < .8d0) then
    inactive_histogram(9) = 1
  else if (r1 >= .6d0 .and. r1 < .7d0) then
    inactive_histogram(8) = 1
  else if (r1 >= .5d0 .and. r1 < .6d0) then
    inactive_histogram(7) = 1
  else if (r1 >= .4d0 .and. r1 < .5d0) then
    inactive_histogram(6) = 1
  else if (r1 >= .3d0 .and. r1 < .4d0) then
    inactive_histogram(5) = 1
  else if (r1 >= .2d0 .and. r1 < .3d0) then
    inactive_histogram(4) = 1
  else if (r1 >= .1d0 .and. r1 < .2d0) then
    inactive_histogram(3) = 1
  else if (r1 > 1.d-20 .and. r1 < .1d0) then
    inactive_histogram(2) = 1
  else if (r1 < 1.d-20) then
    inactive_histogram(1) = 1
  endif
  
  call MPI_Allreduce(inactive_histogram,temp_int_out,TWELVE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! why I cannot use *100, I do not know....geh
  inactive_percentages = dble(temp_int_out)/dble(option%mycommsize)*10.d0
  inactive_percentages = inactive_percentages+1.d-8

  r1 = 0.d0
  do i1 = 1, 12
    r1 = r1 + inactive_percentages(i1)
  enddo
                                
  i1 = -999
  i2 = -999
  i3 = -999
  if (associated(grid%structured_grid)) then
    i1 = grid%structured_grid%npx_final
    i2 = grid%structured_grid%npy_final
    i3 = grid%structured_grid%npz_final
  endif
  if (OptionPrintToScreen(option)) then
    write(*,'(/," Grid Stats:",/, &
              & "                       Global # cells: ",i12,/, &
              & "                Global # active cells: ",i12,/, &
              & "                              # cores: ",i12,/, &
              & "         Processor core decomposition: ",3i6,/, &
              & "               Maximum # cells / core: ",i12,/, &
              & "               Minimum # cells / core: ",i12,/, &
              & "               Average # cells / core: ",1pe12.4,/, &
              & "               Std Dev # cells / core: ",1pe12.4,/, &
              & "        Maximum # active cells / core: ",i12,/, &
              & "        Minimum # active cells / core: ",i12,/, &
              & "        Average # active cells / core: ",1pe12.4,/, &
              & "        Std Dev # active cells / core: ",1pe12.4,/,/, &
              & "        % cores with % active cells =       0%: ",1f7.2,/, &
              & "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
              & "        % cores with % active cells =   10-20%: ",1f7.2,/, &
              & "        % cores with % active cells =   20-30%: ",1f7.2,/, &
              & "        % cores with % active cells =   30-40%: ",1f7.2,/, &
              & "        % cores with % active cells =   40-50%: ",1f7.2,/, &
              & "        % cores with % active cells =   50-60%: ",1f7.2,/, &
              & "        % cores with % active cells =   60-70%: ",1f7.2,/, &
              & "        % cores with % active cells =   70-80%: ",1f7.2,/, &
              & "        % cores with % active cells =   80-90%: ",1f7.2,/, &
              & "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
              & "        % cores with % active cells =     100%: ",1f7.2,/, &
              & "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'(/," Grid Stats:",/, &
               & "                       Global # cells: ",i12,/, &
               & "                Global # active cells: ",i12,/, &
               & "                              # cores: ",i12,/, &
               & "         Processor core decomposition: ",3i6,/, &
               & "               Maximum # cells / core: ",i12,/, &
               & "               Minimum # cells / core: ",i12,/, &
               & "               Average # cells / core: ",1pe12.4,/, &
               & "               Std Dev # cells / core: ",1pe12.4,/, &
               & "        Maximum # active cells / core: ",i12,/, &
               & "        Minimum # active cells / core: ",i12,/, &
               & "        Average # active cells / core: ",1pe12.4,/, &
               & "        Std Dev # active cells / core: ",1pe12.4,/,/, &
               & "        % cores with % active cells =       0%: ",1f7.2,/, &
               & "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
               & "        % cores with % active cells =   10-20%: ",1f7.2,/, &
               & "        % cores with % active cells =   20-30%: ",1f7.2,/, &
               & "        % cores with % active cells =   30-40%: ",1f7.2,/, &
               & "        % cores with % active cells =   40-50%: ",1f7.2,/, &
               & "        % cores with % active cells =   50-60%: ",1f7.2,/, &
               & "        % cores with % active cells =   60-70%: ",1f7.2,/, &
               & "        % cores with % active cells =   70-80%: ",1f7.2,/, &
               & "        % cores with % active cells =   80-90%: ",1f7.2,/, &
               & "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
               & "        % cores with % active cells =     100%: ",1f7.2,/, &
               & "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif

end subroutine RealizationPrintGridStatistics

! ************************************************************************** !
!
! RealizationCalculateCFL1Timestep: Calculates largest time step that  
!                                   preserves a CFL # of 1 in a realization
! author: Glenn Hammond
! date: 10/07/11
!
! ************************************************************************** !
subroutine RealizationCalculateCFL1Timestep(realization,max_dt_cfl_1)

  implicit none

  type(realization_type) realization
  PetscReal :: max_dt_cfl_1
  
  type(patch_type), pointer :: patch
  PetscReal :: max_dt_cfl_1_patch
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  
  max_dt_cfl_1 = 1.d20
  patch => realization%patch
  call PatchCalculateCFL1Timestep(patch,realization%option, &
                                  max_dt_cfl_1_patch)
  max_dt_cfl_1 = min(max_dt_cfl_1,max_dt_cfl_1_patch)

  ! get the minimum across all cores
  call MPI_Allreduce(max_dt_cfl_1,tempreal,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN, &
                     realization%option%mycomm,ierr)
  max_dt_cfl_1 = tempreal

end subroutine RealizationCalculateCFL1Timestep

! ************************************************************************** !
!
! RealizationDestroy: Deallocates a realization
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RealizationDestroy(realization)

  implicit none
  
  type(realization_type), pointer :: realization
  
  if (.not.associated(realization)) return
    
  call FieldDestroy(realization%field)

!  call OptionDestroy(realization%option) !geh it will be destroy externally
  call OutputOptionDestroy(realization%output_option)
  call RegionDestroyList(realization%regions)
  
  call FlowConditionDestroyList(realization%flow_conditions)
  call TranConditionDestroyList(realization%transport_conditions)
  call TranConstraintDestroyList(realization%transport_constraints)

  call LevelDestroyList(realization%level_list)

  if (associated(realization%debug)) deallocate(realization%debug)
  nullify(realization%debug)
  
  if (associated(realization%fluid_property_array)) &
    deallocate(realization%fluid_property_array)
  nullify(realization%fluid_property_array)
  call FluidPropertyDestroy(realization%fluid_properties)
  
  if (associated(realization%material_property_array)) &
    deallocate(realization%material_property_array)
  nullify(realization%material_property_array)
  call MaterialPropertyDestroy(realization%material_properties)

  if (associated(realization%saturation_function_array)) &
    deallocate(realization%saturation_function_array)
  nullify(realization%saturation_function_array)
  call SaturationFunctionDestroy(realization%saturation_functions)

  call DatasetDestroy(realization%datasets)
  
  call VelocityDatasetDestroy(realization%velocity_dataset)
  
  call DiscretizationDestroy(realization%discretization)
  
  call ReactionDestroy(realization%reaction)
  
end subroutine RealizationDestroy

end module Realization_module
