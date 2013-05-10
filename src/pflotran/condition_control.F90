module Condition_Control_module
  
  ! This module store routines that operate on conditions from a level above 
  ! that of the realization_module.  This is necessary to access capability
  ! such as HDF5 which is unavailable from within the realization object
  ! and below.  Routines in this module will loop over realization, levels,
  ! and patches without calling underlying level/patch versions of the
  ! subroutines, which is common in realization.F90 - GEH
 
  implicit none

  private

#include "definitions.h"

  public :: CondControlAssignFlowInitCond, &
            CondControlAssignTranInitCond, &
#ifdef SURFACE_FLOW
            CondControlAssignFlowInitCondSurface, &
#endif
            CondControlScaleSourceSink
 
contains

! ************************************************************************** !
!
! CondControlAssignFlowInitCond: Assigns flow initial conditions to model
! author: Glenn Hammond
! date: 11/02/07, 10/18/11
!
! ************************************************************************** !
subroutine CondControlAssignFlowInitCond(realization)

  use Realization_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Dataset_Aux_module
  use Grid_module
  use Level_module
  use Patch_module
  use Water_EOS_module
#ifdef DASVYAT
  use MFD_module, only : MFDInitializeMassMatrices
#endif

  use Global_Aux_module
  use General_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof, iface
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:), xx_faces_p(:)
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(flow_general_condition_type), pointer :: general
  type(dataset_type), pointer :: dataset
  type(global_auxvar_type) :: global_aux
  type(general_auxvar_type) :: general_aux
  PetscBool :: use_dataset
  PetscBool :: dataset_flag(realization%option%nflowdof)
  PetscInt :: num_connections
  PetscInt, pointer :: conn_id_ptr(:)
  PetscInt :: ghosted_offset
  PetscReal :: x(realization%option%nflowdof)
  PetscReal :: temperature, p_sat

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  if (option%iflowmode == G_MODE) then
    call GlobalAuxVarInit(global_aux,option)
    call GeneralAuxVarInit(general_aux,option)
  endif
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid

      select case(option%iflowmode)
      
        case(G_MODE) ! general phase mode

          call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
          call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
      
          xx_p = -999.d0
      
          initial_condition => cur_patch%initial_conditions%first
          do
      
            if (.not.associated(initial_condition)) exit

            if (.not.associated(initial_condition%flow_aux_real_var)) then
              if (.not.associated(initial_condition%flow_condition)) then
                option%io_buffer = 'Flow condition is NULL in initial condition'
                call printErrMsg(option)
              endif
              
              general => initial_condition%flow_condition%general
              
              string = 'in flow condition "' // &
                trim(initial_condition%flow_condition%name) // &
                '" within initial condition "' // &
                trim(initial_condition%flow_condition%name) // &
                '" must be of type Dirichlet or Hydrostatic'
              ! error checking.  the data must match the state
              select case(initial_condition%flow_condition%iphase)
                case(TWO_PHASE_STATE)  
                  if (.not. &
                      (general%gas_pressure%itype == DIRICHLET_BC .or. &
                       general%gas_pressure%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Gas pressure ' // trim(string)
                    call printErrMsg(option)
                  endif
                  if (.not. &
                      (general%gas_saturation%itype == DIRICHLET_BC .or. &
                       general%gas_saturation%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Gas saturation ' // trim(string)
                    call printErrMsg(option)
                  endif
                case(LIQUID_STATE)
                  if (.not. &
                      (general%liquid_pressure%itype == DIRICHLET_BC .or. &
                       general%liquid_pressure%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Liquid pressure ' // trim(string)
                    call printErrMsg(option)
                  endif
                  if (.not. &
                      (general%mole_fraction%itype == DIRICHLET_BC .or. &
                       general%mole_fraction%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Mole fraction ' // trim(string)
                    call printErrMsg(option)
                  endif
                case(GAS_STATE)
                  if (.not. &
                      (general%gas_pressure%itype == DIRICHLET_BC .or. &
                       general%gas_pressure%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Gas pressure ' // trim(string)
                    call printErrMsg(option)
                  endif
                  if (.not. &
                      (general%mole_fraction%itype == DIRICHLET_BC .or. &
                       general%mole_fraction%itype == HYDROSTATIC_BC)) then
                    option%io_buffer = 'Gas saturation ' // trim(string)
                    call printErrMsg(option)
                  endif
              end select
              if (.not. &
                  (general%temperature%itype == DIRICHLET_BC .or. &
                   general%temperature%itype == HYDROSTATIC_BC)) then
                option%io_buffer = 'Temperature ' // trim(string)
                call printErrMsg(option)
              endif                              
              
              
              do icell=1,initial_condition%region%num_cells
                local_id = initial_condition%region%cell_ids(icell)
                ghosted_id = grid%nL2G(local_id)
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                if (cur_patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  iphase_loc_p(ghosted_id) = 0
                  cycle
                endif
                ! decrement ibegin to give a local offset of 0
                ibegin = ibegin - 1
                select case(initial_condition%flow_condition%iphase)
                  case(TWO_PHASE_STATE)
                    xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                      general%gas_pressure%flow_dataset%time_series%cur_value(1)
                    xx_p(ibegin+GENERAL_GAS_SATURATION_DOF) = &
                      general%gas_saturation%flow_dataset%time_series%cur_value(1)
                    temperature = general%temperature%flow_dataset%time_series%cur_value(1)
                    call psat(temperature,p_sat,ierr)
                    ! p_a = p_g - p_s(T)
                    xx_p(ibegin+GENERAL_AIR_PRESSURE_DOF) = &
                      general%gas_pressure%flow_dataset%time_series%cur_value(1) - &
                      p_sat
                  case(LIQUID_STATE)
                    xx_p(ibegin+GENERAL_LIQUID_PRESSURE_DOF) = &
                      general%liquid_pressure%flow_dataset%time_series%cur_value(1)
                    xx_p(ibegin+GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF) = &
                      general%mole_fraction%flow_dataset%time_series%cur_value(1)
                    xx_p(ibegin+GENERAL_LIQUID_STATE_TEMPERATURE_DOF) = &
                      general%temperature%flow_dataset%time_series%cur_value(1)
                  case(GAS_STATE)
                    xx_p(ibegin+GENERAL_GAS_PRESSURE_DOF) = &
                      general%gas_pressure%flow_dataset%time_series%cur_value(1)
                    xx_p(ibegin+GENERAL_AIR_PRESSURE_DOF) = &
                      general%gas_pressure%flow_dataset%time_series%cur_value(1) * &
                      general%mole_fraction%flow_dataset%time_series%cur_value(1)
                    xx_p(ibegin+GENERAL_GAS_STATE_TEMPERATURE_DOF) = &
                      general%temperature%flow_dataset%time_series%cur_value(1)
                end select
                iphase_loc_p(ghosted_id) = initial_condition%flow_condition%iphase
                cur_patch%aux%Global%aux_vars(ghosted_id)%istate = &
                  initial_condition%flow_condition%iphase
              enddo
            else
              do iconn=1,initial_condition%connection_set%num_connections
                local_id = initial_condition%connection_set%id_dn(iconn)
                ghosted_id = grid%nL2G(local_id)
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                if (cur_patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  iphase_loc_p(ghosted_id) = 0
                  cycle
                endif
                xx_p(ibegin:iend) = &
                  initial_condition%flow_aux_real_var(1:option%nflowdof,iconn)
                iphase_loc_p(ghosted_id) = initial_condition%flow_condition%iphase
                cur_patch%aux%Global%aux_vars(ghosted_id)%istate = &
                  initial_condition%flow_condition%iphase
              enddo
            endif
            initial_condition => initial_condition%next
          enddo
     
          call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
          call GridVecRestoreArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
              
        case default
          ! assign initial conditions values to domain
          if (discretization%itype == STRUCTURED_GRID_MIMETIC.or. &
              discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then
            call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
            call VecGetArrayF90(field%flow_xx_faces, xx_faces_p, ierr); CHKERRQ(ierr)        
          else
            call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
          end if
          call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
      
          xx_p = -999.d0
      
          initial_condition => cur_patch%initial_conditions%first
          do
      
            if (.not.associated(initial_condition)) exit

            if (discretization%itype == STRUCTURED_GRID_MIMETIC.or. &
                discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then
#ifdef DASVYAT
              use_dataset = PETSC_FALSE
              dataset_flag = PETSC_FALSE
              do idof = 1, option%nflowdof
                dataset =>  initial_condition%flow_condition% &
                                 sub_condition_ptr(idof)%ptr% &
                                 flow_dataset%dataset
                if (associated(dataset)) then
                  if (dataset%is_cell_indexed) then
                    use_dataset = PETSC_TRUE
                    dataset_flag(idof) = PETSC_TRUE
                    call ConditionControlMapDatasetToVec(realization, &
                           initial_condition%flow_condition% &
                             sub_condition_ptr(idof)%ptr% &
                             flow_dataset%dataset,idof,field%flow_xx,GLOBAL)
                  endif
                endif
              enddo
              if (.not.associated(initial_condition%flow_aux_real_var)) then
                conn_id_ptr => initial_condition%region%cell_ids

                do iconn=1,initial_condition%region%num_cells
                  local_id = conn_id_ptr(iconn)
                  ghosted_id = grid%nL2G(local_id)
                  iend = local_id*option%nflowdof
                  ibegin = iend-option%nflowdof+1
                  do idof = 1, option%nflowdof
                    if (.not.dataset_flag(idof)) then
                      xx_p(ibegin+idof-1) = &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%flow_dataset% &
                          time_series%cur_value(1)
                    endif
                  enddo
                  xx_faces_p(ibegin:iend) = xx_p(ibegin:iend) ! for LP -formulation
                  xx_faces_p(grid%nlmax_faces + &
                            ibegin:grid%nlmax_faces + iend) = &
                                xx_p(ibegin:iend) ! for LP -formulation
                enddo

                do icell=1,initial_condition%region%num_cells
                  local_id = initial_condition%region%cell_ids(icell)
                  ghosted_id = grid%nL2G(local_id)
                  if (cur_patch%imat(ghosted_id) <= 0) then
                    iphase_loc_p(ghosted_id) = 0
                    cycle
                  endif
                  iphase_loc_p(ghosted_id)=initial_condition%flow_condition%iphase
                 enddo
              else
                do iface=1,initial_condition%numfaces_set
                  ghosted_id = initial_condition%faces_set(iface)
                  local_id = grid%fG2L(ghosted_id)
                  if (local_id > 0) then
                    iend = local_id*option%nflowdof
                    ibegin = iend-option%nflowdof+1
                    xx_faces_p(ibegin:iend) = &
                    initial_condition%flow_aux_real_var(1:option%nflowdof,iface)
                  endif
                enddo
                do iconn=1,initial_condition%connection_set%num_connections
                  local_id = initial_condition%region%cell_ids(iconn)
                  ghosted_id = grid%nL2G(local_id)
                  iend = local_id*option%nflowdof
                  ibegin = iend-option%nflowdof+1
                  if (cur_patch%imat(ghosted_id) <= 0) then
                    xx_p(ibegin:iend) = 0.d0
                    iphase_loc_p(ghosted_id) = 0
                    cycle
                  endif
                  xx_p(ibegin:iend) = &
                    initial_condition%flow_aux_real_var(1:option%nflowdof, &
                                                        iconn + &
                                                        initial_condition%numfaces_set)
                  xx_faces_p(grid%nlmax_faces + &
                            ibegin:grid%nlmax_faces + iend) = &
                                xx_p(ibegin:iend) ! for LP -formulation
                  iphase_loc_p(ghosted_id) = &
                    initial_condition%flow_aux_int_var(1,iconn + &
                                                      initial_condition%numfaces_set)
                enddo
              endif
#endif
            else 
              use_dataset = PETSC_FALSE
              dataset_flag = PETSC_FALSE
              do idof = 1, option%nflowdof
                dataset =>  initial_condition%flow_condition% &
                                 sub_condition_ptr(idof)%ptr% &
                                 flow_dataset%dataset
                if (associated(dataset)) then
                  if (dataset%is_cell_indexed) then
                    use_dataset = PETSC_TRUE
                    dataset_flag(idof) = PETSC_TRUE
                    call ConditionControlMapDatasetToVec(realization, &
                           initial_condition%flow_condition% &
                             sub_condition_ptr(idof)%ptr% &
                             flow_dataset%dataset,idof,field%flow_xx,GLOBAL)
                  endif
                endif
              enddo            
              if (.not.associated(initial_condition%flow_aux_real_var) .and. &
                  .not.associated(initial_condition%flow_condition)) then
                option%io_buffer = 'Flow condition is NULL in initial condition'
                call printErrMsg(option)
              endif
              if (associated(initial_condition%flow_aux_real_var)) then
                num_connections = &
                  initial_condition%connection_set%num_connections
                conn_id_ptr => initial_condition%connection_set%id_dn
              else
                num_connections = initial_condition%region%num_cells
                conn_id_ptr => initial_condition%region%cell_ids
              endif
              do iconn=1, num_connections
                local_id = conn_id_ptr(iconn)
                ghosted_id = grid%nL2G(local_id)
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                if (cur_patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  iphase_loc_p(ghosted_id) = 0
                  cycle
                endif
                if (associated(initial_condition%flow_aux_real_var)) then
                  do idof = 1, option%nflowdof
                    if (.not.dataset_flag(idof)) then
                      xx_p(ibegin+idof-1) =  &
                        initial_condition%flow_aux_real_var(idof,iconn)
                    endif
                  enddo
                else
                  do idof = 1, option%nflowdof
                    if (.not.dataset_flag(idof)) then
                      xx_p(ibegin+idof-1) = &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%flow_dataset% &
                          time_series%cur_value(1)
                    endif
                  enddo
                endif
                iphase_loc_p(ghosted_id) = &
                  initial_condition%flow_condition%iphase
                if (option%iflowmode == G_MODE) then
                  cur_patch%aux%Global%aux_vars(ghosted_id)%istate = &
                    int(iphase_loc_p(ghosted_id))
                endif
              enddo
            end if
            initial_condition => initial_condition%next
          enddo
     
          call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)

      end select 
   
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  if (option%iflowmode == G_MODE) then
    call GlobalAuxVarStrip(global_aux)
    call GeneralAuxVarStrip(general_aux)
  endif  
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx,field%flow_xx_loc,NFLOWDOF)  
  

  call VecCopy(field%flow_xx, field%flow_yy, ierr)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)  
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_old_loc,ONEDOF)

#ifdef DASVYAT
  if (discretization%itype == STRUCTURED_GRID_MIMETIC.or. &
      discretization%itype == UNSTRUCTURED_GRID_MIMETIC) then

    call VecRestoreArrayF90(field%flow_xx_faces,xx_faces_p, ierr)
    call RealizationSetUpBC4Faces(realization)

    !call DiscretizationGlobalToLocalFaces(discretization, field%flow_xx_faces, field%flow_xx_loc_faces, NFLOWDOF)
    call DiscretizationGlobalToLocalLP(discretization, field%flow_xx_faces, field%flow_xx_loc_faces, NFLOWDOF)
    call VecCopy(field%flow_xx_faces, field%flow_yy_faces, ierr)
    call MFDInitializeMassMatrices(realization%discretization%grid,&
                                  realization%field, &
                                  realization%discretization%MFD, realization%option)
    patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_TRUE

  endif
#endif
!  stop
end subroutine CondControlAssignFlowInitCond

! ************************************************************************** !
!
! CondControlAssignTranInitCond: Assigns transport initial conditions to model
! author: Glenn Hammond
! date: 11/02/07, 10/18/11
!
! ************************************************************************** !
subroutine CondControlAssignTranInitCond(realization)

  use Realization_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Constraint_module
  use Grid_module
  use Dataset_Aux_module
  use Level_module
  use Patch_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Reaction_module
  use HDF5_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof, isub_condition, temp_int, iimmobile
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscInt :: irxn, isite, imnrl, ikinrxn
  PetscReal, pointer :: xx_p(:), xx_loc_p(:), porosity_loc(:), vec_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(tran_constraint_coupler_type), pointer :: constraint_coupler

  PetscInt :: iphase
  PetscInt :: offset
  PetscBool :: re_equilibrate_at_each_cell
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(dataset_type), pointer :: dataset
  PetscInt :: aq_dataset_to_idof(realization%reaction%naqcomp)
  PetscInt :: iaqdataset, num_aq_datasets
  PetscBool :: use_aq_dataset
  PetscReal :: ave_num_iterations
  PetscReal :: tempreal
  PetscLogDouble :: tstart, tend
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  
  iphase = 1
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      rt_aux_vars => cur_patch%aux%RT%aux_vars
      global_aux_vars => cur_patch%aux%Global%aux_vars

      ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid,field%tran_xx,xx_p,ierr)
      call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc,ierr)
      
      xx_p = -999.d0
      
      initial_condition => cur_patch%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit
        
        constraint_coupler => initial_condition%tran_condition%cur_constraint_coupler

        re_equilibrate_at_each_cell = PETSC_FALSE
        use_aq_dataset = PETSC_FALSE
        num_aq_datasets = 0
        aq_dataset_to_idof = 0
        do idof = 1, reaction%naqcomp ! primary aqueous concentrations
          if (constraint_coupler%aqueous_species%external_dataset(idof)) then
            num_aq_datasets = num_aq_datasets + 1
            aq_dataset_to_idof(num_aq_datasets) = idof
            re_equilibrate_at_each_cell = PETSC_TRUE
            use_aq_dataset = PETSC_TRUE
            string = 'constraint ' // trim(constraint_coupler%constraint_name)
            dataset => DatasetGetPointer(realization%datasets, &
                         constraint_coupler%aqueous_species%constraint_aux_string(idof), &
                         string,option)
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 field%tran_xx_loc,LOCAL)
          endif
        enddo

        ! read in heterogeneous mineral volume fractions
        if (associated(constraint_coupler%minerals)) then
          do imnrl = 1, reaction%mineral%nkinmnrl
            if (constraint_coupler%minerals%external_dataset(imnrl)) then
              re_equilibrate_at_each_cell = PETSC_TRUE
              string = 'constraint ' // trim(constraint_coupler%constraint_name)
              dataset => DatasetGetPointer(realization%datasets, &
                           constraint_coupler%minerals%constraint_aux_string(imnrl), &
                           string,option)
              string = '' ! group name
              string2 = dataset%h5_dataset_name ! dataset name
              call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                                dataset%filename, &
                                                string,string2, &
                                                dataset%realization_dependent)
              call DiscretizationGlobalToLocal(discretization,field%work,field%work_loc,ONEDOF)
              call GridVecGetArrayF90(grid,field%work_loc,vec_p,ierr)
              do icell=1,initial_condition%region%num_cells
                local_id = initial_condition%region%cell_ids(icell)
                ghosted_id = grid%nL2G(local_id)
                rt_aux_vars(ghosted_id)%mnrl_volfrac0(imnrl) = vec_p(ghosted_id)
                rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = vec_p(ghosted_id)
              enddo
              call GridVecRestoreArrayF90(grid,field%work_loc,vec_p,ierr)
            endif
          enddo
        endif
          
        ! read in heterogeneous immobile
        if (associated(constraint_coupler%immobile_species)) then
          do iimmobile = 1, reaction%immobile%nimmobile
            if (constraint_coupler%immobile_species%external_dataset(iimmobile)) then
              ! no need to requilibrate at each cell
              string = 'constraint ' // trim(constraint_coupler%constraint_name)
              dataset => DatasetGetPointer(realization%datasets, &
                  constraint_coupler%immobile_species%constraint_aux_string(iimmobile), &
                  string,option)
              string = '' ! group name
              string2 = dataset%h5_dataset_name ! dataset name
              call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                                dataset%filename, &
                                                string,string2, &
                                                dataset%realization_dependent)
              call DiscretizationGlobalToLocal(discretization,field%work,field%work_loc,ONEDOF)
              call GridVecGetArrayF90(grid,field%work_loc,vec_p,ierr)
              do icell=1,initial_condition%region%num_cells
                local_id = initial_condition%region%cell_ids(icell)
                ghosted_id = grid%nL2G(local_id)
                rt_aux_vars(ghosted_id)%immobile(iimmobile) = vec_p(ghosted_id)
              enddo
              call GridVecRestoreArrayF90(grid,field%work_loc,vec_p,ierr)
            endif
          enddo
        endif
          
        if (.not.option%use_isothermal) then
          re_equilibrate_at_each_cell = PETSC_TRUE
        endif
        
        if (use_aq_dataset) then
          call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p,ierr); CHKERRQ(ierr)
          call PetscTime(tstart,ierr) 
        endif
        
        ave_num_iterations = 0.d0
        do icell=1,initial_condition%region%num_cells
          local_id = initial_condition%region%cell_ids(icell)
          ghosted_id = grid%nL2G(local_id)
          iend = local_id*option%ntrandof
          ibegin = iend-option%ntrandof+1
          if (cur_patch%imat(ghosted_id) <= 0) then
            xx_p(ibegin:iend) = 1.d-200
            cycle
          endif
          if (re_equilibrate_at_each_cell) then
            if (use_aq_dataset) then
              offset = (ghosted_id-1)*option%ntrandof
              do iaqdataset = 1, num_aq_datasets
                ! remember that xx_loc_p holds the data set values that were read in
                temp_int = aq_dataset_to_idof(iaqdataset)
                constraint_coupler%aqueous_species%constraint_conc(temp_int) = &
                  xx_loc_p(offset+temp_int)
              enddo
            endif
            option%iflag = grid%nG2A(grid%nL2G(local_id))
            if (icell == 1) then
              call ReactionEquilibrateConstraint(rt_aux_vars(ghosted_id), &
                global_aux_vars(ghosted_id),reaction, &
                constraint_coupler%constraint_name, &
                constraint_coupler%aqueous_species, &
                constraint_coupler%minerals, &
                constraint_coupler%surface_complexes, &
                constraint_coupler%colloids, &
                constraint_coupler%immobile_species, &
                porosity_loc(ghosted_id), &
                constraint_coupler%num_iterations, &
                PETSC_FALSE,option)
            else
!geh              call RTAuxVarCopy(rt_aux_vars(ghosted_id), &
!geh                rt_aux_vars(grid%nL2G(initial_condition%region%cell_ids(icell-1))), &
!geh                option)
              rt_aux_vars(ghosted_id)%pri_molal = &
                rt_aux_vars(grid%nL2G(initial_condition%region%cell_ids(icell-1)))%pri_molal
              call ReactionEquilibrateConstraint(rt_aux_vars(ghosted_id), &
                global_aux_vars(ghosted_id),reaction, &
                constraint_coupler%constraint_name, &
                constraint_coupler%aqueous_species, &
                constraint_coupler%minerals, &
                constraint_coupler%surface_complexes, &
                constraint_coupler%colloids, &
                constraint_coupler%immobile_species, &
                porosity_loc(ghosted_id), &
                constraint_coupler%num_iterations, &
                PETSC_TRUE,option)
            endif
            option%iflag = 0
            ave_num_iterations = ave_num_iterations + &
              constraint_coupler%num_iterations
          endif
          ! ibegin is the local non-ghosted offset: (local_id-1)*option%ntrandof+1
          offset = ibegin + reaction%offset_aqueous - 1
          ! primary aqueous concentrations
          do idof = 1, reaction%naqcomp 
            xx_p(offset+idof) = &
              constraint_coupler%aqueous_species%basis_molarity(idof) / &
              global_aux_vars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
          enddo
          ! mineral volume fractions
          if (associated(constraint_coupler%minerals)) then
            do imnrl = 1, reaction%mineral%nkinmnrl
              ! if read from a dataset, the vol frac was set above.  Don't want to
              ! overwrite
              if (.not.constraint_coupler%minerals%external_dataset(imnrl)) then
                rt_aux_vars(ghosted_id)%mnrl_volfrac0(imnrl) = &
                  constraint_coupler%minerals%constraint_vol_frac(imnrl)
                rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = &
                  constraint_coupler%minerals%constraint_vol_frac(imnrl)
              endif
              rt_aux_vars(ghosted_id)%mnrl_area0(imnrl) = &
                constraint_coupler%minerals%constraint_area(imnrl)
              rt_aux_vars(ghosted_id)%mnrl_area(imnrl) = &
                constraint_coupler%minerals%constraint_area(imnrl)
            enddo
          endif
          ! kinetic surface complexes
          if (associated(constraint_coupler%surface_complexes)) then
            do idof = 1, reaction%surface_complexation%nkinsrfcplx
              rt_aux_vars(ghosted_id)%kinsrfcplx_conc(idof,-1) = & !geh: to catch bug
                constraint_coupler%surface_complexes%constraint_conc(idof)
            enddo
            do ikinrxn = 1, reaction%surface_complexation%nkinsrfcplxrxn
              irxn = reaction%surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
              isite = reaction%surface_complexation%srfcplxrxn_to_surf(irxn)
              rt_aux_vars(ghosted_id)%kinsrfcplx_free_site_conc(isite) = &
                constraint_coupler%surface_complexes%basis_free_site_conc(isite)
            enddo
          endif
          ! this is for the multi-rate surface complexation model
          if (reaction%surface_complexation%nkinmrsrfcplxrxn > 0 .and. &
            ! geh: if we re-equilibrate at each grid cell, we do not want to
            ! overwrite the reequilibrated values with those from the constraint
              .not. re_equilibrate_at_each_cell) then
            ! copy over total sorbed concentration
            rt_aux_vars(ghosted_id)%kinmr_total_sorb = &
              constraint_coupler%rt_auxvar%kinmr_total_sorb
            ! copy over free site concentration
            rt_aux_vars(ghosted_id)%srfcplxrxn_free_site_conc = &
              constraint_coupler%rt_auxvar%srfcplxrxn_free_site_conc
          endif
          ! colloids fractions
          if (associated(constraint_coupler%colloids)) then
            offset = ibegin + reaction%offset_colloid - 1
            do idof = 1, reaction%ncoll ! primary aqueous concentrations
              xx_p(offset+idof) = &
                constraint_coupler%colloids%basis_conc_mob(idof) / &
                global_aux_vars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
              rt_aux_vars(ghosted_id)%colloid%conc_imb(idof) = &
                constraint_coupler%colloids%basis_conc_imb(idof)
            enddo
          endif
          ! immobile
          if (associated(constraint_coupler%immobile_species)) then
            offset = ibegin + reaction%offset_immobile - 1
            do iimmobile = 1, reaction%immobile%nimmobile
              if (constraint_coupler%immobile_species%external_dataset(iimmobile)) then
                ! already read into rt_aux_vars above.
                xx_p(offset+iimmobile) = &
                  rt_aux_vars(ghosted_id)%immobile(iimmobile)
              else
                xx_p(offset+iimmobile) = &
                  constraint_coupler%immobile_species%constraint_conc(iimmobile)
                rt_aux_vars(ghosted_id)%immobile(iimmobile) = &
                  constraint_coupler%immobile_species%constraint_conc(iimmobile)
              endif
            enddo
          endif
        enddo ! icell=1,initial_condition%region%num_cells
        if (use_aq_dataset) then
          call PetscTime(tend,ierr) 
          call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p,ierr)
          ave_num_iterations = ave_num_iterations / &
            initial_condition%region%num_cells
          write(option%io_buffer,&
                '("Average number of iterations in ReactionEquilibrateConstraint():", &
                & f5.1)') ave_num_iterations
          call printMsg(option)
          write(option%io_buffer,'(f10.2," Seconds to equilibrate constraints")') &
            tend-tstart
          call printMsg(option)
        endif
        initial_condition => initial_condition%next
      enddo
      
      call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
      call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc,ierr)

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  ! check to ensure that minimum concentration is not less than or equal
  ! to zero
  call VecMin(field%tran_xx,PETSC_NULL_INTEGER,tempreal,ierr)
  if (tempreal <= 0.d0) then
    option%io_buffer = 'ERROR: Zero concentrations found in initial ' // &
      'transport solution.'
    call printMsg(option)
    ! now figure out which species have zero concentrations
    do idof = 1, option%ntrandof
      call VecStrideMin(field%tran_xx,idof-1,offset,tempreal,ierr)
      if (tempreal <= 0.d0) then
        write(string,*) tempreal
        if (idof <= reaction%naqcomp) then
          string2 = '  Aqueous species "' // &
            trim(reaction%primary_species_names(idof))
        else
          string2 = '  Immobile species "' // &
            trim(reaction%immobile%names(idof-reaction%offset_immobile))
        endif
          string2 = trim(string2) // &
            '" has zero concentration (' // &
            trim(adjustl(string)) // ').'
        call printMsg(option)
      endif
    enddo
    option%io_buffer = 'Free ion concentations must be positive.  Try ' // &
      'using a small value such as 1.e-20 or 1.e-40 instead of zero.'
    call printErrMsg(option)
  endif
  
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr)

end subroutine CondControlAssignTranInitCond

! ************************************************************************** !
!
! ConditionControlMapDatasetToVec: maps an external dataset to a PETSc vec
!                                  representing values at each grid cell
! author: Glenn Hammond
! date: 03/23/12
!
! ************************************************************************** !
subroutine ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                           mdof_vec,vec_type)
  use Realization_class
  use Option_module
  use Field_module
  use Dataset_Aux_module
  use HDF5_module
  use Discretization_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  
  
  type(realization_type) :: realization
  type(dataset_type), pointer :: dataset
  PetscInt :: idof
  Vec :: mdof_vec
  PetscInt :: vec_type

  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
  
  if (associated(dataset)) then
    string = '' ! group name
    ! have to copy to string2 due to mismatch in string size
    string2 = dataset%h5_dataset_name
    call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                      dataset%filename, &
                                      string,string2, &
                                      dataset%realization_dependent)
    if (vec_type == GLOBAL) then
      call VecStrideScatter(field%work,idof-1,mdof_vec, &
                            INSERT_VALUES,ierr)    
    else
      call DiscretizationGlobalToLocal(realization%discretization, &
                                       field%work, &
                                       field%work_loc,ONEDOF)
      call VecStrideScatter(field%work_loc,idof-1,mdof_vec, &
                            INSERT_VALUES,ierr)    
    endif
  endif

end subroutine ConditionControlMapDatasetToVec

! ************************************************************************** !
!
! CondControlScaleSourceSink: Scales select source/sinks based on perms
! author: Glenn Hammond
! date: 09/03/08, 10/18/11
!
! ************************************************************************** !
subroutine CondControlScaleSourceSink(realization)

  use Realization_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Level_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdmda.h"

  
  type(realization_type) :: realization
  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: perm_loc_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount
  PetscInt :: x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  
  PetscInt :: ghosted_neighbors(0:27)
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  ! GB: grid was uninitialized
  grid => patch%grid
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_loc_ptr,ierr)

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      ! BIG-TIME warning here.  I assume that all source/sink cells are within 
      ! a single patch - geh

      grid => cur_patch%grid

      cur_source_sink => cur_patch%source_sinks%first
      do
        if (.not.associated(cur_source_sink)) exit

        call VecZeroEntries(field%work,ierr)
        call GridVecGetArrayF90(grid,field%work,vec_ptr,ierr)

        cur_connection_set => cur_source_sink%connection_set
    
        do iconn = 1, cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          select case(option%iflowmode)
            case(RICHARDS_MODE,G_MODE)
               call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                            x_width,y_width,z_width, &
                                            x_count,y_count,z_count, &
                                            ghosted_neighbors,option)
               ! ghosted neighbors is ordered first in x, then, y, then z
               icount = 0
               sum = 0.d0
               ! x-direction
               do while (icount < x_count)
                 icount = icount + 1
                 neighbor_ghosted_id = ghosted_neighbors(icount)
                 sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                             grid%structured_grid%dy(neighbor_ghosted_id)* &
                             grid%structured_grid%dz(neighbor_ghosted_id)
                 
               enddo
               ! y-direction
               do while (icount < x_count + y_count)
                 icount = icount + 1
                 neighbor_ghosted_id = ghosted_neighbors(icount)                 
                 sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                             grid%structured_grid%dx(neighbor_ghosted_id)* &
                             grid%structured_grid%dz(neighbor_ghosted_id)
                 
               enddo
               ! z-direction
               do while (icount < x_count + y_count + z_count)
                 icount = icount + 1
                 neighbor_ghosted_id = ghosted_neighbors(icount)                 
                 sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                             grid%structured_grid%dx(neighbor_ghosted_id)* &
                             grid%structured_grid%dy(neighbor_ghosted_id)
               enddo
               vec_ptr(local_id) = vec_ptr(local_id) + sum
            case(TH_MODE)
            case(THC_MODE)
            case(THMC_MODE)
            case(MPH_MODE)
            case(IMS_MODE)
            case(MIS_MODE)
            case(FLASH2_MODE)
          end select

        enddo
        
        call GridVecRestoreArrayF90(grid,field%work,vec_ptr,ierr)
        call VecNorm(field%work,NORM_1,scale,ierr)
        scale = 1.d0/scale
        call VecScale(field%work,scale,ierr)

        call GridVecGetArrayF90(grid,field%work,vec_ptr, ierr)
        do iconn = 1, cur_connection_set%num_connections      
          local_id = cur_connection_set%id_dn(iconn)
          select case(option%iflowmode)
            case(RICHARDS_MODE,G_MODE)
              cur_source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
                vec_ptr(local_id)
            case(TH_MODE)
            case(THC_MODE)
            case(THMC_MODE)
            case(MPH_MODE)
            case(IMS_MODE)
            case(MIS_MODE)
            case(FLASH2_MODE)
          end select

        enddo
        call GridVecRestoreArrayF90(grid,field%work,vec_ptr,ierr)
        
        cur_source_sink => cur_source_sink%next
      enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_loc_ptr, ierr)
   
end subroutine CondControlScaleSourceSink

! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine CondControlAssignFlowInitCondSurface(surf_realization)

  use Surface_Realization_class
  use Discretization_module
  use Region_module
  use Option_module
  use Surface_Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  use Level_module
  use Patch_module
  use Water_EOS_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(surface_realization_type) :: surf_realization
  
  PetscInt :: icell, iconn, idof, iface
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)!, iphase_loc_p(:), xx_faces_p(:)
  PetscErrorCode :: ierr
  
  PetscReal :: temperature, p_sat
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(flow_general_condition_type), pointer :: general

  option => surf_realization%option
  discretization => surf_realization%discretization
  surf_field => surf_realization%surf_field
  patch => surf_realization%patch


  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid

      select case(option%iflowmode)
      
        case(G_MODE) ! general phase mode

        case (RICHARDS_MODE,TH_MODE)
          ! assign initial conditions values to domain
          call GridVecGetArrayF90(grid,surf_field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
    
          xx_p = -999.d0
      
          initial_condition => cur_patch%initial_conditions%first
          do
      
            if (.not.associated(initial_condition)) exit

              if (.not.associated(initial_condition%flow_aux_real_var)) then
                if (.not.associated(initial_condition%flow_condition)) then
                  option%io_buffer = 'Flow condition is NULL in initial condition'
                  call printErrMsg(option)
                endif
                do icell=1,initial_condition%region%num_cells
                  local_id = initial_condition%region%cell_ids(icell)
                  ghosted_id = grid%nL2G(local_id)
                  iend = local_id*option%nflowdof
                  ibegin = iend-option%nflowdof+1
                  if (cur_patch%imat(ghosted_id) <= 0) then
                    xx_p(ibegin:iend) = 0.d0
                    cycle
                  endif
                  do idof = 1, option%nflowdof
                    xx_p(ibegin+idof-1) = &
                          initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%flow_dataset%time_series%cur_value(1)
                    !TODO(GB): Correct the initialization of surface flow condition
                    if (idof == 1) xx_p(ibegin+idof-1) = 0.d0
                    !if (idof == 1.and.option%iflowmode==TH_MODE) then
                    !  xx_p(ibegin+idof-1) = 0.d0+option%reference_pressure
                    !endif
                  enddo
                enddo
              else
                do iconn=1,initial_condition%connection_set%num_connections
                  local_id = initial_condition%connection_set%id_dn(iconn)
                  ghosted_id = grid%nL2G(local_id)
                  iend = local_id*option%nflowdof
                  ibegin = iend-option%nflowdof+1
                  if (cur_patch%imat(ghosted_id) <= 0) then
                    xx_p(ibegin:iend) = 0.d0
                   cycle
                  endif
                  xx_p(ibegin:iend) = &
                        initial_condition%flow_aux_real_var(1:option%nflowdof,iconn)
                  !TODO(gb): Correct the initialization of surface flow condition
                  xx_p(ibegin) = 0.d0
                  !if (idof == 1.and.option%iflowmode==TH_MODE) then
                  !  xx_p(ibegin) = 0.d0 + option%reference_pressure
                  !endif
                enddo
              endif
            initial_condition => initial_condition%next
          enddo
     
          call GridVecRestoreArrayF90(grid,surf_field%flow_xx,xx_p, ierr)
        case default
          option%io_buffer = 'CondControlAssignFlowInitCondSurface not ' // &
            'for this mode'
          call printErrMsg(option)
      end select 
   
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

end subroutine CondControlAssignFlowInitCondSurface
#endif
! SURFACE_FLOW

end module Condition_Control_module
