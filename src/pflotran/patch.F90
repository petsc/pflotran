module Patch_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Breakthrough_module
  use Strata_module
  use Region_module
  use Material_module
  
  use Auxilliary_module

  implicit none

  private

#include "definitions.h"

  type, public :: patch_type 
    
    PetscInt :: id
    
    ! thiese arrays will be used by all modes, mode-specific arrays should
    ! go in the auxilliary data stucture for that mode
    PetscInt, pointer :: imat(:)
    type(material_ptr_type), pointer :: material_array(:)
    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    
    type(grid_type), pointer :: grid

    type(region_list_type), pointer :: regions

    type(coupler_list_type), pointer :: boundary_conditions
    type(coupler_list_type), pointer :: initial_conditions
    type(coupler_list_type), pointer :: source_sinks

    type(strata_list_type), pointer :: strata
    type(breakthrough_list_type), pointer :: breakthrough
    
    type(auxilliary_type) :: aux
    
    type(patch_type), pointer :: next

  end type patch_type

  ! pointer data structure required for making an array of patch pointers in F90
  type, public :: patch_ptr_type
    type(patch_type), pointer :: ptr           ! pointer to the patch_type
  end type patch_ptr_type 

  type, public :: patch_list_type
    PetscInt :: num_patch_objects
    type(patch_type), pointer :: first
    type(patch_type), pointer :: last
    type(patch_ptr_type), pointer :: array(:)
  end type patch_list_type
    
  public :: PatchCreate, PatchDestroy, PatchCreateList, PatchDestroyList, &
            PatchAddToList, PatchConvertListToArray, PatchProcessCouplers, &
            PatchUpdateAllCouplerAuxVars, PatchInitAllCouplerAuxVars, &
            PatchLocalizeRegions, PatchAssignUniformVelocity, &
            PatchGetDataset, PatchGetDatasetValueAtCell, &
            PatchSetDataset, &
            PatchInitConstraints

contains

! ************************************************************************** !
!
! PatchCreate: Allocates and initializes a new Patch object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function PatchCreate()

  implicit none
  
  type(patch_type), pointer :: PatchCreate
  
  type(patch_type), pointer :: patch
  
  allocate(patch)

  patch%id = 0
  nullify(patch%imat)
  nullify(patch%material_array)
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)

  nullify(patch%grid)

  allocate(patch%regions)
  call RegionInitList(patch%regions)
  
  allocate(patch%boundary_conditions)
  call CouplerInitList(patch%boundary_conditions)
  allocate(patch%initial_conditions)
  call CouplerInitList(patch%initial_conditions)
  allocate(patch%source_sinks)
  call CouplerInitList(patch%source_sinks)

  allocate(patch%breakthrough)
  call BreakthroughInitList(patch%breakthrough)

  allocate(patch%strata)
  call StrataInitList(patch%strata)
  
  call AuxInit(patch%aux)
  
  nullify(patch%next)
  
  PatchCreate => patch
  
end function PatchCreate

! ************************************************************************** !
!
! PatchListCreate: Creates a patch list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function PatchCreateList()

  implicit none

  type(patch_list_type), pointer :: PatchCreateList

  type(patch_list_type), pointer :: patch_list
  
  allocate(patch_list)
  nullify(patch_list%first)
  nullify(patch_list%last)
  nullify(patch_list%array)
  patch_list%num_patch_objects = 0

  PatchCreateList => patch_list

end function PatchCreateList

! ************************************************************************** !
!
! PatchAddToList: Adds a new patch to list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchAddToList(new_patch,patch_list)

  implicit none
  
  type(patch_type), pointer :: new_patch
  type(patch_list_type) :: patch_list
  
  if(associated(new_patch)) then
     patch_list%num_patch_objects = patch_list%num_patch_objects + 1
     new_patch%id = patch_list%num_patch_objects
     if (.not.associated(patch_list%first)) patch_list%first => new_patch
     if (associated(patch_list%last)) patch_list%last%next => new_patch
     patch_list%last => new_patch
  end if
end subroutine PatchAddToList

! ************************************************************************** !
!
! PatchConvertListToArray: Creates an array of pointers to the 
!                               patchs in the patch list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchConvertListToArray(patch_list)

  implicit none
  
  type(patch_list_type) :: patch_list
    
  PetscInt :: count
  type(patch_type), pointer :: cur_patch
  
  
  allocate(patch_list%array(patch_list%num_patch_objects))
  
  cur_patch => patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    patch_list%array(cur_patch%id)%ptr => cur_patch
    cur_patch => cur_patch%next
  enddo

end subroutine PatchConvertListToArray

! ************************************************************************** !
!
! PatchLocalizeRegions: Localizes regions within each patch
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchLocalizeRegions(patch,regions,option)

  use Option_module
  use Region_module

  implicit none
  
  type(patch_type) :: patch
  type(region_list_type) :: regions
  type(option_type) :: option
  
  type(region_type), pointer :: cur_region
  type(region_type), pointer :: patch_region
  
  cur_region => regions%first
  do
    if (.not.associated(cur_region)) exit
    patch_region => RegionCreate(cur_region)
    call RegionAddToList(patch_region,patch%regions)
    cur_region => cur_region%next
  enddo
  
  call GridLocalizeRegions(patch%grid,patch%regions,option)
 
end subroutine PatchLocalizeRegions

! ************************************************************************** !
!
! PatchProcessCouplers: Assigns conditions and regions to couplers
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
                                materials,option)

  use Option_module
  use Material_module
  use Condition_module
  use Connection_module

  implicit none
  
  type(patch_type) :: patch
  type(material_type), pointer :: materials
  type(condition_list_type) :: flow_conditions
  type(tran_condition_list_type) :: transport_conditions
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list 
  type(strata_type), pointer :: strata
  type(breakthrough_type), pointer :: breakthrough, next_breakthrough
  
  PetscInt :: temp_int
  
  call MaterialConvertListToArray(materials,patch%material_array)
  
  ! boundary conditions
  coupler => patch%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in boundary condition list'
      call printErrMsg(option,string)
    endif
    ! pointer to flow condition
    if (len_trim(coupler%flow_condition_name) > 0) then
      coupler%flow_condition => &
        FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in boundary condition list'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        TranConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in boundary condition list'
        call printErrMsg(option,string)
      endif
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => patch%initial_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in initial condition list'
      call printErrMsg(option,string)
    endif
    ! pointer to flow condition
    if (len_trim(coupler%flow_condition_name) > 0) then
      coupler%flow_condition => &
        FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in initial condition list'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        TranConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in initial condition list'
        call printErrMsg(option,string)
      endif
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => patch%source_sinks%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in source/sink list'
      call printErrMsg(option,string)
    endif
    ! pointer to flow condition
    if (len_trim(coupler%flow_condition_name) > 0) then
      coupler%flow_condition => &
        FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in source/sink condition list'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        TranConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in source/sink condition list'
        call printErrMsg(option,string)
      endif
    endif
    coupler => coupler%next
  enddo

!----------------------------  
! AUX  
    
  ! strata
  ! connect pointers from strata to regions
  strata => patch%strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 1) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  patch%regions)
      if (.not.associated(strata%region)) then
        string = 'Region ' // trim(strata%region_name) // &
                 ' not found in region list'
        call printErrMsg(option,string)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material => MaterialGetPtrFromArray(strata%material_name, &
                                                   patch%material_array)
        if (.not.associated(strata%material)) then
          string = 'Material ' // trim(strata%material_name) // &
                   ' not found in material list'
          call printErrMsg(option,string)
        endif
      endif
    else
      nullify(strata%region)
      nullify(strata%material)
    endif
    strata => strata%next
  enddo 

  ! breakthrough
  breakthrough => patch%breakthrough%first
  do
    if (.not.associated(breakthrough)) exit
    next_breakthrough => breakthrough%next
    ! pointer to region
    breakthrough%region => RegionGetPtrFromList(breakthrough%region_name, &
                                                patch%regions)
    if (.not.associated(breakthrough%region)) then
      string = 'Region ' // trim(breakthrough%region_name) // &
               ' not found in region list'
      call printErrMsg(option,string)
    endif
    if (breakthrough%region%num_cells == 0) then
      ! remove the breakthrough object
      call BreakthroughRemoveFromList(breakthrough,patch%breakthrough)
    endif
    breakthrough => next_breakthrough
  enddo
 
  ! connectivity between initial conditions, boundary conditions, srcs/sinks, etc and grid
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%initial_conditions)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%boundary_conditions)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%source_sinks)
                                     
  allocate(patch%internal_velocities(option%nphase, &
           ConnectionGetNumberInList(patch%grid%internal_connection_set_list)))
  patch%internal_velocities = 0.d0
  temp_int = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  allocate(patch%boundary_velocities(option%nphase,temp_int)) 
  patch%boundary_velocities = 0.d0          

end subroutine PatchProcessCouplers

! ************************************************************************** !
!
! PatchInitAllCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchInitAllCouplerAuxVars(patch,reaction,option)

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  type(patch_type) :: patch
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscTruth :: force_update_flag = PETSC_TRUE
  
  call PatchInitCouplerAuxVars(patch%initial_conditions,reaction, &
                               option)
  call PatchInitCouplerAuxVars(patch%boundary_conditions,reaction, &
                               option)
  call PatchInitCouplerAuxVars(patch%source_sinks,reaction, &
                               option)
  
  call PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

end subroutine PatchInitAllCouplerAuxVars

! ************************************************************************** !
!
! PatchInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchInitCouplerAuxVars(coupler_list,reaction,option)

  use Option_module
  use Connection_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Condition_module
  
  implicit none
  
  type(coupler_list_type), pointer :: coupler_list
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: num_connections
  PetscTruth :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  
  if (.not.associated(coupler_list)) return
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%connection_set)) then
      num_connections = coupler%connection_set%num_connections

      ! FLOW
      if (associated(coupler%flow_condition) .and. &
          (coupler%itype == INITIAL_COUPLER_TYPE .or. &
           coupler%itype == BOUNDARY_COUPLER_TYPE)) then

        if (associated(coupler%flow_condition%pressure)) then

          ! allocate arrays that match the number of connections
          select case(option%iflowmode)

            case(THC_MODE,RICHARDS_MODE)
           
              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0

            case(MPH_MODE)

              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
                
            case default
          end select
      
        endif
      
      endif
      
    endif

    ! TRANSPORT   
    if (associated(coupler%tran_condition)) then
      cur_constraint_coupler => &
        coupler%tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        allocate(cur_constraint_coupler%global_auxvar)
        allocate(cur_constraint_coupler%rt_auxvar)
        call GlobalAuxVarInit(cur_constraint_coupler%global_auxvar,option)
        call RTAuxVarInit(cur_constraint_coupler%rt_auxvar,reaction,option)
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
      
    coupler => coupler%next
  enddo
  
end subroutine PatchInitCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateAllCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

  use Option_module
  
  implicit none
  
  type(patch_type) :: patch
  PetscTruth :: force_update_flag
  type(option_type) :: option
  
  call PatchUpdateCouplerAuxVars(patch,patch%initial_conditions, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%boundary_conditions, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%source_sinks, &
                                 force_update_flag,option)

end subroutine PatchUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in list
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVars(patch,coupler_list,force_update_flag, &
                                     option)

  use Option_module
  use Condition_module
  use Hydrostatic_module

  implicit none
  
  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  PetscTruth :: force_update_flag
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  PetscTruth :: update
  
  PetscInt :: idof, num_connections
  
  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first
  
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
        
      num_connections = coupler%connection_set%num_connections

      flow_condition => coupler%flow_condition

      update = PETSC_FALSE
      select case(option%iflowmode)
        case(THC_MODE,MPH_MODE)
          if (force_update_flag .or. &
              flow_condition%pressure%dataset%is_transient .or. &
              flow_condition%pressure%gradient%is_transient .or. &
              flow_condition%pressure%datum%is_transient .or. &
              flow_condition%temperature%dataset%is_transient .or. &
              flow_condition%temperature%gradient%is_transient .or. &
              flow_condition%temperature%datum%is_transient .or. &
              flow_condition%concentration%dataset%is_transient .or. &
              flow_condition%concentration%gradient%is_transient .or. &
              flow_condition%concentration%datum%is_transient) then
            update = PETSC_TRUE
          endif
        case(RICHARDS_MODE)
          if (force_update_flag .or. &
              flow_condition%pressure%dataset%is_transient .or. &
              flow_condition%pressure%gradient%is_transient .or. &
              flow_condition%pressure%datum%is_transient) then
            update = PETSC_TRUE
          endif
      end select
      
      if (update) then
        if (associated(flow_condition%pressure)) then
          select case(flow_condition%pressure%itype)
            case(DIRICHLET_BC,NEUMANN_BC,MASS_RATE_SS,ZERO_GRADIENT_BC, &
                 VOLUMETRIC_RATE_SS)
!              do idof = 1, condition%num_sub_conditions
!                if (associated(condition%sub_condition_ptr(idof)%ptr)) then
!                  coupler%flow_aux_real_var(idof,1:num_connections) = &
!                    condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
!                endif
!              enddo
              select case(option%iflowmode)
                case(MPH_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)  ! <-- Chuan Fix
                   coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(THC_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)
                  coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(RICHARDS_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)
              end select
            case(HYDROSTATIC_BC,SEEPAGE_BC)
              call HydrostaticUpdateCoupler(coupler,option,patch%grid)
          end select
        endif
      endif
     
    endif
      
    ! TRANSPORT
    ! nothing for transport at this point in time

    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !
!
! PatchInitConstraints: Initializes constraint concentrations
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine PatchInitConstraints(patch,reaction,option)

  use Reaction_Aux_module
    
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  type(reaction_type), pointer :: reaction

  call PatchInitCouplerConstraints(patch%initial_conditions, &
                                   reaction,option)
  call PatchInitCouplerConstraints(patch%boundary_conditions, &
                                   reaction,option)
  call PatchInitCouplerConstraints(patch%source_sinks, &
                                   reaction,option)

end subroutine PatchInitConstraints

! ************************************************************************** !
!
! PatchInitCouplerConstraints: Initializes constraint concentrations
!                              for a given coupler
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine PatchInitCouplerConstraints(coupler_list,reaction,option)

  use Reaction_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Condition_module
  use water_eos_module
    
  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(option_type) :: option
  type(reaction_type), pointer :: reaction

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(coupler_type), pointer :: cur_coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  PetscReal :: r1, r2, r3, r4, r5, r6
  PetscErrorCode :: ierr
  
  cur_coupler => coupler_list%first
  do
    if (.not.associated(cur_coupler)) exit

    cur_constraint_coupler => &
      cur_coupler%tran_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      global_auxvar => cur_constraint_coupler%global_auxvar
      rt_auxvar => cur_constraint_coupler%rt_auxvar
      if (associated(cur_coupler%flow_condition)) then
        if (associated(cur_coupler%flow_condition%pressure)) then
          global_auxvar%pres = &
            cur_coupler%flow_condition%pressure%dataset%cur_value(1)
        else
          global_auxvar%pres = option%reference_pressure
        endif
        if (associated(cur_coupler%flow_condition%temperature)) then
          global_auxvar%temp = &
            cur_coupler%flow_condition%temperature%dataset%cur_value(1)
        else
          global_auxvar%temp = option%reference_temperature
        endif
        call wateos(global_auxvar%temp(1),global_auxvar%pres(1), &
                    global_auxvar%den_kg(1),r1,r2,r3,r4,r5,r6, &
                    option%scale,ierr) 
      else
        global_auxvar%pres = option%reference_pressure
        global_auxvar%temp = option%reference_temperature
        global_auxvar%den_kg = option%reference_density
      endif     
      global_auxvar%sat = option%reference_saturation  
      call ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                            reaction,cur_constraint_coupler%constraint_name, &
                            cur_constraint_coupler%aqueous_species, &
                            cur_constraint_coupler%num_iterations,option)
      ! turn on flag indicating constraint has not yet been used
      cur_constraint_coupler%iflag = ONE_INTEGER
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    cur_coupler => cur_coupler%next
  enddo

end subroutine PatchInitCouplerConstraints

! ************************************************************************** !
!
! PatchAssignUniformVelocity: Assigns uniform velocity in connection list
!                        darcy velocities
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine PatchAssignUniformVelocity(patch,option)

  use Option_module
  use Coupler_module
  use Condition_module
  use Connection_module
  
  implicit none
  
  type(patch_type), pointer :: patch   
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: vdarcy

  grid => patch%grid
    
  ! Internal Flux Terms -----------------------------------
  cur_connection_set => grid%internal_connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = OptionDotProduct(option%uniform_velocity, &
                                cur_connection_set%dist(1:3,iconn))
      patch%internal_velocities(1,sum_connection) = vdarcy
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
      vdarcy = OptionDotProduct(option%uniform_velocity, &
                                cur_connection_set%dist(1:3,iconn))
      patch%boundary_velocities(1,sum_connection) = vdarcy
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PatchAssignUniformVelocity

! ************************************************************************** !
!
! PatchAuxVarsUpToDate: Checks to see if aux vars are up to date
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function PatchAuxVarsUpToDate(patch)

  use Grid_module
  use Option_module
  use Field_module
  
  use Mphase_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module  
  
  type(patch_type) :: patch
  
  PetscTruth :: PatchAuxVarsUpToDate
  PetscTruth :: flow_up_to_date
  PetscTruth :: transport_up_to_date
  
  if (associated(patch%aux%THC)) then
    flow_up_to_date = patch%aux%THC%aux_vars_up_to_date
  else if (associated(patch%aux%Richards)) then
    flow_up_to_date = patch%aux%Richards%aux_vars_up_to_date
  else if (associated(patch%aux%Mphase)) then
    flow_up_to_date = patch%aux%Mphase%aux_vars_up_to_date
  endif

  if (associated(patch%aux%RT)) then
    transport_up_to_date = patch%aux%RT%aux_vars_up_to_date
  endif
  
  PatchAuxVarsUpToDate = flow_up_to_date .or. transport_up_to_date
  
end function PatchAuxVarsUpToDate

! ************************************************************************** !
!
! PatchGetDataset: Extracts variables indexed by ivar and isubvar from a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchGetDataset(patch,field,option,vec,ivar,isubvar)

  use Grid_module
  use Option_module
  use Field_module
  
  use Mphase_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module  
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid

  call GridVecGetArrayF90(grid,vec,vec_ptr,ierr)

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%temp
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%pres
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_kg
            enddo
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = 0.d0
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%u
            enddo
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%temp
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%pres
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(1)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(2)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(2)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1)
            enddo
        end select
      endif
      
    case(PH,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY,TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_RATE,MINERAL_VOLUME_FRACTION,SURFACE_CMPLX)
         
      select case(ivar)
        case(PH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              -log10(patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)* &
                     patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))
          enddo
        case(PRIMARY_MOLALITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%pri_molal(isubvar)
          enddo
        case(PRIMARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar)* &
              (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(SECONDARY_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)
          enddo
        case(SECONDARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)* &
              (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(TOTAL_MOLALITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)/ &
              (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(TOTAL_MOLARITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%total(isubvar,iphase)
          enddo
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_volfrac(isubvar)
          enddo
        case(MINERAL_RATE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_rate(isubvar)
          enddo
        case(SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%eqsurfcmplx_conc(isubvar)
          enddo
      end select
    case(PHASE)
      call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
    case default
      call printErrMsg(option,'IVAR not found in OutputGetVarFromArray')      
  end select

  call GridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)
  
end subroutine PatchGetDataset

! ************************************************************************** !
!
! PatchGetDatasetValueAtCell: Returns variables indexed by ivar,
!                             isubvar, local id from Reactive Transport type
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function PatchGetDatasetValueAtCell(patch,field,option,ivar,isubvar, &
                                    ghosted_id)

  use Grid_module
  use Option_module
  use Field_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscReal :: PatchGetDatasetValueAtCell
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase
  PetscInt :: ghosted_id

  PetscReal :: value
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr2(:)  
  PetscErrorCode :: ierr

  grid => patch%grid
  
  value = -999.d0

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%THC%aux_vars(ghosted_id)%temp
          case(PRESSURE)
            value = patch%aux%THC%aux_vars(ghosted_id)%pres
          case(LIQUID_SATURATION)
            value = patch%aux%THC%aux_vars(ghosted_id)%sat
          case(LIQUID_DENSITY)
            value = patch%aux%THC%aux_vars(ghosted_id)%den_kg
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            value = 0.d0
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%THC%aux_vars(ghosted_id)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%THC%aux_vars(ghosted_id)%u
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%temp
          case(PRESSURE)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%pres
          case(LIQUID_SATURATION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(1)
          case(GAS_SATURATION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%sat(2)
          case(GAS_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar)
          case(GAS_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(2)
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
        end select
      endif
      
    case(PH,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY, &
         TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_VOLUME_FRACTION,MINERAL_RATE,SURFACE_CMPLX)
         
      select case(ivar)
        case(PH)
          value = -log10(patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)* &
                         patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))
        case(PRIMARY_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar)
        case(PRIMARY_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar)* &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0
        case(SECONDARY_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)
        case(SECONDARY_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)* &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0
        case(TOTAL_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)/ &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)*1000.d0
        case(TOTAL_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)
        case(MINERAL_VOLUME_FRACTION)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(isubvar)
        case(MINERAL_RATE)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_rate(isubvar)
        case(SURFACE_CMPLX)
          value = patch%aux%RT%aux_vars(ghosted_id)%eqsurfcmplx_conc(isubvar)
      end select
    case(PHASE)
      call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      value = patch%imat(ghosted_id)
  end select

  PatchGetDatasetValueAtCell = value
 
end function PatchGetDatasetValueAtCell

! ************************************************************************** !
!
! PatchSetDataset: Sets variables indexed by ivar and isubvar within a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchSetDataset(patch,field,option,vec,ivar,isubvar)

  use Grid_module
  use Option_module
  use Field_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase

  PetscInt :: local_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid

  call GridVecGetArrayF90(grid,vec,vec_ptr,ierr)

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY)
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%temp = vec_ptr(local_id)
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%pres = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_kg = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%xmol(isubvar) = vec_ptr(local_id)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%THC%aux_vars(grid%nL2G(local_id))%u = vec_ptr(local_id)
            enddo
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1) = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1) = vec_ptr(local_id)
            enddo
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%temp = vec_ptr(local_id)
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%pres = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(1) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(1) = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(2) = vec_ptr(local_id)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar) = vec_ptr(local_id)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2) = vec_ptr(local_id)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(2) = vec_ptr(local_id)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar) = vec_ptr(local_id)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1) = vec_ptr(local_id)
            enddo
        end select
      endif
    case(PRIMARY_MOLALITY,TOTAL_MOLARITY,MINERAL_VOLUME_FRACTION)
      select case(ivar)
        case(PRIMARY_MOLALITY)
          do local_id=1,grid%nlmax
            patch%aux%RT%aux_vars(grid%nL2G(local_id))%pri_molal(isubvar) = vec_ptr(local_id)
          enddo
        case(TOTAL_MOLARITY)
          do local_id=1,grid%nlmax
            patch%aux%RT%aux_vars(grid%nL2G(local_id))%total(isubvar,iphase) = vec_ptr(local_id)
          enddo
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_volfrac(isubvar) = vec_ptr(local_id)
          enddo
      end select
    case(PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY,TOTAL_MOLALITY)
      select case(ivar)
        case(PRIMARY_MOLARITY)
          call printErrMsg(option,'Setting of primary molarity at grid cell not supported.')
        case(SECONDARY_MOLALITY)
          call printErrMsg(option,'Setting of secondary molality at grid cell not supported.')
        case(SECONDARY_MOLARITY)
          call printErrMsg(option,'Setting of secondary molarity at grid cell not supported.')
        case(TOTAL_MOLALITY)
          call printErrMsg(option,'Setting of total molality at grid cell not supported.')
      end select
    case(PHASE)
      call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
      enddo
      call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        patch%imat(grid%nL2G(local_id)) = vec_ptr(local_id)
      enddo
  end select

  call GridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)
  
end subroutine PatchSetDataset

! ************************************************************************** !
!
! PatchDestroyList: Deallocates a patch list and array of patches
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine PatchDestroyList(patch_list)

  implicit none
  
  type(patch_list_type), pointer :: patch_list
    
  type(patch_type), pointer :: cur_patch, prev_patch
  
  if (.not.associated(patch_list)) return
  
  if (associated(patch_list%array)) deallocate(patch_list%array)
  nullify(patch_list%array)
  
  cur_patch => patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    prev_patch => cur_patch
    cur_patch => cur_patch%next
    call PatchDestroy(prev_patch)
  enddo
  
  nullify(patch_list%first)
  nullify(patch_list%last)
  patch_list%num_patch_objects = 0
  
  deallocate(patch_list)
  nullify(patch_list)

end subroutine PatchDestroyList

! ************************************************************************** !
!
! PatchDestroy: Deallocates a patch object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchDestroy(patch)

  implicit none
  
  type(patch_type), pointer :: patch
  
  if (associated(patch%imat)) deallocate(patch%imat)
  nullify(patch%imat)
  if (associated(patch%internal_velocities)) deallocate(patch%internal_velocities)
  nullify(patch%internal_velocities)
  if (associated(patch%boundary_velocities)) deallocate(patch%boundary_velocities)
  nullify(patch%boundary_velocities)

  call GridDestroy(patch%grid)
  call RegionDestroyList(patch%regions)
  call CouplerDestroyList(patch%boundary_conditions)
  call CouplerDestroyList(patch%initial_conditions)
  call CouplerDestroyList(patch%source_sinks)
  
  call BreakthroughDestroyList(patch%breakthrough)
  call StrataDestroyList(patch%strata)
  
  call AuxDestroy(patch%aux)
  
  call BreakthroughDestroyList(patch%breakthrough)
  
  deallocate(patch)
  nullify(patch)
  
end subroutine PatchDestroy

end module Patch_module
