module Patch_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Breakthrough_module
  use Strata_module
  use Region_module
  
  use Auxilliary_module

  implicit none

  private

#include "definitions.h"

  type, public :: patch_type 
    
    PetscInt :: id
    
    ! thiese arrays will be used by all modes, mode-specific arrays should
    ! go in the auxilliary data stucture for that mode
    PetscInt, pointer :: imat(:)
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
            PatchLocalizeRegions, PatchAssignUniformVelocity

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
  
  patch_list%num_patch_objects = patch_list%num_patch_objects + 1
  new_patch%id = patch_list%num_patch_objects
  if (.not.associated(patch_list%first)) patch_list%first => new_patch
  if (associated(patch_list%last)) patch_list%last%next => new_patch
  patch_list%last => new_patch
  
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
  type(condition_list_type) :: transport_conditions
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list 
  type(strata_type), pointer :: strata
  type(breakthrough_type), pointer :: breakthrough, next_breakthrough
  
  PetscInt :: temp_int
  
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
        ConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in boundary condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%flow_condition%iclass /= FLOW_CLASS) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' is not of the FLOW class'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        ConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in boundary condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%tran_condition%iclass /= TRANSPORT_CLASS) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' is not of the TRAN class'
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
        ConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in initial condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%flow_condition%iclass /= FLOW_CLASS) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' is not of the FLOW class'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        ConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in initial condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%tran_condition%iclass /= TRANSPORT_CLASS) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' is not of the TRAN class'
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
        ConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
      if (.not.associated(coupler%flow_condition)) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' not found in source/sink condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%flow_condition%iclass /= FLOW_CLASS) then
        string = 'Condition ' // trim(coupler%flow_condition_name) // &
                 ' is not of the FLOW class'
        call printErrMsg(option,string)
      endif
    endif
    ! pointer to transport condition
    if (len_trim(coupler%tran_condition_name) > 0) then
      coupler%tran_condition => &
        ConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
      if (.not.associated(coupler%tran_condition)) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' not found in source/sink condition list'
        call printErrMsg(option,string)
      endif
      if (coupler%tran_condition%iclass /= TRANSPORT_CLASS) then
        string = 'Condition ' // trim(coupler%tran_condition_name) // &
                 ' is not of the TRAN class'
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
        strata%material => MaterialGetPtrFromList(strata%material_name, &
                                                  materials)
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
subroutine PatchInitAllCouplerAuxVars(patch,option)

  use Option_module
  
  implicit none
  
  type(patch_type) :: patch
  type(option_type) :: option
  
  logical :: force_update_flag = .true.
  
  call PatchInitCouplerAuxVars(patch,patch%boundary_conditions,option)
  call PatchInitCouplerAuxVars(patch,patch%initial_conditions,option)
  
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
subroutine PatchInitCouplerAuxVars(patch,coupler_list,option)

  use Option_module
  use Connection_module
  
  implicit none
  
  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  type(option_type) :: option
  
  PetscInt :: num_connections
  logical :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  
  if (.not.associated(coupler_list)) return
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%connection_set)) then
      num_connections = coupler%connection_set%num_connections

      ! FLOW
      if (associated(coupler%flow_condition)) then

        if (associated(coupler%flow_condition%pressure)) then

          ! allocate arrays that match the number of connections
          select case(option%iflowmode)

            case(RICHARDS_MODE,RICHARDS_LITE_MODE)
           
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
      
      ! TRANSPORT   
      if (associated(coupler%tran_condition)) then

        allocate(coupler%tran_aux_real_var(option%ntrandof,num_connections))
        coupler%tran_aux_real_var = 0.d0

      endif
      
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
  logical :: force_update_flag
  type(option_type) :: option
  
  call PatchUpdateCouplerAuxVars(patch,patch%boundary_conditions, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%initial_conditions, &
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
  logical :: force_update_flag
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(condition_type), pointer :: condition
  logical :: update
  
  PetscInt :: idof, num_connections
  
  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first
  
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
        
      num_connections = coupler%connection_set%num_connections

      condition => coupler%flow_condition

      update = .false.
      select case(option%iflowmode)
        case(RICHARDS_MODE,MPH_MODE)
          if (force_update_flag .or. &
              condition%pressure%dataset%is_transient .or. &
              condition%pressure%gradient%is_transient .or. &
              condition%pressure%datum%is_transient .or. &
              condition%temperature%dataset%is_transient .or. &
              condition%temperature%gradient%is_transient .or. &
              condition%temperature%datum%is_transient .or. &
              condition%concentration%dataset%is_transient .or. &
              condition%concentration%gradient%is_transient .or. &
              condition%concentration%datum%is_transient) then
            update = .true.
          endif
        case(RICHARDS_LITE_MODE)
          if (force_update_flag .or. &
              condition%pressure%dataset%is_transient .or. &
              condition%pressure%gradient%is_transient .or. &
              condition%pressure%datum%is_transient) then
            update = .true.
          endif
      end select
      
      if (update) then
        if (associated(condition%pressure)) then
          select case(condition%pressure%itype)
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
                    condition%pressure%dataset%cur_value(1)  ! <-- Chuan Fix
                  coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    condition%temperature%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    condition%concentration%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    condition%iphase
                case(RICHARDS_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    condition%pressure%dataset%cur_value(1)
                  coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    condition%temperature%dataset%cur_value(1)
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    condition%concentration%dataset%cur_value(1)
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    condition%iphase
                case(RICHARDS_LITE_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    condition%pressure%dataset%cur_value(1)
              end select
            case(HYDROSTATIC_BC,SEEPAGE_BC)
    !          call HydrostaticUpdateCoupler(coupler,patch%option,patch%grid)
              call HydrostaticUpdateCouplerBetter(coupler,option,patch%grid)
          end select
        endif
      endif
     
    endif
      
    ! TRANSPORT
    if (associated(coupler%tran_aux_real_var)) then
    
      num_connections = coupler%connection_set%num_connections

      condition => coupler%tran_condition

      update = .false.
      if (force_update_flag) then
        update = .true.
      else 
        do idof = 1, condition%num_sub_conditions
          if (condition%sub_condition_ptr(idof)%ptr%dataset%is_transient .or. &
              condition%sub_condition_ptr(idof)%ptr%gradient%is_transient .or. &
              condition%sub_condition_ptr(idof)%ptr%datum%is_transient) then
            update = .true.
            exit
          endif
        enddo
      endif
      
      if (update) then ! for now, everything transport is dirichlet-type
        do idof = 1, condition%num_sub_conditions
          coupler%tran_aux_real_var(idof,1:num_connections) = &
            condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
        enddo
      endif
      
    endif

    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

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
