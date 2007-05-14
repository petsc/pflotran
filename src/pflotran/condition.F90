module Condition_module

 use pflow_gridtype_module

 implicit none
 
 public
 save

#include "definitions.h"

  type condition_type
    real*8, pointer :: times(:)
    real*8, pointer :: values(:)
    real*8 :: cur_value
    real*8 :: datum
    integer :: id
    integer :: itype
    character(len=MAXWORDLENGTH) :: ctype
    integer :: cur_time_index, max_time_index
    type(condition_type), pointer :: next
  end type condition_type

  type condition_type_ptr
    type(condition_type), pointer :: ptr
  end type condition_type_ptr

  type(condition_type), pointer :: condition_list
  type(condition_type_ptr), pointer :: condition_array(:)

  integer, private :: num_conditions

  type(condition_type), pointer, private :: last_condition

  public InitializeConditionList, AddCondition, UpdateConditions, &
         UpdateBoundaryConditions
  
contains
  
! ************************************************************************** !
!
! InitializeConditionList: Initializes the condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine InitializeConditionList()

  implicit none

  nullify(condition_list)
  nullify(last_condition)
  num_conditions = 0

end subroutine InitializeConditionList

! ************************************************************************** !
!
! AddCondition: Add condition to condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine AddCondition(new_condition)

  implicit none

  type(condition_type), pointer :: new_condition

  new_condition%cur_time_index = 1

  num_conditions = num_conditions + 1
  new_condition%id =  num_conditions
  if (.not.associated(condition_list)) condition_list => new_condition
  if (associated(last_condition)) last_condition%next => new_condition
  last_condition => new_condition

end subroutine AddCondition

! ************************************************************************** !
!
! ConvertListToArray: Converts the linked list of conditions to an array
!                     of conditions
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine ConvertListToArray()

  implicit none

  type(condition_type), pointer :: cur_condition

  allocate(condition_array(num_conditions))

  cur_condition => condition_list
  do
    if (.not.associated(cur_condition)) exit
    condition_array(cur_condition%id)%ptr => cur_condition
    cur_condition => cur_condition%next
  enddo

end subroutine ConvertListToArray

! ************************************************************************** !
!
! UpdateConditions: Updates the value of the boundary/initial condition based
!                   on the current time
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine UpdateConditions(time)

  implicit none

  real*8 :: time

  type(condition_type), pointer :: cur_condition

  integer :: cur_time_index, next_time_index
  real*8 :: time_fraction

  cur_condition => condition_list
  do

    if (.not.associated(cur_condition)) exit

    cur_time_index = cur_condition%cur_time_index
    next_time_index = min(cur_condition%cur_time_index+1, &
                          cur_condition%max_time_index)

    ! ensure that condition has started
    if (time >= cur_condition%times(cur_time_index) .or. &
        abs(time-cur_condition%times(cur_time_index)) < 1.d-40) then

      ! find appropriate time interval
      do
        if (time < cur_condition%times(next_time_index) .or. &
            cur_time_index == next_time_index) &
          exit
        cur_time_index = next_time_index
        ! ensure that time index does not go beyond end of array
        if (next_time_index < cur_condition%max_time_index) &
          next_time_index = next_time_index + 1
      enddo

      ! interpolate value based on time
      if (cur_time_index < cur_condition%max_time_index) then
        time_fraction = (time-cur_condition%times(cur_time_index)) / &
                        (cur_condition%times(next_time_index) - &
                         cur_condition%times(cur_time_index))
        cur_condition%cur_value = cur_condition%values(cur_time_index) + &
                                  time_fraction * &
                                  (cur_condition%values(next_time_index) - &
                                   cur_condition%values(cur_time_index))
      else
        cur_condition%cur_value =  &
                             cur_condition%values(cur_condition%max_time_index) 
      endif 
    endif 

    cur_condition => cur_condition%next
    
  enddo

end subroutine UpdateConditions

! ************************************************************************** !
!
! InitializeBoundaryConditions: Initializes boundary condition type
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine InitializeBoundaryConditions(grid)

  use pflow_gridtype_module
  use fileio_module

  implicit none

  type(pflowGrid) :: grid

  integer :: icond
  character(len=MAXWORDLENGTH) :: word
  type(condition_type), pointer :: cur_condition

  call ConvertListToArray()

  do icond = 1, num_conditions
    word = condition_array(icond)%ptr%ctype 
    call fiWordToLower(word)
    if (fiStringCompare("dirichlet",word,9)) then
      condition_array(icond)%ptr%itype = 1  ! integer pointer to condition
    else if (fiStringCompare("neumann",word,7)) then
      condition_array(icond)%ptr%itype = 2
    else if (fiStringCompare("hydraulic gradient",word,18)) then
      condition_array(icond)%ptr%itype = 3
    else if (fiStringCompare("seepage face",word,12)) then
      condition_array(icond)%ptr%itype = 4
    else
      if (grid%myrank == 0) print *, 'Condition type: ', word, ' not supported'
      stop
    endif

  enddo

  ! set up native pflotran integer pointer to boundary condition type
  allocate(grid%ibndtyp(num_conditions)) ! condition (boundary) type
  do icond = 1, num_conditions
    if (condition_array(icond)%ptr%itype == 2) then
      grid%ibndtyp(icond) = 2  ! neumann
    else
      grid%ibndtyp(icond) = 3  ! pressure based dirichlet
    endif
  enddo

  ! updated to initial boundary condition values
  call UpdateBoundaryConditions(grid)

end subroutine InitializeBoundaryConditions

! ************************************************************************** !
!
! UpdateBoundaryConditions: Updates the value of the boundary condition based
!                           on the current time
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine UpdateBoundaryConditions(grid)

  use pflow_gridtype_module

  implicit none

  type(pflowGrid) :: grid

  logical :: use_eos
  integer :: iconnbc, icond, itype
  real*8 :: value, datum, delz, delp
  real*8 :: patm = 101325.d0, tref, pref
  real*8 :: datum_coord(3), cell_coord(3)

  call UpdateConditions(grid%t)

  use_eos = .true.
  
  tref = 25.d0
  
  do iconnbc = 1, grid%nconnbc
    icond = grid%ibconn(iconnbc)
    itype = condition_array(icond)%ptr%itype
    
    datum_coord = 0.d0
    datum_coord(3) = condition_array(icond)%ptr%datum
    
    cell_coord = 0.d0
    cell_coord(3) = grid%z(grid%nL2A(grid%mblkbc(iconnbc))+1)
    
    pref = condition_array(icond)%ptr%cur_value
    value = pref
    if (itype == 3 .or. itype == 4) then ! correct for pressure gradient
#if 1
      value = ComputeHydrostaticPressure(grid,use_eos, &
                                         datum_coord,cell_coord, &
                                         pref,0.d0,0.d0,9.81d0*998.32, &
                                         tref,0.d0,0.d0,0.d0)
#else
      delz = grid%z(grid%nL2A(grid%mblkbc(iconnbc))+1)-datum
      delp = delz*9.81d0*998.32
      value = pref - delp
#endif
    endif
    if (itype == 1) then ! dirichlet
      grid%xxbc(1,iconnbc) = value
      grid%xxbc(2,iconnbc) = tref  ! currently hardwired temperature
      grid%xxbc(3,iconnbc) = 1.d-6  ! currently hardwired solute concentration
    elseif (itype == 2) then ! neumann
      grid%velocitybc(1:grid%nphase,iconnbc) = value  ! all xxbc dofs get over-
     ! grid%xxbc(1:grid%nphase,iconnbc) = patm        ! written later in 
    else                                              ! pflow_vadose_ResJac.F90
      grid%xxbc(1,iconnbc) = value   ! dofs 2+ get overwritten later 
    endif                            ! in pflow_vadose_ResJac.F90
  enddo

end subroutine UpdateBoundaryConditions


! ************************************************************************** !
!
! ComputeInitialCondition: Sets up initial condition for primary variables
! author: Glenn Hammond
! date: 5/8/07
!
! ************************************************************************** !
subroutine ComputeInitialCondition(grid, icondition)

  use pflow_gridtype_module
  use translator_vad_module
  use Vadose_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  type(pflowGrid) :: grid
  integer :: icondition
  
  integer :: iln, na, ierr, icond, itype
  real*8 :: cell_coord(3), datum_coord(3)
  real*8 :: pref, tref, sref, value
  real*8 :: patm = 101325.d0
  PetscScalar, pointer :: xx_p(:), iphase_p(:)
  
  do icond = 1, num_conditions
    if (condition_array(icond)%ptr%id == icondition) then
      exit
    endif
  enddo

  if (icond > num_conditions) then
    if (grid%myrank == 0) &
      print *, 'Condition: ', icondition, ' not found in list of conditions'
    stop
  endif
  
  itype = condition_array(icond)%ptr%itype
  if (itype == 2) then
    if (grid%myrank == 0) &
      print *, 'Neumann bc not supported for initial condition.'
    stop
  endif
  
  pref = condition_array(icond)%ptr%cur_value
  tref = 25.d0
  sref = 1.d-6
  
  datum_coord(1) = 0.d0
  datum_coord(2) = 0.d0
  datum_coord(3) = condition_array(icond)%ptr%datum
  
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)
  
  do iln=1, grid%nlmax
  
    na = grid%nL2A(iln)+1
    cell_coord(1) = grid%x(na)
    cell_coord(2) = grid%y(na)
    cell_coord(3) = grid%z(na)
    
    if (itype == 3 .or. itype == 4) then ! correct for pressure gradient
      value = ComputeHydrostaticPressure(grid,.true., &
                                         datum_coord,cell_coord, &
                                         pref,0.d0,0.d0,0.d0, &
                                         tref,0.d0,0.d0,0.d0)
    else
      value = pref
    endif     

!geh    iphase_p(iln) = grid%iphas_ini(ir)
    if (value > patm) then
      iphase_p(iln) = 1
    else
      iphase_p(iln) = 3
    endif

    xx_p(1+(iln-1)*grid%ndof) = value  
    xx_p(2+(iln-1)*grid%ndof) = tref
    xx_p(3+(iln-1)*grid%ndof) = sref
    
  enddo
              
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)

  call pflow_update_vadose(grid)
  call Translator_Vadose_Switching(grid%xx,grid,0,ierr)

end subroutine ComputeInitialCondition

! ************************************************************************** !
!
! ComputeHydrostaticPressure: Computes the hydrostatic pressure at a given
!                             cell
! author: Glenn Hammond
! date: 5/8/07
!
! ************************************************************************** !
real*8 function ComputeHydrostaticPressure(grid, use_eos, datum_coordinate, &
                                           cell_coordinate, &
                                           reference_pressure, &
                                           pressure_gradient_X, &
                                           pressure_gradient_Y, &
                                           pressure_gradient_Z, &
                                           reference_temperature, &
                                           temperature_gradient_X, &
                                           temperature_gradient_Y, &
                                           temperature_gradient_Z)

  use pflow_gridtype_module
  use water_eos_module

  implicit none

  type(pflowGrid) :: grid
  logical :: use_eos
  real*8 :: datum_coordinate(3), cell_coordinate(3)
  real*8 :: reference_pressure, pressure_gradient_X, pressure_gradient_Y,  &
            pressure_gradient_Z
  real*8 :: reference_temperature, temperature_gradient_X, &
            temperature_gradient_Y, temperature_gradient_Z

  integer :: i, num_increments, ierr
  real*8 :: dx, dy, dz, z, dum, rho, dw_mol, dwp
  real*8 :: temperature_at_xy, temperature_at_xyz
  real*8 :: pressure_at_xy, pressure_at_xyz
  real*8 :: increment, final_increment

  dx = cell_coordinate(1) - datum_coordinate(1)
  dy = cell_coordinate(2) - datum_coordinate(2)
  dz = cell_coordinate(3) - datum_coordinate(3)

  pressure_at_xy = reference_pressure + dx*pressure_gradient_X + &
                                        dy*pressure_gradient_Y
  if (use_eos) then
    temperature_at_xy = reference_temperature + dx*temperature_gradient_X + &
                                                dy*temperature_gradient_Y

    num_increments = int(abs(dz))  ! 1m increments
    final_increment = sign(abs(dz)-num_increments*1.d0,-dz)
    increment = 1.d0
    if (datum_coordinate(3) < cell_coordinate(3)) increment = -1.d0

    z = datum_coordinate(3)
    pressure_at_xyz = pressure_at_xy
    do i=1, num_increments
      z = z - increment
      temperature_at_xyz = temperature_at_xy + (z-datum_coordinate(3))* &
                                                temperature_gradient_Z
      call wateos(temperature_at_xyz,pressure_at_xyz,rho,dw_mol,dwp, &
                  dum,dum,dum,dum,grid%scale,ierr)
      pressure_at_xyz = pressure_at_xyz + increment*rho*grid%gravity  
    enddo
    if (final_increment > 0.d0) then
      z = z - final_increment
      temperature_at_xyz = temperature_at_xy + (z-datum_coordinate(3))* &
                                               temperature_gradient_Z
      call wateos(temperature_at_xyz,pressure_at_xyz,rho,dw_mol,dwp, &
                  dum,dum,dum,dum,grid%scale,ierr)
      pressure_at_xyz = pressure_at_xyz + final_increment*rho*grid%gravity  
    endif
  else
    pressure_at_xyz = pressure_at_xy - dz*pressure_gradient_Z
  endif

  computeHydrostaticPressure = pressure_at_xyz
    
end function ComputeHydrostaticPressure

end module Condition_module
