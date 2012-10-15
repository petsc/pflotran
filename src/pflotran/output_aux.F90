module Output_Aux_module

  implicit none

  private

#include "definitions.h"


  type, public :: output_option_type

    character(len=2) :: tunit
    PetscReal :: tconv

    PetscBool :: print_initial
    PetscBool :: print_final
  
    PetscBool :: print_hdf5
    PetscBool :: print_hdf5_velocities
    PetscBool :: print_hdf5_flux_velocities
    PetscBool :: print_single_h5_file

    PetscBool :: print_tecplot 
    PetscInt :: tecplot_format
    PetscBool :: print_tecplot_velocities
    PetscBool :: print_tecplot_flux_velocities
    
    PetscBool :: print_vtk 
    PetscBool :: print_vtk_velocities

    PetscBool :: print_observation 
    PetscBool :: print_column_ids

    PetscBool :: print_mad 

    PetscInt :: screen_imod
    PetscInt :: output_file_imod
    
    PetscInt :: periodic_output_ts_imod
    PetscInt :: periodic_tr_output_ts_imod
    
    PetscReal :: periodic_output_time_incr
    PetscReal :: periodic_tr_output_time_incr
    
    PetscBool :: print_permeability
    PetscBool :: print_porosity
    
    PetscInt :: plot_number
    character(len=MAXWORDLENGTH) :: plot_name

#ifdef SURFACE_FLOW
    PetscBool :: print_hydrograph
#endif

  end type output_option_type
  
  type, public :: output_variable_list_type
    type(output_variable_type), pointer :: first
    type(output_variable_type), pointer :: last
  end type output_variable_list_type
  
  type, public :: output_variable_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: ivar
    PetscInt :: isubvar
    PetscInt :: isubsubvar
    type(output_variable_type), pointer :: next
  end type output_variable_type
    
  public :: OutputOptionCreate, &
            OutputVariableCreate, &
            OutputVariableListCreate, &
            OutputOptionDestroy, &
            OutputVariableListDestroy

contains

! ************************************************************************** !
!
! OutputOptionCreate: Creates output options object
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function OutputOptionCreate()

  implicit none
  
  type(output_option_type), pointer :: OutputOptionCreate

  type(output_option_type), pointer :: output_option
  
  allocate(output_option)
  output_option%print_hdf5 = PETSC_FALSE
  output_option%print_hdf5_velocities = PETSC_FALSE
  output_option%print_hdf5_flux_velocities = PETSC_FALSE
  output_option%print_single_h5_file = PETSC_TRUE
  output_option%print_tecplot = PETSC_FALSE
  output_option%tecplot_format = 0
  output_option%print_tecplot_velocities = PETSC_FALSE
  output_option%print_tecplot_flux_velocities = PETSC_FALSE
  output_option%print_vtk = PETSC_FALSE
  output_option%print_vtk_velocities = PETSC_FALSE
  output_option%print_observation = PETSC_FALSE
  output_option%print_column_ids = PETSC_FALSE
  output_option%print_mad = PETSC_FALSE
  output_option%print_initial = PETSC_TRUE
  output_option%print_final = PETSC_TRUE
  output_option%plot_number = 0
  output_option%screen_imod = 1
  output_option%output_file_imod = 1
  output_option%periodic_output_ts_imod  = 100000000
  output_option%periodic_output_time_incr = 0.d0
  output_option%periodic_tr_output_ts_imod = 100000000
  output_option%periodic_tr_output_time_incr = 0.d0
  output_option%plot_name = ""
  output_option%print_permeability = PETSC_FALSE
  output_option%print_porosity = PETSC_FALSE
  
  output_option%tconv = 1.d0
  output_option%tunit = 's'
  
#ifdef SURFACE_FLOW
  output_option%print_hydrograph = PETSC_FALSE
#endif

  OutputOptionCreate => output_option
  
end function OutputOptionCreate

! ************************************************************************** !
!
! OutputVariableCreate: initializes output variable object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableCreate()

  implicit none
  
  type(output_variable_type), pointer :: OutputVariableCreate
  
  type(output_variable_type), pointer :: output_variable
  
  allocate(output_variable)
  output_variable%name = ''
  output_variable%ivar = 0
  output_variable%isubvar = 0
  output_variable%isubsubvar = 0
  nullify(output_variable%next)
  
  OutputVariableCreate => output_variable
  
end function OutputVariableCreate

! ************************************************************************** !
!
! OutputVariableListCreate: initializes output variable list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableListCreate()

  implicit none
  
  type(output_variable_list_type), pointer :: OutputVariableListCreate
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  allocate(output_variable_list)
  nullify(output_variable_list%first)
  nullify(output_variable_list%last)
  
  OutputVariableListCreate => output_variable_list
  
end function OutputVariableListCreate

! ************************************************************************** !
!
! OutputVariableAddToList: initializes output variable list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine OutputVariableAddToList(list,variable)

  implicit none
  
  type(output_variable_list_type), pointer :: list
  type(output_variable_type), pointer :: variable
  
  list%last%next => variable
  list%last => variable
  
end subroutine OutputVariableAddToList

! ************************************************************************** !
!
! OutputVariableListDestroy: Deallocates an output variable list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine OutputVariableListDestroy(output_variable_list)

  implicit none
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  nullify(output_variable_list%last)
  call OutputVariableDestroy(output_variable_list%first)
  
  deallocate(output_variable_list)
  nullify(output_variable_list)
  
end subroutine OutputVariableListDestroy

! ************************************************************************** !
!
! OutputVariableDestroy: Deallocates an output variable object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
recursive subroutine OutputVariableDestroy(output_variable)

  implicit none
  
  type(output_variable_type), pointer :: output_variable
  
  if (.not.associated(output_variable)) return
  
  call OutputVariableDestroy(output_variable%next)
  
  deallocate(output_variable)
  nullify(output_variable)
  
end subroutine OutputVariableDestroy

! ************************************************************************** !
!
! OutputOptionDestroy: Deallocates an output option
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
subroutine OutputOptionDestroy(output_option)

  implicit none
  
  type(output_option_type), pointer :: output_option
  

  deallocate(output_option)
  nullify(output_option)
  
end subroutine OutputOptionDestroy

end module Output_Aux_module
