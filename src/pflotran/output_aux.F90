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

    PetscInt :: xmf_vert_len
    
    type(output_variable_list_type), pointer :: output_variable_list
    
    PetscInt :: plot_number
    character(len=MAXWORDLENGTH) :: plot_name

#ifdef SURFACE_FLOW
    PetscBool :: print_hydrograph
    PetscInt  :: surf_xmf_vert_len
#endif

  end type output_option_type
  
  type, public :: output_variable_list_type
    type(output_variable_type), pointer :: first
    type(output_variable_type), pointer :: last
  end type output_variable_list_type
  
  type, public :: output_variable_type
    character(len=MAXWORDLENGTH) :: name   ! string that appears in hdf5 file
    character(len=MAXWORDLENGTH) :: units
    PetscBool :: plot_only
    PetscInt :: iformat   ! 0 = for REAL values; 1 = for INTEGER values
    PetscInt :: icategory ! category for variable-specific regression testing
    PetscInt :: ivar
    PetscInt :: isubvar
    PetscInt :: isubsubvar
    PetscBool :: itaveg ! PETSC_FALSE for instantaneous values; 
                        ! PETSC_TRUE for temporally averaged values
    type(output_variable_type), pointer :: next
  end type output_variable_type
  
  interface OutputVariableCreate
    module procedure OutputVariableCreate1
    module procedure OutputVariableCreate2
    module procedure OutputVariableCreate3
  end interface OutputVariableCreate
  
  interface OutputVariableAddToList
    module procedure OutputVariableAddToList1
    module procedure OutputVariableAddToList2
  end interface OutputVariableAddToList
  
  ! Output categories
  PetscInt, parameter, public :: OUTPUT_GENERIC = 0
  PetscInt, parameter, public :: OUTPUT_PRESSURE = 1
  PetscInt, parameter, public :: OUTPUT_SATURATION = 2
  PetscInt, parameter, public :: OUTPUT_CONCENTRATION = 3
  PetscInt, parameter, public :: OUTPUT_RATE = 4
  PetscInt, parameter, public :: OUTPUT_VOLUME_FRACTION = 5
  PetscInt, parameter, public :: OUTPUT_DISCRETE = 6
  
  public :: OutputOptionCreate, &
            OutputVariableCreate, &
            OutputVariableListCreate, &
            OutputVariableListDuplicate, &
            OutputVariableAddToList, &
            OutputAppendToHeader, &
            OutputVariableListToHeader, &
            OutputVariableToCategoryString, &
            OutputVariableRead, &
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
  
  nullify(output_option%output_variable_list)
  
  output_option%tconv = 1.d0
  output_option%tunit = 's'
  
#ifdef SURFACE_FLOW
  output_option%print_hydrograph = PETSC_FALSE
#endif

  OutputOptionCreate => output_option
  
end function OutputOptionCreate

! ************************************************************************** !
!
! OutputVariableCreate1: initializes output variable object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableCreate1()

  implicit none
  
  type(output_variable_type), pointer :: OutputVariableCreate1
  
  type(output_variable_type), pointer :: output_variable
  
  allocate(output_variable)
  output_variable%name = ''
  output_variable%units = ''
  output_variable%plot_only = PETSC_FALSE
  output_variable%iformat = 0
  output_variable%icategory = OUTPUT_GENERIC
  output_variable%ivar = 0
  output_variable%isubvar = 0
  output_variable%isubsubvar = 0
  output_variable%itaveg = PETSC_FALSE
  nullify(output_variable%next)
  
  OutputVariableCreate1 => output_variable
  
end function OutputVariableCreate1

! ************************************************************************** !
!
! OutputVariableCreate2: initializes output variable object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableCreate2(name,icategory,units,ivar,isubvar,isubsubvar)

  implicit none
  
  character(len=*) :: name
  PetscInt :: icategory ! note that I tuck it inbetween the strings to avoid
                        ! errors
  character(len=*) :: units
  PetscInt :: ivar
  PetscInt, intent(in), optional :: isubvar
  PetscInt, intent(in), optional :: isubsubvar

  type(output_variable_type), pointer :: OutputVariableCreate2
  
  type(output_variable_type), pointer :: output_variable
  
  output_variable => OutputVariableCreate()
  output_variable%name = trim(adjustl(name))
  output_variable%icategory = icategory
  output_variable%units = trim(adjustl(units))
  output_variable%ivar = ivar
  if (present(isubvar)) then
    output_variable%isubvar = isubvar
  endif
  if (present(isubsubvar)) then
    output_variable%isubsubvar = isubsubvar
  endif
  nullify(output_variable%next)
  
  OutputVariableCreate2 => output_variable
  
end function OutputVariableCreate2

! ************************************************************************** !
!
! OutputVariableCreate3: initializes output variable object from an existing
!                        output variabl object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableCreate3(output_variable)

  implicit none
  
  type(output_variable_type), pointer :: output_variable

  type(output_variable_type), pointer :: OutputVariableCreate3
  
  type(output_variable_type), pointer :: new_output_variable
  
  allocate(new_output_variable)
  new_output_variable%name = output_variable%name
  new_output_variable%units = output_variable%units
  new_output_variable%plot_only = output_variable%plot_only
  new_output_variable%iformat = output_variable%iformat
  new_output_variable%icategory = output_variable%icategory
  new_output_variable%ivar = output_variable%ivar
  new_output_variable%isubvar = output_variable%isubvar
  new_output_variable%isubsubvar = output_variable%isubsubvar
  new_output_variable%itaveg = output_variable%itaveg
  nullify(new_output_variable%next)
  
  OutputVariableCreate3 => new_output_variable
  
end function OutputVariableCreate3

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
! OutputVariableListDuplicate: initializes output variable list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableListDuplicate(old_list,new_list)

  implicit none
  
  type(output_variable_list_type) :: old_list
  
  type(output_variable_list_type), pointer :: OutputVariableListDuplicate
  
  type(output_variable_list_type), pointer :: new_list
  type(output_variable_type), pointer :: cur_variable
  
  allocate(new_list)
  nullify(new_list%first)
  nullify(new_list%last)
  
  cur_variable => old_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputVariableAddToList(new_list,OutputVariableCreate(cur_variable))
    cur_variable => cur_variable%next
  enddo

  OutputVariableListDuplicate => new_list
  
end function OutputVariableListDuplicate

! ************************************************************************** !
!
! OutputVariableAddToList1: adds variable to list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine OutputVariableAddToList1(list,variable)

  implicit none
  
  type(output_variable_list_type) :: list
  type(output_variable_type), pointer :: variable
  
  if (.not. associated(list%first)) then
    list%first => variable
  else
    list%last%next => variable
  endif
  list%last => variable
  
end subroutine OutputVariableAddToList1

! ************************************************************************** !
!
! OutputVariableAddToList2: creates variable and adds to list object
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine OutputVariableAddToList2(list,name,icategory,units,ivar, &
                                    isubvar,isubsubvar)

  implicit none
  
  type(output_variable_list_type) :: list
  character(len=*) :: name
  character(len=*) :: units
  PetscInt :: icategory
  PetscInt :: ivar
  PetscInt, intent(in), optional :: isubvar
  PetscInt, intent(in), optional :: isubsubvar
  
  type(output_variable_type), pointer :: variable
  
  if (present(isubvar)) then
    if (present(isubsubvar)) then
      variable => OutputVariableCreate(name,icategory,units, &
                                       ivar,isubvar,isubsubvar)
    else
      variable => OutputVariableCreate(name,icategory,units, &
                                       ivar,isubvar)
    endif
  else
    variable => OutputVariableCreate(name,icategory,units,ivar)
  endif
  call OutputVariableAddToList1(list,variable)
  
end subroutine OutputVariableAddToList2

! ************************************************************************** !
!
! OutputVariableListToHeader: Converts a variable list to a header string
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableListToHeader(variable_list,cell_string,icolumn, &
                                    plot_file)
  
  implicit none
  
  type(output_variable_list_type) :: variable_list
  character(len=*) :: cell_string
  PetscInt :: icolumn
  PetscBool :: plot_file
  
  character(len=MAXHEADERLENGTH) :: OutputVariableListToHeader

  character(len=MAXHEADERLENGTH) :: header
  type(output_variable_type), pointer :: cur_variable
  character(len=MAXWORDLENGTH) :: variable_name, units
  
  header = ''
  
  cur_variable => variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    if (.not. plot_file .and. cur_variable%plot_only) then
      cur_variable => cur_variable%next
      cycle
    endif
    variable_name = cur_variable%name
    units = cur_variable%units
    call OutputAppendToHeader(header,variable_name,units,cell_string,icolumn)
    cur_variable => cur_variable%next
  enddo
  
  OutputVariableListToHeader = header
  
end function OutputVariableListToHeader

! ************************************************************************** !
!
! OutputAppendToHeader: Appends formatted strings to header string
! author: Glenn Hammond
! date: 10/27/11
!
! ************************************************************************** !
subroutine OutputAppendToHeader(header,variable_string,units_string, &
                                cell_string, icolumn)

  implicit none

  character(len=MAXHEADERLENGTH) :: header
  character(len=*) :: variable_string, units_string, cell_string
  character(len=MAXWORDLENGTH) :: column_string
  character(len=MAXWORDLENGTH) :: variable_string_adj, units_string_adj
  character(len=MAXSTRINGLENGTH) :: cell_string_adj
  PetscInt :: icolumn, len_cell_string, len_units

  character(len=MAXSTRINGLENGTH) :: string

  variable_string_adj = variable_string
  units_string_adj = units_string
  cell_string_adj = cell_string

  !geh: Shift to left.  Cannot perform on same string since len=*
  variable_string_adj = adjustl(variable_string_adj)
  units_string_adj = adjustl(units_string_adj)
  cell_string_adj = adjustl(cell_string_adj)

  if (icolumn > 0) then
    icolumn = icolumn + 1
    write(column_string,'(i4,''-'')') icolumn
    column_string = trim(adjustl(column_string))
  else
    column_string = ''
  endif

  !geh: this is all to remove the lousy spaces
  len_units = len_trim(units_string)
  len_cell_string = len_trim(cell_string)
  if (len_units > 0 .and. len_cell_string > 0) then
    write(string,'('',"'',a,a,'' ['',a,''] '',a,''"'')') trim(column_string), &
          trim(variable_string_adj), trim(units_string_adj), &
          trim(cell_string_adj)
  else if (len_units > 0 .or. len_cell_string > 0) then
    if (len_units > 0) then
      write(string,'('',"'',a,a,'' ['',a,'']"'')') trim(column_string), &
            trim(variable_string_adj), trim(units_string_adj)
    else
      write(string,'('',"'',a,a,'' '',a,''"'')') trim(column_string), &
            trim(variable_string_adj), trim(cell_string_adj)
    endif
  else
    write(string,'('',"'',a,a,''"'')') trim(column_string), &
          trim(variable_string_adj)
  endif
  header = trim(header) // trim(string)

end subroutine OutputAppendToHeader

! ************************************************************************** !
!
! OutputVariableToCategoryString: returns a string associated with an 
!                                 output variable category
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
function OutputVariableToCategoryString(icategory)

  implicit none
  
  PetscInt :: icategory
  
  character(len=MAXWORDLENGTH) :: OutputVariableToCategoryString
  
  character(len=MAXWORDLENGTH) :: string
  
  select case(icategory)
    case(OUTPUT_GENERIC)
      string = 'GENERIC'
    case(OUTPUT_PRESSURE)
      string = 'PRESSURE'
    case(OUTPUT_SATURATION)
      string = 'SATURATION'
    case(OUTPUT_CONCENTRATION)
      string = 'CONCENTRATION'
    case(OUTPUT_RATE)
      string = 'RATE'
    case(OUTPUT_VOLUME_FRACTION)
      string = 'VOLUME_FRACTION'
    case(OUTPUT_DISCRETE)
      string = 'DISCRETE'
    case default
      string = 'GENERIC'
  end select

  OutputVariableToCategoryString = string

end function OutputVariableToCategoryString

! ************************************************************************** !
!> This routine reads variable from input file.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/21/12
! ************************************************************************** !
subroutine OutputVariableRead(input,option,output_variable_list)

  use Option_module
  use Input_module
  use String_module
  use Variables_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(output_variable_list_type), pointer :: output_variable_list
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name, units

  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','VARIABLES')
    call StringToUpper(word)
    
    select case(trim(word))
      case ('LIQUID_PRESSURE')
        name = 'Liquid Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name,OUTPUT_PRESSURE,units, &
                                     LIQUID_PRESSURE)

      case ('LIQUID_SATURATION')
        name = 'Liquid Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)

!      case ('LIQUID_VELOCITIY_AT_CELL_CENTER')
!        name = 'Liquid Velocity at Cell Center'
!        units = 'm/s'
!        call OutputVariableAddToList(output_variable_list,name,OUTPUT_SATURATION,units, &
!                               LIQUID_VELOCITY_CELL_CENT)
!
!      case ('LIQUID_VELOCITIY_AT_CELL_FACE')
!        name = 'Liquid Velocity at Cell Face'
!        units = 'm/s'
!        call OutputVariableAddToList(output_variable_list,name,OUTPUT_SATURATION,units, &
!                               LIQUID_VELOCITY_CELL_FACE)
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in VARIABLES.'
    end select

    call InputReadWord(input,option,word,PETSC_TRUE)
    if (len_trim(word) > 1) then
      call StringToUpper(word)
      select case(trim(word))
        case('AVEG')
          output_variable_list%last%itaveg=PETSC_TRUE
        case('INST')
          output_variable_list%last%itaveg=PETSC_FALSE
        case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in VARIABLES.'
      end select
    endif
  enddo

end subroutine OutputVariableRead

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
  
  if (.not.associated(output_option)) return
  
  deallocate(output_option)
  nullify(output_option)
  
end subroutine OutputOptionDestroy

end module Output_Aux_module
