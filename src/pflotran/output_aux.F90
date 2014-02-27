module Output_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  PetscInt, parameter, public :: INSTANTANEOUS_VARS = 1
  PetscInt, parameter, public :: AVERAGED_VARS = 2

  type, public :: output_option_type

    character(len=2) :: tunit
    PetscReal :: tconv

    PetscBool :: print_initial
    PetscBool :: print_final
  
    PetscBool :: print_hdf5
    PetscBool :: print_hdf5_velocities
    PetscBool :: print_hdf5_flux_velocities
    PetscBool :: print_single_h5_file
    PetscInt  :: times_per_h5_file
    PetscBool :: print_hdf5_mass_flowrate
    PetscBool :: print_hdf5_energy_flowrate
    PetscBool :: print_hdf5_aveg_mass_flowrate
    PetscBool :: print_hdf5_aveg_energy_flowrate
    PetscBool :: print_explicit_flowrate

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
    PetscReal :: periodic_checkpoint_time_incr

    PetscBool :: print_permeability
    PetscBool :: print_porosity
    PetscBool :: print_iproc
    PetscBool :: print_volume
    PetscBool :: print_tortuosity

    PetscInt :: xmf_vert_len
    
    type(output_variable_list_type), pointer :: output_variable_list
    type(output_variable_list_type), pointer :: aveg_output_variable_list

    PetscReal :: aveg_var_time
    PetscReal :: aveg_var_dtime
    
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
    PetscInt :: nvars
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
    type(output_variable_type), pointer :: next
  end type output_variable_type

!  type, public, EXTENDS (output_variable_type) :: aveg_output_variable_type
!    PetscReal :: time_interval
!  end type aveg_output_variable_type
  
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

function OutputOptionCreate()
  ! 
  ! Creates output options object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(output_option_type), pointer :: OutputOptionCreate

  type(output_option_type), pointer :: output_option
  
  allocate(output_option)
  output_option%print_hdf5 = PETSC_FALSE
  output_option%print_hdf5_velocities = PETSC_FALSE
  output_option%print_hdf5_flux_velocities = PETSC_FALSE
  output_option%print_single_h5_file = PETSC_TRUE
  output_option%times_per_h5_file = 0
  output_option%print_hdf5_mass_flowrate = PETSC_FALSE
  output_option%print_hdf5_energy_flowrate = PETSC_FALSE
  output_option%print_hdf5_aveg_mass_flowrate = PETSC_FALSE
  output_option%print_hdf5_aveg_energy_flowrate = PETSC_FALSE
  output_option%print_explicit_flowrate = PETSC_FALSE
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
  output_option%print_iproc = PETSC_FALSE
  output_option%print_volume = PETSC_FALSE
  output_option%print_tortuosity = PETSC_FALSE
  output_option%aveg_var_time = 0.d0
  output_option%aveg_var_dtime = 0.d0
  output_option%periodic_checkpoint_time_incr = 0.d0
  
  nullify(output_option%output_variable_list)
  nullify(output_option%aveg_output_variable_list)
  
  output_option%tconv = 1.d0
  output_option%tunit = 's'
  
#ifdef SURFACE_FLOW
  output_option%print_hydrograph = PETSC_FALSE
#endif

  OutputOptionCreate => output_option
  
end function OutputOptionCreate

! ************************************************************************** !

function OutputVariableCreate1()
  ! 
  ! initializes output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

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
  nullify(output_variable%next)
  
  OutputVariableCreate1 => output_variable
  
end function OutputVariableCreate1

! ************************************************************************** !

function OutputVariableCreate2(name,icategory,units,ivar,isubvar,isubsubvar)
  ! 
  ! initializes output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

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

function OutputVariableCreate3(output_variable)
  ! 
  ! initializes output variable object from an existing
  ! output variabl object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

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
  nullify(new_output_variable%next)
  
  OutputVariableCreate3 => new_output_variable
  
end function OutputVariableCreate3

! ************************************************************************** !

function OutputVariableListCreate()
  ! 
  ! initializes output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type), pointer :: OutputVariableListCreate
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  allocate(output_variable_list)
  nullify(output_variable_list%first)
  nullify(output_variable_list%last)
  output_variable_list%nvars = 0
  
  OutputVariableListCreate => output_variable_list
  
end function OutputVariableListCreate

! ************************************************************************** !

function OutputVariableListDuplicate(old_list,new_list)
  ! 
  ! initializes output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type) :: old_list
  
  type(output_variable_list_type), pointer :: OutputVariableListDuplicate
  
  type(output_variable_list_type), pointer :: new_list
  type(output_variable_type), pointer :: cur_variable
  
  allocate(new_list)
  nullify(new_list%first)
  nullify(new_list%last)
  new_list%nvars = old_list%nvars
  
  cur_variable => old_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputVariableAddToList(new_list,OutputVariableCreate(cur_variable))
    cur_variable => cur_variable%next
  enddo

  OutputVariableListDuplicate => new_list
  
end function OutputVariableListDuplicate

! ************************************************************************** !

subroutine OutputVariableAddToList1(list,variable)
  ! 
  ! adds variable to list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type) :: list
  type(output_variable_type), pointer :: variable
  
  if (.not. associated(list%first)) then
    list%first => variable
  else
    list%last%next => variable
  endif
  list%last => variable
  
  list%nvars = list%nvars+1
  
end subroutine OutputVariableAddToList1

! ************************************************************************** !

subroutine OutputVariableAddToList2(list,name,icategory,units,ivar, &
                                    isubvar,isubsubvar)
  ! 
  ! creates variable and adds to list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

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

function OutputVariableListToHeader(variable_list,cell_string,icolumn, &
                                    plot_file)
  ! 
  ! Converts a variable list to a header string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
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

subroutine OutputAppendToHeader(header,variable_string,units_string, &
                                cell_string, icolumn)
  ! 
  ! Appends formatted strings to header string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/27/11
  ! 

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

function OutputVariableToCategoryString(icategory)
  ! 
  ! returns a string associated with an
  ! output variable category
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

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

subroutine OutputVariableRead(input,option,output_variable_list)
  ! 
  ! This routine reads variable from input file.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/21/12
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Variables_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(output_variable_list_type), pointer :: output_variable_list
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable
  PetscInt :: temp_int

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','VARIABLES')
    call StringToUpper(word)
    
    select case(trim(word))
      case ('LIQUID_PRESSURE')
        name = 'Liquid Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     LIQUID_PRESSURE)
      case ('LIQUID_SATURATION')
        name = 'Liquid Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     LIQUID_SATURATION)
      case ('LIQUID_DENSITY')
        name = 'Liquid Density'
        units = 'kg/m^3'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_DENSITY)
      case ('LIQUID_MOBILITY')
        name = 'Liquid Mobility'
        units = '1/Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOBILITY)
      case ('LIQUID_ENERGY')
        name = 'Liquid Energy'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'PER_VOLUME')) then
            units = 'MJ/m^3'
            temp_int = ONE_INTEGER
          else
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,LIQUID_ENERGY')
          endif
        else
          units = 'MJ/kmol'
          temp_int = ZERO_INTEGER
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_ENERGY,temp_int)
    
      case ('GAS_PRESSURE')
        name = 'Gas Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     GAS_PRESSURE)
      case ('GAS_SATURATION')
        name = 'Gas Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     GAS_SATURATION)
      case ('GAS_DENSITY')
        name = 'Gas Density'
        units = 'kg/m^3'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_DENSITY)
      case ('GAS_MOBILITY')
        name = 'Gas Mobility'
        units = '1/Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOBILITY)
      case ('GAS_ENERGY')
        name = 'Gas Energy'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'PER_VOLUME')) then
            units = 'MJ/m^3'
            temp_int = ONE_INTEGER
          else
            input%ierr = 1
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,GAS_ENERGY')
          endif
        else
          units = 'MJ/kmol'
          temp_int = ZERO_INTEGER
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_ENERGY,temp_int)
      case ('LIQUID_MOLE_FRACTIONS')
        name = 'X_g^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOLE_FRACTION, &
                                     option%air_id)
        name = 'X_l^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOLE_FRACTION, &
                                     option%water_id)
      case ('GAS_MOLE_FRACTIONS')
        name = 'X_g^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOLE_FRACTION, &
                                     option%air_id)
        name = 'X_l^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOLE_FRACTION, &
                                     option%water_id)
      case ('AIR_PRESSURE')
        name = 'Air Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     AIR_PRESSURE)
      case ('CAPILLARY_PRESSURE')
        name = 'Capillary Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     CAPILLARY_PRESSURE)
      case('THERMODYNAMIC_STATE')
        name = 'Thermodynamic State'
         units = ''
         output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                 units,STATE)
         ! toggle output off for observation
         output_variable%plot_only = PETSC_TRUE 
         output_variable%iformat = 1 ! integer
         call OutputVariableAddToList(output_variable_list,output_variable)
         nullify(output_variable)
      case ('TEMPERATURE')
        name = 'Temperature'
        units = 'C'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     TEMPERATURE)
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in VARIABLES.'
        call printErrMsg(option)
    end select

  enddo

end subroutine OutputVariableRead

! ************************************************************************** !

subroutine OutputVariableListDestroy(output_variable_list)
  ! 
  ! Deallocates an output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  nullify(output_variable_list%last)
  call OutputVariableDestroy(output_variable_list%first)
  
  deallocate(output_variable_list)
  nullify(output_variable_list)
  
end subroutine OutputVariableListDestroy

! ************************************************************************** !

recursive subroutine OutputVariableDestroy(output_variable)
  ! 
  ! Deallocates an output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_type), pointer :: output_variable
  
  if (.not.associated(output_variable)) return
  
  call OutputVariableDestroy(output_variable%next)
  
  deallocate(output_variable)
  nullify(output_variable)
  
end subroutine OutputVariableDestroy

! ************************************************************************** !

subroutine OutputOptionDestroy(output_option)
  ! 
  ! Deallocates an output option
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(output_option_type), pointer :: output_option
  
  if (.not.associated(output_option)) return
  
  call OutputVariableListDestroy(output_option%output_variable_list)
  call OutputVariableListDestroy(output_option%aveg_output_variable_list)
  
  deallocate(output_option)
  nullify(output_option)
  
end subroutine OutputOptionDestroy

end module Output_Aux_module
