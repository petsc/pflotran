module EOS_module
 
  use PFLOTRAN_Constants_module
  use EOS_Water_module
  
  implicit none

  private

#include "finclude/petscsys.h"
 
  public :: EOSInit, &
            EOSRead

contains

! ************************************************************************** !

subroutine EOSInit()

  implicit none
  
  call EOSWaterInit()
  
end subroutine EOSInit

! ************************************************************************** !

subroutine EOSRead(input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal, tempreal2, tempreal3
  PetscErrorCode :: ierr

  input%ierr = 0

  call InputReadWord(input,option,keyword,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','EOS')
  call StringToUpper(keyword)   
      
  select case(trim(keyword))
    case('WATER')
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit  
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,WATER')
        call StringToUpper(keyword)   
        select case(trim(keyword))
          case('DENSITY') 
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,WATER')
            call StringToUpper(word)   
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,DENSITY,CONSTANT')
                call EOSWaterSetDensityConstant(tempreal)
              case('EXPONENTIAL')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,tempreal2)
                call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,tempreal3)
                call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call EOSWaterSetDensityExponential(tempreal,tempreal2, &
                                                   tempreal3)
              case('IFC67','DEFAULT')
                call EOSWaterSetDensityIFC67()
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                                   ' not recognized in EOS,WATER,DENSITY'    
                call printErrMsg(option)
            end select
          case('ENTHALPY') 
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,WATER')
            call StringToUpper(word)   
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,ENTHALPY,CONSTANT')
                call EOSWaterSetEnthalpyConstant(tempreal)
              case('IFC67','DEFAULT')
                call EOSWaterSetEnthalpyIFC67()
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                                   ' not recognized in EOS,WATER,ENTHALPY'    
                call printErrMsg(option)
            end select
          case('VISCOSITY') 
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,WATER')
            call StringToUpper(word)   
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,VISCOSITY,CONSTANT')
                call EOSWaterSetViscosityConstant(tempreal)
              case('DEFAULT')
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                                   ' not recognized in EOS,WATER,VISCOSITY'    
                call printErrMsg(option)
            end select
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
                               ' not recognized in EOS,WATER'    
            call printErrMsg(option)
        end select
      enddo
      string = ''
      call EOSWaterVerify(ierr,string)
      if (ierr /= 0) then
        option%io_buffer = 'Error in Water EOS'    
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif
    case default
      option%io_buffer = 'Keyword: ' // trim(keyword) // &
                         ' not recognized in EOS'    
      call printErrMsg(option)
  end select
  
end subroutine EOSRead

end module EOS_module
