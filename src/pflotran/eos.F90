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
  PetscReal :: tempreal
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
              case('IFC67')
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
              case('IFC67')
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
      call EOSWaterVerify(ierr)
      if (ierr /= 0) then
        option%io_buffer = 'Error in Water EOS'    
      call printErrMsg(option)
      endif
    case default
      option%io_buffer = 'Keyword: ' // trim(keyword) // &
                         ' not recognized in EOS'    
      call printErrMsg(option)
  end select
  
end subroutine EOSRead

end module EOS_module
