module SrcSink_Sandbox_Downreg_class

! Source/sink Sandbox for downregulating source or sink terms to avoid 
! overpressurization (too big a pressure, + or -)
! the source is regulated by multiplying it by 
!   pressure_max/(pressure + pressure_max)
! the sink is regulated by multiplying it by 
!   pressure_min/(pressure + pressure_min)
! extended from srcsink_sandbox_mass_rate
! currently work in RICHARDS mode
  
  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_downreg_type
    PetscReal, pointer :: rate(:)
    PetscReal :: pressure_min       ! Pa
    PetscReal :: pressure_max       ! Pa
  contains
    procedure, public :: ReadInput => DownregRead
    procedure, public :: Setup => DownregSetup
    procedure, public :: Evaluate => DownregSrcSink
    procedure, public :: Destroy => DownregDestroy
  end type srcsink_sandbox_downreg_type

  public :: DownregCreate

contains

! ************************************************************************** !

function DownregCreate()
  ! 
  ! Allocates mass rate src/sink object.
  ! 
  ! Author: Guoping Tang
  ! Date: 06/03/14

  implicit none
  
  class(srcsink_sandbox_downreg_type), pointer :: DownregCreate

  allocate(DownregCreate)
  call SSSandboxBaseInit(DownregCreate)
  nullify(DownregCreate%rate)
  DownregCreate%pressure_min = -1.0d9
  DownregCreate%pressure_max = 1.0d9
  
end function DownregCreate

! ************************************************************************** !

subroutine DownregRead(this,input,option)
  ! 
  ! Reads input deck for mass rate src/sink parameters
  ! 
  ! Author: Guoping Tang
  ! Date: 06/03/14

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(srcsink_sandbox_downreg_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'SOURCE_SINK_SANDBOX,DOWNREG')
    call StringToUpper(word)   

    ! reads the REGION
    call SSSandboxBaseRead(this,input,option,word,found)
    if (found) cycle
    
    select case(trim(word))

      case('RATE')
        allocate(this%rate(option%nflowdof))
        do i = 1, option%nflowdof
          word = ''
          select case(i)
            case(ONE_INTEGER)
              word = 'Liquid Component Rate'
            case(TWO_INTEGER)
              word = 'Gas Component Rate'
            case(THREE_INTEGER)
              word = 'Energy Rate'
            case default
              write(word,*) i
              option%io_buffer = 'Unknown dof #' // trim(adjustl(word)) // &
                                 ' in DownregRead.'
              call printErrMsg(option)
          end select
          call InputReadDouble(input,option,this%rate(i))
          call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,DOWNREG')
        enddo
      case('POSITIVE_REG_PRESSURE')
        call InputReadDouble(input,option,this%pressure_max)
        call InputErrorMsg(input,option,'maximum pressure (Pa)', &
          'SOURCE_SINK_SANDBOX,DOWNREG')

      case('NEGATIVE_REG_PRESSURE')
        call InputReadDouble(input,option,this%pressure_min)
        call InputErrorMsg(input,option,'minimum pressure (Pa)', &
          'SOURCE_SINK_SANDBOX,DOWNREG')
      case default
        option%io_buffer = 'SRCSINK_SANDBOX,DOWNREG keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine DownregRead

! ************************************************************************** !

subroutine DownregSetup(this,region_list,option)
  ! 
  ! Sets up the mass rate src/sink
  ! 
  ! Author: Guoping Tang
  ! Date: 06/03/14

  use Option_module
  use Region_module
  use General_Aux_module, only : general_fmw_com => fmw_comp

  implicit none
  
  class(srcsink_sandbox_downreg_type) :: this
  type(region_list_type) :: region_list
  type(option_type) :: option
  
  PetscInt :: i
  
  call SSSandboxBaseSetup(this,region_list,option)
  ! convert rate from kg/s to mol/s
  select case(option%iflowmode)
    case(RICHARDS_MODE)
      this%rate(1) = this%rate(1) / FMWH2O
    case(G_MODE)
      this%rate(1) = this%rate(1) / general_fmw_com(1)
      this%rate(2) = this%rate(2) / general_fmw_com(2)
    case default
      option%io_buffer = 'Rate conversion not set up for flow mode in ' // &
                         'DownregSetup'
      call printErrMsg(option)
  end select

end subroutine DownregSetup 

! ************************************************************************** !

subroutine DownregSrcSink(this,Residual,Jacobian,compute_derivative, &
                           material_auxvar,aux_real,option)
  ! 
  ! Evaluates src/sink storing residual and/or Jacobian
  ! 
  ! Author: Guoping Tang
  ! Date: 06/04/14
  ! Glenn suggested to use reference pressure for RICHARDS mode 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_downreg_type) :: this  
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
  PetscReal :: pressure
  PetscReal :: rate_regulator  
  PetscReal :: drate_regulator  
  PetscReal :: temp_real

  PetscInt :: idof
  
  do idof = 1, option%nflowdof
    if (option%iflowmode == RICHARDS_MODE .and. idof == ONE_INTEGER) then
      ! regulate liquid pressure in Richards mode
      pressure = aux_real(1)
      if (this%rate(idof) > 0.0d0) then
        ! source (try to avoid this%pressure_max + pressure = 0)
        ! regulate the source under wet (not dry) condition 
        ! not the case with too much infiltration to dry soil
        if (pressure > option%reference_pressure) then
          rate_regulator = this%pressure_max / &
            (this%pressure_max + pressure - option%reference_pressure)
        else
          rate_regulator = 1.0d0 
        endif
          Residual(idof) = this%rate(idof) * rate_regulator
      else
        ! sink (try to avoid this%pressure_min + pressure = 0)
        ! regulate the sink under dry (not wet) condition, 
        if (pressure < option%reference_pressure) then
          rate_regulator = this%pressure_min / &
            (this%pressure_min + pressure - option%reference_pressure)
        else
          rate_regulator = 1.0d0
        endif
        Residual(idof) = this%rate(idof) * rate_regulator
      endif
    else
        Residual(idof) = this%rate(idof)
    endif
  enddo
  
  if (compute_derivative) then
    
    do idof = 1, option%nflowdof
      if (option%iflowmode == RICHARDS_MODE .and. idof == ONE_INTEGER) then
        ! regulate liquid pressure in Richards mode
        pressure = aux_real(1)
        if (this%rate(idof) > 0.0d0) then
          ! source
          if (pressure > option%reference_pressure) then
            temp_real = this%pressure_max - option%reference_pressure 
            drate_regulator = (-1.0d0) * this%pressure_max / &
              (temp_real + pressure) / (temp_real + pressure)
          else
            drate_regulator = 0.0d0
          endif

          Jacobian(idof,idof) = -1.0d0 * this%rate(idof) * drate_regulator
        else
          ! sink           
          if (pressure < option%reference_pressure) then
            temp_real = this%pressure_min - option%reference_pressure 
            drate_regulator = -1.0d0 * this%pressure_min / &
              (temp_real + pressure) / (temp_real + pressure)
          else
            drate_regulator = 0.0d0
          endif
          Jacobian(idof,idof) = -1.0d0 * this%rate(idof) * drate_regulator
        endif
      else
        ! since the rates are constant, there is no derivative
      endif
    enddo
  endif
  
end subroutine DownregSrcSink

! ************************************************************************** !

subroutine DownregDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this module
  ! 
  ! Author: Guoping Tang
  ! Date: 06/04/14

  implicit none
  
  class(srcsink_sandbox_downreg_type) :: this
  
  call SSSandboxBaseDestroy(this)  
  deallocate(this%rate)
  nullify(this%rate)

end subroutine DownregDestroy

end module SrcSink_Sandbox_Downreg_class
