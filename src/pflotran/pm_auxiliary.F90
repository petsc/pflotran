module PM_Auxiliary_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(pm_base_type) :: pm_auxiliary_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    character(len=MAXWORDLENGTH) :: ctype
    procedure(PMAuxliaryEvaluate), pointer :: Evaluate => null()
  contains
    procedure, public :: InitializeRun => PMAuxiliaryInitializeRun
    procedure, public :: Destroy => PMAuxiliaryDestroy
  end type pm_auxiliary_type
  
  ! interface blocks
  interface
    subroutine PMAuxliaryEvaluate(this,time,ierr)
      import :: pm_auxiliary_type
      implicit none
      class(pm_auxiliary_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr      
    end subroutine PMAuxliaryEvaluate
  end interface
  
  public :: PMAuxiliaryCreate, &
            PMAuxiliaryInit, &
            PMAuxiliaryCast, &
            PMAuxiliarySetFunctionPointer
  
contains

! ************************************************************************** !

function PMAuxiliaryCreate()
  ! 
  ! Creates reactive transport process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_auxiliary_type), pointer :: PMAuxiliaryCreate

  class(pm_auxiliary_type), pointer :: pm

  allocate(pm)
  call PMAuxiliaryInit(pm)
  
  PMAuxiliaryCreate => pm

end function PMAuxiliaryCreate
  
! ************************************************************************** !

subroutine PMAuxiliaryInit(this)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this  

  nullify(this%realization)
  nullify(this%comm1)
  this%ctype = ''
  
  call PMBaseInit(this)
  
end subroutine PMAuxiliaryInit

 
! ************************************************************************** !

function PMAuxiliaryCast(this)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_base_type), pointer :: this  

  class(pm_auxiliary_type), pointer :: PMAuxiliaryCast  
  
  nullify(PMAuxiliaryCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_auxiliary_type)
      PMAuxiliaryCast => this
    class default
      !geh: have default here to pass a null pointer if not of type ascii
  end select
  
end function PMAuxiliaryCast

! ************************************************************************** !

subroutine PMAuxiliarySetFunctionPointer(this,string)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Option_module
  
  implicit none
  
  class(pm_auxiliary_type) :: this
  character(len=MAXSTRINGLENGTH) :: string

  this%ctype = trim(string)
  select case(string)
    case('EVOLVING_STRATA')
      this%Evaluate => PMAuxiliaryEvolvingStrata
    case('SALINITY')
      this%Evaluate => PMAuxiliarySalinity
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMAuxiliarySetFunctionPointer.'
      call printErrMsg(this%option)
  end select
  
end subroutine PMAuxiliarySetFunctionPointer

! ************************************************************************** !

recursive subroutine PMAuxiliaryInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16 

  implicit none

  class(pm_auxiliary_type) :: this
  
  PetscReal :: time
  PetscErrorCode :: ierr
  
  time = 0.d0
  select case(this%ctype)
    case('EVOLVING_STRATA')
    case('SALINITY')
      
      call this%Evaluate(time,ierr)
  end select  

end subroutine PMAuxiliaryInitializeRun

! ************************************************************************** !

subroutine PMAuxiliaryEvolvingStrata(this,time,ierr)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Init_Subsurface_module

  implicit none
  
  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  call InitSubsurfAssignMatIDsToRegns(this%realization)
  call InitSubsurfAssignMatProperties(this%realization)
  
end subroutine PMAuxiliaryEvolvingStrata

! ************************************************************************** !

subroutine PMAuxiliarySalinity(this,time,ierr)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16
  !
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none
  
  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  PetscReal, parameter :: FMWNA = 22.989769d0
  PetscReal, parameter :: FMWCL = 35.4527d0
  PetscReal, parameter :: FMWNACL = FMWNA + FMWCL
  PetscInt i, j, nacl_id, na_id, cl_id
  PetscReal :: M_na, M_cl, M_h2o, xnacl
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt, parameter :: iphase = 1
    
  na_id = 0
  cl_id = 0
  nacl_id = 0
  M_na = 0.d0
  M_cl = 0.d0
  if (associated(this%realization%reaction%species_idx)) then
    na_id = this%realization%reaction%species_idx%na_ion_id
    cl_id = this%realization%reaction%species_idx%cl_ion_id
  else
    nacl_id = 1
  endif
  do j = 1, 2
    if (j == 1) then
      rt_auxvars => this%realization%patch%aux%RT%auxvars
      global_auxvars => this%realization%patch%aux%Global%auxvars
    else
      rt_auxvars => this%realization%patch%aux%RT%auxvars_bc
      global_auxvars => this%realization%patch%aux%Global%auxvars_bc
    endif
    do i = 1, size(rt_auxvars)
                                                   ! mol/L * g/mol = g/L and
      if (nacl_id > 0) then                        !   g/L => kg/m^3
        M_na = rt_auxvars(i)%total(nacl_id,iphase)*FMWNACL 
      else
        M_na = rt_auxvars(i)%total(na_id,iphase)*FMWNA
        M_cl = rt_auxvars(i)%total(cl_id,iphase)*FMWCL
      endif
      M_h2o = global_auxvars(i)%den_kg(iphase)  ! kg/m^3
      xnacl = (M_na + M_cl) / (M_na + M_cl + M_h2o)
      global_auxvars(i)%m_nacl(iphase) = xnacl
    enddo
  enddo
  
end subroutine PMAuxiliarySalinity

! ************************************************************************** !

subroutine PMAuxiliaryDestroy(this)
  ! 
  ! Destroys auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this
  
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  
end subroutine PMAuxiliaryDestroy

end module PM_Auxiliary_class
