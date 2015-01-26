module PM_Waste_Form_class

  use PM_Base_class
  use Realization_class
  use Option_module
  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type :: mpm_type
    PetscInt :: local_id
    type(point3d_type) :: coordinate
    PetscReal :: dose_rate
    PetscReal :: temperature
    PetscReal :: oxygen
    PetscReal :: carbonate
    PetscReal :: hydrogen
    PetscReal :: iron
  end type mpm_type
  
  type, public, extends(pm_base_type) :: pm_mpm_type
    class(realization_type), pointer :: realization
    type(mpm_type), pointer :: waste_form(:)
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: ioxygen
    PetscInt :: icarbonate
    PetscInt :: ihydrogen
    PetscInt :: iiron
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Init => PMWasteFormInit
     procedure, public :: Read => PMWasteFormRead
!    procedure, public :: SetupSolvers => PMWasteFormSetupSolvers
    procedure, public :: PMWasteFormSetRealization
!    procedure, public :: InitializeRun => PMWasteFormInitializeRun
!!    procedure, public :: FinalizeRun => PMWasteFormFinalizeRun
!    procedure, public :: InitializeTimestep => PMWasteFormInitializeTimestep
!    procedure, public :: FinalizeTimestep => PMWasteFormFinalizeTimestep
!    procedure, public :: PreSolve => PMWasteFormPreSolve
    procedure, public :: Solve => PMWasteFormSolve
!    procedure, public :: PostSolve => PMWasteFormPostSolve
!    procedure, public :: AcceptSolution => PMWasteFormAcceptSolution
!    procedure, public :: TimeCut => PMWasteFormTimeCut
!    procedure, public :: UpdateSolution => PMWasteFormUpdateSolution
!    procedure, public :: UpdateAuxvars => PMWasteFormUpdateAuxvars
!    procedure, public :: Checkpoint => PMWasteFormCheckpoint    
!    procedure, public :: Restart => PMWasteFormRestart  
!    procedure, public :: Destroy => PMWasteFormDestroy
  end type pm_mpm_type
  
  public :: PMWasteFormCreate, &
            PMWasteFormInit !, &
!            PMWasteFormSetupSolvers, &
!            PMWasteFormInitializeTimestepA, &
!            PMWasteFormInitializeTimestepB, &
!            PMWasteFormInitializeRun, &
!            PMWasteFormUpdateSolution, &
!            PMWasteFormUpdatePropertiesNI, &
!            PMWasteFormTimeCut, &
!            PMWasteFormDestroy
  
contains

! ************************************************************************** !

function PMWasteFormCreate()
  ! 
  ! Intializes shared members of subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type), pointer :: PMWasteFormCreate
  
  allocate(PMWasteFormCreate)
  nullify(PMWasteFormCreate%realization)
  nullify(PMWasteFormCreate%waste_form)
  nullify(PMWasteFormCreate%coordinates)
  PMWasteFormCreate%ioxygen = UNINITIALIZED_INTEGER
  PMWasteFormCreate%icarbonate = UNINITIALIZED_INTEGER
  PMWasteFormCreate%ihydrogen = UNINITIALIZED_INTEGER
  PMWasteFormCreate%iiron = UNINITIALIZED_INTEGER

  call PMBaseCreate(PMWasteFormCreate)

end function PMWasteFormCreate

! ************************************************************************** !

subroutine PMWasteFormRead(this,input,option)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  use Input_Aux_module
  use String_module
  use Utility_module
  
  implicit none
  
  class(pm_mpm_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  type(point3d_type) :: coordinate
  PetscReal, pointer :: x(:), y(:), z(:)
  PetscInt :: coordinate_size = 100
  PetscInt :: num_coordinates, dummy_int
  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'MPM'
  
  allocate(x(coordinate_size))
  x = 0.d0
  allocate(y(coordinate_size))
  y = 0.d0
  allocate(z(coordinate_size))
  z = 0.d0
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
    
      case('COORDINATES')
        num_coordinates = 0
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          num_coordinates = num_coordinates + 1
          call GeometryReadCoordinate(input,option,coordinate,error_string)
          if (num_coordinates > coordinate_size) then
            ! use dummy_int on the first two since the value is doubled on
            ! every call.
            dummy_int = coordinate_size
            call reallocateRealArray(x,dummy_int)
            dummy_int = coordinate_size
            call reallocateRealArray(y,dummy_int)
            call reallocateRealArray(z,coordinate_size)
          endif
          x(num_coordinates) = coordinate%x
          y(num_coordinates) = coordinate%y
          z(num_coordinates) = coordinate%z
        enddo
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  allocate(this%coordinates(num_coordinates))
  do i = 1, num_coordinates
    this%coordinates(i)%x = x(i)
    this%coordinates(i)%y = y(i)
    this%coordinates(i)%z = z(i)
  enddo
  deallocate(x,y,z)
  nullify(x,y,z)  
  
end subroutine PMWasteFormRead

! ************************************************************************** !

subroutine PMWasteFormInit(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Reaction_Aux_module
  use Option_module

  implicit none
  
  class(pm_mpm_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  PetscInt :: num_local_coordinates
  PetscInt, allocatable :: cell_ids(:,:)
  PetscInt :: icoord, i, j, k, local_id
  PetscErrorCode :: ierr
  
  grid => this%realization%patch%grid
  option => this%realization%option
  
  allocate(cell_ids(2,size(this%coordinates)))
  num_local_coordinates = 0
  do icoord = 1, size(this%coordinates)
    local_id = -1
    select case(grid%itype)
        case(STRUCTURED_GRID)
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                              this%coordinates(icoord)%x, &
                                              this%coordinates(icoord)%y, &
                                              this%coordinates(icoord)%z, &
                                              i,j,k)
          if (i > 0 .and. j > 0 .and. k > 0) then
            local_id = i + (j-1)*grid%structured_grid%nlx + &
                       (k-1)*grid%structured_grid%nlxy
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call UGridGetCellFromPoint(this%coordinates(icoord)%x, &
                                     this%coordinates(icoord)%y, &
                                     this%coordinates(icoord)%z, &
                                     grid%unstructured_grid,option,local_id)
        case default
           option%io_buffer = 'Only STRUCTURED_GRID and ' // &
             'IMPLICIT_UNSTRUCTURED_GRID types supported in PMWasteForm.'
           call printErrMsg(option)
      end select
      if (local_id > 0) then
        num_local_coordinates = num_local_coordinates + 1
        cell_ids(1,num_local_coordinates) = i
        cell_ids(2,num_local_coordinates) = local_id
      endif
  enddo
  allocate(this%waste_form(num_local_coordinates))
  do i = 1, num_local_coordinates
    this%waste_form(i)%local_id = cell_ids(2,i)
    this%waste_form(i)%coordinate%x = this%coordinates(cell_ids(1,i))%x
    this%waste_form(i)%coordinate%y = this%coordinates(cell_ids(1,i))%y
    this%waste_form(i)%coordinate%z = this%coordinates(cell_ids(1,i))%z
  enddo
  deallocate(cell_ids)
  deallocate(this%coordinates)
  nullify(this%coordinates)
  
  ! set up indexing of solute concentrations
  species_name = 'O2(aq)'
  this%ioxygen = GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'HCO3-'
  this%icarbonate = GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'H2(aq)'
  this%ihydrogen = GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'Fe++'
  this%iiron = GetPrimarySpeciesIDFromName(species_name,reaction,option)
  
end subroutine PMWasteFormInit

! ************************************************************************** !

subroutine PMWasteFormSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Realization_class

  implicit none
  
  class(pm_mpm_type) :: this
  class(realization_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWasteFormSetRealization

! ************************************************************************** !

recursive subroutine PMWasteFormInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15 

  implicit none
  
  class(pm_mpm_type) :: this

end subroutine PMWasteFormInitializeRun

! ************************************************************************** !

subroutine PMWasteFormInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Global_module
  
  implicit none
  
  class(pm_mpm_type) :: this


end subroutine PMWasteFormInitializeTimestep

! ************************************************************************** !

subroutine PMWasteFormPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Grid_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_mpm_type) :: this
  
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: i
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
  do i = 1, size(this%waste_form)
    local_id = this%waste_form(i)%local_id
    ghosted_id = grid%nL2G(local_id)
    this%waste_form(i)%temperature = global_auxvars(ghosted_id)%temp
    this%waste_form(i)%oxygen = rt_auxvars(ghosted_id)%total(this%ioxygen,1)
    this%waste_form(i)%carbonate = &
      rt_auxvars(ghosted_id)%total(this%icarbonate,1)
    this%waste_form(i)%hydrogen = &
      rt_auxvars(ghosted_id)%total(this%ihydrogen,1)
    this%waste_form(i)%iron = rt_auxvars(ghosted_id)%total(this%iiron,1)
  enddo
  
end subroutine PMWasteFormPreSolve

! ************************************************************************** !

subroutine PMWasteFormSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  !
  use Argonne_Mixed_Potential_module
  
  implicit none
  
  class(pm_mpm_type) :: this
  PetscReal :: time
  PetscInt :: ierr
  
  PetscReal :: conc(4)
  
  PetscInt :: i
  
  do i = 1, size(this%waste_form)
    call AMPRun(time,conc)
  enddo

end subroutine PMWasteFormSolve

! ************************************************************************** !

subroutine PMWasteFormPostSolve(this)
  ! 
  ! PMWasteFormUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_mpm_type) :: this
  
  this%option%io_buffer = 'PMWasteFormPostSolve() must be extended.'
  call printErrMsg(this%option)  
  
end subroutine PMWasteFormPostSolve

! ************************************************************************** !

function PMWasteFormAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this
  
  PetscBool :: PMWasteFormAcceptSolution
  
  ! do nothing
  PMWasteFormAcceptSolution = PETSC_TRUE
  
end function PMWasteFormAcceptSolution

! ************************************************************************** !

subroutine PMWasteFormUpdatePropertiesTS(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this
  
!  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMWasteFormUpdatePropertiesTS

! ************************************************************************** !

subroutine PMWasteFormTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15 
  
  implicit none
  
  class(pm_mpm_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMWasteFormTimeCut

! ************************************************************************** !

subroutine PMWasteFormFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this

end subroutine PMWasteFormFinalizeTimestep

! ************************************************************************** !

subroutine PMWasteFormUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMWasteFormUpdateSolution  

! ************************************************************************** !

subroutine PMWasteFormUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this

  this%option%io_buffer = 'PMWasteFormUpdateAuxvars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMWasteFormUpdateAuxvars   

! ************************************************************************** !

subroutine PMWasteFormCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_mpm_type) :: this
  PetscViewer :: viewer
  
end subroutine PMWasteFormCheckpoint

! ************************************************************************** !

subroutine PMWasteFormRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_mpm_type) :: this
  PetscViewer :: viewer
  
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMWasteFormRestart

! ************************************************************************** !

recursive subroutine PMWasteFormFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWasteFormFinalizeRun

! ************************************************************************** !

subroutine PMWasteFormDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_mpm_type) :: this
  
end subroutine PMWasteFormDestroy
  
end module PM_Waste_Form_class
