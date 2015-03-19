module PM_Waste_Form_class

  use PM_Base_class
  use Realization_class
  use Option_module
  use Geometry_module
  use Mass_Transfer2_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type :: mpm_type
    PetscInt :: local_id
    type(point3d_type) :: coordinate
    PetscReal :: temperature
    PetscReal, pointer :: concentration(:,:)
    PetscReal, pointer :: flux(:,:)
  end type mpm_type
  
  type, public, extends(pm_base_type) :: pm_mpm_type
    class(realization_type), pointer :: realization
    type(mpm_type), pointer :: waste_form(:)
    type(point3d_type), pointer :: coordinate(:)
    PetscInt :: num_grid_cells_in_waste_form
    PetscInt :: num_waste_forms
    ! mapping of mpm species into mpm concentration array
    PetscInt, pointer :: mapping_mpm(:)
    ! mapping of species in mpm concentration array to pflotran
    PetscInt, pointer :: mapping_mpm_to_pflotran(:)
    PetscInt :: iUO2_2p
    PetscInt :: iUCO3_2n
    PetscInt :: iUO2
    PetscInt :: iCO3_2n
    PetscInt :: iO2
    PetscInt :: iH2O2
    PetscInt :: iFe_2p
    PetscInt :: iH2
    PetscInt :: iUO2_sld
    PetscInt :: iUO3_sld
    PetscInt :: iUO4_sld
    PetscInt :: num_concentrations
    type(mass_transfer2_type), pointer :: mass_transfer
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Init => PMWasteFormInit
    procedure, public :: Read => PMWasteFormRead
!    procedure, public :: SetupSolvers => PMWasteFormSetupSolvers
    procedure, public :: PMWasteFormSetRealization
    procedure, public :: InitializeRun => PMWasteFormInitializeRun
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
    procedure, public :: Checkpoint => PMWasteFormCheckpoint    
    procedure, public :: Restart => PMWasteFormRestart  
    procedure, public :: Destroy => PMWasteFormDestroy
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
  nullify(PMWasteFormCreate%coordinate)
  PMWasteFormCreate%num_grid_cells_in_waste_form = UNINITIALIZED_INTEGER
  PMWasteFormCreate%num_waste_forms = UNINITIALIZED_INTEGER
  PMWasteFormCreate%num_concentrations = 11
  PMWasteFormCreate%iUO2_2p = 1
  PMWasteFormCreate%iUCO3_2n = 2
  PMWasteFormCreate%iUO2 = 3
  PMWasteFormCreate%iCO3_2n = 4
  PMWasteFormCreate%iO2 = 5
  PMWasteFormCreate%iH2O2 = 6
  PMWasteFormCreate%iFe_2p = 7
  PMWasteFormCreate%iH2 = 8
  PMWasteFormCreate%iUO2_sld = 9
  PMWasteFormCreate%iUO3_sld = 10
  PMWasteFormCreate%iUO4_sld = 11
  nullify(PMWasteFormCreate%mass_transfer)
  allocate(PMWasteFormCreate%mapping_mpm_to_pflotran( &
             PMWasteFormCreate%num_concentrations))
  PMWasteFormCreate%mapping_mpm_to_pflotran = UNINITIALIZED_INTEGER
  allocate(PMWasteFormCreate%mapping_mpm(4))
  PMWasteFormCreate%mapping_mpm = [PMWasteFormCreate%iO2, &
                                   PMWasteFormCreate%iCO3_2n, &
                                   PMWasteFormCreate%iH2, &
                                   PMWasteFormCreate%iFe_2p]
  
  call PMBaseCreate(PMWasteFormCreate)

end function PMWasteFormCreate

! ************************************************************************** !

subroutine PMWasteFormRead(this,input)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  
  implicit none
  
  class(pm_mpm_type) :: this
  type(input_type) :: input
  
  type(point3d_type) :: coordinate
  PetscReal, pointer :: x(:), y(:), z(:)
  PetscInt :: coordinate_size = 100
  PetscInt :: num_coordinates, dummy_int
  PetscInt :: i
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => this%option
  
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
    
      case('NUM_GRID_CELLS')
        call InputReadInt(input,option,this%num_grid_cells_in_waste_form)
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
  this%num_waste_forms = num_coordinates
  allocate(this%coordinate(num_coordinates))
  do i = 1, num_coordinates
    this%coordinate(i)%x = x(i)
    this%coordinate(i)%y = y(i)
    this%coordinate(i)%z = z(i)
  enddo
  deallocate(x,y,z)
  nullify(x,y,z)
  
  if (this%num_grid_cells_in_waste_form == UNINITIALIZED_INTEGER) then
    option%io_buffer = &
      'NUM_GRID_CELLS must be specified for mixed potential model.'
    call printErrMsg(option)
  endif
  
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
  reaction => this%realization%reaction
  
  allocate(cell_ids(2,size(this%coordinate)))
  num_local_coordinates = 0
  do icoord = 1, size(this%coordinate)
    local_id = -1
    select case(grid%itype)
        case(STRUCTURED_GRID)
          call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                              this%coordinate(icoord)%x, &
                                              this%coordinate(icoord)%y, &
                                              this%coordinate(icoord)%z, &
                                              i,j,k)
          if (i > 0 .and. j > 0 .and. k > 0) then
            local_id = i + (j-1)*grid%structured_grid%nlx + &
                       (k-1)*grid%structured_grid%nlxy
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call UGridGetCellFromPoint(this%coordinate(icoord)%x, &
                                     this%coordinate(icoord)%y, &
                                     this%coordinate(icoord)%z, &
                                     grid%unstructured_grid,option,local_id)
        case default
           option%io_buffer = 'Only STRUCTURED_GRID and ' // &
             'IMPLICIT_UNSTRUCTURED_GRID types supported in PMWasteForm.'
           call printErrMsg(option)
      end select
      if (local_id > 0) then
        num_local_coordinates = num_local_coordinates + 1
        cell_ids(1,num_local_coordinates) = icoord
        cell_ids(2,num_local_coordinates) = local_id
      endif
  enddo
  allocate(this%waste_form(num_local_coordinates))
  do i = 1, num_local_coordinates
    allocate(this%waste_form(i)%concentration(this%num_concentrations, &
                                            this%num_grid_cells_in_waste_form))
                                              
!geh    this%waste_form(i)%concentration = UNINITIALIZED_DOUBLE
    this%waste_form(i)%concentration = 1.d-20
    
    allocate(this%waste_form(i)%flux(this%num_concentrations,1))
    this%waste_form(i)%flux = UNINITIALIZED_DOUBLE
    this%waste_form(i)%local_id = cell_ids(2,i)
    this%waste_form(i)%coordinate%x = this%coordinate(cell_ids(1,i))%x
    this%waste_form(i)%coordinate%y = this%coordinate(cell_ids(1,i))%y
    this%waste_form(i)%coordinate%z = this%coordinate(cell_ids(1,i))%z
  enddo
  deallocate(cell_ids)
  deallocate(this%coordinate)
  nullify(this%coordinate)
  
  ! set up indexing of solute concentrations
  species_name = 'O2(aq)'
  this%mapping_mpm_to_pflotran(this%iO2) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'HCO3-'
  this%mapping_mpm_to_pflotran(this%iCO3_2n) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'H2(aq)'
  this%mapping_mpm_to_pflotran(this%iH2) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'Fe++'
  this%mapping_mpm_to_pflotran(this%iFe_2p) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  
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
  use Communicator_Base_module
  
  implicit none
  
  class(pm_mpm_type) :: this
  
  PetscInt :: num_waste_form_cells
  PetscInt :: i
  PetscInt, allocatable :: waste_form_cell_ids(:)
  PetscErrorCode :: ierr
  
#ifdef PM_WP_DEBUG  
  call printMsg(this%option,'PMRT%InitializeRun()')
#endif
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif
  
  ! set up mass transfer
  this%mass_transfer => MassTransfer2Create()
  ! create a Vec sized by # waste packages * # primary dofs influenced by 
  ! waste package
  num_waste_form_cells = size(this%waste_form)
  call VecCreateMPI(this%option%mycomm,num_waste_form_cells, &
                    this%num_waste_forms,this%mass_transfer%vec, &
                    ierr);CHKERRQ(ierr)

  allocate(waste_form_cell_ids(num_waste_form_cells))
  waste_form_cell_ids = 0
  do i = 1, num_waste_form_cells
    waste_form_cell_ids(i) = this%waste_form(i)%local_id
  enddo
  call this%realization%comm1%AONaturalToPetsc(waste_form_cell_ids)
  call VecScatterCreate(this%mass_transfer%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%work,waste_form_cell_ids, &
                        this%mass_transfer%scatter_ctx,ierr);CHKERRQ(ierr)
  deallocate(waste_form_cell_ids)

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
  PetscInt :: iwasteform
  PetscInt :: i
  PetscInt :: icomp_mpm
  PetscInt :: icomp_pflotran
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
  do iwasteform = 1, size(this%waste_form)
    local_id = this%waste_form(iwasteform)%local_id
    ghosted_id = grid%nL2G(local_id)
    this%waste_form(iwasteform)%temperature = global_auxvars(ghosted_id)%temp
    ! overwrite the components in this%mapping_pflotran array
    do i = 1, size(this%mapping_mpm)
      icomp_mpm = this%mapping_mpm(i)
      icomp_pflotran = this%mapping_mpm_to_pflotran(icomp_mpm)
      this%waste_form(iwasteform)%concentration(icomp_mpm,1) = &
        ! the 1 in the second index if for the liquid phase
        rt_auxvars(ghosted_id)%total(icomp_pflotran,1)
    enddo
  enddo
  
end subroutine PMWasteFormPreSolve

! ************************************************************************** !

subroutine PMWasteFormSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  !
!  use Argonne_Mixed_Potential_module
  
  implicit none

  interface
    subroutine AMP_step ( sTme, conc, initialRun, flux, status )
      real ( kind = 8), intent( in )  :: sTme   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      real ( kind = 8), intent(out), dimension (:,:) :: flux
      integer ( kind = 4), intent(out) :: status
    end subroutine
  end interface  

  class(pm_mpm_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  integer ( kind = 4) :: status
  logical ( kind = 4) :: initialRun = PETSC_FALSE
  
  PetscInt :: i
  
  ierr = 0
  call PMWasteFormPreSolve(this)
  do i = 1, size(this%waste_form)
#ifdef MPM_MODEL  
    call AMP_step(time, this%waste_form(i)%concentration, initialRun, &
                  this%waste_form(i)%flux, status)
#endif
    if (status == 0) then
      ierr = 1
      exit
    endif
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
  
  PetscReal, pointer :: vec_p(:)
  PetscErrorCode :: ierr
  PetscInt :: i
  
  call VecGetArrayF90(this%mass_transfer%vec,vec_p,ierr);CHKERRQ(ierr)
  do i = 1, size(this%waste_form)
    ! do something here.
  enddo
  call VecRestoreArrayF90(this%mass_transfer%vec,vec_p,ierr);CHKERRQ(ierr)
  
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
  use Utility_module

  implicit none
  
  class(pm_mpm_type) :: this
  
  PetscInt :: i
  
  do i = 1, size(this%waste_form)
    call DeallocateArray(this%waste_form(i)%concentration)
    call DeallocateArray(this%waste_form(i)%flux)
  enddo
  call DeallocateArray(this%mapping_mpm_to_pflotran)
  call DeallocateArray(this%mapping_mpm)
  call MassTransfer2Destroy(this%mass_transfer)
  deallocate(this%waste_form)
  nullify(this%waste_form)
  
end subroutine PMWasteFormDestroy
  
end module PM_Waste_Form_class
