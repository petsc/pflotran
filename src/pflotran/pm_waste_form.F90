module PM_Waste_Form_class

  use PM_Base_class
  use Realization_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  use Dataset_Base_class
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  PetscBool, public :: bypass_warning_message = PETSC_FALSE

  type :: fmdm_type
    PetscInt :: local_id
    type(point3d_type) :: coordinate
    PetscReal :: temperature
    PetscReal, pointer :: concentration(:,:)
    PetscReal :: fuel_dissolution_rate
    PetscReal :: specific_surface_area
    PetscReal :: volume
    PetscReal :: burnup
    type(fmdm_type), pointer :: next
  end type fmdm_type
  
  type, public, extends(pm_base_type) :: pm_waste_form_type
    class(realization_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: data_mediator_species
    class(data_mediator_vec_type), pointer :: data_mediator
  contains
    procedure, public :: PMWasteFormSetRealization
  end type pm_waste_form_type
  
  type :: glass_type
    PetscInt :: local_id
    type(point3d_type) :: coordinate
    PetscReal :: exposure_factor
    PetscReal :: volume
    PetscReal :: fuel_dissolution_rate
    type(glass_type), pointer :: next
  end type glass_type
  
  type, public, extends(pm_waste_form_type) :: pm_glass_type
    type(glass_type), pointer :: waste_form_list
    PetscReal :: specific_surface_area
    PetscReal :: formula_weight
    class(dataset_base_type), pointer ::  mass_fraction_dataset
    PetscReal :: glass_density
  contains
    procedure, public :: Setup => PMGlassSetup
    procedure, public :: Read => PMGlassRead
    procedure, public :: InitializeRun => PMGlassInitializeRun
    procedure, public :: InitializeTimestep => PMGlassInitializeTimestep
    procedure, public :: FinalizeTimestep => PMGlassFinalizeTimestep
    procedure, public :: UpdateSolution => PMGlassUpdateSolution
    procedure, public :: Solve => PMGlassSolve
    procedure, public :: Checkpoint => PMGlassCheckpoint    
    procedure, public :: Restart => PMGlassRestart  
    procedure, public :: Destroy => PMGlassDestroy
  end type pm_glass_type
  
  type, public, extends(pm_waste_form_type) :: pm_fmdm_type
    type(fmdm_type), pointer :: waste_form_list
    PetscInt :: num_grid_cells_in_waste_form
    ! mapping of fmdm species into fmdm concentration array
    PetscInt, pointer :: mapping_fmdm(:)
    ! mapping of species in fmdm concentration array to pflotran
    PetscInt, pointer :: mapping_fmdm_to_pflotran(:)
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
    PetscBool :: initialized
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Setup => PMFMDMSetup
    procedure, public :: Read => PMFMDMRead
!    procedure, public :: SetupSolvers => PMFMDMSetupSolvers
    procedure, public :: InitializeRun => PMFMDMInitializeRun
!!    procedure, public :: FinalizeRun => PMFMDMFinalizeRun
    procedure, public :: InitializeTimestep => PMFMDMInitializeTimestep
    procedure, public :: FinalizeTimestep => PMFMDMFinalizeTimestep
!    procedure, public :: PreSolve => PMFMDMPreSolve
    procedure, public :: Solve => PMFMDMSolve
!    procedure, public :: PostSolve => PMFMDMPostSolve
!    procedure, public :: AcceptSolution => PMFMDMAcceptSolution
!    procedure, public :: TimeCut => PMFMDMTimeCut
!    procedure, public :: UpdateSolution => PMFMDMUpdateSolution
!    procedure, public :: UpdateAuxVars => PMFMDMUpdateAuxVars
    procedure, public :: Checkpoint => PMFMDMCheckpoint    
    procedure, public :: Restart => PMFMDMRestart  
    procedure, public :: Destroy => PMFMDMDestroy
  end type pm_fmdm_type
  
  public :: PMFMDMCreate, &
            PMFMDMSetup, &
            PMGlassCreate, &
            PMGlassSetup
  
contains

! ************************************************************************** !

subroutine PMWasteFormInit(this)
  ! 
  ! Creates the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15, 07/20/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
  call PMBaseInit(this)
  nullify(this%realization)
  nullify(this%data_mediator)
  this%data_mediator_species = ''

end subroutine PMWasteFormInit

! ************************************************************************** !

subroutine PMWasteFormReadSelectCase(this,input,keyword,found,error_string, &
                                     option)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15, 07/20/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type) :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option  

  found = PETSC_TRUE
  select case(trim(keyword))
    
    case('DATA_MEDIATOR_SPECIES')
      call InputReadWord(input,option,this%data_mediator_species,PETSC_TRUE)
      call InputErrorMsg(input,option,'data_mediator_species',error_string)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMWasteFormReadSelectCase

! ************************************************************************** !

subroutine PMWasteFormSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Realization_class

  implicit none
  
  class(pm_waste_form_type) :: this
  class(realization_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWasteFormSetRealization

! ************************************************************************** !

subroutine PMWasteFormStrip(this)
  ! 
  ! Destroys a waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_waste_form_type) :: this

  nullify(this%realization)
  nullify(this%data_mediator)
  
end subroutine PMWasteFormStrip

! ************************************************************************** !

subroutine GlassInit(glass)
  ! 
  ! Initializes the fuel matrix degradation model waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  type(glass_type) :: glass

  glass%local_id = UNINITIALIZED_INTEGER
  glass%coordinate%x = UNINITIALIZED_DOUBLE
  glass%coordinate%y = UNINITIALIZED_DOUBLE
  glass%coordinate%z = UNINITIALIZED_DOUBLE
  glass%exposure_factor = UNINITIALIZED_DOUBLE
  glass%fuel_dissolution_rate = UNINITIALIZED_DOUBLE
  glass%volume = UNINITIALIZED_DOUBLE
  nullify(glass%next)

end subroutine GlassInit

! ************************************************************************** !

function PMGlassCreate()
  ! 
  ! Creates the Glass waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  class(pm_glass_type), pointer :: PMGlassCreate
  
  allocate(PMGlassCreate)
  call PMWasteFormInit(PMGlassCreate)
  nullify(PMGlassCreate%waste_form_list)
  PMGlassCreate%specific_surface_area = UNINITIALIZED_DOUBLE
  PMGlassCreate%formula_weight = UNINITIALIZED_DOUBLE
  nullify(PMGlassCreate%mass_fraction_dataset)
  PMGlassCreate%glass_density = 2.65d3 ! kg/m^3

end function PMGlassCreate

! ************************************************************************** !

subroutine PMGlassRead(this,input)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class
  
  implicit none
  
  class(pm_glass_type) :: this
  type(input_type) :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  type(glass_type), pointer :: new_waste_form, prev_waste_form
  class(dataset_ascii_type), pointer :: dataset_ascii

  option => this%option
  
  error_string = 'GLASS'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    call PMWasteFormReadSelectCase(this,input,word,found,error_string,option)
    
    if (found) cycle
    
    select case(trim(word))
    
      case('WASTE_FORM')
        error_string = 'GLASS,WASTE_FORM'
        allocate(new_waste_form)
        call GlassInit(new_waste_form)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
            case('EXPOSURE_FACTOR')
              call InputReadDouble(input,option,new_waste_form%exposure_factor)
              call InputErrorMsg(input,option,'exposure factor',error_string)
            case('VOLUME')
              call InputReadDouble(input,option,new_waste_form%volume)
              call InputErrorMsg(input,option,'volume',error_string)
            case('COORDINATE')
              call GeometryReadCoordinate(input,option, &
                                          new_waste_form%coordinate, &
                                          error_string)
          end select
        enddo
        if (Uninitialized(new_waste_form%volume)) then
          option%io_buffer = &
            'VOLUME must be specified for all glass ' // &
            'waste packages.'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_waste_form%exposure_factor)) then
          option%io_buffer = &
            'EXPOSURE_FACTOR must be specified for all glass ' // &
            'waste packages.'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_waste_form%coordinate%z)) then
          option%io_buffer = &
            'COORDINATE must be specified for all glass ' // &
            'waste packages.'
          call printErrMsg(option)
        endif
        if (.not.associated(this%waste_form_list)) then
          this%waste_form_list => new_waste_form
          prev_waste_form => new_waste_form
        else
          prev_waste_form%next => new_waste_form
          prev_waste_form => new_waste_form
        endif
        nullify(new_waste_form)
        error_string = 'GLASS'
      case('SPECIFIC_SURFACE_AREA')
        call InputReadDouble(input,option,this%specific_surface_area)
        call InputErrorMsg(input,option,'specific surface area',error_string)
      case('FORMULA_WEIGHT')
        call InputReadDouble(input,option,this%formula_weight)
        call InputErrorMsg(input,option,'formula_weight',error_string)
      case('MASS_FRACTION')
        dataset_ascii => DatasetAsciiCreate()
        dataset_ascii%array_width = 1 ! does not include time
        this%mass_fraction_dataset => dataset_ascii
        dataset_ascii%data_type = DATASET_REAL
        call ConditionReadValues(input,option,word, &
                                 this%mass_fraction_dataset,units)
        if (associated(dataset_ascii%time_storage)) then
          ! default time interpolation is linear
          if (dataset_ascii%time_storage%time_interpolation_method == &
              INTERPOLATION_NULL) then
            dataset_ascii%time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
          endif
        endif
      case('GLASS_DENSITY')
        call InputReadDouble(input,option,this%glass_density)
        call InputErrorMsg(input,option,'glass density',error_string)
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
  if (Uninitialized(this%specific_surface_area)) then
    option%io_buffer = 'SPECIFIC_SURFACE_AREA must be specified in ' // &
      trim(error_string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%formula_weight)) then
    option%io_buffer = 'FORMULA_WEIGHT must be specified in ' // &
      trim(error_string)
    call printErrMsg(option)
  endif
  if (.not.associated(this%mass_fraction_dataset)) then
    option%io_buffer = 'MASS_FRACTION must be specified in ' // &
      trim(error_string)
    call printErrMsg(option)
  endif
  
  if (len_trim(this%data_mediator_species) == 0) then
    option%io_buffer = &
      'DATA_MEDIATOR_SPECIES must be specified must be specified in ' // &
      trim(error_string)
  endif
  
end subroutine PMGlassRead

! ************************************************************************** !

subroutine PMGlassSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Reaction_Aux_module
  use Option_module

  implicit none
  
  class(pm_glass_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  type(glass_type), pointer :: cur_waste_form, prev_waste_form, next_waste_form
  PetscInt :: i, j, k, local_id
  PetscErrorCode :: ierr
  
  grid => this%realization%patch%grid
  option => this%realization%option
  reaction => this%realization%reaction
  
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    local_id = -1
    select case(grid%itype)
      case(STRUCTURED_GRID)
        call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                            cur_waste_form%coordinate%x, &
                                            cur_waste_form%coordinate%y, &
                                            cur_waste_form%coordinate%z, &
                                            i,j,k)
        if (i > 0 .and. j > 0 .and. k > 0) then
          local_id = i + (j-1)*grid%structured_grid%nlx + &
                      (k-1)*grid%structured_grid%nlxy
        endif
      case(IMPLICIT_UNSTRUCTURED_GRID)
        call UGridGetCellFromPoint(cur_waste_form%coordinate%x, &
                                   cur_waste_form%coordinate%y, &
                                   cur_waste_form%coordinate%z, &
                                   grid%unstructured_grid,option,local_id)
      case default
          option%io_buffer = 'Only STRUCTURED_GRID and ' // &
            'IMPLICIT_UNSTRUCTURED_GRID types supported in PMGlass.'
          call printErrMsg(option)
    end select
    if (local_id > 0) then
      cur_waste_form%local_id = local_id
      prev_waste_form => cur_waste_form
      cur_waste_form => cur_waste_form%next
    else
      ! remove waste form
      next_waste_form => cur_waste_form%next
      if (associated(prev_waste_form)) then
        prev_waste_form%next => next_waste_form
      else
        this%waste_form_list => next_waste_form
      endif
      deallocate(cur_waste_form)
      cur_waste_form => next_waste_form
    endif
  enddo
  
end subroutine PMGlassSetup

! ************************************************************************** !

recursive subroutine PMGlassInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  use Communicator_Base_module
  use Reaction_Aux_module
  use Realization_Base_class
  use Data_Mediator_Vec_class
  use Dataset_module
  use Time_Storage_module
  
  implicit none

#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pm_glass_type) :: this
  
  IS :: is
  type(glass_type), pointer :: cur_waste_form
  PetscInt :: num_waste_form_cells
  PetscInt :: i
  PetscInt :: data_mediator_species_id
  PetscInt, allocatable :: waste_form_cell_ids(:)
  PetscReal :: time
  type(time_storage_type), pointer :: null_time_storage
  PetscErrorCode :: ierr
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif

  data_mediator_species_id = &
    GetPrimarySpeciesIDFromName(this%data_mediator_species, &
                                this%realization%reaction,this%option)
  nullify(null_time_storage)
  call DatasetVerify(this%mass_fraction_dataset,null_time_storage,this%option)
  ! set up mass transfer
  call RealizCreateTranMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%tran_data_mediator_list)
  ! create a Vec sized by # waste packages * # primary dofs influenced by 
  ! waste package
  ! count of waste form cells
  cur_waste_form => this%waste_form_list
  num_waste_form_cells = 0
  do
    if (.not.associated(cur_waste_form)) exit
    num_waste_form_cells = num_waste_form_cells + 1
    cur_waste_form => cur_waste_form%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,num_waste_form_cells, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(waste_form_cell_ids(num_waste_form_cells))
    waste_form_cell_ids = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      i = i + 1
      waste_form_cell_ids(i) = cur_waste_form%local_id
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    waste_form_cell_ids(:) = waste_form_cell_ids(:) - 1
    waste_form_cell_ids(:) = waste_form_cell_ids(:) + &
                             this%realization%patch%grid%global_offset
    waste_form_cell_ids(:) = waste_form_cell_ids(:) * this%option%ntrandof
    waste_form_cell_ids(:) = waste_form_cell_ids(:)  + &
                             data_mediator_species_id - 1
  endif
  call ISCreateGeneral(this%option%mycomm,num_waste_form_cells, &
                       waste_form_cell_ids,PETSC_COPY_VALUES,is, &
                       ierr);CHKERRQ(ierr)
  if (allocated(waste_form_cell_ids)) deallocate(waste_form_cell_ids)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  
  call PMGlassSolve(this,0.d0,ierr)
  
end subroutine PMGlassInitializeRun

! ************************************************************************** !

subroutine PMGlassInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15
  use Dataset_module
  
  implicit none
  
  class(pm_glass_type) :: this
  
  type(glass_type), pointer :: cur_waste_form
  PetscReal :: rate
  PetscReal :: dV
  PetscReal, parameter :: conversion = 1.d0/(24.d0*3600.d0)
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GLASS MODEL ",65("="))')
  endif
  
  ! update mass balances after transport step
  cur_waste_form => this%waste_form_list
  do 
    if (.not.associated(cur_waste_form)) exit
            ! kg glass/sec
    rate = cur_waste_form%fuel_dissolution_rate * &  ! kg glass/m^2/day
            this%specific_surface_area * &            ! m^2/kg glass
            cur_waste_form%exposure_factor * &        ! [-]
            cur_waste_form%volume * &                 ! m^3 glass
            this%glass_density * &                    ! kg glass/m^3 glass
            conversion                                ! 1/day -> 1/sec  
         ! kg glass/sec / kg glass/m^3 glass = m^3 glass
    dV = rate / this%glass_density * this%option%tran_dt
    cur_waste_form%volume = cur_waste_form%volume - dV
    cur_waste_form => cur_waste_form%next
  enddo
  
  call DatasetUpdate(this%mass_fraction_dataset,this%option%time,this%option)

end subroutine PMGlassInitializeTimestep

! ************************************************************************** !

subroutine PMGlassPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  class(pm_glass_type) :: this
  
end subroutine PMGlassPreSolve

! ************************************************************************** !

subroutine PMGlassSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15
  !
  use Grid_module
  use Global_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pm_glass_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(glass_type), pointer :: cur_waste_form
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id  
  PetscInt :: i
  PetscReal, pointer :: vec_p(:)       ! kmol/day -> mol/sec
  PetscReal, parameter :: conversion = 1.d3/(24.d0*3600.d0)

  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  cur_waste_form => this%waste_form_list
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    i = i + 1
    local_id = cur_waste_form%local_id
    ghosted_id = grid%nL2G(local_id)
    if (cur_waste_form%volume > 0.d0) then
      cur_waste_form%fuel_dissolution_rate = & ! kg glass/m^2/day
        560.d0*exp(-7397.d0/(global_auxvars(ghosted_id)%temp+273.15d0))
      ! mol/sec
      vec_p(i) = cur_waste_form%fuel_dissolution_rate * &  ! kg glass/m^2/day
                 this%formula_weight * &                   ! kmol radionuclide/kg radionuclide
                 this%mass_fraction_dataset%rarray(1) * &  ! kg radionuclude/kg glass
                 this%specific_surface_area * &            ! m^2/kg glass
                 cur_waste_form%exposure_factor * &        ! [-]
                 cur_waste_form%volume * &                 ! m^3 glass
                 this%glass_density * &                    ! kg glass/m^3 glass
                 conversion                                ! kmol/day -> mol/sec
    else
      cur_waste_form%fuel_dissolution_rate = 0.d0
    endif
    cur_waste_form => cur_waste_form%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMGlassSolve

! ************************************************************************** !

subroutine PMGlassFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  class(pm_glass_type) :: this
  
end subroutine PMGlassFinalizeTimestep

! ************************************************************************** !

subroutine PMGlassUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  class(pm_glass_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update glass mass here

end subroutine PMGlassUpdateSolution  

! ************************************************************************** !

subroutine PMGlassUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
  
  class(pm_glass_type) :: this

  this%option%io_buffer = 'PMGlassUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMGlassUpdateAuxVars   

! ************************************************************************** !

subroutine PMGlassCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_glass_type) :: this
  PetscViewer :: viewer
  
end subroutine PMGlassCheckpoint

! ************************************************************************** !

subroutine PMGlassRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_glass_type) :: this
  PetscViewer :: viewer
  
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMGlassRestart

! ************************************************************************** !

recursive subroutine PMGlassFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/20/15
  
  implicit none
  
  class(pm_glass_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMGlassFinalizeRun

! ************************************************************************** !

subroutine PMGlassStrip(this)
  ! 
  ! Destroys strips Glass process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  use Utility_module, only : DeallocateArray
  use Dataset_module

  implicit none
  
  class(pm_glass_type) :: this
  type(glass_type), pointer :: cur_waste_form, prev_waste_form
  
  PetscInt :: i
  
  call PMWasteFormStrip(this)
  call DatasetDestroy(this%mass_fraction_dataset)

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => cur_waste_form%next
    deallocate(prev_waste_form)
    nullify(prev_waste_form)
  enddo
  nullify(this%waste_form_list)
  
end subroutine PMGlassStrip
  
! ************************************************************************** !

subroutine PMGlassDestroy(this)
  ! 
  ! Destroys Glass process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_glass_type) :: this
  
  call PMGlassStrip(this)
  
end subroutine PMGlassDestroy
  
! ************************************************************************** !

subroutine FMDMInit(fmdm)
  ! 
  ! Initializes the fuel matrix degradation model waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/05/2015

  implicit none
  
  type(fmdm_type) :: fmdm

  fmdm%local_id = UNINITIALIZED_INTEGER
  fmdm%coordinate%x = UNINITIALIZED_DOUBLE
  fmdm%coordinate%y = UNINITIALIZED_DOUBLE
  fmdm%coordinate%z = UNINITIALIZED_DOUBLE
  fmdm%temperature = UNINITIALIZED_DOUBLE
  fmdm%fuel_dissolution_rate = UNINITIALIZED_DOUBLE
  fmdm%specific_surface_area = UNINITIALIZED_DOUBLE
  fmdm%volume = UNINITIALIZED_DOUBLE
  fmdm%burnup = UNINITIALIZED_DOUBLE
  nullify(fmdm%concentration)
  nullify(fmdm%next)
  
end subroutine FMDMInit

! ************************************************************************** !

function PMFMDMCreate()
  ! 
  ! Creates the FMDM waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type), pointer :: PMFMDMCreate
  
  allocate(PMFMDMCreate)
  call PMWasteFormInit(PMFMDMCreate)
  nullify(PMFMDMCreate%waste_form_list)
  PMFMDMCreate%num_grid_cells_in_waste_form = UNINITIALIZED_INTEGER
  PMFMDMCreate%num_concentrations = 11
  PMFMDMCreate%iUO2_2p = 1
  PMFMDMCreate%iUCO3_2n = 2
  PMFMDMCreate%iUO2 = 3
  PMFMDMCreate%iCO3_2n = 4
  PMFMDMCreate%iO2 = 5
  PMFMDMCreate%iH2O2 = 6
  PMFMDMCreate%iFe_2p = 7
  PMFMDMCreate%iH2 = 8
  PMFMDMCreate%iUO2_sld = 9
  PMFMDMCreate%iUO3_sld = 10
  PMFMDMCreate%iUO4_sld = 11
  allocate(PMFMDMCreate%mapping_fmdm_to_pflotran( &
             PMFMDMCreate%num_concentrations))
  PMFMDMCreate%mapping_fmdm_to_pflotran = UNINITIALIZED_INTEGER
  allocate(PMFMDMCreate%mapping_fmdm(4))
  PMFMDMCreate%mapping_fmdm = [PMFMDMCreate%iO2, &
                               PMFMDMCreate%iCO3_2n, &
                               PMFMDMCreate%iH2, &
                               PMFMDMCreate%iFe_2p]
  PMFMDMCreate%initialized = PETSC_FALSE

end function PMFMDMCreate

! ************************************************************************** !

subroutine PMFMDMRead(this,input)
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
  
  class(pm_fmdm_type) :: this
  type(input_type) :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  type(fmdm_type), pointer :: new_waste_form, prev_waste_form

  option => this%option
  
  error_string = 'FMDM'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    call PMWasteFormReadSelectCase(this,input,word,found,error_string,option)
    
    if (found) cycle
    
    select case(trim(word))
    
      case('NUM_GRID_CELLS')
        call InputReadInt(input,option,this%num_grid_cells_in_waste_form)
        call InputErrorMsg(input,option,'num_grid_cells',error_string)
        if (this%num_grid_cells_in_waste_form /= 40) then
          option%io_buffer = 'The FMDM model is currently hardwired to 40 ' // &
            'grid cells.  Please set NUM_GRID_CELLs to 40. This will be ' // &
            'fixed in the future.'
          call printErrMsg(option)
        endif
      case('WASTE_FORM')
        error_string = 'FMDM,WASTE_FORM'
        allocate(new_waste_form)
        call FMDMInit(new_waste_form)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
            case('SPECIFIC_SURFACE_AREA')
              call InputReadDouble(input,option, &
                                   new_waste_form%specific_surface_area)
              call InputErrorMsg(input,option,'specific surface area', &
                                 error_string)
            case('VOLUME')
              call InputReadDouble(input,option,new_waste_form%volume)
              call InputErrorMsg(input,option,'volume',error_string)
            case('BURNUP')
              call InputReadDouble(input,option,new_waste_form%burnup)
              call InputErrorMsg(input,option,'burnup',error_string)
            case('COORDINATE')
              call GeometryReadCoordinate(input,option, &
                                          new_waste_form%coordinate, &
                                          error_string)
          end select
        enddo
        if (Uninitialized(new_waste_form%specific_surface_area)) then
          option%io_buffer = &
            'SPECIFIC_SURFACE_AREA must be specified for all fuel matrix ' // &
            'degradation model waste packages.'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_waste_form%volume)) then
          option%io_buffer = &
            'VOLUME must be specified for all fuel matrix ' // &
            'degradation model waste packages.'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_waste_form%burnup)) then
          option%io_buffer = &
            'BURNUP must be specified for all fuel matrix ' // &
            'degradation model waste packages.'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_waste_form%coordinate%z)) then
          option%io_buffer = &
            'COORDINATE must be specified for all fuel matrix ' // &
            'degradation model waste packages.'
          call printErrMsg(option)
        endif
        if (.not.associated(this%waste_form_list)) then
          this%waste_form_list => new_waste_form
          prev_waste_form => new_waste_form
        else
          prev_waste_form%next => new_waste_form
          prev_waste_form => new_waste_form
        endif
        nullify(new_waste_form)
        error_string = 'FMDM'
      case('BYPASS_WARNING_MESSAGE')
        bypass_warning_message = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
  if (Uninitialized(this%num_grid_cells_in_waste_form)) then
    option%io_buffer = &
      'NUM_GRID_CELLS must be specified for fuel matrix degradation model.'
    call printErrMsg(option)
  endif
  if (len_trim(this%data_mediator_species) == 0) then
    option%io_buffer = &
      'DATA_MEDIATOR_SPECIES must be specified for fuel matrix ' // &
      'degradation model.'
    call printErrMsg(option)
  endif
  
end subroutine PMFMDMRead

! ************************************************************************** !

subroutine PMFMDMSetup(this)
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
  
  class(pm_fmdm_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  type(fmdm_type), pointer :: cur_waste_form, prev_waste_form, next_waste_form
  PetscInt :: i, j, k, local_id
  PetscErrorCode :: ierr
  
  grid => this%realization%patch%grid
  option => this%realization%option
  reaction => this%realization%reaction
  
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    local_id = -1
    select case(grid%itype)
      case(STRUCTURED_GRID)
        call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                            cur_waste_form%coordinate%x, &
                                            cur_waste_form%coordinate%y, &
                                            cur_waste_form%coordinate%z, &
                                            i,j,k)
        if (i > 0 .and. j > 0 .and. k > 0) then
          local_id = i + (j-1)*grid%structured_grid%nlx + &
                      (k-1)*grid%structured_grid%nlxy
        endif
      case(IMPLICIT_UNSTRUCTURED_GRID)
        call UGridGetCellFromPoint(cur_waste_form%coordinate%x, &
                                   cur_waste_form%coordinate%y, &
                                   cur_waste_form%coordinate%z, &
                                   grid%unstructured_grid,option,local_id)
      case default
          option%io_buffer = 'Only STRUCTURED_GRID and ' // &
            'IMPLICIT_UNSTRUCTURED_GRID types supported in PMFMDM.'
          call printErrMsg(option)
    end select
    if (local_id > 0) then
      cur_waste_form%local_id = local_id
      ! allocate concentration array
      allocate(cur_waste_form%concentration(this%num_concentrations, &
                                            this%num_grid_cells_in_waste_form))
      cur_waste_form%concentration = 1.d-20
      prev_waste_form => cur_waste_form
      cur_waste_form => cur_waste_form%next
    else
      ! remove waste form
      next_waste_form => cur_waste_form%next
      if (associated(prev_waste_form)) then
        prev_waste_form%next => next_waste_form
      else
        this%waste_form_list => next_waste_form
      endif
      deallocate(cur_waste_form)
      cur_waste_form => next_waste_form
    endif
  enddo
  
  ! set up indexing of solute concentrations
  species_name = 'O2(aq)'
  this%mapping_fmdm_to_pflotran(this%iO2) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'HCO3-'
  this%mapping_fmdm_to_pflotran(this%iCO3_2n) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'H2(aq)'
  this%mapping_fmdm_to_pflotran(this%iH2) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  species_name = 'Fe++'
  this%mapping_fmdm_to_pflotran(this%iFe_2p) = &
    GetPrimarySpeciesIDFromName(species_name,reaction,option)
  
end subroutine PMFMDMSetup

! ************************************************************************** !

recursive subroutine PMFMDMInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15 
  use Communicator_Base_module
  use Reaction_Aux_module
  use Realization_Base_class
  use Data_Mediator_Vec_class
  
  implicit none

#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pm_fmdm_type) :: this
  
  IS :: is
  type(fmdm_type), pointer :: cur_waste_form
  PetscInt :: num_waste_form_cells
  PetscInt :: i
  PetscInt :: data_mediator_species_id
  PetscInt, allocatable :: waste_form_cell_ids(:)
  PetscReal :: time
  PetscErrorCode :: ierr
  
#ifdef PM_WP_DEBUG  
  call printMsg(this%option,'PMRT%InitializeRun()')
#endif

#ifndef FMDM_MODEL
  this%option%io_buffer = 'Preprocessing statement FMDM_MODEL must be ' // &
    'defined and the ANL FMDM library must be linked to PFLOTRAN to ' // &
    'employ the fuel matrix degradation model.'
  if (.not.bypass_warning_message) then
    call printErrMsg(this%option)
  endif
#endif

  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif

  data_mediator_species_id = &
    GetPrimarySpeciesIDFromName(this%data_mediator_species, &
                                this%realization%reaction,this%option)
  
  ! set up mass transfer
  call RealizCreateTranMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%tran_data_mediator_list)
  ! create a Vec sized by # waste packages * # primary dofs influenced by 
  ! waste package
  ! count of waste form cells
  cur_waste_form => this%waste_form_list
  num_waste_form_cells = 0
  do
    if (.not.associated(cur_waste_form)) exit
    num_waste_form_cells = num_waste_form_cells + 1
    cur_waste_form => cur_waste_form%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,num_waste_form_cells, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(waste_form_cell_ids(num_waste_form_cells))
    waste_form_cell_ids = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      i = i + 1
      waste_form_cell_ids(i) = cur_waste_form%local_id
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    waste_form_cell_ids(:) = waste_form_cell_ids(:) - 1
    waste_form_cell_ids(:) = waste_form_cell_ids(:) + &
                             this%realization%patch%grid%global_offset
    waste_form_cell_ids(:) = waste_form_cell_ids(:) * this%option%ntrandof
    waste_form_cell_ids(:) = waste_form_cell_ids(:)  + &
                             data_mediator_species_id - 1
  endif
  call ISCreateGeneral(this%option%mycomm,num_waste_form_cells, &
                       waste_form_cell_ids,PETSC_COPY_VALUES,is, &
                       ierr);CHKERRQ(ierr)
  if (allocated(waste_form_cell_ids)) deallocate(waste_form_cell_ids)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  
  time = 0.d0
  call PMFMDMSolve(this,time,ierr)  

end subroutine PMFMDMInitializeRun

! ************************************************************************** !

subroutine PMFMDMInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Global_module
  
  implicit none
  
  class(pm_fmdm_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," FUEL MATRIX DEGRADATION MODEL ",47("="))')
  endif

end subroutine PMFMDMInitializeTimestep

! ************************************************************************** !

subroutine PMFMDMPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Grid_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_fmdm_type) :: this
  
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(fmdm_type), pointer :: cur_waste_form
  PetscInt :: i
  PetscInt :: icomp_fmdm
  PetscInt :: icomp_pflotran
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
  cur_waste_form => this%waste_form_list
  do 
    if (.not.associated(cur_waste_form)) exit
    local_id = cur_waste_form%local_id
    ghosted_id = grid%nL2G(local_id)
    cur_waste_form%temperature = global_auxvars(ghosted_id)%temp
    ! overwrite the components in this%mapping_pflotran array
    do i = 1, size(this%mapping_fmdm)
      icomp_fmdm = this%mapping_fmdm(i)
      icomp_pflotran = this%mapping_fmdm_to_pflotran(icomp_fmdm)
      cur_waste_form%concentration(icomp_fmdm,1) = &
        ! the 1 in the second index if for the liquid phase
        rt_auxvars(ghosted_id)%total(icomp_pflotran,1)
    enddo
    cur_waste_form => cur_waste_form%next
  enddo
  
end subroutine PMFMDMPreSolve

! ************************************************************************** !

subroutine PMFMDMSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  !
!  use Argonne_Mixed_Potential_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  interface
    subroutine AMP_step ( burnup, sTme, temperature_C, conc, initialRun, &
                          fuelDisRate, success )
      real ( kind = 8), intent( in )  :: burnup   
      real ( kind = 8), intent( in )  :: sTme   
      real ( kind = 8), intent( in )  :: temperature_C   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      real ( kind = 8), intent(out) :: fuelDisRate
      integer ( kind = 4), intent(out) :: success
    end subroutine
  end interface  

  class(pm_fmdm_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(fmdm_type), pointer :: cur_waste_form
  PetscInt :: i
  PetscReal, pointer :: vec_p(:)       ! g(U)/m^2/yr -> mol(U)/m^2/sec
  PetscReal, parameter :: conversion = 1.d0/238.d0/(365.d0*24.d0*3600.d0)
  
  integer ( kind = 4) :: success
  logical ( kind = 4) :: initialRun
  
  if (this%initialized) then
    initialRun = PETSC_FALSE
  else
    initialRun = PETSC_TRUE
    this%initialized = PETSC_TRUE
  endif

  ierr = 0
  call PMFMDMPreSolve(this)
  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  cur_waste_form => this%waste_form_list
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    i = i + 1
#ifdef FMDM_MODEL  
    call AMP_step(cur_waste_form%burnup, time, &
                  cur_waste_form%temperature, &
                  cur_waste_form%concentration, initialRun, &
                  cur_waste_form%fuel_dissolution_rate, success)
#else
    success = 1
    cur_waste_form%fuel_dissolution_rate = cur_waste_form%burnup
#endif
    if (success == 0) then
      ierr = 1
      exit
    endif      ! mol(U)/sec
    vec_p(i) = cur_waste_form%fuel_dissolution_rate * & ! g/m^2/yr
               cur_waste_form%specific_surface_area * & ! m^2/m^3 waste
               cur_waste_form%volume * &                ! m^3 waste
               conversion                               ! g(U)/yr -> mol(U)/sec
    cur_waste_form => cur_waste_form%next
  enddo
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMFMDMSolve

! ************************************************************************** !

subroutine PMFMDMPostSolve(this)
  ! 
  ! PMFMDMUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  implicit none
  
  class(pm_fmdm_type) :: this
  
end subroutine PMFMDMPostSolve

! ************************************************************************** !

function PMFMDMAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
  PetscBool :: PMFMDMAcceptSolution
  
  ! do nothing
  PMFMDMAcceptSolution = PETSC_TRUE
  
end function PMFMDMAcceptSolution

! ************************************************************************** !

subroutine PMFMDMUpdatePropertiesTS(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
!  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMFMDMUpdatePropertiesTS

! ************************************************************************** !

subroutine PMFMDMTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15 
  
  implicit none
  
  class(pm_fmdm_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMFMDMTimeCut

! ************************************************************************** !

subroutine PMFMDMFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
end subroutine PMFMDMFinalizeTimestep

! ************************************************************************** !

subroutine PMFMDMUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMFMDMUpdateSolution  

! ************************************************************************** !

subroutine PMFMDMUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this

  this%option%io_buffer = 'PMFMDMUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMFMDMUpdateAuxVars   

! ************************************************************************** !

subroutine PMFMDMCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_fmdm_type) :: this
  PetscViewer :: viewer
  
end subroutine PMFMDMCheckpoint

! ************************************************************************** !

subroutine PMFMDMRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
#include "finclude/petscviewer.h"      

  class(pm_fmdm_type) :: this
  PetscViewer :: viewer
  
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMFMDMRestart

! ************************************************************************** !

recursive subroutine PMFMDMFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  
  implicit none
  
  class(pm_fmdm_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMFMDMFinalizeRun

! ************************************************************************** !

subroutine FMDMDestroy()
  ! 
  ! Initializes the fuel matrix degradation model waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/05/2015
  use Utility_module, only : DeallocateArray

  implicit none
  
  type(fmdm_type) :: fmdm

  call DeallocateArray(fmdm%concentration)
  
end subroutine FMDMDestroy

! ************************************************************************** !

subroutine PMFMDMStrip(this)
  ! 
  ! Destroys strips FMDM process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_fmdm_type) :: this
  type(fmdm_type), pointer :: cur_waste_form, prev_waste_form
  
  PetscInt :: i
  
  call PMWasteFormStrip(this)
  
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => cur_waste_form%next
    call DeallocateArray(prev_waste_form%concentration)
    deallocate(prev_waste_form)
    nullify(prev_waste_form)
  enddo
  call DeallocateArray(this%mapping_fmdm_to_pflotran)
  call DeallocateArray(this%mapping_fmdm)
  ! this is solely a pointer
!  call DataMediatorVecDestroy(this%data_mediator)
  nullify(this%waste_form_list)
  
end subroutine PMFMDMStrip
  
! ************************************************************************** !

subroutine PMFMDMDestroy(this)
  ! 
  ! Destroys FMDM process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
  call PMFMDMStrip(this)
  
end subroutine PMFMDMDestroy
  
end module PM_Waste_Form_class
