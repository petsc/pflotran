module PM_Waste_Form_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  use Dataset_Base_class
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscBool, public :: bypass_warning_message = PETSC_FALSE

  type, public :: wf_species_type
   PetscReal, allocatable :: formula_weight(:)
   PetscInt, allocatable :: column_id(:)
   PetscInt, allocatable :: ispecies(:)
   PetscInt :: num_species
   character(len=MAXWORDLENGTH), allocatable :: name(:)
  end type wf_species_type

  type, public :: fmdm_species_type
   PetscReal, allocatable :: formula_weight(:)
   PetscInt, allocatable :: column_id(:)
   PetscInt :: num_species
   character(len=MAXWORDLENGTH), allocatable :: name(:)
  end type fmdm_species_type

  type :: waste_form_base_type
    PetscInt :: id
    PetscInt :: local_cell_id
    type(point3d_type) :: coordinate
    PetscReal :: volume
    PetscReal, pointer :: instantaneous_mass_rate(:)    ! mol/sec
    PetscReal, pointer :: cumulative_mass(:)            ! mol
    PetscBool :: canister_degradation_flag
    PetscReal :: canister_vitality
    PetscReal :: canister_vitality_rate
    class(waste_form_base_type), pointer :: next
  end type waste_form_base_type
  
  type, extends(waste_form_base_type) :: waste_form_fmdm_type
    PetscReal, pointer :: concentration(:,:)
    PetscReal :: specific_surface_area
    PetscReal :: burnup
  end type waste_form_fmdm_type
  
  type, extends(waste_form_base_type) :: waste_form_glass_type
    PetscReal :: exposure_factor
    PetscReal :: glass_dissolution_rate  ! kg / sec
  end type waste_form_glass_type
  
  type, public, extends(pm_base_type) :: pm_waste_form_type
    class(realization_subsurface_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: data_mediator_species
    class(data_mediator_vec_type), pointer :: data_mediator
    class(waste_form_base_type), pointer :: waste_form_list
    type(wf_species_type) :: wf_species
    class(dataset_base_type), pointer ::  mass_fraction_dataset
    PetscBool :: print_mass_balance
    PetscBool :: canister_degradation_model
    PetscReal :: vitality_rate_mean
    PetscReal :: vitality_rate_stdev
    PetscReal :: vitality_rate_trunc
  contains
    procedure, public :: PMWasteFormSetRealization
  end type pm_waste_form_type
  
  type, public, extends(pm_waste_form_type) :: pm_waste_form_glass_type
    PetscReal :: specific_surface_area
    PetscReal :: glass_density
  contains
    procedure, public :: Setup => PMGlassSetup
    procedure, public :: Read => PMGlassRead
    procedure, public :: InitializeRun => PMWFGlassInitializeRun
    procedure, public :: InitializeTimestep => PMWFGlassInitializeTimestep
    procedure, public :: FinalizeTimestep => PMGlassFinalizeTimestep
    procedure, public :: UpdateSolution => PMGlassUpdateSolution
    procedure, public :: Solve => PMGlassSolve
    procedure, public :: Checkpoint => PMGlassCheckpoint    
    procedure, public :: Restart => PMGlassRestart  
    procedure, public :: Destroy => PMGlassDestroy
  end type pm_waste_form_glass_type
  
  type, public, extends(pm_waste_form_type) :: pm_waste_form_fmdm_type
    type(fmdm_species_type) :: fmdm_species
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
    procedure, public :: CheckpointBinary => PMFMDMCheckpointBinary
    procedure, public :: RestartBinary => PMFMDMRestartBinary
    procedure, public :: Destroy => PMFMDMDestroy
  end type pm_waste_form_fmdm_type
  
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
  nullify(this%mass_fraction_dataset)
  this%wf_species%num_species = 0
 !------- canister degradation model --------------
  this%canister_degradation_model = PETSC_FALSE
  this%vitality_rate_mean = UNINITIALIZED_DOUBLE
  this%vitality_rate_stdev = UNINITIALIZED_DOUBLE
  this%vitality_rate_trunc = UNINITIALIZED_DOUBLE
 !-------------------------------------------------
!geh: only initialize if a pointer, instead of allocatable
!  nullify(this%wf_species%name)
!  nullify(this%wf_species%formula_weight)
!  nullify(this%wf_species%column_id)
!  nullify(this%wf_species%ispecies)

end subroutine PMWasteFormInit

! ************************************************************************** !

subroutine PMWasteFormReadSelectCase(this,input,keyword,found,error_string, &
                                     option)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Updated/modified by Jenn Frederick, 2/3/2016

  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class 
  use String_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, units
  character(len=MAXSTRINGLENGTH) :: temp_buf, string, units_category
  character(len=MAXSTRINGLENGTH) :: species_name_buf
  character(len=MAXSTRINGLENGTH) :: species_formula_wt_buf
  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscInt :: k, icol

  found = PETSC_TRUE
  species_name_buf = ''
  species_formula_wt_buf = ''
  k = 0

  select case(trim(keyword))
!-------------------------------------
    case('SPECIES')
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        k = k + 1
        call InputReadWord(input,option,word,PETSC_TRUE)
        species_name_buf = trim(species_name_buf) // ' ' // trim(word)
        call InputReadWord(input,option,word,PETSC_TRUE)
        species_formula_wt_buf = trim(species_formula_wt_buf) // ' ' &
                                 // trim(word)
        this%wf_species%num_species = k
      enddo
      if (k == 0) then
        option%io_buffer = 'At least one species name and formula weight &
                           &must be given in the ' // trim(error_string) &
                           // ' block.'
        call printErrMsg(option)
      endif
      allocate(this%wf_species%name(k))  
      allocate(this%wf_species%formula_weight(k))
      allocate(this%wf_species%column_id(k))
      allocate(this%wf_species%ispecies(k))
      input%buf = species_name_buf
      k = 0
      do while (k < this%wf_species%num_species)
        k = k + 1
        call InputReadWord(input,option,this%wf_species%name(k),PETSC_TRUE)
        call InputErrorMsg(input,option,'species name',error_string)
      enddo
      input%buf = species_formula_wt_buf
      k = 0
      do while (k < this%wf_species%num_species)
        k = k + 1
        call InputReadDouble(input,option,this%wf_species%formula_weight(k))
        call InputErrorMsg(input,option,'species formula weight',error_string)
        this%wf_species%column_id(k) = UNINITIALIZED_INTEGER
        this%wf_species%ispecies(k) = UNINITIALIZED_INTEGER
      enddo
!-------------------------------------
    case('MASS_FRACTION')
      temp_buf = input%buf
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'mass fraction file/list',error_string)
      call StringToUpper(word)
      select case(trim(word))
        case('FILE') ! OK format, do now throw error
        case('LIST') ! OK format, do not throw error
        case default
          option%io_buffer = 'A mass fraction FILE or LIST must be given in &
                             &the ' // trim(error_string) // ' block. Only &
                             &FILE or LIST supported in the ' &
                             // trim(error_string) // ' block.'
        call printErrMsg(option)
      end select
      input%buf = temp_buf
      dataset_ascii => DatasetAsciiCreate()
      this%mass_fraction_dataset => dataset_ascii
      dataset_ascii%data_type = DATASET_REAL
      units_category = 'unitless'
      call ConditionReadValues(input,option,word,this%mass_fraction_dataset, &
                               units,units_category)
      if (associated(dataset_ascii%time_storage)) then
        ! default time interpolation is linear
        if (dataset_ascii%time_storage%time_interpolation_method == &
             INTERPOLATION_NULL) then
          dataset_ascii%time_storage%time_interpolation_method = &
               INTERPOLATION_LINEAR
        endif
      endif
      if (dataset_ascii%header == '') then
        option%io_buffer = 'A HEADER must be specified in the ' // &
                           trim(error_string) // ' mass fraction file/list.'
        call printErrMsg(option)
      endif
!-------------------------------------
    case('CANISTER_DEGRADATION_MODEL')
      this%canister_degradation_model = PETSC_TRUE
      option%io_buffer = 'canister degradation model'
      call printMsg(option)
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call StringToUpper(word)
        select case(trim(word))
          case('VITALITY_LOG10_MEAN')
            call InputReadDouble(input,option,this%vitality_rate_mean)
            call InputErrorMsg(input,option,'canister vitality log-10 &
                               &mean value',error_string)
          case('VITALITY_LOG10_STDEV')
            call InputReadDouble(input,option,this%vitality_rate_stdev)
            call InputErrorMsg(input,option,'canister vitality log-10 &
                               &st. dev. value',error_string)
          case('VITALITY_UPPER_TRUNCATION')
            call InputReadDouble(input,option,this%vitality_rate_trunc)
            call InputErrorMsg(input,option,'canister vitality log-10 &
                               &upper truncation value',error_string)
          case default
            option%io_buffer = 'Keyword ' // trim(word) // ' not recognized &
                               &in the ' // trim(error_string) // &
                               ' CANISTER_DEGRADATION_MODEL block.'
            call printErrMsg(option)
        end select
      enddo
!-------------------------------------
    case('PRINT_MASS_BALANCE')
      this%print_mass_balance = PETSC_TRUE
!-------------------------------------    
    case default
      found = PETSC_FALSE
!-------------------------------------
  end select

end subroutine PMWasteFormReadSelectCase

! ************************************************************************** !

subroutine PMWFReadError(this,input,option,error_string)
  ! 
  ! Checks for input deck reading errors for the waste form process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/12/2016

  use Option_module
  use Input_Aux_module

  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  if (.not.associated(this%mass_fraction_dataset)) then
    option%io_buffer = 'MASS_FRACTION must be specified in the ' // &
                       trim(error_string) // ' block.'
    call printErrMsg(option)
  endif
  if (.not.associated(this%waste_form_list)) then
    option%io_buffer = 'At least one WASTE_FORM must be specified in the ' // &
                       trim(error_string) // ' block.'
    call printErrMsg(option)
  endif
  if (.not.allocated(this%wf_species%name)) then
    option%io_buffer = 'At least one SPECIES NAME and FORMULA_WEIGHT must be &
                       &specified in the ' // trim(error_string) // ' block.'
    call printErrMsg(option)
  endif
  
  if (this%canister_degradation_model) then
    if (uninitialized(this%vitality_rate_mean)) then
      option%io_buffer = 'VITALITY_LOG10_MEAN must be given in the '&
                         // trim(error_string) // &
                         ', CANISTER_DEGRADATION_MODEL block.'
      call printErrMsg(option)
    endif
    if (uninitialized(this%vitality_rate_stdev)) then
      option%io_buffer = 'VITALITY_LOG10_STDEV must be given in the '&
                         // trim(error_string) // &
                         ', CANISTER_DEGRADATION_MODEL block.'
      call printErrMsg(option)
    endif
    if (uninitialized(this%vitality_rate_trunc)) then
      option%io_buffer = 'VITALITY_UPPER_TRUNCATION must be given in the '&
                         // trim(error_string) // &
                         ', CANISTER_DEGRADATION_MODEL block.'
      call printErrMsg(option)
    endif
  endif

end subroutine PMWFReadError

! ************************************************************************** !

subroutine PMWFAssignColIdsFromHeader(this,input,option,error_string)
  ! 
  ! Reads the mass fraction file header and assigns a column id to each species
  ! in the waste form process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/12/2016

  use Option_module
  use String_module, only : StringToUpper
  use Input_Aux_module

  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscInt :: icol, k
  character(len=MAXWORDLENGTH) :: word

  input%buf = this%mass_fraction_dataset%header
  input%ierr = 0
  icol = 0
  call InputReadWord(input,option,word,PETSC_TRUE)
  call StringToUpper(word)
  if ((input%ierr == 0) .and. (trim(word) == 'HEADER')) then
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'while reading the mass fraction file/&
                       &list header',error_string)
    call StringToUpper(word)
    if (trim(word) == 'TIME') then
      do
        k = 0
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          icol = icol + 1
              
          do while (k < this%wf_species%num_species)
            k = k + 1
            if (trim(word) == trim(this%wf_species%name(k))) then
              this%wf_species%column_id(k) = icol
              exit
            endif
          enddo ! k loop

        else
          exit
        endif
      enddo ! icol loop
      k = 0
      do while (k < this%wf_species%num_species)
        k = k + 1
        if (Uninitialized(this%wf_species%column_id(k))) then
          option%io_buffer = 'Mismatch between species in the ' &
                              // trim(error_string) // ' mass fraction file/&
                              &list header and those listed in the ' &
                              // trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
      enddo
    else
      option%io_buffer = 'The first column in the ' // trim(error_string) &
                         // ' mass fraction file/list must be TIME.'
      call printErrMsg(option)
    endif
  else
    option%io_buffer = 'A HEADER must be specified in the ' // &
                       trim(error_string) // ' mass fraction file/list.'
    call printErrMsg(option)
  endif

end subroutine PMWFAssignColIdsFromHeader

! ************************************************************************** !

subroutine PMWasteFormSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Realization_Subsurface_class

  implicit none
  
  class(pm_waste_form_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWasteFormSetRealization

! ************************************************************************** !

 subroutine PMWFInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/25/15
  use Reaction_Aux_module
  use Realization_Base_class
  use Utility_module, only : GetRndNumFromNormalDist
  
  implicit none

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_type) :: this
  
  IS :: is
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: num_waste_form_cells
  PetscInt :: i, j
  PetscInt :: data_mediator_species_id
  PetscInt, allocatable :: species_indices_in_residual(:)
  PetscErrorCode :: ierr

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    allocate(cur_waste_form%instantaneous_mass_rate&
                                  (this%wf_species%num_species))
    allocate(cur_waste_form%cumulative_mass&
                                  (this%wf_species%num_species))
    cur_waste_form%instantaneous_mass_rate = UNINITIALIZED_DOUBLE
    cur_waste_form%cumulative_mass = UNINITIALIZED_DOUBLE
    if (this%canister_degradation_model) then
      cur_waste_form%canister_degradation_flag = PETSC_TRUE
      cur_waste_form%canister_vitality = 1.d0
      call GetRndNumFromNormalDist(this%vitality_rate_mean, &
                                   this%vitality_rate_stdev, &
                                   cur_waste_form%canister_vitality_rate)
      if (cur_waste_form%canister_vitality_rate > this%vitality_rate_trunc) then
        cur_waste_form%canister_vitality_rate = this%vitality_rate_trunc
      endif
      ! Given rates are in units of log-10/yr, so convert to 1/yr:
      cur_waste_form%canister_vitality_rate = &
                                10.0**(cur_waste_form%canister_vitality_rate)
    endif
    cur_waste_form => cur_waste_form%next
  enddo
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif
  
  if (.not.this%option%restart_flag .and. this%print_mass_balance) then
    call PMWFOutputHeader(this)
    call PMWFOutput(this)
  endif

  do j = 1, this%wf_species%num_species
    this%wf_species%ispecies(j) = &
      GetPrimarySpeciesIDFromName(this%wf_species%name(j), &
                                  this%realization%reaction,this%option)
  enddo
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
  call VecCreateSeq(PETSC_COMM_SELF,num_waste_form_cells* &
                                    this%wf_species%num_species, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(species_indices_in_residual(num_waste_form_cells* &
                                         this%wf_species%num_species))
    species_indices_in_residual = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      do j = 1, this%wf_species%num_species
        i = i + 1
        species_indices_in_residual(i) = &
          (cur_waste_form%local_cell_id-1)*this%option%ntrandof + &
          this%wf_species%ispecies(j)
      enddo
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    species_indices_in_residual(:) = species_indices_in_residual(:) - 1
    ! set to global petsc index
    species_indices_in_residual(:) = species_indices_in_residual(:) + &
      this%realization%patch%grid%global_offset*this%option%ntrandof
  endif
  call ISCreateGeneral(this%option%mycomm,num_waste_form_cells* &
                                          this%wf_species%num_species, &
                       species_indices_in_residual,PETSC_COPY_VALUES,is, &
                       ierr);CHKERRQ(ierr)
  if (allocated(species_indices_in_residual)) &
    deallocate(species_indices_in_residual)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  
end subroutine PMWFInitializeRun

! ************************************************************************** !

subroutine PMWasteFormStrip(this)
  ! 
  ! Destroys a waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this

  nullify(this%realization)
  nullify(this%data_mediator)
  deallocate(this%wf_species%name)  
  deallocate(this%wf_species%formula_weight)
  deallocate(this%wf_species%column_id)
  deallocate(this%wf_species%ispecies)
  
  
end subroutine PMWasteFormStrip

! ************************************************************************** !

subroutine WasteFormBaseInit(base)
  ! 
  ! Initializes the base waste form data
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(waste_form_base_type) :: base

  base%id = UNINITIALIZED_INTEGER
  base%local_cell_id = UNINITIALIZED_INTEGER
  base%coordinate%x = UNINITIALIZED_DOUBLE
  base%coordinate%y = UNINITIALIZED_DOUBLE
  base%coordinate%z = UNINITIALIZED_DOUBLE
  nullify(base%instantaneous_mass_rate)
  nullify(base%cumulative_mass)
  base%volume = UNINITIALIZED_DOUBLE
 !------- canister degradation model -----------------
  base%canister_degradation_flag = PETSC_FALSE
  base%canister_vitality = 0d0
  base%canister_vitality_rate = UNINITIALIZED_DOUBLE
 !----------------------------------------------------

end subroutine WasteFormBaseInit

! ************************************************************************** !

subroutine PMWFBaseSetup(this)
  ! 
  ! Maps waste forms to grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Option_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  class(waste_form_base_type), pointer :: cur_waste_form, prev_waste_form, &
                                          next_waste_form
  PetscInt :: i, j, k, local_id
  PetscInt :: waste_form_id
  PetscErrorCode :: ierr
  
  grid => this%realization%patch%grid
  option => this%realization%option
  
  waste_form_id = 0
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    waste_form_id = waste_form_id + 1
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
      cur_waste_form%id = waste_form_id
      cur_waste_form%local_cell_id = local_id
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
  
end subroutine PMWFBaseSetup

! ************************************************************************** !

subroutine PMWFOutput(this)
  ! 
  ! Maps waste forms to grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module
  use Output_Aux_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  
  if (.not.associated(this%waste_form_list)) return
  
100 format(100es16.8)

  option => this%realization%option
  output_option => this%realization%output_option
  
  fid = 86
  filename = PMWFOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    do i = 1, this%wf_species%num_species
      write(fid,100,advance="no") cur_waste_form%cumulative_mass(i), &
                                  cur_waste_form%instantaneous_mass_rate(i) * &
                                  output_option%tconv
    enddo
    select type(cur_waste_form)
      class is (waste_form_glass_type)
        write(fid,100,advance="no") cur_waste_form%volume, &
                                    cur_waste_form%glass_dissolution_rate * &
                                    output_option%tconv, &
                                    cur_waste_form%canister_vitality*100.0
      class is (waste_form_fmdm_type)
        write(fid,100,advance="no") cur_waste_form%canister_vitality*100.0
    end select
    cur_waste_form => cur_waste_form%next
  enddo
  close(fid)
  
end subroutine PMWFOutput

! ************************************************************************** !

function PMWFOutputFilename(option)
  ! 
  ! Generates filename for waste form output
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: PMWFOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMWFOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-wf_mass-' // trim(adjustl(word)) // '.dat'
  
end function PMWFOutputFilename  

! ************************************************************************** !

subroutine PMWFOutputHeader(this)
  ! 
  ! Writes header for waste form output file
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Output_Aux_module
  use Grid_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(output_option_type), pointer :: output_option
  type(grid_type), pointer :: grid
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  
  fid = 86
  filename = PMWFOutputFilename(this%option)
  open(unit=fid,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    ! cell natural id
    write(cell_string,*) grid%nG2A(grid%nL2G(cur_waste_form%local_cell_id))
    cell_string = ' (' // trim(adjustl(cell_string)) // ')'
    ! coordinate of waste form
    x_string = BestFloat(cur_waste_form%coordinate%x,1.d4,1.d-2)
    y_string = BestFloat(cur_waste_form%coordinate%y,1.d4,1.d-2)
    z_string = BestFloat(cur_waste_form%coordinate%z,1.d4,1.d-2)
    cell_string = trim(cell_string) // &
             ' (' // trim(adjustl(x_string)) // &
             ' ' // trim(adjustl(y_string)) // &
             ' ' // trim(adjustl(z_string)) // ')'
    do i = 1, this%wf_species%num_species
      variable_string = trim(this%wf_species%name(i)) // &
        ' Mass Flux'
      ! cumulative
      units_string = 'mol'
      call OutputWriteToHeader(fid,variable_string,units_string, &
                               cell_string,icolumn)
      ! instantaneous
      units_string = 'mol/' // trim(adjustl(output_option%tunit))
      call OutputWriteToHeader(fid,variable_string,units_string, &
                               cell_string,icolumn)
    enddo
    select type(cur_waste_form)
      class is (waste_form_glass_type)
        variable_string = 'WF Volume'
        units_string = 'm^3'
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'WF Dissolution Rate'
        units_string = 'kg/' // trim(adjustl(output_option%tunit))
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = 'WF Canister Vitality'
        units_string = '%' 
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
      class is (waste_form_fmdm_type)
        variable_string = 'WF Canister Vitality'
        units_string = '%' 
        call OutputWriteToHeader(fid,variable_string,units_string, &
                                 cell_string,icolumn)
    end select
    cur_waste_form => cur_waste_form%next
  enddo
  
  close(fid)
  
end subroutine PMWFOutputHeader

! ************************************************************************** !

subroutine PMWFInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  use Dataset_module
  
  implicit none
  
  class(pm_waste_form_type) :: this

  class(waste_form_base_type), pointer :: cur_waste_form
  PetscReal :: dt

  dt = this%option%tran_dt
  ! update cumulative values
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    cur_waste_form%cumulative_mass = cur_waste_form%cumulative_mass + &
                                   cur_waste_form%instantaneous_mass_rate * dt
    cur_waste_form => cur_waste_form%next
  enddo  

  if (this%print_mass_balance) then
    call PMWFOutput(this)
  endif

end subroutine PMWFInitializeTimestep

! ************************************************************************** !

subroutine WFGlassInit(glass)
  ! 
  ! Initializes the glass waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  implicit none
  
  class(waste_form_glass_type) :: glass

  call WasteFormBaseInit(glass)
  glass%exposure_factor = UNINITIALIZED_DOUBLE
  glass%glass_dissolution_rate = UNINITIALIZED_DOUBLE
  nullify(glass%next)

end subroutine WFGlassInit

! ************************************************************************** !

function WFGlassCast(this)
  ! 
  ! Casts waste_form_base_type to waste_form_glass_type
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(waste_form_base_type), pointer :: this
  
  class(waste_form_glass_type), pointer :: WFGlassCast
  
  nullify(WFGlassCast)
  if (associated(this)) then
    select type(this)
      class is(waste_form_glass_type)
        WFGlassCast => this
      class default
        print *, 'Wrong class in WFGlassCast'
        stop
    end select
  endif

end function WFGlassCast

! ************************************************************************** !

function WFFMDMCast(this)
  ! 
  ! Casts waste_form_base_type to waste_form_glass_type
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(waste_form_base_type), pointer :: this
  
  class(waste_form_fmdm_type), pointer :: WFFMDMCast
  
  nullify(WFFMDMCast)
  if (associated(this)) then
    select type(this)
      class is(waste_form_fmdm_type)
        WFFMDMCast => this
      class default
        print *, 'Wrong class in WFGlassCast'
        stop
    end select
  endif

end function WFFMDMCast

! ************************************************************************** !

function PMGlassCreate()
  ! 
  ! Creates the Glass waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type), pointer :: PMGlassCreate
  
  allocate(PMGlassCreate)
  call PMWasteFormInit(PMGlassCreate)
  nullify(PMGlassCreate%waste_form_list)
  PMGlassCreate%specific_surface_area = UNINITIALIZED_DOUBLE
  PMGlassCreate%glass_density = 2.65d3 ! kg/m^3

end function PMGlassCreate

! ************************************************************************** !

subroutine PMGlassRead(this,input)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Input_Aux_module
  use Utility_module
  use Option_module
  use String_module
  use Units_module
  
  implicit none
  
  class(pm_waste_form_glass_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  class(waste_form_glass_type), pointer :: new_waste_form, prev_waste_form
  PetscBool :: found

  option => this%option
  error_string = 'GLASS'
  input%ierr = 0

  option%io_buffer = 'pflotran card:: ' // trim(error_string)
  call printMsg(option)

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
    !-------------------------------------
      case('WASTE_FORM')
        error_string = 'GLASS,WASTE_FORM'
        allocate(new_waste_form)
        call WFGlassInit(new_waste_form)
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
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                new_waste_form%volume = UnitsConvertToInternal(word,'volume', &
                                        option) * new_waste_form%volume
              endif
            case('COORDINATE')
              call GeometryReadCoordinate(input,option, &
                                          new_waste_form%coordinate, &
                                          error_string)
            case default
              call InputKeywordUnrecognized(word,error_string,option)
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
    !-------------------------------------
      case('SPECIFIC_SURFACE_AREA')
        call InputReadDouble(input,option,this%specific_surface_area)
        call InputErrorMsg(input,option,'specific surface area',error_string)
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          this%specific_surface_area = UnitsConvertToInternal(word, &
                            'area/volume',option) * this%specific_surface_area
        endif
    !-------------------------------------
      case('GLASS_DENSITY')
        call InputReadDouble(input,option,this%glass_density)
        call InputErrorMsg(input,option,'glass density',error_string)
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-------------------------------------
    end select
  enddo

  call PMWFReadError(this,input,option,error_string)
  call PMWFAssignColIdsFromHeader(this,input,option,error_string)

  if (Uninitialized(this%specific_surface_area)) then
    option%io_buffer = 'SPECIFIC_SURFACE_AREA must be specified in ' // &
      trim(error_string)
    call printErrMsg(option)
  endif
    
end subroutine PMGlassRead

! ************************************************************************** !

subroutine PMGlassSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this

  call PMWFBaseSetup(this)
  
end subroutine PMGlassSetup

! ************************************************************************** !

recursive subroutine PMWFGlassInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Dataset_module
  use Time_Storage_module
  
  implicit none

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_glass_type) :: this
  
  type(time_storage_type), pointer :: null_time_storage
  PetscErrorCode :: ierr
  
  call PMWFInitializeRun(this)
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif

  nullify(null_time_storage)
  call DatasetVerify(this%mass_fraction_dataset,null_time_storage,this%option)

  call PMGlassSolve(this,0.d0,ierr)
  
end subroutine PMWFGlassInitializeRun

! ************************************************************************** !

subroutine PMWFGlassInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  use Dataset_module
  
  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
  class(waste_form_glass_type), pointer :: cur_waste_form
  PetscReal :: rate
  PetscReal :: dV
  PetscReal, parameter :: conversion = 1.d0/(24.d0*3600.d0)
  
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GLASS MODEL ",65("="))')
  endif
  
  ! due to output witin, must be called prior to update
  call PMWFInitializeTimestep(this)

  ! update mass balances after transport step
  cur_waste_form => WFGlassCast(this%waste_form_list)
  do 
    if (.not.associated(cur_waste_form)) exit
    ! m^3 glass
    dV = cur_waste_form%glass_dissolution_rate / & ! kg glass/sec
         this%glass_density * &                    ! kg glass/m^3 glass
         this%option%tran_dt                       ! sec
    cur_waste_form%volume = cur_waste_form%volume - dV
    cur_waste_form => WFGlassCast(cur_waste_form%next)
  enddo
  
  call DatasetUpdate(this%mass_fraction_dataset,this%option%time,this%option)

end subroutine PMWFGlassInitializeTimestep

! ************************************************************************** !

subroutine PMGlassPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
end subroutine PMGlassPreSolve

! ************************************************************************** !

subroutine PMGlassSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  use Grid_module
  use Global_Aux_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_glass_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  class(waste_form_glass_type), pointer :: cur_waste_form
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: i, j
  PetscReal, pointer :: vec_p(:)            ! 1/day -> 1/sec
  PetscReal, parameter :: time_conversion = 1.d0/(24.d0*3600.d0)
  PetscReal :: fuel_dissolution_rate

  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  cur_waste_form => WFGlassCast(this%waste_form_list)
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    if (cur_waste_form%canister_degradation_flag) then
!     ---------------- Vitality degradation function --------------------------
      cur_waste_form%canister_vitality = cur_waste_form%canister_vitality &
                        - ( cur_waste_form%canister_vitality_rate * & ! [1/yr]
                            this%option%tran_dt * &                   ! [sec]
                            (1.0/(365.0*24.0*3600.0)) )               ! [yr/sec]
      if (cur_waste_form%canister_vitality < 1.d-3) then
        cur_waste_form%canister_vitality = 0.d0
      endif
!     -------------------------------------------------------------------------
    endif
    if ((cur_waste_form%volume > 0.d0) .and. &
        (cur_waste_form%canister_vitality == 0.d0)) then
      fuel_dissolution_rate = & ! kg glass/m^2/day
        560.d0*exp(-7397.d0/ &
            (global_auxvars(grid%nL2G(cur_waste_form%local_cell_id))%temp+273.15d0))
      ! kg glass / sec
      cur_waste_form%glass_dissolution_rate = &
        fuel_dissolution_rate * &          ! kg glass (dissolving)/m^2/day
        this%specific_surface_area * &     ! m^2/kg glass
        cur_waste_form%exposure_factor * & ! [-]
        cur_waste_form%volume * &          ! m^3 glass
        this%glass_density * &             ! kg glass/m^3 glass
        time_conversion                    ! 1/day -> 1/sec
      ! mol/sec
      do j = 1, this%wf_species%num_species
        i = i + 1
        cur_waste_form%instantaneous_mass_rate(j) = &
          cur_waste_form%glass_dissolution_rate * & ! kg glass / sec
          this%wf_species%formula_weight(j) * &  ! kmol radnuclide/kg radnuclide
          this%mass_fraction_dataset% &
          rarray(this%wf_species%column_id(j)) * &  ! kg radionuclide/kg glass
          1.d3                                      ! kmol -> mol
        vec_p(i) = cur_waste_form%instantaneous_mass_rate(j)    ! mol/sec
      enddo
    else
      cur_waste_form%glass_dissolution_rate = 0.d0
      cur_waste_form%instantaneous_mass_rate = 0.d0
    endif
    cur_waste_form => WFGlassCast(cur_waste_form%next)
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMGlassSolve

! ************************************************************************** !

subroutine PMGlassFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
end subroutine PMGlassFinalizeTimestep

! ************************************************************************** !

subroutine PMGlassUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update glass mass here

end subroutine PMGlassUpdateSolution  

! ************************************************************************** !

subroutine PMGlassUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this

  this%option%io_buffer = 'PMGlassUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMGlassUpdateAuxVars   

! ************************************************************************** !

subroutine PMGlassCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  use Option_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_glass_type) :: this
  PetscViewer :: viewer
  
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: maximum_waste_form_id
  PetscInt :: local_waste_form_count
  PetscInt :: temp_int
  
  Vec :: local, global
  PetscErrorCode :: ierr
  
  this%option%io_buffer = 'PMGlassCheckpoint not implemented.'
  call printErrMsg(this%option)
  
  ! calculate maximum waste form id
  maximum_waste_form_id = 0
  local_waste_form_count = 0
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    local_waste_form_count = local_waste_form_count + 1
    maximum_waste_form_id = max(maximum_waste_form_id,cur_waste_form%id)
    cur_waste_form => cur_waste_form%next
  enddo
  call MPI_Allreduce(maximum_waste_form_id,temp_int,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,this%option%mycomm,ierr)
!  call VecCreateMPI(this%option%mycomm,local_waste_form_count,PETSC_DETERMINE,ierr)
  
                     
end subroutine PMGlassCheckpoint

! ************************************************************************** !

subroutine PMGlassRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_glass_type) :: this
  PetscViewer :: viewer
  
  this%option%io_buffer = 'PMGlassRestart not implemented.'
  call printErrMsg(this%option)
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
  ! Date: 08/26/15
  
  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
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
  ! Date: 08/26/15

  use Dataset_module

  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
  class(waste_form_glass_type), pointer :: cur_waste_form, prev_waste_form
  
  PetscInt :: i
  
  call PMWasteFormStrip(this)
  call DatasetDestroy(this%mass_fraction_dataset)

  cur_waste_form => WFGlassCast(this%waste_form_list)
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => WFGlassCast(cur_waste_form%next)
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
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_glass_type) :: this
  
  call PMGlassStrip(this)
  
end subroutine PMGlassDestroy
  
! ************************************************************************** !

subroutine WFFMDMInit(fmdm)
  ! 
  ! Initializes the fuel matrix degradation model waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(waste_form_fmdm_type) :: fmdm

  call WasteFormBaseInit(fmdm)  
  fmdm%specific_surface_area = UNINITIALIZED_DOUBLE
  fmdm%burnup = UNINITIALIZED_DOUBLE
  nullify(fmdm%concentration)
  nullify(fmdm%next)
  
end subroutine WFFMDMInit

! ************************************************************************** !

function PMFMDMCreate()
  ! 
  ! Creates the FMDM waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type), pointer :: PMFMDMCreate
  
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
  ! Date: 08/26/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module
  
  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  class(waste_form_fmdm_type), pointer :: new_waste_form, prev_waste_form
  PetscBool :: found

  option => this%option
  error_string = 'FMDM'
  input%ierr = 0

  option%io_buffer = 'pflotran card:: ' // trim(error_string)
  call printMsg(option)

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
    !-------------------------------------
      case('NUM_GRID_CELLS')
        call InputReadInt(input,option,this%num_grid_cells_in_waste_form)
        call InputErrorMsg(input,option,'num_grid_cells',error_string)
        if (this%num_grid_cells_in_waste_form /= 40) then
          option%io_buffer = 'The FMDM model is currently hardwired to 40 ' // &
            'grid cells.  Please set NUM_GRID_CELLs to 40. This will be ' // &
            'fixed in the future.'
          call printErrMsg(option)
        endif
      !-------------------------------------
      case('WASTE_FORM')
        error_string = 'FMDM,WASTE_FORM'
        allocate(new_waste_form)
        call WFFMDMInit(new_waste_form)
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
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                new_waste_form%specific_surface_area = &
                           UnitsConvertToInternal(word,'area/volume',option) * &
                           new_waste_form%specific_surface_area
              endif
            case('VOLUME')
              call InputReadDouble(input,option,new_waste_form%volume)
              call InputErrorMsg(input,option,'volume',error_string)
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                new_waste_form%volume = UnitsConvertToInternal(word,'volume', &
                                        option) * new_waste_form%volume
              endif
            case('BURNUP')
              call InputReadDouble(input,option,new_waste_form%burnup)
              call InputErrorMsg(input,option,'burnup',error_string)
            case('COORDINATE')
              call GeometryReadCoordinate(input,option, &
                                          new_waste_form%coordinate, &
                                          error_string)
            case default
              call InputKeywordUnrecognized(word,error_string,option)
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
    !-------------------------------------
      case('BYPASS_WARNING_MESSAGE')
        bypass_warning_message = PETSC_TRUE
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-------------------------------------
    end select
  enddo

  call PMWFReadError(this,input,option,error_string)
  call PMWFAssignColIdsFromHeader(this,input,option,error_string)
  
  if (Uninitialized(this%num_grid_cells_in_waste_form)) then
    option%io_buffer = &
      'NUM_GRID_CELLS must be specified for fuel matrix degradation model.'
    call printErrMsg(option)
  endif

end subroutine PMFMDMRead

! ************************************************************************** !

subroutine PMFMDMSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Reaction_Aux_module
  use Option_module

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  class(waste_form_fmdm_type), pointer :: cur_waste_form
  
  option => this%realization%option
  reaction => this%realization%reaction
  
  call PMWFBaseSetup(this)

  cur_waste_form => WFFMDMCast(this%waste_form_list)
  do 
    if (.not.associated(cur_waste_form)) exit
    ! allocate concentration array
    allocate(cur_waste_form%concentration(this%num_concentrations, &
                                          this%num_grid_cells_in_waste_form))
    cur_waste_form%concentration = 1.d-20
    allocate(cur_waste_form%instantaneous_mass_rate(1))
    cur_waste_form%instantaneous_mass_rate = 0.d0
    allocate(cur_waste_form%cumulative_mass(1))
    cur_waste_form%cumulative_mass = 0.d0
    cur_waste_form => WFFMDMCast(cur_waste_form%next)
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
  ! Initializes the process model for the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  
  implicit none

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_fmdm_type) :: this
  
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
  
  call PMWFInitializeRun(this)

  time = 0.d0
  call PMFMDMSolve(this,time,ierr)  

end subroutine PMFMDMInitializeRun

! ************************************************************************** !

subroutine PMFMDMInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," FUEL MATRIX DEGRADATION MODEL ",47("="))')
  endif

  ! due to output witin, must be called prior to update
  call PMWFInitializeTimestep(this)

end subroutine PMFMDMInitializeTimestep

! ************************************************************************** !

subroutine PMFMDMPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Grid_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  class(waste_form_fmdm_type), pointer :: cur_waste_form
  PetscInt :: i
  PetscInt :: icomp_fmdm
  PetscInt :: icomp_pflotran
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
  cur_waste_form => WFFMDMCast(this%waste_form_list)
  do 
    if (.not.associated(cur_waste_form)) exit
    ghosted_id = grid%nL2G(cur_waste_form%local_cell_id)
    ! overwrite the components in this%mapping_pflotran array
    do i = 1, size(this%mapping_fmdm)
      icomp_fmdm = this%mapping_fmdm(i)
      icomp_pflotran = this%mapping_fmdm_to_pflotran(icomp_fmdm)
      cur_waste_form%concentration(icomp_fmdm,1) = &
        ! the 1 in the second index if for the liquid phase
        rt_auxvars(ghosted_id)%total(icomp_pflotran,1)
    enddo
    cur_waste_form => WFFMDMCast(cur_waste_form%next)
  enddo
  
end subroutine PMFMDMPreSolve

! ************************************************************************** !

subroutine PMFMDMSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
!  use Argonne_Mixed_Potential_module
  use Grid_module
  use Global_Aux_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  interface
    subroutine AMP_step ( burnup, sTme, temperature_C, conc, initialRun, &
                          fuelDisRate, success )
      real ( kind = 8), intent( in ) :: burnup   
      real ( kind = 8), intent( in ) :: sTme   
      real ( kind = 8), intent( in ) :: temperature_C   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      real ( kind = 8), intent(out) :: fuelDisRate
      integer ( kind = 4), intent(out) :: success
    end subroutine
  end interface  

  class(pm_waste_form_fmdm_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  class(waste_form_fmdm_type), pointer :: cur_waste_form
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(grid_type), pointer :: grid
   PetscInt :: i
  PetscReal, pointer :: vec_p(:)       ! g(U)/m^2/yr -> mol(U)/m^2/sec
  PetscReal, parameter :: conversion = 1.d0/238.d0/(365.d0*24.d0*3600.d0)
  PetscReal :: fuel_dissolution_rate   ! g/m^2/yr
  
  integer ( kind = 4) :: success
  logical ( kind = 4) :: initialRun
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
 
  if (this%initialized) then
    initialRun = PETSC_FALSE
  else
    initialRun = PETSC_TRUE
    this%initialized = PETSC_TRUE
  endif

  ierr = 0
  call PMFMDMPreSolve(this)
  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  cur_waste_form => WFFMDMCast(this%waste_form_list)
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    i = i + 1
#ifdef FMDM_MODEL  
    call AMP_step(cur_waste_form%burnup, time, &
                  global_auxvars(grid%nL2G(cur_waste_form%local_cell_id))%temp, &
                  cur_waste_form%concentration, initialRun, &
                  fuel_dissolution_rate, success)
#else
    success = 1
    fuel_dissolution_rate = cur_waste_form%burnup
#endif
    if (success == 0) then
      ierr = 1
      exit
    endif      ! mol(U)/sec
    cur_waste_form%instantaneous_mass_rate = &
               fuel_dissolution_rate * &                ! g/m^2/yr
               cur_waste_form%specific_surface_area * & ! m^2/m^3 waste
               cur_waste_form%volume * &                ! m^3 waste
               conversion                               ! g(U)/yr -> mol(U)/sec
    vec_p(i) = cur_waste_form%instantaneous_mass_rate(1)
    cur_waste_form => WFFMDMCast(cur_waste_form%next)
  enddo
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMFMDMSolve

! ************************************************************************** !

subroutine PMFMDMPostSolve(this)
  ! 
  ! PMFMDMUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! 
  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
end subroutine PMFMDMPostSolve

! ************************************************************************** !

function PMFMDMAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
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
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
!  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMFMDMUpdatePropertiesTS

! ************************************************************************** !

subroutine PMFMDMTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  
  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMFMDMTimeCut

! ************************************************************************** !

subroutine PMFMDMFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
end subroutine PMFMDMFinalizeTimestep

! ************************************************************************** !

subroutine PMFMDMUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMFMDMUpdateSolution  

! ************************************************************************** !

subroutine PMFMDMUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this

  this%option%io_buffer = 'PMFMDMUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMFMDMUpdateAuxVars   

! ************************************************************************** !

subroutine PMFMDMCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_fmdm_type) :: this
  PetscViewer :: viewer
  
end subroutine PMFMDMCheckpointBinary

! ************************************************************************** !

subroutine PMFMDMRestartBinary(this,viewer)
  !
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_fmdm_type) :: this
  PetscViewer :: viewer
  
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMFMDMRestartBinary

! ************************************************************************** !

recursive subroutine PMFMDMFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  
  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMFMDMFinalizeRun

! ************************************************************************** !

subroutine WFFMDMStrip(fmdm)
  ! 
  ! Initializes the fuel matrix degradation model waste form
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(waste_form_fmdm_type) :: fmdm

  call DeallocateArray(fmdm%concentration)
  
end subroutine WFFMDMStrip

! ************************************************************************** !

subroutine PMFMDMStrip(this)
  ! 
  ! Destroys strips FMDM process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  class(waste_form_fmdm_type), pointer :: cur_waste_form, prev_waste_form
  
  PetscInt :: i
  
  call PMWasteFormStrip(this)
  
  cur_waste_form => WFFMDMCast(this%waste_form_list)
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => WFFMDMCast(cur_waste_form%next)
    call WFFMDMStrip(prev_waste_form)
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
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_fmdm_type) :: this
  
  call PMFMDMStrip(this)
  
end subroutine PMFMDMDestroy
  
end module PM_Waste_Form_class
