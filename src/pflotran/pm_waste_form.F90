module PM_Waste_Form_class

  use PM_Base_class
  use Realization_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type :: fmdm_type
    PetscInt :: local_id
    type(point3d_type) :: coordinate
    PetscReal :: temperature
    PetscReal, pointer :: concentration(:,:)
    PetscReal :: fuel_dissolution_rate
    PetscReal :: specific_surface_area
    PetscReal :: burnup
    type(fmdm_type), pointer :: next
  end type fmdm_type
  
  type, public, extends(pm_base_type) :: pm_fmdm_type
    class(realization_type), pointer :: realization
    type(fmdm_type), pointer :: waste_form_list
    type(point3d_type), pointer :: coordinate(:)
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
    character(len=MAXWORDLENGTH) :: data_mediator_species
    class(data_mediator_vec_type), pointer :: data_mediator
    PetscBool :: initialized
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Init => PMWasteFormInit
    procedure, public :: Read => PMWasteFormRead
!    procedure, public :: SetupSolvers => PMWasteFormSetupSolvers
    procedure, public :: PMWasteFormSetRealization
    procedure, public :: InitializeRun => PMWasteFormInitializeRun
!!    procedure, public :: FinalizeRun => PMWasteFormFinalizeRun
    procedure, public :: InitializeTimestep => PMWasteFormInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWasteFormFinalizeTimestep
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
  end type pm_fmdm_type
  
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
  fmdm%burnup = UNINITIALIZED_DOUBLE
  nullify(fmdm%concentration)
  nullify(fmdm%next)
  
end subroutine FMDMInit

! ************************************************************************** !

function PMWasteFormCreate()
  ! 
  ! Creates the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type), pointer :: PMWasteFormCreate
  
  allocate(PMWasteFormCreate)
  nullify(PMWasteFormCreate%realization)
  nullify(PMWasteFormCreate%waste_form_list)
  nullify(PMWasteFormCreate%coordinate)
  PMWasteFormCreate%num_grid_cells_in_waste_form = UNINITIALIZED_INTEGER
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
  nullify(PMWasteFormCreate%data_mediator)
  PMWasteFormCreate%data_mediator_species = ''
  allocate(PMWasteFormCreate%mapping_fmdm_to_pflotran( &
             PMWasteFormCreate%num_concentrations))
  PMWasteFormCreate%mapping_fmdm_to_pflotran = UNINITIALIZED_INTEGER
  allocate(PMWasteFormCreate%mapping_fmdm(4))
  PMWasteFormCreate%mapping_fmdm = [PMWasteFormCreate%iO2, &
                                   PMWasteFormCreate%iCO3_2n, &
                                   PMWasteFormCreate%iH2, &
                                   PMWasteFormCreate%iFe_2p]
  PMWasteFormCreate%initialized = PETSC_FALSE
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
  
  class(pm_fmdm_type) :: this
  type(input_type) :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
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
    
    select case(trim(word))
    
      case('NUM_GRID_CELLS')
        call InputReadInt(input,option,this%num_grid_cells_in_waste_form)
        call InputErrorMsg(input,option,'num_grid_cells',error_string)
      case('DATA_MEDIATOR_SPECIES')
        call InputReadWord(input,option,this%data_mediator_species,PETSC_TRUE)
        call InputErrorMsg(input,option,'data_mediator_species',error_string)
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
            'IMPLICIT_UNSTRUCTURED_GRID types supported in PMWasteForm.'
          call printErrMsg(option)
    end select
    if (local_id > 0) then
      cur_waste_form%local_id = local_id
      ! allocate concentrationa rray
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
  
end subroutine PMWasteFormInit

! ************************************************************************** !

subroutine PMWasteFormSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Realization_class

  implicit none
  
  class(pm_fmdm_type) :: this
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
  call printErrMsg(this%option)
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
  call PMWasteFormSolve(this,time,ierr)  

end subroutine PMWasteFormInitializeRun

! ************************************************************************** !

subroutine PMWasteFormInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  use Global_module
  
  implicit none
  
  class(pm_fmdm_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," FUEL MATRIX DEGRADATION MODEL ",47("="))')
  endif

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
  
end subroutine PMWasteFormPreSolve

! ************************************************************************** !

subroutine PMWasteFormSolve(this,time,ierr)
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
                          fuelDisRate, status )
      real ( kind = 8), intent( in )  :: burnup   
      real ( kind = 8), intent( in )  :: sTme   
      real ( kind = 8), intent( in )  :: temperature_C   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      real ( kind = 8), intent(out) :: fuelDisRate
      integer ( kind = 4), intent(out) :: status
    end subroutine
  end interface  

  class(pm_fmdm_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(fmdm_type), pointer :: cur_waste_form
  PetscInt :: i
  PetscReal, pointer :: vec_p(:)       ! g(U)/m^2/yr -> mol(U)/m^2/sec
  PetscReal, parameter :: conversion = 1.d0/238.d0/(365.d0*24.d0*3600.d0)
  
  integer ( kind = 4) :: status
  logical ( kind = 4) :: initialRun
  
  if (this%initialized) then
    initialRun = PETSC_FALSE
  else
    initialRun = PETSC_TRUE
    this%initialized = PETSC_TRUE
  endif

  ierr = 0
  call PMWasteFormPreSolve(this)
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
                  cur_waste_form%fuel_dissolution_rate, status)
#endif
    if (status == 0) then
      ierr = 1
      exit
    endif
    vec_p(i) = cur_waste_form%fuel_dissolution_rate * & ! g/m^2/yr
               cur_waste_form%specific_surface_area * &
               conversion
    cur_waste_form => cur_waste_form%next
  enddo
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWasteFormSolve

! ************************************************************************** !

subroutine PMWasteFormPostSolve(this)
  ! 
  ! PMWasteFormUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  implicit none
  
  class(pm_fmdm_type) :: this
  
end subroutine PMWasteFormPostSolve

! ************************************************************************** !

function PMWasteFormAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
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
  
  class(pm_fmdm_type) :: this
  
!  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMWasteFormUpdatePropertiesTS

! ************************************************************************** !

subroutine PMWasteFormTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15 
  
  implicit none
  
  class(pm_fmdm_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMWasteFormTimeCut

! ************************************************************************** !

subroutine PMWasteFormFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
end subroutine PMWasteFormFinalizeTimestep

! ************************************************************************** !

subroutine PMWasteFormUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMWasteFormUpdateSolution  

! ************************************************************************** !

subroutine PMWasteFormUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15

  implicit none
  
  class(pm_fmdm_type) :: this

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

  class(pm_fmdm_type) :: this
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

  class(pm_fmdm_type) :: this
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
  
  class(pm_fmdm_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWasteFormFinalizeRun

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

subroutine PMWasteFormDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_fmdm_type) :: this
  type(fmdm_type), pointer :: cur_waste_form, prev_waste_form
  
  PetscInt :: i
  
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
  
end subroutine PMWasteFormDestroy
  
end module PM_Waste_Form_class
