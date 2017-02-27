module PM_WIPP_SrcSink_class

  use PM_Base_class
  use Region_module
  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Data_Mediator_Vec_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: chem_species_type
    PetscReal :: in_panel                      ! [mol/m3]
    PetscReal, pointer :: initial_conc_mol(:)  ! [mol/m3]
    PetscReal, pointer :: inst_rate(:)         ! [mol/m3/sec]
    PetscReal, pointer :: current_conc_mol(:)  ! [mol/m3]
    PetscReal, pointer :: current_conc_kg(:)   ! [kg/m3] per BRAGFLO U.M.
  end type chem_species_type
  
  type, public :: inventory_type
    character(len=MAXWORDLENGTH) :: name
    type(chem_species_type) :: Fe_s 
    type(chem_species_type) :: FeOH2_s
    type(chem_species_type) :: C6H10O5_s
    type(chem_species_type) :: H_ion_aq
    type(chem_species_type) :: NO3_minus_aq
    type(chem_species_type) :: CO2_g
    type(chem_species_type) :: N2_g
    type(chem_species_type) :: SO42_minus_aq
    type(chem_species_type) :: H2S_g
    type(chem_species_type) :: FeS_s
    type(chem_species_type) :: MgO_s
    type(chem_species_type) :: MgOH2_s
    type(chem_species_type) :: Mg5CO34OH24H2_s
    type(chem_species_type) :: MgCO3_s
    type(inventory_type), pointer :: next
  end type inventory_type
  
  type, public :: srcsink_panel_type
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    type(inventory_type), pointer :: inventory
    character(len=MAXWORDLENGTH) :: inventory_name
    PetscReal, pointer :: scaling_factor(:)
    PetscReal, pointer :: gas_generation_rate(:)
    PetscReal, pointer :: brine_generation_rate(:)
    PetscInt :: id
    PetscMPIInt :: myMPIgroup_id
    PetscMPIInt :: myMPIcomm
    type(srcsink_panel_type), pointer :: next
  end type srcsink_panel_type

  type, public, extends(pm_base_type) :: pm_wipp_srcsink_type
    PetscReal :: alpharxn
    PetscReal :: smin
    PetscReal :: satwick
    PetscReal :: inundated_corrosion_rate
    PetscReal :: humid_corrosion_factor
    PetscReal :: inundated_biodeg_rate
    PetscReal :: humid_biodeg_factor
    PetscReal :: inundated_brucite_rate
    PetscReal :: humid_brucite_rate
    PetscReal :: RXH2S_factor
    PetscReal :: RXCO2_factor
    PetscReal :: hymagcon_rate
    type(srcsink_panel_type), pointer :: waste_panel_list
    type(inventory_type), pointer :: inventory_list
    class(data_mediator_vec_type), pointer :: data_mediator
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: PMWSSSetRealization
    procedure, public :: Setup => PMWSSSetup
    procedure, public :: Read => PMWSSRead
    procedure, public :: InitializeRun => PMWSSInitializeRun
    procedure, public :: InitializeTimestep => PMWSSInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWSSFinalizeTimestep
    procedure, public :: Solve => PMWSSSolve
    procedure, public :: Output => PMWSSOutput
    procedure, public :: InputRecord => PMWSSInputRecord
    procedure, public :: Destroy => PMWSSDestroy
  end type pm_wipp_srcsink_type

  public :: PMWSSCreate, &
            WastePanelCreate, &
            InventoryCreate

contains

! *************************************************************************** !

function PMWSSCreate()
  !
  ! Creates and initializes the WIPP source/sink process model.
  !
  ! Author: Jenn Frederick
  ! Date: 1/31/2017
  !
  
  implicit none
  
  class(pm_wipp_srcsink_type), pointer :: PMWSSCreate
  
  allocate(PMWSSCreate)
  nullify(PMWSSCreate%waste_panel_list)
  nullify(PMWSSCreate%inventory_list)
  nullify(PMWSSCreate%data_mediator)
  nullify(PMWSSCreate%realization)
  PMWSSCreate%name = 'wipp source sink'
  PMWSSCreate%alpharxn = UNINITIALIZED_DOUBLE
  PMWSSCreate%smin = UNINITIALIZED_DOUBLE
  PMWSSCreate%satwick = UNINITIALIZED_DOUBLE
  PMWSSCreate%inundated_biodeg_rate = UNINITIALIZED_DOUBLE
  PMWSSCreate%inundated_brucite_rate = UNINITIALIZED_DOUBLE
  PMWSSCreate%inundated_corrosion_rate = UNINITIALIZED_DOUBLE
  PMWSSCreate%humid_biodeg_factor = UNINITIALIZED_DOUBLE
  PMWSSCreate%humid_brucite_rate = UNINITIALIZED_DOUBLE
  PMWSSCreate%humid_corrosion_factor = UNINITIALIZED_DOUBLE
  PMWSSCreate%RXH2S_factor = UNINITIALIZED_DOUBLE
  PMWSSCreate%RXCO2_factor = UNINITIALIZED_DOUBLE
  PMWSSCreate%hymagcon_rate = UNINITIALIZED_DOUBLE
  
  call PMBaseInit(PMWSSCreate)
  
end function PMWSSCreate

! *************************************************************************** !

function WastePanelCreate()
  !
  ! Creates and initializes a waste panel type.
  !
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
  type(srcsink_panel_type), pointer :: WastePanelCreate
  
  allocate(WastePanelCreate)
  
  nullify(WastePanelCreate%next)
  nullify(WastePanelCreate%region)
  nullify(WastePanelCreate%inventory)
  nullify(WastePanelCreate%scaling_factor)
  nullify(WastePanelCreate%gas_generation_rate)
  nullify(WastePanelCreate%brine_generation_rate)
  WastePanelCreate%name = ''
  WastePanelCreate%region_name = ''
  WastePanelCreate%inventory_name = ''
  WastePanelCreate%id = 0
  WastePanelCreate%myMPIgroup_id = 0
  WastePanelCreate%myMPIcomm = 0

end function WastePanelCreate

! *************************************************************************** !

function InventoryCreate()
  !
  ! Creates and initializes a waste panel inventory.
  !
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
  type(inventory_type), pointer :: InventoryCreate
  
  allocate(InventoryCreate)
  
  nullify(InventoryCreate%next)
  InventoryCreate%name = ''
  
  call InitChemSpecies(InventoryCreate%Fe_s)  ! iron
  call InitChemSpecies(InventoryCreate%FeOH2_s)
  call InitChemSpecies(InventoryCreate%C6H10O5_s)  ! cellulose
  call InitChemSpecies(InventoryCreate%H_ion_aq)  ! h+
  call InitChemSpecies(InventoryCreate%NO3_minus_aq)  ! nitrate  
  call InitChemSpecies(InventoryCreate%CO2_g)  ! carbon dioxide gas                     
  call InitChemSpecies(InventoryCreate%N2_g)  ! nitrogen gas          
  call InitChemSpecies(InventoryCreate%SO42_minus_aq)  ! sulfate
  call InitChemSpecies(InventoryCreate%H2S_g)  
  call InitChemSpecies(InventoryCreate%FeS_s)                    
  call InitChemSpecies(InventoryCreate%MgO_s)  ! magnesia
  call InitChemSpecies(InventoryCreate%MgOH2_s)                 
  call InitChemSpecies(InventoryCreate%Mg5CO34OH24H2_s)  ! hydromagnesite
  call InitChemSpecies(InventoryCreate%MgCO3_s)

end function InventoryCreate

! *************************************************************************** !

subroutine InitChemSpecies(chem_species)
  !
  ! Initializes a waste panel inventory's chemical species.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  implicit none
  
  type(chem_species_type) :: chem_species
  
  nullify(chem_species%current_conc_mol)        ! [mol/m3]
  nullify(chem_species%current_conc_kg)         ! [kg/m3]
  nullify(chem_species%initial_conc_mol)        ! [mol/m3]
  nullify(chem_species%inst_rate)               ! [mol/m3/sec]
  chem_species%in_panel = UNINITIALIZED_DOUBLE  ! [mol/m3]
  
end subroutine

! *************************************************************************** !

subroutine PMWSSSetRealization(this,realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/02/2017
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWSSSetRealization

! *************************************************************************** !

subroutine PMWSSAssociateRegion(this,region_list)
  ! 
  ! Associates the waste panel to its assigned region via the REGION keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Region_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  PetscBool :: matched
  
  option => this%option
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
      cur_region => region_list%first     
      do
        if (.not.associated(cur_region)) exit
        matched = PETSC_FALSE
        if (StringCompare(trim(cur_region%name), &
                          trim(cur_waste_panel%region_name))) then
          cur_waste_panel%region => cur_region
          matched = PETSC_TRUE
        endif
        if (matched) exit
        cur_region => cur_region%next
      enddo      
      if (.not.associated(cur_waste_panel%region)) then
        option%io_buffer = 'WASTE_PANEL REGION ' // &
                           trim(cur_waste_panel%region_name) // ' not found.'
        call printErrMsg(option)
      endif
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine PMWSSAssociateRegion

! *************************************************************************** !

subroutine PMWSSAssociateInventory(this)
  ! 
  ! Associates the waste panel to its assigned inventory via INVENTORY keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Option_module
  use String_module
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(inventory_type), pointer :: cur_inventory
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  PetscBool :: matched
  
  option => this%option
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
      cur_inventory => this%inventory_list     
      do
        if (.not.associated(cur_inventory)) exit
        matched = PETSC_FALSE
        if (StringCompare(trim(cur_inventory%name), &
                          trim(cur_waste_panel%inventory_name))) then
          cur_waste_panel%inventory => cur_inventory
          matched = PETSC_TRUE
        endif
        if (matched) exit
        cur_inventory => cur_inventory%next
      enddo      
      if (.not.associated(cur_waste_panel%inventory)) then
        option%io_buffer = 'WASTE_PANEL INVENTORY ' // &
                           trim(cur_waste_panel%inventory_name) // ' not found.'
        call printErrMsg(option)
      endif
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine PMWSSAssociateInventory

! *************************************************************************** !

subroutine PMWSSSetRegionScaling(this,waste_panel)
  ! 
  ! Calculates and sets the scaling factor vector for each of the waste panels
  ! that have assigned regions. It assumes the volume of the cells that make up
  ! the region do not change over the course of the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !

  use Material_Aux_class
  use Grid_module

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  type(srcsink_panel_type), pointer :: waste_panel
  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: k, cell_id
  PetscReal :: total_volume_local, total_volume_global
  PetscErrorCode :: ierr
  
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  allocate(waste_panel%scaling_factor(waste_panel%region%num_cells))
  total_volume_local = 0.d0
  total_volume_global = 0.d0
  
  ! scale by cell volume
  do k = 1,waste_panel%region%num_cells
    cell_id = grid%nL2G(waste_panel%region%cell_ids(k))
    waste_panel%scaling_factor(k) = material_auxvars(cell_id)%volume ! [m^3]
    total_volume_local = total_volume_local &
                         + material_auxvars(cell_id)%volume  ! [m^3]
  enddo
  call MPI_Allreduce(total_volume_local,total_volume_global,ONE_INTEGER_MPI, &
              MPI_DOUBLE_PRECISION,MPI_SUM,waste_panel%myMPIcomm,ierr)
  waste_panel%scaling_factor = waste_panel%scaling_factor/total_volume_global 
  
end subroutine PMWSSSetRegionScaling

! *************************************************************************** !

subroutine PMWSSRead(this,input)
  !
  ! Reads input file parameters for the WIPP source/sink process model.
  !
  ! Author: Jenn Frederick
  ! Date: 1/31/2017
  !
  
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: double
  character(len=MAXSTRINGLENGTH) :: error_string
  type(srcsink_panel_type), pointer :: new_waste_panel
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(inventory_type), pointer :: new_inventory
  type(inventory_type), pointer :: cur_inventory
  PetscBool :: added
  
  option => this%option
  input%ierr = 0
  option%io_buffer = 'pflotran card:: WIPP_SOURCE_SINK'
  call printMsg(option)
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    error_string = 'WIPP_SOURCE_SINK'
    call StringToUpper(word)
    select case(trim(word))
    !-----------------------------------------
      case('ALPHARXN')
        call InputReadDouble(input,option,this%alpharxn)
        call InputErrorMsg(input,option,'ALPHARXN',error_string)
      case('SOCMIN')
        call InputReadDouble(input,option,this%smin)
        call InputErrorMsg(input,option,'SOCMIN',error_string)
      case('SATWICK')
        call InputReadDouble(input,option,this%satwick)
        call InputErrorMsg(input,option,'SATWICK',error_string)
      case('INUNDATED_CORROSION_RATE')
        call InputReadDouble(input,option,this%inundated_corrosion_rate)
        call InputErrorMsg(input,option,'INUNDATED_CORROSION_RATE',error_string)
      case('HUMID_CORROSION_FACTOR')
        call InputReadDouble(input,option,this%humid_corrosion_factor)
        call InputErrorMsg(input,option,'HUMID_CORROSION_FACTOR',error_string)
      case('INUNDATED_BIODEGRADATION_RATE')
        call InputReadDouble(input,option,this%inundated_biodeg_rate)
        call InputErrorMsg(input,option,'INUNDATED_BIODEGRADATION_RATE', &
                           error_string)
      case('HUMID_BIODEGRADATION_FACTOR')
        call InputReadDouble(input,option,this%humid_biodeg_factor)
        call InputErrorMsg(input,option,'HUMID_BIODEGRADATION_FACTOR', &
                           error_string)
      case('INUNDATED_BRUCITE_RATE')
        call InputReadDouble(input,option,this%inundated_brucite_rate)
        call InputErrorMsg(input,option,'INUNDATED_BRUCITE_RATE',error_string)
      case('HUMID_BRUCITE_RATE')
        call InputReadDouble(input,option,this%humid_brucite_rate)
        call InputErrorMsg(input,option,'HUMID_BRUCITE_RATE',error_string)
      case('RXH2S_FACTOR')
        call InputReadDouble(input,option,this%RXH2S_factor)
        call InputErrorMsg(input,option,'RXH2S_FACTOR',error_string)
      case('RXCO2_FACTOR')
        call InputReadDouble(input,option,this%RXCO2_factor)
        call InputErrorMsg(input,option,'RXCO2_FACTOR',error_string)
      case('HYDROMAGNESITE_CON_RATE')
        call InputReadDouble(input,option,this%hymagcon_rate)
        call InputErrorMsg(input,option,'HYDROMAGNESITE_CON_RATE',error_string)
    !-----------------------------------------
    !-----------------------------------------
      case('WASTE_PANEL')
        error_string = trim(error_string) // ',WASTE_PANEL'
        allocate(new_waste_panel)
        new_waste_panel => WastePanelCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_waste_panel%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_waste_panel%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('REGION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'region assignment',error_string)
              new_waste_panel%region_name = trim(word)
          !-----------------------------------
            case('INVENTORY')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'inventory assignment',error_string)
              new_waste_panel%inventory_name = trim(word)
          !-----------------------------------    
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages ---------------------
        if (new_waste_panel%region_name == '') then
          option%io_buffer = 'REGION must be specified in the ' // &
                             trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        if (new_waste_panel%inventory_name == '') then
          option%io_buffer = 'INVENTORY must be specified in the ' // &
                             trim(error_string) // ' block.'
          call printErrMsg(option)
        endif
        added = PETSC_FALSE
        if (.not.associated(this%waste_panel_list)) then
          this%waste_panel_list => new_waste_panel
        else
          cur_waste_panel => this%waste_panel_list
          do
            if (.not.associated(cur_waste_panel)) exit
            if (.not.associated(cur_waste_panel%next)) then
              cur_waste_panel%next => new_waste_panel
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_waste_panel => cur_waste_panel%next
          enddo
        endif
        nullify(new_waste_panel)
    !-----------------------------------------
    !-----------------------------------------
      case('INVENTORY')
        error_string = trim(error_string) // ',INVENTORY'
        allocate(new_inventory)
        new_inventory => InventoryCreate()
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        new_inventory%name = adjustl(trim(word))
        error_string = trim(error_string) // ' ' // trim(new_inventory%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
          !-----------------------------------
            case('SOLIDS','SOLID')
              error_string = trim(error_string) // ',SOLIDS'
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !-----------------------------
                  case('FE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial Fe',error_string)
                    new_inventory%Fe_s%in_panel = double
                !-----------------------------
                  case('MGO')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial MgO',error_string)
                    new_inventory%MgO_s%in_panel = double
                !-----------------------------
                  case('CELLULOSE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial cellulose', &
                                       error_string)
                    new_inventory%C6H10O5_s%in_panel = double
                !-----------------------------
                  case default
                    call InputKeywordUnrecognized(word,error_string,option)
                !-----------------------------
                end select
              enddo
          !-----------------------------------
            case('AQUEOUS')
              error_string = trim(error_string) // ',AQUEOUS'
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword',error_string)
                call StringToUpper(word)
                select case(trim(word))
                !-----------------------------
                  case('H+')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial H+',error_string)
                    new_inventory%H_ion_aq%in_panel = double
                !-----------------------------
                  case('NITRATE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial nitrate', &
                                       error_string)
                    new_inventory%NO3_minus_aq%in_panel = double
                !-----------------------------
                  case('SULFATE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial sulfate', &
                                       error_string)
                    new_inventory%SO42_minus_aq%in_panel = double
                !-----------------------------
                  case default
                    call InputKeywordUnrecognized(word,error_string,option)
                !-----------------------------
                end select
              enddo
          !-----------------------------------    
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          !----------------------------------- 
          end select
        enddo
        ! error messages ---------------------
        if (uninitialized(new_inventory%Fe_s%in_panel)) then
          option%io_buffer = 'Initial Fe (solid) inventory must be specified &
                         &using the SOLIDS,Fe keyord in the WIPP_SOURCE_SINK &
                         &block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_inventory%MgO_s%in_panel)) then
          option%io_buffer = 'Initial MgO (solid) inventory must be specified &
                        &using the SOLIDS,MgO keyord in the WIPP_SOURCE_SINK &
                        &block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_inventory%C6H10O5_s%in_panel)) then
          option%io_buffer = 'Initial cellulose (solid) inventory must be &
                        &specified using the SOLIDS,CELLULOSE keyord in the &
                        &WIPP_SOURCE_SINK block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_inventory%H_ion_aq%in_panel)) then
          option%io_buffer = 'Initial H+ (aqueous) inventory must be specified &
                        &using the AQUEOUS,H+ keyord in the WIPP_SOURCE_SINK &
                        &block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_inventory%NO3_minus_aq%in_panel)) then
          option%io_buffer = 'Initial nitrate (aqueous) inventory must be &
                        &specified using the AQUEOUS,NITRATE keyord in the &
                        &WIPP_SOURCE_SINK block.'
          call printErrMsg(option)
        endif
        if (uninitialized(new_inventory%SO42_minus_aq%in_panel)) then
          option%io_buffer = 'Initial sulfate (aqueous) inventory must be &
                        &specified using the AQUEOUS,SULFATE keyord in the &
                        &WIPP_SOURCE_SINK block.'
          call printErrMsg(option)
        endif
        added = PETSC_FALSE
        if (.not.associated(this%inventory_list)) then
          this%inventory_list => new_inventory
        else
          cur_inventory => this%inventory_list
          do
            if (.not.associated(cur_inventory)) exit
            if (.not.associated(cur_inventory%next)) then
              cur_inventory%next => new_inventory
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_inventory => cur_inventory%next
          enddo
        endif
        nullify(new_inventory)
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(word,'WIPP_SOURCE_SINK',option)
    !-----------------------------------------
    end select  
  enddo
  
  if (.not.associated(this%waste_panel_list)) then
    option%io_buffer = 'At least one WASTE_PANEL must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (.not.associated(this%inventory_list)) then
    option%io_buffer = 'At least one INVENTORY must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%alpharxn)) then
    option%io_buffer = 'ALPHARXN must be specified in the WIPP_SOURCE_SINK &
                       &block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%smin)) then
    option%io_buffer = 'SOCMIN must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%satwick)) then
    option%io_buffer = 'SATWICK must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%inundated_biodeg_rate)) then
    option%io_buffer = 'INUNDATED_BIODEGRADATION_RATE must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%inundated_brucite_rate)) then
    option%io_buffer = 'INUNDATED_BRUCITE_RATE must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%inundated_corrosion_rate)) then
    option%io_buffer = 'INUNDATED_CORROSION_RATE must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%humid_biodeg_factor)) then
    option%io_buffer = 'HUMID_BIODEGRADATION_FACTOR must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%humid_brucite_rate)) then
    option%io_buffer = 'HUMID_BRUCITE_RATE must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%humid_corrosion_factor)) then
    option%io_buffer = 'HUMID_CORROSION_FACTOR must be specified in &
                       &the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%RXH2S_factor)) then
    option%io_buffer = 'RXH2S_FACTOR must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%RXCO2_factor)) then
    option%io_buffer = 'RXCO2_FACTOR must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (uninitialized(this%hymagcon_rate)) then
    option%io_buffer = 'HYDROMAGNESITE_CON_RATE must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  
  call PMWSSAssociateInventory(this)
  
end subroutine PMWSSRead

! *************************************************************************** !

subroutine PMWSSSetup(this)
  !
  ! Associates the waste panels to their regions and sets the waste panel id.
  ! Creates an MPI group/communicator for processes that own a waste panel.
  ! Throws out waste panels on processes that do not own the waste panel region.
  !
  ! Author: Jenn Frederick
  ! Date: 2/06/2017
  !
  
  use Option_module

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(option_type), pointer :: option
  type(srcsink_panel_type), pointer :: cur_waste_panel, prev_waste_panel
  type(srcsink_panel_type), pointer :: next_waste_panel
  PetscInt :: waste_panel_id
  PetscBool :: local
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm
  
  option => this%realization%option
  
  ! point the waste panel region to the desired region 
  call PMWSSAssociateRegion(this,this%realization%patch%region_list)
  
  waste_panel_id = 0
  nullify(prev_waste_panel)
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    waste_panel_id = waste_panel_id + 1
    local = PETSC_FALSE
    if (associated(cur_waste_panel%region)) then
      if (cur_waste_panel%region%num_cells > 0) then
          local = PETSC_TRUE
      endif
    endif
    if (local) then
      ! assign the waste panel ID number and MPI communicator
      cur_waste_panel%id = waste_panel_id
      cur_waste_panel%myMPIgroup_id = waste_panel_id
      call MPI_Comm_split(option%mycomm,cur_waste_panel%myMPIgroup_id, &
                          option%myrank,newcomm,ierr)
    else
      cur_waste_panel%id = 0
      cur_waste_panel%myMPIgroup_id = 0
      call MPI_Comm_split(option%mycomm,MPI_UNDEFINED,option%myrank, &
                          newcomm,ierr)
    endif
    cur_waste_panel%myMPIcomm = newcomm
    if (local) then
      call PMWSSSetRegionScaling(this,cur_waste_panel)
      call InventoryAllocate(cur_waste_panel%inventory, &
                             cur_waste_panel%region%num_cells)
      allocate(cur_waste_panel%gas_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%gas_generation_rate(:) = 0.d0
      allocate(cur_waste_panel%brine_generation_rate( &
               cur_waste_panel%region%num_cells))
      cur_waste_panel%brine_generation_rate(:) = 0.d0
      prev_waste_panel => cur_waste_panel
      cur_waste_panel => cur_waste_panel%next
    else 
      ! remove waste panel because it is not local
      next_waste_panel => cur_waste_panel%next
      if (associated(prev_waste_panel)) then
        prev_waste_panel%next => next_waste_panel
      else
        this%waste_panel_list => next_waste_panel
      endif
      deallocate(cur_waste_panel)
      cur_waste_panel => next_waste_panel
    endif
  enddo
  
end subroutine PMWSSSetup

! ************************************************************************** !

subroutine InventoryAllocate(inventory,num_cells)
  ! 
  ! Allocates the size of the chemical species arrays within the inventory to 
  ! the number of cells in the waste panel region.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/23/2017
  !
  
  implicit none
  
  type(inventory_type) :: inventory
  PetscInt :: num_cells
  
  call ChemSpeciesAllocate(num_cells,inventory%Fe_s)
  call ChemSpeciesAllocate(num_cells,inventory%FeOH2_s)
  call ChemSpeciesAllocate(num_cells,inventory%C6H10O5_s)
  call ChemSpeciesAllocate(num_cells,inventory%H_ion_aq)
  call ChemSpeciesAllocate(num_cells,inventory%NO3_minus_aq)
  call ChemSpeciesAllocate(num_cells,inventory%CO2_g)
  call ChemSpeciesAllocate(num_cells,inventory%N2_g)
  call ChemSpeciesAllocate(num_cells,inventory%SO42_minus_aq)
  call ChemSpeciesAllocate(num_cells,inventory%H2S_g)
  call ChemSpeciesAllocate(num_cells,inventory%FeS_s)
  call ChemSpeciesAllocate(num_cells,inventory%MgO_s)
  call ChemSpeciesAllocate(num_cells,inventory%MgOH2_s)
  call ChemSpeciesAllocate(num_cells,inventory%Mg5CO34OH24H2_s)
  call ChemSpeciesAllocate(num_cells,inventory%MgCO3_s)
  
end subroutine InventoryAllocate

! ************************************************************************** !

subroutine ChemSpeciesAllocate(num_cells,chem_species)
  ! 
  ! Allocates the size of the chemical species arrays within the inventory to 
  ! the number of cells in the waste panel region.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/23/2017
  !
  
  implicit none
  
  PetscInt :: num_cells
  type(chem_species_type) :: chem_species
  
  allocate(chem_species%initial_conc_mol(num_cells))
  allocate(chem_species%current_conc_mol(num_cells))
  allocate(chem_species%current_conc_kg(num_cells))
  allocate(chem_species%inst_rate(num_cells))
  
  if (uninitialized(chem_species%in_panel)) then
    chem_species%in_panel = 0.d0
  endif
  
  chem_species%initial_conc_mol(:) = chem_species%in_panel
  chem_species%current_conc_mol(:) = chem_species%in_panel
  chem_species%current_conc_kg(:) = 0.d0  ! just a placeholder
  chem_species%inst_rate(:) = 0.d0        ! just a placeholder
  
end subroutine ChemSpeciesAllocate

! ************************************************************************** !

subroutine PMWSSInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/14/2017
  !
  
  use Data_Mediator_Vec_class
  use Realization_Base_class
  
  implicit none
  
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_wipp_srcsink_type) :: this
  
  IS :: is
  PetscInt :: i, j, k
  PetscErrorCode :: ierr
  PetscReal :: size_of_vec
  PetscInt, allocatable :: dofs_in_residual(:)
  type(srcsink_panel_type), pointer :: cur_waste_panel
  
  ierr = 0
  
  ! set up mass transfer vector
  call RealizCreateFlowMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%flow_data_mediator_list)
  
  ! create a Vec sized by # waste panels * # waste panel cells in region *
  ! # src/sinks (water [mol/s], gas [mol/s], energy [MJ/s])
  cur_waste_panel => this%waste_panel_list
  size_of_vec = 0
  do
    if (.not.associated(cur_waste_panel)) exit
    size_of_vec = size_of_vec + (cur_waste_panel%region%num_cells * &
                                 this%option%nflowdof)
    cur_waste_panel => cur_waste_panel%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)
  
  if (size_of_vec > 0) then
    allocate(dofs_in_residual(size_of_vec))
    dofs_in_residual = 0
    i = 0
    cur_waste_panel => this%waste_panel_list
    do
      if (.not.associated(cur_waste_panel)) exit
      do k = 1,cur_waste_panel%region%num_cells
        do j = 1,this%option%nflowdof
          i = i + 1
          dofs_in_residual(i) = &
            (cur_waste_panel%region%cell_ids(k)-1)*this%option%nflowdof + j
        enddo
      enddo
      cur_waste_panel => cur_waste_panel%next
    enddo
    ! zero-based indexing
    dofs_in_residual(:) = dofs_in_residual(:) - 1
    ! index to global petsc ordering
    dofs_in_residual(:) = dofs_in_residual(:) + &
               this%realization%patch%grid%global_offset * this%option%nflowdof
  endif
  ! create the index set (IS)
  call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                       dofs_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  if (allocated(dofs_in_residual)) deallocate(dofs_in_residual)
  ! load the data mediator vec scatter context with the IS
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%flow_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  
  call PMWSSSolve(this,0.d0,ierr)
  
end subroutine PMWSSInitializeRun

! *************************************************************************** !

subroutine PMWSSInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  use Global_Aux_module
  use Material_Aux_class
  use Option_module
  use Grid_module
  
  implicit none

  class(pm_wipp_srcsink_type) :: this
  
  type(srcsink_panel_type), pointer :: cur_waste_panel
  PetscReal :: dt

  dt = this%option%flow_dt
  
    if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," WIPP SRC/SINK PANEL MODEL ",51("="))')
  endif
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    call UpdateInventory(cur_waste_panel%inventory, &
                         cur_waste_panel%region%num_cells,dt)
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine PMWSSInitializeTimestep

! *************************************************************************** !

subroutine UpdateInventory(inventory,num_cells,dt)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 02/10/2017
  !
 
  implicit none
 
  type(inventory_type) :: inventory
  PetscInt :: num_cells
  PetscReal :: dt
 
  PetscReal :: molar_mass  ! [kg/mol]
 
  molar_mass = (55.845d-3) ! [Fe kg/mol]
  call UpdateChemSpecies(inventory%Fe_s,molar_mass,num_cells,dt)
  
  molar_mass = (55.845d-3 + 2.d0*15.9994d-3 + 2.d0*1.01d-3) ! [FeOH2 kg/mol]
  call UpdateChemSpecies(inventory%FeOH2_s,molar_mass,num_cells,dt)
  
  molar_mass = (6.d0*12.0107d-3 + 10.d0*1.01d-3 + &
                5.d0*15.9994d-3) ! [C6H10O5 kg/mol]
  call UpdateChemSpecies(inventory%C6H10O5_s,molar_mass,num_cells,dt)
  
  molar_mass = (1.01d-3) ! [H+ kg/mol] 
  call UpdateChemSpecies(inventory%H_ion_aq,molar_mass,num_cells,dt)
  
  molar_mass = (14.0067d-3 + 3.d0*15.9994d-3) ! [NO3- kg/mol]
  call UpdateChemSpecies(inventory%NO3_minus_aq,molar_mass,num_cells,dt)
  
  molar_mass = (2.d0*14.0067d-3) ! [N2 kg/mol]
  call UpdateChemSpecies(inventory%N2_g,molar_mass,num_cells,dt)
  
  molar_mass = (12.0107d-3 + 2.d0*15.9994d-3) ! [CO2 kg/mol]
  call UpdateChemSpecies(inventory%CO2_g,molar_mass,num_cells,dt)
  
  molar_mass = (32.065d-3 + 4.d0*15.9994d-3) ! [SO42- kg/mol]
  call UpdateChemSpecies(inventory%SO42_minus_aq,molar_mass,num_cells,dt)
  
  molar_mass = (2.d0*1.01d-3 + 32.065d-3) ! [H2S kg/mol]
  call UpdateChemSpecies(inventory%H2S_g,molar_mass,num_cells,dt)
  
  molar_mass = (55.846d-3 + 32.065d-3) ! [FeS kg/mol]
  call UpdateChemSpecies(inventory%FeS_s,molar_mass,num_cells,dt)
  
  molar_mass = (24.305d-3 + 15.9994d-3) ! [MgO kg/mol]
  call UpdateChemSpecies(inventory%MgO_s,molar_mass,num_cells,dt)
  
  molar_mass = (24.305d-3 + 2.d0*15.9994d-3 + 2.d0*1.01d-3) ! [Mg(OH)2 kg/mol]
  call UpdateChemSpecies(inventory%MgOH2_s,molar_mass,num_cells,dt)
  
  molar_mass = (5.d0*24.305d-3 + 4.d0*12.0107d-3 + 4.d0*3.d0*15.9994d-3 + &
                2.d0*15.9994d-3 + 2.d0*1.01d-3 + 8.d0*1.01d-3 + &
                4.d0*15.9994d-3) ! [Mg5CO34OH24H2 kg/mol]
  call UpdateChemSpecies(inventory%Mg5CO34OH24H2_s,molar_mass,num_cells,dt)
                
  molar_mass = (24.305d-3 + 12.0107d-3 + 3.d0*15.9994d-3) ! [MgCO3 kg/mol]
  call UpdateChemSpecies(inventory%MgCO3_s,molar_mass,num_cells,dt)
                                      
 end subroutine UpdateInventory
 
! *************************************************************************** !

subroutine UpdateChemSpecies(chem_species,molar_mass,num_cells,dt)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  implicit none
  
  type(chem_species_type) :: chem_species
  PetscReal :: molar_mass  ! [kg/mol]
  PetscInt :: num_cells
  PetscReal :: dt          ! [sec]
  
  PetscInt :: k
  
  do k = 1,num_cells
    ! [mol/m3]
    chem_species%current_conc_mol(k) = &
                 chem_species%current_conc_mol(k) + &   ! [mol/m3]
                 chem_species%inst_rate(k) * &          ! [mol/m3/sec]
                 dt                                     ! [sec]
    ! [kg/m3]
    chem_species%current_conc_kg(k) = &
                 chem_species%current_conc_mol(k) * &   ! [mol/m3]
                 molar_mass                             ! [kg/mol]
  enddo
  
end subroutine UpdateChemSpecies

! *************************************************************************** !

 subroutine PMWSSSolve(this,time,ierr)
  ! 
  ! .
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  use Option_module
  use Grid_module
  use General_Aux_module
  use Material_Aux_class
  use Global_Aux_module
  use EOS_Gas_module
  use EOS_Water_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_wipp_srcsink_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(general_auxvar_type), pointer :: gen_auxvar(:,:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_p(:)
  type(srcsink_panel_type), pointer :: cur_waste_panel
  PetscInt :: i, j, k
  PetscInt :: cell_id
  PetscInt :: num_cells
  ! brine/gas generation variable
  PetscReal, allocatable :: water_saturation(:)
  PetscReal, allocatable :: s_eff(:)
  PetscReal, allocatable :: rxnrate_corrosion(:)
  PetscReal, allocatable :: rxnrate_biodeg(:)
  PetscReal, allocatable :: rxnrate_FeS(:)
  PetscReal, allocatable :: rxnrate_mgoh2(:)
  PetscReal, allocatable :: rxnrate_hydromag(:)
  PetscReal, allocatable :: rxnrate_hymagcon(:)
  ! enthalpy calculation variables
  PetscReal :: temperature
  PetscReal :: pressure_liq
  PetscReal :: pressure_gas
  PetscReal :: H_liq
  PetscReal :: H_gas
  PetscReal :: U_gas
  PetscReal :: gas_energy
  PetscReal :: brine_energy
  
  option => this%realization%option
  grid => this%realization%patch%grid
  gen_auxvar => this%realization%patch%aux%General%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  global_auxvar => this%realization%patch%aux%Global%auxvars
  
  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
  j = 0  ! (j indexes the data mediator vec)
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    num_cells = cur_waste_panel%region%num_cells
    !-----allocate-arrays----------------------------------------------------
    allocate(water_saturation(num_cells))
    water_saturation(:) = 0.d0
    allocate(s_eff(num_cells))
    s_eff(:) = 0.d0
    allocate(rxnrate_corrosion(num_cells))
    rxnrate_corrosion(:) = 0.d0
    allocate(rxnrate_biodeg(num_cells))
    rxnrate_biodeg(:) = 0.d0
    allocate(rxnrate_FeS(num_cells))
    rxnrate_FeS(:) = 0.d0
    allocate(rxnrate_mgoh2(num_cells))
    rxnrate_mgoh2(:) = 0.d0
    allocate(rxnrate_hydromag(num_cells))
    rxnrate_hydromag(:) = 0.d0
    allocate(rxnrate_hymagcon(num_cells))
    rxnrate_hymagcon(:) = 0.d0

    do i = 1,num_cells
    !-----effective-brine-saturation------------------------------------------
      water_saturation(i) = &
        gen_auxvar(ZERO_INTEGER,grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
        sat(option%liquid_phase)
      s_eff(i) = water_saturation(i) - this%smin + this%satwick*(1.d0 - &
                 exp(200.d0*this%alpharxn* &
                 (max((water_saturation(i)-this%smin),0.d0))**2.d0))
    !-----anoxic-iron-corrosion-----------------------------------------------
      rxnrate_corrosion(i) = (this%inundated_corrosion_rate*s_eff(i)) + &
        ((this%humid_corrosion_factor*this%inundated_corrosion_rate) * &
        (1.d0-s_eff(i)))
    !-----biodegradation------------------------------------------------------
      rxnrate_biodeg(i) = (this%inundated_biodeg_rate*s_eff(i)) + &
        ((this%humid_biodeg_factor*this%inundated_biodeg_rate) * &
        (1.d0-s_eff(i)))
    !-----iron-sulfidation----------------------------------------------------
      rxnrate_FeS(i) = rxnrate_biodeg(i)*this%RXH2S_factor
    !-----MgO-hydration-------------------------------------------------------
      rxnrate_mgoh2(i) = (this%inundated_brucite_rate*s_eff(i)) + &
        ((this%humid_brucite_rate)*(1.d0-s_eff(i)))
    !-----hydromagnesite------------------------------------------------------
      rxnrate_hydromag(i) = rxnrate_biodeg(i)*this%RXCO2_factor
    !-----hydromagnesite-conversion-------------------------------------------
      rxnrate_hymagcon(i) = this%hymagcon_rate* &
        cur_waste_panel%inventory%Mg5CO34OH24H2_s%current_conc_kg(i)
    !-----tracked-species-[mol/m3/sec]----------------------------------------
      cur_waste_panel%inventory%FeOH2_s%inst_rate(i) = &
          1.d0*rxnrate_corrosion(i) + (-1.d0*rxnrate_FeS(i))
      cur_waste_panel%inventory%Fe_s%inst_rate(i) = &
          (-1.d0*rxnrate_corrosion(i)) + (-1.d0*rxnrate_FeS(i)) 
      cur_waste_panel%inventory%FeS_s%inst_rate(i) = &
          1.d0*rxnrate_FeS(i) + 1.d0*rxnrate_FeS(i)
      cur_waste_panel%inventory%C6H10O5_s%inst_rate(i) = &
          (-1.d0*rxnrate_biodeg(i)) + (-1.d0*rxnrate_biodeg(i))
      cur_waste_panel%inventory%H_ion_aq%inst_rate(i) = &
          (-4.8d0*rxnrate_biodeg(i)) + (-6.d0*rxnrate_biodeg(i))
      cur_waste_panel%inventory%NO3_minus_aq%inst_rate(i) = &
          (-4.8d0*rxnrate_biodeg(i))
      cur_waste_panel%inventory%CO2_g%inst_rate(i) = &
          6.d0*rxnrate_biodeg(i) + 6.d0*rxnrate_biodeg(i) + &
          (-4.d0*rxnrate_hydromag(i))
      cur_waste_panel%inventory%N2_g%inst_rate(i) = &
          2.4d0*rxnrate_biodeg(i)
      cur_waste_panel%inventory%SO42_minus_aq%inst_rate(i) = &
          (-3.d0*rxnrate_biodeg(i))
      cur_waste_panel%inventory%H2S_g%inst_rate(i) = &
          3.d0*rxnrate_biodeg(i) + (-1.d0*rxnrate_FeS(i)) + &
          (-1.d0*rxnrate_FeS(i))
      cur_waste_panel%inventory%MgO_s%inst_rate(i) = &
          (-1.d0*rxnrate_mgoh2(i))
      cur_waste_panel%inventory%MgOH2_s%inst_rate(i) = &
          1.d0*rxnrate_mgoh2(i) + (-5.d0*rxnrate_hydromag(i)) + &
          1.d0*rxnrate_hymagcon(i)
      cur_waste_panel%inventory%Mg5CO34OH24H2_s%inst_rate(i) = &
          1.d0*rxnrate_hydromag(i) + (-1.d0*rxnrate_hymagcon(i)) 
      cur_waste_panel%inventory%MgCO3_s%inst_rate(i) = &
          4.d0*rxnrate_hymagcon(i) 
    !-------------------------------------------------------------------------
    ! adjust rates in case a reactant will get used up
    ! is this even possible, because it requires knowing the next time step
    !-----gas-generation-[molH2/m3/sec]---------------------------------------
      cur_waste_panel%gas_generation_rate(i) = &
          1.d0*rxnrate_corrosion(i) + 1.d0*rxnrate_FeS(i)
    !-----brine-generation-[molH2O/m3/sec]------------------------------------
      cur_waste_panel%brine_generation_rate(i) = &
          (-2.d0*rxnrate_corrosion(i)) + 7.4d0*rxnrate_biodeg(i) + &
          5.d0*rxnrate_biodeg(i) + 2.d0*rxnrate_FeS(i) + &
          (-1.d0*rxnrate_mgoh2(i)) + 4.d0*rxnrate_hymagcon(i)
    !------source-term-calculation--------------------------------------------
      cell_id = grid%nL2G(cur_waste_panel%region%cell_ids(i))
      j = j + 1
      ! liquid source term [mol/sec]
      vec_p(j) = cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
                 material_auxvars(cell_id)%volume              ! [m3] 
      j = j + 1
      ! gas source term [mol/sec]
      vec_p(j) = cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
                 material_auxvars(cell_id)%volume              ! [m3]
      j = j + 1
      ! energy source term [MJ/sec]; H from EOS [J/kmol]
      brine_energy = 0.d0
      gas_energy = 0.d0
      temperature = gen_auxvar(ZERO_INTEGER, &
          grid%nL2G(cur_waste_panel%region%cell_ids(i)))%temp
      select case(global_auxvar(cell_id)%istate)
        case(GAS_STATE) !------------------------------------------------------
          pressure_gas = gen_auxvar(ZERO_INTEGER, &
              grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
              pres(option%gas_phase)
          call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
          gas_energy = &
              cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
        case(LIQUID_STATE) !---------------------------------------------------
          pressure_liq = gen_auxvar(ZERO_INTEGER, &
              grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
              pres(option%liquid_phase)
          call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
          brine_energy = &
              cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3] 
              H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
        case(TWO_PHASE_STATE) !------------------------------------------------
          pressure_liq = gen_auxvar(ZERO_INTEGER, &
              grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
              pres(option%liquid_phase)
          pressure_gas = gen_auxvar(ZERO_INTEGER, &
              grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
              pres(option%gas_phase)
          call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
          call EOSGasEnergy(temperature,pressure_gas,H_gas,U_gas,ierr)
          brine_energy = &
              cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3] 
              H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
          gas_energy = &
              cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
      end select
      vec_p(j) = brine_energy + gas_energy
    enddo
    !-----deallocate-arrays---------------------------------------------------
    deallocate(water_saturation)
    deallocate(s_eff)
    deallocate(rxnrate_corrosion)
    deallocate(rxnrate_biodeg)
    deallocate(rxnrate_FeS)
    deallocate(rxnrate_mgoh2)
    deallocate(rxnrate_hydromag)
    deallocate(rxnrate_hymagcon)
    !-------------------------------------------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWSSSolve

! ************************************************************************** !

subroutine PMWSSFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/22/2017

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
end subroutine PMWSSFinalizeTimestep

! *************************************************************************** !

 subroutine PMWSSOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  
  implicit none

  class(pm_wipp_srcsink_type) :: this
  
end subroutine PMWSSOutput

! *************************************************************************** !

subroutine PMWSSInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  ! 
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMWSSInputRecord

! *************************************************************************** !

subroutine PMWSSDestroy(this)
  ! 
  ! Destroys the WIPP source/sink process model
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  !call PMWSSStrip(this)
  !deallocate gas_generation_rate, brine_generation_rate, scaling_factor
  
end subroutine PMWSSDestroy

! *************************************************************************** !

end module PM_WIPP_SrcSink_class