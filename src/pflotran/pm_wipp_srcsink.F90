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
    PetscReal, pointer :: initial_conc_mol(:)  ! [mol/m3-bulk]
    PetscReal, pointer :: inst_rate(:)         ! [mol/m3-bulk/sec]
    PetscReal, pointer :: current_conc_mol(:)  ! [mol/m3-bulk]
    PetscReal, pointer :: current_conc_kg(:)   ! [kg/m3-bulk] per BRAGFLO U.M.
    PetscReal :: molar_mass                    ! [kg/mol]
    PetscReal :: tot_mass_in_panel             ! [kg/panel-volume]
  end type chem_species_type
  
  type, public :: inventory_type
    character(len=MAXWORDLENGTH) :: name
    type(chem_species_type) :: Fe_s 
    type(chem_species_type) :: FeOH2_s
    type(chem_species_type) :: C6H10O5_s
    type(chem_species_type) :: RuPl_s
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
    PetscReal :: num_drums_packing
    type(pre_inventory_type), pointer :: preinventory
  end type inventory_type
  
  type, public :: pre_inventory_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: Fe_in_panel          ! [total initial kg in waste panel]
    PetscReal :: MgO_in_panel         ! [total initial kg in waste panel]
    PetscReal :: Cellulose_in_panel   ! [total initial kg in waste panel]
    PetscReal :: RubberPlas_in_panel  ! [total initial kg in waste panel]
    PetscReal :: H_ion_in_panel       ! [total initial kg in waste panel]
    PetscReal :: Nitrate_in_panel     ! [total initial kg in waste panel]
    PetscReal :: Sulfate_in_panel     ! [total initial kg in waste panel]
    PetscReal :: num_drums_packing
    type(pre_inventory_type), pointer :: next
  end type pre_inventory_type
  
  type, public :: srcsink_panel_type
    character(len=MAXWORDLENGTH) :: name
    type(region_type), pointer :: region
    character(len=MAXWORDLENGTH) :: region_name
    type(inventory_type) :: inventory
    character(len=MAXWORDLENGTH) :: inventory_name
    PetscReal, pointer :: scaling_factor(:)        ! [-]
    PetscReal, pointer :: gas_generation_rate(:)   ! [mol/m3-bulk/sec]
    PetscReal, pointer :: brine_generation_rate(:) ! [mol/m3-bulk/sec]
    PetscReal :: inundated_corrosion_rate          ! [mol/m3-bulk/sec]
    PetscReal :: humid_corrosion_rate              ! [mol/m3-bulk/sec], [-]
    PetscReal :: inundated_biodeg_rate             ! [mol/m3-bulk/sec]
    PetscReal :: humid_biodeg_rate                 ! [mol/m3-bulk/sec]
    PetscReal :: inundated_brucite_rate            ! [mol/m3-bulk/sec]
    PetscReal :: humid_brucite_rate                ! [mol/m3-bulk/sec]
    PetscReal :: RXH2S_factor                      ! [-]
    PetscReal :: volume                            ! [m3]
    PetscInt :: id
    PetscMPIInt :: myMPIcomm
    PetscMPIInt :: myMPIgroup
    PetscInt, pointer :: rank_list(:)
    type(srcsink_panel_type), pointer :: next
  end type srcsink_panel_type

  type, public, extends(pm_base_type) :: pm_wipp_srcsink_type
    PetscReal :: alpharxn           ! [-] 
    PetscReal :: smin               ! [-]
    PetscReal :: satwick            ! [-]
    PetscReal :: corrmco2           ! [m/s]
    PetscReal :: humcorr            ! [m/s]
    PetscReal :: gratmici           ! [mol/kg/sec]
    PetscReal :: gratmich           ! [mol/kg/sec]
    PetscReal :: brucitei           ! [mol/kg/sec]
    PetscReal :: bruciteh           ! [mol/kg/sec]
    PetscReal :: RXCO2_factor
    PetscReal :: hymagcon_rate      ! [mol/kg/sec]
    PetscReal :: drum_surface_area  ! [m2/drum]
    PetscReal :: biogenfc           ! [-]
    type(srcsink_panel_type), pointer :: waste_panel_list
    type(pre_inventory_type), pointer :: pre_inventory_list
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
  
  interface TaperRxnrate
    module procedure TaperRxnrate1
    module procedure TaperRxnrate2
  end interface

  public :: PMWSSCreate, &
            WastePanelCreate, &
            PreInventoryCreate

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
  nullify(PMWSSCreate%pre_inventory_list)
  nullify(PMWSSCreate%data_mediator)
  nullify(PMWSSCreate%realization)
  PMWSSCreate%name = 'wipp source sink'
  PMWSSCreate%alpharxn = UNINITIALIZED_DOUBLE
  PMWSSCreate%smin = UNINITIALIZED_DOUBLE
  PMWSSCreate%satwick = UNINITIALIZED_DOUBLE
  PMWSSCreate%corrmco2 = UNINITIALIZED_DOUBLE
  PMWSSCreate%humcorr = UNINITIALIZED_DOUBLE
  PMWSSCreate%gratmici = UNINITIALIZED_DOUBLE
  PMWSSCreate%gratmich = UNINITIALIZED_DOUBLE
  PMWSSCreate%brucitei = UNINITIALIZED_DOUBLE
  PMWSSCreate%bruciteh = UNINITIALIZED_DOUBLE
  PMWSSCreate%RXCO2_factor = UNINITIALIZED_DOUBLE
  PMWSSCreate%hymagcon_rate = UNINITIALIZED_DOUBLE
  PMWSSCreate%drum_surface_area = UNINITIALIZED_DOUBLE
  PMWSSCreate%biogenfc = UNINITIALIZED_DOUBLE
  
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
  nullify(WastePanelCreate%scaling_factor)
  nullify(WastePanelCreate%gas_generation_rate)
  nullify(WastePanelCreate%brine_generation_rate)
  nullify(WastePanelCreate%rank_list)
  call InventoryInit(WastePanelCreate%inventory)
  WastePanelCreate%name = ''
  WastePanelCreate%region_name = ''
  WastePanelCreate%inventory_name = ''
  WastePanelCreate%volume = UNINITIALIZED_DOUBLE
  WastePanelCreate%inundated_corrosion_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%humid_corrosion_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%inundated_biodeg_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%humid_biodeg_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%inundated_brucite_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%humid_brucite_rate = UNINITIALIZED_DOUBLE
  WastePanelCreate%RXH2S_factor = UNINITIALIZED_DOUBLE
  WastePanelCreate%id = 0
  WastePanelCreate%myMPIgroup = 0
  WastePanelCreate%myMPIcomm = 0

end function WastePanelCreate

! *************************************************************************** !

function PreInventoryCreate()
  !
  ! Creates and initializes a waste panel pre-inventory.
  !
  ! Author: Jenn Frederick
  ! Date: 3/01/2017
  !
  
  implicit none
  
  type(pre_inventory_type), pointer :: PreInventoryCreate
  
  allocate(PreInventoryCreate)
  nullify(PreInventoryCreate%next)

  PreInventoryCreate%name = ''
  PreInventoryCreate%Fe_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%MgO_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%Cellulose_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%RubberPlas_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%H_ion_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%Nitrate_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%Sulfate_in_panel = UNINITIALIZED_DOUBLE
  PreInventoryCreate%num_drums_packing = UNINITIALIZED_DOUBLE

end function PreInventoryCreate

! *************************************************************************** !

subroutine InventoryInit(inventory)
  !
  ! Initializes a waste panel inventory object.
  !
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
  type(inventory_type) :: inventory
  
  PetscReal :: molar_mass
  
  nullify(inventory%preinventory)
  inventory%name = ''
  
  molar_mass = (0.055847) ! [kg/mol MW_FE]
  call InitChemSpecies(inventory%Fe_s,molar_mass)  ! iron
  
  molar_mass = (0.08986) ! [kg/mol MW_FEOH2]
  call InitChemSpecies(inventory%FeOH2_s,molar_mass) ! iron hydroxide
  
  molar_mass = (0.027023) ! [kg/mol MW_CELL]
  call InitChemSpecies(inventory%C6H10O5_s,molar_mass)  ! cellulose
  
  molar_mass = (0.027023) ! [rubber/plastic kg/mol, same as MW_CELL]
  call InitChemSpecies(inventory%RuPl_s,molar_mass)  ! rubber/plastics
  
  molar_mass = (1.01d-3) ! [H+ kg/mol] 
  call InitChemSpecies(inventory%H_ion_aq,molar_mass)  ! h+
  
  molar_mass = (14.0067d-3 + 3.d0*15.9994d-3) ! [NO3- kg/mol]
  call InitChemSpecies(inventory%NO3_minus_aq,molar_mass)  ! nitrate
  
  molar_mass = (0.0440098) ! [kg/mol MW_CO2]
  call InitChemSpecies(inventory%CO2_g,molar_mass)  ! carbon dioxide gas 
  
  molar_mass = (0.02801348) ! [kg/mol MW_N2]
  call InitChemSpecies(inventory%N2_g,molar_mass)  ! nitrogen gas  
  
  molar_mass = (32.065d-3 + 4.d0*15.9994d-3) ! [SO42- kg/mol]
  call InitChemSpecies(inventory%SO42_minus_aq,molar_mass)  ! sulfate
  
  molar_mass = (0.03408188) ! [kg/mol MW_H2S]
  call InitChemSpecies(inventory%H2S_g,molar_mass)  
  
  molar_mass = (0.087911) ! [kg/mol MW_FES]
  call InitChemSpecies(inventory%FeS_s,molar_mass)  ! iron sulfide 
  
  molar_mass = (0.040304) ! [kg/mol MW_MGO]
  call InitChemSpecies(inventory%MgO_s,molar_mass)  ! magnesium oxide
  
  molar_mass = (0.05832) ! [kg/mol MW_MGOH2]
  call InitChemSpecies(inventory%MgOH2_s,molar_mass)  ! magnesium hydroxide              
  
  molar_mass = (5.d0*24.305d-3 + 4.d0*12.0107d-3 + 4.d0*3.d0*15.9994d-3 + &
                2.d0*15.9994d-3 + 2.d0*1.01d-3 + 8.d0*1.01d-3 + &
                4.d0*15.9994d-3) ! [Mg5CO34OH24H2 kg/mol]
  call InitChemSpecies(inventory%Mg5CO34OH24H2_s,molar_mass)  ! hydromagnesite
  
  molar_mass = (0.084314) ! [kg/mol MW_MGCO3]
  call InitChemSpecies(inventory%MgCO3_s,molar_mass)  ! magnesium carbonate

end subroutine InventoryInit

! *************************************************************************** !

subroutine InitChemSpecies(chem_species,molar_mass)
  !
  ! Initializes a waste panel inventory's chemical species.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  implicit none
  
  type(chem_species_type) :: chem_species
  PetscReal :: molar_mass
  
  nullify(chem_species%current_conc_mol)        ! [mol/m3]
  nullify(chem_species%current_conc_kg)         ! [kg/m3]
  nullify(chem_species%initial_conc_mol)        ! [mol/m3]
  nullify(chem_species%inst_rate)               ! [mol/m3/sec]
  chem_species%molar_mass = molar_mass          ! [kg/mol]
  chem_species%tot_mass_in_panel = 0.d0         ! [kg/panel-volume]
  
end subroutine InitChemSpecies

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

subroutine AssociateRegion(this,region_list)
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
  
end subroutine AssociateRegion

! *************************************************************************** !

subroutine AssociateInventory(this)
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
  
  type(pre_inventory_type), pointer :: cur_preinventory
  class(srcsink_panel_type), pointer :: cur_waste_panel
  type(option_type), pointer :: option
  PetscBool :: matched
  
  option => this%option
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
      cur_preinventory => this%pre_inventory_list     
      do
        if (.not.associated(cur_preinventory)) exit
        matched = PETSC_FALSE
        if (StringCompare(trim(cur_preinventory%name), &
                          trim(cur_waste_panel%inventory_name))) then
          call CopyPreInvToInv(cur_preinventory,cur_waste_panel%inventory)
          matched = PETSC_TRUE
        endif
        if (matched) exit
        cur_preinventory => cur_preinventory%next
      enddo      
      if (.not.associated(cur_waste_panel%inventory%preinventory)) then
        option%io_buffer = 'WASTE_PANEL INVENTORY ' // &
                           trim(cur_waste_panel%inventory_name) // ' not found.'
        call printErrMsg(option)
      endif
    cur_waste_panel => cur_waste_panel%next
  enddo
  
end subroutine AssociateInventory

! *************************************************************************** !

subroutine CopyPreInvToInv(preinventory,inventory)
  ! 
  ! Copies information from a pre-inventory to an inventory object.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/02/2017
  !
  
  implicit none
  
  type(pre_inventory_type), pointer :: preinventory
  type(inventory_type) :: inventory
  
  inventory%name = preinventory%name
  inventory%num_drums_packing = preinventory%num_drums_packing
  inventory%preinventory => preinventory
  
end subroutine CopyPreInvToInv

! *************************************************************************** !

subroutine SetRegionScaling(this,waste_panel)
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
  waste_panel%volume = total_volume_global
  
end subroutine SetRegionScaling

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
  type(pre_inventory_type), pointer :: new_inventory
  type(pre_inventory_type), pointer :: cur_preinventory
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
      case('SAMPLED_SAT_WICK')
        call InputReadDouble(input,option,this%satwick)
        call InputErrorMsg(input,option,'wicking saturation parameter &
                           &(SAT_WICK)',error_string)
      case('SAMPLED_CORRMCO2','CORRMCO2')
        call InputReadDouble(input,option,this%corrmco2)
        call InputErrorMsg(input,option,'inundated steel corrosion rate &
                           &(CORRMCO2)',error_string)
      case('HUMCORR')
        call InputReadDouble(input,option,this%humcorr)
        call InputErrorMsg(input,option,'humid steel corrosion rate &
                           &(HUMCORR)',error_string)
      case('SAMPLED_GRATMICI','GRATMICI')
        call InputReadDouble(input,option,this%gratmici)
        call InputErrorMsg(input,option,'inundated biodegradation rate for &
                           &cellulose (GRATMICI)',error_string)
      case('SAMPLED_GRATMICH','GRATMICH')
        call InputReadDouble(input,option,this%gratmich)
        call InputErrorMsg(input,option,'humid diodegradation rate for &
                           &cellulose (GRATMICH)',error_string)
      case('SAMPLED_BRUCITEI','BRUCITEI')
        call InputReadDouble(input,option,this%brucitei)
        call InputErrorMsg(input,option,'MgO inundated hydration rate in &
                           &brine (BRUCITEI)',error_string)
      case('SAMPLED_BRUCITEH','BRUCITEH')
        call InputReadDouble(input,option,this%bruciteh)
        call InputErrorMsg(input,option,'MgO humid hydration rate (BRUCITEH)', &
                           error_string)
      case('SAMPLED_HYMAGCON','HYMAGCON')
        call InputReadDouble(input,option,this%hymagcon_rate)
        call InputErrorMsg(input,option,'hydromagnesite to magnesite &
                           &conversion rate (HYMAGCON)',error_string)
      case('ASDRUM')
        call InputReadDouble(input,option,this%drum_surface_area)
        call InputErrorMsg(input,option,'surface area of corrodable metal &
                           &per drum (ASDRUM)',error_string)
      case('BIOGENFC')
        call InputReadDouble(input,option,this%biogenfc)
        call InputErrorMsg(input,option,'probability of attaining sampled &
                           &microbial gas generation rates (BIOGENFC)', &
                           error_string)
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
                 trim(error_string) // ' block. WASTE_PANEL name "' // &
                 trim(new_waste_panel%name) // '".'
          call printErrMsg(option)
        endif
        if (new_waste_panel%inventory_name == '') then
          option%io_buffer = 'INVENTORY must be specified in the ' // &
                 trim(error_string) // ' block. WASTE_PANEL name "' // &
                 trim(new_waste_panel%name) // '".'
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
        new_inventory => PreInventoryCreate()
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
                  case('WTFETOT')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial Fe mass &
                                       &(WTFETOT)',error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string) // ',initial Fe mass (WTFETOT) &
                          &units',option)
                    new_inventory%Fe_in_panel = double
                !-----------------------------
                  case('WTMGOTOT')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial MgO mass &
                                       &(WTMGOTOT)',error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string) // ',initial MgO mass units &
                         &(WTMGOTOT)',option)
                    new_inventory%MgO_in_panel = double
                !-----------------------------
                  case('WTCELTOT')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial cellulose mass &
                                       &(WTCELTOT)',error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string) // ',initial cellulose mass &
                         &(WTCELTOT) units',option)
                    new_inventory%Cellulose_in_panel = double
                !-----------------------------
                  case('WTRPLTOT')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial rubber and &
                                       &plastics mass (WTRPLTOT)',error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                         trim(error_string) // ',initial rubber and plastics &
                         &mass (WTRPLTOT) units',option)
                    new_inventory%RubberPlas_in_panel = double
                !-----------------------------------
                  case('DRROOM')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'number of metal drums per &
                                 &panel in ideal packing (DRROOM)',error_string)
                    new_inventory%num_drums_packing = double
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
                    call InputErrorMsg(input,option,'initial H+ mass', &
                                       error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                          trim(error_string) // ',initial H+ mass units',option)
                    new_inventory%H_ion_in_panel = double
                !-----------------------------
                  case('NITRATE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial nitrate mass', &
                                       error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                      trim(error_string) // ',initial nitrate mass units', &
                      option)
                    new_inventory%Nitrate_in_panel = double
                !-----------------------------
                  case('SULFATE')
                    call InputReadDouble(input,option,double)
                    call InputErrorMsg(input,option,'initial sulfate mass', &
                                       error_string)
                    call InputReadAndConvertUnits(input,double,'kg', &
                      trim(error_string) // ',initial sulfate mass units', &
                      option)
                    new_inventory%Sulfate_in_panel = double
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
        if (Uninitialized(new_inventory%Fe_in_panel)) then
          option%io_buffer = 'Initial Fe (solid) inventory must be specified &
                        &using the SOLIDS,WTFETOT keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%MgO_in_panel)) then
          option%io_buffer = 'Initial MgO (solid) inventory must be specified &
                        &using the SOLIDS,WTMGOTOT keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%Cellulose_in_panel)) then
          option%io_buffer = 'Initial cellulose (solid) inventory must be &
                        &specified using the SOLIDS,WTCELTOT keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%RubberPlas_in_panel)) then
          option%io_buffer = 'Initial rubber/plastic (solid) inventory must be &
                        &specified using the SOLIDS,WTRPLTOT keyword in &
                        &the WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%H_ion_in_panel)) then
          option%io_buffer = 'Initial H+ (aqueous) inventory must be specified &
                        &using the AQUEOUS,H+ keyword in the WIPP_SOURCE_SINK &
                        &block. Inventory name "' // trim(new_inventory%name) &
                         // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%Nitrate_in_panel)) then
          option%io_buffer = 'Initial nitrate (aqueous) inventory must be &
                        &specified using the AQUEOUS,NITRATE keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%Sulfate_in_panel)) then
          option%io_buffer = 'Initial sulfate (aqueous) inventory must be &
                        &specified using the AQUEOUS,SULFATE keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        if (Uninitialized(new_inventory%num_drums_packing)) then
          option%io_buffer = 'Number of metal drums must be specified in the &
                        &specified using the SOLIDS,DRROOM keyword in the &
                        &WIPP_SOURCE_SINK block. Inventory name "' // &
                        trim(new_inventory%name) // '".'
          call printErrMsg(option)
        endif
        added = PETSC_FALSE
        if (.not.associated(this%pre_inventory_list)) then
          this%pre_inventory_list => new_inventory
        else
          cur_preinventory => this%pre_inventory_list
          do
            if (.not.associated(cur_preinventory)) exit
            if (.not.associated(cur_preinventory%next)) then
              cur_preinventory%next => new_inventory
              added = PETSC_TRUE
            endif
            if (added) exit
            cur_preinventory => cur_preinventory%next
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
  if (.not.associated(this%pre_inventory_list)) then
    option%io_buffer = 'At least one INVENTORY must be specified in the &
                       &WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%alpharxn)) then
    option%io_buffer = 'ALPHARXN must be specified in the WIPP_SOURCE_SINK &
                       &block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%smin)) then
    option%io_buffer = 'SOCMIN must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%satwick)) then
    option%io_buffer = 'SAT_WICK (wicking saturation parameter) must be &
                       &specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%gratmici)) then
    option%io_buffer = 'GRATMICI (inundated biodegradation rate for cellulose) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%brucitei)) then
    option%io_buffer = 'BRUCITEI (MgO inundated hydration rate in brine) must &
                       &be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%corrmco2)) then
    option%io_buffer = 'CORRMCO2 (inundated steel corrosion rate) must be &
                       &specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%gratmich)) then
    option%io_buffer = 'GRATMICH (humid biodegradation rate for cellulose) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%bruciteh)) then
    option%io_buffer = 'BRUCITEH (MgO humid hydration rate) must be specified &
                       &in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%humcorr)) then
    option%io_buffer = 'HUMCORR (humid steel corrosion rate) must be specified &
                       &in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%hymagcon_rate)) then
    option%io_buffer = 'HYMAGCON (hydromagnesite to magnesite conversion rate) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%drum_surface_area)) then
    option%io_buffer = 'ASDRUM (metal drum surface area) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%biogenfc)) then
    option%io_buffer = 'BIOGENFC (microbial gas generation probability) &
                       &must be specified in the WIPP_SOURCE_SINK block.'
    call printErrMsg(option)
  endif
  
  call AssociateInventory(this)
  
end subroutine PMWSSRead

! *************************************************************************** !

subroutine ProcessAfterRead(this)
  !
  ! After reading input parameters, ALGEBRA processing is done to get final 
  ! input parameters required.
  !
  ! Author: Jenn Frederick
  ! Date: 3/10/2017
  !

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(srcsink_panel_type), pointer :: cur_waste_panel
  type(pre_inventory_type), pointer :: preinventory
  PetscReal, parameter :: DN_FE = 7870.d0   ! [kg/m3] density of iron
  PetscReal, parameter :: MW_FE = 5.5847d-2 ! [kg/mol] mol weight of iron
  PetscReal, parameter :: MW_MGO = 4.03d-2  ! [kg/mol] mol weight of MgO
  PetscReal, parameter :: MW_C = 2.70d-2    ! [kg/mol] mol weight of cellulosics
  PetscReal :: D_c                          ! [kg/m3] mass conc biodegradables
  PetscReal :: D_m                          ! [kg/m3] mass conc MgO
  PetscReal :: D_s                          ! [m2/m3] area conc iron steel
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    preinventory => cur_waste_panel%inventory%preinventory
    !-----mass-concentrations----------------------------------units---------
    D_c = (preinventory%Cellulose_in_panel + &               ! [kg]
           preinventory%RubberPlas_in_panel) / &             ! [kg]
          cur_waste_panel%volume                             ! [m3]
          ! D_c = (m_c + m_r + 1.7*m_p)/volume ???
    D_s = this%drum_surface_area * &                         ! [m2]
          cur_waste_panel%inventory%num_drums_packing / &    ! [-]
          cur_waste_panel%volume                             ! [m3]
    D_m = 1.2d0*D_c * &                                      ! [kg/m3]
          MW_MGO / &                                         ! [kg/mol] 
          MW_C                                               ! [kg/mol]
    !-------------------------------------------------------------------------
    !-----anoxic-iron-corrosion--------------------------------units----------
    cur_waste_panel%inundated_corrosion_rate = &             ! [mol-Fe/m3/sec]
          this%corrmco2 * &                                  ! [m/s]
          D_s * &                                            ! [m2/m3]
          DN_FE / &                                          ! [kg/m3]
          MW_FE                                              ! [kg/mol]
    cur_waste_panel%humid_corrosion_rate = &                 ! [mol-Fe/m3/sec]
          this%humcorr * &                                   ! [m/s]
          D_s * &                                            ! [m2/m3]
          DN_FE / &                                          ! [kg/m3]
          MW_FE                                              ! [kg/mol]
    cur_waste_panel%humid_corrosion_rate = &                 ! [-]
          cur_waste_panel%humid_corrosion_rate / &
          cur_waste_panel%inundated_corrosion_rate
    !-----biodegradation------------------------------------units-------------
    cur_waste_panel%inundated_biodeg_rate = &             ! [mol-cell/m3/sec]
          this%gratmici * &                               ! [mol-cell/kg/sec]
          D_c * &                                         ! [kg/m3]
          this%biogenfc                                   ! [-]
    cur_waste_panel%humid_biodeg_rate = &                 ! [mol-cell/m3/sec]
          this%gratmich * &                               ! [mol-cell/kg/sec]
          D_c * &                                         ! [kg/m3]
          this%biogenfc                                   ! [-]            
    !-----iron-sulfidation----------------------------------------------------
    cur_waste_panel%RXH2S_factor = 1.11d-2  ! fill in later
           ! its 1 - (some ratio of initial nitrate to carbon loss via biodeg)
    !-----MgO-hydration-------------------------------------units-------------
    cur_waste_panel%inundated_brucite_rate = &            ! [mol-bruc/m3/sec]
          this%brucitei * &                               ! [mol-bruc/kg/sec]
          D_m                                             ! [kg/m3]
    cur_waste_panel%humid_brucite_rate = &                ! [mol-bruc/m3/sec]
          this%bruciteh * &                               ! [mol-bruc/kg/sec]
          D_m                                             ! [kg/m3]
    !-------------------------------------------------------------------------
    this%RXCO2_factor = 1.d0  ! BRAGFLO User's Manual Eq. 155, based on
                              ! Eqs. 145 & 146 stoichiometry 
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo

end subroutine ProcessAfterRead

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
  PetscInt :: i, j
  PetscBool :: local
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm_size
  PetscInt, pointer :: ranks(:)
  
  option => this%realization%option
  
  ! point the waste panel region to the desired region 
  call AssociateRegion(this,this%realization%patch%region_list)
  
  allocate(ranks(option%mycommsize))
  
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
    ranks(:) = 0
    newcomm_size = 0
    if (local) then
      cur_waste_panel%id = waste_panel_id
      ranks(option%myrank+1) = 1
    else
      cur_waste_panel%id = 0
      ranks(option%myrank+1) = 0
    endif
    ! count the number of processes that own the waste panel
    call MPI_Allreduce(MPI_IN_PLACE,ranks,option%mycommsize,MPI_INTEGER, &
                       MPI_SUM,option%mycomm,ierr)
    newcomm_size = sum(ranks)
    allocate(cur_waste_panel%rank_list(newcomm_size))
    j = 0
    do i = 1,option%mycommsize
      if (ranks(i) == 1) then
        j = j + 1
        cur_waste_panel%rank_list(j) = (i - 1)
      endif
    enddo
    ! create an MPI group and communicator for each waste panel
    call MPI_Group_incl(option%mygroup,newcomm_size,cur_waste_panel%rank_list, &
                        cur_waste_panel%myMPIgroup,ierr)
    call MPI_Comm_create(option%mycomm,cur_waste_panel%myMPIgroup, &
                         cur_waste_panel%myMPIcomm,ierr)
    if (local) then
      call SetRegionScaling(this,cur_waste_panel)
      call InventoryAllocate(cur_waste_panel%inventory, &
                             cur_waste_panel%region%num_cells, &
                             cur_waste_panel%volume)
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
  
  deallocate(ranks)
  
end subroutine PMWSSSetup

! ************************************************************************** !

subroutine InventoryAllocate(inventory,num_cells,volume)
  ! 
  ! Allocates the size of the chemical species arrays within the inventory to 
  ! the number of cells in the waste panel region, and assigns the initial
  ! values.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/23/2017
  !
  
  implicit none
  
  type(inventory_type) :: inventory
  PetscInt :: num_cells
  PetscReal :: volume
  
  call ChemSpeciesAllocate(num_cells,inventory%Fe_s, &
                           inventory%preinventory%Fe_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%MgO_s, &
                           inventory%preinventory%MgO_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%RuPl_s, &
                           inventory%preinventory%RubberPlas_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%C6H10O5_s, &
                           inventory%preinventory%Cellulose_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%SO42_minus_aq, &
                           inventory%preinventory%Sulfate_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%NO3_minus_aq, &
                           inventory%preinventory%Nitrate_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%H_ion_aq, &
                           inventory%preinventory%H_ion_in_panel,volume)
  call ChemSpeciesAllocate(num_cells,inventory%FeOH2_s,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%CO2_g,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%N2_g,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%H2S_g,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%FeS_s,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%MgOH2_s,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%Mg5CO34OH24H2_s,0.d0,volume)
  call ChemSpeciesAllocate(num_cells,inventory%MgCO3_s,0.d0,volume)
  
end subroutine InventoryAllocate

! ************************************************************************** !

subroutine ChemSpeciesAllocate(num_cells,chem_species,initial_mass,volume)
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
  PetscReal :: initial_mass     ! mass species in waste panel [kg]
  PetscReal :: volume           ! volume of waste panel [m3]
  
  allocate(chem_species%initial_conc_mol(num_cells))
  allocate(chem_species%current_conc_mol(num_cells))
  allocate(chem_species%current_conc_kg(num_cells))
  allocate(chem_species%inst_rate(num_cells))
  
  chem_species%current_conc_kg(:) = initial_mass / volume             ! [kg/m3]
  chem_species%initial_conc_mol(:) = initial_mass / volume / &
                                     chem_species%molar_mass
  chem_species%current_conc_mol(:) = chem_species%initial_conc_mol(:) ! [mol/m3]
  chem_species%inst_rate(:) = 0.d0        
  chem_species%tot_mass_in_panel = initial_mass
  
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
  PetscInt :: size_of_vec
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
  
  call ProcessAfterRead(this)
  
  if (.not.this%option%restart_flag) then
    call PMWSSOutputHeader(this)
    call PMWSSOutput(this)
  endif
  
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
    call UpdateInventory(cur_waste_panel,dt,this%option)
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  call PMWSSOutput(this)
  
end subroutine PMWSSInitializeTimestep

! *************************************************************************** !

subroutine UpdateInventory(waste_panel,dt,option)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 02/10/2017
  !
  
  use Option_module
 
  implicit none
  
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: dt ! [sec; flow_dt]
  type(option_type) :: option
 
  call UpdateChemSpecies(waste_panel%inventory%Fe_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%FeOH2_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%C6H10O5_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%RuPl_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%H_ion_aq,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%NO3_minus_aq,waste_panel,dt, &
                         option)
  call UpdateChemSpecies(waste_panel%inventory%N2_g,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%CO2_g,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%SO42_minus_aq,waste_panel,dt, &
                         option)
  call UpdateChemSpecies(waste_panel%inventory%H2S_g,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%FeS_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%MgO_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%MgOH2_s,waste_panel,dt,option)
  call UpdateChemSpecies(waste_panel%inventory%Mg5CO34OH24H2_s,waste_panel,dt, &
                         option)
  call UpdateChemSpecies(waste_panel%inventory%MgCO3_s,waste_panel,dt,option)
                                      
 end subroutine UpdateInventory
 
! *************************************************************************** !

subroutine UpdateChemSpecies(chem_species,waste_panel,dt,option)
  !
  ! Updates the waste panel tracked species inventory concentrations.
  !
  ! Author: Jenn Frederick
  ! Date: 2/23/2017
  !
  
  use Option_module
  
  implicit none
  
  type(chem_species_type) :: chem_species
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: dt       ! [sec; flow_dt]
  type(option_type) :: option
  
  PetscInt :: k
  PetscInt :: num_cells
  PetscReal :: local_conc_kg, global_conc_kg
  
  num_cells = waste_panel%region%num_cells
  local_conc_kg = 0.d0  
  do k = 1,num_cells
    ! [mol/m3]
    chem_species%current_conc_mol(k) = &
                 chem_species%current_conc_mol(k) + &   ! [mol/m3]
                 chem_species%inst_rate(k) * &          ! [mol/m3/sec]
                 dt                                     ! [sec]
    ! [kg/m3]
    chem_species%current_conc_kg(k) = &
                 chem_species%current_conc_mol(k) * &   ! [mol/m3]
                 chem_species%molar_mass                ! [kg/mol]
    ! [kg/m3]             
    local_conc_kg = local_conc_kg + &
                (chem_species%current_conc_kg(k)*waste_panel%scaling_factor(k))
  enddo
  ! [kg]
  call CalcParallelSUM(option,waste_panel,local_conc_kg,global_conc_kg)
  chem_species%tot_mass_in_panel = global_conc_kg * &   ! [kg/m3]
                                   waste_panel%volume   ! [m3]
                                   
end subroutine UpdateChemSpecies

! *************************************************************************** !

 subroutine PMWSSSolve(this,time,ierr)
  ! 
  ! Calculates reaction rates, gas generation rate, and brine generation rate.
  ! Sets the fluid and energy source terms.
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
  PetscInt :: i, j
  PetscInt :: cell_id
  PetscInt :: num_cells
  PetscReal :: conc_ratio
  ! brine/gas generation variable
  PetscReal :: water_saturation
  PetscReal :: s_eff
  PetscReal :: rxnrate_corrosion
  PetscReal :: rxnrate_biodeg_nitrate
  PetscReal :: rxnrate_biodeg_sulfate
  PetscReal :: rxnrate_FeS_Fe
  PetscReal :: rxnrate_FeS_FeOH2
  PetscReal :: rxnrate_mgoh2
  PetscReal :: rxnrate_hydromag
  PetscReal :: rxnrate_hymagcon
  PetscReal :: rxnrate1, rxnrate2, rxnrate3, rxnrate4
  ! enthalpy calculation variables
  PetscReal :: temperature
  PetscReal :: pressure_liq
  PetscReal :: pressure_gas
  PetscReal :: H_liq
  PetscReal :: H_gas
  PetscReal :: U_gas
  PetscReal :: gas_energy
  PetscReal :: brine_energy
  
  !return
  
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
    !-----zero-out-reaction-rates---------------------------------------------
    water_saturation = 0.d0
    s_eff = 0.d0
    rxnrate_corrosion = 0.d0
    rxnrate_biodeg_nitrate = 0.d0
    rxnrate_biodeg_sulfate = 0.d0
    rxnrate_FeS_Fe = 0.d0
    rxnrate_FeS_FeOH2 = 0.d0
    rxnrate_mgoh2 = 0.d0
    rxnrate_hydromag = 0.d0
    rxnrate_hymagcon = 0.d0

    do i = 1,num_cells
    !-----effective-brine-saturation------------------------------------------
      water_saturation = &
        gen_auxvar(ZERO_INTEGER,grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
        sat(option%liquid_phase)
      s_eff = water_saturation - this%smin + this%satwick*(1.d0 - &
        exp(200.d0*this%alpharxn*(max((water_saturation-this%smin),0.d0))**2.d0))
    !-----anoxic-iron-corrosion-[mol-Fe/m3/sec]-------------------------------
      rxnrate_corrosion = (cur_waste_panel%inundated_corrosion_rate*s_eff) + &
                          (cur_waste_panel%humid_corrosion_rate*(1.d0-s_eff))
      call SmoothRxnrate(rxnrate_corrosion,i,cur_waste_panel%inventory%Fe_s, &
                         this%alpharxn) 
      call TaperRxnrate(rxnrate_corrosion,i,cur_waste_panel%inventory%Fe_s)
    !-----biodegradation-[mol-cell/m3/sec]------------------------------------
      rxnrate_biodeg_nitrate = (cur_waste_panel%inundated_biodeg_rate*s_eff) + &
                               (cur_waste_panel%humid_biodeg_rate*(1.d0-s_eff))
      rxnrate_biodeg_sulfate = (cur_waste_panel%inundated_biodeg_rate*s_eff) + &
                               (cur_waste_panel%humid_biodeg_rate*(1.d0-s_eff))
      call SmoothRxnrate(rxnrate_biodeg_nitrate,i, &
                         cur_waste_panel%inventory%C6H10O5_s,this%alpharxn)
      call SmoothRxnrate(rxnrate_biodeg_sulfate,i, &
                         cur_waste_panel%inventory%C6H10O5_s,this%alpharxn)
      call TaperRxnrate(rxnrate_biodeg_nitrate,i, &
                        cur_waste_panel%inventory%C6H10O5_s, &
                        cur_waste_panel%inventory%NO3_minus_aq)
      call TaperRxnrate(rxnrate_biodeg_sulfate,i, &
                        cur_waste_panel%inventory%C6H10O5_s, &
                        cur_waste_panel%inventory%SO42_minus_aq)
    !-----iron-sulfidation-[mol-H2S/m3/sec]-----------------------------------
      rxnrate_FeS_Fe = rxnrate_biodeg_sulfate*cur_waste_panel%RXH2S_factor
      rxnrate_FeS_FeOH2 = rxnrate_biodeg_sulfate*cur_waste_panel%RXH2S_factor
      call SmoothRxnrate(rxnrate_FeS_Fe,i,cur_waste_panel%inventory%Fe_s, &
                         this%alpharxn) 
      call SmoothRxnrate(rxnrate_FeS_FeOH2,i,cur_waste_panel%inventory%Fe_s, &
                         this%alpharxn)
      call TaperRxnrate(rxnrate_FeS_Fe,i, &
                        cur_waste_panel%inventory%Fe_s, &
                        cur_waste_panel%inventory%H2S_g)
      call TaperRxnrate(rxnrate_FeS_FeOH2,i, &
                        cur_waste_panel%inventory%FeOH2_s, &
                        cur_waste_panel%inventory%H2S_g)
    !-----MgO-hydration-[mol-MgO/m3/sec]--------------------------------------
      rxnrate_mgoh2 = (cur_waste_panel%inundated_brucite_rate*s_eff) + &
                      ((cur_waste_panel%humid_brucite_rate)*(1.d0-s_eff))
      call SmoothRxnrate(rxnrate_mgoh2,i,cur_waste_panel%inventory%MgO_s, &
                         this%alpharxn)
      call TaperRxnrate(rxnrate_mgoh2,i,cur_waste_panel%inventory%MgO_s)
    !-----hydromagnesite-[mol/m3-bulk/sec]------------------------------------
      rxnrate_hydromag = max(rxnrate_biodeg_nitrate,rxnrate_biodeg_sulfate) * &
                         this%RXCO2_factor
      call SmoothRxnrate(rxnrate_hydromag,i,cur_waste_panel%inventory%MgO_s, &
                         this%alpharxn)
      call TaperRxnrate(rxnrate_hydromag,i,cur_waste_panel%inventory%MgOH2_s)
    !-----hydromagnesite-conversion-[mol/m3-bulk/sec]-------------------------
      rxnrate_hymagcon = this%hymagcon_rate* &
                  cur_waste_panel%inventory%Mg5CO34OH24H2_s%current_conc_kg(i)
      call SmoothRxnrate(rxnrate_hymagcon,i,cur_waste_panel%inventory%MgO_s, &
                         this%alpharxn)
      call TaperRxnrate(rxnrate_hydromag,i, &
                        cur_waste_panel%inventory%Mg5CO34OH24H2_s)
    !-----tracked-species-[mol-species/m3-bulk/sec]---------------------------
      cur_waste_panel%inventory%FeOH2_s%inst_rate(i) = &
          1.d0*rxnrate_corrosion + (-1.d0*rxnrate_FeS_FeOH2)
      cur_waste_panel%inventory%Fe_s%inst_rate(i) = &
          (-1.d0*rxnrate_corrosion) + (-1.d0*rxnrate_FeS_Fe) 
      cur_waste_panel%inventory%FeS_s%inst_rate(i) = &
          1.d0*rxnrate_FeS_Fe + 1.d0*rxnrate_FeS_FeOH2
      cur_waste_panel%inventory%C6H10O5_s%inst_rate(i) = &
          (-1.d0*rxnrate_biodeg_nitrate) + (-1.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%RuPl_s%inst_rate(i) = &
          (-1.d0*rxnrate_biodeg_nitrate) + (-1.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%H_ion_aq%inst_rate(i) = &
          (-4.8d0*rxnrate_biodeg_nitrate) + (-6.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%NO3_minus_aq%inst_rate(i) = &
          (-4.8d0*rxnrate_biodeg_nitrate)
      cur_waste_panel%inventory%CO2_g%inst_rate(i) = &
          6.d0*rxnrate_biodeg_sulfate + 6.d0*rxnrate_biodeg_nitrate + &
          (-4.d0*rxnrate_hydromag)
      cur_waste_panel%inventory%N2_g%inst_rate(i) = &
          2.4d0*rxnrate_biodeg_nitrate
      cur_waste_panel%inventory%SO42_minus_aq%inst_rate(i) = &
          (-3.d0*rxnrate_biodeg_sulfate)
      cur_waste_panel%inventory%H2S_g%inst_rate(i) = &
          3.d0*rxnrate_biodeg_sulfate + (-1.d0*rxnrate_FeS_Fe) + &
          (-1.d0*rxnrate_FeS_FeOH2)
      cur_waste_panel%inventory%MgO_s%inst_rate(i) = &
          (-1.d0*rxnrate_mgoh2)
      cur_waste_panel%inventory%MgOH2_s%inst_rate(i) = &
          1.d0*rxnrate_mgoh2 + (-5.d0*rxnrate_hydromag) + 1.d0*rxnrate_hymagcon
      cur_waste_panel%inventory%Mg5CO34OH24H2_s%inst_rate(i) = &
          1.d0*rxnrate_hydromag + (-1.d0*rxnrate_hymagcon) 
      cur_waste_panel%inventory%MgCO3_s%inst_rate(i) = &
          4.d0*rxnrate_hymagcon 
    !-----gas-generation-[mol-H2/m3-bulk/sec]---------------------------------
      cur_waste_panel%gas_generation_rate(i) = &
          1.d0*rxnrate_corrosion + 1.d0*rxnrate_FeS_Fe + &
          2.4d0*rxnrate_biodeg_nitrate
    !-----brine-generation-[mol-H2O/m3-bulk/sec]------------------------------
      cur_waste_panel%brine_generation_rate(i) = &
          (-2.d0*rxnrate_corrosion) + 7.4d0*rxnrate_biodeg_nitrate + &
          5.d0*rxnrate_biodeg_sulfate + 2.d0*rxnrate_FeS_FeOH2 + &
          (-1.d0*rxnrate_mgoh2) + 4.d0*rxnrate_hymagcon
    !------source-term-calculation--------------------------------------------
      cell_id = grid%nL2G(cur_waste_panel%region%cell_ids(i))
      j = j + 1
      ! liquid source term [kmol/sec]
      vec_p(j) = cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
                 material_auxvars(cell_id)%volume / &          ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
      j = j + 1
      ! gas source term [kmol/sec]
      vec_p(j) = cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
                 material_auxvars(cell_id)%volume / &          ! [m3-bulk]
                 1.d3                                          ! [mol -> kmol]
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
              material_auxvars(cell_id)%volume * &          ! [m3-bulk] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
        case(LIQUID_STATE) !---------------------------------------------------
          pressure_liq = gen_auxvar(ZERO_INTEGER, &
              grid%nL2G(cur_waste_panel%region%cell_ids(i)))% &
              pres(option%liquid_phase)
          call EOSWaterEnthalpy(temperature,pressure_liq,H_liq,ierr)
          brine_energy = &
              cur_waste_panel%brine_generation_rate(i) * &  ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3-bulk] 
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
              material_auxvars(cell_id)%volume * &          ! [m3-bulk] 
              H_liq * 1.d-3 * 1.d-6                         ! [MJ/mol]
          gas_energy = &
              cur_waste_panel%gas_generation_rate(i) * &    ! [mol/m3/sec]
              material_auxvars(cell_id)%volume * &          ! [m3-bulk] 
              H_gas * 1.d-3 * 1.d-6                         ! [MJ/mol]
      end select
      vec_p(j) = brine_energy + gas_energy
    enddo
    !-------------------------------------------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWSSSolve

! ************************************************************************** !

subroutine SmoothRxnrate(rxnrate,cell_num,limiting_species,alpharxn)
  !
  ! Smooths the reaction rate near the point where the reaction runs out of a
  ! limiting relevant reactant/species. This implements Eq. 158 in the BRAGFLO
  ! User's Manual.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species
  PetscReal :: alpharxn
  
  PetscReal :: conc_ratio
  
  conc_ratio = ( limiting_species%initial_conc_mol(cell_num) / &
                 limiting_species%current_conc_mol(cell_num) ) 
  ! K_smoothed = K * (1.0 - exp(A*C/Ci)  BRAGFLO User's Manual Eq. 158
  rxnrate = rxnrate * (1.d0 - exp(alpharxn*conc_ratio))
  
end subroutine SmoothRxnrate

! ************************************************************************** !

subroutine TaperRxnrate1(rxnrate,cell_num,limiting_species1)
  !
  ! Tapers the reaction rate if the reaction runs out of a single
  ! limiting relevant reactant/species. The limiting reactant/species is
  ! chosen according to the equations in the BRAGFLO User's Manual, 
  ! Section 14.13.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  
  if (limiting_species1%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate = 0.d0
  else
    rxnrate = rxnrate
  endif
  
end subroutine TaperRxnrate1

! ************************************************************************** !

subroutine TaperRxnrate2(rxnrate,cell_num,limiting_species1,limiting_species2)
  !
  ! Tapers the reaction rate if the reaction runs out of one of two possible
  ! limiting relevant reactants/species. The limiting reactants/species are
  ! chosen according to the equations in the BRAGFLO User's Manual, 
  ! Section 14.13.
  !
  ! Author: Jennifer Frederick
  ! Date: 03/27/3017
  !
  
  implicit none
  
  PetscReal :: rxnrate
  PetscInt :: cell_num
  type(chem_species_type) :: limiting_species1
  type(chem_species_type) :: limiting_species2
  
  PetscReal :: rxnrate1
  PetscReal :: rxnrate2
  
  if (limiting_species1%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate1 = 0.d0
  else
    rxnrate1 = rxnrate
  endif
  if (limiting_species2%current_conc_mol(cell_num) <= 0.d0) then
    rxnrate2 = 0.d0
  else
    rxnrate2 = rxnrate
  endif
  rxnrate = min(rxnrate1,rxnrate2)
  
end subroutine TaperRxnrate2

! ************************************************************************** !

subroutine PMWSSFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/22/2017
  !

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
end subroutine PMWSSFinalizeTimestep

! ************************************************************************** !

subroutine CalcParallelSUM(option,waste_panel,local_val,global_sum)
  ! 
  ! Calculates global sum for a MPI_DOUBLE_PRECISION number over a
  ! waste panel region. This function uses only MPI_Send and MPI_Recv functions
  ! and does not need a communicator object. It reduces communication to the
  ! processes that are in the waste panel's rank_list object rather than using
  ! a call to MPI_Allreduce.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/23/17
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(srcsink_panel_type) :: waste_panel
  PetscReal :: local_val
  PetscReal :: global_sum

  PetscReal, pointer :: temp_array(:)
  PetscInt :: num_ranks
  PetscInt :: m
  PetscInt :: TAG
  PetscErrorCode :: ierr
  
  num_ranks = size(waste_panel%rank_list)
  allocate(temp_array(num_ranks))
  temp_array = 0.d0
  TAG = 0
  
  if (num_ranks > 1) then
  !------------------------------------------
    if (option%myrank .ne. waste_panel%rank_list(1)) then
      call MPI_Send(local_val,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                    waste_panel%rank_list(1),TAG,option%mycomm,ierr)
    else
      temp_array(1) = local_val
      do m = 2,num_ranks
        call MPI_Recv(local_val,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      waste_panel%rank_list(m),TAG,option%mycomm, &
                      MPI_STATUS_IGNORE,ierr)
        temp_array(m) = local_val
      enddo
      global_sum = sum(temp_array)
    endif
    if (option%myrank == waste_panel%rank_list(1)) then
      do m = 2,num_ranks
        call MPI_Send(global_sum,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      waste_panel%rank_list(m),TAG,option%mycomm,ierr)
      enddo
    else
      call MPI_Recv(global_sum,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                    waste_panel%rank_list(1),TAG,option%mycomm, &
                    MPI_STATUS_IGNORE,ierr)
    endif             
  !------------------------------------------        
  else 
    global_sum = local_val
  endif
  
  deallocate(temp_array)

end subroutine CalcParallelSUM

! *************************************************************************** !

 subroutine PMWSSOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/27/2017
  !
  
  use Option_module
  use Output_Aux_module
  
  implicit none

  class(pm_wipp_srcsink_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(srcsink_panel_type), pointer :: cur_waste_panel
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  PetscReal :: local_gas_rate, global_gas_rate
  PetscReal :: local_brine_rate, global_brine_rate
  
  if (.not.associated(this%waste_panel_list)) return
  
100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option
  
  fid = 88
  filename = PMWSSOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")
       
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    ! pre-calculations
    local_gas_rate = 0.d0 
    local_brine_rate = 0.d0
    do i = 1,cur_waste_panel%region%num_cells            
      local_gas_rate = local_gas_rate + &    ! [mol/m3-bulk/sec]
                   (cur_waste_panel%gas_generation_rate(i) * &
                    cur_waste_panel%scaling_factor(i))
      local_brine_rate = local_brine_rate + &    ! [mol/m3-bulk/sec]
                   (cur_waste_panel%brine_generation_rate(i) * &
                    cur_waste_panel%scaling_factor(i))
    enddo
    call CalcParallelSUM(option,cur_waste_panel,local_gas_rate,global_gas_rate)
    call CalcParallelSUM(option,cur_waste_panel,local_brine_rate, &
                         global_brine_rate)
    write(fid,101,advance="no") cur_waste_panel%id
    write(fid,100,advance="no") &
      cur_waste_panel%inventory%Fe_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%MgO_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%C6H10O5_s%tot_mass_in_panel, &
      cur_waste_panel%inventory%RuPl_s%tot_mass_in_panel, &
      global_gas_rate, &
      global_brine_rate
  !----------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  close(fid)
  
end subroutine PMWSSOutput

! ************************************************************************** !

function PMWSSOutputFilename(option)
  ! 
  ! Generates filename for waste panel output.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 03/27/17

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: PMWSSOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMWSSOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-' // trim(adjustl(word)) // '.pnl'
  
end function PMWSSOutputFilename  

! ************************************************************************** !

subroutine PMWSSOutputHeader(this)
  ! 
  ! Writes header for waste panel output file.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 03/27/17

  use Output_Aux_module
    
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  class(srcsink_panel_type), pointer :: cur_waste_panel
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units_string 
  character(len=MAXWORDLENGTH) :: variable_string
  character(len=MAXWORDLENGTH) :: cell_string
  PetscInt :: fid
  PetscInt :: icolumn
  
  if (.not.associated(this%waste_panel_list)) return
  
  output_option => this%realization%output_option
  
  fid = 88
  filename = PMWSSOutputFilename(this%option)
  open(unit=fid,file=filename,action="write",status="replace")
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'
  
  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    cell_string = trim(cur_waste_panel%region_name)
    variable_string = 'WP ID#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Fe(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'MgO(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Cellulosics(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Rubber/Plastics(s) mass'
    units_string = 'kg'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Avg. gas gen. rate'
    units_string = 'mol/m3-bulk/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'Avg. brine gen. rate'
    units_string = 'mol/m3-bulk/sec'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
  !----------------------------------------
    cur_waste_panel => cur_waste_panel%next
  enddo
  
  close(fid)
  
end subroutine PMWSSOutputHeader

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
  ! Strips and destroys the WIPP source/sink process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !

  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  call PMWSSStrip(this)
  
end subroutine PMWSSDestroy

! ************************************************************************** !

subroutine PMWSSStrip(this)
  ! 
  ! Strips the WIPP source/sink process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/31/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_wipp_srcsink_type) :: this
  
  type(srcsink_panel_type), pointer :: cur_waste_panel, prev_waste_panel
  type(pre_inventory_type), pointer :: cur_preinventory, prev_preinventory

  cur_waste_panel => this%waste_panel_list
  do
    if (.not.associated(cur_waste_panel)) exit
    prev_waste_panel => cur_waste_panel
    cur_waste_panel => cur_waste_panel%next
    call InventoryDeallocate(prev_waste_panel%inventory)
    call DeallocateArray(prev_waste_panel%scaling_factor)
    call DeallocateArray(prev_waste_panel%gas_generation_rate)
    call DeallocateArray(prev_waste_panel%brine_generation_rate)
    call DeallocateArray(prev_waste_panel%rank_list)
    deallocate(prev_waste_panel)
    nullify(prev_waste_panel)
  enddo
  nullify(this%waste_panel_list)
  
  cur_preinventory => this%pre_inventory_list
  do
    if (.not.associated(cur_preinventory)) exit
    prev_preinventory => cur_preinventory
    cur_preinventory => cur_preinventory%next
    deallocate(prev_preinventory)
    nullify(prev_preinventory)
  enddo
  nullify(this%pre_inventory_list)

end subroutine PMWSSStrip

! ************************************************************************** !

subroutine InventoryDeallocate(inventory)
  ! 
  ! Deallocates the inventory's chemical species.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/01/2017
  !
  
  implicit none
  
  type(inventory_type) :: inventory
  
  call ChemSpeciesDeallocate(inventory%Fe_s)
  call ChemSpeciesDeallocate(inventory%FeOH2_s)
  call ChemSpeciesDeallocate(inventory%C6H10O5_s)
  call ChemSpeciesDeallocate(inventory%RuPl_s)
  call ChemSpeciesDeallocate(inventory%H_ion_aq)
  call ChemSpeciesDeallocate(inventory%NO3_minus_aq)
  call ChemSpeciesDeallocate(inventory%CO2_g)
  call ChemSpeciesDeallocate(inventory%N2_g)
  call ChemSpeciesDeallocate(inventory%SO42_minus_aq)
  call ChemSpeciesDeallocate(inventory%H2S_g)
  call ChemSpeciesDeallocate(inventory%FeS_s)
  call ChemSpeciesDeallocate(inventory%MgO_s)
  call ChemSpeciesDeallocate(inventory%MgOH2_s)
  call ChemSpeciesDeallocate(inventory%Mg5CO34OH24H2_s)
  call ChemSpeciesDeallocate(inventory%MgCO3_s)
  nullify(inventory%preinventory)
  
end subroutine InventoryDeallocate

! ************************************************************************** !

subroutine ChemSpeciesDeallocate(chem_species)
  ! 
  ! Deallocates the inventory's chemical species.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/01/2017
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(chem_species_type) :: chem_species
  
  call DeallocateArray(chem_species%initial_conc_mol)
  call DeallocateArray(chem_species%current_conc_mol)
  call DeallocateArray(chem_species%current_conc_kg)
  call DeallocateArray(chem_species%inst_rate)
  
end subroutine ChemSpeciesDeallocate

! *************************************************************************** !

end module PM_WIPP_SrcSink_class
