module PM_UFD_Decay_class

  use PM_Base_class
  use Realization_class
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  
  type, public, extends(pm_base_type) :: pm_ufd_decay_type
    class(realization_type), pointer :: realization
!    character(len=MAXWORDLENGTH) :: species
    PetscInt, pointer :: component(:)
    PetscInt, pointer :: isotope(:,:)
    PetscReal, pointer :: decay_rate(:,:)
    PetscInt, pointer :: daughter(:,:,:)
    PetscReal, pointer :: solubility(:)
    PetscReal, pointer :: molar_volume(:)
    PetscReal, pointer :: Kd(:)
    PetscReal, pointer :: mole_fraction(:,:,:)
    PetscReal, pointer :: mass0(:,:)
    PetscInt :: num_components = 5
    PetscInt :: max_num_isotopes = 3
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Setup => PMUFDDecayInit
    procedure, public :: Read => PMUFDDecayRead
!    procedure, public :: SetupSolvers => PMUFDDecaySetupSolvers
    procedure, public :: PMUFDDecaySetRealization
    procedure, public :: InitializeRun => PMUFDDecayInitializeRun
!!    procedure, public :: FinalizeRun => PMUFDDecayFinalizeRun
    procedure, public :: InitializeTimestep => PMUFDDecayInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDDecayFinalizeTimestep
!    procedure, public :: PreSolve => PMUFDDecayPreSolve
    procedure, public :: Solve => PMUFDDecaySolve
!    procedure, public :: PostSolve => PMUFDDecayPostSolve
!    procedure, public :: AcceptSolution => PMUFDDecayAcceptSolution
!    procedure, public :: TimeCut => PMUFDDecayTimeCut
!    procedure, public :: UpdateSolution => PMUFDDecayUpdateSolution
!    procedure, public :: UpdateAuxvars => PMUFDDecayUpdateAuxvars
!    procedure, public :: Checkpoint => PMUFDDecayCheckpoint    
!    procedure, public :: Restart => PMUFDDecayRestart  
    procedure, public :: Destroy => PMUFDDecayDestroy
  end type pm_ufd_decay_type
  
  public :: PMUFDDecayCreate, &
            PMUFDDecayInit !, &
!            PMUFDDecaySetupSolvers, &
!            PMUFDDecayInitializeTimestepA, &
!            PMUFDDecayInitializeTimestepB, &
!            PMUFDDecayInitializeRun, &
!            PMUFDDecayUpdateSolution, &
!            PMUFDDecayUpdatePropertiesNI, &
!            PMUFDDecayTimeCut, &
!            PMUFDDecayDestroy
  
contains

! ************************************************************************** !

function PMUFDDecayCreate()
  ! 
  ! Creates the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type), pointer :: PMUFDDecayCreate
  
  allocate(PMUFDDecayCreate)
  nullify(PMUFDDecayCreate%realization)

  call PMBaseInit(PMUFDDecayCreate)

end function PMUFDDecayCreate

! ************************************************************************** !

subroutine PMUFDDecayRead(this,input)
  ! 
  ! Reads input file parameters associated with the ufd decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  type(input_type) :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => this%option
  
  error_string = 'UFD Decay'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMUFDDecayRead

! ************************************************************************** !

subroutine PMUFDDecayInit(this)
  ! 
  ! Initializes variables associated with the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Reaction_Aux_module
  use Option_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  
  PetscInt, parameter :: max_daughters = 1
  PetscInt, parameter :: max_num_isotopes = 1
  
  grid => this%realization%patch%grid
  option => this%realization%option
  reaction => this%realization%reaction

#if 0
  this%num_components = 5
  allocate(this%num_isotopes(this%num_components))
  this%num_isotopes = 1
  max_num_isotopes = maxval(this%num_isotopes)
  allocate(this%isotope(max_num_isotopes,this%num_components))
  this%isotopes(1,:) = [1, 2, 3, 4, 5]
  allocate(this%decay_rate(max_num_isotopes,this%num_components))
  this%decay_rate(1,:) = [1.29d-15, 5.08d-11, 1.03d-14, 1.38d-13, 2.78d-12]
  allocate(this%daughters(0:max_daughter,this%max_num_isotopes))
  this%daughter = 0
  this%daughter(0,2:4) = 1
  this%daughter(1,2) = 3
  this%daughter(1,3) = 4
  this%daughter(1,4) = 5
  allocate(this%solubility(this%num_components))
  this%solubility = 0.d0
#endif

end subroutine PMUFDDecayInit

! ************************************************************************** !

subroutine PMUFDDecaySetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Realization_class

  implicit none
  
  class(pm_ufd_decay_type) :: this
  class(realization_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDDecaySetRealization

! ************************************************************************** !

recursive subroutine PMUFDDecayInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none

  class(pm_ufd_decay_type) :: this
  

end subroutine PMUFDDecayInitializeRun

! ************************************************************************** !

subroutine PMUFDDecayInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Global_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," USED FUEL DISPOSITION DECAY MODEL ",43("="))')
  endif

end subroutine PMUFDDecayInitializeTimestep

! ************************************************************************** !

subroutine PMUFDDecayPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Grid_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
end subroutine PMUFDDecayPreSolve

! ************************************************************************** !

subroutine PMUFDDecaySolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  !
  use Option_module
  use Reaction_Aux_module
  use Patch_module
  use Grid_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none

  class(pm_ufd_decay_type) :: this
  
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: icomp, i, idaughter
  PetscReal :: Aq0, Sorb0, Ppt0
  PetscReal :: Aq1, Sorb1, Ppt1
  PetscReal :: Mass_Aq0, Mass_Sorb0, Mass_Ppt0
  PetscReal :: Mass_Aq1, Mass_Sorb1, Mass_Ppt1
  PetscReal :: Aq_tot1
  PetscReal :: Mass_tot0(this%num_components), Mass_tot1(this%num_components)
  PetscReal :: Delta_Mass_tot(this%num_components)
  PetscReal :: kd_kgw_m3b
  PetscReal :: vol, den_w_kg, por, sat, vps, dt
  PetscReal :: mol_fraction(this%max_num_isotopes)
  PetscReal :: mass(this%max_num_isotopes,this%num_components,this%realization%patch%grid%nlmax)
  PetscReal, pointer :: vec_loc_p(:)
  
  ierr = 0
  
  option => this%realization%option
  reaction => this%realization%reaction
  patch => this%realization%patch
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  dt = option%dt

#if 0  
  do i = 1, this%num_components
    do j = 1, this%num_isotopes(i)
      call VecGetArrayF90(field%work_loc, vec_loc_p, ierr);CHKERRQ(ierr)
      vec_loc_p(:) = this%mol_fraction(j,i,:)
      call VecRestoreArrayF90(field%work_loc, vec_loc_p, ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,field%work_loc, &
                                       field%work_loc,ONEDOF)   
      call VecGetArrayReadF90(field%work_loc, vec_loc_p, ierr);CHKERRQ(ierr)
      this%mol_fraction(j,i,:) = vec_loc_p(:)
      call VecRestoreArrayReadF90(field%work_loc, vec_loc_p, ierr);CHKERRQ(ierr)
    enddo
  enddo
  
  mass = 0.d0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    vol = material_auxvars(ghosted_id)%volume
    por = material_auxvars(ghosted_id)%porosity
    sat = global_auxvars(ghosted_id)%sat(1)
    vps = vol * por * sat ! m^3 water
    do i = 1, this%num_components
      icomp = this%component(i)
      total = rt_auxvars(ghosted_id)%total(icomp,1)* &
              vps*1.d3 ! mol/L * m^3 water * 1000 L /m^3 = mol
      do j = 1, this%num_isotopes(i)
        mass(j,i,local_id) = total*this%mole_fraction(j,i,ghosted_id)
      enddo
    enddo
  enddo  
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
!      if (patch%imat(ghosted_id_up) <= 0 .or.  &
!          patch%imat(ghosted_id_dn) <= 0) cycle
      do i = 1, this%num_components
        icomp = this%component(i)
        flux = patch%internal_tran_fluxes(icomp,iconn) * dt
        do j = 1, this%num_isotopes(i)
          iiso = this%isotope(j,i)
          if (flux > 0.d0) then
            mass(j,i,ghosted_id_up) = mass(j,i,local_id_up) - flux*this%mol_fraction(j,i,ghosted_id_up)
            mass(j,i,ghosted_id_dn) = mass(j,i,local_id_dn) + flux*this%mol_fraction(j,i,ghosted_id_up)
          else
            mass(j,i,ghosted_id_up) = mass(j,i,local_id_up) + flux*this%mol_fraction(j,i,ghosted_id_dn)
            mass(j,i,ghosted_id_dn) = mass(j,i,local_id_dn) - flux*this%mol_fraction(j,i,ghosted_id_dn)
          endif
        enddo
      enddo
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
    
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
!      if (patch%imat(ghosted_id) <= 0) cycle
      do i = 1, this%num_components
        icomp = this%component(i)
        flux = patch%boundary_tran_fluxes(icomp,iconn) * dt
        do j = 1, this%num_isotopes(i)
          iiso = this%isotope(j,i)
          mass(j,i,ghosted_id) = mass(j,i,ghosted_id) + flux*this%mol_fraction(j,i,ghosted_id)
        enddo
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    do i = 1, this%num_components
      do j = 1, this%num_isotopes(i)
        this%mol_fraction(j,i,ghosted_id) = mass0(i,ghosted_id)* &
                                            this%mol_fraction(j,i,ghosted_id) / &
                                            (mass(
      enddo  
    enddo
  enddo
  
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    vol = material_auxvars(ghosted_id)%volume
    den_w_kg = global_auxvars(ghosted_id)%den_kg(1)
    por = material_auxvars(ghosted_id)%porosity
    sat = global_auxvars(ghosted_id)%sat(1)
    vps = vol * por * sat ! m^3 water
    do icomp = 1, this%num_components
      do i = 1, this%num_isotopes(icomp)
        iiso = this%isotope(i,icomp)
        Aq0 = rt_auxvars(ghosted_id)%total(iiso,1) ! mol/L
        Sorb0 = rt_auxvars(ghosted_id)%total_sorb_eq(iiso) ! mol/m^3 bulk
        Ppt0 = rt_auxvars(ghosted_id)%mnrl_volfrac(iiso)
        Mass_Aq0 = Aq0*vps*1.d3 ! mol/L * m^3 water * 1000 L /m^3 = mol
        Mass_Sorb0 = Sorb0 * vol ! mol/m^3 bulk * m^3 bulk = mol
        Mass_Ppt0 = Ppt0 * vol / &  ! m^3 mnrl/m^3 bulk * m^3 bulk / (m^3 mnrl/mol mnrl) = mol
                    this%molar_volume(icomp)
        Mass_tot0 = Mass_Aq0 + Mass_Sorb0 + Mass_Ppt0
        Mass_tot1(iiso) = Mass_tot0 * exp(this%decay_rate(iiso)*dt)
        Delta_Mass_tot(iiso) = Mass_tot1(iiso) - Mass_tot0
      enddo
    enddo 
    
    do icomp = 1, this%num_components
      do i = 1, this%num_isotopes(icomp)
        iiso = this%isotope(i,icomp)
        do i = 1, this%daughter(0,iiso)
          idaughter = this%daughter(i,iiso)
          ! Delta is always negative for parent, positive for child
          Mass_tot1(idaughter) = Mass_tot1(idaughter) - Delta_Mass_tot(iiso)
        enddo
      enddo
    enddo
    
    do icomp = 1, this%num_components
      ! calculate mole fractions
      Mass_total = 0.d0
      mol_fraction = 0.d0
      do i = 1, this%num_isotopes(icomp)
        iiso = this%isotope(i,icomp)
        Mass_total = Mass_total + Mass_tot1(iiso)
      enddo
      do i = 1, this%num_isotopes(icomp)
        mol_fraction(i) = Mass_tot1(iiso) / Mass_total
      enddo
      ! split mass between phases
      kd_kgw_m3b = this%KD(icomp)
      Aq_tot1 = Mass_total /  (1.d0+kd_kgw_m3b/(den_w_kg*por*sat)) / &
              (vps*1.d3)
      if (Aq_tot1 > this%solubility(icomp)) then
        Aq1 = this%solubility(icomp)
      else
        Aq1 = Aq_tot1
      endif
      Sorb1 = Aq_tot1 / den_w_kg * 1.d3 * kd_kgw_m3b
      Mass_Aq1 = Aq_tot1*vps*1.d3
      Mass_Sorb1 = Sorb0 * vol
      Mass_Ppt1 = max(Mass_total - Mass_Aq0 - Mass_Sorb0,0.d0)
      Ppt1 = Mass_Ppt1 * this%molar_volume(icomp) / vol
      ! store mass in data structures
      do i = 1, this%num_isotopes(icomp)
        iiso = this%isotope(i,icomp)
        rt_auxvars(ghosted_id)%total(iiso,1) = Aq1 * mole_fraction(i)
        rt_auxvars(ghosted_id)%free_ion(iiso) = Aq1 / den_w_kg * 1.d3 * mole_fraction(i)
        rt_auxvars(ghosted_id)%total_sorb_eq(iiso) = Sorb1 * mole_fraction(i)
        rt_auxvars(ghosted_id)%mnrl_volfrac(iiso) = Ppt1 * mole_fraction(i)
      enddo
    enddo      
  enddo
#endif
  
end subroutine PMUFDDecaySolve

! ************************************************************************** !

subroutine PMUFDDecayPostSolve(this)
  ! 
  ! PMUFDDecayUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  ! 
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayPostSolve

! ************************************************************************** !

function PMUFDDecayAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscBool :: PMUFDDecayAcceptSolution
  
  ! do nothing
  PMUFDDecayAcceptSolution = PETSC_TRUE
  
end function PMUFDDecayAcceptSolution

! ************************************************************************** !

subroutine PMUFDDecayUpdatePropertiesTS(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayUpdatePropertiesTS

! ************************************************************************** !

subroutine PMUFDDecayTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMUFDDecayTimeCut

! ************************************************************************** !

subroutine PMUFDDecayFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayFinalizeTimestep

! ************************************************************************** !

subroutine PMUFDDecayUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMUFDDecayUpdateSolution  

! ************************************************************************** !

subroutine PMUFDDecayUpdateAuxvars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this

  this%option%io_buffer = 'PMUFDDecayUpdateAuxvars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMUFDDecayUpdateAuxvars   

! ************************************************************************** !

subroutine PMUFDDecayCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayCheckpoint

! ************************************************************************** !

subroutine PMUFDDecayRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayRestart

! ************************************************************************** !

recursive subroutine PMUFDDecayFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMUFDDecayFinalizeRun

! ************************************************************************** !

subroutine PMUFDDecayDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayDestroy
  
end module PM_UFD_Decay_class
