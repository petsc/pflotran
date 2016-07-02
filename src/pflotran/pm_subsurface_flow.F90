module PM_Subsurface_Flow_class

  use PM_Base_class
!geh: using Init_Subsurface_module here fails with gfortran (internal compiler error)
!  use Init_Subsurface_module
  use Realization_Subsurface_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"

  type, public, extends(pm_base_type) :: pm_subsurface_flow_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    PetscBool :: transient_permeability
    PetscBool :: store_porosity_for_ts_cut
    PetscBool :: store_porosity_for_transport
    PetscBool :: check_post_convergence
    ! these govern the size of subsequent time steps
    PetscReal :: max_pressure_change
    PetscReal :: max_temperature_change
    PetscReal :: max_saturation_change
    PetscReal :: max_xmol_change
    PetscReal :: pressure_change_governor
    PetscReal :: temperature_change_governor
    PetscReal :: saturation_change_governor
    PetscReal :: xmol_change_governor
    PetscReal :: cfl_governor
    ! these limit (truncate) the maximum change in a Newton iteration
    ! truncation occurs within PMXXXCheckUpdatePre
    PetscReal :: pressure_dampening_factor
    PetscReal :: saturation_change_limit
    PetscReal :: pressure_change_limit
    PetscReal :: temperature_change_limit
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Setup => PMSubsurfaceFlowSetup
    procedure, public :: SetupSolvers => PMSubsurfaceFlowSetupSolvers
    procedure, public :: PMSubsurfaceFlowSetRealization
    procedure, public :: InitializeRun => PMSubsurfaceFlowInitializeRun
!    procedure, public :: FinalizeRun => PMSubsurfaceFlowFinalizeRun
!    procedure, public :: InitializeTimestep => PMSubsurfaceFlowInitializeTimestep
    procedure, public :: FinalizeTimestep => PMSubsurfaceFlowFinalizeTimestep
    procedure, public :: PreSolve => PMSubsurfaceFlowPreSolve
    procedure, public :: PostSolve => PMSubsurfaceFlowPostSolve
    procedure, public :: AcceptSolution => PMSubsurfaceFlowAcceptSolution
!    procedure, public :: TimeCut => PMSubsurfaceFlowTimeCut
!    procedure, public :: UpdateSolution => PMSubsurfaceFlowUpdateSolution
    procedure, public :: UpdateAuxvars => PMSubsurfaceFlowUpdateAuxvars
    procedure, public :: CheckpointBinary => PMSubsurfaceFlowCheckpointBinary
    procedure, public :: RestartBinary => PMSubsurfaceFlowRestartBinary
#if defined(PETSC_HAVE_HDF5)
    procedure, public :: CheckpointHDF5 => PMSubsurfaceFlowCheckpointHDF5
    procedure, public :: RestartHDF5 => PMSubsurfaceFlowRestartHDF5
#endif
    procedure, public :: InputRecord => PMSubsurfaceFlowInputRecord
    procedure  :: AllWellsInit
    procedure :: AllWellsUpdate
!    procedure, public :: Destroy => PMSubsurfaceFlowDestroy
  end type pm_subsurface_flow_type
  
  public :: PMSubsurfaceFlowCreate, &
            PMSubsurfaceFlowSetup, &
            PMSubsurfaceFlowSetupSolvers, &
            PMSubsurfaceFlowInitializeTimestepA, &
            PMSubsurfaceFlowInitializeTimestepB, &
            PMSubsurfaceFlowInitializeRun, &
            PMSubsurfaceFlowUpdateSolution, &
            PMSubsurfaceFlowUpdatePropertiesNI, &
            PMSubsurfaceFlowTimeCut, &
            PMSubsurfaceFlowLimitDTByCFL, &
            PMSubsurfaceFlowCheckpointBinary, &
            PMSubsurfaceFlowRestartBinary, &
            PMSubsurfaceFlowReadSelectCase, &
            PMSubsurfaceFlowDestroy
  
contains

! ************************************************************************** !

subroutine PMSubsurfaceFlowCreate(this)
  ! 
  ! Intializes shared members of subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  nullify(this%realization)
  nullify(this%comm1)
  this%transient_permeability = PETSC_FALSE
  this%store_porosity_for_ts_cut = PETSC_FALSE
  this%store_porosity_for_transport = PETSC_FALSE
  this%check_post_convergence = PETSC_FALSE
  
  ! defaults
  this%max_pressure_change = 0.d0
  this%max_temperature_change = 0.d0
  this%max_saturation_change = 0.d0
  this%max_xmol_change = 0.d0
  this%pressure_change_governor = 5.d5
  this%temperature_change_governor = 5.d0
  this%saturation_change_governor = 0.5d0
  this%xmol_change_governor = 1.d0
  this%cfl_governor = UNINITIALIZED_DOUBLE
  this%pressure_dampening_factor = UNINITIALIZED_DOUBLE
  this%saturation_change_limit = UNINITIALIZED_DOUBLE
  this%pressure_change_limit = UNINITIALIZED_DOUBLE
  this%temperature_change_limit = UNINITIALIZED_DOUBLE
  
  call PMBaseInit(this)

end subroutine PMSubsurfaceFlowCreate

! ************************************************************************** !

subroutine PMSubsurfaceFlowReadSelectCase(this,input,keyword,found,option)
  ! 
  ! Reads input file parameters associated with the subsurface flow process 
  !       model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/05/16

  use Input_Aux_module
  use String_module
  use Option_module
 
  implicit none
  
  class(pm_subsurface_flow_type) :: this
  type(input_type) :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  type(option_type) :: option

  found = PETSC_TRUE
  select case(trim(keyword))
  
    case('MAX_PRESSURE_CHANGE')
      call InputReadDouble(input,option,this%pressure_change_governor)
      call InputDefaultMsg(input,option,'dpmxe')

    case('MAX_TEMPERATURE_CHANGE')
      call InputReadDouble(input,option,this%temperature_change_governor)
      call InputDefaultMsg(input,option,'dtmpmxe')
  
    case('MAX_CONCENTRATION_CHANGE')
      call InputReadDouble(input,option,this%xmol_change_governor)
      call InputDefaultMsg(input,option,'dcmxe')

    case('MAX_SATURATION_CHANGE')
      call InputReadDouble(input,option,this%saturation_change_governor)
      call InputDefaultMsg(input,option,'dsmxe')

    case('PRESSURE_DAMPENING_FACTOR')
      call InputReadDouble(input,option,this%pressure_dampening_factor)
      call InputErrorMsg(input,option,'PRESSURE_DAMPENING_FACTOR', &
                         'SUBSURFACE_FLOW OPTIONS')

    case('SATURATION_CHANGE_LIMIT')
      call InputReadDouble(input,option,this%saturation_change_limit)
      call InputErrorMsg(input,option,'SATURATION_CHANGE_LIMIT', &
                          'SUBSURFACE_FLOW OPTIONS')
                           
    case('PRESSURE_CHANGE_LIMIT')
      call InputReadDouble(input,option,this%pressure_change_limit)
      call InputErrorMsg(input,option,'PRESSURE_CHANGE_LIMIT', &
                          'SUBSURFACE_FLOW OPTIONS')
                           
    case('TEMPERATURE_CHANGE_LIMIT')
      call InputReadDouble(input,option,this%temperature_change_limit)
      call InputErrorMsg(input,option,'TEMPERATURE_CHANGE_LIMIT', &
                          'SUBSURFACE_FLOW OPTIONS')

    case('MAX_CFL')
      call InputReadDouble(input,option,this%cfl_governor)
      call InputErrorMsg(input,option,'MAX_CFL', &
                          'SUBSURFACE_FLOW OPTIONS')

    case default
      found = PETSC_FALSE
  end select  
  
end subroutine PMSubsurfaceFlowReadSelectCase

! ************************************************************************** !

subroutine PMSubsurfaceFlowSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  PetscErrorCode :: ierr

  ! set the communicator
  this%comm1 => this%realization%comm1
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mnrl_surf_with_porosity) then
      this%store_porosity_for_ts_cut = PETSC_TRUE
      this%store_porosity_for_transport = PETSC_TRUE
    endif
  endif
  if (this%option%flow%transient_porosity) then
    this%store_porosity_for_ts_cut = PETSC_TRUE
    if (this%option%ntrandof > 0) then
      this%store_porosity_for_transport = PETSC_TRUE
    endif
  endif
  
end subroutine PMSubsurfaceFlowSetup

! ************************************************************************** !

subroutine PMSubsurfaceFlowSetupSolvers(this,solver)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/03/14

  use Solver_module
  
  implicit none
  
  class(pm_subsurface_flow_type) :: this
  type(solver_type) :: solver
  
  PetscErrorCode :: ierr
  
end subroutine PMSubsurfaceFlowSetupSolvers

! ************************************************************************** !

subroutine PMSubsurfaceFlowSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Realization_Subsurface_class
  use Grid_module

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

  this%solution_vec = realization%field%flow_xx
  this%residual_vec = realization%field%flow_r
  
end subroutine PMSubsurfaceFlowSetRealization

! ************************************************************************** !

recursive subroutine PMSubsurfaceFlowInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 
  use Condition_Control_module
  use Material_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_MINERAL, POROSITY_CURRENT

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  PetscBool :: update_initial_porosity

  ! overridden in pm_general only
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity) then
      call RealizationCalcMineralPorosity(this%realization)
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   POROSITY,POROSITY_MINERAL)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity0)
    endif
  endif
  
  ! restart
  if (this%option%restart_flag .and. this%option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(this%realization)
!geh: for testing only.  In general, we only revert parameter, not flow.
!    call CondControlAssignFlowInitCond(this%realization)
!    call this%UpdateAuxVars()
  endif
  ! update material properties that are a function of mineral vol fracs
  update_initial_porosity = PETSC_TRUE
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mineral_surface_area) then
      call RealizationUpdatePropertiesTS(this%realization)
      update_initial_porosity = PETSC_FALSE
    endif
  endif
  if (update_initial_porosity) then
    call this%comm1%GlobalToLocal(this%realization%field%porosity0, &
                                  this%realization%field%work_loc)
    ! push values to porosity_base
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_MINERAL)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_CURRENT)
  endif  

  call this%PreSolve()
  call this%UpdateAuxVars()
  call this%UpdateSolution() 
  call this%AllWellsInit() !does nothing if no well exist
    
end subroutine PMSubsurfaceFlowInitializeRun

! ************************************************************************** !

subroutine AllWellsInit(this)
  !
  ! Initialise all wells - does nothing if no well exist
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/16

  !use Well_Base_class
  use Coupler_module
  implicit none

  class(pm_subsurface_flow_type) :: this

  type(coupler_type), pointer :: source_sink

  PetscMPIInt :: cur_w_myrank
  character(len=MAXSTRINGLENGTH) :: wfile_name
  PetscInt :: beg_cpl_conns, end_cpl_conns
  PetscInt :: ierr 

  source_sink => this%realization%patch%source_sink_list%first

  beg_cpl_conns = 1
  do
    if (.not.associated(source_sink)) exit
    if( associated(source_sink%well) ) then
      !exlude empty wells - not included in well comms
      if(source_sink%connection_set%num_connections > 0) then
        !TO DO - move all of this chunk of code into (well_base.F90) 
        !        well%InitRun(Material%auxvars,grid,option,source_sink%name) 
  
        end_cpl_conns = beg_cpl_conns + &
                        source_sink%connection_set%num_connections - 1  

#ifdef WELL_DEBUG
        print *,"AllWellsInit - cntrl_lcell_id", source_sink%well%cntrl_lcell_id
        print *, "perm_yy", this%realization%patch%aux%Material% &
                                  auxvars(1)%permeability(2)
#endif   
        
        call source_sink%well%WellFactorUpdate(this%realization%patch%grid, &
                                        source_sink%connection_set, &
                                this%realization%patch%aux%Material%auxvars, &
                                           this%realization%option)
        !need to initialize pressure and well densities for injectectors
        call source_sink%well%InitDensity(this%realization%patch%grid, &
                                          this%realization%option )
        call source_sink%well%ExplUpdate(this%realization%patch%grid, &
                      this%realization%patch% &
                      ss_flow_vol_fluxes(:,beg_cpl_conns:end_cpl_conns), &
                                         this%realization%option)
        call source_sink%well%TempUpdate(this%realization%patch%grid, &
                             this%realization%option)
        call source_sink%well%HydroCorrUpdates(this%realization%patch%grid, &
                      this%realization%patch% &
                      ss_flow_vol_fluxes(:,beg_cpl_conns:end_cpl_conns), &
                      this%realization%option)
        !update the pressure again after H correction, 
        ! only to print the right value at t=0
        call source_sink%well%ExplUpdate(this%realization%patch%grid, &
                      this%realization%patch% &
                      ss_flow_vol_fluxes(:,beg_cpl_conns:end_cpl_conns), &
                                         this%realization%option)

        
       ! create well outputfile - should be moved into a well class
       ! For now open files to print the well variables by default 
       ! TODO: add to well_spec user options to control well printing
       call MPI_Comm_rank(source_sink%well%comm, cur_w_myrank, ierr )  
       if(source_sink%well%cntr_rank == cur_w_myrank ) then
         !w_file_id = source_sink%id 
         wfile_name = trim(this%realization%option%global_prefix) // "_" // &
                      trim(source_sink%name) // ".tec" 
         !open(unit=w_file_id,file=wfile_name)
         !can close the file, and reopen it in appending mode later??
         open(unit=IUNIT_TEMP,file=wfile_name)
         call source_sink%well%PrintOutputHeader(this%realization% &
                                                 output_option,IUNIT_TEMP)
         close(unit=IUNIT_TEMP)
       end if

      end if
    end if
    source_sink => source_sink%next 
  end do 

end subroutine AllWellsInit

! ************************************************************************** !

subroutine PMSubsurfaceFlowInitializeTimestepA(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  use Material_module
  use Material_Aux_class, only : POROSITY_MINERAL
  
  implicit none
  
  class(pm_subsurface_flow_type) :: this

  this%option%flow_dt = this%option%dt

  if (this%transient_permeability) then
  !geh:remove
    call MaterialAuxVarCommunicate(this%comm1, &
                                   this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   PERMEABILITY_X,ZERO_INTEGER)
    call MaterialAuxVarCommunicate(this%comm1, &
                                   this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   PERMEABILITY_Y,ZERO_INTEGER)
    call MaterialAuxVarCommunicate(this%comm1, &
                                   this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   PERMEABILITY_Z,ZERO_INTEGER)
  endif

  if (this%store_porosity_for_ts_cut) then
    ! store base properties for reverting at time step cut
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY, &
                                 POROSITY_MINERAL)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity_base_store)
  endif

end subroutine PMSubsurfaceFlowInitializeTimestepA

! ************************************************************************** !

subroutine PMSubsurfaceFlowInitializeTimestepB(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
  use Material_module
  use Material_Aux_class, only : POROSITY_CURRENT
  
  implicit none
  
  class(pm_subsurface_flow_type) :: this

  if (this%option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(this%realization,TIME_T,this%option%time)
    if (this%store_porosity_for_transport) then
      ! store time t properties for transport
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc,POROSITY, &
                                   POROSITY_CURRENT)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity_t)
    endif
  endif 
  
  ! update porosity for time t+dt
  if (associated(this%realization%reaction)) then
    if (this%realization%reaction%update_porosity .or. &
        this%realization%reaction%update_tortuosity .or. &
        this%realization%reaction%update_permeability .or. &
        this%realization%reaction%update_mineral_surface_area) then
      call RealizationUpdatePropertiesTS(this%realization)
    endif
  endif

  call this%AllWellsUpdate()
  
end subroutine PMSubsurfaceFlowInitializeTimestepB

! ************************************************************************** !

subroutine AllWellsUpdate(this)
  !
  ! Update all wells at the beginning of each time step
  !  - is permeability changes updates well factor 
  !  - update hydrostatic corrections
  !  - 
  ! 
  ! Author: Paolo Orsini
  ! Date: 06/06/16

  use Coupler_module
  implicit none

  class(pm_subsurface_flow_type) :: this

  type(coupler_type), pointer :: source_sink

  PetscInt :: beg_cpl_conns, end_cpl_conns

  source_sink => this%realization%patch%source_sink_list%first

  beg_cpl_conns = 1
  do
    if (.not.associated(source_sink)) exit
    if( associated(source_sink%well) ) then
      !exlude empty wells - not included in well comms
      if(source_sink%connection_set%num_connections > 0) then
        !TO DO: move this chunk of code in well_base class
        !       well%InitTimeStep() 
        if(this%realization%option%update_flow_perm) then
          call source_sink%well%WellFactorUpdate(this%realization%patch%grid, &
                                         source_sink%connection_set, &
                                 this%realization%patch%aux%Material%auxvars, &
                                            this%realization%option)
        end if

        end_cpl_conns = beg_cpl_conns + &
                        source_sink%connection_set%num_connections - 1  

        ! for fully explicit well - or extra coupling
        !call source_sink%well%ExplUpdate(this%realization%patch%grid, &
        !              this%realization%patch% &
        !             ss_flow_vol_fluxes(:,beg_cpl_conns:end_cpl_conns), &
        !                                 this%realization%option)
        !commenting this freezes the hydrostatic corrections at the intial solution 
        call source_sink%well%TempUpdate(this%realization%patch%grid, &
                             this%realization%option)
        call source_sink%well%HydroCorrUpdates(this%realization%patch%grid, &
                      this%realization%patch% &
                      ss_flow_vol_fluxes(:,beg_cpl_conns:end_cpl_conns), &
                      this%realization%option)

      end if
    end if
    source_sink => source_sink%next 
  end do 

end subroutine AllWellsUpdate

! ************************************************************************** !

subroutine PMSubsurfaceFlowPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Global_module

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  this%option%io_buffer = 'PMSubsurfaceFlowPreSolve() must be extended.'
  call printErrMsg(this%option)  

end subroutine PMSubsurfaceFlowPreSolve

! ************************************************************************** !

subroutine PMSubsurfaceFlowPostSolve(this)
  ! 
  ! PMSubsurfaceFlowUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  this%option%io_buffer = 'PMSubsurfaceFlowPostSolve() must be extended.'
  call printErrMsg(this%option)  
  
end subroutine PMSubsurfaceFlowPostSolve

! ************************************************************************** !

function PMSubsurfaceFlowAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  PetscBool :: PMSubsurfaceFlowAcceptSolution
  
  ! do nothing
  PMSubsurfaceFlowAcceptSolution = PETSC_TRUE
  
end function PMSubsurfaceFlowAcceptSolution

! ************************************************************************** !

subroutine PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  call RealizationUpdatePropertiesNI(this%realization)

end subroutine PMSubsurfaceFlowUpdatePropertiesNI

! ************************************************************************** !

subroutine PMSubsurfaceFlowTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 
  use Material_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_MINERAL
  
  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  PetscErrorCode :: ierr
  
  this%option%flow_dt = this%option%dt
  call VecCopy(this%realization%field%flow_yy, &
               this%realization%field%flow_xx,ierr);CHKERRQ(ierr)
  if (this%store_porosity_for_transport) then
    ! store base properties for reverting at time step cut
    call this%comm1%GlobalToLocal(this%realization%field%porosity_base_store, &
                                  this%realization%field%work_loc)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,POROSITY, &
                                 POROSITY_MINERAL)
  endif             

end subroutine PMSubsurfaceFlowTimeCut

! ************************************************************************** !

subroutine PMSubsurfaceFlowLimitDTByCFL(this,dt)
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/16 
  !
  use Option_module
  use Output_Aux_module

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  PetscReal :: dt

  PetscReal :: max_dt_cfl_1
  PetscReal :: prev_dt
  type(output_option_type), pointer :: output_option
  
  if (Initialized(this%cfl_governor)) then
    call RealizationCalculateCFL1Timestep(this%realization,max_dt_cfl_1) 
    if (dt/this%cfl_governor > max_dt_cfl_1) then
      prev_dt = dt
      dt = max_dt_cfl_1*this%cfl_governor
      output_option => this%realization%output_option
      if (OptionPrintToScreen(this%option)) then
        write(*, &
          '(" CFL Limiting (",f4.1,"): ",1pe12.4," -> ",1pe12.4," [",a,"]")') &
              this%cfl_governor,prev_dt/output_option%tconv, &
              dt/output_option%tconv,trim(output_option%tunit)
      endif
      if (OptionPrintToFile(this%option)) then
        write(this%option%fid_out, &
          '(" CFL Limiting (",f4.1,"): ",1pe12.4," -> ",1pe12.4," [",a,"]")') &
              this%cfl_governor,prev_dt/output_option%tconv, &
              dt/output_option%tconv,trim(output_option%tunit)
      endif
    endif
  endif

end subroutine PMSubsurfaceFlowLimitDTByCFL

! ************************************************************************** !

subroutine PMSubsurfaceFlowFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use Material_module
  use Global_module
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_CURRENT

  implicit none
  
  class(pm_subsurface_flow_type) :: this

  if (this%option%ntrandof > 0) then 
    ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(this%realization,TIME_TpDT,this%option%time)
    if (this%store_porosity_for_transport) then
      ! store time t properties for transport
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc,POROSITY, &
                                   POROSITY_CURRENT)
      call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                    this%realization%field%porosity_tpdt)
    endif    
  endif
  
  call this%MaxChange()
  
end subroutine PMSubsurfaceFlowFinalizeTimestep

! ************************************************************************** !

subroutine PMSubsurfaceFlowUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Condition_module
  use Integral_Flux_module
  use SrcSink_Sandbox_module

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  PetscBool :: force_update_flag = PETSC_FALSE
  PetscErrorCode :: ierr

  call VecCopy(this%realization%field%flow_xx, &
               this%realization%field%flow_yy,ierr);CHKERRQ(ierr)
  
  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%realization%flow_conditions, &
                           this%realization%option, &
                           this%realization%option%time)
  call SSSandboxUpdate(ss_sandbox_list,this%realization%option%time, &
                       this%realization%option,this%realization%output_option)
  ! right now, RealizUpdateAllCouplerAuxVars only updates flow
  call RealizUpdateAllCouplerAuxVars(this%realization,force_update_flag)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif
  call IntegralFluxUpdate(this%realization%patch%integral_flux_list, &
                          this%realization%patch%internal_flow_fluxes, &
                          this%realization%patch%boundary_flow_fluxes, &
                          INTEGRATE_FLOW,this%option)
  ! end from RealizationUpdate()

end subroutine PMSubsurfaceFlowUpdateSolution  

! ************************************************************************** !

subroutine PMSubsurfaceFlowUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this

  this%option%io_buffer = 'PMSubsurfaceFlowUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMSubsurfaceFlowUpdateAuxVars   

! ************************************************************************** !

subroutine PMSubsurfaceFlowCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_subsurface_flow_type) :: this
  PetscViewer :: viewer
  
  call CheckpointFlowProcessModelBinary(viewer,this%realization) 
  
end subroutine PMSubsurfaceFlowCheckpointBinary

! ************************************************************************** !

subroutine PMSubsurfaceFlowRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_subsurface_flow_type) :: this
  PetscViewer :: viewer
  
  call RestartFlowProcessModelBinary(viewer,this%realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()
  
end subroutine PMSubsurfaceFlowRestartBinary

! ************************************************************************** !

#if defined(PETSC_HAVE_HDF5)
subroutine PMSubsurfaceFlowCheckpointHDF5(this, pm_grp_id)
  !
  ! Checkpoints data associated with Subsurface PM
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15

  use Checkpoint_module
  use hdf5

  implicit none

  class(pm_subsurface_flow_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

  call CheckpointFlowProcessModelHDF5(pm_grp_id, this%realization)

end subroutine PMSubsurfaceFlowCheckpointHDF5

! ************************************************************************** !

subroutine PMSubsurfaceFlowRestartHDF5(this, pm_grp_id)
  !
  ! Checkpoints data associated with Subsurface PM
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15

  use Checkpoint_module
  use hdf5

  implicit none

  class(pm_subsurface_flow_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

  call RestartFlowProcessModelHDF5(pm_grp_id, this%realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()


end subroutine PMSubsurfaceFlowRestartHDF5
#endif

! ************************************************************************** !

recursive subroutine PMSubsurfaceFlowFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMSubsurfaceFlowFinalizeRun

! ************************************************************************** !

subroutine PMSubsurfaceFlowInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_subsurface_flow_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMSubsurfaceFlowInputRecord

! ************************************************************************** !

subroutine PMSubsurfaceFlowDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  implicit none
  
  class(pm_subsurface_flow_type) :: this
  
  ! destroyed in realization
  nullify(this%comm1)
  nullify(this%option)
  nullify(this%output_option)
  
end subroutine PMSubsurfaceFlowDestroy
  
end module PM_Subsurface_Flow_class
