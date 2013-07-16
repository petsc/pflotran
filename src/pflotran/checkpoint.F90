! The header for the binary checkpoint files.
! This needs to go in its own module for some reason, because otherwise 
! some compilers won't like the INTERFACE block needed to allow us to 
! use the PetscBagGetData() routine.
! RTM: This is pretty makeshift.  We need to think about what should 
! go into this header and how it should be organized.

module Checkpoint_Header_module
  implicit none
  private
  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  ! Currently:
  !  PetscReal: 8
  !  PetscInt:  19
  !  Total: 27 * 8 = 216
! IMPORTANT: If you change the contents of the header, you MUST update 
! 'bagsize' or risk corrupting memory.
#ifdef PetscSizeT
  PetscSizeT, parameter :: bagsize = 216
#else
  ! PETSC_SIZEOF_SIZE_T isn't defined, so we just have to assume that it 
  ! is 8 bytes.  This is dangerous, but what can we do?
  integer*8, parameter :: bagsize = 216
#endif
  public :: bagsize
  type, public :: checkpoint_header_type
    integer*8 :: revision_number  ! increment this every time there is a change
    integer*8 :: plot_number      ! in the checkpoint file format
    integer*8 :: match_waypoint_flag

    integer*8 :: grid_discretization_type

    integer*8 :: nflowdof
    real*8 :: flow_time
    real*8 :: flow_dt
    real*8 :: flow_prev_dt
    integer*8 :: flow_time_steps
    integer*8 :: flow_cumulative_newton_iterations
    integer*8 :: flow_cumulative_time_step_cuts
    integer*8 :: flow_cumulative_linear_iterations
    integer*8 :: flow_num_constant_time_steps
    integer*8 :: flow_num_newton_iterations
    real*8 :: flow_cumulative_solver_time  ! don't implement yet; will screw up restarts

    integer*8 :: ntrandof
    real*8 :: tran_time
    real*8 :: tran_dt
    real*8 :: tran_prev_dt
    integer*8 :: tran_time_steps
    integer*8 :: tran_cumulative_newton_iterations
    integer*8 :: tran_cumulative_time_step_cuts
    integer*8 :: tran_cumulative_linear_iterations
    integer*8 :: tran_num_constant_time_steps
    integer*8 :: tran_num_newton_iterations
    real*8 :: tran_cumulative_solver_time  ! don't implement yet; will screw up restarts
    integer*8 :: checkpoint_activity_coefs
  end type checkpoint_header_type
end module Checkpoint_Header_module

module Checkpoint_module

  use Checkpoint_Header_module

  implicit none
  
  private

  public :: Checkpoint, Restart

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdef.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

Interface PetscBagGetData
Subroutine PetscBagGetData(bag,ctx,ierr)
      use Checkpoint_Header_module
      PetscBag bag
      type(checkpoint_header_type), pointer :: ctx
      PetscErrorCode ierr
End Subroutine
End Interface PetscBagGetData

contains

subroutine Checkpoint(realization, &
                      flow_time_steps,flow_cumulative_newton_iterations, &
                      flow_cumulative_time_step_cuts, &
                      flow_cumulative_linear_iterations, &
                      flow_num_constant_time_steps, &
                      flow_num_newton_iterations, &
                      flow_cumulative_solver_time, &
                      flow_prev_dt, &
                      tran_time_steps,tran_cumulative_newton_iterations, &
                      tran_cumulative_time_step_cuts, &
                      tran_cumulative_linear_iterations, &
                      tran_num_constant_time_steps, &
                      tran_num_newton_iterations, &
                      tran_cumulative_solver_time, &
                      tran_prev_dt, &
                      id_string)

  use Realization_class
  use Realization_Base_class, only : RealizationGetDataset
  use Reaction_Aux_module
  use Discretization_module
  use Option_module
  use Output_Aux_module
  use Field_module
  use Logging_module
  use Grid_module
  
  use Flash2_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION

  use Reactive_Transport_module, only : RTCheckpointKineticSorption
  use String_module, only : StringNull

  implicit none

  type(realization_type) :: realization
  PetscInt :: grid_discretization_type
  PetscInt :: flow_num_constant_time_steps
  PetscInt :: flow_num_newton_iterations
  PetscInt :: flow_time_steps
  PetscInt :: flow_cumulative_newton_iterations
  PetscInt :: flow_cumulative_time_step_cuts
  PetscInt :: flow_cumulative_linear_iterations
  PetscReal :: flow_cumulative_solver_time
  PetscReal :: flow_prev_dt
  PetscInt :: tran_num_constant_time_steps
  PetscInt :: tran_num_newton_iterations
  PetscInt :: tran_time_steps
  PetscInt :: tran_cumulative_newton_iterations
  PetscInt :: tran_cumulative_time_step_cuts
  PetscInt :: tran_cumulative_linear_iterations
  PetscReal :: tran_cumulative_solver_time
  PetscReal :: tran_prev_dt
  character(len=MAXWORDLENGTH), intent(in) :: id_string ! should not be altered
  
!  PetscInt :: checkpoint_activity_coefs
!  PetscInt :: match_waypoint_flag

  character(len=MAXSTRINGLENGTH) :: filename
  PetscViewer :: viewer
  PetscBag :: bag
  type(checkpoint_header_type), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend  
  
  Vec :: global_vec
  PetscInt :: int_flag
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option
  PetscInt :: i, j, k

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)
  call PetscLogEventBegin(logging%event_checkpoint,ierr)  
    
  field => realization%field
  option => realization%option
  discretization => realization%discretization
  output_option => realization%output_option
  grid => discretization%grid 

  global_vec = 0
  ! Open the checkpoint file.
  call PetscTime(tstart,ierr)   
  if (StringNull(id_string)) then
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-restart.chk'
  else 
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-' // trim(adjustl(id_string)) // '.chk'
  endif
  !geh: To skip .info file, need to split PetscViewerBinaryOpen() 
  !     into the routines it calls so that PetscViewerBinarySkipInfo()
  !     can be called after PetscViewerSetType(), but before
  !     PetscViewerFileSetName().  See note in PETSc docs.
  !call PetscViewerBinaryOpen(option%mycomm, filename, FILE_MODE_WRITE, &
  !                           viewer, ierr)
  call PetscViewerCreate(option%mycomm,viewer,ierr)
  call PetscViewerSetType(viewer,PETSCVIEWERBINARY,ierr)
  call PetscViewerFileSetMode(viewer,FILE_MODE_WRITE,ierr)
  call PetscViewerBinarySkipInfo(viewer,ierr)
  call PetscViewerFileSetName(viewer,filename,ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  call PetscBagCreate(option%mycomm,bagsize, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  call CheckpointRegisterBagHeader(bag,header)
  ! Revision # register in PetscBagRegister since it is default.  All other 
  ! header entities default to 0 or 0.d0
  header%plot_number = output_option%plot_number
  header%match_waypoint_flag = ZERO_INTEGER
  if (option%match_waypoint) then
    header%match_waypoint_flag = ONE_INTEGER
  endif
  header%grid_discretization_type = grid%itype

  ! FLOW
  header%nflowdof = option%nflowdof
  header%flow_num_newton_iterations = flow_num_newton_iterations
  header%flow_num_constant_time_steps = flow_num_constant_time_steps

  ! TRANSPORT
  header%tran_num_newton_iterations = tran_num_newton_iterations
  header%tran_num_constant_time_steps = tran_num_constant_time_steps
  
  ! Register relevant components of the stepper.
  ! FLOW
  header%flow_time = option%flow_time
  header%flow_dt = option%flow_dt
  header%flow_prev_dt = flow_prev_dt
                            
  header%flow_time_steps = flow_time_steps
  header%flow_cumulative_newton_iterations = flow_cumulative_newton_iterations
  header%flow_cumulative_time_step_cuts = flow_cumulative_time_step_cuts
  header%flow_cumulative_linear_iterations = flow_cumulative_linear_iterations
!!$  header%flow_cumulative_solver_time = flow_cumulative_solver_time
  ! NOTE(bja, 2013-06) : zero out wall clock time so restart files
  ! will be bit for bit identical
  header%flow_cumulative_solver_time = 0.d0
                                                        
  ! TRANSPORT
  header%ntrandof = option%ntrandof
  header%tran_time = option%tran_time
  header%tran_dt = option%tran_dt
  header%tran_prev_dt = tran_prev_dt
                            
  header%tran_time_steps = tran_time_steps
  header%tran_cumulative_newton_iterations = tran_cumulative_newton_iterations
  header%tran_cumulative_time_step_cuts = tran_cumulative_time_step_cuts
  header%tran_cumulative_linear_iterations = tran_cumulative_linear_iterations
!!$ header%tran_cumulative_solver_time = tran_cumulative_solver_time
  ! NOTE(bja, 2013-06) : zero out wall clock time so restart files
  ! will be bit for bit identical
  header%tran_cumulative_solver_time = 0.d0

  if (associated(realization%reaction)) then
    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      header%checkpoint_activity_coefs = ONE_INTEGER
    else
      header%checkpoint_activity_coefs = ZERO_INTEGER
    endif
  else
    header%checkpoint_activity_coefs = ZERO_INTEGER
  endif

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! grid%flow_xx is the vector into which all of the primary variables are 
    ! packed for the SNESSolve().
    call VecView(field%flow_xx, viewer, ierr)


    ! If we are running with multiple phases, we need to dump the vector 
    ! that indicates what phases are present, as well as the 'var' vector 
    ! that holds variables derived from the primary ones via the translator.
    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,MIS_MODE, &
           FLASH2_MODE,G_MODE)
        call DiscretizationLocalToGlobal(realization%discretization, &
                                         field%iphas_loc,global_vec,ONEDOF)
        call VecView(global_vec, viewer, ierr)
       case default
    end select 

    ! Porosity and permeability.
    ! (We only write diagonal terms of the permeability tensor for now, 
    ! since we have yet to add the full-tensor formulation.)
    call DiscretizationLocalToGlobal(discretization,field%porosity_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr)
    call DiscretizationLocalToGlobal(discretization,field%perm_xx_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr)
    call DiscretizationLocalToGlobal(discretization,field%perm_yy_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr)
    call DiscretizationLocalToGlobal(discretization,field%perm_zz_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr)


    if (grid%itype == STRUCTURED_GRID_MIMETIC) then 
      call DiscretizationLocalToGlobal(discretization,field%perm_xz_loc, &
                                        global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr)

      call DiscretizationLocalToGlobal(discretization,field%perm_xy_loc, &
                                        global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr)

      call DiscretizationLocalToGlobal(discretization,field%perm_yz_loc, &
                                        global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr)

      call VecGetSize(field%flow_xx_faces, j, ierr)
      call VecGetBlockSize(field%flow_xx_faces, k, ierr)
!              write(*,*) "Size", j, "block", k
      call VecView(field%flow_xx_faces, viewer, ierr) 
    end if

  endif

  if (option%ntrandof > 0) then
    call VecView(field%tran_xx, viewer, ierr)
    ! create a global vec for writing below 
    if (global_vec == 0) then
      call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
    endif
    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      ! allocated vector
      do i = 1, realization%reaction%naqcomp
        call RealizationGetDataset(realization,global_vec, &
                                   PRIMARY_ACTIVITY_COEF,i)
        call VecView(global_vec,viewer,ierr)
      enddo
      do i = 1, realization%reaction%neqcplx
        call RealizationGetDataset(realization,global_vec, &
                                   SECONDARY_ACTIVITY_COEF,i)
        call VecView(global_vec,viewer,ierr)
      enddo
    endif
    ! mineral volume fractions for kinetic minerals
    if (realization%reaction%mineral%nkinmnrl > 0) then
      do i = 1, realization%reaction%mineral%nkinmnrl
        call RealizationGetDataset(realization,global_vec, &
                                   MINERAL_VOLUME_FRACTION,i)
        call VecView(global_vec,viewer,ierr)
      enddo
    endif
    ! sorbed concentrations for multirate kinetic sorption
    if (realization%reaction%surface_complexation%nkinmrsrfcplxrxn > 0 .and. &
        .not.option%no_checkpoint_kinetic_sorption) then
      ! PETSC_TRUE flag indicates write to file
      call RTCheckpointKineticSorption(realization,viewer,PETSC_TRUE)
    endif
  endif

  if (global_vec /= 0) then
    call VecDestroy(global_vec,ierr)
  endif
 
  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

  write(option%io_buffer,'(" --> Dump checkpoint file: ", a16)') trim(filename)
  call printMsg(option)

  call PetscTime(tend,ierr) 
  write(option%io_buffer, &
        '("      Seconds to write to checkpoint file: ", f10.2)') tend-tstart
  call printMsg(option)

  call PetscLogEventEnd(logging%event_checkpoint,ierr)  
  call PetscLogStagePop(ierr)

end subroutine Checkpoint

subroutine Restart(realization, &
                   flow_time_steps, &
                   flow_cumulative_newton_iterations, &
                   flow_cumulative_time_step_cuts, &
                   flow_cumulative_linear_iterations, &
                   flow_num_constant_time_steps, &
                   flow_num_newton_iterations, &
                   flow_cumulative_solver_time, &
                   flow_prev_dt, &
                   tran_time_steps, &
                   tran_cumulative_newton_iterations, &
                   tran_cumulative_time_step_cuts, &
                   tran_cumulative_linear_iterations, &
                   tran_num_constant_time_steps, &
                   tran_num_newton_iterations, &
                   tran_cumulative_solver_time, &
                   tran_prev_dt, &
                   flow_read, &
                   transport_read, &
                   activity_coefs_read)

  use Realization_class
  use Realization_Base_class, only : RealizationSetDataset
  use Discretization_module
  use Option_module
  use Output_Aux_module
  use Field_module
  use Logging_module
  use Grid_module

  use Flash2_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION
  
  use Reactive_Transport_module, only: RTCheckpointKineticSorption

  implicit none

  type(realization_type) :: realization
  PetscInt :: flow_num_constant_time_steps
  PetscInt :: flow_num_newton_iterations
  PetscInt :: flow_time_steps, flow_cumulative_newton_iterations
  PetscInt :: flow_cumulative_time_step_cuts, flow_cumulative_linear_iterations
  PetscReal :: flow_cumulative_solver_time
  PetscReal :: flow_prev_dt
  PetscInt :: tran_num_constant_time_steps
  PetscInt :: tran_num_newton_iterations
  PetscInt :: tran_time_steps, tran_cumulative_newton_iterations
  PetscInt :: tran_cumulative_time_step_cuts, tran_cumulative_linear_iterations
  PetscReal :: tran_cumulative_solver_time
  PetscReal :: tran_prev_dt
  PetscBool :: activity_coefs_read, flow_read, transport_read

  PetscViewer viewer
  PetscBag bag
  type(checkpoint_header_type), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend

  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: int_flag
  PetscInt :: i,j,k
  PetscInt :: read_activity_coefs
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: string
  
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  call PetscLogEventBegin(logging%event_restart,ierr)  
  
  field => realization%field
  option => realization%option
  discretization => realization%discretization
  output_option => realization%output_option
  grid => discretization%grid
  
  global_vec = 0
  local_vec = 0

  call PetscTime(tstart,ierr)
  option%io_buffer = '--> Open checkpoint file: ' // &
                     trim(option%restart_filename)
  call printMsg(option)
  call PetscViewerBinaryOpen(option%mycomm,option%restart_filename, &
                             FILE_MODE_READ,viewer,ierr)
  ! skip reading info file when loading, but not working
  call PetscViewerBinarySetSkipOptions(viewer,PETSC_TRUE,ierr)
  activity_coefs_read = PETSC_FALSE
  
  ! Get the header data.
  call PetscBagCreate(option%mycomm, bagsize, bag, ierr)
  call PetscBagGetData(bag, header, ierr)
  call CheckpointRegisterBagHeader(bag,header)
  call PetscBagLoad(viewer, bag, ierr)
  
  if (header%revision_number /= CHECKPOINT_REVISION_NUMBER) then
    write(string,*) header%revision_number
    option%io_buffer = 'The revision number # of checkpoint file (' // &
                       trim(option%restart_filename) // ', rev=' // &
                       trim(adjustl(string)) // &
                       ') does not match the current revision number' // &
                       ' of PFLOTRAN checkpoint files ('
    write(string,*) CHECKPOINT_REVISION_NUMBER
    option%io_buffer = trim(option%io_buffer) // trim(adjustl(string)) // ').'
    call printErrMsg(option)
  endif
  
  output_option%plot_number = header%plot_number
  option%match_waypoint = (header%match_waypoint_flag == ONE_INTEGER)

   if (header%grid_discretization_type /= grid%itype) then
     write(string,*) header%grid_discretization_type
     option%io_buffer = 'The discretization of checkpoint file (' // &
                       trim(option%restart_filename) // ', grid_type=' // &
                       trim(adjustl(string)) // &
                       ') does not match the discretization of the current problem' // &
                       ' grid_type= ('
    write(string,*) grid%itype
    option%io_buffer = trim(option%io_buffer) // trim(adjustl(string)) // ').'
    call printErrMsg(option)
  endif

  
  ! FLOW
  if (option%nflowdof > 0 .and. option%nflowdof == header%nflowdof) then
    option%flow_time = header%flow_time
    option%flow_dt = header%flow_dt
    flow_num_newton_iterations = header%flow_num_newton_iterations
    flow_num_constant_time_steps = header%flow_num_constant_time_steps
    flow_time_steps = header%flow_time_steps
    flow_cumulative_newton_iterations = header%flow_cumulative_newton_iterations
    flow_cumulative_time_step_cuts = header%flow_cumulative_time_step_cuts
    flow_cumulative_linear_iterations = header%flow_cumulative_linear_iterations
    flow_cumulative_solver_time = header%flow_cumulative_solver_time
    flow_cumulative_solver_time = 0.d0
    flow_prev_dt = header%flow_prev_dt
    flow_read = PETSC_TRUE
  endif
  ! TRANSPORT
  if (option%ntrandof > 0 .and. option%ntrandof == header%ntrandof) then
    option%tran_time = header%tran_time
    option%tran_dt = header%tran_dt
    tran_num_newton_iterations = header%tran_num_newton_iterations
    tran_num_constant_time_steps = header%tran_num_constant_time_steps
    tran_time_steps = header%tran_time_steps
    tran_cumulative_newton_iterations = header%tran_cumulative_newton_iterations
    tran_cumulative_time_step_cuts = header%tran_cumulative_time_step_cuts
    tran_cumulative_linear_iterations = header%tran_cumulative_linear_iterations
    tran_cumulative_solver_time = header%tran_cumulative_solver_time
    tran_cumulative_solver_time = 0.d0
    tran_prev_dt = header%tran_prev_dt
    read_activity_coefs = header%checkpoint_activity_coefs
    transport_read = PETSC_TRUE
  else if (option%ntrandof /= header%ntrandof) then
    write(string,*) header%ntrandof
    option%io_buffer = 'Number of transport dofs in restart file (' // &
                       trim(adjustl(string)) // &
           ') does not match the number of transport dofs in the input file ('
    write(string,*) option%ntrandof
    option%io_buffer = trim(option%io_buffer) // string // ')'
    call printWrnMsg(option)
  endif

  if (flow_read) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
  ! Load the PETSc vectors.
    call VecLoad(field%flow_xx,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr)  

    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,THC_MODE,THMC_MODE,RICHARDS_MODE,IMS_MODE,MIS_MODE, &
           FLASH2_MODE,G_MODE)
        call VecLoad(global_vec,viewer,ierr)      
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         field%iphas_loc,ONEDOF)
        call VecCopy(field%iphas_loc,field%iphas_old_loc,ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                        field%iphas_old_loc,ONEDOF)
        if (option%iflowmode == MPH_MODE) then
        ! set vardof vec in mphase
        endif
        if (option%iflowmode == IMS_MODE) then
        ! set vardof vec in mphase
        endif
        if (option%iflowmode == FLASH2_MODE) then
        ! set vardof vec in mphase
        endif
 
      case default
    end select
    
    call VecLoad(global_vec,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%porosity_loc,ONEDOF)
    call VecLoad(global_vec,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_xx_loc,ONEDOF)
    call VecLoad(global_vec,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_yy_loc,ONEDOF)
    call VecLoad(global_vec,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_zz_loc,ONEDOF)
    
    if (grid%itype == STRUCTURED_GRID_MIMETIC) then
      call VecLoad(global_vec,viewer,ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                        field%perm_xz_loc,ONEDOF)
      call VecLoad(global_vec,viewer,ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                        field%perm_xy_loc,ONEDOF)
      call VecLoad(global_vec,viewer,ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_yz_loc,ONEDOF)

      call VecGetSize(field%flow_xx_faces, j, ierr)
      call VecGetBlockSize(field%flow_xx_faces, k, ierr)
!              write(*,*) "Size", j, "block", k
      call VecLoad(field%flow_xx_faces, viewer,ierr)
      call DiscretizationGlobalToLocalLP(discretization, field%flow_xx_faces, field%flow_xx_loc_faces, NFLOWDOF)
      call VecCopy(field%flow_xx_faces,field%flow_yy_faces,ierr) 
    end if
    
  endif
  
  if (transport_read) then
    call VecLoad(field%tran_xx,viewer,ierr)
    call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                     field%tran_xx_loc,NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr)

    if (global_vec == 0) then
      call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
    endif    
    if (read_activity_coefs == ONE_INTEGER) then
      activity_coefs_read = PETSC_TRUE
      call DiscretizationCreateVector(discretization,ONEDOF,local_vec, &
                                      LOCAL,option)
      do i = 1, realization%reaction%naqcomp
        call VecLoad(global_vec,viewer,ierr)
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         local_vec,ONEDOF)
        call RealizationSetDataset(realization,local_vec,LOCAL, &
                                   PRIMARY_ACTIVITY_COEF,i)
      enddo
      do i = 1, realization%reaction%neqcplx
        call VecLoad(global_vec,viewer,ierr)
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         local_vec,ONEDOF)
        call RealizationSetDataset(realization,local_vec,LOCAL, &
                                   SECONDARY_ACTIVITY_COEF,i)
      enddo
    endif
    ! mineral volume fractions for kinetic minerals
    if (realization%reaction%mineral%nkinmnrl > 0) then
      do i = 1, realization%reaction%mineral%nkinmnrl
        ! have to load the vecs no matter what
        call VecLoad(global_vec,viewer,ierr)
        if (.not.option%no_restart_mineral_vol_frac) then
          call RealizationSetDataset(realization,global_vec,GLOBAL, &
                                     MINERAL_VOLUME_FRACTION,i)
        endif
      enddo
    endif
    ! sorbed concentrations for multirate kinetic sorption
    if (realization%reaction%surface_complexation%nkinmrsrfcplxrxn > 0 .and. &
        .not.option%no_checkpoint_kinetic_sorption .and. &
        ! we need to fix this.  We need something to skip over the reading
        ! of sorbed concentrations altogether if they do not exist in the
        ! checkpoint file
        .not.option%no_restart_kinetic_sorption) then
      ! PETSC_FALSE flag indicates read from file
      call RTCheckpointKineticSorption(realization,viewer,PETSC_FALSE)
    endif
  endif
    
  ! We are finished, so clean up.
  if (global_vec /= 0) then
    call VecDestroy(global_vec,ierr)
  endif
  if (local_vec /= 0) then
    call VecDestroy(local_vec,ierr)
  endif
  call PetscViewerDestroy(viewer, ierr)
  call PetscTime(tend,ierr) 

  call PetscBagDestroy(bag, ierr)

  write(option%io_buffer, &
        '("      Seconds to read to checkpoint file: ", f6.2)') tend-tstart
  call printMsg(option)

  call PetscLogEventEnd(logging%event_restart,ierr) 

  
end subroutine Restart

! ************************************************************************** !
!
! CheckpointRegisterBagHeader: Registers entities within the PETSc bag to
!                              header
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine CheckpointRegisterBagHeader(bag,header)

  implicit none
  
  PetscBag :: bag
  type(checkpoint_header_type), pointer :: header

  PetscInt :: i
  PetscErrorCode :: ierr
  
  i = CHECKPOINT_REVISION_NUMBER
  call PetscBagRegisterInt(bag,header%revision_number,i, &
                           "revision_number", &
                           "revision_number", &
                           ierr)
  ! Register variables that are passed into timestepper().
  call PetscBagRegisterInt(bag,header%plot_number,0, &
                           "plot_number", &
                           "plot_number", &
                           ierr)
  call PetscBagRegisterInt(bag,header%match_waypoint_flag,0, &
                           "match_waypoint_flag","match_waypoint_flag",ierr)
  call PetscBagRegisterInt(bag,header%grid_discretization_type, 0, &
                           "grid_discretization_type", &
                           "grid_discretization_type", &
                           ierr) 
  ! FLOW
  call PetscBagRegisterInt(bag,header%nflowdof,0, &
                           "nflowdof","Number of flow degrees of freedom",ierr)
  call PetscBagRegisterInt(bag,header%flow_num_newton_iterations,0, &
                           "flow_num_newton_iterations", &
                           "Number of flow Newton iterations in last SNES solve", &
                           ierr)
  call PetscBagRegisterInt(bag,header%flow_num_constant_time_steps,0, &
                           "flow_num_constant_time_steps", &
                           "flow_num_constant_time_steps", &
                           ierr)

  ! TRANSPORT
  call PetscBagRegisterInt(bag,header%tran_num_newton_iterations,0, &
                           "tran_num_newton_iterations", &
                           "Number of transport Newton iterations in last SNES solve", &
                           ierr)
  call PetscBagRegisterInt(bag,header%tran_num_constant_time_steps,0, &
                           "tran_num_constant_time_steps", &
                           "tran_num_constant_time_steps", &
                           ierr)
  
  ! Register relevant components of the stepper.
  ! FLOW
  call PetscBagRegisterReal(bag,header%flow_time,0.d0, &
                            "flow_time", &
                            "Flow Simulation time (seconds)", &
                            ierr)
  call PetscBagRegisterReal(bag,header%flow_dt,0.d0, &
                            "flow_dt", &
                            "Current size of flow timestep (seconds)", &
                            ierr)
  call PetscBagRegisterReal(bag,header%flow_prev_dt,0.d0, &
                            "flow_prev_dt", &
                            "Previous size of flow timestep (seconds)", &
                            ierr)
  call PetscBagRegisterInt(bag,header%flow_time_steps,0, &
                           "flow_steps", &
                           "Total number of flow steps taken", &
                           ierr)
  call PetscBagRegisterInt(bag,header%flow_cumulative_newton_iterations,0, &
                           "flow_cumulative_newton_iterations", &
                            "Total number of flow Newton steps taken", &
                           ierr)
  call PetscBagRegisterInt(bag,header%flow_cumulative_time_step_cuts,0, &
                           "flow_cumulative_time_step_cuts", &
                            "Total number of flow time step cuts", &
                           ierr)
  call PetscBagRegisterInt(bag,header%flow_cumulative_linear_iterations,0, &
                           "flow_cumulative_linear_iterations", &
                            "Total number of flow linear iterations", &
                           ierr)
  call PetscBagRegisterReal(bag,header%flow_cumulative_solver_time,0.d0, &
                            "flow_cumulative_solver_time", &
                            "flow_cumulative_solver_time", &
                            ierr)
                                                        
  ! TRANSPORT
  call PetscBagRegisterInt(bag,header%ntrandof,0, &
                           "ntrandof", &
                           "Number of transport degrees of freedom", &
                           ierr)
  call PetscBagRegisterReal(bag,header%tran_time,0.d0, &
                            "tran_time", &
                            "Transport Simulation time (seconds)", &
                            ierr)
  call PetscBagRegisterReal(bag,header%tran_dt,0.d0, &
                            "tran_dt", &
                            "Current size of transport timestep (seconds)", &
                            ierr)
  call PetscBagRegisterReal(bag,header%tran_prev_dt,0.d0, &
                            "tran_prev_dt", &
                            "Previous size of transport timestep (seconds)", &
                            ierr)
                            
  call PetscBagRegisterInt(bag,header%tran_time_steps,0, &
                           "tran_steps", &
                           "Total number of transport steps taken", &
                           ierr)
  call PetscBagRegisterInt(bag,header%tran_cumulative_newton_iterations,0, &
                           "tran_cumulative_newton_iterations", &
                           "Total number of transport Newton steps taken", &
                           ierr)
  call PetscBagRegisterInt(bag,header%tran_cumulative_time_step_cuts,0, &
                           "tran_cumulative_time_step_cuts", &
                           "Total number of transport time step cuts", &
                           ierr)
  call PetscBagRegisterInt(bag,header%tran_cumulative_linear_iterations,0, &
                           "tran_cumulative_linear_iterations", &
                           "Total number of transport linear iterations", &
                           ierr)
  call PetscBagRegisterReal(bag,header%tran_cumulative_solver_time,0.d0, &
                            "tran_cumulative_solver_time", &
                            "tran_cumulative_solver_time", &
                            ierr)

  call PetscBagRegisterInt(bag,header%checkpoint_activity_coefs,0, &
                           "checkpoint_activity_coefs", &
                           "Flag indicating whether activity coefficients were checkpointed", &
                           ierr)                            

end subroutine CheckpointRegisterBagHeader

end module Checkpoint_module
