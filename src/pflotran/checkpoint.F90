! The header for the binary checkpoint files.
! This needs to go in its own module for some reason, because otherwise 
! some compilers won't like the INTERFACE block needed to allow us to 
! use the PetscBagGetData() routine.
! RTM: This is pretty makeshift.  We need to think about what should 
! go into this header and how it should be organized.

module Checkpoint_Header_module
  implicit none
  private
  type, public :: checkpoint_header_type
    real*8 :: flow_time
    real*8 :: flow_dt
    integer*8 :: flow_steps
    integer*8 :: flow_newtcum
    integer*8 :: flow_icutcum
    integer*8 :: flow_linear_cum
    integer*8 :: flow_num_const_timesteps
    integer*8 :: flow_num_newton_iterations
    real*8 :: tran_time
    real*8 :: tran_dt
    integer*8 :: tran_steps
    integer*8 :: tran_newtcum
    integer*8 :: tran_icutcum
    integer*8 :: tran_linear_cum
    integer*8 :: tran_num_const_timesteps
    integer*8 :: tran_num_newton_iterations
    integer*8 :: plot_number
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
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petscdef.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
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
                      flow_steps,flow_newtcum,flow_icutcum, &
                      flow_linear_cum, &
                      flow_num_const_timesteps, &
                      flow_num_newton_iterations, &
                      tran_steps,tran_newtcum,tran_icutcum, &
                      tran_linear_cum, &
                      tran_num_const_timesteps, &
                      tran_num_newton_iterations, &
                      id)

  use Realization_module
  use Discretization_module
  use Option_module
  use Field_module
  use Logging_module
  
  use MPHASE_module

  implicit none

  type(realization_type) :: realization
  PetscInt :: flow_num_const_timesteps
  PetscInt :: flow_num_newton_iterations
  PetscInt :: flow_steps, flow_newtcum, flow_icutcum, flow_linear_cum
  PetscInt :: tran_num_const_timesteps
  PetscInt :: tran_num_newton_iterations
  PetscInt :: tran_steps, tran_newtcum, tran_icutcum, tran_linear_cum
  PetscInt :: id
#ifdef PetscSizeT
  PetscSizeT :: bagsize
#else
  ! PETSC_SIZEOF_SIZE_T isn't defined, so we just have to assume that it 
  ! is 8 bytes.  This is dangerous, but what can we do?
  integer*8 :: bagsize
#endif

  character(len=MAXSTRINGLENGTH) :: fname
  PetscViewer :: viewer
  PetscBag :: bag
  type(checkpoint_header_type), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend  
  
  Vec :: global_vec, global_var
  PetscInt :: int_flag
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)
  call PetscLogEventBegin(logging%event_checkpoint, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)  
    
  field => realization%field
  option => realization%option
  discretization => realization%discretization
  output_option => realization%output_option

  ! Open the checkpoint file.
  call PetscGetTime(tstart,ierr)   
  if (id < 0) then
    fname = 'restart.chk'
  else if (id < 10) then
    write(fname, '(a12,i1)') 'pflotran.chk', id
  else if (id < 100) then
    write(fname, '(a12,i2)') 'pflotran.chk', id
  else if (id < 1000) then
    write(fname, '(a12,i3)') 'pflotran.chk', id
  else if (id < 10000) then
    write(fname, '(a12,i4)') 'pflotran.chk', id
  else if (id < 100000) then
    write(fname, '(a12,i5)') 'pflotran.chk', id
  else if (id < 1000000) then
    write(fname, '(a12,i6)') 'pflotran.chk', id
  endif
  call PetscViewerBinaryOpen(option%comm, fname, FILE_MODE_WRITE, &
                             viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  bagsize = 136
  call PetscBagCreate(option%comm, bagsize, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  ! Register variables that are passed into timestepper().
  call PetscBagRegisterInt(bag,header%plot_number,output_option%plot_number, &
                           "plot_number","plot_number",ierr)
  ! FLOW
  call PetscBagRegisterInt(bag,header%flow_num_newton_iterations, &
                           flow_num_newton_iterations, &
                           "flow_num_newton_iterations", &
                           "Number of flow Newton iterations in last SNES solve", &
                           ierr)
  call PetscBagRegisterInt(bag,header%flow_num_const_timesteps, &
                           flow_num_const_timesteps, &
                           "flow_num_const_timesteps", &
                           "flow_num_const_timesteps",ierr)
  ! TRANSPORT
  call PetscBagRegisterInt(bag,header%tran_num_newton_iterations, &
                           tran_num_newton_iterations, &
                           "tran_num_newton_iterations", &
                           "Number of transport Newton iterations in last SNES solve", &
                           ierr)
  call PetscBagRegisterInt(bag,header%tran_num_const_timesteps, &
                           tran_num_const_timesteps, &
                           "tran_num_const_timesteps", &
                           "tran_num_const_timesteps",ierr)
  
  ! Register relevant components of the stepper.
  ! FLOW
  call PetscBagRegisterReal(bag,header%flow_time,option%flow_time,"flow_time", &
                            "Flow Simulation time (seconds)",ierr)
  call PetscBagRegisterReal(bag,header%flow_dt,option%flow_dt,"flow_dt", &
                            "Current size of flow timestep (seconds)",ierr)
  call PetscBagRegisterInt(bag,header%flow_steps,flow_steps,"flow_steps", &
                            "Total number of flow steps taken",ierr)
  call PetscBagRegisterInt(bag,header%flow_newtcum,flow_newtcum,"flow_newtcum", &
                            "Total number of flow Newton steps taken",ierr)
  call PetscBagRegisterInt(bag,header%flow_icutcum,flow_icutcum,"flow_icutcum", &
                            "Total number of flow time step cuts",ierr)
  call PetscBagRegisterInt(bag,header%flow_linear_cum,flow_linear_cum,"flow_linear_cum", &
                            "Total number of flow linear iterations",ierr)
                            
  ! TRANSPORT
  call PetscBagRegisterReal(bag,header%tran_time,option%tran_time,"tran_time", &
                            "Transport Simulation time (seconds)",ierr)
  call PetscBagRegisterReal(bag,header%tran_dt,option%tran_dt,"tran_dt", &
                            "Current size of transport timestep (years)",ierr)
                            
  call PetscBagRegisterInt(bag,header%tran_steps,tran_steps,"tran_steps", &
                            "Total number of transport steps taken",ierr)
  call PetscBagRegisterInt(bag,header%tran_newtcum,tran_newtcum,"tran_newtcum", &
                            "Total number of transport Newton steps taken",ierr)
  call PetscBagRegisterInt(bag,header%tran_icutcum,tran_icutcum,"tran_icutcum", &
                            "Total number of transport time step cuts",ierr)
  call PetscBagRegisterInt(bag,header%tran_linear_cum,tran_linear_cum,"tran_linear_cum", &
                            "Total number of transport linear iterations",ierr)                            

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  if (option%nflowdof > 0) then
    ! grid%flow_xx is the vector into which all of the primary variables are 
    ! packed for the SNESSolve().
    call VecView(field%flow_xx, viewer, ierr)

    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! If we are running with multiple phases, we need to dump the vector 
    ! that indicates what phases are present, as well as the 'var' vector 
    ! that holds variables derived from the primary ones via the translator.
    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,RICHARDS_MODE)
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

    call VecDestroy(global_vec,ierr)
  
  endif

  if (option%ntrandof > 0) then
    call VecView(field%tran_xx, viewer, ierr)
  endif

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

  write(option%io_buffer,'(" --> Dump checkpoint file: ", a16)') trim(fname)
  call printMsg(option)

  call PetscGetTime(tend,ierr) 
  write(option%io_buffer, &
        '("      Seconds to write to checkpoint file: ", f6.4)') tend-tstart
  call printMsg(option)

  call PetscLogEventEnd(logging%event_checkpoint, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)  
  call PetscLogStagePop(ierr)

end subroutine Checkpoint

subroutine Restart(realization, &
                   flow_steps,flow_newtcum,flow_icutcum, &
                   flow_linear_cum, &
                   flow_num_const_timesteps, &
                   flow_num_newton_iterations, &
                   tran_steps,tran_newtcum,tran_icutcum, &
                   tran_linear_cum, &
                   tran_num_const_timesteps, &
                   tran_num_newton_iterations)

  use Realization_module
  use Discretization_module
  use Option_module
  use Field_module
  use Logging_module

  use MPHASE_module

  implicit none

  type(realization_type) :: realization
  PetscInt :: flow_num_const_timesteps
  PetscInt :: flow_num_newton_iterations
  PetscInt :: flow_steps, flow_newtcum, flow_icutcum, flow_linear_cum
  PetscInt :: tran_num_const_timesteps
  PetscInt :: tran_num_newton_iterations
  PetscInt :: tran_steps, tran_newtcum, tran_icutcum, tran_linear_cum

  PetscViewer viewer
  PetscBag bag
  type(checkpoint_header_type), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend

  Vec :: global_vec, global_var
  PetscInt :: int_flag
  
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  call PetscLogEventBegin(logging%event_restart, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)  
  
  field => realization%field
  option => realization%option
  discretization => realization%discretization
  output_option => realization%output_option
  
  call PetscGetTime(tstart,ierr)
  option%io_buffer = '--> Open checkpoint file: ' // trim(option%restart_file)
  call printMsg(option)
  call PetscViewerBinaryOpen(option%comm,option%restart_file, &
                             FILE_MODE_READ,viewer,ierr)
 
  ! Get the header data.
  call PetscBagLoad(viewer, bag, ierr)
  call PetscBagGetData(bag, header, ierr)
  output_option%plot_number = header%plot_number
  ! FLOW
  option%flow_time = header%flow_time
  option%flow_dt = header%flow_dt
  flow_num_newton_iterations = header%flow_num_newton_iterations
  flow_num_const_timesteps = header%flow_num_const_timesteps
  flow_steps = header%flow_steps
  flow_newtcum = header%flow_newtcum
  flow_icutcum = header%flow_icutcum
  flow_linear_cum = header%flow_linear_cum
  ! TRANSPORT
  option%tran_time = header%tran_time
  option%tran_dt = header%tran_dt
  tran_num_newton_iterations = header%tran_num_newton_iterations
  tran_num_const_timesteps = header%tran_num_const_timesteps
  tran_steps = header%tran_steps
  tran_newtcum = header%tran_newtcum
  tran_icutcum = header%tran_icutcum
  tran_linear_cum = header%tran_linear_cum
  call PetscBagDestroy(bag, ierr)
  
  ! Load the PETSc vectors.
  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)

    call VecLoadIntoVector(viewer,field%flow_xx,ierr)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr)
    
    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,RICHARDS_MODE)
        call VecLoadIntoVector(viewer, global_vec, ierr)      
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         field%iphas_loc,ONEDOF)
        call VecCopy(field%iphas_loc,field%iphas_old_loc,ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                        field%iphas_old_loc,ONEDOF)
        if (option%iflowmode == MPH_MODE) then
        ! set vardof vec in mphase
        endif
      case default
    end select
    
    call VecLoadIntoVector(viewer, global_vec, ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%porosity_loc,ONEDOF)
    call VecLoadIntoVector(viewer,global_vec,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_xx_loc,ONEDOF)
    call VecLoadIntoVector(viewer,global_vec,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_yy_loc,ONEDOF)
    call VecLoadIntoVector(viewer,global_vec,ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     field%perm_zz_loc,ONEDOF)
    
    call VecDestroy(global_vec,ierr)
  
  endif
  
  if (option%ntrandof > 0) then
    call VecLoadIntoVector(viewer,field%tran_xx,ierr)
    call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                     field%tran_xx_loc,NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr)
  endif
  
  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)
  call PetscGetTime(tend,ierr) 

  write(option%io_buffer, &
        '("      Seconds to read to checkpoint file: ", f6.4)') tend-tstart
  call printMsg(option)

  call PetscLogEventEnd(logging%event_restart, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)  
  
end subroutine Restart

end module Checkpoint_module
