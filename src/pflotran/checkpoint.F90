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
    real*8 :: tran_time
    real*8 :: tran_dt
    integer*8 :: steps
    integer*8 :: newtcum
    integer*8 :: icutcum
    integer*8 :: num_const_timesteps
    integer*8 :: num_newton_iterations
    integer*8 :: plot_number
  end type checkpoint_header_type
end module Checkpoint_Header_module

module Checkpoint_module

  use Checkpoint_Header_module

  implicit none
  
  private

  public :: Checkpoint, Restart

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"
#if (PETSC_VERSION_RELEASE == 0 || PETSC_VERSION_SUBMINOR == 3)
#include "include/finclude/petscbag.h"

  Interface PetscBagGetData
    Subroutine PetscBagGetData(bag,ctx,ierr)
      use Checkpoint_Header_module
      PetscBag bag
      type(checkpoint_header_type), pointer :: ctx
      PetscErrorCode ierr
    End Subroutine
  End Interface PetscBagGetData

#endif

contains

#if (PETSC_VERSION_RELEASE == 1 && PETSC_VERSION_SUBMINOR < 3)

subroutine Checkpoint(realization,steps,newtcum,icutcum, &
                               timestep_cut_flag, &
                               num_newton_iterations,id)
  use Realization_module
  use Option_module

  implicit none
  
  type(realization_type) :: realization
  logical :: timestep_cut_flag
  PetscInt :: num_newton_iterations
  PetscInt :: id, steps, newtcum, icutcum

  if(realization%option%myrank == 0) then
    print *, "Warning: Checkpoint() not supported with PETSc 2.3.2."
  endif
end subroutine Checkpoint

#else

subroutine Checkpoint(realization,steps,newtcum,icutcum, &
                               num_const_timesteps, &
                               num_newton_iterations,id)

  use Realization_module
  use Discretization_module
  use Option_module
  use Field_module
  
  use MPHASE_module

  implicit none

  type(realization_type) :: realization
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations
  PetscInt :: id, steps, newtcum, icutcum
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
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_WRITE, &
                             viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  bagsize = 80
  call PetscBagCreate(PETSC_COMM_WORLD, bagsize, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  ! Register variables that are passed into timestepper().
  call PetscBagRegisterInt(bag, header%num_newton_iterations, num_newton_iterations, &
                           "num_newton_iterations", &
                           "Number of Newton iterations in last SNES solve", ierr)
  call PetscBagRegisterInt(bag, header%plot_number, output_option%plot_number, &
                           "plot_number","plot_number", ierr)
  call PetscBagRegisterInt(bag, header%num_const_timesteps, num_const_timesteps, &
                           "num_const_timesteps","num_const_timesteps", ierr)
  
  ! Register relevant components of the stepper.
  call PetscBagRegisterReal(bag, header%flow_time, option%flow_time, "time", &
                            "Flow Simulation time (seconds)", ierr)
  call PetscBagRegisterReal(bag, header%flow_dt, option%flow_dt, "dt", &
                            "Current size of flow timestep (seconds)", ierr)
  call PetscBagRegisterReal(bag, header%tran_time, option%tran_time, "time", &
                            "Transport Simulation time (seconds)", ierr)
  call PetscBagRegisterReal(bag, header%tran_dt, option%tran_dt, "dt", &
                            "Current size of transport timestep (years)", ierr)
  call PetscBagRegisterInt(bag, header%steps, steps, "steps", &
                            "Total number of flow steps taken", ierr)
  call PetscBagRegisterInt(bag, header%newtcum, newtcum, "newtcum", &
                            "Total number of Newton steps taken", ierr)
  call PetscBagRegisterInt(bag, header%icutcum, icutcum, "icutcum", &
                            "Total number of time step cuts", ierr)

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  ! grid%flow_xx is the vector into which all of the primary variables are 
  ! packed for the SNESSolve().
  call VecView(field%flow_xx, viewer, ierr)

  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL)
  ! If we are running with multiple phases, we need to dump the vector 
  ! that indicates what phases are present, as well as the 'var' vector 
  ! that holds variables derived from the primary ones via the translator.
  select case(option%iflowmode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
      call DiscretizationLocalToGlobal(realization%discretization,field%iphas_loc, &
                             global_vec,ONEDOF)
      call VecView(global_vec, viewer, ierr)
      if (option%iflowmode == MPH_MODE) then
      ! get vardof vec from mphase
        call MphaseCheckpointWrite(realization%discretization,viewer)
      endif
    case default
  end select 

  ! Porosity and permeability.
  ! (We only write diagonal terms of the permeability tensor for now, 
  ! since we have yet to add the full-tensor formulation.)
  call DiscretizationLocalToGlobal(discretization,field%porosity_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call DiscretizationLocalToGlobal(discretization,field%perm_xx_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call DiscretizationLocalToGlobal(discretization,field%perm_yy_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call DiscretizationLocalToGlobal(discretization,field%perm_zz_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)

  call VecDestroy(global_vec,ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

  if(option%myrank == 0) write(*, '(" --> Dump checkpoint file: ", a16)') trim(fname)
  call PetscGetTime(tend,ierr) 
  if (realization%option%myrank == 0) &
    print *, '      Seconds to write to checkpoint file: ', (tend-tstart)

end subroutine Checkpoint

#endif

#if (PETSC_VERSION_RELEASE == 1 && PETSC_VERSION_SUBMINOR < 3)

subroutine Restart(realization,steps,newtcum,icutcum, &
                            num_const_timesteps, &
                            num_newton_iterations)
  use Realization_module
  use Option_module

  character(len=MAXSTRINGLENGTH) :: fname
  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations
  PetscInt :: steps, newtcum, icutcum

  if(realization%option%myrank == 0) then
    print *, "Warning: Restart() not supported with PETSc 2.3.2."
  endif
end subroutine Restart

#else

subroutine Restart(realization,steps,newtcum,icutcum, &
                            num_const_timesteps, &
                            num_newton_iterations)
  use Realization_module
  use Discretization_module
  use Option_module
  use Field_module
  
  use MPHASE_module

  implicit none

  type(realization_type) :: realization
  PetscInt :: num_const_timesteps
  PetscInt :: num_newton_iterations
  PetscInt :: steps, newtcum, icutcum

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
  
  field => realization%field
  option => realization%option
  discretization => realization%discretization
  output_option => realization%output_option
  
  call PetscGetTime(tstart,ierr)   
  if (option%myrank == 0) print *,'--> Open checkpoint file: ', trim(option%restart_file)
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, option%restart_file, FILE_MODE_READ, &
                             viewer, ierr)
 
  ! Get the header data.
  call PetscBagLoad(viewer, bag, ierr)
  call PetscBagGetData(bag, header, ierr)
  num_newton_iterations = header%num_newton_iterations
  output_option%plot_number = header%plot_number
  num_const_timesteps = header%num_const_timesteps
  option%flow_time = header%flow_time
  option%flow_dt = header%flow_dt
  option%tran_time = header%tran_time
  option%tran_dt = header%tran_dt
  steps = header%steps
  newtcum = header%newtcum
  icutcum = header%icutcum
  call PetscBagDestroy(bag, ierr)
  
  ! Load the PETSc vectors.
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL)

  call VecLoadIntoVector(viewer, field%flow_xx, ierr)
  call VecCopy(field%flow_xx, field%flow_yy, ierr)
  
  select case(option%iflowmode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
      call VecLoadIntoVector(viewer, global_vec, ierr)      
      call DiscretizationGlobalToLocal(discretization,global_vec,field%iphas_loc,ONEDOF)
      call VecCopy(field%iphas_loc, field%iphas_old_loc, ierr)
      call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_old_loc,ONEDOF)
      if (option%iflowmode == MPH_MODE) then
      ! set vardof vec in mphase
      endif
    case default
  end select
  
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call DiscretizationGlobalToLocal(discretization,global_vec,field%porosity_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call DiscretizationGlobalToLocal(discretization,global_vec,field%perm_xx_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call DiscretizationGlobalToLocal(discretization,global_vec,field%perm_yy_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call DiscretizationGlobalToLocal(discretization,global_vec,field%perm_zz_loc,ONEDOF)
  
  call VecDestroy(global_vec,ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)
  call PetscGetTime(tend,ierr) 
  if (realization%option%myrank == 0) &
    print *, '      Seconds to read checkpoint file: ', (tend-tstart)
  
end subroutine Restart

#endif

end module Checkpoint_module
