! The header for the binary checkpoint files.
! This needs to go in its own module for some reason, because otherwise 
! some compilers won't like the INTERFACE block needed to allow us to 
! use the PetscBagGetData() routine.
! RTM: This is pretty makeshift.  We need to think about what should 
! go into this header and how it should be organized.

module pflow_chkptheader
  implicit none
  private
  type, public :: pflowChkPtHeader
    real*8 :: time
    real*8 :: dt
    integer*8 :: flowsteps
    integer*8 :: newtcum
    integer*8 :: icutcum
    integer*8 :: timestep_cut_flag
    integer*8 :: num_timestep_cuts
    integer*8 :: num_newton_iterations
    integer*8 :: plot_number
  end type pflowChkPtHeader
end module pflow_chkptheader

module pflow_checkpoint
  use pflow_chkptheader

  implicit none
  
  private

  public :: pflowGridCheckpoint, pflowGridRestart

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
      use pflow_chkptheader
      PetscBag bag
      type(pflowChkPtHeader), pointer :: ctx
      PetscErrorCode ierr
    End Subroutine
  End Interface PetscBagGetData

#endif

contains

#if (PETSC_VERSION_RELEASE == 1 && PETSC_VERSION_SUBMINOR < 3)

subroutine pflowGridCheckpoint(realization,flowsteps,newtcum,icutcum, &
                               timestep_cut_flag,num_timestep_cuts, &
                               num_newton_iterations,id)
  use Realization_module
  use Option_module

  implicit none
  
  type(realization_type) :: realization
  logical :: timestep_cut_flag
  PetscInt :: num_timestep_cuts, num_newton_iterations
  PetscInt :: id, flowsteps, newtcum, icutcum

  if(realization%option%myrank == 0) then
    print *, "Warning: pflowGridCheckpoint() not supported with PETSc 2.3.2."
  endif
end subroutine pflowGridCheckpoint

#else

subroutine pflowGridCheckpoint(realization,flowsteps,newtcum,icutcum, &
                               timestep_cut_flag,num_timestep_cuts, &
                               num_newton_iterations,id)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none

  type(realization_type) :: realization
  logical :: timestep_cut_flag
  PetscInt :: num_timestep_cuts, num_newton_iterations
  PetscInt :: id, flowsteps, newtcum, icutcum
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
  type(pflowChkPtHeader), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend  
  
  Vec :: global_vec, global_var
  PetscInt :: int_flag
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  
  field => realization%field
  option => realization%option
  grid => realization%grid
  output_option => realization%output_option

  ! Open the checkpoint file.
  call PetscGetTime(tstart,ierr)   
  if (id < 0) then
    fname = 'restart.chk'
  else if (id < 10) then
    write(fname, '(a9,i1)') 'pflow.chk', id
  else if (id < 100) then
    write(fname, '(a9,i2)') 'pflow.chk', id
  else if (id < 1000) then
    write(fname, '(a9,i3)') 'pflow.chk', id
  else if (id < 10000) then
    write(fname, '(a9,i4)') 'pflow.chk', id
  else if (id < 100000) then
    write(fname, '(a9,i5)') 'pflow.chk', id
  else if (id < 1000000) then
    write(fname, '(a9,i6)') 'pflow.chk', id
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
  bagsize = 72
  call PetscBagCreate(PETSC_COMM_WORLD, bagsize, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  ! Register variables that are passed into pflowGrid_step().
  call PetscBagRegisterInt(bag, header%num_newton_iterations, num_newton_iterations, &
                           "num_newton_iterations", &
                           "Number of Newton iterations in last SNES solve", ierr)
  call PetscBagRegisterInt(bag, header%plot_number, output_option%plot_number, &
                           "plot_number","plot_number", ierr)
  int_flag = 0
  if (timestep_cut_flag) int_flag = 1
  call PetscBagRegisterInt(bag, header%timestep_cut_flag, int_flag, &
                           "timestep_cut_flag","timestep_cut_flag", ierr)
  call PetscBagRegisterInt(bag, header%num_timestep_cuts, num_timestep_cuts, &
                           "num_timestep_cuts","num_timestep_cuts", ierr)
  
  ! Register relevant components of the pflowGrid.
  call PetscBagRegisterReal(bag, header%time, option%time, "time", &
                            "Simulation time (years)", ierr)
  call PetscBagRegisterReal(bag, header%dt, option%dt, "dt", &
                            "Current size of timestep (years)", ierr)
  call PetscBagRegisterInt(bag, header%flowsteps, flowsteps, "flowsteps", &
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

  ! grid%xx is the vector into which all of the primary variables are 
  ! packed for the SNESSolve().
  call VecView(field%xx, viewer, ierr)

  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)
  ! If we are running with multiple phases, we need to dump the vector 
  ! that indicates what phases are present, as well as the 'var' vector 
  ! that holds variables derived from the primary ones via the translator.
  select case(option%imode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
      call GridLocalToGlobal(grid,field%iphas_loc,global_vec,ONEDOF)
      call VecView(global_vec, viewer, ierr)
#ifdef RICHARDS_ANALYTICAL
      if (option%imode /= RICHARDS_MODE .and. &
          option%imode /= RICHARDS_LITE_MODE) then
#endif
        call GridCreateVector(grid,VARDOF,global_var,GLOBAL)
        call GridLocalToGlobal(grid,field%var_loc,global_var,VARDOF)
        call VecView(global_var, viewer, ierr)
        call VecDestroy(global_var,ierr)
#ifdef RICHARDS_ANALYTICAL
      endif
#endif
    case default
      call VecView(field%hh, viewer, ierr)
      call VecView(field%ddensity, viewer, ierr)
  end select 

  ! solid volume fraction
  if (option%rk > 0.d0) then
    call VecView(field%phis, viewer, ierr)
  endif

  ! Porosity and permeability.
  ! (We only write diagonal terms of the permeability tensor for now, 
  ! since we have yet to add the full-tensor formulation.)
  call GridLocalToGlobal(grid,field%porosity_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call GridLocalToGlobal(grid,field%perm_xx_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call GridLocalToGlobal(grid,field%perm_yy_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)
  call GridLocalToGlobal(grid,field%perm_zz_loc,global_vec,ONEDOF)
  call VecView(global_vec,viewer,ierr)

  call VecDestroy(global_vec,ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

  if(option%myrank == 0) write(*, '(" --> Dump checkpoint file: ", a16)') trim(fname)
  call PetscGetTime(tend,ierr) 
  if (realization%option%myrank == 0) &
    print *, '      Seconds to write to checkpoint file: ', (tend-tstart)

end subroutine pflowGridCheckpoint

#endif

#if (PETSC_VERSION_RELEASE == 1 && PETSC_VERSION_SUBMINOR < 3)

subroutine pflowGridRestart(realization,flowsteps,newtcum,icutcum, &
                            timestep_cut_flag,num_timestep_cuts, &
                            num_newton_iterations)
  use Realization_module
  use Option_module

  character(len=MAXSTRINGLENGTH) :: fname
  type(realization_type) :: realization
  type(stepper_type) :: stepper
  logical :: timestep_cut_flag
  PetscInt :: num_timestep_cuts, num_newton_iterations
  PetscInt :: flowsteps, newtcum, icutcum

  if(realization%option%myrank == 0) then
    print *, "Warning: pflowGridRestart() not supported with PETSc 2.3.2."
  endif
end subroutine pflowGridRestart

#else

subroutine pflowGridRestart(realization,flowsteps,newtcum,icutcum, &
                            timestep_cut_flag,num_timestep_cuts, &
                            num_newton_iterations)
  use Realization_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none

  type(realization_type) :: realization
  logical :: timestep_cut_flag
  PetscInt :: num_timestep_cuts, num_newton_iterations
  PetscInt :: flowsteps, newtcum, icutcum

  PetscViewer viewer
  PetscBag bag
  type(pflowChkPtHeader), pointer :: header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend

  Vec :: global_vec, global_var
  PetscInt :: int_flag
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  
  field => realization%field
  option => realization%option
  grid => realization%grid
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
  timestep_cut_flag = header%timestep_cut_flag
  num_timestep_cuts = header%num_timestep_cuts
  option%time = header%time
  option%dt = header%dt
  flowsteps = header%flowsteps
  newtcum = header%newtcum
  icutcum = header%icutcum
  call PetscBagDestroy(bag, ierr)
  
  ! Load the PETSc vectors.
  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)

  call VecLoadIntoVector(viewer, field%xx, ierr)
  call VecCopy(field%xx, field%yy, ierr)
  
  select case(option%imode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
      call VecLoadIntoVector(viewer, global_vec, ierr)      
      call GridGlobalToLocal(grid,global_vec,field%iphas_loc,ONEDOF)
      call VecCopy(field%iphas_loc, field%iphas_old_loc, ierr)
      call GridLocalToLocal(grid,field%iphas_loc,field%iphas_old_loc,ONEDOF)
#ifdef RICHARDS_ANALYTICAL
      if (option%imode /= RICHARDS_MODE .and. &
          option%imode /= RICHARDS_LITE_MODE) then
#endif      
        call GridCreateVector(grid,VARDOF,global_var,GLOBAL)
        call VecLoadIntoVector(viewer, global_var, ierr)
        call GridGlobalToLocal(grid,global_var,field%var_loc,VARDOF)
        call VecDestroy(global_var,ierr)
#ifdef RICHARDS_ANALYTICAL
      endif
#endif
    case default
      call VecLoadIntoVector(viewer, field%hh, ierr)
      call VecCopy(field%hh, field%h, ierr)
      call VecLoadIntoVector(viewer, field%ddensity, ierr)
      call VecCopy(field%ddensity, field%density, ierr)
  end select
  
  if (option%rk > 0.d0) then
    call VecLoadIntoVector(viewer, field%phis, ierr)
  endif
  
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call GridGlobalToLocal(grid,global_vec,field%porosity_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call GridGlobalToLocal(grid,global_vec,field%perm_xx_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call GridGlobalToLocal(grid,global_vec,field%perm_yy_loc,ONEDOF)
  call VecLoadIntoVector(viewer, global_vec, ierr)
  call GridGlobalToLocal(grid,global_vec,field%perm_zz_loc,ONEDOF)
  
  call VecDestroy(global_vec,ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)
  call PetscGetTime(tend,ierr) 
  if (realization%option%myrank == 0) &
    print *, '      Seconds to read checkpoint file: ', (tend-tstart)
  
end subroutine pflowGridRestart

#endif

#if 0
subroutine pflowGridTHCBinaryOut(grid, kplt)

  implicit none

  type(pflowGrid), intent(inout) :: grid
  PetscInt, intent(inout) :: kplt

  character(len=MAXSTRINGLENGTH) :: fname
  PetscViewer viewer
  PetscInt :: ierr

  ! Open the output file.
  write(fname, '(a9,i2)') 'pflow.chk', kplt
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_WRITE, &
                        viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  ! RTM: I need to write this code!
  ! Members of 'grid' that have to be dumped: t, dt, flowsteps, kplot
  ! What else?

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  call VecView(grid%conc, viewer, ierr)
  call VecView(grid%vl, viewer, ierr)
  call VecView(grid%vg, viewer, ierr)
    ! RTM: When do we actually need to dump the gas velocities vg?
  call VecView(grid%pressure, viewer, ierr)
  call VecView(grid%temp, viewer, ierr)
  call VecView(grid%sat, viewer, ierr)

  ! primary variables
  if(grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                .or. grid%use_vadose == PETSC_TRUE &
                                .or. grid%use_flash == PETSC_TRUE &
                                 .or. grid%use_richards == PETSC_TRUE) then
    call VecView(grid%xmol, viewer, ierr)
  endif  
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_richards == PETSC_TRUE ) then
    call VecView(grid%iphas, viewer, ierr)
  endif  

  ! solid volume fraction
  if (grid%rk > 0.d0) then
    call VecView(grid%phis, viewer, ierr)
  endif

  ! Porosity and permeability.
  ! (We only write diagonal terms of the permeability tensor for now, 
  ! since we have yet to add the full-tensor formulation.)
  call VecView(grid%porosity, viewer, ierr)
  call VecView(grid%perm_xx, viewer, ierr)
  call VecView(grid%perm_yy, viewer, ierr)
  call VecView(grid%perm_zz, viewer, ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

end subroutine pflowGridTHCBinaryOut
#endif
end module pflow_checkpoint
