module pflow_checkpoint

  private
  
  public :: pflowGridCheckpoint, pflowGridRestart

contains

subroutine pflowGridCheckpoint(grid, kplt)

  use pflow_gridtype_module
  use TTPHASE_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(inout) :: kplt

  character(len=256) :: fname
  PetscViewer viewer
  integer ierr

  ! Open the checkpoint file.
  write(fname, '(a10,i2)') 'pflow.chk.', kplt
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
                                .or. grid%use_flash == PETSC_TRUE) then
    call VecView(grid%xmol, viewer, ierr)
  endif  
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE) then
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

end subroutine pflow_checkpoint(grid)

module pflow_checkpoint
