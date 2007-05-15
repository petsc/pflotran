! The header for the binary checkpoint files.
! This needs to go in its own module for some reason, because otherwise 
! some compilers won't like the INTERFACE block needed to allow us to 
! use the PetscBagGetData() routine.
! RTM: This is pretty makeshift.  We need to think about what should 
! go into this header and how it should be organized.
module pflow_chkptheader
  type, public :: pflowChkPtHeader
    real*8 :: t
    real*8 :: dt
    integer :: flowsteps 
    integer :: kplot
    integer :: newtcum
  end type pflowChkPtHeader
end module pflow_chkptheader

module pflow_checkpoint
  use pflow_chkptheader

  private

  public :: pflowGridCheckpoint, pflowGridRestart

#include "include/finclude/petsc.h"
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
#include "include/finclude/petscbag.h"

  Interface PetscBagGetData
    Subroutine PetscBagGetData(bag,ctx,ierr)
      use pflow_chkptheader
      PetscBag bag
      type(pflowChkPtHeader), pointer :: ctx
      PetscErrorCode ierr
    End Subroutine
  End Interface PetscBagGetData

contains

subroutine pflowGridCheckpoint(grid, id)

  use pflow_gridtype_module
  use TTPHASE_module

  implicit none

#include "definitions.h"

  type(pflowGrid) :: grid
  integer, intent(in) :: id

  character(len=256) :: fname
  PetscViewer viewer
  PetscBag bag
  type(pflowChkPtHeader), pointer :: header
  integer ierr

  ! Open the checkpoint file.
  write(fname, '(a10,i2)') 'pflow.chk.', id
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_WRITE, &
                        viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  call PetscBagCreate(PETSC_COMM_WORLD,40,bag,ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)
  call PetscBagRegisterReal(bag, header%t, grid%t, "t", &
                            "Simulation time (years)", ierr)
  call PetscBagRegisterReal(bag, header%dt, grid%dt, "dt", &
                            "Current size of timestep (years)", ierr)
  call PetscBagRegisterInt(bag, header%flowsteps, grid%flowsteps, "flowsteps", &
                            "Total number of flow steps taken", ierr)
  call PetscBagRegisterInt(bag, header%kplot, grid%kplot, "kplot", &
                            "Printout steps", ierr)
  call PetscBagRegisterInt(bag, header%newtcum, grid%newtcum, "newtcum", &
                            "Total number of Newton steps taken", ierr)
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  ! grid%xx is the vector into which all of the primary variables are 
  ! packed for the SNESSolve().
  call VecView(grid%xx, viewer, ierr)

  ! If we are running with multiple phases, we need to dump the vector 
  ! that indicates what phases are present. 
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE .or. &
     grid%use_flash == PETSC_TRUE .or. grid%use_2ph == PETSC_TRUE) then
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

end subroutine pflowGridCheckpoint


subroutine pflowGridTHCBinaryOut(grid, kplt)

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

end subroutine pflowGridTHCBinaryOut

end module pflow_checkpoint
