! The header for the binary checkpoint files.
! This needs to go in its own module for some reason, because otherwise 
! some compilers won't like the INTERFACE block needed to allow us to 
! use the PetscBagGetData() routine.
! RTM: This is pretty makeshift.  We need to think about what should 
! go into this header and how it should be organized.
module pflow_chkptheader
  type, public :: pflowChkPtHeader
    ! Paramters passed into pflowGrid_step().
    integer :: ntstep, kplt, iplot, iflgcut, ihalcnt, its

    ! Necessary components of the pflowGrid (excluding PETSc Vecs).
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
#if (PETSC_VERSION_RELEASE == 0)
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

#if (PETSC_VERSION_RELEASE == 1)

subroutine pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                               its, id)
  use pflow_gridtype_module

  type(pflowGrid), intent(inout) :: grid
  integer, intent(in) :: ntstep, kplt, iplot, iflgcut, ihalcnt, its
  integer, intent(in) :: id

  print *, "Warning: pflowGridCheckpoint() not supported with PETSc 2.3.2."
end subroutine pflowGridCheckpoint

#else

subroutine pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                               its, id)

  use pflow_gridtype_module
  use TTPHASE_module

  implicit none

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(in) :: ntstep, kplt, iplot, iflgcut, ihalcnt, its
  integer, intent(in) :: id

  character(len=256) :: fname
  PetscViewer viewer
  PetscBag bag
  type(pflowChkPtHeader), pointer :: header
  integer ierr

  ! Open the checkpoint file.
  write(fname, '(a10,i6.6)') 'pflow.chk.', id
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_WRITE, &
                             viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------

  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  call PetscBagCreate(PETSC_COMM_WORLD, 88, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  ! Register variables that are passed into pflowGrid_step().
  call PetscBagRegisterInt(bag, header%ntstep, ntstep, "ntstep", &
                           "ntstep", ierr)
  call PetscBagRegisterInt(bag, header%kplt, kplt, "kplt", &
                           "kplt", ierr)
  call PetscBagRegisterInt(bag, header%iplot, iplot, "iplot", &
                           "iplot", ierr)
  call PetscBagRegisterInt(bag, header%iflgcut, iflgcut, "iflgcut", &
                           "iflgcut", ierr)
  call PetscBagRegisterInt(bag, header%ihalcnt, ihalcnt, "ihalcnt", &
                           "ihalcnt", ierr)
  call PetscBagRegisterInt(bag, header%its, its, "its", &
                           "its", ierr)
  
  ! Register relevant components of the pflowGrid.
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

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  ! grid%xx is the vector into which all of the primary variables are 
  ! packed for the SNESSolve().
  call VecView(grid%xx, viewer, ierr)
  call VecView(grid%hh, viewer, ierr)
  call VecView(grid%ddensity, viewer, ierr)

  ! If we are running with multiple phases, we need to dump the vector 
  ! that indicates what phases are present, as well as the 'var' vector 
  ! that holds variables derived from the primary ones via the translator.
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE .or. &
     grid%use_flash == PETSC_TRUE .or. grid%use_2ph == PETSC_TRUE &
      .or. grid%use_richard == PETSC_TRUE ) then
    call VecView(grid%iphas, viewer, ierr)
    call VecView(grid%var, viewer, ierr)
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

  write(*, '(a23, a16)') "Dumped checkpoint file ", fname

end subroutine pflowGridCheckpoint

#endif

#if (PETSC_VERSION_RELEASE == 1)

subroutine pflowGridRestart(grid, fname, ntstep, kplt, iplot, iflgcut, &
                            ihalcnt, its)
  use pflow_gridtype_module

  type(pflowGrid), intent(inout) :: grid
  character(len=256) :: fname
  integer, intent(out) :: ntstep, kplt, iplot, iflgcut, ihalcnt, its

  print *, "Warning: pflowGridRestart() not supported with PETSc 2.3.2."
end subroutine pflowGridRestart

#else

subroutine pflowGridRestart(grid, fname, ntstep, kplt, iplot, iflgcut, &
                            ihalcnt, its)
  use pflow_gridtype_module
  use TTPHASE_module

  implicit none

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  character(len=256) :: fname
  integer, intent(out) :: ntstep, kplt, iplot, iflgcut, ihalcnt, its

  PetscViewer viewer
  PetscBag bag
  type(pflowChkPtHeader), pointer :: header
  integer ierr

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_READ, &
                             viewer, ierr)
 
  ! Get the header data.
  call PetscBagLoad(viewer, bag, ierr)
  call PetscBagGetData(bag, header, ierr)
  ntstep = header%ntstep
  kplt = header%kplt
  iplot = header%iplot
  iflgcut = header%iflgcut
  ihalcnt = header%ihalcnt
  its = header%its
  grid%t = header%t
  grid%dt = header%dt
  grid%flowsteps = header%flowsteps
  grid%kplot = header%kplot
  grid%newtcum = header%newtcum
  call PetscBagDestroy(bag, ierr)
  
  ! Load the PETSc vectors.
  call VecLoadIntoVector(viewer, grid%xx, ierr)
  call VecCopy(grid%xx, grid%yy, ierr)
  
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE .or. &
     grid%use_flash == PETSC_TRUE .or. grid%use_2ph == PETSC_TRUE .or. &
     grid%use_richard == PETSC_TRUE ) then
    call VecLoadIntoVector(viewer, grid%iphas, ierr)
    call VecCopy(grid%iphas, grid%iphas_old, ierr)
    call VecLoadIntoVector(viewer, grid%var, ierr)
  else
    call VecLoadIntoVector(viewer, grid%hh, ierr)
    call VecCopy(grid%hh, grid%h, ierr)
    call VecLoadIntoVector(viewer, grid%ddensity, ierr)
    call VecCopy(grid%ddensity, grid%density, ierr)
  endif
  
  if (grid%rk > 0.d0) then
    call VecLoadIntoVector(viewer, grid%phis, ierr)
  endif
  call VecLoadIntoVector(viewer, grid%porosity, ierr)
  call VecLoadIntoVector(viewer, grid%perm_xx, ierr)
  call VecLoadIntoVector(viewer, grid%perm_yy, ierr)
  call VecLoadIntoVector(viewer, grid%perm_zz, ierr)

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)
end subroutine pflowGridRestart

#endif


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

  ! Open the output file.
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
                                .or. grid%use_flash == PETSC_TRUE &
                                 .or. grid%use_richard == PETSC_TRUE) then
    call VecView(grid%xmol, viewer, ierr)
  endif  
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_richard == PETSC_TRUE ) then
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
