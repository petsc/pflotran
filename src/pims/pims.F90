  program pflow
  use pflow_gridtype_module
  use pflow_grid_module
  use pflow_read_module
  use pflow_output_module
  use pflow_step
  
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscpc.h"
!#include "include/finclude/petscsles.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  integer :: ierr, ihalcnt
  type(pflowGrid) :: grid
  integer :: igeom
  integer*4 :: nx, ny, nz
  integer*4 :: npx, npy, npz
  integer :: nphase, ndof, icouple, idcdm, itable
  integer :: nspec,npricomp
  integer :: kplt, iplot, iflgcut, its, ntstep
  integer :: myid
  integer*4 :: steps

! Initialize Startup Time
 ! call PetscGetCPUTime(timex(1), ierr)
 ! call PetscGetTime(timex_wall(1), ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_Comm_rank(PETSC_COMM_WORLD, myid, ierr)

  call pflow_read_gridsize("pflow.in", igeom, nx, ny, nz, npx, npy, npz, &
  nphase,  ierr)
  
  icouple = 0
  ntstep=1
  iflgcut = 0

  if(ierr /= 0) then
    if(myid == 0) then
      print *, "Error reading gridsize--exiting."
    endif 
    call PetscFinalize(ierr)
    stop
  endif

! set up structure constructor
! npx = PETSC_DECIDE; npy = PETSC_DECIDE; npz = PETSC_DECIDE
  grid = pflowGrid_new(igeom, nx, ny, nz, npx, npy, npz, nphase)

 
  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

  call pflowGrid_setup(grid, "pflow.in")
  CHKMEMQ
  kplt = 0
  iplot = 1

! call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)

      call pflow_var_output(grid,kplt,iplot)
  
  	  
  if(myid == 0) print *, ""
  if(myid == 0) print *, ""

  kplt = 1
  iplot = 0
  ihalcnt = 0
  grid%isrc1 = 2

  do steps = 1, grid%stepmax

    call pflowGrid_step(grid,ntstep,kplt,iplot,iflgcut,ihalcnt,its)

!   update field variables
    call pflowGrid_update(grid)
  ! call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   comment out for now-need in pflotran (no longer needed-pcl)
!   call pflowgrid_setvel(grid, grid%vl, grid%vlbc, ibndtyp)
      call pflow_var_output(grid,kplt,iplot)


    if (iflgcut == 0) call pflowgrid_update_dt(grid,its)

    if (kplt .gt. grid%kplot) exit
  enddo

!  call pflowgrid_destroy(grid)

  call PetscFinalize (ierr)

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  call PetscGetTime(timex_wall(2), ierr)
  if (myid == 0) then
!       write(*,'("steps: ",13x,i6,/, &
!   &   "total newton iters: ",i6,/, &
!   &   "total PETSc gmres:  ",i6,/, &
!   &   "total cuts: ",11x,i3)') steps,newtcum,itpetscum,icutcum
        
        write(*,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
        steps-1,grid%newtcum,grid%icutcum

        write(IUNIT2,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
        steps-1,grid%newtcum,grid%icutcum

        write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0

        write(IUNIT2,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(IUNIT2,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0
  endif
      
  close(IUNIT2)
  if (grid%ibrkcrv > 0) close(IUNIT4)
  
  end program pflow
