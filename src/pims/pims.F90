  program pflow
#include "include/finclude/petscsnes.h"
  use petscsnes
  use pflow_gridtype_module
  use pflow_grid_module
  use pflow_read_module
  use pflow_output_module
  use pflow_step
  
  implicit none

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  integer :: ierr, ihalcnt
  type(pflowGrid), target :: grid
  type(pflowGridParameters) :: gridparameters
  type(pflow_solver_context) :: pflowsolv
  type(time_stepping_context), target :: timestep
  
!  type(pflow_localpatch_info) :: gridpatch
    
  integer :: igeom
  integer*4 :: nx, ny, nz
  integer*4 :: npx, npy, npz
  integer :: nphase, ndof, idcdm, itable
  integer :: nspec,npricomp
  integer :: kplt, iplot, iflgcut, its, ntstep
  integer :: myid
  integer*4 :: steps

  type(pflow_localpatch_info),pointer :: locpat

! Initialize Startup Time
 ! call PetscGetCPUTime(timex(1), ierr)
 ! call PetscGetTime(timex_wall(1), ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call MPI_Comm_rank(PETSC_COMM_WORLD, myid, ierr)
 ! call MPI_Comm_size(PETSC_COMM_WORLD, grid%commsize, ierr)
 
 ! For current pims petsc version, every processor only have one patch on one level
  print *,' Begin reading'
 
   call pflow_read_gridsize("pflow.in", igeom, nx, ny, nz, npx, npy, npz, &
  nphase,  ierr)
  
  ntstep=1
  iflgcut = 0

  if(ierr /= 0) then
    if(myid == 0) then
      print *, "Error reading gridsize--exiting."
    endif 
    call PetscFinalize(ierr)
    stop
  endif

  gridparameters%grid => grid
  gridparameters%timestep => timestep
  gridparameters%igeom = igeom
  gridparameters%nx    = nx
  gridparameters%nx    = ny
  gridparameters%nx    = nz
  gridparameters%npx   = npx
  gridparameters%npy   = npy
  gridparameters%npz   = npz
  gridparameters%nphase= nphase

! set up structure constructor
! npx = PETSC_DECIDE; npy = PETSC_DECIDE; npz = PETSC_DECIDE
  call pflowGrid_new(grid, timestep,igeom, nx, ny, nz, npx, npy, npz, nphase)
  !  call pflowGrid_newpm(gridparameters)

  locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr

  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

  call pflowGrid_setup(grid, timestep,locpat,"pflow.in")

  !-----------------------------------------------------------------------
  ! Set up the Jacobian matrix.  We do this here instead of in 
  ! pflowgrid_new() because we may have to parse the input file to 
  ! determine how we want to do the Jacobian (matrix vs. matrix-free, for
  ! example).
  !-----------------------------------------------------------------------
  
  call pflowGrid_Setup_SNES(grid, pflowsolv )

  CHKMEMQ
  kplt = 0
  iplot = 1

! call VecView(grid%porosity,PETSC_VIEWER_STDOUT_WORLD,ierr)

   print *,' Begin output'
      call pflow_var_output(grid,timestep,kplt,iplot)
    print *,' End output'
  	  
  if(myid == 0) print *, ""
  if(myid == 0) print *, ""

  kplt = 1
  iplot = 0
  ihalcnt = 0
  grid%isrc1 = 2

  do steps = 1, grid%stepmax

    call pflowGrid_step(grid,pflowsolv,timestep,ntstep,kplt,iplot,iflgcut,ihalcnt,its)

!   update field variables
    call pflowGrid_update(grid)
  ! call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   comment out for now-need in pflotran (no longer needed-pcl)
!   call pflowgrid_setvel(grid, grid%vl, grid%vlbc, ibndtyp)
      call pflow_var_output(grid,timestep,kplt,iplot)


    if (iflgcut == 0) call pflowgrid_update_dt(grid,timestep,its)

    if (kplt .gt. timestep%kplot) exit
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
