!=======================================================================
! PFLOTRAN v1.0 LA-CC-06-093
!=======================================================================

! The software titled "PFLOTRAN v 1.0" has been assigned LA-CC-06-093. 
! The software is unclassified and does not contain Unclassified 
! Controlled Nuclear Information (UCNI).

! The software is under review to be released as open source, which 
! requires review and approval by the appropriate DOE Program Office. 
! If DOE declines to release this software as publicly available open 
! source, it would be subject to export control under Department of 
! Commerce regulations, and classified as ECCN 
! (Export Control Classification Number) EAR99. 
 
! If released as open source, the software will be publicly available and 
! not subject to export control. Until DOE approval is obtained, the 
! software should be treated as export controlled. DOE reserves the 
! right to release or deny release of the software as open source.

! Send all bug reports/questions/comments to:
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-6, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM
!=======================================================================
  program pflow
  use pflow_gridtype_module
  use pflow_grid_module
  use pflow_read_gridsize_module
  use pflow_output_module
  use pflow_output_module_new
  use span_wagner_module
  use pflow_checkpoint
  use readfield, only: Read_init_field
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  integer :: ierr, ihalcnt
  type(pflowGrid) :: grid
  integer :: igeom
  integer :: nx, ny, nz
  integer :: npx, npy, npz
  integer :: nphase, ndof, icouple, idcdm, itable, mcomp, mphas
  integer :: nspec,npricomp
  integer :: kplt, iplot, iflgcut, its, ntstep
  integer :: myid
  integer :: steps
  integer :: idx, ista=0
  PetscInt :: stage(10)
  PetscTruth :: restartflag
  character(len=MAXSTRINGLENGTH) :: restartfile
  PetscInt :: chkptfreq  
    ! The number of flow steps between checkpoints.
  PetscTruth :: chkptflag
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflowin
  real*8 dt_cur
  real*8, pointer :: dxdt(:)

! Initialize Startup Time
 ! call PetscGetCPUTime(timex(1), ierr)
 ! call PetscGetTime(timex_wall(1), ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

! Register stages for profiling.
! We identify three stages here: setup, output, and cleanup.
! Note that we do NOT identify a separate stage for the calculations -- 
! we include those as part of the default "Main" stage (so we need 
! to pop all other stages off when calculations are proceeding).
  call PetscLogStageRegister(stage(1), "PFLOW setup", ierr)
  call PetscLogStageRegister(stage(2), "PFLOW output", ierr)
  call PetscLogStageRegister(stage(3), "PFLOW cleanup", ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, myid, ierr)

! Starting PFLOW setup -- push this onto the log stack.
  call PetscLogStagePush(stage(1), ierr)
  
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflowin", &
    pflowin, option_found, ierr)
  if(option_found /= PETSC_TRUE) pflowin = "pflow.in"
  
  call pflow_read_gridsize(pflowin, igeom, nx, ny, nz, npx, npy, npz, &
  nphase, nspec, npricomp, ndof, idcdm, itable, mcomp, mphas,ierr)
  
  allocate(dxdt(1:ndof))
  
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
  grid = pflowGrid_new(igeom, nx, ny, nz, npx, npy, npz, nphase, nspec, &
         npricomp, ndof, icouple, idcdm, itable, mcomp, mphas)

  if(grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE) &
  call initialize_span_wagner(itable,grid%myrank)

  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

  call pflowGrid_setup(grid, pflowin)
  CHKMEMQ
  kplt = 0
  iplot = 1

! Done with PFLOW setup.
  call PetscLogStagePop(ierr)

! Now we do some initial output, so push this onto the log stack.
  call PetscLogStagePush(stage(2), ierr)
! call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)

  if(grid%use_owg/=PETSC_TRUE) then
    call pflow_output_new(grid,kplt,iplot)
  else
    call pflow_var_output(grid,kplt,iplot)
  endif

  if(myid == 0) print *, ""
  if(myid == 0) print *, ""

  kplt = 1
  iplot = 0
  ihalcnt = 0
  grid%isrc1 = 2

  

  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', restartfile, &
                             restartflag, ierr)
  if(restartflag == PETSC_TRUE) then
    call pflowGridRestart(grid, restartfile, ntstep, kplt, iplot, iflgcut, &
                          ihalcnt,its)
    call pflowGridInitAccum(grid)
  endif
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-chkptfreq', chkptfreq, &
                          chkptflag, ierr)
                             
  call PetscLogStagePop(ierr)

 if (grid%iread_init==2 .and. restartflag == PETSC_FALSE)then
      call Read_init_field(grid, kplt)
      call pflowGridInitAccum(grid)
      print *, 'Restart from ASCII file: pflow_init.dat'
      print *, 't, dt, kplt :: ', grid%t,grid%dt, kplt, grid%flowsteps 
  endif


           
  do steps = grid%flowsteps+1, grid%stepmax

    call pflowGrid_step(grid,ntstep,kplt,iplot,iflgcut,ihalcnt,its)

!   update field variables
    call pflowGrid_update(grid)

    dt_cur = grid%dt 
   
    
    if(grid%use_thc == PETSC_TRUE)then
      dxdt(1)=grid%dpmax/dt_cur
      dxdt(2)=grid%dtmpmax/dt_cur
      dxdt(3)=grid%dcmax/dt_cur
    endif  

    call PetscLogStagePush(stage(2), ierr)
    if(grid%use_owg/=PETSC_TRUE) then
      call pflow_output_new(grid,kplt,iplot)
     ! print *,'XX ::...........'; call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else
      call pflow_var_output(grid,kplt,iplot)
    endif
    call PetscLogStagePop(ierr)
  
  
    if (iflgcut == 0) call pflowgrid_update_dt(grid,its)


    call PetscLogStagePush(stage(2), ierr)
    if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) == 0) then
      call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                               its, steps)
    endif
    call PetscLogStagePop(ierr)
    
    ista=0
    if(grid%use_thc == PETSC_TRUE)then
      do idx = 1, grid%ndof
        if(dxdt(idx) < grid%steady_eps(idx)) ista=ista+1
      enddo 
      
      if(ista >= grid%ndof)then
        kplt=grid%kplot; iplot=1     
      endif
    endif   
!   comment out for now-need in pflotran (no longer needed-pcl)
!   call pflowgrid_setvel(grid, grid%vl, grid%vlbc, ibndtyp)



!    call PetscLogStagePush(stage(2), ierr)
!    if(grid%use_owg/=PETSC_TRUE) then
!      call pflow_output_new(grid,kplt,iplot)
     ! print *,'XX ::...........'; call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    else
!      call pflow_var_output(grid,kplt,iplot)
!    endif
!    call PetscLogStagePop(ierr)
  
  
!    if (iflgcut == 0) call pflowgrid_update_dt(grid,its)
      
    
    if (kplt .gt. grid%kplot) exit
  enddo

  if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) /= 0) then
    call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                             its, steps)
  endif

! Clean things up.
  call PetscLogStagePush(stage(3), ierr)
  call pflowgrid_destroy(grid)
  call PetscLogStagePop(ierr)

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
