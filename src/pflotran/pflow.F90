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
  use Output_module
  use span_wagner_module
  
  use Simulation_module
  use Solution_module
  use Solver_module
  use Grid_module
  use Option_module
  use Init_module
  use Timestepper_module
  
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
  integer :: kplt, iplot, iflgcut, its, ntstep
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
  
  type(simulation_type), pointer :: simulation
  type(solver_type), pointer :: solver
  type(solution_type), pointer :: solution
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  simulation => createSimulation()
  simulation%stepper => createTimestepper()
  solution => simulation%solution
  option => solution%option

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
  call MPI_Comm_rank(PETSC_COMM_WORLD,option%myrank, ierr)

! Starting PFLOW setup -- push this onto the log stack.
  call PetscLogStagePush(stage(1), ierr)
  
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflowin", &
    pflowin, option_found, ierr)
  if(option_found /= PETSC_TRUE) pflowin = "pflow.in"
  
  ntstep=1
  iflgcut = 0

  if(ierr /= 0) then
    if(option%myrank == 0) then
      print *, "Error reading gridsize--exiting."
    endif 
    call PetscFinalize(ierr)
    stop
  endif

  if(option%use_mph == PETSC_TRUE .or. &
     option%use_owg == PETSC_TRUE .or. &
     option%use_flash == PETSC_TRUE) &
    call initialize_span_wagner(option%itable,option%myrank)

  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

  call initPFLOW(simulation,pflowin)
  allocate(dxdt(1:option%ndof))

  kplt = 0
  iplot = 1

! Done with PFLOW setup.
  call PetscLogStagePop(ierr)

! Now we do some initial output, so push this onto the log stack.
  call PetscLogStagePush(stage(2), ierr)
! call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)

  if(option%use_owg/=PETSC_TRUE) then
    call Output(solution,kplt,iplot)
  else
 !   call pflow_var_output(grid,kplt,iplot)
  endif

  if(option%myrank == 0) print *, ""
  if(option%myrank == 0) print *, ""

  kplt = 1
  iplot = 0
  ihalcnt = 0
  option%isrc1 = 2

  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', restartfile, &
                             restartflag, ierr)
#if 0                             
  if(restartflag == PETSC_TRUE) then
    call pflowGridRestart(grid, restartfile, ntstep, kplt, iplot, iflgcut, &
                          ihalcnt,its)
    call pflowGridInitAccum(grid)
    option%dt_max = grid%dtstep(ntstep)
  endif
#endif  
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-chkptfreq', chkptfreq, &
                          chkptflag, ierr)
                             
  call PetscLogStagePop(ierr)

#if 0
 if (grid%iread_init==2 .and. restartflag == PETSC_FALSE)then
      call Read_init_field(grid, kplt)
      call pflowGridInitAccum(grid)
      print *, 'Restart from ASCII file: pflow_init.dat'
    !  print *, 't, dt, kplt :: ',grid%t,grid%dt,kplt,grid%flowsteps,iplot,ntstep
  endif
#endif  
           
  do steps = option%flowsteps+1, option%stepmax

!    call pflowGrid_step(solution,ntstep,kplt,iplot,iflgcut,ihalcnt,its)
    call step(solution,solver,ntstep,kplt,iplot,iflgcut,ihalcnt,its)

!   update field variables
!    call pflowGrid_update(solution)
    call update(solution)

    dt_cur = option%dt 
   
    
    if(option%use_thc == PETSC_TRUE)then
      dxdt(1)=option%dpmax/dt_cur
      dxdt(2)=option%dtmpmax/dt_cur
      dxdt(3)=option%dcmax/dt_cur
    endif  

    call PetscLogStagePush(stage(2), ierr)
    if(option%use_owg/=PETSC_TRUE) then
      call Output(solution,kplt,iplot)
     ! print *,'XX ::...........'; call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else
 !     call pflow_var_output(grid,kplt,iplot)
    endif
    call PetscLogStagePop(ierr)
  
  
!    if (iflgcut == 0) call pflowgrid_update_dt(solution,its)
    if (iflgcut == 0) call updateDT(option,its)

#if 0
    call PetscLogStagePush(stage(2), ierr)
    if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) == 0) then
      call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                               its, steps)
    endif
    call PetscLogStagePop(ierr)
#endif
    
    ista=0
    if(option%use_thc == PETSC_TRUE)then
      do idx = 1, option%ndof
        if(dxdt(idx) < option%steady_eps(idx)) ista=ista+1
      enddo 
      
      if(ista >= option%ndof)then
        kplt=option%kplot; iplot=1     
      endif
    endif   
!   comment out for now-need in pflotran (no longer needed-pcl)
!   call pflowgrid_setvel(grid, grid%vl, grid%vlbc, ibndtyp)



!    call PetscLogStagePush(stage(2), ierr)
!    if(option%use_owg/=PETSC_TRUE) then
!      call pflow_output_new(grid,kplt,iplot)
     ! print *,'XX ::...........'; call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    else
!      call pflow_var_output(grid,kplt,iplot)
!    endif
!    call PetscLogStagePop(ierr)
  
  
!    if (iflgcut == 0) call pflowgrid_update_dt(grid,its)
      
    
    if (kplt .gt. option%kplot) exit
  enddo

#if 0
  if(chkptflag == PETSC_TRUE .and. mod(steps, chkptfreq) /= 0) then
    call pflowGridCheckpoint(grid, ntstep, kplt, iplot, iflgcut, ihalcnt, &
                             its, steps)
  endif
#endif  

#if 0
! Clean things up.
  call PetscLogStagePush(stage(3), ierr)
  call pflowgrid_destroy(grid)
  call PetscLogStagePop(ierr)
#endif
  call PetscFinalize (ierr)

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  call PetscGetTime(timex_wall(2), ierr)
  if (option%myrank == 0) then
!       write(*,'("steps: ",13x,i6,/, &
!   &   "total newton iters: ",i6,/, &
!   &   "total PETSc gmres:  ",i6,/, &
!   &   "total cuts: ",11x,i3)') steps,newtcum,itpetscum,icutcum
        
        write(*,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
        steps-1,option%newtcum,option%icutcum

        write(IUNIT2,'(/," PFLOW steps = ",i6," newton = ",i6," cuts = ",i6)') &
        steps-1,option%newtcum,option%icutcum

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
  if (option%ibrkcrv > 0) close(IUNIT4)

  end program pflow
