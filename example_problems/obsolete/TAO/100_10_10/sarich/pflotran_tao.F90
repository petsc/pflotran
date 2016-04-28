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
  program pflotran
  implicit none
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/tao_solver.h"
      
  
  external FormFunction
      
  TAO_SOLVER :: tao
  TAO_APPLICATION :: taoapp
  PetscInt :: reason

  Vec :: parameters,xl,xu
  PetscInt :: par_indices(6)
  PetscScalar :: par_array(6)
  PetscInt :: ierr
  PetscInt :: MAXTIMESTEPS
  parameter (MAXTIMESTEPS=3000)
  PetscScalar :: f

  PetscScalar :: reference_pressure(MAXTIMESTEPS,8)
  PetscScalar :: reference_times(MAXTIMESTEPS)
  integer reference_size      
  logical usetransform
  common /appctx/ usetransform,reference_size,reference_pressure,reference_times

  data par_indices/0,1,2,3,4,5/
      
  reason=0
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CHKERRQ(ierr)
  call TaoInitialize(PETSC_NULL_CHARACTER, ierr)      
  CHKERRQ(ierr)

! Check for transform
      
  call PetscOptionsHasName(TAO_NULL_CHARACTER,'-logtransform',usetransform,ierr)
  CHKERRQ(ierr)

  if (usetransform) then
    par_array(1)=-12d0
    par_array(2)=-12d0
    par_array(3)=-12d0
    par_array(4)=-12d0
    par_array(5)=-12d0
    par_array(6)=-12d0
  else
    par_array(1)=1.01d-13
    par_array(2)=1.01d-13
    par_array(3)=1.01d-13
    par_array(4)=0.99d-14
    par_array(5)=1.01d-14
    par_array(6)=1.01d-14
  endif
  call VecCreateSeq(PETSC_COMM_SELF,6,parameters,ierr)
  CHKERRQ(ierr)
  call VecSetValues(parameters,6,par_indices,par_array,INSERT_VALUES,ierr)
  CHKERRQ(ierr)
  call VecAssemblyBegin(parameters,ierr)
  CHKERRQ(ierr)
  call VecAssemblyEnd(parameters,ierr)
  CHKERRQ(ierr)

  call VecDuplicate(parameters,xl,ierr)
  CHKERRQ(ierr)
  call VecSet(xl,1.0d-20,ierr)
  CHKERRQ(ierr)

  call VecDuplicate(parameters,xu,ierr)
  CHKERRQ(ierr)
  call VecSet(xu,1.0d-10,ierr)
  CHKERRQ(ierr)


  call TaoCreate(PETSC_COMM_SELF,'tao_nm',tao,ierr)
  CHKERRQ(ierr)

  if (usetransform .eqv. PETSC_TRUE) then
    call PetscOptionsSetValue('-tao_nm_lamda','1d0',ierr)
    CHKERRQ(ierr)
  else
    print *,'set tao_nm_lamda'
!    call PetscOptionsSetValue('-tao_nm_lamda','1d-13',ierr)
!    CHKERRQ(ierr)
  endif	
  call TaoApplicationCreate(PETSC_COMM_SELF,taoapp,ierr)
  CHKERRQ(ierr)
  call TaoAppSetInitialSolutionVec(taoapp,parameters,ierr)
  CHKERRQ(ierr)
  call TaoAppSetObjectiveRoutine(taoapp,FormFunction,TAO_NULL_OBJECT,ierr)
  CHKERRQ(ierr)  
  call TaoAppSetVariableBounds(taoapp,xl,xu,ierr)
  CHKERRQ(ierr)  
  
  call extractPressure(50,'breakthrough_ref.tec',reference_size,reference_times,reference_pressure)

    

  call TaoSetOptions(taoapp,tao,ierr)
  CHKERRQ(ierr)
  call TaoSolveApplication(taoapp,tao,ierr)
  CHKERRQ(ierr)

  call TaoGetTerminationReason(tao,reason,ierr)
  CHKERRQ(ierr)
  print *,'Termination reason: ',reason
  call VecView(parameters,PETSC_VIEWER_STDOUT_SELF,ierr)
  CHKERRQ(ierr)
  
  call TaoApplicationDestroy(taoapp,ierr)
  CHKERRQ(ierr)
  call TaoDestroy(tao,ierr)
  CHKERRQ(ierr)

  call VecDestroy(xl,ierr)
  CHKERRQ(ierr)
  call VecDestroy(xu,ierr)
  CHKERRQ(ierr)
  call VecDestroy(parameters,ierr)
  CHKERRQ(ierr)
  
  call TaoFinalize (ierr)
  CHKERRQ(ierr)
  call PetscFinalize (ierr)
  CHKERRQ(ierr)

  end program pflotran



!*********************************************************

  subroutine FormFunction(taoapp, X, f, dummy, ierr)
      
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module
  
  implicit none
#include "definitions.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/tao_solver.h"
  TAO_APPLICATION :: taoapp
  Vec :: X
  PetscScalar :: f
  PetscInt :: dummy
  PetscErrorCode :: ierr

  PetscScalar, pointer :: x_val(:)
  integer :: MAXLINELENGTH
  parameter (MAXLINELENGTH=128)
  integer :: MAXTIMESTEPS
  parameter (MAXTIMESTEPS=3000)      
  character*(MAXLINELENGTH) :: buffer


  PetscLogDouble :: timex(4), timex_wall(4)

  PetscInt :: stage(10)
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflotranin

  PetscMPIInt :: myrank, commsize
  
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option

  PetscScalar :: cur_pressure(MAXTIMESTEPS,8)
  integer :: i,j,fd1,fd2,cur_size
  PetscScalar :: cur_times(MAXTIMESTEPS)
  PetscScalar :: reference_pressure(MAXTIMESTEPS,8)
  PetscScalar :: reference_times(MAXTIMESTEPS)  
  integer reference_size
  logical usetransform
  common /appctx/ usetransform,reference_size,reference_pressure,reference_times
  character(len=MAXWORDLENGTH) :: inputfilename
  
  double precision :: calculate_norm

  call LoggingCreate()

  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank, ierr)
  CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)
  CHKERRQ(ierr)

  simulation => SimulationCreate()
  realization => simulation%realization
  option => realization%option

  option%myrank = myrank
  option%commsize = commsize

  

!  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflotranin", &
!                             pflotranin, option_found, ierr)

  option_found = .false.
  if(.not.option_found) pflotranin = "pflotran_tao.in"

  call VecGetArrayF90(X,x_val,ierr)
  CHKERRQ(ierr)
  fd1=25
  fd2=26
  do i=1,6
      if (usetransform) then
         if (x_val(i).le.-25d0 .or. x_val(i) .gt. -10d0) then
           print *,'Parameter ',i,'(',x_val(i),') is out of range'
           f=1.0e30
           call VecRestoreArrayF90(X,x_val,ierr)
           CHKERRQ(ierr)
           return
         endif
      else if (x_val(i).le.0 .or. x_val(i).gt.1.e-10) then
         print *,'Parameter ',i,'(',x_val(i),') is out of range'
         f=1.0e30
         call VecRestoreArrayF90(X,x_val,ierr)
         CHKERRQ(ierr)
         return
      endif
  enddo

  

! delete old breakthrough_0.tec file
  open(unit=fd1,file='breakthrough_0.tec',status='unknown')
  write(fd1,*),''
  close(fd1)

  inputfilename="pflotran_tao.in"
  open(unit=fd1,status='unknown',file=inputfilename)
  open(unit=fd2,status='old',file='pflotranin.start')
  do i=1,35
    read (fd2,50),buffer
    write (fd1,75),buffer
  enddo
  close(fd2)
  write (fd1,75),'MATERIALS'
  write (fd1,75),':name id icap ithm por  tau permx  permy  permz  permpwr'
!'
  if (usetransform) then
    write (fd1,100),1,1,1,1,0.25d0,0.5d0,10.0**x_val(1),10.0**x_val(2),10.0**x_val(3)
    write (fd1,100),2,2,1,1,0.25d0,0.5d0,10.0**x_val(4),10.0**x_val(5),10.0**x_val(6)
    write (*,100),1,1,1,1,0.25d0,0.5d0,10.0**x_val(1),10.0**x_val(2),10.0**x_val(3)
    write (*,100),2,2,1,1,0.25d0,0.5d0,10.0**x_val(4),10.0**x_val(5),10.0**x_val(6)
  else
    write (fd1,100),1,1,1,1,0.25d0,0.5d0,x_val(1),x_val(2),x_val(3)
    write (fd1,100),2,2,1,1,0.25d0,0.5d0,x_val(4),x_val(5),x_val(6)
    write (*,100),1,1,1,1,0.25d0,0.5d0,x_val(1),x_val(2),x_val(3)
    write (*,100),2,2,1,1,0.25d0,0.5d0,x_val(4),x_val(5),x_val(6)
  endif
  write (fd1,75),'END'
  open(unit=fd2,status='old',file='pflotranin.end')
  do i=1,140
    read (fd2,50),buffer
    write (fd1,75),buffer
  enddo
  close(fd1)
  close(fd2)
 50   format(a128)
 75   format(TL1,a)
 100  format(TL1,'soil',I1,' ',I3,' ',I4,' ',I4,' ',F6.2,' ',F6.2,' ', 3(E10.2,' '),'1.d0')
  call VecRestoreArrayF90(X,x_val,ierr)
  CHKERRQ(ierr)

  call PetscGetCPUTime(timex(1), ierr)
  CHKERRQ(ierr)
  call PetscGetTime(timex_wall(1), ierr)
  CHKERRQ(ierr)
  option%start_time = timex_wall(1)

  call OptionCheckCommandLine(option)

  call Init(simulation,inputfilename)

  call StepperRun(simulation%realization,simulation%flow_stepper,simulation%tran_stepper)

! Clean things up.
  print *,'Calling SimulationDestroy'
  call SimulationDestroy(simulation)
  print *,'Done Calling SimulationDestroy'

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  CHKERRQ(ierr)
  call PetscGetTime(timex_wall(2), ierr)
  CHKERRQ(ierr)
  
  if (myrank == 0) then
  
    write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
      (timex(2)-timex(1))/3600.d0

    write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
      (timex_wall(2)-timex_wall(1))/3600.d0

    write(IUNIT2,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
      (timex(2)-timex(1))/3600.d0

    write(IUNIT2,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
      (timex_wall(2)-timex_wall(1))/3600.d0
  endif

  close(IUNIT2)

!  call LoggingDestroy()
  call extractPressure(MAXTIMESTEPS,'breakthrough_0.tec',cur_size,cur_times,cur_pressure)
  f = calculate_norm(cur_size,cur_times,cur_pressure)

  ierr=0
  end subroutine FormFunction      



      subroutine extractPressure(maxn,filename,size,times,pressure)
      implicit none
      integer maxn,size
      character*(128) filename
      PetscScalar :: pressure(3000,8)
      PetscScalar :: times(3000)

      character dummy
      integer fd
      integer i,io
      PetscScalar :: d(8)
      PetscScalar :: p(8)
      
      pressure(1,1)=7.0
      fd=29
      
      open(unit=fd,file=filename,status='old',iostat=io)
      read (fd,*) dummy
      if (dummy .eq. 'T') then
            print *,'TAO IGNORING LINE: ',dummy
      else
        close(fd)
	open(unit=fd,file=filename,status='old',iostat=io)
      endif	
	
      size=0
      do i=1,maxn
         read (fd,*,iostat=io) d(1),p(1),d(2),p(2),d(3),p(3),           &
     &                         d(4),p(4),d(5),p(5),d(6),p(6),           &
     &                         d(7),p(7),d(8),p(8)  
         if (io .eq. 0) then
            pressure(i,1)=p(1)
            pressure(i,2)=p(2)
            pressure(i,3)=p(3)
            pressure(i,4)=p(4)
            pressure(i,5)=p(5)
            pressure(i,6)=p(6)
            pressure(i,7)=p(7)
            pressure(i,8)=p(8)
            times(i)=d(1)
            size=size+1
         else
            exit
         endif
            
      enddo
      close(fd)
      end subroutine extractPressure

      
      

      function calculate_norm(size,times,pressure)
      implicit none
      double precision :: calculate_norm
      integer size,index
      integer :: MAXTIMESTEPS
      parameter (MAXTIMESTEPS=3000)      
      double precision :: interpolated_pressure(MAXTIMESTEPS,8)
      double precision :: pressure(MAXTIMESTEPS,8),times(MAXTIMESTEPS)
      integer :: reference_size
      double precision :: reference_pressure(MAXTIMESTEPS,8)
      double precision :: reference_times(MAXTIMESTEPS)
      logical usetransform
      common /appctx/ usetransform,reference_size,reference_pressure,reference_times


      double precision :: f
      integer i,j

      f=0.0d0
      call interpolate(size,times,pressure,interpolated_pressure)
      do i=1,reference_size
         do j=1,8
            f = f + (interpolated_pressure(i,j)-reference_pressure(i,j))**2
         enddo
      enddo



      calculate_norm=sqrt(f)
      return 
      end function calculate_norm

      subroutine interpolate(n,times,pressure,interpolated_pressure)
      implicit none

      integer :: n
      double precision :: times(3000)
      double precision :: pressure(3000,8)
      integer :: MAXTIMESTEPS
      parameter (MAXTIMESTEPS=3000)      
      double precision :: interpolated_pressure(MAXTIMESTEPS,8)

      integer :: reference_size,i
      double precision :: reference_pressure(MAXTIMESTEPS,8)
      double precision :: reference_times(MAXTIMESTEPS)
      logical usetransform
      common /appctx/ usetransform,reference_size,reference_pressure,reference_times

      integer :: t1_index,t2_index,ref_index
      double precision :: t1,t2,tr
      logical stopme

      call pressuresort(n,times,pressure)
      t1_index = 1
      t2_index = 2
      t1=times(1)
      t2=times(2)
      
      ! For each time in the reference output ...
      do ref_index=1,reference_size
         tr = reference_times(ref_index)
         ! ... find the entries in current output that bracket that time ...
	 stopme = .false.
         do while (stopme .eqv. .false.)
            if (t2 .lt. tr) then
               t1_index=t1_index+1
               t2_index=t2_index+1
               t1 = t2
               t2 = times(t2_index)
	    else 
	      stopme = .true.
            endif

            if (t2_index .gt. n) then
               print *,'Error interpolating times'
               print *,'n=',n
               print *,'t1_index=',t1_index
               print *,'t2_index=',t2_index
               print *,'t1=',t1
               print *,'t2=',t2
               print *,'ref_index=',ref_index
               print *,'tr=',tr
               stop
            endif
         enddo
         ! ... and do the interpolation for each pressure reading
         do i=1,8
            if (t2 .eq. tr) then
               interpolated_pressure(ref_index,i) = pressure(t2_index,i)
            elseif (t1 .eq. tr) then
               interpolated_pressure(ref_index,i) = pressure(t1_index,i)
            else				  
               interpolated_pressure(ref_index,i) = pressure(t1_index,i)&
     &      +(pressure(t2_index,i)-pressure(t1_index,i))/(t2-t1)*(tr-t1)

            endif
         enddo
      enddo
      end subroutine interpolate

!     bubble sort -- should be okay algorithm since mostly in order
      subroutine pressuresort(n,times,pressure)
      implicit none
      integer :: n
      double precision :: times(n), pressure(n,8)
      double precision :: tmp_t, tmp_p(8)
      
      integer :: i,j
      logical :: bubbled
      
      bubbled = .true.
      do while (bubbled)
        bubbled = .false.
        do i=1,n-1
	  if (times(i) .gt. times(i+1)) then
             bubbled = .true.
             tmp_t=times(i)
             do j=1,8
                tmp_p(j) = pressure(i,j)
             enddo

             times(i)=times(i+1)
             do j=1,8
                pressure(i,j)=tmp_p(j)
             enddo

             times(i+1)=tmp_t
             do j=1,8
                pressure(i+1,j)=tmp_p(j)
             enddo
          endif
       enddo
      enddo

      end subroutine pressuresort
