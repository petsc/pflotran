module Utility_module

#include "definitions.h"

  interface DotProduct
    module procedure DotProduct2
    module procedure DotProduct3
  end interface
  
contains

function rnd()

!  common/rndnr/iseed
  integer*8 iseed

  iseed = iseed*125
  iseed = iseed - (iseed/2796203) * 2796203
  rnd   = iseed/2796203.0
  return
end function rnd

function ran1(idum)
      
  implicit none
      
  save

!-----returns a random number in the range (0,1). Set idum to neg.
!     value to initialize

  PetscReal :: ran1
  PetscReal :: r(97),rm1,rm2
  PetscInt :: idum,iff,ix1,ix2,ix3,j,m1,ia1,ic1,m2,ia2,ic2,m3,ia3,ic3

  parameter (M1  = 259200)
  parameter (IA1 = 7141)
  parameter (IC1 = 54773)
  parameter (RM1 = 1.0/M1)
  parameter (M2  = 134456)
  parameter (IA2 = 8121)
  parameter (IC2 = 28411)
  parameter (RM2 = 1.0/M2)
  parameter (M3  = 243000)
  parameter (IA3 = 4561)
  parameter (IC3 = 51349)

  data iff/0/

  if (idum.lt.0 .or. iff.eq.0) then
    iff=1
    ix1=mod(IC1-idum,M1)
    ix1=mod(IA1*ix1+IC1,M1)
    ix2=mod(ix1,M2)
    ix1=mod(IA1*ix1+IC1,M1)
    ix3=mod(ix1,M3)
    do j=1,97
      ix1=mod(IA1*ix1+IC1,M1)
      ix2=mod(IA2*ix2+IC2,M2)
      r(j)=(float(ix1)+float(ix2)*RM2)*RM1
    enddo
    idum=1
  endif
  ix1=mod(IA1*ix1+IC1,M1)
  ix2=mod(IA2*ix2+IC2,M2)
  ix3=mod(IA3*ix3+IC3,M3)
  j=1+(97*ix3)/M3
!  if (j.gt.97 .or. j.lt.1) pause

  ran1=r(j)
  r(j)=(float(ix1)+float(ix2)*RM2)*RM1
      
  return
end function ran1

function ran2(idum)
      
!  implicit none

!-----Minimal random number generator of Park and Miller 
!     in the range (0,1)

  parameter (IA = 16807)
  parameter (IM = 2147483647)
  parameter (AM = 1.0/IM)
  parameter (IQ = 127773)
  parameter (IR = 2836)
  parameter (NTAB = 32)
  parameter (NDIV = 1+(IM-1)/NTAB)
  parameter (EPS  = 1.2e-7)
  parameter (RNMX = 1.0-EPS)

  dimension iv(NTAB)

  iy = 0
  if (idum.le.0 .or. iy.eq.0) then
    if (-idum .lt. 1) then
      idum = 1
    else 
      idum = -idum
    endif
    do j = NTAB+7,0,-1
      k = idum/IQ
      idum = IA*(idum-k*IQ)-IR*k
      if (idum .lt. 0) idum = idum+IM
      if (j .lt. NTAB) iv(j) = idum
    enddo
    iy = iv(1)
  endif
  k = idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum .lt. 0) idum = idum+IM
  j= iy/NDIV
  iy = iv(j)
  iv(j) = idum
  temp = AM*iy
  if (temp .gt. RNMX) then
    ran2 = RNMX
  else 
    ran2 = temp
  endif

  return
end function ran2
      
      
subroutine Natural2LocalIndex(ir, nl, llist, llength)
  implicit none
  PetscInt :: nl, ir,na, l_search, itt, llength
  PetscInt :: llist(*)
  
  PetscInt ::  nori0, nori1, nori
  
  
  nl=-1
  l_search = llength
  
  na = ir!-1 
  itt=0
  nori0 =1
  nori1 = llength
  if(na>=llist(1) .and. na <= llist(llength))then
  do while(l_search > 1 .and.itt<=50)
  
    itt=itt+1
    if(na == llist(nori0))then
      nl = nori0
      exit
    elseif(na == llist(nori1))then
       nl = nori1
      exit
    endif   
     
     ! nori = int((real(nori0 + nori1))/ 2.) + mod ( nori0 + nori1,2 )
    nori =  int(floor(real(nori0+nori1)/2D0 + .75D0))
    if( na > llist(nori)) then
      nori0 = nori
    elseif(na < llist(nori))then
      nori1 = nori
    else
      if(na == llist(nori))then
        nl = nori
        exit
      else
        print *, 'wrong index', na, nori, llist(nori); stop
      endif  
    endif
    l_search = nori1-nori0
    if (itt>=40)then
      print *, na, nori0,nori1,nori, llist(nori0), llist(nori1)
      if (itt>=50) stop
    endif
  enddo  
 endif         
          
end subroutine Natural2LocalIndex

! ************************************************************************** !
!
! reallocateIntArray: Reallocates an integer array to a larger size and copies
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine reallocateIntArray(array,size)

  implicit none

  PetscInt, pointer :: array(:)
  PetscInt :: size
  
  PetscInt, allocatable :: array2(:)
  
  allocate(array2(size))
  array2(1:size) = array(1:size)
  deallocate(array)
  allocate(array(2*size))
  array = 0
  array(1:size) = array2(1:size)
  size = 2*size
  deallocate(array2)

end subroutine reallocateIntArray

! ************************************************************************** !
!
! reallocateRealArray: Reallocates a real array to a larger size and copies
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine reallocateRealArray(array,size)

  implicit none

  PetscReal, pointer :: array(:)
  PetscInt :: size
  
  PetscReal, allocatable :: array2(:)
  
  allocate(array2(size))
  array2(1:size) = array(1:size)
  deallocate(array)
  allocate(array(2*size))
  array = 0.d0
  array(1:size) = array2(1:size)
  size = 2*size
  deallocate(array2)

end subroutine reallocateRealArray

!* Given an NxN matrix A, with physical dimension NP, this routine replaces it
!* by the LU decomposition of a rowwise permutation of itself.
!* A and N are input. A is output; INDX is output vector which records the
!* row permutation effected by the partial pivoting; D id output as +1 or -1
!* depending on whether the number of row interchanges was odd or even,
!* respectively. This routine is used in combination with lubksb to solve
!* linear equations or invert a matrix.
subroutine ludcmp(A,N,INDX,D)

  implicit none

  integer :: N
  real*8, parameter :: tiny=1.0d-20
  real*8 :: A(N,N),VV(N)
  integer :: INDX(N)
  integer :: D

  integer :: i, j, k, imax, ierr, rank
  real*8 :: aamax, sum, dum

  D=1
  do i=1,N
    aamax=0
    do j=1,N
      if (abs(A(i,j)).gt.aamax) aamax=abs(A(i,j))
    enddo
    if (aamax.eq.0) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      print *, "ERROR: Singular value encountered in ludcmp() on process", rank
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      call MPI_Finalize(ierr)
      stop
    endif
    VV(i)=1./aamax
  enddo
  do j=1,N
    do i=1,j-1
      sum=A(i,j)
      do k=1,i-1
        sum=sum-A(i,k)*A(k,j)
      enddo
      A(i,j)=sum
    enddo
    aamax=0
    do i=j,N
      sum=A(i,j)
      do k=1,j-1
        sum=sum-A(i,k)*A(k,j)
      enddo
      A(i,j)=sum
      dum=VV(i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j.ne.imax) then
      do k=1,N
        dum=A(imax,k)
        A(imax,k)=A(j,k)
        A(j,k)=dum
      enddo
      D=-D
      VV(imax)=VV(j)
    endif
    INDX(j)=imax
    if (A(j,j).eq.0.) A(j,j)=tiny
    if (j.ne.N) then
      dum=1./A(j,j)
      do i=j+1,N
        A(i,j)=A(i,j)*dum
      enddo
    endif
  enddo
  return

end subroutine ludcmp

!* Solves the set of N linear equations A.X=D. Here A is input, not as a matrix
!* A but rather as its LU decomposition. INDX is the input as the permutation
!* vector returned bu ludcmp. B is input as the right-hand side vector B, and
!* returns with the solution vector X.
subroutine lubksb(A,N,INDX,B)

  implicit none

  integer :: N
  real*8 :: A(N,N),B(N)
  integer :: INDX(N)

  integer :: i, j, ii, ll
  real*8 :: sum


  ii=0
  do i=1,N
    ll=INDX(i)
    sum=B(ll)
    B(ll)=B(i)
    if (ii.ne.0) then
      do j=ii,i-1
        sum=sum-A(i,j)*B(j)
      enddo
    else if (sum.ne.0) then
      ii=i
    endif
    B(i)=sum
  enddo
  do i=N,1,-1
    sum=B(i)
    if (i.lt.N) then
      do j=i+1,N
        sum=sum-A(i,j)*B(j)
      enddo
    endif
    B(i)=sum/A(i,i)
  enddo
  return

end subroutine lubksb

! ************************************************************************** !
!
! Interpolate: Interpolates values between two reference values
! author: Glenn Hammond
! date: 02/09/09
!
! ************************************************************************** !
subroutine Interpolate(x_high,x_low,x,y_high,y_low,y)

  implicit none

  PetscReal :: x_high, x_low, x
  PetscReal :: y_high, y_low, y
  
  PetscReal :: weight
  
  if (dabs(x_high-x_low) < 1.d-10) then
    y = x
  else
    weight = (x-x_low)/(x_high-x_low)
    y = y_low + weight*(y_high-y_low)
  endif

end subroutine Interpolate

! ************************************************************************** !
!
! DotProduct1: Computes the dot product between two 3d vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function DotProduct1(v1,v2)

  implicit none
  
  PetscReal :: v1(3), v2(3)
  
  PetscReal :: DotProduct1
  
  DotProduct1 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

end function DotProduct1

! ************************************************************************** !
!
! DotProduct2: Computes the dot product between two 3d vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function DotProduct2(v1,v2x,v2y,v2z)

  implicit none
  
  PetscReal :: v1(3), v2x, v2y, v2z
  
  PetscReal :: DotProduct2
  
  DotProduct2 = v1(1)*v2x+v1(2)*v2y+v1(3)*v2z

end function DotProduct2

! ************************************************************************** !
!
! DotProduct3: Computes the dot product between components of two 3d 
!                    vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function DotProduct3(v1x,v1y,v1z,v2x,v2y,v2z)

  implicit none
  
  PetscReal :: v1x, v1y, v1z, v2x, v2y, v2z
  
  PetscReal :: DotProduct3
  
  DotProduct3 = v1x*v2x+v1y*v2y+v1z*v2z

end function DotProduct3

! ************************************************************************** !
!
! Erf: Computes an approximate to erf(x)
! author: Glenn Hammond
! date: 05/20/09
! from: http://jin.ece.uiuc.edu/routines/merror.for
!
! ************************************************************************** !
function Erf(x)

  implicit none
  
  PetscReal :: x
  
  PetscReal :: Erf

  PetscReal, parameter :: EPS = 1.d-15
  PetscReal, parameter :: PI=3.141592653589793d0
  PetscReal :: x2, er, r, co
  PetscInt :: k
  
  x2=x*x
  if (dabs(x) < 3.5d0) then
    er=1.d0
    r=1.d0
    do k = 1, 50
      r=r*x2/(dble(k)+0.5d0)
      er=er+r
      if (dabs(r) < dabs(er)*EPS) exit
    enddo
    co=2.d0/sqrt(PI)*x*exp(-x2)
    Erf=co*er
  else
    er=1.d0
    r=1.d0
    do k = 1, 12
      r=-r*(dble(k)-0.5d0)/x2
      er=er+r
    enddo
    co=exp(-x2)/(dabs(x)*sqrt(PI))
    Erf=1.d0-co*er
    if (x < 0.d0) Erf=-Erf
  endif

end function Erf

! ************************************************************************** !
!
! Erf: Computes an approximate to erf(x)
! author: Glenn Hammond
! date: 05/20/09
!
!/* adapted from 
!#
!# Lower tail quantile for standard normal distribution function.
!#
!# This function returns an approximation of the inverse cumulative
!# standard normal distribution function.  I.e., given P, it returns
!# an approximation to the X satisfying P = Pr{Z <= X} where Z is a
!# random variable from the standard normal distribution.
!#
!# The algorithm uses a minimax approximation by rational functions
!# and the result has a relative error whose absolute value is less
!# than 1.15e-9.
!#
!# Author:      Peter J. Acklam
!# Time-stamp:  2000-07-19 18:26:14
!# E-mail:      pjacklam@online.no
!# WWW URL:     http://home.online.no/~pjacklam
!*/
!
! ************************************************************************** !
function InverseErf(p)

  implicit none
  
  PetscReal :: p
  
  PetscReal :: InverseErf
  
 ! Coefficients in rational approximations.
  PetscReal, parameter :: A(6) = (/-3.969683028665376d+1,2.209460984245205d+2, &
                                  -2.759285104469687d+2,1.383577518672690d+2, &
                                  -3.066479806614716d+1,2.506628277459239d+0/)
  PetscReal, parameter :: B(5) = (/-5.447609879822406d+1,1.615858368580409d+2, &
                                  -1.556989798598866d+2,6.680131188771972d+1, &
                                  -1.328068155288572d+1/)
  PetscReal, parameter :: C(6) = (/-7.784894002430293d-3,-3.223964580411365d-1, &
                                  -2.400758277161838d+0,-2.549732539343734d+0, &
                                  4.374664141464968d+0,2.938163982698783d+0/)
  PetscReal, parameter :: D(4) = (/7.784695709041462d-03,  3.224671290700398d-01, &
                                  2.445134137142996d+00,  3.754408661907416d+0/)

  ! Define break-points.
  PetscReal, parameter :: PLOW  = 0.02425d0;
  PetscReal, parameter :: PHIGH = 0.97575d0 ! 1 - PLOW;
  PetscReal :: q, r

  ! Rational approximation for lower region:
  if (p < PLOW) then
    q = sqrt(-2.d0*log(p))
    InverseErf = (((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6)) / &
                  ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.d0)
  ! Rational approximation for upper region:
  elseif (PHIGH < p) then
    q = sqrt(-2.d0*log(1.d0-p))
    InverseErf = -(((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6)) / &
                   ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.d0)
  ! Rational approximation for central region:
  else
    q = p - 0.5d0;
    r = q*q;
    InverseErf = (((((A(1)*r+A(2))*r+A(3))*r+A(4))*r+A(5))*r+A(6))*q / &
                 (((((B(1)*r+B(2))*r+B(3))*r+B(4))*r+B(5))*r+1.d0)
  endif

end function InverseErf

! ************************************************************************** !
!
! UtilityReadArray: Reads an array of double precision numbers from the  
!                   input file
! author: Glenn Hammond
! date: 05/21/09
!
! ************************************************************************** !
subroutine UtilityReadArray(array,array_size,comment,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type), target :: input
  character(len=MAXSTRINGLENGTH) :: comment
  PetscInt :: array_size
  PetscReal, pointer :: array(:)
  
  PetscInt :: i, num_values, count
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word, word2, word3
  character(len=1) :: backslash
  PetscTruth :: continuation_flag
  PetscReal :: value
  PetscReal, pointer :: temp_array(:)
  PetscInt :: max_size
  PetscErrorCode :: ierr

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  max_size = 1000
  if (array_size > 0) then
    max_size = array_size
  endif
  allocate(temp_array(max_size))
  temp_array = 0.d0
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToLower(word)
  if (StringCompare(word,'file',FOUR_INTEGER)) then
    call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
    input%err_buf = 'filename'
    input%err_buf2 = comment
    call InputErrorMsg(input,option)
    input2 => InputCreate(input%fid + 1,string2)
  else
    input2 => input
    input%buf = string2
  endif
  
  if (len_trim(input2%buf) > 1) then
    continuation_flag = PETSC_FALSE
  else
    continuation_flag = PETSC_TRUE
  endif
  
  count = 0
  do
    
    if (count >= array_size .and. array_size > 0) exit
    
    if (.not.continuation_flag .and. count /= 1 .and. array_size > 0) then
      write(string,*) count
      write(string2,*) array_size
      if (len_trim(comment) < 1) then
        option%io_buffer = 'Within call to UtilityReadArray(), ' // &
                           'insufficient values read: ' // trim(string) // &
                           ' of ' // trim(string2) // '.'
      else
        option%io_buffer = 'Within call to UtilityReadArray() in ' // &
                           trim(comment) // &
                           'insufficient values read: ' // trim(string) // &
                           ' of ' // trim(string2) // '.'
      endif
      call printErrMsg(option)
    else if (count == 1) then
      temp_array = temp_array(count)
      exit
    else if (.not.continuation_flag .and. array_size <= 0 .and. count /= 0) then
      exit
    endif
    
    if (continuation_flag) then
      call InputReadFlotranString(input2,option)
      call InputReadStringErrorMsg(input2,option,'DXYZ')
    endif

    continuation_flag = PETSC_FALSE
    if (index(input2%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input2,option,word,PETSC_TRUE)
      if (InputError(input2) .or. StringCompare(word,backslash,1)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input2%ierr)
        call InputErrorMsg(input2,option,'# values','UtilityReadArray')
        string2 = word3
        call InputReadDouble(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value','UtilityReadArray')
        do while (count+num_values > max_size)
          ! careful.  reallocateRealArray double max_size every time.
          call reallocateRealArray(temp_array,max_size) 
        enddo
        do i=1, num_values
          count = count + 1
          temp_array(count) = value
        enddo
      else
        string2 = word
        call InputReadDouble(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value','UtilityReadArray')
        count = count + 1
        if (count > max_size) then
          ! careful.  reallocateRealArray double max_size every time.
          call reallocateRealArray(temp_array,max_size) 
        endif
        temp_array(count) = value
      endif
    enddo
  enddo
  
  if (array_size > 0 .and. count > array_size) then
    count = array_size
  endif
  
  if (.not.associated(input2,input)) call InputDestroy(input2)
  nullify(input2)
  
  if (associated(array)) deallocate(array)
  allocate(array(count))
  array(1:count) = temp_array(1:count)
  deallocate(temp_array)
  nullify(temp_array)

end subroutine UtilityReadArray

end module Utility_module
