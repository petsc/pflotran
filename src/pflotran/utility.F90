module Utility_module

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  interface DotProduct
    module procedure DotProduct1
    module procedure DotProduct2
    module procedure DotProduct3
  end interface
  
  interface CrossProduct
    module procedure CrossProduct1
  end interface
  
  interface UtilityReadArray
    module procedure UtilityReadIntArray
    module procedure UtilityReadRealArray
  end interface

  interface DeallocateArray
    ! TODO(geh) replace deallocations with the below
    module procedure DeallocateArray1DInteger
    module procedure DeallocateArray2DInteger
    module procedure DeallocateArray3DInteger
    module procedure DeallocateArray1DReal
    module procedure DeallocateArray2DReal
    module procedure DeallocateArray3DReal
    module procedure DeallocateArray1DLogical
    module procedure DeallocateArray2DLogical
    module procedure DeallocateArray3DLogical
    module procedure DeallocateArray1DString
    module procedure DeallocateArray2DString
  end interface
  
contains

function rnd()

  implicit none

!  common/rndnr/iseed
  integer*8 iseed
  PetscReal :: rnd

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
      
  implicit none

!-----Minimal random number generator of Park and Miller 
!     in the range (0,1)

  PetscReal :: ran2, AM, EPS, RNMX, temp
  PetscInt :: IA, IM, IQ, IR, NTAB, idum, iy, j, k, iv(32), NDIV

  parameter (IA = 16807)
  parameter (IM = 2147483647)
  parameter (AM = 1.0/IM)
  parameter (IQ = 127773)
  parameter (IR = 2836)
  parameter (NTAB = 32)
  parameter (NDIV = 1+(IM-1)/NTAB)
  parameter (EPS  = 1.2e-7)
  parameter (RNMX = 1.0-EPS)

  !dimension iv(NTAB)

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

  PetscInt :: N
  PetscReal, parameter :: tiny=1.0d-20
  PetscReal :: A(N,N),VV(N)
  PetscInt :: INDX(N)
  PetscInt :: D

  PetscInt :: i, j, k, imax
  PetscReal :: aamax, sum, dum
  PetscMPIInt ::  rank
  PetscErrorCode :: ierr

  D=1
  do i=1,N
    aamax=0.d0
    do j=1,N
      if (abs(A(i,j)).gt.aamax) aamax=abs(A(i,j))
    enddo
    if (aamax <= 0.d0) then
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      print *, "ERROR: Singular value encountered in ludcmp() on processor: ", rank, ' aamax = ',aamax,' row = ',i
      do k = 1, N
        print *, "Jacobian: ",k,(j,A(k,j),j=1,N)
      enddo
      call MPI_Abort(MPI_COMM_WORLD,ONE_INTEGER_MPI,ierr)
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
      dum=1.d0/A(j,j)
      do i=j+1,N
        A(i,j)=A(i,j)*dum
      enddo
    endif
  enddo
  return

end subroutine ludcmp

!* Solves the set of N linear equations A.X=D. Here A is input, not as a matrix
!* A but rather as its LU decomposition. INDX is the input as the permutation
!* vector returned by ludcmp. B is input as the right-hand side vector B, and
!* returns with the solution vector X.
subroutine lubksb(A,N,INDX,B)

  implicit none

  PetscInt :: N
  PetscReal :: A(N,N),B(N)
  PetscInt :: INDX(N)

  PetscInt :: i, j, ii, ll
  PetscReal :: sum


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

!* Given an NxN matrix A, with physical dimension NP, this routine replaces it
!* by the LU decomposition of a rowwise permutation of itself.
!* A and N are input. A is output; INDX is output vector which records the
!* row permutation effected by the partial pivoting; D id output as +1 or -1
!* depending on whether the number of row interchanges was odd or even,
!* respectively. This routine is used in combination with lubksb to solve
!* linear equations or invert a matrix.
subroutine ludcmp_chunk(A,N,INDX,D,chunk_size,ithread,num_threads)

  implicit none

  PetscInt :: N
  PetscInt :: chunk_size
  PetscInt :: num_threads
  PetscReal, parameter :: tiny=1.0d-20
  PetscReal :: A(chunk_size,num_threads,N,N),VV(chunk_size,num_threads,N)
  PetscInt :: INDX(chunk_size,num_threads,N)
  PetscInt :: D(chunk_size,num_threads)
  PetscInt :: ithread

  PetscInt :: i, j, k, imax
  PetscReal :: aamax, sum, dum
  PetscMPIInt ::  rank
  PetscErrorCode :: ierr
  
  PetscInt :: ichunk

  do ichunk = 1, chunk_size

  D(ichunk,ithread)=1
  do i=1,N
    aamax=0
    do j=1,N
      if (abs(A(ichunk,ithread,i,j)).gt.aamax) aamax=abs(A(ichunk,ithread,i,j))
    enddo
    if (aamax.eq.0) then
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      print *, "ERROR: Singular value encountered in ludcmp() on processor", rank, ichunk,ithread
      call MPI_Abort(MPI_COMM_WORLD,ONE_INTEGER_MPI,ierr)
      call MPI_Finalize(ierr)
      stop
    endif
    VV(ichunk,ithread,i)=1./aamax
  enddo
  do j=1,N
    do i=1,j-1
      sum=A(ichunk,ithread,i,j)
      do k=1,i-1
        sum=sum-A(ichunk,ithread,i,k)*A(ichunk,ithread,k,j)
      enddo
      A(ichunk,ithread,i,j)=sum
    enddo
    aamax=0
    do i=j,N
      sum=A(ichunk,ithread,i,j)
      do k=1,j-1
        sum=sum-A(ichunk,ithread,i,k)*A(ichunk,ithread,k,j)
      enddo
      A(ichunk,ithread,i,j)=sum
      dum=VV(ichunk,ithread,i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j.ne.imax) then
      do k=1,N
        dum=A(ichunk,ithread,imax,k)
        A(ichunk,ithread,imax,k)=A(ichunk,ithread,j,k)
        A(ichunk,ithread,j,k)=dum
      enddo
      D(ichunk,ithread)=-D(ichunk,ithread)
      VV(ichunk,ithread,imax)=VV(ichunk,ithread,j)
    endif
    INDX(ichunk,ithread,j)=imax
    if (A(ichunk,ithread,j,j).eq.0.) A(ichunk,ithread,j,j)=tiny
    if (j.ne.N) then
      dum=1./A(ichunk,ithread,j,j)
      do i=j+1,N
        A(ichunk,ithread,i,j)=A(ichunk,ithread,i,j)*dum
      enddo
    endif
  enddo
  
  enddo ! chunk loop
  
  return

end subroutine ludcmp_chunk

!* Solves the set of N linear equations A.X=D. Here A is input, not as a matrix
!* A but rather as its LU decomposition. INDX is the input as the permutation
!* vector returned bu ludcmp. B is input as the right-hand side vector B, and
!* returns with the solution vector X.
subroutine lubksb_chunk(A,N,INDX,B,chunk_size,ithread,num_threads)

  implicit none

  PetscInt :: N
  PetscInt :: chunk_size
  PetscInt :: num_threads
  PetscReal :: A(chunk_size,num_threads,N,N),B(chunk_size,num_threads,N)
  PetscInt :: INDX(chunk_size,num_threads,N)
  PetscInt :: ithread

  PetscInt :: i, j, ii, ll
  PetscReal :: sum

  PetscInt :: ichunk

  do ichunk = 1, chunk_size
  
  ii=0
  do i=1,N
    ll=INDX(ichunk,ithread,i)
    sum=B(ichunk,ithread,ll)
    B(ichunk,ithread,ll)=B(ichunk,ithread,i)
    if (ii.ne.0) then
      do j=ii,i-1
        sum=sum-A(ichunk,ithread,i,j)*B(ichunk,ithread,j)
      enddo
    else if (sum.ne.0) then
      ii=i
    endif
    B(ichunk,ithread,i)=sum
  enddo
  do i=N,1,-1
    sum=B(ichunk,ithread,i)
    if (i.lt.N) then
      do j=i+1,N
        sum=sum-A(ichunk,ithread,i,j)*B(ichunk,ithread,j)
      enddo
    endif
    B(ichunk,ithread,i)=sum/A(ichunk,ithread,i,i)
  enddo
  
  enddo ! chunk loop
  
  return

end subroutine lubksb_chunk

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
    y = y_low
  else
    weight = (x-x_low)/(x_high-x_low)
    y = y_low + weight*(y_high-y_low)
  endif

end subroutine Interpolate

! ************************************************************************** !
!
! InterpolateBilinear: Interpolates values between four reference values in 2D
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
function InterpolateBilinear(x,y,x1,x2,y1,y2,z1,z2,z3,z4)

  implicit none

  PetscReal :: x,y,x1,x2,y1,y2,z1,z2,z3,z4
  PetscReal :: InterpolateBilinear
  
  
  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2
  
  
  InterpolateBilinear = (z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+ &
                         z3*(x2-x)*(y-y1)+z4*(x-x1)*(y-y1))/ &
                        ((x2-x1)*(y2-y1))

end function InterpolateBilinear

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
! CrossProduct1: Computes the cross product between two 3d vectors
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
function CrossProduct1(v1,v2)

  implicit none
  
  PetscReal :: v1(3), v2(3)
  
  PetscReal :: CrossProduct1(3)
  
  CrossProduct1(1) = v1(2)*v2(3)-v1(3)*v2(2)
  CrossProduct1(2) = v1(3)*v2(1)-v1(1)*v2(3)
  CrossProduct1(3) = v1(1)*v2(2)-v1(2)*v2(1)

end function CrossProduct1

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
! UtilityReadIntArray: Reads an array of integers from an input file
! author: Glenn Hammond
! date: 11/30/11
!
! ************************************************************************** !
subroutine UtilityReadIntArray(array,array_size,comment,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type), target :: input
  character(len=MAXSTRINGLENGTH) :: comment
  PetscInt :: array_size
  PetscInt, pointer :: array(:)
  
  PetscInt :: i, num_values, count
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word, word2, word3
  character(len=1) :: backslash
  PetscBool :: continuation_flag
  PetscInt :: value
  PetscInt, pointer :: temp_array(:)
  PetscInt :: max_size
  PetscErrorCode :: ierr

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  max_size = 1000
  if (array_size > 0) then
    max_size = array_size
  endif
  allocate(temp_array(max_size))
  temp_array = 0
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','UtilityReadIntArray')
  call StringToLower(word)
  if (StringCompare(word,'file',FOUR_INTEGER)) then
    call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
    input%err_buf = 'filename'
    input%err_buf2 = comment
    call InputErrorMsg(input,option)
    input2 => InputCreate(input%fid + 1,string2,option)
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
        option%io_buffer = 'Within call to UtilityReadIntArray(), ' // &
                           'insufficient values read: ' // trim(string) // &
                           ' of ' // trim(string2) // '.'
      else
        option%io_buffer = 'Within call to UtilityReadIntArray() in ' // &
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
      call InputReadStringErrorMsg(input2,option,comment)
    endif

    continuation_flag = PETSC_FALSE
    if (index(input2%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input2,option,word,PETSC_TRUE)
      if (InputError(input2) .or. StringCompare(word,backslash,ONE_INTEGER)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input2%ierr)
        call InputErrorMsg(input2,option,'# values','UtilityReadIntArray')
        string2 = word3
        call InputReadInt(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value','UtilityReadIntArray')
        do while (count+num_values > max_size)
          ! careful.  reallocateRealArray double max_size every time.
          call reallocateIntArray(temp_array,max_size) 
        enddo
        do i=1, num_values
          count = count + 1
          temp_array(count) = value
        enddo
      else
        string2 = word
        call InputReadInt(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value','UtilityReadIntArray')
        count = count + 1
        if (count > max_size) then
          ! careful.  reallocateRealArray double max_size every time.
          call reallocateIntArray(temp_array,max_size) 
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

end subroutine UtilityReadIntArray

! ************************************************************************** !
!
! UtilityReadRealArray: Reads an array of double precision numbers from the  
!                       input file
! author: Glenn Hammond
! date: 05/21/09
!
! ************************************************************************** !
subroutine UtilityReadRealArray(array,array_size,comment,input,option)

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
  PetscBool :: continuation_flag
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
    input2 => InputCreate(input%fid + 1,string2,option)
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
        option%io_buffer = 'Within call to UtilityReadRealArray(), ' // &
                           'insufficient values read: ' // trim(string) // &
                           ' of ' // trim(string2) // '.'
      else
        option%io_buffer = 'Within call to UtilityReadRealArray() in ' // &
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
      call InputReadStringErrorMsg(input2,option,comment)
    endif

    continuation_flag = PETSC_FALSE
    if (index(input2%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input2,option,word,PETSC_TRUE)
      if (InputError(input2) .or. StringCompare(word,backslash,ONE_INTEGER)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input2%ierr)
        call InputErrorMsg(input2,option,'# values','UtilityReadRealArray')
        string2 = word3
        call InputReadDouble(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value','UtilityReadRealArray')
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

end subroutine UtilityReadRealArray

! ************************************************************************** !
!
! SearchOrderedArray: Locates an integer value in an ordered array and
!                     returned the index
! author: Glenn Hammond
! date: 10/21/09
!
! ************************************************************************** !
function SearchOrderedArray(array,array_length,int_value)

  implicit none

  PetscInt :: array_length 
  PetscInt :: array(array_length)
  PetscInt :: int_value

  PetscInt :: SearchOrderedArray
  PetscInt :: i
  PetscInt :: array_value
  PetscInt :: upper_bound, lower_bound

  SearchOrderedArray = -1

  upper_bound = array_length
  lower_bound = 1

  i = array_length/2
  if (i == 0) i = 1

  do 
    array_value = array(i)
    if (array_value == int_value) then
      SearchOrderedArray = i
      return
    endif
    if (array_value > int_value) then
      upper_bound = i 
    else
      lower_bound = i 
    endif
    i = lower_bound + (upper_bound-lower_bound) / 2
    if (i == lower_bound) then
      if (array(lower_bound) == int_value) SearchOrderedArray = lower_bound
      if (array(upper_bound) == int_value) SearchOrderedArray = upper_bound
      return
    endif
  enddo

end function SearchOrderedArray

! ************************************************************************** !
!
! FileExists: Returns PETSC_TRUE if file exists
! author: Glenn Hammond
! date: 04/27/11
!
! ************************************************************************** !
function FileExists(filename)

  implicit none
  
  PetscBool FileExists

  character(len=*) :: filename
  PetscInt :: fid
  PetscInt :: ios
  
  ios = 0
  fid = 86
  open(unit=fid,file=filename,action="read",status="old",iostat=ios)
  if (ios == 0) then
    close(fid)
    FileExists = PETSC_TRUE
  else
    FileExists = PETSC_FALSE
  endif

end function FileExists

! ************************************************************************** !
!
! Equal: Returns PETSC_TRUE if values are equal
! author: Glenn Hammond
! date: 04/27/11
!
! ************************************************************************** !
function Equal(value1, value2)

  implicit none
  
  PetscBool :: Equal

  PetscReal :: value1, value2

  Equal = PETSC_FALSE
  if (dabs(value1 - value2) <= 1.d-14 * dabs(value1))  Equal = PETSC_TRUE
  
end function Equal

! ************************************************************************** !
!
! BestFloat: Returns the best format for a floating point number
! author: Glenn Hammond
! date: 11/21/11
!
! ************************************************************************** !
function BestFloat(float,upper_bound,lower_bound)

  implicit none
  
  PetscReal :: float
  PetscReal :: upper_bound
  PetscReal :: lower_bound

  character(len=MAXWORDLENGTH) :: BestFloat
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  
100 format(f12.3)
101 format(es12.2)
102 format(es12.4)

  if (dabs(float) <= upper_bound .and. dabs(float) >= lower_bound) then
    write(word,100) float
    word = adjustl(word)
    do i = len_trim(word), 1, -1
      if (word(i:i) == '0') then
        word(i:i) = ' '
      else
        exit
      endif
    enddo
  else if (dabs(float) < lower_bound) then
    write(word,101) float
  else
    write(word,102) float
  endif
  
  BestFloat = adjustl(word)
  
end function BestFloat

! ************************************************************************** !
!
! CubicPolynomialSetup: Sets up a cubic polynomial for smoothing 
!                       discontinuous functions
! author: Glenn Hammond
! date: 03/12/12
!
! ************************************************************************** !
subroutine CubicPolynomialSetup(upper_value,lower_value,coefficients)

  implicit none

  PetscReal :: upper_value
  PetscReal :: lower_value
  PetscReal :: coefficients(4)
  
  PetscReal :: A(4,4)
  PetscInt :: indx(4)
  PetscInt :: d

  A(1,1) = 1.d0
  A(2,1) = 1.d0
  A(3,1) = 0.d0
  A(4,1) = 0.d0
  
  A(1,2) = upper_value
  A(2,2) = lower_value
  A(3,2) = 1.d0
  A(4,2) = 1.d0
  
  A(1,3) = upper_value**2.d0
  A(2,3) = lower_value**2.d0
  A(3,3) = 2.d0*upper_value
  A(4,3) = 2.d0*lower_value
  
  A(1,4) = upper_value**3.d0
  A(2,4) = lower_value**3.d0
  A(3,4) = 3.d0*upper_value**2.d0
  A(4,4) = 3.d0*lower_value**2.d0
  
  ! coefficients(1): value at upper_value
  ! coefficients(2): value at lower_value
  ! coefficients(3): derivative at upper_value
  ! coefficients(4): derivative at lower_value
  
  call ludcmp(A,FOUR_INTEGER,indx,d)
  call lubksb(A,FOUR_INTEGER,indx,coefficients)

end subroutine CubicPolynomialSetup

! ************************************************************************** !
!
! CubicPolynomialEvaluate: Evaluates value in cubic polynomial
! author: Glenn Hammond
! date: 03/12/12
!
! ************************************************************************** !
subroutine CubicPolynomialEvaluate(coefficients,x,f,df_dx)

  implicit none

  PetscReal :: coefficients(4)
  PetscReal :: x
  PetscReal :: f
  PetscReal :: df_dx

  PetscReal :: x_squared
  
  x_squared = x*x
  
  f = coefficients(1) + &
      coefficients(2)*x + &
      coefficients(3)*x_squared + &
      coefficients(4)*x_squared*x
  
  df_dx = coefficients(2) + &
          coefficients(3)*2.d0*x + &
          coefficients(4)*3.d0*x_squared
  
end subroutine CubicPolynomialEvaluate

! ************************************************************************** !
!
! DeallocateArray1DInteger: Deallocates a 1D integer array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray1DInteger(array)

  implicit none
  
  PetscInt, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DInteger

! ************************************************************************** !
!
! DeallocateArray2DInteger: Deallocates a 2D integer array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray2DInteger(array)

  implicit none
  
  PetscInt, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DInteger

! ************************************************************************** !
!
! DeallocateArray3DInteger: Deallocates a 3D integer array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray3DInteger(array)

  implicit none
  
  PetscInt, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DInteger

! ************************************************************************** !
!
! DeallocateArray1DReal: Deallocates a 1D real array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray1DReal(array)

  implicit none
  
  PetscReal, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DReal

! ************************************************************************** !
!
! DeallocateArray2DReal: Deallocates a 2D real array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray2DReal(array)

  implicit none
  
  PetscReal, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DReal

! ************************************************************************** !
!
! DeallocateArray3DReal: Deallocates a 3D real array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray3DReal(array)

  implicit none
  
  PetscReal, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DReal

! ************************************************************************** !
!
! DeallocateArray1DLogical: Deallocates a 1D logical array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray1DLogical(array)

  implicit none
  
  PetscBool, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DLogical

! ************************************************************************** !
!
! DeallocateArray2DLogical: Deallocates a 2D logical array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray2DLogical(array)

  implicit none
  
  PetscBool, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DLogical

! ************************************************************************** !
!
! DeallocateArray3DLogical: Deallocates a 3D logical array
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray3DLogical(array)

  implicit none
  
  PetscBool, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DLogical

! ************************************************************************** !
!
! DeallocateArray1DString: Deallocates a 1D array of character strings
! author: Glenn Hammond
! date: 03/13/12
!
! ************************************************************************** !
subroutine DeallocateArray1DString(array)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DString

! ************************************************************************** !
!
! DeallocateArray2DString: Deallocates a 2D array of character strings
! author: Glenn Hammond
! date: 10/30/12
!
! ************************************************************************** !
subroutine DeallocateArray2DString(array)

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DString

end module Utility_module
