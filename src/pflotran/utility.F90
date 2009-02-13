module Utility_module

#include "definitions.h"

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

  integer :: i, j, k, imax
  real*8 :: aamax, sum, dum

  D=1
  do i=1,N
    aamax=0
    do j=1,N
      if (abs(A(i,j)).gt.aamax) aamax=abs(A(i,j))
    enddo
    if (aamax.eq.0) &
      pause 'Singular Matrix.'
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
! Interpolate: Interpolation log K function: temp - temperature [C]
!                             b - fit coefficients determined from fit(...)
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !

double precision function flogk(b,temp)

  PetscReal :: b(5), temp

  flogk = b(1)*log(temp) &
        + b(2)           &
        + b(3)*temp      &
        + b(4)/temp      &
        + b(5)/(temp*temp)
     
end function flogk


! ************************************************************************** !
!
! Interpolate: Least squares fit to log K over database temperature range
! author: P.C. Lichtner
! date: 02/13/09
!
! ************************************************************************** !
#if 0
subroutine fit (n0,nbasis,ntemp,iflgint,alogk0,w,bvec,vec,indx)

  integer :: i,j,n0,ndimmx,ntemp,ntmpmx,iflgint
  
!     include 'impl.h'
!     include 'paramtrs.h'

  integer :: int(ntmpmx),indx(ndimmx)
  PetscReal :: alogk0(ntmpmx),w(ndimmx,ndimmx),vec(5,*),bvec(5)

  iflgint = 0
  do j = 1, nbasis
    bvec(j) = 0.d0
    do i = 1, ntemp
      if (alogk0(i) .ne. 500.) then
        bvec(j) = bvec(j) + alogk0(i)*vec(j,i)
        int(i) = 1
      else if (alogk0(i) .eq. 500.) then
        iflgint = 1
        int(i) = 0
      else if (alogk0(i) .gt. 500.) then
        write(*,*) 'error in fit: log K .gt. 500---stop!'
        stop
      endif
    enddo
  enddo

  do j = 1, nbasis
    do k = j, nbasis
      w(j,k) = 0.d0
      do i = 1, ntemp
        if (int(i) .eq. 1) then
          w(j,k) = w(j,k) + vec(j,i)*vec(k,i)
        endif
      enddo
      if (j .ne. k) w(k,j) = w(j,k)
    enddo
  enddo
  call ludcmp(w,nbasis,n0,indx,dd)
  call lubksb(w,nbasis,n0,indx,bvec)

end subroutine fit
#endif

end module Utility_module
