module Gauss_module

use Unstructured_Cell_module

  implicit none

#include "definitions.h"

  private
  
  PetscInt, parameter, public :: LINE_TYPE          = 7
  
  type, public :: gauss_type
    PetscInt :: dim                 ! dimension
    PetscInt :: EleType             ! Element type
    PetscInt :: NGPTS               ! Number of gauss points
    PetscReal, pointer :: r(:,:)    ! location of points
    PetscReal, pointer :: w(:)      ! weights
  end type gauss_type
    
  public :: GaussCalculatePoints, GaussDestroy
  
  contains

! ************************************************************************** !
!
! GaussCalculatePoints: Calculates Gauss points 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine GaussCalculatePoints(gauss)

  type(gauss_type) :: gauss
  
  select case(gauss%dim)
    case(ONE_DIM_GRID)
      call Gauss1D(gauss%EleType,gauss%NGPTS,gauss%r,gauss%w)
    case(TWO_DIM_GRID)
      call Gauss2D(gauss%EleType,gauss%NGPTS,gauss%r,gauss%w)
    case(THREE_DIM_GRID)
      call Gauss3D(gauss%EleType,gauss%NGPTS,gauss%r,gauss%w)
    case default
      print *, 'Error: Invalid dimension for Gauss point calculation'
      stop
    end select  

end subroutine GaussCalculatePoints  


! ************************************************************************** !
!
! Gauss1D: Calculates Gauss points for 1D elements 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine Gauss1D(EleType,NGPTS,r,w)

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  allocate(r(NGPTS,1))
  allocate(w(NGPTS))
  
  if (EleType /= LINE_TYPE) then
    print *, 'Error: in Element type. Only L2 ' // &
             '(line type) can be used for 1D Gauss quadrature.'
  endif
    
  select case(NGPTS)
    case(1)

    !------------------------------
    ! No of Gauss Points = 1
    !------------------------------
        
      r(1,1) = 0.0
      !
      w(1) = 2.0

    !------------------------------
    ! No of Gauss Points = 2
    !------------------------------

    case(2)
        
      r(1,1) = -1.0/sqrt(3.0)
      r(2,1) = -r(1,1) 
      !
      w(1) = 1.0
      w(2) = 1.0

    !-------------------------------
    ! No of Gauss Points = 3
    !-------------------------------

    case(3)
        
      r(1,1) = -sqrt(0.6)
      r(2,1) = 0.0
      r(3,1) = -r(1,1)
      !
      w(1) = 5.0/9.0
      w(2) = 8.0/9.0
      w(3) = w(1)

    !--------------------------------
    ! No of Gauss Points = 4
    !--------------------------------

    case(4)
        
      r(1,1) = -0.861136311594053
      r(2,1) = -0.339981043584856
      r(3,1) =  0.339981043584856
      r(4,1) =  0.861136311594053
      ! 
      w(1) = 0.347854845137454
      w(2) = 0.652145154862546
      w(3) = 0.652145154862546
      w(4) = 0.347854845137454

    !----------------------------------
    ! No of Gauss Points = 5
    !----------------------------------

    case(5)
        
      r(1,1) = -0.906179845938664
      r(2,1) = -0.538469310105683
      r(3,1) =  0.000000000000000
      r(4,1) =  0.538469310105683
      r(5,1) =  0.906179845938664
      !
      w(1) =  0.236926885056189
      w(2) =  0.478628670499366 
      w(3) =  0.568888888888889
      w(4) =  0.478628670499366
      w(5) =  0.236926885056189

    !----------------------------------
    ! No of Gauss Points = 6
    !----------------------------------
   
    case(6)
        
      r(1,1) = -0.932469514203152
      r(2,1) = -0.661209386466265
      r(3,1) = -0.238619186083197
      r(4,1) =  0.238619186083197
      r(5,1) =  0.661209386466265
      r(6,1) =  0.932469514203152
      !
      w(1) =  0.171324492379170
      w(2) =  0.360761573048139
      w(3) =  0.467913934572691
      w(4) =  0.467913934572691
      w(5) =  0.360761573048139
      w(6) =  0.171324492379170

    !------------------------------------
    ! No of Gauss Points = 7
    !------------------------------------

    case(7)
        
      r(1,1) = -0.949107912342759
      r(2,1) = -0.741531185599394
      r(3,1) = -0.405845151377397
      r(4,1) =  0.000000000000000
      r(5,1) =  0.405845151377397
      r(6,1) =  0.741531185599394
      r(7,1) =  0.949107912342759
      !
      w(1) =  0.129484966168870
      w(2) =  0.279705391489277
      w(3) =  0.381830050505119
      w(4) =  0.417959183673469
      w(5) =  0.381830050505119
      w(6) =  0.279705391489277
      w(7) =  0.129484966168870

    !------------------------------------
    ! No of Gauss Points = 8
    !------------------------------------
    
    case(8)
        
      r(1,1) = -0.960289856497536
      r(2,1) = -0.796666477413627
      r(3,1) = -0.525532409916329
      r(4,1) = -0.183434642495650
      r(5,1) =  0.183434642495650
      r(6,1) =  0.525532409916329
      r(7,1) =  0.796666477413627
      r(8,1) =  0.960289856497536
      !
      w(1) =  0.101228536290376
      w(2) =  0.222381034453374
      w(3) =  0.313706645877887
      w(4) =  0.362683783378362
      w(5) =  0.362683783378362
      w(6) =  0.313706645877887
      w(7) =  0.222381034453374
      w(8) =  0.101228536290376

    !------------------------------------
    ! No of Gauss Points = 9
    !------------------------------------
 
    case(9)
        
      r(1,1) = -0.968160239507626
      r(2,1) = -0.836031170326636
      r(3,1) = -0.613371432700590
      r(4,1) = -0.324253423403809
      r(5,1) =  0.000000000000000
      r(6,1) =  0.324253423403809
      r(7,1) =  0.613371432700590
      r(8,1) =  0.836031107326636
      r(9,1) =  0.968160239507626

      w(1) =  0.081274388361574
      w(2) =  0.180648160694857
      w(3) =  0.260610696402935
      w(4) =  0.312347077040003
      w(5) =  0.330239355001260
      w(6) =  0.312347077040003
      w(7) =  0.260610696402935
      w(8) =  0.180648160694857
      w(9) =  0.081274388361574

   case default
     print *, 'Error in NGPTS for 1D Gauss quadrature'
     stop
   end select

end subroutine Gauss1D  

! ************************************************************************** !
!
! Gauss2D: Calculates Gauss points for 2D elements 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine Gauss2D(EleType,NGPTS,r,w)

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
    
  select case(EleType)
    case(QUAD_TYPE)
      call GaussSquare(NGPTS,r,w)
!    case(T3)
!      call GaussTriangle(NGPTS,r,w)
    case default
      print *, 'Error: Only T3 and Q4 elements available for 2D.'
      stop
  end select

end subroutine Gauss2D

! ************************************************************************** !
!
! GaussSquare: Calculates Gauss points for Q4 element 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine GaussSquare(NGPTS,r,w)

  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter,i,j 
  
  allocate(r(NGPTS,2))
  allocate(w(NGPTS))
  allocate(l(NGPTS,1))
  allocate(m(NGPTS))

  ! 1D Gauss points are stored in l vector and weights are stored in m vector

  call Gauss1D(LINE_TYPE,NGPTS,l,m)
  
  ! Generate the Q4 Gauss points and weights using for loops
  counter = 1
  do i = 1, NGPTS
    do j = 1, NGPTS
      r(counter,1) = l(i,1)
      r(counter,2) = l(j,1)
      w(counter) = m(i)*m(j)
      counter = counter + 1
    enddo
  enddo

  deallocate(l)
  deallocate(m)

end subroutine GaussSquare

! ************************************************************************** !
!
! GaussTriangle: Calculates Gauss points for T3 elements 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
!subroutine GaussTriangle(NGPTS,r,w)
!
!  PetscInt :: NGPTS
!  
!  
!end subroutine GaussTriangle

! ************************************************************************** !
!
! Gauss3D: Calculates Gauss points for 3D element 
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine Gauss3D(EleType,NGPTS,r,w)

  PetscInt :: EleType
  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  
  select case(EleType)
    case(HEX_TYPE)
      call GaussBrick(NGPTS,r,w)
!    case(TET4)
!      call GaussTetrahedra(NGPTS,r,w)
    case default
      print *, 'Error: Only B8 and TET4 elements available for 3D.'
      stop
  end select
  
end subroutine Gauss3D

! ************************************************************************** !
!
! GaussBrick: Calculates Gauss points for B8 element
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !  
subroutine GaussBrick(NGPTS,r,w)

  PetscInt :: NGPTS
  PetscReal, pointer :: r(:,:)
  PetscReal, pointer :: w(:)
  PetscReal, pointer :: l(:,:)
  PetscReal, pointer :: m(:)
  PetscInt :: counter, i, j, k
  
  allocate(r(NGPTS,3))
  allocate(w(NGPTS))
  allocate(l(NGPTS,1))
  allocate(m(NGPTS))
  
  call Gauss1D(LINE_TYPE,NGPTS,l,m)
  
  ! Generate the B8 Gauss points and weights using for loops

  counter = 1
  do i = 1, NGPTS
    do j = 1, NGPTS
      do k = 1, NGPTS
        r(counter,1) = l(i,1)
        r(counter,2) = l(j,1)
        r(counter,3) = l(k,1)
        w(counter) = m(i)*m(j)*m(k)
        counter = counter + 1
      enddo
    enddo
  enddo
  
end subroutine GaussBrick


! ************************************************************************** !
!
! GaussDestroy: Deallocate gauss type
! author: Satish Karra, LANL
! date: 5/17/2013
!
! ************************************************************************** !
subroutine GaussDestroy(gauss)

  type(gauss_type) :: gauss
  
  deallocate(gauss%r)
  nullify(gauss%r)
  deallocate(gauss%w)
  nullify(gauss%w)

end subroutine GaussDestroy

#if 0
subroutine ConvertMatrixToVector(A,vecA)

  PetscReal :: A(:,:)
  PetscReal, pointer :: vecA(:)
  PetscInt :: m, n, i, j
  
  m = size(A,1)
  n = size(A,2)
  
  allocate(vecA(m*n))
  
  do i = 1, m
    do j = 1, n
      vecA(j+(i-1)*n) = A(i,j)
    enddo
  enddo

end subroutine ConvertMatrixToVector
#endif
     
end module Gauss_module

