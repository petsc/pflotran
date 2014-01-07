module e4d_mat_inv
! Updated 10/24/2001.
!
!cccccccccccccccccccccccc     Program 4.4     cccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Please Note:                                                         c
!                                                                      c
! (1) This computer program is part of the book, "An Introduction to   c
!     Computational Physics," written by Tao Pang and published and    c
!     copyrighted by Cambridge University Press in 1997.               c
!                                                                      c
! (2) No warranties, express or implied, are made for this program.    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

  public :: MIGS

contains

SUBROUTINE MIGS(A,N,X,INDX)
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.
  !gehadd
  implicit none
  integer :: N
  integer :: INDX(N)
  real :: A(N,N)
  real :: X(N,N)
!
  integer :: I, J, K
  real :: B(N,N)
!
!geh      DIMENSION A(N,N),X(N,N),INDX(N),B(N,N)
!

!geh      DO      20 I = 1, N
      DO I = 1, N
!geh        DO    10 J = 1, N
        DO J = 1, N
          B(I,J) = 0.0
!geh   10   CONTINUE
        ENDDO
!geh   20 CONTINUE
      ENDDO
      
!geh      DO      30 I = 1, N
      DO I = 1, N
          B(I,I) = 1.0
!geh   30 CONTINUE
      ENDDO
!
      CALL ELGS(A,N,INDX)
!
!geh      DO     100 I = 1, N-1
!geh        DO    90 J = I+1, N
!geh          DO  80 K = 1, N
      DO I = 1, N-1
        DO J = I+1, N
          DO K = 1, N
            B(INDX(J),K) = B(INDX(J),K) &
                          -A(INDX(J),I)*B(INDX(I),K)
!geh   80     CONTINUE
!geh   90   CONTINUE
!geh  100 CONTINUE
          ENDDO
        ENDDO
      ENDDO
!
!geh      DO     200 I = 1, N
      DO I = 1, N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
!geh        DO   190 J = N-1, 1, -1
        DO J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
!geh          DO 180 K = J+1, N
          DO K = J+1, N
            X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
!geh  180     CONTINUE
          ENDDO
          X(J,I) =  X(J,I)/A(INDX(J),J)
!geh  190   CONTINUE
        ENDDO
!geh  200 CONTINUE
      ENDDO
!
      RETURN
END SUBROUTINE
!
SUBROUTINE ELGS(A,N,INDX)
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.
  !gehadd
  implicit none
  integer :: N
  integer :: INDX(N)
  real :: A(N,N)

  integer :: I, J, K, ITMP
  real :: C(N)
  real :: C1
  real :: PI
  real :: PI1
  real :: PJ

!
!      DIMENSION A(N,N),INDX(N),C(N)
!
! Initialize the index
!
!geh      DO     50    I = 1, N
      DO I = 1, N
        INDX(I) = I
!geh   50 CONTINUE
      ENDDO
!
! Find the rescaling factors, one from each row
!
!geh        DO     100   I = 1, N
        DO I = 1, N
          C1= 0.0
!geh          DO    90   J = 1, N
          DO J = 1, N
            C1 = AMAX1(C1,ABS(A(I,J)))
!geh   90     CONTINUE
          ENDDO
          C(I) = C1
!geh  100   CONTINUE
        ENDDO
!
! Search the pivoting (largest) element from each column
!
!geh      DO     200   J = 1, N-1
      DO J = 1, N-1
        PI1 = 0.0
!geh        DO   150   I = J, N
        DO I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
!geh  150   CONTINUE
        ENDDO
!
! Interchange the rows via INDX(N) to record pivoting order
!
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
!geh        DO   170   I = J+1, N
        DO I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
!geh          DO 160   K = J+1, N
          DO K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
!geh  160     CONTINUE
          ENDDO
!geh  170   CONTINUE
        ENDDO
!geh  200 CONTINUE
      ENDDO
!
      RETURN
END SUBROUTINE
end module e4d_mat_inv
