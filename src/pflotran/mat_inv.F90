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
SUBROUTINE MIGS(A,N,X,INDX)
!
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.

!
      DIMENSION A(N,N),X(N,N),INDX(N),B(N,N)
!

      DO      20 I = 1, N
        DO    10 J = 1, N
          B(I,J) = 0.0
   10   CONTINUE
   20 CONTINUE
      DO      30 I = 1, N
          B(I,I) = 1.0
   30 CONTINUE
!
      CALL ELGS(A,N,INDX)
!
      DO     100 I = 1, N-1
        DO    90 J = I+1, N
          DO  80 K = 1, N
            B(INDX(J),K) = B(INDX(J),K) &
                          -A(INDX(J),I)*B(INDX(I),K)
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
!
      DO     200 I = 1, N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
        DO   190 J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
          DO 180 K = J+1, N
            X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
  180     CONTINUE
          X(J,I) =  X(J,I)/A(INDX(J),J)
  190   CONTINUE
  200 CONTINUE
!
      RETURN
END SUBROUTINE
!
SUBROUTINE ELGS(A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.

!
      DIMENSION A(N,N),INDX(N),C(N)
!
! Initialize the index
!
      DO     50    I = 1, N
        INDX(I) = I
   50 CONTINUE
!
! Find the rescaling factors, one from each row
!
        DO     100   I = 1, N
          C1= 0.0
          DO    90   J = 1, N
            C1 = AMAX1(C1,ABS(A(I,J)))
   90     CONTINUE
          C(I) = C1
  100   CONTINUE
!
! Search the pivoting (largest) element from each column
!
      DO     200   J = 1, N-1
        PI1 = 0.0
        DO   150   I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
  150   CONTINUE
!
! Interchange the rows via INDX(N) to record pivoting order
!
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO   170   I = J+1, N
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
          DO 160   K = J+1, N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
  160     CONTINUE
  170   CONTINUE
  200 CONTINUE
!
      RETURN
END SUBROUTINE
end module e4d_mat_inv
