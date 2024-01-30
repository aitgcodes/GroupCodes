C=======================================================================
C==            COMPUTES THE INVERSE OF A MATRIX 3*3                   ==
C=======================================================================

      SUBROUTINE INV3(A,C,DEN)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      REAL*8 A(3,3),B(3,3),C(3,3)

      DEN = 0.D0
      S   = 1.D0
      I   = 1
      J   = 2
      K   = 3

    1 DO IPERM=1,3
        DEN = DEN + S*A(1,I)*A(2,J)*A(3,K)
        L   = I
        I   = J
        J   = K
        K   = L
      END DO

      I = 2
      J = 1
      K = 3
      S = - S
      IF(S.LT.0.D0) GO TO 1

      I = 1
      J = 2
      K = 3

      DO IR=1,3
        B(IR,1) = (A(2,J)*A(3,K) - A(2,K)*A(3,J)) / DEN
        B(IR,2) = (A(3,J)*A(1,K) - A(3,K)*A(1,J)) / DEN
        B(IR,3) = (A(1,J)*A(2,K) - A(1,K)*A(2,J)) / DEN
        L = I
        I = J
        J = K
        K = L
      END DO

      DO KK=1,9
        C(KK,1) = B(KK,1)
      END DO

      RETURN
      END
