* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVD3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric
* 3x3 matrix A using Cuppen's Divide & Conquer algorithm.
* The function accesses only the diagonal and upper triangular parts of
* A. The access is read-only. 
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
* Dependencies:
*   DSYEV2(), SLVSEC3(), DSYTRD3()
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

*     .. Local Variables ..
      DOUBLE PRECISION R(3,3)
      DOUBLE PRECISION P(3,3)
      DOUBLE PRECISION D(3), E(2), Z(3)
      DOUBLE PRECISION C, S, T
      INTEGER          I, J, K

*     .. External Functions ..
      EXTERNAL         DSYEV2, SLVSEC3, DSYTRD3
      
*     Initialize Q
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 I = 1, N
        DO 11 J = 1, N
          Q(I, J) = 0.0D0
   11   CONTINUE
   10 CONTINUE
  
*     Transform A to real tridiagonal form by the Householder method
      CALL DSYTRD3(A, R, W, E)

      
*     "Divide"
*     -------------------------------
  
*     Detect matrices that factorize to avoid multiple eigenvalues in the Divide/Conquer algorithm
      DO 20 I = 1, N-1
        T = ABS(W(I)) + ABS(W(I+1))
        IF (ABS(E(I)) .LE. 8.0D0 * EPS * T) THEN
          IF (I .EQ. 1) THEN
            CALL DSYEV2(W(2), E(2), W(3), D(2), D(3), C, S)
            W(2) = D(2)
            W(3) = D(3)
*           --- The rest of this IF-branch can be omitted if only the eigenvalues are desired ---
            Q(1, 1) = 1.0D0
            DO 30 J = 2, N
              Q(J, 2) = S * R(J, 3) + C * R(J, 2)
              Q(J, 3) = C * R(J, 3) - S * R(J, 2)
   30       CONTINUE
          ELSE
            CALL DSYEV2(W(1), E(1), W(2), D(1), D(2), C, S)
            W(1)    = D(1)
            W(2)    = D(2)
*           --- The rest of this ELSE-branch can be omitted if only the eigenvalues are desired ---
            Q(1, 1) = C
            Q(1, 2) = -S
            Q(2, 1) = R(2, 2) * S
            Q(2, 2) = R(2, 2) * C
            Q(2, 3) = R(2, 3)
            Q(3, 1) = R(3, 2) * S
            Q(3, 2) = R(3, 2) * C
            Q(3, 3) = R(3, 3)
          END IF
          RETURN
        END IF
   20 CONTINUE

*     Calculate eigenvalues and eigenvectors of 2x2 block
      CALL DSYEV2(W(2) - E(1), E(2), W(3), D(2), D(3), C, S)
      D(1) = W(1) - E(1)

      
*     "Conquer"
*     -------------------------------

*     Determine coefficients of secular equation
      Z(1) = E(1)
      Z(2) = E(1) * C**2
      Z(3) = E(1) * S**2

*     Call SLVSEC3 with D sorted in ascending order. We make use of the
*     fact that DSYEV2 guarantees D[1] >= D[2].
      IF (D(1) .LT. D(3)) THEN
        CALL SLVSEC3(D, Z, W, P, 1, 3, 2)
      ELSE IF (D(1) .LT. D(2)) THEN
        CALL SLVSEC3(D, Z, W, P, 3, 1, 2)
      ELSE
        CALL SLVSEC3(D, Z, W, P, 3, 2, 1)
      END IF

*     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

*     Calculate eigenvectors of matrix D + BETA * Z * Z^T and store them
*     in the columns of P
      Z(1) = SQRT(ABS(E(1)))
      Z(2) = C * Z(1)
      Z(3) = -S * Z(1)

*     Detect duplicate elements in D to avoid division by zero
      T = 8.0D0 * EPS * ( ABS(D(1)) + ABS(D(2)) + ABS(D(3)) )
      IF (ABS(D(2) - D(1)) .LE. T) THEN
        DO 40 J = 1, N
          IF (P(1, J) * P(2, J) .LE. 0.0D0) THEN
            P(1, J) = Z(2)
            P(2, J) = -Z(1)
            P(3, J) = 0.0D0
          ELSE
            DO 45 I = 1, N
              P(I, J) = Z(I) / P(I, J)
   45       CONTINUE
          END IF
   40   CONTINUE
      ELSE IF (ABS(D(3) - D(1)) .LE. T) THEN
        DO 50 J = 1, N
          IF (P(1, J) * P(3, J) .LE. 0.0D0) THEN
            P(1, J) = Z(3)
            P(2, J) = 0.0D0
            P(3, J) = -Z(1)
          ELSE
            DO 55 I = 1, N
              P(I, J) = Z(I) / P(I, J)
   55       CONTINUE
          END IF
   50   CONTINUE
      ELSE
        DO 60 J = 1, N
          DO 61 I = 1, N
            IF (P(I, J) .EQ. 0.0) THEN
              P(I, J) = 1.0
              P(1 + MOD(I, N), J)   = 0.0
              P(1 + MOD(I+1, N), J) = 0.0
              GO TO 60
            ELSE
              P(I, J) = Z(I) / P(I, J)
            END IF
   61     CONTINUE
   60   CONTINUE
      END IF

*     Normalize eigenvectors of D + BETA * Z * Z^T
      DO 70 J = 1, N
        T = P(1, J)**2 + P(2, J)**2 + P(3, J)**2
        T = 1.0D0 / SQRT(T)
        DO 75 I = 1, N
          P(I, J) = P(I, J) * T
   75   CONTINUE
   70 CONTINUE
  
*     Undo diagonalization of 2x2 block
      DO 80 J = 1, N
        T       = P(2, J)
        P(2, J) = C * T - S * P(3, J)
        P(3, J) = S * T + C * P(3, J)
   80 CONTINUE

*     Undo Householder transformation
      DO 90 J = 1, N
        DO 91 K = 1, N
          T       = P(K, J)
          DO 95 I = 1, N
            Q(I, J) = Q(I, J) + T * R(I, K)
   95     CONTINUE
   91   CONTINUE
   90 CONTINUE

      END SUBROUTINE
* End of subroutine DSYEVD3

