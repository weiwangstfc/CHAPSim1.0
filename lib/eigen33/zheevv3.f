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
      SUBROUTINE ZHEEVV3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
* matrix A using Cardano's method for the eigenvalues and an analytical
* method based on vector cross products for the eigenvectors.
* Only the diagonal and upper triangular parts of A need to contain
* meaningful values. However, all of A may be used as temporary storage
* and may hence be destroyed.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The hermitian input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
* Dependencies:
*   SQRABS(), ZHEEVC3()
* ----------------------------------------------------------------------------
*     .. Arguments ..
      COMPLEX*16       A(3,3)
      COMPLEX*16       Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

*     .. Local Variables ..
      DOUBLE PRECISION NORM, N1, N2, N1TMP, N2TMP
      DOUBLE PRECISION THRESH, ERROR, WMAX, T
      COMPLEX*16       F
      INTEGER          I, J

*     .. External Functions ..
      DOUBLE PRECISION SQRABS
      EXTERNAL         SQRABS, ZHEEVC3

*     Calculate eigenvalues
      CALL ZHEEVC3(A, W)

*     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

      WMAX   = MAX(ABS(W(1)), ABS(W(2)), ABS(W(3)))
      THRESH = (8.0D0 * EPS * WMAX)**2

*     Prepare calculation of eigenvectors
      N1TMP   = SQRABS(A(1, 2)) + SQRABS(A(1, 3))
      N2TMP   = SQRABS(A(1, 2)) + SQRABS(A(2, 3))
      Q(1, 1) = A(1, 2) * A(2, 3) - A(1, 3) * DREAL(A(2, 2))
      Q(1, 2) = Q(1, 1)
      Q(2, 1) = A(1, 3) * DCONJG(A(1, 2)) - A(2, 3) * DREAL(A(1, 1))
      Q(2, 2) = Q(2, 1)
      Q(3, 2) = SQRABS(A(1, 2))

*     Calculate first eigenvector by the formula
*       v[0] = conj( (A - lambda[0]).e1 x (A - lambda[0]).e2 )
      A(1, 1) = A(1, 1) - W(1)
      A(2, 2) = A(2, 2) - W(1)
      Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
      Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
      Q(3, 1) = DREAL(A(1, 1)) * DREAL(A(2, 2)) - Q(3, 2)
      NORM    = SQRABS(Q(1, 1)) + SQRABS(Q(2, 1)) + SQRABS(Q(3, 1))
      N1      = N1TMP + DREAL(A(1, 1))**2
      N2      = N2TMP + DREAL(A(2, 2))**2
      ERROR   = N1 * N2

*     If the first column is zero, then (1, 0, 0) is an eigenvector
      IF (N1 .LE. THRESH) THEN
        Q(1, 1) = 1.0D0
        Q(2, 1) = 0.0D0
        Q(3, 1) = 0.0D0
*     If the second column is zero, then (0, 1, 0) is an eigenvector
      ELSE IF (N2 .LE. THRESH) THEN
        Q(1, 1) = 0.0D0
        Q(2, 1) = 1.0D0
        Q(3, 1) = 0.0D0
*     If angle between A(*,1) and A(*,2) is too small, don't use
*     cross product, but calculate v ~ (1, -A0/A1, 0)
      ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
        T = SQRABS(A(1, 2))
        F = -A(1, 1) / A(1, 2)
        IF (SQRABS(A(2, 2)) .GT. T) THEN
          T = SQRABS(A(2, 2))
          F = -DCONJG(A(1, 2)) / A(2, 2)
        END IF
        IF (SQRABS(A(2, 3)) .GT. T) THEN
          F = -DCONJG(A(1, 3) / A(2, 3))
        END IF
        NORM    = 1.0D0 / SQRT(1.0D0 + SQRABS(F))
        Q(1, 1) = NORM
        Q(2, 1) = F * NORM
        Q(3, 1) = 0.0D0
*     This is the standard branch
      ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 20, J = 1, 3
          Q(J, 1) = Q(J, 1) * NORM
   20   CONTINUE
      END IF
 
*     Prepare calculation of second eigenvector     
      T = W(1) - W(2)

*     Is this eigenvalue degenerate?
      IF (ABS(T) .GT. 8.0D0 * EPS * WMAX) THEN
*       For non-degenerate eigenvalue, calculate second eigenvector by
*       the formula
*         v[1] = conj( (A - lambda[1]).e1 x (A - lambda[1]).e2 )
        A(1, 1) = A(1, 1) + T
        A(2, 2) = A(2, 2) + T
        Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
        Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
        Q(3, 2) = DREAL(A(1, 1)) * DREAL(A(2, 2)) - DREAL(Q(3, 2))
        NORM    = SQRABS(Q(1, 2)) + SQRABS(Q(2, 2)) + DREAL(Q(3, 2))**2
        N1      = N1TMP + DREAL(A(1, 1))**2
        N2      = N2TMP + DREAL(A(2, 2))**2
        ERROR   = N1 * N2

        IF (N1 .LE. THRESH) THEN
          Q(1, 2) = 1.0D0
          Q(2, 2) = 0.0D0
          Q(3, 2) = 0.0D0
        ELSE IF (N2 .LE. THRESH) THEN
          Q(1, 2) = 0.0D0
          Q(2, 2) = 1.0D0
          Q(3, 2) = 0.0D0
        ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
          T = SQRABS(A(1, 2))
          F = -A(1, 1) / A(1, 2)
          IF (SQRABS(A(2, 2)) .GT. T) THEN
            T = SQRABS(A(2, 2))
            F = -DCONJG(A(1, 2)) / A(2, 2)
          END IF
          IF (SQRABS(A(2, 3)) .GT. T) THEN
            F = -DCONJG( A(1, 3) / A(2, 3) )
          END IF
          NORM    = 1.0D0 / SQRT(1.0D0 + SQRABS(F))
          Q(1, 2) = NORM
          Q(2, 2) = F * NORM
          Q(3, 2) = 0.0D0
        ELSE
          NORM = SQRT(1.0D0 / NORM)
          DO 40, J = 1, 3 
            Q(J, 2) = Q(J, 2) * NORM
   40     CONTINUE
        END IF
      ELSE
*       For degenerate eigenvalue, calculate second eigenvector according to
*         v[1] = conj( v[0] x (A - lambda[1]).e[i] )
*   
*       This would really get to complicated if we could not assume all of A to
*       contain meaningful values.
        A(2, 1) = DCONJG(A(1, 2))
        A(3, 1) = DCONJG(A(1, 3))
        A(3, 2) = DCONJG(A(2, 3))
        A(1, 1) = A(1, 1) + W(1)
        A(2, 2) = A(2, 2) + W(1)
        DO 50 I = 1, 3
          A(I, I) = A(I, I) - W(2)
          N1      = SQRABS(A(1, I)) + SQRABS(A(2, I)) + SQRABS(A(3, I))
          IF (N1 .GT. THRESH) THEN
            Q(1, 2) = DCONJG( Q(2, 1) * A(3, I) - Q(3, 1) * A(2, I) )
            Q(2, 2) = DCONJG( Q(3, 1) * A(1, I) - Q(1, 1) * A(3, I) )
            Q(3, 2) = DCONJG( Q(1, 1) * A(2, I) - Q(2, 1) * A(1, I) )
            NORM    = SQRABS(Q(1,2)) + SQRABS(Q(2,2)) + SQRABS(Q(3,2))
            IF (NORM .GT. (256.0D0 * EPS)**2 * N1) THEN
              NORM = SQRT(1.0D0 / NORM)
              DO 55 J = 1, 3
                Q(J, 2) = Q(J, 2) * NORM
   55         CONTINUE
              GO TO 60
            END IF
          END IF
   50   CONTINUE
   
*       This means that any vector orthogonal to v[0] is an EV.
   60   IF (I .EQ. 4) THEN
          DO 70 J = 1, 3
*           Find nonzero element of v[0] and swap it with the next one
            IF (DREAL(Q(J, 1)).NE.0.0D0.OR.DIMAG(Q(J, 1)).NE.0.0D0) THEN
              NORM = SQRABS(Q(J, 1)) + SQRABS(Q(1 + MOD(J,3), 1))
              NORM = 1.0D0 / SQRT(NORM)
              Q(J, 2)              = DCONJG(Q(1 + MOD(J,3), 1)) * NORM
              Q(1 + MOD(J,3), 2)   = -DCONJG(Q(J, 1)) * NORM
              Q(1 + MOD(J+1,3), 2) = 0.0D0
              GO TO 80
            END IF
   70     CONTINUE
        END IF
      END IF

*     Calculate third eigenvector according to
*       v[2] = conj(v[0] x v[1])
   80 Q(1, 3) = DCONJG( Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2) )
      Q(2, 3) = DCONJG( Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2) )
      Q(3, 3) = DCONJG( Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2) )

      END SUBROUTINE
* End of subroutine ZHEEVV3

