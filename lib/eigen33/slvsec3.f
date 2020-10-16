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
      SUBROUTINE SLVSEC3(D, Z, W, R, I1, I2, I3)
* ----------------------------------------------------------------------------
* Finds the three roots lambda_j of the secular equation
*   f(W_j) = 1 + Sum[ Z_i / (D_i - W_j) ]  ==  0.
* It is assumed that D_0 <= D_1 <= D_2, and that all Z_i have the same sign.
* The arrays R will contain the information required for the calculation
* of the eigenvectors:
*   R_ij = D_i - W_j.
* These differences can be obtained with better accuracy from intermediate
* results.
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION D(3), Z(3), W(3)
      DOUBLE PRECISION R(3,3)
      INTEGER          I1, I2, I3

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
      DOUBLE PRECISION SQRT3
      PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

*     .. Local Variables ..
      DOUBLE PRECISION A(4)
      DOUBLE PRECISION DELTA
      DOUBLE PRECISION DD(3)
      DOUBLE PRECISION XL, XH, X, X0(3)
      DOUBLE PRECISION DX, DXOLD
      DOUBLE PRECISION F, DF
      DOUBLE PRECISION ERROR
      DOUBLE PRECISION T(3)
      DOUBLE PRECISION ALPHA, BETA, GAMMA
      DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
      INTEGER          I, J, NITER

*     Determine intervals which must contain the roots
      IF (Z(1) .GT. 0.0D0) THEN
        A(1) = D(I1)
        A(2) = D(I2)
        A(3) = D(I3)
        A(4) = ABS(D(1) + 3.0D0*Z(1))
     $           + ABS(D(2) + 3.0D0*Z(2))
     &           + ABS(D(3) + 3.0D0*Z(3))
      ELSE
        A(1) = -ABS(D(1) + 3.0D0*Z(1))
     $           - ABS(D(2) + 3.0D0*Z(2))
     &           - ABS(D(3) + 3.0D0*Z(3))
        A(2) = D(I1)
        A(3) = D(I2)
        A(4) = D(I3)
      END IF

*     Calculate roots of f(x) = 0 analytically (analogous to ZHEEVC3)
      T(1)  = D(2) * D(3)
      T(2)  = D(1) * D(3)
      T(3)  = D(1) * D(2)
      GAMMA = T(1) * D(1) + ( Z(1)*T(1) + Z(2)*T(2) + Z(3)*T(3) )
      BETA  = ( Z(1) * (D(2) + D(3) ) + Z(2) * (D(1) + D(3) )
     $            + Z(3) * (D(1) + D(2)) )
     $            + ( T(1) + T(2) + T(3) )
      ALPHA = ( Z(1) + Z(2) + Z(3) ) + ( D(1) + D(2) + D(3) )

      P     = ALPHA**2 - 3.0D0 * BETA
      Q     = ALPHA * (P - (3.0D0/2.0D0) * BETA) + (27.0D0/2.0D0)*GAMMA
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0D0 * ( 0.25D0 * BETA**2 * (P - BETA)
     $                       - GAMMA * (Q - (27.0D0/4.0D0) * GAMMA) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)
      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * ABS(SIN(PHI))

*     Make sure the roots are in ascending order
      X0(1) = (1.0D0/3.0D0) * (ALPHA - C)
      X0(2) = X0(1)
      X0(3) = X0(1)
      IF (C .GT. S) THEN
        X0(1) = X0(1) - S
        X0(2) = X0(2) + S
        X0(3) = X0(3) + C
      ELSE IF (C .LT. -S) THEN
        X0(1) = X0(1) + C
        X0(2) = X0(2) - S
        X0(3) = X0(3) + S
      ELSE
        X0(1) = X0(1) - S
        X0(2) = X0(2) + C
        X0(3) = X0(3) + S
      END IF

*     Refine roots with a combined Bisection/Newton-Raphson method
      DO 10 I = 1, N
        XL    = A(I)
        XH    = A(I+1)
        DX    = 0.5D0 * (XH - XL)
        DXOLD = DX

*       Make sure that XL != XH
        IF (DX .EQ. 0.0D0) THEN
          W(I) = XL
          DO 15 J = 1, N
            R(J, I) = D(J) - XL
   15     CONTINUE
          GO TO 10
        END IF

*       Shift the root close to zero to achieve better accuracy
        IF (X0(I) .GE. XH) THEN
          DELTA = XH
          X     = -DX
          DO 20 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J) - X
   20     CONTINUE
        ELSE IF (X0(I) .LE. XL) THEN
          DELTA = XL
          X     = DX
          DO 30 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J) - X
   30     CONTINUE
        ELSE
          DELTA = X0(I)
          X     = 0.0D0
          DO 40 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J)
   40     CONTINUE
        END IF
        XL = XL - DELTA
        XH = XH - DELTA

*       Make sure that f(XL) < 0 and f(XH) > 0
        IF (Z(1) .LT. 0.0D0) THEN
          P  = XH
          XH = XL
          XL = P
        END IF

*       Main iteration loop
        DO 50 NITER = 1, 500
*         Evaluate f and f', and calculate an error estimate
          F     = 1.0D0
          DF    = 0.0D0
          ERROR = 1.0D0
          DO 60 J = 1, N
            T(1)  = 1.0D0 * (1.0 / R(J, I))
            T(2)  = Z(J) * T(1)
            T(3)  = T(2) * T(1)
            F     = F + T(2)
            ERROR = ERROR + ABS(T(2))
            DF    = DF + T(3)
   60     CONTINUE

*         Check for convergence and leave loop if applicable
          IF (ABS(F) .LE. EPS * (8.0D0 * ERROR + ABS(X * DF))) THEN
            GO TO 70
          END IF

*         Adjust interval boundaries
          IF (F .LT. 0.0D0) THEN
            XL   = X
          ELSE
            XH   = X
          END IF

*         Check, whether Newton-Raphson would converge fast enough.
*         If so, give it a try. If not, or if it would run out of
*         bounds, use bisection.
          IF ( ABS(2.0D0 * F) .LT. ABS(DXOLD * DF) ) THEN
            DXOLD = DX
            DX    = F * (1.0 / DF)
            X     = X - DX
            IF ( (X - XH) * (X - XL) .GE. 0.0D0 ) THEN
              DX = 0.5D0 * (XH - XL)
              X  = XL + DX
            END IF
          ELSE
            DX = 0.5D0 * (XH - XL)
            X  = XL + DX
          END IF
            
*         Prepare next iteration
          DO 80 J = 1, N
            R(J, I) = DD(J) - X
   80     CONTINUE
   50   CONTINUE

*       Un-shift result
   70   W(I) = X + DELTA
   10 CONTINUE

      END SUBROUTINE
* End of subroutine SLVSEC3

 

