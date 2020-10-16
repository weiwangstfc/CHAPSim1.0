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
      SUBROUTINE ZHEEVC3(A, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues of a hermitian 3x3 matrix A using Cardano's
* analytical algorithm.
* Only the diagonal and upper triangular parts of A are accessed. The access
* is read-only.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The hermitian input matrix
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
* Dependencies:
*   SQRABS()
* ---------------------------------------------------------------------------
*     .. Arguments ..
      COMPLEX*16       A(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      DOUBLE PRECISION SQRT3
      PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )

*     .. Local Variables ..
      DOUBLE PRECISION M, C1, C0
      COMPLEX*16       DE
      DOUBLE PRECISION DD, EE, FF
      DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
  
*     .. External Functions ..
      DOUBLE PRECISION SQRABS
      EXTERNAL         SQRABS

*     Determine coefficients of characteristic poynomial. We write
*           | A   D   F  |
*      A =  | D*  B   E  |
*           | F*  E*  C  |
      DE    = A(1,2) * A(2,3)
      DD    = SQRABS(A(1,2))
      EE    = SQRABS(A(2,3))
      FF    = SQRABS(A(1,3))
      M     = DREAL(A(1,1)) + DREAL(A(2,2)) + DREAL(A(3,3))
      C1    = ( DREAL(A(1,1)) * DREAL(A(2,2))
     $         + DREAL(A(1,1)) * DREAL(A(3,3))
     $         + DREAL(A(2,2)) * DREAL(A(3,3)) )
     $         - (DD + EE + FF)
      C0    = DREAL(A(3,3))*DD + DREAL(A(1,1))*EE + DREAL(A(2,2))*FF
     $         - DREAL(A(1,1))*DREAL(A(2,2))*DREAL(A(3,3))
     $         - 2.0D0*(DREAL(A(1,3))*DREAL(DE)+DIMAG(A(1,3))*DIMAG(DE))

      P     = M**2 - 3.0D0 * C1
      Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1)
     $          + C0 * (Q + (27.0D0/4.0D0)*C0) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0D0/3.0D0) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE
* End of subroutine ZHEEVC3

