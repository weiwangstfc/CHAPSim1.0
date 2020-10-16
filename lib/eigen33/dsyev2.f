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
      SUBROUTINE DSYEV2(A, B, C, RT1, RT2, CS, SN)
* ----------------------------------------------------------------------------
* Calculates the eigensystem of a real symmetric 2x2 matrix
*    [ A  B ]
*    [ B  C ]
* in the form
*    [ A  B ]  =  [ CS  -SN ] [ RT1   0  ] [  CS  SN ]
*    [ B  C ]     [ SN   CS ] [  0   RT2 ] [ -SN  CS ]
* where RT1 >= RT2. Note that this convention is different from the one used
* in the LAPACK routine DLAEV2, where |RT1| >= |RT2|.
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A, B, C, RT1, RT2, CS, SN

*     .. Local Variables ..
      DOUBLE PRECISION SM, DF, RT, T

      SM = A + C
      DF = A - C
      RT = SQRT(DF**2 + 4.0D0 * B**2)

*     Calculate eigenvalues
      IF (SM .GT. 0.0D0) THEN
        RT1 = 0.5D0 * (SM + RT)
        T   = 1.0D0 / RT1
        RT2 = (A*T)*C - (B*T)*B
      ELSE IF (SM .LT. 0.0D0) THEN
        RT2 = 0.5D0 * (SM - RT)
        T   = 1.0D0 / RT2
        RT1 = (A*T)*C - (B*T)*B
      ELSE
*       This case needs to be treated separately to avoid DIV by 0
        RT1 = 0.5D0 * RT
        RT2 = -0.5D0 * RT
      END IF
      
*     Calculate eigenvectors
      IF (DF .GT. 0.0D0) THEN
        CS = DF + RT
      ELSE
        CS = DF - RT
      END IF

      IF (ABS(CS) .GT. 2.0D0 * ABS(B)) THEN
        T  = -2.0D0 * B / CS
        SN = 1.0D0 / SQRT(1.0D0 + T**2)
        CS = T * SN
      ELSE IF (ABS(B) .EQ. 0.0D0) THEN
        CS = 1.0D0
        SN = 0.0D0
      ELSE
        T  = -0.5D0 * CS / B
        CS = 1.0D0 / SQRT(1.0D0 + T**2)
        SN = T * CS
      END IF

      IF (DF .GT. 0.0D0) THEN
        T  = CS
        CS = -SN
        SN = T
      END IF
      
      END SUBROUTINE
* End of subroutine DSYEV2

