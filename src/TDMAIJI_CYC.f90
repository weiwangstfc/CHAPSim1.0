!***********************************************************************
!> @author
!> Wei Wang, The University of Sheffield.
!
! DESCRIPTION:
!> Brief description of SUBROUTINE.
!> @details
!>  SOLVER FOR REAL- DEFINITE, CYCLIC- TRIDIAGONAL, CONSTANT COEFFICIENTS
!>  LINEAR SYSTEMS OF EQUATIONS.
!>  THE SOLVER IS WRITTEN IN FORTRAN 90 LANGUAGE AND ADOPTS ONLY
!>  DOUBLE - PRECISION (8 BYTES) REAL variables.
!>********************************************************************************
!>  @par MethoDOlogy                                                             *
!>  THE METHOD SOLVES FIRSTLY FOR THE X(1) UNKNOWN, WHICH MAKES THE              *
!>  COEFFICIENT MATRIX CYCLIC: THIS IS THE INITIALIZATION PHASE, PERFORMED       *
!>  IN THE PRESENT SUBROUTINE. NOW THE SYSTEM IS A (N - 1) BY (N - 1) TRIDIAGONAL    *
!>  (NOT CYCLIC) SYSTEM, AND IS SOLVED with THE STANDARD (AND FAST) LU           *
!>  FACTORIZATION TECHNIQUE (ROUTINES LU_PRE AND LU_SOLV).                       *
!>                                                                               *
!>  COEFFICIENTS ON main AND SECONDARY DIAGONALS MUST BE CONSTANT, I.E. WE are   *
!>  DEALING with SYSTEMS OF THE following KIND:                                  *
!>                                                                               *
!>                      | A  B  0  0  B | |X1|   |D1|                            *
!>                      | B  A  B  0  0 | |X2|   |D2|                            *
!>                      | 0  B  A  B  0 | |X3| = |D3|                            *
!>                      | 0  0  B  A  B | |X4|   |D4|                            *
!>                      | B  0  0  B  A | |X5|   |D5|                            *
!>                                                                               *
!>  OF COURSE, THE IMPROUVEMENT IN RESPECT OF USING THE CLASSIC AND efficient    *
!>  CTDMA SOLVER OF THOMAS ARISES WHEN ONE NEEDS TO SOLVE MANY CYCLIC TRIDIAGONAL *
!>  SYSTEMS with THE SAME COEFFICIENT MATRIX.                                    *
!>                                                                               *
!>  NOTE! IN THE INITIALIZATION PHASE WE PERFORM ALSO THE LU FACTORIZATION       *
!>  OF THE n - 1 BY n - 1 TRIDIAGONAL COEFFICIENT MATRIX.                            *
!>                                                                               *
!>  THE NOn - UNITY DIAGONALS OF MATRICES L & U are STORED IN THE ARRAYS            *
!>  DIAG_PRIN(1 : n - 1) & DIAG_SUPE(1 : n -2) (IN THIS ORDER)                           *
!>                                                                                *
!>  THE RIGHT - HAND SIDE IS STORED IN THE ARRAY RHS (N ELEMENTS). IN               *
!>  OUTPUT FROM SUB. ALG4_SOL RHS CONTAINS THE SOLUTION.                          *
!>                                                                                *
!>  Z IS A REAL *8 ARRAY OF N ELEMENTS, AND doesN'T NEED TO BE INITIALIZED.        *
!>                                                                                *
!> @par USAGE:                                                                    *
!>      CALL TDMAIJI_CYC(A,B,C,R, IS, ISZ, JS, JSZ)                                   *
!>      TRIDIAGNAL MARCH in I direction in the plane of I *J
!>
!> @par INPUT                                                                     *
!>      IS:  INTEGER, The Starting index of I direction in the I *J plane          *
!>      ISZ: INTEGER, The total size of points in I direction                     *
!>      JS:  INTEGER, The Starting index of J direction in the I *J plane          *
!>      JSZ: INTEGER, The total size of points in J direction                     *
!>      A(IS:1:IS + ISZ - 1, JS:1:JS + JSZ - 1) : Matrix ISZ*JSZ, Coef of u_{I - 1, J}       *
!>      B(IS:1:IS + ISZ - 1, JS:1:JS + JSZ - 1) : Matrix ISZ*JSZ, Coef of u_{i, J}         *
!>      C(IS:1:IS + ISZ - 1, JS:1:JS + JSZ - 1) : Matrix ISZ*JSZ, Coef of u_{I + 1, J}       *
!>      R(IS:1:IS + ISZ - 1, JS:1:JS + JSZ - 1) : Matrix ISZ*JSZ, RHS of the original     *
!>                                        tridiagnal sys                          *
!>
!>par   OUTPUT
!>      R(IS:1:IS + ISZ - 1, JS:1:JS + JSZ - 1) : Matrix ISZ*JSZ, Results u_{i, J}         *
!> @todo
!> Nothing left to DO.
!
! REVISION HISTORY:
! 30/01/2013- Initial Version, by Wei Wang
!**********************************************************************************************************************************
SUBROUTINE TDMAIJI_CYC(A,B,C,R, IS, ISZ, JS, JSZ)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IS
    INTEGER(4), INTENT(IN) :: ISZ
    INTEGER(4), INTENT(IN) :: JS
    INTEGER(4), INTENT(IN) :: JSZ
    REAL(WP), INTENT(IN) :: A(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(IN) :: B(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(IN) :: C(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(INOUT) :: R(IS : ISZ + IS - 1, JS : JSZ + JS - 1)

    INTEGER(4) :: I, IEND
    INTEGER(4) :: J, JEND
    REAL(WP) :: PP
    REAL(WP) :: QQ1, QQ2
    REAL(WP) :: H0(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP) :: G1(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP) :: G2(IS : ISZ + IS - 1, JS : JSZ + JS - 1)

    H0 = 0.0_WP
    G1 = 0.0_WP
    G2 = 0.0_WP

    !>    @note Set up the ending index in I and J directions.
    IEND = IS + ISZ - 1
    JEND = JS + JSZ - 1

    !>    @note Set up the initial H1, H2 and G1, G2, in whICh H1 and G1 IS
    !>          For the first Tridiagnal System and H2 and G2 for the second
    !>          One in the sepARated system. (H1 = H2 = H0)

    I = IS
    DO J = JS, JEND
        H0(I, J) = C(I, J) / B(I, J)
        G1(I, J) = R(I, J) / B(I, J)
        G2(I, J) = - A(I, J) / B(I, J)
    END DO

    !>    @note forward substitution to get H, G1 and G2, from I = 2 to I = n -2
    DO I = IS + 1, IEND - 2
        DO J = JS, JEND
            PP = 1.0_WP / (B(I, J) - H0(I - 1, J) * A(I, J) )
            QQ1 = G1(I - 1, J) * A(I, J)
            QQ2 = G2(I - 1, J) * A(I, J)
            H0(I, J) = C(I, J) * PP
            G1(I, J) = ( R(I, J) - QQ1 ) * PP
            G2(I, J) = -QQ2 * PP
        END DO
    END DO

    !>    @note For the last point H0, G1, G2 at n - 1
    I = IEND - 1
    DO J = JS, JEND
        PP = 1.0_WP / (B(I, J) - H0(I - 1, J) * A(I, J) )
        QQ1 = G1(I - 1, J) * A(I, J)
        QQ2 = G2(I - 1, J) * A(I, J)
        G1(I, J) = (  R(I, J) - QQ1 ) * PP
        G2(I, J) = ( -C(I, J) - QQ2 ) * PP
    END DO
    !

    !>    @note backward substitution to calcute U1, U2 from I = n - 1 to 1
    !>          G1 and G2 IS overWRitten by the results U1 and U2
    DO I = IEND - 2, IS, -1
        DO J = JS, JEND
            G1(I, J) = - H0(I, J) * G1(I + 1, J) + G1(I, J)
            G2(I, J) = - H0(I, J) * G2(I + 1, J) + G2(I, J)
        END DO
    END DO

    !>    @note Calcuate the original unknowns, x_n
    !>         R IS rE -WRitten by the results U0.
    I = IEND
    DO J = JS, JEND
        PP = R(I, J) - C(I, J) * G1(IS, J) - A(I, J) * G1(I - 1, J)
        QQ1 = B(I, J) + C(I, J) * G2(IS, J) + A(I, J) * G2(I - 1, J)
        R(I, J) = PP / QQ1
    END DO

    !>    @note backward substitusion from I = n - 1 to 1
    !>         R IS rE -WRitten by the results U0.
    DO I = IEND - 1, IS, -1
        DO J = JS, JEND
            R(I, J) = G1(I, J) + G2(I, J) * R(IEND, J)
        END DO
    END DO

    RETURN
END SUBROUTINE
