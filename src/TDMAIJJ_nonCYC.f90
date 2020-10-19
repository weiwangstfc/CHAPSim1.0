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
SUBROUTINE TDMAIJJ_nonCYC(A, B, C, R, BCJ, JS, JSZ, IS, ISZ)  !(A,B,C,R, N,UUU, M)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: JS, JSZ
    INTEGER(4), INTENT(IN) :: IS, ISZ
    REAL(WP), INTENT(IN) :: A(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(IN) :: B(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(IN) :: C(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(INOUT) :: R(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP), INTENT(IN) :: BCJ(IS : ISZ + IS - 1, 2)

    REAL(WP) :: H(IS : ISZ + IS - 1, JS : JSZ + JS - 1)
    REAL(WP) :: G(IS : ISZ + IS - 1, JS : JSZ + JS - 1)

    REAL(WP) :: PP, QQ
    INTEGER(4) :: I, IEND
    INTEGER(4) :: J, JEND

    H = 0.0_WP
    G = 0.0_WP
    PP = 0.0_WP
    QQ = 0.0_WP
    IEND = IS + ISZ - 1
    JEND = JS + JSZ - 1

    J = JS
    DO I = IS, IEND
        H(I, J) = C(I, J) / B(I, J)
        G(I, J) = ( R(I, J) - A(I, J) * BCJ(I, 1) ) / B(I, J)
    END DO


    DO J = JS + 1, JEND - 1
        DO I = IS, IEND
            PP = 1.0_WP / (B(I, J) - H(I, J - 1) * A(I, J) )
            H(I, J) = C(I, J) * PP
            G(I, J) = ( R(I, J) - A(I, J) * G(I, J - 1) ) * PP
        END DO
    END DO

    J = JEND
    DO I = IS, IEND
        PP = 1.0_WP / (B(I, J) - H(I, J - 1) * A(I, J) )
        G(I, J) = ( R(I, J) - C(I, J) * BCJ(I, 2) - A(I, J) * G(I, J - 1) ) * PP
        R(I, J) = G(I, J)
    END DO


    DO J = JEND - 1, JS, -1
        DO I = IS, IEND
            R(I, J) = G(I, J) - H(I, J) * R(I, J + 1)
        END DO
    END DO

    RETURN
END SUBROUTINE
