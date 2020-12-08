!**********************************************************************************************************************************
!> @brief
!>        functions to calculate thermo-properties from table.
!> @details
!> function: spline_interpolation_TD
!> function: spline_interpolation_TM
!> function: spline_interpolation_TK
!> function: spline_interpolation_TB
!> function: spline_interpolation_TCp
!> function: spline_interpolation_HT
!> function: spline_interpolation_HD
!> function: spline_interpolation_HM
!> function: spline_interpolation_HK
!> function: spline_interpolation_HB
!> function: spline_interpolation_HCp
!> function: spline_interpolation_DHH
!> function: spline_interpolation_DHT
!> SUBROUTINE: BiSection_LIST_H
!> SUBROUTINE: BiSection_LIST_T
!> SUBROUTINE: SPLINE_COEFFICIENTS_FOR_DH
!> SUBROUTINE: SPLINE_COEFFICIENTS_FOR_H
!> SUBROUTINE: SPLINE_COEFFICIENTS_FOR_T
!> SUBROUTINE: CUBIC_SPLINE
!> SUBROUTINE: FIND_LOCAL_MAXIMA
!> SUBROUTINE: Property_Table_Nondimensionalization
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
! 09/2020 - Added more fluid types and optimized, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
FUNCTION spline_interpolation_TH(T)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TH

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_H(1) - LIST_H(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_H(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_H(N) - LIST_H(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_H(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_H(I) + DX * (SplineCoeff_TH_B(I) + DX * (SplineCoeff_TH_C(I) + DX * SplineCoeff_TH_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 0.0_WP !test
    spline_interpolation_TH = EVAL
END FUNCTION

!**********************************************************************************************************************************
FUNCTION spline_interpolation_TD(T)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TD

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_D(1) - LIST_D(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_D(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_D(N) - LIST_D(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_D(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_D(I) + DX * (SplineCoeff_TD_B(I) + DX * (SplineCoeff_TD_C(I) + DX * SplineCoeff_TD_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_TD = EVAL
END FUNCTION

!**********************************************************************************************************************************
FUNCTION spline_interpolation_TM(T)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TM

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_M(1) - LIST_M(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_M(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_M(N) - LIST_M(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_M(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_M(I) + DX * (SplineCoeff_TM_B(I) + DX * (SplineCoeff_TM_C(I) + DX * SplineCoeff_TM_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_TM = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_TK(T)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TK

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_K(1) - LIST_K(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_K(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_K(N) - LIST_K(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_K(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_K(I) + DX * (SplineCoeff_TK_B(I) + DX * (SplineCoeff_TK_C(I) + DX * SplineCoeff_TK_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_TK = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_TB(T)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TB

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_B(1) - LIST_B(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_B(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_B(N) - LIST_B(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_B(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_B(I) + DX * (SplineCoeff_TB_B(I) + DX * (SplineCoeff_TB_C(I) + DX * SplineCoeff_TB_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_TB = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_TCp(T)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: T, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_TCp

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF(T < LIST_T(1)) THEN
        !EVAL = ( LIST_Cp(1) - LIST_Cp(2) ) / ( LIST_T(1) - LIST_T(2) ) * (T - LIST_T(1)) + LIST_Cp(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its minimum value in the given NIST table.', LIST_T(1), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmin =, T = ", MYID, LIST_T(1), T)
        STOP
    END IF
    IF(T > LIST_T(N)) THEN
        !EVAL = ( LIST_Cp(N) - LIST_Cp(N - 1) ) / ( LIST_T(N) - LIST_T(N - 1) ) * (T - LIST_T(N)) + LIST_Cp(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: T reaches its maximum value in the given NIST table.', LIST_T(N), T
        CALL CHK2RLHDL("Error! T reaches its limit. Tmax =, T = ", MYID, LIST_T(N), T)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(T < LIST_T(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = T - LIST_T(I)
    EVAL = LIST_Cp(I) + DX * (SplineCoeff_TCp_B(I) + DX * (SplineCoeff_TCp_C(I) + DX * SplineCoeff_TCp_D(I)))
    !IF(DABS(T - 1.0_WP) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_TCp = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HT(H)
    !======================================================================
    ! FUNCTION ispline evaluates the CUBIC spline interpolation at point z
    ! ispline = y(I) +b(I) * (u-x(I)) +c(I) * (u-x(I))**2 + D(I) * (u-x(I))**3
    ! where  x(I) <= u <= x(I + 1)
    !----------------------------------------------------------------------
    ! input..
    ! u     = the abscISsa at whICh the spline IS to be evaluated
    ! x, y  = the ARrays of given data points
    ! b, c, d = ARrays of spline coefficients computed by spline
    ! n     = the number of data points
    ! output:
    ! ispline = INTerpolated value at point u
    !=======================================================================
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HT

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_T(1) - LIST_T(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_T(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_T(N) - LIST_T(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_T(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(h < LIST_H(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = h - LIST_H(I)
    EVAL = LIST_T(I) + DX * (SplineCoeff_HT_B(I) + DX * (SplineCoeff_HT_C(I) + DX * SplineCoeff_HT_D(I)))

    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HT = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HD(H)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HD

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_D(1) - LIST_D(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_D(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_D(N) - LIST_D(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_D(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(h < LIST_H(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = h - LIST_H(I)
    EVAL = LIST_D(I) + DX * (SplineCoeff_HD_B(I) + DX * (SplineCoeff_HD_C(I) + DX * SplineCoeff_HD_D(I)))
    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HD = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HM(H)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HM

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_M(1) - LIST_M(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_M(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_M(N) - LIST_M(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_M(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(h < LIST_H(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = h - LIST_H(I)
    EVAL = LIST_M(I) + DX * (SplineCoeff_HM_B(I) + DX * (SplineCoeff_HM_C(I) + DX * SplineCoeff_HM_D(I)))
    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HM = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HK(H)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HK

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_K(1) - LIST_K(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_K(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_K(N) - LIST_K(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_K(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(h < LIST_H(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = h - LIST_H(I)
    EVAL = LIST_K(I) + DX * (SplineCoeff_HK_B(I) + DX * (SplineCoeff_HK_C(I) + DX * SplineCoeff_HK_D(I)))
    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HK = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HB(H)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HB

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_B(1) - LIST_B(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_B(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_B(N) - LIST_B(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_B(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (j > I + 1)
        k = (I + J) / 2
        IF(h < LIST_H(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = h - LIST_H(I)
    EVAL = LIST_B(I) + DX * (SplineCoeff_HB_B(I) + DX * (SplineCoeff_HB_C(I) + DX * SplineCoeff_HB_D(I)))
    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HB = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_HCp(H)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: H, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_HCp

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (H - LIST_H(IMIN_H)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !EVAL = ( LIST_Cp(1) - LIST_Cp(2) ) / ( LIST_H(1) - LIST_H(2) ) * (H- LIST_H(1)) + LIST_Cp(1)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its minimum value in the given NIST table.', LIST_H(1), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmin =, H = ", MYID, LIST_H(1), H)
        STOP
    END IF
    IF( (H - LIST_H(IMAX_H)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !EVAL = ( LIST_Cp(N) - LIST_Cp(N - 1) ) / ( LIST_H(N) - LIST_H(N - 1) ) * (H- LIST_H(N)) + LIST_Cp(N)
        !WRITE(logflg_io, '(A, 2F12.5)') '# WARning: H reaches its maximum value in the given NIST table.', LIST_H(N), H
        CALL CHK2RLHDL("Error! H reaches its limit. Hmax =, H = ", MYID, LIST_H(N), H)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (J > I + 1)
        K = (I + J) / 2
        IF(H < LIST_H(K)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = H - LIST_H(I)
    EVAL = LIST_Cp(I) + DX * (SplineCoeff_HCp_B(I) + DX * (SplineCoeff_HCp_C(I) + DX * SplineCoeff_HCp_D(I)))
    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_HCp = EVAL
END FUNCTION
!**********************************************************************************************************************************
FUNCTION spline_interpolation_DHH(RH)!by Junjie, 2017/03 / 26
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: RH, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_DHH

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (RH - LIST_DH(IMIN_DH)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !CALL CHKHDL('min RH, MYID = ', MYID)
        !EVAL = ( LIST_DH(1) - LIST_DH(2) ) / ( LIST_RH(1) - LIST_RH(2) ) * (RH- LIST_RH(1)) + LIST_DH(1)
        CALL CHK2RLHDL("Error! DH reaches its limit in spline_interpolation_DHH. DHmin =, DH = ", MYID, LIST_DH(IMIN_DH), RH)
        STOP
    END IF
    IF( (RH - LIST_DH(IMAX_DH)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !CALL CHKHDL('max RH, MYID = ', MYID)
        !EVAL = ( LIST_DH(N) - LIST_DH(N - 1) ) / ( LIST_RH(N) - LIST_RH(N - 1) ) * (RH- LIST_RH(N)) + LIST_DH(N)
        CALL CHK2RLHDL("Error! DH reaches its limit in spline_interpolation_DHH. DHmax =, DH = ", MYID, LIST_DH(IMAX_DH), RH)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (J > I + 1)
        K = (I + J) / 2
        IF(RH < LIST_DH(K)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = RH - LIST_DH(I)
    EVAL = LIST_H(I) + DX * (SplineCoeff_DHH_B(I) + DX * (SplineCoeff_DHH_C(I) + DX * SplineCoeff_DHH_D(I)))

    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_DHH = EVAL
END FUNCTION

!**********************************************************************************************************************************
FUNCTION spline_interpolation_DHT(RH)
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    REAL(WP) :: RH, DX, EVAL
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: spline_interpolation_DHT

    N = N_LIST
    ! IF u IS ouside the x() interval take a boundary value (left or right)
    IF( (RH - LIST_DH(IMIN_DH)) < -1.0E-8_WP ) THEN
        !EVAL = LIST_T(1)
        !CALL CHKHDL('min RH, MYID = ', MYID)
        !EVAL = ( LIST_DH(1) - LIST_DH(2) ) / ( LIST_RH(1) - LIST_RH(2) ) * (RH- LIST_RH(1)) + LIST_DH(1)
        CALL CHK2RLHDL("Error! DH reaches its limit in spline_interpolation_DHT. DHmin =, DH = ", MYID, LIST_DH(IMIN_DH), RH)
        STOP
    END IF
    IF( (RH - LIST_DH(IMAX_DH)) > 1.0E-8_WP ) THEN
        !EVAL = LIST_T(N)
        !CALL CHKHDL('max RH, MYID = ', MYID)
        !EVAL = ( LIST_DH(N) - LIST_DH(N - 1) ) / ( LIST_RH(N) - LIST_RH(N - 1) ) * (RH- LIST_RH(N)) + LIST_DH(N)
        CALL CHK2RLHDL("Error! DH reaches its limit in spline_interpolation_DHT. DHmax =, DH = ", MYID, LIST_DH(IMAX_DH), RH)
        STOP
    END IF

    !*
    !  binary search for for i, such that x(I) <= u <= x(I + 1)
    !*
    I = 1
    J = N + 1
    DO WHILE (J > I + 1)
        K = (I + J) / 2
        IF(RH < LIST_DH(k)) THEN
            J = K
        ELSE
            I = K
        END IF
    END DO
    !*
    !  evaluate spline interpolation
    !*
    DX = RH - LIST_DH(I)
    EVAL = LIST_T(I) + DX * (SplineCoeff_DHT_B(I) + DX * (SplineCoeff_DHT_C(I) + DX * SplineCoeff_DHT_D(I)))

    !IF(DABS(H) < 1.0E-10_WP) EVAL = 1.0_WP !test
    spline_interpolation_DHT = EVAL
END FUNCTION

!**********************************************************************************************************************************
SUBROUTINE BiSection_LIST_H(ENTH, IE, IS)
    USE thermal_info
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: ENTH
    INTEGER(4), INTENT(OUT) :: IE, IS
    REAL(WP) :: VS, VM
    INTEGER(4) :: IM

    !      IF( ( ENTH < LIST_H(IMIN_H) ).OR. ( ENTH  > LIST_H(N_LIST) ) ) THEN
    !        WRITE(*, '(A, 2ES13.5)') 'RANGE OF LIST_H  = ', LIST_H(1), LIST_H(N_LIST)
    !        WRITE(*, '(A, 1ES13.5)') 'current GIVEN ENTHALPY = ', ENTH
    !        CALL ERRHDL(' ENTHOPHY IS OUT OF GIVEN NIST TABLE RANGE.', MYID)
    !      END IF

    IS = 1
    IE = N_LIST
    DO WHILE ( (IE - IS) > 1 )
        IM = IS + (IE -IS) / 2
        VS = LIST_H(IS) - ENTH
        VM = LIST_H(IM) - ENTH
        IF( VS * VM  >   0.0_WP) THEN
            IS = IM
        ELSE
            IE = IM
        END IF
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BiSection_LIST_T(TEMP, IE, IS)
    !>    1) TABLE USED HERE IS DIMENSIONAL.
    !>    2) AS H IS SORTED FROM SMALL TO LARGE NUMBERS, AND
    !>          T IS A MONOTONE FUNCTION OF H, THUS
    !>          T IS ALSO FROM SMALL TO LARGE IN THE TABLE.

    USE thermal_info
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: TEMP
    INTEGER(4), INTENT(OUT) :: IE, IS
    REAL(WP) :: VS, VM
    INTEGER(4) :: IM

    !IF( ( TEMP < LIST_T(1) ).OR. ( TEMP  > LIST_T(N_LIST) ) ) &
    !CALL ERRHDL(' Reference TEMPERATURE IS OUT OF GIVEN NIST TABLE RANGE.', MYID)

    IS = 1
    IE = N_LIST
    DO WHILE ( (IE - IS) > 1 )
        IM = IS + (IE -IS) / 2
        VS = LIST_T(IS) - TEMP
        VM = LIST_T(IM) - TEMP
        IF( VS * VM  >   0.0_WP) THEN
            IS = IM
        ELSE
            IE = IM
        END IF
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CUBIC_SPLINE (N, X, Y, B, C, D)
!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE CoefficientS B, C, D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCREET FONCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X, Y    VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     B,C, D   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!     Function:
!     Y =  X
!     Reference:
!     FORSYTHE, G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE - HALL, INC.
!---------------------------------------------------------------------
    USE WPRECISION
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: N
    REAL(WP), INTENT(IN) :: X(N), Y(N)
    REAL(WP), INTENT(OUT) :: B(N), C(N), D(N)

    INTEGER(4) :: NM1, I, L
    REAL(WP) :: T


    IF (N < 2) RETURN
    IF (N < 3) THEN
        B(1) = (Y(2) - Y(1)) / (X(2) - X(1))
        C(1) = 0.0_WP
        D(1) = 0.0_WP
        B(2) = B(1)
        C(2) = 0.0_WP
        D(2) = 0.0_WP
        RETURN
    END IF

    ! Step 1: preparation
    !        BUILD THE TRIDIAGONAL SYSTEM
    !        B (DIAGONAL), D (UPPERDIAGONAL), C (SECOND MEMBER)
    NM1 = N - 1

    D(1) = X(2) - X(1)
    C(2) = (Y(2) - Y(1)) / D(1)
    DO I = 2, NM1
        D(I) = X(I + 1) - X(I)
        B(I) = 2.0_WP * (D(I - 1) + D(I))
        C(I + 1) = (Y(I + 1) - Y(I)) / D(I)
        C(I) = C(I + 1) - C(I)
    END DO

    ! step 2: end conditions
    !     CONDITIONS AT LIMITS
    !     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES
    B(1) = - D(1)
    B(N) = - D(N - 1)
    C(1) = 0.0_WP
    C(N) = 0.0_WP

    IF(N /= 3) THEN
        C(1) = C(3) / (X(4) - X(2)) - C(2) / (X(3) - X(1))
        C(N) = C(N - 1) / (X(N) - X(N - 2)) - C(N - 2) / (X(N - 1) - X(N - 3))
        C(1) = C(1) * D(1) * D(1) / (X(4) - X(1))
        C(N) = - C(N) * D(N - 1)**2 / (X(N) - X(N - 3))
    END IF

    ! step 3:     forward ELIMINATION
    DO I = 2, N
        T = D(I - 1) / B(I - 1)
        B(I) = B(I) - T * D(I - 1)
        C(I) = C(I) - T * C(I - 1)
    END DO

    !step 4:     BACK SUBSTITUTION
    C(N) = C(N) / B(N)
    DO  L = 1, NM1
        I = N - L
        C(I) = (C(I) - D(I) * C(I + 1)) / B(I)
    END DO

    !step 5: CoefficientS OF 3RD DEGREE POLYNOMIAL
    B(N) = (Y(N) - Y(NM1)) / D(NM1) + D(NM1) * (C(NM1) + 2.0_WP * C(N))
    DO  I = 1, NM1
        B(I) = (Y(I + 1) - Y(I)) / D(I) - D(I) * (C(I + 1) + 2.0_WP * C(I))
        D(I) = (C(I + 1) -C(I)) / D(I)
        C(I) = 3.0_WP * C(I)
    END DO
    C(N) = 3.0_WP * C(N)
    D(N) = D(NM1)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE FIND_LOCAL_MAXIMA(N, A, JMAX, NMAXL, AMAXL, JMIN, NMINL, AMINL)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N
    REAL(WP), INTENT (IN) :: A(N)
    INTEGER(4), INTENT(OUT) :: JMAX, JMIN           ! how many peaks
    !===== Max. 10 peaks max. and 10 peaks of mini.
    REAL(WP), INTENT (OUT) :: AMAXL(10), AMINL(10) ! peak values
    INTEGER(4), INTENT(OUT) :: NMAXL(10), NMINL(10) ! peaks location

    REAL(WP) :: MULTP, DL, DR
    INTEGER(4) :: I

    JMAX = 0
    JMIN = 0
    DO I = 2, n - 1
        DL = A(I) - A(I - 1)
        DR = A(I) - A(I + 1)
        MULTP = DL * DR
        IF(MULTP >  0.0_WP) THEN
            IF(DL  >   0.0_WP) THEN
                JMAX = JMAX + 1
                AMAXL(JMAX) = A(I)
                NMAXL(JMAX) = I
            END IF
            IF(DL  <  0.0_WP) THEN
                JMIN = JMIN + 1
                AMINL(JMIN) = A(I)
                NMINL(JMIN) = I
            END IF
        END IF
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE SPLINE_COEFFICIENTS_FOR_DH
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    INTEGER(4) :: I
    REAL(WP) :: H_tmp, T_tmp
    REAL(WP) :: DH_tmp
    REAL(WP) :: spline_interpolation_DHH
    REAL(WP) :: spline_interpolation_DHT
    REAL(WP) :: dDHdt1, dDHdt2, dDHdh1, dDHdh2

    IF(MYID /= 0) RETURN

    ALLOCATE(SplineCoeff_DHH_B(N_LIST)); SplineCoeff_DHH_B = 0.0_WP
    ALLOCATE(SplineCoeff_DHH_C(N_LIST)); SplineCoeff_DHH_C = 0.0_WP
    ALLOCATE(SplineCoeff_DHH_D(N_LIST)); SplineCoeff_DHH_D = 0.0_WP
    CALL CUBIC_SPLINE (N_LIST, LIST_DH, LIST_H, SplineCoeff_DHH_B, SplineCoeff_DHH_C, SplineCoeff_DHH_D)

    !/*below is used for equations only*/
    ALLOCATE(SplineCoeff_DHT_B(N_LIST)); SplineCoeff_DHT_B = 0.0_WP
    ALLOCATE(SplineCoeff_DHT_C(N_LIST)); SplineCoeff_DHT_C = 0.0_WP
    ALLOCATE(SplineCoeff_DHT_D(N_LIST)); SplineCoeff_DHT_D = 0.0_WP
    CALL CUBIC_SPLINE (N_LIST, LIST_DH, LIST_T, SplineCoeff_DHT_B, SplineCoeff_DHT_C, SplineCoeff_DHT_D)

    OPEN(20, FILE = TRIM(FilePath0) // 'CHK_LIST_TABLE_NONDIM_SPLINE_DH.dat')
    WRITE(20, '(A)') '# DH  H  T'
    WRITE(20, *) '# ', N_LIST * 3 + 1

    DO I = 1, N_LIST * 3 + 1
        DH_tmp = LIST_DH(IMIN_DH) + DBLE(I - 1) * (LIST_DH(IMAX_DH) - LIST_DH(IMIN_DH)) / DBLE(N_LIST * 3)
        H_tmp = spline_interpolation_DHH(DH_tmp)
        T_tmp = spline_interpolation_DHT(DH_tmp)
        WRITE(20, '(3ES13.5)') DH_tmp, H_tmp, T_tmp
        ! check monotonicity
    END DO
    CLOSE(20)

    ! check monotonicity of rH =f(T) and rH =f(h)
    DO I = 2, N_LIST - 1
        dDHdt1 = ( LIST_DH(I) - LIST_DH(I - 1) ) /  ( LIST_T(I) - LIST_T(I - 1) )
        dDHdt2 = ( LIST_DH(I + 1) - LIST_DH(I) ) /  ( LIST_T(I + 1) - LIST_T(I) )
        IF ( dDHdt1 * dDHdt2 < 0.0_WP) CALL ERRHDL(' Error. The relation (rho * h) = FUNCTION (T) is not monotonicity.')
        dDHdh1 = ( LIST_DH(I) - LIST_DH(I - 1) ) /  ( LIST_H(I) - LIST_H(I - 1) )
        dDHdh2 = ( LIST_DH(I + 1) - LIST_DH(I) ) /  ( LIST_H(I + 1) - LIST_H(I) )
        IF ( dDHdh1 * dDHdh2 < 0.0_WP) CALL ERRHDL(' Error. The relation (rho * h) = FUNCTION (h) is not monotonicity.')
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE SPLINE_COEFFICIENTS_FOR_H
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    INTEGER(4) :: I
    REAL(WP) :: H_tmp
    REAL(WP) :: T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp, dDDH_tmp
    REAL(WP) :: D0TMP, H0TMP, dDdH0TMP
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: spline_interpolation_HB

    IF(MYID /= 0) RETURN

    ALLOCATE(SplineCoeff_HT_B(N_LIST)); SplineCoeff_HT_B = 0.0_WP
    ALLOCATE(SplineCoeff_HT_C(N_LIST)); SplineCoeff_HT_C = 0.0_WP
    ALLOCATE(SplineCoeff_HT_D(N_LIST)); SplineCoeff_HT_D = 0.0_WP

    ALLOCATE(SplineCoeff_HD_B(N_LIST)); SplineCoeff_HD_B = 0.0_WP
    ALLOCATE(SplineCoeff_HD_C(N_LIST)); SplineCoeff_HD_C = 0.0_WP
    ALLOCATE(SplineCoeff_HD_D(N_LIST)); SplineCoeff_HD_D = 0.0_WP

    ALLOCATE(SplineCoeff_HM_B(N_LIST)); SplineCoeff_HM_B = 0.0_WP
    ALLOCATE(SplineCoeff_HM_C(N_LIST)); SplineCoeff_HM_C = 0.0_WP
    ALLOCATE(SplineCoeff_HM_D(N_LIST)); SplineCoeff_HM_D = 0.0_WP

    ALLOCATE(SplineCoeff_HK_B(N_LIST)); SplineCoeff_HK_B = 0.0_WP
    ALLOCATE(SplineCoeff_HK_C(N_LIST)); SplineCoeff_HK_C = 0.0_WP
    ALLOCATE(SplineCoeff_HK_D(N_LIST)); SplineCoeff_HK_D = 0.0_WP

    ALLOCATE(SplineCoeff_HCp_B(N_LIST)); SplineCoeff_HCp_B = 0.0_WP
    ALLOCATE(SplineCoeff_HCp_C(N_LIST)); SplineCoeff_HCp_C = 0.0_WP
    ALLOCATE(SplineCoeff_HCp_D(N_LIST)); SplineCoeff_HCp_D = 0.0_WP

    ALLOCATE(SplineCoeff_HB_B(N_LIST)); SplineCoeff_HB_B = 0.0_WP
    ALLOCATE(SplineCoeff_HB_C(N_LIST)); SplineCoeff_HB_C = 0.0_WP
    ALLOCATE(SplineCoeff_HB_D(N_LIST)); SplineCoeff_HB_D = 0.0_WP

    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_T, SplineCoeff_HT_B, SplineCoeff_HT_C, SplineCoeff_HT_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_D, SplineCoeff_HD_B, SplineCoeff_HD_C, SplineCoeff_HD_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_M, SplineCoeff_HM_B, SplineCoeff_HM_C, SplineCoeff_HM_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_K, SplineCoeff_HK_B, SplineCoeff_HK_C, SplineCoeff_HK_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_CP, SplineCoeff_HCp_B, SplineCoeff_HCp_C, SplineCoeff_HCp_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_H, LIST_B, SplineCoeff_HB_B, SplineCoeff_HB_C, SplineCoeff_HB_D)

    IF(MYID == 0) THEN
        OPEN(10, FILE = TRIM(FilePath0) // 'CHK_LIST_TABLE_NONDIM_SPLINE_H.dat')
        WRITE(10, '(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA  dDdH'
        WRITE(10, *) '# ', N_LIST* 3 + 1

        DHmax = -1E-14_WP
        DHmin = + 1E+14_WP
        CpMaX = -1E-14_WP
        DO I = 2, N_LIST*3

            !D0TMP = D0
            !H0TMP = H0

            H_tmp = LIST_H(1) + DBLE(I - 1) * (LIST_H(N_LIST) - LIST_H(1)) / DBLE(N_LIST*3)
            T_tmp = spline_interpolation_HT(H_tmp)
            D_tmp = spline_interpolation_HD(H_tmp)
            M_tmp = spline_interpolation_HM(H_tmp)
            K_tmp = spline_interpolation_HK(H_tmp)
            Cp_tmp = spline_interpolation_HCp(H_tmp)
            B_tmp = spline_interpolation_HB(H_tmp)

            !dDdH0TMP = (D_tmp - D0TMP) / (H_tmp - H0TMP)

            WRITE(10, '(8ES13.5)') P0, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp!, dDdH0TMP
        END DO

        DHmax = LIST_DH(IMAX_DH)
        DHmin = LIST_DH(IMIN_DH)
        T4MaxDH = LIST_T(IMAX_DH)
        T4MinDH = LIST_T(IMIN_DH)
        WRITE(10, *) '# Max DH =               ', DHmax
        WRITE(10, *) '# Corresponding Tempeture = ', T4MaxDH
        WRITE(10, *) '# Min DH =               ', DHmin
        WRITE(10, *) '# Corresponding Tempeture = ', T4MinDH
        CLOSE(10)

        CALL CHKRLHDL  ('Max DH (undim) = ', MYID, DHmax)
        CALL CHKRLHDL  (' At Temperature(undim) = ', MYID, T4MaxDH)
        CALL CHKRLHDL  (' At Temperature(dim)   = ', MYID, T4MaxDH* T0)
        CALL CHKRLHDL  ('Min DH (undim): = ', MYID, DHmin)
        CALL CHKRLHDL  (' At Temperature(undim) = ', MYID, T4MinDH)
        CALL CHKRLHDL  (' At Temperature(dim)   = ', MYID, T4MinDH* T0)

        CALL CHKRLHDL  ('Max Cp (undim) = ', MYID, CpMax)
        CALL CHKRLHDL  (' Tpc = At Temperature(undim) = ', MYID, T4CpMax)
        CALL CHKRLHDL  (' Tpc = At Temperature(dim)   = ', MYID, T4CpMaX * T0)



    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE SPLINE_COEFFICIENTS_FOR_T
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    INTEGER(4) :: I
    REAL(WP) :: H_tmp
    REAL(WP) :: T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp
    REAL(WP) :: spline_interpolation_TH
    REAL(WP) :: spline_interpolation_TD
    REAL(WP) :: spline_interpolation_TM
    REAL(WP) :: spline_interpolation_TK
    REAL(WP) :: spline_interpolation_TCp
    REAL(WP) :: spline_interpolation_TB

    IF(MYID /= 0) RETURN


    ALLOCATE(SplineCoeff_TH_B(N_LIST)); SplineCoeff_TH_B = 0.0_WP
    ALLOCATE(SplineCoeff_TH_C(N_LIST)); SplineCoeff_TH_C = 0.0_WP
    ALLOCATE(SplineCoeff_TH_D(N_LIST)); SplineCoeff_TH_D = 0.0_WP

    ALLOCATE(SplineCoeff_TD_B(N_LIST)); SplineCoeff_TD_B = 0.0_WP
    ALLOCATE(SplineCoeff_TD_C(N_LIST)); SplineCoeff_TD_C = 0.0_WP
    ALLOCATE(SplineCoeff_TD_D(N_LIST)); SplineCoeff_TD_D = 0.0_WP

    ALLOCATE(SplineCoeff_TM_B(N_LIST)); SplineCoeff_TM_B = 0.0_WP
    ALLOCATE(SplineCoeff_TM_C(N_LIST)); SplineCoeff_TM_C = 0.0_WP
    ALLOCATE(SplineCoeff_TM_D(N_LIST)); SplineCoeff_TM_D = 0.0_WP

    ALLOCATE(SplineCoeff_TK_B(N_LIST)); SplineCoeff_TK_B = 0.0_WP
    ALLOCATE(SplineCoeff_TK_C(N_LIST)); SplineCoeff_TK_C = 0.0_WP
    ALLOCATE(SplineCoeff_TK_D(N_LIST)); SplineCoeff_TK_D = 0.0_WP

    ALLOCATE(SplineCoeff_TCp_B(N_LIST)); SplineCoeff_TCp_B = 0.0_WP
    ALLOCATE(SplineCoeff_TCp_C(N_LIST)); SplineCoeff_TCp_C = 0.0_WP
    ALLOCATE(SplineCoeff_TCp_D(N_LIST)); SplineCoeff_TCp_D = 0.0_WP

    ALLOCATE(SplineCoeff_TB_B(N_LIST)); SplineCoeff_TB_B = 0.0_WP
    ALLOCATE(SplineCoeff_TB_C(N_LIST)); SplineCoeff_TB_C = 0.0_WP
    ALLOCATE(SplineCoeff_TB_D(N_LIST)); SplineCoeff_TB_D = 0.0_WP

    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_H, SplineCoeff_TH_B, SplineCoeff_TH_C, SplineCoeff_TH_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_D, SplineCoeff_TD_B, SplineCoeff_TD_C, SplineCoeff_TD_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_M, SplineCoeff_TM_B, SplineCoeff_TM_C, SplineCoeff_TM_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_K, SplineCoeff_TK_B, SplineCoeff_TK_C, SplineCoeff_TK_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_CP, SplineCoeff_TCp_B, SplineCoeff_TCp_C, SplineCoeff_TCp_D)
    CALL CUBIC_SPLINE (N_LIST, LIST_T, LIST_B, SplineCoeff_TB_B, SplineCoeff_TB_C, SplineCoeff_TB_D)


    IF(MYID == 0) THEN
        OPEN(10, FILE = TRIM(FilePath0) // 'CHK_LIST_TABLE_NONDIM_SPLINE_T.dat')
        WRITE(10, '(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA'
        WRITE(10, *) '# ', N_LIST*3

        DO I = 2, N_LIST*3
            T_tmp = LIST_T(1) + DBLE(I - 1) * (LIST_T(N_LIST) - LIST_T(1)) / DBLE(N_LIST*3)
            H_tmp = spline_interpolation_TH(T_tmp)
            D_tmp = spline_interpolation_TD(T_tmp)
            M_tmp = spline_interpolation_TM(T_tmp)
            K_tmp = spline_interpolation_TK(T_tmp)
            Cp_tmp = spline_interpolation_TCp(T_tmp)
            B_tmp = spline_interpolation_TB(T_tmp)

            WRITE(10, '(8ES13.5)') P0, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp
        END DO
        CLOSE(10)

    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE Property_Table_Nondimensionalization
!>    THE FORMAT OF NIST TABLE USED HERE SHUOLD BE
!>            # P H  T  \RHO \MU \KAPPA CP
!>            NUMBER OF LINES
!>            P  H  T  D  M  K CP
!>            ................
!>
!>    ALL H VALUSE are SORTED FROM SMALL TO BIG NUMBERS.

!>    As H has been sorted, the Newton BISection method can be used for searching.
!>    IF the H data has not been sorted, binary tree method IS recommended.

        USE thermal_info
        USE mesh_info
        USE flow_info
        USE INIT_INFO
        IMPLICIT NONE

        INTEGER(4) :: IOS = 0
        INTEGER(4) :: Fflg = 15           !file I /O id
        CHARACTER(64) :: STMP
        REAL(WP) :: RTMP

        REAL(WP) :: SFE, SFS
        INTEGER(4) :: IE, IS, I
        INTEGER(4) :: K
        REAL(WP) :: BUF_T, BUF_H, BUF_D, BUF_M, BUF_K, BUF_B, BUF_Cp

        !INTEGER(4):: JMAX = 0, JMIN = 0           ! how many peaks
        !===== Max. 10 peaks max. and 10 peaks of mini.
        !REAL(WP) :: AMAXL(10) = 0.0_WP, AMINL(10) = 0.0_WP ! peak values
        !INTEGER(4):: NMAXL(10) = 0, NMINL(10) = 0 ! peaks location


        IF(MYID /= 0) RETURN

        !================ READ IN NIST TABLE (DIMENSIONAL) =====================
        OPEN(FFLG, FILE = TRIM(NISTFLNM), STATUS = 'old', IOSTAT = IOS)
        IF(IOS  /= 0)  &
        CALL ERRHDL(' File ' // TRIM(NISTFLNM) // ' cannot be found.', MYID)

        READ(FFLG, *) STMP
        READ(FFLG, *) N_LIST

        ALLOCATE( LIST_H (N_LIST) ) ;  LIST_H  = 0.0_WP
        ALLOCATE( LIST_T (N_LIST) ) ;  LIST_T  = 0.0_WP
        ALLOCATE( LIST_D (N_LIST) ) ;  LIST_D  = 0.0_WP
        ALLOCATE( LIST_M (N_LIST) ) ;  LIST_M  = 0.0_WP
        ALLOCATE( LIST_K (N_LIST) ) ;  LIST_K  = 0.0_WP
        ALLOCATE( LIST_B (N_LIST) ) ;  LIST_B  = 0.0_WP
        ALLOCATE( LIST_CP(N_LIST) ) ;  LIST_CP = 0.0_WP
        ALLOCATE( LIST_DH(N_LIST) ) ;  LIST_DH = 0.0_WP!updated by Junjie, 2017/03 /13
        MEMPC_Byte = MEMPC_Byte + N_LIST* 8 *8

        DO I = 1, N_LIST
            READ(FFLG, *) RTMP, &
                        LIST_H(I), LIST_T(I), LIST_D(I), &
                        LIST_M(I), LIST_K(I), LIST_CP(I), LIST_B(I)
        END DO

        ! Data sorting from small to big, based on LIST_T, idiot-proof
        DO I = 1, N_LIST
            K = MINLOC( LIST_T(I : N_LIST), DIM = 1 ) + I - 1

            BUF_T = LIST_T(I)
            BUF_D = LIST_D(I)
            BUF_H = LIST_H(I)
            BUF_M = LIST_M(I)
            BUF_K = LIST_K(I)
            BUF_B = LIST_B(I)
            BUF_Cp = LIST_Cp(I)

            LIST_T(I) = LIST_T(K)
            LIST_D(I) = LIST_D(K)
            LIST_H(I) = LIST_H(K)
            LIST_M(I) = LIST_M(K)
            LIST_K(I) = LIST_K(K)
            LIST_B(I) = LIST_B(K)
            LIST_Cp(I) = LIST_Cp(K)

            LIST_T(K) = BUF_T
            LIST_D(K) = BUF_D
            LIST_H(K) = BUF_H
            LIST_M(K) = BUF_M
            LIST_K(K) = BUF_K
            LIST_B(K) = BUF_B
            LIST_Cp(K) = BUF_Cp
        END DO


        !============== Reference statment ==== (dimensional) ======================
        CALL BiSection_LIST_T(T0, IE, IS)
        SFS = ( LIST_T(IE) - T0 ) / ( LIST_T(IE) - LIST_T(IS) )
        SFE = 1.0_WP - SFS

        H0  = SFE * LIST_H(IE)  + SFS * LIST_H(IS)
        D0  = SFE * LIST_D(IE)  + SFS * LIST_D(IS)
        K0  = SFE * LIST_K(IE)  + SFS * LIST_K(IS)
        M0  = SFE * LIST_M(IE)  + SFS * LIST_M(IS)
        CP0 = SFE * LIST_CP(IE) + SFS * LIST_CP(IS) ! updated
        B0  = SFE * LIST_B(IE)  + SFS * LIST_B(IS)

        !================ Local peakS =================================== LIST_SPLINE ========
        ! CALL FIND_LOCAL_MAXIMA(N_LIST, LIST_K, JMAX, NMAXL, AMAXL, JMIN, NMINL, AMINL)
        ! DO I = 1, JMAX
        !     CALL CHK2RLHDL  ('    Local Max. of K = ', MYID, AMAXL(I), AMAXL(I) / K0 )
        !     CALL CHK2RLHDL  ('          Locate at T = ', MYID, LIST_T(NMAXL(I)), LIST_T(NMAXL(I)) / T0)
        ! END DO
        ! DO I = 1, JMIN
        !     CALL CHK2RLHDL  ('    Local Min. of K = ', MYID, AMINL(I), AMINL(I) / K0 )
        !     CALL CHK2RLHDL  ('          Locate at T = ', MYID, LIST_T(NMINL(I)), LIST_T(NMINL(I)) / T0)
        ! END DO

        !=============== Scale the NIST table based on the Reference statE == (dimensionlessization) ===============
        DO I = 1, N_LIST
            LIST_H(I) = (LIST_H(I) - H0) / T0 / CP0
            LIST_T(I) = LIST_T(I) / T0
            LIST_D(I) = LIST_D(I) / D0
            LIST_M(I) = LIST_M(I) / M0
            LIST_K(I) = LIST_K(I) / K0
            LIST_CP(I) = LIST_CP(I) / CP0
            LIST_B(I) = LIST_B(I) / B0
            LIST_DH(I) = LIST_D(I) * LIST_H(I)!updated by Junjie, 2017/03 /13
        END DO

        IMAX_H = MAXLOC(LIST_H, DIM = 1)
        IMIN_H = MINLOC(LIST_H, DIM = 1)

        IMAX_T = MAXLOC(LIST_H, DIM = 1)
        IMIN_T = MINLOC(LIST_H, DIM = 1)

        IMAX_DH = MAXLOC(LIST_DH, DIM = 1)
        IMIN_DH = MINLOC(LIST_DH, DIM = 1)

        IMAX_Cp = MAXLOC(LIST_H, DIM = 1)
        IMIN_Cp = MINLOC(LIST_H, DIM = 1)



        RETURN
END SUBROUTINE
