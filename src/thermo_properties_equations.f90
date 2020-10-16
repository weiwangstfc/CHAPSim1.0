!**********************************************************************************************************************************
!> @brief
!>        functions to calculate thermo-properties from equations.
!> @details
!> SUBROUTINE: thermophysical_function_Check_Trange (in MYID = all)
!> function: thermophysical_function_TD (in MYID = all)
!> function: thermophysical_function_TCp (in MYID = all)
!> function: thermophysical_function_TM
!> function: thermophysical_function_TK
!> function: thermophysical_function_TB
!> function: thermophysical_function_TH
!> SUBROUTINE:thermophysical_function_DH_preparation
!> @note
!> @toDO
! REVISION HISTORY:
! 09/2020 - Created, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE thermophysical_function_Check_Trange(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: T_tmp0
    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: T_tmp
    ! dim = 0 :: dimensionless
    ! dim = 1 :: dimensional
    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), D(kg/ M3)
        CASE (iLiquidSodium)
            IF (T_tmp < Tm0_Na .OR. T_tmp > Tb0_Na) THEN
                CALL CHK3RLHDL("Tm0_Na, Tb0_Na, T_tmp (K) = ", MYID, Tm0_Na, Tb0_Na, T_tmp)
                CALL ERRHDL("The calcualted temperature exceeds its limits in a liquid state.", MYID)
            END IF
        CASE (iLiquidLead)
            IF (T_tmp < Tm0_Pb .OR. T_tmp > Tb0_Pb) THEN
                CALL CHK3RLHDL("Tm0_Pb, Tb0_Pb, T_tmp (K) = ", MYID, Tm0_Pb, Tb0_Pb, T_tmp)
                CALL ERRHDL("The calcualted temperature exceeds its limits in a liquid state.", MYID)
            END IF
        CASE (iLiquidBismuth)
            IF (T_tmp < Tm0_Bi .OR. T_tmp > Tb0_Bi) THEN
                CALL CHK3RLHDL("Tm0_Bi, Tb0_Bi, T_tmp (K) = ", MYID, Tm0_Bi, Tb0_Bi, T_tmp)
                CALL ERRHDL("The calcualted temperature exceeds its limits in a liquid state.", MYID)
            END IF
        CASE (iLiquidLBE)
            IF (T_tmp < Tm0_LBE .OR. T_tmp > Tb0_LBE) THEN
                CALL CHK3RLHDL("Tm0_LBE, Tb0_LBE, T_tmp (K) = ", MYID, Tm0_LBE, Tb0_LBE, T_tmp)
                CALL ERRHDL("The calcualted temperature exceeds its limits in a liquid state.", MYID)
            END IF
        CASE DEFAULT
            IF (T_tmp < Tm0_Na .OR. T_tmp > Tb0_Na) THEN
                CALL CHK3RLHDL("Tm0_Na, Tb0_Na, T_tmp (K) = ", MYID, Tm0_Na, Tb0_Na, T_tmp)
                CALL ERRHDL("The calcualted temperature exceeds its limits in a liquid state.", MYID)
            END IF
    END SELECT

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TD(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: thermophysical_function_TD
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: D_tmp, T_tmp

    ! dim = 0 :: dimensionless
    ! dim = 1 :: dimensional
    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), D(kg/ M3)
        CASE (iLiquidSodium)
            D_tmp = CoD_Na(0) + CoD_Na(1) * T_tmp
        CASE (iLiquidLead)
            D_tmp = CoD_Pb(0) + CoD_Pb(1) * T_tmp
        CASE (iLiquidBismuth)
            D_tmp = CoD_Bi(0) + CoD_Bi(1) * T_tmp
        CASE (iLiquidLBE)
            D_tmp = CoD_LBE(0) + CoD_LBE(1) * T_tmp
        CASE DEFAULT
            D_tmp = CoD_Na(0) + CoD_Na(1) * T_tmp
    END SELECT

    IF (dim == 0) D_tmp = D_tmp / D0 ! convert dim to umdim

    thermophysical_function_TD = D_tmp

END FUNCTION

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TCp(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: thermophysical_function_TCp
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: Cp_tmp, T_tmp

    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), Cp(J / Kg/ K)
        CASE (iLiquidSodium)
            Cp_tmp = CoCp_Na(-2) * T_tmp**(-2) + CoCp_Na(-1) * T_tmp**(-1) + CoCp_Na(0) + CoCp_Na(1) * T_tmp + CoCp_Na(2) * T_tmp**2
        CASE (iLiquidLead)
            Cp_tmp = CoCp_Pb(-2) * T_tmp**(-2) + CoCp_Pb(-1) * T_tmp**(-1) + CoCp_Pb(0) + CoCp_Pb(1) * T_tmp + CoCp_Pb(2) * T_tmp**2
        CASE ( iLiquidBismuth)
            Cp_tmp = CoCp_Bi(-2) * T_tmp**(-2) + CoCp_Bi(-1) * T_tmp**(-1) + CoCp_Bi(0) + CoCp_Bi(1) * T_tmp + CoCp_Bi(2) * T_tmp**2
        CASE (iLiquidLBE)
            Cp_tmp = CoCp_LBE(-2) * T_tmp**(-2) + CoCp_LBE(-1) * T_tmp**(-1) + CoCp_LBE(0) + &
                     CoCp_LBE(1) * T_tmp + CoCp_LBE(2) * T_tmp**2
        CASE DEFAULT
            Cp_tmp = CoCp_Na(-2) * T_tmp**(-2) + CoCp_Na(-1) * T_tmp**(-1) + CoCp_Na(0) + CoCp_Na(1) * T_tmp + CoCp_Na(2) * T_tmp**2
    END SELECT

    IF (dim == 0) Cp_tmp = Cp_tmp / Cp0 ! convert dim to umdim

    thermophysical_function_TCp = Cp_tmp

END FUNCTION

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TM(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    REAL(WP) :: thermophysical_function_TM
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: M_tmp, T_tmp

    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), M(Pa S)
        CASE (iLiquidSodium)
            M_tmp = EXP( CoM_Na(-1) / T_tmp + CoM_Na(0) + CoM_Na(1) * LOG(T_tmp) )
        CASE (iLiquidLead)
            M_tmp = CoM_Pb(0) * EXP (CoM_Pb(-1) / T_tmp)
        CASE (iLiquidBismuth)
            M_tmp = CoM_Bi(0) * EXP (CoM_Bi(-1) / T_tmp)
        CASE (iLiquidLBE)
            M_tmp = CoM_LBE(0) * EXP (CoM_LBE(-1) / T_tmp)
        CASE DEFAULT
            M_tmp = EXP( CoM_Na(-1) / T_tmp + CoM_Na(0) + CoM_Na(1) * LOG(T_tmp) )
    END SELECT

    IF (dim == 0) M_tmp = M_tmp / M0 ! convert dim to umdim
    thermophysical_function_TM = M_tmp

END FUNCTION

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TK(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: thermophysical_function_TK
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: K_tmp, T_tmp

    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), K(W/ M / K)
        CASE (iLiquidSodium)
            K_tmp = CoK_Na(0) + CoK_Na(1) * T_tmp + CoK_Na(2) * T_tmp**2
        CASE (iLiquidLead)
            K_tmp = CoK_Pb(0) + CoK_Pb(1) * T_tmp + CoK_Pb(2) * T_tmp**2
        CASE (iLiquidBismuth)
            K_tmp = CoK_Bi(0) + CoK_Bi(1) * T_tmp + CoK_Bi(2) * T_tmp**2
        CASE (iLiquidLBE)
            K_tmp = CoK_LBE(0) + CoK_LBE(1) * T_tmp + CoK_LBE(2) * T_tmp**2
        CASE DEFAULT
            K_tmp = CoK_Na(0) + CoK_Na(1) * T_tmp + CoK_Na(2) * T_tmp**2
    END SELECT

    IF (dim == 0) K_tmp = K_tmp / K0 ! convert dim to umdim
    thermophysical_function_TK = K_tmp

END FUNCTION

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TB(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: thermophysical_function_TB
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: B_tmp, T_tmp

    ! dim = 0 :: dimensionless
    ! dim = 1 :: dimensional
    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), B(1 / K)
        CASE (iLiquidSodium)
            B_tmp = 1.0_WP / (CoB_Na - T_tmp)
        CASE (iLiquidLead)
            B_tmp = 1.0_WP / (CoB_Pb - T_tmp)
        CASE (iLiquidBismuth)
            B_tmp = 1.0_WP / (CoB_Bi - T_tmp)
        CASE (iLiquidLBE)
            B_tmp = 1.0_WP / (CoB_LBE - T_tmp)
        CASE DEFAULT
            B_tmp = 1.0_WP / (CoB_Na - T_tmp)  ! default IS Na.
    END SELECT

    IF (dim == 0) B_tmp = B_tmp / B0 ! convert dim to umdim

    thermophysical_function_TB = B_tmp

END FUNCTION

!**********************************************************************************************************************************
FUNCTION thermophysical_function_TH(T_tmp0, dim)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: thermophysical_function_TH
    REAL(WP), INTENT(IN) :: T_tmp0

    INTEGER(4), INTENT(IN) :: dim
    REAL(WP) :: H_tmp, T_tmp

    ! dim = 0 :: dimensionless
    ! dim = 1 :: dimensional
    IF (dim == 0) THEN
        T_tmp = T_tmp0 * T0 ! convert undim to dim
    ELSE
        T_tmp = T_tmp0
    END IF

    SELECT CASE (iFluidMedia)
        ! unit: T(Kelvin), H(J / Kg)
        CASE (iLiquidSodium)
            H_tmp = Hm0_Na + &
                    CoH_Na(-1) * (1.0_WP / T_tmp - 1.0_WP / Tm0_Na) + &
                    CoH_Na(0) + &
                    CoH_Na(1) * (T_tmp - Tm0_Na) + &
                    CoH_Na(2) * (T_tmp**2 - Tm0_Na**2) + &
                    CoH_Na(3) * (T_tmp**3 - Tm0_Na**3)
        CASE (iLiquidLead)
            H_tmp = Hm0_Pb + &
                    CoH_Pb(-1) * (1.0_WP / T_tmp - 1.0_WP / Tm0_Pb) + &
                    CoH_Pb(0) + &
                    CoH_Pb(1) * (T_tmp - Tm0_Pb) + &
                    CoH_Pb(2) * (T_tmp**2 - Tm0_Pb**2) + &
                    CoH_Pb(3) * (T_tmp**3 - Tm0_Pb**3)
        CASE (iLiquidBismuth)
            H_tmp = Hm0_Bi + &
                    CoH_Bi(-1) * (1.0_WP / T_tmp - 1.0_WP / Tm0_Bi) + &
                    CoH_Bi(0) + &
                    CoH_Bi(1) * (T_tmp - Tm0_Bi) + &
                    CoH_Bi(2) * (T_tmp**2 - Tm0_BI**2) + &
                    CoH_Bi(3) * (T_tmp**3 - Tm0_BI**3)
        CASE (iLiquidLBE)
            H_tmp = Hm0_LBE + &
                    CoH_LBE(-1) * (1.0_WP / T_tmp - 1.0_WP / Tm0_LBE) + &
                    CoH_LBE(0) + &
                    CoH_LBE(1) * (T_tmp - Tm0_LBE) + &
                    CoH_LBE(2) * (T_tmp**2 - Tm0_LBE**2) + &
                    CoH_LBE(3) * (T_tmp**3 - Tm0_LBE**3)
        CASE DEFAULT
            H_tmp = Hm0_Na + &
                    CoH_Na(-1) * (1.0_WP / T_tmp - 1.0_WP / Tm0_Na) + &
                    CoH_Na(0) + &
                    CoH_Na(1) * (T_tmp - Tm0_Na) + &
                    CoH_Na(2) * (T_tmp**2 - Tm0_Na**2) + &
                    CoH_Na(3) * (T_tmp**3 - Tm0_Na**3)
    END SELECT

    IF (dim == 0) H_tmp = (H_tmp - H0) / (Cp0 * T0) ! convert dim to umdim

    thermophysical_function_TH = H_tmp

END FUNCTION

!**********************************************************************************************************************************
SUBROUTINE thermophysical_function_DH_preparation
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    INTEGER(4) :: dim
    REAL(WP) :: T, D, H
    INTEGER(4) :: I, K
    REAL(WP) :: BUF_T, BUF_H, BUF_DH
    REAL(WP) :: thermophysical_function_TD
    REAL(WP) :: thermophysical_function_TH

    N_LIST = 1024

    ALLOCATE( LIST_T (N_LIST) ) ;  LIST_T = 0.0_WP
    ALLOCATE( LIST_H (N_LIST) ) ;  LIST_T = 0.0_WP
    ALLOCATE( LIST_DH(N_LIST) ) ;  LIST_DH = 0.0_WP

    ! for dimensionless
    dim = 0
    SELECT CASE (iFluidMedia)
    CASE (iLiquidSodium)
        DO I = 1, N_LIST
            T = ( Tm0_Na + (Tb0_Na - Tm0_Na) / DBLE(N_LIST) * DBLE(I) ) / T0
            D = thermophysical_function_TD(T, dim)
            H = thermophysical_function_TH(T, dim)
            LIST_T(I) = T
            LIST_H(I) = H
            LIST_DH(I) = D * H
        END DO
    CASE (iLiquidLead)
        DO I = 1, N_LIST
            T = ( Tm0_Pb + (Tb0_Pb - Tm0_Pb) / DBLE(N_LIST) * DBLE(I) ) / T0
            D = thermophysical_function_TD(T, dim)
            H = thermophysical_function_TH(T, dim)
            LIST_T(I) = T
            LIST_H(I) = H
            LIST_DH(I) = D * H
        END DO
    CASE (iLiquidBismuth)
        DO I = 1, N_LIST
            T = ( Tm0_Bi + (Tb0_Bi - Tm0_Bi) / DBLE(N_LIST) * DBLE(I) ) / T0
            D = thermophysical_function_TD(T, dim)
            H = thermophysical_function_TH(T, dim)
            LIST_T(I) = T
            LIST_H(I) = H
            LIST_DH(I) = D * H
        END DO
    CASE (iLiquidLBE)
        DO I = 1, N_LIST
            T = ( Tm0_LBE + (Tb0_LBE - Tm0_LBE) / DBLE(N_LIST) * DBLE(I) ) / T0
            D = thermophysical_function_TD(T, dim)
            H = thermophysical_function_TH(T, dim)
            LIST_T(I) = T
            LIST_H(I) = H
            LIST_DH(I) = D * H
        END DO
    CASE DEFAULT ! Na
        DO I = 1, N_LIST
            T = ( Tm0_Na + (Tb0_Na - Tm0_Na) / DBLE(N_LIST) * DBLE(I) ) / T0
            D = thermophysical_function_TD(T, dim)
            H = thermophysical_function_TH(T, dim)
            LIST_T(I) = T
            LIST_H(I) = H
            LIST_DH(I) = D * H
        END DO
    END SELECT

    ! Data sorting from small to big, based on LIST_T, idiot-proof
    DO I = 1, N_LIST
        K = MINLOC( LIST_T(I : N_LIST), DIM = 1 ) + I - 1
        BUF_T = LIST_T(I)
        BUF_H = LIST_H(I)
        BUF_DH = LIST_DH(I)

        LIST_T(I) = LIST_T(K)
        LIST_H(I) = LIST_H(K)
        LIST_DH(I) = LIST_DH(K)

        LIST_T(K) = BUF_T
        LIST_H(K) = BUF_H
        LIST_DH(K) = BUF_DH
    END DO


    IMAX_T = MAXLOC(LIST_T, DIM = 1)
    IMIN_T = MINLOC(LIST_T, DIM = 1)

    IMAX_H = MAXLOC(LIST_H, DIM = 1)
    IMIN_H = MINLOC(LIST_H, DIM = 1)

    IMAX_DH = MAXLOC(LIST_DH, DIM = 1)
    IMIN_DH = MINLOC(LIST_DH, DIM = 1)

    CALL SPLINE_COEFFICIENTS_FOR_DH

    RETURN

END SUBROUTINE
