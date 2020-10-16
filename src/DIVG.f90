!**********************************************************************************************************************************
!> @brief
!>       divergence of the velocity field.
!> @details
!> SUBROUTINE: DIVG_tg (in MYID = all)
!> SUBROUTINE: DIVG_io(NS) (in MYID = all)
!> SUBROUTINE: DIVG_U_io(in MYID = all)
!> SUBROUTINE: DENSITY_TIME_DERIVATION (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 12 / 2013- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE DIVG_tg(NS)
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS

    INTEGER(4) :: KC, KP
    INTEGER(4) :: JC, JP, JJ
    INTEGER(4) :: IC, IP
    REAL(WP) :: DVIGVELO
    REAL(WP) :: COE0, COE1, COE2

    COE0 = 1.0_WP / (DT * TALP(NS))
    RHSLLPHI_tg = 0.0_WP

    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        COE1 = DXI / RCCI2(JJ)
        COE2 = DYFI(JJ) / RCCI1(JJ)
        DO KC = 1, NCL3
            KP = KPV(KC)
            DO IC = 1, NCL1_tg
                IP = IPV_tg(IC)
                DVIGVELO= (Q_tg(IP, JC, KC, 1) - Q_tg(IC, JC, KC, 1)) * COE1 +  &
                          (Q_tg(IC, JP, KC, 2) - Q_tg(IC, JC, KC, 2)) * COE2 +  &
                          (Q_tg(IC, JC, KP, 3) - Q_tg(IC, JC, KC, 3)) * DZI
                RHSLLPHI_tg(IC, JC, KC) = DVIGVELO * COE0
            END DO
        END DO
    END DO

    RETURN

END SUBROUTINE

!*******************************************************************************************
SUBROUTINE DIVG_io(NS)
    !>    @NOTE
    !>    1) IF the velocity field IS calculated based on timE -fixed DENSITY field,
    !>       the DENSITY time rate IS excluded from the continuity equation.
    USE thermal_info
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS

    INTEGER(4) :: IC, JC, KC
    INTEGER(4) :: IP, JP, KP, JJ
    REAL(WP) :: DIVX
    REAL(WP) :: DIVY
    REAL(WP) :: DIVZ
    REAL(WP) :: DTRHO, DQCAP
    REAL(WP) :: COE0, COE1, COE2, COE4

    RHSLLPHI_io = 0.0_WP
    DTRHO = 0.0_WP
    DQCAP = 0.0_WP

    IF(iWeightedPre == 1) THEN
        COE0 = 1.0_WP / DT / TALP(NS) / (0.25_WP + pres_epslon)
    ELSE
        COE0 = 1.0_WP / DT / TALP(NS)
    END IF

    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        COE1 = DXI / RCCI2(JJ)
        COE2 = DYFI(JJ) / RCCI1(JJ)
        COE4 = 1.0_WP / DT / RCCI2(JJ)
        DO KC = 1, NCL3
            KP = KPV(KC)
            !================== MAX. DIV IN THE main DOMAIN =================
            DO IC = 1, NCL1_io
                IP = IPV_io(IC)

                !========== D(\rho u) / DX at (i, J, K) ===========================
                DIVX  = ( G_io(IP, JC, KC, 1) - G_io(IC, JC, KC, 1) ) * COE1

                !========== D(\rho v) / Dy at (i, J, K) ===========================
                DIVY  = ( G_io(IC, JP, KC, 2) - G_io(IC, JC, KC, 2) ) * COE2

                !========== D(\rho w) / Dz at (i, J, K) ===========================
                DIVZ  = ( G_io(IC, JC, KP, 3) - G_io(IC, JC, KC, 3) ) * DZI

                !========== D \RHO / DT ========================================
                !DTRHO = ( DENSITY(IC, JC, KC) - DENSITYP(IC, JC, KC) ) * COE4
                DTRHO = DrhoDtP(IC, JC, KC) / RCCI2(JJ)

                !=========== TOTAL ============================================
                DQCAP = DIVX + DIVY + DIVZ + DTRHO

                RHSLLPHI_io(IC, JC, KC) = DQCAp * COE0

                !                    IF(JJ == 1)    THEN
                !                        RHSLLPHI_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC) + AMPH0 * DPDYWAL(IC, KC, 1) / DYFI(1)
                !                        !WRITE(*, *) JJ, DQCAp * COE0, AMPH0 * DPDYWAL(IC, KC, 1) * DYFI(1), RHSLLPHI_io(IC, JC, KC)
                !                    END IF
                !                    IF(JJ == NCL2) THEN
                !                        RHSLLPHI_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC) - APPH0 * DPDYWAL(IC, KC, 2) / DYFI(NCL2)
                !                        !WRITE(*, *) JJ, DQCAp * COE0, AMPH0 * DPDYWAL(IC, KC, 2) * DYFI(NCL2), RHSLLPHI_io(IC, JC, KC)
                !                    END IF

                !                    IF(MYID == 0) &
                !                    WRITE(*, '(A, 4I3.1,7ES13.5)') 'divgoy', MYID, JC, KC, IC, &
                !                    DIVX / COE1, DIVY / COE2, DIVZ/ DZI, &
                !                    DENSITY(IC, JC, KC), DENSITY0(IC, JC, KC), &
                !                    DTRHO / COE4, DQCAP

            END DO
        END DO
    END DO



    RETURN

END SUBROUTINE

!************************************************************************************************************
SUBROUTINE DIVG_U_io
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: IC, JC, KC!, NYI
    INTEGER(4) :: IP, JP, KP, JJ
    REAL(WP) :: DIVX
    REAL(WP) :: DIVY
    REAL(WP) :: DIVZ
    REAL(WP) :: COE2, COE3


    DivU_io = 0.0_WP


    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        COE2 = DYFI(JJ) * RCCI1(JJ)
        COE3 = DZI * RCCI2(JJ)
        DO KC = 1, NCL3
            KP = KPV(KC)
            !================== MAX. DIV IN THE main DOMAIN =================
            DO IC = 1, NCL1_io
                IP = IPV_io(IC)

                !========== D(\rho u) / DX at (i, J, K) ===========================
                DIVX  = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
                !========== D(\rho v) / Dy at (i, J, K) ===========================
                DIVY  = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * COE2
                !========== D(\rho w) / Dz at (i, J, K) ===========================
                DIVZ  = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * COE3

                !=========== TOTAL ============================================
                DivU_io(IC, JC, KC) = (DIVX + DIVY + DIVZ) / 3.0_WP

                !WRITE(*, '(A, 3I4.1, 4ES13.5)') 'divxyz', JJ, IC, KC, DIVX, DIVY, DIVZ, DivU_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !        !===========================below bottom wall JC = 0 IS not USEd============================
    !        IF(MYID == 0) THEN
    !            JC = 0
    !            JP = 1
    !            JJ = 1
    !            COE2 = DYFI(JJ) * RCCI1(JJ)
    !            COE3 = DZI * RCCI2(JJ)
    !            DO KC = 1, NCL3
    !                KP = KPV(KC)
    !                !================== MAX. DIV IN THE main DOMAIN =================
    !                DO IC = 1, NCL1_io
    !                    IP = IPV_io(IC)

    !                    !========== D(\rho u) / DX at (i, J, K) ===========================
    !                    DIVX  = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
    !                    !========== D(\rho v) / Dy at (i, J, K) ===========================
    !                    DIVY  = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * COE2
    !                    !========== D(\rho w) / Dz at (i, J, K) ===========================
    !                    DIVZ  = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * COE3

    !                    !=========== TOTAL ============================================
    !                    DivU_io(IC, JC, KC) = (DIVX + DIVY + DIVZ) / 3.0_WP

    !                    !WRITE(*, '(A, 3I4.1, 4ES13.5)') 'divxyz', JJ, IC, KC, DIVX, DIVY, DIVZ, DivU_io(IC, JC, KC)
    !                END DO
    !            END DO
    !        END IF

    !        !===========================below top wall JC = 0 IS not USEd============================
    !        IF(MYID == NPSLV) THEN
    !            JC = N2DO(MYID) + 1
    !            JP = N2DO(MYID) + 2
    !            JJ = NND2
    !            COE2 = DYFI(NCL2) * RCCI1(JJ)
    !            COE3 = DZI * RCCI2(JJ)
    !            DO KC = 1, NCL3
    !                KP = KPV(KC)
    !                !================== MAX. DIV IN THE main DOMAIN =================
    !                DO IC = 1, NCL1_io
    !                    IP = IPV_io(IC)

    !                    !========== D(\rho u) / DX at (i, J, K) ===========================
    !                    DIVX  = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
    !                    !========== D(\rho v) / Dy at (i, J, K) ===========================
    !                    DIVY  = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * COE2
    !                    !========== D(\rho w) / Dz at (i, J, K) ===========================
    !                    DIVZ  = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * COE3

    !                    !=========== TOTAL ============================================
    !                    DivU_io(IC, JC, KC) = (DIVX + DIVY + DIVZ) / 3.0_WP

    !                    !WRITE(*, '(A, 3I4.1, 4ES13.5)') 'divxyz', JJ, IC, KC, DIVX, DIVY, DIVZ, DivU_io(IC, JC, KC)
    !                END DO
    !            END DO
    !        END IF

    !===========================================================================================

    !CALL INTFC_ALL_DIVU_io
    CALL INTFC_VARS1(1, NCL1_io, NCL1S, NCL1E, DivU_io)

    RETURN

END SUBROUTINE

!********************************************************************************
SUBROUTINE DENSITY_TIME_DERIVATION(NS)
    USE thermal_info
    USE MESH_INFO
    USE flow_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, JC, KC, JP, JJ, IP, KP, NS
    REAL(WP) :: DIVX, DIVY, DIVZ, COE2, COE3, DIVX0, DIVY0, DIVZ0

    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        COE2 =  DYFI(JJ) * RCCI1(JJ)
        COE3 =  DZI * RCCI2(JJ)
        DO IC = 1, NCL1_io
            IP = IPV_io(IC)
            DO KC = 1, NCL3
                KP = KPV(KC)

                !                    !METHOD 1 ==================================
                !                    !========== D(\rho u) / DX at (i, J, K) ===
                !                    DIVX  = ( G_io(IP, JC, KC, 1) - G_io(IC, JC, KC, 1) ) * DXI

                !                    !========== D(\rho v) / Dy at (i, J, K) ===
                !                    DIVY  = ( G_io(IC, JP, KC, 2) - G_io(IC, JC, KC, 2) ) * COE2

                !                    !========== D(\rho w) / Dz at (i, J, K) ===
                !                    DIVZ  = ( G_io(IC, JC, KP, 3) - G_io(IC, JC, KC, 3) ) * COE3
                !                    !(D - Dp) / DT + GV = 0  dP - D=GV * Dt
                !                    !DENSITYP(IC, JC, KC) = DENSITY(IC, JC, KC) + DT* (DIVX + DIVY+ DIVZ)
                !                    DrhoDtP(IC, JC, KC) = -1.0_WP * (DIVX + DIVY+ DIVZ)

                !METHOD 2 ======GOOD==checK ======================
                !DrhoDtP(IC, JC, KC) = (DENSITY(IC, JC, KC) - DENSITY0(IC, JC, KC)) / DT / TALP(NS)

                !METHOD 4 ====== AdaM -Bashforth method===========
                !DrhoDtP(IC, JC, KC) = (3.0_WP * DENSITY(IC, JC, KC) -4.0_WP * DENSITY0(IC, JC, KC) + DENSITY1(IC, JC, KC) ) * 0.5_WP / DT



                !METHOD 3 =========== RK == (divegence free) = GOOD==== Ref: NICoud1999JCP ============
                !========== D(\rho u) / DX at (i, J, K) ===
                DIVX  = ( G_io(IP, JC, KC, 1) - G_io(IC, JC, KC, 1) ) * DXI
                DIVX0 = ( G0_io(IP, JC, KC, 1) - G0_io(IC, JC, KC, 1) ) * DXI

                !========== D(\rho v) / Dy at (i, J, K) ===
                DIVY  = ( G_io(IC, JP, KC, 2) - G_io(IC, JC, KC, 2) ) * COE2
                DIVY0 = ( G0_io(IC, JP, KC, 2) - G0_io(IC, JC, KC, 2) ) * COE2

                !========== D(\rho w) / Dz at (i, J, K) ===
                DIVZ  = ( G_io(IC, JC, KP, 3) - G_io(IC, JC, KC, 3) ) * COE3
                DIVZ0 = ( G0_io(IC, JC, KP, 3) - G0_io(IC, JC, KC, 3) ) * COE3

                !DrhoDtP(IC, JC, KC) = ( (TGAM(NS) + TROH(NS)) * (DIVX + DIVY + DIVZ) &
                !                            + TROH(NS) * (DIVX0 + DIVY0 + DIVZ0) ) / (TGAM(NS) - 2.0_WP)

                DrhoDtP(IC, JC, KC) = - TGAM(NS) * (DIVX + DIVY + DIVZ) - TROH(NS) * (DIVX0 + DIVY0 + DIVZ0)

            END DO
        END DO
    END DO

    !THIS DENSITY IS USED ONLY IN CONTINUITY EQUATION.


    RETURN

END SUBROUTINE
