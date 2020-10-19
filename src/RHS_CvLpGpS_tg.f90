!**********************************************************************************************************************************
!> @brief
!>       to calculate RHS of the momentum equation. Eq(A6)
!> @details
!> SUBROUTINE: RHS_CvLpGpS_tg (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE RHS_CvLpGpS_tg(NS, IDR)
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4), INTENT(IN) :: IDR

    REAL(WP) :: INTGRHSY
    REAL(WP) :: INTGRHSY_WORK
    !REAL(WP) :: INTGU
    !REAL(WP) :: INTGU_WORK
    REAL(WP) :: DPGRNS
    !REAL(WP) :: INTGVOL
    !REAL(WP) :: INTGVOL_WORK
    !REAL(WP) :: RHSCCC
    REAL(WP) :: SUCACJ
    REAL(WP) :: DERQ
    REAL(WP) :: RHSC, RHSL
    !REAL(WP) :: VOLTMP
    REAL(WP) :: PGM
    INTEGER(4) :: NII
    INTEGER(4) :: I, IC, IM, IP
    INTEGER(4) :: J, JC, JM, JP, JJ
    INTEGER(4) :: K, KC, KM, KP
    REAL(WP) :: RMC2(N2DO(0))
    REAL(WP) :: CONVH_tg(NCL1_tg, N2DO(0), NCL3)
    REAL(WP) :: COE1,COE2

    !>      @note Setup the initial y values
    NII = 1
    IF((IDR == 2) .AND. (MYID == 0)) NII = 2

    RMC2 = 0.0_WP
    IF(IDR == 2)THEN
        DO JC = NII, N2DO(MYID)
            JJ = JCL2G(JC)
            RMC2(JC) = RNDI2(JJ)
        END DO
    ELSE
        DO JC = NII, N2DO(MYID)
            JJ = JCL2G(JC)
            RMC2(JC) = RCCI2(JJ)
        END DO
    END IF

    !======================construct the convectiont term ===============================
    CONVH_tg = 0.0_WP
    IF(IDR == 1) THEN

        DO I = 1, NCL1_tg
            DO J = NII, N2DO(MYID)
                DO K = 1, NCL3
                    CONVH_tg(I, J, K) = QTMP_tg(I, J, K)
                END DO
            END DO
        END DO

    ELSE IF(IDR == 2) THEN

        DO I = 1, NCL1_tg
            DO J = NII, N2DO(MYID)
                DO K = 1, NCL3
                    CONVH_tg(I, J, K) = DPH_tg(I, J, K)
                END DO
            END DO
        END DO

    ELSE IF(IDR == 3) THEN

        DO I = 1, NCL1_tg
            DO J = NII, N2DO(MYID)
                DO K = 1, NCL3
                    CONVH_tg(I, J, K) = RHSLLPHI_tg(I, J, K)
                END DO
            END DO
        END DO

    ELSE
    END IF

    !======================construct viscous term and summing up the convection /viscous terms ===============================
    COE1 = TALP(NS) * CVISC * M_inlet / D_inlet!modIFied by Junjie, 2017/05/ 22
    RHSC = 0.0_WP
    RHSL = 0.0_WP

    DO KC = 1, NCL3
        KM = KMV(KC)
        KP = KPV(KC)
        DO JC = NII, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            DO IC = 1, NCL1_tg
                IP = IPV_tg(IC)
                IM = IMV_tg(IC)
                RHSC  = TGAM(NS) * CONVH_tg(IC, JC, KC) + TROH(NS) * CONVH0_tg(IC, JC, KC, IDR)
                RHSL  = COE1 * ( &
                        (Q_tg(IP, JC, KC, IDR) - 2.0_WP * Q_tg(IC, JC, KC, IDR) + Q_tg(IM, JC, KC, IDR)) * DXQI +  &
                        (Q_tg(IC, JC, KP, IDR) - 2.0_WP * Q_tg(IC, JC, KC, IDR) + Q_tg(IC, JC, KM, IDR)) * DZQI * RMC2(JC) + &
                        (Q_tg(IC, JP, KC, IDR) * APVR(JJ, IDR) + &
                         Q_tg(IC, JC, KC, IDR) * ACVR(JJ, IDR) + &
                         Q_tg(IC, JM, KC, IDR) * AMVR(JJ, IDR) ) )
                CONVH0_tg(IC, JC, KC, IDR) = CONVH_tg(IC, JC, KC)
                RHS_tg(IC, JC, KC) = (RHSC + RHSL) * DT
            END DO
        END DO
    END DO

    !=========================pressure gradient terms ==================================================
    PGM = 0.0_WP
    COE2 = TALP(NS) * DT / D_inlet!modIFied by Junjie, 2017/05/ 22
    IF (IDR == 1) THEN
        DO KC = 1, NCL3
            DO JC = NII, N2DO(MYID)
                DO IC = 1, NCL1_tg
                    IM = IMV_tg(IC)
                    PGM = (PR_tg(IC, JC, KC) - PR_tg(IM, JC, KC)) * DXI * COE2
                    RHS_tg(IC, JC, KC) = RHS_tg(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE IF (IDR == 2) THEN
        DO KC = 1, NCL3
            DO JC = NII, N2DO(MYID)
                JM = JLMV(JC)
                JJ = JCL2G(JC)
                SUCACJ = DYCI(JJ) / RNDI1(JJ)
                DO IC = 1, NCL1_tg
                    PGM = (PR_tg(IC, JC, KC) - PR_tg(IC, JM, KC)) * SUCACJ * COE2
                    RHS_tg(IC, JC, KC) = RHS_tg(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE IF (IDR == 3) THEN
        DO KC = 1, NCL3
            KM = KMV(KC)
            DO JC = NII, N2DO(MYID)
                DO IC = 1, NCL1_tg
                    PGM = (PR_tg(IC, JC, KC) - PR_tg(IC, JC, KM)) * DZI * COE2
                    RHS_tg(IC, JC, KC) = RHS_tg(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE
    ENDIF

    !====================flow drive terms (source terms) in periodic Streamwise flow===========================
    DPGRNS = 0.0_WP
    DERQ = 0.0_WP
    IF (IDR == NFLOW) THEN
        !=============constant mass flow ratE =============================
        IF(iFlowDriven == 1) THEN
            INTGRHSY = 0.0_WP

            DO IC = 1, NCL1_tg
                DO JC = 1, N2DO(MYID)
                    DO KC = 1, NCL3
                        JJ = JCL2G(JC)
                        INTGRHSY = INTGRHSY + RHS_tg(IC, JC, KC) / DYFI(JJ) / RCCI1(JJ)
                    END DO
                END DO
            END DO

            CALL MPI_BARRIER(ICOMM, IERROR)
            CALL MPI_ALLREDUCE(INTGRHSY, INTGRHSY_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
            IF(MYID == 0) DPGRNS = INTGRHSY_WORK / VOLM_tg


        ELSE IF(iFlowDriven == 2) THEN
            IF(MYID == 0) DPGRNS = -0.50_WP * Cf_Given * COE2 !dimensionless based on \Delta and U_m
            !constant pressure gradient, dimensionless based on \Delta and U_m
            !DPGRNS = -2.0_WP      ! DIMENSIONLESS BASED ON U_TAU
        ELSE
        END IF

        CALL MPI_BCAST( DPGRNS, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        DO K = 1, NCL3
            DO I = 1, NCL1_tg
                DO J = NII, N2DO(MYID)
                    RHS_tg(I, J, K) = RHS_tg(I, J, K) - DPGRNS
                END DO
            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE
