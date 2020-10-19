!**********************************************************************************************************************************
!> @brief
!>       Factorization Approximation of the momentume equation. Eq A4a
!> @details
!> SUBROUTINE: MOMFA_tg (in MYID = all)
!> SUBROUTINE: MOMFA1_X_tg (in MYID = all)
!> SUBROUTINE: MOMFA2_Y_tg (in MYID = all)
!> SUBROUTINE: MOMFA3_Z_tg (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE MOMFA_tg(NS, IDR)
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4), INTENT(IN) :: IDR

    INTEGER(4) :: N2I
    INTEGER(4) :: I, J, K


    FACOE_tg = 0.50_WP * TALP(NS) * DT * CVISC * M_inlet / D_inlet     != Alpha* Dt/ (2Re), Junjie

    N2I = 1
    IF ((MYID == 0) .AND. (IDR == 2))  N2I = 2

    CALL MOMFA1_X_tg(N2I, IDR)
    CALL MOMFA2_Z_tg(N2I)
    CALL MOMFA3_Y_tg(IDR)

    DO K = 1, NCL3
        DO J = N2I, N2DO(MYID)
            DO I = 1, NCL1_tg
                Q_tg(I, J, K, IDR) = RHS_tg(I, J, K) + Q_tg(I, J, K, IDR)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE MOMFA1_X_tg(N2I, IDR)
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N2I
    INTEGER(4), INTENT(IN) :: IDR

    REAL(WP) :: AMIV(NCL1_tg, N2DO(0))
    REAL(WP) :: ACIV(NCL1_tg, N2DO(0))
    REAL(WP) :: APIV(NCL1_tg, N2DO(0))
    REAL(WP) :: FI  (NCL1_tg, N2DO(0))
    REAL(WP) :: COE

    INTEGER(4) :: I, J, K, JSZ

    AMIV = 0.0_WP
    ACIV = 0.0_WP
    APIV = 0.0_WP
    FI  = 0.0_WP

    COE = - FACOE_tg * DXQI
    DO K = 1, NCL3


        DO J = N2I, N2DO(MYID)
            DO I = 1, NCL1_tg
                APIV(I, J) = COE     !ci
                AMIV(I, J) = COE      !AI
                ACIV(I, J) = 1.0_WP - APIV(I, J) - AMIV(I, J)             !bi
                FI(I, J) = RHS_tg(I, J, K)       !RHSi
            END DO
        END DO
        !             CALL TRIPVI(AMI, ACI, API, FI, 1, N1M, N2I, N2DO)
        JSZ = N2DO(MYID) - N2I + 1
        CALL TDMAIJI_CYC(AMIV(1 : NCL1_tg, N2I : N2DO(MYID)), &
                         ACIV(1 : NCL1_tg, N2I : N2DO(MYID)), &
                         APIV(1 : NCL1_tg, N2I : N2DO(MYID)), &
                           FI(1 : NCL1_tg, N2I : N2DO(MYID)), &
                         1, NCL1_tg, N2I, JSZ)
        DO J = N2I, N2DO(MYID)
            DO I = 1, NCL1_tg
                RHS_tg(I, J, K) = FI(I, J)
            END DO
        END DO

        IF(IoFlowFlg) THEN
            DO J = N2I, N2DO(MYID)
                BC_U_SSTAR(J, K, IDR) = RHS_tg(NCL1_tg, J, K)
            END DO
        END IF

    END DO

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE MOMFA2_Z_tg(N2I)
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N2I

    REAL(WP) :: AMKV(NCL3, N2DO(0))
    REAL(WP) :: ACKV(NCL3, N2DO(0))
    REAL(WP) :: APKV(NCL3, N2DO(0))
    REAL(WP) :: FK  (NCL3, N2DO(0))

    REAL(WP) :: RMC2(N2DO(MYID))

    INTEGER(4) :: I, J, K, JJ, JSZ

    AMKV = 0.0_WP
    ACKV = 0.0_WP
    APKV = 0.0_WP
    FK  = 0.0_WP

    RMC2 = 0.0_WP

    DO J = N2I, N2DO(MYID)
        JJ = JCL2G(J)
        IF (N2I == 2) THEN
            RMC2(J) = RNDI2(JJ) * (-FACOE_tg * DZQI)
        ELSE
            RMC2(J) = RCCI2(JJ) * (-FACOE_tg * DZQI)
        END IF
    END DO

    DO I = 1, NCL1_tg

        DO K = 1, NCL3
            DO J = N2I, N2DO(MYID)
                APKV(K, J) = RMC2(J)    !@
                AMKV(K, J) = RMC2(J)    !@
                ACKV(K, J) = 1.0_WP - APKV(K, J) - AMKV(K, J)
                FK(K, J) = RHS_tg(I, J, K)
            END DO
        END DO

        !         CALL  TRVPJK(AMK, ACK, APK, FK, 1, N3M, N2I, N2DO )
        JSZ = N2DO(MYID) - N2I + 1
        CALL TDMAIJI_CYC(AMKV(1 : NCL3, N2I : N2DO(MYID)), &
        ACKV(1 : NCL3, N2I : N2DO(MYID)), &
        APKV(1 : NCL3, N2I : N2DO(MYID)), &
        FK(1 : NCL3, N2I : N2DO(MYID)), &
        1, NCL3, N2I, JSZ)
        DO K = 1, NCL3
            DO J = N2I, N2DO(MYID)
                RHS_tg(I, J, K) = FK(K, J)
            END DO
        END DO

    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE MOMFA3_Y_tg(IDR)
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IDR

    INTEGER(4) :: I, J, K
    REAL(WP) :: FJJ (NCL1_tg, NND2)
    REAL(WP) :: AMJV(NCL1_tg, NND2)
    REAL(WP) :: ACJV(NCL1_tg, NND2)
    REAL(WP) :: APJV(NCL1_tg, NND2)
    REAL(WP) :: BCJJ(NCL1_tg, 2)

    REAL(WP) :: F(NCL1_TG, NCL2, N3DO(0) )

    INTEGER(4) :: NYS, NYE, JSZ


    !IF(IDR == 2) THEN
    !NYS = 1
    !NYE = NND2
    !ELSE
    NYS = 1
    NYE = NCL2
    !END IF
    FJJ = 0.0_WP
    AMJV = 0.0_WP
    APJV = 0.0_WP
    ACJV = 0.0_WP
    BCJJ = 0.0_WP
    F  = 0.0_WP

    !CALL TRASP23L2G_RHS
    CALL TRASP23_Y2Z(NCL1_tg, 1, N2DO(0), RHS_TG, F)

    DO K = 1, N3DO(MYID)

        DO I = 1, NCL1_tg
            DO J = NYS, NYE
                FJJ(I, J) = F(I, J, K)
                ACJV(I, J) = 1.0_WP - FACOE_tg * ACVR(J, IDR)
                APJV(I, J) = - FACOE_tg * APVR(J, IDR)
                AMJV(I, J) = - FACOE_tg * AMVR(J, IDR)
            END DO
            BCJJ(I, :) = 0.0_WP
        END DO

        JSZ = NYE - NYS + 1
        CALL TDMAIJJ_nonCYC (AMJV(1 : NCL1_tg, NYS : NYE), &
        ACJV(1 : NCL1_tg, NYS : NYE), &
        APJV(1 : NCL1_tg, NYS : NYE), &
        FJJ(1 : NCL1_tg, NYS : NYE), &
        BCJJ(1 : NCL1_tg, 1:2), &
        NYS, JSZ, 1, NCL1_tg)

        DO I = 1, NCL1_tg
            DO J = NYS, NYE
                F(I, J, K) = FJJ(I, J)
            END DO
        END DO
    END DO

    !CALL TRASP23G2L_RHS
    CALL TRASP23_Z2Y(NCL1_TG, 1, N2DO(0), RHS_TG, F)

    RETURN
END SUBROUTINE
