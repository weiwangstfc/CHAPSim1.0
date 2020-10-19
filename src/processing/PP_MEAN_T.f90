!**********************************************************************************************************************************
!> @brief
!>        consider below two methods to do a smooth average
!>        -moving average
!>        - Exponential smoothing
!> @details
!> SUBROUTINE: PP_TMEAN_INI_TG (in MYID = all)
!> SUBROUTINE: PP_TMEAN_INI_ZERO_TG
!> SUBROUTINE: PP_TMEAN_INI_io
!> SUBROUTINE: PP_TMEAN_INI_ZERO_io
!> SUBROUTINE: PP_MEAN_T_TG
!> SUBROUTINE: PP_MEAN_T_nonXperiodic_io
!> SUBROUTINE: PP_MEAN_T_Xperiodic_io
!> @note
!> @todo
! REVISION HISTORY:
! 04/2015 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PP_TMEAN_INI_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE

    IF(tRunAve1 > tRunAve_Reset)  tRunAve_Reset = tRunAve1
    IF(iIniField_tg == IniField_reStart .AND. PhyTIME_TG  >   tRunAve1) THEN ! REStart

        CALL ReStart_AVERAGE_VARS_TG
        CALL WRT_AVERAGE_PPED_TG

    ELSE

        CALL PP_TMEAN_INI_ZERO_TG

    END IF

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE PP_TMEAN_INI_ZERO_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE

    NSTATIS_tg = 0

    U1xztL_tg = 0.0_WP
    UPxztL_tg = 0.0_WP

    U2xztL_tg = 0.0_WP
    U3xztL_tg = 0.0_WP

    DVDL1xztL_tg = 0.0_WP
    DVDLPxztL_tg = 0.0_WP
    DVDL2xztL_tg = 0.0_WP
    IF(iPPQuadrants == 1)  THEN
        QUADUVxztL_tg = 0.0_WP
        QUADVzxztL_tg = 0.0_WP
        QUADTKxztL_tg = 0.0_WP
        QUADDRxztL_tg = 0.0_WP
    END IF

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE PP_TMEAN_INI_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    ! iIniField_io == IniField_reStart : reStarting from the last step
    ! PhyTIME_io  >   tRunAve1 : to DO time averaging
    ! iIniFieldType == 0 : ReStart both, flow and thermal fields
    ! iIniFieldTime == 0 :  ReStart following previous time
    IF( (iIniField_io == IniField_reStart) .AND. &
        (PhyTIME_io  >   tRunAve1) .AND. &
        (iIniFieldType == 0) .AND. &
        (iIniFieldTime == 0) ) THEN ! REStart
        IF(TgFlowFlg) THEN
            !=============pp space averaged=======
            !CALL PP_MEAN_Z_FLOW_nonXperiodic_io
            !IF(iThermoDynamics == 1) &
            !CALL PP_MEAN_Z_THEML_nonXperiodic_io
            !============ Read in space and time averaged=========
            CALL ReStart_AVERAGE_VARS_nonXperiodic_io
            !============WRITE out space and time averaged=======
            CALL WRT_AVERAGE_PPED_nonXperiodic_io
            CALL MPI_BARRIER(ICOMM, IERROR)
            !CALL PP_MONITOR_nonXperiodic_io
        ELSE
            !=============pp space averaged=======
            !CALL PP_MEAN_ZX_FLOW_Xperiodic_io
            !(ThIS will be replaced by below)
            !============ Read in space and time averaged=========
            CALL ReStart_AVERAGE_VARS_Xperiodic_io
            IF(tRunAve_Reset < PhyTIME_io) THEN
                PhyTIME = TimeReStart_io
                CALL WRT_AVERAGE_VARS_Xperiodic_io
            END IF
            !============WRITE out space and time averaged=======
            CALL WRT_AVERAGE_PPED_Xperiodic_io
            CALL PP_MEAN_ZX_FLOW_Xperiodic_io
            IF(iPPSpectra == 1) CALL PP_SPECOSPEC
            CALL MPI_BARRIER(ICOMM, IERROR)
            !CALL PP_MONITOR_Xperiodic_io
        END IF

    ELSE

        CALL PP_TMEAN_INI_ZERO_io

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_TMEAN_INI_ZERO_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE


    NSTATIS_io = 0

    IF(TgFlowFlg) THEN
        U1ztL_io = 0.0_WP
        G1ztL_io = 0.0_WP

        UPztL_io = 0.0_WP

        U2ztL_io = 0.0_WP
        UGztL_io = 0.0_WP
        UGUztL_io = 0.0_WP

        DVDL1ztL_io = 0.0_WP
        DVDLPztL_io = 0.0_WP
        DVDL2ztL_io = 0.0_WP

        IF(iThermoDynamics == 1) THEN
            T1ztL_io = 0.0_WP
            D1ztL_io = 0.0_WP
            H1ztL_io = 0.0_WP
            K1ztL_io = 0.0_WP
            M1ztL_io = 0.0_WP

            T2ztL_io = 0.0_WP
            D2ztL_io = 0.0_WP
            H2ztL_io = 0.0_WP
            K2ztL_io = 0.0_WP
            M2ztL_io = 0.0_WP

            UHztL_io = 0.0_WP
            GHztL_io = 0.0_WP
        END IF


    ELSE

        U1xztL_io = 0.0_WP
        G1xztL_io = 0.0_WP
        UPxztL_io = 0.0_WP

        U2xztL_io = 0.0_WP
        UGxztL_io = 0.0_WP
        UGUxztL_io = 0.0_WP
        U3xztL_io = 0.0_WP

        DVDL1xztL_io = 0.0_WP
        DVDLPxztL_io = 0.0_WP
        DVDL2xztL_io = 0.0_WP
        IF(iPPQuadrants == 1)  THEN
            QUADUVxztL_io = 0.0_WP
            QUADVzxztL_io = 0.0_WP
            QUADTKxztL_io = 0.0_WP
            QUADDRxztL_io = 0.0_WP

            QUADDUV1xztL_io = 0.0_WP
            QUADDUV2xztL_io = 0.0_WP
        END IF

        FUxztL_io = 0.0_WP


        IF(iThermoDynamics == 1) THEN
            T1xztL_io = 0.0_WP
            D1xztL_io = 0.0_WP
            H1xztL_io = 0.0_WP
            M1xztL_io = 0.0_WP

            T2xztL_io = 0.0_WP
            D2xztL_io = 0.0_WP
            H2xztL_io = 0.0_WP

            DHxztL_io = 0.0_WP
            PHxztL_io = 0.0_WP

            DVDL1MxztL_io = 0.0_WP
            DVDL1MHxztL_io = 0.0_WP
            DVDL1MUxztL_io = 0.0_WP
            DVDL2MxztL_io = 0.0_WP

            UHxztL_io = 0.0_WP
            GHxztL_io = 0.0_WP
            U2DHxztL_io = 0.0_WP

            DhDL1xztL_io = 0.0_WP
            DhDLPxztL_io = 0.0_WP
            DTDLKxztL_io = 0.0_WP
            DTDLKUxztL_io = 0.0_WP
            DTDLKDVDLxztL_io = 0.0_WP
            DHDLMDVDLxztL_io = 0.0_WP
        END IF

        IF(MYID == 0) THEN
            !=========== Spectra ======================
            !================ Velocity =====================
            IF(iPPSpectra == 1) THEN
                R11X1_xztLA = 0.0_WP
                R22X1_xztLA = 0.0_WP
                R33X1_xztLA = 0.0_WP
                R12X1_xztLA = 0.0_WP
                R13X1_xztLA = 0.0_WP
                R23X1_xztLA = 0.0_WP

                R11X3_xztLA = 0.0_WP
                R22X3_xztLA = 0.0_WP
                R33X3_xztLA = 0.0_WP
                R12X3_xztLA = 0.0_WP
                R13X3_xztLA = 0.0_WP
                R23X3_xztLA = 0.0_WP

                ENE11T_xztLA = 0.0_WP
                ENE22T_xztLA = 0.0_WP
                ENE33T_xztLA = 0.0_WP
                ENE12T_xztLA = 0.0_WP
                ENE13T_xztLA = 0.0_WP
                ENE23T_xztLA = 0.0_WP

                ENE11Z_xztLA = 0.0_WP
                ENE22Z_xztLA = 0.0_WP
                ENE33Z_xztLA = 0.0_WP
                ENE12Z_xztLA = 0.0_WP
                ENE13Z_xztLA = 0.0_WP
                ENE23Z_xztLA = 0.0_WP

                !================ VoritICity ====================
                V11X1_xztLA = 0.0_WP
                V22X1_xztLA = 0.0_WP
                V33X1_xztLA = 0.0_WP
                V12X1_xztLA = 0.0_WP
                V13X1_xztLA = 0.0_WP
                V23X1_xztLA = 0.0_WP

                V11X3_xztLA = 0.0_WP
                V22X3_xztLA = 0.0_WP
                V33X3_xztLA = 0.0_WP
                V12X3_xztLA = 0.0_WP
                V13X3_xztLA = 0.0_WP
                V23X3_xztLA = 0.0_WP

                ENV11T_xztLA = 0.0_WP
                ENV22T_xztLA = 0.0_WP
                ENV33T_xztLA = 0.0_WP
                ENV12T_xztLA = 0.0_WP
                ENV13T_xztLA = 0.0_WP
                ENV23T_xztLA = 0.0_WP

                ENV11Z_xztLA = 0.0_WP
                ENV22Z_xztLA = 0.0_WP
                ENV33Z_xztLA = 0.0_WP
                ENV12Z_xztLA = 0.0_WP
                ENV13Z_xztLA = 0.0_WP
                ENV23Z_xztLA = 0.0_WP


                !=============== VoritICity & VelocitY =========================
                VO11X1_xztLA = 0.0_WP
                VO12X1_xztLA = 0.0_WP
                VO13X1_xztLA = 0.0_WP

                VO21X1_xztLA = 0.0_WP
                VO22X1_xztLA = 0.0_WP
                VO23X1_xztLA = 0.0_WP

                VO31X1_xztLA = 0.0_WP
                VO32X1_xztLA = 0.0_WP
                VO33X1_xztLA = 0.0_WP

                VO11X3_xztLA = 0.0_WP
                VO12X3_xztLA = 0.0_WP
                VO13X3_xztLA = 0.0_WP

                VO21X3_xztLA = 0.0_WP
                VO22X3_xztLA = 0.0_WP
                VO23X3_xztLA = 0.0_WP

                VO31X3_xztLA = 0.0_WP
                VO32X3_xztLA = 0.0_WP
                VO33X3_xztLA = 0.0_WP

                EVO11T_xztLA = 0.0_WP
                EVO12T_xztLA = 0.0_WP
                EVO13T_xztLA = 0.0_WP

                EVO21T_xztLA = 0.0_WP
                EVO22T_xztLA = 0.0_WP
                EVO23T_xztLA = 0.0_WP

                EVO31T_xztLA = 0.0_WP
                EVO32T_xztLA = 0.0_WP
                EVO33T_xztLA = 0.0_WP

                EVO11Z_xztLA = 0.0_WP
                EVO12Z_xztLA = 0.0_WP
                EVO13Z_xztLA = 0.0_WP

                EVO21Z_xztLA = 0.0_WP
                EVO22Z_xztLA = 0.0_WP
                EVO23Z_xztLA = 0.0_WP

                EVO31Z_xztLA = 0.0_WP
                EVO32Z_xztLA = 0.0_WP
                EVO33Z_xztLA = 0.0_WP
            END IF

        END IF
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_T_TG(TP)
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: TP
    REAL(WP) :: COE1, COE2, COE3

    NSTATIS_TG = NSTATIS_TG + 1

    IF(TP == T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
        !COE1 = DBLE(NSTATIS_io - 1)
        COE2 = 1.0_WP / DBLE(NSTATIS_TG)
        COE3 = DBLE(NSTATIS_TG - 1) / DBLE(NSTATIS_TG)
    ELSE IF(TP == T_Summing_average) THEN ! just adding all them together
        !COE1 = 1.0_WP /
        COE2 = 1.0_WP / DBLE(pp_instn_sz)
        COE3 = 1.0_WP
    ELSE
    END IF


    !        NSTATIS_TG = NSTATIS_TG + 1

    !        COE1 = DBLE(NSTATIS_TG - 1)
    !        COE2 = 1.0_WP / DBLE(NSTATIS_TG)
    !        COE3 = DBLE(NSTATIS_TG - 1) / DBLE(NSTATIS_TG)


    U1xztL_tg = COE3 * U1xztL_tg + COE2 * U1xzL_tg
    UPxztL_tg = COE3 * UPxztL_tg + COE2 * UPxzL_tg

    U2xztL_tg = COE3 * U2xztL_tg + COE2 * U2xzL_tg
    U3xztL_tg = COE3 * U3xztL_tg + COE2 * U3xzL_tg

    DVDL1xztL_tg = COE3 * DVDL1xztL_tg + COE2 * DVDL1xzL_tg
    DVDLPxztL_tg = COE3 * DVDLPxztL_tg + COE2 * DVDLPxzL_tg
    DVDL2xztL_tg = COE3 * DVDL2xztL_tg + COE2 * DVDL2xzL_tg
    IF(iPPQuadrants == 1)  THEN
        QUADUVxztL_tg = COE3 * QUADUVxztL_tg + COE2 * QUADUVxzL_tg
        QUADVzxztL_tg = COE3 * QUADVzxztL_tg + COE2 * QUADVzxzL_tg
        QUADTKxztL_tg = COE3 * QUADTKxztL_tg + COE2 * QUADTKxzL_tg
        QUADDRxztL_tg = COE3 * QUADDRxztL_tg + COE2 * QUADDRxzL_tg
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_T_nonXperiodic_io(TP)
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: TP
    REAL(WP) :: COE1, COE2, COE3


    NSTATIS_io = NSTATIS_io + 1

    IF(TP == T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
        !COE1 = DBLE(NSTATIS_io - 1)
        COE2 = 1.0_WP / DBLE(NSTATIS_io)
        COE3 = DBLE(NSTATIS_io - 1) / DBLE(NSTATIS_io)
    ELSE IF(TP == T_Summing_average) THEN ! just adding all them together
        !COE1 = 1.0_WP /
        COE2 = 1.0_WP / DBLE(pp_instn_sz)
        COE3 = 1.0_WP
    ELSE
    END IF


    !COE1 = DBLE(NSTATIS_io - 1)
    !COE2 = 1.0_WP / DBLE(NSTATIS_io)
    !COE3 = DBLE(NSTATIS_io - 1) / DBLE(NSTATIS_io)

    U1ztL_io = COE3 * U1ztL_io + COE2 * U1zL_io
    G1ztL_io = COE3 * G1ztL_io + COE2 * G1zL_io

    UPztL_io = COE3 * UPztL_io + COE2 * UPzL_io

    U2ztL_io = COE3 * U2ztL_io  + COE2 * U2zL_io
    UGztL_io = COE3 * UGztL_io  + COE2 * UGzL_io
    UGUztL_io = COE3 * UGUztL_io + COE2 * UGUzL_io

    DVDL1ztL_io = COE3 * DVDL1ztL_io + COE2 * DVDL1zL_io
    DVDLPztL_io = COE3 * DVDLPztL_io + COE2 * DVDLPzL_io
    DVDL2ztL_io = COE3 * DVDL2ztL_io + COE2 * DVDL2zL_io

    IF(iThermoDynamics == 1) THEN
        T1ztL_io = COE3 * T1ztL_io + COE2 * T1zL_io
        D1ztL_io = COE3 * D1ztL_io + COE2 * D1zL_io
        H1ztL_io = COE3 * H1ztL_io + COE2 * H1zL_io
        K1ztL_io = COE3 * K1ztL_io + COE2 * K1zL_io
        M1ztL_io = COE3 * M1ztL_io + COE2 * M1zL_io

        T2ztL_io = COE3 * T2ztL_io + COE2 * T2zL_io
        D2ztL_io = COE3 * D2ztL_io + COE2 * D2zL_io
        H2ztL_io = COE3 * H2ztL_io + COE2 * H2zL_io
        K2ztL_io = COE3 * K2ztL_io + COE2 * K2zL_io
        M2ztL_io = COE3 * M2ztL_io + COE2 * M2zL_io

        DHztL_io = COE3 * DHztL_io + COE2 * DHzL_io

        UHztL_io = COE3 * UHztL_io + COE2 * UHzL_io
        GHztL_io = COE3 * GHztL_io + COE2 * GHzL_io
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_T_Xperiodic_io(TP)
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: TP
    REAL(WP) :: COE2, COE3


    NSTATIS_io = NSTATIS_io + 1

    IF(TP == T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
        !COE1 = DBLE(NSTATIS_io - 1)
        COE2 = 1.0_WP / DBLE(NSTATIS_io)
        COE3 = DBLE(NSTATIS_io - 1) / DBLE(NSTATIS_io)
    ELSE IF(TP == T_Summing_average) THEN ! just adding all them together
        !COE1 = 1.0_WP /
        COE2 = 1.0_WP / DBLE(pp_instn_sz)
        COE3 = 1.0_WP
    ELSE
    END IF

    U1xztL_io = COE3 * U1xztL_io + COE2 * U1xzL_io
    G1xztL_io = COE3 * G1xztL_io + COE2 * G1xzL_io
    UPxztL_io = COE3 * UPxztL_io + COE2 * UPxzL_io

    U2xztL_io = COE3 * U2xztL_io + COE2 * U2xzL_io
    UGxztL_io = COE3 * UGxztL_io + COE2 * UGxzL_io

    U3xztL_io = COE3 * U3xztL_io + COE2 * U3xzL_io
    UGUxztL_io = COE3 * UGUxztL_io + COE2 * UGUxzL_io

    DVDL1xztL_io = COE3 * DVDL1xztL_io + COE2 * DVDL1xzL_io
    DVDLPxztL_io = COE3 * DVDLPxztL_io + COE2 * DVDLPxzL_io
    DVDL2xztL_io = COE3 * DVDL2xztL_io + COE2 * DVDL2xzL_io
    IF(iPPQuadrants == 1)  THEN
        QUADUVxztL_io = COE3 * QUADUVxztL_io + COE2 * QUADUVxzL_io
        QUADVzxztL_io = COE3 * QUADVzxztL_io + COE2 * QUADVzxzL_io
        QUADTKxztL_io = COE3 * QUADTKxztL_io + COE2 * QUADTKxzL_io
        QUADDRxztL_io = COE3 * QUADDRxztL_io + COE2 * QUADDRxzL_io
        QUADDUV1xztL_io = COE3 * QUADDUV1xztL_io + COE2 * QUADDUV1xzL_io
        QUADDUV2xztL_io = COE3 * QUADDUV2xztL_io + COE2 * QUADDUV2xzL_io

        OCTDUVxztL_io = COE3 * OCTDUVxztL_io + COE2 * OCTDUVxzL_io
        OCTDVzxztL_io = COE3 * OCTDVzxztL_io + COE2 * OCTDVzxzL_io
        OCTDTKxztL_io = COE3 * OCTDTKxztL_io + COE2 * OCTDTKxzL_io
        OCTDDRxztL_io = COE3 * OCTDDRxztL_io + COE2 * OCTDDRxzL_io
        OCTDDUV1xztL_io = COE3 * OCTDDUV1xztL_io + COE2 * OCTDDUV1xzL_io
        OCTDDUV2xztL_io = COE3 * OCTDDUV2xztL_io + COE2 * OCTDDUV2xzL_io

        OCTTUVxztL_io = COE3 * OCTTUVxztL_io + COE2 * OCTTUVxzL_io
        OCTTVzxztL_io = COE3 * OCTTVzxztL_io + COE2 * OCTTVzxzL_io
        OCTTTKxztL_io = COE3 * OCTTTKxztL_io + COE2 * OCTTTKxzL_io
        OCTTDRxztL_io = COE3 * OCTTDRxztL_io + COE2 * OCTTDRxzL_io
        OCTTDUV1xztL_io = COE3 * OCTTDUV1xztL_io + COE2 * OCTTDUV1xzL_io
        OCTTDUV2xztL_io = COE3 * OCTTDUV2xztL_io + COE2 * OCTTDUV2xzL_io
    END IF

    FUxztL_io    = COE3 * FUxztL_io    + COE2 * FUxzL_io


    IF(iThermoDynamics == 1) THEN
        T1xztL_io = COE3 * T1xztL_io + COE2 * T1xzL_io
        D1xztL_io = COE3 * D1xztL_io + COE2 * D1xzL_io
        H1xztL_io = COE3 * H1xztL_io + COE2 * H1xzL_io
        M1xztL_io = COE3 * M1xztL_io + COE2 * M1xzL_io

        T2xztL_io = COE3 * T2xztL_io + COE2 * T2xzL_io
        D2xztL_io = COE3 * D2xztL_io + COE2 * D2xzL_io
        H2xztL_io = COE3 * H2xztL_io + COE2 * H2xzL_io

        DHxztL_io = COE3 * DHxztL_io + COE2 * DHxzL_io
        PHxztL_io = COE3 * PHxztL_io + COE2 * PHxzL_io

        DVDL1MxztL_io = COE3 * DVDL1MxztL_io + COE2 * DVDL1MxzL_io
        DVDL1MHxztL_io = COE3 * DVDL1MHxztL_io + COE2 * DVDL1MHxzL_io
        DVDL1MUxztL_io = COE3 * DVDL1MUxztL_io + COE2 * DVDL1MUxzL_io
        DVDL2MxztL_io = COE3 * DVDL2MxztL_io + COE2 * DVDL2MxzL_io

        UHxztL_io = COE3 * UHxztL_io + COE2 * UHxzL_io
        GHxztL_io = COE3 * GHxztL_io + COE2 * GHxzL_io
        U2DHxztL_io = COE3 * U2DHxztL_io + COE2 * U2DHxzL_io

        DhDL1xztL_io = COE3 * DhDL1xztL_io + COE2 * DhDL1xzL_io
        DhDLPxztL_io = COE3 * DhDLPxztL_io + COE2 * DhDLPxzL_io
        DTDLKxztL_io = COE3 * DTDLKxztL_io + COE2 * DTDLKxzL_io
        DTDLKUxztL_io = COE3 * DTDLKUxztL_io + COE2 * DTDLKUxzL_io

        DTDLKDVDLxztL_io = COE3 * DTDLKDVDLxztL_io + COE2 * DTDLKDVDLxzL_io
        DHDLMDVDLxztL_io = COE3 * DHDLMDVDLxztL_io + COE2 * DHDLMDVDLxzL_io
    END IF

    IF(MYID == 0) THEN
        !=========== Spectra ======================
        !================ Velocity =====================
        IF(iPPSpectra == 1) THEN
            R11X1_xztLa = COE3 * R11X1_xztLa + COE2 * R11X1_xzLa
            R22X1_xztLa = COE3 * R22X1_xztLa + COE2 * R22X1_xzLa
            R33X1_xztLa = COE3 * R33X1_xztLa + COE2 * R33X1_xzLa
            R12X1_xztLa = COE3 * R12X1_xztLa + COE2 * R12X1_xzLa
            R13X1_xztLa = COE3 * R13X1_xztLa + COE2 * R13X1_xzLa
            R23X1_xztLa = COE3 * R23X1_xztLa + COE2 * R23X1_xzLa

            R11X3_xztLa = COE3 * R11X3_xztLa + COE2 * R11X3_xzLa
            R22X3_xztLa = COE3 * R22X3_xztLa + COE2 * R22X3_xzLa
            R33X3_xztLa = COE3 * R33X3_xztLa + COE2 * R33X3_xzLa
            R12X3_xztLa = COE3 * R12X3_xztLa + COE2 * R12X3_xzLa
            R13X3_xztLa = COE3 * R13X3_xztLa + COE2 * R13X3_xzLa
            R23X3_xztLa = COE3 * R23X3_xztLa + COE2 * R23X3_xzLa

            ENE11T_xztLa = COE3 * ENE11T_xztLa + COE2 * ENE11T_xzLa
            ENE22T_xztLa = COE3 * ENE22T_xztLa + COE2 * ENE22T_xzLa
            ENE33T_xztLa = COE3 * ENE33T_xztLa + COE2 * ENE33T_xzLa
            ENE12T_xztLa = COE3 * ENE12T_xztLa + COE2 * ENE12T_xzLa
            ENE13T_xztLa = COE3 * ENE13T_xztLa + COE2 * ENE13T_xzLa
            ENE23T_xztLa = COE3 * ENE23T_xztLa + COE2 * ENE23T_xzLa

            ENE11Z_xztLa = COE3 * ENE11Z_xztLa + COE2 * ENE11Z_xzLa
            ENE22Z_xztLa = COE3 * ENE22Z_xztLa + COE2 * ENE22Z_xzLa
            ENE33Z_xztLa = COE3 * ENE33Z_xztLa + COE2 * ENE33Z_xzLa
            ENE12Z_xztLa = COE3 * ENE12Z_xztLa + COE2 * ENE12Z_xzLa
            ENE13Z_xztLa = COE3 * ENE13Z_xztLa + COE2 * ENE13Z_xzLa
            ENE23Z_xztLa = COE3 * ENE23Z_xztLa + COE2 * ENE23Z_xzLa

            !================ VoritICity ====================
            V11X1_xztLa = COE3 * V11X1_xztLa + COE2 * V11X1_xzLa
            V22X1_xztLa = COE3 * V22X1_xztLa + COE2 * V22X1_xzLa
            V33X1_xztLa = COE3 * V33X1_xztLa + COE2 * V33X1_xzLa
            V12X1_xztLa = COE3 * V12X1_xztLa + COE2 * V12X1_xzLa
            V13X1_xztLa = COE3 * V13X1_xztLa + COE2 * V13X1_xzLa
            V23X1_xztLa = COE3 * V23X1_xztLa + COE2 * V23X1_xzLa

            V11X3_xztLa = COE3 * V11X3_xztLa + COE2 * V11X3_xzLa
            V22X3_xztLa = COE3 * V22X3_xztLa + COE2 * V22X3_xzLa
            V33X3_xztLa = COE3 * V33X3_xztLa + COE2 * V33X3_xzLa
            V12X3_xztLa = COE3 * V12X3_xztLa + COE2 * V12X3_xzLa
            V13X3_xztLa = COE3 * V13X3_xztLa + COE2 * V13X3_xzLa
            V23X3_xztLa = COE3 * V23X3_xztLa + COE2 * V23X3_xzLa

            ENV11T_xztLa = COE3 * ENV11T_xztLa + COE2 * ENV11T_xzLa
            ENV22T_xztLa = COE3 * ENV22T_xztLa + COE2 * ENV22T_xzLa
            ENV33T_xztLa = COE3 * ENV33T_xztLa + COE2 * ENV33T_xzLa
            ENV12T_xztLa = COE3 * ENV12T_xztLa + COE2 * ENV12T_xzLa
            ENV13T_xztLa = COE3 * ENV13T_xztLa + COE2 * ENV13T_xzLa
            ENV23T_xztLa = COE3 * ENV23T_xztLa + COE2 * ENV23T_xzLa

            ENV11Z_xztLa = COE3 * ENV11Z_xztLa + COE2 * ENV11Z_xzLa
            ENV22Z_xztLa = COE3 * ENV22Z_xztLa + COE2 * ENV22Z_xzLa
            ENV33Z_xztLa = COE3 * ENV33Z_xztLa + COE2 * ENV33Z_xzLa
            ENV12Z_xztLa = COE3 * ENV12Z_xztLa + COE2 * ENV12Z_xzLa
            ENV13Z_xztLa = COE3 * ENV13Z_xztLa + COE2 * ENV13Z_xzLa
            ENV23Z_xztLa = COE3 * ENV23Z_xztLa + COE2 * ENV23Z_xzLa


            !=============== VoritICity & VelocitY =========================
            VO11X1_xztLa = COE3 * VO11X1_xztLa + COE2 * VO11X1_xzLa
            VO12X1_xztLa = COE3 * VO12X1_xztLa + COE2 * VO12X1_xzLa
            VO13X1_xztLa = COE3 * VO13X1_xztLa + COE2 * VO13X1_xzLa

            VO21X1_xztLa = COE3 * VO21X1_xztLa + COE2 * VO21X1_xzLa
            VO22X1_xztLa = COE3 * VO22X1_xztLa + COE2 * VO22X1_xzLa
            VO23X1_xztLa = COE3 * VO23X1_xztLa + COE2 * VO23X1_xzLa

            VO31X1_xztLa = COE3 * VO31X1_xztLa + COE2 * VO31X1_xzLa
            VO32X1_xztLa = COE3 * VO32X1_xztLa + COE2 * VO32X1_xzLa
            VO33X1_xztLa = COE3 * VO33X1_xztLa + COE2 * VO33X1_xzLa

            VO11X3_xztLa = COE3 * VO11X3_xztLa + COE2 * VO11X3_xzLa
            VO12X3_xztLa = COE3 * VO12X3_xztLa + COE2 * VO12X3_xzLa
            VO13X3_xztLa = COE3 * VO13X3_xztLa + COE2 * VO13X3_xzLa

            VO21X3_xztLa = COE3 * VO21X3_xztLa + COE2 * VO21X3_xzLa
            VO22X3_xztLa = COE3 * VO22X3_xztLa + COE2 * VO22X3_xzLa
            VO23X3_xztLa = COE3 * VO23X3_xztLa + COE2 * VO23X3_xzLa

            VO31X3_xztLa = COE3 * VO31X3_xztLa + COE2 * VO31X3_xzLa
            VO32X3_xztLa = COE3 * VO32X3_xztLa + COE2 * VO32X3_xzLa
            VO33X3_xztLa = COE3 * VO33X3_xztLa + COE2 * VO33X3_xzLa

            EVO11T_xztLa = COE3 * EVO11T_xztLa + COE2 * EVO11T_xzLa
            EVO12T_xztLa = COE3 * EVO12T_xztLa + COE2 * EVO12T_xzLa
            EVO13T_xztLa = COE3 * EVO13T_xztLa + COE2 * EVO13T_xzLa

            EVO21T_xztLa = COE3 * EVO21T_xztLa + COE2 * EVO21T_xzLa
            EVO22T_xztLa = COE3 * EVO22T_xztLa + COE2 * EVO22T_xzLa
            EVO23T_xztLa = COE3 * EVO23T_xztLa + COE2 * EVO23T_xzLa

            EVO31T_xztLa = COE3 * EVO31T_xztLa + COE2 * EVO31T_xzLa
            EVO32T_xztLa = COE3 * EVO32T_xztLa + COE2 * EVO32T_xzLa
            EVO33T_xztLa = COE3 * EVO33T_xztLa + COE2 * EVO33T_xzLa

            EVO11Z_xztLa = COE3 * EVO11Z_xztLa + COE2 * EVO11Z_xzLa
            EVO12Z_xztLa = COE3 * EVO12Z_xztLa + COE2 * EVO12Z_xzLa
            EVO13Z_xztLa = COE3 * EVO13Z_xztLa + COE2 * EVO13Z_xzLa

            EVO21Z_xztLa = COE3 * EVO21Z_xztLa + COE2 * EVO21Z_xzLa
            EVO22Z_xztLa = COE3 * EVO22Z_xztLa + COE2 * EVO22Z_xzLa
            EVO23Z_xztLa = COE3 * EVO23Z_xztLa + COE2 * EVO23Z_xzLa

            EVO31Z_xztLa = COE3 * EVO31Z_xztLa + COE2 * EVO31Z_xzLa
            EVO32Z_xztLa = COE3 * EVO32Z_xztLa + COE2 * EVO32Z_xzLa
            EVO33Z_xztLa = COE3 * EVO33Z_xztLa + COE2 * EVO33Z_xzLa
        END IF
    END IF

    RETURN
END SUBROUTINE
