!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE: POSTPROCESS_TG (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE POSTPROCESS_TG
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE


    !======================== MONITORING THE CALCULATION RUNNING PROPERLY ========================
    CALL DIVGCK_tg

    !=============== CALL average every time step after the setting tsav1 =========================
    IF(PhyTIME > tRunAve1) THEN
        CALL PP_MEAN_ZX_FLOW_TG
        CALL PP_MEAN_T_TG(T_Asymptotic_Average)
    END IF

    !=============== Screen printing for monitoring the codE ======================================
    IF (ITERG == 1 .OR. MOD(ITERG, iterMonitor) == 0)  THEN
        IF(PhyTIME <= tRunAve1) THEN
            CALL PP_MEAN_ZX_FLOW_TG
        END IF
        CALL PP_MONITOR_TG
    END IF

    !=================WRITE out data============================================================

    IF( DMOD(PhyTIME, dtSave) < DT) THEN
        !================WRITE out basIC variables =================
        CALL WRT_INSTANT_VARS_TG
        !================WRITE out time averaged variables ===========
        IF(PhyTIME > tRunAve1) THEN
            CALL WRT_AVERAGE_VARS_TG
            IF(DMOD(PhyTIME, dtAveView)  <  DT) THEN
                CALL WRT_AVERAGE_PPED_TG
            END IF
        END IF

    END IF

    RETURN
END SUBROUTINE

!*******************************************************************

SUBROUTINE POSTPROCESS_io
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE


    !============ instantanous saving===================================
    IF( DMOD(PhyTIME, dtSave)  < DT) CALL WRT_INSTANT_VARS_io

    !=============== CALL average every time step after the setting tsav1 =========================
    IF(PhyTIME > tRunAve1) THEN

        !=================== instantanous processing for each steP ==========================
        IF(TgFlowFlg) THEN
            CALL PP_MEAN_Z_FLOW_nonXperiodic_io
            IF(iThermoDynamics == 1) &
            CALL PP_MEAN_Z_THEML_nonXperiodic_io
            CALL PP_MEAN_T_nonXperiodic_io(T_Asymptotic_Average)
        ELSE
            CALL PP_MEAN_ZX_FLOW_Xperiodic_io
            IF(iPPSpectra == 1) CALL PP_SPECOSPEC
            CALL PP_MEAN_T_Xperiodic_io(T_Asymptotic_Average)
        END IF

        !==================WRiting ouT ========================
        IF( DMOD(PhyTIME, dtSave)  < DT) THEN

            IF(TgFlowFlg) THEN
                CALL WRT_AVERAGE_VARS_nonXperiodic_io
                IF(DMOD(PhyTIME, dtAveView)  <  DT) THEN
                    CALL WRT_AVERAGE_PPED_nonXperiodic_io
                END IF
            ELSE
                CALL WRT_AVERAGE_VARS_Xperiodic_io
                IF(DMOD(PhyTIME, dtAveView)  <  DT) THEN
                    CALL WRT_AVERAGE_PPED_Xperiodic_io
                    IF(iPPSpectra == 1) CALL PP_SPECOSPEC
                END IF
            END IF
        END IF

    END IF

    !=============== Screen printing for monitoring the codE ======================================
    IF ( (ITERG == 1) .OR. MOD(ITERG, iterMonitor) == 0 )  THEN
        !======================== MONITORING THE CALCULATION RUNNING PROPERLY ==
        CALL DIVGCK_io

        IF(TgFlowFlg) THEN

            IF(PhyTIME <= tRunAve1) THEN
                CALL PP_MEAN_Z_FLOW_nonXperiodic_io
                IF(iThermoDynamics == 1) &
                CALL PP_MEAN_Z_THEML_nonXperiodic_io
            END IF
            CALL PP_MONITOR_nonXperiodic_io

        ELSE

            IF(PhyTIME <= tRunAve1) THEN
                CALL PP_MEAN_ZX_FLOW_Xperiodic_io
            END IF
            CALL PP_MONITOR_Xperiodic_io

        END IF

    END IF

    RETURN
END SUBROUTINE

!*******************************************************************
