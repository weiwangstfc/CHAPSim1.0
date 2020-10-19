!**********************************************************************************************************************************
!> @brief
!>        to postprocess all given instantanous results
!> @details
!> SUBROUTINE: POSTPROCESS_INTEGRAL_INSTANS (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE POSTPROCESS_INTEGRAL_INSTANS
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE


    INTEGER(4) :: N
    REAL(WP) :: TI

    IF(iPostProcess /= 2) RETURN

    IF(MYID == 0) CALL CHKHDL('IO: rE -postprocessing all instantanous flows', MYID)

    !===========FOR TG PART ==========================
    IF(TgFlowFlg) THEN
        !== initialize the averaged valueS ======
        CALL PP_TMEAN_INI_ZERO_TG

        !== Read in all data=======
        DO N = 1, pp_instn_sz
            ti = pp_instn_tim(N)
            CALL ReStart_INSTANT_VARS_TG(ti)
            CALL PP_MEAN_ZX_FLOW_TG
            CALL PP_MEAN_T_TG(T_Summing_average)
        END DO

        !====WRITE data ouT ===
        CALL WRT_AVERAGE_VARS_TG
        CALL WRT_AVERAGE_PPED_TG
        IF(iPPSpectra == 1) CALL PP_SPECOSPEC

    END IF

    !==========FOR main PART ========================

    IF(IoFlowFlg) THEN
        !== initialize the averaged valueS ======
        CALL PP_TMEAN_INI_ZERO_io

        IF(TgFlowFlg) THEN
            ! to add....

        ELSE  !===periodic main io domain =======================
            !== Read in all data=======
            DO N = 1, pp_instn_sz
                tI = pp_instn_tim(N)
                CALL ReStart_INSTANT_VARS_io(ti)
                CALL PP_MEAN_ZX_FLOW_Xperiodic_io
                IF(iPPSpectra == 1) CALL PP_SPECOSPEC
                CALL PP_MEAN_T_Xperiodic_io(T_Summing_average)
            END DO

            PhyTIME = PhyTIME_io
            ITERG0  = ITERG0_io

            !====WRITE data ouT ===
            CALL WRT_AVERAGE_VARS_Xperiodic_io
            CALL WRT_AVERAGE_PPED_Xperiodic_io
            !IF(iPPSpectra == 1) CALL PP_SPECOSPEC

        END IF

    END IF


    RETURN
END SUBROUTINE
