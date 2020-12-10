!**********************************************************************************************************************************
!> @brief
!>        to calculate the initial fields
!> @details
!> SUBROUTINE: FLOWStart (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE FLOWStart
    USE init_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    !======== iNITIAL FLOW FIELD and CONDITIONS ===FOR TURBULENCE GENERATOR ==========================================
    IF(TgFlowFlg) THEN

        IF (iIniField_tg == IniField_random) THEN
            IF (MYID == 0) &
            CALL CHKHDL('14. TG: Flow initialization from random velocity field', MYID)
            CALL random_FL_FLD_TG
            CALL CALC_INITIALIZATION_tg
        ELSE IF (iIniField_tg == IniField_extrapolation) THEN
            IF (MYID == 0) &
            CALL CHKHDL('14. TG: Flow initialization from previous coarse mesh to a fine mesh using extrapolation.', MYID)
            CALL INITIAL_INTERP_tg
            CALL WRT_INSTANT_VARS_TG
        ELSE IF (iIniField_tg == IniField_reStart) THEN
            IF (MYID == 0)  &
            CALL CHKHDL('14. TG: Flow initialization from last step using the same mesh.', MYID)
            CALL ReStart_INSTANT_VARS_TG(TimeReStart_tg)
        ENDIF

        CALL INTFC_VARS1(1, NCL1_tg, 1, NCL1_tg,PR_tg)
        CALL INTFC_VARS3(1, NCL1_tg, 1, NCL1_tg, Q_tg)
        CALL BC_WALL_Q_tg
        CALL BC_WALL_PR_TG
        CALL DIVGCK_tg
        IF(MYID == 0) THEN
            CALL CHKHDL   ('(6-1) TG: Max divergence of the velocity in final initial flow field', MYID)
            CALL CHKRLHDL('          Div(Velocity) = ', MYID, MAXDIVGV_tg(2))
            IF(MAXDIVGV_tg(2) > 1.0E-6_WP) CALL ERRHDL('Large divergence of the velocity field.', MYID)
        END IF
        MAXDIVGV_tg(2) = MAXDIVGV_tg(1)
        CALL PP_TMEAN_INI_TG
        tRunAve1 = DMAX1(tRunAve1, tRunAve_Reset)

    END IF


    !======== iNITIAL FLOW FIELD and CONDITIONS ===FOR THE main DOMAIN ============================================
    IF(IoFlowFlg) THEN

        !================ initial flow field===========================
        IF (iIniField_io == IniField_random) THEN  !
            IF (MYID == 0) CALL CHKHDL('14. IO: Flow initialization from random velocity field', MYID)

            CALL random_FL_THEML_FLD_io
            CALL CALC_INITIALIZATION_io
        ELSE IF (iIniField_io == IniField_extrapolation) THEN  !
            !==== NREAD= 1 : from coarse mesh to fine mesh using extrapolatioN ======
            IF (MYID == 0) &
            CALL CHKHDL('14.IO: Flow initialization from previous coarse mesh to a fine mesh using extrapolation.', MYID)
            CALL INITIAL_INTERP_io
            CALL WRT_INSTANT_VARS_io

        ELSE IF (iIniField_io == IniField_reStart) THEN
            !==== NREAD= 2 reStarting from the same mesH ======
            IF (MYID == 0) THEN
                CALL CHKHDL('14. IO: Flow initialization from last step using the same mesh.', MYID)
            END IF

            CALL ReStart_INSTANT_VARS_io(TimeReStart_io)

            CALL DIVGCK_io
            IF(MYID == 0) THEN
                CALL CHKHDL   ('IO: Max divergence of initial flow field. Main Domain; Inlet; Outlet.', MYID)
                CALL CHK3RLHDL('    Div(Velocity) = ', MYID, MAXDIVGV_io(1), MAXDIVGV_io(2), MAXDIVGV_io(3))
            END IF
            IF(iWeightedPre == 1) THEN
                PR0_io(:, :, :, 1) = PR_io(:, :, :)
                PR0_io(:, :, :, 2) = PR_io(:, :, :)
            END IF
        ENDIF

        !=========== After initial flow field==============================
        CALL PP_TMEAN_INI_io
        tRunAve1 = DMAX1(tRunAve1, tRunAve_Reset)

        !==========below IS only a repeat for securitY ====================
        IF(iThermoDynamics == 1) THEN
            CALL BC_WALL_THERMAL(IALLDM)
            IF(TgFlowFlg) CALL BC_TINLET_THERML
            CALL THERM_PROP_UPDATE(IALLDM)
            IF(TgFlowFlg) THEN
                CALL INTFC_OUL_THERMAL_io
                CALL INTFC_INL_THERMAL_io
            END IF
            CALL INTFC_MFD_THERMAL_io
        END IF

        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
        CALL BC_WALL_Q_io
        CALL BC_WALL_G_io
        CALL BC_WALL_PR_io

    END IF

    CALL CALL_TEC360

    IF(iPostProcess == 1) THEN
        CALL MPI_BARRIER(ICOMM, IERROR)
        IF(MYID == 0) CALL CHKHDL('<===Only postprocessed results, now the code stops...==>', MYID)
        STOP "Finished postprocessing=1"
    END IF



    RETURN
END SUBROUTINE
