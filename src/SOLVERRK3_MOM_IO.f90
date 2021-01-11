!**********************************************************************************************************************************
!> @brief
!>
!> @details
!> SUBROUTINE: SOLVERRK3_MOM_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE SOLVERRK3_MOM_io(NS)
    USE cparam
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4) :: IDR

    !INTEGER(4) :: JC, JJ, KC, IC

    !============for time stagereed in energy equatioN =====================
    G0_io(:, :, :, :) = G_io(:, :, :, :)
    !======= STEP6: CALCULATE VELOCITY =================
    CALL VELOCITY_CALC_io

    !======= STEP1: CALCULATE OUTLET B.C.=================================
    IF(TgFlowFlg) THEN
        CALL BC_TINLET_FLOW
        CALL BC_COUTLET_MOM_RK3(NS)  ! OUTLET

        IF (iVisScheme == VisImplicit) CALL BC_CBC_TDMA(NS) !!to DO...
    END IF


    !======= STEP1: CALCULATE CONVECTION TERMS ===============================
    CALL CONVECTION_X_io
    CALL CONVECTION_Y_io
    CALL CONVECTION_Z_io

    !======= STEP2: CALCULATE THE WHOLE RHS FOR X, Y,Z ==========================
    IF(iVisScheme == VisExplicit) THEN
        CALL DIVG_U_io

        IDR = 1
        CALL VISCOUS_ALL_EXPLT_X_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)

        IDR = 2
        CALL VISCOUS_ALL_EXPLT_Y_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)

        IDR = 3
        CALL VISCOUS_ALL_EXPLT_Z_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)

        CALL MASSFLUX_CALC_io

    ELSE IF (iVisScheme == VisImplicit) THEN
        IDR = 1
        CALL VISCOUS_PAR_EXPLT_X_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)
        CALL MOMFA_io(NS, IDR)

        IDR = 2
        CALL VISCOUS_PAR_EXPLT_Y_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)
        CALL MOMFA_io(NS, IDR)

        IDR = 3
        CALL VISCOUS_PAR_EXPLT_Z_io
        CALL RHS_MOM_EXPLICIT_io(NS, IDR)
        CALL MOMFA_io(NS, IDR)

    ELSE
        CALL ERRHDL(' No Such Scheme for vIScous term', MYID)
    END IF

    CALL INTFC_VARS3(1, NCL1_io, NCL1S, NCL1E, G_io)
    CALL BC_WALL_G_io

    !======= STEP4: CONSTRUCTING AND SOLVING POISSION EQ.=====================
    CALL DIVG_io(NS)

    !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io, 1, N2DO(MYID), 'divg') !test

    IF(TgFlowFlg) THEN
        CALL FISHPACK_POIS3D_SIMPLE
    ELSE
      if(is_FFT_FFT99) CALL FFT99_POIS3D_periodicxz(IIO) !Method One
      if(is_FFT_FISHPACK) CALL FISHPACK_POIS3D_SIMPLE ! Method Two,  good
    END IF

    CALL INTFC_VARS1(1, NCL1_io, NCL1S, NCL1E, DPH_io)
    CALL BC_WALL_DPH_io

    !CALL CHECK_FFT_SOLVER!test
    !CALL DEBUG_WRT_LOCAL(DPH_io, 0, N2DO(MYID) + 1, 'dphi') !test

    !======= STEP5: CALCULATE pressure =================
    CALL PRCALC_io(NS)

    !======= STEP6: CALCULATE mass fluX =================
    CALL MASSFLUX_UPDATE_io(NS)

    !======= STEP6: CALCULATE VELOCITY =================
    CALL VELOCITY_CALC_io

    !==========checK =============
    !CALL CHK_MassConsv_io
    !IF(MYID == 0) WRITE(*, *) 'mass conservation(final): NS ', NS, CHK_MASS_CONSV0
    !IF(MYID == 0)     WRITE(*, *) 'pressure1', &
    !PR_io(4, 0, 4),         PR_io(4, 1, 4),           PR_io(4, N2DO(MYID), 4), PR_io(4, N2DO(MYID) + 1, 4)
    !IF(MYID == NPSLV) WRITE(*, *) 'pressure2', &
    !PR_io(4, N2DO(MYID), 4), PR_io(4, N2DO(MYID) + 1, 4), PR_io(4, 0, 4),         PR_io(4, 1, 4)

    RETURN
END SUBROUTINE
