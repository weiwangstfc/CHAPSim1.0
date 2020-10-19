!**********************************************************************************************************************************
!> @brief
!>
!> @details
!> SUBROUTINE: SOLVERRK3_ENG_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE SOLVERRK3_ENG_io(NS)
    USE init_info
    USE mesh_info
    USE thermal_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4) :: I, J, K
    REAL(WP) :: COE1, COE0


    Viscousity0(:, :, :) = Viscousity(:, :, :) ! for time stagger in momentum equation.

    !DENSITY1(:, :, :) = DENSITY0(:, :, :)
    DENSITY0(:, :, :) = DENSITY(:, :, :)
    DH0(:, :, :)    = DH(:, :, :)

    !======== STEP 2: main DOMAIN  ENERGY EQUATION ==\DH =================
    COE1 = TGAM(NS) * DT
    COE0 = TROH(NS) * DT

    CALL RHS_ENERGY_EXPLICIT

    DO J = 1, N2DO(MYID)
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                DH(I, J, K) = DH(I, J, K) + COE1 * RHS_ENERGY (I, J, K) + COE0 * RHS_ENERGY0(I, J, K)
                IF(iThermoProperty == search_table) THEN
                    DH(I, J, K) = DMIN1(DHmax - REALMIN, DH(I, J, K)) ! limiter!
                    DH(I, J, K) = DMAX1(DHmin + REALMIN, DH(I, J, K)) ! limiter!
                END IF

                RHS_ENERGY0(I, J, K) = RHS_ENERGY(I, J, K)
            END DO
        END DO
    END DO

    !==================checK =============
    !IF(NS == 3) CALL CHK_EnegConsv_io
    !IF(MYID == 0) WRITE(*, *) 'eneg conservation(inter): NS ', NS, CHK_ENEG_CONSV0

    !======== STEP 3:  PROVIsoNAL ESTIMATE FOR H ================
    !CALL ENTHALPY_UPDATE(IMAIND) ! DH/ Rho

    !======== STEP 4:  DENSITY UPDATE ===========================
    !CALL DENSITY_UPDATE(IMAIND)  ! h->T ->rho

    !======== STEP 5:  ENTHALPY UPDATE ==========================
    !CALL ENTHALPY_UPDATE(IMAIND) ! DH/ Rho

    !======== STEP 6:  THERMAL_PROPERTY_UPDATE ==================
    CALL THERM_PROP_UPDATE(IMAIND) ! no wall b.c.

    !======== STEP 7: SET UP THERMAL PROPERTIES IN INTERFACES ============
    CALL INTFC_MFD_THERMAL_io

    !======== STEP 8:  b.c. ENERGY EQUATION ==\DH =================
    IF(TgFlowFlg) THEN
        CALL BC_TINLET_THERML
        CALL BC_COUTLET_ENEG_RK3(NS)
    END IF

    CALL BC_WALL_THERMAL(IMAIND)

    !CALL DEBUG_WRT_THERMAL
    !==================checK =============
    !CALL CHK_EnegConsv_io
    !IF(MYID == 0) WRITE(*, *) 'eneg conservation(final): NS ', NS, CHK_ENEG_CONSV0

    RETURN
END SUBROUTINE
