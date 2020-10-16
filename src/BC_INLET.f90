!**********************************************************************************************************************************
!> @brief
!>        inlet condition for io developing flow
!> @details
!> SUBROUTINE: BC_TINLET_FLOW (in MYID = all)
!> SUBROUTINE: BC_TINLET_THERML (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 04/ 2014- Initial version, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE BC_TINLET_FLOW
    !>     @note
    !>     Inlet bc will USE the data from the adjacent cells in the turbulence
    !>     generator.
    !>     For TG, energy equation may be solved with no heat flux fed in, therefore
    !>     It IS still like incompressible flow with constant thermal properties.
    !>     Thus, the inlet thermal properties USE constant values.

    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4) :: J, K

    IF(.NOT. TgFlowFlg) RETURN

    !===================FOR FLOW FIELD=========================
    DO J = 0, N2DO(MYID) + 1
        DO K = 1, NCL3

            !=====U=========
            Q_io(0, J, K, 1) = Q_tg(NCL1_tg, J, K, 1)
            Q_io(1, J, K, 1) = Q_tg(1,       J, K, 1)

            G_io(0, J, K, 1) = Q_io(0, J, K, 1) * DENSITY0(0, J, K)
            !G_io(1, J, K, 1) = Q_io(1, J, K, 1) * ( DENSITY(0, J, K) + DENSITY(1, J, K) ) * XND2CL
            G_io(1, J, K, 1) = Q_io(1, J, K, 1) * DENSITY0(0, J, K) ! check for RK iteration

            !===== V =========
            Q_io(0, J, K, 2) = Q_tg(NCL1_tg, J, K, 2)
            G_io(0, J, K, 2) = Q_io(0, J, K, 2) * DENSITY0(0, J, K)

            !=====W=========
            Q_io(0, J, K, 3) = Q_tg(NCL1_tg, J, K, 3)
            G_io(0, J, K, 3) = Q_io(0, J, K, 3) * DENSITY0(0, J, K)

            !==========P ===========
            PR_io(0, J, K) = PR_tg(NCL1_tg, J, K)
        END DO
    END DO
    DPH_io(0, :, :) = 0.0_WP

    !== As J = 0~n2DO+ 1, the intfc IS not necessary =====
    !CALL INTFC_INL_G_io
    !CALL INTFC_INL_Q_io
    !CALL INTFC_INL_P_io

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_TINLET_THERML
    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4) :: J, K

    IF(.NOT. TgFlowFlg) RETURN

    DO J = 0, N2DO(MYID) + 1, 1
        DO K = 1, NCL3
            TEMPERATURE(0, J, K) = T_inlet
            THERMCONDT(0, J, K) = K_inlet
            DENSITY   (0, J, K) = D_inlet
            !DENSITY0  (0, J, K) = 1.0_WP
            Viscousity(0, J, K) = M_inlet
            HEATCAP   (0, J, K) = CP_inlet

            DH      (0, J, K) = DH_inlet
            ENTHALPY  (0, J, K) = H_inlet

        END DO
    END DO

    !        !================FOR THERMAL FIELD============================
    !        !====ONLY TO BE CALLED BEFORE main SOLVER ONCE =======
    !        IF(iThermalWallType == BC_Fixed_Heat_Flux) THEN
    !            DO J = 0, N2DO(MYID) + 1, 1
    !                DO K = 1, NCL3
    !                    TEMPRATURE(0, J, K) = 1.0_WP
    !                    THERMCONDT(0, J, K) = 1.0_WP
    !                    DENSITY   (0, J, K) = 1.0_WP
    !                    !DENSITY0  (0, J, K) = 1.0_WP
    !                    Viscousity(0, J, K) = 1.0_WP

    !                    DH      (0, J, K) = 0.0_WP
    !                    ENTHALPY  (0, J, K) = 0.0_WP

    !                END DO
    !            END DO
    !        END IF

    !        IF(iThermalWallType == BC_Fixed_Temperature) THEN
    !            DO J = 0, N2DO(MYID) + 1, 1
    !                IF(MYID == 0     .AND. J == 0             ) cycle
    !                IF(MYID == NPSLV .AND. J == (N2DO(MYID) + 1)) cycle
    !                DO K = 1, NCL3
    !                    TEMPRATURE(0, J, K) = 2.0_WP * TEMPRATURE(1, J, K) - TEMPRATURE(2, J, K)
    !                    THERMCONDT(0, J, K) = 2.0_WP * THERMCONDT(1, J, K) - THERMCONDT(2, J, K)
    !                    DENSITY   (0, J, K) = 2.0_WP * DENSITY   (1, J, K) - DENSITY   (2, J, K)
    !                    !DENSITY0  (0, J, K) = 2.0_WP * DENSITY0  (1, J, K) - DENSITY0  (2, J, K)
    !                    Viscousity(0, J, K) = 2.0_WP * Viscousity(1, J, K) - Viscousity(2, J, K)
    !                    ENTHALPY  (0, J, K) = 2.0_WP * ENTHALPY  (1, J, K) - ENTHALPY  (2, J, K)
    !                    DH      (0, J, K) = 2.0_WP * DH      (1, J, K) - DH      (2, J, K)
    !                END DO
    !            END DO
    !            CALL INTFC_INL_THERMAL_io
    !        END IF

    RETURN
END SUBROUTINE
