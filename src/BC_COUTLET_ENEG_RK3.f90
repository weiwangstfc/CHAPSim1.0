!**********************************************************************************************************************************
!> @brief
!>        convective outlet b.c. for the energy equation
!> @details
!> SUBROUTINE: BC_COUTLET_ENEG_RK3 (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 04/ 2014- Initial version, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE BC_COUTLET_ENEG_RK3(NS)
    USE mesh_info
    USE flow_info
    USE init_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS

    INTEGER(4) :: J, K
    REAL(WP) :: BC_RHS
    REAL(WP) :: BC_CONV


    IF(.NOT. TgFlowFlg) RETURN

    CALL CONVCTION_OUTLET_U

    DO K = 1, NCL3
        DO J = 1, N2DO(MYID)
            !U_OUTLET =     !Q_io(NCL1_io + 1, J, K, 1) whICh may be negative.
            BC_CONV = - DXI * U_OUTLET * &
            ( DH(NCL1_io + 1, J, K) - DH(NCL1_io, J, K) )
            BC_RHS = TGAM(NS) * BC_CONV + TROH(NS) * BC_CONV0_ENG(J, K)
            BC_CONV0_ENG(J, K) = BC_CONV

            DH(NCL1_io + 1, J, K) = DH(NCL1_io + 1, J, K) + DT * BC_RHS
        END DO
    END DO

    ! limiter!
    IF(iThermoProperty == search_table) THEN
        DO K = 1, NCL3
            DO J = 1, N2DO(MYID)
                DH(NCL1_io + 1, J, K) = DMIN1(DHmax - REALMIN, DH(NCL1_io + 1, J, K))
                DH(NCL1_io + 1, J, K) = DMAX1(DHmin + REALMIN, DH(NCL1_io + 1, J, K))
            END DO
        END DO
    END IF

    !CALL ENTHALPY_UPDATE(IOULET)
    !CALL DENSITY_UPDATE(IOULET)
    !CALL ENTHALPY_UPDATE(IOULET)

    CALL THERM_PROP_UPDATE(IOULET)
    CALL INTFC_OUL_THERMAL_io
    CALL BC_WALL_THERMAL(IOULET)

    RETURN
END SUBROUTINE
