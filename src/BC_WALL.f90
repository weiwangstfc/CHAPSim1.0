!**********************************************************************************************************************************
!> @brief
!>        The Wall B.C.
!> @details
!> SUBROUTINE: BC_PIPE_CENTRE_Q (in MYID = all)
!> SUBROUTINE: BC_PIPE_CENTRE_PR (in MYID = all)
!> SUBROUTINE: BC_PIPE_CENTRE_DPH (in MYID = all)
!> SUBROUTINE: BC_PIPE_CENTRE_G (in MYID = all)
!> SUBROUTINE: BC_PIPE_CENTRE_THERMAL (in MYID = all)
!> SUBROUTINE: BC_WALL_Q (in MYID = all)
!> SUBROUTINE: BC_WALL_PR (in MYID = all)
!> SUBROUTINE: BC_WALL_DPH (in MYID = all)
!> SUBROUTINE: BC_WALL_G (in MYID = all)
!> SUBROUTINE: BC_WALL_IsoTHERMAL (in MYID = all)
!> SUBROUTINE: BC_WALL_Isoflux (in MYID = all)
!> SUBROUTINE: BC_WALL_THERMAL (in MYID = all)
!> SUBROUTINE: BC_WALL_Q_tg (in MYID = all)
!> SUBROUTINE: BC_WALL_Q_io(in MYID = all)
!> SUBROUTINE: BC_WALL_G_tg (in MYID = all)
!> SUBROUTINE: BC_WALL_G_io(in MYID = all)
!> SUBROUTINE: BC_WALL_PR_tg (in MYID = all)
!> SUBROUTINE: BC_WALL_PR_io(in MYID = all)
!> SUBROUTINE: BC_WALL_DPH_tg (in MYID = all)
!> SUBROUTINE: BC_WALL_DPH_io(in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 12 / 2013- Initial Version, by Wei Wang
!**********************************************************************************************************************************
SUBROUTINE BC_PIPE_CENTRE_Q(iDomain_tmp)
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iDomain_tmp
    INTEGER(4) :: I, K, KS

    IF(iCase /= IPIPEC) RETURN
    IF(MYID /= 0) RETURN

    IF(iDomain_tmp == ITG) THEN

        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                KS = KSYM(K)
                Q_tg(I, 0, K, 1) = Q_tg(I, 1, KS, 1)
                Q_tg(I, 0, K, 3) = Q_tg(I, 1, KS, 3)
                Q_tg(I, 0, K, 2) = Q_tg(I, 2, KS, 2)
                Q_tg(I, 1, K, 2) = 0.0_WP
            END DO
        END DO

    ELSE IF(iDomain_tmp == IIO) THEN

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                KS = KSYM(K)
                Q_io(I, 0, K, 1) = Q_io(I, 1, KS, 1)
                Q_io(I, 0, K, 3) = Q_io(I, 1, KS, 3)
                Q_io(I, 0, K, 2) = Q_io(I, 2, KS, 2)
                Q_io(I, 1, K, 2) = 0.0_WP
            END DO
        END DO

    ELSE

    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_PIPE_CENTRE_PR(iDomain_tmp)
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iDomain_tmp
    INTEGER(4) :: I, K, KS

    IF(iCase /= IPIPEC) RETURN
    IF(MYID /= 0) RETURN

    IF(iDomain_tmp == ITG) THEN

        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                KS = KSYM(K)
                PR_TG(I, 0, K) = PR_TG(I, 1, KS)
            END DO
        END DO

    ELSE IF(iDomain_tmp == IIO) THEN

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                KS = KSYM(K)
                PR_io(I, 0, K) = PR_io(I, 1, KS)
            END DO
        END DO

    ELSE

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_PIPE_CENTRE_DPH(iDomain_tmp)
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iDomain_tmp
    INTEGER(4) :: I, K, KS

    IF(iCase /= IPIPEC) RETURN
    IF(MYID /= 0) RETURN

    IF(iDomain_tmp == ITG) THEN

        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                KS = KSYM(K)
                DPH_TG(I, 0, K) = DPH_TG(I, 1, KS)
            END DO
        END DO

    ELSE IF(iDomain_tmp == IIO) THEN

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                KS = KSYM(K)
                DPH_io(I, 0, K) = DPH_io(I, 1, KS)
            END DO
        END DO

    ELSE

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_PIPE_CENTRE_G(iDomain_tmp)
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iDomain_tmp
    INTEGER(4) :: I, K, KS

    IF(iCase /= IPIPEC) RETURN
    IF(iDomain_tmp /= IIO) RETURN
    IF(MYID /= 0) RETURN

    DO I = NCL1S, NCL1E
        DO K = 1, NCL3
            KS = KSYM(K)
            G_io(I, 0, K, 1) = G_io(I, 1, KS, 1)
            G_io(I, 0, K, 3) = G_io(I, 1, KS, 3)
            G_io(I, 0, K, 2) = G_io(I, 2, KS, 2)
            G_io(I, 1, K, 2) = 0.0_WP
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_PIPE_CENTRE_THERMAL
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: I, K, KS

    IF(iCase /= IPIPEC) RETURN
    IF(MYID /= 0) RETURN

    DO I = NCL1S, NCL1E
        DO K = 1, NCL3
            KS = KSYM(K)
            ENTHALPY  (I, 0, K) = ENTHALPY  (I, 1, KS)
            DENSITY   (I, 0, K) = DENSITY   (I, 1, KS)
            TEMPERATURE(I, 0, K) = TEMPERATURE(I, 1, KS)
            Viscousity(I, 0, K) = Viscousity(I, 1, KS)
            THERMCONDT(I, 0, K) = THERMCONDT(I, 1, KS)
            HEATCAP   (I, 0, K) = HEATCAP   (I, 1, KS)
            DH      (I, 0, K) = ENTHALPY  (I, 1, KS) *  &
            DENSITY   (I, 1, KS)
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_Q(N, iDomain_tmp)

    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE


    INTEGER(4), INTENT(IN) :: N
    INTEGER(4), INTENT(IN) :: iDomain_tmp

    INTEGER(4) :: I, K

    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF(iDomain_tmp == ITG) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    Q_tg(I, 0, K, 1) = 0.0_WP!-Q_tg(I, 1, K, N)
                    Q_tg(I, 0, K, 3) = 0.0_WP!-Q_tg(I, 3, K, N)
                    Q_tg(I, 1, K, 2) = 0.0_WP
                    Q_tg(I, 0, K, 2) = 0.0_WP
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    Q_tg(I, N2DO(MYID) + 1, K, 1) = 0.0_WP
                    Q_tg(I, N2DO(MYID) + 1, K, 3) = 0.0_WP
                    Q_tg(I, N2DO(MYID) + 1, K, 2) = 0.0_WP
                END DO
            END DO

        END IF

    ELSE IF(iDomain_tmp == IIO) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    Q_io(I, 0, K, 1) = 0.0_WP!-Q_io(I, 1, K, N)
                    Q_io(I, 0, K, 3) = 0.0_WP!-Q_io(I, 3, K, N)
                    Q_io(I, 1, K, 2) = 0.0_WP
                    Q_io(I, 0, K, 2) = 0.0_WP
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    Q_io(I, N2DO(MYID) + 1, K, 1) = 0.0_WP
                    Q_io(I, N2DO(MYID) + 1, K, 3) = 0.0_WP
                    Q_io(I, N2DO(MYID) + 1, K, 2) = 0.0_WP
                END DO
            END DO

        END IF

    ELSE

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_PR(N, iDomain_tmp)

    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE


    INTEGER(4), INTENT(IN) :: N
    INTEGER(4), INTENT(IN) :: iDomain_tmp

    INTEGER(4) :: I, K

    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF(iDomain_tmp == ITG) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    PR_tg(I, 0, K) = PR_tg(I, 1, K)
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    PR_tg(I, N2DO(MYID) + 1, K) = PR_tg(I, N2DO(MYID), K)
                END DO
            END DO

        END IF

    ELSE IF(iDomain_tmp == IIO) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    PR_io(I, 0, K) = PR_io(I, 1, K)
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    PR_io(I, N2DO(MYID) + 1, K) = PR_io(I, N2DO(MYID), K)
                END DO
            END DO

        END IF

    ELSE

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_DPH(N, iDomain_tmp)

    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE


    INTEGER(4), INTENT(IN) :: N
    INTEGER(4), INTENT(IN) :: iDomain_tmp

    INTEGER(4) :: I, K

    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF(iDomain_tmp == ITG) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    DPH_tg(I, 0, K) = DPH_tg(I, 1, K)
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    DPH_tg(I, N2DO(MYID) + 1, K) = DPH_tg(I, N2DO(MYID), K)
                END DO
            END DO

        END IF

    ELSE IF(iDomain_tmp == IIO) THEN

        IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    DPH_io(I, 0, K) = DPH_io(I, 1, K)
                END DO
            END DO

        END IF

        IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    DPH_io(I, N2DO(MYID) + 1, K) = DPH_io(I, N2DO(MYID), K)
                END DO
            END DO

        END IF

    ELSE

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_G(N, iDomain_tmp)

    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE


    INTEGER(4), INTENT(IN) :: N
    INTEGER(4), INTENT(IN) :: iDomain_tmp

    INTEGER(4) :: I, K

    IF(iDomain_tmp /= IIO) RETURN
    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF (MYID == 0 .AND. N == iBotWall) THEN !=========bottom walL ==================

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                G_io(I, 0, K, 1) = 0.0_WP!-Q_io(I, 1, K, N)
                G_io(I, 0, K, 3) = 0.0_WP!-Q_io(I, 3, K, N)
                G_io(I, 1, K, 2) = 0.0_WP
                G_io(I, 0, K, 2) = 0.0_WP
            END DO
        END DO

    END IF

    IF (MYID == NPSLV .AND. N == iTopWall) THEN !========= Top walL ==================

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                G_io(I, N2DO(MYID) + 1, K, 1) = 0.0_WP
                G_io(I, N2DO(MYID) + 1, K, 3) = 0.0_WP
                G_io(I, N2DO(MYID) + 1, K, 2) = 0.0_WP
            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_IsoTHERMAL(N)
    !>      @note
    !>      All thermal variables at wall are located on the wall,
    !>      Rather than the ghost cells in exterior
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N
    INTEGER(4) :: I, J, K

    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    !WRITE(*, *) 'MYID', MYID, N !test
    IF(MYID == 0 .AND. N == iBotWall) THEN
        !WRITE(*, *) 'MYIDsss',  MYID, N !test
        J = 0
        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                TEMPERATURE(I, J, K) = T_WAL_GV(I, N)
                ENTHALPY  (I, J, K) = H_WAL_GV(I, N)
                DENSITY   (I, J, K) = D_WAL_GV(I, N)
                Viscousity(I, J, K) = M_WAL_GV(I, N)
                THERMCONDT(I, J, K) = K_WAL_GV(I, N)
                DH      (I, J, K) = ENTHALPY(I, J, K) * DENSITY(I, J, K)
            END DO
        END DO
        !WRITE(*, *) 'N, J, T', N, J, T_WAL_GV(1, N), T_WAL_GV(NCL1_io / 2, N)!test
    END IF

    IF(MYID == NPSLV .AND. N == iTopWall) THEN
        !WRITE(*, *) 'MYIDsss',  MYID, N !test
        J = N2DO(MYID) + 1
        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                TEMPERATURE(I, J, K) = T_WAL_GV(I, N)
                ENTHALPY  (I, J, K) = H_WAL_GV(I, N)
                DENSITY   (I, J, K) = D_WAL_GV(I, N)
                Viscousity(I, J, K) = M_WAL_GV(I, N)
                THERMCONDT(I, J, K) = K_WAL_GV(I, N)
                DH      (I, J, K) = ENTHALPY(I, J, K) * DENSITY(I, J, K)
            END DO
        END DO
        !WRITE(*, *) 'N, J, T', N, J, T_WAL_GV(1, N), T_WAL_GV(NCL1_io / 2, N)!test
    END IF



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_Isoflux(N, IREGION)
    USE flow_info
    USE MESH_INFO
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IREGION, N
    INTEGER(4) :: J, K, I1, I2, I!, JS, KS
    REAL(WP) :: H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp


    IF(N == iBotWall .AND. iCase == IPIPEC) RETURN

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF(IREGION == IOULET) THEN  ! outlet
        I1 = NCL1_io + 1
        I2 = NCL1_io + 1
    ELSE IF(IREGION == IINLET) THEN !INLET
        I1 = 0
        I2 = 0
    ELSE IF(IREGION == IMAIND) THEN
        I1 = 1
        I2 = NCL1_io
    ELSE IF(IREGION == IALLDM) THEN
        I1 = NCL1S
        I2 = NCL1E
    ELSE
        I1 = 1
        I2 = NCL1_io
    END IF


    IF(MYID == 0 .OR. MYID == NPSLV) THEN

        IF(MYID == 0 .AND. N == iBotWall)     THEN
            J = 0
            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    DH(I, J, K) = 2.0_WP * DH(I, J + 1, K) &
                    - YCL2ND_WFF(J + 2) * DH(I, J + 2, K) &
                    - YCL2ND_WFB(J + 2) * DH(I, J + 1, K)
                    CALL THERM_PROP_UPDATE_FROM_DH(DH(I, J, K), H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
                    ENTHALPY(I, J, K) = H_tmp
                    TEMPERATURE(I, J, K) = T_tmp
                    DENSITY(I, J, K) = D_tmp
                    Viscousity(I, J, K) = M_tmp
                    THERMCONDT(I, J, K) = K_tmp
                    HEATCAP(I, J, K) = Cp_tmp
                END DO
            END DO

        END IF

        IF(MYID == NPSLV .AND. N == iTopWall) THEN
            J = N2DO(MYID) + 1
            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    DH(I, J, K) = 2.0_WP * DH(I, J - 1, K) &
                    - YCL2ND_WFF(J - 1) * DH(I, J - 1, K) &
                    - YCL2ND_WFB(J - 1) * DH(I, J -2, K)
                    CALL THERM_PROP_UPDATE_FROM_DH(DH(I, J, K), H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
                    ENTHALPY(I, J, K) = H_tmp
                    TEMPERATURE(I, J, K) = T_tmp
                    DENSITY(I, J, K) = D_tmp
                    Viscousity(I, J, K) = M_tmp
                    THERMCONDT(I, J, K) = K_tmp
                    HEATCAP(I, J, K) = Cp_tmp
                END DO
            END DO

        END IF
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_THERMAL(IREGION)
    USE init_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IREGION

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux)    CALL BC_WALL_Isoflux   (iTopWall, IREGION)
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature) CALL BC_WALL_IsoTHERMAL(iTopWall)

    IF(iCase /= Ipipec) THEN
        IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux)    CALL BC_WALL_Isoflux   (iBotWall, IREGION)
        IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature) CALL BC_WALL_IsoTHERMAL(iBotWall)
    ELSE
        CALL BC_PIPE_CENTRE_THERMAL
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_Q_tg
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_Q(iTopWall, ITG)
    IF(iCase /= Ipipec) CALL BC_WALL_Q(iBotWall, ITG)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_Q(ITG)
    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_Q_io
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_Q(iTopWall, IIO)
    IF(iCase /= Ipipec) CALL BC_WALL_Q(iBotWall, IIO)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_Q(IIO)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_G_io
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_G(iTopWall, IIO)
    IF(iCase /= Ipipec) CALL BC_WALL_G(iBotWall, IIO)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_G(IIO)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_PR_TG
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_PR(iTopWall, ITG)
    IF(iCase /= Ipipec) CALL BC_WALL_PR(iBotWall, ITG)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_PR(ITG)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_PR_io
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_PR(iTopWall, IIO)
    IF(iCase /= Ipipec) CALL BC_WALL_PR(iBotWall, IIO)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_PR(IIO)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_DPH_TG
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_DPH(iTopWall, ITG)
    IF(iCase /= Ipipec) CALL BC_WALL_DPH(iBotWall, ITG)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_DPH(ITG)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BC_WALL_DPH_io
    USE init_info
    IMPLICIT NONE

    IF(BCY(1) /= IBCWALL) RETURN
    IF(BCY(2) /= IBCWALL) RETURN

    CALL BC_WALL_DPH(iTopWall, IIO)
    IF(iCase /= Ipipec) CALL BC_WALL_DPH(iBotWall, IIO)
    IF(iCase == Ipipec) CALL BC_PIPE_CENTRE_DPH(IIO)

    RETURN
END SUBROUTINE
