!**********************************************************************************************************************************
!> @brief
!>        Flow field initialization of the tg domAIn
!> @details
!> SUBROUTINE: random_FL_FLD_TG (in MYID = all)
!> SUBROUTINE: CALC_INITIALIZATION_tg (in MYID = all)
!> SUBROUTINE: BULK_VELOCITY_TG (in MYID = all)
!> SUBROUTINE: IniField_FLOW_tg (in MYID = all)
!>             - to calculate initial velocity and pressure flow field
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 12 / 2013- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE random_FL_FLD_TG
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, L
    INTEGER(4) :: J, JJ
    INTEGER(4) :: K
    REAL(WP) :: Umean1_tg


    !>===========Generate scaled random Q and PR = 1 ==========================
    !       Turblence Generator***
    IF(MYID == 0) CALL CHKHDL('(1-1) TG: Generating random velocity ...', MYID)
    CALL IniField_FLOW_tg
    CALL VMAV_tg
    IF(MYID == 0) THEN
        CALL CHKHDL   ('TG: VMV(1:3) in generated random velocity field ', MYID)
        CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_tg(1), VMAX_tg(2), VMAX_tg(3))
        CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_tg(1), VMIN_tg(2), VMIN_tg(3))
    END IF

    !>================ Q = Q-QxzmeaN =========================================
    !        Turblence Generator***
    CALL INTFC_VARS3(1, NCL1_tg, 1, NCL1_tg, Q_tg)
    CALL INITIAL_MEAN_TG
    DO L = 1, 3
        DO J = 1, N2DO(MYID)
            Q_tg(:, J, :,L) = Q_tg(:, J, :,L) - UU(J, L, 1)
        END DO
    END DO
    CALL VMAV_tg
    IF(MYID == 0) THEN
        CALL CHKHDL('(2-1) TG: VMV(1:3) in random velocity field after subtracting meanXZ ', MYID)
        CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_tg(1), VMAX_tg(2), VMAX_tg(3))
        CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_tg(1), VMIN_tg(2), VMIN_tg(3))
    END IF

    !>===========update Q in the incoming flow direction ====================
    !        Turblence Generator***
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                Q_tg(I, J, K, NFLOW) = Q_tg(I, J, K, NFLOW) + Vini(JJ)
            END DO
        END DO
    END DO

    !IF(iCase /= 1) CALL VELO2RVELO_tg

    !>       @note Set Q(I, J, K, 2) = v at Wall b.c. be zero.
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        IF ( (JJ == 1) .OR. (JJ == (NND2)) ) THEN
            Q_tg(:, J, :, 2) = 0.0_WP
        ENDIF
    END DO

    !==============CORRECT THE MASS FLOW RATE = 1 =========================================
    CALL BULK_VELOCITY_TG(Umean1_tg)
    IF(MYID == 0) CALL CHKRLHDL  ('TG: The bulk velocity (original) = ', MYID, Umean1_tg)

    DO J = 1, N2DO(MYID)
        Q_tg(:, J, :, NFLOW) = Q_tg(:, J, :, NFLOW) / Umean1_tg
    END DO

    CALL BULK_VELOCITY_TG(Umean1_tg)
    IF(MYID == 0) CALL CHKRLHDL  ('TG: The bulk velocity (corrected) = ', MYID, Umean1_tg)

    !>============ SCALING DUE TO NOn -INLET Reference ==JunjiE ====================================
    Q_tg(:, :, :, :) = Q_tg(:, :, :, :) * M_inlet / D_inlet
    CALL BULK_VELOCITY_TG(Umean1_tg)
    IF(MYID == 0) CALL CHKRLHDL  ('TG: The bulk velocity (scaled) = ', MYID, Umean1_tg)

    !>============CHECK DIVERGENCE ======================================
    CALL INTFC_VARS3(1, NCL1_tg, 1, NCL1_tg, Q_tg)
    CALL DIVGCK_tg
    IF(MYID == 0) THEN
        CALL CHKHDL  ('(3-1) TG: Max divergence of initial flow field', MYID)
        CALL CHKRLHDL('          Div(Velocity) = ', MYID, MAXDIVGV_tg(2))
    END IF

    CALL CALL_TEC360

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CALC_INITIALIZATION_tg
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: L

    !======Calculate the QP fields =========================================
    !        Turblence Generator***
    IF(MYID == 0) CALL CHKHDL('(4-1) TG: Initial flow field, calc QP...', MYID)

    CALL SOLVERRK3_MOM_tg(0)

    !==============CHECK INFO========================================
    CALL VMAV_tg
    IF(MYID == 0) THEN
        CALL CHKHDL('(5-1) TG: Initial flow field with first calculated velocity, VMAX(1:3)', MYID)
        CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_tg(1), VMAX_tg(2), VMAX_tg(3))
        CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_tg(1), VMIN_tg(2), VMIN_tg(3))
    END IF

    CALL DIVGCK_tg
    IF(MYID == 0) THEN
        CALL CHKHDL  ('(6-1) TG: Max divergence of the velocity in final initial flow field', MYID)
        CALL CHKRLHDL('          Div(Velocity) = ', MYID, MAXDIVGV_tg(2))
        IF(MAXDIVGV_tg(2) > 1.0E-6_WP) CALL ERRHDL('Large divergence of the velocity field.', MYID)
    END IF
    MAXDIVGV_tg(2) = MAXDIVGV_tg(1)

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE BULK_VELOCITY_TG(UUMEAN_WORK)
    USE init_info
    USE mesh_info
    USE flow_info
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    INTEGER(4) :: I
    INTEGER(4) :: J, JJ
    INTEGER(4) :: K
    REAL(WP) :: UUMEAN_WORK, UUMEAN

    UUMEAN = 0.0_WP
    DO J = 1, N2DO(MYID)  !@
        JJ = JCL2G(J)
        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                UUMEAN = UUMEAN + Q_tg(I, J, K, NFLOW) / DYFI(JJ) / RCCI1(JJ)
            END DO
        END DO
    END DO
    UUMEAN = UUMEAN / DZI

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(UUMEAN, UUMEAN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    UUMEAN_WORK = UUMEAN_WORK / (Area_inlet * DBLE(NCL1_tg))

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE IniField_FLOW_tg
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE


    INTEGER(4) :: I
    INTEGER(4) :: J, JJ
    INTEGER(4) :: K
    INTEGER(4) :: RDCOUNT, seed
    !       INTEGER(4) :: IDUM
    REAL(WP) :: VPERQ
    REAL(WP) :: XRDM(3)
    REAL(WP) :: XRD(3)

    Q_tg = 0.0_WP
    PR_tg = 0.0_WP

    IF(iRandomType /= flag_random_sine) THEN
        XRDM = 0.0_WP
        RDCOUNT = 0
        SEED = 0
        !>@note Generate the random flow field u, V,w.
        DO J = 1, N2DO(MYID)           ! only local y_Cell no.
            JJ = JCL2G(J)
            !>@pARam VPERQ = Scaled magnitude of perturbation based on Up (pARabolIC max. velo.)
            !>@note for 1 /4 neAR wall region, perturbation ratio IS decreased by 25%.
            VPERQ = VPERG
            IF((1.0_WP - DABS(YCC(JJ))) < 0.250_WP) VPERQ = SVPERG

            DO I = 1, NCL1_tg
                DO K = 1, NCL3
                    RDCOUNT = I + K + JJ  ! make it independent of processor number.
                    IF(iRandomType == flag_random_real) seed = 0
                    IF(iRandomType == flag_random_fixed) seed = 1973 + RDCOUNT
                    CALL random_initialize ( seed )
                    CALL rvec_random ( - 1.0_WP, 1.0_WP, 3, XRD )
                    !>@note the scaled random perturbations, in three directions.
                    Q_tg(I, J, K, 1) = VPERQ * XRD(1)
                    Q_tg(I, J, K, 2) = VPERQ * XRD(2) / RNDI1(JJ)
                    Q_tg(I, J, K, 3) = VPERQ * XRD(3) / RCCI1(JJ)
                    !WRITE(*, '(3I4.2, 3ES13.5)') I, J, K,XRD(1), XRD(2), XRD(3) !test

                END DO
            END DO
            IF(JJ == 1) Q_tg(:, J, :, 2) = 0.0_WP
        END DO

        PR_tg(:, :, :) = 0.0_WP ! check 1 or 0
    END IF

    IF(iRandomType == flag_random_sine) THEN
        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                DO J = 1, N2DO(MYID)
                    JJ = JCL2G(J)
                    Q_tg(I, J, K, 1) = DSIN(XND_tg(I)) * DCOS(YCC(JJ)) * DCOS(ZCC(K))
                    Q_tg(I, J, K, 2) = - DCOS(XCC_tg(I)) * DSIN(YND(JJ)) * DCOS(ZCC(K))
                    Q_tg(I, J, K, 3) = 0.0_WP
                    PR_tg(I, J, K) = (DCOS(2.0_WP * XCC_tg(I)) + DCOS(2.0_WP * YCC(JJ))) * &
                    (DCOS(2.0_WP * ZCC(K) + 2.0_WP)) / 16.0_WP
                END DO
            END DO
        END DO
    END IF



    RETURN
END SUBROUTINE

!***********************************************************************************
!SUBROUTINE VELO2RVELO_TG
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        INTEGER(4) :: I, J, K, NYI, JJ

!        NYI = 1
!        IF(MYID == 0) THEN
!            J = 1
!            DO I = 1, NCL1_TG
!                DO K = 1, NCL3
!                    Q_tg(I, J, K, 1) = Q_tg(I, J, K, 1)
!                    Q_tg(I, J, K, 2) = 0.0_WP
!                    Q_tg(I, J, K, 3) = Q_tg(I, J, K, 3) / RCCI1(JJ)
!                END DO
!            END DO
!        NYI = 2
!        END IF

!        DO J = NYI, N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I = 1, NCL1_TG
!                DO K = 1, NCL3
!                    Q_tg(I, J, K, 1) = Q_tg(I, J, K, 1)
!                    Q_tg(I, J, K, 2) = Q_tg(I, J, K, 2) / RNDI1(JJ)
!                    Q_tg(I, J, K, 3) = Q_tg(I, J, K, 3) / RCCI1(JJ)
!                END DO
!            END DO
!        END DO

!        RETURN
!    END SUBROUTINE
