!**********************************************************************************************************************************
!> @brief
!>        Flow field initialization of the tg domain
!> @details
!> SUBROUTINE: random_FL_THEML_FLD_io (in MYID = all)
!> SUBROUTINE: CALC_INITIALIZATION_io (in MYID = all)
!> SUBROUTINE: IniField_FLOW_io (in MYID = all)
!> SUBROUTINE: IniField_THERMAL_io (in MYID = all)
!> SUBROUTINE: BULK_VELOCITY_io
!> SUBROUTINE: BULK_MASSFLUX_io
!> SUBROUTINE: BULK_H_io
!> SUBROUTINE: VELO2MASSFLUX
!> SUBROUTINE: Unified_MassFlux
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010- Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013- Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE random_FL_THEML_FLD_io
        USE init_info
        USE mesh_info
        USE flow_info
        USE thermal_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE

        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        INTEGER(4) :: L
        REAL(WP) :: Umean1_io

        IF(.not.IoFlowFlg) RETURN

        !========Generate scaled random Q and PR = 1 ==========================
        !    INFLOW/OUTFLOW DOmAIN***
        IF(MYID == 0) CALL CHKHDL('(1-2) IO: Generating random velocity', MYID)
        IF(iThermoDynamics == 1) CALL IniField_THERMAL_io
        CALL IniField_FLOW_io
        CALL VMAV_io
        IF(MYID == 0) THEN
            CALL CHKHDL   ('IO: VMV(1:3) In random velocity field', MYID)
            CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_io(1), VMAX_io(2), VMAX_io(3))
            CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_io(1), VMIN_io(2), VMIN_io(3))
        END IF

        !>================ Q = Q-QxzmeaN =========================================
        !>    INFLOW/OUTFLOW DOmAIN***
        IF(iThermoDynamics == 1) CALL INTFC_MFD_THERMAL_io
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL BC_WALL_Q_io
        CALL INITIAL_MEAN_io
        DO L = 1, 3
            DO J = 1, N2DO(MYID)
                Q_io(:, J, :,L) = Q_io(:, J, :,L) - UU(J, L, 2)
            END DO
        END DO

        CALL VMAV_io
        IF(MYID == 0) THEN
            CALL CHKHDL('(2-2) IO: VMV(1:3) in random velocity after subtracting meanXZ ', MYID)
            CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_io(1), VMAX_io(2), VMAX_io(3))
            CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_io(1), VMIN_io(2), VMIN_io(3))
        END IF

        IF(iCase /= IBox3P) THEN !!!added
        !>===========update Q in the incoming flow direction ====================
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            DO I = 1, NCL1E
                DO K = 1, NCL3
                    Q_io(I, J, K, NFLOW) = Q_io(I, J, K, NFLOW) + Vini(JJ)
                END DO
            END DO
        END DO

        !IF(iCase /= ICHANL) CALL VELO2RVELO_io

        !>    @note Set Q(I, J, K, 2) = v at Wall b.c. be zero.
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL BC_WALL_Q_io

        !==============CORRECT THE MASS FLOW RATE = 1 =========================================
        CALL BULK_VELOCITY_io(Umean1_io)
        IF(MYID == 0) CALL CHKRLHDL  ('IO: The bulk velocity (original) = ', MYID, Umean1_io)

        DO J = 1, N2DO(MYID)
            Q_io(:, J, :, NFLOW) = Q_io(:, J, :, NFLOW) / Umean1_io
        END DO

        CALL BULK_VELOCITY_io(Umean1_io)
        IF(MYID == 0) CALL CHKRLHDL  ('IO: The bulk velocity (corrected) = ', MYID, Umean1_io)
        END IF

        !>============ SCALING DUE TO NOn -INLET Reference ======================================
        Q_io(:, :, :, :) = Q_io(:, :, :, :) * M_inlet / D_inlet
        CALL BULK_VELOCITY_io(Umean1_io)
        IF(MYID == 0) CALL CHKRLHDL  ('IO: The bulk velocity (scaled) = ', MYID, Umean1_io)
        !>============CHECK DIVERGENCE ======================================
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL BC_WALL_Q_io

        CALL VELO2MASSFLUX
        IF(TgFlowFlg) CALL BC_TINLET_FLOW

        CALL Unified_MassFlux

        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL BC_WALL_Q_io
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
        CALL BC_WALL_G_io
        CALL VMAV_io
        IF(MYID == 0) THEN
            CALL CHKHDL   ('(3-2) IO: VMV(1:3) in the real initial flow field ', MYID)
            CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_io(1), VMAX_io(2), VMAX_io(3))
            CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_io(1), VMIN_io(2), VMIN_io(3))
        END IF

        CALL DIVGCK_io
        IF(MYID == 0) THEN
            CALL CHKHDL   ('(4-2) IO: Max divergence of initial flow field. main domain, inlet, outlet', MYID)
            CALL CHK3RLHDL('      Div(Velocity) = ', MYID, MAXDIVGV_io(1), MAXDIVGV_io(2), MAXDIVGV_io(3))
        END IF

        !CALL CALL_TEC360

        RETURN

    END SUBROUTINE

!*********************************************************************************************************************
SUBROUTINE CALC_INITIALIZATION_io ! initial thermal field IS not calculated for random reStart...
        USE init_info
        USE mesh_info
        USE flow_info
        USE thermal_info
        IMPLICIT NONE

        INTEGER(4) :: L

        IF(.not.IoFlowFlg) RETURN


        !======Calculate the QP fields ================================
        IF(MYID == 0) CALL CHKHDL('(4-2) IO: Initial flow field, calc QP...', MYID)
        CALL VMAV_io
        IF(MYID == 0) THEN
            CALL CHKHDL('(4-3) IO: Initial flow field before first calcuted velocity, VMAX(1:3)', MYID)
            CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_io(1), VMAX_io(2), VMAX_io(3))
            CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_io(1), VMIN_io(2), VMIN_io(3))
        END IF
        IF(iThermoDynamics == 1) THEN
            CALL BC_WALL_THERMAL(IALLDM)

            !CALL SOLVERRK3_ENG_io(0)! Junjie, IF H_ref IS not as inlet, energy equation may caUSE unphysICal behavior

            CALL DENSITY_Staggered
            CALL MU_Staggered
        END IF
        CALL SOLVERRK3_MOM_io(0)

        !==============CHECK INFO========================================
        CALL VMAV_io
        IF(MYID == 0) THEN
            CALL CHKHDL('(5-2) IO: Initial flow field in the first calcuted velocity, VMAX(1:3)', MYID)
            CALL CHK3RLHDL('      Max U, V, W = ', MYID, VMAX_io(1), VMAX_io(2), VMAX_io(3))
            CALL CHK3RLHDL('      Min U, V, W = ', MYID, VMIN_io(1), VMIN_io(2), VMIN_io(3))
        END IF

        CALL DIVGCK_io
        IF(MYID == 0) THEN
            CALL CHKHDL('(6-2) IO: Max divergence of final flow field. main domain, inlet, outlet', MYID)
            CALL CHK3RLHDL('      Div(Velocity) = ', MYID, MAXDIVGV_io(1), MAXDIVGV_io(2), MAXDIVGV_io(3))
        END IF
        !MAXDIVGV_io(1) = MAXDIVGV_io(2)

        RETURN
    END SUBROUTINE
!*********************************************************************************************************************
SUBROUTINE IniField_FLOW_io
        USE init_info
        USE mesh_info
        USE flow_info
        IMPLICIT NONE


        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        INTEGER(4) :: RDCOUNT
        INTEGER(4) :: seed
!       INTEGER(4) :: IDUM
        REAL(WP) :: VPERQ
        !REAL(WP) :: XRDM(3)
        REAL(WP) :: XRD(3)

        Q_io = 0.0_WP
        PR_io = 0.0_WP
        !******************************************************************
        IF(iRandomType /= flag_random_sine) THEN
            XRD = 0.0_WP
            RDCOUNT = 0
            SEED = 0
            !>     @note Generate the random flow field u, V,w.
            DO J = 1, N2DO(MYID)           ! only local y_Cell no.
                JJ = JCL2G(J)
            !>          @pARam VPERQ = Scaled magnitude of perturbation based on Up (pARabolIC max. velo.)
            !>          @note for 1 /4 neAR wall region, perturbation ratio IS decreased by 25%.
                VPERQ = VPERG
                IF((1.0_WP - DABS(YCC(JJ))) < 0.250_WP) VPERQ = SVPERG

                DO I = 1, NCL1E
                    DO K = 1, NCL3
                        !RDCOUNT = RDCOUNT + 1
                        RDCOUNT = I + K + JJ  ! make it independent of processor number.
                        IF(iRandomType == flag_random_real) THEN
                          seed = 0
                        ELSE IF(iRandomType == flag_random_fixed) THEN
                          seed = 1973 + RDCOUNT + (NCL1_TG * NCL3 * NCL2) + 1024
                        ELSE
                          seed = RDCOUNT
                        END IF
                        CALL random_initialize ( seed )
                        CALL rvec_random ( - 1.0_WP, 1.0_WP, 3, XRD )
            !>                   @note the scaled random perturbations, in three directions.
                        Q_io(I, J, K, 1) = VPERQ * XRD(1)
                        Q_io(I, J, K, 2) = VPERQ * XRD(2) / RNDI1(JJ)
                        Q_io(I, J, K, 3) = VPERQ * XRD(3) / RCCI1(JJ)
                        !WRITE(*, '(3I4.2, 3ES13.5)') I, J, K,XRD(1), XRD(2), XRD(3) !test
            !>               @note No perturbulation IS added to pressure
                    END DO
                END DO
                IF(JJ == 1) Q_io(:, J, :, 2) = 0.0_WP
            END DO
            PR_io(:, :, :) = 0.0_WP ! check 0 or 1
        END IF

        !******************************************************************
        IF(iRandomType == flag_random_sine) THEN
            DO I = 1, NCL1E
                DO K = 1, NCL3
                    DO J = 1, N2DO(MYID)
                        JJ = JCL2G(J)
                        Q_io(I, J, K, 1) = DSIN(XND_io(I)) * DCOS(YCC(JJ)) * DCOS(ZCC(K))
                        Q_io(I, J, K, 2) = - DCOS(XCC_io(I)) * DSIN(YND(JJ)) * DCOS(ZCC(K))
                        Q_io(I, J, K, 3) = 0.0_WP
                        PR_io(I, J, K) = (DCOS(2.0_WP * XCC_io(I)) + DCOS(2.0_WP * YCC(JJ))) * &
                                       (DCOS(2.0_WP * ZCC(K) + 2.0_WP)) / 16.0_WP
                    END DO
                END DO
            END DO

        END IF
        !*****************************************************************************



        IF(iWeightedPre == 1) THEN
            PR0_io(:, :, :, 1) = PR_io(:, :, :)
            PR0_io(:, :, :, 2) = PR_io(:, :, :)
        END IF


        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
        CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)

        ! TEST ===========

        !

        CALL VELO2MASSFLUX

        !CALL DEBUG_WRT_QGP_io

        RETURN
    END SUBROUTINE


!*********************************************************************************************************************
SUBROUTINE IniField_THERMAL_io
!>    @note
!>          Set up the initial thermal field.
!>    Known:
!>           P_0, T_0
!>           P_0 IS only USEd to creat the NIST table, and not used in the code.
!>    To set up:
!>           T, \rho, \mu, \kappa, h, h\rho
        USE MESH_INFO
        USE thermal_info
        USE init_info
        IMPLICIT NONE

        !INTEGER(4) :: I, J, K, JJ, IE, IS
        !REAL(WP) :: T_tmp, SFS, SFE, H_tmp, D_tmp, K_tmp, M_tmp

        NTHERMAL = 7


        DH (:, :, :)     = DH_inlet
        DH0(:, :, :)     = DH_inlet

        ENTHALPY   (:, :, :) = H_inlet
        TEMPERATURE(:, :, :) = T_inlet
        THERMCONDT (:, :, :) = K_inlet
        HEATCAP    (:, :, :) = CP_inlet

        Viscousity (:, :, :)    = M_inlet
        Viscousity0(:, :, :)    = M_inlet
        MU_STG     (:, :, :, :) = M_inlet

        DENSITY (:, :, :)    = D_inlet
        DENSITY0(:, :, :)    = D_inlet
        D_STG   (:, :, :, :) = D_inlet


        RETURN
    END SUBROUTINE


    !*********************************************************************************************************************
SUBROUTINE BULK_VELOCITY_io(UUMEAN_WORK)
        USE init_info
        USE mesh_info
        USE flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE

        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP) :: UUMEAN, UUMEAN_WORK

        UUMEAN = 0.0_WP
        DO J = 1, N2DO(MYID)  !@
            JJ = JCL2G(J)
            DO I = 1, NCL1_io
                DO K = 1, NCL3
                    UUMEAN = UUMEAN + Q_io(I, J, K, NFLOW) / DYFI(JJ) / RCCI1(JJ)
                END DO
            END DO
        END DO
        UUMEAN = UUMEAN / DZI
        CALL MPI_BARRIER(ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UUMEAN, UUMEAN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

        UUMEAN_WORK = UUMEAN_WORK / (Area_inlet * DBLE(NCL1_io))

        RETURN
    END SUBROUTINE

!*********************************************************************************************************************
SUBROUTINE BULK_MASSFLUX_io(GMEAN_WORK)
        USE init_info
        USE mesh_info
        USE flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE

        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP) :: GMEAN, GMEAN_WORK

        GMEAN = 0.0_WP
        DO J = 1, N2DO(MYID)  !@
            JJ = JCL2G(J)
            DO I = 1, NCL1_io
                DO K = 1, NCL3
                    GMEAN = GMEAN + G_io(I, J, K, NFLOW) / DYFI(JJ) / RCCI1(JJ)
                END DO
            END DO
        END DO
        GMEAN = GMEAN / DZI
        CALL MPI_BARRIER(ICOMM, IERROR)
        CALL MPI_ALLREDUCE(GMEAN, GMEAN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

        GMEAN_WORK = GMEAN_WORK / (Area_inlet * DBLE(NCL1_io))

        RETURN
    END SUBROUTINE

!*********************************************************************************************************************
SUBROUTINE BULK_H_io(HMEAN_WORK)
        USE init_info
        USE mesh_info
        USE flow_info
        USE thermal_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE

        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP) :: HMEAN, HMEAN_WORK

        HMEAN = 0.0_WP
        DO J = 1, N2DO(MYID)  !@
            JJ = JCL2G(J)
            DO I = 1, NCL1_io
                DO K = 1, NCL3
                    HMEAN = HMEAN + ENTHALPY(I, J, K) / DYFI(JJ) / RCCI1(JJ)
                END DO
            END DO
        END DO
        HMEAN = HMEAN / DZI
        CALL MPI_BARRIER(ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HMEAN, HMEAN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

        HMEAN_WORK = HMEAN_WORK / (Area_inlet * DBLE(NCL1_io))

        RETURN
    END SUBROUTINE



!*********************************************************************************************************
SUBROUTINE VELO2MASSFLUX
        USE MESH_INFO
        USE THERMAL_INFO
        USE FLOW_INFO
        USE init_info
        IMPLICIT NONE
        INTEGER(4) :: IDR, NII, IC, JC, KC, IM, JJ, KM, JM, KS
        !REAL(WP) :: RHO


!>      !===================================== X =============================
        IDR = 1
        DO IC = 1, NCL1_io
            IM = IMV_io(IC)
            DO JC = 1, N2DO(MYID)
                DO KC = 1, NCL3
                    !RHO = ( DENSITY(IC, JC, KC) + DENSITY(IM, JC, KC) ) * XND2CL
                    G_io(IC, JC, KC, IDR) = Q_io(IC, JC, KC, IDR) * D_STG(IC, JC, KC, IDR)!* RHO
                END DO
            END DO
        END DO

!>      !===================================== Z =============================
        IDR = 3
        DO KC = 1, NCL3
            KM = KMV(KC)
            DO JC = 1, N2DO(MYID)
                DO IC = 1, NCL1_io
                    !RHO = ( DENSITY(IC, JC, KC) + DENSITY(IC, JC, KM) ) * ZND2CL
                    G_io(IC, JC, KC, IDR) = Q_io(IC, JC, KC, IDR) * D_STG(IC, JC, KC, IDR)!* RHO
                END DO
            END DO
        END DO

!>      !=====================================Y =============================
        IDR = 2
        NII = 1
        IF(MYID == 0 .AND. iCase /= IBox3P) NII = 2
        DO JC = NII, N2DO(MYID)
            JJ = JCL2G(JC)
            JM = JLMV(JC)
            DO KC = 1, NCL3
                DO IC = 1, NCL1_io
                    !RHO = YCL2ND_WFF(JJ) * DENSITY(IC, JC, KC) + &
                    !      YCL2ND_WFB(JJ) * DENSITY(IC, JM, KC)
                    G_io(IC, JC, KC, IDR) = Q_io(IC, JC, KC, IDR) * D_STG(IC, JC, KC, IDR)!* RHO
                END DO
            END DO
        END DO

        IF (MYID == 0 .AND. iCase /= IBox3P) THEN
            DO KC = 1, NCL3
                KS = KSYM(KC)
                DO IC = NCL1S, NCL1E
                    IF(iCase == ICHANL) THEN
                        G_io(IC, 1, KC, IDR) = 0.0_WP
                        G_io(IC, 0, KC, IDR) = 0.0_WP ! not used!
                    END IF
                    IF(iCase == IPIPEC) THEN
                        G_io(IC, 1, KC, IDR) = 0.0_WP
                        G_io(IC, 0, KC, IDR) = G_io(IC, 1, KS, IDR)
                    END IF
                END DO
            END DO
        ENDIF

        IF (MYID == NPSLV) THEN
            DO KC = 1, NCL3
                DO IC = NCL1S, NCL1E
                    IF(iCase == IBox3P) THEN
                        G_io(IC, N2DO(MYID) + 1, KC, IDR) = Q_io(IC, N2DO(MYID) + 1, KC, IDR) * D_STG(IC, N2DO(MYID) + 1, KC, IDR)
                    ELSE
                        G_io(IC, N2DO(MYID) + 1, KC, IDR) = 0.0_WP
                    END IF
                END DO
            END DO
        ENDIF

        RETURN
    END SUBROUTINE


!SUBROUTINE INTFC_G_WALL
!        USE MESH_INFO
!        USE THERMAL_INFO
!        USE FLOW_INFO
!        USE init_info
!        IMPLICIT NONE
!        INTEGER(4) :: IDR, NII, IC, JC, KC, IM, JJ, KM, JM, KS, JP
!        REAL(WP) :: RHO1,RHO2, RHO0
!        REAL(WP) :: RHO_ND2, RHO_ND0, RHO_WAL, RHO_WL0

!        IF(MYID == 0) THEN
!            JC = 1
!            JP = 0
!        ELSE IF(MYID == NPSLV) THEN
!            JC = N2DO(MYID)
!            JP = N2DO(MYID) + 1
!        ELSE
!            RETURN
!        END IF

!        !===================================== X =============================
!        IDR = 1
!        DO IC = NCL1S, NCL1E
!            IF(IC == 0) THEN
!                IM = 0
!            ELSE
!                IM = IMV_io(IC)
!            END IF

!            DO KC = 1, NCL3
!                RHO1 = ( DENSITY(IC, JC, KC) + DENSITY(IM, JC, KC) ) * XND2CL
!                RHO2 = ( DENSITY(IC, JP, KC) + DENSITY(IM, JP, KC) ) * XND2CL
!                Q_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) / RHO1
!                G_io(IC, JP, KC, IDR) = Q_io(IC, JC, KC, IDR) * RHO2 * (-1.0_WP)
!            END DO

!            !Q_io(:, JP, :, IDR) = -Q_io(:, JC, :, IDR) !test
!            !G_io(:, JP, :, IDR) = - G_io(:, JC, :, IDR) !test

!        END DO

!!>      !===================================== Z =============================
!        IDR = 3
!        DO KC = 1, NCL3
!            KM = KMV(KC)
!            DO IC = NCL1S, NCL1E
!                RHO1 = ( DENSITY(IC, JC, KC) + DENSITY(IC, JC, KM) ) * ZND2CL
!                RHO2 = ( DENSITY(IC, JP, KC) + DENSITY(IC, JP, KM) ) * ZND2CL
!                Q_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) / RHO1
!                G_io(IC, JP, KC, IDR) = Q_io(IC, JC, KC, IDR) * RHO2 * (-1.0_WP)
!            END DO

!            !Q_io(:, JP, :, IDR) = -Q_io(:, JC, :, IDR) !test
!            !G_io(:, JP, :, IDR) = - G_io(:, JC, :, IDR) !test

!        END DO

!!>      !=====================================Y =============================
!        IF(MYID == 0) THEN
!            IDR = 2
!            DO KC = 1, NCL3
!                KM = KMV(KC)
!                DO IC = NCL1S, NCL1E

!                    RHO_WAL = ( DENSITY (IC, 0, KC) + DENSITY(IC, 1, KC)  ) * 0.5_WP
!                    RHO_WL0 = ( DENSITY0(IC, 0, KC) + DENSITY0(IC, 1, KC) ) * 0.5_WP
!                    RHO_ND2 =  DENSITY(IC, 1, KC) * YCL2ND_WFF(2) + DENSITY(IC, 2, KC) * YCL2ND_WFF(2)
!					RHO_ND0 = 2.0_WP * RHO_WAL- RHO_ND2

!                    Q_io(IC, 2, KC, IDR) = G_io(IC, 2, KC, IDR) / RHO_ND2
!                    Q_io(IC, 0, KC, IDR) = (-1.0_WP) * Q_io(IC, 2, KC, IDR)!-2.0_WP / DT* (RHO_WAL- RHO_WL0) / RHO_WAL
!                    G_io(IC, 0, KC, IDR) = Q_io(IC, 0, KC, IDR) * RHO_ND0
!                END DO
!            END DO

!            !Q_io(:, 0:1, :, IDR) = 0.0_WP !test
!            !G_io(:, 0:1, :, IDR) = 0.0_WP !test

!        END IF

!        IF(MYID == NPSLV) THEN
!            IDR = 2
!            Q_io(:, N2DO(MYID) + 1, :, IDR) = 0.0_WP
!            G_io(:, N2DO(MYID) + 1, :, IDR) = 0.0_WP
!        END IF

!        RETURN
!    END SUBROUTINE

!*********************************************************************
SUBROUTINE Unified_MassFlux
        USE MESH_INFO
        USE THERMAL_INFO
        USE FLOW_INFO
        USE init_info
        IMPLICIT NONE

        REAL(WP) :: GMEAN1, GMEAN2, Umean1_io

        IF(iCase == IBox3P) RETURN

        GMEAN1 = 0.0_WP
        CALL BULK_MASSFLUX_io(GMEAN1)
        IF(MYID == 0) CALL CHKRLHDL('The initial mass flux = ', MYID, GMEAN1)


        Q_io(:, :, :, 1:3) = Q_io(:, :, :, 1:3) /GMEAN1
        CALL BULK_VELOCITY_io(Umean1_io)
        IF(MYID == 0) CALL CHKRLHDL('The unified Umean = ', MYID, Umean1_io)


        CALL VELO2MASSFLUX

        CALL MPI_BARRIER(ICOMM, IERROR)
        GMEAN2 = 0.0_WP
        CALL BULK_MASSFLUX_io(GMEAN2)
        IF(MYID == 0) CALL CHKRLHDL('The unified mass flux = ', MYID, GMEAN2)

        RETURN
    END SUBROUTINE

!!***********************************************************************************
!SUBROUTINE VELO2RVELO_io
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        INTEGER(4) :: I, J, K, NYI, JJ

!        NYI = 1
!        IF(MYID == 0) THEN
!            J = 1
!            DO I = 1, NCL1_io
!                DO K = 1, NCL3
!                    Q_io(I, J, K, 1) = Q_io(I, J, K, 1)
!                    Q_io(I, J, K, 2) = 0.0_WP
!                    Q_io(I, J, K, 3) = Q_io(I, J, K, 3) / RCCI1(JJ)
!                END DO
!            END DO
!        NYI = 2
!        END IF

!        DO J = NYI, N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I = 1, NCL1_io
!                DO K = 1, NCL3
!                    Q_io(I, J, K, 1) = Q_io(I, J, K, 1)
!                    Q_io(I, J, K, 2) = Q_io(I, J, K, 2) / RNDI1(JJ)
!                    Q_io(I, J, K, 3) = Q_io(I, J, K, 3) / RCCI1(JJ)
!                END DO
!            END DO
!        END DO

!        RETURN
!    END SUBROUTINE
