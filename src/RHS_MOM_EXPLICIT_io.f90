!**********************************************************************************************************************************
!> @brief
!>
!> @details
!> SUBROUTINE: RHS_MOM_EXPLICIT_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE RHS_MOM_EXPLICIT_io(NS, IDR) ! not using other Gs or Qs in the current step
    USE FLOW_INFO
    USE THERMAL_INFO
    USE MESH_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4), INTENT(IN) :: IDR

    REAL(WP) :: EXPLT_io(NCL1_io, N2DO(0), NCL3)
    REAL(WP) :: DEN0
    REAL(WP) :: COE2
    REAL(WP) :: RHSC, RHSL
    REAL(WP) :: PGM
    REAL(WP) :: DPGRNS, INTGRHSY, INTGRHSY_WORK
    INTEGER(4) :: NXI, NYI
    INTEGER(4) :: IC, IM
    INTEGER(4) :: JC, JM, JJ
    INTEGER(4) :: KC, KM


    !=============== SET UP INDEX withOUT B.C.===========================
    NXI = 1
    IF(TgFlowFlg .AND. (IDR == 1) )  NXI = 2
    NYI = 1
    IF(IDR == 2 .AND. MYID == 0 .AND. iCase /= IBox3P) NYI = 2

    !=============== SET UP THE CONVECTION TERM INTO ONE variable ========
    EXPLT_io = 0.0_WP
    IF(IDR == 1) THEN

        DO IC = NXI, NCL1_io
            DO JC = NYI, N2DO(MYID)
                DO KC = 1, NCL3
                    EXPLT_io(IC, JC, KC) = QTMP_io(IC, JC, KC)
                END DO
            END DO
        END DO

    ELSE IF(IDR == 2) THEN

        DO IC = NXI, NCL1_io
            DO JC = NYI, N2DO(MYID)
                DO KC = 1, NCL3
                    EXPLT_io(IC, JC, KC) = DPH_io(IC, JC, KC)
                END DO
            END DO
        END DO

    ELSE IF(IDR == 3) THEN

        DO IC = NXI, NCL1_io
            DO JC = NYI, N2DO(MYID)
                DO KC = 1, NCL3
                    EXPLT_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC)
                END DO
            END DO
        END DO

    ELSE
    END IF

    !============ SET UP THE RHS with CONVECTION AND VISCOUS TERMS ================
    IF(iVisScheme == VisExplicit) THEN
        DO IC = NXI, NCL1_io
            DO JC = NYI, N2DO(MYID)
                DO KC = 1, NCL3
                    RHSC = EXPLT_io(IC, JC, KC)    !current
                    RHSL = EXPLT0_io(IC, JC, KC, IDR) ! LAST
                    RHS_io(IC, JC, KC) = (TGAM(NS) * RHSC + TROH(NS) * RHSL) * DT
                    EXPLT0_io(IC, JC, KC, IDR) = RHSC
                END DO
            END DO
        END DO
    ELSE IF(iVisScheme == VisImplicit) THEN
        COE2 = TALP(NS) * DT
        DO IC = NXI, NCL1_io
            DO JC = NYI, N2DO(MYID)
                DO KC = 1, NCL3
                    RHSC = EXPLT_io(IC, JC, KC)    !current
                    RHSL = EXPLT0_io(IC, JC, KC, IDR) ! LAST
                    RHS_io(IC, JC, KC) = (TGAM(NS) * RHSC + TROH(NS) * RHSL) * DT + &
                    COE2 * RHS_io(IC, JC, KC) + &
                    0.5_WP * COE2 * DIVU_io(IC, JC, KC)
                    EXPLT0_io(IC, JC, KC, IDR) = RHSC
                END DO
            END DO
        END DO
    ELSE
    END IF

    !=========== SET UP pressure GRADIENT =========================================
    PGM = 0.0_WP
    IF (IDR == 1) THEN
        COE2 = TALP(NS) * DT* DXI
        DO IC = NXI, NCL1_io
            IM = IMV_io(IC)
            DO KC = 1, NCL3
                DO JC = NYI, N2DO(MYID)
                    IF(iWeightedPre == 1) THEN
                        PGM = ( (PR0_io(IC, JC, KC, 2) - PR0_io(IM, JC, KC, 2)) * (0.25_WP - pres_epslon) + &
                        (PR0_io(IC, JC, KC, 1) - PR0_io(IM, JC, KC, 1)) * 0.5_WP + &
                        (PR_io(IC, JC, KC)   - PR_io(IM, JC, KC)) * (0.25_WP + pres_epslon)) * COE2
                    ELSE
                        PGM = ( PR_io(IC, JC, KC) - PR_io(IM, JC, KC) ) * COE2
                    END IF
                    RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE IF (IDR == 2) THEN

        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JJ = JCL2G(JC)
            COE2 = TALP(NS) * DT* DYCI(JJ) / RNDI1(JJ)
            DO IC = NXI, NCL1_io
                DO KC = 1, NCL3
                    IF(iWeightedPre == 1) THEN
                        PGM = ( (PR0_io(IC, JC, KC, 2) - PR0_io(IC, JM, KC, 2)) * (0.25_WP - pres_epslon) + &
                        (PR0_io(IC, JC, KC, 1) - PR0_io(IC, JM, KC, 1)) * 0.5_WP + &
                        (PR_io(IC, JC, KC)   - PR_io(IC, JM, KC))    * (0.25_WP + pres_epslon)) * COE2
                    ELSE
                        PGM = (PR_io(IC, JC, KC) - PR_io(IC, JM, KC)) * COE2
                    END IF
                    RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE IF (IDR == 3) THEN
        COE2 = TALP(NS) * DT* DZI
        DO KC = 1, NCL3
            KM = KMV(KC)
            DO JC = NYI, N2DO(MYID)
                DO IC = NXI, NCL1_io
                    IF(iWeightedPre == 1) THEN
                        PGM = ( (PR0_io(IC, JC, KC, 2) - PR0_io(IC, JC, KM, 2)) * (0.25_WP - pres_epslon) + &
                        (PR0_io(IC, JC, KC, 1) - PR0_io(IC, JC, KM, 1)) * 0.5_WP + &
                        (PR_io(IC, JC, KC)   - PR_io(IC, JC, KM))    * (0.25_WP + pres_epslon)) * COE2
                    ELSE
                        PGM = (PR_io(IC, JC, KC) - PR_io(IC, JC, KM)) * COE2
                    END IF
                    RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) - PGM
                END DO
            END DO
        END DO
    ELSE
    ENDIF

    !========= SET UP THE BUOYANCY FORCE =======================
    IF(ABS(iGravity) == IDR) THEN

        IF(ABS(iGravity) == 1) THEN ! vertical downwards/ upwards
            COE2 = TALP(NS) * DT * F_A
            DO IC = NXI, NCL1_io
                IM = IMV_io(IC)
                DO KC = 1, NCL3
                    DO JC = NYI, N2DO(MYID)
                        !DEN0 = (DENSITY(IC, JC, KC) + DENSITY(IM, JC, KC)) * XND2CL
                        RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) + D_STG(IC, JC, KC, IDR) * COE2
                        !WRITE(*, *) MYID, IC, KC, JC, DEN0, D_STG(IC, JC, KC, IDR), DEN0- D_STG(IC, JC, KC, IDR) !test
                    END DO
                END DO
            END DO
        END IF

        IF(ABS(iGravity) == 2) THEN ! horizontal flow
            DO JC = NYI, N2DO(MYID)
                JJ = JCL2G(JC)
                JM = JLMV(JC)
                DO KC = 1, NCL3
                    IF (iCase == 1 .OR. iCase == 4) COE2 = TALP(NS) * DT * F_A / RNDI1(JJ)
                    IF (iCase == 2 .OR. iCase == 3) COE2 = TALP(NS) * DT * F_A / RNDI1(JJ) * DCOS((DBLE(KC) - 0.5_WP) * DZ)
                    DO IC = NXI, NCL1_io
                        !DEN0 = YCL2ND_WFF(JJ) * DENSITY(IC, JC, KC) + &
                        !        YCL2ND_WFB(JJ) * DENSITY(IC, JM, KC)
                        RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) + D_STG(IC, JC, KC, IDR) * COE2
                    END DO
                END DO
            END DO
        END IF

        IF(ABS(iGravity) == 3) THEN
            COE2 = TALP(NS) * DT

            DO JC = NYI, N2DO(MYID)
                JJ = JCL2G(JC)
                COE2 = TALP(NS) * DT * F_A / RCCI1(JJ)
                DO KC = 1, NCL3
                    KM = KMV(KC)
                    DO IC = NXI, NCL1_io
                        !DEN0 = DENSITY(IC, JC, KC) + DENSITY(IC, JC, KM)
                        RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) + D_STG(IC, JC, KC, IDR) * COE2
                    END DO
                END DO
            END DO
        END IF

    END IF

    IF (iCase == 2 .OR. iCase == 3) THEN
        IF(ABS(iGravity) == 2 .AND. IDR == 3) THEN
            DO JC = NYI, N2DO(MYID)
                JJ = JCL2G(JC)
                DO KC = 1, NCL3
                    COE2 = TALP(NS) * DT * F_A / RNDI1(JJ) * (-DSIN((DBLE(KC - 1)) * DZ))
                    KM = KMV(KC)
                    DO IC = NXI, NCL1_io                 !TREATMENT OF CIRCUMFERENTIAL GRAVITY IN PIPE OR Annular FLOW
                        !DEN0 = DENSITY(IC, JC, KC) + DENSITY(IC, JC, KM)
                        RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) + D_STG(IC, JC, KC, IDR) * COE2
                    END DO
                END DO
            END DO
        END IF
    END IF

    !====================flow drive terms (source terms) in periodic Streamwise flow===========================
    IF ( (.NOT.TgFlowFlg) .AND. (IDR == NFLOW) .AND. (iCase /= IBox3P)) THEN

        DPGRNS = 0.0_WP
        !=============constant mass flow ratE =============================
        IF(iFlowDriven == 1) THEN
            INTGRHSY = 0.0_WP

            DO JC = 1, N2DO(MYID)
                JJ = JCL2G(JC)
                DO IC = 1, NCL1_io
                    DO KC = 1, NCL3
                        INTGRHSY = INTGRHSY + RHS_io(IC, JC, KC) / DYFI(JJ) / RCCI1(JJ)
                    END DO
                END DO
            END DO

            CALL MPI_BARRIER(ICOMM, IERROR)
            CALL MPI_ALLREDUCE(INTGRHSY, INTGRHSY_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
            IF(MYID == 0) DPGRNS = INTGRHSY_WORK / VOLM_io
        ELSE IF(iFlowDriven == 2) THEN
            !constant pressure gradient, dimensionless based on \Delta and U_m
            !DPGRNS = -2.0_WP      ! DIMENSIONLESS BASED ON U_TAU
            COE2 = TALP(NS) * DT
            IF(MYID == 0) DPGRNS = -0.50_WP * Cf_Given * COE2 !dimensionless based on \Delta and U_m
        ELSE
        END IF

        CALL MPI_BCAST( DPGRNS, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        DO KC = 1, NCL3
            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    RHS_io(IC, JC, KC) = RHS_io(IC, JC, KC) - DPGRNS
                END DO
            END DO
        END DO
    END IF

    !===== SAVE DATA BACK ===========================================
    IF(iVisScheme == VisExplicit) THEN
        IF(IDR == 1) THEN

            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    DO KC = 1, NCL3
                        QTMP_io(IC, JC, KC) = RHS_io(IC, JC, KC)
                    END DO
                END DO
            END DO
            !CALL DEBUG_WRT_LOCAL(QTMP_io, 1, N2DO(MYID), 'conx') ! test
        ELSE IF(IDR == 2) THEN

            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    DO KC = 1, NCL3
                        DPH_io(IC, JC, KC) = RHS_io(IC, JC, KC)
                    END DO
                END DO
            END DO
            !CALL DEBUG_WRT_LOCAL(DPH_io, 1, N2DO(MYID), 'cony') ! test
        ELSE IF(IDR == 3) THEN

            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    DO KC = 1, NCL3
                        RHSLLPHI_io(IC, JC, KC) = RHS_io(IC, JC, KC)
                    END DO
                END DO
            END DO
            !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io, 1, N2DO(MYID), 'conz') ! test
        ELSE
        END IF
    END IF

    RETURN
END SUBROUTINE
