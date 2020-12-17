!**********************************************************************************************************************************
!> @brief
!>        The main CFD solver
!> @details
!> SUBROUTINE: SOLVE (in MYID = all)
!> SUBROUTINE: BCAST_COMM_STEP (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE SOLVE
    USE init_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    REAL(WP) :: ENDTIME(2)
    REAL(WP) :: StartTIME
    REAL(WP) :: RENTMP

    INTEGER(4) :: NS
    INTEGER(4) :: ITERL

    !========= INITIAL TIME/STEP SETTING UP =======================================
    ITERL = 0
    IF(TgFlowFlg) THEN
        ITERG0_TG  = 0
        PhyTIME_TG = 0.0_WP

        CONVH0_tg = 0.0_WP
    END IF

    IF(IoFlowFlg) THEN
        ITERG0_io = 0
        PhyTIME_io = 0.0_WP

        EXPLT0_io   = 0.0_WP
        IF(iThermoDynamics == 1) RHS_ENERGY0 = 0.0_WP
    END IF


    !========FLOW INITIALIZATION, EITHER FROM RANDOM FIELDS OR FROM RESTART===========
    CALL FLOWStart
    RENTMP = REN
    DT0 = DT
    IF(IoFlowFlg) THEN
        IF(iIniField_io == IniField_reStart .AND. iIniFieldTime == 1) THEN
            ITERG0_io = 0
            PhyTIME_io = 0.0_WP
        END IF
        PhyTIME = PhyTIME_io
        ITERG0  = ITERG0_io
    ELSE
        PhyTIME = PhyTIME_TG
        ITERG0  = ITERG0_TG
    END IF
    dtSave  = dtSave1

    !======== SCREEN PRINGITNG=====================================================
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL PP_MONITOR_INI
    CALL MPI_BARRIER(ICOMM, IERROR)

    !========= main STEPS ADVANCING========= main SOLVE =============================
    !NTSTF = IDINT(TSTOP / DT) *10000
    NTSTF = MIN(NTSTF, IDINT(TSTOP / DT))
    IF(MYID == 0) THEN
        CALL CHKINTHDL(' The total iterations to calculation (given) = ', MYID, NTSTF)
        CALL CHKINTHDL(' The current iterations = ', MYID, ITERG0)
        CALL CHKHDL('*** End oF Code Configuration ***', MYID)
    END IF
    
    DO ITERG= ITERG0 + 1, NTSTF
        StartTIME = MPI_WTIME()

        !=========== SET UP CFL ===============
        IF(TgFlowFlg) CALL CFL_tg
        IF(IoFlowFlg) CALL CFL_io
        IF(IoFlowFlg .AND. TgFlowFlg) THEN
            CFLMM = DMAX1(CFLMM_tg, CFLMM_io)
        ELSE
            IF(TgFlowFlg)  CFLMM = CFLMM_tg
            IF(IoFlowFlg)  CFLMM = CFLMM_io
        END IF

        !========== ShARe the same common info for approaching=====
        IF ( PhyTIME < TLgRe ) THEN
            REN  = ReIni
        ELSE
            REN  = RENTMP
        ENDIF
        CALL BCAST_COMM_STEP

        !============ STEP AND TIME CONTROL =======================
        ITERL    = ITERL    + 1
        PhyTIME  = PhyTIME  + DT
        IF(TgFlowFlg) PhyTIME_tg = PhyTIME
        IF(IoFlowFlg) PhyTIME_io = PhyTIME

        !            !=========== SET UP RECORDING TIME INTERVAL ===============
        IF(PhyTIME < tRunAve1) THEN
            dtSave  = dtSave1 * 2.0_WP
        ELSE
            dtSave  = dtSave1
        END IF


        !============= RK ===== The main pART ===========================
        ! about how to solve 'n -coupled' dIFfeRENtial equations using RK method:
        ! ref: https://www.researchgate.net/post/How_DO_you_solve_n -coupled_diffeRENtial_equations_using_RK4
        ! ref: https:// Math.stackexchange.coM /questions/721076/helP -with- Using- ThE - RungE - Kutta-4th-order-methoD -on - A-systeM -oF -2 -firsT -order-od
        DO NS = 1, NSST
            IF(TgFlowFlg) CALL SOLVERRK3_MOM_tg(NS)


            ! in each stage of RK, the energy equation and the momentum equation IS decoupled.
            IF(IoFlowFlg .AND. (iThermoDynamics == 1) .AND. (PhyTIME > timeThermoStart)) THEN

                CALL SOLVERRK3_ENG_io(NS)
                CALL DENSITY_TIME_DERIVATION(NS)
                CALL DENSITY_Staggered
                CALL MU_Staggered
                IF(iVisScheme == VisImplicit) THEN
                    CALL DENSITY_IMPLICIT_add
                END IF
            END IF
            IF(IoFlowFlg .AND. (PhyTIME > timeFlowStart)) THEN
                !======= STEP6: CALCULATE VELOCITY =================
                !CALL VELOCITY_CALC_io
                CALL SOLVERRK3_MOM_io(NS)
            END IF
        END DO

        !===================================================
        ENDTIME(1) = MPI_WTIME()
        CPUTIME_tmp =ENDTIME(1) - StartTIME


        !=======POSTPROCESS =============================
        IF( DMOD(PhyTIME, dtRawView) < DT) &
        CALL CALL_TEC360
        IF(TgFlowFlg) CALL POSTPROCESS_tg
        IF(IoFlowFlg) CALL POSTPROCESS_io

    END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BCAST_COMM_STEP
    USE init_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    REAL(WP) :: COMMINFO(5)
    REAL(WP) :: DT0HALF

    COMMINFO = 0.0_WP
    IF(MYID == 0) THEN
        !============ REYNOLDS NUMBER AND PRT0 =====================

        CVISC  = 1.0_WP / REN
        IF(IoFlowFlg .AND. (iThermoDynamics == 1)) CTHECD = 1.0_WP / REN /PRT0

        !==========CALCULATE DT FROM GIVEN DT AND CFL ============
        IF(PhyTIME < tRunAve1) THEN
            DT0HALF = DT0
        ELSE
            DT0HALF = DT0 !* 0.5_WP
        END IF

        DT = DMIN1(DT0HALF, 5.0_WP * DT,CFLGV/ CFLMM) !DT can not exceed five times of the last step
        DT = DMAX1(DT, DTMIN)

        COMMINFO(1) = REN
        COMMINFO(2) = CVISC
        COMMINFO(3) = CTHECD
        COMMINFO(4) = CFLMM
        COMMINFO(5) = DT
    END IF

    CALL MPI_BCAST( COMMINFO, 5, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    IF(MYID /= 0) THEN
        REN  = COMMINFO(1)
        CVISC  = COMMINFO(2)
        CTHECD = COMMINFO(3)
        CFLMM  = COMMINFO(4)
        DT   = COMMINFO(5)
    END IF


    ! IN CASE OF ANY SMALL DT DUE TO DIVERGENCE
    IF((DT/ DT0) < 0.01_WP) THEN

        CALL CALL_TEC360
        IF(TgFlowFlg) CALL POSTPROCESS_tg
        IF(IoFlowFlg) CALL POSTPROCESS_io

        CALL MPI_BARRIER(ICOMM, IERROR)
        STOP
    END IF


    RETURN
END SUBROUTINE
