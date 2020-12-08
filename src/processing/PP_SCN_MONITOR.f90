!**********************************************************************************************************************************
!> @brief
!>        monitor
!> @details
!> SUBROUTINE: PP_MONITOR_TG
!> SUBROUTINE: PP_MONITOR_Xperiodic_io
!> SUBROUTINE: PP_MONITOR_nonXperiodic_io
!> SUBROUTINE: PP_wall_thermal_shear
!> SUBROUTINE: PP_Wall_thermal_properties
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PP_MONITOR_INI
    USE init_info
    USE mesh_info
    USE WRT_INFO
    USE thermal_info
    IMPLICIT NONE

    REAL(WP) :: dxPlus, yMinPlus, dzPlus, L2, HY, dyMin, dyMax, yMaxPlus

    IF(MYID /= 0) RETURN

    !======================== MESH INFO==================================================================
    IF(iCase == ICHANL .OR. iCase == IBox3P) THEN
        dyMin = 1.0_WP / DYFI(1)
        dyMax = 1.0_WP / DYFI(NCL2 / 2)
        HY  = 2.0_WP * ALX2
    ELSE
        dyMin = 1.0_WP / DYFI(NCL2)
        dyMax = 1.0_WP / DYFI(1)
        HY  = 2.0_WP * PI * ALX2
    END IF
    dxPlus = DX * REN * DSQRT(0.5_WP * Cf_Given)
    dzPlus = DZ * REN * DSQRT(0.5_WP * Cf_Given)
    yMinPlus = dyMin * REN * DSQRT(0.5_WP * Cf_Given)
    yMaxPlus = dyMax * REN * DSQRT(0.5_WP * Cf_Given)


    WRITE(logflg_pg, '(A)') '#********** Mesh and flow details based on initial flow field******************************'
    IF(TgFlowFlg) &
    WRITE(logflg_pg, '(A, I5, 3(2X, A, F9.5))') '#   NCL1_tg = ', NCL1_tg, 'HX_TG = ', ALX1(1), 'DX =  ', DX,    'DX += ', dxPlus
    IF(IoFlowFlg) &
    WRITE(logflg_pg, '(A, I5, 3(2X, A, F9.5))') '#   NCL1_io = ', NCL1_io, 'HX_io = ', ALX1(2), 'DX =  ', DX,    'DX += ', dxPlus
    WRITE(logflg_pg, '(A, I5, 3(2X, A, F9.5))') '#   NCL3 = ', NCL3,    'HZ = ', ALX3,    'DZ =  ', DZ,    'DZ += ', dzPlus
    WRITE(logflg_pg, '(A, I5, 5(2X, A, F9.5))') '#   NCL2 = ', NCL2,    'HY = ', HY,      'DY1 = ', dyMin, 'Y1 += ', yMinPlus
    WRITE(logflg_pg, '(A, I5, 5(2X, A, F9.5))') '#   NCL2 = ', NCL2,    'HY = ', HY,      'DYC = ', dyMax, 'YC += ', yMaxPlus

    IF(TgFlowFlg) &
    WRITE(logflg_pg, '(A, I12)'         ) '#   MESH_SIZE_TG    = ', NCL1_tg * NCL2 * NCL3
    IF(IoFlowFlg) &
    WRITE(logflg_pg, '(A, I12)'         ) '#   MESH_SIZE_IO    = ', NCL1_io * NCL2 * NCL3
    IF(IoFlowFlg .AND. TgFlowFlg) &
    WRITE(logflg_pg, '(A, I12)'         ) '#   MESH_SIZE_TOTAL = ', (NCL1_io + NCL1_tg) * NCL2 * NCL3
    WRITE(logflg_pg, '(A)') '#****************************************************************************************'
    WRITE(logflg_pg, '(A, 2X, F18.5)')    '#   CONSTANT MEMORY PER CORE (Mb) = ', DBLE(MEMPC_Byte) / 1024.0_WP / 1024.0_WP
    WRITE(logflg_pg, '(A)') '#****************************************************************************************'

    !========================== iNDICATORS ===============================================================

    ! IF(IoFlowFlg .AND. TgFlowFlg)  THEN
    !     CALL date_and_time(DATE = Date, TIME = Time)
    !     fllog= Date(1:4) // '.' // Date(5:8) // '.' // Time(1:4) // '.log'
    !     logflg_tg = 110
    !     OPEN(logflg_tg, FILE = TRIM(FilePath0) // 'hIStory.periodicxz.' //fllog)
    !     logflg_io = 6
    ! ELSE
    !     IF(TgFlowFlg) logflg_tg = 6
    !     IF(IoFlowFlg) logflg_io = 6
    ! END IF


    IF(TgFlowFlg)  THEN
        WRITE(logflg_tg, '(27A10)') &
        '#', '1STEP', '2PhyT', '3CpuT', '4CFL', '5DT', '6MaxDiv', &
        '07U', '08Ut', '09V', '10Vt', '11W', '12Wt', &
        '13UU', '14UUt', '15UV', '16UVt', '17VV', '18VVt', '19WW', '20WWT', &
        '21Cf_L', '22Cf_Lt', '23Cf_U', '24Cf_Ut', '25Umean', '26Umaxx'
    END IF

    IF(IoFlowFlg) THEN
        IF(TgFlowFlg) THEN ! tg + io
            IF(iThermoDynamics == 1) THEN
                WRITE(logflg_io, '(13A)')&
                '# "1STEP","2PhyT", "3CpuT", "4CFL", "5DT", "6DivMFD", "7DivINL", "8DivOUL", ', &
                '"09UU","10UUt","11UV","12UVt", "13VV","14VVt","15WW","16WWT", ', &
                '"17DD","18DDt","19TT","20TTt", "21HH", "22HHt", ', &
                '"23Cf_L","24Cf_Lt","25Cf_U","26Cf_Ut", ', &
                '"27UU","28UUt","29UV","30UVt", "31VV","32VVt","33WW","34WWT", ', &
                '"35DD","36DDt","37TT","38TTt", "39HH", "40HHt", ', &
                '"41Cf_L","42Cf_Lt","43Cf_U","44Cf_Ut", ', &
                '"45UU","46UUt","47UV","48UVt", "49VV","50VVt","51WW","52WWT", ', &
                '"53DD","54DDt","55TT","56TTt", "57HH", "58HHt", ', &
                '"59Cf_L","60Cf_Lt","61Cf_U","62Cf_Ut", ', &
                '"63Gmean", "64Gmaxx","65MassConvs", "66EnegConvs", "67EnegTotal", ', &
                '"68TOTAL_ENERGY_RA", "69TOTAL_ENERGY_FA", "70TOTAL_ENSTROPHY", ', &
                '"71UUV","72UUVt","73UVV","74UVVt","75Fcdrv"'
            ELSE
                WRITE(logflg_io, '(6A)') &
                '# "1STEP","2PhyT", "3CpuT", "4CFL", "5DT", "6DivMFD", "7DivINL", "8DivOUL", ', &
                '"09UU","10UUt","11UV","12UVt", "13VV","14VVt","15WW","16WWT","17Cf_L","18Cf_Lt","19Cf_U","20Cf_Ut", ', &
                '"21UU","22UUt","23UV","24UVt", "25VV","26VVt","27WW","28WWT","29Cf_L","30Cf_Lt","31Cf_U","32Cf_Ut", ', &
                '"33UU","34UUt","35UV","36UVt", "37VV","38VVt","39WW","40WWT","41Cf_L","42Cf_Lt","43Cf_U","44Cf_Ut", ', &
                '"45Gmean", "46Gmaxx","47MassConvs", "48TOTAL_ENERGY_RA", "49TOTAL_ENERGY_FA", "50TOTAL_ENSTROPHY", ', &
                '"51UUV","52UUVt","53UVV","54UVVt","55Fcdrv"'
            END IF
        ELSE ! io only
            IF(iThermoDynamics == 1) THEN
                WRITE(logflg_io, '(A)')'#TITLE = " instantanous xZ -periodic with thermal "'
                WRITE(logflg_io, '(10A)') &
                '#variables = ', &
                '"1STEP", "2PhyT", "3CpuT", "4CFL", "5DT", "6MaxDiv", ', &
                '"07Uc", "08Uct", "09UUc", "10UUct", "11UVC", "12UVCt", ', &
                '"13Cf_B", "14Cf_Bt", "15Cf_T", "16Cf_Tt", ', &
                '"17Gmean", "18Gmaxx", "19MassConvs", ', &
                '"20Tc", "21Tct", "22Trms_B", "23Trmst_Bt","24Trms_T", "25Trmst_Tt", ', &
                '"26Qwflux_B", "27Qwflux_Bt", "28Qwflux_T", "29Qwflux_Tt", ', &
                '"30Tbulk", "31Tmaxx", "32EnegConvs", "33EnegTotal", "34QWRatio", "35QWRatio_t", ', &
                '"36TOTAL_ENERGY_RA", "37TOTAL_ENERGY_FA", "38TOTAL_ENSTROPHY" ', &
                '"39UUV","40UUVt","41UVV","42UVVt","43Fcdrv"'

                WRITE(logflg_io, '(A)')'#ZONE T = "tracking"'
            ELSE
                WRITE(logflg_io, '(A)')'#TITLE = " instantanous xz-periodic without thermal "'
                WRITE(logflg_io, '(7A)')'#variables = ', &
                '"1STEP", "2PhyT", "3CpuT", "4CFL", "5DT", "6MaxDiv", ', &
                '"07Uc", "08Uct", "09UUc", "10UUct", "11UVC", "12UVCt", ', &
                '"13Cf_B", "14Cf_Bt", "15Cf_T", "16Cf_Tt", ', &
                '"17Gmean", "18Gmaxx", "19MassConvs", ', &
                '"20TOTAL_ENERGY_RA","21tke_energy_ave", ', &
                '"22TOTAL_ENSTROPHY","23enstrophy_ave","24UUV","25UUVt","26UVV","27UVVt","28Fcdrv"'
                WRITE(logflg_io, '(A)')'#ZONE T = "tracking"'


            END IF
        END IF

    END IF

    RETURN
END SUBROUTINE


!******************************************************************************************************************************
SUBROUTINE PP_MONITOR_TG
    USE flow_info
    USE init_info
    USE mesh_info
    USE postprocess_info
    USE WRT_INFO
    IMPLICIT NONE

    REAL(WP) :: CF_XZT(4)  !1 = LW, 2 = T_LW, 3 = UW,  4 = T_UW
    REAL(WP) :: CF_XZT_WORK(4)
    REAL(WP) :: DUDY_xz, DUDY_xzt
    REAL(WP) :: UU_XZT(8), U1_XZT(6)
    REAL(WP) :: UU_XZT_WORK(8), U1_XZT_WORK(6)
    INTEGER(4) :: JID, J, JJ, JJSCN

    !CALL MPI_BARRIER(ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(CPUTIME_tmp, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

    !========CF ============================================
    CF_XZT(:) = 0.0_WP
    IF(iCase /= IPIPEC) THEN
        IF(MYID == 0) THEN
            JID = 1
            JJ = JCL2G(JID)
            DUDY_xz = (U1xzL_tg(JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
            !DUDY_xz = (Q_tg(NCL1_tg / 2, JID, NCL3 / 2, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
            CF_XZT(1) = DUDY_xz  * 2.0_WP / REN

            DUDY_xzt = (U1xztL_tg(JID, 1) - 0.0_WP) * DYFI(JID) * 2.0_WP
            CF_XZT(2) = DUDY_xzt * 2.0_WP / REN
        END IF
    END IF

    IF(MYID == NPSLV) THEN
        JID = N2DO(MYID)
        JJ = JCL2G(JID)
        DUDY_xz = (U1xzL_tg(JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
        !DUDY_xz = (Q_tg(NCL1_tg / 2, JID, NCL3 / 2, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
        CF_XZT(3) = DUDY_xz  * 2.0_WP / REN

        DUDY_xzt = (U1xztL_tg(JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
        CF_XZT(4) = DUDY_xzt * 2.0_WP / REN
    END IF

    !=========UU, VV,WW,UV ==================================
    IF(iCase == iPIPEC) THEN
        JJSCN = 1
    ELSE
        JJSCN = NCL2 / 2
    END IF

    U1_XZT(:) = 0.0_WP
    U1_XZT_WORK(:) = 0.0_WP
    UU_XZT(:) = 0.0_WP
    UU_XZT_WORK(:) = 0.0_WP

    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        IF(JJ == JJSCN) THEN
            U1_XZT(1) = U1xzL_tg(J, 1)
            U1_XZT(3) = U1xzL_tg(J, 2)
            U1_XZT(5) = U1xzL_tg(J, 3)

            U1_XZT(2) = U1xztL_tg(J, 1)
            U1_XZT(4) = U1xztL_tg(J, 2)
            U1_XZT(6) = U1xztL_tg(J, 3)

            UU_XZT(1) = U2xzL_tg(J, 1)  - U1xzL_tg(J, 1) * U1xzL_tg(J, 1)!UUXZ
            UU_XZT(2) = U2xztL_tg(J, 1) - U1xztL_tg(J, 1) * U1xztL_tg(J, 1)

            UU_XZT(3) = U2xzL_tg(J, 2)  - U1xzL_tg(J, 1) * U1xzL_tg(J, 2)!UVXZ
            UU_XZT(4) = U2xztL_tg(J, 2) - U1xztL_tg(J, 1) * U1xztL_tg(J, 2)

            UU_XZT(5) = U2xzL_tg(J, 4)  - U1xzL_tg(J, 2) * U1xzL_tg(J, 2)!VVXZ
            UU_XZT(6) = U2xztL_tg(J, 4) - U1xztL_tg(J, 2) * U1xztL_tg(J, 2)

            UU_XZT(7) = U2xzL_tg(J,6)  - U1xzL_tg(J, 3) * U1xzL_tg(J, 3)!WWXZ
            UU_XZT(8) = U2xztL_tg(J,6) - U1xztL_tg(J, 3) * U1xztL_tg(J, 3)
        END IF
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(CF_XZT, CF_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(UU_XZT, UU_XZT_WORK, 8, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(U1_XZT, U1_XZT_WORK, 6, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)


    !================PRINT DATA ON SCREEN FOR MONITORING==============
    IF (MYID == 0) THEN
        WRITE(logflg_tg, '(1I14, 1F10.3, 2F7.3, 1F10.5, 1ES18.10, 18ES13.5, 2F17.9)') &
        ITERG, PhyTIME, CPUTIME(1), CFLMM * DT, DT, MAXDIVGV_TG(2), &
        U1_XZT_WORK(1),U1_XZT_WORK(2), &
        U1_XZT_WORK(3),U1_XZT_WORK(4), &
        U1_XZT_WORK(5),U1_XZT_WORK(6), &
        UU_XZT_WORK(1),UU_XZT_WORK(2), &
        UU_XZT_WORK(3),UU_XZT_WORK(4), &
        UU_XZT_WORK(5),UU_XZT_WORK(6), &
        UU_XZT_WORK(7),UU_XZT_WORK(8), &
        CF_XZT_WORK(1),CF_XZT_WORK(2), &
        CF_XZT_WORK(3),CF_XZT_WORK(4), &
        U1mean_WORK_TG, U1maxx_WORK_TG
    ENDIF

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE PP_MONITOR_Xperiodic_io
    USE flow_info
    USE init_info
    USE mesh_info
    USE postprocess_info
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE
    !=== Top and bottoM ==================
    !REAL(WP) :: CFw_XZT(4), CFw_XZT_WORK(4)  !1 = LW, 2 = T_LW, 3 = UW,  4 = T_UW
    REAL(WP) :: Cfw_io(4)
    REAL(WP) :: Trmsw_XZT(4), Trmsw_XZT_WORK(4)
    REAL(WP) :: WHw(4), Hw, Tw, WH_ratio(2)

    !========centrE ======================
    REAL(WP) :: UGc_XZT(4), UGc_XZT_WORK(4)
    REAL(WP) :: UGU_XZT(4),UGU_XZT_WORK(4)
    REAL(WP) :: Uc_XZT(2), Uc_XZT_WORK(2)
    REAL(WP) :: Tc_XZT(2), Tc_XZT_WORK(2)
    INTEGER(4) :: I, K, IP, KP, JP, IM, KM, JM
    INTEGER(4) :: JID, J, JJ, JJSCN, N, JJA, JJB, N1, N2
    INTEGER(4) :: H, M, LMNH
    REAL(WP) :: TMP, TMP1, TMP2
    REAL(WP) :: TOTAL_ENERGY_FA, TOTAL_ENERGY_FA_WORK
    REAL(WP) :: TOTAL_ENERGY_RA, TOTAL_ENERGY_RA_WORK
    REAL(WP) :: TOTAL_ENERGY_M2, TOTAL_ENERGY_M2_WORK
    REAL(WP) :: TOTAL_ENSTROPHY, TOTAL_ENSTROPHY_WORK
    REAL(WP) :: TOTAL_ENSPHY_M2, TOTAL_ENSPHY_M2_WORK

    !REAL(WP) :: aveg_energy_RA, aveg_energy_RA_WORK
    !REAL(WP) :: aveg_enstrophy, aveg_enstrophy_WORK


    !CALL MPI_BARRIER(ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(CPUTIME_tmp, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

    CALL CHK_MassConsv_io
    IF(iThermoDynamics == 1) CALL CHK_EnegConsv_io

    Cfw_io(:)      = 0.0_WP
    WHw(:)         = 0.0_WP
    WH_ratio(:)    = 0.0_WP
    !============================ info on the walL =========================================
    CALL PP_wall_thermal_shear(flgxz)
    Cfw_io(1) = 2.0_WP * Tauw_io(1)
    Cfw_io(3) = 2.0_WP * Tauw_io(2)

    IF(iThermoDynamics == 1) THEN
        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
            WHw(1) = Twal(1) * T0 ! Tw_D
            WHw(3) = Twal(2) * T0 ! Tw_D
        END IF
        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            WHw(1) = Qw(1) * D0 * U0 * CP0 * T0
            WHw(3) = Qw(2) * D0 * U0 * CP0 * T0
        END IF
    END IF

    IF(PhyTIME > tRunAve1) THEN
        CALL PP_wall_thermal_shear(flgxzt)
        Cfw_io(2) = 2.0_WP * Tauw_io(1)
        Cfw_io(4) = 2.0_WP * Tauw_io(2)
        IF(iThermoDynamics == 1) THEN
            IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
                WHw(2) = Twal(1) * T0 ! Tw_D
                WHw(4) = Twal(2) * T0 ! Tw_D
            END IF
            IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
                WHw(2) = Qw(1) * D0 * U0 * CP0 * T0
                WHw(4) = Qw(2) * D0 * U0 * CP0 * T0
            END IF
        END IF
    ELSE
        Cfw_io(2) = Cfw_io(1)
        Cfw_io(4) = Cfw_io(3)
        IF(iThermoDynamics == 1) THEN
            WHw(2) = WHw(1)
            WHw(4) = WHw(3)
        ENDIF
    END IF

    !================================================================================
    Trmsw_XZT(:)   = 0.0_WP
    Trmsw_XZT_WORK(:) = 0.0_WP
    IF( (iThermoDynamics == 1) .AND. ( (iCase /= IPIPEC .AND. MYID == 0) .OR. MYID == NPSLV ) ) THEN
        IF(iCase /= IPIPEC .AND. MYID == 0) THEN
            JID = 1
            JJ  = JCL2G(JID)
            N1  = 1
            N2  = 2
        END IF
        IF(MYID == NPSLV) THEN
            JID = N2DO(MYID)
            JJ  = JCL2G(JID)
            N1  = 3
            N2  = 4
        END IF

        !================== TrmS ======================
        DO N = N1, N2
            TMP2 = T2xzL_io(JID)
            TMP1 = T1xzL_io(JID)
            IF(N == N2) THEN
                IF(PhyTIME > tRunAve1) THEN
                    TMP2 = T2xztL_io(JID)
                    TMP1 = T1xztL_io(JID)

                ELSE
                    Trmsw_XZT(N) = Trmsw_XZT(N1)
                    EXIT
                END IF
            END IF
            Trmsw_XZT(N) = DSQRT(DABS(TMP2 - TMP1 * TMP1))
        END DO

    END IF

    !============================ info on y-centre =========================================
    IF(iCase == iPIPEC) THEN
        JJSCN = 1
    ELSE
        JJSCN = NCL2 / 2
    END IF

    UGc_XZT(:) = 0.0_WP
    UGU_XZT(:) = 0.0_WP
    Uc_XZT(:)   = 0.0_WP
    Tc_XZT(:)   = 0.0_WP

    UGc_XZT_WORK(:) = 0.0_WP
    UGU_XZT_WORK(:) = 0.0_WP
    Uc_XZT_WORK(:)   = 0.0_WP
    Tc_XZT_WORK(:)   = 0.0_WP
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        IF(JJ == JJSCN) THEN
            UGc_XZT(1) = UGxzL_io(J, 1) - G1xzL_io(J, 1) * G1xzL_io(J, 1) / D1xzL_io(J)!UUXZ
            UGc_XZT(3) = UGxzL_io(J, 2) - G1xzL_io(J, 1) * G1xzL_io(J, 2) / D1xzL_io(J)!UVXZ
            Uc_XZT(1)  = G1xzL_io(J, 1) / D1xzL_io(J)

            M = 1
            N = 1
            H = 2
            LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
            UGU_XZT(1) = UGUxzL_io(J, LMNH) / D1xzL_io(J) &
            - G1xzL_io(J, M) / D1xzL_io(J) * UGxzL_io(J, (N * (7-N)) / 2 + H- 3) / D1xzL_io(J) &
            - G1xzL_io(J, N) / D1xzL_io(J) * UGxzL_io(J, (M * (7-M)) / 2 + H- 3) / D1xzL_io(J) &
            - G1xzL_io(J,H) / D1xzL_io(J) * UGxzL_io(J, (M * (7-M)) / 2 + N - 3) / D1xzL_io(J) &
            + 2.0_WP * G1xzL_io(J, M) / D1xzL_io(J) * &
            G1xzL_io(J, N) / D1xzL_io(J) * &
            G1xzL_io(J,H) / D1xzL_io(J)

            M = 1
            N = 2
            H = 2
            LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
            UGU_XZT(3) = UGUxzL_io(J, LMNH) / D1xzL_io(J) &
            - G1xzL_io(J, M) / D1xzL_io(J) * UGxzL_io(J, (N * (7-N)) / 2 + H- 3) / D1xzL_io(J) &
            - G1xzL_io(J, N) / D1xzL_io(J) * UGxzL_io(J, (M * (7-M)) / 2 + H- 3) / D1xzL_io(J) &
            - G1xzL_io(J,H) / D1xzL_io(J) * UGxzL_io(J, (M * (7-M)) / 2 + N - 3) / D1xzL_io(J) &
            + 2.0_WP * G1xzL_io(J, M) / D1xzL_io(J) * &
            G1xzL_io(J, N) / D1xzL_io(J) * &
            G1xzL_io(J,H) / D1xzL_io(J)

            IF(iThermoDynamics == 1) Tc_XZT(1) = T1xzL_io(J)
        END IF
    END DO

    IF(PhyTIME > tRunAve1) THEN
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            IF(JJ == JJSCN) THEN
                UGc_XZT(2) = UGxztL_io(J, 1) - G1xztL_io(J, 1) * G1xztL_io(J, 1) / D1xztL_io(J)
                UGc_XZT(4) = UGxztL_io(J, 2) - G1xztL_io(J, 1) * G1xztL_io(J, 2) / D1xztL_io(J)
                Uc_XZT(2) = G1xztL_io(J, 1) / D1xztL_io(J)

                M = 1
                N = 1
                H = 2
                LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
                UGU_XZT(2) = UGUxztL_io(J, LMNH) / D1xztL_io(J) &
                - G1xztL_io(J, M) / D1xztL_io(J) * UGxztL_io(J, (N * (7-N)) / 2 + H- 3) / D1xztL_io(J) &
                - G1xztL_io(J, N) / D1xztL_io(J) * UGxztL_io(J, (M * (7-M)) / 2 + H- 3) / D1xztL_io(J) &
                - G1xztL_io(J,H) / D1xztL_io(J) * UGxztL_io(J, (M * (7-M)) / 2 + N - 3) / D1xztL_io(J) &
                + 2.0_WP * G1xztL_io(J, M) / D1xztL_io(J) * &
                G1xztL_io(J, N) / D1xztL_io(J) * &
                G1xztL_io(J,H) / D1xztL_io(J)

                M = 1
                N = 2
                H = 2
                LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
                UGU_XZT(4) = UGUxztL_io(J, LMNH) / D1xztL_io(J) &
                - G1xztL_io(J, M) / D1xztL_io(J) * UGxztL_io(J, (N * (7-N)) / 2 + H- 3) / D1xztL_io(J) &
                - G1xztL_io(J, N) / D1xztL_io(J) * UGxztL_io(J, (M * (7-M)) / 2 + H- 3) / D1xztL_io(J) &
                - G1xztL_io(J,H) / D1xztL_io(J) * UGxztL_io(J, (M * (7-M)) / 2 + N - 3) / D1xztL_io(J) &
                + 2.0_WP * G1xztL_io(J, M) / D1xztL_io(J) * &
                G1xztL_io(J, N) / D1xztL_io(J) * &
                G1xztL_io(J,H) / D1xztL_io(J)

                IF(iThermoDynamics == 1) Tc_XZT(2) = T1xztL_io(J)
            END IF
        END DO
    ELSE
        UGc_XZT(2) = UGc_XZT(1)
        UGc_XZT(4) = UGc_XZT(3)

        UGU_XZT(2) = UGU_XZT(1)
        UGU_XZT(4) = UGU_XZT(3)

        Uc_XZT(2) = Uc_XZT(1)
        Tc_XZT(2) = Tc_XZT(1)
    END IF

    !==== Total energy and enstrophY ==============
    TOTAL_ENERGY_RA = 0.0_WP
    TOTAL_ENERGY_FA = 0.0_WP
    TOTAL_ENSTROPHY = 0.0_WP
    DO J = 1, N2DO(MYID)
        ! check, Jundi version without the geometry info.
        JJ = JCL2G(J)
        !TOTAL_ENERGY_RA =  TOTAL_ENERGY_RA + ( &
        !                    U2xzL_io(J, 1) - U1xzL_io(J, 1) * U1xzL_io(J, 1) + &
        !                    U2xzL_io(J, 2) - U1xzL_io(J, 2) * U1xzL_io(J, 2) + &
        !                    U2xzL_io(J, 3) - U1xzL_io(J, 3) * U1xzL_io(J, 3) ) * 0.5_WP * (1.0W_P / DYFI(JJ))

        TOTAL_ENERGY_RA =  TOTAL_ENERGY_RA + ( &
        U2xzL_io(J, 1)  + &
        U2xzL_io(J, 2)  + &
        U2xzL_io(J, 3)  ) * 0.5_WP * (1.0_WP / DYFI(JJ))

        TOTAL_ENERGY_FA =  TOTAL_ENERGY_FA + ( &
        UGxzL_io(J, 1) - G1xzL_io(J, 1) * G1xzL_io(J, 1) / D1xzL_io(J) + &
        UGxzL_io(J, 2) - G1xzL_io(J, 2) * G1xzL_io(J, 2) / D1xzL_io(J) + &
        UGxzL_io(J, 3) - G1xzL_io(J, 3) * G1xzL_io(J, 3) / D1xzL_io(J) ) * 0.5_WP * (1.0_WP / DYFI(JJ))



        TOTAL_ENSTROPHY =  TOTAL_ENSTROPHY + ( &
        (DVDL2xzL_io(J, (3- 1) * NDV + 2, (3- 1) * NDV + 2) - 2.0_WP * &
        DVDL2xzL_io(J, (3- 1) * NDV + 2, (2 - 1) * NDV +3) +         &
        DVDL2xzL_io(J, (2 - 1) * NDV +3, (2 - 1) * NDV +3) ) -       &
        (DVDL1xzL_io(J, 3, 2) * DVDL1xzL_io(J, 3, 2) - 2.0_WP * &
        DVDL1xzL_io(J, 3, 2) * DVDL1xzL_io(J, 2, 3)  +         &
        DVDL1xzL_io(J, 2, 3) * DVDL1xzL_io(J, 2, 3) ) +        &
        (DVDL2xzL_io(J, (1- 1) * NDV +3, (1- 1) * NDV +3) - 2.0_WP * &
        DVDL2xzL_io(J, (1- 1) * NDV +3, (3- 1) * NDV + 1) +         &
        DVDL2xzL_io(J, (3- 1) * NDV + 1, (3- 1) * NDV + 1) ) -       &
        (DVDL1xzL_io(J, 1, 3) * DVDL1xzL_io(J, 1, 3) - 2.0_WP * &
        DVDL1xzL_io(J, 1, 3) * DVDL1xzL_io(J, 3, 1) +          &
        DVDL1xzL_io(J, 3, 1) * DVDL1xzL_io(J, 3, 1) ) +        &
        (DVDL2xzL_io(J, (2 - 1) * NDV + 1, (2 - 1) * NDV + 1) - 2.0_WP * &
        DVDL2xzL_io(J, (2 - 1) * NDV + 1, (1- 1) * NDV + 2) +         &
        DVDL2xzL_io(J, (1- 1) * NDV + 2, (1- 1) * NDV + 2) ) -       &
        (DVDL1xzL_io(J, 2, 1) * DVDL1xzL_io(J, 2, 1) - 2.0_WP * &
        DVDL1xzL_io(J, 2, 1) * DVDL1xzL_io(J, 1, 2) +          &
        DVDL1xzL_io(J, 1, 2) * DVDL1xzL_io(J, 1, 2) ) ) * 0.5_WP * (1.0_WP / DYFI(JJ))
    END DO


    ! Below IS for TGV
    TOTAL_ENERGY_M2 = 0.0_WP
    DO I = 1, NCL1_io
        IP = IPV_io(I)
        DO K = 1, NCL3
            KP = KPV(K)
            DO J = 1, N2DO(MYID)
                JP = JLPV(J)
                JJ = JCL2G(J)
                TOTAL_ENERGY_M2 = TOTAL_ENERGY_M2 + ( &
                ((Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * 0.5_WP)**2 + &
                ((Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * 0.5_WP)**2 + &
                ((Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * 0.5_WP)**2 ) * 0.5_WP * DX * DZ/ DYFI(JJ)
            END DO
        END DO
    END DO


    TOTAL_ENSPHY_M2 = 0.0_WP
    DO I = 1, NCL1_io
        IP = IPV_io(I)
        IM = IMV_io(I)
        DO K = 1, NCL3
            KP = KPV(K)
            KM = KMV(K)
            DO J = 1, N2DO(MYID)
                JP = JLPV(J)
                JM = JLMV(J)
                JJ = JCL2G(J)

                TOTAL_ENSPHY_M2 = TOTAL_ENSPHY_M2 + ( &
                ((( ( (Q_io(I, JP, K, 3) + Q_io(I, JP, KP, 3)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * 0.5_WP ) + &
                ( (Q_io(I, JM, K, 3) + Q_io(I, JM, KP, 3)) * 0.5_WP ) ) * 0.5_WP) * DYFI(JJ) - &
                (( ( (Q_io(I, J, KP, 2) + Q_io(I, JP, KP, 2)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * 0.5_WP ) + &
                ( (Q_io(I, J, KM, 2) + Q_io(I, JP, KM, 2)) * 0.5_WP ) ) * 0.5_WP) * DZI )**2 +  &
                ((( ( (Q_io(I, J, KP, 1) + Q_io(IP, J, KP, 1)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * 0.5_WP ) + &
                ( (Q_io(I, J, KM, 1) + Q_io(IP, J, KM, 1)) * 0.5_WP ) ) * 0.5_WP) * DZI - &
                (( ( (Q_io(IP, J, K, 3) + Q_io(IP, J, KP, 3)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * 0.5_WP ) + &
                ( (Q_io(IM, J, K, 3) + Q_io(IM, J, KP, 3)) * 0.5_WP ) ) * 0.5_WP) * DXI )**2 +  &
                ((( ( (Q_io(IP, J, K, 2) + Q_io(IP, JP, K, 2)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * 0.5_WP ) + &
                ( (Q_io(IM, J, K, 2) + Q_io(IM, JP, K, 2)) * 0.5_WP ) ) * 0.5_WP) * DXI - &
                (( ( (Q_io(I, JP, K, 1) + Q_io(IP, JP, K, 1)) * 0.5_WP ) + &
                ( (Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * 0.5_WP ) ) * 0.5_WP - &
                ( ( (Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * 0.5_WP ) + &
                ( (Q_io(I, JM, K, 1) + Q_io(IP, JM, K, 1)) * 0.5_WP ) ) * 0.5_WP) * DYFI(JJ) )**2 ) * DX * DZ/ DYFI(JJ)
            END DO
        END DO
    END DO

    !============================ info on y-centre =========================================
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(UGc_XZT, UGc_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(UGU_XZT, UGU_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(Uc_XZT,  Uc_XZT_WORK,  2, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    CALL MPI_ALLREDUCE(TOTAL_ENERGY_RA,  TOTAL_ENERGY_RA_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(TOTAL_ENERGY_FA,  TOTAL_ENERGY_FA_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(TOTAL_ENSTROPHY,  TOTAL_ENSTROPHY_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    CALL MPI_ALLREDUCE(TOTAL_ENERGY_M2,  TOTAL_ENERGY_M2_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(TOTAL_ENSPHY_M2,  TOTAL_ENSPHY_M2_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    TOTAL_ENERGY_RA_WORK  = TOTAL_ENERGY_RA_WORK / ((HYT - HYB))
    TOTAL_ENERGY_FA_WORK  = TOTAL_ENERGY_FA_WORK / ((HYT - HYB))
    TOTAL_ENSTROPHY_WORK  = TOTAL_ENSTROPHY_WORK / ((HYT - HYB))

    TOTAL_ENERGY_M2_WORK  = TOTAL_ENERGY_M2_WORK / (HX_io * HZ* (HYT - HYB))
    TOTAL_ENSPHY_M2_WORK  = TOTAL_ENSPHY_M2_WORK / (HX_io * HZ* (HYT - HYB))

    !IF(MYID == 0) THEN
    !    WRITE(*, *) '##',PhyTIME, TOTAL_ENERGY_RA_WORK, TOTAL_ENERGY_M2_WORK, TOTAL_ENSTROPHY_WORK, TOTAL_ENSPHY_M2_WORK
    ! END IF
    IF(iThermoDynamics == 1) THEN
        CALL MPI_ALLREDUCE(Trmsw_XZT, Trmsw_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(Tc_XZT,    Tc_XZT_WORK,    2, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    END IF

    !=====================ON THE TOP WALL =================================
    IF(MYID == 0) THEN
        IF(iThermoDynamics == 1) THEN
            IF(iCase /= IPIPEC) THEN
                WH_ratio(1) = ( DABS(WHw(1)) - DABS(WHw(3)) ) / ( (DABS(WHw(1)) + DABS(WHw(3))) )
                WH_ratio(2) = ( DABS(WHw(2)) - DABS(WHw(4)) ) / ( (DABS(WHw(2)) + DABS(WHw(4))) )
            END IF

            WRITE(logflg_io, '(1I14, 1F12.5, 3F9.5, 1ES11.3, 10ES13.5, 2F9.5, 1ES13.5, 10ES13.5, 2F9.5, 12ES13.5)') &
            ITERG, PhyTIME, CPUTIME(1), CFLMM * DT, DT, MAXDIVGV_io(1), &
            Uc_XZT_WORK(1:2), UGc_XZT_WORK(1:4),   Cfw_io(1:4), G1BULK_WORK_io, G1maxx_WORK_io, CHK_MASS_CONSV0, &
            Tc_XZT_WORK(1:2), Trmsw_XZT_WORK(1:4), WHw(1:4),    T1BULK_WORK_io, T1maxx_WORK_io, CHK_ENEG_CONSV0, &
            CHK_ENEG_TOTAL, WH_ratio(2), TOTAL_ENERGY_RA_WORK, TOTAL_ENERGY_FA_WORK, TOTAL_ENSTROPHY_WORK, &
            UGU_XZT_WORK(1:4), FcDrv_io

        ELSE
            WRITE(logflg_io, '(1I14, 1F12.5, 3F9.5, 1ES11.3, 10ES13.5, 2F9.5, 10ES13.5)') &
            ITERG, PhyTIME, CPUTIME(1), CFLMM * DT, DT, MAXDIVGV_io(1), &
            Uc_XZT_WORK(1:2), UGc_XZT_WORK(1:4),  Cfw_io(1:4), G1BULK_WORK_io, G1maxx_WORK_io, CHK_MASS_CONSV0, &
            TOTAL_ENERGY_RA_WORK, TOTAL_ENERGY_M2_WORK, TOTAL_ENSTROPHY_WORK, TOTAL_ENSPHY_M2_WORK, &
            UGU_XZT_WORK(1:4), FcDrv_io

        END IF
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_MONITOR_nonXperiodic_io
    !>  @note: three point at x direction.
    USE flow_info
    USE init_info
    USE mesh_info
    USE postprocess_info
    USE WRT_INFO
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4) :: IID(NDV) ! NDV = 3 ONLY SHOW NUMBERS, NOT MEANING.
    REAL(WP) :: CF_ZT(4, NDV)  !1 = LW, 2 = T_LW, 3 = UW,  4 = T_UW
    REAL(WP) :: CF_ZT_WORK(4, NDV)
    REAL(WP) :: UG_ZT(8, NDV)
    REAL(WP) :: UG_ZT_WORK(8, NDV)

    REAL(WP) :: D2_ZT(2, NDV)
    REAL(WP) :: H2_ZT(2, NDV)
    REAL(WP) :: T2_ZT(2, NDV)
    REAL(WP) :: K2_ZT(2, NDV)
    REAL(WP) :: M2_ZT(2, NDV)

    REAL(WP) :: DUDY_z, DUDY_zt
    INTEGER(4) :: N, JID, J, JJ, JJSCN


    !CALL MPI_BARRIER(ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(CPUTIME_tmp, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

    !CALL CHK_MassConsv_io
    IF(iThermoDynamics == 1) CALL CHK_EnegConsv_io


    DO N = 1, NDV
        IID(N) = NCL1_io / (NDV + 1) * N
    END DO

    !========CF ============================================
    CF_ZT(:, :) = 0.0_WP

    IF(iCase /= IPIPEC) THEN
        IF(MYID == 0) THEN
            JID = 1
            JJ = JCL2G(JID)
            DO N = 1, NDV
                DUDY_z = (U1zL_io(IID(N), JID, 1) - 0.0_WP) * DYFI(JID) * 2.0_WP
                !DUDY_z = (Q_io(IID(N), JID, NCL3 / 2, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
                CF_ZT(1, N) = DUDY_z  * 2.0_WP / REN

                DUDY_zt = (U1ztL_io(IID(N), JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
                CF_ZT(2, N) = DUDY_zt * 2.0_WP / REN
            END DO
        END IF
    END IF

    IF(MYID == NPSLV) THEN
        JID = N2DO(MYID)
        JJ = JCL2G(JID)
        DO N = 1, NDV
            DUDY_z = (U1zL_io(IID(N), JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
            !DUDY_z = (Q_io(IID(N), JID, NCL3 / 2, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
            CF_ZT(3, N) = DUDY_z  * 2.0_WP / REN
            DUDY_zt = (U1ztL_io(IID(N), JID, 1) - 0.0_WP) * DYFI(JJ) * 2.0_WP
            CF_ZT(4, N) = DUDY_zt * 2.0_WP / REN
        END DO
    END IF

    !=========UU, VV,WW,UV ==================================
    IF(iCase == iPIPEC) THEN
        JJSCN = 1
    ELSE
        JJSCN = NCL2 / 2
    END IF

    UG_ZT(:, :) = 0.0_WP
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        IF(JJ == JJSCN) THEN
            DO N = 1, NDV
                UG_ZT(1, N) = UGzL_io(IID(N), J, 1) - G1zL_io(IID(N), J, 1) * G1zL_io(IID(N), J, 1) / D1zL_io(IID(N), J)!UUXZ
                UG_ZT(3, N) = UGzL_io(IID(N), J, 2) - G1zL_io(IID(N), J, 1) * G1zL_io(IID(N), J, 2) / D1zL_io(IID(N), J)!UVXZ
                UG_ZT(5, N) = UGzL_io(IID(N), J, 4) - G1zL_io(IID(N), J, 2) * G1zL_io(IID(N), J, 2) / D1zL_io(IID(N), J)!VVXZ
                UG_ZT(7, N) = UGzL_io(IID(N), J,6) - G1zL_io(IID(N), J, 3) * G1zL_io(IID(N), J, 3) / D1zL_io(IID(N), J)!WWXZ
            END DO
        END IF
    END DO

    IF(PhyTIME > tRunAve1) THEN
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            IF(JJ == JJSCN) THEN
                DO N = 1, NDV
                    UG_ZT(2, N) = UGztL_io(IID(N), J, 1) - G1ztL_io(IID(N), J, 1) * G1ztL_io(IID(N), J, 1) / D1ztL_io(IID(N), J)
                    UG_ZT(4, N) = UGztL_io(IID(N), J, 2) - G1ztL_io(IID(N), J, 1) * G1ztL_io(IID(N), J, 2) / D1ztL_io(IID(N), J)
                    UG_ZT(6, N) = UGztL_io(IID(N), J, 4) - G1ztL_io(IID(N), J, 2) * G1ztL_io(IID(N), J, 2) / D1ztL_io(IID(N), J)
                    UG_ZT(8, N) = UGztL_io(IID(N), J,6) - G1ztL_io(IID(N), J, 3) * G1ztL_io(IID(N), J, 3) / D1ztL_io(IID(N), J)
                END DO
            END IF
        END DO
    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(CF_ZT, CF_ZT_WORK, 4* NDV, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(UG_ZT, UG_ZT_WORK, 8* NDV, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    !=====================ON THE TOP WALL =================================
    IF(MYID == 0) THEN


        IF(iThermoDynamics == 1) THEN

            J = N2DO(MYID)
            DO N = 1, NDV
                D2_ZT(1, N) = D2zL_io(IID(N), J) - D1zL_io(IID(N), J) * D1zL_io(IID(N), J)
                T2_ZT(1, N) = T2zL_io(IID(N), J) - T1zL_io(IID(N), J) * T1zL_io(IID(N), J)
                H2_ZT(1, N) = H2zL_io(IID(N), J) - H1zL_io(IID(N), J) * H1zL_io(IID(N), J)
                K2_ZT(1, N) = K2zL_io(IID(N), J) - K1zL_io(IID(N), J) * K1zL_io(IID(N), J)
                M2_ZT(1, N) = M2zL_io(IID(N), J) - M1zL_io(IID(N), J) * M1zL_io(IID(N), J)

                D2_ZT(2, N) = D2zTL_io(IID(N), J) - D1zTL_io(IID(N), J) * D1zTL_io(IID(N), J)
                T2_ZT(2, N) = T2zTL_io(IID(N), J) - T1zTL_io(IID(N), J) * T1zTL_io(IID(N), J)
                H2_ZT(2, N) = H2zTL_io(IID(N), J) - H1zTL_io(IID(N), J) * H1zTL_io(IID(N), J)
                K2_ZT(2, N) = K2zTL_io(IID(N), J) - K1zTL_io(IID(N), J) * K1zTL_io(IID(N), J)
                M2_ZT(2, N) = M2zTL_io(IID(N), J) - M1zTL_io(IID(N), J) * M1zTL_io(IID(N), J)

            END DO
            WRITE(logflg_io, '(1I14, 1F10.3, 2F7.3, 1F10.5, 3ES10.2, 54ES13.5, 2F13.5, 3ES13.5)') &
            ITERG, PhyTIME, CPUTIME(1), CFLMM * DT,  DT, MAXDIVGV_io(1:3), &
            ( UG_ZT_WORK(1, N),UG_ZT_WORK(2, N),UG_ZT_WORK(3, N),UG_ZT_WORK(4, N), &
            UG_ZT_WORK(5, N),UG_ZT_WORK(6, N),UG_ZT_WORK(7, N),UG_ZT_WORK(8, N), &
            CF_ZT_WORK(1, N),CF_ZT_WORK(2, N),CF_ZT_WORK(3, N),CF_ZT_WORK(4, N), &
            D2_ZT(1, N), D2_ZT(2, N), T2_ZT(1, N), T2_ZT(2, N), H2_ZT(1, N), H2_ZT(2, N), &
            N = 1, NDV), G1BULK_WORK_io, G1maxx_WORK_io, CHK_MASS_CONSV0, CHK_ENEG_CONSV0, &
            CHK_ENEG_TOTAL

        ELSE
            WRITE(logflg_io, '(1I14, 1F10.3, 2F7.3, 1F10.5, 3ES10.2, 36ES13.5, 1F13.5, 2ES13.5)') &
            ITERG, PhyTIME, CPUTIME(1), CFLMM * DT,  DT, MAXDIVGV_io(1:3), &
            ( UG_ZT_WORK(1, N),UG_ZT_WORK(2, N),UG_ZT_WORK(3, N),UG_ZT_WORK(4, N), &
            UG_ZT_WORK(5, N),UG_ZT_WORK(6, N),UG_ZT_WORK(7, N),UG_ZT_WORK(8, N), &
            CF_ZT_WORK(1, N),CF_ZT_WORK(2, N),CF_ZT_WORK(3, N),CF_ZT_WORK(4, N), &
            N = 1, NDV), G1BULK_WORK_io, G1maxx_WORK_io, CHK_MASS_CONSV0

        END IF
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_wall_thermal_shear(flg_xzt)
    USE flow_info
    USE init_info
    USE mesh_info
    USE postprocess_info
    USE thermal_info
    USE WRT_INFO
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: FLG_XZT
    INTEGER(4) :: JID, JJ, JJA, N1, NSZ, J
    REAL(WP) :: DUMMY(7, 2)
    REAL(WP) :: DUMMY_WORK(7, 2)
    REAL(WP) :: DUDY
    REAL(WP) :: DEN, YL, VIS
    REAL(WP) :: Tau_diff, Tau_avag


    DEN = 1.0_WP
    VIS = 1.0_WP
    YL  = 1.0_WP
    IF(iThermoDynamics == 1) then 
        CALL PP_Wall_thermal_properties(flg_xzt)

    !================ ARithmetIC mean of DENSITY ==================
        DEN = 0.0_WP
        VIS = 0.0_WP
        YL  = 0.0_WP
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            YL = YL + 1.0_WP / DYFI(JJ)
            IF(flg_xzT == flgxz) THEN
                DEN = DEN + D1xzL_io(J) / DYFI(JJ)
                VIS = VIS + M1xzL_io(J) / DYFI(JJ)
            ELSE IF (flg_xzT == flgxzt)  THEN
                DEN = DEN + D1xztL_io(J) / DYFI(JJ)
                VIS = VIS + M1xztL_io(J) / DYFI(JJ)
            ELSE
            END IF
        END DO
    END IF

    !================ ARithmetIC mean of DENSITY ==================
    Utaw_io(1:2) = 0.0_WP
    Tauw_io(1:2) = 0.0_WP
    Utaw_D_io(1:2) = 0.0_WP
    Tauw_D_io(1:2) = 0.0_WP
    Ret_io(1:2) = 0.0_WP
    DUDY = 0.0_WP
    IF((iCase /= IPIPEC .AND. MYID == 0) .OR. MYID == NPSLV ) THEN ! FOR TOP OR BOTTOM WALL RANKID


        IF(iCase /= IPIPEC .AND. MYID == 0) THEN
            JID = 1
            JJ  = JCL2G(JID) !1
            JJA = 1
            N1  = 1
        END IF
        IF(MYID == NPSLV) THEN
            JID = N2DO(MYID) !LOCAL
            JJ  = JCL2G(JID) !GLOCAL
            JJA = NND2       !GLOBAL
            N1  = 2
        END IF

        !================== DU / DY ======================
        IF(flg_xzT == flgxz) THEN
            DUDY = ( U1xzL_io(JID, 1) - 0.0_WP) / (YCC(JJ) - YND(JJA))
        ELSE IF (flg_xzT == flgxzt)  THEN
            DUDY = ( U1xztL_io(JID, 1) - 0.0_WP) / (YCC(JJ) - YND(JJA))
        ELSE
        END IF

        !================= Skin variables ==================
        IF(iThermoDynamics == 0) THEN

            Tauw_io(N1) = DUDY / REN               ! undim tau_w = tau_w/ (\rhO * U2)
            Utaw_io(N1) = DSQRT(DABS(Tauw_io(N1)))
            Ret_io(N1) = REN * Utaw_io(N1)

        ELSE IF (iThermoDynamics == 1) THEN

            Tauw_io(N1) = DUDY * Mwal(N1) / REN
            Utaw_io(N1) = DSQRT( DABS(Tauw_io(N1) / Dwal(N1)) )
            Ret_io(N1) = REN * Utaw_io(N1) * Dwal(N1) / Mwal(N1)

            Tauw_D_io(N1) = Tauw_io(N1) * U0 * U0 * D0 !dimensional
            Utaw_D_io(N1) = Utaw_io(N1) * U0           !dimensional

        ELSE
        END IF

    END IF

    DUMMY(1, 1:2) = Utaw_io  (1:2)
    DUMMY(2, 1:2) = Tauw_io  (1:2)
    DUMMY(3, 1:2) = Ret_io   (1:2)
    DUMMY(4, 1:2) = Utaw_D_io(1:2)
    DUMMY(5, 1:2) = Tauw_D_io(1:2)
    DUMMY(6, 1) = YL
    DUMMY(6, 2) = DEN
    DUMMY(7, 1) = VIS
    DUMMY(7, 2) = 0.0_WP
    NSZ = 7 * 2

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(DUMMY(1, 1), DUMMY_WORK(1, 1), NSZ, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    Utaw_io(1:2) = DUMMY_WORK(1, 1:2)
    Tauw_io(1:2) = DUMMY_WORK(2, 1:2)
    Ret_io   (1:2) = DUMMY_WORK(3, 1:2)
    Utaw_D_io(1:2) = DUMMY_WORK(4, 1:2)
    Tauw_D_io(1:2) = DUMMY_WORK(5, 1:2)

    YL = DUMMY_WORK(6, 1)
    DEN = DUMMY_WORK(6, 2)
    VIS = DUMMY_WORK(7, 1)

    Tauw_io(1:2) = DABS(Tauw_io(1:2))

    IF(iThermoDynamics == 1) THEN
        Tau_diff = Tauw_io(2) - Tauw_io(1)!;  WRITE(*, *) 'Tau_diff', Tau_diff, Tauw_io(2), Tauw_io(1)
        Tau_avag = Tauw_io(2) + Tauw_io(1)!;  WRITE(*, *) 'Tau_avag', Tau_avag
        Ldist_io(1) = DABS(-1.0_WP + DABS(Tau_diff / Tau_avag))
        Ldist_io(2) = DABS( 1.0_WP + DABS(Tau_diff / Tau_avag))
        !WRITE(*, *) 'test',Ret_io(1:2), Ldist_io(1:2)
        Ret_io(1:2) = Ret_io(1:2) * Ldist_io(1:2)
    ELSE
        Ldist_io = 1.0_WP
    END IF

    !============== Averaged tauw based variables ====================
    Tauw_ave_io = 0.5_WP * (DABS(Tauw_io(1)) + DABS(Tauw_io(2)))
    DenAvew = DEN /YL
    VisAvew = VIS/YL

    !WRITE(*, *) 'DEN, YL', DEN, YL
    Utaw_ave_io = DSQRT( Tauw_ave_iO / DenAvew)
    Ret_ave_io = REN * Utaw_ave_io * DenAvew / VisAvew

    Tauw_D_ave_io = Tauw_ave_io * U0 * U0 * D0
    Utaw_D_ave_io = Utaw_ave_io * U0

    !WRITE(*, *) '# Ret_io', Ret_io(1), Ret_io(2) !test

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_Wall_thermal_properties(flg_xzt)
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: flg_xzt
    INTEGER(4) :: J, JP, JM, NSZ, JJ, JJP, JJM
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: DUMMY(8, 2)
    REAL(WP) :: DUMMY_WORK(8, 2)
    REAL(WP) :: B_tmp

    IF(iThermoDynamics /= 1) RETURN

    Hwal_RA(1:2) = 0.0_WP
    Hwal_FA(1:2) = 0.0_WP
    DHwal_RA(1:2) = 0.0_WP
    Twal(1:2) = 0.0_WP
    Dwal(1:2) = 0.0_WP
    Mwal(1:2) = 0.0_WP
    Kwal(1:2) = 0.0_WP
    Cpwal(1:2) = 0.0_WP
    Qw(1:2) = 0.0_WP

    !=============wall PARAMETERS ==undiM ==================================
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. MYID == NPSLV) THEN
        HWAL_RA(iTopWall) = H_WAL_GV(NCL1_io / 2, iTopWall)
        Dwal(iTopWall) = D_WAL_GV(NCL1_io / 2, iTopWall)
        Twal(iTopWall) = T_WAL_GV(NCL1_io / 2, iTopWall)
        Mwal(iTopWall) = M_WAL_GV(NCL1_io / 2, iTopWall)
        Kwal(iTopWall) = K_WAL_GV(NCL1_io / 2, iTopWall)
        Cpwal(iTopWall) = Cp_WAL_GV(NCL1_io / 2, iTopWall)

        HWAL_FA(iTopWall) = HWAL_RA(iTopWall)

        J  = N2DO(MYID)
        JJ = JCL2G(J)
        IF(flg_xzT == flgxz) THEN
            Qw(iTopWall) = -1.0_WP * Kwal(iTopWall) * (Twal(iTopWall) - T1xzL_io(J)) / ( YND(NND2) - YCC(NCL2) ) * CTHECD
        ELSE IF (flg_xzT == flgxzt) THEN
            Qw(iTopWall) = -1.0_WP * Kwal(iTopWall) * (Twal(iTopWall) - T1xztL_io(J)) / ( YND(NND2) - YCC(NCL2) ) * CTHECD
        ELSE
        END IF

    END IF

    IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature .AND. MYID == 0) THEN
        HWAL_RA(iBotWall) = H_WAL_GV(NCL1_io / 2, iBotWall)
        Dwal(iBotWall) = D_WAL_GV(NCL1_io / 2, iBotWall)
        Twal(iBotWall) = T_WAL_GV(NCL1_io / 2, iBotWall)
        Mwal(iBotWall) = M_WAL_GV(NCL1_io / 2, iBotWall)
        Kwal(iBotWall) = K_WAL_GV(NCL1_io / 2, iBotWall)
        Cpwal(iBotWall) = Cp_WAL_GV(NCL1_io / 2, iBotWall)

        HWAL_FA(iBotWall) = HWAL_RA(iBotWall)

        IF(flg_xzT == flgxz) THEN
            Qw(iBotWall) = -1.0_WP * Kwal(iBotWall) * (T1xzL_io(1) - Twal(iBotWall)) / (YCC(1) - YND(1)) * CTHECD
        ELSE IF (flg_xzT == flgxzt) THEN
            Qw(iBotWall) = -1.0_WP * Kwal(iBotWall) * (T1xztL_io(1) - Twal(iBotWall)) / (YCC(1) - YND(1)) * CTHECD
        ELSE
        END IF

    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. MYID == NPSLV) THEN

        J = N2DO(MYID)
        JP = J
        JM = N2DO(MYID) - 1

        JJ = JCL2G(J)
        JJP = JJ
        JJM = JJ - 1
        IF(flg_xzT == flgxz) THEN
            DHWAL_RA(iTopWall) = 2.0_WP * DHxzL_io(J) - YCL2ND_WFF(JJP) * DHxzL_io(JP) - YCL2ND_WFB(JJP) * DHxzL_io(JM)
            HWAL_RA(iTopWall) = 2.0_WP * H1xzL_io(J) - YCL2ND_WFF(JJP) * H1xzL_io(JP) - YCL2ND_WFB(JJP) * H1xzL_io(JM)
            HWAL_FA(iTopWall) = (2.0_WP * DHxzL_io(J) - YCL2ND_WFF(JJP) * DHxzL_io(JP) - YCL2ND_WFB(JJP) * DHxzL_io(JM)) /&
            (2.0_WP * D1xzL_io(J) - YCL2ND_WFF(JJP) * D1xzL_io(JP) - YCL2ND_WFB(JJP) * D1xzL_io(JM))
        ELSE IF (flg_xzT == flgxzt) THEN
            DHWAL_RA(iTopWall) = 2.0_WP * DHxztL_io(J) - YCL2ND_WFF(JJP) * DHxztL_io(JP) - YCL2ND_WFB(JJP) * DHxztL_io(JM)
            HWAL_RA(iTopWall) = 2.0_WP * H1xztL_io(J) - YCL2ND_WFF(JJP) * H1xztL_io(JP) - YCL2ND_WFB(JJP) * H1xztL_io(JM)
            HWAL_FA(iTopWall) = (2.0_WP * DHxztL_io(J) - YCL2ND_WFF(JJP) * DHxztL_io(JP) - YCL2ND_WFB(JJP) * DHxztL_io(JM)) /&
            (2.0_WP * D1xztL_io(J) - YCL2ND_WFF(JJP) * D1xztL_io(JP) - YCL2ND_WFB(JJP) * D1xztL_io(JM))
        ELSE
        END IF
        ! check Hwal_RA two methods...
        CALL THERM_PROP_UPDATE_FROM_DH(DHWAL_RA(iTopWall), Hwal_RA(iTopWall), Twal(iTopWall), Dwal(iTopWall), &
        Mwal(iTopWall), Kwal(iTopWall), Cpwal(iTopWall), B_tmp)

        IF(flg_xzT == flgxz) THEN
            Qw(iTopWall) = -1.0_WP * Kwal(iTopWall) * (Twal(iTopWall) - T1xzL_io(J)) / ( YND(NND2) - YCC(NCL2) ) * CTHECD
        ELSE IF (flg_xzT == flgxzt) THEN
            Qw(iTopWall) = -1.0_WP * Kwal(iTopWall) * (Twal(iTopWall) - T1xztL_io(J)) / ( YND(NND2) - YCC(NCL2) ) * CTHECD
        ELSE
        END IF

    END IF


    IF((iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux .OR. iCase == iPIPEC ) .AND. MYID == 0) THEN

        J = 1
        JP = 2
        JM = 1

        IF(flg_xzT == flgxz) THEN
            DHWAL_RA(iBotWall) = 2.0_WP * DHxzL_io(J) - YCL2ND_WFF(JP) * DHxzL_io(JP) - YCL2ND_WFB(JP) * DHxzL_io(JM)
            HWAL_RA(iBotWall) = 2.0_WP * H1xzL_io(J) - YCL2ND_WFF(JP) * H1xzL_io(JP) - YCL2ND_WFB(JP) * H1xzL_io(JM)
            HWAL_FA(iBotWall) = (2.0_WP * DHxzL_io(J) - YCL2ND_WFF(JP) * DHxzL_io(JP) - YCL2ND_WFB(JP) * DHxzL_io(JM)) / &
            2.0_WP * D1xzL_io(J) - YCL2ND_WFF(JP) * D1xzL_io(JP) - YCL2ND_WFB(JP) * D1xzL_io(JM)
        ELSE IF (flg_xzT == flgxzt) THEN
            DHWAL_RA(iBotWall) = 2.0_WP * DHxztL_io(J) - YCL2ND_WFF(JP) * DHxztL_io(JP) - YCL2ND_WFB(JP) * DHxztL_io(JM)
            HWAL_RA(iBotWall) = 2.0_WP * H1xztL_io(J) - YCL2ND_WFF(JP) * H1xztL_io(JP) - YCL2ND_WFB(JP) * H1xztL_io(JM)
            HWAL_FA(iBotWall) = (2.0_WP * DHxztL_io(J) - YCL2ND_WFF(JP) * DHxztL_io(JP) - YCL2ND_WFB(JP) * DHxztL_io(JM)) / &
            2.0_WP * D1xztL_io(J) - YCL2ND_WFF(JP) * D1xztL_io(JP) - YCL2ND_WFB(JP) * D1xztL_io(JM)
        ELSE
        END IF

        ! check Hwal_RA two methods...
        CALL THERM_PROP_UPDATE_FROM_DH(DHWAL_RA(iBotWall), Hwal_RA(iBotWall), Twal(iBotWall), Dwal(iBotWall), &
        Mwal(iBotWall), Kwal(iBotWall), Cpwal(iBotWall), B_tmp)

        IF(flg_xzT == flgxz) THEN
            Qw(iBotWall) = -1.0_WP * Kwal(iBotWall) * (T1xzL_io(1) - Twal(iBotWall)) / (YCC(1) - YND(1)) * CTHECD
        ELSE IF (flg_xzT == flgxzt) THEN
            Qw(iBotWall) = -1.0_WP * Kwal(iBotWall) * (T1xztL_io(1) - Twal(iBotWall)) / (YCC(1) - YND(1)) * CTHECD
        ELSE
        END IF
    END IF


    DUMMY(1, 1:2) = Hwal_RA(1:2)
    DUMMY(2, 1:2) = Hwal_FA(1:2)
    DUMMY(3, 1:2) = Twal(1:2)
    DUMMY(4, 1:2) = Dwal(1:2)
    DUMMY(5, 1:2) = Mwal(1:2)
    DUMMY(6, 1:2) = Kwal(1:2)
    DUMMY(7, 1:2) = Cpwal(1:2)
    DUMMY(8, 1:2) = Qw(1:2)
    NSZ = 8 * 2

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(DUMMY, DUMMY_WORK, NSZ, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    Hwal_RA(1:2) = DUMMY_WORK(1, 1:2)
    Hwal_FA(1:2) = DUMMY_WORK(2, 1:2)
    Twal(1:2) = DUMMY_WORK(3, 1:2)
    Dwal(1:2) = DUMMY_WORK(4, 1:2)
    Mwal(1:2) = DUMMY_WORK(5, 1:2)
    Kwal(1:2) = DUMMY_WORK(6, 1:2)
    Cpwal(1:2) = DUMMY_WORK(7, 1:2)
    Qw(1:2)    = DUMMY_WORK(8, 1:2)

    RETURN
END SUBROUTINE
