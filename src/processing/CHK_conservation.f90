!**********************************************************************************************************************************
!> @brief
!>        to check mass conservation
!> @details
!> SUBROUTINE: CHK_MassConsv_io(in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CHK_MassConsv_io
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    REAL(WP) :: FluxIntg_INL,  FluxIntg_INL_WORK
    REAL(WP) :: FluxIntg_OUT,  FluxIntg_OUT_WORK
    REAL(WP) :: MASS_RATE_INTG, MASS_RATE_INTG_WORK
    REAL(WP) :: MASS_RATE, COE2
    REAL(WP) :: DUMMY(3), DUMMY_WORK(3)

    INTEGER(4) :: IC, JC, KC, JJ

    !======== inlet /OUTLET ===========================
    FluxIntg_INL = 0.0_WP
    FluxIntg_OUT = 0.0_WP
    IF(TgFlowFlg) THEN
        DO JC = 1, N2DO(MYID)
            JJ = JCL2G(JC)
            COE2 = 1.0_WP / DYFI(JJ) / RCCI1(JJ)
            DO KC = 1, NCL3
                FluxIntg_INL = FluxIntg_INL + G_io(1,    JC, KC, NFLOW) * COE2
                FluxIntg_OUT = FluxIntg_OUT + G_io(NCL1E, JC, KC, NFLOW) * COE2
            END DO
        END DO
        FluxIntg_INL = FluxIntg_INL / DZI
        FluxIntg_OUT = FluxIntg_OUT / DZI
    END IF

    !======== INSIDE ===========================
    MASS_RATE_INTG = 0.0_WP
    DO JC = 1, N2DO(MYID)
        JJ = JCL2G(JC)
        COE2 = 1.0_WP / DYFI(JJ) / RCCI1(JJ)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_io
                !MASS_RATE    = ( DENSITY(IC, JC, KC) - DENSITYP(IC, JC, KC) ) * COE2
                MASS_RATE    = DrhoDtP(IC, JC, KC) * COE2
                MASS_RATE_INTG = MASS_RATE + MASS_RATE_INTG
            END DO
        END DO
    END DO
    MASS_RATE_INTG = MASS_RATE_INTG / DZI / DXI

    !=========== ADD TOGETHER FROM ALL RANKS ===============================
    DUMMY(1) = FluxIntg_INL
    DUMMY(2) = FluxIntg_OUT
    DUMMY(3) = MASS_RATE_INTG
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(DUMMY, DUMMY_WORK, 3, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    FluxIntg_INL_WORK = DUMMY_WORK(1)
    FluxIntg_OUT_WORK = DUMMY_WORK(2)
    MASS_RATE_INTG_WORK = DUMMY_WORK(3)

    !! method 1
    IF(TgFlowFlg) &
    CHK_Mass_CONSV0_SCALING = (FluxIntg_INL_WORK - MASS_RATE_INTG_WORK) / FluxIntg_OUT_WORK
    !CHK_MASS_CONSV0       = DABS(CHK_Mass_CONSV0_SCALING - 1.0_WP)

    !! method 2
    CHK_MASS_CONSV0 = FluxIntg_INL_WORK - MASS_RATE_INTG_WORK - FluxIntg_OUT_WORK

    !IF(MYID == 0) WRITE(*, '(A, I3.1, 5ES13.5)') 'mass coe = ', MYID, CHK_MASS_CONSV0,CHK_Mass_CONSV0_SCALING, &
    !FluxIntg_INL_WORK, MASS_RATE_INTG_WORK, FluxIntg_OUT_WORK


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CHK_EnegConsv_io
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE

    REAL(WP) :: FluxIntg_INL,  FluxIntg_INL_WORK
    REAL(WP) :: FluxIntg_OUT,  FluxIntg_OUT_WORK
    REAL(WP) :: ENEG_RATE_INTG, ENEG_RATE_INTG_WORK, ENEG_TOTAL, ENEG_TOTAL_WORK

    REAL(WP) :: DHV_BTW, KDT_BTW
    REAL(WP) :: DHV_TPW, KDT_TPW
    REAL(WP) :: FluxIntg_BTW, FluxIntg_BTW_WORK
    REAL(WP) :: FluxIntg_TPW, FluxIntg_TPW_WORK
    REAL(WP) :: DHU_INL, DHU_OUT, KDT_INL, KDT_OUT!, WHEATFLUX
    INTEGER(4) :: IC, JC, KC, JJ, JJP, JM, JP, IC1, IC2, IM1, IP2
    REAL(WP) :: COE1, COE2, COE221, COE222
    REAL(WP) :: DUMMY(6), DUMMY_WORK(6)



    !CALL CHK_EnegConsv_each_Cell_io ! CHECK

    IF(iThermoDynamics /= 1) RETURN

    !============== inlet /OUTLET AT Y-Z PLANE =======================
    FluxIntg_INL = 0.0_WP
    FluxIntg_OUT = 0.0_WP
    IF(TgFlowFlg) THEN
        IC1 = 1
        IM1 = IMV_io(IC1)

        IC2 = NCL1_io
        IP2 = IPV_io(IC2)

        DO JC = 1, N2DO(MYID)
            JJ = JCL2G(JC)
            COE1 = XND2CL / DYFI(JJ) / RCCI1(JJ)
            COE2 = 1.0_WP / DYFI(JJ) / RCCI1(JJ)
            DO KC = 1, NCL3
                DHU_INL = ( ENTHALPY(IM1, JC, KC) + ENTHALPY(IC1, JC, KC) ) * G_io(IC1, JC, KC, 1)
                KDT_INL = ( THERMCONDT(IC1, JC, KC) + THERMCONDT(IM1, JC, KC) ) * &
                ( TEMPERATURE(IC1, JC, KC) - TEMPERATURE(IM1, JC, KC) ) * DXI * CTHECD

                DHU_OUT = ( ENTHALPY(IC2, JC, KC) + ENTHALPY(IP2, JC, KC) ) * G_io(IP2, JC, KC, 1)
                KDT_OUT = ( THERMCONDT(IP2, JC, KC) + THERMCONDT(IC2, JC, KC) ) * &
                ( TEMPERATURE(IP2, JC, KC) - TEMPERATURE(IC2, JC, KC) ) * DXI * CTHECD

                FluxIntg_INL  = FluxIntg_INL + (DHU_INL- KDT_INL) * COE1
                FluxIntg_OUT  = FluxIntg_OUT + (DHU_OUT - KDT_OUT) * COE1

            END DO
        END DO
        FluxIntg_INL = FluxIntg_INL / DZI
        FluxIntg_OUT = -FluxIntg_OUT / DZI
    END IF

    !================= INSIDE =================================
    ENEG_RATE_INTG = 0.0_WP
    ENEG_TOTAL = 0.0_WP
    DO JC = 1, N2DO(MYID)
        JJ = JCL2G(JC)
        COE2 = 1.0_WP / DYFI(JJ) / RCCI1(JJ)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_io
                ENEG_RATE_INTG = (DH(IC, JC, KC) - DH0(IC, JC, KC)) * COE2 + ENEG_RATE_INTG
                ENEG_TOTAL = ENEG_TOTAL + DH(IC, JC, KC) * COE2
            END DO
        END DO
    END DO
    ENEG_RATE_INTG = ENEG_RATE_INTG / DZI / DXI / DT
    ENEG_TOTAL = ENEG_TOTAL / DZI / DXI

    !==============BOTTOM WALL AT X -Z PLANE =======================
    FluxIntg_BTW = 0.0_WP
    IF(MYID == 0) THEN
        JC = 1
        JP = JLPV(JC)
        JM = JLMV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)

        IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux    ) THEN
            IF(iCase == iPIPEC) THEN
                COE222 = 0.0_WP
            ELSE
                COE222 = 1.0_WP / RNDI1(JJ) * CTHECD
            END IF
        END IF

        IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature    ) THEN
            IF(iCase == iPIPEC) THEN
                COE222 = 0.0_WP
            ELSE
                COE222 = DYCI(JJ) / RNDI1(JJ) * CTHECD
            END IF
        END IF

        DO IC = 1, NCL1_io
            DO KC = 1, NCL3
                DHV_BTW = 0.0_WP
                !KDT_BTW = WHEATFLUX * CTHECD

                IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux    ) KDT_BTW = -WALLFLUX(IC, 1, KC) * COE222
                IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature ) THEN
                    KDT_BTW = ( YCL2ND_WFF(JJ) * THERMCONDT(IC, JC, KC) +   &
                    YCL2ND_WFB(JJ) * THERMCONDT(IC, JM, KC) ) * &
                    ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IC, JM, KC) ) * COE222        !DYCI(JJ) / RNDI(JJ)
                END IF

                FluxIntg_BTW  = FluxIntg_BTW + (DHV_BTW - KDT_BTW)
            END DO
        END DO
        FluxIntg_BTW = FluxIntg_BTW/ DXI / DZI
        IF(iCase == iPIPEC) FluxIntg_BTW = 0.0_WP
    END IF

    !============== TOP WALL AT X -Z PLANE =======================
    FluxIntg_TPW = 0.0_WP
    IF(MYID == NPSLV) THEN
        JC = N2DO(MYID)
        JP = JLPV(JC)
        JM = JLMV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)

        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux    ) THEN
            COE221 = 1.0_WP / RNDI1(JJP) * CTHECD
        END IF
        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature    ) THEN
            COE221 = DYCI(JJP) / RNDI1(JJP) * CTHECD
        END IF

        DO IC = 1, NCL1_io

            DO KC = 1, NCL3
                DHV_TPW = 0.0_WP
                IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux    ) KDT_TPW = -WALLFLUX(IC, 2, KC) * COE221
                IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature ) THEN
                    KDT_TPW = ( YCL2ND_WFF(JJP) * THERMCONDT(IC, JP, KC) +   &
                    YCL2ND_WFB(JJP) * THERMCONDT(IC, JC, KC) ) * &
                    ( TEMPERATURE(IC, JP, KC) - TEMPERATURE(IC, JC, KC) ) * COE221   !DYCI(JJP) / RNDI(JJP)
                END IF

                FluxIntg_TPW  = FluxIntg_TPW + (DHV_TPW + KDT_TPW) !positive
            END DO
        END DO
        FluxIntg_TPW = FluxIntg_TPW / DXI / DZI / RNDI1(NND2)
    END IF

    DUMMY(1) = FluxIntg_INL
    DUMMY(2) = FluxIntg_OUT
    DUMMY(3) = ENEG_RATE_INTG
    DUMMY(4) = FluxIntg_BTW
    DUMMY(5) = FluxIntg_TPW
    DUMMY(6) = ENEG_TOTAL
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(DUMMY,  DUMMY_WORK,  6, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    FluxIntg_INL_WORK  = DUMMY_WORK(1)
    FluxIntg_OUT_WORK  = DUMMY_WORK(2)
    ENEG_RATE_INTG_WORK  = DUMMY_WORK(3)
    FluxIntg_BTW_WORK  = DUMMY_WORK(4)
    FluxIntg_TPW_WORK  = DUMMY_WORK(5)
    ENEG_TOTAL_WORK    = DUMMY_WORK(6)

    !CALL MPI_ALLREDUCE(FluxIntg_INL,  FluxIntg_INL_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(FluxIntg_OUT,  FluxIntg_OUT_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(ENEG_RATE_INTG, ENEG_RATE_INTG_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(FluxIntg_BTW,  FluxIntg_BTW_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(FluxIntg_TPW,  FluxIntg_TPW_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    !CALL MPI_ALLREDUCE(ENEG_TOTAL,ENEG_TOTAL_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    CHK_ENEG_CONSV0 = ENEG_RATE_INTG_WORK + FluxIntg_INL_WORK + FluxIntg_OUT_WORK + FluxIntg_BTW_WORK + FluxIntg_TPW_WORK
    CHK_ENEG_TOTAL  = ENEG_TOTAL_WORK
    !IF(MYID == 0) WRITE(*, *) 'eng coe = ',CHK_ENEG_CONSV0,ENEG_RATE_INTG_WORK, &
    !-FluxIntg_BTW_WORK, FluxIntg_TPW_WORK, FluxIntg_TPW_WORK -FluxIntg_BTW_WORK

    RETURN
END SUBROUTINE


!**********************************************************************
SUBROUTINE CHK_EnegConsv_each_Cell_io
    USE init_info
    USE mesh_info
    USE flow_info
    USE thermal_info
    IMPLICIT NONE



    INTEGER(4) :: IC, IP, IM
    INTEGER(4) :: KC, KP, KM
    INTEGER(4) :: JC, JP, JM, JJ, JJP

    REAL(WP) :: TOTLENG, ENEG_RATE_INTG, G1HP, G1HN, G2HP, G2HN, G3HP, G3HN, KXTP, KXTN, KZTP, KZTN, KYTP, KYTN
    REAL(WP) :: G1H, G2H, G3H, KXT, KYT, KZT

    INTEGER(4) :: FLID
    CHARACTER(4) :: NNN


    IF(iThermoDynamics /= 1) RETURN

    WRITE(NNN, '(1I4.4)') MYID
    FLID = MYID + 100
    OPEN(FLID, FILE = 'OUT_VARS' // NNN // '.debug',POSITION = "APPEND")
    IF(MYID == 0) WRITE(FLID, *) '#===check energy equatioN === '

    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JM = JLMV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)
        DO IC = 1, NCL1_io
            IP = IPV_io(IC)
            IM = IMV_io(IC)
            DO KC = 1, NCL3
                KP = KPV(KC)
                KM = KMV(KC)

                !=== DH =============================
                ENEG_RATE_INTG = (DH(IC, JC, KC) - DH0(IC, JC, KC)) / DYFI(JJ) / DZI / DXI / DT

                !===G1H ===============
                G1HP = ( ENTHALPY(IP, JC, KC) + ENTHALPY(IC, JC, KC) ) * G_io(IP, JC, KC, 1) / DYFI(JJ) / DZI
                G1HN = ( ENTHALPY(IC, JC, KC) + ENTHALPY(IM, JC, KC) ) * G_io(IC, JC, KC, 1) / DYFI(JJ) / DZI

                !===G2H ===============
                G2HP = ( YCL2ND_WFF(JJP) * ENTHALPY(IC, JP, KC) + &
                YCL2ND_WFB(JJP) * ENTHALPY(IC, JC, KC) ) * G_io(IC, JP, KC, 2) / DZI / DXI
                G2HN = ( YCL2ND_WFF(JJ) * ENTHALPY(IC, JC, KC) +                         &
                YCL2ND_WFB(JJ) * ENTHALPY(IC, JM, KC) ) * G_io(IC, JC, KC, 2) / DZI / DXI

                !===G3H ===============
                G3HP = ( ENTHALPY(IC, JC, KP) + ENTHALPY(IC, JC, KC) ) * G_io(IC, JC, KP, 3) / DYFI(JJ) / DXI
                G3HN = ( ENTHALPY(IC, JC, KC) + ENTHALPY(IC, JC, KM) ) * G_io(IC, JC, KC, 3) / DYFI(JJ) / DXI

                !======== KDTX =========================
                KXTP = ( THERMCONDT(IP, JC, KC) + THERMCONDT(IC, JC, KC) ) * 0.5_WP * &
                ( TEMPERATURE(IP, JC, KC) - TEMPERATURE(IC, JC, KC) ) * DXI / DYFI(JJ) / DZI * CTHECD
                KXTN = ( THERMCONDT(IC, JC, KC) + THERMCONDT(IM, JC, KC) ) * 0.5_WP * &
                ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IM, JC, KC) ) * DXI / DYFI(JJ) / DZI * CTHECD

                !======== KDTZ =========================
                KZTP = ( THERMCONDT(IC, JC, KP) + THERMCONDT(IC, JC, KC) ) * 0.5_WP * &
                ( TEMPERATURE(IC, JC, KP) - TEMPERATURE(IC, JC, KC) ) * DZI / DYFI(JJ) / DXI * CTHECD
                KZTN = ( THERMCONDT(IC, JC, KM) + THERMCONDT(IC, JC, KC) ) * 0.5_WP * &
                ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IC, JC, KM) ) * DZI / DYFI(JJ) / DXI * CTHECD

                !======== KDTY =========================
                KYTP = ( YCL2ND_WFF(JJP) * THERMCONDT(IC, JP, KC) +   &
                YCL2ND_WFB(JJP) * THERMCONDT(IC, JC, KC) ) * &
                ( TEMPERATURE(IC, JP, KC) - TEMPERATURE(IC, JC, KC) ) * DYCI(JJP) / DZI / DXI * CTHECD
                KYTN = ( YCL2ND_WFF(JJ) * THERMCONDT(IC, JC, KC) +   &
                YCL2ND_WFB(JJ) * THERMCONDT(IC, JM, KC) ) * &
                ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IC, JM, KC) ) * DYCI(JJ) / DZI / DXI * CTHECD
                IF(JJ == 1 .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
                    KYTN = -WALLFLUX(IC, iBotWall, KC) / DZI / DXI * CTHECD
                ELSE IF(JJ == NCL2 .AND. iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
                    KYTP = -WALLFLUX(IC, iTopWall, KC) / DZI / DXI * CTHECD
                ELSE
                END IF

                G1H = G1Hn - G1HP
                G2H = G2Hn - G2HP
                G3H = G3Hn - G3HP
                KXT = KXTP - KXTN
                KYT = KYTP - KYTN
                KZT = KZTP - KZTN

                TOTLENG = - ENEG_RATE_INTG + G1H + G2H + G3H + KXT + KZT + KYT
                IF( DABS(TOTLENG) >  1.0E-5_WP) THEN
                    WRITE(FLID, *) JJ, IC, KC, &
                    TOTLENG, - ENEG_RATE_INTG,G1H + G2H + G3H, KXT + KyT + KzT
                    WRITE(FLID, *) G1HP, G1HN, G1H
                    WRITE(FLID, *) G2HP, G2HN, G2H
                    WRITE(FLID, *) G3HP, G3HN, G3H
                    WRITE(FLID, *) KXTP, KXTN, KXT
                    WRITE(FLID, *) KYTP, KYTN, KYT
                    WRITE(FLID, *) KZTP, KZTN, KZT
                END IF
            END DO
        END DO
    END DO


    RETURN
END SUBROUTINE
