!**********************************************************************************************************************************
!> @brief
!>        to calculate the viscous term
!> @details
!> SUBROUTINE: VISCOUS_ALL_EXPLT_X_io(in MYID = all)
!> SUBROUTINE: VISCOUS_PAR_EXPLT_X_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
! 09/2020 - Added more fluid types and optimized, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE VISCOUS_ALL_EXPLT_X_io
    USE MESH_INFO
    USE FLOW_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP, NXI
    INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
    INTEGER(4) :: KC, KM, KP

    REAL(WP) :: VIS12C, VIS12P
    REAL(WP) :: VIS13C, VIS13P
    REAL(WP) :: DUDX, DUDY, DVDX, DUDZ, DWDX
    REAL(WP) :: TAU11F, TAU11B, DTAU11DX
    REAL(WP) :: TAU12F, TAU12B, DTAU12DY
    REAL(WP) :: TAU13F, TAU13B, DTAU13DZ
    REAL(WP) :: COE1, COE2, COE3,  COE32 !, TESTVAL, TESTVAL1, TESTVAL2, TESTVAL3

    RHS_io = 0.0_WP

    COE1 = 2.0_WP * DXI * CVISC
    !COE3 = DZI * XND2CL * ZND2CL * CVISC
    COE3 = DZI * CVISC

    NXI = 1
    IF(TgFlowFlg) NXI  = 2

    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJM = JGMV(JJ)
        JJP = JGPV(JJ)
        !COE2 = DYFI(JJ) * XND2CL * CVISC * RCCI1(JJ)
        COE2 = DYFI(JJ) * CVISC * RCCI1(JJ)
        COE32 = COE3 * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)

            DO IC = NXI, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)

                !================ DX_TAU_11 =====================================
                ! at (i, J, K)
                DUDX = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
                TAU11F = 0.5_WP * ( Viscousity0(IC, JC, KC) + Viscousity0(IC, JC, KC) ) * ( DUDX - DivU_io(IC, JC, KC) )

                ! at (I - 1, J, K)
                DUDX = ( Q_io(IC, JC, KC, 1) - Q_io(IM, JC, KC, 1) ) * DXI
                TAU11B = 0.5_WP * ( Viscousity0(IM, JC, KC) + Viscousity0(IM, JC, KC) ) * ( DUDX - DivU_io(IM, JC, KC) )

                ! at (i', J, K)
                DTAU11DX = (TAU11F - TAU11B) * COE1

                !================ DY_TAU_12 =====================================
                ! at (i', J'+ 1, K)
                !VIS12P = ( Viscousity(IC, JP, KC) + Viscousity(IM, JP, KC) ) * YCL2ND_WFF(JJP) + &
                !         ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) * YCL2ND_WFB(JJP)
                VIS12P = MU_STG(IC, JP, KC, 1)
                DUDY = ( Q_io(IC, JP, KC, 1) - Q_io(IC, JC, KC, 1) ) * DYCI(JJP) / RNDI1(JJP)
                DVDX = ( Q_io(IC, JP, KC, 2) - Q_io(IM, JP, KC, 2) ) * DXI
                TAU12F = VIS12P * ( DUDY + DVDX )

                ! at (i', J', K)
                !VIS12C = ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) * YCL2ND_WFF(JJ) + &
                !         ( Viscousity(IC, JM, KC) + Viscousity(IM, JM, KC) ) * YCL2ND_WFB(JJ)
                VIS12C = MU_STG(IC, JC, KC, 1)
                IF(iCase == IPIPEC .AND. JJ == 1) THEN
                    DUDY = 0.0_WP
                ELSE
                    DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ) / RNDI1(JJ)
                END IF
                DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI
                TAU12B = VIS12C * ( DUDY + DVDX )

                ! at (i', J, K)
                DTAU12DY = (TAU12F - TAU12B) * COE2
                !IF(MYID == 0) WRITE(*, *) 'vIScx', JC, KC, IC, &
                !TAU12F, TAU12B, VIS12P, VIS12C

                !================ DZ_TAU_13 =====================================
                ! at (i', J, K'+ 1)
                !VIS13P = ( Viscousity(IC, JC, KP) + Viscousity(IM, JC, KP) ) + &
                !         ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) )
                VIS13P = MU_STG(IC, JC, KP, 2)
                DUDZ = ( Q_io(IC, JC, KP, 1) - Q_io(IC, JC, KC, 1) ) * DZI
                DWDX = ( Q_io(IC, JC, KP, 3) - Q_io(IM, JC, KP, 3) ) * DXI
                TAU13F = VIS13P * ( DUDZ + DWDX )

                ! at (i', J, K')
                !VIS13C = ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) + &
                !         ( Viscousity(IC, JC, KM) + Viscousity(IM, JC, KM) )
                VIS13C = MU_STG(IC, JC, KC, 2)
                DUDZ = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JC, KM, 1) ) * DZI
                DWDX = ( Q_io(IC, JC, KC, 3) - Q_io(IM, JC, KC, 3) ) * DXI
                TAU13B = VIS13C * ( DUDZ + DWDX )

                ! at (i', J, K)
                DTAU13DZ = (TAU13F - TAU13B) * COE32

                !================ D_TAU_X direction =================================
                QTMP_io(IC, JC, KC) = QTMP_io(IC, JC, KC) + DTAU11DX + DTAU12DY + DTAU13DZ

                !IF(MYID == 0) WRITE(*, *) 'vIScx', JC, KC, IC, &
                !DTAU11DX, DTAU12DY, DTAU13DZ, QTMP_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !CALL DEBUG_WRT_LOCAL(QTMP_io, 1, N2DO(MYID), 'vISx') ! test

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE VISCOUS_PAR_EXPLT_X_io
    ! Refer to NICoud (2000) for the split of vIScous terms.
    USE MESH_INFO
    USE FLOW_INFO
    USE THERMAL_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP, NXI
    INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
    INTEGER(4) :: KC, KM, KP

    REAL(WP) :: VIS12C, VIS12P
    REAL(WP) :: VIS13C, VIS13P
    REAL(WP) :: DUDX, DUDY, DUDZ, DVDX, DVDY, DWDX, DWDZ
    REAL(WP) :: TAU11F_EX, TAU11B_EX, DTAU11DX_EX
    REAL(WP) :: TAU12F_EX, TAU12B_EX, DTAU12DY_EX
    REAL(WP) :: TAU13F_EX, TAU13B_EX, DTAU13DZ_EX
    REAL(WP) :: TAU11F_IM, TAU11B_IM, DTAU11DX_IM
    REAL(WP) :: TAU12F_IM, TAU12B_IM, DTAU12DY_IM
    REAL(WP) :: TAU13F_IM, TAU13B_IM, DTAU13DZ_IM
    REAL(WP) :: TAU11F_AD, TAU11B_AD, DTAU11DX_AD
    REAL(WP) :: TAU12F_AD, TAU12B_AD, DTAU12DY_AD
    REAL(WP) :: TAU13F_AD, TAU13B_AD, DTAU13DZ_AD
    REAL(WP) :: COE1, COE2, COE3,  COE32 !, TESTVAL, TESTVAL1, TESTVAL2, TESTVAL3

    RHS_io = 0.0_WP

    COE1 = -2.0_WP / 3.0_WP * DXI * CVISC
    !COE3 = DZI * XND2CL * ZND2CL * CVISC
    COE3 = DZI * CVISC

    NXI = 1
    IF(TgFlowFlg) NXI  = 2

    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJM = JGMV(JJ)
        JJP = JGPV(JJ)
        !COE2 = DYFI(JJ) * XND2CL * CVISC * RCCI1(JJ)
        COE2 = DYFI(JJ) * CVISC * RCCI1(JJ)
        COE32 = COE3 * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)

            DO IC = NXI, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)

                !================ DX_TAU_11 = (pARtly) ====================================
                ! at (i, J, K)
                DUDX = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
                DVDY = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * DYFI(JJ) * RCCI1(JJ)
                DWDZ = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * DZI * RCCI2(JJ)
                TAU11F_EX = Viscousity(IC, JC, KC) * ( DVDY + DWDZ )
                TAU11F_IM = Viscousity(IC, JC, KC) * DUDX

                DUDX = ( G_io(IP, JC, KC, 1) * DRHOI_STG(IP, JC, KC, 1) - &
                G_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) ) * DXI
                TAU11F_AD = Viscousity(IC, JC, KC) * DUDX

                ! at (I - 1, J, K)
                DUDX = ( Q_io(IC, JC, KC, 1) - Q_io(IM, JC, KC, 1) ) * DXI
                DVDY = ( Q_io(IM, JP, KC, 2) - Q_io(IM, JC, KC, 2) ) * DYFI(JJ) * RCCI1(JJ)
                DWDZ = ( Q_io(IM, JC, KP, 3) - Q_io(IM, JC, KC, 3) ) * DZI * RCCI2(JJ)
                TAU11B_EX = Viscousity(IM, JC, KC) * ( DVDY + DWDZ )
                TAU11B_IM = Viscousity(IM, JC, KC) * DUDX

                DUDX = ( G_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) - &
                G_io(IM, JC, KC, 1) * DRHOI_STG(IM, JC, KC, 1) ) * DXI
                TAU11B_AD = Viscousity(IM, JC, KC) * DUDX

                ! at (i', J, K)
                DTAU11DX_EX = (TAU11F_EX - TAU11B_EX) * COE1
                DTAU11DX_IM = (TAU11F_IM - TAU11B_IM) * COE1 * (-2.0_WP)

                DTAU11DX_AD = (TAU11F_AD - TAU11B_AD) * COE1 * (-2.0_WP)

                !================ DY_TAU_12 =====================================
                ! at (i', J'+ 1, K)
                !VIS12P = ( Viscousity(IC, JP, KC) + Viscousity(IM, JP, KC) ) * YCL2ND_WFF(JJP) + &
                !         ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) * YCL2ND_WFB(JJP)
                VIS12P = MU_STG(IC, JP, KC, 1)
                DVDX = ( Q_io(IC, JP, KC, 2) - Q_io(IM, JP, KC, 2) ) * DXI
                DUDY = ( Q_io(IC, JP, KC, 1) - Q_io(IC, JC, KC, 1) ) * DYCI(JJP) / RNDI1(JJP)
                TAU12F_EX = VIS12P * ( DVDX )
                TAU12F_IM = VIS12P * ( DUDY )

                DUDY = ( G_io(IC, JP, KC, 1) * DRHOI_STG(IC, JP, KC, 1) - &
                G_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) ) * DYCI(JJP) / RNDI1(JJP)
                TAU12F_AD = VIS12P * ( DUDY )
                ! at (i', J', K)
                !VIS12C = ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) * YCL2ND_WFF(JJ) + &
                !         ( Viscousity(IC, JM, KC) + Viscousity(IM, JM, KC) ) * YCL2ND_WFB(JJ)
                VIS12C = MU_STG(IC, JC, KC, 1)
                IF(iCase == IPIPEC .AND. JJ == 1) THEN
                    DUDY = 0.0_WP
                ELSE
                    DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ) / RNDI1(JJ)
                END IF
                DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI
                TAU12B_EX = VIS12C * ( DVDX )
                TAU12B_IM = VIS12C * ( DUDY )

                DUDY = ( G_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) - &
                G_io(IC, JM, KC, 1) * DRHOI_STG(IC, JM, KC, 1) ) * DYCI(JJ) / RNDI1(JJ)
                TAU12B_AD = VIS12C * ( DUDY )

                ! at (i', J, K)
                DTAU12DY_EX = (TAU12F_EX - TAU12B_EX) * COE2
                DTAU12DY_IM = (TAU12F_IM - TAU12B_IM) * COE2

                DTAU12DY_AD = (TAU12F_AD - TAU12B_AD) * COE2

                !================ DZ_TAU_13 =====================================
                ! at (i', J, K'+ 1)
                !VIS13P = ( Viscousity(IC, JC, KP) + Viscousity(IM, JC, KP) ) + &
                !         ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) )
                VIS13P = MU_STG(IC, JC, KP, 2)
                DUDZ = ( Q_io(IC, JC, KP, 1) - Q_io(IC, JC, KC, 1) ) * DZI
                DWDX = ( Q_io(IC, JC, KP, 3) - Q_io(IM, JC, KP, 3) ) * DXI
                TAU13F_EX = VIS13P * ( DWDX )
                TAU13F_IM = VIS13P * ( DUDZ )

                DUDZ = ( Q_io(IC, JC, KP, 1) * DRHOI_STG(IC, JC, KP, 1) - &
                Q_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) ) * DZI
                TAU13F_AD = VIS13P * ( DUDZ )

                ! at (i', J, K')
                !VIS13C = ( Viscousity(IC, JC, KC) + Viscousity(IM, JC, KC) ) + &
                !         ( Viscousity(IC, JC, KM) + Viscousity(IM, JC, KM) )
                VIS13C = MU_STG(IC, JC, KC, 2)
                DUDZ = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JC, KM, 1) ) * DZI
                DWDX = ( Q_io(IC, JC, KC, 3) - Q_io(IM, JC, KC, 3) ) * DXI
                TAU13B_EX = VIS13C * ( DWDX )
                TAU13B_IM = VIS13C * ( DUDZ )

                DUDZ = ( Q_io(IC, JC, KC, 1) * DRHOI_STG(IC, JC, KC, 1) - &
                Q_io(IC, JC, KM, 1) * DRHOI_STG(IC, JC, KM, 1) ) * DZI
                TAU13B_AD = VIS13C * ( DUDZ )

                ! at (i', J, K)
                DTAU13DZ_EX = (TAU13F_EX - TAU13B_EX) * COE32
                DTAU13DZ_IM = (TAU13F_IM - TAU13B_IM) * COE32

                DTAU13DZ_AD = (TAU13F_AD - TAU13B_AD) * COE32

                !================ D_TAU_X direction =================================
                QTMP_io(IC, JC, KC) = DTAU11DX_EX + DTAU12DY_EX + DTAU13DZ_EX + QTMP_io(IC, JC, KC)
                RHS_io(IC, JC, KC) = DTAU11DX_IM + DTAU12DY_IM + DTAU13DZ_IM
                DIVU_io(IC, JC, KC) = DTAU11DX_AD + DTAU12DY_AD + DTAU13DZ_AD
                !WRITE(*, *) DTAU11DX_AD, DTAU12DY_AD, DTAU13DZ_AD
                !IF(JJ == 2 .AND. IC == 1 .AND. KC == 1) WRITE(*, '(A, 4ES13.5)') 'vIScx', &
                !DTAU11DX, DTAU12DY, DTAU13DZ,RHS_io(IC, JC, KC)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
