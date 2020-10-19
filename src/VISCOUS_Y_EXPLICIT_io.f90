!**********************************************************************************************************************************
!> @brief
!>        to calculate the viscous term
!> @details
!> SUBROUTINE: VISCOUS_ALL_EXPLT_Y_io(in MYID = all)
!> SUBROUTINE: VISCOUS_PAR_EXPLT_Y_io(in MYID = all)
!> SUBROUTINE: CHECK_GradP_ON_WALL
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
! 09/2020 - Added more fluid types and optimized, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE VISCOUS_ALL_EXPLT_Y_io
    USE init_info
    USE MESH_INFO
    USE FLOW_INFO
    USE THERMAL_INFO
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP, KS
    INTEGER(4) :: NYI

    REAL(WP) :: VIS21C, VIS21P
    REAL(WP) :: VIS23C, VIS23P
    REAL(WP) :: DVDY, DVDX, DUDY, DVDZ, DWDY
    REAL(WP) :: TAU21F, TAU21B, DTAU21DX
    REAL(WP) :: TAU22F, TAU22B, DTAU22DY
    REAL(WP) :: TAU23F, TAU23B, DTAU23DZ
    REAL(WP) :: DTAU22DD, DWDZ, QR_R2, QT_R2, UR1, TAU24B, TAU24F
    REAL(WP) :: COE1, COE2, COE3 !, TESTVAL, TESTVAL1, TESTVAL2, TESTVAL3

    !REAL(WP) :: q1e, q1w, rhsl
    !REAL(WP) :: d11q2e
    !INTEGER(4) :: idr

    IF(iCase == ICHANL .OR. iCase == IBox3P) THEN
        RHS_io = 0.0_WP
        NYI = 1
        IF (MYID == 0 .AND. iCase == ICHANL) NYI = 2

        !COE1 = DXI * XND2CL * CVISC
        !COE3 = DZI * ZND2CL * CVISC

        COE1 = DXI * CVISC
        COE3 = DZI * CVISC

        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)
            COE2 = 2.0_WP * DYCI(JJ) * CVISC

            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)

                DO IC = 1, NCL1_io
                    IP = IPV_io(IC)
                    IM = IMV_io(IC)

                    !================ DY_TAU_22 =====================================
                    ! at (i, J,  k)
                    DVDY = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * DYFI(JJ)
                    TAU22F = 0.5_WP * ( Viscousity0(IC, JC, KC) + Viscousity0(IC, JC, KC) ) * ( DVDY - DivU_io(IC, JC, KC) )
                    ! at (i, J - 1, K)
                    DVDY = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JM, KC, 2) ) * DYFI(JJM)
                    TAU22B = 0.5_WP * ( Viscousity0(IC, JM, KC) + Viscousity0(IC, JM, KC) ) * ( DVDY - DivU_io(IC, JM, KC) )
                    ! at (i, J', k)
                    DTAU22DY = (TAU22F - TAU22B) * COE2

                    !================ DX_TAU_21 =====================================
                    ! at (i'+ 1, J', k)
                    !VIS21P = ( Viscousity(IC, JC, KC) + Viscousity(IP, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IP, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21P = MU_STG(IP, JC, KC, 1)
                    DVDX = ( Q_io(IP, JC, KC, 2) - Q_io(IC, JC, KC, 2) ) * DXI
                    DUDY = ( Q_io(IP, JC, KC, 1) - Q_io(IP, JM, KC, 1) ) * DYCI(JJ)
                    TAU21F = VIS21P * ( DVDX + DUDY )

                    ! at (i', J', k)
                    !VIS21C = ( Viscousity(IM, JC, KC) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IM, JM, KC) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21C = MU_STG(IC, JC, KC, 1)
                    DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI
                    DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ)
                    TAU21B = VIS21C * ( DVDX + DUDY )

                    ! at (i,   j', k)
                    DTAU21DX = (TAU21F - TAU21B) * COE1

                    !WRITE(*, *) 'G, J, I, K', JJ, KC, IC, Viscousity(IM, JM, KC), Viscousity(IC, JM, KC)

                    !================ DZ_TAU_23 =====================================
                    ! at (i, j', k'+ 1)
                    !VIS23P = ( Viscousity(IC, JC, KC) + Viscousity(IC, JC, KP) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IC, JM, KP) ) * YCL2ND_WFB(JJ)
                    VIS23P = MU_STG(IC, JC, KP, 3)
                    DVDZ = ( Q_io(IC, JC, KP, 2) - Q_io(IC, JC, KC, 2) ) * DZI
                    DWDY = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JM, KP, 3) ) * DYCI(JJ)
                    TAU23F = VIS23P * ( DVDZ + DWDY )

                    ! at (i, j', k'  )
                    !VIS23C = ( Viscousity(IC, JC, KM) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KM) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS23C = MU_STG(IC, JC, KC, 3)
                    DVDZ = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KM, 2) ) * DZI
                    DWDY = ( Q_io(IC, JC, KC, 3) - Q_io(IC, JM, KC, 3) ) * DYCI(JJ)
                    TAU23B = VIS23C * ( DVDZ + DWDY )

                    ! at (i, j', k   )
                    DTAU23DZ = (TAU23F - TAU23B) * COE3

                    !================ D_TAU_Y direction =================================
                    DPH_io(IC, JC, KC) = DPH_io(IC, JC, KC) + DTAU21DX + DTAU22DY + DTAU23DZ
                    !IF(JJ == 2 .AND. IC == 1 .AND. KC == 1) WRITE(*, '(A, 4ES13.5)') 'vIScy', &
                    !DTAU21DX, DTAU22DY, DTAU23DZ,RHS_io(IC, JC, KC)
                END DO
            END DO
        END DO
    END IF


    IF(iCase == iPIPEC .OR. iCase == iANNUL) THEN
        RHS_io = 0.0_WP
        NYI = 1
        IF (MYID == 0) NYI = 2

        !COE1 = DXI * XND2CL * CVISC
        !COE3 = DZI * ZND2CL * CVISC
        COE1 = DXI * CVISC
        COE3 = DZI * CVISC

        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)
            COE2 = 2.0_WP * DYCI(JJ) * CVISC

            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)
                KS = KSYM(KC)
                DO IC = 1, NCL1_io
                    IP = IPV_io(IC)
                    IM = IMV_io(IC)

                    !================ DY_TAU_22 =====================================
                    ! at (i, J,  k)
                    DVDY = ( Q_io(IC, JP, KC, 2) * RNDI1(JJP) - Q_io(IC, JC, KC, 2) * RNDI1(JJ) ) * DYFI(JJ)
                    TAU22F = 0.5_WP * ( Viscousity0(IC, JC, KC) + Viscousity0(IC, JC, KC) ) * &
                            ( DVDY - DivU_io(IC, JC, KC) ) / RCCI1(JJ)
                    ! at (i, J - 1, K)
                    IF(iCase == iPIPEC .AND. JJ == 2) THEN
                        Ur1  = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KS, 2) ) * 0.50_WP * RNDI1(JJ)
                        DVDY = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ)  - Ur1 ) * DYFI(JJM)
                    ELSE
                        DVDY = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ)  - Q_io(IC, JM, KC, 2) * RNDI1(JJM) ) * DYFI(JJM)
                    END IF
                    TAU22B = 0.5_WP * ( Viscousity0(IC, JM, KC) + Viscousity0(IC, JM, KC) ) * &
                             ( DVDY - DivU_io(IC, JM, KC) ) / RCCI1(JJM)
                    ! at (i, J', k)
                    DTAU22DY = (TAU22F - TAU22B) * COE2
                    !WRITE(*, '(3I3.1, 3ES13.5)') JJ, KC, IC,  TAU22F, TAU22B, DTAU22DY

                    !================ DX_TAU_21 =====================================
                    ! at (i'+ 1, J', k)
                    !VIS21P = ( Viscousity(IC, JC, KC) + Viscousity(IP, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IP, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21P = MU_STG(IP, JC, KC, 1)
                    DVDX = ( Q_io(IP, JC, KC, 2) - Q_io(IC, JC, KC, 2) ) * DXI * RNDI1(JJ)
                    DUDY = ( Q_io(IP, JC, KC, 1) - Q_io(IP, JM, KC, 1) ) * DYCI(JJ)
                    TAU21F = VIS21P * ( DVDX + DUDY )
                    ! at (i', J', k)
                    !VIS21C = ( Viscousity(IM, JC, KC) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IM, JM, KC) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21C = MU_STG(IC, JC, KC, 1)
                    DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI * RNDI1(JJ)
                    DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ)
                    TAU21B = VIS21C * ( DVDX + DUDY )

                    ! at (i,   j', k)
                    DTAU21DX = (TAU21F - TAU21B) * COE1 / RNDI1(JJ)

                    !================ DZ_TAU_23 =====================================
                    ! at (i, j', k'+ 1)
                    !VIS23P = ( Viscousity(IC, JC, KC) + Viscousity(IC, JC, KP) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IC, JM, KP) ) * YCL2ND_WFB(JJ)
                    VIS23P = MU_STG(IC, JC, KP, 3)
                    DVDZ = ( Q_io(IC, JC, KP, 2) - Q_io(IC, JC, KC, 2) ) * DZI * RNDI2(JJ)
                    DWDY = ( Q_io(IC, JC, KP, 3) * RCCI1(JJ) - Q_io(IC, JM, KP, 3) * RCCI1(JJM) ) * DYCI(JJ)
                    Qt_R2  = ( Q_io(IC, JC, KP, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + &
                    Q_io(IC, JM, KP, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ) ) * RNDI1(JJ)
                    TAU23F = VIS23P * ( DVDZ + DWDY - Qt_R2)

                    ! at (i, j', k'  )
                    !VIS23C = ( Viscousity(IC, JC, KM) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KM) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS23C = MU_STG(IC, JC, KC, 3)
                    DVDZ = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KM, 2) ) * DZI * RNDI2(JJ)
                    DWDY = ( Q_io(IC, JC, KC, 3) * RCCI1(JJ) - Q_io(IC, JM, KC, 3) * RCCI1(JJM) ) * DYCI(JJ)
                    Qt_R2  = ( Q_io(IC, JC, KC, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + &
                    Q_io(IC, JM, KC, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ) ) * RNDI1(JJ)
                    TAU23B = VIS23C * ( DVDZ + DWDY - Qt_R2)

                    ! at (i, j', k   )
                    DTAU23DZ = (TAU23F - TAU23B) * COE3

                    !================= ADDITIONAL ======================================
                    ! at (i, J,  k)
                    DWDZ = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * DZI * RCCI2(JJ)
                    QR_R2  = ( Q_io(IC, JP, KC, 2) * RNDI1(JJP) + Q_io(IC, JC, KC, 2) * RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
                    TAU24F = 0.5_WP * (Viscousity0(IC, JC, KC) + Viscousity0(IC, JC, KC)) * ( DWDZ + QR_R2 - DivU_io(IC, JC, KC) )

                    ! at (i, J - 1, k)
                    DWDZ = ( Q_io(IC, JM, KP, 3) - Q_io(IC, JM, KC, 3) ) * DZI * RCCI2(JJM)
                    IF(iCase == iPIPEC .AND. JJ == 2) THEN
                        Ur1  = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KS, 2) ) * 0.50_WP * RNDI1(JJ)
                        QR_R2  = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Ur1 ) * YND2CL * RCCI1(JJM)
                    ELSE
                        QR_R2  = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Q_io(IC, JM, KC, 2) * RNDI1(JJM) ) * YND2CL * RCCI1(JJM)
                    END IF
                    TAU24B = 0.5_WP * ( Viscousity0(IC, JM, KC) + Viscousity0(IC, JM, KC) ) * &
                             ( DWDZ + QR_R2 - DivU_io(IC, JM, KC) )

                    ! at (i, J',  k)
                    DTAU22DD = (TAU24F * YCL2ND_WFF(JJ) + TAU24B * YCL2ND_WFB(JJ)) * 2.0_WP * CVISC

                    !================ D_TAU_Y direction =================================
                    DPH_io(IC, JC, KC) = DPH_io(IC, JC, KC) + DTAU21DX + DTAU22DY + DTAU23DZ - DTAU22DD

                    !WRITE(*, '(3I3.1, 5ES13.5)') JJ, KC, IC,  RHS_io(IC, JC, KC), DTAU21DX, DTAU22DY, DTAU23DZ, - DTAU22DD
                END DO
            END DO
        END DO

    END IF

    !CALL DEBUG_WRT_LOCAL(DPH_io, 1, N2DO(MYID), 'vISy') ! test

    !        IDR = 2
    !        DO KC = 1, NCL3
    !            KM = KMV(KC)
    !            KP = KPV(KC)
    !            DO JC = 2, N2DO(MYID)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJ = JCL2G(JC)
    !                DO IC = 1, NCL1_io
    !                    IP = IPV_io(IC)
    !                    IM = IMV_io(IC)
    !                    RHSL  = CVISC * ( (Q_io(IP, JC, KC, IDR) - 2.0_WP * Q_io(IC, JC, KC, IDR) + Q_io(IM, JC, KC, IDR)) * DXQI +              &
    !                                    (Q_io(IC, JC, KP, IDR) - 2.0_WP * Q_io(IC, JC, KC, IDR) + Q_io(IC, JC, KM, IDR)) * DZQI * RNDI2(JJ) +      &
    !                                    (Q_io(IC, JP, KC, IDR) * APVR(JJ, IDR) +             &
    !                                     Q_io(IC, JC, KC, IDR) * ACVR(JJ, IDR) +             &
    !                                     Q_io(IC, JM, KC, IDR) * AMVR(JJ, IDR) )            &
    !                                 )
    !                    q1E = Q_io(IC, JC, KP, 3) * RCCI1(JJ) + Q_io(IC, JM, KP, 3) * RCCI1(JJM)
    !                    q1w= Q_io(IC, JC, KC, 3) * RCCI1(JJ) + Q_io(IC, JM, KC, 3) * RCCI1(JJM)
    !                    d11q2E = -(q1E -q1w) * DZI * RNDI1(JJ)
    !                    RHSL = d11q2e/ REN + RHSL
    !                    !IF(MYID == Npslv) WRITE(*, '(a, 3I3.1, 3ES13.5)') 'vIS_y', JJ, KC, IC, RHSL, RHS_io(IC, JC, KC), RHSL- RHS_io(IC, JC, KC)
    !                END DO
    !            END DO
    !        END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE VISCOUS_PAR_EXPLT_Y_io
    USE init_info
    USE MESH_INFO
    USE FLOW_INFO
    USE THERMAL_INFO
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP, KS
    INTEGER(4) :: NYI

    REAL(WP) :: VIS21C, VIS21P
    REAL(WP) :: VIS23C, VIS23P
    REAL(WP) :: DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ
    REAL(WP) :: TAU21F_EX, TAU21B_EX, DTAU21DX_EX
    REAL(WP) :: TAU22F_EX, TAU22B_EX, DTAU22DY_EX
    REAL(WP) :: TAU23F_EX, TAU23B_EX, DTAU23DZ_EX
    REAL(WP) :: TAU21F_IM, TAU21B_IM, DTAU21DX_IM
    REAL(WP) :: TAU22F_IM, TAU22B_IM, DTAU22DY_IM
    REAL(WP) :: TAU23F_IM, TAU23B_IM, DTAU23DZ_IM
    REAL(WP) :: TAU21F_AD, TAU21B_AD, DTAU21DX_AD
    REAL(WP) :: TAU22F_AD, TAU22B_AD, DTAU22DY_AD
    REAL(WP) :: TAU23F_AD, TAU23B_AD, DTAU23DZ_AD
    REAL(WP) :: DTAU22DD, QR_R2, QT_R2, UR1, TAU24B, TAU24F
    REAL(WP) :: COE1, COE2, COE3 !, TESTVAL, TESTVAL1, TESTVAL2, TESTVAL3

    !REAL(WP) :: q1e, q1w, rhsl
    !REAL(WP) :: d11q2e
    !INTEGER(4) :: idr

    IF(iCase == ICHANL) THEN
        RHS_io = 0.0_WP
        NYI = 1
        IF (MYID == 0) NYI = 2

        !COE1 = DXI * XND2CL * CVISC
        !COE3 = DZI * ZND2CL * CVISC

        COE1 = DXI * CVISC
        COE3 = DZI * CVISC

        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)
            COE2 = -2.0_WP / 3.0_WP * DYCI(JJ) * CVISC

            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)

                DO IC = 1, NCL1_io
                    IP = IPV_io(IC)
                    IM = IMV_io(IC)

                    !================ DY_TAU_22 === (PARTLY) ================================
                    ! at (i, J,  k)
                    DUDX = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
                    DVDY = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * DYFI(JJ)
                    DWDZ = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * DZI
                    TAU22F_EX = Viscousity(IC, JC, KC) * (DUDX + DWDZ)
                    TAU22F_IM = Viscousity(IC, JC, KC) * (DVDY)

                    DVDY = ( G_io(IC, JP, KC, 2) * DRHOI_STG(IC, JP, KC, 2) - &
                             G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) ) * DYFI(JJ)
                    TAU22F_AD = Viscousity(IC, JC, KC) * (DVDY)

                    ! at (i, J - 1,  k)
                    DUDX = ( Q_io(IP, JM, KC, 1) - Q_io(IC, JM, KC, 1) ) * DXI
                    DVDY = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JM, KC, 2) ) * DYFI(JJM)
                    DWDZ = ( Q_io(IC, JM, KP, 3) - Q_io(IC, JM, KC, 3) ) * DZI
                    TAU22B_EX = Viscousity(IC, JM, KC) * (DUDX + DWDZ)
                    TAU22B_IM = Viscousity(IC, JM, KC) * (DVDY)

                    DVDY = ( G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) - &
                             G_io(IC, JM, KC, 2) * DRHOI_STG(IC, JM, KC, 2) ) * DYFI(JJM)
                    TAU22B_AD = Viscousity(IC, JM, KC) * (DVDY)

                    ! at (i, J', k)
                    DTAU22DY_EX = (TAU22F_EX - TAU22B_EX) * COE2
                    DTAU22DY_IM = (TAU22F_IM - TAU22B_IM) * COE2 * (-2.0_WP)

                    DTAU22DY_AD = (TAU22F_AD - TAU22B_AD) * COE2 * (-2.0_WP)

                    !================ DX_TAU_21 =====================================
                    ! at (i'+ 1, J', k)
                    !VIS21P = ( Viscousity(IC, JC, KC) + Viscousity(IP, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IP, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21P = MU_STG(IP, JC, KC, 1)
                    DVDX = ( Q_io(IP, JC, KC, 2) - Q_io(IC, JC, KC, 2) ) * DXI
                    DUDY = ( Q_io(IP, JC, KC, 1) - Q_io(IP, JM, KC, 1) ) * DYCI(JJ)
                    TAU21F_EX = VIS21P * ( DUDY )
                    TAU21F_IM = VIS21P * ( DVDX )

                    DVDX = ( G_io(IP, JC, KC, 2) * DRHOI_STG(IP, JC, KC, 2) - &
                             G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) ) * DXI
                    TAU21F_AD = VIS21P * ( DVDX )

                    ! at (i', J', k)
                    !VIS21C = ( Viscousity(IM, JC, KC) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IM, JM, KC) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS21C = MU_STG(IC, JC, KC, 1)
                    DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI
                    DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ)
                    TAU21B_EX = VIS21C * ( DUDY )
                    TAU21B_IM = VIS21C * ( DVDX )

                    DVDX = ( G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) - &
                             G_io(IM, JC, KC, 2) * DRHOI_STG(IM, JC, KC, 2) ) * DXI
                    TAU21B_AD = VIS21C * ( DVDX )

                    ! at (i,   j', k)
                    DTAU21DX_EX = (TAU21F_EX - TAU21B_EX) * COE1
                    DTAU21DX_IM = (TAU21F_IM - TAU21B_IM) * COE1

                    DTAU21DX_AD = (TAU21F_AD - TAU21B_AD) * COE1

                    !WRITE(*, *) 'G, J, I, K', JJ, KC, IC, Viscousity(IM, JM, KC), Viscousity(IC, JM, KC)

                    !================ DZ_TAU_23 =====================================
                    ! at (i, j', k'+ 1)
                    !VIS23P = ( Viscousity(IC, JC, KC) + Viscousity(IC, JC, KP) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KC) + Viscousity(IC, JM, KP) ) * YCL2ND_WFB(JJ)
                    VIS23P = MU_STG(IC, JC, KP, 3)
                    DVDZ = ( Q_io(IC, JC, KP, 2) - Q_io(IC, JC, KC, 2) ) * DZI
                    DWDY = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JM, KP, 3) ) * DYCI(JJ)
                    TAU23F_EX = VIS23P * ( DWDY )
                    TAU23F_IM = VIS23P * ( DVDZ )

                    DVDZ = ( G_io(IC, JC, KP, 2) * DRHOI_STG(IC, JC, KP, 2) - &
                             G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) ) * DZI
                    TAU23F_AD = VIS23P * ( DVDZ )

                    ! at (i, j', k'  )
                    !VIS23C = ( Viscousity(IC, JC, KM) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( Viscousity(IC, JM, KM) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
                    VIS23C = MU_STG(IC, JC, KC, 3)
                    DVDZ = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KM, 2) ) * DZI
                    DWDY = ( Q_io(IC, JC, KC, 3) - Q_io(IC, JM, KC, 3) ) * DYCI(JJ)
                    TAU23B_EX = VIS23C * ( DWDY )
                    TAU23B_IM = VIS23C * ( DVDZ )

                    DVDZ = ( G_io(IC, JC, KC, 2) * DRHOI_STG(IC, JC, KC, 2) - &
                             G_io(IC, JC, KM, 2) * DRHOI_STG(IC, JC, KM, 2) ) * DZI
                    TAU23B_AD = VIS23C * ( DVDZ )

                    ! at (i, j', k   )
                    DTAU23DZ_EX = (TAU23F_EX - TAU23B_EX) * COE3
                    DTAU23DZ_IM = (TAU23F_IM - TAU23B_IM) * COE3

                    DTAU23DZ_AD = (TAU23F_AD - TAU23B_AD) * COE3

                    !================ D_TAU_Y direction =================================
                    DPH_io(IC, JC, KC) = DTAU21DX_EX + DTAU22DY_EX + DTAU23DZ_EX + DPH_io(IC, JC, KC)
                    RHS_io(IC, JC, KC) = DTAU21DX_IM + DTAU22DY_IM + DTAU23DZ_IM
                    DIVU_io(IC, JC, KC) = DTAU21DX_AD + DTAU22DY_AD + DTAU23DZ_AD
                    !WRITE(*, *) DTAU21DX_AD, DTAU22DY_AD, DTAU23DZ_AD
                    !IF(JJ == 2 .AND. IC == 1 .AND. KC == 1) WRITE(*, '(A, 4ES13.5)') 'vIScy', &
                    !DTAU21DX, DTAU22DY, DTAU23DZ,RHS_io(IC, JC, KC)
                END DO
            END DO
        END DO
    END IF

    ! for cylinderICal coordinates, to DO....
    !        IF(iCase == iPIPEC .OR. iCase == iANNUL) THEN
    !            RHS_io = 0.0_WP
    !            NYI = 1
    !            IF (MYID == 0) NYI = 2

    !            COE1 = DXI * XND2CL * CVISC
    !            COE3 = DZI * ZND2CL * CVISC

    !            DO JC = NYI, N2DO(MYID)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJ = JCL2G(JC)
    !                JJM = JGMV(JJ)
    !                JJP = JGPV(JJ)
    !                COE2 = 2.0_WP * DYCI(JJ) * CVISC

    !                DO KC = 1, NCL3
    !                    KM = KMV(KC)
    !                    KP = KPV(KC)
    !                    KS = KSYM(KC)
    !                    DO IC = 1, NCL1_io
    !                        IP = IPV_io(IC)
    !                        IM = IMV_io(IC)

    !                        !================ DY_TAU_22 == (PARTLY) ================================
    !                        ! at (i, J,  k)
    !                        DVDY = ( Q_io(IC, JP, KC, 2) * RNDI1(JJP) - Q_io(IC, JC, KC, 2) * RNDI1(JJ) ) * DYFI(JJ)
    !                        TAU22F = Viscousity(IC, JC, KC) * ( DVDY - DivU_io(IC, JC, KC) ) / RCCI1(JJ)
    !                        ! at (i, J - 1, K)
    !                        IF(iCase == iPIPEC .AND. JJ == 2) THEN
    !                            Ur1  = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KS, 2) ) * 0.50_WP * RNDI1(JJ)
    !                            DVDY = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ)  - Ur1 ) * DYFI(JJM)
    !                        ELSE
    !                            DVDY = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ)  - Q_io(IC, JM, KC, 2) * RNDI1(JJM) ) * DYFI(JJM)
    !                        END IF
    !                        TAU22B = Viscousity(IC, JM, KC) * ( DVDY - DivU_io(IC, JM, KC) ) / RCCI1(JJM)
    !                        ! at (i, J', k)
    !                        DTAU22DY = (TAU22F - TAU22B) * COE2
    !                        !WRITE(*, '(3I3.1, 3ES13.5)') JJ, KC, IC,  TAU22F, TAU22B, DTAU22DY

    !                        ! at (i, J,  k)
    !                        DUDX = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
    !                        DWDZ = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * DZI
    !                        TAU22F = Viscousity(IC, JC, KC) * (DUDX + DWDZ)

    !                        ! at (i, J - 1,  k)
    !                        DUDX = ( Q_io(IP, JM, KC, 1) - Q_io(IC, JM, KC, 1) ) * DXI
    !                        DWDZ = ( Q_io(IC, JM, KP, 3) - Q_io(IC, JM, KC, 3) ) * DZI
    !                        TAU22B = Viscousity(IC, JM, KC) * (DUDX + DWDZ)

    !                        ! at (i, J', k)
    !                        DTAU22DY = (TAU22F - TAU22B) * COE2

    !                        !================ DX_TAU_21 =====================================
    !                        ! at (i'+ 1, J', k)
    !                        VIS21P = ( Viscousity(IC, JC, KC) + Viscousity(IP, JC, KC) ) * YCL2ND_WFF(JJ) + &
    !                                 ( Viscousity(IC, JM, KC) + Viscousity(IP, JM, KC) ) * YCL2ND_WFB(JJ)
    !                        DVDX = ( Q_io(IP, JC, KC, 2) - Q_io(IC, JC, KC, 2) ) * DXI * RNDI1(JJ)
    !                        DUDY = ( Q_io(IP, JC, KC, 1) - Q_io(IP, JM, KC, 1) ) * DYCI(JJ)
    !                        TAU21F = VIS21P * ( DVDX + DUDY )
    !                        ! at (i', J', k)
    !                        VIS21C = ( Viscousity(IM, JC, KC) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
    !                                 ( Viscousity(IM, JM, KC) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
    !                        DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI * RNDI1(JJ)
    !                        DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ)
    !                        TAU21B = VIS21C * ( DVDX + DUDY )

    !                        ! at (i,   j', k)
    !                        DTAU21DX = (TAU21F - TAU21B) * COE1 / RNDI1(JJ)

    !                        !================ DZ_TAU_23 =====================================
    !                        ! at (i, j', k'+ 1)
    !                        VIS23P = ( Viscousity(IC, JC, KC) + Viscousity(IC, JC, KP) ) * YCL2ND_WFF(JJ) + &
    !                                 ( Viscousity(IC, JM, KC) + Viscousity(IC, JM, KP) ) * YCL2ND_WFB(JJ)
    !                        DVDZ = ( Q_io(IC, JC, KP, 2) - Q_io(IC, JC, KC, 2) ) * DZI * RNDI2(JJ)
    !                        DWDY = ( Q_io(IC, JC, KP, 3) * RCCI1(JJ) - Q_io(IC, JM, KP, 3) * RCCI1(JJM) ) * DYCI(JJ)
    !                        Qt_R2  = ( Q_io(IC, JC, KP, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + &
    !                                   Q_io(IC, JM, KP, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ) ) * RNDI1(JJ)
    !                        TAU23F = VIS23P * ( DVDZ + DWDY - Qt_R2)

    !                        ! at (i, j', k'  )
    !                        VIS23C = ( Viscousity(IC, JC, KM) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
    !                                 ( Viscousity(IC, JM, KM) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
    !                        DVDZ = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KM, 2) ) * DZI * RNDI2(JJ)
    !                        DWDY = ( Q_io(IC, JC, KC, 3) * RCCI1(JJ) - Q_io(IC, JM, KC, 3) * RCCI1(JJM) ) * DYCI(JJ)
    !                        Qt_R2  = ( Q_io(IC, JC, KC, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + &
    !                                   Q_io(IC, JM, KC, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ) ) * RNDI1(JJ)
    !                        TAU23B = VIS23C * ( DVDZ + DWDY - Qt_R2)

    !                        ! at (i, j', k   )
    !                        DTAU23DZ = (TAU23F - TAU23B) * COE3

    !                        !================= ADDITIONAL ======================================
    !                        ! at (i, J,  k)
    !                        DWDZ = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * DZI * RCCI2(JJ)
    !                        QR_R2  = ( Q_io(IC, JP, KC, 2) * RNDI1(JJP) + Q_io(IC, JC, KC, 2) * RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
    !                        TAU24F = Viscousity(IC, JC, KC) * ( DWDZ + QR_R2 - DivU_io(IC, JC, KC) )

    !                        ! at (i, J - 1, k)
    !                        DWDZ = ( Q_io(IC, JM, KP, 3) - Q_io(IC, JM, KC, 3) ) * DZI * RCCI2(JJM)
    !                        IF(iCase == iPIPEC .AND. JJ == 2) THEN
    !                            Ur1  = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KS, 2) ) * 0.50_WP * RNDI1(JJ)
    !                            QR_R2  = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Ur1 ) * YND2CL * RCCI1(JJM)
    !                        ELSE
    !                            QR_R2  = ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Q_io(IC, JM, KC, 2) * RNDI1(JJM) ) * YND2CL * RCCI1(JJM)
    !                        END IF
    !                        TAU24B = Viscousity(IC, JM, KC) * ( DWDZ + QR_R2 - DivU_io(IC, JM, KC) )

    !                        ! at (i, J',  k)
    !                        DTAU22DD = (TAU24F * YCL2ND_WFF(JJ) + TAU24B * YCL2ND_WFB(JJ)) * 2.0_WP * CVISC

    !                        !================ D_TAU_Y direction =================================
    !                        RHS_io(IC, JC, KC) = DTAU21DX + DTAU22DY + DTAU23DZ - DTAU22DD
    !                        !RHS_io(IC, JC, KC) = DTAU22DY
    !                        !WRITE(*, '(3I3.1, 5ES13.5)') JJ, KC, IC,  RHS_io(IC, JC, KC), DTAU21DX, DTAU22DY, DTAU23DZ, - DTAU22DD
    !                    END DO
    !                END DO
    !            END DO

    !        END IF

    !!        IDR = 2
    !!        DO KC = 1, NCL3
    !!            KM = KMV(KC)
    !!            KP = KPV(KC)
    !!            DO JC = 2, N2DO(MYID)
    !!                JM = JLMV(JC)
    !!                JP = JLPV(JC)
    !!                JJ = JCL2G(JC)
    !!                DO IC = 1, NCL1_io
    !!                    IP = IPV_io(IC)
    !!                    IM = IMV_io(IC)
    !!                    RHSL  = CVISC * ( (Q_io(IP, JC, KC, IDR) - 2.0_WP * Q_io(IC, JC, KC, IDR) + Q_io(IM, JC, KC, IDR)) * DXQI +              &
    !!                                    (Q_io(IC, JC, KP, IDR) - 2.0_WP * Q_io(IC, JC, KC, IDR) + Q_io(IC, JC, KM, IDR)) * DZQI * RNDI2(JJ) +      &
    !!                                    (Q_io(IC, JP, KC, IDR) * APVR(JJ, IDR) +             &
    !!                                     Q_io(IC, JC, KC, IDR) * ACVR(JJ, IDR) +             &
    !!                                     Q_io(IC, JM, KC, IDR) * AMVR(JJ, IDR) )            &
    !!                                 )
    !!                    q1E = Q_io(IC, JC, KP, 3) * RCCI1(JJ) + Q_io(IC, JM, KP, 3) * RCCI1(JJM)
    !!                    q1w= Q_io(IC, JC, KC, 3) * RCCI1(JJ) + Q_io(IC, JM, KC, 3) * RCCI1(JJM)
    !!                    d11q2E = -(q1E -q1w) * DZI * RNDI1(JJ)
    !!                    RHSL = d11q2e/ REN + RHSL
    !!                    !IF(MYID == Npslv) WRITE(*, '(a, 3I3.1, 3ES13.5)') 'vIS_y', JJ, KC, IC, RHSL, RHS_io(IC, JC, KC), RHSL- RHS_io(IC, JC, KC)
    !!                END DO
    !!            END DO
    !!        END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CHECK_GradP_ON_WALL
    USE init_info
    USE MESH_INFO
    USE FLOW_INFO
    USE THERMAL_INFO
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP, KS
    INTEGER(4) :: NYI

    REAL(WP) :: VIS21C, VIS21P
    REAL(WP) :: VIS23C, VIS23P
    REAL(WP) :: DVDY, DVDX, DUDY, DVDZ, DWDY
    REAL(WP) :: TAU21F, TAU21B, DTAU21DX
    REAL(WP) :: TAU22F, TAU22B, DTAU22DY
    REAL(WP) :: TAU23F, TAU23B, DTAU23DZ
    REAL(WP) :: DTAU22DD, DWDZ, QR_R2, QT_R2, UR1, TAU24B, TAU24F, VISYALL
    REAL(WP) :: COE1, COE2, COE3 !, TESTVAL, TESTVAL1, TESTVAL2, TESTVAL3

    !REAL(WP) :: q1e, q1w, rhsl
    !REAL(WP) :: d11q2e
    !INTEGER(4) :: idr

    IF(MYID /= 0) RETURN

    COE1 = DXI * XND2CL * CVISC
    COE3 = DZI * ZND2CL * CVISC

    JC = 1
    JM = JLMV(JC)
    JP = JLPV(JC)
    JJ = JCL2G(JC)
    JJM = JGMV(JJ)
    JJP = JGPV(JJ)
    COE2 = 2.0_WP * DYCI(JJ) * CVISC

    DO KC = 1, NCL3
        KM = KMV(KC)
        KP = KPV(KC)

        DO IC = 1, NCL1_io
            IP = IPV_io(IC)
            IM = IMV_io(IC)

            !================ DY_TAU_22 =====================================
            ! at (i, J,  k)
            DVDY = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * DYFI(JJ)
            TAU22F = Viscousity(IC, JC, KC) * ( DVDY - DivU_io(IC, JC, KC) )
            ! at (i, J - 1, K)
            DVDY = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JM, KC, 2) ) * DYFI(JJM)
            TAU22B = Viscousity(IC, JM, KC) * ( DVDY - DivU_io(IC, JM, KC) )
            ! at (i, J', k)
            DTAU22DY = (TAU22F - TAU22B) * COE2

            !================ DX_TAU_21 =====================================
            ! at (i'+ 1, J', k)
            VIS21P = ( Viscousity(IC, JC, KC) + Viscousity(IP, JC, KC) ) * YCL2ND_WFF(JJ) + &
            ( Viscousity(IC, JM, KC) + Viscousity(IP, JM, KC) ) * YCL2ND_WFB(JJ)
            DVDX = ( Q_io(IP, JC, KC, 2) - Q_io(IC, JC, KC, 2) ) * DXI
            DUDY = ( Q_io(IP, JC, KC, 1) - Q_io(IP, JM, KC, 1) ) * DYCI(JJ)
            TAU21F = VIS21P * ( DVDX + DUDY )

            ! at (i', J', k)
            VIS21C = ( Viscousity(IM, JC, KC) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
            ( Viscousity(IM, JM, KC) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
            DVDX = ( Q_io(IC, JC, KC, 2) - Q_io(IM, JC, KC, 2) ) * DXI
            DUDY = ( Q_io(IC, JC, KC, 1) - Q_io(IC, JM, KC, 1) ) * DYCI(JJ)
            TAU21B = VIS21C * ( DVDX + DUDY )

            ! at (i,   j', k)
            DTAU21DX = (TAU21F - TAU21B) * COE1

            !WRITE(*, *) 'G, J, I, K', JJ, KC, IC, Viscousity(IM, JM, KC), Viscousity(IC, JM, KC)

            !================ DZ_TAU_23 =====================================
            ! at (i, j', k'+ 1)
            !VIS23P = ( Viscousity(IC, JC, KC) + Viscousity(IC, JC, KP) ) * YCL2ND_WFF(JJ) + &
            !         ( Viscousity(IC, JM, KC) + Viscousity(IC, JM, KP) ) * YCL2ND_WFB(JJ)
            VIS23P = MU_STG(IC, JC, KP, 3)
            DVDZ = ( Q_io(IC, JC, KP, 2) - Q_io(IC, JC, KC, 2) ) * DZI
            DWDY = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JM, KP, 3) ) * DYCI(JJ)
            TAU23F = VIS23P * ( DVDZ + DWDY )

            ! at (i, j', k'  )
            !VIS23C = ( Viscousity(IC, JC, KM) + Viscousity(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
            !         ( Viscousity(IC, JM, KM) + Viscousity(IC, JM, KC) ) * YCL2ND_WFB(JJ)
            VIS23C = MU_STG(IC, JC, KC, 3)
            DVDZ = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KM, 2) ) * DZI
            DWDY = ( Q_io(IC, JC, KC, 3) - Q_io(IC, JM, KC, 3) ) * DYCI(JJ)
            TAU23B = VIS23C * ( DVDZ + DWDY )

            ! at (i, j', k   )
            DTAU23DZ = (TAU23F - TAU23B) * COE3

            !================ D_TAU_Y direction =================================
            VISYALL = DTAU22DY + DTAU21DX + DTAU23DZ
            WRITE(logflg_pg, '(A, 2I3.1, 4ES18.10)') 'G(P)_WALL', KC, IC, DTAU22DY, DTAU21DX, DTAU23DZ, VISYALL
        END DO
    END DO




    RETURN
END SUBROUTINE
