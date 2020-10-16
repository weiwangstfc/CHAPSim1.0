!**********************************************************************************************************************************
!> @brief
!>        Calculate the convection term in Y direction.
!> @details
!> SUBROUTINE: CONVECTION_Y_tg (in MYID = all)
!> SUBROUTINE: CONVECTION_Y_io(in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 04/ 2014- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CONVECTION_Y_tg
    USE flow_info
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP
    INTEGER(4) :: KS
    INTEGER(4) :: NYI
    REAL(WP) :: H21, H22, H23
    REAL(WP) :: q2s1
    REAL(WP) :: h23n
    REAL(WP) :: q1e, q1w
    REAL(WP) :: d11q2e
    REAL(WP) :: COE1,COE2, COE30, COE31

    DPH_tg(:, :, :) = 0.0_WP  ! to store the convection term in y direction.

    NYI = 1
    IF (MYID == 0) NYI = 2

    COE1  = DXI * XND2CL
    COE30 = DZI * ZND2CL
    DO JC = NYI, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJM = JGMV(JJ)
        JJP = JGPV(JJ)
        COE2  = DYCI(JJ) * YND2CL * YND2CL
        COE31 = COE30 * RNDI1(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)
            DO IC = 1, NCL1_tg

                IP = IPV_tg(IC)
                IM = IMV_tg(IC)
                ! I, J', K = (I'+ 1, J', K)  vs  (I', J', K)
                H21 = ( ( YCL2ND_WFF(JJ) * Q_tg(IP, JC, KC, 1) +                      &
                YCL2ND_WFB(JJ) * Q_tg(IP, JM, KC, 1) ) *                   &
                ( Q_tg(IP, JC, KC, 2) + Q_tg(IC, JC, KC, 2) ) -                     &
                ( YCL2ND_WFF(JJ) * Q_tg(IC, JC, KC, 1) +                      &
                YCL2ND_WFB(JJ) * Q_tg(IC, JM, KC, 1) ) *                  &
                ( Q_tg(IC, JC, KC, 2) + Q_tg(IM, JC, KC, 2) ) ) * COE1
                IF (iCase == IPIPEC .AND. JJ == 2) THEN
                    KS = KSYM(KC)
                    q2s1 = (Q_tg(IC, JC, KC, 2) - Q_tg(IC, JC, KS, 2)) * 0.50_WP * RNDI1(JJ)
                    h22 = ((Q_tg(IC, JP, KC, 2) * RNDI1(JJP) + Q_tg(IC, JC, KC, 2) * RNDI1(JJ)) *     &
                    (Q_tg(IC, JP, KC, 2) + Q_tg(IC, JC, KC, 2)) -                  &
                    (Q_tg(IC, JC, KC, 2) * RNDI1(JJ) + Q2s1) *                   &
                    (Q_tg(IC, JC, KC, 2) + Q_tg(IC, JM, KC, 2)) ) * COE2
                ELSE
                    ! I, J', K = (I, J, K)  vs  (I, J - 1, K)
                    H22 = ((Q_tg(IC, JP, KC, 2) * RNDI1(JJP) + Q_tg(IC, JC, KC, 2) * RNDI1(JJ)) * &
                    (Q_tg(IC, JP, KC, 2) + Q_tg(IC, JC, KC, 2)) -                  &
                    (Q_tg(IC, JC, KC, 2) * RNDI1(JJ) + Q_tg(IC, JM, KC, 2) * RNDI1(JJM)) * &
                    (Q_tg(IC, JC, KC, 2) + Q_tg(IC, JM, KC, 2)) ) * COE2
                END IF
                ! I, J', K = (I, J', K'+ 1)  vs  (I, J', K')
                H23 = ( ( YCL2ND_WFF(JJ) * Q_tg(IC, JC, KP, 3) * RCCI1(JJ) +      &
                YCL2ND_WFB(JJ) * Q_tg(IC, JM, KP, 3) * RCCI1(JJM) ) * &
                ( Q_tg(IC, JC, KP, 2) + Q_tg(IC, JC, KC, 2) ) - &
                ( YCL2ND_WFF(JJ) * Q_tg(IC, JC, KC, 3) * RCCI1(JJ) +      &
                YCL2ND_WFB(JJ) * Q_tg(IC, JM, KC, 3) * RCCI1(JJM) ) * &
                ( Q_tg(IC, JC, KC, 2) + Q_tg(IC, JC, KM, 2) ) ) * COE31        !@
                DPH_tg(IC, JC, KC) = -(H21 + H22 + H23)   !H for momentum y direction.
            END DO
        END DO
    END DO


    IF (iCase == IPIPEC .OR. iCase == IANNUL) THEN            !@   Centripetal force
        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)
            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)
                DO IC = 1, NCL1_tg
                    !                        q1E = Q_tg(IC, JC, KP, 3) * RCCI1(JJ) + Q_tg(IC, JM, KP, 3) * RCCI1(JJM)
                    !                        q1w= Q_tg(IC, JC, KC, 3) * RCCI1(JJ) + Q_tg(IC, JM, KC, 3) * RCCI1(JJM)
                    !                        h23N = ( (q1E + Q1w) * 0.250_WP )**2
                    !                        d11q2E = -(q1E -q1w) * DZI * RNDI1(JJ)
                    !                        DPH_tg(IC, JC, KC) = DPH_tg(IC, JC, KC) + H23N + D11q2e/ REN
                    q1e = Q_tg(IC, JC, KP, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + Q_tg(IC, JM, KP, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)
                    q1w = Q_tg(IC, JC, KC, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + Q_tg(IC, JM, KC, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)
                    h23n = ( (q1e + Q1w) * ZND2CL )**2
                    d11q2e = -2.0_WP * (q1e - q1w) * DZI * RNDI1(JJ)
                    DPH_tg(IC, JC, KC) = DPH_tg(IC, JC, KC) + H23n + d11q2e / REN
                END DO
            END DO
        END DO
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CONVECTION_Y_io
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP, KS
    INTEGER(4) :: NYI

    REAL(WP) :: H21, H22, H23, H23N
    REAL(WP) :: H21F
    REAL(WP) :: H21B
    REAL(WP) :: H22F
    REAL(WP) :: H22B
    REAL(WP) :: H23F
    REAL(WP) :: H23B
    !REAL(WP) :: H24, H24F, H24B, H24F1, H24F2, H24B1, H24B2
    REAL(WP) :: COE1, COE2, COE30, COE31!, COE41, COE42
    REAL(WP) :: Ur1, q1e, q1w, g1e, g1w

    DPH_io(:, :, :) = 0.0_WP  ! to store the convection term in y direction.

    NYI = 1
    IF (MYID == 0 .AND. iCase /= IBox3P) NYI = 2

    COE1  = DXI * XND2CL
    COE30 = DZI * ZND2CL
    DO JC = NYI, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJM = JGMV(JJ)
        JJP = JGPV(JJ)
        COE2 = DYCI(JJ) * YND2CL * YND2CL
        COE31 = COE30 * RNDI1(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)

            DO IC = 1, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)

                ! \frac{\pARtial {\rho u v}}{\pARtial x}_{i, J', K}
                !{i'+ 1, J', K}
                H21F = ( YCL2ND_WFF(JJ) * G_io(IP, JC, KC, 1) +   &
                YCL2ND_WFB(JJ) * G_io(IP, JM, KC, 1) ) * &
                ( Q_io(IP, JC, KC, 2) + Q_io(IC, JC, KC, 2) )
                !{i', J', K}
                H21B = ( YCL2ND_WFF(JJ) * G_io(IC, JC, KC, 1) +   &
                YCL2ND_WFB(JJ) * G_io(IC, JM, KC, 1) ) * &
                ( Q_io(IC, JC, KC, 2) + Q_io(IM, JC, KC, 2) )
                !{i, J', K}
                H21  = ( H21F - H21B) * COE1

                ! \frac{\pARtial {\rho v v}}{\pARtial y}_{i, J', K}
                !{i, J, K}
                H22F = ( G_io(IC, JP, KC, 2) + G_io(IC, JC, KC, 2) ) * &
                ( Q_io(IC, JP, KC, 2) * RNDI1(JJP) + Q_io(IC, JC, KC, 2) * RNDI1(JJ) )
                !{i, J - 1, K}
                IF(iCase == IPIPEC .AND. JJ == 2) THEN
                    KS = KSYM(KC)
                    Ur1  = ( Q_io(IC, JC, KC, 2) - Q_io(IC, JC, KS, 2) ) * 0.50_WP * RNDI1(JJ)
                    H22B = ( G_io(IC, JC, KC, 2) + G_io(IC, JM, KC, 2) ) * &
                    ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Ur1 )
                ELSE
                    H22B = ( G_io(IC, JC, KC, 2) + G_io(IC, JM, KC, 2) ) * &
                    ( Q_io(IC, JC, KC, 2) * RNDI1(JJ) + Q_io(IC, JM, KC, 2) * RNDI1(JJM) )
                END IF
                !{i, J', K}
                H22 = ( H22F - H22B) * COE2

                ! \frac{\pARtial {\rho w v}}{\pARtial x}_{i, J', K}
                !{i, J', K'+ 1}
                H23F = ( YCL2ND_WFF(JJ) * G_io(IC, JC, KP, 3) * RCCI1(JJ) +   &
                YCL2ND_WFB(JJ) * G_io(IC, JM, KP, 3) * RCCI1(JJM) ) * &
                ( Q_io(IC, JC, KP, 2) + Q_io(IC, JC, KC, 2) ) !* RNDI1(JJ)
                !{i, J', K'}
                H23B = ( YCL2ND_WFF(JJ) * G_io(IC, JC, KC, 3) * RCCI1(JJ) +   &
                YCL2ND_WFB(JJ) * G_io(IC, JM, KC, 3) * RCCI1(JJM) ) * &
                ( Q_io(IC, JC, KC, 2) + Q_io(IC, JC, KM, 2) ) !* RNDI1(JJ)
                !{i, J', K}
                H23 = ( H23F - H23B) * COE31

                DPH_io(IC, JC, KC) = -(H21 + H22 + H23)   !H for momentum y direction.
                !IF(JJ == 2 .AND. IC == 1 .AND. KC == 1) WRITE(*, '(A, 4ES13.5)') 'convy',H21,H22,H23, DPH_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !@   Centripetal force
    IF (iCase == IPIPEC .OR. iCase == IANNUL) THEN            !@   Centripetal force
        DO JC = NYI, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)
            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)
                DO IC = 1, NCL1_io
                    q1e = Q_io(IC, JC, KP, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + Q_io(IC, JM, KP, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)
                    q1w = Q_io(IC, JC, KC, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + Q_io(IC, JM, KC, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)
                    g1e = G_io(IC, JC, KP, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + G_io(IC, JM, KP, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)
                    g1w = G_io(IC, JC, KC, 3) * RCCI1(JJ) * YCL2ND_WFF(JJ) + G_io(IC, JM, KC, 3) * RCCI1(JJM) * YCL2ND_WFB(JJ)

                    !h23N = ( (q1E + Q1w) * 0.250_WP )**2
                    h23n = (q1e + Q1w) * ZND2CL * (g1e + g1w) * ZND2CL
                    DPH_io(IC, JC, KC) = DPH_io(IC, JC, KC) + H23n
                END DO
            END DO
        END DO
    END IF

    !        !@   Centripetal force
    !        IF (iCase == IPIPEC .OR. iCase == IANNUL) THEN

    !            DO JC = NYI, N2DO(MYID)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJ = JCL2G(JC)
    !                JJM = JGMV(JJ)
    !                JJP = JGPV(JJ)
    !                COE41 = ZND2CL * ZND2CL * RCCI2(JJ)
    !                COE42 = ZND2CL * ZND2CL * RCCI2(JJM)
    !                DO KC = 1, NCL3
    !                    KM = KMV(KC)
    !                    KP = KPV(KC)

    !                    DO IC = 1, NCL1_io
    !                        ! AT I, J,  K
    !                        H24F1 = ( G_io(IC, JC, KP, 3) + G_io(IC, JC, KC, 3) ) !* ZND2CL * RCCI1(JJ)
    !                        H24F2 = ( Q_io(IC, JC, KP, 3) + Q_io(IC, JC, KC, 3) ) !* ZND2CL * RCCI1(JJ)
    !                        H24F  = H24F1 * H24F2 * COE41
    !                        ! AT I, J - 1, K
    !                        H24B1 = ( G_io(IC, JM, KP, 3) + G_io(IC, JM, KC, 3) ) !* ZND2CL * RCCI1(JJM)
    !                        H24B2 = ( Q_io(IC, JM, KP, 3) + Q_io(IC, JM, KC, 3) ) !* ZND2CL * RCCI1(JJM)
    !                        H24B  = H24B1 * H24B2 * COE42

    !                        ! AT I, J', K
    !                        H24 = YCL2ND_WFF(JJ) * H24F + YCL2ND_WFB(JJ) * H24B

    !                        DPH_io(IC, JC, KC) = DPH_io(IC, JC, KC) + H24
    !                    END DO
    !                END DO
    !            END DO
    !        END IF

    !CALL DEBUG_WRT_LOCAL(DPH_io, 1, N2DO(MYID), 'cony') ! test

    RETURN
END SUBROUTINE
