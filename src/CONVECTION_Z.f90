!**********************************************************************************************************************************
!> @brief
!>        Calculate the convection term in Z direction.
!> @details
!> SUBROUTINE: CONVECTION_Y_tg (in MYID = all)
!> SUBROUTINE: CONVECTION_Y_io(in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 04/ 2014- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CONVECTION_Z_tg
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
    INTEGER(4) :: KC, KM, KP, KS, KSM
    REAL(WP) :: COE1, COE2, COE3
    REAL(WP) :: H31, H32, H33
    REAL(WP) :: q2s1, q2s2
    REAL(WP) :: h32n
    REAL(WP) :: q2e, q2w
    REAL(WP) :: d11q1e


    RHSLLPHI_tg(:, :, :) = 0.0_WP ! to store convection term in z direction.
    COE1 = XND2CL * ZND2CL * DXI
    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJM = JGMV(JJ)
        JJP = JGPV(JJ)
        COE2 = DYFI(JJ) * ZND2CL
        COE3 = ZND2CL * ZND2CL * DZI * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)
            DO IC = 1, NCL1_tg
                IP = IPV_tg(IC)
                IM = IMV_tg(IC)

                H31 = ((Q_tg(IP, JC, KC, 1) + Q_tg(IP, JC, KM, 1)) *           &
                (Q_tg(IP, JC, KC, 3) + Q_tg(IC, JC, KC, 3)) -                  &
                (Q_tg(IC, JC, KC, 1) + Q_tg(IC, JC, KM, 1)) *                &
                (Q_tg(IC, JC, KC, 3) + Q_tg(IM, JC, KC, 3)) ) * COE1

                H32 = ( ( Q_tg(IC, JP, KC, 2) + Q_tg(IC, JP, KM, 2) ) *         &
                ( YCL2ND_WFF(JJP) * Q_tg(IC, JP, KC, 3) * RCCI1(JJP) +       &
                YCL2ND_WFB(JJP) * Q_tg(IC, JC, KC, 3) * RCCI1(JJ)   ) -      &
                ( Q_tg(IC, JC, KC, 2) + Q_tg(IC, JC, KM, 2) ) *               &
                ( YCL2ND_WFF(JJ) * Q_tg(IC, JC, KC, 3) * RCCI1(JJ) +        &
                YCL2ND_WFB(JJ) * Q_tg(IC, JM, KC, 3) * RCCI1(JJM) ) ) * COE2

                H33 = ((Q_tg(IC, JC, KP, 3) + Q_tg(IC, JC, KC, 3)) *            &
                (Q_tg(IC, JC, KP, 3) + Q_tg(IC, JC, KC, 3)) -                   &
                (Q_tg(IC, JC, KC, 3) + Q_tg(IC, JC, KM, 3)) *                 &
                (Q_tg(IC, JC, KC, 3) + Q_tg(IC, JC, KM, 3)) ) * COE3

                RHSLLPHI_tg(IC, JC, KC) = -(H31 + H32 + H33)
            END DO
        END DO
    END DO

    !>       @note below IS only for pipe/ CoriolIS force!!!!!!!!!!!!!!!      for athmut
    IF (iCase == IPIPEC .OR. iCase == IANNUL) THEN
        DO JC = 1, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)      !@
            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)

                DO IC = 1, NCL1_tg
                    IF (JJ == 1) THEN
                        KS  = KSYM(KC)
                        KSM = KSYM(KM)
                        q2s1 = (Q_tg(IC, JP, KC, 2) - Q_tg(IC, JP, KS, 2)) * 0.50_WP * RNDI1(JJP)
                        q2s2 = (Q_tg(IC, JP, KM, 2) - Q_tg(IC, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP)
                        q2e = Q_tg(IC, JP, KC, 2) * RNDI1(JJP) + q2s1
                        q2w = Q_tg(IC, JP, KM, 2) * RNDI1(JJP) + q2s2
                    ELSE
                        q2e = Q_tg(IC, JP, KC, 2) * RNDI1(JJP) + Q_tg(IC, JC, KC, 2) * RNDI1(JJ)
                        q2w = Q_tg(IC, JP, KM, 2) * RNDI1(JJP) + Q_tg(IC, JC, KM, 2) * RNDI1(JJ)
                    END IF
                    h32n = Q_tg(IC, JC, KC, 3) * (q2E + Q2w) * 0.250_WP * RCCI1(JJ)
                    d11q1e = (q2E -q2w) * DZI * RCCI1(JJ)
                    RHSLLPHI_tg(IC, JC, KC) = - H32n + d11q1e / REn + RHSLLPHI_tg(IC, JC, KC)
                END DO
            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE

!********************************************************************************************
SUBROUTINE CONVECTION_Z_io
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
    INTEGER(4) :: KC, KM, KP, KS, KSM!, NYI

    REAL(WP) :: H31, H32, H33
    REAL(WP) :: H31F, H31B
    REAL(WP) :: H32F, H32B
    REAL(WP) :: H33F, H33B
    !REAL(WP) :: H34, H34F, H34B, H34F1, H34F2, H34B1, H34B2
    REAL(WP) :: COE1, COE2, COE30, COE31
    !REAL(WP) :: Ur1z1, Ur1z2
    REAL(WP) :: q2s1, q2s2, h32n, q2e, q2w


    RHSLLPHI_io(:, :, :) = 0.0_WP ! to store convection term in z direction.

    COE1  = DXI * XND2CL * ZND2CL
    COE30 = DZI * ZND2CL * ZND2CL
    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)
        JJM = JGMV(JJ)
        COE2 = DYFI(JJ) * ZND2CL
        COE31 =COE30 * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)
            DO IC = 1, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)

                ! \frac{\pARtial {\rho u w}}{\pARtial x}_{i, J, K'}
                !{i'+ 1, J, K'}
                H31F = ( G_io(IP, JC, KC, 1) + G_io(IP, JC, KM, 1) ) * &
                ( Q_io(IP, JC, KC, 3) + Q_io(IC, JC, KC, 3) )
                !{i', J, K'}
                H31B = ( G_io(IC, JC, KC, 1) + G_io(IC, JC, KM, 1) ) * &
                ( Q_io(IC, JC, KC, 3) + Q_io(IM, JC, KC, 3) )
                !{i, J, K'}
                H31  = ( H31F - H31B) * COE1

                ! \frac{\pARtial {\rho v w}}{\pARtial y}_{i, J, K'}
                !{i, J'+ 1, K'}
                H32F = ( G_io(IC, JP, KC, 2) + G_io(IC, JP, KM, 2) ) * &
                ( YCL2ND_WFF(JJP) * Q_io(IC, JP, KC, 3) * RCCI1(JJP) +   &
                YCL2ND_WFB(JJP) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) )
                !{i, J', K'}
                H32B = ( G_io(IC, JC, KC, 2) + G_io(IC, JC, KM, 2) ) * &
                ( YCL2ND_WFF(JJ) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) +   &
                YCL2ND_WFB(JJ) * Q_io(IC, JM, KC, 3) * RCCI1(JJM) )
                !{i, J, K'}
                H32 = ( H32F - H32B) * COE2

                ! \frac{\pARtial {\rho w w}}{\pARtial z}_{i, J, K'}
                !(I, J, K)
                H33F = ( G_io(IC, JC, KP, 3) + G_io(IC, JC, KC, 3) ) *  &
                ( Q_io(IC, JC, KP, 3) + Q_io(IC, JC, KC, 3) )
                !(I, J, K - 1)
                H33B = ( G_io(IC, JC, KC, 3) + G_io(IC, JC, KM, 3) ) *  &
                ( Q_io(IC, JC, KC, 3) + Q_io(IC, JC, KM, 3) )
                !(I, J, K')
                H33 = ( H33F - H33B) * COE31

                RHSLLPHI_io(IC, JC, KC) = -(H31 + H32 + H33)     !H for momentum Z direction.

                !IF(JJ == 1 .AND. IC == 1 .AND. KC == 1) WRITE(*, '(A, 4ES13.5)') 'convz',H31,H32,H33,RHSLLPHI_io(IC, JC, KC)
            END DO

        END DO

    END DO


    !below IS only for pipe/ CoriolIS force!!!!!!!!!!!!!!!
    IF (iCase == IPIPEC .OR. iCase == IANNUL ) THEN
        DO JC = 1, N2DO(MYID)
            JM = JLMV(JC)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            JJM = JGMV(JJ)
            JJP = JGPV(JJ)      !@
            DO KC = 1, NCL3
                KM = KMV(KC)
                KP = KPV(KC)

                DO IC = 1, NCL1_io
                    IF (JJ == 1) THEN
                        KS  = KSYM(KC)
                        KSM = KSYM(KM)
                        q2s1 = (G_io(IC, JP, KC, 2) - G_io(IC, JP, KS, 2)) * 0.50_WP * RNDI1(JJP)
                        q2s2 = (G_io(IC, JP, KM, 2) - G_io(IC, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP)
                        q2e = G_io(IC, JP, KC, 2) * RNDI1(JJP) + q2s1
                        q2w = G_io(IC, JP, KM, 2) * RNDI1(JJP) + q2s2
                    ELSE
                        q2e = G_io(IC, JP, KC, 2) * RNDI1(JJP) + G_io(IC, JC, KC, 2) * RNDI1(JJ)
                        q2w = G_io(IC, JP, KM, 2) * RNDI1(JJP) + G_io(IC, JC, KM, 2) * RNDI1(JJ)
                    END IF
                    h32N = Q_io(IC, JC, KC, 3) * (q2e + Q2w) * YND2CL * ZND2CL * RCCI1(JJ)
                    RHSLLPHI_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC) - H32n
                END DO
            END DO
        END DO

    END IF


    !        IF(iCase == iPIPEC .OR. iCase == iANNUL) THEN
    !            NYI = 1

    !            IF(MYID == 0 .AND. iCase == IPIPEC) THEN
    !                !( Q_io(IC, JC, KC, 2) * RNDI1(JJ) ==>(Q_io(IC, JP, KC, 2) - Q_io(IC, JP, KS, 2)) * RNDI1(JJP) * 0.5_WP
    !                JC = 1
    !                JJ = JCL2G(JC)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJP = JGPV(JJ) !JPV(JJ)
    !                JJM = JGMV(JJ) !JMV(JJ)
    !                COE1 = ZND2CL * RNDI1(JJP)

    !                DO KC = 1, NCL3
    !                    KP = KPV(KC)
    !                    KM = KMV(KC)
    !                    KS = KSYM(KC)
    !                    KSM = KSYM(KM)
    !                    DO IC = 1, NCL1_io
    !                        ! AT I, J'+ 1, K'
    !                        H34F1 = ( G_io(IC, JP, KC, 2) + G_io(IC, JP, KM, 2) ) * COE1     !Gr average over theta
    !                        H34F2 = ( YCL2ND_WFF(JJP) * Q_io(IC, JP, KC, 3) * RCCI1(JJP) + &
    !                                  YCL2ND_WFB(JJP) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) ) !G_theta average over r
    !                        H34F = H34F1 * H34F2

    !                        ! AT I, J', K'
    !                        Ur1z1  = ( G_io(IC, JP, KC, 2) - G_io(IC, JP, KS, 2)) * COE1!* RNDI1(JJP) * 0.5_WP
    !                        Ur1z2  = ( G_io(IC, JP, KM, 2) - G_io(IC, JP, KSM, 2)) * COE1!* RNDI1(JJP) * 0.5_WP
    !                        H34B1 =  (Ur1z1 + Ur1z2) * ZND2CL     !Gr average over theta
    !                        H34B2 = ( YCL2ND_WFF(JJ) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) +   &
    !                                  YCL2ND_WFB(JJ) * Q_io(IC, JM, KC, 3) * RCCI1(JJM) ) !G_theta average over r
    !                        H34B = H34B1 * H34B2

    !                        ! AT I, J, K'
    !                        H34 = (H34F + H34B) * YND2CL

    !                        RHSLLPHI_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC) - H34
    !                    END DO
    !                END DO

    !                NYI = 2
    !            END IF

    !            DO JC = NYI, N2DO(MYID)
    !                JJ = JCL2G(JC)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJP = JGPV(JJ)
    !                JJM = JGMV(JJ)
    !                COE1 = ZND2CL * RNDI1(JJP)
    !                COE2 = ZND2CL * RNDI1(JJ)

    !                DO KC = 1, NCL3
    !                    KP = KPV(KC)
    !                    KM = KMV(KC)
    !                    DO IC = 1, NCL1_io
    !                        ! AT I, J'+ 1, K'
    !                        H34F1 = ( G_io(IC, JP, KC, 2) + G_io(IC, JP, KM, 2) ) * COE1    !Gr average over theta
    !                        H34F2 = ( YCL2ND_WFF(JJP) * Q_io(IC, JP, KC, 3) * RCCI1(JJP) + &
    !                                  YCL2ND_WFB(JJP) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) ) !G_theta average over r
    !                        H34F = H34F1 * H34F2

    !                        ! AT I, J', K'
    !                        H34B1 = ( G_io(IC, JC, KC, 2) + G_io(IC, JC, KM, 2) ) * COE2     !Gr average over theta
    !                        H34B2 = ( YCL2ND_WFF(JJ) * Q_io(IC, JC, KC, 3) * RCCI1(JJ) +   &
    !                                  YCL2ND_WFB(JJ) * Q_io(IC, JM, KC, 3) * RCCI1(JJM) ) !G_theta average over r
    !                        H34B = H34B1 * H34B2

    !                        ! AT I, J, K'
    !                        H34 = (H34F + H34B) * YND2CL

    !                        RHSLLPHI_io(IC, JC, KC) = RHSLLPHI_io(IC, JC, KC) - H34
    !                    END DO
    !                END DO
    !            END DO


    !        END IF

    !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io, 1, N2DO(MYID), 'conz') ! test

    RETURN
END SUBROUTINE
