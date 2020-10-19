!**********************************************************************************************************************************
!> @brief
!>        To calculate the mean velocity in the Streamwise
!> @details
!> SUBROUTINE: INITIAL_MEAN_TG (in MYID = all)
!> SUBROUTINE: INITIAL_MEAN_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE INITIAL_MEAN_TG
    USE init_info
    USE mesh_info
    USE flow_info
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    INTEGER(4) :: IDR, I, J, K, IP, KP, JJ
    REAL(WP) :: UUJJ, VELENTER

    UU(:, :, 1) = 0.0_WP
    DO IDR = 1, NDV

        IF (IDR == 1) THEN

            DO J = 1, N2DO(MYID)
                JJ = JCL2G(J)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP

                DO I = 1, NCL1_TG
                    IP = IPV_TG(I)
                    DO K = 1, NCL3
                        VELENTER = ( Q_tg(IP, J, K, IDR) + Q_tg(I, J, K, IDR) ) * 0.5_WP
                        UUJJ = VELENTER + UUJJ
                    END DO

                END DO
                UU(J, IDR, 1) = UUJJ / DBLE(NCL1_Tg * NCL3)
            END DO

        ELSE IF (IDR == 3) THEN

            DO J = 1, N2DO(MYID)
                JJ = JCL2G(J)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP

                DO I = 1, NCL1_TG

                    DO K = 1, NCL3
                        KP = KPV(K)
                        !VELENTER = ( Q_tg(I, J, KP, IDR) * RCCI1(JJ) + Q_tg(I, J, K, IDR) * RCCI1(JJ) ) * 0.5_WP
                        VELENTER = ( Q_tg(I, J, KP, IDR) + Q_tg(I, J, K, IDR) ) * 0.5_WP
                        UUJJ = VELENTER + UUJJ
                    END DO

                END DO
                UU(J, IDR, 1) = UUJJ / DBLE(NCL1_Tg * NCL3)
            END DO

        ELSE IF (IDR == 2) THEN
            !                NYI = 1
            !                IF(MYID == 0 .AND. iCase == iPIPEC) THEN
            !                    J = 1
            !                    JJ = JCL2G(J)
            !                    JP = JLPV(J)
            !                    JJP = JGPV(JJ)
            !                    UUJJ  = 0.0_WP
            !                    VELENTER = 0.0_WP
            !                    DO I = 1, NCL1_TG
            !                        DO K = 1, NCL3
            !                            KS = KSYM(K)
            !                            !V1 = (Q_tg(I, JP, K, IDR) - Q_tg(I, JP, KS, IDR)) * RNDI1(JJP) / 2.0_WP
            !                            !VELENTER = ( Q_tg(I, JP, K, IDR) * RNDI1(JJP) + V1 ) * 0.5_WP

            !                            UUJJ = VELENTER + UUJJ
            !                        END DO
            !                    END DO
            !                    UU(J, IDR, 1) = UUJJ / DBLE(NCL1_Tg * NCL3)
            !                    NYI = 2
            !                END IF

            !                DO J = NYI, N2DO(MYID)
            !                    JJ = JCL2G(J)
            !                    JP = JLPV(J)
            !                    UUJJ  = 0.0_WP
            !                    VELENTER = 0.0_WP
            !                    DO I = 1, NCL1_TG
            !                        DO K = 1, NCL3
            !                            VELENTER = ( Q_tg(I, JP, K, IDR) * RNDI1(JJ) + Q_tg(I, J, K, IDR) * RNDI1(JJ) ) * 0.5_WP
            !                            UUJJ = VELENTER + UUJJ
            !                        END DO
            !                    END DO
            !                    UU(J, IDR, 1) = UUJJ / DBLE(NCL1_Tg * NCL3)
            !                END DO

            DO J = 1, N2DO(MYID)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP
                DO I = 1, NCL1_TG
                    DO K = 1, NCL3
                        VELENTER = Q_tg(I, J, K, IDR)
                        UUJJ = VELENTER + UUJJ
                    END DO
                END DO
                UU(J, IDR, 1) = UUJJ / DBLE(NCL1_Tg * NCL3)
            END DO

        ELSE
            CALL ERRHDL('directionS SHOULD BE 1, 2 or 3, and no others.', MYID)
        END IF
    END DO


    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INITIAL_MEAN_io
    USE init_info
    USE mesh_info
    USE flow_info
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    INTEGER(4) :: IDR, I, J, K, IP, KP, JJ
    REAL(WP) :: UUJJ, VELENTER

    UU(:, :, 2) = 0.0_WP
    DO IDR = 1, NDV

        IF (IDR == 1) THEN

            DO J = 1, N2DO(MYID)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP

                DO I = 1, NCL1_io
                    IP = IPV_io(I)
                    DO K = 1, NCL3
                        VELENTER = ( Q_io(IP, J, K, IDR) + Q_io(I, J, K, IDR) ) * 0.5_WP
                        UUJJ = VELENTER + UUJJ
                    END DO

                END DO
                UU(J, IDR, 2) = UUJJ / DBLE(NCL1_io * NCL3)
            END DO

        ELSE IF (IDR == 3) THEN

            DO J = 1, N2DO(MYID)
                JJ = JCL2G(J)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP
                DO I = 1, NCL1_io
                    DO K = 1, NCL3
                        KP = KPV(K)
                        !VELENTER = ( Q_io(I, J, KP, IDR) * RCCI1(JJ) + Q_io(I, J, K, IDR) * RCCI1(JJ) ) * 0.5_WP
                        VELENTER = ( Q_io(I, J, KP, IDR) + Q_io(I, J, K, IDR) ) * 0.5_WP
                        UUJJ = VELENTER + UUJJ
                    END DO

                END DO
                UU(J, IDR, 2) = UUJJ / DBLE(NCL1_io * NCL3)
            END DO

        ELSE IF (IDR == 2) THEN

            !                NYI = 1

            !                IF(MYID == 0 .AND. iCase == iPIPEC) THEN
            !                    J = 1
            !                    JJ = JCL2G(J)
            !                    JP = JLPV(J)
            !                    JJP = JGPV(JJ)
            !                    UUJJ  = 0.0_WP
            !                    VELENTER = 0.0_WP
            !                    DO I = 1, NCL1_io
            !                        DO K = 1, NCL3
            !                            KS = KSYM(K)
            !                            V1 = (Q_io(I, JP, K, IDR) - Q_io(I, JP, KS, IDR)) * RNDI1(JJP) / 2.0_WP
            !                            VELENTER = ( Q_io(I, JP, K, IDR) * RNDI1(JJP) + V1 ) * 0.5_WP
            !                            UUJJ = VELENTER + UUJJ
            !                        END DO
            !                    END DO
            !                    UU(J, IDR, 2) = UUJJ / DBLE(NCL1_io * NCL3)
            !                    NYI = 2
            !                END IF

            !                DO J = NYI, N2DO(MYID)
            !                    JJ = JCL2G(J)
            !                    JP =JLPV(J)
            !                    UUJJ  = 0.0_WP
            !                    VELENTER = 0.0_WP
            !                    DO I = 1, NCL1_io
            !                        DO K = 1, NCL3
            !                            VELENTER = ( Q_io(I, JP, K, IDR) * RNDI1(JJ) + Q_io(I, J, K, IDR) * RNDI1(JJ) ) * 0.5_WP
            !                            UUJJ = VELENTER + UUJJ
            !                        END DO
            !                    END DO
            !                    UU(J, IDR, 2) = UUJJ / DBLE(NCL1_io * NCL3)
            !                END DO

            DO J = 1, N2DO(MYID)
                UUJJ  = 0.0_WP
                VELENTER = 0.0_WP
                DO I = 1, NCL1_io
                    DO K = 1, NCL3
                        VELENTER = Q_io(I, J, K, IDR)
                        UUJJ = VELENTER + UUJJ
                    END DO
                END DO
                UU(J, IDR, 2) = UUJJ / DBLE(NCL1_io * NCL3)
            END DO



        ELSE
            CALL ERRHDL('directionS SHOULD BE 1, 2 or 3, and no others.', MYID)
        END IF
    END DO


    RETURN

END SUBROUTINE
