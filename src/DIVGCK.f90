!**********************************************************************************************************************************
!> @brief
!>        To calculate divergence of velocity and check continuty
!> @details
!> SUBROUTINE: DIVGCK_tg (in MYID = all)
!> SUBROUTINE: DIVGCK_io(in MYID = all)
!> SUBROUTINE: DIVGCK_Comm_io(in MYID = all)
!> SUBROUTINE: DIVGCK_Q_Comm_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010- Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013- Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE DIVGCK_tg
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    REAL(WP) :: DIVX
    REAL(WP) :: QMA
    REAL(WP) :: QMAX_WORK
    REAL(WP) :: DQCAP
    INTEGER(4) :: K, KP
    INTEGER(4) :: J, JP, JJ
    INTEGER(4) :: I, IP

    DIVX    = 0.0_WP
    QMA     = 0.0_WP
    QMAX_WORK = 0.0_WP
    DQCAP   = 0.0_WP
    DO K = 1, NCL3
        KP = KPV(K)
        DO J = 1, N2DO(MYID)
            JP = JLPV(J)
            JJ = JCL2G(J)
            DO I = 1, NCL1_tg
                IP = IPV_tg(I)
                DQCAP = (Q_tg(IP, J, K, 1) - Q_tg(I, J, K, 1)) * DXI       &
              + (Q_tg(I, JP, K, 2) - Q_tg(I, J, K, 2)) * DYFI(JJ) * RCCI1(JJ)     &
              + (Q_tg(I, J, KP, 3) - Q_tg(I, J, K, 3)) * DZI * RCCI2(JJ)
                DIVX = DMAX1(DABS(DQCAP), DIVX)
            END DO
        END DO
    END DO

    QMA = DIVX

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(QMA, QMAX_WORK, 1, MPI_DOUBLE_PRECISION,  &
    MPI_MAX, ICOMM, IERROR)

    MAXDIVGV_tg(1) = MAXDIVGV_tg(2)
    MAXDIVGV_tg(2) = QMAX_WORK

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DIVGCK_io
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info !test
    IMPLICIT NONE

    !REAL(WP) :: DIVX1, DIVX2, DIVX3
    REAL(WP) :: QMA1, QMA2, QMA3
    REAL(WP) :: QMAX_WORK1, QMAX_WORK2, QMAX_WORK3

    QMA1       = 0.0_WP
    QMAX_WORK1 = 0.0_WP

    QMA2       = 0.0_WP
    QMAX_WORK2 = 0.0_WP

    QMA3       = 0.0_WP
    QMAX_WORK3 = 0.0_WP

    IF(TgFlowFlg) THEN
        CALL DIVGCK_Comm_io(2, NCL1_io - 1,    QMA1)
        CALL DIVGCK_Comm_io(0, 1,              QMA2)
        CALL DIVGCK_Comm_io(NCL1_io, NCL1_io,  QMA3)
    ELSE
        CALL DIVGCK_Comm_io(1, NCL1_io,    QMA1)
        QMA2 = QMA1
        QMA3 = QMA1
    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(QMA1, QMAX_WORK1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(QMA2, QMAX_WORK2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(QMA3, QMAX_WORK3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

    MAXDIVGV_io(1) = QMAX_WORK1
    MAXDIVGV_io(2) = QMAX_WORK2
    MAXDIVGV_io(3) = QMAX_WORK3

    !IF(MYID == 0) WRITE(*, '(A, 3ES15.7)') '#Continuity Eq. = ', QMAX_WORK1, QMAX_WORK2, QMAX_WORK3

    ! !============= Test below===================
    ! IF(TgFlowFlg) THEN
    !     CALL DIVGCK_Q_Comm_io(2, NCL1_io - 1,    QMA1)
    !     CALL DIVGCK_Q_Comm_io(0, 1,              QMA2)
    !     CALL DIVGCK_Q_Comm_io(NCL1_io, NCL1_io,  QMA3)
    ! ELSE
    !     CALL DIVGCK_Q_Comm_io(1, NCL1_io,    QMA1)
    !     QMA2 = QMA1
    !     QMA3 = QMA1
    ! END IF
    !
    ! CALL MPI_BARRIER(ICOMM, IERROR)
    ! CALL MPI_ALLREDUCE(QMA1, QMAX_WORK1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    ! CALL MPI_ALLREDUCE(QMA2, QMAX_WORK2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    ! CALL MPI_ALLREDUCE(QMA3, QMAX_WORK3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    !
    ! !IF(MYID == 0) WRITE(*, '(A, 3ES15.7)') '#Divergence of V = ', QMAX_WORK1, QMAX_WORK2, QMAX_WORK3
    ! !============= Test abovE ==================

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE DIVGCK_Comm_io(IS, IE, DIVMAX)
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info !test
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IS
    INTEGER(4), INTENT(IN) :: IE
    REAL(WP), INTENT(OUT) :: DIVMAX

    INTEGER(4) :: KC, KP
    INTEGER(4) :: JC, JP, JJ
    INTEGER(4) :: IC, IP
    REAL(WP) :: DIVX, DIVY, DIVZ
    REAL(WP) :: DTRHO, DQCAP
    REAL(WP) :: COE2, COE3, COE4

    DIVMAX   = 0.0_WP
    !        IF(MYID == 0) WRITE(*, *) 'check:'

    COE4 = 1.0_WP / DT
    DO JC = 1, N2DO(MYID)
        JP  = JLPV(JC)
        JJ  = JCL2G(JC)
        COE2 = DYFI(JJ) * RCCI1(JJ)
        COE3 = DZI * RCCI2(JJ)
        DO KC = 1, NCL3
            KP = KPV(KC)
            !================== MAX. DIV IN THE main DOMAIN =================
            DO IC = IS, IE
                IP = IPV_io(IC)
                !========== D(\rho u) / DX at (i, J, K) ===========================
                DIVX  = ( G_io(IP, JC, KC, 1) - G_io(IC, JC, KC, 1) ) * DXI
                !========== D(\rho v) / Dy at (i, J, K) ===========================
                DIVY  = ( G_io(IC, JP, KC, 2) - G_io(IC, JC, KC, 2) ) * COE2
                !========== D(\rho w) / Dz at (i, J, K) ===========================
                DIVZ  = ( G_io(IC, JC, KP, 3) - G_io(IC, JC, KC, 3) ) * COE3
                !========== D \RHO / DT ========================================
                !DTRHO = ( DENSITY(IC, JC, KC) - DENSITYP(IC, JC, KC) ) * COE4
                DTRHO = DrhoDtP(IC, JC, KC)
                !=========== TOTAL ============================================
                DQCAP = DIVX + DIVY + DIVZ + DTRHO
                DIVMAX = DMAX1(DABS(DQCAP), DIVMAX)

                !                    IF(DABS(DQCAP) >  1.0E-7_WP) &
                !                    !IF(JJ == 1 .OR. JJ == NCL2) &!
                !                    WRITE(*, '(A, 4I3.1,8ES13.5)') 'divgck', MYID, JC, KC, IC, &
                !                    DIVX / DXI, DIVY / COE2, DIVZ/ COE3, &
                !                    DIVX / DXI + DIVY / COE2 + DIVZ/ COE3, &
                !                    DTRHO / COE4, DQCAP, &
                !                    G_io(IC, JP, KC, 2), G_io(IC, JC, KC, 2)


            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DIVGCK_Q_Comm_io(IS, IE, DIVMAX)
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE init_info !test
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IS
    INTEGER(4), INTENT(IN) :: IE
    REAL(WP), INTENT(OUT) :: DIVMAX

    INTEGER(4) :: KC, KP
    INTEGER(4) :: JC, JP, JJ
    INTEGER(4) :: IC, IP
    REAL(WP) :: DIVX, DIVY, DIVZ
    REAL(WP) :: DQCAP
    REAL(WP) :: COE2, COE3, COE4

    DIVMAX   = 0.0_WP
    !        IF(MYID == 0) WRITE(*, *) 'check:'

    COE4 = 1.0_WP / DT
    DO JC = 1, N2DO(MYID)
        JP  = JLPV(JC)
        JJ  = JCL2G(JC)
        COE2 = DYFI(JJ) * RCCI1(JJ)
        COE3 = DZI * RCCI2(JJ)
        DO KC = 1, NCL3
            KP = KPV(KC)
            !================== MAX. DIV IN THE main DOMAIN =================
            DO IC = IS, IE
                IP = IPV_io(IC)
                !========== D(\rho u) / DX at (i, J, K) ===========================
                DIVX  = ( Q_io(IP, JC, KC, 1) - Q_io(IC, JC, KC, 1) ) * DXI
                !========== D(\rho v) / Dy at (i, J, K) ===========================
                DIVY  = ( Q_io(IC, JP, KC, 2) - Q_io(IC, JC, KC, 2) ) * COE2
                !========== D(\rho w) / Dz at (i, J, K) ===========================
                DIVZ  = ( Q_io(IC, JC, KP, 3) - Q_io(IC, JC, KC, 3) ) * COE3

                DQCAP = DIVX + DIVY + DIVZ
                DIVMAX = DMAX1(DABS(DQCAP), DIVMAX)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
