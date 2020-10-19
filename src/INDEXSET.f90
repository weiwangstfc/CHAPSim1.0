!**********************************************************************************************************************************
!> @brief
!>       Set up the plus and minus index for I, J, K
!> @details
!> SUBROUTINE: INDEXSET (in MYID = all)
!>             Set up the plus and minus index for I, J, K
!> @note
!> allow IPV_tg(NCL1) = NCL1 + 1, outlet BC
!>       IMV_tg(1) = 1- 1 = 0,  inlet B.C.
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE INDEXSET
    USE mesh_info
    USE cparam
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IC
    INTEGER(4) :: JC
    INTEGER(4) :: K, KC

    !== I ===
    IF(TgFlowFlg) THEN
        DO IC = 1, NCL1_tg
            IPV_tg(IC) = IC + 1
            IMV_tg(IC) = IC- 1
            IF(IC == NCL1_tg) IPV_tg(IC) = 1
            IF(IC == 1)       IMV_tg(IC) = NCL1_tg
        END DO
    END IF

    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        DO IC = 0, NCL1_io
            IPV_io(IC) = IC + 1
        END DO

        DO IC = 1, NCL1_io + 1
            IMV_io(IC) = IC- 1
        END DO
    ELSE
        IF(IoFlowFlg) THEN
            DO IC = 1, NCL1_io
                IPV_io(IC) = IC + 1
                IMV_io(IC) = IC- 1
                IF(IC == NCL1_io) IPV_io(IC) = 1
                IF(IC == 1)       IMV_io(IC) = NCL1_io
            END DO

        END IF
    END IF

    !== K ===
    !>      !===== KPV:    2, 3, 4, ...., NCL3,   1 =============
    !===== KC :    1, 2, 3, ...., NCL3- 1, NCL3 =============
    !===== KMV: NCL3, 1, 2, ...., NCL3-2, NCL3- 1 ============
    DO KC = 1, NCL3
        KPV(KC) = KC + 1
        KMV(KC) = KC - 1
        IF(KC == NCL3) KPV(KC) = 1
        IF(KC == 1)    KMV(KC) = NCL3
    END DO

    !>      @note: For symmetric boundary Condition, only in 3 direction.
    !>           Move forward half size of total cells.
    DO K = 1, NCL3                                  !@
        KSYM(k) = k + NCL3 / 2                          !@
        IF(KSYM(k) >  NCL3) KSYM(k) = KSYM(k) - NCL3  !@
    END DO
    !@
    !>      @note: for pipe radial direction in J /Y / R direction
    !>           Repeat the ending/Starting cell rather than mapping to the Starting/ending ones.
    !=====JPV:    2, 3, 4, ...., NCL2,   NCL2 =============
    !=====JC :    1, 2, 3, ...., NCL2 - 1, NCL2 =============
    !=====JMV:    1, 1, 2, ...., NCL2 -2, NCL2 - 1 ============
    !         DO JC = 1, NCL2
    !             JMv(JC) = JC - 1
    !             JPv(JC) = JC + 1
    !             IF(JC == 1)    JMv(JC) = JC
    !             IF(JC == NCL2) JPv(JC) = JC
    !         END DO

    !== J GLOBAL ===
    DO JC = 1, NCL2
        JGMV(JC) = JC - 1 ! 0~NCL2 - 1
        JGPV(JC) = JC + 1 ! 2~NND2
    END DO
    IF(iCase == IBox3P) THEN
        JGMV(1) = NCL2
        JGPV(NCL2) = 1
    ELSE
        JGMV(NND2) = NCL2
        JGPV(0) = 1
    END IF

    !== J LOCAL ==
    DO JC = 1, N2DO(MYID)
        JLMV(JC) = JC - 1 ! 0~N2DO(MYID) - 1
        JLPV(JC) = JC + 1 ! 2~N2DO(MYID) + 1
    END DO
    JLMV(N2DO(MYID) + 1) = N2DO(MYID)
    JLPV(0) = 1


    IF(MYID == 0) THEN

        OPEN(10, FILE = TRIM(FilePath0) // 'CHK_INDEX_X_Y_Z.dat')

        IF(IoFlowFlg) THEN
            WRITE(10, '(A)') 'IC, IPV_io, IMV_io'
            DO IC = 1, NCL1_io
                WRITE(10, '(3I5.1)') IC, IPV_io(IC), IMV_io(IC)
            END DO
        END IF

        IF(TgFlowFlg) THEN
            WRITE(10, '(A)') 'IC, IPV_tg, IMV_tg'
            DO IC = 1, NCL1_tg
                WRITE(10, '(3I5.1)') IC, IPV_tg(IC), IMV_tg(IC)
            END DO
        END IF


        WRITE(10, '(A)') 'KC, KPV, KMV'
        DO KC = 1, NCL3
            WRITE(10, '(3I5.1)') KC, KPV(KC), KMV(KC)
        END DO

        WRITE(10, '(A)') 'JC, JGPV, JGMV'
        DO JC = 1, NCL2
            WRITE(10, '(3I5.1)') JC, JGPV(JC), JGMV(JC)
        END DO

        WRITE(10, '(A)') 'JC, JLPV, JLMV, JJ, JGPV, JGMV'
        DO JC = 1, N2DO(MYID)
            WRITE(10, '(6I5.1)') JC, JLPV(JC), JLMV(JC), JCL2G(JC), JGPV(JCL2G(JC)), JGMV(JCL2G(JC))
        END DO


        CLOSE(10)
    END IF


    RETURN

END SUBROUTINE
