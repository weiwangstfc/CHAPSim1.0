!**********************************************************************************************************************************
!> @brief
!>        Calculate the velocity field and rhO * U field
!>        Eq.(A1c) in Mehdi paper or Eq.(4.61) in Mehdi thesis.
!> @details
!> SUBROUTINE: VELOUPDT_tg (in MYID = all)
!> SUBROUTINE: MASSFLUX_UPDATE_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE VELOUPDT_tg(NS)

    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS

    INTEGER(4) :: IC, IM
    INTEGER(4) :: JC, JM, JJ
    INTEGER(4) :: KC, KM
    INTEGER(4) :: NII
    REAL(WP) :: DFX11, DFX22, DFX33
    REAL(WP) :: COE1


    !===================update u===================================
    COE1 = DT * TALP(NS) * DXI
    DO KC = 1, NCL3
        DO JC = 1, N2DO(MYID)
            DO IC = 1, NCL1_tg
                IM = IMV_tg(IC)
                DFX11 = DPH_tg(IC, JC, KC) - DPH_tg(IM, JC, KC)
                Q_tg(IC, JC, KC, 1) = Q_tg(IC, JC, KC, 1) - DFX11 * COE1
            END DO
        END DO
    END DO

    !===================update W===================================
    COE1 = DT * TALP(NS) * DZI
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = 1, N2DO(MYID)
            DO IC = 1, NCL1_tg
                DFX33 = DPH_tg(IC, JC, KC) - DPH_tg(IC, JC, KM)
                Q_tg(IC, JC, KC, 3) = Q_tg(IC, JC, KC, 3) - DFX33 * COE1
            END DO
        END DO
    END DO


    !===================update V ===================================
    NII = 1
    IF(MYID == 0) NII = 2

    DO JC = NII, N2DO(MYID) !====== main DOMAIN ================
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        COE1 = DT * TALP(NS) * DYCI(JJ) / RNDI1(JJ)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_tg
                DFX22 = DPH_tg(IC, JC, KC) - DPH_tg(IC, JM, KC)
                Q_tg(IC, JC, KC, 2) = Q_tg(IC, JC, KC, 2) - DFX22 * COE1
            END DO
        END DO
    END DO

    !        IF (MYID == 0) THEN          !======BOTTOM ========== Repeat in the interfacE ======
    !            IF (iCase == iPIPEC) THEN
    !                DO KC = 1, NCL3
    !                    KS = KSYM(KC)
    !                    DO IC = 1, NCL1_tg
    !                        Q_tg(IC, 1, KC, 2) = 0.0_WP
    !                        Q_tg(IC, 0, KC, 2) = Q_tg(IC, 2, KS, 2)
    !                    END DO
    !                END DO
    !            ELSE
    !                DO KC = 1, NCL3
    !                    DO IC = 1, NCL1_tg
    !                        Q_tg(IC, 1, KC, 2) = 0.0_WP
    !                        Q_tg(IC, 0, KC, 2) = 0.0_WP
    !                    END DO
    !                END DO
    !            END IF
    !        ENDIF

    !        IF (MYID == NPSLV) THEN      !====== TOP ================
    !            DO KC = 1, NCL3
    !                DO IC = 1, NCL1_tg
    !                    Q_tg(IC, N2DO(MYID) + 1, KC, 2) = 0.0_WP
    !                END DO
    !            END DO
    !        ENDIF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE MASSFLUX_UPDATE_io(NS)
    USE thermal_info
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS

    INTEGER(4) :: IC, IM
    INTEGER(4) :: JC, JM, JJ
    INTEGER(4) :: KC, KM
    INTEGER(4) :: NXI, NYI, IDR, NXE
    REAL(WP) :: DFX11, DFX22, DFX33
    REAL(WP) :: COE1



    !===================================== X =============================
    IDR = 1
    NXI = 1
    NXE = NCL1_io
    IF(TgFlowFlg) THEN
        NXI = 2
        NXE = NCL1_io + 1
    END IF

    NYI = 1
    COE1 = DT * TALP(NS) * DXI
    DO IC = NXI, NXE
        IM = IMV_io(IC)
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                DFX11 = DPH_io(IC, JC, KC) - DPH_io(IM, JC, KC)    !d\phi / DX
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) - DFX11 * COE1
            END DO
        END DO
    END DO

    !===================================== Z =============================
    IDR = 3
    NXI = 1
    NYI = 1
    COE1 = DT * TALP(NS) * DZI
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = NYI, N2DO(MYID)
            DO IC = NXI, NCL1_io
                DFX33 = DPH_io(IC, JC, KC) - DPH_io(IC, JC, KM)
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) - DFX33 * COE1
            END DO
        END DO
    END DO

    !=====================================Y =============================
    IDR = 2
    NXI = 1
    NYI = 1
    IF(MYID == 0 .AND. iCase /= IBox3P) NYI = 2
    DO JC = NYI, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        COE1 = DT * TALP(NS) * DYCI(JJ) / RNDI1(JJ)
        DO IC = NXI, NCL1_io
            DO KC = 1, NCL3
                DFX22 = DPH_io(IC, JC, KC) - DPH_io(IC, JM, KC)
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) - DFX22 * COE1
            END DO
        END DO
    END DO


    CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
    CALL BC_WALL_G_io

    RETURN
END SUBROUTINE
