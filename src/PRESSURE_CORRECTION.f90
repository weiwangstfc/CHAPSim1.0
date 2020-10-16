!**********************************************************************************************************************************
!> @brief
!>        Calculate the pressure field
!>        Eq.(A1c) in Mehdi paper or Eq.(4.61) in Mehdi thesIS.
!> @details
!> SUBROUTINE: PRCALC_tg (in MYID = all)
!> SUBROUTINE: PRCALC_tg (in MYID = all)
!> SUBROUTINE: CHECK_FFT_SOLVER (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 04/ 2014- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PRCALC_tg(NS)
    USE init_info
    USE flow_info
    USE thermal_info
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    REAL(WP) :: BE, LLPHI !,  COE21, COE22, COE3
    INTEGER(4) :: IC, JC, KC !, KM, KP, IM, IP, JM, JP, JJ, JJP

    BE = 0.50_WP * TALP(NS) * DT * CVISC * M_inlet / D_inlet ! \alpha_1 * dt / 2 / Re in Eq.(A1d) of Mehdi thesIS.


    DO KC = 1, NCL3
        DO JC = 1, N2DO(MYID)
            DO IC = 1, NCL1_TG
                LLphI = RHSLLPHI_tg(IC, JC, KC)
                PR_TG(IC, JC, KC) = PR_TG(IC, JC, KC) + DPH_TG(IC, JC, KC) - BE * LLphi
            END DO
        END DO
    END DO


    !        IF(iCase == ICHANL) THEN
    !            DO KC = 1, NCL3
    !                KP = KPV(KC)
    !                KM = KMV(KC)
    !                DO JC = 1, N2DO(MYID)
    !                    JM = JLMV(JC)
    !                    JP = JLPV(JC)
    !                    JJ = JCL2G(JC)
    !                    DO IC = 1, NCL1_TG
    !                        IP = IPV_TG(IC)
    !                        IM = IMV_TG(IC)
    !                        LLphI = (       DPH_TG(IP, JC, KC) - &
    !                                 2.0_WP * DPH_TG(IC, JC, KC) + &
    !                                        DPH_TG(IM, JC, KC) ) * DXQI + &
    !                                ( DPH_TG(IC, JP, KC) * APPH(JJ) + &
    !                                  DPH_TG(IC, JC, KC) * ACPH(JJ) + &
    !                                  DPH_TG(IC, JM, KC) * AMPH(JJ) ) + &
    !                                (       DPH_TG(IC, JC, KP) - &
    !                                 2.0_WP * DPH_TG(IC, JC, KC) + &
    !                                        DPH_TG(IC, JC, KM) ) * DZQI
    !                        !or
    !                        !LLphI = RHSLLPHI_tg(IC, JC, KC)
    !                        PR_TG(IC, JC, KC) = PR_TG(IC, JC, KC) + DPH_TG(IC, JC, KC) - BE * LLphi
    !                    END DO
    !                END DO
    !            END DO
    !        END IF

    !        IF(iCase /= ICHANL) THEN
    !            DO JC = 1, N2DO(MYID)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJ = JCL2G(JC)
    !                JJP = JGPV(JJ)
    !                COE3 = DZQI * RCCI2(JJ)
    !                COE21 = 1.0_WP / RNDI1(JJP) * DYCI(JJP) * DYFI(JJ) * RCCI1(JJ)
    !                COE22 = 1.0_WP / RNDI1(JJ) * DYCI(JJ) * DYFI(JJ) * RCCI1(JJ)
    !                DO KC = 1, NCL3
    !                    KP = KPV(KC)
    !                    KM = KMV(KC)

    !                    DO IC = 1, NCL1_TG
    !                        IP = IPV_TG(IC)
    !                        IM = IMV_TG(IC)
    !                        LLphI = (         DPH_TG(IP, JC, KC) -  &
    !                                  2.0_WP * DPH_TG(IC, JC, KC) +  &
    !                                         DPH_TG(IM, JC, KC) ) * DXQI +  &
    !                               ( (DPH_TG(IC, JP, KC) - DPH_TG(IC, JC, KC) ) * COE21-      &
    !                                 (DPH_TG(IC, JC, KC) - DPH_TG(IC, JM, KC) ) * COE22  ) +  &
    !                               (         DPH_TG(IC, JC, KP) -                         &
    !                                  2.0_WP * DPH_TG(IC, JC, KC) +                         &
    !                                         DPH_TG(IC, JC, KM) ) * COE3
    !                        PR_TG(IC, JC, KC) = PR_TG(IC, JC, KC) + DPH_TG(IC, JC, KC) - BE * LLphi
    !                    END DO
    !                END DO
    !            END DO
    !        END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PRCALC_io(NS)
    USE init_info
    USE flow_info
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NS
    INTEGER(4) :: IC
    INTEGER(4) :: JC
    INTEGER(4) :: KC
    REAL(WP) :: BE,LLphi

    IF(iVisScheme == VisImplicit) THEN

        BE = 0.50_WP * TALP(NS) * DT * CVISC ! \alpha_1 * dt / 2 / Re in Eq.(A1d) of Mehdi thesIS.
        DO KC = 1, NCL3
            DO JC = 1, N2DO(MYID)
                DO IC = 1, NCL1_io
                    LLphI = RHSLLPHI_io(IC, JC, KC)
                    PR_io(IC, JC, KC) = PR_io(IC, JC, KC) + DPH_io(IC, JC, KC) - BE * LLphi
                END DO
            END DO
        END DO

    ELSE IF(iVisScheme == VisExplicit) THEN

        DO KC = 1, NCL3
            DO JC = 1, N2DO(MYID)
                DO IC = 1, NCL1_io
                    PR_io(IC, JC, KC) = PR_io(IC, JC, KC) + DPH_io(IC, JC, KC)
                END DO
            END DO
        END DO

    ELSE
    END IF

    !BELOW IS FOR IMPLICIT VISCOUS AND iCase = 2
    !        DO KC = 1, NCL3
    !            KP = KPV(KC)
    !            KM = KMV(KC)
    !            DO JC = 1, N2DO(MYID)
    !                JM = JLMV(JC)
    !                JP = JLPV(JC)
    !                JJ = JCL2G(JC)
    !                JJP = JGPV(JJ)
    !                COE2 = DYFI(JJ) * RCCI1(JJ)
    !                COE3 = DZQI * RCCI2(JJ)
    !                DO IC = 1, NCL1_io
    !                    IP = IPV(IC)
    !                    IM = IMV(IC)
    !                    LLphI = (         DPH(IP, JC, KC) -  &
    !                              2.0_WP * DPH(IC, JC, KC) +  &
    !                                     DPH(IM, JC, KC) ) * DXQI +  &
    !                           ( (DPH(IC, JP, KC) - DPH(IC, JC, KC) ) * Rc(JPp) * DYCI(JJP) -      &
    !                             (DPH(IC, JC, KC) - DPH(IC, JM, KC) ) * Rc(JJ) * DYCI(JJ)  ) * COE2 + &
    !                           (         DPH(IC, JC, KP) -                         &
    !                              2.0_WP * DPH(IC, JC, KC) +                         &
    !                                     DPH(IC, JC, KM) ) * COE3
    !                    PR(IC, JC, KC) = PR(IC, JC, KC) + DPH(IC, JC, KC) - BE * LLphi
    !                END DO
    !            END DO
    !        END DO

    CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
    CALL BC_WALL_PR_io

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CHECK_FFT_SOLVER
    USE init_info
    USE flow_info
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4) :: IC, JC, KC, KM, KP, IM, IP, JM, JP, JJ
    REAL(WP) :: LLphiX, LLphiY, LLphiZ, LLPHI, COE0


    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        DO KC = 1, NCL3
            KP = KPV(KC)
            KM = KMV(KC)
            DO IC = 1, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)
                LLphiX = (       DPH_io(IP, JC, KC) - &
                2.0_WP * DPH_io(IC, JC, KC) + &
                DPH_io(IM, JC, KC) ) * DXQI

                LLphiY = DPH_io(IC, JP, KC) * APPH(JJ) + &
                DPH_io(IC, JC, KC) * ACPH(JJ) + &
                DPH_io(IC, JM, KC) * AMPH(JJ)

                LLphiZ = (       DPH_io(IC, JC, KP) - &
                2.0_WP * DPH_io(IC, JC, KC) + &
                DPH_io(IC, JC, KM) ) * DZQI

                LLphI = LLphiX + LLphiY + LLphiZ
                COE0  = DABS( LLphI - RHSLLPhi_io(IC, JC, KC) )

                IF(COE0 >  1.0E-8_WP ) &
                WRITE(logflg_pg, '(A, 4I4.1, 3ES22.14)') '# FFT -comp', MYID, JC, KC, IC,LLphi, &
                RHSLLPhi_io(IC, JC, KC), COE0
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE



!SUBROUTINE GRAD_P_WALL_NONZERO
!        USE init_info
!        USE flow_info
!        USE mesh_info
!        USE thermal_info
!        IMPLICIT NONE

!        INTEGER(4) :: IC, JC, JJ, KC, N
!        REAL(WP) :: COE1, COE2, COE3


!        DPDYWAL(:, :, :) = 0.0_WP
!        IF(iThermoDynamics == 0) RETURN
!        IF(iThermalWallType /= BC_Fixed_Temperature )  RETURN
!        IF(MYID >  0 .AND. MYID < NPSLV) RETURN

!        IF(MYID == 0) THEN
!            JC = 1
!            N  = 1
!            JJ = 1
!            COE1 = 2.0_WP / 3.0_WP * DYFI(JJ) / DT/ REN

!            COE3 = AMPH0
!        END IF

!        IF(MYID == NPSLV) THEN
!            JC = N2DO(MYID)
!            N  = 2
!            JJ = NCL2
!            COE1 = -2.0_WP / 3.0_WP * DYFI(JJ) / DT/ REN

!            COE3 = APPH0
!        END IF



!        DO IC = 1, NCL1_io
!            COE2 = M_WAL_GV(IC, N) * D_WAL_GV(IC, N)
!            DO KC = 1, NCL3
!                DPDYWAL(IC, KC, N) = COE1 * COE2 / DENSITY(IC, JC, KC) / (2.0_WP * D_WAL_GV(IC, N) - DENSITY(IC, JC, KC)) * &
!                (DENSITY(IC, JC, KC) - DENSITY0(IC, JC, KC))
!                !WRITE(*, *) 'DPDY', MYID, IC,COE1,COE2, DPDYWAL(IC, KC, N), DPDYWAL(IC, KC, N) * COE3 / DYFI(JJ)
!            END DO

!        END DO


!        RETURN
!    END SUBROUTINE



!SUBROUTINE REVERSE_DPH_GHOST
!        USE init_info
!        USE flow_info
!        USE mesh_info
!        IMPLICIT NONE

!        REAL(WP) :: COE1, BE, LLPHI,   COE21, COE22, COE3
!        INTEGER(4) :: IC, JC, KC, KM, KP, IM, IP, JM, JP, JJ, JJP

!        IF(MYID == 0) THEN
!            DO JC = 1, N2DO(MYID)
!                JM = JLMV(JC)
!                JP = JLPV(JC)
!                JJ = JCL2G(JC)
!                DO KC = 1, NCL3
!                KP = KPV(KC)
!                KM = KMV(KC)
!                    DO IC = 1, NCL1_io
!                        IP = IPV_io(IC)
!                        IM = IMV_io(IC)
!                        LLphI = (       DPH_io(IP, JC, KC) - &
!                                 2.0_WP * DPH_io(IC, JC, KC) + &
!                                        DPH_io(IM, JC, KC) ) * DXQI + &
!                                ( DPH_io(IC, JP, KC) * APPH(JJ) + &
!                                  DPH_io(IC, JC, KC) * ACPH(JJ) + &
!                                  0.0_WP                   ) + &
!                                (       DPH_io(IC, JC, KP) - &
!                                 2.0_WP * DPH_io(IC, JC, KC) + &
!                                        DPH_io(IC, JC, KM) ) * DZQI
!                         DPH_io(IC, JM, KC) = (RHSLLPhi_io(IC, JC, KC) -LLphi) /AMPH(JJ) !INFINITY
!                         WRITE(*, *) MYID, JC, KC, IC, RHSLLPhi_io(IC, JC, KC), LLphi, (RHSLLPhi_io(IC, JC, KC) -LLphi), AMPH(JJ)
!                    END DO
!                END DO
!            END DO
!        END IF

!        IF(MYID == NPSLV) THEN
!            DO JC = 1, N2DO(MYID)
!                JM = JLMV(JC)
!                JP = JLPV(JC)
!                JJ = JCL2G(JC)
!                DO KC = 1, NCL3
!                KP = KPV(KC)
!                KM = KMV(KC)
!                    DO IC = 1, NCL1_io
!                        IP = IPV_io(IC)
!                        IM = IMV_io(IC)
!                        LLphI = (       DPH_io(IP, JC, KC) - &
!                                 2.0_WP * DPH_io(IC, JC, KC) + &
!                                        DPH_io(IM, JC, KC) ) * DXQI + &
!                                ( 0.0_WP + &
!                                  DPH_io(IC, JC, KC) * ACPH(JJ) + &
!                                  DPH_io(IC, JM, KC) * AMPH(JJ) ) + &
!                                (       DPH_io(IC, JC, KP) - &
!                                 2.0_WP * DPH_io(IC, JC, KC) + &
!                                        DPH_io(IC, JC, KM) ) * DZQI
!                        DPH_io(IC, JP, KC) = (RHSLLPhi_io(IC, JC, KC) -LLphi) /APPH(JJ) !INFINITY
!                        WRITE(*, *) MYID, JC, KC, IC, RHSLLPhi_io(IC, JC, KC), LLphi, (RHSLLPhi_io(IC, JC, KC) -LLphi), APPH(JJ)
!                    END DO
!                END DO
!            END DO
!        END IF


!        RETURN
!    END SUBROUTINE
