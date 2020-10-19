!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE: PP_MEAN_ZX_FLOW_TG
!> SUBROUTINE: PP_MEAN_ZX_FLOW_Xperiodic_io
!> SUBROUTINE: INST_Qio_Gradient_CellCentred
!> SUBROUTINE: PP_DRIVEN_FORCE
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_ZX_FLOW_TG
    !>  @NOTE :  I, J,  L
    !>           1, 1,  1
    !>           1, 2,  2
    !>           1, 3,  3
    !>           2, 2,  4
    !>           2, 3,  5
    !>           3, 3,  6
    !>  @NOTE :  I, J, K,  L
    !>           1, 1, 1,  1
    !>           1, 1, 2,  2
    !>           1, 1, 3,  3
    !>           1, 2, 2,  4
    !>           1, 2, 3,  5
    !>           1, 3, 3,  6
    !>           2, 2, 2,  7
    !>           2, 2, 3,  8
    !>           2, 3, 3,  9
    !>           3, 3, 3,  10

    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, KS, KSM, KSP
    INTEGER(4) :: M, N, H, L, P, L1, L2
    REAL(WP) :: U_CCT(NDV)
    REAL(WP) :: U_CCT_IP(NDV), U_CCT_JP(NDV), U_CCT_KP(NDV)
    REAL(WP) :: U_CCT_IM(NDV), U_CCT_JM(NDV), U_CCT_KM(NDV)
    REAL(WP) :: DVDL_CCT(NDV, NDV)
    REAL(WP) :: COE0, COE1, COE2, COE3
    REAL(WP) :: COG1(NDV, NDV), COG2(NDV * NDV, NDV * NDV)

    U_CCT = 0.0_WP
    U_CCT_IP = 0.0_WP
    U_CCT_JP = 0.0_WP
    U_CCT_KP = 0.0_WP
    U_CCT_IM = 0.0_WP
    U_CCT_JM = 0.0_WP
    U_CCT_KM = 0.0_WP

    COE0 = VL1313_tg
    COE1 = VL1313_tg * 0.5_WP
    COE2 = VL1313_tg * (0.5_WP**2)
    COE3 = VL1313_tg * (0.5_WP**3)


    COG1(1, 1) = VL1313_tg * DXI
    COG1(1, 2) = VL1313_tg * 0.25_WP
    COG1(1, 3) = VL1313_tg * 0.25_WP * DZI

    COG1(2, 1) = VL1313_tg * 0.25_WP * DXI
    COG1(2, 2) = VL1313_tg
    COG1(2, 3) = VL1313_tg * 0.25_WP * DZI

    COG1(3, 1) = VL1313_tg * 0.25_WP * DXI
    COG1(3, 2) = VL1313_tg * 0.25_WP
    COG1(3, 3) = VL1313_tg * DZI


    !        DO M = 1, NDV
    !            DO N = 1, NDV
    !                DO H = 1, NDV
    !                    DO P = 1, NDV
    !                        L1 = (M - 1) * 3 + H
    !                        L2 = (N - 1) * 3 + P
    !                        COG2(L1, L2) = COG1(M,H) * COG1(N,P)
    !                    END DO
    !                END DO
    !            END DO
    !        END DO
    DO M = 1, NDV
        DO H = 1, NDV
            L1 = (M - 1) * NDV + H
            DO N = 1, NDV
                DO P = 1, NDV
                    L2 = (N - 1) * NDV + P
                    COG2(L1, L2) = COG1(M,H) * COG1(N,P)
                END DO
            END DO
        END DO
    END DO


    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        JM = JLMV(J)
        JP = JLPV(J)
        JJP = JGPV(JJ)
        JJM = JGMV(JJ)

        U1xzL_tg(J, :) = 0.0_WP
        UPxzL_tg(J, :) = 0.0_WP

        U2xzL_tg(J, :) = 0.0_WP
        U3xzL_tg(J, :) = 0.0_WP

        DVDL1xzL_tg(J, :, :) = 0.0_WP
        DVDLPxzL_tg(J, :, :) = 0.0_WP
        DVDL2xzL_tg(J, :, :) = 0.0_WP

        DO I = 1, NCL1_tg
            IP = IPV_TG(I)
            IM = IMV_TG(I)
            DO K = 1, NCL3

                KP = KPV(K)
                KM = KMV(K)
                KS = KSYM(K)
                KSM = KSYM(KM)
                KSP = KSYM(KP)

                U_CCT(1) = Q_tg(I, J, K, 1)  + Q_tg(IP, J, K, 1)            ! U at i, j, k
                U_CCT_JP(1) = Q_tg(I, JP, K, 1) + Q_tg(IP, JP, K, 1)           ! U at i, J + 1, k
                U_CCT_JM(1) = Q_tg(I, JM, K, 1) + Q_tg(IP, JM, K, 1)           ! U at i, J - 1, k
                U_CCT_KP(1) = Q_tg(I, J, KP, 1) + Q_tg(IP, J, KP, 1)           ! U at i, j, k+ 1
                U_CCT_KM(1) = Q_tg(I, J, KM, 1) + Q_tg(IP, J, KM, 1)           ! U at i, j, K - 1

                U_CCT(3) = (Q_tg(I, J, K, 3)  + Q_tg(I, J, KP, 3)) * RCCI1(JJ)   ! W at i, j, k
                U_CCT_IP(3) = (Q_tg(IP, J, K, 3) + Q_tg(IP, J, KP, 3)) * RCCI1(JJ)  ! W at I + 1, j, k
                U_CCT_IM(3) = (Q_tg(IM, J, K, 3) + Q_tg(IM, J, KP, 3)) * RCCI1(JJ)  ! W at I - 1, j, k
                U_CCT_JP(3) = (Q_tg(I, JP, K, 3) + Q_tg(I, JP, KP, 3)) * RCCI1(JJP) ! W at i, J + 1, k
                U_CCT_JM(3) = (Q_tg(I, JM, K, 3) + Q_tg(I, JM, KP, 3)) * RCCI1(JJM) ! W at i, J - 1, k

                IF(JJ == 1 .AND. iCase == iPIPEC) THEN
                    U_CCT(2) = (Q_tg(I, JP, K, 2) - Q_tg(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_tg(I, JP, K, 2) * RNDI1(JJP)  ! V at i, j, k
                    U_CCT_IP(2) = (Q_tg(IP, JP, K, 2) - Q_tg(IP, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_tg(IP, JP, K, 2) * RNDI1(JJP)  ! V at I + 1, j, k
                    U_CCT_IM(2) = (Q_tg(IM, JP, K, 2) - Q_tg(IM, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_tg(IM, JP, K, 2) * RNDI1(JJP)  ! V at I - 1, j, k
                    U_CCT_KP(2) = (Q_tg(I, JP, KP, 2) - Q_tg(I, JP, KSP, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_tg(I, JP, KP, 2) * RNDI1(JJP)  ! V at i, j, k+ 1
                    U_CCT_KM(2) = (Q_tg(I, JP, KM, 2) - Q_tg(I, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_tg(I, JP, KM, 2) * RNDI1(JJP)  ! V at i, j, K - 1
                ELSE
                    U_CCT(2) = (Q_tg(I, J, K, 2)  + Q_tg(I, JP, K, 2) ) * RCCI1(JJ) ! V at i, j, k
                    U_CCT_IP(2) = (Q_tg(IP, J, K, 2) + Q_tg(IP, JP, K, 2)) * RCCI1(JJ) ! V at I + 1, j, k
                    U_CCT_IM(2) = (Q_tg(IM, J, K, 2) + Q_tg(IM, JP, K, 2)) * RCCI1(JJ) ! V at I - 1, j, k
                    U_CCT_KP(2) = (Q_tg(I, J, KP, 2) + Q_tg(I, JP, KP, 2)) * RCCI1(JJ) ! V at i, j, k+ 1
                    U_CCT_KM(2) = (Q_tg(I, J, KM, 2) + Q_tg(I, JP, KM, 2)) * RCCI1(JJ) ! V at i, j, K - 1
                END IF

                !================= MEAN U, V,W,P =========================
                DO M = 1, NDV
                    U1xzL_tg(J, M) = U1xzL_tg(J, M) + U_CCT(M)
                    UPxzL_tg(J, M) = UPxzL_tg(J, M) + U_CCT(M) * PR_tg(I, J, K)
                END DO
                U1xzL_tg(J, NDV + 1) = U1xzL_tg(J, NDV + 1) + PR_tg(I, J, K)

                !==== MEAN UU,UV,UW,ETC...===============================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        L = (M * (7-M)) / 2 + N - 3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                        U2xzL_tg(J, L) = U2xzL_tg(J, L) +  U_CCT(M) * U_CCT(N)
                    END DO
                END DO

                !====== MEAN UUU, UUV, UUW, UVV, UVW ETC...================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        DO H = 1, NDV
                            IF(N >  H) CYCLE
                            L = M * (6-M) + (N * (7-N)) / 2 + H-8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                            U3xzL_tg(J, L) = U3xzL_tg(J, L) +  U_CCT(M) * U_CCT(N) * U_CCT(H)
                        END DO
                    END DO
                END DO


                !======================== DU / DX.DY.DZ ======================================
                DVDL_CCT(1, 1) =  Q_tg(IP, J, K, 1) - Q_tg(I, J, K, 1)      ! * DXI
                DVDL_CCT(1, 2) = ( U_CCT_JP(1) - U_CCT(1) ) * DYCI(JJP) + &
                ( U_CCT(1) - U_CCT_JM(1) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                DVDL_CCT(1, 3) =  U_CCT_KP(1) - U_CCT_KM(1)           ! * 0.5_WP * 0.5_WP * DZI

                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(2, 1) =  U_CCT_IP(2) - U_CCT_IM(2)           ! * 0.5_WP * 0.5_WP * DXI
                IF(JJ == 1 .AND. iCase == IPIPEC) THEN
                    DVDL_CCT(2, 2) = (   Q_tg(I, JP, K, 2) * RNDI1(JJP) - &
                    ( Q_tg(I, JP, K, 2) - Q_tg(I, JP, KS, 2) ) * 0.50_WP * RNDI1(JJP) ) * DYFI(JJ)
                ELSE
                    DVDL_CCT(2, 2) = ( Q_tg(I, JP, K, 2) * RNDI1(JJP) - Q_tg(I, J, K, 2) * RNDI1(JJ) ) * DYFI(JJ)
                END IF
                DVDL_CCT(2, 3) =  U_CCT_KP(2) - U_CCT_KM(2)           ! * 0.5_WP * 0.5_WP * DZI

                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(3, 1) =  U_CCT_IP(3) - U_CCT_IM(3)           ! * 0.5_WP * 0.5_WP * DXI
                DVDL_CCT(3, 2) = ( U_CCT_JP(3) - U_CCT(3) ) * DYCI(JJP) + &
                ( U_CCT(3) - U_CCT_JM(3) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                DVDL_CCT(3, 3) =  Q_tg(I, J, KP, 3) * RCCI1(JJ) - Q_tg(I, J, K, 3) * RCCI1(JJ)      ! * DZI

                !=========== Mean dU / DX,  P dU / DX =======================================
                DO M = 1, NDV
                    DO N = 1, NDV
                        DVDL1xzL_tg(J, M, N) = DVDL1xzL_tg(J, M, N) + DVDL_CCT(M, N)
                        DVDLPxzL_tg(J, M, N) = DVDLPxzL_tg(J, M, N) + DVDL_CCT(M, N) * PR_tg(I, J, K)
                    END DO
                END DO

                !=========== Mean dUI / DX * dUJ / DX ==========================================
                !                    DO M = 1, NDV
                !                        DO N = 1, NDV
                !                            DO H = 1, NDV
                !                                DO P = 1, NDV
                !                                    L1 = (M - 1) * 3 + H
                !                                    L2 = (N - 1) * 3 + P
                !                                    DVDL2xzL_tg(J, L1, L2) = DVDL2xzL_tg(J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)
                !                                END DO
                !                            END DO
                !                        END DO
                !                    END DO

                DO M = 1, NDV
                    DO H = 1, NDV
                        L1 = (M - 1) * NDV + H
                        DO N = 1, NDV
                            DO P = 1, NDV
                                L2 = (N - 1) * NDV + P
                                DVDL2xzL_tg (J, L1, L2) =  DVDL2xzL_tg (J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)
                            END DO
                        END DO
                    END DO
                END DO


            END DO
        END DO
        U1xzL_tg(J, NDV + 1) = U1xzL_tg(J, NDV + 1) * COE0
        U1xzL_tg(J, 1 : NDV) = U1xzL_tg(J, 1 : NDV) * COE1

        UPxzL_tg(J, :) = UPxzL_tg(J, :) * COE1
        U2xzL_tg(J, :) = U2xzL_tg(J, :) * COE2
        U3xzL_tg(J, :) = U3xzL_tg(J, :) * COE3


        DO M = 1, NDV
            DO N = 1, NDV
                DVDL1xzL_tg(J, M, N) = DVDL1xzL_tg(J, M, N) * COG1(M, N)
                DVDLPxzL_tg(J, M, N) = DVDLPxzL_tg(J, M, N) * COG1(M, N)
            END DO
        END DO

        !            DO M = 1, NDV
        !                DO N = 1, NDV
        !                    DO H = 1, NDV
        !                        DO P = 1, NDV
        !                            L1 = (M - 1) * 3 + H
        !                            L2 = (N - 1) * 3 + P
        !                            DVDL2xzL_tg(J, L1, L2) = DVDL2xzL_tg(J, L1, L2) * COG2(L1, L2)
        !                        END DO
        !                    END DO
        !                END DO
        !            END DO

        DO M = 1, NDV
            DO H = 1, NDV
                L1 = (M - 1) * NDV + H
                DO N = 1, NDV
                    DO P = 1, NDV
                        L2 = (N - 1) * NDV + P
                        DVDL2xzL_tg(J, L1, L2) = DVDL2xzL_tg(J, L1, L2) * COG2(L1, L2)
                    END DO
                END DO
            END DO
        END DO


    END DO


    IF(iPPQuadrants == 1) CALL PP_QUADRANTANALYSIS_xzL(ITG)

    !============================================

    U1mean_tg = 0.0_WP
    U1maxx_tg = 0.0_WP
    DO J = 1, N2DO(MYID)  !@
        JJ = JCL2G(J)
        U1mean_tg = U1mean_tg + U1xzL_tg(J, 1) / DYFI(JJ) / RCCI1(JJ)
        U1maxx_tg = DMAX1(U1maxx_tg, DABS(U1xzL_tg(J, 1)))
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(U1mean_tg, U1mean_WORK_tg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(U1maxx_tg, U1maxx_WORK_tg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    U1mean_WORK_tg = U1mean_WORK_tg / DZI * DBLE(NCL3) /Area_inlet

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_ZX_FLOW_Xperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, M, N, H, L, KS, KSM, KSP,P, L1, L2
    REAL(WP) :: U_CCT(NDV), G_CCT(NDV)
    REAL(WP) :: U_CCT_IP(NDV), U_CCT_JP(NDV), U_CCT_KP(NDV)
    REAL(WP) :: U_CCT_IM(NDV), U_CCT_JM(NDV), U_CCT_KM(NDV)
    REAL(WP) :: DVDL_CCT(NDV, NDV)
    REAL(WP) :: dHdL_CCT(NDV), dTdL_CCT(NDV)
    REAL(WP) :: RTMP
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_DHT

    !IF (MYID == 0) CALL CHKHDL('==> CALL PP_MEAN_ZX_FLOW_Xperiodic_io', MYID)
    DO J = 1, N2DO(MYID)
        JP = JLPV(J)
        JM = JLMV(J)
        JJ = JCL2G(J)
        JJP = JGPV(JJ)
        JJM = JGMV(JJ)

        U1xzL_io(J, :) = 0.0_WP
        G1xzL_io(J, :) = 0.0_WP
        UPxzL_io(J, :) = 0.0_WP

        U2xzL_io(J, :) = 0.0_WP
        UGxzL_io(J, :) = 0.0_WP

        U3xzL_io(J, :) = 0.0_WP
        UGUxzL_io(J, :) = 0.0_WP

        DVDL1xzL_io(J, :, :) = 0.0_WP
        DVDLPxzL_io(J, :, :) = 0.0_WP
        DVDL2xzL_io(J, :, :) = 0.0_WP


        IF(iThermoDynamics == 1) THEN
            T1xzL_io(J) = 0.0_WP
            D1xzL_io(J) = 0.0_WP
            H1xzL_io(J) = 0.0_WP
            M1xzL_io(J) = 0.0_WP

            T2xzL_io(J) = 0.0_WP
            D2xzL_io(J) = 0.0_WP
            H2xzL_io(J) = 0.0_WP

            DHxzL_io(J) = 0.0_WP
            PHxzL_io(J) = 0.0_WP

            DVDL1MxzL_io(J, :,:  ) = 0.0_WP
            DVDL1MHxzL_io(J, :,:  ) = 0.0_WP
            DVDL1MUxzL_io(J, :, :, :) = 0.0_WP
            DVDL2MxzL_io(J, :,:  ) = 0.0_WP

            UHxzL_io(J, :) = 0.0_WP
            GHxzL_io(J, :) = 0.0_WP
            U2DHxzL_io(J, :) = 0.0_WP

            DhDL1xzL_io(J, :) = 0.0_WP
            DhDLPxzL_io(J, :) = 0.0_WP
            DTDLKxzL_io(J, :) = 0.0_WP
            DTDLKUxzL_io(J, :, :) = 0.0_WP
            DTDLKDVDLxzL_io(J, :, :, :) = 0.0_WP
            DHDLMDVDLxzL_io(J, :, :, :) = 0.0_WP
        END IF


        DO I = 1, NCL1_io
            IP = IPV_io(I)
            IM = IMV_io(I)

            DO K = 1, NCL3
                KP = KPV(K)
                KM = KMV(K)
                KS = KSYM(K)
                KSP = KSYM(KP)
                KSM = KSYM(KM)

                G_CCT(1) = ( G_io(I, J, K, 1)  + G_io(IP, J, K, 1)  ) * 0.5_WP ! RHO * U at i, j, k
                U_CCT(1) = ( Q_io(I, J, K, 1)  + Q_io(IP, J, K, 1)  ) * 0.5_WP ! U at i, j, k
                U_CCT_JP(1) = ( Q_io(I, JP, K, 1) + Q_io(IP, JP, K, 1) ) * 0.5_WP ! U at i, J + 1, k
                U_CCT_JM(1) = ( Q_io(I, JM, K, 1) + Q_io(IP, JM, K, 1) ) * 0.5_WP ! U at i, J - 1, k
                U_CCT_KP(1) = ( Q_io(I, J, KP, 1) + Q_io(IP, J, KP, 1) ) * 0.5_WP ! U at i, j, k+ 1
                U_CCT_KM(1) = ( Q_io(I, J, KM, 1) + Q_io(IP, J, KM, 1) ) * 0.5_WP ! U at i, j, K - 1

                G_CCT(3) = ( (G_io(I, J, K, 3)  + G_io(I, J, KP, 3)) * RCCI1(JJ)   ) * 0.5_WP  ! RHO * U at i, j, k
                U_CCT(3) = ( (Q_io(I, J, K, 3)  + Q_io(I, J, KP, 3)) * RCCI1(JJ)   ) * 0.5_WP  ! W at i, j, k
                U_CCT_IP(3) = ( (Q_io(IP, J, K, 3) + Q_io(IP, J, KP, 3)) * RCCI1(JJ)  ) * 0.5_WP  ! W at I + 1, j, k
                U_CCT_IM(3) = ( (Q_io(IM, J, K, 3) + Q_io(IM, J, KP, 3)) * RCCI1(JJ)  ) * 0.5_WP  ! W at I - 1, j, k
                U_CCT_JP(3) = ( (Q_io(I, JP, K, 3) + Q_io(I, JP, KP, 3)) * RCCI1(JJP) ) * 0.5_WP  ! W at i, J + 1, k
                U_CCT_JM(3) = ( (Q_io(I, JM, K, 3) + Q_io(I, JM, KP, 3)) * RCCI1(JJM) ) * 0.5_WP  ! W at i, J - 1, k


                IF(JJ == 1 .AND. iCase == iPIPEC) THEN
                    G_CCT(2) = ( (G_io(I, JP, K, 2) - G_io(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    G_io(I, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, k
                    U_CCT(2) = ( (Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, k

                    U_CCT_IP(2) = ( (Q_io(IP, JP, K, 2) - Q_io(IP, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(IP, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at I + 1, j, k
                    U_CCT_IM(2) = ( (Q_io(IM, JP, K, 2) - Q_io(IM, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(IM, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at I - 1, j, k
                    U_CCT_KP(2) = ( (Q_io(I, JP, KP, 2) - Q_io(I, JP, KSP, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, KP, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, k+ 1
                    U_CCT_KM(2) = ( (Q_io(I, JP, KM, 2) - Q_io(I, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, KM, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, K - 1
                ELSE
                    G_CCT(2) = ( (G_io(I, J, K, 2) + G_io(I, JP, K, 2)) * RCCI1(JJ)   ) * 0.5_WP  ! RHO * U at i, j, k
                    U_CCT(2) = ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * RCCI1(JJ)   ) * 0.5_WP  ! V at i, j, k
                    U_CCT_IP(2) = ( (Q_io(IP, J, K, 2) + Q_io(IP, JP, K, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at I + 1, j, k
                    U_CCT_IM(2) = ( (Q_io(IM, J, K, 2) + Q_io(IM, JP, K, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at I - 1, j, k
                    U_CCT_KP(2) = ( (Q_io(I, J, KP, 2) + Q_io(I, JP, KP, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at i, j, k+ 1
                    U_CCT_KM(2) = ( (Q_io(I, J, KM, 2) + Q_io(I, JP, KM, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at i, j, K - 1
                END IF


                !================U, V,W,P, G1,G2,G3 ========================
                DO M = 1, NDV
                    U1xzL_io(J, M) = U1xzL_io(J, M) + U_CCT(M)
                    G1xzL_io(J, M) = G1xzL_io(J, M) + G_CCT(M)                 !++Method1 ++
                    !G1xzL_io(J, M) = G1xzL_io(J, M) + U_CCT(M) * DENSITY(I, J, K)!++Method2 ++  for rE -pp will introduce errors
                    UPxzL_io(J, M) = UPxzL_io(J, M) + U_CCT(M) * PR_io(I, J, K)
                END DO
                U1xzL_io(J, NDV + 1) = U1xzL_io(J, NDV + 1) + PR_io(I, J, K)

                !====U* FLUX ===============================================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        L = (M * (7-M)) / 2 + N - 3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                        UGxzL_io(J, L) = UGxzL_io(J, L) +  U_CCT(M) * G_CCT(N)                 !++Method1 ++
                        !UGxzL_io(J, L) = UGxzL_io(J, L) +  U_CCT(M) * U_CCT(N) * DENSITY(I, J, K)!++Method2 ++
                        U2xzL_io(J, L) = U2xzL_io(J, L) +  U_CCT(M) * U_CCT(N)
                    END DO
                END DO

                !====U* FLUX * U============================================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        DO H = 1, NDV
                            IF(N >  H) CYCLE
                            L = M * (6 - M) + (N * (7 - N)) / 2 + H - 8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                            UGUxzL_io(J, L) = UGUxzL_io(J, L) +  U_CCT(M) * G_CCT(N) * U_CCT(H)                 !++Method1 ++
                            !UGUxzL_io(J, L) = UGUxzL_io(J, L) +  U_CCT(M) * U_CCT(N) * U_CCT(H) * DENSITY(I, J, K)!++Method2 ++
                        END DO
                    END DO
                END DO

                !====== MEAN UUU, UUV, UUW, UVV, UVW ETC...================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        DO H = 1, NDV
                            IF(N >  H) CYCLE
                            L = M * (6 - M) + (N * (7 - N)) / 2 + H - 8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                            U3xzL_io(J, L) = U3xzL_io(J, L) +  U_CCT(M) * U_CCT(N) * U_CCT(H)
                        END DO
                    END DO
                END DO
                ! U3(1) = U1 U1 U1
                ! U3(2) = U1 U1 U2
                ! U3(3) = U1 U1 U3
                ! U3(4) = U1 U2 U2
                ! U3(5) = U1 U2 U3
                ! U3(6) = U1 U3 U3
                ! U3(7) = U2 U2 U2
                ! U3(8) = U2 U2 U3
                ! U3(9) = U2 U3 U3
                ! U3(10) = U3 U3 U3



                !======================== DU / DX.DY.DZ ======================================
                DVDL_CCT(1, 1) = ( Q_io(IP, J, K, 1) - Q_io(I, J, K, 1) ) * DXI
                DVDL_CCT(1, 3) = ( U_CCT_KP(1) - U_CCT_KM(1) ) * DZI * 0.5_WP
                DVDL_CCT(1, 2) = ( ( U_CCT_JP(1) - U_CCT(1)    ) * DYCI(JJP) + &
                ( U_CCT(1)    - U_CCT_JM(1) ) * DYCI(JJ) ) * 0.5_WP


                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(2, 1) = ( U_CCT_IP(2) - U_CCT_IM(2) ) * DXI * 0.5_WP
                DVDL_CCT(2, 3) = ( U_CCT_KP(2) - U_CCT_KM(2) ) * DZI * 0.5_WP
                IF(JJ == 1 .AND. iCase == IPIPEC) THEN
                    DVDL_CCT(2, 2) = ( Q_io(I, JP, K, 2) * RNDI1(JJP) - &
                    ( Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2) ) * 0.50_WP * RNDI1(JJP) ) * DYFI(JJ)
                ELSE
                    DVDL_CCT(2, 2) = ( Q_io(I, JP, K, 2) * RNDI1(JJP) - Q_io(I, J, K, 2) * RNDI1(JJ) ) * DYFI(JJ)
                END IF


                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(3, 3) = ( Q_io(I, J, KP, 3) * RCCI1(JJ) - Q_io(I, J, K, 3) * RCCI1(JJ) ) * DZI
                DVDL_CCT(3, 1) = ( U_CCT_IP(3) - U_CCT_IM(3) ) * DXI * 0.5_WP
                DVDL_CCT(3, 2) = ( ( U_CCT_JP(3) - U_CCT(3)    ) * DYCI(JJP) + &
                ( U_CCT(3)    - U_CCT_JM(3) ) * DYCI(JJ) ) * 0.5_WP

                !=============== D(u_m) / D(x_n) ===============================
                ! Eq. DVDL1xzL_io(J, M, N) = d(u_m) / D(x_n)
                !     DVDLPxzL_io(J, M, N) = P * D(u_m) / D(x_n)
                DO M = 1, NDV
                    DO N = 1, NDV
                        DVDL1xzL_io(J, M, N) = DVDL1xzL_io(J, M, N) + DVDL_CCT(M, N)
                        DVDLPxzL_io(J, M, N) = DVDLPxzL_io(J, M, N) + DVDL_CCT(M, N) * PR_io(I, J, K)
                    END DO
                END DO

                !====== D(u_m) / D(x_h) * d(u_n) / D(x_p) =========================
                ! Eq. DVDL2xzL_io(J, (M - 1) * 3 + H, (N - 1) * 3 + P) = d(u_m) / D(x_h) * d(u_n) / D(x_p)
                !                    DO M = 1, NDV
                !                        DO N = 1, NDV
                !                            DO H = 1, NDV
                !                                DO P = 1, NDV
                !                                    L1 = (M - 1) * 3 + H ! 1, 2, 3, 4, 5,6,7,8,9
                !                                    L2 = (N - 1) * 3 + P ! 1, 2, 3, 4, 5,6,7,8,9
                !                                    DVDL2xzL_io(J, L1, L2) =  DVDL2xzL_io(J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)
                !                                END DO
                !                            END DO
                !                        END DO
                !                    END DO
                !du_M / DX_h * du_N / DX_p
                DO M = 1, NDV
                    DO H = 1, NDV
                        L1 = (M - 1) * NDV + H
                        DO N = 1, NDV
                            DO P = 1, NDV
                                L2 = (N - 1) * NDV + P
                                DVDL2xzL_io(J, L1, L2) =  DVDL2xzL_io(J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)

                                ! M = 1, H = 1, N = 1, M = 1:3
                                ! DVDL2xzL_io(J, 1, 1) = dU / DX * dU / DX  (J, (1-C1) * NDV + 1, (1-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 1, 2) = dU / DX * dU / Dy  (J, (1-C1) * NDV + 1, (1-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 1, 3) = dU / DX * dU / Dz  (J, (1-C1) * NDV + 1, (1-C1) * NDV +3 )
                                ! M = 1, H = 1, N = 2, M = 1:3
                                ! DVDL2xzL_io(J, 1, 4) = dU / DX * dv/ DX  (J, (1-C1) * NDV + 1, (2 -C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 1, 5) = dU / DX * dv/ Dy  (J, (1-C1) * NDV + 1, (2 -C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 1,6) = dU / DX * dv/ Dz  (J, (1-C1) * NDV + 1, (2 -C1) * NDV +3 )
                                ! M = 1, H = 1, N = 3, M = 1:3
                                ! DVDL2xzL_io(J, 1,7) = dU / DX * dw/ DX  (J, (1-C1) * NDV + 1, (3-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 1,8) = dU / DX * dw/ Dy  (J, (1-C1) * NDV + 1, (3-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 1,9) = dU / DX * dw/ Dz  (J, (1-C1) * NDV + 1, (3-C1) * NDV +3 )

                                ! M = 1, H = 2, N = 1, M = 1:3
                                ! DVDL2xzL_io(J, 2, 1) = dU / Dy * dU / DX  (J, (1-C1) * NDV + 2, (1-C1) * NDV + 1 ) = DVDL2xzL_io(J, 1, 2)
                                ! DVDL2xzL_io(J, 2, 2) = dU / Dy * dU / Dy  (J, (1-C1) * NDV + 2, (1-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 2, 3) = dU / Dy * dU / Dz  (J, (1-C1) * NDV + 2, (1-C1) * NDV +3 )
                                ! M = 1, H = 2, N = 2, M = 1:3
                                ! DVDL2xzL_io(J, 2, 4) = dU / Dy * dv/ DX  (J, (1-C1) * NDV + 2, (2 -C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 2, 5) = dU / Dy * dv/ Dy  (J, (1-C1) * NDV + 2, (2 -C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 2,6) = dU / Dy * dv/ Dz  (J, (1-C1) * NDV + 2, (2 -C1) * NDV +3 )
                                ! M = 1, H = 2, N = 3, M = 1:3
                                ! DVDL2xzL_io(J, 2,7) = dU / Dy * dw/ DX  (J, (1-C1) * NDV + 2, (3-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 2,8) = dU / Dy * dw/ Dy  (J, (1-C1) * NDV + 2, (3-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 2,9) = dU / Dy * dw/ Dz  (J, (1-C1) * NDV + 2, (3-C1) * NDV +3 )

                                ! M = 1, H = 3, N = 1, M = 1:3
                                ! DVDL2xzL_io(J, 3, 1) = dU / Dz * dU / DX  (J, (1-C1) * NDV +3, (1-C1) * NDV + 1 ) = DVDL2xzL_io(J, 1, 3)
                                ! DVDL2xzL_io(J, 3, 2) = dU / Dz * dU / Dy  (J, (1-C1) * NDV +3, (1-C1) * NDV + 2 ) = DVDL2xzL_io(J, 2, 3)
                                ! DVDL2xzL_io(J, 3, 3) = dU / Dz * dU / Dz  (J, (1-C1) * NDV +3, (1-C1) * NDV +3 )
                                ! M = 1, H = 3, N = 2, M = 1:3
                                ! DVDL2xzL_io(J, 3, 4) = dU / Dz * dv/ DX  (J, (1-C1) * NDV +3, (2 -C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 3, 5) = dU / Dz * dv/ Dy  (J, (1-C1) * NDV +3, (2 -C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 3,6) = dU / Dz * dv/ Dz  (J, (1-C1) * NDV +3, (2 -C1) * NDV +3 )
                                ! M = 1, H = 3, N = 3, M = 1:3
                                ! DVDL2xzL_io(J, 3,7) = dU / Dz * dw/ DX  (J, (1-C1) * NDV +3, (3-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 3,8) = dU / Dz * dw/ Dy  (J, (1-C1) * NDV +3, (3-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 3,9) = dU / Dz * dw/ Dz  (J, (1-C1) * NDV +3, (3-C1) * NDV +3 )

                                ! M = 2, H = 1, N = 1, M = 1:3
                                ! DVDL2xzL_io(J, 4, 1) = dv/ DX * dU / DX  (J, (2 -C1) * NDV + 1, (1-C1) * NDV + 1 ) = DVDL2xzL_io(J, 1, 4)
                                ! DVDL2xzL_io(J, 4, 2) = dv/ DX * dU / Dy  (J, (2 -C1) * NDV + 1, (1-C1) * NDV + 2 ) = DVDL2xzL_io(J, 2, 4)
                                ! DVDL2xzL_io(J, 4, 3) = dv/ DX * dU / Dz  (J, (2 -C1) * NDV + 1, (1-C1) * NDV +3 ) = DVDL2xzL_io(J, 3, 4)
                                ! M = 2, H = 1, N = 2, M = 1:3
                                ! DVDL2xzL_io(J, 4, 4) = dv/ DX * dv/ DX  (J, (2 -C1) * NDV + 1, (2 -C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 4, 5) = dv/ DX * dv/ Dy  (J, (2 -C1) * NDV + 1, (2 -C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 4,6) = dv/ DX * dv/ Dz  (J, (2 -C1) * NDV + 1, (2 -C1) * NDV +3 )
                                ! M = 2, H = 1, N = 3, M = 1:3
                                ! DVDL2xzL_io(J, 4,7) = dv/ DX * dw/ DX  (J, (2 -C1) * NDV + 1, (3-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 4,8) = dv/ DX * dw/ Dy  (J, (2 -C1) * NDV + 1, (3-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 4,9) = dv/ DX * dw/ Dz  (J, (2 -C1) * NDV + 1, (3-C1) * NDV +3 )

                                ! M = 2, H = 2, N = 1, M = 1:3
                                ! DVDL2xzL_io(J, 5, 1) = dv/ Dy * dU / DX  (J, (2 -C1) * NDV + 2, (1-C1) * NDV + 1 ) = DVDL2xzL_io(J, 1, 5)
                                ! DVDL2xzL_io(J, 5, 2) = dv/ Dy * dU / Dy  (J, (2 -C1) * NDV + 2, (1-C1) * NDV + 2 ) = DVDL2xzL_io(J, 2, 5)
                                ! DVDL2xzL_io(J, 5, 3) = dv/ Dy * dU / Dz  (J, (2 -C1) * NDV + 2, (1-C1) * NDV +3 ) = DVDL2xzL_io(J, 3, 5)
                                ! M = 2, H = 2, N = 2, M = 1:3
                                ! DVDL2xzL_io(J, 5, 4) = dv/ Dy * dv/ DX  (J, (2 -C1) * NDV + 2, (2 -C1) * NDV + 1 ) = DVDL2xzL_io(J, 4, 5)
                                ! DVDL2xzL_io(J, 5, 5) = dv/ Dy * dv/ Dy  (J, (2 -C1) * NDV + 2, (2 -C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 5,6) = dv/ Dy * dv/ Dz  (J, (2 -C1) * NDV + 2, (2 -C1) * NDV +3 )
                                ! M = 2, H = 2, N = 3, M = 1:3
                                ! DVDL2xzL_io(J, 5,7) = dv/ Dy * dw/ DX  (J, (2 -C1) * NDV + 2, (3-C1) * NDV + 1 )
                                ! DVDL2xzL_io(J, 5,8) = dv/ Dy * dw/ Dy  (J, (2 -C1) * NDV + 2, (3-C1) * NDV + 2 )
                                ! DVDL2xzL_io(J, 5,9) = dv/ Dy * dw/ Dz  (J, (2 -C1) * NDV + 2, (3-C1) * NDV +3 )

                                ! M = 2, H = 3, N = 1, M = 1:3
                                ! DVDL2xzL_io(J,6, 1) = dv/ Dz * dU / DX = DVDL2xzL_io(J, 1,6)
                                ! DVDL2xzL_io(J,6, 2) = dv/ Dz * dU / Dy = DVDL2xzL_io(J, 2,6)
                                ! DVDL2xzL_io(J,6, 3) = dv/ Dz * dU / Dz = DVDL2xzL_io(J, 3,6)
                                ! M = 2, H = 3, N = 2, M = 1:3
                                ! DVDL2xzL_io(J,6, 4) = dv/ Dz * dv/ DX = DVDL2xzL_io(J, 4,6)
                                ! DVDL2xzL_io(J,6, 5) = dv/ Dz * dv/ Dy = DVDL2xzL_io(J, 5,6)
                                ! DVDL2xzL_io(J,6,6) = dv/ Dz * dv/ Dz
                                ! M = 2, H = 3, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,6,7) = dv/ Dz * dw/ DX
                                ! DVDL2xzL_io(J,6,8) = dv/ Dz * dw/ Dy
                                ! DVDL2xzL_io(J,6,9) = dv/ Dz * dw/ Dz

                                ! M = 3, H = 1, N = 1, M = 1:3
                                ! DVDL2xzL_io(J,7, 1) = dw/ DX * dU / DX = DVDL2xzL_io(J, 1,7)
                                ! DVDL2xzL_io(J,7, 2) = dw/ DX * dU / Dy = DVDL2xzL_io(J, 2,7)
                                ! DVDL2xzL_io(J,7, 3) = dw/ DX * dU / Dz = DVDL2xzL_io(J, 3,7)
                                ! M = 3, H = 1, N = 2, M = 1:3
                                ! DVDL2xzL_io(J,7, 4) = dw/ DX * dv/ DX = DVDL2xzL_io(J, 4,7)
                                ! DVDL2xzL_io(J,7, 5) = dw/ DX * dv/ Dy = DVDL2xzL_io(J, 5,7)
                                ! DVDL2xzL_io(J,7,6) = dw/ DX * dv/ Dz = DVDL2xzL_io(J,6,7)
                                ! M = 3, H = 1, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,7,7) = dw/ DX * dw/ DX
                                ! DVDL2xzL_io(J,7,8) = dw/ DX * dw/ Dy
                                ! DVDL2xzL_io(J,7,9) = dw/ DX * dw/ Dz

                                ! M = 3, H = 2, N = 1, M = 1:3
                                ! DVDL2xzL_io(J,8, 1) = dw/ Dy * dU / DX = DVDL2xzL_io(J, 1,8)
                                ! DVDL2xzL_io(J,8, 2) = dw/ Dy * dU / Dy = DVDL2xzL_io(J, 2,8)
                                ! DVDL2xzL_io(J,8, 3) = dw/ Dy * dU / Dz = DVDL2xzL_io(J, 3,8)
                                ! M = 3, H = 2, N = 2, M = 1:3
                                ! DVDL2xzL_io(J,8, 4) = dw/ Dy * dv/ DX = DVDL2xzL_io(J, 4,8)
                                ! DVDL2xzL_io(J,8, 5) = dw/ Dy * dv/ Dy = DVDL2xzL_io(J, 5,8)
                                ! DVDL2xzL_io(J,8,6) = dw/ Dy * dv/ Dz = DVDL2xzL_io(J,6,8)
                                ! M = 3, H = 2, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,8,7) = dw/ Dy * dw/ DX = DVDL2xzL_io(J,7,8)
                                ! DVDL2xzL_io(J,8,8) = dw/ Dy * dw/ Dy
                                ! DVDL2xzL_io(J,8,9) = dw/ Dy * dw/ Dz

                                ! M = 3, H = 3, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,9, 1) = dw/ Dz * dU / DX = DVDL2xzL_io(J, 1,9)
                                ! DVDL2xzL_io(J,9, 2) = dw/ Dz * dU / Dy = DVDL2xzL_io(J, 2,9)
                                ! DVDL2xzL_io(J,9, 3) = dw/ Dz * dU / Dz = DVDL2xzL_io(J, 3,9)
                                ! M = 3, H = 3, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,9, 4) = dw/ Dz * dv/ DX = DVDL2xzL_io(J, 4,9)
                                ! DVDL2xzL_io(J,9, 5) = dw/ Dz * dv/ Dy = DVDL2xzL_io(J, 5,9)
                                ! DVDL2xzL_io(J,9,6) = dw/ Dz * dv/ Dz = DVDL2xzL_io(J,6,9)
                                ! M = 3, H = 3, N = 3, M = 1:3
                                ! DVDL2xzL_io(J,9,7) = dw/ Dz * dw/ DX = DVDL2xzL_io(J,7,9)
                                ! DVDL2xzL_io(J,9,8) = dw/ Dz * dw/ Dy = DVDL2xzL_io(J,8,9)
                                ! DVDL2xzL_io(J,9,9) = dw/ Dz * dw/ Dz
                            END DO
                        END DO
                    END DO
                END DO



                IF(iThermoDynamics == 1) THEN

                    T1xzL_io(J) = T1xzL_io(J) + TEMPERATURE(I, J, K)
                    D1xzL_io(J) = D1xzL_io(J) + DENSITY    (I, J, K)
                    H1xzL_io(J) = H1xzL_io(J) + ENTHALPY   (I, J, K)
                    M1xzL_io(J) = M1xzL_io(J) + Viscousity (I, J, K)

                    T2xzL_io(J) = T2xzL_io(J) + TEMPERATURE(I, J, K) * TEMPERATURE(I, J, K)
                    D2xzL_io(J) = D2xzL_io(J) + DENSITY    (I, J, K) * DENSITY    (I, J, K)
                    H2xzL_io(J) = H2xzL_io(J) + ENTHALPY   (I, J, K) * ENTHALPY   (I, J, K)

                    DHxzL_io(J) = DHxzL_io(J) + ENTHALPY(I, J, K) * DENSITY(I, J, K)
                    PHxzL_io(J) = PHxzL_io(J) + ENTHALPY(I, J, K) * PR_io(I, J, K)

                    !===== D(h) / D(x_i) ==============================================
                    DhDL_CCT(1) = ( ENTHALPY(IP, J, K) - ENTHALPY(IM, J, K) ) * DXI * 0.5_WP
                    DhDL_CCT(3) = ( ENTHALPY(I, J, KP) - ENTHALPY(I, J, KM) ) * DZI * 0.5_WP
                    DhDL_CCT(2) = ( ( ENTHALPY(I, JP, K) - ENTHALPY(I, J, K ) ) * DYCI(JJP) + &
                    ( ENTHALPY(I, J, K) - ENTHALPY(I, JM, K ) ) * DYCI(JJ) ) * 0.5_WP

                    !===== D(T) / D(x_i) ==============================================
                    DTDL_CCT(1) = ( TEMPERATURE(IP, J, K) - TEMPERATURE(IM, J, K) ) * DXI * 0.5_WP
                    DTDL_CCT(3) = ( TEMPERATURE(I, J, KP) - TEMPERATURE(I, J, KM) ) * DZI * 0.5_WP
                    DTDL_CCT(2) = ( ( TEMPERATURE(I, JP, K) - TEMPERATURE(I, J, K ) ) * DYCI(JJP) + &
                    ( TEMPERATURE(I, J, K) - TEMPERATURE(I, JM, K ) ) * DYCI(JJ) ) * 0.5_WP

                    !=============== D(u_m) / D(x_n) ===============================
                    ! Eq.
                    !     DVDL1MxzL_io(J, M, N) = M * D(u_m) / D(x_n)
                    !     DVDL1MhxzL_io(J, M, N) = M * D(u_m) / D(x_n) * H
                    !     DVDL1MUxzL_io(J, M, N,H) = M * D(u_m) / D(x_n) * u_h
                    DO M = 1, NDV
                        DO N = 1, NDV
                            DVDL1MxzL_io(J, M, N) = DVDL1MxzL_io(J, M, N) + DVDL_CCT(M, N) * Viscousity(I, J, K)
                            DVDL1MHxzL_io(J, M, N) = DVDL1MHxzL_io(J, M, N) + &
                                                     DVDL_CCT(M, N) * Viscousity(I, J, K) * ENTHALPY(I, J, K)
                            DO H = 1, NDV
                                DVDL1MUxzL_io(J, M, N,H) = DVDL1MUxzL_io(J, M, N,H) + &
                                                           DVDL_CCT(M, N) * Viscousity(I, J, K) * U_CCT(H)
                            END DO
                        END DO
                    END DO

                    !====== D(u_m) / D(x_h) * d(u_n) / D(x_p) =========================
                    ! Eq.
                    !     DVDL2MxzL_io(J, (M - 1) * 3 + H, (N - 1) * 3 + P) = d(u_m) / D(x_h) * d(u_n) / D(x_p) * M
                    !                        DO M = 1, NDV
                    !                            DO N = 1, NDV
                    !                                DO H = 1, NDV
                    !                                    DO P = 1, NDV
                    !                                        L1 = (M - 1) * 3 + H ! 1, 2, 3, 4, 5,6,7,8,9
                    !                                        L2 = (N - 1) * 3 + P ! 1, 2, 3, 4, 5,6,7,8,9
                    !                                        DVDL2MxzL_io(J, L1, L2) =  DVDL2MxzL_io(J, L1, L2) + &
                    !                                                       DVDL_CCT(M,H) * DVDL_CCT(N,P) * Viscousity(I, J, K)
                    !                                    END DO
                    !                                END DO
                    !                            END DO
                    !                        END DO

                    DO M = 1, NDV
                        DO H = 1, NDV
                            L1 = (M - 1) * NDV + H
                            DO N = 1, NDV
                                DO P = 1, NDV
                                    L2 = (N - 1) * NDV + P
                                    DVDL2MxzL_io(J, L1, L2) =  DVDL2MxzL_io(J, L1, L2) + &
                                    DVDL_CCT(M,H) * DVDL_CCT(N,P) * Viscousity(I, J, K)
                                END DO
                            END DO
                        END DO
                    END DO


                    !=========================================
                    !Eq.  UHxzL_io(J, M) = u_m * h
                    !     GHxzL_io(J, M) = g_m * h
                    !     DhDL1xzL_io(J, M)      = d(h) / D(x_m)
                    !     DhDLPxzL_io(J, M)      = d(h) / D(x_m) * p
                    !     DhDLKxzL_io(J, M)      = d(h) / D(x_m) * k
                    !     DTDLKUxzL_io(J, M, N)   = d(h) / D(x_m) * k * u_n
                    !     DTDLKDVDLxzL_io(J, M, N,H) = d(h) / D(x_m) * k * d(u_n) / D(x_h)
                    !     DHDLMDVDLxzL_io(J, M, N,H) = d(h) / D(x_m) * mu* d(u_n) / D(x_h)
                    DO M = 1, NDV

                        UHxzL_io(J, M) = UHxzL_io(J, M) + ENTHALPY(I, J, K) * U_CCT(M)
                        GHxzL_io(J, M) = GHxzL_io(J, M) + ENTHALPY(I, J, K) * G_CCT(M)  !++method1 +++
                        !GHxzL_io(J, M) = GHxzL_io(J, M) + ENTHALPY(I, J, K) * U_CCT(M) * DENSITY(I, J, K) ! ++method2 +++

                        DhDL1xzL_io(J, M) = DhDL1xzL_io(J, M) + DhDL_CCT(M)
                        DhDLPxzL_io(J, M) = DhDLPxzL_io(J, M) + DhDL_CCT(M) * PR_io(I, J, K)
                        DTDLKxzL_io(J, M) = DTDLKxzL_io(J, M) + DTDL_CCT(M) * THERMCONDT(I, J, K)

                        DO N = 1, NDV
                            DTDLKUxzL_io(J, M, N) = DTDLKUxzL_io(J, M, N) + DTDL_CCT(M) * THERMCONDT(I, J, K) * U_CCT(N)
                            DO H = 1, NDV
                                DTDLKDVDLxzL_io(J, M, N,H) = DTDLKDVDLxzL_io(J, M, N,H) + &
                                DTDL_CCT(M) * THERMCONDT(I, J, K) * DVDL_CCT(N,H)
                                DHDLMDVDLxzL_io(J, M, N,H) = DHDLMDVDLxzL_io(J, M, N,H) + &
                                DHDL_CCT(M) * Viscousity(I, J, K) * DVDL_CCT(N,H)
                            END DO
                        END DO
                    END DO

                    !====\rho h * U* U============================================
                    !U2DHxzL_io(J, L) = \rho * h * u_m * u_n
                    DO M = 1, NDV
                        DO N = 1, NDV
                            IF(M >  N) CYCLE
                            L = (M * (7 - M)) / 2 + N - 3
                            U2DHxzL_io(J, L) = U2DHxzL_io(J, L) + DENSITY(I, J, K) * ENTHALPY(I, J, K) * U_CCT(M) * U_CCT(N)
                        END DO
                    END DO

                END IF


            END DO
        END DO

        U1xzL_io(J, :) = U1xzL_io(J, :) * VL1313_io
        G1xzL_io(J, :) = G1xzL_io(J, :) * VL1313_io
        UPxzL_io(J, :) = UPxzL_io(J, :) * VL1313_io

        U2xzL_io(J, :) = U2xzL_io(J, :) * VL1313_io
        UGxzL_io(J, :) = UGxzL_io(J, :) * VL1313_io
        U3xzL_io(J, :) = U3xzL_io(J, :) * VL1313_io
        UGUxzL_io(J, :) = UGUxzL_io(J, :) * VL1313_io

        DVDL1xzL_io(J, :, :) = DVDL1xzL_io(J, :, :) * VL1313_io
        DVDLPxzL_io(J, :, :) = DVDLPxzL_io(J, :, :) * VL1313_io
        DVDL2xzL_io(J, :, :) = DVDL2xzL_io(J, :, :) * VL1313_io

        !WRITE(*, *) J, U1xzL_io(J, 2), G1xzL_io(J, 2) !test
        IF(iThermoDynamics == 1) THEN
            T1xzL_io(J) = T1xzL_io(J) * VL1313_io
            D1xzL_io(J) = D1xzL_io(J) * VL1313_io
            H1xzL_io(J) = H1xzL_io(J) * VL1313_io
            M1xzL_io(J) = M1xzL_io(J) * VL1313_io


            T2xzL_io(J) = T2xzL_io(J) * VL1313_io
            D2xzL_io(J) = D2xzL_io(J) * VL1313_io
            H2xzL_io(J) = H2xzL_io(J) * VL1313_io

            DHxzL_io(J) = DHxzL_io(J) * VL1313_io
            PHxzL_io(J) = PHxzL_io(J) * VL1313_io

            DVDL1MxzL_io(J, :, :) = DVDL1MxzL_io(J, :, :) * VL1313_io
            DVDL1MHxzL_io(J, :, :) = DVDL1MHxzL_io(J, :, :) * VL1313_io
            DVDL1MUxzL_io(J, :, :, :) = DVDL1MUxzL_io(J, :, :, :) * VL1313_io
            DVDL2MxzL_io(J, :, :) = DVDL2MxzL_io(J, :, :) * VL1313_io

            !IF (MYID == 0) THEN
            !                DO L1 = 1, (3- 1) * NDV +3
            !                    DO L2 = 1, (3- 1) * NDV +3
            !                        WRITE(*, *) J, L1, L2, DVDL2xzL_io(J, L1, L2), DVDL2MxzL_io(J, L1, L2)
            !                    END DO
            !                END DO
            !WRITE(*, *) J, DVDLPxzL_io(J, 1, 1), DVDL1xzL_io(J, 1, 1), U1xzL_io(J, 4)
            !END IF

            UHxzL_io(J, :) = UHxzL_io(J, :) * VL1313_io
            GHxzL_io(J, :) = GHxzL_io(J, :) * VL1313_io
            U2DHxzL_io(J, :) = U2DHxzL_io(J, :) * VL1313_io

            DhDL1xzL_io(J, :) = DhDL1xzL_io(J, :) * VL1313_io
            DhDLPxzL_io(J, :) = DhDLPxzL_io(J, :) * VL1313_io
            DTDLKxzL_io(J, :) = DTDLKxzL_io(J, :) * VL1313_io
            DTDLKUxzL_io(J, :, :)    = DTDLKUxzL_io(J, :, :) * VL1313_io
            DTDLKDVDLxzL_io(J, :, :, :) = DTDLKDVDLxzL_io(J, :, :, :) * VL1313_io
            DHDLMDVDLxzL_io(J, :, :, :) = DHDLMDVDLxzL_io(J, :, :, :) * VL1313_io
        END IF

    END DO


    IF(iPPQuadrants == 1) CALL PP_QUADRANTANALYSIS_xzL(IIO)

    CALL PP_DRIVEN_FORCE
    !==============================

    G1rate_io = 0.0_WP
    G1maxx_io = 0.0_WP
    DO J = 1, N2DO(MYID)  !@
        JJ = JCL2G(J)
        G1rate_io = G1rate_io + G1xzL_io(J, 1) / DYFI(JJ) / RCCI1(JJ)
        G1maxx_io = DMAX1(G1maxx_io, DABS(G1xzL_io(J, 1)))
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(G1rate_io, G1rate_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(G1maxx_io, G1maxx_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    G1BULK_WORK_io = G1rate_WORK_iO / DZI * DBLE(NCL3) /Area_inlet

    IF(iThermoDynamics == 1) THEN
        DH1rate_io = 0.0_WP
        T1maxx_io = 0.0_WP
        DO J = 1, N2DO(MYID)  !@
            JJ = JCL2G(J)
            DH1rate_io = DH1rate_io + DHxzL_io(J) / DYFI(JJ) / RCCI1(JJ)
            T1maxx_io = DMAX1(T1maxx_io, DABS(T1xzL_io(J)))
        END DO
        CALL MPI_BARRIER(ICOMM, IERROR)
        CALL MPI_ALLREDUCE(DH1rate_io, DH1rate_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(T1maxx_io, T1maxx_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

        !H1bulk_WORK_io = H1RATE_WORK_iO /G1rate_WORK_io

        T1bulk_WORK_io = spline_interpolation_DHT(DH1bulk_WORK_io)

    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INST_Qio_Gradient_CellCentred(I, J, K, DVDL_CCT)
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: I, J, K
    REAL(WP), INTENT(OUT) :: DVDL_CCT(NDV, NDV)

    INTEGER(4) :: JJ, JM, JP, JJP, JJM, IP, IM, KP, KM, KS, KSP, KSM
    REAL(WP) :: U_CCT(NDV)
    REAL(WP) :: U_CCT_IP(NDV), U_CCT_JP(NDV), U_CCT_KP(NDV)
    REAL(WP) :: U_CCT_IM(NDV), U_CCT_JM(NDV), U_CCT_KM(NDV)


    JJ = JCL2G(J)
    JM = JLMV(J)
    JP = JLPV(J)
    JJP = JGPV(JJ)
    JJM = JGMV(JJ)

    IP = IPV_io(I)
    IM = IMV_io(I)
    KP = KPV(K)
    KM = KMV(K)
    KS = KSYM(K)
    KSP = KSYM(KP)
    KSM = KSYM(KM)

    U_CCT(1) = ( Q_io(I, J, K, 1)  + Q_io(IP, J, K, 1)  ) * 0.5_WP ! U at i, j, k
    U_CCT_JP(1) = ( Q_io(I, JP, K, 1) + Q_io(IP, JP, K, 1) ) * 0.5_WP ! U at i, J + 1, k
    U_CCT_JM(1) = ( Q_io(I, JM, K, 1) + Q_io(IP, JM, K, 1) ) * 0.5_WP ! U at i, J - 1, k
    U_CCT_KP(1) = ( Q_io(I, J, KP, 1) + Q_io(IP, J, KP, 1) ) * 0.5_WP ! U at i, j, k+ 1
    U_CCT_KM(1) = ( Q_io(I, J, KM, 1) + Q_io(IP, J, KM, 1) ) * 0.5_WP ! U at i, j, K - 1

    U_CCT(3) = ( (Q_io(I, J, K, 3)  + Q_io(I, J, KP, 3)) * RCCI1(JJ)   ) * 0.5_WP  ! W at i, j, k
    U_CCT_IP(3) = ( (Q_io(IP, J, K, 3) + Q_io(IP, J, KP, 3)) * RCCI1(JJ)  ) * 0.5_WP  ! W at I + 1, j, k
    U_CCT_IM(3) = ( (Q_io(IM, J, K, 3) + Q_io(IM, J, KP, 3)) * RCCI1(JJ)  ) * 0.5_WP  ! W at I - 1, j, k
    U_CCT_JP(3) = ( (Q_io(I, JP, K, 3) + Q_io(I, JP, KP, 3)) * RCCI1(JJP) ) * 0.5_WP  ! W at i, J + 1, k
    U_CCT_JM(3) = ( (Q_io(I, JM, K, 3) + Q_io(I, JM, KP, 3)) * RCCI1(JJM) ) * 0.5_WP  ! W at i, J - 1, k


    IF(JJ == 1 .AND. iCase == iPIPEC) THEN
        U_CCT(2) = ( (Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
        Q_io(I, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, k

        U_CCT_IP(2) = ( (Q_io(IP, JP, K, 2) - Q_io(IP, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
        Q_io(IP, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at I + 1, j, k
        U_CCT_IM(2) = ( (Q_io(IM, JP, K, 2) - Q_io(IM, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
        Q_io(IM, JP, K, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at I - 1, j, k
        U_CCT_KP(2) = ( (Q_io(I, JP, KP, 2) - Q_io(I, JP, KSP, 2)) * 0.50_WP * RNDI1(JJP) + &
        Q_io(I, JP, KP, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, k+ 1
        U_CCT_KM(2) = ( (Q_io(I, JP, KM, 2) - Q_io(I, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP) + &
        Q_io(I, JP, KM, 2) * RNDI1(JJP)  ) * 0.5_WP  ! V at i, j, K - 1
    ELSE
        U_CCT(2) = ( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * RCCI1(JJ)   ) * 0.5_WP  ! V at i, j, k
        U_CCT_IP(2) = ( (Q_io(IP, J, K, 2) + Q_io(IP, JP, K, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at I + 1, j, k
        U_CCT_IM(2) = ( (Q_io(IM, J, K, 2) + Q_io(IM, JP, K, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at I - 1, j, k
        U_CCT_KP(2) = ( (Q_io(I, J, KP, 2) + Q_io(I, JP, KP, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at i, j, k+ 1
        U_CCT_KM(2) = ( (Q_io(I, J, KM, 2) + Q_io(I, JP, KM, 2)) * RCCI1(JJ) ) * 0.5_WP  ! V at i, j, K - 1
    END IF


    !======================== DU / DX.DY.DZ ======================================
    DVDL_CCT(1, 1) = ( Q_io(IP, J, K, 1) - Q_io(I, J, K, 1) ) * DXI
    DVDL_CCT(1, 3) = ( U_CCT_KP(1) - U_CCT_KM(1) ) * DZI * 0.5_WP
    DVDL_CCT(1, 2) = ( ( U_CCT_JP(1) - U_CCT(1)    ) * DYCI(JJP) + &
    ( U_CCT(1)    - U_CCT_JM(1) ) * DYCI(JJ) ) * 0.5_WP


    !======================== DV/ DX.DY.DZ ======================================
    DVDL_CCT(2, 1) = ( U_CCT_IP(2) - U_CCT_IM(2) ) * DXI * 0.5_WP
    DVDL_CCT(2, 3) = ( U_CCT_KP(2) - U_CCT_KM(2) ) * DZI * 0.5_WP
    IF(JJ == 1 .AND. iCase == IPIPEC) THEN
        DVDL_CCT(2, 2) = (   Q_io(I, JP, K, 2) * RNDI1(JJP) - &
        ( Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2) ) * 0.50_WP * RNDI1(JJP) ) * DYFI(JJ)
    ELSE
        DVDL_CCT(2, 2) = ( Q_io(I, JP, K, 2) * RNDI1(JJP) - Q_io(I, J, K, 2) * RNDI1(JJ) ) * DYFI(JJ)
    END IF


    !======================== DW/ DX.DY.DZ ======================================
    DVDL_CCT(3, 3) = ( Q_io(I, J, KP, 3) * RCCI1(JJ) - Q_io(I, J, K, 3) * RCCI1(JJ) ) * DZI
    DVDL_CCT(3, 1) = ( U_CCT_IP(3) - U_CCT_IM(3) ) * DXI * 0.5_WP
    DVDL_CCT(3, 2) = ( ( U_CCT_JP(3) - U_CCT(3)    ) * DYCI(JJP) + &
    ( U_CCT(3)    - U_CCT_JM(3) ) * DYCI(JJ) ) * 0.5_WP


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_DRIVEN_FORCE
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE THERMAL_INFO
    IMPLICIT NONE

    REAL(WP) :: VisForcWall1, VisForcWall1_WORK
    REAL(WP) :: VisForcWall2, VisForcWall2_WORK
    REAL(WP) :: BuoFc, BuoFc_WORK
    INTEGER(4) :: IC, KC, M, JC, JJ
    REAL(WP) :: B_tmp

    IF(iFlowDriven /= 1) RETURN


    ! VISCOUS FORCE on the bottom wall
    VisForcWall1 = 0.0_WP
    VisForcWall2 = 0.0_WP
    IF(MYID == 0) THEN
        DO IC = 1, NCL1_io
            DO KC = 1, NCL3
                IF(iThermoDynamics == 1) THEN
                    IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
                        Mwal(iBotWall) = M_WAL_GV(NCL1_io / 2, iBotWall)
                    END IF
                    IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
                        CALL THERM_PROP_UPDATE_FROM_DH(DH(IC, 1, KC), Hwal_RA(iBotWall), Twal(iBotWall), Dwal(iBotWall), &
                    Mwal(iBotWall), Kwal(iBotWall), Cpwal(iBotWall), B_tmp)
                    END IF

                    VisForcWall1 = VisForcWall1 + Mwal(iBotWall) * (Q_io(IC, 1, KC, 1) - 0.0_WP) / &
                    (YCC(1) - YND(1)) * CVISC * DX * DZ
                ELSE
                    VisForcWall1 = VisForcWall1 + (Q_io(IC, 1, KC, 1) - 0.0_WP) / (YCC(1) - YND(1)) * CVISC * DX * DZ
                END IF
            END DO
        END DO
        !WRITE(*, *) 'Wall1 Friction Force:', VisForcWall1
    END IF

    IF(MYID == NPSLV) THEN

        DO IC = 1, NCL1_io
            DO KC = 1, NCL3
                IF(iThermoDynamics == 1) THEN
                    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature) THEN
                        Mwal(iTopWall) = M_WAL_GV(NCL1_io / 2, iTopWall)
                    END IF
                    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
                        CALL THERM_PROP_UPDATE_FROM_DH(DH(IC, N2DO(MYID), KC), Hwal_RA(iTopWall), Twal(iTopWall), Dwal(iTopWall), &
                    Mwal(iTopWall), Kwal(iTopWall), Cpwal(iTopWall), B_tmp)
                    END IF

                    VisForcWall2 = VisForcWall2 + &
                    Mwal(iTopWall) * (Q_io(IC, N2DO(MYID), KC, 1) - 0.0_WP) / (YCC(NCL2) - YND(NND2)) * CVISC * DX * DZ
                ELSE
                    VisForcWall2 = VisForcWall2 + (Q_io(IC, N2DO(MYID), KC, 1) - 0.0_WP) / (YCC(NCL2) - YND(NND2)) * CVISC * DX * DZ
                END IF
            END DO
        END DO
        !WRITE(*, *) 'Wall2 Friction Force:', VisForcWall2
    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VisForcWall1, VisForcWall1_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VisForcWall2, VisForcWall2_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

    ! Gravity force
    BuoFc = 0.0_WP
    BuoFc_WORK = 0.0_WP
    IF ( iGravity == 1 ) THEN
        DO JC = 1, N2DO(MYID)
            JJ = JCL2G(JC)
            DO IC = 1, NCL1_io
                DO KC = 1, NCL3
                    BuoFc = BuoFc + DENSITY(IC, JC, KC) * F_A / DYFI(JJ) * DX * DZ
                END DO
            END DO
        END DO
        !WRITE(*, *) 'Gravity in MYID', BuoFc, MYID
        CALL MPI_BARRIER(ICOMM, IERROR)
        CALL MPI_ALLREDUCE(BuoFc,BuoFc_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    END IF

    IF(MYID == 0) THEN
        FcDrv_io = ( DABS(VisForcWall1_WORK) + DABS(VisForcWall2_WORK) -BuoFc_WORK) / 2.0_WP / DBLE(NCL1_io * NCL3) / DX / DZ

        !WRITE(*, *) 'FRC', DABS(VisForcWall1_WORK) / DBLE(NCL1_io * NCL3) / DX / DZ, &
        !                  DABS(VisForcWall2_WORK) / DBLE(NCL1_io * NCL3) / DX / DZ, &
        !                  BuoFc_WORK / DBLE(NCL1_io * NCL3) / DX / DZ, FcDrv_io
    END IF

    CALL MPI_BCAST( FcDrv_io, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    DO JC = 1, N2DO(MYID)
        DO M = 1, NDV
            FUxzL_io(JC, M) = FcDrv_io * U1xzL_io(JC, M)
        END DO
        FUxzL_io(JC, NDV + 1) = FcDrv_io
    END DO



    RETURN
END SUBROUTINE
