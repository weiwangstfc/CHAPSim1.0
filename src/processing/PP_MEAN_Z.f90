!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE: PP_MEAN_Z_FLOW_nonXperiodic_io
!> SUBROUTINE: PP_MEAN_Z_THEML_nonXperiodic_io
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_Z_FLOW_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, M, N, H, L, KS, KSM, KSP, P, L1, L2
    REAL(WP) :: U_CCT(NDV), G_CCT(NDV)
    REAL(WP) :: U_CCT_IP(NDV), U_CCT_JP(NDV), U_CCT_KP(NDV)
    REAL(WP) :: U_CCT_IM(NDV), U_CCT_JM(NDV), U_CCT_KM(NDV)
    REAL(WP) :: DVDL_CCT(NDV, NDV)
    REAL(WP) :: COG1(NDV, NDV), COG2(NDV * NDV, NDV * NDV)
    REAL(WP) :: COE0, COE1, COE2, COE3
    REAL(WP) :: RTMP


    COE0 = 1.0_WP / DBLE(NCL3)
    COE1 = 1.0_WP / DBLE(NCL3) * 0.5_WP
    COE2 = 1.0_WP / DBLE(NCL3) * (0.5_WP**2)
    COE3 = 1.0_WP / DBLE(NCL3) * (0.5_WP**3)

    COG1(1, 1) = 1.0_WP / DBLE(NCL3) * DXI
    COG1(1, 2) = 1.0_WP / DBLE(NCL3) * 0.25_WP
    COG1(1, 3) = 1.0_WP / DBLE(NCL3) * 0.25_WP * DZI

    COG1(2, 1) = 1.0_WP / DBLE(NCL3) * 0.25_WP * DXI
    COG1(2, 2) = 1.0_WP / DBLE(NCL3)
    COG1(2, 3) = 1.0_WP / DBLE(NCL3) * 0.25_WP * DZI

    COG1(3, 1) = 1.0_WP / DBLE(NCL3) * 0.25_WP * DXI
    COG1(3, 2) = 1.0_WP / DBLE(NCL3) * 0.25_WP
    COG1(3, 3) = 1.0_WP / DBLE(NCL3) * DZI


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
                    COG2(L1, L2) = COG1(M,H) * COG1(N,P) * DBLE(NCL3) !Junjie
                END DO
            END DO
        END DO
    END DO




    DO J = 1, N2DO(MYID)
        JP = JLPV(J)
        JM = JLMV(J)
        JJ = JCL2G(J)
        JJP = JGPV(JJ)
        JJM = JGMV(JJ)
        DO I = 1, NCL1_io
            IP = IPV_io(I)
            IM = IMV_io(I)
            U1zL_io(I, J, :) = 0.0_WP
            G1zL_io(I, J, :) = 0.0_WP
            UPzL_io(I, J, :) = 0.0_WP

            U2zL_io(I, J, :) = 0.0_WP
            UGzL_io(I, J, :) = 0.0_WP
            UGUzL_io(I, J, :) = 0.0_WP

            DVDL1zL_io(I, J, :, :) = 0.0_WP
            DVDLPzL_io(I, J, :, :) = 0.0_WP
            DVDL2zL_io(I, J, :, :) = 0.0_WP

            DO K = 1, NCL3
                KP = KPV(K)
                KM = KMV(K)
                KS = KSYM(K)
                KSP = KSYM(KP)
                KSM = KSYM(KM)
                U_CCT(1) = Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)              ! U at i, j, k
                G_CCT(1) = G_io(I, J, K, 1) + G_io(IP, J, K, 1)             ! RHO * U at i, j, k
                U_CCT_JP(1) = Q_io(I, JP, K, 1) + Q_io(IP, JP, K, 1)  ! U at i, J + 1, k
                U_CCT_JM(1) = Q_io(I, JM, K, 1) + Q_io(IP, JM, K, 1)  ! U at i, J - 1, k
                U_CCT_KP(1) = Q_io(I, J, KP, 1) + Q_io(IP, J, KP, 1)  ! U at i, j, k+ 1
                U_CCT_KM(1) = Q_io(I, J, KM, 1) + Q_io(IP, J, KM, 1)  ! U at i, j, K - 1


                U_CCT(3) = (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * RCCI1(JJ)   ! W at i, j, k
                G_CCT(3) = (G_io(I, J, K, 3) + G_io(I, J, KP, 3)) * RCCI1(JJ)  ! RHO * U at i, j, k
                U_CCT_IP(3) = (Q_io(IP, J, K, 3) + Q_io(IP, J, KP, 3)) * RCCI1(JJ)  ! W at I + 1, j, k
                U_CCT_IM(3) = (Q_io(IM, J, K, 3) + Q_io(IM, J, KP, 3)) * RCCI1(JJ)  ! W at I - 1, j, k
                U_CCT_JP(3) = (Q_io(I, JP, K, 3) + Q_io(I, JP, KP, 3)) * RCCI1(JJP) ! W at i, J + 1, k
                U_CCT_JM(3) = (Q_io(I, JM, K, 3) + Q_io(I, JM, KP, 3)) * RCCI1(JJM) ! W at i, J - 1, k


                IF(JJ == 1 .AND. iCase == iPIPEC) THEN
                    U_CCT(2) = (Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, K, 2) * RNDI1(JJP)  ! V at i, j, k
                    G_CCT(2) = (G_io(I, JP, K, 2) - G_io(I, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    G_io(I, JP, K, 2) * RNDI1(JJP)  ! V at i, j, k
                    U_CCT_IP(2) = (Q_io(IP, JP, K, 2) - Q_io(IP, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(IP, JP, K, 2) * RNDI1(JJP)  ! V at I + 1, j, k
                    U_CCT_IM(2) = (Q_io(IM, JP, K, 2) - Q_io(IM, JP, KS, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(IM, JP, K, 2) * RNDI1(JJP)  ! V at I - 1, j, k
                    U_CCT_KP(2) = (Q_io(I, JP, KP, 2) - Q_io(I, JP, KSP, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, KP, 2) * RNDI1(JJP)  ! V at i, j, k+ 1
                    U_CCT_KM(2) = (Q_io(I, JP, KM, 2) - Q_io(I, JP, KSM, 2)) * 0.50_WP * RNDI1(JJP) + &
                    Q_io(I, JP, KM, 2) * RNDI1(JJP)  ! V at i, j, K - 1
                ELSE
                    U_CCT(2) = (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * RCCI1(JJ)   ! V at i, j, k
                    G_CCT(2) = (G_io(I, J, K, 2) + G_io(I, JP, K, 2)) * RCCI1(JJ)  ! RHO * U at i, j, k
                    U_CCT_IP(2) = (Q_io(IP, J, K, 2) + Q_io(IP, JP, K, 2)) * RCCI1(JJ)  ! V at I + 1, j, k
                    U_CCT_IM(2) = (Q_io(IM, J, K, 2) + Q_io(IM, JP, K, 2)) * RCCI1(JJ)  ! V at I - 1, j, k
                    U_CCT_KP(2) = (Q_io(I, J, KP, 2) + Q_io(I, JP, KP, 2)) * RCCI1(JJ)  ! V at i, j, k+ 1
                    U_CCT_KM(2) = (Q_io(I, J, KM, 2) + Q_io(I, JP, KM, 2)) * RCCI1(JJ)  ! V at i, j, K - 1
                END IF


                !================U, V,W,P, G1,G2,G3 ========================
                DO M = 1, NDV
                    U1zL_io(I, J, M) = U1zL_io(I, J, M) + U_CCT(M)
                    G1zL_io(I, J, M) = G1zL_io(I, J, M) + G_CCT(M)
                    UPzL_io(I, J, M) = UPzL_io(I, J, M) + U_CCT(M) * PR_io(I, J, K)
                END DO
                U1zL_io(I, J, NDV + 1) = U1zL_io(I, J, NDV + 1) + PR_io(I, J, K)

                !====U* FLUX ===============================================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        L = (M * (7-M)) / 2 + N - 3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                        UGzL_io(I, J, L) = UGzL_io(I, J, L) +  U_CCT(M) * G_CCT(N)
                        U2zL_io(I, J, L) = U2zL_io(I, J, L) +  U_CCT(M) * U_CCT(N)
                    END DO
                END DO

                !====U* FLUX * U============================================
                DO M = 1, NDV
                    DO N = 1, NDV
                        IF(M >  N) CYCLE
                        DO H = 1, NDV
                            IF(N >  H) CYCLE
                            L = M * (6 - M) + (N * (7-N)) / 2 + H - 8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER tri-matrix
                            UGUzL_io(I, J, L) = UGUzL_io(I, J, L) +  U_CCT(M) * G_CCT(N) * U_CCT(H)
                        END DO
                    END DO
                END DO


                !======================== DU / DX.DY.DZ ======================================
                DVDL_CCT(1, 1) =  Q_io(IP, J, K, 1) - Q_io(I, J, K, 1)    ! * DXI
                DVDL_CCT(1, 2) = ( U_CCT_JP(1) - U_CCT(1) ) * DYCI(JJP) + &
                ( U_CCT(1) - U_CCT_JM(1) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                DVDL_CCT(1, 3) =  U_CCT_KP(1) - U_CCT_KM(1)           ! * 0.5_WP * 0.5_WP * DZI

                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(2, 1) =  U_CCT_IP(2) - U_CCT_IM(2)           ! * 0.5_WP * 0.5_WP * DXI
                IF(JJ == 1 .AND. iCase == IPIPEC) THEN
                    DVDL_CCT(2, 2) = (   Q_io(I, JP, K, 2) * RNDI1(JJP) - &
                    ( Q_io(I, JP, K, 2) - Q_io(I, JP, KS, 2) ) * 0.50_WP * RNDI1(JJP) ) * DYFI(JJ)
                ELSE
                    DVDL_CCT(2, 2) = ( Q_io(I, JP, K, 2) * RNDI1(JJP) - Q_io(I, J, K, 2) * RNDI1(JJ) ) * DYFI(JJ)
                END IF
                DVDL_CCT(2, 3) =  U_CCT_KP(2) - U_CCT_KM(2)           ! * 0.5_WP * 0.5_WP * DZI

                !======================== DV/ DX.DY.DZ ======================================
                DVDL_CCT(3, 1) =  U_CCT_IP(3) - U_CCT_IM(3)           ! * 0.5_WP * 0.5_WP * DXI
                DVDL_CCT(3, 2) = ( U_CCT_JP(3) - U_CCT(3) ) * DYCI(JJP) + &
                ( U_CCT(3) - U_CCT_JM(3) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                DVDL_CCT(3, 3) =  Q_io(I, J, KP, 3) * RCCI1(JJ) - Q_io(I, J, K, 3) * RCCI1(JJ)      ! * DZI


                DO M = 1, NDV
                    DO N = 1, NDV
                        DVDL1zL_io(I, J, M, N) = DVDL1zL_io(I, J, M, N) + DVDL_CCT(M, N)
                        DVDLPzL_io(I, J, M, N) = DVDLPzL_io(I, J, M, N) + DVDL_CCT(M, N) * PR_io(I, J, K)
                    END DO
                END DO

                !                     DO M = 1, NDV
                !                        DO N = 1, NDV
                !                            DO H = 1, NDV
                !                                DO P = 1, NDV
                !                                    L1 = (M - 1) * 3 + H
                !                                    L2 = (N - 1) * 3 + P
                !                                    DVDL2zL_io(I, J, L1, L2) = DVDL2zL_io(I, J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)
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
                                DVDL2zL_io(I, J, L1, L2) =  DVDL2zL_io(I, J, L1, L2) + DVDL_CCT(M,H) * DVDL_CCT(N,P)
                            END DO
                        END DO
                    END DO
                END DO

            END DO

            U1zL_io(I, J, NDV + 1) = U1zL_io(I, J, NDV + 1) * COE0
            U1zL_io(I, J, 1 : NDV) = U1zL_io(I, J, 1 : NDV) * COE1
            G1zL_io(I, J, 1 : NDV) = G1zL_io(I, J, 1 : NDV) * COE1
            UPzL_io(I, J, :) = UPzL_io(I, J, :) * COE1

            U2zL_io(I, J, :) = U2zL_io(I, J, :) * COE2
            UGzL_io(I, J, :) = UGzL_io(I, J, :) * COE2

            UGUzL_io(I, J, :) = UGUzL_io(I, J, :) * COE3

            DO M = 1, NDV
                DO N = 1, NDV
                    DVDL1zL_io(I, J, M, N) = DVDL1zL_io(I, J, M, N) * COG1(M, N)
                    DVDLPzL_io(I, J, M, N) = DVDLPzL_io(I, J, M, N) * COG1(M, N)
                END DO
            END DO

            !                DO M = 1, NDV
            !                    DO N = 1, NDV
            !                        DO H = 1, NDV
            !                            DO P = 1, NDV
            !                                L1 = (M - 1) * 3 + H
            !                                L2 = (N - 1) * 3 + P
            !                                DVDL2zL_io(I, J, L1, L2) = DVDL2zL_io(I, J, L1, L2) * COG2(L1, L2)
            !                            END DO
            !                        END DO
            !                    END DO
            !                END DO

            DO M = 1, NDV
                DO H = 1, NDV
                    L1 = (M - 1) * NDV + H
                    DO N = 1, NDV
                        DO P = 1, NDV
                            L2 = (N - 1) * NDV + P
                            DVDL2zL_io(I, J, L1, L2) = DVDL2zL_io(I, J, L1, L2) * COG2(L1, L2)
                        END DO
                    END DO
                END DO
            END DO


        END DO
    END DO

    G1rate_io = 0.0_WP
    G1maxx_io = 0.0_WP
    DO J = 1, N2DO(MYID)  !@
        JJ = JCL2G(J)
        DO I = 1, NCL1_io
            G1rate_io =G1rate_io + G1zL_io(I, J, 1) / DYFI(JJ) / RCCI1(JJ)
            G1maxx_io = DMAX1(G1maxx_io, DABS(G1zL_io(I, J, 1)))
        END DO
    END DO
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(G1rate_io, G1rate_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(G1maxx_io, G1maxx_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    G1BULK_WORK_io = G1rate_WORK_io * DBLE(NCL3) / VOLM_io


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_MEAN_Z_THEML_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE postprocess_info
    USE thermal_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IP, KP, JP, N, JJ
    REAL(WP) :: COE0, COE1
    REAL(WP) :: U_CCT(NDV), G_CCT(NDV)
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_DHT


    COE0 = 1.0_WP / DBLE(NCL3)
    COE1 = 1.0_WP / DBLE(NCL3) * 0.5_WP

    DO J = 1, N2DO(MYID)
        JP = JLPV(J)
        JJ = JCL2G(J)
        DO I = 1, NCL1_io
            IP = IPV_io(I)

            T1zL_io(I, J) = 0.0_WP
            D1zL_io(I, J) = 0.0_WP
            H1zL_io(I, J) = 0.0_WP
            K1zL_io(I, J) = 0.0_WP
            M1zL_io(I, J) = 0.0_WP

            T2zL_io(I, J) = 0.0_WP
            D2zL_io(I, J) = 0.0_WP
            H2zL_io(I, J) = 0.0_WP
            K2zL_io(I, J) = 0.0_WP
            M2zL_io(I, J) = 0.0_WP

            DHzL_io(I, J) = 0.0_WP

            UHzL_io(I, J, :) = 0.0_WP
            GHzL_io(I, J, :) = 0.0_WP

            DO K = 1, NCL3
                KP = KPV(K)

                T1zL_io(I, J) = T1zL_io(I, J) + TEMPERATURE(I, J, K)
                D1zL_io(I, J) = D1zL_io(I, J) + DENSITY(I, J, K)
                H1zL_io(I, J) = H1zL_io(I, J) + ENTHALPY(I, J, K)
                K1zL_io(I, J) = K1zL_io(I, J) + THERMCONDT(I, J, K)
                M1zL_io(I, J) = M1zL_io(I, J) + Viscousity(I, J, K)

                T2zL_io(I, J) = T2zL_io(I, J) + TEMPERATURE(I, J, K) * TEMPERATURE(I, J, K)
                D2zL_io(I, J) = D2zL_io(I, J) + DENSITY(I, J, K)   * DENSITY(I, J, K)
                H2zL_io(I, J) = H2zL_io(I, J) + ENTHALPY(I, J, K) * ENTHALPY(I, J, K)
                K2zL_io(I, J) = K2zL_io(I, J) + THERMCONDT(I, J, K) * THERMCONDT(I, J, K)
                M2zL_io(I, J) = M2zL_io(I, J) + Viscousity(I, J, K) * Viscousity(I, J, K)

                DHzL_io(I, J) = DHzL_io(I, J) + ENTHALPY(I, J, K) * DENSITY(I, J, K)

                U_CCT(1) = Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)              ! U at i, j, k
                U_CCT(2) = (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * RCCI1(JJ)   ! V at i, j, k
                U_CCT(3) = (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * RCCI1(JJ)   ! W at i, j, k

                G_CCT(1) = G_io(I, J, K, 1) + G_io(IP, J, K, 1)             ! RHO * U at i, j, k
                G_CCT(2) = (G_io(I, J, K, 2) + G_io(I, JP, K, 2)) * RCCI1(JJ)  ! RHO * U at i, j, k
                G_CCT(3) = (G_io(I, J, K, 3) + G_io(I, J, KP, 3)) * RCCI1(JJ)  ! RHO * U at i, j, k


                DO N = 1, NDV
                    UHzL_io(I, J, N) = UHzL_io(I, J, N) + ENTHALPY(I, J, K) * U_CCT(N)
                    GHzL_io(I, J, N) = GHzL_io(I, J, N) + ENTHALPY(I, J, K) * G_CCT(N)
                END DO


            END DO

            T1zL_io(I, J) = T1zL_io(I, J) * COE0
            D1zL_io(I, J) = D1zL_io(I, J) * COE0
            H1zL_io(I, J) = H1zL_io(I, J) * COE0
            K1zL_io(I, J) = K1zL_io(I, J) * COE0
            M1zL_io(I, J) = M1zL_io(I, J) * COE0

            T2zL_io(I, J) = T2zL_io(I, J) * COE0
            D2zL_io(I, J) = D2zL_io(I, J) * COE0
            H2zL_io(I, J) = H2zL_io(I, J) * COE0
            K2zL_io(I, J) = K2zL_io(I, J) * COE0
            M2zL_io(I, J) = M2zL_io(I, J) * COE0

            DHzL_io(I, J) = DHzL_io(I, J) * COE0

            UHzL_io(I, J, :) = UHzL_io(I, J, :) * COE1
            GHzL_io(I, J, :) = GHzL_io(I, J, :) * COE1

        END DO
    END DO


    DH1rate_io = 0.0_WP
    T1maxx_io = 0.0_WP
    DO J = 1, N2DO(MYID)  !@
        JJ = JCL2G(J)
        DO I = 1, NCL1_io
            DH1rate_io = DH1rate_io + DHzL_io(I, J) / DYFI(JJ) / RCCI1(JJ)
            T1maxx_io = DMAX1(T1maxx_io, DABS(T1zL_io(I, J)))
        END DO
    END DO
    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(DH1rate_io, DH1rate_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(T1maxx_io, T1maxx_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

    !H1bulk_WORK_io = H1RATE_WORK_iO /G1rate_WORK_io
    T1bulk_WORK_io = spline_interpolation_DHT(DH1bulk_WORK_io)

    RETURN
END SUBROUTINE
