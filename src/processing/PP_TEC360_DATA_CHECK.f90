!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE:
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
!============================= TG+ IO==================================
MODULE TEC360_INFO
    USE init_info
    USE flow_info
    USE postprocess_info
    USE mesh_info
    USE WRT_INFO
    USE thermal_info
    !======================== tg =========================
    REAL(WP), ALLOCATABLE :: U_F0_tg(:, :, :, :)
    REAL(WP), ALLOCATABLE :: U1xzL_F0_tg(:, :)
    REAL(WP), ALLOCATABLE :: U1xzL_INTP_tg(:, :)

    REAL(WP), ALLOCATABLE :: U_INTP_tg(:, :, :, :)
    REAL(WP), ALLOCATABLE :: Uprime_tg(:, :, :, :)

    REAL(WP), ALLOCATABLE :: Qcr_tg(:, :, :)
    REAL(WP), ALLOCATABLE :: Vor_tg(:, :, :, :)
    REAL(WP), ALLOCATABLE :: Delta_tg(:, :, :)
    REAL(WP), ALLOCATABLE :: Lambda2_tg(:, :, :)
    REAL(WP), ALLOCATABLE :: SWIrlStrength_tg(:, :, :, :)

    !======================== iO=========================
    REAL(WP), ALLOCATABLE :: U_F0_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: G_F0_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: H_F0_io(:, :, :)
    REAL(WP), ALLOCATABLE :: T_F0_io(:, :, :)
    REAL(WP), ALLOCATABLE :: D_F0_io(:, :, :)
    REAL(WP), ALLOCATABLE :: M_F0_io(:, :, :)


    REAL(WP), ALLOCATABLE :: U1xzL_F0_io(:, :)
    REAL(WP), ALLOCATABLE :: G1xzL_F0_io(:, :)
    REAL(WP), ALLOCATABLE :: H1xzL_F0_io(:)
    REAL(WP), ALLOCATABLE :: T1xzL_F0_io(:)
    REAL(WP), ALLOCATABLE :: D1xzL_F0_io(:)
    REAL(WP), ALLOCATABLE :: M1xzL_F0_io(:)

    REAL(WP), ALLOCATABLE :: U1zL_F0_io(:, :, :)
    REAL(WP), ALLOCATABLE :: G1zL_F0_io(:, :, :)
    REAL(WP), ALLOCATABLE :: H1zL_F0_io(:, :)
    REAL(WP), ALLOCATABLE :: T1zL_F0_io(:, :)
    REAL(WP), ALLOCATABLE :: D1zL_F0_io(:, :)
    REAL(WP), ALLOCATABLE :: M1zL_F0_io(:, :)


    REAL(WP), ALLOCATABLE :: U_INTP_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: G_INTP_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: H_INTP_io(:, :, :)
    REAL(WP), ALLOCATABLE :: T_INTP_io(:, :, :)
    REAL(WP), ALLOCATABLE :: D_INTP_io(:, :, :)
    REAL(WP), ALLOCATABLE :: M_INTP_io(:, :, :)


    REAL(WP), ALLOCATABLE :: U1xzL_INTP_io(:, :)
    REAL(WP), ALLOCATABLE :: G1xzL_INTP_io(:, :)
    REAL(WP), ALLOCATABLE :: H1xzL_INTP_io(:)
    REAL(WP), ALLOCATABLE :: T1xzL_INTP_io(:)
    REAL(WP), ALLOCATABLE :: D1xzL_INTP_io(:)
    REAL(WP), ALLOCATABLE :: M1xzL_INTP_io(:)

    REAL(WP), ALLOCATABLE :: U1zL_INTP_io(:, :, :)
    REAL(WP), ALLOCATABLE :: G1zL_INTP_io(:, :, :)
    REAL(WP), ALLOCATABLE :: H1zL_INTP_io(:, :)
    REAL(WP), ALLOCATABLE :: T1zL_INTP_io(:, :)
    REAL(WP), ALLOCATABLE :: D1zL_INTP_io(:, :)
    REAL(WP), ALLOCATABLE :: M1zL_INTP_io(:, :)

    REAL(WP), ALLOCATABLE :: Uprime_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: Gprime_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: Tprime_io(:, :, :)
    REAL(WP), ALLOCATABLE :: Dprime_io(:, :, :)
    REAL(WP), ALLOCATABLE :: Mprime_io(:, :, :)
    REAL(WP), ALLOCATABLE :: UDprime_io(:, :, :, :)

    REAL(WP), ALLOCATABLE :: Qcr_io(:, :, :)
    REAL(WP), ALLOCATABLE :: Vor_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: Delta_io(:, :, :)
    REAL(WP), ALLOCATABLE :: Lambda2_io(:, :, :)
    REAL(WP), ALLOCATABLE :: SWIrlStrength_io(:, :, :, :)

    INTEGER(4) :: NI
    !INTEGER(4) :: NJ
    INTEGER(4) :: NK
    INTEGER(4), ALLOCATABLE :: IID(:)
    !INTEGER(4), ALLOCATABLE :: JID(:)
    INTEGER(4), ALLOCATABLE :: KID(:)
    INTEGER(4) :: NCOUNT(5) = 0

END MODULE TEC360_INFO
!**********************************************************************************************************************************
SUBROUTINE DEALLO_TEC
    USE TEC360_INFO
    IMPLICIT NONE

    !======================== tg =========================
    IF ( ALLOCATED(U_F0_tg) )        DEALLOCATE(U_F0_tg)
    IF ( ALLOCATED(U1xzL_F0_tg) )    DEALLOCATE(U1xzL_F0_tg)
    IF ( ALLOCATED(U1xzL_INTP_tg) )  DEALLOCATE(U1xzL_INTP_tg)


    IF ( ALLOCATED(Uprime_tg) )  DEALLOCATE(Uprime_tg)
    IF ( ALLOCATED(U_INTP_tg) )  DEALLOCATE(U_INTP_tg)

    IF ( ALLOCATED(Qcr_tg) )   DEALLOCATE(Qcr_tg)
    IF ( ALLOCATED(Vor_tg) )   DEALLOCATE(Vor_tg)
    IF ( ALLOCATED(Delta_tg) ) DEALLOCATE(Delta_tg)
    IF ( ALLOCATED(Lambda2_tg) ) DEALLOCATE(Lambda2_tg)
    IF ( ALLOCATED(SWIrlStrength_tg) )    DEALLOCATE(SWIrlStrength_tg)


    !======================== iO=========================
    IF ( ALLOCATED(U_F0_io) ) DEALLOCATE(U_F0_io)
    IF ( ALLOCATED(G_F0_io) ) DEALLOCATE(G_F0_io)
    IF ( ALLOCATED(H_F0_io) ) DEALLOCATE(H_F0_io)
    IF ( ALLOCATED(T_F0_io) ) DEALLOCATE(T_F0_io)
    IF ( ALLOCATED(D_F0_io) ) DEALLOCATE(D_F0_io)
    IF ( ALLOCATED(M_F0_io) ) DEALLOCATE(M_F0_io)



    IF ( ALLOCATED(U1xzL_F0_io) ) DEALLOCATE(U1xzL_F0_io)
    IF ( ALLOCATED(G1xzL_F0_io) ) DEALLOCATE(G1xzL_F0_io)
    IF ( ALLOCATED(T1xzL_F0_io) ) DEALLOCATE(H1xzL_F0_io)
    IF ( ALLOCATED(T1xzL_F0_io) ) DEALLOCATE(T1xzL_F0_io)
    IF ( ALLOCATED(D1xzL_F0_io) ) DEALLOCATE(D1xzL_F0_io)
    IF ( ALLOCATED(M1xzL_F0_io) ) DEALLOCATE(M1xzL_F0_io)

    IF ( ALLOCATED(U1zL_F0_io) )  DEALLOCATE(U1zL_F0_io)
    IF ( ALLOCATED(G1zL_F0_io) )  DEALLOCATE(G1zL_F0_io)
    IF ( ALLOCATED(T1zL_F0_io) )  DEALLOCATE(H1zL_F0_io)
    IF ( ALLOCATED(T1zL_F0_io) )  DEALLOCATE(T1zL_F0_io)
    IF ( ALLOCATED(D1zL_F0_io) )  DEALLOCATE(D1zL_F0_io)
    IF ( ALLOCATED(M1zL_F0_io) )  DEALLOCATE(M1zL_F0_io)


    IF ( ALLOCATED(U_INTP_io) ) DEALLOCATE(U_INTP_io)
    IF ( ALLOCATED(G_INTP_io) ) DEALLOCATE(G_INTP_io)
    IF ( ALLOCATED(H_INTP_io) ) DEALLOCATE(H_INTP_io)
    IF ( ALLOCATED(T_INTP_io) ) DEALLOCATE(T_INTP_io)
    IF ( ALLOCATED(D_INTP_io) ) DEALLOCATE(D_INTP_io)
    IF ( ALLOCATED(M_INTP_io) ) DEALLOCATE(M_INTP_io)


    IF ( ALLOCATED(U1xzL_INTP_io) ) DEALLOCATE(U1xzL_INTP_io)
    IF ( ALLOCATED(G1xzL_INTP_io) ) DEALLOCATE(G1xzL_INTP_io)
    IF ( ALLOCATED(T1xzL_INTP_io) ) DEALLOCATE(H1xzL_INTP_io)
    IF ( ALLOCATED(T1xzL_INTP_io) ) DEALLOCATE(T1xzL_INTP_io)
    IF ( ALLOCATED(D1xzL_INTP_io) ) DEALLOCATE(D1xzL_INTP_io)
    IF ( ALLOCATED(M1xzL_INTP_io) ) DEALLOCATE(M1xzL_INTP_io)

    IF ( ALLOCATED(U1zL_INTP_io) )  DEALLOCATE(U1zL_INTP_io)
    IF ( ALLOCATED(G1zL_INTP_io) )  DEALLOCATE(G1zL_INTP_io)
    IF ( ALLOCATED(T1zL_INTP_io) )  DEALLOCATE(H1zL_INTP_io)
    IF ( ALLOCATED(T1zL_INTP_io) )  DEALLOCATE(T1zL_INTP_io)
    IF ( ALLOCATED(D1zL_INTP_io) )  DEALLOCATE(D1zL_INTP_io)
    IF ( ALLOCATED(M1zL_INTP_io) )  DEALLOCATE(M1zL_INTP_io)

    IF ( ALLOCATED(Uprime_io) )   DEALLOCATE(Uprime_io)
    IF ( ALLOCATED(Gprime_io) )   DEALLOCATE(Gprime_io)
    IF ( ALLOCATED(Tprime_io) )   DEALLOCATE(Tprime_io)
    IF ( ALLOCATED(Dprime_io) )   DEALLOCATE(Dprime_io)
    IF ( ALLOCATED(Mprime_io) )   DEALLOCATE(Mprime_io)
    IF ( ALLOCATED(UDprime_io) )  DEALLOCATE(UDprime_io)

    IF ( ALLOCATED(Qcr_io) )        DEALLOCATE(Qcr_io)
    IF ( ALLOCATED(Vor_io) )        DEALLOCATE(Vor_io)
    IF ( ALLOCATED(Delta_io) )      DEALLOCATE(Delta_io)
    IF ( ALLOCATED(Lambda2_io) )    DEALLOCATE(Lambda2_io)
    IF ( ALLOCATED(SWIrlStrength_io) )    DEALLOCATE(SWIrlStrength_io)

    IF ( ALLOCATED(IID) ) DEALLOCATE(IID)
    !IF ( ALLOCATED(JID) ) DEALLOCATE(JID)
    IF ( ALLOCATED(KID) ) DEALLOCATE(KID)


    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
!********************** WRITE INITIAL CALCUATED RESULTS OUT*******************************************
SUBROUTINE CALL_TEC360
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: N

    IF(iPostProcess == 1) THEN
        IF(TgFlowFlg) THEN
            !!IF(IoFlowFlg) add ... more....
        ELSE
            CALL PP_MEAN_ZX_FLOW_Xperiodic_io
        END IF

    END IF

    !============WRITE INITIAL CALCUATED RESULTS OUT ======================================
    IF(TgFlowFlg .AND. IoFlowFlg) THEN

        CALL GATHERING_TG
        CALL GATHERING_io

        IF(MYID == 0) THEN
            NI = 4
            !NJ = 4
            NK = 2
            ALLOCATE( IID(NI) )
            !ALLOCATE( JID(NJ) )
            ALLOCATE( KID(NK) )

            CALL TEC360_EXTRPLAT2ND_MASTER_TG
            CALL TEC360_EXTRPLAT2ND_MASTER_io
        END IF

    ELSE
        IF(MYID == 0) THEN
            NI = 1
            !NJ = 8
            NK = 1
            ALLOCATE( IID(NI) )
            !ALLOCATE( JID(NJ) )
            ALLOCATE( KID(NK) )
        END IF

        IF (TgFlowFlg) THEN
            CALL GATHERING_TG
            IF(MYID == 0) CALL TEC360_EXTRPLAT2ND_MASTER_TG
        END IF

        IF(IoFlowFlg) THEN
            CALL GATHERING_io
            IF(MYID == 0) CALL TEC360_EXTRPLAT2ND_MASTER_io
        END IF
    END IF

    IF(MYID == 0) THEN
        CALL TEC360_INSTANT_VORTEX_CRITERIA
        CALL TEC360_INSTANT_Uprime

        IF(iPostProcess == 1 .AND. iPPInst == 1) THEN
            CALL TEC360_ALL_NODES
        END IF

        !===================WRITE X -SLICE ==================================

        DO N = 1, NI
            IID(N) = N * (NND1_io + 1) / (2 * NI) + 1
            CALL TEC360_XSLICE(N)
        END DO

        !===================WRITE Y-SLICE ==================================
        !DO N = 1, NJ
        !JID(N) = N * (NND2 + 1) / (2 * NJ) + 1
        DO N = 1, MGRID
            CALL TEC360_YSLICE(N)
        END DO
        !===================WRITE Z -SLICE ==================================
        DO N = 1, NK
            KID(N) = N * (NND3 + 1) / (2 * NK) + 1
            CALL TEC360_ZSLICE(N)
        END DO
    END IF

    CALL DEALLO_TEC

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE GATHERING_TG
    USE TEC360_INFO
    IMPLICIT NONE
    REAL(WP), ALLOCATABLE :: VARLoc1(:, :, :, :)
    REAL(WP), ALLOCATABLE :: VARGlb1(:, :, :, :)
    REAL(WP), ALLOCATABLE :: VARLoc2(:, :)
    REAL(WP), ALLOCATABLE :: VARGlb2(:, :)
    INTEGER(4) :: N, I, JJ, K, KS

    !=======================global =============================
    ALLOCATE ( U_F0_tg     (NCL1_tg, NCL2, NCL3, NDV + 1) )
    ALLOCATE ( U1xzL_F0_tg (NCL2, NDV + 1 ) )

    !=========================gather U,P ==========================
    N = 4
    ALLOCATE ( VARLoc1(NCL1_tg, 0:(N2DO(0) + 1), NCL3, N))
    ALLOCATE ( VARGlb1(NCL1_tg, NCL2,         NCL3, N))

    VARLoc1(1 : NCL1_tg, 0:(N2DO(0) + 1), 1 : NCL3, 1:3) = Q_tg (1 : NCL1_tg, 0:(N2DO(0) + 1), 1 : NCL3, 1:3)
    VARLoc1(1 : NCL1_tg, 0:(N2DO(0) + 1), 1 : NCL3, 4) = PR_tg(1 : NCL1_tg, 0:(N2DO(0) + 1), 1 : NCL3    )

    CALL GATHERING_xyzn(VARLoc1, VARGlb1, 1, NCL1_tg, N)

    IF(MYID == 0) THEN
        DO I = 1, NCL1_tg
            DO K = 1, NCL3
                DO JJ = 1, NCL2
                    U_F0_tg(I, JJ, K, 1) = VARGlb1(I, JJ, K, 1)
                    U_F0_tg(I, JJ, K, 3) = VARGlb1(I, JJ, K, 3) * RCCI1(JJ)
                    U_F0_tg(I, JJ, K, 4) = VARGlb1(I, JJ, K, 4)
                    IF(JJ == 1 .AND. iCase == ipipec) THEN
                        KS = KSYM(K)
                        U_F0_tg(I, JJ, K, 2) = (VARGlb1(I, 2, K, 2) - VARGlb1(I, 2, KS, 2) ) * 0.50_WP * RNDI1(2)
                    ELSE
                        U_F0_tg(I, JJ, K, 2) = VARGlb1(I, JJ, K, 2) * RNDI1(JJ)
                    END IF
                END DO
            END DO
        END DO
    END IF
    DEALLOCATE (VARLoc1)
    DEALLOCATE (VARGlb1)

    !=========================gather averaged==========================
    N = 4
    ALLOCATE ( VARLoc2(1 : N2DO(0), N))
    ALLOCATE ( VARGlb2(NCL2,     N))

    VARLoc2(1 : N2DO(0), 1:4) = U1xzL_tg (1 : N2DO(0), 1:4)   ! should space and time averaged? or just space averaged?

    CALL GATHERING_yn(VARLoc2, VARGlb2, N)

    IF(MYID == 0) THEN
        U1xzL_F0_tg(1 : NCL2, 1:4) = VARGlb2(1 : NCL2, 1:4)
    END IF

    DEALLOCATE (VARLoc2)
    DEALLOCATE (VARGlb2)


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE GATHERING_io
    USE TEC360_INFO
    IMPLICIT NONE
    REAL(WP), ALLOCATABLE :: VARLoc1(:, :, :, :)
    REAL(WP), ALLOCATABLE :: VARGlb1(:, :, :, :)
    REAL(WP), ALLOCATABLE :: VARLoc2(:, :)
    REAL(WP), ALLOCATABLE :: VARGlb2(:, :)
    REAL(WP), ALLOCATABLE :: VARLoc3(:, :, :)
    REAL(WP), ALLOCATABLE :: VARGlb3(:, :, :)
    INTEGER(4) :: N, I, K, JJ, KS

    !======================globaL ===================================
    ALLOCATE ( U_F0_io(NCL1S : NCL1E, NCL2, NCL3, NDV + 1) )

    IF(iThermoDynamics == 1) THEN
        ALLOCATE ( G_F0_io(NCL1S : NCL1E, NCL2, NCL3, NDV) )
        ALLOCATE ( T_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
        ALLOCATE ( D_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
        ALLOCATE ( M_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
        ALLOCATE ( H_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
    END IF

    IF(TgFlowFlg) THEN
        ALLOCATE ( U1zL_F0_io(NCL1S : NCL1E, NCL2, NDV + 1 ) )
        IF(iThermoDynamics == 1) THEN
            ALLOCATE ( G1zL_F0_io(NCL1S : NCL1E, NCL2, NDV ) )
            ALLOCATE ( H1zL_F0_io(NCL1S : NCL1E, NCL2      ) )
            ALLOCATE ( T1zL_F0_io(NCL1S : NCL1E, NCL2      ) )
            ALLOCATE ( D1zL_F0_io(NCL1S : NCL1E, NCL2      ) )
            ALLOCATE ( M1zL_F0_io(NCL1S : NCL1E, NCL2      ) )

        END IF
    ELSE
        ALLOCATE ( U1xzL_F0_io(NCL2, NDV + 1 ) )
        IF(iThermoDynamics == 1) THEN
            ALLOCATE ( G1xzL_F0_io(NCL2, NDV ) )
            ALLOCATE ( H1xzL_F0_io(NCL2      ) )
            ALLOCATE ( T1xzL_F0_io(NCL2      ) )
            ALLOCATE ( D1xzL_F0_io(NCL2      ) )
            ALLOCATE ( M1xzL_F0_io(NCL2      ) )
        END IF
    END IF

    !=========================gather U,P ==========================
    IF(iThermoDynamics == 1) THEN
        N = 11
    ELSE
        N = 4
    END IF
    ALLOCATE ( VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), NCL3, N))
    ALLOCATE ( VARGlb1(NCL1S : NCL1E, NCL2,         NCL3, N))

    VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 1:3) = Q_io(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 1:3)
    VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 4) = PR_io(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3    )
    IF(iThermoDynamics == 1) THEN
        VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 5:7) = G_io       (NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 1:3)
        VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3,8) = ENTHALPY   (NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3)
        VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3,9) = TEMPERATURE(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3)
        VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 10) = DENSITY    (NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3)
        VARLoc1(NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3, 11) = Viscousity (NCL1S : NCL1E, 0:(N2DO(0) + 1), 1 : NCL3)
    END IF

    CALL GATHERING_xyzn(VARLoc1, VARGlb1, NCL1S, NCL1E, N)

    IF(MYID == 0) THEN

        DO I = NCL1S, NCL1E
            DO K = 1, NCL3
                DO JJ = 1, NCL2
                    U_F0_io(I, JJ, K, 1) = VARGlb1(I, JJ, K, 1)
                    U_F0_io(I, JJ, K, 3) = VARGlb1(I, JJ, K, 3) * RCCI1(JJ)
                    U_F0_io(I, JJ, K, 4) = VARGlb1(I, JJ, K, 4)
                    IF(JJ == 1 .AND. iCase == ipipec) THEN
                        KS = KSYM(K)
                        U_F0_io(I, JJ, K, 2) = (VARGlb1(I, 2, K, 2) - VARGlb1(I, 2, KS, 2) ) * 0.50_WP * RNDI1(2)
                    ELSE
                        U_F0_io(I, JJ, K, 2) = VARGlb1(I, JJ, K, 2) * RNDI1(JJ)
                    END IF
                END DO
            END DO
        END DO

        IF(iThermoDynamics == 1) THEN
            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    DO JJ = 1, NCL2
                        G_F0_io(I, JJ, K, 1) = VARGlb1(I, JJ, K, 5)
                        G_F0_io(I, JJ, K, 3) = VARGlb1(I, JJ, K, 7) * RCCI1(JJ)
                        IF(JJ == 1 .AND. iCase == ipipec) THEN
                            KS = KSYM(K)
                            G_F0_io(I, JJ, K, 2) = (VARGlb1(I, 2, K,6) - VARGlb1(I, 2, KS,6) ) * 0.50_WP * RNDI1(2)
                        ELSE
                            G_F0_io(I, JJ, K, 2) = VARGlb1(I, JJ, K, 6) * RNDI1(JJ)
                        END IF
                    END DO
                END DO
            END DO
            H_F0_io(NCL1S : NCL1E, 1 : NCL2, 1 : NCL3) = VARGlb1(NCL1S : NCL1E, 1 : NCL2,  1 : NCL3, 8)
            T_F0_io(NCL1S : NCL1E, 1 : NCL2, 1 : NCL3) = VARGlb1(NCL1S : NCL1E, 1 : NCL2,  1 : NCL3, 9)
            D_F0_io(NCL1S : NCL1E, 1 : NCL2, 1 : NCL3) = VARGlb1(NCL1S : NCL1E, 1 : NCL2,  1 : NCL3, 10)
            M_F0_io(NCL1S : NCL1E, 1 : NCL2, 1 : NCL3) = VARGlb1(NCL1S : NCL1E, 1 : NCL2,  1 : NCL3, 11)
        END IF
    END IF

    DEALLOCATE (VARLoc1)
    DEALLOCATE (VARGlb1)

    !=========================gather averaged==========================
    IF(TgFlowFlg) THEN

        N = 4
        ALLOCATE ( VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), N))
        ALLOCATE ( VARGlb3(NCL1S : NCL1E, NCL2,     N))

        VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), 1:4) = U1zL_io(NCL1S : NCL1E, 1 : N2DO(0), 1:4)

        CALL GATHERING_xyn(VARLoc3, VARGlb3, NCL1S, NCL1E, N)

        IF(MYID == 0) THEN
            U1zL_F0_io(NCL1S : NCL1E, 1 : NCL2, 1:4) = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 1:4)
        END IF
        DEALLOCATE (VARLoc3)
        DEALLOCATE (VARGlb3)

        IF(iThermoDynamics == 1) THEN
            N = 7

            ALLOCATE ( VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), N))
            ALLOCATE ( VARGlb3(NCL1S : NCL1E, NCL2,     N))

            VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), 1:3) = G1zL_io(NCL1S : NCL1E, 1 : N2DO(0), 1:3)
            VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), 4) = H1zL_io(NCL1S : NCL1E, 1 : N2DO(0)    )
            VARLoc3(NCL1S : NCL1E, 1 : N2DO(0), 5) = T1zL_io(NCL1S : NCL1E, 1 : N2DO(0)    )
            VARLoc3(NCL1S : NCL1E, 1 : N2DO(0),6) = D1zL_io(NCL1S : NCL1E, 1 : N2DO(0)    )
            VARLoc3(NCL1S : NCL1E, 1 : N2DO(0),7) = M1zL_io(NCL1S : NCL1E, 1 : N2DO(0)    )

            CALL GATHERING_xyn(VARLoc3, VARGlb3, NCL1S, NCL1E, N)

            IF(MYID == 0) THEN
                G1zL_F0_io(NCL1S : NCL1E, 1 : NCL2, 1:3) = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 1:3)
                H1zL_F0_io(NCL1S : NCL1E, 1 : NCL2)    = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 4  )
                T1zL_F0_io(NCL1S : NCL1E, 1 : NCL2)    = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 5  )
                D1zL_F0_io(NCL1S : NCL1E, 1 : NCL2)    = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 6  )
                M1zL_F0_io(NCL1S : NCL1E, 1 : NCL2)    = VARGlb3(NCL1S : NCL1E, 1 : NCL2, 7  )
            END IF

            DEALLOCATE (VARLoc3)
            DEALLOCATE (VARGlb3)
        END IF
    ELSE
        N = 4
        ALLOCATE ( VARLoc2(1 : N2DO(0), N))
        ALLOCATE ( VARGlb2(NCL2,     N))

        VARLoc2(1 : N2DO(0), 1:4) = U1xzL_io(1 : N2DO(0), 1:4)

        CALL GATHERING_yn(VARLoc2, VARGlb2, N)

        IF(MYID == 0) THEN
            U1xzL_F0_io(1 : NCL2, 1:4) = VARGlb2(1 : NCL2, 1:4)
        END IF

        DEALLOCATE (VARLoc2)
        DEALLOCATE (VARGlb2)

        IF(iThermoDynamics == 1) THEN
            N = 7

            ALLOCATE ( VARLoc2(1 : N2DO(0), N))
            ALLOCATE ( VARGlb2(NCL2,     N))

            VARLoc2(1 : N2DO(0), 1:3) = G1xzL_io(1 : N2DO(0), 1:3)
            VARLoc2(1 : N2DO(0), 4) = H1xzL_io(1 : N2DO(0)    )
            VARLoc2(1 : N2DO(0), 5) = T1xzL_io(1 : N2DO(0)    )
            VARLoc2(1 : N2DO(0),6) = D1xzL_io(1 : N2DO(0)    )
            VARLoc2(1 : N2DO(0),7) = M1xzL_io(1 : N2DO(0)    )


            CALL GATHERING_yn(VARLoc2, VARGlb2, N)

            IF(MYID == 0) THEN
                G1xzL_F0_io(1 : NCL2, 1:3) = VARGlb2(1 : NCL2, 1:3)
                H1xzL_F0_io(1 : NCL2)    = VARGlb2(1 : NCL2, 4  )
                T1xzL_F0_io(1 : NCL2)    = VARGlb2(1 : NCL2, 5  )
                D1xzL_F0_io(1 : NCL2)    = VARGlb2(1 : NCL2, 6  )
                M1xzL_F0_io(1 : NCL2)    = VARGlb2(1 : NCL2, 7 )
            END IF

            DEALLOCATE (VARLoc2)
            DEALLOCATE (VARGlb2)
        END IF

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE GATHERING_xyzn(VARLoc, VARGlb, NSX1, NSX2, NSN)
    USE TEC360_INFO
    IMPLICIT none
    INTEGER(4) :: NSX1, NSX2, NSN
    REAL(WP), INTENT(IN) :: VARLoc(NSX1 : NSX2, 0 : N2DO(0) + 1, NCL3, NSN)
    REAL(WP), INTENT(OUT) :: VARGlb(NSX1 : NSX2, 1 : NCL2,     NCL3, NSN)
    REAL(WP) :: AUX  (NSX1 : NSX2, 0 : N2DO(0) + 1, NCL3, NSN, 1 : NPTOT)
    INTEGER(4) :: INN, N, KK, J, N2DOID, I, JJ, K

    INN = (NSX2 - NSX1 + 1) * (N2DO(0) + 2) * NCL3 * NSN
    CALL MPI_GATHER( VARLoc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !============== RE - ARrange gathered data========================
    IF (MYID == 0) THEN
        DO N = 1, NSN
            DO KK = 0, NPSLV
                N2DOID=JDEWT(KK) - JDSWT(KK) + 1
                DO  J = 1, N2DOID
                    JJ = JDSWT(KK) - 1 + J
                    DO I = NSX1, NSX2
                        DO K = 1, NCL3
                            VARGlb(I, JJ, K, N) = AUX(I, J, K, N, KK+ 1)
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE GATHERING_yn(VARLoc, VARGlb, NSN)
    USE tec360_info
    IMPLICIT none
    INTEGER(4) :: NSN
    REAL(WP), INTENT(IN) :: VARLoc(N2DO(MYID),   NSN)
    REAL(WP), INTENT(OUT) :: VARGlb(NCL2,         NSN)
    REAL(WP) :: AUX  (N2DO(MYID),   NSN, 1 : NPTOT)
    INTEGER(4) :: INN, N, IP, J, N2DOID, I, JJ

    INN = N2DO(0) * NSN
    CALL MPI_GATHER( VARLoc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !============== RE - ARrange gathered data========================
    IF (MYID == 0) THEN

        DO IP = 0, NPSLV
            N2DOID=JDEWT(IP) - JDSWT(IP) + 1
            DO  J = 1, N2DOID
                JJ = JDSWT(IP) - 1 + J

                DO N = 1, nsn
                    VARGlb(JJ, N) = AUX(J, N, IP + 1)
                END DO

            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE GATHERING_xyn(VARLoc, VARGlb, NSX1, NSX2, NSN)
    USE tec360_info
    IMPLICIT none
    INTEGER(4) :: NSX1, NSX2, NSN
    REAL(WP), INTENT(IN) :: VARLoc(NSX1 : NSX2, 1 : N2DO(0), NSN)
    REAL(WP), INTENT(OUT) :: VARGlb(NSX1 : NSX2, 1 : NCL2,    NSN)
    REAL(WP) :: AUX  (NSX1 : NSX2, 1 : N2DO(0), NSN, 1 : NPTOT)
    INTEGER(4) :: INN, N, IP, J, N2DOID, I, JJ

    INN = (NSX2 -NSX1 + 1) * (N2DO(0)) * NSN
    CALL MPI_GATHER( VARLoc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !============== RE - ARrange gathered data========================
    IF (MYID == 0) THEN

        DO IP = 0, NPSLV
            N2DOID=JDEWT(IP) - JDSWT(IP) + 1
            DO  J = 1, N2DOID
                JJ = JDSWT(IP) - 1 + J
                DO I = NSX1, NSX2
                    DO N = 1, nsn
                        VARGlb(I, JJ, N) = AUX(I, J, N, IP + 1)
                    END DO
                END DO
            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE TEC360_EXTRPLAT2ND_MASTER_TG
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IM, JM, KM, KS


    IF (MYID /= 0) RETURN

    ALLOCATE ( U1xzL_INTP_tg(NND2, NDV + 1 ) )

    DO J = 2, NCL2
        JM = J - 1
        U1xzL_INTP_tg(J, 1:4) = (U1xzL_F0_tg(JM, 1:4) + U1xzL_F0_tg(J, 1:4)) * 0.5_WP
    END DO
    U1xzL_INTP_tg(1,   1:3) = 0.0_WP
    U1xzL_INTP_tg(NND2, 1:3) = 0.0_WP
    U1xzL_INTP_tg(1,   4) = U1xzL_F0_tg(1,   4)
    U1xzL_INTP_tg(NND2, 4) = U1xzL_F0_tg(NCL2, 4)


    ALLOCATE ( U_INTP_tg(NND1_tg, NND2, NND3, NDV + 1) )

    !>============= INTERPOLATION ALL VALUES TO POINTS FOR X, Z periodic B.C.==========

    DO I = 1, NCL1_tg
        IM = IMV_tg(I)
        DO K = 1, NCL3
            KM = KMV(K)
            KS = KSYM(K)
            DO J = 2, NCL2
                JM = JGMV(J)
                U_INTP_tg(I, J, K, 1) = 0.250_WP * ( U_F0_tg(I, JM, KM, 1) + &
                U_F0_tg(I, JM, K, 1)  + &
                U_F0_tg(I, J, KM, 1) + &
                U_F0_tg(I, J, K,  1)    )
                U_INTP_tg(I, J, K, 3) = 0.250_WP * ( U_F0_tg(IM, JM, K, 3) + &
                U_F0_tg(IM, J, K, 3)  + &
                U_F0_tg(I, JM, K, 3)  + &
                U_F0_tg(I, J, K,  3) )
                U_INTP_tg(I, J, K, 2) = 0.250_WP * ( U_F0_tg(IM, J, KM, 2) + &
                U_F0_tg(IM, J, K, 2)  + &
                U_F0_tg(I, J, KM, 2)  + &
                U_F0_tg(I, J, K, 2) )
                U_INTP_tg(I, J, K, 4) = 0.1250_WP * (U_F0_tg(IM, JM, KM, 4) + U_F0_tg(IM, JM, K, 4) + &
                U_F0_tg(IM, J, KM, 4) + U_F0_tg(IM, J, K, 4) + &
                U_F0_tg(I, JM, KM, 4) + U_F0_tg(I, JM, K, 4) + &
                U_F0_tg(I, J, KM, 4) + U_F0_tg(I, J, K, 4) )
            END DO
            IF(iCase == iPIPEC) THEN
                U_INTP_tg(I, 1, K, 1:4) = 0.5_WP * ( U_F0_tg(I, 1, K, 1:4) + U_F0_tg(I, 1, KS, 1:4))
            ELSE
                U_INTP_tg(I, 1, K, 1:3) = 0.0_WP
                U_INTP_tg(I, 1, K, 4) = U_F0_tg(I, 1, K, 4)
            END IF

        END DO
    END DO
    U_INTP_tg(NND1_tg, :, :, :) = U_INTP_tg(1, :, :, :)

    U_INTP_tg(:, :, NND3, :) = U_INTP_tg(:, :, 1, :)

    U_INTP_tg(:, NND2, :, 1:3) = 0.0_WP
    U_INTP_tg(:, NND2, :, 4)   = U_INTP_tg(:, NCL2, :, 4)


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TEC360_EXTRPLAT2ND_MASTER_io
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, IM, JM, KM, KS, IE

    IF(.not.IoFlowFlg) RETURN
    IF (MYID /= 0) RETURN


    !=== instantanous variables ==u, V,w,P =======================================
    ALLOCATE ( U_INTP_io(NND1_io, NND2, NND3, NDV + 1) )
    !============ main DOmAIn without b.c. points ==========
    DO I = 1, NCL1_io
        IM = IMV_io(I)
        DO K = 1, NCL3
            KM = KMV(K)
            KS = KSYM(K)
            DO J = 2, NCL2
                JM = JGMV(J)
                U_INTP_io(I, J, K, 1) = 0.250_WP * ( U_F0_io(I, JM, KM, 1)  + &
                U_F0_io(I, JM, K, 1)  + &
                U_F0_io(I, J, KM, 1)  + &
                U_F0_io(I, J, K, 1) )
                U_INTP_io(I, J, K, 2) = 0.250_WP * ( U_F0_io(IM, J, KM, 2) + &
                U_F0_io(IM, J, K, 2)  + &
                U_F0_io(I, J, KM, 2)  + &
                U_F0_io(I, J, K, 2) )
                U_INTP_io(I, J, K, 3) = 0.250_WP * ( U_F0_io(IM, JM, K, 3) + &
                U_F0_io(IM, J, K, 3)  + &
                U_F0_io(I, JM, K, 3)  + &
                U_F0_io(I, J, K, 3) )
                U_INTP_io(I, J, K, 4) = 0.1250_WP * ( U_F0_io(IM, JM, KM, 4) + U_F0_io(IM, JM, K, 4) + &
                U_F0_io(IM, J, KM, 4) + U_F0_io(IM, J, K, 4) + &
                U_F0_io(I, JM, KM, 4) + U_F0_io(I, JM, K, 4) + &
                U_F0_io(I, J, KM, 4) + U_F0_io(I, J, K, 4) )
            END DO

            IF(iCase == iPIPEC) THEN
                U_INTP_io(I, 1, K, 1:4) = 0.5_WP * ( U_F0_io(I, 1, K, 1:4) + U_F0_io(I, 1, KS, 1:4))
            ELSE IF(iCase == IBox3P) THEN
                U_INTP_io(I, 1, K, 1:4) = 0.5_WP * ( U_F0_io(I, 1, K, 1:4) + U_F0_io(I, NCL2, K, 1:4))
            ELSE
                U_INTP_io(I, 1, K, 1:3) = 0.0_WP
                U_INTP_io(I, 1, K, 4) = U_F0_io(I, 1, K, 4)
            END IF
        END DO
    END DO
    !========y b.c. pointS ============
    IF(iCase == IBox3P) THEN
        U_INTP_io(:, NND2, :, 1:3) = U_INTP_io(:, 1, :, 1:3)
    ELSE
        U_INTP_io(:, NND2, :, 1:3) = 0.0_WP
    END IF
    U_INTP_io(:, NND2, :, 4) = U_INTP_io(:, NCL2, :, 4)
    !======== X b.c. pointS ============
    IF(TgFlowFlg) THEN
        U_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3, :) = U_F0_io(NND1_io, 1 : NCL2, 1 : NCL3, :)
    ELSE
        U_INTP_io(NND1_io, :, :, :) = U_INTP_io(1, :, :, :)
    END IF
    !======== Z b.c. pointS ============
    U_INTP_io(:, :, NND3, :) = U_INTP_io(:, :, 1, :)


    !=== instantanous variables \rho u,\rho v,\rho w, h, t, \rho, mu=======================================
    IF(iThermoDynamics == 1) THEN
        ALLOCATE ( G_INTP_io(NND1_io, NND2, NND3, NDV) )
        ALLOCATE ( H_INTP_io(NND1_io, NND2, NND3) )
        ALLOCATE ( T_INTP_io(NND1_io, NND2, NND3) )
        ALLOCATE ( D_INTP_io(NND1_io, NND2, NND3) )
        ALLOCATE ( M_INTP_io(NND1_io, NND2, NND3) )
        !============ main DOmAIn without b.c. points ==========
        DO I = 1, NCL1_io
            IM = IMV_io(I)

            DO K = 1, NCL3
                KM = KMV(K)
                KS = KSYM(K)
                DO J = 2, NCL2
                    JM = JGMV(J)
                    G_INTP_io(I, J, K, 1) = 0.250_WP * ( G_F0_io(I, JM, KM, 1) + &
                    G_F0_io(I, JM, K, 1)  + &
                    G_F0_io(I, J, KM, 1)  + &
                    G_F0_io(I, J, K, 1) )
                    G_INTP_io(I, J, K, 2) = 0.250_WP * ( G_F0_io(IM, J, KM, 2) + &
                    G_F0_io(IM, J, K, 2)  + &
                    G_F0_io(I, J, KM, 2)  + &
                    G_F0_io(I, J, K, 2) )
                    G_INTP_io(I, J, K, 3) = 0.250_WP * ( G_F0_io(IM, JM, K, 3) + &
                    G_F0_io(IM, J, K, 3)  + &
                    G_F0_io(I, JM, K, 3)  + &
                    G_F0_io(I, J, K, 3) )

                    H_INTP_io(I, J, K) = 0.1250_WP * ( H_F0_io(IM, JM, KM) + H_F0_io(IM, JM, K) + &
                    H_F0_io(IM, J, KM) + H_F0_io(IM, J, K) + &
                    H_F0_io(I, JM, KM) + H_F0_io(I, JM, K) + &
                    H_F0_io(I, J, KM) + H_F0_io(I, J, K) )

                    T_INTP_io(I, J, K) = 0.1250_WP *  ( T_F0_io(IM, JM, KM) + T_F0_io(IM, JM, K) + &
                    T_F0_io(IM, J, KM) + T_F0_io(IM, J, K) + &
                    T_F0_io(I, JM, KM) + T_F0_io(I, JM, K) + &
                    T_F0_io(I, J, KM) + T_F0_io(I, J, K) )

                    D_INTP_io(I, J, K) = 0.1250_WP * ( D_F0_io(IM, JM, KM) + D_F0_io(IM, JM, K) + &
                    D_F0_io(IM, J, KM) + D_F0_io(IM, J, K) + &
                    D_F0_io(I, JM, KM) + D_F0_io(I, JM, K) + &
                    D_F0_io(I, J, KM) + D_F0_io(I, J, K) )

                    M_INTP_io(I, J, K) = 0.1250_WP * ( M_F0_io(IM, JM, KM) + M_F0_io(IM, JM, K) + &
                    M_F0_io(IM, J, KM) + M_F0_io(IM, J, K) + &
                    M_F0_io(I, JM, KM) + M_F0_io(I, JM, K) + &
                    M_F0_io(I, J, KM) + M_F0_io(I, J, K) )


                END DO
                IF(iCase == iPIPEC) THEN
                    G_INTP_io(I, 1, K, :) = 0.5_WP * ( G_F0_io(I, 1, K, :) + G_F0_io(I, 1, KS, :))
                ELSE
                    G_INTP_io(I, 1, K, :) = 0.0_WP
                END IF
            END DO
        END DO
        !========y b.c. pointS ============
        G_INTP_io(:, NND2, :, :) = 0.0_WP
        !======== X b.c. pointS ============
        IF(TgFlowFlg) THEN
            G_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3, 1:3)       = G_F0_io(NND1_io, 1 : NCL2, 1 : NCL3, 1:3)
            H_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = H_F0_io(NND1_io, 1 : NCL2, 1 : NCL3)
            T_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = T_F0_io(NND1_io, 1 : NCL2, 1 : NCL3)
            D_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = D_F0_io(NND1_io, 1 : NCL2, 1 : NCL3)
            M_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = M_F0_io(NND1_io, 1 : NCL2, 1 : NCL3)
        ELSE
            G_INTP_io(NND1_io, :, :, :)       = G_INTP_io(1, :, :, :)
            H_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = H_INTP_io(1, 1 : NCL2, 1 : NCL3)
            T_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = T_INTP_io(1, 1 : NCL2, 1 : NCL3)
            D_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = D_INTP_io(1, 1 : NCL2, 1 : NCL3)
            M_INTP_io(NND1_io, 1 : NCL2, 1 : NCL3) = M_INTP_io(1, 1 : NCL2, 1 : NCL3)
        END IF
        !======== Z b.c. pointS ============
        G_INTP_io(:, :, NND3, :) = G_INTP_io(:, :, 1, :)
        H_INTP_io(:, :, NND3) = H_INTP_io(:, :, 1)
        T_INTP_io(:, :, NND3) = T_INTP_io(:, :, 1)
        D_INTP_io(:, :, NND3) = D_INTP_io(:, :, 1)
        M_INTP_io(:, :, NND3) = M_INTP_io(:, :, 1)




        IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux .AND. iCase == iPIPEC) THEN
            H_INTP_io(:, 1, :) = 2.0_WP * H_INTP_io(:, 2, :) - H_INTP_io(:, 3, :)
            T_INTP_io(:, 1, :) = 2.0_WP * T_INTP_io(:, 2, :) - T_INTP_io(:, 3, :)
            D_INTP_io(:, 1, :) = 2.0_WP * D_INTP_io(:, 2, :) - D_INTP_io(:, 3, :)
            M_INTP_io(:, 1, :) = 2.0_WP * M_INTP_io(:, 2, :) -M_INTP_io(:, 3, :)
        ELSE IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            H_INTP_io(:, 1, :) = H_WAL_GV(1, iBotWall)
            T_INTP_io(:, 1, :) = T_WAL_GV(1, iBotWall)
            D_INTP_io(:, 1, :) = D_WAL_GV(1, iBotWall)
            M_INTP_io(:, 1, :) = M_WAL_GV(1, iBotWall)
        ELSE
        END IF


        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
            H_INTP_io(:, NND2, :) = 2.0_WP * H_INTP_io(:, NCL2, :) - H_INTP_io(:, NCL2 - 1, :)
            T_INTP_io(:, NND2, :) = 2.0_WP * T_INTP_io(:, NCL2, :) - T_INTP_io(:, NCL2 - 1, :)
            D_INTP_io(:, NND2, :) = 2.0_WP * D_INTP_io(:, NCL2, :) - D_INTP_io(:, NCL2 - 1, :)
            M_INTP_io(:, NND2, :) = 2.0_WP * M_INTP_io(:, NCL2, :) -M_INTP_io(:, NCL2 - 1, :)
        END IF
        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature) THEN
            H_INTP_io(:, NND2, :) = H_WAL_GV(1, iTopWall)
            T_INTP_io(:, NND2, :) = T_WAL_GV(1, iTopWall)
            D_INTP_io(:, NND2, :) = D_WAL_GV(1, iTopWall)
            M_INTP_io(:, NND2, :) = M_WAL_GV(1, iTopWall)
        END IF

    END IF

    !============ INTERPOLATION ALL VALUES TO POINTS FOR Z periodic B.C. X inlet /OUTLET ============================
    IF(TgFlowFlg) THEN

        !============= Z averaged valueS ===================
        ALLOCATE ( U1zL_INTP_io(NND1_io, NND2, NDV + 1 ) )
        DO I = 1, NND1_io
            IM = IMV_io(I)
            DO J = 2, NCL2
                JM = JGMV(J)
                U1zL_INTP_io(I, J, 1:4) = ( U1zL_F0_io(IM, JM, 1:4) + U1zL_F0_io(IM, J, 1:4) + &
                U1zL_F0_io(I, JM, 1:4) + U1zL_F0_io(I, J, 1:4) ) * 0.25_WP
            END DO
            U1zL_INTP_io(I, 1,   1:3) = 0.0_WP
            U1zL_INTP_io(I, 1,   4) = U1zL_INTP_io(I, 2,   4)
            U1zL_INTP_io(I, NND2, 1:3) = 0.0_WP
            U1zL_INTP_io(I, NND2, 4) = U1zL_INTP_io(I, NCL2, 4)
        END DO

        IF(iThermoDynamics == 1) THEN
            ALLOCATE ( G1zL_INTP_io(NND1_io, NND2, NDV ) )
            ALLOCATE ( H1zL_INTP_io(NND1_io, NND2) )
            ALLOCATE ( T1zL_INTP_io(NND1_io, NND2) )
            ALLOCATE ( D1zL_INTP_io(NND1_io, NND2) )
            ALLOCATE ( M1zL_INTP_io(NND1_io, NND2) )

            DO I = 1, NND1_io
                IM = IMV_io(I)
                DO J = 2, NCL2
                    JM = JGMV(J)
                    G1zL_INTP_io(I, J, 1:3) = ( G1zL_F0_io(IM, JM, 1:3) + G1zL_F0_io(IM, J, 1:3) + &
                    G1zL_F0_io(I, JM, 1:3) + G1zL_F0_io(I, J, 1:3) ) * 0.25_WP
                    H1zL_INTP_io(I, J    ) = ( H1zL_F0_io(IM, JM    ) + H1zL_F0_io(IM, J    ) + &
                    H1zL_F0_io(I, JM    ) + H1zL_F0_io(I, J    ) ) * 0.25_WP
                    T1zL_INTP_io(I, J    ) = ( T1zL_F0_io(IM, JM    ) + T1zL_F0_io(IM, J    ) + &
                    T1zL_F0_io(I, JM    ) + T1zL_F0_io(I, J    ) ) * 0.25_WP
                    D1zL_INTP_io(I, J    ) = ( D1zL_F0_io(IM, JM    ) + D1zL_F0_io(IM, J    ) + &
                    D1zL_F0_io(I, JM    ) + D1zL_F0_io(I, J    ) ) * 0.25_WP
                    M1zL_INTP_io(I, J    ) = ( M1zL_F0_io(IM, JM    ) + M1zL_F0_io(IM, J    ) + &
                    M1zL_F0_io(I, JM    ) + M1zL_F0_io(I, J    ) ) * 0.25_WP

                END DO
                G1zL_INTP_io(I, 1,   1:3) = 0.0_WP
                G1zL_INTP_io(I, NND2, 1:3) = 0.0_WP
                H1zL_INTP_io(I, 1)    = H1zL_F0_io(I, 1   )
                H1zL_INTP_io(I, NND2) = H1zL_F0_io(I, NCL2)
                T1zL_INTP_io(I, 1)    = T1zL_F0_io(I, 1   )
                T1zL_INTP_io(I, NND2) = T1zL_F0_io(I, NCL2)
                D1zL_INTP_io(I, 1)    = D1zL_F0_io(I, 1   )
                D1zL_INTP_io(I, NND2) = D1zL_F0_io(I, NCL2)
                M1zL_INTP_io(I, 1)    = M1zL_F0_io(I, 1   )
                M1zL_INTP_io(I, NND2) = M1zL_F0_io(I, NCL2)

            END DO
        END IF

    ELSE
        ALLOCATE ( U1xzL_INTP_io(NND2, NDV + 1 ) )
        DO J = 2, NCL2
            JM = J - 1
            U1xzL_INTP_io(J, 1:4) = (U1xzL_F0_io(JM, 1:4) + U1xzL_F0_io(J, 1:4)) * 0.5_WP
        END DO
        U1xzL_INTP_io(1,   1:3) = 0.0_WP
        U1xzL_INTP_io(1,     4) = U1xzL_INTP_io(2,     4)
        U1xzL_INTP_io(NND2, 1:3) = 0.0_WP
        U1xzL_INTP_io(NND2,  4) = U1xzL_INTP_io(NCL2,  4)

        IF(iThermoDynamics == 1) THEN
            ALLOCATE ( G1xzL_INTP_io(NND2, NDV ) )
            ALLOCATE ( H1xzL_INTP_io(NND2) )
            ALLOCATE ( T1xzL_INTP_io(NND2) )
            ALLOCATE ( D1xzL_INTP_io(NND2) )
            ALLOCATE ( M1xzL_INTP_io(NND2) )

            DO J = 2, NCL2
                JM = J - 1
                G1xzL_INTP_io(J, 1:3) = (G1xzL_F0_io(JM, 1:3) + G1xzL_F0_io(J, 1:3)) * 0.5_WP
                H1xzL_INTP_io(J    ) = (H1xzL_F0_io(JM    ) + H1xzL_F0_io(J    )) * 0.5_WP
                T1xzL_INTP_io(J    ) = (T1xzL_F0_io(JM    ) + T1xzL_F0_io(J    )) * 0.5_WP
                D1xzL_INTP_io(J    ) = (D1xzL_F0_io(JM    ) + D1xzL_F0_io(J    )) * 0.5_WP
                M1xzL_INTP_io(J    ) = (M1xzL_F0_io(JM    ) + M1xzL_F0_io(J    )) * 0.5_WP
            END DO
            G1xzL_INTP_io(1,   1:3) = 0.0_WP
            H1xzL_INTP_io(1) = HWAL_RA(iBotWall)
            T1xzL_INTP_io(1) = TWAL(iBotWall)
            D1xzL_INTP_io(1) = DWAL(iBotWall)
            M1xzL_INTP_io(1) = MWAL(iBotWall)

            G1xzL_INTP_io(NND2, 1:3) = 0.0_WP
            H1xzL_INTP_io(NND2) = HWAL_RA(iTopWall)
            T1xzL_INTP_io(NND2) = TWAL(iTopWall)
            D1xzL_INTP_io(NND2) = DWAL(iTopWall)
            M1xzL_INTP_io(NND2) = MWAL(iTopWall)
        END IF

    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TEC360_ALL_NODES
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, N1
    INTEGER(4) :: TECFLG = 202
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists
    REAL(WP) :: Uprime(3)

    ! Below IS testing....
    FLNM = 'RESULT.UVWP.Centreline.plt'

    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist = File_exists)
    IF(File_exists) THEN
        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)), POSITION = 'APPEND')
    ELSE
        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)))
    END IF
    WRITE(TECFLG, '(A, 1E13.5)') '#variables = "X", "Y", "Z", "U", "V","W", "P"', PhyTIME
    DO I = NCL1_io / 2, NCL1_io / 2
        DO K = NCL3 / 2, NCL3 / 2
            DO J = 1, NCL2
                WRITE(TECFLG, '(7ES15.7)') XCC_io(I), YCC(J),ZCC(K),U_F0_io(I, J, K, 1:4)
            END DO
        END DO
    END DO

    CLOSE(TECFLG)

    ! .....








    FLNM = TRIM(FilePath5) // 'RESULT.UVWP.ALLPOINTS.plt'

    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist =File_exists)

    IF(IoFlowFlg .AND. TgFlowFlg) THEN
        N1 = NND1_tg + NND1_io
    ELSE
        IF(TgFlowFlg) N1 = NND1_tg
        IF(IoFlowFlg) N1 = NND1_io
    END IF

    IF(File_exists) THEN
        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)), POSITION = 'APPEND')

        WRITE(TECFLG, '(A, 1I13.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') &
        'ZONE T = " ', ITERG,PhyTIME, ' ", I = ', &
        N1, ', J = ', NND2, ', K = ', NND3, ', F =POINT'
        WRITE(TECFLG, '(A)') 'VARSHARELIST = ([1-3]= 1)'
        DO K = 1, NND3
            DO J = 1, NND2
                IF(TgFlowFlg) THEN
                    DO I = 1, NND1_tg
                        WRITE(TECFLG, '(30ES15.7)') U_INTP_tg(I, J, K, 1),U_INTP_tg(I, J, K, 2),U_INTP_tg(I, J, K, 3), &
                        U_INTP_tg(I, J, K, 4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                        Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 3), &
                        Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                    END DO
                END IF
                IF(IoFlowFlg) THEN
                    DO I = 1, NND1_io
                        IF(iThermoDynamics == 1) THEN
                            WRITE(TECFLG, '(30ES15.7)') &
                            U_INTP_io(I, J, K, 1:4), &
                            T_INTP_io(I, J, K), D_INTP_io(I, J, K), H_INTP_io(I, J, K), &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                            Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                        ELSE
                            WRITE(TECFLG, '(30ES15.7)') &
                            U_INTP_io(I, J, K, 1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP,  &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                            0.0_WP, 0.0_WP, UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                        END IF
                    END DO
                END IF
            END DO
        END DO

    ELSE
        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        WRITE(TECFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Vorx","Vory","Vorz","VorM","Qcr","Delta","Lambda2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"SWIrlStrengthR1","SWIrlStrengthR2","SWIrlStrengthI2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime", '
        WRITE(TECFLG, '(A)')              '"Tprime", "Dprime","Mprime","uvprime", "UDprime", "vDprime"'

        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') &
        'ZONE T = " ', ITERG,PhyTIME, ' ", I = ', &
        N1, ', J = ', NND2, ', K = ', NND3, ', F =POINT'
        DO K = 1, NND3
            DO J = 1, NND2
                IF(TgFlowFlg) THEN
                    DO I = 1, NND1_tg
                        WRITE(TECFLG, '(33ES15.7)') XND_tg(I), YND(J),ZND(K), &
                        U_INTP_tg(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                        Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 1:3), &
                        Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                    END DO
                END IF
                IF(IoFlowFlg) THEN
                    DO I = 1, NND1_io
                        IF(iThermoDynamics == 1) THEN
                            WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                            U_INTP_io(I, J, K, 1:4), &
                            T_INTP_io(I, J, K), D_INTP_io(I, J, K), H_INTP_io(I, J, K), &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                            Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                        ELSE
                            WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                            U_INTP_io(I, J, K, 1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                            0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                        END IF
                    END DO
                END IF
            END DO
        END DO

    END IF


    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        WRITE(TECFLG, '(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID, ' // &
        'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK, ' // &
        'F = POINT, S = GLOBAL'
        WRITE(TECFLG, '(I2.1)') 1
        WRITE(TECFLG, '(I2.1)') 2
        WRITE(TECFLG, '(2ES15.7)') XND_io(1), YND(1)
        WRITE(TECFLG, '(2ES15.7)') XND_io(1), YND(NND2)
    END IF


    CLOSE(TECFLG)


END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TEC360_XSLICE(N)
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, N
    CHARACTER(4) :: PNTIM
    INTEGER(4) :: TECFLG = 202
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists
    REAL(WP) :: Uprime(3)



    IF(.NOT.IoFlowFlg) RETURN
    TECFLG = TECFLG + IID(N)

    WRITE(PNTIM, '(I4.4)') IID(N)

    FLNM = TRIM(FilePath5) // 'RESULT.UVWP.XSLICE.' //PNTIM // '.plt'

    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist = File_exists)

    IF(File_exists) THEN

        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)), POSITION = 'APPEND')
        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG,PhyTIME, &
        ' ", I = ', 1, ', J = ', NND2, ', K = ', NND3, ', F =POINT'
        WRITE(TECFLG, '(A)') 'VARSHARELIST = ([1-3]= 1)'

        DO K = 1, NND3
            DO J = 1, NND2
                DO I = iID(N), IID(N)
                    IF(iThermoDynamics == 1) THEN
                        WRITE(TECFLG, '(30ES15.7)') &
                        U_INTP_io(I, J, K, 1:4), &
                        T_INTP_io(I, J, K), D_INTP_io(I, J, K), H_INTP_io(I, J, K), &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                        Dprime_io(I, J, K), Mprime_io(I, J, K), UDprime_io(I, J, K, 1:3)
                    ELSE
                        WRITE(TECFLG, '(30ES15.7)') &
                        U_INTP_io(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                        0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                    END IF
                END DO
            END DO
        END DO


    ELSE


        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)) )
        WRITE(TECFLG, '(A)') 'TITLE = "DNS FLOW X -SLICE"'
        WRITE(TECFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Vorx","Vory","Vorz","VorM","Qcr","Delta","Lambda2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"SWIrlStrengthR1","SWIrlStrengthR2","SWIrlStrengthI2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime", '
        WRITE(TECFLG, '(A)')              '"Tprime", "Dprime","Mprime","uvprime", "UDprime", "vDprime"'


        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG, PhyTIME, &
        ' ", I = ', 1, ', J = ', NND2, ', K = ', NND3, ', F =POINT'

        DO K = 1, NND3
            DO J = 1, NND2
                DO I = iID(N), IID(N)
                    IF(iThermoDynamics == 1) THEN
                        WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                        U_INTP_io(I, J, K, 1:4), &
                        T_INTP_io(I, J, K), D_INTP_io(I, J, K),H_INTP_io(I, J, K), &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                        Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                    ELSE
                        WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                        U_INTP_io(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                        0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                    END IF
                END DO
            END DO
        END DO

    END IF

    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        WRITE(TECFLG, '(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID, ' // &
        'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK, ' // &
        'F = POINT, S = GLOBAL'
        WRITE(TECFLG, '(I2.1)') 1
        WRITE(TECFLG, '(I2.1)') 2
        WRITE(TECFLG, '(2ES15.7)') ZND(1), YND(1)
        WRITE(TECFLG, '(2ES15.7)') ZND(1), YND(NND2)

    END IF

    CLOSE(TECFLG)


END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE TEC360_YSLICE(N)
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, N1, N
    CHARACTER(4) :: PNTIM
    INTEGER(4) :: TECFLG = 202
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists
    REAL(WP) :: Uprime(3)

    IF(N < (MGRID/ 2)) THEN
        J = JGMOV(N) + 1
    ELSE
        J = JGMOV(N)
    END IF
    TECFLG = TECFLG + J
    WRITE(PNTIM, '(I4.4)') J
    FLNM = TRIM(FilePath5) // 'RESULT.UVWP.YSLICE.' //PNTIM // '.plt'

    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist = File_exists)

    IF(IoFlowFlg .AND. TgFlowFlg) THEN
        N1 = NND1_tg + NND1_io
    ELSE
        IF(TgFlowFlg) N1 = NND1_tg
        IF(IoFlowFlg) N1 = NND1_io
    END IF

    IF(File_exists) THEN

        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)), POSITION = 'APPEND')
        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG, PhyTIME, &
        ' ", I = ', N1, ', J = ', 1, ', K = ', NND3, ', F = POINT'
        WRITE(TECFLG, '(A)') 'VARSHARELIST = ([1-3]= 1)'

        DO K = 1, NND3
            IF(TgFlowFlg) THEN
                DO I = 1, NND1_tg
                    WRITE(TECFLG, '(30ES15.7)') U_INTP_tg(I, J, K, 1:4), &
                    1.0_WP, 1.0_WP, 0.0_WP, &
                    Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                    Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 1:3), &
                    Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                    0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                END DO
            END IF
            IF(IoFlowFlg) THEN
                DO I = 1, NND1_io
                    IF(iThermoDynamics == 1) THEN
                        WRITE(TECFLG, '(30ES15.7)') &
                        U_INTP_io(I, J, K, 1:4), &
                        T_INTP_io(I, J, K), D_INTP_io(I, J, K),H_INTP_io(I, J, K), &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                        Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                    ELSE
                        WRITE(TECFLG, '(30ES15.7)') &
                        U_INTP_io(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                        0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                    END IF
                END DO
            END IF
            !END DO
        END DO


    ELSE

        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)) )
        WRITE(TECFLG, '(A)') 'TITLE = "DNS FLOW Y-SLICE"'
        WRITE(TECFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Vorx","Vory","Vorz","VorM","Qcr","Delta","Lambda2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"SWIrlStrengthR1","SWIrlStrengthR2","SWIrlStrengthI2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime", '
        WRITE(TECFLG, '(A)')              '"Tprime", "Dprime","Mprime","uvprime", "UDprime", "vDprime"'

        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG, PhyTIME, &
        ' ", I = ', N1, ', J = ', 1, ', K = ', NND3, ', F =POINT'

        DO K = 1, NND3
            IF(TgFlowFlg) THEN
                DO I = 1, NND1_tg
                    WRITE(TECFLG, '(33ES15.7)') XND_tg(I), YND(J),ZND(K), &
                    U_INTP_tg(I, J, K, 1:4), &
                    1.0_WP, 1.0_WP, 0.0_WP, &
                    Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                    Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 1:3), &
                    Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                END DO
            END IF
            IF(IoFlowFlg) THEN
                DO I = 1, NND1_io
                    IF(iThermoDynamics == 1) THEN
                        WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                        U_INTP_io(I, J, K, 1:4), &
                        T_INTP_io(I, J, K), D_INTP_io(I, J, K),H_INTP_io(I, J, K), &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                        Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                    ELSE
                        WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                        U_INTP_io(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                        Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                        Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                        0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                    END IF
                END DO
            END IF
            !END DO
        END DO


    END IF

    IF(IoFlowFlg .AND. TgFlowFlg) THEN
        WRITE(TECFLG, '(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID, ' // &
        'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK, ' // &
        'F = POINT, S = GLOBAL'
        WRITE(TECFLG, '(I2.1)') 1
        WRITE(TECFLG, '(I2.1)') 2
        WRITE(TECFLG, '(2ES15.7)') XND_io(1),ZND(1)
        WRITE(TECFLG, '(2ES15.7)') XND_io(1),ZND(NND3)

    END IF

    CLOSE(TECFLG)


END SUBROUTINE TEC360_YSLICE

!**********************************************************************************************************************************
SUBROUTINE TEC360_ZSLICE(N)
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, N1, N
    CHARACTER(4) :: PNTIM
    INTEGER(4) :: TECFLG = 202
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists
    REAL(WP) :: Uprime(3)

    TECFLG = TECFLG + KID(N)
    IF(IoFlowFlg .AND. TgFlowFlg) THEN
        N1 = NND1_tg + NND1_io
    ELSE
        IF(TgFlowFlg) N1 = NND1_tg
        IF(IoFlowFlg) N1 = NND1_io
    END IF
    WRITE(PNTIM, '(I4.4)') KID(N)

    FLNM = TRIM(FilePath5) // 'RESULT.UVWP.ZSLICE.' //PNTIM // '.plt'

    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist = File_exists)

    IF(File_exists) THEN

        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)), POSITION = 'APPEND')
        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG,PhyTIME, &
        ' ", I = ', N1, ', J = ', NND2, ', K = ', 1, ', F =POINT'
        WRITE(TECFLG, '(A)') 'VARSHARELIST = ([1-3]= 1)'

        DO K = KID(N), KID(N)
            DO J = 1, NND2
                IF(TgFlowFlg) THEN
                    DO I = 1, NND1_tg
                        WRITE(TECFLG, '(30ES15.7)') U_INTP_tg(I, J, K, 1),U_INTP_tg(I, J, K, 2), U_INTP_tg(I, J, K, 3), &
                        U_INTP_tg(I, J, K, 4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                        Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 1:3), &
                        Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                    END DO
                END IF
                IF(IoFlowFlg) THEN
                    DO I = 1, NND1_io

                        IF(iThermoDynamics == 1) THEN
                            WRITE(TECFLG, '(30ES15.7)') &
                            U_INTP_io(I, J, K, 1:4), &
                            T_INTP_io(I, J, K), D_INTP_io(I, J, K),H_INTP_io(I, J, K), &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                            Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                        ELSE
                            WRITE(TECFLG, '(30ES15.7)') &
                            U_INTP_io(I, J, K, 1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                            0.0_WP, 0.0_WP,UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                        END IF
                    END DO
                END IF
            END DO
        END DO


    ELSE

        OPEN(TECFLG, FILE = TRIM(ADJUSTL(FLNM)) )
        WRITE(TECFLG, '(A)') 'TITLE = "DNS FLOW Z -SLICE"'
        WRITE(TECFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Vorx","Vory","Vorz","VorM","Qcr","Delta","Lambda2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"SWIrlStrengthR1","SWIrlStrengthR2","SWIrlStrengthI2", '
        WRITE(TECFLG, '(A)', AdvancE = "no") '"Uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime", '
        WRITE(TECFLG, '(A)')              '"Tprime", "Dprime","Mprime","uvprime", "UDprime", "vDprime"'

        WRITE(TECFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG,PhyTIME, &
        ' ", I = ', N1, ', J = ', NND2, ', K = ', 1, ', F =POINT'


        DO K = KID(N), KID(N)
            DO J = 1, NND2
                IF(TgFlowFlg) THEN
                    DO I = 1, NND1_tg
                        WRITE(TECFLG, '(33ES15.7)') XND_tg(I), YND(J),ZND(K), &
                        U_INTP_tg(I, J, K, 1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_tg(I, J, K, 1:4), Qcr_tg(I, J, K), Delta_tg(I, J, K), &
                        Lambda2_tg(I, J, K), SWIrlStrength_tg(I, J, K, 1:3), &
                        Uprime_tg(I, J, K, 1:4), 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP
                    END DO
                END IF
                IF(IoFlowFlg) THEN
                    DO I = 1, NND1_io
                        IF(iThermoDynamics == 1) THEN
                            WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                            U_INTP_io(I, J, K, 1:4), &
                            T_INTP_io(I, J, K), D_INTP_io(I, J, K),H_INTP_io(I, J, K), &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), Tprime_io(I, J, K), &
                            Dprime_io(I, J, K), Mprime_io(I, J, K),UDprime_io(I, J, K, 1:3)
                        ELSE
                            WRITE(TECFLG, '(33ES15.7)') XND_io(I), YND(J),ZND(K), &
                            U_INTP_io(I, J, K, 1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I, J, K, 1:4), Qcr_io(I, J, K), Delta_io(I, J, K), &
                            Lambda2_io(I, J, K), SWIrlStrength_io(I, J, K, 1:3), &
                            Uprime_io(I, J, K, 1:4), Gprime_io(I, J, K, 1:3), 0.0_WP, &
                            0.0_WP, 0.0_WP, UDprime_io(I, J, K, 1), 0.0_WP, 0.0_WP
                        END IF
                    END DO
                END IF
            END DO
        END DO

    END IF

    IF(IoFlowFlg .AND. TgFlowFlg) THEN
        WRITE(TECFLG, '(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID, ' // &
        'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK, ' // &
        'F = POINT, S = GLOBAL'
        WRITE(TECFLG, '(I2.1)') 1
        WRITE(TECFLG, '(I2.1)') 2
        WRITE(TECFLG, '(2ES15.7)') XND_io(1), YND(1)
        WRITE(TECFLG, '(2ES15.7)') XND_io(1), YND(NND2)

    END IF

    CLOSE(TECFLG)


END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TEC360_INSTANT_Uprime
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) I, J, K

    IF(TgFlowFlg) THEN
        ALLOCATE ( Uprime_tg(NND1_tg, NND2, NND3, NDV + 1) )
        DO I = 1, NND1_tg
            DO K = 1, NND3
                DO J = 1, NND2
                    Uprime_tg(I, J, K, 1:4) = U_INTP_tg(I, J, K, 1:4) - U1xzL_INTP_tg(J, 1:4)
                END DO
            END DO
        END DO
    END IF

    IF(IoFlowFlg) THEN
        ALLOCATE ( Uprime_io(NND1_io, NND2, NND3, NDV + 1) )
        ALLOCATE ( UDprime_io(NND1_io, NND2, NND3, NDV) )

        IF(TgFlowFlg) THEN
            DO I = 1, NND1_io
                DO K = 1, NND3
                    DO J = 1, NND2
                        Uprime_io(I, J, K, 1:4) = U_INTP_io(I, J, K, 1:4) - U1zL_INTP_io(I, J, 1:4)
                    END DO
                END DO
            END DO
            UDprime_io(:, :, :, 1) = Uprime_io(:, :, :, 1) * Uprime_io(:, :, :, 2) !u'v'

            ALLOCATE ( Gprime_io(NND1_io, NND2, NND3, NDV) )
            IF(iThermoDynamics == 1) THEN
                ALLOCATE ( Tprime_io(NND1_io, NND2, NND3) )
                ALLOCATE ( Dprime_io(NND1_io, NND2, NND3) )
                ALLOCATE ( Mprime_io(NND1_io, NND2, NND3) )

                DO I = 1, NND1_io
                    DO K = 1, NND3
                        DO J = 1, NND2
                            !Gprime_io(I, J, K, 1:3) = U_INTP_io(I, J, K, 1:3) - G1zL_INTP_io(I, J, 1:3) / D1zL_INTP_io(I, J)
                            Gprime_io(I, J, K, 1:3) = G_INTP_io(I, J, K, 1:3) - G1zL_INTP_io(I, J, 1:3)
                            Tprime_io(I, J, K)    = T_INTP_io(I, J, K) -     T1zL_INTP_io(I, J)
                            Dprime_io(I, J, K)    = D_INTP_io(I, J, K) -     D1zL_INTP_io(I, J)
                        END DO
                    END DO
                END DO
                UDprime_io(:, :, :, 2) = Uprime_io(:, :, :, 1) * Dprime_io(:, :, :) !u'd'
                UDprime_io(:, :, :, 3) = Uprime_io(:, :, :, 2) * Dprime_io(:, :, :) !v'd'
            ELSE
                Gprime_io = Uprime_io
            END IF

        ELSE
            DO I = 1, NND1_io
                DO K = 1, NND3
                    DO J = 1, NND2
                        Uprime_io(I, J, K, 1:4) = U_INTP_io(I, J, K, 1:4) - U1xzL_INTP_io(J, 1:4)
                    END DO
                END DO
            END DO
            UDprime_io(:, :, :, 1) = Uprime_io(:, :, :, 1) * Uprime_io(:, :, :, 2) !u'v'

            ALLOCATE ( Gprime_io(NND1_io, NND2, NND3, NDV) )
            IF(iThermoDynamics == 1) THEN
                ALLOCATE ( Tprime_io(NND1_io, NND2, NND3) )
                ALLOCATE ( Dprime_io(NND1_io, NND2, NND3) )
                ALLOCATE ( Mprime_io(NND1_io, NND2, NND3) )
                DO I = 1, NND1_io
                    DO K = 1, NND3
                        DO J = 1, NND2
                            !Gprime_io(I, J, K, 1:3) = U_INTP_io(I, J, K, 1:3) - G1xzL_INTP_io(J, 1:3) / D1xzL_INTP_io(J)
                            Gprime_io(I, J, K, 1:3) = G_INTP_io(I, J, K, 1:3) - G1xzL_INTP_io(J, 1:3)
                            Tprime_io(I, J, K)    = T_INTP_io(I, J, K) -     T1xzL_INTP_io(J)
                            Dprime_io(I, J, K)    = D_INTP_io(I, J, K) -     D1xzL_INTP_io(J)
                            Mprime_io(I, J, K)    = M_INTP_io(I, J, K) -     M1xzL_INTP_io(J)
                        END DO
                    END DO
                END DO
                UDprime_io(:, :, :, 2) = Uprime_io(:, :, :, 1) * Dprime_io(:, :, :) !u'd'
                UDprime_io(:, :, :, 3) = Uprime_io(:, :, :, 2) * Dprime_io(:, :, :) !v'd'
            ELSE
                Gprime_io = Uprime_io
            END IF
        END IF

    END IF



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TEC360_INSTANT_VORTEX_CRITERIA
    USE TEC360_INFO
    IMPLICIT NONE

    INTEGER(4) :: IC, JC, KC
    INTEGER(4) :: IM, JM, KM
    INTEGER(4) :: IP, JP, KP
    REAL(WP) :: DUDX1, DUDX2, DUDX3, DUDX4, DUDX5, DUDX6, DUDX7, DUDX8
    REAL(WP) :: DVDY1, DVDY2, DVDY3, DVDY4, DVDY5, DVDY6, DVDY7, DVDY8
    REAL(WP) :: DWDZ1, DWDZ2, DWDZ3, DWDZ4, DWDZ5, DWDZ6, DWDZ7, DWDZ8
    REAL(WP) :: DUDY1, DUDY2, DUDZ1, DUDZ2
    REAL(WP) :: DWDY1, DWDY2, DVDZ1, DVDZ2
    REAL(WP) :: DWDX1, DWDX2, DVDX1, DVDX2
    REAL(WP) :: StrSsym(3, 3)  !symmetrIC components of velcocity gradient
    REAL(WP) :: StrAsym(3, 3)  !antISymmetrIC components of velcocity gradient
    REAL(WP) :: gradU(3, 3)
    REAL(WP) :: S2plusO2(3, 3)
    REAL(WP) :: EIG(3)
    REAL(WP) :: detGradU
    REAL(WP) :: TraceGradU
    REAL(WP) :: STSSYM
    REAL(WP) :: STASYM
    REAL(WP), ALLOCATABLE :: dQ_io(:, :, :, :, :)
    REAL(WP), ALLOCATABLE :: dQ_tg(:, :, :, :, :)

    REAL(WP) :: AR(3, 3), AI(3, 3)
    REAL(WP) :: WR(3), WI(3)
    REAL(WP) :: ZR(3, 3), ZI(3, 3)
    INTEGER(4) :: NM, NN
    INTEGER(4) :: MATZ
    INTEGER(4) :: IERR
    REAL(WP) :: FV1(3), FV2(3), FV3(3)


    IF(IoFlowFlg) THEN
        allocate (Qcr_io    (NND1_io, NND2, NND3)       )
        allocate (Vor_io    (NND1_io, NND2, NND3, 4)     )
        allocate (Delta_io(NND1_io, NND2, NND3)       )
        allocate (Lambda2_io(NND1_io, NND2, NND3)       )
        allocate (SWIrlStrength_io(NND1_io, NND2, NND3, 3) )

        allocate (dQ_io    (NND1_io, NND2, NND3, 3, 3)   )

        !====================GET ALL du_I / DX_i on nodeS ===================================================
        DO IC = 1, NCL1_io
            IM = IMV_io(IC)
            IP = IPV_io(IC)
            DO JC = 2, NCL2
                JM = JGMV(JC)
                JP = JGPV(JC)
                DO KC = 1, NCL3
                    KP = KPV(KC)
                    KM = KMV(KC)
                    !===============================================================================
                    !=== DU / DX at i', j', k'============================
                    DUDX1 = (U_F0_io(IC, JC, KM, 1) - U_F0_io(IM, JC, KM, 1)) * DXI       ! (I - 1, j, K - 1)
                    DUDX2 = (U_F0_io(IP, JC, KM, 1) - U_F0_io(IC, JC, KM, 1)) * DXI       ! (i,   j, K - 1)

                    DUDX3 = (U_F0_io(IC, JC, KC, 1) - U_F0_io(IM, JC, KC, 1)) * DXI       ! (I - 1, j, k)
                    DUDX4 = (U_F0_io(IP, JC, KC, 1) - U_F0_io(IC, JC, KC, 1)) * DXI       ! (i,   j, k)

                    DUDX5 = (U_F0_io(IC, JM, KM, 1) - U_F0_io(IM, JM, KM, 1)) * DXI       ! (I - 1, J - 1, K - 1)
                    DUDX6 = (U_F0_io(IP, JM, KM, 1) - U_F0_io(IC, JM, KM, 1)) * DXI       ! (i,   J - 1, K - 1)

                    DUDX7 = (U_F0_io(IC, JM, KC, 1) - U_F0_io(IM, JM, KC, 1)) * DXI       ! (I - 1, J - 1, k)
                    DUDX8 = (U_F0_io(IP, JM, KC, 1) - U_F0_io(IC, JM, KC, 1)) * DXI       ! (i,   J - 1, k)

                    dQ_io(IC, JC, KC, 1, 1) = ( (DUDX1 + DUDX2) * 0.5_WP + &
                    (DUDX3 + DUDX4) * 0.5_WP ) * 0.5_WP * YCL2ND_WFF(JC) + &
                    ( (DUDX3 + DUDX4) * 0.5_WP + (DUDX5+ DUDX6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== DU / Dy at i', J', K'============================
                    IF(TgFlowFlg) THEN
                        DUDY1 = ( (U_F0_io(IC, JC, KM, 1) - U_F0_io(IC, JM, KM, 1)) - &
                        (U1zL_F0_io(IC, JC, 1) - U1zL_F0_io(IC, JM, 1)) ) * DYCI(JC)  ! (i', J', K - 1)
                        DUDY2 = ( (U_F0_io(IC, JC, KC, 1) - U_F0_io(IC, JM, KC, 1)) - &
                        (U1zL_F0_io(IC, JC, 1) - U1zL_F0_io(IC, JM, 1)) ) * DYCI(JC)  ! (i', J', K)

                    ELSE
                        DUDY1 = ( (U_F0_io(IC, JC, KM, 1) - U_F0_io(IC, JM, KM, 1)) - &
                        (U1xzL_F0_io(JC, 1) - U1xzL_F0_io(JM, 1)) ) * DYCI(JC)  ! (i', J', K - 1)
                        DUDY2 = ( (U_F0_io(IC, JC, KC, 1) - U_F0_io(IC, JM, KC, 1)) - &
                        (U1xzL_F0_io(JC, 1) - U1xzL_F0_io(JM, 1)) ) * DYCI(JC)  ! (i', J', K)
                    END IF

                    dQ_io(IC, JC, KC, 1, 2) = 0.5_WP * (DUDY1 + DUDY2)                 ! (i', J', K')

                    !=== DU / Dz at i', J', K'============================
                    DUDZ1 = (U_F0_io(IC, JM, KC, 1) - U_F0_io(IC, JM, KM, 1)) * DZI       ! (i', J - 1, K')
                    DUDZ2 = (U_F0_io(IC, JC, KC, 1) - U_F0_io(IC, JC, KM, 1)) * DZI       ! (i', J,  k')

                    dQ_io(IC, JC, KC, 1, 3) = DUDZ1 * YCL2ND_WFB(JC) + DUDZ2 * YCL2ND_WFF(JC)

                    !===============================================================================
                    !=== Dv/ Dy at i', j', k'============================
                    IF(TgFlowFlg) THEN
                        IF(JC == NCL2) THEN
                            DVDY1 = ( (0.0_WP - U_F0_io(IM, JC, KM, 2)) - &
                            (0.0_WP - U1zL_F0_io(IC, JC, 2)) ) * DYFI(JC)       ! (I - 1, j, K - 1)
                            DVDY2 = ( (0.0_WP - U_F0_io(IC, JC, KM, 2)) - &
                            (0.0_WP - U1zL_F0_io(IC, JC, 2)) ) * DYFI(JC)       ! (i,   j, K - 1)
                            DVDY3 = (0.0_WP - U_F0_io(IM, JC, KC, 2)) * DYFI(JC)       ! (I - 1, j, k)
                            DVDY4 = (0.0_WP - U_F0_io(IC, JC, KC, 2)) * DYFI(JC)       ! (i,   j, k)
                        ELSE
                            DVDY1 = ( (U_F0_io(IM, JP, KM, 2) - U_F0_io(IM, JC, KM, 2)) - &
                            (U1zL_F0_io(IC, JP, 2) - U1zL_F0_io(IC, JC, 2)) )&
                            * DYFI(JC)       ! (I - 1, j, K - 1)
                            DVDY2 = ( (U_F0_io(IC, JP, KM, 2) - U_F0_io(IC, JC, KM, 2)) - &
                            (U1zL_F0_io(IC, JP, 2) - U1zL_F0_io(IC, JC, 2)) )&
                            * DYFI(JC)       ! (i,   j, K - 1)

                            DVDY3 = ( (U_F0_io(IM, JP, KC, 2) - U_F0_io(IM, JC, KC, 2)) - &
                            (U1zL_F0_io(IC, JP, 2) - U1zL_F0_io(IC, JC, 2)) )&
                            * DYFI(JC)       ! (I - 1, j, k)
                            DVDY4 = ( (U_F0_io(IC, JP, KC, 2) - U_F0_io(IC, JC, KC, 2)) - &
                            (U1zL_F0_io(IC, JP, 2) - U1zL_F0_io(IC, JC, 2)) )&
                            * DYFI(JC)       ! (i,   j, k)
                        END IF
                        DVDY5= ( (U_F0_io(IM, JC, KM, 2) - U_F0_io(IM, JM, KM, 2)) - &
                        (U1zL_F0_io(IC, JC, 2) - U1zL_F0_io(IC, JM, 2)) ) * DYFI(JM)      ! (I - 1, J - 1, K - 1)
                        DVDY6= ( (U_F0_io(IC, JC, KM, 2) - U_F0_io(IC, JM, KM, 2)) - &
                        (U1zL_F0_io(IC, JC, 2) - U1zL_F0_io(IC, JM, 2)) ) * DYFI(JM)      ! (i,   J - 1, K - 1)

                        DVDY7= ( (U_F0_io(IM, JC, KC, 2) - U_F0_io(IM, JM, KC, 2)) - &
                        (U1zL_F0_io(IC, JC, 2) - U1zL_F0_io(IC, JM, 2)) ) * DYFI(JM)      ! (I - 1, J - 1, k)
                        DVDY8= ( (U_F0_io(IC, JC, KC, 2) - U_F0_io(IC, JM, KC, 2)) - &
                        (U1zL_F0_io(IC, JC, 2) - U1zL_F0_io(IC, JM, 2)) ) * DYFI(JM)      ! (i,   J - 1, k)
                    ELSE
                        IF(JC == NCL2) THEN
                            DVDY1 = ( (0.0_WP - U_F0_io(IM, JC, KM, 2)) - (0.0_WP - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (I - 1, j, K - 1)
                            DVDY2 = ( (0.0_WP - U_F0_io(IC, JC, KM, 2)) - (0.0_WP - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (i,   j, K - 1)
                            DVDY3 = (0.0_WP - U_F0_io(IM, JC, KC, 2)) * DYFI(JC)       ! (I - 1, j, k)
                            DVDY4 = (0.0_WP - U_F0_io(IC, JC, KC, 2)) * DYFI(JC)       ! (i,   j, k)
                        ELSE
                            DVDY1 = ( (U_F0_io(IM, JP, KM, 2) - U_F0_io(IM, JC, KM, 2)) - &
                            (U1xzL_F0_io(JP, 2) - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (I - 1, j, K - 1)
                            DVDY2 = ( (U_F0_io(IC, JP, KM, 2) - U_F0_io(IC, JC, KM, 2)) - &
                            (U1xzL_F0_io(JP, 2) - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (i,   j, K - 1)

                            DVDY3 = ( (U_F0_io(IM, JP, KC, 2) - U_F0_io(IM, JC, KC, 2)) - &
                            (U1xzL_F0_io(JP, 2) - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (I - 1, j, k)
                            DVDY4 = ( (U_F0_io(IC, JP, KC, 2) - U_F0_io(IC, JC, KC, 2)) - &
                            (U1xzL_F0_io(JP, 2) - U1xzL_F0_io(JC, 2)) ) * DYFI(JC)       ! (i,   j, k)
                        END IF
                        DVDY5 = ( (U_F0_io(IM, JC, KM, 2) - U_F0_io(IM, JM, KM, 2)) - &
                        (U1xzL_F0_io(JC, 2) - U1xzL_F0_io(JM, 2)) ) * DYFI(JM)      ! (I - 1, J - 1, K - 1)
                        DVDY6 = ( (U_F0_io(IC, JC, KM, 2) - U_F0_io(IC, JM, KM, 2)) - &
                        (U1xzL_F0_io(JC, 2) - U1xzL_F0_io(JM, 2)) ) * DYFI(JM)      ! (i,   J - 1, K - 1)

                        DVDY7 = ( (U_F0_io(IM, JC, KC, 2) - U_F0_io(IM, JM, KC, 2)) - &
                        (U1xzL_F0_io(JC, 2) - U1xzL_F0_io(JM, 2)) ) * DYFI(JM)      ! (I - 1, J - 1, k)
                        DVDY8 = ( (U_F0_io(IC, JC, KC, 2) - U_F0_io(IC, JM, KC, 2)) - &
                        (U1xzL_F0_io(JC, 2) - U1xzL_F0_io(JM, 2)) ) * DYFI(JM)      ! (i,   J - 1, k)

                    END IF

                    dQ_io(IC, JC, KC, 2, 2) = ( (DVDY1 + DVDY2) * 0.5_WP + (DVDY3 + DVDY4) * 0.5_WP ) * &
                    0.5_WP * YCL2ND_WFF(JC) + &
                    ( (DVDY3 + DVDY4) * 0.5_WP + (DVDY5+ DVDY6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== Dv/ DX at i', J', K'============================
                    DVDX1 = (U_F0_io(IC, JC, KM, 2) - U_F0_io(IM, JC, KM, 2)) * DXI  ! (i', J', K - 1)
                    DVDX2 = (U_F0_io(IC, JC, KC, 2) - U_F0_io(IM, JC, KC, 2)) * DXI  ! (i', J', K)

                    dQ_io(IC, JC, KC, 2, 1) = 0.5_WP * (DVDX1 + DVDX2)                 ! (i', J', K')

                    !=== Dv/ Dz at i', J', K'============================
                    DVDZ1 = (U_F0_io(IM, JC, KC, 2) - U_F0_io(IM, JC, KM, 2)) * DZI       ! (I - 1, J', K')
                    DVDZ2 = (U_F0_io(IC, JC, KC, 2) - U_F0_io(IC, JC, KM, 2)) * DZI       ! (i, J', K')

                    dQ_io(IC, JC, KC, 2, 3) = (DVDZ1 + DVDZ2) * 0.5_WP

                    !===============================================================================
                    !=== Dw/ Dz at i', j', k'============================
                    DWDZ1 = (U_F0_io(IM, JC, KC, 3) - U_F0_io(IM, JC, KM, 3)) * DXI       ! (I - 1, j, K - 1)
                    DWDZ2 = (U_F0_io(IC, JC, KC, 3) - U_F0_io(IC, JC, KM, 3)) * DXI       ! (i,   j, K - 1)

                    DWDZ3 = (U_F0_io(IM, JC, KP, 3) - U_F0_io(IM, JC, KC, 3)) * DXI       ! (I - 1, j, k)
                    DWDZ4 = (U_F0_io(IC, JC, KP, 3) - U_F0_io(IC, JC, KC, 3)) * DXI       ! (i,   j, k)

                    DWDZ5 = (U_F0_io(IM, JM, KC, 3) - U_F0_io(IM, JM, KM, 3)) * DXI       ! (I - 1, J - 1, K - 1)
                    DWDZ6 = (U_F0_io(IC, JM, KC, 3) - U_F0_io(IC, JM, KM, 3)) * DXI       ! (i,   J - 1, K - 1)

                    DWDZ7 = (U_F0_io(IM, JM, KP, 3) - U_F0_io(IM, JM, KC, 3)) * DXI       ! (I - 1, J - 1, k)
                    DWDZ8 = (U_F0_io(IC, JM, KP, 3) - U_F0_io(IC, JM, KC, 3)) * DXI       ! (i,   J - 1, k)

                    dQ_io(IC, JC, KC, 3, 3) = ( (DWDZ1 + DWDZ2) * 0.5_WP + (DWDZ3 + DWDZ4) * 0.5_WP ) * 0.5_WP * &
                    YCL2ND_WFF(JC) + &
                    ( (DWDZ3 + DWDZ4) * 0.5_WP + (DWDZ5+ DWDZ6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== Dw/ Dy at i', J', K'============================
                    IF(TgFlowFlg) THEN
                        DWDY1 = ( (U_F0_io(IM, JC, KC, 3) - U_F0_io(IM, JM, KC, 3)) - &
                        (U1zL_F0_io(IC, JC, 3) - U1zL_F0_io(IC, JM, 3)) ) * DYCI(JC)  ! (I - 1, J', K')
                        DWDY2 = ( (U_F0_io(IC, JC, KC, 3) - U_F0_io(IC, JM, KC, 3)) - &
                        (U1zL_F0_io(IC, JC, 3) - U1zL_F0_io(IC, JM, 3)) ) * DYCI(JC)  ! (i, J', K)
                    ELSE
                        DWDY1 = ( (U_F0_io(IM, JC, KC, 3) - U_F0_io(IM, JM, KC, 3)) - &
                        (U1xzL_F0_io(JC, 3) - U1xzL_F0_io(JM, 3)) ) * DYCI(JC)  ! (I - 1, J', K')
                        DWDY2 = ( (U_F0_io(IC, JC, KC, 3) - U_F0_io(IC, JM, KC, 3)) - &
                        (U1xzL_F0_io(JC, 3) - U1xzL_F0_io(JM, 3)) ) * DYCI(JC)  ! (i, J', K)
                    END IF

                    dQ_io(IC, JC, KC, 3, 2) = 0.5_WP * (DWDY1 + DWDY2)                 ! (i', J', K')

                    !=== Dw/ DX at i', J', K'============================
                    DWDX1 = (U_F0_io(IC, JM, KC, 3) - U_F0_io(IM, JM, KC, 3)) * DZI       ! (i', J - 1, K')
                    DWDX2 = (U_F0_io(IC, JC, KC, 3) - U_F0_io(IM, JC, KC, 3)) * DZI       ! (i', J,  k')

                    dQ_io(IC, JC, KC, 3, 1) = DWDX1 * YCL2ND_WFB(JC) + DWDX2 * YCL2ND_WFF(JC)
                END DO
            END DO
        END DO
        dQ_io(:, NND2, :, :, :) = 2.0_WP * DQ_io(:, NCL2, :, :, :) - DQ_io(:, NCL2 - 1, :, :, :)
        dQ_io(:, 1,   :, :, :) = 2.0_WP * DQ_io(:, 2,   :, :, :) - DQ_io(:, 3,     :, :, :)
        IF(.NOT.TgFlowFlg) dQ_io(NND1_io, :,:,   :, :) = dQ_io(1, :, :, :, :)
        dQ_io(:,      :, NND3, :, :) = dQ_io(:, :, 1, :, :)



        DO IC = 1, NND1_io
            DO JC = 1, NND2
                DO KC = 1, NND3
                    !=================== VortICitY ===============================
                    Vor_io(IC, JC, KC, 1) = dQ_io(IC, JC, KC, 3, 2) - dQ_io(IC, JC, KC, 2, 3) ! Vor_io-x
                    Vor_io(IC, JC, KC, 2) = dQ_io(IC, JC, KC, 1, 3) - dQ_io(IC, JC, KC, 3, 1) ! Vor_io- Y
                    Vor_io(IC, JC, KC, 3) = dQ_io(IC, JC, KC, 2, 1) - dQ_io(IC, JC, KC, 1, 2) ! Vor_io-z
                    Vor_io(IC, JC, KC, 4) = DSQRT( DOT_PRODUCT( Vor_io(IC, JC, KC, 1:3), Vor_io(IC, JC, KC, 1:3) ) )!Vor_io-mag

                    !==========================================================================================
                    !===================omega = antISymetrIC components of graduate U===============
                    StrAsym(1, 1) = dQ_io(IC, JC, KC, 1, 1) - dQ_io(IC, JC, KC, 1, 1)  !dU / DX - dU / DX
                    StrAsym(1, 2) = dQ_io(IC, JC, KC, 2, 1) - dQ_io(IC, JC, KC, 1, 2)  !dV/ DX - dU / DY
                    StrAsym(1, 3) = dQ_io(IC, JC, KC, 3, 1) - dQ_io(IC, JC, KC, 1, 3)  !dW/ DX - dU / DZ

                    StrAsym(2, 1) = dQ_io(IC, JC, KC, 1, 2) - dQ_io(IC, JC, KC, 2, 1)  !dU / DY - dV/ DX
                    StrAsym(2, 2) = dQ_io(IC, JC, KC, 2, 2) - dQ_io(IC, JC, KC, 2, 2)  !dV/ DY - dV/ DY
                    StrAsym(2, 3) = dQ_io(IC, JC, KC, 3, 2) - dQ_io(IC, JC, KC, 2, 3)  !dW/ DY - dV/ DZ

                    StrAsym(3, 1) = dQ_io(IC, JC, KC, 1, 3) - dQ_io(IC, JC, KC, 3, 1)  !dU / DZ - dW/ DX
                    StrAsym(3, 2) = dQ_io(IC, JC, KC, 2, 3) - dQ_io(IC, JC, KC, 3, 2)  !dV/ DZ - dW/ DY
                    StrAsym(3, 3) = dQ_io(IC, JC, KC, 3, 3) - dQ_io(IC, JC, KC, 3, 3)  !dW/ DZ - dW/ DZ
                    !===================omega = symetrIC components of graduate U===============
                    StrSsym(1, 1) = dQ_io(IC, JC, KC, 1, 1) + dQ_io(IC, JC, KC, 1, 1)  !dU / DX + dU / DX
                    StrSsym(1, 2) = dQ_io(IC, JC, KC, 2, 1) + dQ_io(IC, JC, KC, 1, 2)  !dV/ DX + dU / DY
                    StrSsym(1, 3) = dQ_io(IC, JC, KC, 3, 1) + dQ_io(IC, JC, KC, 1, 3)  !dW/ DX + dU / DZ

                    StrSsym(2, 1) = dQ_io(IC, JC, KC, 1, 2) + dQ_io(IC, JC, KC, 2, 1)  !dU / DY + dV/ DX
                    StrSsym(2, 2) = dQ_io(IC, JC, KC, 2, 2) + dQ_io(IC, JC, KC, 2, 2)  !dV/ DY + dV/ DY
                    StrSsym(2, 3) = dQ_io(IC, JC, KC, 3, 2) + dQ_io(IC, JC, KC, 2, 3)  !dW/ DY + dV/ DZ

                    StrSsym(3, 1) = dQ_io(IC, JC, KC, 1, 3) + dQ_io(IC, JC, KC, 3, 1)  !dU / DZ + dW/ DX
                    StrSsym(3, 2) = dQ_io(IC, JC, KC, 2, 3) + dQ_io(IC, JC, KC, 3, 2)  !dV/ DZ + dW/ DY
                    StrSsym(3, 3) = dQ_io(IC, JC, KC, 3, 3) + dQ_io(IC, JC, KC, 3, 3)  !dW/ DZ + dW/ DZ

                    StrAsym  = StrAsym  * 0.5_WP  ! Sij
                    StrSsym  = StrSsym  * 0.5_WP  ! WIj

                    ! trace of (STSYM *STSYM^T)
                    STSSYM = DOT_PRODUCT(StrSsym(:, 1), StrSsym(:, 1)) + &
                    DOT_PRODUCT(StrSsym(:, 2), StrSsym(:, 2)) + &
                    DOT_PRODUCT(StrSsym(:, 3), StrSsym(:, 3))    ! SiJ *Sij
                    ! trace of (STASYM *STASYM^T)
                    STASYM = DOT_PRODUCT(StrAsym(:, 1), StrAsym(:, 1)) + &
                    DOT_PRODUCT(StrAsym(:, 2), StrAsym(:, 2)) + &
                    DOT_PRODUCT(StrAsym(:, 3), StrAsym(:, 3))    ! WIJ * WIj

                    !================================ Q cretiria===================================================
                    QCR_io(IC, JC, KC) = 0.5_WP * (STASYM -STSSYM) ! Q


                    !================================ Delta_io criteria===============================================
                    ! grad U
                    gradU = StrAsym + StrSsym

                    !  negtive trace of velocity gradient  !P
                    TraceGradU = (StrAsym(1, 1) +StrAsym(2, 2) +StrAsym(3, 3) +&
                    StrSsym(1, 1) +StrSsym(2, 2) +StrSsym(3, 3)) * (-1.0_WP)

                    ! det(gradU) !R
                    detGradU = (gradU(1, 1) * ( gradU(2, 2) * GradU(3, 3) - gradU(2, 3) * GradU(3, 2) ) - &
                    gradU(1, 2) * ( gradU(2, 1) * GradU(3, 3) - gradU(3, 1) * GradU(2, 3) ) + &
                    gradU(1, 3) * ( gradU(2, 1) * GradU(3, 2) - gradU(2, 2) * GradU(3, 1) )) * (-1.0_WP)
                    ! Delta for incompressible flow
                    !Delta_io(IC, JC, KC) = (QCR_io(IC, JC, KC) / 3.0_WP)**3 + (detGradU / 2.0_WP)**2

                    ! Delta for compressible flow
                    Delta_io(IC, JC, KC) = 27.0_WP * DetGradU* DetGradU + &
                    (4.0_WP * TraceGradU* TraceGradU* TraceGradU- 18.0_WP * TraceGradU* QCR_io(IC, JC, KC)) * DetGradU + &
                    4.0_WP * QCR_io(IC, JC, KC) * QCR_io(IC, JC, KC) * QCR_io(IC, JC, KC) - &
                    TraceGradU* TraceGradU* QCR_io(IC, JC, KC) * QCR_io(IC, JC, KC)

                    !get rid of possible infinity
                    IF (Delta_io(IC, JC, KC) > 1.0E30_WP) Delta_io(IC, JC, KC) = 1.0E30_WP
                    IF (Delta_io(IC, JC, KC) < -1.0E30_WP) Delta_io(IC, JC, KC) = -1.0E30_WP

                    !============== SWIrl stRENgth ==================================
                    IF (Delta_io(IC, JC, KC) > 0.0_WP) THEN

                        MATZ = 0
                        NM = 3
                        NN = 3
                        AR(:, :) = StrAsym(:, :) + StrSsym(:, :)
                        AI(:, :) = 0.0_WP

                        CALL cg(NM, NN, AR, AI, WR, WI, MATZ, ZR, ZI, FV1, FV2, FV3, IERR)

                        IF(DABS(WI(1)) < (DABS(WI(2)) + 1.0E-14) .AND. DABS(WI(1)) < (DABS(WI(3)) + 1.0E-14) ) THEN
                            SWIrlStrength_io(IC, JC, KC, 1) = WR(1)
                            SWIrlStrength_io(IC, JC, KC, 2) = WR(2)
                            SWIrlStrength_io(IC, JC, KC, 3) = WR(2)
                        ELSE IF (DABS(WI(2)) < (DABS(WI(1)) + 1.0E-14) .AND. DABS(WI(2)) < (DABS(WI(3)) + 1.0E-14) ) THEN
                            SWIrlStrength_io(IC, JC, KC, 1) = WR(2)
                            SWIrlStrength_io(IC, JC, KC, 2) = WR(1)
                            SWIrlStrength_io(IC, JC, KC, 3) = WR(1)
                        ELSE IF (DABS(WI(3)) < (DABS(WI(1)) + 1.0E-14) .AND. DABS(WI(3)) < (DABS(WI(2)) + 1.0E-14) ) THEN
                            SWIrlStrength_io(IC, JC, KC, 1) = WR(3)
                            SWIrlStrength_io(IC, JC, KC, 2) = WR(1)
                            SWIrlStrength_io(IC, JC, KC, 3) = WR(1)
                        ELSE
                        END IF

                    END IF


                    !************* Lamda2 criteria*************************
                    S2plusO2(1, 1) = StrSsym(1, 1) * StrSsym(1, 1) + &
                    StrSsym(1, 2) * StrSsym(2, 1) + &
                    StrSsym(1, 3) * StrSsym(3, 1) + &
                    StrAsym(1, 1) * StrAsym(1, 1) + &
                    StrAsym(1, 2) * StrAsym(2, 1) + &
                    StrAsym(1, 3) * StrAsym(3, 1)
                    S2plusO2(1, 2) = StrSsym(1, 1) * StrSsym(1, 2) + &
                    StrSsym(1, 2) * StrSsym(2, 2) + &
                    StrSsym(1, 3) * StrSsym(3, 2) + &
                    StrAsym(1, 1) * StrAsym(1, 2) + &
                    StrAsym(1, 2) * StrAsym(2, 2) + &
                    StrAsym(1, 3) * StrAsym(3, 2)
                    S2plusO2(1, 3) = StrSsym(1, 1) * StrSsym(1, 3) + &
                    StrSsym(1, 2) * StrSsym(2, 3) + &
                    StrSsym(1, 3) * StrSsym(3, 3) + &
                    StrAsym(1, 1) * StrAsym(1, 3) + &
                    StrAsym(1, 2) * StrAsym(2, 3) + &
                    StrAsym(1, 3) * StrAsym(3, 3)

                    S2plusO2(2, 1) = StrSsym(2, 1) * StrSsym(1, 1) + &
                    StrSsym(2, 2) * StrSsym(2, 1) + &
                    StrSsym(2, 3) * StrSsym(3, 1) + &
                    StrAsym(2, 1) * StrAsym(1, 1) + &
                    StrAsym(2, 2) * StrAsym(2, 1) + &
                    StrAsym(2, 3) * StrAsym(3, 1)
                    S2plusO2(2, 2) = StrSsym(2, 1) * StrSsym(1, 2) + &
                    StrSsym(2, 2) * StrSsym(2, 2) + &
                    StrSsym(2, 3) * StrSsym(3, 2) + &
                    StrAsym(2, 1) * StrAsym(1, 2) + &
                    StrAsym(2, 2) * StrAsym(2, 2) + &
                    StrAsym(2, 3) * StrAsym(3, 2)
                    S2plusO2(2, 3) = StrSsym(2, 1) * StrSsym(1, 3) + &
                    StrSsym(2, 2) * StrSsym(2, 3) + &
                    StrSsym(2, 3) * StrSsym(3, 3) + &
                    StrAsym(2, 1) * StrAsym(1, 3) + &
                    StrAsym(2, 2) * StrAsym(2, 3) + &
                    StrAsym(2, 3) * StrAsym(3, 3)

                    S2plusO2(3, 1) = StrSsym(3, 1) * StrSsym(1, 1) + &
                    StrSsym(3, 2) * StrSsym(2, 1) + &
                    StrSsym(3, 3) * StrSsym(3, 1) + &
                    StrAsym(3, 1) * StrAsym(1, 1) + &
                    StrAsym(3, 2) * StrAsym(2, 1) + &
                    StrAsym(3, 3) * StrAsym(3, 1)
                    S2plusO2(3, 2) = StrSsym(3, 1) * StrSsym(1, 2) + &
                    StrSsym(3, 2) * StrSsym(2, 2) + &
                    StrSsym(3, 3) * StrSsym(3, 2) + &
                    StrAsym(3, 1) * StrAsym(1, 2) + &
                    StrAsym(3, 2) * StrAsym(2, 2) + &
                    StrAsym(3, 3) * StrAsym(3, 2)
                    S2plusO2(3, 3) = StrSsym(3, 1) * &
                    StrSsym(1, 3) + StrSsym(3, 2) * StrSsym(2, 3) + &
                    StrSsym(3, 3) * StrSsym(3, 3) + &
                    StrAsym(3, 1) * StrAsym(1, 3) + &
                    StrAsym(3, 2) * StrAsym(2, 3) + &
                    StrAsym(3, 3) * StrAsym(3, 3)

                    CALL DSYEVC3(S2plusO2,EIG)
                    IF (((EIG(1) - EIG(2)) * (EIG(1) - EIG(3))) < 0.0_WP) Lambda2_io(IC, JC, KC) = EIG(1)
                    IF (((EIG(2) - EIG(1)) * (EIG(2) - EIG(3))) < 0.0_WP) Lambda2_io(IC, JC, KC) = EIG(2)
                    IF (((EIG(3) - EIG(1)) * (EIG(3) - EIG(2))) < 0.0_WP) Lambda2_io(IC, JC, KC) = EIG(3)


                    !get rid of possible infinity
                    IF (Lambda2_io(IC, JC, KC) > 1.0E30_WP) Lambda2_io(IC, JC, KC) = 1.0E30_WP
                    IF (Lambda2_io(IC, JC, KC) < -1.0E30_WP) Lambda2_io(IC, JC, KC) = -1.0E30_WP
                    !***********************************************************************************
                END DO
            END DO
        END DO

        DEALLOCATE (dQ_io)
    END IF

    IF(TgFlowFlg)THEN
        allocate (Qcr_tg    (NND1_tg, NND2, NND3)       )
        allocate (Vor_tg    (NND1_tg, NND2, NND3, 4)     )
        allocate (Delta_tg  (NND1_tg, NND2, NND3)       )
        allocate (Lambda2_tg(NND1_tg, NND2, NND3)       )
        allocate (SWIrlStrength_tg(NND1_tg, NND2, NND3, 3) )
        allocate (dQ_tg     (NND1_tg, NND2, NND3, 3, 3)   )


        !====================GET ALL du_I / DX_i on nodeS ===================================================
        DO IC = 1, NCL1_tg
            IM = IMV_tg(IC)
            IP = IPV_tg(IC)
            DO JC = 2, NCL2
                JM = JGMV(JC)
                JP = JGPV(JC)
                DO KC = 1, NCL3
                    KP = KPV(KC)
                    KM = KMV(KC)
                    !===============================================================================
                    !=== DU / DX at i', j', k'============================
                    DUDX1 = (U_F0_tg(IC, JC, KM, 1) - U_F0_tg(IM, JC, KM, 1)) * DXI       ! (I - 1, j, K - 1)
                    DUDX2 = (U_F0_tg(IP, JC, KM, 1) - U_F0_tg(IC, JC, KM, 1)) * DXI       ! (i,   j, K - 1)

                    DUDX3 = (U_F0_tg(IC, JC, KC, 1) - U_F0_tg(IM, JC, KC, 1)) * DXI       ! (I - 1, j, k)
                    DUDX4 = (U_F0_tg(IP, JC, KC, 1) - U_F0_tg(IC, JC, KC, 1)) * DXI       ! (i,   j, k)

                    DUDX5 = (U_F0_tg(IC, JM, KM, 1) - U_F0_tg(IM, JM, KM, 1)) * DXI       ! (I - 1, J - 1, K - 1)
                    DUDX6 = (U_F0_tg(IP, JM, KM, 1) - U_F0_tg(IC, JM, KM, 1)) * DXI       ! (i,   J - 1, K - 1)

                    DUDX7 = (U_F0_tg(IC, JM, KC, 1) - U_F0_tg(IM, JM, KC, 1)) * DXI       ! (I - 1, J - 1, k)
                    DUDX8 = (U_F0_tg(IP, JM, KC, 1) - U_F0_tg(IC, JM, KC, 1)) * DXI       ! (i,   J - 1, k)

                    dQ_tg(IC, JC, KC, 1, 1) = ( (DUDX1 + DUDX2) * 0.5_WP + (DUDX3 + DUDX4) * 0.5_WP ) * 0.5_WP * YCL2ND_WFF(JC) + &
                    ( (DUDX3 + DUDX4) * 0.5_WP + (DUDX5+ DUDX6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== DU / Dy at i', J', K'============================
                    DUDY1 = (U_F0_tg(IC, JC, KM, 1) - U_F0_tg(IC, JM, KM, 1)) * DYCI(JC)  ! (i', J', K - 1)
                    DUDY2 = (U_F0_tg(IC, JC, KC, 1) - U_F0_tg(IC, JM, KC, 1)) * DYCI(JC)  ! (i', J', K)

                    dQ_tg(IC, JC, KC, 1, 2) = 0.5_WP * (DUDY1 + DUDY2)                 ! (i', J', K')

                    !=== DU / Dz at i', J', K'============================
                    DUDZ1 = (U_F0_tg(IC, JM, KC, 1) - U_F0_tg(IC, JM, KM, 1)) * DZI       ! (i', J - 1, K')
                    DUDZ2 = (U_F0_tg(IC, JC, KC, 1) - U_F0_tg(IC, JC, KM, 1)) * DZI       ! (i', J,  k')

                    dQ_tg(IC, JC, KC, 1, 3) = DUDZ1 * YCL2ND_WFB(JC) + DUDZ2 * YCL2ND_WFF(JC)

                    !===============================================================================
                    !=== Dv/ Dy at i', j', k'============================
                    IF(JC == NCL2) THEN
                        DVDY1 = (0.0_WP - U_F0_tg(IM, JC, KM, 2)) * DYFI(JC)       ! (I - 1, j, K - 1)
                        DVDY2 = (0.0_WP - U_F0_tg(IC, JC, KM, 2)) * DYFI(JC)       ! (i,   j, K - 1)

                        DVDY3 = (0.0_WP - U_F0_tg(IM, JC, KC, 2)) * DYFI(JC)       ! (I - 1, j, k)
                        DVDY4 = (0.0_WP - U_F0_tg(IC, JC, KC, 2)) * DYFI(JC)       ! (i,   j, k)
                    ELSE
                        DVDY1 = (U_F0_tg(IM, JP, KM, 2) - U_F0_tg(IM, JC, KM, 2)) * DYFI(JC)       ! (I - 1, j, K - 1)
                        DVDY2 = (U_F0_tg(IC, JP, KM, 2) - U_F0_tg(IC, JC, KM, 2)) * DYFI(JC)       ! (i,   j, K - 1)

                        DVDY3 = (U_F0_tg(IM, JP, KC, 2) - U_F0_tg(IM, JC, KC, 2)) * DYFI(JC)       ! (I - 1, j, k)
                        DVDY4 = (U_F0_tg(IC, JP, KC, 2) - U_F0_tg(IC, JC, KC, 2)) * DYFI(JC)       ! (i,   j, k)
                    END IF

                    DVDY5 = (U_F0_tg(IM, JC, KM, 2) - U_F0_tg(IM, JM, KM, 2)) * DYFI(JM)      ! (I - 1, J - 1, K - 1)
                    DVDY6 = (U_F0_tg(IC, JC, KM, 2) - U_F0_tg(IC, JM, KM, 2)) * DYFI(JM)      ! (i,   J - 1, K - 1)

                    DVDY7 = (U_F0_tg(IM, JC, KC, 2) - U_F0_tg(IM, JM, KC, 2)) * DYFI(JM)      ! (I - 1, J - 1, k)
                    DVDY8 = (U_F0_tg(IC, JC, KC, 2) - U_F0_tg(IC, JM, KC, 2)) * DYFI(JM)      ! (i,   J - 1, k)

                    dQ_tg(IC, JC, KC, 2, 2) = ( (DVDY1 + DVDY2) * 0.5_WP + (DVDY3 + DVDY4) * 0.5_WP ) * 0.5_WP * YCL2ND_WFF(JC) + &
                    ( (DVDY3 + DVDY4) * 0.5_WP + (DVDY5+ DVDY6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== Dv/ DX at i', J', K'============================
                    DVDX1 = (U_F0_tg(IC, JC, KM, 2) - U_F0_tg(IM, JC, KM, 2)) * DXI  ! (i', J', K - 1)
                    DVDX2 = (U_F0_tg(IC, JC, KC, 2) - U_F0_tg(IM, JC, KC, 2)) * DXI  ! (i', J', K)

                    dQ_tg(IC, JC, KC, 2, 1) = 0.5_WP * (DVDX1 + DVDX2)                 ! (i', J', K')

                    !=== Dv/ Dz at i', J', K'============================
                    DVDZ1 = (U_F0_tg(IM, JC, KC, 2) - U_F0_tg(IM, JC, KM, 2)) * DZI       ! (I - 1, J', K')
                    DVDZ2 = (U_F0_tg(IC, JC, KC, 2) - U_F0_tg(IC, JC, KM, 2)) * DZI       ! (i, J', K')

                    dQ_tg(IC, JC, KC, 2, 3) = (DVDZ1 + DVDZ2) * 0.5_WP

                    !===============================================================================
                    !=== Dw/ Dz at i', j', k'============================
                    DWDZ1 = (U_F0_tg(IM, JC, KC, 3) - U_F0_tg(IM, JC, KM, 3)) * DXI       ! (I - 1, j, K - 1)
                    DWDZ2 = (U_F0_tg(IC, JC, KC, 3) - U_F0_tg(IC, JC, KM, 3)) * DXI       ! (i,   j, K - 1)

                    DWDZ3 = (U_F0_tg(IM, JC, KP, 3) - U_F0_tg(IM, JC, KC, 3)) * DXI       ! (I - 1, j, k)
                    DWDZ4 = (U_F0_tg(IC, JC, KP, 3) - U_F0_tg(IC, JC, KC, 3)) * DXI       ! (i,   j, k)

                    DWDZ5 = (U_F0_tg(IM, JM, KC, 3) - U_F0_tg(IM, JM, KM, 3)) * DXI       ! (I - 1, J - 1, K - 1)
                    DWDZ6 = (U_F0_tg(IC, JM, KC, 3) - U_F0_tg(IC, JM, KM, 3)) * DXI       ! (i,   J - 1, K - 1)

                    DWDZ7 = (U_F0_tg(IM, JM, KP, 3) - U_F0_tg(IM, JM, KC, 3)) * DXI       ! (I - 1, J - 1, k)
                    DWDZ8 = (U_F0_tg(IC, JM, KP, 3) - U_F0_tg(IC, JM, KC, 3)) * DXI       ! (i,   J - 1, k)

                    dQ_tg(IC, JC, KC, 3, 3) = ( (DWDZ1 + DWDZ2) * 0.5_WP + (DWDZ3 + DWDZ4) * 0.5_WP ) * 0.5_WP * YCL2ND_WFF(JC) + &
                    ( (DWDZ3 + DWDZ4) * 0.5_WP + (DWDZ5+ DWDZ6) * 0.5_WP ) * 0.5_WP * YCL2ND_WFB(JC)

                    !=== Dw/ Dy at i', J', K'============================
                    DWDY1 = (U_F0_tg(IM, JC, KC, 3) - U_F0_tg(IM, JM, KC, 3)) * DYCI(JC)  ! (I - 1, J', K')
                    DWDY2 = (U_F0_tg(IC, JC, KC, 3) - U_F0_tg(IC, JM, KC, 3)) * DYCI(JC)  ! (i, J', K)

                    dQ_tg(IC, JC, KC, 3, 2) = 0.5_WP * (DWDY1 + DWDY2)                 ! (i', J', K')

                    !=== Dw/ DX at i', J', K'============================
                    DWDX1 = (U_F0_tg(IC, JM, KC, 3) - U_F0_tg(IM, JM, KC, 3)) * DZI       ! (i', J - 1, K')
                    DWDX2 = (U_F0_tg(IC, JC, KC, 3) - U_F0_tg(IM, JC, KC, 3)) * DZI       ! (i', J,  k')

                    dQ_tg(IC, JC, KC, 3, 1) = DWDX1 * YCL2ND_WFB(JC) + DWDX2 * YCL2ND_WFF(JC)
                END DO
            END DO
        END DO

        dQ_tg(:, NND2, :, :, :) = 2.0_WP * DQ_tg(:, NCL2, :, :, :) - DQ_tg(:, NCL2 - 1, :, :, :)
        dQ_tg(:, 1,   :, :, :) = 2.0_WP * DQ_tg(:, 2,   :, :, :) - DQ_tg(:, 3,     :, :, :)
        dQ_tg(NND1_tg, :,:,   :, :) = dQ_tg(1, :, :, :, :)
        dQ_tg(:,      :, NND3, :, :) = dQ_tg(:, :, 1, :, :)


        DO IC = 1, NND1_tg
            DO JC = 1, NND2
                DO KC = 1, NND3
                    !=================== VortICitY ===============================
                    Vor_tg(IC, JC, KC, 1) = dQ_tg(IC, JC, KC, 3, 2) - dQ_tg(IC, JC, KC, 2, 3) ! Vor_tg-x
                    Vor_tg(IC, JC, KC, 2) = dQ_tg(IC, JC, KC, 1, 3) - dQ_tg(IC, JC, KC, 3, 1) ! Vor_tg- Y
                    Vor_tg(IC, JC, KC, 3) = dQ_tg(IC, JC, KC, 2, 1) - dQ_tg(IC, JC, KC, 1, 2) ! Vor_tg-z
                    Vor_tg(IC, JC, KC, 4) = DSQRT( DOT_PRODUCT( Vor_tg(IC, JC, KC, 1:3), Vor_tg(IC, JC, KC, 1:3) ) )!Vor_tg-mag

                    !==========================================================================================
                    !===================omega = antISymetrIC components of graduate U===============
                    StrAsym(1, 1) = dQ_tg(IC, JC, KC, 1, 1) - dQ_tg(IC, JC, KC, 1, 1)  !dU / DX - dU / DX
                    StrAsym(1, 2) = dQ_tg(IC, JC, KC, 2, 1) - dQ_tg(IC, JC, KC, 1, 2)  !dV/ DX - dU / DY
                    StrAsym(1, 3) = dQ_tg(IC, JC, KC, 3, 1) - dQ_tg(IC, JC, KC, 1, 3)  !dW/ DX - dU / DZ

                    StrAsym(2, 1) = dQ_tg(IC, JC, KC, 1, 2) - dQ_tg(IC, JC, KC, 2, 1)  !dU / DY - dV/ DX
                    StrAsym(2, 2) = dQ_tg(IC, JC, KC, 2, 2) - dQ_tg(IC, JC, KC, 2, 2)  !dV/ DY - dV/ DY
                    StrAsym(2, 3) = dQ_tg(IC, JC, KC, 3, 2) - dQ_tg(IC, JC, KC, 2, 3)  !dW/ DY - dV/ DZ

                    StrAsym(3, 1) = dQ_tg(IC, JC, KC, 1, 3) - dQ_tg(IC, JC, KC, 3, 1)  !dU / DZ - dW/ DX
                    StrAsym(3, 2) = dQ_tg(IC, JC, KC, 2, 3) - dQ_tg(IC, JC, KC, 3, 2)  !dV/ DZ - dW/ DY
                    StrAsym(3, 3) = dQ_tg(IC, JC, KC, 3, 3) - dQ_tg(IC, JC, KC, 3, 3)  !dW/ DZ - dW/ DZ
                    !===================omega = symetrIC components of graduate U===============
                    StrSsym(1, 1) = dQ_tg(IC, JC, KC, 1, 1) + dQ_tg(IC, JC, KC, 1, 1)  !dU / DX + dU / DX
                    StrSsym(1, 2) = dQ_tg(IC, JC, KC, 2, 1) + dQ_tg(IC, JC, KC, 1, 2)  !dV/ DX + dU / DY
                    StrSsym(1, 3) = dQ_tg(IC, JC, KC, 3, 1) + dQ_tg(IC, JC, KC, 1, 3)  !dW/ DX + dU / DZ

                    StrSsym(2, 1) = dQ_tg(IC, JC, KC, 1, 2) + dQ_tg(IC, JC, KC, 2, 1)  !dU / DY + dV/ DX
                    StrSsym(2, 2) = dQ_tg(IC, JC, KC, 2, 2) + dQ_tg(IC, JC, KC, 2, 2)  !dV/ DY + dV/ DY
                    StrSsym(2, 3) = dQ_tg(IC, JC, KC, 3, 2) + dQ_tg(IC, JC, KC, 2, 3)  !dW/ DY + dV/ DZ

                    StrSsym(3, 1) = dQ_tg(IC, JC, KC, 1, 3) + dQ_tg(IC, JC, KC, 3, 1)  !dU / DZ + dW/ DX
                    StrSsym(3, 2) = dQ_tg(IC, JC, KC, 2, 3) + dQ_tg(IC, JC, KC, 3, 2)  !dV/ DZ + dW/ DY
                    StrSsym(3, 3) = dQ_tg(IC, JC, KC, 3, 3) + dQ_tg(IC, JC, KC, 3, 3)  !dW/ DZ + dW/ DZ

                    StrAsym  = StrAsym  * 0.5_WP  ! Sij
                    StrSsym  = StrSsym  * 0.5_WP  ! WIj

                    ! trace of (STSYM *STSYM^T)
                    STSSYM = DOT_PRODUCT(StrSsym(:, 1), StrSsym(:, 1)) + &
                    DOT_PRODUCT(StrSsym(:, 2), StrSsym(:, 2)) + &
                    DOT_PRODUCT(StrSsym(:, 3), StrSsym(:, 3))    ! SiJ *Sij
                    ! trace of (STASYM *STASYM^T)
                    STASYM = DOT_PRODUCT(StrAsym(:, 1), StrAsym(:, 1)) + &
                    DOT_PRODUCT(StrAsym(:, 2), StrAsym(:, 2)) + &
                    DOT_PRODUCT(StrAsym(:, 3), StrAsym(:, 3))    ! WIJ * WIj

                    !================================ Q cretiria===================================================
                    QCR_tg(IC, JC, KC) = 0.5_WP * (STASYM -STSSYM)


                    !================================ Delta_io criteria===============================================
                    ! grad U
                    gradU = StrAsym + StrSsym

                    !  negtive trace of velocity gradient  !P
                    TraceGradU = (StrAsym(1, 1) +StrAsym(2, 2) +StrAsym(3, 3) +&
                    StrSsym(1, 1) +StrSsym(2, 2) +StrSsym(3, 3)) * (-1.0_WP)

                    ! det(gradU) !R
                    detGradU = (gradU(1, 1) * ( gradU(2, 2) * GradU(3, 3) - gradU(2, 3) * GradU(3, 2) ) - &
                    gradU(1, 2) * ( gradU(2, 1) * GradU(3, 3) - gradU(3, 1) * GradU(2, 3) ) + &
                    gradU(1, 3) * ( gradU(2, 1) * GradU(3, 2) - gradU(2, 2) * GradU(3, 1) )) * (-1.0_WP)
                    ! Delta for incompressible flow
                    !Delta_io(IC, JC, KC) = (QCR_io(IC, JC, KC) / 3.0_WP)**3 + (detGradU / 2.0_WP)**2

                    ! Delta for compressible flow
                    Delta_tg(IC, JC, KC) = 27.0_WP * DetGradU* DetGradU + &
                    (4.0_WP * TraceGradU* TraceGradU* TraceGradU- 18.0_WP * TraceGradU* QCR_tg(IC, JC, KC)) * DetGradU + &
                    4.0_WP * QCR_tg(IC, JC, KC) * QCR_tg(IC, JC, KC) * QCR_tg(IC, JC, KC) - &
                    TraceGradU* TraceGradU* QCR_tg(IC, JC, KC) * QCR_tg(IC, JC, KC)

                    !get rid of possible infinity
                    IF (Delta_tg(IC, JC, KC) > 1.0E30_WP) Delta_tg(IC, JC, KC) = 1.0E30_WP
                    IF (Delta_tg(IC, JC, KC) < -1.0E30_WP) Delta_tg(IC, JC, KC) = -1.0E30_WP

                    !============== SWIrl stRENgth ==================================
                    IF (Delta_tg(IC, JC, KC) > 0.0_WP) THEN

                        MATZ = 0
                        NM = 3
                        NN = 3
                        AR(:, :) = StrAsym(:, :) + StrSsym(:, :)
                        AI(:, :) = 0.0

                        CALL cg(NM, NN, AR, AI,WR, WI, MATZ, ZR, ZI, FV1, FV2, FV3, IERR)

                        IF(DABS(WI(1)) < (DABS(WI(2)) + 1.0E-14) .AND. DABS(WI(1)) < (DABS(WI(3)) + 1.0E-14) ) THEN
                            SWIrlStrength_tg(IC, JC, KC, 1) = WR(1)
                            SWIrlStrength_tg(IC, JC, KC, 2) = WR(2)
                            SWIrlStrength_tg(IC, JC, KC, 3) = WR(2)
                        ELSE IF (DABS(WI(2)) < (DABS(WI(1)) + 1.0E-14) .AND. DABS(WI(2)) < (DABS(WI(3)) + 1.0E-14) ) THEN
                            SWIrlStrength_tg(IC, JC, KC, 1) = WR(2)
                            SWIrlStrength_tg(IC, JC, KC, 2) = WR(1)
                            SWIrlStrength_tg(IC, JC, KC, 3) = WR(1)
                        ELSE IF (DABS(WI(3)) < (DABS(WI(1)) + 1.0E-14) .AND. DABS(WI(3)) < (DABS(WI(2)) + 1.0E-14) ) THEN
                            SWIrlStrength_tg(IC, JC, KC, 1) = WR(3)
                            SWIrlStrength_tg(IC, JC, KC, 2) = WR(1)
                            SWIrlStrength_tg(IC, JC, KC, 3) = WR(1)
                        ELSE
                        END IF

                    END IF



                    !************* Lamda2 criteria*************************
                    S2plusO2(1, 1) = StrSsym(1, 1) * StrSsym(1, 1) + &
                    StrSsym(1, 2) * StrSsym(2, 1) + &
                    StrSsym(1, 3) * StrSsym(3, 1) + &
                    StrAsym(1, 1) * StrAsym(1, 1) + &
                    StrAsym(1, 2) * StrAsym(2, 1) + &
                    StrAsym(1, 3) * StrAsym(3, 1)
                    S2plusO2(1, 2) = StrSsym(1, 1) * StrSsym(1, 2) + &
                    StrSsym(1, 2) * StrSsym(2, 2) + &
                    StrSsym(1, 3) * StrSsym(3, 2) + &
                    StrAsym(1, 1) * StrAsym(1, 2) + &
                    StrAsym(1, 2) * StrAsym(2, 2) + &
                    StrAsym(1, 3) * StrAsym(3, 2)
                    S2plusO2(1, 3) = StrSsym(1, 1) * StrSsym(1, 3) + &
                    StrSsym(1, 2) * StrSsym(2, 3) + &
                    StrSsym(1, 3) * StrSsym(3, 3) + &
                    StrAsym(1, 1) * StrAsym(1, 3) + &
                    StrAsym(1, 2) * StrAsym(2, 3) + &
                    StrAsym(1, 3) * StrAsym(3, 3)

                    S2plusO2(2, 1) = StrSsym(2, 1) * StrSsym(1, 1) + &
                    StrSsym(2, 2) * StrSsym(2, 1) + &
                    StrSsym(2, 3) * StrSsym(3, 1) + &
                    StrAsym(2, 1) * StrAsym(1, 1) + &
                    StrAsym(2, 2) * StrAsym(2, 1) + &
                    StrAsym(2, 3) * StrAsym(3, 1)
                    S2plusO2(2, 2) = StrSsym(2, 1) * StrSsym(1, 2) + &
                    StrSsym(2, 2) * StrSsym(2, 2) + &
                    StrSsym(2, 3) * StrSsym(3, 2) + &
                    StrAsym(2, 1) * StrAsym(1, 2) + &
                    StrAsym(2, 2) * StrAsym(2, 2) + &
                    StrAsym(2, 3) * StrAsym(3, 2)
                    S2plusO2(2, 3) = StrSsym(2, 1) * StrSsym(1, 3) + &
                    StrSsym(2, 2) * StrSsym(2, 3) + &
                    StrSsym(2, 3) * StrSsym(3, 3) + &
                    StrAsym(2, 1) * StrAsym(1, 3) + &
                    StrAsym(2, 2) * StrAsym(2, 3) + &
                    StrAsym(2, 3) * StrAsym(3, 3)

                    S2plusO2(3, 1) = StrSsym(3, 1) * StrSsym(1, 1) + &
                    StrSsym(3, 2) * StrSsym(2, 1) + &
                    StrSsym(3, 3) * StrSsym(3, 1) + &
                    StrAsym(3, 1) * StrAsym(1, 1) + &
                    StrAsym(3, 2) * StrAsym(2, 1) + &
                    StrAsym(3, 3) * StrAsym(3, 1)
                    S2plusO2(3, 2) = StrSsym(3, 1) * StrSsym(1, 2) + &
                    StrSsym(3, 2) * StrSsym(2, 2) + &
                    StrSsym(3, 3) * StrSsym(3, 2) + &
                    StrAsym(3, 1) * StrAsym(1, 2) + &
                    StrAsym(3, 2) * StrAsym(2, 2) + &
                    StrAsym(3, 3) * StrAsym(3, 2)
                    S2plusO2(3, 3) = StrSsym(3, 1) * StrSsym(1, 3) + &
                    StrSsym(3, 2) * StrSsym(2, 3) + &
                    StrSsym(3, 3) * StrSsym(3, 3) + &
                    StrAsym(3, 1) * StrAsym(1, 3) + &
                    StrAsym(3, 2) * StrAsym(2, 3) + StrAsym(3, 3) * StrAsym(3, 3)

                    CALL DSYEVC3(S2plusO2,EIG)
                    IF (((EIG(1) - EIG(2)) * (EIG(1) - EIG(3))) < 0.0_WP) Lambda2_tg(IC, JC, KC) = EIG(1)
                    IF (((EIG(2) - EIG(1)) * (EIG(2) - EIG(3))) < 0.0_WP) Lambda2_tg(IC, JC, KC) = EIG(2)
                    IF (((EIG(3) - EIG(1)) * (EIG(3) - EIG(2))) < 0.0_WP) Lambda2_tg(IC, JC, KC) = EIG(3)


                    !get rid of possible infinity
                    IF (Lambda2_tg(IC, JC, KC) > 1.0E30_WP) Lambda2_tg(IC, JC, KC) = 1.0E30_WP
                    IF (Lambda2_tg(IC, JC, KC) < -1.0E30_WP) Lambda2_tg(IC, JC, KC) = -1.0E30_WP
                    !***********************************************************************************
                END DO
            END DO
        END DO
        DEALLOCATE(dQ_tg)
    END IF



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
!SUBROUTINE TEC360_PROFILE_TEST
!        USE QPGATHERING_INFO
!        USE TEC360_INFO
!        USE mesh_info
!        USE flow_info
!        USE init_info
!        USE thermal_info
!        IMPLICIT NONE

!        INTEGER(4) :: I, J, K
!        INTEGER(4) :: TECFLG = 202

!        LOGICAL :: File_exists

!        INQUIRE(FILE = 'RESULT.Centerline.plt', exist =File_exists)

!        IF(File_exists) THEN
!            OPEN(TECFLG, FILE = 'RESULT.Centerline.plt', POSITION = 'APPEND')
!        ELSE
!            OPEN(TECFLG, FILE = 'RESULT.Centerline.plt')
!        END IF

!        WRITE(TECFLG, '(A)') 'TITLE = "DNS Centerline"'
!        WRITE(TECFLG, '(A)') 'variables = "X", "U", "V", "W", "P"'

!        WRITE(TECFLG, '(A, 1I6.1, 1ES13.5, A)') 'ZONE T = " ', ITERG,PhyTIME, '"'

!        J = NCL2 / 2
!        K = NCL3 / 2

!        DO I = 0, NCL1_io
!           WRITE(TECFLG, '(2ES15.7)') 0.5_WP * (XND_io(I) + XND_io(I + 1)), P_F0_io(I, J, K)
!        END DO
!        close(tecflg)

!        RETURN
!     END SUBROUTINE

!    !***********************************************************************************************
!SUBROUTINE BULK_PARAMETER_CHECK
!    USE TEC360_INFO
!    USE mesh_info
!    USE flow_info
!    USE init_info
!    IMPLICIT NONE

!    INTEGER(4) :: IM, IC, JC, KC
!    INTEGER(4) :: TECFLG = 202

!    REAL(WP) :: MASS_RATE(NCL1_io)
!    REAL(WP) :: MASS_FLUX(NCL1_io)

!    REAL(WP) :: ENTH_RATE
!    REAL(WP) :: ENTH_FLUX(NCL1_io)

!    REAL(WP) :: RHOU, RHO, RHOU_Z
!    REAL(WP) :: DH, DHU,DHU_Z


!    !=============BULK MASS FLUX ========================
!    DO IC = 1, NCL1_io
!        IM = IMV_io(IC)

!        MASS_RATE(IC) = 0.0_WP

!        DO JC = 1, NCL2
!            !==== AVERAGE IN HOMOGENOUS direction
!            RHOU = 0.0_WP
!            DO KC = 1, NCL3
!                RHO  = D_F0_io(IM, JC, KC) + D_F0_io(IC, JC, KC)
!                RHOU = RHOU+ U_F0_io(IC, JC, KC) * RHO
!            END DO
!            RHOU_Z = RHOU / NCL3 / 2.0_WP

!            MASS_RATE(IC) = MASS_RATE(IC) + RHOU_Z* DYFI(JC) / DZI
!        END DO
!        MASS_FLUX(IC) = MASS_RATE(IC) /Area_inlet

!    END DO


!    !=============BULK ENTHALPY ========================
!    DO IC = 1, NCL1_io
!        IM = IMV_io(IC)

!        ENTH_RATE = 0.0_WP

!        DO JC = 1, NCL2
!            !==== AVERAGE IN HOMOGENOUS direction
!            DHU = 0.0_WP
!            DO KC = 1, NCL3
!                DH  = D_F0_io(IM, JC, KC) * H_F0_io(IM, JC, KC) + &
!                        D_F0_io(IC, JC, KC) * H_F0_io(IC, JC, KC)

!                DHU = DHU+ U_F0_io(IC, JC, KC) * DH
!            END DO
!            DHU_Z = DHU / NCL3 / 2.0_WP

!            ENTH_RATE = ENTH_RATE + DHU_Z* DYFI(JC) / DZI
!        END DO
!        ENTH_FLUX(IC) = ENTH_RATE/ MASS_RATE(IC)
!    END DO

!    !===============WRITE OUT ================================

!    IF(NCOUNT(4) == 0) THEN
!       OPEN(TECFLG, FILE = 'CHECK_CONSERVATION.plt')
!    ELSE
!       OPEN(TECFLG, FILE = 'CHECK_CONSERVATION.plt', POSITION = 'APPEND')
!    END IF
!    NCOUNT(4) = NCOUNT(4) + 1
!    WRITE(TECFLG, '(A)') 'TITLE = "DNS FLOW CONSERVATION"'
!    WRITE(TECFLG, '(A)') 'variables = "X", "Gb", "Hb"'
!    WRITE(TECFLG, '(A, 1I6.1, 1ES13.5)') 'ZONE T = " ', ITERG,PhyTIME

!    DO IC = 1, NCL1
!       WRITE(*, *) XND_io(IC), MASS_FLUX(IC),ENTH_FLUX(IC)
!    END DO

!    RETURN
!    END SUBROUTINE





!===============CONSERVATION CHECK =============================================================================
!         !=============BULK MASS FLUX ========================
!        DO IC = 1, NCL1_io
!            IM = IMV_io(IC)

!            MASS_RATE(IC) = 0.0_WP

!            DO JC = 1, NCL2
!                !==== AVERAGE IN HOMOGENOUS direction
!                RHOUt = 0.0_WP
!                DO KC = 1, NCL3
!                    RHOt  = D_F0_io(IM, JC, KC) + D_F0_io(IC, JC, KC)
!                    RHOUt = RHOUT + U_F0_io(IC, JC, KC) * RHOt
!                END DO
!                RHOU_Z = RHOUt/ NCL3 / 2.0_WP

!                MASS_RATE(IC) = MASS_RATE(IC) + RHOU_Z/ DYFI(JC) / DZI
!            END DO
!            MASS_FLUX(IC) = MASS_RATE(IC) /Area_inlet

!        END DO


!        !=============BULK ENTHALPY ========================
!        DO IC = 1, NCL1_io
!            IM = IMV_io(IC)

!            ENTH_RATE = 0.0_WP

!            DO JC = 1, NCL2
!                !==== AVERAGE IN HOMOGENOUS direction
!                DHUt = 0.0_WP
!                DO KC = 1, NCL3
!                    DHt  = D_F0_io(IM, JC, KC) * H_F0_io(IM, JC, KC) + &
!                             D_F0_io(IC, JC, KC) * H_F0_io(IC, JC, KC)

!                    DHUt = DHUT + U_F0_io(IC, JC, KC) * DHt
!                END DO
!                DHU_Z = DHUt/ NCL3 / 2.0_WP

!                ENTH_RATE = ENTH_RATE + DHU_Z/ DYFI(JC) / DZI
!            END DO
!            ENTH_FLUX(IC) = ENTH_RATE/ MASS_RATE(IC)
!        END DO

!        !===============WRITE OUT ================================

!        IF(NCOUNT(5) == 0) THEN
!           OPEN(TECFLG, FILE = 'CHECK_CONSERVATION.plt')
!        ELSE
!           OPEN(TECFLG, FILE = 'CHECK_CONSERVATION.plt', POSITION = 'APPEND')
!        END IF
!        NCOUNT(5) = NCOUNT(5) + 1
!        WRITE(TECFLG, '(A)') 'TITLE = "DNS FLOW CONSERVATION"'
!        WRITE(TECFLG, '(A)') 'variables = "X", "Gb", "Hb"'
!        WRITE(TECFLG, '(A, 1I6.1, 1ES13.5)') 'ZONE T = " ', ITERG,PhyTIME

!        DO IC = 1, NCL1_io
!           WRITE(TECFLG, '(3ES15.7)') XND_io(IC), MASS_FLUX(IC),ENTH_FLUX(IC)
!        END DO


!<=================================================================
!=============================== tg =========================

!=========================================================================================
!SUBROUTINE QPGATHERING_TG
!!>    @WARning: 1) MPI_GATHER requires the data account from each pARtition
!!>             sould be the same.
!        USE QPGATHERING_INFO
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        USE postprocess_info
!        IMPLICIT NONE

!        INTEGER(4) :: INN1, INN2, INN3
!        INTEGER(4) :: I,L
!        INTEGER(4) :: J, JJ, JJP, JP
!        INTEGER(4) :: K, KK, KS
!        INTEGER(4) :: N2DOID


!        IF( .NOT.TgFlowFlg) RETURN

!        !================== ALLOCATTE DATA=============================
!        ALLOCATE ( U_F0_tg(NCL1_tg, 0 : NND2, NCL3) )
!        ALLOCATE ( U_F0_tg(NCL1_tg, 0 : NND2, NCL3) )
!        ALLOCATE ( U_F0_tg(NCL1_tg, 0 : NND2, NCL3) )
!        ALLOCATE ( P_F0_tg(NCL1_tg, NCL2, NCL3) )

!        ALLOCATE ( QAUX_tg(NCL1_tg, 0 : N2DO(0) + 1, NCL3, NDV, 1 : NPTOT) )
!        ALLOCATE ( PAUX_tg(NCL1_tg, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )

!        ALLOCATE ( D1AUX (N2DO(MYID), NDV + 1, 1 : NPTOT) )
!        ALLOCATE ( U1xzL_F0_tg( NCL2, NDV + 1 ) )



!        !=============gather data to masteR ===========================
!        INN1 = NCL1_tg * (N2DO(0) + 2) * NCL3 * NDV
!        CALL MPI_GATHER( Q_tg, INN1, MPI_DOUBLE_PRECISION, QAUX_tg, INN1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM, IERROR)

!        INN2 = NCL1_tg * (N2DO(0) + 2) * NCL3
!        CALL MPI_GATHER( PR_tg, INN2, MPI_DOUBLE_PRECISION, PAUX_tg, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM, IERROR)

!        INN3 = N2DO(MYID) * (NDV + 1)
!        CALL MPI_GATHER( U1xzL_tg, INN3, MPI_DOUBLE_PRECISION, D1AUX, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM, IERROR)


!        !============== RE - ARrange gathered data========================
!        IF (MYID == 0) THEN

!            DO KK = 0, NPSLV
!                N2DOID=JDEWT(KK) - JDSWT(KK) + 1
!                DO J = 1, N2DOID
!                    JJ = JDSWT(KK) - 1 + J
!                    DO L = 1, NDV + 1
!                        U1xzL_F0_tg(JJ, L) = D1AUX(J, L, KK+ 1)
!                    END DO
!                END DO
!            END DO

!            DO KK = 0, NPSLV
!                N2DOID=JDEWT(KK) - JDSWT(KK) + 1
!                DO  J = 1, N2DOID
!                    JJ = JDSWT(KK) - 1 + J
!                    JP = J + 1
!                    JJP = JJ + 1
!                    DO I = 1, NCL1_tg
!                        DO K = 1, NCL3
!                            KS = KSYM(K)

!                            U_F0_tg(I, JJ, K) = (QAUX_tg(I, J, K, 1, KK+ 1))
!                            U_F0_tg(I, JJ, K) = (QAUX_tg(I, J, K, 3, KK+ 1) * RCCI1(JJ))
!                            IF(JJ == 1 .AND. iCase == iPIPEC) THEN
!                                U_F0_tg(I, JJ, K) = (QAUX_tg(I, JP, K, 2, KK+ 1) - QAUX_tg(I, JP, KS, 2, KK+ 1) ) * 0.50_WP * RNDI1(JJP)
!                            ELSE
!                                U_F0_tg(I, JJ, K) = (QAUX_tg(I, J, K, 2, KK+ 1) * RNDI1(JJ))
!                            END IF
!                            P_F0_tg(I, JJ, K) = (PAUX_tg(I, J, K, KK+ 1))
!                        END DO
!                    END DO
!                END DO
!            END DO

!            DO I = 1, NCL1_tg
!                DO K = 1, NCL3
!                    U_F0_tg(I, NND2, K) = (QAUX_tg(I, N2DO(MYID) + 1, K, 2, NPTOT))
!                END DO
!            END DO

!            DO I = 1, NCL1_tg
!                DO K = 1, NCL3
!                    U_F0_tg(I, 0,   K) = U_F0_tg(I, 1,   K) * (-1.0_WP)
!                    U_F0_tg(I, NND2, K) = U_F0_tg(I, NND2, K) * (-1.0_WP)
!                    U_F0_tg(I, 0,   K) = U_F0_tg(I, 1,   K) * (-1.0_WP)
!                    U_F0_tg(I, NND2, K) = U_F0_tg(I, NND2, K) * (-1.0_WP)
!                    U_F0_tg(I, 0,   K) = 0.0_WP
!                END DO
!            END DO

!        END IF

!        RETURN
!    END SUBROUTINE

!!========================= iO===============================================
!SUBROUTINE QPGATHERING_io
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        USE THERMAL_INFO
!        USE QPGATHERING_INFO
!        USE postprocess_info
!        IMPLICIT NONE

!        INTEGER(4) :: INN1, INN2, INN3
!        INTEGER(4) :: I, L
!        INTEGER(4) :: J, JJ, JJP, JP
!        INTEGER(4) :: K, KK, KS
!        INTEGER(4) :: N2DOID



!        IF( .NOT.IoFlowFlg) RETURN

!        !================== ALLOCATTE DATA=============================

!        ALLOCATE ( U_F0_io(NCL1S : NCL1E, 0 : NND2, NCL3) )
!        ALLOCATE ( U_F0_io(NCL1S : NCL1E, 0 : NND2, NCL3) )
!        ALLOCATE ( U_F0_io(NCL1S : NCL1E, 0 : NND2, NCL3) )
!        ALLOCATE ( QAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3, NDV, 1 : NPTOT) )
!        ALLOCATE ( PAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!        ALLOCATE ( P_F0_io(NCL1S : NCL1E, NCL2, NCL3) )

!        IF(iThermoDynamics == 1) THEN
!            ALLOCATE ( T_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
!            ALLOCATE ( D_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
!            ALLOCATE ( H_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
!            ALLOCATE ( K_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
!            ALLOCATE ( M_F0_io(NCL1S : NCL1E, NCL2, NCL3) )
!            ALLOCATE ( TAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!            ALLOCATE ( DAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!            ALLOCATE ( HAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!            ALLOCATE ( KAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!            ALLOCATE ( MAUX_io(NCL1S : NCL1E, 0 : N2DO(0) + 1, NCL3,    1 : NPTOT) )
!        END IF

!        ALLOCATE ( D1AUXX (NCL1_io, N2DO(MYID), NDV + 1, 1 : NPTOT) )
!        ALLOCATE ( U1zL_F0_io(NCL1_io, NCL2, NDV + 1 ) )

!        ALLOCATE ( D1AUXZ (N2DO(MYID), NDV + 1, 1 : NPTOT) )
!        ALLOCATE ( U1xzL_F0_io(NCL2, NDV + 1 ) )

!        !IF(MYID == 0) THEN
!          !WRITE(*, *) 'U, MYID', Q_io(Ncl1_tg / 2, :, Ncl3 / 2, 1)
!        !END IF

!        !=============gather data to masteR ===========================
!        INN1 = (NCL1E - NCL1S + 1) * (N2DO(0) + 2) * NCL3 * NDV
!        CALL MPI_GATHER( Q_io, INN1, MPI_DOUBLE_PRECISION, QAUX_io, INN1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM, IERROR)

!        INN2 = (NCL1E - NCL1S + 1) * (N2DO(0) + 2) * NCL3
!        CALL MPI_GATHER( PR_io, INN2, MPI_DOUBLE_PRECISION, PAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM, IERROR)

!        IF(iThermoDynamics == 1) THEN
!            CALL MPI_GATHER( TEMPERATURE, INN2, MPI_DOUBLE_PRECISION, TAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)

!            CALL MPI_GATHER( DENSITY, INN2, MPI_DOUBLE_PRECISION, DAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)

!            CALL MPI_GATHER( ENTHALPY, INN2, MPI_DOUBLE_PRECISION, HAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)

!            CALL MPI_GATHER( Viscousity, INN2, MPI_DOUBLE_PRECISION, MAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)

!            CALL MPI_GATHER( THERMCONDT, INN2, MPI_DOUBLE_PRECISION, KAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)
!        END IF

!        IF(TgFlowFlg) THEN
!            INN3  = NCL1_io * N2DO(MYID) * (NDV + 1)
!            CALL MPI_GATHER( U1zL_io, INN3, MPI_DOUBLE_PRECISION, D1AUXX, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)
!        ELSE
!            INN3 = N2DO(MYID) * (NDV + 1)
!            CALL MPI_GATHER( U1xzL_io, INN3, MPI_DOUBLE_PRECISION, D1AUXZ, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!            CALL MPI_BARRIER(ICOMM, IERROR)
!        END IF

!        !============== RE - ARrange gathered data=============================
!        IF (MYID == 0) THEN

!            IF(TgFlowFlg) THEN
!                DO KK = 0, NPSLV
!                    N2DOID=JDEWT(KK) - JDSWT(KK) + 1
!                    DO I = 1, NCL1_io
!                        DO J = 1, N2DOID
!                            JJ = JDSWT(KK) - 1 + J
!                            DO L = 1, NDV + 1
!                                !U1xtL_F0_io(I, JJ, L) = D1AUXX(I, J, L, KK+ 1)
!                            END DO
!                        END DO
!                    END DO
!                END DO
!            ELSE
!                DO KK = 0, NPSLV
!                    N2DOID=JDEWT(KK) - JDSWT(KK) + 1
!                    DO J = 1, N2DOID
!                        JJ = JDSWT(KK) - 1 + J
!                        DO L = 1, NDV + 1
!                            U1xzL_F0_io(JJ, L) = D1AUXZ(J, L, KK+ 1)
!                        END DO
!                    END DO
!                END DO
!            END IF



!            DO KK = 0, NPSLV
!                N2DOID=JDEWT(KK) - JDSWT(KK) + 1
!                DO J = 1, N2DOID
!                    JJ = JDSWT(KK) - 1 + J
!                    JP = J + 1
!                    JJP = JJ + 1
!                    DO I = NCL1S, NCL1E
!                        DO K = 1, NCL3
!                            KS = KSYM(K)
!                            U_F0_io(I, JJ, K) = QAUX_io(I, J, K, 1, KK+ 1)
!                            U_F0_io(I, JJ, K) = QAUX_io(I, J, K, 3, KK+ 1) * RCCI1(JJ)
!                            IF(JJ == 1 .AND. iCase == iPIPEC) THEN
!                                U_F0_io(I, JJ, K) = (QAUX_io(I, JP, K, 2, KK+ 1) - QAUX_io(I, JP, KS, 2, KK+ 1) ) * 0.50_WP * RNDI1(JJP)
!                            ELSE
!                                U_F0_io(I, JJ, K) = (QAUX_io(I, J, K, 2, KK+ 1) * RNDI1(JJ))
!                            END IF
!                            P_F0_io(I, JJ, K) = PAUX_io(I, J, K, KK+ 1)
!                            IF(iThermoDynamics == 1) THEN
!                                T_F0_io(I, JJ, K) = TAUX_io(I, J, K, KK+ 1)
!                                D_F0_io(I, JJ, K) = DAUX_io(I, J, K, KK+ 1)
!                                H_F0_io(I, JJ, K) = HAUX_io(I, J, K, KK+ 1)
!                                K_F0_io(I, JJ, K) = KAUX_io(I, J, K, KK+ 1)
!                                M_F0_io(I, JJ, K) = MAUX_io(I, J, K, KK+ 1)
!                            END IF
!                        END DO
!                    END DO
!                END DO
!            END DO

!            DO I = NCL1S, NCL1E
!                DO K = 1, NCL3
!                    U_F0_io(I, NND2, K) = QAUX_io(I, N2DO(MYID) + 1, K, 2, NPTOT)
!                END DO
!            END DO

!            DO I = NCL1S, NCL1E
!                DO K = 1, NCL3
!                    U_F0_io(I, 0,   K) = U_F0_io(I, 1,   K) * (-1.0_WP)
!                    U_F0_io(I, NND2, K) = U_F0_io(I, NND2, K) * (-1.0_WP)
!                    U_F0_io(I, 0,   K) = U_F0_io(I, 1,   K) * (-1.0_WP)
!                    U_F0_io(I, NND2, K) = U_F0_io(I, NND2, K) * (-1.0_WP)
!                    U_F0_io(I, 0,   K) = 0.0_WP
!                END DO
!            END DO

!        END IF


!        RETURN
!    END SUBROUTINE
