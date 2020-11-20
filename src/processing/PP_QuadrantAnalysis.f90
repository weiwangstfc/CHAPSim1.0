!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> SUBROUTINE: PP_QUADRANTANALYSIS_xzL
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE PP_QUADRANTANALYSIS_xzL(iDomain_tmp)
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iDomain_tmp

    INTEGER(4) :: NCL1
    REAL(WP) :: VL1313
    REAL(WP), ALLOCATABLE :: Q(:, :, :, :)
    REAL(WP) :: U1xzL(N2DO(MYID), 4)
    REAL(WP) :: DVDL1xzL(N2DO(MYID), 3, 3)
    INTEGER(4), ALLOCATABLE :: IPV(:), IMV(:)

    INTEGER(4) :: IC, JC, KC, JJ, JM, JP, IP, IM, KP, KM, JJP
    INTEGER(4) :: N, M
    REAL(WP) :: Ucc(3), Ucc_JP, Ucc_JM, Vcc_IP, Vcc_IM, VorzJ
    REAL(WP), ALLOCATABLE ::  UVper(:, :), Vorz(:, :), TKE(:, :)
    REAL(WP), ALLOCATABLE ::  Uper(:, :), Vper(:, :), Wper(:, :), Tper(:, :), Dper(:, :)
    REAL(WP), ALLOCATABLE ::  DUVper (:, :)
    REAL(WP), ALLOCATABLE ::  DUVpper(:, :), Upper(:, :), Vpper(:, :),Wpper(:, :)

    REAL(WP) :: UVper_xzL  (N2DO(MYID)), Vorz_xzL(N2DO(MYID)), TKE_xzL(N2DO(MYID))
    REAL(WP) :: DUVper_xzL (N2DO(MYID))
    REAL(WP) :: DUVpper_xzL(N2DO(MYID))

    REAL(WP) :: QUADUVJ  (4, QUADHN), QUADVZJ(4, QUADHN), QUADTKJ(4, QUADHN), QUADDRJ(4, QUADHN)
    REAL(WP) :: QUADDUV1J(4, QUADHN)
    REAL(WP) :: QUADDUV2J(4, QUADHN)

    REAL(WP) :: OCTTUVJ  (8, QUADHN), OCTTVZJ(8, QUADHN), OCTTTKJ(8, QUADHN), OCTTDRJ(8, QUADHN)
    REAL(WP) :: OCTTDUV1J(8, QUADHN)
    REAL(WP) :: OCTTDUV2J(8, QUADHN)

    REAL(WP) :: OCTDUVJ  (8, QUADHN), OCTDVZJ(8, QUADHN), OCTDTKJ(8, QUADHN), OCTDDRJ(8, QUADHN)
    REAL(WP) :: OCTDDUV1J(8, QUADHN)
    REAL(WP) :: OCTDDUV2J(8, QUADHN)

    REAL(WP) :: UVABS, DUV1ABS, DUV2ABS

    LOGICAL :: COND_TH

    ! Attention
    ! after examinations, <\rho u'v'>  and <u'v'> shows almost the same figures....
    !

    !============== SpecIFy io or tg ==========
    IF(iDomain_tmp == ITG) THEN

        NCL1 = NCL1_TG
        VL1313 = VL1313_tg

        ALLOCATE(Q(NCL1_tg, 0 : N2DO(MYID) + 1, NCL3, NDV)); Q = 0.0_WP
        ALLOCATE(IPV(NCL1)) ; IPV = 0
        ALLOCATE(IMV(NCL1)) ; IMV = 0

        Q = Q_tg
        IPV = IPV_tg
        IMV = IMV_tg

        U1xzL = U1xzL_tg
        DVDL1xzL = DVDL1xzL_tg

    ELSE IF(iDomain_tmp == IIO) THEN

        NCL1 = NCL1_io
        VL1313 = VL1313_io

        ALLOCATE(Q(NCL1_io, 0 : N2DO(MYID) + 1, NCL3, NDV)) ; Q = 0.0_WP
        ALLOCATE(IPV(NCL1)) ; IPV = 0
        ALLOCATE(IMV(NCL1)) ; IMV = 0

        Q = Q_io
        IPV = IPV_io
        IMV = IMV_io

        U1xzL = U1xzL_io
        DVDL1xzL = DVDL1xzL_io
    ELSE
        NCL1 = NCL1_io
        VL1313 = VL1313_io
    END IF


    ALLOCATE ( Uper (NCL1, NCL3) ); Upper = 0.0_WP
    ALLOCATE ( Vper (NCL1, NCL3) ); Vpper = 0.0_WP
    ALLOCATE ( Wper (NCL1, NCL3) ); Wpper = 0.0_WP

    ALLOCATE ( UVper  (NCL1, NCL3) ) ; UVper = 0.0_WP
    ALLOCATE ( Vorz (NCL1, NCL3) )   ; Vorz = 0.0_WP
    ALLOCATE ( TKE(NCL1, NCL3) )     ; TKE = 0.0_WP

    IF(iDomain_tmp == IIO) THEN
        ALLOCATE ( Upper (NCL1, NCL3) ) ; Upper = 0.0_WP
        ALLOCATE ( Vpper (NCL1, NCL3) ) ; Vpper = 0.0_WP
        ALLOCATE ( Wpper (NCL1, NCL3) ) ; Wpper = 0.0_WP

        ALLOCATE ( Tper (NCL1, NCL3) ) ; Tper = 0.0_WP
        ALLOCATE ( Dper (NCL1, NCL3) ) ; Dper = 0.0_WP

        ALLOCATE ( DUVper (NCL1, NCL3) ) ; DUVper = 0.0_WP
        ALLOCATE ( DUVpper(NCL1, NCL3) ) ; DUVpper = 0.0_WP
    END IF

    !=========CALCULATE QUADRANT ===========================================================================
    DO  JC = 1, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJP = JGPV(JJ)

        !============ Space Averaged perturbatioN ============================================
        UVper_xzL(JC) = 0.0_WP
        Vorz_xzL(JC) = 0.0_WP
        TKE_xzL(JC)   = 0.0_WP

        IF(iDomain_tmp == IIO) THEN
            DUVper_xzL(JC) = 0.0_WP
            DUVpper_xzL(JC) = 0.0_WP
        END IF

        DO IC = 1, NCL1
            IP = IPV(IC)
            IM = IMV(IC)

            DO KC = 1, NCL3
                KP = KPV(KC)
                KM = KMV(KC)

                !======= instantanous velocity in the cell centrE =====
                Ucc(1) = ( Q(IP, JC, KC, 1) + Q(IC, JC, KC, 1) ) * 0.5_WP
                Ucc(2) = ( Q(IC, JP, KC, 2) + Q(IC, JC, KC, 2) ) * 0.5_WP
                Ucc(3) = ( Q(IC, JC, KP, 3) + Q(IC, JC, KC, 3) ) * 0.5_WP

                Ucc_JP = ( Q(IP, JP, KC, 1) + Q(IC, JP, KC, 1) ) * 0.5_WP
                Ucc_JM = ( Q(IP, JM, KC, 1) + Q(IC, JM, KC, 1) ) * 0.5_WP
                Vcc_IP = ( Q(IP, JP, KC, 2) + Q(IP, JC, KC, 2) ) * 0.5_WP
                Vcc_IM = ( Q(IM, JP, KC, 2) + Q(IM, JC, KC, 2) ) * 0.5_WP

                !========perturbationS ==============================
                Uper(IC, KC) = Ucc(1) - U1xzL(JC, 1)      ! u'
                Vper(IC, KC) = Ucc(2) - U1xzL(JC, 2)      ! v'
                Wper(IC, KC) = Ucc(3) - U1xzL(JC, 3)      ! w'

                !========= M1: {u'v'}================================
                UVper(IC, KC) = Uper(IC, KC) * Vper(IC, KC)   ! u'v'
                UVper_xzL(JC) = UVper_xzL(JC) + UVper(IC, KC)! \Simga^{xz}{u'v'}

                !======== VortICitY = Z = dv'/ DX - du'/ DY ==============
                VorzJ = 0.5_WP * ( ( Ucc(2) + Vcc_IP ) - (Ucc(2) + Vcc_IM) ) * DXI - &
                0.5_WP * ( ( Ucc_JP - Ucc(1) ) * DYCI(JJP) + ( Ucc(1) - Ucc_JM ) * DYCI(JJ) ) -&
                ( DVDL1xzL(JC, 2, 1) - DVDL1xzL(JC, 1, 2) )
                Vorz(IC, KC) = DSQRT(VorzJ * VorzJ)
                Vorz_xzL(JC) = Vorz_xzL(JC) + Vorz(IC, KC)

                !======= TKE =========================================
                TKE(IC, KC) = 0.5_WP * (Uper(IC, KC) * Uper(IC, KC) + Vper(IC, KC) * Vper(IC, KC) + Wper(IC, KC) * Wper(IC, KC))
                TKE_xzL(JC) = TKE_xzL  (JC) + TKE(IC, KC)

                IF(iDomain_tmp == IIO) THEN
                    IF(iThermoDynamics == 1) THEN
                        Tper(IC, KC) = TEMPERATURE(IC, JC, KC) - T1xzL_io(JC)        ! T'
                        Dper(IC, KC) = DENSITY    (IC, JC, KC) - D1xzL_io(JC)        ! D'
                    END IF

                    !========perturbationS ==============================
                    Upper(IC, KC) = Ucc(1) - G1xzL_io(JC, 1) / D1xzL_io(JC)      ! u"
                    Vpper(IC, KC) = Ucc(2) - G1xzL_io(JC, 2) / D1xzL_io(JC)      ! v"
                    Wpper(IC, KC) = Ucc(3) - G1xzL_io(JC, 3) / D1xzL_io(JC)      ! w"

                    !========= M2: {\rho u'v'}===========================
                    DUVper(IC, KC) = D1xzL_io(JC) * Uper(IC, KC) * Vper(IC, KC)   ! \rho * u'v'
                    DUVper_xzL(JC) = DUVper_xzL(JC) + DUVper(IC, KC)           ! \Simga^{xz}{\rho * u'v'}

                    !========= M3: {\rho u" v"}==========================
                    DUVpper(IC, KC) = D1xzL_io(JC) * Upper(IC, KC) * Vpper(IC, KC) ! \rho * u" * v"
                    DUVpper_xzL(JC) = DUVpper_xzL(JC) + DUVpper(IC, KC)         ! \Simga^{xz}{\rho * u" * v"}
                END IF

            END DO
        END DO
        UVper_xzL  (JC) = UVper_xzL(JC) * VL1313   ! <u'v'>_{xz}
        Vorz_xzL(JC) = Vorz_xzL(JC) * VL1313
        TKE_xzL(JC)   = TKE_xzL(JC)   * VL1313
        IF(iDomain_tmp == IIO) THEN
            DUVper_xzL (JC) = DUVper_xzL(JC) * VL1313   ! <\rho u'v'>_{xz}
            DUVpper_xzL(JC) = DUVpper_xzL(JC) * VL1313   ! <\rho u"v">_{xz}
        END IF

        !=============== Space sum of QUADRANT ================================================
        QUADDRJ(:, :) = 0.0_WP
        QUADUVJ(:, :) = 0.0_WP
        QUADVzJ(:, :) = 0.0_WP
        QUADTKJ(:, :) = 0.0_WP

        OCTTDRJ(:, :) = 0.0_WP
        OCTTUVJ(:, :) = 0.0_WP
        OCTTVzJ(:, :) = 0.0_WP
        OCTTTKJ(:, :) = 0.0_WP

        OCTDDRJ(:, :) = 0.0_WP
        OCTDUVJ(:, :) = 0.0_WP
        OCTDVzJ(:, :) = 0.0_WP
        OCTDTKJ(:, :) = 0.0_WP

        IF(iDomain_tmp == IIO) THEN
            QUADDUV1J(:, :) = 0.0_WP
            QUADDUV2J(:, :) = 0.0_WP

            OCTTDUV1J(:, :) = 0.0_WP
            OCTTDUV2J(:, :) = 0.0_WP

            OCTDDUV1J(:, :) = 0.0_WP
            OCTDDUV2J(:, :) = 0.0_WP
        END IF

        DO IC = 1, NCL1
            IP = IPV(IC)
            IM = IMV(IC)

            DO KC = 1, NCL3
                KP = KPV(KC)
                KM = KMV(KC)

                !======= QUADRANT FOR QUADHN TH levelS =================
                UVABS  = DABS(UVper  (IC, KC) / UVper_xzL  (JC))  ! |u'v'|/|<u'v'>_{xz}|

                IF(iDomain_tmp == IIO) THEN
                    DUV2ABS = DABS(DUVpper(IC, KC) / DUVpper_xzL(JC))  ! |\rho u"v"|/|<\rho u"v">_{xz}|
                    DUV1ABS = DABS(DUVper (IC, KC) / DUVper_xzL (JC))  ! |\rho u'v'|/|<\rho u'v'>_{xz}|
                END IF

                !======= movement direction based on u' not u" ===================
                DO M = 1, 4   ! M = quadrant phase 1, 2, 3, 4
                    COND_TH = .FALSE.

                    !==below IS based on assumption of zero-sheAR-stress at Channel centre...
                    !                        IF (M == 1) THEN
                    !                            IF(YCC(JJ) < 0.0_WP) THEN
                    !                                COND_TH = ( (Uper(IC, KC) >   0.0_WP) .AND. (Vper(IC, KC) >  0.0_WP) )
                    !                            ELSE
                    !                                COND_TH = ( (Uper(IC, KC) >   0.0_WP) .AND. (Vper(IC, KC) < 0.0_WP) )
                    !                            END IF
                    !                        ELSE IF(M == 2) THEN
                    !                            IF(YCC(JJ) < 0.0_WP) THEN
                    !                                COND_TH = ( (Uper(IC, KC) <= 0.0_WP) .AND. (Vper(IC, KC) >= 0.0_WP) )
                    !                            ELSE
                    !                                COND_TH = ( (Uper(IC, KC) <= 0.0_WP) .AND. (Vper(IC, KC) <= 0.0_WP) )
                    !                            END IF
                    !                        ELSE IF(M == 3) THEN
                    !                            IF(YCC(JJ) < 0.0_WP) THEN
                    !                                COND_TH = ( (Uper(IC, KC)  <  0.0_WP) .AND. (Vper(IC, KC) < 0.0_WP) )
                    !                            ELSE
                    !                                COND_TH = ( (Uper(IC, KC)  <  0.0_WP) .AND. (Vper(IC, KC) >  0.0_WP) )
                    !                            END IF
                    !                        ELSE IF(M == 4) THEN
                    !                            IF(YCC(JJ) < 0.0_WP) THEN
                    !                                COND_TH = ( (Uper(IC, KC) >= 0.0_WP) .AND. (Vper(IC, KC) <= 0.0_WP) )
                    !                            ELSE
                    !                                COND_TH = ( (Uper(IC, KC) >= 0.0_WP) .AND. (Vper(IC, KC) >= 0.0_WP) )
                    !                            END IF
                    !                        ELSE
                    !                            COND_TH = .FALSE.
                    !                        END IF

                    !===== Two wall sides using the same signs of uv
                    IF (M == 1) THEN
                        COND_TH = ( (Uper(IC, KC) >   0.0_WP) .AND. (Vper(IC, KC) >  0.0_WP) )
                    ELSE IF(M == 2) THEN
                        COND_TH = ( (Uper(IC, KC) <= 0.0_WP) .AND. (Vper(IC, KC) >= 0.0_WP) )
                    ELSE IF(M == 3) THEN
                        COND_TH = ( (Uper(IC, KC)  <  0.0_WP) .AND. (Vper(IC, KC) < 0.0_WP) )
                    ELSE IF(M == 4) THEN
                        COND_TH = ( (Uper(IC, KC) >= 0.0_WP) .AND. (Vper(IC, KC) <= 0.0_WP) )
                    ELSE
                        COND_TH = .FALSE.
                    END IF

                    DO N = 1, QUADHN    ! N = number of TH levels
                        !=======|u'v'|===============
                        IF(COND_TH .AND. (UVABS>=QUADHV(N)) ) THEN

                            QUADUVJ(M, N) = QUADUVJ(M, N) + UVper(IC, KC)

                            !====below three has the same condition as u'v'==
                            QUADVzJ(M, N) = QUADVzJ(M, N) + Vorz (IC, KC)
                            QUADTKJ(M, N) = QUADTKJ(M, N) + TKE(IC, KC)
                            QUADDRJ(M, N) = QUADDRJ(M, N) + 1.0_WP

                            !====for octant analysIS ==========================================
                            IF( (iDomain_tmp == IIO) .AND. (iThermoDynamics == 1) ) THEN
                                IF(Tper(IC, KC) >= 0.0_WP) THEN
                                    OCTTUVJ(M, N) = OCTTUVJ(M, N) + UVper(IC, KC)

                                    !====below three has the same condition as u'v'==
                                    OCTTVzJ(M, N) = OCTTVzJ(M, N) + Vorz (IC, KC)
                                    OCTTTKJ(M, N) = OCTTTKJ(M, N) + TKE(IC, KC)
                                    OCTTDRJ(M, N) = OCTTDRJ(M, N) + 1.0_WP
                                ELSE
                                    OCTTUVJ(M + 4, N) = OCTTUVJ(M + 4, N) + UVper(IC, KC)

                                    !====below three has the same condition as u'v'==
                                    OCTTVzJ(M + 4, N) = OCTTVzJ(M + 4, N) + Vorz (IC, KC)
                                    OCTTTKJ(M + 4, N) = OCTTTKJ(M + 4, N) + TKE(IC, KC)
                                    OCTTDRJ(M + 4, N) = OCTTDRJ(M + 4, N) + 1.0_WP
                                END IF

                                IF(Dper(IC, KC) >= 0.0_WP) THEN
                                    OCTDUVJ(M, N) = OCTDUVJ(M, N) + UVper(IC, KC)

                                    !====below three has the same condition as u'v'==
                                    OCTDVzJ(M, N) = OCTDVzJ(M, N) + Vorz (IC, KC)
                                    OCTDTKJ(M, N) = OCTDTKJ(M, N) + TKE(IC, KC)
                                    OCTDDRJ(M, N) = OCTDDRJ(M, N) + 1.0_WP
                                ELSE
                                    OCTDUVJ(M + 4, N) = OCTDUVJ(M + 4, N) + UVper(IC, KC)

                                    !====below three has the same condition as u'v'==
                                    OCTDVzJ(M + 4, N) = OCTDVzJ(M + 4, N) + Vorz (IC, KC)
                                    OCTDTKJ(M + 4, N) = OCTDTKJ(M + 4, N) + TKE(IC, KC)
                                    OCTDDRJ(M + 4, N) = OCTDDRJ(M + 4, N) + 1.0_WP
                                END IF
                            END IF

                        END IF

                        IF((iDomain_tmp == IIO) .AND. (iThermoDynamics == 1)) THEN
                            !=======|\rho u'v'|===============
                            IF(COND_TH .AND. (DUV1ABS>=QUADHV(N)) ) THEN
                                QUADDUV1J(M, N) = QUADDUV1J(M, N) + DUVper(IC, KC)
                                !====for octant analysIS ==========================================
                                IF(Tper(IC, KC) >= 0.0_WP) THEN
                                    OCTTDUV1J(M, N) = OCTTDUV1J(M, N) + DUVper(IC, KC)
                                ELSE
                                    OCTTDUV1J(M + 4, N) = OCTTDUV1J(M + 4, N) + DUVper(IC, KC)
                                END IF

                                IF(Dper(IC, KC) >= 0.0_WP) THEN
                                    OCTDDUV1J(M, N) = OCTDDUV1J(M, N) + DUVper(IC, KC)
                                ELSE
                                    OCTDDUV1J(M + 4, N) = OCTDDUV1J(M + 4, N) + DUVper(IC, KC)
                                END IF
                            END IF
                        END IF
                    END DO

                END DO

                IF(iDomain_tmp == IIO) THEN
                    !======= movement direction based on  u" ===================
                    DO M = 1, 4   ! M = quadrant phase 1, 2, 3, 4
                        COND_TH = .FALSE.
                        !==below IS based on assumption of zero-sheAR-stress at Channel centre...

                        !                            IF (M == 1) THEN
                        !                                IF(YCC(JJ) < 0.0_WP) THEN
                        !                                    COND_TH = ( (Upper(IC, KC) >   0.0_WP) .AND. (Vpper(IC, KC) >  0.0_WP) )
                        !                                ELSE
                        !                                    COND_TH = ( (Upper(IC, KC) >   0.0_WP) .AND. (Vpper(IC, KC) < 0.0_WP) )
                        !                                END IF
                        !                            ELSE IF(M == 2) THEN
                        !                                IF(YCC(JJ) < 0.0_WP) THEN
                        !                                    COND_TH = ( (Upper(IC, KC) <= 0.0_WP) .AND. (Vpper(IC, KC) >= 0.0_WP) )
                        !                                ELSE
                        !                                    COND_TH = ( (Upper(IC, KC) <= 0.0_WP) .AND. (Vpper(IC, KC) <= 0.0_WP) )
                        !                                END IF
                        !                            ELSE IF(M == 3) THEN
                        !                                IF(YCC(JJ) < 0.0_WP) THEN
                        !                                    COND_TH = ( (Upper(IC, KC)  <  0.0_WP) .AND. (Vpper(IC, KC) < 0.0_WP) )
                        !                                ELSE
                        !                                    COND_TH = ( (Upper(IC, KC)  <  0.0_WP) .AND. (Vpper(IC, KC) >  0.0_WP) )
                        !                                END IF
                        !                            ELSE IF(M == 4) THEN
                        !                                IF(YCC(JJ) < 0.0_WP) THEN
                        !                                    COND_TH = ( (Upper(IC, KC) >= 0.0_WP) .AND. (Vpper(IC, KC) <= 0.0_WP) )
                        !                                ELSE
                        !                                    COND_TH = ( (Upper(IC, KC) >= 0.0_WP) .AND. (Vpper(IC, KC) >= 0.0_WP) )
                        !                                END IF
                        !                            ELSE
                        !                                COND_TH = .FALSE.
                        !                            END IF

                        !===== Two wall sides using the same signs of uv
                        IF (M == 1) THEN
                            COND_TH = ( (Upper(IC, KC) >   0.0_WP) .AND. (Vpper(IC, KC) >  0.0_WP) )
                        ELSE IF(M == 2) THEN
                            COND_TH = ( (Upper(IC, KC) <= 0.0_WP) .AND. (Vpper(IC, KC) >= 0.0_WP) )
                        ELSE IF(M == 3) THEN
                            COND_TH = ( (Upper(IC, KC)  <  0.0_WP) .AND. (Vpper(IC, KC) < 0.0_WP) )
                        ELSE IF(M == 4) THEN
                            COND_TH = ( (Upper(IC, KC) >= 0.0_WP) .AND. (Vpper(IC, KC) <= 0.0_WP) )
                        ELSE
                            COND_TH = .FALSE.
                        END IF

                        DO N = 1, QUADHN    ! N = number of TH levels
                            !=======|\rho u"v"|===============
                            IF(COND_TH .AND. (DUV2ABS>=QUADHV(N)) ) THEN
                                QUADDUV2J(M, N) = QUADDUV2J(M, N) + DUVpper(IC, KC)
                                !====for octant analysIS ==========================================
                                IF(Tper(IC, KC) >= 0.0_WP) THEN
                                    OCTTDUV2J(M, N) = OCTTDUV2J(M, N) + DUVpper(IC, KC)
                                ELSE
                                    OCTTDUV2J(M + 4, N) = OCTTDUV2J(M + 4, N) + DUVpper(IC, KC)
                                END IF

                                IF(Dper(IC, KC) >= 0.0_WP) THEN
                                    OCTDDUV2J(M, N) = OCTDDUV2J(M, N) + DUVpper(IC, KC)
                                ELSE
                                    OCTDDUV2J(M + 4, N) = OCTDDUV2J(M + 4, N) + DUVpper(IC, KC)
                                END IF
                            END IF
                        END DO

                    END DO
                END IF


            END DO
        END DO

        !            IF(MYID == 0) THEN
        !                WRITE(*, '(1I2.1,6I8.1)')   JC, QUADNSJ(1:4, 1), QUADNSJ(1, 1) + QUADNSJ(2, 1) + QUADNSJ(3, 1) + QUADNSJ(4, 1), NCL1 * NCL3
        !                WRITE(*, '(1I2.1,6ES13.5)') JC, QUADUVJ(1:4, 1) / DBLE(QUADNSJ(1:4, 1)), &
        !                               QUADUVJ(1, 1) / DBLE(NCL1 * NCL3) + QUADUVJ(2, 1) / DBLE(NCL1 * NCL3) + QUADUVJ(3, 1) / DBLE(NCL1 * NCL3) +&
        !                               QUADUVJ(4, 1) / DBLE(NCL1 * NCL3),UVper_xzL(JC)
        !            END IF

        !========= Space averaged  QUADRANT ==================
        DO N = 1, QUADHN  ! TH levels
            !WRITE(*, *) N, M, QUADNSJ(M, N), QUADUVJ(M, N), UVper_xzL(JC),UVper_xzL(JC), TKE_xzL(JC), QUADDRJ(M, N)
            IF(iDomain_tmp == ITG) THEN

                QUADUVxzL_TG(JC, :, N) = QUADUVJ(:, N) * VL1313 / (UVper_xzL(JC))

                QUADVzxzL_TG(JC, :, N) = QUADVzJ(:, N) * VL1313 / (Vorz_xzL(JC))
                QUADTKxzL_TG(JC, :, N) = QUADTKJ(:, N) * VL1313 / (TKE_xzL(JC))
                QUADDRxzL_TG(JC, :, N) = QUADDRJ(:, N) * VL1313

            ELSE IF(iDomain_tmp == IIO) THEN

                QUADUVxzL_io(JC, :, N) = QUADUVJ(:, N) * VL1313 / (UVper_xzL(JC))
                QUADVzxzL_io(JC, :, N) = QUADVzJ(:, N) * VL1313 / (Vorz_xzL(JC))
                QUADTKxzL_io(JC, :, N) = QUADUVJ(:, N) * VL1313 !for test ! not ratio, but abosulute value
                !QUADTKxzL_io(JC, :, N) = QUADTKJ(:, N) * VL1313 / (TKE_xzL(JC))


                QUADDRxzL_io(JC, :, N) = QUADDRJ(:, N) * VL1313

                IF(iThermoDynamics == 1) THEN
                    QUADDUV1xzL_io(JC, :, N) = QUADDUV1J(:, N) * VL1313 / (DUVper_xzL(JC))
                    QUADDUV2xzL_io(JC, :, N) = QUADDUV2J(:, N) * VL1313 / (DUVpper_xzL(JC))


                    OCTTUVxzL_io(JC, :, N) = OCTTUVJ(:, N) * VL1313 / (UVper_xzL(JC))
                    OCTTVzxzL_io(JC, :, N) = OCTTVzJ(:, N) * VL1313 / (Vorz_xzL(JC))
                    OCTTTKxzL_io(JC, :, N) = OCTTTKJ(:, N) * VL1313 / (TKE_xzL(JC))
                    OCTTDRxzL_io(JC, :, N) = OCTTDRJ(:, N) * VL1313
                    OCTTDUV1xzL_io(JC, :, N) = OCTTDUV1J(:, N) * VL1313 / (DUVper_xzL(JC))
                    OCTTDUV2xzL_io(JC, :, N) = OCTTDUV2J(:, N) * VL1313 / (DUVpper_xzL(JC))

                    OCTDUVxzL_io(JC, :, N) = OCTDUVJ(:, N) * VL1313 / (UVper_xzL(JC))
                    OCTDVzxzL_io(JC, :, N) = OCTDVzJ(:, N) * VL1313 / (Vorz_xzL(JC))
                    OCTDTKxzL_io(JC, :, N) = OCTDTKJ(:, N) * VL1313 / (TKE_xzL(JC))
                    OCTDDRxzL_io(JC, :, N) = OCTDDRJ(:, N) * VL1313
                    OCTDDUV1xzL_io(JC, :, N) = OCTDDUV1J(:, N) * VL1313 / (DUVper_xzL(JC))
                    OCTDDUV2xzL_io(JC, :, N) = OCTDDUV2J(:, N) * VL1313 / (DUVpper_xzL(JC))
                END IF
            ELSE
            END IF

            !IF(MYID == 0) WRITE(*, '(2I3.1, 5ES13.5)') N, JC, QUADUVxzL_io(JC, 1:4, N), &
            !            QUADUVxzL_io(JC, 1, N) + QUADUVxzL_io(JC, 2, N) + QUADUVxzL_io(JC, 3, N) + QUADUVxzL_io(JC, 4, N)
        END DO


        !WRITE(*, '(1I2.1, 5ES13.5)') JC, QUADUVxzL_io(JC, 1, 1), QUADUVxzL_io(JC, 2, 1), QUADUVxzL_io(JC, 3, 1), QUADUVxzL_io(JC, 4, 1), &
        !              QUADUVxzL_io(JC, 1, 1) + QUADUVxzL_io(JC, 2, 1) + QUADUVxzL_io(JC, 3, 1) + QUADUVxzL_io(JC, 4, 1)

    END DO

    DEALLOCATE ( Q     )
    DEALLOCATE ( IPV   )
    DEALLOCATE ( IMV   )

    DEALLOCATE ( Uper  )
    DEALLOCATE ( Vper  )
    DEALLOCATE ( Wper  )

    DEALLOCATE ( UVper   )

    DEALLOCATE ( Vorz  )
    DEALLOCATE ( TKE   )

    IF(iDomain_tmp == IIO) THEN
        DEALLOCATE ( DUVper  )
        DEALLOCATE ( DUVpper )
        DEALLOCATE ( Upper  )
        DEALLOCATE ( Vpper  )
        DEALLOCATE ( Wpper  )
    END IF


    RETURN
END SUBROUTINE
