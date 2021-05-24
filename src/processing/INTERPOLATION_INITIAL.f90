!**********************************************************************************************************************************
!> @brief
!>        data mapping form one mesh to another.
!>        step 1: repeating/cutting the original domain to match the current domain.
!>        step 2: mapping
!> @details
!> module: INTP_VARS
!> SUBROUTINE: INITIAL_INTERP_MESH (in MYID = all)
!> SUBROUTINE: INITIAL_INTERP_tg
!> SUBROUTINE: INITIAL_INTERP_io

!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE INTP_VARS
    USE WPRECISION
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info

    ! the initial mesh info
    INTEGER(4) :: NCLOO1, NCLOO3
    REAL(WP) :: HXOO, HZOO

    ! the real used (processed old) mesh for data mapping
    INTEGER(4) :: NCLO1, NCLO2, NCLO3
    REAL(WP) :: HXO, HZO
    INTEGER(4), ALLOCATABLE :: N2DOO(:)
    INTEGER(4), ALLOCATABLE :: N3DOO(:)
    REAL(WP), ALLOCATABLE :: XNDO(:)
    REAL(WP), ALLOCATABLE :: YNDO(:)
    REAL(WP), ALLOCATABLE :: ZNDO(:)
    REAL(WP), ALLOCATABLE :: XCCO(:)
    REAL(WP), ALLOCATABLE :: YCCO(:)
    REAL(WP), ALLOCATABLE :: ZCCO(:)

    ! the current mesh info
    INTEGER(4) :: NCL1
    REAL(WP) :: HX
    !REAL(WP), ALLOCATABLE :: ZCC(:)
    REAL(WP), ALLOCATABLE :: XND(:)
    REAL(WP), ALLOCATABLE :: XCC(:)
END MODULE

!=========================================================================================================
SUBROUTINE INITIAL_INTERP_MESH
    USE INTP_VARS
    IMPLICIT NONE

    CHARACTER(128) :: SECT          ! Section names
    CHARACTER(128) :: SKEY          ! key names
    INTEGER(4) :: IOS = 0
    INTEGER(4) :: INI = 13             !file I /O id
    INTEGER(4) :: LENS
    !REAL(WP) :: RTMP
    !CHARACTER(256) :: STMP
    INTEGER(4) :: I, J, K


    !******* Read in old y mesh*******************************************************
    IF(MYID == 0) CALL CHKHDL('    Reading in intp y coordinates', MYID)

    OPEN(INI, FILE = 'INTP_COORD_YND.ini', STATUS = 'old', IOSTAT = IOS)
    IF(IOS  /= 0)  &
    CALL ERRHDL(' File INTP_COORD_YND.ini cannot be found.', MYID)

    READ(INI, *) SECT
    LENS = LEN_TRIM(SECT)
    IF (MYID == 0) CALL CHKHDL('      Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[origeo]')  &
    CALL ERRHDL('Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, HXOO
    READ(INI, *) SKEY, HZOO
    READ(INI, *) SKEY, NCLOO1
    READ(INI, *) SKEY, NCLO2
    READ(INI, *) SKEY, NCLOO3

    IF(MYID == 0) THEN
        IF((NCLO2 / NPTOT) < 1) CALL ERRHDL(' Interpolation from coarse mesh with a CPU number lARger than NCL2', MYID)
    END IF

    !******* Allocate variables*************************************************

    ALLOCATE (YNDO(NCLO2 + 1));  YNDO = 0.0_WP
    READ(INI, *) SECT
    LENS = LEN_TRIM(SECT)
    IF (MYID == 0) CALL CHKHDL('      Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[ynd]')  &
    CALL ERRHDL('Reading fails: ' // SECT(1:LENS), MYID)


    DO J = 1, NCLO2 + 1
        READ(INI, *) LENS, YNDO(J)
    END DO
    CLOSE(INI)

    !******* Repeating (Not stretching) to the real DOmAIN******************
    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        WRITE(logflg_pg, *) 'MODIFY CODE, PLEASE'
    ELSE
        IF(TgFlowFlg) THEN
            HX  = HX_TG
            NCL1  = NCL1_TG
        END IF
        IF(IoFlowFlg) THEN
            HX  = HX_io
            NCL1  = NCL1_io
        END IF
        ALLOCATE (XND(NCL1 + 1)) ;  XND = 0.0_WP
        ALLOCATE (XCC(NCL1)) ;  XCC = 0.0_WP
        IF(TgFlowFlg) THEN
            XND = XND_tg
            XCC = XCC_tg
        END IF
        IF(IoFlowFlg) THEN
            XND = XND_io
            XCC = XCC_io
        END IF
    END IF

    NCLO1 = CEILING(HX / (HXOO / DBLE(NCLOO1)))
    NCLO3 = CEILING(HZ/ (HZOO / DBLE(NCLOO3)))
    HXO = HX
    HZO = HZ

    IF(MYID == 0) THEN
        CALL CHKHDL    ('      Repeating/cutting the originl domain to the real domain', MYID)
        CALL CHK2RLHDL ('         Lx00->Lx0', MYID, HXOO, HXO)
        CALL CHK2RLHDL ('         Lz00->Lz0', MYID, HZOO, HZO)
        CALL CHK2INTHDL('         Nx00->Nx0', MYID, NCLOO1, NCLO1)
        CALL CHK2INTHDL('         Nz00->Nz0', MYID, NCLOO3, NCLO3)

        CALL CHKHDL    ('      Mapping domain from the given to the current', MYID)
        CALL CHK2RLHDL ('         Lx0->Lx', MYID, HXO, HX)
        CALL CHK2RLHDL ('         Lz0->Lz', MYID, HZO, HZ)
        CALL CHK2INTHDL('         Nx0->Nx', MYID, NCLO1, NCL1)
        CALL CHK2INTHDL('         Nz0->Nz', MYID, NCLO3, NCL3)
    END IF



    ALLOCATE (N3DOO(0 : NPSLV)) ;  N3DOO = 0
    ALLOCATE (N2DOO(0 : NPSLV)) ;  N2DOO = 0
    N2DOO(0 : NPSLV) = INT(NCLO2 / NPTOT)
    N3DOO(0 : NPSLV) = INT(NCLO3 / NPTOT)
    !********* OLD AND NEW COORDINATES******************************************************************
    ALLOCATE (YCCO(NCLO2) ) ;  YCCO = 0.0_WP
    ALLOCATE (XNDO(NCLO1 + 1));  XNDO = 0.0_WP
    ALLOCATE (ZNDO(NCLO3 + 1));  ZNDO = 0.0_WP
    ALLOCATE (XCCO(NCLO1));  XCCO = 0.0_WP
    ALLOCATE (ZCCO(NCLO3));  ZCCO = 0.0_WP
    !ALLOCATE (ZCC(NCL3))

    DO I = 1, (NCLO1 + 1)
        XNDO(I) = DBLE(I - 1) * HXO / DBLE(NCLO1)
    END DO

    DO  k = 1, (NCLO3 + 1)
        ZNDO(K) = DBLE(K - 1) * HZO / DBLE(NCLO3)
    END DO

    DO I = 1, NCLO1
        XCCO(I) = 0.5_WP * (XNDO(I) + XNDO(I + 1))
    END DO

    DO K = 1, NCLO3
        ZCCO(K) = 0.5_WP * (ZNDO(K) + ZNDO(K+ 1))
    END DO

    DO J = 1, NCLO2
        YCCO(J) = 0.5_WP * ( YNDO(J) + YNDO(J + 1) )
    END DO

    !DO K = 1, NCL3
    !ZCC(K) = 0.5_WP * ( ZND(K) + ZND(K+ 1) )
    !END DO
    CALL MPI_BARRIER(ICOMM, IERROR)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INITIAL_INTERP_tg
    USE INTP_VARS
    IMPLICIT NONE

    CHARACTER(LEN = 256) :: WRT_RST_FNM
    CHARACTER(15) :: PNTIM
    INTEGER(4) :: I, J, K, N, KK, II
    REAL(WP) :: Umean1_tg
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
    REAL(WP), ALLOCATABLE :: QO   (:, :, :)
    REAL(WP), ALLOCATABLE :: QN   (:, :, :)

    if(iIniField_tg /= 1) return

    IF(MYID == 0) CALL CHKHDL(' TG: Interpolating from one mesh to another', MYID)

    !========= STEP 1, GENERATE MESH INFO==================
    CALL INITIAL_INTERP_MESH

    !========= STEP 2 ALLOCATE VARS ===================
    ALLOCATE( DUMMY (NCLOO1, N2DOO(MYID), NCLOO3) ); DUMMY = 0.0_WP
    ALLOCATE( QO    (NCLO1, N2DOO(MYID), NCLO3)  ); QO = 0.0_WP
    ALLOCATE( QN    (NCL1, N2DO(MYID), NCL3 )  ); QN = 0.0_WP

    !========== STEP 3 READ IN OLD DATA AND INTP INTO NEW ONES ===============
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg
    DO N = 1, NDV + 1
        !========== (1) FILE NAME ===========
        IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_U.D'
        IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_V.D'
        IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_W.D'
        IF(N == 4)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_P.D'

        !========== (2) READ IN DATA===========
        CALL READ_3D_VARS(NCLOO1, NCLO2, NCLOO3, N2DOO(MYID), N2DOO(0) * MYID, ITERG0_TG,PhyTIME_TG, DUMMY, WRT_RST_FNM)

        !========== (3) EXPAND READ -IN TO REAL LENGTH ===========
        DO I = 1, NCLO1
            DO J = 1, N2DOO(MYID)
                DO K = 1, NCLO3

                    II = MOD(I, NCLOO1)
                    IF(II == 0) II = NCLOO1

                    KK = MOD(K, NCLOO3)
                    IF(KK == 0) KK = NCLOO3

                    QO(I, J, K) = DUMMY(II, J, KK)
                END DO
            END DO
        END DO

        !========== (4) INTEP DATA===========
        CALL INTP_3D_VARS(QO, QN, N)

        DO J = 1, N2DO(MYID)
            DO I = 1, NCL1
                DO K = 1, NCL3
                    IF(N == 1) Q_tg(I, J, K, 1) = QN(I, J, K) !U
                    IF(N == 2) Q_tg(I, J, K, 2) = QN(I, J, K) !V
                    IF(N == 3) Q_tg(I, J, K, 3) = QN(I, J, K) !W
                    IF(N == 4) PR_TG(I, J, K) = QN(I, J, K) !P
                END DO
            END DO
        END DO

    END DO

    PhyTIME_tg = 0.0_WP
    ITERG0_tg= 0
    PhyTIME = PhyTIME_tg
    ITERG0 = ITERG0_tg



    DEALLOCATE (DUMMY)
    DEALLOCATE (QO)
    DEALLOCATE (QN)
    DEALLOCATE  (N2DOO)
    DEALLOCATE  (N3DOO)
    DEALLOCATE  (XNDO)
    DEALLOCATE  (YNDO)
    DEALLOCATE  (ZNDO)
    DEALLOCATE  (XCCO)
    DEALLOCATE  (YCCO)
    DEALLOCATE  (ZCCO)
    !DEALLOCATE  (ZCC)
    DEALLOCATE  (XCC)
    DEALLOCATE  (XND)

    CALL CALL_TEC360

    CALL BULK_VELOCITY_TG(Umean1_tg)
    IF(MYID == 0) CALL CHKRLHDL  ('         TG: The bulk velocity (orignal) = ', MYID, Umean1_tg)

    DO J = 1, N2DO(MYID)
        Q_tg(:, J, :, NFLOW) = Q_tg(:, J, :, NFLOW) / Umean1_tg
    END DO

    CALL BULK_VELOCITY_TG(Umean1_tg)
    IF(MYID == 0) CALL CHKRLHDL  ('         TG: The bulk velocity (corcted) = ', MYID, Umean1_tg)


    CALL MPI_BARRIER(ICOMM, IERROR)


    IF(MYID == 0) CALL CHKRLHDL('    TG: END of interplating from a coarse mesh at time', MYID, PhyTIME_TG)



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INITIAL_INTERP_io
    USE INTP_VARS
    IMPLICIT NONE

    CHARACTER(LEN = 256) :: WRT_RST_FNM
    CHARACTER(15) :: PNTIM
    INTEGER(4) :: I, J, K, N, KK, II, JJ
    REAL(WP) :: Umean1_io, Gmean1_io, HMEAN, HMEAN_WORK, HMEANO_WORK, VOLOO, VOLOO_WORK
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
    REAL(WP), ALLOCATABLE :: QO   (:, :, :)
    REAL(WP), ALLOCATABLE :: QN   (:, :, :)

    if(iIniField_io /= 1) return

    IF(MYID == 0) CALL CHKHDL(' IO: Interpolating from one mesh to another', MYID)

    !========= STEP 1, GENERATE MESH INFO==================
    CALL INITIAL_INTERP_MESH

    !========= STEP 2 ALLOCATE VARS ===================
    ALLOCATE( DUMMY (NCLOO1, N2DOO(MYID), NCLOO3) ); DUMMY = 0.0_WP
    ALLOCATE( QO    (NCLO1, N2DOO(MYID), NCLO3)  ); QO = 0.0_WP
    ALLOCATE( QN    (NCL1, N2DO(MYID), NCL3 )  ); QN = 0.0_WP

    !========== STEP 3 READ IN OLD DATA AND INTP INTO NEW ONES ===============
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    DO N = 1, NDV + 1
        !========== (1) FILE NAME ===========
        IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_U.D'
        IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_V.D'
        IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_W.D'
        IF(N == 4)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_P.D'

        !========== (2) READ IN DATA===========
        CALL READ_3D_VARS(NCLOO1, NCLO2, NCLOO3, N2DOO(MYID), N2DOO(0) * MYID, ITERG0_io, PhyTIME_IO, DUMMY, WRT_RST_FNM)

        !========== (3) EXPAND READ -IN TO REAL LENGTH ===========
        DO I = 1, NCLO1
            DO J = 1, N2DOO(MYID)
                DO K = 1, NCLO3

                    II = MOD(I, NCLOO1)
                    IF(II == 0) II = NCLOO1

                    KK = MOD(K, NCLOO3)
                    IF(KK == 0) KK = NCLOO3

                    QO(I, J, K) = DUMMY(II, J, KK)
                END DO
            END DO
        END DO

        !========== (4) INTEP DATA===========
        CALL INTP_3D_VARS(QO, QN, N)

        IF(iIniFieldType == 0) THEN
            DO I = 1, NCL1E
                DO J = 1, N2DO(MYID)
                    DO K = 1, NCL3
                        IF(N == 1) G_io(I, J, K, 1) = QN(I, J, K) !RHO U
                        IF(N == 2) G_io(I, J, K, 2) = QN(I, J, K) !RHO V
                        IF(N == 3) G_io(I, J, K, 3) = QN(I, J, K) !RHO W
                        IF(N == 4) PR_io(I, J, K) = QN(I, J, K) !P
                    END DO
                END DO
            END DO

            IF(N == NFLOW) THEN
                CALL BULK_MASSFLUX_io( Gmean1_io )
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk G (original) =      ', MYID, Gmean1_io)


                DO J = 1, N2DO(MYID)
                    G_io(:, J, :, NFLOW) = G_io(:, J, :, NFLOW) /Gmean1_io
                END DO

                CALL BULK_MASSFLUX_io( Gmean1_io )
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk G (corrected) =     ', MYID, Gmean1_io)

            END IF

        END IF

        IF(iIniFieldType == 1) THEN

            DO I = 1, NCL1E
                DO J = 1, N2DO(MYID)
                    DO K = 1, NCL3
                        IF(N == 1) Q_io(I, J, K, 1) = QN(I, J, K) !U
                        IF(N == 2) Q_io(I, J, K, 2) = QN(I, J, K) !V
                        IF(N == 3) Q_io(I, J, K, 3) = QN(I, J, K) !W
                        IF(N == 4) PR_io(I, J, K) = QN(I, J, K) !P
                    END DO
                END DO
            END DO

            IF(N == NFLOW) THEN
                CALL BULK_VELOCITY_io( Umean1_io )
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk velocity(original) = ', MYID, Umean1_io)


                DO J = 1, N2DO(MYID)
                    Q_io(:, J, :, NFLOW) = Q_io(:, J, :, NFLOW) / Umean1_io
                END DO

                CALL BULK_VELOCITY_io( Umean1_io )
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk velocity (coreted) = ', MYID, Umean1_io)

            END IF

        END IF

    END DO

    IF(iThermoDynamics == 1 .AND. iIniFieldType == 0 ) THEN
        DO N = 1, 3
            !========== (1) FILE NAME ===========
            IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_T.D'
            IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_D.D'
            IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath3) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_H.D'

            !========== (2) READ IN DATA===========
            CALL READ_3D_VARS(NCLOO1, NCLO2, NCLOO3, N2DOO(MYID), N2DOO(0) * MYID, ITERG0_IO, PhyTIME_IO, DUMMY, WRT_RST_FNM)

            !========== (3) EXPAND READ -IN TO REAL LENGTH ===========
            DO I = 1, NCLO1
                DO J = 1, N2DOO(MYID)
                    DO K = 1, NCLO3

                        II = MOD(I, NCLOO1)
                        IF(II == 0) II = NCLOO1

                        KK = MOD(K, NCLOO3)
                        IF(KK == 0) KK = NCLOO3

                        QO(I, J, K) = DUMMY(II, J, KK)
                    END DO
                END DO
            END DO

            !======================================
            IF(N == 3) THEN
                DO J = 1, N2DOO(MYID)
                    JJ = MYID * N2DOO(0) + J
                    DO I = 1, NCLO1
                        DO K = 1, NCLO3
                            HMEAN = HMEAN + QO(I, J, K) * (YNDO(JJ + 1) - YNDO(JJ)) * (HXO / DBLE(NCLO1)) * (HZO / DBLE(NCLO3))
                            VOLOO = VOLOO + (YNDO(JJ + 1) - YNDO(JJ)) * (HXO / DBLE(NCLO1)) * (HZO / DBLE(NCLO3))
                        END DO
                    END DO
                END DO
                CALL MPI_BARRIER(ICOMM, IERROR)
                CALL MPI_ALLREDUCE(HMEAN, HMEANO_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
                CALL MPI_ALLREDUCE(VOLOO, VOLOO_WORK, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)

                HMEANO_WORK = HMEANO_WORK / VOLOO_WORK

                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk H orig-coarse mesh = ', MYID, HMEANO_WORK)

            END IF

            !========== (4) INTEP DATA===========
            CALL INTP_3D_VARS(QO, QN, 0)

            DO J = 1, N2DO(MYID)
                DO I = 1, NCL1
                    DO K = 1, NCL3
                        IF(N == 1) TEMPERATURE(I, J, K) = QN(I, J, K) !T
                        IF(N == 2) DENSITY    (I, J, K) = QN(I, J, K) !D
                        IF(N == 3) ENTHALPY   (I, J, K) = QN(I, J, K) !E
                    END DO
                END DO
            END DO
            !================================
            IF(N == 3) THEN
                CALL BULK_H_io(HMEAN_WORK)
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk H orig-fine mesh = ', MYID, HMEAN_WORK)

            END IF

        END DO

        !================= Scaled to originaL ==================
        DO J = 1, N2DO(MYID)
            ENTHALPY(:, J, :) = ENTHALPY(:, J, :) * HMEANO_WORK /HMEAN_WORK
        END DO

        IF(N == 3) THEN
            CALL BULK_H_io(HMEAN_WORK)
            IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk H-scaled fine mesH =       ', MYID, HMEAN_WORK)
        END IF


        !=============BUILD UP DH ===========================================
        DO J = 1, N2DO(MYID)
            DO I = 1, NCL1E
                DO K = 1, NCL3
                    DH(I, J, K) = DENSITY(I, J, K) * ENTHALPY(I, J, K)
                END DO
            END DO
        END DO

        !===========build up thermal BC and other thermal variables =======
        CALL BC_WALL_THERMAL(IALLDM)
        IF(TgFlowFlg) CALL BC_TINLET_THERML

        CALL THERM_PROP_UPDATE(IALLDM)
        IF(TgFlowFlg) THEN
            CALL INTFC_OUL_THERMAL_io
            CALL INTFC_INL_THERMAL_io
        END IF
        CALL INTFC_MFD_THERMAL_io


        !=============update velocity based on the given G and DENSITY ===================
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
        CALL BC_WALL_G_io
        CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
        CALL BC_WALL_Pr_io
        IF(TgFlowFlg) THEN
            CALL BC_TINLET_FLOW
        END IF
        CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
        CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
        CALL BC_WALL_G_io
        CALL BC_WALL_PR_io

        CALL DENSITY_Staggered
        CALL MU_Staggered
        CALL VELOCITY_CALC_io


    ELSE
        !=======build up thermal field=================
        !CALL IniField_THERMAL_io
        IF(iThermoDynamics == 1) THEN
            CALL BC_WALL_THERMAL(IALLDM)
            IF(TgFlowFlg) CALL BC_TINLET_THERML
            CALL THERM_PROP_UPDATE(IALLDM)
            IF(TgFlowFlg) THEN
                CALL INTFC_OUL_THERMAL_io
                CALL INTFC_INL_THERMAL_io
            END IF
            CALL INTFC_MFD_THERMAL_io
        END IF

        !==========flow field===================================
        IF(iIniFieldType == 1) THEN
            !====================
            IF(TgFlowFlg) THEN
                CALL BC_TINLET_FLOW
            END IF
            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
            CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
            CALL BC_WALL_PR_io
            CALL BC_WALL_Q_io

            CALL VELO2MASSFLUX

            CALL Unified_MassFlux

            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
            CALL BC_WALL_G_io
            CALL BC_WALL_Q_io

            PhyTIME_io = 0.0_WP
            ITERG0_io = 0

        END IF

        IF(iIniFieldType == 0) THEN
            IF(TgFlowFlg) THEN
                CALL BC_TINLET_FLOW
            END IF
            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
            CALL INTFC_VARS1(NCL1S, NCL1E, NCL1S, NCL1E,PR_io)
            CALL BC_WALL_G_io
            CALL BC_WALL_PR_io

            CALL DENSITY_Staggered
            CALL MU_Staggered
            CALL VELOCITY_CALC_io

            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, G_io)
            CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
            CALL BC_WALL_G_io
            CALL BC_WALL_Q_io
        END IF

    END IF

    PhyTIME_io = 0.0_WP
    ITERG0_io = 0
    PhyTIME = PhyTIME_io
    ITERG0 = ITERG0_io

    DEALLOCATE (DUMMY)
    DEALLOCATE (QO)
    DEALLOCATE (QN)

    DEALLOCATE  (N2DOO)
    DEALLOCATE  (N3DOO)
    DEALLOCATE  (XNDO)
    DEALLOCATE  (YNDO)
    DEALLOCATE  (ZNDO)
    DEALLOCATE  (XCCO)
    DEALLOCATE  (YCCO)
    DEALLOCATE  (ZCCO)
    !DEALLOCATE  (ZCC)

    CALL CALL_TEC360

    CALL MPI_BARRIER(ICOMM, IERROR)
    IF(MYID == 0) CALL CHKRLHDL('    IO: END of interplating from a coarse mesh at time', MYID, PhyTIME_io)

    RETURN
END SUBROUTINE
