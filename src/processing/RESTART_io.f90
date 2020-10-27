!**********************************************************************************************************************************
!> @brief
!>        reStart
!> @details
!> subroutine: ReStart_INSTANT_VARS_io
!> subroutine: ReStart_AVERAGE_VARS_nonXperiodic_io
!> subroutine: ReStart_AVERAGE_VARS_THERMAL_nonXperiodic_io
!> subroutine: ReStart_AVERAGE_VARS_Xperiodic_io
!> subroutine: ReStart_AVERAGE_VARS_THERMAL_Xperiodic_io
!> subroutine: REStart_AVERAGE_SPECTRA
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE ReStart_INSTANT_VARS_io(TFM)  ! FOR BOTH KINDS OF DOMAINS....
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: TFM
    CHARACTER(LEN = 256) :: WRT_RST_FNM
    CHARACTER(LEN =64) :: FLNAME
    CHARACTER(15) :: PNTIM
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
    INTEGER(4) :: I, J, K, N
    REAL(WP) :: Gmean1_io, Umean1_io, Hmean1



    IF(MYID == 0) CALL CHKHDL('IO: ReStart instantanous flow field', MYID)

    ALLOCATE(DUMMY(NCL1E, N2DO(MYID), NCL3))
    WRITE(PNTIM, '(1ES15.9)') TFM !TimeReStart_io

    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        FLNAME  = 'DNS_perioz_INSTANT_T'
    ELSE
        FLNAME  = 'DNS_perixz_INSTANT_T'
    END IF

    DO N = 1, NDV + 1
        IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_U.D'
        IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_V.D'
        IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_W.D'
        IF(N == 4)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_P.D'

        CALL READ_3D_VARS(NCL1E, NCL2, NCL3, N2DO(MYID), JCL2G(1) - 1, ITERG0_io,PhyTIME_io, DUMMY, WRT_RST_FNM)

        IF(iIniFieldType == 0) THEN
            DO I = 1, NCL1E
                DO J = 1, N2DO(MYID)
                    DO K = 1, NCL3
                        IF(N == 1) G_io(I, J, K, 1) = DUMMY(I, J, K) !RHO U
                        IF(N == 2) G_io(I, J, K, 2) = DUMMY(I, J, K) !RHO V
                        IF(N == 3) G_io(I, J, K, 3) = DUMMY(I, J, K) !RHO W
                        IF(N == 4) PR_io(I, J, K) = DUMMY(I, J, K) !P
                    END DO
                END DO
            END DO

            IF(N == NFLOW) THEN
                CALL BULK_MASSFLUX_io(Gmean1_io)
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk massflux (reStart) = ', MYID, Gmean1_io)
            END IF

        END IF

        IF(iIniFieldType == 1) THEN
            DO I = 1, NCL1E
                DO J = 1, N2DO(MYID)
                    DO K = 1, NCL3
                        IF(N == 1) Q_io(I, J, K, 1) = DUMMY(I, J, K) !U
                        IF(N == 2) Q_io(I, J, K, 2) = DUMMY(I, J, K) !V
                        IF(N == 3) Q_io(I, J, K, 3) = DUMMY(I, J, K) !W
                        IF(N == 4) PR_io(I, J, K) = DUMMY(I, J, K) !P
                    END DO
                END DO
            END DO
            IF(N == NFLOW) THEN
                CALL BULK_VELOCITY_io(Umean1_io)
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk velocity (reStart) = ', MYID, Umean1_io)
            END IF
        END IF

    END DO


    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting instantanous flow field', MYID, PhyTIME_io)
    !============================ Thermal field==============================================================
    IF(iThermoDynamics == 1 .AND. iIniFieldType == 0 ) THEN
        IF(MYID == 0) CALL CHKHDL(' IO: ReStart instantanous thermal field', MYID)
        DO N = 1, 3
            IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_T.D'
            IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_D.D'
            IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // TRIM(FLNAME) // TRIM(PNTIM) // '_E.D'

            CALL READ_3D_VARS(NCL1E, NCL2, NCL3, N2DO(MYID), JCL2G(1) - 1, ITERG0_io,PhyTIME_io, DUMMY, WRT_RST_FNM)

            DO I = 1, NCL1E
                DO J = 1, N2DO(MYID)
                    DO K = 1, NCL3
                        IF(N == 1) TEMPERATURE(I, J, K) = DUMMY(I, J, K)
                        IF(N == 2) DENSITY    (I, J, K) = DUMMY(I, J, K)
                        IF(N == 3) ENTHALPY   (I, J, K) = DUMMY(I, J, K)
                    END DO
                END DO
            END DO

            IF(N == 3) THEN
                CALL BULK_H_io(Hmean1)
                IF(MYID == 0) CALL CHKRLHDL  ('         IO: The bulk H (reStart) =     ', MYID, Hmean1)
            END IF

        END DO
        IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting instantanous thermal field', MYID, PhyTIME_io)

        !=============BUILD UP DH ===========================================
        DO J = 1, N2DO(MYID)
            DO I = 1, NCL1E
                DO K = 1, NCL3
                    DH(I, J, K) = DENSITY(I, J, K) * ENTHALPY(I, J, K)
                END DO
            END DO
        END DO

        !============PrepARe DENSITY0 =====================================
        DENSITY0      (:, :, :) = DENSITY      (:, :, :)

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
            CALL IniField_THERMAL_io
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

    DEALLOCATE (DUMMY)

    IF(MYID == 0) CALL CHKRLHDL('===========PhyTIME_io ================== ', MYID, PhyTIME_io)
    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE ReStart_AVERAGE_VARS_nonXperiodic_io
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG, NARRAY
    INTEGER(4) :: IERR,SIZES_ARRAY(4), NEWTYPE,SUBSIZES(4),StartS(4)
    INTEGER(4) :: INISIZE, IRLSIZE, INTMPI(4), INTBYTE, DBLBYTE
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, I, J, L, L1, L2, N, M, H, P, NSZ3, NSZ4
    REAL(WP) :: RLEMPI(3)
    REAL(WP) :: RENL, DTL
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :, :)
    LOGICAL :: File_exists

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged flow field', MYID)

    NARRAY = 4

    NSZ3 = NDV + 1 + NDV + NDV + NDV + NDV * (7 - NDV) / 2 + NDV - 3 + &
    NDV * (7 - NDV) / 2 + NDV - 3 + NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8 + &
    NDV + NDV + (NDV - 1) * NDV + NDV
    NSZ4 = 1 + 1 + 1 + 1 + 1 + 1 + NDV +  NDV + (NDV - 1) * NDV + NDV

    ALLOCATE( DUMMY(1 : NCL1_io, 1 : N2DO(MYID), NSZ3, NSZ4 ) )
    DUMMY = 0.0_WP

    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perioz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

    SIZES_ARRAY(1) = NCL1_io
    SIZES_ARRAY(2) = NCL2
    SIZES_ARRAY(3) = NSZ3
    SIZES_ARRAY(4) = NSZ4

    SUBSIZES(1) = SIZES_ARRAY(1)
    SUBSIZES(2) = N2DO(MYID)
    SUBSIZES(3) = SIZES_ARRAY(3)
    SUBSIZES(4) = SIZES_ARRAY(4)

    StartS(1) = 0
    StartS(2) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    StartS(3) = 0
    StartS(4) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3) * SUBSIZES(4)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

    INISIZE = 4
    IRLSIZE = 3
    OFFSET = 0_MPI_OFFSET_KIND
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)

    CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)
    OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ_ALL(MYFILE, DUMMY,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_CLOSE(MYFILE, IERR)
    CALL MPI_TYPE_FREE(NEWTYPE, IERR)


    N1ML   = INTMPI(1)
    N2ML   = INTMPI(2)
    ITERG0_io = INTMPI(3)
    NSTATIS_io = INTMPI(4)

    PhyTIME_io = RLEMPI(1)
    RENL    = RLEMPI(2)
    DTL     = RLEMPI(3)

    IF(MYID == 0) THEN
        CALL CHKINTHDL('     I IN AVERAGED variables (I, J, L, M) ', MYID, N1ML)
        CALL CHKINTHDL('     J IN AVERAGED variables (I, J, L, M) ', MYID, N2ML)
        CALL CHKINTHDL('     ITERG0_io IN AVERAGED variables ', MYID, ITERG0_io)
        CALL CHKINTHDL('     NSTATIS_io IN AVERAGED variables', MYID, NSTATIS_io)

        CALL CHKRLHDL ('     PhyTIME_io IN AVERAGED variables', MYID, PhyTIME_io)
        CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RENL)
        CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, DTL)
    END IF

    L1 = NDV + 1
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                U1ztL_io(I, J, L) = DUMMY(I, J, L, 1)
            END DO
        END DO
    END DO

    L1 = NDV
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                G1ztL_io(I, J, L) = DUMMY(I, J, L, 2)
            END DO
        END DO
    END DO

    L1 = NDV
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                UPztL_io(I, J, L) = DUMMY(I, J, L, 3)
            END DO
        END DO
    END DO

    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                U2ztL_io(I, J, L) = DUMMY(I, J, L, 4)
            END DO
        END DO
    END DO

    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                UGztL_io(I, J, L) = DUMMY(I, J, L, 5)
            END DO
        END DO
    END DO

    L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                UGUztL_io(I, J, L) = DUMMY(I, J, L,6)
            END DO
        END DO
    END DO


    L1 = NDV
    L2 = NDV
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                DO M = 1, L2
                    DVDL1ztL_io(I, J, L, M) = DUMMY(I, J, L, M +6)
                END DO
            END DO
        END DO
    END DO

    L1 = NDV
    L2 = NDV
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                DO M = 1, L2
                    DVDLPztL_io(I, J, L, M) = DUMMY(I, J, L, M +9)
                END DO
            END DO
        END DO
    END DO

    DO M = 1, NDV
        DO H = 1, NDV
            L1 = (M - 1) * NDV + H
            DO N = 1, NDV
                DO P = 1, NDV
                    L2 = (N - 1) * NDV + P !
                    DO J = 1, N2DO(MYID)
                        DO I = 1, NCL1_io
                            DVDL2ztL_io(I, J, L1, L2) = DUMMY(I, J, L1, L2 + 12)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO

    DEALLOCATE (DUMMY)

    IF(iIniFieldType == 1) THEN
        PhyTIME_io = 0.0_WP
        ITERG0_io = 0
    END IF


    IF(iThermoDynamics == 1 .AND. iIniFieldType == 0 ) CALL ReStart_AVERAGE_VARS_THERMAL_nonXperiodic_io

    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting averaged flow field', MYID, PhyTIME_io)

    RETURN

END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE ReStart_AVERAGE_VARS_THERMAL_nonXperiodic_io
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG, NARRAY
    INTEGER(4) :: IERR,SIZES_ARRAY(3), NEWTYPE,SUBSIZES(3),StartS(3)
    INTEGER(4) :: INISIZE, IRLSIZE, INTMPI(4), INTBYTE, DBLBYTE
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, I, J, L1, L, N
    REAL(WP) :: RLEMPI(3)
    REAL(WP) :: RENL, DTL
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
    LOGICAL :: File_exists

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged thermal field', MYID)

    NARRAY = 3
    ALLOCATE( DUMMY(1 : NCL1_io, 1 : N2DO(MYID), 7+ 2 * NDV + 4) )
    DUMMY = 0.0_WP

    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perioz_AVERAGD_T' // TRIM(PNTIM) // '_THEL.D'

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)


    SIZES_ARRAY(1) = NCL1_io
    SIZES_ARRAY(2) = NCL2
    SIZES_ARRAY(3) = 7+ 2 * NDV + 4

    SUBSIZES(1) = SIZES_ARRAY(1)
    SUBSIZES(2) = N2DO(MYID)
    SUBSIZES(3) = SIZES_ARRAY(3)

    StartS(1) = 0
    StartS(2) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    StartS(3) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

    INISIZE = 4
    IRLSIZE = 3
    OFFSET = 0_MPI_OFFSET_KIND
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

    OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ_ALL(MYFILE, DUMMY,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_CLOSE(MYFILE, IERR)
    CALL MPI_TYPE_FREE(NEWTYPE, IERR)


    N1ML   = INTMPI(1)
    N2ML   = INTMPI(2)
    ITERG0_io = INTMPI(3)
    NSTATIS_io = INTMPI(4)

    PhyTIME_io = RLEMPI(1)
    RENL    = RLEMPI(2)
    DTL     = RLEMPI(3)

    N  = 1
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            T1ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO


    N  = 2
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            D1ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO

    N  = 3
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            H1ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO

    N  = 4
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            T2ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO


    N  = 5
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            D2ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO

    N  = 6
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            H2ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO

    N  = 7
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DHztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO

    L1 = NDV
    N  = 7
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                UHztL_io(I, J, L) = DUMMY(I, J, L+ N)
            END DO
        END DO
    END DO

    L1 = NDV
    N  = 7+ NDV
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            DO L = 1, L1
                GHztL_io(I, J, L) = DUMMY(I, J, L+ N)
            END DO
        END DO
    END DO


    N  = 14
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            K1ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO
    N  = 15
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            K2ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO
    N  = 16
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            M1ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO
    N  = 17
    DO I = 1, NCL1_io
        DO J = 1, N2DO(MYID)
            M2ztL_io(I, J) = DUMMY(I, J, N)
        END DO
    END DO


    DEALLOCATE (DUMMY)

    IF(iIniFieldType == 1) THEN
        PhyTIME_io = 0.0_WP
        ITERG0_io = 0
    END IF

    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting averaged thermal field', MYID, PhyTIME_io)

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE ReStart_AVERAGE_VARS_Xperiodic_io
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG
    INTEGER(4) :: IERR, NEWTYPE
    INTEGER(4) :: INISIZE, IRLSIZE, INTBYTE, DBLBYTE
    INTEGER(4) :: NARRAY, NSZ
    INTEGER(4), ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),StartS(:)
    INTEGER(4), ALLOCATABLE :: INTMPI(:)
    REAL(WP), ALLOCATABLE :: RLEMPI(:)
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, J, L, L1, L2, N, NLLL, NNNL, M, K
    REAL(WP) :: RENL, DTL
    REAL(WP), ALLOCATABLE :: DUMMY(:, :), DUMMY1(:, :), DUMMY2(:, :)
    LOGICAL :: File_exists
    INTEGER(4) :: NSTATIS_io1, NSTATIS_io2

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged flow field', MYID)

    !============================
    NARRAY = 2
    ALLOCATE ( SIZES_ARRAY(NARRAY) )
    ALLOCATE ( SUBSIZES   (NARRAY) )
    ALLOCATE ( StartS     (NARRAY) )
    INISIZE = 4
    IRLSIZE = 3
    ALLOCATE ( INTMPI(INISIZE)       )
    ALLOCATE ( RLEMPI(IRLSIZE)       )

    ! original
    NSZ = NDV + 1 + 2 * NDV + 2 * (NDV * (7 - NDV) / 2 + NDV - 3) + &
    (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) + &
    2 * NDV * NDV + ((NDV - 1) * 3 + NDV) * ((NDV - 1) * 3 + NDV)

    ! IF with U3
    NSZ = NSZ + (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)

    ! IF with quadrant
    IF(iPPQuadrants == 1) NSZ = NSZ + 6* (4* QUADHN)

    ! IF with driven force
    NSZ = NSZ + NDV + 1

    ! IF with octant
    !NSZ = NSZ + 12 * (8* QUADHN)


    ALLOCATE( DUMMY (1 : N2DO(MYID), NSZ ) )
    ALLOCATE( DUMMY1(1 : N2DO(MYID), NSZ ) )
    DUMMY  = 0.0_WP
    DUMMY1 = 0.0_WP


    !=======prepARe sizes and rangeS ===================================================
    SIZES_ARRAY(1) = NCL2
    SIZES_ARRAY(2) = NSZ

    SUBSIZES(1) = N2DO(MYID)
    SUBSIZES(2) = SIZES_ARRAY(2)

    StartS(1) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    StartS(2) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2)

    !=========== Read data1 ==============================================================
    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged flow field ' // TRIM(WRT_AVE_FNM_io), MYID)

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

    OFFSET = 0_MPI_OFFSET_KIND
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)

    CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)
    OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ_ALL(MYFILE, DUMMY1,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_CLOSE(MYFILE, IERR)
    CALL MPI_TYPE_FREE(NEWTYPE, IERR)


    N2ML     = INTMPI(1)
    NNNL     = INTMPI(2)
    ITERG0_io = INTMPI(3)
    NSTATIS_io1 = INTMPI(4)

    PhyTIME_io = RLEMPI(1)
    RENL     = RLEMPI(2)
    DTL      = RLEMPI(3)

    IF(MYID == 0) THEN
        CALL CHK2INTHDL('     Expected and reaD -in DUMMY Sizes ARe', MYID, NSZ, NNNL)
        IF(NSZ /= NNNL) CALL ERRHDL('ReStarting averaged data formats are not consistent! ', MYID)
        CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, N2ML)
        CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, NNNL)
        CALL CHKINTHDL('     ITERG0_io IN AVERAGED variables ', MYID, ITERG0_io)
        CALL CHKINTHDL('     NSTATIS_io IN AVERAGED variables', MYID, NSTATIS_io1)

        CALL CHKRLHDL ('     PhyTIME_io IN AVERAGED variables', MYID, PhyTIME_io)
        CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RENL)
        CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, DTL)
    END IF

    NSTATIS_io = NSTATIS_io1
    DUMMY = DUMMY1

    !=========== Read data2 ==============================================================
    IF(tRunAve_Reset < PhyTIME_io .AND. tRunAve_Reset > tRunAve1) THEN
        ALLOCATE( DUMMY2(1 : N2DO(MYID), NSZ ) )
        DUMMY2 = 0.0_WP

        DFLG = 101
        WRITE(PNTIM, '(1ES15.9)') tRunAve_Reset
        WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

        IF(MYID == 0) CALL CHKHDL(' IO: Resetting averaged flow field ' // TRIM(WRT_AVE_FNM_io), MYID)

        INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
        IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY,SUBSIZES,StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

        OFFSET = 0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)

        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)
        OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_READ_ALL(MYFILE, DUMMY2,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_CLOSE(MYFILE, IERR)
        CALL MPI_TYPE_FREE(NEWTYPE, IERR)


        NSTATIS_io2 = INTMPI(4)
        IF(MYID == 0) THEN
            IF(NSZ /= INTMPI(2)) CALL ERRHDL('ReStarting averaged data formats are not consistent! ', MYID)
            CALL CHK2INTHDL('     Expected and reaD -in DUMMY Sizes ARe', MYID, NSZ, INTMPI(2))
            CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, INTMPI(1))
            CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, INTMPI(2))
            CALL CHKINTHDL('     ITERG0_io IN AVERAGED variables ', MYID, INTMPI(3))
            CALL CHKINTHDL('     NSTATIS_io IN AVERAGED variables', MYID, NSTATIS_io2)

            CALL CHKRLHDL ('     PhyTIME_io IN AVERAGED variables', MYID, RLEMPI(1))
            CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RLEMPI(2))
            CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, RLEMPI(3))
        END IF

        DUMMY = (DUMMY1 * DBLE(NSTATIS_io1) - DUMMY2 * DBLE(NSTATIS_io2)) / DBLE(NSTATIS_io1-NSTATIS_io2)
        NSTATIS_io = NSTATIS_io1-NSTATIS_io2
        DEALLOCATE (DUMMY2)
    END IF

    !==========================================
    N = 0
    L1 = NDV + 1
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U1xztL_io(J, L) = DUMMY(J, N +L) !1-4
        END DO
    END DO
    N = N + L1 !4
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in U1xztL_io', MYID, L1, N)

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            G1xztL_io(J, L) = DUMMY(J, N +L)!5-7
        END DO
    END DO
    N = N + L1 !7
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in G1xztL_io', MYID, L1, N)

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            UPxztL_io(J, L) = DUMMY(J, N +L) !8- 10
        END DO
    END DO
    N = N + L1 !10
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in UPxztL_io', MYID, L1, N)

    !==========================================
    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U2xztL_io(J, L) = DUMMY(J, N +L) !11- 16
        END DO
    END DO
    N = N + L1 !16
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in U2xztL_io', MYID, L1, N)

    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            UGxztL_io(J, L) = DUMMY(J, N +L)
        END DO
    END DO
    N = N + L1 !22
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in UGxztL_io', MYID, L1, N)


    L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8 !10
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            UGUxztL_io(J, L) = DUMMY(J, N +L)
        END DO
    END DO
    N = N + L1 !32
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in UGUxztL_io', MYID, L1, N)

    !!!!        ! Commented below IS no U3 ===========================
    L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8 !10
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U3xztL_io(J, L) = DUMMY(J, N +L)
        END DO
    END DO
    N = N + L1 !32
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in U3xztL_io', MYID, L1, N)

    !==========================================
    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDL1xztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in DVDL1xztL_io', MYID, L1, N)


    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDLPxztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in DVDLPxztL_io', MYID, L1, N)

    L1 = (NDV - 1) * 3 + NDV !9
    L2 = (NDV - 1) * 3 + NDV !9
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDL2xztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in DVDL2xztL_io', MYID, L1, N)

    !        !        ! Commented below IS no quadrant
    !============= QuadranT ===============
    IF(iPPQuadrants == 1)  THEN
        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADUVxztL_io(J, L, K) = DUMMY(J, M)
                    !WRITE(*, *) QUADUVxztL_io(J, L, K) !!!test
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADUVxztL_io', MYID, L1, N)

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADVzxztL_io(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADVzxztL_io', MYID, L1, N)

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADTKxztL_io(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADTKxztL_io', MYID, L1, N)

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADDRxztL_io(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADDRxztL_io', MYID, L1, N)

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADDUV1xztL_io(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADDUV1xztL_io', MYID, L1, N)

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADDUV2xztL_io(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in QUADDUV2xztL_io', MYID, L1, N)
    END IF
    !=========================================================
    ! with driven force
    L1 = NDV + 1
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            FUxztL_io(J, L) = DUMMY(J, N +L) !1-4
        END DO
    END DO
    N = N + L1 !4
    !IF(MYID == 0) CALL CHK2INTHDL('     Reading in FUxztL_io', MYID, L1, N)

    !        !=======================================================
    !        !=============octanT ==\rho=============
    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDUVxztL_io(J, L, K) = DUMMY(J, M)
    !                    !WRITE(*, *) OCTDUVxztL_io(J, L, K) !!!test
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDUVxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDVzxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDVzxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDTKxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDTKxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDDRxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDDRxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDDUV1xztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDDUV1xztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTDDUV2xztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTDDUV2xztL_io', MYID, L1, N)

    !        !=============octanT == T =============
    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTUVxztL_io(J, L, K) = DUMMY(J, M)
    !                    !WRITE(*, *) OCTTUVxztL_io(J, L, K) !!!test
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTUVxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTVzxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTVzxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTTKxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTTKxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTDRxztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTDRxztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTDUV1xztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTDUV1xztL_io', MYID, L1, N)

    !        L1 = 8
    !        L2 = QUADHN
    !        DO J = 1, N2DO(MYID)
    !            DO L = 1, L1 !
    !                DO K = 1, L2
    !                    M = N + (L- 1) * L2 + K
    !                    OCTTDUV2xztL_io(J, L, K) = DUMMY(J, M)
    !                END DO
    !            END DO
    !        END DO
    !        N = N + L1 * L2
    !        !IF(MYID == 0) CALL CHK2INTHDL('     Reading in OCTTDUV2xztL_io', MYID, L1, N)



    IF(N /= NSZ) THEN
        CALL ERRHDL('#ERROR in ReStart_AVERAGE_VARS_Xperiodic_io', MYID)
    END IF

    DEALLOCATE (DUMMY)
    DEALLOCATE (DUMMY1)
    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting averaged flow    field', MYID, PhyTIME_io)

    IF(iThermoDynamics == 1 .AND. iIniFieldType == 0 ) CALL ReStart_AVERAGE_VARS_THERMAL_Xperiodic_io

    IF(iIniFieldType == 0 )  CALL REStart_AVERAGE_SPECTRA

    RETURN

END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE ReStart_AVERAGE_VARS_THERMAL_Xperiodic_io
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG
    INTEGER(4) :: IERR, NEWTYPE
    INTEGER(4) :: INISIZE, IRLSIZE, INTBYTE, DBLBYTE
    INTEGER(4) :: NARRAY, NSZ
    INTEGER(4), ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),StartS(:)
    INTEGER(4), ALLOCATABLE :: INTMPI(:)
    REAL(WP), ALLOCATABLE :: RLEMPI(:)
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, J, L1, L2, L3, L, N, NNNL, K, S, M
    REAL(WP) :: RENL, DTL
    REAL(WP), ALLOCATABLE :: DUMMY(:, :), DUMMY1(:, :), DUMMY2(:, :)
    LOGICAL :: File_exists
    INTEGER(4) :: NSTATIS_io1, NSTATIS_io2




    !============================
    NARRAY = 2
    ALLOCATE ( SIZES_ARRAY(NARRAY) )
    ALLOCATE ( SUBSIZES   (NARRAY) )
    ALLOCATE ( StartS     (NARRAY) )
    INISIZE = 4
    IRLSIZE = 3
    ALLOCATE ( INTMPI(INISIZE)       )
    ALLOCATE ( RLEMPI(IRLSIZE)       )


    !!NSZ = 9+5* NDV + ((NDV * (7 - NDV)) / 2 + NDV - 3) + 4* NDV * NDV +3 * NDV * NDV * NDV
    NSZ = 9 + 5 * NDV + ((NDV * (7 - NDV)) / 2 + NDV - 3) +3 * NDV * NDV +3 * NDV * NDV * NDV + &
            ((NDV - 1) * NDV + NDV) * ((NDV - 1) * NDV + NDV)
    ALLOCATE( DUMMY (1 : N2DO(MYID), NSZ ) )
    ALLOCATE( DUMMY1(1 : N2DO(MYID), NSZ ) )
    DUMMY  = 0.0_WP
    DUMMY1 = 0.0_WP

    !============================
    SIZES_ARRAY(1) = NCL2
    SIZES_ARRAY(2) = NSZ

    SUBSIZES(1) = N2DO(MYID)
    SUBSIZES(2) = SIZES_ARRAY(2)

    StartS(1) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    StartS(2) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2)

    INISIZE = 4
    IRLSIZE = 3
    !============================

    !!=========== Read data1 ==============================================================
    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_THEL.D'

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged thermal field' //WRT_AVE_FNM_io, MYID)

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)


    OFFSET = 0_MPI_OFFSET_KIND
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

    OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_READ_ALL(MYFILE, DUMMY1,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_CLOSE(MYFILE, IERR)
    CALL MPI_TYPE_FREE(NEWTYPE, IERR)


    N2ML     = INTMPI(1)
    NNNL     = INTMPI(2)
    ITERG0_io = INTMPI(3)
    NSTATIS_io1 = INTMPI(4)

    PhyTIME_io = RLEMPI(1)
    RENL    = RLEMPI(2)
    DTL     = RLEMPI(3)

    IF(MYID == 0) THEN
        CALL CHK2INTHDL('     Expected and reaD -in DUMMY Sizes ARe', MYID, NSZ, NNNL)
        CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, N2ML)
        CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, NNNL)
        CALL CHKINTHDL('     ITERG0_io IN AVERAGED variables ', MYID, ITERG0_io)
        CALL CHKINTHDL('     NSTATIS_io IN AVERAGED variables', MYID, NSTATIS_io1)

        CALL CHKRLHDL ('     PhyTIME_io IN AVERAGED variables', MYID, PhyTIME_io)
        CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RENL)
        CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, DTL)
    END IF

    NSTATIS_io = NSTATIS_io1
    DUMMY = DUMMY1

    !!=========== Read data2 ==============================================================
    IF(tRunAve_Reset < PhyTIME_io) THEN
        ALLOCATE( DUMMY2 (1 : N2DO(MYID), NSZ ) )
        DUMMY2 = 0.0_WP

        DFLG = 101
        WRITE(PNTIM, '(1ES15.9)') tRunAve_Reset
        WRITE(WRT_AVE_FNM_io, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_THEL.D'

        IF(MYID == 0) CALL CHKHDL(' IO: Resetting averaged thermal field' //WRT_AVE_FNM_io, MYID)

        INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
        IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY,SUBSIZES,StartS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)


        OFFSET = 0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

        OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_READ_ALL(MYFILE, DUMMY2,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_CLOSE(MYFILE, IERR)
        CALL MPI_TYPE_FREE(NEWTYPE, IERR)

        NSTATIS_io2 = INTMPI(4)

        IF(MYID == 0) THEN
            CALL CHK2INTHDL('     Expected and reaD -in DUMMY Sizes ARe', MYID, NSZ, INTMPI(2))
            CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, INTMPI(1))
            CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, INTMPI(2))
            CALL CHKINTHDL('     ITERG0_io IN AVERAGED variables ', MYID, INTMPI(3))
            CALL CHKINTHDL('     NSTATIS_io IN AVERAGED variables', MYID, NSTATIS_io2)

            CALL CHKRLHDL ('     PhyTIME_io IN AVERAGED variables', MYID, RLEMPI(1))
            CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RLEMPI(2))
            CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, RLEMPI(3))
        END IF

        DUMMY = (DUMMY1 * DBLE(NSTATIS_io1) - DUMMY2 * DBLE(NSTATIS_io2)) / DBLE(NSTATIS_io1-NSTATIS_io2)
        NSTATIS_io = NSTATIS_io1-NSTATIS_io2
        DEALLOCATE (DUMMY2)
    END IF

    !==================================================

    N  = 1
    DO J = 1, N2DO(MYID)
        T1xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !2

    DO J = 1, N2DO(MYID)
        D1xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !3

    DO J = 1, N2DO(MYID)
        H1xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !4

    DO J = 1, N2DO(MYID)
        M1xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !5

    !===============================
    DO J = 1, N2DO(MYID)
        T2xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !6

    DO J = 1, N2DO(MYID)
        D2xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !7

    DO J = 1, N2DO(MYID)
        H2xztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !8

    !=============================
    DO J = 1, N2DO(MYID)
        DHxztL_io(J) = DUMMY(J, N)
    END DO
    N  = N + 1 !9

    DO J = 1, N2DO(MYID)
        PHxztL_io(J) = DUMMY(J, N)
    END DO
    N  = N +0 !9

    !===============================
    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            UHxztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 12

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            GHxztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 15

    L1 = (NDV * (7 - NDV)) / 2 + NDV - 3 !6
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 ! 16, 17, 18, 19, 20, 21
            U2DHxztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 21

    !==========================================
    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 ! 22, 23, 24
            DhDL1xztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 24

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 ! 25, 26, 27
            DhDLPxztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 27

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 ! 28, 29, 30
            DTDLKxztL_io(J, L) = DUMMY(J, L+ N)
        END DO
    END DO
    N  = N +L1 ! 30

    !==========================================
    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K  !31- 33, 34- 36, 37- 39
                DVDL1MxztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N +L1 * L2 !39

    L1 = (NDV - 1) * NDV + NDV
    L2 = (NDV - 1) * NDV + NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K  !40-42, 43-45, 46-48
                DVDL2MxztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N +L1 * L2 !48

    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K  !49- 51, 52 - 54, 55- 57
                DVDL1MHxztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N +L1 * L2 !57

    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K  !58-60,61-63,64-66
                DTDLKUxztL_io(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N +L1 * L2 !66

    !==========================================
    L1 = NDV
    L2 = NDV
    L3 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                DO S = 1, L3
                    M = N + (L- 1) * L2 * L3 + (K - 1) * L3 +S !67-75,76-84,85-93
                    DVDL1MUxztL_io(J, L, K, S) = DUMMY(J, M)
                END DO
            END DO
        END DO
    END DO
    N = N +L1 * L2 * L3 !93

    L1 = NDV
    L2 = NDV
    L3 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                DO S = 1, L3
                    M = N + (L- 1) * L2 * L3 + (K - 1) * L3 +S !
                    DTDLKDVDLxztL_io(J, L, K, S) = DUMMY(J, M)
                END DO
            END DO
        END DO
    END DO
    N = N +L1 * L2 * L3 !120

    L1 = NDV
    L2 = NDV
    L3 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                DO S = 1, L3
                    M = N + (L- 1) * L2 * L3 + (K - 1) * L3 +S !67-75,76-84,85-93
                    DHDLMDVDLxztL_io(J, L, K, S) = DUMMY(J, M)
                END DO
            END DO
        END DO
    END DO
    N = N +L1 * L2 * L3 !219

    IF(N /= NSZ) THEN
        CALL ERRHDL('ERROR in ReStart_AVERAGE_VARS_THERMAL_Xperiodic_io', MYID)
    END IF

    DEALLOCATE (DUMMY)
    DEALLOCATE (DUMMY1)

    IF(iIniFieldType == 1) THEN
        PhyTIME_io = 0.0_WP
        ITERG0_io = 0
    END IF

    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting averaged thermal field', MYID, PhyTIME_io)

    RETURN

END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE REStart_AVERAGE_SPECTRA
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    CHARACTER(64) :: FIL1
    LOGICAL :: File_exists
    INTEGER(4) :: DFLG = 10
    INTEGER(4) :: N1ML, N2DOL, N3ML


    REAL(WP) :: RENL
    REAL(WP) :: DTL
    INTEGER(4) :: ITEMP
    INTEGER(4) :: IOS

    !IF(iPPSpectra  /= 1) RETURN
    IF(MYID  /= 0)  RETURN


    WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
    WRITE(WRT_AVE_FNM_io, '(A, I4.4, A4)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_SPEC.D'

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_iO // ' does not exist.', 0)

    !        !==============OPEN FILE =============================================
    OPEN(DFLG, FILE = TRIM(ADJUSTL(WRT_AVE_FNM_io)), FORM = 'UNFORMATTED',STATUS = 'OLD', IOSTAT = IOS)

    READ(DFLG) ITEMP, ITEMP
    READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_io
    READ(DFLG) NSTATIS_io
    READ(DFLG) PhyTIME_io,RENL, DTL


    !================ Velocity =====================
    READ(DFLG) R11X1_xztLa
    READ(DFLG) R22X1_xztLa
    READ(DFLG) R33X1_xztLa
    READ(DFLG) R12X1_xztLa
    READ(DFLG) R13X1_xztLa
    READ(DFLG) R23X1_xztLa

    READ(DFLG) R11X3_xztLa
    READ(DFLG) R22X3_xztLa
    READ(DFLG) R33X3_xztLa
    READ(DFLG) R12X3_xztLa
    READ(DFLG) R13X3_xztLa
    READ(DFLG) R23X3_xztLa

    READ(DFLG) ENE11T_xztLa
    READ(DFLG) ENE22T_xztLa
    READ(DFLG) ENE33T_xztLa
    READ(DFLG) ENE12T_xztLa
    READ(DFLG) ENE13T_xztLa
    READ(DFLG) ENE23T_xztLa

    READ(DFLG) ENE11Z_xztLa
    READ(DFLG) ENE22Z_xztLa
    READ(DFLG) ENE33Z_xztLa
    READ(DFLG) ENE12Z_xztLa
    READ(DFLG) ENE13Z_xztLa
    READ(DFLG) ENE23Z_xztLa

    !================ VoritICity ====================
    READ(DFLG) V11X1_xztLa
    READ(DFLG) V22X1_xztLa
    READ(DFLG) V33X1_xztLa
    READ(DFLG) V12X1_xztLa
    READ(DFLG) V13X1_xztLa
    READ(DFLG) V23X1_xztLa

    READ(DFLG) V11X3_xztLa
    READ(DFLG) V22X3_xztLa
    READ(DFLG) V33X3_xztLa
    READ(DFLG) V12X3_xztLa
    READ(DFLG) V13X3_xztLa
    READ(DFLG) V23X3_xztLa

    READ(DFLG) ENV11T_xztLa
    READ(DFLG) ENV22T_xztLa
    READ(DFLG) ENV33T_xztLa
    READ(DFLG) ENV12T_xztLa
    READ(DFLG) ENV13T_xztLa
    READ(DFLG) ENV23T_xztLa

    READ(DFLG) ENV11Z_xztLa
    READ(DFLG) ENV22Z_xztLa
    READ(DFLG) ENV33Z_xztLa
    READ(DFLG) ENV12Z_xztLa
    READ(DFLG) ENV13Z_xztLa
    READ(DFLG) ENV23Z_xztLa


    !=============== VoritICity & VelocitY =========================
    READ(DFLG) VO11X1_xztLa
    READ(DFLG) VO12X1_xztLa
    READ(DFLG) VO13X1_xztLa

    READ(DFLG) VO21X1_xztLa
    READ(DFLG) VO22X1_xztLa
    READ(DFLG) VO23X1_xztLa

    READ(DFLG) VO31X1_xztLa
    READ(DFLG) VO32X1_xztLa
    READ(DFLG) VO33X1_xztLa

    READ(DFLG) VO11X3_xztLa
    READ(DFLG) VO12X3_xztLa
    READ(DFLG) VO13X3_xztLa

    READ(DFLG) VO21X3_xztLa
    READ(DFLG) VO22X3_xztLa
    READ(DFLG) VO23X3_xztLa

    READ(DFLG) VO31X3_xztLa
    READ(DFLG) VO32X3_xztLa
    READ(DFLG) VO33X3_xztLa

    READ(DFLG) EVO11T_xztLa
    READ(DFLG) EVO12T_xztLa
    READ(DFLG) EVO13T_xztLa

    READ(DFLG) EVO21T_xztLa
    READ(DFLG) EVO22T_xztLa
    READ(DFLG) EVO23T_xztLa

    READ(DFLG) EVO31T_xztLa
    READ(DFLG) EVO32T_xztLa
    READ(DFLG) EVO33T_xztLa

    READ(DFLG) EVO11Z_xztLa
    READ(DFLG) EVO12Z_xztLa
    READ(DFLG) EVO13Z_xztLa

    READ(DFLG) EVO21Z_xztLa
    READ(DFLG) EVO22Z_xztLa
    READ(DFLG) EVO23Z_xztLa

    READ(DFLG) EVO31Z_xztLa
    READ(DFLG) EVO32Z_xztLa
    READ(DFLG) EVO33Z_xztLa

    CLOSE(DFLG)
    IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(WRT_AVE_FNM_io), MYID)



    RETURN
END SUBROUTINE


!!***********************************************************************************
!SUBROUTINE REStart_INSTANT_GLOBAL_io
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE thermal_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1

!        !INTEGER(4) :: ITIM
!        INTEGER(4) :: I, J, K, JJ
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML

!        INTEGER(4) :: IOS

!        REAL(WP), ALLOCATABLE :: U_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: V_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: W_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: P_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: T_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: D_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: H_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: D0_F0_io(:, :, :)
!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_io

!        ALLOCATE( U_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!        ALLOCATE( W_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!        ALLOCATE( V_F0_io(0 : NCL1_io + 1, 1 : NND2, 1 : NCL3) )
!        ALLOCATE( P_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )


!        IF(iIniFieldType == 0)THEN
!            ALLOCATE( T_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!            ALLOCATE( D_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!            ALLOCATE( H_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!            ALLOCATE( D0_F0_io(0 : NCL1_io + 1, 1 : NCL2, 1 : NCL3) )
!        END IF

!        !=======================================U=============================================
!        FIL1 = 'DNS_inouDOmAIn_RST_' // TRIM(PNTIM) // '_FLOW_INST.bin'
!        IF(MYID == 0) CALL CHKHDL( '    Start reading data from ' //FIL1, MYID)

!        OPEN(10, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF

!        READ(10) N1ML, N2DOL, N3ML, ITERG0_io
!        READ(10) PhyTIME_io,RENL, DTL
!        IF(MYID == 0) THEN
!            CALL CHKRLHDL   ('        UVW field TimE = ', MYID, PhyTIME_io)
!            CALL CHKINTHDL  ('        UVW field ITERG=   ', MYID, ITERG0_io)
!        END IF

!        READ(10) U_F0_io
!        READ(10) V_F0_io
!        READ(10) W_F0_io
!        READ(10) P_F0_io
!        CLOSE(10)
!        IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(FIL1), MYID)
!        IF(DABS(PhyTIME_io- TimeReStart_io) >TSAVE1) CALL CHKHDL('WARning: ReStarting from an ealier time step!', MYID)

!!====================================================================================
!        DO J = 1, N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I = 0, NCL1_io + 1, 1
!                DO K = 1, NCL3
!                    G_io(I, J, K, 1) = (U_F0_io(I, JJ, K))
!                    G_io(I, J, K, 2) = (V_F0_io(I, JJ, K))
!                    G_io(I, J, K, 3) = (W_F0_io(I, JJ, K))
!                    PR_io(I, J, K)   = (P_F0_io(I, JJ, K))
!                END DO
!            END DO
!        END DO

!        IF (MYID == NPSLV) THEN
!            DO I = 0, NCL1_io + 1, 1
!                DO K = 1, NCL3
!                    !Q_io(I, N2DO(MYID) + 1, K, 2) = (V_F0_io(I, NND2, K))
!                    G_io(I, N2DO(MYID) + 1, K, 2) = 0.0_WP
!                END DO
!            END DO
!        ENDIF

!        IF (MYID == 0) THEN
!            DO I = 0, NCL1_io + 1, 1
!                DO K = 1, NCL3
!                    G_io(I, 1, K, 2) = 0.0_WP
!                    G_io(I, 0, K, 2) = 0.0_WP
!                END DO
!            END DO
!        ENDIF
!        !====================================================================================
!        IF( iIniFieldType == 0 ) THEN  !reStart thermal field as well
!            FIL1 = 'DNS_inouDOmAIn_RST_' // TRIM(PNTIM) // '_THEL_INST.bin'
!            IF(MYID == 0) CALL CHKHDL( '    Start reading ' //FIL1, MYID)

!            OPEN(10, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!            IF(IOS/= 0) THEN
!                CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!            END IF

!            READ(10) N1ML, N2DOL, N3ML, ITERG0_io
!            READ(10) PhyTIME_io,RENL, DTL
!            IF(MYID == 0) THEN
!                CALL CHKRLHDL   ('        TDH field TimE = ', MYID, PhyTIME_io)
!                CALL CHKINTHDL  ('        TDH field ITERG=   ', MYID, ITERG0_io)
!            END IF
!            READ(10) T_F0_io
!            READ(10) D_F0_io
!            READ(10) H_F0_io
!            READ(10) D0_F0_io
!            CLOSE(10)
!            IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(FIL1), MYID)
!            IF(DABS(PhyTIME_io- TimeReStart_io) >TSAVE1) CALL CHKHDL('WARning: ReStarting from an ealier time step!', MYID)

!            !==============================================================================
!            DO J = 1, N2DO(MYID)
!                JJ = JCL2G(J)
!                DO I = 0, NCL1_io + 1, 1
!                    DO K = 1, NCL3
!                        TEMPERATURE(I, J, K) = (T_F0_io(I, JJ, K))
!                        DENSITY(I, J, K) = (D_F0_io(I, JJ, K))
!                        ENTHALPY(I, J, K) = (H_F0_io(I, JJ, K))
!                        DENSITY0(I, J, K) = (D0_F0_io(I, JJ, K))
!                    END DO
!                END DO
!            END DO

!            !=============BUILD UP DH ============
!            DO J = 1, N2DO(MYID)
!                DO I = 0, NCL1_io + 1, 1
!                    DO K = 1, NCL3
!                        DH(I, J, K) = DENSITY(I, J, K) * ENTHALPY(I, J, K)
!                    END DO
!                END DO
!            END DO

!            !===========build up other variables =======
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML
!            CALL THERM_PROP_UPDATE(IALLDM)
!            CALL INTFC_OUL_THERMAL_io
!            CALL INTFC_MFD_THERMAL_io
!            IF(iThermalWallType == BC_Fixed_Heat_Flux ) CALL BC_WALL_Isoflux(IALLDM)
!            CALL VELOCITY_CALC_io
!            !===========build up mass fluX =========
!            !CALL Massflux_reStart_io
!        END IF

!        IF(iIniFieldType == 1) THEN

!            CALL IniField_THERMAL_io
!            Q_io = G_io   ! based on constant unit DENSITY
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML

!            IF(iThermoDynamics == 1)THEN
!                CALL SOLVERRK3_ENG_io(0)
!                PhyTIME_io = 0.0_WP
!                ITERG0_io = 0
!            END IF

!        END IF

!!        IF(MYID == 0) THEN
!!           DO J = 1, N2DO(MYID)
!!              DO I = 0, NCL1_io + 1, 1
!!                 DO K = 1, ncl3
!!                    WRITE(*, '(3I5.1,7ES13.5)') J, I, K, U_F0_io(I, J, K), V_F0_io(I, J, K),w_F0_io(I, J, K), P_F0_io(I, J, K), &
!!                    T_F0_io(I, J, K), D_F0_io(I, J, K),h_F0_io(I, J, K)
!!                 END DO
!!              END DO
!!            END DO
!!        END IF

!!        DT = DTL
!!        REN = RENL

!        DEALLOCATE (U_F0_io,W_F0_io,P_F0_io, V_F0_io)
!        IF( iIniFieldType == 0 ) DEALLOCATE (T_F0_io, D_F0_io, H_F0_io)

!        CALL MPI_BARRIER(ICOMM, IERROR)

!        RETURN
!    END SUBROUTINE

!!***********************************************************************************
!SUBROUTINE REStart_INSTANT_LOCAL_io
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE thermal_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: I, J, K
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML
!        INTEGER(4) :: DFLG = 10
!        INTEGER(4) :: MYID_tmp, NPTOT_tmp

!        INTEGER(4) :: IOS

!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        !===============CREAT FILE NAME ================================
!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
!        WRITE(STRIP, '(1I3.3)') MYID

!        IF(MYID < NPTOT - 1) THEN
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_inouDOmAIn_INSTANT_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'C.bin'
!        ELSE
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_inouDOmAIn_INSTANT_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'E.bin'
!        END IF

!        !==============OPEN FILE =============================================
!        OPEN(DFLG, FILE = TRIM(FIL1), FORM = 'UNFORMATTED',STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF
!        CALL CHKHDL( '#    READING FILE FROM: ' // TRIM(FIL1), MYID)

!        !=============== READ IN PROCESSOR INFO.====================
!        READ(DFLG) MYID_tmp, NPTOT_tmp

!        IF(MYID /= MYID_tmp) THEN
!            WRITE(*, *) 'MYID /= iP, WHICH ARE', MYID, MYID_tmp, ' IN FILE ' // TRIM(FIL1)
!        END IF
!        IF(NPTOT_tmp /= NPTOT) THEN
!            WRITE(*, *) 'NPTOT_tmp /= NPTOT, WHICH ARE', NPTOT_tmp, NPTOT, ' IN FILE ' // TRIM(FIL1)
!        END IF

!        !===================== READ IN DATA============================

!        READ(10) N1ML, N2DOL, N3ML, ITERG0_io
!        READ(10) PhyTIME_io,RENL, DTL

!        READ(DFLG) G_io
!        READ(DFLG) PR_io

!        IF( iIniFieldType == 0 ) THEN

!            READ(DFLG) TEMPERATURE
!            READ(DFLG) DENSITY
!            READ(DFLG) ENTHALPY
!            READ(DFLG) DENSITY0

!            !=============BUILD UP DH ============
!            DO J = 1, N2DO(MYID)
!                DO I = 0, NCL1_io + 1, 1
!                    DO K = 1, NCL3
!                        DH(I, J, K) = DENSITY(I, J, K) * ENTHALPY(I, J, K)
!                    END DO
!                END DO
!            END DO

!            !===========build up other variables =======
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML
!            CALL THERM_PROP_UPDATE(IALLDM)
!            CALL INTFC_OUL_THERMAL_io
!            CALL INTFC_MFD_THERMAL_io
!            IF(iThermalWallType == BC_Fixed_Heat_Flux ) CALL BC_WALL_Isoflux(IALLDM)
!            CALL VELOCITY_CALC_io
!            !===========build up mass fluX =========
!            !CALL Massflux_reStart_io

!        END IF
!        CLOSE(DFLG)

!        IF(iIniFieldType == 1) THEN

!            CALL IniField_THERMAL_io
!            Q_io = G_io   ! based on constant unit DENSITY
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML

!            IF(iThermoDynamics == 1)THEN
!                CALL SOLVERRK3_ENG_io(0)
!                PhyTIME_io = 0.0_WP
!                ITERG0_io = 0
!            END IF

!        END IF

!        RETURN
!    END SUBROUTINE

!!*************************************************************************
!SUBROUTINE REStart_AVERAGE_GLOBAL_io
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE thermal_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: DFLG = 10
!        INTEGER(4) :: I, J, JJ
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML

!        INTEGER(4) :: IOS

!        REAL(WP), ALLOCATABLE :: U1ztL_F0_io(:, :, :)!Ui
!        REAL(WP), ALLOCATABLE :: G1ztL_F0_io(:, :, :)!Gi
!        REAL(WP), ALLOCATABLE :: UPztL_F0_io(:, :, :)!UiP
!        REAL(WP), ALLOCATABLE :: U2ztL_F0_io(:, :, :)!UiUj
!        REAL(WP), ALLOCATABLE :: UGztL_F0_io(:, :, :)!UiGj
!        REAL(WP), ALLOCATABLE :: UGUztL_F0_io(:, :, :)!UiGjUk

!        REAL(WP), ALLOCATABLE :: DVDL1ztL_F0_io(:, :, :, :) ! dUI / DXj
!        REAL(WP), ALLOCATABLE :: DVDLPztL_F0_io(:, :, :, :) !PdUI / DXj
!        REAL(WP), ALLOCATABLE :: DVDL2ztL_F0_io(:, :, :, :) ! dUI / DXj * dUI / DXj

!        REAL(WP), ALLOCATABLE :: T1ztL_F0_io(:, :)
!        REAL(WP), ALLOCATABLE :: D1ztL_F0_io(:, :)
!        REAL(WP), ALLOCATABLE :: H1ztL_F0_io(:, :)
!        REAL(WP), ALLOCATABLE :: T2ztL_F0_io(:, :)
!        REAL(WP), ALLOCATABLE :: D2ztL_F0_io(:, :)
!        REAL(WP), ALLOCATABLE :: H2ztL_F0_io(:, :)

!        REAL(WP), ALLOCATABLE :: DHztL_F0_io(:, :)

!        REAL(WP), ALLOCATABLE :: UHztL_F0_io(:, :, :)
!        REAL(WP), ALLOCATABLE :: GHztL_F0_io(:, :, :)

!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_io

!        ALLOCATE ( U1ztL_F0_io( NCL1_io, NCL2, NDV + 1 ) )
!        ALLOCATE ( G1ztL_F0_io( NCL1_io, NCL2, NDV   ) )
!        ALLOCATE ( UPztL_F0_io( NCL1_io, NCL2, NDV   ) )
!        ALLOCATE ( U2ztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3) )
!        ALLOCATE ( UGztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3) )
!        ALLOCATE ( UGUztL_F0_io(NCL1_io, NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) )

!        ALLOCATE ( DVDL1ztL_F0_io( NCL1_io, NCL2, NDV, NDV  ) )
!        ALLOCATE ( DVDLPztL_F0_io( NCL1_io, NCL2, NDV, NDV  ) )
!        ALLOCATE ( DVDL2ztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3, NDV  ) )

!        IF(iIniFieldType == 0)THEN
!            ALLOCATE ( T1ztL_F0_io( NCL1_io, NCL2 ) )
!            ALLOCATE ( D1ztL_F0_io( NCL1_io, NCL2 ) )
!            ALLOCATE ( H1ztL_F0_io( NCL1_io, NCL2 ) )

!            ALLOCATE ( T2ztL_F0_io( NCL1_io, NCL2 ) )
!            ALLOCATE ( D2ztL_F0_io( NCL1_io, NCL2 ) )
!            ALLOCATE ( H2ztL_F0_io( NCL1_io, NCL2 ) )

!            ALLOCATE ( UHztL_F0_io( NCL1_io, NCL2, NDV ) )
!            ALLOCATE ( GHztL_F0_io( NCL1_io, NCL2, NDV ) )
!        END IF
!        !=======================================U=============================================
!        FIL1 = 'DNS_inouDOmAIn_RST_' // TRIM(PNTIM) // '_FLOW_AVEG.bin'
!        IF(MYID == 0) CALL CHKHDL( '    Start reading ' //FIL1, MYID)

!        OPEN(DFLG, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF

!        READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_io
!        READ(DFLG) NSTATIS_io
!        READ(DFLG) PhyTIME_io,RENL, DTL
!        IF(MYID == 0) THEN
!            CALL CHKRLHDL   ('        FLOW_AVEG field TimE = ', MYID, PhyTIME_io)
!            CALL CHKINTHDL  ('        FLOW_AVEG field ITERG=   ', MYID, ITERG0_io)
!            CALL CHKINTHDL  ('        FLOW_AVEG field INSTATIS_io =   ', MYID, NSTATIS_io)
!        END IF

!        READ(DFLG) U1ztL_F0_io
!        READ(DFLG) G1ztL_F0_io
!        READ(DFLG) UPztL_F0_io

!        READ(DFLG) U2ztL_F0_io
!        READ(DFLG) UGztL_F0_io
!        READ(DFLG) UGUztL_F0_io

!        READ(DFLG) DVDL1ztL_F0_io
!        READ(DFLG) DVDLPztL_F0_io
!        READ(DFLG) DVDL2ztL_F0_io

!        CLOSE(DFLG)
!        IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(FIL1), MYID)

!        DO J = 1, N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I = 1, NCL1_io
!                U1ztL_io( I, J,: ) = U1ztL_F0_io( I, JJ,: )
!                G1ztL_io( I, J,: ) = G1ztL_F0_io( I, JJ,: )
!                UPztL_io( I, J,: ) = UPztL_F0_io( I, JJ,: )

!                U2ztL_io( I, J,: ) = U2ztL_F0_io( I, JJ,: )
!                UGztL_io( I, J,: ) = UGztL_F0_io( I, JJ,: )

!                UGUztL_io( I, J,: ) = UGUztL_F0_io( I, JJ,: )

!                DVDL1ztL_io( I, J, :,: ) = DVDL1ztL_F0_io( I, JJ, :,: )
!                DVDLPztL_io( I, J, :,: ) = DVDLPztL_F0_io( I, JJ, :,: )
!                DVDL2ztL_io( I, J, :,: ) = DVDL2ztL_F0_io( I, JJ, :,: )
!            END DO
!        END DO

!        IF(iIniFieldType == 0) THEN

!            FIL1 = 'DNS_inouDOmAIn_RST_' // TRIM(PNTIM) // '_THEL_AVEG.bin'
!            IF(MYID == 0) CALL CHKHDL( '    Start reading data from ' //FIL1, MYID)

!            OPEN(DFLG, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!            IF(IOS/= 0) THEN
!                CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!            END IF

!            READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_io
!            READ(DFLG) NSTATIS_io
!            READ(DFLG) PhyTIME_io,RENL, DTL
!            IF(MYID == 0) THEN
!                CALL CHKRLHDL   ('        THEL_AVEG field TimE = ', MYID, PhyTIME_io)
!                CALL CHKINTHDL  ('        THEL_AVEG field ITERG=   ', MYID, ITERG0_io)
!            END IF

!            READ(DFLG) T1ztL_F0_io
!            READ(DFLG) D1ztL_F0_io
!            READ(DFLG) H1ztL_F0_io

!            READ(DFLG) T2ztL_F0_io
!            READ(DFLG) D2ztL_F0_io
!            READ(DFLG) H2ztL_F0_io

!            READ(DFLG) DHztL_F0_io

!            READ(DFLG) UHztL_F0_io
!            READ(DFLG) GHztL_F0_io

!            CLOSE(DFLG)
!            IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(FIL1), MYID)

!            DO J = 1, N2DO(MYID)
!                JJ = JCL2G(J)
!                DO I = 1, NCL1_io
!                    T1ztL_io( I, J ) = T1ztL_F0_io( I, JJ )
!                    D1ztL_io( I, J ) = D1ztL_F0_io( I, JJ )
!                    H1ztL_io( I, J ) = H1ztL_F0_io( I, JJ )

!                    T2ztL_io( I, J ) = T2ztL_F0_io( I, JJ )
!                    D2ztL_io( I, J ) = D2ztL_F0_io( I, JJ )
!                    H2ztL_io( I, J ) = H2ztL_F0_io( I, JJ )

!                    UHztL_io( I, J,: ) = UHztL_F0_io( I, JJ,: )
!                    GHztL_io( I, J,: ) = GHztL_F0_io( I, JJ,: )
!                END DO
!            END DO

!        END IF


!        DEALLOCATE (U1ztL_F0_io, G1ztL_F0_io, UPztL_F0_io, U2ztL_F0_io, UGztL_F0_io, &
!        UGUztL_F0_io, DVDL1ztL_F0_io, DVDLPztL_F0_io, DVDL2ztL_F0_io)
!        IF( iIniFieldType == 0 ) DEALLOCATE (T1ztL_F0_io, D1ztL_F0_io, H1ztL_F0_io, &
!        T2ztL_F0_io, D2ztL_F0_io, H2ztL_F0_io, UHztL_F0_io, GHztL_F0_io)

!        CALL MPI_BARRIER(ICOMM, IERROR)

!        RETURN
!    END SUBROUTINE
!!********************************************************************************
!SUBROUTINE REStart_AVERAGE_LOCAL_io
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE thermal_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: DFLG = 10
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML
!        INTEGER(4) :: MYID_tmp, NPTOT_tmp

!        INTEGER(4) :: IOS

!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        !================CREAT FILE NAME ====================================
!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_io
!        WRITE(STRIP, '(1I3.3)') MYID
!        IF(MYID < NPTOT - 1) THEN
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_inouDOmAIn_AVERAGE_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'C.bin'
!        ELSE
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_inouDOmAIn_AVERAGE_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'E.bin'
!        END IF

!        !==============OPEN FILE =============================================
!        OPEN(DFLG, FILE = TRIM(FIL1), FORM = 'UNFORMATTED',STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF
!        CALL CHKHDL( '#    READING FILE FROM: ' // TRIM(FIL1), MYID)

!        !============== Read in processor info.================================
!        READ(DFLG) MYID_tmp, NPTOT_tmp

!        IF(MYID /= MYID_tmp) THEN
!            WRITE(*, *) 'MYID /= iP, WHICH ARE', MYID, MYID_tmp, ' IN FILE ' // TRIM(FIL1)
!        END IF
!        IF(NPTOT_tmp /= NPTOT) THEN
!            WRITE(*, *) 'NPTOT_tmp /= NPTOT, WHICH ARE', NPTOT_tmp, NPTOT, ' IN FILE ' // TRIM(FIL1)
!        END IF

!        !================= Read in data================================
!        READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_io
!        READ(DFLG) NSTATIS_io
!        READ(DFLG) PhyTIME_io,RENL, DTL

!        READ(DFLG) U1ztL_io
!        READ(DFLG) G1ztL_io
!        READ(DFLG) UPztL_io

!        READ(DFLG) U2ztL_io
!        READ(DFLG) UGztL_io

!        READ(DFLG) UGUztL_io

!        READ(DFLG) DVDL1ztL_io
!        READ(DFLG) DVDLPztL_io
!        READ(DFLG) DVDL2ztL_io

!        IF(iIniFieldType == 0) THEN

!            READ(DFLG) T1ztL_io
!            READ(DFLG) D1ztL_io
!            READ(DFLG) H1ztL_io

!            READ(DFLG) T2ztL_io
!            READ(DFLG) D2ztL_io
!            READ(DFLG) H2ztL_io

!            READ(DFLG) DHztL_io

!            READ(DFLG) UHztL_io
!            READ(DFLG) GHztL_io

!        END IF

!        CLOSE(DFLG)

!        RETURN
!    END SUBROUTINE
