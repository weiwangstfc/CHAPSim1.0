!**********************************************************************************************************************************
!> @brief
!>        reStart
!> @details
!> subroutine: ReStart_INSTANT_VARS_TG
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE ReStart_INSTANT_VARS_TG(TFM)
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    USE thermal_info
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: TFM
    CHARACTER(LEN = 256) :: WRT_RST_FNM
    CHARACTER(15) :: PNTIM
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
    INTEGER(4) :: I, J, K, N

    IF(MYID == 0) CALL CHKHDL('TG: ReStart instantanous flow field', MYID)

    ALLOCATE(DUMMY(NCL1_TG, N2DO(MYID), NCL3))
    WRITE(PNTIM, '(1ES15.9)') TFM !TimeReStart_tg
    DO N = 1, NDV + 1
        IF(N == 1)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_U.D'
        IF(N == 2)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_V.D'
        IF(N == 3)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_W.D'
        IF(N == 4)  WRITE(WRT_RST_FNM, '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_P.D'

        CALL READ_3D_VARS(NCL1_TG, NCL2, NCL3, N2DO(MYID), JCL2G(1) - 1, ITERG0_TG, PhyTIME_TG, DUMMY, WRT_RST_FNM)

        DO I = 1, NCL1_TG
            DO J = 1, N2DO(MYID)
                DO K = 1, NCL3
                    IF(N == 1) Q_tg(I, J, K, 1) = DUMMY(I, J, K) !U
                    IF(N == 2) Q_tg(I, J, K, 2) = DUMMY(I, J, K) !V
                    IF(N == 3) Q_tg(I, J, K, 3) = DUMMY(I, J, K) !W
                    IF(N == 4) PR_TG(I, J, K) = DUMMY(I, J, K) !P
                END DO
            END DO
        END DO

    END DO

    DEALLOCATE(DUMMY)

    IF(MYID == 0) CALL CHKRLHDL('===========PhyTIME_tg ================== ', MYID, PhyTIME_TG)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE ReStart_AVERAGE_VARS_TG
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG
    INTEGER(4) :: IERR, NEWTYPE
    INTEGER(4) :: INISIZE, IRLSIZE, INTBYTE, DBLBYTE
    INTEGER(4) :: NARRAY, NSZ
    INTEGER(4), ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
    INTEGER(4), ALLOCATABLE :: INTMPI(:)
    REAL(WP), ALLOCATABLE :: RLEMPI(:)
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, J, L, L1, L2, N, NLLL, NNNL, M, K
    REAL(WP) :: RENL, DTL
    REAL(WP), ALLOCATABLE :: DUMMY(:, :), DUMMY1(:, :), DUMMY2(:, :)
    LOGICAL :: File_exists
    INTEGER(4) :: NSTATIS_TG1, NSTATIS_TG2

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged flow field', MYID)

    !============================
    NARRAY = 2
    ALLOCATE ( SIZES_ARRAY(NARRAY) )
    ALLOCATE ( SUBSIZES   (NARRAY) )
    ALLOCATE ( STARTS     (NARRAY) )
    INISIZE = 4
    IRLSIZE = 3
    ALLOCATE ( INTMPI(INISIZE)       )
    ALLOCATE ( RLEMPI(IRLSIZE)       )

    NSZ = NDV + 1 + NDV + NDV * (7 - NDV) / 2 + NDV - 3 +  NDV * (6 - NDV) + &
    (NDV * (7 - NDV)) / 2 + NDV - 8  + NDV * NDV  + NDV * NDV +  &
    ((NDV - 1) * 3 + NDV) * ((NDV - 1) * 3 + NDV) !+ 4* (4* QUADHN)
    IF(iPPQuadrants == 1)  NSZ = NSZ + 4* (4* QUADHN)

    ALLOCATE( DUMMY (1 : N2DO(MYID), NSZ ) )
    ALLOCATE( DUMMY1(1 : N2DO(MYID), NSZ ) )
    DUMMY  = 0.0_WP
    DUMMY1 = 0.0_WP


    !=======prepARe sizes and rangeS ===================================================
    SIZES_ARRAY(1) = NCL2
    SIZES_ARRAY(2) = NSZ

    SUBSIZES(1) = N2DO(MYID)
    SUBSIZES(2) = SIZES_ARRAY(2)

    STARTS(1) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    STARTS(2) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2)

    !=========== Read data1 ==============================================================
    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg
    WRITE(WRT_AVE_FNM_tg, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

    IF(MYID == 0) CALL CHKHDL(' IO: ReStart averaged flow field ' // TRIM(WRT_AVE_FNM_tg), MYID)

    INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_tg)), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_tg // ' does not exist.', 0)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_tg, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

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
    ITERG0_TG  = INTMPI(3)
    NSTATIS_TG1 = INTMPI(4)

    PhyTIME_TG = RLEMPI(1)
    RENL     = RLEMPI(2)
    DTL      = RLEMPI(3)

    IF(MYID == 0) THEN
        CALL CHK2INTHDL('     Expected and read in DUMMY Sizes are', MYID, NSZ, NNNL)
        IF(NSZ /= NNNL) CALL ERRHDL('ReStarting averaged data formats are not consistent! ', MYID)
        CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, N2ML)
        CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, NNNL)
        CALL CHKINTHDL('     ITERG0_tg IN AVERAGED variables ', MYID, ITERG0_TG)
        CALL CHKINTHDL('     NSTATIS_tg IN AVERAGED variables', MYID, NSTATIS_TG1)

        CALL CHKRLHDL ('     PhyTIME_tg IN AVERAGED variables', MYID, PhyTIME_TG)
        CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RENL)
        CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, DTL)
    END IF

    NSTATIS_TG = NSTATIS_TG1
    DUMMY = DUMMY1

    !=========== Read data2 ==============================================================
    IF(tRunAve_Reset < PhyTIME_TG .AND. tRunAve_Reset > tRunAve1) THEN
        ALLOCATE( DUMMY2(1 : N2DO(MYID), NSZ ) )
        DUMMY2 = 0.0_WP

        DFLG = 101
        WRITE(PNTIM, '(1ES15.9)') tRunAve_Reset
        WRITE(WRT_AVE_FNM_TG, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

        IF(MYID == 0) CALL CHKHDL(' IO: Resetting averaged flow field ' // TRIM(WRT_AVE_FNM_tg), MYID)

        INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_tg)), exist =File_exists)
        IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_tg // ' does not exist.', 0)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY,SUBSIZES,STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_tg, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

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


        NSTATIS_tg2 = INTMPI(4)
        IF(MYID == 0) THEN
            IF(NSZ /= INTMPI(2)) CALL ERRHDL('ReStarting averaged data formats are not consistent! ', MYID)
            CALL CHK2INTHDL('     Expected and reaD -in DUMMY Sizes ARe', MYID, NSZ, INTMPI(2))
            CALL CHKINTHDL('     J IN AVERAGED variables (J, M) ', MYID, INTMPI(1))
            CALL CHKINTHDL('     M IN AVERAGED variables (J, M) ', MYID, INTMPI(2))
            CALL CHKINTHDL('     ITERG0_tg IN AVERAGED variables ', MYID, INTMPI(3))
            CALL CHKINTHDL('     NSTATIS_tg IN AVERAGED variables', MYID, NSTATIS_tg2)

            CALL CHKRLHDL ('     PhyTIME_tg IN AVERAGED variables', MYID, RLEMPI(1))
            CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RLEMPI(2))
            CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, RLEMPI(3))
        END IF

        DUMMY = (DUMMY1 * DBLE(NSTATIS_tg1) - DUMMY2 * DBLE(NSTATIS_tg2)) / DBLE(NSTATIS_tg1-NSTATIS_tg2)
        NSTATIS_tg = NSTATIS_tg1-NSTATIS_tg2
        DEALLOCATE (DUMMY2)
    END IF

    !==========================================
    N = 0
    L1 = NDV + 1
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U1xztL_tg(J, L) = DUMMY(J, N +L) !1-4
        END DO
    END DO
    N = N + L1 !4


    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            UPxztL_tg(J, L) = DUMMY(J, N +L) !8- 10
        END DO
    END DO
    N = N + L1 !10

    !==========================================
    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U2xztL_tg(J, L) = DUMMY(J, N +L) !11- 16
        END DO
    END DO
    N = N + L1 !16


    !        ! Commented below IS no U3
    L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8 !10
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            U3xztL_tg(J, L) = DUMMY(J, N +L)
        END DO
    END DO
    N = N + L1 !32

    !==========================================
    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDL1xztL_tg(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2


    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDLPxztL_tg(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2

    L1 = (NDV - 1) * 3 + NDV !9
    L2 = (NDV - 1) * 3 + NDV !9
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DVDL2xztL_tg(J, L, K) = DUMMY(J, M)
            END DO
        END DO
    END DO
    N = N + L1 * L2

    !============= QuadranT ===============
    IF(iPPQuadrants == 1)  THEN
        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADUVxztL_tg(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADVzxztL_tg(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADTKxztL_tg(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2

        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    QUADDRxztL_tg(J, L, K) = DUMMY(J, M)
                END DO
            END DO
        END DO
        N = N + L1 * L2
    END IF

    IF(N /= NSZ) THEN
        CALL ERRHDL('#ERROR in ReStart_AVERAGE_VARS_Xperiodic_tg', MYID)
    END IF

    DEALLOCATE (DUMMY)
    DEALLOCATE (DUMMY1)
    IF(MYID == 0) CALL CHKRLHDL('    IO: End of reStarting averaged flow    field', MYID, PhyTIME_tg)

    RETURN

END SUBROUTINE



!SUBROUTINE ReStart_AVERAGE_VARS_TG
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        USE WRT_INFO
!        USE POSTPROCESS_INFO
!        USE thermal_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        INTEGER(4) ::  DFLG, NARRAY, L1, N, L2, J, L, M, NLLL, NNNL, N2ML, H, P
!        INTEGER(4) :: IERR,SIZES_ARRAY(3), NEWTYPE,SUBSIZES(3),STARTS(3)
!        INTEGER(4) :: INISIZE, IRLSIZE, INTMPI(4), INTBYTE, DBLBYTE
!        INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
!        INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
!        REAL(WP) :: RLEMPI(3)
!        REAL(WP) :: RENL, DTL
!        REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)
!        LOGICAL :: File_exists

!        IF(MYID == 0) CALL CHKHDL(' TG: ReStart averaged flow field', MYID)

!        ALLOCATE(DUMMY(1 : N2DO(MYID), 1:(NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8), 1:19) )
!        DUMMY = 0.0_WP

!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg


!      ==================GENERATE FILE NAMES ===============================================
!        DFLG = 100
!        WRITE(WRT_AVE_FNM_TG, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

!        INQUIRE(FILE = TRIM(ADJUSTL(WRT_AVE_FNM_TG)), exist =File_exists)
!        IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //WRT_AVE_FNM_tg // ' does not exist.', 0)

!        SIZES_ARRAY(1) = NCL2
!        SIZES_ARRAY(2) = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8
!        SIZES_ARRAY(3) = 19

!        SUBSIZES(1) = N2DO(MYID)
!        SUBSIZES(2) = SIZES_ARRAY(2)
!        SUBSIZES(3) = SIZES_ARRAY(3)

!        STARTS(1) = JCL2G(1) - 1
!        STARTS(2) = 0
!        STARTS(3) = 0

!        BUFSIZE  = SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)

!        CALL MPI_TYPE_CREATE_SUBARRAY(3, SIZES_ARRAY, SUBSIZES, STARTS, &
!               MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
!        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
!        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_TG, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

!        INISIZE = 4
!        IRLSIZE = 3

!        OFFSET = 0_MPI_OFFSET_KIND
!        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
!        CALL MPI_FILE_READ(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
!        CALL MPI_FILE_READ(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
!        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
!        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

!        OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE
!        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
!        CALL MPI_FILE_READ_ALL(MYFILE, DUMMY,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
!        CALL MPI_FILE_CLOSE(MYFILE, IERR)
!        CALL MPI_TYPE_FREE(NEWTYPE, IERR)

!        N2ML    = INTMPI(1)
!        NNNL    = INTMPI(2)
!        ITERG0_TG = INTMPI(3)
!        NSTATIS_tg = INTMPI(4)

!        PhyTIME_tg = RLEMPI(1)
!        RENL    = RLEMPI(2)
!        DTL     = RLEMPI(3)

!        IF(MYID == 0) THEN
!            CALL CHKINTHDL('     J IN AVERAGED variables (J, L, M) ', MYID, N2ML)
!            CALL CHKINTHDL('     L IN AVERAGED variables (J, L, M) ', MYID, NNNL)
!            CALL CHKINTHDL('     ITERG0_TG IN AVERAGED variables ', MYID, ITERG0_TG)
!            CALL CHKINTHDL('     NSTATIS_tg IN AVERAGED variables', MYID, NSTATIS_tg)

!            CALL CHKRLHDL ('     PhyTIME_TG IN AVERAGED variables', MYID, PhyTIME_TG)
!            CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RENL)
!            CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, DTL)
!        END IF


!        DO J = 1, N2DO(MYID)
!            DO L = 1, NDV + 1
!                U1xztL_tg(J, L) = DUMMY(J, L, 1)
!            END DO
!        END DO

!        DO J = 1, N2DO(MYID)
!            DO L = 1, NDV
!                UPxztL_tg(J, L) = DUMMY(J, L, 2)
!            END DO
!        END DO

!        L1 = NDV * (7 - NDV) / 2 + NDV - 3
!        DO J = 1, N2DO(MYID)
!            DO L = 1, L1
!                U2xztL_tg(J, L) = DUMMY(J, L, 3)
!            END DO
!        END DO

!        L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8
!        DO J = 1, N2DO(MYID)
!            DO L = 1, L1
!                U3xztL_tg(J, L) = DUMMY(J, L, 4)
!            END DO
!        END DO

!        DO M = 1, NDV
!            DO J = 1, N2DO(MYID)
!                DO L = 1, NDV
!                    DVDL1xztL_tg(J, L, M) = DUMMY(J, L, M + 4)
!                END DO
!            END DO
!        END DO

!        DO M = 1, NDV
!            DO J = 1, N2DO(MYID)
!            DO L = 1, NDV
!                    DVDLPxztL_tg(J, L, M) = DUMMY(J, L, M +7)
!                END DO
!            END DO
!        END DO

!        DO M = 1, NDV
!            DO H = 1, NDV
!                DO N = 1, NDV
!                    DO P = 1, NDV
!                        L1 = (M - 1) * 3 + H
!                        L2 = (N - 1) * 3 + P
!                        DO J = 1, N2DO(MYID)
!                            DVDL2xzL_tg(J, L1, L2) = DUMMY(J, L1, L2 + 10)
!                        END DO
!                    END DO
!                END DO
!            END DO
!        END DO

!        !test
!        DO J = 1, N2DO(MYID)
!            m = JCl2g(J)
!            WRITE(*, '(1i4.1, 4es13.5)') m, U1xztL_tg(J, 1:4)
!        END DO

!        DO J = 1, N2DO(MYID)
!            m = JCl2g(J)
!            DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
!                WRITE(*, '(2i4.1, 13es13.5)') m, l, DUMMY(J, L, 1:(4 + NDV *3))
!            END DO
!        END DO


!        DEALLOCATE (DUMMY)

!        IF(MYID == 0) CALL CHKRLHDL('    TG: End of ReStarting averaged flow field', MYID, PhyTIME_TG)


!        RETURN
!    END SUBROUTINE


!!**********************************************************************************************************
!SUBROUTINE REStart_INSTANT_GLOBAL_tg
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: IOS
!        INTEGER(4) :: I, J, K, JJ
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML

!        REAL(WP), ALLOCATABLE :: U_F0_tg(:, :, :)
!        REAL(WP), ALLOCATABLE :: V_F0_tg(:, :, :)
!        REAL(WP), ALLOCATABLE :: W_F0_tg(:, :, :)
!        REAL(WP), ALLOCATABLE :: P_F0_tg(:, :, :)
!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg

!        ALLOCATE( U_F0_tg(1 : NCL1_tg, 1 : NCL2,  1 : NCL3) )
!        ALLOCATE( W_F0_tg(1 : NCL1_tg, 1 : NCL2,  1 : NCL3) )
!        ALLOCATE( V_F0_tg(1 : NCL1_tg, 1 : NCL2 + 1, 1 : NCL3) )
!        ALLOCATE( P_F0_tg(1 : NCL1_tg, 1 : NCL2,  1 : NCL3) )

!        !=======================================UVW=============================================
!        FIL1 = 'DNS_periodicxz_RST_' // TRIM(PNTIM) // '_FLOW_INST.bin'

!        IF(MYID == 0) CALL CHKHDL( '    Start reading ' //FIL1, MYID)

!        OPEN(10, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF

!        READ(10) N1ML, N2DOL, N3ML, ITERG0_TG
!        READ(10) PhyTIME_Tg, rENL, DTL

!        IF(MYID == 0) THEN
!            CALL CHKRLHDL   ('        UVWP field TimE = ', MYID, PhyTIME_TG)
!            CALL CHKINTHDL  ('        UVWP field ITERG=   ', MYID, ITERG0_TG)
!        END IF

!        READ(10) U_F0_tg
!        READ(10) V_F0_tg
!        READ(10) W_F0_tg
!        READ(10) P_F0_tg
!        CLOSE(10)

!        IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' // TRIM(FIL1), MYID)
!        IF(DABS(PhyTIME_TG- TimeReStart_tg) >TSAVE1) CALL CHKHDL('WARning: ReStarting from an ealier time step!', MYID)

!        !=======================================================================================
!        DO J = 1, N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I = 1, NCL1_tg
!                DO K = 1, NCL3
!                    Q_tg(I, J, K, 1) = (U_F0_tg(I, JJ, K))
!                    Q_tg(I, J, K, 2) = (V_F0_tg(I, JJ, K))
!                    Q_tg(I, J, K, 3) = (W_F0_tg(I, JJ, K))
!                    PR_tg(I, J, K) = (P_F0_tg(I, JJ, K))
!                END DO
!            END DO
!        END DO

!        IF (MYID == NPSLV) THEN
!            J = N2DO(MYID) + 1
!            DO I = 1, NCL1_tg
!                DO K = 1, NCL3
!                    Q_tg(I, J, K, 2) = (V_F0_tg(I, NND2, K))
!                END DO
!            END DO
!        ENDIF

!        IF (MYID == 0) THEN
!            DO I = 1, NCL1_tg
!                DO K = 1, NCL3
!                    Q_tg(I, 1, K, 2) = 0.0_WP
!                    Q_tg(I, 0, K, 2) = 0.0_WP
!                END DO
!            END DO
!        ENDIF

!        DEALLOCATE (U_F0_tg,W_F0_tg,P_F0_tg, V_F0_tg)
!        CALL MPI_BARRIER(ICOMM, IERROR)

!        RETURN
!    END SUBROUTINE

!!**********************************************************************************************************
!SUBROUTINE REStart_INSTANT_LOCAL_tg
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: IOS
!        INTEGER(4) :: DFLG = 10
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML
!        INTEGER(4) :: MYID_tmp, NPTOT_tmp

!        REAL(WP) :: RENL
!        REAL(WP) :: DTL

!        !=============creat file namE ========================================
!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg
!        WRITE(STRIP, '(1I3.3)')   MYID

!        IF(MYID < NPTOT - 1) THEN
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_periodicxz_INSTANT_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'C.bin'
!        ELSE
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_periodicxz_INSTANT_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'E.bin'
!        END IF

!        !==============OPEN FILE =============================================
!        OPEN(DFLG, FILE = TRIM(FIL1), FORM = 'UNFORMATTED',STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF
!        CALL CHKHDL( '#    READING FILE FROM: ' // TRIM(FIL1), MYID)

!        !================= Read in processor info.================================
!        READ(DFLG) MYID_tmp, NPTOT_tmp

!        IF(MYID /= MYID_tmp) THEN
!            WRITE(*, *) 'MYID /= iP, WHICH ARE', MYID, MYID_tmp, ' IN FILE ' // TRIM(FIL1)
!        END IF
!        IF(NPTOT_tmp /= NPTOT) THEN
!            WRITE(*, *) 'NPTOT_tmp /= NPTOT, WHICH ARE', NPTOT_tmp, NPTOT, ' IN FILE ' // TRIM(FIL1)
!        END IF

!        !================= Read in data================================
!        READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_TG
!        READ(DFLG) PhyTIME_Tg, rENL, DTL

!        READ(DFLG) Q_tg
!        READ(DFLG) PR_TG

!        CLOSE(DFLG)

!        RETURN
!    END SUBROUTINE

!***********************************************************************************
!SUBROUTINE REStart_AVERAGE_GLOBAL_tg
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: IOS
!        INTEGER(4) :: J, JJ
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML

!        REAL(WP), ALLOCATABLE :: U1xztL_F0_tg(:, :)
!        REAL(WP), ALLOCATABLE :: UPxztL_F0_tg(:, :)
!        REAL(WP), ALLOCATABLE :: U2xztL_F0_tg(:, :)
!        REAL(WP), ALLOCATABLE :: U3xztL_F0_tg(:, :)

!        REAL(WP), ALLOCATABLE :: DVDL1xztL_F0_tg(:, :, :)
!        REAL(WP), ALLOCATABLE :: DVDLPxztL_F0_tg(:, :, :)
!        REAL(WP), ALLOCATABLE :: DVDL2xztL_F0_tg(:, :, :)
!        REAL(WP) :: RENL
!        REAL(WP) :: DTL


!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg

!        ALLOCATE ( U1xztL_F0_tg( NCL2, NDV + 1 ) )
!        ALLOCATE ( UPxztL_F0_tg( NCL2, NDV   ) )
!        ALLOCATE ( U2xztL_F0_tg( NCL2, NDV * (7 - NDV) / 2 + NDV - 3                ) )
!        ALLOCATE ( U3xztL_F0_tg( NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8  ) )

!        ALLOCATE ( DVDL1xztL_F0_tg( NCL2, NDV, NDV  ) )
!        ALLOCATE ( DVDLPxztL_F0_tg( NCL2, NDV, NDV  ) )
!        ALLOCATE ( DVDL2xztL_F0_tg( NCL2, NDV * (7 - NDV) / 2 + NDV - 3, NDV  ) )

!        !=======================================UVW=============================================
!        FIL1 = 'DNS_periodicxz_RST_' // TRIM(PNTIM) // '_FLOW_AVEG.bin'

!        IF(MYID == 0) CALL CHKHDL( '    Start reading ' //FIL1, MYID)

!        OPEN(10, FILE =FIL1, FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOS)
!        IF(IOS/= 0) THEN
!            CALL ERRHDL('# cannot OPEN FILE: ' // TRIM(FIL1), MYID)
!        END IF

!        READ(10) N1ML, N2DOL, N3ML, ITERG0_TG
!        READ(10) NSTATIS_tg
!        READ(10) PhyTIME_Tg, rENL, DTL

!        IF(MYID == 0) THEN
!            CALL CHKRLHDL   ('        UVWP field TimE = ', MYID, PhyTIME_TG)
!            CALL CHKINTHDL  ('        UVWP field ITERG=   ', MYID, ITERG0_TG)
!        END IF

!        READ(10) U1xztL_F0_tg
!        READ(10) UPxztL_F0_tg
!        READ(10) U2xztL_F0_tg
!        READ(10) U3xztL_F0_tg

!        READ(10) DVDL1xztL_F0_tg
!        READ(10) DVDLPxztL_F0_tg
!        READ(10) DVDL2xztL_F0_tg

!        CLOSE(10)

!       IF(MYID == 0) CALL CHKHDL( '    FinISh reading ' //FIL1, MYID)
!        IF(DABS(PhyTIME_TG- TimeReStart_tg) >TSAVE1) CALL CHKHDL('WARning: ReStarting from an ealier time step!', MYID)
!        !=======================================================================================
!        DO J = 1, N2DO(MYID)
!            JJ = JCL2G(J)

!            U1xztL_tg(J, :) = U1xztL_F0_tg(JJ, :)
!            UPxztL_tg(J, :) = UPxztL_F0_tg(JJ, :)
!            U2xztL_tg(J, :) = U2xztL_F0_tg(JJ, :)
!            U3xztL_tg(J, :) = U3xztL_F0_tg(JJ, :)

!            DVDL1xztL_tg(J, :, :) = DVDL1xztL_F0_tg(JJ, :, :)
!            DVDLPxztL_tg(J, :, :) = DVDLPxztL_F0_tg(JJ, :, :)
!            DVDL2xztL_tg(J, :, :) = DVDL2xztL_F0_tg(JJ, :, :)

!        END DO

!        DEALLOCATE (U1xztL_F0_tg, UPxztL_F0_tg, U2xztL_F0_tg, U3xztL_F0_tg, &
!                     DVDL1xztL_F0_tg, DVDLPxztL_F0_tg, DVDL2xztL_F0_tg)

!        CALL MPI_BARRIER(ICOMM, IERROR)

!        RETURN
!    END SUBROUTINE


!!***********************************************************************************
!SUBROUTINE REStart_AVERAGE_LOCAL_tg
!        USE init_info
!        USE mesh_info
!        USE flow_info
!        USE postprocess_info
!        IMPLICIT NONE

!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1

!        INTEGER(4) :: DFLG = 10

!        INTEGER(4) :: IOS
!        INTEGER(4) :: N1ML
!        INTEGER(4) :: N2DOL
!        INTEGER(4) :: N3ML
!        INTEGER(4) :: MYID_tmp, NPTOT_tmp
!        REAL(WP) :: RENL
!        REAL(WP) :: DTL

!        !==================create file namE =================================
!        WRITE(PNTIM, '(1ES15.9)') TimeReStart_tg
!        WRITE(STRIP, '(1I3.3)') MYID

!        IF(MYID < NPTOT - 1) THEN
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_periodicxz_AVERAGE_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'C.bin'
!        ELSE
!            WRITE(FIL1, '(A, I4.4, A4)') &
!                'DNS_periodicxz_AVERAGE_RST_' // TRIM(PNTIM) // '_IP' // TRIM(STRIP) // 'E.bin'
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
!        READ(DFLG) N1ML, N2DOL, N3ML, ITERG0_TG
!        READ(DFLG) NSTATIS_tg
!        READ(DFLG) PhyTIME_Tg, rENL, DTL

!        READ(DFLG) U1xztL_tg
!        READ(DFLG) UPxztL_tg

!        READ(DFLG) U2xztL_tg
!        READ(DFLG) U3xztL_tg

!        READ(DFLG) DVDL1xztL_tg
!        READ(DFLG) DVDLPxztL_tg
!        READ(DFLG) DVDL2xztL_tg

!        CLOSE(DFLG)

!        RETURN
!    END SUBROUTINE
