!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> subroutine: WRT_INSTANT_VARS_TG
!> subroutine: WRT_AVERAGE_VARS_TG
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE WRT_INSTANT_VARS_TG
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: DFLG(4), IVEL, I, J, K
    INTEGER(4) :: IERR,SIZES_ARRAY(3)
    INTEGER(4) :: NEWTYPE,SUBSIZES(3),STARTS(3)
    INTEGER(4) :: IRLSIZE, INTMPI(4), INTBYTE, DBLBYTE, INISIZE
    INTEGER(4) ::  BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND)::  OFFSET
    REAL(WP) :: RLEMPI(3)
    REAL(WP), ALLOCATABLE :: DUMMY(:, :, :)

    ALLOCATE(DUMMY(1 : NCL1_TG, 1 : N2DO(MYID), 1 : NCL3))
    DUMMY = 0.0_WP

    DO IVEL = 1, NDV + 1
        !==================GENERATE FILE NAMES ===============================================
        DFLG(IVEL) = 100 + IVEL
        WRITE(PNTIM, '(1ES15.9)') PhyTIME
        IF(IVEL == 1) WRITE(WRT_RST_FNM_TG(IVEL), '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_U.D'
        IF(IVEL == 2) WRITE(WRT_RST_FNM_TG(IVEL), '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_V.D'
        IF(IVEL == 3) WRITE(WRT_RST_FNM_TG(IVEL), '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_W.D'
        IF(IVEL == 4) WRITE(WRT_RST_FNM_TG(IVEL), '(A)') TRIM(FilePath1) // 'DNS_perixz_INSTANT_T' // TRIM(PNTIM) // '_P.D'

        !======================== RE - ASSIGN variables =======================================
        DO I = 1, NCL1_tg
            DO J = 1, N2DO(MYID)
                DO K = 1, NCL3
                    IF(IVEL == 1) DUMMY(I, J, K) = Q_tg(I, J, K, 1) !U
                    IF(IVEL == 2) DUMMY(I, J, K) = Q_tg(I, J, K, 2) !V
                    IF(IVEL == 3) DUMMY(I, J, K) = Q_tg(I, J, K, 3) !W
                    IF(IVEL == 4) DUMMY(I, J, K) = PR_TG(I, J, K)  !P
                END DO
            END DO
        END DO

        SIZES_ARRAY(1) = NCL1_tg
        SIZES_ARRAY(2) = NCL2
        SIZES_ARRAY(3) = NCL3

        SUBSIZES(1) = NCL1_tg
        SUBSIZES(2) = N2DO(MYID)
        SUBSIZES(3) = NCL3

        STARTS(1) = 0
        STARTS(2) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
        STARTS(3) = 0
        BUFSIZE  = SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)

        CALL MPI_TYPE_CREATE_SUBARRAY(3,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_RST_FNM_TG(IVEL), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, MYFILE, IERR)

        INISIZE = 4
        IRLSIZE = 3

        INTMPI(1) = NCL1_tg
        INTMPI(2) = NCL2
        INTMPI(3) = NCL3
        INTMPI(4) = ITERG

        RLEMPI(1) = PhyTIME
        RLEMPI(2) = REN
        RLEMPI(3) = DT


        OFFSET = 0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_WRITE(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_WRITE(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

        OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE

        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
        CALL MPI_FILE_WRITE_ALL(MYFILE, DUMMY,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
        CALL MPI_FILE_CLOSE(MYFILE, IERR)
        CALL MPI_TYPE_FREE(NEWTYPE, IERR)
        CALL MPI_BARRIER(ICOMM, IERROR)
    END DO
    DEALLOCATE (DUMMY)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_VARS_TG
    USE MESH_INFO
    USE FLOW_INFO
    USE INIT_INFO
    USE WRT_INFO
    USE POSTPROCESS_INFO
    IMPLICIT NONE

    CHARACTER(15) :: PNTIM
    INTEGER(4) :: L1, L2, NSZ
    INTEGER(4) :: DFLG, J, L, M, N, K
    INTEGER(4) :: IERR, NEWTYPE
    INTEGER(4) :: NARRAY
    INTEGER(4), ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
    INTEGER(4), ALLOCATABLE :: INTMPI(:)
    REAL(WP), ALLOCATABLE :: RLEMPI(:)
    INTEGER(4) :: IRLSIZE, INTBYTE, DBLBYTE, INISIZE
    INTEGER(4) ::  BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND)::  OFFSET
    REAL(WP), ALLOCATABLE :: DUMMY(:, :)


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
    IF(iPPQuadrants == 1) NSZ =  NSZ + 4* (4* QUADHN)

    ALLOCATE( DUMMY(1 : N2DO(MYID), NSZ ) )
    DUMMY = 0.0_WP
    !============================

    SIZES_ARRAY(1) = NCL2
    SIZES_ARRAY(2) = NSZ

    SUBSIZES(1) = N2DO(MYID)
    SUBSIZES(2) = SIZES_ARRAY(2)

    STARTS(1) = JCL2G(1) - 1  !SUBSIZES(2) * MYID
    STARTS(2) = 0
    BUFSIZE  = SUBSIZES(1) * SUBSIZES(2)


    INTMPI(1) = NCL2
    INTMPI(2) = NSZ
    INTMPI(3) = ITERG
    INTMPI(4) = NSTATIS_tg

    RLEMPI(1) = PhyTIME
    RLEMPI(2) = REN
    RLEMPI(3) = DT

    !============================================================
    N = 0
    L1 = NDV + 1
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            DUMMY(J, N +L) = U1xztL_tg(J, L) !1-4
        END DO
    END DO
    N = N + L1 !4

    L1 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            DUMMY(J, N +L) = UPxztL_tg(J, L) !5-7
        END DO
    END DO
    N = N + L1 !10

    L1 = NDV * (7 - NDV) / 2 + NDV - 3
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            DUMMY(J, N +L) = U2xztL_tg(J, L)
        END DO
    END DO
    N = N + L1

    L1 = NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8
    DO J = 1, N2DO(MYID)
        DO L = 1, L1
            DUMMY(J, N +L) = U3xztL_tg(J, L)
        END DO
    END DO
    N = N + L1


    !==========================================
    L1 = NDV
    L2 = NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DUMMY(J, M) = DVDL1xztL_tg(J, L, K)
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
                DUMMY(J, M) = DVDLPxztL_tg(J, L, K)
            END DO
        END DO
    END DO
    N = N + L1 * L2

    L1 = (NDV - 1) * 3 + NDV
    L2 = (NDV - 1) * 3 + NDV
    DO J = 1, N2DO(MYID)
        DO L = 1, L1 !
            DO K = 1, L2
                M = N + (L- 1) * L2 + K
                DUMMY(J, M) = DVDL2xztL_tg(J, L, K)
            END DO
        END DO
    END DO
    N = N + L1 * L2

    !=======for quadranT ================
    IF(iPPQuadrants == 1)  THEN
        L1 = 4
        L2 = QUADHN
        DO J = 1, N2DO(MYID)
            DO L = 1, L1 !
                DO K = 1, L2
                    M = N + (L- 1) * L2 + K
                    DUMMY(J, M) = QUADUVxztL_tg(J, L, K)
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
                    DUMMY(J, M) = QUADVzxztL_tg(J, L, K)
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
                    DUMMY(J, M) = QUADTKxztL_tg(J, L, K)
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
                    DUMMY(J, M) = QUADDRxztL_tg(J, L, K)
                END DO
            END DO
        END DO
        N = N + L1 * L2
    END IF

    IF(N /= NSZ) THEN
        CALL ERRHDL('ERROR in WRT_AVERAGE_VARS_Xperiodic_TG', MYID)
    END IF



    !!==========================================

    DFLG = 100
    WRITE(PNTIM, '(1ES15.9)') PhyTIME
    WRITE(WRT_AVE_FNM_TG, '(A)') TRIM(FilePath2) // 'DNS_perixz_AVERAGD_T' // TRIM(PNTIM) // '_FLOW.D'

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_tg, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, MYFILE, IERR)

    !        IF(MYID == 0) THEN
    !            CALL CHKINTHDL('     J IN AVERAGED variables (J, L, M) ', MYID, INTMPI(1))
    !            CALL CHKINTHDL('     L IN AVERAGED variables (J, L, M) ', MYID, INTMPI(2))
    !            CALL CHKINTHDL('     M IN AVERAGED variables (J, L, M) ', MYID, INTMPI(3))
    !            CALL CHKINTHDL('     ITERG0_TG IN AVERAGED variables ', MYID, INTMPI(4))
    !            CALL CHKINTHDL('     NSTATIS_tg IN AVERAGED variables', MYID, INTMPI(5))

    !            CALL CHKRLHDL ('     PhyTIME_TG IN AVERAGED variables', MYID, RLEMPI(1))
    !            CALL CHKRLHDL ('     RENL       IN AVERAGED variables', MYID, RLEMPI(2))
    !            CALL CHKRLHDL ('     DTL        IN AVERAGED variables', MYID, RLEMPI(3))
    !        END IF

    OFFSET = 0_MPI_OFFSET_KIND
    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_WRITE(MYFILE, INTMPI, INISIZE, MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_WRITE(MYFILE,RLEMPI, IRLSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE, IERR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE, IERR)

    OFFSET = 0_MPI_OFFSET_KIND + INISIZE * INTBYTE + IRLSIZE * DBLBYTE

    CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET, MPI_DOUBLE_PRECISION, NEWTYPE,"NATIVE", MPI_INFO_NULL, IERR)
    CALL MPI_FILE_WRITE_ALL(MYFILE, DUMMY,BUFSIZE, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, IERR)
    CALL MPI_FILE_CLOSE(MYFILE, IERR)
    CALL MPI_TYPE_FREE(NEWTYPE, IERR)
    CALL MPI_BARRIER(ICOMM, IERROR)


    DEALLOCATE (DUMMY)

    DEALLOCATE ( SIZES_ARRAY )
    DEALLOCATE ( SUBSIZES    )
    DEALLOCATE ( STARTS      )
    DEALLOCATE ( INTMPI      )
    DEALLOCATE ( RLEMPI      )

    RETURN
END SUBROUTINE
