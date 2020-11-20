!**********************************************************************************************************************************
!> @brief
!>        read data for rStart/initialization
!> @details
!> subroutine: READ_3D_VARS
!> subroutine: INTP_3D_VARS
!> subroutine: INTRPRI
!> subroutine: INTRPRK
!> subroutine: INTRPRJ
!> subroutine: INTRPRJG
!> subroutine: SPLINE
!> subroutine: SPLINT
!> subroutine: TRASPO23_Y2Z
!> subroutine: TRASPO23_Z2Y
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE READ_3D_VARS(N1, N2, N3, N2L, N2P, IterIn, TimeIn, DUMMY, FLNM)
    USE WPRECISION
    USE mpi_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N1, N2, N3, N2L, N2P
    CHARACTER(LEN = 256), INTENT(IN) :: FLNM
    INTEGER(4), INTENT(OUT) :: IterIn
    REAL(WP), INTENT(OUT) :: TimeIn
    REAL(WP), INTENT(OUT) :: DUMMY(N1, N2L, N3)

    INTEGER(4) :: NARRAY
    INTEGER(4) :: IERR, SIZES_ARRAY(3), NEWTYPE,SUBSIZES(3),STARTS(3)
    INTEGER(4) :: INISIZE, IRLSIZE, INTMPI(4), INTBYTE, DBLBYTE
    INTEGER(4) :: BUFSIZE, MYFILE,STATUS(MPI_STATUS_SIZE)
    INTEGER(KIND = MPI_OFFSET_KIND) OFFSET
    INTEGER(4) :: N1ML, N2ML, N3ML, I, J, K, IVEL
    REAL(WP) :: RLEMPI(3)
    REAL(WP) :: RENL, DTL
    LOGICAL :: File_exists

    IF(MYID == 0) CALL CHKHDL('     Reading instantanous flow field in ' // TRIM(FLNM), MYID)
    NARRAY = 3

    INQUIRE( FILE = TRIM(ADJUSTL(FLNM) ), exist =File_exists)
    IF(.not.File_exists .AND. MYID == 0) CALL ERRHDL('File ' //FLNM // ' does not exist.', 0)

    SIZES_ARRAY(1) = N1
    SIZES_ARRAY(2) = N2
    SIZES_ARRAY(3) = N3

    SUBSIZES(1) = SIZES_ARRAY(1)
    SUBSIZES(2) = N2L
    SUBSIZES(3) = SIZES_ARRAY(3)

    STARTS(1) = 0
    STARTS(2) = N2P
    STARTS(3) = 0

    BUFSIZE = SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)

    CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
    CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)
    CALL MPI_FILE_OPEN(ICOMM, FLNM, MPI_MODE_RDONLY, MPI_INFO_NULL, MYFILE, IERR)

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
    N3ML   = INTMPI(3)
    IterIn = INTMPI(4)

    TimeIn = RLEMPI(1)
    RENL   = RLEMPI(2)
    DTL    = RLEMPI(3)

    IF(MYID == 0) THEN
        CALL CHKHDL   ('     variables IN   ' // Trim(FLNM), MYID)
        CALL CHKINTHDL('         I IN AVERAGED variables (I, J, K) ', MYID, N1ML)
        CALL CHKINTHDL('         J IN AVERAGED variables (I, J, K) ', MYID, N2ML)
        CALL CHKINTHDL('         K IN AVERAGED variables (I, J, K) ', MYID, N3ML)
        CALL CHKINTHDL('         IterIn    IN AVERAGED variables ', MYID, IterIn)

        CALL CHKRLHDL ('         TimeIn IN AVERAGED variables     ', MYID, TimeIn)
        CALL CHKRLHDL ('         RENL       IN AVERAGED variables ', MYID, RENL)
        CALL CHKRLHDL ('         DTL        IN AVERAGED variables ', MYID, DTL)
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INTP_3D_VARS(QO, QN, NDIR)
    USE INTP_VARS
    INTEGER(4), INTENT(IN) :: NDIR
    REAL(WP), INTENT(IN) :: QO(NCLO1, N2DOO(0), NCLO3)
    REAL(WP), INTENT(OUT) :: QN(NCL1, N2DO(0), NCL3 )

    REAL(WP) :: QIOKO( NCLO1, NCLO3 )
    REAL(WP) :: QINKO( NCL1, NCLO3 )
    REAL(WP) :: QINKN( NCL1, NCL3  )
    REAL(WP) :: QINJO( NCL1, N2DOO(0) )
    REAL(WP) :: QINJN( NCL1, N2DO(0)  )
    REAL(WP) :: QINJOKN( NCL1, N2DOO(0), NCL3 )
    REAL(WP) :: QINJNKN( NCL1, N2DO(0), NCL3 )

    INTEGER(4) :: I, J, K

    DO J = 1, N2DOO(MYID)

        DO I = 1, NCLO1
            DO K = 1, NCLO3
                QIOKO(I, K) = QO(I, J, K)
            END DO
        END DO

        IF(NDIR == 1) THEN
            CALL INTRPRI(QIOKO, QINKO, XNDO, XND, NCLO1, NCL1, NCLO3)
        ELSE
            CALL INTRPRI(QIOKO, QINKO, XCCO, XCC, NCLO1, NCL1, NCLO3)
        END IF

        IF(NDIR == 3) THEN
            CALL INTRPRK(QINKO, QINKN, ZNDO, ZND, NCLO3, NCL3, NCL1)
        ELSE
            CALL INTRPRK(QINKO, QINKN, ZCCO, ZCC, NCLO3, NCL3, NCL1)
        END IF

        DO I = 1, NCL1
            DO K = 1, NCL3
                QINJOKN(I, J, K) = QINKN(I, K)
            END DO
        END DO

    END DO

    IF(NDIR == 2) THEN
        CALL INTRPRJG( QINJOKN, QINJNKN, YNDO, YND, N2DOO, N2DO, NCLO2, NCL2, NCL1, NCL3)
    ELSE
        CALL INTRPRJG( QINJOKN, QINJNKN, YCCO, YCC, N2DOO, N2DO, NCLO2, NCL2, NCL1, NCL3)
    END IF

    DO J = 1, N2DO(MYID)
        DO I = 1, NCL1
            DO K = 1, NCL3
                QN(I, J, K) = QINJNKN(I, J, K)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!C *************  SUBROUTINE INTRPRI ************
!C
!C  THIS ROUTINE INTERPOLATES ALONG X1 COORDINATES AT THE
!C  OLD POSITIONS OF X3  AND X2
!C
SUBROUTINE INTRPRI(QO, QN,XOLD,XNEW, NIOLD, NINEW, KO)
    USE WPRECISION
    !USE mpi_info !test
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NIOLD, NINEW, KO
    REAL(WP), INTENT(IN) :: QO(NIOLD, KO) ! ORIGINAL coarse
    REAL(WP), INTENT(OUT) :: QN(NINEW, KO) ! NEW FINE
    REAL(WP), INTENT(IN) :: XOLD(NIOLD)
    REAL(WP), INTENT(IN) :: XNEW(NINEW)

    REAL(WP) :: FOLD(NIOLD)
    REAL(WP) :: YU(NIOLD)
    REAL(WP) :: DPN2, DP1
    REAL(WP) :: XX, YY
    INTEGER(4) :: I, K

    DO K = 1, KO
        DO I = 1, NIOLD
            FOLD(I) = QO(I, K)
        END DO
        DPN2 = 0.0_WP
        DP1 = 0.0_WP

        CALL SPLINE(XOLD, FOLD, NIOLD, DP1, DPN2, YU)
        DO I = 1, NINEW  !!1, NINEW
            XX = XNEW(I)
            CALL SPLINT(XOLD, FOLD, YU, NIOLD,XX, YY)
            QN(I, K) = YY
            !WRITE(*, *) MYID, K, I, XNEW(I), QN(I, K) !TEST
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
!C
!C IN THIS SUBR EVALUATE THE VELOCITY, READ FROM A CONTINUATION
!C FILE, BY AN CUBIC SPLINE INTERPOLATION IN THE NEW GRID
!C
!C
!C
SUBROUTINE INTRPRK(QO, QN,XOLD,XNEW, NKOLD, NKNEW, IN)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NKOLD, NKNEW, IN
    REAL(WP), INTENT(IN) :: QO(IN, NKOLD)
    REAL(WP), INTENT(IN) :: XOLD(NKOLD)
    REAL(WP), INTENT(IN) :: XNEW(NKNEW)
    REAL(WP), INTENT(out) :: QN(IN, NKNEW)

    INTEGER(4) :: I, K
    REAL(WP) :: FOLD(NKOLD)
    REAL(WP) :: DPN22, DP12, XX, YY
    REAL(WP) :: YU(NKOLD)

    DO I = 1, IN
        DO K = 1, NKOLD
            FOLD(K) = QO(I, K)
        END DO
        DPN22 = 0.0_WP
        DP12 = 0.0_WP

        CALL SPLINE(XOLD, FOLD, NKOLD, DP12, DPN22, YU)

        DO K = 1, NKNEW
            XX = XNEW(K)
            CALL SPLINT(XOLD, FOLD, YU, NKOLD,XX, YY)
            QN(I, K) = YY
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INTRPRJ(QO, QN,XOLD,XNEW, NJOLD, NJNEW, IN)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NJOLD, NJNEW, IN
    REAL(WP), INTENT(IN) :: QO(IN, NJOLD)
    REAL(WP), INTENT(IN) :: XOLD(NJOLD)
    REAL(WP), INTENT(IN) :: XNEW(NJNEW)
    REAL(WP), INTENT(out) :: QN(IN, NJNEW)

    INTEGER(4) :: I, J
    REAL(WP) :: FOLD(NJOLD)
    REAL(WP) :: DPN22, DP12, XX, YY
    REAL(WP) :: YU(NJOLD)

    DO I = 1, IN
        DO J = 1, NJOLD
            FOLD(J) = QO(I, J)
        END DO
        DPN22 = 0.0_WP
        DP12 = 0.0_WP

        CALL SPLINE(XOLD, FOLD, NJOLD, DP12, DPN22, YU)

        DO J = 1, NJNEW
            XX = XNEW(J)
            CALL SPLINT(XOLD, FOLD, YU, NJOLD,XX, YY)
            QN(I, J) = YY
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
!C *************  SUBROUTINE INTRPRJG ************
!C
!C  THIS ROUTINE INTERPOLATES ALONG X2 COORDINATES AT THE
!C  NEW POSITIONS OF X3 AND AT THE NEW POSITIONS OF X1
!C
SUBROUTINE INTRPRJG(QO, QN,XOLD,XNEW, N2DOO, N2DO, NJOLD, NJNEW, IN, KN)
    USE MPI_INFO
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NJOLD, NJNEW, IN, KN
    INTEGER(4), INTENT(IN) :: N2DOO(0 : NPSLV), N2DO(0 : NPSLV)
    REAL(WP), INTENT(IN) :: QO(IN, N2DOO(0), KN)
    REAL(WP), INTENT(OUT) :: QN(IN, N2DO(0), KN)
    REAL(WP), INTENT(IN) :: XOLD(NJOLD)
    REAL(WP), INTENT(IN) :: XNEW(NJNEW)


    INTEGER(4) :: N3DO(0 : NPSLV)
    REAL(WP) :: TF1(IN, NJOLD, KN / NPTOT)
    REAL(WP) :: TF2(IN, NJNEW, KN / NPTOT)
    REAL(WP) :: FOLD(NJOLD)
    REAL(WP) :: YU  (NJOLD)
    REAL(WP) :: DPN2, DP1, XX, YY
    INTEGER(4) :: K, I, J

    N3DO(:) = KN / NPTOT
    CALL TRASPO23_Y2Z(IN, NJOLD, KN, 1, N2DOO(0), N2DOO, N3DO, QO, TF1)

    DO K = 1, N3DO(MYID)
        DO I = 1, IN

            DO J = 1, NJOLD
                FOLD(J) = TF1(I, J, K)
            END DO
            DPN2 = 0.0_WP
            DP1 = 0.0_WP

            CALL SPLINE(XOLD, FOLD, NJOLD, DP1, DPN2, YU)

            DO J = 1, NJNEW
                XX = XNEW(J)
                CALL SPLINT(XOLD, FOLD, YU, NJOLD,XX, YY)
                TF2(I, J, K) = YY
            END DO

        END DO
    END DO

    CALL TRASPO23_Z2Y(IN, NJNEW, KN, 1, N2DO(0), N2DO, N3DO, QN, TF2)

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE SPLINE(XIN, Y, N,RP1,RPN, Y2)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N
    REAL(WP), INTENT(IN) :: XIN(N)
    REAL(WP), INTENT(IN) :: Y(N)
    REAL(WP), INTENT(IN) :: RP1, RPN
    REAL(WP), INTENT(OUT):: y2(N)

    INTEGER(4) :: I, K
    REAL(WP) :: SIG, PNN, QN, UN
    REAL(WP) :: UU1(N)


    IF (RP1 >  .99E30) THEN
        Y2(1) = 0.0_WP
        UU1(1) = 0.0_WP
    ELSE
        Y2(1) = -0.50_WP
        UU1(1) = (3.0_WP / (XIN(2) -XIN(1))) * ((Y(2) - Y(1)) / (XIN(2) -XIN(1)) - RP1)
    ENDIF

    DO I = 2, n - 1
        SIG = (XIN(I) -XIN(I - 1)) / (XIN(I + 1) -XIN(I - 1))
        PNN = SIG* Y2(I - 1) + 2.0_WP
        Y2(I) = (SIG - 1.0_WP) /PNN
        UU1(I) = (6.0_WP * ((Y(I + 1) - Y(I)) / (XIN(I + 1) - XIN(I)) - &
        (Y(I) - Y(I - 1)) / (XIN(I) - XIN(I - 1))) / (XIN(I + 1) - XIN(I - 1)) -SIG* UU1(I - 1)) / PNN
    END DO

    IF (RPN >  .99E30) THEN
        QN = 0.0_WP
        UN = 0.0_WP
    ELSE
        QN = 0.50_WP
        UN = (3.0_WP / (XIN(N) - XIN(N - 1))) * (RPn -(Y(N) - Y(N - 1)) / (XIN(N) - XIN(N - 1)))
    ENDIF

    Y2(N) = (Un - QN * UU1(N - 1)) / (QN * Y2(N - 1) + 1.0_WP)
    DO K = n - 1, 1, -1
        Y2(K) = Y2(K) * Y2(K + 1) + UU1(K)
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE SPLINT(XAIN, YA, Y2A, N,X11, Y11)
    USE WPRECISION
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: N
    REAL(WP), INTENT(IN) :: XAIN(N)
    REAL(WP), INTENT(IN) :: YA(N)
    REAL(WP), INTENT(IN) :: Y2A(N)
    REAL(WP), INTENT(IN) :: X11
    REAL(WP), INTENT(OUT) :: Y11

    INTEGER(4) :: KLO, KHI, K
    REAL(WP) :: H1, A1, B1

    KLO= 1
    KHI = N
    DO WHILE(KHI - KLO >  1)
        K = (KHI + KLO) / 2
        IF(XAIN(K) >  X11)THEN
            KHI = K
        ELSE
            KLO= K
        ENDIF
    END DO

    !1       IF (KHI - KLO >  1) THEN
    !            K = (KHI + KLO) / 2
    !            IF(XAIN(K) >  X11)THEN
    !                KHI = K
    !            ELSE
    !                KLO= K
    !            ENDIF
    !            GOTO 1
    !        ENDIF

    H1 = XAIN(KHI) -XAIN(KLO)
    !IF (H1 == 0.) PAUSE 'BAD XA INPUT.'
    A1 = (XAIN(KHI) -X11) /H1
    B1 = (X11-XAIN(KLO)) /H1
    Y11 = A1 * YA(KLO) +B1 * YA(KHI) + ((A1**3- A1) * Y2A(KLO) + (B1**3-B1) * Y2A(KHI)) * (H1**2) /6.0_WP

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TRASPO23_Y2Z(NCL1, NCL2, NCL3, NYS, NYE, N2DO, N3DO, TRHS, TF)
    USE WPRECISION
    USE mpi_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NCL1, NCL2, NCL3, NYS, NYE
    INTEGER(4), INTENT(IN) :: N2DO(0 : NPSLV), N3DO(0 : NPSLV)
    REAL(WP),   INTENT(IN) :: TRHS(NCL1, NYS : NYE, NCL3    )
    REAL(WP),   INTENT(OUT) :: TF  (NCL1, NCL2,   N3DO(0) )

    REAL(WP) :: SENBLK( NCL1, N2DO(0), N3DO(0) )
    REAL(WP) :: RECBLK( NCL1, N2DO(0), N3DO(0) )
    INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NTOL
    INTEGER(4) :: JL, JG
    INTEGER(4) :: IG
    INTEGER(4) :: KL, KG
    INTEGER(4) :: IP
    INTEGER(4) :: ISEND
    INTEGER(4) :: IRECV
    INTEGER(4) :: ITAG
    INTEGER(4) :: NNTR

    NNTR = NCL1

    DO IP = 0, NPSLV
        IF(MYID == IP) THEN

            DO JL = 1, N2DO(MYID)
                JG=JL+ Ip * N2DO(0)
                DO KL = 1, N3DO(MYID)
                    KG = KL+ Ip * N3DO(0)
                    DO IG = 1, NNTR
                        TF(IG, JG, KL) = TRHS(IG, JL, KG)
                    END DO
                END DO
            END DO

        ELSE

            NTOL = NNTR * N2DO(0) * N3DO(0)
            DO KL = 1, N3DO(IP)
                KG = KL+ Ip * N3DO(0)
                DO IG = 1, NNTR
                    DO JL = 1, N2DO(MYID)
                        SENBLK(IG, JL, KL) = TRHS(IG, JL, KG)
                    END DO
                END DO
            END DO

            ISEND = IP
            IRECV = IP
            ITAG = 0
            CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
            IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
            ISEND, ITAG, ICOMM, TRAS_STS, IERROR)

            DO JL = 1, N2DO(IP)
                JG=JL+ Ip * N2DO(0)
                DO IG = 1, NNTR
                    DO KL = 1, N3DO(MYID)
                        TF(IG, JG, KL) = RECBLK(IG, JL, KL)
                    END DO
                END DO
            END DO

        END IF
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE TRASPO23_Z2Y(NCL1, NCL2, NCL3, NYS, NYE, N2DO, N3DO, TRHS, TF)
    USE WPRECISION
    USE mpi_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NCL1, NCL2, NCL3, NYS, NYE
    INTEGER(4), INTENT(IN) :: N2DO(0 : NPSLV), N3DO(0 : NPSLV)
    REAL(WP),   INTENT(OUT) :: TRHS(NCL1, NYS : NYE, NCL3    )
    REAL(WP),   INTENT(IN) :: TF  (NCL1, NCL2,   N3DO(0) )

    REAL(WP) :: SENBLK(NCL1, N2DO(0), N3DO(0))
    REAL(WP) :: RECBLK(NCL1, N2DO(0), N3DO(0))
    INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NTOL
    INTEGER(4) :: JL, JG
    INTEGER(4) :: IG
    INTEGER(4) :: KL, KG
    INTEGER(4) :: IP
    INTEGER(4) :: ISEND
    INTEGER(4) :: IRECV
    INTEGER(4) :: ITAG
    INTEGER(4) :: NNTR



    NNTR = NCL1
    DO IP = 0, NPSLV
        IF(MYID == IP) THEN

            DO KL = 1, N3DO(MYID)
                KG = KL+ Ip * N3DO(0)
                DO JL = 1, N2DO(MYID)
                    JG=JL+ Ip * N2DO(0)
                    DO IG = 1, NNTR
                        TRHS(IG, JL, KG) = TF(IG, JG, KL)
                    END DO
                END DO
            END DO

        ELSE
            NTOL = NNTR * N3DO(0) * N2DO(0)
            DO JL = 1, N2DO(IP)
                JG=JL+ Ip * N2DO(0)
                DO IG = 1, NNTR
                    DO KL = 1, N3DO(MYID)
                        SENBLK(IG, JL, KL) = TF(IG, JG, KL)
                    END DO
                END DO
            END DO
            ISEND = IP
            IRECV = IP
            ITAG = 0
            CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
            IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
            ISEND, ITAG, ICOMM, TRAS_STS, IERROR)


            DO KL = 1, N3DO(IP)
                KG = KL+ Ip * N3DO(0)
                DO IG = 1, NNTR
                    DO JL = 1, N2DO(MYID)
                        TRHS(IG, JL, KG) = RECBLK(IG, JL, KL)
                    END DO
                END DO
            END DO

        END IF
    END DO

    RETURN
END SUBROUTINE
