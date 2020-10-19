!**********************************************************************************************************************************
!> @brief
!>       Transformation from Y /Z to Z/Y
!> @details
!> SUBROUTINE: TRASP23_Y2Z (in MYID = all)
!> SUBROUTINE: TRASP23_Z2Y (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE TRASP23_Y2Z(NCL1, NYS, NYE, TRHS, TF)
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NCL1, NYS, NYE
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
                JG=JCL2G(JL)
                DO KL = 1, N3DO(MYID)
                    KG = KDSWT(MYID) - 1 + KL
                    DO IG = 1, NNTR
                        TF(IG, JG, KL) = TRHS(IG, JL, KG)
                    END DO
                END DO
            END DO

        ELSE

            NTOL = NNTR * N2DO(0) * N3DO(0)
            DO KL = 1, N3DO(IP)
                KG = KDSWT(IP) - 1 + KL
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
                JG=JDSWT(IP) + JL- 1
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
SUBROUTINE TRASP23_Z2Y(NCL1, NYS, NYE, TRHS, TF)
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NCL1, NYS, NYE
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
                KG = KDSWT(MYID) - 1 + KL
                DO JL = 1, N2DO(MYID)
                    JG=JDSWT(MYID) - 1 + JL
                    DO IG = 1, NNTR
                        TRHS(IG, JL, KG) = TF(IG, JG, KL)
                    END DO
                END DO
            END DO

        ELSE
            NTOL = NNTR * N3DO(0) * N2DO(0)
            DO JL = 1, N2DO(IP)
                JG=JDSWT(IP) - 1 + JL
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
                KG = KDSWT(IP) + KL- 1
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
