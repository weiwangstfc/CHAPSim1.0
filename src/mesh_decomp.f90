!**********************************************************************************************************************************
!> @brief
!>       mesh decomPOSITION in y and z directions.
!> @details
!> SUBROUTINE: mesh_Ydecomp (in MYID = 0)
!>             Decomposite the mesh in y- direction
!> SUBROUTINE: mesh_Zdecomp (in MYID = 0)
!>             Decomposite the mesh in Z - direction
!> SUBROUTINE: mesh_Decomp_Bcast (in MYID = all)
!>             Broadcast the mesh decomPOSITION information
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 12 / 2013- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE mesh_Ydecomp
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4) :: NTOT1
    INTEGER(4) :: NLOCAL
    INTEGER(4) :: DEFICIT
    INTEGER(4) :: IP

    !>     @note GET THE SIZE OF THE PROBLEM = NTOT1
    !>           WE are CHOOSING 1- D DECOMPOSITION IN Y direction.
    !>           DecomPOSITION of mesh IS only in Y (or r) direction

    ALLOCATE ( JDSWT(0 : NPSLV) )
    ALLOCATE ( JDEWT(0 : NPSLV) )
    ALLOCATE ( N2DO (0 : NPSLV) )
    JDSWT = 0
    JDEWT = 0
    !>       @wARning: y velocity IS stored on y cell surface, whICh are N2M + 1
    NTOT1 = NCL2

    DO IP = 0, NPSLV

        NLOCAL  = NTOT1 / NPTOT
        DEFICIT = MOD(NTOT1, NPTOT)

        JDSWT(IP) = IP * NLOCAL + 1
        JDSWT(IP) = JDSWT(IP) + MIN(IP, DEFICIT)
        !
        IF (IP  <  DEFICIT) THEN
            NLOCAL = NLOCAL + 1
        ENDIF
        !
        JDEWT(IP) = JDSWT(IP) + NLOCAL - 1
        !
        IF ( (JDEWT(IP) > NTOT1) .OR. ( IP == NPSLV ) ) THEN
            JDEWT(IP) = NTOT1
        END IF
    END DO

    !CALL CHKHDL('Ydecomp,      RANK,  Start,  End,  Interval', MYID)
    DO IP = 0, NPSLV
        N2DO(IP) = JDEWT(IP) - JDSWT(IP) + 1
        !WRITE(*, '(A, 4I8.1)') '# MYID   0        Ydecomp, ', Ip, JDSWT(IP), JDEWT(IP), N2DO(IP)
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE mesh_Zdecomp
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4) :: NTOT1
    INTEGER(4) :: NLOCAL
    INTEGER(4) :: DEFICIT
    INTEGER(4) :: IP

    !>     @note GET THE NPTOT OF THE PROBLEM = NTOT1
    !>           WE are CHOOSING 1- D DECOMPOSITION IN Y direction.
    !>           DecomPOSITION of mesh IS only in Y (or r) direction

    ALLOCATE ( KDSWT(0 : NPSLV) )
    ALLOCATE ( KDEWT(0 : NPSLV) )
    ALLOCATE ( N3DO (0 : NPSLV) )
    KDSWT = 0
    KDEWT = 0

    !>       @wARning: y velocity IS stored on y cell surface, whICh are N2M + 1
    NTOT1 = NCL3

    DO IP = 0, NPSLV

        NLOCAL  = NTOT1 / NPTOT
        DEFICIT = MOD(NTOT1, NPTOT)

        KDSWT(IP) = IP * NLOCAL + 1
        KDSWT(IP) = KDSWT(IP) + MIN(IP, DEFICIT)
        !
        IF (IP  <  DEFICIT) THEN
            NLOCAL = NLOCAL + 1
        ENDIF
        !
        KDEWT(IP) = KDSWT(IP) + NLOCAL - 1
        !
        IF ( (KDEWT(IP) > NTOT1) .OR. ( IP == NPSLV ) ) THEN
            KDEWT(IP) = NTOT1
        END IF
    END DO

    !CALL CHKHDL('Zdecomp,      RANK,  Start,  End,  Interval', MYID)
    DO IP = 0, NPSLV
        N3DO(IP) = KDEWT(IP) - KDSWT(IP) + 1
        !WRITE(*, '(A, 4I8.1)') '# MYID   0        Zdecomp, ', Ip, KDSWT(IP), KDEWT(IP), N3DO(IP)
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE mesh_Decomp_Bcast
    USE mesh_info
    IMPLICIT NONE

    IF(MYID /= 0) THEN
        ALLOCATE ( JDSWT(0 : NPSLV) )
        ALLOCATE ( JDEWT(0 : NPSLV) )
        ALLOCATE ( KDSWT(0 : NPSLV) )
        ALLOCATE ( KDEWT(0 : NPSLV) )
        ALLOCATE ( N2DO (0 : NPSLV) )
        ALLOCATE ( N3DO (0 : NPSLV) )
        KDSWT = 0
        KDEWT = 0
        JDSWT = 0
        JDEWT = 0
        N2DO = 0
        N3DO = 0
    END IF

    CALL MPI_BCAST(JDSWT(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)
    CALL MPI_BCAST(JDEWT(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)

    CALL MPI_BCAST(KDSWT(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)
    CALL MPI_BCAST(KDEWT(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)

    CALL MPI_BCAST(N2DO(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)
    CALL MPI_BCAST(N3DO(0), NPTOT, MPI_INTEGER4, 0, ICOMM, IERROR)

    RETURN
END SUBROUTINE
