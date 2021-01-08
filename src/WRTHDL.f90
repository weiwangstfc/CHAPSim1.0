!**********************************************************************************************************************************
!> @brief
!>        to write information to the log/screen file
!> @details
!> SUBROUTINE: MKDIR
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE CHKHDL_FORMAT
    CHARACTER(LEN = 32) :: FORMAT_CHAR  = '(4X, A)'
    CHARACTER(LEN = 32) :: FORMAT_ERROR = '(4X, A, I3.1, 2X, A)'

    CHARACTER(LEN = 32) :: FORMAT_1INT_SHORT = '(A32, 2X, I11.1)'
    CHARACTER(LEN = 32) :: FORMAT_1INT_LONG  = '(A64, 2X, I11.1)'

    CHARACTER(LEN = 32) :: FORMAT_2INT_SHORT = '(A32, 2X, 2I11.1)'
    CHARACTER(LEN = 32) :: FORMAT_2INT_LONG  = '(A64, 2X, 2I11.1)'

    CHARACTER(LEN = 32) :: FORMAT_1FLT_SHORT = '(A32, 2X, ES17.9)'
    CHARACTER(LEN = 32) :: FORMAT_1FLT_LONG  = '(A64, 2X, ES17.9)'

    CHARACTER(LEN = 32) :: FORMAT_2FLT_SHORT = '(A32, 2X, 2ES17.9)'
    CHARACTER(LEN = 32) :: FORMAT_2FLT_LONG  = '(A64, 2X, 2ES17.9)'

    CHARACTER(LEN = 32) :: FORMAT_3FLT_SHORT = '(A32, 2X, 3ES17.9)'
    CHARACTER(LEN = 32) :: FORMAT_3FLT_LONG  = '(A64, 2X, 3ES17.9)'

END MODULE
!**********************************************************************************************************************************
SUBROUTINE MKDIR
    USE WRT_INFO
    IMPLICIT none

    CALL SYSTEM('mkdir -p ' // TRIM(DIR0))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR1))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR2))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR3))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR4))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR5))
    CALL SYSTEM('mkdir -p ' // TRIM(DIR6))

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE ERRHDL(MSG, RANK)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'
    WRITE(logflg_pg, TRIM(FORMAT_ERROR)) '# Error in MYID ', RANK, MSG

    STOP

END SUBROUTINE ERRHDL
!**********************************************************************************************************************************
!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     WRITE out the check message
!>**********************************************************************
SUBROUTINE CHKHDL(MSG, RANK)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A)') '# MYID ', RANK, MSG
    WRITE(logflg_pg, TRIM(FORMAT_CHAR)) MSG

END SUBROUTINE
!**********************************************************************************************************************************
!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     WRITE out the check message
!>**********************************************************************
SUBROUTINE CHKINTHDL(MSG, RANK, N)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    INTEGER :: N
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A, 5X, I11.1)') '# MYID ', RANK, MSG, N
    IF(LEN(TRIM(MSG)) < 32) THEN
        WRITE(logflg_pg, TRIM(FORMAT_1INT_SHORT)) MSG, N
    ELSE
        WRITE(logflg_pg, TRIM(FORMAT_1INT_LONG)) MSG, N
    END IF
    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CHK2INTHDL(MSG, RANK, N1, N2)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    INTEGER :: N1, N2
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A, 5X, 2I11.1)') '# MYID ', RANK, MSG, N1, N2
    IF(LEN(TRIM(MSG)) < 32) THEN
        WRITE(logflg_pg, TRIM(FORMAT_2INT_SHORT)) MSG, N1, N2
    ELSE
        WRITE(logflg_pg, TRIM(FORMAT_2INT_LONG)) MSG, N1, N2
    END IF
    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     WRITE out the check message
!>**********************************************************************
SUBROUTINE CHKRLHDL(MSG, RANK, A)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    DOUBLE PRECISION :: A
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A, 5X, ES17.9)') '# MYID ', RANK, MSG, A

    IF(LEN(TRIM(MSG)) < 32) THEN
        WRITE(logflg_pg, TRIM(FORMAT_1FLT_SHORT)) MSG, A
    ELSE
        WRITE(logflg_pg, TRIM(FORMAT_1FLT_LONG)) MSG, A
    END IF
    RETURN

END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE CHK2RLHDL(MSG, RANK, A1, A2)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    DOUBLE PRECISION :: A1, A2
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A, 5X, 2ES17.9)') '# MYID ', RANK, MSG, A1, a2
    IF(LEN(TRIM(MSG)) < 32) THEN
        WRITE(logflg_pg, TRIM(FORMAT_2FLT_SHORT)) MSG, A1, A2
    ELSE
        WRITE(logflg_pg, TRIM(FORMAT_2FLT_LONG)) MSG, A1, A2
    END IF
    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CHK3RLHDL(MSG, RANK, A1, A2, A3)
    USE WRT_INFO
    USE CHKHDL_FORMAT
    IMPLICIT NONE

    CHARACTER(*) :: MSG
    INTEGER :: RANK
    DOUBLE PRECISION :: A1, A2, A3
    !LOGICAL :: ex1, ex2
    !INQUIRE(FILE=TRIM(LOGFIL_PG), EXIST = ex1, OPENED = ex2)
    !if(.not. ex1) write(*, *) 'The FILE does not EXIST.'
    !if(.not. ex2) write(*, *) 'The FILE/UNIT cannot be connected.'

    !WRITE(logflg_pg, '(A, I3.1, 5X, A, 5X, 2ES17.9)') '# MYID ', RANK, MSG, A1, a2
    IF(LEN(TRIM(MSG)) < 32) THEN
        WRITE(logflg_pg, TRIM(FORMAT_3FLT_SHORT)) MSG, A1, A2, A3
    ELSE
        WRITE(logflg_pg, TRIM(FORMAT_3FLT_LONG)) MSG, A1, A2, A3
    END IF
    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_QP_tg
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    IMPLICIT NONE


    INTEGER(4) :: I, J, K, JJ, FLID
    CHARACTER(4) :: NNN

    WRITE(NNN, '(1I4.4)') MYID
    FLID = MYID + 20

    OPEN(FLID, FILE = 'OUT_UP' // NNN // '.debug',POSITION = "APPEND")
    WRITE(FLID, '(A)') '#=====================checK ======================= '
    DO J = 0, N2DO(MYID) + 1
        DO K = 1, NCL3
            DO I = 1, NCL1_TG
                WRITE(FLID, *) 'UGP', K, I, Q_tg(I, J, K, 1:3), PR_tg(I, J, K)
            END DO
        END DO
    END DO

    CLOSE(FLID)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_UWP_io
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    IMPLICIT NONE


    INTEGER(4) :: I, J, K, JJ, FLID
    CHARACTER(4) :: NNN

    CALL MPI_BARRIER(ICOMM, IERROR)
    WRITE(NNN, '(1I4.4)') MYID

    FLID = MYID + 200
    OPEN(FLID, FILE = 'OUT_UGP' // NNN // '.debug',POSITION = "APPEND")

    WRITE(FLID, '(A)') '#=====================checK ======================= '
    IF(MYID < (NPTOT/ 2)) THEN
        DO J = 1, N2DO(MYID)
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    WRITE(FLID, *) 'UGP', K, I, Q_io(I, J, K, 1:3), G_io(I, J, K, 1:3), PR_io(I, J, K)
                END DO
            END DO
        END DO
    ELSE
        DO J = N2DO(MYID), 1
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    WRITE(FLID, *) 'UGP', K, I, Q_io(I, J, K, 1:3), G_io(I, J, K, 1:3), PR_io(I, J, K)
                END DO
            END DO
        END DO

    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)
    CLOSE(FLID)

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_LOCAL(VAR, NJS, NJE,STR)
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: NJS, NJE
    REAL(WP), INTENT(IN) :: VAR(NCL1S : NCL1E, NJS : NJE, 1 : NCL3)
    CHARACTER(4) :: STR
    INTEGER(4) :: FLID

    INTEGER(4) :: I, J, K, JJ
    CHARACTER(4) :: NNN

    CALL MPI_BARRIER(ICOMM, IERROR)
    WRITE(NNN, '(1I4.4)') MYID

    FLID = MYID + 100
    OPEN(FLID, FILE = 'OUT_VARS' // NNN // '.debug',POSITION = "APPEND")

    WRITE(FLID, '(A)') '#=====================checK ======================= '
    IF(MYID < (NPTOT/ 2)) THEN
        DO J = NJS, NJE
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    !WRITE(MYID + 100, '(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I, J, K)
                    WRITE(FLID, *) STR, K, I, VAR(I, J, K)
                END DO
            END DO
        END DO
    ELSE
        DO J = NJE, NJS, -1
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    !WRITE(MYID + 100, '(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I, J, K)
                    WRITE(FLID, *) STR, K, I, VAR(I, J, K)
                END DO
            END DO
        END DO

    END IF


    CALL MPI_BARRIER(ICOMM, IERROR)
    CLOSE(FLID)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_DENSITY
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    USE thermal_info
    IMPLICIT NONE


    INTEGER(4) :: I, J, K, JJ
    CHARACTER(4) :: NNN

    CALL MPI_BARRIER(ICOMM, IERROR)
    WRITE(NNN, '(1I4.4)') MYID
    OPEN(MYID + 100, FILE = 'OUT_' // NNN // '.debug',POSITION = "APPEND")

    WRITE(MYID + 100, '(A, 1ES13.5)') 'check: t = ',PhyTIME
    DO J = 0, N2DO(MYID) + 1
        DO K = 1, NCL3
            DO I = 1, NCL1_io
                WRITE(MYID + 100, *) 'DENSITY', MYID, J, K, I, DENSITY(I, J, K), DENSITY0(I, J, K)
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_VG2_io
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    IMPLICIT NONE


    INTEGER(4) :: I, J, K, JJ
    CHARACTER(4) :: NNN


    CALL MPI_BARRIER(ICOMM, IERROR)

    WRITE(NNN, '(1I4.4)') MYID
    OPEN(MYID + 200, FILE = 'OUT_VG2A' // NNN // '.debug',POSITION = "APPEND")

    WRITE(MYID + 200, '(A)') '=====================checK ======================= '


    IF(MYID < (NPTOT/ 2)) THEN
        DO J = 1, N2DO(MYID) + 1, +1
            WRITE(MYID + 200, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    WRITE(MYID + 200, '(A, 2I4.1, 2ES15.7)') 'VG2', K, I, Q_io(I, J, K, 2), G_io(I, J, K, 2)
                END DO
            END DO
        END DO
    ELSE
        DO J = N2DO(MYID) + 1, 1, -1
            WRITE(MYID + 200, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    WRITE(MYID + 200, '(A, 2I4.1, 2ES15.7)') 'VG2', K, I, Q_io(I, J, K, 2), G_io(I, J, K, 2)
                END DO
            END DO
        END DO

    END IF

    CALL MPI_BARRIER(ICOMM, IERROR)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DEBUG_WRT_THERMAL
    USE FLOW_INFO
    USE INIT_INFO
    USE MESH_INFO
    USE thermal_info
    IMPLICIT NONE


    INTEGER(4) :: I, J, K, JJ, FLID
    CHARACTER(4) :: NNN
    CHARACTER(5) :: STR

    CALL MPI_BARRIER(ICOMM, IERROR)
    WRITE(NNN, '(1I4.4)') MYID
    STR = 'HDTMK'

    FLID = MYID + 100
    OPEN(FLID, FILE = 'OUT_VARS' // NNN // '.debug',POSITION = "APPEND")

    WRITE(FLID, '(A)') '#=====================checK ======================= '
    IF(MYID < (NPTOT / 2)) THEN
        DO J = 0, N2DO(MYID) + 1
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    !WRITE(MYID + 100, '(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I, J, K)
                    WRITE(FLID, *) STR, K, I, ENTHALPY(I, J, K), DENSITY(I, J, K), TEMPERATURE(I, J, K), &
                    Viscousity(I, J, K), THERMCONDT(I, J, K)
                END DO
            END DO
        END DO
    ELSE
        DO J = N2DO(MYID) + 1, 0, -1
            WRITE(FLID, '(A, I4.1)') 'J = ', J
            DO K = 1, NCL3
                DO I = NCL1S, NCL1E
                    !WRITE(MYID + 100, '(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I, J, K)
                    WRITE(FLID, *) STR, K, I, ENTHALPY(I, J, K), DENSITY(I, J, K), TEMPERATURE(I, J, K), &
                    Viscousity(I, J, K), THERMCONDT(I, J, K)
                END DO
            END DO
        END DO

    END IF


    CALL MPI_BARRIER(ICOMM, IERROR)
    CLOSE(FLID)

    RETURN
END SUBROUTINE
