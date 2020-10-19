!**********************************************************************************************************************************
!> @brief
!>        To calculate the maximum velocity in three directions
!> @details
!> SUBROUTINE: VMAV_tg (in MYID = all)
!> SUBROUTINE: VMAV_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 12/2013 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE VMAV_tg

    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, JJ, NYI
    REAL(WP) :: VMA, VMI
    REAL(WP) :: VMAX_WORK, VMIN_WORK

    !>       @note Max. Q(i, J, K, 1)
    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    DO K = 1, NCL3
        DO J = 1, N2DO(MYID)
            DO I = 1, NCL1_tg
                VMA = DMAX1(VMA, DABS(Q_tg(I, J, K, 1)))
                VMI = DMIN1(VMI, DABS(Q_tg(I, J, K, 1)))
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)
    VMAX_tg(1) = VMAX_WORK
    VMIN_tg(1) = VMIN_WORK


    !>       @note Max. Q(i, J, K, 3)
    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    DO K = 1, NCL3
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            DO I = 1, NCL1_tg
                VMA = DMAX1( VMA, DABS( Q_tg(I, J, K, 3) * RCCI1(JJ) )   )
                VMI = DMIN1( VMI, DABS( Q_tg(I, J, K, 3) * RCCI1(JJ) )   )
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)
    VMAX_tg(3) = VMAX_WORK
    VMIN_tg(3) = VMIN_WORK

    !>        @note Max. Q(i, J, K, 2)
    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO K = 1, NCL3
        DO J = NYI, N2DO(MYID)
            JJ = JCL2G(J)
            DO I = 1, NCL1_tg
                VMA = DMAX1(VMA, DABS( Q_tg(I, J, K, 2) * RNDI1(JJ)  ) )
                VMI = DMIN1(VMI, DABS( Q_tg(I, J, K, 2) * RNDI1(JJ)  ) )
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)
    VMAX_tg(2) = VMAX_WORK
    VMIN_tg(2) = VMIN_WORK

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE VMAV_io
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, JJ, NYI
    REAL(WP) :: VMA, VMI
    REAL(WP) :: VMAX_WORK, VMIN_WORK


    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    DO K = 1, NCL3
        DO J = 1, N2DO(MYID)
            DO I = NCL1S, NCL1E
                VMA = DMAX1(VMA, DABS(Q_io(I, J, K, 1)))
                VMI = DMIN1(VMI, DABS(Q_io(I, J, K, 1)))
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
    VMAX_io(1) = VMAX_WORK
    VMIN_io(1) = VMIN_WORK


    !>      @note Max. Q(i, J, K, 3)
    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    DO K = 1, NCL3
        DO J = 1, N2DO(MYID)
            JJ = JCL2G(J)
            DO I = NCL1S, NCL1E
                VMA = DMAX1(  VMA, DABS(Q_io(I, J, K, 3) * RCCI1(JJ))  )
                VMI = DMIN1(  VMI, DABS(Q_io(I, J, K, 3) * RCCI1(JJ))  )
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
    VMAX_io(3) = VMAX_WORK
    VMIN_io(3) = VMIN_WORK

    !>     @note Max. Q(i, J, K, 2)
    VMA = -1.0e20_WP
    VMI = 1.0e20_WP
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO K = 1, NCL3
        DO J = NYI, N2DO(MYID)
            JJ = JCL2G(J)
            DO I = NCL1S, NCL1E
                VMA = DMAX1(VMA, DABS(Q_io(I, J, K, 2) * RNDI1(JJ)))
                VMI = DMIN1(VMA, DABS(Q_io(I, J, K, 2) * RNDI1(JJ)))
            END DO
        END DO
    END DO


    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
    CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
    VMAX_io(2) = VMAX_WORK
    VMIN_io(2) = VMIN_WORK

    RETURN

END SUBROUTINE
