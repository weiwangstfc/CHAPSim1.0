!**********************************************************************************************************************************
!> @brief
!>        Initialize MPI
!> @details
!> SUBROUTINE: SOLVE (in myid =  all)
!> SUBROUTINE: BCAST_COMM_STEP (in myid =  all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE Start_mpi
    USE mpi_info
    IMPLICIT NONE

    INTEGER(4) :: NDIM11
    INTEGER(4) :: DIMS(1)
    INTEGER(4) :: ICOORDS(1)
    INTEGER(4) :: NMINUS
    INTEGER(4) :: NPLUS
    INTEGER(4) :: IDIMS
    LOGICAL     PERIODS(1), REORDER

    CALL MPI_INIT(IERROR)
    IF(IERROR /= 0) THEN
        CALL ERRHDL('mpi_init fails!', 0)
    END IF

    !>       @note FIND GLOBAL RANKING AND NUMBER OF PROCESSORS
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERROR)

    !*   PARAMETERS USED TO DEFINE periodic Cartesian PROCESSOR:
    NDIM11 = 1
    DIMS(1) = SIZE
    !
    !*     ALLOWS TO FORCE periodic NO. FOR PROCESSORS!!!(SO PROCESSOR AFTER
    !*     (SIZE - 1) IS 0, AND BEFORE 0 IS (SIZE - 1)
    PERIODS(1) = .TRUE.
    !*     ALLOWS THE HARDWARE TO OPTIMIZE PROCESSOR ARRANGEMENT
    REORDER = .TRUE.
    !
    !>       @note HARDWARE OPTIMIZED PROCESSOR CONFIGURATOR AND NEW RANKING
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIM11, DIMS, PERIODS,  &
    REORDER, ICOMM, IERROR)
    !>       @note ABOVE replace MPI_COMM_WORLD by ICOMM

    !>       GET MY POSITION IN THIS COMMUNICATOR AND MY NEIGHBORS
    CALL MPI_COMM_RANK(ICOMM, MYID, IERROR)
    CALL MPI_CART_SHIFT(ICOMM, 0, 1, NMINUS, NPLUS, IERROR)
    CALL MPI_CART_GET(ICOMM, 1, IDIMS,PERIODS, ICOORDS, IERROR)

    NPTOT = SIZE
    NPSLV = SIZE - 1

    RETURN
END SUBROUTINE
