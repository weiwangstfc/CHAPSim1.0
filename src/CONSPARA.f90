!**********************************************************************************************************************************
!> @brief
!>       setup constant PARAMETERs for time dIScretization method, like the RK method
!> @details
!> SUBROUTINE: CONSPARA (in MYID = 0)
!>             Calculate below information
!>             - Kronecker_Delta          \n
!>             - mesh spacing             \n
!>             - quadrant analysIS        \n
!> SUBROUTINE: BCAST_CONSPARA (in MYID = all)
!>             Broadcast the common information
!> @note
!> @toDO
! REVISION HISTORY:
! 05/ 2010- Initial Version (tg domAIn only), by Mehdi Seddighi
! 12 / 2013- added io domAIn, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CONSPARA
    USE init_info
    USE mesh_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    REAL(WP) :: DX1, DX2
    INTEGER(4) :: I, J

    IF(MYID /= 0) RETURN

    !=================
    DO I = 1, NDV
        DO J = 1, NDV
            IF(I ==J) THEN
                Kronecker_Delta(I, J) = 1
            ELSE
                Kronecker_Delta(I, J) = 0
            END IF
        END DO
    END DO
    !=================
    ALX1(1) = HX_tg
    ALX1(2) = HX_io

    IF(TgFlowFlg .AND. IoFlowFlg) THEN
        DX1 = ALX1(1) / DBLE(NCL1_TG)
        DX2 = ALX1(2) / DBLE(NCL1_io)
        IF(DABS(DX1 / DX2 - 1.0_WP) > 1.0E-8_WP) THEN
            CALL ERRHDL('# DX_tg /= DX_io', MYID)
        ELSE
            DX = DX1
        END IF
    ELSE
        IF(TgFlowFlg) THEN
            DX = ALX1(1) / DBLE(NCL1_tg)  !check toDO for clustering grids in x direction.
        END IF
        IF(IoFlowFlg) THEN
            DX = ALX1(2) / DBLE(NCL1_io)
        END IF
    END IF

    DXI  = 1.0_WP / DX
    DXQI = DXI * DXI


    IF(iCase == ICHANL .OR. iCase == IBox3P) THEN
        ALX3 = HZ
    ELSE IF (iCase == IPIPEC .OR. iCase == IANNUL) THEN
        ALX3 = 2.0_WP * PI
    ELSE
    endIF
    DZ = ALX3 / DBLE(NCL3)
    DZI  = 1.0_WP / DZ
    DZQI = DZI * DZI

    ALX2 = (HYT - HYB) / 2.0_WP   ! HALF Channel HEIGHT

    CVISC = 1.0_WP / REN ! updated by ReIni later.

    IF(TgFlowFlg) THEN
        VL1313_tg = 1.0_WP / DBLE(NCL1_tg * NCL3)
    END IF

    IF(IoFlowFlg) THEN
        VL1313_io = 1.0_WP / DBLE(NCL1_io * NCL3)
    END IF

    !========= Define threshold level H for quadrant analysIS ===============
    QUADHV(1) = 0.00_WP
    QUADHV(2) = 0.25_WP
    QUADHV(3) = 0.50_WP
    QUADHV(4) = 0.75_WP
    QUADHV(5) = 1.00_WP
    QUADHV(6) = 2.00_WP
    QUADHV(7) = 3.00_WP
    QUADHV(8) = 4.00_WP
    QUADHV(9) = 5.00_WP
    !


    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE BCAST_CONSPARA
    USE init_info
    USE mesh_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CALL MPI_BCAST( PI,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( ALX1, 2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( ALX2, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( ALX3, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( DX,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DXI,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DXQI, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( DZ,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DZI,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DZQI, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( CVISC, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( VL1313_tg, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( VL1313_io, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( Kronecker_Delta, 6, MPI_INTEGER4, 0, ICOMM, IERROR )

    CALL MPI_BCAST( QUADHV, QUADHN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    RETURN
END SUBROUTINE
