!**********************************************************************************************************************************
!> @brief
!>        To calculate the CFL number
!> @details
!> SUBROUTINE: CFL_tg (in MYID = all)
!> SUBROUTINE: CFL_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010- Initial Version (tg only), by Mehdi Seddighi
! 12/2013- Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CFL_tg
    USE flow_info
    USE init_info
    USE mesh_info

    IMPLICIT NONE

    REAL(WP) :: CFLMA
    REAL(WP) :: CFLM_WORK
    REAL(WP) :: QCF!, COE1
    INTEGER(4) :: I, IP
    INTEGER(4) :: J, JP, JJ, JC
    INTEGER(4) :: K, KP

    CFLMA = 0.0_WP
    CFLM_WORK = 0.0_WP
    QCF = 0.0_WP

    DO K = 1, NCL3
        KP = KPV(K)
        DO J = 1, N2DO(MYID)
            JP = JLPV(J)
            JJ = JCL2G(J)
            JC = JJ !NCL2 / 2
            DO I = 1, NCL1_tg
                IP = IPV_tg(I)
                QCF = ( DABS( (Q_tg(I, J, K, 1) + Q_tg(IP, J, K, 1)) * DXI                ) +        &
                DABS( (Q_tg(I, J, K, 2) + Q_tg(I, JP, K, 2)) * DYFI(JC) * RCCI1(JC) ) +        &
                DABS( (Q_tg(I, J, K, 3) + Q_tg(I, J, KP, 3)) * DZI * RCCI2(JC)      ) ) * 0.50_WP
                CFLMA = DMAX1(CFLMA, QCF)
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CFLMM_tg = CFLM_WORK

    RETURN
END SUBROUTINE

!**********************************************************************
SUBROUTINE CFL_io
    !>    @note
    !>    1) Energy equation has the same convection velocity as the
    !>       momentum equation, thus, no additional CFL is required
    !>       for energy equation.
    !>    2) For conservative variables, the convection variable are Q

    USE flow_info
    USE init_info
    USE mesh_info
    USE thermal_info

    IMPLICIT NONE

    REAL(WP) :: CFLMA
    REAL(WP) :: CFLM_WORK
    REAL(WP) :: QCF   !, COE1
    REAL(WP) :: VCF   !, COE1
    INTEGER(4) :: I, IP
    INTEGER(4) :: J, JP, JJ!, JC
    INTEGER(4) :: K, KP


    CFLMA = 0.0_WP
    CFLM_WORK = 0.0_WP
    QCF = 0.0_WP

    DO K = 1, NCL3
        KP = KPV(K)
        DO J = 1, N2DO(MYID)
            JP = JLPV(J)
            JJ = JCL2G(J)
            DO I = 1, NCL1_io
                IP = IPV_io(I)
                QCF = ( DABS( (Q_io(I, J, K, 1) + Q_io(IP, J, K, 1)) * DXI                 ) +        &
                DABS( (Q_io(I, J, K, 2) + Q_io(I, JP, K, 2)) * DYFI(JJ) * RCCI1(JJ)  ) +        &
                DABS( (Q_io(I, J, K, 3) + Q_io(I, J, KP, 3)) * DZI * RCCI2(JJ)       )   ) * 0.50_WP

                VCF = 2.0_WP * CVISC * (DXQI + DZQI * RCCI2(JJ) + DYFI(JJ) * DYFI(JJ) ) * Viscousity(I, J, K) / DENSITY(I, J, K)

                QCF = DMAX1(QCF, VCF)
                CFLMA = DMAX1(CFLMA, QCF)
                !>                  @warning Why rm is involved in 3rd direction?
            END DO
        END DO
    END DO

    CALL MPI_BARRIER(ICOMM, IERROR)
    CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
    CFLMM_io = CFLM_WORK

    RETURN
END SUBROUTINE

!**********************************************************************
SUBROUTINE CFL_VISCOUS ! not used any more, merged above to CFL_io
    USE flow_info
    USE init_info
    USE mesh_info
    USE thermal_info

    IMPLICIT NONE


    INTEGER(4) :: J, JJ, I, K
    REAL(WP) :: TEMP

    CFLVIS = 0.0_WP
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        DO K = 1, NCL3
            DO I = 1, NCL1_io
                TEMP = 2.0_WP * CVISC * (DXQI + DZQI * RCCI2(JJ) + DYFI(JJ) * DYFI(JJ) ) * Viscousity(I, J, K) / DENSITY(I, J, K)
                CFLVIS(J) = DMAX1(CFLVIS(J), TEMP)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
