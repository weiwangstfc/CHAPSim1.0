!**********************************************************************************************************************************
!> @brief
!>       To calculate coefficients of the PoISsion equation (Tridiagonal System of Equation)
!> @details
!> SUBROUTINE: TDMA_COEF (in MYID = 0)
!>             To calculate below parameters:
!> @pARam AMPH(1 : NND2M)   only USEd in  SUBROUTINE PRCALC
!> @pARam APPH(1 : NND2M)   only USEd in  SUBROUTINE PRCALC
!> @pARam ACPH(1 : NND2M)   only USEd in  SUBROUTINE PRCALC
!> @pARam AMVR(1 : NND2M)   \phi_{I - 1} Coefficient of Possion Eq. in UnIForm and non - UnIForm Eq.3.13 Mehdi thesis
!> @pARam APVR(1 : NND2M)   \phi_{i}   Coefficient of Possion Eq. in UnIForm and non - UnIForm Eq.3.13 Mehdi thesis
!> @pARam ACVR(1 : NND2M)   \phi_{I + 1} Coefficient of Possion Eq. in UnIForm and non - UnIForm Eq.3.13 Mehdi thesis
!> SUBROUTINE: BCAST_TDMA_COEF (in MYID = all)
!>             To broadcast the common information
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 02/2014 - Added io domain, optimized the code structure in f90 and code structure was optimized, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE TDMA_COEF
    !>      TEST the same as Ke Kui version in both Channel and pipe.
    !>       For WALL b.c.: u, v, w   on the wall
    !>                    : thermal properties on the wall
    !>                    : pressure on the ghost cells
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: JC, JM, JP, IDR
    INTEGER(4) :: NQ
    INTEGER(4) :: NFIL
    REAL(WP) :: A22ICC
    REAL(WP) :: A22ICP
    REAL(WP) :: A22
    REAL(WP) :: A22DP
    REAL(WP) :: UCAJ11, UCAJ
    REAL(WP) :: UGMMM
    REAL(WP) :: AC2

    IF(MYID /= 0) RETURN

    !================Cartesian ============================================
    IF(iCase == ICHANL) THEN
        DO JC = 1, NCL2
            JP = JC + 1

            AMPH(JC) = DYFI(JC) * DYCI(JC)
            APPH(JC) = DYFI(JC) * DYCI(JP)

            IF(JC == 1)    AMPH(JC) = DYFI(JC) * DYCI(JC) * 0.5_WP ! for b.c. phi in ghost cell
            IF(JP == NND2) APPH(JC) = DYFI(JC) * DYCI(JP) * 0.5_WP

        END DO
        ! set b.c. for phI =====
        AMPH0 = AMPH(1)
        APPH0 = APPH(NCL2)
        ACPH(1 : NCL2) = -( AMPH(1 : NCL2) + APPH(1 : NCL2) )
        !================
        !pressure gradient equal to zero at wall
        ACPH(1) = AMPH(1) + ACPH(1)
        AMPH(1) = 0.0_WP

        ACPH(NCL2) = APPH(NCL2) + ACPH(NCL2)
        APPH(NCL2) = 0.0_WP

        !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
        DO IDR = 1, 3
            IF( (IDR == 1) .OR. (IDR == 3) ) THEN
                DO JC = 2, NCL2 - 1
                    JP = JC + 1
                    AMVR(JC, IDR) = DYFI(JC) * DYCI(JC)
                    APVR(JC, IDR) = DYFI(JC) * DYCI(JP)
                    ACVR(JC, IDR) = -( AMVR(JC, IDR) + APVR(JC, IDR) )
                END DO
                !>              @note : FOR JC = 1   ! to DO for u,w on the wall!!! check
                UCAJ = 4.0_WP / ( 2.0_WP / DYCI(2) + 1.0_WP / DYFI(1) )
                AMVR(1, IDR) =   0.0_WP
                APVR(1, IDR) =   UCAJ * DYCI(2)
                ACVR(1, IDR) =  - UCAJ * ( DYCI(2) + 2.0_WP * DYFI(1) )
                AMVR1 = UCAJ * (2.0_WP * DYFI(1))
                !>               @note : FOR JC = NND2M ! to DO for u,w on the wall!!!  check
                UCAJ = 4.0_WP / ( 2.0_WP / DYCI(NCL2) + 1.0_WP / DYFI(NCL2) )
                AMVR(NCL2, IDR) = UCAJ * DYCI(NCL2)
                APVR(NCL2, IDR) = 0.0_WP
                ACVR(NCL2, IDR) = -UCAJ * ( 2.0_WP * DYFI(NCL2) + DYCI(NCL2) )
                APVRN = UCAJ * (2.0_WP * DYFI(NCL2))

            ELSE
                !>          @note: clustered grids in non - Homogenous direction (y), the same methods as above.
                !>                 coefficients in Eq.3.13 Page 44 of Mehdi thesis.
                DO JC = 2, NCL2
                    JM = JC - 1
                    AMVR(JC, IDR) =  DYFI(JM) * DYCI(JC)
                    APVR(JC, IDR) =  DYFI(JC) * DYCI(JC)
                    ACVR(JC, IDR) = -( AMVR(JC, IDR) + APVR(JC, IDR) )
                END DO
                AMVR(1, IDR) = 0.0_WP
                APVR(1, IDR) = 0.0_WP
                ACVR(1, IDR) = 1.0_WP
                AMVR(NND2, IDR) = 0.0_WP
                APVR(NND2, IDR) = 0.0_WP
                ACVR(NND2, IDR) = 1.0_WP

            ENDIF

        END DO
        !END IF

    END IF


    !================cylindrical ============================================
    IF(iCase == iPIPEC .OR. iCase == iANNUL) THEN
        ! for pressure correctioN ==========
        DO JC = 2, NCL2 - 1
            JP = JC + 1
            A22ICC = DYCI(JC) / RNDI1(JC) !rc(JC)
            A22ICP = DYCI(JP) / RNDI1(JP) !rc(JP)
            AC2 = -(A22ICC + A22ICP)
            UGMMM = DYFI(JC) / RCCI1(JC)       !rm(JC) *
            AMPH(JC) = (A22ICC) * UGMMM
            APPH(JC) = (A22ICP) * UGMMM
            ACPH(JC) = AC2 * UGMMM
        END DO

        JC = 1
        JP = JC + 1
        A22ICP = DBLE(JP - JC) * DYCI(JP) / RNDI1(JP)     !rc(JP) *
        UGMMM = DYFI(JC) / RCCI1(JC) !rm(JC) *
        AMPH(JC) = 0.0_WP
        APPH(JC) = (A22ICP) * UGMMM
        ACPH(JC) = -(A22ICP) * UGMMM

        JC = NCL2
        A22ICC = DYCI(JC) / RNDI1(JC)  !rc(JC) *
        UGMMM = DYFI(JC) / RCCI1(JC)  !rm(JC) *
        AMPH(JC) = (A22ICC) * UGMMM
        APPH(JC) = 0.0_WP
        ACPH(JC) = -(A22ICC) * UGMMM

        !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
        DO NQ = 1, 3
            IF(NQ == 3) THEN
                DO JC = 2, NCL2 - 1
                    JP = JC + 1
                    JM = JC - 1     !@
                    A22ICC = DYCI(JC) / RNDI1(JC)          !rc(JC) *
                    A22ICP = DYCI(JP) / RNDI1(JP)          !rc(JP) *
                    AMVR(JC, NQ) =  A22ICC * DYFI(JC) * RCCI1(JM)      !/ Rm(JM)
                    APVR(JC, NQ) =  A22ICP * DYFI(JC) * RCCI1(JP)      !/ Rm(JP)
                    ACVR(JC, NQ) = -(A22ICC + A22ICP) * DYFI(JC) * RCCI1(JC) - RCCI2(JC)!1.0_WP / Rm(JC)**2
                END DO

                JC = 1
                JP = JC + 1
                UCAJ11 = DYFI(JC) * DYCI(JP) / RNDI1(JP)      !* Rc(JP)
                AMVR(JC, NQ) = 0.0_WP
                APVR(JC, NQ) = UCAJ11 * RCCI1(JP)!/ Rm(JP)
                ACVR(JC, NQ) = -UCAJ11 * RCCI1(JC) - RCCI2(JC)!1.0_WP / Rm(JC)**2

                JC = NCL2
                JP = NND2
                JM = JC - 1
                UGMMM = DYFI(JC) / RNDI1(JC)!* Rm(JC)

                AMVR(JC, NQ) = UGMMM * DYCI(JC) * RNDI1(JC)    !/ Rc(JC)
                APVR(JC, NQ) = 0. !- UGMMM / Rc(JP) / Cac(JP) * 2.0_WP
                !ACVR(JC, NQ) = -UGMMM * DYCI(JC) * RNDI1(JC) - UGMMM * DYCI(JP) * 2.0_WP * RNDI1(JP)   ! W on the ghost cell
                ACVR(JC, NQ) = -UGMMM * DYCI(JC) * RNDI1(JC) - UGMMM * DYCI(JP) * RNDI1(JP)     ! W on the wall.
            endIF

            IF(NQ == 1) THEN
                DO JC = 2, NCL2 - 1
                    JP = JC + 1
                    UGMMM = DYFI(JC) * RCCI1(JC)!/ Rm(JC)
                    AMVR(JC, NQ) = + UGMMM * DYCI(JC) / RNDI1(JC)   !* Rc(JC)
                    APVR(JC, NQ) = + UGMMM * DYCI(JP) / RNDI1(JP)   !* Rc(JP)
                    ACVR(JC, NQ) = -(AMVR(JC, NQ) + APVR(JC, NQ))
                END DO

                JC = 1
                JP = JC + 1
                UGMMM = DYFI(JC) * RCCI1(JC)!/ Rm(JC)
                AMVR(JC, NQ) = 0.0_WP
                APVR(JC, NQ) = UGMMM * (DYCI(JP)) / RNDI1(JP)    !rc(JP) *
                ACVR(JC, NQ) = -UGMMM * (DYCI(JP)) / RNDI1(JP)      !rc(JP) *

                JC = NCL2
                JP = NND2
                UGMMM = DYFI(JC) * RCCI1(JC)!/ Rm(JC)
                AMVR(JC, NQ) = + UGMMM * DYCI(JC) / RNDI1(JC)      !rc(JC) *
                APVR(JC, NQ) = 0.0_WP !UGMMM * Rc(JP) / Cac(JP) * 2.0_WP
                !ACVR(JC, NQ) = -UGMMM * DYCI(JC) / RNDI1(JC) - UGMMM * DYCI(JP) * 2.0_WP / RNDI1(JP)      ! U on the ghost cell
                ACVR(JC, NQ) = -UGMMM * DYCI(JC) / RNDI1(JC) - UGMMM * DYCI(JP) / RNDI1(JP)  ! U on the wall.
            endIF

            IF(NQ == 2) THEN
                DO JC = 2, NCL2
                    JM = JC - 1
                    A22  = DYCI(JC)
                    A22DP = DYCI(JC) / (2.0_WP / RNDI1(JC))
                    APVR(JC, NQ) =  A22 * DYFI(JC) - A22DP
                    ACVR(JC, NQ) = -(A22 * DYFI(JC) + A22 * DYFI(JM))
                    AMVR(JC, NQ) =  A22 * DYFI(JM) + A22DP
                END DO

                DO JC = 1, NND2, NCL2
                    APVR(JC, NQ) = 0.0_WP
                    ACVR(JC, NQ) = 0.0_WP
                    AMVR(JC, NQ) = 0.0_WP
                END DO
            endIF

        END DO
        !END IF
    END IF



    !================Cartesian ============================================
    IF(iCase == IBox3P) THEN

        DO JC = 1, NCL2
            JP = JC + 1

            AMPH(JC) = DYFI(JC) * DYCI(JC)
            APPH(JC) = DYFI(JC) * DYCI(JP)

            IF(JC == 1)    AMPH(JC) = DYFI(JC) * DYCI(JC) ! DYCI IS whole length, not half at wall
            IF(JP == NND2) APPH(JC) = DYFI(JC) * DYCI(JP) ! DYCI IS whole length, not half at wall

        END DO
        ! set b.c. for phI =====
        AMPH0 = AMPH(1)
        APPH0 = APPH(NCL2)
        ACPH(1 : NCL2) = -( AMPH(1 : NCL2) + APPH(1 : NCL2) )


        !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
        DO IDR = 1, 3

            IF(IDR == 1 .OR. IDR == 3) THEN
                DO JC = 2, NCL2 - 1
                    JP = JC + 1
                    AMVR(JC, IDR) = DYFI(JC) * DYCI(JC)
                    APVR(JC, IDR) = DYFI(JC) * DYCI(JP)
                    ACVR(JC, IDR) = -( AMVR(JC, IDR) + APVR(JC, IDR) )
                END DO
                JC = 1
                AMVR(1, IDR) = DYFI(1) * DYCI(1)
                APVR(1, IDR) = DYFI(1) * DYCI(2)
                ACVR(1, IDR) = -( AMVR(1, IDR) + APVR(1, IDR) )
                JC = NCL2
                AMVR(NCL2, IDR) = DYFI(NCL2) * DYCI(NCL2)
                APVR(NCL2, IDR) = DYFI(NCL2) * DYCI(NND2)
                ACVR(NCL2, IDR) = -( AMVR(NCL2, IDR) + APVR(NCL2, IDR) )
            ELSE
                DO JC = 2, NCL2
                    JM = JC - 1
                    AMVR(JC, IDR) =  DYFI(JM) * DYCI(JC)
                    APVR(JC, IDR) =  DYFI(JC) * DYCI(JC)
                    ACVR(JC, IDR) = -( AMVR(JC, IDR) + APVR(JC, IDR) )
                END DO
                JC = 1
                AMVR(1, IDR) = DYFI(NCL2) * DYCI(1)
                APVR(1, IDR) = DYFI(1) * DYCI(1)
                ACVR(1, IDR) = -( AMVR(1, IDR) + APVR(1, IDR) )
                JC = NND2
                AMVR(NND2, IDR) =  DYFI(NCL2) * DYCI(NND2)
                APVR(NND2, IDR) =  DYFI(NCL2) * DYCI(NND2)
                ACVR(NND2, IDR) = -( AMVR(1, IDR) + APVR(1, IDR) )
            END IF


        END DO
        !END IF


        !            ACPH(:) = ACVR(:, 1)
        !            AMPH(:) = AMVR(:, 1)
        !            APPH(:) = APVR(:, 1)

    END IF

    !!**********************************************************************
    NFIL = 16
    OPEN(NFIL, FILE = TRIM(FilePath0) // 'CHK_COEF_TDMA.dat')
    !REWIND NFIL
    WRITE(NFIL, '(A)') &
    '#AMPH(JC),  ACPH(JC), APPH(JC)'
    DO JC = 1, NCL2
        WRITE(NFIL, '(3ES15.7)') &
        AMPH(JC),  ACPH(JC), APPH(JC)
    END DO

    !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
    WRITE(NFIL, '(A)') &
    '#AMVR(JC, 1), ACVR(JC, 1), APVR(JC, 1)'
    DO JC = 1, NCL2
        WRITE(NFIL, '(3ES15.7)') &
        AMVR(JC, 1), ACVR(JC, 1), APVR(JC, 1)
    END DO

    WRITE(NFIL, '(A)') &
    '#AMVR(JC, 2), ACVR(JC, 2), APVR(JC, 2)'
    DO JC = 1, NND2
        WRITE(NFIL, '(3ES15.7)') &
        AMVR(JC, 2), ACVR(JC, 2), APVR(JC, 2)
    END DO

    WRITE(NFIL, '(A)') &
    '#AMVR(JC, 3), ACVR(JC, 3), APVR(JC, 3)'
    DO JC = 1, NCL2
        WRITE(NFIL, '(3ES15.7)') &
        AMVR(JC, 3), ACVR(JC, 3), APVR(JC, 3)
    END DO

    !            !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
    !            WRITE(NFIL, '(A)') &
    !                   '#AMVR(JC, 1), ACVR(JC, 1), APVR(JC, 1), ',  &
    !                   '#AMVR(JC, 2), ACVR(JC, 2), APVR(JC, 2), ',  &
    !                   '#AMVR(JC, 3), ACVR(JC, 3), APVR(JC, 3)  '
    !            DO JC = 1, NCL2
    !                WRITE(NFIL, '(9ES15.7)') &
    !                    AMVR(JC, 1), ACVR(JC, 1), APVR(JC, 1),  &
    !                    AMVR(JC, 2), ACVR(JC, 2), APVR(JC, 2),  &
    !                    AMVR(JC, 3), ACVR(JC, 3), APVR(JC, 3)
    !            END DO
    !END IF
    CLOSE(NFIL)

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BCAST_TDMA_COEF
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    CALL MPI_BCAST( AMPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( ACPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( APPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !IF(TgFlowFlg .OR. (iVisScheme == VisImplicit)) THEN
    CALL MPI_BCAST( AMVR, NND2 * 3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( ACVR, NND2 * 3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( APVR, NND2 * 3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( AMPH0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( APPH0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    !END IF


    RETURN
END SUBROUTINE
