!**********************************************************************************************************************************
!> @brief
!>       To create the velocity profile (laminar Poiseuille Profile) U = f(y)
!> @details
!> SUBROUTINE: LAMPOISLPROF (in MYID = 0)
!>             To create the velocity profile (laminar Poiseuille Profile) U = f(y)
!>             To satISfy U_Bulk = 1.0.
!> SUBROUTINE: BCAST_LAMPOISLPROF (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 06/2011 - Added Pipe flow profile, by Kui He
! 12/2013 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
! 04/2016 - Added Annular and TGV-box meshing, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE LAMPOISLPROF
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: J, FLLG, k
    REAL(WP) :: Y1, Y2, PCOE, a, b, U1mean, U1maxx

    IF(MYID /= 0) RETURN
    IF(iCase == IBox3P) RETURN

    Y1 = HYB
    Y2 = HYT
    a = (HYT - HYB) / 2.0_WP
    b = 0.0_WP
    PCOE = 3.0_WP / 2.0_WP
    IF(iCase  == Ipipec) THEN
        PCOE = 2.0_WP
        a = 2.0_WP * A
    ELSE IF(iCase  == IANNUL) THEN
        PCOE = 2.0_WP
        b  = a + HYB
    ELSE
    END IF

    Vini(:) = 0.0_WP
    DO J = 1, NCL2
        IF ( YCC(J) < Y2 .AND. YCC(J) >  Y1 ) THEN
            Vini(J) = (1.0_WP - ((YCC(J) - b)**2) / a / a) * PCOE
        ELSE
            Vini(J) = 0.0_WP
        ENDIF
    END DO

    U1mean = 0.0_WP
    U1maxx = 0.0_WP
    DO J = 1, NCL2
        DO K = 1, NCL3
            U1mean = U1mean + Vini(J) / RCCI1(J) / DYFI(J) / DZI
            U1maxx = DMAX1(U1maxx, Vini(J))
        END DO
    END DO
    U1mean = U1mean / Area_inlet


    FLLG = 50
    OPEN(FLLG, FILE = TRIM(FilePath0) // 'CHK_INIL_Poiseuille.dat')
    WRITE(FLLG, '(A)') '#  YCC, Uini'
    DO J = 1, NCL2
        WRITE(FLLG, '(I8, 2ES15.7)') J, YCC(J), Vini(J)
    END DO
    CLOSE(FLLG)

    CALL CHKRLHDL  ('The given lamimar Poiseuille flow bulk velocitY =       ', MYID, U1mean)
    CALL CHKRLHDL  ('The given lamimar Poiseuille flow maxi velocitY =       ', MYID, U1maxx)


    RETURN

END  SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE BCAST_LAMPOISLPROF
    USE init_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    CALL MPI_BCAST( Vini, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    RETURN
END SUBROUTINE
