!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> module: VARS_AVERAGED_nonXperiodic_io
!> subroutine: MEMO_ALLOCT_AVERAGE_nonXperiodic_io
!> subroutine: MEMO_DEALLT_AVERAGE_nonXperiodic_io
!> subroutine: MEMO_ALLOCT_INTP_nonXperiodic_io
!> subroutine: MEMO_DEALLT_INTP_nonXperiodic_io
!> subroutine: WRT_AVERAGE_PPED_nonXperiodic_io
!> subroutine: WRT_Reynolds_Average_Contour_nonXperiodic_io
!> subroutine: WRT_FARve_Average_Contour_nonXperiodic_io
!> subroutine: WRT_HeatTransfer_XProfile_nonXperiodic_io
!> subroutine: WRT_Cf_XProfile_nonXperiodic_io
!> subroutine: POSITIONS_nonXperiodic_io
!> subroutine: CL2ND_INTP_AVEG_nonXperiodic_io
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE VARS_AVERAGED_nonXperiodic_io
    USE thermal_info
    USE WRT_INFO

    CHARACTER(15) :: PNTIM

    INTEGER(4), PARAMETER :: NX = 5

    !============= Averaged global data in each processoR ======================
    REAL(WP), ALLOCATABLE :: U1ztL_F0_io( :, :, : )
    REAL(WP), ALLOCATABLE :: G1ztL_F0_io( :, :, : )

    REAL(WP), ALLOCATABLE :: UPztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: U2ztL_F0_io( :, :, : )
    REAL(WP), ALLOCATABLE :: UGztL_F0_io( :, :, : )

    REAL(WP), ALLOCATABLE :: UGUztL_F0_io(:, :, :)

    REAL(WP), ALLOCATABLE :: DVDL1ztL_F0_io( :, :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDLPztL_F0_io( :, :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDL2ztL_F0_io( :, :, :, :  )

    !=============== INTerp =====
    REAL(WP), ALLOCATABLE :: U1ztL_F0_INTP_io( :, :, : )
    REAL(WP), ALLOCATABLE :: G1ztL_F0_INTP_io( :, :, : )
    REAL(WP), ALLOCATABLE :: UPztL_F0_INTP_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: U2ztL_F0_INTP_io( :, :, : )
    REAL(WP), ALLOCATABLE :: UGztL_F0_INTP_io( :, :, : )

    REAL(WP), ALLOCATABLE :: UGUztL_F0_INTP_io(:, :, :)

    !============= Averaged global data in each processor thermaL ======================
    REAL(WP), ALLOCATABLE :: T1ztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: D1ztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: H1ztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: T2ztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: D2ztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: H2ztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: DHztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: UHztL_F0_io( :, :, : )
    REAL(WP), ALLOCATABLE :: GHztL_F0_io( :, :, : )
    !==== INTerP ==========
    REAL(WP), ALLOCATABLE :: T1ztL_F0_INTP_io( :, : )
    REAL(WP), ALLOCATABLE :: D1ztL_F0_INTP_io( :, : )
    REAL(WP), ALLOCATABLE :: H1ztL_F0_INTP_io( :, : )

    REAL(WP), ALLOCATABLE :: T2ztL_F0_INTP_io( :, : )
    REAL(WP), ALLOCATABLE :: D2ztL_F0_INTP_io( :, : )
    REAL(WP), ALLOCATABLE :: H2ztL_F0_INTP_io( :, : )

    REAL(WP), ALLOCATABLE :: DHztL_F0_INTP_io( :, : )

    REAL(WP), ALLOCATABLE :: UHztL_F0_INTP_io( :, :, : )
    REAL(WP), ALLOCATABLE :: GHztL_F0_INTP_io( :, :, : )

    !=====================================================================


    REAL(WP), ALLOCATABLE :: Cf_LW_io(:),    Cf_UW_io(:),     Cf_ave_io(:)
    REAL(WP), ALLOCATABLE :: U_tau_LW_io(:), Re_tau_LW_io(:), U_tau_ave_io(:)
    REAL(WP), ALLOCATABLE :: U_tau_UW_io(:), Re_tau_UW_io(:), Re_tau_ave_io(:)

    !========================= ThermaL ===================================
    REAL(WP), ALLOCATABLE :: MDOt(:)
    REAL(WP), ALLOCATABLE :: Gbuk(:)
    REAL(WP), ALLOCATABLE :: DHDOt(:)
    REAL(WP), ALLOCATABLE :: Hbuk(:)
    REAL(WP), ALLOCATABLE :: DHbuk(:)

    REAL(WP), ALLOCATABLE :: Tbuk(:)
    REAL(WP), ALLOCATABLE :: Dbuk(:)
    REAL(WP), ALLOCATABLE :: Mbuk(:)
    REAL(WP), ALLOCATABLE :: Kbuk(:)
    REAL(WP), ALLOCATABLE :: Bbuk(:)
    REAL(WP), ALLOCATABLE :: Cpbk(:)

    REAL(WP), ALLOCATABLE :: DHwalz_RA(:, :)
    REAL(WP), ALLOCATABLE :: Hwalz_RA(:, :)
    REAL(WP), ALLOCATABLE :: Hwalz_FA(:, :)
    REAL(WP), ALLOCATABLE :: Twalz(:, :)
    REAL(WP), ALLOCATABLE :: Dwalz(:, :)
    REAL(WP), ALLOCATABLE :: Mwalz(:, :)
    REAL(WP), ALLOCATABLE :: Kwalz(:, :)
    REAL(WP), ALLOCATABLE :: Cpwalz(:, :)

    REAL(WP), ALLOCATABLE :: Ubuk(:)
    REAL(WP), ALLOCATABLE :: Rebk(:)
    REAL(WP), ALLOCATABLE :: Prbk(:)
    REAL(WP), ALLOCATABLE :: Nubk(:, :)
    REAL(WP), ALLOCATABLE :: Bobk(:)
    REAL(WP), ALLOCATABLE :: Grbk(:, :)

    REAL(WP), ALLOCATABLE :: qw_D(:, :)
    REAL(WP), ALLOCATABLE :: hc_D(:, :)
    REAL(WP), ALLOCATABLE :: Hwalz_RA_D(:, :)
    REAL(WP), ALLOCATABLE :: Twalz_D(:, :)
    REAL(WP), ALLOCATABLE :: Dwalz_D(:, :)
    REAL(WP), ALLOCATABLE :: Mwalz_D(:, :)
    REAL(WP), ALLOCATABLE :: Kwalz_D(:, :)
    REAL(WP), ALLOCATABLE :: Cpwalz_D(:, :)

    REAL(WP), ALLOCATABLE :: Ubuk_D(:)
    REAL(WP), ALLOCATABLE :: Gbuk_D(:)
    REAL(WP), ALLOCATABLE :: Hbuk_D(:)
    REAL(WP), ALLOCATABLE :: Tbuk_D(:)
    REAL(WP), ALLOCATABLE :: Dbuk_D(:)
    REAL(WP), ALLOCATABLE :: Mbuk_D(:)
    REAL(WP), ALLOCATABLE :: Kbuk_D(:)
    REAL(WP), ALLOCATABLE :: Bbuk_D(:)
    REAL(WP), ALLOCATABLE :: Cpbk_D(:)


    REAL(WP), ALLOCATABLE :: Cfbk_U(:)
    REAL(WP), ALLOCATABLE :: Cfbk_L(:)
    REAL(WP), ALLOCATABLE :: Cfbk_A(:)

END MODULE


!**********************************************************************************************************************************
SUBROUTINE MEMO_ALLOCT_AVERAGE_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    ALLOCATE( U1ztL_F0_io( NCL1_io, NCL2, NDV + 1 ) )
    ALLOCATE( G1ztL_F0_io( NCL1_io, NCL2, NDV   ) )
    ALLOCATE( UPztL_F0_io( NCL1_io, NCL2, NDV   ) )

    ALLOCATE( U2ztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3 ) )
    ALLOCATE( UGztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3 ) )

    ALLOCATE( UGUztL_F0_io(NCL1_io, NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) )

    ALLOCATE( DVDL1ztL_F0_io( NCL1_io, NCL2, NDV, NDV  ) )
    ALLOCATE( DVDLPztL_F0_io( NCL1_io, NCL2, NDV, NDV  ) )
    ALLOCATE( DVDL2ztL_F0_io( NCL1_io, NCL2, NDV * (7 - NDV) / 2 + NDV - 3, NDV  ) )

    ALLOCATE( D1ztL_F0_io( NCL1_io, NCL2 ) ) ;  D1ztL_F0_io = 1.0_WP

    IF(iThermoDynamics == 1) THEN

        ALLOCATE( T1ztL_F0_io( NCL1_io, NCL2 ) )

        ALLOCATE( H1ztL_F0_io( NCL1_io, NCL2 ) )

        ALLOCATE( T2ztL_F0_io( NCL1_io, NCL2 ) )
        ALLOCATE( D2ztL_F0_io( NCL1_io, NCL2 ) )
        ALLOCATE( H2ztL_F0_io( NCL1_io, NCL2 ) )

        ALLOCATE( DHztL_F0_io( NCL1_io, NCL2 ) )

        ALLOCATE( UHztL_F0_io( NCL1_io, NCL2, NDV ) )
        ALLOCATE( GHztL_F0_io( NCL1_io, NCL2, NDV ) )

        ALLOCATE( MDOt(NCL1_io) )
        ALLOCATE( Gbuk(NCL1_io) )
        ALLOCATE( DHDOt(NCL1_io) )
        ALLOCATE( Hbuk(NCL1_io) )
        ALLOCATE( DHbuk(NCL1_io) )

        ALLOCATE( Tbuk(NCL1_io) )
        ALLOCATE( Dbuk(NCL1_io) )
        ALLOCATE( Mbuk(NCL1_io) )
        ALLOCATE( Kbuk(NCL1_io) )
        ALLOCATE(Bbuk(NCL1_io) )
        ALLOCATE( Cpbk(NCL1_io) )

        ALLOCATE( Hwalz_FA(NCL1_io, 2) )
        ALLOCATE( Hwalz_RA(NCL1_io, 2) )
        ALLOCATE( DHwalz_RA(NCL1_io, 2) )
        ALLOCATE( Twalz(NCL1_io, 2) )
        ALLOCATE( Dwalz(NCL1_io, 2) )
        ALLOCATE( Mwalz(NCL1_io, 2) )
        ALLOCATE( Kwalz(NCL1_io, 2) )
        ALLOCATE( Cpwalz(NCL1_io, 2) )

        ALLOCATE( qw_D(NCL1_io, 2) )
        ALLOCATE( hc_D(NCL1_io, 2) )
        ALLOCATE( Hwalz_RA_D(NCL1_io, 2) )
        ALLOCATE( Twalz_D(NCL1_io, 2) )
        ALLOCATE( Dwalz_D(NCL1_io, 2) )
        ALLOCATE( Mwalz_D(NCL1_io, 2) )
        ALLOCATE( Kwalz_D(NCL1_io, 2) )
        ALLOCATE( Cpwalz_D(NCL1_io, 2) )

        ALLOCATE( Gbuk_D(NCL1_io) )
        ALLOCATE( Ubuk_D(NCL1_io) )
        ALLOCATE( Hbuk_D(NCL1_io) )
        ALLOCATE( Tbuk_D(NCL1_io) )
        ALLOCATE( Dbuk_D(NCL1_io) )
        ALLOCATE( Mbuk_D(NCL1_io) )
        ALLOCATE( Kbuk_D(NCL1_io) )
        ALLOCATE(Bbuk_D(NCL1_io) )
        ALLOCATE( Cpbk_D(NCL1_io) )



        ALLOCATE( Ubuk(NCL1_io) )
        ALLOCATE( Rebk(NCL1_io) )
        ALLOCATE( Prbk(NCL1_io) )
        ALLOCATE( Grbk(NCL1_io, 2) )

        ALLOCATE( Nubk(NCL1_io, 2) )
        ALLOCATE(Bobk(NCL1_io) )

        ALLOCATE( Cfbk_A(NCL1_io) )
        ALLOCATE( Cfbk_L(NCL1_io) )
        ALLOCATE( Cfbk_U(NCL1_io) )
    END IF

    ALLOCATE( Cf_LW_io(NCL1_io) )
    ALLOCATE( Cf_UW_io(NCL1_io) )
    ALLOCATE( Cf_AVE_io(NCL1_io) )

    ALLOCATE( U_tau_LW_io(NCL1_io) )
    ALLOCATE( U_tau_UW_io(NCL1_io) )
    ALLOCATE( U_tau_AVE_io(NCL1_io) )

    ALLOCATE( Re_tau_LW_io(NCL1_io) )
    ALLOCATE( Re_tau_UW_io(NCL1_io) )
    ALLOCATE( Re_tau_AVE_io(NCL1_io) )


    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE MEMO_DEALLT_AVERAGE_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    IMPLICIT NONE


    DEALLOCATE( U1ztL_F0_io )
    DEALLOCATE( G1ztL_F0_io )
    DEALLOCATE( UPztL_F0_io )

    DEALLOCATE( U2ztL_F0_io )
    DEALLOCATE( UGztL_F0_io )

    DEALLOCATE( UGUztL_F0_io )

    DEALLOCATE( DVDL1ztL_F0_io )
    DEALLOCATE( DVDLPztL_F0_io )
    DEALLOCATE( DVDL2ztL_F0_io )

    DEALLOCATE( D1ztL_F0_io )

    IF(iThermoDynamics == 1) THEN

        DEALLOCATE( T1ztL_F0_io )

        DEALLOCATE( H1ztL_F0_io )

        DEALLOCATE( T2ztL_F0_io )
        DEALLOCATE( D2ztL_F0_io )
        DEALLOCATE( H2ztL_F0_io )

        DEALLOCATE( DHztL_F0_io )

        DEALLOCATE( UHztL_F0_io )
        DEALLOCATE( GHztL_F0_io )

        DEALLOCATE( MDOt )
        DEALLOCATE( Gbuk )
        DEALLOCATE( DHDOt )
        DEALLOCATE( Hbuk )
        DEALLOCATE( DHbuk )

        DEALLOCATE( Tbuk )
        DEALLOCATE( Dbuk )
        DEALLOCATE( Mbuk )
        DEALLOCATE( Kbuk )
        DEALLOCATE(Bbuk )
        DEALLOCATE( Cpbk )

        DEALLOCATE( Hwalz_FA )
        DEALLOCATE( Hwalz_RA )
        DEALLOCATE( DHwalz_RA )
        DEALLOCATE( Twalz )
        DEALLOCATE( Dwalz )
        DEALLOCATE( Mwalz )
        DEALLOCATE( Kwalz )
        DEALLOCATE( Cpwalz )

        DEALLOCATE( qw_D )
        DEALLOCATE( hc_D )
        DEALLOCATE( Hwalz_RA_D )
        DEALLOCATE( Twalz_D )
        DEALLOCATE( Dwalz_D )
        DEALLOCATE( Mwalz_D )
        DEALLOCATE( Kwalz_D )
        DEALLOCATE( Cpwalz_D )

        DEALLOCATE( Gbuk_D )
        DEALLOCATE( Ubuk_D )
        DEALLOCATE( Hbuk_D )
        DEALLOCATE( Tbuk_D )
        DEALLOCATE( Dbuk_D )
        DEALLOCATE( Mbuk_D )
        DEALLOCATE( Kbuk_D )
        DEALLOCATE(Bbuk_D )
        DEALLOCATE( Cpbk_D )



        DEALLOCATE( Ubuk )
        DEALLOCATE( Rebk )
        DEALLOCATE( Prbk )
        DEALLOCATE( Grbk )

        DEALLOCATE( Nubk )
        DEALLOCATE(Bobk )

        DEALLOCATE( Cfbk_A )
        DEALLOCATE( Cfbk_L )
        DEALLOCATE( Cfbk_U )
    END IF

    DEALLOCATE( Cf_LW_io )
    DEALLOCATE( Cf_UW_io )
    DEALLOCATE( Cf_AVE_io )

    DEALLOCATE( U_tau_LW_io )
    DEALLOCATE( U_tau_UW_io )
    DEALLOCATE( U_tau_AVE_io )

    DEALLOCATE( Re_tau_LW_io )
    DEALLOCATE( Re_tau_UW_io )
    DEALLOCATE( Re_tau_AVE_io )

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE MEMO_ALLOCT_INTP_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    ALLOCATE( U1ztL_F0_INTP_io( NND1_io, NND2, NDV + 1 ) )
    ALLOCATE( G1ztL_F0_INTP_io( NND1_io, NND2, NDV   ) )
    ALLOCATE( UPztL_F0_INTP_io( NND1_io, NND2, NDV  ) )

    ALLOCATE( U2ztL_F0_INTP_io( NND1_io, NND2, NDV * (7 - NDV) / 2 + NDV - 3 ) )
    ALLOCATE( UGztL_F0_INTP_io( NND1_io, NND2, NDV * (7 - NDV) / 2 + NDV - 3 ) )

    ALLOCATE( UGUztL_F0_INTP_io(NND1_io, NND2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) )

    ALLOCATE( D1ztL_F0_INTP_io( NND1_io, NND2 ) )  ; D1ztL_F0_INTP_io = 1.0_WP

    IF(iThermoDynamics == 1) THEN

        ALLOCATE( T1ztL_F0_INTP_io( NND1_io, NND2 ) )

        ALLOCATE( H1ztL_F0_INTP_io( NND1_io, NND2 ) )

        ALLOCATE( T2ztL_F0_INTP_io( NND1_io, NND2 ) )
        ALLOCATE( D2ztL_F0_INTP_io( NND1_io, NND2 ) )
        ALLOCATE( H2ztL_F0_INTP_io( NND1_io, NND2 ) )

        ALLOCATE( DHztL_F0_INTP_io( NND1_io, NND2 ) )

        ALLOCATE( UHztL_F0_INTP_io( NND1_io, NND2, NDV ) )
        ALLOCATE( GHztL_F0_INTP_io( NND1_io, NND2, NDV ) )
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE MEMO_DEALLT_INTP_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    DEALLOCATE( U1ztL_F0_INTP_io )
    DEALLOCATE( G1ztL_F0_INTP_io )
    DEALLOCATE( UPztL_F0_INTP_io )

    DEALLOCATE( U2ztL_F0_INTP_io )
    DEALLOCATE( UGztL_F0_INTP_io )

    DEALLOCATE( UGUztL_F0_INTP_io )

    DEALLOCATE( D1ztL_F0_INTP_io )

    IF(iThermoDynamics == 1) THEN
        DEALLOCATE( T1ztL_F0_INTP_io )


        DEALLOCATE( H1ztL_F0_INTP_io )

        DEALLOCATE( T2ztL_F0_INTP_io )
        DEALLOCATE( D2ztL_F0_INTP_io )
        DEALLOCATE( H2ztL_F0_INTP_io )

        DEALLOCATE( DHztL_F0_INTP_io )

        DEALLOCATE( UHztL_F0_INTP_io )
        DEALLOCATE( GHztL_F0_INTP_io )
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: INN
    INTEGER(4) :: J, JJ, I
    INTEGER(4) :: L, IP, M
    INTEGER(4) :: N2DOID

    CHARACTER(15) :: NXSTR
    REAL(WP) :: urms, vrms, WRms, uv, uw, vw
    REAL(WP) :: COE1, COE2,COE
    INTEGER(4) :: TECFLG1, TECFLG2, TECFLG3

    REAL(WP) :: D1AUX (NCL1_io, N2DO(MYID), NDV + 1,                               1 : NPTOT)
    REAL(WP) :: D2AUX (NCL1_io, N2DO(MYID), NDV,                                 1 : NPTOT)
    REAL(WP) :: D3AUX (NCL1_io, N2DO(MYID), NDV,                                 1 : NPTOT)
    REAL(WP) :: D4AUX (NCL1_io, N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3),                1 : NPTOT)
    REAL(WP) :: D5AUX (NCL1_io, N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3),                1 : NPTOT)
    REAL(WP) :: D6AUX (NCL1_io, N2DO(MYID), (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8),  1 : NPTOT)

    REAL(WP) :: T1AUX (NCL1_io, N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: T2AUX (NCL1_io, N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: T3AUX (NCL1_io, N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3), NDV, 1 : NPTOT)

    REAL(WP) :: H1AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H2AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H3AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H4AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H5AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H6AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H7AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H8AUX (NCL1_io, N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: H9AUX (NCL1_io, N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: H10AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H11AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H12AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)
    REAL(WP) :: H13AUX (NCL1_io, N2DO(MYID),     1 : NPTOT)


    INN = NCL1_io * N2DO(MYID) * (NDV + 1)
    CALL MPI_GATHER( U1ztL_io, INN, MPI_DOUBLE_PRECISION, D1AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * NDV
    CALL MPI_GATHER( G1ztL_io, INN, MPI_DOUBLE_PRECISION, D2AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * NDV
    CALL MPI_GATHER( UPztL_io, INN, MPI_DOUBLE_PRECISION, D3AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
    CALL MPI_GATHER( U2ztL_io, INN, MPI_DOUBLE_PRECISION, D4AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
    CALL MPI_GATHER( UGztL_io, INN, MPI_DOUBLE_PRECISION, D5AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
    CALL MPI_GATHER( UGUztL_io, INN, MPI_DOUBLE_PRECISION, D6AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDL1ztL_io, INN, MPI_DOUBLE_PRECISION, T1AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDLPztL_io, INN, MPI_DOUBLE_PRECISION, T2AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = NCL1_io * N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3) * NDV
    CALL MPI_GATHER( DVDL2ztL_io, INN, MPI_DOUBLE_PRECISION, T3AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    IF(iThermoDynamics == 1) THEN

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( T1ztL_io, INN, MPI_DOUBLE_PRECISION, H1AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( D1ztL_io, INN, MPI_DOUBLE_PRECISION, H2AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( H1ztL_io, INN, MPI_DOUBLE_PRECISION, H3AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( T2ztL_io, INN, MPI_DOUBLE_PRECISION, H4AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( D2ztL_io, INN, MPI_DOUBLE_PRECISION, H5AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( H2ztL_io, INN, MPI_DOUBLE_PRECISION, H6AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( DHztL_io, INN, MPI_DOUBLE_PRECISION, H7AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID) * NDV
        CALL MPI_GATHER( UHztL_io, INN, MPI_DOUBLE_PRECISION, H8AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID) * NDV
        CALL MPI_GATHER( GHztL_io, INN, MPI_DOUBLE_PRECISION, H9AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)


        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( K1ztL_io, INN, MPI_DOUBLE_PRECISION, H10AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( K2ztL_io, INN, MPI_DOUBLE_PRECISION, H11AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( M1ztL_io, INN, MPI_DOUBLE_PRECISION, H12AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = NCL1_io * N2DO(MYID)
        CALL MPI_GATHER( M2ztL_io, INN, MPI_DOUBLE_PRECISION, H13AUX, INN,  &
        MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

    END IF


    IF(MYID == 0) THEN
        !==================================================================
        CALL MEMO_ALLOCT_AVERAGE_nonXperiodic_io

        DO IP = 0, NPSLV
            N2DOID=JDEWT(IP) - JDSWT(IP) + 1
            DO J = 1, N2DOID
                JJ = JDSWT(IP) - 1 + J
                DO I = 1, NCL1_io

                    DO L = 1, NDV + 1
                        U1ztL_F0_io(I, JJ, L) = D1AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, NDV
                        G1ztL_F0_io(I, JJ, L) = D2AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, NDV
                        UPztL_F0_io(I, JJ, L) = D3AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                        U2ztL_F0_io(I, JJ, L) = D4AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                        UGztL_F0_io(I, JJ, L) = D5AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
                        UGUztL_F0_io(I, JJ, L) = D6AUX(I, J, L, IP + 1)
                    END DO

                    DO L = 1, NDV
                        DO M = 1, NDV
                            DVDL1ztL_F0_io(I, JJ, L, M) = T1AUX(I, J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, NDV
                        DO M = 1, NDV
                            DVDLPztL_F0_io(I, JJ, L, M) = T2AUX(I, J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                        DO M = 1, NDV
                            DVDL2ztL_F0_io(I, JJ, L, M) = T3AUX(I, J, L, M, IP + 1)
                        END DO
                    END DO

                    IF(iThermoDynamics == 1) THEN

                        T1ztL_F0_io(I, JJ) = H1AUX(I, J, IP + 1)
                        D1ztL_F0_io(I, JJ) = H2AUX(I, J, IP + 1)
                        H1ztL_F0_io(I, JJ) = H3AUX(I, J, IP + 1)

                        T2ztL_F0_io(I, JJ) = H4AUX(I, J, IP + 1)
                        D2ztL_F0_io(I, JJ) = H5AUX(I, J, IP + 1)
                        H2ztL_F0_io(I, JJ) = H6AUX(I, J, IP + 1)

                        DHztL_F0_io(I, JJ) = H7AUX(I, J, IP + 1)

                        DO L = 1, NDV
                            UHztL_F0_io(I, JJ, L) = H8AUX(I, J, L, IP + 1)
                        END DO

                        DO L = 1, NDV
                            GHztL_F0_io(I, JJ, L) = H9AUX(I, J, L, IP + 1)
                        END DO
                    END IF

                END DO
            END DO
        END DO

        CALL MEMO_ALLOCT_INTP_nonXperiodic_io
        CALL CL2ND_INTP_AVEG_nonXperiodic_io

        IF(iThermoDynamics == 1) THEN
            CALL WRT_FARve_Average_Contour_nonXperiodic_io
            CALL WRT_HeatTransfer_XProfile_nonXperiodic_io
        END IF

        CALL WRT_Cf_XProfile_nonXperiodic_io
        CALL WRT_Reynolds_Average_Contour_nonXperiodic_io

        CALL MEMO_DEALLT_INTP_nonXperiodic_io
        CALL MEMO_DEALLT_AVERAGE_nonXperiodic_io
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_Reynolds_Average_Contour_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CHARACTER(128) :: FLNM
    REAL(WP) :: ux, uy, uz, p, &
    tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
    den, tem, enh, drms, trms, hrms, Thfx, thfy, thfz
    INTEGER(4) :: I, J, TECFLG_ReyAG

    TECFLG_ReyAG = 200
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath2) // 'Result.IO.Reynd.Averaged.Flow.' // TRIM(PNTIM) // '.StAR.Contour.plt'

    !==================== Reynolds Avderaged ContourS ================================================================
    OPEN(TECFLG_ReyAG, FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_ReyAG, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    IF(iThermoDynamics == 0) THEN
        WRITE(TECFLG_ReyAG, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
        ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw"'
    ELSE IF(iThermoDynamics == 1) THEN
        WRITE(TECFLG_ReyAG, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
        ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
        ' "D", "T", "H", "Drms", "Trms", "Hrms", ', &
        ' "thfx", "thfy", "thfz" '
    ELSE
    END IF

    WRITE(TECFLG_ReyAG, '(A, 1ES13.5, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', &
    PhyTIME, ' ", I = ', NND1_io, ', J = ', NND2, ', F =POINT'

    DO J = 1, NND2
        DO I = 1, NND1_io

            ux  = U1ztL_F0_INTP_io(I, J, 1)
            uy  = U1ztL_F0_INTP_io(I, J, 2)
            uz  = U1ztL_F0_INTP_io(I, J, 3)

            p = U1ztL_F0_INTP_io(I, J, 4)

            Ruu = ( U2ztL_F0_INTP_io(I, J, 1) - U1ztL_F0_INTP_io(I, J, 1) * U1ztL_F0_INTP_io(I, J, 1) ) * D1ztL_F0_INTP_io(I, J)
            Ruv = ( U2ztL_F0_INTP_io(I, J, 2) - U1ztL_F0_INTP_io(I, J, 1) * U1ztL_F0_INTP_io(I, J, 2) ) * D1ztL_F0_INTP_io(I, J)
            Ruw = ( U2ztL_F0_INTP_io(I, J, 3) - U1ztL_F0_INTP_io(I, J, 1) * U1ztL_F0_INTP_io(I, J, 3) ) * D1ztL_F0_INTP_io(I, J)
            Rvv = ( U2ztL_F0_INTP_io(I, J, 4) - U1ztL_F0_INTP_io(I, J, 2) * U1ztL_F0_INTP_io(I, J, 2) ) * D1ztL_F0_INTP_io(I, J)
            Rvw = ( U2ztL_F0_INTP_io(I, J, 5) - U1ztL_F0_INTP_io(I, J, 2) * U1ztL_F0_INTP_io(I, J, 3) ) * D1ztL_F0_INTP_io(I, J)
            Rww = ( U2ztL_F0_INTP_io(I, J,6) - U1ztL_F0_INTP_io(I, J, 3) * U1ztL_F0_INTP_io(I, J, 3) ) * D1ztL_F0_INTP_io(I, J)

            tke = 0.5_WP * (Ruu+ Rvv+ Rww)

            IF(iThermoDynamics == 1) THEN
                den = D1ztL_F0_INTP_io(I, J)
                tem = T1ztL_F0_INTP_io(I, J)
                enh = H1ztL_F0_INTP_io(I, J)

                drms = DSQRT(DABS( D2ztL_F0_INTP_io(I, J) - D1ztL_F0_INTP_io(I, J) * D1ztL_F0_INTP_io(I, J) ) )
                trms = DSQRT(DABS( T2ztL_F0_INTP_io(I, J) - T1ztL_F0_INTP_io(I, J) * T1ztL_F0_INTP_io(I, J) ) )
                hrms = DSQRT(DABS( H2ztL_F0_INTP_io(I, J) - H1ztL_F0_INTP_io(I, J) * H1ztL_F0_INTP_io(I, J) ) )

                thfx = ( UHztL_F0_INTP_io(I, J, 1) - U1ztL_F0_INTP_io(I, J, 1) * H1ztL_F0_INTP_io(I, J) ) * D1ztL_F0_INTP_io(I, J)
                thfy = ( UHztL_F0_INTP_io(I, J, 2) - U1ztL_F0_INTP_io(I, J, 2) * H1ztL_F0_INTP_io(I, J) ) * D1ztL_F0_INTP_io(I, J)
                thfz = ( UHztL_F0_INTP_io(I, J, 3) - U1ztL_F0_INTP_io(I, J, 3) * H1ztL_F0_INTP_io(I, J) ) * D1ztL_F0_INTP_io(I, J)
            END IF

            IF(iThermoDynamics == 0) THEN
                WRITE(TECFLG_ReyAG, '(13ES17.9)') XND_io(I), YND(J), ux, uy, uz, p, &
                tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw
            ELSE IF(iThermoDynamics == 1) THEN
                WRITE(TECFLG_ReyAG, '(22ES17.9)') XND_io(I), YND(J), ux, uy, uz, p, &
                tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
                den, tem, enh, drms, trms, hrms, &
                thfx, thfy, thfz
            ELSE
            END IF

        END DO
    END DO

    CLOSE(TECFLG_ReyAG)

    CALL POSITIONS_nonXperiodic_io

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_FARve_Average_Contour_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CHARACTER(128) :: FLNM
    REAL(WP) :: ux, uy, uz, p, &
    tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
    den, tem, enh, drms, trms, hrms, thfx, thfy, thfz
    INTEGER(4) :: I, J, TECFLG_FavAG

    IF(iThermoDynamics /= 1) RETURN

    TECFLG_FavAG = 200

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath2) // 'Result.IO.Favre.Averaged.Flow.' // TRIM(PNTIM) // '.StAR.Contour.plt'

    !====================FARve Avderaged ContourS ================================================================
    OPEN(TECFLG_FavAG, FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_FavAG, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    WRITE(TECFLG_FavAG, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
    ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
    ' "D", "T", "H", "Drms", "Trms", "Hrms", ', &
    ' "thfx", "thfy", "thfz" '
    WRITE(TECFLG_FavAG, '(A, 1ES13.5, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', &
    PhyTIME, ' ", I = ', NND1_io, ', J = ', NND2, ', F =POINT'

    DO J = 1, NND2
        DO I = 1, NND1_io

            ux  = G1ztL_F0_INTP_io(I, J, 1) / D1ztL_F0_INTP_io(I, J)
            uy  = G1ztL_F0_INTP_io(I, J, 2) / D1ztL_F0_INTP_io(I, J)
            uz  = G1ztL_F0_INTP_io(I, J, 3) / D1ztL_F0_INTP_io(I, J)

            p = U1ztL_F0_INTP_io(I, J, 4) ! same

            Ruu = UGztL_F0_INTP_io(I, J, 1) - G1ztL_F0_INTP_io(I, J, 1) * G1ztL_F0_INTP_io(I, J, 1) / D1ztL_F0_INTP_io(I, J)
            Ruv = UGztL_F0_INTP_io(I, J, 2) - G1ztL_F0_INTP_io(I, J, 1) * G1ztL_F0_INTP_io(I, J, 2) / D1ztL_F0_INTP_io(I, J)
            Ruw = UGztL_F0_INTP_io(I, J, 3) - G1ztL_F0_INTP_io(I, J, 1) * G1ztL_F0_INTP_io(I, J, 3) / D1ztL_F0_INTP_io(I, J)
            Rvv = UGztL_F0_INTP_io(I, J, 4) - G1ztL_F0_INTP_io(I, J, 2) * G1ztL_F0_INTP_io(I, J, 2) / D1ztL_F0_INTP_io(I, J)
            Rvw = UGztL_F0_INTP_io(I, J, 5) - G1ztL_F0_INTP_io(I, J, 2) * G1ztL_F0_INTP_io(I, J, 3) / D1ztL_F0_INTP_io(I, J)
            Rww = UGztL_F0_INTP_io(I, J,6) - G1ztL_F0_INTP_io(I, J, 3) * G1ztL_F0_INTP_io(I, J, 3) / D1ztL_F0_INTP_io(I, J)

            tke = 0.5_WP * (Ruu+ Rvv+ Rww)

            den = D1ztL_F0_INTP_io(I, J) !same
            tem = T1ztL_F0_INTP_io(I, J) !same
            enh = DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J)

            drms = DSQRT(DABS( D2ztL_F0_INTP_io(I, J) - D1ztL_F0_INTP_io(I, J) * D1ztL_F0_INTP_io(I, J) ) ) ! same
            trms = DSQRT(DABS( T2ztL_F0_INTP_io(I, J) - T1ztL_F0_INTP_io(I, J) * T1ztL_F0_INTP_io(I, J) ) ) ! same
            hrms = DSQRT(DABS( H2ztL_F0_INTP_io(I, J) - DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J) * &
            (2.0_WP * H1ztL_F0_INTP_io(I, J) - DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J) ) ) )

            thfx = GHztL_F0_INTP_io(I, J, 1) - G1ztL_F0_INTP_io(I, J, 1) * DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J)
            thfy = GHztL_F0_INTP_io(I, J, 2) - G1ztL_F0_INTP_io(I, J, 2) * DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J)
            thfz = GHztL_F0_INTP_io(I, J, 3) - G1ztL_F0_INTP_io(I, J, 3) * DHztL_F0_INTP_io(I, J) / D1ztL_F0_INTP_io(I, J)

            WRITE(TECFLG_FavAG, '(22ES17.9)') XND_io(I), YND(J), ux, uy, uz, p, &
            tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
            den, tem, enh, drms, trms, hrms, &
            thfx, thfy, thfz
        END DO
    END DO
    CLOSE(TECFLG_FavAG)


    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE WRT_HeatTransfer_XProfile_nonXperiodic_io
    !=================================================================
    !====== Defination based on Joong Hun Bae, 2005===========
    !=== MDOt = mass flow rate
    !=== Mbuk = bulk mass flux
    !=== HDOt = enthalpy flow rate
    !=== Hbuk = bulk enthalpy
    !==
    !=== Tbuk = bulk temperature
    !=== Dbuk = bulk DENSITY
    !=== Mbuk = bulk Viscousity
    !=== Kbuk = bulk thermal conductivity
    !===Cpbk = bulk Cp
    !===
    !===Ubuk = bulk velocity
    !=== Rebk = bulk Reynolds no.
    !===Prbk = bulk Prandtl no.
    !=== Nubk = bulk Nusselt number
    !===Bobk = Boyancy PARAMETER
    !===Cfbk = local Friction coefficient

    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CHARACTER(128) :: FLNM
    INTEGER(4) :: I, J, TECFLG_PROF, IP
    REAL(WP) :: Hw, Tw, Dw, Cpw, Mw, Kw, AREA, Tw1, Tw2, AREA_io
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: B_tmp


    !=============wall PARAMETERS ==undiM ==================================
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature ) THEN
        DO I = 1, NCL1_io
            Hwalz_RA(I, iTopWall) = H_WAL_GV(NCL1_io / 2, iTopWall)
            Dwalz(I, iTopWall) = D_WAL_GV(NCL1_io / 2, iTopWall)
            Twalz(I, iTopWall) = T_WAL_GV(NCL1_io / 2, iTopWall)
            Mwalz(I, iTopWall) = M_WAL_GV(NCL1_io / 2, iTopWall)
            Kwalz(I, iTopWall) = K_WAL_GV(NCL1_io / 2, iTopWall)
            Cpwalz(I, iTopWall) = Cp_WAL_GV(NCL1_io / 2, iTopWall)
            Hwalz_FA(I, iTopWall) = Hwalz_RA(I, iTopWall)
        END DO
    END IF

    IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature ) THEN
        DO I = 1, NCL1_io
            Hwalz_RA(I, iBotWall) = H_WAL_GV(NCL1_io / 2, iBotWall)
            Dwalz(I, iBotWall) = D_WAL_GV(NCL1_io / 2, iBotWall)
            Twalz(I, iBotWall) = T_WAL_GV(NCL1_io / 2, iBotWall)
            Mwalz(I, iBotWall) = M_WAL_GV(NCL1_io / 2, iBotWall)
            Kwalz(I, iBotWall) = K_WAL_GV(NCL1_io / 2, iBotWall)
            Cpwalz(I, iBotWall) = Cp_WAL_GV(NCL1_io / 2, iBotWall)
            Hwalz_FA(I, iBotWall) = Hwalz_RA(I, iTopWall)
        END DO
    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux ) THEN
        J = NND2
        DO I = 1, NCL1_io
            IP = I + 1
            DHwalz_RA(I, iTopWall) = 0.5_WP * ( DHztL_F0_INTP_io(I, J) + DHztL_F0_INTP_io(IP, J) )
            Hwalz_FA(I, iTopWall) = (0.5_WP * ( DHztL_F0_INTP_io(I, J) + DHztL_F0_INTP_io(IP, J) )) / &
            0.5_WP * ( D1ztL_F0_INTP_io(I, J) + D1ztL_F0_INTP_io(IP, J) )
            CALL THERM_PROP_UPDATE_FROM_DH(DHwalz_RA(I, iTopWall), Hwalz_RA(I, iTopWall), Twalz(I, iTopWall), &
            Dwalz(I, iTopWall), Mwalz(I, iTopWall), Kwalz(I, iTopWall), Cpwalz(I, iTopWall), B_tmp)

        END DO
    END IF

    IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux ) THEN
        J = 1
        DO I = 1, NCL1_io
            IP = I + 1
            DHwalz_RA(I, iBotWall) = 0.5_WP * ( DHztL_F0_INTP_io(I, J) + DHztL_F0_INTP_io(IP, J) )
            Hwalz_FA(I, iBotWall) = (0.5_WP * ( DHztL_F0_INTP_io(I, J) + DHztL_F0_INTP_io(IP, J) )) / &
            0.5_WP * ( D1ztL_F0_INTP_io(I, J) + D1ztL_F0_INTP_io(IP, J) )
            CALL THERM_PROP_UPDATE_FROM_DH(DHwalz_RA(I, iBotWall), Hwalz_RA(I, iBotWall), Twalz(I, iBotWall), &
            Dwalz(I, iBotWall), Mwalz(I, iBotWall), Kwalz(I, iBotWall), Cpwalz(I, iBotWall), B_tmp)
        END DO

    END IF

    Hwalz_RA_D(:, :) = Hwalz_RA(:, :) * T0 * CP0 + H0
    Twalz_D(:, :) = Twalz(:, :) * T0
    Dwalz_D(:, :) = Dwalz(:, :) * D0
    Mwalz_D(:, :) = Mwalz(:, :) * M0
    Kwalz_D(:, :) = Kwalz(:, :) * K0
    Cpwalz_D(:, :) = Cpwalz(:, :) * CP0

    !=======================bulk
    AREA_io = 0.0_WP
    DO J = 1, NCL2
        AREA_io = AREA_io + 1.0_WP / RCCI1(J) / DYFI(J)
    END DO

    !===============bulk valueS =======================================================================
    !==================calculate bulk mass flux, and bulk enthalpY =========
    DO I = 1, NCL1_io
        MDOt(I) = 0.0_WP
        DHDOt(I) = 0.0_WP
        DO J = 1, NCL2
            MDOt(I) = MDOt(I) + G1ztL_F0_io(I, J, 1) / RCCI1(J) / DYFI(J)
            DHDOt(I) = DHDOt(I) + DHztL_F0_io(I, J) / RCCI1(J) / DYFI(J)
        END DO
        Gbuk(I) = MDOt(I) /AREA_io
        DHbuk(I) = DHDOt(I) /AREA_io
        Gbuk_D(I) = Gbuk(I) * D0 * U0
        !Hbuk_D(I) = Hbuk(I) * T0 * CP0 + H0
    END DO

    DO I = 1, NCL1_io

        CALL THERM_PROP_UPDATE_FROM_DH(DHbuk(I), Hbuk(I), Tbuk(I), Dbuk(I), Mbuk(I), Kbuk(I), Cpbk(I), Bbuk(I))

        Dbuk_D(I) = Dbuk(I) * D0
        Tbuk_D(I) = Tbuk(I) * T0
        Mbuk_D(I) = Mbuk(I) * M0
        Kbuk_D(I) = Kbuk(I) * K0
        Cpbk_D(I) = CPbk(I) * CP0
        Bbuk_D(I) = Bbuk(I) * B0
        Hbuk(I) = Hbuk(I) * T0 * CP0 + H0

    END DO

    !====================================================================
    DO I = 1, NCL1_io
        Ubuk(I) = Gbuk(I) / Dbuk(I)
        Ubuk_D(I) = Ubuk(I) * U0

        Rebk(I) = Dbuk(I) * Ubuk(I) * (HYT - HYB) / Mbuk(I)
        Prbk(I) = Mbuk(I) * Cpbk(I) / Kbuk(I) * Prt0

    END DO


    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
        DO I = 1, NCL1_io
            qw_D(I, 1:2) = thermalWallBC_Dim(1:2)
            !====convective heat transfer ratE ==== HC = Qw/ (Tw- Tb) ====
            hc_D(I, 1:2) = qw_D(I, 1:2) / (Twalz_D(I, 1:2) - Tbuk_D(I))
            !====bulk Grashof numbeR ===========

            Grbk(I, 1) = G_A*Bbuk_D(I) * Qw_D(I, 1) * ((2.0_WP * L0)**4) * Dbuk_D(I) * Dbuk_D(I) / Mbuk_D(I) / Mbuk_D(I) / Kbuk_D(I)
            Grbk(I, 2) = G_A*Bbuk_D(I) * Qw_D(I, 2) * ((2.0_WP * L0)**4) * Dbuk_D(I) * Dbuk_D(I) / Mbuk_D(I) / Mbuk_D(I) / Kbuk_D(I)
        END DO
    END IF
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        !====wall heat fluX ========== here IS not corret, as bAR over (kdT/ Dy) /= kdbART/ Dy
        DO I = 1, NCL1_io
            qw_D(I, 1) = Kwalz_D(I, 1) * (T1ztL_F0_io(I, 1) * T0- Twalz_D(I, 1)) / ((YCC(1) - YND(1)) * L0)
            qw_D(I, 2) = Kwalz_D(I, 2) * (Twalz_D(I, 2) - T1ztL_F0_io(I, NCL2) * T0) / ((YND(NND2) - YCC(NCL2)) * L0)
            !====convective heat transfer ratE ==== HC = Qw/ (Tw- Tb) ====
            hc_D(I, 1:2) = qw_D(I, 1:2) / (Twalz_D(I, 1:2) - Tbuk_D(I))
            !====bulk Grashof numbeR ===========
            Grbk(I, 1) = G_A*Bbuk_D(I) * DABS(Twalz_D(I, 1) - &
            Twalz_D(I, 2)) * ((2.0_WP * L0)**3) * Dbuk_D(I) * Dbuk_D(I) / Mbuk_D(I) / Mbuk_D(I)
        END DO
    END IF
    !Added by Junjie Yan, 2017/02 /16, for pipe flow.
    !OtherWISe qw_D, hc_D and Grbk won't be calculated and THEN Bobk and Nubk will be 0
    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == 0) THEN
        DO I = 1, NCL1_io
            qw_D(I, 1:2) = thermalWallBC_Dim(1:2)
            !====convective heat transfer ratE ==== HC = Qw/ (Tw- Tb) ====
            hc_D(I, 1:2) = qw_D(I, 1:2) / (Twalz_D(I, 1:2) - Tbuk_D(I))
            !====bulk Grashof numbeR ===========

            Grbk(I, 1) = G_A*Bbuk_D(I) * Qw_D(I, 1) * ((2.0_WP * L0)**4) * Dbuk_D(I) * Dbuk_D(I) / Mbuk_D(I) / Mbuk_D(I) / Kbuk_D(I)
            Grbk(I, 2) = G_A*Bbuk_D(I) * Qw_D(I, 2) * ((2.0_WP * L0)**4) * Dbuk_D(I) * Dbuk_D(I) / Mbuk_D(I) / Mbuk_D(I) / Kbuk_D(I)
        END DO
    END IF

    !=============== Nusselt numbeR ==============
    DO I = 1, NCL1_io
        Nubk(I, 1) = 2.0_WP * L0 * Hc_D(I, 1) / Kbuk_D (I)
        Nubk(I, 2) = 2.0_WP * L0 * Hc_D(I, 2) / Kbuk_D (I)
        Bobk(I) = Grbk(I, 2) / (Rebk(I)**3.425_WP) / (Prbk(I)**0.8_WP)
    END DO


    !===============WRITE data ouT =============================================================
    TECFLG_PROF = 200
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath2) // 'Result.IO.X.DIStributed.HeatTransfer.' // TRIM(PNTIM) // '.StAR.Profile.plt'

    !==================== Reynolds Avderaged ContourS ================================================================
    OPEN(TECFLG_PROF, FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_PROF, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    WRITE(TECFLG_PROF, '(A, A, A)') 'variables = "1X", "2Gb", "3Ub", "4Hb", "5Db", "6Tb", "7Mb", "8Kb", ', &
    '"9Cpb", "10Reb", "11Prb", "12Grb", "13Bob", "14Nub", ', &
    '"15Hw", "16Tw", "17Dw", "18Cpw", "19Mw", "20Kw"'
    WRITE(TECFLG_PROF, '(A, 1ES13.5, A)') 'ZONE T = " ',PhyTIME, ' "'
    DO I = 1, NCL1_io

        Tw  = Twalz(I, 2)   !            TW  = 0.5_WP * ( Twalz(I, 1) + Twalz(I, 2) )
        Hw  = Hwalz_RA(I, 2)!            Hw  = 0.5_WP * ( Hwalz_RA(I, 1) + Hwalz_RA(I, 2) )
        Dw  = Dwalz(I, 2)   !            Dw  = 0.5_WP * ( Dwalz(I, 1) + Dwalz(I, 2) )
        Mw  = Mwalz(I, 2)   !            Mw  = 0.5_WP * ( Mwalz(I, 1) +Mwalz(I, 2) )
        Kw  = Kwalz(I, 2)   !            Kw  = 0.5_WP * ( Kwalz(I, 1) + Kwalz(I, 2) )
        Cpw = Cpwalz(I, 2)  !            Cpw  = 0.5_WP * ( Cpwalz(I, 1) +Cpwalz(I, 2) )



        WRITE(TECFLG_PROF, '(18ES17.9)') XCC_io(I), &
        Gbuk(I), Ubuk(I), Hbuk(I), Dbuk(I), Tbuk(I), Mbuk(I), &
        Kbuk(I), Cpbk(I), Rebk(I), Prbk(I), Grbk(I, 2), Bobk(I), Nubk(I, 2), &
        Hw, Tw, Dw, Cpw, Mw, Kw

    END DO
    close(TECFLG_PROF)


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_Cf_XProfile_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: TECFLG = 200
    INTEGER(4) :: I
    REAL(WP) :: DUDYL, DUDYU


    DO I = 1, NCL1_io

        ! method 1 and 2 give slightly dIFfeRENt results.
        ! method 2
        !DUDYL = DVDL1ztL_F0_io(I, 1,   1, 2)
        !DUDYU = DVDL1ztL_F0_io(I, NCL2, 1, 2)
        DUDYU = -(U1ztL_F0_io(I, NCL2, 1) - 0.0_WP) / (YCC(NCL2) - YND(NND2))
        IF(iCase == iPIPEC) THEN
            DUDYL = DUDYU
        ELSE
            DUDYL = (U1ztL_F0_io(I, 1,   1) - 0.0_WP) / (YCC(1)   - YND(1)   )
        END IF

        Cf_LW_io(I) = 2.0_WP * DUDYL / REN
        Cf_UW_io(I) = 2.0_WP * DUDYU / REN



        U_tau_LW_io(I) = DSQRT(DABS(Cf_LW_io(I)) * 0.5_WP)
        U_tau_UW_io(I) = DSQRT(DABS(Cf_UW_io(I)) * 0.5_WP)

        Re_tau_LW_io(I) = REN * U_tau_LW_io(I)
        Re_tau_UW_io(I) = REN * U_tau_UW_io(I)

        Cf_ave_io(I) = 0.5_WP * (DABS(Cf_LW_io(I)) + DABS(Cf_UW_io(I)))
        U_tau_ave_io(I) = DSQRT(DABS(Cf_ave_io(I)) * 0.5_WP)
        Re_tau_ave_io(I) = REN * U_tau_ave_io(I)

        IF(iThermoDynamics == 1) THEN
            Cfbk_L(I) = 2.0_WP * DUDYL / REN * Mwalz(I, 1) / Dbuk(I) / Ubuk(I) / Dbuk(I)
            Cfbk_U(I) = 2.0_WP * DUDYU / REN * Mwalz(I, 2) / Dbuk(I) / Ubuk(I) / Dbuk(I)
            Cfbk_A(I) = 0.5_WP * (Cfbk_L(I) + Cfbk_U(I))
        END IF

    END DO

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    OPEN(TECFLG, FILE = TRIM(FilePath2) // 'Result.IO.X.DIStributed.Cf.' // TRIM(PNTIM) // '.Profile.plt')
    WRITE(TECFLG, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    IF(iThermoDynamics == 0) THEN
        WRITE(TECFLG, '(A, A, A)') &
        'variables = "X", "Cf_LW", "Cf_UW", "Cf_AVE", ', &
        '"U_tau_LW", "U_tau_UW", "U_tau_AVE", ', &
        '"Re_tau_LW", "Re_tau_UW", "Re_tau_AVE"'
    ELSE IF(iThermoDynamics == 1) THEN
        WRITE(TECFLG, '(A, A, A, A)') &
        'variables = "X", "Cf_LW", "Cf_UW", "Cf_AVE", ', &
        '"U_tau_LW", "U_tau_UW", "U_tau_AVE", ', &
        '"Re_tau_LW", "Re_tau_UW", "Re_tau_AVE", ', &
        '"Cfbk_LW", "Cfbk_UW", "Cfbk_Ave"'
    ELSE
    END IF

    WRITE(TECFLG, '(A)') 'ZONE T = "Time and Z averaged" '

    IF(iThermoDynamics == 0) THEN

        DO I = 1, NCL1_io
            WRITE(TECFLG, '(10ES15.7)') XCC_io(I), Cf_LW_io(I), Cf_UW_io(I), Cf_ave_io(I),    &
            U_tau_LW_io(I), U_tau_UW_io(I), U_tau_ave_io(I), &
            Re_tau_LW_io(I), Re_tau_UW_io(I), Re_tau_ave_io(I)

        END DO

    ELSE IF(iThermoDynamics == 1) THEN

        DO I = 1, NCL1_io
            WRITE(TECFLG, '(13ES15.7)') XCC_io(I), Cf_LW_io(I), Cf_UW_io(I), Cf_ave_io(I),    &
            U_tau_LW_io(I), U_tau_UW_io(I), U_tau_ave_io(I), &
            Re_tau_LW_io(I), Re_tau_UW_io(I), Re_tau_ave_io(I), &
            Cfbk_L(I), Cfbk_U(I), Cfbk_A(I)
        END DO

    ELSE
    END IF


    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!==========================PROFILES AT 10 POSTION =========================================================================
SUBROUTINE POSITIONS_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CHARACTER(15) :: NXSTR
    REAL(WP) :: urms, vrms, WRms, ux, uy, uz
    REAL(WP) :: p, &
    tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
    den, tem, enh, drms, trms, hrms, Thfx, thfy, thfz
    REAL(WP) :: COE1, COE2,COE
    INTEGER(4) :: I, J, K
    INTEGER(4) :: TECFLG1 = 101
    INTEGER(4) :: TECFLG2 = 102


    !==========================pluS ==================================================================
    DO K = 1, NX

        IF(K == 1) THEN
            I = NCL1_io / (NX - 1) / 2
        ELSE IF(K == NX) THEN
            I = NCL1_io- 1
        ELSE
            I = NCL1_io / (NX - 1) * (K - 1)
        END IF

        COE1 = DSQRT(DABS(Cf_LW_io(I) * 0.5_WP))
        COE2 = DSQRT(DABS(Cf_UW_io(I) * 0.5_WP))

        WRITE(NXSTR, '(1ES15.9)') XND_io(I)

        OPEN(TECFLG1, FILE = TRIM(FilePath2) // 'Result.IO.Reynd.Averaged.Flow.' // 'Plus.Profile.X' &
        // TRIM(ADJUSTL(NXSTR)) // '.plt')
        WRITE(TECFLG1, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        WRITE(TECFLG1, '(A)') 'variables = "Y+", "UX +", "Uy+", "UZ +", "UrmS +", "VrmS +", "WRmS +","uv+","uw+","vw+" '
        WRITE(TECFLG1, '(A, 1ES13.5, A)') 'ZONE T = " ',XND_io(I), ' "'


        OPEN(TECFLG2, FILE = TRIM(FilePath2) // 'Result.IO.Favre.Averaged.Flow.' // 'Plus.Profile.X' &
        // TRIM(ADJUSTL(NXSTR)) // '.plt')
        WRITE(TECFLG2, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        WRITE(TECFLG2, '(A)') 'variables = "Y+", "UX +", "Uy+", "UZ +", "UrmS +", "VrmS +", "WRmS +","uv+","uw+","vw+" '
        WRITE(TECFLG2, '(A, 1ES13.5, A)') 'ZONE T = " ',XND_io(I), ' "'

        DO J = 1, NCL2

            IF(YND(J) < 0.0_WP) THEN
                COE = COE1
            ELSE
                COE = COE2
            END IF

            ux  = U1ztL_F0_io(I, J, 1)
            uy  = U1ztL_F0_io(I, J, 2)
            uz  = U1ztL_F0_io(I, J, 3)

            Ruu = ( U2ztL_F0_io(I, J, 1) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 1) ) * D1ztL_F0_io(I, J)
            Ruv = ( U2ztL_F0_io(I, J, 2) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 2) ) * D1ztL_F0_io(I, J)
            Ruw = ( U2ztL_F0_io(I, J, 3) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)
            Rvv = ( U2ztL_F0_io(I, J, 4) - U1ztL_F0_io(I, J, 2) * U1ztL_F0_io(I, J, 2) ) * D1ztL_F0_io(I, J)
            Rvw = ( U2ztL_F0_io(I, J, 5) - U1ztL_F0_io(I, J, 2) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)
            Rww = ( U2ztL_F0_io(I, J,6) - U1ztL_F0_io(I, J, 3) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)

            urms = DSQRT( DABS( Ruu ) )
            vrms = DSQRT( DABS( Rvv ) )
            WRms = DSQRT( DABS( Rww ) )

            WRITE(TECFLG1, '(10ES15.7)')  (1.0_WP - DABS(YCC(J))) * REN * COE, uX / COE, uY / COE, uz/ COE, &
            urms/ COE, vrms/ COE, WRms/ COE, &
            - 1.0_WP * DABS(Ruv) / COE/ COE, Ruw/ COE/ COE,Rvw/ COE/ COE
        END DO
        CLOSE(TECFLG1)


        IF(iThermoDynamics == 1) THEN
            OPEN(TECFLG2, FILE = TRIM(FilePath2) // 'Result.IO.Favre.Averaged.Flow.' &
            // 'Plus.Profile.X' // TRIM(ADJUSTL(NXSTR)) // '.plt')
            WRITE(TECFLG2, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
            WRITE(TECFLG2, '(A)') 'variables = "Y+", "UX +", "Uy+", "UZ +", "UrmS +", "VrmS +", "WRmS +","uv+","uw+","vw+" '
            WRITE(TECFLG2, '(A, 1ES13.5, A)') 'ZONE T = " ',XND_io(I), ' "'


            DO J = 1, NCL2

                IF(YND(J) < 0.0_WP) THEN
                    COE = COE1
                ELSE
                    COE = COE2
                END IF

                ux  = G1ztL_F0_io(I, J, 1) / D1ztL_F0_io(I, J)
                uy  = G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                uz  = G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)

                Ruu = UGztL_F0_io(I, J, 1) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 1) / D1ztL_F0_io(I, J)
                Ruv = UGztL_F0_io(I, J, 2) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                Ruw = UGztL_F0_io(I, J, 3) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)
                Rvv = UGztL_F0_io(I, J, 4) - G1ztL_F0_io(I, J, 2) * G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                Rvw = UGztL_F0_io(I, J, 5) - G1ztL_F0_io(I, J, 2) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)
                Rww = UGztL_F0_io(I, J,6) - G1ztL_F0_io(I, J, 3) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)

                urms = DSQRT( DABS( Ruu ) )
                vrms = DSQRT( DABS( Rvv ) )
                WRms = DSQRT( DABS( Rww ) )

                WRITE(TECFLG2, '(10ES15.7)')  (1.0_WP - DABS(Ycc(J))) * REN * COE, uX / COE, uY / COE, uz/ COE, &
                urms/ COE, vrms/ COE, WRms/ COE, &
                - 1.0_WP * DABS(Ruv) / COE/ COE, Ruw/ COE/ COE,Rvw/ COE/ COE
            END DO
            CLOSE(TECFLG2)
        END IF
    END DO

    !============================= Reynolds averaged========================================
    DO K = 1, NX

        IF(K == 1) THEN
            I = NCL1_io / (NX - 1) / 2
        ELSE IF(K == NX) THEN
            I = NCL1_io- 1
        ELSE
            I = NCL1_io / (NX - 1) * (K - 1)
        END IF


        WRITE(NXSTR, '(1ES15.9)') XND_io(I)
        OPEN(TECFLG1, FILE = TRIM(FilePath2) // 'Result.IO.Reynd.Averaged.Flow.' &
        // 'StAR.Profile.X' // TRIM(ADJUSTL(NXSTR)) // '.plt')
        WRITE(TECFLG1, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        IF(iThermoDynamics == 0) THEN
            WRITE(TECFLG1, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
            ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw"'
        ELSE IF(iThermoDynamics == 1) THEN
            WRITE(TECFLG1, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
            ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
            ' "D", "T", "H", "Drms", "Trms", "Hrms", ', &
            ' "thfx", "thfy", "thfz" '
        END IF
        WRITE(TECFLG1, '(A, 1ES13.5, A)') 'ZONE T = " ',XCC_io(I), ' "'


        DO J = 1, NCL2

            ux  = U1ztL_F0_io(I, J, 1)
            uy  = U1ztL_F0_io(I, J, 2)
            uz  = U1ztL_F0_io(I, J, 3)

            p = U1ztL_F0_io(I, J, 4)

            Ruu = ( U2ztL_F0_io(I, J, 1) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 1) ) * D1ztL_F0_io(I, J)
            Ruv = ( U2ztL_F0_io(I, J, 2) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 2) ) * D1ztL_F0_io(I, J)
            Ruw = ( U2ztL_F0_io(I, J, 3) - U1ztL_F0_io(I, J, 1) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)
            Rvv = ( U2ztL_F0_io(I, J, 4) - U1ztL_F0_io(I, J, 2) * U1ztL_F0_io(I, J, 2) ) * D1ztL_F0_io(I, J)
            Rvw = ( U2ztL_F0_io(I, J, 5) - U1ztL_F0_io(I, J, 2) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)
            Rww = ( U2ztL_F0_io(I, J,6) - U1ztL_F0_io(I, J, 3) * U1ztL_F0_io(I, J, 3) ) * D1ztL_F0_io(I, J)

            tke = 0.5_WP * (Ruu+ Rvv+ Rww)

            IF(iThermoDynamics == 1) THEN
                den = D1ztL_F0_io(I, J)
                tem = T1ztL_F0_io(I, J)
                enh = H1ztL_F0_io(I, J)

                drms = DSQRT(DABS( D2ztL_F0_io(I, J) - D1ztL_F0_io(I, J) * D1ztL_F0_io(I, J) ) )
                trms = DSQRT(DABS( T2ztL_F0_io(I, J) - T1ztL_F0_io(I, J) * T1ztL_F0_io(I, J) ) )
                hrms = DSQRT(DABS( H2ztL_F0_io(I, J) - H1ztL_F0_io(I, J) * H1ztL_F0_io(I, J) ) )

                thfx = ( UHztL_F0_io(I, J, 1) - U1ztL_F0_io(I, J, 1) * H1ztL_F0_io(I, J) ) * D1ztL_F0_io(I, J)
                thfy = ( UHztL_F0_io(I, J, 2) - U1ztL_F0_io(I, J, 2) * H1ztL_F0_io(I, J) ) * D1ztL_F0_io(I, J)
                thfz = ( UHztL_F0_io(I, J, 3) - U1ztL_F0_io(I, J, 3) * H1ztL_F0_io(I, J) ) * D1ztL_F0_io(I, J)
            END IF

            IF(iThermoDynamics == 0) THEN
                WRITE(TECFLG1, '(13ES17.9)') XCC_io(I), YCC(J), ux, uy, uz, p, &
                tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw
            ELSE IF(iThermoDynamics == 1) THEN
                WRITE(TECFLG1, '(22ES17.9)') XCC_io(I), YCC(J), ux, uy, uz, p, &
                tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
                den, tem, enh, drms, trms, hrms, &
                thfx, thfy, thfz
            ELSE
            END IF

        END DO

        CLOSE(TECFLG1)
    END DO


    !=================fARve averaged==============================================================================
    IF(iThermoDynamics == 1) THEN
        DO K = 1, NX

            IF(K == 1) THEN
                I = NCL1_io / (NX - 1) / 2
            ELSE IF(K == NX) THEN
                I = NCL1_io- 1
            ELSE
                I = NCL1_io / (NX - 1) * (K - 1)
            END IF

            WRITE(NXSTR, '(1ES15.9)') XCC_io(I)
            OPEN(TECFLG2, FILE = TRIM(FilePath2) // 'Result.IO.Favre.Averaged.Flow.' &
            // 'StAR.Profile.X' // TRIM(ADJUSTL(NXSTR)) // '.plt')
            WRITE(TECFLG2, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
            WRITE(TECFLG2, '(A, A, A)') 'variables = "X", "Y", "Ux", "Uy", "Uz", "P", ', &
            ' "TKE", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
            ' "D", "T", "H", "Drms", "Trms", "Hrms", ', &
            ' "thfx", "thfy", "thfz" '
            WRITE(TECFLG2, '(A, 1ES13.5, A)') 'ZONE T = " ',XCC_io(I), ' "'



            DO J = 1, NCL2

                ux  = G1ztL_F0_io(I, J, 1) / D1ztL_F0_io(I, J)
                uy  = G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                uz  = G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)

                p = U1ztL_F0_io(I, J, 4) ! same

                Ruu = UGztL_F0_io(I, J, 1) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 1) / D1ztL_F0_io(I, J)
                Ruv = UGztL_F0_io(I, J, 2) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                Ruw = UGztL_F0_io(I, J, 3) - G1ztL_F0_io(I, J, 1) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)
                Rvv = UGztL_F0_io(I, J, 4) - G1ztL_F0_io(I, J, 2) * G1ztL_F0_io(I, J, 2) / D1ztL_F0_io(I, J)
                Rvw = UGztL_F0_io(I, J, 5) - G1ztL_F0_io(I, J, 2) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)
                Rww = UGztL_F0_io(I, J,6) - G1ztL_F0_io(I, J, 3) * G1ztL_F0_io(I, J, 3) / D1ztL_F0_io(I, J)

                tke = 0.5_WP * (Ruu+ Rvv+ Rww)

                den = D1ztL_F0_io(I, J) !same
                tem = T1ztL_F0_io(I, J) !same
                enh = DHztL_F0_io(I, J) / D1ztL_F0_io(I, J)

                drms = DSQRT(DABS( D2ztL_F0_io(I, J) - D1ztL_F0_io(I, J) * D1ztL_F0_io(I, J) ) ) ! same
                trms = DSQRT(DABS( T2ztL_F0_io(I, J) - T1ztL_F0_io(I, J) * T1ztL_F0_io(I, J) ) ) ! same
                hrms = DSQRT(DABS( H2ztL_F0_io(I, J) - DHztL_F0_io(I, J) / D1ztL_F0_io(I, J) * &
                (2.0_WP * H1ztL_F0_io(I, J) - DHztL_F0_io(I, J) / D1ztL_F0_io(I, J) ) ) )

                thfx = GHztL_F0_io(I, J, 1) - G1ztL_F0_io(I, J, 1) * DHztL_F0_io(I, J) / D1ztL_F0_io(I, J)
                thfy = GHztL_F0_io(I, J, 2) - G1ztL_F0_io(I, J, 2) * DHztL_F0_io(I, J) / D1ztL_F0_io(I, J)
                thfz = GHztL_F0_io(I, J, 3) - G1ztL_F0_io(I, J, 3) * DHztL_F0_io(I, J) / D1ztL_F0_io(I, J)

                WRITE(TECFLG2, '(22ES17.9)') XCC_io(I), YCC(J), ux, uy, uz, p, &
                tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
                den, tem, enh, drms, trms, hrms, &
                thfx, thfy, thfz
            END DO
            CLOSE(TECFLG2)
        END DO
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CL2ND_INTP_AVEG_nonXperiodic_io
    USE VARS_AVERAGED_nonXperiodic_io
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: IC, JC, IM, JM, N, M,H,L1, L2, IP, JP, IS, JS
    REAL(WP) :: tlow, tupp
    REAL(WP) :: side1(NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8), side2(NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)

    !============ main DOMAIN ==FOR FLOW============================================
    DO JC = 2, NCL2
        JM = JC - 1
        DO IC = 2, NCL1_io
            IM = IC- 1

            DO N = 1, NDV + 1
                tlow = 0.5_WP * ( U1ztL_F0_io(IM, JM, N) + U1ztL_F0_io(IC, JM, N) )
                tupp = 0.5_WP * ( U1ztL_F0_io(IM, JC, N) + U1ztL_F0_io(IC, JC, N) )
                U1ztL_F0_INTP_io(IC, JC, N) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp
            END DO

            DO N = 1, NDV
                tlow = 0.5_WP * ( G1ztL_F0_io(IM, JM, N) + G1ztL_F0_io(IC, JM, N) )
                tupp = 0.5_WP * ( G1ztL_F0_io(IM, JC, N) + G1ztL_F0_io(IC, JC, N) )
                G1ztL_F0_INTP_io(IC, JC, N) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                tlow = 0.5_WP * ( UPztL_F0_io(IM, JM, N) + UPztL_F0_io(IC, JM, N) )
                tupp = 0.5_WP * ( UPztL_F0_io(IM, JC, N) + UPztL_F0_io(IC, JC, N) )
                UPztL_F0_INTP_io(IC, JC, N) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    L1 = (M * (7-M)) / 2 + N - 3

                    tlow = 0.5_WP * ( U2ztL_F0_io(IM, JM,L1) + U2ztL_F0_io(IC, JM,L1) )
                    tupp = 0.5_WP * ( U2ztL_F0_io(IM, JC,L1) + U2ztL_F0_io(IC, JC,L1) )
                    U2ztL_F0_INTP_io(IC, JC,L1) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                    tlow = 0.5_WP * ( UGztL_F0_io(IM, JM,L1) + UGztL_F0_io(IC, JM,L1) )
                    tupp = 0.5_WP * ( UGztL_F0_io(IM, JC,L1) + UGztL_F0_io(IC, JC,L1) )
                    UGztL_F0_INTP_io(IC, JC,L1) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                    DO H = 1, NDV
                        L2 = M * (6-M) + (N * (7-N)) / 2 + H-8
                        tlow = 0.5_WP * ( UGUztL_F0_io(IM, JM,L2) + UGUztL_F0_io(IC, JM,L2) )
                        tupp = 0.5_WP * ( UGUztL_F0_io(IM, JC,L2) + UGUztL_F0_io(IC, JC,L2) )
                        UGUztL_F0_INTP_io(IC, JC,L2) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp
                    END DO

                END DO
            END DO

        END DO
    END DO


    !============FOR ToP / Bottom. Except four cornerS ==================================
    DO IC = 2, NCL1_io

        DO JC = 1, NND2, NND2 - 1
            IF(JC == 1) THEN
                JP = JC + 1
                JS = JC + 2
            ELSE IF(JC == NND2) THEN
                JP = JC - 1
                JS = JC - 2
            ELSE
                cycle
            END IF

            U1ztL_F0_INTP_io(IC, JC, 4) = 2.0_WP * U1ztL_F0_INTP_io(IC, JP, 4) - 1.0_WP * U1ztL_F0_INTP_io(IC, JS, 4)

            U1ztL_F0_INTP_io(IC, JC, 1:3) = 0.0_WP
            G1ztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP
            UPztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP

            U2ztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP
            UGztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP

            UGUztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP

            IF(iCase == iPIPEC .AND. JC == 1) THEN

                JP = JC + 1
                JS = JC + 2
                U1ztL_F0_INTP_io(IC, JC, 4) = 2.0_WP * U1ztL_F0_INTP_io(IC, JP, 4) - 1.0_WP * U1ztL_F0_INTP_io(IC, JS, 4)

                U1ztL_F0_INTP_io(IC, JC, 1:3) = 2.0_WP * U1ztL_F0_INTP_io(IC, JP, 1:3) - 1.0_WP * U1ztL_F0_INTP_io(IC, JS, 1:3)
                G1ztL_F0_INTP_io(IC, JC,:  ) = 2.0_WP * G1ztL_F0_INTP_io(IC, JC,:  ) - 1.0_WP * G1ztL_F0_INTP_io(IC, JC,:  )
                UPztL_F0_INTP_io(IC, JC,:  ) = 2.0_WP * UPztL_F0_INTP_io(IC, JC,:  ) - 1.0_WP * UPztL_F0_INTP_io(IC, JC,:  )

                U2ztL_F0_INTP_io(IC, JC,:  ) = 2.0_WP * U2ztL_F0_INTP_io(IC, JC,:  ) - 1.0_WP * U2ztL_F0_INTP_io(IC, JC,:  )
                UGztL_F0_INTP_io(IC, JC,:  ) = 2.0_WP * UGztL_F0_INTP_io(IC, JC,:  ) - 1.0_WP * UGztL_F0_INTP_io(IC, JC,:  )

                UGUztL_F0_INTP_io(IC, JC,:  ) = 2.0_WP * UGUztL_F0_INTP_io(IC, JC,:  ) - 1.0_WP * UGUztL_F0_INTP_io(IC, JC,:  )
            END IF
        END DO

    END DO
    !============FOR inlet /Outlet. Except four cornerS ==================================
    DO JC = 2, NCL2
        DO IC = 1, NCL1_io + 1, NCL1_io
            IF(IC == 1) THEN
                IP = IC + 1
                IS = IC + 2
            ELSE IF(IC == NCL1_io + 1) THEN
                IP = IC - 1
                IS = IC - 2
            ELSE
                cycle
            END IF
            U1ztL_F0_INTP_io(IC, JC, :) = 2.0_WP * U1ztL_F0_INTP_io(IP, JC, :)   - 1.0_WP * U1ztL_F0_INTP_io(IS, JC, :)
            G1ztL_F0_INTP_io(IC, JC, :) = 2.0_WP * G1ztL_F0_INTP_io(IP, JC, :)   - 1.0_WP * G1ztL_F0_INTP_io(IS, JC, :)
            UPztL_F0_INTP_io(IC, JC, :) = 2.0_WP * UPztL_F0_INTP_io(IP, JC, :)   - 1.0_WP * UPztL_F0_INTP_io(IS, JC, :)

            U2ztL_F0_INTP_io(IC, JC, :) = 2.0_WP * U2ztL_F0_INTP_io(IP, JC, :)   - 1.0_WP * U2ztL_F0_INTP_io(IS, JC, :)
            UGztL_F0_INTP_io(IC, JC, :) = 2.0_WP * UGztL_F0_INTP_io(IP, JC, :)   - 1.0_WP * UGztL_F0_INTP_io(IS, JC, :)

            UGUztL_F0_INTP_io(IC, JC, :) = 2.0_WP * UGUztL_F0_INTP_io(IP, JC, :)  - 1.0_WP * UGUztL_F0_INTP_io(IS, JC, :)

        END DO
    END DO

    !==============For Four coNNerS ======================================================
    DO IC = 1, NCL1_io + 1, NCL1_io
        IF(IC == 1) THEN
            IP = IC + 1
            IS = IC + 2
        ELSE IF(IC == NCL1_io + 1) THEN
            IP = IC - 1
            IS = IC - 2
        ELSE
            cycle
        END IF
        DO JC = 1, NCL2 + 1, NCL2
            IF(JC == 1) THEN
                JP = JC + 1
                JS = JC + 2
            ELSE IF(JC == NND2) THEN
                JP = JC - 1
                JS = JC - 2
            ELSE
                cycle
            END IF

            U1ztL_F0_INTP_io(IC, JC, 1:3) = 0.0_WP
            G1ztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP
            UPztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP

            U2ztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP
            UGztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP

            UGUztL_F0_INTP_io(IC, JC,:  ) = 0.0_WP


            IF(iCase == iPIPEC .AND. JC == 1) THEN
                JP = JC + 1
                JS = JC + 2
                SIDE1(1:3) = 2.0_WP * U1ztL_F0_INTP_io(IC, JP, 1:3) - 1.0_WP * U1ztL_F0_INTP_io(IC, JS, 1:3)
                SIDE2(1:3) = 2.0_WP * U1ztL_F0_INTP_io(IP, JC, 1:3) - 1.0_WP * U1ztL_F0_INTP_io(IS, JC, 1:3)
                U1ztL_F0_INTP_io(IC, JC, 1:3) = 0.5_WP * (SIDE1(1:3) + SIDE2(1:3))

                SIDE1(1:3) = 2.0_WP * G1ztL_F0_INTP_io(IC, JP, :) - 1.0_WP * G1ztL_F0_INTP_io(IC, JS, :)
                SIDE2(1:3) = 2.0_WP * G1ztL_F0_INTP_io(IP, JC, :) - 1.0_WP * G1ztL_F0_INTP_io(IS, JC, :)
                G1ztL_F0_INTP_io(IC, JC, 1:3) = 0.5_WP * (SIDE1(1:3) + SIDE2(1:3))

                SIDE1(1:3) = 2.0_WP * UPztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UPztL_F0_INTP_io(IC, JS, :)
                SIDE2(1:3) = 2.0_WP * UPztL_F0_INTP_io(IP, JC, :) - 1.0_WP * UPztL_F0_INTP_io(IS, JC, :)
                UPztL_F0_INTP_io(IC, JC, 1:3) = 0.5_WP * (SIDE1(1:3) + SIDE2(1:3))

                SIDE1(1:6) = 2.0_WP * U2ztL_F0_INTP_io(IC, JP, :) - 1.0_WP * U2ztL_F0_INTP_io(IC, JS, :)
                SIDE2(1:6) = 2.0_WP * U2ztL_F0_INTP_io(IP, JC, :) - 1.0_WP * U2ztL_F0_INTP_io(IS, JC, :)
                U2ztL_F0_INTP_io(IC, JC, 1:6) = 0.5_WP * (SIDE1(1:6) + SIDE2(1:6))

                SIDE1(1:6) = 2.0_WP * UGztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UGztL_F0_INTP_io(IC, JS, :)
                SIDE2(1:6) = 2.0_WP * UGztL_F0_INTP_io(IP, JC, :) - 1.0_WP * UGztL_F0_INTP_io(IS, JC, :)
                UGztL_F0_INTP_io(IC, JC, 1:6) = 0.5_WP * (SIDE1(1:6) + SIDE2(1:6))

                SIDE1(1:10) = 2.0_WP * UGUztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UGUztL_F0_INTP_io(IC, JS, :)
                SIDE2(1:10) = 2.0_WP * UGUztL_F0_INTP_io(IP, JC, :) - 1.0_WP * UGUztL_F0_INTP_io(IS, JC, :)
                UGUztL_F0_INTP_io(IC, JC, 1:10) = 0.5_WP * (SIDE1(1:10) + SIDE2(1:10))

            END IF

            SIDE1(1) = 2.0_WP * U1ztL_F0_INTP_io(IC, JP, 4) - 1.0_WP * U1ztL_F0_INTP_io(IC, JS, 4)
            SIDE2(1) = 2.0_WP * U1ztL_F0_INTP_io(IP, JC, 4) - 1.0_WP * U1ztL_F0_INTP_io(IS, JC, 4)
            U1ztL_F0_INTP_io(IC, JC, 4) = 0.5_WP * (SIDE1(1) + SIDE2(1))
        END DO
    END DO

    !======================= main DOMAIN FOR THERMAL FIELD========================================
    IF(iThermoDynamics == 1) THEN

        DO JC = 2, NCL2
            JM = JC - 1
            DO IC = 2, NCL1_io
                IM = IC- 1
                !============== D1 AND D2 ===========================================
                tlow = 0.5_WP * ( D1ztL_F0_io(IM, JM) + D1ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( D1ztL_F0_io(IM, JC) + D1ztL_F0_io(IC, JC) )
                D1ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                tlow = 0.5_WP * ( D2ztL_F0_io(IM, JM) + D2ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( D2ztL_F0_io(IM, JC) + D2ztL_F0_io(IC, JC) )
                D2ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                !============== T1 AND T2 ===========================================
                tlow = 0.5_WP * ( T1ztL_F0_io(IM, JM) + T1ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( T1ztL_F0_io(IM, JC) + T1ztL_F0_io(IC, JC) )
                T1ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                tlow = 0.5_WP * ( T2ztL_F0_io(IM, JM) + T2ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( T2ztL_F0_io(IM, JC) + T2ztL_F0_io(IC, JC) )
                T2ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                !============== H1 AND H2 ===========================================
                tlow = 0.5_WP * ( H1ztL_F0_io(IM, JM) + H1ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( H1ztL_F0_io(IM, JC) + H1ztL_F0_io(IC, JC) )
                H1ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                tlow = 0.5_WP * ( H2ztL_F0_io(IM, JM) + H2ztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( H2ztL_F0_io(IM, JC) + H2ztL_F0_io(IC, JC) )
                H2ztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                !============== DH ==================================================
                tlow = 0.5_WP * ( DHztL_F0_io(IM, JM) + DHztL_F0_io(IC, JM) )
                tupp = 0.5_WP * ( DHztL_F0_io(IM, JC) + DHztL_F0_io(IC, JC) )
                DHztL_F0_INTP_io(IC, JC) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                !==============UH AND GH ===========================================
                DO N = 1, NDV
                    tlow = 0.5_WP * ( UHztL_F0_io(IM, JM, N) + UHztL_F0_io(IC, JM, N) )
                    tupp = 0.5_WP * ( UHztL_F0_io(IM, JC, N) + UHztL_F0_io(IC, JC, N) )
                    UHztL_F0_INTP_io(IC, JC, N) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp

                    tlow = 0.5_WP * ( GHztL_F0_io(IM, JM, N) + GHztL_F0_io(IC, JM, N) )
                    tupp = 0.5_WP * ( GHztL_F0_io(IM, JC, N) + GHztL_F0_io(IC, JC, N) )
                    GHztL_F0_INTP_io(IC, JC, N) = YCL2ND_WFB(JC) * Tlow + YCL2ND_WFF(JC) * Tupp
                END DO

            END DO
        END DO

        !============FOR ToP / Bottom. Except four cornerS ==================================
        DO IC = 2, NCL1_io
            DO JC = 1, NND2, NND2 - 1
                IF(JC == 1) THEN
                    JP = JC + 1
                    JS = JC + 2
                ELSE IF(JC == NND2) THEN
                    JP = JC - 1
                    JS = JC - 2
                ELSE
                    cycle
                END IF

                D1ztL_F0_INTP_io(IC, JC) = 2.0_WP * D1ztL_F0_INTP_io(IC, JP) - 1.0_WP * D1ztL_F0_INTP_io(IC, JS)
                D2ztL_F0_INTP_io(IC, JC) = 2.0_WP * D2ztL_F0_INTP_io(IC, JP) - 1.0_WP * D2ztL_F0_INTP_io(IC, JS)
                H1ztL_F0_INTP_io(IC, JC) = 2.0_WP * H1ztL_F0_INTP_io(IC, JP) - 1.0_WP * H1ztL_F0_INTP_io(IC, JS)
                H2ztL_F0_INTP_io(IC, JC) = 2.0_WP * H2ztL_F0_INTP_io(IC, JP) - 1.0_WP * H2ztL_F0_INTP_io(IC, JS)
                T1ztL_F0_INTP_io(IC, JC) = 2.0_WP * T1ztL_F0_INTP_io(IC, JP) - 1.0_WP * T1ztL_F0_INTP_io(IC, JS)
                T2ztL_F0_INTP_io(IC, JC) = 2.0_WP * T2ztL_F0_INTP_io(IC, JP) - 1.0_WP * T2ztL_F0_INTP_io(IC, JS)
                DHztL_F0_INTP_io(IC, JC) = 2.0_WP * DHztL_F0_INTP_io(IC, JP) - 1.0_WP * DHztL_F0_INTP_io(IC, JS)
                UHztL_F0_INTP_io(IC, JC, :) = 0.0_WP
                GHztL_F0_INTP_io(IC, JC, :) = 0.0_WP
                IF(iCase == iPIPEC .AND. JC == 1) THEN
                    JP = JC + 1
                    JS = JC + 2
                    UHztL_F0_INTP_io(IC, JC, :) = 2.0_WP * UHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UHztL_F0_INTP_io(IC, JS, :)
                    GHztL_F0_INTP_io(IC, JC, :) = 2.0_WP * GHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * GHztL_F0_INTP_io(IC, JS, :)
                END IF

            END DO
        END DO
        !============FOR inlet /Outlet. Except four cornerS ==================================
        DO JC = 2, NCL2
            DO IC = 1, NCL1_io + 1, NCL1_io
                IF(IC == 1) THEN
                    IP = IC + 1
                    IS = IC + 2
                ELSE IF(IC == NCL1_io + 1) THEN
                    IP = IC - 1
                    IS = IC - 2
                ELSE
                    cycle
                END IF
                D1ztL_F0_INTP_io(IC, JC) = 2.0_WP * D1ztL_F0_INTP_io(IP, JC) - 1.0_WP * D1ztL_F0_INTP_io(IS, JC)
                D2ztL_F0_INTP_io(IC, JC) = 2.0_WP * D2ztL_F0_INTP_io(IP, JC) - 1.0_WP * D2ztL_F0_INTP_io(IS, JC)
                H1ztL_F0_INTP_io(IC, JC) = 2.0_WP * H1ztL_F0_INTP_io(IP, JC) - 1.0_WP * H1ztL_F0_INTP_io(IS, JC)
                H2ztL_F0_INTP_io(IC, JC) = 2.0_WP * H2ztL_F0_INTP_io(IP, JC) - 1.0_WP * H2ztL_F0_INTP_io(IS, JC)
                T1ztL_F0_INTP_io(IC, JC) = 2.0_WP * T1ztL_F0_INTP_io(IP, JC) - 1.0_WP * T1ztL_F0_INTP_io(IS, JC)
                T2ztL_F0_INTP_io(IC, JC) = 2.0_WP * T2ztL_F0_INTP_io(IP, JC) - 1.0_WP * T2ztL_F0_INTP_io(IS, JC)
                DHztL_F0_INTP_io(IC, JC) = 2.0_WP * DHztL_F0_INTP_io(IP, JC) - 1.0_WP * DHztL_F0_INTP_io(IS, JC)
                UHztL_F0_INTP_io(IC, JC, :) = 0.0_WP
                GHztL_F0_INTP_io(IC, JC, :) = 0.0_WP
                IF(iCase == iPIPEC .AND. JC == 1) THEN
                    JP = JC + 1
                    JS = JC + 2
                    UHztL_F0_INTP_io(IC, JC, :) = 2.0_WP * UHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UHztL_F0_INTP_io(IC, JS, :)
                    GHztL_F0_INTP_io(IC, JC, :) = 2.0_WP * GHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * GHztL_F0_INTP_io(IC, JS, :)
                END IF

            END DO
        END DO

        !==============For Four coNNerS ======================================================
        DO IC = 1, NCL1_io + 1, NCL1_io
            IF(IC == 1) THEN
                IP = IC + 1
                IS = IC + 2
            ELSE IF(IC == NCL1_io + 1) THEN
                IP = IC - 1
                IS = IC - 2
            ELSE
                cycle
            END IF
            DO JC = 1, NCL2 + 1, NCL2
                IF(JC == 1) THEN
                    JP = JC + 1
                    JS = JC + 2
                ELSE IF(JC == NND2) THEN
                    JP = JC - 1
                    JS = JC - 2
                ELSE
                    cycle
                END IF



                SIDE1(1) = 2.0_WP * D1ztL_F0_INTP_io(IC, JP) - 1.0_WP * D1ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * D1ztL_F0_INTP_io(IP, JC) - 1.0_WP * D1ztL_F0_INTP_io(IS, JC)
                D1ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * D2ztL_F0_INTP_io(IC, JP) - 1.0_WP * D2ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * D2ztL_F0_INTP_io(IP, JC) - 1.0_WP * D2ztL_F0_INTP_io(IS, JC)
                D2ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * H1ztL_F0_INTP_io(IC, JP) - 1.0_WP * H1ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * H1ztL_F0_INTP_io(IP, JC) - 1.0_WP * H1ztL_F0_INTP_io(IS, JC)
                H1ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * H2ztL_F0_INTP_io(IC, JP) - 1.0_WP * H2ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * H2ztL_F0_INTP_io(IP, JC) - 1.0_WP * H2ztL_F0_INTP_io(IS, JC)
                H2ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * T1ztL_F0_INTP_io(IC, JP) - 1.0_WP * T1ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * T1ztL_F0_INTP_io(IP, JC) - 1.0_WP * T1ztL_F0_INTP_io(IS, JC)
                T1ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * T2ztL_F0_INTP_io(IC, JP) - 1.0_WP * T2ztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * T2ztL_F0_INTP_io(IP, JC) - 1.0_WP * T2ztL_F0_INTP_io(IS, JC)
                T2ztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                SIDE1(1) = 2.0_WP * DHztL_F0_INTP_io(IC, JP) - 1.0_WP * DHztL_F0_INTP_io(IC, JS)
                SIDE2(1) = 2.0_WP * DHztL_F0_INTP_io(IP, JC) - 1.0_WP * DHztL_F0_INTP_io(IS, JC)
                DHztL_F0_INTP_io(IC, JC) = 0.5_WP * (SIDE1(1) + SIDE2(1))

                UHztL_F0_INTP_io(IC, JC, :) = 0.0_WP
                GHztL_F0_INTP_io(IC, JC, :) = 0.0_WP

                IF(iCase == iPIPEC .AND. JC == 1) THEN
                    JP = JC + 1
                    JS = JC + 2
                    SIDE1(1:3) = 2.0_WP * UHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * UHztL_F0_INTP_io(IC, JS, :)
                    SIDE2(1:3) = 2.0_WP * UHztL_F0_INTP_io(IP, JC, :) - 1.0_WP * UHztL_F0_INTP_io(IS, JC, :)
                    UHztL_F0_INTP_io(IC, JC, :) = 0.5_WP * (SIDE1(1:3) + SIDE2(1:3))

                    SIDE1(1:3) = 2.0_WP * GHztL_F0_INTP_io(IC, JP, :) - 1.0_WP * GHztL_F0_INTP_io(IC, JS, :)
                    SIDE2(1:3) = 2.0_WP * GHztL_F0_INTP_io(IP, JC, :) - 1.0_WP * GHztL_F0_INTP_io(IS, JC, :)
                    GHztL_F0_INTP_io(IC, JC, :) = 0.5_WP * (SIDE1(1:3) + SIDE2(1:3))
                END IF

            END DO
        END DO


    END IF


    RETURN
END SUBROUTINE
