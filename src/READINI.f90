!**********************************************************************************************************************************
!> @brief
!>      Read in the user-specified parameters from readdata.ini
!> @details
!> SUBROUTINE: READINI
!>             Read in the user-specified parameters from readdata.ini
!> SUBROUTINE: BCAST_READINI
!>             Broadcast the read information
!> @note
!> The reading subsequence in this subroutine should be consistent with the readdata.ini.
!> @todo
! REVISION HISTORY:
! 06/12 /2013- Initial Version, by Wei Wang (wei.wang@sheffield.ac.uk)
! 29/07/ 2020- revised by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE READINI
    USE init_info
    USE mpi_info
    USE mesh_info
    USE thermal_info
    USE CHT_info
    IMPLICIT NONE

    CHARACTER(len = 80) :: SECT
    CHARACTER(len = 128) :: SKEY
    INTEGER(4) :: IOS = 0
    INTEGER(4) :: INI = 13
    INTEGER(4) :: LENS
    INTEGER(4) :: iTMP
    INTEGER(4) :: N
    REAL(WP) :: RTMP
    CHARACTER(256) :: STMP
    LOGICAL :: File_exists = .FALSE.


    PI = 2.0_WP * (DASIN(1.0_WP))

    IF(MYID /= 0) RETURN

    OPEN(INI, FILE = 'readdata.ini', STATUS = 'old', IOSTAT = IOS)
    IF(IOS /= 0)  &
    CALL ERRHDL(' File readdata.ini cannot be found.', MYID)

    !========== Readin Section [flowtype]======================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[flowtype]')  &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, iCase
    READ(INI, *) SKEY, iDomain
    READ(INI, *) SKEY, iThermoDynamics
    READ(INI, *) SKEY, iCHT

    IF(iCase == ICHANL) THEN
        CALL CHKHDL   ('iCase = Plane Channel (Cartesian Cooridnates) ', MYID)
        zoneNameView = 'CHANNEL'
    ELSE IF(iCase == IBox3P) THEN
        CALL CHKHDL   ('iCase = Box with 3 periodic b.c (Cartesian Cooridnates) ', MYID)
        zoneNameView = 'BOX'
    ELSE IF(iCase == IPIPEC) THEN
        CALL CHKHDL   ('iCase = Circular Tube (cylindrical Cooridnates)', MYID)
        zoneNameView = 'PIPE'
    ELSE IF(iCase == IANNUL) THEN
        CALL CHKHDL   ('iCase = Annular  Tube (cylindrical Cooridnates)', MYID)
        zoneNameView = 'ANNULAR'
    ELSE
        iCase = ICHANL
        zoneNameView = 'CHANNEL'
        CALL CHKHDL   ('iCase = Plane Channel (Cartesian Cooridnates) ', MYID)
    END IF

    IF(iDomain == ITG) THEN
        CALL CHKHDL   ('iDomain = Turbulent flow generator only with periodic streamwise direction.', MYID)
        TgFlowFlg = .TRUE.
        IoFlowFlg = .FALSE.
    ELSE IF(iDomain == IIO) THEN
        CALL CHKHDL   ('iDomain = The main flow domain only with periodic streamwise direction. ', MYID)
        TgFlowFlg = .FALSE.
        IoFlowFlg = .TRUE.
    ELSE IF(iDomain == ITGIO) THEN
        CALL CHKHDL   ('iDomain = Turbulent flow generator and the main flow domain with inlet and outlet.', MYID)
        TgFlowFlg = .TRUE.
        IoFlowFlg = .TRUE.
    ELSE
        iDomain = IIO
        CALL CHKHDL   ('iDomain = The main flow domain only with periodic streamwise direction.', MYID)
        TgFlowFlg = .FALSE.
        IoFlowFlg = .TRUE.
    END IF

    IF(TgFlowFlg .AND. iCase == IBox3P) &
    CALL ERRHDL('The Box Domain with 3 periodic BCs can only be in the main flow.', MYID)

    IF(TgFlowFlg) OPEN(logflg_tg, FILE = TRIM(FilePath0) // 'monitor.tg.' //fllog)
    IF(IoFlowFlg) OPEN(logflg_io, FILE = TRIM(FilePath0) // 'monitor.io.' //fllog)

    IF(iThermoDynamics == 0) THEN
        CALL CHKHDL    ('iThermoDynamics = only flow field, no thermodynamics.', MYID)
        IF(iCHT /= 0) iCHT = 0
    ELSE
        iThermoDynamics = 1
        CALL CHKHDL    ('iThermoDynamics = flow and thermal fields.', MYID)
    END IF

    IF(iCHT == 0) THEN
        CALL CHKHDL   ('iCHT = conjugate Heat Transfer is not considered.', MYID)
    ELSE
        iCHT = 1
        CALL CHKHDL   ('iCHT = conjugate Heat Transfer is considered.', MYID)
    END IF

    !========== Readin Section [geometry]======================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[geometry]')  &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, HX_tg
    READ(INI, *) SKEY, Hx_io
    READ(INI, *) SKEY, HZ
    READ(INI, *) SKEY, HYB
    READ(INI, *) SKEY, HYT

    IF(iCase == ICHANL) THEN
        HYB = -1.0_WP
        HYT = 1.0_WP
    END IF

    IF(iCase == IPIPEC) THEN
        HZ = 2.0_WP * PI
        HYB = 0.0_WP
        HYT = 1.0_WP
    END IF

    IF(iCase == IANNUL) THEN
        HZ = 2.0_WP * PI
        HYT = 1.0_WP
    END IF

    IF(iCase == IBox3P) THEN
        HX_io = 2.0_WP * PI
        HZ = 2.0_WP * PI
        HYB = -PI
        HYT = PI
    END IF

    !IF(NFLOW== 1) CALL CHKHDL ('NFLOW=                X Streamwise Flow', MYID)
    !IF(NFLOW== 2) CALL CHKHDL ('NFLOW=                Y Streamwise Flow', MYID)
    !IF(NFLOW== 3) CALL CHKHDL ('NFLOW=                Z Streamwise Flow', MYID)

    IF(TgFlowFlg) CALL CHKRLHDL  ('HX_tg = ', MYID, HX_tg)
    IF(IoFlowFlg) CALL CHKRLHDL  ('HX_io = ', MYID, HX_io)
    CALL CHKRLHDL  ('HZ = ', MYID, HZ)
    CALL CHKRLHDL  ('Y/R (bottom) = ', MYID, HYB)
    CALL CHKRLHDL  ('Y/R (top) = ', MYID, HYT)

    !========== Readin Section [mesh]======================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[mesh]')       &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, NCL1_tg
    READ(INI, *) SKEY, NCL1_io
    READ(INI, *) SKEY, NCL3
    READ(INI, *) SKEY, NCL2
    READ(INI, *) SKEY, STR2

    IF(TgFlowFlg) CALL CHKINTHDL ('NCL1(for TG) = ', MYID, NCL1_tg )
    IF(IoFlowFlg) CALL CHKINTHDL ('NCL1(for IO) = ', MYID, NCL1_io )
    CALL CHKINTHDL ('NCL3 = ', MYID, NCL3)
    CALL CHKINTHDL ('NCL2 = ', MYID, NCL2)
    CALL CHKRLHDL  ('STR2 (mesh stretching) = ', MYID, STR2)

    IF(iDomain == ITG) THEN
        IF(NCL1_tg < 2 .OR. HX_tg < 0.0) CALL ERRHDL(' Error in NCL1_tg or HX_tg ', MYID)
    END IF

    IF(iDomain == IIO) THEN
        IF(NCL1_io < 2 .OR. HX_io < 0.0) CALL ERRHDL(' Error in NCL1_io or HX_io ', MYID)
    END IF

    IF(iDomain == ITGIO) THEN
        IF(NCL1_io < 2 .OR. HX_io < 0.0) CALL ERRHDL(' Error in NCL1_io or HX_io ', MYID)
        IF(NCL1_tg < 2 .OR. HX_tg < 0.0) CALL ERRHDL(' Error in NCL1_tg or HX_tg ', MYID)
    END IF

    !========== Readin Section [boundary]======================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[boundary]')       &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, BCX_tg(1), BCX_tg(2)
    READ(INI, *) SKEY, BCX_io(1), BCX_io(2)
    READ(INI, *) SKEY, BCZ(1), BCZ(2)
    READ(INI, *) SKEY, BCY(1), BCY(2)

    IF (TgFlowFlg) THEN
        IF(BCX_tg(1) == 3 ) BCX_tg(2) = 3
        IF(BCX_tg(2) == 3 ) BCX_tg(1) = 3
        CALL CHECK_BC(BCX_tg(1), ' TG Inlet')
        CALL CHECK_BC(BCX_tg(2), ' TG Outlet')
    END IF

    IF (IoFlowFlg ) THEN
        IF(BCX_io(1) == 3 ) BCX_io(2) = 3
        IF(BCX_io(2) == 3 ) BCX_io(1) = 3
        CALL CHECK_BC(BCX_io(1), ' Main Domain Inlet')
        CALL CHECK_BC(BCX_io(2), ' Main Domain Outlet')
    END IF

    IF(BCZ(1) == 3 ) BCZ(2) = 3
    IF(BCZ(2) == 3 ) BCZ(1) = 3

    IF(BCZ(1) /= 3) THEN
        CALL ERRHDL('The spanwise direction MUST BE periodic', MYID)
    ELSE
        CALL CHKHDL(' Spanwise BC = Periodic', MYID)
    END IF

    CALL CHECK_BC(BCY(1), ' Y/Wall-normal')

    !========== Readin Section [fluid]=====================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[fluid]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, REN
    READ(INI, *) SKEY, ReIni
    READ(INI, *) SKEY, TLgRe
    READ(INI, *) SKEY, iFlowDriven
    READ(INI, *) SKEY, Cf_Given

    CALL CHKRLHDL('Reynolds number = ', MYID, REN)
    CALL CHKRLHDL('Initial Re No.  = ', MYID, ReIni)
    CALL CHKRLHDL('Time for Ini R  = ', MYID, TLgRe)
    IF(iFlowDriven == 1) CALL CHKHDL ('iFlowDriven = Constant mass flux driven', MYID)
    IF(iFlowDriven == 2) THEN
        CALL CHKRLHDL('Given Cf = ', MYID, Cf_Given)
        CALL CHKHDL ('iFlowDriven = Constant pressure gradient driven', MYID)
    END IF

    !========== Readin Section [thermohydraulics]=====================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[thermohydraulics]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, iFluidMedia
    READ(INI, *) SKEY, iGravity
    READ(INI, *) SKEY, L0
    READ(INI, *) SKEY, T0
    READ(INI, *) SKEY, Ti
    READ(INI, *) SKEY, iThermalWallType(iBotWall), iThermalWallType(iTopWall)
    READ(INI, *) SKEY, thermalWallBC_Dim(iBotWall), thermalWallBC_Dim(iTopWall)

    IF (iThermoDynamics == 1) THEN
        IF(iFluidMedia == iScpWater)      CALL CHKHDL ('iFluidMedia = Supercritical Water', MYID)
        IF(iFluidMedia == iScpCO2)        CALL CHKHDL ('iFluidMedia = Supercritical CO2', MYID)
        IF(iFluidMedia == iLiquidSodium)  CALL CHKHDL ('iFluidMedia = Liquid Sodium', MYID)
        IF(iFluidMedia == iLiquidLead)    CALL CHKHDL ('iFluidMedia = Liquid Lead', MYID)
        IF(iFluidMedia == iLiquidBismuth) CALL CHKHDL ('iFluidMedia = Liquid Bismuth', MYID)
        IF(iFluidMedia == iLiquidLBE)     CALL CHKHDL ('iFluidMedia = Liquid LBE', MYID)

        SELECT CASE (iFluidMedia)
        CASE (iScpWater)
            iThermoProperty = search_table
            NISTFLNM = 'NIST_WATER_23.5MP.DAT'
            INQUIRE(FILE = TRIM(NISTFLNM), exist = File_exists)
            IF(File_exists) THEN
                CALL CHKHDL('PROPREF FILE = ' // TRIM(NISTFLNM), MYID)
            ELSE
                CALL ERRHDL('File does not exist. Please copy ' // TRIM(NISTFLNM) // ' to the current directory.', MYID)
            END IF
        CASE (iScpCO2)
            iThermoProperty = search_table
            NISTFLNM = 'NIST_CO2_8MP.DAT'
            INQUIRE(FILE = TRIM(NISTFLNM), exist = File_exists)
            IF(File_exists) THEN
                CALL CHKHDL('PROPREF FILE = ' // TRIM(NISTFLNM), MYID)
            ELSE
                CALL ERRHDL('File does not exist. Please copy ' // TRIM(NISTFLNM) // ' to the current directory.', MYID)
            END IF
        CASE DEFAULT
            iThermoProperty = properties_functions
        END SELECT

        IF(iThermoProperty == search_table) CALL CHKHDL('iThermoProperty = from searching table', MYID)
        IF(iThermoProperty == properties_functions) &
        CALL CHKHDL('iThermoProperty = Thermodynamic properties by functions of Temperature', MYID)

        IF(iGravity == 0)  THEN
            CALL CHKHDL ('Gravity = Not considered', MYID)
        ELSE IF(iGravity == 1) THEN
            CALL CHKHDL ('Gravity = +X direction = vertical downwards flow.', MYID)
        ELSE IF(iGravity == -1) THEN
            CALL CHKHDL ('Gravity = -X direction = vertical upwards flow.', MYID)
        ELSE IF(iGravity == -2) THEN
            CALL CHKHDL ('Gravity = -Y direction = horizontal flow', MYID)
        ELSE
            CALL ERRHDL ('The specified gravity direction is not valid. Gravity is NOT considered', MYID)
            iGravity = 0
        END IF

        CALL CHKRLHDL('Ref: Length (m) = ', MYID, L0)
        CALL CHKRLHDL('Ref: Temperature (K) T0 = ', MYID, T0)
        CALL CHKRLHDL('Inlet: Temperature (K) Ti =  ', MYID, Ti)


        IF(iCase == IPIPEC) THEN
            iThermalWallType(iBotWall) = 0
            thermalWallBC_Dim(iBotWall) = 0.0_WP
        END IF

        IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
            CALL CHKHDL    ('iThermalWallType (bottom wall) = Constant Wall Heat Flux', MYID)
            CALL CHKRLHDL  ('Wall Heat flux (W/M2) on the bottom walL = ', MYID, thermalWallBC_Dim(iBotWall))
        ELSE IF (iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            CALL CHKHDL    ('iThermalWallType (bottomwall) = Constant Wall Temperature', MYID)
            CALL CHKRLHDL  ('Wall Temperature (K) on the bottom walL = ', MYID, thermalWallBC_Dim(iBotWall))
        ELSE
        END IF

        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
            CALL CHKHDL    ('iThermalWallType (topwall) = Constant Wall Heat Flux', MYID)
            CALL CHKRLHDL  ('Wall Heat flux (W/M2) on the top walL = ', MYID, thermalWallBC_Dim(iTopWall))
        ELSE IF (iThermalWallType(iTopWall) == BC_Fixed_Temperature) THEN
            CALL CHKHDL    ('iThermalWallType (topmwall) = Constant Wall Temperature', MYID)
            CALL CHKRLHDL  ('Wall Temperature (K) on the top walL = ', MYID, thermalWallBC_Dim(iTopWall))
        ELSE
        END IF
    ELSE
        CALL CHKHDL( SECT(1:LENS)//' is not considered.', MYID)
    END IF

    !========== Readin Section [conjugateHeatTransfer]=====================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[conjugateHeatTransfer]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, HY_Solid(1), HY_Solid(2)
    READ(INI, *) SKEY, HXs_Solid(1), HXs_Solid(2)
    READ(INI, *) SKEY, HXe_Solid(1), HXe_Solid(2)
    READ(INI, *) SKEY, NCL2_Solid(1), NCL2_Solid(2)
    READ(INI, *) SKEY, Cp_Solid(1), Cp_Solid(2)
    READ(INI, *) SKEY, D_Solid(1), D_Solid(2)
    READ(INI, *) SKEY, K_Solid(1), K_Solid(2)
    IF(iCHT == 1) THEN
        CALL CHKINTHDL ('NCL2_Solid_Top = ', MYID, NCL2_Solid(2))
        CALL CHKINTHDL ('NCL2_Solid_Bot = ', MYID, NCL2_Solid(1))
        CALL CHKRLHDL  ('HY_Solid_Top   = ', MYID, HY_Solid(2))
        CALL CHKRLHDL  ('HY_Solid_Bot   = ', MYID, HY_Solid(1))
        CALL CHKRLHDL  ('HXs_Solid_Top  = ', MYID, HXs_Solid(2))
        CALL CHKRLHDL  ('HXs_Solid_Bot  = ', MYID, HXs_Solid(1))
        CALL CHKRLHDL  ('HXe_Solid_Top  = ', MYID, HXe_Solid(2))
        CALL CHKRLHDL  ('HXe_Solid_Bot  = ', MYID, HXe_Solid(1))
        CALL CHKRLHDL  ('Cp_Solid_Top   = ', MYID, Cp_Solid(2))
        CALL CHKRLHDL  ('Cp_Solid_Bot   = ', MYID, Cp_Solid(1))
        CALL CHKRLHDL  ('D_Solid_Top    = ', MYID, D_Solid(2))
        CALL CHKRLHDL  ('D_Solid_Bot    = ', MYID, D_Solid(1))
        CALL CHKRLHDL  ('K_Solid_Top    = ', MYID, K_Solid(2))
        CALL CHKRLHDL  ('K_Solid_Bot    = ', MYID, K_Solid(1))
    ELSE
        CALL CHKHDL( SECT(1:LENS)//' is not considered.', MYID)
    END IF

    !========== Readin Section [initialisation]=====================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[initialisation]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, iIniField_tg
    READ(INI, *) SKEY, iIniField_io
    READ(INI, *) SKEY, iRandomType
    READ(INI, *) SKEY, VPERG
    READ(INI, *) SKEY, SVPERG
    READ(INI, *) SKEY, TimeReStart_tg
    READ(INI, *) SKEY, TimeReStart_io
    READ(INI, *) SKEY, iIniFieldType
    READ(INI, *) SKEY, iIniFieldTime

    IF(iCase == IBox3P) iRandomType = flag_random_sine

    IF(TgFlowFlg) THEN
        IF(iIniField_tg == IniField_random)        CALL CHKHDL ('iIniField_tg = TG: Start from random flow field', MYID)
        IF(iIniField_tg == IniField_extrapolation) &
        CALL CHKHDL ('iIniField_tg = TG: Start from interpolation of a coarse mesh', MYID)
        IF(iIniField_tg == IniField_reStart)       CALL CHKHDL ('iIniField_tg = TG: ReStart from last step', MYID)
        IF(iIniField_tg == IniField_random) THEN
            CALL CHKRLHDL  ('VPERG (initial velocity perturbation)  = ', MYID, VPERG)
            CALL CHKRLHDL  ('SVPERG (initial velocity perturbuatio near 1/4 Y-wall) = ', MYID, SVPERG)
        END IF
    END IF

    IF(IoFlowFlg) THEN
        IF(iIniField_io == IniField_random)        CALL CHKHDL ('iIniField_io = IO: Start from random flow field', MYID)
        IF(iIniField_io == IniField_extrapolation) &
        CALL CHKHDL ('iIniField_io = IO: Start from interpolation of a coarse mesh', MYID)
        IF(iIniField_io == IniField_reStart)       CALL CHKHDL ('iIniField_io = IO: ReStart from last step', MYID)
        IF(iIniField_io == IniField_extrapolation .OR. iIniField_io == IniField_reStart) THEN
            IF(iThermoDynamics == 0) iIniFieldType = 1
            IF(iIniFieldType == 0)  CALL CHKHDL ('iIniFieldType = IO: ReStart both flow and thermal fields', MYID)
            IF(iIniFieldType == 1)  CALL CHKHDL ('iIniFieldType = IO: ReStart only flow field, no thermal field', MYID)
            IF(iIniFieldTime == 0)  CALL CHKHDL ('iIniFieldTime = IO: ReStart following previous time.', MYID)
            IF(iIniFieldTime == 1)  CALL CHKHDL ('iIniFieldTime = IO: ReStart with re-setting time to zero', MYID)
        END IF
        IF(iIniField_io == IniField_random) THEN
            CALL CHKRLHDL  ('VPERG (initial velocity perturbation)  = ', MYID, VPERG)
            CALL CHKRLHDL  ('SVPERG (initial velocity perturbuatio near 1/4 Y-wall) = ', MYID, SVPERG)
        END IF
    END IF

    !========== Reading Section [numerics]===================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[numerics]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, NTSTF
    READ(INI, *) SKEY, TSTOP
    READ(INI, *) SKEY, timeFlowStart
    READ(INI, *) SKEY, timeThermoStart
    READ(INI, *) SKEY, DT
    READ(INI, *) SKEY, DTMIN
    READ(INI, *) SKEY, CFLGV

    CALL CHKINTHDL ('Iteration No. for STOP = ', MYID, NTSTF)
    CALL CHKRLHDL  ('Physical Time for STOP = ', MYID, TStop)
    CALL CHKRLHDL  ('Time_thermosolver_Start = ', MYID, timeThermoStart)
    CALL CHKRLHDL  ('Time_flowsolver_Start = ', MYID, timeFlowStart)
    CALL CHKRLHDL  ('Time Step = ', MYID, DT)
    CALL CHKRLHDL  ('Minimum Time Step = ', MYID, DTMIN)
    CALL CHKRLHDL  ('CFL No. = ', MYID, CFLGV)

    !========== Reading Section [methods]===================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[methods]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, iVisScheme
    READ(INI, *) SKEY, iWeightedPre

    IF(iVisScheme == VisImplicit) CALL CHKHDL ('iVisScheme = Implicit Viscous Discretization', MYID)
    IF(iVisScheme == VisExplicit) CALL CHKHDL ('iVisScheme = Explicit Viscous Discretization', MYID)
    IF(iWeightedPre == 1) CALL CHKHDL ('iWeightedPre = Weighted Pressure', MYID)
    IF(iWeightedPre == 0) CALL CHKHDL ('iWeightedPre = Not Weighted Pressure', MYID)

    !========== Readin Section [statistics]===============
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[statistics]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, dtSave1
    READ(INI, *) SKEY, tRunAve1
    READ(INI, *) SKEY, tRunAve_Reset
    READ(INI, *) SKEY, dtAveView
    READ(INI, *) SKEY, dtRawView
    READ(INI, *) SKEY, iterMonitor
    READ(INI, *) SKEY, MGRID
    READ(INI, *) SKEY, JINI

    CALL CHKRLHDL  ('Time Interval For Raw Data Saving = ', MYID, dtSave1)
    CALL CHKRLHDL  ('Time From When For Running Average = ', MYID, tRunAve1)
    CALL CHKRLHDL  ('Reset Time From When For Running Average = ', MYID, tRunAve_Reset)
    CALL CHKRLHDL  ('Time Interval For Averaged Data Visualization = ', MYID, dtAveView)
    CALL CHKRLHDL  ('Time Interval For Raw Data Visualization = ', MYID, dtRawView)
    CALL CHKHDL    ('Tecplot Zone Name For Data Visualization =  ' // TRIM(zoneNameView), MYID)
    CALL CHKINTHDL ('Iteration Interval For Data Monitoring = ', MYID, iterMonitor)
    CALL CHKINTHDL ('MGRID = ', MYID, MGRID)
    CALL CHKINTHDL ('JINI = ', MYID, JINI)

    !========== Readin Section [postprocess]=====================
    READ(INI, *) SECT
    DO WHILE((SECT(1:1) == ';') .OR. (SECT(1:1) == '#') .OR. (SECT(1:1) == ' '))
        READ(INI, *) SECT   ! skipping the comments lines.
    END DO
    LENS = LEN_TRIM(SECT)
    CALL CHKHDL('--------------------------------', MYID)
    CALL CHKHDL('Read in ' // SECT(1:LENS), MYID)
    IF(SECT(1:LENS) /= '[postprocess]') &
    CALL ERRHDL(' Reading fails: ' // SECT(1:LENS), MYID)

    READ(INI, *) SKEY, iPostProcess
    READ(INI, *) SKEY, iPPInst
    READ(INI, *) SKEY, iPPSpectra
    READ(INI, *) SKEY, iPPQuadrants
    READ(INI, *) SKEY, iPPDimension
    READ(INI, *) SKEY, pp_instn_sz

    IF(iPostProcess == 0) CALL CHKHDL ('iPostProcess = normal reStart', MYID)
    IF(iPostProcess == 1) CALL CHKHDL ('iPostProcess = postprocessing data only', MYID)
    IF(iPostProcess == 2) CALL CHKHDL ('iPostProcess = postprocessing data from given instantanous flow field', MYID)

    IF(iPPInst == 1) CALL CHKHDL ('iPPInst = postprocessing instantanous isosurface', MYID)
    IF(iPPInst == 0) CALL CHKHDL ('iPPInst = not postprocessing instantanous isosurface', MYID)

    IF(iPPQuadrants == 1) CALL CHKHDL ('iPPQuadrants = postprocessing quadrants', MYID)
    IF(iPPQuadrants == 0) CALL CHKHDL ('iPPQuadrants = not postprocessing quadrants', MYID)

    IF(iPPSpectra == 1) CALL CHKHDL ('iPPSpectra = postprocessing energy spectra and correlations', MYID)
    IF(iPPSpectra == 0) CALL CHKHDL ('iPPSpectra = not postprocessing energy spectra and correlations', MYID)

    IF(iPPDimension == 1) CALL CHKHDL ('iPPDimension = output both dim and undim results', MYID)
    IF(iPPDimension == 0) CALL CHKHDL ('iPPDimension = output only undim results', MYID)

    IF(iPostProcess == 2) THEN
        CALL CHKINTHDL ('Steps of given instantanous flow field  = ', MYID, pp_instn_sz)
        IF(iPostProcess == 2 .AND. pp_instn_sz < 1) THEN
            CALL ERRHDL(' iPPInstansz is less than 1 for postprocessing instantanous flow fields', MYID)
        END IF

        ALLOCATE (pp_instn_tim(pp_instn_sz)); pp_instn_tiM = 0.0_WP
        DO N = 1, pp_instn_sz
            READ(INI, *) pp_instn_tim(N)
            CALL CHKRLHDL  ('pp_instn_tim  = ', MYID, pp_instn_tim(N))
        END DO
    ELSE
        DO N = 1, pp_instn_sz
            READ(INI, *) RTMP
        END DO
    END IF

    CLOSE(INI)
    !================end of reading in data================================

    IF(MOD(MGRID, 2) == 0) MGRID = MGRID + 1

    NND1_tg = NCL1_tg + 1
    NND1_io = NCL1_io + 1
    NND2 = NCL2 + 1
    NND3 = NCL3 + 1

    IF(IoFlowFlg) THEN
        IF(TgFlowFlg) THEN
            NCL1S = 0
            NCL1E = NCL1_io + 1
        ELSE
            NCL1S = 1
            NCL1E = NCL1_io
        END IF
    END IF


    RETURN

END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE BCAST_READINI
    USE init_info
    USE mpi_info
    USE mesh_info
    USE thermal_info
    USE CHT_info
    IMPLICIT NONE

    !================[flowtype]===================
    CALL MPI_BCAST( iCase, 1, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iDomain, 1, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iThermoDynamics, 1, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iCHT,   1, MPI_INTEGER4, 0, ICOMM, IERROR )

    !================[geometry]===================
    CALL MPI_BCAST( HX_tg, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HX_io, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HZ,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HYB,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HYT,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !================[mesh]===================
    CALL MPI_BCAST( NCL1_tg, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NCL1_io, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NCL2,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NCL3,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( STR2,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !================[boundary]===================
    CALL MPI_BCAST(BCX_tg,      2, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST(BCX_io,      2, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST(BCZ,         2, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST(BCY,         2, MPI_INTEGER4,         0, ICOMM, IERROR )

    !================[fluid]===================
    CALL MPI_BCAST( REN,       1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( ReIni,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TLgRe,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iFlowDriven, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( Cf_Given,      1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )


    !================[thermohydraulics]===================
    CALL MPI_BCAST( iFluidMedia,     1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iGravity,     1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( L0,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( T0,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( Ti,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iThermalWallType,  2, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( thermalWallBC_Dim,  2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !================[conjugateHeatTransfer]===================
    CALL MPI_BCAST( NCL2_Solid,  2, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( HY_Solid,    2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HXs_Solid,   2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( HXe_Solid,   2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( Cp_Solid,    2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( D_Solid,     2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( K_Solid,     2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !================[initialisation]===================
    CALL MPI_BCAST( iIniField_tg,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iIniField_io,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iRandomType,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( VPERG,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( SVPERG,       1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TimeReStart_tg,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TimeReStart_io,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( iIniFieldType, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iIniFieldTime, 1, MPI_INTEGER4,      0, ICOMM, IERROR )

    !================[numerics]===================
    CALL MPI_BCAST( NTSTF,        1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( TSTOP,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( timeFlowStart, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( timeThermoStart, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DT,           1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( CFLGV,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DTmin,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    !================[methods]===================
    CALL MPI_BCAST( iWeightedPre, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iVisScheme,   1, MPI_INTEGER4,         0, ICOMM, IERROR )

    !================[statistics]===================
    CALL MPI_BCAST( dtSave1,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( tRunAve1,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( tRunAve_Reset, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( dtAveView,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( dtRawView,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( zoneNameView, 64, MPI_CHARACTER,        0, ICOMM, IERROR )
    CALL MPI_BCAST( iterMonitor,      1, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( MGRID,      1, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( JINI,       1, MPI_INTEGER4, 0, ICOMM, IERROR )

    !================[postprocess]===================
    CALL MPI_BCAST( iPostProcess, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iPPInst,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iPPSpectra,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iPPQuadrants,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( iPPDimension,        1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( pp_instn_sz,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
    IF(iPostProcess == 2) THEN
        IF(MYID /= 0) THEN
            allocate(pp_instn_tim(pp_instn_sz)); pp_instn_tim = 0.0_WP
        END IF
        CALL MPI_BCAST( pp_instn_tim, pp_instn_sz, MPI_DOUBLE_PRECISION,         0, ICOMM, IERROR )
    END IF

    !================ All otherS =================================================
    CALL MPI_BCAST( IoFlowFlg,  1, MPI_LOGICAL,          0, ICOMM, IERROR )
    CALL MPI_BCAST( TgFlowFlg,  1, MPI_LOGICAL,          0, ICOMM, IERROR )
    CALL MPI_BCAST( iThermoProperty,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NISTFLNM,   64, MPI_CHARACTER,        0, ICOMM, IERROR )
    CALL MPI_BCAST( NND1_tg,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NND1_io,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NND2,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NND3,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NCL1S,      1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( NCL1E,      1, MPI_INTEGER4,         0, ICOMM, IERROR )

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE CHECK_BC(BC_TYPE, STR)
    USE mpi_info
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: BC_TYPE
    CHARACTER(*), INTENT(IN) :: STR

    SELECT CASE(BC_TYPE)
    CASE(1)
        CALL CHKHDL(TRIM(STR) // ' BC = Dirichlet (specified value)', MYID)
    CASE(2)
        CALL CHKHDL(TRIM(STR) // ' BC = Neumann (specified the normal derivation)', MYID)
    CASE(3)
        CALL CHKHDL(TRIM(STR) // ' BC = Periodic', MYID)
    CASE DEFAULT
        CALL ERRHDL('No Such BC Defined For' // TRIM(STR), MYID)
    END SELECT
    RETURN
END SUBROUTINE
