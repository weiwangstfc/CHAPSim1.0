!**********************************************************************************************************************************
!> @brief
!>        to prepARe PARAMETERs for thermo- dynamics calculation
!> @details
!> SUBROUTINE: thermal_init (in MYID = 0)
!> SUBROUTINE: thermal_gravity (in MYID = 0)
!> SUBROUTINE: thermal_wall_Boundaries_noCHT (in MYID = 0)
!> SUBROUTINE: thermal_init_write (in MYID = 0)
!> SUBROUTINE: BCAST_THERMAL_INIT (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
! 09/ 2020- Added more fluid types and optimized, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE thermal_init ! master only
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: spline_interpolation_TH
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: thermophysical_function_TD
    REAL(WP) :: thermophysical_function_TK
    REAL(WP) :: thermophysical_function_TM
    REAL(WP) :: thermophysical_function_TCp
    REAL(WP) :: thermophysical_function_TB
    REAL(WP) :: thermophysical_function_TH
    INTEGER(4) :: dim

    ! =============build up ref and thermo- PRoperty (undim) relations ============
    IF(iThermoProperty == search_table) THEN
        !(get ref and undim table)
        CALL Property_Table_Nondimensionalization

        !build up relation of DH- H, H~other, T~others
        CALL SPLINE_COEFFICIENTS_FOR_DH
        CALL SPLINE_COEFFICIENTS_FOR_H
        CALL SPLINE_COEFFICIENTS_FOR_T

        ! get the dimensionless inlet Temperature
        T_inlet = Ti / T0
        H_inlet = spline_interpolation_TH(T_inlet)
        D_inlet = spline_interpolation_HD(H_inlet)
        K_inlet = spline_interpolation_HK(H_inlet)
        M_inlet = spline_interpolation_HM(H_inlet)
        CP_inlet = spline_interpolation_HCp(H_inlet)
        DH_inlet = D_inlet * H_inlet

    ELSE IF (iThermoProperty == properties_functions) THEN
        !== get ref values (dimensional) =====
        dim = 1
        D0  = thermophysical_function_TD(T0, dim)
        K0  = thermophysical_function_TK(T0, dim)
        M0  = thermophysical_function_TM(T0, dim)
        CP0 = thermophysical_function_TCp(T0, dim)
        B0  = thermophysical_function_TB(T0, dim)
        H0  = thermophysical_function_TH(T0, dim)

        ! DH- T, DH- H relations are always using table -searching (undim only)...
        CALL thermophysical_function_DH_preparation

        !==undim inlet info =============
        dim = 0
        T_inlet = Ti / T0
        D_inlet = thermophysical_function_TD(T_inlet, dim)
        K_inlet = thermophysical_function_TK(T_inlet, dim)
        M_inlet = thermophysical_function_TM(T_inlet, dim)
        CP_inlet = thermophysical_function_TCp(T_inlet, dim)
        H_inlet = thermophysical_function_TH(T_inlet, dim)
        DH_inlet = D_inlet * H_inlet

    ELSE
    END IF

    PRT0 = M0 * CP0 / K0
    CTHECD = 1.0_WP / REN / PRT0
    U0   = REN * M0 / D0 / L0
    RHOU20 = D0 * U0 * U0

    !========gravitY ==============
    CALL thermal_gravity

    !======= Applied wall b.c.==================
    CALL thermal_wall_Boundaries_noCHT

    !=======check PARAMETERS ==================
    CALL thermal_init_write

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE thermal_gravity
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    IBuoF = 0
    IF(iGravity == 0) THEN       ! no gravity considered
        G_A = 0.0_WP
        F_A = 0.0_WP
    ELSE IF(iGravity == 1) THEN  ! vertical downwards flow, buoyancy = x direction
        G_A = 9.80665_WP
        F_A = L0/ U0/ U0 * G_A
        IBuoF(1) = 1
    ELSE IF(iGravity == -1) THEN  ! vertical upwards flow, buoyancy = x direction
        G_A = 9.80665_WP
        F_A = -L0/ U0/ U0 * G_A
        IBuoF(1) = 1
    ELSE IF(iGravity == -2) THEN  ! horizontal flow, buoyancy = y direction
        G_A = 9.80665_WP
        F_A = -L0/ U0/ U0 * G_A
        IBuoF(2) = 1
    ELSE IF(iGravity == 3) THEN  ! horizontal flow, should not be USEd.
        G_A = 9.80665_WP
        F_A = L0/ U0/ U0 * G_A
        IBuoF(3) = 1
    ELSE
        G_A = 0.0_WP
        F_A = 0.0_WP
    END IF
    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE thermal_wall_Boundaries_noCHT
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    INTEGER(4) :: NIN1, NIN2, I, K

    !===================== Isoflux ========================
    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .OR. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN

        ALLOCATE (WALLFLUX(NCL1S : NCL1E, 2, 1 : NCL3)); WALLFLUX = 0.0_WP

        !cARefully in ini file
        ! bottom wall heating (heat flux in ), qw = - KDT/ Dy = positive
        ! bottom wall cooling (heat flux out), qw = - KDT/ Dy = negtive
        ! top    wall heating (heat flux in ), qw = - KDT/ Dy = negtive
        ! top    wall cooling (heat flux out), qw = - KDT/ Dy = positive

        !===== TOP WALL CONSTANT HEAT FLUX====
        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN

            thermalWall_nondim(iTopWall) = thermalWallBC_Dim(iTopWall) * L0 / K0 / T0
            NIN2 = 0
            DO I = NCL1S, NCL1E
                !===== TOP WALL HEATING====
                DO K = 1, NCL3
                    IF(I < NIN2) THEN
                        WALLFLUX(I, iTopWall, K) = 0.0_WP
                    ELSE
                        WALLFLUX(I, iTopWall, K) = -thermalWall_nondim(iTopWall)
                    END IF
                END DO
            END DO
        END IF

        !=====BOT WALL  CONSTANT HEAT FLUX====
        IF(iCase /= Ipipec .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN

            thermalWall_nondim(iBotWall) = thermalWallBC_Dim(iBotWall) * L0/ K0/ T0
            NIN1 = 0
            DO I = NCL1S, NCL1E
                DO K = 1, NCL3
                    IF(I < NIN1) THEN
                        WALLFLUX(I, iBotWall, K) = 0.0_WP
                    ELSE
                        WALLFLUX(I, iBotWall, K) = thermalWall_nondim(iBotWall)
                    END IF
                END DO
            END DO
        END IF

        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
            Gr0(1:2) = G_A * B0 * thermalWallBC_Dim(1:2) * (L0**4) / K0 / ((M0 / D0)**2)
            BO0(1, 1:2) = Gr0(1:2) / REN**(2.7_WP)
            BO0(2, 1:2) = Gr0(1:2) / (REN**3.425_WP) / (PRT0**0.8_WP)
        END IF
    END IF


    !===================== IsothermaL ========================
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .OR. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        ALLOCATE (T_WAL_GV(NCL1S : NCL1E, 2)  ); T_WAL_GV = 0.0_WP
        ALLOCATE (H_WAL_GV(NCL1S : NCL1E, 2)  ); H_WAL_GV = 0.0_WP
        ALLOCATE (D_WAL_GV(NCL1S : NCL1E, 2)  ); D_WAL_GV = 0.0_WP
        ALLOCATE (M_WAL_GV(NCL1S : NCL1E, 2)  ); M_WAL_GV = 0.0_WP
        ALLOCATE (K_WAL_GV(NCL1S : NCL1E, 2)  ); K_WAL_GV = 0.0_WP
        ALLOCATE (Cp_WAL_GV(NCL1S : NCL1E, 2) ); Cp_WAL_GV = 0.0_WP
        ALLOCATE (B_WAL_GV(NCL1S : NCL1E, 2) ); B_WAL_GV = 0.0_WP
        ALLOCATE (DH_WAL_GV(NCL1S : NCL1E, 2) ); DH_WAL_GV = 0.0_WP
        ! top wall
        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature) THEN
            thermalWall_nondim(iTopWall) = thermalWallBC_Dim(iTopWall) / T0
            NIN2 = 0
            DO I = NCL1S, NCL1E
                IF(I < NIN2) THEN
                    T_WAL_GV(I, iTopWall) = thermalWall_nondim(iTopWall)
                ELSE
                    T_WAL_GV(I, iTopWall) = thermalWall_nondim(iTopWall)
                END IF

                CALL THERM_PROP_UPDATE_FROM_T(DH_WAL_GV(I, iTopWall), H_WAL_GV(I, iTopWall), T_WAL_GV(I, iTopWall), &
                     D_WAL_GV(I, iTopWall), M_WAL_GV(I, iTopWall), K_WAL_GV(I, iTopWall), Cp_WAL_GV(I, iTopWall), &
                     B_WAL_GV(I, iTopWall))
            END DO
        END IF
        ! bottom wall
        IF(iCase /= Ipipec .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            thermalWall_nondim(iBotWall) = thermalWallBC_Dim(iBotWall) / T0
            NIN1 = 0
            DO I = NCL1S, NCL1E
                IF(I < NIN1) THEN
                    T_WAL_GV(I, iBotWall) = thermalWall_nondim(iBotWall)
                ELSE
                    T_WAL_GV(I, iBotWall) = thermalWall_nondim(iBotWall)
                END IF

                CALL THERM_PROP_UPDATE_FROM_T(DH_WAL_GV(I, iBotWall), H_WAL_GV(I, iBotWall), T_WAL_GV(I, iBotWall), &
                     D_WAL_GV(I, iBotWall), M_WAL_GV(I, iBotWall), K_WAL_GV(I, iBotWall), Cp_WAL_GV(I, iBotWall), &
                     B_WAL_GV(I, iBotWall))
            END DO
        END IF

        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            Gr0(1:2) = G_A*B0 * DABS(thermalWallBC_Dim(1) -thermalWallBC_Dim(2)) * ((2.0_WP * L0)**3) / ((M0/ D0)**2)
            BO0(1, 1:2) = Gr0(1:2) / (2.0_WP * REN)**(2.7_WP)
            BO0(2, 1:2) = Gr0(1:2) / ((2.0_WP * REN)**3.425_WP) / (PRT0**0.8_WP)
        END IF
    END IF

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE thermal_init_write
    USE THERMAL_INFO
    USE MESH_INFO
    USE INIT_INFO
    IMPLICIT NONE
    REAL(WP) :: massflux
    INTEGER(4) :: I
    INTEGER(4) :: N
    REAL(WP) :: DH_tmp, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp

    massflux = REN * M0 / L0

    IF(MYID /= 0) RETURN

    IF(iThermoProperty == search_table) THEN
        OPEN(10, FILE = TRIM(FilePath0) // 'CHK_LIST_TABLE_NONDIM.dat')
        WRITE(10, '(A, 1ES13.5)') '# REF T0(K) =             ', T0
        WRITE(10, '(A, 1ES13.5)') '# REF H0(J) =             ', H0
        WRITE(10, '(A, 1ES13.5)') '# REF D0(Kg/ M3) = ', D0
        WRITE(10, '(A, 1ES13.5)') '# REF M0(Pa-s) =          ', M0
        WRITE(10, '(A, 1ES13.5)') '# REF K0(W/ M - K) = ', K0
        WRITE(10, '(A, 1ES13.5)') '# REF B0(1 / K) =           ', B0
        WRITE(10, '(A, 1ES13.5)') '# REF massflux(Kg/ M2s) =  ', massflux
        WRITE(10, '(A, 1ES13.5)') '# REF CP0 =              ', CP0
        WRITE(10, '(A, 1ES13.5)') '# REF H0/ CP0/ T0 =        ', H0/ CP0/ T0

        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            WRITE(10, '(A, A, A)')     '#         ', 'BOTTOM WALL', ' TOP WALL'
            WRITE(10, '(A, 2ES13.5)') '# T =      ', T_WAL_GV(1, 1), T_WAL_GV(1, 2)
            WRITE(10, '(A, 2ES13.5)') '# H =      ', H_WAL_GV(1, 1), H_WAL_GV(1, 2)
            WRITE(10, '(A, 2ES13.5)') '# D=       ', D_WAL_GV(1, 1), D_WAL_GV(1, 2)
            WRITE(10, '(A, 2ES13.5)') '# M =      ', M_WAL_GV(1, 1), M_WAL_GV(1, 2)
            WRITE(10, '(A, 2ES13.5)') '# K =      ', K_WAL_GV(1, 1), K_WAL_GV(1, 2)
        END IF

        WRITE(10, '(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA'
        WRITE(10, *) '# ', N_LIST

        DO I = 1, N_LIST
            WRITE(10, '(8ES13.5)') P0, &
            LIST_H(I), LIST_T(I), LIST_D(I), &
            LIST_M(I), LIST_K(I), LIST_CP(I), LIST_B(I)
        END DO
        CLOSE(10)
    END IF

    IF(iThermoProperty == properties_functions)  THEN
        N = 512
        OPEN(10, FILE = TRIM(FilePath0) // 'CHK_properties_functions_NONDIM_T.dat')
        WRITE(10, '(A)')      '#1DH    2H     3T    4D    5M     6K      7CP      8BETA'
        WRITE(10, '(A, I3.1)') '# ', N
        DO I = 1, N
            SELECT CASE (iFluidMedia)
                CASE (iLiquidSodium)
                    T_tmp = Tm0_Na + (Tb0_Na - Tm0_Na) / DBLE(N) * DBLE(I)
                CASE (iLiquidLead)
                    T_tmp = Tm0_Pb + (Tb0_Pb - Tm0_Pb) / DBLE(N) * DBLE(I)
                CASE (iLiquidBismuth)
                    T_tmp = Tm0_Bi + (Tb0_Bi - Tm0_Bi) / DBLE(N) * DBLE(I)
                CASE (iLiquidLBE)
                    T_tmp = Tm0_LBE + (Tb0_LBE - Tm0_LBE) / DBLE(N) * DBLE(I)
                CASE DEFAULT
                    T_tmp = Tm0_Na + (Tb0_Na - Tm0_Na) / DBLE(N) * DBLE(I)
            END SELECT

            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T_tmp / T0, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)

            WRITE(10, '(8ES13.5)') DH_tmp, H_tmp, T_tmp / T0, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp
        END DO
        CLOSE(10)
    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        CALL CHKHDL    ('Wall States (dimensionaless) at (Y = -1) & Y = (+ 1)', MYID)
        CALL CHK2RLHDL ('Wall T = ', MYID, T_WAL_GV(1, 1), T_WAL_GV(1, 2))
        CALL CHK2RLHDL ('Wall H = ', MYID, H_WAL_GV(1, 1), H_WAL_GV(1, 2))
        CALL CHK2RLHDL ('Wall D = ', MYID, D_WAL_GV(1, 1), D_WAL_GV(1, 2))
        CALL CHK2RLHDL ('Wall M = ', MYID, M_WAL_GV(1, 1), M_WAL_GV(1, 2))
        CALL CHK2RLHDL ('Wall K = ', MYID, K_WAL_GV(1, 1), K_WAL_GV(1, 2))
    END IF

    !========WRITE ===========
    CALL CHKHDL    ('Reference States (dimensional)', MYID)
    CALL CHKRLHDL  ('REF P0 (Pa) = ', MYID, P0)
    CALL CHKRLHDL  ('REF RHOU20 (Pa) = ', MYID, RHOU20)
    CALL CHKRLHDL  ('REF L0 (m) = ', MYID, L0)
    CALL CHKRLHDL  ('REF U0 (m/s) = ', MYID, U0)
    CALL CHKRLHDL  ('REF DU (Kg/m2s) = ', MYID, U0 * D0)
    CALL CHKRLHDL  ('REF Mdot (Kg/m2s) = ', MYID, Massflux)
    CALL CHKRLHDL  ('REF T0 (K) = ', MYID, T0)
    CALL CHKRLHDL  ('REF H0 (J) = ', MYID, H0)
    CALL CHKRLHDL  ('REF D0 (Kg/M3) = ', MYID, D0)
    CALL CHKRLHDL  ('REF Cp (J/Kg/K) = ', MYID, CP0)
    CALL CHKRLHDL  ('REF K0 (W/M-K) = ', MYID, K0)
    CALL CHKRLHDL  ('REF M0 (Pa-s) = ', MYID, M0)
    CALL CHKRLHDL  ('REF B0 (1/K) = ', MYID, B0)

    CALL CHKRLHDL  ('Initial T_inlet = ', MYID, T_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial H_inlet = ', MYID, H_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial D_inlet = ', MYID, D_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial K_inlet = ', MYID, K_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial M_inlet = ', MYID, M_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial CP_inlet = ', MYID, CP_inlet)    !updated by Junjie, 2017/03 / 26
    CALL CHKRLHDL  ('Initial DH_inlet = ', MYID, DH_inlet) !updated by Junjie, 2017/03 / 26

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
        CALL CHKRLHDL  ('REF WALHFLX TOPWALL (W/M2) = ', MYID, thermalWallBC_Dim(1))
        CALL CHKRLHDL  ('REF WALHFLX BOTWALL (W/M2) = ', MYID, thermalWallBC_Dim(2))
    END IF
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        CALL CHKRLHDL  ('REF T1 TOPWALL (K) = ', MYID, thermalWallBC_Dim(1))
        CALL CHKRLHDL  ('REF T2 BOTWALL (K) = ', MYID, thermalWallBC_Dim(2))
    END IF


    CALL CHKHDL    ('Reference States (dimensionless)', MYID)
    CALL CHKRLHDL  ('REF CP0 = ', MYID, CP0)

    CALL CHKRLHDL  ('REF RE = ', MYID, REN)
    CALL CHKRLHDL  ('REF PRT0 = ', MYID, PRT0)
    CALL CHKRLHDL  ('REF 1/Fr2 = ', MYID, F_A)

    !CALL CHKRLHDL  ('Wall HEAT (BOT) FROM X / L0 = ', MYID,XND_io(NIN1 + 1))
    !CALL CHKRLHDL  ('Wall HEAT (TOP) FROM X / L0 = ', MYID,XND_io(NIN2 + 1))

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
        CALL CHKRLHDL  ('REF Gr0(Bottom Wall) = ', MYID, Gr0(1))
        CALL CHKRLHDL  ('REF Gr0(Top Wall) = ', MYID, Gr0(2))

        CALL CHKRLHDL  ('REF Gr0/ RE^2.7(Bottom Wall) = ', MYID, BO0(1, 1))
        CALL CHKRLHDL  ('REF Gr0/ RE^2.7(Top Wall) = ', MYID, BO0(1, 2))

        CALL CHKRLHDL  ('REF Gr0/ RE^3.425/Pr^0.8(Bottom Wall) = ', MYID, BO0(2, 1))
        CALL CHKRLHDL  ('REF Gr0/ RE^3.425/Pr^0.8(Top Wall) = ', MYID, BO0(2, 2))

        CALL CHKRLHDL  ('REF Q+ (WHFLUX0) (Bottom Wall) = ', MYID, thermalWall_nondim(1))
        CALL CHKRLHDL  ('REF Q+ (WHFLUX0) (Top Wall) = ', MYID, thermalWall_nondim(2))
    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        CALL CHKRLHDL  ('REF Gr0(estimated Bottom Wall) = ', MYID, Gr0(1))
        CALL CHKRLHDL  ('REF Gr0(estimated Top Wall) = ', MYID, Gr0(2))

        CALL CHKRLHDL  ('REF Gr0/ RE^2.7(Bottom Wall) = ', MYID, BO0(1, 1))
        CALL CHKRLHDL  ('REF Gr0/ RE^2.7(Top Wall) = ', MYID, BO0(1, 2))

        CALL CHKRLHDL  ('REF Gr0/ RE^3.425/Pr^0.8(Bottom Wall) = ', MYID, BO0(2, 1))
        CALL CHKRLHDL  ('REF Gr0/ RE^3.425/Pr^0.8(Top Wall) = ', MYID, BO0(2, 2))

        CALL CHKRLHDL  ('REF T + (WALLTEMP) (Bottom Wall) = ', MYID, thermalWall_nondim(1))
        CALL CHKRLHDL  ('REF T + (WALLTEMP) (Top Wall) = ', MYID, thermalWall_nondim(2))
    END IF


END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE BCAST_THERMAL_INIT
    USE thermal_info
    USE mesh_info
    USE flow_info
    USE INIT_INFO
    IMPLICIT NONE

    ! =common variables for tablE -searching and functions =
    CALL MPI_BCAST( H0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( D0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( K0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( M0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( B0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( Cp0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( PRT0,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( CTHECD, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( U0,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( RHOU20, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( G_A,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( F_A,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IBuoF,  3, MPI_INTEGER4,         0, ICOMM, IERROR )

    CALL MPI_BCAST( T_inlet,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( H_inlet,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( D_inlet,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( M_inlet,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( K_inlet,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( CP_inlet,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( DH_inlet, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( thermalWall_nondim, 2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    !CALL MPI_BCAST( WHFLUX_UND,   2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( Gr0,          2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( BO0,          4, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( IMAX_DH,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMIN_DH,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMAX_H,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMIN_H,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMAX_T,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMIN_T,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMAX_Cp,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
    CALL MPI_BCAST( IMIN_Cp,  1, MPI_INTEGER4,         0, ICOMM, IERROR )

    IF(MYID /= 0) THEN
        IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .OR. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
            ALLOCATE (WALLFLUX(NCL1S : NCL1E, 2, NCL3)  )
        END IF

        IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .OR. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
            ALLOCATE (T_WAL_GV(NCL1S : NCL1E, 2)  )
            ALLOCATE (H_WAL_GV(NCL1S : NCL1E, 2)  )
            ALLOCATE (D_WAL_GV(NCL1S : NCL1E, 2)  )
            ALLOCATE (M_WAL_GV(NCL1S : NCL1E, 2)  )
            ALLOCATE (K_WAL_GV(NCL1S : NCL1E, 2)  )
            ALLOCATE (Cp_WAL_GV(NCL1S : NCL1E, 2)  )
        END IF
    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .OR. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
        CALL MPI_BCAST( WALLFLUX, 2 * (NCL1E - NCL1S + 1) * NCL3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    END IF

    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .OR. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        CALL MPI_BCAST( T_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( H_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( D_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( M_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( K_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST(Cp_WAL_GV, 2 * (NCL1E - NCL1S + 1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    END IF

    CALL MPI_BCAST( N_LIST, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    IF(MYID /= 0) THEN
        ALLOCATE( LIST_H (N_LIST) ) ;  LIST_H  = 0.0_WP
        ALLOCATE( LIST_T (N_LIST) ) ;  LIST_T  = 0.0_WP
        ALLOCATE( LIST_DH(N_LIST) ) ;  LIST_DH = 0.0_WP

        ALLOCATE(SplineCoeff_DHH_B(N_LIST)); SplineCoeff_DHH_B = 0.0_WP
        ALLOCATE(SplineCoeff_DHH_C(N_LIST)); SplineCoeff_DHH_C = 0.0_WP
        ALLOCATE(SplineCoeff_DHH_D(N_LIST)); SplineCoeff_DHH_D = 0.0_WP!updated by Junjie, 2017/03 / 26

        ALLOCATE(SplineCoeff_DHT_B(N_LIST)); SplineCoeff_DHT_B = 0.0_WP
        ALLOCATE(SplineCoeff_DHT_C(N_LIST)); SplineCoeff_DHT_C = 0.0_WP
        ALLOCATE(SplineCoeff_DHT_D(N_LIST)); SplineCoeff_DHT_D = 0.0_WP
    END IF

    CALL MPI_BCAST( LIST_H, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( LIST_T, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( LIST_DH, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( SplineCoeff_DHH_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( SplineCoeff_DHH_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( SplineCoeff_DHH_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( SplineCoeff_DHT_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( SplineCoeff_DHT_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( SplineCoeff_DHT_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

   ! = Special variables for tablE -searching =
    IF(iThermoProperty == search_table) THEN

        CALL MPI_BCAST( DHmax,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DHmin,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T4MaxDH, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T4MinDH, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( Cpmax,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T4CpMax,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        IF(MYID /= 0) THEN
            ALLOCATE( LIST_D (N_LIST) ) ;  LIST_D  = 0.0_WP
            ALLOCATE( LIST_M (N_LIST) ) ;  LIST_M  = 0.0_WP
            ALLOCATE( LIST_K (N_LIST) ) ;  LIST_K  = 0.0_WP
            ALLOCATE( LIST_B (N_LIST) ) ;  LIST_B  = 0.0_WP
            ALLOCATE( LIST_CP(N_LIST) ) ;  LIST_CP = 0.0_WP

            ALLOCATE(SplineCoeff_HT_B(N_LIST)); SplineCoeff_HT_B = 0.0_WP
            ALLOCATE(SplineCoeff_HT_C(N_LIST)); SplineCoeff_HT_C = 0.0_WP
            ALLOCATE(SplineCoeff_HT_D(N_LIST)); SplineCoeff_HT_D = 0.0_WP

            ALLOCATE(SplineCoeff_HD_B(N_LIST)); SplineCoeff_HD_B = 0.0_WP
            ALLOCATE(SplineCoeff_HD_C(N_LIST)); SplineCoeff_HD_C = 0.0_WP
            ALLOCATE(SplineCoeff_HD_D(N_LIST)); SplineCoeff_HD_D = 0.0_WP

            ALLOCATE(SplineCoeff_HM_B(N_LIST)); SplineCoeff_HM_B = 0.0_WP
            ALLOCATE(SplineCoeff_HM_C(N_LIST)); SplineCoeff_HM_C = 0.0_WP
            ALLOCATE(SplineCoeff_HM_D(N_LIST)); SplineCoeff_HM_D = 0.0_WP

            ALLOCATE(SplineCoeff_HK_B(N_LIST)); SplineCoeff_HK_B = 0.0_WP
            ALLOCATE(SplineCoeff_HK_C(N_LIST)); SplineCoeff_HK_C = 0.0_WP
            ALLOCATE(SplineCoeff_HK_D(N_LIST)); SplineCoeff_HK_D = 0.0_WP

            ALLOCATE(SplineCoeff_HCp_B(N_LIST)); SplineCoeff_HCp_B = 0.0_WP
            ALLOCATE(SplineCoeff_HCp_C(N_LIST)); SplineCoeff_HCp_C = 0.0_WP
            ALLOCATE(SplineCoeff_HCp_D(N_LIST)); SplineCoeff_HCp_D = 0.0_WP

            ALLOCATE(SplineCoeff_HB_B(N_LIST)); SplineCoeff_HB_B = 0.0_WP
            ALLOCATE(SplineCoeff_HB_C(N_LIST)); SplineCoeff_HB_C = 0.0_WP
            ALLOCATE(SplineCoeff_HB_D(N_LIST)); SplineCoeff_HB_D = 0.0_WP


            ALLOCATE(SplineCoeff_TH_B(N_LIST)); SplineCoeff_TH_B = 0.0_WP
            ALLOCATE(SplineCoeff_TH_C(N_LIST)); SplineCoeff_TH_C = 0.0_WP
            ALLOCATE(SplineCoeff_TH_D(N_LIST)); SplineCoeff_TH_D = 0.0_WP

            ALLOCATE(SplineCoeff_TD_B(N_LIST)); SplineCoeff_TD_B = 0.0_WP
            ALLOCATE(SplineCoeff_TD_C(N_LIST)); SplineCoeff_TD_C = 0.0_WP
            ALLOCATE(SplineCoeff_TD_D(N_LIST)); SplineCoeff_TD_D = 0.0_WP

            ALLOCATE(SplineCoeff_TM_B(N_LIST)); SplineCoeff_TM_B = 0.0_WP
            ALLOCATE(SplineCoeff_TM_C(N_LIST)); SplineCoeff_TM_C = 0.0_WP
            ALLOCATE(SplineCoeff_TM_D(N_LIST)); SplineCoeff_TM_D = 0.0_WP

            ALLOCATE(SplineCoeff_TK_B(N_LIST)); SplineCoeff_TK_B = 0.0_WP
            ALLOCATE(SplineCoeff_TK_C(N_LIST)); SplineCoeff_TK_C = 0.0_WP
            ALLOCATE(SplineCoeff_TK_D(N_LIST)); SplineCoeff_TK_D = 0.0_WP

            ALLOCATE(SplineCoeff_TCp_B(N_LIST)); SplineCoeff_TCp_B = 0.0_WP
            ALLOCATE(SplineCoeff_TCp_C(N_LIST)); SplineCoeff_TCp_C = 0.0_WP
            ALLOCATE(SplineCoeff_TCp_D(N_LIST)); SplineCoeff_TCp_D = 0.0_WP

            ALLOCATE(SplineCoeff_TB_B(N_LIST)); SplineCoeff_TB_B = 0.0_WP
            ALLOCATE(SplineCoeff_TB_C(N_LIST)); SplineCoeff_TB_C = 0.0_WP
            ALLOCATE(SplineCoeff_TB_D(N_LIST)); SplineCoeff_TB_D = 0.0_WP

        END IF

        CALL MPI_BCAST( LIST_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( LIST_M, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( LIST_K, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( LIST_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( LIST_CP, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )


        CALL MPI_BCAST( SplineCoeff_HT_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HT_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HT_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_HD_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HD_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HD_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_HM_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HM_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HM_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_HK_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HK_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HK_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_HB_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HB_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HB_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_HCp_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HCp_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_HCp_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )


        CALL MPI_BCAST( SplineCoeff_TH_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TH_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TH_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_TD_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TD_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TD_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_TM_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TM_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TM_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_TK_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TK_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TK_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_TB_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TB_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TB_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

        CALL MPI_BCAST( SplineCoeff_TCp_B, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TCp_C, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SplineCoeff_TCp_D, N_LIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    END IF


    RETURN
END SUBROUTINE
