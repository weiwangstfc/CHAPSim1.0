!**********************************************************************************************************************************
!> @brief
!>        to updata thermo-properties
!> @details
!> SUBROUTINE: THERM_PROP_UPDATE_FROM_DH (in MYID = all)
!> SUBROUTINE: THERM_PROP_UPDATE_FROM_T (in MYID = all)
!> SUBROUTINE: THERM_PROP_UPDATE (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
! 09/ 2020- Added more fluid types and optimized, by Wei Wang (wei.wang@stfc.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE THERM_PROP_UPDATE_FROM_DH(DH_tmp, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: DH_tmp
    REAL(WP), INTENT(OUT) :: H_tmp
    REAL(WP), INTENT(OUT) :: T_tmp
    REAL(WP), INTENT(OUT) :: D_tmp
    REAL(WP), INTENT(OUT) :: M_tmp
    REAL(WP), INTENT(OUT) :: K_tmp
    REAL(WP), INTENT(OUT) :: Cp_tmp
    REAL(WP), INTENT(OUT) :: B_tmp

    REAL(WP) :: spline_interpolation_DHT
    REAL(WP) :: spline_interpolation_DHH
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: spline_interpolation_HB

    REAL(WP) :: thermophysical_function_TD
    REAL(WP) :: thermophysical_function_TK
    REAL(WP) :: thermophysical_function_TM
    REAL(WP) :: thermophysical_function_TCp
    REAL(WP) :: thermophysical_function_TB
    REAL(WP) :: thermophysical_function_TH

    INTEGER(4) :: dim = 0

    IF(iThermoProperty == search_table) THEN
        H_tmp = spline_interpolation_DHH( DH_tmp )
        T_tmp = spline_interpolation_HT( H_tmp )
        D_tmp = spline_interpolation_HM( H_tmp )
        M_tmp = spline_interpolation_HK( H_tmp )
        K_tmp = spline_interpolation_HD( H_tmp )
        Cp_tmp = spline_interpolation_HCp( H_tmp )
        B_tmp = spline_interpolation_HB( H_tmp )
    ELSE IF(iThermoProperty == properties_functions) THEN
        T_tmp = spline_interpolation_DHT( DH_tmp )
        CALL thermophysical_function_Check_Trange(T_tmp, 0)
        H_tmp = thermophysical_function_TH( T_tmp, dim )
        D_tmp = thermophysical_function_TD( T_tmp, dim )
        M_tmp = thermophysical_function_TM( T_tmp, dim )
        K_tmp = thermophysical_function_TK( T_tmp, dim )
        Cp_tmp = thermophysical_function_TCp( T_tmp, dim )
        B_tmp = thermophysical_function_TB( T_tmp, dim )
    ELSE
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
    USE THERMAL_INFO
    USE INIT_INFO
    IMPLICIT NONE

    REAL(WP), INTENT(IN) :: T_tmp
    REAL(WP), INTENT(OUT) :: DH_tmp
    REAL(WP), INTENT(OUT) :: H_tmp
    REAL(WP), INTENT(OUT) :: D_tmp
    REAL(WP), INTENT(OUT) :: M_tmp
    REAL(WP), INTENT(OUT) :: K_tmp
    REAL(WP), INTENT(OUT) :: Cp_tmp
    REAL(WP), INTENT(OUT) :: B_tmp

    REAL(WP) :: spline_interpolation_TH
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: spline_interpolation_HB

    REAL(WP) :: thermophysical_function_TD
    REAL(WP) :: thermophysical_function_TK
    REAL(WP) :: thermophysical_function_TM
    REAL(WP) :: thermophysical_function_TCp
    REAL(WP) :: thermophysical_function_TB
    REAL(WP) :: thermophysical_function_TH

    INTEGER(4) :: dim = 0

    IF(iThermoProperty == search_table) THEN
        H_tmp = spline_interpolation_TH( T_tmp )
        D_tmp = spline_interpolation_HM( H_tmp )
        M_tmp = spline_interpolation_HK( H_tmp )
        K_tmp = spline_interpolation_HD( H_tmp )
        Cp_tmp = spline_interpolation_HCp( H_tmp )
        B_tmp = spline_interpolation_HB( H_tmp )
        DH_tmp = D_tmp * H_tmp
    ELSE IF(iThermoProperty == properties_functions) THEN
        CALL thermophysical_function_Check_Trange(T_tmp, 0)
        H_tmp = thermophysical_function_TH( T_tmp, dim )
        D_tmp = thermophysical_function_TD( T_tmp, dim )
        M_tmp = thermophysical_function_TM( T_tmp, dim )
        K_tmp = thermophysical_function_TK( T_tmp, dim )
        Cp_tmp = thermophysical_function_TCp( T_tmp, dim )
        B_tmp = thermophysical_function_TB( T_tmp, dim )
        DH_tmp = D_tmp * H_tmp
    ELSE
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE THERM_PROP_UPDATE(IREGION)
    !>    @NOTE: WHY TO USE "H" TO UPDATE OTHER PROPERTIES, NOT "H* D"
    !>           Answer: FOR "HD", FOR EACH VALUE OF "HD", THERE IS MORE THAN
    !>                   ONE NUMBER FOUND FOR OTHER PROPERTIES.
    !>                   AS "H" AND "D" HAS OPPOPSITE CHANGING TREND AS "T" CHANGES.
    USE init_info
    USE thermal_info
    USE FLOW_INFO
    USE mesh_info
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: IREGION
    INTEGER(4) :: I, J, K
    INTEGER(4) :: I1, I2
    REAL(WP) :: B_tmp

    IF(IREGION == IOULET) THEN  ! outlet
        I1 = NCL1_io + 1
        I2 = NCL1_io + 1
    ELSE IF(IREGION == IINLET) THEN !INLET
        I1 = 0
        I2 = 0
    ELSE IF(IREGION == IMAIND) THEN
        I1 = 1
        I2 = NCL1_io
    ELSE IF(IREGION == IALLDM) THEN
        I1 = NCL1S
        I2 = NCL1E
    ELSE
        I1 = 0
        I2 = NCL1_io + 1
    END IF

    DO J = 1, N2DO(MYID)
        DO I = i1, I2
            DO K = 1, NCL3
                CALL THERM_PROP_UPDATE_FROM_DH(DH(I, J, K), ENTHALPY(I, J, K), TEMPERATURE(I, J, K), &
                DENSITY(I, J, K), Viscousity(I, J, K), THERMCONDT(I, J, K), HEATCAP(I, J, K), B_tmp)
            END DO
        END DO
    END DO


    RETURN
END SUBROUTINE


! SUBROUTINE ENTHALPY_UPDATE(IREGION)
!     USE THERMAL_INFO
!     USE MESH_INFO
!     IMPLICIT NONE
!     INTEGER(4), INTENT(IN) :: IREGION
!     INTEGER(4) :: I, J, K, I1, I2
!
!     IF(IREGION == IOULET) THEN  ! outlet
!         I1 = NCL1_io + 1
!         I2 = NCL1_io + 1
!     ELSE IF(IREGION == IINLET) THEN !INLET
!         I1 = 0
!         I2 = 0
!     ELSE IF(IREGION == IMAIND) THEN
!         I1 = 1
!         I2 = NCL1_io
!     ELSE IF(IREGION == IALLDM) THEN
!         I1 = NCL1S
!         I2 = NCL1E
!     ELSE
!         I1 = 0
!         I2 = NCL1_io + 1
!     END IF
!
!     DO J = 1, N2DO(MYID)
!         DO I = i1, I2
!             DO K = 1, NCL3
!                 ENTHALPY(I, J, K) = DH(I, J, K) / DENSITY(I, J, K)
!             END DO
!         END DO
!     END DO
!     RETURN
! END SUBROUTINE


!=================================================================================
! SUBROUTINE DH_UPDATE(IREGION)
!     USE THERMAL_INFO
!     USE MESH_INFO
!     IMPLICIT NONE
!     INTEGER(4), INTENT(IN) :: IREGION
!     INTEGER(4) :: I, J, K, I1, I2
!
!     IF(IREGION == IOULET) THEN  ! outlet
!         I1 = NCL1_io + 1
!         I2 = NCL1_io + 1
!     ELSE IF(IREGION == IINLET) THEN !INLET
!         I1 = 0
!         I2 = 0
!     ELSE IF(IREGION == IMAIND) THEN
!         I1 = 1
!         I2 = NCL1_io
!     ELSE IF(IREGION == IALLDM) THEN
!         I1 = NCL1S
!         I2 = NCL1E
!     ELSE
!         I1 = 0
!         I2 = NCL1_io + 1
!     END IF
!
!     DO J = 1, N2DO(MYID)
!         DO I = i1, I2
!             DO K = 1, NCL3
!                 DH(I, J, K) = ENTHALPY(I, J, K) * DENSITY(I, J, K)
!             END DO
!         END DO
!     END DO
!     RETURN
! END SUBROUTINE
