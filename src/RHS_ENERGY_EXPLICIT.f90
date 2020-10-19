!**********************************************************************************************************************************
!> @brief
!>
!> @details
!> SUBROUTINE: RHS_ENERGY_EXPLICIT (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE RHS_ENERGY_EXPLICIT
    !>   @note
    !>    1) the calculation of heat flux is based on a constant velocity field.
    !>    2) \rho h u is treated as (\rho u)h, not (\rho h)u
    !>    3) wall bc is not included, which will be dealed with separately.

    USE thermal_info
    USE mesh_info
    USE flow_info
    USE INIT_INFO
    IMPLICIT NONE

    INTEGER(4) :: IC, IP, IM
    INTEGER(4) :: KC, KP, KM
    INTEGER(4) :: JC, JP, JM, JJ, JJP, JS, JE
    REAL(WP) :: DX_DHU, DY_DHV, DZ_DHW
    REAL(WP) :: DX_KDT, DY_KDT, DZ_KDT, KDTy1, KDTy2
    REAL(WP) :: COE1, COE2, COE3, COE11, COE22, COE221, COE222, COE33, COE331

    !LOGICAL,external :: ISNAN1, ISinf1

    RHS_ENERGY = 0.0_WP

    COE1 = DXI  * XND2CL * 0.5_WP
    COE11 = DXQI * XND2CL * CTHECD
    DO JC = 1, N2DO(MYID)
        JP = JLPV(JC)
        JM = JLMV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)

        !=============coefficientS =============================
        !COEFFICIENT FOR CONVECTIONS
        COE2  = DYFI(JJ)     * RCCI1(JJ) * 0.5_WP
        COE3  = DZI * ZND2CL * RCCI2(JJ) * 0.5_WP

        !COEFFICIENTS FOR VISCOUS TERMS
        COE33 = DZQI * ZND2CL * RCCI2(JJ) * CTHECD
        COE22 = DYFI(JJ) * RCCI1(JJ) * CTHECD

        IF(JJ == 1) THEN
            COE221 = DYCI(JJP) / RNDI1(JJP)
            IF(iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux    ) THEN
                IF(iCase == iPIPEC) THEN
                    COE222 = 0.0_WP
                ELSE
                    COE222 = 1.0_WP / RNDI1(JJ)
                END IF
            END IF

            IF(iThermalWallType(iBotWall) == BC_Fixed_Temperature    ) THEN
                IF(iCase == iPIPEC) THEN
                    COE222 = 0.0_WP
                ELSE
                    COE222 = DYCI(JJ) / RNDI1(JJ)
                END IF
            END IF

        ELSE IF(JJ == NCL2) THEN
            COE222 = DYCI(JJ) / RNDI1(JJ)
            IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux    ) THEN
                COE221 = 1.0_WP / RNDI1(JJP)
            END IF
            IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature    ) THEN
                COE221 = DYCI(JJP) / RNDI1(JJP)
            END IF
        ELSE
            COE221 = DYCI(JJP) / RNDI1(JJP)
            COE222 = DYCI(JJ)  / RNDI1(JJ)
        END IF
        !=============coefficientS =============================

        DO IC = 1, NCL1_io
            IP = IPV_io(IC)
            IM = IMV_io(IC)
            DO KC = 1, NCL3
                KP = KPV(KC)
                KM = KMV(KC)

                !=====================CONVECTION ===================
                ! $ \Delta_x [ (h)^x \cDOt G_X ] $
                ! (I, J, K) = (I'+ 1, J, K)  - (I', J, K)   DX
                DX_DHU = ( ENTHALPY(IP, JC, KC) + ENTHALPY(IC, JC, KC) ) * ( G_io(IP, JC, KC, 1) + G_io(IP, JC, KC, 1) ) - &
                ( ENTHALPY(IC, JC, KC) + ENTHALPY(IM, JC, KC) ) * ( G_io(IC, JC, KC, 1) + G_io(IC, JC, KC, 1) )
                DX_DHU =  DX_DHU * COE1  ! * XND2CL * DXI * 0.5_WP


                ! $ \Delta_y [ (h)^y \cDOt G_Y ] $
                ! (I, J, K) = (I, J'+ 1, K)  - (I, J', K)   DY
                DY_DHV = ( YCL2ND_WFF(JJP) * ENTHALPY(IC, JP, KC) +                         &
                YCL2ND_WFB(JJP) * ENTHALPY(IC, JC, KC) ) * ( G_io(IC, JP, KC, 2) + G_io(IC, JP, KC, 2) ) -    &
                ( YCL2ND_WFF(JJ) * ENTHALPY(IC, JC, KC) +                         &
                YCL2ND_WFB(JJ) * ENTHALPY(IC, JM, KC) ) * ( G_io(IC, JC, KC, 2) + G_io(IC, JC, KC, 2) )
                DY_DHV = DY_DHV * COE2 !* RCCI1(JJ) * DYFI(JJ) * 0.5_WP

                ! $ \Delta_z [ (\rho h)^z \cDOt w ] $
                ! (I, J, K) = (I, J, K'+ 1)  - (I, J, K')   DZ
                DZ_DHW = ( ENTHALPY(IC, JC, KP) + ENTHALPY(IC, JC, KC) ) * ( G_io(IC, JC, KP, 3) + G_io(IC, JC, KP, 3) ) - &
                ( ENTHALPY(IC, JC, KC) + ENTHALPY(IC, JC, KM) ) * ( G_io(IC, JC, KC, 3) + G_io(IC, JC, KC, 3) )
                DZ_DHW = DZ_DHW * COE3 !* ZND2CL * DZI * RCCI2(JJ) * 0.5_WP


                !==========CONDUCTION =================================
                ! $ \Delta_x [k^x (\Delta_x T)] $
                ! (I, J, K) = (I'+ 1, J, K)  - (I', J, K)   DX
                DX_KDT = ( THERMCONDT(IP, JC, KC)  + THERMCONDT(IC, JC, KC) ) * &
                         ( TEMPERATURE(IP, JC, KC) - TEMPERATURE(IC, JC, KC) ) - &
                         ( THERMCONDT(IC, JC, KC)  + THERMCONDT(IM, JC, KC) ) * &
                         ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IM, JC, KC) )
                DX_KDT = DX_KDT * COE11 !XND2CL * DXI * DXI * CTHECD

                ! $ \Delta_y [k^y (\Delta_y T)] $
                ! (I, J, K) = (I, J'+ 1, K)  - (I, J', K)   DY
                KDTy1 = ( YCL2ND_WFF(JJP) * THERMCONDT(IC, JP, KC) +   &
                YCL2ND_WFB(JJP) * THERMCONDT(IC, JC, KC) ) * &
                ( TEMPERATURE(IC, JP, KC) - TEMPERATURE(IC, JC, KC) ) * COE221        !DYCI(JJP) / RNDI(JJP)
                KDTy2 = ( YCL2ND_WFF(JJ) * THERMCONDT(IC, JC, KC) +   &
                YCL2ND_WFB(JJ) * THERMCONDT(IC, JM, KC) ) * &
                ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IC, JM, KC) ) * COE222        !DYCI(JJ) / RNDI(JJ)
                IF(JJ == 1 .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
                    KDTy2 = -WALLFLUX(IC, iBotWall, KC) * COE222
                ELSE IF(JJ == NCL2 .AND. iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux) THEN
                    KDTy1 = -WALLFLUX(IC, iTopWall, KC) * COE221
                ELSE
                END IF
                DY_KDT = (KDTy1- KDTy2) * COE22  !DYFI(JJ) * RCCI1(JJ) * CTHECD

                ! $ \Delta_z [k^z (\Delta_z T)] $
                DZ_KDT = ( THERMCONDT(IC, JC, KP)  + THERMCONDT(IC, JC, KC) ) * &
                         ( TEMPERATURE(IC, JC, KP) - TEMPERATURE(IC, JC, KC) ) - &
                         ( THERMCONDT(IC, JC, KM)  + THERMCONDT(IC, JC, KC) ) * &
                         ( TEMPERATURE(IC, JC, KC) - TEMPERATURE(IC, JC, KM) )
                DZ_KDT = DZ_KDT * COE33  !ZND2CL * DZI * DZI * RCCI2(JJ) * CTHECD

                !=======EXPLICIT RHS OF ENERGY EQUATION ====================
                RHS_ENERGY(IC, JC, KC) = - DX_DHU - DY_DHV - DZ_DHW + DX_KDT + DY_KDT + DZ_KDT

            END DO
        END DO
    END DO

    !CALL DEBUG_WRT_LOCAL(RHS_ENERGY, 1, N2DO(MYID)) !test

    RETURN
END SUBROUTINE
