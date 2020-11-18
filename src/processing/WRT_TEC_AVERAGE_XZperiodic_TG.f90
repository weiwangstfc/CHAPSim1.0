!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> module: VARS_AVERAGED_TG
!> subroutine: MEMO_ALLOCT_AVERAGE_TG
!> subroutine: MEMO_DEALLT_AVERAGE_TG
!> subroutine: WRT_AVERAGE_PPED_TG
!> subroutine: WRT_AVERAGE_PPED_TG_WRT_TEC
!> subroutine: WRT_AVERAGE_PPED_TG_CALC_RSTE_BUDG
!> subroutine: Cf_Utau_TG
!> subroutine: WRT_AVERAGE_PPED_TG_GATHER
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE VARS_AVERAGED_TG
    USE WRT_INFO

    REAL(WP) :: Cf_LW_TG, Cf_UW_TG, Cf_ave_TG
    REAL(WP) :: U_tau_LW_TG, Re_tau_LW_TG, U_tau_ave_TG
    REAL(WP) :: U_tau_UW_TG, Re_tau_UW_TG, Re_tau_ave_TG

    CHARACTER(15) :: PNTIM

    !===============GLOABL DATA===================
    REAL(WP), ALLOCATABLE :: U1xztL_F0_tg( :, : )
    REAL(WP), ALLOCATABLE :: UPxztL_F0_tg( :, : )
    REAL(WP), ALLOCATABLE :: U2xztL_F0_tg( :, : )
    REAL(WP), ALLOCATABLE :: U3xztL_F0_tg( :, : )

    REAL(WP), ALLOCATABLE :: DVDL1xztL_F0_tg( :, :, : )
    REAL(WP), ALLOCATABLE :: DVDLPxztL_F0_tg( :, :, : )
    REAL(WP), ALLOCATABLE :: DVDL2xztL_F0_tg( :, :, : )

    REAL(WP), ALLOCATABLE :: U2PER( :, :, : )
    REAL(WP), ALLOCATABLE :: U3PER( :, :, :, :)
    REAL(WP), ALLOCATABLE :: VORper2( :, :)
    REAL(WP), ALLOCATABLE :: DUDX1(:, :, :)
    REAL(WP), ALLOCATABLE :: Skewness(:, :)

    REAL(WP), ALLOCATABLE :: Budg_productn( :, :)
    REAL(WP), ALLOCATABLE :: Budg_DISsipat( :, :)
    REAL(WP), ALLOCATABLE :: Budg_Tpr_diff( :, :)
    REAL(WP), ALLOCATABLE :: Budg_VIS_diff( :, :)
    REAL(WP), ALLOCATABLE :: Budg_VPG_diff( :, :)
    REAL(WP), ALLOCATABLE :: Budg_VPG_stra( :, :)


END MODULE

!**********************************************************************************************************************************
SUBROUTINE MEMO_ALLOCT_AVERAGE_TG
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    ALLOCATE( U1xztL_F0_tg( NCL2, NDV + 1 ) ); U1xztL_F0_tg = 0.0_WP
    ALLOCATE( UPxztL_F0_tg( NCL2, NDV   )  ); UPxztL_F0_tg = 0.0_WP
    ALLOCATE( U2xztL_F0_tg( NCL2, NDV * (7 - NDV) / 2 + NDV - 3 )   ); U2xztL_F0_tg = 0.0_WP
    ALLOCATE( U3xztL_F0_tg( NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8  ) ); U3xztL_F0_tg = 0.0_WP

    ALLOCATE( DVDL1xztL_F0_tg( NCL2, NDV, NDV  ) ); DVDL1xztL_F0_tg = 0.0_WP
    ALLOCATE( DVDLPxztL_F0_tg( NCL2, NDV, NDV  ) ); DVDLPxztL_F0_tg = 0.0_WP
    ALLOCATE( DVDL2xztL_F0_tg( NCL2, NDV * NDV, NDV * NDV  ) ); DVDL2xztL_F0_tg = 0.0_WP

    ALLOCATE( U2PER(NCL2, NDV, NDV)                     ); U2PER = 0.0_WP
    ALLOCATE( U3PER(NCL2, NDV, NDV, NDV)            ); U3PER = 0.0_WP
    ALLOCATE( VORper2(NCL2, NDV)                     ); VORper2 = 0.0_WP
    ALLOCATE( DUDX1(NCL2, NDV, NDV)                     ); DUDX1 = 0.0_WP
    ALLOCATE( Skewness(NCL2, NDV)                    ); Skewness = 0.0_WP


    ALLOCATE(Budg_productn( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_productn = 0.0_WP
    ALLOCATE(Budg_Tpr_diff( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_Tpr_diff = 0.0_WP
    ALLOCATE(Budg_VIS_diff( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_VIS_diff = 0.0_WP
    ALLOCATE(Budg_VPG_diff( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_VPG_diff = 0.0_WP
    ALLOCATE(Budg_VPG_stra( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_VPG_stra = 0.0_WP
    ALLOCATE(Budg_DISsipat( NCL2, (NDV * (7 - NDV)) / 2 + NDV - 3 ) ); Budg_DISsipat = 0.0_WP


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE MEMO_DEALLT_AVERAGE_TG
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    IMPLICIT NONE

    DEALLOCATE(U1xztL_F0_tg)
    DEALLOCATE(UPxztL_F0_tg)
    DEALLOCATE(U2xztL_F0_tg)
    DEALLOCATE(U3xztL_F0_tg)

    DEALLOCATE(DVDL1xztL_F0_tg)
    DEALLOCATE(DVDLPxztL_F0_tg)
    DEALLOCATE(DVDL2xztL_F0_tg)

    DEALLOCATE( U2PER     )
    DEALLOCATE( U3PER     )
    DEALLOCATE( VORper2    )
    DEALLOCATE( DUDX1     )
    DEALLOCATE( Skewness  )

    DEALLOCATE(Budg_productn )
    DEALLOCATE(Budg_Tpr_diff )
    DEALLOCATE(Budg_VIS_diff )
    DEALLOCATE(Budg_VPG_diff )
    DEALLOCATE(Budg_VPG_stra )
    DEALLOCATE(Budg_DISsipat )

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_TG
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    CALL WRT_AVERAGE_PPED_TG_GATHER
    IF(MYID == 0) THEN

        CALL Cf_Utau_TG
        CALL WRT_AVERAGE_PPED_TG_CALC_RSTE_BUDG

        CALL WRT_AVERAGE_PPED_TG_WRT_TEC

        CALL MEMO_DEALLT_AVERAGE_TG
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_TG_WRT_TEC
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: TECFLG1, TECFLG2
    REAL(WP) :: COE1, COE2, COE
    REAL(WP) :: urms, vrms, WRms, uv, uw, vw, yplus

    INTEGER(4) :: INN
    INTEGER(4) :: J, JJ
    INTEGER(4) :: L, IP, M
    INTEGER(4) :: N2DOID



    IF(MYID /= 0) RETURN

    !================WRITE DATA OUT =======================

    TECFLG1 = 200

    COE1 = DSQRT(DABS(Cf_LW_tg * 0.5_WP))
    COE2 = DSQRT(DABS(Cf_UW_tg * 0.5_WP))

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_tg

    OPEN(TECFLG1, FILE = TRIM(FilePath4) // 'Result.TG.Reynolds.Averaged.Flow.' // TRIM(PNTIM) // '.Profile.plt')
    WRITE(TECFLG1, '(A)') 'TITLE = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    WRITE(TECFLG1, '(A)') &
    'variables = "Y", "Y+", "Ut", "Ux", "Uy", "Uz", "P", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
    '"vorxper2", "voryper2", "vorzper2", "Suuu", "Svvv", "Swww", ', &
    '"uu_Prodc", "uu_TdIFf", "uu_VdIFf", "uu_VPGdIFf", "uu_VPG-stAIn", "uu_DISsip", ', &
    '"vv_Prodc", "vv_TdIFf", "vv_VdIFf", "vv_VPGdIFf", "vv_VPG-stAIn", "vv_DISsip", ', &
    '"ww_Prodc", "ww_TdIFf", "ww_VdIFf", "ww_VPGdIFf", "ww_VPG-stAIn", "ww_DISsip", ', &
    '"uv_Prodc", "uv_TdIFf", "uv_VdIFf", "uv_VPGdIFf", "uv_VPG-stAIn", "uv_DISsip", ', &
    '"uw_Prodc", "uw_TdIFf", "uw_VdIFf", "uw_VPGdIFf", "uw_VPG-stAIn", "uw_DISsip", ', &
    '"vw_Prodc", "vw_TdIFf", "vw_VdIFf", "vw_VPGdIFf", "vw_VPG-stAIn", "vw_DISsip"  ', &
    '"dmeanUdy", "dmeanWdy"'
    WRITE(TECFLG1, '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '


    DO J = 1, NCL2
        IF(YND(J) < 0.0_WP) THEN
            COE = COE1
        ELSE
            COE = COE2
        END IF
        urms = DSQRT( DABS(U2PER(J, 1, 1) ) )
        vrms = DSQRT( DABS(U2PER(J, 2, 2) ) )
        WRms = DSQRT( DABS(U2PER(J, 3, 3) ) )
        uv = U2PER(J, 1, 2)
        uw = U2PER(J, 1, 3)
        vw = U2PER(J, 2, 3)
        ypluS = (1.0_WP - Abs(YCC(J))) * REN * COE


        WRITE(TECFLG1, '(57ES22.14)') &
        YCC(J), yplus, COE, &
        U1xztL_F0_tg(J, 1:4), &
        urms, vrms, WRms, uv, uw, vw, &
        VORper2(J, 1), VORper2(J, 2), VORper2(J, 3), &
        Skewness(J, 1:3), &
        Budg_productn(J, 1),Budg_Tpr_diff(J, 1),Budg_VIS_diff(J, 1), &
        Budg_VPG_diff(J, 1),Budg_VPG_stra(J, 1),Budg_DISsipat(J, 1), &
        Budg_productn(J, 4),Budg_Tpr_diff(J, 4),Budg_VIS_diff(J, 4), &
        Budg_VPG_diff(J, 4),Budg_VPG_stra(J, 4),Budg_DISsipat(J, 4), &
        Budg_productn(J,6),Budg_Tpr_diff(J,6),Budg_VIS_diff(J,6), &
        Budg_VPG_diff(J,6),Budg_VPG_stra(J,6),Budg_DISsipat(J,6), &
        Budg_productn(J, 2),Budg_Tpr_diff(J, 2),Budg_VIS_diff(J, 2), &
        Budg_VPG_diff(J, 2),Budg_VPG_stra(J, 2),Budg_DISsipat(J, 2), &
        Budg_productn(J, 3),Budg_Tpr_diff(J, 3),Budg_VIS_diff(J, 3), &
        Budg_VPG_diff(J, 3),Budg_VPG_stra(J, 3),Budg_DISsipat(J, 3), &
        Budg_productn(J, 5),Budg_Tpr_diff(J, 5),Budg_VIS_diff(J, 5), &
        Budg_VPG_diff(J, 5),Budg_VPG_stra(J, 5),Budg_DISsipat(J, 5), &
        DUDX1(J, 1, 2), DUDX1(J, 3, 2)

        !WRITE(*, *) YCC(J), VORper2(J, 1) / COE1 / COE1, VORper2(J, 2) / COE1 / COE1, VORper2(J, 3) / COE1 / COE1


    END DO
    CLOSE(TECFLG1)

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_TG_CALC_RSTE_BUDG
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, H, LMN, LMH, LNH, LMNH, L
    INTEGER(4) :: H1(3), P1(3), L1(3), L2(3)


    DO J = 1, NCL2

        !====u'u', u'v', u'w', v'w', w'w'============
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                LMN = (M * (7-M)) / 2 + N - 3
                U2PER(J, M, N) = U2xztL_F0_tg(J, LMN) - U1xztL_F0_tg(J, M) * U1xztL_F0_tg(J, N)
            END DO
        END DO
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) THEN
                    U2PER(J, M, N) = U2PER(J, N, M)
                END IF
            END DO
        END DO

        !== Above tested OK ======

        !=====u'u'u', u'u'v', u'u'w', u'v'v', u'v'w', u'w'w', v'v'v', v'v'w', v'w'w', w'w'w ===
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                DO H = 1, NDV
                    IF(N >  H) CYCLE
                    LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
                    LMN = (M * (7-M)) / 2 + N - 3
                    LMH = (M * (7-M)) / 2 + H- 3
                    LNH = (N * (7-N)) / 2 + H- 3

                    U3PER(J, M, N,H) = U3xztL_F0_tg(J, LMNH)  &
                    - U1xztL_F0_tg(J, M) * U2xztL_F0_tg(J, LNH) &
                    - U1xztL_F0_tg(J, N) * U2xztL_F0_tg(J, LMH) &
                    - U1xztL_F0_tg(J,H) * U2xztL_F0_tg(J, LMN)  &
                    + 2.0_WP * U1xztL_F0_tg(J, M) * U1xztL_F0_tg(J, N) * U1xztL_F0_tg(J,H)
                END DO
            END DO
        END DO

        !{111} {112} {113}
        ![121] {122} {123}
        ![131] [132] {133}
        ![211] [212] [213]
        ![221] {222} {223}
        ![231] [232] {233}
        ![311] [312] [313]
        ![321] [322] [323]
        ![331] [332] {333}


        U3PER(J, 1, 2, 1) = U3PER(J, 1, 1, 2)
        U3PER(J, 1, 3, 1) = U3PER(J, 1, 1, 3)
        U3PER(J, 1, 3, 2) = U3PER(J, 1, 2, 3)

        U3PER(J, 2, 1, 1) = U3PER(J, 1, 1, 2)
        U3PER(J, 2, 1, 2) = U3PER(J, 1, 2, 2)
        U3PER(J, 2, 1, 3) = U3PER(J, 1, 2, 3)

        U3PER(J, 2, 2, 1) = U3PER(J, 1, 2, 2)
        U3PER(J, 2, 3, 1) = U3PER(J, 1, 2, 3)
        U3PER(J, 2, 3, 2) = U3PER(J, 2, 2, 3)

        U3PER(J, 3, 1, 1) = U3PER(J, 1, 1, 3)
        U3PER(J, 3, 1, 2) = U3PER(J, 1, 2, 3)
        U3PER(J, 3, 1, 3) = U3PER(J, 1, 3, 3)

        U3PER(J, 3, 2, 1) = U3PER(J, 1, 2, 3)
        U3PER(J, 3, 2, 2) = U3PER(J, 2, 2, 3)
        U3PER(J, 3, 2, 3) = U3PER(J, 2, 3, 3)

        U3PER(J, 3, 3, 1) = U3PER(J, 1, 3, 3)
        U3PER(J, 3, 3, 2) = U3PER(J, 2, 3, 3)

        !==========
        DO M = 1, NDV
            DO N = 1, NDV
                IF(N == 2) THEN
                    IF(J == 1) THEN
                        DUDX1(J, M, N) = ( ( YCL2ND_WFB(J + 1) * U1xztL_F0_tg(J, M) + &
                        YCL2ND_WFF(J + 1) * U1xztL_F0_tg(J + 1, M) ) - &
                        0.0_WP  ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        DUDX1(J, M, N) = (   0.0_WP - &
                        ( YCL2ND_WFF(J) * U1xztL_F0_tg(J, M) + YCL2ND_WFB(J) * U1xztL_F0_tg(J - 1, M) )  ) * DYFI(J)
                    ELSE
                        DUDX1(J, M, N) = ( ( YCL2ND_WFB(J + 1) * U1xztL_F0_tg(J, M) + &
                        YCL2ND_WFF(J + 1) * U1xztL_F0_tg(J + 1, M) ) - &
                        ( YCL2ND_WFF(J) * U1xztL_F0_tg(J, M) + YCL2ND_WFB(J) * U1xztL_F0_tg(J - 1, M) ) ) * DYFI(J)
                    END IF
                ELSE
                    DUDX1(J, M, N) = 0.0_WP
                END IF

            END DO
        END DO

    END DO


    DO J = 1, NCL2
        DO M = 1, NDV
            Skewness(J, M) = U3PER(J, M, M, M) / ( U2PER(J, M, M)**(3.0_WP / 2.0_WP) )
        END DO
    END DO


    !=== Refer to note on 19/11 /2014
    DO J = 1, NCL2
        VORper2(J, 1) = DVDL2xztL_F0_tg(J,8,8) + DVDL2xztL_F0_tg(J,6,6) - 2.0_WP * DVDL2xztL_F0_tg(J,8,6)  &
        + DUDX1(J, 3, 2)**2 + DUDX1(J, 2, 3)**2  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 3, 2) * DUDX1(J, 3, 2)  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 2, 3) * DUDX1(J, 2, 3)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 3, 2) * DUDX1(J, 2, 3)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 2, 3) * DUDX1(J, 3, 2)  &
        -2.0_WP * DUDX1(J, 3, 2) * DUDX1(J, 2, 3)
        VORper2(J, 2) = DVDL2xztL_F0_tg(J, 3, 3) + DVDL2xztL_F0_tg(J,7,7) - 2.0_WP * DVDL2xztL_F0_tg(J, 3,7)  &
        + DUDX1(J, 1, 3)**2 + DUDX1(J, 3, 1)**2  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 1, 3) * DUDX1(J, 1, 3)  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 3, 1) * DUDX1(J, 3, 1)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 1, 3) * DUDX1(J, 3, 1)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 3, 1) * DUDX1(J, 1, 3)  &
        -2.0_WP * DUDX1(J, 1, 3) * DUDX1(J, 3, 1)
        VORper2(J, 3) = DVDL2xztL_F0_tg(J, 4, 4) + DVDL2xztL_F0_tg(J, 2, 2) - 2.0_WP * DVDL2xztL_F0_tg(J, 4, 2)  &
        + DUDX1(J, 2, 1)**2 + DUDX1(J, 1, 2)**2  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 2, 1) * DUDX1(J, 2, 1)  &
        -2.0_WP * DVDL1xztL_F0_tg(J, 1, 2) * DUDX1(J, 1, 2)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 2, 1) * DUDX1(J, 1, 2)  &
        + 2.0_WP * DVDL1xztL_F0_tg(J, 1, 2) * DUDX1(J, 2, 1)  &
        -2.0_WP * DUDX1(J, 2, 1) * DUDX1(J, 1, 2)
        VORper2(J, 1) = DSQRT(DABS(VORper2(J, 1)))
        VORper2(J, 2) = DSQRT(DABS(VORper2(J, 2)))
        VORper2(J, 3) = DSQRT(DABS(VORper2(J, 3)))

    END DO

    !============================================================================================================
    DO J = 1, NCL2
        !============================BUDGET TERMS =========================================================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                !================PRODUCTUON TERMS =======================================
                Budg_productn(J, L) = U2PER(J, 1, M) * DUDX1(J, N, 1) + &
                U2PER(J, 2, M) * DUDX1(J, N, 2) + &
                U2PER(J, 3, M) * DUDX1(J, N, 3) + &
                U2PER(J, 1, N) * DUDX1(J, M, 1) + &
                U2PER(J, 2, N) * DUDX1(J, M, 2) + &
                U2PER(J, 3, N) * DUDX1(J, M, 3)
                Budg_productn(J, L) = Budg_productn(J, L) * (-1.0_WP)

                !================ DIFFUSION TERMS ========================================

                !================= TURBULENCE TRANSPORT RATE == (T) ==========
                H = 2
                IF(J == 1) THEN
                    Budg_Tpr_diff(J, L) = ( ( YCL2ND_WFF(J + 1) * U3PER(J + 1, M, N,H) + &
                    YCL2ND_WFB(J + 1) * U3PER(J,  M, N,H) ) - &
                    0.0_WP ) * DYFI(J)
                ELSE IF (J == NCL2) THEN
                    Budg_Tpr_diff(J, L) = ( 0.0_WP - &
                    ( YCL2ND_WFF(J) * U3PER(J,  M, N,H) + &
                    YCL2ND_WFB(J) * U3PER(J - 1, M, N,H) ) ) * DYFI(J)
                ELSE
                    Budg_Tpr_diff(J, L) = ( ( YCL2ND_WFF(J + 1) * U3PER(J + 1, M, N,H) + &
                    YCL2ND_WFB(J + 1) * U3PER(J,  M, N,H) ) - &
                    ( YCL2ND_WFF(J) * U3PER(J,  M, N,H) + &
                    YCL2ND_WFB(J) * U3PER(J - 1, M, N,H) ) ) * DYFI(J)
                END IF
                Budg_Tpr_diff(J, L) = -1.0_WP *Budg_Tpr_diff(J, L)


                !=============== VISCOUS DIFFUSION TERM == (D) ===============

                IF(J == 1) THEN
                    Budg_VIS_diff(J, L) = U2PER(J + 1, M, N) * APVR(J, 1) + U2PER(J, M, N) * ACVR(J, 1) + 0.0_WP * AMVR(J, 1)
                ELSE IF(J == NCL2) THEN
                    Budg_VIS_diff(J, L) = 0.0_WP * APVR(J, 1) + U2PER(J, M, N) * ACVR(J, 1) + U2PER(J - 1, M, N) * AMVR(J, 1)
                ELSE
                    Budg_VIS_diff(J, L) = U2PER(J + 1, M, N) * APVR(J, 1) + &
                    U2PER(J,  M, N) * ACVR(J, 1) + U2PER(J - 1, M, N) * AMVR(J, 1)
                END IF
                Budg_VIS_diff(J, L) = Budg_VIS_diff(J, L) * CVISC

                !============= VELOCITY- PRessure GRADIENT (DIFFUSION) ======
                IF(M == 2 .AND. N == 2) THEN
                    IF(J == 1) THEN
                        Budg_VPG_diff(J, L) = &
                        ((  YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, M) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, M) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * 0.0_WP   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)  &
                        + (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, N) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, N) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * 0.0_WP   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        Budg_VPG_diff(J, L) = &
                        ((  YCL2ND_WFF(J + 1) * 0.0_WP   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, M) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, M) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)  &
                        + (( YCL2ND_WFF(J + 1) * 0.0_WP   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, N) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, N) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)
                    ELSE
                        Budg_VPG_diff(J, L) =  &
                        (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, M) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, M) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, M) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, M) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)  &
                        + (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, N) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, N) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, N) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, N) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)
                    END IF
                END IF

                IF((N /= 2) .AND. (m == 2)) THEN
                    IF(J == 1) THEN
                        Budg_VPG_diff(J, L) = &
                        (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, N) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, N) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * 0.0_WP   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        Budg_VPG_diff(J, L) = &
                        (( YCL2ND_WFF(J + 1) * 0.0_WP   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, N) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, N) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)
                    ELSE
                        Budg_VPG_diff(J, L) =  &
                        (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, N) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, N) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, N) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, N) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J, N) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J, N) ) ) &
                        ) * DYFI(J)

                    END IF
                END IF

                IF((M /= 2) .AND. (N == 2)) THEN
                    IF(J == 1) THEN
                        Budg_VPG_diff(J, L) = &
                        ((  YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, M) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, M) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * 0.0_WP   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        Budg_VPG_diff(J, L) = &
                        ((  YCL2ND_WFF(J + 1) * 0.0_WP   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, M) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, M) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)
                    ELSE
                        Budg_VPG_diff(J, L) =  &
                        (( YCL2ND_WFF(J + 1) * ( UPxztL_F0_tg(J + 1, M) - U1xztL_F0_tg(J + 1, 4) * U1xztL_F0_tg(J + 1, M) )   &
                        + YCL2ND_WFB(J + 1) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        -( YCL2ND_WFF(J) * ( UPxztL_F0_tg(J - 1, M) - U1xztL_F0_tg(J - 1, 4) * U1xztL_F0_tg(J - 1, M) )   &
                        + YCL2ND_WFB(J) * ( UPxztL_F0_tg(J,  M) - U1xztL_F0_tg(J,  4) * U1xztL_F0_tg(J,  M) ) ) &
                        ) * DYFI(J)
                    END IF
                END IF

                IF (.not.(M == 2 .OR. N == 2)) THEN
                    Budg_VPG_diff(J, L) = 0.0_WP
                END IF

                Budg_VPG_diff(J, L) = -1.0_WP *Budg_VPG_diff(J, L)


                !================ VELOCITY- PRessure GRADIENT (STRAIN) ==============================
                Budg_VPG_stra(J, L) = DVDLPxztL_F0_tg(J, M, N) - U1xztL_F0_tg(J, 4) * DVDL1xztL_F0_tg(J, M, N) + &
                DVDLPxztL_F0_tg(J, N, M) - U1xztL_F0_tg(J, 4) * DVDL1xztL_F0_tg(J, N, M)


                !================ DISSIPATION RATE =================================================
                H1(1) = 1
                H1(2) = 2
                H1(3) = 3
                P1(1) = 1
                P1(2) = 2
                P1(3) = 3
                L1(1:3) = (M - 1) * 3 + H1(1:3)
                L2(1:3) = (N - 1) * 3 + P1(1:3)
                Budg_DISsipat(J, L) =  DVDL2xztL_F0_tg(J, L1(1), L2(1)) &
                + DVDL2xztL_F0_tg(J, L1(2), L2(2)) &
                + DVDL2xztL_F0_tg(J, L1(3), L2(3)) &
                - DVDL1xztL_F0_tg(J, N, 1) * DUDX1(J, M, 1) &
                - DVDL1xztL_F0_tg(J, N, 2) * DUDX1(J, M, 2) &
                - DVDL1xztL_F0_tg(J, N, 3) * DUDX1(J, M, 3) &
                - DVDL1xztL_F0_tg(J, M, 1) * DUDX1(J, N, 1) &
                - DVDL1xztL_F0_tg(J, M, 2) * DUDX1(J, N, 2) &
                - DVDL1xztL_F0_tg(J, M, 3) * DUDX1(J, N, 3) &
                + DUDX1(J, M, 1) * DUDX1(J, N, 1) &
                + DUDX1(J, M, 2) * DUDX1(J, N, 2) &
                + DUDX1(J, M, 3) * DUDX1(J, N, 3)
                Budg_DISsipat(J, L) = Budg_DISsipat(J, L) * CVISC * 2.0_WP
            END DO
        END DO

    END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE Cf_Utau_TG
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: TECFLG = 200
    INTEGER(4) :: J
    REAL(WP) :: MDOt, ARea_tg

    ! method 1 and 2 give slightly dIFfeRENt results.
    ! method 2
    !Cf_LW_TG = 2.0_WP * (DVDL1xztL_F0_tg(1,   1, 2)) / REN
    !Cf_UW_TG = -2.0_WP * (DVDL1xztL_F0_tg(NCL2, 1, 2)) / REN

    ! method 1
    Cf_UW_TG = -2.0_WP * (U1xztL_F0_tg(NCL2, 1) - 0.0_WP) / (YCC(NCL2) - YND(NND2)) / REN
    IF(iCase == iPIPEC) THEN
        Cf_LW_TG = Cf_UW_TG
    ELSE
        Cf_LW_TG = 2.0_WP * (U1xztL_F0_tg(1, 1) - 0.0_WP) / (YCC(1) - YND(1)) / REN
    END IF

    U_tau_LW_TG = DSQRT(DABS(Cf_LW_TG) * 0.5_WP)
    U_tau_UW_TG = DSQRT(DABS(Cf_UW_TG) * 0.5_WP)

    Re_tau_LW_tg = REN * U_tau_LW_TG
    Re_tau_UW_tg = REN * U_tau_UW_TG


    Cf_ave_tg = 0.5* (DABS(Cf_LW_TG) + DABS(Cf_UW_TG))
    U_tau_ave_TG = DSQRT(DABS(Cf_ave_TG) * 0.5_WP)
    Re_tau_ave_tg = REN * U_tau_ave_TG


    MDOt = 0.0_WP
    ARea_TG = 0.0_WP
    DO J = 1, NCL2
        MDOt = MDOt + U1xztL_F0_tg(J, 1) / DYFI(J) / RCCI1(J)
        ARea_TG = ARea_TG + 1.0_WP / DYFI(J) / RCCI1(J)
    END DO
    MDOt = MDOt/ARea_TG



    OPEN(TECFLG, FILE = TRIM(FilePath4) // 'Result.TG.Reynolds.Averaged.Cf.Utau.table.plt')

    IF(iCase /= IPIPEC)  &
    WRITE(TECFLG, '(A, 1ES15.8)') 'Cf on the lower walL = ', Cf_LW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'Cf on the upper walL = ', Cf_UW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'Cf average of walls = ', Cf_ave_TG

    WRITE(TECFLG, '(A)') '   '

    IF(iCase /= IPIPEC)  &
    WRITE(TECFLG, '(A, 1ES15.8)') 'U_tau on the lower walL = ', U_tau_LW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'U_tau on the upper walL = ', U_tau_UW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'U_tau average of walls = ', U_tau_ave_TG

    WRITE(TECFLG, '(A)') '   '

    IF(iCase /= IPIPEC) &
    WRITE(TECFLG, '(A, 1ES15.8)') 'Re_tau on the lower walL = ', Re_tau_LW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'Re_tau on the upper walL = ', Re_tau_UW_TG
    WRITE(TECFLG, '(A, 1ES15.8)') 'Re_tau average of walls = ', Re_tau_ave_TG

    WRITE(TECFLG, '(A)') '   '
    WRITE(TECFLG, '(A, 1ES15.8)') 'Bulk Mass Flux = ', MDOt

    CLOSE(TECFLG)



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_TG_GATHER
    USE VARS_AVERAGED_TG
    USE mesh_info
    USE init_info
    USE flow_info
    USE thermal_info
    USE postprocess_info
    IMPLICIT NONE

    INTEGER(4) :: TECFLG1, TECFLG2
    REAL(WP) :: COE1, COE2, COE
    !REAL(WP) :: urms, vrms, WRms, uv, uw, vw

    INTEGER(4) :: INN
    INTEGER(4) :: J, JJ
    INTEGER(4) :: L, IP, M, N, H, P, L1, L2
    INTEGER(4) :: N2DOID
    REAL(WP) :: D1AUX (N2DO(MYID), NDV + 1,                               1 : NPTOT)
    REAL(WP) :: D2AUX (N2DO(MYID), NDV,                                 1 : NPTOT)
    REAL(WP) :: D3AUX (N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3),                1 : NPTOT)
    REAL(WP) :: D4AUX (N2DO(MYID), (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8),  1 : NPTOT)

    REAL(WP) :: T1AUX (N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: T2AUX (N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: T3AUX (N2DO(MYID), NDV * NDV, NDV * NDV, 1 : NPTOT)



    INN = N2DO(MYID) * (NDV + 1)
    CALL MPI_GATHER( U1xztL_tg, INN, MPI_DOUBLE_PRECISION, D1AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * NDV
    CALL MPI_GATHER( UPxztL_tg, INN, MPI_DOUBLE_PRECISION, D2AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
    CALL MPI_GATHER( U2xztL_tg, INN, MPI_DOUBLE_PRECISION, D3AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)


    INN = N2DO(MYID) * (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
    CALL MPI_GATHER( U3xztL_tg, INN, MPI_DOUBLE_PRECISION, D4AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)


    INN = N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDL1xztL_tg, INN, MPI_DOUBLE_PRECISION, T1AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDLPxztL_tg, INN, MPI_DOUBLE_PRECISION, T2AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * (NDV * NDV) * NDV * NDV
    CALL MPI_GATHER( DVDL2xztL_tg, INN, MPI_DOUBLE_PRECISION, T3AUX, INN,  &
    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)


    IF(MYID == 0) THEN
        !==================================================================
        CALL MEMO_ALLOCT_AVERAGE_TG
        DO IP = 0, NPSLV
            N2DOID=JDEWT(IP) - JDSWT(IP) + 1
            DO J = 1, N2DOID
                JJ = JDSWT(IP) - 1 + J

                DO L = 1, NDV + 1
                    U1xztL_F0_tg(JJ, L) = D1AUX(J, L, IP + 1)
                END DO

                DO L = 1, NDV
                    UPxztL_F0_tg(JJ, L) = D2AUX(J, L, IP + 1)
                END DO

                DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                    U2xztL_F0_tg(JJ, L) = D3AUX(J, L, IP + 1)
                END DO

                DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
                    U3xztL_F0_tg(JJ, L) = D4AUX(J, L, IP + 1)
                END DO

                DO L = 1, NDV
                    DO M = 1, NDV
                        DVDL1xztL_F0_tg(JJ, L, M) = T1AUX(J, L, M, IP + 1)
                    END DO
                END DO

                DO L = 1, NDV
                    DO M = 1, NDV
                        DVDLPxztL_F0_tg(JJ, L, M) = T2AUX(J, L, M, IP + 1)
                    END DO
                END DO

                DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                    DO M = 1, NDV
                        DVDL2xztL_F0_tg(JJ, L, M) = T3AUX(J, L, M, IP + 1)
                    END DO
                END DO

                DO M = 1, NDV
                    DO N = 1, NDV
                        DO H = 1, NDV
                            DO P = 1, NDV
                                L1 = (M - 1) * 3 + H
                                L2 = (N - 1) * 3 + P
                                DVDL2xztL_F0_tg(JJ, L1, L2) = T3AUX(J, L1, L2, IP + 1)
                            END DO
                        END DO
                    END DO
                END DO

            END DO
        END DO

    END IF

    RETURN
END SUBROUTINE
