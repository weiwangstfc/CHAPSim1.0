!**********************************************************************************************************************************
!> @brief
!>        to calculate mass flux
!> @details
!> SUBROUTINE: MASSFLUX_CALC_io(in MYID = all)
!> SUBROUTINE: VELOCITY_CALC_io(in MYID = all)
!> SUBROUTINE: MU_Staggered (in MYID = all)
!> SUBROUTINE: DENSITY_Staggered (in MYID = all)
!> SUBROUTINE: DENSITY_IMPLICIT_add (in MYID = all)
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE MASSFLUX_CALC_io
    USE FLOW_INFO
    USE THERMAL_INFO
    USE MESH_INFO
    USE init_info
    IMPLICIT NONE
    INTEGER(4) :: NXI, NYI, IDR
    INTEGER(4) :: IC, JC, KC, KS

    !===================================== X =============================
    IDR = 1
    NXI = 1
    IF(TgFlowFlg) NXI = 2
    NYI = 1
    DO IC = NXI, NCL1_io
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) + QTMP_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !===================================== Z =============================
    IDR = 3
    NXI = 1
    NYI = 1
    DO IC = NXI, NCL1_io
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) + RHSLLPHI_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !=====================================Y =============================
    IDR = 2
    NXI = 1
    NYI = 1
    IF(MYID == 0 .AND. iCase /= IBox3P) NYI = 2
    DO IC = NXI, NCL1_io
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                G_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) + DPH_io(IC, JC, KC)
            END DO
        END DO
    END DO

    !==========below pARt IS repeating interfaceS ===========
    !        IF(MYID == 0) THENMASSFLUX_CALC_io
    !            IF(iCase == IPIPEC) THEN
    !                DO KC = 1, NCL3
    !                    KS = KSYM(KC)
    !                    DO IC = NXI, NCL1_io
    !                        G_io(IC, 1, KC, IDR) = 0.0_WP
    !                        G_io(IC, 0, KC, IDR) = G_io(IC, 2, KS, IDR)
    !                    END DO
    !                END DO
    !            ELSE
    !                DO KC = 1, NCL3
    !                    DO IC = NXI, NCL1_io
    !                        G_io(IC, 1, KC, IDR) = 0.0_WP
    !                        G_io(IC, 0, KC, IDR) = 0.0_WP
    !                    END DOMASSFLUX_CALC_io
    !                END DO
    !            END IF
    !        END IF

    !        IF(MYID == NPSLV) THEN
    !            DO KC = 1, NCL3
    !                DO IC = NXI, NCL1_io
    !                    G_io(IC, N2DO(MYID) + 1, KC, IDR) = 0.0_WP
    !                END DO
    !            END DO
    !        END IF

    !CALL INTFC_G_WALL

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE VELOCITY_CALC_io
    USE MESH_INFO
    USE THERMAL_INFO
    USE FLOW_INFO
    USE init_info
    IMPLICIT NONE
    INTEGER(4) :: IDR, NYI, IC, JC, KC, IM, JJ, KM, JM!, KS
    !REAL(WP) :: RHO

    !===================================== X =============================
    IDR = 1
    NYI = 1
    DO IC = 1, NCL1_io
        IM = IMV_io(IC)
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                !RHO = ( DENSITY(IC, JC, KC) + DENSITY(IM, JC, KC) ) * XND2CL
                Q_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) / D_STG(IC, JC, KC, IDR)
            END DO
        END DO
    END DO

    !===================================== Z =============================
    IDR = 3
    NYI = 1
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = NYI, N2DO(MYID)
            DO IC = 1, NCL1_io
                !RHO = ( DENSITY(IC, JC, KC) + DENSITY(IC, JC, KM) ) * ZND2CL
                Q_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) / D_STG(IC, JC, KC, IDR)
            END DO
        END DO
    END DO


    !=====================================Y =============================
    IDR = 2
    NYI = 1
    IF(MYID == 0 .AND. iCase /= IBox3P) NYI = 2
    DO JC = NYI, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_io
                !RHO = YCL2ND_WFF(JJ) * DENSITY(IC, JC, KC) + &
                !      YCL2ND_WFB(JJ) * DENSITY(IC, JM, KC)
                Q_io(IC, JC, KC, IDR) = G_io(IC, JC, KC, IDR) / D_STG(IC, JC, KC, IDR)
            END DO
        END DO
    END DO

    CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, Q_io)
    CALL BC_WALL_Q_io

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE MU_Staggered
    USE FLOW_INFO
    USE THERMAL_INFO
    USE MESH_INFO
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IDR, NYI, IC, JC, KC, IM, JM, KM, JJ


    IF(iThermoDynamics /= 1) THEN
        MU_STG = M_inlet
        RETURN
    END IF

    !===================================== X - Y = (i', J', K) ============================
    IDR = 1
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO IC = 1, NCL1_io
        IM = IMV_io(IC)
        DO JC = NYI, N2DO(MYID)
            JJ = JCL2G(JC)
            JM = JLMV(JC)
            DO KC = 1, NCL3
                MU_STG(IC, JC, KC, IDR) = &
                ( ( Viscousity0(IC, JC, KC)  + Viscousity0(IM, JC, KC) ) * YCL2ND_WFF(JJ) + &
                ( Viscousity0(IC, JM, KC)  + Viscousity0(IM, JM, KC) ) * YCL2ND_WFB(JJ) ) * XND2CL !Junjie
                !MU_STG(IC, JC, KC, IDR) = &
                !( ( Viscousity(IC, JC, KC)  + Viscousity(IM, JC, KC) + &
                !    Viscousity0(IC, JC, KC) + Viscousity0(IM, JC, KC)) * YCL2ND_WFF(JJ) + &
                !  ( Viscousity(IC, JM, KC)  + Viscousity(IM, JM, KC) + &
                !    Viscousity0(IC, JM, KC) + Viscousity0(IM, JM, KC)) * YCL2ND_WFB(JJ) ) * XND2CL * 0.5_WP
            END DO
        END DO
    END DO

    !===================================== X -Z = (i', J, K') =============================
    ! NotICe: the brackets in below equation will introduce dIFfeRENces in the order of 1E-14.
    IDR = 2
    NYI = 1
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = NYI, N2DO(MYID)
            DO IC = 1, NCL1_io
                IM = IMV_io(IC)
                MU_STG(IC, JC, KC, IDR) = &
                ( Viscousity0(IC, JC, KC)  + Viscousity0(IM, JC, KC) + &
                Viscousity0(IC, JC, KM)  + Viscousity0(IM, JC, KM) ) * XND2CL * ZND2CL !Junjie
                !MU_STG(IC, JC, KC, IDR) = &
                !      ( Viscousity(IC, JC, KC)  + Viscousity(IM, JC, KC)  + &
                !        Viscousity0(IC, JC, KC) + Viscousity0(IM, JC, KC) + &
                !        Viscousity(IC, JC, KM)  + Viscousity(IM, JC, KM)  + &
                !        Viscousity0(IC, JC, KM) + Viscousity0(IM, JC, KM)  ) * XND2CL * ZND2CL * 0.5_WP
            END DO
        END DO
    END DO


    !=====================================Y-Z = (i, J', K') ============================
    IDR = 3
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO JC = NYI, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        DO KC = 1, NCL3
            KM = KMV(KC)
            DO IC = 1, NCL1_io
                MU_STG(IC, JC, KC, IDR) = &
                (( Viscousity0(IC, JC, KM)  + Viscousity0(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                ( Viscousity0(IC, JM, KM)  + Viscousity0(IC, JM, KC) ) * YCL2ND_WFB(JJ) ) * ZND2CL !Junjie

                !MU_STG(IC, JC, KC, IDR) = &
                !    (( Viscousity(IC, JC, KM)  + Viscousity(IC, JC, KC) + &
                !       Viscousity0(IC, JC, KM) + Viscousity0(IC, JC, KC) ) * YCL2ND_WFF(JJ) + &
                !     ( Viscousity(IC, JM, KM)  + Viscousity(IC, JM, KC)  + &
                !       Viscousity0(IC, JM, KM) + Viscousity0(IC, JM, KC) ) * YCL2ND_WFB(JJ) ) * ZND2CL * 0.5_WP
            END DO
        END DO
    END DO

    !=======================================================================================
    CALL INTFC_VARS3(NCL1S, NCL1E, NCL1S, NCL1E, MU_STG)

    IF(MYID == 0) THEN
        !MU_STG(:, 1, :, 1) = 0.5_WP * ( Viscousity0(:, 0, :) + Viscousity(:, 0, :) )
        !MU_STG(:, 1, :, 3) = 0.5_WP * ( Viscousity0(:, 0, :) + Viscousity(:, 0, :) )
        MU_STG(:, 1, :, 1) = Viscousity0(:, 0, :) !junjie
        MU_STG(:, 1, :, 3) = Viscousity0(:, 0, :)
    END IF

    IF(MYID == NPSLV) THEN
        !MU_STG(:, N2DO(MYID) + 1, :, 1) = 0.5_WP * ( Viscousity0(:, N2DO(MYID) + 1, :) + Viscousity(:, N2DO(MYID) + 1, :) )
        !MU_STG(:, N2DO(MYID) + 1, :, 3) = 0.5_WP * ( Viscousity0(:, N2DO(MYID) + 1, :) + Viscousity(:, N2DO(MYID) + 1, :) )
        MU_STG(:, N2DO(MYID) + 1, :, 1) = Viscousity0(:, N2DO(MYID) + 1, :) !Junjie
        MU_STG(:, N2DO(MYID) + 1, :, 3) = Viscousity0(:, N2DO(MYID) + 1, :) !Junjie
    END IF

    !CALL DEBUG_WRT_LOCAL(MU_STG(:, :, :, 1), 0, N2DO(MYID) + 1, 'MUG1')
    !CALL DEBUG_WRT_LOCAL(MU_STG(:, :, :, 2), 0, N2DO(MYID) + 1, 'MUG2')
    !CALL DEBUG_WRT_LOCAL(MU_STG(:, :, :, 3), 0, N2DO(MYID) + 1, 'MUG3')

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DENSITY_Staggered
    USE FLOW_INFO
    USE THERMAL_INFO
    USE MESH_INFO
    USE init_info
    IMPLICIT NONE
    INTEGER(4) :: IDR, NYI, IC, JC, KC, IM, JJ, KM, JM

    IF(iThermoDynamics /= 1) THEN
        D_STG = D_inlet
        RETURN
    END IF

    !===================================== X = (i', J, K) ==========================
    IDR = 1
    NYI = 1
    DO IC = 1, NCL1_io
        IM = IMV_io(IC)
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                D_STG(IC, JC, KC, IDR) = &
                ( DENSITY0(IC, JC, KC)  + DENSITY0(IM, JC, KC) ) * XND2CL
            END DO
        END DO
    END DO

    !===================================== Z === X = (i, J, K') =======================
    IDR = 3
    NYI = 1
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = NYI, N2DO(MYID)
            DO IC = 1, NCL1_io
                D_STG(IC, JC, KC, IDR) = &
                ( DENSITY0(IC, JC, KC)  + DENSITY0(IC, JC, KM) ) * ZND2CL
            END DO
        END DO
    END DO


    !=====================================Y ==== X = (i, J', K) =======================
    IDR = 2
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO JC = NYI, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_io
                D_STG(IC, JC, KC, IDR) = &
                YCL2ND_WFF(JJ) * DENSITY0(IC, JC, KC) + &
                YCL2ND_WFB(JJ) * DENSITY0(IC, JM, KC)
            END DO
        END DO
    END DO

    CALL INTFC_VARS3(1, NCL1_io, NCL1S, NCL1E,D_STG)

    IF(MYID == 0) THEN
        D_STG(:, 1, :, 2)           = DENSITY0(:, 0, :)
    END IF

    IF(MYID == NPSLV) THEN
        D_STG(:, N2DO(MYID) + 1, :, 2) = DENSITY0(:, N2DO(MYID) + 1, :)
    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE DENSITY_IMPLICIT_add
    USE FLOW_INFO
    USE THERMAL_INFO
    USE MESH_INFO
    USE init_info
    IMPLICIT NONE
    INTEGER(4) :: IDR, NYI, IC, JC, KC, IM, JJ, KM, JM
    REAL(WP) :: DEN_STG

    IF(iThermoDynamics /= 1) THEN
        DRHOI_STG = 0.0_WP
        RETURN
    END IF
    !===================================== X =============================
    IDR = 1
    NYI = 1
    DO IC = 1, NCL1_io
        IM = IMV_io(IC)
        DO JC = NYI, N2DO(MYID)
            DO KC = 1, NCL3
                DEN_STG = ( DENSITY0(IC, JC, KC) + DENSITY0(IM, JC, KC) ) * XND2CL
                DRHOI_STG(IC, JC, KC, IDR) = 1.0_WP / D_STG(IC, JC, KC, IDR) - 1.0_WP / DEN_STG
            END DO
        END DO
    END DO

    !===================================== Z =============================
    IDR = 3
    NYI = 1
    DO KC = 1, NCL3
        KM = KMV(KC)
        DO JC = NYI, N2DO(MYID)
            DO IC = 1, NCL1_io
                DEN_STG = ( DENSITY0(IC, JC, KC) + DENSITY0(IC, JC, KM) ) * ZND2CL
                DRHOI_STG(IC, JC, KC, IDR) = 1.0_WP / D_STG(IC, JC, KC, IDR) - 1.0_WP / DEN_STG
            END DO
        END DO
    END DO


    !=====================================Y =============================
    IDR = 2
    NYI = 1
    IF(MYID == 0) NYI = 2
    DO JC = NYI, N2DO(MYID)
        JJ = JCL2G(JC)
        JM = JLMV(JC)
        DO KC = 1, NCL3
            DO IC = 1, NCL1_io
                DEN_STG = YCL2ND_WFF(JJ) * DENSITY0(IC, JC, KC) + &
                YCL2ND_WFB(JJ) * DENSITY0(IC, JM, KC)
                DRHOI_STG(IC, JC, KC, IDR) = 1.0_WP / D_STG(IC, JC, KC, IDR) - 1.0_WP / DEN_STG
            END DO
        END DO
    END DO

    CALL INTFC_VARS3(1, NCL1_io, NCL1S, NCL1E, DRHOI_STG)


    IF(MYID == 0) THEN
        DO IC = 1, NCL1_io
            DRHOI_STG(IC, 1, :, 2) = 0.0_WP
        END DO
    END IF

    IF(MYID == NPSLV) THEN
        DO IC = 1, NCL1_io
            DRHOI_STG(IC, N2DO(MYID) + 1, :, 2) = 0.0_WP
        END DO
    END IF


    !test
    !        DO JC = 1, N2DO(MYID)
    !            DO IC = 1, NCL1_io
    !                DO KC = 1, Ncl3
    !                    IF(  DABS(DRHOI_STG(IC, JC, KC, 1)) > 1.0E-10_WP) &!
    !                    WRITE(*, *) MYID, JC, IC, KC, DRHOI_STG(IC, JC, KC, 1 : Ndv)
    !                END DO
    !            END DO
    !        END DO
    !CALL DEBUG_WRT_LOCAL(DRHOI_STG, 0, N2DO(MYID) + 1, 'dstg')

    RETURN
END SUBROUTINE
