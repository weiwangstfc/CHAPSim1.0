!**********************************************************************************************************************************
!> @brief
!>        to build up interface for thermo- PRoperties
!> @details
!> SUBROUTINE: INTFC_MFD_THERMAL_io(in MYID = all)
!> SUBROUTINE:INTFC_OUL_THERMAL_io(in MYID = all)
!> SUBROUTINE:INTFC_INL_THERMAL_io(in MYID = all)
!> SUBROUTINE:INTFC_MFD_DENSITY_io(in MYID = all), not used.
!> @note
!> @toDO
! REVISION HISTORY:
! 06/ 2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE INTFC_MFD_THERMAL_io
    !>      @NOTE:
    !>         TOP AND BOTTOM WALL HERE IS USING ZERO GRADIENT CONDITION.
    !>        THIS TREATMENT IS NOT ACCURATE/PHYSICALLY REASONABLE.
    !>           A B.C. SUBROUTINE IS USED TO TREAT THE TOP AND BOTTOM WALL.
    USE thermal_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I, IK
    INTEGER(4) :: K
    INTEGER(4) :: ITAG
    INTEGER(4) :: IDESF
    INTEGER(4) :: IDESB
    INTEGER(4) :: TRC_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NSZ
    REAL(WP) :: BSEN_F_io(NCL1_io * NCL3, NTHERMAL)
    REAL(WP) :: BSEN_L_io(NCL1_io * NCL3, NTHERMAL)
    REAL(WP) :: BREC_F_io(NCL1_io * NCL3, NTHERMAL)
    REAL(WP) :: BREC_L_io(NCL1_io * NCL3, NTHERMAL)

    BSEN_F_io = 0.0_WP
    BSEN_L_io = 0.0_WP
    BREC_F_io = 0.0_WP
    BREC_L_io = 0.0_WP
    NSZ = NCL1_io * NCL3 * NTHERMAL

    !================== 1 ENTHALPY ===================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 1) = ENTHALPY(I, 1,         K) !Y = Local Floor
            BSEN_L_io(IK, 1) = ENTHALPY(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !================== 2 DENSITY ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 2) = DENSITY(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK, 2) = DENSITY(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !================== 3 TEMPERATURE ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 3) = TEMPERATURE(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK, 3) = TEMPERATURE(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO


    !================== 4 Viscousity ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 4) = Viscousity(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK, 4) = Viscousity(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !==================5 THERMAL CONDUCTIVITY =========================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 5) = THERMCONDT(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK, 5) = THERMCONDT(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !==================6 DH ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK,6) = DH(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK,6) = DH(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !==================7 CP ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK,7) = HEATCAP(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK,7) = HEATCAP(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO


    !======= TRANSFERING DATA IN THE main PARTITIONS ====================================
    ITAG = 0
    IF ( (MYID >  0) .AND. (MYID < NPSLV) ) THEN
        IDESF = MYID + 1    !next MYID
        IDESB = MYID - 1    !last MYID

        !SEND Floor B.C. plane p,u, V,w from current MYID to IDESB= MYID - 1
        !RECEIVE Floor B.C. plane p,u, V,w from IDESF = MYID + 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)

        !SEND TOP B.C. plane p,u, V,w from current MYID to IDESF = MYID + 1
        !RECEIVE TOP B.C. plane p,u, V,w from IDESB= MYID - 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                   &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !Constructing of ghost cells y_id= 0 and N2NO+ 1, received from adjacent pARtitions.

        !================== 1 ENTHALPY ======================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                ENTHALPY(I, 0,           K) = BREC_L_io(IK, 1)
                ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 1)
            END DO
        END DO

        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DENSITY(I, 0,           K) = BREC_L_io(IK, 2)
                DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 2)
            END DO
        END DO

        !================== 3 TEMPERATURE ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                TEMPERATURE(I, 0,           K) = BREC_L_io(IK, 3)
                TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 3)
            END DO
        END DO

        !================== 4 Viscousity ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                Viscousity(I, 0,           K) = BREC_L_io(IK, 4)
                Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 4)
            END DO
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                THERMCONDT(I, 0,           K) = BREC_L_io(IK, 5)
                THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 5)
            END DO
        END DO

        !==================6 DH ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DH(I, 0, K)         = BREC_L_io(IK,6)
                DH(I, N2DO(MYID) + 1, K) = BREC_F_io(IK,6)
            END DO
        END DO

        !==================7 CP ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                HEATCAP(I, 0, K)         = BREC_L_io(IK,7)
                HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(IK,7)
            END DO
        END DO


    ENDIF

    !======= TRANSFERING DATA IN THE TOP WALL PARTITION ========================================
    IF(MYID == NPSLV) THEN
        IDESF = 0
        IDESB = MYID - 1

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,               &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !@note vAR(I, N2DO(MYID) + 1, K) not USEd actually.
        !================== 1 ENTHALPY ===========================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                ENTHALPY(I, 0, K)         = BREC_L_io(IK, 1)
                !ENTHALPY(I, N2DO(MYID) + 1, K) = 2.0_WP * ENTHALPY(I, N2DO(MYID),  K) - ENTHALPY(I, N2DO(MYID) - 1, K)
            END DO
        END DO

        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DENSITY(I, 0, K)         = BREC_L_io(IK, 2)
                !DENSITY(I, N2DO(MYID) + 1, K) = 2.0_WP * DENSITY(I, N2DO(MYID),  K) - DENSITY(I, N2DO(MYID) - 1, K)
            END DO
        END DO

        !================== 3 TEMPERATURE ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                TEMPERATURE(I, 0, K)         = BREC_L_io(IK, 3)
                !TEMPERATURE(I, N2DO(MYID) + 1, K) = 2.0_WP * TEMPERATURE(I, N2DO(MYID),  K) - TEMPERATURE(I, N2DO(MYID) - 1, K)
            END DO
        END DO

        !================== 4 Viscousity ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                Viscousity(I, 0, K)         = BREC_L_io(IK, 4)
                !Viscousity(I, N2DO(MYID) + 1, K) = 2.0_WP * Viscousity(I, N2DO(MYID),  K) - Viscousity(I, N2DO(MYID) - 1, K)
            END DO
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                THERMCONDT(I, 0, K)         = BREC_L_io(IK, 5)
                !THERMCONDT(I, N2DO(MYID) + 1, K) = 2.0_WP * THERMCONDT(I, N2DO(MYID),  K) - THERMCONDT(I, N2DO(MYID) - 1, K)
            END DO
        END DO


        !==================6 DH ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DH(I, 0, K)         = BREC_L_io(IK,6)
                !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
            END DO
        END DO

        !==================7 CP ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                HEATCAP(I, 0, K)         = BREC_L_io(IK,7)
                !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
            END DO
        END DO

    ENDIF

    !======= TRANSFERING DATA IN THE BOTTOM WALL PARTITION ========================================
    IF (MYID == 0) THEN
        IDESF = MYID + 1
        IDESB = NPSLV

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)

        !================== 1 ENTHALPY ===========================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !ENTHALPY(I, 0, K)         = 2.0_WP * ENTHALPY(I, 1, K) - ENTHALPY(I, 2, K)
                ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 1)
            END DO
        END DO


        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !DENSITY(I, 0, K)         = 2.0_WP * DENSITY(I, 1, K) - DENSITY(I, 2, K)
                DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 2)
            END DO
        END DO

        !================== 3 TEMPERATURE ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !TEMPERATURE(I, 0, K)         = 2.0_WP * TEMPERATURE(I, 1, K) - TEMPERATURE(I, 2, K)
                TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 3)
            END DO
        END DO

        !================== 4 Viscousity ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !Viscousity(I, 0, K)         = 2.0_WP * Viscousity(I, 1, K) - Viscousity(I, 2, K)
                Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 4)
            END DO
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !THERMCONDT(I, 0, K)         = 2.0_WP * THERMCONDT(I, 1, K) - THERMCONDT(I, 2, K)
                THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 5)
            END DO
        END DO

        !==================6 DH ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
                DH(I, N2DO(MYID) + 1, K) = BREC_F_io(IK,6)
            END DO
        END DO

        !==================7 CP ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
                HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(IK,7)
            END DO
        END DO

    ENDIF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INTFC_OUL_THERMAL_io
    USE thermal_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I
    INTEGER(4) :: K
    INTEGER(4) :: ITAG
    INTEGER(4) :: IDESF
    INTEGER(4) :: IDESB
    INTEGER(4) :: TRC_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NSZ
    REAL(WP) :: BSEN_F_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BSEN_L_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BREC_F_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BREC_L_io(1 * NCL3, NTHERMAL)

    !REAL(WP) :: H_tmp, M_tmp, K_tmp, D_tmp, T_tmp

    BSEN_F_io = 0.0_WP
    BSEN_L_io = 0.0_WP
    BREC_F_io = 0.0_WP
    BREC_L_io = 0.0_WP
    NSZ = 1 * NCL3 * NTHERMAL

    !================== 1 ENTHALPY ===================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K, 1) = ENTHALPY(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 1) = ENTHALPY(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !================== 2 DENSITY ==================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K, 2) = DENSITY(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 2) = DENSITY(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !================== 3 TEMPERATURE ==================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K, 3) = TEMPERATURE(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 3) = TEMPERATURE(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO


    !================== 4 Viscousity ==================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K, 4) = Viscousity(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 4) = Viscousity(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================5 THERMAL CONDUCTIVITY =========================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K, 5) = THERMCONDT(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 5) = THERMCONDT(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================6 DH ==================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K,6) = DH(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K,6) = DH(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================7 CP ==================================
    I = NCL1_io + 1
    DO K = 1, NCL3
        BSEN_F_io(K,7) = HEATCAP(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K,7) = HEATCAP(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO


    !======= TRANSFERING DATA IN THE main PARTITIONS ====================================
    ITAG = 0
    IF ( (MYID >  0) .AND. (MYID < NPSLV) ) THEN
        IDESF = MYID + 1    !next MYID
        IDESB = MYID - 1    !last MYID

        !SEND Floor B.C. plane p,u, V,w from current MYID to IDESB= MYID - 1
        !RECEIVE Floor B.C. plane p,u, V,w from IDESF = MYID + 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)

        !SEND TOP B.C. plane p,u, V,w from current MYID to IDESF = MYID + 1
        !RECEIVE TOP B.C. plane p,u, V,w from IDESB= MYID - 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                   &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !Constructing of ghost cells y_id= 0 and N2NO+ 1, received from adjacent pARtitions.

        !================== 1 ENTHALPY ======================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            ENTHALPY(I, 0,           K) = BREC_L_io(K, 1)
            ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 1)
        END DO

        !================== 2 DENSITY ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            DENSITY(I, 0,           K) = BREC_L_io(K, 2)
            DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 2)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            TEMPERATURE(I, 0,           K) = BREC_L_io(K, 3)
            TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 3)
        END DO

        !================== 4 Viscousity ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            Viscousity(I, 0,           K) = BREC_L_io(K, 4)
            Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 4)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = NCL1_io + 1
        DO K = 1, NCL3
            THERMCONDT(I, 0,           K) = BREC_L_io(K, 5)
            THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 5)
        END DO

        !==================6 DH ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            DH(I, 0, K)         = BREC_L_io(K,6)
            DH(I, N2DO(MYID) + 1, K) = BREC_F_io(K,6)
        END DO

        !==================7 CP ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            HEATCAP(I, 0, K)         = BREC_L_io(K,7)
            HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(K,7)
        END DO


    ENDIF

    !======= TRANSFERING DATA IN THE TOP WALL PARTITION ========================================
    IF(MYID == NPSLV) THEN
        IDESF = 0
        IDESB = MYID - 1

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,               &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !@note vAR(I, N2DO(MYID) + 1, K) not USEd actually.
        !================== 1 ENTHALPY ===========================
        I = NCL1_io + 1
        DO K = 1, NCL3
            ENTHALPY(I, 0, K)         = BREC_L_io(K, 1)
            !ENTHALPY(I, N2DO(MYID) + 1, K) = 2.0_WP * ENTHALPY(I, N2DO(MYID),  K) - ENTHALPY(I, N2DO(MYID) - 1, K)
        END DO

        !================== 2 DENSITY ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            DENSITY(I, 0, K)         = BREC_L_io(K, 2)
            !DENSITY(I, N2DO(MYID) + 1, K) = 2.0_WP * DENSITY(I, N2DO(MYID),  K) - DENSITY(I, N2DO(MYID) - 1, K)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            TEMPERATURE(I, 0, K)         = BREC_L_io(K, 3)
            !TEMPERATURE(I, N2DO(MYID) + 1, K) = 2.0_WP * TEMPERATURE(I, N2DO(MYID),  K) - TEMPERATURE(I, N2DO(MYID) - 1, K)
        END DO

        !================== 4 Viscousity ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            Viscousity(I, 0, K)         = BREC_L_io(K, 4)
            !Viscousity(I, N2DO(MYID) + 1, K) = 2.0_WP * Viscousity(I, N2DO(MYID),  K) - Viscousity(I, N2DO(MYID) - 1, K)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = NCL1_io + 1
        DO K = 1, NCL3
            THERMCONDT(I, 0, K)         = BREC_L_io(K, 5)
            !THERMCONDT(I, N2DO(MYID) + 1, K) = 2.0_WP * THERMCONDT(I, N2DO(MYID),  K) - THERMCONDT(I, N2DO(MYID) - 1, K)
        END DO


        !==================6 DH ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            DH(I, 0, K)         = BREC_L_io(K,6)
            !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
        END DO

        !==================7 CP ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            HEATCAP(I, 0, K)         = BREC_L_io(K,7)
            !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
        END DO

    ENDIF

    !======= TRANSFERING DATA IN THE BOTTOM WALL PARTITION ========================================
    IF (MYID == 0) THEN
        IDESF = MYID + 1
        IDESB = NPSLV

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)

        !================== 1 ENTHALPY ===========================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !ENTHALPY(I, 0, K)         = 2.0_WP * ENTHALPY(I, 1, K) - ENTHALPY(I, 2, K)
            ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 1)
        END DO


        !================== 2 DENSITY ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !DENSITY(I, 0, K)         = 2.0_WP * DENSITY(I, 1, K) - DENSITY(I, 2, K)
            DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 2)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !TEMPERATURE(I, 0, K)         = 2.0_WP * TEMPERATURE(I, 1, K) - TEMPERATURE(I, 2, K)
            TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 3)
        END DO

        !================== 4 Viscousity ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !Viscousity(I, 0, K)         = 2.0_WP * Viscousity(I, 1, K) - Viscousity(I, 2, K)
            Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 4)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !THERMCONDT(I, 0, K)         = 2.0_WP * THERMCONDT(I, 1, K) - THERMCONDT(I, 2, K)
            THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 5)
        END DO

        !==================6 DH ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
            DH(I, N2DO(MYID) + 1, K) = BREC_F_io(K,6)
        END DO

        !==================7 CP ==================================
        I = NCL1_io + 1
        DO K = 1, NCL3
            !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
            HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(K,7)
        END DO

    ENDIF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INTFC_INL_THERMAL_io
    USE thermal_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I
    INTEGER(4) :: K
    INTEGER(4) :: ITAG
    INTEGER(4) :: IDESF
    INTEGER(4) :: IDESB
    INTEGER(4) :: TRC_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NSZ
    REAL(WP) :: BSEN_F_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BSEN_L_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BREC_F_io(1 * NCL3, NTHERMAL)
    REAL(WP) :: BREC_L_io(1 * NCL3, NTHERMAL)

    !REAL(WP) :: H_tmp, M_tmp, K_tmp, D_tmp, T_tmp

    BSEN_F_io = 0.0_WP
    BSEN_L_io = 0.0_WP
    BREC_F_io = 0.0_WP
    BREC_L_io = 0.0_WP
    NSZ = 1 * NCL3 * NTHERMAL

    !================== 1 ENTHALPY ===================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K, 1) = ENTHALPY(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 1) = ENTHALPY(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !================== 2 DENSITY ==================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K, 2) = DENSITY(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 2) = DENSITY(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !================== 3 TEMPERATURE ==================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K, 3) = TEMPERATURE(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 3) = TEMPERATURE(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO


    !================== 4 Viscousity ==================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K, 4) = Viscousity(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 4) = Viscousity(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================5 THERMAL CONDUCTIVITY =========================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K, 5) = THERMCONDT(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K, 5) = THERMCONDT(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================6 DH ==================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K,6) = DH(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K,6) = DH(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO

    !==================7 CP ==================================
    I = 0
    DO K = 1, NCL3
        BSEN_F_io(K,7) = HEATCAP(I, 1,         K)  !Y = Local Floor
        BSEN_L_io(K,7) = HEATCAP(I, N2DO(MYID), K)  !Y = Local ROOF
    END DO


    !======= TRANSFERING DATA IN THE main PARTITIONS ====================================
    ITAG = 0
    IF ( (MYID >  0) .AND. (MYID < NPSLV) ) THEN
        IDESF = MYID + 1    !next MYID
        IDESB = MYID - 1    !last MYID

        !SEND Floor B.C. plane p,u, V,w from current MYID to IDESB= MYID - 1
        !RECEIVE Floor B.C. plane p,u, V,w from IDESF = MYID + 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)

        !SEND TOP B.C. plane p,u, V,w from current MYID to IDESF = MYID + 1
        !RECEIVE TOP B.C. plane p,u, V,w from IDESB= MYID - 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                   &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !Constructing of ghost cells y_id= 0 and N2NO+ 1, received from adjacent pARtitions.

        !================== 1 ENTHALPY ======================================
        I = 0
        DO K = 1, NCL3
            ENTHALPY(I, 0,           K) = BREC_L_io(K, 1)
            ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 1)
        END DO

        !================== 2 DENSITY ==================================
        I = 0
        DO K = 1, NCL3
            DENSITY(I, 0,           K) = BREC_L_io(K, 2)
            DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 2)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = 0
        DO K = 1, NCL3
            TEMPERATURE(I, 0,           K) = BREC_L_io(K, 3)
            TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 3)
        END DO

        !================== 4 Viscousity ==================================
        I = 0
        DO K = 1, NCL3
            Viscousity(I, 0,           K) = BREC_L_io(K, 4)
            Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 4)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = 0
        DO K = 1, NCL3
            THERMCONDT(I, 0,           K) = BREC_L_io(K, 5)
            THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 5)
        END DO

        !==================6 DH ==================================
        I = 0
        DO K = 1, NCL3
            DH(I, 0, K)         = BREC_L_io(K,6)
            DH(I, N2DO(MYID) + 1, K) = BREC_F_io(K,6)
        END DO

        !==================7 CP ==================================
        I = 0
        DO K = 1, NCL3
            HEATCAP(I, 0, K)         = BREC_L_io(K,7)
            HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(K,7)
        END DO


    ENDIF

    !======= TRANSFERING DATA IN THE TOP WALL PARTITION ========================================
    IF(MYID == NPSLV) THEN
        IDESF = 0
        IDESB = MYID - 1

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,               &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !@note vAR(I, N2DO(MYID) + 1, K) not USEd actually.
        !================== 1 ENTHALPY ===========================
        I = 0
        DO K = 1, NCL3
            ENTHALPY(I, 0, K)         = BREC_L_io(K, 1)
            !ENTHALPY(I, N2DO(MYID) + 1, K) = 2.0_WP * ENTHALPY(I, N2DO(MYID),  K) - ENTHALPY(I, N2DO(MYID) - 1, K)
        END DO

        !================== 2 DENSITY ==================================
        I = 0
        DO K = 1, NCL3
            DENSITY(I, 0, K)         = BREC_L_io(K, 2)
            !DENSITY(I, N2DO(MYID) + 1, K) = 2.0_WP * DENSITY(I, N2DO(MYID),  K) - DENSITY(I, N2DO(MYID) - 1, K)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = 0
        DO K = 1, NCL3
            TEMPERATURE(I, 0, K)         = BREC_L_io(K, 3)
            !TEMPERATURE(I, N2DO(MYID) + 1, K) = 2.0_WP * TEMPERATURE(I, N2DO(MYID),  K) - TEMPERATURE(I, N2DO(MYID) - 1, K)
        END DO

        !================== 4 Viscousity ==================================
        I = 0
        DO K = 1, NCL3
            Viscousity(I, 0, K)         = BREC_L_io(K, 4)
            !Viscousity(I, N2DO(MYID) + 1, K) = 2.0_WP * Viscousity(I, N2DO(MYID),  K) - Viscousity(I, N2DO(MYID) - 1, K)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = 0
        DO K = 1, NCL3
            THERMCONDT(I, 0, K)         = BREC_L_io(K, 5)
            !THERMCONDT(I, N2DO(MYID) + 1, K) = 2.0_WP * THERMCONDT(I, N2DO(MYID),  K) - THERMCONDT(I, N2DO(MYID) - 1, K)
        END DO


        !==================6 DH ==================================
        I = 0
        DO K = 1, NCL3
            DH(I, 0, K)         = BREC_L_io(K,6)
            !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
        END DO

        !==================7 CP ==================================
        I = 0
        DO K = 1, NCL3
            HEATCAP(I, 0, K)         = BREC_L_io(K,7)
            !DH(I, N2DO(MYID) + 1, K) = 2.0_WP * DH(I, N2DO(MYID),  K) - DH(I, N2DO(MYID) - 1, K)
        END DO

    ENDIF

    !======= TRANSFERING DATA IN THE BOTTOM WALL PARTITION ========================================
    IF (MYID == 0) THEN
        IDESF = MYID + 1
        IDESB = NPSLV

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)

        !================== 1 ENTHALPY ===========================
        I = 0
        DO K = 1, NCL3
            !ENTHALPY(I, 0, K)         = 2.0_WP * ENTHALPY(I, 1, K) - ENTHALPY(I, 2, K)
            ENTHALPY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 1)
        END DO


        !================== 2 DENSITY ==================================
        I = 0
        DO K = 1, NCL3
            !DENSITY(I, 0, K)         = 2.0_WP * DENSITY(I, 1, K) - DENSITY(I, 2, K)
            DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 2)
        END DO

        !================== 3 TEMPERATURE ==================================
        I = 0
        DO K = 1, NCL3
            !TEMPERATURE(I, 0, K)         = 2.0_WP * TEMPERATURE(I, 1, K) - TEMPERATURE(I, 2, K)
            TEMPERATURE(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 3)
        END DO

        !================== 4 Viscousity ==================================
        I = 0
        DO K = 1, NCL3
            !Viscousity(I, 0, K)         = 2.0_WP * Viscousity(I, 1, K) - Viscousity(I, 2, K)
            Viscousity(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 4)
        END DO

        !==================5 THERMAL CONDUCTIVITY ===========================
        I = 0
        DO K = 1, NCL3
            !THERMCONDT(I, 0, K)         = 2.0_WP * THERMCONDT(I, 1, K) - THERMCONDT(I, 2, K)
            THERMCONDT(I, N2DO(MYID) + 1, K) = BREC_F_io(K, 5)
        END DO

        !==================6 DH ==================================
        I = 0
        DO K = 1, NCL3
            !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
            DH(I, N2DO(MYID) + 1, K) = BREC_F_io(K,6)
        END DO

        !==================7 CP ==================================
        I = 0
        DO K = 1, NCL3
            !DH(I, 0, K)         = 2.0_WP * DH(I, 1, K) - DH(I, 2, K)
            HEATCAP(I, N2DO(MYID) + 1, K) = BREC_F_io(K,7)
        END DO

    ENDIF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE INTFC_MFD_DENSITY_io ! not used
    !>      @NOTE:
    !>         TOP AND BOTTOM WALL HERE IS USING ZERO GRADIENT CONDITION.
    !>        THIS TREATMENT IS NOT ACCURATE/PHYSICALLY REASONABLE.
    !>           A B.C. SUBROUTINE IS USED TO TREAT THE TOP AND BOTTOM WALL.
    USE thermal_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: I, IK
    INTEGER(4) :: K
    INTEGER(4) :: ITAG
    INTEGER(4) :: IDESF
    INTEGER(4) :: IDESB
    INTEGER(4) :: TRC_STS(MPI_STATUS_SIZE)
    INTEGER(4) :: NSZ
    REAL(WP) :: BSEN_F_io(NCL1_io * NCL3, 1)
    REAL(WP) :: BSEN_L_io(NCL1_io * NCL3, 1)
    REAL(WP) :: BREC_F_io(NCL1_io * NCL3, 1)
    REAL(WP) :: BREC_L_io(NCL1_io * NCL3, 1)

    BSEN_F_io = 0.0_WP
    BSEN_L_io = 0.0_WP
    BREC_F_io = 0.0_WP
    BREC_L_io = 0.0_WP
    NSZ = NCL1_io * NCL3 * 1


    !================== 2 DENSITY ==================================
    DO I = 1, NCL1_io
        DO K = 1, NCL3
            IK = (I - 1) * NCL3 + K
            BSEN_F_io(IK, 1) = DENSITY(I, 1,         K)  !Y = Local Floor
            BSEN_L_io(IK, 1) = DENSITY(I, N2DO(MYID), K)  !Y = Local ROOF
        END DO
    END DO

    !======= TRANSFERING DATA IN THE main PARTITIONS ====================================
    ITAG = 0
    IF ( (MYID >  0) .AND. (MYID < NPSLV) ) THEN
        IDESF = MYID + 1    !next MYID
        IDESB = MYID - 1    !last MYID

        !SEND Floor B.C. plane p,u, V,w from current MYID to IDESB= MYID - 1
        !RECEIVE Floor B.C. plane p,u, V,w from IDESF = MYID + 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)

        !SEND TOP B.C. plane p,u, V,w from current MYID to IDESF = MYID + 1
        !RECEIVE TOP B.C. plane p,u, V,w from IDESB= MYID - 1 to current MYID
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                   &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !Constructing of ghost cells y_id= 0 and N2NO+ 1, received from adjacent pARtitions.
        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DENSITY(I, 0,           K) = BREC_L_io(IK, 1)
                DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 1)
            END DO
        END DO

    ENDIF

    !======= TRANSFERING DATA IN THE TOP WALL PARTITION ========================================
    IF(MYID == NPSLV) THEN
        IDESF = 0
        IDESB = MYID - 1

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,               &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ, &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !@note vAR(I, N2DO(MYID) + 1, K) not USEd actually.
        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DENSITY(I, 0, K)         = BREC_L_io(IK, 1)
                DENSITY(I, N2DO(MYID) + 1, K) = 2.0_WP * DENSITY(I, N2DO(MYID),  K) - DENSITY(I, N2DO(MYID) - 1, K)
            END DO
        END DO

    ENDIF

    !======= TRANSFERING DATA IN THE BOTTOM WALL PARTITION ========================================
    IF (MYID == 0) THEN
        IDESF = MYID + 1
        IDESB = NPSLV

        CALL  MPI_SENDRECV(BSEN_F_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, BREC_F_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, ICOMM, TRC_STS, IERROR)
        CALL  MPI_SENDRECV(BSEN_L_io(1, 1), NSZ,                 &
        MPI_DOUBLE_PRECISION, IDESF, ITAG, BREC_L_io(1, 1), NSZ,  &
        MPI_DOUBLE_PRECISION, IDESB, ITAG, ICOMM, TRC_STS, IERROR)


        !================== 2 DENSITY ==================================
        DO I = 1, NCL1_io
            DO K = 1, NCL3
                IK = (I - 1) * NCL3 + K
                DENSITY(I, 0, K)         = 2.0_WP * DENSITY(I, 1, K) - DENSITY(I, 2, K)
                DENSITY(I, N2DO(MYID) + 1, K) = BREC_F_io(IK, 1)
            END DO
        END DO

    ENDIF

    RETURN
END SUBROUTINE
