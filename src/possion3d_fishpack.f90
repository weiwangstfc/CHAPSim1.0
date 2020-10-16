!**********************************************************************************************************************************
!> @brief
!>        to prepARe information required to USE FISHPACK librARy
!> @details
!> module: FISHPACK_POIS3D_INFO
!>         variables for coefficients
!> SUBROUTINE: FISHPACK_POIS3D_INIT (in MYID = all)
!>             - CALLed by the main solver
!> SUBROUTINE: FFTPACK_ROOT (in MYID = all)
!>             - CALLed by FISHPACK_POIS3D_INIT
!> SUBROUTINE: FISHPACK_POIS3D_SIMPLE (in MYID = all)
!>             - CALLed by the main solver SOLVERRK3_MOM_io
!> SUBROUTINE: FFTPACK_XZ (in MYID = all)
!>             - CALLed by FISHPACK_POIS3D_SIMPLE
!> SUBROUTINE: TRID0 (in MYID = all)
!>             - CALLed by FISHPACK_POIS3D_SIMPLE
!> @note
!> @toDO
! REVISION HISTORY:
! 04/ 2014- created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE FISHPACK_POIS3D_INFO
    USE WPRECISION

    INTEGER(4) :: LPEROD, MPEROD, NPEROD
    INTEGER(4) :: LP, MP, NP
    INTEGER(4) :: L, M, NG, NL, ML
    INTEGER(4) :: LDIMF, MDIMF
    INTEGER(4) :: WSZ

    REAL(WP) :: C1, C2
    REAL(WP) :: SCALX, SCALY
    REAL(WP), ALLOCATABLE :: A(:), B(:), C(:), D(:), BB(:)
    REAL(WP), ALLOCATABLE :: XRT(:), YRT(:), WX(:), WY(:)
    REAL(WP), ALLOCATABLE :: FR(:, :, :), FK(:, :, :), T(:)
    REAL(WP), ALLOCATABLE :: F_io   (:, :, :)

END MODULE

!**********************************************************************************************************************************
SUBROUTINE FISHPACK_POIS3D_INIT
    !    SHOULD BE CALLED IN SLAVES BEFORE FIRSTLY
    USE mesh_info
    USE flow_info
    USE init_info
    USE FISHPACK_POIS3D_INFO
    IMPLICIT NONE

    INTEGER(4) :: K


    IF(BCX_io(1) == 3 .AND. BCX_io(2) == 3 )  LPEROD = 0
    IF(BCX_io(1) == 1 .AND. BCX_io(2) == 1 )  LPEROD = 1
    IF(BCX_io(1) == 1 .AND. BCX_io(2) == 2 )  LPEROD = 2
    IF(BCX_io(1) == 2 .AND. BCX_io(2) == 2 )  LPEROD = 3
    IF(BCX_io(1) == 2 .AND. BCX_io(2) == 1 )  LPEROD = 4

    IF(BCZ(1) == 3 .AND. BCZ(2) == 3 )  MPEROD = 0
    IF(BCZ(1) == 1 .AND. BCZ(2) == 1 )  MPEROD = 1
    IF(BCZ(1) == 1 .AND. BCZ(2) == 2 )  MPEROD = 2
    IF(BCZ(1) == 2 .AND. BCZ(2) == 2 )  MPEROD = 3
    IF(BCZ(1) == 2 .AND. BCZ(2) == 1 )  MPEROD = 4

    NPEROD = 1

    LP = LPEROD + 1
    MP = MPEROD + 1
    NP = NPEROD + 1

    L  = NCL1_io
    M  = NCL3
    NG = NCL2
    NL = N2DO(MYID)
    ML = N3DO(MYID)

    IF( (L <= 3) .OR. (M <= 3) .OR. (NG <= 3) ) &
    CALL ERRHDL('Dimensions in poISson solver should be lARger than 3!', MYID)

    C1 = DXQI
    C2 = DZQI

    WSZ = 30 + L + M + 2 * NG + MAX(L, M, NG) + &
    7* (INT((L+ 1) / 2) + INT((M + 1) / 2)) + 128

    ALLOCATE (A (NG) ) ; A = 0.0_WP
    ALLOCATE (B (NG) ) ; B = 0.0_WP
    ALLOCATE (C (NG) ) ; C = 0.0_WP
    ALLOCATE (D (NG) ) ; D = 0.0_WP
    ALLOCATE (BB(NG) ) ; BB = 0.0_WP

    MEMPC_Byte = MEMPC_Byte + NG*10 *8

    ALLOCATE (XRT(L))  ; XRT = 0.0_WP
    ALLOCATE (YRT(M))  ; YRT = 0.0_WP
    ALLOCATE (WX(WSZ)) ; WX = 0.0_WP
    ALLOCATE (WY(WSZ)) ; WY = 0.0_WP
    ALLOCATE (FR(L, M, NL) ) ; FR = 0.0_WP
    ALLOCATE (FK(L, ML, NG) ) ; FK = 0.0_WP
    ALLOCATE (T(MAX0(L, M, NG)) ) ; T = 0.0_WP

    ALLOCATE ( F_io   (NCL1_io, NCL2, N3DO(0) )     )       ;  F_io = 0.0_WP

    MEMPC_Byte = MEMPC_Byte + (L+M + WSZ*2 +L * M * NL+L * ML * NG+MAX0(L, M, NG)) *8

    DO K = 1, NG
        A(K) = AMPH(K)
        B(K) = ACPH(K)
        C(K) = APPH(K)
    END DO

    CALL FFTPACK_ROOT

    RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE FFTPACK_ROOT
!       CALLED IN SLAVES AFTER FFT_INITIALIZATION

        USE FISHPACK_POIS3D_INFO
        IMPLICIT NONE

        REAL(WP) :: PI
        REAL(WP) :: DUM
        REAL(WP) :: DY, DJ
        REAL(WP),external :: PIMACH
        REAL(WP) :: DX, DI

        INTEGER(4) :: I, J
        INTEGER(4) :: LRDEL
        INTEGER(4) :: MRDEL
        INTEGER(4) :: LR, MR

        PI = PIMACH(DUM)
        LR = L
        MR = M
!C
!C     GENERATE TRANSFORM ROOTS FOR X direction
!C
        LRDEL = ((LP - 1) * (LP - 3) * (LP - 5)) / 3
        SCALX = DBLE(LR + LRDEL)
        DX = PI / (2.0_WP * SCALX)
        GO TO (108, 103, 101, 102, 101), LP
  101   DI = 0.50_WP
        SCALX = 2.0_WP * SCALX
        GO TO 104
  102   DI = 1.00_WP
        GO TO 104
  103   DI = 0.00_WP
  104   DO 105 I = 1, LR
            XRT(I) = -4.0_WP * C1 * (SIN((DBLE(I) - DI) * DX))**2
  105   CONTINUE
        SCALX = 2.0_WP * SCALX
        GO TO (112, 106, 110, 107, 111), LP
  106   CALL SINTI (LR, WX)
        GO TO 112
  107   CALL COSTI (LR, WX)
        GO TO 112
  108   XRT(1) = 0.0_WP
        XRT(LR) = -4.0_WP * C1
        DO 109 I = 3,LR, 2
            XRT(I - 1) = -4.0_WP * C1 * (SIN(DBLE((I - 1)) * DX))**2
            XRT(I) = XRT(I - 1)
  109   CONTINUE
        CALL RFFTI (LR, WX)
        GO TO 112
  110   CALL SINQI (LR, WX)
        GO TO 112
  111   CALL COSQI (LR, WX)
  112   CONTINUE

!C
!C     GENERATE TRANSFORM ROOTS FOR Y direction (Z IN REAL)
!C
      MRDEL = ((MP - 1) * (MP - 3) * (MP - 5)) / 3
      SCALY = DBLE(MR + MRDEL)
      DY = PI / (2.0_WP * SCALY)
      GO TO (120, 115, 113, 114, 113), MP
  113 DJ = 0.50_WP
      SCALY = 2.0_WP * SCALY
      GO TO 116
  114 DJ = 1.00_WP
      GO TO 116
  115 DJ = 0.00_WP
  116 DO 117 J = 1, MR
         YRT(J) = -4.0_WP * C2 * (SIN((DBLE(J) - DJ) * DY))**2
  117 CONTINUE
      SCALY = 2.0_WP * SCALY
      GO TO (124, 118, 122, 119, 123), MP
  118 CALL SINTI (MR, WY)
      GO TO 124
  119 CALL COSTI (MR, WY)
      GO TO 124
  120 YRT(1) = 0.0_WP
      YRT(MR) = -4.0_WP * C2
      DO 121 J = 3, MR, 2
         YRT(J - 1) = -4.0_WP * C2 * (SIN(DBLE((J - 1)) * DY))**2
         YRT(J) = YRT(J - 1)
  121 CONTINUE
      CALL RFFTI (MR, WY)
      GO TO 124
  122 CALL SINQI (MR, WY)
      GO TO 124
  123 CALL COSQI (MR, WY)
  124 CONTINUE

     RETURN

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE FISHPACK_POIS3D_SIMPLE
    !    CALLED IN SLAVES EVERY RK STAGE TO CALCULATE pressure CORRECTION TERMS.

    USE FISHPACK_POIS3D_INFO
    USE mesh_info
    USE flow_info
    USE init_info
    IMPLICIT NONE

    INTEGER(4) :: IFWRD
    INTEGER(4) :: I, J, K, JJ

    !   ======== RECONSTRUCT RHS TO FIT POIS3D.========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                FR(I, K, J) = RHSLLPHI_io(I, J, K)
            END DO
        END DO
    END DO

    !   ========forward FFT IN X AND Z direction ========
    IFWRD = 1
    CALL FFTPACK_XZ(IFWRD)

    !   ======== RESTORE RHSLLPHI ========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                RHSLLPHI_io(I, J, K) = FR(I, K, J)
            END DO
        END DO
    END DO

    !   ======== TRANSPORT Y- DECOMP TO K - DECOMP ========
    !CALL TRASP_Y2Z_RHSLLPHI_io
    CALL TRASP23_Y2Z(NCL1_io, 1, N2DO(0), RHSLLPHI_io, F_io)

    !   ========CONSTRUCT DATA FOR TDMA========
    DO I = 1, L
        DO K = 1, N3DO(MYID)
            DO J = 1, NG
                FK(I, K, J) = F_io(I, J, K)
            END DO
        END DO
    END DO

    !   ======== TDMA IN Y direction FOR PART OF FR(:,PART, :) ========
    DO I = 1, L
        DO J = 1, N3DO(MYID)
            JJ = KCL2G(J)
            DO K = 1, NG
                BB(K) = B(K) + XRT(I) / RCCI2(K) + YRT(JJ)
                T(K) = FK(I, J, K)
            END DO

            CALL TRID0

            DO K = 1, NG
                FK(I, J, K) = T(K)
            END DO

        END DO
    END DO

    !      ======== RE -CONSTRUCT DATA BACK ========
    DO I = 1, L
        DO K = 1, N3DO(MYID)
            DO J = 1, NG
                F_io(I, J, K) = FK(I, K, J)
            END DO
        END DO
    END DO

    !   ======== TRANSPORT Z - DECOMP TO Y- DECOMP ========
    !CALL TRASP_Z2Y_RHSLLPHI_io
    CALL TRASP23_Z2Y(NCL1_io, 1, N2DO(0), RHSLLPHI_io, F_io)

    !   ======== RESTORE RHSLLPHI ========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                FR(I, K, J) = RHSLLPHI_io(I, J, K)
            END DO
        END DO
    END DO

    !   ========backward FFT IN Z AND X directionS ========
    IFWRD = 2
    CALL FFTPACK_XZ(IFWRD)

    !   ======== SCALE THE CALCULATED VALUE ========
    DO I = 1, L
        DO J = 1, M
            DO K = 1, NL
                FR(I, J, K) = FR(I, J, K) / (SCALX * SCALY)
            END DO
        END DO
    END DO


    !   ======== RE -STORE AND ASSIGN DATA TO DPH ========
    DPH_io = 0.0_WP
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                DPH_io(I, J, K) = FR(I, K, J)
            END DO
        END DO
    END DO

    RETURN

END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE FFTPACK_XZ(IFWRD)
!     FOR FFT TRANSFORM FROM SPACE TO WAVENUMBER, IFWRD = 1, IS = 1
!     FOR FFT TRANSFORM FROM WAVENUMBER TO SPACE, IFWRD = 2, IS = -1
      USE FISHPACK_POIS3D_INFO
      IMPLICIT NONE
      INTEGER(4) :: IFWRD
      INTEGER(4) :: I, J, K
      INTEGER(4) :: LR, MR, NR



      LR = L
      MR = M
      NR = NL

      GO TO(125, 142) IFWRD

  125 CONTINUE
!
!     TRANSFORM X
!
      DO 141 J = 1, MR
         DO 140 k = 1, NR
            DO 126 I = 1, LR
               T(I) = FR(I, J, K)
  126       CONTINUE
            GO TO (127, 130, 131, 134, 135), LP
  127       GO TO (128, 129), IFWRD
  128       CALL RFFTF (LR, T,WX)
            GO TO 138
  129       CALL RFFTB (LR, T,WX)
            GO TO 138
  130       CALL SINT (LR, T,WX)
            GO TO 138
  131       GO TO (132, 133), IFWRD
  132       CALL SINQF (LR, T,WX)
            GO TO 138
  133       CALL SINQB (LR, T,WX)
            GO TO 138
  134       CALL COST (LR, T,WX)
            GO TO 138
  135       GO TO (136, 137), IFWRD
  136       CALL COSQF (LR, T,WX)
            GO TO 138
  137       CALL COSQB (LR, T,WX)
  138       CONTINUE
            DO 139 I = 1, LR
               FR(I, J, K) = T(I)
  139       CONTINUE
  140    CONTINUE
  141 CONTINUE
      GO TO (142, 159), IFWRD


!C
!C     TRANSFORM Y
!C
  142 CONTINUE
      DO 158 I = 1, LR
         DO 157 k = 1, NR
            DO 143 J = 1, MR
               T(J) = FR(I, J, K)
  143       CONTINUE
            GO TO (144, 147, 148, 151, 152), MP
  144       GO TO (145, 146), IFWRD
  145       CALL RFFTF (MR, T,WY)
            GO TO 155
  146       CALL RFFTB (MR, T,WY)
            GO TO 155
  147       CALL SINT (MR, T,WY)
            GO TO 155
  148       GO TO (149, 150), IFWRD
  149       CALL SINQF (MR, T,WY)
            GO TO 155
  150       CALL SINQB (MR, T,WY)
            GO TO 155
  151       CALL COST (MR, T,WY)
            GO TO 155
  152       GO TO (153, 154), IFWRD
  153       CALL COSQF (MR, T,WY)
            GO TO 155
  154       CALL COSQB (MR, T,WY)
  155       CONTINUE
            DO 156 J = 1, MR
               FR(I, J, K) = T(J)
  156       CONTINUE
  157    CONTINUE
  158 CONTINUE
      GO TO (159, 125), IFWRD

  159 CONTINUE


    RETURN

END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE TRID0
      USE FISHPACK_POIS3D_INFO
      IMPLICIT NONE

      INTEGER(4) :: NR, MM1
      REAL(WP) :: Z
      INTEGER(4) :: I, IP

      NR = NG
      MM1 = NR- 1
      Z = 1.0_WP / BB(1)
      D(1) = C(1) * Z
      T(1) = T(1) * Z
      DO 101 I = 2, MM1
         Z = 1.0_WP / (BB(I) - A(I) * D(I - 1))
         D(I) = C(I) * Z
         T(I) = (T(I) - A(I) * T(I - 1)) * Z
  101 CONTINUE
      Z = BB(NR) - A(NR) * D(MM1)
      IF (Z /= 0.0_WP) GO TO 102
      T(NR) = 0.0_WP
      GO TO 103
  102 T(NR) = (T(NR) - A(NR) * T(MM1)) /Z
  103 CONTINUE
      DO 104 IP = 1, MM1
         I = NR-IP
         T(I) = T(I) - D(I) * T(I + 1)
  104 CONTINUE
      RETURN

END SUBROUTINE
