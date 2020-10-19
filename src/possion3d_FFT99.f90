!**********************************************************************************************************************************
!> @brief
!>        to prepare information required to USE FFT99 library
!> @details
!> module: FFT99_info1
!>         variables for coefficients
!> module: FFT99_info2
!>         variables for parameters
!> SUBROUTINE: FFT99_POIS3D_INIT (in MYID = all)
!> SUBROUTINE: FFT99_ROOT (in MYID = 0)
!> SUBROUTINE: BCAST_FFT99_ROOT (in MYID = all)
!> SUBROUTINE: TDMA_PHI_COEF (in MYID = all)
!> SUBROUTINE: FFT99_POIS3D_periodicxz (in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010- Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014- Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE FFT99_info1
    USE WPRECISION

    REAL(WP), ALLOCATABLE :: AK1(:)
    REAL(WP), ALLOCATABLE :: AK3(:)

    REAL(WP), ALLOCATABLE :: TRIGXX1(:)
    REAL(WP), ALLOCATABLE :: TRIGXX2(:)
    REAL(WP), ALLOCATABLE :: TRIGXX3(:)

    REAL(WP), ALLOCATABLE :: TRIGXC1(:)
    REAL(WP), ALLOCATABLE :: TRIGXC2(:)
    REAL(WP), ALLOCATABLE :: TRIGXC3(:)

    INTEGER(4) :: IFXX1(13)
    INTEGER(4) :: IFXX2(13)
    INTEGER(4) :: IFXX3(13)
    INTEGER(4) :: IFXC3(13)
    INTEGER(4) :: IFXC1(13)
    INTEGER(4) :: IFXC2(13)

END MODULE

!**********************************************************************************************************************************
MODULE FFT99_info2
    USE WPRECISION

    REAL(WP), ALLOCATABLE :: XR  (:, :)
    REAL(WP), ALLOCATABLE :: WORK(:, :)
    COMPLEX(WP), ALLOCATABLE :: XA  (:, :)
    COMPLEX(WP), ALLOCATABLE :: WOR (:, :)

    INTEGER(4) :: LHFP

    INTEGER(4) :: L, M, N, ML, NL, MLmax, NLmax

    REAL(WP) :: DXQI0, DZQI0

    REAL(WP), ALLOCATABLE :: AMJP(:, :, :)
    REAL(WP), ALLOCATABLE :: ACJP(:, :, :)
    REAL(WP), ALLOCATABLE :: APJP(:, :, :)

    REAL(WP), ALLOCATABLE :: RHSLLPHIRe(:, :, :)
    REAL(WP), ALLOCATABLE :: RHSLLPHIIm(:, :, :)

    REAL(WP), ALLOCATABLE :: FJ(:, :)
    REAL(WP), ALLOCATABLE :: BCJ(:, :)
    REAL(WP), ALLOCATABLE :: F(:, :, :)

END MODULE

!**********************************************************************************************************************************
SUBROUTINE FFT99_POIS3D_INIT(iDomain_tmp)
    USE FFT99_info1
    USE FFT99_info2
    USE mesh_info
    USE mpi_info
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: iDomain_tmp

    IF(iDomain_tmp == ITG) L = NCL1_TG
    IF(iDomain_tmp == IIO) L = NCL1_io
    M = NCL2
    N = NCL3

    ML = N2DO(MYID)
    NL = N3DO(MYID)

    MLmax = N2DO(0)
    NLmax = N3DO(0)

    LHFP = L / 2 + 1

    DXQI0 = DXQI
    DZQI0 = DZQI

    ALLOCATE ( AK1(L) ) ; AK1 = 0.0_WP
    ALLOCATE ( AK3(N) ) ; AK3 = 0.0_WP

    ALLOCATE ( TRIGXX1(3 * L / 2 + 1) ) ; TRIGXX1 = 0.0_WP
    ALLOCATE ( TRIGXX2(3 * M / 2 + 1) ) ; TRIGXX2 = 0.0_WP
    ALLOCATE ( TRIGXX3(3 * N / 2 + 1) ) ; TRIGXX3 = 0.0_WP

    ALLOCATE ( TRIGXC1(2 * L) ) ; TRIGXC1 = 0.0_WP
    ALLOCATE ( TRIGXC2(2 * M) ) ; TRIGXC2 = 0.0_WP
    ALLOCATE ( TRIGXC3(2 * N) ) ; TRIGXC3 = 0.0_WP

    ALLOCATE (  XR  (L + 2, N)  )
    ALLOCATE (  WORK(L + 2, N)  )
    ALLOCATE (  XA  (N, L)  )
    ALLOCATE (  WOR (N, L)  )

    ALLOCATE ( AMJP(LHFP, M, NLmax ) ) ; AMJP = 0.0_WP
    ALLOCATE ( ACJP(LHFP, M, NLmax ) ) ; ACJP = 0.0_WP
    ALLOCATE ( APJP(LHFP, M, NLmax ) ) ; APJP = 0.0_WP

    ALLOCATE ( RHSLLPHIRe(LHFP, MLmax, N) ) ; RHSLLPHIRe = 0.0_WP
    ALLOCATE ( RHSLLPHIIm(LHFP, MLmax, N) ) ; RHSLLPHIIm = 0.0_WP

    ALLOCATE ( FJ (LHFP, M) ) ; FJ  = 0.0_WP
    ALLOCATE ( BCJ(LHFP, 2) ); BCJ = 0.0_WP

    ALLOCATE ( F (LHFP, NCL2, N3DO(0) ) )       ;  F = 0.0_WP

    MEMPC_Byte = MEMPC_Byte  + (L+ N +3 * L / 2 + 1 +3 * M / 2 + 1 +3 * N / 2 + 1 + 2 * L + 2 * M + 2 * N + &
    (2 * L+3) * N + N * L+5* LHFp * M * NLmaX +LHFp * (M + 2)) * 8




    IF(MYID == 0) &
    CALL FFT99_ROOT
    CALL BCAST_FFT99_ROOT

    CALL TDMA_PHI_COEF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE FFT99_ROOT
    USE FFT99_info1
    USE FFT99_info2
    USE mpi_info
    IMPLICIT NONE

    INTEGER(4) :: I, K
    INTEGER(4) :: N1MH, N3MH
    INTEGER(4) :: N1MP, N3MP

    REAL(WP), ALLOCATABLE :: AP(:)
    REAL(WP), ALLOCATABLE :: AN(:)
    REAL(WP) :: PI

    IF(MYID /= 0) RETURN

    !>      @note: set up coefficients IFXX1, TRIGXX1 in x and z directions
    CALL FFTFAX(L, IFXX1, TRIGXX1)
    CALL FFTFAX(N, IFXX3, TRIGXX3)
    !>       @note: set up coefficients IFXC1, TRIGXC1 in x and z directions
    CALL CFTFAX(N, IFXC3, TRIGXC3)
    CALL CFTFAX(L, IFXC1, TRIGXC1)


    N1MH = L / 2
    N3MH = N / 2
    N1MP = N1MH + 1
    N3MP = N3MH + 1

    PI = 2.0_WP * (DASIN(1.0_WP))

    !=============calculate (2 - 2cos(2PI / Nk)) / Dz^2 for z direction =============.
    ALLOCATE( AN(N) )
    AN = 0.0_WP
    DO K = 1, N3MH
        AN(K) = DBLE(K - 1) * 2.0_WP * PI
    END DO
    DO K = N3MP, N
        AN(K) = -DBLE(N - K + 1) * 2.0_WP * PI
    END DO

    DO K = 1, N
        AK3(K) = 2.0_WP * (1.0_WP - DCOS(AN(K) / DBLE(N))) * DZQI0
    END DO
    DEALLOCATE(AN)

    !=============calculate (2 -2cos(2PI / Nk)) / DX^2 for x direction =============.
    ALLOCATE( AP(L) )
    AP = 0.0_WP
    DO I = 1, N1MH
        AP(I) = DBLE(I - 1) * 2.0_WP * PI
    END DO
    DO I = N1MP,L
        AP(I) = -DBLE(L - I + 1) * 2.0_WP * PI
    END DO

    DO I = 1, L
        AK1(I) = 2.0_WP * (1.0_WP - DCOS(AP(I) / DBLE(L))) * DXQI0
    END DO
    DEALLOCATE(AP)

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE TDMA_PHI_COEF
    USE mesh_info
    USE FFT99_info1
    USE FFT99_info2
    IMPLICIT NONE

    INTEGER(4) :: I, J, K, KK

    ! For TDMA of Phi
    DO K = 1, NL
        KK = KCL2G(K)
        !>           @note calcuate FJ(I, J), ACJP(I, J), APJP(I, J), AMJP(I, J).
        DO I = 1, LHFP
            DO J = 1, M
                !ACC = 1.0_WP / (ACPH(J) + (-AK1(I) * (1.0_WP / RCCI2(J)) - AK3(KK)))
                APJP(I, J, K) = APPH(J)!* ACC
                AMJP(I, J, K) = AMPH(J)!* ACC
                ACJP(I, J, K) = ACPH(J) + (-AK1(I) / RCCI2(J) - AK3(KK))
            END DO

            !IF (I == 1 .AND. KK == 1 .AND. iCase /= IBox3P) THEN
            IF (I == 1 .AND. KK == 1) THEN
                !WRITE(*, *) '111, ', AMJP(1, 1, K), ACJP(1, 1, K), APJP(1, 1, K)
                ACJP(1, 1, 1) = ACPH(1) + (-AK1(1) / RCCI2(1) - AK3(1))
                APJP(1, 1, 1) = 0.0_WP
                AMJP(1, 1, 1) = 0.0_WP
                !WRITE(*, *) '111, ', AMJP(1, 1, K), ACJP(1, 1, K), APJP(1, 1, K)
            ENDIF
        END DO


    END DO

    RETURN

END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE BCAST_FFT99_ROOT
    USE FFT99_info1
    USE FFT99_info2
    USE mesh_info
    USE mpi_info
    IMPLICIT NONE

    CALL MPI_BCAST( AK1, L, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( AK3, N, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( TRIGXX1, 3 * L / 2 + 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TRIGXX2, 3 * M / 2 + 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TRIGXX3, 3 * N / 2 + 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

    CALL MPI_BCAST( TRIGXC1, 2 * L, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TRIGXC2, 2 * M, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
    CALL MPI_BCAST( TRIGXC3, 2 * N, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )


    CALL MPI_BCAST( IFXX1, 13, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IFXX2, 13, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IFXX3, 13, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IFXC1, 13, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IFXC2, 13, MPI_INTEGER4, 0, ICOMM, IERROR )
    CALL MPI_BCAST( IFXC3, 13, MPI_INTEGER4, 0, ICOMM, IERROR )

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE FFT99_POIS3D_periodicxz(iDomain_tmp)
    USE mpi_info
    USE FFT99_info1
    USE FFT99_info2
    USE mesh_info
    USE flow_info
    IMPLICIT none
    INTEGER(4), Intent(in) :: iDomain_tmp
    INTEGER(4) :: I, J, K, KK
    REAL(WP) :: UN3M

    RHSLLPHIRE = 0.0_WP
    RHSLLPHIIM = 0.0_WP

    ! FFT for x direction, and then CFFT for z direction.
    UN3M = 1.0_WP / DBLE(N)
    DO J = 1, ML

        IF(iDomain_tmp == ITG) THEN
            DO K = 1, N
                XR(1,  K) = RHSLLPHI_tg(L, J, K)
                XR(L + 2, K) = RHSLLPHI_tg(1, J, K)
            END DO
            DO I = 1, L
                DO K = 1, N
                    XR (I + 1, K) = RHSLLPHI_tg(I, J, K) !2 to NCL1 + 1
                    WORK(I, K ) = 0.0_WP
                END DO
            END DO
        END IF

        IF(iDomain_tmp == IIO) THEN
            DO K = 1, N
                XR(1,  K) = RHSLLPHI_io(L, J, K)
                XR(L + 2, K) = RHSLLPHI_io(1, J, K)
            END DO
            DO I = 1, L
                DO K = 1, N
                    XR (I + 1, K) = RHSLLPHI_io(I, J, K) !2 to NCL1 + 1
                    WORK(I, K ) = 0.0_WP
                END DO
            END DO
        END IF

        !WORK(:, :) = 0.0_WP
        CALL FFT99(XR, WORK, TRIGXX1, IFXX1, 1, L + 2, L, N, -1)

        DO K = 1, N
            DO I = 1, LHFP
                XA(K, I) = DCMPLX(XR(2 * I - 1, K), XR(2 * I, K))
                XR(2 * I - 1, K) = 0.0_WP
                XR(2 * I,  K) = 0.0_WP
            END DO
        END DO

        !XR(:, :) = 0.0_WP
        CALL CFFT99(XA, WOR, TRIGXC3, IFXC3, 1, N, N, LHFP, -1)

        DO I = 1, LHFP
            DO K = 1, N
                RHSLLPHIRe(I, J, K) = DREAL(XA(K, I) * UN3M)
                RHSLLPHIIm(I, J, K) = DIMAG(XA(K, I) * UN3M)
            END DO
        END DO

    END DO

    !       Correction for b.c.
    !WRITE(*, *) 'fftcheck',RHSLLPHIRe(1, 1, 1), RHSLLPHIIm(1, 1, 1)
    !IF (MYID == 0 .AND. iCase /= IBox3P) THEN !
    IF (MYID == 0 ) THEN
        RHSLLPHIRe(1, 1, 1) = 0.0_WP
        RHSLLPHIIm(1, 1, 1) = 0.0_WP
    ENDIF

    !       TDMA for the real part of FFTxz
    FJ  = 0.0_WP
    BCJ = 0.0_WP
    F = 0.0_WP
    !CALL TRASP23L2G_PHIRe
    CALL TRASP23_Y2Z(LHFP, 1, MLmax, RHSLLPHIRe, F)



    IF(iCase == IBox3P) THEN
        DO K = 1, NL
            KK = KCL2G(K)
            DO I = 1, LHFP
                DO J = 1, M
                    FJ(I, J) = F(I, J, K)
                END DO
            END DO

            IF (KK == 1) THEN
                FJ(1, 1) = 0.0_WP
            ENDIF

            CALL TDMAIJJ_CYC(AMJP(:, :, K), ACJP(:, :, K), APJP(:, :, K), FJ(:, :), 1, M, 1, LHFP)
            DO I = 1, LHFP
                DO J = 1, M
                    F(I, J, K) = FJ(I, J)
                END DO
            END DO

        END DO

    ELSE
        DO K = 1, NL
            KK = KCL2G(K)
            DO I = 1, LHFP
                DO J = 1, M
                    FJ(I, J) = F(I, J, K)
                END DO
                BCJ(I, :) = 0.0_WP
            END DO

            IF (KK == 1) THEN
                FJ(1, 1) = 0.0_WP
            ENDIF

            CALL TDMAIJJ_nonCYC(AMJP(:, :, K), ACJP(:, :, K), APJP(:, :, K), &
            FJ, BCJ, 1, M, 1, LHFP)
            DO I = 1, LHFP
                DO J = 1, M
                    F(I, J, K) = FJ(I, J)
                END DO
            END DO

        END DO
    END IF
    !CALL TRASP23G2L_PHIRe
    CALL TRASP23_Z2Y(LHFP, 1, MLmax, RHSLLPHIRe, F)

    !        TDMA for the imaginARy pARt of FFTxz
    FJ  = 0.0_WP
    BCJ = 0.0_WP
    F = 0.0_WP
    !CALL TRASP23L2G_PHIIm
    CALL TRASP23_Y2Z(LHFP,  1, MLmax, RHSLLPHIIm, F)

    IF(iCase == IBox3P) THEN

        DO K = 1, NL
            KK = KCL2G(K)
            DO I = 1, LHFP
                DO J = 1, M
                    FJ(I, J) = F(I, J, K)
                END DO
            END DO

            IF (KK == 1) THEN
                FJ(1, 1) = 0.0_WP
            ENDIF

            CALL TDMAIJJ_CYC(AMJP(:, :, K), ACJP(:, :, K), APJP(:, :, K), FJ(:, :), 1, M, 1, LHFP)

            DO I = 1, LHFP
                DO J = 1, M
                    F(I, J, K) = FJ(I, J)
                END DO
            END DO

        END DO
    ELSE
        DO K = 1, NL
            KK = KCL2G(K)
            DO I = 1, LHFP
                DO J = 1, M
                    FJ(I, J) = F(I, J, K)
                END DO
                BCJ(I, :) = 0.0_WP
            END DO

            IF (KK == 1) THEN
                FJ(1, 1) = 0.0_WP
            ENDIF

            CALL TDMAIJJ_nonCYC(AMJP(:, :, K), ACJP(:, :, K), APJP(:, :, K), &
            FJ, BCJ, 1, M, 1, LHFP)

            DO I = 1, LHFP
                DO J = 1, M
                    F(I, J, K) = FJ(I, J)
                END DO
            END DO

        END DO
    END IF
    !CALL TRASP23G2L_PHIIm
    CALL TRASP23_Z2Y(LHFP, 1, MLmax, RHSLLPHIIm, F)

    !        backward CFFT for z direction, and THEN backward FFT for x direction.
    IF(iDomain_tmp == ITG) DPH_tg = 0.0_WP
    IF(iDomain_tmp == IIO) DPH_io = 0.0_WP
    DO J = 1, ML

        DO K = 1, N
            DO I = 1, LHFP
                XA(K, I) = DCMPLX(RHSLLPHIRe(I, J, K), RHSLLPHIIm(I, J, K))
            END DO
        END DO

        CALL CFFT99(XA, WOR, TRIGXC3, IFXC3, 1, N, N, LHFP, +1)

        DO I = 1, LHFP
            DO  k = 1, N
                XR(2 * I - 1, K) = DREAL(XA(K, I))
                XR(2 * I,  K) = DIMAG(XA(K, I))
            END DO
        END DO

        CALL FFT99(XR, WORK, TRIGXX1, IFXX1, 1, L + 2, L, N, +1)

        IF(iDomain_tmp == ITG) THEN
            DO I = 1, L
                DO K = 1, N
                    DPH_tg (I, J, K) = XR(I + 1, K)
                    WORK(I, K  ) = 0.0_WP
                END DO
            END DO
        END IF

        IF(iDomain_tmp == IIO) THEN
            DO I = 1, L
                DO K = 1, N
                    DPH_io(I, J, K) = XR(I + 1, K)
                    WORK(I, K  ) = 0.0_WP
                END DO
            END DO
        END IF

    END DO

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!!========================================================================================================================
!SUBROUTINE TRASP23L2G_PHIRe   !I23 = 1
!        USE FFT99_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        REAL(WP) :: SENBLK( LHFP, N2DO(0), N3DO(0) )
!        REAL(WP) :: RECBLK( LHFP, N2DO(0), N3DO(0) )
!        INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
!        INTEGER(4) :: NTOL
!        INTEGER(4) :: JL, JG
!        INTEGER(4) :: IG
!        INTEGER(4) :: KL, KG
!        INTEGER(4) :: IP
!        INTEGER(4) :: ISEND
!        INTEGER(4) :: IRECV
!        INTEGER(4) :: ITAG
!        INTEGER(4) :: NNTR

!        NNTR = LHFP

!        DO IP = 0, NPSLV
!            IF(MYID == IP) THEN

!                DO JL = 1, N2DO(MYID)
!                    JG=JCL2G(JL)
!                    DO KL = 1, N3DO(MYID)
!                        KG = KDSWT(MYID) - 1 + KL
!                        DO IG = 1, NNTR
!                            F(IG, JG, KL) = RHSLLPHIRe(IG, JL, KG)
!                        END DO
!                    END DO
!                END DO

!            ELSE

!                NTOL = NNTR * N2DO(0) * N3DO(0)
!                DO KL = 1, N3DO(IP)
!                   KG = KDSWT(IP) - 1 + KL
!                   DO IG = 1, NNTR
!                      DO JL = 1, N2DO(MYID)
!                         SENBLK(IG, JL, KL) = RHSLLPHIRe(IG, JL, KG)
!                      END DO
!                   END DO
!                END DO

!                ISEND = IP
!                IRECV = IP
!                ITAG = 0
!                CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
!                     IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
!                     ISEND, ITAG, ICOMM, TRAS_STS, IERROR)

!                DO JL = 1, N2DO(IP)
!                   JG=JDSWT(IP) + JL- 1
!                    DO IG = 1, NNTR
!                        DO KL = 1, N3DO(MYID)
!                            F(IG, JG, KL) = RECBLK(IG, JL, KL)
!                        END DO
!                    END DO
!                END DO

!            END IF
!        END DO

!        RETURN
!    END SUBROUTINE

!!>*********************************************************************
!SUBROUTINE TRASP23G2L_PHIRe  !I23 = -1
!        USE FFT99_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        REAL(WP) :: SENBLK(LHFP, N2DO(0), N3DO(0))
!        REAL(WP) :: RECBLK(LHFP, N2DO(0), N3DO(0))
!        INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
!        INTEGER(4) :: NTOL
!        INTEGER(4) :: JL, JG
!        INTEGER(4) :: IG
!        INTEGER(4) :: KL, KG
!        INTEGER(4) :: IP
!        INTEGER(4) :: ISEND
!        INTEGER(4) :: IRECV
!        INTEGER(4) :: ITAG
!        INTEGER(4) :: NNTR



!        NNTR = LHFP
!        DO IP = 0, NPSLV
!            IF(MYID == IP) THEN

!                DO KL = 1, N3DO(MYID)
!                    KG = KDSWT(MYID) - 1 + KL
!                    DO JL = 1, N2DO(MYID)
!                        JG=JDSWT(MYID) - 1 + JL
!                        DO IG = 1, NNTR
!                            RHSLLPHIRe(IG, JL, KG) = F(IG, JG, KL)
!                        END DO
!                    END DO
!                END DO

!            ELSE
!                NTOL = NNTR * N3DO(0) * N2DO(0)
!                DO JL = 1, N2DO(IP)
!                    JG=JDSWT(IP) - 1 + JL
!                    DO IG = 1, NNTR
!                        DO KL = 1, N3DO(MYID)
!                            SENBLK(IG, JL, KL) = F(IG, JG, KL)
!                        END DO
!                    END DO
!                END DO
!                ISEND = IP
!                IRECV = IP
!                ITAG = 0
!                CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
!                    IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
!                    ISEND, ITAG, ICOMM, TRAS_STS, IERROR)


!                DO KL = 1, N3DO(IP)
!                    KG = KDSWT(IP) + KL- 1
!                    DO IG = 1, NNTR
!                        DO JL = 1, N2DO(MYID)
!                            RHSLLPHIRe(IG, JL, KG) = RECBLK(IG, JL, KL)
!                        END DO
!                    END DO
!                END DO

!            END IF
!        END DO

!        RETURN
!    END SUBROUTINE
!!===================================================================
!SUBROUTINE TRASP23L2G_PHIIm   !I23 = 1
!        USE FFT99_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        REAL(WP) :: SENBLK( LHFP, N2DO(0), N3DO(0) )
!        REAL(WP) :: RECBLK( LHFP, N2DO(0), N3DO(0) )
!        INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
!        INTEGER(4) :: NTOL
!        INTEGER(4) :: JL, JG
!        INTEGER(4) :: IG
!        INTEGER(4) :: KL, KG
!        INTEGER(4) :: IP
!        INTEGER(4) :: ISEND
!        INTEGER(4) :: IRECV
!        INTEGER(4) :: ITAG
!        INTEGER(4) :: NNTR

!        NNTR = LHFP

!        DO IP = 0, NPSLV
!            IF(MYID == IP) THEN

!                DO JL = 1, N2DO(MYID)
!                    JG=JCL2G(JL)
!                    DO KL = 1, N3DO(MYID)
!                        KG = KDSWT(MYID) - 1 + KL
!                        DO IG = 1, NNTR
!                            F(IG, JG, KL) = RHSLLPHIIm(IG, JL, KG)
!                        END DO
!                    END DO
!                END DO

!            ELSE

!                NTOL = NNTR * N2DO(0) * N3DO(0)
!                DO KL = 1, N3DO(IP)
!                    KG = KDSWT(IP) - 1 + KL
!                    DO IG = 1, NNTR
!                        DO JL = 1, N2DO(MYID)
!                            SENBLK(IG, JL, KL) = RHSLLPHIIm(IG, JL, KG)
!                        END DO
!                    END DO
!                END DO

!                ISEND = IP
!                IRECV = IP
!                ITAG = 0
!                CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
!                    IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
!                    ISEND, ITAG, ICOMM, TRAS_STS, IERROR)

!                DO JL = 1, N2DO(IP)
!                    JG=JDSWT(IP) + JL- 1
!                    DO IG = 1, NNTR
!                        DO KL = 1, N3DO(MYID)
!                            F(IG, JG, KL) = RECBLK(IG, JL, KL)
!                        END DO
!                    END DO
!                END DO

!            END IF
!        END DO

!        RETURN
!    END SUBROUTINE

!!>*********************************************************************
!SUBROUTINE TRASP23G2L_PHIIm  !I23 = -1
!        USE FFT99_info
!        USE mesh_info
!        USE flow_info
!        IMPLICIT NONE

!        REAL(WP) :: SENBLK(LHFP, N2DO(0), N3DO(0))
!        REAL(WP) :: RECBLK(LHFP, N2DO(0), N3DO(0))
!        INTEGER(4) :: TRAS_STS(MPI_STATUS_SIZE)
!        INTEGER(4) :: NTOL
!        INTEGER(4) :: JL, JG
!        INTEGER(4) :: IG
!        INTEGER(4) :: KL, KG
!        INTEGER(4) :: IP
!        INTEGER(4) :: ISEND
!        INTEGER(4) :: IRECV
!        INTEGER(4) :: ITAG
!        INTEGER(4) :: NNTR



!        NNTR = LHFP
!        DO IP = 0, NPSLV
!            IF(MYID == IP) THEN

!                DO KL = 1, N3DO(MYID)
!                    KG = KDSWT(MYID) - 1 + KL
!                    DO JL = 1, N2DO(MYID)
!                        JG=JDSWT(MYID) - 1 + JL
!                        DO IG = 1, NNTR
!                            RHSLLPHIIm(IG, JL, KG) = F(IG, JG, KL)
!                        END DO
!                    END DO
!                END DO

!            ELSE
!                NTOL = NNTR * N3DO(0) * N2DO(0)
!                DO JL = 1, N2DO(IP)
!                    JG=JDSWT(IP) - 1 + JL
!                    DO IG = 1, NNTR
!                        DO KL = 1, N3DO(MYID)
!                            SENBLK(IG, JL, KL) = F(IG, JG, KL)
!                        END DO
!                    END DO
!                END DO
!                ISEND = IP
!                IRECV = IP
!                ITAG = 0
!                CALL MPI_SENDRECV(SENBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION, &
!                    IRECV, ITAg, rECBLK(1, 1, 1), NTOL, MPI_DOUBLE_PRECISION,   &
!                    ISEND, ITAG, ICOMM, TRAS_STS, IERROR)


!                DO KL = 1, N3DO(IP)
!                    KG = KDSWT(IP) + KL- 1
!                    DO IG = 1, NNTR
!                        DO JL = 1, N2DO(MYID)
!                            RHSLLPHIIm(IG, JL, KG) = RECBLK(IG, JL, KL)
!                        END DO
!                    END DO
!                END DO

!            END IF
!        END DO

!        RETURN
!    END SUBROUTINE
