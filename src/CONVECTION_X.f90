!**********************************************************************************************************************************
!> @brief
!>        Calculate the convection term in x direction.
!> @details
!> SUBROUTINE: CONVECTION_X_tg (in MYID = all)
!> SUBROUTINE: CONVECTION_X_io(in MYID = all)
!> @note
!> @todo
! REVISION HISTORY:
! 05/2010 - Initial Version (tg domain only), by Mehdi Seddighi
! 04/2014 - Added io domain, optimized the code structure in f90, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
SUBROUTINE CONVECTION_X_tg
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJP
    INTEGER(4) :: KC, KM, KP
    REAL(WP) :: COE1, COE2, COE3
    REAL(WP) :: H11, H12, H13

    QTMP_tg(:, :, :) = 0.0_WP

    COE1 = XND2CL * XND2CL * DXI
    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)
        COE2 = DYFI(JJ) * XND2CL * RCCI1(JJ)
        COE3 = XND2CL  * ZND2CL * DZI * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)

            DO IC = 1, NCL1_tg
                IP = IPV_tg(IC)
                IM = IMV_tg(IC)
                H11 = ( (Q_tg(IP, JC, KC, 1) + Q_tg(IC, JC, KC, 1)) * &
                (Q_tg(IP, JC, KC, 1) + Q_tg(IC, JC, KC, 1))     -   &
                (Q_tg(IM, JC, KC, 1) + Q_tg(IC, JC, KC, 1)) *     &
                (Q_tg(IM, JC, KC, 1) + Q_tg(IC, JC, KC, 1)) ) * COE1

                H12 = ( (Q_tg(IC, JP, KC, 2) + Q_tg(IM, JP, KC, 2)) * &
                (YCL2ND_WFF(JJP) * Q_tg(IC, JP, KC, 1) + &
                YCL2ND_WFB(JJP) * Q_tg(IC, JC, KC, 1) ) - &
                (Q_tg(IC, JC, KC, 2) + Q_tg(IM, JC, KC, 2)) * &
                (YCL2ND_WFF(JJ) * Q_tg(IC, JC, KC, 1) + &
                YCL2ND_WFB(JJ) * Q_tg(IC, JM, KC, 1))  ) * COE2

                H13 = ( (Q_tg(IC, JC, KP, 3) + Q_tg(IM, JC, KP, 3)) *  &
                (Q_tg(IC, JC, KP, 1) + Q_tg(IC, JC, KC, 1))    -   &
                (Q_tg(IC, JC, KC, 3) + Q_tg(IM, JC, KC, 3)) *    &
                (Q_tg(IC, JC, KC, 1) + Q_tg(IC, JC, KM, 1)) ) * COE3

                QTMP_tg(IC, JC, KC) = -(H11 + H12 + H13)     !H for momentum x direction.
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE CONVECTION_X_io
    !>    The last time \rho IS USEd. No combination of N + 1 and n step of \rho

    USE thermal_info
    USE mesh_info
    USE flow_info
    IMPLICIT NONE

    INTEGER(4) :: IC, IM, IP
    INTEGER(4) :: JC, JM, JP, JJ, JJP, NXI
    INTEGER(4) :: KC, KM, KP
    REAL(WP) :: H11, H12, H13
    REAL(WP) :: H11F
    REAL(WP) :: H11B
    REAL(WP) :: H12F
    REAL(WP) :: H12B
    REAL(WP) :: H13F
    REAL(WP) :: H13B
    REAL(WP) :: COE1, COE2, COE3, COE32

    QTMP_io(:, :, :) = 0.0_WP
    COE1 = DXI * XND2CL * XND2CL
    COE3 = DZI * XND2CL * ZND2CL

    NXI = 1
    IF(TgFlowFlg) NXI = 2

    DO JC = 1, N2DO(MYID)
        JM = JLMV(JC)
        JP = JLPV(JC)
        JJ = JCL2G(JC)
        JJP = JGPV(JJ)
        COE2 = DYFI(JJ) * XND2CL * RCCI1(JJ)
        COE32 = COE3 * RCCI2(JJ)
        DO KC = 1, NCL3
            KM = KMV(KC)
            KP = KPV(KC)

            DO IC = NXI, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)

                ! \frac{\partial {\rho u u}}{\partial x}_{i', J, K}
                !(I, J, K)
                H11F = ( G_io(IP, JC, KC, 1) + G_io(IC, JC, KC, 1) ) *  &
                ( Q_io(IP, JC, KC, 1) + Q_io(IC, JC, KC, 1) )
                !(I - 1, J, K)
                H11B = ( G_io(IC, JC, KC, 1) + G_io(IM, JC, KC, 1) ) *  &
                ( Q_io(IC, JC, KC, 1) + Q_io(IM, JC, KC, 1) )
                !(I', J, K)
                H11  = (H11F - H11B) * COE1

                ! \frac{\partial {\rho v u}}{\partial y}_{i', J, K}
                !(I', J'+ 1, K)
                H12F = ( G_io(IC, JP, KC, 2) +  G_io(IM, JP, KC, 2) ) *  &
                ( YCL2ND_WFF(JJP) * Q_io(IC, JP, KC, 1) +      &
                YCL2ND_WFB(JJP) * Q_io(IC, JC, KC, 1) )
                !(I', J',  K)
                H12B = ( G_io(IC, JC, KC, 2) +  G_io(IM, JC, KC, 2) ) *  &
                ( YCL2ND_WFF(JJ) * Q_io(IC, JC, KC, 1) +    &
                YCL2ND_WFB(JJ) * Q_io(IC, JM, KC, 1) )
                !(I', J,  K)
                H12  = ( H12F - H12B) * COE2

                ! \frac{\partial {\rho w u}}{\partial z}_{i', J, K}
                !(I', J, K'+ 1)
                H13F = ( G_io(IC, JC, KP, 3) + G_io(IM, JC, KP, 3) ) *  &
                ( Q_io(IC, JC, KP, 1) + Q_io(IC, JC, KC, 1) )
                !(I', J, K')
                H13B = ( G_io(IC, JC, KC, 3) + G_io(IM, JC, KC, 3) ) *  &
                ( Q_io(IC, JC, KC, 1) + Q_io(IC, JC, KM, 1) )
                !(I', J, K)
                H13 = ( H13F - H13B) * COE32

                QTMP_io(IC, JC, KC) = -(H11 + H12 + H13)     !H for momentum x direction.

                !IF(MYID == 0 .AND. JC == N2DO(0)) WRITE(*, '(A, 3I4.1, 5ES13.5)') 'convx', &
                !JC, KC, IC,H12F,G_io(IC, JP, KC, 2), G_io(IM, JP, KC, 2), Q_io(IC, JP, KC, 1), Q_io(IC, JC, KC, 1)

            END DO

        END DO

    END DO

    !CALL DEBUG_WRT_UWP_io ! test
    !CALL DEBUG_WRT_LOCAL(QTMP_io, 1, N2DO(MYID), 'conx') ! test

    RETURN
END SUBROUTINE
