!**********************************************************************************************************************************
!> @brief
!>        common variables
!> @details
!> module: WPRECISION
!> module: mpi_info
!> module: cparam
!> module: WRT_INFO
!> module: mesh_info
!> module: BC_INLOUT_info
!> module: init_info
!> module: flow_info
!> module: CHT_info
!> module: thermal_info
!> module: postprocess_info
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014- Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
MODULE WPRECISION
    INTEGER, PARAMETER :: WP = 8 !KIND(0.0D0) !WORKING PRECESION
END MODULE WPRECISION
!**********************************************************************************************************************************
MODULE mpi_info
    INCLUDE 'mpif.h'

    INTEGER(4) ::  MYID
    INTEGER(4) ::  IERROR
    INTEGER(4) ::  SIZE
    INTEGER(4) ::  ICOMM
    INTEGER(4) ::  STS(MPI_STATUS_SIZE)

    INTEGER(4) ::  NPSLV
    INTEGER(4) ::  NPTOT

END MODULE mpi_info

!**********************************************************************************************************************************
MODULE cparam
    USE WPRECISION
    USE mpi_info

    INTEGER(4), PARAMETER :: NFLOW = 1
    INTEGER(4), PARAMETER ::  NDV = 3
    !>      @pARam NDV IS three directions

    INTEGER(4), PARAMETER :: IniField_random = 0
    INTEGER(4), PARAMETER :: IniField_extrapolation = 1
    INTEGER(4), PARAMETER :: IniField_reStart = 2

    REAL(WP), PARAMETER :: REALMIN = 1.0E-14_WP
    REAL(WP), PARAMETER :: REALMAX = 1.0E+14_WP
    REAL(WP), PARAMETER :: Cmu0 = 0.09_WP
    REAL(WP), PARAMETER :: TK2C = 273.15_WP

    INTEGER(4), PARAMETER :: IINLET = -1
    INTEGER(4), PARAMETER :: IOULET = 1
    INTEGER(4), PARAMETER :: IMAIND = 0
    INTEGER(4), PARAMETER :: IALLDM = 2

    INTEGER(4), PARAMETER :: IBCWALL = 1
    INTEGER(4), PARAMETER :: IBCPERI = 3
    INTEGER(4), PARAMETER :: IBCSYMM = 2


    INTEGER(4), PARAMETER :: iBotWall = 1
    INTEGER(4), PARAMETER :: iTopWall = 2

    INTEGER(4), PARAMETER ::  T_Asymptotic_Average = 1
    INTEGER(4), PARAMETER ::  T_Summing_average = 2

    INTEGER(4), PARAMETER ::  VisImplicit = 1
    INTEGER(4), PARAMETER ::  VisExplicit = 0

    INTEGER(4), PARAMETER ::  search_table = 1
    INTEGER(4), PARAMETER ::  properties_functions = 2

    INTEGER(4), PARAMETER ::  flgxz  = 1
    INTEGER(4), PARAMETER ::  flgxzt = 2

    INTEGER(4), PARAMETER ::  monatomic_gas = 1
    INTEGER(4), PARAMETER ::  diatomic_gas  = 2
    INTEGER(4), PARAMETER ::  trivalence_gaS = 3

    INTEGER(4), PARAMETER :: powerlawflg = 1
    INTEGER(4), PARAMETER :: sutherlandflg = 2

    INTEGER(4) :: iCase
    INTEGER(4), PARAMETER :: ICHANL = 1
    INTEGER(4), PARAMETER :: IPIPEC = 2
    INTEGER(4), PARAMETER :: IANNUL = 3
    INTEGER(4), PARAMETER :: IBox3P = 4

    INTEGER(4) :: iDomain
    INTEGER(4), PARAMETER :: ITG = 1
    INTEGER(4), PARAMETER :: IIO = 2
    INTEGER(4), PARAMETER :: ITGIO = 3

    INTEGER(4) :: iCHT

    INTEGER(4) :: iThermoDynamics
    INTEGER(4), PARAMETER :: iScpWater = 1
    INTEGER(4), PARAMETER :: iScpCO2 = 2
    INTEGER(4), PARAMETER :: iLiquidSodium = 3
    INTEGER(4), PARAMETER :: iLiquidLead = 4
    INTEGER(4), PARAMETER :: iLiquidBismuth = 5
    INTEGER(4), PARAMETER :: iLiquidLBE = 6

    REAL(WP), PARAMETER :: Tm0_Na = 371.0 ! unit: K, melting temperature at 1 atm for Na
    REAL(WP), PARAMETER :: Tm0_Pb = 600.6 ! unit: K, melting temperature at 1 atm for Lead
    REAL(WP), PARAMETER :: Tm0_BI = 544.6 ! unit: K, melting temperature at 1 atm for Bismuth
    REAL(WP), PARAMETER :: Tm0_LBE = 398.0 ! unit: K, melting temperature at 1 atm for LBE

    REAL(WP), PARAMETER :: Tb0_Na = 1155.0 ! unit: K, boling temperature at 1 atm for Na
    REAL(WP), PARAMETER :: Tb0_Pb = 2021.0 ! unit: K, boling temperature at 1 atm for Lead
    REAL(WP), PARAMETER :: Tb0_BI = 1831.0 ! unit: K, boling temperature at 1 atm for Bismuth
    REAL(WP), PARAMETER :: Tb0_LBE = 1927.0 ! unit: K, boling temperature at 1 atm for LBE

    REAL(WP), PARAMETER :: Hm0_Na = 113.0e3 ! unit: J / Kg, latent melting heat, enthalpy
    REAL(WP), PARAMETER :: Hm0_Pb = 23.07e3 ! unit: J / Kg, latent melting heat, enthalpy
    REAL(WP), PARAMETER :: Hm0_BI = 53.3e3 ! unit: J / Kg, latent melting heat, enthalpy
    REAL(WP), PARAMETER :: Hm0_LBE = 38.6e3 ! unit: J / Kg, latent melting heat, enthalpy
    ! D = CoD(0) + CoD(1) * T
    REAL(WP), PARAMETER :: CoD_Na(0:1) = (/1014.0, -0.235/)
    REAL(WP), PARAMETER :: CoD_Pb(0:1) = (/11441.0, -1.2795/)
    REAL(WP), PARAMETER :: CoD_Bi(0:1) = (/10725.0, -1.22 /)
    REAL(WP), PARAMETER :: CoD_LBE(0:1) = (/11065.0, 1.293 /)
    ! K = CoK(0) + CoK(1) * T + CoK(2) * T^2
    REAL(WP), PARAMETER :: CoK_Na(0:2) = (/104.0, -0.047, 0.0/)
    REAL(WP), PARAMETER :: CoK_Pb(0:2) = (/9.2, 0.011, 0.0/)
    REAL(WP), PARAMETER :: CoK_Bi(0:2) = (/7.34, 9.5E-3, 0.0/)
    REAL(WP), PARAMETER :: CoK_LBE(0:2) = (/ 3.284, 1.617E-2, -2.305E-6/)
    ! B = 1 / (CoB - T)
    REAL(WP), PARAMETER :: CoB_Na = 4316.0
    REAL(WP), PARAMETER :: CoB_Pb = 8942.0
    REAL(WP), PARAMETER :: CoB_BI = 8791.0
    REAL(WP), PARAMETER :: CoB_LBE = 8558.0
    ! Cp = CoCp(-2) * T^(-2) + CoCp(-1) * T^(-1) + CoCp(0) + CoCp(1) * T + CoCp(2) * T^2
    REAL(WP), PARAMETER :: CoCp_Na(-2:2) = (/- 3.001e6, 0.0, 1658.0, -0.8479, 4.454E-4/)
    REAL(WP), PARAMETER :: CoCp_Pb(-2:2) = (/- 1.524e6, 0.0, 176.2, -4.923E-2, 1.544E-5/)
    REAL(WP), PARAMETER :: CoCp_Bi(-2:2) = (/7.183e6, 0.0, 118.2, 5.934E-3, 0.0/)
    REAL(WP), PARAMETER :: CoCp_LBE(-2:2) = (/-4.56e5, 0.0, 164.8, - 3.94E-2, 1.25E-5/)
    ! H = Hm0 + CoH(-1) * (1 / T - 1 / Tm0) + CoH(0) + CoH(1) * (T - Tm0) +  CoH(2) * (T^2 - Tm0^2) +  CoH(3) * (T^3- Tm0^3)
    REAL(WP), PARAMETER :: CoH_Na(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)
    REAL(WP), PARAMETER :: CoH_Pb(-1:3) = (/1.524e6, 0.0, 176.2, -2.4615E-2, 5.147E-6/)
    REAL(WP), PARAMETER :: CoH_Bi(-1:3) = (/-7.183e6, 0.0, 118.2, 2.967E-3, 0.0/)
    REAL(WP), PARAMETER :: CoH_LBE(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)! check, WRong from literature.
    ! M = vARies
    REAL(WP), PARAMETER :: CoM_Na(-1:1) = (/556.835, -6.4406, -0.3958/) ! M = exp ( CoM(-1) / T + CoM(0) + CoM(1) * ln(T) )
    REAL(WP), PARAMETER :: CoM_Pb(-1:1) = (/1069.0, 4.55E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
    REAL(WP), PARAMETER :: CoM_Bi(-1:1) = (/780.0, 4.456E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
    REAL(WP), PARAMETER :: CoM_LBE(-1:1) = (/754.1, 4.94E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)

    INTEGER(4) ::  MGRID, JINI
    INTEGER(4), ALLOCATABLE ::  JGMOV(:)

    INTEGER(4) ::  NCL1_tg
    INTEGER(4) ::  NCL1_io
    INTEGER(4) ::  NCL2
    INTEGER(4) ::  NCL3

    INTEGER(4) :: NCL1S
    INTEGER(4) :: NCL1E

    INTEGER(4) ::  NND1_tg
    INTEGER(4) ::  NND1_io
    INTEGER(4) ::  NND2
    INTEGER(4) ::  NND3

    INTEGER(4) ::  Kronecker_Delta(NDV, NDV)


    LOGICAL ::  TgFlowFlg, IoFlowFlg

    INTEGER(4) :: flag_random_real = 1
    INTEGER(4) :: flag_random_fixed = 2
    INTEGER(4) :: flag_random_sine = 3
    INTEGER(4) :: iRandomType

    INTEGER(4) :: MEMPC_Byte = 0  ! memory per core in the unit of byte

END MODULE cparam

!**********************************************************************************************************************************
MODULE WRT_INFO
    USE cparam
    CHARACTER(8) :: date
    CHARACTER(10) :: time

    CHARACTER(18) :: fllog
    INTEGER(4) :: logflg_pg = 200
    INTEGER(4) :: logflg_tg = 201
    INTEGER(4) :: logflg_io = 202


    CHARACTER(LEN = 256) :: WRT_RST_FNM_TG(NDV + 1)
    CHARACTER(LEN = 256) :: WRT_RST_FNM_io(NDV + 1)
    CHARACTER(LEN = 256) :: WRT_AVE_FNM_TG
    CHARACTER(LEN = 256) :: WRT_AVE_FNM_io

    CHARACTER(LEN = 64) :: DIR0 = '0_log_monitors'
    CHARACTER(LEN = 64) :: DIR1 = '1_instant_rawdata'
    CHARACTER(LEN = 64) :: DIR2 = '2_averaged_rawdata'
    CHARACTER(LEN = 64) :: DIR3 = '3_interpolation_rawdata'
    CHARACTER(LEN = 64) :: DIR4 = '4_averaged_pltdata'
    CHARACTER(LEN = 64) :: DIR5 = '5_instant_pltdata'
    CHARACTER(LEN = 64) :: DIR6 = '6_src_Backup'

    CHARACTER(LEN = 256) :: FilePath0 = '0_log_monitors/'
    CHARACTER(LEN = 256) :: FilePath1 = '1_instant_rawdata/'
    CHARACTER(LEN = 256) :: FilePath2 = '2_averaged_rawdata/'
    CHARACTER(LEN = 256) :: FilePath3 = '3_interpolation_rawdata/'
    CHARACTER(LEN = 256) :: FilePath4 = '4_averaged_pltdata/'
    CHARACTER(LEN = 256) :: FilePath5 = '5_instant_pltdata/'


END MODULE WRT_INFO

!**********************************************************************************************************************************
MODULE mesh_info
    USE WRT_INFO
    USE cparam

    INTEGER(4), ALLOCATABLE :: N2DO(:)
    INTEGER(4), ALLOCATABLE :: N3DO(:)

    INTEGER(4), ALLOCATABLE :: JDSWT(:)
    INTEGER(4), ALLOCATABLE :: JDEWT(:)
    INTEGER(4), ALLOCATABLE :: KDSWT(:)
    INTEGER(4), ALLOCATABLE :: KDEWT(:)
    INTEGER(4), ALLOCATABLE :: JCL2G(:)
    INTEGER(4), ALLOCATABLE :: KCL2G(:)
    INTEGER(4), ALLOCATABLE :: JG2IP(:)
    INTEGER(4), ALLOCATABLE :: JG2LC(:)

    INTEGER(4), ALLOCATABLE :: IPV_tg(:)
    INTEGER(4), ALLOCATABLE :: IMV_tg(:)
    INTEGER(4), ALLOCATABLE :: IPV_io(:)
    INTEGER(4), ALLOCATABLE :: IMV_io(:)
    INTEGER(4), ALLOCATABLE :: KPV(:)
    INTEGER(4), ALLOCATABLE :: KMV(:)
    INTEGER(4), ALLOCATABLE :: JLPV(:)
    INTEGER(4), ALLOCATABLE :: JLMV(:)
    INTEGER(4), ALLOCATABLE :: JGPV(:)
    INTEGER(4), ALLOCATABLE :: JGMV(:)
    INTEGER(4), ALLOCATABLE :: KSYM(:)


    REAL(WP) :: HX_tg, HX_io, HZ, HYB, HYT
    REAL(WP) :: ALX1(2), ALX2, ALX3
    REAL(WP) :: DX, DXI, DXQI
    REAL(WP) :: DZ, DZI, DZQI
    REAL(WP) :: VL1313_tg, VL1313_io

    REAL(WP), ALLOCATABLE :: CFLVIS(:)

    REAL(WP), ALLOCATABLE :: AMPH(:)
    REAL(WP), ALLOCATABLE :: ACPH(:)
    REAL(WP), ALLOCATABLE :: APPH(:)

    REAL(WP) :: AMPH0
    REAL(WP) :: APPH0
    !REAL(WP), ALLOCATABLE :: DPDYWAL(:, :, :)

    REAL(WP) :: AMVR1, APVRN

    REAL(WP), ALLOCATABLE :: AMVR(:, :)
    REAL(WP), ALLOCATABLE :: ACVR(:, :)
    REAL(WP), ALLOCATABLE :: APVR(:, :)

    REAL(WP), ALLOCATABLE :: XND_tg(:)
    REAL(WP), ALLOCATABLE :: XND_io(:)
    REAL(WP), ALLOCATABLE :: YND(:)
    REAL(WP), ALLOCATABLE :: ZND(:)
    REAL(WP), ALLOCATABLE :: ZCC(:)

    REAL(WP), ALLOCATABLE :: RCCI1(:)
    REAL(WP), ALLOCATABLE :: RNDI1(:)
    REAL(WP), ALLOCATABLE :: RCCI2(:)
    REAL(WP), ALLOCATABLE :: RNDI2(:)

    REAL(WP), ALLOCATABLE :: YCC(:)

    REAL(WP), ALLOCATABLE :: DYFI(:)
    REAL(WP), ALLOCATABLE :: DYCI(:)

    !===============WEIGHTING FACTORS FOR X, Y, Z AVERAGE (NODE VS CELL CENTRE) ==
    REAL(WP) :: XND2CL
    REAL(WP) :: YND2CL
    REAL(WP) :: ZND2CL

    REAL(WP), ALLOCATABLE :: Xcc_tg(:)
    REAL(WP), ALLOCATABLE :: Xcc_io(:)

    REAL(WP), ALLOCATABLE :: YCL2ND_WFF(:)
    REAL(WP), ALLOCATABLE :: YCL2ND_WFB(:)

END MODULE mesh_info

!**********************************************************************************************************************************
MODULE BC_INLOUT_info
    USE cparam

    REAL(WP) :: U_OUTLET
    REAL(WP), ALLOCATABLE :: BC_CONV0(:, :, :)
    REAL(WP), ALLOCATABLE :: BC_CONV0_ENG(:, :)
    REAL(WP), ALLOCATABLE :: BC_TDMA (:, :, :, :)
    REAL(WP), ALLOCATABLE :: BC_U_SSTAR(:, :, :)

END MODULE BC_INLOUT_info

!**********************************************************************************************************************************
MODULE init_info
    USE WRT_INFO
    USE cparam

    INTEGER(4) :: ITERG0_io
    INTEGER(4) :: ITERG0_TG
    INTEGER(4) :: ITERG0

    INTEGER(4) :: ITERG

    REAL(WP) :: PhyTIME_io
    REAL(WP) :: PhyTIME_TG
    REAL(WP) :: PhyTIME


    INTEGER(4) :: NTSTF
    INTEGER(4) :: NSST
    INTEGER(4) :: iIniField_tg, iIniField_io, iIniFieldType, iIniFieldTime
    INTEGER(4) :: iPostProcess, iPPInst, iPPSpectra, iPPDimension, pp_instn_sz, iPPQuadrants
    REAL(WP), ALLOCATABLE :: pp_instn_tim(:)

    !INTEGER(4) :: ISTR2
    CHARACTER(64) :: zoneNameView

    INTEGER(4) :: iFlowDriven

    INTEGER(4) ::  iWeightedPre
    INTEGER(4) ::  iVisScheme
    INTEGER(4) ::  iThermoProperty
    INTEGER(4) ::  iFluidMedia

    REAL(WP) :: PI
    REAL(WP) :: DT     ! real time step
    REAL(WP) :: DT0    ! given time step
    REAL(WP) :: DTMIN  ! the minimum time step

    REAL(WP) :: CPUTIME, CPUTIME_tmp
    REAL(WP) :: CFLGV


    REAL(WP) :: TimeReStart_tg
    REAL(WP) :: TimeReStart_io
    REAL(WP) :: TSTOP
    REAL(WP) :: timeThermoStart
    REAL(WP) :: timeFlowStart
    REAL(WP) :: TLgRe

    REAL(WP) :: dtSave
    REAL(WP) :: dtSave1
    REAL(WP) :: dtAveView
    REAL(WP) :: dtRawView
    INTEGER(4) :: iterMonitor

    REAL(WP) :: tRunAve1, tRunAve_Reset

    REAL(WP) :: Cf_Given

    REAL(WP) :: ReIni
    REAL(WP) :: REN
    REAL(WP) :: PRT0

    REAL(WP) :: Area_inlet, VOLM_tg, VOLM_io
    REAL(WP) :: CVISC

    REAL(WP) :: VPERG
    REAL(WP) :: SVPERG

    REAL(WP) :: STR2
    REAL(WP) :: TGAM(0:3), TROH(0:3), TALP(0:3)

    REAL(WP), ALLOCATABLE :: Vini(:)


    INTEGER(4) :: BCX_tg(2)
    INTEGER(4) :: BCZ(2)
    INTEGER(4) :: BCY(2)

    INTEGER(4) :: BCX_io(2)



END MODULE init_info

!**********************************************************************************************************************************
MODULE flow_info
    USE cparam
    USE BC_INLOUT_info

    REAL(WP) :: FACOE_TG, FACOE_io
    REAL(WP) :: U1mean_tg,  U1maxx_tg, U1mean_WORK_tg,  U1maxx_WORK_tg

    REAL(WP) :: G1maxx_io, G1maxx_WORK_io
    REAL(WP) :: G1rate_io, G1rate_WORK_io, G1BULK_WORK_io
    REAL(WP) :: T1maxx_io, T1maxx_WORK_io
    REAL(WP) :: DH1rate_io, DH1rate_WORK_io, DH1BULK_WORK_io, T1BULK_WORK_io


    REAL(WP) :: CHK_Mass_CONSV0, CHK_Mass_CONSV0_SCALING
    REAL(WP) :: CHK_ENEG_CONSV0, CHK_ENEG_TOTAL

    !=============== TURBULENCE ==GENERATOR =================================
    REAL(WP) :: MAXDIVGV_tg(2)
    REAL(WP) :: CFLMM_tg
    REAL(WP) :: CFLMM
    REAL(WP) :: VMAX_tg(3), VMIN_tg(3)

    !=========PRIMARY == variables == DIMENSIONLESS ==========
    REAL(WP), ALLOCATABLE :: Q_tg   (:, :, :, :)
    REAL(WP), ALLOCATABLE :: PR_tg  (:, :, :)

    !========= variables == iN ==CALCULATION =================
    REAL(WP), ALLOCATABLE :: RHS_tg     (:, :, :)
    REAL(WP), ALLOCATABLE :: CONVH0_tg  (:, :, :, :)
    REAL(WP), ALLOCATABLE :: DPH_tg     (:, :, :)
    REAL(WP), ALLOCATABLE :: RHSLLPHI_tg(:, :, :)

    !======== ASSISTANT == variables == iN ==CALCULATION =======
    !REAL(WP), ALLOCATABLE :: F   (:, :, :)
    REAL(WP), ALLOCATABLE :: QTMP_tg(:, :, :)

    !=========== main == DOMAIN =================================================
    REAL(WP) :: MAXDIVGV_io(3) ! 1formAInDOmAIn, 2formAInDOmAIn, 3forb.c.outlet.
    REAL(WP) :: VMAX_io(3), VMIN_io(3)
    REAL(WP) :: CFLMM_io

    !=========PRIMARY == variables == DIMENSIONLESS ==========
    REAL(WP), ALLOCATABLE :: Q_io   (:, :, :, :)
    REAL(WP), ALLOCATABLE :: G_io   (:, :, :, :)
    REAL(WP), ALLOCATABLE :: G0_io   (:, :, :, :)
    REAL(WP), ALLOCATABLE :: PR_io(:, :, :)
    REAL(WP), ALLOCATABLE :: PR0_io(:, :, :, :)
    REAL(WP), PARAMETER :: pres_epslon = 0.001_WP

    !========= variables == iN ==CALCULATION =================
    REAL(WP), ALLOCATABLE :: RHS_io     (:, :, :)

    !REAL(WP), ALLOCATABLE :: viscs0_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: EXPLT0_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DPH_io     (:, :, :)
    REAL(WP), ALLOCATABLE :: RHSLLPHI_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DivU_io(:, :, :)


    !======== ASSISTANT == variables == iN ==CALCULATION =======
    REAL(WP), ALLOCATABLE :: QTMP_io(:, :, :)

    REAL(WP), ALLOCATABLE :: UU(:, :, :)

END MODULE flow_info

!**********************************************************************************************************************************
MODULE  CHT_info
    USE cparam

    INTEGER(4) :: NCL2_Solid(2)
    INTEGER(4) :: NCL1s_CHT(2)
    INTEGER(4) :: NCL1e_CHT(2)

    REAL(WP) :: HY_Solid(2)
    REAL(WP) :: HXs_Solid(2)
    REAL(WP) :: HXe_Solid(2)
    REAL(WP) :: Cp_Solid(2)
    REAL(WP) :: D_Solid(2)
    REAL(WP) :: K_Solid(2)
    REAL(WP) :: HEAT_SRC_SOLI(2)
END MODULE CHT_info

!**********************************************************************************************************************************
MODULE thermal_info
    USE cparam
    !============== iNI GIVEN ================
    INTEGER(4) :: FLOWDIR
    INTEGER(4) :: iGravity
    CHARACTER(64) :: NISTFLNM
    INTEGER(4), PARAMETER :: BC_Fixed_Heat_Flux = 1
    INTEGER(4), PARAMETER :: BC_Fixed_Temperature = 2
    INTEGER(4), PARAMETER :: IsoPowerDENSITY = 3
    INTEGER(4) :: iThermalWallType(2)

    REAL(WP) :: CTHECD



    REAL(WP) :: G_A
    REAL(WP) :: F_A
    INTEGER(4) :: IBuoF(3)

    REAL(WP) :: Gr0(2)
    REAL(WP) :: BO0(2, 2)

    !===========for perfect gaS ==============
    INTEGER(4) :: IdealGasType
    INTEGER(4) :: fstatetype ! flow state equation.
    REAL(WP) :: R0
    REAL(WP) :: CV0
    REAL(WP) :: gamma0
    REAL(WP) :: Cvhat0
    REAL(WP) :: RHOU20

    REAL(WP) :: PL_MT  = 0.67_WP
    REAL(WP) :: PL_CpT = 0.095_WP
    REAL(WP) :: PL_HT  = 1.095_WP
    REAL(WP) :: PL_KT  = 0.805_WP

    !============== iNI FOR REFEENCE STATE == DIMENSIONAL ====
    REAL(WP) :: L0  !  (m)
    REAL(WP) :: U0  !  (M /s)
    REAL(WP) :: P0  !  (Pa)
    REAL(WP) :: T0  !  (K)
    REAL(WP) :: Ti  !  (K)   updated by Junjie, 2017/03 / 26
    REAL(WP) :: H0  !  (Pa)
    REAL(WP) :: D0  !  (J / Kg)
    REAL(WP) :: M0  !  (Pa-s)
    REAL(WP) :: K0  !  (W/ M - K)
    REAL(WP) :: CP0 ! (J / Kg- K)
    REAL(WP) :: B0  !  (1 / K)
    REAL(WP) :: thermalWallBC_Dim(2) ! (W/ M2) given heat flux rate on the wall or ! (K) given constant wall temperature

    REAL(WP) :: thermalWall_nondim(2)

    !REAL(WP) :: WHFLUX_UND(2)
    REAL(WP) :: DHmax, T4MaxDH
    REAL(WP) :: DHmin, T4MinDH
    REAL(WP) :: CpMax, T4CpMax
    INTEGER(4) :: IMAX_DH, IMIN_DH
    INTEGER(4) :: IMAX_H,  IMIN_H
    INTEGER(4) :: IMAX_T,  IMIN_T
    INTEGER(4) :: IMAX_Cp, IMIN_Cp

    !============== iNI FOR INLETSTATE ==dimensionless ====
    REAL(WP) :: T_inlet    !updated by Junjie, 2017/03 / 26
    REAL(WP) :: H_inlet    !updated by Junjie, 2017/03 / 26
    REAL(WP) :: D_inlet = 1.0  !updated by Junjie, 2017/03 / 26
    REAL(WP) :: DH_inlet !updated by Junjie, 2017/03 / 26
    REAL(WP) :: M_inlet = 1.0    !updated by Junjie, 2017/03 / 26
    REAL(WP) :: K_inlet    !updated by Junjie, 2017/03 / 26
    REAL(WP) :: CP_inlet

    !=============== variables ==== DIMENSIONALESS ===============
    REAL(WP), ALLOCATABLE :: DH(:, :, :)
    REAL(WP), ALLOCATABLE :: DENSITY(:, :, :)
    REAL(WP), ALLOCATABLE :: TEMPERATURE(:, :, :)
    REAL(WP), ALLOCATABLE :: Viscousity(:, :, :)
    REAL(WP), ALLOCATABLE :: Viscousity0(:, :, :)
    REAL(WP), ALLOCATABLE :: THERMCONDT(:, :, :)
    REAL(WP), ALLOCATABLE :: HEATCAP(:, :, :)
    REAL(WP), ALLOCATABLE :: ENTHALPY(:, :, :)
    REAL(WP), ALLOCATABLE :: D_STG(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DRHOI_STG(:, :, :, :)
    REAL(WP), ALLOCATABLE :: MU_STG(:, :, :, :)

    !============== variables ==== iN =CALCULATION ================
    REAL(WP), ALLOCATABLE :: RHS_ENERGY (:, :, :)
    REAL(WP), ALLOCATABLE :: RHS_ENERGY0(:, :, :)
    REAL(WP), ALLOCATABLE :: DENSITY0   (:, :, :)
    REAL(WP), ALLOCATABLE :: DrhoDtP    (:, :, :)
    REAL(WP), ALLOCATABLE :: DH0   (:, :, :)
    REAL(WP), ALLOCATABLE :: WALLFLUX   (:, :, :)

    !REAL(WP), ALLOCATABLE :: DENSITY1   (:, :, :)

    !============== NIST TABLE INFO=====================
    INTEGER(4) :: N_LIST
    REAL(WP), ALLOCATABLE :: LIST_H(:)
    REAL(WP), ALLOCATABLE :: LIST_T(:)
    REAL(WP), ALLOCATABLE :: LIST_D(:)
    REAL(WP), ALLOCATABLE :: LIST_M(:)
    REAL(WP), ALLOCATABLE :: LIST_K(:)
    REAL(WP), ALLOCATABLE :: LIST_CP(:)
    REAL(WP), ALLOCATABLE :: LIST_B(:)
    REAL(WP), ALLOCATABLE :: LIST_DH(:)!updated by Junjie, 2017/03 /13


    REAL(WP), ALLOCATABLE :: SplineCoeff_HT_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HD_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HM_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HK_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HCp_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HB_B(:)

    REAL(WP), ALLOCATABLE :: SplineCoeff_DHH_B(:)!updated by Junjie, 2017/03 / 26
    REAL(WP), ALLOCATABLE :: SplineCoeff_DHH_C(:)!updated by Junjie, 2017/03 / 26
    REAL(WP), ALLOCATABLE :: SplineCoeff_DHH_D(:)!updated by Junjie, 2017/03 / 26

    REAL(WP), ALLOCATABLE :: SplineCoeff_DHT_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_DHT_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_DHT_D(:)

    REAL(WP), ALLOCATABLE :: SplineCoeff_HT_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HD_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HM_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HK_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HCp_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HB_C(:)


    REAL(WP), ALLOCATABLE :: SplineCoeff_HT_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HD_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HM_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HK_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HCp_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_HB_D(:)

    REAL(WP), ALLOCATABLE :: SplineCoeff_TH_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TD_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TM_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TK_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TCp_B(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TB_B(:)

    REAL(WP), ALLOCATABLE :: SplineCoeff_TH_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TD_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TM_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TK_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TCp_C(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TB_C(:)

    REAL(WP), ALLOCATABLE :: SplineCoeff_TH_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TD_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TM_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TK_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TCp_D(:)
    REAL(WP), ALLOCATABLE :: SplineCoeff_TB_D(:)


    !============== SOME CONSTANT ==========================
    INTEGER(4) :: NTHERMAL = 7

    !========================
    REAL(WP), ALLOCATABLE :: T_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: H_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: D_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: M_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: K_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: Cp_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: DH_WAL_GV(:, :)
    REAL(WP), ALLOCATABLE :: B_WAL_GV(:, :)

END MODULE

!**********************************************************************************************************************************


!    MODULE WALLSTRESS_info
!        USE cparam

!        REAL(WP), ALLOCATABLE :: Cf_LW_io(:)
!        REAL(WP), ALLOCATABLE :: Cf_UW_io(:)
!        REAL(WP), ALLOCATABLE :: Cf_LW_tg(:)
!        REAL(WP), ALLOCATABLE :: Cf_UW_tg(:)
!        INTEGER(4) :: NCOUNT1 = 0
!        INTEGER(4) :: TECFLG
!    END MODULE WALLSTRESS_info

!**********************************************************************************************************************************
MODULE postprocess_info
    USE cparam

    INTEGER(4) :: NSTATIS_tg
    INTEGER(4) :: NSTATIS_io

    !=============== TG pART ==== Space averaged===========================
    REAL(WP), ALLOCATABLE :: U1xzL_tg(:, :), U1xztL_tg(:, :) !Ui
    REAL(WP), ALLOCATABLE :: UPxzL_tg(:, :), UPxztL_tg(:, :) !UiP
    REAL(WP), ALLOCATABLE :: U2xzL_tg(:, :), U2xztL_tg(:, :) !UiUj
    REAL(WP), ALLOCATABLE :: U3xzL_tg(:, :), U3xztL_tg(:, :) !UiUjUk

    REAL(WP), ALLOCATABLE :: DVDL1xzL_tg(:, :, :), DVDL1xztL_tg(:, :, :) ! dUI / DXj
    REAL(WP), ALLOCATABLE :: DVDLPxzL_tg(:, :, :), DVDLPxztL_tg(:, :, :) !PdUI / DXj
    REAL(WP), ALLOCATABLE :: DVDL2xzL_tg(:, :, :), DVDL2xztL_tg(:, :, :) ! dUI / DXj * dUI / DXj

    !=============== iO pARt (Non -periodic) ==== Space averaged================
    REAL(WP), ALLOCATABLE :: U1zL_io(:, :, :),  U1ztL_io(:, :, :)!Ui
    REAL(WP), ALLOCATABLE :: G1zL_io(:, :, :),  G1ztL_io(:, :, :)!Gi
    REAL(WP), ALLOCATABLE :: UPzL_io(:, :, :),  UPztL_io(:, :, :)!UiP

    REAL(WP), ALLOCATABLE :: U2zL_io(:, :, :),  U2ztL_io(:, :, :)!UiUj
    REAL(WP), ALLOCATABLE :: UGzL_io(:, :, :),  UGztL_io(:, :, :)!UiGj

    REAL(WP), ALLOCATABLE :: U3zL_io(:, :, :),  U3ztL_io(:, :, :)!UiGjUk
    REAL(WP), ALLOCATABLE :: UGUzL_io(:, :, :), UGUztL_io(:, :, :)!UiGjUk

    REAL(WP), ALLOCATABLE :: DVDL1zL_io(:, :, :, :), DVDL1ztL_io(:, :, :, :) ! dUI / DXj
    REAL(WP), ALLOCATABLE :: DVDLPzL_io(:, :, :, :), DVDLPztL_io(:, :, :, :) !PdUI / DXj
    REAL(WP), ALLOCATABLE :: DVDL2zL_io(:, :, :, :), DVDL2ztL_io(:, :, :, :) ! dUI / DXj * dUI / DXj

    REAL(WP), ALLOCATABLE :: T1zL_io(:, :), T1ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: D1zL_io(:, :), D1ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: H1zL_io(:, :), H1ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: K1zL_io(:, :), K1ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: M1zL_io(:, :), M1ztL_io(:, :)

    REAL(WP), ALLOCATABLE :: T2zL_io(:, :), T2ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: D2zL_io(:, :), D2ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: H2zL_io(:, :), H2ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: K2zL_io(:, :), K2ztL_io(:, :)
    REAL(WP), ALLOCATABLE :: M2zL_io(:, :), M2ztL_io(:, :)

    REAL(WP), ALLOCATABLE :: DHzL_io(:, :), DHztL_io(:, :)
    REAL(WP), ALLOCATABLE :: PHzL_io(:, :), PHztL_io(:, :)

    REAL(WP), ALLOCATABLE :: DVDL1MzL_io(:, :, :, :),    DVDL1MztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DVDL1MHzL_io(:, :, :, :),   DVDL1MHztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DVDL1MUzL_io(:, :, :, :, :), DVDL1MUztL_io(:, :, :, :, :)
    REAL(WP), ALLOCATABLE :: DVDL2MzL_io(:, :, :, :),    DVDL2MztL_io(:, :, :, :)


    REAL(WP), ALLOCATABLE :: UHzL_io(:, :, :), UHztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: GHzL_io(:, :, :), GHztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: U2DHzL_io(:, :, :), U2DHztL_io(:, :, :)

    REAL(WP), ALLOCATABLE :: DhDL1zL_io(:, :, :),         DhDL1ztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DhDLPzL_io(:, :, :),         DhDLPztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DTDLKzL_io(:, :, :),         DTDLKztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DTDLKUzL_io(:, :, :, :),      DTDLKUztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DTDLKDVDLzL_io(:, :, :, :, :), DTDLKDVDLztL_io(:, :, :, :, :)
    REAL(WP), ALLOCATABLE :: DhDLMDVDLzL_io(:, :, :, :, :), DhDLMDVDLztL_io(:, :, :, :, :)


    !=============== iO pARt (periodic) ==== Space averaged================
    REAL(WP), ALLOCATABLE :: U1xzL_io(:, :),  U1xztL_io(:, :) !Ui
    REAL(WP), ALLOCATABLE :: G1xzL_io(:, :),  G1xztL_io(:, :)!Gi
    REAL(WP), ALLOCATABLE :: UPxzL_io(:, :),  UPxztL_io(:, :)!UiP

    REAL(WP), ALLOCATABLE :: U2xzL_io(:, :),  U2xztL_io(:, :)!UiUj
    REAL(WP), ALLOCATABLE :: UGxzL_io(:, :),  UGxztL_io(:, :)!UiGj

    REAL(WP), ALLOCATABLE :: U3xzL_io(:, :),  U3xztL_io(:, :)!UiUjUk
    REAL(WP), ALLOCATABLE :: UGUxzL_io(:, :), UGUxztL_io(:, :)!UiGjUk

    REAL(WP), ALLOCATABLE :: DVDL1xzL_io(:, :, :), DVDL1xztL_io(:, :, :) ! dUI / DXj
    REAL(WP), ALLOCATABLE :: DVDLPxzL_io(:, :, :), DVDLPxztL_io(:, :, :) ! PdUI / DXj
    REAL(WP), ALLOCATABLE :: DVDL2xzL_io(:, :, :), DVDL2xztL_io(:, :, :) ! dUI / DXj * dUI / DXj

    REAL(WP), ALLOCATABLE :: T1xzL_io(:), T1xztL_io(:)
    REAL(WP), ALLOCATABLE :: D1xzL_io(:), D1xztL_io(:)
    REAL(WP), ALLOCATABLE :: H1xzL_io(:), H1xztL_io(:)
    REAL(WP), ALLOCATABLE :: M1xzL_io(:), M1xztL_io(:)

    REAL(WP), ALLOCATABLE :: T2xzL_io(:), T2xztL_io(:)
    REAL(WP), ALLOCATABLE :: D2xzL_io(:), D2xztL_io(:)
    REAL(WP), ALLOCATABLE :: H2xzL_io(:), H2xztL_io(:)

    REAL(WP), ALLOCATABLE :: DHxzL_io(:), DHxztL_io(:)
    REAL(WP), ALLOCATABLE :: PHxzL_io(:), PHxztL_io(:)

    REAL(WP), ALLOCATABLE :: DVDL1MxzL_io(:, :, :),    DVDL1MxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DVDL1MHxzL_io(:, :, :),   DVDL1MHxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DVDL1MUxzL_io(:, :, :, :), DVDL1MUxztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DVDL2MxzL_io(:, :, :),    DVDL2MxztL_io(:, :, :)

    REAL(WP), ALLOCATABLE :: UHxzL_io(:, :), UHxztL_io(:, :)
    REAL(WP), ALLOCATABLE :: GHxzL_io(:, :), GHxztL_io(:, :)
    REAL(WP), ALLOCATABLE :: U2DHxzL_io(:, :), U2DHxztL_io(:, :)

    REAL(WP), ALLOCATABLE :: DhDL1xzL_io(:, :),         DhDL1xztL_io(:, :)
    REAL(WP), ALLOCATABLE :: DhDLPxzL_io(:, :),         DhDLPxztL_io(:, :)
    REAL(WP), ALLOCATABLE :: DTDLKxzL_io(:, :),         DTDLKxztL_io(:, :)
    REAL(WP), ALLOCATABLE :: DTDLKUxzL_io(:, :, :),      DTDLKUxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: DTDLKDVDLxzL_io(:, :, :, :), DTDLKDVDLxztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: DhDLMDVDLxzL_io(:, :, :, :), DhDLMDVDLxztL_io(:, :, :, :)

    !================ variables ==================================================
    REAL(WP) :: DHwal_RA(2)
    REAL(WP) :: Hwal_RA(2)
    REAL(WP) :: Hwal_FA(2)
    REAL(WP) :: Twal(2)
    REAL(WP) :: Dwal(2) = 1.0_WP
    REAL(WP) :: Mwal(2) = 1.0_WP
    REAL(WP) :: Kwal(2)
    REAL(WP) :: Cpwal(2)
    REAL(WP) :: Qw(2)

    !======below values are non - DimensionaL ===============

    REAL(WP) :: Utaw_D_io(2), Utaw_io(2), Utaw_ave_io, Utaw_D_ave_io
    REAL(WP) :: Tauw_D_io(2), Tauw_io(2), Tauw_ave_io, Tauw_D_ave_io
    REAL(WP) :: Ret_io(2), Ldist_io(2)
    REAL(WP) :: Ret_ave_io, DenAvew, VisAvew


    !===below for quadrant analysIS =========
    INTEGER(4), PARAMETER :: QUADHN =9
    REAL(WP) :: QUADHV(QUADHN)

    !========== tg =================================
    REAL(WP), ALLOCATABLE :: QUADUVxzL_TG(:, :, :), QUADUVxztL_TG(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADVzxzL_TG(:, :, :), QUADVzxztL_TG(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADTKxzL_TG(:, :, :), QUADTKxztL_TG(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADDRxzL_TG(:, :, :), QUADDRxztL_TG(:, :, :)

    !========== Developping io======================
    REAL(WP), ALLOCATABLE :: QUADUVzL_io(:, :, :, :), QUADUVztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: QUADVzzL_io(:, :, :, :), QUADVzztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: QUADTKzL_io(:, :, :, :), QUADTKztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: QUADDRzL_io(:, :, :, :), QUADDRztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: QUADDUV1zL_io(:, :, :, :), QUADDUV1ztL_io(:, :, :, :)
    REAL(WP), ALLOCATABLE :: QUADDUV2zL_io(:, :, :, :), QUADDUV2ztL_io(:, :, :, :)

    !==========periodic io=========================
    REAL(WP), ALLOCATABLE :: QUADUVxzL_io(:, :, :), QUADUVxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADVzxzL_io(:, :, :), QUADVzxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADTKxzL_io(:, :, :), QUADTKxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADDRxzL_io(:, :, :), QUADDRxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADDUV1xzL_io(:, :, :), QUADDUV1xztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: QUADDUV2xzL_io(:, :, :), QUADDUV2xztL_io(:, :, :)

    REAL(WP), ALLOCATABLE :: OCTDUVxzL_io(:, :, :), OCTDUVxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTDVzxzL_io(:, :, :), OCTDVzxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTDTKxzL_io(:, :, :), OCTDTKxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTDDRxzL_io(:, :, :),  OCTDDRxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTDDUV1xzL_io(:, :, :), OCTDDUV1xztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTDDUV2xzL_io(:, :, :), OCTDDUV2xztL_io(:, :, :)

    REAL(WP), ALLOCATABLE :: OCTTUVxzL_io(:, :, :), OCTTUVxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTTVzxzL_io(:, :, :), OCTTVzxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTTTKxzL_io(:, :, :), OCTTTKxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTTDRxzL_io(:, :, :), OCTTDRxztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTTDUV1xzL_io(:, :, :), OCTTDUV1xztL_io(:, :, :)
    REAL(WP), ALLOCATABLE :: OCTTDUV2xzL_io(:, :, :), OCTTDUV2xztL_io(:, :, :)

    !========= DRIVEN FORCE ============================
    REAL(WP) :: FcDrv_io

    REAL(WP), ALLOCATABLE :: FUxzL_io(:, :)
    REAL(WP), ALLOCATABLE :: FUxztL_io(:, :)



    !=========== Spectra ======================
    !================ Velocity =====================
    REAL(WP), ALLOCATABLE :: R11X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R22X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R33X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R12X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R13X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R23X1_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: R11X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R22X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R33X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R12X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R13X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R23X3_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENE11T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE22T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE33T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE12T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE13T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE23T_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENE11Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE22Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE33Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE12Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE13Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE23Z_xzLa(:, :, :)

    !================ VoritICity ====================
    REAL(WP), ALLOCATABLE :: V11X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V22X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V33X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V12X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V13X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V23X1_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: V11X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V22X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V33X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V12X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V13X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V23X3_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENV11T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV22T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV33T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV12T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV13T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV23T_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENV11Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV22Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV33Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV12Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV13Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV23Z_xzLa(:, :, :)


    !=============== VoritICity & VelocitY =========================
    REAL(WP), ALLOCATABLE :: VO11X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO12X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO13X1_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO21X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO22X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO23X1_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO31X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO32X1_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO33X1_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO11X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO12X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO13X3_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO21X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO22X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO23X3_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO31X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO32X3_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO33X3_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO11T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO12T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO13T_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO21T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO22T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO23T_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO31T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO32T_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO33T_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO11Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO12Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO13Z_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO21Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO22Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO23Z_xzLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO31Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO32Z_xzLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO33Z_xzLa(:, :, :)

    !=========== Spectra ======================
    !================ Velocity =====================
    REAL(WP), ALLOCATABLE :: R11X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R22X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R33X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R12X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R13X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R23X1_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: R11X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R22X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R33X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R12X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R13X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: R23X3_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENE11T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE22T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE33T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE12T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE13T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE23T_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENE11Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE22Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE33Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE12Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE13Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENE23Z_xztLa(:, :, :)

    !================ VoritICity ====================
    REAL(WP), ALLOCATABLE :: V11X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V22X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V33X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V12X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V13X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V23X1_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: V11X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V22X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V33X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V12X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V13X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: V23X3_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENV11T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV22T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV33T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV12T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV13T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV23T_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: ENV11Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV22Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV33Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV12Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV13Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: ENV23Z_xztLa(:, :, :)


    !=============== VoritICity & VelocitY =========================
    REAL(WP), ALLOCATABLE :: VO11X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO12X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO13X1_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO21X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO22X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO23X1_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO31X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO32X1_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO33X1_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO11X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO12X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO13X3_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO21X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO22X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO23X3_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: VO31X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO32X3_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: VO33X3_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO11T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO12T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO13T_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO21T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO22T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO23T_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO31T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO32T_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO33T_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO11Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO12Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO13Z_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO21Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO22Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO23Z_xztLa(:, :, :)

    REAL(WP), ALLOCATABLE :: EVO31Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO32Z_xztLa(:, :, :)
    REAL(WP), ALLOCATABLE :: EVO33Z_xztLa(:, :, :)

END MODULE postprocess_info

!**********************************************************************************************************************************
! MODULE Flow_State_Constant
!     USE cparam
!
!     REAL(WP), PARAMETER :: SL_MU0 = 1.716E-5_WP  ! Kg/ M /s
!     REAL(WP), PARAMETER :: SL_T0  = 273.15_WP    ! K
!     REAL(WP), PARAMETER :: SL_S = 110.4_WP     ! K
!     REAL(WP), PARAMETER :: SL_C1  = 1.458E-6_WP  ! kg/ M /s/sqrt(K)
!
!     REAL(WP), PARAMETER :: SL_K0  = 0.023961_WP  ! W/ MK
!
! END MODULE Flow_State_Constant
