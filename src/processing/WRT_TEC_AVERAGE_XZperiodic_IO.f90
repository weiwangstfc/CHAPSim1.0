!**********************************************************************************************************************************
!> @brief
!>        postprocessing
!> @details
!> module: VARS_AVERAGED_XZ_io
!> subroutine: MEMO_ALLOCT_AVERAGE_XZ_io
!> subroutine: MEMO_DEALLT_AVERAGE_XZ_io
!> subroutine: WRT_AVERAGE_PPED_Xperiodic_io
!> subroutine: WRT_AVERAGE_PPED_XZ_io_GATHER
!> subroutine: PP_FLOW_BASIC_VARS_XZ_io
!> subroutine: PP_Budg_INIT
!> subroutine: PP_FLOW_FA_RSTE_Budg_XZ_io
!> subroutine: PP_FLOW_RA_noDen_RSTE_Budg_XZ_io
!> subroutine: PP_HEAT_BASIC_VARS_XZ_io
!> subroutine: PP_HEAT_FA_RSTE_Budg_XZ_io
!> subroutine: WRT_FLOW_FA_Profile_XZ_io
!> subroutine: WRT_FLOW_RA_Profile_XZ_io
!> subroutine: WRT_HEAT_FA_Profile_XZ_io
!> subroutine: WRT_HeatTransfer_Table_XZ_io
!> subroutine: WRT_Cf_Table_XZ_io
!> subroutine: WRT_Checking_TABLE_XZ_io
!> subroutine: PP_SSzero_SIDED
!> subroutine: PP_J4TbukTsd
!> subroutine: PP_Umax_SIDED
!> subroutine: WRT_FLOW_Budgets_Profile_XZ_io
!> subroutine: WRITE_SPECO_AVE_PROFILE
!> subroutine: WRITE_SPECO_AVE_Contour
!> @note
!> @todo
! REVISION HISTORY:
! 06/2014 - Created, by Wei Wang (wei.wang@sheffield.ac.uk)
!**********************************************************************************************************************************
!================================periodic BELOW========================================================
MODULE VARS_AVERAGED_XZ_io
    USE thermal_info
    USE WRT_INFO
    USE postprocess_info
    USE mesh_info
    USE init_info
    USE flow_info
    CHARACTER(15) :: PNTIM

    INTEGER(4), PARAMETER :: NX = 5

    !============= Averaged global Dwta in each processoR ======================
    REAL(WP), ALLOCATABLE :: U1xztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: G1xztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: UPxztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: U2xztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: UGxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: UGUxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: U3xztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: DVDL1xztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDLPxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDL2xztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: QuadUVxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: QuadVzxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: QuadTKxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: QuadDRxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: QuadDUV1xztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: QuadDUV2xztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: OctDUVxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctDVzxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctDTKxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctDDRxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctDDUV1xztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctDDUV2xztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: OctTUVxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctTVzxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctTTKxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctTDRxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctTDUV1xztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: OctTDUV2xztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: FUxztL_F0_io( :, :)

    !============= Averaged global Dwta in each processor thermaL ======================
    REAL(WP), ALLOCATABLE :: T1xztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: D1xztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: H1xztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: M1xztL_F0_io( : )

    REAL(WP), ALLOCATABLE :: T2xztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: D2xztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: H2xztL_F0_io( : )

    REAL(WP), ALLOCATABLE :: DHxztL_F0_io( : )
    REAL(WP), ALLOCATABLE :: PHxztL_F0_io( : )

    REAL(WP), ALLOCATABLE :: DVDL1MxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDL1MHxztL_F0_io( :, :, :  )
    REAL(WP), ALLOCATABLE :: DVDL1MUxztL_F0_io( :, :, :, :)
    REAL(WP), ALLOCATABLE :: DVDL2MxztL_F0_io( :, :, :  )

    REAL(WP), ALLOCATABLE :: UHxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: GHxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: U2DHxztL_F0_io( :, : )

    REAL(WP), ALLOCATABLE :: DhDL1xztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: DhDLPxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: DTDLKxztL_F0_io( :, : )
    REAL(WP), ALLOCATABLE :: DTDLKUxztL_F0_io( :, :, :)
    REAL(WP), ALLOCATABLE :: DTDLKDVDLxztL_F0_io( :, :, :, :)
    REAL(WP), ALLOCATABLE :: DHDLMDVDLxztL_F0_io( :, :, :, :)

    !===========================================
    REAL(WP), ALLOCATABLE :: TauwSD(:)
    REAL(WP), ALLOCATABLE :: DensSD(:)
    REAL(WP), ALLOCATABLE :: YWdISD(:)
    REAL(WP), ALLOCATABLE :: viscsD(:)

    REAL(WP), ALLOCATABLE :: CpSD(:)
    REAL(WP), ALLOCATABLE :: QwSD(:)
    REAL(WP), ALLOCATABLE :: TwSD(:)
    REAL(WP), ALLOCATABLE :: HwSD(:)

    !========================pp intermediate variables ===============================
    !==== RA related==============
    REAL(WP), ALLOCATABLE :: dPDX_RA   (:, :)          !d<p>/ DX_i
    REAL(WP), ALLOCATABLE :: ufpf_RA   (:, :)          !<p' u'_i> = <p' u"_i>
    REAL(WP), ALLOCATABLE :: uf2_RA    (:, :, :)       !<u'_i u'_j>
    REAL(WP), ALLOCATABLE :: uf2d_RA   (:, :, :)       !<\rho> <u'_i u'_j>
    REAL(WP), ALLOCATABLE :: uf3_RA    (:, :, :, :)    !<u'_i u'_j u'_k>
    REAL(WP), ALLOCATABLE :: uf3d_RA   (:, :, :, :)    !<\rho> * <u'_i u'_j u'_k>
    REAL(WP), ALLOCATABLE :: UU_RA      (:, :, :)    !{u_i u_j}

    REAL(WP), ALLOCATABLE :: dUiDXi    (:)             ! dUI / DXI
    REAL(WP), ALLOCATABLE :: StrAInTensor(:, :, :)     ! 0.5_WP * (dUI / DXJ + DUJ / DXI)
    REAL(WP), ALLOCATABLE :: VortcyTensor(:, :, :)     ! 0.5_WP * (dUI / DXJ - DUJ / DXI)
    REAL(WP), ALLOCATABLE :: Skewness_RA(:, :)

    REAL(WP), ALLOCATABLE :: MKE_RA(:)
    REAL(WP), ALLOCATABLE :: TKE_RA(:)
    REAL(WP), ALLOCATABLE :: ufTKEfd_RA(:)
    REAL(WP), ALLOCATABLE :: ufMKEfd_RA(:)

    REAL(WP), ALLOCATABLE :: Omega2_RA(:, :)
    REAL(WP), ALLOCATABLE :: Omega_RA2(:, :)
    REAL(WP), ALLOCATABLE :: Omega_rms(:, :)

    REAL(WP), ALLOCATABLE :: AnIsotropy_RA(:, :, :)
    REAL(WP), ALLOCATABLE :: Anistpinva_RA(:, :)
    REAL(WP), ALLOCATABLE :: LumleyAxis_RA(:, :)

    !====FA related==============
    REAL(WP), ALLOCATABLE :: DrivenForce (:)         !

    REAL(WP), ALLOCATABLE :: U_FA       (:, :)       !{u_i}
    REAL(WP), ALLOCATABLE :: UU_FA      (:, :, :)    !{u_i u_j}
    REAL(WP), ALLOCATABLE :: dUDX_FA    (:, :, :)    !d{u_i}/ DX_j

    REAL(WP), ALLOCATABLE :: uff_RA     (:, :)       ! <u"_i>
    REAL(WP), ALLOCATABLE :: uff2_FA    (:, :, :)    !{u"_i u"_j}
    REAL(WP), ALLOCATABLE :: uff2d_FA   (:, :, :)    !<\rho>*{u"_i u"_j}
    REAL(WP), ALLOCATABLE :: uff3_FA    (:, :, :, :) !{u"_i u"_j u"_k}
    REAL(WP), ALLOCATABLE :: uff3d_FA   (:, :, :, :) !<\rho>*{u"_i u"_j u"_k}
    REAL(WP), ALLOCATABLE :: TDIFU_FA   (:, :, :)    !
    REAL(WP), ALLOCATABLE :: dUiDXiM      (:)        ! <\mu dUI / DXI>
    REAL(WP), ALLOCATABLE :: StrAInTensorM(:, :, :)  ! < 0.5_WP * (dUI / DXJ + DUJ / DXI) *\mu>
    REAL(WP), ALLOCATABLE :: VortcyTensorM(:, :, :)  ! < 0.5_WP * (dUI / DXJ - DUJ / DXI) *\mu>
    REAL(WP), ALLOCATABLE :: Skewness_FA(:, :)

    REAL(WP), ALLOCATABLE :: MKE_FA(:)
    REAL(WP), ALLOCATABLE :: TKE_FA(:)
    REAL(WP), ALLOCATABLE :: uffTKEffd_FA(:)
    REAL(WP), ALLOCATABLE :: uffMKEffd_FA(:)

    REAL(WP), ALLOCATABLE :: Tau_Mean_RA (:, :, :)       !<tau_ij>
    REAL(WP), ALLOCATABLE :: Tau_meaU_RA (:, :, :)       !<tau_ij>(<S>,<\mu>)
    REAL(WP), ALLOCATABLE :: dTaudy_RA   (:, :, :)       !d<tau_ij>/ Dy
    REAL(WP), ALLOCATABLE :: dTSSdy_RA   (:, :, :)       !d<R_ij>/ Dy
    REAL(WP), ALLOCATABLE :: TauU_RA     (:, :, :, :)    !< u_h* \tau_mn >
    REAL(WP), ALLOCATABLE :: Taufuf_RA   (:, :, :, :)    !< u'_h* \tau'_mn >
    REAL(WP), ALLOCATABLE :: TauDvDL_RA  (:, :, :, :, :) !< d(u_m) / D(x_n) * \tau_hp >
    !REAL(WP), ALLOCATABLE :: Tau_ik_Du_jDX_i_RA  (:, :, :) !< d(u_m) / D(x_n) * \tau_hp >


    REAL(WP), ALLOCATABLE :: AnIsotropy_FA(:, :, :)
    REAL(WP), ALLOCATABLE :: Anistpinva_FA(:, :)
    REAL(WP), ALLOCATABLE :: LumleyAxis_FA(:, :)

    !==================Checking==============================================
    REAL(WP) :: NSFbalt_RA(NDV)
    REAL(WP) :: NSFbalt_FA(NDV)
    REAL(WP), ALLOCATABLE :: NSFbal_RA (:, :)
    REAL(WP), ALLOCATABLE :: NSFbal_FA (:, :)
    REAL(WP) :: BuoyForceTT
    REAL(WP), ALLOCATABLE :: ENEbal_FA (:)
    REAL(WP) :: ENEbalt_FA

    REAL(WP), ALLOCATABLE :: bdfcintg (:)
    REAL(WP), ALLOCATABLE :: DensIntg (:)

    !===============budgetS ================================
    !================== RuV =============================
    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_pduDX_stran_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_DpuDX_diffu_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_Duiuj(:, :)

    REAL(WP), ALLOCATABLE :: Budg_press_accl1_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc1_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_Duiuj(:, :)

    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc2_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc2_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_turss_accl2_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Balance2_Duiuj(:, :)

    REAL(WP), ALLOCATABLE :: Budg_pressure3_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_vistress3_Duiuj(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Balance3_Duiuj(:, :)

    !================== RuV == Sum along Y =======================
    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_pduDX_stran_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_DpuDX_diffu_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_Duiuj_ysum(:)

    REAL(WP), ALLOCATABLE :: Budg_press_accl1_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc1_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_Duiuj_ysum(:)

    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc2_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc2_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_turss_accl2_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance2_Duiuj_ysum(:)

    REAL(WP), ALLOCATABLE :: Budg_pressure3_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_vistress3_Duiuj_ysum(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance3_Duiuj_ysum(:)

    !================== TKE =============================
    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_pduDX_stran_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_DpuDX_diffu_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_TKE(:)

    REAL(WP), ALLOCATABLE :: Budg_press_accl1_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc1_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_TKE(:)

    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc2_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc2_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_turss_accl2_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance2_TKE(:)

    REAL(WP), ALLOCATABLE :: Budg_pressure3_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_vistress3_TKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance3_TKE(:)

    !================== TKE sum along Y =============================
    REAL(WP) :: Budg_prodc_stres_TKE_ysum
    REAL(WP) :: Budg_viscs_dissp_TKE_ysum
    REAL(WP) :: Budg_pduDX_stran_TKE_ysum
    REAL(WP) :: Budg_Turbu_diffu_TKE_ysum
    REAL(WP) :: Budg_DpuDX_diffu_TKE_ysum
    REAL(WP) :: Budg_viscs_diffu_TKE_ysum

    REAL(WP) :: Budg_press_accl1_TKE_ysum
    REAL(WP) :: Budg_viscs_accl1_TKE_ysum
    REAL(WP) :: Budg_prodc_Dvfc1_TKE_ysum
    REAL(WP) :: Budg_Balance1_TKE_ysum

    REAL(WP) :: Budg_prodc_gvfc2_TKE_ysum
    REAL(WP) :: Budg_prodc_Dvfc2_TKE_ysum
    REAL(WP) :: Budg_turss_accl2_TKE_ysum
    REAL(WP) :: Budg_Balance2_TKE_ysum

    REAL(WP) :: Budg_pressure3_TKE_ysum
    REAL(WP) :: Budg_vistress3_TKE_ysum
    REAL(WP) :: Budg_Balance3_TKE_ysum

    !================== MKE =============================
    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_pduDX_stran_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_DpuDX_diffu_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_MKE(:)

    REAL(WP), ALLOCATABLE :: Budg_press_accl1_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc1_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc1_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_MKE(:)

    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc2_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_Dvfc2_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_turss_accl2_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance2_MKE(:)

    REAL(WP), ALLOCATABLE :: Budg_pressure3_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_vistress3_MKE(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance3_MKE(:)

    !================== MKE sum along Y =============================
    REAL(WP) :: Budg_prodc_stres_MKE_ysum
    REAL(WP) :: Budg_viscs_dissp_MKE_ysum
    REAL(WP) :: Budg_pduDX_stran_MKE_ysum
    REAL(WP) :: Budg_Turbu_diffu_MKE_ysum
    REAL(WP) :: Budg_DpuDX_diffu_MKE_ysum
    REAL(WP) :: Budg_viscs_diffu_MKE_ysum

    REAL(WP) :: Budg_press_accl1_MKE_ysum
    REAL(WP) :: Budg_viscs_accl1_MKE_ysum
    REAL(WP) :: Budg_prodc_Dvfc1_MKE_ysum
    REAL(WP) :: Budg_prodc_gvfc1_MKE_ysum
    REAL(WP) :: Budg_Balance1_MKE_ysum

    REAL(WP) :: Budg_prodc_gvfc2_MKE_ysum
    REAL(WP) :: Budg_prodc_Dvfc2_MKE_ysum
    REAL(WP) :: Budg_turss_accl2_MKE_ysum
    REAL(WP) :: Budg_Balance2_MKE_ysum

    REAL(WP) :: Budg_pressure3_MKE_ysum
    REAL(WP) :: Budg_vistress3_MKE_ysum
    REAL(WP) :: Budg_Balance3_MKE_ysum


    !==============For RANS =================================
    REAL(WP), ALLOCATABLE :: RANS_Mut(:)

    !========================pp variables 2 ===============================
    REAL(WP), ALLOCATABLE :: H_FA(:)
    REAL(WP), ALLOCATABLE :: hff_RA(:)
    REAL(WP), ALLOCATABLE :: hfpf_RA(:)
    REAL(WP), ALLOCATABLE :: dTDX(:, :)
    REAL(WP), ALLOCATABLE :: dDDX(:, :)
    REAL(WP), ALLOCATABLE :: dHDX_RA(:, :)
    REAL(WP), ALLOCATABLE :: dHDX_FA(:, :)
    REAL(WP), ALLOCATABLE :: UH_FA(:, :)
    REAL(WP), ALLOCATABLE :: uff2hffd_FA(:, :, :)

    REAL(WP), ALLOCATABLE :: uffhffd_FA(:, :)
    REAL(WP), ALLOCATABLE :: ufhfd_RA(:, :)
    REAL(WP), ALLOCATABLE :: viscstressEnth_RA(:, :, :)
    REAL(WP), ALLOCATABLE :: viscstressEnthGrad_RA(:, :, :, :)

    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_prodc_enthg_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_press_accl1_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_DphDX_diffu_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_pdhDX_stran_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_accel_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_diffu_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_dissp_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_thf(:, :)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_thf(:, :)

    REAL(WP), ALLOCATABLE :: Budg_prodc_stres_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_prodc_enthg_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_Turbu_diffu_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_press_accl1_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_DphDX_diffu_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_pdhDX_stran_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_accel_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_diffu_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_ConHF_dissp_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_accl1_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_diffu_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_viscs_dissp_IEN(:)
    REAL(WP), ALLOCATABLE :: Budg_Balance1_IEN(:)

    !========================body forcE =================================

    REAL(WP), ALLOCATABLE :: Budg_prodc_gvfc2_thf(:)

    !========================= ThermaL ===================================
    INTEGER(4) :: J4SS0, J4Tpc, J4Tbk, J4MaxU, J4TbukSsd(2)

    !======on the wall ==================
    REAL(WP) :: Hwal_RA_D(2)
    REAL(WP) :: Twal_D(2)
    REAL(WP) :: Dwal_D(2)
    REAL(WP) :: Mwal_D(2)
    REAL(WP) :: Kwal_D(2)
    REAL(WP) :: Cpwal_D(2)
    REAL(WP) :: Bowal_D(2)
    REAL(WP) :: qw_D(2), qw_ave, qw_D_ave

    !=======coeffecieinT ================
    REAL(WP) :: Cf0_io(2),        Cf0_ave_io
    REAL(WP) :: Cfbk_io(2),       Cfbk_ave_io
    REAL(WP) :: CfbkSsd_io(2),    CfbkSsd_ave_io

    REAL(WP) :: Rebk, RebkSsd(2), RebkTsd(2)
    REAL(WP) :: Prbk, PrbkSsd(2), PrbkTsd(2)
    REAL(WP) :: Grbk(2), Grbk_Drho, GrbkSsd(2)
    REAL(WP) :: Bobk(2)
    REAL(WP) :: Nubk(2), Nubk_DTw(2), NubkSsd(2), NubkTsd(2), Nu_int

    REAL(WP) :: hc_D(2),  hc_DTw_D(2), hc_DTw_D_ave, hcsd_D(2)


    REAL(WP) :: L4TbkSsd(2), L4TbkTsd(2), L4Tbk(2)


    REAL(WP) :: RIChARdsoNNo(2), Ttau(2), Ttau_D(2)
    REAL(WP) :: MDOt
    REAL(WP) :: DHDOt
    REAL(WP) :: Gbuk, Gbuk_D, GbukSsd(2), GbukTsd(2)
    REAL(WP) :: Hbuk, Hbuk_D, HbukSsd(2), HbukTsd(2)
    REAL(WP) :: DHbuk, DHbuk_D, DHbukSsd(2), DHbukTsd(2)
    REAL(WP) :: Ubuk, Ubuk_D, UbukSsd(2), UbukTsd(2)
    REAL(WP) :: Tbuk, Tbuk_D, TbukSsd(2), TbukTsd(2)
    REAL(WP) :: Mbuk, Mbuk_D, MbukSsd(2), MbukTsd(2)
    REAL(WP) :: Kbuk, Kbuk_D, KbukSsd(2), KbukTsd(2)
    REAL(WP) :: Cpbk, Cpbk_D, CpbkSsd(2), CpbkTsd(2)
    REAL(WP) :: Bbuk, Bbuk_D, BbukSsd(2), BbukTsd(2)
    REAL(WP) :: Dbuk, Dbuk_D, DbukSsd(2), DbukTsd(2)
    REAL(WP) :: D_int, D_int_D, Dsd_int(2)


    !REAL(WP), ALLOCATABLE :: Nuy(:, :)



END MODULE

!**********************************************************************************************************************************
SUBROUTINE MEMO_ALLOCT_AVERAGE_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    !========================================
    ALLOCATE( U1xztL_F0_io( NCL2, NDV + 1 ) ) ;  U1xztL_F0_io = 0.0_WP
    ALLOCATE( G1xztL_F0_io( NCL2, NDV   ) ) ; G1xztL_io = 0.0_WP
    ALLOCATE( UPxztL_F0_io( NCL2, NDV   ) ) ; UPxztL_io = 0.0_WP

    ALLOCATE( U2xztL_F0_io( NCL2, NDV * (7 - NDV) / 2 + NDV - 3 ) ) ; U2xztL_F0_io = 0.0_WP
    ALLOCATE( UGxztL_F0_io( NCL2, NDV * (7 - NDV) / 2 + NDV - 3 ) ) ; UGxztL_F0_io = 0.0_WP
    ALLOCATE( UGUxztL_F0_io(NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) ) ; UGUxztL_F0_io = 0.0_WP
    ALLOCATE( U3xztL_F0_io(NCL2, NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8) ) ;  U3xztL_F0_io = 0.0_WP

    ALLOCATE( DVDL1xztL_F0_io( NCL2, NDV, NDV  ) ) ; DVDL1xztL_F0_io = 0.0_WP
    ALLOCATE( DVDLPxztL_F0_io( NCL2, NDV, NDV  ) ) ; DVDLPxztL_F0_io = 0.0_WP
    ALLOCATE( DVDL2xztL_F0_io( NCL2, (NDV - 1) * 3 + NDV, (NDV - 1) * 3 + NDV  ) ) ; DVDL2xztL_F0_io = 0.0_WP
    IF(iPPQuadrants == 1)  THEN
        ALLOCATE( QuadUVxztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadUVxztL_F0_io = 0.0_WP
        ALLOCATE( QuadVzxztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadVzxztL_F0_io = 0.0_WP
        ALLOCATE( QuadTKxztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadTKxztL_F0_io = 0.0_WP
        ALLOCATE( QuadDRxztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadDRxztL_F0_io = 0.0_WP
        ALLOCATE( QuadDUV1xztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadDUV1xztL_F0_io = 0.0_WP
        ALLOCATE( QuadDUV2xztL_F0_io( NCL2, 4, QUADHN  ) ) ; QuadDUV2xztL_F0_io = 0.0_WP

        ALLOCATE( OctDUVxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDUVxztL_F0_io = 0.0_WP
        ALLOCATE( OctDVzxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDVzxztL_F0_io = 0.0_WP
        ALLOCATE( OctDTKxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDTKxztL_F0_io = 0.0_WP
        ALLOCATE( OctDDRxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDDRxztL_F0_io = 0.0_WP
        ALLOCATE( OctDDUV1xztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDDUV1xztL_F0_io = 0.0_WP
        ALLOCATE( OctDDUV2xztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctDDUV2xztL_F0_io = 0.0_WP

        ALLOCATE( OctTUVxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTUVxztL_F0_io = 0.0_WP
        ALLOCATE( OctTVzxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTVzxztL_F0_io = 0.0_WP
        ALLOCATE( OctTTKxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTTKxztL_F0_io = 0.0_WP
        ALLOCATE( OctTDRxztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTDRxztL_F0_io = 0.0_WP
        ALLOCATE( OctTDUV1xztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTDUV1xztL_F0_io = 0.0_WP
        ALLOCATE( OctTDUV2xztL_F0_io( NCL2, 8, QUADHN  ) ) ; OctTDUV2xztL_F0_io = 0.0_WP
    END IF
    ALLOCATE( FUxztL_F0_io( NCL2, NDV + 1 ) ) ; FUxztL_F0_io = 0.0_WP

    !============= RA =========================
    ALLOCATE( dPDX_RA   (NCL2, NDV)    ) ; dPDX_RA = 0.0_WP
    ALLOCATE( ufpf_RA   (NCL2, NDV)    ) ; ufpf_RA = 0.0_WP
    ALLOCATE( uf2_RA    (NCL2, NDV, NDV)    ) ; uf2_RA = 0.0_WP
    ALLOCATE( uf2d_RA   (NCL2, NDV, NDV)    ) ; uf2d_RA = 0.0_WP
    ALLOCATE( uf3_RA    (NCL2, NDV, NDV, NDV) ) ; uf3_RA = 0.0_WP
    ALLOCATE( uf3d_RA   (NCL2, NDV, NDV, NDV) ) ; uf3d_RA = 0.0_WP
    ALLOCATE( UU_RA         (NCL2, NDV, NDV) ) ; UU_RA = 0.0_WP

    ALLOCATE( dUiDXi        (NCL2)            ) ; dUiDXi = 0.0_WP
    ALLOCATE( StrAInTensor  (NCL2, NDV, NDV)  ) ; StrAInTensor = 0.0_WP
    ALLOCATE( VortcyTensor  (NCL2, NDV, NDV)  ) ; VortcyTensor = 0.0_WP
    ALLOCATE( Skewness_RA   (NCL2, NDV)       ) ; Skewness_RA = 0.0_WP

    ALLOCATE( MKE_RA        (NCL2)            ) ; MKE_RA = 0.0_WP
    ALLOCATE( TKE_RA        (NCL2)            ) ; TKE_RA = 0.0_WP
    ALLOCATE( ufTKEfd_RA    (NCL2)            ) ; ufTKEfd_RA = 0.0_WP
    ALLOCATE( ufMKEfd_RA    (NCL2)            ) ; ufMKEfd_RA = 0.0_WP


    ALLOCATE( Omega2_RA(NCL2, NDV) ) ; Omega2_RA = 0.0_WP
    ALLOCATE( Omega_RA2(NCL2, NDV) ) ; Omega_RA2 = 0.0_WP
    ALLOCATE( Omega_rms(NCL2, NDV) ) ; Omega_rms = 0.0_WP

    ALLOCATE( AnIsotropy_RA(NCL2, NDV, NDV) ) ; AnIsotropy_RA = 0.0_WP
    ALLOCATE( Anistpinva_RA(NCL2, NDV) ) ; Anistpinva_RA = 0.0_WP
    ALLOCATE( LumleyAxis_RA(NCL2, 2) ) ; LumleyAxis_RA = 0.0_WP



    !=============FA =========================
    IF(iThermoDynamics == 1) THEN
        ALLOCATE( DrivenForce   (NCL2)           ) ; DrivenForce = 0.0_WP
        ALLOCATE( U_FA          (NCL2, NDV)      ) ; U_FA = 0.0_WP
        ALLOCATE( UU_FA         (NCL2, NDV, NDV) ) ; UU_FA = 0.0_WP
        ALLOCATE( dUDX_FA       (NCL2, NDV, NDV) ) ; dUDX_FA = 0.0_WP


        ALLOCATE( uff_RA        (NCL2, NDV)      ) ; uff_RA = 0.0_WP
        ALLOCATE( uff2_FA       (NCL2, NDV, NDV) ) ; uff2_FA = 0.0_WP
        ALLOCATE( uff2d_FA      (NCL2, NDV, NDV) ) ; uff2d_FA = 0.0_WP
        ALLOCATE( uff3_FA       (NCL2, NDV, NDV, NDV) ) ; uff3_FA = 0.0_WP
        ALLOCATE( uff3d_FA      (NCL2, NDV, NDV, NDV) ) ; uff3d_FA = 0.0_WP
        ALLOCATE( TDIFU_FA      (NCL2, NDV, NDV) ) ; TDIFU_FA = 0.0_WP

        ALLOCATE( dUiDXiM       (NCL2)            ) ; dUiDXiM = 0.0_WP
        ALLOCATE( StrAInTensorM (NCL2, NDV, NDV)  ) ; StrAInTensorM = 0.0_WP
        ALLOCATE( VortcyTensorM (NCL2, NDV, NDV)  ) ; VortcyTensorM = 0.0_WP
        ALLOCATE( Skewness_FA   (NCL2, NDV)       ) ; Skewness_FA = 0.0_WP

        ALLOCATE( MKE_FA        (NCL2)            ) ; MKE_FA = 0.0_WP
        ALLOCATE( TKE_FA        (NCL2)            ) ; TKE_FA = 0.0_WP
        ALLOCATE( uffTKEffd_FA  (NCL2)            ) ; uffTKEffd_FA = 0.0_WP
        ALLOCATE( uffMKEffd_FA  (NCL2)            ) ; uffMKEffd_FA = 0.0_WP

        ALLOCATE( AnIsotropy_FA(NCL2, NDV, NDV) ) ; AnIsotropy_FA = 0.0_WP
        ALLOCATE( Anistpinva_FA(NCL2, NDV) ) ; Anistpinva_FA = 0.0_WP
        ALLOCATE( LumleyAxis_FA(NCL2, 2) ) ; LumleyAxis_FA = 0.0_WP
    END IF

    ALLOCATE( Tau_Mean_RA   (NCL2, NDV, NDV) ) ; Tau_Mean_RA = 0.0_WP
    ALLOCATE( Tau_meaU_RA   (NCL2, NDV, NDV) ) ; Tau_meaU_RA = 0.0_WP
    ALLOCATE( dTaudy_RA     (NCL2, NDV, NDV) ) ; dTaudy_RA = 0.0_WP
    ALLOCATE( dTSSdy_RA     (NCL2, NDV, NDV) ) ; dTSSdy_RA = 0.0_WP
    ALLOCATE( TauU_RA       (NCL2, NDV, NDV, NDV) ) ; TauU_RA = 0.0_WP
    ALLOCATE( Taufuf_RA     (NCL2, NDV, NDV, NDV     ) ) ; Taufuf_RA = 0.0_WP
    ALLOCATE( TauDvDL_RA    (NCL2, NDV, NDV, NDV, NDV) ) ; TauDvDL_RA = 0.0_WP
    !ALLOCATE( Tau_ik_Du_jDX_i_RA(NCL2, NDV, NDV))

    ALLOCATE( NSFbal_FA(NCL2, NDV)) ; NSFbal_FA = 0.0_WP
    ALLOCATE( NSFbal_RA(NCL2, NDV)) ; NSFbal_RA = 0.0_WP
    ALLOCATE( ENEbal_FA(NCL2)) ; ENEbal_FA = 0.0_WP

    ALLOCATE( Bdfcintg(NCL2)) ; bdfciNTG = 0.0_WP
    ALLOCATE( DensIntg(NCL2)) ; DensIntg = 0.0_WP


    !========================================
    !================== RuV =============================
    ALLOCATE (Budg_prodc_stres_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_stres_Duiuj = 0.0_WP
    ALLOCATE (Budg_viscs_dissp_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_dissp_Duiuj = 0.0_WP
    ALLOCATE (Budg_pduDX_stran_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_pduDX_stran_Duiuj = 0.0_WP
    ALLOCATE (Budg_Turbu_diffu_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Turbu_diffu_Duiuj = 0.0_WP
    ALLOCATE (Budg_DpuDX_diffu_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_DpuDX_diffu_Duiuj = 0.0_WP
    ALLOCATE (Budg_viscs_diffu_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_diffu_Duiuj = 0.0_WP

    ALLOCATE (Budg_press_accl1_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_press_accl1_Duiuj = 0.0_WP
    ALLOCATE (Budg_viscs_accl1_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_accl1_Duiuj = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc1_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_Dvfc1_Duiuj = 0.0_WP
    ALLOCATE (Budg_Balance1_Duiuj   (NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance1_Duiuj = 0.0_WP

    ALLOCATE (Budg_prodc_gvfc2_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_gvfc2_Duiuj = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc2_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_Dvfc2_Duiuj = 0.0_WP
    ALLOCATE (Budg_turss_accl2_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_turss_accl2_Duiuj = 0.0_WP
    ALLOCATE (Budg_Balance2_Duiuj   (NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance2_Duiuj = 0.0_WP

    ALLOCATE (Budg_pressure3_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_pressure3_Duiuj = 0.0_WP
    ALLOCATE (Budg_vistress3_Duiuj(NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_vistress3_Duiuj = 0.0_WP
    ALLOCATE (Budg_Balance3_Duiuj (NCL2, NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance3_Duiuj = 0.0_WP

    !================== RuV == Sum along Y =======================
    ALLOCATE (Budg_prodc_stres_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_stres_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_viscs_dissp_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_dissp_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_pduDX_stran_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_pduDX_stran_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_Turbu_diffu_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Turbu_diffu_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_DpuDX_diffu_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_DpuDX_diffu_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_viscs_diffu_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_diffu_Duiuj_ysum = 0.0_WP

    ALLOCATE (Budg_press_accl1_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_press_accl1_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_viscs_accl1_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_viscs_accl1_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc1_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_Dvfc1_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_Balance1_Duiuj_ysum   (NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance1_Duiuj_ysum = 0.0_WP

    ALLOCATE (Budg_prodc_gvfc2_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_gvfc2_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc2_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_prodc_Dvfc2_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_turss_accl2_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_turss_accl2_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_Balance2_Duiuj_ysum   (NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance2_Duiuj_ysum = 0.0_WP

    ALLOCATE (Budg_pressure3_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_pressure3_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_vistress3_Duiuj_ysum(NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_vistress3_Duiuj_ysum = 0.0_WP
    ALLOCATE (Budg_Balance3_Duiuj_ysum (NDV * (7 - NDV) / 2 + NDV - 3)) ; Budg_Balance3_Duiuj_ysum = 0.0_WP

    !================== TKE =============================
    ALLOCATE (Budg_prodc_stres_TKE(NCL2) ) ; Budg_prodc_stres_TKE = 0.0_WP
    ALLOCATE (Budg_viscs_dissp_TKE(NCL2) ) ; Budg_viscs_dissp_TKE = 0.0_WP
    ALLOCATE (Budg_pduDX_stran_TKE(NCL2) ) ; Budg_pduDX_stran_TKE = 0.0_WP
    ALLOCATE (Budg_Turbu_diffu_TKE(NCL2) ) ; Budg_Turbu_diffu_TKE = 0.0_WP
    ALLOCATE (Budg_DpuDX_diffu_TKE(NCL2) ) ; Budg_DpuDX_diffu_TKE = 0.0_WP
    ALLOCATE (Budg_viscs_diffu_TKE(NCL2) ) ; Budg_viscs_diffu_TKE = 0.0_WP

    ALLOCATE (Budg_press_accl1_TKE(NCL2) ) ; Budg_press_accl1_TKE = 0.0_WP
    ALLOCATE (Budg_viscs_accl1_TKE(NCL2) ) ; Budg_viscs_accl1_TKE = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc1_TKE(NCL2) ) ; Budg_prodc_Dvfc1_TKE = 0.0_WP
    ALLOCATE (Budg_Balance1_TKE   (NCL2) ) ; Budg_Balance1_TKE = 0.0_WP

    ALLOCATE (Budg_prodc_gvfc2_TKE(NCL2) ) ; Budg_prodc_gvfc2_TKE = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc2_TKE(NCL2) ) ; Budg_prodc_Dvfc2_TKE = 0.0_WP
    ALLOCATE (Budg_turss_accl2_TKE(NCL2) ) ; Budg_turss_accl2_TKE = 0.0_WP
    ALLOCATE (Budg_Balance2_TKE   (NCL2) ) ; Budg_Balance2_TKE = 0.0_WP

    ALLOCATE (Budg_pressure3_TKE(NCL2) ) ; Budg_pressure3_TKE = 0.0_WP
    ALLOCATE (Budg_vistress3_TKE(NCL2) ) ; Budg_vistress3_TKE = 0.0_WP
    ALLOCATE (Budg_Balance3_TKE (NCL2) ) ; Budg_Balance3_TKE = 0.0_WP

    !================== MKE =============================
    ALLOCATE (Budg_prodc_stres_MKE(NCL2) ) ; Budg_prodc_stres_MKE = 0.0_WP
    ALLOCATE (Budg_viscs_dissp_MKE(NCL2) ) ; Budg_viscs_dissp_MKE = 0.0_WP
    ALLOCATE (Budg_pduDX_stran_MKE(NCL2) ) ; Budg_pduDX_stran_MKE = 0.0_WP
    ALLOCATE (Budg_Turbu_diffu_MKE(NCL2) ) ; Budg_Turbu_diffu_MKE = 0.0_WP
    ALLOCATE (Budg_DpuDX_diffu_MKE(NCL2) ) ; Budg_DpuDX_diffu_MKE = 0.0_WP
    ALLOCATE (Budg_viscs_diffu_MKE(NCL2) ) ; Budg_viscs_diffu_MKE = 0.0_WP

    ALLOCATE (Budg_press_accl1_MKE(NCL2) ) ; Budg_press_accl1_MKE = 0.0_WP
    ALLOCATE (Budg_viscs_accl1_MKE(NCL2) ) ; Budg_viscs_accl1_MKE = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc1_MKE(NCL2) ) ; Budg_prodc_Dvfc1_MKE = 0.0_WP
    ALLOCATE (Budg_prodc_gvfc1_MKE(NCL2) ) ; Budg_prodc_gvfc1_MKE = 0.0_WP
    ALLOCATE (Budg_Balance1_MKE   (NCL2) ) ; Budg_Balance1_MKE = 0.0_WP

    ALLOCATE (Budg_prodc_gvfc2_MKE(NCL2) ) ; Budg_prodc_gvfc2_MKE = 0.0_WP
    ALLOCATE (Budg_prodc_Dvfc2_MKE(NCL2) ) ; Budg_prodc_Dvfc2_MKE = 0.0_WP
    ALLOCATE (Budg_turss_accl2_MKE(NCL2) ) ; Budg_turss_accl2_MKE = 0.0_WP
    ALLOCATE (Budg_Balance2_MKE   (NCL2) ) ; Budg_Balance2_MKE = 0.0_WP

    ALLOCATE (Budg_pressure3_MKE(NCL2) ) ; Budg_pressure3_MKE = 0.0_WP
    ALLOCATE (Budg_vistress3_MKE(NCL2) ) ; Budg_vistress3_MKE = 0.0_WP
    ALLOCATE (Budg_Balance3_MKE (NCL2) ) ; Budg_Balance3_MKE = 0.0_WP


    !===========================================
    ALLOCATE ( RANS_Mut(NCL2) ) ; RANS_Mut = 0.0_WP
    !========================================
    ALLOCATE( TauwSD(NCL2)         ) ; TauwSD = 0.0_WP
    ALLOCATE( DensSD(NCL2)         ) ; DensSD = 0.0_WP
    ALLOCATE( viscsD(NCL2)         ) ; viscsD = 0.0_WP
    ALLOCATE( YWdISD(NCL2)         ) ; YWdISD = 0.0_WP

    ALLOCATE( D1xztL_F0_io( NCL2 ) ) ;  D1xztL_F0_io = 1.0_WP
    ALLOCATE( M1xztL_F0_io( NCL2 ) ) ;  M1xztL_F0_io = 1.0_WP

    ALLOCATE( DVDL1MxztL_F0_io( NCL2, NDV, NDV  ) )  ; DVDL1MxztL_F0_io = 0.0_WP
    ALLOCATE( DVDL1MUxztL_F0_io( NCL2, NDV, NDV, NDV  ) ) ; DVDL1MUxztL_F0_io = 0.0_WP
    ALLOCATE( DVDL2MxztL_F0_io( NCL2, (NDV - 1) * 3 + NDV, (NDV - 1) * 3 + NDV  ) ) ; DVDL2MxztL_F0_io = 0.0_WP

    IF(iThermoDynamics == 1) THEN
        !========================================
        ALLOCATE( T1xztL_F0_io( NCL2 ) ) ; T1xztL_F0_io = 0.0_WP
        ALLOCATE( H1xztL_F0_io( NCL2 ) ) ; H1xztL_F0_io = 0.0_WP

        ALLOCATE( DVDL1MHxztL_F0_io( NCL2, NDV, NDV  ) )  ; DVDL1MHxztL_F0_io = 0.0_WP
        ALLOCATE( T2xztL_F0_io( NCL2 ) ) ; T2xztL_F0_io = 0.0_WP
        ALLOCATE( D2xztL_F0_io( NCL2 ) ) ; D2xztL_F0_io = 0.0_WP
        ALLOCATE( H2xztL_F0_io( NCL2 ) ) ; H2xztL_F0_io = 0.0_WP

        ALLOCATE( DHxztL_F0_io( NCL2 ) ) ; DHxztL_F0_io = 0.0_WP
        ALLOCATE( PHxztL_F0_io( NCL2 ) ) ; PHxztL_F0_io = 0.0_WP

        !            ALLOCATE( DVDL1MxztL_F0_io( NCL2, NDV, NDV  )                     )
        !            ALLOCATE( DVDL1MHxztL_F0_io( NCL2, NDV, NDV  )                     )
        !            ALLOCATE( DVDL1MUxztL_F0_io( NCL2, NDV, NDV, NDV  )                )
        !            ALLOCATE( DVDL2MxztL_F0_io( NCL2, (NDV - 1) * 3 + NDV, (NDV - 1) * 3 + NDV  ) )

        ALLOCATE( UHxztL_F0_io( NCL2, NDV )                 ) ; UHxztL_F0_io = 0.0_WP
        ALLOCATE( GHxztL_F0_io( NCL2, NDV )                 ) ; GHxztL_F0_io = 0.0_WP
        ALLOCATE( U2DHxztL_F0_io( NCL2, NDV * (7 - NDV) / 2 + NDV - 3 ) ) ; U2DHxztL_F0_io = 0.0_WP

        ALLOCATE( DhDL1xztL_F0_io( NCL2, NDV )     ) ; DhDL1xztL_F0_io = 0.0_WP
        ALLOCATE( DhDLPxztL_F0_io( NCL2, NDV )     ) ; DhDLPxztL_F0_io = 0.0_WP
        ALLOCATE( DTDLKxztL_F0_io( NCL2, NDV )     ) ; DTDLKxztL_F0_io = 0.0_WP
        ALLOCATE( DTDLKUxztL_F0_io(NCL2, NDV, NDV ) ) ; DTDLKUxztL_F0_io = 0.0_WP

        ALLOCATE( DTDLKDVDLxztL_F0_io(NCL2, NDV, NDV, NDV ) ) ; DTDLKDVDLxztL_F0_io = 0.0_WP
        ALLOCATE( DHDLMDVDLxztL_F0_io(NCL2, NDV, NDV, NDV ) ) ; DHDLMDVDLxztL_F0_io = 0.0_WP


        ALLOCATE( CpSD  (NCL2)         ) ; CpSD = 0.0_WP
        ALLOCATE( QwSD  (NCL2)         ) ; QwSD = 0.0_WP
        ALLOCATE( TwSD  (NCL2)         ) ; TwSD = 0.0_WP
        ALLOCATE( HwSD  (NCL2)         ) ; HwSD = 0.0_WP

        !============================
        !ALLOCATE( Nuy(NCL2, 2))

        !========================================
        ALLOCATE( H_FA      (0 : NND2)         ) ; H_FA = 0.0_WP
        ALLOCATE( hff_RA    (NCL2)         ) ; hff_RA = 0.0_WP
        ALLOCATE( hfpf_RA (NCL2)         ) ; hfpf_RA = 0.0_WP
        ALLOCATE( dTDX      (NCL2, NDV)     ) ; dTDX = 0.0_WP
        ALLOCATE( dDDX      (NCL2, NDV)     ) ; dDDX = 0.0_WP
        ALLOCATE( dHDX_RA   (NCL2, NDV)     ) ; dHDX_RA = 0.0_WP
        ALLOCATE( dHDX_FA   (NCL2, NDV)     ) ; dHDX_FA = 0.0_WP
        ALLOCATE( UH_FA     (NCL2, NDV)     ) ; UH_FA = 0.0_WP
        ALLOCATE( uff2hffd_FA(NCL2, NDV, NDV) ) ; uff2hffd_FA = 0.0_WP

        ALLOCATE( uffhffd_FA        (NCL2, NDV)      ) ; uffhffd_FA = 0.0_WP
        ALLOCATE( ufhfd_RA        (NCL2, NDV)      ) ; ufhfd_RA = 0.0_WP
        ALLOCATE( viscstressEnth_RA    (NCL2, NDV, NDV) ) ; viscstressEnth_RA = 0.0_WP
        ALLOCATE( viscstressEnthGrad_RA(NCL2, NDV, NDV, NDV) ) ; viscstressEnthGrad_RA = 0.0_WP
        !========================================

        ALLOCATE(Budg_prodc_stres_thf(NCL2, NDV) ) ; Budg_prodc_stres_thF = 0.0_WP
        ALLOCATE(Budg_prodc_enthg_thf(NCL2, NDV) ) ; Budg_prodc_enthg_thF = 0.0_WP
        ALLOCATE(Budg_Turbu_diffu_thf(NCL2, NDV) ) ; Budg_Turbu_diffu_thF = 0.0_WP
        ALLOCATE(Budg_press_accl1_thf(NCL2, NDV) ) ; Budg_press_accl1_thF = 0.0_WP
        ALLOCATE(Budg_DphDX_diffu_thf(NCL2, NDV) ) ; Budg_DphDX_diffu_thF = 0.0_WP
        ALLOCATE(Budg_pdhDX_stran_thf(NCL2, NDV) ) ; Budg_pdhDX_stran_thF = 0.0_WP
        ALLOCATE(Budg_ConHF_accel_thf(NCL2, NDV) ) ; Budg_ConHF_accel_thF = 0.0_WP
        ALLOCATE(Budg_ConHF_diffu_thf(NCL2, NDV) ) ; Budg_ConHF_diffu_thF = 0.0_WP
        ALLOCATE(Budg_ConHF_dissp_thf(NCL2, NDV) ) ; Budg_ConHF_dissp_thF = 0.0_WP
        ALLOCATE(Budg_viscs_accl1_thf(NCL2, NDV) ) ; Budg_viscs_accl1_thF = 0.0_WP
        ALLOCATE(Budg_viscs_diffu_thf(NCL2, NDV) ) ; Budg_viscs_diffu_thF = 0.0_WP
        ALLOCATE(Budg_viscs_dissp_thf(NCL2, NDV) ) ; Budg_viscs_dissp_thF = 0.0_WP
        ALLOCATE(Budg_Balance1_thf    (NCL2, NDV) ) ; Budg_Balance1_thf  = 0.0_WP
        
        ALLOCATE(Budg_prodc_stres_IEN(NCL2) ) ; Budg_prodc_stres_IEN = 0.0_WP
        ALLOCATE(Budg_prodc_enthg_IEN(NCL2) ) ; Budg_prodc_enthg_IEN = 0.0_WP
        ALLOCATE(Budg_Turbu_diffu_IEN(NCL2) ) ; Budg_Turbu_diffu_IEN = 0.0_WP
        ALLOCATE(Budg_press_accl1_IEN(NCL2) ) ; Budg_press_accl1_IEN = 0.0_WP
        ALLOCATE(Budg_DphDX_diffu_IEN(NCL2) ) ; Budg_DphDX_diffu_IEN = 0.0_WP
        ALLOCATE(Budg_pdhDX_stran_IEN(NCL2) ) ; Budg_pdhDX_stran_IEN = 0.0_WP
        ALLOCATE(Budg_ConHF_accel_IEN(NCL2) ) ; Budg_ConHF_accel_IEN = 0.0_WP
        ALLOCATE(Budg_ConHF_diffu_IEN(NCL2) ) ; Budg_ConHF_diffu_IEN = 0.0_WP
        ALLOCATE(Budg_ConHF_dissp_IEN(NCL2) ) ; Budg_ConHF_dissp_IEN = 0.0_WP
        ALLOCATE(Budg_viscs_accl1_IEN(NCL2) ) ; Budg_viscs_accl1_IEN = 0.0_WP
        ALLOCATE(Budg_viscs_diffu_IEN(NCL2) ) ; Budg_viscs_diffu_IEN = 0.0_WP
        ALLOCATE(Budg_viscs_dissp_IEN(NCL2) ) ; Budg_viscs_dissp_IEN = 0.0_WP
        ALLOCATE(Budg_Balance1_IEN    (NCL2) ) ; Budg_Balance1_IEN  = 0.0_WP
!========================body forcE =================================
        ALLOCATE(Budg_prodc_gvfc2_thf(NCL2)     ) ; Budg_prodc_gvfc2_thF = 0.0_WP

    END IF

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE MEMO_DEALLT_AVERAGE_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE


    !========================================
    DEALLOCATE( U1xztL_F0_io )
    DEALLOCATE( G1xztL_F0_io )
    DEALLOCATE( UPxztL_F0_io )

    DEALLOCATE( U2xztL_F0_io  )
    DEALLOCATE( UGxztL_F0_io  )
    DEALLOCATE( UGUxztL_F0_io )
    DEALLOCATE( U3xztL_F0_io )

    DEALLOCATE( DVDL1xztL_F0_io )
    DEALLOCATE( DVDLPxztL_F0_io )
    DEALLOCATE( DVDL2xztL_F0_io )
    IF(iPPQuadrants == 1)  THEN
    DEALLOCATE( QuadUVxztL_F0_io )
    DEALLOCATE( QuadVzxztL_F0_io )
    DEALLOCATE( QuadTKxztL_F0_io )
    DEALLOCATE( QuadDRxztL_F0_io )
    DEALLOCATE( QuadDUV1xztL_F0_io )
    DEALLOCATE( QuadDUV2xztL_F0_io )

    DEALLOCATE( OctTUVxztL_F0_io )
    DEALLOCATE( OctTVzxztL_F0_io )
    DEALLOCATE( OctTTKxztL_F0_io )
    DEALLOCATE( OctTDRxztL_F0_io )
    DEALLOCATE( OctTDUV1xztL_F0_io )
    DEALLOCATE( OctTDUV2xztL_F0_io )

    DEALLOCATE( OctDUVxztL_F0_io )
    DEALLOCATE( OctDVzxztL_F0_io )
    DEALLOCATE( OctDTKxztL_F0_io )
    DEALLOCATE( OctDDRxztL_F0_io )
    DEALLOCATE( OctDDUV1xztL_F0_io )
    DEALLOCATE( OctDDUV2xztL_F0_io )
    END IF

    DEALLOCATE( FUxztL_F0_io )
    !========================================
    DEALLOCATE( dPDX_RA    )
    DEALLOCATE( ufpf_RA    )
    DEALLOCATE( uf2_RA     )
    DEALLOCATE( uf2d_RA    )
    DEALLOCATE( uf3_RA     )
    DEALLOCATE( uf3d_RA    )
    DEALLOCATE( UU_RA      )

    DEALLOCATE( dUiDXi     )
    DEALLOCATE( StrAInTensor)
    DEALLOCATE( VortcyTensor)
    DEALLOCATE( Skewness_RA)

    DEALLOCATE( MKE_RA     )
    DEALLOCATE( TKE_RA     )
    DEALLOCATE( ufTKEfd_RA )
    DEALLOCATE( ufMKEfd_RA )

    DEALLOCATE( Omega2_RA              )
    DEALLOCATE( Omega_RA2              )
    DEALLOCATE( Omega_rms              )

    DEALLOCATE( AnIsotropy_RA          )
    DEALLOCATE( Anistpinva_RA          )
    DEALLOCATE( LumleyAxis_RA          )

    !========================================
    IF(iThermoDynamics == 1) THEN
        DEALLOCATE( DrivenForce)
        DEALLOCATE( U_FA       )
        DEALLOCATE( UU_FA      )
        DEALLOCATE( dUDX_FA    )

        DEALLOCATE( uff_RA     )
        DEALLOCATE( uff2_FA    )
        DEALLOCATE( uff2d_FA   )
        DEALLOCATE( uff3_FA    )
        DEALLOCATE( uff3d_FA   )
        DEALLOCATE( TDIFU_FA   )

        DEALLOCATE( dUiDXiM      )
        DEALLOCATE( StrAInTensorM)
        DEALLOCATE( VortcyTensorM)
        DEALLOCATE( Skewness_FA  )

        DEALLOCATE( MKE_FA       )
        DEALLOCATE( TKE_FA       )
        DEALLOCATE( uffTKEffd_FA )
        DEALLOCATE( uffMKEffd_FA )

        DEALLOCATE( AnIsotropy_FA          )
        DEALLOCATE( Anistpinva_FA          )
        DEALLOCATE( LumleyAxis_FA          )
    END IF

    DEALLOCATE( Tau_Mean_RA  )
    DEALLOCATE( Tau_meaU_RA  )
    DEALLOCATE( dTaudy_RA    )
    DEALLOCATE( dTSSdy_RA    )
    DEALLOCATE( TauU_RA      )
    DEALLOCATE( Taufuf_RA    )
    DEALLOCATE( TauDvDL_RA   )
    !DEALLOCATE( Tau_ik_Du_jDX_i_RA)


    DEALLOCATE( NSFbal_FA)
    DEALLOCATE( NSFbal_RA)
    DEALLOCATE( ENEbal_FA)

    DEALLOCATE(Bdfcintg)
    DEALLOCATE( DensIntg)
    !========================================
    !================== RuV =============================
    DEALLOCATE (  Budg_prodc_stres_Duiuj )
    DEALLOCATE (  Budg_viscs_dissp_Duiuj )
    DEALLOCATE (  Budg_pduDX_stran_Duiuj )
    DEALLOCATE (  Budg_Turbu_diffu_Duiuj )
    DEALLOCATE (  Budg_DpuDX_diffu_Duiuj )
    DEALLOCATE (  Budg_viscs_diffu_Duiuj )

    DEALLOCATE (  Budg_press_accl1_Duiuj )
    DEALLOCATE (  Budg_viscs_accl1_Duiuj )
    DEALLOCATE (  Budg_prodc_Dvfc1_Duiuj )
    DEALLOCATE (  Budg_Balance1_Duiuj    )

    DEALLOCATE (  Budg_prodc_gvfc2_Duiuj )
    DEALLOCATE (  Budg_prodc_Dvfc2_Duiuj )
    DEALLOCATE (  Budg_turss_accl2_Duiuj )
    DEALLOCATE (  Budg_Balance2_Duiuj    )

    DEALLOCATE (  Budg_pressure3_Duiuj )
    DEALLOCATE (  Budg_vistress3_Duiuj )
    DEALLOCATE (  Budg_Balance3_Duiuj  )

    !================== RuV == Sum along Y =======================
    DEALLOCATE (  Budg_prodc_stres_Duiuj_ysum )
    DEALLOCATE (  Budg_viscs_dissp_Duiuj_ysum )
    DEALLOCATE (  Budg_pduDX_stran_Duiuj_ysum )
    DEALLOCATE (  Budg_Turbu_diffu_Duiuj_ysum )
    DEALLOCATE (  Budg_DpuDX_diffu_Duiuj_ysum )
    DEALLOCATE (  Budg_viscs_diffu_Duiuj_ysum )

    DEALLOCATE (  Budg_press_accl1_Duiuj_ysum )
    DEALLOCATE (  Budg_viscs_accl1_Duiuj_ysum )
    DEALLOCATE (  Budg_prodc_Dvfc1_Duiuj_ysum )
    DEALLOCATE (  Budg_Balance1_Duiuj_ysum    )

    DEALLOCATE (  Budg_prodc_gvfc2_Duiuj_ysum )
    DEALLOCATE (  Budg_prodc_Dvfc2_Duiuj_ysum )
    DEALLOCATE (  Budg_turss_accl2_Duiuj_ysum )
    DEALLOCATE (  Budg_Balance2_Duiuj_ysum    )

    DEALLOCATE (  Budg_pressure3_Duiuj_ysum )
    DEALLOCATE (  Budg_vistress3_Duiuj_ysum )
    DEALLOCATE (  Budg_Balance3_Duiuj_ysum  )

    !================== TKE =============================
    DEALLOCATE (  Budg_prodc_stres_TKE )
    DEALLOCATE (  Budg_viscs_dissp_TKE )
    DEALLOCATE (  Budg_pduDX_stran_TKE )
    DEALLOCATE (  Budg_Turbu_diffu_TKE )
    DEALLOCATE (  Budg_DpuDX_diffu_TKE )
    DEALLOCATE (  Budg_viscs_diffu_TKE )

    DEALLOCATE (  Budg_press_accl1_TKE )
    DEALLOCATE (  Budg_viscs_accl1_TKE )
    DEALLOCATE (  Budg_prodc_Dvfc1_TKE )
    DEALLOCATE (  Budg_Balance1_TKE    )

    DEALLOCATE (  Budg_prodc_gvfc2_TKE )
    DEALLOCATE (  Budg_prodc_Dvfc2_TKE )
    DEALLOCATE (  Budg_turss_accl2_TKE )
    DEALLOCATE (  Budg_Balance2_TKE    )

    DEALLOCATE (  Budg_pressure3_TKE )
    DEALLOCATE (  Budg_vistress3_TKE )
    DEALLOCATE (  Budg_Balance3_TKE  )

    !================== MKE =============================
    DEALLOCATE (  Budg_prodc_stres_MKE )
    DEALLOCATE (  Budg_viscs_dissp_MKE )
    DEALLOCATE (  Budg_pduDX_stran_MKE )
    DEALLOCATE (  Budg_Turbu_diffu_MKE )
    DEALLOCATE (  Budg_DpuDX_diffu_MKE )
    DEALLOCATE (  Budg_viscs_diffu_MKE )

    DEALLOCATE (  Budg_press_accl1_MKE )
    DEALLOCATE (  Budg_viscs_accl1_MKE )
    DEALLOCATE (  Budg_prodc_Dvfc1_MKE )
    DEALLOCATE (  Budg_prodc_gvfc1_MKE )
    DEALLOCATE (  Budg_Balance1_MKE    )

    DEALLOCATE (  Budg_prodc_gvfc2_MKE )
    DEALLOCATE (  Budg_prodc_Dvfc2_MKE )
    DEALLOCATE (  Budg_turss_accl2_MKE )
    DEALLOCATE (  Budg_Balance2_MKE    )

    DEALLOCATE (  Budg_pressure3_MKE )
    DEALLOCATE (  Budg_vistress3_MKE )
    DEALLOCATE (  Budg_Balance3_MKE  )


    !======================================
    DEALLOCATE( RANS_Mut)

    DEALLOCATE( YWdISD )
    DEALLOCATE( TauwSD )
    DEALLOCATE( DensSD )
    DEALLOCATE( viscsD )


    DEALLOCATE( D1xztL_F0_io )
    DEALLOCATE( M1xztL_F0_io )

    DEALLOCATE( DVDL1MxztL_F0_io  )
    DEALLOCATE( DVDL1MUxztL_F0_io )
    DEALLOCATE( DVDL2MxztL_F0_io  )

    IF(iThermoDynamics == 1) THEN
        !========================================
        DEALLOCATE( T1xztL_F0_io )
        DEALLOCATE( H1xztL_F0_io )

        DEALLOCATE( DVDL1MHxztL_F0_io )
        DEALLOCATE( T2xztL_F0_io )
        DEALLOCATE( D2xztL_F0_io )
        DEALLOCATE( H2xztL_F0_io )

        DEALLOCATE( DHxztL_F0_io )
        DEALLOCATE( PHxztL_F0_io )

        DEALLOCATE( UHxztL_F0_io   )
        DEALLOCATE( GHxztL_F0_io   )
        DEALLOCATE( U2DHxztL_F0_io )

        DEALLOCATE( DhDL1xztL_F0_io      )
        DEALLOCATE( DhDLPxztL_F0_io      )
        DEALLOCATE( DTDLKxztL_F0_io      )
        DEALLOCATE( DTDLKUxztL_F0_io     )

        DEALLOCATE( DTDLKDVDLxztL_F0_io )
        DEALLOCATE( DHDLMDVDLxztL_F0_io )

        !==============================
        !DEALLOCATE( Nuy  )
        !==================

        DEALLOCATE( CpSD )
        DEALLOCATE( QwSD )
        DEALLOCATE( TwSD )
        DEALLOCATE( HwSD )

        !========================================
        DEALLOCATE( H_FA       )
        DEALLOCATE( hff_RA     )
        DEALLOCATE( hfpf_RA  )
        DEALLOCATE( dTDX       )
        DEALLOCATE( dDDX       )
        DEALLOCATE( dHDX_RA    )
        DEALLOCATE( dHDX_FA    )
        DEALLOCATE( UH_FA      )
        DEALLOCATE( uff2hffd_FA )

        DEALLOCATE( uffhffd_FA         )
        DEALLOCATE( ufhfd_RA         )
        DEALLOCATE( viscstressEnth_RA     )
        DEALLOCATE( viscstressEnthGrad_RA )
        !========================================

        DEALLOCATE(Budg_prodc_stres_thf )
        DEALLOCATE(Budg_prodc_enthg_thf )
        DEALLOCATE(Budg_Turbu_diffu_thf )
        DEALLOCATE(Budg_press_accl1_thf )
        DEALLOCATE(Budg_DphDX_diffu_thf )
        DEALLOCATE(Budg_pdhDX_stran_thf )
        DEALLOCATE(Budg_ConHF_accel_thf )
        DEALLOCATE(Budg_ConHF_diffu_thf )
        DEALLOCATE(Budg_ConHF_dissp_thf )
        DEALLOCATE(Budg_viscs_accl1_thf )
        DEALLOCATE(Budg_viscs_diffu_thf )
        DEALLOCATE(Budg_viscs_dissp_thf )
        DEALLOCATE(Budg_Balance1_thf    )

        DEALLOCATE(Budg_prodc_stres_IEN )
        DEALLOCATE(Budg_prodc_enthg_IEN )
        DEALLOCATE(Budg_Turbu_diffu_IEN )
        DEALLOCATE(Budg_press_accl1_IEN )
        DEALLOCATE(Budg_DphDX_diffu_IEN )
        DEALLOCATE(Budg_pdhDX_stran_IEN )
        DEALLOCATE(Budg_ConHF_accel_IEN )
        DEALLOCATE(Budg_ConHF_diffu_IEN )
        DEALLOCATE(Budg_ConHF_dissp_IEN )
        DEALLOCATE(Budg_viscs_accl1_IEN )
        DEALLOCATE(Budg_viscs_diffu_IEN )
        DEALLOCATE(Budg_viscs_dissp_IEN )
        DEALLOCATE(Budg_Balance1_IEN    )
        !========================body forcE =================================

        DEALLOCATE(Budg_prodc_gvfc2_thf   )

    END IF

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_Xperiodic_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    IF(MYID == 0) CALL CHKRLHDL('15.IO: postprocessing general data at', MYID, PhyTIME_io)

    !==========Gather Dwta================
    CALL WRT_AVERAGE_PPED_XZ_io_GATHER

    !==========Calcuate wall info================
    CALL PP_wall_thermal_shear(flgxzt)

    IF(MYID /= 0) RETURN

    !=========Calculate data=======================
    CALL PP_FLOW_BASIC_VARS_XZ_io
    CALL PP_FLOW_FA_RSTE_Budg_XZ_io
    IF(iThermoDynamics == 1) CALL PP_SSzero_SIDED(1)
    !==========WRITE out wall info=================
    IF(iThermoDynamics == 1) THEN
        CALL WRT_HeatTransfer_Table_XZ_io
    END IF
    CALL WRT_Cf_Table_XZ_io ! must be after heat
    !CALL PP_Umax_SIDED
    CALL WRT_FLOW_FA_Profile_XZ_io
    IF(iThermoDynamics == 1) THEN
        CALL PP_HEAT_BASIC_VARS_XZ_io
        CALL PP_HEAT_FA_RSTE_Budg_XZ_io
        CALL WRT_HEAT_FA_Profile_XZ_io
    END IF

    CALL PP_FLOW_RA_noDen_RSTE_Budg_XZ_io
    IF(iThermoDynamics  /= 1) CALL PP_SSzero_SIDED(2)
    CALL WRT_FLOW_RA_Profile_XZ_io

    CALL WRT_Checking_TABLE_XZ_io

    IF(iPPSpectra == 1) THEN
        CALL WRITE_SPECO_AVE_PROFILE('FLOW')
        CALL WRITE_SPECO_AVE_Contour('FLOW')
        IF(iThermoDynamics == 1) THEN
            CALL WRITE_SPECO_AVE_PROFILE('HEAT')
            CALL WRITE_SPECO_AVE_Contour('HEAT')
        END IF
    END IF
    !CALL MEMO_DEALLT_INTP_XZ_io
    CALL MEMO_DEALLT_AVERAGE_XZ_io

    IF(MYID == 0) CALL CHKRLHDL('IO: finished postprocessing general data at', MYID, PhyTIME_io)

END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_AVERAGE_PPED_XZ_io_GATHER
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: INN
    INTEGER(4) :: J, JJ, I
    INTEGER(4) :: L, IP, M, H, N
    INTEGER(4) :: N2DOID

    CHARACTER(15) :: NXSTR
    REAL(WP) :: urms, vrms, WRms, uv, uw, vw
    REAL(WP) ::  COE, Pwall
    INTEGER(4) :: TECFLG1, TECFLG2, TECFLG3

    REAL(WP) :: FU_AUX (N2DO(MYID), NDV + 1,                               1 : NPTOT)
    REAL(WP) :: U1_AUX (N2DO(MYID), NDV + 1,                               1 : NPTOT)
    REAL(WP) :: G1_AUX (N2DO(MYID), NDV,                                 1 : NPTOT)
    REAL(WP) :: UP_AUX (N2DO(MYID), NDV,                                 1 : NPTOT)

    REAL(WP) :: U2_AUX (N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3),                1 : NPTOT)
    REAL(WP) :: UG_AUX (N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3),                1 : NPTOT)
    REAL(WP) :: UGU_AUX(N2DO(MYID), (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8),  1 : NPTOT)
    REAL(WP) :: U3_AUX (N2DO(MYID), (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8),  1 : NPTOT)

    REAL(WP) :: DVDL1_AUX (N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: DVDLP_AUX (N2DO(MYID), NDV,                 NDV, 1 : NPTOT)
    REAL(WP) :: DVDL2_AUX (N2DO(MYID), (NDV - 1) * 3 + NDV, (NDV - 1) * 3 + NDV, 1 : NPTOT)

    REAL(WP) :: QuadUV_AUX (N2DO(MYID),   4, QUADHN, 1 : NPTOT)
    REAL(WP) :: QuadVz_AUX (N2DO(MYID),   4, QUADHN, 1 : NPTOT)
    REAL(WP) :: QuadTK_AUX (N2DO(MYID),   4, QUADHN, 1 : NPTOT)
    REAL(WP) :: QuadDR_AUX (N2DO(MYID),   4, QUADHN, 1 : NPTOT)
    REAL(WP) :: QuadDUV1_AUX (N2DO(MYID), 4, QUADHN, 1 : NPTOT)
    REAL(WP) :: QuadDUV2_AUX (N2DO(MYID), 4, QUADHN, 1 : NPTOT)

    REAL(WP) :: OctDUV_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctDVz_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctDTK_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctDDR_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctDDUV1_AUX (N2DO(MYID), 8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctDDUV2_AUX (N2DO(MYID), 8, QUADHN, 1 : NPTOT)

    REAL(WP) :: OctTUV_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctTVz_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctTTK_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctTDR_AUX (N2DO(MYID),   8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctTDUV1_AUX (N2DO(MYID), 8, QUADHN, 1 : NPTOT)
    REAL(WP) :: OctTDUV2_AUX (N2DO(MYID), 8, QUADHN, 1 : NPTOT)

    REAL(WP) :: T1_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: D1_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: H1_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: M1_AUX (N2DO(MYID), 1 : NPTOT)

    REAL(WP) :: T2_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: D2_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: H2_AUX (N2DO(MYID), 1 : NPTOT)

    REAL(WP) :: DH_AUX (N2DO(MYID), 1 : NPTOT)
    REAL(WP) :: PH_AUX (N2DO(MYID), 1 : NPTOT)

    REAL(WP) :: DVDL1M_AUX  (N2DO(MYID), NDV,                 NDV,     1 : NPTOT)
    REAL(WP) :: DVDL1MH_AUX (N2DO(MYID), NDV,                 NDV,     1 : NPTOT)
    REAL(WP) :: DVDL1MU_AUX (N2DO(MYID), NDV, NDV, NDV,                1 : NPTOT)
    REAL(WP) :: DVDL2M_AUX  (N2DO(MYID), (NDV - 1) * 3 + NDV, (NDV - 1) * 3 + NDV, 1 : NPTOT)

    REAL(WP) :: UH_AUX   (N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: GH_AUX   (N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: U2DH_AUX (N2DO(MYID), (NDV * (7 - NDV) / 2 + NDV - 3), 1 : NPTOT)

    REAL(WP) :: DhDL1_AUX (N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: DhDLP_AUX (N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: DTDLK_AUX     (N2DO(MYID), NDV, 1 : NPTOT)
    REAL(WP) :: DTDLKU_AUX    (N2DO(MYID), NDV, NDV, 1 : NPTOT)
    REAL(WP) :: DTDLKDVDL_AUX (N2DO(MYID), NDV, NDV, NDV, 1 : NPTOT)
    REAL(WP) :: DHDLMDVDL_AUX (N2DO(MYID), NDV, NDV, NDV, 1 : NPTOT)

    !======================================
    INN = N2DO(MYID) * (NDV + 1)
    CALL MPI_GATHER( U1xztL_io, INN, MPI_DOUBLE_PRECISION, U1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * NDV
    CALL MPI_GATHER( G1xztL_io, INN, MPI_DOUBLE_PRECISION, G1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * NDV
    CALL MPI_GATHER( UPxztL_io, INN, MPI_DOUBLE_PRECISION, UP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !======================================
    INN = N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
    CALL MPI_GATHER( U2xztL_io, INN, MPI_DOUBLE_PRECISION, U2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
    CALL MPI_GATHER( UGxztL_io, INN, MPI_DOUBLE_PRECISION, UG_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
    CALL MPI_GATHER( UGUxztL_io, INN, MPI_DOUBLE_PRECISION, UGU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
    CALL MPI_GATHER( U3xztL_io, INN, MPI_DOUBLE_PRECISION, U3_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !======================================
    INN = N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDL1xztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * NDV * NDV
    CALL MPI_GATHER( DVDLPxztL_io, INN, MPI_DOUBLE_PRECISION, DVDLP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    INN = N2DO(MYID) * ((NDV - 1) * 3 + NDV) * ((NDV - 1) * 3 + NDV)
    CALL MPI_GATHER( DVDL2xztL_io, INN, MPI_DOUBLE_PRECISION, DVDL2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)
    IF(iPPQuadrants == 1)  THEN
        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadUVxztL_io, INN, MPI_DOUBLE_PRECISION, QuadUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadVzxztL_io, INN, MPI_DOUBLE_PRECISION, QuadVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadTKxztL_io, INN, MPI_DOUBLE_PRECISION, QuadTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadDRxztL_io, INN, MPI_DOUBLE_PRECISION, QuadDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, QuadDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 4 * QUADHN
        CALL MPI_GATHER( QuadDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, QuadDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)
    END IF
    !=======================================================
    INN = N2DO(MYID) * (NDV + 1)
    CALL MPI_GATHER( FUxztL_io, INN, MPI_DOUBLE_PRECISION, FU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
    CALL MPI_BARRIER(ICOMM, IERROR)

    !=======================================================
    IF(iPPQuadrants == 1)  THEN
        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDUVxztL_io, INN, MPI_DOUBLE_PRECISION, OctDUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDVzxztL_io, INN, MPI_DOUBLE_PRECISION, OctDVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDTKxztL_io, INN, MPI_DOUBLE_PRECISION, OctDTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDDRxztL_io, INN, MPI_DOUBLE_PRECISION, OctDDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, OctDDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctDDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, OctDDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !=======================================================
        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTUVxztL_io, INN, MPI_DOUBLE_PRECISION, OctTUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTVzxztL_io, INN, MPI_DOUBLE_PRECISION, OctTVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTTKxztL_io, INN, MPI_DOUBLE_PRECISION, OctTTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTDRxztL_io, INN, MPI_DOUBLE_PRECISION, OctTDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, OctTDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * 8 * QUADHN
        CALL MPI_GATHER( OctTDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, OctTDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)
    END IF

    IF(iThermoDynamics == 1) THEN
        !===============================
        INN = N2DO(MYID)
        CALL MPI_GATHER( T1xztL_io, INN, MPI_DOUBLE_PRECISION, T1_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( D1xztL_io, INN, MPI_DOUBLE_PRECISION, D1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( H1xztL_io, INN, MPI_DOUBLE_PRECISION, H1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( M1xztL_io, INN, MPI_DOUBLE_PRECISION, M1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !===============================
        INN = N2DO(MYID)
        CALL MPI_GATHER( T2xztL_io, INN, MPI_DOUBLE_PRECISION, T2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( D2xztL_io, INN, MPI_DOUBLE_PRECISION, D2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( H2xztL_io, INN, MPI_DOUBLE_PRECISION, H2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !================================
        INN = N2DO(MYID)
        CALL MPI_GATHER( DHxztL_io, INN, MPI_DOUBLE_PRECISION, DH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID)
        CALL MPI_GATHER( PHxztL_io, INN, MPI_DOUBLE_PRECISION, PH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !================================
        INN = N2DO(MYID) * NDV * NDV
        CALL MPI_GATHER( DVDL1MxztL_io,  INN, MPI_DOUBLE_PRECISION, DVDL1M_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV * NDV
        CALL MPI_GATHER( DVDL1MHxztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1MH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV * NDV * NDV
        CALL MPI_GATHER( DVDL1MUxztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1MU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * ((NDV - 1) * 3 + NDV) * ((NDV - 1) * 3 + NDV)
        CALL MPI_GATHER( DVDL2MxztL_io,  INN, MPI_DOUBLE_PRECISION, DVDL2M_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !            !====== TesT =============
        !            IF (MYID == 0) THEN
        !            DO J = 1, N2DO(MYID)
        !                DO L = 1, (NDV - 1) * 3 + NDV
        !                    DO M = 1, (NDV - 1) * 3 + NDV
        !                        !WRITE(*, *) L, M, DVDL2xztL_io(J, L, M) * M1xztL_io(J), DVDL2MxztL_io(J, L, M)
        !                    END DO
        !                END DO
        !            END DO
        !            END IF
        !            !====== TesT =============

        !============================
        INN = N2DO(MYID) * NDV
        CALL MPI_GATHER( UHxztL_io, INN, MPI_DOUBLE_PRECISION, UH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV
        CALL MPI_GATHER( GHxztL_io, INN, MPI_DOUBLE_PRECISION, GH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * (NDV * (7 - NDV) / 2 + NDV - 3)
        CALL MPI_GATHER( U2DHxztL_io, INN, MPI_DOUBLE_PRECISION, U2DH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        !============================
        INN = N2DO(MYID) * NDV
        CALL MPI_GATHER( DhDL1xztL_io, INN, MPI_DOUBLE_PRECISION, DhDL1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV
        CALL MPI_GATHER( DhDLPxztL_io, INN, MPI_DOUBLE_PRECISION, DhDLP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV
        CALL MPI_GATHER( DTDLKxztL_io, INN, MPI_DOUBLE_PRECISION, DTDLK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV * NDV
        CALL MPI_GATHER( DTDLKUxztL_io, INN, MPI_DOUBLE_PRECISION, DTDLKU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV * NDV * NDV
        CALL MPI_GATHER( DTDLKDVDLxztL_io, INN, MPI_DOUBLE_PRECISION, DTDLKDVDL_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

        INN = N2DO(MYID) * NDV * NDV * NDV
        CALL MPI_GATHER( DHDLMDVDLxztL_io, INN, MPI_DOUBLE_PRECISION, DHDLMDVDL_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM, IERROR)

    END IF


    IF(MYID == 0) THEN
        !==================================================================
        CALL MEMO_ALLOCT_AVERAGE_XZ_io

        DO IP = 0, NPSLV
            N2DOID=JDEWT(IP) - JDSWT(IP) + 1
            DO J = 1, N2DOID
                JJ = JDSWT(IP) - 1 + J

                !===============================
                DO L = 1, NDV + 1
                    U1xztL_F0_io(JJ, L) = U1_AUX(J, L, IP + 1)
                END DO

                DO L = 1, NDV
                    G1xztL_F0_io(JJ, L) = G1_AUX(J, L, IP + 1)
                END DO

                DO L = 1, NDV
                    UPxztL_F0_io(JJ, L) = UP_AUX(J, L, IP + 1)
                END DO

                !===============================
                DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                    U2xztL_F0_io(JJ, L) = U2_AUX(J, L, IP + 1)
                END DO

                DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                    UGxztL_F0_io(JJ, L) = UG_AUX(J, L, IP + 1)
                END DO

                DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
                    UGUxztL_F0_io(JJ, L) = UGU_AUX(J, L, IP + 1)
                END DO

                DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
                    U3xztL_F0_io(JJ, L) = U3_AUX(J, L, IP + 1)
                END DO

                !===============================
                DO L = 1, NDV
                    DO M = 1, NDV
                        DVDL1xztL_F0_io(JJ, L, M) = DVDL1_AUX(J, L, M, IP + 1)
                    END DO
                END DO

                DO L = 1, NDV
                    DO M = 1, NDV
                        DVDLPxztL_F0_io(JJ, L, M) = DVDLP_AUX(J, L, M, IP + 1)
                    END DO
                END DO

                DO L = 1, (NDV - 1) * 3 + NDV
                    DO M = 1, (NDV - 1) * 3 + NDV
                        DVDL2xztL_F0_io(JJ, L, M) = DVDL2_AUX(J, L, M, IP + 1)
                    END DO
                END DO

                !========= QuadranT ======================
                IF(iPPQuadrants == 1)  THEN
                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadUVxztL_F0_io(JJ, L, M) = QuadUV_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadVzxztL_F0_io(JJ, L, M) = DSQRT(QuadVz_AUX(J, L, M, IP + 1))
                        END DO
                    END DO

                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadTKxztL_F0_io(JJ, L, M) = QuadTK_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadDRxztL_F0_io(JJ, L, M) = QuadDR_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadDUV1xztL_F0_io(JJ, L, M) = QuadDUV1_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, 4
                        DO M = 1, QUADHN
                            QuadDUV2xztL_F0_io(JJ, L, M) = QuadDUV2_AUX(J, L, M, IP + 1)
                        END DO
                    END DO
                END IF
                !================================
                DO L = 1, NDV + 1
                    FUxztL_F0_io(JJ, L) = FU_AUX(J, L, IP + 1)
                END DO

                !=========OctanT ==\rho==================
                IF(iPPQuadrants == 1)  THEN
                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDUVxztL_F0_io(JJ, L, M) = OctDUV_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDVzxztL_F0_io(JJ, L, M) = DSQRT(OctDVz_AUX(J, L, M, IP + 1))
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDTKxztL_F0_io(JJ, L, M) = OctDTK_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDDRxztL_F0_io(JJ, L, M) = OctDDR_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDDUV1xztL_F0_io(JJ, L, M) = OctDDUV1_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctDDUV2xztL_F0_io(JJ, L, M) = OctDDUV2_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    !=========OctanT == T ===================
                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTUVxztL_F0_io(JJ, L, M) = OctTUV_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTVzxztL_F0_io(JJ, L, M) = DSQRT(OctTVz_AUX(J, L, M, IP + 1))
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTTKxztL_F0_io(JJ, L, M) = OctTTK_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTDRxztL_F0_io(JJ, L, M) = OctTDR_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTDUV1xztL_F0_io(JJ, L, M) = OctTDUV1_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1,8
                        DO M = 1, QUADHN
                            OctTDUV2xztL_F0_io(JJ, L, M) = OctTDUV2_AUX(J, L, M, IP + 1)
                        END DO
                    END DO
                END IF
                !============================
                IF(iThermoDynamics == 1) THEN
                    !===============================
                    T1xztL_F0_io(JJ) = T1_AUX(J, IP + 1)
                    D1xztL_F0_io(JJ) = D1_AUX(J, IP + 1)
                    H1xztL_F0_io(JJ) = H1_AUX(J, IP + 1)
                    M1xztL_F0_io(JJ) = M1_AUX(J, IP + 1)

                    !===============================
                    T2xztL_F0_io(JJ) = T2_AUX(J, IP + 1)
                    D2xztL_F0_io(JJ) = D2_AUX(J, IP + 1)
                    H2xztL_F0_io(JJ) = H2_AUX(J, IP + 1)

                    !===============================
                    DHxztL_F0_io(JJ) = DH_AUX(J, IP + 1)
                    PHxztL_F0_io(JJ) = PH_AUX(J, IP + 1)
                    !===============================

                    DO L = 1, NDV
                        DO M = 1, NDV
                            DVDL1MxztL_F0_io(JJ, L, M) = DVDL1M_AUX (J, L, M, IP + 1)
                            DVDL1MHxztL_F0_io(JJ, L, M) = DVDL1MH_AUX(J, L, M, IP + 1)
                            DO N = 1, NDV
                                DVDL1MUxztL_F0_io(JJ, L, M, N) = DVDL1MU_AUX(J, L, M, N, IP + 1)
                            END DO
                        END DO
                    END DO

                    DO L = 1, (NDV - 1) * 3 + NDV
                        DO M = 1, (NDV - 1) * 3 + NDV
                            DVDL2MxztL_F0_io(JJ, L, M) = DVDL2M_AUX(J, L, M, IP + 1)
                        END DO
                    END DO

                    DO L = 1, NDV
                        UHxztL_F0_io(JJ, L) = UH_AUX(J, L, IP + 1)
                        GHxztL_F0_io(JJ, L) =GH_AUX(J, L, IP + 1)
                    END DO

                    DO L = 1, (NDV * (7 - NDV) / 2 + NDV - 3)
                        U2DHxztL_F0_io(JJ, L) = U2DH_AUX(J, L, IP + 1)
                    END DO

                    DO L = 1, NDV
                        DhDL1xztL_F0_io(JJ, L) = DhDL1_AUX(J, L, IP + 1)
                        DhDLPxztL_F0_io(JJ, L) = DhDLP_AUX(J, L, IP + 1)
                        DTDLKxztL_F0_io(JJ, L) = DTDLK_AUX(J, L, IP + 1)
                        DO M = 1, NDV
                            DTDLKUxztL_F0_io(JJ, L, M) = DTDLKU_AUX(J, L, M, IP + 1)
                            DO N = 1, NDV
                                DTDLKDVDLxztL_F0_io(JJ, L, M, N) = DTDLKDVDL_AUX(J, L, M, N, IP + 1)
                                DHDLMDVDLxztL_F0_io(JJ, L, M, N) = DHDLMDVDL_AUX(J, L, M, N, IP + 1)
                            END DO
                        END DO
                    END DO

                ELSE
                    ! below IS to test FA recovering to RA...
                    CALL CHKHDL('NOTICE: ThIS IS only for Isothermal flow!', MYID)
                    D1xztL_F0_io = 1.0_WP
                    M1xztL_F0_io = 1.0_WP
                    DO L = 1, NDV
                        DO M = 1, NDV
                            DVDL1MxztL_F0_io(JJ, L, M) = DVDL1xztL_F0_io(JJ, L, M)
                            DO N = 1, NDV
                                DVDL1MUxztL_F0_io(JJ, L, M, N) = 0.0_WP ! to get !!!?
                            END DO
                        END DO
                    END DO


                    DO L = 1, (NDV - 1) * 3 + NDV
                        DO M = 1, (NDV - 1) * 3 + NDV
                            DVDL2MxztL_F0_io(JJ, L, M) = DVDL2xztL_F0_io(JJ, L, M)
                        END DO
                    END DO


                END IF
            END DO
        END DO

        !================check data=============================
        !            DO L = 1, (NDV * (6 - NDV) + (NDV * (7 - NDV)) / 2 + NDV - 8)
        !            DO J = 1, NCL2
        !                !WRITE(*, *) L, J, UGUxztL_F0_io(J, L) / D1xztL_F0_io(J), U3xztL_F0_io(J, L), &
        !                                UGUxztL_F0_io(J, L) / D1xztL_F0_io(J) - U3xztL_F0_io(J, L) !test
        !            END DO
        !            END DO

        !================= Adjust any pressure term to make wall pressure zero=======================
        PwalL = U1xztL_F0_io(1, 4)
        DO J = 1, NCL2

            U1xztL_F0_io(J, 4) = U1xztL_F0_io(J, 4)  - Pwall

            DO M = 1, NDV
                UPxztL_F0_io(J, M) = UPxztL_F0_io(J, M) - PwalL * U1xztL_F0_io(J, M)
            END DO

            DO L = 1, NDV
                DO M = 1, NDV
                    DVDLPxztL_F0_io(J, L, M) = DVDLPxztL_F0_io(J, L, M) -PwalL * DVDL1xztL_F0_io(J, L, M)
                END DO
            END DO

            IF(iThermoDynamics == 1) THEN
                DO L = 1, NDV
                    DhDLPxztL_F0_io(J, L) = DhDLPxztL_F0_io(J, L) -PwalL * DhDL1xztL_F0_io(J, L)
                END DO
            END IF
            !!WRITE(*, *) YCC(J), U1xztL_F0_io(J, 4) !test
        END DO
        !================== Adjust any pressure term to make wall pressure zero==============
        !            CALL CHKHDL(' ==> In WRT..DVDL2 VS DVDL2M.', MYID)
        !            DO J = 1, NCL2
        !                DO L = 1, (NDV - 1) * 3 + NDV
        !                    DO M = 1, (NDV - 1) * 3 + NDV
        !                        !WRITE(*, *) J, L, M, DVDL2xztL_F0_io(J, L, M), DVDL2MxztL_F0_io(J, L, M), &
        !                                            DVDL2MxztL_F0_io(J, L, M) / DVDL2xztL_F0_io(J, L, M)
        !                    END DO
        !                END DO
        !            END DO


    END IF

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_FLOW_BASIC_VARS_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, H, P, LMNH, LMH, LMN, LNH
    REAL(WP) :: Matrix(3, 3)
    REAL(WP) :: EIG(3)
    REAL(WP) :: FF, TT
    REAL(WP) :: DrivenFCTT1, DrivenFCTT2
    REAL(WP) :: DrivenFCTTU1, DrivenFCTTU2
    INTEGER(4) :: TECflg = 200
    REAL(WP) :: DENtemp, ddenintg, DeNMintg

    !====================================================================
    !======= INTroductioN =================================================
    !====================================================================
    ! {*} = FA = Favred Averaged Mean.
    ! ''  = ff = Favred Averaged Fluctuation
    ! <>  = RA = Reynolds Averaged Mean
    ! ' = f  = Reynolds Averaged Fluctuation.

    !====================================================================
    !====================== RA based======================================
    !====================================================================

    CALL CHKHDL('   Calculating basic variables', MYID)
    !============== DPDX_RA(CL, M) = d<p>/ DX_m ========================
    dPDX_RA(:, :) = 0.0_WP
    N = 2
    DO J = 1, NCL2
        IF(J == 1) THEN
            dPDX_RA(J, N) = ( ( YCL2ND_WFB(J + 1) * U1xztL_F0_io(J,  4) + &
            YCL2ND_WFF(J + 1) * U1xztL_F0_io(J + 1, 4) ) - U1xztL_F0_io(J, 4)  ) * DYFI(J)

        ELSE IF (J == NCL2) THEN
            dPDX_RA(J, N) = (  U1xztL_F0_io(J, 4) - &
            ( YCL2ND_WFF(J) * U1xztL_F0_io(J,  4) + &
            YCL2ND_WFB(J) * U1xztL_F0_io(J - 1, 4) )  ) * DYFI(J)

        ELSE
            dPDX_RA(J, N) = ( ( YCL2ND_WFB(J + 1) * U1xztL_F0_io(J,  4) + &
            YCL2ND_WFF(J + 1) * U1xztL_F0_io(J + 1, 4) ) - &
            ( YCL2ND_WFF(J) * U1xztL_F0_io(J,  4) + &
            YCL2ND_WFB(J) * U1xztL_F0_io(J - 1, 4) ) ) * DYFI(J)
        END IF
    END DO
    CALL CHKHDL(' ==>Calculated d<P>/ Dy', MYID)

    !====<p' u'_i>=<p' u"_i>=<p u_i> - <p>*<u_i>==ufpf_RA(J, I) =============
    DO J = 1, NCL2
        DO M = 1, NDV
            ufpf_RA(J, M) =  UPxztL_F0_io(J, M) - U1xztL_F0_io(J, 4) * U1xztL_F0_io(J, M)
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated <u`_i p`>', MYID)

    !=======================================
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                LMN = (M * (7-M)) / 2 + N - 3
                UU_RA(J, M, N) = U2xztL_F0_io(J, LMN)
            END DO
        END DO

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) THEN
                    UU_RA(J, M, N) = UU_RA(J, N, M)
                END IF
            END DO
        END DO
        !WRITE(*, '(4ES13.5)') YCC(J), UU_FA(J, 1:3, 1:3) !test
    END DO
    CALL CHKHDL(' ==>Calculated <UU>', MYID)

    !==============<u'_i u'_j> ========================
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                !LMN = (M * (7-M)) / 2 + N - 3
                uf2_RA(J, M, N) = UU_RA(J, M, N) - U1xztL_F0_io(J, M) * U1xztL_F0_io(J, N)
            END DO
        END DO

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) THEN
                    uf2_RA(J, M, N) = Uf2_RA(J, N, M)
                END IF
            END DO
        END DO
        !WRITE(*, '(10ES13.5)') YCC(J), uf2_RA(J, 1, 1), uf2_RA(J, 2, 2),uf2_RA(J, 3, 3), &
        !            uf2_RA(J, 1, 2), uf2_RA(J, 2, 1),  &
        !            uf2_RA(J, 1, 3), uf2_RA(J, 3, 1),  &
        !            uf2_RA(J, 2, 3), uf2_RA(J, 3, 2)  ! test
    END DO
    CALL CHKHDL(' ==>Calculated <u`_i u`_j>', MYID)

    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                uf2d_RA(J, M, N) = uf2_RA(J, M, N) * D1xztL_F0_io(J)
            END DO
            !WRITE(*, '(10ES13.5)') YCC(J), uf2d_RA(J, 1, 1), uf2d_RA(J, 2, 2),uf2d_RA(J, 3, 3), &
            !        uf2d_RA(J, 1, 2), uf2d_RA(J, 2, 1),  &
            !        uf2d_RA(J, 1, 3), uf2d_RA(J, 3, 1),  &
            !        uf2d_RA(J, 2, 3), uf2d_RA(J, 3, 2)  ! test
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated <\rho>*<u`_i u`_j>', MYID)

    !=====u'u'u', u'u'v', u'u'w', u'v'v', u'v'w', u'w'w', v'v'v', v'v'w', v'w'w', w'w'w ===
    !{111} {112} {113}
    ![121] {122} {123}
    ![131] [132] {133}
    ![211] [212] [213]
    ![221] {222} {223}
    ![231] [232] {233}
    ![311] [312] [313]
    ![321] [322] [323]
    ![331] [332] {333}
    !        ! Below IS for test
    !        DO J = 1, NCL2
    !            DO M = 1, NDV
    !                DO N = 1, NDV
    !                    IF(M >  N) CYCLE
    !                    DO H = 1, NDV
    !                        IF(N >  H) CYCLE
    !                        LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
    !                        LMN = (M * (7-M)) / 2 + N - 3
    !                        LMH = (M * (7-M)) / 2 + H- 3
    !                        LNH = (N * (7-N)) / 2 + H- 3

    !                        U3xztL_F0_io(J, LMNH) = ( UGUxztL_F0_io(J, LMNH)  &
    !                          +3.0_WP * D1xztL_F0_io(J) * U1xztL_F0_io(J, M) * U1xztL_F0_io(J, N) * U1xztL_F0_io(J,H) &
    !                            - U1xztL_F0_io(J, M) * UGxztL_F0_io(J, LNH) &
    !                            - U1xztL_F0_io(J, N) * UGxztL_F0_io(J, LMH) &
    !                            - U1xztL_F0_io(J,H) * UGxztL_F0_io(J, LMN) &
    !                          + U1xztL_F0_io(J, M) * U1xztL_F0_io(J, N) * G1xztL_F0_io(J,H) &
    !                          + U1xztL_F0_io(J, M) * U1xztL_F0_io(J,H) * G1xztL_F0_io(J, N) &
    !                          + U1xztL_F0_io(J, N) * U1xztL_F0_io(J,H) * G1xztL_F0_io(J, M) &
    !                          + D1xztL_F0_io(J) * U1xztL_F0_io(J, M) * U2xztL_F0_io(J, LNH) &
    !                          + D1xztL_F0_io(J) * U1xztL_F0_io(J, N) * U2xztL_F0_io(J, LMH) &
    !                          + D1xztL_F0_io(J) * U1xztL_F0_io(J,H) * U2xztL_F0_io(J, LMN) ) / D1xztL_F0_io(J)
    !                    END DO
    !                END DO
    !            END DO
    !        END DO
    !        ! above treatment IS an approximat
    !   As third order introduces more dIFfeRENces than the second order, thus above approximation IS too rough to be used.

    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                DO H = 1, NDV
                    IF(N >  H) CYCLE
                    LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
                    LMN = (M * (7-M)) / 2 + N - 3
                    LMH = (M * (7-M)) / 2 + H- 3
                    LNH = (N * (7-N)) / 2 + H- 3

                    uf3_RA(J, M, N,H) = U3xztL_F0_io(J, LMNH)  &
                    - U1xztL_F0_io(J, M) * UU_RA(J, N,H) &
                    - U1xztL_F0_io(J, N) * UU_RA(J, M,H) &
                    - U1xztL_F0_io(J,H) * UU_RA(J, M, N)  &
                    + 2.0_WP * U1xztL_F0_io(J, M) * U1xztL_F0_io(J, N) * U1xztL_F0_io(J,H)
                END DO
            END DO
        END DO

        uf3_RA(J, 1, 2, 1) = uf3_RA(J, 1, 1, 2)
        uf3_RA(J, 1, 3, 1) = uf3_RA(J, 1, 1, 3)
        uf3_RA(J, 1, 3, 2) = uf3_RA(J, 1, 2, 3)

        uf3_RA(J, 2, 1, 1) = uf3_RA(J, 1, 1, 2)
        uf3_RA(J, 2, 1, 2) = uf3_RA(J, 1, 2, 2)
        uf3_RA(J, 2, 1, 3) = uf3_RA(J, 1, 2, 3)

        uf3_RA(J, 2, 2, 1) = uf3_RA(J, 1, 2, 2)
        uf3_RA(J, 2, 3, 1) = uf3_RA(J, 1, 2, 3)
        uf3_RA(J, 2, 3, 2) = uf3_RA(J, 2, 2, 3)

        uf3_RA(J, 3, 1, 1) = uf3_RA(J, 1, 1, 3)
        uf3_RA(J, 3, 1, 2) = uf3_RA(J, 1, 2, 3)
        uf3_RA(J, 3, 1, 3) = uf3_RA(J, 1, 3, 3)

        uf3_RA(J, 3, 2, 1) = uf3_RA(J, 1, 2, 3)
        uf3_RA(J, 3, 2, 2) = uf3_RA(J, 2, 2, 3)
        uf3_RA(J, 3, 2, 3) = uf3_RA(J, 2, 3, 3)

        uf3_RA(J, 3, 3, 1) = uf3_RA(J, 1, 3, 3)
        uf3_RA(J, 3, 3, 2) = uf3_RA(J, 2, 3, 3)

        !WRITE(*, '(28ES13.5)') YCC(J), uf3_RA(J, 1:3, 1:3, 1:3)  ! test

    END DO
    CALL CHKHDL(' ==>Calculated <u`_i u`_j u`_k>', MYID)

    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                DO H = 1, NDV
                    uf3d_RA(J, M, N,H) = uf3_RA(J, M, N,H) * D1xztL_F0_io(J)
                END DO
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated <\rho>*<u`_i u`_j u`_k>', MYID)


    !============= Dialation of each volumE ======================
    DO J = 1, NCL2
        dUiDXi(J) = DVDL1xztL_F0_io(J, 1, 1) + &
        DVDL1xztL_F0_io(J, 2, 2) + &
        DVDL1xztL_F0_io(J, 3, 3)
        !WRITE(*, '(2ES13.5)') YCC(J), dUiDXi(J) !test
    END DO
    CALL CHKHDL(' ==>Calculated d<u>/ DX + D<v>/ Dy+ D<w>/ Dz ', MYID)

    !============= Mean strAIn rate and vortICity tensor ==============
    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                StrAInTensor(J, M, N) = 0.5_WP * ( DVDL1xztL_F0_io(J, M, N) + DVDL1xztL_F0_io(J, N, M) )
                VortcyTensor(J, M, N) = 0.5_WP * ( DVDL1xztL_F0_io(J, M, N) - DVDL1xztL_F0_io(J, N, M) )
            END DO
        END DO
        !WRITE(*, '(7ES13.5)') YCC(J), StrAInTensor(J, 1:3, 1:3) !test
    END DO
    CALL CHKHDL(' ==>Calculated semi and asymetrIC pARts of d<u_i>/ DX_j ', MYID)

    !============== SkewnesS ========================
    DO J = 1, NCL2
        DO M = 1, NDV
            Skewness_RA(J, M) = uf3_RA(J, M, M, M) / ( DABS(uf2_RA(J, M, M))**(3.0_WP / 2.0_WP) + RealMin)
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Skewness', MYID)

    DO J = 1, NCL2
        MKE_RA(J) = 0.5_WP * ( U1xztL_F0_io(J, 1) * U1xztL_F0_io(J, 1) + &
        U1xztL_F0_io(J, 2) * U1xztL_F0_io(J, 2) + &
        U1xztL_F0_io(J, 3) * U1xztL_F0_io(J, 3) ) * D1xztL_F0_io(J)
        TKE_RA(J) = 0.5_WP * ( uf2_RA(J, 1, 1) + &
        uf2_RA(J, 2, 2) + &
        uf2_RA(J, 3, 3) )      * D1xztL_F0_io(J)
    END DO
    CALL CHKHDL(' ==>Calculated MKE_RA=<\rho>*<u_i>*<u_i>/ 2 and TKE_RA=<\rho>*<u`_I * U`_i>/ 2', MYID)


    !==== VortICitY ==================================
    ! Omega1 = \partial u_3 / \partial x_2 -  \partial u_2 / \partial x_3
    ! Omega2 = \partial u_1 / \partial x_3 -  \partial u_3 / \partial x_1
    ! Omega3 = \partial u_2 / \partial x_1 -  \partial u_1 / \partial x_2
    ! Omega_I = <Omega_i> + Omega'_i (or omega)
    ! <Omega'_i * Omega'_i> = <Omega_i * Omega_i> - <Omega_i><Omega_i>
    DO J = 1, NCL2
        ! Omega2_RA(J, 1) = (dw/ Dy * dw/ Dy) -2 * Dw/ Dy* Dv/ Dz + (dv/ Dz* Dv/ Dz)
        Omega2_RA(J, 1) = DVDL2xztL_F0_io(J, (3- 1) * NDV + 2, (3- 1) * NDV + 2) - 2.0_WP * &
        DVDL2xztL_F0_io(J, (3- 1) * NDV + 2, (2 - 1) * NDV +3) + &
        DVDL2xztL_F0_io(J, (2 - 1) * NDV +3, (2 - 1) * NDV +3)

        Omega2_RA(J, 2) = DVDL2xztL_F0_io(J, (1- 1) * NDV +3, (1- 1) * NDV +3) - 2.0_WP * &
        DVDL2xztL_F0_io(J, (1- 1) * NDV +3, (3- 1) * NDV + 1) + &
        DVDL2xztL_F0_io(J, (3- 1) * NDV + 1, (3- 1) * NDV + 1)

        Omega2_RA(J, 3) = DVDL2xztL_F0_io(J, (2 - 1) * NDV + 1, (2 - 1) * NDV + 1) - 2.0_WP * &
        DVDL2xztL_F0_io(J, (2 - 1) * NDV + 1, (1- 1) * NDV + 2) + &
        DVDL2xztL_F0_io(J, (1- 1) * NDV + 2, (1- 1) * NDV + 2)

        Omega_RA2(J, 1) = DVDL1xztL_F0_io(J, 3, 2) * DVDL1xztL_F0_io(J, 3, 2) - 2.0_WP * &
        DVDL1xztL_F0_io(J, 3, 2) * DVDL1xztL_F0_io(J, 2, 3) + &
        DVDL1xztL_F0_io(J, 2, 3) * DVDL1xztL_F0_io(J, 2, 3)

        Omega_RA2(J, 2) = DVDL1xztL_F0_io(J, 1, 3) * DVDL1xztL_F0_io(J, 1, 3) - 2.0_WP * &
        DVDL1xztL_F0_io(J, 1, 3) * DVDL1xztL_F0_io(J, 3, 1) + &
        DVDL1xztL_F0_io(J, 3, 1) * DVDL1xztL_F0_io(J, 3, 1)

        Omega_RA2(J, 3) = DVDL1xztL_F0_io(J, 2, 1) * DVDL1xztL_F0_io(J, 2, 1) - 2.0_WP * &
        DVDL1xztL_F0_io(J, 2, 1) * DVDL1xztL_F0_io(J, 1, 2) + &
        DVDL1xztL_F0_io(J, 1, 2) * DVDL1xztL_F0_io(J, 1, 2)

        omega_rms(J, 1) = DSQRT( DABS(Omega2_RA(J, 1) - Omega_RA2(J, 1)) )
        omega_rms(J, 2) = DSQRT( DABS(Omega2_RA(J, 2) - Omega_RA2(J, 2)) )
        omega_rms(J, 3) = DSQRT( DABS(Omega2_RA(J, 3) - Omega_RA2(J, 3)) )
    END DO
    CALL CHKHDL(' ==>Calculated omega_rms', MYID)

    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                AnIsotropy_RA(J, M, N) = uf2_RA(J, M, N) / (uf2_RA(J, 1, 1) + Uf2_RA(J, 2, 2) + Uf2_RA(J, 3, 3)) - &
                DBLE(Kronecker_Delta(M, N)) / 3.0_WP
            END DO
        END DO
        !            !WRITE(*, *) '  --  '
        !            !WRITE(*, *) AnIsotropy_RA(J, 1, 1:3)
        !            !WRITE(*, *) AnIsotropy_RA(J, 2, 1:3)
        !            !WRITE(*, *) AnIsotropy_RA(J, 3, 1:3)
        !            !WRITE(*, *) '  --  '
    END DO

    DO J = 1, NCL2
        ! below two methods are the same. Checked Good!
        !            Anistpinva_RA(J, 1) = 0.0_WP
        !            Anistpinva_RA(J, 2) = -( AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 1) + &
        !                                    AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 1, 2) + &
        !                                    AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 1, 3) + &
        !                                    AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 2, 1) + &
        !                                    AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 2) + &
        !                                    AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 2, 3) + &
        !                                    AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 3, 1) + &
        !                                    AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 3, 2) + &
        !                                    AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 3) ) / 2.0_WP
        !            Anistpinva_RA(J, 3) = (  AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 1) + &
        !                                    AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 1) + &
        !                                    AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 1) + &
        !                                    AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 2) + &
        !                                    AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 2) + &
        !                                    AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 2) + &
        !                                    AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 3) + &
        !                                    AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 3) + &
        !                                    AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 3) + &
        !                                    AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 1) + &
        !                                    AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 1) + &
        !                                    AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 1) + &
        !                                    AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 2) + &
        !                                    AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 2) + &
        !                                    AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 2) + &
        !                                    AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 3) + &
        !                                    AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 3) + &
        !                                    AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 3) + &
        !                                    AnIsotropy_RA(J, 1, 1) * AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 1) + &
        !                                    AnIsotropy_RA(J, 1, 2) * AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 1) + &
        !                                    AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 1) + &
        !                                    AnIsotropy_RA(J, 2, 1) * AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 2) + &
        !                                    AnIsotropy_RA(J, 2, 2) * AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 2) + &
        !                                    AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 2) + &
        !                                    AnIsotropy_RA(J, 3, 1) * AnIsotropy_RA(J, 1, 3) * AnIsotropy_RA(J, 3, 3) + &
        !                                    AnIsotropy_RA(J, 3, 2) * AnIsotropy_RA(J, 2, 3) * AnIsotropy_RA(J, 3, 3) + &
        !                                    AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 3) * AnIsotropy_RA(J, 3, 3) ) / 3.0_WP

        !            LumleyAxis_RA(J, 1) = DSQRT(-Anistpinva_RA(J, 2) / 3.0_WP)
        !            LumleyAxis_RA(J, 2) = (Anistpinva_RA(J, 3) / 2.0_WP)**(1.0_WP / 3.0_WP)
        !            !WRITE(*, *) 'invARs', J, Anistpinva_RA(J, 1:3), LumleyAxis_RA(J, 1:2)

        Matrix(1:3, 1:3) = AnIsotropy_RA(J, 1:3, 1:3)
        CALL DSYEVC3(MatriX, EIG)
        Anistpinva_RA(J, 1) = 0.0_WP
        Anistpinva_RA(J, 2) = (EIG(1) * EIG(1) + EIG(1) * EIG(2) + EIG(2) * EIG(2)) * (-1.0_WP)
        Anistpinva_RA(J, 3) = EIG(1) * EIG(2) * (EIG(1) + EIG(2)) * (-1.0_WP)
        
        !LumleyAxis_RA(J, 1) = DSQRT( DABS(-Anistpinva_RA(J, 2) / 3.0_WP))
        !LumleyAxis_RA(J, 2) = (Anistpinva_RA(J, 3) / 2.0_WP)**(1.0_WP / 3.0_WP)

        !!WRITE(*, *) 'invARs', J, Anistpinva_RA(J, 1:3), LumleyAxis_RA(J, 1:2)
        !!WRITE(*, *) 'EIGenv', J, EIG(1:3), - EIG(1) - EIG(2), EIG(1) + EIG(2) + EIG(3)

    END DO


    !====================================================================
    !======================FA based======================================
    !====================================================================
    IF(iThermoDynamics == 1) THEN
        !===============U_FA(J, M) = {u_M}==========================
        ! Eq. U_FA(J, M) = {u_M} = <\rho u_M>/ <\rho>
        !     uff_RA(J, M) = <u"_M>=<u_M>-{u_M}
        U_FA = 0.0_WP
        uff_RA = 0.0_WP
        DO J = 1, NCL2
            DO M = 1, NDV
                U_FA(J, M) = G1xztL_F0_io(J, M) / D1xztL_F0_io(J)
                uff_RA(J, M) = U1xztL_F0_io(J, M) - U_FA(J, M)
            END DO
            !WRITE(*, '(3ES13.5)') YCC(J), U_FA(J, 1:3), uff_RA(J, 1:3) !test
        END DO
        CALL CHKHDL(' ==>Calculated {U} and <u">', MYID)


        !============== DUDX_FA(CL, M, N) = d{u_m}/ DX_N ========================
        dUDX_FA = 0.0_WP
        N = 2
        DO M = 1, NDV      ! u-components

            DO J = 1, NCL2
                IF(J == 1) THEN
                    dUDX_FA(J, M, N) = ( ( YCL2ND_WFB(J + 1) * U_FA(J,  M) + &
                    YCL2ND_WFF(J + 1) * U_FA(J + 1, M) ) - 0.0_WP  ) * DYFI(J)
                ELSE IF (J == NCL2) THEN
                    dUDX_FA(J, M, N) = (  0.0_WP - &
                    ( YCL2ND_WFF(J) * U_FA(J,  M) + &
                    YCL2ND_WFB(J) * U_FA(J - 1, M) )  ) * DYFI(J)
                ELSE
                    dUDX_FA(J, M, N) = ( ( YCL2ND_WFB(J + 1) * U_FA(J,  M) + &
                    YCL2ND_WFF(J + 1) * U_FA(J + 1, M) ) - &
                    ( YCL2ND_WFF(J) * U_FA(J,  M) + &
                    YCL2ND_WFB(J) * U_FA(J - 1, M) ) ) * DYFI(J)
                END IF
            END DO
            !WRITE(*, '(4ES13.5)') YCC(J), dUDX_FA(J, 1:3, 2) !test
        END DO
        CALL CHKHDL(' ==>Calculated d{U}/ DX', MYID)

        !======UU_FA(CL, I, J) = {u_i u_j} = <\rho u_i u_j>/<\rho>  Not Equal to <u"_u u"_j>=======
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    LMN = (M * (7-M)) / 2 + N - 3
                    UU_FA(J, M, N) = UGxztL_F0_io(J, LMN) / D1xztL_F0_io(J)
                END DO
            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) THEN
                        UU_FA(J, M, N) = UU_FA(J, N, M)
                    END IF
                END DO
            END DO
            !WRITE(*, '(4ES13.5)') YCC(J), UU_FA(J, 1:3, 1:3) !test
        END DO
        CALL CHKHDL(' ==>Calculated {UU}', MYID)

        !=============={u"_i u"_j} ========================
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    uff2_FA(J, M, N) = UU_FA(J, M, N)  - U_FA(J, M) * U_FA(J, N)
                END DO
            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) THEN
                        uff2_FA(J, M, N) = uff2_FA(J, N, M)
                    END IF
                END DO
            END DO
            !WRITE(*, '(4ES13.5)') YCC(J), uff2_FA(J, 1:3, 1:3)
        END DO
        CALL CHKHDL(' ==>Calculated {u"_i u"_j }', MYID)

        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    uff2d_FA(J, M, N) = uff2_FA(J, M, N) * D1xztL_F0_io(J)
                END DO
            END DO
        END DO
        CALL CHKHDL(' ==>Calculated <\rho>*{u"_i u"_j }', MYID)

        !======{\rho u_2 u"_i u"_j}===for checking turb dIFfu, and exactly the same as the below one.
        !        DO J = 1, NCL2

        !            DO M = 1, NDV
        !                DO N = 1, NDV
        !                    IF(M >  N) CYCLE
        !                    IF(M == 1 .AND. N == 1)   LMNH =  2
        !                    IF(M == 1 .AND. N == 2)   LMNH =  4
        !                    IF(M == 1 .AND. N == 3)   LMNH =  5
        !                    IF(M == 2 .AND. N == 2)   LMNH =  7
        !                    IF(M == 2 .AND. N == 3)   LMNH =  8
        !                    IF(M == 3 .AND. N == 3)   LMNH =  9

        !                    TDIFU_FA(J, M, N) = UGUxztL_F0_io(J, LMNH) &
        !                                        - U_FA(J, M) * UU_FA(J, N, 2) * D1xztL_F0_io(J) &
        !                                        - U_FA(J, N) * UU_FA(J, M, 2) * D1xztL_F0_io(J) &
        !                                      + U_FA(J, M) * U_FA(J, N) * U_FA(J, 2) * D1xztL_F0_io(J)
        !                END DO
        !            END DO
        !            DO M = 1, NDV
        !                DO N = 1, NDV
        !                    IF(M >  N) THEN
        !                        TDIFU_FA(J, M, N) = TDIFU_FA(J, N, M)
        !                    END IF
        !                END DO
        !            END DO
        !        END DO


        !==={u"_k u"_i u"_j}==============================
        !=====u"u"u", u"u"v'', u"u"w'', u"v''v'', u"v''w'' =
        !=====u"w''w'', v''v''v'', v''v''w'', v''w''w'', w''w''w'' ===
        !{111} {112} {113}
        ![121] {122} {123}
        ![131] [132] {133}
        ![211] [212] [213]
        ![221] {222} {223}
        ![231] [232] {233}
        ![311] [312] [313]
        ![321] [322] [323]
        ![331] [332] {333}
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    DO H = 1, NDV
                        IF(N >  H) CYCLE
                        LMNH = M * (6-M) + (N * (7-N)) / 2 + H-8
                        uff3_FA(J, M, N,H) = UGUxztL_F0_io(J, LMNH) / D1xztL_F0_io(J) &
                        - U_FA(J, M) * UU_FA(J, N,H) &
                        - U_FA(J, N) * UU_FA(J, M,H) &
                        - U_FA(J,H) * UU_FA(J, M, N) &
                        + 2.0_WP * U_FA(J, M) * U_FA(J, N) * U_FA(J,H)
                        !!WRITE(*, *) M, N, H, uff3_FA(J, M, N,H),uf3_RA(J, M, N,H),uff3_FA(J, M, N,H) - Uf3_RA(J, M, N,H)
                    END DO
                END DO
            END DO

            !uff3_FA(J, 1, 1, 1)  ! #1
            !uff3_FA(J, 1, 1, 2)  ! #2
            !uff3_FA(J, 1, 1, 3)  ! #3

            uff3_FA(J, 1, 2, 1) = uff3_FA(J, 1, 1, 2) !2
            !uff3_FA(J, 1, 2, 2)  ! #4
            !uff3_FA(J, 1, 2, 3)  ! #5

            uff3_FA(J, 1, 3, 1) = uff3_FA(J, 1, 1, 3) !3
            uff3_FA(J, 1, 3, 2) = uff3_FA(J, 1, 2, 3) !5
            !uff3_FA(J, 1, 3, 3)  ! #6

            uff3_FA(J, 2, 1, 1) = uff3_FA(J, 1, 1, 2) !2
            uff3_FA(J, 2, 1, 2) = uff3_FA(J, 1, 2, 2) !4
            uff3_FA(J, 2, 1, 3) = uff3_FA(J, 1, 2, 3) !5

            uff3_FA(J, 2, 2, 1) = uff3_FA(J, 1, 2, 2) !4
            !uff3_FA(J, 2, 2, 2)  ! #7
            !uff3_FA(J, 2, 2, 3)  ! #8

            uff3_FA(J, 2, 3, 1) = uff3_FA(J, 1, 2, 3) !5
            uff3_FA(J, 2, 3, 2) = uff3_FA(J, 2, 2, 3) !8
            !uff3_FA(J, 2, 3, 3)  ! #9

            uff3_FA(J, 3, 1, 1) = uff3_FA(J, 1, 1, 3) !3
            uff3_FA(J, 3, 1, 2) = uff3_FA(J, 1, 2, 3) !5
            uff3_FA(J, 3, 1, 3) = uff3_FA(J, 1, 3, 3) !6

            uff3_FA(J, 3, 2, 1) = uff3_FA(J, 1, 2, 3) !5
            uff3_FA(J, 3, 2, 2) = uff3_FA(J, 2, 2, 3) !8
            uff3_FA(J, 3, 2, 3) = uff3_FA(J, 2, 3, 3) !9

            uff3_FA(J, 3, 3, 1) = uff3_FA(J, 1, 3, 3) !6
            uff3_FA(J, 3, 3, 2) = uff3_FA(J, 2, 3, 3) !9
            !uff3_FA(J, 3, 3, 3)  ! #10
            !WRITE(*, '(28ES13.5)') YCC(J), uff3_FA(J, 1:3, 1:3, 1:3)  ! test
        END DO

        !        !WRITE(*, *) 'Check UUU'
        !        DO J = 1, NCL2
        !!            !WRITE(*, *) TDIFU_FA(J, 1, 1), TDIFU_FA(J, 1, 2), TDIFU_FA(J, 1, 3), TDIFU_FA(J, 2, 2), TDIFU_FA(J, 2, 3), TDIFU_FA(J, 3, 3)
        !            WRITE(*, '(9ES13.5)') YCC(J), uff3_FA (J, 1, 1, 1),uf3_RA (J, 1, 1, 1), &
        !            UGUxztL_F0_io(J, 1) / D1xztL_F0_io(J), U3xztL_F0_io(J, 1), &
        !            U_FA(J, 1), U1xztL_F0_io(J, 1), &
        !            UU_FA(J, 1, 1), UU_RA(J, 1, 1)
        !        END DO
        ! Conclusion: Above IS correct, U3 indeed introduces lARge variations from UGU / D

        CALL CHKHDL(' ==>Calculated {u"_i u"_j u"_k}', MYID)

        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    DO H = 1, NDV
                        uff3d_FA(J, M, N,H) = uff3_FA(J, M, N,H) * D1xztL_F0_io(J)
                    END DO
                END DO
            END DO
        END DO
        CALL CHKHDL(' ==>Calculated <\rho> * {u"_i u"_j u"_k}', MYID)

        !============= Dialation of each volumE ======================
        DO J = 1, NCL2
            dUiDXiM(J) = (   DVDL1MxztL_F0_io(J, 1, 1) + &
            DVDL1MxztL_F0_io(J, 2, 2) + &
            DVDL1MxztL_F0_io(J, 3, 3)    ) * CVISC
        END DO
        CALL CHKHDL(' ==>Calculated <\mu dU / DX> + <\mu dV/ Dy> + <\mu dW/ Dz>', MYID)

        !============= Mean strAIn rate and vortICity tensor ==============
        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    StrAInTensorM(J, M, N) = 0.5_WP * ( DVDL1MxztL_F0_io(J, M, N) + DVDL1MxztL_F0_io(J, N, M) ) * CVISC
                    VortcyTensorM(J, M, N) = 0.5_WP * ( DVDL1MxztL_F0_io(J, M, N) - DVDL1MxztL_F0_io(J, N, M) ) * CVISC
                END DO
            END DO
        END DO
        CALL CHKHDL(' ==>Calculated \mu S and \mu\Omega', MYID)

        !============== SkewnesS ========================

        DO M = 1, NDV
            DO J = 1, NCL2
                Skewness_FA(J, M) = uff3_FA(J, M, M, M) / ( DABS(uff2_FA(J, M, M))**(3.0_WP / 2.0_WP) + RealMin)
                !WRITE(*, '(2I3.1, 3(3ES13.5, 2X))') M, J, Skewness_FA(J, M),Skewness_RA(J, M), &
                !                (Skewness_FA(J, M) -Skewness_RA(J, M)) /Skewness_FA(J, M), &
                !                 uff3_FA(J, M, M, M),uf3_RA(J, M, M, M), (uff3_FA(J, M, M, M) - Uf3_RA(J, M, M, M)) / Uff3_FA(J, M, M, M), &
                !                 uff2_FA(J, M, M),uf2_RA(J, M, M), (uff2_FA(J, M, M) - Uf2_RA(J, M, M)) / Uff2_FA(J, M, M)
            END DO
        END DO
        CALL CHKHDL(' ==>Calculated Skewness_FA', MYID)


        DO J = 1, NCL2
            MKE_FA(J) = 0.5_WP * ( U_FA(J, 1) * U_FA(J, 1) + &
            U_FA(J, 2) * U_FA(J, 2) + &
            U_FA(J, 3) * U_FA(J, 3) ) * D1xztL_F0_io(J)
            TKE_FA(J) = 0.5_WP * ( uff2_FA(J, 1, 1) + &
            uff2_FA(J, 2, 2) + &
            uff2_FA(J, 3, 3) )      * D1xztL_F0_io(J)
        END DO
        CALL CHKHDL(' ==>Calculated MKE_FA and TKE_FA', MYID)
    

        !=========== VIScous sheAR stresS =====================
        !=========<tau_mn>(<u>) AND <tau_mn>(u') ============
        ! Eq. Tau_Mean_RA(J, M, N) = viscstress_Tau_Umea + viscstress_Tau_Uper
        !     <tau_mn>(<u>) * REN = <mu> [ (\partial <u_m>) / (\partial x_n) +
        !                               (\partial <u_n>) / (\partial x_m) ] -
        !                     2 / 3<mu> [ (\partial <u_l>) / (\partial x_l) \Delta_mn]
        !     <tau_mn>(u') * REN = < mu' [ (\partial u'_m) / (\partial x_n) +
        !                               (\partial u'_n) / (\partial x_m) ] > -
        !                     2 / 3<mu' [ (\partial u'_l) / (\partial x_l) \Delta_mn]
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    !==Eq.2.20 in Huang1995====
                    Tau_Mean_RA(J, M, N) = 2.0_WP * &
                    ( StrAInTensorM(J, M, N) - dUiDXiM(J) * DBLE(Kronecker_Delta(M, N)) / 3.0_WP  )
                    Tau_meaU_RA(J, M, N) = 2.0_WP * M1xztL_F0_io(J) * CVISC * &
                    ( StrAInTensor (J, M, N) - dUiDXi(J) * DBLE(Kronecker_Delta(M, N)) / 3.0_WP  )
                END DO
            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) THEN
                        Tau_Mean_RA(J, M, N) = Tau_Mean_RA(J, N, M)
                        Tau_meaU_RA(J, M, N) = Tau_meaU_RA(J, N, M)
                    END IF
                END DO
            END DO

        END DO
    
        CALL CHKHDL(' ==>Calculated <\mu S*>, <\mu>*<S*>, <\mu` S`>', MYID)


        !======== TauU_RA(J, M, N, H) = <u_h \tau_mn>======================
        !Eq.<u_h \tau_mn> = <\mu u_h \partial{u_m}/\partial{u_n}>+
        !                  <\mu u_h \partial{u_n}/\partial{u_m}>-
        !               2 / 3<\mu u_h \partial{u_l}/\partial{u_l}>Delta_mn
        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    DO H = 1, NDV
                        TauU_RA(J, M, N, H) = ( &
                        DVDL1MUxztL_F0_io(J, M, N,H) + &
                        DVDL1MUxztL_F0_io(J, N, M,H) - &
                        2.0_WP / 3.0_WP * ( DVDL1MUxztL_F0_io(J, 1, 1,H) + &
                        DVDL1MUxztL_F0_io(J, 2, 2,H) + &
                        DVDL1MUxztL_F0_io(J, 3, 3,H) ) * DBLE(Kronecker_Delta(M, N)) ) * CVISC
                        Taufuf_RA(J, M, N, H) = TauU_RA(J, M, N, H) - &
                        U1xztL_F0_io(J,H) * Tau_Mean_RA(J, M, N)
                    END DO
                END DO
            END DO

            !===below for tesT == test OK ===
            !!WRITE(*, *) TauU_RA(J, 1, 2, 1:3) - TauU_RA(J, 2, 1, 1:3)
            !!WRITE(*, *) TauU_RA(J, 1, 3, 1:3) - TauU_RA(J, 3, 1, 1:3)
            !!WRITE(*, *) TauU_RA(J, 2, 3, 1:3) - TauU_RA(J, 3, 2, 1:3)
            !test

        END DO
        CALL CHKHDL(' ==>Calculated TauU_RA', MYID)

        !========= D<tau_mn>(<u>,u') / DY = DTaudy_RA(J, M, N) ============
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    IF(J == 1) THEN
                        IF (M == 1 .AND. N == 2) THEN
                            FF = DABS(Tauw_io(1))
                        ELSE
                            FF = 0.0_WP
                        END IF
                        dTaudy_RA(J, M, N) = &
                        ( ( YCL2ND_WFF(J + 1) * Tau_Mean_RA(J + 1, M, N) + &
                        YCL2ND_WFB(J + 1) * Tau_Mean_RA(J,  M, N) ) -  &
                        FF ) * DYFI(J)

                    ELSE IF(J == NCL2) THEN
                        IF (M == 1 .AND. N == 2) THEN
                            FF = - DABS(Tauw_io(2))
                        ELSE
                            FF = 0.0_WP
                        END IF
                        dTaudy_RA(J, M, N) = &
                        ( FF -  &
                        ( YCL2ND_WFF(J) * Tau_Mean_RA(J,  M, N) + &
                        YCL2ND_WFB(J) * Tau_Mean_RA(J - 1, M, N) ) ) * DYFI(J)

                    ELSE
                        dTaudy_RA(J, M, N) = &
                        ( ( YCL2ND_WFF(J + 1) * Tau_Mean_RA(J + 1, M, N) + &
                        YCL2ND_WFB(J + 1) * Tau_Mean_RA(J,  M, N) ) -  &
                        ( YCL2ND_WFF(J) * Tau_Mean_RA(J,  M, N) + &
                        YCL2ND_WFB(J) * Tau_Mean_RA(J - 1, M, N) ) ) * DYFI(J)
                    END IF
                END DO
            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) THEN
                        dTaudy_RA(J, M, N) = dTaudy_RA(J, N, M)
                    END IF
                END DO
            END DO

        END DO
        CALL CHKHDL(' ==>Calculated dTaudy_RA', MYID)


        !========= D<R_mn>/ DY = DTSSdy_RA(J, M, N) ============
        !========= D<\rho u"_m u"_n>/ Dy = d (<\rho> {u"_m u"_n}) / Dy
        DO J = 1, NCL2

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) CYCLE
                    IF(J == 1) THEN
                        dTSSdy_RA(J, M, N) = &
                        ( ( YCL2ND_WFF(J + 1) * Uff2d_FA(J + 1, M, N) + &
                        YCL2ND_WFB(J + 1) * Uff2d_FA(J,  M, N) ) -  &
                        0.0_WP ) * DYFI(J)

                    ELSE IF(J == NCL2) THEN
                        dTSSdy_RA(J, M, N) = &
                        ( 0.0_WP -  &
                        ( YCL2ND_WFF(J) * Uff2d_FA(J,  M, N) + &
                        YCL2ND_WFB(J) * Uff2d_FA(J - 1, M, N) ) ) * DYFI(J)

                    ELSE
                        dTSSdy_RA(J, M, N) = &
                        ( ( YCL2ND_WFF(J + 1) * Uff2d_FA(J + 1, M, N) + &
                        YCL2ND_WFB(J + 1) * Uff2d_FA(J,  M, N) ) -  &
                        ( YCL2ND_WFF(J) * Uff2d_FA(J,  M, N) + &
                        YCL2ND_WFB(J) * Uff2d_FA(J - 1, M, N) ) ) * DYFI(J)
                    END IF
                END DO
            END DO

            DO M = 1, NDV
                DO N = 1, NDV
                    IF(M >  N) THEN
                        dTSSdy_RA(J, M, N) = dTSSdy_RA(J, N, M)
                    END IF
                END DO
            END DO

        END DO
        CALL CHKHDL(' ==>Calculated dTSSdy_RA', MYID)


        !========<du_M / DX_n * \tau_hp >======================
        ! Eq. = <d(u_m) / D(x_n) * mu * d(u_h) / D(x_p) > +
        !       <d(u_m) / D(x_n) * mu * d(u_p) / D(x_h) > - 2 / 3 * Delta_hp
        !       <d(u_m) / D(x_n) * mu * d(u_l) / D(x_l) >
        ! (M - 1) * NDV + H
        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    DO H = 1, NDV
                        DO P = 1, NDV
                            TauDvDL_RA(J, M, N, H, P) = (     &
                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + N, (H- 1) * 3 + P ) + &
                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + N, (P - 1) * 3 + H ) - &
                            2.0_WP / 3.0_WP * DBLE(Kronecker_Delta(H,P)) * (&
                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + N, (1- 1) * 3 + 1 ) + &
                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + N, (2 - 1) * 3 + 2 ) + &
                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + N, (3- 1) * 3 +3 ) ) ) * CVISC

                        END DO
                    END DO
                END DO
            END DO
        END DO
        CALL CHKHDL(' ==>Calculated TauDvDL_RA', MYID)


        !        DO J = 1, NCL2
        !            DO M = 1, NDV
        !                DO N = 1, NDV
        !                    Tau_ik_Du_jDX_i_RA(J, M, N) = ( &
        !                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + 1, (N - 1) * 3 + 1 ) + &
        !                            DVDL2MxztL_F0_io(J, (M - 1) * 3 + 2, (N - 1) * 3 + 2 ) + &
        !                            DVDL2MxztL_F0_io(J, (M - 1) * 3 +3, (N - 1) * 3 +3 ) + &
        !                            DVDL2MxztL_F0_io(J, (1- 1) * 3 +M, (N - 1) * 3 + 1 ) + &
        !                            DVDL2MxztL_F0_io(J, (2 - 1) * 3 +M, (N - 1) * 3 + 2 ) + &
        !                            DVDL2MxztL_F0_io(J, (3- 1) * 3 +M, (N - 1) * 3 +3 ) - &
        !                            2.0_WP / 3.0_WP * ( &
        !                            DVDL2MxztL_F0_io(J, (1- 1) * 3 + 1, (N - 1) * 3 + 1 ) * DBLE(Kronecker_Delta(M, 1)) + &
        !                            DVDL2MxztL_F0_io(J, (1- 1) * 3 + 1, (N - 1) * 3 + 2 ) * DBLE(Kronecker_Delta(M, 2)) + &
        !                            DVDL2MxztL_F0_io(J, (1- 1) * 3 + 1, (N - 1) * 3 +3 ) * DBLE(Kronecker_Delta(M, 3)) + &
        !                            DVDL2MxztL_F0_io(J, (2 - 1) * 3 + 2, (N - 1) * 3 + 1 ) * DBLE(Kronecker_Delta(M, 1)) + &
        !                            DVDL2MxztL_F0_io(J, (2 - 1) * 3 + 2, (N - 1) * 3 + 2 ) * DBLE(Kronecker_Delta(M, 2)) + &
        !                            DVDL2MxztL_F0_io(J, (2 - 1) * 3 + 2, (N - 1) * 3 +3 ) * DBLE(Kronecker_Delta(M, 3)) + &
        !                            DVDL2MxztL_F0_io(J, (3- 1) * 3 +3, (N - 1) * 3 + 1 ) * DBLE(Kronecker_Delta(M, 1)) + &
        !                            DVDL2MxztL_F0_io(J, (3- 1) * 3 +3, (N - 1) * 3 + 2 ) * DBLE(Kronecker_Delta(M, 2)) + &
        !                            DVDL2MxztL_F0_io(J, (3- 1) * 3 +3, (N - 1) * 3 +3 ) * DBLE(Kronecker_Delta(M, 3)) ) &
        !                            ) * CVISC
        !                END DO
        !            END DO
        !        END DO

        DO J = 1, NCL2
            DO M = 1, NDV
                DO N = 1, NDV
                    AnIsotropy_FA(J, M, N) = uff2_FA(J, M, N) / ( uff2_FA(J, 1, 1) + uff2_FA(J, 2, 2) + uff2_FA(J, 3, 3) ) -&
                    DBLE(Kronecker_Delta(M, N)) / 3.0_WP
                END DO
            END DO
            !            !WRITE(*, *) '  --  '
            !            !WRITE(*, *) AnIsotropy_FA(J, 1, 1:3)
            !            !WRITE(*, *) AnIsotropy_FA(J, 2, 1:3)
            !            !WRITE(*, *) AnIsotropy_FA(J, 3, 1:3)
            !            !WRITE(*, *) AnIsotropy_FA(J, 1, 1) + AnIsotropy_FA(J, 2, 2) + AnIsotropy_FA(J, 3, 3)
            !            !WRITE(*, *) '  --  '
        END DO
        CALL CHKHDL(' ==>AnIsotropy_FA', MYID)

        DO J = 1, NCL2
            ! below two methods are the same, checked good!
            !            Anistpinva_FA(J, 1) = AnIsotropy_FA(J, 1, 1) + AnIsotropy_FA(J, 2, 2) + AnIsotropy_FA(J, 3, 3)
            !            Anistpinva_FA(J, 2) = -( AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 1) + &
            !                                    AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 1, 2) + &
            !                                    AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 1, 3) + &
            !                                    AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 2, 1) + &
            !                                    AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 2) + &
            !                                    AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 2, 3) + &
            !                                    AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 3, 1) + &
            !                                    AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 3, 2) + &
            !                                    AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 3) ) / 2.0_WP
            !            Anistpinva_FA(J, 3) = (  AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 1) + &
            !                                    AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 1) + &
            !                                    AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 1) + &
            !                                    AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 2) + &
            !                                    AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 2) + &
            !                                    AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 2) + &
            !                                    AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 3) + &
            !                                    AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 3) + &
            !                                    AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 3) + &
            !                                    AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 1) + &
            !                                    AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 1) + &
            !                                    AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 1) + &
            !                                    AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 2) + &
            !                                    AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 2) + &
            !                                    AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 2) + &
            !                                    AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 3) + &
            !                                    AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 3) + &
            !                                    AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 3) + &
            !                                    AnIsotropy_FA(J, 1, 1) * AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 1) + &
            !                                    AnIsotropy_FA(J, 1, 2) * AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 1) + &
            !                                    AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 1) + &
            !                                    AnIsotropy_FA(J, 2, 1) * AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 2) + &
            !                                    AnIsotropy_FA(J, 2, 2) * AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 2) + &
            !                                    AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 2) + &
            !                                    AnIsotropy_FA(J, 3, 1) * AnIsotropy_FA(J, 1, 3) * AnIsotropy_FA(J, 3, 3) + &
            !                                    AnIsotropy_FA(J, 3, 2) * AnIsotropy_FA(J, 2, 3) * AnIsotropy_FA(J, 3, 3) + &
            !                                    AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 3) * AnIsotropy_FA(J, 3, 3) ) / 3.0_WP

            !            LumleyAxis_FA(J, 1) = DSQRT(-Anistpinva_FA(J, 2) / 3.0_WP)
            !            LumleyAxis_FA(J, 2) = (Anistpinva_FA(J, 3) / 2.0_WP)**(1.0_WP / 3.0_WP)

            Matrix(1:3, 1:3) = AnIsotropy_FA(J, 1:3, 1:3)
            CALL DSYEVC3(MatriX, EIG)

            Anistpinva_FA(J, 1) = EIG(1) + EIG(2) + EIG(3)
            Anistpinva_FA(J, 2) = (EIG(1) * EIG(1) + EIG(1) * EIG(2) + EIG(2) * EIG(2)) * (-1.0_WP)
            Anistpinva_FA(J, 3) = EIG(1) * EIG(2) * (EIG(1) + EIG(2)) * (-1.0_WP)
            LumleyAxis_FA(J, 1) = DSQRT(DABS(-Anistpinva_FA(J, 2) / 3.0_WP))            ! \eta
            LumleyAxis_FA(J, 2) = SIGN( DABS(Anistpinva_FA(J, 3) / 2.0_WP)**(1.0_WP / 3.0_WP),  Anistpinva_FA(J, 3) / 2.0_WP) ! \xi

            !            !WRITE(*, *) 'invARs', J, Anistpinva_FA(J, 1:3), LumleyAxis_FA(J, 1:2)
            !            !WRITE(*, *) 'EIGenv', J, EIG(1:3), - EIG(1) - EIG(2), EIG(1) + EIG(2) + EIG(3)

        END DO
        CALL CHKHDL(' ==>LumleyAxis_FA', MYID)
    

        ! calculate DrivenForce in Streamwise direction for periodic directions.
        DO J = 1, NCL2

            IF(J == 1) THEN
                FF = ( ( YCL2ND_WFF(J + 1) * (Tau_Mean_RA(J + 1, 1, 2) - Uff2d_FA(J + 1, 1, 2)) + &
                YCL2ND_WFB(J + 1) * (Tau_Mean_RA(J,  1, 2) - Uff2d_FA(J,  1, 2)) ) -  &
                DABS(Tauw_io(1)) ) * DYFI(J)


            ELSE IF(J == NCL2) THEN
                FF = ( - DABS(Tauw_io(2)) -  &
                ( YCL2ND_WFF(J) * (Tau_Mean_RA(J,  1, 2) - Uff2d_FA(J,  1, 2)) + &
                YCL2ND_WFB(J) * (Tau_Mean_RA(J - 1, 1, 2) - Uff2d_FA(J - 1, 1, 2)) ) ) * DYFI(J)


            ELSE
                FF = ( ( YCL2ND_WFF(J + 1) * (Tau_Mean_RA(J + 1, 1, 2) - Uff2d_FA(J + 1, 1, 2)) + &
                YCL2ND_WFB(J + 1) * (Tau_Mean_RA(J,  1, 2) - Uff2d_FA(J,  1, 2)) ) -  &
                ( YCL2ND_WFF(J) * (Tau_Mean_RA(J,  1, 2) - Uff2d_FA(J,  1, 2)) + &
                YCL2ND_WFB(J) * (Tau_Mean_RA(J - 1, 1, 2) - Uff2d_FA(J - 1, 1, 2)) ) ) * DYFI(J)

            END IF

            DrivenForce(J) = (FF+F_A* IBuoF(1) * D1xztL_F0_io(J)) * (-1.0_WP)
            !!WRITE(*, *) YCC(J), DrivenForce(J), FF, F_A* IBuoF(1) * D1xztL_F0_io(J)
        END DO

        DrivenFCTT1  = 0.0_WP
        DrivenFCTT2  = 0.0_WP
        DrivenFCTTU1 = 0.0_WP
        DrivenFCTTU2 = 0.0_WP
        DO J = 1, NCL2
            DrivenFCTT1 = DrivenFCTT1 + DrivenForce(J) / DYFI(J)
            DrivenFCTTU1 = DrivenFCTTU1 + DrivenForce(J) * U_FA(J, 1) / DYFI(J)

            DrivenFCTT2 = DrivenFCTT2 +FUxztL_F0_io(J, 4) / DYFI(J)
            DrivenFCTTU2 = DrivenFCTTU2 +FUxztL_F0_io(J, 4) * U_FA(J, 1) / DYFI(J)
        END DO

        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        OPEN(TECFLG, FILE = TRIM(FilePath4) // 'Result.IO.Table.WallandBulk.' // TRIM(PNTIM) // '.plt', POSITION = 'APPEND')
        WRITE(TECFLG, '(A)') '==================================== '
        WRITE(TECFLG, '(A, 2ES20.7)') 'Driven Force (Constant) & its MKE production = ', DrivenFCTT2, DrivenFCTTU2
        WRITE(TECFLG, '(A, 2ES20.7)') 'Driven Force (inversed) & its MKE production = ', DrivenFCTT1, DrivenFCTTU1
        CLOSE(TECFLG)

        CALL CHKHDL(' ==>DrivenForce', MYID)


        !=============================================================
        ddenintg = 0.0_WP
        deNMintg = 0.0_WP
        bdfciNTG = 0.0_WP
        DensIntg = 0.0_WP
        DO J = 1, NCL2
            ! second order intgeral
            !            IF(J == 1) THEN
            !                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
            !                                  Dwal(1) )
            !            ELSE IF(J == NCL2) THEN
            !                DENtemp = 0.5_WP * ( Dwal(2) +             &
            !                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
            !            ELSE
            !                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
            !                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
            !            END IF
            !first order integral
            DENtemp = D1xztL_F0_io(J)

            deNMintg = deNMintg+ DENtemP / DYFI(J)
            DensIntg(J) = deNMintg !; !WRITE(*, *) J, DENtemp, DensIntg(J)

            ddenintg = ddenintg + (DENtemP - DenAvew) / DYFI(J)
            bdfcintg(J) = F_A* Ddenintg
        END DO
        CALL CHKHDL(' ==>Calculated bodyforce distribution', MYID)
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE PP_Budg_INIT
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    !================== RuV =============================
    Budg_prodc_stres_DuiuJ = 0.0_WP
    Budg_viscs_dissp_DuiuJ = 0.0_WP
    Budg_pduDX_stran_DuiuJ = 0.0_WP
    Budg_Turbu_diffu_DuiuJ = 0.0_WP
    Budg_DpuDX_diffu_DuiuJ = 0.0_WP
    Budg_viscs_diffu_DuiuJ = 0.0_WP

    Budg_press_accl1_DuiuJ = 0.0_WP
    Budg_viscs_accl1_DuiuJ = 0.0_WP
    Budg_prodc_Dvfc1_DuiuJ = 0.0_WP
    Budg_Balance1_Duiuj  = 0.0_WP

    Budg_prodc_gvfc2_DuiuJ = 0.0_WP
    Budg_prodc_Dvfc2_DuiuJ = 0.0_WP
    Budg_turss_accl2_DuiuJ = 0.0_WP
    Budg_Balance2_Duiuj  = 0.0_WP

    Budg_pressure3_DuiuJ = 0.0_WP
    Budg_vistress3_DuiuJ = 0.0_WP
    Budg_Balance3_Duiuj  = 0.0_WP

    !================== RuV == Sum along Y =======================
    Budg_prodc_stres_Duiuj_ysum = 0.0_WP
    Budg_viscs_dissp_Duiuj_ysum = 0.0_WP
    Budg_pduDX_stran_Duiuj_ysum = 0.0_WP
    Budg_Turbu_diffu_Duiuj_ysum = 0.0_WP
    Budg_DpuDX_diffu_Duiuj_ysum = 0.0_WP
    Budg_viscs_diffu_Duiuj_ysum = 0.0_WP

    Budg_press_accl1_Duiuj_ysum = 0.0_WP
    Budg_viscs_accl1_Duiuj_ysum = 0.0_WP
    Budg_prodc_Dvfc1_Duiuj_ysum = 0.0_WP
    Budg_Balance1_Duiuj_ysum = 0.0_WP

    Budg_prodc_gvfc2_Duiuj_ysum = 0.0_WP
    Budg_prodc_Dvfc2_Duiuj_ysum = 0.0_WP
    Budg_turss_accl2_Duiuj_ysum = 0.0_WP
    Budg_Balance2_Duiuj_ysum  = 0.0_WP

    Budg_pressure3_Duiuj_ysum = 0.0_WP
    Budg_vistress3_Duiuj_ysum = 0.0_WP
    Budg_Balance3_Duiuj_ysum  = 0.0_WP

    !================== TKE =============================
    Budg_prodc_stres_TKE = 0.0_WP
    Budg_viscs_dissp_TKE = 0.0_WP
    Budg_pduDX_stran_TKE = 0.0_WP
    Budg_Turbu_diffu_TKE = 0.0_WP
    Budg_DpuDX_diffu_TKE = 0.0_WP
    Budg_viscs_diffu_TKE = 0.0_WP

    Budg_press_accl1_TKE = 0.0_WP
    Budg_viscs_accl1_TKE = 0.0_WP
    Budg_prodc_Dvfc1_TKE = 0.0_WP
    Budg_Balance1_TKE  = 0.0_WP

    Budg_prodc_gvfc2_TKE = 0.0_WP
    Budg_prodc_Dvfc2_TKE = 0.0_WP
    Budg_turss_accl2_TKE = 0.0_WP
    Budg_Balance2_TKE  = 0.0_WP

    Budg_pressure3_TKE = 0.0_WP
    Budg_vistress3_TKE = 0.0_WP
    Budg_Balance3_TKE  = 0.0_WP

    !================== TKE sum along Y =============================
    Budg_prodc_stres_TKE_ysum = 0.0_WP
    Budg_viscs_dissp_TKE_ysum = 0.0_WP
    Budg_pduDX_stran_TKE_ysum = 0.0_WP
    Budg_Turbu_diffu_TKE_ysum = 0.0_WP
    Budg_DpuDX_diffu_TKE_ysum = 0.0_WP
    Budg_viscs_diffu_TKE_ysum = 0.0_WP

    Budg_press_accl1_TKE_ysum = 0.0_WP
    Budg_viscs_accl1_TKE_ysum = 0.0_WP
    Budg_prodc_Dvfc1_TKE_ysum = 0.0_WP
    Budg_Balance1_TKE_ysum  = 0.0_WP

    Budg_prodc_gvfc2_TKE_ysum = 0.0_WP
    Budg_prodc_Dvfc2_TKE_ysum = 0.0_WP
    Budg_turss_accl2_TKE_ysum = 0.0_WP
    Budg_Balance2_TKE_ysum  = 0.0_WP

    Budg_pressure3_TKE_ysum = 0.0_WP
    Budg_vistress3_TKE_ysum = 0.0_WP
    Budg_Balance3_TKE_ysum  = 0.0_WP

    !================== MKE =============================
    Budg_prodc_stres_MKE = 0.0_WP
    Budg_viscs_dissp_MKE = 0.0_WP
    Budg_pduDX_stran_MKE = 0.0_WP
    Budg_Turbu_diffu_MKE = 0.0_WP
    Budg_DpuDX_diffu_MKE = 0.0_WP
    Budg_viscs_diffu_MKE = 0.0_WP

    Budg_press_accl1_MKE = 0.0_WP
    Budg_viscs_accl1_MKE = 0.0_WP
    Budg_prodc_Dvfc1_MKE = 0.0_WP
    Budg_Balance1_MKE  = 0.0_WP

    Budg_prodc_gvfc1_MKE = 0.0_WP
    Budg_prodc_Dvfc2_MKE = 0.0_WP
    Budg_turss_accl2_MKE = 0.0_WP
    Budg_Balance2_MKE  = 0.0_WP

    Budg_pressure3_MKE = 0.0_WP
    Budg_vistress3_MKE = 0.0_WP
    Budg_Balance3_MKE  = 0.0_WP

    !================== MKE sum along Y =============================
    Budg_prodc_stres_MKE_ysum = 0.0_WP
    Budg_viscs_dissp_MKE_ysum = 0.0_WP
    Budg_pduDX_stran_MKE_ysum = 0.0_WP
    Budg_Turbu_diffu_MKE_ysum = 0.0_WP
    Budg_DpuDX_diffu_MKE_ysum = 0.0_WP
    Budg_viscs_diffu_MKE_ysum = 0.0_WP

    Budg_press_accl1_MKE_ysum = 0.0_WP
    Budg_viscs_accl1_MKE_ysum = 0.0_WP
    Budg_prodc_Dvfc1_MKE_ysum = 0.0_WP
    Budg_prodc_gvfc1_MKE_ysum = 0.0_WP
    Budg_Balance1_MKE_ysum = 0.0_WP

    Budg_prodc_gvfc2_MKE_ysum = 0.0_WP
    Budg_prodc_Dvfc2_MKE_ysum = 0.0_WP
    Budg_turss_accl2_MKE_ysum = 0.0_WP
    Budg_Balance2_MKE_ysum = 0.0_WP

    Budg_pressure3_MKE_ysum = 0.0_WP
    Budg_vistress3_MKE_ysum = 0.0_WP
    Budg_Balance3_MKE_ysum = 0.0_WP

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!###########################  FA   ##############################################################
SUBROUTINE PP_FLOW_FA_RSTE_Budg_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, L, k
    INTEGER(4) :: TECflg = 200
    REAL(WP) :: FF, TT
    REAL(WP) :: Budg_prod_Pij1, Budg_prod_Pij2, coe
    REAL(WP) :: vISc_DISsipation_Duiduj1, vISc_DISsipation_Duiduj2
    REAL(WP) :: sum_p_related_mke
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists

    IF(iThermoDynamics /= 1) RETURN
    CALL CHKHDL('=====Calculating FA BudegtS ===== ', MYID)
    CALL PP_Budg_INIT

    !===FA====PRODUCTION TERMS ==P ==================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !       Eq: Budg_prodc_stres_Duiuj(J, L) = p *_{ij}=P_{ij}-2 / 3P\Delta_{ij}
    !       Eq: P_{ij} = -<\rho u"_i u"_k> (\partial {u_j} / \partial x_k) +
    !                    -<\rho u"_j u"_k> (\partial {u_i} / \partial x_k)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_prodc_stres_MKE(J) = uff2d_FA(J, 1, 2) * DUDX_FA(J, 1, 2) + &
        uff2d_FA(J, 2, 2) * DUDX_FA(J, 2, 2) + &
        uff2d_FA(J, 3, 2) * DUDX_FA(J, 3, 2)
        Budg_prodc_stres_TKE(J) = -1.0_WP *Budg_prodc_stres_MKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                Budg_prod_Pij1 = uff2d_FA(J, M, 2) * DUDX_FA(J, N, 2)
                Budg_prod_Pij2 = uff2d_FA(J, N, 2) * DUDX_FA(J, M, 2)
                Budg_prodc_stres_Duiuj(J, L) = -1.0_WP * (Budg_prod_Pij1 + Budg_prod_Pij2)

            END DO
        END DO

    END DO
    CALL CHKHDL(' ==>Calculated Budg_prodc_stres', MYID)

    !==FA===== VIScous ENERGY dISsipation term ===========================================================
    !       TauDvDL_RA(J, M, N,H,P) = <\partial(u_m) /\partial x_n \tau_hp>
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_viscs_dissp_TKE(J) = &
        TauDvDL_RA(J, 1, 1, 1, 1) - Tau_Mean_RA(J, 1, 1) * DVDL1xztL_F0_io(J, 1, 1) + &
        TauDvDL_RA(J, 2, 1, 1, 2) - Tau_Mean_RA(J, 1, 2) * DVDL1xztL_F0_io(J, 2, 1) + &
        TauDvDL_RA(J, 3, 1, 1, 3) - Tau_Mean_RA(J, 1, 3) * DVDL1xztL_F0_io(J, 3, 1) + &
        TauDvDL_RA(J, 1, 2, 2, 1) - Tau_Mean_RA(J, 2, 1) * DVDL1xztL_F0_io(J, 1, 2) + &
        TauDvDL_RA(J, 2, 2, 2, 2) - Tau_Mean_RA(J, 2, 2) * DVDL1xztL_F0_io(J, 2, 2) + &
        TauDvDL_RA(J, 3, 2, 2, 3) - Tau_Mean_RA(J, 2, 3) * DVDL1xztL_F0_io(J, 3, 2) + &
        TauDvDL_RA(J, 1, 3, 3, 1) - Tau_Mean_RA(J, 3, 1) * DVDL1xztL_F0_io(J, 1, 3) + &
        TauDvDL_RA(J, 2, 3, 3, 2) - Tau_Mean_RA(J, 3, 2) * DVDL1xztL_F0_io(J, 2, 3) + &
        TauDvDL_RA(J, 3, 3, 3, 3) - Tau_Mean_RA(J, 3, 3) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_viscs_dissp_TKE(J) = -1.0_WP *Budg_viscs_dissp_TKE(J)

        Budg_viscs_dissp_MKE(J) = &
        Tau_Mean_RA(J, 1, 1) * DVDL1xztL_F0_io(J, 1, 1) + &
        Tau_Mean_RA(J, 2, 1) * DVDL1xztL_F0_io(J, 2, 1) + &
        Tau_Mean_RA(J, 3, 1) * DVDL1xztL_F0_io(J, 3, 1) + &
        Tau_Mean_RA(J, 1, 2) * DVDL1xztL_F0_io(J, 1, 2) + &
        Tau_Mean_RA(J, 2, 2) * DVDL1xztL_F0_io(J, 2, 2) + &
        Tau_Mean_RA(J, 3, 2) * DVDL1xztL_F0_io(J, 3, 2) + &
        Tau_Mean_RA(J, 1, 3) * DVDL1xztL_F0_io(J, 1, 3) + &
        Tau_Mean_RA(J, 2, 3) * DVDL1xztL_F0_io(J, 2, 3) + &
        Tau_Mean_RA(J, 3, 3) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_viscs_dissp_MKE(J) = -1.0_WP *Budg_viscs_dissp_MKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3  !Tau_ik_Du_jDX_i_RA(J, M, N) + &Tau_ik_Du_jDX_i_RA(J, N, M) +  &
                vISc_DISsipation_Duiduj1 = &
                TauDvDL_RA(J, N, 1, M, 1) + &
                TauDvDL_RA(J, N, 2, M, 2) + &
                TauDvDL_RA(J, N, 3, M, 3) - &
                Tau_Mean_RA(J, M, 1) * DVDL1xztL_F0_io(J, N, 1) - &
                Tau_Mean_RA(J, M, 2) * DVDL1xztL_F0_io(J, N, 2) - &
                Tau_Mean_RA(J, M, 3) * DVDL1xztL_F0_io(J, N, 3)
                vISc_DISsipation_Duiduj2 = &
                TauDvDL_RA(J, M, 1, N, 1) + &
                TauDvDL_RA(J, M, 2, N, 2) + &
                TauDvDL_RA(J, M, 3, N, 3) - &
                Tau_Mean_RA(J, N, 1) * DVDL1xztL_F0_io(J, M, 1) - &
                Tau_Mean_RA(J, N, 2) * DVDL1xztL_F0_io(J, M, 2) - &
                Tau_Mean_RA(J, N, 3) * DVDL1xztL_F0_io(J, M, 3)
                Budg_viscs_dissp_Duiuj(J, L) = -1.0_WP * (vISc_DISsipation_Duiduj1 + VISc_DISsipation_Duiduj2)
                ! DO not need times DENSITY, as it goes into the equation. See Huang1995.
                !!WRITE(*, *) L, vISc_DISsipation_Duiduj1, vISc_DISsipation_Duiduj2
            END DO
        END DO

    END DO
    CALL CHKHDL(' ==>Calculated Budg_viscs_dissp_Duiuj', MYID)

    !==FA===== Velocity- PRessure gradient pressure strAIn terM ======= (StrAIn ) =======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !       Eq = <p' \partial { u"_i}/\partial{x_j} > +
    !            <p' \partial { u"_j}/\partial{x_i} >
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_pduDX_stran_TKE(J) =  DVDLPxztL_F0_io(J, 1, 1) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 1, 1) + &
        DVDLPxztL_F0_io(J, 2, 2) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 2, 2) + &
        DVDLPxztL_F0_io(J, 3, 3) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_pduDX_stran_MKE(J) =  U1xztL_F0_io(J, 4) * dUiDXi(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_pduDX_stran_Duiuj(J, L) = DVDLPxztL_F0_io(J, M, N) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, M, N) + &
                DVDLPxztL_F0_io(J, N, M) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, N, M)
            END DO
        END DO
        !!WRITE(*, *) J, YCC(J), DVDLPxztL_F0_io(J, 1, 1),U1xztL_F0_io(J, 4), DVDL1xztL_F0_io(J, 1, 1), Budg_pduDX_stran_Duiuj(J, 1) !test
    END DO
    CALL CHKHDL(' ==>Calculated Budg_pduDX_stran_Duiuj', MYID)


    !==FA===== TURBUELCEN DIFFUSION TERMS = ( Turb. DIFfusion) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !       Eq.  Budg_TdIFf_Duiuj(J, L) = ( \partial <\rho u"_i u"_j u"_k > ) / (\partial x_k )
    DO J = 1, NCL2
        uffMKEffd_FA(J) = U_FA(J, 1) * Uff2d_FA(J, 1, 2) + &
        U_FA(J, 2) * Uff2d_FA(J, 2, 2) + &
        U_FA(J, 3) * Uff2d_FA(J, 3, 2)

        uffTKEffd_FA(J) = 0.5_WP * ( uff3d_FA(J, 1, 1, 2) + uff3d_FA(J, 2, 2, 2) + uff3d_FA(J, 3, 3, 2) )
    END DO

    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_Turbu_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UffTKEffd_FA(J + 1) + &
            YCL2ND_WFB(J + 1) * UffTKEffd_FA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UffMKEffd_FA(J + 1) + &
            YCL2ND_WFB(J + 1) * UffMKEffd_FA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
        ELSE IF (J == NCL2) THEN
            Budg_Turbu_diffu_TKE(J) = &
            ( 0.0_WP - &
            ( YCL2ND_WFF(J) * UffTKEffd_FA(J  ) + &
            YCL2ND_WFB(J) * UffTKEffd_FA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( 0.0_WP - &
            ( YCL2ND_WFF(J) * UffMKEffd_FA(J  ) + &
            YCL2ND_WFB(J) * UffMKEffd_FA(J - 1) ) ) * DYFI(J) * (-1.0_WP)

        ELSE
            Budg_Turbu_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UffTKEffd_FA(J + 1) + &
            YCL2ND_WFB(J + 1) * UffTKEffd_FA(J  ) ) -         &
            ( YCL2ND_WFF(J) * UffTKEffd_FA(J  ) + &
            YCL2ND_WFB(J) * UffTKEffd_FA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UffMKEffd_FA(J + 1) + &
            YCL2ND_WFB(J + 1) * UffMKEffd_FA(J  ) ) -         &
            ( YCL2ND_WFF(J) * UffMKEffd_FA(J  ) + &
            YCL2ND_WFB(J) * UffMKEffd_FA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
        END IF

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                IF(J == 1) THEN
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * Uff3d_FA(J + 1, M, N, 2) + &
                    YCL2ND_WFB(J + 1) * Uff3d_FA(J,  M, N, 2) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
                ELSE IF (J == NCL2) THEN
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( 0.0_WP - &
                    ( YCL2ND_WFF(J) * Uff3d_FA(J,  M, N, 2) + &
                    YCL2ND_WFB(J) * Uff3d_FA(J - 1, M, N, 2) ) ) * DYFI(J) * (-1.0_WP)

                ELSE
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * Uff3d_FA(J + 1, M, N, 2) + &
                    YCL2ND_WFB(J + 1) * Uff3d_FA(J,  M, N, 2) ) -         &
                    ( YCL2ND_WFF(J) * Uff3d_FA(J,  M, N, 2) + &
                    YCL2ND_WFB(J) * Uff3d_FA(J - 1, M, N, 2) ) ) * DYFI(J) * (-1.0_WP)
                END IF


                !                    !======= TesT ==checked the same as abovE ==================
                !                    FF =Budg_Turbu_diffu_Duiuj(J, L)
                !                    IF(J == 1) THEN
                !                        Budg_Turbu_diffu_Duiuj(J, L) = &
                !                            ( ( YCL2ND_WFF(J + 1) * TDIFU_FA(J + 1, M, N) + &
                !                                YCL2ND_WFB(J + 1) * TDIFU_FA(J,  M, N) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
                !                    ELSE IF (J == NCL2) THEN
                !                        Budg_Turbu_diffu_Duiuj(J, L) = &
                !                            ( 0.0_WP - &
                !                              ( YCL2ND_WFF(J) * TDIFU_FA(J,  M, N) + &
                !                                YCL2ND_WFB(J) * TDIFU_FA(J - 1, M, N) ) ) * DYFI(J) * (-1.0_WP)

                !                    ELSE
                !                        Budg_Turbu_diffu_Duiuj(J, L) = &
                !                            ( ( YCL2ND_WFF(J + 1) * TDIFU_FA(J + 1, M, N) + &
                !                                YCL2ND_WFB(J + 1) * TDIFU_FA(J,  M, N) ) -         &
                !                              ( YCL2ND_WFF(J) * TDIFU_FA(J,  M, N) + &
                !                                YCL2ND_WFB(J) * TDIFU_FA(J - 1, M, N) ) ) * DYFI(J) * (-1.0_WP)
                !                    END IF
                !                    !WRITE(*, *) M, N, J, Budg_Turbu_diffu_Duiuj(J, L) -FF
                !                    !======== Test

            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_Turbu_diffu_Duiuj', MYID)


    !==FA===== Velocity- PRessure gradient dIFfusion terM ======= (pressure dIFfusion) ==================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !       Eq. = - \partial (<p' u"_j>) /\partial (x_i) - \partial (<p' u"_i>) /\partial (x_j)
    !         = - \partial (<p' u'_j>) /\partial (x_i) - \partial (<p' u'_i>) /\partial (x_j)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_DpuDX_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * Ufpf_RA(J + 1, 2) + &
            YCL2ND_WFB(J + 1) * Ufpf_RA(J,  2) ) &
            -  0.0_WP ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * (U1xztL_F0_io(J + 1, 2) * U1xztL_F0_io(J + 1, 4)) + &
            YCL2ND_WFB(J + 1) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) ) &
            -  0.0_WP ) * DYFI(J) * (-1.0_WP)

        ELSE IF(J == NCL2) THEN
            Budg_DpuDX_diffu_TKE(J) = &
            ( 0.0_WP -  &
            ( YCL2ND_WFF(J) * Ufpf_RA(J,  2) + &
            YCL2ND_WFB(J) * Ufpf_RA(J - 1, 2) ) ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( 0.0_WP -  &
            ( YCL2ND_WFF(J) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) + &
            YCL2ND_WFB(J) * (U1xztL_F0_io(J - 1, 2) * U1xztL_F0_io(J - 1, 4)) )) * DYFI(J) * (-1.0_WP)
        ELSE

            Budg_DpuDX_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * Ufpf_RA(J + 1, 2) + &
            YCL2ND_WFB(J + 1) * Ufpf_RA(J,  2) ) -  &
            ( YCL2ND_WFF(J) * Ufpf_RA(J,  2) + &
            YCL2ND_WFB(J) * Ufpf_RA(J - 1, 2) ) ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * (U1xztL_F0_io(J + 1, 2) * U1xztL_F0_io(J + 1, 4)) + &
            YCL2ND_WFB(J + 1) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) ) -  &
            ( YCL2ND_WFF(J) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) + &
            YCL2ND_WFB(J) * (U1xztL_F0_io(J - 1, 2) * U1xztL_F0_io(J - 1, 4)) ) ) * DYFI(J) * (-1.0_WP)
        END IF


        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                IF(J == 1) THEN
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * ( ufpf_RA(J + 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J + 1, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J + 1) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) ) - &
                    0.0_WP ) * DYFI(J)
                ELSE IF(J == NCL2) THEN
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    (  0.0_WP - &
                    ( YCL2ND_WFF(J) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J) * ( ufpf_RA(J - 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J - 1, N) * DBLE(Kronecker_Delta(M, 2))   ) ) ) * DYFI(J)
                ELSE
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * ( ufpf_RA(J + 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J + 1, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J + 1) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) ) - &
                    ( YCL2ND_WFF(J) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J) * ( ufpf_RA(J - 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J - 1, N) * DBLE(Kronecker_Delta(M, 2))   ) ) ) * DYFI(J)
                END IF


                Budg_DpuDX_diffu_Duiuj(J, L) = Budg_DpuDX_diffu_Duiuj(J, L) * (-1.0_WP)
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_DpuDX_diffu_Duiuj', MYID)


    !===FA==== VIScous dIFfusion term ===========================================================
    !       VIScous stress term IS based on RA decomPOSITION like Huang, not FA.
    !       Eq. = \partial <u"_j tau'_ki> / pARtial (x_k) +
    !             \partial <u"_i tau'_kj> / pARtial (x_k)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_viscs_diffu_TKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J + 1) * ( Taufuf_RA(J + 1, 1, 2, 1) + &
            Taufuf_RA(J + 1, 2, 2, 2) + &
            Taufuf_RA(J + 1, 3, 2, 3) ) ) -&
            0.0_WP &
            ) * DYFI(J)
            Budg_viscs_diffu_MKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J + 1) * ( U1xztL_F0_io(J + 1, 1) * Tau_Mean_RA(J + 1, 1, 2) + &
            U1xztL_F0_io(J + 1, 2) * Tau_Mean_RA(J + 1, 2, 2) + &
            U1xztL_F0_io(J + 1, 3) * Tau_Mean_RA(J + 1, 3, 2) ) ) -&
            0.0_WP &
            ) * DYFI(J)
        ELSE IF (J == NCL2) THEN
            Budg_viscs_diffu_TKE(J) = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J) * ( Taufuf_RA(J - 1, 1, 2, 1) + &
            Taufuf_RA(J - 1, 2, 2, 2) + &
            Taufuf_RA(J - 1, 3, 2, 3) ) ) &
            ) * DYFI(J)

            Budg_viscs_diffu_MKE(J) = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J, 1,  2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J, 2,  2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J, 3,  2) ) &
            + YCL2ND_WFB(J) * ( U1xztL_F0_io(J - 1, 1) * Tau_Mean_RA(J - 1, 1, 2) + &
            U1xztL_F0_io(J - 1, 2) * Tau_Mean_RA(J - 1, 2, 2) + &
            U1xztL_F0_io(J - 1, 3) * Tau_Mean_RA(J - 1, 3, 2) ) ) &
            ) * DYFI(J)
        ELSE
            Budg_viscs_diffu_TKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J + 1) * ( Taufuf_RA(J + 1, 1, 2, 1) + &
            Taufuf_RA(J + 1, 2, 2, 2) + &
            Taufuf_RA(J + 1, 3, 2, 3) ) ) -&
            ( YCL2ND_WFF(J) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J) * ( Taufuf_RA(J - 1, 1, 2, 1) + &
            Taufuf_RA(J - 1, 2, 2, 2) + &
            Taufuf_RA(J - 1, 3, 2, 3) ) ) &
            ) * DYFI(J)

            Budg_viscs_diffu_MKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J + 1) * ( U1xztL_F0_io(J + 1, 1) * Tau_Mean_RA(J + 1, 1, 2) + &
            U1xztL_F0_io(J + 1, 2) * Tau_Mean_RA(J + 1, 2, 2) + &
            U1xztL_F0_io(J + 1, 3) * Tau_Mean_RA(J + 1, 3, 2) ) ) -&
            ( YCL2ND_WFF(J) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J) * ( U1xztL_F0_io(J - 1, 1) * Tau_Mean_RA(J - 1, 1, 2) + &
            U1xztL_F0_io(J - 1, 2) * Tau_Mean_RA(J - 1, 2, 2) + &
            U1xztL_F0_io(J - 1, 3) * Tau_Mean_RA(J - 1, 3, 2) ) ) &
            ) * DYFI(J)
        END IF

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                IF(J == 1) THEN
                    Budg_viscs_diffu_Duiuj(J, L) = ( &
                    ( YCL2ND_WFF(J + 1) * Taufuf_RA(J,  2, M, N)   &
                    + YCL2ND_WFB(J + 1) * Taufuf_RA(J + 1, 2, M, N) ) -&
                    0.0_WP &
                    ) * DYFI(J) + (&
                    ( YCL2ND_WFF(J + 1) * Taufuf_RA(J,  2, N, M)   &
                    + YCL2ND_WFB(J + 1) * Taufuf_RA(J + 1, 2, N, M) ) -&
                    0.0_WP &
                    ) * DYFI(J)
                ELSE IF (J == NCL2) THEN
                    Budg_viscs_diffu_Duiuj(J, L) = ( &
                    0.0_WP -&
                    ( YCL2ND_WFF(J) * Taufuf_RA(J,  2, M, N)   &
                    + YCL2ND_WFB(J) * Taufuf_RA(J - 1, 2, M, N) ) &
                    ) * DYFI(J) + (&
                    0.0_WP -&
                    ( YCL2ND_WFF(J) * Taufuf_RA(J,  2, N, M)   &
                    + YCL2ND_WFB(J) * Taufuf_RA(J - 1, 2, N, M) ) &
                    ) * DYFI(J)
                ELSE
                    Budg_viscs_diffu_Duiuj(J, L) = ( &
                    ( YCL2ND_WFF(J + 1) * Taufuf_RA(J,  2, M, N)   &
                    + YCL2ND_WFB(J + 1) * Taufuf_RA(J + 1, 2, M, N) ) -&
                    ( YCL2ND_WFF(J) * Taufuf_RA(J,  2, M, N)   &
                    + YCL2ND_WFB(J) * Taufuf_RA(J - 1, 2, M, N) ) &
                    ) * DYFI(J) + (&
                    ( YCL2ND_WFF(J + 1) * Taufuf_RA(J,  2, N, M)   &
                    + YCL2ND_WFB(J + 1) * Taufuf_RA(J + 1, 2, N, M) ) -&
                    ( YCL2ND_WFF(J) * Taufuf_RA(J,  2, N, M)   &
                    + YCL2ND_WFB(J) * Taufuf_RA(J - 1, 2, N, M) ) &
                    ) * DYFI(J)
                END IF
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_viscs_diffu_Duiuj', MYID)
    IF(iThermoDynamics == 1) THEN
    !====FA==========below 1 first decomPOSITION method ===========================================
    !==FA=====pressure acceleration term, including alL ======== (FA) ===========>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_press_accl1_MKE(J) = uff_RA(J, 1) * DPDX_RA(J, 1) + uff_RA(J, 2) * DPDX_RA(J, 2) + uff_RA(J, 3) * DPDX_RA(J, 3)
        Budg_press_accl1_TKE(J) = -1.0_WP *Budg_press_accl1_MKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_press_accl1_Duiuj(J, L) = - DPDX_RA(J, M) * Uff_RA(J, N) &
                - DPDX_RA(J, N) * Uff_RA(J, M)
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_press_accl1_Duiuj', MYID)

    !==FA===== VIScous acceleration, including alL ======== (based on tau_RA) ======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !       Eq. <u"_j> (\partial <tau_ki>) / (\partial  x_k) +
    !           <u"_i> (\partial <tau_kj>) / (\partial  x_k)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_viscs_accl1_TKE(J) = dTaudy_RA(J, 2, 1) * Uff_RA(J, 1) + &
        dTaudy_RA(J, 2, 2) * Uff_RA(J, 2) + &
        dTaudy_RA(J, 2, 3) * Uff_RA(J, 3)
        Budg_viscs_accl1_MKE(J) = -1.0_WP * Budg_viscs_accl1_TKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_viscs_accl1_Duiuj(J, L) = dTaudy_RA(J, 2, M) * Uff_RA(J, N) + &
                dTaudy_RA(J, 2, N) * Uff_RA(J, M)
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_viscs_accl1_Duiuj', MYID)

    !==FA=====Production by driven force method 1 ======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Budg_prodc_Dvfc1_DuiuJ = 0.0_WP
    DO J = 1, NCL2

        !============exact method============================
        Budg_prodc_Dvfc1_Duiuj(J, 1) = ( FUxztL_F0_io(J, 1) -FUxztL_F0_io(J, 4) * U_FA(J, 1)  ) * 2.0_WP
        Budg_prodc_Dvfc1_Duiuj(J, 2) =  FUxztL_F0_io(J, 2) -FUxztL_F0_io(J, 4) * U_FA(J, 2)
        Budg_prodc_Dvfc1_Duiuj(J, 3) =  FUxztL_F0_io(J, 3) -FUxztL_F0_io(J, 4) * U_FA(J, 3)
        Budg_prodc_Dvfc1_TKE(J)   = FUxztL_F0_io(J, 1) - FUxztL_F0_io(J, 4) * U_FA(J, 1)
        !            Budg_prodc_Dvfc1_MKE(J) = FUxztL_F0_io(J, 4) * U_FA(J, 1)

        !============== Rough method==================
        !            Budg_prodc_Dvfc1_Duiuj(J, 1) =  DrivenForce(J) * Uff_RA(J, 1) * 2.0_WP
        !            Budg_prodc_Dvfc1_Duiuj(J, 2) =  DrivenForce(J) * Uff_RA(J, 2)
        !            Budg_prodc_Dvfc1_Duiuj(J, 3) =  DrivenForce(J) * Uff_RA(J, 3)
        !            Budg_prodc_Dvfc1_TKE(J)   =  DrivenForce(J) * Uff_RA(J, 1)
        Budg_prodc_Dvfc1_MKE(J)   =  DrivenForce(J) * U_FA(J, 1)   ! Good

    END DO

    !=========FA====== Method 1 decomPOSITION, no gravity production for turb. =================
    DO J = 1, NCL2
        Budg_prodc_gvfc1_MKE(J) = F_A* D1xztL_F0_io(J) * U1xztL_F0_io(J, 1)
    END DO

    !====FA==========below 2 second decomPOSITION method ===========================================
    !==FA=====Production by driven force method 2 ==
    DO J = 1, NCL2
        Budg_prodc_Dvfc2_Duiuj(J, 1) = ( FUxztL_F0_io(J, 1) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 1)  ) * 2.0_WP
        Budg_prodc_Dvfc2_Duiuj(J, 2) =  FUxztL_F0_io(J, 2) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 2)
        Budg_prodc_Dvfc2_Duiuj(J, 3) =  FUxztL_F0_io(J, 3) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 3)
        Budg_prodc_Dvfc2_TKE(J)   =  FUxztL_F0_io(J, 1) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 1)
        Budg_prodc_Dvfc2_MKE(J)   =  DrivenForce(J) * U1xztL_F0_io(J, 1)
    END DO

    !==FA=====Production by gravity force method 2 ==
    Budg_prodc_gvfc2_DuiuJ = 0.0_WP
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                IF(M /= ABS(iGravity) .AND. N /= ABS(iGravity)) CYCLE

                L = (M * (7-M)) / 2 + N - 3

                IF(M == ABS(iGravity) .AND. N == ABS(iGravity)) THEN
                    COE = 2.0_WP
                    K = ABS(iGravity)
                ELSE IF (M == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = N
                ELSE IF (N == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = M
                ELSE
                    COE = 0.0_WP
                    K = -1 !(WHICH will LEAD TO ERROR!)
                END IF
                ! F_A INCLUDEs a postive or negtive sign
                Budg_prodc_gvfc2_Duiuj(J, L) = F_A * COE * ( G1xztL_F0_io(J, K) - D1xztL_F0_io(J) * U1xztL_F0_io(J, K) )

                !du= ( G1xztL_F0_io(J, 1) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 1) ) * F_A
                !dV = ( G1xztL_F0_io(J, 2) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 2) ) * F_A
                !dw= ( G1xztL_F0_io(J, 3) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 3) ) * F_A

            END DO
        END DO
        !==============for TKE and MKE ==================
        Budg_prodc_gvfc2_TKE(J) = 0.5_WP * (Budg_prodc_gvfc2_Duiuj(J, 1) +Budg_prodc_gvfc2_Duiuj(J, 4) +Budg_prodc_gvfc2_Duiuj(J,6))

        Budg_prodc_gvfc2_MKE(J) = F_A* D1xztL_F0_io(J) * U_FA(J, 1)
    END DO

    !==FA=====Production by gravity force method 2 ==
    Budg_turss_accl2_DuiuJ = 0.0_WP
    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_turss_accl2_Duiuj(J, L) = uff_RA(J, M) * dTSSdy_RA(J, 2, N) + uff_RA(J, N) * dTSSdy_RA(J, 2, M)
            END DO
        END DO

        Budg_turss_accl2_TKE(J) =  uff_RA(J, 1) * dTSSdy_RA(J, 2, 1) + &
        uff_RA(J, 2) * dTSSdy_RA(J, 2, 2) + &
        uff_RA(J, 3) * dTSSdy_RA(J, 2, 3)
        Budg_turss_accl2_MKE(J) = -1.0_WP * Budg_turss_accl2_TKE(J)
    END DO

    !====FA==========below IS the 3 third  sum of some termS ===========================
    !========pressure related===================
    DO J = 1, NCL2
        Budg_pressure3_MKE(J) = -U_FA(J, 2) * DPDX_RA(J, 2)
        FF = 0.0_WP
        IF(J == 1) THEN
            FF = ( &
            ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J,  2)   &
            + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, 2) ) -&
            0.0_WP) * DYFI(J)

        ELSE IF (J == NCL2) THEN
            FF = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * UPxztL_F0_io(J,  2)   &
            + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, 2) ) &
            ) * DYFI(J)
        ELSE
            FF = ( &
            ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J,  2)   &
            + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, 2) ) -&
            ( YCL2ND_WFF(J) * UPxztL_F0_io(J,  2)   &
            + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, 2) ) &
            ) * DYFI(J)
        END IF
        Budg_pressure3_TKE(J) = FF - DVDLPxztL_F0_io(J, 1, 1) - DVDLPxztL_F0_io(J, 2, 2) - DVDLPxztL_F0_io(J, 3, 3) &
        - U_FA(J, 2) * DPDX_RA(J, 2)
        Budg_pressure3_TKE(J) = -1.0_WP * Budg_pressure3_TKE(J)


        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                FF = 0.0_WP
                IF(N == 2) THEN
                    IF(J == 1) THEN
                        FF = FF+ ( &
                        ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J,  M)   &
                        + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, M) ) -&
                        0.0_WP &
                        ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        FF = FF+ ( &
                        0.0_WP -&
                        ( YCL2ND_WFF(J) * UPxztL_F0_io(J,  M)   &
                        + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, M) ) &
                        ) * DYFI(J)
                    ELSE
                        FF = FF+ ( &
                        ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J,  M)   &
                        + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, M) ) -&
                        ( YCL2ND_WFF(J) * UPxztL_F0_io(J,  M)   &
                        + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, M) ) &
                        ) * DYFI(J)
                    END IF
                END IF

                IF(M == 2) THEN
                    IF(J == 1) THEN
                        FF = FF+  ( &
                        ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J, N)   &
                        + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, N) ) -&
                        0.0_WP &
                        ) * DYFI(J)
                    ELSE IF (J == NCL2) THEN
                        FF = FF+ ( &
                        0.0_WP -&
                        ( YCL2ND_WFF(J) * UPxztL_F0_io(J, N)   &
                        + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, N) ) &
                        ) * DYFI(J)
                    ELSE
                        FF = FF+  ( &
                        ( YCL2ND_WFF(J + 1) * UPxztL_F0_io(J, N)   &
                        + YCL2ND_WFB(J + 1) * UPxztL_F0_io(J + 1, N) ) -&
                        ( YCL2ND_WFF(J) * UPxztL_F0_io(J, N)   &
                        + YCL2ND_WFB(J) * UPxztL_F0_io(J - 1, N) ) &
                        ) * DYFI(J)
                    END IF
                END IF

                Budg_pressure3_Duiuj(J, L) = FF - DVDLPxztL_F0_io(J, M, N) - DVDLPxztL_F0_io(J, N, M) &
                - U_FA(J, M) * DPDX_RA(J, N) - U_FA(J, N) * DPDX_RA(J, M)
                Budg_pressure3_Duiuj(J, L) = -1.0_WP * Budg_pressure3_Duiuj(J, L)
            END DO
        END DO

    END DO
END IF


    !======== Stress related===================
    DO J = 1, NCL2
        Budg_vistress3_MKE(J) = U_FA(J, 1) * DTaudy_RA(J, 2, 1) + &
        U_FA(J, 2) * DTaudy_RA(J, 2, 2) + &
        U_FA(J, 3) * DTaudy_RA(J, 2, 3)
        FF = 0.0_WP
        IF(J == 1) THEN
            FF = ( &
            ( YCL2ND_WFF(J + 1) * (TauU_RA(J,  1, 2, 1) + TauU_RA(J,  2, 2, 2) + TauU_RA(J,  3, 2, 3))   &
            + YCL2ND_WFB(J + 1) * (TauU_RA(J + 1, 1, 2, 1) + TauU_RA(J + 1, 2, 2, 2) + TauU_RA(J + 1, 3, 2, 3)) ) -&
            0.0_WP &
            ) * DYFI(J)
        ELSE IF (J == NCL2) THEN
            FF = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * (TauU_RA(J,  1, 2, 1) + TauU_RA(J,  2, 2, 2) + TauU_RA(J,  3, 2, 3))   &
            + YCL2ND_WFB(J) * (TauU_RA(J - 1, 1, 2, 1) + TauU_RA(J - 1, 2, 2, 2) + TauU_RA(J - 1, 3, 2, 3)) ) &
            ) * DYFI(J)
        ELSE
            FF = ( &
            ( YCL2ND_WFF(J + 1) * (TauU_RA(J,  1, 2, 1) + TauU_RA(J,  2, 2, 2) + TauU_RA(J,  3, 2, 3))   &
            + YCL2ND_WFB(J + 1) * (TauU_RA(J + 1, 1, 2, 1) + TauU_RA(J + 1, 2, 2, 2) + TauU_RA(J + 1, 3, 2, 3)) ) -&
            ( YCL2ND_WFF(J) * (TauU_RA(J,  1, 2, 1) + TauU_RA(J,  2, 2, 2) + TauU_RA(J,  3, 2, 3))   &
            + YCL2ND_WFB(J) * (TauU_RA(J - 1, 1, 2, 1) + TauU_RA(J - 1, 2, 2, 2) + TauU_RA(J - 1, 3, 2, 3)) ) &
            ) * DYFI(J)
        END IF
        Budg_vistress3_TKE(J) = FF &
        - TauDvDL_RA(J, 1, 1, 1, 1) - TauDvDL_RA(J, 2, 1, 2, 1) - TauDvDL_RA(J, 3, 1, 3, 1) &
        - TauDvDL_RA(J, 1, 2, 1, 2) - TauDvDL_RA(J, 2, 2, 2, 2) - TauDvDL_RA(J, 3, 2, 3, 2) &
        - TauDvDL_RA(J, 1, 3, 1, 3) - TauDvDL_RA(J, 2, 3, 2, 3) - TauDvDL_RA(J, 3, 3, 3, 3) &
        - U_FA(J, 1) * DTaudy_RA(J, 2, 1) &
        - U_FA(J, 2) * DTaudy_RA(J, 2, 2) &
        - U_FA(J, 3) * DTaudy_RA(J, 2, 3)


        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                FF = 0.0_WP

                IF(J == 1) THEN
                    FF = ( &
                    ( YCL2ND_WFF(J + 1) * TauU_RA(J,  M, 2, N)   &
                    + YCL2ND_WFB(J + 1) * TauU_RA(J + 1, M, 2, N) ) -&
                    0.0_WP &
                    ) * DYFI(J) + ( &
                    ( YCL2ND_WFF(J + 1) * TauU_RA(J, N, 2, M)   &
                    + YCL2ND_WFB(J + 1) * TauU_RA(J + 1, N, 2, M) ) -&
                    0.0_WP &
                    ) * DYFI(J)
                ELSE IF (J == NCL2) THEN
                    FF = ( &
                    0.0_WP -&
                    ( YCL2ND_WFF(J) * TauU_RA(J,  M, 2, N)   &
                    + YCL2ND_WFB(J) * TauU_RA(J - 1, M, 2, N) ) &
                    ) * DYFI(J) + ( &
                    0.0_WP -&
                    ( YCL2ND_WFF(J) * TauU_RA(J, N, 2, M)   &
                    + YCL2ND_WFB(J) * TauU_RA(J - 1, N, 2, M) ) &
                    ) * DYFI(J)
                ELSE
                    FF = ( &
                    ( YCL2ND_WFF(J + 1) * TauU_RA(J,  M, 2, N)   &
                    + YCL2ND_WFB(J + 1) * TauU_RA(J + 1, M, 2, N) ) -&
                    ( YCL2ND_WFF(J) * TauU_RA(J,  M, 2, N)   &
                    + YCL2ND_WFB(J) * TauU_RA(J - 1, M, 2, N) ) &
                    ) * DYFI(J) +  ( &
                    ( YCL2ND_WFF(J + 1) * TauU_RA(J, N, 2, M)   &
                    + YCL2ND_WFB(J + 1) * TauU_RA(J + 1, N, 2, M) ) -&
                    ( YCL2ND_WFF(J) * TauU_RA(J, N, 2, M)   &
                    + YCL2ND_WFB(J) * TauU_RA(J - 1, N, 2, M) ) &
                    ) * DYFI(J)
                END IF




                Budg_vistress3_Duiuj(J, L) = FF &
                - TauDvDL_RA(J, M, 1, N, 1) - TauDvDL_RA(J, M, 2, N, 2) - TauDvDL_RA(J, M, 3, N, 3) &
                - TauDvDL_RA(J, N, 1, M, 1) - TauDvDL_RA(J, N, 2, M, 2) - TauDvDL_RA(J, N, 3, M, 3) &
                - U_FA(J, M) * DTaudy_RA(J, 2, N) - U_FA(J, N) * DTaudy_RA(J, 2, M)

                !IF (L == 1) THEN
                !!WRITE(*, *) FF, - TauDvDL_RA(J, M, 1, N, 1) - TauDvDL_RA(J, M, 2, N, 2) - TauDvDL_RA(J, M, 3, N, 3), &
                !- TauDvDL_RA(J, N, 1, M, 1) - TauDvDL_RA(J, N, 2, M, 2) - TauDvDL_RA(J, N, 3, M, 3), &
                !- U_FA(J, M) * DTaudy_RA(J, 2, N), - U_FA(J, N) * DTaudy_RA(J, 2, M)
                !!WRITE(*, *) - U_FA(J, M) * DTaudy_RA(J, 2, N), - U_FA(J, N) * DTaudy_RA(J, 2, M), U_FA(J, M), dTaudy_RA(J, 2, M)
                !END IF

            END DO
        END DO

    END DO

    !====FA===========BALANCE ===================
    !==============buoyancy force DIRect production IS zero======================
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_Balance1_TKE(J) =   Budg_prodc_stres_TKE(J) + &
        Budg_viscs_dissp_TKE(J) + &
        Budg_pduDX_stran_TKE(J) + &
        Budg_Turbu_diffu_TKE(J) + &
        Budg_DpuDX_diffu_TKE(J) + &
        Budg_viscs_diffu_TKE(J) + &
        Budg_press_accl1_TKE(J) + &
        Budg_viscs_accl1_TKE(J) + &
        Budg_prodc_Dvfc1_TKE(J)

        Budg_Balance2_TKE(J) =   Budg_prodc_stres_TKE(J) + &
        Budg_viscs_dissp_TKE(J) + &
        Budg_pduDX_stran_TKE(J) + &
        Budg_Turbu_diffu_TKE(J) + &
        Budg_DpuDX_diffu_TKE(J) + &
        Budg_viscs_diffu_TKE(J) + &
        Budg_turss_accl2_TKE(J) + &
        Budg_prodc_Dvfc2_TKE(J) + &
        Budg_prodc_gvfc2_TKE(J)

        Budg_Balance3_TKE(J) =   Budg_prodc_stres_TKE(J) + &
        Budg_Turbu_diffu_TKE(J) + &
        Budg_prodc_Dvfc1_TKE(J) + &
        Budg_pressure3_TKE(J) + &
        Budg_vistress3_TKE(J)

        Budg_Balance1_MKE(J) =   Budg_prodc_stres_MKE(J) + &
        Budg_viscs_dissp_MKE(J) + &
        Budg_pduDX_stran_MKE(J) + &
        Budg_Turbu_diffu_MKE(J) + &
        Budg_DpuDX_diffu_MKE(J) + &
        Budg_viscs_diffu_MKE(J) + &
        Budg_press_accl1_MKE(J) + &
        Budg_viscs_accl1_MKE(J) + &
        Budg_prodc_gvfc1_MKE(J) + &
        Budg_prodc_Dvfc1_MKE(J)

        Budg_Balance2_MKE(J) =   Budg_prodc_stres_MKE(J) + &
        Budg_viscs_dissp_MKE(J) + &
        Budg_pduDX_stran_MKE(J) + &
        Budg_Turbu_diffu_MKE(J) + &
        Budg_DpuDX_diffu_MKE(J) + &
        Budg_viscs_diffu_MKE(J) + &
        Budg_turss_accl2_MKE(J) + &
        Budg_prodc_gvfc2_MKE(J) + &
        Budg_prodc_Dvfc2_MKE(J)

        Budg_Balance3_MKE(J) =   Budg_prodc_stres_MKE(J) + &
        Budg_Turbu_diffu_MKE(J) + &
        Budg_prodc_Dvfc1_MKE(J) + &
        Budg_pressure3_MKE(J)  + &
        Budg_vistress3_MKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_Balance1_Duiuj(J, L) =Budg_prodc_stres_Duiuj(J, L) + &
                Budg_viscs_dissp_Duiuj(J, L) + &
                Budg_pduDX_stran_Duiuj(J, L) + &
                Budg_Turbu_diffu_Duiuj(J, L) + &
                Budg_DpuDX_diffu_Duiuj(J, L) + &
                Budg_viscs_diffu_Duiuj(J, L) + &
                Budg_press_accl1_Duiuj(J, L) + &
                Budg_viscs_accl1_Duiuj(J, L) + &
                Budg_prodc_Dvfc1_Duiuj(J, L)

                Budg_Balance2_Duiuj(J, L) =Budg_prodc_stres_Duiuj(J, L) + &
                Budg_viscs_dissp_Duiuj(J, L) + &
                Budg_pduDX_stran_Duiuj(J, L) + &
                Budg_Turbu_diffu_Duiuj(J, L) + &
                Budg_DpuDX_diffu_Duiuj(J, L) + &
                Budg_viscs_diffu_Duiuj(J, L) + &
                Budg_turss_accl2_Duiuj(J, L) + &
                Budg_prodc_gvfc2_Duiuj(J, L) + &
                Budg_prodc_Dvfc2_Duiuj(J, L)

                Budg_Balance3_Duiuj(J, L) =Budg_prodc_stres_Duiuj(J, L) + &
                Budg_Turbu_diffu_Duiuj(J, L) + &
                Budg_prodc_Dvfc1_Duiuj(J, L) + &
                Budg_pressure3_Duiuj(J, L) + &
                Budg_vistress3_Duiuj(J, L)


            END DO
        END DO
    END DO


    !================ INTegral of each terms along Y ==================================
    DO J = 1, NCL2
        !============ TKE ======================================================================
        Budg_prodc_stres_TKE_ysum = Budg_prodc_stres_TKE_ysum + Budg_prodc_stres_TKE(J) / DYFI(J)
        Budg_viscs_dissp_TKE_ysum = Budg_viscs_dissp_TKE_ysum + Budg_viscs_dissp_TKE(J) / DYFI(J)
        Budg_pduDX_stran_TKE_ysum = Budg_pduDX_stran_TKE_ysum + Budg_pduDX_stran_TKE(J) / DYFI(J)
        Budg_Turbu_diffu_TKE_ysum = Budg_Turbu_diffu_TKE_ysum + Budg_Turbu_diffu_TKE(J) / DYFI(J)
        Budg_DpuDX_diffu_TKE_ysum = Budg_DpuDX_diffu_TKE_ysum + Budg_DpuDX_diffu_TKE(J) / DYFI(J)
        Budg_viscs_diffu_TKE_ysum = Budg_viscs_diffu_TKE_ysum + Budg_viscs_diffu_TKE(J) / DYFI(J)

        Budg_press_accl1_TKE_ysum = Budg_press_accl1_TKE_ysum + Budg_press_accl1_TKE(J) / DYFI(J)
        Budg_viscs_accl1_TKE_ysum = Budg_viscs_accl1_TKE_ysum + Budg_viscs_accl1_TKE(J) / DYFI(J)
        Budg_prodc_Dvfc1_TKE_ysum = Budg_prodc_Dvfc1_TKE_ysum + Budg_prodc_Dvfc1_TKE(J) / DYFI(J)
        Budg_Balance1_TKE_ysum  = Budg_Balance1_TKE_ysum  + Budg_Balance1_TKE(J) / DYFI(J)

        Budg_turss_accl2_TKE_ysum = Budg_turss_accl2_TKE_ysum + Budg_turss_accl2_TKE(J) / DYFI(J)
        Budg_prodc_gvfc2_TKE_ysum = Budg_prodc_gvfc2_TKE_ysum + Budg_prodc_gvfc2_TKE(J) / DYFI(J)
        Budg_prodc_Dvfc2_TKE_ysum = Budg_prodc_Dvfc2_TKE_ysum + Budg_prodc_Dvfc2_TKE(J) / DYFI(J)
        Budg_Balance2_TKE_ysum  = Budg_Balance2_TKE_ysum  + Budg_Balance2_TKE(J) / DYFI(J)


        Budg_pressure3_TKE_ysum = Budg_pressure3_TKE_ysum + Budg_pressure3_TKE(J) / DYFI(J)
        Budg_vistress3_TKE_ysum = Budg_vistress3_TKE_ysum + Budg_vistress3_TKE(J) / DYFI(J)
        Budg_Balance3_TKE_ysum  = Budg_Balance3_TKE_ysum  + Budg_Balance3_TKE(J) / DYFI(J)


        !=========== MKE =======================================================================
        Budg_prodc_stres_MKE_ysum = Budg_prodc_stres_MKE_ysum + Budg_prodc_stres_MKE(J) / DYFI(J)
        Budg_viscs_dissp_MKE_ysum = Budg_viscs_dissp_MKE_ysum + Budg_viscs_dissp_MKE(J) / DYFI(J)
        Budg_pduDX_stran_MKE_ysum = Budg_pduDX_stran_MKE_ysum + Budg_pduDX_stran_MKE(J) / DYFI(J)
        Budg_Turbu_diffu_MKE_ysum = Budg_Turbu_diffu_MKE_ysum + Budg_Turbu_diffu_MKE(J) / DYFI(J)
        Budg_DpuDX_diffu_MKE_ysum = Budg_DpuDX_diffu_MKE_ysum + Budg_DpuDX_diffu_MKE(J) / DYFI(J)
        Budg_viscs_diffu_MKE_ysum = Budg_viscs_diffu_MKE_ysum + Budg_viscs_diffu_MKE(J) / DYFI(J)

        Budg_press_accl1_MKE_ysum = Budg_press_accl1_MKE_ysum + Budg_press_accl1_MKE(J) / DYFI(J)
        Budg_viscs_accl1_MKE_ysum = Budg_viscs_accl1_MKE_ysum + Budg_viscs_accl1_MKE(J) / DYFI(J)
        Budg_prodc_Dvfc1_MKE_ysum = Budg_prodc_Dvfc1_MKE_ysum + Budg_prodc_Dvfc1_MKE(J) / DYFI(J)
        Budg_prodc_gvfc1_MKE_ysum = Budg_prodc_gvfc1_MKE_ysum + Budg_prodc_gvfc1_MKE(J) / DYFI(J)
        Budg_Balance1_MKE_ysum  = Budg_Balance1_MKE_ysum  + Budg_Balance1_MKE(J) / DYFI(J)

        Budg_turss_accl2_MKE_ysum = Budg_turss_accl2_MKE_ysum + Budg_turss_accl2_MKE(J) / DYFI(J)
        Budg_prodc_gvfc2_MKE_ysum = Budg_prodc_gvfc2_MKE_ysum + Budg_prodc_gvfc2_MKE(J) / DYFI(J)
        Budg_prodc_Dvfc2_MKE_ysum = Budg_prodc_Dvfc2_MKE_ysum + Budg_prodc_Dvfc2_MKE(J) / DYFI(J)
        Budg_Balance2_MKE_ysum  = Budg_Balance2_MKE_ysum  + Budg_Balance2_MKE(J) / DYFI(J)


        Budg_pressure3_MKE_ysum = Budg_pressure3_MKE_ysum + Budg_pressure3_MKE(J) / DYFI(J)
        Budg_vistress3_MKE_ysum = Budg_vistress3_MKE_ysum + Budg_vistress3_MKE(J) / DYFI(J)
        Budg_Balance3_MKE_ysum  = Budg_Balance3_MKE_ysum  + Budg_Balance3_MKE(J) / DYFI(J)


        DO L = 1, (NDV * (7 - NDV)) / 2 + NDV - 3
            Budg_prodc_stres_Duiuj_ysum(L) = Budg_prodc_stres_Duiuj_ysum(L) + Budg_prodc_stres_Duiuj(J, L) / DYFI(J)
            Budg_viscs_dissp_Duiuj_ysum(L) = Budg_viscs_dissp_Duiuj_ysum(L) + Budg_viscs_dissp_Duiuj(J, L) / DYFI(J)
            Budg_pduDX_stran_Duiuj_ysum(L) = Budg_pduDX_stran_Duiuj_ysum(L) + Budg_pduDX_stran_Duiuj(J, L) / DYFI(J)
            Budg_Turbu_diffu_Duiuj_ysum(L) = Budg_Turbu_diffu_Duiuj_ysum(L) + Budg_Turbu_diffu_Duiuj(J, L) / DYFI(J)
            Budg_DpuDX_diffu_Duiuj_ysum(L) = Budg_DpuDX_diffu_Duiuj_ysum(L) + Budg_DpuDX_diffu_Duiuj(J, L) / DYFI(J)
            Budg_viscs_diffu_Duiuj_ysum(L) = Budg_viscs_diffu_Duiuj_ysum(L) + Budg_viscs_diffu_Duiuj(J, L) / DYFI(J)

            Budg_press_accl1_Duiuj_ysum(L) = Budg_press_accl1_Duiuj_ysum(L) + Budg_press_accl1_Duiuj(J, L) / DYFI(J)
            Budg_viscs_accl1_Duiuj_ysum(L) = Budg_viscs_accl1_Duiuj_ysum(L) + Budg_viscs_accl1_Duiuj(J, L) / DYFI(J)
            Budg_prodc_Dvfc1_Duiuj_ysum(L) = Budg_prodc_Dvfc1_Duiuj_ysum(L) + Budg_prodc_Dvfc1_Duiuj(J, L) / DYFI(J)
            Budg_Balance1_Duiuj_ysum(L) = Budg_Balance1_Duiuj_ysum(L)  + Budg_Balance1_Duiuj(J, L) / DYFI(J)

            Budg_turss_accl2_Duiuj_ysum(L) = Budg_turss_accl2_Duiuj_ysum(L) + Budg_turss_accl2_Duiuj(J, L) / DYFI(J)
            Budg_prodc_gvfc2_Duiuj_ysum(L) = Budg_prodc_gvfc2_Duiuj_ysum(L) + Budg_prodc_gvfc2_Duiuj(J, L) / DYFI(J)
            Budg_prodc_Dvfc2_Duiuj_ysum(L) = Budg_prodc_Dvfc2_Duiuj_ysum(L) + Budg_prodc_Dvfc2_Duiuj(J, L) / DYFI(J)
            Budg_Balance2_Duiuj_ysum(L) = Budg_Balance2_Duiuj_ysum(L)  + Budg_Balance2_Duiuj(J, L) / DYFI(J)


            Budg_pressure3_Duiuj_ysum(L) = Budg_pressure3_Duiuj_ysum(L) + Budg_pressure3_Duiuj(J, L) / DYFI(J)
            Budg_vistress3_Duiuj_ysum(L) = Budg_vistress3_Duiuj_ysum(L) + Budg_vistress3_Duiuj(J, L) / DYFI(J)
            Budg_Balance3_Duiuj_ysum(L) = Budg_Balance3_Duiuj_ysum(L)  + Budg_Balance3_Duiuj(J, L) / DYFI(J)
        END DO

    END DO

    !        !=============below for validatioN ======== (theoretICALLy zero) ===================
    !        sum_p_related_mkE = 0.0_WP
    !        DO J = 1, NCL2
    !            sum_p_related_mke = Budg_pduDX_stran_MKE(J) +Budg_DpuDX_diffu_MKE(J) +Budg_press_accl1_MKE(J)
    !            !WRITE(*, *) '##Sum of P -MKE:', J, sum_p_related_mke
    !        END DO
    !        !================

    !======== Sum of the pressure term related in MKE IS zero=========================
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.budget.check.' // TRIM(PNTIM) // '.plt'
    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist =File_exists)
    IF(File_exists) THEN
        OPEN(TECFLG, FILE =FLNM, POSITION = 'APPEND')
    ELSE
        OPEN(TECFLG, FILE =FLNM)
    END IF

    WRITE(TECFLG, '(A)') '=======FA=====uu,uv,uw, Vv, Vw,ww========== '
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_stres_Duiuj_ysum(1:6) = ', Budg_prodc_stres_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_dissp_Duiuj_ysum(1:6) = ', Budg_viscs_dissp_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_pduDX_stran_Duiuj_ysum(1:6) = ', Budg_pduDX_stran_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Turbu_diffu_Duiuj_ysum(1:6) = ', Budg_Turbu_diffu_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_DpuDX_diffu_Duiuj_ysum(1:6) = ', Budg_DpuDX_diffu_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_diffu_Duiuj_ysum(1:6) = ', Budg_viscs_diffu_Duiuj_ysum(1:6)

    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_press_accl1_Duiuj_ysum(1:6) = ', Budg_press_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_accl1_Duiuj_ysum(1:6) = ', Budg_viscs_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_Dvfc1_Duiuj_ysum(1:6) = ', Budg_prodc_Dvfc1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Balance1_Duiuj_ysum(1:6) = ', Budg_Balance1_Duiuj_ysum(1:6)

    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_turss_accl2_Duiuj_ysum(1:6) = ', Budg_turss_accl2_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_gvfc2_Duiuj_ysum(1:6) = ', Budg_prodc_gvfc2_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_Dvfc2_Duiuj_ysum(1:6) = ', Budg_prodc_Dvfc2_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Balance2_Duiuj_ysum(1:6) = ', Budg_Balance2_Duiuj_ysum(1:6)


    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_pressure3_Duiuj_ysum(1:6)(calc) = ',  Budg_pressure3_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_pressure3_Duiuj_ysum(1:6)(addt) = ',  Budg_pduDX_stran_Duiuj_ysum(1:6) + &
    Budg_DpuDX_diffu_Duiuj_ysum(1:6) + &
    Budg_press_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_vistress3_Duiuj_ysum(1:6)(calc) = ',  Budg_vistress3_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_vistress3_Duiuj_ysum(1:6)(addt) = ',  Budg_viscs_dissp_Duiuj_ysum(1:6) + &
    Budg_viscs_diffu_Duiuj_ysum(1:6) + &
    Budg_viscs_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Balance3_Duiuj_ysum(1:6)      = ',  Budg_Balance3_Duiuj_ysum(1:6)


    WRITE(TECFLG, '(A)') '=======FA===== TKE, MKE ========== '
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_prodc_stres_ysum_TKE_and_MKE = ', Budg_prodc_stres_TKE_ysum, Budg_prodc_stres_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_viscs_dissp_ysum_TKE_and_MKE = ', Budg_viscs_dissp_TKE_ysum, Budg_viscs_dissp_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_pduDX_stran_ysum_TKE_and_MKE = ', Budg_pduDX_stran_TKE_ysum, Budg_pduDX_stran_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_Turbu_diffu_ysum_TKE_and_MKE = ', Budg_Turbu_diffu_TKE_ysum, Budg_Turbu_diffu_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_DpuDX_diffu_ysum_TKE_and_MKE = ', Budg_DpuDX_diffu_TKE_ysum, Budg_DpuDX_diffu_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_viscs_diffu_ysum_TKE_and_MKE = ', Budg_viscs_diffu_TKE_ysum, Budg_viscs_diffu_MKE_ysum

    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_press_accl1_ysum_TKE_and_MKE = ', Budg_press_accl1_TKE_ysum, Budg_press_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_viscs_accl1_ysum_TKE_and_MKE = ', Budg_viscs_accl1_TKE_ysum, Budg_viscs_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_prodc_Dvfc1_ysum_TKE_and_MKE = ', Budg_prodc_Dvfc1_TKE_ysum, Budg_prodc_Dvfc1_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_Balance1_ysum_TKE_and_MKE  = ', Budg_Balance1_TKE_ysum,    Budg_Balance1_MKE_ysum

    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_turss_accl2_ysum_TKE_and_MKE = ', Budg_turss_accl2_TKE_ysum, Budg_turss_accl2_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_prodc_gvfc2_ysum_TKE_and_MKE = ', Budg_prodc_gvfc2_TKE_ysum, Budg_prodc_gvfc2_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_prodc_Dvfc2_ysum_TKE_and_MKE = ', Budg_prodc_Dvfc2_TKE_ysum, Budg_prodc_Dvfc2_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_Balance2_ysum_TKE_and_MKE  = ', Budg_Balance2_TKE_ysum,    Budg_Balance2_MKE_ysum


    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_pressure3_ysum_TKE_and_MKE(calc) = ', Budg_pressure3_TKE_ysum,Budg_pressure3_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_pressure3_ysum_TKE_and_MKE(addt) = ', Budg_pduDX_stran_TKE_ysum + &
    Budg_DpuDX_diffu_TKE_ysum + &
    Budg_press_accl1_TKE_ysum,  &
    Budg_pduDX_stran_MKE_ysum + &
    Budg_DpuDX_diffu_MKE_ysum + &
    Budg_press_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_vistress3_ysum_TKE_and_MKE(calc) = ', Budg_vistress3_TKE_ysum, Budg_vistress3_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_vistress3_ysum_TKE_and_MKE(addt) = ', Budg_viscs_dissp_TKE_ysum + &
    Budg_viscs_diffu_TKE_ysum + &
    Budg_viscs_accl1_TKE_ysum, &
    Budg_viscs_dissp_MKE_ysum + &
    Budg_viscs_diffu_MKE_ysum + &
    Budg_viscs_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 2ES20.7)') 'Budg_Balance3_ysum_TKE_and_MKE      = ', Budg_Balance3_TKE_ysum, Budg_Balance3_MKE_ysum

    CLOSE(TECFLG)



    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
!###########################  RA   ##############################################################
SUBROUTINE PP_FLOW_RA_noDen_RSTE_Budg_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, H, L, K
    REAL(WP) :: Budg_prod_Pij1, Budg_prod_Pij2, COE
    INTEGER(4) :: TECflg = 200
    CHARACTER(128) :: FLNM
    LOGICAL :: File_exists


    CALL CHKHDL('=====Calculating RA BudegtS ===== ', MYID)
    CALL PP_Budg_INIT
    !=======PRODUCTION TERMS due to mean sheAR === (RS) =========================================================
    !       Eq: Budg_prodc_stres_Duiuj(J, L) = p *_{ij}=P_{ij}-2 / 3P\Delta_{ij}
    !       Eq: P_{ij} = -<\rho u''_i u''_k> (\partial {u_j} / \partial x_k) +
    !                    -<\rho u''_j u''_k> (\partial {u_i} / \partial x_k)
    Budg_prodc_stres_DuiuJ = 0.0_WP
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_prodc_stres_MKE(J) = uf2_RA(J, 1, 1) * DVDL1xztL_F0_io(J, 1, 1) + &
        uf2_RA(J, 1, 2) * DVDL1xztL_F0_io(J, 1, 2) + &
        uf2_RA(J, 1, 3) * DVDL1xztL_F0_io(J, 1, 3) + &
        uf2_RA(J, 2, 1) * DVDL1xztL_F0_io(J, 2, 1) + &
        uf2_RA(J, 2, 2) * DVDL1xztL_F0_io(J, 2, 2) + &
        uf2_RA(J, 2, 3) * DVDL1xztL_F0_io(J, 2, 3) + &
        uf2_RA(J, 3, 1) * DVDL1xztL_F0_io(J, 3, 1) + &
        uf2_RA(J, 3, 2) * DVDL1xztL_F0_io(J, 3, 2) + &
        uf2_RA(J, 3, 3) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_prodc_stres_TKE(J) = -1.0_WP *Budg_prodc_stres_MKE(J)

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                Budg_prod_Pij1 = uf2_RA(J, 1, M) * DVDL1xztL_F0_io(J, N, 1) + &
                uf2_RA(J, 2, M) * DVDL1xztL_F0_io(J, N, 2) + &
                uf2_RA(J, 3, M) * DVDL1xztL_F0_io(J, N, 3)
                Budg_prod_Pij2 = uf2_RA(J, 1, N) * DVDL1xztL_F0_io(J, M, 1) + &
                uf2_RA(J, 2, N) * DVDL1xztL_F0_io(J, M, 2) + &
                uf2_RA(J, 3, N) * DVDL1xztL_F0_io(J, M, 3)
                Budg_prodc_stres_Duiuj(J, L) = -1.0_WP * (Budg_prod_Pij1 + Budg_prod_Pij2)

            END DO
        END DO

    END DO
    CALL CHKHDL(' ==>Calculated Budg_prodc_stres_Duiuj', MYID)

    !== RA===== VIScous ENERGY dISsipation term ===========================================================
    !       viscstressVeloGrad_RA(J, M, N,H,P) = <\partial(u_m) /\partial x_n \tau_hp>
    Budg_viscs_dissp_TKE = 0.0_WP
    Budg_viscs_dissp_MKE = 0.0_WP
    Budg_viscs_dissp_DuiuJ = 0.0_WP
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_viscs_dissp_TKE(J) = &
        DVDL2xztL_F0_io(J, (1- 1) * 3 + 1, (1- 1) * 3 + 1) - DVDL1xztL_F0_io(J, 1, 1) * DVDL1xztL_F0_io(J, 1, 1) &
        + DVDL2xztL_F0_io(J, (1- 1) * 3 + 2, (1- 1) * 3 + 2) - DVDL1xztL_F0_io(J, 1, 2) * DVDL1xztL_F0_io(J, 1, 2) &
        + DVDL2xztL_F0_io(J, (1- 1) * 3 +3, (1- 1) * 3 +3) - DVDL1xztL_F0_io(J, 1, 3) * DVDL1xztL_F0_io(J, 1, 3) &
        + DVDL2xztL_F0_io(J, (2 - 1) * 3 + 1, (2 - 1) * 3 + 1) - DVDL1xztL_F0_io(J, 2, 1) * DVDL1xztL_F0_io(J, 2, 1) &
        + DVDL2xztL_F0_io(J, (2 - 1) * 3 + 2, (2 - 1) * 3 + 2) - DVDL1xztL_F0_io(J, 2, 2) * DVDL1xztL_F0_io(J, 2, 2) &
        + DVDL2xztL_F0_io(J, (2 - 1) * 3 +3, (2 - 1) * 3 +3) - DVDL1xztL_F0_io(J, 2, 3) * DVDL1xztL_F0_io(J, 2, 3) &
        + DVDL2xztL_F0_io(J, (3- 1) * 3 + 1, (3- 1) * 3 + 1) - DVDL1xztL_F0_io(J, 3, 1) * DVDL1xztL_F0_io(J, 3, 1) &
        + DVDL2xztL_F0_io(J, (3- 1) * 3 + 2, (3- 1) * 3 + 2) - DVDL1xztL_F0_io(J, 3, 2) * DVDL1xztL_F0_io(J, 3, 2) &
        + DVDL2xztL_F0_io(J, (3- 1) * 3 +3, (3- 1) * 3 +3) - DVDL1xztL_F0_io(J, 3, 3) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_viscs_dissp_TKE(J) = Budg_viscs_dissp_TKE(J) * CVISC * (-1.0_WP)

        Budg_viscs_dissp_MKE(J) = &
        Tau_Mean_RA(J, 1, 1) * DVDL1xztL_F0_io(J, 1, 1) + &
        Tau_Mean_RA(J, 2, 1) * DVDL1xztL_F0_io(J, 2, 1) + &
        Tau_Mean_RA(J, 3, 1) * DVDL1xztL_F0_io(J, 3, 1) + &
        Tau_Mean_RA(J, 1, 2) * DVDL1xztL_F0_io(J, 1, 2) + &
        Tau_Mean_RA(J, 2, 2) * DVDL1xztL_F0_io(J, 2, 2) + &
        Tau_Mean_RA(J, 3, 2) * DVDL1xztL_F0_io(J, 3, 2) + &
        Tau_Mean_RA(J, 1, 3) * DVDL1xztL_F0_io(J, 1, 3) + &
        Tau_Mean_RA(J, 2, 3) * DVDL1xztL_F0_io(J, 2, 3) + &
        Tau_Mean_RA(J, 3, 3) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_viscs_dissp_MKE(J) = -1.0_WP *Budg_viscs_dissp_MKE(J)

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                Budg_viscs_dissp_Duiuj(J, L) =  &
                DVDL2xztL_F0_io(J, (M - 1) * 3 + 1, (N - 1) * 3 + 1) - &
                DVDL1xztL_F0_io(J, M, 1) * DVDL1xztL_F0_io(J, N, 1)             &
                + DVDL2xztL_F0_io(J, (M - 1) * 3 + 2, (N - 1) * 3 + 2) - &
                DVDL1xztL_F0_io(J, M, 2) * DVDL1xztL_F0_io(J, N, 2)             &
                + DVDL2xztL_F0_io(J, (M - 1) * 3 +3, (N - 1) * 3 +3) - &
                DVDL1xztL_F0_io(J, M, 3) * DVDL1xztL_F0_io(J, N, 3)

                Budg_viscs_dissp_Duiuj(J, L) = Budg_viscs_dissp_Duiuj(J, L) * CVISC * (-2.0_WP)
                !IF(DABS(dUDX_RA(J, M, N) - DVDL1xztL_F0_io(J, M, N)) >  1.0E-12_WP) &
                !WRITE(*, '(3I4.1, 3ES13.5)') J, M, N, dUDX_RA(J, M, N), &
                !DVDL1xztL_F0_io(J, M, N), DUDX_RA(J, M, N) - DVDL1xztL_F0_io(J, M, N)

            END DO
        END DO
    END DO

    CALL CHKHDL(' ==>Calculated Budg_viscs_dissp_Duiuj', MYID)


    !== RA===== Velocity- PRessure gradient pressure strAIn terM ======= (RS) ==============================
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_pduDX_stran_TKE(J) =  DVDLPxztL_F0_io(J, 1, 1) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 1, 1) + &
        DVDLPxztL_F0_io(J, 2, 2) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 2, 2) + &
        DVDLPxztL_F0_io(J, 3, 3) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, 3, 3)
        Budg_pduDX_stran_MKE(J) =  U1xztL_F0_io(J, 4) * dUiDXi(J)

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_pduDX_stran_Duiuj(J, L) = DVDLPxztL_F0_io(J, M, N) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, M, N) + &
                DVDLPxztL_F0_io(J, N, M) - U1xztL_F0_io(J, 4) * DVDL1xztL_F0_io(J, N, M)
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_pduDX_stran_Duiuj', MYID)


    !== RA===== TURBUELCEN DIFFUSION TERMS = (turblence tranport rate) === (RS) ========================
    !       Eq.  Budg_TdIFf_Duiuj(J, L) = ( \partial <\rho u''_i u''_j u''_k > ) / (\partial x_k )
    DO J = 1, NCL2
        ufMKEfd_RA(J) =  0.5_WP * ( U3xztL_F0_io(J, 2) + U3xztL_F0_io(J,7) + U3xztL_F0_io(J,9) )

        ufTKEfd_RA(J) = 0.5_WP * ( uf3_RA(J, 1, 1, 2) + uf3_RA(J, 2, 2, 2) + uf3_RA(J, 3, 3, 2) )

    END DO


    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_Turbu_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UfTKEfd_RA(J + 1) + &
            YCL2ND_WFB(J + 1) * UfTKEfd_RA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UfMKEfd_RA(J + 1) + &
            YCL2ND_WFB(J + 1) * UfMKEfd_RA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
        ELSE IF (J == NCL2) THEN
            Budg_Turbu_diffu_TKE(J) = &
            ( 0.0_WP - &
            ( YCL2ND_WFF(J) * UfTKEfd_RA(J  ) + &
            YCL2ND_WFB(J) * UfTKEfd_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( 0.0_WP - &
            ( YCL2ND_WFF(J) * UfMKEfd_RA(J  ) + &
            YCL2ND_WFB(J) * UfMKEfd_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)

        ELSE
            Budg_Turbu_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UfTKEfd_RA(J + 1) + &
            YCL2ND_WFB(J + 1) * UfTKEfd_RA(J  ) ) -         &
            ( YCL2ND_WFF(J) * UfTKEfd_RA(J  ) + &
            YCL2ND_WFB(J) * UfTKEfd_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
            Budg_Turbu_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * UfMKEfd_RA(J + 1) + &
            YCL2ND_WFB(J + 1) * UfMKEfd_RA(J  ) ) -         &
            ( YCL2ND_WFF(J) * UfMKEfd_RA(J  ) + &
            YCL2ND_WFB(J) * UfMKEfd_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
        END IF


        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                IF(J == 1) THEN
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * Uf3_RA(J + 1, M, N, 2) + &
                    YCL2ND_WFB(J + 1) * Uf3_RA(J,  M, N, 2) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
                ELSE IF (J == NCL2) THEN
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( 0.0_WP - &
                    ( YCL2ND_WFF(J) * Uf3_RA(J,  M, N, 2) + &
                    YCL2ND_WFB(J) * Uf3_RA(J - 1, M, N, 2) ) ) * DYFI(J) * (-1.0_WP)

                ELSE
                    Budg_Turbu_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * Uf3_RA(J + 1, M, N, 2) + &
                    YCL2ND_WFB(J + 1) * Uf3_RA(J,  M, N, 2) ) -         &
                    ( YCL2ND_WFF(J) * Uf3_RA(J,  M, N, 2) + &
                    YCL2ND_WFB(J) * Uf3_RA(J - 1, M, N, 2) ) ) * DYFI(J) * (-1.0_WP)
                END IF

            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_Turbu_diffu_Duiuj', MYID)

    !== RA===== Velocity- PRessure gradient dIFfusion terM ======= (RS) ==============================
    !       Eq. = - \partial (<p' u''_j>) /\partial (x_i) - \partial (<p' u''_i>) /\partial (x_j)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_DpuDX_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * Ufpf_RA(J + 1, 2) + &
            YCL2ND_WFB(J + 1) * Ufpf_RA(J,  2) ) &
            -  0.0_WP ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * (U1xztL_F0_io(J + 1, 2) * U1xztL_F0_io(J + 1, 4)) + &
            YCL2ND_WFB(J + 1) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) ) &
            -  0.0_WP ) * DYFI(J) * (-1.0_WP)

        ELSE IF(J == NCL2) THEN
            Budg_DpuDX_diffu_TKE(J) = &
            ( 0.0_WP -  &
            ( YCL2ND_WFF(J) * Ufpf_RA(J,  2) + &
            YCL2ND_WFB(J) * Ufpf_RA(J - 1, 2) ) ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( 0.0_WP -  &
            ( YCL2ND_WFF(J) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) + &
            YCL2ND_WFB(J) * (U1xztL_F0_io(J - 1, 2) * U1xztL_F0_io(J - 1, 4)) )) * DYFI(J) * (-1.0_WP)
        ELSE

            Budg_DpuDX_diffu_TKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * Ufpf_RA(J + 1, 2) + &
            YCL2ND_WFB(J + 1) * Ufpf_RA(J,  2) ) -  &
            ( YCL2ND_WFF(J) * Ufpf_RA(J,  2) + &
            YCL2ND_WFB(J) * Ufpf_RA(J - 1, 2) ) ) * DYFI(J) * (-1.0_WP)

            Budg_DpuDX_diffu_MKE(J) = &
            ( ( YCL2ND_WFF(J + 1) * (U1xztL_F0_io(J + 1, 2) * U1xztL_F0_io(J + 1, 4)) + &
            YCL2ND_WFB(J + 1) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) ) -  &
            ( YCL2ND_WFF(J) * (U1xztL_F0_io(J,  2) * U1xztL_F0_io(J,  4)) + &
            YCL2ND_WFB(J) * (U1xztL_F0_io(J - 1, 2) * U1xztL_F0_io(J - 1, 4)) ) ) * DYFI(J) * (-1.0_WP)
        END IF

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                IF(J == 1) THEN
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * ( ufpf_RA(J + 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J + 1, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J + 1) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) ) - &
                    0.0_WP ) * DYFI(J)
                ELSE IF(J == NCL2) THEN
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    (  0.0_WP - &
                    ( YCL2ND_WFF(J) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J) * ( ufpf_RA(J - 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J - 1, N) * DBLE(Kronecker_Delta(M, 2))   ) ) ) * DYFI(J)
                ELSE
                    Budg_DpuDX_diffu_Duiuj(J, L) = &
                    ( ( YCL2ND_WFF(J + 1) * ( ufpf_RA(J + 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J + 1, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J + 1) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) ) - &
                    ( YCL2ND_WFF(J) * ( ufpf_RA(J,  M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J, N) * DBLE(Kronecker_Delta(M, 2))   ) +   &
                    YCL2ND_WFB(J) * ( ufpf_RA(J - 1, M) * DBLE(Kronecker_Delta(N, 2)) +      &
                    ufpf_RA(J - 1, N) * DBLE(Kronecker_Delta(M, 2))   ) ) ) * DYFI(J)
                END IF


                Budg_DpuDX_diffu_Duiuj(J, L) = Budg_DpuDX_diffu_Duiuj(J, L) * (-1.0_WP)
            END DO
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated Budg_DpuDX_diffu_Duiuj', MYID)


    !== RA===== VIScous dIFfusion term ===========================================================
    !       Eq. = \partial <u''_j tau'_ki> / pARtial (x_k) +
    !             \partial <u''_i tau'_kj> / pARtial (x_k)
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        IF(J == 1) THEN
            Budg_viscs_diffu_TKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J + 1) * ( Taufuf_RA(J + 1, 1, 2, 1) + &
            Taufuf_RA(J + 1, 2, 2, 2) + &
            Taufuf_RA(J + 1, 3, 2, 3) ) ) -&
            0.0_WP &
            ) * DYFI(J)
            Budg_viscs_diffu_MKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J + 1) * ( U1xztL_F0_io(J + 1, 1) * Tau_Mean_RA(J + 1, 1, 2) + &
            U1xztL_F0_io(J + 1, 2) * Tau_Mean_RA(J + 1, 2, 2) + &
            U1xztL_F0_io(J + 1, 3) * Tau_Mean_RA(J + 1, 3, 2) ) ) -&
            0.0_WP &
            ) * DYFI(J)
        ELSE IF (J == NCL2) THEN
            Budg_viscs_diffu_TKE(J) = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J) * ( Taufuf_RA(J - 1, 1, 2, 1) + &
            Taufuf_RA(J - 1, 2, 2, 2) + &
            Taufuf_RA(J - 1, 3, 2, 3) ) ) &
            ) * DYFI(J)

            Budg_viscs_diffu_MKE(J) = ( &
            0.0_WP -&
            ( YCL2ND_WFF(J) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J, 1,  2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J, 2,  2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J, 3,  2) ) &
            + YCL2ND_WFB(J) * ( U1xztL_F0_io(J - 1, 1) * Tau_Mean_RA(J - 1, 1, 2) + &
            U1xztL_F0_io(J - 1, 2) * Tau_Mean_RA(J - 1, 2, 2) + &
            U1xztL_F0_io(J - 1, 3) * Tau_Mean_RA(J - 1, 3, 2) ) ) &
            ) * DYFI(J)
        ELSE
            Budg_viscs_diffu_TKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J + 1) * ( Taufuf_RA(J + 1, 1, 2, 1) + &
            Taufuf_RA(J + 1, 2, 2, 2) + &
            Taufuf_RA(J + 1, 3, 2, 3) ) ) -&
            ( YCL2ND_WFF(J) * ( Taufuf_RA(J,  1, 2, 1) + &
            Taufuf_RA(J,  2, 2, 2) + &
            Taufuf_RA(J,  3, 2, 3) ) &
            + YCL2ND_WFB(J) * ( Taufuf_RA(J - 1, 1, 2, 1) + &
            Taufuf_RA(J - 1, 2, 2, 2) + &
            Taufuf_RA(J - 1, 3, 2, 3) ) ) &
            ) * DYFI(J)

            Budg_viscs_diffu_MKE(J) = ( &
            ( YCL2ND_WFF(J + 1) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J + 1) * ( U1xztL_F0_io(J + 1, 1) * Tau_Mean_RA(J + 1, 1, 2) + &
            U1xztL_F0_io(J + 1, 2) * Tau_Mean_RA(J + 1, 2, 2) + &
            U1xztL_F0_io(J + 1, 3) * Tau_Mean_RA(J + 1, 3, 2) ) ) -&
            ( YCL2ND_WFF(J) * ( U1xztL_F0_io(J,  1) * Tau_Mean_RA(J,  1, 2) + &
            U1xztL_F0_io(J,  2) * Tau_Mean_RA(J,  2, 2) + &
            U1xztL_F0_io(J,  3) * Tau_Mean_RA(J,  3, 2) ) &
            + YCL2ND_WFB(J) * ( U1xztL_F0_io(J - 1, 1) * Tau_Mean_RA(J - 1, 1, 2) + &
            U1xztL_F0_io(J - 1, 2) * Tau_Mean_RA(J - 1, 2, 2) + &
            U1xztL_F0_io(J - 1, 3) * Tau_Mean_RA(J - 1, 3, 2) ) ) &
            ) * DYFI(J)
        END IF

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3

                IF(J == 1) THEN
                    Budg_viscs_diffu_Duiuj(J, L) = uf2_RA(J + 1, M, N) * APVR(J, 1) + &
                    uf2_RA(J,  M, N) * ACVR(J, 1) + 0.0_WP * AMVR(J, 1)
                ELSE IF(J == NCL2) THEN
                    Budg_viscs_diffu_Duiuj(J, L) = 0.0_WP * APVR(J, 1) + &
                    uf2_RA(J,  M, N) * ACVR(J, 1) + &
                    uf2_RA(J - 1, M, N) * AMVR(J, 1)
                ELSE
                    Budg_viscs_diffu_Duiuj(J, L) = uf2_RA(J + 1, M, N) * APVR(J, 1) + &
                    uf2_RA(J,  M, N) * ACVR(J, 1) + &
                    uf2_RA(J - 1, M, N) * AMVR(J, 1)
                END IF
                Budg_viscs_diffu_Duiuj(J, L) = Budg_viscs_diffu_Duiuj(J, L) * CVISC

            END DO
        END DO

        Budg_viscs_diffu_TKE(J) = 0.5_WP * (Budg_viscs_diffu_Duiuj(J, 1) + &
        Budg_viscs_diffu_Duiuj(J, 4) + &
        Budg_viscs_diffu_Duiuj(J,6) )

    END DO
    CALL CHKHDL(' ==>Calculated Budg_viscs_diffu_Duiuj', MYID)
    IF(iThermoDynamics == 1) THEN
    !================BODY FORCE/ BUOYANCY PRODUCTION, whICh IS pARt of pressure acceleration terM ===========
    !       ! not an independe contribution, but belongs to pARt of Budg_press_accl1_Duiuj
    Budg_prodc_gvfc2_DuiuJ = 0.0_WP
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                IF(M /= ABS(iGravity) .AND. N /= ABS(iGravity)) CYCLE

                L = (M * (7-M)) / 2 + N - 3

                IF(M == ABS(iGravity) .AND. N == ABS(iGravity)) THEN
                    COE = 2.0_WP
                    K = ABS(iGravity)
                ELSE IF (M == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = N
                ELSE IF (N == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = M
                ELSE
                    COE = 0.0_WP
                    K = -1 !(WHICH will LEAD TO ERROR!)
                END IF
                ! F_A INCLUDEs a postive or negtive sign
                Budg_prodc_gvfc2_Duiuj(J, L) = F_A * COE * ( G1xztL_F0_io(J, K) - D1xztL_F0_io(J) * U1xztL_F0_io(J, K) )

                !du= ( G1xztL_F0_io(J, 1) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 1) ) * F_A
                !dV = ( G1xztL_F0_io(J, 2) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 2) ) * F_A
                !dw= ( G1xztL_F0_io(J, 3) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 3) ) * F_A

            END DO
        END DO

    END DO
    Budg_press_accl1_Duiuj(:, :) = 0.0_WP
    Budg_viscs_accl1_Duiuj(:, :) = 0.0_WP


    Budg_prodc_Dvfc1_DuiuJ = 0.0_WP
    DO J = 1, NCL2
        ! f_1 * u"_1 + f_1 * u"_1
        Budg_prodc_Dvfc1_Duiuj(J, 1) = ( FUxztL_F0_io(J, 1) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 1)  ) * 2.0_WP
        ! f_1 * u"_2
        Budg_prodc_Dvfc1_Duiuj(J, 2) =  FUxztL_F0_io(J, 2) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 2)
        ! f_1 * u"_3
        Budg_prodc_Dvfc1_Duiuj(J, 3) =  FUxztL_F0_io(J, 3) -FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 3)

        !==============for TKE and MKE ==================
        !Budg_prodc_Dvfc1_TKE(J) = DrivenForce(J) * Uff_RA(J, 1)
        !Budg_prodc_Dvfc1_MKE(J) = DrivenForce(J) * U_FA(J, 1)

        !!WRITE(*, *) 'Driven force', DrivenForce(J), FUxztL_F0_io(J, 4)
        !!WRITE(*, *) 'FU', FUxztL_F0_io(J, 1:4) !test

        Budg_prodc_Dvfc1_TKE(J) = FUxztL_F0_io(J, 1) - FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 1)
        !Budg_prodc_Dvfc1_MKE(J) = FUxztL_F0_io(J, 4) * U1xztL_F0_io(J, 1)
        Budg_prodc_Dvfc1_MKE(J) = DrivenForce(J) * U_FA(J, 1)

    END DO
END IF
    !==== RA===========BALANCE ===================
    DO J = 1, NCL2
        !==============for TKE and MKE ==================
        Budg_Balance1_TKE(J) =    Budg_prodc_stres_TKE(J) + &
        Budg_viscs_dissp_TKE(J) + &
        Budg_pduDX_stran_TKE(J) + &
        Budg_Turbu_diffu_TKE(J) + &
        Budg_DpuDX_diffu_TKE(J) + &
        Budg_viscs_diffu_TKE(J) + &
        Budg_press_accl1_TKE(J) + &
        Budg_viscs_accl1_TKE(J) + &
        Budg_prodc_gvfc2_TKE(J) + &
        Budg_prodc_Dvfc1_TKE(J)

        Budg_Balance1_MKE(J) =    Budg_prodc_stres_MKE(J) + &
        Budg_viscs_dissp_MKE(J) + &
        Budg_pduDX_stran_MKE(J) + &
        Budg_Turbu_diffu_MKE(J) + &
        Budg_DpuDX_diffu_MKE(J) + &
        Budg_viscs_diffu_MKE(J) + &
        Budg_press_accl1_MKE(J) + &
        Budg_viscs_accl1_MKE(J) + &
        Budg_prodc_gvfc1_MKE(J) + &
        Budg_prodc_Dvfc1_MKE(J)

        !==========for each RuV =========================
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                Budg_Balance1_Duiuj(J, L) = Budg_prodc_stres_Duiuj(J, L) + &
                Budg_viscs_dissp_Duiuj(J, L) + &
                Budg_pduDX_stran_Duiuj(J, L) + &
                Budg_Turbu_diffu_Duiuj(J, L) + &
                Budg_DpuDX_diffu_Duiuj(J, L) + &
                Budg_pduDX_stran_Duiuj(J, L) + &
                Budg_viscs_diffu_Duiuj(J, L) + &
                Budg_press_accl1_Duiuj(J, L) + &
                Budg_viscs_accl1_Duiuj(J, L) + &
                Budg_prodc_gvfc2_Duiuj(J, L) + &
                Budg_prodc_Dvfc1_Duiuj(J, L)


            END DO
        END DO
    END DO
    !================ INTegral of each terms along Y ==================================
    DO J = 1, NCL2
        Budg_prodc_stres_TKE_ysum = Budg_prodc_stres_TKE_ysum + Budg_prodc_stres_TKE(J) / DYFI(J)
        Budg_viscs_dissp_TKE_ysum = Budg_viscs_dissp_TKE_ysum + Budg_viscs_dissp_TKE(J) / DYFI(J)
        Budg_pduDX_stran_TKE_ysum = Budg_pduDX_stran_TKE_ysum + Budg_pduDX_stran_TKE(J) / DYFI(J)
        Budg_Turbu_diffu_TKE_ysum = Budg_Turbu_diffu_TKE_ysum + Budg_Turbu_diffu_TKE(J) / DYFI(J)
        Budg_DpuDX_diffu_TKE_ysum = Budg_DpuDX_diffu_TKE_ysum + Budg_DpuDX_diffu_TKE(J) / DYFI(J)
        Budg_viscs_diffu_TKE_ysum = Budg_viscs_diffu_TKE_ysum + Budg_viscs_diffu_TKE(J) / DYFI(J)
        Budg_press_accl1_TKE_ysum = Budg_press_accl1_TKE_ysum + Budg_press_accl1_TKE(J) / DYFI(J)
        Budg_viscs_accl1_TKE_ysum = Budg_viscs_accl1_TKE_ysum + Budg_viscs_accl1_TKE(J) / DYFI(J)
        Budg_prodc_gvfc2_TKE_ysum = Budg_prodc_gvfc2_TKE_ysum + Budg_prodc_gvfc2_TKE(J) / DYFI(J)
        Budg_prodc_Dvfc1_TKE_ysum = Budg_prodc_Dvfc1_TKE_ysum + Budg_prodc_Dvfc1_TKE(J) / DYFI(J)
        Budg_Balance1_TKE_ysum   = Budg_Balance1_TKE_ysum   + Budg_Balance1_TKE(J) / DYFI(J)

        Budg_prodc_stres_MKE_ysum = Budg_prodc_stres_MKE_ysum + Budg_prodc_stres_MKE(J) / DYFI(J)
        Budg_viscs_dissp_MKE_ysum = Budg_viscs_dissp_MKE_ysum + Budg_viscs_dissp_MKE(J) / DYFI(J)
        Budg_pduDX_stran_MKE_ysum = Budg_pduDX_stran_MKE_ysum + Budg_pduDX_stran_MKE(J) / DYFI(J)
        Budg_Turbu_diffu_MKE_ysum = Budg_Turbu_diffu_MKE_ysum + Budg_Turbu_diffu_MKE(J) / DYFI(J)
        Budg_DpuDX_diffu_MKE_ysum = Budg_DpuDX_diffu_MKE_ysum + Budg_DpuDX_diffu_MKE(J) / DYFI(J)
        Budg_viscs_diffu_MKE_ysum = Budg_viscs_diffu_MKE_ysum + Budg_viscs_diffu_MKE(J) / DYFI(J)
        Budg_press_accl1_MKE_ysum = Budg_press_accl1_MKE_ysum + Budg_press_accl1_MKE(J) / DYFI(J)
        Budg_viscs_accl1_MKE_ysum = Budg_viscs_accl1_MKE_ysum + Budg_viscs_accl1_MKE(J) / DYFI(J)
        Budg_prodc_gvfc1_MKE_ysum = Budg_prodc_gvfc1_MKE_ysum + Budg_prodc_gvfc1_MKE(J) / DYFI(J)
        Budg_prodc_Dvfc1_MKE_ysum = Budg_prodc_Dvfc1_MKE_ysum + Budg_prodc_Dvfc1_MKE(J) / DYFI(J)
        Budg_Balance1_MKE_ysum   = Budg_Balance1_MKE_ysum   + Budg_Balance1_MKE(J) / DYFI(J)

        DO L = 1, (NDV * (7 - NDV)) / 2 + NDV - 3
            Budg_prodc_stres_Duiuj_ysum(L) = Budg_prodc_stres_Duiuj_ysum(L) + Budg_prodc_stres_Duiuj(J, L) / DYFI(J)
            Budg_viscs_dissp_Duiuj_ysum(L) = Budg_viscs_dissp_Duiuj_ysum(L) + Budg_viscs_dissp_Duiuj(J, L) / DYFI(J)
            Budg_pduDX_stran_Duiuj_ysum(L) = Budg_pduDX_stran_Duiuj_ysum(L) + Budg_pduDX_stran_Duiuj(J, L) / DYFI(J)
            Budg_Turbu_diffu_Duiuj_ysum(L) = Budg_Turbu_diffu_Duiuj_ysum(L) + Budg_Turbu_diffu_Duiuj(J, L) / DYFI(J)
            Budg_DpuDX_diffu_Duiuj_ysum(L) = Budg_DpuDX_diffu_Duiuj_ysum(L) + Budg_DpuDX_diffu_Duiuj(J, L) / DYFI(J)
            Budg_viscs_diffu_Duiuj_ysum(L) = Budg_viscs_diffu_Duiuj_ysum(L) + Budg_viscs_diffu_Duiuj(J, L) / DYFI(J)
            Budg_press_accl1_Duiuj_ysum(L) = Budg_press_accl1_Duiuj_ysum(L) + Budg_press_accl1_Duiuj(J, L) / DYFI(J)
            Budg_viscs_accl1_Duiuj_ysum(L) = Budg_viscs_accl1_Duiuj_ysum(L) + Budg_viscs_accl1_Duiuj(J, L) / DYFI(J)
            Budg_prodc_gvfc2_Duiuj_ysum(L) = Budg_prodc_gvfc2_Duiuj_ysum(L) + Budg_prodc_gvfc2_Duiuj(J, L) / DYFI(J)
            Budg_prodc_Dvfc1_Duiuj_ysum(L) = Budg_prodc_Dvfc1_Duiuj_ysum(L) + Budg_prodc_Dvfc1_Duiuj(J, L) / DYFI(J)
            Budg_Balance1_Duiuj_ysum(L)   = Budg_Balance1_Duiuj_ysum(L)   + Budg_Balance1_Duiuj(J, L) / DYFI(J)
        END DO

    END DO

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.budget.check.' // TRIM(PNTIM) // '.plt'
    INQUIRE(FILE = TRIM(ADJUSTL(FLNM)), exist =File_exists)
    IF(File_exists) THEN
        OPEN(TECFLG, FILE =FLNM, POSITION = 'APPEND')
    ELSE
        OPEN(TECFLG, FILE =FLNM)
    END IF
    WRITE(TECFLG, '(A)') '======= RA=====uu,uv,uw, Vv, Vw,ww========== '
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_stres_Duiuj_ysum(1:6) = ', Budg_prodc_stres_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_dissp_Duiuj_ysum(1:6) = ', Budg_viscs_dissp_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_pduDX_stran_Duiuj_ysum(1:6) = ', Budg_pduDX_stran_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Turbu_diffu_Duiuj_ysum(1:6) = ', Budg_Turbu_diffu_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_DpuDX_diffu_Duiuj_ysum(1:6) = ', Budg_DpuDX_diffu_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_diffu_Duiuj_ysum(1:6) = ', Budg_viscs_diffu_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_press_accl1_Duiuj_ysum(1:6) = ', Budg_press_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_viscs_accl1_Duiuj_ysum(1:6) = ', Budg_viscs_accl1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_gvfc2_Duiuj_ysum(1:6) = ', Budg_prodc_gvfc2_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_prodc_Dvfc1_Duiuj_ysum(1:6) = ', Budg_prodc_Dvfc1_Duiuj_ysum(1:6)
    WRITE(TECFLG, '(A,6ES20.7)') 'Budg_Balance1_Duiuj_ysum(1:6) = ', Budg_Balance1_Duiuj_ysum(1:6)

    WRITE(TECFLG, '(A)') '======= RA===== TKE ========== '
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_stres_TKE_ysum = ', Budg_prodc_stres_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_dissp_TKE_ysum = ', Budg_viscs_dissp_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_pduDX_stran_TKE_ysum = ', Budg_pduDX_stran_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_Turbu_diffu_TKE_ysum = ', Budg_Turbu_diffu_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_DpuDX_diffu_TKE_ysum = ', Budg_DpuDX_diffu_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_diffu_TKE_ysum = ', Budg_viscs_diffu_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_press_accl1_TKE_ysum = ', Budg_press_accl1_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_accl1_TKE_ysum = ', Budg_viscs_accl1_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_gvfc2_TKE_ysum = ', Budg_prodc_gvfc2_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_Dvfc1_TKE_ysum = ', Budg_prodc_Dvfc1_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_Balance1_TKE_ysum = ', Budg_Balance1_TKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Sum of all terms = ', Budg_prodc_stres_TKE_ysum + Budg_viscs_dissp_TKE_ysum + &
    Budg_pduDX_stran_TKE_ysum + Budg_Turbu_diffu_TKE_ysum + Budg_DpuDX_diffu_TKE_ysum + Budg_viscs_diffu_TKE_ysum &
    + Budg_press_accl1_TKE_ysum + Budg_viscs_accl1_TKE_ysum + Budg_prodc_gvfc2_TKE_ysum + Budg_prodc_Dvfc1_TKE_ysum
    WRITE(TECFLG, '(A)') '====== RA====== MKE ========== '
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_stres_MKE_ysum = ', Budg_prodc_stres_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_dissp_MKE_ysum = ', Budg_viscs_dissp_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_pduDX_stran_MKE_ysum = ', Budg_pduDX_stran_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_Turbu_diffu_MKE_ysum = ', Budg_Turbu_diffu_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_DpuDX_diffu_MKE_ysum = ', Budg_DpuDX_diffu_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_diffu_MKE_ysum = ', Budg_viscs_diffu_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_press_accl1_MKE_ysum = ', Budg_press_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_viscs_accl1_MKE_ysum = ', Budg_viscs_accl1_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_gvfc1_MKE_ysum = ', Budg_prodc_gvfc1_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_prodc_Dvfc1_MKE_ysum = ', Budg_prodc_Dvfc1_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Budg_Balance1_MKE_ysum = ', Budg_Balance1_MKE_ysum
    WRITE(TECFLG, '(A, 1ES20.7)') 'Sum of all terms = ', Budg_prodc_stres_MKE_ysum + Budg_viscs_dissp_MKE_ysum + &
    Budg_pduDX_stran_MKE_ysum + Budg_Turbu_diffu_MKE_ysum + Budg_DpuDX_diffu_MKE_ysum + Budg_viscs_diffu_MKE_ysum &
    + Budg_press_accl1_MKE_ysum + Budg_viscs_accl1_MKE_ysum + Budg_prodc_gvfc1_MKE_ysum + Budg_prodc_Dvfc1_MKE_ysum

    CLOSE(TECFLG)

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE PP_HEAT_BASIC_VARS_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, H, P, L, K
    REAL(WP) :: COE
    REAL(WP) :: htauJC, htauJP

    !=============={h}===================================
    ! Eq. {h} =<\rho h>/<\rho>
    H_FA = 0.0_WP
    DO J = 1, NCL2
        H_FA(J) = DHxztL_F0_io(J) / D1xztL_F0_io(J)
        hff_RA(J) = H1xztL_F0_io(J) - H_FA(J)
    END DO
    H_FA(0) = HWAL_FA(iBotWall)
    H_FA(NND2) = HWAL_FA(iTopWall)

    CALL CHKHDL(' ==>Calculated {h} and <h">', MYID)

    !====<p' h'>=<p' h''>=<ph> - <p>*<h>== Hp_per_RA(J, I) =============
    DO J = 1, NCL2
        hfpf_RA(J) =  PHxztL_F0_io(J) - U1xztL_F0_io(J,  4) * H1xztL_F0_io(J)
    END DO
    CALL CHKHDL(' ==>Calculated <h`p`>', MYID)

    !============== D{h}/ DX_M = dHDX_FA(CL, M) ========================
    !============== D<h>/ DX_M = dHDX_RA(CL, M) ========================
    !============== D<T>/ DX_M = dTDX(CL, M) ========================
    M = 2
    dHDX_FA = 0.0_WP
    dHDX_RA = 0.0_WP
    dTDX  = 0.0_WP
    dDDX  = 0.0_WP
    DO J = 1, NCL2

        IF(J == 1) THEN
            dHDX_RA(J, M) = ( ( YCL2ND_WFB(J + 1) * H1xztL_F0_io(J)  + &
            YCL2ND_WFF(J + 1) * H1xztL_F0_io(J + 1) ) - &
            HWAL_RA(iBotWall)  ) * DYFI(J)

            dHDX_FA(J, M) = ( ( YCL2ND_WFB(J + 1) * H_FA(J)  + &
            YCL2ND_WFF(J + 1) * H_FA(J + 1) ) - &
            H_FA(0)  ) * DYFI(J)

            dTDX(J, M) = ( ( YCL2ND_WFB(J + 1) * T1xztL_F0_io(J)  + &
            YCL2ND_WFF(J + 1) * T1xztL_F0_io(J + 1) ) - &
            TWAL(iBotWall)  ) * DYFI(J)
            dDDX(J, M) = ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J)  + &
            YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) - &
            DWAL(iBotWall)  ) * DYFI(J)

        ELSE IF (J == NCL2) THEN
            dHDX_RA(J, M) = (  HWAL_RA(iTopWall) - &
            ( YCL2ND_WFF(J) * H1xztL_F0_io(J) + &
            YCL2ND_WFB(J) * H1xztL_F0_io(J - 1) )  ) * DYFI(J)

            dHDX_FA(J, M) = (  H_FA(NND2) - &
            ( YCL2ND_WFF(J) * H_FA(J) + &
            YCL2ND_WFB(J) * H_FA(J - 1) )  ) * DYFI(J)

            dTDX(J, M) = (  TWAL(iTopWall) - &
            ( YCL2ND_WFF(J) * T1xztL_F0_io(J) + &
            YCL2ND_WFB(J) * T1xztL_F0_io(J - 1) )  ) * DYFI(J)

            dDDX(J, M) = (  DWAL(iTopWall) - &
            ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + &
            YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) )  ) * DYFI(J)

        ELSE
            dHDX_RA(J, M) = ( ( YCL2ND_WFB(J + 1) * H1xztL_F0_io(J  ) + &
            YCL2ND_WFF(J + 1) * H1xztL_F0_io(J + 1) ) - &
            ( YCL2ND_WFF(J) * H1xztL_F0_io(J  ) + &
            YCL2ND_WFB(J) * H1xztL_F0_io(J - 1) ) ) * DYFI(J)

            dHDX_FA(J, M) = ( ( YCL2ND_WFB(J + 1) * H_FA(J  ) + &
            YCL2ND_WFF(J + 1) * H_FA(J + 1) ) - &
            ( YCL2ND_WFF(J) * H_FA(J  ) + &
            YCL2ND_WFB(J) * H_FA(J - 1) ) ) * DYFI(J)

            dTDX(J, M) = ( ( YCL2ND_WFB(J + 1) * T1xztL_F0_io(J  ) + &
            YCL2ND_WFF(J + 1) * T1xztL_F0_io(J + 1) ) - &
            ( YCL2ND_WFF(J) * T1xztL_F0_io(J  ) + &
            YCL2ND_WFB(J) * T1xztL_F0_io(J - 1) ) ) * DYFI(J)

            dDDX(J, M) = ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J  ) + &
            YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) - &
            ( YCL2ND_WFF(J) * D1xztL_F0_io(J  ) + &
            YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) ) * DYFI(J)

        END IF

    END DO
    CALL CHKHDL(' ==>Calculated d<h>/ Dy, d<T>/ Dy and d{h}/ Dy', MYID)


    !==============<\rho>{h'' u''_m} = <\rho h'' u''_m> = uffhffd_FA(CL, M) ========================
    ! Eq: uffhffd_FA(J, M) = <\rho>{h'' u''_m} = <\rho u_m h> - <\rho u_m> * <\rho h>/ <\rho>
    !     UH_FA(J, M) = {u_m h} = <\rho u_m h>/ <\rho>
    DO J = 1, NCL2
        DO M = 1, NDV
            uffhffd_FA(J, M) = GHxztL_F0_io(J, M) - G1xztL_F0_io(J, M) * DHxztL_F0_io(J) / D1xztL_F0_io(J)
            ufhfd_RA(J, M) = (UHxztL_F0_io(J, M) - U1xztL_F0_io(J, M) * H1xztL_F0_io(J)) * D1xztL_F0_io(J)
            UH_FA(J, M)       = GHxztL_F0_io(J, M) / D1xztL_F0_io(J)
        END DO
    END DO
    CALL CHKHDL(' ==>Calculated <uh>, uffhffd_FA, ufhfd_RA', MYID)

    !===<\rho u''_k u''_i h''> = <\rho> {u''_k u''_i h''}==============================
    ! Eq. :
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                L = (M * (7-M)) / 2 + N - 3
                uff2hffd_FA(J, M, N) = U2DHxztL_F0_io(J, L)                        &
                - D1xztL_F0_io(J) * H_FA(J)   * UU_FA(J, M, N) &
                - D1xztL_F0_io(J) * U_FA(J, N) * UH_FA(J, M) &
                - D1xztL_F0_io(J) * U_FA(J, M) * UH_FA(J, N) &
                + 2.0_WP * D1xztL_F0_io(J) * H_FA(J) * U_FA(J, N) * U_FA(J, M)
            END DO
        END DO

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) THEN
                    uff2hffd_FA(J, M, N) = uff2hffd_FA(J, N, M)
                END IF
            END DO
        END DO

    END DO

    !======== viscstressEnth_RA(J, M, N) = <h \tau_mn>======================
    !Eq.
    DO J = 1, NCL2
        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                viscstressEnth_RA(J, M, N) = ( &
                DVDL1MHxztL_F0_io(J, M, N) + &
                DVDL1MHxztL_F0_io(J, N, M) - &
                2.0_WP / 3.0_WP * ( DVDL1MHxztL_F0_io(J, 1, 1) + &
                DVDL1MHxztL_F0_io(J, 2, 2) + &
                DVDL1MHxztL_F0_io(J, 3, 3) ) * DBLE(Kronecker_Delta(M, N)) ) / REN
            END DO
        END DO

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) THEN
                    viscstressEnth_RA(J, M, N) = viscstressEnth_RA(J, N, M)
                END IF
            END DO
        END DO

    END DO

    !========<\partial(h) /\partial x_m \tau_hp>======================
    ! Eq. viscstressEnthGrad_RA(J, M, H, P)
    !   = <mu * d(u_h) / D(x_p) * d(h) / D(x_m) > +
    !       <mu * d(u_p) / D(x_h) * d(h) / D(x_m) > - 2 / 3 *
    !       <mu * d(u_l) / D(x_l) * d(h) / D(x_m) >* Delta_hp
    DO J = 1, NCL2
        DO M = 1, NDV

            DO H = 1, NDV
                DO P = 1, NDV
                    IF(H >  P) CYCLE
                    viscstressEnthGrad_RA(J, M, H, P) = (     &
                    DHDLMDVDLxztL_F0_io(J, M, H, P ) + &
                    DHDLMDVDLxztL_F0_io(J, M, P, H ) - &
                    2.0_WP / 3.0_WP * ( &
                    DHDLMDVDLxztL_F0_io(J, M, 1, 1 ) + &
                    DHDLMDVDLxztL_F0_io(J, M, 2, 2 ) + &
                    DHDLMDVDLxztL_F0_io(J, M, 3, 3 )   &
                    ) * DBLE(Kronecker_Delta(H,P)) ) / REN
                END DO
            END DO

            DO H = 1, NDV
                DO P = 1, NDV
                    IF(H < P) CYCLE
                    viscstressEnthGrad_RA(J, M, H, P) = viscstressEnthGrad_RA(J, M, P, H)
                END DO
            END DO

        END DO
    END DO

    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE PP_HEAT_FA_RSTE_Budg_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    INTEGER(4) :: M, N, H, L, K
    REAL(WP) :: COE
    REAL(WP) :: htauJC, htauJP

    !==PRODUCTION TERMS due to mean sheAR and mean enthalpy gradienT === (THF(turbulent heat flux)) ========
    !Eq: Budg_prod_thf_stres(J, L) = -<\rho h''u''_k> (\partial {u_i} / \partial x_k)
    !Eq: Budg_prod_thf_entpg(J, L) = -<\rho u_i'' u''_k> (\partial {h}   / \partial x_k)
    Budg_prodc_stres_thf = 0.0_WP
    Budg_prodc_enthg_thf = 0.0_WP
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_prodc_stres_thf(J, M) = ( uffhffd_FA(J, 1) * dUDX_FA(J, M, 1) + &
            uffhffd_FA(J, 2) * dUDX_FA(J, M, 2) + &
            uffhffd_FA(J, 3) * dUDX_FA(J, M, 3) ) * (-1.0_WP)
            Budg_prodc_enthg_thf(J, M) = ( uff2d_FA(J, 1, M) * dHDX_FA(J, 1) + &
            uff2d_FA(J, 2, M) * dHDX_FA(J, 2) + &
            uff2d_FA(J, 3, M) * dHDX_FA(J, 3) ) * (-1.0_WP)
        END DO
    END DO

    !== TURBUELCEN DIFFUSION TERMS = (turblence tranport rate) ================================
    !Eq.  Budg_TdIFf_Duiuj(J, L) = ( \partial <\rho u''_i u''_j u''_k > ) / (\partial x_k )
    DO J = 1, NCL2
        DO M = 1, NDV
            IF(J == 1) THEN
                Budg_Turbu_diffu_thf(J, M) = &
                ( ( YCL2ND_WFF(J + 1) * Uff2hffd_FA(J + 1, M, 2) + &
                YCL2ND_WFB(J + 1) * Uff2hffd_FA(J,  M, 2) ) - 0.0_WP ) * DYFI(J) * (-1.0_WP)
            ELSE IF (J == NCL2) THEN
                Budg_Turbu_diffu_thf(J, M) = &
                ( 0.0_WP - &
                ( YCL2ND_WFF(J) * Uff2hffd_FA(J,  M, 2) + &
                YCL2ND_WFB(J) * Uff2hffd_FA(J - 1, M, 2) ) ) * DYFI(J) * (-1.0_WP)

            ELSE
                Budg_Turbu_diffu_thf(J, M) = &
                ( ( YCL2ND_WFF(J + 1) * Uff2hffd_FA(J + 1, M, 2) + &
                YCL2ND_WFB(J + 1) * Uff2hffd_FA(J,  M, 2) ) - &
                ( YCL2ND_WFF(J) * Uff2hffd_FA(J,  M, 2) + &
                YCL2ND_WFB(J) * Uff2hffd_FA(J - 1, M, 2) ) ) * DYFI(J) * (-1.0_WP)
            END IF
        END DO
    END DO

    !==================pressure acceleration terM ===============
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_press_accl1_thf(J, M) = - DPDX_RA(J, M) * Hff_RA(J)
        END DO
    END DO

    !==================enthalpy- PRessure gradient dIFfusioN ===============
    DO J = 1, NCL2
        DO M = 1, NDV
            IF (M == 2) THEN
                IF(J == 1) THEN
                    Budg_DphDX_diffu_thf(J, M) = &
                    ( ( YCL2ND_WFF(J + 1) * Hfpf_RA(J + 1) + &
                    YCL2ND_WFB(J + 1) * Hfpf_RA(J) ) -  &
                    0.0_WP ) * DYFI(J) * (-1.0_WP)
                ELSE IF(J == NCL2) THEN
                    Budg_DphDX_diffu_thf(J, M) = &
                    (  0.0_WP -  &
                    ( YCL2ND_WFF(J) * Hfpf_RA(J) + &
                    YCL2ND_WFB(J) * Hfpf_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
                ELSE
                    Budg_DphDX_diffu_thf(J, M) = &
                    ( ( YCL2ND_WFF(J + 1) * Hfpf_RA(J + 1) + &
                    YCL2ND_WFB(J + 1) * Hfpf_RA(J) ) -  &
                    ( YCL2ND_WFF(J) * Hfpf_RA(J) + &
                    YCL2ND_WFB(J) * Hfpf_RA(J - 1) ) ) * DYFI(J) * (-1.0_WP)
                END IF
            ELSE
                Budg_DphDX_diffu_thf(J, M) = 0.0_WP
            END IF
        END DO
    END DO

    !==================enthalpy- PRessure strAIN ===============
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_pdhDX_stran_thf(J, M) = DHDLPxztL_F0_io(J, M) - U1xztL_F0_io(J, 4) * DHDL1xztL_F0_io(J, M)
        END DO
    END DO

    !================conductive heat flux acceleratioN =============
    DO J = 1, NCL2
        DO M = 1, NDV
            IF(J == 1) THEN
                Budg_ConHF_accel_thf(J, M) = uff_RA(J, M) * &
                ( ( YCL2ND_WFF(J + 1) * DTDLKxztL_F0_io(J + 1, 2) + &
                YCL2ND_WFB(J + 1) * DTDLKxztL_F0_io(J,  2) ) -  &
                Qw(iBotWall) ) * DYFI(J) * CTHECD

            ELSE IF(J == NCL2) THEN
                Budg_ConHF_accel_thf(J, M) = uff_RA(J, M) * &
                (  Qw(iTopWall) -  &
                ( YCL2ND_WFF(J) * DTDLKxztL_F0_io(J,  2) + &
                YCL2ND_WFB(J) * DTDLKxztL_F0_io(J - 1, 2) ) ) * DYFI(J) * CTHECD
            ELSE
                Budg_ConHF_accel_thf(J, M) = uff_RA(J, M) * &
                ( ( YCL2ND_WFF(J + 1) * DTDLKxztL_F0_io(J + 1, 2) + &
                YCL2ND_WFB(J + 1) * DTDLKxztL_F0_io(J,  2) ) -  &
                ( YCL2ND_WFF(J) * DTDLKxztL_F0_io(J,  2) + &
                YCL2ND_WFB(J) * DTDLKxztL_F0_io(J - 1, 2) ) ) * DYFI(J) * CTHECD
            END IF

        END DO
    END DO

    !================conductive heat flux dIFfusioN =============
    DO J = 1, NCL2
        DO M = 1, NDV
            IF(J == 1) THEN
                Budg_ConHF_diffu_thf(J, M) = &
                ( ( YCL2ND_WFF(J + 1) * ( DTDLKUxztL_F0_io(J + 1, 2, M) - U1xztL_F0_io(J + 1, M) * DTDLKxztL_F0_io(J + 1, 2) ) + &
                YCL2ND_WFB(J + 1) * ( DTDLKUxztL_F0_io(J, 2, M) - U1xztL_F0_io(J, M) * DTDLKxztL_F0_io(J,  2) ) ) -&
                0.0_WP )&
                * DYFI(J) * CTHECD
            ELSE IF(J == NCL2) THEN
                Budg_ConHF_diffu_thf(J, M) = &
                ( 0.0_WP -&
                ( YCL2ND_WFF(J) * ( DTDLKUxztL_F0_io(J, 2, M) - U1xztL_F0_io(J, M) * DTDLKxztL_F0_io(J,  2) ) + &
                YCL2ND_WFB(J) * ( DTDLKUxztL_F0_io(J - 1, 2, M) - U1xztL_F0_io(J - 1, M) * DTDLKxztL_F0_io(J - 1, 2) ) ) )&
                * DYFI(J) * CTHECD
            ELSE
                Budg_ConHF_diffu_thf(J, M) = &
                ( ( YCL2ND_WFF(J + 1) * ( DTDLKUxztL_F0_io(J + 1, 2, M) - U1xztL_F0_io(J + 1, M) * DTDLKxztL_F0_io(J + 1, 2) ) + &
                YCL2ND_WFB(J + 1) * ( DTDLKUxztL_F0_io(J, 2, M) - U1xztL_F0_io(J, M) * DTDLKxztL_F0_io(J,  2) ) ) -&
                ( YCL2ND_WFF(J) * ( DTDLKUxztL_F0_io(J, 2, M) - U1xztL_F0_io(J, M) * DTDLKxztL_F0_io(J,  2) ) + &
                YCL2ND_WFB(J) * ( DTDLKUxztL_F0_io(J - 1, 2, M) - U1xztL_F0_io(J - 1, M) * DTDLKxztL_F0_io(J - 1, 2) ) ) )&
                * DYFI(J) * CTHECD
            END IF
        END DO
    END DO

    !============conductive heat flux dISsipatioN ==================
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_ConHF_dissp_thf(J, M) = -1.0_WP * ( dTdLKdVdLxztL_F0_io(J, 1, M, 1) + &
            dTdLKdVdLxztL_F0_io(J, 2, M, 2) + &
            dTdLKdVdLxztL_F0_io(J, 3, M, 3) ) * CTHECD + &
            ( DTDLKxztL_F0_io(J, 1) * DVDL1xztL_F0_io(J, M, 1) + &
            DTDLKxztL_F0_io(J, 2) * DVDL1xztL_F0_io(J, M, 2) + &
            DTDLKxztL_F0_io(J, 3) * DVDL1xztL_F0_io(J, M, 3) ) * CTHECD
        END DO
    END DO

    !======= VIScous acceleration, including alL ======== (based on tau_RA) ===================================
    !       Eq. <h''> (\partial <tau_ki>) / (\partial  x_k)
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_viscs_accl1_thf(J, M) = dTaudy_RA(J, 2, M) * Hff_RA(J)
        END DO
    END DO

    !======= VIScous dIFfusion term ===========================================================
    !       Eq. = \partial <h'' tau'_ki> / pARtial (x_k)
    DO J = 1, NCL2
        DO M = 1, NDV
            IF(J == 1) THEN
                htauJP = YCL2ND_WFF(J + 1) * ( viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)) +&
                YCL2ND_WFB(J + 1) * ( viscstressEnth_RA(J + 1, 2, M) - H1xztL_F0_io(J + 1) * Tau_Mean_RA(J + 1, 2, M))
                htauJC = viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)
                Budg_viscs_diffu_thf(J, M) = (htauJP - HtauJC) * DYFI(J) * 2.0_WP
            ELSE IF(J == NCL2) THEN
                htauJP = viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)
                htauJC = YCL2ND_WFF(J) * ( viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)) +&
                YCL2ND_WFB(J) * ( viscstressEnth_RA(J - 1, 2, M) - H1xztL_F0_io(J - 1) * Tau_Mean_RA(J - 1, 2, M))
                Budg_viscs_diffu_thf(J, M) = (htauJP - HtauJC) * DYFI(J) * 2.0_WP
            ELSE
                htauJP = YCL2ND_WFF(J + 1) * ( viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)) +&
                YCL2ND_WFB(J + 1) * ( viscstressEnth_RA(J + 1, 2, M) - H1xztL_F0_io(J + 1) * Tau_Mean_RA(J + 1, 2, M))
                htauJC = YCL2ND_WFF(J) * ( viscstressEnth_RA(J,  2, M) - H1xztL_F0_io(J  ) * Tau_Mean_RA(J,  2, M)) +&
                YCL2ND_WFB(J) * ( viscstressEnth_RA(J - 1, 2, M) - H1xztL_F0_io(J - 1) * Tau_Mean_RA(J - 1, 2, M))
                Budg_viscs_diffu_thf(J, M) = (htauJP - HtauJC) * DYFI(J)
            END IF

        END DO
    END DO

    !======= VIScous ENERGY dISsipation term ===========================================================
    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_viscs_dissp_thf(J, M) = -1.0_WP * ( &
            viscstressEnthGrad_RA(J, 1, 1, M) + &
            viscstressEnthGrad_RA(J, 2, 2, M) + &
            viscstressEnthGrad_RA(J, 3, 3, M) - &
            Tau_Mean_RA(J, 1, M) * DHDX_RA(J, 1) - &
            Tau_Mean_RA(J, 2, M) * DHDX_RA(J, 2) - &
            Tau_Mean_RA(J, 3, M) * DHDX_RA(J, 3) )
        END DO
    END DO

    DO J = 1, NCL2
        DO M = 1, NDV
            Budg_Balance1_thf(J, M) = Budg_prodc_stres_thf(J, M) +Budg_prodc_enthg_thf(J, M) + &
            Budg_Turbu_diffu_thf(J, M) +Budg_press_accl1_thf(J, M) + &
            Budg_DphDX_diffu_thf(J, M) +Budg_pdhDX_stran_thf(J, M) + &
            Budg_ConHF_accel_thf(J, M) +Budg_ConHF_diffu_thf(J, M) + &
            Budg_ConHF_dissp_thf(J, M) +Budg_viscs_accl1_thf(J, M) + &
            Budg_viscs_diffu_thf(J, M) +Budg_viscs_dissp_thf(J, M)
        END DO
    END DO


    !================BODY FORCE/ BUOYANCY PRODUCTION, whICh IS pARt of pressure acceleration terM ===========
    !       ! not an independe contribution, but belongs to pARt of Budg_press_accl1_Duiuj
    Budg_prodc_gvfc2_DuiuJ = 0.0_WP
    DO J = 1, NCL2

        DO M = 1, NDV
            DO N = 1, NDV
                IF(M >  N) CYCLE
                IF(M /= ABS(iGravity) .AND. N /= ABS(iGravity)) CYCLE

                L = (M * (7-M)) / 2 + N - 3

                IF(M == ABS(iGravity) .AND. N == ABS(iGravity)) THEN
                    COE = 2.0_WP
                    K = ABS(iGravity)
                ELSE IF (M == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = N
                ELSE IF (N == ABS(iGravity)) THEN
                    COE = 1.0_WP
                    K = M
                ELSE
                    COE = 0.0_WP
                    K = -1 !(WHICH will LEAD TO ERROR!)
                END IF
                ! F_A INCLUDEs a postive or negtive sign
                Budg_prodc_gvfc2_Duiuj(J, L) = F_A * COE * ( G1xztL_F0_io(J, K) - D1xztL_F0_io(J) * U1xztL_F0_io(J, K) )

                !du= ( G1xztL_F0_io(J, 1) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 1) ) * F_A
                !dV = ( G1xztL_F0_io(J, 2) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 2) ) * F_A
                !dw= ( G1xztL_F0_io(J, 3) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 3) ) * F_A

            END DO
        END DO

        IF(ABS(iGravity) == 0) THEN
            Budg_prodc_gvfc2_thf(J) = 0.0_WP
        ELSE
            Budg_prodc_gvfc2_thf(J) = F_A * ( DHxztL_F0_io(J) - D1xztL_F0_io(J) * H1xztL_F0_io(J) )
        END IF

    END DO

    RETURN
END SUBROUTINE



!**********************************************************************************************************************************
SUBROUTINE WRT_FLOW_FA_Profile_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    CHARACTER(128) :: FLNM, FLNN
    CHARACTER(5) :: STDIM(2)
    CHARACTER(4) :: STFASTRESS(8)
    CHARACTER(2) :: SJ2
    REAL(WP) :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
    REAL(WP) :: tke,Ruv_vIS
    INTEGER(4) :: I, J, N, L, TECFLG_FavAG(2), TECFLG_FavAG3, NMax
    REAL(WP) ::   COE, tem, M_tmp
    REAL(WP) ::   TKE2, EPPSI, Fmu
    REAL(WP) ::   Cmu, K2De

    REAL(WP) :: FCT(0 : NND2, NDV)
    REAL(WP) :: COE1, COE2, DENtemp, Dintg, intgbdfc

    IF(iThermoDynamics /= 1) RETURN

    !====================Favre Avderaged profileS ===============================================================
    NMax = 1
    STDIM(1) = 'undim'
    TECFLG_FavAG(1) = 51
    IF(iPPDimension == 1) THEN
        NMax = 2
        STDIM(2) = 'dimen'
        TECFLG_FavAG(2) = 52
    END IF


    !scaling1 = U0 * U0 * D0
    DO N = 1, NMax
        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(N)) // '.Profile.Flow.Favre.' // TRIM(PNTIM) // '.plt'
        OPEN (TECFLG_FavAG(N), FILE = TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_FavAG(N), '(A)') 'TITLE = " Favre Averged Flow (35 variables)" '
        J = 0;                          WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") 'variables = '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ux", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Uy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Uz", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'P", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rvv", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rww", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruv", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rvw", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RuvvIS", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DUDY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dilation", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DensIntg", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'bdfcintg", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsX", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsZ", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoUU", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoVV", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoWW", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MKE", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'TKE", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUpp1_FA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUpp2_FA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUpp3_FA"'
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUpp112FA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUpp332FA", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_U_FA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_V_FA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_W_FA"'

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Anistpinva2", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Anistpinva3", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'LumleyX", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'LumleyY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)'          ) '"' // SJ2 // 'drvFC" '

        WRITE(TECFLG_FavAG(N), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    END DO

    DO J = 1, NCL2

        WRITE(TECFLG_FavAG(1), '(46ES20.12)') &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J), &
        U_FA(J, 1), U_FA(J, 2), U_FA(J, 3), U1xztL_F0_io(J, 4), &
        uff2_FA(J, 1, 1), &
        uff2_FA(J, 2, 2), &
        uff2_FA(J, 3, 3), &
        uff2_FA(J, 1, 2), &
        uff2_FA(J, 1, 3), &
        uff2_FA(J, 2, 3), &
        Tau_Mean_RA(J, 1, 2), &
        dUDX_FA(J, 1, 2), &
        dUiDXi(J), &
        DensIntg(J), &
        bdfcintg(J), &
        Omega_rms(J, 1:3), &
        UGxztL_F0_io(J, 1), &
        UGxztL_F0_io(J, 4), &
        UGxztL_F0_io(J,6), &
        MKE_FA(J), &
        TKE_FA(J), &
        uff3_FA(J, 1, 1, 1), &
        uff3_FA(J, 2, 2, 2), &
        uff3_FA(J, 3, 3, 3), &
        uff3_FA(J, 1, 1, 2), &
        uff3_FA(J, 3, 3, 2), &
        skewness_FA(J, 1:3), &
        Anistpinva_FA(J, 2:3), &
        LumleyAxis_FA(J, 1:2), &
        FUxztL_F0_io(J, 4)

    END DO
    DO N = 1, NMax
        CLOSE(TECFLG_FavAG(N))
    END DO

    CALL WRT_FLOW_Budgets_Profile_XZ_io('Favre')

    !======================================================================================================
    N = 1
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.VALIDATION.Profile.Flow.FANS.' // TRIM(PNTIM) // '.plt'

    OPEN (TECFLG_FavAG(N), FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_FavAG(N), '(A)') 'TITLE = " Favre Averged Flow (20 variables)" '
    J = 0;                           WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") 'variables = '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'X - VTau12 - 1", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'X - VTau12 -2", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'X - TTauT12", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'X -bdfc", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'X - Total", '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y- PRessure", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y- VTau22 - 1", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y- VTau22 -2", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y- TTauT22", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y-bdfc", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Y- Total", '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Z - VTau23- 1", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Z - VTau23-2", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'z-tTauT23", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Z -bdfc", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)')              '"' // SJ2 // 'z-total"'

    WRITE(TECFLG_FavAG(N), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    DO J = 1, NCL2
        WRITE(TECFLG_FavAG(N), '(26ES20.12)')     &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J), &
        Tau_meaU_RA(J, 1, 2), Tau_Mean_RA(J, 1, 2) - Tau_meaU_RA(J, 1, 2), - Uff2d_FA(J, 1, 2), F_A* IBuoF(1) * D1xztL_F0_io(J), &
        Tau_Mean_RA(J, 1, 2) - Uff2d_FA(J, 1, 2) +F_A* IBuoF(1) * D1xztL_F0_io(J), &
        - U1xztL_F0_io(J, 4), Tau_meaU_RA(J, 2, 2), Tau_Mean_RA(J, 2, 2) - Tau_meaU_RA(J, 2, 2), - Uff2d_FA(J, 2, 2), &
        F_A* IBuoF(2) * D1xztL_F0_io(J), &
        - U1xztL_F0_io(J, 4) + Tau_Mean_RA(J, 2, 2) - Uff2d_FA(J, 2, 2) +F_A* IBuoF(2) * D1xztL_F0_io(J), &
        Tau_meaU_RA(J, 2, 3), &
        Tau_Mean_RA(J, 2, 3) - Tau_meaU_RA(J, 2, 3), &
        - Uff2d_FA(J, 2, 3), &
        F_A* IBuoF(3) * D1xztL_F0_io(J), &
        Tau_Mean_RA(J, 2, 3) - Uff2d_FA(J, 2, 3) +F_A* IBuoF(3) * D1xztL_F0_io(J)
    END DO

    !========================================================================================
    ! show RA FA relations
    !        N = 1
    !        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    !        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(N)) // '.Profile.Flow.FA_RA_Comp.' // TRIM(PNTIM) // '.plt'
    !        OPEN (TECFLG_FavAG(N), FILE = TRIM(ADJUSTL(FLNM)))
    !        WRITE(TECFLG_FavAG(N), '(A)') 'TITLE = " Favre Averged Flow (29 variables)" '
    !        J = 0;                          WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") 'variables = '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Utw", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Utave", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dave", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'upp", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'vpp", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'wpp", '
    !        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)'          ) '"' // SJ2 // 'Mut"'
    !        WRITE(TECFLG_FavAG(N), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    !        DO J = 1, NCL2

    !            Ruv_vIS = Tau_Mean_RA(J, 1, 2)!
    !            tke = 0.5_WP * (uff2d_FA(J, 1, 1) + Uff2d_FA(J, 2, 2) + Uff2d_FA(J, 3, 3))
    !            EPPSI = 0.5_WP * (Budg_viscs_dissp_Duiuj(J, 1) +Budg_viscs_dissp_Duiuj(J, 4) +Budg_viscs_dissp_Duiuj(J,6)) / D1xztL_F0_io(J)
    !            !http://www.cfD -online.coM /WIkI / Turbulence_DISsipation_rate

    !            WRITE(TECFLG_FavAG(N), '(29ES20.12)') &
    !                UtauSD(J), DensSD(J), viscsD(J), Utaw_ave_io, DenAvew, VisAvew, YCC(J), YWdISD(J), &
    !                D1xztL_F0_io(J), M1xztL_F0_io(J), &
    !                U1xztL_F0_io(J, 1) - U_FA(J, 1), U1xztL_F0_io(J, 2) - U_FA(J, 2), U1xztL_F0_io(J, 3) - U_FA(J, 3)
    !        END DO
    !        CLOSE(TECFLG_FavAG(N))



    RETURN
END SUBROUTINE


!**********************************************************************************************************************************
SUBROUTINE WRT_FLOW_RA_Profile_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    CHARACTER(128) :: FLNM, FLNN
    CHARACTER(5) :: STDIM(2)
    CHARACTER(4) :: STFASTRESS(8)
    CHARACTER(2) :: SJ2
    CHARACTER(3) :: SJ3
    REAL(WP) :: ux, uy, uz, p, &
    tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, vorx, vory, vorz
    REAL(WP) ::   COE, tem, M_tmp,Ruv_vIS, DUDY
    REAL(WP) :: scaling1, scaling2
    INTEGER(4) :: I, J, L
    INTEGER(4) :: TECFLG_ReyAG(2), TECFLG_ReyAG3
    INTEGER(4) :: N, NMax, M
    REAL(WP) :: TKE2, EPPSI, Fmu
    REAL(WP) :: Cmu, K2De
    REAL(WP) :: FCT(0 : NCL2 + 1, NDV)
    REAL(WP) :: COE1, COE2
    REAL(WP) :: Yplus

    !========================== TitlE =========================================
    NMax = 1
    STDIM(1) = 'undim'
    TECFLG_ReyAG(1) = 51
    IF(iPPDimension == 1) THEN
        NMax = 2
        STDIM(2) = 'dimen'
        TECFLG_ReyAG(2) = 102
    END IF

    !scaling1 = U0 * U0 * D0

    DO N = 1, NMax
        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(N)) // '.Profile.Flow.Reynd.' // TRIM(PNTIM) // '.plt'
        OPEN (TECFLG_ReyAG(N), FILE = TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_ReyAG(N), '(A)')     'TITLE = " Reynods Averaged Flow (41 variables)" '
        J = 0;                          WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") 'variables = '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ux", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Uy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Uz", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'P", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rvv", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rww", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruv", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ruw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Rvw", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RuvvIS", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DUDY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dialation", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DensIntg", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'bdfcintg", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsX", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'OmegRmsZ", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoUU", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoVV", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'RhoWW", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MKE", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'TKE", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUp1_RA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUp2_RA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUp3_RA"'
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUp112RA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'UUUp332RA"'

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_U_RA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_V_RA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'skewness_W_RA"'


        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Anistpinva2", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Anistpinva3", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'LumleyX", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'LumleyY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(N), '(A)'          ) '"' // SJ2 // 'drvFC" '



        WRITE(TECFLG_ReyAG(N), '(A)')     'ZONE T = " ' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    END DO

    !==========================WRITE Flow info=========================================
    DO J = 1, NCL2

        WRITE(TECFLG_ReyAG(1), '(46ES20.12)') &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J), &
        U1xztL_F0_io(J, 1), &
        U1xztL_F0_io(J, 2), &
        U1xztL_F0_io(J, 3), &
        U1xztL_F0_io(J, 4), &
        uf2_RA(J, 1, 1), &
        uf2_RA(J, 2, 2), &
        uf2_RA(J, 3, 3), &
        uf2_RA(J, 1, 2), &
        uf2_RA(J, 1, 3), &
        uf2_RA(J, 2, 3), &
        Tau_Mean_RA(J, 1, 2), &
        DVDL1xztL_F0_io(J, 1, 2), &
        dUiDXi(J), &
        DensIntg(J), &
        bdfcintg(J), &
        Omega_rms(J, 1:3), &
        UU_RA(J, 1, 1) * D1xztL_F0_io(J), &
        UU_RA(J, 2, 2) * D1xztL_F0_io(J), &
        UU_RA(J, 3, 3) * D1xztL_F0_io(J), &
        MKE_RA(J), &
        TKE_RA(J), &
        uf3_RA(J, 1, 1, 1), &
        uf3_RA(J, 2, 2, 2), &
        uf3_RA(J, 3, 3, 3), &
        uf3_RA(J, 1, 1, 2), &
        uf3_RA(J, 3, 3, 2), &
        skewness_RA(J, 1:3), &
        Anistpinva_RA(J, 2:3), &
        LumleyAxis_RA(J, 1:2), &
        FUxztL_F0_io(J, 4)

    END DO
    DO N = 1, NMax
        CLOSE(TECFLG_ReyAG(N))
    END DO

    CALL WRT_FLOW_Budgets_Profile_XZ_io('Reynd')

    !======================= Quadrant termS =======================
if(iPPQuadrants == 1) then
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.Quad.Reynd.' // TRIM(PNTIM) // '.plt'

    OPEN (TECFLG_ReyAG(1), FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_ReyAG(1), '(A)') 'TITLE = " Reynolds Averged Flow (10 + 16*9 variables)" '
    J = 0;                          WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") 'variables = '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '
    IF(iPPQuadrants == 1)  THEN
        DO N = 1, QUADHN
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADHV", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDR1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDR2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDR3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDR4", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADUV1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADUV2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADUV3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADUV4", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV11", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV12", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV13", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV14", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV21", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV22", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV23", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADDUV24", '


            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADVz1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADVz2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADVz3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADVz4", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADTK1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADTK2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADTK3", '
            IF (N == QUADHN) THEN
                J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)'          ) '"' // SJ3 // 'QUADTK4", '
            ELSE
                J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'QUADTK4", '
            END IF
        END DO
    END IF

    WRITE(TECFLG_ReyAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    DO J = 1, NCL2

        WRITE(TECFLG_ReyAG(1), '(10ES20.12)', AdvancE = "no") &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J)
        IF(iPPQuadrants == 1)  THEN
            DO N = 1, QUADHn - 1
                WRITE(TECFLG_ReyAG(1), '(25ES20.12)', AdvancE = "no") &
                QUADHV(N), &
                QUADDRxztL_F0_io(J, 1:4, N), &
                QUADUVxztL_F0_io(J, 1:4, N), &
                QUADDUV1xztL_F0_io(J, 1:4, N), &
                QUADDUV2xztL_F0_io(J, 1:4, N), &
                QUADVzxztL_F0_io(J, 1:4, N), &
                QUADTKxztL_F0_io(J, 1:4, N)
            END DO
            N = QUADHN
            WRITE(TECFLG_ReyAG(1), '(25ES20.12)') &
            QUADHV(N), &
            QUADDRxztL_F0_io(J, 1:4, N), &
            QUADUVxztL_F0_io(J, 1:4, N), &
            QUADDUV1xztL_F0_io(J, 1:4, N), &
            QUADDUV2xztL_F0_io(J, 1:4, N), &
            QUADVzxztL_F0_io(J, 1:4, N), &
            QUADTKxztL_F0_io(J, 1:4, N)
        END IF

    END DO
    CLOSE(TECFLG_ReyAG(1))

    !=======================octant termS ====\rho=================

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.OctD.Reynd.' // TRIM(PNTIM) // '.plt'

    OPEN (TECFLG_ReyAG(1), FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_ReyAG(1), '(A)') 'TITLE = " Reynolds Averged Flow (10 + 16*9 variables)" '
    J = 0;                          WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") 'variables = '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

    DO N = 1, QUADHN
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDHV", '

        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR1", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR2", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR3", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR4", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR5", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR6", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR7", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDR8", '

        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV1", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV2", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV3", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV4", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV5", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV6", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV7", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDUV8", '

        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV11", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV12", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV13", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV14", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV15", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV16", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV17", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV18", '

        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV21", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV22", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV23", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV24", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV25", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV26", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV27", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDDUV28", '


        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz1", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz2", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz3", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz4", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz5", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz6", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz7", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDVz8", '

        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK1", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK2", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK3", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK4", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK5", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK6", '
        J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK7", '
        IF (N == QUADHN) THEN
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)'          ) '"' // SJ3 // 'OctDTK8", '
        ELSE
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctDTK8", '
        END IF
    END DO

    WRITE(TECFLG_ReyAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    DO J = 1, NCL2

        WRITE(TECFLG_ReyAG(1), '(10ES20.12)', AdvancE = "no") &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J)

        DO N = 1, QUADHn - 1
            WRITE(TECFLG_ReyAG(1), '(49ES20.12)', AdvancE = "no") &
            QuadHV(N), &
            OctDDRxztL_F0_io(J, 1:8, N), &
            OctDUVxztL_F0_io(J, 1:8, N), &
            OctDDUV1xztL_F0_io(J, 1:8, N), &
            OctDDUV2xztL_F0_io(J, 1:8, N), &
            OctDVzxztL_F0_io(J, 1:8, N), &
            OctDTKxztL_F0_io(J, 1:8, N)
        END DO
        N = QUADHN
        WRITE(TECFLG_ReyAG(1), '(49ES20.12)') &
        QuadHV(N), &
        OctDDRxztL_F0_io(J, 1:8, N), &
        OctDUVxztL_F0_io(J, 1:8, N), &
        OctDDUV1xztL_F0_io(J, 1:8, N), &
        OctDDUV2xztL_F0_io(J, 1:8, N), &
        OctDVzxztL_F0_io(J, 1:8, N), &
        OctDTKxztL_F0_io(J, 1:8, N)


    END DO
    CLOSE(TECFLG_ReyAG(1))

    !=======================octant termS ==== T =================

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.OctT.Reynd.' // TRIM(PNTIM) // '.plt'

    OPEN (TECFLG_ReyAG(1), FILE = TRIM(ADJUSTL(FLNM)))
    WRITE(TECFLG_ReyAG(1), '(A)') 'TITLE = " Reynolds Averged Flow (10 + 16*9 variables)" '
    J = 0;                          WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") 'variables = '

    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
    J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '
    IF(iPPQuadrants == 1)  THEN
        DO N = 1, QUADHN
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTHV", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR4", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR5", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR6", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR7", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDR8", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV4", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV5", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV6", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV7", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTUV8", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV11", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV12", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV13", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV14", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV15", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV16", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV17", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV18", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV21", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV22", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV23", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV24", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV25", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV26", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV27", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTDUV28", '


            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz4", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz5", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz6", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz7", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTVz8", '

            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK1", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK2", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK3", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK4", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK5", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK6", '
            J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK7", '
            IF (N == QUADHN) THEN
                J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)'          ) '"' // SJ3 // 'OctTTK8", '
            ELSE
                J = J + 1; WRITE(SJ3, '(1I3.3)') J; WRITE(TECFLG_ReyAG(1), '(A)', AdvancE = "no") '"' // SJ3 // 'OctTTK8", '
            END IF
        END DO
    END IF

    WRITE(TECFLG_ReyAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

    DO J = 1, NCL2

        WRITE(TECFLG_ReyAG(1), '(10ES20.12)', AdvancE = "no") &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J)
        IF(iPPQuadrants == 1)  THEN
            DO N = 1, QUADHn - 1
                WRITE(TECFLG_ReyAG(1), '(49ES20.12)', AdvancE = "no") &
                QuadHV(N), &
                OctTDRxztL_F0_io(J, 1:8, N), &
                OctTUVxztL_F0_io(J, 1:8, N), &
                OctTDUV1xztL_F0_io(J, 1:8, N), &
                OctTDUV2xztL_F0_io(J, 1:8, N), &
                OctTVzxztL_F0_io(J, 1:8, N), &
                OctTTKxztL_F0_io(J, 1:8, N)
            END DO
            N = QUADHN
            WRITE(TECFLG_ReyAG(1), '(49ES20.12)') &
            QuadHV(N), &
            OctTDRxztL_F0_io(J, 1:8, N), &
            OctTUVxztL_F0_io(J, 1:8, N), &
            OctTDUV1xztL_F0_io(J, 1:8, N), &
            OctTDUV2xztL_F0_io(J, 1:8, N), &
            OctTVzxztL_F0_io(J, 1:8, N), &
            OctTTKxztL_F0_io(J, 1:8, N)
        END IF


    END DO
    CLOSE(TECFLG_ReyAG(1))

end if

    RETURN
END SUBROUTINE



!**********************************************************************************************************************************
SUBROUTINE WRT_HEAT_FA_Profile_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    CHARACTER(128) :: FLNM, FLNN
    CHARACTER(5) :: STDIM(2)
    CHARACTER(4) :: STFASTRESS(3)
    CHARACTER(2) :: SJ2
    REAL(WP) :: ux, uy, uz, p, &
    tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
    den, tem, enh_ra, enh_fa, drms, trms, &
    hrms_ra, hrms_fa, Thfx_ra, thfy_ra, thfz_ra, &
    thfx_fa, thfy_fa, thfz_fa
    REAL(WP) :: CpT(NCL2), Cp_eval
    INTEGER(4) :: I, J, TECFLG_FavAG(2), L, N, NMax
    REAL(WP) :: COE2, COE3
    REAL(WP) :: H_tmp, Ptemp, D_tmp, T_tmp, M_tmp, K_tmp, Cp_tmp,  B_tmp, DH_tmp, PRTMP, Qflux, Du, Dv, Dw, Dhe
    REAL(WP) :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
    REAL(WP) :: du_per, dv_per, dw_per, dhe_per
    REAL(WP) :: MutOverPrt(3), Pruv(3)
    REAL(WP) :: FCT(0 : NCL2 + 1)
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCp

    NMax = 1
    STDIM(1) = 'undim'
    TECFLG_FavAG(1) = 51
    IF(iPPDimension == 1) THEN
        NMax = 2
        STDIM(2) = 'dimen'
        TECFLG_FavAG(2) = 102
    END IF
    !==================profileS ================================================
    DO N = 1, NMax
        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(N)) // '.Profile.Heat.Transfer.' // TRIM(PNTIM) // '.plt'

        OPEN (TECFLG_FavAG(N), FILE = TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_FavAG(N), '(A)') 'TITLE = " Thermal variables (29 variables)" '
        J = 0;                          WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") 'variables = '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Hw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Qw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Cpw", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'T", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'HRA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'HFA", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D(T)", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M(T)", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'K(T)", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Cp(T)", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Pr(T)", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Drms", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Trms", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'HRARms", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'HFARms", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfRAx", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfRAy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfRAz", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfFAx", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfFAy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'ThfFAz", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'QfluxX", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'QfluxY", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'QfluxZ", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'du_per", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dv_per", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dw_per", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dh_per", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dTdy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dHRAdy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dHFAdy", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(N), '(A)'          ) '"' // SJ2 // 'dDdy", '
        WRITE(TECFLG_FavAG(N), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '
    END DO

    DO J = 1, NCL2

        drms  = DSQRT(DABS( D2xztL_F0_io(J) - D1xztL_F0_io(J) * D1xztL_F0_io(J) ) ) ! same
        trms  = DSQRT(DABS( T2xztL_F0_io(J) - T1xztL_F0_io(J) * T1xztL_F0_io(J) ) ) ! same
        hrms_ra = DSQRT(DABS( H2xztL_F0_io(J) - H1xztL_F0_io(J) * H1xztL_F0_io(J) ) )
        hrms_fa = DSQRT(DABS( H2xztL_F0_io(J) - DHxztL_F0_io(J) / D1xztL_F0_io(J) * &
        (2.0_WP * H1xztL_F0_io(J) - DHxztL_F0_io(J) / D1xztL_F0_io(J) ) ) )

        IF(ABS(iGravity) == 0) THEN
            du_peR = 0.0_WP
            dv_peR = 0.0_WP
            dw_peR = 0.0_WP
            dhe_peR = 0.0_WP
        ELSE
            du_peR = ( G1xztL_F0_io(J, 1) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 1) ) * F_A
            dv_peR = ( G1xztL_F0_io(J, 2) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 2) ) * F_A
            dw_peR = ( G1xztL_F0_io(J, 3) - D1xztL_F0_io(J) * U1xztL_F0_io(J, 3) ) * F_A
            dhe_peR = ( DHxztL_F0_io(J)   - D1xztL_F0_io(J) * H1xztL_F0_io(J)   ) * F_A
        END IF

        T_tmp = T1xztL_F0_io(J)
        CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T_tmp, D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)

        PRTMP = Cp_tmp * M_tmp / K_tmp
        CpT(J) =    Cp_tmp

        WRITE(TECFLG_FavAG(1), '(43ES20.12)') &
        TauwSD(J), &
        DensSD(J), DenAvew, D1xztL_F0_io(J), &
        viscsD(J), VisAvew, M1xztL_F0_io(J), &
        F_A, YCC(J), YWdISD(J), &
        TwSD(J), HwSD(J), &
        QwSD(J), CpSD(J), &
        T1xztL_F0_io(J), &
        H1xztL_F0_io(J), &
        H_FA(J), &
        D_tmp, M_tmp, K_tmp, Cp_tmp, PRTMP, &
        drms, trms, hrms_ra, hrms_fa, &
        ufhfd_RA(J, 1:3), &
        uffhffd_FA(J, 1:3), &
        - DTDLKxztL_F0_io(J, 1:3) * CTHECD, &
        du_per, Dv_per, Dw_per, Dhe_per, &
        dTDX(J, 2), dHDX_RA(J, 2), dHDX_FA(J, 2), dDDX(J, 2)
        !            IF(iPPDimension == 1) THEN
        !                WRITE(TECFLG_FavAG(2), '(29ES20.12)') YCC(J) * L0, (1.0_WP - DABS(YCC(J))) * REN * COE, &
        !                        D1xztL_F0_io(J) * D0, T1xztL_F0_io(J) * T0, &
        !                        H1xztL_F0_io(J) * CP0 * T0 + H0, &
        !                        H_FA(J) * CP0 * T0 + H0, &
        !                        M1xztL_F0_io(J) * M0, &
        !                        drms* D0, trms* T0, &
        !                        hrms_ra* CP0 * T0 + H0, hrms_fa* CP0 * T0 + H0, &
        !                        D_tmp * D0, M_tmp * M0, K_tmp * K0, Cp_tmp * CP0, PRTMP, &
        !                        ufhfd_RA(J, 1:3) * U0 * D0 * CP0 * T0 + U0 * D0 * H0, uffhffd_FA(J, 1:3) * U0 * D0 * CP0 * T0 + U0 * D0 * H0, &
        !                        - DTDLKxztL_F0_io(J, 1:3) * CTHECD* T0/ L0 * K0, &
        !                        du_per* D0 * U0, Dv_per* D0 * U0, Dw_per* D0 * U0, Dh_per* D0 * U0
        !            END IF
    END DO
    DO N = 1, NMax
        CLOSE(TECFLG_FavAG(N))
    END DO

    !=========================budgetS =============================================================
    STFASTRESS(1) = 'THFx'
    STFASTRESS(2) = 'THFy'
    STFASTRESS(3) = 'THFz'
    !        scaling1 = CP0 * T0 * D0 * U0 * U0/ L0
    !        scaling2 = H0 * D0 * U0 * U0/ L0
    !        scaling3 = U0 * K0 * T0/ L0/ L0
    !        scaling4 = CP0 * T0 * M0 * U0/ L0/ L0
    !        scaling5 = H0 * M0 * U0/ L0/ L0
    !        scaling6 = CP0 * T0 * D0/ L0
    !        scaling7 = H0 * D0/ L0
    DO L = 1, 3

        TECFLG_FavAG(1) = TECFLG_FavAG(1) + 1

        !DO N = 1, NMax
        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.BUDGETS.Favre.' // TRIM(STFASTRESS(L)) //  &
        '.' // TRIM(PNTIM) // '.plt'
        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (26 variables)" '

        J = 0;                          WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") 'variables = '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Tw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Hw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Qw", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Cpw", '

        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_stres", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_enthg", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'Turbu_diffu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'dphDX_diffu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'pdhDX_stran", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'ConHF_diffu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'ConHF_dissp", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_diffu", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_dissp", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'press_accel", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'ConHF_accel", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_accel", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_Dgbdf", '
        J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG_FavAG(1), '(A)'          ) '"' // SJ2 // 'balance"'
        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        !END DO

        DO J = 1, NCL2

            WRITE(TECFLG_FavAG(1), '(28ES20.12)') &
            TauwSD(J), &
            DensSD(J), DenAvew, D1xztL_F0_io(J), &
            viscsD(J), VisAvew, M1xztL_F0_io(J), &
            F_A, YCC(J), YWdISD(J), &
            TwSD(J), HwSD(J), &
            QwSD(J), CpSD(J), &
            Budg_prodc_stres_thf(J, L), Budg_prodc_enthg_thf(J, L), &
            Budg_Turbu_diffu_thf(J, L), Budg_DphDX_diffu_thf(J, L), &
            Budg_pdhDX_stran_thf(J, L), Budg_ConHF_diffu_thf(J, L), &
            Budg_ConHF_dissp_thf(J, L), Budg_viscs_diffu_thf(J, L), &
            Budg_viscs_dissp_thf(J, L), Budg_press_accl1_thf(J, L), &
            Budg_ConHF_accel_thf(J, L), Budg_viscs_accl1_thf(J, L), &
            Budg_prodc_gvfc2_thf(J),   Budg_Balance1_thf(J, L)
            !                IF(iPPDimension == 1 ) THEN
            !                    WRITE(TECFLG_FavAG(2), '(17ES20.12)') YCC(J) * L0, (1.0_WP - DABS(YCC(J))) * REN * COE, COE* U0, &
            !                        Budg_prodc_stres_thf(J, L) * SCaling1 +scaling2, Budg_prodc_enthg_thf(J, L) * SCaling1 +scaling2, &
            !                        Budg_Turbu_diffu_thf(J, L) * SCaling1 +scaling2, Budg_press_accl1_thf(J, L) * SCaling1 +scaling2, &
            !                        Budg_DphDX_diffu_thf(J, L) * SCaling1 +scaling2, Budg_pdhDX_stran_thf(J, L) * SCaling1 +scaling2, &
            !                        Budg_ConHF_accel_thf(J, L) * SCaling3, &
            !                        Budg_ConHF_diffu_thf(J, L) * SCaling3, &
            !                        Budg_ConHF_dissp_thf(J, L) * SCaling3, &
            !                        Budg_viscs_accl1_thf(J, L) * SCaling4+scaling5, Budg_viscs_diffu_thf(J, L) * SCaling4+scaling5, &
            !                        Budg_viscs_dissp_thf(J, L) * SCaling4+scaling5, &
            !                        Budg_prodc_gvfc2_thf(J) * SCaling6+scaling7, Budg_Balance1_thf(J, L) * SCaling4+scaling5
            !                END IF
        END DO
        CLOSE(TECFLG_FavAG(1))
    END DO



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_HeatTransfer_Table_XZ_io
    !=================================================================
    !====== Defination based on Joong Hun Bae, 2005===========
    !=== MDOt = mass flow rate
    !=== Mbuk = bulk mass flux
    !=== HDOt = enthalpy flow rate
    !=== Hbuk = bulk enthalpy
    !==
    !=== Tbuk = bulk temperature
    !=== Dbuk = bulk DENSITY
    !=== Mbuk = bulk Viscousity
    !=== Kbuk = bulk thermal conductivity
    !===Cpbk = bulk Cp
    !===
    !===Ubuk = bulk velocity
    !=== Rebk = bulk Reynolds no.
    !===Prbk = bulk Prandtl no.
    !=== Nubk = bulk Nusselt number
    !===Bobk = Boyancy PARAMETER
    !===Cfbk = local Friction coefficient

    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    CHARACTER(128) :: FLNM
    INTEGER(4) :: I, J, IP, JP, JM, K, N
    REAL(WP) :: Hw, Tw, Dw, Cpw, Mw, Kw
    REAL(WP) :: H_tmp, D_tmp, D_tmp1, D_tmp2, T_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp, Ptemp, DH_tmp
    REAL(WP) :: T_tmp1, T_tmp2, K_tmp1, K_tmp2, Ltemp(2)
    REAL(WP) :: Buoya_1, Buoya_2, hc_Dwve(2), NubkTsdave(2), Nubk_DTw_ave
    INTEGER(4) :: TECFLG = 200
    REAL(WP) :: DUDYL, DUDYU
    REAL(WP) :: RTMPmin, TempDiff
    REAL(WP) :: CMWORK, CMWORK_D ! work DOne by effective pressure gradient
    REAL(WP) :: dqw_D, dqw
    REAL(WP) :: Buoyasd_1(2), Buoyasd_2(2)
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HD
    REAL(WP) :: spline_interpolation_HM
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: spline_interpolation_HB


    !===================================================================
    !++++++++++++++++sided bulk valueS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL PP_J4TbukTsd
    !WRITE(*, *) '# L4TbkSsd',L4TbkSsd(1:2)
    !WRITE(*, *) '# Ldist_io',Ldist_io(1:2)
    !WRITE(*, *) '# DbukSsd  ', DbukSsd(1:2)
    !WRITE(*, *) '# MbukSsd  ', MbukSsd(1:2)

    !=========================================================================
    TECFLG = 200
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    OPEN(TECFLG, FILE = TRIM(FilePath4) // 'Result.IO.Table.WallandBulk.' // TRIM(PNTIM) // '.plt')!, POSITION = 'APPEND')

    WRITE(TECFLG, '(A)         ') '%######## Reference States (dimensionaless)####################'
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF Re = ', REN
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF PRT0 = ', PRT0
    WRITE(TECFLG, '(A)') ''
    WRITE(TECFLG, '(A)         ') '%######## Reference States (dimensional)  #####################'
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF L0(m) = ', L0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF U0(M/s) = ', U0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF MDOt(Kg/m2s) = ', U0 * D0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF P0(Pa) = ', P0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF T0(K) = ', T0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF T0(C) = ', T0- TK2C
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF H0(J) = ', H0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF D0(Kg/M3) = ', D0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF K0(W/M-K) = ', K0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF M0(Pa-s) = ', M0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF Cp0(J/Kg/K) = ', CP0
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF B0(1/K) = ', B0
    WRITE(TECFLG, '(A, 1ES17.9)') 'Gravity/U0^2 = ', F_A
    WRITE(TECFLG, '(A)') ''

    !==========wall PARAMETERS ==================================================
    HWAL_RA_D(:) = HWAL_RA(:) * T0 * CP0 + H0
    Dwal_D(:) = Dwal(:) * D0
    Twal_D(:) = Twal(:) * T0
    Mwal_D(:) = Mwal(:) * M0
    Kwal_D(:) = Kwal(:) * K0
    Cpwal_D(:) = CPwal(:) * CP0

    Ttau(1:2)       = DABS(Qw(1:2) / (Dwal(1:2) * Cpwal(1:2) * Utaw_io(1:2)))!Ttau_D(1:2) / T0
    Ttau_D(1:2)     = DABS(qw_D(1:2) / (Dwal_D(1:2) * Cpwal_D(1:2) * Utaw_io(1:2) * U0))

    WRITE(TECFLG, '(A)         ') '%########variables (Wall)######################################'
    WRITE(TECFLG, '(A, 24X, 2(A,6X), 2(A, 4X))') '%#', 'dimensional', 'dimensional', 'dimensionless', 'dimensionless'
    WRITE(TECFLG, '(A, 24X, 4(A,6X))') '%#', 'wall(Y = -1) ', 'wall(Y =+ 1) ', 'wall(Y = -1) ', 'wall(Y =+ 1) '
    WRITE(TECFLG, '(A, 4ES17.9)') 'Tw: Temperature (K) = ', Twal_D(1:2), Twal(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') 'Tw: Temperature (C) = ', Twal_D(1:2) - TK2C
    WRITE(TECFLG, '(A, 4ES17.9)') 'Tt: Friction Ttau (K) = ', Ttau_D(1:2), Ttau(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') 'Dw: DENSITY (Kg/ M3)  = ', Dwal_D(1:2), Dwal(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') 'Hw: Enthalpy(J / Kg)  = ', Hwal_RA_D(1:2), Hwal_RA(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') 'Mw: Viscousity(Pa-s)  = ', Mwal_D(1:2), Mwal(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') 'Kw: Conductivity(W/M-K) = ', Kwal_D(1:2), Kwal(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') 'Cpw: Specific heat(Cp) = ', Cpwal_D(1:2), Cpwal(1:2)
    WRITE(TECFLG, '(A)') ''

    !============find the location of TcP =========================
    RTMPmin = 1.0E+10_WP
    IF(DMAX1(Twal(1), Twal(2)) >  T4CpMax) THEN
        DO J = 1, NCL2
            TempDiff = DABS(T1xztL_F0_io(J) - T4CpMax)
            IF(TempDiff < RTMPmin) THEN
                RTMPmin = TempDiff
                J4Tpc = J
            END IF
        END DO
    ELSE
        J4Tpc = 0
    END IF

    !===============bulk valueS =======================================================================
    !==================calculate bulk mass flux, and bulk enthalpY =========
    MDOt = 0.0_WP
    DHDOt = 0.0_WP
    DO J = 1, NCL2
        DO K = 1, NCL3
            MDOt = MDOt + G1xztL_F0_io(J, 1) / RCCI1(J) / DYFI(J) / DZI
            DHDOt = DHDOt + DHxztL_F0_io(J) / RCCI1(J) / DYFI(J) / DZI
        END DO
    END DO

    Gbuk = MDOt /Area_inlet
    DHbuk = DHDOt /Area_inlet

    CALL THERM_PROP_UPDATE_FROM_DH(DHbuk, Hbuk, Tbuk, Dbuk, Mbuk, Kbuk, Cpbk, Bbuk)

    Ubuk  = Gbuk/ Dbuk

    Gbuk_D = Gbuk * D0 * U0
    Hbuk_D = Hbuk * T0 * CP0 + H0
    Ubuk_D = Ubuk* U0

    Dbuk_D  = Dbuk* D0
    Tbuk_D  = Tbuk* T0
    Mbuk_D  = Mbuk* M0
    Kbuk_D  = Kbuk* K0
    Cpbk_D  = CPbk* CP0
    Bbuk_D  = Bbuk*B0

    WRITE(TECFLG, '(A)         ') '%########Case Set up##########################################'
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        WRITE(TECFLG, '(A, A)      ') 'Wall thermal conditions        = ', ' Constant Wall temperature(K), (C)'
        WRITE(TECFLG, '(A, 2ES17.9)') '  Wall temperature (K) & (C) at (Y = -1) = ', Twal_D(1), Twal_D(1) - TK2C
        WRITE(TECFLG, '(A, 2ES17.9)') '  Wall temperature (K) & (C) at (Y =+ 1) = ', Twal_D(2), Twal_D(2) - TK2C
        WRITE(TECFLG, '(A, 1ES17.9)') 'Constant Mass Flux (Kg/ M2 -s)       = ', Gbuk_D
        WRITE(TECFLG, '(A)') ''
    END IF

    WRITE(TECFLG, '(A)         ') '%########variables (Bulk)#####################################'
    WRITE(TECFLG, '(A, 24X, 4(A,6X))') '%#', 'dimensional', 'dimensionless'
    WRITE(TECFLG, '(A, 2ES17.9)') 'MDOt:Mass flux(Kg/ M2 -s)     = ', Gbuk_D, Gbuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Ub  :Velocity (M /s)         = ', Ubuk_D, Ubuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Hb  :Enthalpy(J / Kg)         = ', Hbuk_D, Hbuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Tb: Temperature (K)         = ', Tbuk_D, Tbuk
    WRITE(TECFLG, '(A, 1ES17.9)') 'Tb: Temperature (C)         = ', Tbuk_D - TK2C
    WRITE(TECFLG, '(A, 2ES17.9)') 'Db: DENSITY (Kg/ M3)         = ', Dbuk_D, Dbuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Mub:Viscousity(Pa-s)        = ', Mbuk_D, Mbuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Kb:Conductivity(W/ M - K)      = ', Kbuk_D, Kbuk
    WRITE(TECFLG, '(A, 2ES17.9)') 'Cpb:Specific Heat(Cp)       = ', Cpbk_D, Cpbk
    WRITE(TECFLG, '(A)') ''


    !=============== search === TablE ==========
    !================================
    Rebk = Dbuk* Ubuk* 2.0_WP * ALX2 / Mbuk * REN ! based on diameter or the whole hEIGht of Channel.
    Prbk = Mbuk* Cpbk/ Kbuk      * Prt0

    WRITE(TECFLG, '(A)         ') '%########variables (Bulk), SSzero Y = -1 side, Y = 1 side dimensionless######################'
    WRITE(TECFLG, '(A, 3ES17.9)') 'Reynolds No                 = ', Rebk, RebkSsd(1:2)
    WRITE(TECFLG, '(A, 3ES17.9)') 'Prandlt No                  = ', Prbk, PrbkSsd(1:2)
    !====================================================================

    IF(iThermalWallType(iTopWall) == BC_Fixed_Heat_Flux .AND. iThermalWallType(iBotWall) == BC_Fixed_Heat_Flux) THEN
        qw_D(1:2) = thermalWallBC_Dim(1:2)
        !====convective heat transfer ratE ==== HC = Qw/ (Tw- Tb) ====
        hc_D(1:2) = qw_D(1:2) / (Twal_D(1:2) - Tbuk_D)
        hc_Dwve(1:2) = 0.5_WP * (DABS(qw_D(1)) + DABS(qw_D(2))) / (Twal_D(1:2) - Tbuk_D)
        !====bulk Grashof numbeR ===========

        Grbk(1:2) = G_A*Bbuk_D* Qw_D(1:2) * (L0**4) * Dbuk_D* Dbuk_D/ Mbuk_D/ Mbuk_D/ Kbuk_D
        !WRITE(*, *) '# Grbk', Grbk(1:2)
    END IF
    IF(iThermalWallType(iTopWall) == BC_Fixed_Temperature .AND. iThermalWallType(iBotWall) == BC_Fixed_Temperature) THEN
        !====wall heat fluX ========== here IS not corret, as bAR over (kdT/ Dy) /= kdbART/ Dy
        qw_D(1) = Qw(1) * D0 * U0 * CP0 * T0 !Kwal_D(1) * (Twal_D(1) - T1xztL_F0_io(1) * T0   ) / ((YCC(1) - YND(1)) * L0)
        qw_D(2) = Qw(2) * D0 * U0 * CP0 * T0 !Kwal_D(2) * (Twal_D(2) - T1xztL_F0_io(NCL2) * T0) / ((YND(NND2) - YCC(NCL2)) * L0)

        !test
        !!WRITE(*, *) 'qwc', qw_D(1), Kwal_D(1) * (Twal_D(1) - T1xztL_F0_io(1) * T0   ) / ((YCC(1) - YND(1)) * L0)       !Same
        !!WRITE(*, *) 'qwh', qw_D(2), Kwal_D(2) * (Twal_D(2) - T1xztL_F0_io(NCL2) * T0) / ((YND(NND2) - YCC(NCL2)) * L0) !Same

        qw_D_ave = 0.5_WP * (DABS(qw_D(1)) + DABS(qw_D(2)))
        !====convective heat transfer ratE ==== HC = Qw/ (Tw- Tb) ====
        hc_DTw_D(1:2) = qw_D(1:2) / DABS(Twal_D(1) - Twal_D(2))
        hc_DTw_D_ave = qw_D_ave / DABS(Twal_D(1) - Twal_D(2))

        hc_D(1:2) = qw_D(1:2) / DABS(Twal_D(1:2) - Tbuk_D)
        hc_Dwve(1:2) = DABS(qw_D_ave / DABS(Twal_D(1:2) - Tbuk_D))



        !====bulk Grashof numbeR ===========
        Grbk(1) = G_A*Bbuk_D* DABS(Twal_D(1) - Twal_D(2)) * ((2.0_WP * L0)**3) * Dbuk_D* Dbuk_D/ Mbuk_D/ Mbuk_D
        !WRITE(*, *) '# Grbk', Grbk(1)
    END IF

    !=============== Nasselt numbeR ==============
    Nubk_DTw(1:2)   = 2.0_WP * L0 * Hc_DTw_D(1:2) / Kbuk_D  ! based on q_w(h,c) / (Tc- Th)
    Nubk_DTw_ave    = 2.0_WP * L0 * Hc_DTw_D_ave/ Kbuk_D   ! based on averaged q_w(H +c) / (Tc- Th)

    NubkTsd(1:2)       = L4Tbk(1:2) * Hc_D(1:2) / KbukTsd(1:2) / K0    ! based on q_w(h,c) / (Tw- Tb)
    NubkTsdave(1:2)    = L4Tbk(1:2) * Hc_Dwve(1:2) / KbukTsd(1:2) / K0  ! based on averaged q_w(H +c) / (Tw- Tb)



    !============work DOne by sheAR stresS =======================
    !CMWORK_D = DABS( (2.0_WP * Tauw_D_ave_io) / (DenAvew* D0) * (U0 * D0) )
    CMWORK_D = Tauw_D_ave_io * Ubuk_D
    dqw_D  = DABS( DABS(qw_D(1)) - DABS(qw_D(2)))

    !CMWORK = DABS( (2.0_WP * Tauw_ave_io) / DenAvew * Gbuk)
    CMWORK = Tauw_ave_io * Ubuk
    dqw    = DABS( DABS(Qw(1)) - DABS(Qw(2))) * (T0 * Cp0/ U0/ U0)

    WRITE(TECFLG, '(A)         ') '%########Energy Conservation#############################'
    WRITE(TECFLG, '(A, 24X, 4(A,6X))') '%#', 'dimensionless', 'dimensional'
    WRITE(TECFLG, '(A, 2ES19.9)') 'Net Heat Flux =          ', dqw,    dqw_D
    WRITE(TECFLG, '(A, 2ES19.9)') 'Work DOne by (dP / DX)_eff  = ', CMWORK, CMWORK_D
    WRITE(TECFLG, '(A, 2ES19.9)') 'Above two dIFf = ', dqw-CMWORK, dqw_D -CMWORK_D
    WRITE(TECFLG, '(A)') '  '

    WRITE(TECFLG, '(A)         ') '%########Heat transfer parameters#############################'
    WRITE(TECFLG, '(A, 30X, 2(A,7X), A, 1X, A)') '%#', 'wall(Y = -1)', 'wall(Y =+ 1)', &
    'Two-wall average', 'Error(\%) = (W1-W2) / (W1 + W2) *100'
    WRITE(TECFLG, '(A, 4ES19.9)') 'Qw: Wall heat flux(undim)   = ', Qw(1:2), 0.5_WP * (Qw(1) + Qw(2)), &
    DABS(DABS(Qw(1)) - DABS(Qw(2))) / (Qw(1) + Qw(2)) *100.0_WP
    WRITE(TECFLG, '(A, 4ES19.9)') 'Qw: Wall heat flux(W/ M2)    = ', qw_D(1:2),    qw_D_ave, &
    0.5_WP * DABS(DABS(qw_D(1)) - DABS(qw_D(2))) /qw_D_ave*100.0_WP
    WRITE(TECFLG, '(A)') '  '

    !=========== INTegrated DENSITY = And Nussult numbeR =========
    D_inT = 0.0_WP
    Nu_inT = 0.0_WP
    DO J = 1, NND2
        IF (J == 1) THEN
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            D_int  = D_int  + ( T1xztL_F0_io(J) - Twal(1) ) * D_tmp
            Nu_int = Nu_int + YCC(J) / (K_tmp + Kwal(1)) * (YCC(J) - 0.0_WP)
        ELSE IF(J == NND2) THEN
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J - 1), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            D_tmp = Dwal(2)
            D_int = D_int + ( Twal(2) - T1xztL_F0_io(J - 1) ) * D_tmp
            Nu_int = Nu_int + (YCC(J - 1) + YND(J)) / (K_tmp + Kwal(2)) * (YND(J) - YCC(J - 1))
        ELSE
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J - 1), D_tmp1, M_tmp, K_tmp1, Cp_tmp, B_tmp)
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp2, M_tmp, K_tmp2, Cp_tmp, B_tmp)
            D_tmp = D_tmp2
            D_int  = D_int + ( T1xztL_F0_io(J) - T1xztL_F0_io(J - 1) ) * D_tmp
            Nu_int = Nu_int + YND(J) / ((YCL2ND_WFF(J) * K_tmp2 + YCL2ND_WFB(J) * K_tmp1)) * ((YCC(J) - YCC(J - 1)))
        END IF
    END DO
    Nu_int = Nu_int/ (2.0_WP) * L0/ K0 * Qw_D_ave/ ((Twal(2) - Twal(1)) * T0)
    !WRITE(*, *) '# Nu_int', Nu_int

    IF(DABS(Twal(2) - Twal(1)) < 1.0E-12_WP) THEN
        D_int = D_tmp
    ELSE
        D_int = D_int/ (Twal(2) - Twal(1))
    END IF
    !WRITE(*, *) '# D_int', D_int
    D_int_D = D_int * D0

    Grbk_Drho = (Dbuk_D - D_int_D) / Dbuk_D * G_A * (L0**3) * Dbuk_D * Dbuk_D / Mbuk_D / Mbuk_D
    !WRITE(*, *) '# Grbk_Drho', Grbk_Drho

    !====================================================================
    Buoya_1 = Grbk(1) / (Rebk**3.425_WP) * (Prbk**0.8_WP)
    Buoya_2 = Grbk(1) / (Rebk**2.7_WP)
    !WRITE(*, *) '# Buoya_12', Buoya_1, Buoya_2



    Dsd_inT = 0.0_WP
    DO J = 1, J4TbukSsd(1)
        IF (J == 1) THEN
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            Dsd_int(1) = Dsd_int(1)  + ( T1xztL_F0_io(J) - Twal(1) ) * D_tmp
        ELSE IF(J == NND2) THEN
            D_tmp = Dwal(2)
            Dsd_int(1) = Dsd_int(1) + ( Twal(2) - T1xztL_F0_io(J - 1) ) * D_tmp
        ELSE
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            Dsd_int(1) = Dsd_int(1) + ( T1xztL_F0_io(J) - T1xztL_F0_io(J - 1) ) * D_tmp
        END IF
    END DO

    IF(DABS(Twal(2) - Twal(1)) < 1.0E-12_WP) THEN
        Dsd_int(1) = D_tmp
    ELSE
        Dsd_int(1) = Dsd_int(1) / (T1xztL_F0_io(J4TbukSsd(1)) - Twal(1))
    END IF


    DO J = J4TbukSsd(2), NND2
        IF (J == 1) THEN
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            Dsd_int(2) = Dsd_int(2)  + ( T1xztL_F0_io(J) - Twal(1) ) * D_tmp
        ELSE IF(J == NND2) THEN
            D_tmp = Dwal(2)
            Dsd_int(2) = Dsd_int(2) + ( Twal(2) - T1xztL_F0_io(J - 1) ) * D_tmp
        ELSE
            CALL THERM_PROP_UPDATE_FROM_T(DH_tmp, H_tmp, T1xztL_F0_io(J), D_tmp, M_tmp, K_tmp, Cp_tmp, B_tmp)
            Dsd_int(2) = Dsd_int(2) + ( T1xztL_F0_io(J) - T1xztL_F0_io(J - 1) ) * D_tmp
        END IF
    END DO

    IF(DABS(Twal(2) - Twal(1)) < 1.0E-12_WP) THEN
        Dsd_int(2) = D_tmp
    ELSE
        Dsd_int(2) = Dsd_int(2) / (Twal(2) - T1xztL_F0_io(J4TbukSsd(2)))
    END IF

    !WRITE(*, *) '# Dsd_int', Dsd_int(1:2)

    !Ltemp(1:2) = Ldist_io(1:2)
    Ltemp(1:2) = Ldist_io(1:2)

    GrbkSsd(1:2) = G_A * ((Ltemp(1:2) * L0)**3) * (DbukSsd(1:2) * D0 * DbukSsd(1:2) * D0) / &
    (MbukSsd(1:2) * M0 * MbukSsd(1:2) * M0) * &
    (DbukSsd(1:2) - Dsd_int(1:2)) / DbukSsd(1:2)
    !!WRITE(*, *) '# GrbkSsd', GrbkSsd(1:2)

    hcsd_D(1:2) = qw_D(1:2) / (Twal_D(1:2) - TbukSsd(1:2) * T0)
    NubkSsd(1:2) = Ltemp(1:2) * L0 * Hcsd_D(1:2) / (KbukSsd(1:2) * K0)    ! based on q_w(h,c) / (Tw- Tb)
    !!WRITE(*, *) '# NubkSsd', NubkSsd(1:2)

    Buoyasd_1(1:2) = GrbkSsd(1:2) / ((RebkSsd(1:2) / Ldist_io(1:2) * Ltemp(1:2))**3.425_WP) * (PrbkSsd(1:2)**0.8_WP)
    Buoyasd_2(1:2) = GrbkSsd(1:2) / ((RebkSsd(1:2) / Ldist_io(1:2) * Ltemp(1:2))**2.7_WP)
    !WRITE(*, *) '# Buoyasd_1', Buoyasd_1(1:2)
    !WRITE(*, *) '# Buoyasd_2', Buoyasd_2(1:2)
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    WRITE(TECFLG, '(A)         ') '%########integral, dimensionless(c,h), dimensinoal(c,h) ######################'
    WRITE(TECFLG, '(A, 2ES17.9)') 'Nu:      Integral               = ', Nu_int
    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A)         ') '%########sided bulk average(Zero sheAR stress balance)#######'
    WRITE(TECFLG, '(A, 2X, I3.1, 2ES17.9)')'Jindex, Ycc, = ', J4SS0, YCC(J4SS0)
    WRITE(TECFLG, '(A, 15X, 2(A,7X), A)')   '%#', 'side(Y = -1)', 'side(Y =+ 1)'
    WRITE(TECFLG, '(A, 3ES17.9)')       'Y+max       = ', (YCC(J4SS0) + 1.0_WP) * Ret_io(1) / Ldist_io(1), &
    (1.0_WP - YCC(J4SS0)) * Ret_io(2) / Ldist_io(2)
    WRITE(TECFLG, '(A, 2ES17.9)')       'Reynolds tau  each side:', Ret_io(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)')       'Reynolds bulk each side:', RebkSsd(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)')       'Velocity bulk each side:', UbukSsd(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)')       'ChARac L bulk each side:', Ldist_io(1:2)
    WRITE(TECFLG, '(A)') '  '

    WRITE(TECFLG, '(A)         ') '%########dimensional (1)on the cold wall, (2)on the hot wall, (3)averaged (cold, hot)###'
    WRITE(TECFLG, '(A, 3ES17.9)') 'Qw: Wall heat flux(W/ M2)          = ', qw_D(1:2),    qw_D_ave
    WRITE(TECFLG, '(A         )') '1:==based on whole domain L === '
    WRITE(TECFLG, '(A, 2ES17.9)') '1:Dint:    Integral DENSITY(Kg/ M3) = ', D_int, D_int_D
    WRITE(TECFLG, '(A, 3ES17.9)') '1:hC = Qw/ (Th- Tc) (W/ M2 / K)          = ', hc_DTw_D(1:2), hc_DTw_D_ave
    WRITE(TECFLG, '(A, 3ES17.9)') '1 : Nu= HC *2R/ Kb                     = ', Nubk_DTw(1:2), Nubk_DTw_ave
    WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb_Drho                        = ', Grbk_Drho
    WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb/ Re2.7                       = ', Buoya_2
    WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb/ Re3.4/Pr0.8                 = ', Buoya_1
    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A         )') '2:==based on Tbulk  sided=====   '
    WRITE(TECFLG, '(A, I3.1, 1ES17.9)')'2:Jindex, Ycc                      = ', J4Tbk, YCC(J4Tbk)
    WRITE(TECFLG, '(A, 2ES17.9)') '2:Ldist:                          = ', L4TbkTsd(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') '2:hC = Qw/ (T_{c/h}- Tb) (W/ M2 / K)     = ', hc_D(1:2),    hc_Dwve(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') '2:Nu= HC * Lt/ Kbt(Tbulk sep)         = ', NubkTsd(1:2), NubkTsdave(1:2)
    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A         )') '3:==based on SSzero sided=====   '
    WRITE(TECFLG, '(A, I3.1, 1ES17.9)')'3:J4SS0, Ycc                       = ', J4SS0, YCC(J4SS0)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Ldist:                          = ', Ldist_io(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Ldist to Tbulk each side        = ', L4TbkSsd(1:2)
    WRITE(TECFLG, '(A, 2I3.1, 2ES17.9)')'3:J4TbukSsd, Ycc                   = ', J4TbukSsd(1:2), YCC(J4TbukSsd(1:2))
    WRITE(TECFLG, '(A, 2ES17.9)') '3:TbukSsd:                        = ', TbukSsd(1:2)
    WRITE(TECFLG, '(A, 4ES17.9)') '3:DintSsd: Integral DENSITY(Kg/ M3) = ', Dsd_int(1:2), Dsd_int(1:2) * D0
    WRITE(TECFLG, '(A, 2ES17.9)') '3:hC = Qw/ (T_{wc/h}- T{bc/h)  (W/ M2 / K) = ', hcsd_D(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Nu= HC * Lt/ Kbt(SS0 sep)           = ', NubkSsd(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb_Drho(lc,rhobc)              = ', GrbkSsd(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb/ Re2.7                       = ', Buoyasd_2(1:2)
    WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb/ Re3.4/Pr0.8                 = ', Buoyasd_1(1:2)


    WRITE(TECFLG, '(A)') '  '



    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_Cf_Table_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: TECFLG  = 200
    INTEGER(4) :: TECFLG2 = 202
    INTEGER(4) :: I, J
    REAL(WP) :: DUDY(2)
    REAL(WP) :: dxPlus(2), dzPlus(2), yMinPlus(2), yMaxPlus(2), dyMin, dyMax
    REAL(WP) :: DENtemp, ddenintg, DeNMintg

    REAL(WP), ALLOCATABLE :: dxPlusj(:, :), dzPlusj(:, :), dyplusj(:, :)
    REAL(WP) :: INNerL(2)
    REAL(WP) :: Tau_diff, Tau_avag, Llocal(2)

    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    OPEN(TECFLG, FILE = TRIM(FilePath4) // 'Result.IO.Table.WallandBulk.' // TRIM(PNTIM) // '.plt', POSITION = 'APPEND')

    Tauw_io(1:2) = DABS(Tauw_io(1:2))
    Cf0_io(1:2) = 2.0_WP * DABS(Tauw_io(1:2))
    IF(iThermoDynamics == 0) THEN
        Cfbk_io(1:2) = Cf0_io(1:2)
    END IF
    IF(iThermoDynamics == 1) THEN
        Cfbk_io(1:2) = 2.0_WP * Tauw_io(1:2) / Dbuk/ Ubuk/ Ubuk
    END IF

    Cf0_ave_io = 0.5_WP * (Cf0_io(1)   + Cf0_io(2) )

    INNerL(1:2) = 1.0_WP / Ret_io(1:2)

    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A)         ') '%########Wall Coefficient ####################################'
    WRITE(TECFLG, '(A, 21X, 2(A,7X), A, 1X, A)') '%#', 'wall(Y = -1)', 'wall(Y =+ 1)', 'Two walls average', &
    'Error(\%) = (W1-W2) / (W1 + W2) *100'
    WRITE(TECFLG, '(A, 4ES17.9) ') 'Cf0                = ',Cf0_io(1),  Cf0_io(2),   Cf0_ave_io, &
    DABS(Cf0_io(1) -Cf0_io(2)) / (Cf0_io(1) +Cf0_io(2)) *100.0_WP
    WRITE(TECFLG, '(A, 3ES17.9) ') 'Cfbk               = ',Cfbk_io(1), Cfbk_io(2),  Cfbk_ave_io
    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A)         ') '%########Wall SheAR, dimensionless ###########################'
    WRITE(TECFLG, '(A, 21X, 2(A,7X), A)') '%#', 'wall(Y = -1)', 'wall(Y =+ 1)'
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Tau                = ', Tauw_io(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Utau               = ',Utaw_io(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Re_tau             = ',Ret_io(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Delta_nu(iNNer L) = ', INNerL(1:2)
    WRITE(TECFLG, '(A)') '  '

    IF(iCase == ICHANL) THEN
        dyMin = 1.0_WP / DYFI(1)
        dyMax = 1.0_WP / DYFI(NCL2 / 2)
    ELSE
        dyMin = 1.0_WP / DYFI(NCL2)
        dyMax = 1.0_WP / DYFI(1)
    END IF

    dxPlus(1:2) = DX * Ret_io(1:2) / Ldist_io(1:2)
    dzPlus(1:2) = DZ * Ret_io(1:2) / Ldist_io(1:2)

    yMinPlus(1:2) = dyMin * Ret_io(1:2) / Ldist_io(1:2)
    yMaxPlus(1:2) = dyMax * Ret_io(1:2) / Ldist_io(1:2)


    WRITE(TECFLG, '(A)         ') '%########Real mesh resolution ###############################'
    WRITE(TECFLG, '(A, 4X, 2(A,7X), A, 1X, A)') '%#Based on Utau at', ' wall(Y = -1)', 'wall(Y =+ 1), and original'
    WRITE(TECFLG, '(A, 2ES17.9) ') 'DX +                = ', dxPlus(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'dZ +                = ', dzPlus(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'dy1 +               = ', yMinPlus(1:2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'dy1maX +            = ', yMaxPlus(1:2)

    IF(iThermoDynamics == 1) THEN

        ALLOCATE (dxPlusj(NCL2, 3))
        ALLOCATE (dyplusj(NCL2, 3))
        ALLOCATE (dzPlusj(NCL2, 3))
        OPEN(TECFLG2, FILE = TRIM(FilePath4) // 'Result.IO.SemiLocal.Yplus.' // TRIM(PNTIM) // '.plt')
        WRITE(TECFLG2, '(A)') 'TITLE = " Local Yplus (26 variables)" '

        WRITE(TECFLG2, '(A)') 'variables = "YCC", "DXN +", "dyN +", "dzN +", "DXP +", "dyhP +", "dzP +", "DXa+", "dya+", "dza+" '
        WRITE(TECFLG2, '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        DO J = 1, NCL2

            dxPlusj(J, 1) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(1) ) / M1xztL_F0_io(J) * DX
            dyplusj(J, 1) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(1) ) / M1xztL_F0_io(J) / DYFI(J)
            dzPlusj(J, 1) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(1) ) / M1xztL_F0_io(J) * DZ

            dxPlusj(J, 2) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(2) ) / M1xztL_F0_io(J) * DX
            dyplusj(J, 2) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(2) ) / M1xztL_F0_io(J) / DYFI(J)
            dzPlusj(J, 2) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_io(2) ) / M1xztL_F0_io(J) * DZ

            dxPlusj(J, 3) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_ave_io ) / M1xztL_F0_io(J) * DX
            dyplusj(J, 3) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_ave_io ) / M1xztL_F0_io(J) / DYFI(J)
            dzPlusj(J, 3) = REN * DSQRT(D1xztL_F0_io(J) * Tauw_ave_io ) / M1xztL_F0_io(J) * DZ
            WRITE(TECFLG2, '(10ES17.9)') YCC(J), dxPlusj(J, 1), dyplusj(J, 1), dzPlusj(J, 1), &
            dxPlusj(J, 2), dyplusj(J, 2), dzPlusj(J, 2), &
            dxPlusj(J, 3), dyplusj(J, 3), dzPlusj(J, 3)

        END DO
        CLOSE(TECFLG2)

        !            dxPlus(1:2) = REN * DSQRT(DWAL(1:2) * Tauw_io(1:2) ) / MWal(1:2) * DX
        !            dzPlus(1:2) = REN * DSQRT(DWAL(1:2) * Tauw_io(1:2) ) / MWal(1:2) * DZ

        !            yMinPlus(1:2) = dyMin * REN * DSQRT(DWAL(1:2) * Tauw_io(1:2) ) / MWal(1:2)
        !            yMaxPlus(1:2) = dyMax * REN * DSQRT(DWAL(1:2) * Tauw_io(1:2) ) / MWal(1:2)

        !            WRITE(TECFLG, '(A)         ') '%########Real mesh resolution ###############################'
        !            WRITE(TECFLG, '(A, 4X, 2(A,7X), A, 1X, A)') '%#Based on Utau at', ' wall(Y = -1)', 'wall(Y =+ 1), and original'
        !            WRITE(TECFLG, '(A, 2ES17.9) ') 'DX +                = ', dxPlus(1:2)
        !            WRITE(TECFLG, '(A, 2ES17.9) ') 'dZ +                = ', dzPlus(1:2)
        !            WRITE(TECFLG, '(A, 2ES17.9) ') 'dy1 +               = ', yMinPlus(1:2)
        !            WRITE(TECFLG, '(A, 2ES17.9) ') 'dy1maX +            = ', yMaxPlus(1:2)


        WRITE(TECFLG, '(A, 4X, 2(A,7X), A, 1X, A)') '%#Based on Utau_w at', ' wall(Y = -1)', 'wall(Y =+ 1), and semI -local'
        WRITE(TECFLG, '(A)'               )              'min.    max.     min.  max.'
        WRITE(TECFLG, '(A, 4ES17.9) ') 'DX +                = ', MINVAL(dxPlusj(:, 1)), MAXVAL(dxPlusj(:, 1)), &
        MINVAL(dxPlusj(:, 2)), MAXVAL(dxPlusj(:, 2))
        WRITE(TECFLG, '(A, 4ES17.9) ') 'dZ +                = ', MINVAL(dzPlusj(:, 1)), MAXVAL(dzPlusj(:, 1)), &
        MINVAL(dzPlusj(:, 2)), MAXVAL(dzPlusj(:, 2))
        WRITE(TECFLG, '(A, 4ES17.9) ') 'dy+                = ', MINVAL(dyplusj(:, 1)), MAXVAL(dyplusj(:, 1)), &
        MINVAL(dyplusj(:, 2)), MAXVAL(dyplusj(:, 2))

        WRITE(TECFLG, '(A, 4X, 2(A,7X), A, 1X, A)') '%#Based on Utau_ave at', ' wall(Y = -1)', 'wall(Y =+ 1), and semI -local'
        WRITE(TECFLG, '(A)'               )              'min.    max.     min.  max.'
        WRITE(TECFLG, '(A, 2ES17.9) ') 'DX +                = ', MINVAL(dxPlusj(:, 3)), MAXVAL(dxPlusj(:, 3))
        WRITE(TECFLG, '(A, 2ES17.9) ') 'dZ +                = ', MINVAL(dzPlusj(:, 3)), MAXVAL(dzPlusj(:, 3))
        WRITE(TECFLG, '(A, 2ES17.9) ') 'dy+                = ', MINVAL(dyplusj(:, 3)), MAXVAL(dyplusj(:, 3))

    END IF


    !=============================================================
    ddenintg = 0.0_WP
    deNMintg = 0.0_WP
    bdfciNTG = 0.0_WP
    DensIntg = 0.0_WP
    DO J = 1, NCL2
        ! second order intgeral
        !            IF(J == 1) THEN
        !                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
        !                                  Dwal(1) )
        !            ELSE IF(J == NCL2) THEN
        !                DENtemp = 0.5_WP * ( Dwal(2) +             &
        !                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
        !            ELSE
        !                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
        !                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
        !            END IF
        !first order integral
        DENtemp = D1xztL_F0_io(J)

        deNMintg = deNMintg+ DENtemP / DYFI(J)
        DensIntg(J) = deNMintg !; !WRITE(*, *) J, DENtemp, DensIntg(J)

        ddenintg = ddenintg + (DENtemP - DenAvew) / DYFI(J)
        bdfcintg(J) = F_A* Ddenintg
    END DO
    CALL CHKHDL(' ==>Calculated bodyforce distribution', MYID)


    !        WRITE(TECFLG, '(A)') '  '
    !        IF(iThermoDynamics == 1) THEN
    !            WRITE(TECFLG, '(A)         ') '%########Wall SheAR, dimensional ########################'
    !            WRITE(TECFLG, '(A, 21X, 2(A,7X), A)') '%#', 'wall(Y = -1)', 'wall(Y =+ 1)'
    !            WRITE(TECFLG, '(A, 2ES17.9) ') 'Tau(Kg/ M /s^2)      = ', Tauw_D_io(1), Tauw_D_io(2)
    !            WRITE(TECFLG, '(A, 2ES17.9) ') 'Utau(M /s)          = ',Utaw_D_io(1), Utaw_D_io(2)
    !        END IF

    CLOSE(TECFLG)

    IF(iThermoDynamics == 1) THEN
        DEALLOCATE (dxPlusj)
        DEALLOCATE (dyplusj)
        DEALLOCATE (dzPlusj)
    END IF


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_Checking_TABLE_XZ_io
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: TECFLG = 200
    INTEGER(4) :: I



    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    OPEN(TECFLG, FILE = TRIM(FilePath4) // 'Result.IO.Table.WallandBulk.' // TRIM(PNTIM) // '.plt', POSITION = 'APPEND')

    WRITE(TECFLG, '(A)         ') '%########Reference States (dimensionaless)###################'
    WRITE(TECFLG, '(A, 1ES17.9)') 'REF Re                      = ', REN
    WRITE(TECFLG, '(A)') '  '
    WRITE(TECFLG, '(A)         ') '%########Checking force balance (wall unit scaled) ##########'
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Total Force in X direction (RA and FA) = ', NSFbalt_RA(1), NSFbalt_FA(1)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Total Force in Y direction (RA and FA) = ', NSFbalt_RA(2), NSFbalt_FA(2)
    WRITE(TECFLG, '(A, 2ES17.9) ') 'Total Force in Z direction (RA and FA) = ', NSFbalt_RA(3), NSFbalt_FA(3)
    WRITE(TECFLG, '(A)') '  '


    IF( iThermoDynamics == 1 ) THEN
        WRITE(TECFLG, '(A)         ') '%########Checking energy balance in #########################'
        WRITE(TECFLG, '(A, 2ES17.9) ') 'Total energy in(FA)                 = ',ENEbalt_FA
        WRITE(TECFLG, '(A)') '  '

        WRITE(TECFLG, '(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########Location of Tbk ######################################'
        WRITE(TECFLG, '(A, 1I5.1  )') 'J index  for Tbk      = ', J4Tbk
        WRITE(TECFLG, '(A, 1ES17.9)') 'y  for Tpc            = ', YCC(J4Tbk)
        WRITE(TECFLG, '(A)') ''
        WRITE(TECFLG, '(A)         ') '%########Location of Tpc ######################################'
        WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(K)                = ', T4CpMaX * T0
        WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(C)                = ', T4CpMaX * T0- TK2C
        WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(undim)            = ', T4CpMax
        IF(DMAX1(Twal(1), Twal(2)) >  T4CpMax) THEN
            WRITE(TECFLG, '(A, 1I5.1  )') 'J index  for Tpc      = ', J4Tpc
            WRITE(TECFLG, '(A, 1ES17.9)') 'y  for Tpc            = ', YCC(J4Tpc)
            IF(YCC(J4Tpc) >  0.0_WP) THEN
                WRITE(TECFLG, '(A, 1ES17.9)') 'y+                    = ', (1.0_WP - DABS(YCC(J4Tpc))) * Ret_io(2) / Ldist_io(2)
            ELSE
                WRITE(TECFLG, '(A, 1ES17.9)') 'y+                    = ', (1.0_WP - DABS(YCC(J4Tpc))) * Ret_io(1) / Ldist_io(1)
            END IF
        ELSE
            WRITE(TECFLG, '(A)'      ) 'The current Temperature ange does not cover Tpc.'
        END IF
        WRITE(TECFLG, '(A)') '  '

        WRITE(TECFLG, '(A)         ') '%########sided bulk average(Tbuk seperation)###################'
        WRITE(TECFLG, '(A)          ') '%#            (based on sided bulk T, K and Qw (Kasagi1996))#'
        WRITE(TECFLG, '(A, 2X, I3.1, 2ES17.9)')'Jindex, Ycc, = ', J4Tbk, YCC(J4Tbk)
        WRITE(TECFLG, '(A, 15X, 2(A,7X), A)')   '%#', 'side(Y = -1)', 'side(Y =+ 1)', 'gloabl'
        WRITE(TECFLG, '(A, 3ES17.9)')       'Y+max       = ', (YCC(J4Tbk) + 1.0_WP) * Ret_io(1) / Ldist_io(1), &
        (1.0_WP - YCC(J4Tbk)) * Ret_io(2) / Ldist_io(2)
        WRITE(TECFLG, '(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########based on Int{\rho}dy (ARithmetIC mean DENSITY) #####'
        WRITE(TECFLG, '(A, 1ES17.9) ') 'Tau   = 0.5_WP * (Tauw1 + Tauw2)       = ', Tauw_ave_io
        WRITE(TECFLG, '(A, 1ES17.9) ') 'GravityF = 2.0 * F_A* Int{RHO}dy      = ', DenAvew* 2.0_WP * F_A
        WRITE(TECFLG, '(A, 1ES17.9) ') 'DENSITY = 1 / Ly* Int{RHO}dy         = ', DenAvew
        WRITE(TECFLG, '(A, 1ES17.9) ') 'VIScity = 1 / Ly* Int{Mu}dy          = ', VisAvew
        WRITE(TECFLG, '(A, 1ES17.9) ') 'Utau  = SQRT(Tau_ave/ DenAve)    = ',Utaw_ave_io
        WRITE(TECFLG, '(A, 1ES17.9) ') 'Re_tau  = REN * Utawave* DenAve/ MuAve  = ',Ret_ave_io
    END IF

    CLOSE(TECFLG)

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE PP_SSzero_SIDED(FLG)
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE
    INTEGER(4) :: FLG
    INTEGER(4) :: J
    REAL(WP) :: SSmax, SSFA, Umax, UFA
    REAL(WP) :: MDOtt(2), DHDOtt(2), Ltt(2), Gbukt(2), Hbukt(2), DHbukt(2)
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCP
    REAL(WP) :: spline_interpolation_HB


    !============ Sided by zero stress balancE ===================================
    !=========Find the zero stress balance value.==================================
    SSmax = 1.0E14_WP
    DO J = 1, NCL2
        IF(FLg == 1 .AND. iThermoDynamics == 1) THEN
            SSFA = Tau_Mean_RA(J, 1, 2) - Uff2d_FA(J, 1, 2) +bdfcintg(J)
        ELSE IF (flg == 2 .AND. iThermoDynamics /= 1) THEN
            SSFA = Tau_Mean_RA(J, 1, 2) - Uf2d_RA(J, 1, 2) +bdfcintg(J)
        END IF
        IF (DABS(SSFA) < SSmax) THEN
            SSmax = SSFA
            J4SS0 = J
        END IF

    END DO
    !WRITE(*, *) '# zero ss', J4SS0, Ycc(J4SS0)

    DO J = 1, NCL2
        IF(J < J4SS0) THEN ! neAR Y = -1
            TauwSD(J) = Tauw_io(1)
            DensSD(J) = Dwal(1)
            viscsD(J) = Mwal(1)
            YWdISD(J) = YCC(J) + 1.0_WP
        ELSE                ! neAR Y =+ 1
            TauwSD(J) = Tauw_io(2)
            DensSD(J) = Dwal(2)
            viscsD(J) = Mwal(2)
            YWdISD(J) = 1.0_WP - YCC(J)
        END IF
    END DO

    IF(iThermoDynamics == 1) THEN
        DO J = 1, NCL2
            IF(J < J4SS0) THEN ! neAR Y = -1
                CpSD(J) = Cpwal(1)
                QwSD(J) = Qw(1)
                TwSD(J) = Twal(1)
                HwSD(J) = Hwal_RA(1)
            ELSE                ! neAR Y =+ 1
                CpSD(J) = Cpwal(2)
                QwSD(J) = Qw(2)
                TwSD(J) = Twal(2)
                HwSD(J) = Hwal_RA(2)
            END IF
        END DO

        !===== Sided bulk g and H =============
        MDOtt  = 0.0_WP
        DHDOtt  = 0.0_WP
        UbukSsd = 0.0_WP
        MbukSsd = 0.0_WP
        DbukSsd = 0.0_WP
        Ltt  = 0.0_WP

        DO J = 1, J4SS0
            MDOtt(1) = MDOtt(1) + G1xztL_F0_io(J, 1) / DYFI(J)
            DHDOtt(1) = DHDOtt(1) + DHxztL_F0_io(J) / DYFI(J)
            UbukSsd(1) = UbukSsd(1) + U1xztL_F0_io(J, 1) / DYFI(J)
            MbukSsd(1) = MbukSsd(1) + M1xztL_F0_io(J) / DYFI(J)
            DbukSsd(1) = DbukSsd(1) + D1xztL_F0_io(J) / DYFI(J)
            Ltt(1) = Ltt(1) + 1.0_WP / DYFI(J)
        END DO
        DO J = J4SS0, NCL2
            MDOtt(2) = MDOtt(2) + G1xztL_F0_io(J, 1) / DYFI(J)
            DHDOtt(2) = DHDOtt(2) + DHxztL_F0_io(J) / DYFI(J)
            UbukSsd(2) = UbukSsd(2) + U1xztL_F0_io(J, 1) / DYFI(J)
            MbukSsd(2) = MbukSsd(2) + M1xztL_F0_io(J) / DYFI(J)
            DbukSsd(2) = DbukSsd(2) + D1xztL_F0_io(J) / DYFI(J)
            Ltt(2) = Ltt(2) + 1.0_WP / DYFI(J)
        END DO

        Gbukt(1:2) = MDOtt(1:2) / Ltt(1:2)
        DHbukt(1:2) = DHDOtt(1:2) / Ltt(1:2)
        UbukSsd(1:2) = UbukSsd(1:2) / Ltt(1:2)
        MbukSsd(1:2) = MbukSsd(1:2) / Ltt(1:2)
        DbukSsd(1:2) = DbukSsd(1:2) / Ltt(1:2)

        !WRITE(*, '(A, 4ES13.5)') '# MDOt, HDOt, Gbuk, Hbuk', MDOt, HDOt, Gbuk, Hbuk
        !=============== search === TablE ==========
        DO J = 1,  2

            CALL THERM_PROP_UPDATE_FROM_DH(DHbukt(J), Hbukt(J), TbukSsd(J), DbukSsd(J), MbukSsd(J), &
            KbukSsd(J), CpbkSsd(J), BbukSsd(J))
            !UbukSsd(J) = Gbukt(J) / D_tmp
        END DO

        RebkSsd(1:2) = DbukSsd(1:2) * UbukSsd(1:2) * Ldist_io(1:2) / MbukSsd(1:2) * REN ! based on diameter or the whole hEIGht of Channel.
        PrbkSsd(1:2) = MbukSsd(1:2) * CpbkSsd(1:2) / KbukSsd(1:2) * Prt0

        !WRITE(*, *) '# RebkSsd',RebkSsd(1:2)
        !WRITE(*, *) '# PrbkSsd',PrbkSsd(1:2)

    END IF

    RETURN
END SUBROUTINE
!**********************************************************************************************************************************
SUBROUTINE PP_J4TbukTsd
    USE VARS_AVERAGED_XZ_io

    IMPLICIT NONE
    INTEGER(4) :: J
    REAL(WP) :: TTmax, TTdIF
    INTEGER(4) :: FLG
    REAL(WP) :: SSmax, SSFA, Umax, UFA
    REAL(WP) :: MDOtt(2), DHDOtt(2), Ltt(2), Gbukt(2), Hbukt(2), DHbukt(2)
    REAL(WP) :: spline_interpolation_HT
    REAL(WP) :: spline_interpolation_HK
    REAL(WP) :: spline_interpolation_HCp
    REAL(WP) :: spline_interpolation_HB


    !=========Find the TbulK ==================================
    TTmax = 1.0E14_WP
    DO J = 1, NCL2
        TTdIF = DABS(T1xztL_F0_io(J) - Tbuk)
        IF (TTdIF < TTmax) THEN
            TTmax = TTdIF
            J4Tbk = J
        END IF
    END DO

    L4Tbk(1) = YCC(J4Tbk) - (-1.0_WP)
    L4Tbk(2) = 1.0_WP - YCC(J4Tbk)
    !WRITE(*, *) '# J4Tbk', J, L4Tbk(1), L4Tbk(2)

    MDOtt  = 0.0_WP
    DHDOtt  = 0.0_WP
    UbukTsd = 0.0_WP
    MbukTsd = 0.0_WP
    DbukTsd = 0.0_WP
    Ltt  = 0.0_WP

    DO J = 1, J4Tbk
        MDOtt(1) = MDOtt(1) + G1xztL_F0_io(J, 1) / DYFI(J)
        DHDOtt(1) = DHDOtt(1) + DHxztL_F0_io(J) / DYFI(J)
        UbukTsd(1) = UbukTsd(1) + U1xztL_F0_io(J, 1) / DYFI(J)
        MbukTsd(1) = MbukTsd(1) + M1xztL_F0_io(J) / DYFI(J)
        DbukTsd(1) = DbukTsd(1) + D1xztL_F0_io(J) / DYFI(J)
        Ltt(1) = Ltt(1) + 1.0_WP / DYFI(J)
    END DO
    DO J = J4TbK, NCL2
        MDOtt(2) = MDOtt(2) + G1xztL_F0_io(J, 1) / DYFI(J)
        DHDOtt(2) = DHDOtt(2) + DHxztL_F0_io(J) / DYFI(J)
        UbukTsd(2) = UbukTsd(2) + U1xztL_F0_io(J, 1) / DYFI(J)
        MbukTsd(2) = MbukTsd(2) + M1xztL_F0_io(J) / DYFI(J)
        DbukTsd(2) = DbukTsd(2) + D1xztL_F0_io(J) / DYFI(J)
        Ltt(2) = Ltt(2) + 1.0_WP / DYFI(J)
    END DO

    Gbukt(1:2) = MDOtt(1:2) / Ltt(1:2)
    DHbukt(1:2) = DHDOtt(1:2) / Ltt(1:2)
    UbukTsd(1:2) = UbukTsd(1:2) / Ltt(1:2)
    MbukTsd(1:2) = MbukTsd(1:2) / Ltt(1:2)
    DbukTsd(1:2) = DbukTsd(1:2) / Ltt(1:2)

    !WRITE(*, '(A, 4ES13.5)') '# MDOt, HDOt, Gbuk, Hbuk', MDOt, HDOt, Gbuk, Hbuk
    !=============== search === TablE ==========
    DO J = 1,  2
        CALL THERM_PROP_UPDATE_FROM_DH(DHbukt(J), Hbukt(J), TbukTsd(J), DbukTsd(J), MbukTsd(J), &
        KbukTsd(J), CpbkTsd(J), BbukTsd(J))
    END DO
    !============ Sided by zero stress balancE ===================================
    !=========Find the zero stress balance value.==================================
    TTmax = 1.0E14_WP
    ! find where TbukSsd J
    DO J = 1, J4SS0
        TTdIF = TbukSsd(1) - T1xztL_F0_io(J)
        IF (DABS(TTdIF) < TTmax) THEN
            TTmax = TTdIF
            J4TbukSsd(1) = J
        END IF
    END DO
    TTmax = 1.0E14_WP
    DO J = J4SS0, NCL2
        TTdIF = TbukSsd(2) - T1xztL_F0_io(J)
        IF (DABS(TTdIF) < TTmax) THEN
            TTmax = TTdIF
            J4TbukSsd(2) = J
        END IF
    END DO
    !WRITE(*, *) '# J4TbukSsd', J4TbukSsd(1:2), Ycc(J4TbukSsd(1:2)), TbukSsd(1:2)

    L4TbkSsd(1) = YCC(J4TbukSsd(1)) - (-1.0_WP)
    L4TbkSsd(2) = 1.0_WP - YCC(J4TbukSsd(2))


    RETURN
END SUBROUTINE



!**********************************************************************************************************************************
!===============================================================================
SUBROUTINE PP_Umax_SIDED
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    INTEGER(4) :: J
    REAL(WP) :: SSmax, SSFA, Umax, UFA
    REAL(WP) :: MDOtt(2), HDOtt(2), Ltt(2), Gbukt(2), Hbukt(2), DHbukt(2)
    REAL(WP) :: H_tmp, T_tmp, D_tmp, K_tmp, Ptemp

    !=============== Sd by max. velocitY ========================
    !=========Find the max. VelocitY ==================================
    Umax = 0.0_WP
    DO J = 1, NCL2
        UFA = G1xztL_F0_io(J, 1) / D1xztL_F0_io(J)
        IF (UFA > Umax) THEN
            Umax = UFA
            J4MaxU = J
        END IF
    END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRT_FLOW_Budgets_Profile_XZ_io(STR)
    USE VARS_AVERAGED_XZ_io
    IMPLICIT NONE

    CHARACTER(5), INTENT(IN) :: STR  ! Favre or Reynd
    CHARACTER(128) :: FLNM, FLNN
    CHARACTER(5) :: STDIM(2)
    CHARACTER(4) :: STFASTRESS(8)
    CHARACTER(2) :: SJ2
    REAL(WP) :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
    REAL(WP) :: tke,Ruv_vIS
    INTEGER(4) :: I, J, N, L, TECFLG(2), TECFLG3, NMax
    REAL(WP) ::   COE, tem, M_tmp
    REAL(WP) ::   TKE2, EPPSI, Fmu
    REAL(WP) ::   Cmu, K2De

    REAL(WP) :: FCT(0 : NND2, NDV)
    REAL(WP) :: COE1, COE2, DENtemp, Dintg, intgbdfc

    IF(MYID == 0) CALL CHKRLHDL(' **>WRT_OUT_BUDGET ' // STR// ' at ', MYID, PhyTIME_io)
    !============================budgetS =============================================
    !============================budgetS =============================================
    IF(TRIM(STR) == 'Favre') THEN
        STFASTRESS(1) = 'FAuu'
        STFASTRESS(2) = 'FAuv'
        STFASTRESS(3) = 'FAuw'
        STFASTRESS(4) = 'FAvv'
        STFASTRESS(5) = 'FAvw'
        STFASTRESS(6) = 'FAww'
        STFASTRESS(7) = 'FATK'
        STFASTRESS(8) = 'FAMK'
    ELSE IF (TRIM(STR) == 'Reynd') THEN
        STFASTRESS(1) = 'RAuu'
        STFASTRESS(2) = 'RAuv'
        STFASTRESS(3) = 'RAuw'
        STFASTRESS(4) = 'RAvv'
        STFASTRESS(5) = 'RAvw'
        STFASTRESS(6) = 'RAww'
        STFASTRESS(7) = 'RATK'
        STFASTRESS(8) = 'RAMK'
    ELSE
    END IF

    NMax = 1
    STDIM(1) = 'undim'
    TECFLG(1) = 51
    IF(iPPDimension == 1) THEN
        NMax = 2
        STDIM(2) = 'dimen'
        TECFLG(2) = 102
    END IF

    DO L = 1, 8

        DO N = 1, NMax
            WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
            FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(N)) // '.Profile.BUDGETS.' &
            // STR// '.' // TRIM(STFASTRESS(L)) //  &
            '.' // TRIM(PNTIM) // '.plt'
            OPEN (TECFLG(N), FILE = TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG(N), '(A)') 'TITLE = "' // STR// ' Averged Flow (28 variables)" '
            J = 0;                          WRITE(TECFLG(N), '(A)', AdvancE = "no") 'variables = '

            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Tauw", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Dw", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'DInt", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'D", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Muw", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'MuInt", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'M", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'F_A", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'YCC", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Ywd", '

            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_stres", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_dissp", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'pduDX_stran", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'Turbu_diffu", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'dpuDX_diffu", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_diffu", '

            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'press_accel1", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'viscs_accel1", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_Dgbdf1", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_Drvfc1", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'balance1"'

            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'turss_accl2", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_Dgbdf2", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'prodc_Drvfc2", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'balance2"'

            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'pressure3", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)', AdvancE = "no") '"' // SJ2 // 'vistress3", '
            J = J + 1; WRITE(SJ2, '(1I2.2)') J; WRITE(TECFLG(N), '(A)'          ) '"' // SJ2 // 'balance3"'

            WRITE(TECFLG(N), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '
        END DO

        DO J = 1, NCL2
            IF (L ==7) THEN
                WRITE(TECFLG(1), '(28ES20.12)') &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                viscsD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdISD(J), &
                Budg_prodc_stres_tke(J), Budg_viscs_dissp_tke(J), Budg_pduDX_stran_tke(J), &
                Budg_Turbu_diffu_tke(J), Budg_DpuDX_diffu_tke(J), Budg_viscs_diffu_tke(J), &
                Budg_press_accl1_tke(J), Budg_viscs_accl1_tke(J), 0.0_WP, Budg_prodc_Dvfc1_tke(J), Budg_Balance1_tke(J), &
                Budg_turss_accl2_tke(J), Budg_prodc_gvfc2_tke(J),         Budg_prodc_Dvfc1_tke(J), Budg_Balance2_tke(J), &
                Budg_pressure3_tke(J),   Budg_vistress3_tke(J),   Budg_Balance3_tke(J)

            ELSE IF (L == 8) THEN
                WRITE(TECFLG(1), '(28ES20.12)') &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                viscsD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdISD(J), &
                Budg_prodc_stres_mke(J), Budg_viscs_dissp_mke(J), Budg_pduDX_stran_mke(J), &
                Budg_Turbu_diffu_mke(J), Budg_DpuDX_diffu_mke(J), Budg_viscs_diffu_mke(J), &
                Budg_press_accl1_mke(J), Budg_viscs_accl1_mke(J), 0.0_WP, Budg_prodc_Dvfc1_mke(J), Budg_Balance1_mke(J), &
                Budg_turss_accl2_mke(J), Budg_prodc_gvfc2_mke(J),         Budg_prodc_Dvfc1_mke(J), Budg_Balance2_mke(J), &
                Budg_pressure3_mke(J),   Budg_vistress3_mke(J),   Budg_Balance3_mke(J)
            ELSE
                WRITE(TECFLG(1), '(28ES20.12)') &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                viscsD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdISD(J), &
                Budg_prodc_stres_Duiuj(J, L), Budg_viscs_dissp_Duiuj(J, L), Budg_pduDX_stran_Duiuj(J, L), &
                Budg_Turbu_diffu_Duiuj(J, L), Budg_DpuDX_diffu_Duiuj(J, L), Budg_viscs_diffu_Duiuj(J, L), &
                Budg_press_accl1_Duiuj(J, L), Budg_viscs_accl1_Duiuj(J, L), 0.0_WP, &
                Budg_prodc_Dvfc1_Duiuj(J, L), Budg_Balance1_Duiuj(J, L), &
                Budg_turss_accl2_Duiuj(J, L), Budg_prodc_gvfc2_Duiuj(J, L),         &
                Budg_prodc_Dvfc1_Duiuj(J, L), Budg_Balance2_Duiuj(J, L), &
                Budg_pressure3_Duiuj(J, L),   Budg_vistress3_Duiuj(J, L),   Budg_Balance3_Duiuj(J, L)
            END IF
            !                IF(iPPDimension == 1) THEN
            !                    WRITE(TECFLG(2), '(13ES20.12)') YCC(J) * L0, (1.0_WP - DABS(YCC(J))) * REN * COE, COE* U0, &
            !                        Budg_prodc_stres_Duiuj(J, L) * SCaling1, Budg_Turbu_diffu_Duiuj(J, L) * SCaling1, &
            !                        Budg_DpuDX_diffu_Duiuj(J, L) * SCaling1, Budg_pduDX_stran_Duiuj(J, L) * SCaling1, &
            !                        Budg_press_accl1_Duiuj(J, L) * SCaling1, Budg_viscs_accl1_Duiuj(J, L) * SCaling2, &
            !                        Budg_viscs_diffu_Duiuj(J, L) * SCaling2, Budg_viscs_dissp_Duiuj(J, L) * SCaling2, &
            !                        Budg_prodc_gvfc2_Duiuj(J, L) * D0 * U0, Budg_Balance1_Duiuj(J, L) * SCaling2
            !                END IF
        END DO

        DO N = 1, NMax
            CLOSE(TECFLG(N))
        END DO
    END DO

    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRITE_SPECO_AVE_PROFILE(STR)
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
    USE mesh_info
    USE init_info
    USE SPECO_info
    USE postprocess_info
    USE WRT_INFO

    IMPLICIT NONE
    CHARACTER(4), INTENT(IN) :: STR

    CHARACTER(15) :: PNTIM
    CHARACTER(15) :: PNLOC
    INTEGER(4) :: DFLG(MGRID)
    CHARACTER(256) :: FLNAME
    !REAL(WP) :: Ret_ave, U_tau_ave
    REAL(WP) :: AKE, WDS(2)

    INTEGER(4) :: N, JJ, L, KC, IC, M

    IF(MYID /= 0) RETURN
    IF(iPPSpectra /= 1) RETURN

    IF(TRIM(STR) == 'FLOW') M = 1
    IF(TRIM(STR) == 'HEAT') M = 2

    N3MH = NCL3 / 2 + 1
    N3MD= NCL3 + 2
    N1MH = NCL1_io / 2 + 1
    N1MD= NCL1_io + 2
    !CALL PP_wall_thermal_shear(flgxzt)

    !========== Test for a ASCII outpuT ==============
    !Ret_ave  = 0.5_WP * (Ret_io(1) + Ret_io(2)) !;!WRITE(*, *) Ret_ave, Ret_io(1), Ret_io(2) !test
    !U_tau_avE = 0.5_WP * (Utaw_io(1) + Utaw_io(2))
    OPEN(100, FILE = TRIM(FilePath4) // 'CHK_PROBE_for_spectra_instant.dat')
    WRITE(100, '(A)') '## MGRID, JJ, YCC, Yplus1, Yplus2, Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2)'
    DO N = 1, MGRID
        JJ = JGMOV(N)
        DFLG(N) = 100 + N

        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
        WRITE(PNLOC, '(1I3.3)')   N

        !==============correlationS ======================================
        DO L = 1, 2
            IF(L == 1) THEN
                FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.2PCorrelation.X.T'  &
                // TRIM(PNTIM) // '.YLC' // TRIM(PNLOC) // '.plt'
                OPEN(DFLG(N), FILE = TRIM(ADJUSTL(FLNAME)))
                WRITE(DFLG(N), '(A)') 'TITLE = " Correlation (Streamwise sepARation)" '
                WRITE(DFLG(N), '(A)', AdvancE = "no") 'variables = "1X", "U11", "U22", "U33", "U12", "U13", "U23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"O11", "O22", "O33", "O12", "O13", "O23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"VO11", "VO12", "VO13", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"VO21", "VO22", "VO23", '
                WRITE(DFLG(N), '(A)'         ) '"VO31", "VO32", "VO33"  '

                WRITE(DFLG(N), '(A, 1ES13.5, A, 2ES13.5, A, 3ES13.5, A)') &
                'ZONE T = "Rkk(X) at Y / Delta= ', YCC(JJ), ' Utauw12 = ',Utaw_io(1:2), ' Ret12 = ',Ret_io(1:2), Ret_ave_io, ' " '

                WDS(1) = YCC(JJ) - (-1.0_WP)
                WDS(2) = 1.0_WP - YCC(JJ)
                WRITE(100, '(2I4.2, 10ES13.5)') N, JJ, YCC(JJ), WDS(1:2) * Ret_io(1:2) / Ldist_io(1:2),  &
                Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2)

                DO IC = 1, N1MH
                    WRITE(DFLG(N), '(22ES13.5)') 0.5_WP * ( XND_io(IC) + XND_io(IC + 1) ), &
                    R11X1_xztLa (JJ, IC, M), R22X1_xztLa (JJ, IC, M), R33X1_xztLa (JJ, IC, M), &
                    R12X1_xztLa (JJ, IC, M), R13X1_xztLa (JJ, IC, M), R23X1_xztLa (JJ, IC, M), &
                    V11X1_xztLa (JJ, IC, M), V22X1_xztLa (JJ, IC, M), V33X1_xztLa (JJ, IC, M), &
                    V12X1_xztLa (JJ, IC, M), V13X1_xztLa (JJ, IC, M), V23X1_xztLa (JJ, IC, M), &
                    VO11X1_xztLa(JJ, IC, M), VO12X1_xztLa(JJ, IC, M), VO13X1_xztLa(JJ, IC, M), &
                    VO21X1_xztLa(JJ, IC, M), VO22X1_xztLa(JJ, IC, M), VO23X1_xztLa(JJ, IC, M), &
                    VO31X1_xztLa(JJ, IC, M), VO32X1_xztLa(JJ, IC, M), VO33X1_xztLa(JJ, IC, M)
                END DO
                CLOSE(DFLG(N))
            END IF

            IF(L == 2) THEN
                FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.2PCorrelation.Z.T'  &
                // TRIM(PNTIM) // '.YLC' // TRIM(PNLOC) // '.plt'
                OPEN(DFLG(N), FILE = TRIM(ADJUSTL(FLNAME)))
                WRITE(DFLG(N), '(A)') 'TITLE = " Correlation (SpanWISe sepARation)" '
                WRITE(DFLG(N), '(A)', AdvancE = "no") 'variables = "3Z", "U11", "U22", "U33", "U12", "U13", "U23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"O11", "O22", "O33", "O12", "O13", "O23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"VO11", "VO12", "VO13", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"VO21", "VO22", "VO23", '
                WRITE(DFLG(N), '(A)'         ) '"VO31", "VO32", "VO33"  '
                WRITE(DFLG(N), '(A, 1ES13.5, A, 2ES13.5, A, 3ES13.5, A)') &
                'ZONE T = "Rkk(Z) at Y / Delta= ', YCC(JJ), ' Utauw12 = ',Utaw_io(1:2), ' Ret12 = ',Ret_io(1:2), Ret_ave_io, ' " '

                DO KC = 1, N3MH
                    WRITE(DFLG(N), '(22ES13.5)') 0.5_WP * ( ZND(KC) + ZND(KC + 1) ), &
                    R11X3_xztLa (JJ, KC, M), R22X3_xztLa (JJ, KC, M), R33X3_xztLa (JJ, KC, M), &
                    R12X3_xztLa (JJ, KC, M), R13X3_xztLa (JJ, KC, M), R23X3_xztLa (JJ, KC, M), &
                    V11X3_xztLa (JJ, KC, M), V22X3_xztLa (JJ, KC, M), V33X3_xztLa (JJ, KC, M), &
                    V12X3_xztLa (JJ, KC, M), V13X3_xztLa (JJ, KC, M), V23X3_xztLa (JJ, KC, M), &
                    VO11X3_xztLa(JJ, KC, M), VO12X3_xztLa(JJ, KC, M), VO13X3_xztLa(JJ, KC, M), &
                    VO21X3_xztLa(JJ, KC, M), VO22X3_xztLa(JJ, KC, M), VO23X3_xztLa(JJ, KC, M), &
                    VO31X3_xztLa(JJ, KC, M), VO32X3_xztLa(JJ, KC, M), VO33X3_xztLa(JJ, KC, M)
                END DO
                CLOSE(DFLG(N))

            END IF
        END DO

        !============== ENErgy espectra ======================================
        DO L = 1, 2
            IF(L == 1) THEN
                FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.energy.Xwavenumber.T'  &
                // TRIM(PNTIM) // '.YLC' // TRIM(PNLOC) // '.plt'
                OPEN(DFLG(N), FILE = TRIM(ADJUSTL(FLNAME)))
                WRITE(DFLG(N), '(A)') 'TITLE = " energy.spectra (Streamwise)" '
                WRITE(DFLG(N), '(A)', AdvancE = "no") 'variables = "1IC", "2WaveNumberX"'
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"E11", "E22", "E33", "E12", "E13", "E23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"H11", "H22", "H33", "H12", "H13", "H23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"EH11", "EH12", "EH13", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"EH21", "EH22", "EH23", '
                WRITE(DFLG(N), '(A)'         ) '"EH31", "EH32", "EH33"  '

                WRITE(DFLG(N), '(A, 1ES13.5, A, 2ES13.5, A, 3ES13.5, A)') &
                'ZONE T = "Ekk(K1) at Y / Delta= ', YCC(JJ), ' Utauw12 = ',Utaw_io(1:2), ' Ret12 = ',Ret_io(1:2), Ret_ave_io, ' " '
                DO IC = 1, N1MH
                    AKE = ( DBLE(IC- 1) * 2.0_WP * PI /HX_io )
                    WRITE(DFLG(N), '(1I8.1, 22ES13.5)') IC, AKE, &
                    ENE11T_xztLa(JJ, IC, M), ENE22T_xztLa(JJ, IC, M), ENE33T_xztLa(JJ, IC, M), &
                    ENE12T_xztLa(JJ, IC, M), ENE13T_xztLa(JJ, IC, M), ENE23T_xztLa(JJ, IC, M), &
                    ENV11T_xztLa(JJ, IC, M), ENV22T_xztLa(JJ, IC, M), ENV33T_xztLa(JJ, IC, M), &
                    ENV12T_xztLa(JJ, IC, M), ENV13T_xztLa(JJ, IC, M), ENV23T_xztLa(JJ, IC, M), &
                    EVO11T_xztLa(JJ, IC, M), EVO12T_xztLa(JJ, IC, M), EVO13T_xztLa(JJ, IC, M), &
                    EVO21T_xztLa(JJ, IC, M), EVO22T_xztLa(JJ, IC, M), EVO23T_xztLa(JJ, IC, M), &
                    EVO31T_xztLa(JJ, IC, M), EVO32T_xztLa(JJ, IC, M), EVO33T_xztLa(JJ, IC, M)
                END DO
                CLOSE(DFLG(N))
            END IF

            IF(L == 2) THEN
                FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.energy.Zwavenumber.T'  &
                // TRIM(PNTIM) // '.YLC' // TRIM(PNLOC) // '.plt'
                OPEN(DFLG(N), FILE = TRIM(ADJUSTL(FLNAME)))
                WRITE(DFLG(N), '(A)') 'TITLE = " energy.spectra (spanWISe)" '
                WRITE(DFLG(N), '(A)', AdvancE = "no") 'variables = "3KC", "2WaveNumberZ"'
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"E11", "E22", "E33", "E12", "E13", "E23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"H11", "H22", "H33", "H12", "H13", "H23", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"EH11", "EH12", "EH13", '
                WRITE(DFLG(N), '(A)', AdvancE = "no") '"EH21", "EH22", "EH23", '
                WRITE(DFLG(N), '(A)'         ) '"EH31", "EH32", "EH33"  '
                WRITE(DFLG(N), '(A, 1ES13.5, A, 2ES13.5, A, 3ES13.5, A)') &
                'ZONE T = "Ekk(K3) at Y / Delta= ', YCC(JJ), ' Utauw12 = ',Utaw_io(1:2), ' Ret12 = ',Ret_io(1:2), Ret_ave_io, ' " '
                DO KC = 1, N3MH
                    AKE = ( DBLE(KC - 1) * 2.0_WP * PI /HZ )
                    WRITE(DFLG(N), '(1I8.1, 22ES13.5)') KC, AKE, &
                    ENE11Z_xztLa(JJ, KC, M), ENE22Z_xztLa(JJ, KC, M), ENE33Z_xztLa(JJ, KC, M), &
                    ENE12Z_xztLa(JJ, KC, M), ENE13Z_xztLa(JJ, KC, M), ENE23Z_xztLa(JJ, KC, M), &
                    ENV11Z_xztLa(JJ, KC, M), ENV22Z_xztLa(JJ, KC, M), ENV33Z_xztLa(JJ, KC, M), &
                    ENV12Z_xztLa(JJ, KC, M), ENV13Z_xztLa(JJ, KC, M), ENV23Z_xztLa(JJ, KC, M), &
                    EVO11Z_xztLa(JJ, KC, M), EVO12Z_xztLa(JJ, KC, M), EVO13Z_xztLa(JJ, KC, M), &
                    EVO21Z_xztLa(JJ, KC, M), EVO22Z_xztLa(JJ, KC, M), EVO23Z_xztLa(JJ, KC, M), &
                    EVO31Z_xztLa(JJ, KC, M), EVO32Z_xztLa(JJ, KC, M), EVO33Z_xztLa(JJ, KC, M)
                END DO
                CLOSE(DFLG(N))
            END IF
        END DO

    END DO


    RETURN
END SUBROUTINE

!**********************************************************************************************************************************
SUBROUTINE WRITE_SPECO_AVE_Contour(STR)
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
    USE mesh_info
    USE init_info
    USE SPECO_info
    USE postprocess_info
    USE WRT_INFO

    IMPLICIT NONE
    CHARACTER(4), INTENT(IN) :: STR

    CHARACTER(15) :: PNTIM
    CHARACTER(15) :: PNLOC
    INTEGER(4) :: Dflg = 51
    CHARACTER(256) :: FLNAME
    !REAL(WP) :: Ret_ave, U_tau_ave
    REAL(WP) :: AKE
    REAL(WP) :: yplus

    INTEGER(4) :: N, JJ, L, KC, IC, JJM, JJC, M


    IF(MYID /= 0) RETURN

    IF(TRIM(STR) == 'FLOW') M = 1
    IF(TRIM(STR) == 'HEAT') M = 2

    N3MH = NCL3 / 2 + 1
    N3MD= NCL3 + 2
    N1MH = NCL1_io / 2 + 1
    N1MD= NCL1_io + 2

    !========== Test for a ASCII outpuT ==============
    !Ret_ave  = 0.5_WP * (Ret_io(1) + Ret_io(2)) !;!WRITE(*, *) Ret_ave, Ret_io(1), Ret_io(2) !test
    !U_tau_avE = 0.5_WP * (Utaw_io(1) + Utaw_io(2))

    !===============plane y-Z =====================
    WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
    FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.Contours.yz.' // TRIM(PNTIM) // '.plt'
    OPEN (DFLG, FILE = TRIM(ADJUSTL(FLNAME)))

    WRITE(DFLG, '(A)') 'TITLE = "DNS FLOW YZ -plane"'
    WRITE(DFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "Y+", "WaveNo3", ' !5
    WRITE(DFLG, '(A)', AdvancE = "no") '"U11", "U22", "U33", "U12", "U13", "U23", '!6
    WRITE(DFLG, '(A)', AdvancE = "no") '"O11", "O22", "O33", "O12", "O13", "O23", '!6
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO11", "VO12", "VO13", ' !3
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO21", "VO22", "VO23", ' !3
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO31", "VO32", "VO33"  ' !3
    WRITE(DFLG, '(A)', AdvancE = "no") '"E11", "E22", "E33", "E12", "E13", "E23", ' !6
    WRITE(DFLG, '(A)', AdvancE = "no") '"H11", "H22", "H33", "H12", "H13", "H23", ' !6
    WRITE(DFLG, '(A)', AdvancE = "no") '"EH11", "EH12", "EH13", ' !3
    WRITE(DFLG, '(A)', AdvancE = "no") '"EH21", "EH22", "EH23", ' !3
    WRITE(DFLG, '(A)'         ) '"EH31", "EH32", "EH33"  ' !3

    WRITE(DFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG, PhyTIME, &
    ' ", I = ', 1, ', J = ', NND2, ', K = ', N3MH, ', F =POINT'


    DO KC = 1, N3MH
        AKE = ( DBLE(KC - 1) * 2.0_WP * PI /HZ ) / Ret_ave_io
        DO JJ = 1, NND2
            JJM = JGMV(JJ)
            JJC = JJ
            ypluS = (1.0_WP - DABS(YND(JJ))) * REN * Utaw_ave_io

            IF(JJ == 1)    JJM = 1
            IF(JJ == NND2) JJC = NCL2
            WRITE(DFLG, '(47ES15.7)') 0.0_WP, YND(JJ),ZND(KC), yplus, AKE, &
            0.5_WP * ( R11X3_xztLa(JJC, KC, M) + R11X3_xztLa (JJM, KC, M) ), &
            0.5_WP * ( R22X3_xztLa(JJC, KC, M) + R22X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( R33X3_xztLa(JJC, KC, M) + R33X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( R12X3_xztLa(JJC, KC, M) + R12X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( R13X3_xztLa(JJC, KC, M) + R13X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( R23X3_xztLa(JJC, KC, M) + R23X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V11X3_xztLa(JJC, KC, M) + V11X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V22X3_xztLa(JJC, KC, M) + V22X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V33X3_xztLa(JJC, KC, M) + V33X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V12X3_xztLa(JJC, KC, M) + V12X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V13X3_xztLa(JJC, KC, M) + V13X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( V23X3_xztLa(JJC, KC, M) + V23X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO11X3_xztLa(JJC, KC, M) + VO11X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO12X3_xztLa(JJC, KC, M) + VO12X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO13X3_xztLa(JJC, KC, M) + VO13X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO21X3_xztLa(JJC, KC, M) + VO21X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO22X3_xztLa(JJC, KC, M) + VO22X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO23X3_xztLa(JJC, KC, M) + VO23X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO31X3_xztLa(JJC, KC, M) + VO31X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO32X3_xztLa(JJC, KC, M) + VO32X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( VO33X3_xztLa(JJC, KC, M) + VO33X3_xztLa(JJM, KC, M) ), &
            0.5_WP * ( ENE11Z_xztLa(JJC, KC, M) + ENE11Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE22Z_xztLa(JJC, KC, M) + ENE22Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE33Z_xztLa(JJC, KC, M) + ENE33Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE12Z_xztLa(JJC, KC, M) + ENE12Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE13Z_xztLa(JJC, KC, M) + ENE13Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE23Z_xztLa(JJC, KC, M) + ENE23Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV11Z_xztLa(JJC, KC, M) + ENV11Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV22Z_xztLa(JJC, KC, M) + ENV22Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV33Z_xztLa(JJC, KC, M) + ENV33Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV12Z_xztLa(JJC, KC, M) + ENV12Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV13Z_xztLa(JJC, KC, M) + ENV13Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV23Z_xztLa(JJC, KC, M) + ENV23Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO11Z_xztLa(JJC, KC, M) + EVO11Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO12Z_xztLa(JJC, KC, M) + EVO12Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO13Z_xztLa(JJC, KC, M) + EVO13Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO21Z_xztLa(JJC, KC, M) + EVO21Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO22Z_xztLa(JJC, KC, M) + EVO22Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO23Z_xztLa(JJC, KC, M) + EVO23Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO31Z_xztLa(JJC, KC, M) + EVO31Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO32Z_xztLa(JJC, KC, M) + EVO32Z_xztLa(JJM, KC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO33Z_xztLa(JJC, KC, M) + EVO33Z_xztLa(JJM, KC, M) ) * Ret_ave_io

        END DO
    END DO
    CLOSE(DFLG)
    !===============plane X - Y =====================
    FLNAME = TRIM(FilePath4) // 'Result.IO.Spectral.' // TRIM(STR) // '.Contours.yx.' // TRIM(PNTIM) // '.plt'
    OPEN (DFLG, FILE = TRIM(ADJUSTL(FLNAME)))

    WRITE(DFLG, '(A)') 'TITLE = "DNS FLOW XY-plane"'
    WRITE(DFLG, '(A)', AdvancE = "no") 'variables = "X", "Y", "Z", "Y+", "WaveNo3", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"U11", "U22", "U33", "U12", "U13", "U23", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"O11", "O22", "O33", "O12", "O13", "O23", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO11", "VO12", "VO13", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO21", "VO22", "VO23", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"VO31", "VO32", "VO33"  '
    WRITE(DFLG, '(A)', AdvancE = "no") '"E11", "E22", "E33", "E12", "E13", "E23", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"H11", "H22", "H33", "H12", "H13", "H23", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"EH11", "EH12", "EH13", '
    WRITE(DFLG, '(A)', AdvancE = "no") '"EH21", "EH22", "EH23", '
    WRITE(DFLG, '(A)'         ) '"EH31", "EH32", "EH33"  '

    WRITE(DFLG, '(A, 1I11.1, 1ES13.5, A, 1I4.1, A, 1I4.1, A, 1I4.1, A)') 'ZONE T = " ', ITERG,PhyTIME, &
    ' ", I = ', N1MH, ', J = ', NND2, ', K = ', 1, ', F =POINT'
    DO JJ = 1, NND2
        ypluS = (1.0_WP - DABS(YND(JJ))) * REN * Utaw_ave_io
        JJM = JGMV(JJ)
        JJC = JJ
        IF(JJ == 1)    JJM = 1
        IF(JJ == NND2) JJC = NCL2
        DO IC = 1, N1MH
            AKE = ( DBLE(IC- 1) * 2.0_WP * PI /HX_io ) / Ret_ave_io

            WRITE(DFLG, '(47ES15.7)') XND_io(IC), YND(JJ), 0.0_WP, yplus, AKE, &
            0.5_WP * ( R11X1_xztLa(JJC, IC, M) + R11X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( R22X1_xztLa(JJC, IC, M) + R22X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( R33X1_xztLa(JJC, IC, M) + R33X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( R12X1_xztLa(JJC, IC, M) + R12X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( R13X1_xztLa(JJC, IC, M) + R13X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( R23X1_xztLa(JJC, IC, M) + R23X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V11X1_xztLa(JJC, IC, M) + V11X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V22X1_xztLa(JJC, IC, M) + V22X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V33X1_xztLa(JJC, IC, M) + V33X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V12X1_xztLa(JJC, IC, M) + V12X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V13X1_xztLa(JJC, IC, M) + V13X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( V23X1_xztLa(JJC, IC, M) + V23X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO11X1_xztLa(JJC, IC, M) + VO11X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO12X1_xztLa(JJC, IC, M) + VO12X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO13X1_xztLa(JJC, IC, M) + VO13X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO21X1_xztLa(JJC, IC, M) + VO21X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO22X1_xztLa(JJC, IC, M) + VO22X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO23X1_xztLa(JJC, IC, M) + VO23X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO31X1_xztLa(JJC, IC, M) + VO31X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO32X1_xztLa(JJC, IC, M) + VO32X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( VO33X1_xztLa(JJC, IC, M) + VO33X1_xztLa(JJM, IC, M) ), &
            0.5_WP * ( ENE11T_xztLa(JJC, IC, M) + ENE11T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE22T_xztLa(JJC, IC, M) + ENE22T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE33T_xztLa(JJC, IC, M) + ENE33T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE12T_xztLa(JJC, IC, M) + ENE12T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE13T_xztLa(JJC, IC, M) + ENE13T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENE23T_xztLa(JJC, IC, M) + ENE23T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV11T_xztLa(JJC, IC, M) + ENV11T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV22T_xztLa(JJC, IC, M) + ENV22T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV33T_xztLa(JJC, IC, M) + ENV33T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV12T_xztLa(JJC, IC, M) + ENV12T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV13T_xztLa(JJC, IC, M) + ENV13T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( ENV23T_xztLa(JJC, IC, M) + ENV23T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO11T_xztLa(JJC, IC, M) + EVO11T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO12T_xztLa(JJC, IC, M) + EVO12T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO13T_xztLa(JJC, IC, M) + EVO13T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO21T_xztLa(JJC, IC, M) + EVO21T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO22T_xztLa(JJC, IC, M) + EVO22T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO23T_xztLa(JJC, IC, M) + EVO23T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO31T_xztLa(JJC, IC, M) + EVO31T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO32T_xztLa(JJC, IC, M) + EVO32T_xztLa(JJM, IC, M) ) * Ret_ave_io, &
            0.5_WP * ( EVO33T_xztLa(JJC, IC, M) + EVO33T_xztLa(JJM, IC, M) ) * Ret_ave_io
        END DO
    END DO

    CLOSE(DFLG)


    RETURN
END SUBROUTINE
!======================checking...==========================================================
!        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
!        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.Checking.Favre.' // TRIM(PNTIM) // '.plt'
!        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (22 variables)" '
!        WRITE(TECFLG_FavAG(1), '(7A)') &
!            'variables = "1Y", "2Y+", "3Utau", "4Dwal", ', &
!            '"5viscstress_Tau_Mean_uu", "6viscstress_Tau_Umea_uu", "7viscstress_Tau_Uper_uu"', &
!            '"8viscstress_Tau_Mean_uv",  "9viscstress_Tau_Umea_uv", "10viscstress_Tau_Uper_uv"', &
!            '"11viscstress_Tau_Mean_uw", "12viscstress_Tau_Umea_uw", "13viscstress_Tau_Uper_uw',  &
!            '"14viscstress_Tau_Mean_vv", "15viscstress_Tau_Umea_vv", "16viscstress_Tau_Uper_vv"', &
!            '"17viscstress_Tau_Mean_vw", "18viscstress_Tau_Umea_vw", "19viscstress_Tau_Uper_vw"', &
!            '"20viscstress_Tau_Mean_ww", "21viscstress_Tau_Umea_ww", "22viscstress_Tau_Uper_ww"'
!        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

!        DO J = 1, NCL2
!            WRITE(TECFLG_FavAG(1), '(22ES20.12)') YCC(J), YpluSD(J), UtauSD(J), DensSD(J) &
!                              Tau_Mean_RA(J, 1, 1), viscstress_Tau_Umea(J, 1, 1), viscstress_Tau_Uper(J, 1, 1), &
!                              Tau_Mean_RA(J, 1, 2), viscstress_Tau_Umea(J, 1, 2), viscstress_Tau_Uper(J, 1, 2), &
!                              Tau_Mean_RA(J, 1, 3), viscstress_Tau_Umea(J, 1, 3), viscstress_Tau_Uper(J, 1, 3), &
!                              Tau_Mean_RA(J, 2, 2), viscstress_Tau_Umea(J, 2, 2), viscstress_Tau_Uper(J, 2, 2), &
!                              Tau_Mean_RA(J, 2, 3), viscstress_Tau_Umea(J, 2, 3), viscstress_Tau_Uper(J, 2, 3), &
!                              Tau_Mean_RA(J, 3, 3), viscstress_Tau_Umea(J, 3, 3), viscstress_Tau_Uper(J, 3, 3)

!        END DO
!        CLOSE(TECFLG_FavAG(1))


!        !======================checking== Momentum Equation IS xyz direction ======================================================
!        BuoyForceTT = 0.0_WP
!        FCT = 0.0_WP
!        DO J = 1, NCL2
!            FCT(J, 1) = - D1xztL_F0_io(J) * U_FA(J, 1) * U_FA(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew + &
!                        Tau_Mean_RA(J, 1, 2)      * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew- &
!                        uff2d_FA  (J, 1, 2)      / Utaw_ave_iO / Utaw_ave_iO / DenAvew

!            FCT(J, 2) = - D1xztL_F0_io(J) * U_FA(J, 2) * U_FA(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew + &
!                        Tau_Mean_RA(J, 2, 2)     * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew- &
!                        uff2d_FA  (J, 2, 2)     / Utaw_ave_iO / Utaw_ave_iO / DenAvew-&
!                        dPDX_RA(J, 2)                      / Utaw_ave_iO / Utaw_ave_iO / D1xztL_F0_io(J)!

!            FCT(J, 3) = - D1xztL_F0_io(J) * U_FA(J, 3) * U_FA(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew + &
!                        Tau_Mean_RA(J, 2, 3)     * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew- &
!                        uff2d_FA  (J, 2, 3)      / Utaw_ave_iO / Utaw_ave_iO / DenAvew

!            BuoyForceTT = BuoyForceTT + F_A* D1xztL_F0_io(J) / DYFI(J)
!        END DO
!        !=== At wall surfaceS ====
!        FCT(0, 1) = Tauw_io(1) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew
!        FCT(0, 2) = 0.0_WP
!        FCT(0, 3) = 0.0_WP

!        FCT(NND2, 1) = Tauw_io(2) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew
!        FCT(NND2, 2) = 0.0_WP
!        FCT(NND2, 3) = 0.0_WP

!        !===gradient at cell centrE ====
!        NSFbal_FA = 0.0_WP
!        DO J = 2, NCL2 - 1
!            DO N = 1, NDV
!                NSFbal_FA(J, N) = ( ( YCL2ND_WFB(J + 1) * FCT(J, N) + YCL2ND_WFF(J + 1) * FCT(J + 1, N) ) -             &
!                                  ( YCL2ND_WFF(J) * FCT(J, N) + YCL2ND_WFB(J) * FCT(J - 1, N) ) ) * DYFI(J) + &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * D1xztL_F0_io(J) / DenAvew* DBLE(IBuoF(N)) - &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!!                !WRITE(*, *) N, J, ( ( YCL2ND_WFB(J + 1) * FCT(J, N) + YCL2ND_WFF(J + 1) * FCT(J + 1, N) ) -             &
!!                                  ( YCL2ND_WFF(J) * FCT(J, N) + YCL2ND_WFB(J) * FCT(J - 1, N) ) ) * DYFI(J), &
!!                                  F_A / Utaw_ave_iO / Utaw_ave_io * D1xztL_F0_io(J) / DenAvew* DBLE(IBuoF(N)), &
!!                                  F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!            END DO
!        END DO
!        DO N = 1, NDV
!            NSFbal_FA(1,   N) = ( ( YCL2ND_WFB(1 + 1) * FCT(1, N) +                                           &
!                                   YCL2ND_WFF(1 + 1) * FCT(1 + 1, N) ) - FCT(0, N) ) * DYFI(1) +               &
!                                 F_A / Utaw_ave_iO / Utaw_ave_io * D1xztL_F0_io(1) / DenAvew* DBLE(IBuoF(N)) - &
!                                 F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))

!            NSFbal_FA(NCL2, N) = ( FCT(NND2, N) - ( YCL2ND_WFF(NCL2) * FCT(NCL2, N) +                       &
!                                                 YCL2ND_WFB(NCL2) * FCT(NCL2 - 1, N) ) ) * DYFI(NCL2) +      &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * D1xztL_F0_io(NCL2) / DenAvew* DBLE(IBuoF(N)) - &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!        END DO

!        !=== INTegral along y direction ===============================
!        NSFbalt_FA = 0.0_WP
!        DO J = 1, NCL2
!            DO N = 1, NDV
!                NSFbalt_FA(N) = NSFbalt_FA(N) + NSFbal_FA(J, N) / DYFI(J)
!            END DO
!            !!WRITE(*, *) J, NSFbalt_FA(1:3)
!        END DO


!        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
!        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.Checking.NSForceBalance.FA.' // TRIM(PNTIM) // '.plt'
!        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (21 variables)" '
!        WRITE(TECFLG_FavAG(1), '(6A)') &
!            'variables = "1Y", "2Y+", "3Utau", "4FC_x", "5FC_y", "6FC_z","7FCGRD_x", "8FCGRD_y", "9FCGRD_z"'
!        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

!        Dintg = 0.0_WP
!        DO J = 1, NCL2

!            !IF(YND(J) < 0.0_WP) THEN
!            IF(J < J4SS0) THEN
!                COE = Utaw_io(1)
!            ELSE
!                COE = Utaw_io(2)
!            END IF
!            IF(J == 1) THEN
!                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
!                                  Dwal(1) )
!            ELSE IF(J == NCL2) THEN
!                DENtemp = 0.5_WP * ( Dwal(2) +             &
!                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
!            ELSE
!                DENtemp = 0.5_WP * ( ( YCL2ND_WFB(J + 1) * D1xztL_F0_io(J) + YCL2ND_WFF(J + 1) * D1xztL_F0_io(J + 1) ) +             &
!                                  ( YCL2ND_WFF(J) * D1xztL_F0_io(J) + YCL2ND_WFB(J) * D1xztL_F0_io(J - 1) ) )
!            END IF
!            Dintg = Dintg + (DENtemP / DenAvew- 1.0_WP) / DYFI(J)
!            WRITE(TECFLG_FavAG(1), '(13ES20.12)') YCC(J), (1.0_WP - DABS(YCC(J))) * REN * COE, COE, &
!                              FCT(J, 1:3), NSFbal_FA(J, 1:3), &
!                              Tau_Mean_RA(J, 1, 2)     * REN / Ret_ave_iO / Utaw_ave_io, &
!                             - Uff2d_FA  (J, 1, 2)     / Utaw_ave_iO / Utaw_ave_iO / DenAvew, &
!                              F_A / Utaw_ave_iO / Utaw_ave_io * Dintg, &
!                             -F_A / Utaw_ave_iO / Utaw_ave_iO / DYFI(J)
!        END DO
!        CLOSE(TECFLG_FavAG(1))



!!======================checking== Momentum Equation IS xyz direction ======================================================
!        FCT = 0.0_WP
!        DO J = 1, NCL2

!            FCT(J, 1) = - D1xztL_F0_io(J) * U1xztL_F0_io(J, 1) * U1xztL_F0_io(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew+ &
!                        Tau_Mean_RA(J, 1, 2) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew-&
!                        uf2d_RA  (J, 1, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew

!            FCT(J, 2) = - D1xztL_F0_io(J) * U1xztL_F0_io(J, 2) * U1xztL_F0_io(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew+ &
!                        Tau_Mean_RA(J, 2, 2) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew-&
!                        uf2d_RA  (J, 2, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew-&
!                        dPDX_RA(J, 2) / Utaw_ave_iO / Utaw_ave_iO / D1xztL_F0_io(J)!

!            FCT(J, 3) = - D1xztL_F0_io(J) * U1xztL_F0_io(J, 3) * U1xztL_F0_io(J, 2) / Utaw_ave_iO / Utaw_ave_iO / DenAvew+ &
!                        Tau_Mean_RA(J, 2, 3) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew-&
!                        uf2d_RA  (J, 2, 3) / Utaw_ave_iO / Utaw_ave_iO / DenAvew
!        END DO
!        !=== At wall surfaceS ====
!        FCT(0, 1) = Tauw_io(1) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew
!        FCT(0, 2) = 0.0_WP
!        FCT(0, 3) = 0.0_WP

!        FCT(NND2, 1) = Tauw_io(2) * REN / Ret_ave_iO / Utaw_ave_iO /VisAvew
!        FCT(NND2, 2) = 0.0_WP
!        FCT(NND2, 3) = 0.0_WP

!        !===gradient at cell centrE ====
!        NSFbal_RA = 0.0_WP
!        DO J = 2, NCL2 - 1
!            DO N = 1, NDV
!                NSFbal_RA(J, N) = ( ( YCL2ND_WFB(J + 1) * FCT(J, N) + YCL2ND_WFF(J + 1) * FCT(J + 1, N) ) -             &
!                                  ( YCL2ND_WFF(J) * FCT(J, N) + YCL2ND_WFB(J) * FCT(J - 1, N) ) ) * DYFI(J) + &
!                                F_A* D1xztL_F0_io(J) / DenAvew/ Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N)) - &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!                !!WRITE(*, *) N, J, ( ( YCL2ND_WFB(J + 1) * FCT(J, N) + YCL2ND_WFF(J + 1) * FCT(J + 1, N) ) -             &
!                !                  ( YCL2ND_WFF(J) * FCT(J, N) + YCL2ND_WFB(J) * FCT(J - 1, N) ) ) * DYFI(J), &
!                !                  F_A* D1xztL_F0_io(J) / DenAvew/ Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N)), &
!                !                  F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!            END DO
!        END DO
!        DO N = 1, NDV
!            NSFbal_RA(1,   N) = ( ( YCL2ND_WFB(1 + 1) * FCT(1, N) +                                           &
!                                   YCL2ND_WFF(1 + 1) * FCT(1 + 1, N) ) - FCT(0, N) ) * DYFI(1) +               &
!                                 F_A* D1xztL_F0_io(1) / DenAvew/ Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N)) - &
!                                 F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))

!            NSFbal_RA(NCL2, N) = ( FCT(NND2, N) - ( YCL2ND_WFF(NCL2) * FCT(NCL2, N) +                       &
!                                                 YCL2ND_WFB(NCL2) * FCT(NCL2 - 1, N) ) ) * DYFI(NCL2) +      &
!                                F_A* D1xztL_F0_io(NCL2) / DenAvew/ Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N)) - &
!                                F_A / Utaw_ave_iO / Utaw_ave_io * DBLE(IBuoF(N))
!        END DO

!        !=== INTegral along y direction ===============================
!        NSFbalt_RA = 0.0_WP
!        DO J = 1, NCL2
!            DO N = 1, NDV
!                NSFbalt_RA(N) = NSFbalt_RA(N) + NSFbal_RA(J, N) / RCCI1(J) / DYFI(J)
!            END DO
!            !!WRITE(*, *) J, NSFbalt_RA(1:3)
!        END DO
!!======================checking== Momentum Equation IS xyz direction ======================================================

!!====================Checking=================================
!        FCT = 0.0_WP
!        DO J = 1, NCL2
!            FCT(J) = - D1xztL_F0_io(J) * U_FA(J, 2) * H_FA(J) + &
!                      DTDLKxztL_F0_io(J, 2) * CTHECD - &
!                      uffhffd_FA  (J, 2)
!        END DO
!        !=== At wall surfaceS ====
!        FCT(0) = Qw(1)
!        FCT(NND2) = Qw(2)

!        !===gradient at cell centrE ====
!        ENEbal_FA = 0.0_WP
!        DO J = 2, NCL2 - 1

!            ENEbal_FA(J) = ( ( YCL2ND_WFB(J + 1) * FCT(J) + YCL2ND_WFF(J + 1) * FCT(J + 1) ) - &
!                            ( YCL2ND_WFF(J) * FCT(J) + YCL2ND_WFB(J) * FCT(J - 1) ) ) * DYFI(J) * (-1.0_WP)

!        END DO
!        ENEbal_FA(1) = ( ( YCL2ND_WFB(1 + 1) * FCT(1)  + YCL2ND_WFF(1 + 1) * FCT(1 + 1)    ) - FCT(0)   ) * DYFI(1) * (-1.0_WP)
!        ENEbal_FA(NCL2) = ( ( YCL2ND_WFF(NCL2) * FCT(NCL2) + YCL2ND_WFB(NCL2) * FCT(NCL2 - 1) ) - FCT(NND2)) * DYFI(NCL2) * ((-1.0_WP)**2)

!        !=== INTegral along y direction ===============================
!        ENEbalt_FA = 0.0_WP
!        DO J = 1, NCL2
!            ENEbalt_FA = ENEbalt_FA + ENEbal_FA(J) / RCCI1(J) / DYFI(J)
!            !!WRITE(*, *) J, NSFbalt_FA(1:3)
!        END DO

!!====================calcuate variables for RANS ============================================
!        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
!        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.ForRANS.Heat.Transfer.' // TRIM(PNTIM) // '.plt'
!        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (9 variables)" '
!        WRITE(TECFLG_FavAG(1), '(A, A)') &
!            'variables = "1Y", "2Y+", "3Utau", "4MutOverPrt(x)", "5MutOverPrt(y)", "6MutOverPrt(z)", ', &
!            ' "7Prt(x)", "8Prt(y)", "9Prt(z)"'
!        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

!        DO J = 1, NCL2

!            MutOverPrt(1:3) = uffhffd_FA(J, 1:3) / (dHDX_FA(J, 1:3) + REALMIN)

!            Pruv(1:3) = RANS_Mut(J, 1, 2) / MutOverPrt(1:3)

!            WRITE(TECFLG_FavAG(1), '(9ES20.12)') YCC(J), YpluSD(J), UtauSD(J), &
!                              MutOverPrt(1:3), Pruv(1:3)
!        END DO
!        CLOSE(TECFLG_FavAG(1))


!!======================checking...==========================================================
!        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
!        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.Checking.Heat.Transfer.' // TRIM(PNTIM) // '.plt'
!        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (21 variables)" '
!        WRITE(TECFLG_FavAG(1), '(A)') &
!            'variables = "1Y", "2Y+", "3Utau", "4Cp(<T>)", "5Cp(dh/ DT)"'
!        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

!        DO J = 1, NCL2
!            Cp_eval = dHDX_FA(J, 2) / DTDX(J, 2)

!            WRITE(TECFLG_FavAG(1), '(5ES20.12)') YCC(J), YpluSD(J), UtauSD(J), &
!                              CpT(J), Cp_eval

!        END DO
!        CLOSE(TECFLG_FavAG(1))
!======================checking...==========================================================



!        WRITE(PNTIM, '(1ES15.9)') PhyTIME_io
!        FLNM = TRIM(FilePath4) // 'Result.IO.' // TRIM(STDIM(1)) // '.Profile.Checking.EnergyBalance.FA.' // TRIM(PNTIM) // '.plt'
!        OPEN (TECFLG_FavAG(1), FILE = TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1), '(A)') 'TITLE = " Favre Averged Flow (4 variables)" '
!        WRITE(TECFLG_FavAG(1), '(6A)') &
!            'variables = "1Y", "2Y+", "3Utau", "4ENEbal"'
!        WRITE(TECFLG_FavAG(1), '(A)') 'ZONE T = "' // TRIM(ADJUSTL(zoneNameView)) // ' " '

!        DO J = 1, NCL2
!            WRITE(TECFLG_FavAG(1), '(4ES20.12)') YCC(J), YpluSD(J), UtauSD(J), ENEbal_FA(J)
!        END DO
!        CLOSE(TECFLG_FavAG(1))
