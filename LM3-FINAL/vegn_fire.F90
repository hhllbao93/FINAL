module vegn_fire_mod

#include "../shared/debug.inc"

! This stuff is boilerplate for making/reading namelists.
use mpp_mod, only: mpp_pe, mpp_root_pe
#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use constants_mod,   only: PI, VONKARM
use time_manager_mod, only : time_type, get_date, days_in_month
use fms_mod, only : &
     write_version_number, file_exist, check_nml_error, error_mesg, &
     close_file, stdlog, stdout, WARNING, FATAL

use land_io_mod,     only : external_ts_type, init_external_ts, del_external_ts, &
     read_external_ts, read_field
use land_debug_mod,  only : check_conservation, check_var_range, is_watch_point, &
                            carbon_cons_tol, &
                            is_watch_cell, force_watch_cell_burning
use land_data_mod,   only : land_state_type, lnd
use land_tile_mod,   only : land_tile_type, land_tile_carbon, &
                            land_tile_list_type, land_tile_enum_type, &
                            first_elmt, tail_elmt, next_elmt, &
                            operator(/=), operator(==), current_tile
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type, &
     OP_SUM

use vegn_data_mod,   only : agf_bs, fsc_liv, fsc_wood, fsc_froot, &
                            SP_C4GRASS, SP_C3GRASS, SP_TEMPDEC, SP_TROPICAL, SP_EVERGR, &
                            LU_CROP, LU_PAST, LU_NTRL, LU_SCND, &
                            TROP_SAV, TROP_FOR, thresh_trop_SavFor, &   ! SSR20150831
                            TROP_SHR, thresh_trop_ShrSav, &   ! SSR20160211
                            split_past_tiles   ! SSR20151118
use vegn_tile_mod,   only : vegn_tile_type
use soil_tile_mod,   only : soil_tile_type, soil_ave_theta1, soil_ave_theta2
use sphum_mod,       only : qscomp
use vegn_cohort_mod, only : vegn_cohort_type
use soil_mod,        only : add_root_litter
use soil_carbon_mod, only : add_litter, soil_carbon_option, poolTotalCarbon, &
                            remove_carbon_fraction_from_pool, &
                            SOILC_CENTURY, SOILC_CENTURY_BY_LAYER, SOILC_CORPSE, &
                            n_c_types
use vegn_harvesting_mod, only : adj_nppPrevDay

implicit none
private

! ==== public interfaces =====================================================
public  ::  vegn_fire_init
public  ::  update_fire_fast
public  ::  update_fire_agri
public  ::  update_fire_data ! reads external data for fire model
public  ::  update_fire_Fk
public  ::  update_fire_agb
public  ::  vegn_burn
public  ::  vegn_fire_sendtiledata_Cburned
public  ::  send_tile_data_BABF_forAgri
public  ::  fire_fragmentation
public  ::  defo_annCalcs

integer, public  ::  fire_option = 0
integer, public, parameter :: FIRE_NONE = 0, FIRE_LM3 = 1, FIRE_UNPACKED = 2
integer, public  ::  fire_option_past = 0
integer, public, parameter :: FIRE_PASTLI = 0, FIRE_PASTFP = 1
integer, public  ::  fire_windType = 0
integer, public, parameter :: FIRE_WIND_CANTOP = 0, FIRE_WIND_10MTMP = 1, FIRE_WIND_10MSHEFFIELD = 2
integer, public :: fire_option_fAGB = 0
integer, public :: fire_option_fRH = 0
integer, public :: fire_option_fTheta = 0
integer, public, parameter :: FIRE_AGB_LI2012 = 0, FIRE_AGB_LOGISTIC = 1, FIRE_AGB_GOMPERTZ = 2
integer, public, parameter :: FIRE_RH_LI2012 = 0, FIRE_RH_LOGISTIC = 1, FIRE_RH_GOMPERTZ = 2
integer, public, parameter :: FIRE_THETA_LI2012 = 0, FIRE_THETA_LOGISTIC = 1, FIRE_THETA_GOMPERTZ = 2
integer, public  ::  order_LU_vegn_fire_option = 0
integer, public, parameter :: ORDER_FLV = 0, ORDER_LVF = 1, ORDER_VFL = 2
integer, public, parameter :: NOW_DO_FIRELU = 0, NOW_DO_FIRE = 1, NOW_DO_LU = 2

! =====end of public interfaces ==============================================

! ==== constants =============================================================
character(len=*), parameter   :: &
     version     = '$Id: $', &
     tagname     = '$Name: $', &
     module_name = 'vegn_fire_mod', &
     diag_mod_name = 'vegn'

! ==== variables =============================================================

!---- namelist ---------------------------------------------------------------

! Overarching settings
character(32) :: fire_to_use = 'unpacked'
    ! 'LM3' for old once-a-year fire, same as always existed in LM3
    ! 'unpacked' for modern fire parameterization
character(32) :: fire_for_past = 'pastfp'
    ! 'pastfp' for unpacked fire based on Fpast
    ! 'pastli' for fire using NTRL/SCND model (after Li)
logical :: do_calc_derivs = .FALSE.
logical :: zero_tmf = .TRUE.    ! SSR20150921: If we have to reduce burned area it
                                ! exceeded the unburned area in the grid cell, or because
                                ! fire size was larger than that allowed by fragmentation,
                                ! then if TRUE, set all derivatives to zero.
character(3) :: order_LU_vegn_fire = 'flv'   ! Refers to order of Land-use transitions, 
                                             ! Fire, and update_Vegn slow as called in 
                                             ! update_land_model_slow_0d. Other valid
                                             ! options: 'lvf' and 'vfl'.
real   ::   magic_scalar = 1.0   ! SSR20160222 on advice of SLP.

! Create fire tiles?
logical, public :: do_fire_tiling = .TRUE.
real, public    :: min_BA_to_split = 0.0   ! Km2. If BA < than this, fire tile will not be split.
                                           ! Instead, combustion completeness and mortality will
                                           ! be multiplied by burned_frac and applied to the
                                           ! existing tile.
                                           
! Incorporate fragmentation effects?
logical, public  ::  do_fire_fragmentation = .FALSE.
real             ::  max_fire_size_min = 1e5   ! m2. Default taken from Pfeiffer et al. (2013).
logical          ::  frag_incl_PAST = .TRUE.     ! Include pasture as fragmenting? (When using unpacked PAST)
logical          ::  frag_incl_nonveg = .TRUE.   ! Include non-vegetated tiles as fragmenting?

! RH
logical :: use_Frh = .TRUE.
character(32) :: f_rh_style = 'li2012'   ! Or 'logistic' or 'gompertz'
real :: rh_up = 0.7   ! From Li et al. (2012)
real :: rh_lo = 0.3   ! From Li et al. (2012)
real :: rh_psi2 = -13.2522     ! For logistic-style fRH. Default value determined by MATLAB fitting.
real :: rh_psi3 = 0.5          ! For logistic-style fRH. Default value determined by MATLAB fitting.
real :: rh_gom2 = 0.0062     ! For Gompertz-style fRH. Default value determined by MATLAB fitting.
real :: rh_gom3 = -9.1912    ! For Gompertz-style fRH. Default value determined by MATLAB fitting.

! Soil moisture
character(32) :: f_theta_style = 'li2012'   ! Or 'logistic' or 'gompertz'. SSR20160203
logical :: use_Ftheta = .TRUE.
logical :: use_Cm = .TRUE.
real :: theta_extinction = 0.69 ! Previously 0.829583324318996
real :: theta_psi2 = -9.3182     ! For logistic-style fTheta. Default value determined by MATLAB fitting.
real :: theta_psi3 = 0.3338      ! For logistic-style fTheta. Default value determined by MATLAB fitting.
real :: theta_gom2 = 0.0750     ! For Gompertz-style fTheta. Default value determined by MATLAB fitting.
real :: theta_gom3 = -6.3741    ! For Gompertz-style fTheta. Default value determined by MATLAB fitting.
logical :: thetaE_already_squared = .FALSE.   ! If .TRUE., input theta_extinction is
                                              ! assumed to be squared already.
real :: C_beta_threshUP = 0.7   ! For calculating effect of soil moisture on fire ROS
real :: C_beta_threshLO = 0.3   ! For calculating effect of soil moisture on fire ROS
logical :: C_beta_params_likeRH = .FALSE.   ! If .TRUE., sets C_beta parameters equal to
                                            ! whatever RH_up and RH_lo are. This will then
                                            ! affect derivative w/r/t RH parameters.
logical :: theta_ROSeffect_asFnTheta = .FALSE.   ! SSR20151216. If .TRUE., overrides
                                                 ! C_beta_params_likeRH.
logical :: theta_include_ice = .TRUE.   ! FALSE reproduces original behavior
real    :: depth_for_theta = 0.05   ! -1.0 for original behavior. Anything positive for average down to that depth (m).

! Temperature (from Li et al. [2013])
logical :: use_Ftemp = .TRUE.
real :: T_up_celsius = 0.
real :: T_lo_celsius = -10.

! Biomass
character(32) :: f_agb_style = 'li2012'   ! Or 'logistic' or 'gompertz'
real :: fire_biomass_threshold = -1   ! -1 means calculate F_b, otherwise no burning if < this
real :: agb_up = 1.05   ! For li2012-style fAGB. Above this, biomass availability doesn't matter. (kg)
real :: agb_lo = 0.155   ! For li2012-style fAGB. Below this, no fire allowed. (kg)
real :: agb_psi2 = 5.8864  ! For logistic-style fAGB. Default value determined by MATLAB fitting.
real :: agb_psi3 = 0.6021  ! For logistic-style fAGB. Default value determined by MATLAB fitting.
real :: agb_gom2 = 7.3157  ! For Gompertz-style fAGB. Default value determined by MATLAB fitting.
real :: agb_gom3 = 4.1100  ! For Gompertz-style fAGB. Default value determined by MATLAB fitting.

! Fire shape (ellipse length:breadth)
real :: LBratio_a = 1.
real :: LBratio_b = 10.
real :: LBratio_c = 0.06

! Population density: Ignitions (from Li et al. [2012]).
real :: Ia_param1 = 6.8
real :: Ia_param2 = 0.6

! Population density: Suppression
logical :: use_FpopD_nf = .TRUE.
logical :: use_FpopD_ba = .TRUE.
real :: popD_supp_eps1 = 0.99   ! From Li et al. (2012)
real :: popD_supp_eps2 = 0.98   ! From Li et al. (2012)
real :: popD_supp_eps3 = 0.025  ! From Li et al. (2012)

! GDP suppression
logical :: use_Fgdp_nf = .TRUE.
logical :: use_Fgdp_ba = .TRUE.

! Can turn off lightning and/or human ignitions using these two parameters
real :: In_c2g_ign_eff = 0.25       ! From Li et al. (2012) + Corrigendum
real :: Ia_alpha_monthly = 0.0035   ! From Li et al. (2013). Ignitions/person/month; 
                                    ! converted to daily in fire_init.

! Fire duration. If species-specific value not specified in namelist,
! vegn_fire_init will set it to fire_duration_ave.
real :: fire_duration_ave = 60.*60.*24.   ! 24 hours. From Li et al. (2012)
real :: fire_duration_ave_c4 = -1.   ! C4 grass (vegn%species==0)
real :: fire_duration_ave_c3 = -1.   ! C3 grass (vegn%species==1)
real :: fire_duration_ave_dt = -1.   ! Deciduous tree (vegn%species==2)
real :: fire_duration_ave_tt = -1.   ! Tropical tree (vegn%species==3)
real :: fire_duration_ave_et = -1.   ! Evergreen tree (vegn%species==4)

! Wind
logical :: use_Fwind = .TRUE.   ! If FALSE, ROS_max even with wind=0 (although C_m still matters)
character(32) :: wind_to_use = '10m_sheffield'
    ! 'canopy_top' for original, top-of-canopy
    ! '10m_tmp' for kludgey 10-m wind (uses neutral stability assumption)
    ! '10m_sheffield' for Sheffield's 10-m wind
real    :: constant_wind_forFire = -1   ! -1 means use wind_forFire from canopy_air.F90. Anything else (x) sets wind_forFire to x m/s.

! Maximum rate of spread for each species. From Li et al. (2012); Doubled according to Corrigendum.
logical :: lock_ROSmaxC3_to_ROSmaxC4 = .FALSE.
real    :: ROS_max_C4GRASS  = 0.4
real    :: ROS_max_C3GRASS  = 0.4
real    :: ROS_max_TEMPDEC  = 0.22
! real    :: ROS_max_TROPICAL = 0.22
real    :: ROS_max_TROPSHR = -1      ! SSR20160211. vegn_fire_ROS will fall back on
                                     ! ROS_max_TROPICAL if ROS_max_TROPSHR is <=0.
real    :: ROS_max_TROPSAV = -1      ! SSR20150831. vegn_fire_ROS will fall back on
                                     ! ROS_max_TROPICAL if ROS_max_TROPSAV is <=0.
real    :: ROS_max_TROPICAL = 0.22   ! SSR20150831
real    :: ROS_max_EVERGR   = 0.3

! Dealing with low or high burning
real :: min_fire_size = 1.E-6   ! Minimum fire size (i.e., BAperFire): km2
real :: min_combined_ba = 1.E-6   ! Minimum total tile BA (km2)
logical :: print_min_fire_violations = .FALSE.

! Combustion completeness (SSR20160223)
logical :: burn_cwl_like_stem = .FALSE.   ! SSR20160410
real    :: CC_GRASS_leaf    = 0.85
real    :: CC_GRASS_stem    = 1.00   ! Li2012 has no value, but LM3 actually can have some wood in grassland.
real    :: CC_GRASS_litter    = 0.85
real    :: CC_TEMPDEC_leaf  = 0.70
real    :: CC_TEMPDEC_stem  = 0.10
real    :: CC_TEMPDEC_litter  = 0.45
real    :: CC_TROPSHRSAV_leaf = 0.70
real    :: CC_TROPSHRSAV_stem = 0.10
real    :: CC_TROPSHRSAV_litter = 0.50
real    :: CC_TROPICAL_leaf = 0.70
real    :: CC_TROPICAL_stem = 0.15
real    :: CC_TROPICAL_litter = 0.50
real    :: CC_EVERGR_leaf   = 0.75
real    :: CC_EVERGR_stem   = 0.20
real    :: CC_EVERGR_litter   = 0.55


! Post-fire mortality
real    :: fireMort_GRASS_leaf = 0.85
real    :: fireMort_GRASS_stem = 0.00      ! Li2012 has no value, but LM3 actually can have some wood in grassland.
real    :: fireMort_GRASS_root = 0.20
real    :: fireMort_TEMPDEC_leaf = 0.70
real    :: fireMort_TEMPDEC_stem = 0.55
real    :: fireMort_TEMPDEC_root = 0.07
real    :: fireMort_TROPSHRSAV_leaf = 0.70   ! SSR20160223
real    :: fireMort_TROPSHRSAV_stem = 0.55   ! SSR20160223
real    :: fireMort_TROPSHRSAV_root = 0.07   ! SSR20160223
real    :: fireMort_TROPICAL_leaf = 0.70
real    :: fireMort_TROPICAL_stem = 0.60
real    :: fireMort_TROPICAL_root = 0.10
real    :: fireMort_EVERGR_leaf = 0.75
real    :: fireMort_EVERGR_stem = 0.65
real    :: fireMort_EVERGR_root = 0.13

! For CLMCN-style deforestation (Kloster et al., 2010)
real    :: fs_min = 0.2
real    :: fs_max = 0.8
real    :: fp_lo = 0.01
real    :: fp_hi = 0.3

logical :: minimal_fire_diagnostics = .FALSE.   ! SSR20150812

namelist /fire_nml/ fire_to_use, fire_for_past, &
                    f_rh_style, rh_up, rh_lo, rh_psi2, rh_psi3, rh_gom2, rh_gom3, &
                    theta_extinction, &
                    thetaE_already_squared, &   ! SSR20151009
                    use_Cm, & ! SSR20160217
                    T_up_celsius, T_lo_celsius, &
                    f_agb_style, agb_up, agb_lo, agb_psi2, agb_psi3, agb_gom2, agb_gom3, &
                    LBratio_a, LBratio_b, LBratio_c, &
                    C_beta_threshUP, C_beta_threshLO, &
                    C_beta_params_likeRH, &   ! SSR20151009
                    popD_supp_eps1, popD_supp_eps2, popD_supp_eps3, &
                    In_c2g_ign_eff, &
                    Ia_alpha_monthly, Ia_param1, Ia_param2, &
                    fire_duration_ave, fire_duration_ave_c4, fire_duration_ave_c3, &
                    fire_duration_ave_dt, fire_duration_ave_tt, fire_duration_ave_et, &
                    min_fire_size, min_combined_ba, &
                    print_min_fire_violations, &
                    fire_biomass_threshold, wind_to_use, &
                    use_Ftheta, use_Frh, use_Ftemp, use_Fwind, &
                    use_FpopD_nf, use_FpopD_ba, use_Fgdp_nf, use_Fgdp_ba, &
                    theta_include_ice, depth_for_theta, &
                    fireMort_GRASS_leaf, fireMort_GRASS_stem, fireMort_GRASS_root, &
                    fireMort_TEMPDEC_leaf, fireMort_TEMPDEC_stem, fireMort_TEMPDEC_root, &
                    fireMort_TROPSHRSAV_leaf, fireMort_TROPSHRSAV_stem, fireMort_TROPSHRSAV_root, &
                    fireMort_TROPICAL_leaf, fireMort_TROPICAL_stem, fireMort_TROPICAL_root, &
                    fireMort_EVERGR_leaf, fireMort_EVERGR_stem, fireMort_EVERGR_root, &
                    CC_GRASS_leaf, CC_GRASS_stem, CC_GRASS_litter, &   ! SSR20160223
                    CC_TEMPDEC_leaf, CC_TEMPDEC_stem, CC_TEMPDEC_litter, &   ! SSR20160223
                    CC_TROPSHRSAV_leaf, CC_TROPSHRSAV_stem, CC_TROPSHRSAV_litter, &   ! SSR20160223
                    CC_TROPICAL_leaf, CC_TROPICAL_stem, CC_TROPICAL_litter, &   ! SSR20160223
                    CC_EVERGR_leaf, CC_EVERGR_stem, CC_EVERGR_litter, &   ! SSR20160223
                    constant_wind_forFire, &
                    ROS_max_C4GRASS, ROS_max_C3GRASS, ROS_max_TEMPDEC, ROS_max_TROPICAL, ROS_max_EVERGR, &
                    ROS_max_TROPSHR, &   ! SSR20160211
                    ROS_max_TROPSAV, &   ! SSR20150831
                    do_calc_derivs, do_fire_tiling, &
                    order_LU_vegn_fire, &
                    min_BA_to_split, &
                    minimal_fire_diagnostics, &
                    do_fire_fragmentation, max_fire_size_min, &
                    frag_incl_PAST, frag_incl_nonveg, &
                    fs_min, fs_max, fp_lo, fp_hi, &   ! SSR20150903
                    zero_tmf, &   ! SSR20150921
                    theta_ROSeffect_asFnTheta, &   ! SSR20151216
                    lock_ROSmaxC3_to_ROSmaxC4, &   ! SSR20160131
                    f_theta_style, theta_psi2, theta_psi3, theta_gom2, theta_gom3, &   ! SSR20160203
                    magic_scalar, &   ! SSR20160222
                    burn_cwl_like_stem   ! SSR20160413

!---- end of namelist --------------------------------------------------------
real :: dt_fast      ! fast time step, s
real :: dt_fast_yr   ! fast time step, yr
real :: dt_slow      ! slow time step, s
real :: dt_slow_yr   ! slow time step, yr
real :: seconds_per_year, days_per_year

real :: LB_max = -1 ! Maximum length:breadth ratio
real :: HB_max = -1 ! Maximum head:back ratio
real :: Ia_alpha_daily = -1 ! Daily per-person ignition rate

real, dimension(0:4) :: fire_duration_array

! Placeholders for derivatives
real :: Ia_DERIVwrt_alphaM            = -1.0e+20
real :: Ia_DERIVwrt_IaParam1          = -1.0e+20
real :: Ia_DERIVwrt_IaParam2          = -1.0e+20
real :: fire_fn_theta_DERIVwrt_thetaE = -1.0e+20
real :: fire_TOTALfn_theta_DERIVwrt_thetaE = -1.0e+20
real :: fire_TOTALfn_theta_DERIVwrt_param1 = -1.0e+20
real :: fire_TOTALfn_theta_DERIVwrt_param2 = -1.0e+20
real :: fire_fn_agb_DERIVwrt_param1   = -1.0e+20
real :: fire_fn_agb_DERIVwrt_param2   = -1.0e+20
real :: fire_fn_rh_DERIVwrt_param1    = -1.0e+20
real :: fire_fn_rh_DERIVwrt_param2    = -1.0e+20
real :: fire_TOTALfn_rh_DERIVwrt_param1 = -1.0e+20
real :: fire_TOTALfn_rh_DERIVwrt_param2 = -1.0e+20
real :: popDsupp_NF_DERIVwrt_eps1     = -1.0e+20
real :: popDsupp_NF_DERIVwrt_eps2     = -1.0e+20
real :: popDsupp_NF_DERIVwrt_eps3     = -1.0e+20
real :: ROS_DERIVwrt_ROSmax           = -1.0e+20
real :: BAperFire0_DERIVwrt_ROSmax    = -1.0e+20
real :: BAperFire0_DERIVwrt_fireDur   = -1.0e+20

type(external_ts_type) :: &
    population_ts, & ! input time series of population density
!SSR20151124    lightning_ts, &  ! input time series of lightning
    GDPpc_billion_ts, &        ! input time series of GDP
    Fc_ts, &         ! input time series of Fc
    Fp_ts            ! input time series of Fp
real, allocatable :: &
    population_in(:,:), & ! input buffer for population data
    lightning_in(:,:),  & ! input buffer for lightning data
    GDPpc_billion_in(:,:), &! input buffer for GDP data
    Fc_in(:,:), &         ! input buffer for Fc
    Fp_in(:,:), &         ! input buffer for Fp
    crop_burn_rate_in(:,:,:),& ! input buffer for monthly Fc values
    past_burn_rate_in(:,:,:),&   ! input buffer for monthly Fp values
    lightning_in_v2(:,:,:)

integer :: & ! diag field IDs
    id_population, id_lightning, id_GDPpc, &
    id_fire_fn_rh, id_fire_fn_theta, id_fire_fn_Tca, id_fire_fn_agb, &
    id_ROS, id_LB, id_HB, id_gW, &
    id_BAperFire_0, id_BAperFire_1, id_BAperFire_2, id_BAperFire_3, &   ! SSR20150812
    id_BAperFire_max, &   ! SSR20150727
    id_fragmenting_frac, id_burnable_frac, &   ! SSR20150811
    id_fire_duration_ave, &
    id_fire_fn_popD_NF, id_fire_fn_popD_BA, id_fire_fn_GDPpc_NF, id_fire_fn_GDPpc_BA, &
    id_Ia, id_In, id_Nfire_perKm2, id_Nfire_rate, &
    id_Fcrop, id_Fpast, &
    id_burn_Cemit, id_burn_Ckill, &
    id_burn_Cemit_noCWL, &   ! SSR20151227
    id_burn_Cemit_leaf,id_burn_Cemit_stem,id_burn_Cemit_litter, &   ! SSR20160224
    id_burn_Cemit_litter_lf,id_burn_Cemit_litter_fw,id_burn_Cemit_litter_cw, &              ! SSR20160224
    id_burn_Ckill_leaf,id_burn_Ckill_stem,id_burn_Ckill_root, &   ! SSR20160224
    id_fire_wind_forFire, &   ! Used for troubleshooting
    id_fire_agb, &   ! Used for troubleshooting
    id_fire_q, id_fire_rh, id_fire_theta, id_fire_Tca, &   ! Used for troubleshooting
    id_fire_depth_ave, &   ! Used for troubleshooting
    id_BF_rate, id_BA_rate, &
    id_BA_DERIVwrt_alphaM, id_BA_DERIVwrt_thetaE, &
    id_BA_DERIVwrt_IaParam1, id_BA_DERIVwrt_IaParam2, &
    id_BA_DERIVwrt_AGBparam1, id_BA_DERIVwrt_AGBparam2, &
    id_BA_DERIVwrt_RHparam1, id_BA_DERIVwrt_RHparam2, &
    id_BA_DERIVwrt_popDsupp_NF_eps1, id_BA_DERIVwrt_popDsupp_NF_eps2, id_BA_DERIVwrt_popDsupp_NF_eps3, &
    id_BA_DERIVwrt_fireDur_c4, id_BA_DERIVwrt_fireDur_c3, &
    id_BA_DERIVwrt_fireDur_dt, id_BA_DERIVwrt_fireDur_tt, id_BA_DERIVwrt_fireDur_et, &
    id_BA_DERIVwrt_ROSmax_c4, id_BA_DERIVwrt_ROSmax_c3, &
    id_BA_DERIVwrt_ROSmax_dt, id_BA_DERIVwrt_ROSmax_tt, id_BA_DERIVwrt_ROSmax_et, &
!    id_BA_DERIVwrt_ROSmax_ts, &   ! SSR20150831: Tropical savanna
    id_BA_DERIVwrt_ROSmax_tshr, id_BA_DERIVwrt_ROSmax_tsav, &   ! SSR20160211
    id_Ia_DERIVwrt_alphaM, id_fire_fn_theta_DERIVwrt_thetaE, &
    id_BA_DERIVwrt_THETAparam1, id_BA_DERIVwrt_THETAparam2, &   ! SSR20160203
    id_BA_DERIVwrt_magicScalar, &   ! SSR20160222
    id_Ia_DERIVwrt_IaParam1, id_Ia_DERIVwrt_IaParam2, &
    id_fire_fn_agb_DERIVwrt_param1, id_fire_fn_agb_DERIVwrt_param2, &
    id_popDsupp_NF_DERIVwrt_eps1, id_popDsupp_NF_DERIVwrt_eps2, id_popDsupp_NF_DERIVwrt_eps3, &
    id_fire_fn_rh_DERIVwrt_param1, id_fire_fn_rh_DERIVwrt_param2, &
    id_BAperFire0_DERIVwrt_fireDur, id_BAperFire0_DERIVwrt_ROSmax, &
    id_tropType   ! SSR20160211

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ==============================================================================
subroutine vegn_fire_init(id_lon, id_lat, dt_fast_in, dt_fast_yr_in, dt_slow_in, dt_slow_yr_in, time)
  integer,         intent(in) :: id_lon, id_lat ! ids of diagnostic axes
  real,            intent(in) :: dt_fast_in        ! fast time step, seconds
  real,            intent(in) :: dt_fast_yr_in     ! fast time step, years
  real,            intent(in) :: dt_slow_in        ! slow time step, seconds
  real,            intent(in) :: dt_slow_yr_in     ! slow time step, years
  type(time_type), intent(in) :: time           ! current time
  
  integer :: io, ierr, unit
  integer :: i, j   ! SSR
  logical, allocatable :: burn_rate_mask(:,:,:)
  logical, allocatable :: lightning_mask_v2(:,:,:)   ! SSR20151124

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=fire_nml, iostat=io)
  ierr = check_nml_error(io, 'fire_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=fire_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'fire_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=fire_nml)
  endif
  
  ! get dt_*
  dt_fast    = dt_fast_in
  dt_fast_yr = dt_fast_yr_in
  dt_slow    = dt_slow_in
  dt_slow_yr = dt_slow_yr_in
  seconds_per_year = dt_fast / dt_fast_yr
  days_per_year = seconds_per_year / 86400.

  ! parse the fire model options for efficiency.
  select case (trim(fire_to_use))
  case ('none')
     fire_option = FIRE_NONE
  case ('lm3')
     fire_option = FIRE_LM3
  case ('unpacked')
     fire_option = FIRE_UNPACKED
  case default
     call error_mesg('vegn_fire_init',&
        'option fire_to_use="'//trim(fire_to_use)//'" is invalid, use "lm3", "unpacked", or "none"', &
        FATAL)
  end select
  
  select case (trim(fire_for_past))
  case ('pastfp')
     fire_option_past = FIRE_PASTFP
  case ('pastli')
     fire_option_past = FIRE_PASTLI
  case default
     call error_mesg('vegn_fire_init',&
        'option fire_for_past="'//trim(fire_for_past)//'" is invalid, use "pastfp" or "pastli"', &
        FATAL)
  end select
  
  select case (trim(wind_to_use))
  case ('canopy_top')
     fire_windType = FIRE_WIND_CANTOP
  case ('10m_tmp')
     fire_windType = FIRE_WIND_10MTMP
  case ('10m_sheffield')
     fire_windType = FIRE_WIND_10MSHEFFIELD
  case default
     call error_mesg('vegn_fire_init',&
        'option wind_to_use="'//trim(wind_to_use)//'" is invalid, use "canopy_top", "10m_tmp", or "10m_sheffield"', &
        FATAL)
  end select
      
  select case (trim(f_agb_style))
  case ('li2012')
     fire_option_fAGB = FIRE_AGB_LI2012
  case ('logistic')
     fire_option_fAGB = FIRE_AGB_LOGISTIC
  case ('gompertz')
     fire_option_fAGB = FIRE_AGB_GOMPERTZ
  case default
     call error_mesg('vegn_fire_init',&
        'option fire_option_fAGB="'//trim(f_agb_style)//'" is invalid, use "li2012", "logistic", or "gompertz".', &
        FATAL)
  end select
  
  select case (trim(f_rh_style))
  case ('li2012')
     fire_option_fRH = FIRE_RH_LI2012
  case ('logistic')
     fire_option_fRH = FIRE_RH_LOGISTIC
  case ('gompertz')
     fire_option_fRH = FIRE_RH_GOMPERTZ
  case default
     call error_mesg('vegn_fire_init',&
        'option fire_option_fRH="'//trim(f_rh_style)//'" is invalid, use "li2012", "logistic", or "gompertz"', &
        FATAL)
  end select
  
  select case (trim(f_theta_style))
  case ('li2012')
     fire_option_fTheta = FIRE_THETA_LI2012
  case ('logistic')
     fire_option_fTheta = FIRE_THETA_LOGISTIC
  case ('gompertz')
     fire_option_fTheta = FIRE_THETA_GOMPERTZ
  case default
     call error_mesg('vegn_fire_init',&
        'option fire_option_fTheta="'//trim(f_theta_style)//'" is invalid, use "li2012", "logistic", or "gompertz"', &
        FATAL)
  end select
  
  ! parse order_LU_vegn_fire for efficiency
  select case (trim(order_LU_vegn_fire))
  case ('flv')
     order_LU_vegn_fire_option = ORDER_FLV
  case ('lvf')
     order_LU_vegn_fire_option = ORDER_LVF
  case ('vfl')
     order_LU_vegn_fire_option = ORDER_VFL
  case default
     call error_mesg('vegn_fire_init',&
        'option order_LU_vegn_fire="'//trim(order_LU_vegn_fire)//'" is invalid, use "flv", "lvf", or "vfl"', &
        FATAL)
  end select
  
  ! SSR20151009
  if (C_beta_params_likeRH) write(*,*) 'C_beta_params_likeRH is .TRUE., so ignoring setting of C_beta_threshUP and C_beta_threshLO.'
  
  ! parse species-specific values for fire_duration
  if (fire_duration_ave_c4 == -1.) fire_duration_ave_c4 = fire_duration_ave
  if (fire_duration_ave_c3 == -1.) fire_duration_ave_c3 = fire_duration_ave
  if (fire_duration_ave_dt == -1.) fire_duration_ave_dt = fire_duration_ave
  if (fire_duration_ave_tt == -1.) fire_duration_ave_tt = fire_duration_ave
  if (fire_duration_ave_et == -1.) fire_duration_ave_et = fire_duration_ave
  fire_duration_array(SP_C4GRASS) = fire_duration_ave_c4
  fire_duration_array(SP_C3GRASS) = fire_duration_ave_c3
  fire_duration_array(SP_TEMPDEC) = fire_duration_ave_dt
  fire_duration_array(SP_TROPICAL) = fire_duration_ave_tt
  fire_duration_array(SP_EVERGR) = fire_duration_ave_et
    
  ! fire_biomass_threshold can only be -1 or >=0
  if (fire_biomass_threshold /= -1) then
     call check_var_range(fire_biomass_threshold, 0.0,10.**37,'vegn_fire_init', 'fire_biomass_threshold', lnd%time, FATAL)
  endif
  
  ! depth_for_theta can only be -1 or positive
  if (depth_for_theta /= -1.0) then
     call check_var_range(depth_for_theta, 10.**-37., 10.**37., 'vegn_fire_init', 'depth_for_theta', lnd%time, FATAL)
  endif
  
  ! constant_wind_forFire can only be -1 or >=0
  if (constant_wind_forFire /= -1.0) then
     call check_var_range(constant_wind_forFire, 0.0, 10.**37., 'vegn_fire_init', 'constant_wind_forFire', lnd%time, FATAL)
  endif
    
  if (fire_option /= FIRE_UNPACKED) return ! nothing more to do
  
  ! For rate-of-spread calculation
  LB_max = LBratio_a + LBratio_b
  HB_max = (LB_max + sqrt(LB_max**2. - 1.)) / (LB_max - sqrt(LB_max**2. - 1.))
  
  ! Find daily per-person ignition rate
  Ia_alpha_daily = Ia_alpha_monthly * 12./days_per_year
  
  ! SSR20160131
  if (lock_ROSmaxC3_to_ROSmaxC4) then
     ROS_max_C3GRASS = ROS_max_C4GRASS
  endif
  
  ! initialize external fields
  ! SSR: Does horizontal interpolation
  if (use_FpopD_nf .OR. use_FpopD_ba .OR. Ia_alpha_monthly>0.0) then
     call init_external_ts(population_ts, 'INPUT/population.nc', 'pop_density',&
          lnd%lonb, lnd%latb, 'conservative', lnd%domain)
  endif
  if (use_Fgdp_nf .OR. use_Fgdp_ba) then
     call init_external_ts(GDPpc_billion_ts, 'INPUT/GDP.nc', 'GDPPC', &
          lnd%lonb, lnd%latb, 'conservative', lnd%domain)
  endif

  ! allocate buffers for external data
  allocate(lightning_in(lnd%is:lnd%ie,lnd%js:lnd%je))
!  if (use_FpopD_nf .OR. use_FpopD_ba) allocate(population_in(lnd%is:lnd%ie,lnd%js:lnd%je))
!  if (use_Fgdp_nf .OR. use_Fgdp_ba) allocate(GDPpc_billion_in(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate(population_in(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate(GDPpc_billion_in(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate(Fc_in(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate(Fp_in(lnd%is:lnd%ie,lnd%js:lnd%je))
  ! monthly buffers for the input burn rate data
  allocate(crop_burn_rate_in(lnd%is:lnd%ie,lnd%js:lnd%je,12))
  allocate(past_burn_rate_in(lnd%is:lnd%ie,lnd%js:lnd%je,12))
  allocate(burn_rate_mask(lnd%is:lnd%ie,lnd%js:lnd%je,12))
  
  ! SSR20151124: Monthly buffers for lightning
  allocate(lightning_in_v2(lnd%is:lnd%ie,lnd%js:lnd%je,12))
  allocate(lightning_mask_v2(lnd%is:lnd%ie,lnd%js:lnd%je,12))
  call read_field('INPUT/lightning.nc', 'LRMTS_COM_FR', lnd%lonb, lnd%latb, lightning_in_v2, &
                  interp = 'conservative', mask=lightning_mask_v2)
  where (.not.lightning_mask_v2) lightning_in_v2 = 0.0
  deallocate(lightning_mask_v2)

  ! read the monthly data
  call read_field('INPUT/Fk.nc', 'Fcrop', lnd%lonb, lnd%latb, crop_burn_rate_in, &
                  interp = 'conservative', mask=burn_rate_mask)
  where (.not.burn_rate_mask) crop_burn_rate_in = 0.0
  call read_field('INPUT/Fk.nc', 'Fpast', lnd%lonb, lnd%latb, past_burn_rate_in, &
                  interp = 'conservative', mask=burn_rate_mask)
  where (.not.burn_rate_mask) past_burn_rate_in = 0.0
  deallocate(burn_rate_mask)
  
  ! get the data for current date
  call update_fire_data(time,lnd%is,lnd%ie,lnd%js,lnd%je)
  
  ! ---- initialize the diagnostics --------------------------------------------
  
  ! Yearly (although calculated and/or saved daily)
!  if (use_FpopD_nf .OR. use_FpopD_ba) then
     id_population = register_tiled_diag_field (diag_mod_name, 'population', (/id_lon,id_lat/), &
          time, 'population density', 'person/km2', missing_value=-999.0)
!  endif
  id_fire_fn_popD_NF = register_tiled_diag_field (diag_mod_name, 'fire_fn_popD_NF',(/id_lon,id_lat/), &
       time, 'Fraction of fires suppressed by popD', 'dimensionless', &
       missing_value=-1.0)
  id_fire_fn_popD_BA = register_tiled_diag_field (diag_mod_name, 'fire_fn_popD_BA',(/id_lon,id_lat/), &
       time, 'Fraction of potential BA realized via popD', 'dimensionless', &
       missing_value=-1.0)
!  if (use_Fgdp_nf .OR. use_Fgdp_ba) then
     id_GDPpc = register_tiled_diag_field (diag_mod_name, 'GDPpc', (/id_lon,id_lat/), &
          time, 'gross domestic product per capita', 'dollars/person', missing_value=-999.0)
!  endif
  id_fire_fn_GDPpc_NF = register_tiled_diag_field (diag_mod_name, 'fire_fn_GDPpc_NF',(/id_lon,id_lat/), &
       time, 'Fraction of ignitions becoming fires: GDP per capita', 'dimensionless', &
       missing_value=-1.0)
  id_fire_fn_GDPpc_BA = register_tiled_diag_field (diag_mod_name, 'fire_fn_GDPpc_BA',(/id_lon,id_lat/), &
       time, 'Fraction of potential BA realized via GDP per capita', 'dimensionless', &
       missing_value=-1.0)
  id_Ia = register_tiled_diag_field (diag_mod_name, 'Ia',(/id_lon,id_lat/), &
       time, 'Number of ignitions: Anthropogenic', 'Ignitions/km2', &
       missing_value=-1.0)

  ! Daily (saved daily)
  id_lightning = register_tiled_diag_field (diag_mod_name, 'lightning', (/id_lon,id_lat/), &
       time, 'lightning strike density', 'flashes/km2/day', missing_value=-999.0)
  id_In = register_tiled_diag_field (diag_mod_name, 'In',(/id_lon,id_lat/), &
       time, 'Number of ignitions: Natural', 'Ignitions/km2', &
       missing_value=-1.0)
!  id_BAperFire = register_tiled_diag_field (diag_mod_name, 'BAperFire',(/id_lon,id_lat/), &
!       time, 'Burned area per fire', 'km2/fire', &
!       missing_value=-1.0)
!  id_BAperFire_noAnthro = register_tiled_diag_field (diag_mod_name, 'BAperFire_noAnthro',(/id_lon,id_lat/), &
!       time, 'Average burned area per fire, without popD/GDP effects', 'km2/fire', &
!       missing_value=-1.0)
  id_BAperFire_0 = register_tiled_diag_field (diag_mod_name, 'BAperFire_0',(/id_lon,id_lat/), &
       time, 'Burned area per fire, before adjusting for tile size', 'km2/fire', &
       missing_value=-1.0)
  id_BAperFire_1 = register_tiled_diag_field (diag_mod_name, 'BAperFire_1',(/id_lon,id_lat/), &
       time, 'Burned area per fire, after adjusting for tile size', 'km2/fire', &
       missing_value=-1.0)
  id_BAperFire_2 = register_tiled_diag_field (diag_mod_name, 'BAperFire_2',(/id_lon,id_lat/), &
       time, 'Burned area per fire, after max_fire_size limitation', 'km2/fire', &
       missing_value=-1.0)
  id_BAperFire_3 = register_tiled_diag_field (diag_mod_name, 'BAperFire_3',(/id_lon,id_lat/), &
       time, 'Burned area per fire, after suppression', 'km2/fire', &
       missing_value=-1.0)
  id_BAperFire_max = register_tiled_diag_field (diag_mod_name, 'BAperFire_max',(/id_lon,id_lat/), &
       time, 'maximum fire size', 'km2', &
       missing_value=-1.0)
  id_fragmenting_frac = register_tiled_diag_field (diag_mod_name, 'fragmenting_frac',(/id_lon,id_lat/), &
       time, 'fraction of grid cell contributing to fragmentation', 'unitless', &
       missing_value=-1.0)
  id_burnable_frac = register_tiled_diag_field (diag_mod_name, 'burnable_frac',(/id_lon,id_lat/), &
       time, 'fraction of grid cell that could burn but is getting fragmented', 'unitless', &
       missing_value=-1.0)
  id_fire_agb = register_tiled_diag_field (diag_mod_name, 'fire_agb',(/id_lon,id_lat/), &
       time, 'Average aboveground biomass for fire calc', 'kgC/m2', &
       missing_value=-1.0)
  id_Nfire_perKm2 = register_tiled_diag_field (diag_mod_name, 'Nfire_perKm2',(/id_lon,id_lat/), &
       time, 'Number of fires per km2', 'Fires/km2', &
       missing_value=-1.0)
  id_Nfire_rate = register_tiled_diag_field (diag_mod_name, 'Nfire_rate',(/id_lon,id_lat/), &
       time, 'Number of fires in tile', 'Fires/day', &
       missing_value=-1.0, op=OP_SUM)
  id_BF_rate = register_tiled_diag_field (diag_mod_name, 'BF_rate',(/id_lon,id_lat/), &
     time, 'Burned fraction', '/day', &
     missing_value=-1.0)
  id_BA_rate = register_tiled_diag_field (diag_mod_name, 'BA_rate',(/id_lon,id_lat/), &
     time, 'Burned area', 'km2/day', &
     missing_value=-1.0, op=OP_SUM)
  id_Fcrop = register_tiled_diag_field (diag_mod_name, 'Fcrop',(/id_lon,id_lat/), &
       time, 'Fraction of crop that burns', '/day', &
       missing_value=-1.0)
  id_Fpast = register_tiled_diag_field (diag_mod_name, 'Fpast',(/id_lon,id_lat/), &
       time, 'Fraction of pasture that burns', '/day', &
       missing_value=-1.0)
  id_burn_Cemit = register_tiled_diag_field (diag_mod_name, 'burn_Cemit',(/id_lon,id_lat/), &
       time, 'Biomass combusted by fire', 'kg C/day', &
       missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_noCWL = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_noCWL',(/id_lon,id_lat/), &
       time, 'Biomass combusted by fire, excluding coarse woody litter', 'kg C/day', &
       missing_value=-1.0, op=OP_SUM)
  id_burn_Ckill = register_tiled_diag_field (diag_mod_name, 'burn_Ckill',(/id_lon,id_lat/), &
       time, 'Biomass killed by fire', 'kg C/day', &
       missing_value=-1.0, op=OP_SUM)
!  id_tropType = register_tiled_diag_field (diag_mod_name, 'tropType',(/id_lon,id_lat/), &
!       time, 'SP_TROPICAL, but what kind?', 'integer', &
!       missing_value=-1.0)

  ! SSR20160224
  id_burn_Cemit_leaf = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_leaf',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Leaf', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_stem = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_stem',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Stem', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_litter = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_litter',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Litter', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_litter_lf = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_litter_lf',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Leaf litter', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_litter_fw = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_litter_fw',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Fine woody litter', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Cemit_litter_cw = register_tiled_diag_field (diag_mod_name, 'burn_Cemit_litter_cw',(/id_lon,id_lat/), &
     time, 'Biomass combusted by fire: Coarse woody litter', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Ckill_leaf = register_tiled_diag_field (diag_mod_name, 'burn_Ckill_leaf',(/id_lon,id_lat/), &
     time, 'Biomass killed by fire: Leaf', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Ckill_stem = register_tiled_diag_field (diag_mod_name, 'burn_Ckill_stem',(/id_lon,id_lat/), &
     time, 'Biomass killed by fire: Stem', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  id_burn_Ckill_root = register_tiled_diag_field (diag_mod_name, 'burn_Ckill_root',(/id_lon,id_lat/), &
     time, 'Biomass killed by fire: Fine root', 'kg C/day', &
     missing_value=-1.0, op=OP_SUM)
  
  ! Sub-daily, instant values (saved sub-daily)
  id_fire_q = register_tiled_diag_field (diag_mod_name, 'fire_q',(/id_lon,id_lat/), &
       time, 'Specific humidity for fire (canopy)', 'kg/kg', &
       missing_value=-1.0)
  id_fire_rh = register_tiled_diag_field (diag_mod_name, 'fire_rh',(/id_lon,id_lat/), &
       time, 'Relative humidity for fire', 'dimensionless', &
       missing_value=-1.0)
  id_fire_theta = register_tiled_diag_field (diag_mod_name, 'fire_theta',(/id_lon,id_lat/), &
       time, 'Soil moisture for fire', 'dimensionless', &
       missing_value=-1.0)
  id_fire_Tca = register_tiled_diag_field (diag_mod_name, 'fire_Tca',(/id_lon,id_lat/), &
       time, 'Temperature for fire', 'degC', &
       missing_value=-1.0)
  id_fire_fn_theta = register_tiled_diag_field (diag_mod_name, 'fire_fn_theta',(/id_lon,id_lat/), &
       time, 'Fraction of ignitions becoming fires: Theta', 'dimensionless', &
       missing_value=-1.0)
  id_fire_fn_rh = register_tiled_diag_field (diag_mod_name, 'fire_fn_rh',(/id_lon,id_lat/), &
       time, 'Fraction of ignitions becoming fires: RH', 'dimensionless', &
       missing_value=-1.0)
  id_fire_fn_Tca = register_tiled_diag_field (diag_mod_name, 'fire_fn_Tca',(/id_lon,id_lat/), &
       time, 'Fraction of ignitions becoming fires: Canopy air temp.', 'dimensionless', &
       missing_value=-1.0)
  id_fire_fn_agb = register_tiled_diag_field (diag_mod_name, 'fire_fn_agb',(/id_lon,id_lat/), &
       time, 'Fraction of ignitions becoming fires: AGB', 'dimensionless', &
       missing_value=-1.0)
  id_fire_wind_forFire = register_tiled_diag_field (diag_mod_name, 'wind_forFire',(/id_lon,id_lat/), &
       time, 'Wind speed for use in fire model', 'm/s', &
       missing_value=-1.0)
  id_ROS = register_tiled_diag_field (diag_mod_name, 'ROS',(/id_lon,id_lat/), &
       time, 'Fire rate of spread', 'm/s', &
       missing_value=-1.0)
  id_fire_duration_ave = register_tiled_diag_field (diag_mod_name, 'fire_duration_ave',(/id_lon,id_lat/), &
       time, 'Average fire duration', 's', &
       missing_value=-1.0)
  id_LB = register_tiled_diag_field (diag_mod_name, 'LB',(/id_lon,id_lat/), &
       time, 'Length:breadth ratio', 'dimensionless', &
       missing_value=-1.0)
  id_HB = register_tiled_diag_field (diag_mod_name, 'HB',(/id_lon,id_lat/), &
       time, 'Head:back ratio', 'dimensionless', &
       missing_value=-1.0)
  id_gW = register_tiled_diag_field (diag_mod_name, 'gW',(/id_lon,id_lat/), &
       time, 'g(W)', 'dimensionless', &
       missing_value=-1.0)
  id_fire_depth_ave = register_tiled_diag_field (diag_mod_name, 'fire_depth_ave',(/id_lon,id_lat/), &
       time, 'depth_ave for theta calculation for fire', 'm', &
       missing_value=-1.0)
  id_BA_DERIVwrt_alphaM = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_alphaM',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t alpha_m', '(km2/day)/(alpha_monthly_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_thetaE = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_thetaE',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t theta_extinction', '(km2/day)/(theta_extinction_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
     
  ! SSR20160203
  id_BA_DERIVwrt_THETAparam1 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_THETAparam1',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t theta_psi2 or theta_gom2', '(km2/day)/(THETAparam1_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_THETAparam2 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_THETAparam2',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t theta_psi3 or theta_gom3', '(km2/day)/(THETAparam2_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
     
  id_BA_DERIVwrt_AGBparam1 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_AGBparam1',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t agb_lo, agb_psi2, or agb_gom2', '(km2/day)/(AGBparam1_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_AGBparam2 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_AGBparam2',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t agb_up, agb_psi3, or agb_gom3', '(km2/day)/(AGBparam2_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_RHparam1 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_RHparam1',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t rh_lo, rh_psi2, or rh_gom2', '(km2/day)/(RHparam1_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_RHparam2 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_RHparam2',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t rh_up, rh_psi3, or rh_gom3', '(km2/day)/(RHparam2_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_IaParam1 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_IaParam1',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t Ia param1', '(km2/day)/(IaParam1_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_IaParam2 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_IaParam2',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t Ia param2', '(km2/day)/(IaParam2_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_popDsupp_NF_eps1 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_popDsupp_NF_eps1',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t popDsupp_NF eps1', '(km2/day)/(eps1_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_popDsupp_NF_eps2 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_popDsupp_NF_eps2',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t popDsupp_NF eps2', '(km2/day)/(eps2_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_popDsupp_NF_eps3 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_popDsupp_NF_eps3',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t popDsupp_NF eps3', '(km2/day)/(eps3_unit)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_fireDur_c4 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_fireDur_c4',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t mean fire duration for C4 grass', '(km2/day)/second', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_fireDur_c3 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_fireDur_c3',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t mean fire duration for C3 grass', '(km2/day)/second', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_fireDur_dt = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_fireDur_dt',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t mean fire duration for deciduous trees', '(km2/day)/second', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_fireDur_tt = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_fireDur_tt',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t mean fire duration for tropical trees', '(km2/day)/second', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_fireDur_et = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_fireDur_et',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t mean fire duration for evergreen trees', '(km2/day)/second', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_c4 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_c4',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for C4 grass', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_c3 = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_c3',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for C3 grass', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_dt = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_dt',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for deciduous trees', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_tt = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_tt',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for tropical forest', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
!  id_BA_DERIVwrt_ROSmax_ts = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_ts',(/id_lon,id_lat/), &
!     time, 'Partial derivative of BA w/r/t max. rate of spread for tropical savanna', '(km2/day)/(m/s)', &
!     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_tshr = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_tshr',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for tropical shrubland', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_tsav = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_tsav',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for tropical savanna', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_ROSmax_et = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_ROSmax_et',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t max. rate of spread for evergreen trees', '(km2/day)/(m/s)', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BA_DERIVwrt_magicScalar = register_tiled_diag_field (diag_mod_name, 'BA_DERIVwrt_magicScalar',(/id_lon,id_lat/), &
     time, 'Partial derivative of BA w/r/t magic scalar', 'km2/day', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_Ia_DERIVwrt_alphaM = register_tiled_diag_field (diag_mod_name, 'Ia_DERIVwrt_alphaM',(/id_lon,id_lat/), &
     time, 'Partial derivative of Ia w/r/t alphaM', '(ignitions/km2)/alphaM_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_Ia_DERIVwrt_IaParam1 = register_tiled_diag_field (diag_mod_name, 'Ia_DERIVwrt_IaParam1',(/id_lon,id_lat/), &
     time, 'Partial derivative of Ia w/r/t IaParam1', '(ignitions/km2)/IaParam1_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_Ia_DERIVwrt_IaParam2 = register_tiled_diag_field (diag_mod_name, 'Ia_DERIVwrt_IaParam2',(/id_lon,id_lat/), &
     time, 'Partial derivative of Ia w/r/t IaParam2', '(ignitions/km2)/IaParam2_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_fire_fn_theta_DERIVwrt_thetaE = register_tiled_diag_field (diag_mod_name, 'fire_fn_theta_DERIVwrt_thetaE',(/id_lon,id_lat/), &
     time, 'Partial derivative of fTheta w/r/t thetaE', 'fTheta_units/thetaE_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_fire_fn_agb_DERIVwrt_param1 = register_tiled_diag_field (diag_mod_name, 'fire_fn_agb_DERIVwrt_param1',(/id_lon,id_lat/), &
     time, 'Partial derivative of fAGB w/r/t AGBparam1', 'fAGB_units/AGBparam1_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_fire_fn_agb_DERIVwrt_param2 = register_tiled_diag_field (diag_mod_name, 'fire_fn_agb_DERIVwrt_param2',(/id_lon,id_lat/), &
     time, 'Partial derivative of fAGB w/r/t AGBparam2', 'fAGB_units/AGBparam2_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_popDsupp_NF_DERIVwrt_eps1 = register_tiled_diag_field (diag_mod_name, 'popDsupp_NF_DERIVwrt_eps1',(/id_lon,id_lat/), &
     time, 'Partial derivative of popDsupp_NF w/r/t eps1', 'popDsupp_NF_units/eps1_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_popDsupp_NF_DERIVwrt_eps2 = register_tiled_diag_field (diag_mod_name, 'popDsupp_NF_DERIVwrt_eps2',(/id_lon,id_lat/), &
     time, 'Partial derivative of popDsupp_NF w/r/t eps2', 'popDsupp_NF_units/eps2_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_popDsupp_NF_DERIVwrt_eps3 = register_tiled_diag_field (diag_mod_name, 'popDsupp_NF_DERIVwrt_eps3',(/id_lon,id_lat/), &
     time, 'Partial derivative of popDsupp_NF w/r/t eps3', 'popDsupp_NF_units/eps3_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_fire_fn_rh_DERIVwrt_param1 = register_tiled_diag_field (diag_mod_name, 'fire_fn_rh_DERIVwrt_param1',(/id_lon,id_lat/), &
     time, 'Partial derivative of fRH w/r/t RHparam1', 'fRH_units/RHparam1_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_fire_fn_rh_DERIVwrt_param2 = register_tiled_diag_field (diag_mod_name, 'fire_fn_rh_DERIVwrt_param2',(/id_lon,id_lat/), &
     time, 'Partial derivative of fRH w/r/t RHparam2', 'fRH_units/RHparam2_unit', &
     missing_value=-1.0e+20, op=OP_SUM)
  id_BAperFire0_DERIVwrt_fireDur = register_tiled_diag_field (diag_mod_name, 'BAperFire0_DERIVwrt_fireDur',(/id_lon,id_lat/), &
       time, 'Partial derivative of BAperFire_0 w/r/t fire duration', 'm2/s', &
       missing_value=-1.0e+20, op=OP_SUM)
     
end subroutine vegn_fire_init


! ==============================================================================
subroutine vegn_fire_end()
  if (allocated(population_in)) deallocate(population_in)
  if (allocated(GDPpc_billion_in))        deallocate(GDPpc_billion_in)
  if (allocated(lightning_in))  deallocate(lightning_in)
  if (allocated(lightning_in_v2))  deallocate(lightning_in_v2)   ! SSR20151124
  if (allocated(Fc_in))         deallocate(Fc_in)
  if (allocated(Fp_in))         deallocate(Fp_in)
  if (allocated(crop_burn_rate_in)) deallocate(crop_burn_rate_in)
  if (allocated(past_burn_rate_in)) deallocate(past_burn_rate_in)
end subroutine vegn_fire_end


! ==============================================================================
! reads external data for the fire model
subroutine update_fire_data(time,is,ie,js,je)
  type(time_type), intent(in) :: time
  logical mask(size(GDPpc_billion_in,1),size(GDPpc_billion_in,2))
  integer :: year,month,day,hour,minute,second
  integer, intent(in) :: is, ie, js, je
  integer :: i, j

  if (fire_option /= FIRE_UNPACKED) return ! we don't need do do anything
  
  ! read_external_ts interpolates data conservatively in space (see init)
  ! and linearly in time.
  ! SSR: Really just does time interpolation
  if (use_FpopD_nf .OR. use_FpopD_ba .OR. Ia_alpha_monthly>0.0) then
     call read_external_ts(population_ts,time,population_in,mask=mask)
     where (.not.mask) population_in = 0.0 ! fill points where there is no data
  else
     population_in = 0.0
  endif
  if (use_Fgdp_nf .OR. use_Fgdp_ba) then
     call read_external_ts(GDPpc_billion_ts,time,GDPpc_billion_in,mask=mask)
     where (.not.mask) GDPpc_billion_in = 0.0 ! fill points where there is no data
  else
     GDPpc_billion_in = 0.0
  endif
! SSR20151124  call read_external_ts(lightning_ts,time,lightning_in,mask=mask)
! SSR20151124  where (.not.mask) lightning_in = 0.0
  
  call get_date(time,year,month,day,hour,minute,second)
  Fc_in = crop_burn_rate_in(:,:,month)
  Fp_in = past_burn_rate_in(:,:,month)
  lightning_in = lightning_in_v2(:,:,month)
  
  ! SSR: Check lightning data. (Sergey's suggestion.)
  do i = lnd%is, lnd%ie
      do j = lnd%js, lnd%je
          call check_var_range(lightning_in(i,j), 0.0, 10.**37, 'update_fire_data', 'lightning', lnd%time, FATAL)
      end do
  end do
  
  if (is_watch_point()) then
     write(*,*) '#### checkpoint update_fire_data #####'
     if (use_Fgdp_nf .OR. use_Fgdp_ba) __DEBUG2__(minval(GDPpc_billion_in),maxval(GDPpc_billion_in))
     if (use_FpopD_nf .OR. use_FpopD_ba .OR. Ia_alpha_monthly>0.0) __DEBUG2__(minval(population_in),maxval(population_in))
     __DEBUG2__(minval(lightning_in),maxval(lightning_in))
     write(*,*) '######################################'
  endif

end subroutine update_fire_data


! ==============================================================================
subroutine update_fire_fast(vegn,soil,diag, &
                            Tca,q,p_surf,cplr2land_wind, &
                            i,j,tile_area, &
                            latitude)
    type(vegn_tile_type), intent(inout) :: vegn
    type(soil_tile_type), intent(in) :: soil
    type(diag_buff_type), intent(inout) :: diag
    real, intent(in) :: q
    real, intent(in) :: Tca   ! Kelvin
    real, intent(in) :: p_surf
    real, intent(in) :: cplr2land_wind   ! Sheffield's 10-m wind
    integer, intent(in)  :: i,j   ! Coordinates of current point, for fire data
    real, intent(in)     :: tile_area   ! Area of tile (m2)
    real, intent(in)     :: latitude
    real   ::   lightning  ! Lightning flash density (flashes/km2/day)
    real   ::   popD       ! Population density (people/km2)
    real   ::   GDPpc ! GDP per capita ($/person)
    real   ::   fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb
    real   ::   qsat, rh
    real   ::   percentile = 0.95
    real   ::   depth_ave, theta
    real   ::   wind_forFire
    real   ::   ROS, LB, HB, gW, BAperFire_0
    real   ::   ROSmax, fire_dur, C_beta   ! SSR20151009
    integer::   vegn_cohort_1_species    
        
    call qscomp(Tca, p_surf, qsat)   ! (For RH) Calculates qsat, which is saturation specific humidity
    rh = q / qsat
    
    if (fire_windType == FIRE_WIND_CANTOP) then
       wind_forFire = vegn%wind_canaTop
    elseif (fire_windType == FIRE_WIND_10MTMP) then
       wind_forFire = vegn%wind_10mTmp
    elseif (fire_windType == FIRE_WIND_10MSHEFFIELD) then
       wind_forFire = cplr2land_wind
    endif
    
    if (depth_for_theta == -1.0) then
       depth_ave = -log(1.-percentile)*vegn%cohorts(1)%root_zeta   ! (For theta) Note that this is for the one-cohort version of the model.
    elseif (depth_for_theta > 0.0) then
       depth_ave = depth_for_theta
    endif
    if (theta_include_ice) then
       theta = soil_ave_theta2(soil,depth_ave)
    else
       theta = soil_ave_theta1(soil,depth_ave)
    endif
    
    if (use_Ftheta) then
       call vegn_fire_fn_theta(theta,fire_fn_theta)
    else
       fire_fn_theta = 1.0
    endif
    vegn%fireFtheta_av = vegn%fireFtheta_av + fire_fn_theta   ! SSR20150903
    
    if (use_Frh) then
       call vegn_fire_fn_rh(rh,fire_fn_rh)
    else
       fire_fn_rh = 1.0
    endif
    if (use_Ftemp) then
       call vegn_fire_fn_Tca(Tca,fire_fn_Tca)
    else
       fire_fn_Tca = 1.0
    endif
    call vegn_fire_fn_agb(vegn,soil,fire_fn_agb)
!    call vegn_fire_ROS(vegn,fire_fn_rh,theta,wind_forFire,ROS,LB,HB,gW)
!    call vegn_fire_ROS(vegn,fire_fn_rh,theta,wind_forFire,ROS,LB,HB,gW,ROSmax,C_beta)   ! SSR20151009
    call vegn_fire_ROS(vegn,fire_fn_rh,theta,fire_fn_theta,wind_forFire,ROS,LB,HB,gW,ROSmax,C_beta)   ! SSR20151216
    
    vegn_cohort_1_species = vegn%cohorts(1)%species
    call vegn_fire_BAperFire_noAnthro(ROS,LB,HB,vegn_cohort_1_species, &
                                      BAperFire_0, &
                                      fire_dur)   ! SSR20151009
    
    ! Could speed things up by only changing these monthly (or even yearly, for
    ! popD and GDPpc).
    lightning = lightning_in(i,j)
    if (use_FpopD_nf .OR. use_FpopD_ba .OR. Ia_alpha_monthly>0.0) popD =  population_in(i,j)
    if (use_Fgdp_nf .OR. use_Fgdp_ba)   GDPpc = GDPpc_billion_in(i,j) * 1.e9   ! Converting to dollars/person.
    
    if (is_watch_point()) then
       write(*,*) '#### checkpoint update_fire_fast #####'
       if (use_Fgdp_nf .OR. use_Fgdp_ba)   __DEBUG4__(minval(GDPpc_billion_in),maxval(GDPpc_billion_in),GDPpc_billion_in(i,j),GDPpc)
       if (use_FpopD_nf .OR. use_FpopD_ba .OR. Ia_alpha_monthly>0.0) __DEBUG4__(minval(population_in),maxval(population_in),population_in(i,j),popD)
       __DEBUG3__(minval(lightning_in),maxval(lightning_in),lightning_in(i,j))
       write(*,*) '######################################'
    endif
    
    ! SSR20160128: If update_fire_fast is getting called, make sure that this
    ! tile has a valid value.
    if (do_fire_fragmentation) then
       call check_var_range(vegn%max_fire_size, max_fire_size_min, 10.**37, 'update_fire_fast', 'vegn%max_fire_size', lnd%time, FATAL)
    endif
    
    call update_Nfire_BA_fast(diag,i,j,tile_area, &
                              fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                              BAperFire_0, &
                              vegn_cohort_1_species, &
                              vegn%trop_code, &   ! SSR20150831
                              lightning, popD, GDPpc, &
                              vegn%max_fire_size, &   ! SSR20150727
                              vegn%burned_frac, & ! inout variable. All others just in.
                              latitude, &
                              ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta)   ! SSR20151009

    ! Make sure that vegn%burned_frac is within good range
    call check_var_range(vegn%burned_frac, 0.0, 1.0, 'update_fire_fast', 'vegn%burned_frac', lnd%time, FATAL)
    
    if (.NOT. minimal_fire_diagnostics) then
       ! Save state variables for diagnostics
       
       
       call send_tile_data(id_lightning,      lightning,           diag)
       call send_tile_data(id_population, popD,  diag)
!       if (id_population>0) call send_tile_data(id_population, popD,  diag)
       call send_tile_data(id_GDPpc,      GDPpc, diag)
!       if (id_GDPpc>0)      call send_tile_data(id_GDPpc,      GDPpc, diag)
    
       ! Save fire stuff for diagnostics
       call send_tile_data(id_fire_q,           q,                diag)
       call send_tile_data(id_fire_rh,          rh,               diag)
       call send_tile_data(id_fire_theta,       theta,            diag)
       call send_tile_data(id_fire_Tca,         Tca-273.15,       diag)
       call send_tile_data(id_fire_fn_theta,    fire_fn_theta,    diag)
       call send_tile_data(id_fire_fn_rh,       fire_fn_rh,       diag)
       call send_tile_data(id_fire_fn_Tca,      fire_fn_Tca,      diag)
       call send_tile_data(id_fire_fn_agb,      fire_fn_agb,      diag)
       call send_tile_data(id_fire_agb,         vegn%fire_agb,         diag)
       call send_tile_data(id_BAperFire_0, BAperFire_0,        diag)
       call send_tile_data(id_fire_wind_forFire, wind_forFire, diag)
       call send_tile_data(id_ROS,              ROS,              diag)
       call send_tile_data(id_fire_duration_ave, fire_duration_array(vegn_cohort_1_species), diag)
       call send_tile_data(id_LB,               LB,               diag)
       call send_tile_data(id_HB,               HB,               diag)
       call send_tile_data(id_gW,               gW,               diag)
       call send_tile_data(id_fire_depth_ave,   depth_ave, diag)
       call send_tile_data(id_tropType,    1.0*vegn%trop_code, diag)
    endif
    
    ! Print diagnostics if really small BAperFire_0
    if (0.<BAperFire_0 .AND. BAperFire_0<min_fire_size) then
       if (print_min_fire_violations) then
          write(*,*) '######## BAperFire_0<min_fire_size ########'
          write(*,*) 'min_fire_size', min_fire_size
          write(*,*) 'BAperFire_0', BAperFire_0
          write(*,*) 'depth_ave', depth_ave
          write(*,*) 'fire_fn_theta', fire_fn_theta
          write(*,*) 'theta', theta
          write(*,*) 'fire_fn_rh', fire_fn_rh
          write(*,*) 'rh', rh
          write(*,*) 'fire_fn_Tca', fire_fn_Tca
          write(*,*) 'Tca (degC)', Tca-273.15
          write(*,*) 'wind_forFire', wind_forFire
          write(*,*) 'ROS', ROS
          write(*,*) 'LB', LB
          write(*,*) 'HB', HB
       endif
    endif
            
end subroutine update_fire_fast


subroutine update_fire_agri(vegn,Time,tile_area_km2,BF_mth,BA_mth)
  type(vegn_tile_type), intent(inout) :: vegn
  type(time_type), intent(in)  :: Time
  real, intent(in) :: tile_area_km2
  real, intent(out):: BF_mth, BA_mth   ! Total burned area (km2) or burned fraction of tile
  
  ! Calculate burned fraction
  call vegn_fire_BA_agri(vegn,Time,tile_area_km2,BA_mth,BF_mth)
  
  ! accumulate burned fraction since last burn (SSR: after vegn_disturbance.F90)
  vegn%burned_frac = vegn%burned_frac + BF_mth
  
  ! vegn%burned_frac shouldn't be >1 after having restricted BF_mth in vegn_fire_BA_agri
  call check_var_range(vegn%burned_frac, 0.0, 1.0, 'update_fire_agri', 'vegn%burned_frac', lnd%time, FATAL)
  
end subroutine update_fire_agri


subroutine update_fire_Fk(vegn,diag,i,j)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  integer, intent(in)  :: i,j   ! Coordinates of current point, for fire data
  
  vegn%Fcrop = Fc_in(i,j)
  vegn%Fpast = Fp_in(i,j)
  
  call send_tile_data(id_Fcrop, vegn%Fcrop, diag)
  call send_tile_data(id_Fpast, vegn%Fpast, diag)
  
end subroutine update_fire_Fk


! ==============================================================================


subroutine vegn_fire_fn_theta(theta,fire_fn_theta)
    real, intent(in) :: theta
    real, intent(out) :: fire_fn_theta
        
    if (fire_option_fTheta==FIRE_THETA_LI2012) then
       if (.NOT. thetaE_already_squared) then
          fire_fn_theta = exp(-pi*((theta/theta_extinction)**2.))
       else
          fire_fn_theta = exp(-pi*((theta**2.)/theta_extinction))
       endif
    elseif (fire_option_fTheta==FIRE_THETA_LOGISTIC) then
       fire_fn_theta = 1. / (1. + exp(-1.*theta_psi2*(theta-theta_psi3)))
    elseif (fire_option_fTheta==FIRE_THETA_GOMPERTZ) then
       fire_fn_theta = exp(-theta_gom2*exp(-theta_gom3*theta))
    endif
!     if (do_calc_derivs) then
!        if (.NOT. thetaE_already_squared) then
!           fire_fn_theta_DERIVwrt_thetaE = 2.*pi*(theta**2.)*fire_fn_theta/(theta_extinction**3.)
!        else
!           fire_fn_theta_DERIVwrt_thetaE = fire_fn_theta*pi*(theta**2.)/(theta_extinction**2.)
!        endif
!        call check_var_range(fire_fn_theta_DERIVwrt_thetaE, -10.**37, 10.**37, 'vegn_fire_fn_theta', 'fire_fn_theta_DERIVwrt_thetaE', lnd%time, FATAL)
!     endif
    
    ! SSR20160203
    ! Avoid craziness resulting from super-tiny values associated with logistic
    ! or Gompertz functions
    if (fire_fn_theta<1e-9 .AND. fire_fn_theta>0.0) then
       fire_fn_theta = 1e-9
    endif
    
    ! SSR20160202; SSR20160203
    if (do_calc_derivs) then
       if (fire_option_fTheta==FIRE_THETA_LI2012) then
          if (theta_ROSeffect_asFnTheta) then
             if (.NOT. thetaE_already_squared) then
                fire_TOTALfn_theta_DERIVwrt_thetaE = (6.*pi*(theta**2.)*exp(-(3.*pi*(theta**2.))/(theta_extinction**2.)))/(theta_extinction**3.)
             else
                fire_TOTALfn_theta_DERIVwrt_thetaE = (3.*pi*(theta**2.)*exp(-(3.*pi*(theta**2.))/theta_extinction))/(theta_extinction**2.)
             endif
          else
             call error_mesg('vegn_fire_fn_theta',&
                             'Add code to get theta derivative for when theta_ROSeffect_asFnTheta==FALSE.', &
                             FATAL)
          endif
          call check_var_range(fire_TOTALfn_theta_DERIVwrt_thetaE, -10.**37, 10.**37, 'vegn_fire_fn_theta', 'fire_fn_theta_DERIVwrt_thetaE', lnd%time, FATAL)
       elseif (fire_option_fTheta==FIRE_THETA_LOGISTIC) then
          if (theta_ROSeffect_asFnTheta) then
             fire_TOTALfn_theta_DERIVwrt_param1 = (3.*exp(-theta_psi2*(theta - theta_psi3))*(theta - theta_psi3))/((exp(-theta_psi2*(theta - theta_psi3)) + 1)**4.)
             fire_TOTALfn_theta_DERIVwrt_param2 = -(3.*theta_psi2*exp(-theta_psi2*(theta - theta_psi3)))/((exp(-theta_psi2*(theta - theta_psi3)) + 1)**4.)
          else
             call error_mesg('vegn_fire_fn_theta',&
                             'Add code to get theta derivatives for when theta_ROSeffect_asFnTheta==FALSE.', &
                             FATAL)
          endif
       elseif (fire_option_fTheta==FIRE_THETA_GOMPERTZ) then
          if (theta_ROSeffect_asFnTheta) then
             fire_TOTALfn_theta_DERIVwrt_param1 =  -3.*exp(-theta_gom3*theta)*exp(-3.*theta_gom2*exp(-theta_gom3*theta))
             fire_TOTALfn_theta_DERIVwrt_param2 = 3.*theta_gom2*theta*exp(-theta_gom3*theta)*exp(-3.*theta_gom2*exp(-theta_gom3*theta))
          else
             call error_mesg('vegn_fire_fn_theta',&
                             'Add code to get theta derivatives for when theta_ROSeffect_asFnTheta==FALSE.', &
                             FATAL)
          endif
       endif
    endif
    
    if (is_watch_point()) then
       write(*,*) '######## checkpoint vegn_fire_fn_theta ########'
       write(*,*) 'pi', pi
       write(*,*) 'theta', theta
       if (.NOT. thetaE_already_squared) then
          write(*,*) 'theta_extinction', theta_extinction
       else
          write(*,*) 'theta_extinction', sqrt(theta_extinction)
       endif
       write(*,*) 'fire_fn_theta', fire_fn_theta
       if (do_calc_derivs) then
          write(*,*) 'fire_fn_theta_DERIVwrt_thetaE', fire_fn_theta_DERIVwrt_thetaE
       endif
       write(*,*) '#######################################################'
    endif
    
    call check_var_range(fire_fn_theta, 0.0, 1.0, 'vegn_fire_fn_theta', 'fire_fn_theta', lnd%time, FATAL)
    
end subroutine vegn_fire_fn_theta



subroutine vegn_fire_fn_rh(rh,fire_fn_rh)
    real, intent(in)  :: rh
    real, intent(out) :: fire_fn_rh
        
    if (fire_option_fRH==FIRE_RH_LI2012) then
       fire_fn_rh = max(0.,min(1., (rh_up-rh) / (rh_up-rh_lo)))
       if (do_calc_derivs) then
          if (rh <= rh_lo) then
             ! Actually not differentiable at rh==rh_lo! So it's fudged.
             fire_fn_rh_DERIVwrt_param1 = 0.
             fire_fn_rh_DERIVwrt_param2 = 0.
          elseif (rh > rh_lo .AND. rh < rh_up) then
!             fire_fn_rh_DERIVwrt_param1 = (rh_up - rh)/(rh_lo - rh_up)**2.   ! dFrh/d_rh_lo
!             fire_fn_rh_DERIVwrt_param2 = -1./(rh_lo - rh_up) - (rh_up - rh)/(rh_lo - rh_up)**2.   ! dFrh/d_rh_up
              ! SSR20160202
              if ((.NOT. C_beta_params_likeRH) .OR. theta_ROSeffect_asFnTheta) then
                 fire_TOTALfn_rh_DERIVwrt_param1 = -(3.*((rh - rh_up)**3.))/((rh_lo - rh_up)**4.)
                 fire_TOTALfn_rh_DERIVwrt_param2 = (3.*((rh - rh_up)**3.))/((rh_lo - rh_up)**4.) - (3.*((rh - rh_up)**2.))/((rh_lo - rh_up)**3.) ;
              else
                 call error_mesg('vegn_fire_fn_rh',&
                                 'Add code to get RH derivatives for when C_beta_params_likeRH==TRUE and theta_ROSeffect_asFnTheta==FALSE.', &
                                  FATAL)
              endif
          elseif (rh >= rh_up) then
             ! Actually not differentiable at rh==rh_up! So it's fudged.
             fire_fn_rh_DERIVwrt_param1 = 0.
             fire_fn_rh_DERIVwrt_param2 = 0.
          endif
       endif
    elseif (fire_option_fRH==FIRE_RH_LOGISTIC) then
       fire_fn_rh = 1. / (1. + exp(-1.*rh_psi2*(rh-rh_psi3)))
       if (do_calc_derivs) then
!          fire_fn_rh_DERIVwrt_param1 = (exp(-rh_psi2*(rh - rh_psi3))*(rh - rh_psi3)) &
!                                     /(exp(-rh_psi2*(rh - rh_psi3)) + 1.)**2.
!          fire_fn_rh_DERIVwrt_param2 = -(rh_psi2*exp(-rh_psi2*(rh - rh_psi3))) &
!                                      /(exp(-rh_psi2*(rh - rh_psi3)) + 1)**2.
          ! SSR20160202
          if ((.NOT. C_beta_params_likeRH) .OR. theta_ROSeffect_asFnTheta) then
             fire_TOTALfn_rh_DERIVwrt_param1 = (3.*exp(-rh_psi2*(rh - rh_psi3))*(rh - rh_psi3))/((exp(-rh_psi2*(rh - rh_psi3)) + 1)**4.)
             fire_TOTALfn_rh_DERIVwrt_param2 = -(3.*rh_psi2*exp(-rh_psi2*(rh - rh_psi3)))/((exp(-rh_psi2*(rh - rh_psi3)) + 1)**4.)
          else
             call error_mesg('vegn_fire_fn_rh',&
                             'Add code to get RH derivatives for when C_beta_params_likeRH==TRUE and theta_ROSeffect_asFnTheta==FALSE.', &
                             FATAL)
          endif
       endif
    elseif (fire_option_fRH==FIRE_RH_GOMPERTZ) then
       fire_fn_rh = exp(-rh_gom2*exp(-rh_gom3*rh))
       if (do_calc_derivs) then
!          fire_fn_rh_DERIVwrt_param1 = -exp(-rh_gom3*rh)*exp(-rh_gom2*exp(-rh_gom3*rh))
!          fire_fn_rh_DERIVwrt_param2 = rh_gom2*rh*exp(-rh_gom3*rh)*exp(-rh_gom2*exp(-rh_gom3*rh))
           ! SSR20160202
           if ((.NOT. C_beta_params_likeRH) .OR. theta_ROSeffect_asFnTheta) then
              fire_TOTALfn_rh_DERIVwrt_param1 = -3.*exp(-rh_gom3*rh)*exp(-3.*rh_gom2*exp(-rh_gom3*rh))
              fire_TOTALfn_rh_DERIVwrt_param2 = 3.*rh_gom2*rh*exp(-rh_gom3*rh)*exp(-3.*rh_gom2*exp(-rh_gom3*rh))
              ! SSR20160217: Kludgey fix for when rh_gom3 is extreme
              if (.NOT. fire_TOTALfn_rh_DERIVwrt_param1>-10e37 .AND. .NOT. fire_TOTALfn_rh_DERIVwrt_param1<10e37) then
                 fire_TOTALfn_rh_DERIVwrt_param1 = 0.0
              endif
              if (.NOT. fire_TOTALfn_rh_DERIVwrt_param2>-10e37 .AND. .NOT. fire_TOTALfn_rh_DERIVwrt_param2<10e37) then
                 fire_TOTALfn_rh_DERIVwrt_param2 = 0.0
              endif
           else
             call error_mesg('vegn_fire_fn_rh',&
                             'Add code to get RH derivatives for when C_beta_params_likeRH==TRUE and theta_ROSeffect_asFnTheta==FALSE.', &
                             FATAL)
           endif
       endif
    endif
    
    ! Avoid craziness resulting from super-tiny values associated with logistic
    ! or Gompertz functions
    if (fire_fn_rh<1e-9 .AND. fire_fn_rh>0.0) then
       fire_fn_rh = 1e-9
    endif
    
    if (is_watch_point()) then
          write(*,*) '######## checkpoint vegn_fire_fn_rh ########'
          write(*,*) 'rh', rh
          write(*,*) 'fire_fn_rh', fire_fn_rh
          write(*,*) 'f_rh_style', f_rh_style
          if (fire_option_fRH==FIRE_RH_LI2012) then
             write(*,*) 'rh_up', rh_up
             write(*,*) 'rh_lo', rh_lo
          elseif (fire_option_fRH==FIRE_RH_LOGISTIC) then
             write(*,*) 'rh_psi2', rh_psi2
             write(*,*) 'rh_psi3', rh_psi3
          elseif (fire_option_fRH==FIRE_RH_GOMPERTZ) then
             write(*,*) 'rh_gom2', rh_gom2
             write(*,*) 'rh_gom3', rh_gom3
          endif
          if (do_calc_derivs) then
             write(*,*) 'fire_fn_rh_DERIVwrt_param1', fire_fn_rh_DERIVwrt_param1
             write(*,*) 'fire_fn_rh_DERIVwrt_param2', fire_fn_rh_DERIVwrt_param2
          endif
          write(*,*) '#######################################################'
    endif
        
    call check_var_range(fire_fn_rh, 0.0, 1.0, 'vegn_fire_fn_rh', 'fire_fn_rh', lnd%time, FATAL)
    if (do_calc_derivs) then
       call check_var_range(fire_fn_rh_DERIVwrt_param1, -10.**37, 10.**37, 'vegn_fire_fn_rh', 'fire_fn_rh_DERIVwrt_param1', lnd%time, FATAL)
       call check_var_range(fire_fn_rh_DERIVwrt_param2, -10.**37, 10.**37, 'vegn_fire_fn_rh', 'fire_fn_rh_DERIVwrt_param2', lnd%time, FATAL)
       call check_var_range(fire_TOTALfn_rh_DERIVwrt_param1, -10.**37, 10.**37, 'vegn_fire_fn_rh', 'fire_TOTALfn_rh_DERIVwrt_param1', lnd%time, FATAL)
       call check_var_range(fire_TOTALfn_rh_DERIVwrt_param2, -10.**37, 10.**37, 'vegn_fire_fn_rh', 'fire_TOTALfn_rh_DERIVwrt_param2', lnd%time, FATAL)
    endif
    
end subroutine vegn_fire_fn_rh



subroutine vegn_fire_fn_Tca(Tca,fire_fn_Tca)
    real, intent(in) :: Tca
    real, intent(out) :: fire_fn_Tca
    
    real :: Tca_celsius
    
    Tca_celsius = Tca - 273.15
    fire_fn_Tca = max(0., min(1., (Tca_celsius - T_lo_celsius)/(T_up_celsius-T_lo_celsius)))
    
    if (is_watch_point()) then
       write(*,*) '######## checkpoint vegn_fire_fn_Tca ########'
       write(*,*) 'Tca_celsius', Tca_celsius
       write(*,*) 'T_lo_celsius', T_lo_celsius
       write(*,*) 'T_up_celsius', T_up_celsius
       write(*,*) 'fire_fn_Tca', fire_fn_Tca
       write(*,*) '#######################################################'
    endif
        
    call check_var_range(fire_fn_Tca, 0.0, 1.0, 'vegn_fire_fn_Tca', 'fire_fn_Tca', lnd%time, FATAL)
    
end subroutine vegn_fire_fn_Tca



subroutine vegn_fire_fn_agb(vegn,soil,fire_fn_agb)
    type(vegn_tile_type), intent(inout) :: vegn
    type(soil_tile_type), intent(in) :: soil
    real, intent(out) :: fire_fn_agb
    type(vegn_cohort_type), pointer :: cc    ! current cohort
    
   ! vegn%fire_agb is updated in update_land_model_fast_0d to ensure that it's done for 
   ! all tiles, since fire_agb is what's used to determine whether tiles can merge when 
   ! using new fire model. 
    
    if (fire_biomass_threshold==-1) then
       if (fire_option_fAGB==FIRE_AGB_LI2012) then
          fire_fn_agb = max(0.,min(1., (vegn%fire_agb-agb_lo) / (agb_up-agb_lo)))
          if (do_calc_derivs) then
             if (vegn%fire_agb <= agb_lo) then
                ! Actually not differentiable at agb==agb_lo! So it's fudged.
                fire_fn_agb_DERIVwrt_param1 = 0.
                fire_fn_agb_DERIVwrt_param2 = 0.
             elseif (vegn%fire_agb > agb_lo .AND. vegn%fire_agb < agb_up) then
                fire_fn_agb_DERIVwrt_param1 = 1./(agb_lo - agb_up) - (agb_lo - vegn%fire_agb)/(agb_lo - agb_up)**2.   ! dFagb/d_agb_lo
                fire_fn_agb_DERIVwrt_param2 = (agb_lo - vegn%fire_agb)/(agb_lo - agb_up)**2.   ! dFagb/d_agb_up
             elseif (vegn%fire_agb >= agb_up) then
                ! Actually not differentiable at agb==agb_up! So it's fudged.
                fire_fn_agb_DERIVwrt_param1 = 0.
                fire_fn_agb_DERIVwrt_param2 = 0.
             endif
          endif
       elseif (fire_option_fAGB==FIRE_AGB_LOGISTIC) then
          fire_fn_agb = 1. / (1. + exp(-1.*agb_psi2*(vegn%fire_agb-agb_psi3)))
          if (do_calc_derivs) then
             fire_fn_agb_DERIVwrt_param1 = (exp(-agb_psi2*(vegn%fire_agb - agb_psi3))*(vegn%fire_agb - agb_psi3)) &
                                        /(exp(-agb_psi2*(vegn%fire_agb - agb_psi3)) + 1.)**2.
             fire_fn_agb_DERIVwrt_param2 = -(agb_psi2*exp(-agb_psi2*(vegn%fire_agb - agb_psi3))) &
                                         /(exp(-agb_psi2*(vegn%fire_agb - agb_psi3)) + 1)**2.
          endif
       elseif (fire_option_fAGB==FIRE_AGB_GOMPERTZ) then
          fire_fn_agb = exp(-agb_gom2*exp(-agb_gom3*vegn%fire_agb))
          if (do_calc_derivs) then
             fire_fn_agb_DERIVwrt_param1 = -exp(-agb_gom3*vegn%fire_agb)*exp(-agb_gom2*exp(-agb_gom3*vegn%fire_agb))   ! dFagb/d_agb_gom2
             fire_fn_agb_DERIVwrt_param2 = agb_gom2*vegn%fire_agb*exp(-agb_gom3*vegn%fire_agb)*exp(-agb_gom2*exp(-agb_gom3*vegn%fire_agb))   ! dFagb/d_agb_gom3
          endif
       endif
    elseif (fire_biomass_threshold>=0) then
       if (vegn%fire_agb > fire_biomass_threshold) then
          fire_fn_agb = 1.0
       else
          fire_fn_agb = 0.0
       endif
    endif
    
    ! Avoid craziness resulting from super-tiny values associated with logistic
    ! or Gompertz functions
    if (fire_fn_agb<1e-9 .AND. fire_fn_agb>0.0) then
       fire_fn_agb = 1e-9
    endif
    
    if (is_watch_point()) then
          write(*,*) '######## checkpoint vegn_fire_fn_agb ########'
          write(*,*) 'agb', vegn%fire_agb
          write(*,*) 'fire_fn_agb', fire_fn_agb
          write(*,*) 'f_agb_style', f_agb_style
          if (fire_option_fAGB==FIRE_AGB_LI2012) then
             write(*,*) 'agb_up', agb_up
             write(*,*) 'agb_lo', agb_lo
          elseif (fire_option_fAGB==FIRE_AGB_LOGISTIC) then
             write(*,*) 'agb_psi2', agb_psi2
             write(*,*) 'agb_psi3', agb_psi3
          elseif (fire_option_fAGB==FIRE_AGB_GOMPERTZ) then
             write(*,*) 'agb_gom2', agb_gom2
             write(*,*) 'agb_gom3', agb_gom3
          endif
          if (do_calc_derivs) then
             write(*,*) 'fire_fn_agb_DERIVwrt_param1', fire_fn_agb_DERIVwrt_param1
             write(*,*) 'fire_fn_agb_DERIVwrt_param2', fire_fn_agb_DERIVwrt_param2
          endif
          write(*,*) '#######################################################'
    endif
        
    call check_var_range(fire_fn_agb, 0.0, 1.0, 'vegn_fire_fn_agb', 'fire_fn_agb', lnd%time, FATAL)
    if (do_calc_derivs) then
       call check_var_range(fire_fn_agb_DERIVwrt_param1, -10.**37, 10.**37, 'vegn_fire_fn_agb', 'fire_fn_agb_DERIVwrt_param1', lnd%time, FATAL)
       call check_var_range(fire_fn_agb_DERIVwrt_param2, -10.**37, 10.**37, 'vegn_fire_fn_agb', 'fire_fn_agb_DERIVwrt_param2', lnd%time, FATAL)
    endif
    
end subroutine vegn_fire_fn_agb


subroutine vegn_fire_fn_popD(vegn_cohort_1_species,popD,fire_fn_popD_NF, fire_fn_popD_BA)
    integer, intent(in) ::  vegn_cohort_1_species
    real, intent(in)    ::  popD
    real, intent(out)   ::  fire_fn_popD_NF     ! Fraction of potential fires suppressed by popD
    real, intent(out)   ::  fire_fn_popD_BA     ! Fraction of potential BA realized via popD
        
    ! Number of fires
    if (use_FpopD_nf==.TRUE. .AND. popD>0.0) then
       fire_fn_popD_NF = popD_supp_eps1 - popD_supp_eps2*exp(-popD_supp_eps3*popD)
       if (do_calc_derivs) then
          popDsupp_NF_DERIVwrt_eps1 = 1.0
          popDsupp_NF_DERIVwrt_eps2 = -1.0 * exp(-popD_supp_eps3*popD)
          popDsupp_NF_DERIVwrt_eps3 = popD_supp_eps2 * popD * exp(-popD_supp_eps3*popD)
       endif
    else
       fire_fn_popD_NF = 0.0
       if (do_calc_derivs) then
          popDsupp_NF_DERIVwrt_eps1 = 0.0
          popDsupp_NF_DERIVwrt_eps2 = 0.0
          popDsupp_NF_DERIVwrt_eps3 = 0.0
       endif
    endif
    
    ! Burned area per fire
    if (use_FpopD_ba==.TRUE.) then
       if (popD <= 0.1) then
           fire_fn_popD_BA = 1.
       else if (vegn_cohort_1_species==SP_C4GRASS .OR. vegn_cohort_1_species==SP_C3GRASS) then
           fire_fn_popD_BA = 0.2 + 0.8*exp(-pi*sqrt(popD/450.))
       else
           fire_fn_popD_BA = 0.4 + 0.6*exp(-pi*popD/125.)
       endif
    else
       fire_fn_popD_BA = 1.0
    endif
    
    if (is_watch_point()) then
       write(*,*) '######## checkpoint vegn_fire_fn_popD ########'
       write(*,*) 'popD_supp_eps1', popD_supp_eps1
       write(*,*) 'popD_supp_eps2', popD_supp_eps2
       write(*,*) 'popD_supp_eps3', popD_supp_eps3
       write(*,*) 'popD', popD
       write(*,*) 'vegn_cohort_1_species', vegn_cohort_1_species
       write(*,*) 'fire_fn_popD_NF', fire_fn_popD_NF
       write(*,*) 'fire_fn_popD_BA', fire_fn_popD_BA
       if (do_calc_derivs) then
          write(*,*) 'popDsupp_NF_DERIVwrt_eps1', popDsupp_NF_DERIVwrt_eps1
          write(*,*) 'popDsupp_NF_DERIVwrt_eps2', popDsupp_NF_DERIVwrt_eps2
          write(*,*) 'popDsupp_NF_DERIVwrt_eps3', popDsupp_NF_DERIVwrt_eps3
       endif
       write(*,*) '#######################################################'
    endif
    
    call check_var_range(fire_fn_popD_NF, 0.0, 1.0, 'vegn_fire_fn_popD_NF', 'fire_fn_popD_NF', lnd%time, FATAL)
    call check_var_range(fire_fn_popD_BA, 0.0, 1.0, 'vegn_fire_fn_popD_BA', 'fire_fn_popD_BA', lnd%time, FATAL)
    
end subroutine vegn_fire_fn_popD


subroutine vegn_fire_fn_GDPpc(vegn_cohort_1_species,GDPpc,popD,fire_fn_GDPpc_NF,fire_fn_GDPpc_BA)
    ! Note that Li et al. (2013) don't apply this to "tropical closed forests," because
    ! they have a separate model for burning there. I may also want to change this so
    ! that the "tree-dominated" ecosystems are better delineated, e.g., through Olson
    ! classifications.
    
    integer, intent(in) :: vegn_cohort_1_species
    real, intent(in)    :: GDPpc   ! $/person
    real, intent(in)    :: popD    ! People/km2
    real, intent(out)   :: fire_fn_GDPpc_NF     ! Fraction of potential fires realized via GDPpc
    real, intent(out)   :: fire_fn_GDPpc_BA     ! Fraction of potential BA/fire realized via GDPpc  
    
    real                :: GDPpc_k   ! k$/person
    
    if (use_Fgdp_nf==.TRUE. .OR. use_Fgdp_ba==.TRUE.) then
       ! Convert GDPpc from $/person to k$/person
       GDPpc_k = GDPpc / 1000.
    endif
    
    ! Number of fires
    if (use_Fgdp_nf==.TRUE.) then
       if (popD <= 0.1) then
           fire_fn_GDPpc_NF = 1.
       else if (vegn_cohort_1_species==SP_C4GRASS .OR. vegn_cohort_1_species==SP_C3GRASS) then
           fire_fn_GDPpc_NF = 0.1 + 0.9*exp(-pi*sqrt(GDPpc_k/8.))
       else
           if (GDPpc_k < 20.) then
               fire_fn_GDPpc_NF = 1.
           else
               fire_fn_GDPpc_NF = 0.39
           endif
       end if
       if (is_watch_point() .AND. (fire_fn_GDPpc_NF<0. .OR. fire_fn_GDPpc_NF>1.)) then
           write(*,*) '######## BAD fire_fn_GDPpc_NF ########'
           write(*,*) 'GDPpc_k', GDPpc_k
       endif
    else
       fire_fn_GDPpc_NF = 1.0
    endif
    
    ! Burned area per fire
    if (use_Fgdp_ba==.TRUE.) then
       if (popD <= 0.1) then
           fire_fn_GDPpc_BA = 1.
       else if (vegn_cohort_1_species==SP_C4GRASS .OR. vegn_cohort_1_species==SP_C3GRASS) then
           fire_fn_GDPpc_BA = 0.2 + 0.8*exp(-pi*GDPpc_k/7.) ;
       else
           if (GDPpc_k <= 8.) then
               fire_fn_GDPpc_BA = 1.
           else if (GDPpc_k > 8. .AND. GDPpc_k <= 20.) then
               fire_fn_GDPpc_BA = 0.83
           else
               fire_fn_GDPpc_BA = 0.62
           end if
       end if
       if (is_watch_point() .AND. (fire_fn_GDPpc_BA<0. .OR. fire_fn_GDPpc_BA>1.)) then
           write(*,*) '######## BAD fire_fn_GDPpc_BA ########'
           write(*,*) 'GDPpc_k', GDPpc_k
       endif
    else
       fire_fn_GDPpc_BA = 1.0
    endif
    
    if (is_watch_point()) then
       write(*,*) '######## checkpoint vegn_fire_fn_GDPpc ########'
       write(*,*) 'GDPpc_k', GDPpc_k
       write(*,*) 'popD', popD
       write(*,*) 'vegn_cohort_1_species', vegn_cohort_1_species
       write(*,*) 'fire_fn_GDPpc_NF', fire_fn_GDPpc_NF
       write(*,*) 'fire_fn_GDPpc_BA', fire_fn_GDPpc_BA
       write(*,*) '#######################################################'
    endif
    
    call check_var_range(fire_fn_GDPpc_NF, 0.0, 1.0, 'vegn_fire_fn_GDPpc_NF', 'fire_fn_GDPpc_NF', lnd%time, FATAL)
    call check_var_range(fire_fn_GDPpc_BA, 0.0, 1.0, 'vegn_fire_fn_GDPpc_BA', 'fire_fn_GDPpc_BA', lnd%time, FATAL)
    
end subroutine vegn_fire_fn_GDPpc


subroutine vegn_fire_In(latitude,lightning,In)
    real, intent(in)    ::  latitude    ! Radians
    real, intent(in)    ::  lightning   ! Flashes/km2/day (cloud-cloud + cloud-ground)
    real, intent(out)   ::  In          ! Potential ignitions, natural (#/km2)

    real :: cloud2ground_frac
    
    cloud2ground_frac = 1. / (5.16 + 2.16*cos(latitude))
    In = lightning * cloud2ground_frac * In_c2g_ign_eff
    In = In * dt_fast/86400.

    if (is_watch_point()) then
       write(*,*) '#### checkpoint vegn_fire_In #####'
       write(*,*) 'latitude', latitude
       write(*,*) 'cloud2ground_frac', cloud2ground_frac
       write(*,*) 'lightning', lightning
       write(*,*) 'dt_fast', dt_fast
       write(*,*) 'In_c2g_ign_eff', In_c2g_ign_eff
       write(*,*) 'In', In
       write(*,*) '##################################'
    endif
    
    call check_var_range(lightning, 0.0, 10.**37, 'vegn_fire_In', 'lightning', lnd%time, FATAL)
    call check_var_range(cloud2ground_frac, 0.0, 1.0, 'vegn_fire_In', 'cloud2ground_frac', lnd%time, FATAL)
    call check_var_range(In, 0.0, lightning, 'vegn_fire_In', 'In', lnd%time, FATAL)
end subroutine vegn_fire_In


subroutine vegn_fire_Ia(popD,Ia)
    real, intent(in)    :: popD   ! People/km2
    real, intent(out)   :: Ia     ! Potential ignitions, anthropogenic (#/km2)
    
    if (popD > 0.0) then
       Ia = Ia_alpha_daily * Ia_param1 * (popD**(1.0-Ia_param2)) * dt_fast/86400.
    else
       Ia = 0.0   ! To avoid 0 in denominator when Ia_param2 > 1
    endif
    
    if (do_calc_derivs) then
       if (popD > 0.0) then
          Ia_DERIVwrt_alphaM = Ia_param1 * (popD**(1.0-Ia_param2)) * (12./days_per_year) * dt_fast/86400.
          Ia_DERIVwrt_IaParam1 = Ia_alpha_daily * (popD**(1.0-Ia_param2)) * dt_fast/86400.
          Ia_DERIVwrt_IaParam2 = -1.0 * Ia_alpha_daily * Ia_param1 * (popD**(1.0-Ia_param2)) * log(popD) * dt_fast/86400.
       else
          Ia_DERIVwrt_alphaM = 0.0   ! To avoid 0 in denominator when Ia_param2 > 1
          Ia_DERIVwrt_IaParam1 = 0.0   ! To avoid 0 in denominator when Ia_param2 > 1
          Ia_DERIVwrt_IaParam2 = 0.0   ! To avoid log(0)
       endif
    endif
    
    if (is_watch_point()) then
       write(*,*) '#### checkpoint vegn_fire_Ia #####'
       write(*,*) 'popD', popD
       write(*,*) 'Ia_alpha_daily', Ia_alpha_daily
       write(*,*) 'Ia_param1', Ia_param1
       write(*,*) 'Ia_param2', Ia_param2
       write(*,*) 'dt_fast', dt_fast
       write(*,*) 'Ia', Ia
       write(*,*) '##################################'
    endif
    call check_var_range(Ia, 0.0, 10.**37, 'vegn_fire_Ia', 'Ia', lnd%time, FATAL)
    
end subroutine vegn_fire_Ia


subroutine vegn_fire_Nfire(In, Ia, &
                               fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                               fire_fn_popD_NF, fire_fn_GDPpc_NF, &
                               Nfire_perKm2)
    real, intent(in)    ::  In, Ia, &
                            fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                            fire_fn_popD_NF, fire_fn_GDPpc_NF
    real, intent(out)   ::  Nfire_perKm2   ! Number of fires per km2
    
    Nfire_perKm2 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca * fire_fn_agb &
                      * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF
    
    if(is_watch_point()) then
        write(*,*) '######## checkpoint vegn_fire_NFire ########'
        write(*,*) 'In', In
        write(*,*) 'Ia', Ia
        write(*,*) 'fire_fn_theta', fire_fn_theta
        write(*,*) 'fire_fn_rh', fire_fn_rh
        write(*,*) 'fire_fn_Tca', fire_fn_Tca
        write(*,*) 'fire_fn_agb', fire_fn_agb
        write(*,*) 'fire_fn_popD_NF', fire_fn_popD_NF
        write(*,*) 'fire_fn_GDPpc_NF', fire_fn_GDPpc_NF
        write(*,*) 'Nfire_perKm2', Nfire_perKm2
        write(*,*) '############################################'
    endif
    
    call check_var_range(Nfire_perKm2, 0.0, 10.**37, 'vegn_fire_Nfire', 'Nfire_perKm2', lnd%time, FATAL)
end subroutine vegn_fire_Nfire


subroutine vegn_fire_ROS(vegn,fire_fn_rh,theta,&
                         fire_fn_theta,&   ! SSR20151216
                         wind_forFire,ROS,LB,HB,gW, &
                         ROS_max, C_beta)   ! SSR20151009
    type(vegn_tile_type), intent(inout) :: vegn
    real, intent(in)    :: fire_fn_rh
    real, intent(in)    :: theta    ! Calculated in subroutine vegn_fire_fn_theta
    real, intent(in)    :: fire_fn_theta    ! SSR20151216
    real                :: wind_forFire
    real, intent(out)   :: ROS      ! Downwind spread rate, m/s
    real, intent(out)   :: LB, HB   ! Length:breadth and head:back ratios
    real, intent(out)   :: gW
!    real                :: C_beta, C_m
    real                :: C_m   ! SSR20151009
    real, parameter     ::  vegn_min_height = 0.1
    real                ::  vegn_height_limited
    real    ::  g0
!    real    ::  ROS_max
    real, intent(out)    ::  ROS_max, C_beta   ! SSR20151009
    
    if (constant_wind_forFire >= 0.) then
       wind_forFire = constant_wind_forFire
    endif
    call check_var_range(wind_forFire, 0.0, 10.**37., 'vegn_fire_ROS', 'wind_forFire', lnd%time, FATAL)
    
    ! Calculate length:breadth and head:back ratios
    !!! LB_max and HB_max calculated in fire_init
    LB = LBratio_a + LBratio_b*(1. - exp(-LBratio_c*wind_forFire))
    call check_var_range(LB, LBratio_a, LB_max, 'vegn_fire_ROS', 'LB', lnd%time, FATAL)
    HB = (LB + sqrt(LB**2. - 1.)) / (LB - sqrt(LB**2. - 1.)) 
    
    ! Calculate wind multiplier effect
    g0 = (1. + 1./HB_max) / (2.*LB_max) ;
    if (use_Fwind) then
       gW = ((2.*LB) / (1. + 1./HB)) * g0 ;
    else
       gW = 1.0
    endif
    
    ! Choose maximum rate of spread, based on species
    ! SSR: At some point, should move this into vegn_data.F90 because these are species-level parameters
    if (vegn%cohorts(1)%species==SP_C4GRASS) then
        ROS_max = ROS_max_C4GRASS
    elseif (vegn%cohorts(1)%species==SP_C3GRASS) then
        ROS_max = ROS_max_C3GRASS
    elseif (vegn%cohorts(1)%species==SP_TEMPDEC) then
        ROS_max = ROS_max_TEMPDEC
    elseif (vegn%cohorts(1)%species==SP_TROPICAL) then
        if (vegn%trop_code == TROP_FOR) then
           ROS_max = ROS_max_TROPICAL
!        elseif (ROS_max_TROPSAV > 0.0) then
!           ROS_max = ROS_max_TROPSAV
!        else   ! If ROS_max_TROPSAV is nonpositive, then fall back on ROS_max_TROPICAL.
!           ROS_max = ROS_max_TROPICAL
!        endif
        ! SSR20160211
        elseif (vegn%trop_code == TROP_SHR) then
           if (ROS_max_TROPSHR > 0.0) then
              ROS_max = ROS_max_TROPSHR
           else
              ROS_max = ROS_max_TROPICAL
           endif
        elseif (vegn%trop_code == TROP_SAV) then
           if (ROS_max_TROPSAV > 0.0) then
              ROS_max = ROS_max_TROPSAV
           else
              ROS_max = ROS_max_TROPICAL
           endif
        endif
    elseif (vegn%cohorts(1)%species==SP_EVERGR) then
        ROS_max = ROS_max_EVERGR
    endif
    
    ! Calculate effects of RH and theta on fire ROS
    if (use_Cm) then
       if (theta_ROSeffect_asFnTheta) then
          C_beta = fire_fn_theta
       else
          if (.NOT. C_beta_params_likeRH) then
             C_beta = min(1.,max(0.,(C_beta_threshUP - theta) / (C_beta_threshUP - C_beta_threshLO)))
          else
             call vegn_fire_fn_rh(theta,C_beta)
          endif
       endif
       C_m = C_beta * fire_fn_rh
    else
       C_m = 1.0
    endif
    
    ! Calculate downwind rate of spread
    ROS = ROS_max * C_m * gW
    
    if (do_calc_derivs) then
       ROS_DERIVwrt_ROSmax = C_m * gW
    endif
    
    if(is_watch_point()) then
        write(*,*) '######## fire_checkpoint_3 ########'
        write(*,*) 'wind_forFire', wind_forFire
        write(*,*) 'LB', LB
        write(*,*) 'LB_max', LB_max
        write(*,*) 'HB', HB
        write(*,*) 'HB_max', HB_max
        write(*,*) 'g0', g0
        write(*,*) 'gW', gW
        write(*,*) 'C_beta', C_beta
        write(*,*) 'fire_fn_rh', fire_fn_rh
        write(*,*) 'C_m', C_m
        write(*,*) 'ROS', ROS
        if (do_calc_derivs) then
           write(*,*) 'ROS_DERIVwrt_ROSmax', ROS_DERIVwrt_ROSmax
        endif
    endif
    
    call check_var_range(C_m, 0.0, 1.0, 'vegn_fire_ROS', 'C_m', lnd%time, FATAL)
    call check_var_range(gW, 0.0, 1.0, 'vegn_fire_ROS', 'gW', lnd%time, FATAL)
    call check_var_range(ROS, 0.0, 10.**37, 'vegn_fire_ROS', 'ROS', lnd%time, FATAL)
    if (do_calc_derivs) then
       call check_var_range(ROS_DERIVwrt_ROSmax, -10.**37, 10.**37, 'vegn_fire_ROS', 'ROS_DERIVwrt_ROSmax', lnd%time, FATAL)
    endif
    
end subroutine vegn_fire_ROS


subroutine vegn_fire_BAperFire_noAnthro(ROS,LB,HB,vegn_cohort_1_species,BAperFire_0, &
                                        fire_dur)
    real, intent(in)    :: ROS,LB,HB
    integer, intent(in) :: vegn_cohort_1_species
    real, intent(out)   :: BAperFire_0   ! km2
    real, intent(out)   :: fire_dur   ! SSR20151009 (s)
    
    fire_dur = fire_duration_array(vegn_cohort_1_species)
    BAperFire_0 = (ROS**2.) &
                         * (fire_dur**2.) &
                         * ((1. + 1./HB)**2.) &
                         * pi &
                         / (4. * LB * (10.**6.))
                         
    if (do_calc_derivs) then
       BAperFire0_DERIVwrt_fireDur = BAperFire_0 &
                                     * (2.*fire_duration_array(vegn_cohort_1_species)) &
                                     / (fire_duration_array(vegn_cohort_1_species)**2.)
       if (ROS > 0.0) then
          BAperFire0_DERIVwrt_ROSmax = BAperFire_0 &
                                       * (2.*ROS*ROS_DERIVwrt_ROSmax) &
                                       / (ROS**2.)
       else
          BAperFire0_DERIVwrt_ROSmax = 0.0
       endif
    endif
                         
    if(is_watch_point()) then
       write(*,*) '### vegn_fire_BAperFire_0 ###'
       write(*,*) 'ROS', ROS
       write(*,*) 'LB ', LB
       write(*,*) 'HB ', HB
       write(*,*) 'vegn_cohort_1_species', vegn_cohort_1_species
       write(*,*) 'duration', fire_duration_array(vegn_cohort_1_species)
       write(*,*) 'BAperFire_0', BAperFire_0
       if (do_calc_derivs) then
          write(*,*) 'BAperFire0_DERIVwrt_fireDur', BAperFire0_DERIVwrt_fireDur
          write(*,*) 'BAperFire0_DERIVwrt_ROSmax', BAperFire0_DERIVwrt_ROSmax
       endif
       write(*,*) '####################################'
    endif
    
    call check_var_range(BAperFire_0, 0.0, 10.**37, 'vegn_fire_BAperFire_noAnthro', 'BAperFire_0', lnd%time, FATAL)
    if (do_calc_derivs) then
       call check_var_range(BAperFire0_DERIVwrt_fireDur, -10.**37, 10.**37, 'vegn_fire_BAperFire_noAnthro', 'BAperFire0_DERIVwrt_fireDur', lnd%time, FATAL)
       call check_var_range(BAperFire0_DERIVwrt_ROSmax, -10.**37, 10.**37, 'vegn_fire_BAperFire_noAnthro', 'BAperFire0_DERIVwrt_ROSmax', lnd%time, FATAL)
    endif
    
end subroutine vegn_fire_BAperFire_noAnthro


subroutine vegn_fire_BA_ntrlscnd(Nfire_perKm2,tile_area_km2, &
                                 max_fire_size, &   ! SSR20150727
                                 vegn_burned_frac, &
                                 BAperFire_0, BAperFire_1, BAperFire_2, BAperFire_3, &
                                 fire_fn_popD_BA,fire_fn_GDPpc_BA, &
                                 BA,BF,Nfire,BA_rate,BF_rate,Nfire_rate, &
                                 BA_reduction, vegn_unburned_area)
    real, intent(in)    ::  Nfire_perKm2, &
                            BAperFire_0, &
                            fire_fn_popD_BA, fire_fn_GDPpc_BA, &
                            tile_area_km2, &
                            max_fire_size   ! m2
    real, intent(inout) ::  vegn_burned_frac     ! Fraction of tile burned so far today
    real, intent(out)   ::  BAperFire_1, BAperFire_2, BAperFire_3, &
                            BA, BF, Nfire
    real, intent(out)   ::  BA_rate, BF_rate, Nfire_rate   ! Rates (per day)
    real, intent(out)   ::  BA_reduction, vegn_unburned_area   ! BA_reduction maybe not needed
    
    real   ::  vegn_unburned_frac, BA_tmp
    
    
    if (vegn_burned_frac < 1.0) then
    
       ! # fires, knowing that fire can't happen on land that already burned today
       vegn_unburned_frac = 1.0 - vegn_burned_frac
       vegn_unburned_area = tile_area_km2 * vegn_unburned_frac
       Nfire = Nfire_perKm2 * vegn_unburned_area   ! Number of fires per timestep
       
       ! Restrict this time step's burned fraction to amount of unburned area in tile
       BA_tmp = Nfire * BAperFire_0
       if (BA_tmp > vegn_unburned_area) then
          BA_reduction = (BA_tmp - vegn_unburned_area) / BA_tmp
          BAperFire_1 = vegn_unburned_area / Nfire
          BAperFire0_DERIVwrt_fireDur = 0.0
          if (is_watch_point()) then
             write(*,*) '######## BA > vegn_unburned_area. Reducing. ########'
             write(*,*) 'BA (before)', BA_tmp
             write(*,*) 'vegn_unburned_area', vegn_unburned_area
             write(*,*) 'BA_reduction', BA_reduction
          endif
       else
          BAperFire_1 = BAperFire_0
          BA_reduction = 0.0
       endif
       
       ! Restrict based on max fire size
       if (do_fire_fragmentation .AND. BAperFire_1 > max_fire_size*1e-6) then
          BAperFire_2 = max_fire_size*1e-6
          BA_reduction = (BAperFire_0 - max_fire_size*1e-6) / BAperFire_0
       else
          BAperFire_2 = BAperFire_1
       endif
    
       ! Incorporate human suppression
       BAperFire_3 = BAperFire_2 * fire_fn_popD_BA * fire_fn_GDPpc_BA
       
       ! Calculate BA
       BA = Nfire * BAperFire_3 * magic_scalar
       BF = BA / tile_area_km2
       
       ! Calculate per-day rates
       BF_rate =    BF    * (86400./dt_fast)
       Nfire_rate = Nfire * (86400./dt_fast)
       BA_rate =    BA    * (86400./dt_fast)
       
       if (0.<BA .AND. BA<min_combined_ba) then
          if (print_min_fire_violations .OR. is_watch_point()) then
             write(*,*) '######## BA<min_combined_ba ########'
             write(*,*) 'min_combined_ba', min_combined_ba
             write(*,*) 'BA', BA
             write(*,*) 'BF', BF
             write(*,*) 'Nfire', Nfire
          endif
          BF = 0.
          Nfire = 0.
          BA = 0.
          BAperFire_3 = 0.
       endif
       
       ! Update vegn%burned_frac. Done this way to avoid rounding error resulting in
       ! vegn%burned_frac > 1.0 by <1e-15.
       if (BA_reduction == 0.0) then
          vegn_burned_frac = vegn_burned_frac + BF
       else
          vegn_burned_frac = 1.0
       endif
       
    else ! (So vegn_burned_frac==1.0)
       BF = 0.0
       Nfire = 0.0
       BA = 0.0
       BF_rate = 0.0
       Nfire_rate = 0.0
       BA_rate = 0.0
       BA_reduction = 1.0
       BAperFire_1 = 0.0
       BAperFire_2 = 0.0
       BAperFire_3 = 0.0
       vegn_unburned_area = 0.0
    endif ! (vegn_burned_frac < 1.0)
    
    if(is_watch_point()) then
        write(*,*) '######## fire_checkpoint_4 ########'
        write(*,*) 'fire_fn_popD_BA', fire_fn_popD_BA
        write(*,*) 'fire_fn_GDPpc_BA', fire_fn_GDPpc_BA
        write(*,*) 'BAperFire_0', BAperFire_0
        write(*,*) 'BAperFire_1', BAperFire_1
        write(*,*) 'BAperFire_2', BAperFire_2
        write(*,*) 'BAperFire_3', BAperFire_3
        if (do_fire_fragmentation) then
           write(*,*) 'max_fire_size', max_fire_size
        endif
        write(*,*) 'tile_area_km2', tile_area_km2
        write(*,*) 'vegn_unburned_frac', vegn_unburned_frac
        write(*,*) 'vegn_unburned_area', vegn_unburned_area
        write(*,*) 'BF', BF
        write(*,*) 'BF_rate', BF_rate
        write(*,*) 'Nfire', Nfire
        write(*,*) 'Nfire_rate', Nfire_rate
        write(*,*) 'BA', BA
        write(*,*) 'BA_rate', BA_rate
    endif
    
    ! SSR20151208: Add some tolerance for BF at high end
    if (BF>1.0 .AND. BF<1.0+10.**-12) then
       BF = 1.0
    endif
    
    ! Error check
    call check_var_range(BAperFire_0, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'BAperFire_0', lnd%time, FATAL)
    call check_var_range(BAperFire_1, 0.0, BAperFire_0, 'vegn_fire_BA_ntrlscnd', 'BAperFire_1', lnd%time, FATAL)
    call check_var_range(BAperFire_2, 0.0, BAperFire_1, 'vegn_fire_BA_ntrlscnd', 'BAperFire_2', lnd%time, FATAL)
    call check_var_range(BAperFire_3, 0.0, BAperFire_2, 'vegn_fire_BA_ntrlscnd', 'BAperFire_3', lnd%time, FATAL)
    call check_var_range(BF, 0.0, 1.0, 'vegn_fire_BA_ntrlscnd', 'BF', lnd%time, FATAL)
    call check_var_range(Nfire, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'Nfire', lnd%time, FATAL)
    call check_var_range(BA, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'BA', lnd%time, FATAL)
    call check_var_range(BF_rate, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'BF_rate', lnd%time, FATAL)
    call check_var_range(Nfire_rate, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'Nfire_rate', lnd%time, FATAL)
    call check_var_range(BA_rate, 0.0, 10.**37, 'vegn_fire_BA_ntrlscnd', 'BA_rate', lnd%time, FATAL)

end subroutine vegn_fire_BA_ntrlscnd


subroutine vegn_burn(vegn,soil,tile_area_m2,force_watch_point)
    type(vegn_tile_type), intent(inout) :: vegn
    type(soil_tile_type), intent(inout) :: soil
    type(vegn_cohort_type), pointer :: cc    ! current cohort
    
    real, intent(in) :: tile_area_m2   ! Area of land in tile, m2
    logical, intent(in) :: force_watch_point
    
    real    ::    CC_leaf, CC_stem, CC_litter, CC_living   ! What fraction of tissues get combusted?
    real    ::    CC_litter_inclBF   ! SSR20151217
    real    ::    burned_total, killed_total, &   ! Sum of all burned (i.e., combusted), fire-killed biomass
                  burned_total_noCWL   ! SSR20151227
    real    ::    burned_bl_cc, burned_bwood_cc, burned_bsw_cc, burned_blv_cc   ! How much of each tissue got combusted in the current cohort?
    real    ::    burned_bl,    burned_bwood,    burned_bsw,    burned_blv      ! How much of each tissue got combusted across all cohorts?
    real    ::    burned_leafLitter, burned_fineWoodLitter, burned_coarseWoodLitter   ! How much of litter gets combusted across all litter cohorts?
    real    ::    burnedLitter_1_lf(n_c_types), burnedLitter_2_lf(n_c_types), burnedLitter_3_lf   ! For leaf litter, how much gets burned in: "Litter," protected, and microbial C?
    real    ::    burnedLitter_1_fw(n_c_types), burnedLitter_2_fw(n_c_types), burnedLitter_3_fw   ! For fine wood litter, how much gets burned in: "Litter," protected, and microbial C?
    real    ::    burnedLitter_1_cw(n_c_types), burnedLitter_2_cw(n_c_types), burnedLitter_3_cw   ! For coarse wood, how much gets burned in: "Litter," protected, and microbial C?
    real    ::    fireMort_leaf, fireMort_stem, fireMort_root, fireMort_living   ! What fraction of non-combusted tissues get killed by fire?
    real    ::    killed_bl_cc, killed_bwood_cc, killed_bsw_cc, killed_blv_cc, killed_br_cc   ! How much of non-combusted tissue got killed in the current cohort?
    real    ::    killed_bl,    killed_bwood,    killed_bsw,    killed_blv,    killed_br      ! How much of non-combusted tissue got killed across all cohorts?
    real    ::    killed_wood, killed_live   ! From standing biomass, partitions wood from leaves/roots/virtual pools
    real    ::    new_fast_C_ag, new_slow_C_ag, new_fast_C_bg, new_slow_C_bg
    real    ::    burned_living_cc, burned_Cgain_cc, burned_Wgain_cc, burned_Cgain, burned_Wgain, &
                  killed_living_cc, killed_Cgain_cc, killed_Wgain_cc, killed_Cgain, killed_Wgain
    real    ::    bliving_0, bwood_0, npp_prevDay_0, npp_prevDayTmp_0, &
                  bliving_1, bwood_1, npp_prevDay_1, npp_prevDayTmp_1
    real    ::    fire_agb_0, fire_agb_1
    integer ::    i
    real    ::    frac_nppPrevDay_toRemove
    real    ::    burned_frac
    real    ::    bl_inBurn_cc, bwood_inBurn_cc, bsw_inBurn_cc, blv_inBurn_cc, &
                  br_inBurn_cc, bliving_inBurn_cc, Cgain_inBurn_cc, Wgain_inBurn_cc, &
                  npp_prevDay_inBurn_cc, npp_prevDayTmp_inBurn_cc
    real    ::    bl0, blv0, bwood0, bsw0, br0, Cgain0, Wgain0, &
                  bl1, blv1, bwood1, bsw1, br1, Cgain1, Wgain1
    real    ::    seconds_per_year, days_per_year
    
    ! SSR20151217
    real    ::    lfl_unpr_fast, lfl_unpr_slow, lfl_unpr_dmic, &
                  lfl_prot_fast, lfl_prot_slow, lfl_prot_dmic, &
                  cwl_unpr_fast, cwl_unpr_slow, cwl_unpr_dmic, &
                  cwl_prot_fast, cwl_prot_slow, cwl_prot_dmic
                  
    
    if (vegn%cohorts(1)%species==SP_C4GRASS .OR. vegn%cohorts(1)%species==SP_C3GRASS) then
        ! Use Li2012 values for C4 = C3 Non-Arctic = C3 Arctic
        CC_leaf = CC_GRASS_leaf
        CC_stem = CC_GRASS_stem
        CC_litter = CC_GRASS_litter
        fireMort_leaf = fireMort_GRASS_leaf
        fireMort_stem = fireMort_GRASS_stem
        fireMort_root = fireMort_GRASS_root
    else if (vegn%cohorts(1)%species==SP_TEMPDEC) then
        ! Use Li2012 values for Broadleaf Deciduous Trees, temperate
        CC_leaf = CC_TEMPDEC_leaf
        CC_stem = CC_TEMPDEC_stem
        CC_litter = CC_TEMPDEC_litter
        fireMort_leaf = fireMort_TEMPDEC_leaf
        fireMort_stem = fireMort_TEMPDEC_stem
        fireMort_root = fireMort_TEMPDEC_root
    else if (vegn%cohorts(1)%species==SP_TROPICAL) then
!        ! Use Li2012 values for Broadleaf Evergreen Trees, tropical
!        CC_leaf = 0.70
!        CC_stem = 0.15
!        CC_litter = 0.50
!        fireMort_leaf = fireMort_TROPICAL_leaf
!        fireMort_stem = fireMort_TROPICAL_stem
!        fireMort_root = fireMort_TROPICAL_root
       if (vegn%trop_code == TROP_SHR .OR. vegn%trop_code == TROP_SAV) then
          ! Use Li2012 values for Broadleaf Deciduous Trees, tropical
          CC_leaf = CC_TROPSHRSAV_leaf
          CC_stem = CC_TROPSHRSAV_stem
          CC_litter = CC_TROPSHRSAV_litter
          fireMort_leaf = fireMort_TROPSHRSAV_leaf
          fireMort_stem = fireMort_TROPSHRSAV_stem
          fireMort_root = fireMort_TROPSHRSAV_root
       elseif (vegn%trop_code == TROP_FOR) then
          ! Use Li2012 values for Broadleaf Evergreen Trees, tropical
          CC_leaf = CC_TROPICAL_leaf
          CC_stem = CC_TROPICAL_stem
          CC_litter = CC_TROPICAL_litter
          fireMort_leaf = fireMort_TROPICAL_leaf
          fireMort_stem = fireMort_TROPICAL_stem
          fireMort_root = fireMort_TROPICAL_root
       endif
    else if (vegn%cohorts(1)%species==SP_EVERGR) then
        ! Use Li2012 values for Needleleaf Evergreen Trees, boreal
        CC_leaf = CC_EVERGR_leaf
        CC_stem = CC_EVERGR_stem
        CC_litter = CC_EVERGR_litter
        fireMort_leaf = fireMort_EVERGR_leaf
        fireMort_stem = fireMort_EVERGR_stem
        fireMort_root = fireMort_EVERGR_root
    end if
    
    if (do_fire_tiling==.FALSE. &
    .OR. vegn%landuse==LU_CROP &
    .OR. (vegn%landuse==LU_PAST .AND..NOT. split_past_tiles) &   ! SSR20151118
    .OR. vegn%burned_frac*tile_area_m2*1e-6 < min_BA_to_split) then
       burned_frac = vegn%burned_frac
    else
       burned_frac = 1.0
    endif
    
    ! Burn litter
    if (soil_carbon_option==SOILC_CORPSE) then
       
       ! SSR20151217
       call poolTotalCarbon(soil%leafLitter, &
                            fastC=lfl_unpr_fast,slowC=lfl_unpr_slow,deadMicrobeC=lfl_unpr_dmic, &
                            fast_protectedC=lfl_prot_fast,slow_protectedC=lfl_prot_slow,deadmic_protectedC=lfl_prot_dmic)
       call poolTotalCarbon(soil%coarseWoodLitter, &
                            fastC=cwl_unpr_fast,slowC=cwl_unpr_slow,deadMicrobeC=cwl_unpr_dmic, &
                            fast_protectedC=cwl_prot_fast,slow_protectedC=cwl_prot_slow,deadmic_protectedC=cwl_prot_dmic)
       
       CC_litter_inclBF = CC_litter * burned_frac
       call remove_carbon_fraction_from_pool(soil%leafLitter,CC_litter_inclBF,litter_removed=burnedLitter_1_lf, &
                                                                              protected_removed=burnedLitter_2_lf, &
                                                                              liveMicrobe_removed=burnedLitter_3_lf)
       burned_leafLitter = sum(burnedLitter_1_lf) + sum(burnedLitter_2_lf) + burnedLitter_3_lf
       call remove_carbon_fraction_from_pool(soil%fineWoodLitter,CC_litter_inclBF,litter_removed=burnedLitter_1_fw, &
                                                                                  protected_removed=burnedLitter_2_fw, &
                                                                                  liveMicrobe_removed=burnedLitter_3_fw)
       burned_fineWoodLitter = sum(burnedLitter_1_fw) + sum(burnedLitter_2_fw) + burnedLitter_3_fw
       if (burn_cwl_like_stem) then
          call remove_carbon_fraction_from_pool(soil%coarseWoodLitter,CC_stem*burned_frac,litter_removed=burnedLitter_1_cw, &
                                                                                       protected_removed=burnedLitter_2_cw, &
                                                                                       liveMicrobe_removed=burnedLitter_3_cw)
       else
          call remove_carbon_fraction_from_pool(soil%coarseWoodLitter,CC_litter_inclBF,litter_removed=burnedLitter_1_cw, &
                                                                                       protected_removed=burnedLitter_2_cw, &
                                                                                       liveMicrobe_removed=burnedLitter_3_cw)
       endif
       burned_coarseWoodLitter = sum(burnedLitter_1_cw) + sum(burnedLitter_2_cw) + burnedLitter_3_cw
       
       ! SSR20151217
       if (lfl_unpr_fast > 0.0) soil%ssr_lfl_unpr_fast_Ko_accum = soil%ssr_lfl_unpr_fast_Ko_accum + CC_litter_inclBF
       if (lfl_unpr_slow > 0.0) soil%ssr_lfl_unpr_slow_Ko_accum = soil%ssr_lfl_unpr_slow_Ko_accum + CC_litter_inclBF
       if (lfl_unpr_dmic > 0.0) soil%ssr_lfl_unpr_dmic_Ko_accum = soil%ssr_lfl_unpr_dmic_Ko_accum + CC_litter_inclBF
       if (lfl_prot_fast > 0.0) soil%ssr_lfl_prot_fast_Ko_accum = soil%ssr_lfl_prot_fast_Ko_accum + CC_litter_inclBF
       if (lfl_prot_slow > 0.0) soil%ssr_lfl_prot_slow_Ko_accum = soil%ssr_lfl_prot_slow_Ko_accum + CC_litter_inclBF
       if (lfl_prot_dmic > 0.0) soil%ssr_lfl_prot_dmic_Ko_accum = soil%ssr_lfl_prot_dmic_Ko_accum + CC_litter_inclBF
       if (cwl_unpr_fast > 0.0) soil%ssr_cwl_unpr_fast_Ko_accum = soil%ssr_cwl_unpr_fast_Ko_accum + CC_litter_inclBF
       if (cwl_unpr_slow > 0.0) soil%ssr_cwl_unpr_slow_Ko_accum = soil%ssr_cwl_unpr_slow_Ko_accum + CC_litter_inclBF
       if (cwl_unpr_dmic > 0.0) soil%ssr_cwl_unpr_dmic_Ko_accum = soil%ssr_cwl_unpr_dmic_Ko_accum + CC_litter_inclBF
       if (cwl_prot_fast > 0.0) soil%ssr_cwl_prot_fast_Ko_accum = soil%ssr_cwl_prot_fast_Ko_accum + CC_litter_inclBF
       if (cwl_prot_slow > 0.0) soil%ssr_cwl_prot_slow_Ko_accum = soil%ssr_cwl_prot_slow_Ko_accum + CC_litter_inclBF
       if (cwl_prot_dmic > 0.0) soil%ssr_cwl_prot_dmic_Ko_accum = soil%ssr_cwl_prot_dmic_Ko_accum + CC_litter_inclBF
    endif
    
    fire_agb_0 = vegn%fire_agb
    
    ! Burn/kill living vegetation (including "dead" wood)
    !!!!! Note that all species die from fire (i.e., this ignores spdata(cc%species)%mortality_kills_balive)
    burned_bl=0.0 ; burned_bwood=0.0 ; burned_bsw=0.0 ; burned_blv=0.0 ; burned_Cgain = 0.0; burned_Wgain = 0.0;
    killed_bl=0.0 ; killed_bwood=0.0 ; killed_bsw=0.0 ; killed_blv=0.0 ; killed_br=0.0 ; killed_Cgain = 0.0; killed_Wgain = 0.0;
    if(is_watch_point() .OR. force_watch_point .OR. (is_watch_cell() .AND. force_watch_cell_burning)) then
       bliving_0 = 0.0; bwood_0 = 0.0; npp_prevDay_0 = 0.0 ;
       bliving_1 = 0.0; bwood_1 = 0.0; npp_prevDay_1 = 0.0 ;
       bl0 = 0.0; blv0 = 0.0; bwood0 = 0.0; bsw0 = 0.0; br0 = 0.0;
       bl1 = 0.0; blv1 = 0.0; bwood1 = 0.0; bsw1 = 0.0; br1 = 0.0;
    endif
    do i = 1,vegn%n_cohorts
        cc => vegn%cohorts(i)
        
        ! Get initial biomass
        bl_inBurn_cc = cc%bl
        bwood_inBurn_cc = cc%bwood
        bsw_inBurn_cc = cc%bsw
        blv_inBurn_cc = cc%blv
        br_inBurn_cc = cc%br
        bliving_inBurn_cc = cc%bliving
        Cgain_inBurn_cc = cc%carbon_gain
        Wgain_inBurn_cc = cc%bwood_gain
        npp_prevDay_inBurn_cc = cc%npp_previous_day
        npp_prevDayTmp_inBurn_cc = cc%npp_previous_day_tmp
        if(is_watch_point() .OR. force_watch_point .OR. (is_watch_cell() .AND. force_watch_cell_burning)) then
           bl0 = bl0 + cc%bl
           blv0 = blv0 + cc%blv
           bwood0 = bwood0 + cc%bwood
           bsw0 = bsw0 + cc%bsw
           br0 = br0 + cc%br
           Cgain0 = Cgain0 + cc%carbon_gain
           Wgain0 = Wgain0 + cc%bwood_gain
           bliving_0 = bliving_0 + cc%bliving
           bwood_0 = bwood_0 + cc%bwood
           npp_prevDay_0 = npp_prevDay_0 + cc%npp_previous_day
           npp_prevDayTmp_0 = npp_prevDayTmp_0 + cc%npp_previous_day_tmp
        endif
        
        ! Combustion of leaves
        burned_bl_cc = CC_leaf * cc%bl
        burned_bl = burned_bl + burned_bl_cc
        
        ! Combustion of wood
        burned_bwood_cc = CC_stem * cc%bwood
        burned_bwood = burned_bwood + burned_bwood_cc
        
        ! Combustion of sapwood
        burned_bsw_cc = CC_stem * cc%bsw
        burned_bsw = burned_bsw + burned_bsw_cc
        
        ! Combustion of virtual leaf biomass
        burned_blv_cc = CC_stem * cc%blv
        burned_blv = burned_blv + burned_blv_cc
        
        ! Combustion of carbon_gain and bwood_gain
        burned_living_cc = burned_bl_cc + burned_bsw_cc + burned_blv_cc
        CC_living = burned_living_cc / bliving_inBurn_cc
        if (bliving_inBurn_cc==0.0) then
           CC_living = 1.0
        endif
        burned_Cgain_cc = Cgain_inBurn_cc * CC_living
        burned_Cgain = burned_Cgain + burned_Cgain_cc
        burned_Wgain_cc = Wgain_inBurn_cc * CC_living
        burned_Wgain = burned_Wgain + burned_WGain_cc
                
        ! Update biomass after combustion
        bl_inBurn_cc = bl_inBurn_cc - burned_bl_cc
        bwood_inBurn_cc = bwood_inBurn_cc - burned_bwood_cc
        bsw_inBurn_cc = bsw_inBurn_cc - burned_bsw_cc
        blv_inBurn_cc = blv_inBurn_cc - burned_blv_cc
        bliving_inBurn_cc = bsw_inBurn_cc + bl_inBurn_cc + br_inBurn_cc + blv_inBurn_cc;
        Cgain_inBurn_cc = Cgain_inBurn_cc - burned_Cgain_cc
        Wgain_inBurn_cc = Wgain_inBurn_cc - burned_Wgain_cc
        
        ! Mortality of remaining leaves
        killed_bl_cc = fireMort_leaf * bl_inBurn_cc
        killed_bl = killed_bl + killed_bl_cc
        
        ! Mortality of remaining wood
        killed_bwood_cc = fireMort_stem * bwood_inBurn_cc
        killed_bwood = killed_bwood + killed_bwood_cc
        
        ! Mortality of remaining sapwood
        killed_bsw_cc = fireMort_stem * bsw_inBurn_cc
        killed_bsw = killed_bsw + killed_bsw_cc
        
        ! Mortality of remaining virtual leaf biomass
        killed_blv_cc = fireMort_stem * blv_inBurn_cc
        killed_blv = killed_blv + killed_blv_cc
        
        ! Mortality of roots
        killed_br_cc = fireMort_root * br_inBurn_cc
        killed_br = killed_br + killed_br_cc
        
        ! Mortality of carbon_gain and bwood_gain
        killed_living_cc = killed_bl_cc + killed_bsw_cc + killed_blv_cc + killed_br_cc
        fireMort_living = killed_living_cc / bliving_inBurn_cc
        if (bliving_inBurn_cc==0.0) then
           fireMort_living = 1.0
        endif
        killed_Cgain_cc = Cgain_inBurn_cc * fireMort_living
        killed_Cgain = killed_Cgain + killed_Cgain_cc
        killed_Wgain_cc = Wgain_inBurn_cc * fireMort_living
        killed_Wgain = killed_Wgain + killed_Wgain_cc
        
        ! Update biomass after mortality
        bl_inBurn_cc = bl_inBurn_cc - killed_bl_cc
        bwood_inBurn_cc = bwood_inBurn_cc - killed_bwood_cc
        bsw_inBurn_cc = bsw_inBurn_cc - killed_bsw_cc
        blv_inBurn_cc = blv_inBurn_cc - killed_blv_cc
        br_inBurn_cc = br_inBurn_cc - killed_br_cc
        bliving_inBurn_cc = bsw_inBurn_cc + bl_inBurn_cc + br_inBurn_cc + blv_inBurn_cc;
        Cgain_inBurn_cc = Cgain_inBurn_cc - killed_Cgain_cc
        Wgain_inBurn_cc = Wgain_inBurn_cc - killed_Wgain_cc
        
        ! Update biomass after combustion AND mortality
        ! Weights based on values inside and outside burn perimeter.
        cc%bl =          (1.0-burned_frac)*cc%bl          + burned_frac*bl_inBurn_cc
        cc%bwood =       (1.0-burned_frac)*cc%bwood       + burned_frac*bwood_inBurn_cc
        cc%bsw =         (1.0-burned_frac)*cc%bsw         + burned_frac*bsw_inBurn_cc
        cc%blv =         (1.0-burned_frac)*cc%blv         + burned_frac*blv_inBurn_cc
        cc%br =          (1.0-burned_frac)*cc%br          + burned_frac*br_inBurn_cc
        cc%carbon_gain = (1.0-burned_frac)*cc%carbon_gain + burned_frac*Cgain_inBurn_cc
        cc%bwood_gain =  (1.0-burned_frac)*cc%bwood_gain  + burned_frac*Wgain_inBurn_cc
        cc%bliving =     (1.0-burned_frac)*cc%bliving     + burned_frac*(bl_inBurn_cc+blv_inBurn_cc+bsw_inBurn_cc+br_inBurn_cc)
        
        if (adj_nppPrevDay /= 0) then
           if (adj_nppPrevDay==1) then
              frac_nppPrevDay_toRemove = (burned_bl_cc+killed_bl_cc)/(cc%bl+burned_bl_cc+killed_bl_cc)
              if (cc%bl+burned_bl_cc+killed_bl_cc == 0) then
                 frac_nppPrevDay_toRemove = 1.0
              endif
           elseif (adj_nppPrevDay==2) then
              frac_nppPrevDay_toRemove = (burned_bsw_cc + burned_bl_cc + burned_blv_cc + killed_bsw_cc + killed_bl_cc + killed_br_cc + killed_blv_cc) / (cc%bliving + burned_bsw_cc + burned_bl_cc + burned_blv_cc + killed_bsw_cc + killed_bl_cc + killed_br_cc + killed_blv_cc)
              if (cc%bliving + burned_bsw_cc + burned_bl_cc + burned_blv_cc + killed_bsw_cc + killed_bl_cc + killed_br_cc + killed_blv_cc == 0) then
                 frac_nppPrevDay_toRemove = 1.0
              endif
           elseif (adj_nppPrevDay==3) then
              ! Elena's idea to just say "no NPP, because it burned." Kinda sketchy.
              frac_nppPrevDay_toRemove = 1.0
           else
              call error_mesg('vegn_burn','adj_nppPrevDay option is invalid.', FATAL)
           endif
           cc%npp_previous_day     = (1.0-burned_frac)*cc%npp_previous_day     + burned_frac*((1.0-frac_nppPrevDay_toRemove)*cc%npp_previous_day)
           cc%npp_previous_day_tmp = (1.0-burned_frac)*cc%npp_previous_day_tmp + burned_frac*((1.0-frac_nppPrevDay_toRemove)*cc%npp_previous_day_tmp)
        endif
        
        if(is_watch_point() .OR. force_watch_point .OR. (is_watch_cell() .AND. force_watch_cell_burning)) then
           bliving_1 = bliving_1 + cc%bliving
           bwood_1 = bwood_1 + cc%bwood
           npp_prevDay_1 = npp_prevDay_1 + cc%npp_previous_day
           npp_prevDayTmp_1 = npp_prevDayTmp_1 + cc%npp_previous_day_tmp
           bl1 = bl1 + cc%bl
           blv1 = blv1 + cc%blv
           bwood1 = bwood1 + cc%bwood
           bsw1 = bsw1 + cc%bsw
           br1 = br1 + cc%br
           Cgain1 = Cgain1 + cc%carbon_gain
           Wgain1 = Wgain1 + cc%bwood_gain
        endif
    end do
    
    if (burned_frac < 1.0) then
       burned_bl = burned_bl*burned_frac
       burned_bwood = burned_bwood*burned_frac
       burned_bsw = burned_bsw*burned_frac
       burned_blv = burned_blv*burned_frac
       burned_Cgain = burned_Cgain*burned_frac
       burned_Wgain = burned_Wgain*burned_frac
       killed_bl = killed_bl*burned_frac
       killed_bwood = killed_bwood*burned_frac
       killed_bsw = killed_bsw*burned_frac
       killed_blv = killed_blv*burned_frac
       killed_br = killed_br*burned_frac
       killed_Cgain = killed_Cgain*burned_frac
       killed_Wgain = killed_Wgain*burned_frac
    endif
    
    ! Get total combusted, killed
    burned_total = burned_bl + burned_bwood + burned_bsw + burned_blv + burned_Cgain + burned_Wgain
    if (soil_carbon_option==SOILC_CORPSE) then
       burned_total       = burned_total + burned_leafLitter + burned_fineWoodLitter + burned_coarseWoodLitter
       burned_total_noCWL = burned_total - burned_coarseWoodLitter
    else
       burned_total_noCWL = burned_total
    endif
    killed_wood = killed_bwood+killed_bsw   ! Should killed_bsw actually be in killed_live?
    killed_live = killed_bl+killed_br+killed_blv+killed_Cgain+killed_Wgain
    killed_total = killed_wood + killed_live
    ! Deal with tiny negative values
    if (burned_total<0 .AND. -1.0*burned_total<carbon_cons_tol/2.0) then
       burned_total = 0.0
    endif
    if (burned_total_noCWL<0 .AND. -1.0*burned_total_noCWL<carbon_cons_tol/2.0) then
       burned_total_noCWL = 0.0
    endif
    if (killed_live<0 .AND. -1.0*killed_total<carbon_cons_tol/2.0) then
       killed_live = 0.0
       killed_total = killed_wood
    endif
    
    ! Transfer killed C to soil
    select case (soil_carbon_option)
       case (SOILC_CENTURY, SOILC_CENTURY_BY_LAYER)
          soil%slow_soil_C(1) = soil%slow_soil_C(1) + (1.-fsc_liv) *killed_live &   ! "Alive"
                                                    + (1.-fsc_wood)*killed_wood   ! "Wood"
          soil%fast_soil_C(1) = soil%fast_soil_C(1) +     fsc_liv  *killed_live &   ! "Alive"
                                                    +     fsc_wood *killed_wood   ! "Dead"
!          soil%fsc_in(1) = soil%fsc_in(1) + killed_live
!          soil%ssc_in(1) = soil%ssc_in(1) + killed_wood
       case (SOILC_CORPSE)
          ! Dead biomass. (Does NOT include BLV as before. Instead, now like vegn_disturbance from CORPSE.)
          !!! EVERY TIME YOU CALL ADD_LITTER, IT CREATES A NEW COHORT. SO HERE, I'M CREATING TWO ABOVEGROUND LITTER
          !!! COHORTS AND TWO BELOWGROUND LITTER COHORTS. INSTEAD, COMBINE SO YOU CAN ADD ALL ABOVEGROUND (BELOWGROUND)
          !!! IN ONE CALL.
          call add_litter(soil%coarseWoodLitter, &
                          (/     fsc_wood *killed_wood, &   ! Cellulose
                             (1.-fsc_wood)*killed_wood, &   ! Lignin
                             0.0 &                                        ! Microbial products
                          /)*agf_bs &
                         )
          soil%coarsewoodlitter_fsc_in = soil%coarsewoodlitter_fsc_in +     fsc_wood *killed_wood
          soil%coarsewoodlitter_ssc_in = soil%coarsewoodlitter_ssc_in + (1.-fsc_wood)*killed_wood
          
          ! SSR20151217
          soil%ssr_cwl_unpr_fast_Ki_accum = soil%ssr_cwl_unpr_fast_Ki_accum +     fsc_wood *killed_wood
          soil%ssr_cwl_unpr_slow_Ki_accum = soil%ssr_cwl_unpr_slow_Ki_accum + (1.-fsc_wood)*killed_wood
          
          call add_root_litter(soil, vegn, &
                               (/     fsc_wood *killed_wood, &   ! Cellulose
                                  (1.-fsc_wood)*killed_wood, &   ! Lignin
                                  0.0 &                                        ! Microbial products
                               /)*(1-agf_bs) &
                              )   ! fsc_in and ssc_in are updated in add_root_litter
          ! Living biomass, aboveground: Leaves and virtual pool
          new_fast_C_ag =      fsc_liv *(killed_bl + killed_blv)
          new_slow_C_ag = (1.0-fsc_liv)*(killed_bl + killed_blv)
          call add_litter(soil%leafLitter,(/new_fast_C_ag,new_slow_C_ag,0.0/))
          soil%leaflitter_fsc_in = soil%leaflitter_fsc_in + new_fast_C_ag
          soil%leaflitter_ssc_in = soil%leaflitter_ssc_in + new_slow_C_ag
          
          ! SSR20151218
          soil%ssr_lfl_unpr_fast_Ki_accum = soil%ssr_lfl_unpr_fast_Ki_accum + new_fast_C_ag
          soil%ssr_lfl_unpr_slow_Ki_accum = soil%ssr_lfl_unpr_slow_Ki_accum + new_slow_C_ag
          
          ! Living biomass, belowground: Roots
          new_fast_C_bg =      fsc_froot *killed_br
          new_slow_C_bg = (1.0-fsc_froot)*killed_br
          call add_root_litter(soil,vegn,(/new_fast_C_bg,new_slow_C_bg,0.0/)) !fsc_in and ssc_in updated in add_root_litter
       case default
          call error_mesg('vegn_burn','soil_carbon_option is invalid. This should never happen. Contact developer.', FATAL)
    end select
    
    ! Accumulate smoke pool
    !!!!! The idea is to set csmoke_rate such that csmoke_pool is completely depleted via
    !!!!! constant emissions at each fast time step throughout the day.
    ! Deal with tiny negative values of csmoke_pool (before)
    if (vegn%csmoke_pool<0 .AND. -1.0*vegn%csmoke_pool<carbon_cons_tol/2.0) then
       vegn%csmoke_pool = 0.0
    endif
    vegn%csmoke_pool = vegn%csmoke_pool + burned_total   ! vegn%csmoke_pool SHOULD be zero
                                                         ! at the end of the day, but this
                                                         ! ensures no weirdness. (Maybe a
                                                         ! weird restart, e.g.)
    vegn%csmoke_rate_daily = vegn%csmoke_pool   ! kg C/(m2 day)
    vegn%csmoke_rate = vegn%csmoke_rate_daily * days_per_year ! kg C/(m2 yr)
    
    ! Biomass burned diagnostics
    vegn%burn_Cemit       = burned_total      *tile_area_m2
    vegn%burn_Cemit_noCWL = burned_total_noCWL*tile_area_m2
    vegn%burn_Ckill = killed_total*tile_area_m2

    ! More biomass burned diagnostics (2016-03-14)
    vegn%burn_Cemit_leaf = burned_bl
    vegn%burn_Cemit_stem = burned_bwood + burned_bsw
    vegn%burn_Cemit_litter = burned_leafLitter + burned_fineWoodLitter + burned_coarseWoodLitter
    vegn%burn_Cemit_litter_lf = burned_leafLitter
    vegn%burn_Cemit_litter_fw = burned_fineWoodLitter
    vegn%burn_Cemit_litter_cw = burned_coarseWoodLitter
    vegn%burn_Ckill_leaf = killed_bl
    vegn%burn_Ckill_stem = killed_wood
    vegn%burn_Ckill_root = killed_br
    
    call update_fire_agb(vegn,soil)
    fire_agb_1 = vegn%fire_agb
    
    if(is_watch_point() .OR. force_watch_point .OR. (is_watch_cell() .AND. force_watch_cell_burning)) then
        write(*,*) '######## After vegn_burn ########'
        write(*,*) 'burned_bl', burned_bl
        write(*,*) 'burned_bwood', burned_bwood
        write(*,*) 'burned_bsw', burned_bsw
        write(*,*) 'burned_blv', burned_blv
        write(*,*) 'burned_Cgain', burned_Cgain
        write(*,*) 'burned_Wgain', burned_Wgain
        if (soil_carbon_option==SOILC_CORPSE) then
           write(*,*) 'burned_leafLitter         ', burned_leafLitter
           write(*,*) 'burned_fineWoodLitter     ', burned_fineWoodLitter
           write(*,*) 'burned_coarseWoodLitter   ', burned_coarseWoodLitter
        endif
        write(*,*) 'burned_total'        , burned_total
        write(*,*) 'burned_total_noCWL', burned_total_noCWL
        write(*,*) 'vegn%csmoke_pool', vegn%csmoke_pool
        write(*,*) 'vegn%csmoke_rate', vegn%csmoke_rate
        write(*,*) ' '
        write(*,*) 'killed_bl', killed_bl
        write(*,*) 'killed_bwood', killed_bwood
        write(*,*) 'killed_bsw', killed_bsw
        write(*,*) 'killed_blv', killed_blv
        write(*,*) 'killed_br', killed_br
        write(*,*) 'killed_Cgain', killed_Cgain
        write(*,*) 'killed_Wgain', killed_Wgain
        write(*,*) 'killed_total', killed_total
        write(*,*) ' '
        write(*,*) 'bliving_0', bliving_0
        write(*,*) 'bliving_1', bliving_1
        write(*,*) 'fire_agb_0', fire_agb_0
        write(*,*) 'fire_agb_1', fire_agb_1
        write(*,*) ' '
        write(*,*) '(Each pair below should be equivalent.)'
        write(*,*) 'bl0 - bl1            ', bl0-bl1
        write(*,*) 'burned_bl + killed_bl', burned_bl + killed_bl
        write(*,*) 'bwood0 - bwood1            ', bwood0-bwood1
        write(*,*) 'burned_bwood + killed_bwood', burned_bwood + killed_bwood
        write(*,*) 'bsw0 - bsw1            ', bsw0-bsw1
        write(*,*) 'burned_bsw + killed_bsw', burned_bsw + killed_bsw
        write(*,*) 'blv0 - blv1            ', blv0-blv1
        write(*,*) 'burned_blv + killed_blv', burned_blv + killed_blv
        write(*,*) 'br0 - br1', br0-br1
        write(*,*) 'killed_br', killed_br
        write(*,*) 'Cgain0 - Cgain1            ', Cgain0-Cgain1
        write(*,*) 'burned_Cgain + killed_Cgain', burned_Cgain + killed_Cgain
        write(*,*) 'Wgain0 - Wgain1            ', Wgain0-Wgain1
        write(*,*) 'burned_Wgain + killed_Wgain', burned_Wgain + killed_Wgain
        write(*,*) ' '
        write(*,*) 'npp_prevDay_0', npp_prevDay_0
        write(*,*) 'npp_prevDay_1', npp_prevDay_1
        write(*,*) 'npp_prevDayTmp_0', npp_prevDayTmp_0
        write(*,*) 'npp_prevDayTmp_1', npp_prevDayTmp_1
        write(*,*) ' '
        write(*,*) 'burn_Cemit'      , vegn%burn_Cemit
        write(*,*) 'burn_Cemit_noCWL', vegn%burn_Cemit_noCWL
        write(*,*) 'burn_Ckill', vegn%burn_Ckill
        
        write(*,*) '##########################################################'
    endif
    
    ! Check values
    call check_var_range(vegn%burn_Cemit, 0.0, 10.**37, 'vegn_burn', 'burn_Cemit', lnd%time, FATAL)
    call check_var_range(vegn%burn_Ckill, 0.0, 10.**37, 'vegn_burn', 'burn_Ckill', lnd%time, FATAL)
    call check_var_range(vegn%csmoke_pool, 0.0, 10.**37, 'vegn_burn', 'vegn%csmoke_pool', lnd%time, FATAL)
    call check_var_range(vegn%csmoke_rate, 0.0, 10.**37, 'vegn_burn', 'vegn%csmoke_rate', lnd%time, FATAL)
    
end subroutine vegn_burn


subroutine vegn_fire_sendtiledata_Cburned(diag,tile_vegn)
    type(diag_buff_type), intent(inout) :: diag
    type(vegn_tile_type), intent(inout),optional :: tile_vegn
    real :: num_days
    
    if (present(tile_vegn)) then
       call send_tile_data(id_burn_Cemit,           tile_vegn%burn_Cemit,           diag)
       call send_tile_data(id_burn_Cemit_noCWL,     tile_vegn%burn_Cemit_noCWL,     diag)   ! SSR20151227
       call send_tile_data(id_burn_Ckill,           tile_vegn%burn_Ckill,           diag)
       call send_tile_data(id_burn_Cemit_leaf,      tile_vegn%burn_Cemit_leaf,      diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_stem,      tile_vegn%burn_Cemit_stem,      diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter,    tile_vegn%burn_Cemit_litter,    diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_lf, tile_vegn%burn_Cemit_litter_lf, diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_fw, tile_vegn%burn_Cemit_litter_fw, diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_cw, tile_vegn%burn_Cemit_litter_cw, diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_leaf,      tile_vegn%burn_Ckill_leaf,      diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_stem,      tile_vegn%burn_Ckill_stem,      diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_root,      tile_vegn%burn_Ckill_root,      diag)   ! SSR20160224
    else
       call send_tile_data(id_burn_Cemit,           0.0,diag)
       call send_tile_data(id_burn_Cemit_noCWL,     0.0,diag)   ! SSR20151227
       call send_tile_data(id_burn_Ckill,           0.0,diag)
       call send_tile_data(id_burn_Cemit_leaf,      0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_stem,      0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter,    0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_lf, 0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_fw, 0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Cemit_litter_cw, 0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_leaf,      0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_stem,      0.0,diag)   ! SSR20160224
       call send_tile_data(id_burn_Ckill_root,      0.0,diag)   ! SSR20160224
    endif
    
end subroutine vegn_fire_sendtiledata_Cburned


subroutine vegn_fire_BA_agri(vegn,Time,tile_area_km2,BA_mth,BF_mth)
  type(vegn_tile_type), intent(inout) :: vegn
  type(time_type), intent(in)  :: Time
  real, intent(in) :: tile_area_km2
  real, intent(out)   ::  BA_mth, BF_mth
  real :: num_days

  num_days = days_in_month(Time)
  
  ! Fcrop and Fpast are per-day rates. Since they're just being called once per month,
  ! they need to be multiplied by the number of days in the month.
  if (vegn%landuse==LU_CROP) then
     BF_mth = vegn%Fcrop * num_days
  elseif (vegn%landuse==LU_PAST) then
     BF_mth = vegn%Fpast * num_days
  endif
  
  ! BF can't be > 1!
  if (BF_mth > 1.0) then
     call check_var_range(BF_mth, 0.0, 1.0, 'vegn_fire_BA_agri', 'BF_mth', lnd%time, WARNING)
     if (vegn%landuse==LU_CROP) then
        write(*,*) 'Setting BF_mth for this CROP tile to 1.0.'
     elseif (vegn%landuse==LU_PAST) then
        write(*,*) 'Setting BF_mth for this PAST tile to 1.0.'
     endif
     BF_mth = 1.0
  endif
  
  BA_mth = BF_mth * tile_area_km2
  
  if (is_watch_point()) then
     write(*,*) '##### checkpoint vegn_fire_BA_agri #####'
     if (vegn%landuse==LU_CROP) then
        write(*,*) 'CROP tile. Fcrop = ', vegn%Fcrop
     elseif (vegn%landuse==LU_PAST) then
        write(*,*) 'PAST tile. Fpast = ', vegn%Fpast
     endif
     write(*,*) 'BF_mth', BF_mth
     write(*,*) '########################################'
  endif
  
end subroutine vegn_fire_BA_agri


subroutine update_Nfire_BA_fast(diag,i,j,tile_area, &
                                fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                                BAperFire_0, &
                                vegn_cohort_1_species, &
                                vegn_trop_code, &   ! SSR20150831
                                lightning, popD, GDPpc, &
                                max_fire_size, &   ! SSR20150727
                                vegn_burned_frac, &
                                latitude, &
                                ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta)   ! SSR20151009
  type(diag_buff_type), intent(inout) :: diag
  integer, intent(in) :: i,j   ! Coordinates of current point, for fire data
  real, intent(in)    :: tile_area   ! Area of tile (m2)
  real, intent(in)    :: fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb
  real, intent(in)    :: BAperFire_0
  integer, intent(in) :: vegn_cohort_1_species
  integer, intent(in) :: vegn_trop_code   ! SSR20150831
  real, intent(in)    :: lightning ! Lightning flash density (flashes/km2/day)
  real, intent(in)    :: popD      ! Population density (people/km2)
  real, intent(in)    :: GDPpc     ! GDP per capita ($/person)
  real, intent(in)    :: max_fire_size   ! SSR20150727
  real, intent(inout) :: vegn_burned_frac     ! Fraction of tile burned so far today
  
  real                :: tile_area_km2   ! Area of tile (km2)
  real, intent(in)    :: latitude   ! Latitude (radians)
  real, intent(in)    :: ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta   ! SSR20151009
  real                :: fire_fn_popD_NF, fire_fn_popD_BA, &
                         fire_fn_GDPpc_NF, fire_fn_GDPpc_BA
  real                :: In, Ia
  real                :: Nfire_perKm2
  real                :: BA, Nfire, BA_rate, BF, BF_rate, Nfire_rate
  real                :: BA_DERIVwrt_alphaM, &
                         BA_DERIVwrt_thetaE, &
                         BA_DERIVwrt_THETAparam1, BA_DERIVwrt_THETAparam2, &   ! SSR20160203
                         BA_DERIVwrt_AGBparam1, BA_DERIVwrt_AGBparam2, &
                         BA_DERIVwrt_RHparam1, BA_DERIVwrt_RHparam2, &
                         BA_DERIVwrt_IaParam1, BA_DERIVwrt_IaParam2, &
                         BA_DERIVwrt_popDsupp_NF_eps1, BA_DERIVwrt_popDsupp_NF_eps2, BA_DERIVwrt_popDsupp_NF_eps3, &
                         BA_DERIVwrt_fireDur_c4, BA_DERIVwrt_fireDur_c3, &
                         BA_DERIVwrt_fireDur_dt, BA_DERIVwrt_fireDur_tt, BA_DERIVwrt_fireDur_et, &
                         BA_DERIVwrt_ROSmax_c4, BA_DERIVwrt_ROSmax_c3, &
                         BA_DERIVwrt_ROSmax_dt, BA_DERIVwrt_ROSmax_tt, BA_DERIVwrt_ROSmax_et, &
!                         BA_DERIVwrt_ROSmax_ts   ! SSR20150831
                         BA_DERIVwrt_ROSmax_tshr, BA_DERIVwrt_ROSmax_tsav, &   ! SSR20160211
                         BA_DERIVwrt_magicScalar   ! SSR20160222
  real                :: num_days   ! Number of days in month
  real                :: BA_reduction, vegn_unburned_area_before
  real                :: BAperFire_1, BAperFire_2, BAperFire_3
  logical             :: tooMuchFire


  tile_area_km2 = tile_area / 1000000.   ! Convert m2 to km2
  
  call vegn_fire_fn_popD(vegn_cohort_1_species,popD,fire_fn_popD_NF,fire_fn_popD_BA)
  call vegn_fire_fn_GDPpc(vegn_cohort_1_species,GDPpc,popD,fire_fn_GDPpc_NF,fire_fn_GDPpc_BA)
  call vegn_fire_In(latitude,lightning,In)
  call vegn_fire_Ia(popD,Ia)
  call vegn_fire_Nfire(In, Ia, &
                       fire_fn_theta, fire_fn_rh, &
                       fire_fn_Tca, fire_fn_agb, &
                       fire_fn_popD_NF, fire_fn_GDPpc_NF, &
                       Nfire_perKm2)
  call vegn_fire_BA_ntrlscnd(Nfire_perKm2, tile_area_km2, &
                             max_fire_size, &
                             vegn_burned_frac, &
                             BAperFire_0, BAperFire_1, BAperFire_2, BAperFire_3, &
                             fire_fn_popD_BA, fire_fn_GDPpc_BA, &
                             BA,BF,Nfire,BA_rate,BF_rate,Nfire_rate, &
                             BA_reduction, vegn_unburned_area_before)
                             
  ! Did we have to reduce burned area it exceeded the unburned area in the grid cell,
  ! or because fire size was larger than that allowed by fragmentation? (SSR20150921)
  if (BA_reduction > 0.0) then
     tooMuchFire = .TRUE.
  else
     tooMuchFire = .FALSE.
  endif

  if (do_calc_derivs) then
     if (vegn_unburned_area_before > 0.0 &
     .AND..NOT.(zero_tmf .AND. tooMuchFire) ) then   ! SSR20150921
        call calc_fire_derivs(& ! input
                                BA_rate, BAperFire_0, BA_reduction, &
                                In, Ia, &
                                fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                                fire_fn_popD_NF, fire_fn_GDPpc_NF, &
                                fire_fn_popD_BA, fire_fn_GDPpc_BA, &
                                vegn_cohort_1_species, vegn_unburned_area_before, &
                                vegn_trop_code, &   ! SSR20150831
                                ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta, &   ! SSR20151009
                                ! output
                                BA_DERIVwrt_alphaM, BA_DERIVwrt_IaParam1, BA_DERIVwrt_IaParam2, &
                                BA_DERIVwrt_AGBparam1, BA_DERIVwrt_AGBparam2, &
                                BA_DERIVwrt_RHparam1, BA_DERIVwrt_RHparam2, &
                                BA_DERIVwrt_thetaE, &
                                BA_DERIVwrt_THETAparam1, BA_DERIVwrt_THETAparam2, &   ! SSR20160203
                                BA_DERIVwrt_popDsupp_NF_eps1, BA_DERIVwrt_popDsupp_NF_eps2, BA_DERIVwrt_popDsupp_NF_eps3, &
                                BA_DERIVwrt_fireDur_c4, BA_DERIVwrt_fireDur_c3, &
                                BA_DERIVwrt_fireDur_dt, BA_DERIVwrt_fireDur_tt, BA_DERIVwrt_fireDur_et, &
                                BA_DERIVwrt_ROSmax_c4, BA_DERIVwrt_ROSmax_c3, &
                                BA_DERIVwrt_ROSmax_dt, BA_DERIVwrt_ROSmax_tt, BA_DERIVwrt_ROSmax_et, &
!                                BA_DERIVwrt_ROSmax_ts)   ! SSR20150831
                                BA_DERIVwrt_ROSmax_tshr, BA_DERIVwrt_ROSmax_tsav)   ! SSR20160211
        BA_DERIVwrt_magicScalar = BA_rate / magic_scalar   ! SSR20160222
     else
        BA_DERIVwrt_alphaM           = 0.0
        BA_DERIVwrt_IaParam1         = 0.0
        BA_DERIVwrt_IaParam2         = 0.0
        BA_DERIVwrt_AGBparam1        = 0.0
        BA_DERIVwrt_AGBparam2        = 0.0
        BA_DERIVwrt_RHparam1         = 0.0
        BA_DERIVwrt_RHparam2         = 0.0
        BA_DERIVwrt_thetaE           = 0.0
        BA_DERIVwrt_popDsupp_NF_eps1 = 0.0
        BA_DERIVwrt_popDsupp_NF_eps2 = 0.0
        BA_DERIVwrt_popDsupp_NF_eps3 = 0.0
        BA_DERIVwrt_fireDur_c4       = 0.0
        BA_DERIVwrt_fireDur_c3       = 0.0
        BA_DERIVwrt_fireDur_dt       = 0.0
        BA_DERIVwrt_fireDur_tt       = 0.0
        BA_DERIVwrt_fireDur_et       = 0.0
        BA_DERIVwrt_ROSmax_c4        = 0.0
        BA_DERIVwrt_ROSmax_c3        = 0.0
        BA_DERIVwrt_ROSmax_dt        = 0.0
        BA_DERIVwrt_ROSmax_tt        = 0.0
        BA_DERIVwrt_ROSmax_et        = 0.0
!        BA_DERIVwrt_ROSmax_ts        = 0.0   ! SSR20150831
        BA_DERIVwrt_ROSmax_tshr      = 0.0   ! SSR20160211
        BA_DERIVwrt_ROSmax_tsav      = 0.0   ! SSR20160211
        BA_DERIVwrt_magicScalar      = 0.0   ! SSR20160222
     endif
  endif
  
  ! Save fire stuff for diagnostics
  if (.NOT. minimal_fire_diagnostics) then
     call send_tile_data(id_BAperFire_1,        BAperFire_1,        diag)
     call send_tile_data(id_BAperFire_2,        BAperFire_2,        diag)
     call send_tile_data(id_BAperFire_3,        BAperFire_3,        diag)
     call send_tile_data(id_BAperFire_max,    max_fire_size*1e-6, diag)
     call send_tile_data(id_fire_fn_popD_NF,  fire_fn_popD_NF,  diag)
     call send_tile_data(id_fire_fn_popD_BA,  fire_fn_popD_BA,  diag)
     call send_tile_data(id_fire_fn_GDPpc_NF, fire_fn_GDPpc_NF, diag)
     call send_tile_data(id_fire_fn_GDPpc_BA, fire_fn_GDPpc_BA, diag)
     call send_tile_data(id_In,               In,               diag)
     call send_tile_data(id_Ia,               Ia,               diag)
     call send_tile_data(id_Nfire_perKm2,     Nfire_perKm2,     diag)
  endif
  call send_tile_data(id_Nfire_rate,       Nfire_rate,       diag)
  call send_tile_data(id_BF_rate, BF_rate, diag)
  call send_tile_data(id_BA_rate, BA_rate, diag)

  if (do_calc_derivs) then
     call send_tile_data(id_BA_DERIVwrt_alphaM, BA_DERIVwrt_alphaM, diag)
     call send_tile_data(id_BA_DERIVwrt_thetaE, BA_DERIVwrt_thetaE, diag)
     
     ! SSR20160203
     call send_tile_data(id_BA_DERIVwrt_THETAparam1, BA_DERIVwrt_THETAparam1, diag)
     call send_tile_data(id_BA_DERIVwrt_THETAparam2, BA_DERIVwrt_THETAparam2, diag)
     
     call send_tile_data(id_BA_DERIVwrt_AGBparam1, BA_DERIVwrt_AGBparam1, diag)
     call send_tile_data(id_BA_DERIVwrt_AGBparam2, BA_DERIVwrt_AGBparam2, diag)
     call send_tile_data(id_BA_DERIVwrt_RHparam1, BA_DERIVwrt_RHparam1, diag)
     call send_tile_data(id_BA_DERIVwrt_RHparam2, BA_DERIVwrt_RHparam2, diag)
     call send_tile_data(id_BA_DERIVwrt_IaParam1, BA_DERIVwrt_IaParam1, diag)
     call send_tile_data(id_BA_DERIVwrt_IaParam2, BA_DERIVwrt_IaParam2, diag)
     call send_tile_data(id_BA_DERIVwrt_popDsupp_NF_eps1, BA_DERIVwrt_popDsupp_NF_eps1, diag)
     call send_tile_data(id_BA_DERIVwrt_popDsupp_NF_eps2, BA_DERIVwrt_popDsupp_NF_eps2, diag)
     call send_tile_data(id_BA_DERIVwrt_popDsupp_NF_eps3, BA_DERIVwrt_popDsupp_NF_eps3, diag)
     call send_tile_data(id_BA_DERIVwrt_popDsupp_NF_eps3, BA_DERIVwrt_popDsupp_NF_eps3, diag)
     call send_tile_data(id_BA_DERIVwrt_fireDur_c4, BA_DERIVwrt_fireDur_c4, diag)
     call send_tile_data(id_BA_DERIVwrt_fireDur_c3, BA_DERIVwrt_fireDur_c3, diag)
     call send_tile_data(id_BA_DERIVwrt_fireDur_dt, BA_DERIVwrt_fireDur_dt, diag)
     call send_tile_data(id_BA_DERIVwrt_fireDur_tt, BA_DERIVwrt_fireDur_tt, diag)
     call send_tile_data(id_BA_DERIVwrt_fireDur_et, BA_DERIVwrt_fireDur_et, diag)
     call send_tile_data(id_BA_DERIVwrt_ROSmax_c4, BA_DERIVwrt_ROSmax_c4, diag)
     call send_tile_data(id_BA_DERIVwrt_ROSmax_c3, BA_DERIVwrt_ROSmax_c3, diag)
     call send_tile_data(id_BA_DERIVwrt_ROSmax_dt, BA_DERIVwrt_ROSmax_dt, diag)
     call send_tile_data(id_BA_DERIVwrt_ROSmax_tt, BA_DERIVwrt_ROSmax_tt, diag)
     call send_tile_data(id_BA_DERIVwrt_ROSmax_et, BA_DERIVwrt_ROSmax_et, diag)
!     call send_tile_data(id_BA_DERIVwrt_ROSmax_ts, BA_DERIVwrt_ROSmax_et, diag)   ! SSR20150831
     call send_tile_data(id_BA_DERIVwrt_ROSmax_tshr, BA_DERIVwrt_ROSmax_tshr, diag)   ! SSR20160211
     call send_tile_data(id_BA_DERIVwrt_ROSmax_tsav, BA_DERIVwrt_ROSmax_tsav, diag)   ! SSR20160211
     call send_tile_data(id_BA_DERIVwrt_magicScalar, BA_DERIVwrt_magicScalar, diag)   ! SSR20160222
     call send_tile_data(id_Ia_DERIVwrt_alphaM, Ia_DERIVwrt_alphaM, diag)
     call send_tile_data(id_fire_fn_theta_DERIVwrt_thetaE, fire_fn_theta_DERIVwrt_thetaE, diag)
     call send_tile_data(id_Ia_DERIVwrt_IaParam1, Ia_DERIVwrt_IaParam1, diag)
     call send_tile_data(id_Ia_DERIVwrt_IaParam2, Ia_DERIVwrt_IaParam2, diag)
     call send_tile_data(id_fire_fn_agb_DERIVwrt_param1, fire_fn_agb_DERIVwrt_param1, diag)
     call send_tile_data(id_fire_fn_agb_DERIVwrt_param2, fire_fn_agb_DERIVwrt_param2, diag)
     call send_tile_data(id_fire_fn_rh_DERIVwrt_param1, fire_fn_rh_DERIVwrt_param1, diag)
     call send_tile_data(id_fire_fn_rh_DERIVwrt_param2, fire_fn_rh_DERIVwrt_param2, diag)
     call send_tile_data(id_popDsupp_NF_DERIVwrt_eps1, popDsupp_NF_DERIVwrt_eps1, diag)
     call send_tile_data(id_popDsupp_NF_DERIVwrt_eps2, popDsupp_NF_DERIVwrt_eps2, diag)
     call send_tile_data(id_popDsupp_NF_DERIVwrt_eps3, popDsupp_NF_DERIVwrt_eps3, diag)
     call send_tile_data(id_BAperFire0_DERIVwrt_fireDur, BAperFire0_DERIVwrt_fireDur, diag)
     call send_tile_data(id_BAperFire0_DERIVwrt_ROSmax, BAperFire0_DERIVwrt_ROSmax, diag)
  endif
end subroutine update_Nfire_BA_fast


subroutine send_tile_data_BABF_forAgri(diag,BF_mth,BA_mth)
   type(diag_buff_type), intent(inout) :: diag
   real, intent(in)    :: BF_mth, BA_mth
   
   call send_tile_data(id_BF_rate, BF_mth, diag)   ! Remember, sending BF_mth when it's calculated and 0 otherwise results in a good per-day rate for the month as a whole
   call send_tile_data(id_BA_rate, BA_mth, diag)   ! Remember, sending BA_mth when it's calculated and 0 otherwise results in a good per-day rate for the month as a whole

end subroutine send_tile_data_BABF_forAgri


subroutine calc_fire_derivs(&
                            ! input
                            BA_rate, BAperFire_0, BA_reduction, &
                            In, Ia, &
                            fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                            fire_fn_popD_NF, fire_fn_GDPpc_NF, &
                            fire_fn_popD_BA, fire_fn_GDPpc_BA, &
                            vegn_cohort_1_species, vegn_unburned_area, &
                            vegn_trop_code, &   ! SSR20150831
                            ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta, &   ! SSR20151009
                            ! output
                            BA_DERIVwrt_alphaM, BA_DERIVwrt_IaParam1, BA_DERIVwrt_IaParam2, &
                            BA_DERIVwrt_AGBparam1, BA_DERIVwrt_AGBparam2, &
                            BA_DERIVwrt_RHparam1, BA_DERIVwrt_RHparam2, &
                            BA_DERIVwrt_thetaE, &
                            BA_DERIVwrt_THETAparam1, BA_DERIVwrt_THETAparam2, &   ! SSR20160203
                            BA_DERIVwrt_popDsupp_NF_eps1, BA_DERIVwrt_popDsupp_NF_eps2, BA_DERIVwrt_popDsupp_NF_eps3, &
                            BA_DERIVwrt_fireDur_c4, BA_DERIVwrt_fireDur_c3, &
                            BA_DERIVwrt_fireDur_dt, BA_DERIVwrt_fireDur_tt, BA_DERIVwrt_fireDur_et, &
                            BA_DERIVwrt_ROSmax_c4, BA_DERIVwrt_ROSmax_c3, &
                            BA_DERIVwrt_ROSmax_dt, BA_DERIVwrt_ROSmax_tt, BA_DERIVwrt_ROSmax_et, &
!                            BA_DERIVwrt_ROSmax_ts)   ! SSR20150831
                            BA_DERIVwrt_ROSmax_tshr, BA_DERIVwrt_ROSmax_tsav)   ! SSR20160211
   real,    intent(in)   ::   BA_rate, BAperFire_0, BA_reduction, In, Ia, &
                              fire_fn_popD_NF, fire_fn_GDPpc_NF, &
                              fire_fn_popD_BA, fire_fn_GDPpc_BA, &
                              fire_fn_theta, fire_fn_rh, fire_fn_Tca, fire_fn_agb, &
                              vegn_unburned_area, &
                              ROSmax, gW, fire_dur, HB, LB, rh, theta, C_beta   ! SSR20151009
   integer, intent(in)   ::   vegn_cohort_1_species, &
                              vegn_trop_code   ! SSR20150831
   real,    intent(out)  ::   BA_DERIVwrt_alphaM, BA_DERIVwrt_IaParam1, BA_DERIVwrt_IaParam2, &
                              BA_DERIVwrt_AGBparam1, BA_DERIVwrt_AGBparam2, &
                              BA_DERIVwrt_RHparam1, BA_DERIVwrt_RHparam2, &
                              BA_DERIVwrt_thetaE, &
                              BA_DERIVwrt_THETAparam1, BA_DERIVwrt_THETAparam2, &   ! SSR20160203
                              BA_DERIVwrt_popDsupp_NF_eps1, BA_DERIVwrt_popDsupp_NF_eps2, BA_DERIVwrt_popDsupp_NF_eps3, &
                              BA_DERIVwrt_fireDur_c4, BA_DERIVwrt_fireDur_c3, &
                              BA_DERIVwrt_fireDur_dt, BA_DERIVwrt_fireDur_tt, BA_DERIVwrt_fireDur_et, &
                              BA_DERIVwrt_ROSmax_c4, BA_DERIVwrt_ROSmax_c3, &
                              BA_DERIVwrt_ROSmax_dt, BA_DERIVwrt_ROSmax_tt, BA_DERIVwrt_ROSmax_et, &
                              !BA_DERIVwrt_ROSmax_ts   ! SSR20150831
                              BA_DERIVwrt_ROSmax_tshr, BA_DERIVwrt_ROSmax_tsav   ! SSR20160211
   real  ::  BA_DERIVwrt_fireDur_TMP, BA_DERIVwrt_ROSmax_TMP, fast_to_daily
   real  ::  RHderiv_X, RHderiv_Z   ! SSR20151009
   real  ::  THETAderiv_X, THETAderiv_Z   ! SSR20151216

!!!! NOTE:
!    BA_rate = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF
!       * vegn_unburned_area * 86400. / dt_fast &
!       * (BAperFire_0 * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   
   fast_to_daily = 86400. / dt_fast
   
   ! Number of ignitions
   BA_DERIVwrt_alphaM = Ia_DERIVwrt_alphaM * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   BA_DERIVwrt_IaParam1 = Ia_DERIVwrt_IaParam1 * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   BA_DERIVwrt_IaParam2 = Ia_DERIVwrt_IaParam2 * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)

   ! Fraction of ignitions becoming fires
   BA_DERIVwrt_AGBparam1 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca & &
                       * fire_fn_agb_DERIVwrt_param1 * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   BA_DERIVwrt_AGBparam2 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb_DERIVwrt_param2 * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
                         
                         
                         
!   BA_DERIVwrt_RHparam1 = (In + Ia) * fire_fn_theta * fire_fn_rh_DERIVwrt_param1 * fire_fn_Tca &
!                       * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                       * vegn_unburned_area * fast_to_daily &
!                       * (BAperFire_0 * (1.0 - BA_reduction) &
!                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
!   BA_DERIVwrt_RHparam2 = (In + Ia) * fire_fn_theta * fire_fn_rh_DERIVwrt_param2 * fire_fn_Tca &
!                       * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                       * vegn_unburned_area * fast_to_daily &
!                       * (BAperFire_0 * (1.0 - BA_reduction) &
!                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)

!!    New method of calculating derivative w/r/t RH parameters. Ignores
!!    fire_fn_rh_DERIVwrt_param[1,2], instead doing everything here.
!    RHderiv_X = (In + Ia) * fire_fn_theta * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                * vegn_unburned_area * fast_to_daily
!    if ((.NOT. C_beta_params_likeRH) .OR. theta_ROSeffect_asFnTheta) then
!       RHderiv_Z = ((ROSmax * gW * C_beta * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6)) &
!                 * (BAperFire_0 * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA)
!       BA_DERIVwrt_RHparam1 = -(3.*RHderiv_X*RHderiv_Z*((rh - RH_up)**3.))/((RH_lo - RH_up)**4.)
!       BA_DERIVwrt_RHparam2 = (3.*RHderiv_X*RHderiv_Z*((rh - RH_up)**3.))/((RH_lo - RH_up)**4.) - (3.*RHderiv_X*RHderiv_Z*((rh - RH_up)**2.))/((RH_lo - RH_up)**3.)
!    else
!       RHderiv_Z = (ROSmax * gW * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6) &
!                 * (BAperFire_0 * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA)
!       BA_DERIVwrt_RHparam1 = -(5.*RHderiv_X*RHderiv_Z*((rh - RH_up)**3.)*((theta - RH_up)**2.))/((RH_lo - RH_up)**6.)
!       BA_DERIVwrt_RHparam2 = (5.*RHderiv_X*RHderiv_Z*((rh - RH_up)**3.)*((theta - RH_up)**2.))/((RH_lo - RH_up)**6.) - (3.*RHderiv_X*RHderiv_Z*((rh - RH_up)**2.)*((theta - RH_up)**2.))/((RH_lo - RH_up)**5.) - (RHderiv_X*RHderiv_Z*(2.*theta - 2.*RH_up)*((rh - RH_up)**3.))/((RH_lo - RH_up)**5.)
!    endif

   ! SSR20160202
   RHderiv_X = (In + Ia) * fire_fn_theta * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
               * vegn_unburned_area * fast_to_daily &
               * ((ROSmax * gW * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6)) &
               * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA
   if ((.NOT. C_beta_params_likeRH) .OR. theta_ROSeffect_asFnTheta) then
      RHderiv_X = RHderiv_X * (C_beta**2.)
   else
      call error_mesg('calc_fire_derivs',&
                      'Add code to get RH derivatives for when C_beta_params_likeRH==TRUE and theta_ROSeffect_asFnTheta==FALSE.', &
                      FATAL)
   endif
   BA_DERIVwrt_RHparam1 = RHderiv_X * fire_TOTALfn_rh_DERIVwrt_param1
   BA_DERIVwrt_RHparam2 = RHderiv_X * fire_TOTALfn_rh_DERIVwrt_param2
   
!    if (theta_ROSeffect_asFnTheta) then
!       ! New method of calculating derivative w/r/t theta_extinction. Ignores
!       ! fire_fn_theta_DERIVwrt_thetaE, instead doing everything here.
!       THETAderiv_X = (In + Ia) * fire_fn_rh * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                    * vegn_unburned_area * fast_to_daily
!       THETAderiv_Z = (ROSmax * gW * fire_fn_rh * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6) &
!                    * (BAperFire_0 * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA)
!       if (thetaE_already_squared) then
!          BA_DERIVwrt_thetaE = (3*pi*THETAderiv_X*THETAderiv_Z*(theta**2)*exp(-(3*pi*(theta**2))/theta_extinction))/(theta_extinction**2)
!       else
!          BA_DERIVwrt_thetaE = (6*pi*THETAderiv_X*THETAderiv_Z*(theta**2)*exp(-(3*pi*(theta**2))/(theta_extinction**2)))/(theta_extinction**3)
!       endif
!    else
!       BA_DERIVwrt_thetaE = (In + Ia) * fire_fn_theta_DERIVwrt_thetaE * fire_fn_rh * fire_fn_Tca &
!                           * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                           * vegn_unburned_area * fast_to_daily &
!                           * (BAperFire_0 * (1.0 - BA_reduction) &
!                             * fire_fn_popD_BA * fire_fn_GDPpc_BA)
!    endif
   
!   ! SSR20160202
!   if (theta_ROSeffect_asFnTheta) then
!      THETAderiv_X = (In + Ia) * fire_fn_rh * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
!                     * vegn_unburned_area * fast_to_daily &
!                     * ((ROSmax * gW * fire_fn_rh * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6)) &
!                     * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA
!   else
!      call error_mesg('calc_fire_derivs',&
!                      'Add code to get RH derivatives for when theta_ROSeffect_asFnTheta==FALSE.', &
!                      FATAL)
!   endif
   
   ! SSR20160203
   if (theta_ROSeffect_asFnTheta) then
      THETAderiv_X = (In + Ia) * fire_fn_rh * fire_fn_Tca * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                     * vegn_unburned_area * fast_to_daily &
                     * ((ROSmax * gW * fire_fn_rh * fire_dur * (1 + 1/HB))**2. * pi / (4 * LB * 10**6)) &
                     * (1.0 - BA_reduction) * fire_fn_popD_BA * fire_fn_GDPpc_BA
      if (fire_option_fTheta==FIRE_THETA_LI2012) then
         BA_DERIVwrt_thetaE = THETAderiv_X * fire_TOTALfn_theta_DERIVwrt_thetaE
      else
         BA_DERIVwrt_THETAparam1 = THETAderiv_X * fire_TOTALfn_theta_DERIVwrt_param1
         BA_DERIVwrt_THETAparam2 = THETAderiv_X * fire_TOTALfn_theta_DERIVwrt_param2
      endif
   else
      call error_mesg('calc_fire_derivs',&
                      'Add code to get theta derivatives for when theta_ROSeffect_asFnTheta==FALSE.', &
                      FATAL)
   endif
   
   BA_DERIVwrt_popDsupp_NF_eps1 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (0.-popDsupp_NF_DERIVwrt_eps1) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   BA_DERIVwrt_popDsupp_NF_eps2 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (0.-popDsupp_NF_DERIVwrt_eps2) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
   BA_DERIVwrt_popDsupp_NF_eps3 = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                       * fire_fn_agb * (0.-popDsupp_NF_DERIVwrt_eps3) * fire_fn_GDPpc_NF &
                       * vegn_unburned_area * fast_to_daily &
                       * (BAperFire_0 * (1.0 - BA_reduction) &
                         * fire_fn_popD_BA * fire_fn_GDPpc_BA)
                         
   ! Fire duration
   BA_DERIVwrt_fireDur_c4 = 0.0
   BA_DERIVwrt_fireDur_c3 = 0.0
   BA_DERIVwrt_fireDur_dt = 0.0
   BA_DERIVwrt_fireDur_tt = 0.0
   BA_DERIVwrt_fireDur_et = 0.0
   if (BA_reduction==0.0) then   ! Otherwise, the slope of BA w/r/t fire duration should
                                 ! be zero because BAperFire0 is too big anyway.
      BA_DERIVwrt_fireDur_TMP = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                          * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                          * vegn_unburned_area * fast_to_daily &
                          * (BAperFire0_DERIVwrt_fireDur * (1.0 - BA_reduction) &
                            * fire_fn_popD_BA * fire_fn_GDPpc_BA)
      if (vegn_cohort_1_species==SP_C4GRASS)  BA_DERIVwrt_fireDur_c4 = BA_DERIVwrt_fireDur_TMP
      if (vegn_cohort_1_species==SP_C3GRASS)  BA_DERIVwrt_fireDur_c3 = BA_DERIVwrt_fireDur_TMP
      if (vegn_cohort_1_species==SP_TEMPDEC)  BA_DERIVwrt_fireDur_dt = BA_DERIVwrt_fireDur_TMP
      if (vegn_cohort_1_species==SP_TROPICAL) BA_DERIVwrt_fireDur_tt = BA_DERIVwrt_fireDur_TMP
      if (vegn_cohort_1_species==SP_EVERGR)   BA_DERIVwrt_fireDur_et = BA_DERIVwrt_fireDur_TMP
   endif
   
   ! Rate of spread
   BA_DERIVwrt_ROSmax_c4 = 0.0
   BA_DERIVwrt_ROSmax_c3 = 0.0
   BA_DERIVwrt_ROSmax_dt = 0.0
   BA_DERIVwrt_ROSmax_tt = 0.0
   BA_DERIVwrt_ROSmax_et = 0.0
!   BA_DERIVwrt_ROSmax_ts = 0.0   ! SSR20150831
   BA_DERIVwrt_ROSmax_tshr = 0.0   ! SSR20160211
   BA_DERIVwrt_ROSmax_tsav = 0.0   ! SSR20160211
   if (BA_reduction==0.0) then   ! Otherwise, the slope of BA w/r/t fire duration should
                                 ! be zero because BAperFire0 is too big anyway.
      BA_DERIVwrt_ROSmax_TMP = (In + Ia) * fire_fn_theta * fire_fn_rh * fire_fn_Tca &
                          * fire_fn_agb * (1.-fire_fn_popD_NF) * fire_fn_GDPpc_NF &
                          * vegn_unburned_area * fast_to_daily &
                          * (BAperFire0_DERIVwrt_ROSmax * (1.0 - BA_reduction) &
                            * fire_fn_popD_BA * fire_fn_GDPpc_BA)
      if (vegn_cohort_1_species==SP_C4GRASS)  BA_DERIVwrt_ROSmax_c4 = BA_DERIVwrt_ROSmax_TMP
      if (vegn_cohort_1_species==SP_C3GRASS) then
         if (lock_ROSmaxC3_to_ROSmaxC4) then
            BA_DERIVwrt_ROSmax_c4 = BA_DERIVwrt_ROSmax_TMP
         else
            BA_DERIVwrt_ROSmax_c3 = BA_DERIVwrt_ROSmax_TMP
         endif
      endif
      if (vegn_cohort_1_species==SP_TEMPDEC)  BA_DERIVwrt_ROSmax_dt = BA_DERIVwrt_ROSmax_TMP
      if (vegn_cohort_1_species==SP_TROPICAL) then
         ! SSR20150831
         if (vegn_trop_code==TROP_FOR &
             .OR. (vegn_trop_code==TROP_SHR .AND. ROS_max_TROPSHR <= 0.0) &
             .OR. (vegn_trop_code==TROP_SAV .AND. ROS_max_TROPSAV <= 0.0)) then   ! Because in this case, we fell back on ROS_max_TROPFOR
            BA_DERIVwrt_ROSmax_tt = BA_DERIVwrt_ROSmax_TMP
!         else
!            BA_DERIVwrt_ROSmax_ts = BA_DERIVwrt_ROSmax_TMP
!         endif
         ! SSR20160211
         elseif (vegn_trop_code==TROP_SHR) then
            BA_DERIVwrt_ROSmax_tshr = BA_DERIVwrt_ROSmax_TMP
         elseif (vegn_trop_code==TROP_SAV) then
            BA_DERIVwrt_ROSmax_tsav = BA_DERIVwrt_ROSmax_TMP
         endif
      endif
      if (vegn_cohort_1_species==SP_EVERGR)   BA_DERIVwrt_ROSmax_et = BA_DERIVwrt_ROSmax_TMP
   endif
   
   if (is_watch_point()) then
      write(*,*) '### calc_fire_derivs_v2 ###'
      write(*,*) 'BA_rate', BA_rate
      write(*,*) 'BAperFire_0', BAperFire_0
      write(*,*) 'BA_reduction', BA_reduction
      write(*,*) 'dt_fast', dt_fast
      write(*,*) 'fast_to_daily', fast_to_daily
      __DEBUG2__(In, Ia)
      __DEBUG2__(fire_fn_popD_NF,fire_fn_GDPpc_NF)
      __DEBUG2__(fire_fn_popD_BA,fire_fn_GDPpc_BA)
      __DEBUG2__(fire_fn_theta,fire_fn_rh)
      __DEBUG2__(fire_fn_Tca,fire_fn_agb)
      write(*,*) 'vegn_unburned_area', vegn_unburned_area
      write(*,*) 'vegn_cohort_1_species', vegn_cohort_1_species
      __DEBUG2__(Ia_DERIVwrt_alphaM,BA_DERIVwrt_alphaM)
      __DEBUG2__(Ia_DERIVwrt_IaParam1,BA_DERIVwrt_IaParam1)
      __DEBUG2__(Ia_DERIVwrt_IaParam2,BA_DERIVwrt_IaParam2)
      __DEBUG2__(fire_fn_agb_DERIVwrt_param1,BA_DERIVwrt_AGBparam1)
      __DEBUG2__(fire_fn_agb_DERIVwrt_param2,BA_DERIVwrt_AGBparam2)
!      __DEBUG2__(fire_fn_rh_DERIVwrt_param1,BA_DERIVwrt_RHparam1)
!      __DEBUG2__(fire_fn_rh_DERIVwrt_param2,BA_DERIVwrt_RHparam2)
      __DEBUG2__(fire_TOTALfn_rh_DERIVwrt_param1,BA_DERIVwrt_RHparam1)
      __DEBUG2__(fire_TOTALfn_rh_DERIVwrt_param2,BA_DERIVwrt_RHparam2)
      __DEBUG2__(fire_fn_theta_DERIVwrt_thetaE,BA_DERIVwrt_thetaE)
      __DEBUG2__(popDsupp_NF_DERIVwrt_eps1,BA_DERIVwrt_popDsupp_NF_eps1)
      __DEBUG2__(popDsupp_NF_DERIVwrt_eps2,BA_DERIVwrt_popDsupp_NF_eps2)
      __DEBUG2__(popDsupp_NF_DERIVwrt_eps3,BA_DERIVwrt_popDsupp_NF_eps3)
      write(*,*) 'BAperFire0_DERIVwrt_fireDur', BAperFire0_DERIVwrt_fireDur
      __DEBUG2__(BA_DERIVwrt_fireDur_c4,BA_DERIVwrt_fireDur_c3)
      __DEBUG3__(BA_DERIVwrt_fireDur_dt,BA_DERIVwrt_fireDur_tt,BA_DERIVwrt_fireDur_et)
      write(*,*) 'BAperFire0_DERIVwrt_ROSmax', BAperFire0_DERIVwrt_ROSmax
      __DEBUG3__(BA_DERIVwrt_ROSmax_c4,BA_DERIVwrt_ROSmax_c3,BA_DERIVwrt_ROSmax_dt)
!      __DEBUG3__(BA_DERIVwrt_ROSmax_tt,BA_DERIVwrt_ROSmax_ts,BA_DERIVwrt_ROSmax_et)
      __DEBUG2__(BA_DERIVwrt_ROSmax_tt,BA_DERIVwrt_ROSmax_et)   ! SSR20160211
      __DEBUG2__(BA_DERIVwrt_ROSmax_tshr,BA_DERIVwrt_ROSmax_tsav)   ! SSR20160211
      write(*,*) '###########################'
   endif
   
   ! Range checks
   call check_var_range(BA_DERIVwrt_alphaM, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_alphaM', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_IaParam1, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_IaParam1', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_IaParam2, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_IaParam2', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_AGBparam1, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_AGBparam1', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_AGBparam2, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_AGBparam2', lnd%time, FATAL)
   call check_var_range(RHderiv_X, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'RHderiv_X', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_RHparam1, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_RHparam1', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_RHparam2, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_RHparam2', lnd%time, FATAL)
   call check_var_range(THETAderiv_X, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'THETAderiv_X', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_thetaE, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_thetaE', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_THETAparam1, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_THETAparam1', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_THETAparam2, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_THETAparam2', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_popDsupp_NF_eps1, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_popDsupp_NF_eps1', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_popDsupp_NF_eps2, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_popDsupp_NF_eps2', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_popDsupp_NF_eps3, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_popDsupp_NF_eps3', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_fireDur_c4, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_fireDur_c4', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_fireDur_c3, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_fireDur_c3', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_fireDur_dt, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_fireDur_dt', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_fireDur_tt, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_fireDur_tt', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_fireDur_et, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_fireDur_et', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_ROSmax_c4, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_c4', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_ROSmax_c3, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_c3', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_ROSmax_dt, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_dt', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_ROSmax_tt, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_tt', lnd%time, FATAL)
!   call check_var_range(BA_DERIVwrt_ROSmax_ts, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_ts', lnd%time, FATAL)
   call check_var_range(BA_DERIVwrt_ROSmax_tshr, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_tshr', lnd%time, FATAL)   ! SSR20160211
   call check_var_range(BA_DERIVwrt_ROSmax_tsav, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_tsav', lnd%time, FATAL)   ! SSR20160211
   call check_var_range(BA_DERIVwrt_ROSmax_et, -10.**37., 10.**37., 'calc_fire_derivs_v2', 'BA_DERIVwrt_ROSmax_et', lnd%time, FATAL)

end subroutine calc_fire_derivs


subroutine update_fire_agb(vegn,soil)
   type(vegn_tile_type), intent(inout) :: vegn
   type(soil_tile_type), intent(in)    :: soil
   integer    ::    nv ! Number of vegetation cohorts
   real       ::    leafLitter_total_C, fineWoodLitter_total_C
   
   nv = vegn%n_cohorts
   
   ! When using Ensheng's multi-cohort model, the stats for each cohort are per individual, and
   ! then I would need to multiply by # individuals.
   vegn%fire_agb = sum(vegn%cohorts(1:nv)%bl    &
                      +vegn%cohorts(1:nv)%blv ) &
                   + agf_bs*sum(vegn%cohorts(1:nv)%bsw     &
                              + vegn%cohorts(1:nv)%bwood )
   if (soil_carbon_option==SOILC_CORPSE) then
      ! Calculate litter carbon, ignoring coarseWoodLitter, which should not contribute to spread
      call poolTotalCarbon(soil%leafLitter,totalCarbon=leafLitter_total_C)
      call poolTotalCarbon(soil%fineWoodLitter,totalCarbon=fineWoodLitter_total_C)
      vegn%fire_agb = vegn%fire_agb + leafLitter_total_C + fineWoodLitter_total_C
   endif
   
end subroutine update_fire_agb

subroutine fire_fragmentation(tile_list,land_area_m2,recalc,fragmenting_frac, burnable_frac)
   type(land_tile_list_type), intent(inout)  ::  tile_list
   real, intent(in)     ::  land_area_m2
   real, intent(inout)  ::  fragmenting_frac, burnable_frac
   logical, intent(in)  ::  recalc   ! If .TRUE. then actually recalculate. Otherwise just send diagnostics.
   type(land_tile_type), pointer :: ptr
   type(land_tile_enum_type) :: ts, te
   integer   ::   k   ! To keep track of tile number for printing watch_cell output.
   
   if (recalc) then
      call get_fragburn_fracs(tile_list,fragmenting_frac,burnable_frac)
      call check_var_range(fragmenting_frac, 0.0, 1.0,'fire_fragmentation', 'fragmenting_frac', lnd%time, FATAL)
      call check_var_range(burnable_frac,    0.0, 1.0,'fire_fragmentation', 'burnable_frac',    lnd%time, FATAL)
      call get_max_fire_size(tile_list,land_area_m2,fragmenting_frac,burnable_frac)
   endif
   
   if (is_watch_cell()) then
      write(*,*) '### checkpoint fire_fragmentation ###'
      write(*,*) 'fragmenting_frac', fragmenting_frac
      write(*,*) 'burnable_frac', burnable_frac
      write(*,*) 'land_area_m2', land_area_m2
      write(*,*) '#####################################'
   endif
   
   ts = first_elmt(tile_list) ; te=tail_elmt(tile_list) ; k = 0
   do while (ts /= te)
      ptr=>current_tile(ts); ts=next_elmt(ts)
      if (.not.associated(ptr%vegn)) cycle
      call send_tile_data(id_fragmenting_frac, fragmenting_frac, ptr%diag)
      call send_tile_data(id_burnable_frac, burnable_frac, ptr%diag)
   enddo
   
end subroutine fire_fragmentation


subroutine get_fragburn_fracs(tile_list,fragmenting_frac,burnable_frac)
   type(land_tile_list_type), intent(in)   ::   tile_list
   real, intent(out)  :: fragmenting_frac, burnable_frac
   type(land_tile_type), pointer :: ptr
   type(land_tile_enum_type) :: ts, te
   
   fragmenting_frac = 0.0
   burnable_frac = 0.0
   ts = first_elmt(tile_list) ; te=tail_elmt(tile_list)
   do while (ts /= te)
      ptr=>current_tile(ts); ts=next_elmt(ts)
      if (.not.associated(ptr%vegn)) then
         if (frag_incl_nonveg) fragmenting_frac = fragmenting_frac + ptr%frac
         cycle
      endif
      if (ptr%vegn%landuse==LU_CROP .OR. (frag_incl_PAST .AND. ptr%vegn%landuse==LU_PAST .AND. fire_option_past==FIRE_PASTFP)) then
         fragmenting_frac = fragmenting_frac + ptr%frac
      elseif  (ptr%vegn%landuse==LU_NTRL .OR. ptr%vegn%landuse==LU_SCND &
              .OR. (ptr%vegn%landuse==LU_PAST .AND. fire_option_past==FIRE_PASTLI)) then
         burnable_frac = burnable_frac + ptr%frac
      endif
   enddo
   
   ! Allow some wiggle room for burnable_frac
   if (burnable_frac>1.0 .AND. burnable_frac<1.000001) then
      burnable_frac = 1.0
   endif
   
end subroutine get_fragburn_fracs


subroutine get_max_fire_size(tile_list,lnd_area_m2,fragmenting_frac,burnable_frac)
   type(land_tile_list_type), intent(inout)   ::   tile_list
   real, intent(in)  :: lnd_area_m2
   real, intent(in)  :: fragmenting_frac, burnable_frac
   type(land_tile_type), pointer :: ptr
   type(land_tile_enum_type) :: ts, te
   real   ::   fnat   ! How to calculate the variable that Pfeiffer et al. (2013) call "fnat"
   
   ! Doesn't need to consider vegn%burned_frac because this only ever happens at the end
   ! of the day, AFTER resetting vegn%burned_frac to zero.
   
   ts = first_elmt(tile_list) ; te=tail_elmt(tile_list)
   do while (ts /= te)
      ptr=>current_tile(ts); ts=next_elmt(ts)
   
      if (.not.associated(ptr%vegn)) then
!         call send_tile_data(id_BAperFire_max, -1.0, ptr%diag)
         cycle
      endif
   
      ! This only needs to be done for tiles where fire size is actually calculated/used.
      if (ptr%vegn%landuse==LU_NTRL .OR. ptr%vegn%landuse==LU_SCND &
      .OR. (ptr%vegn%landuse==LU_PAST .AND. fire_option_past==FIRE_PASTLI)) then
         ! Based on Pfeiffer et al. (2013), eq. 33.
         ptr%vegn%max_fire_size = ((1.003 + exp(16.607-41.503*burnable_frac))**-2.169)*(ptr%frac*lnd_area_m2)
         if (is_watch_cell()) then
            write(*,*) '### get_max_fire_size output ###'
            if (is_watch_point()) then
               write(*,*) '(This is watch_point.)'
            endif
            write(*,*) 'fragmenting_frac', fragmenting_frac
            write(*,*) 'burnable_frac', burnable_frac
            write(*,*) 'lnd_area_m2', lnd_area_m2
            write(*,*) 'ptr%frac', ptr%frac
            write(*,*) 'fnat', fnat
            write(*,*) 'ptr%vegn%max_fire_size', ptr%vegn%max_fire_size
            write(*,*) 'max_fire_size_min', max_fire_size_min
            write(*,*) '################################'
         endif
         if (ptr%vegn%max_fire_size < max_fire_size_min) then
            ptr%vegn%max_fire_size = max_fire_size_min
         endif
      else
         ptr%vegn%max_fire_size = -1.0
         if (is_watch_cell()) then
            write(*,*) '### get_max_fire_size output ###'
            if (is_watch_point()) then
               write(*,*) '(This is watch_point.)'
            endif
            write(*,*) 'ptr%vegn%max_fire_size', ptr%vegn%max_fire_size
            write(*,*) '################################'
         endif
      endif
      
      
   enddo
end subroutine get_max_fire_size


subroutine defo_annCalcs(fireFtheta_ann_acm,nmn_acm, &
                                 fireFtheta_ann,defoFS)
   real, intent(in)    ::   fireFtheta_ann_acm
   integer, intent(in) ::   nmn_acm
   real, intent(out)   ::   fireFtheta_ann, defoFS
   
   fireFtheta_ann  = fireFtheta_ann_acm/nmn_acm
   defoFS = fs_min + max(0.0,min(1.0,(fireFtheta_ann-fp_lo)/(fp_hi-fp_lo)))*(fs_max-fs_min)

end subroutine defo_annCalcs


end module vegn_fire_mod
