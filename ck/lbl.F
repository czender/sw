c     $Id$

      program lbl

#if ( defined DOUBLE )
#define COMPUTATIONAL_PRECISION double precision
c#define NF_REAL NF_DOUBLE
#else
#define COMPUTATIONAL_PRECISION real
c#define NF_DOUBLE NF_REAL
#endif

c     Purpose: Line-by-line RT model employing HITRAN line information.
c     Monochromatic optical depth from gaseous line absorption is combined  
c     with continuum absorption and scattering cross sections and fed to DISORT. 
c     All I/O uses netCDF interface.

c     NB: DISORT uses fortran unit 99 for its own purposes.

c     Note on nomenclature:
c     lbl makes a distinction betweeen lines (ln), bands (bnd), intervals (ntv), and regions (rgn).
c     Thus there are four different wavelength/wavenumber grids in lbl.
c     The line is the fundamental unit of lbl.
c     lbl stores line parameters from HITRAN in arrays like ln_ctr(), ln_str(), and HWHM_air(). 
c     DISORT works with intervals, which in general may be composed of many lines (up to 10000). 
c     In lbl, however, intervals must be smaller than the HWHM of a particular line in order to resolve the line.
c     Furthermore, the spacing (resolution) of intervals varies from line to line, 
c     and depends on the line parameters themselves.
c     The average wavenumber resolution of all the intervals in the SW or LW is order .001 cm-1.
c     The RT computations in lbl scale as the number of intervals.
c     Because there are so many intervals, lbl groups them together in bands for output. 
c     Bands are just summations of one or more intervals.
c     The boundaries of bands are normally chosen to represent physically meaningful bands,
c     such as the 15 um band system of CO2, or the 9.6 um band system of O3.
c     The region is the entire spectral region for the RT calculation.
c     For example, the region might be the SW (i.e., 2000--20000 cm-1), the LW, or both.
c     The code ensures that band boundaries coincide with interval boundaries.
c     The code ensures that region boundaries coincide with band boundaries.

c     If wvn_min_rgn and wvn_max_rgn are not specified on the command line then they are obtained from the HITRAN line data file.

c     A main purpose of lbl is to provide a benchmark against which to parameterize c-k models.
c     Thus the bands in lbl are usually chosen to correspond to bands in the c-k model.
c     The output of lbl includes parameterized k-distributions intended for use in c-k models.

c     Key band indices:

c     Debugging usage:
c     Perform full setup for all bands but only call DISORT() for #1593:
c     lbl -D 1 -E -e 1593 -d foo.nc 

c     Production usage:

c     Abbreviations:
c     htrn = HITRAN
c     rtl = rotational
c     dpd = dependence
c     tpt = temperature
c     xpn = exponent

c     Code proceeds as follows:
c     Section 1: Initialization
c     Section 2: Main computation loop over spectral bands
c     Section 3: Postprocessing

c     Begin section 1: Initialization
c     IEEE setup
c     initialization of program defaults
c     processing of command line switches
c     netCDF input of atmospheric thermodynamic and cloud profile
c     netCDF input of gaseous absorption and aerosol scattering and absorption characteristics
c     override netCDF input data according to command line switches
c     band-independent DISORT() initialization
c     band loops for figuring out where to splice input data together
c     initialization of  band quantities that are level independent
c     initialization of level quantities that are band independent
c     a band-level loop for initializing band-level quantities

      implicit none
c     Parameters
      character*(*) 
     $     rcs_Id
      integer
     $     max_ln_nbr,
     $     max_ntv_nbr,
     $     max_nbr_cst,
     $     max_nbr_htrn_bnd,
     $     max_nbr_azi,         ! # user azimuthal angles
     $     bnd_nbr_max,
     $     max_nbr_lev,
     $     max_nbr_levp,
     $     max_nbr_nst_bnd,
     $     max_nbr_plr,         ! # user polar angles
     $     max_nbr_str,         ! # computational polar angles
     $     max_nbr_tau          ! # user optical depths

      parameter(max_ln_nbr=44228, ! O3 from .2-5 microns
     $     max_ntv_nbr=10000,   ! # of line in region
     $     max_nbr_cst=31,      ! # of gases in HITRAN96
     $     max_nbr_htrn_bnd=1590, ! 
     $     max_nbr_azi=16,      ! arbitrary
     $     bnd_nbr_max=1690,    ! WMO spliced onto 10 cm-1 resolution from .56--5 microns
     $     max_nbr_lev=110,     ! roughly 10 mb resolution from 1010 mb to TOA
     $     max_nbr_levp=max_nbr_lev+1,
     $     max_nbr_nst_bnd=69,  ! TSBR resolution
     $     max_nbr_plr=16,      ! arbitrary
     $     max_nbr_str=16,      ! arbitrary
     $     max_nbr_tau=max_nbr_lev+1,
     $     rcs_Id='$Id$')

c     Commons
#include <fnd_cst.com>
#include <hitran.com>
#include <drv_cst.com>
#include <netcdf.inc>
c     Input Arguments
c     Input/Output Arguments
c     Output Arguments
c     Local workspace
#include <time_ini.com> /* Initialize system-dependent date/time functions */
#include <ieee_ini.com> /* Initialize system-dependent IEEE functions */

      character
     $     argv*80,
     $     aer_sng*100,
     $     azi_sng*100,
     $     bga_sng*100,
     $     char_foo*1,
     $     fl_H2OH2O*80,
     $     fl_htrn*80,
     $     fl_NO2*80,
     $     fl_O2*80,
     $     fl_O3*80,
     $     fl_O2O2*80,
     $     fl_aer*80,
     $     fl_bga*80,
     $     fl_clm*80,
     $     fl_ice*80,
     $     fl_liq*80,
     $     fl_nst*80,
     $     fl_out*80,
     $     fl_slr*80,
     $     sng*80,
     $     opt_dep_sng*100,
     $     plr_sng*100,
     $     prf_sng*100,
     $     rcs_Date*28,
     $     rcs_Revision*16,
     $     str_sng*100,
     $     stt_CO2*100,
     $     stt_H2O*100,
     $     stt_H2OH2O*100,
     $     stt_Herzberg*100,
     $     stt_NO2*100,
     $     stt_O2*100,
     $     stt_O3*100,
     $     stt_O2O2*100,
     $     stt_OH*100,
     $     stt_Planck*100,
     $     stt_Rayleigh*100,
     $     stt_aer*100,
     $     stt_bga*100,
     $     stt_ice*100,
     $     stt_liq*100,
     $     stt_slr*100,
     $     stt_flt_nst*100,
     $     stt_top_lvl*100

      integer
     $     arg,
     $     dbg_lvl,
     $     bnd_dbg,
     $     exit_status,          ! program exit status
     $     int_foo,
     $     nbr_arg,
     $     rcd                  ! return success code

      logical
     $     cmd_ln_mpc_CWP,
     $     cmd_ln_odxc_obs_aer,
     $     cmd_ln_alb,
     $     cmd_ln_slr_zen_ngl_cos,
     $     cmd_ln_slr_cst,
     $     flg_CO2,
     $     flg_H2OH2O,
     $     flg_H2O,
     $     flg_Herzberg,
     $     flg_NO2,
     $     flg_O2,
     $     flg_O3,
     $     flg_O2O2,
     $     flg_OH,
     $     flg_Planck,
     $     flg_Rayleigh,
     $     flg_aer,
     $     flg_bga,
     $     flg_ice,
     $     flg_liq,
     $     flt_nst,
     $     force_ice_phz,
     $     force_liq_phz,
     $     logical_foo,
     $     single_bnd_computation,
     $     sv_cmp_opt_dep,
     $     sv_cmp_plr_ngl,
     $     sv_ntn,
     $     top_lvl,
     $     tst_case_HG,
     $     tst_case_Rayleigh

      integer
     $     bnd_nbr,             ! dimension size
     $     nbr_azi,             ! dimension size
     $     nbr_bnd,             ! dimension size
     $     nbr_lev,             ! dimension size
     $     nbr_levp,            ! dimension size
     $     ln_nbr,              ! dimension size
     $     nbr_mom,             ! dimension size
     $     nbr_nst_bnd,         ! dimension size
     $     nbr_plr,             ! dimension size
     $     nbr_str,             ! dimension size
     $     nbr_tau              ! dimension size
      
      integer
     $     azi_idx,             ! counting index
     $     bnd_idx,             ! counting index
     $     bnd_idx_nst,         ! counting index
     $     lev_idx,             ! counting index
     $     levp_idx,            ! counting index
     $     ln_idx,              ! counting index
     $     mom_idx,             ! counting index
     $     plr_idx,             ! counting index
     $     str_idx,             ! counting index
     $     tau_idx              ! counting index
      
      integer
     $     azi_dim_id,          ! dimension ID for azi
     $     bnd_dim_id,          ! dimension ID for bnd
     $     cnt_azi_plr_bnd_levp(4),
     $     cnt_azi_plr_levp(3),
     $     cnt_bnd_lev(2),
     $     cnt_bnd_levp(2),
     $     cnt_plr_bnd(2),
     $     cnt_plr_bnd_levp(3),
     $     cnt_plr_levp(2),
     $     dim_azi_plr_bnd_levp(4),
     $     dim_azi_plr_levp(3),
     $     dim_bnd_lev(2),
     $     dim_bnd_levp(2),
     $     dim_plr_bnd(2),
     $     dim_plr_bnd_levp(3),
     $     dim_plr_levp(2),
     $     lev_dim_id,          ! dimension ID for lev
     $     levp_dim_id,         ! dimension ID for levp
     $     ln_ctr_dim_id,       ! dimension ID for ln
     $     nc_id,               ! file handle
     $     plr_dim_id,          ! dimension ID for plr
     $     srt(4),              ! starting offsets (long *) 
     $     tau_dim_id           ! dimension ID for tau
      
c     netCDF output variables not contained in input files
      integer
     $     abs_bb_SAS_id,
     $     abs_bb_atm_id,
     $     abs_bb_sfc_id,
     $     abs_nst_SAS_id,
     $     abs_nst_atm_id,
     $     abs_nst_sfc_id,
     $     abs_spc_SAS_id,
     $     abs_spc_atm_id,
     $     abs_spc_sfc_id,
     $     alb_sfc_id,
     $     azi_dgr_id,          
     $     azi_id,              ! coordinate ID
     $     bnd_id,              ! coordinate ID
     $     flx_abs_atm_rdr_id,
     $     flx_bb_abs_atm_id,
     $     flx_bb_abs_id,
     $     flx_bb_abs_sfc_id,
     $     flx_bb_abs_ttl_id,
     $     flx_bb_dwn_TOA_id,
     $     flx_bb_dwn_dff_id,
     $     flx_bb_dwn_drc_id,
     $     flx_bb_dwn_id,
     $     flx_bb_dwn_sfc_id,
     $     flx_bb_net_id,
     $     flx_bb_up_id,
     $     flx_nst_abs_atm_id,
     $     flx_nst_abs_id,
     $     flx_nst_abs_sfc_id,
     $     flx_nst_abs_ttl_id,
     $     flx_nst_dwn_TOA_id,
     $     flx_nst_dwn_id,
     $     flx_nst_dwn_sfc_id,
     $     flx_nst_net_id,
     $     flx_nst_up_id,
     $     flx_spc_abs_SAS_id,
     $     flx_spc_abs_atm_id,
     $     flx_spc_abs_id,
     $     flx_spc_abs_sfc_id,
     $     flx_spc_dwn_TOA_id,
     $     flx_spc_dwn_dff_id,
     $     flx_spc_dwn_drc_id,
     $     flx_spc_dwn_id,
     $     flx_spc_dwn_sfc_id,
     $     flx_spc_up_id,
     $     flx_slr_frc_id,
     $     htg_rate_bb_id,
     $     lev_id,              ! coordinate ID
     $     levp_id,             ! coordinate ID
     $     ntn_bb_aa_id,
     $     ntn_bb_mean_id,
     $     ntn_spc_aa_ndr_id,
     $     ntn_spc_aa_ndr_sfc_id,
     $     ntn_spc_aa_sfc_id,
     $     ntn_spc_aa_zen_id,
     $     ntn_spc_aa_zen_sfc_id,
     $     ntn_spc_chn_id,
     $     ntn_spc_mean_id,
     $     odxc_spc_ttl_id,
     $     plr_cos_id,          
     $     plr_dgr_id,          
     $     plr_id,              ! coordinate ID
     $     rfl_bb_SAS_id,
     $     rfl_bb_sfc_id,
     $     rfl_nst_SAS_id,
     $     rfl_nst_sfc_id,
     $     rfl_spc_SAS_id,
     $     tau_id,              ! coordinate ID
     $     tau_prs_id,          
     $     trn_bb_atm_id,
     $     trn_nst_atm_id,
     $     trn_spc_atm_ttl_id,
     $     wvl_ctr_id,
     $     wvl_max_id,
     $     wvl_min_id,
     $     wvl_sz_id,
     $     wvn_ctr_id,
     $     wvn_max_id,
     $     wvn_min_id,
     $     wvn_sz_id,
     $     z_id,
     $     z_int_id

c     Instrument input variables
      integer
     $     nst_SRF_id,
     $     wvl_ctr_nst_id,
     $     wvl_max_nst_id,
     $     wvl_min_nst_id

c     CLM input variables
      integer
     $     RH_liq_id,
     $     alb_NIR_str_zen_ngl_id,
     $     alb_NIR_weak_zen_ngl_id,
     $     alb_vis_str_zen_ngl_id,
     $     alb_vis_weak_zen_ngl_id,
     $     slr_zen_ngl_cos_id,
     $     prs_dlt_id,
     $     frc_ice_id,
     $     frc_ice_ttl_id,
     $     grv_id,
     $     lat_dgr_id,
     $     lcl_time_hr_id,
     $     lcl_yr_day_id,
     $     mmw_mst_air_id,
     $     mpc_CWP_id,
     $     mpl_CO2_id,
     $     mpl_CWP_id,
     $     mpl_H2O_id,
     $     mpl_IWP_id,
     $     mpl_LWP_id,
     $     mpl_O2_id,
     $     mpl_OH_id,
     $     mpl_aer_id,
     $     mpl_bga_id,
     $     mpl_mst_air_id,
     $     npl_NO2_id,
     $     npl_O2_id,
     $     npl_O3_id,
     $     npl_O2O2_id,
     $     odxc_obs_aer_id,
     $     odxc_obs_bga_id,
     $     odxl_obs_aer_id,
     $     odxl_obs_bga_id,
     $     p_id,
     $     p_int_id,
     $     q_CO2_id,
     $     q_H2O_id,
     $     q_O2_id,
     $     q_OH_id,
     $     spec_heat_mst_air_id,
     $     t_id,
     $     t_int_id,
     $     t_skn_id,
     $     xnt_fac_id,
     $     wvl_obs_aer_id,
     $     wvl_obs_bga_id,
     $     z_cld_btm_id,
     $     z_cld_thick_id

c     netCDF output variables
      real
     $     abs_bb_SAS,
     $     abs_bb_atm,
     $     abs_bb_sfc,
     $     abs_nst_SAS,
     $     abs_nst_atm,
     $     abs_nst_sfc,
     $     abs_spc_SAS(bnd_nbr_max),
     $     abs_spc_atm(bnd_nbr_max),
     $     abs_spc_sfc(bnd_nbr_max),
     $     alb_sfc,
     $     asm_prm_Mie_ttl(bnd_nbr_max,max_nbr_lev),
     $     azi(max_nbr_azi),    ! coordinate variable
     $     azi_dgr(max_nbr_azi),
     $     bnd(bnd_nbr_max),  ! coordinate variable
     $     flx_abs_atm_rdr(bnd_nbr_max),
     $     flx_bb_abs(max_nbr_lev),
     $     flx_bb_abs_atm,
     $     flx_bb_abs_sfc,
     $     flx_bb_abs_ttl,
     $     flx_bb_dwn(max_nbr_levp),
     $     flx_bb_dwn_TOA,
     $     flx_bb_dwn_dff(max_nbr_levp),
     $     flx_bb_dwn_drc(max_nbr_levp),
     $     flx_bb_dwn_sfc,
     $     flx_bb_net(max_nbr_levp),
     $     flx_bb_up(max_nbr_levp),
     $     flx_nst_abs(max_nbr_lev),
     $     flx_nst_abs_atm,
     $     flx_nst_abs_sfc,
     $     flx_nst_abs_ttl,
     $     flx_nst_dwn(max_nbr_levp),
     $     flx_nst_dwn_TOA,
     $     flx_nst_dwn_sfc,
     $     flx_nst_net(max_nbr_levp),
     $     flx_nst_up(max_nbr_levp),
     $     flx_spc_abs(bnd_nbr_max,max_nbr_lev),
     $     flx_spc_abs_SAS(bnd_nbr_max),
     $     flx_spc_abs_atm(bnd_nbr_max),
     $     flx_spc_abs_sfc(bnd_nbr_max),
     $     flx_spc_dwn(bnd_nbr_max,max_nbr_levp),
     $     flx_spc_dwn_TOA(bnd_nbr_max),
     $     flx_spc_dwn_dff(bnd_nbr_max,max_nbr_levp),
     $     flx_spc_dwn_drc(bnd_nbr_max,max_nbr_levp),
     $     flx_spc_dwn_sfc(bnd_nbr_max),
     $     flx_spc_up(bnd_nbr_max,max_nbr_levp),
     $     flx_slr_frc(bnd_nbr_max),
     $     htg_rate_bb(max_nbr_lev),
     $     lev(max_nbr_lev),    ! coordinate variable
     $     levp(max_nbr_levp),  ! coordinate variable
     $     ntn_bb_aa(max_nbr_plr,max_nbr_levp),
     $     ntn_bb_mean(max_nbr_lev),
     $     ntn_spc_aa_ndr(bnd_nbr_max,max_nbr_levp),
     $     ntn_spc_aa_ndr_sfc(bnd_nbr_max),
     $     ntn_spc_aa_sfc(max_nbr_plr,bnd_nbr_max),
     $     ntn_spc_aa_zen(bnd_nbr_max,max_nbr_levp),
     $     ntn_spc_aa_zen_sfc(bnd_nbr_max),
     $     ntn_spc_chn(max_nbr_azi,max_nbr_plr,max_nbr_levp),
     $     ntn_spc_mean(bnd_nbr_max,max_nbr_lev),
     $     odxc_spc_ttl(bnd_nbr_max),
     $     odxl_spc_ttl(bnd_nbr_max,max_nbr_lev),
     $     plr(max_nbr_plr),    ! coordinate variable
     $     plr_cos(max_nbr_plr),
     $     plr_dgr(max_nbr_plr),
     $     rfl_bb_SAS,
     $     rfl_bb_sfc,
     $     rfl_nst_SAS,
     $     rfl_nst_sfc,
     $     rfl_spc_SAS(bnd_nbr_max),
     $     ss_alb_fct(bnd_nbr_max,max_nbr_lev),
     $     tau(max_nbr_tau),    ! coordinate variable
     $     tau_prs(max_nbr_tau),
     $     trn_bb_atm,
     $     trn_nst_atm,
     $     trn_spc_atm_ttl(bnd_nbr_max),
     $     wvl_ctr(bnd_nbr_max),
     $     wvl_max(bnd_nbr_max),
     $     wvl_min(bnd_nbr_max),
     $     wvl_sz(bnd_nbr_max),
     $     wvn_ctr(bnd_nbr_max),
     $     wvn_max(bnd_nbr_max),
     $     wvn_min(bnd_nbr_max),
     $     wvn_sz(bnd_nbr_max)

c     Instrument input variables
      real
     $     nst_SRF(max_nbr_nst_bnd),
     $     wvl_ctr_nst(max_nbr_nst_bnd),
     $     wvl_max_nst(max_nbr_nst_bnd),
     $     wvl_min_nst(max_nbr_nst_bnd)

c     CLM input variables
      double precision
     $     lat_dgr,
     $     lcl_time_hr,
     $     lcl_yr_day,
     $     slr_zen_ngl_cos
      real
     $     RH_liq(max_nbr_lev),
     $     alb_NIR_str_zen_ngl,
     $     alb_NIR_weak_zen_ngl,
     $     alb_vis_str_zen_ngl,
     $     alb_vis_weak_zen_ngl,
     $     prs_dlt(max_nbr_lev),
     $     frc_ice(max_nbr_lev),
     $     frc_ice_ttl,
     $     grv(max_nbr_lev),
     $     mmw_mst_air(max_nbr_lev),
     $     mpc_CWP,
     $     mpl_CO2(max_nbr_lev),
     $     mpl_CWP(max_nbr_lev),
     $     mpl_H2O(max_nbr_lev),
     $     mpl_IWP(max_nbr_lev),
     $     mpl_LWP(max_nbr_lev),
     $     mpl_O2(max_nbr_lev),
     $     mpl_OH(max_nbr_lev),
     $     mpl_aer(max_nbr_lev),
     $     mpl_bga(max_nbr_lev),
     $     mpl_mst_air(max_nbr_lev),
     $     npl_NO2(max_nbr_lev),
     $     npl_O2(max_nbr_lev),
     $     npl_O3(max_nbr_lev),
     $     npl_O2O2(max_nbr_lev),
     $     odxc_obs_aer,
     $     odxc_obs_bga,
     $     odxl_obs_aer(max_nbr_lev),
     $     odxl_obs_bga(max_nbr_lev),
     $     p(max_nbr_lev),
     $     p_int(max_nbr_levp),
     $     q_CO2(max_nbr_lev),
     $     q_H2O(max_nbr_lev),
     $     q_O2(max_nbr_lev),
     $     q_OH(max_nbr_lev),
     $     spec_heat_mst_air(max_nbr_lev),
     $     t(max_nbr_lev),
     $     t_int(max_nbr_levp),
     $     t_skn,
     $     wvl_obs_aer,
     $     wvl_obs_bga,
     $     xnt_fac,
     $     z(max_nbr_lev),
     $     z_cld_btm,
     $     z_cld_thick,
     $     z_int(max_nbr_levp)

c     Local arrays
      real
     $     asm_prm_Mie_aer(bnd_nbr_max),
     $     asm_prm_Mie_bga(bnd_nbr_max),
     $     asm_prm_Mie_ice(bnd_nbr_max),
     $     asm_prm_Mie_liq(bnd_nbr_max),
     $     bnd_wgt(bnd_nbr_max),
     $     flx_spc_net(bnd_nbr_max,max_nbr_levp),
     $     idx_rfr_air_STP(bnd_nbr_max),
     $     odal_spc_ttl(bnd_nbr_max,max_nbr_lev),
     $     odsl_spc_ttl(bnd_nbr_max,max_nbr_lev),
     $     p_bar(max_nbr_levp),
     $     phi_wgt(max_nbr_lev),
     $     psi_wgt(max_nbr_lev),
     $     sca_cff_mss_Ray(max_nbr_lev),
     $     sca_frc_Mie(max_nbr_lev),
     $     u_bar(max_nbr_levp)

c     Local scalars
      double precision
     $     slr_zen_ngl_cos_cmd_ln
      integer
     $     bnd_obs_aer,
     $     bnd_obs_bga,
     $     aer_lvl_nbr,
     $     cld_lvl_nbr
      real
     $     alb_cmd_ln,
     $     float_foo,
     $     mpc_CWP_cmd_ln,
     $     mpc_IWP,
     $     odxc_obs_aer_cmd_ln,
     $     pi,
     $     slr_cst,
     $     slr_cst_cmd_ln,
     $     slr_cst_xnt_fac,
     $     wvl_Planck

c     LBL variables
c     HITRAN line parameter input variables
      integer
     $     HWHM_air_id,
     $     HWHM_tpt_dpn_xpn_id,
     $     iso_idx,
     $     ln_ctr_id,
     $     ln_iso_id,
     $     ln_mlc_id,
     $     ln_nrg_lwr_id,
     $     ln_str_id,
     $     mlc_idx

      integer*2                 ! NF_SHORT
     $     iso_nbr,
     $     ln_iso(max_ln_nbr),
     $     ln_mlc(max_ln_nbr),
     $     mlc_nbr

      COMPUTATIONAL_PRECISION
     $     HWHM_air(max_ln_nbr),
     $     HWHM_tpt_dpn_xpn(max_ln_nbr),
     $     ln_ctr(max_ln_nbr),
     $     ln_nrg_lwr(max_ln_nbr),
     $     ln_str(max_ln_nbr)

      external Voigt_Hum82     
      complex Voigt_Hum82     

      complex Voigt_arg_z       ! Argument to complex error function for Voigt profiles

      integer ln_nbr_ngh        ! Number of lines with wng_sz_wvn of current spectral interval
      integer ln_idx_end        ! Index (into ln array) of righmost neighbor 
      integer ln_idx_srt        ! Index (into ln array) of leftmost neighbor 
      integer map(iso_per_mlc_nbr_max_htrn,mlc_nbr_max_htrn) ! map of mlc,iso indices to isopotomer index

      double precision HWHM_Lorentz ! HWHM of pressure broadened line
      double precision HWeM_Doppler ! HWHM of Doppler broadened line (half-width at e-1 of peak)
      double precision ln_dst_ctr ! Distance of current spectral interval to line center
      double precision ln_prf   ! Normalized profile of spectral line in current spectral interval
      double precision wng_sz_wvn ! cm-1 Distance from line center to wing tips
      double precision rsn_wvn  ! cm-1 Resolution of spectral grid in wavenumbers
      double precision sqrt_ln2 ! Constant for Voigt function
      double precision sqrt_rcp_pi ! Constant for Voigt function
      double precision sqrt_ln2_rcp_pi ! Constant for Voigt function
      double precision t0_rcp_t(max_nbr_lev) ! temperature ratio for Lorentz HWHM computation
      double precision p_rcp_p0(max_nbr_lev) ! pressure ratio for Lorentz HWHM computation
      double precision wvn_max_rgn ! cm-1 Maximum wavenumber in region
      double precision wvn_min_rgn ! cm-1 Minimum wavenumber in region
      double precision Voigt_arg_x ! Argument to complex error function for Voigt profiles
      double precision Voigt_arg_y ! Argument to complex error function for Voigt profiles
      double precision Voigt_fac ! Normalization factor for complex error function for Voigt profiles
      double precision xpn_mlc(mlc_nbr_max_htrn) ! Exponent defining temperature dependence of rotational partition function
      double precision mmw_iso(iso_nbr_max_htrn) ! kg mol-1 mean molecular weight of all HITRAN isopotomers
      
c     The following line is from the DISORT() subroutine. 
c     parameter ( mxcly = 85, mxulv = 85, mxcmu = 49, mxumu = 16,
c     $     maxpphi = 3)
c     The parameters declared below should match the values from
c     the DISORT() subroutine.
      integer 
     $     maxcly,
     $     maxcmu,
     $     maxphi,
     $     maxulv,
     $     maxumu
      parameter(maxcly=max_nbr_tau,
     $     maxcmu = max_nbr_str,
     $     maxphi = max_nbr_azi,
     $     maxulv = max_nbr_tau,
     $     maxumu = max_nbr_plr)

c     DISORT() input variables:
      character  header*127
      logical  deltam, lamber, plank, onlyfl, prnt(7), usrang, usrtau
      integer  ibcnd, nlyr, numu, nstr, nphi, ntau
      real     accur, albedo, btemp, dtauc( maxcly ), fbeam, fisot,
     $     hl( 0:maxcmu ), phi( maxphi ), pmom( 0:maxcmu, maxcly ),
     $     phi0, ssalb( maxcly ), temper( 0:maxcly ), temis, ttemp,
     $     wvnmlo, wvnmhi, umu( maxumu ), umu0, utau( maxulv )

c     DISORT() output variables:
      real     rfldir( maxulv ), rfldn( maxulv ), flup( maxulv ),
     $     dfdt( maxulv ), uavg( maxulv ), u0u( maxumu, maxulv ),
     $     uu( maxumu, maxulv, maxphi ), albmed( maxumu ),
     $     trnmed( maxumu )

c     Externals
      integer 
     $     iargc
c     Main code

#include <time_srt.com> /* Set timers, print date */
#include <ieee_srt.com> /* Set IEEE handler options */

c     Initialize default values
#if ( defined RS6K )
      fl_H2OH2O='odal_H2OH2O_Chy97()'
      fl_htrn='/usr/tmp/zender/lbl_CO2_1931_1939.nc'
      fl_O2='/usr/tmp/zender/swnb_O2.nc'
      fl_O3='/usr/tmp/zender/abs_xsx_O3.nc'
      fl_O2O2='/usr/tmp/zender/abs_xsx_O2O2.nc'
      fl_NO2='/usr/tmp/zender/abs_xsx_NO2.nc'
      fl_clm='/usr/tmp/zender/mls_cld.nc'
      fl_ice='/usr/tmp/zender/ice_20.nc'
      fl_liq='/usr/tmp/zender/liq_10.nc'
      fl_nst='/usr/tmp/zender/nst_FSBR.nc'
      fl_aer='/usr/tmp/zender/aer_mineral_dust.nc'
      fl_bga='/usr/tmp/zender/aer_h2so4_215K.nc'
      fl_out='/usr/tmp/zender/dsrt.nc'
      fl_slr='/usr/tmp/zender/spc_Kur95_01wvn.nc'
#else 
      fl_H2OH2O='odal_H2OH2O_Chy97()'
      fl_htrn='/data2/zender/hitran/lbl_CO2_1931_1939.nc'
      fl_O2='/data/zender/aca/swnb_O2.nc'
      fl_O3='/data/zender/aca/abs_xsx_O3.nc'
      fl_O2O2='/data/zender/aca/abs_xsx_O2O2.nc'
      fl_NO2='/data/zender/aca/abs_xsx_NO2.nc'
      fl_clm='/data/zender/aca/mls_cld.nc'
      fl_ice='/data/zender/aca/ice_20.nc'
      fl_liq='/data/zender/aca/liq_10.nc'
      fl_nst='/data/zender/aca/nst_FSBR.nc'
      fl_aer='/data/zender/aca/aer_mineral_dust.nc'
      fl_bga='/data/zender/aca/aer_h2so4_215K.nc'
      fl_out='/data/zender/aca/dsrt.nc'
      fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
#endif
      bnd_dbg=1603
      cmd_ln_mpc_CWP=.false.
      cmd_ln_odxc_obs_aer=.false.
      cmd_ln_alb=.false.
      cmd_ln_slr_zen_ngl_cos=.false.
      cmd_ln_slr_cst=.false.
      dbg_lvl=0
      exit_status=0
      flg_CO2=.false.
      flg_H2OH2O=.false.
      flg_H2O=.true.
      flg_Herzberg=.false.
      flg_NO2=.false.
      flg_O2=.false.
      flg_O3=.false.
      flg_O2O2=.false.
      flg_OH=.false.
      flg_Planck=.false.
      flg_Rayleigh=.false.
      flg_aer=.false.
      flg_bga=.false.
      flg_ice=.false.
      flg_liq=.false.
      float_foo=0.
      flt_nst=.true.
      force_ice_phz=.false.
      force_liq_phz=.false.
      logical_foo=.true.
      nbr_azi=max_nbr_azi
      nbr_H2OH2O_bnd=max_nbr_H2OH2O_bnd
      nbr_plr=2
      nbr_str=4
      pi=4.*atan(1.)
      rcs_Date='$Date$'
      rcs_Revision='$Revision$'
      rsn_wvn=0.01              ! cm-1
      single_bnd_computation=.false.
      slr_cst=slr_cst_CCM
      sqrt_ln2=sqrt(log(2.))
      sqrt_ln2_rcp_pi=sqrt(log(2.)/pi)
      sqrt_rcp_pi=sqrt(1./pi)
      srt(1)=1
      sv_cmp_opt_dep=.true.
      sv_cmp_plr_ngl=.true.
      sv_ntn=.false.
      top_lvl=.false.
      tst_case_HG=.false.
      tst_case_Rayleigh=.false.
      wvl_Planck=2.e-6
      wvn_max_rgn=1939.
      wvn_min_rgn=1931.

c     Retrieve command line arguments
      write (6,'(a)') 'lbl.F version ' // rcs_Revision(11:15) // 'dated ' // rcs_Date(8:26) // ' GMT' 
      nbr_arg=iargc()
      do arg=1,nbr_arg
         call getarg(arg,argv)

         if (argv(1:2).eq.'-A') then
            flg_aer=.not.flg_aer
         endif
         if (argv(1:2).eq.'-a') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_aer
         endif
         if (argv(1:2).eq.'-B') then
            flg_bga=.not.flg_bga
         endif
         if (argv(1:2).eq.'-b') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_bga
         endif
         if (argv(1:2).eq.'-C') then
            flg_CO2=.not.flg_CO2
         endif
         if (argv(1:2).eq.'-c') then
         endif
         if (argv(1:2).eq.'-D') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') dbg_lvl
         endif
         if (argv(1:2).eq.'-d') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_out
         endif
         if (argv(1:2).eq.'-E') then
            single_bnd_computation=.not.single_bnd_computation
         endif
         if (argv(1:2).eq.'-e') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') bnd_dbg
         endif
         if (argv(1:2).eq.'-F') then
            force_ice_phz=.not.force_ice_phz
         endif
         if (argv(1:2).eq.'-f') then
            force_liq_phz=.not.force_liq_phz
         endif
         if (argv(1:2).eq.'-H') then
            flg_H2O=.not.flg_H2O
         endif
         if (argv(1:2).eq.'-h') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_htrn
         endif
         if (argv(1:2).eq.'-I') then
            flg_ice=.not.flg_ice
         endif
         if (argv(1:2).eq.'-i') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_ice
         endif
         if (argv(1:2).eq.'-J') then
            flg_Planck=.not.flg_Planck
         endif
         if (argv(1:2).eq.'-j') then
         endif
         if (argv(1:2).eq.'-K') then
            flg_O2O2=.not.flg_O2O2
         endif
         if (argv(1:2).eq.'-k') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O2O2
         endif
         if (argv(1:2).eq.'-L') then
            flg_liq=.not.flg_liq
         endif
         if (argv(1:2).eq.'-l') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_liq
         endif
         if (argv(1:2).eq.'-M') then
            cmd_ln_odxc_obs_aer=.not.cmd_ln_odxc_obs_aer
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') odxc_obs_aer_cmd_ln
         endif
         if (argv(1:2).eq.'-m') then
            cmd_ln_mpc_CWP=.not.cmd_ln_mpc_CWP
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') mpc_CWP_cmd_ln
         endif
         if (argv(1:2).eq.'-N') then
            flt_nst=.not.flt_nst
         endif
         if (argv(1:2).eq.'-n') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_nst
         endif
         if (argv(1:2).eq.'-O') then
            flg_O2=.not.flg_O2
         endif
         if (argv(1:2).eq.'-o') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O2
         endif
         if (argv(1:2).eq.'-P') then
            tst_case_HG=.not.tst_case_HG
         endif
         if (argv(1:2).eq.'-p') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_clm
         endif
         if (argv(1:2).eq.'-Q') then
            flg_H2OH2O=.not.flg_H2OH2O
         endif
         if (argv(1:2).eq.'-q') then
         endif
         if (argv(1:2).eq.'-R') then
            flg_Rayleigh=.not.flg_Rayleigh
         endif
         if (argv(1:2).eq.'-r') then
            cmd_ln_alb=.not.cmd_ln_alb
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') alb_cmd_ln
         endif
         if (argv(1:2).eq.'-r') then
         endif
         if (argv(1:2).eq.'-S') then
            cmd_ln_slr_cst=.not.cmd_ln_slr_cst
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') slr_cst_cmd_ln
         endif
         if (argv(1:2).eq.'-s') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') nbr_str
         endif
         if (argv(1:2).eq.'-T') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_slr
         endif
         if (argv(1:2).eq.'-t') then
            tst_case_Rayleigh=.not.tst_case_Rayleigh
         endif
         if (argv(1:2).eq.'-U') then
         endif
         if (argv(1:2).eq.'-u') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') nbr_plr
         endif
         if (argv(1:2).eq.'-V') then
            write (6,'(a)') rcs_Id
            goto 1000
         endif
         if (argv(1:2).eq.'-v') then
         endif
         if (argv(1:2).eq.'-W') then
            flg_O3=.not.flg_O3
         endif
         if (argv(1:2).eq.'-w') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O3
         endif
         if (argv(1:2).eq.'-X') then
            flg_NO2=.not.flg_NO2
         endif
         if (argv(1:2).eq.'-x') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_NO2
         endif
         if (argv(1:2).eq.'-Y') then
            flg_OH=.not.flg_OH
         endif
         if (argv(1:2).eq.'-y') then
         endif
         if (argv(1:2).eq.'-Z') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') nbr_azi
         endif
         if (argv(1:2).eq.'-z') then
            cmd_ln_slr_zen_ngl_cos=.not.cmd_ln_slr_zen_ngl_cos
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') slr_zen_ngl_cos_cmd_ln
         endif
      end do

c     Compute any quantities that might depend on command line input.
      flg_Herzberg=flg_O3
      nbr_mom=nbr_str           ! # moments always equals # streams
      if (sv_cmp_plr_ngl) then
         nbr_plr=nbr_str
      endif

      if (top_lvl) then
         stt_top_lvl='Pure absorbing layer to space at model top not yet implemented'
      else
         stt_top_lvl='Pure absorbing layer to space at model top not yet implemented'
      endif
      if (flg_CO2) then
         stt_CO2='CO2 line absorption: monochromatic Voigt line profile'
      else
         stt_CO2='CO2 line absorption: Off'
      endif
      if (flg_H2O) then
         stt_H2O='H2O line absorption: monochromatic Voigt line profile'
      else
         stt_H2O='H2O line absorption: Off'
      endif
      if (flg_H2OH2O) then
         stt_H2OH2O='H2O dimer absorption: Broadband average dimer continuum absorption from Chy97'
      else
         stt_H2OH2O='H2O dimer absorption: Off'
      endif
      if (flg_OH) then
         stt_OH='OH line absorption: monochromatic Voigt line profile'
      else
         stt_OH='OH line absorption: Off'
      endif
      if (flg_O2) then
         stt_O2='O2 line absorption: monochromatic Voigt line profile'
      else
         stt_O2='O2 line absorption: Off'
      endif
      if (flg_Herzberg) then
         stt_Herzberg='O2 Herzberg bands: continuum absorption cross sections from ' // fl_O3
      else
         stt_Herzberg='O2 Herzberg bands: Off'
      endif
      if (flg_O3) then
         stt_O3='O3 Hartley, Huggins, and Chappuis bands: continuum absorption cross sections from ' // fl_O3
      else
         stt_O3='O3 Hartley, Huggins, and Chappuis bands: Off'
      endif
      if (flg_O2O2) then
         stt_O2O2='O2-O2 collision-induced absorption: continuum absorption cross sections from ' // fl_O2O2
      else
         stt_O2O2='O2-O2 collision-induced absorption: Off'
      endif
      if (flg_NO2) then
         stt_NO2='NO2 absorption: continuum absorption cross sections from ' // fl_NO2
      else
         stt_NO2='NO2 absorption: Off'
      endif
      if (flg_ice) then
         stt_ice='Ice water crystal scattering and absorption: Mie theory from ' // fl_ice
      else
         stt_ice='Ice water crystal scattering and absorption: Off'
      endif
      if (flg_liq) then
         stt_liq='Liquid water droplet scattering and absorption: Mie theory from ' // fl_liq
      else
         stt_liq='Liquid water droplet scattering and absorption: Off'
      endif
      if (flg_Planck) then
         stt_Planck='Thermal emission of atmosphere: Included for lambda > 2 microns'
      else
         stt_Planck='Thermal emission of atmosphere: Off'
      endif
      if (flg_Rayleigh) then
         stt_Rayleigh='Rayleigh scattering: method of Lenoble (1993)'
      else
         stt_Rayleigh='Rayleigh scattering: Off'
      endif
      if (flg_aer) then
         stt_aer='Aerosol scattering and absorption: Mie theory from ' // fl_aer
      else
         stt_aer='Aerosol scattering and absorption: Off'
      endif
      if (flg_bga) then
         stt_bga='Background aerosol scattering and absorption: Mie theory from ' // fl_bga
      else
         stt_bga='Background aerosol scattering and absorption: Off'
      endif
      if (flt_nst) then
         stt_flt_nst='Instrument filter spectral response function: instrument from ' // fl_nst
      else
         stt_flt_nst='Instrument filter spectral response function: Off'
      endif
      stt_slr='TOA solar spectrum from ' // fl_slr

c     Set netCDF error handling
      call ncpopt(NCVERBOS+NCFATAL)

c     Ingest fl_clm
      rcd=nf_open(fl_clm,NCNOWRIT,nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_clm)

c     Get global attributes
c      call ncainq(nc_id,NF_GLOBAL,'',int_foo,int_foo,rcd)
      call ncagtc(nc_id,NF_GLOBAL,'prf_sng',prf_sng,len(prf_sng),rcd)

c     Get dimension IDs
      lev_dim_id=ncdid(nc_id,'lev',rcd)
      levp_dim_id=ncdid(nc_id,'levp',rcd)

c     Get dimension sizes
      call ncdinq(nc_id,lev_dim_id,sng,nbr_lev,rcd)
      if (nbr_lev.gt.max_nbr_lev) stop 'nbr_lev.gt.max_nbr_lev'
      call ncdinq(nc_id,levp_dim_id,sng,nbr_levp,rcd)
      if (nbr_levp.gt.max_nbr_levp) stop 'nbr_levp.gt.max_nbr_levp'

c     Get variable IDs
      RH_liq_id=ncvid(nc_id,'RH_liq',rcd)
      alb_NIR_str_zen_ngl_id=ncvid(nc_id,'alb_NIR_str_zen_ngl',rcd)
      alb_NIR_weak_zen_ngl_id=ncvid(nc_id,'alb_NIR_weak_zen_ngl',rcd)
      alb_vis_str_zen_ngl_id=ncvid(nc_id,'alb_vis_str_zen_ngl',rcd)
      alb_vis_weak_zen_ngl_id=ncvid(nc_id,'alb_vis_weak_zen_ngl',rcd)
      slr_zen_ngl_cos_id=ncvid(nc_id,'slr_zen_ngl_cos',rcd)
      prs_dlt_id=ncvid(nc_id,'prs_dlt',rcd)
      frc_ice_id=ncvid(nc_id,'frc_ice',rcd)
      frc_ice_ttl_id=ncvid(nc_id,'frc_ice_ttl',rcd)
      grv_id=ncvid(nc_id,'grv',rcd)
      lat_dgr_id=ncvid(nc_id,'lat_dgr',rcd)
      lcl_time_hr_id=ncvid(nc_id,'lcl_time_hr',rcd)
      lcl_yr_day_id=ncvid(nc_id,'lcl_yr_day',rcd)
      lev_id=ncvid(nc_id,'lev',rcd)
      levp_id=ncvid(nc_id,'levp',rcd)
      mmw_mst_air_id=ncvid(nc_id,'mmw_mst_air',rcd)
      mpc_CWP_id=ncvid(nc_id,'mpc_CWP',rcd)
      mpl_CO2_id=ncvid(nc_id,'mpl_CO2',rcd)
      mpl_CWP_id=ncvid(nc_id,'mpl_CWP',rcd)
      mpl_H2O_id=ncvid(nc_id,'mpl_H2O',rcd)
      mpl_IWP_id=ncvid(nc_id,'mpl_IWP',rcd)
      mpl_LWP_id=ncvid(nc_id,'mpl_LWP',rcd)
      mpl_O2_id=ncvid(nc_id,'mpl_O2',rcd)
      mpl_OH_id=ncvid(nc_id,'mpl_OH',rcd)
      mpl_aer_id=ncvid(nc_id,'mpl_aer',rcd)
      mpl_bga_id=ncvid(nc_id,'mpl_bga',rcd)
      mpl_mst_air_id=ncvid(nc_id,'mpl_mst_air',rcd)
      npl_NO2_id=ncvid(nc_id,'npl_NO2',rcd)
      npl_O2_id=ncvid(nc_id,'npl_O2',rcd)
      npl_O3_id=ncvid(nc_id,'npl_O3',rcd)
      npl_O2O2_id=ncvid(nc_id,'npl_O2O2',rcd)
      odxc_obs_aer_id=ncvid(nc_id,'odxc_obs_aer',rcd)
      odxc_obs_bga_id=ncvid(nc_id,'odxc_obs_bga',rcd)
      odxl_obs_aer_id=ncvid(nc_id,'odxl_obs_aer',rcd)
      odxl_obs_bga_id=ncvid(nc_id,'odxl_obs_bga',rcd)
      p_id=ncvid(nc_id,'p',rcd)
      p_int_id=ncvid(nc_id,'p_int',rcd)
      q_CO2_id=ncvid(nc_id,'q_CO2',rcd)
      q_H2O_id=ncvid(nc_id,'q_H2O',rcd)
      q_O2_id=ncvid(nc_id,'q_O2',rcd)
      q_OH_id=ncvid(nc_id,'q_OH',rcd)
      spec_heat_mst_air_id=ncvid(nc_id,'spec_heat_mst_air',rcd)
      t_id=ncvid(nc_id,'t',rcd)
      t_int_id=ncvid(nc_id,'t_int',rcd)
      t_skn_id=ncvid(nc_id,'t_skn',rcd)
      xnt_fac_id=ncvid(nc_id,'xnt_fac',rcd)
      wvl_obs_aer_id=ncvid(nc_id,'wvl_obs_aer',rcd)
      wvl_obs_bga_id=ncvid(nc_id,'wvl_obs_bga',rcd)
      z_cld_btm_id=ncvid(nc_id,'z_cld_btm',rcd)
      z_cld_thick_id=ncvid(nc_id,'z_cld_thick',rcd)
      z_id=ncvid(nc_id,'z',rcd)
      z_int_id=ncvid(nc_id,'z_int',rcd)

c     Get data
      call ncvgt(nc_id,RH_liq_id,srt,nbr_lev,RH_liq,rcd)
      call ncvgt(nc_id,prs_dlt_id,srt,nbr_lev,prs_dlt,rcd)
      call ncvgt(nc_id,frc_ice_id,srt,nbr_lev,frc_ice,rcd)
      call ncvgt(nc_id,grv_id,srt,nbr_lev,grv,rcd)
      call ncvgt(nc_id,lev_id,srt,nbr_lev,lev,rcd)
      call ncvgt(nc_id,levp_id,srt,nbr_levp,levp,rcd)
      call ncvgt(nc_id,mmw_mst_air_id,srt,nbr_lev,mmw_mst_air,rcd)
      call ncvgt(nc_id,mpl_CO2_id,srt,nbr_lev,mpl_CO2,rcd)
      call ncvgt(nc_id,mpl_CWP_id,srt,nbr_lev,mpl_CWP,rcd)
      call ncvgt(nc_id,mpl_H2O_id,srt,nbr_lev,mpl_H2O,rcd)
      call ncvgt(nc_id,mpl_IWP_id,srt,nbr_lev,mpl_IWP,rcd)
      call ncvgt(nc_id,mpl_LWP_id,srt,nbr_lev,mpl_LWP,rcd)
      call ncvgt(nc_id,mpl_O2_id,srt,nbr_lev,mpl_O2,rcd)
      call ncvgt(nc_id,mpl_OH_id,srt,nbr_lev,mpl_OH,rcd)
      call ncvgt(nc_id,mpl_aer_id,srt,nbr_lev,mpl_aer,rcd)
      call ncvgt(nc_id,mpl_bga_id,srt,nbr_lev,mpl_bga,rcd)
      call ncvgt(nc_id,mpl_mst_air_id,srt,nbr_lev,mpl_mst_air,rcd)
      call ncvgt(nc_id,npl_NO2_id,srt,nbr_lev,npl_NO2,rcd)
      call ncvgt(nc_id,npl_O2_id,srt,nbr_lev,npl_O2,rcd)
      call ncvgt(nc_id,npl_O3_id,srt,nbr_lev,npl_O3,rcd)
      call ncvgt(nc_id,npl_O2O2_id,srt,nbr_lev,npl_O2O2,rcd)
      call ncvgt(nc_id,odxl_obs_aer_id,srt,nbr_lev,odxl_obs_aer,rcd)
      call ncvgt(nc_id,odxl_obs_bga_id,srt,nbr_lev,odxl_obs_bga,rcd)
      call ncvgt(nc_id,p_id,srt,nbr_lev,p,rcd)
      call ncvgt(nc_id,p_int_id,srt,nbr_levp,p_int,rcd)
      call ncvgt(nc_id,q_CO2_id,srt,nbr_lev,q_CO2,rcd)
      call ncvgt(nc_id,q_H2O_id,srt,nbr_lev,q_H2O,rcd)
      call ncvgt(nc_id,q_O2_id,srt,nbr_lev,q_O2,rcd)
      call ncvgt(nc_id,q_OH_id,srt,nbr_lev,q_OH,rcd)
      call ncvgt(nc_id,spec_heat_mst_air_id,srt,nbr_lev,spec_heat_mst_air,rcd)
      call ncvgt(nc_id,t_id,srt,nbr_lev,t,rcd)
      call ncvgt(nc_id,t_int_id,srt,nbr_levp,t_int,rcd)
      call ncvgt(nc_id,z_id,srt,nbr_lev,z,rcd)
      call ncvgt(nc_id,z_int_id,srt,nbr_levp,z_int,rcd)
      call ncvgt1(nc_id,alb_NIR_str_zen_ngl_id,0,alb_NIR_str_zen_ngl,rcd)
      call ncvgt1(nc_id,alb_NIR_weak_zen_ngl_id,0,alb_NIR_weak_zen_ngl,rcd)
      call ncvgt1(nc_id,alb_vis_str_zen_ngl_id,0,alb_vis_str_zen_ngl,rcd)
      call ncvgt1(nc_id,alb_vis_weak_zen_ngl_id,0,alb_vis_weak_zen_ngl,rcd)
      call ncvgt1(nc_id,slr_zen_ngl_cos_id,0,slr_zen_ngl_cos,rcd)
      call ncvgt1(nc_id,frc_ice_ttl_id,0,frc_ice_ttl,rcd)
      call ncvgt1(nc_id,lat_dgr_id,0,lat_dgr,rcd)
      call ncvgt1(nc_id,lcl_time_hr_id,0,lcl_time_hr,rcd)
      call ncvgt1(nc_id,lcl_yr_day_id,0,lcl_yr_day,rcd)
      call ncvgt1(nc_id,mpc_CWP_id,0,mpc_CWP,rcd)
      call ncvgt1(nc_id,odxc_obs_aer_id,0,odxc_obs_aer,rcd)
      call ncvgt1(nc_id,odxc_obs_bga_id,0,odxc_obs_bga,rcd)
      call ncvgt1(nc_id,t_skn_id,0,t_skn,rcd)
      call ncvgt1(nc_id,xnt_fac_id,0,xnt_fac,rcd)
      call ncvgt1(nc_id,wvl_obs_aer_id,0,wvl_obs_aer,rcd)
      call ncvgt1(nc_id,wvl_obs_bga_id,0,wvl_obs_bga,rcd)
      call ncvgt1(nc_id,z_cld_btm_id,0,z_cld_btm,rcd)
      call ncvgt1(nc_id,z_cld_thick_id,0,z_cld_thick,rcd)

      rcd=nf_close(nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_clm)
      write (6,'(a20,1x,a40)') 'Ingested netCDF file',fl_clm

c     Ingest fl_htrn
      rcd=nf_open(fl_htrn,NCNOWRIT,nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_htrn)
      
c     Get attributes
      rcd=nf_get_att_double(nc_id,NF_GLOBAL,'wvn_min_rgn',wvn_min_rgn)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR wvn_min_rgn')
      rcd=nf_get_att_double(nc_id,NF_GLOBAL,'wvn_max_rgn',wvn_max_rgn)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR wvn_max_rgn')

c     Get dimension IDs
      rcd=nf_inq_dimid(nc_id,'ln_ctr',ln_ctr_dim_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_ctr_dim_id')

c     Get dimension sizes
      rcd=nf_inq_dimlen(nc_id,ln_ctr_dim_id,ln_nbr)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_nbr')
      if (ln_nbr.gt.max_ln_nbr) stop 'ln_nbr.gt.max_ln_nbr'

c     Get variable IDs
      rcd=nf_inq_varid(nc_id,'HWHM_air',HWHM_air_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR HWHM_air_id')
      rcd=nf_inq_varid(nc_id,'HWHM_tpt_dpn_xpn',HWHM_tpt_dpn_xpn_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR HWHM_tpt_dpn_xpn_id')
      rcd=nf_inq_varid(nc_id,'ln_ctr',ln_ctr_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_ctr_id')
      rcd=nf_inq_varid(nc_id,'ln_iso',ln_iso_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_iso_id')
      rcd=nf_inq_varid(nc_id,'ln_mlc',ln_mlc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_mlc_id')
      rcd=nf_inq_varid(nc_id,'ln_nrg_lwr',ln_nrg_lwr_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_nrg_lwr_id')
      rcd=nf_inq_varid(nc_id,'ln_str',ln_str_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_str_id')

c     Get data
      srt(1)=1
      
      rcd=nf_get_vara_double(nc_id,HWHM_air_id,srt,ln_nbr,HWHM_air)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR HWHM_air')
      rcd=nf_get_vara_double(nc_id,HWHM_tpt_dpn_xpn_id,srt,ln_nbr,HWHM_tpt_dpn_xpn)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR HWHM_tpt_dpn_xpn')
      rcd=nf_get_vara_double(nc_id,ln_ctr_id,srt,ln_nbr,ln_ctr)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_ctr')
      rcd=nf_get_vara_int2(nc_id,ln_iso_id,srt,ln_nbr,ln_iso)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_iso')
      rcd=nf_get_vara_int2(nc_id,ln_mlc_id,srt,ln_nbr,ln_mlc)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_mlc')
      rcd=nf_get_vara_double(nc_id,ln_nrg_lwr_id,srt,ln_nbr,ln_nrg_lwr)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_nrg_lwr')
      rcd=nf_get_vara_double(nc_id,ln_str_id,srt,ln_nbr,ln_str)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR ln_str')

      rcd=nf_close(nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_htrn)
      write (6,'(a20,1x,a40)') 'Ingested netCDF file',fl_htrn
      
c     Ingest fl_nst
      rcd=nf_open(fl_nst,NCNOWRIT,nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_nst)
      
c     Get dimension IDs
      bnd_dim_id=ncdid(nc_id,'wvl',rcd)
      
c     Get dimension sizes
      call ncdinq(nc_id,bnd_dim_id,sng,nbr_nst_bnd,rcd)
      if (nbr_nst_bnd.gt.max_nbr_nst_bnd) stop 'nbr_nst_bnd.gt.max_nbr_nst_bnd'
      
c     Get variable IDs
      nst_SRF_id=ncvid(nc_id,'nst_SRF',rcd)
      wvl_max_nst_id=ncvid(nc_id,'wvl_max',rcd)
      wvl_min_nst_id=ncvid(nc_id,'wvl_min',rcd)
      wvl_ctr_nst_id=ncvid(nc_id,'wvl',rcd)

c     Get data
      call ncvgt(nc_id,nst_SRF_id,srt,nbr_nst_bnd,nst_SRF,rcd)
      call ncvgt(nc_id,wvl_max_nst_id,srt,nbr_nst_bnd,wvl_max_nst,rcd)
      call ncvgt(nc_id,wvl_min_nst_id,srt,nbr_nst_bnd,wvl_min_nst,rcd)
      call ncvgt(nc_id,wvl_ctr_nst_id,srt,nbr_nst_bnd,wvl_ctr_nst,rcd)

      rcd=nf_close(nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_nst)
      write (6,'(a20,1x,a40)') 'Ingested netCDF file',fl_nst

c     All necessary input data has been read
c     Initialize input data that may overridden by command line switches
      if (cmd_ln_odxc_obs_aer) then
         aer_lvl_nbr=0
         do lev_idx=1,nbr_lev
            if (odxl_obs_aer(lev_idx).gt.0.) aer_lvl_nbr=aer_lvl_nbr+1
         enddo                  ! end loop over lev
         if (aer_lvl_nbr.eq.0) stop 'aer_lvl_nbr.eq.0'
         do lev_idx=1,nbr_lev
            if (odxl_obs_aer(lev_idx).gt.0.) odxl_obs_aer(lev_idx)=odxc_obs_aer_cmd_ln/aer_lvl_nbr
         enddo                  ! end loop over lev
c     Recompute column path
         odxc_obs_aer=0.
         do lev_idx=1,nbr_lev
               odxc_obs_aer=odxc_obs_aer+odxl_obs_aer(lev_idx)
         enddo                  ! end loop over lev
      endif                     ! end if overriding CLM profile odxc_obs_aer

      if (cmd_ln_mpc_CWP) then
         cld_lvl_nbr=0
         do lev_idx=1,nbr_lev
            if (mpl_CWP(lev_idx).gt.0.) cld_lvl_nbr=cld_lvl_nbr+1
         enddo                  ! end loop over lev
         if (cld_lvl_nbr.eq.0) stop 'cld_lvl_nbr.eq.0'
         do lev_idx=1,nbr_lev
            if (mpl_CWP(lev_idx).gt.0.) mpl_CWP(lev_idx)=mpc_CWP_cmd_ln/cld_lvl_nbr
         enddo                  ! end loop over lev
      endif                     ! end if overriding CLM profile mpl_CWP
      if (force_ice_phz) then
         do lev_idx=1,nbr_lev
            frc_ice(lev_idx)=1.
         enddo                  ! end loop over lev
      endif                     ! end if forcing ice phase mpl_CWP
      if (force_liq_phz) then
         do lev_idx=1,nbr_lev
            frc_ice(lev_idx)=0.
         enddo                  ! end loop over lev
      endif                     ! end if forcing liquid phase mpl_CWP
      if (cmd_ln_mpc_CWP.or.force_ice_phz.or.force_liq_phz) then
         mpc_CWP=0.
         do lev_idx=1,nbr_lev
               mpl_IWP(lev_idx)=frc_ice(lev_idx)*mpl_CWP(lev_idx)
               mpl_LWP(lev_idx)=max(0.,mpl_CWP(lev_idx)-mpl_IWP(lev_idx))
               mpc_CWP=mpc_CWP+mpl_CWP(lev_idx)
               mpc_IWP=mpc_IWP+mpl_IWP(lev_idx)
         enddo                  ! end loop over lev
         if (mpc_CWP.ne.0.) then
            frc_ice_ttl=mpc_IWP/mpc_CWP
         else                   ! end if column is cloudy
            frc_ice_ttl=0.
            z_cld_btm=0.
            z_cld_thick=0.
         endif                  ! end if column is clear
      endif                     ! end if recomputing mpl_CWP, mpl_IWP, mpl_LWP
      if (cmd_ln_slr_zen_ngl_cos) then
         slr_zen_ngl_cos=slr_zen_ngl_cos_cmd_ln
         xnt_fac=1.
      endif                     ! end if overriding CLM profile zen ang
      if (cmd_ln_alb) then
         alb_NIR_str_zen_ngl=alb_cmd_ln
         alb_NIR_weak_zen_ngl=alb_cmd_ln
         alb_vis_str_zen_ngl=alb_cmd_ln
         alb_vis_weak_zen_ngl=alb_cmd_ln
         alb_sfc=alb_cmd_ln
      else
         alb_sfc=.5*(alb_vis_str_zen_ngl+alb_NIR_str_zen_ngl)
      endif                     ! end if overriding CLM profile albedo
      if (cmd_ln_slr_cst) then
         slr_cst=slr_cst_cmd_ln
      endif                     ! end if overriding solar constant

c     Many DISORT arguments only need to be set once. 
c     Set them here, outside band loop.
      nstr=nbr_str
      nlyr=nbr_lev

c     Temperatures are specified at interfaces.
      do lev_idx=0,nbr_lev
         temper(lev_idx)=t_int(lev_idx+1)
      end do

c     NB: Setting usrtau to .true. has grave consequences on 
c     lengths of arrays returned by DISORT. In particular,
c     it is safer to keep usrtau = .false. because then all 
c     radiant quantities are returned at boundaries to 
c     computational layers.
      usrtau=.not.sv_cmp_opt_dep
      if (sv_cmp_opt_dep) then 
         nbr_tau=nbr_levp
      else
         nbr_tau=nbr_levp
      endif                     ! end if no user defined levels
      ntau=nbr_tau
      do tau_idx=1,nbr_tau
         tau_prs(tau_idx)=p_int(tau_idx)
      end do                    ! end loop over tau

c     Set angles at which to report intensities.
c     Array must be in increasing order of cosine polar angle.
      usrang=.not.sv_cmp_plr_ngl
      if (usrang) then
         plr(1)=pi
         do plr_idx=2,nbr_plr-1
            plr(plr_idx)=pi-(plr_idx-1)*pi/(nbr_plr-1)
         enddo                  ! end loop over plr
         plr(nbr_plr)=0.

         do plr_idx=1,nbr_plr
            plr_dgr(plr_idx)=180.*plr(plr_idx)/pi
         enddo                  ! end loop over plr

         plr_cos(1)=-1.0        ! 180 degrees, i.e., towards nadir
         do plr_idx=2,nbr_plr-1
            plr_cos(plr_idx)=cos(plr(plr_idx))
         enddo                  ! end loop over plr
         plr_cos(nbr_plr)=1.0   ! 0 degrees, i.e., towards zenith

         numu=nbr_plr
         do plr_idx=1,nbr_plr
            umu(plr_idx)=plr_cos(plr_idx)
         enddo                  ! end loop over plr
      endif                     ! endif if looking at user angles

c     Set azimuthal angles at which to report intensities.
c     nphi=0 is valid only when onlyfl=.true.
      nphi=nbr_azi
      do azi_idx=1,nbr_azi
         azi(azi_idx)=(azi_idx-1)*2.*pi/nbr_azi
         azi_dgr(azi_idx)=(azi_idx-1)*360./nbr_azi
         phi(azi_idx)=azi_dgr(azi_idx)
      enddo                     ! end loop over azimuthal angles

c     Set boundary conditions
c     0 includes case of beam illumination from top.
      ibcnd=0

c     Set cosine of the solar zenith angle. 
c     Beware of using .5 when nbr_streams/2 is odd---it's a 
c     quadrature point and can crash DISORT.
      umu0=slr_zen_ngl_cos

c     Set azimuth angle of incident sunlight.
      phi0=0.

c     Set incident isotropic illumination at top boundary.
      fisot=0.

c     Specify bottom boundary temperature btemp. 
c     Btemp is skin temperature in CCM parlance.
c     Bottom boundary emissivity is derived from albedo above. 
c     Top boundary temperature and emissivity must also be specified.
c     Of course these are only used when plank is .true..
      btemp=t_skn
      ttemp=t_int(1)
      temis=0.                  ! i have no idea what this should be.

c     Set control flags:
c     Set deltam=.true. unless looking at radiances within 10 degrees of forward peak.
      deltam=.true.

c     Sanity check
      if (slr_zen_ngl_cos.lt.0..and..not.flg_Planck) then 
         write (6,*) 'Sun beneath horizon for pure solar calculation'
         call abort
      endif

c     Set onlyfl=.false. when looking at user angle radiances.
c     When onlyfl=.true., azimuthally averaged radiances at 
c     computational angles--NOT user angles--will be reported 
c     in u0u whenever maxumu .ge. nstr.
      onlyfl=.false.

c     accur is maximum relative error in last three terms 
c     in azimuthal series, and it determines convergence.
c     Value of accur does not seem to contribute to any problems
c     with single precision.
c     Set 0. < accur < .01
      accur=0.

c     Set printing flags for DISORT:
c     prnt(1)=.true.: print input variables (except pmom)
      prnt(1)=.false.

c     prnt(2)=.true.: print fluxes
      prnt(2)=.false.

c     prnt(3)=.true.: print azimuthally-averaged intensities at user
c     levels and computational polar angles.
      prnt(3)=.false.

c     prnt(4)=.true.: print azimuthally-averaged intensities at user
c     levels and user polar angles.
      prnt(4)=.false.

c     prnt(5)=.true.: print intensities at user levels and angles.
      prnt(5)=.false.

c     prnt(6)=.true.: print planar transmissivity and planar albedo 
c     as a function of solar zenith angle (for ibcnd=1)
      prnt(6)=.false.

c     prnt(7)=.true.: print phase function moments 'pmom' for each
c     layer (only if prnt(1)=.true., and only for layers with scattering)
      prnt(7)=.false.

c     Band-independent DISORT() initialization is now complete.
c     Remaining DISORT arguments need to be set inside main band loop.

c     Given spectral range of valid lines, divide by resolution to obtain number of monochromatic intervals.     
      bnd_nbr=(wvn_max_rgn-wvn_min_rgn)/rsn_wvn

c     Initialize level-independent arrays that depend on bnd_nbr.
      do bnd_idx=1,bnd_nbr
         wvn_min(bnd_idx)=wvn_min_rgn+(bnd_idx-1)*rsn_wvn
         wvn_max(bnd_idx)=wvn_min_rgn+bnd_idx*rsn_wvn
         wvn_ctr(bnd_idx)=
     $        .5*(wvn_min(bnd_idx)+wvn_max(bnd_idx))
         wvl_ctr(bnd_idx)=1./(100.*wvn_ctr(bnd_idx))
         wvl_min(bnd_idx)=1./(100.*wvn_max(bnd_idx))
         wvl_max(bnd_idx)=1./(100.*wvn_min(bnd_idx))
         wvl_sz(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
      enddo                     ! end loop over bnd

c     Get TOA solar spectrum
      call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc)

c     Let user know where the wavelength chips have fallen.
      write (str_sng,'(a45,i2,a8)') 'Discrete ordinate computation performed with ',nbr_str,' streams'
      if (sv_cmp_plr_ngl) then
         write (plr_sng,'(a24,i2,a28)') 'Intensities reported at ',nbr_plr,' computational polar angles.'
      else
         write (plr_sng,'(a24,i2,a33)') 'Intensities reported at ',nbr_plr,' evenly spaced user polar angles.'
      endif
      write (azi_sng,'(a24,i2,a48,i4,a3,f7.5,a8)') 
     $     'Intensities reported at ',nbr_azi,
     $     ' evenly spaced user azimuthal angles for band = ',
     $     bnd_dbg,' = ',wvl_ctr(bnd_dbg)*1.e6,' microns'
      write (aer_sng,'(a46,i4,a3,f7.5,a8)') 
     $     'Layer aerosol optical depths saved for band = ',
     $     bnd_obs_aer,' = ',wvl_ctr(bnd_obs_aer)*1.e6,' microns'
      write (bga_sng,'(a57,i4,a3,f7.5,a8)') 
     $     'Layer background aerosol optical depths saved for band = ',
     $     bnd_obs_bga,' = ',wvl_ctr(bnd_obs_bga)*1.e6,' microns'
      if (sv_cmp_opt_dep) then
         write (opt_dep_sng,'(a33,i3,a27)') 'Radiative quantities reported at ',nbr_tau,' interface pressure levels.'
      else
         write (opt_dep_sng,'(a33,i3,a27)') 'Radiative quantities reported at ',nbr_tau,' user defined optical depths.'
      endif

      write (6,'(25(a,/))')       
     $     str_sng,
     $     prf_sng,
     $     plr_sng,
     $     opt_dep_sng,
     $     azi_sng,
     $     aer_sng,
     $     bga_sng,
     $     stt_CO2,
     $     stt_H2O,
     $     stt_H2OH2O,
     $     stt_Herzberg,
     $     stt_ice,
     $     stt_liq,
     $     stt_OH,
     $     stt_O2,
     $     stt_O3,
     $     stt_O2O2,
     $     stt_NO2,
     $     stt_Planck,
     $     stt_Rayleigh,
     $     stt_aer,
     $     stt_bga,
     $     stt_flt_nst,
     $     stt_slr,
     $     stt_top_lvl
      if (dbg_lvl.gt.0) then
         write (6,'(
     $        a14,f9.6,a5,
     $        a14,f9.6,a5,
     $        a14,f9.6,a5,
     $        a14,f9.6,a5,
     $        a10,i6,
     $        a10,i6,
     $        a10,i6
     $        )')       
     $        'wvn_min_rgn = ',wvn_min_rgn,' cm-1',
     $        'wvn_max_rgn = ',wvn_max_rgn,' cm-1',
     $        'wvl_min_rgn = ',wvl_min_rgn*1.e6,' um',
     $        'wvl_max_rgn = ',wvl_max_rgn*1.e6,' um',
     $        'ntv_nbr = ',ntv_nbr,
     $        'ln_nbr  = ',ln_nbr,
     $        'bnd_nbr = ',bnd_nbr
      endif                     ! end if dbg

      slr_cst_xnt_fac=slr_cst*xnt_fac
      do bnd_idx=1,nbr_bnd
         odxc_spc_ln(bnd_idx)=0.
         odxc_spc_ttl(bnd_idx)=0.
      enddo                     ! end loop over bnd

c     Initialize band-independent arrays that depend on level.
      do lev_idx=1,nbr_lev
         t0_rcp_t(lev_idx)=tpt_HITRAN/t(lev_idx)
         p_rcp_p0(lev_idx)=p(lev_idx)/p_HITRAN
      enddo                     ! end loop over lev

c     Zero broad-band arrays which will accumulate spectral fluxes.
      do lev_idx=1,nbr_lev
         ntn_bb_mean(lev_idx)=0.
      enddo                     ! end loop over lev
      do lev_idx=1,nbr_levp
         do plr_idx=1,nbr_plr
            ntn_bb_aa(plr_idx,lev_idx)=0.
         enddo                  ! end loop over plr
         flx_bb_dwn_drc(lev_idx)=0.
         flx_bb_dwn_dff(lev_idx)=0.
         flx_bb_up(lev_idx)=0.
      enddo                     ! end loop over lev

c     End section 1: Initialization
c     Begin section 2: Main computation loop over all bands
c     begin main loop over bands
c     assignment of counting indices into input data arrays
c     a level loop for zeroing band quantities that depend on level
c     computation and binning of process-specific optical properties
c     summation of individual optical properties into bulk layer props
c     assignment of bulk properties into actual DISORT() input arrays
c     call DISORT()
c     assignment of DISORT() output to final output storage arrays
c     end main loop over bands

      do bnd_idx=1,nbr_bnd

         write (6,'(a1,$)') '.'

c     O3 data is reverse indexed with an offset relative to 
c     combined band data, so index juggling is necessary.
         bnd_idx_O3=nbr_pure_O3_bnd-(bnd_idx-nbr_H2O_bnd-1)

c     H2O data is indexed identically to the combined band data.
         bnd_idx_H2O=bnd_idx

c     Before anything is actually computed for this band, 
c     zero out all optical depth arrays. This is important
c     both as a failsafe programming technique and for the case
c     when arguments are given to ignore specified radiative
c     processes, e.g., the pure Rayleigh scattering case.
c     Also, at least one optical property, odal_O2, is
c     accumulated between two processes (Herzberg continuum absorption
c     and line absorption near .69, .77, and 1.28 microns in the 
c     case of O2).
         do lev_idx=1,nbr_lev
            odsl_Ray(lev_idx)=0.
            odsl_ice(lev_idx)=0.
            odsl_liq(lev_idx)=0.
            odsl_aer(lev_idx)=0.
            odsl_bga(lev_idx)=0.
            odal_ice(lev_idx)=0.
            odal_liq(lev_idx)=0.
            odal_aer(lev_idx)=0.
            odal_bga(lev_idx)=0.
            odal_OH(lev_idx)=0.
            odal_O2(lev_idx)=0.
            odal_O3(lev_idx)=0.
            odal_O2O2(lev_idx)=0.
            odal_NO2(lev_idx)=0.
            odal_CO2(lev_idx)=0.
            odal_H2OH2O(lev_idx)=0.
            odal_H2O(lev_idx)=0.
         enddo                  ! end loop over lev

c     Convert input data to SI units where necessary. 
c     HITRAN line data is an exception to the rule that I store everything in SI units.
c     Line data is stored exactly as it appears in the HITRAN database, which is basically CGS. 
c     The reason why ln_str should be scaled to SI, even though there is a 
c     possible loss of precision involved for weak lines in single precision, 
c     is that ln_str is used in all sorts of SI algebra involving fundamental constants 
c     and it's easiest to do it all in SI.
         do ln_idx=1,ln_nbr
            ln_str(ln_idx)=ln_str(ln_idx)*1.e-4 ! cm2/mlc -> m2/mlc
            ln_str(ln_idx)=ln_str(ln_idx)*Avagadro ! m2/mlc -> m2/mol
         enddo                  ! end loop over lines
         
c     LBL computation includes contributions from line wings within 25 cm-1.
c     See GoY89 p. 109 for discussion of line wings.
         wng_sz_wvn=25.

c     Fill in ln_idx_srt and ln_idx_end for this line
         call ln_nbr_get(ln_ctr,ln_nbr,wng_sz_wvn,ln_idx_srt,ln_idx_end)

c     Compute the profiles of all these lines in the current spectral interval
         ln_nbr_ngh=ln_idx_end-ln_idx_srt+1

c     Initialize the arrays of molecular properties
         call rtl_fnc_tpt_xpn_get(xpn_mlc)
         call mmw_iso_get(mmw_iso)
         call iso_idx_map_get(map)

         do ln_idx=ln_idx_srt,ln_idx_end
            iso_idx=map(ln_iso(ln_idx),ln_mlc(ln_idx))
            do lev_idx=1,nbr_lev
c     Compute Lorentz half-width at half-max. See, e.g., GoY89 p. 99 (3.51), Lio92 p. 28 (2.2.6)
               HWHM_Lorentz=HWHM_air(ln_idx)*p_rcp_p0(lev_idx)*t0_rcp_t(lev_idx)**HWHM_tpt_dpn_xpn(ln_idx)
c     Compute Doppler half-width at e-1 of max. See, e.g., GoY89 p. 111 (3.82), Lio92 p. 30 (2.2.9)
               HWeM_Doppler=(ln_ctr(ln_idx)/speed_of_light)*sqrt(2.*Boltzmann*t(lev_idx)/mmw_iso(iso_idx))  
            enddo               ! end loop over lev
c     Compute Voigt profile. See, e.g., GoY89 p. 112 (3.83), Lio92 p. 31 (2.2.11)
c     The Voigt profile is symmetric in frequency space, but always give argument as
c     nu-nu0, just in case.
            ln_dst_ctr=wvn_ctr(bnd_idx)-ln_ctr(ln_idx)
            Voigt_arg_x=ln_dst_ctr/HWeM_Doppler
            Voigt_arg_y=HWHM_Lorentz/HWeM_Doppler
            Voigt_arg_z=cmplx(Voigt_arg_x,Voigt_arg_y)
            Voigt_fac=sqrt_rcp_pi/HWeM_Doppler
            ln_prf=Voigt_fac*real(Voigt_Hum82(Voigt_arg_z))
c     Compute line strength. See, e.g., Lio92 p. 33 (2.2.22a)
            do lev_idx=1,nbr_lev
               odal_H2O(lev_idx)=odal_H2O(lev_idx)+ln_str(ln_idx)*ln_prf
            enddo               ! end loop over lev
         enddo                  ! end loop over ln

c     All the line contributions from the current gas to the current bin have been computed.

c     Add on the continuum contributions to the current bin.

c     Ready to invoke DISORT

c     Need to move from English tags to hashed indices for atmospheric constituents

c     Continuum contribution #1: H2O
            
         if ((bnd_idx_O3.gt.nbr_pure_O3_bnd).and.
     $        (bnd_idx.gt.nbr_non_O3_bnd)) then
c     We're in overlap region where there is both H2O and continuum data. 
c     We'll need to interpolate O3,O2 continuum absorption data from
c     O3 bins onto the H2O bins. This is messy.
c     For now we just search for the O3 band which contains the H2O
c     band center, and assign that O3 cross-section to that H2O band.
 120        if ((wvl_ctr(bnd_idx).lt.wvl_min_O3(bnd_idx_tpt_O3)).or.
     $           (wvl_ctr(bnd_idx).gt.wvl_max_O3(bnd_idx_tmp_O3))) then
               bnd_idx_tmp_O3=bnd_idx_tmp_O3-1
               goto 120
            endif

            do lev_idx=1,nbr_lev
               odal_O3(lev_idx)=
     $              abs_xsx_O3(bnd_idx_tmp_O3)*npl_O3(lev_idx)
               odal_O2(lev_idx)=odal_O2(lev_idx)+
     $              abs_xsx_O2(bnd_idx_tmp_O3)*npl_O2(lev_idx)
            enddo               ! end loop over lev

         endif                  ! endif in overlap region

         if (bnd_idx_O3.eq.nbr_pure_O3_bnd) then
c     Special care must be taken with this band as it is the band we truncated. 
c     O3,O2 cross-sections here should be interpolated.
            do lev_idx=1,nbr_lev
               odal_O3(lev_idx)=
     $              abs_xsx_O3(bnd_idx_O3)*npl_O3(lev_idx)
               odal_O2(lev_idx)=odal_O2(lev_idx)+
     $              abs_xsx_O2(bnd_idx_O3)*npl_O2(lev_idx)
            enddo               ! end loop over lev
         endif                  ! endif band is the splice band

         if ((bnd_idx_O3.lt.nbr_pure_O3_bnd).and.
     $        (bnd_idx.le.nbr_bnd)) then
c     We're in pure continuum region. 
c     Absorption will be due solely to O3, O2.
            do lev_idx=1,nbr_lev
               odal_O3(lev_idx)=
     $              abs_xsx_O3(bnd_idx_O3)*npl_O3(lev_idx)
               odal_O2(lev_idx)=odal_O2(lev_idx)+
     $              abs_xsx_O2(bnd_idx_O3)*npl_O2(lev_idx)
            enddo               ! end loop over lev
         endif                  ! end if outside H2O data

c     O2-O2 absorption
c     NB: Possible difficulty with O2-O2 cross-sections is their high resolution (up to .1 nm).
c     Narrow band model resolution, 10 cm-1, is 1 nm at 1.06 um, .3 nm at .55 um, and 
c     roughly 5 nm beneath (bluer than) .5 um. 
c     A professional would average the cross-section data (weighted by insolation, no less!)
c     up to the grid size of the narrow band model.
 190     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_O2O2(1)) then
            bnd_idx_O2O2=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_O2O2(nbr_O2O2_bnd)) then
            bnd_idx_O2O2=nbr_O2O2_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_O2O2(bnd_idx_O2O2)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_O2O2(bnd_idx_O2O2))) then
            bnd_idx_O2O2=bnd_idx_O2O2
         else
            bnd_idx_O2O2=bnd_idx_O2O2-1
            goto 190
         endif

c     H2OH2O absorption
         call odal_H2OH2O_Chy97(wvn_ctr(bnd_idx),nbr_lev,t,mpl_H2O,RH_liq,odal_H2OH2O,dbg_lvl)

c     Compute O2-O2 absorption optical depth.
         if ((wvl_ctr(bnd_idx).ge.wvl_min_O2O2(bnd_idx_O2O2)).and.
     $        (wvl_ctr(bnd_idx).le.wvl_max_O2O2(bnd_idx_O2O2))) then
            do lev_idx=1,nbr_lev
               odal_O2O2(lev_idx)=
     $              abs_xsx_O2O2(bnd_idx_O2O2)*npl_O2O2(lev_idx)
            enddo               ! end loop over lev
         endif

c     NO2 absorption
 160     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_NO2(1)) then
            bnd_idx_NO2=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_NO2(nbr_NO2_bnd)) then
            bnd_idx_NO2=nbr_NO2_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_NO2(bnd_idx_NO2)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_NO2(bnd_idx_NO2))) then
            bnd_idx_NO2=bnd_idx_NO2
         else
            bnd_idx_NO2=bnd_idx_NO2-1
            goto 160
         endif

c     Compute NO2 absorption optical depth.
         if ((wvl_ctr(bnd_idx).ge.wvl_min_NO2(bnd_idx_NO2)).and.
     $        (wvl_ctr(bnd_idx).le.wvl_max_NO2(bnd_idx_NO2))) then
            do lev_idx=1,nbr_lev
               odal_NO2(lev_idx)=
     $              abs_xsx_NO2(bnd_idx_NO2)*npl_NO2(lev_idx)
            enddo               ! end loop over lev
         endif

c     Liquid cloud properties are allowed to affect all spectral bands. 
c     This is easier than interpolating all Mie parameters onto a regular grid.
 130     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_liq(1)) then
            bnd_idx_liq=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_liq(nbr_liq_bnd)) then
            bnd_idx_liq=nbr_liq_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_liq(bnd_idx_liq)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_liq(bnd_idx_liq))) then
            bnd_idx_liq=bnd_idx_liq
         else
            bnd_idx_liq=bnd_idx_liq-1
            goto 130
         endif

c     Compute liquid Mie scattering/absorption optical depths.
         asm_prm_Mie_liq(bnd_idx)=asm_prm_liq(bnd_idx_liq)
         do lev_idx=1,nbr_lev
            odal_liq(lev_idx)=
     $           abs_cff_mss_liq(bnd_idx_liq)*mpl_LWP(lev_idx)
            odsl_liq(lev_idx)=
     $           sca_cff_mss_liq(bnd_idx_liq)*mpl_LWP(lev_idx)
         enddo                  ! end loop over lev

c     Ice cloud properties are allowed to affect all spectral bands. 
c     This is easier than interpolating all Mie parameters onto a regular grid.
 140     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_ice(1)) then
            bnd_idx_ice=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_ice(nbr_ice_bnd)) then
            bnd_idx_ice=nbr_ice_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_ice(bnd_idx_ice)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_ice(bnd_idx_ice))) then
            bnd_idx_ice=bnd_idx_ice
         else
            bnd_idx_ice=bnd_idx_ice-1
            goto 140
         endif

c     Compute ice Mie scattering/absorption optical depths.
         asm_prm_Mie_ice(bnd_idx)=asm_prm_ice(bnd_idx_ice)
         do lev_idx=1,nbr_lev
            odal_ice(lev_idx)=
     $           abs_cff_mss_ice(bnd_idx_ice)*mpl_IWP(lev_idx)
            odsl_ice(lev_idx)=
     $           sca_cff_mss_ice(bnd_idx_ice)*mpl_IWP(lev_idx)
         enddo                  ! end loop over lev

c     Aerosol properties are allowed to affect all spectral bands. 
c     This is easier than interpolating all Mie parameters onto a regular grid.
 150     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_aer(1)) then
            bnd_idx_aer=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_aer(nbr_aer_bnd)) then
            bnd_idx_aer=nbr_aer_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_aer(bnd_idx_aer)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_aer(bnd_idx_aer))) then
            bnd_idx_aer=bnd_idx_aer
         else
            bnd_idx_aer=bnd_idx_aer-1
            goto 150
         endif

c     Compute aerosol Mie scattering/absorption optical depths.
         asm_prm_Mie_aer(bnd_idx)=asm_prm_aer(bnd_idx_aer)
         do lev_idx=1,nbr_lev
            odal_aer(lev_idx)=
     $           abs_cff_mss_aer(bnd_idx_aer)*mpl_aer(lev_idx)
            odsl_aer(lev_idx)=
     $           sca_cff_mss_aer(bnd_idx_aer)*mpl_aer(lev_idx)
         enddo                  ! end loop over lev

c     Background aerosol properties are allowed to affect all spectral bands. 
c     This is easier than interpolating all Mie parameters onto a regular grid.
 210     continue
         if (wvl_ctr(bnd_idx).lt.wvl_min_bga(1)) then
            bnd_idx_bga=1
         else if (wvl_ctr(bnd_idx).gt.wvl_max_bga(nbr_bga_bnd)) then
            bnd_idx_bga=nbr_bga_bnd
         else if ((wvl_ctr(bnd_idx).ge.wvl_min_bga(bnd_idx_bga)).and.
     $           (wvl_ctr(bnd_idx).le.wvl_max_bga(bnd_idx_bga))) then
            bnd_idx_bga=bnd_idx_bga
         else
            bnd_idx_bga=bnd_idx_bga-1
            goto 210
         endif

c     Compute background aerosol Mie scattering/absorption optical depths.
         asm_prm_Mie_bga(bnd_idx)=asm_prm_bga(bnd_idx_bga)
         do lev_idx=1,nbr_lev
            odal_bga(lev_idx)=
     $           abs_cff_mss_bga(bnd_idx_bga)*mpl_bga(lev_idx)
            odsl_bga(lev_idx)=
     $           sca_cff_mss_bga(bnd_idx_bga)*mpl_bga(lev_idx)
         enddo                  ! end loop over lev

c     Figure out Rayleigh scattering optical depth.
         do lev_idx=1,nbr_lev

c     See Len93 p. 154 for details. See also BrS84 p. 107, GoY89 p. 297.
            idx_rfr_air_STP(bnd_idx)=
     $           1.+
     $           1.e-6*(77.46+.459/(1.e12*wvl_ctr(bnd_idx)**2))*
     $           prs_STP*.01/tpt_STP
            float_foo=32.*(idx_rfr_air_STP(bnd_idx)-1)**2/3.

c     Following line can easily cause overflow in single precision if not handled properly. 
c     Moving factor of N_STP around avoids this problem.
            sca_cff_mss_Ray(lev_idx)=
     $           pi**3*float_foo*(Avagadro/N_STP)/
     $           (mmw_mst_air(lev_idx)*N_STP*wvl_ctr(bnd_idx)**4)
            odsl_Ray(lev_idx)=
     $           mpl_mst_air(lev_idx)*sca_cff_mss_Ray(lev_idx)
         enddo                  ! end loop over lev

c     All individual optical properties for each radiative process
c     have been computed. Look to see if user has specified any
c     processes to be ignored, i.e., compute clear-sky case, pure
c     Rayleigh scattering case, etc. Zero out appropriate 
c     optical constants then proceed normally to calculation 
c     of bulk layer properties. 
c     NB: Mie asymmetry parameter never needs to be zeroed 
c     (e.g., in clear sky case) because it is always weighted 
c     by Mie optical depth before being used in a computation.
         if (.not.flg_Rayleigh) then
            do lev_idx=1,nbr_lev
               odsl_Ray(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no Rayleigh processes
         if (.not.flg_ice) then
            do lev_idx=1,nbr_lev
               odal_ice(lev_idx)=0.
               odsl_ice(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no ice processes
         if (.not.flg_liq) then
            do lev_idx=1,nbr_lev
               odal_liq(lev_idx)=0.
               odsl_liq(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no liq processes
         if (.not.flg_aer) then
            do lev_idx=1,nbr_lev
               odal_aer(lev_idx)=0.
               odsl_aer(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no aer processes
         if (.not.flg_bga) then
            do lev_idx=1,nbr_lev
               odal_bga(lev_idx)=0.
               odsl_bga(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no bga processes
         if (.not.flg_H2O) then
            do lev_idx=1,nbr_lev
               odal_H2O(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no H2O processes
         if (.not.flg_H2OH2O) then
            do lev_idx=1,nbr_lev
               odal_H2OH2O(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no H2OH2O processes
         if (.not.flg_CO2) then
            do lev_idx=1,nbr_lev
               odal_CO2(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no CO2 processes
         if (.not.flg_OH) then
            do lev_idx=1,nbr_lev
               odal_OH(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no OH processes
         if (.not.flg_O2) then
            do lev_idx=1,nbr_lev
               odal_O2(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no O2 processes
         if (.not.flg_O3) then
            do lev_idx=1,nbr_lev
               odal_O3(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no O3 processes
         if (.not.flg_O2O2) then
            do lev_idx=1,nbr_lev
               odal_O2O2(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no O2-O2 processes
         if (.not.flg_NO2) then
            do lev_idx=1,nbr_lev
               odal_NO2(lev_idx)=0.
            enddo               ! end loop over lev
         endif                  ! end if no NO2 processes

c     NB: Now that processes have been turned on/off, save diagnostic values
         if (bnd_idx.eq.bnd_obs_aer) then
            odxl_obs_aer(lev_idx)=odal_aer(lev_idx)+odsl_aer(lev_idx)
         endif                  ! end if bnd_obs_aer
         if (bnd_idx.eq.bnd_obs_bga) then
            odxl_obs_bga(lev_idx)=odal_bga(lev_idx)+odsl_bga(lev_idx)
         endif                  ! end if bnd_obs_bga

c     Weight optical parameters as per CCY83.
         do lev_idx=1,nbr_lev
            odsl_spc_ttl(bnd_idx,lev_idx)=
     $           odsl_Ray(lev_idx)+
     $           odsl_liq(lev_idx)+odsl_ice(lev_idx)+
     $           odsl_aer(lev_idx)+odsl_bga(lev_idx)
            odal_spc_ttl(bnd_idx,lev_idx)=
     $           odal_H2O(lev_idx)+odal_CO2(lev_idx)+
     $           odal_O2(lev_idx)+odal_O3(lev_idx)+
     $           odal_NO2(lev_idx)+odal_OH(lev_idx)+
     $           odal_O2O2(lev_idx)+odal_H2OH2O(lev_idx)+
     $           odal_liq(lev_idx)+odal_ice(lev_idx)+
     $           odal_aer(lev_idx)+odal_bga(lev_idx)
            odxl_spc_ttl(bnd_idx,lev_idx)=
     $           odsl_spc_ttl(bnd_idx,lev_idx)+odal_spc_ttl(bnd_idx,lev_idx)
         enddo                  ! end loop over lev

c     In rare, diagnostic cases where all scattering is turned off,
c     the following lines have a zero over zero condition and
c     will crash under many compilers unless protected by conditions.
         do lev_idx=1,nbr_lev
            if (odsl_spc_ttl(bnd_idx,lev_idx).le.0.) then
               ss_alb_fct(bnd_idx,lev_idx)=0.
               sca_frc_Mie(lev_idx)=0.
               asm_prm_Mie_ttl(bnd_idx,lev_idx)=0.
            else                ! endif no scattering whatsoever

c     Single scattering albedo is ill-conditioned when 
c     scattering optical depth is zero and there is no absorption.
c     This may only be a problem in single precision, i'm not sure.
               ss_alb_fct(bnd_idx,lev_idx)=
     $              odsl_spc_ttl(bnd_idx,lev_idx)/
     $              odxl_spc_ttl(bnd_idx,lev_idx)

c     Mie scattering asymmetry parameter and fraction of
c     scattering due to Mie scattering must be saved for all 
c     bands and levels until DISORT() is called. They are needed to
c     correctly weight the moments of the total phase function 
c     between Rayleigh and Henyey-Greenstein components.
               sca_frc_Mie(lev_idx)=
     $              (odsl_liq(lev_idx)+odsl_ice(lev_idx)+
     $              odsl_aer(lev_idx)+odsl_bga(lev_idx))/
     $              odsl_spc_ttl(bnd_idx,lev_idx)

c     See CZP III p. #115 for discussion of effective asymmetry parameter.
               asm_prm_Mie_ttl(bnd_idx,lev_idx)=
     $              (asm_prm_Mie_liq(bnd_idx)*odsl_liq(lev_idx)+
     $              asm_prm_Mie_ice(bnd_idx)*odsl_ice(lev_idx)+
     $              asm_prm_Mie_aer(bnd_idx)*odsl_aer(lev_idx)+
     $              asm_prm_Mie_bga(bnd_idx)*odsl_bga(lev_idx))/
     $              odsl_spc_ttl(bnd_idx,lev_idx)
            endif               ! endif there is scattering
         enddo                  ! end loop over lev

c     Sanity check for faulty ssalb.
c     omega's as large as 1.0000656 can occur (at least under LINUX) for very
c     small Rayleigh scattering optical depths (p <~ 1 Pa) in single precision
         do lev_idx=1,nbr_lev
            if (ss_alb_fct(bnd_idx,lev_idx).gt.1.) then
               write (6,'(a,i4,a,i3,a,f10.7)') 'WARNING: ss_alb_fct(',bnd_idx,',',lev_idx,') = ',ss_alb_fct(bnd_idx,lev_idx)
               ss_alb_fct(bnd_idx,lev_idx)=1.

            endif
         enddo                  ! end loop over lev
            
c     Compute diagnostics
         do lev_idx=1,nbr_lev
            odxc_spc_OH(bnd_idx)=odxc_spc_OH(bnd_idx)+
     $           odal_OH(lev_idx)
            odxc_spc_O2(bnd_idx)=odxc_spc_O2(bnd_idx)+
     $           odal_O2(lev_idx)
            odxc_spc_CO2(bnd_idx)=odxc_spc_CO2(bnd_idx)+
     $           odal_CO2(lev_idx)
            odxc_spc_H2O(bnd_idx)=odxc_spc_H2O(bnd_idx)+
     $           odal_H2O(lev_idx)
            odxc_spc_H2OH2O(bnd_idx)=odxc_spc_H2OH2O(bnd_idx)+
     $           odal_H2OH2O(lev_idx)
            odxc_spc_O3(bnd_idx)=odxc_spc_O3(bnd_idx)+
     $           odal_O3(lev_idx)
            odxc_spc_O2O2(bnd_idx)=odxc_spc_O2O2(bnd_idx)+
     $           odal_O2O2(lev_idx)
            odxc_spc_NO2(bnd_idx)=odxc_spc_NO2(bnd_idx)+
     $           odal_NO2(lev_idx)
            odxc_spc_ice(bnd_idx)=odxc_spc_ice(bnd_idx)+
     $           odal_ice(lev_idx)+odsl_ice(lev_idx)
            odxc_spc_liq(bnd_idx)=odxc_spc_liq(bnd_idx)+
     $           odal_liq(lev_idx)+odsl_liq(lev_idx)
            odxc_spc_aer(bnd_idx)=odxc_spc_aer(bnd_idx)+
     $           odal_aer(lev_idx)+odsl_aer(lev_idx)
            odxc_spc_bga(bnd_idx)=odxc_spc_bga(bnd_idx)+
     $           odal_bga(lev_idx)+odsl_bga(lev_idx)
            odxc_spc_Ray(bnd_idx)=odxc_spc_Ray(bnd_idx)+
     $           odsl_Ray(lev_idx)
         enddo                  ! end loop over lev
         odxc_spc_ttl(bnd_idx)=
     $        odxc_spc_OH(bnd_idx)+
     $        odxc_spc_O2(bnd_idx)+
     $        odxc_spc_CO2(bnd_idx)+
     $        odxc_spc_H2O(bnd_idx)+
     $        odxc_spc_H2OH2O(bnd_idx)+
     $        odxc_spc_O3(bnd_idx)+
     $        odxc_spc_O2O2(bnd_idx)+
     $        odxc_spc_NO2(bnd_idx)+
     $        odxc_spc_ice(bnd_idx)+
     $        odxc_spc_liq(bnd_idx)+
     $        odxc_spc_aer(bnd_idx)+
     $        odxc_spc_bga(bnd_idx)+
     $        odxc_spc_Ray(bnd_idx)

c     All spectral information needed for DISORT and data reduction has now been saved.
c     Massage data for input to DISORT. 
c     DISORT uses a top down numbering scheme like CCM2.
c     Layer indices from HPMM cirrus cloud model must be reversed.

c     Assign optical depth and single scattering albedos
         do lev_idx=1,nbr_lev
            dtauc(lev_idx)=odxl_spc_ttl(bnd_idx,lev_idx)
            ssalb(lev_idx)=ss_alb_fct(bnd_idx,lev_idx)
         end do

c     Compute phase function moments (currently a weighted mixture
c     of Henyey-Greenstein and Rayleigh phase functions).
         do lev_idx=1,nbr_lev
            pmom(0,lev_idx)=1.
         end do

c     All moments of Rayleigh phase function except second are zero.
c     To blend this with HG phase function from Mie
c     scattering we weight HG coefficients by fraction of
c     scattering due to Mie processes.
         do mom_idx=1,nbr_mom
            do lev_idx=1,nbr_lev
               pmom(mom_idx,lev_idx)=
     $              sca_frc_Mie(lev_idx)*
     $              asm_prm_Mie_ttl(bnd_idx,lev_idx)**mom_idx
            end do
         end do

c     For the second moment, we form a weighted average of 
c     HG scattering fraction times the HG moment (g^2), 
c     plus Rayleigh scattered fraction times Rayleigh moment (.1).
         do lev_idx=1,nbr_lev
            pmom(2,lev_idx)=
     $           (1.-sca_frc_Mie(lev_idx))*.1+
     $           sca_frc_Mie(lev_idx)*
     $           asm_prm_Mie_ttl(bnd_idx,lev_idx)**2
         end do
         if (.not.sv_cmp_opt_dep) then 

c     If usrtau.eq.true. then place levels at top and 
c     bottom of atmosphere and at computational levels in between. 
c     This should be same as setting usrtau=.false.
c     There is not much point in setting user levels to anything
c     else, i.e., evenly spaced in optical depth, or at top and
c     bottom of cloud, yet. But there may be soon. In that case
c     all radiant quantities stored in netCDF will need to be 
c     made dimensioned by nbr_tau instead of nbr_levp.
            tau(1)=0.
            float_foo=0.
            do tau_idx=1,nbr_tau
               float_foo=float_foo+odxl_spc_ttl(bnd_idx,tau_idx)
               tau(tau_idx+1)=float_foo
            end do              ! end loop over tau

            ntau=nbr_tau
            do tau_idx=1,nbr_tau
               utau(tau_idx)=tau(tau_idx)
            end do              ! end loop over tau
         endif                  ! end if setting user defined levels

c     Set wavenumbers for this spectral interval.
         wvnmlo=wvn_min(bnd_idx)
         wvnmhi=wvn_max(bnd_idx)

c     Set plank=.true. whenever considering thermal emission.
c     Computing thermal source function generates underflows in single precision when lambda < XXX microns.
         if (wvl_ctr(bnd_idx).gt.wvl_Planck) then
            plank=flg_Planck
         else 
            plank=.false.
         endif

c     Set intensity of incident parallel beam at top boundary.
c     Units are arbitrary but must match fisot in solar case.
c     Must be in W/m^2 if there is any thermal emission (plank=.true.)
c     Recall that fbeam is flux normal to earth-sun path,
c     not normal to ground, i.e., zenith angle correction
c     to incoming flux is made within DISORT() so DO NOT
c     pre-multiply fbeam by cosine of the solar zenith angle.
         fbeam=slr_cst_xnt_fac*flx_slr_frc(bnd_idx)

c     Give details about lower boundary reflectance.
c     When lamber is true, albedo specifies isotropic reflectance.
c     Otherwise hl array must be specified to give BDRF of bottom boundary.
         lamber=.true.
         if (wvl_ctr(bnd_idx).lt..7e-6) then
            if (slr_zen_ngl_cos.gt..5) then
               albedo=alb_vis_str_zen_ngl
            else
               albedo=alb_vis_weak_zen_ngl
            endif
         else
            if (slr_zen_ngl_cos.gt..5) then
               albedo=alb_NIR_str_zen_ngl
            else
               albedo=alb_NIR_weak_zen_ngl
            endif
         endif
      
         if (single_bnd_computation) then

c     Set header which DISORT will use printing results.
c     If header has length greater than zero, DISORT() prints
c     annoying header message on each call.
            write (header,'(a,i4,a,f6.4,a,e8.1,a,e8.1)')
     $           'lbl: i = ',bnd_idx,', lambda = ',
     $           wvl_ctr(bnd_idx)*1.e6,
     $           ', flx_slr_frc = ',flx_slr_frc(bnd_idx),', bandwith = ',
     $           wvl_sz(bnd_idx)

            if (bnd_idx.eq.bnd_dbg) then

c     Hardcode numbers to test radiative code against benchmark values.
               if (tst_case_Rayleigh) then

c     Pure Rayleigh scattering. Answer: Albedo = .796920
                  fbeam=pi
                  albedo=0.
                  umu0=.5
                  do lev_idx=1,nbr_lev
                     dtauc(lev_idx)=4.284/nbr_lev
                     odsl_Ray(lev_idx)=dtauc(lev_idx)
                     odsl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     odxl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                  enddo         ! end loop over lev
               endif            ! end if tst_case_Rayleigh
               if (tst_case_HG) then

c     Pure Henyey-Greenstein scattering. Answer: Albedo = .123420
c     Atmospheric absorptance = .360522
                  fbeam=pi
                  albedo=0.
                  umu0=.5
                  do lev_idx=1,nbr_lev
                     dtauc(lev_idx)=1./nbr_lev
                     ssalb(lev_idx)=.8
                     odsl_Ray(lev_idx)=dtauc(lev_idx)
                     odsl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     odxl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     ss_alb_fct(bnd_idx,lev_idx)=ssalb(lev_idx)
                     sca_frc_Mie(lev_idx)=1.
                     asm_prm_Mie_liq(bnd_idx)=.75
                     asm_prm_Mie_ttl(bnd_idx,lev_idx)=.75
                     odsl_Ray(lev_idx)=0.
                  enddo         ! end loop over lev
                  do lev_idx=1,nbr_lev
                     pmom(0,lev_idx)=1.
                  end do
                  do mom_idx=1,nbr_mom
                     do lev_idx=1,nbr_lev
                        pmom(mom_idx,lev_idx)=.75**mom_idx
                     end do
                  end do
               endif            ! end if tst_case_HG
               write (6,'(/,a,i4,a4,f9.6,a8)')
     $              'Detailed debugging information for bnd(',
     $              bnd_idx,') = ',wvl_ctr(bnd_idx)*1.e6,' microns'
               write (6,'(a4,16(x,a8))')
     $              'lev_idx',
     $              'tau_sca',
     $              'tau_abs',
     $              'tau_ext',
     $              'ss_alb',
     $              'Mie_frc',
     $              'liq_asym',
     $              'tot_asym',
     $              'Ray_msig',
     $              'pth_wet',
     $              'Ray_tau',
     $              'dtauc',
     $              'ssalb',
     $              'pmom(1)',
     $              'pmom(2)',
     $              'pmom(3)',
     $              'pmom(4)'
               do lev_idx=1,nbr_lev
                  write (6,'(i4,16(x,1p,e8.1))')
     $                 lev_idx,
     $                 odsl_spc_ttl(bnd_idx,lev_idx),
     $                 odal_spc_ttl(bnd_idx,lev_idx),
     $                 odxl_spc_ttl(bnd_idx,lev_idx),
     $                 ss_alb_fct(bnd_idx,lev_idx),
     $                 sca_frc_Mie(lev_idx),
     $                 asm_prm_Mie_liq(bnd_idx),
     $                 asm_prm_Mie_ttl(bnd_idx,lev_idx),
     $                 sca_cff_mss_Ray(lev_idx),
     $                 mpl_mst_air(lev_idx),
     $                 odsl_Ray(lev_idx),
     $                 dtauc(lev_idx),
     $                 ssalb(lev_idx),
     $                 pmom(1,lev_idx),
     $                 pmom(2,lev_idx),
     $                 pmom(3,lev_idx),
     $                 pmom(4,lev_idx)
               enddo            ! end loop over lev
            else                ! end if dbg band
               goto 999
            endif               ! end if not dbg band
         endif                  ! end if single_bnd_computation

         call DISORT( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $        wvnmhi, usrtau, ntau, utau, nstr, usrang,
     $        numu, umu, nphi, phi, ibcnd, fbeam, umu0,
     $        phi0, fisot, lamber, albedo, hl, btemp, ttemp,
     $        temis, deltam, plank, onlyfl, accur, prnt,
     $        header, maxcly, maxulv, maxumu, maxcmu,
     $        maxphi, rfldir, rfldn, flup, dfdt, uavg,
     $        uu, u0u, albmed, trnmed )

c     Store spectral fluxes returned from DISORT.
c     DISORT uses a funny indexing system: 
c     Returned flux arrays, which are on layer interfaces, 
c     are always indexed starting from 1 at TOA even though same levels 
c     start with index 0 on input to DISORT, i.e., temper().
c     Make sure you know what you're doing before changing these indices.
         do lev_idx=1,nbr_lev
            ntn_spc_mean(bnd_idx,lev_idx)=
     $           .5*(uavg(lev_idx)+uavg(lev_idx+1))/
     $           wvl_sz(bnd_idx)
         enddo                  ! end loop over lev
         do lev_idx=1,nbr_levp
            flx_spc_dwn_drc(bnd_idx,lev_idx)=rfldir(lev_idx)/
     $           wvl_sz(bnd_idx)
            flx_spc_dwn_dff(bnd_idx,lev_idx)=rfldn(lev_idx)/
     $           wvl_sz(bnd_idx)
            flx_spc_dwn(bnd_idx,lev_idx)=
     $           flx_spc_dwn_drc(bnd_idx,lev_idx)+
     $           flx_spc_dwn_dff(bnd_idx,lev_idx)
            flx_spc_up(bnd_idx,lev_idx)=flup(lev_idx)/
     $           wvl_sz(bnd_idx)
            flx_spc_net(bnd_idx,lev_idx)=
     $           flx_spc_dwn(bnd_idx,lev_idx)-flx_spc_up(bnd_idx,lev_idx)
         enddo                  ! end loop over lev

c     u0u returns azimuthally averaged intensities at user polar angles. 
c     uu returns exact intensities at each user azimuthal angle at each user polar angle.
c     User polar angles are ordered by increasing cosine of polar angle, 
c     i.e., from downwelling radiance to upwelling radiance. 

c     NB: Possible bug in current DISORT() documentation. 
c     Intensities are returned in arrays in order of descending cosine of polar angle.
         ntn_spc_aa_ndr_sfc(bnd_idx)=u0u(1,nbr_levp)/
     $        wvl_sz(bnd_idx)
         ntn_spc_aa_zen_sfc(bnd_idx)=u0u(nbr_plr,nbr_levp)/
     $        wvl_sz(bnd_idx)

         do lev_idx=1,nbr_levp
            do plr_idx=1,nbr_plr
               ntn_bb_aa(plr_idx,lev_idx)=ntn_bb_aa(plr_idx,lev_idx)+
     $              u0u(plr_idx,lev_idx)
            enddo               ! end loop over plr
            ntn_spc_aa_ndr(bnd_idx,lev_idx)=
     $           u0u(1,lev_idx)/wvl_sz(bnd_idx)
            ntn_spc_aa_zen(bnd_idx,lev_idx)=
     $           u0u(nbr_plr,lev_idx)/wvl_sz(bnd_idx)
         enddo                  ! end loop over lev

         do plr_idx=1,nbr_plr
            ntn_spc_aa_sfc(plr_idx,bnd_idx)=
     $           u0u(plr_idx,nbr_levp)/wvl_sz(bnd_idx)
         enddo                  ! end loop over lev

         if (bnd_idx.eq.bnd_dbg) then
            do azi_idx=1,nbr_azi
               do plr_idx=1,nbr_plr
                  do lev_idx=1,nbr_levp
                     ntn_spc_chn(azi_idx,plr_idx,lev_idx)=
     $                    uu(plr_idx,lev_idx,azi_idx)/wvl_sz(bnd_idx)
                  enddo         ! end loop over lev
               enddo            ! end loop over plr
            enddo               ! end loop over azimuthal angles
         endif                  ! end if dbg band

         if (sv_ntn) then
c            do azi_idx=1,nbr_azi
c               do plr_idx=1,nbr_plr
c                  do lev_idx=1,nbr_levp
c                     ntn_spc_aa(plr_idx,bnd_idx,lev_idx)=
c     $                    u0u(plr_idx,lev_idx)/wvl_sz(bnd_idx)
c                  enddo         ! end loop over lev
c               enddo            ! end loop over plr
c            enddo               ! end loop over azimuthal angles
         endif                  ! end if saving full intensity arrays

 999     continue

      enddo                     ! end loop over bnd

c     End section 2: Main computation loop over all bands
c     Begin section 3: Postprocessing
c     computation of diagnostic arrays from output storage arrays
c     summation of spectral diagnostic arrays to integrated arrays
c     reduction of integrated arrays to diagnostic scalars
c     netCDF output

      write (6,'(/)')

c     If we are storing values at computational angles then
c     they should have been returned in umu array.
      if (sv_cmp_plr_ngl) then

         do plr_idx=1,nbr_plr
            plr_cos(plr_idx)=umu(plr_idx)
            plr(plr_idx)=acos(plr_cos(plr_idx))
            plr_dgr(plr_idx)=180.*plr(plr_idx)/pi
         enddo                  ! end loop over plr
      endif                     ! endif looking at computational angles

      do bnd_idx=1,nbr_bnd
         trn_spc_atm_CO2(bnd_idx)=exp(-min(odxc_spc_CO2(bnd_idx),25.))
         trn_spc_atm_H2O(bnd_idx)=exp(-min(odxc_spc_H2O(bnd_idx),25.))
         trn_spc_atm_H2OH2O(bnd_idx)=exp(-min(odxc_spc_H2OH2O(bnd_idx),25.))
         trn_spc_atm_ice(bnd_idx)=exp(-min(odxc_spc_ice(bnd_idx),25.))
         trn_spc_atm_liq(bnd_idx)=exp(-min(odxc_spc_liq(bnd_idx),25.))
         trn_spc_atm_aer(bnd_idx)=exp(-min(odxc_spc_aer(bnd_idx),25.))
         trn_spc_atm_bga(bnd_idx)=exp(-min(odxc_spc_bga(bnd_idx),25.))
         trn_spc_atm_OH(bnd_idx)=exp(-min(odxc_spc_OH(bnd_idx),25.))
         trn_spc_atm_O2(bnd_idx)=exp(-min(odxc_spc_O2(bnd_idx),25.))
         trn_spc_atm_O3(bnd_idx)=exp(-min(odxc_spc_O3(bnd_idx),25.))
         trn_spc_atm_O2O2(bnd_idx)=exp(-min(odxc_spc_O2O2(bnd_idx),25.))
         trn_spc_atm_NO2(bnd_idx)=exp(-min(odxc_spc_NO2(bnd_idx),25.))
         trn_spc_atm_Ray(bnd_idx)=exp(-min(odxc_spc_Ray(bnd_idx),25.))
      enddo                     ! end loop over bnd
      do bnd_idx=1,nbr_bnd
         flx_spc_dwn_TOA(bnd_idx)=flx_spc_dwn(bnd_idx,1)
         flx_spc_dwn_sfc(bnd_idx)=flx_spc_dwn(bnd_idx,nbr_levp)

c     Compute absorbed spectral fluxes.
         flx_spc_abs_SAS(bnd_idx)=flx_spc_net(bnd_idx,1)
         flx_spc_abs_sfc(bnd_idx)=flx_spc_net(bnd_idx,nbr_levp)
         flx_spc_abs_atm(bnd_idx)=
     $        flx_spc_abs_SAS(bnd_idx)-flx_spc_abs_sfc(bnd_idx)

         do lev_idx=1,nbr_lev
            flx_spc_abs(bnd_idx,lev_idx)=
     $           max(0.,
     $           flx_spc_net(bnd_idx,lev_idx)-flx_spc_net(bnd_idx,lev_idx+1))
         enddo                  ! end loop over lev
      enddo                     ! end loop over bnd

c     Compute running sum of atmospheric absorption
c     NB:  flx_abs_atm_rdr is meant to hold total atmospheric absorption
c     occuring redward (at longer wavelengths) of a given band. 
c     This quantity, at certain wavelengths, can be backed out of 
c     fractional spectral instruments, like Valero's FSBR.
      flx_abs_atm_rdr(1)=flx_spc_abs_atm(1)*wvl_sz(1)
      do bnd_idx=2,nbr_bnd
         flx_abs_atm_rdr(bnd_idx)=flx_abs_atm_rdr(bnd_idx-1)+
     $        flx_spc_abs_atm(bnd_idx)*wvl_sz(bnd_idx)
      enddo                     ! end loop over bnd

c     Instrument fluxes are a paradigm useful for any fractional band measurement. 
c     Therefore we compute all instrument diagnostics in one location.
c     The computation of instrument diagnostics is perfectly analogous to 
c     broadband (bb) diagnostics; reuse code by changing bb to nst.
c     All that need be changed to adapt to a new instrument are the
c     lower and upper limits of the (rectangular window) bandpass. 
c     Initialize basic instrument quantities which will be incremented.
      do lev_idx=1,nbr_levp
         flx_nst_dwn(lev_idx)=0.
         flx_nst_up(lev_idx)=0.
      enddo                     ! end loop over lev
c     Initialize a default instrument spectral response function
      do bnd_idx=1,nbr_bnd
         bnd_wgt(bnd_idx)=1.
      enddo                     ! end loop over bnd
      if (flt_nst) then
c     Instrument filter response affects all spectral bands. 
         bnd_idx_nst=nbr_nst_bnd
         do bnd_idx=1,nbr_bnd
 170        continue
            if (wvl_ctr(bnd_idx).lt.wvl_min_nst(1)) then
               bnd_idx_nst=1
            else if (wvl_ctr(bnd_idx).gt.wvl_max_nst(nbr_nst_bnd)) then
               bnd_idx_nst=nbr_nst_bnd
            else if ((wvl_ctr(bnd_idx).ge.wvl_min_nst(bnd_idx_nst)).and.
     $              (wvl_ctr(bnd_idx).le.wvl_max_nst(bnd_idx_nst))) then
               bnd_idx_nst=bnd_idx_nst
            else
               bnd_idx_nst=bnd_idx_nst-1
               goto 170
            endif
            bnd_wgt(bnd_idx)=nst_SRF(bnd_idx_nst)
         enddo                  ! end loop over bnd
      endif                     ! endif flt_nst
      do bnd_idx=1,nbr_bnd
         do lev_idx=1,nbr_levp
            flx_nst_dwn(lev_idx)=flx_nst_dwn(lev_idx)+
     $           flx_spc_dwn(bnd_idx,lev_idx)*wvl_sz(bnd_idx)*bnd_wgt(bnd_idx)
            flx_nst_up(lev_idx)=flx_nst_up(lev_idx)+
     $           flx_spc_up(bnd_idx,lev_idx)*wvl_sz(bnd_idx)*bnd_wgt(bnd_idx)
         enddo                  ! end loop over lev
      enddo                     ! end loop over bnd
      do lev_idx=1,nbr_levp
         flx_nst_net(lev_idx)=flx_nst_dwn(lev_idx)-flx_nst_up(lev_idx)
      enddo                     ! end loop over lev
      do lev_idx=1,nbr_lev
         flx_nst_abs(lev_idx)=flx_nst_net(lev_idx)-flx_nst_net(lev_idx+1)
      enddo                  ! end loop over lev
      if (flx_nst_dwn(1).gt.0.) then 
         abs_nst_SAS=flx_nst_net(1)/flx_nst_dwn(1)
         abs_nst_atm=(flx_nst_net(1)-flx_nst_net(nbr_levp))/flx_nst_dwn(1)
         abs_nst_sfc=flx_nst_net(nbr_levp)/flx_nst_dwn(1)
         rfl_nst_SAS=flx_nst_up(1)/flx_nst_dwn(1)
         rfl_nst_sfc=flx_nst_up(nbr_levp)/flx_nst_dwn(nbr_levp)
         trn_nst_atm=flx_nst_dwn(nbr_levp)/flx_nst_dwn(1)
      else
         abs_nst_SAS=0.
         abs_nst_atm=0.
         abs_nst_sfc=0.
         rfl_nst_SAS=0.
         rfl_nst_sfc=0.
         trn_nst_atm=0.
      endif
      flx_nst_abs_atm=flx_nst_net(1)-flx_nst_net(nbr_levp)
      flx_nst_abs_sfc=flx_nst_net(nbr_levp)
      flx_nst_abs_ttl=flx_nst_net(1)
      flx_nst_dwn_TOA=flx_nst_dwn(1)
      flx_nst_dwn_sfc=flx_nst_dwn(nbr_levp)
c     End instrument computations

c     Define system transmittance, reflectance, and absorptance.
c     Because these quantities are all normalized by the insolation,
c     which may be zero, these definitions must be made in a 
c     conditional clause.
      do bnd_idx=1,nbr_bnd
         if (flx_spc_dwn_TOA(bnd_idx).gt.0) then
            trn_spc_atm_ttl(bnd_idx)=
     $           flx_spc_dwn_sfc(bnd_idx)/
     $           flx_spc_dwn_TOA(bnd_idx)
            rfl_spc_SAS(bnd_idx)=
     $           flx_spc_up(bnd_idx,1)/
     $           flx_spc_dwn_TOA(bnd_idx)

c     Layer absorptance is absorbed flux in a layer normalized by total flux entering layer. 
c     Define absorptance so surface + atmospheric absorptances sum to total SAS absorptance,
c     i.e., as fraction of solar insolation absorbed by atmosphere, surface, and SAS, respectively.
            abs_spc_SAS(bnd_idx)=
     $           flx_spc_net(bnd_idx,1)/flx_spc_dwn_TOA(bnd_idx)
            abs_spc_sfc(bnd_idx)=
     $           flx_spc_net(bnd_idx,nbr_levp)/flx_spc_dwn_TOA(bnd_idx)
         else
            trn_spc_atm_ttl(bnd_idx)=0.
            rfl_spc_SAS(bnd_idx)=0.
            abs_spc_SAS(bnd_idx)=0.
            abs_spc_sfc(bnd_idx)=0.
         endif
      enddo                     ! end loop over bnd
      do bnd_idx=1,nbr_bnd
         abs_spc_atm(bnd_idx)=abs_spc_SAS(bnd_idx)-abs_spc_sfc(bnd_idx)
      enddo                     ! end loop over bnd

      do bnd_idx=1,nbr_bnd

c     Accumulate spectral fluxes into broadband arrays
         do lev_idx=1,nbr_lev
            ntn_bb_mean(lev_idx)=ntn_bb_mean(lev_idx)+
     $           ntn_spc_mean(bnd_idx,lev_idx)*
     $           wvl_sz(bnd_idx)
         enddo                  ! end loop over lev
         do lev_idx=1,nbr_levp
            flx_bb_dwn_drc(lev_idx)=flx_bb_dwn_drc(lev_idx)+
     $           flx_spc_dwn_drc(bnd_idx,lev_idx)*
     $           wvl_sz(bnd_idx)
            flx_bb_dwn_dff(lev_idx)=flx_bb_dwn_dff(lev_idx)+
     $           flx_spc_dwn_dff(bnd_idx,lev_idx)*
     $           wvl_sz(bnd_idx)
            flx_bb_up(lev_idx)=flx_bb_up(lev_idx)+
     $           flx_spc_up(bnd_idx,lev_idx)*
     $           wvl_sz(bnd_idx)
         enddo                  ! end loop over lev
      enddo                  ! end loop over bnd

c     Process broadband fluxes
      do lev_idx=1,nbr_levp
         flx_bb_dwn(lev_idx)=flx_bb_dwn_drc(lev_idx)+flx_bb_dwn_dff(lev_idx)
         flx_bb_net(lev_idx)=flx_bb_dwn(lev_idx)-flx_bb_up(lev_idx)
      enddo                  ! end loop over lev

c     Compute scalar diagnostics
      abs_bb_SAS=flx_bb_net(1)/flx_bb_dwn(1)
      abs_bb_atm=(flx_bb_net(1)-flx_bb_net(nbr_levp))/flx_bb_dwn(1)
      abs_bb_sfc=flx_bb_net(nbr_levp)/flx_bb_dwn(1)
      rfl_bb_SAS=flx_bb_up(1)/flx_bb_dwn(1)
      rfl_bb_sfc=flx_bb_up(nbr_levp)/flx_bb_dwn(nbr_levp)
      trn_bb_atm=flx_bb_dwn(nbr_levp)/flx_bb_dwn(1)
      flx_bb_abs_ttl=flx_bb_net(1)
      flx_bb_abs_sfc=flx_bb_net(nbr_levp)
      flx_bb_abs_atm=flx_bb_net(1)-flx_bb_net(nbr_levp)
      flx_bb_dwn_TOA=flx_bb_dwn(1)
      flx_bb_dwn_sfc=flx_bb_dwn(nbr_levp)

      do lev_idx=1,nbr_lev
         flx_bb_abs(lev_idx)=flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1)
         htg_rate_bb(lev_idx)=
     $        (flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1))*grv(lev_idx)/
     $        (spec_heat_mst_air(lev_idx)*prs_dlt(lev_idx))
      enddo                  ! end loop over lev
c      write (6,'(i4,a,e8.1)') bnd_idx,' flurt flx_spc_dwn(1690,nbr_levp) =  ',flx_spc_dwn(1690,nbr_levp)

c     Begin netCDF output routines
      rcd=nccre(fl_out,NF_CLOBBER,nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_out)

c     Define dimension IDs
      azi_dim_id=ncddef(nc_id,'azi',nbr_azi,rcd)
      bnd_dim_id=ncddef(nc_id,'bnd',nbr_bnd,rcd)
      lev_dim_id=ncddef(nc_id,'lev',nbr_lev,rcd)
      levp_dim_id=ncddef(nc_id,'levp',nbr_levp,rcd)
      plr_dim_id=ncddef(nc_id,'plr',nbr_plr,rcd)
      tau_dim_id=ncddef(nc_id,'tau',nbr_tau,rcd)

c     Assemble ID and count vectors for each multidimensional combination of dimensions
      cnt_bnd_lev(1)=nbr_bnd
      cnt_bnd_lev(2)=nbr_lev
      dim_bnd_lev(1)=bnd_dim_id
      dim_bnd_lev(2)=lev_dim_id

      cnt_bnd_levp(1)=nbr_bnd
      cnt_bnd_levp(2)=nbr_levp
      dim_bnd_levp(1)=bnd_dim_id
      dim_bnd_levp(2)=levp_dim_id

      cnt_plr_levp(1)=nbr_plr
      cnt_plr_levp(2)=nbr_levp
      dim_plr_levp(1)=plr_dim_id
      dim_plr_levp(2)=levp_dim_id

      cnt_plr_bnd(1)=nbr_plr
      cnt_plr_bnd(2)=nbr_bnd
      dim_plr_bnd(1)=plr_dim_id
      dim_plr_bnd(2)=bnd_dim_id

      cnt_plr_bnd_levp(1)=nbr_plr
      cnt_plr_bnd_levp(2)=nbr_bnd
      cnt_plr_bnd_levp(3)=nbr_levp
      dim_plr_bnd_levp(1)=plr_dim_id
      dim_plr_bnd_levp(2)=bnd_dim_id
      dim_plr_bnd_levp(3)=levp_dim_id

      cnt_azi_plr_levp(1)=nbr_azi
      cnt_azi_plr_levp(2)=nbr_plr
      cnt_azi_plr_levp(3)=nbr_levp
      dim_azi_plr_levp(1)=azi_dim_id
      dim_azi_plr_levp(2)=plr_dim_id
      dim_azi_plr_levp(3)=levp_dim_id

      cnt_azi_plr_bnd_levp(1)=nbr_azi
      cnt_azi_plr_bnd_levp(2)=nbr_plr
      cnt_azi_plr_bnd_levp(3)=nbr_bnd
      cnt_azi_plr_bnd_levp(4)=nbr_levp
      dim_azi_plr_bnd_levp(1)=azi_dim_id
      dim_azi_plr_bnd_levp(2)=plr_dim_id
      dim_azi_plr_bnd_levp(3)=bnd_dim_id
      dim_azi_plr_bnd_levp(4)=levp_dim_id

c     Variable definitions.
      abs_bb_SAS_id=ncvdef(nc_id,'abs_bb_SAS',NF_REAL,0,0,rcd)
      abs_bb_atm_id=ncvdef(nc_id,'abs_bb_atm',NF_REAL,0,0,rcd)
      abs_bb_sfc_id=ncvdef(nc_id,'abs_bb_sfc',NF_REAL,0,0,rcd)
      abs_nst_SAS_id=ncvdef(nc_id,'abs_nst_SAS',NF_REAL,0,0,rcd)
      abs_nst_atm_id=ncvdef(nc_id,'abs_nst_atm',NF_REAL,0,0,rcd)
      abs_nst_sfc_id=ncvdef(nc_id,'abs_nst_sfc',NF_REAL,0,0,rcd)
      abs_spc_SAS_id=ncvdef(nc_id,'abs_spc_SAS',NF_REAL,1,bnd_dim_id,rcd)
      abs_spc_atm_id=ncvdef(nc_id,'abs_spc_atm',NF_REAL,1,bnd_dim_id,rcd)
      abs_spc_sfc_id=ncvdef(nc_id,'abs_spc_sfc',NF_REAL,1,bnd_dim_id,rcd)
      alb_sfc_id=ncvdef(nc_id,'alb_sfc',NF_REAL,0,0,rcd)
      azi_dgr_id=ncvdef(nc_id,'azi_dgr',NF_REAL,1,azi_dim_id,rcd)
      azi_id=ncvdef(nc_id,'azi',NF_REAL,1,azi_dim_id,rcd)
      bnd_id=ncvdef(nc_id,'bnd',NF_REAL,1,bnd_dim_id,rcd)
      slr_zen_ngl_cos_id=ncvdef(nc_id,'slr_zen_ngl_cos',NF_DOUBLE,0,0,rcd)
      flx_abs_atm_rdr_id=ncvdef(nc_id,'flx_abs_atm_rdr',NF_REAL,1,bnd_dim_id,rcd)
      flx_bb_abs_atm_id=ncvdef(nc_id,'flx_bb_abs_atm',NF_REAL,0,0,rcd)
      flx_bb_abs_id=ncvdef(nc_id,'flx_bb_abs',NF_REAL,1,lev_dim_id,rcd)
      flx_bb_abs_sfc_id=ncvdef(nc_id,'flx_bb_abs_sfc',NF_REAL,0,0,rcd)
      flx_bb_abs_ttl_id=ncvdef(nc_id,'flx_bb_abs_ttl',NF_REAL,0,0,rcd)
      flx_bb_dwn_TOA_id=ncvdef(nc_id,'flx_bb_dwn_TOA',NF_REAL,0,0,rcd)
      flx_bb_dwn_dff_id=ncvdef(nc_id,'flx_bb_dwn_dff',NF_REAL,1,levp_dim_id,rcd)
      flx_bb_dwn_drc_id=ncvdef(nc_id,'flx_bb_dwn_drc',NF_REAL,1,levp_dim_id,rcd)
      flx_bb_dwn_id=ncvdef(nc_id,'flx_bb_dwn',NF_REAL,1,levp_dim_id,rcd)
      flx_bb_dwn_sfc_id=ncvdef(nc_id,'flx_bb_dwn_sfc',NF_REAL,0,0,rcd)
      flx_bb_net_id=ncvdef(nc_id,'flx_bb_net',NF_REAL,1,levp_dim_id,rcd)
      flx_bb_up_id=ncvdef(nc_id,'flx_bb_up',NF_REAL,1,levp_dim_id,rcd)
      flx_nst_abs_atm_id=ncvdef(nc_id,'flx_nst_abs_atm',NF_REAL,0,0,rcd)
      flx_nst_abs_id=ncvdef(nc_id,'flx_nst_abs',NF_REAL,1,lev_dim_id,rcd)
      flx_nst_abs_sfc_id=ncvdef(nc_id,'flx_nst_abs_sfc',NF_REAL,0,0,rcd)
      flx_nst_abs_ttl_id=ncvdef(nc_id,'flx_nst_abs_ttl',NF_REAL,0,0,rcd)
      flx_nst_dwn_TOA_id=ncvdef(nc_id,'flx_nst_dwn_TOA',NF_REAL,0,0,rcd)
      flx_nst_dwn_id=ncvdef(nc_id,'flx_nst_dwn',NF_REAL,1,levp_dim_id,rcd)
      flx_nst_dwn_sfc_id=ncvdef(nc_id,'flx_nst_dwn_sfc',NF_REAL,0,0,rcd)
      flx_nst_net_id=ncvdef(nc_id,'flx_nst_net',NF_REAL,1,levp_dim_id,rcd)
      flx_nst_up_id=ncvdef(nc_id,'flx_nst_up',NF_REAL,1,levp_dim_id,rcd)
      flx_spc_abs_SAS_id=ncvdef(nc_id,'flx_spc_abs_SAS',NF_REAL,1,bnd_dim_id,rcd)
      flx_spc_abs_atm_id=ncvdef(nc_id,'flx_spc_abs_atm',NF_REAL,1,bnd_dim_id,rcd)
      flx_spc_abs_id=ncvdef(nc_id,'flx_spc_abs',NF_REAL,2,dim_bnd_lev,rcd)
      flx_spc_abs_sfc_id=ncvdef(nc_id,'flx_spc_abs_sfc',NF_REAL,1,bnd_dim_id,rcd)
      flx_spc_dwn_TOA_id=ncvdef(nc_id,'flx_spc_dwn_TOA',NF_REAL,1,bnd_dim_id,rcd)
      flx_spc_dwn_dff_id=ncvdef(nc_id,'flx_spc_dwn_dff',NF_REAL,2,dim_bnd_levp,rcd)
      flx_spc_dwn_drc_id=ncvdef(nc_id,'flx_spc_dwn_drc',NF_REAL,2,dim_bnd_levp,rcd)
      flx_spc_dwn_id=ncvdef(nc_id,'flx_spc_dwn',NF_REAL,2,dim_bnd_levp,rcd)
      flx_spc_dwn_sfc_id=ncvdef(nc_id,'flx_spc_dwn_sfc',NF_REAL,1,bnd_dim_id,rcd)
      flx_spc_up_id=ncvdef(nc_id,'flx_spc_up',NF_REAL,2,dim_bnd_levp,rcd)
      frc_ice_ttl_id=ncvdef(nc_id,'frc_ice_ttl',NF_REAL,0,0,rcd)
      flx_slr_frc_id=ncvdef(nc_id,'flx_slr_frc',NF_REAL,1,bnd_dim_id,rcd)
      htg_rate_bb_id=ncvdef(nc_id,'htg_rate_bb',NF_REAL,1,lev_dim_id,rcd)
      lat_dgr_id=ncvdef(nc_id,'lat_dgr',NF_DOUBLE,0,0,rcd)
      lcl_time_hr_id=ncvdef(nc_id,'lcl_time_hr',NF_DOUBLE,0,0,rcd)
      lcl_yr_day_id=ncvdef(nc_id,'lcl_yr_day',NF_DOUBLE,0,0,rcd)
      lev_id=ncvdef(nc_id,'lev',NF_REAL,1,lev_dim_id,rcd)
      levp_id=ncvdef(nc_id,'levp',NF_REAL,1,levp_dim_id,rcd)
      mpc_CWP_id=ncvdef(nc_id,'mpc_CWP',NF_REAL,0,0,rcd)
      ntn_bb_aa_id=ncvdef(nc_id,'ntn_bb_aa',NF_REAL,2,dim_plr_levp,rcd)
      ntn_bb_mean_id=ncvdef(nc_id,'ntn_bb_mean',NF_REAL,1,lev_dim_id,rcd)
      ntn_spc_aa_ndr_id=ncvdef(nc_id,'ntn_spc_aa_ndr',NF_REAL,2,dim_bnd_levp,rcd)
      ntn_spc_aa_ndr_sfc_id=ncvdef(nc_id,'ntn_spc_aa_ndr_sfc',NF_REAL,1,bnd_dim_id,rcd)
      ntn_spc_aa_sfc_id=ncvdef(nc_id,'ntn_spc_aa_sfc',NF_REAL,2,dim_plr_bnd,rcd)
      ntn_spc_aa_zen_id=ncvdef(nc_id,'ntn_spc_aa_zen',NF_REAL,2,dim_bnd_levp,rcd)
      ntn_spc_aa_zen_sfc_id=ncvdef(nc_id,'ntn_spc_aa_zen_sfc',NF_REAL,1,bnd_dim_id,rcd)
      ntn_spc_chn_id=ncvdef(nc_id,'ntn_spc_chn',NF_REAL,3,dim_azi_plr_levp,rcd)
      ntn_spc_mean_id=ncvdef(nc_id,'ntn_spc_mean',NF_REAL,2,dim_bnd_lev,rcd)
      odxc_obs_aer_id=ncvdef(nc_id,'odxc_obs_aer',NF_REAL,0,0,rcd)
      odxc_obs_bga_id=ncvdef(nc_id,'odxc_obs_bga',NF_REAL,0,0,rcd)
      odxl_obs_aer_id=ncvdef(nc_id,'odxl_obs_aer',NF_REAL,1,lev_dim_id,rcd)
      odxl_obs_bga_id=ncvdef(nc_id,'odxl_obs_bga',NF_REAL,1,lev_dim_id,rcd)
      odxc_spc_CO2_id=ncvdef(nc_id,'odxc_spc_CO2',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_H2OH2O_id=ncvdef(nc_id,'odxc_spc_H2OH2O',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_H2O_id=ncvdef(nc_id,'odxc_spc_H2O',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_NO2_id=ncvdef(nc_id,'odxc_spc_NO2',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_O2_id=ncvdef(nc_id,'odxc_spc_O2',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_O3_id=ncvdef(nc_id,'odxc_spc_O3',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_O2O2_id=ncvdef(nc_id,'odxc_spc_O2O2',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_OH_id=ncvdef(nc_id,'odxc_spc_OH',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_Ray_id=ncvdef(nc_id,'odxc_spc_Ray',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_aer_id=ncvdef(nc_id,'odxc_spc_aer',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_bga_id=ncvdef(nc_id,'odxc_spc_bga',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_ice_id=ncvdef(nc_id,'odxc_spc_ice',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_liq_id=ncvdef(nc_id,'odxc_spc_liq',NF_REAL,1,bnd_dim_id,rcd)
      odxc_spc_ttl_id=ncvdef(nc_id,'odxc_spc_ttl',NF_REAL,1,bnd_dim_id,rcd)
      plr_cos_id=ncvdef(nc_id,'plr_cos',NF_REAL,1,plr_dim_id,rcd)
      plr_dgr_id=ncvdef(nc_id,'plr_dgr',NF_REAL,1,plr_dim_id,rcd)
      plr_id=ncvdef(nc_id,'plr',NF_REAL,1,plr_dim_id,rcd)
      rfl_bb_SAS_id=ncvdef(nc_id,'rfl_bb_SAS',NF_REAL,0,0,rcd)
      rfl_bb_sfc_id=ncvdef(nc_id,'rfl_bb_sfc',NF_REAL,0,0,rcd)
      rfl_nst_SAS_id=ncvdef(nc_id,'rfl_nst_SAS',NF_REAL,0,0,rcd)
      rfl_nst_sfc_id=ncvdef(nc_id,'rfl_nst_sfc',NF_REAL,0,0,rcd)
      rfl_spc_SAS_id=ncvdef(nc_id,'rfl_spc_SAS',NF_REAL,1,bnd_dim_id,rcd)
      t_id=ncvdef(nc_id,'t',NF_REAL,1,lev_dim_id,rcd)
      t_int_id=ncvdef(nc_id,'t_int',NF_REAL,1,levp_dim_id,rcd)
      tau_id=ncvdef(nc_id,'tau',NF_REAL,1,tau_dim_id,rcd)
      tau_prs_id=ncvdef(nc_id,'tau_prs',NF_REAL,1,tau_dim_id,rcd)
      trn_bb_atm_id=ncvdef(nc_id,'trn_bb_atm',NF_REAL,0,0,rcd)
      trn_nst_atm_id=ncvdef(nc_id,'trn_nst_atm',NF_REAL,0,0,rcd)
      trn_spc_atm_CO2_id=ncvdef(nc_id,'trn_spc_atm_CO2',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_H2OH2O_id=ncvdef(nc_id,'trn_spc_atm_H2OH2O',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_H2O_id=ncvdef(nc_id,'trn_spc_atm_H2O',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_NO2_id=ncvdef(nc_id,'trn_spc_atm_NO2',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_O2_id=ncvdef(nc_id,'trn_spc_atm_O2',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_O3_id=ncvdef(nc_id,'trn_spc_atm_O3',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_O2O2_id=ncvdef(nc_id,'trn_spc_atm_O2O2',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_OH_id=ncvdef(nc_id,'trn_spc_atm_OH',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_Ray_id=ncvdef(nc_id,'trn_spc_atm_Ray',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_aer_id=ncvdef(nc_id,'trn_spc_atm_aer',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_bga_id=ncvdef(nc_id,'trn_spc_atm_bga',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_ice_id=ncvdef(nc_id,'trn_spc_atm_ice',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_liq_id=ncvdef(nc_id,'trn_spc_atm_liq',NF_REAL,1,bnd_dim_id,rcd)
      trn_spc_atm_ttl_id=ncvdef(nc_id,'trn_spc_atm_ttl',NF_REAL,1,bnd_dim_id,rcd)
      wvl_ctr_id=ncvdef(nc_id,'wvl_ctr',NF_REAL,1,bnd_dim_id,rcd)
      wvl_max_id=ncvdef(nc_id,'wvl_max',NF_REAL,1,bnd_dim_id,rcd)
      wvl_min_id=ncvdef(nc_id,'wvl_min',NF_REAL,1,bnd_dim_id,rcd)
      wvl_obs_aer_id=ncvdef(nc_id,'wvl_obs_aer',NF_REAL,0,0,rcd)
      wvl_obs_bga_id=ncvdef(nc_id,'wvl_obs_bga',NF_REAL,0,0,rcd)
      wvl_sz_id=ncvdef(nc_id,'wvl_sz',NF_REAL,1,bnd_dim_id,rcd)
      wvn_ctr_id=ncvdef(nc_id,'wvn_ctr',NF_REAL,1,bnd_dim_id,rcd)
      wvn_max_id=ncvdef(nc_id,'wvn_max',NF_REAL,1,bnd_dim_id,rcd)
      wvn_min_id=ncvdef(nc_id,'wvn_min',NF_REAL,1,bnd_dim_id,rcd)
      wvn_sz_id=ncvdef(nc_id,'wvn_sz',NF_REAL,1,bnd_dim_id,rcd)
      z_cld_btm_id=ncvdef(nc_id,'z_cld_btm',NF_REAL,0,0,rcd)
      z_cld_thick_id=ncvdef(nc_id,'z_cld_thick',NF_REAL,0,0,rcd)
      z_id=ncvdef(nc_id,'z',NF_REAL,1,lev_dim_id,rcd)
      z_int_id=ncvdef(nc_id,'z_int',NF_REAL,1,levp_dim_id,rcd)

c     Add global attibutes
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'aer_sng',len(aer_sng),aer_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'azi_sng',len(azi_sng),azi_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'bga_sng',len(bga_sng),bga_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'creation_date',len(lcl_date_time),lcl_date_time)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'opt_dep_sng',len(opt_dep_sng),opt_dep_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'plr_sng',len(plr_sng),plr_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'prf_sng',len(prf_sng),prf_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'rcs_Id',len(rcs_Id),rcs_Id)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'str_sng',len(str_sng),str_sng)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_CO2',len(stt_CO2),stt_CO2)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_H2OH2O',len(stt_H2OH2O),stt_H2OH2O)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_H2O',len(stt_H2O),stt_H2O)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_Herzberg',len(stt_Herzberg),stt_Herzberg)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_NO2',len(stt_NO2),stt_NO2)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_O2',len(stt_O2),stt_O2)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_O3',len(stt_O3),stt_O3)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_O2O2',len(stt_O2O2),stt_O2O2)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_OH',len(stt_OH),stt_OH)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_Planck',len(stt_Planck),stt_Planck)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_Rayleigh',len(stt_Rayleigh),stt_Rayleigh)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_aer',len(stt_aer),stt_aer)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_bga',len(stt_bga),stt_bga)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_ice',len(stt_ice),stt_ice)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_liq',len(stt_liq),stt_liq)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_slr',len(stt_slr),stt_slr)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_flt_nst',len(stt_flt_nst),stt_flt_nst)
      rcd=nf_put_att_text(nc_id,NF_GLOBAL,'stt_top_lvl',len(stt_top_lvl),stt_top_lvl)

c     Add english text descriptions
      rcd=nf_put_att_text(nc_id,abs_bb_SAS_id,'long_name',50,'Broadband absorptance of surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,abs_bb_atm_id,'long_name',32,'Broadband absorptance of surface')
      rcd=nf_put_att_text(nc_id,abs_bb_sfc_id,'long_name',33,'Broadband absorptance of atmosphere')
      rcd=nf_put_att_text(nc_id,abs_nst_SAS_id,'long_name',45,'FSBR absorptance of surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,abs_nst_atm_id,'long_name',27,'FSBR absorptance of surface')
      rcd=nf_put_att_text(nc_id,abs_nst_sfc_id,'long_name',28,'FSBR absorptance of atmosphere')
      rcd=nf_put_att_text(nc_id,abs_spc_SAS_id,'long_name',49,'Spectral absorptance of surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,abs_spc_atm_id,'long_name',31,'Spectral absorptance of surface')
      rcd=nf_put_att_text(nc_id,abs_spc_sfc_id,'long_name',34,'Spectral absorptance of atmosphere')
      rcd=nf_put_att_text(nc_id,alb_sfc_id,'long_name',35,'Specified Lambertian surface albedo')
      rcd=nf_put_att_text(nc_id,azi_dgr_id,'long_name',25,'Azimuthal angle (degrees)')
      rcd=nf_put_att_text(nc_id,azi_id,'long_name',25,'Azimuthal angle (radians)')
      rcd=nf_put_att_text(nc_id,bnd_id,'long_name',19,'Midpoint wavelength')
      rcd=nf_put_att_text(nc_id,slr_zen_ngl_cos_id,'long_name',25,'Cosine solar zenith angle')
      rcd=nf_put_att_text(nc_id,flx_abs_atm_rdr_id,'long_name',49,'Flux absorbed in atmosphere at longer wavelengths')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_atm_id,'long_name',50,'Broadband flux absorbed by atmospheric column only')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_id,'long_name',32,'Broadband flux absorbed by layer')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_sfc_id,'long_name',39,'Broadband flux absorbed by surface only')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_ttl_id,'long_name',52,'Broadband flux absorbed by surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_TOA_id,'long_name',49,'Broadband incoming flux at TOA (total insolation)')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_dff_id,'long_name',34,'Diffuse downwelling broadband flux')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_drc_id,'long_name',33,'Direct downwelling broadband flux')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_id,'long_name',51,'Total downwelling broadband flux (direct + diffuse)')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_sfc_id,'long_name',37,'Broadband downwelling flux at surface')
      rcd=nf_put_att_text(nc_id,flx_bb_net_id,'long_name',44,'Net broadband flux (downwelling - upwelling)')
      rcd=nf_put_att_text(nc_id,flx_bb_up_id,'long_name',24,'Upwelling broadband flux')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_atm_id,'long_name',45,'FSBR flux absorbed by atmospheric column only')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_id,'long_name',27,'FSBR flux absorbed by layer')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_sfc_id,'long_name',34,'FSBR flux absorbed by surface only')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_ttl_id,'long_name',47,'FSBR flux absorbed by surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_TOA_id,'long_name',44,'FSBR incoming flux at TOA (total insolation)')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_id,'long_name',46,'Total downwelling FSBR flux (direct + diffuse)')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_sfc_id,'long_name',32,'FSBR downwelling flux at surface')
      rcd=nf_put_att_text(nc_id,flx_nst_net_id,'long_name',39,'Net FSBR flux (downwelling - upwelling)')
      rcd=nf_put_att_text(nc_id,flx_nst_up_id,'long_name',19,'Upwelling FSBR flux')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_SAS_id,'long_name',51,'Spectral flux absorbed by surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_atm_id,'long_name',49,'Spectral flux absorbed by atmospheric column only')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_id,'long_name',31,'Spectral flux absorbed by layer')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_sfc_id,'long_name',38,'Spectral flux absorbed by surface only')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_TOA_id,'long_name',32,'Spectral solar insolation at TOA')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_dff_id,'long_name',33,'Spectral diffuse downwelling flux')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_drc_id,'long_name',32,'Spectral direct downwelling flux')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_id,'long_name',25,'Spectral downwelling flux')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_sfc_id,'long_name',36,'Spectral solar insolation at surface')
      rcd=nf_put_att_text(nc_id,flx_spc_up_id,'long_name',23,'Spectral upwelling flux')
      rcd=nf_put_att_text(nc_id,frc_ice_ttl_id,'long_name',42,'Fraction of column condensate that is ice')
      rcd=nf_put_att_text(nc_id,flx_slr_frc_id,'long_name',22,'Fraction of solar flux')
      rcd=nf_put_att_text(nc_id,htg_rate_bb_id,'long_name',22,'Broadband heating rate')
      rcd=nf_put_att_text(nc_id,lat_dgr_id,'long_name',18,'Latitude (degrees)')
      rcd=nf_put_att_text(nc_id,lcl_time_hr_id,'long_name',14,'Local day hour')
      rcd=nf_put_att_text(nc_id,lcl_yr_day_id,'long_name',25,'Day of year in local time')
      rcd=nf_put_att_text(nc_id,lev_id,'long_name',14,'Layer pressure')
      rcd=nf_put_att_text(nc_id,levp_id,'long_name',18,'Interface pressure')
      rcd=nf_put_att_text(nc_id,mpc_CWP_id,'long_name',34,'Total column Condensed Water Path')
      rcd=nf_put_att_text(nc_id,ntn_bb_aa_id,'long_name',40,'Broadband azimuthally averaged intensity')
      rcd=nf_put_att_text(nc_id,ntn_bb_mean_id,'long_name',24,'Broadband mean intensity')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_ndr_id,'long_name',37,'Spectral intensity of nadir radiation')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_ndr_sfc_id,'long_name',48,'Spectral intensity of nadir radiation at surface')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_sfc_id,'long_name',42,'Spectral intensity of radiation at surface')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_zen_id,'long_name',38,'Spectral intensity of zenith radiation')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_zen_sfc_id,'long_name',49,'Spectral intensity of zenith radiation at surface')
      rcd=nf_put_att_text(nc_id,ntn_spc_chn_id,'long_name',42,'Full spectral intensity of particular band')
      rcd=nf_put_att_text(nc_id,ntn_spc_mean_id,'long_name',23,'Spectral mean intensity')
      rcd=nf_put_att_text(nc_id,odxc_obs_aer_id,'long_name',39,'Column aerosol extinction optical depth')
      rcd=nf_put_att_text(nc_id,odxc_obs_bga_id,'long_name',50,'Column background aerosol extinction optical depth')
      rcd=nf_put_att_text(nc_id,odxl_obs_aer_id,'long_name',38,'Layer aerosol extinction optical depth')
      rcd=nf_put_att_text(nc_id,odxl_obs_bga_id,'long_name',49,'Layer background aerosol extinction optical depth')
      rcd=nf_put_att_text(nc_id,odxc_spc_CO2_id,'long_name',28,'CO2 optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_H2OH2O_id,'long_name',34,'H2O dimer optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_H2O_id,'long_name',28,'H2O optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_NO2_id,'long_name',28,'NO2 optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_O2_id,'long_name',27,'O2 optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_O3_id,'long_name',27,'O3 optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_O2O2_id,'long_name',27,'O2O2 optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_OH_id,'long_name',27,'OH optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_Ray_id,'long_name',44,'Rayleigh scattering optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_aer_id,'long_name',48,'Aerosol scattering optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_bga_id,'long_name',59,'Background aerosol scattering optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_ice_id,'long_name',45,'Ice water scattering optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_liq_id,'long_name',48,'Liquid water scattering optical depth to surface')
      rcd=nf_put_att_text(nc_id,odxc_spc_ttl_id,'long_name',41,'Total extinction optical depth to surface')
      rcd=nf_put_att_text(nc_id,plr_cos_id,'long_name',28,'Cosine polar angle (degrees)')
      rcd=nf_put_att_text(nc_id,plr_dgr_id,'long_name',21,'Polar angle (degrees)')
      rcd=nf_put_att_text(nc_id,plr_id,'long_name',21,'Polar angle (radians)')
      rcd=nf_put_att_text(nc_id,rfl_bb_SAS_id,'long_name',53,'Broadband albedo of entire surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,rfl_bb_sfc_id,'long_name',27,'Broadband albedo of surface')
      rcd=nf_put_att_text(nc_id,rfl_nst_SAS_id,'long_name',48,'FSBR albedo of entire surface-atmosphere system')
      rcd=nf_put_att_text(nc_id,rfl_nst_sfc_id,'long_name',22,'FSBR albedo of surface')
      rcd=nf_put_att_text(nc_id,rfl_spc_SAS_id,'long_name',35,'Spectral planetary flux reflectance')
      rcd=nf_put_att_text(nc_id,t_id,'long_name',17,'Layer Temperature')
      rcd=nf_put_att_text(nc_id,t_int_id,'long_name',21,'Interface temperature')
      rcd=nf_put_att_text(nc_id,tau_id,'long_name',25,'Optical level (optical depth)')
      rcd=nf_put_att_text(nc_id,tau_prs_id,'long_name',25,'Optical level (pressure)')
      rcd=nf_put_att_text(nc_id,trn_bb_atm_id,'long_name',44,'Broadband transmission of atmospheric column')
      rcd=nf_put_att_text(nc_id,trn_nst_atm_id,'long_name',39,'FSBR transmission of atmospheric column')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_CO2_id,'long_name',41,'Column transmission due to CO2 absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_H2OH2O_id,'long_name',47,'Column transmission due to H2O dimer absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_H2O_id,'long_name',41,'Column transmission due to H2O absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_NO2_id,'long_name',41,'Column transmission due to NO2 absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O2_id,'long_name',40,'Column transmission due to O2 absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O3_id,'long_name',40,'Column transmission due to O3 absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O2O2_id,'long_name',45,'Column transmission due to O2-O2 absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_OH_id,'long_name',40,'Column transmission due to OH absorption')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_Ray_id,'long_name',46,'Column transmission due to Rayleigh scattering')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_aer_id,'long_name',44,'Column transmission due to aerosol extinction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_bga_id,'long_name',55,'Column transmission due to background aerosol extinction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_ice_id,'long_name',41,'Column transmission due to ice extinction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_liq_id,'long_name',44,'Column transmission due to liquid extinction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_ttl_id,'long_name',43,'Spectral flux transmission of entire column')
      rcd=nf_put_att_text(nc_id,wvl_ctr_id,'long_name',27,'Midpoint wavelength in band')
      rcd=nf_put_att_text(nc_id,wvl_max_id,'long_name',26,'Maximum wavelength in band')
      rcd=nf_put_att_text(nc_id,wvl_min_id,'long_name',26,'Minimum wavelength in band')
      rcd=nf_put_att_text(nc_id,wvl_obs_aer_id,'long_name',49,'Wavelength of aerosol optical depth specification')
      rcd=nf_put_att_text(nc_id,wvl_obs_bga_id,'long_name',60,'Wavelength of background aerosol optical depth specification')
      rcd=nf_put_att_text(nc_id,wvl_sz_id,'long_name',13,'Width of band')
      rcd=nf_put_att_text(nc_id,wvn_ctr_id,'long_name',26,'Maximum wavenumber in band')
      rcd=nf_put_att_text(nc_id,wvn_max_id,'long_name',26,'Maximum wavenumber in band')
      rcd=nf_put_att_text(nc_id,wvn_min_id,'long_name',26,'Minimum wavenumber in band')
      rcd=nf_put_att_text(nc_id,wvn_sz_id,'long_name',24,'Bandwidth in wavenumbers')
      rcd=nf_put_att_text(nc_id,z_cld_btm_id,'long_name',46,'Highest interface beneath all clouds in column')
      rcd=nf_put_att_text(nc_id,z_cld_thick_id,'long_name',41,'Thickness of region containing all clouds')
      rcd=nf_put_att_text(nc_id,z_id,'long_name',8,'Altitude')
      rcd=nf_put_att_text(nc_id,z_int_id,'long_name',18,'Interface altitude')

c     Add units
      rcd=nf_put_att_text(nc_id,abs_bb_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_bb_atm_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_bb_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_nst_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_nst_atm_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_nst_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_spc_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_spc_atm_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,abs_spc_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,alb_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,azi_dgr_id,'units',6,'degree')
      rcd=nf_put_att_text(nc_id,azi_id,'units',6,'radian')
      rcd=nf_put_att_text(nc_id,bnd_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,slr_zen_ngl_cos_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,flx_abs_atm_rdr_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_atm_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_sfc_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_abs_ttl_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_TOA_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_dff_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_drc_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_dwn_sfc_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_net_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_bb_up_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_atm_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_sfc_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_abs_ttl_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_TOA_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_dwn_sfc_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_net_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_nst_up_id,'units',12,'watt meter-2')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_SAS_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_atm_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_abs_sfc_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_TOA_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_dff_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_drc_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_dwn_sfc_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,flx_spc_up_id,'units',20,'watt meter-2 meter-1')
      rcd=nf_put_att_text(nc_id,frc_ice_ttl_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,flx_slr_frc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,htg_rate_bb_id,'units',15,'kelvin second-1')
      rcd=nf_put_att_text(nc_id,lat_dgr_id,'units',6,'degree')
      rcd=nf_put_att_text(nc_id,lcl_time_hr_id,'units',4,'hour')
      rcd=nf_put_att_text(nc_id,lcl_yr_day_id,'units',3,'day')
      rcd=nf_put_att_text(nc_id,lev_id,'units',6,'pascal')
      rcd=nf_put_att_text(nc_id,levp_id,'units',6,'pascal')
      rcd=nf_put_att_text(nc_id,mpc_CWP_id,'units',16,'kilogram meter-2')
      rcd=nf_put_att_text(nc_id,ntn_bb_aa_id,'units',25,'watt meter-2 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_bb_mean_id,'units',25,'watt meter-2 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_ndr_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_ndr_sfc_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_sfc_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_zen_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_aa_zen_sfc_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_chn_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,ntn_spc_mean_id,'units',33,'watt meter-2 meter-1 sterradian-1')
      rcd=nf_put_att_text(nc_id,odxc_obs_aer_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_obs_bga_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxl_obs_aer_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxl_obs_bga_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_CO2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_H2OH2O_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_H2O_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_NO2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_O2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_O3_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_O2O2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_OH_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_Ray_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_aer_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_bga_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_ice_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_liq_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,odxc_spc_ttl_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,plr_cos_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,plr_dgr_id,'units',6,'degree')
      rcd=nf_put_att_text(nc_id,plr_id,'units',6,'radian')
      rcd=nf_put_att_text(nc_id,rfl_bb_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,rfl_bb_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,rfl_nst_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,rfl_nst_sfc_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,rfl_spc_SAS_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,t_id,'units',6,'kelvin')
      rcd=nf_put_att_text(nc_id,t_int_id,'units',6,'kelvin')
      rcd=nf_put_att_text(nc_id,tau_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,tau_prs_id,'units',6,'pascal')
      rcd=nf_put_att_text(nc_id,trn_bb_atm_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_nst_atm_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_CO2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_H2OH2O_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_H2O_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_NO2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O3_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_O2O2_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_OH_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_Ray_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_aer_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_bga_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_ice_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_liq_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,trn_spc_atm_ttl_id,'units',8,'fraction')
      rcd=nf_put_att_text(nc_id,wvl_ctr_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvl_max_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvl_min_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvl_obs_aer_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvl_obs_bga_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvl_sz_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,wvn_ctr_id,'units',12,'centimeter-1')
      rcd=nf_put_att_text(nc_id,wvn_max_id,'units',12,'centimeter-1')
      rcd=nf_put_att_text(nc_id,wvn_min_id,'units',12,'centimeter-1')
      rcd=nf_put_att_text(nc_id,wvn_sz_id,'units',12,'centimeter-1')
      rcd=nf_put_att_text(nc_id,z_cld_btm_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,z_cld_thick_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,z_id,'units',5,'meter')
      rcd=nf_put_att_text(nc_id,z_int_id,'units',5,'meter')

c     All dimensions, variables, and attributes have been defined
      rcd=nf_enddef(nc_id)

c     Write out data
      srt(1)=1
      srt(2)=1
      srt(3)=1
      srt(4)=1

      call ncvpt(nc_id,abs_spc_SAS_id,srt,nbr_bnd,abs_spc_SAS,rcd)
      call ncvpt(nc_id,abs_spc_atm_id,srt,nbr_bnd,abs_spc_atm,rcd)
      call ncvpt(nc_id,abs_spc_sfc_id,srt,nbr_bnd,abs_spc_sfc,rcd)
      call ncvpt(nc_id,azi_dgr_id,srt,nbr_azi,azi_dgr,rcd)
      call ncvpt(nc_id,azi_id,srt,nbr_azi,azi,rcd)
      call ncvpt(nc_id,bnd_id,srt,nbr_bnd,bnd,rcd)
      call ncvpt(nc_id,flx_abs_atm_rdr_id,srt,nbr_bnd,flx_abs_atm_rdr,rcd)
      call ncvpt(nc_id,flx_bb_abs_id,srt,nbr_lev,flx_bb_abs,rcd)
      call ncvpt(nc_id,flx_bb_dwn_dff_id,srt,nbr_levp,flx_bb_dwn_dff,rcd)
      call ncvpt(nc_id,flx_bb_dwn_drc_id,srt,nbr_levp,flx_bb_dwn_drc,rcd)
      call ncvpt(nc_id,flx_bb_dwn_id,srt,nbr_levp,flx_bb_dwn,rcd)
      call ncvpt(nc_id,flx_bb_net_id,srt,nbr_levp,flx_bb_net,rcd)
      call ncvpt(nc_id,flx_bb_up_id,srt,nbr_levp,flx_bb_up,rcd)
      call ncvpt(nc_id,flx_nst_abs_id,srt,nbr_lev,flx_nst_abs,rcd)
      call ncvpt(nc_id,flx_nst_dwn_id,srt,nbr_levp,flx_nst_dwn,rcd)
      call ncvpt(nc_id,flx_nst_net_id,srt,nbr_levp,flx_nst_net,rcd)
      call ncvpt(nc_id,flx_nst_up_id,srt,nbr_levp,flx_nst_up,rcd)
      call ncvpt(nc_id,flx_spc_abs_SAS_id,srt,nbr_bnd,flx_spc_abs_SAS,rcd)
      call ncvpt(nc_id,flx_spc_abs_atm_id,srt,nbr_bnd,flx_spc_abs_atm,rcd)
      call ncvpt(nc_id,flx_spc_abs_id,srt,cnt_bnd_lev,flx_spc_abs,rcd)
      call ncvpt(nc_id,flx_spc_abs_sfc_id,srt,nbr_bnd,flx_spc_abs_sfc,rcd)
      call ncvpt(nc_id,flx_spc_dwn_TOA_id,srt,nbr_bnd,flx_spc_dwn_TOA,rcd)
      call ncvpt(nc_id,flx_spc_dwn_dff_id,srt,cnt_bnd_levp,flx_spc_dwn_dff,rcd)
      call ncvpt(nc_id,flx_spc_dwn_drc_id,srt,cnt_bnd_levp,flx_spc_dwn_drc,rcd)
      call ncvpt(nc_id,flx_spc_dwn_id,srt,cnt_bnd_levp,flx_spc_dwn,rcd)
      call ncvpt(nc_id,flx_spc_dwn_sfc_id,srt,nbr_bnd,flx_spc_dwn_sfc,rcd)
      call ncvpt(nc_id,flx_spc_up_id,srt,cnt_bnd_levp,flx_spc_up,rcd)
      call ncvpt(nc_id,flx_slr_frc_id,srt,nbr_bnd,flx_slr_frc,rcd)
      call ncvpt(nc_id,htg_rate_bb_id,srt,nbr_lev,htg_rate_bb,rcd)
      call ncvpt(nc_id,lev_id,srt,nbr_lev,lev,rcd)
      call ncvpt(nc_id,levp_id,srt,nbr_levp,levp,rcd)
      call ncvpt(nc_id,ntn_bb_aa_id,srt,cnt_plr_levp,ntn_bb_aa,rcd)
      call ncvpt(nc_id,ntn_bb_mean_id,srt,nbr_lev,ntn_bb_mean,rcd)
      call ncvpt(nc_id,ntn_spc_aa_ndr_id,srt,cnt_bnd_levp,ntn_spc_aa_ndr,rcd)
      call ncvpt(nc_id,ntn_spc_aa_ndr_sfc_id,srt,nbr_bnd,ntn_spc_aa_ndr_sfc,rcd)
      call ncvpt(nc_id,ntn_spc_aa_sfc_id,srt,cnt_plr_bnd,ntn_spc_aa_sfc,rcd)
      call ncvpt(nc_id,ntn_spc_aa_zen_id,srt,cnt_bnd_levp,ntn_spc_aa_zen,rcd)
      call ncvpt(nc_id,ntn_spc_aa_zen_sfc_id,srt,nbr_bnd,ntn_spc_aa_zen_sfc,rcd)
      call ncvpt(nc_id,ntn_spc_chn_id,srt,cnt_azi_plr_levp,ntn_spc_chn,rcd)
      call ncvpt(nc_id,ntn_spc_mean_id,srt,cnt_bnd_lev,ntn_spc_mean,rcd)
      call ncvpt(nc_id,odxc_spc_CO2_id,srt,nbr_bnd,odxc_spc_CO2,rcd)
      call ncvpt(nc_id,odxc_spc_H2OH2O_id,srt,nbr_bnd,odxc_spc_H2OH2O,rcd)
      call ncvpt(nc_id,odxc_spc_H2O_id,srt,nbr_bnd,odxc_spc_H2O,rcd)
      call ncvpt(nc_id,odxc_spc_NO2_id,srt,nbr_bnd,odxc_spc_NO2,rcd)
      call ncvpt(nc_id,odxc_spc_O2_id,srt,nbr_bnd,odxc_spc_O2,rcd)
      call ncvpt(nc_id,odxc_spc_O3_id,srt,nbr_bnd,odxc_spc_O3,rcd)
      call ncvpt(nc_id,odxc_spc_O2O2_id,srt,nbr_bnd,odxc_spc_O2O2,rcd)
      call ncvpt(nc_id,odxc_spc_OH_id,srt,nbr_bnd,odxc_spc_OH,rcd)
      call ncvpt(nc_id,odxc_spc_Ray_id,srt,nbr_bnd,odxc_spc_Ray,rcd)
      call ncvpt(nc_id,odxc_spc_aer_id,srt,nbr_bnd,odxc_spc_aer,rcd)
      call ncvpt(nc_id,odxc_spc_bga_id,srt,nbr_bnd,odxc_spc_bga,rcd)
      call ncvpt(nc_id,odxc_spc_ice_id,srt,nbr_bnd,odxc_spc_ice,rcd)
      call ncvpt(nc_id,odxc_spc_liq_id,srt,nbr_bnd,odxc_spc_liq,rcd)
      call ncvpt(nc_id,odxc_spc_ttl_id,srt,nbr_bnd,odxc_spc_ttl,rcd)
      call ncvpt(nc_id,odxl_obs_aer_id,srt,nbr_lev,odxl_obs_aer,rcd)
      call ncvpt(nc_id,odxl_obs_bga_id,srt,nbr_lev,odxl_obs_bga,rcd)
      call ncvpt(nc_id,plr_cos_id,srt,nbr_plr,plr_cos,rcd)
      call ncvpt(nc_id,plr_dgr_id,srt,nbr_plr,plr_dgr,rcd)
      call ncvpt(nc_id,plr_id,srt,nbr_plr,plr,rcd)
      call ncvpt(nc_id,rfl_spc_SAS_id,srt,nbr_bnd,rfl_spc_SAS,rcd)
      call ncvpt(nc_id,t_id,srt,nbr_lev,t,rcd)
      call ncvpt(nc_id,t_int_id,srt,nbr_levp,t_int,rcd)
      call ncvpt(nc_id,tau_id,srt,nbr_tau,tau,rcd)
      call ncvpt(nc_id,tau_prs_id,srt,nbr_tau,tau_prs,rcd)
      call ncvpt(nc_id,trn_spc_atm_CO2_id,srt,nbr_bnd,trn_spc_atm_CO2,rcd)
      call ncvpt(nc_id,trn_spc_atm_H2OH2O_id,srt,nbr_bnd,trn_spc_atm_H2OH2O,rcd)
      call ncvpt(nc_id,trn_spc_atm_H2O_id,srt,nbr_bnd,trn_spc_atm_H2O,rcd)
      call ncvpt(nc_id,trn_spc_atm_NO2_id,srt,nbr_bnd,trn_spc_atm_NO2,rcd)
      call ncvpt(nc_id,trn_spc_atm_O2_id,srt,nbr_bnd,trn_spc_atm_O2,rcd)
      call ncvpt(nc_id,trn_spc_atm_O3_id,srt,nbr_bnd,trn_spc_atm_O3,rcd)
      call ncvpt(nc_id,trn_spc_atm_O2O2_id,srt,nbr_bnd,trn_spc_atm_O2O2,rcd)
      call ncvpt(nc_id,trn_spc_atm_OH_id,srt,nbr_bnd,trn_spc_atm_OH,rcd)
      call ncvpt(nc_id,trn_spc_atm_Ray_id,srt,nbr_bnd,trn_spc_atm_Ray,rcd)
      call ncvpt(nc_id,trn_spc_atm_aer_id,srt,nbr_bnd,trn_spc_atm_aer,rcd)
      call ncvpt(nc_id,trn_spc_atm_bga_id,srt,nbr_bnd,trn_spc_atm_bga,rcd)
      call ncvpt(nc_id,trn_spc_atm_ice_id,srt,nbr_bnd,trn_spc_atm_ice,rcd)
      call ncvpt(nc_id,trn_spc_atm_liq_id,srt,nbr_bnd,trn_spc_atm_liq,rcd)
      call ncvpt(nc_id,trn_spc_atm_ttl_id,srt,nbr_bnd,trn_spc_atm_ttl,rcd)
      call ncvpt(nc_id,wvl_ctr_id,srt,nbr_bnd,wvl_ctr,rcd)
      call ncvpt(nc_id,wvl_max_id,srt,nbr_bnd,wvl_max,rcd)
      call ncvpt(nc_id,wvl_min_id,srt,nbr_bnd,wvl_min,rcd)
      call ncvpt(nc_id,wvl_sz_id,srt,nbr_bnd,wvl_sz,rcd)
      call ncvpt(nc_id,wvn_ctr_id,srt,nbr_bnd,wvn_ctr,rcd)
      call ncvpt(nc_id,wvn_max_id,srt,nbr_bnd,wvn_max,rcd)
      call ncvpt(nc_id,wvn_min_id,srt,nbr_bnd,wvn_min,rcd)
      call ncvpt(nc_id,wvn_sz_id,srt,nbr_bnd,wvn_sz,rcd)
      call ncvpt(nc_id,z_id,srt,nbr_lev,z,rcd)
      call ncvpt(nc_id,z_int_id,srt,nbr_levp,z_int,rcd)
      call ncvpt1(nc_id,abs_bb_SAS_id,0,abs_bb_SAS,rcd)
      call ncvpt1(nc_id,abs_bb_atm_id,0,abs_bb_atm,rcd)
      call ncvpt1(nc_id,abs_bb_sfc_id,0,abs_bb_sfc,rcd)
      call ncvpt1(nc_id,abs_nst_SAS_id,0,abs_nst_SAS,rcd)
      call ncvpt1(nc_id,abs_nst_atm_id,0,abs_nst_atm,rcd)
      call ncvpt1(nc_id,abs_nst_sfc_id,0,abs_nst_sfc,rcd)
      call ncvpt1(nc_id,alb_sfc_id,0,alb_sfc,rcd)
      call ncvpt1(nc_id,slr_zen_ngl_cos_id,srt,slr_zen_ngl_cos,rcd)
      call ncvpt1(nc_id,flx_bb_abs_atm_id,0,flx_bb_abs_atm,rcd)
      call ncvpt1(nc_id,flx_bb_abs_sfc_id,0,flx_bb_abs_sfc,rcd)
      call ncvpt1(nc_id,flx_bb_abs_ttl_id,0,flx_bb_abs_ttl,rcd)
      call ncvpt1(nc_id,flx_bb_dwn_TOA_id,0,flx_bb_dwn_TOA,rcd)
      call ncvpt1(nc_id,flx_bb_dwn_sfc_id,0,flx_bb_dwn_sfc,rcd)
      call ncvpt1(nc_id,flx_nst_abs_atm_id,0,flx_nst_abs_atm,rcd)
      call ncvpt1(nc_id,flx_nst_abs_sfc_id,0,flx_nst_abs_sfc,rcd)
      call ncvpt1(nc_id,flx_nst_abs_ttl_id,0,flx_nst_abs_ttl,rcd)
      call ncvpt1(nc_id,flx_nst_dwn_TOA_id,0,flx_nst_dwn_TOA,rcd)
      call ncvpt1(nc_id,flx_nst_dwn_sfc_id,0,flx_nst_dwn_sfc,rcd)
      call ncvpt1(nc_id,frc_ice_ttl_id,0,frc_ice_ttl,rcd)
      call ncvpt1(nc_id,lat_dgr_id,0,lat_dgr,rcd)
      call ncvpt1(nc_id,lcl_time_hr_id,0,lcl_time_hr,rcd)
      call ncvpt1(nc_id,lcl_yr_day_id,0,lcl_yr_day,rcd)
      call ncvpt1(nc_id,mpc_CWP_id,0,mpc_CWP,rcd)
      call ncvpt1(nc_id,odxc_obs_aer_id,0,odxc_obs_aer,rcd)
      call ncvpt1(nc_id,odxc_obs_bga_id,0,odxc_obs_bga,rcd)
      call ncvpt1(nc_id,rfl_bb_SAS_id,0,rfl_bb_SAS,rcd)
      call ncvpt1(nc_id,rfl_bb_sfc_id,0,rfl_bb_sfc,rcd)
      call ncvpt1(nc_id,rfl_nst_SAS_id,0,rfl_nst_SAS,rcd)
      call ncvpt1(nc_id,rfl_nst_sfc_id,0,rfl_nst_sfc,rcd)
      call ncvpt1(nc_id,trn_bb_atm_id,0,trn_bb_atm,rcd)
      call ncvpt1(nc_id,trn_nst_atm_id,0,trn_nst_atm,rcd)
      call ncvpt1(nc_id,wvl_obs_aer_id,srt,wvl_obs_aer,rcd)
      call ncvpt1(nc_id,wvl_obs_bga_id,srt,wvl_obs_bga,rcd)
      call ncvpt1(nc_id,z_cld_btm_id,0,z_cld_btm,rcd)
      call ncvpt1(nc_id,z_cld_thick_id,0,z_cld_thick,rcd)

      if (sv_ntn) then

c     $     ntn_spc_aa_id,
c     $     ntn_spc_aa(max_nbr_plr,bnd_nbr_max,max_nbr_levp),
c         call ncredf(nc_id,rcd)
c         ntn_spc_aa_id=ncvdef(nc_id,'ntn_spc_aa',NF_REAL,3,dim_plr_bnd_levp,rcd)
c         call ncaptc(nc_id,ntn_spc_aa_id,'long_name',NF_CHAR,31,'Spectral intensity of radiation',rcd)
c         call ncaptc(nc_id,ntn_spc_aa_id,'units',NF_CHAR,33,'watt meter-2 meter-1 sterradian-1',rcd)
c         rcd=nf_enddef(nc_id)
c         call ncvpt(nc_id,ntn_spc_aa_id,srt,cnt_plr_bnd_levp,ntn_spc_aa,rcd)
      endif                     ! end if saving full intensity arrays

      rcd=nf_close(nc_id)
      if (rcd.ne.NF_NOERR) call nc_err_exit(rcd,'ERROR '//fl_out)
      write (6,'(a16,1x,a40)') 'Wrote results to',fl_out

 1000 continue

#include <time_end.com> /* Print date, elapsed time */
#include <ieee_end.com> /* Print and reset IEEE errors, if any */

      call exit(exit_status)
      end                       ! end lbl()

