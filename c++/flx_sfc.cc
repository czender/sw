// $Id$ 

// Purpose: Surface flux physics utilities for C++ programs

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <flx_sfc.hh> // Surface flux physics

// Namespaces
using phc::tpt_frz_pnt; // (273.15) [K] Kelvin--Celsius scale offset, see, e.g., Bol80
using phc::ltn_heat_fsn_H2O_std; // (0.3336e06) [J kg-1] Latent heat of fusion of H2O at 0 C, standard CCM:lsm/phyconi.F 
using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant
using phc::eps_H2O_rcp_m1; // (0.60777) [frc] Constant for virtual temperature
using phc::ltn_heat_vpr_H2O_std; // (2.5104e06) [J kg-1] Latent heat of vaporization of H2O, standard CCM:lsm/phyconi.F
using phc::ltn_heat_sbl_H2O_std; // (2.8440e06) [J kg-1] Latent heat of sublimation of H2O, standard CCM:lsm/phyconi.F
using phc::spc_heat_dry_air; // (1004.697) [J kg-1 K-1] IrG81 p. 25
using phc::cst_Stefan_Boltzmann; // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462
using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
using phc::cst_von_krm_rcp; // (2.5) [frc] Reciprocal of Von Karman's constant
using phc::ssh_H2O_sln_crc; // (0.98) [frc] Salinity correction to saturation specific humidity of H2O LaP81 p. 328 CCM:dom/flxoce()
using phc::cp_vpr_rcp_cp_dry_m1; // (0.83745) Constant for moist specific heat IrG81 p. 77

// Soil thermal conductivity
int // O [rcd] Return success code
cnd_trm_soi_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnd_trm_soi, // O [W m-1 K-1] Soil thermal conductivity
 prc_cmp *lvl_dlt_snw, // O [m] Soil layer thickness including snow
 const prc_cmp *mss_frc_cly, // I [frc] Mass fraction clay 
 const prc_cmp *mss_frc_snd, // I [frc] Mass fraction sand
 const prc_cmp *snw_hgt, // I [m] Geometric bulk thickness of snow
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
 const prc_cmp *vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
 const prc_cmp *vwc_sat, // I [m3 m-3] Saturated volumetric water content
 const prc_cmp *vwc_sfc) // I [m3 m-3] Volumetric water content
{
  // Purpose: Thermal properties of soil
  // Dependencies: <phys_cst.hh>, <lsm.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  using phc::cnd_trm_H2O_ice; // (2.2) [W m-1 K-1] Thermal conductivity of ice water CCM:lsm/phyconi
  using phc::cnd_trm_H2O_lqd; // (0.6) [W m-1 K-1] Thermal conductivity of liquid water CCM:lsm/phyconi
  using phc::cnd_trm_snw; // (0.34) [W m-1 K-1] Thermal conductivity of snow CCM:lsm/snoconi
  const prc_cmp lvl_dlt_sfc(0.1); // [m] Soil layer thickness, top layer
  const prc_cmp tpt_dlt(0.5); // [K] Temperature range of mixed phase soil
  std::valarray<prc_cmp> cnd_trm_soi_dry(lon_nbr); // [W m-1 K-1] Thermal conductivity of dry soil
  std::valarray<prc_cmp> cnd_trm_soi_frz(lon_nbr); // [W m-1 K-1] Soil thermal conductivity, frozen
  std::valarray<prc_cmp> cnd_trm_soi_sld(lon_nbr); // [W m-1 K-1] Thermal conductivity of soil solids
  std::valarray<prc_cmp> cnd_trm_soi_wrm(lon_nbr); // [W m-1 K-1] Soil thermal conductivity, unfrozen
  std::valarray<prc_cmp> ltn_heat_fsn_vlm(lon_nbr); // [J m-3] Volumetric latent heat of fusion
  std::valarray<prc_cmp> lvl_dlt(lon_nbr); // [m] Soil layer thickness, no snow
  prc_cmp snw_hgt_bnd; // [m] Bounded geometric bulk thickness of snow
  long lon_idx; // Counting index
  // Main Code
  const std::string sbr_nm("cnd_trm_soi_get");
  
  // Initialize arrays
  vec_set(&lvl_dlt[0],lon_nbr,lvl_dlt_sfc); // [m] Soil layer thickness
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    ltn_heat_fsn_vlm[lon_idx]=vwc_sfc[lon_idx]*ltn_heat_fsn_H2O_std*dns_H2O_lqd_std; // [J m-3] Volumetric latent heat of fusion
    cnd_trm_soi_sld[lon_idx]=(8.80*mss_frc_snd[lon_idx]+2.92*mss_frc_cly[lon_idx])/(mss_frc_snd[lon_idx]+mss_frc_cly[lon_idx]); // [W m-1 K-1] Thermal conductivity of soil solids CCM:lsm/lsmtci() Bon96 p. 77
    cnd_trm_soi_dry[lon_idx]=0.15; // [W m-1 K-1] Thermal conductivity of dry soil CCM:lsm/lsmtci()
    cnd_trm_soi_wrm[lon_idx]= // [W m-1 K-1] Soil thermal conductivity, unfrozen
      +cnd_trm_soi_dry[lon_idx]
      +(std::pow(cnd_trm_soi_sld[lon_idx],(PRC_CMP(1.0)-vwc_sat[lon_idx]))
	*std::pow(cnd_trm_H2O_lqd,vwc_sfc[lon_idx])-cnd_trm_soi_dry[lon_idx])
      *vwc_sfc[lon_idx]/vwc_sat[lon_idx];
    cnd_trm_soi_frz[lon_idx]= // [W m-1 K-1] Soil thermal conductivity, frozen
      +cnd_trm_soi_dry[lon_idx]
      +(std::pow(cnd_trm_soi_sld[lon_idx],(PRC_CMP(1.0)-vwc_sat[lon_idx]))
	*std::pow(cnd_trm_H2O_ice,vwc_sfc[lon_idx])-cnd_trm_soi_dry[lon_idx])
      *vwc_sfc[lon_idx]/vwc_sat[lon_idx];
    if(tpt_soi[lon_idx] < tpt_frz_pnt-tpt_dlt){
      cnd_trm_soi[lon_idx]=cnd_trm_soi_frz[lon_idx]; // [W m-1 K-1] Soil thermal conductivity
    } // endif
    if(tpt_soi[lon_idx] >= tpt_frz_pnt-tpt_dlt && tpt_soi[lon_idx] <= tpt_frz_pnt+tpt_dlt){
      cnd_trm_soi[lon_idx]= // [W m-1 K-1] Soil thermal conductivity
	+cnd_trm_soi_frz[lon_idx]
	+(cnd_trm_soi_frz[lon_idx]-cnd_trm_soi_wrm[lon_idx])
	*(tpt_soi[lon_idx]-tpt_frz_pnt+tpt_dlt)
	/(2.0*tpt_dlt);
    } // endif
    if(tpt_soi[lon_idx] > tpt_frz_pnt+tpt_dlt){
      cnd_trm_soi[lon_idx]=cnd_trm_soi_wrm[lon_idx]; // [W m-1 K-1] Soil thermal conductivity
    } // endif

    // Blend snow into first soil layer
    snw_hgt_bnd=min_cpv(snw_hgt[lon_idx],1.0); // [m] Bounded geometric bulk thickness of snow CCM:lsm/surphy()
    lvl_dlt_snw[lon_idx]=lvl_dlt[lon_idx]+snw_hgt_bnd; // O [m] Soil layer thickness including snow CCM:lsm/surphy()
    cnd_trm_soi[lon_idx]= // [W m-1 K-1] Soil thermal conductivity Bon96 p. 77
      cnd_trm_snw*cnd_trm_soi[lon_idx]*lvl_dlt_snw[lon_idx]
      /(cnd_trm_snw*lvl_dlt[lon_idx]+cnd_trm_soi[lon_idx]*snw_hgt_bnd);
  }  // end loop over lon

  if(vwc_dry[0]+vwc_opt[0]==0.0){;} // CEWU Compiler Error Warning Usage

  return rcd;
} // end cnd_trm_soi_get()

int // O [rcd] Return success code
flx_sfc_lnd
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *vgt, // I [flg] "Vegetated" flag
 const prc_cmp *cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
 const prc_cmp *flx_SW_net_gnd, // I [W m-2] Solar flux absorbed by ground
 const prc_cmp *flx_SW_net_vgt, // I [W m-2] Solar flux absorbed by vegetation
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *lvl_dlt_snw, // I [m] Soil layer thickness including snow
 const prc_cmp *msv_gnd, // I [frc] Bare ground emissivity
 const prc_cmp *msv_snw, // I [frc] Snow emissivity
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 const prc_cmp tm_dlt, // I [s] Timestep
 const long *pnt_typ_idx, // I [idx] Plant type index (1..14)
 const long *soi_typ, // I [idx] LSM soil type (1..5)
 prc_cmp *cff_xch_heat, // O [frc] Exchange coefficient for heat transfer
 prc_cmp *cff_xch_mmn, // O [frc] Exchange coefficient for momentum transfer
 prc_cmp *cff_xch_mmn_ntr, // O [frc] Neutral drag coefficient hgt_mdp to z0m+zpd 
 prc_cmp *cff_xch_vpr, // O [frc] Exchange coefficient for vapor transfer
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm_ttl, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *flx_sns_gnd_ttl, // O [W m-2] Sensible heat flux to soil
 prc_cmp *flx_snw_mlt, // O [W m-2] Snow melt heat flux
 prc_cmp *mno_lng, // I/O [m] Monin-Obukhov length
 prc_cmp *msv_sfc, // O [frc] Surface (bare bare ground+snow) emissivity
 prc_cmp *rss_aer_heat_sfc, // O [s m-1] Aerodynamic resistance to heat transfer
 prc_cmp *rss_aer_mmn_sfc, // O [s m-1] Aerodynamic resistance to momentum transfer
 prc_cmp *rss_aer_vpr_sfc, // O [s m-1] Aerodynamic resistance to vapor transfer
 prc_cmp *tpt_aer, // O [K] "Aerodynamic" temperature at z=zpd+rgh_mmn
 prc_cmp *tpt_ash, // O [K] "Surface" temperature at z=zpd+rgh_heat
 prc_cmp *tpt_ash_p2m, // O [K] "Screen" temperature at z=zpd+rgh_heat+2m
 prc_cmp *tpt_ffc, // O [K] Radiative effective temperature
 prc_cmp *tpt_gnd, // I/O [K] Ground temperature
 prc_cmp *tpt_vgt, // I/O [K] Vegetation temperature
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl) // O [kg m-1 s-2] Zonal wind stress
{
  /* Purpose: Given meteorology at reference level (i.e., GCM layer midpoint),
     compute the surface fluxes. 
     To get the fluxes, we must solve for the surface temperature
     The term "surface temperature" is generic and exactly what it describes 
     depends on the given surface type, i.e., bare ground, snow, ocean, ice, or vegetation.
  */
  /* Variables marked CEWI "Compiler Error Warning Initialization" are 
     initialized to 2.0e36 to prevent spurious 
     "warning: `float foo' might be used uninitialized in this function"
     compiler warnings when -Wunitialized is turned on */
  // Theory and algorithms: Bonan (1996) CCM:lsm/surtem()
  // Dependencies: <phys_cst.hh>,<tdy.hh>
  using lsm::hvt; // [m] Height at top of canopy
  using lsm::dleaf; // [m] Characteristic leaf dimension
 
  /* Notes on algorithm:
     Suffix mdp quantity evaluated at height hgt_mdp
     Suffix gnd quantity evaluated at ground
     Suffix vgt quantity evaluated in vegetation
     Suffix sfc quantity evaluated at surface (i.e., gnd or vgt, depending)
     Currently, we neglect the effects of vegetation */
  
  /* Named surface type indices idx_vgt and idx_gnd enumerate the two layers
     where the surface budget equation is solved iteratively, in the order
     that they are solved. 
     Over bare ground, first iterative loop computes ground to atmosphere exchange
     Second iterative loop is not executed
     Over vegetation, first iterative loop computes vegetation to atmosphere exchange
     Second iterative loop computes ground to vegetation exchange */
  
  // Output
  int rcd(0); // O [rcd] Return success code
  
  // Local
  const prc_cmp eps_dbz(1.0e-6); // [frc] Prevents division by zero
  const prc_cmp rss_sfc_vpr_snw(150.0); // [s m-1] Surface resistance to vapor transfer through snow Bon96 p. 56, Fgr. 17 p. 57
  const prc_cmp svp_H2O_frz(svp_H2O_lqd_PrK78_fst_scl(0.0)); // [Pa] Saturation vapor pressure of H2O at freezing point
  const prc_cmp wnd_cnp_prm_cst(3.0); // [frc] Canopy wind parameter Bon96 p. 63, Bru82 p. 102
  const prc_cmp wnd_mdp_min(1.0); // [m s-1] Minimum surface layer mean wind speed
  const long idx_gnd(1); // [idx] Named index for ground
  const long idx_vgt(0); // [idx] Named index for vegetation
  const long itr_max_gnd(10); // Maximum number of iterations for tpt_gnd loop
  const long itr_max_sfc(20); // Maximum number of iterations for tpt_sfc loop
  const long sgn_chg_ctr_max(4); // Maximum number of sign changes in stability parameter

  prc_cmp cnd_evp_vgt_sfc(CEWI_cpv); // [m s-1] Evaporation conductance vegetation to surface air Bon96 p. 60, Fgr. 17 p. 57
  prc_cmp cnd_heat_gnd_sfc; // [m s-1] Sensible heat conductance ground to surface air
  prc_cmp cnd_heat_sfc_mdp; // [m s-1] Sensible heat conductance surface air to midlayer air
  prc_cmp cnd_heat_ttl; // [m s-1] Sum of conductances
  prc_cmp cnd_heat_vgt_sfc(CEWI_cpv); // [m s-1] Sensible heat conductance vegetation to surface air
  prc_cmp cnd_trn_vgt_sfc(CEWI_cpv); // [m s-1] Transpiration conductance vegetation to surface air Bon96 p. 60, Fgr. 17 p. 57
  std::valarray<prc_cmp> cnd_vpr_cnp(lon_nbr); // [m s-1] Net canopy conductance for vapor
  prc_cmp cnd_vpr_gnd_sfc; // [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
  prc_cmp cnd_vpr_sfc_mdp; // [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57
  prc_cmp cnd_vpr_ttl; // [m s-1] Sum of conductances
  std::valarray<prc_cmp> cst_psych(lon_nbr); // [Pa K-1] Psychrometric constant
  prc_cmp dff_heat_cnp; // [m2 s-1] Eddy diffusivity for heat in canopy 
  std::valarray<prc_cmp> dsvpdt_H2O_gnd(lon_nbr); // [Pa K-1] Derivative of saturation vapor pressure over planar condensed water, ground
  std::valarray<prc_cmp> dsvpdt_H2O_sfc(lon_nbr); // [Pa K-1] Derivative of saturation vapor pressure over planar condensed water, surface
  a2d_cls<prc_cmp> flx_LW_net(lon_nbr,2); // [W m-2] Net longwave flux to atmosphere
  std::valarray<prc_cmp> flx_LW_net_cff_a(lon_nbr); // [W m-2] a in FLWnet = a + b*T^4 Bon96 p. 45
  std::valarray<prc_cmp> flx_LW_net_cff_b(lon_nbr); // [W m-2 K-4] b in FLWnet = a + b*T^4 Bon96 p. 45
  std::valarray<prc_cmp> flx_LW_net_ttl(lon_nbr); // [W m-2] Total net longwave flux to atmosphere
  std::valarray<prc_cmp> flx_LW_msn_sfc(lon_nbr); // [W m-2] Longwave emitted flux at surface
  std::valarray<prc_cmp> flx_LW_rfl_sfc(lon_nbr); // [W m-2] Longwave reflected flux at surface
  a2d_cls<prc_cmp> flx_SW_net(lon_nbr,2); // [W m-2] Solar flux absorbed at surface
  a2d_cls<prc_cmp> flx_sns_gnd(lon_nbr,2); // [W m-2] Sensible heat flux to soil Bon96 p. 64
  std::valarray<prc_cmp> flx_sns_gnd_cff_a(lon_nbr); // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
  std::valarray<prc_cmp> flx_sns_gnd_cff_b(lon_nbr); // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
  a2d_cls<prc_cmp> flx_ltn_evp(lon_nbr,2); // [W m-2] Evaporation flux to atmosphere
  std::valarray<prc_cmp> flx_ltn_evp_cff_a(lon_nbr); // [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> flx_ltn_evp_cff_b(lon_nbr); // [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> flx_ltn_evp_cnp_atm(lon_nbr); // [W m-2] Canopy evaporation flux to atmosphere
  std::valarray<prc_cmp> flx_ltn_evp_gnd_atm(lon_nbr); // [W m-2] Ground evaporation flux to atmosphere
  prc_cmp flx_ltn_fct; // Factor in vapor flux calculations
  a2d_cls<prc_cmp> flx_ltn_trn(lon_nbr,2); // [W m-2] Transpiration flux to atmosphere
  std::valarray<prc_cmp> flx_ltn_trn_cff_a(lon_nbr); // [W m-2] a in LHT = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> flx_ltn_trn_cff_b(lon_nbr); // [W m-2 Pa-1] b in LHT = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> flx_ltn_trn_cnp_atm(lon_nbr); // [W m-2] Canopy transpiration flux to atmosphere
  a2d_cls<prc_cmp> flx_sns_atm(lon_nbr,2); // [W m-2] Sensible heat flux to atmosphere
  std::valarray<prc_cmp> flx_sns_atm_cff_a(lon_nbr); // [W m-2] a in SH = a + b*Ts
  std::valarray<prc_cmp> flx_sns_atm_cff_b(lon_nbr); // [W m-2 K-1] b in SH = a + b*Ts
  std::valarray<prc_cmp> flx_sns_atm_cnp(lon_nbr); // [W m-2] Canopy heat storage (currently defined = 0.0)
  prc_cmp flx_sns_atm_fct; // Factor in sensible heat computation
  prc_cmp flx_sns_atm_tmp(CEWI_cpv); // [W m-2] Temporary sensible heat flux
  prc_cmp flx_sns_atm_vrt_tmp; // [W m-2] Temporary virtual sensible heat flux Bon96 p. 49
  prc_cmp flx_vpr_tmp(CEWI_cpv); // [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
  std::valarray<prc_cmp> hgt_cnp(lon_nbr); // [m] Canopy height 
  std::valarray<prc_cmp> hgt_cnp_bnd(lon_nbr); // [m] Canopy height bounded
  std::valarray<prc_cmp> ltn_heat_trn(lon_nbr); // [J kg-1] Latent heat of sublimation or evaporation
  prc_cmp mno_dnm; // Denominator of Monin-Obukhov length Bon96 p. 49
  std::valarray<prc_cmp> mno_stb_crc_heat(lon_nbr); // [frc] Monin-Obukhov stability correction heat
  prc_cmp mno_stb_crc_heat_crr; // Undamped correction factor heat [frc]
  std::valarray<prc_cmp> mno_stb_crc_mmn(lon_nbr); // [frc] Monin-Obukhov stability correction momentum
  prc_cmp mno_stb_crc_mmn_crr; // Undamped correction factor momentum [frc]
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  prc_cmp mno_stb_crc_tmp2; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp3; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp4; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp5; // Term in stability correction computation
  std::valarray<prc_cmp> mno_stb_prm(lon_nbr); // [frc] Monin-Obukhov stability parameter 
  std::valarray<prc_cmp> mno_stb_prm_old(lon_nbr); // [frc] Monin Obukhov stability parameter old
  std::valarray<prc_cmp> msv_vgt(lon_nbr); // [frc] Vegetation emissivity
  std::valarray<prc_cmp> nrg_bdg(lon_nbr); // [W m-2] Surface energy budget
  prc_cmp nrg_bdg_dlt; // [W m-2 K-1] Temperature derivative of surface energy budget
  std::valarray<prc_cmp> ppr_H2O_cnp(lon_nbr); // [Pa] Canopy vapor pressure of H2O
  std::valarray<prc_cmp> ppr_H2O_cnp_cff_a(lon_nbr); // [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> ppr_H2O_cnp_cff_b(lon_nbr); // [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> ppr_H2O_mdp(lon_nbr); // [Pa] Ambient vapor pressure of H2O
  std::valarray<prc_cmp> rgh_fct(lon_nbr); // [frc] Roughness scale Bon96 p. 55 (binv = B^{-1})
  std::valarray<prc_cmp> rgh_heat(lon_nbr); // [m] Roughness length heat
  std::valarray<prc_cmp> rgh_heat_gnd(lon_nbr); // [m] Roughness length heat ground
  std::valarray<prc_cmp> rgh_mmn_gnd(lon_nbr); // [m] Roughness length momentum ground
  a2d_cls<prc_cmp> rss_aer_heat(lon_nbr,2); // [s m-1] Aerodynamic resistance to heat transfer
  a2d_cls<prc_cmp> rss_aer_mmn(lon_nbr,2); // [s m-1] Aerodynamic resistance to momentum transfer
  a2d_cls<prc_cmp> rss_aer_vpr(lon_nbr,2); // [s m-1] Aerodynamic resistance to vapor transfer
  std::valarray<prc_cmp> rss_aer_fct(lon_nbr); // [s m-1] Term in resistance calculation
  std::valarray<prc_cmp> rss_aer_heat_cnp_fct(lon_nbr); // Term in canopy heat resistance
  prc_cmp rss_aer_heat_cnp_tmp1; // Term in canopy heat resistance
  prc_cmp rss_aer_heat_cnp_tmp2; // Term in canopy heat resistance
  std::valarray<prc_cmp> rss_aer_heat_fct(lon_nbr); // [frc] Term in resistance calculation
  std::valarray<prc_cmp> rss_aer_mmn_fct(lon_nbr); // [frc] Term in resistance calculation
  std::valarray<prc_cmp> rss_leaf_heat_vpr(lon_nbr); // [s m-1] Bulk leaf resistance to heat moisture transfer Bon96 p. 63
  std::valarray<prc_cmp> rss_leaf_heat_vpr_fct(lon_nbr); // Factor to compute bulk leaf resistance
  std::valarray<prc_cmp> rss_sfc_vpr(lon_nbr); // [s m-1] Surface resistance to vapor transfer
  std::valarray<prc_cmp> svp_H2O_gnd(lon_nbr); // [Pa] Saturation vapor pressure over planar condensed water at ground
  std::valarray<prc_cmp> svp_H2O_sfc(lon_nbr); // [Pa] Saturation vapor pressure over planar condensed water at surface
  std::valarray<prc_cmp> tpt_ash_cff_a(lon_nbr); // [K] a in Ta = a + b*Ts
  std::valarray<prc_cmp> tpt_ash_cff_b(lon_nbr); // [frc] b in Ta = a + b*Ts
  prc_cmp tpt_bnd_cls; // [C] Temperature bounded celsius
  std::valarray<prc_cmp> tpt_gnd_dlt(lon_nbr); // [K] Ground temperature adjustment
  prc_cmp tpt_gnd_old; // [K] Previous ground temperature
  std::valarray<prc_cmp> tpt_sfc(lon_nbr); // [K] Surface temperature
  std::valarray<prc_cmp> tpt_sfc_dlt(lon_nbr); // [K] Newton-Raphson temperature adjustment
  prc_cmp tpt_sfc_old; // [K] Previous surface temperature
  std::valarray<prc_cmp> tpt_vgt_dlt(lon_nbr); // [K] Vegetation temperature adjustment
  std::valarray<prc_cmp> tpt_vrt(lon_nbr); // [K] Virtual temperature
  std::valarray<prc_cmp> wnd_cnp_prm(lon_nbr); // [frc] Canopy wind parameter Bon96 p. 63, Bru82 p. 102
  prc_cmp wnd_cnp_top; // [m s-1] Wind speed at canopy top Bon96 p. 63
  std::valarray<prc_cmp> wnd_cnp_top_fct(lon_nbr); // Factor to compute wind at z=ztop 
  std::valarray<prc_cmp> wnd_mdp(lon_nbr); // [m s-1] Surface layer mean wind speed
  std::valarray<prc_cmp> wnd_mdp_bnd(lon_nbr); // [m s-1] Surface layer mean wind speed bounded
  long lon_idx; // Counting index for lon
  long lvl_idx(CEWI_lng); // [idx] Level index, either idx_vgt or idx_gnd
  std::valarray<long> sgn_chg_ctr(lon_nbr); // Number of sign changes in stability parameter
  
  // Main Code
  const std::string sbr_nm("flx_sfc_lnd"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  //  const std::string dbg_sng_lcl(dbg_sng_get()); // [sng] Debugging string
  
  // Initialize variables which normally would be available in LSM
  vec_set(&wnd_cnp_prm[0],lon_nbr,wnd_cnp_prm_cst); // [frc] Canopy wind parameter Bon96 p. 63, Bru82 p. 102
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    ppr_H2O_mdp[lon_idx]=q_H2O_vpr[lon_idx]*prs_mdp[lon_idx]/(eps_H2O+one_mns_eps_H2O*q_H2O_vpr[lon_idx]); // [Pa] Ambient vapor pressure of H2O
    
    hgt_cnp[lon_idx]=hvt[pnt_typ_idx[lon_idx]]; // [m] Canopy height 
    
    // Default ground emissivity
    msv_sfc[lon_idx]=msv_snw[lon_idx]*snw_frc[lon_idx]+
      msv_gnd[lon_idx]*(1.0-snw_frc[lon_idx]); // [frc] Surface (bare ground+snow) emissivity

    msv_vgt[lon_idx]=msv_sfc[lon_idx]; // [frc] Vegetation emissivity fxm: implement plant type dependence

    mno_stb_prm[lon_idx]=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    
    /* Notes on fluxes:
       For vegetated surfaces, the total flux to the atmosphere is the sum of
       the flux from the vegetation and the flux from the ground.
       Bonan makes the assumption that the canopy does not store heat or moisture
       Without this assumption, the fluxes would not be additive, as some of the
       ground flux of moisture, e.g., would be trapped in the canopy before 
       escaping to the atmosphere; likewise with sensible and latent heat.
     */
    /* Notes on roughness lengths:
       The roughness lengths for bare ground are always required
       For vegetated surfaces, rgh_mmn and rgh_heat contain the 
       respective roughness lengths of the corresponding vegetation types 
       rgh_mmn_gnd and rgh_heat_gnd, on the other hand, contain the bare
       ground roughness lengths
       For non-vegetated surfaces, rgh_mmn = rgh_mmn_gnd, rgh_heat = rgh_heat_gnd
    */
    // Bon96 currently uses rgh_fct = 0.0
    rgh_fct[lon_idx]=0.0; // [frc] Roughness scale Bon96 p. 55 (binv = B^{-1})
    // rgh_fct[lon_idx]=cst_von_krm_rcp*log(rgh_mmn[lon_idx]/rgh_heat[lon_idx]); // [frc] Roughness scale Bon96 p. 55 (binv = B^{-1})
    
    //    const prc_cmp rgh_mmn_snw(0.04); // [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
    //    rgh_mmn_gnd[lon_idx]=rgh_mmn_snw*snw_frc[lon_idx]+
    //  rgh_lng[soi_typ[lon_idx]]*(1.0-snw_frc[lon_idx]); // [m] Roughness length momentum ground Bon96 p. 59
    rgh_mmn_gnd[lon_idx]=rgh_mmn[lon_idx]; // [m] Roughness length momentum ground Bon96 p. 59 fxm: implement CCM:lsm/surphy() 
    
    rgh_heat[lon_idx]=rgh_mmn[lon_idx]*std::exp(-cst_von_krm*rgh_fct[lon_idx]); // [m] Bon96 p. 59 CCM:lsm/surphy()
    rgh_heat_gnd[lon_idx]=rgh_mmn_gnd[lon_idx]*std::exp(-cst_von_krm*rgh_fct[lon_idx]); // [m] Bon96 p. 59 CCM:lsm/surphy()
    
    tpt_vrt[lon_idx]=tpt_mdp[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Virtual temperature
    
  }  // end loop over lon
  
  // Initialize variables which are independent of stability iteration
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    // Zero temperature adjustments from last timestep
    tpt_vgt_dlt[lon_idx]=0.0; // [K] Change in vegetation temperature
    tpt_gnd_dlt[lon_idx]=0.0; // [K] Change in ground temperature
    
    // Latent heat of water transformation
    tpt_bnd_cls=tpt_bnd_cls_get(tpt_mdp[lon_idx]); // [C]
    if(tpt_bnd_cls > 0.0){ 
      ltn_heat_trn[lon_idx]=ltn_heat_vpr_H2O_std; // [J kg-1]
    }else{ 
      ltn_heat_trn[lon_idx]=ltn_heat_sbl_H2O_std; // [J kg-1]
    } // endif
    // Psychrometric constant for transformations of water vapor at surface
    cst_psych[lon_idx]=spc_heat_dry_air*prs_mdp[lon_idx]/(eps_H2O*ltn_heat_trn[lon_idx]); // [Pa K-1] 
    
    // Combine separate SW absorbed fluxes into single array
    flx_SW_net(lon_idx,idx_vgt)=flx_SW_net_vgt[lon_idx]; // [W m-2] Solar flux absorbed by vegetation
    flx_SW_net(lon_idx,idx_gnd)=flx_SW_net_gnd[lon_idx]; // [W m-2] Solar flux absorbed by vegetation
    
    // tpt_sfc is the temperature we are solving for, tpt_vgt if vegetated, tpt_gnd otherwise
    if(vgt[lon_idx]) tpt_sfc[lon_idx]=tpt_vgt[lon_idx]; else tpt_sfc[lon_idx]=tpt_gnd[lon_idx];
    
    // Saturation vapor pressure of water at ground temperature
    tpt_bnd_cls=tpt_bnd_cls_get(tpt_gnd[lon_idx]); // [C]
    if(tpt_bnd_cls > 0.0){ 
      svp_H2O_gnd[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
    }else{
      svp_H2O_gnd[lon_idx]=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
    } // endif
    
    // Vegetated fluxes (for non-vegetated surfaces) are set here and only here
    flx_sns_gnd(lon_idx,idx_vgt)=0.0; // [W m-2] Sensible heat flux to soil
    flx_sns_atm(lon_idx,idx_vgt)=0.0; // [W m-2] Sensible heat flux to atmosphere
    flx_ltn_trn(lon_idx,idx_vgt)=0.0; // [W m-2] Transpiration flux to atmosphere
    flx_ltn_evp(lon_idx,idx_vgt)=0.0; // [W m-2] Evaporation flux to atmosphere
    flx_LW_net(lon_idx,idx_vgt)=0.0; // [W m-2] Net longwave flux to atmosphere
  } // end loop over lon
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    // Canopy height: must be greater than roughness length of bare ground
    hgt_cnp_bnd[lon_idx]=max_cpv(hgt_cnp[lon_idx],rgh_mmn_gnd[lon_idx]); // [m] Canopy height bounded 
    
    // Midlayer wind speeds
    wnd_mdp[lon_idx]= // [m s-1] Surface layer mean wind speed
      std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	   wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
    wnd_mdp_bnd[lon_idx]=max_cpv(wnd_mdp[lon_idx],wnd_mdp_min); // [m s-1] Surface layer mean wind speed bounded
    
    // Neutral drag coefficient between z2=hgt_mdp and z1=z0m+zpd 
    cff_xch_mmn_ntr[lon_idx]= // [frc]
      std::pow(cst_von_krm/std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_mmn[lon_idx]),2.0);
    
    // Miscellaneous
    sgn_chg_ctr[lon_idx]=0; // Number of sign changes in stability parameter
    mno_stb_prm_old[lon_idx]=0.0; // [frc] Monin Obukhov stability parameter old
    rss_aer_heat_cnp_tmp1= // Bon96 p. 63
      std::exp(-wnd_cnp_prm[lon_idx]*rgh_heat_gnd[lon_idx]/hgt_cnp_bnd[lon_idx]);
    rss_aer_heat_cnp_tmp2= // Bon96 p. 63
      std::exp(-wnd_cnp_prm[lon_idx]*(rgh_heat[lon_idx]+hgt_zpd[lon_idx])/hgt_cnp_bnd[lon_idx]);
    rss_aer_heat_cnp_fct[lon_idx]= // Term in canopy heat resistance Bon96 p. 63 
      hgt_cnp_bnd[lon_idx]
      *std::exp(wnd_cnp_prm[lon_idx])
      *(rss_aer_heat_cnp_tmp1-rss_aer_heat_cnp_tmp2)
      /wnd_cnp_prm[lon_idx];
    wnd_cnp_top_fct[lon_idx]= // Factor to compute wind at z=ztop Bon96 p. 63 
      std::log((hgt_cnp_bnd[lon_idx]-hgt_zpd[lon_idx])/rgh_mmn[lon_idx])
      /cst_von_krm;
    rss_leaf_heat_vpr_fct[lon_idx]= // Factor to compute bulk canopy resistance Bon96 p. 63 
      wnd_cnp_prm[lon_idx]*50.0
      /(1.0-std::exp(-wnd_cnp_prm[lon_idx]/2.0));
    
    // Stability-independent components of aerodynamic resistance calculations
    rss_aer_fct[lon_idx]=1.0/(cst_von_krm*cst_von_krm*wnd_mdp_bnd[lon_idx]); // [s m-1]
    rss_aer_mmn_fct[lon_idx]=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_mmn[lon_idx]); // [frc]
    rss_aer_heat_fct[lon_idx]=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_heat[lon_idx]); // [frc]
  } // end loop over lon
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    if(vgt[lon_idx]){
      flx_LW_net_cff_a[lon_idx]= // [W m-2] a in FLWnet = a + b*Tg^4 Bon96 p. 45
	-msv_vgt[lon_idx]
	*(1.0+(1.0-msv_vgt[lon_idx])*(1.0-msv_sfc[lon_idx]))
	*flx_LW_dwn_sfc[lon_idx]
	-msv_vgt[lon_idx]*msv_sfc[lon_idx]*cst_Stefan_Boltzmann
	*std::pow(tpt_gnd[lon_idx],PRC_CMP(4.0));
      flx_LW_net_cff_b[lon_idx]= // [W m-2 K-4] b in FLWnet = a + b*Tg^4 Bon96 p. 45
	+(2.0-msv_vgt[lon_idx]*(1.0-msv_sfc[lon_idx]))
	*msv_vgt[lon_idx]
	*cst_Stefan_Boltzmann;
      // Canopy heat flux directly to sub-surface soil is zero, since all heat must pass through ground level first
      flx_sns_gnd_cff_b[lon_idx]=0.0; // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
      flx_sns_gnd_cff_a[lon_idx]=0.0; // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
    }else{
      /* 20210618
	 F(LW net) = LWupw - LWdwn = Reflected up + Emitted up - Absorbed down
	 Assuming absorptivity = emissivity then reflected = 1 - msv_sfc and
	 LWupw = (1.0-msv)*Fdwn + msv*sigma*Tg^4
	 LWdwn = LWdwn
	 F(LW net) = (1.0-msv)*Fdwn + msv*sigma*Tg^4 - Fdwn
	           = Fdwn -msv*Fdwn + msv*sigma*Tg^4 - Fdwn
                   = -msv*Fdwn + msv*sigma*Tg^4 
		   = a + b*Tg^4 Bon96 p. 45
	 where
	   a = -msv*Fdwn
	   b = msv*sigma
	 Formula may look physically incorrect, yet is correct due to cancellation of Fdwn
	 Inverting for Fupw in terms of Fnet,

	 Fupw = Fnet + Fdwn 
	      = -msv*Fdwn + msv*sigma*Tg^4 + Fdwn
	      = (1.0-msv)*Fdwn + msv*sigma*Tg^4
              = Reflected Up + Emitted Up */
      flx_LW_net_cff_a[lon_idx]=-msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] a in FLWnet = a + b*Tg^4 Bon96 p. 45
      flx_LW_net_cff_b[lon_idx]=msv_sfc[lon_idx]*cst_Stefan_Boltzmann; // [W m-2 K-4] b in FLWnet = a + b*Tg^4 Bon96 p. 45
      
      // F(sns heat dwn into soil) = 2*k*(Tg-T1)/dz = a + b*Tg Bon96 p. 64
      flx_sns_gnd_cff_b[lon_idx]=2.0*cnd_trm_soi[lon_idx]/lvl_dlt_snw[lon_idx]; // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
      flx_sns_gnd_cff_a[lon_idx]=-flx_sns_gnd_cff_b[lon_idx]*tpt_soi[lon_idx]; // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
    } // endelse
  }  // end loop over lon
  
  // Iteration loop
  const prc_cmp eps_max(1.0e-5); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  long itr_idx; // Counting index
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    if(dbg_lvl == dbg_crr){
      (void)std::fprintf(stderr,"%s fluxes from %s():\n",(vgt[lon_idx] ? "Vegetation" : "Ground"),sbr_nm.c_str());
      (void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		    "itr","mno_lng","mno_stb","wnd_frc","  Ts   ","LW(msn)"," H(atm)"," H(soi)","   L   ","nrg_bdg","  eps  ");
      (void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		    "    ","   m   ","  frc  "," m s-1 ","   K   "," W m-2 "," W m-2 "," W m-2 "," W m-2 "," W m-2 ","  frc  ");
    } // end if dbg
    
    // This iteration loop solves for tpt_vgt if vegetated, else tpt_gnd
    while(eps_crr > eps_max){
      
      // Save old surface temperature for convergence diagnostic
      if(itr_idx == 1) tpt_sfc_old=tpt_gnd[lon_idx]; else tpt_sfc_old=tpt_sfc[lon_idx]; // [K]
      
      // Stability functions computed as in Bon96 p. 52
      if(mno_stb_prm[lon_idx] < 0.0){
	sml_fnc_mmn_uns_rcp=std::pow(PRC_CMP(1.0)-PRC_CMP(16.0)*mno_stb_prm[lon_idx],PRC_CMP(0.25));
	mno_stb_crc_tmp2=std::log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_tmp3=std::log((1.0+sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_mmn_crr=2.0*mno_stb_crc_tmp3+mno_stb_crc_tmp2-2.0*std::atan(sml_fnc_mmn_uns_rcp)+1.5707963; // [frc]
	mno_stb_crc_heat_crr=2.0*mno_stb_crc_tmp2; // [frc]
      }else{ // not stable
	mno_stb_crc_mmn_crr=-5.0*mno_stb_prm[lon_idx]; // [frc]
	mno_stb_crc_heat_crr=mno_stb_crc_mmn[lon_idx]; // [frc]
      } // not stable
      
      // Filter stability corrections to reduce numerical ping-pong
      if(itr_idx == 1){
	mno_stb_crc_mmn[lon_idx]=mno_stb_crc_mmn_crr; // [frc]
	mno_stb_crc_heat[lon_idx]=mno_stb_crc_heat_crr; // [frc]
      }else{
	mno_stb_crc_mmn[lon_idx]=0.5*(mno_stb_crc_mmn_crr+mno_stb_crc_mmn[lon_idx]); // [frc]
	mno_stb_crc_heat[lon_idx]=0.5*(mno_stb_crc_heat_crr+mno_stb_crc_heat[lon_idx]); // [frc]
      } // endif first iteration
      
      // Aerodynamic resistance between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
      mno_stb_crc_tmp4=rss_aer_mmn_fct[lon_idx]-mno_stb_crc_mmn[lon_idx]; // [frc]
      mno_stb_crc_tmp5=rss_aer_heat_fct[lon_idx]-mno_stb_crc_heat[lon_idx]; // [frc]
      rss_aer_mmn(lon_idx,idx_vgt)=max_cpv(rss_aer_fct[lon_idx]*mno_stb_crc_tmp4*mno_stb_crc_tmp4,1.0); // [s m-1]
      rss_aer_heat(lon_idx,idx_vgt)=max_cpv(rss_aer_fct[lon_idx]*mno_stb_crc_tmp4*mno_stb_crc_tmp5,1.0); // [s m-1]
      // Resistances are equal because rgh_heat = rgh_vpr
      rss_aer_vpr(lon_idx,idx_vgt)=rss_aer_heat(lon_idx,idx_vgt); // [s m-1]
      
      // Exchange coefficients between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
      // Exchange coefficients are dimensionless, inversely proportional to wind speed 
      cff_xch_mmn[lon_idx]=1.0/(rss_aer_mmn(lon_idx,idx_vgt)*wnd_mdp_bnd[lon_idx]); // [frc]
      cff_xch_heat[lon_idx]=1.0/(rss_aer_heat(lon_idx,idx_vgt)*wnd_mdp_bnd[lon_idx]); // [frc]
      cff_xch_vpr[lon_idx]=1.0/(rss_aer_vpr(lon_idx,idx_vgt)*wnd_mdp_bnd[lon_idx]); // [frc]
      
      // Friction velocity
      wnd_frc[lon_idx]=wnd_mdp_bnd[lon_idx]*std::sqrt(cff_xch_mmn[lon_idx]); // [m s-1]
      
      /* Aerodynamic resistances for heat and moisture between z=zpd+z0h and z=z0hg
	 If no vegetation, rss_aer_heat[idx_gnd]=0 because zpd+z0h = z0hg */
      dff_heat_cnp=max_cpv(cst_von_krm*wnd_frc[lon_idx]*(hgt_cnp_bnd[lon_idx]-hgt_zpd[lon_idx]),eps_dbz); // [m2 s-1] Eddy diffusivity for heat in canopy
      rss_aer_mmn(lon_idx,idx_gnd)=0.0; // [s m-1]
      rss_aer_heat(lon_idx,idx_gnd)=rss_aer_heat_cnp_fct[lon_idx]/dff_heat_cnp; // [s m-1] Bon96 p. 63
      // Resistances are equal because rgh_heat = rgh_vpr Bon96 p. 63
      rss_aer_vpr(lon_idx,idx_gnd)=rss_aer_heat(lon_idx,idx_gnd); // [s m-1]
      
      // Bulk leaf resistance to moisture transfer
      wnd_cnp_top=max_cpv(wnd_frc[lon_idx]*wnd_cnp_top_fct[lon_idx],eps_dbz); // [m s-1] Wind speed at canopy top Bon96 p. 63
      rss_leaf_heat_vpr[lon_idx]=rss_leaf_heat_vpr_fct[lon_idx]*std::sqrt(dleaf[pnt_typ_idx[lon_idx]]/wnd_cnp_top); // [s m-1] Bulk leaf resistance to heat and moisture transfer Bon96 p. 63
      
      // Saturation vapor pressure of water at surface temperature
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc[lon_idx]); // [C]
      if(tpt_bnd_cls > 0.0){ 
	svp_H2O_sfc[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	dsvpdt_H2O_sfc[lon_idx]=dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
      }else{
	svp_H2O_sfc[lon_idx]=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	dsvpdt_H2O_sfc[lon_idx]=dsvpdt_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
      } // endif frozen
      
      // Surface resistance to vapor transfer
      if(itr_idx == 1){ // NB: Bon96 only computes this on the first iteration
	if(vgt[lon_idx]){
	  rss_sfc_vpr[lon_idx]= // [s m-1] Surface resistance to vapor transfer Bon96 p. 56, Fgr. 17 p. 57
	    +rss_sfc_vpr_snw*snw_frc[lon_idx]
	    +(1.0-snw_frc[lon_idx])
	    *rss_aer_vpr(lon_idx,idx_gnd)
	    *(1.0-trn_fsh_vpr_soi_atm[lon_idx])/
	    max_cpv(trn_fsh_vpr_soi_atm[lon_idx],eps_dbz);
	}else{
	  rss_sfc_vpr[lon_idx]= // [s m-1] Surface resistance to vapor transfer Bon96 p. 56, Fgr. 17 p. 57
	    +rss_sfc_vpr_snw*snw_frc[lon_idx]
	    +(1.0-snw_frc[lon_idx])
	    *rss_aer_vpr(lon_idx,idx_vgt)
	    *(1.0-trn_fsh_vpr_soi_atm[lon_idx])/
	    max_cpv(trn_fsh_vpr_soi_atm[lon_idx],eps_dbz);
	} // endif vegetated
	// Set minimum rss_sfc_vpr so that in bare ground case es = svp(tg)
	rss_sfc_vpr[lon_idx]=max_cpv(rss_sfc_vpr[lon_idx],eps_dbz);
	
	// Stomatal resistances for sunlit and shaded fractions of canopy fxm: Add stomatal routines
      } // endif first iteration
      
      /* Relationships between various temperatures: Bon96 p. 55 and Fig. 16 p. 57
	 T1 = Skin temperature = Soil temperature of 1st layer (top 10 cm)
	 Tg = "Ground" temperature = Air temperature at z = rgh_heat
	 Ts = "Surface" temperature = Air temperature at z=zpd+rgh_heat
	 Ta = "Aerodynamic" temperature = Air temperature at z=zpd+rgh_mmn
	 Te = Emitting temperature = (FLWup/sigma)**0.25
	 Tgcm = Ambient temperature = Air temperature at z=hgt_mdp
	 
	 For bare ground, Ts = Tg = Air temperature at z=rgh_heat (Bon96 p. 55)
	 For bare ground, Tg = potential temperture defined relative to local surface pressure
	 
	 Heat conductances from ground and canopy to ambient air
	 Conductances are dimensional, exact inverses of resistances
	 Conductances are much more complicated when ground is vegetated than bare
	 Coefficients for temperature at apparent sink for heat Tash = a + b*Ts
	 Coefficients for sensible heat flux SH = a + b*Ts */
      if(vgt[lon_idx]){
	cnd_heat_sfc_mdp=1.0/rss_aer_heat(lon_idx,idx_vgt); // [m s-1] Sensible heat conductance surface air to midlayer air Bon96 p. 60, Fig. 16 p. 57
	//	cnd_heat_vgt_sfc=2.0*vai[pnt_typ_idx[lon_idx]]/rss_leaf_heat_vpr[lon_idx]; // [m s-1] Sensible heat conductance vegetation to surface air Bon96 p. 60, Fig. 16 p. 57 fxm: Add vai to lsm.hh
	cnd_heat_gnd_sfc=1.0/rss_aer_heat(lon_idx,idx_gnd); // [m s-1] Sensible heat conductance ground to surface air Bon96 p. 60, Fig. 16 p. 57
	cnd_heat_ttl=cnd_heat_sfc_mdp+cnd_heat_vgt_sfc+cnd_heat_gnd_sfc; // [m s-1] Sum of conductances
	tpt_ash_cff_a[lon_idx]= // [K] a in Tash = a + b*Ts Bon96 p. 55
	  +(tpt_ptn_mdp[lon_idx]*cnd_heat_sfc_mdp+tpt_gnd[lon_idx]*cnd_heat_gnd_sfc)
	  /cnd_heat_ttl;
	tpt_ash_cff_b[lon_idx]=cnd_heat_vgt_sfc/cnd_heat_ttl; // [frc] b in Tash = a + b*Ts Bon96 p. 55
	flx_sns_atm_fct=dns_mdp[lon_idx]*spc_heat_dry_air*cnd_heat_vgt_sfc; // Bon96 p. 55
	flx_sns_atm_cff_a[lon_idx]=-tpt_ash_cff_a[lon_idx]*flx_sns_atm_fct; // [W m-2] a in SH = a + b*Ts Bon96 p. 55, 69, 70
	flx_sns_atm_cff_b[lon_idx]=(1.0-tpt_ash_cff_b[lon_idx])*flx_sns_atm_fct; // [W m-2 K-1] b in SH = a + b*Ts Bon96 p. 55, 69, 70
      }else{
	cnd_heat_sfc_mdp=1.0/rss_aer_heat(lon_idx,idx_vgt); // [m s-1] Sensible heat conductance surface air to midlayer air Bon96 p. 60, Fig. 16 p. 57
	cnd_heat_vgt_sfc=0.0; // [m s-1] Sensible heat conductance vegetation to surface air Bon96 p. 60, Fig. 16 p. 57
	cnd_heat_gnd_sfc=0.0; // [m s-1] Sensible heat conductance ground to surface air Bon96 p. 60, Fig. 16 p. 57
	cnd_heat_ttl=cnd_heat_sfc_mdp+cnd_heat_vgt_sfc+cnd_heat_gnd_sfc; // [m s-1] Sum of conductances
	// For bare ground, Ts = Tg = Air temperature at z=rgh_heat (Bon96 p. 55)
	tpt_ash_cff_a[lon_idx]=0.0; // [K] a in Tash = a + b*Ts Bon96 p. 55
	tpt_ash_cff_b[lon_idx]=1.0; // [frc] b in Tash = a + b*Ts Bon96 p. 55
	flx_sns_atm_fct=dns_mdp[lon_idx]*spc_heat_dry_air*cnd_heat_sfc_mdp; // Bon96 p. 55
	flx_sns_atm_cff_a[lon_idx]=-tpt_ptn_mdp[lon_idx]*flx_sns_atm_fct; // [W m-2] a in SH = a + b*Ts Bon96 p. 55, 69, 70
	flx_sns_atm_cff_b[lon_idx]=flx_sns_atm_fct; // [W m-2 K-1] b in SH = a + b*Ts Bon96 p. 55, 69, 70
      } // endif vegetated
      
      // Vapor conductances from ground and canopy to ambient air
      flx_ltn_fct=dns_mdp[lon_idx]*spc_heat_dry_air/cst_psych[lon_idx]; // 
      if(vgt[lon_idx]){
      }else{
	cnd_vpr_sfc_mdp=1.0/rss_aer_vpr(lon_idx,idx_vgt); // [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57
	cnd_evp_vgt_sfc=0.0; // [m s-1] Evaporation conductance vegetation to surface air Bon96 p. 60, Fgr. 17 p. 57
	cnd_trn_vgt_sfc=0.0; // [m s-1] Transpiration conductance vegetation to surface air Bon96 p. 60, Fgr. 17 p. 57
	cnd_vpr_gnd_sfc=1.0/rss_sfc_vpr[lon_idx]; // [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
	cnd_vpr_ttl=cnd_vpr_sfc_mdp+cnd_evp_vgt_sfc+cnd_trn_vgt_sfc+cnd_vpr_gnd_sfc; // [m s-1] Sum of conductances
	
	// Coefficients for canopy vapor pressure e(cnp) = a + b*e(Ts)
	ppr_H2O_cnp_cff_a[lon_idx]=ppr_H2O_mdp[lon_idx]*cnd_vpr_sfc_mdp/cnd_vpr_ttl; // [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
	ppr_H2O_cnp_cff_b[lon_idx]=cnd_vpr_gnd_sfc/cnd_vpr_ttl; // [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55
	// Coefficients for evaporation heat flux LHE = a + b*e(Ts) 
	flx_ltn_evp_cff_a[lon_idx]=-flx_ltn_fct*(ppr_H2O_mdp[lon_idx]-ppr_H2O_cnp_cff_a[lon_idx])*cnd_vpr_sfc_mdp; // [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
	flx_ltn_evp_cff_b[lon_idx]=flx_ltn_fct*ppr_H2O_cnp_cff_b[lon_idx]*cnd_vpr_sfc_mdp; // [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
	// Coefficients for transpiration heat flux LHT = a + b*e(Ts)
	flx_ltn_trn_cff_a[lon_idx]=0.0; // [W m-2] a in LHT = a + b*e(Ts) Bon96 p. 55
	flx_ltn_trn_cff_b[lon_idx]=0.0; // [W m-2 Pa-1] b in LHT = a + b*e(Ts) Bon96 p. 55
      } // endif vegetated
      cnd_vpr_cnp[lon_idx]=cnd_evp_vgt_sfc+cnd_trn_vgt_sfc; // [m s-1] Net canopy conductance for vapor
      
      // Evaluate fluxes for current tpt_sfc
      if(vgt[lon_idx]) lvl_idx=idx_vgt; else lvl_idx=idx_gnd;
      flx_sns_gnd(lon_idx,lvl_idx)=flx_sns_gnd_cff_a[lon_idx]+flx_sns_gnd_cff_b[lon_idx]*tpt_sfc[lon_idx]; // [W m-2] Sensible heat flux to soil
      flx_sns_atm(lon_idx,lvl_idx)=flx_sns_atm_cff_a[lon_idx]+flx_sns_atm_cff_b[lon_idx]*tpt_sfc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn_trn(lon_idx,lvl_idx)=flx_ltn_trn_cff_a[lon_idx]+flx_ltn_trn_cff_b[lon_idx]*svp_H2O_sfc[lon_idx]; // [W m-2] Transpiration flux to atmosphere
      flx_ltn_evp(lon_idx,lvl_idx)=flx_ltn_evp_cff_a[lon_idx]+flx_ltn_evp_cff_b[lon_idx]*svp_H2O_sfc[lon_idx]; // [W m-2] Evaporation flux to atmosphere
      flx_LW_net(lon_idx,lvl_idx)=flx_LW_net_cff_a[lon_idx]+flx_LW_net_cff_b[lon_idx]*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Net longwave flux to atmosphere
      nrg_bdg[lon_idx]= // [W m-2] Total energy budget at surface
	+flx_SW_net(lon_idx,lvl_idx)
	-flx_LW_net(lon_idx,lvl_idx)
	-flx_ltn_evp(lon_idx,lvl_idx)
	-flx_ltn_trn(lon_idx,lvl_idx)
	-flx_sns_atm(lon_idx,lvl_idx)
	-flx_sns_gnd(lon_idx,lvl_idx);
      nrg_bdg_dlt= // [W m-2 K-1] Temperature derivative of surface energy budget
	-4.0*flx_LW_net_cff_b[lon_idx]*std::pow(tpt_sfc[lon_idx],PRC_CMP(3.0))
	-(flx_ltn_evp_cff_b[lon_idx]+flx_ltn_trn_cff_b[lon_idx])*dsvpdt_H2O_sfc[lon_idx]
	-flx_sns_atm_cff_b[lon_idx]
	-flx_sns_gnd_cff_b[lon_idx];
      tpt_sfc_dlt[lon_idx]=-nrg_bdg[lon_idx]/nrg_bdg_dlt; // [K] Newton-Raphson temperature adjustment
      
      // Adjust fluxes before adjusting temperatures
      if(vgt[lon_idx]) lvl_idx=idx_vgt; else lvl_idx=idx_gnd;
      flx_sns_gnd(lon_idx,lvl_idx)+=flx_sns_gnd_cff_b[lon_idx]*tpt_sfc_dlt[lon_idx]; // [W m-2] Sensible heat flux to soil
      flx_sns_atm(lon_idx,lvl_idx)+=flx_sns_atm_cff_b[lon_idx]*tpt_sfc_dlt[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn_trn(lon_idx,lvl_idx)+=flx_ltn_trn_cff_b[lon_idx]*dsvpdt_H2O_sfc[lon_idx]*tpt_sfc_dlt[lon_idx]; // [W m-2] Transpiration flux to atmosphere
      flx_ltn_evp(lon_idx,lvl_idx)+=flx_ltn_evp_cff_b[lon_idx]*dsvpdt_H2O_sfc[lon_idx]*tpt_sfc_dlt[lon_idx]; // [W m-2] Evaporation flux to atmosphere
      flx_LW_net(lon_idx,lvl_idx)+=4.0*flx_LW_net_cff_b[lon_idx]*std::pow(tpt_sfc[lon_idx],PRC_CMP(3.0))*tpt_sfc_dlt[lon_idx]; // [W m-2] Net longwave flux to atmosphere
      
      // Adjust temperatures 
      tpt_sfc[lon_idx]+=tpt_sfc_dlt[lon_idx]; // [K]
      tpt_ash[lon_idx]=tpt_ash_cff_a[lon_idx]+tpt_ash_cff_b[lon_idx]*tpt_sfc[lon_idx]; // [K] "Surface" temperature at z=zpd+rgh_heat
      
      // Adjust canopy vapor pressure
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc[lon_idx]); // [C]
      if(tpt_bnd_cls > 0.0){
	ppr_H2O_cnp[lon_idx]=ppr_H2O_cnp_cff_a[lon_idx]+ppr_H2O_cnp_cff_b[lon_idx]*svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      }else{
	ppr_H2O_cnp[lon_idx]=ppr_H2O_cnp_cff_a[lon_idx]+ppr_H2O_cnp_cff_b[lon_idx]*svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      } // endif frozen
      
      // Sensible heat flux
      flx_sns_atm_tmp= // [W m-2] Temporary sensible heat flux
	-(tpt_ptn_mdp[lon_idx]-tpt_ash[lon_idx])
	*dns_mdp[lon_idx]*spc_heat_dry_air
	/rss_aer_heat(lon_idx,idx_vgt);
      // Following step approximates psi_h(z0m/L) = 0 
      // NB: Is there a simpler way to get tpt_aer than this...
      tpt_aer[lon_idx]= // [K] "Aerodynamic" temperature at z=zpd+rgh_mmn Bon96 p. 55
	+tpt_ash[lon_idx]
	-flx_sns_atm_tmp*rgh_fct[lon_idx]
	/(dns_mdp[lon_idx]*spc_heat_dry_air*wnd_frc[lon_idx]);
      
      // Monin-Obukhov stability parameter mno_stb_prm for next iteration
      flx_vpr_tmp= // [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
	// 19990501 fxm: this appears to be latent heat flux, not water vapor flux
	-(ppr_H2O_mdp[lon_idx]-ppr_H2O_cnp[lon_idx])
	*dns_mdp[lon_idx]*spc_heat_dry_air
	/(cst_psych[lon_idx]*rss_aer_vpr(lon_idx,idx_vgt));
      flx_sns_atm_vrt_tmp= // [W m-2] Virtual sensible heat flux Bon96 p. 49
	+flx_sns_atm_tmp
	+eps_H2O_rcp_m1*spc_heat_dry_air*tpt_mdp[lon_idx]*flx_vpr_tmp
	/ltn_heat_trn[lon_idx];
      mno_dnm=  // Denominator of Monin-Obukhov length Bon96 p. 49
	+cst_von_krm
	*(grv_sfc_mean/tpt_vrt[lon_idx])
	*flx_sns_atm_vrt_tmp
	/(dns_mdp[lon_idx]*spc_heat_dry_air);
      // Set denominator of Monin-Obukhov length to minimum value if vapor and heat fluxes equal 0.0
      if(PRC_CMP_ABS(mno_dnm) <= eps_dbz) mno_dnm=eps_dbz;
      mno_lng[lon_idx]=-1.0*std::pow(wnd_frc[lon_idx],PRC_CMP(3.0))/mno_dnm; // [m] Monin-Obukhov length Bon96 p. 49
      // Stability functions only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
      mno_stb_prm[lon_idx]=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc] Monin Obukhov stability parameter
      
      // Accumulate number of times stability parameter changes sign 
      if(mno_stb_prm_old[lon_idx]*mno_stb_prm[lon_idx] < 0.0) sgn_chg_ctr[lon_idx]++;
      // Zero stability parameter if it has changed sign too many times
      if(sgn_chg_ctr[lon_idx] >= sgn_chg_ctr_max){
	wrn_prn(prg_nm,sbr_nm,"Numerical ping-pong...setting mno_stb_prm to 0.0");
	mno_stb_prm[lon_idx]=0.0; // [frc]
	// Zero stability corrections for consistency with stability parameter
	mno_stb_crc_mmn[lon_idx]=0.0; // [frc]
	mno_stb_crc_heat[lon_idx]=0.0; // [frc]
      } // endif
      mno_stb_prm_old[lon_idx]=mno_stb_prm[lon_idx]; // [frc]
      
      eps_crr=PRC_CMP_ABS((tpt_sfc[lon_idx]-tpt_sfc_old)/tpt_sfc[lon_idx]); // Relative convergence
      if(dbg_lvl == dbg_crr){
	// NB: Fluxes reported here have not been adjusted for snowmelt yet
	// NB: LW up flux here assumes non-vegetated fxm
	// 20210618 Emissivity factor was incorrectly applied to Fdwn in below line for ~25 years
	//flx_LW_upw_sfc[lon_idx]=flx_LW_net(lon_idx,lvl_idx)+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44
	flx_LW_upw_sfc[lon_idx]=flx_LW_net(lon_idx,lvl_idx)+flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44
	flx_LW_rfl_sfc[lon_idx]=(1.0-msv_sfc[lon_idx])*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave reflected flux at surface
	flx_LW_msn_sfc[lon_idx]=flx_LW_upw_sfc[lon_idx]-flx_LW_rfl_sfc[lon_idx]; // [W m-2] Longwave emitted flux at surface
	flx_sns_gnd_ttl[lon_idx]=flx_sns_gnd(lon_idx,lvl_idx); // [W m-2] Sensible heat flux to soil
	(void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
		      itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],tpt_sfc[lon_idx],flx_LW_msn_sfc[lon_idx],flx_sns_atm_tmp,flx_sns_gnd_ttl[lon_idx],flx_vpr_tmp,nrg_bdg[lon_idx],eps_crr);
      } // end if dbg
      if(itr_idx > itr_max_sfc){
	std::cerr << "Final: tpt_sfc = " << tpt_sfc[lon_idx] << " K, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Vegetation temperature not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
    
    /* Fill in vegetation and ground temperatures and delta temperatures from ts and dts
       If vegetated: tpt_gnd and tpt_gnd_dlt will be set to correct values later
       If not vegetated: set tpt_vgt and tpt_vgt_dlt to ground values
       Surface can switch from vegetated to non-vegetated if buried by snow 
       Need to set tpt_vgt=tpt_gnd for non-vegetated surfaces so that tpt_vgt is
       available if vegetation is exposed on next time step */
    tpt_vgt[lon_idx]=tpt_sfc[lon_idx]; // [K]
    tpt_vgt_dlt[lon_idx]=tpt_sfc_dlt[lon_idx]; // [K]
    tpt_gnd[lon_idx]=tpt_sfc[lon_idx]; // [K]
    tpt_gnd_dlt[lon_idx]=tpt_sfc_dlt[lon_idx]; // [K]
    
    // Use iterated tpt_vgt to initialize coefficients for fluxes at tpt_gnd
    if(vgt[lon_idx]){
      flx_LW_net_cff_a[lon_idx]= // [W m-2] a in FLWnet = a + b*Tg^4 Bon96 p. 45
	-msv_sfc[lon_idx]*(1.0-msv_vgt[lon_idx])*flx_LW_dwn_sfc[lon_idx] 
	-msv_sfc[lon_idx]*msv_vgt[lon_idx]*cst_Stefan_Boltzmann*std::pow(tpt_vgt[lon_idx],PRC_CMP(4.0));
      flx_LW_net_cff_b[lon_idx]=msv_sfc[lon_idx]*cst_Stefan_Boltzmann; // [W m-2 K-4] b in FLWnet = a + b*Tg^4 Bon96 p. 45
      flx_sns_atm_cff_b[lon_idx]=dns_mdp[lon_idx]*spc_heat_dry_air/rss_aer_heat(lon_idx,idx_gnd); // [W m-2 K-1] b in SH = a + b*Ts Bon96 p. 55, 69, 70
      flx_sns_atm_cff_a[lon_idx]=-tpt_ash[lon_idx]*flx_sns_atm_cff_b[lon_idx]; // [W m-2] a in SH = a + b*Ts Bon96 p. 55, 69, 70
      flx_ltn_evp_cff_b[lon_idx]=dns_mdp[lon_idx]*spc_heat_dry_air/(cst_psych[lon_idx]*(rss_aer_vpr(lon_idx,idx_gnd)+rss_sfc_vpr[lon_idx])); // [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
      flx_ltn_evp_cff_a[lon_idx]=-ppr_H2O_cnp[lon_idx]*flx_ltn_evp_cff_b[lon_idx]; // [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
      flx_sns_gnd_cff_b[lon_idx]=2.0*cnd_trm_soi[lon_idx]/lvl_dlt_snw[lon_idx]; // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
      flx_sns_gnd_cff_a[lon_idx]=-flx_sns_gnd_cff_b[lon_idx]*tpt_soi[lon_idx]; // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
    } // endif vegetated
    
    if(vgt[lon_idx]){
      // Initialize accuracy and counter
      eps_crr=eps_max+1.0; // [frc] Current relative accuracy
      itr_idx=1; // Counting index
      if(dbg_lvl == dbg_crr){
	(void)std::fprintf(stderr,"Ground fluxes from %s():\n",sbr_nm.c_str());
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
			   "itr","mno_lng","mno_stb","wnd_frc","  Ts   ","LW(msn)"," H(atm)"," H(soi)","   L   ","nrg_bdg","  eps  ");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
			   "    ","   m   ","  frc  "," m s-1 ","   K   "," W m-2 "," W m-2 "," W m-2 "," W m-2 "," W m-2 ","  frc  ");
      } // end if dbg
      // This iteration loop solves for tpt_grd if vegetated, else do nothing
      while(eps_crr > eps_max){
	
	// Save old surface temperature for convergence diagnostic
	tpt_gnd_old=tpt_gnd[lon_idx]; // [K]
	
	// Saturation vapor pressure of water at surface temperature
	tpt_bnd_cls=tpt_bnd_cls_get(tpt_gnd[lon_idx]); // [C]
	if(tpt_bnd_cls > 0.0){ 
	  svp_H2O_gnd[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	  dsvpdt_H2O_gnd[lon_idx]=dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
	}else{
	  svp_H2O_gnd[lon_idx]=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	  dsvpdt_H2O_gnd[lon_idx]=dsvpdt_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
	} // endif frozen
	
	// Evaluate fluxes for current Tg
	flx_sns_gnd(lon_idx,idx_gnd)=flx_sns_gnd_cff_a[lon_idx]+flx_sns_gnd_cff_b[lon_idx]*tpt_gnd[lon_idx]; // [W m-2] Sensible heat flux to soil
	flx_sns_atm(lon_idx,idx_gnd)=flx_sns_atm_cff_a[lon_idx]+flx_sns_atm_cff_b[lon_idx]*tpt_gnd[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
	flx_ltn_trn(lon_idx,idx_gnd)=0.0; // [W m-2] Transpiration flux to atmosphere
	flx_ltn_evp(lon_idx,idx_gnd)=flx_ltn_evp_cff_a[lon_idx]+flx_ltn_evp_cff_b[lon_idx]*svp_H2O_gnd[lon_idx]; // [W m-2] Evaporation flux to atmosphere
	flx_LW_net(lon_idx,idx_gnd)=flx_LW_net_cff_a[lon_idx]+flx_LW_net_cff_b[lon_idx]*std::pow(tpt_gnd[lon_idx],PRC_CMP(4.0)); // [W m-2] Net longwave flux to atmosphere
	nrg_bdg[lon_idx]= // [W m-2] Total energy budget at surface
	  +flx_SW_net(lon_idx,idx_gnd)
	  -flx_LW_net(lon_idx,idx_gnd)
	  -flx_ltn_evp(lon_idx,idx_gnd)
	  -flx_ltn_trn(lon_idx,idx_gnd)
	  -flx_sns_atm(lon_idx,idx_gnd)
	  -flx_sns_gnd(lon_idx,idx_gnd);
	nrg_bdg_dlt= // [W m-2 K-1] Temperature derivative of surface energy budget
	  -4.0*flx_LW_net_cff_b[lon_idx]*std::pow(tpt_gnd[lon_idx],PRC_CMP(3.0))
	  -flx_ltn_evp_cff_b[lon_idx]*dsvpdt_H2O_gnd[lon_idx]
	  -flx_sns_atm_cff_b[lon_idx]
	  -flx_sns_gnd_cff_b[lon_idx];
	tpt_gnd_dlt[lon_idx]=-nrg_bdg[lon_idx]/nrg_bdg_dlt; // [K] Newton-Raphson temperature adjustment
	
	eps_crr=PRC_CMP_ABS((tpt_gnd[lon_idx]-tpt_gnd_old)/tpt_gnd[lon_idx]); // Relative convergence
	if(dbg_lvl == dbg_crr){
	  // NB: Fluxes reported here have not been adjusted for snowmelt yet
	  // fxm: LW flux here is ground flux 
	  // 20210618 Emissivity factor was incorrectly applied to Fdwn in below line for ~25 years
	  //flx_LW_upw_sfc[lon_idx]=flx_LW_net(lon_idx,lvl_idx)+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
	  flx_LW_upw_sfc[lon_idx]=flx_LW_net(lon_idx,lvl_idx)+flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
	  flx_LW_rfl_sfc[lon_idx]=(1.0-msv_sfc[lon_idx])*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave reflected flux at surface
	  flx_LW_msn_sfc[lon_idx]=flx_LW_upw_sfc[lon_idx]-flx_LW_rfl_sfc[lon_idx]; // [W m-2] Longwave emitted flux at surface
	  flx_sns_gnd_ttl[lon_idx]=flx_sns_gnd(lon_idx,lvl_idx); // [W m-2] Sensible heat flux to soil
	  (void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
			itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],tpt_sfc[lon_idx],flx_LW_msn_sfc[lon_idx],flx_sns_atm_tmp,flx_sns_gnd_ttl[lon_idx],flx_vpr_tmp,nrg_bdg[lon_idx],eps_crr);
	} // end if dbg
	if(itr_idx > itr_max_gnd){
	  std::cerr << "Final: tpt_gnd = " << tpt_gnd[lon_idx] << " K, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	  wrn_prn(prg_nm,sbr_nm,"Ground temperature not converging, breaking loop...");
	  break;
	} // endif

	tpt_gnd[lon_idx]+=tpt_gnd_dlt[lon_idx]; // [K] Ground temperature fxm: Bon96 does not update the temperature on the final iteration
	itr_idx++;
      } // end loop over itr
      
    } // endif vegetated

    // Adjust fluxes before adjusting temperatures
    if(vgt[lon_idx]){
      flx_sns_gnd(lon_idx,idx_gnd)+=flx_sns_gnd_cff_b[lon_idx]*tpt_gnd_dlt[lon_idx]; // [W m-2] Sensible heat flux to soil
      flx_sns_atm(lon_idx,idx_gnd)+=flx_sns_atm_cff_b[lon_idx]*tpt_gnd_dlt[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn_evp(lon_idx,idx_gnd)+=flx_ltn_evp_cff_b[lon_idx]*dsvpdt_H2O_gnd[lon_idx]*tpt_gnd_dlt[lon_idx]; // [W m-2] Evaporation flux to atmosphere
      flx_LW_net(lon_idx,idx_gnd)+=4.0*flx_LW_net_cff_b[lon_idx]*std::pow(tpt_gnd[lon_idx],PRC_CMP(3.0))*tpt_gnd_dlt[lon_idx]; // [W m-2] Net longwave flux to atmosphere
      tpt_gnd[lon_idx]+=tpt_gnd_dlt[lon_idx]; // [K] Ground temperature
    } // endif vegetated

    /* If snow on ground and tpt_gnd > tpt_frz_pnt: reset tpt_gnd = tpt_frz_pnt and reevaluate ground fluxes
       Condition that snw_hgt_lqd > 0.5 m prevents spurious fluxes
       Procedure described on Bon96 p. 65 */
    if(snw_hgt_lqd[lon_idx] > 0.5 && tpt_gnd[lon_idx] > tpt_frz_pnt){
      tpt_gnd[lon_idx]=tpt_frz_pnt; // [K] 
      flx_LW_net(lon_idx,idx_gnd)=flx_LW_net_cff_a[lon_idx]+flx_LW_net_cff_b[lon_idx]*std::pow(tpt_frz_pnt,4.0); // [W m-2] Net longwave flux to atmosphere
      flx_sns_atm(lon_idx,idx_gnd)=flx_sns_atm_cff_a[lon_idx]+flx_sns_atm_cff_b[lon_idx]*tpt_frz_pnt; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn_evp(lon_idx,idx_gnd)=flx_ltn_evp_cff_a[lon_idx]+flx_ltn_evp_cff_b[lon_idx]*svp_H2O_frz; // [W m-2] Evaporation flux to atmosphere
      flx_ltn_trn(lon_idx,idx_gnd)=0.0; // [W m-2] Transpiration flux to atmosphere
      flx_sns_gnd(lon_idx,idx_gnd)=flx_SW_net(lon_idx,idx_gnd)-(flx_LW_net(lon_idx,idx_gnd)+flx_sns_atm(lon_idx,idx_gnd)+flx_ltn_evp(lon_idx,idx_gnd)+flx_ltn_trn(lon_idx,idx_gnd)); // [W m-2] Sensible heat flux to soil
    } // endif
    /* Energy imbalance used to melt snow and warm ground
       flx_snw_mlt must be less than both flx_sns_gnd and flux required to completely melt snow in one timestep */
    if (snw_hgt_lqd[lon_idx] > 0.0 && tpt_gnd[lon_idx] >= tpt_frz_pnt){
      flx_snw_mlt[lon_idx]=min_cpv(snw_hgt_lqd[lon_idx]*ltn_heat_fsn_H2O_std/tm_dlt,max_cpv(flx_sns_gnd(lon_idx,idx_gnd),0.0)); // [W m-2] Snow melt heat flux Bon96 p. 65
    }else{
      flx_snw_mlt[lon_idx]=0.0; // [W m-2] Snow melt heat flux Bon96 p. 65
    } // endif

    // Final fluxes, radiative temperature, wind stresses
    flx_LW_net_ttl[lon_idx]=flx_LW_net(lon_idx,idx_vgt)+flx_LW_net(lon_idx,idx_gnd); // [W m-2] Total net longwave flux to atmosphere
    flx_sns_atm_ttl[lon_idx]=flx_sns_atm(lon_idx,idx_vgt)+flx_sns_atm(lon_idx,idx_gnd); // [W m-2] Total sensible heat flux to atmosphere
    flx_ltn_evp_cnp_atm[lon_idx]=flx_ltn_evp(lon_idx,idx_vgt); // [W m-2] Canopy evaporation flux to atmosphere
    flx_ltn_trn_cnp_atm[lon_idx]=flx_ltn_trn(lon_idx,idx_vgt); // [W m-2] Canopy transpiration flux to atmosphere
    flx_ltn_evp_gnd_atm[lon_idx]=flx_ltn_evp(lon_idx,idx_gnd); // [W m-2] Ground evaporation flux to atmosphere
    flx_sns_atm_cnp[lon_idx]=flx_sns_gnd(lon_idx,idx_vgt); // [W m-2] Canopy heat storage (currently defined = 0.0)
    flx_sns_gnd_ttl[lon_idx]=flx_sns_gnd(lon_idx,idx_gnd)-flx_snw_mlt[lon_idx]; // [W m-2] Sensible heat flux to soil
    // 20210618 Emissivity factor was incorrectly applied to Fdwn in below line for ~25 years
    //flx_LW_upw_sfc[lon_idx]=flx_LW_net_ttl[lon_idx]+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 fxm: implement correct formulae for vegetated surfaces as well
    flx_LW_upw_sfc[lon_idx]=flx_LW_net_ttl[lon_idx]+flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 fxm: implement correct formulae for vegetated surfaces as well
    flx_LW_rfl_sfc[lon_idx]=(1.0-msv_sfc[lon_idx])*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave reflected flux at surface
    flx_LW_msn_sfc[lon_idx]=flx_LW_upw_sfc[lon_idx]-flx_LW_rfl_sfc[lon_idx]; // [W m-2] Longwave emitted flux at surface
    tpt_ffc[lon_idx]=std::pow(flx_LW_upw_sfc[lon_idx]/cst_Stefan_Boltzmann,0.25); // [K] Radiative effective temperature
    wnd_str_znl[lon_idx]=-dns_mdp[lon_idx]*cff_xch_mmn[lon_idx]*wnd_mdp_bnd[lon_idx]*wnd_znl_mdp[lon_idx]; // [kg m-1 s-2] Zonal wind stress Bon96 p. 54
    wnd_str_mrd[lon_idx]=-dns_mdp[lon_idx]*cff_xch_mmn[lon_idx]*wnd_mdp_bnd[lon_idx]*wnd_mrd_mdp[lon_idx]; // [kg m-1 s-2] Meridional wind stress Bon96 p. 54

    // "Screen" temperature at z=zpd+rgh_heat+2
    // fxm: 19990519: Screen temperature results look strange, since temperature 
    // does not change monotonically from Tsfc to Tgcm, could be a problem? 
    tpt_ash_p2m[lon_idx]= // [K] "Screen" temperature at z=zpd+rgh_heat+2m Bon96 p. 55 
      +tpt_ash[lon_idx] 
      -flx_sns_atm_ttl[lon_idx] 
      /(dns_mdp[lon_idx]*spc_heat_dry_air*wnd_frc[lon_idx]) 
      *cst_von_krm_rcp
      *std::log((2.0+rgh_heat[lon_idx])/rgh_heat[lon_idx]);

    // Non-LSM diagnostics
    flx_ltn[lon_idx]= // [W m-2] Latent heat flux to atmosphere
      +flx_ltn_evp_cnp_atm[lon_idx]+flx_ltn_trn_cnp_atm[lon_idx]+flx_ltn_evp_gnd_atm[lon_idx];
    flx_q_H2O[lon_idx]=flx_vpr_tmp/ltn_heat_trn[lon_idx]; // [kg m-2 s-1] Moisture flux to atmosphere fxm: units? maybe ltn_heat_trn not necessary
    rss_aer_heat_sfc[lon_idx]=rss_aer_heat(lon_idx,idx_vgt); // [s m-1] Aerodynamic resistance to heat transfer
    rss_aer_mmn_sfc[lon_idx]=rss_aer_mmn(lon_idx,idx_vgt); // [s m-1] Aerodynamic resistance to momentum transfer
    rss_aer_vpr_sfc[lon_idx]=rss_aer_vpr(lon_idx,idx_vgt); // [s m-1] Aerodynamic resistance to vapor transfer
      
    if(dbg_lvl == dbg_crr){
      prc_cmp sum_LHS;
      prc_cmp sum_RHS;
      sum_LHS=flx_SW_net_gnd[lon_idx]+flx_SW_net_vgt[lon_idx]+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2]
      sum_RHS=flx_LW_msn_sfc[lon_idx]+flx_sns_atm_ttl[lon_idx]+flx_sns_gnd_ttl[lon_idx]+flx_ltn_evp_gnd_atm[lon_idx]+flx_snw_mlt[lon_idx]; // [W m-2]
      (void)std::fprintf(stderr," LHS     =  SW(abs) +  LW(abs) = SW(abs)  + a*LW(dwn)\n");
      (void)std::fprintf(stderr,"%8.3f = %8.3f + %8.3f\n",sum_LHS,flx_SW_net_gnd[lon_idx]+flx_SW_net_vgt[lon_idx],msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]);
      (void)std::fprintf(stderr," RHS     =  LW(msn) +  H(atm)  +  H(soi)  +  L(atm)  +  M(snw)\n");
      (void)std::fprintf(stderr,"%8.3f = %8.3f + %8.3f + %8.3f + %8.3f + %8.3f\n",sum_RHS,flx_LW_msn_sfc[lon_idx],flx_sns_atm_ttl[lon_idx],flx_sns_gnd_ttl[lon_idx],flx_ltn_evp_gnd_atm[lon_idx],flx_snw_mlt[lon_idx]);
      (void)std::fprintf(stderr," LW(msn) +  LW(rfl) =  LW(upw)\n%8.3f + %8.3f = %8.3f\n",flx_LW_msn_sfc[lon_idx],flx_LW_rfl_sfc[lon_idx],flx_LW_upw_sfc[lon_idx]);
    } // end if dbg

  } // end loop over lon
  
  return rcd;
} // end flx_sfc_lnd()

int // O [rcd] Return success code
blm_mbl
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
 const prc_cmp *flx_SW_net, // I [W m-2] Solar flux absorbed by ground
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *lvl_dlt, // I [m] Soil layer thickness
 const prc_cmp *msv_gnd, // I [frc] Bare ground emissivity
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *mno_lng, // I/O [m] Monin-Obukhov length
 prc_cmp *tpt_gnd, // I/O [K] Ground temperature
 prc_cmp *wnd_frc) // O [m s-1] Surface friction velocity
{
  /* Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
     compute the boundary layer exchange properties
     Routine is optimized for dust source regions: dry, bare, uncovered land
     Theory and algorithms: Bonan (1996) CCM:lsm/surtem() */
  // Dependencies: <phys_cst.hh>,<tdy.hh>
  /* Notes on algorithm:
     Suffix mdp quantity evaluated at height hgt_mdp
     Suffix gnd quantity evaluated at ground */
  /* Relationships between various temperatures: Bon96 p. 55 and Fig. 16 p. 57:
     T1 = Skin temperature = Soil temperature of 1st layer (top 10 cm)
     Tg = "Ground" temperature = Air temperature at z = rgh_heat
     Ts = "Surface" temperature = Air temperature at z=zpd+rgh_heat
     Ta = "Aerodynamic" temperature = Air temperature at z=zpd+rgh_mmn
     Te = Emitting temperature = (FLWup/sigma)**0.25
     Tgcm = Ambient temperature = Air temperature at z=hgt_mdp
     For bare ground, Ts = Tg = Air temperature at z=rgh_heat (Bon96 p. 55)
     For bare ground, Tg = potential temperature defined relative to local surface pressure */
  
  // Output
  int rcd(0); // O [rcd] Return success code

  // Local
  const prc_cmp eps_dbz(1.0e-6); // [frc] Prevents division by zero
  const prc_cmp wnd_min_mbl(1.0); // [m s-1] Minimum windspeed used for mobilization
  const long itr_max_gnd(12); // Maximum number of iterations for tpt_gnd loop
  const long sgn_chg_ctr_max(4); // Maximum number of sign changes in stability parameter

  std::valarray<prc_cmp> cff_xch_mmn(lon_nbr); // [frc] Exchange coefficient for momentum transfer
  prc_cmp cnd_heat_sfc_mdp; // [m s-1] Sensible heat conductance surface air to midlayer air
  prc_cmp cnd_vpr_gnd_sfc; // [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
  prc_cmp cnd_vpr_sfc_mdp; // [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57

  prc_cmp cnd_vpr_ttl; // [m s-1] Sum of conductances
  std::valarray<prc_cmp> cst_psych(lon_nbr); // [Pa K-1] Psychrometric constant
  std::valarray<prc_cmp> dsvpdt_H2O_gnd(lon_nbr); // [Pa K-1] Derivative of saturation vapor pressure over planar condensed water, ground
  std::valarray<prc_cmp> flx_LW_net(lon_nbr); // [W m-2] Net longwave flux to atmosphere
  std::valarray<prc_cmp> flx_LW_net_cff_a(lon_nbr); // [W m-2] a in FLWnet = a + b*T^4 Bon96 p. 45
  std::valarray<prc_cmp> flx_LW_net_cff_b(lon_nbr); // [W m-2 K-4] b in FLWnet = a + b*T^4 Bon96 p. 45
  std::valarray<prc_cmp> flx_LW_upw_sfc(lon_nbr); // [W m-2] Longwave upwelling (emission+reflection) flux at surface
  std::valarray<prc_cmp> flx_LW_msn_sfc(lon_nbr); // [W m-2] Longwave emitted flux at surface
  std::valarray<prc_cmp> flx_LW_rfl_sfc(lon_nbr); // [W m-2] Longwave reflected flux at surface
  std::valarray<prc_cmp> flx_ltn_evp(lon_nbr); // [W m-2] Evaporation flux to atmosphere
  std::valarray<prc_cmp> flx_ltn_evp_cff_a(lon_nbr); // [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> flx_ltn_evp_cff_b(lon_nbr); // [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
  prc_cmp flx_ltn_fct; // Factor in vapor flux calculations
  std::valarray<prc_cmp> flx_sns_atm(lon_nbr); // [W m-2] Sensible heat flux to atmosphere
  std::valarray<prc_cmp> flx_sns_atm_cff_a(lon_nbr); // [W m-2] a in SH = a + b*Ts
  std::valarray<prc_cmp> flx_sns_atm_cff_b(lon_nbr); // [W m-2 K-1] b in SH = a + b*Ts
  prc_cmp flx_sns_atm_fct; // Factor in sensible heat computation
  prc_cmp flx_sns_atm_tmp; // [W m-2] Temporary sensible heat flux
  prc_cmp flx_sns_atm_vrt_tmp; // [W m-2] Temporary virtual sensible heat flux Bon96 p. 49
  std::valarray<prc_cmp> flx_sns_gnd(lon_nbr); // [W m-2] Sensible heat flux to soil
  std::valarray<prc_cmp> flx_sns_gnd_cff_a(lon_nbr); // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
  std::valarray<prc_cmp> flx_sns_gnd_cff_b(lon_nbr); // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
  prc_cmp flx_vpr_tmp; // [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
  std::valarray<prc_cmp> ltn_heat_trn(lon_nbr); // [J kg-1] Latent heat of sublimation or evaporation
  prc_cmp mno_dnm; // Denominator of Monin-Obukhov length Bon96 p. 49
  std::valarray<prc_cmp> mno_stb_crc_heat(lon_nbr); // [frc] Monin-Obukhov stability correction heat
  prc_cmp mno_stb_crc_heat_crr; // Undamped correction factor heat [frc]
  std::valarray<prc_cmp> mno_stb_crc_mmn(lon_nbr); // [frc] Monin-Obukhov stability correction momentum
  prc_cmp mno_stb_crc_mmn_crr; // Undamped correction factor momentum [frc]
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  prc_cmp mno_stb_crc_tmp2; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp3; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp4; // Term in stability correction computation
  prc_cmp mno_stb_crc_tmp5; // Term in stability correction computation
  std::valarray<prc_cmp> mno_stb_prm(lon_nbr); // [frc] Monin-Obukhov stability parameter 
  std::valarray<prc_cmp> mno_stb_prm_old(lon_nbr); // [frc] Monin Obukhov stability parameter old
  std::valarray<prc_cmp> msv_sfc(lon_nbr); // [frc] Surface (bare ground+snow) emissivity
  std::valarray<prc_cmp> nrg_bdg(lon_nbr); // [W m-2] Surface energy budget
  prc_cmp nrg_bdg_dlt; // [W m-2 K-1] Temperature derivative of surface energy budget
  std::valarray<prc_cmp> ppr_H2O_cnp(lon_nbr); // [Pa] Canopy vapor pressure of H2O
  std::valarray<prc_cmp> ppr_H2O_cnp_cff_a(lon_nbr); // [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> ppr_H2O_cnp_cff_b(lon_nbr); // [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55
  std::valarray<prc_cmp> ppr_H2O_mdp(lon_nbr); // [Pa] Ambient vapor pressure of H2O
  std::valarray<prc_cmp> rgh_heat(lon_nbr); // [m] Roughness length heat
  std::valarray<prc_cmp> rss_aer_fct(lon_nbr); // [s m-1] Term in resistance calculation
  std::valarray<prc_cmp> rss_aer_heat(lon_nbr); // [s m-1] Aerodynamic resistance to heat transfer
  std::valarray<prc_cmp> rss_aer_heat_fct(lon_nbr); // [frc] Term in resistance calculation
  std::valarray<prc_cmp> rss_aer_mmn(lon_nbr); // [s m-1] Aerodynamic resistance to momentum transfer
  std::valarray<prc_cmp> rss_aer_mmn_fct(lon_nbr); // [frc] Term in resistance calculation
  std::valarray<prc_cmp> rss_aer_vpr(lon_nbr); // [s m-1] Aerodynamic resistance to vapor transfer
  std::valarray<prc_cmp> rss_sfc_vpr(lon_nbr); // [s m-1] Surface resistance to vapor transfer
  std::valarray<prc_cmp> svp_H2O_gnd(lon_nbr); // [Pa] Saturation vapor pressure over planar condensed water at ground
  prc_cmp tpt_bnd_cls; // [C] Temperature bounded celsius
  std::valarray<prc_cmp> tpt_gnd_dlt(lon_nbr); // [K] Ground temperature adjustment
  prc_cmp tpt_gnd_old; // [K] Previous ground temperature
  std::valarray<prc_cmp> tpt_vrt(lon_nbr); // [K] Virtual temperature
  std::valarray<prc_cmp> wnd_mdp(lon_nbr); // [m s-1] Surface layer mean wind speed
  std::valarray<prc_cmp> wnd_mdp_bnd(lon_nbr); // [m s-1] Surface layer mean wind speed bounded
  long lon_idx; // Counting index for lon
  std::valarray<long> sgn_chg_ctr(lon_nbr); // Number of sign changes in stability parameter
  
  // Main Code
  const std::string sbr_nm("blm_mbl"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  
  // Initialize variables that normally would be available in LSM
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    ppr_H2O_mdp[lon_idx]=q_H2O_vpr[lon_idx]*prs_mdp[lon_idx]/(eps_H2O+one_mns_eps_H2O*q_H2O_vpr[lon_idx]); // [Pa] Ambient vapor pressure of H2O
    
    mno_stb_prm[lon_idx]=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    
    msv_sfc[lon_idx]=msv_gnd[lon_idx]; // [frc] Surface (bare ground+snow) emissivity

    rgh_heat[lon_idx]=rgh_mmn[lon_idx]; // [m] Roughness length heat
    
    tpt_vrt[lon_idx]=tpt_mdp[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Virtual temperature
    
  }  // end loop over lon
  
  // Initialize variables which are independent of stability iteration
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    // Zero temperature adjustments from last timestep
    tpt_gnd_dlt[lon_idx]=0.0; // [K] Change in ground temperature
    
    // Latent heat of water transformation
    tpt_bnd_cls=tpt_bnd_cls_get(tpt_mdp[lon_idx]); // [C]
    if(tpt_bnd_cls > 0.0){ 
      ltn_heat_trn[lon_idx]=ltn_heat_vpr_H2O_std; // [J kg-1]
    }else{ 
      ltn_heat_trn[lon_idx]=ltn_heat_sbl_H2O_std; // [J kg-1]
    } // endif
    // Psychrometric constant for transformations of water vapor at surface
    cst_psych[lon_idx]=spc_heat_dry_air*prs_mdp[lon_idx]/(eps_H2O*ltn_heat_trn[lon_idx]); // [Pa K-1] 
    
  } // end loop over lon
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    // Midlayer wind speeds
    wnd_mdp[lon_idx]= // [m s-1] Surface layer mean wind speed
      std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	   wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
    wnd_mdp_bnd[lon_idx]=max_cpv(wnd_mdp[lon_idx],wnd_min_mbl); // [m s-1] Surface layer mean wind speed bounded
    
    // Miscellaneous
    sgn_chg_ctr[lon_idx]=0; // Number of sign changes in stability parameter
    mno_stb_prm_old[lon_idx]=0.0; // [frc] Monin Obukhov stability parameter old
    
    // Stability-independent components of aerodynamic resistance calculations
    rss_aer_fct[lon_idx]=1.0/(cst_von_krm*cst_von_krm*wnd_mdp_bnd[lon_idx]); // [s m-1]
    rss_aer_mmn_fct[lon_idx]=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_mmn[lon_idx]); // [frc]
    rss_aer_heat_fct[lon_idx]=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_heat[lon_idx]); // [frc]
  } // end loop over lon
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    
    // F(LW net) = -msv*Fdwn + msv*sigma*Tg^4 = a + b*Tg^4 Bon96 p. 45
    // F(sns heat dwn into soil) = 2*k*(Tg-T1)/dz = a + b*Tg Bon96 p. 64
    flx_LW_net_cff_a[lon_idx]=-msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] a in FLWnet = a + b*Tg^4 Bon96 p. 45
    flx_LW_net_cff_b[lon_idx]=msv_sfc[lon_idx]*cst_Stefan_Boltzmann; // [W m-2 K-4] b in FLWnet = a + b*Tg^4 Bon96 p. 45
    
    flx_sns_gnd_cff_b[lon_idx]=2.0*cnd_trm_soi[lon_idx]/lvl_dlt[lon_idx]; // [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
    flx_sns_gnd_cff_a[lon_idx]=-flx_sns_gnd_cff_b[lon_idx]*tpt_soi[lon_idx]; // [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
  }  // end loop over lon
  
  // Iteration loop
  const prc_cmp eps_max(1.0e-5); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  long itr_idx; // Counting index
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    if(dbg_lvl == dbg_sbr){
      (void)std::fprintf(stderr,"Ground fluxes from %s():\n",sbr_nm.c_str());
      (void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
			 "itr","mno_lng","mno_stb","wnd_frc","  Ts   "," LW(up)"," H(atm)"," H(soi)","   L   ","nrg_bdg","  eps  ");
      (void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
			 "    ","   m   ","  frc  "," m s-1 ","   K   "," W m-2 "," W m-2 "," W m-2 "," W m-2 "," W m-2 ","  frc  ");
    } // end if dbg
    
    // This iteration loop solves for tpt_vgt if vegetated, else tpt_gnd
    while(eps_crr > eps_max){
      
      // Save old ground temperature for convergence diagnostic
      if(itr_idx == 1) tpt_gnd_old=tpt_gnd[lon_idx]; else tpt_gnd_old=tpt_gnd[lon_idx]; // [K]
      
      // Stability functions computed as in Bon96 p. 52
      if(mno_stb_prm[lon_idx] < 0.0){
	sml_fnc_mmn_uns_rcp=std::pow(PRC_CMP(1.0)-PRC_CMP(16.0)*mno_stb_prm[lon_idx],PRC_CMP(0.25));
	mno_stb_crc_tmp2=std::log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_tmp3=std::log((1.0+sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_mmn_crr=2.0*mno_stb_crc_tmp3+mno_stb_crc_tmp2-2.0*std::atan(sml_fnc_mmn_uns_rcp)+1.5707963; // [frc]
	mno_stb_crc_heat_crr=2.0*mno_stb_crc_tmp2; // [frc]
      }else{ // not stable
	mno_stb_crc_mmn_crr=-5.0*mno_stb_prm[lon_idx]; // [frc]
	mno_stb_crc_heat_crr=mno_stb_crc_mmn[lon_idx]; // [frc]
      } // not stable
      
      // Filter stability corrections to reduce numerical ping-pong
      if(itr_idx == 1){
	mno_stb_crc_mmn[lon_idx]=mno_stb_crc_mmn_crr; // [frc]
	mno_stb_crc_heat[lon_idx]=mno_stb_crc_heat_crr; // [frc]
      }else{
	mno_stb_crc_mmn[lon_idx]=0.5*(mno_stb_crc_mmn_crr+mno_stb_crc_mmn[lon_idx]); // [frc]
	mno_stb_crc_heat[lon_idx]=0.5*(mno_stb_crc_heat_crr+mno_stb_crc_heat[lon_idx]); // [frc]
      } // endif first iteration
      
      // Aerodynamic resistance between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
      mno_stb_crc_tmp4=rss_aer_mmn_fct[lon_idx]-mno_stb_crc_mmn[lon_idx]; // [frc]
      mno_stb_crc_tmp5=rss_aer_heat_fct[lon_idx]-mno_stb_crc_heat[lon_idx]; // [frc]
      rss_aer_mmn[lon_idx]=max_cpv(rss_aer_fct[lon_idx]*mno_stb_crc_tmp4*mno_stb_crc_tmp4,1.0); // [s m-1]
      rss_aer_heat[lon_idx]=max_cpv(rss_aer_fct[lon_idx]*mno_stb_crc_tmp4*mno_stb_crc_tmp5,1.0); // [s m-1]
      // Resistances are equal because rgh_heat = rgh_vpr
      rss_aer_vpr[lon_idx]=rss_aer_heat[lon_idx]; // [s m-1]
      
      // Exchange coefficients between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
      // Exchange coefficients are dimensionless, inversely proportional to wind speed 
      cff_xch_mmn[lon_idx]=1.0/(rss_aer_mmn[lon_idx]*wnd_mdp_bnd[lon_idx]); // [frc]
      // Friction velocity
      wnd_frc[lon_idx]=wnd_mdp_bnd[lon_idx]*std::sqrt(cff_xch_mmn[lon_idx]); // [m s-1]
      
      // Saturation vapor pressure of water at ground temperature
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_gnd[lon_idx]); // [C]
      if(tpt_bnd_cls > 0.0){ 
	svp_H2O_gnd[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	dsvpdt_H2O_gnd[lon_idx]=dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
      }else{
	svp_H2O_gnd[lon_idx]=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	dsvpdt_H2O_gnd[lon_idx]=dsvpdt_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa K-1]
      } // endif frozen
      
      // Surface resistance to vapor transfer
      if(itr_idx == 1){ // NB: Bon96 only computes this on the first iteration
	rss_sfc_vpr[lon_idx]= // [s m-1] Surface resistance to vapor transfer Bon96 p. 56, Fgr. 17 p. 57
	    rss_aer_vpr[lon_idx]
	    *(1.0-trn_fsh_vpr_soi_atm[lon_idx])/
	    max_cpv(trn_fsh_vpr_soi_atm[lon_idx],eps_dbz);

	// Set minimum rss_sfc_vpr so that in bare ground case es = svp(tg)
	rss_sfc_vpr[lon_idx]=max_cpv(rss_sfc_vpr[lon_idx],eps_dbz);
	
      } // endif first iteration
      
      // Heat conductances from ground to ambient air
      // Conductances are dimensional, exact inverses of resistances
      // Coefficients for sensible heat flux SH = a + b*Ts
      cnd_heat_sfc_mdp=1.0/rss_aer_heat[lon_idx]; // [m s-1] Sensible heat conductance surface air to midlayer air Bon96 p. 60, Fig. 16 p. 57
      flx_sns_atm_fct=dns_mdp[lon_idx]*spc_heat_dry_air*cnd_heat_sfc_mdp; // Bon96 p. 55
      flx_sns_atm_cff_a[lon_idx]=-tpt_ptn_mdp[lon_idx]*flx_sns_atm_fct; // [W m-2] a in SH = a + b*Ts Bon96 p. 55, 69, 70
      flx_sns_atm_cff_b[lon_idx]=flx_sns_atm_fct; // [W m-2 K-1] b in SH = a + b*Ts Bon96 p. 55, 69, 70
      
      // Vapor conductances from ground to ambient air
      flx_ltn_fct=dns_mdp[lon_idx]*spc_heat_dry_air/cst_psych[lon_idx]; // 
      cnd_vpr_sfc_mdp=1.0/rss_aer_vpr[lon_idx]; // [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57
      cnd_vpr_gnd_sfc=1.0/rss_sfc_vpr[lon_idx]; // [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
      cnd_vpr_ttl=cnd_vpr_sfc_mdp+cnd_vpr_gnd_sfc; // [m s-1] Sum of conductances
      
      // Coefficients for canopy vapor pressure e(cnp) = a + b*e(Ts)
      ppr_H2O_cnp_cff_a[lon_idx]=ppr_H2O_mdp[lon_idx]*cnd_vpr_sfc_mdp/cnd_vpr_ttl; // [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
      ppr_H2O_cnp_cff_b[lon_idx]=cnd_vpr_gnd_sfc/cnd_vpr_ttl; // [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55

      // Coefficients for evaporation heat flux LHE = a + b*e(Ts) 
      flx_ltn_evp_cff_a[lon_idx]=-flx_ltn_fct*(ppr_H2O_mdp[lon_idx]-ppr_H2O_cnp_cff_a[lon_idx])*cnd_vpr_sfc_mdp; // [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
      flx_ltn_evp_cff_b[lon_idx]=flx_ltn_fct*ppr_H2O_cnp_cff_b[lon_idx]*cnd_vpr_sfc_mdp; // [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55

      // Evaluate fluxes for current tpt_gnd
      flx_sns_gnd[lon_idx]=flx_sns_gnd_cff_a[lon_idx]+flx_sns_gnd_cff_b[lon_idx]*tpt_gnd[lon_idx]; // [W m-2] Sensible heat flux to soil
      flx_sns_atm[lon_idx]=flx_sns_atm_cff_a[lon_idx]+flx_sns_atm_cff_b[lon_idx]*tpt_gnd[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn_evp[lon_idx]=flx_ltn_evp_cff_a[lon_idx]+flx_ltn_evp_cff_b[lon_idx]*svp_H2O_gnd[lon_idx]; // [W m-2] Evaporation flux to atmosphere
      flx_LW_net[lon_idx]=flx_LW_net_cff_a[lon_idx]+flx_LW_net_cff_b[lon_idx]*std::pow(tpt_gnd[lon_idx],PRC_CMP(4.0)); // [W m-2] Net longwave flux to atmosphere
      nrg_bdg[lon_idx]= // [W m-2] Total energy budget at surface
	+flx_SW_net[lon_idx]
	-flx_LW_net[lon_idx]
	-flx_ltn_evp[lon_idx]
	-flx_sns_atm[lon_idx]
	-flx_sns_gnd[lon_idx];
      nrg_bdg_dlt= // [W m-2 K-1] Temperature derivative of surface energy budget
	-4.0*flx_LW_net_cff_b[lon_idx]*std::pow(tpt_gnd[lon_idx],PRC_CMP(3.0))
	-flx_ltn_evp_cff_b[lon_idx]*dsvpdt_H2O_gnd[lon_idx]
	-flx_sns_atm_cff_b[lon_idx]
	-flx_sns_gnd_cff_b[lon_idx];
      tpt_gnd_dlt[lon_idx]=-nrg_bdg[lon_idx]/nrg_bdg_dlt; // [K] Newton-Raphson temperature adjustment
      
      /* NB:
	 Here the fluxes could be incremented by their "b" coefficient
	 times the temperature differential. flx_sfc_lnd() does, in fact, do this.
	 However, blm_mbl() is optimized for speed so we omit the flux update
	 Instead, we proceed directly to temperature update
	 Since fluxes are diagnostic only, blm_mbl() remains BFB with flx_sfc_lnd() */

      // Adjust temperatures 
      tpt_gnd[lon_idx]+=tpt_gnd_dlt[lon_idx]; // [K]
      
      // Adjust canopy vapor pressure
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_gnd[lon_idx]); // [C]
      if(tpt_bnd_cls > 0.0){
	ppr_H2O_cnp[lon_idx]=ppr_H2O_cnp_cff_a[lon_idx]+ppr_H2O_cnp_cff_b[lon_idx]*svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      }else{
	ppr_H2O_cnp[lon_idx]=ppr_H2O_cnp_cff_a[lon_idx]+ppr_H2O_cnp_cff_b[lon_idx]*svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      } // endif frozen
      
      // Sensible heat flux
      flx_sns_atm_tmp= // [W m-2] Temporary sensible heat flux
	-(tpt_ptn_mdp[lon_idx]-tpt_gnd[lon_idx])
	*dns_mdp[lon_idx]*spc_heat_dry_air
	/rss_aer_heat[lon_idx];
      // Following step approximates psi_h(z0m/L) = 0 
      
      // Monin-Obukhov stability parameter mno_stb_prm for next iteration
      flx_vpr_tmp= // [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
	// 19990501 fxm: this appears to be latent heat flux, not water vapor flux
	-(ppr_H2O_mdp[lon_idx]-ppr_H2O_cnp[lon_idx])
	*dns_mdp[lon_idx]*spc_heat_dry_air
	/(cst_psych[lon_idx]*rss_aer_vpr[lon_idx]);
      flx_sns_atm_vrt_tmp= // [W m-2] Virtual sensible heat flux Bon96 p. 49
	+flx_sns_atm_tmp
	+eps_H2O_rcp_m1*spc_heat_dry_air*tpt_mdp[lon_idx]*flx_vpr_tmp
	/ltn_heat_trn[lon_idx];
      mno_dnm=  // Denominator of Monin-Obukhov length Bon96 p. 49
	+cst_von_krm
	*(grv_sfc_mean/tpt_vrt[lon_idx])
	*flx_sns_atm_vrt_tmp
	/(dns_mdp[lon_idx]*spc_heat_dry_air);
      // Set denominator of Monin-Obukhov length to minimum value if vapor and heat fluxes equal 0.0
      if(PRC_CMP_ABS(mno_dnm) <= eps_dbz) mno_dnm=eps_dbz;
      mno_lng[lon_idx]=-1.0*std::pow(wnd_frc[lon_idx],PRC_CMP(3.0))/mno_dnm; // [m] Monin-Obukhov length Bon96 p. 49
      // Stability functions only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
      mno_stb_prm[lon_idx]=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc] Monin Obukhov stability parameter
      
      // Accumulate number of times stability parameter changes sign 
      if(mno_stb_prm_old[lon_idx]*mno_stb_prm[lon_idx] < 0.0) sgn_chg_ctr[lon_idx]++;
      // Zero stability parameter if it has changed sign too many times
      if(sgn_chg_ctr[lon_idx] >= sgn_chg_ctr_max){
	wrn_prn(prg_nm,sbr_nm,"Numerical ping-pong...setting mno_stb_prm to 0.0");
	mno_stb_prm[lon_idx]=0.0; // [frc]
	// Zero stability corrections for consistency with stability parameter
	mno_stb_crc_mmn[lon_idx]=0.0; // [frc]
	mno_stb_crc_heat[lon_idx]=0.0; // [frc]
      } // endif
      mno_stb_prm_old[lon_idx]=mno_stb_prm[lon_idx]; // [frc]
      
      eps_crr=PRC_CMP_ABS((tpt_gnd[lon_idx]-tpt_gnd_old)/tpt_gnd[lon_idx]); // Relative convergence
      if(dbg_lvl == dbg_sbr){
	/* NB: 
	   Unlike flx_sfc_lnd(), blm_mbl() does not update fluxes after updating 
	   ground temperature each iteration.
	   Thus, these two routines show different flux diagnostics, but do produce
	   identical, BFB results. */
	// 20210618 Emissivity factor was incorrectly applied to Fdwn in below line for ~25 years
	//flx_LW_upw_sfc[lon_idx]=flx_LW_net[lon_idx]+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
	flx_LW_upw_sfc[lon_idx]=flx_LW_net[lon_idx]+flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
	flx_LW_rfl_sfc[lon_idx]=(1.0-msv_sfc[lon_idx])*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave reflected flux at surface
	flx_LW_msn_sfc[lon_idx]=flx_LW_upw_sfc[lon_idx]-flx_LW_rfl_sfc[lon_idx]; // [W m-2] Longwave emitted flux at surface
	(void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
		      itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],tpt_gnd[lon_idx],flx_LW_upw_sfc[lon_idx],flx_sns_atm_tmp,flx_sns_gnd[lon_idx],flx_vpr_tmp,nrg_bdg[lon_idx],eps_crr);
      } // end if dbg
      if(itr_idx > itr_max_gnd){
	std::cerr << "Final: tpt_gnd = " << tpt_gnd[lon_idx] << " K, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Vegetation temperature not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
    
    if(dbg_lvl == dbg_sbr){
      prc_cmp sum_LHS;
      prc_cmp sum_RHS;
      // 20210618 Emissivity factor was incorrectly applied to Fdwn in below line for ~25 years
      //flx_LW_upw_sfc[lon_idx]=flx_LW_net[lon_idx]+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
      flx_LW_upw_sfc[lon_idx]=flx_LW_net[lon_idx]+flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere Bon96 p. 40 Fig. 12, p. 44 
      flx_LW_rfl_sfc[lon_idx]=(1.0-msv_sfc[lon_idx])*flx_LW_dwn_sfc[lon_idx]; // [W m-2] Longwave reflected flux at surface
      flx_LW_msn_sfc[lon_idx]=flx_LW_upw_sfc[lon_idx]-flx_LW_rfl_sfc[lon_idx]; // [W m-2] Longwave emitted flux at surface
      sum_LHS=flx_SW_net[lon_idx]+msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx]; // [W m-2]
      sum_RHS=flx_LW_upw_sfc[lon_idx]+flx_sns_atm[lon_idx]+flx_sns_gnd[lon_idx]+flx_ltn_evp[lon_idx]; // [W m-2]
      (void)std::fprintf(stderr,"SW(abs)  + a*LW(dwn) =  LHS\n");
      (void)std::fprintf(stderr,"%8.3f + %8.3f = %8.3f\n",flx_SW_net[lon_idx],msv_sfc[lon_idx]*flx_LW_dwn_sfc[lon_idx],sum_LHS);
      (void)std::fprintf(stderr," LW(up)  +  H(atm)  +  H(soi)  +  L(atm)  =  RHS\n");
      (void)std::fprintf(stderr,"%8.3f + %8.3f + %8.3f + %8.3f = %8.3f\n",flx_LW_upw_sfc[lon_idx],flx_sns_atm[lon_idx],flx_sns_gnd[lon_idx],flx_ltn_evp[lon_idx],sum_RHS);
    } // end if dbg

  } // end loop over lon
  
  return rcd;
} // end blm_mbl()

int // O [rcd] Return success code
blm_glb // Solve boundary layer meteorology on global scale
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *oro, // I [frc] Orography
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
 prc_cmp *flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *hgt_zpd_dps, // O [m] Zero plane displacement
 prc_cmp *mno_lng_dps, // O [m] Monin-Obukhov length
 prc_cmp *rgh_mmn_dps, // O [m] Roughness length momentum
 prc_cmp *wnd_frc_dps, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_rfr_dps, // O [m s-1] Wind speed at reference height
 prc_cmp *wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl_dps) // O [kg m-1 s-2] Zonal wind stress
{
  /* Purpose: Given vectors of state variables and surface type, 
     divide vectors in land, ocean, and sea-ice components and call the routines 
     to compute the boundary layer exchange properties */
  // Dependencies: 
  
  // Output
  int rcd(0); // O [rcd] Return success code

  // Local
  const prc_cmp hgt_rfr(10.0); // [m] Reference height for deposition processes
  long lon_idx; // Counting index
  // Main Code
  const std::string sbr_nm("blm_glb"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  //  unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  
  // Construct ocean, sea-ice, and land vectors
  std::valarray<bool> flg_lnd(lon_nbr); // [flg] Land flag
  std::valarray<bool> flg_ocn(lon_nbr); // [flg] Ocean flag
  std::valarray<bool> flg_ice(lon_nbr); // [flg] Ice flag
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])) flg_lnd[lon_idx]=true; else flg_lnd[lon_idx]=false;
    if(oro_is_ocn(oro[lon_idx])) flg_ocn[lon_idx]=true; else flg_ocn[lon_idx]=false;
    if(oro_is_ice(oro[lon_idx])) flg_ice[lon_idx]=true; else flg_ice[lon_idx]=false;
    // Sanity check
    if (!flg_lnd[lon_idx] && !flg_ocn[lon_idx] && !flg_ice[lon_idx]) err_prn(prg_nm,sbr_nm,"Invalid surface type in blm_glb()");
  } // end loop over lon
  
  // Midlayer wind speeds
  std::valarray<prc_cmp> wnd_mdp(lon_nbr); // [m s-1] Surface layer mean wind speed
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    wnd_mdp[lon_idx]= // [m s-1] Surface layer mean wind speed
      std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	   wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
  }  // end loop over lon
  
  // Zero-plane displacement
  rcd+=zpd_get
    (lon_nbr, // I [nbr] Size of arrays
     hgt_zpd_dps, // O [m] Zero plane displacement
     oro, // I [frc] Orography
     sfc_typ); // I [idx] LSM surface type (0..28)

  // Fraction of surface covered by snow
  std::valarray<prc_cmp> snw_hgt(lon_nbr); // [m] Geometric bulk thickness of snow
  std::valarray<prc_cmp> snw_frc(lon_nbr); // [frc] Fraction of surface covered by snow
  rcd+=snw_frc_get
    (lon_nbr, // I [nbr] Size of arrays
     &snw_frc[0], // O [frc] Fraction of surface covered by snow
     &snw_hgt[0], // O [m] Geometric bulk thickness of snow
     snw_hgt_lqd); // I [m] Equivalent liquid water snow depth

  // Roughness length
  rcd+=rgh_mmn_get
    (lon_nbr, // I [nbr] Size of arrays
     oro, // I [frc] Orography
     rgh_mmn_dps, // O [m] Roughness length momentum
     sfc_typ, // I [idx] LSM surface type (0..28)
     &snw_frc[0], // I [frc] Fraction of surface covered by snow
     &wnd_mdp[0]); // I [m s-1] Surface layer mean wind speed

  // Land boundary layer
  rcd+=blm_lnd
    (lon_nbr, // I [nbr] Size of arrays
     &flg_lnd[0], // I [flg] Land flag
     dns_mdp, // I [kg m-3] Midlayer density
     hgt_mdp, // I [m] Midlayer height above surface
     prs_mdp, // I [Pa] Pressure
     q_H2O_vpr, // I [kg kg-1] Specific humidity
     rgh_mmn_dps, // I [m] Roughness length momentum
     tpt_mdp, // I [K] Midlayer temperature
     tpt_ptn_mdp, // I [K] Midlayer local potential temperature
     tpt_sfc, // I [K] Surface temperature
     wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
     flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
     flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
     flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
     mno_lng_dps, // O [m] Monin-Obukhov length
     wnd_frc_dps, // O [m s-1] Surface friction velocity
     wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
     wnd_str_znl_dps); // O [kg m-1 s-2] Zonal wind stress

  // Ocean boundary layer
  rcd+=blm_ocn
    (lon_nbr, // I [nbr] Size of arrays
     &flg_ocn[0], // I [flg] Ocean flag
     dns_mdp, // I [kg m-3] Midlayer density
     hgt_mdp, // I [m] Midlayer height above surface
     prs_mdp, // I [Pa] Pressure
     q_H2O_vpr, // I [kg kg-1] Specific humidity
     tpt_mdp, // I [K] Midlayer temperature
     tpt_ptn_mdp, // I [K] Midlayer local potential temperature
     tpt_sfc, // I [K] Surface temperature
     wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
     flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
     flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
     flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
     mno_lng_dps, // O [m] Monin-Obukhov length
     rgh_mmn_dps, // O [m] Roughness length momentum
     wnd_frc_dps, // O [m s-1] Surface friction velocity
     wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
     wnd_str_znl_dps); // O [kg m-1 s-2] Zonal wind stress

  // Sea ice boundary layer
  rcd+=blm_ice
    (lon_nbr, // I [nbr] Size of arrays
     &flg_ice[0], // I [flg] Sea ice flag
     dns_mdp, // I [kg m-3] Midlayer density
     hgt_mdp, // I [m] Midlayer height above surface
     prs_mdp, // I [Pa] Pressure
     q_H2O_vpr, // I [kg kg-1] Specific humidity
     rgh_mmn_dps, // I [m] Roughness length momentum
     tpt_mdp, // I [K] Midlayer temperature
     tpt_ptn_mdp, // I [K] Midlayer local potential temperature
     tpt_sfc, // I [K] Sea surface temperature
     wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
     flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
     flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
     flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
     mno_lng_dps, // O [m] Monin-Obukhov length
     wnd_frc_dps, // O [m s-1] Surface friction velocity
     wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
     wnd_str_znl_dps); // O [kg m-1 s-2] Zonal wind stress

  // 10 m windspeed
  rcd+=wnd_rfr_get
    (lon_nbr, // I [nbr] Size of arrays
     hgt_mdp, // I [m] Midpoint height above surface
     hgt_rfr, // I [m] Reference height for deposition processes
     hgt_zpd_dps, // I [m] Zero plane displacement
     mno_lng_dps, // I [m] Monin-Obukhov length
     wnd_frc_dps, // I [m s-1] Surface friction velocity
     &wnd_mdp[0], // I [m s-1] Surface layer mean wind speed
     wnd_rfr_dps); // O [m s-1] Wind speed at reference height

  return rcd;
} // end blm_glb()

int // O [rcd] Return success code
blm_ice
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_ice, // I [flg] Sea ice flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl) // O [kg m-1 s-2] Zonal wind stress
{
  /* Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
     compute the boundary layer exchange properties over sea ice
     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxsice(), Bonan (1996) CCM:lsm/surtem()
  */
  // Dependencies: <phys_cst.hh>,<tdy.hh>,<blm.hh>
  /* Notes on algorithm:
     Suffix mdp quantity evaluated at height hgt_mdp
     Suffix sfc quantity evaluated at surface temperature
     The roughness length of sea ice depends on many factors
     A globally uniform roughness length for sea ice is foolish but necessary
     CCM1/2/3 used z0m=0.04 m (CCM:dom/parpbl), based on a glacial ice value
     LSM also adopted this value for land ice
     This value is more appropriate for ridged, multi-year ice
     The NCAR Oceanography section used z0m=0.05 m (BKL97 p. F-3)
     These values lead to large drag coefficients and excessive ice extent off of Antarctica
     CSM later adopted z0m=0.0005 m, appropriate for new-formed, seasonal ice
     This value greatly improved sea-ice dynamics in CSM, so we adopt it
     This blm_ice() routine reads roughness length as an input parameter
     for algorithmic reasons, and checks that its value agrees with CCM sea ice
  */
  
  // Output
  int rcd(0); // O [rcd] Return success code

  // Local
  const prc_cmp hgt_rfr_LaP81(10.0); // [m] Reference height for turbulent flux parameterization of LaP81
  const prc_cmp hgt_rfr_tpt(2.0); // [m] Reference height for temperature
#ifndef NDEBUG
  const prc_cmp rgh_mmn_ice_ocn(0.0005); // [m] Roughness length over sea ice BKL97 p. F-3 (updated)
#endif // NDEBUG
  const prc_cmp wnd_min_dps(1.0); // [m s-1] Minimum windspeed used for deposition
  const long itr_max(5); // Maximum number of iterations for surface flux convergence

  prc_cmp hgt_rat_log; // [frc] Log of ratio of local to reference height
  prc_cmp ltn_heat_trn; // [J kg-1] Latent heat of sublimation or evaporation
  prc_cmp mno_stb_crc_heat; // [frc] Monin-Obukhov stability correction heat
  prc_cmp mno_stb_crc_mmn; // [frc] Monin-Obukhov stability correction momentum
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  std::valarray<prc_cmp> mno_stb_prm(lon_nbr); // [frc] Monin-Obukhov stability parameter 
  prc_cmp ntp_fct; // [frc] Interpolation factor in reference height temperature calculation
  prc_cmp q_H2O_vpr_dlt; // [kg kg-1] Humidity change
  prc_cmp spc_heat_mst_air; // [J kg-1 K-1] Specific heat of moist air
  prc_cmp ssh_H2O_sfc; // [kg kg-1] Saturation specific humidity of H2O at surface
  prc_cmp stb_val(CEWI_cpv); // [flg] 1.0 if stable, 0.0 if unstable
  prc_cmp svp_H2O_sfc; // [Pa] Saturation vapor pressure over planar condensed water at surface
  prc_cmp tpt_bnd_cls; // [C] Temperature bounded celsius
  prc_cmp tpt_dlt; // [K] Temperature change
  prc_cmp tpt_ptn_vrt_mdp; // [K] Midlayer virtual potential temperature
  std::valarray<prc_cmp> tpt_rfr(lon_nbr); // [K] Temperature at reference height
  prc_cmp tpt_rfr_fct_0; // [frc] Factor in reference height temperature calculation
  prc_cmp tpt_rfr_fct_3; // [frc] Factor in reference height temperature calculation
  std::valarray<prc_cmp> tpt_scl(lon_nbr); // [K] Temperature scale
  std::valarray<prc_cmp> vpr_scl(lon_nbr); // [kg kg-1] Moisture scale
  prc_cmp wnd_mdp; // [m s-1] Surface layer mean wind speed
  prc_cmp wnd_mdp_bnd; // [m s-1] Surface layer mean wind speed bounded
  std::valarray<prc_cmp> wnd_rfr_ntr(lon_nbr); // [m s-1] Neutral 10 m wind speed
  std::valarray<prc_cmp> wnd_str(lon_nbr); // [kg m-1 s-2] Wind stress
  prc_cmp xch_cff_heat_ntr_sqrt; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
  prc_cmp xch_cff_heat_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Stanton number for heat exchange
  prc_cmp xch_cff_mmn_ntr_sqrt(CEWI_cpv); // [frc] Squareroot of neutral 10 m drag coefficient
  prc_cmp xch_cff_mmn_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer drag coefficient
  prc_cmp xch_cff_vpr_ntr_sqrt; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
  prc_cmp xch_cff_vpr_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Dalton number for vapor exchange
  prc_cmp xpn_heat_fct; // [frc] Factor of heat exchange
  prc_cmp xpn_mmn_fct; // [frc] Factor of momentum exchange
  long lon_idx; // Counting index for lon

  // Main Code
  const std::string sbr_nm("blm_ice"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  prc_cmp wnd_frc_old; // [m s-1] Surface friction velocity old
  
  // Iteration loop
  const prc_cmp eps_max(1.0e-5); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  long itr_idx; // Counting index
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(flg_ice[lon_idx]){ 

      // Sanity check that routine is being called correctly
      assert(rgh_mmn[lon_idx] == rgh_mmn_ice_ocn); // [m] BKL97 p. F-4

      // Initialize variables which are independent of stability iteration
      // Midlayer wind speeds
      wnd_mdp= // [m s-1] Surface layer mean wind speed
	std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	     wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
      wnd_mdp_bnd=max_cpv(wnd_mdp,wnd_min_dps); // [m s-1] Surface layer mean wind speed bounded
      
      tpt_ptn_vrt_mdp=tpt_ptn_mdp[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Midlayer Virtual potential temperature
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc[lon_idx]); // [C]
      // Assume gridpoint is sea ice thus use heat of sublimation not vaporization (contrary to CSM)
      // NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
      svp_H2O_sfc=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      ltn_heat_trn=ltn_heat_sbl_H2O_std; // [J kg-1]
      ssh_H2O_sfc=q_H2O_vpr_fst_scl_get(svp_H2O_sfc,prs_mdp[lon_idx]); // [kg kg-1] Saturation specific humidity of H2O at surface
      // NB: No salinity correction for sea ice (contrary to CSM)
      
      tpt_dlt=tpt_ptn_mdp[lon_idx]-tpt_sfc[lon_idx]; // [K] Temperature change
      q_H2O_vpr_dlt=q_H2O_vpr[lon_idx]-ssh_H2O_sfc; // [kg kg-1] Humidity change
      hgt_rat_log=std::log(hgt_mdp[lon_idx]/hgt_rfr_LaP81); // [frc] Log of ratio of local to reference height
      spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc); // [J kg-1 K-1] Specific heat of moist air
      
      // Initialize accuracy and counter
      eps_crr=eps_max+1.0; // [frc] Current relative accuracy
      itr_idx=1; // Counting index
      if(dbg_lvl == dbg_sbr){
	(void)std::fprintf(stderr,"Sea ice fluxes:\n");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "itr","mno_lng","mno_stb","wnd_frc","  U10N ","  CDN  ","   CD  "," H(atm)","   L   "," LW(up)","  eps  ");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "   ","   m   ","  frc  "," m s-1 "," m s-1 "," x 1000"," x 1000"," W m-2 "," W m-2 "," W m-2 ","  frc  ");
      } // end if dbg

      // This iteration loop solves for convergence of the wind friction speed
      while(eps_crr > eps_max){

	// Save old friction speed for convergence diagnostic
	if(itr_idx == 1) wnd_frc_old=0.0; else wnd_frc_old=wnd_frc[lon_idx]; // [m s-1]
	
	// Roots of neutral, 10 m exchange coefficients
	xch_cff_mmn_ntr_sqrt=cst_von_krm/std::log(hgt_rfr_LaP81/rgh_mmn[lon_idx]); // [frc] Squareroot of neutral 10 m drag coefficient
	xch_cff_heat_ntr_sqrt=xch_cff_mmn_ntr_sqrt; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
	xch_cff_vpr_ntr_sqrt=xch_cff_mmn_ntr_sqrt; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange

	if(itr_idx == 1){
	  // First iteration: Use estimated roots of neutral exchange coefficients
	  wnd_frc[lon_idx]=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_ntr_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	}else{
	  // Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
	  wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	} // endif first iteration
	
	// Compute stability parameter at midlayer and evaluate stability corrections  
	// Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
	mno_stb_prm[lon_idx]= // [frc] Monin-Obukhov stability parameter 
	  cst_von_krm*grv_sfc_mean*hgt_mdp[lon_idx]
	  *(tpt_scl[lon_idx]/tpt_ptn_vrt_mdp+vpr_scl[lon_idx]/(1.0/eps_H2O_rcp_m1+q_H2O_vpr[lon_idx]))
	  /(wnd_frc[lon_idx]*wnd_frc[lon_idx]);
	// Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
	mno_stb_prm[lon_idx]=sign_cpv(min_cpv(PRC_CMP_ABS(mno_stb_prm[lon_idx]),10.0),mno_stb_prm[lon_idx]); // [frc] Monin-Obukhov stability parameter
	stb_val=0.5+sign_cpv(0.5,mno_stb_prm[lon_idx]); // [flg] 1.0 if stable, 0.0 if unstable
	sml_fnc_mmn_uns_rcp=max_cpv(std::sqrt(PRC_CMP_ABS(1.0-16.0*mno_stb_prm[lon_idx])),1.0); // [frc] BKL97 p. F1, LaP81 p. 325
	sml_fnc_mmn_uns_rcp=std::sqrt(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_mmn=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_heat=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	
	// Shift old neutral 10 m exchange coefficients to measurement height and stability
	xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)); // [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
	xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
	xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
	
	eps_crr=PRC_CMP_ABS((wnd_frc[lon_idx]-wnd_frc_old)/wnd_frc[lon_idx]); // Relative convergence
	if(dbg_lvl == dbg_sbr){
	  wnd_rfr_ntr[lon_idx]=wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt; // [m s-1] Neutral 10 m wind speed (by definition)
	  mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
	  wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
	  flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
	  // fxm: units problem again with flx_ltn
	  flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
	  flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere
	  (void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
			itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],wnd_rfr_ntr[lon_idx],1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,flx_sns_atm[lon_idx],flx_ltn[lon_idx],flx_LW_upw_sfc[lon_idx],eps_crr);
	} // end if dbg
	if(itr_idx > itr_max){
	  std::cerr << "Final: wnd_frc = " << wnd_frc[lon_idx] << " m s-1, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	  wrn_prn(prg_nm,sbr_nm,"Wind friction speed not converging, breaking loop...");
	  break;
	} // endif
	itr_idx++;
      } // end loop over itr
      
      // Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
      wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity 
      tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale 
      vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale 
      
      // Compute surface stress components
      wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
      wnd_str_znl[lon_idx]=-wnd_str[lon_idx]*wnd_znl_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Zonal wind stress
      wnd_str_mrd[lon_idx]=-wnd_str[lon_idx]*wnd_mrd_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Meridional wind stress
      
      /* Compute heat flux components at current surface temperature
	 Define positive latent and sensible heat as upwards into atmosphere */
      flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
      flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere

      /* Following Geleyn (1988), interpolate tpt_sfc to fixed height hgt_rfr_tpt
	 Compute function of exchange coefficients
	 Assume that xch_cff_mmn_ntr=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, xch_cff_mmn=xch_cff_mmn_sqrt*xch_cff_mmn_sqrt and xch_cff_heat=xch_cff_heat_sqrt*xch_cff_mmn_sqrt, and therefore 1/sqrt(xch_cff_mmn_ntr)=1/xch_cff_mmn_ntr_sqrt and sqrt(xch_cff_mmn)/xch_cff_heat=1/xch_cff_heat_sqrt  
      */
      xpn_mmn_fct=cst_von_krm/xch_cff_mmn_ntr_sqrt; // [frc] Exponential factor of momentum exchange
      xpn_heat_fct=cst_von_krm/xch_cff_heat_sqrt; // [frc] Exponential factor of heat exchange

      // Interpolation factor for stable and unstable cases
      tpt_rfr_fct_0=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct)-1.0)); // [frc]
      tpt_rfr_fct_3=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct-xpn_heat_fct)-1.0)); // [frc]
      ntp_fct= // [frc] Interpolation factor in reference height temperature calculation
	(tpt_rfr_fct_0-hgt_rfr_tpt/hgt_mdp[lon_idx]*(xpn_mmn_fct-xpn_heat_fct))/xpn_heat_fct*stb_val
	+(tpt_rfr_fct_0-tpt_rfr_fct_3)/xpn_heat_fct*(1.0-stb_val);
      ntp_fct=min_cpv(max_cpv(ntp_fct,0.0),1.0); // [frc]

      // Actual interpolation
      tpt_rfr[lon_idx]=tpt_sfc[lon_idx]+(tpt_mdp[lon_idx]-tpt_sfc[lon_idx])*ntp_fct; // [K]

      // Additional diagnostics
      mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
    } // end if flg_ice
  } // end loop over lon
  
  return rcd;
} // end blm_ice()

int // O [rcd] Return success code
blm_lnd
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_lnd, // I [flg] Land flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl) // O [kg m-1 s-2] Zonal wind stress
{
  /* Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
     compute the boundary layer exchange properties over land
     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxsice(), Bonan (1996) CCM:lsm/surtem()
  */
  // Dependencies: <phys_cst.hh>,<tdy.hh>,<blm.hh>
  /* Notes on algorithm:
     Suffix mdp quantity evaluated at height hgt_mdp
     Suffix sfc quantity evaluated at surface temperature
     Currently, routine is virtually identical to blm_ice() */
  
  // Output
  int rcd(0); // O [rcd] Return success code

  // Local
  const prc_cmp hgt_rfr_LaP81(10.0); // [m] Reference height for turbulent flux parameterization of LaP81
  const prc_cmp hgt_rfr_tpt(2.0); // [m] Reference height for temperature
  const prc_cmp wnd_min_dps(1.0); // [m s-1] Minimum windspeed used for deposition
  const long itr_max(25); // Maximum number of iterations for surface flux convergence

  prc_cmp hgt_rat_log; // [frc] Log of ratio of local to reference height
  prc_cmp ltn_heat_trn; // [J kg-1] Latent heat of sublimation or evaporation
  prc_cmp mno_stb_crc_heat(CEWI_cpv); // [frc] Monin-Obukhov stability correction heat
  prc_cmp mno_stb_crc_mmn(CEWI_cpv); // [frc] Monin-Obukhov stability correction momentum
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  std::valarray<prc_cmp> mno_stb_prm(lon_nbr); // [frc] Monin-Obukhov stability parameter 
  prc_cmp ntp_fct; // [frc] Interpolation factor in reference height temperature calculation
  prc_cmp q_H2O_vpr_dlt; // [kg kg-1] Humidity change
  prc_cmp spc_heat_mst_air; // [J kg-1 K-1] Specific heat of moist air
  prc_cmp ssh_H2O_sfc; // [kg kg-1] Saturation specific humidity of H2O at surface
  prc_cmp stb_val(CEWI_cpv); // [flg] 1.0 if stable, 0.0 if unstable
  prc_cmp svp_H2O_sfc; // [Pa] Saturation vapor pressure over planar condensed water at surface
  prc_cmp tpt_bnd_cls; // [C] Temperature bounded celsius
  prc_cmp tpt_dlt; // [K] Temperature change
  prc_cmp tpt_ptn_vrt_mdp; // [K] Midlayer virtual potential temperature
  std::valarray<prc_cmp> tpt_rfr(lon_nbr); // [K] Temperature at reference height
  prc_cmp tpt_rfr_fct_0; // [frc] Factor in reference height temperature calculation
  prc_cmp tpt_rfr_fct_3; // [frc] Factor in reference height temperature calculation
  std::valarray<prc_cmp> tpt_scl(lon_nbr); // [K] Temperature scale
  std::valarray<prc_cmp> vpr_scl(lon_nbr); // [kg kg-1] Moisture scale
  prc_cmp wnd_mdp; // [m s-1] Surface layer mean wind speed
  prc_cmp wnd_mdp_bnd; // [m s-1] Surface layer mean wind speed bounded
  std::valarray<prc_cmp> wnd_rfr_ntr(lon_nbr); // [m s-1] Neutral 10 m wind speed
  std::valarray<prc_cmp> wnd_str(lon_nbr); // [kg m-1 s-2] Wind stress
  prc_cmp xch_cff_heat_ntr_sqrt; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
  prc_cmp xch_cff_heat_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Stanton number for heat exchange
  prc_cmp xch_cff_mmn_ntr_sqrt(CEWI_cpv); // [frc] Squareroot of neutral 10 m drag coefficient
  prc_cmp xch_cff_mmn_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer drag coefficient
  prc_cmp xch_cff_vpr_ntr_sqrt; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
  prc_cmp xch_cff_vpr_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Dalton number for vapor exchange
  prc_cmp xpn_heat_fct; // [frc] Factor of heat exchange
  prc_cmp xpn_mmn_fct; // [frc] Factor of momentum exchange
  long lon_idx; // Counting index for lon

  // Main Code
  const std::string sbr_nm("blm_lnd"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  prc_cmp wnd_frc_old; // [m s-1] Surface friction velocity old
  
  // Iteration loop
  const prc_cmp eps_max(1.0e-5); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  long itr_idx; // Counting index
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(flg_lnd[lon_idx]){ 

      // Initialize variables which are independent of stability iteration
      // Midlayer wind speeds
      wnd_mdp= // [m s-1] Surface layer mean wind speed
	std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	     wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
      wnd_mdp_bnd=max_cpv(wnd_mdp,wnd_min_dps); // [m s-1] Surface layer mean wind speed bounded
      
      tpt_ptn_vrt_mdp=tpt_ptn_mdp[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Midlayer Virtual potential temperature
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc[lon_idx]); // [C]
      // Saturation vapor pressure of water at ground temperature
      // NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
      if(tpt_bnd_cls > 0.0){ 
	svp_H2O_sfc=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	ltn_heat_trn=ltn_heat_vpr_H2O_std; // [J kg-1]
      }else{
	svp_H2O_sfc=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
	ltn_heat_trn=ltn_heat_sbl_H2O_std; // [J kg-1]
      } // endif
      ssh_H2O_sfc=q_H2O_vpr_fst_scl_get(svp_H2O_sfc,prs_mdp[lon_idx]); // [kg kg-1] Saturation specific humidity of H2O at surface
      // NB: No salinity correction for land
      
      tpt_dlt=tpt_ptn_mdp[lon_idx]-tpt_sfc[lon_idx]; // [K] Temperature change
      q_H2O_vpr_dlt=q_H2O_vpr[lon_idx]-ssh_H2O_sfc; // [kg kg-1] Humidity change
      hgt_rat_log=std::log(hgt_mdp[lon_idx]/hgt_rfr_LaP81); // [frc] Log of ratio of local to reference height
      spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc); // [J kg-1 K-1] Specific heat of moist air
      
      // Initialize accuracy and counter
      eps_crr=eps_max+1.0; // [frc] Current relative accuracy
      itr_idx=1; // Counting index
      if(dbg_lvl == dbg_sbr){
	(void)std::fprintf(stderr,"Land fluxes:\n");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "itr","mno_lng","mno_stb","wnd_frc","  U10N ","  CDN  ","   CD  "," H(atm)","   L   "," LW(up)","  eps  ");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "   ","   m   ","  frc  "," m s-1 "," m s-1 "," x 1000"," x 1000"," W m-2 "," W m-2 "," W m-2 ","  frc  ");
      } // end if dbg

      // This iteration loop solves for convergence of the wind friction speed
      while(eps_crr > eps_max){

	// Save old friction speed for convergence diagnostic
	if(itr_idx == 1) wnd_frc_old=0.0; else wnd_frc_old=wnd_frc[lon_idx]; // [m s-1]
	
	// Roots of neutral, 10 m exchange coefficients
	// dcm: trying limits on rgh_mmn to rgh_mmn_bnd
	const prc_cmp rgh_mmn_max(1.0); // [m] Maximum roughness length momentum
	prc_cmp rgh_mmn_bnd; // [m] Roughness length momentum bounded
	rgh_mmn_bnd=min_cpv(rgh_mmn[lon_idx],rgh_mmn_max); // [m] Roughness length momentum bounded
	if(rgh_mmn[lon_idx] > rgh_mmn_max && itr_idx == 1) std::cerr << "WARNING: "+sbr_nm+" bounding roughness length at 1.0 to facilitate blm convergence" << std::endl;
	rgh_mmn_bnd=min_cpv(rgh_mmn[lon_idx],rgh_mmn_max); // [m] Roughness length momentum bounded

	xch_cff_mmn_ntr_sqrt=cst_von_krm/std::log(hgt_rfr_LaP81/rgh_mmn_bnd); // [frc] Squareroot of neutral 10 m drag coefficient
	xch_cff_heat_ntr_sqrt=xch_cff_mmn_ntr_sqrt; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
	xch_cff_vpr_ntr_sqrt=xch_cff_mmn_ntr_sqrt; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange

	if(itr_idx == 1){
	  // First iteration: Use estimated roots of neutral exchange coefficients
	  wnd_frc[lon_idx]=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_ntr_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	}else{
	  // Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
	  wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	} // endif first iteration
	
	// Compute stability parameter at midlayer and evaluate stability corrections  
	// Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
	mno_stb_prm[lon_idx]= // [frc] Monin-Obukhov stability parameter 
	  cst_von_krm*grv_sfc_mean*hgt_mdp[lon_idx]
	  *(tpt_scl[lon_idx]/tpt_ptn_vrt_mdp+vpr_scl[lon_idx]/(1.0/eps_H2O_rcp_m1+q_H2O_vpr[lon_idx]))
	  /(wnd_frc[lon_idx]*wnd_frc[lon_idx]);
	// Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
	mno_stb_prm[lon_idx]=sign_cpv(min_cpv(PRC_CMP_ABS(mno_stb_prm[lon_idx]),10.0),mno_stb_prm[lon_idx]); // [frc] Monin-Obukhov stability parameter

	stb_val=0.5+sign_cpv(0.5,mno_stb_prm[lon_idx]); // [flg] 1.0 if stable, 0.0 if unstable
	sml_fnc_mmn_uns_rcp=max_cpv(std::sqrt(PRC_CMP_ABS(1.0-16.0*mno_stb_prm[lon_idx])),1.0); // [frc] BKL97 p. F1, LaP81 p. 325
	sml_fnc_mmn_uns_rcp=std::sqrt(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_mmn=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_heat=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	
	// Shift old neutral 10 m exchange coefficients to measurement height and stability
	xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)); // [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
	xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
	xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
	
	eps_crr=PRC_CMP_ABS((wnd_frc[lon_idx]-wnd_frc_old)/wnd_frc[lon_idx]); // Relative convergence
	if(dbg_lvl == dbg_sbr){
	  wnd_rfr_ntr[lon_idx]=wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt; // [m s-1] Neutral 10 m wind speed (by definition)
	  mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
	  wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
	  flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
	  // fxm: units problem again with flx_ltn
	  flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
	  flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere
	  (void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
			itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],wnd_rfr_ntr[lon_idx],1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,flx_sns_atm[lon_idx],flx_ltn[lon_idx],flx_LW_upw_sfc[lon_idx],eps_crr);
	} // end if dbg
	if(itr_idx > itr_max){
	  std::cerr << "Final: wnd_frc = " << wnd_frc[lon_idx] << " m s-1, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	  wrn_prn(prg_nm,sbr_nm,"Wind friction speed not converging, breaking loop...");
	  break;
	} // endif
	itr_idx++;
      } // end loop over itr
      
      // Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
      wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity 
      tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale 
      vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale 
      
      // Compute surface stress components
      wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
      wnd_str_znl[lon_idx]=-wnd_str[lon_idx]*wnd_znl_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Zonal wind stress
      wnd_str_mrd[lon_idx]=-wnd_str[lon_idx]*wnd_mrd_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Meridional wind stress
      
      /* Compute heat flux components at current surface temperature
	 Define positive latent and sensible heat as upwards into atmosphere */
      flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
      flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere

      /* Following Geleyn (1988), interpolate tpt_sfc to fixed height hgt_rfr_tpt
	 Compute function of exchange coefficients
	 Assume that xch_cff_mmn_ntr=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, xch_cff_mmn=xch_cff_mmn_sqrt*xch_cff_mmn_sqrt and xch_cff_heat=xch_cff_heat_sqrt*xch_cff_mmn_sqrt, and therefore 1/sqrt(xch_cff_mmn_ntr)=1/xch_cff_mmn_ntr_sqrt and sqrt(xch_cff_mmn)/xch_cff_heat=1/xch_cff_heat_sqrt  
      */
      xpn_mmn_fct=cst_von_krm/xch_cff_mmn_ntr_sqrt; // [frc] Exponential factor of momentum exchange
      xpn_heat_fct=cst_von_krm/xch_cff_heat_sqrt; // [frc] Exponential factor of heat exchange

      // Interpolation factor for stable and unstable cases
      tpt_rfr_fct_0=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct)-1.0)); // [frc]
      tpt_rfr_fct_3=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct-xpn_heat_fct)-1.0)); // [frc]
      ntp_fct= // [frc] Interpolation factor in reference height temperature calculation
	(tpt_rfr_fct_0-hgt_rfr_tpt/hgt_mdp[lon_idx]*(xpn_mmn_fct-xpn_heat_fct))/xpn_heat_fct*stb_val
	+(tpt_rfr_fct_0-tpt_rfr_fct_3)/xpn_heat_fct*(1.0-stb_val);
      ntp_fct=min_cpv(max_cpv(ntp_fct,0.0),1.0); // [frc]

      // Actual interpolation
      tpt_rfr[lon_idx]=tpt_sfc[lon_idx]+(tpt_mdp[lon_idx]-tpt_sfc[lon_idx])*ntp_fct; // [K]

      // Additional diagnostics
      mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
      // Test that neutral drag coefficient inverts correctly
      //      rgh_mmn[lon_idx]=hgt_rfr_LaP81*std::exp(-cst_von_krm/xch_cff_mmn_ntr_sqrt); // [m] BKL97 p. F-4, LaP81 p. 327 (14) 
    } // end if flg_lnd
  } // end loop over lon
  
  return rcd;
} // end blm_lnd()

int // O [rcd] Return success code
blm_ocn
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_ocn, // I [flg] Ocean flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling (emission+reflection) flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl) // O [kg m-1 s-2] Zonal wind stress
{
  /* Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
     compute the boundary layer exchange properties over open ocean
     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem() 
  */
  // Dependencies: <phys_cst.hh>,<tdy.hh>,<blm.hh>
  /* Validation:
     Attempted to validate implementation of blm_ocn() against LaP81 and LaP82
     Simulated LaP81 p. 333 Fig 7 using
     mie --tpt_sst=277 --tpt_mdp=281 --hgt_mdp=12.5 --wnd_znl_mdp=20
     Agreement is within 15% as wnd_znl_mdp varies from 7.5--20.0 m s-1
  */
  /* Notes on algorithm:
     In this routine surface temperature tpt_sfc = Sea surface temperature SST
     Routine uses sfc suffix rather than SST simply for consistency with land and sea ice routines
     Suffix mdp quantity evaluated at height hgt_mdp
     Suffix sfc quantity evaluated at (sea) surface temperature
  */
  
  // Output
  int rcd(0); // O [rcd] Return success code

  // Local
  const prc_cmp dlt_nbr_ntr_10m_cff(0.0346); // [frc] Coefficient for neutral 10 m Dalton number CCM:dom/flxoce() LaP82 p. 477
  const prc_cmp hgt_rfr_LaP81(10.0); // [m] Reference height for turbulent flux parameterization of LaP81
  const prc_cmp hgt_rfr_tpt(2.0); // [m] Reference height for temperature
  const prc_cmp stn_nbr_ntr_10m_cff_stb(0.0180); // [frc] Coefficient for neutral, stable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
  const prc_cmp stn_nbr_ntr_10m_cff_uns(0.0327); // [frc] Coefficient for neutral, unstable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
  const prc_cmp wnd_min_dps(1.0); // [m s-1] Minimum windspeed used for deposition
  const long itr_max(5); // Maximum number of iterations for surface flux convergence

  prc_cmp hgt_rat_log; // [frc] Log of ratio of local to reference height
  prc_cmp ltn_heat_trn; // [J kg-1] Latent heat of sublimation or evaporation
  prc_cmp mno_stb_crc_heat; // [frc] Monin-Obukhov stability correction heat
  prc_cmp mno_stb_crc_mmn(CEWI_cpv); // [frc] Monin-Obukhov stability correction momentum
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  std::valarray<prc_cmp> mno_stb_prm(lon_nbr); // [frc] Monin-Obukhov stability parameter 
  prc_cmp ntp_fct; // [frc] Interpolation factor in reference height temperature calculation
  prc_cmp q_H2O_vpr_dlt; // [kg kg-1] Humidity change
  prc_cmp spc_heat_mst_air; // [J kg-1 K-1] Specific heat of moist air
  prc_cmp ssh_H2O_sfc; // [kg kg-1] Saturation specific humidity of H2O at surface
  prc_cmp stb_val; // [flg] 1.0 if stable, 0.0 if unstable
  prc_cmp svp_H2O_sfc; // [Pa] Saturation vapor pressure over planar condensed water at surface
  prc_cmp tpt_bnd_cls; // [C] Temperature bounded celsius
  prc_cmp tpt_dlt; // [K] Temperature change
  prc_cmp tpt_ptn_vrt_mdp; // [K] Midlayer virtual potential temperature
  std::valarray<prc_cmp> tpt_rfr(lon_nbr); // [K] Temperature at reference height
  prc_cmp tpt_rfr_fct_0; // [frc] Factor in reference height temperature calculation
  prc_cmp tpt_rfr_fct_3; // [frc] Factor in reference height temperature calculation
  std::valarray<prc_cmp> tpt_scl(lon_nbr); // [K] Temperature scale
  std::valarray<prc_cmp> vpr_scl(lon_nbr); // [kg kg-1] Moisture scale
  prc_cmp wnd_mdp; // [m s-1] Surface layer mean wind speed
  prc_cmp wnd_mdp_bnd; // [m s-1] Surface layer mean wind speed bounded
  std::valarray<prc_cmp> wnd_rfr_ntr(lon_nbr); // [m s-1] Neutral 10 m wind speed
  std::valarray<prc_cmp> wnd_str(lon_nbr); // [kg m-1 s-2] Wind stress
  prc_cmp xch_cff_heat_ntr_sqrt; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
  prc_cmp xch_cff_heat_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Stanton number for heat exchange
  prc_cmp xch_cff_mmn_ntr_sqrt; // [frc] Squareroot of neutral 10 m drag coefficient
  prc_cmp xch_cff_mmn_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer drag coefficient
  prc_cmp xch_cff_vpr_ntr_sqrt; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
  prc_cmp xch_cff_vpr_sqrt(CEWI_cpv); // [frc] Squareroot of mid-layer Dalton number for vapor exchange
  prc_cmp xpn_heat_fct; // [frc] Factor of heat exchange
  prc_cmp xpn_mmn_fct; // [frc] Factor of momentum exchange
  long lon_idx; // Counting index for lon
  std::valarray<prc_cmp> wnd_rfr(lon_nbr); // [m s-1] Wind speed at reference height LaP81 p. 327 (14)
  std::valarray<prc_cmp> xch_cff_mmn(lon_nbr); // [frc] Drag coefficient at mid-layer
  // std::valarray<prc_cmp> xch_cff_heat(lon_nbr); // [frc] Exchange coefficient for heat transfer
  //  std::valarray<prc_cmp> xch_cff_mmn(lon_nbr); // [frc] Exchange coefficient for momentum transfer
  //  std::valarray<prc_cmp> xch_cff_vpr(lon_nbr); // [frc] Exchange coefficient for vapor transfer

  // Main Code
  const std::string sbr_nm("blm_ocn"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  prc_cmp wnd_frc_old; // [m s-1] Surface friction velocity old
  
  // Iteration loop
  const prc_cmp eps_max(1.0e-5); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  long itr_idx; // Counting index
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(flg_ocn[lon_idx]){ 

      // Initialize variables which are independent of stability iteration
      // Midlayer wind speeds
      wnd_mdp= // [m s-1] Surface layer mean wind speed
	std::sqrt(wnd_znl_mdp[lon_idx]*wnd_znl_mdp[lon_idx]+ 
	     wnd_mrd_mdp[lon_idx]*wnd_mrd_mdp[lon_idx]); 
      wnd_mdp_bnd=max_cpv(wnd_mdp,wnd_min_dps); // [m s-1] Surface layer mean wind speed bounded
      
      tpt_ptn_vrt_mdp=tpt_ptn_mdp[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Midlayer Virtual potential temperature
      tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc[lon_idx]); // [C]
      /* Liquid sea water supercools befored freezing to sea ice (CCM:dom/parsst.h/tsice = -1.7999 C)
	 Assume gridpoint liquid sea water thus use heat of vaporization not sublimation */
      svp_H2O_sfc=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls); // [Pa]
      ltn_heat_trn=ltn_heat_vpr_H2O_std; // [J kg-1]
      ssh_H2O_sfc=q_H2O_vpr_fst_scl_get(svp_H2O_sfc,prs_mdp[lon_idx]); // [kg kg-1] Saturation specific humidity of H2O at surface
      // Correct freshwater saturation specific humidity for salinity effects
      ssh_H2O_sfc*=ssh_H2O_sln_crc; // [kg kg-1] Saturation specific humidity of H2O at surface
      
      tpt_dlt=tpt_ptn_mdp[lon_idx]-tpt_sfc[lon_idx]; // [K] Temperature change
      q_H2O_vpr_dlt=q_H2O_vpr[lon_idx]-ssh_H2O_sfc; // [kg kg-1] Humidity change
      hgt_rat_log=std::log(hgt_mdp[lon_idx]/hgt_rfr_LaP81); // [frc] Log of ratio of local to reference height
      spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc); // [J kg-1 K-1] Specific heat of moist air
      
      // Stable if tpt_ptn_mdp > tpt_sfc
      stb_val=0.5+sign_cpv(0.5,tpt_dlt); // [flg] 1.0 if stable, 0.0 if unstable
      // Initial guess for roots of neutral exchange coefficients: z/L=0 and u10n=vmag
      xch_cff_mmn_ntr_sqrt=std::sqrt(xch_cff_mmn_ocn_ntr_get(wnd_mdp_bnd)); // [frc] Squareroot of neutral 10 m drag coefficient
      xch_cff_heat_ntr_sqrt=(1.0-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
      xch_cff_vpr_ntr_sqrt=dlt_nbr_ntr_10m_cff; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange

      // Initialize accuracy and counter
      eps_crr=eps_max+1.0; // [frc] Current relative accuracy
      itr_idx=1; // Counting index
      if(dbg_lvl == dbg_sbr){
	(void)std::fprintf(stderr,"Ocean fluxes:\n");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "itr","mno_lng","mno_stb","wnd_frc","  U10N ","  CDN  ","   CD  "," H(atm)","   L   "," LW(up)","  eps  ");
	(void)std::fprintf(stderr,"%3s %9s %8s %7s %7s %7s %8s %8s %8s %8s %8s\n",
		      "   ","   m   ","  frc  "," m s-1 "," m s-1 "," x 1000"," x 1000"," W m-2 "," W m-2 "," W m-2 ","  frc  ");
      } // end if dbg

      // This iteration loop solves for convergence of the wind friction speed
      while(eps_crr > eps_max){

	// Save old friction speed for convergence diagnostic
	if(itr_idx == 1) wnd_frc_old=0.0; else wnd_frc_old=wnd_frc[lon_idx]; // [m s-1]
	
	if(itr_idx == 1){
	  // First iteration: Use estimated roots of neutral exchange coefficients
	  wnd_frc[lon_idx]=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_ntr_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	}else{
	  // Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
	  wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity
	  tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale
	  vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale
	} // endif first iteration
	
	// Compute stability parameter at midlayer and evaluate stability corrections  
	// Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
	mno_stb_prm[lon_idx]= // [frc] Monin-Obukhov stability parameter 
	  cst_von_krm*grv_sfc_mean*hgt_mdp[lon_idx]
	  *(tpt_scl[lon_idx]/tpt_ptn_vrt_mdp+vpr_scl[lon_idx]/(1.0/eps_H2O_rcp_m1+q_H2O_vpr[lon_idx]))
	  /(wnd_frc[lon_idx]*wnd_frc[lon_idx]);
	// Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
	mno_stb_prm[lon_idx]=sign_cpv(min_cpv(PRC_CMP_ABS(mno_stb_prm[lon_idx]),10.0),mno_stb_prm[lon_idx]); // [frc] Monin-Obukhov stability parameter
	stb_val=0.5+sign_cpv(0.5,mno_stb_prm[lon_idx]); // [flg] 1.0 if stable, 0.0 if unstable
	sml_fnc_mmn_uns_rcp=max_cpv(std::sqrt(PRC_CMP_ABS(1.0-16.0*mno_stb_prm[lon_idx])),1.0); // [frc] BKL97 p. F1, LaP81 p. 325
	sml_fnc_mmn_uns_rcp=std::sqrt(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_mmn=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	mno_stb_crc_heat=-5.0*mno_stb_prm[lon_idx]*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325
	
	// Shift old neutral 10 m exchange coefficient to measurement height and stability 
	xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)); // [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
	
	// Define neutral 10 m wind speed 
	wnd_rfr_ntr[lon_idx]=wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt; // [m s-1] Neutral 10 m wind speed
	
	// Update neutral 10 m exchange coefficients
	xch_cff_mmn_ntr_sqrt=std::sqrt(xch_cff_mmn_ocn_ntr_get(wnd_rfr_ntr[lon_idx])); // [frc] Squareroot of neutral 10 m drag coefficient
	xch_cff_vpr_ntr_sqrt=dlt_nbr_ntr_10m_cff; // [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
	xch_cff_heat_ntr_sqrt=(1.0-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb; // [frc] Squareroot of neutral 10 m Stanton number for heat exchange
	
	// Shift old neutral 10 m exchange coefficients to measurement height and stability
	xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)); // [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
	xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
	xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)); // [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
	
	eps_crr=PRC_CMP_ABS((wnd_frc[lon_idx]-wnd_frc_old)/wnd_frc[lon_idx]); // Relative convergence
	if(dbg_lvl == dbg_sbr){
	  mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
	  wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
	  flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
	  // fxm: units problem again with flx_ltn
	  flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
	  flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere
	  (void)std::fprintf(stderr,"%3ld %9.3f %8.3f %7.4f %7.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.6f\n",
			itr_idx,mno_lng[lon_idx],mno_stb_prm[lon_idx],wnd_frc[lon_idx],wnd_rfr_ntr[lon_idx],1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,flx_sns_atm[lon_idx],flx_ltn[lon_idx],flx_LW_upw_sfc[lon_idx],eps_crr);
	} // end if dbg
	if(itr_idx > itr_max){
	  std::cerr << "Final: wnd_frc = " << wnd_frc[lon_idx] << " m s-1, mno_lng = " << mno_lng[lon_idx] << " m, eps_crr = " << eps_crr << std::endl;
	  wrn_prn(prg_nm,sbr_nm,"Wind friction speed not converging, breaking loop...");
	  break;
	} // endif
	itr_idx++;
      } // end loop over itr
      
      // Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
      wnd_frc[lon_idx]=xch_cff_mmn_sqrt*wnd_mdp_bnd; // [m s-1] Surface friction velocity 
      tpt_scl[lon_idx]=xch_cff_heat_sqrt*tpt_dlt; // [K] Temperature scale 
      vpr_scl[lon_idx]=xch_cff_vpr_sqrt*q_H2O_vpr_dlt; // [kg kg-1] Moisture scale 
      
      // Compute surface stress components
      wnd_str[lon_idx]=dns_mdp[lon_idx]*wnd_frc[lon_idx]*wnd_frc[lon_idx]; // [kg m-1 s-2] Wind stress
      wnd_str_znl[lon_idx]=-wnd_str[lon_idx]*wnd_znl_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Zonal wind stress
      wnd_str_mrd[lon_idx]=-wnd_str[lon_idx]*wnd_mrd_mdp[lon_idx]/wnd_mdp_bnd; // [kg m-1 s-2] Meridional wind stress
      
      /* Compute heat flux components at current surface temperature
	 Define positive latent and sensible heat as upwards into atmosphere */
      flx_sns_atm[lon_idx]=-spc_heat_mst_air*wnd_str[lon_idx]*tpt_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Sensible heat flux to atmosphere
      flx_ltn[lon_idx]=-ltn_heat_trn*wnd_str[lon_idx]*vpr_scl[lon_idx]/wnd_frc[lon_idx]; // [W m-2] Latent heat flux to atmosphere
      flx_LW_upw_sfc[lon_idx]=cst_Stefan_Boltzmann*std::pow(tpt_sfc[lon_idx],PRC_CMP(4.0)); // [W m-2] Longwave upwelling (emission+reflection) flux to atmosphere

      /* Following Geleyn (1988), interpolate tpt_sfc to fixed height hgt_rfr_tpt
	 Compute function of exchange coefficients
	 Assume that xch_cff_mmn_ntr=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, xch_cff_mmn=xch_cff_mmn_sqrt*xch_cff_mmn_sqrt and xch_cff_heat=xch_cff_heat_sqrt*xch_cff_mmn_sqrt, and therefore 1/sqrt(xch_cff_mmn_ntr)=1/xch_cff_mmn_ntr_sqrt and sqrt(xch_cff_mmn)/xch_cff_heat=1/xch_cff_heat_sqrt  
      */
      xpn_mmn_fct=cst_von_krm/xch_cff_mmn_ntr_sqrt; // [frc] Exponential factor of momentum exchange
      xpn_heat_fct=cst_von_krm/xch_cff_heat_sqrt; // [frc] Exponential factor of heat exchange

      // Interpolation factor for stable and unstable cases
      tpt_rfr_fct_0=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct)-1.0)); // [frc]
      tpt_rfr_fct_3=std::log(1.0+(hgt_rfr_tpt/hgt_mdp[lon_idx])*(std::exp(xpn_mmn_fct-xpn_heat_fct)-1.0)); // [frc]
      ntp_fct= // [frc] Interpolation factor in reference height temperature calculation
	(tpt_rfr_fct_0-hgt_rfr_tpt/hgt_mdp[lon_idx]*(xpn_mmn_fct-xpn_heat_fct))/xpn_heat_fct*stb_val
	+(tpt_rfr_fct_0-tpt_rfr_fct_3)/xpn_heat_fct*(1.0-stb_val);
      ntp_fct=min_cpv(max_cpv(ntp_fct,0.0),1.0); // [frc]

      // Actual interpolation
      tpt_rfr[lon_idx]=tpt_sfc[lon_idx]+(tpt_mdp[lon_idx]-tpt_sfc[lon_idx])*ntp_fct; // [K]

      // Additional diagnostics
      mno_lng[lon_idx]=hgt_mdp[lon_idx]/mno_stb_prm[lon_idx]; // [m] Monin-Obukhov length
      rgh_mmn[lon_idx]=hgt_rfr_LaP81*std::exp(-cst_von_krm/xch_cff_mmn_ntr_sqrt); // [m] BKL97 p. F-4, LaP81 p. 327 (14) 

      wnd_rfr[lon_idx]= // [m s-1] Wind speed at reference height LaP81 p. 327 (14)
	+wnd_mdp
	-(wnd_frc[lon_idx]/cst_von_krm)
	*(hgt_rat_log-mno_stb_crc_mmn+mno_stb_crc_mmn_get(hgt_rfr_LaP81/mno_lng[lon_idx]));
      xch_cff_mmn[lon_idx]=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt; // [frc] Drag coefficient at mid-layer

    } // end if flg_ocn
  } // end loop over lon
  
  return rcd;
} // end blm_ocn()
