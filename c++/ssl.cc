// $Id$ 

// Purpose: Sea salt physics utilities for C++ programs

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <ssl.hh> // Sea salt physics

int // O [rcd] Return success code
cst_C1C2C3_And98_get // [fnc] 
(const prc_cmp wnd_10m, // I [m s-1] Wind speed at 10 m reference height
 prc_cmp &cst_C1_And98, // O [frc] Factor in direct production scheme And98
 prc_cmp &cst_C2_And98, // O [frc] Factor in direct production scheme And98
 prc_cmp &cst_C3_And98) // O [frc] Factor in direct production scheme And98
{
  // Purpose: Evaluate C1 constant of And98 based on current wind speed
  // Output
  int rcd(0); // O [rcd] Return success code
  const prc_cmp rds_10p0_mcr(10.0); // [um] Radius
  const prc_cmp rds_37p5_mcr(37.5); // [um] Radius
  const prc_cmp rds_100p0_mcr(100.0); // [um] Radius

  prc_cmp flx_vrt_nbr_ssl_SPC93; // [# m-2 s-1 m-1] Spectral vertical number flux of sea salt particles at 80% RH

  // [fnc] Evaluate number flux of 10 micron sea salt particles using SPC93
  flx_vrt_nbr_ssl_SPC93=flx_vrt_nbr_ssl_10mcr_SPC93_get(wnd_10m); // [# m-2 s-1 m-1]
  // C1 determined by ensuring continuity of And98 and SPC93 at 10 microns
  cst_C1_And98=rds_10p0_mcr*flx_vrt_nbr_ssl_SPC93/wnd_10m; // [frc] And98 p. 2180 (3.5a)
  // [fnc] Evaluate number flux of 37.5 micron sea salt particles using And98
  flx_vrt_nbr_ssl_SPC93=cst_C1_And98*wnd_10m*rds_37p5_mcr; // [# m-2 s-1 m-1] And98 p. 2180 (3.5a)
  // C2 determined by ensuring continuity of And98 at 37.5 microns
  cst_C2_And98=std::pow(rds_37p5_mcr,PRC_CMP(2.8))*flx_vrt_nbr_ssl_SPC93/wnd_10m; // [frc] And98 p. 2180 (3.5b)
  // [fnc] Evaluate number flux of 100.0 micron sea salt particles using And98
  flx_vrt_nbr_ssl_SPC93=cst_C2_And98*wnd_10m*rds_100p0_mcr; // [# m-2 s-1 m-1] And98 p. 2180 (3.5b)
  // C3 determined by ensuring continuity of And98 at 100.0 microns
  cst_C3_And98=std::pow(rds_100p0_mcr,PRC_CMP(8.0))*flx_vrt_nbr_ssl_SPC93/wnd_10m; // [frc] And98 p. 2180 (3.5c)

  return rcd; // O [rcd] Return success code
} // end cst_C1C2C3_And98_get() 

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of sea salt
flx_vrt_nbr_ssl_SPC93_get // [fnc] Sea salt generation function of SPC93
(prc_cmp rds_prt, // I [m] Particle radius at ocean surface = 100% RH
 prc_cmp wnd_10m) // I [m s-1] Wind speed at 10 m reference height
{
  /* Purpose: Return vertical number flux of sea salt according to SPC93
     Input is particle radius at 80% relative humidity and wind speed at 10 m */

  // Set four invariant constants used in SPC93
  const prc_cmp fff_SPC93[2]={3.1,3.3}; // [frc] Factor in small spume drop production SPC93
  const prc_cmp rds_mdl_SPC93_mcr[2]={2.1,9.2}; // [um] Empirical mode radii from SPC93

  // Most sea salt parameterizations, including SPC93, want radius in microns
  const prc_cmp rds_prt_mcr(rds_prt*1.0e6); // [um] Size at bin center
  // Note: SPC93 parameterized in terms of r(RH = 80%)
  // And98 suggest r(RH = 80%) = 0.518*r^0.976
  const prc_cmp rds_RH80pct_mcr(0.518*std::pow(rds_prt_mcr,PRC_CMP(0.976))); // [um] Particle radius at 80% RH And98 p. 2180 (3.4)

  // Convert standard 10 m wind speed to 14 m wind speed needed in SPC93
  const prc_cmp xch_cff_mmn_ocn_ntr_10m(xch_cff_mmn_ocn_ntr_get(wnd_10m)); // [frc] Neutral 10 m drag coefficient over ocean
  const prc_cmp fct_wnd_ntp_10m_14m(1.0+std::log(1.4)*std::sqrt(xch_cff_mmn_ocn_ntr_10m)/phc::cst_von_krm); // [frc] Stability factor to interopolate U(10 m) to U(14 m)
  const prc_cmp wnd_14m(wnd_10m*fct_wnd_ntp_10m_14m); // [m s-1] Wind speed at 14 m

  // Relative strengths of two modes are strongly wind speed dependent
  prc_cmp AAA_SPC93[2]; // [frc] Mode strength coefficients SPC93
  AAA_SPC93[0]=std::exp(0.0676*wnd_14m+2.43); // [frc] Mode strength coefficients SPC93 p. 816 (7), And98 p. 2183 (A2a), VDB01 p. 20227 (4)
  // NB: SPC93 has 0.959 but GMS02 incorrectly says 0.0959
  AAA_SPC93[1]=std::exp(0.959*std::sqrt(wnd_14m))-1.476; // [frc] Mode strength coefficients SPC93 p. 816 (7), And98 p. 2183 (A2b), VDB01 p. 20227 (4)

  prc_cmp fct_log; // [frc] Factor in small spume production SPC93
  prc_cmp flx_nbr_vrt_ssl_SPC93(0.0); // [# m-2 s-1 m-1] Direct sea-salt production SPC93
  for(unsigned short trm_idx=0;trm_idx<=1;trm_idx++){
    fct_log=std::log(rds_RH80pct_mcr/rds_mdl_SPC93_mcr[trm_idx]); // [frc] Factor in small spume production SPC93
    flx_nbr_vrt_ssl_SPC93+=AAA_SPC93[trm_idx]*std::exp(-fff_SPC93[trm_idx]*std::pow(fct_log,PRC_CMP(2.0))); // [# m-2 s-1 m-1] Direct sea-salt production SPC93
  } // end loop over trm

  return flx_nbr_vrt_ssl_SPC93; // [# m-2 s-1 m-1] Direct sea-salt production SPC93
} // end flx_vrt_nbr_ssl_SPC93_get()

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of 10 micron sea salt particles at 80% RH
flx_vrt_nbr_ssl_10mcr_SPC93_get // [fnc] Evaluate number flux of 10 micron sea salt particles using SPC93
(prc_cmp wnd_10m) // I [m s-1] Wind speed at 10 m reference height
{
  /* Purpose: Return spectral vertical number flux of 10 um sea salt according to SPC93
     Input is particle radius at 80% relative humidity and wind speed at 10 m
     This function is to compute C_1 as defined in And98 p. 2180 (3.5)
     Difference with flx_vrt_nbr_ssl_SPC93_get() is that present function
     is for fixed size, and so only wind speed dependent.
     Returned spectral number number flux is based on r=r(80% RH) as in SPC93 */

  // Set four invariant constants used in SPC93
  const prc_cmp fff_SPC93[2]={3.1,3.3}; // [frc] Factor in small spume drop production SPC93
  const prc_cmp rds_mdl_SPC93_mcr[2]={2.1,9.2}; // [um] Empirical mode radii from SPC93

  // Most sea salt parameterizations, including SPC93, want radius in microns
  const prc_cmp rds_RH80pct_mcr(10.0); // [m] Particle radius at 80% RH And98 p. 2180 (3.4)

  // Convert standard 10 m wind speed to 14 m wind speed needed in SPC93
  const prc_cmp xch_cff_mmn_ocn_ntr_10m(xch_cff_mmn_ocn_ntr_get(wnd_10m)); // [frc] Neutral 10 m drag coefficient over ocean
  const prc_cmp fct_wnd_ntp_10m_14m(1.0+std::log(1.4)*std::sqrt(xch_cff_mmn_ocn_ntr_10m)/phc::cst_von_krm); // [frc] Stability factor to interopolate U(10 m) to U(14 m)
  const prc_cmp wnd_14m(wnd_10m*fct_wnd_ntp_10m_14m); // [m s-1] Wind speed at 14 m

  // Relative strengths of two modes are strongly wind speed dependent
  prc_cmp AAA_SPC93[2]; // [frc] Mode strength coefficients SPC93
  AAA_SPC93[0]=std::exp(0.0676*wnd_14m+2.43); // [frc] Mode strength coefficients SPC93 p. 816 (7), And98 p. 2183 (A2a), VDB01 p. 20227 (4)
  // NB: SPC93 has 0.959 but GMS02 incorrectly says 0.0959
  AAA_SPC93[1]=std::exp(0.959*std::sqrt(wnd_14m))-1.476; // [frc] Mode strength coefficients SPC93 p. 816 (7), And98 p. 2183 (A2b), VDB01 p. 20227 (4)

  prc_cmp fct_log; // [frc] Factor in small spume production SPC93
  prc_cmp flx_nbr_vrt_ssl_SPC93(0.0); // [# m-2 s-1 m-1] Direct sea-salt production SPC93
  for(unsigned short trm_idx=0;trm_idx<=1;trm_idx++){
    fct_log=std::log(rds_RH80pct_mcr/rds_mdl_SPC93_mcr[trm_idx]); // [frc] Factor in small spume production SPC93, VDB01 p. 20227 (3)
    flx_nbr_vrt_ssl_SPC93+=AAA_SPC93[trm_idx]*std::exp(-fff_SPC93[trm_idx]*std::pow(fct_log,PRC_CMP(2.0))); // [# m-2 s-1 m-1] Direct sea-salt production SPC93, VDB01 p. 20227 (3)
  } // end loop over 

  return flx_nbr_vrt_ssl_SPC93; // [# m-2 s-1 m-1] Direct sea-salt production SPC93
} // end flx_vrt_nbr_ssl_10mcr_SPC93_get()

int flx_mss_vrt_ssl_get
(const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *sz_ctr, // I [m] Size at bin center
 const prc_cmp *sz_dlt, // I [m] Width of size bin
 const prc_cmp &wnd_10m, // I [m s-1] Wind speed at 10 m reference height
 prc_cmp *flx_mss_vrt_ssl_drc, // O [kg m-2 s-1 m-1] Direct sea-salt production And98
 prc_cmp *flx_mss_vrt_ssl_ndr, // O [kg m-2 s-1 m-1] Indirect sea-salt production MSD86
 prc_cmp *flx_mss_vrt_ssl, // O [kg m-2 s-1 m-1] Vertical mass flux of sea-salt
 prc_cmp *ssl_flx_ndr_drc_rat) // O [m-1] Ratio of direct to indirect sea-salt mass flux
{     
  /* Purpose: Diagnose total vertical mass flux of sea-salt from 10-m wind speed
     Direct production accounts for spume drops in large size mode
     Andreas (1998) (And98) is used for direct production for sizes > 21 um
     Smith et al. (SPC93) is used for direct production for sizes < 21 um
     Indirect production accounts for bursting bubbles and jet drops
     These are generally in the small size mode
     Monahan et al. (1986) (MSD86) is used for indirect production
     Following Andreas (1998) all particle radii are corrected from equilibrium
     values at 80% relative humidity to values at point of departure from ocean
     surface, which corresponds to 100% relative humidity.
     Thus fluxes predicted are numbers of particles of size r per unit r at sea surface
     Mass fluxes are obtained by multiplying by particle volume and density
     Hygroscopic growth formulae must be used to convert particle size to dry particle size
     Theory: 
     Dependencies: <phys_cst.hh> */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using phc::dns_NaCl_std; // (2170.0) [kg m-3] Density of NaCl http://www.crystran.co.uk/nacldata.htm
  const std::string sbr_nm("flx_mss_vrt_ssl_get"); // Subroutine name
  const prc_cmp fff_SPC93[2]={3.1,3.3}; // [frc] Factor in small spume drop production SPC93
  const prc_cmp rds_mdl_SPC93_mcr[2]={2.1,9.2}; // [um] Empirical mode radii from SPC93
  long sz_idx; // Counting index for sz
  prc_cmp cst_C1_And98; // [frc] Factor in direct production scheme And98
  prc_cmp cst_C2_And98; // [frc] Factor in direct production scheme And98
  prc_cmp cst_C3_And98; // [frc] Factor in direct production scheme And98
  prc_cmp bbb_MSD86; // [frc] Factor in indirect production scheme MSD86
  prc_cmp dns_prt_ssl; // [kg m-3] Density of sea salt aerosol
  prc_cmp rds_prt_mcr; // [um] Size at bin center
  prc_cmp rds_RH80pct_mcr; // [um] Particle radius at 80% RH
  prc_cmp fct_log; // [frc] Factor in small spume production SPC93
  prc_cmp vlm_prt_mcr; // [um] Volume at bin center
  std::valarray<prc_cmp> flx_nbr_vrt_ssl_drc(sz_nbr); // [# m-2 s-1 m-1] Direct sea-salt production And98
  std::valarray<prc_cmp> flx_nbr_vrt_ssl_ndr(sz_nbr); // [# m-2 s-1 m-1] Indirect sea-salt production MSD86
  std::valarray<prc_cmp> flx_nbr_vrt_ssl(sz_nbr); // [# m-2 s-1 m-1] Vertical number flux of sea-salt

  // Convert standard 10 m wind speed to 14 m wind speed needed in SPC93
  const prc_cmp xch_cff_mmn_ocn_ntr_10m(xch_cff_mmn_ocn_ntr_get(wnd_10m)); // [frc] Neutral 10 m drag coefficient over ocean
  const prc_cmp fct_wnd_ntp_10m_14m(1.0+std::log(1.4)*std::sqrt(xch_cff_mmn_ocn_ntr_10m)/phc::cst_von_krm); // [frc] Stability factor to interopolate U(10 m) to U(14 m)
  const prc_cmp wnd_14m(wnd_10m*fct_wnd_ntp_10m_14m); // [m s-1] Wind speed at 14 m

  // Relative strengths of two modes are strongly wind speed dependent
  prc_cmp AAA_SPC93[2]; // [frc] Mode strength coefficients SPC93
  AAA_SPC93[0]=std::exp(0.0676*wnd_14m+2.43); // [frc] Mode strength coefficients SPC93 p. 816 (7) And98 p. 2183 (A2a)
  // NB: SPC93 has 0.959 but GMS02 incorrectly says 0.0959
  AAA_SPC93[1]=std::exp(0.959*std::sqrt(wnd_14m))-1.476; // [frc] Mode strength coefficients SPC93 p. 816 (7) And98 p. 2183 (A2b)

  rcd+=cst_C1C2C3_And98_get // [fnc] Evaluate spume production parameters of And98
    (wnd_10m, // I [m s-1] Wind speed at 10 m reference height
     cst_C1_And98, // O [frc] Factor in direct production scheme And98
     cst_C2_And98, // O [frc] Factor in direct production scheme And98
     cst_C3_And98); // O [frc] Factor in direct production scheme And98

  // Initialize variables that will be incremented
  prc_cmp flx_mss_vrt_ssl_ttl(0.0); // [kg m-2 s-1 m-1] Vertical mass flux of sea-salt
  prc_cmp flx_mss_vrt_ssl_drc_ttl(0.0); // [kg m-2 s-1 m-1] Direct sea-salt production And98
  prc_cmp flx_mss_vrt_ssl_ndr_ttl(0.0); // [kg m-2 s-1 m-1] Indirect sea-salt production MSD86
  for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    rds_prt_mcr=sz_ctr[sz_idx]*1.0e6; // [um] Size at bin center
    // fxm: verify r should be in microns in following expressions
    bbb_MSD86=(0.380-std::log(rds_prt_mcr))/0.650; // [frc] Factor in indirect production scheme MSD86 GMS02 (2)
    // NB: Dimensions are spectral number flux not spectral mass flux
    flx_nbr_vrt_ssl_ndr[sz_idx]=std::pow(wnd_10m,PRC_CMP(3.41))*std::pow(rds_prt_mcr,PRC_CMP(3.0))*(1.0+0.057*std::pow(rds_prt_mcr,PRC_CMP(1.05)))*std::pow(PRC_CMP(10.0),PRC_CMP(1.19)*std::exp(-bbb_MSD86*bbb_MSD86)); // [# m-2 s-1 m-1] Indirect sea-salt production MSD86 GMS02 (2)
    // Direct production of sea-salt aerosol smaller than 21 um uses Smith et al. (SPC93)
    if(rds_prt_mcr<21.0){
      // Direct production of sea-salt aerosol smaller than 21 um uses Smith et al. (SPC93)
      flx_nbr_vrt_ssl_drc[sz_idx]=0.0; // [# m-2 s-1 m-1] Direct sea-salt production SPC93
      // Note: 0.518*r^0.976 =  r(RH = 80%) used in SPC93
      rds_RH80pct_mcr=0.518*std::pow(rds_prt_mcr,PRC_CMP(0.976)); // [um] Particle radius at 80% RH And98 p. 2180 (3.4)
      for(unsigned short trm_idx=0;trm_idx<=1;trm_idx++){
	fct_log=std::log(rds_RH80pct_mcr/rds_mdl_SPC93_mcr[trm_idx]); // [frc] Factor in small spume production SPC93
	flx_nbr_vrt_ssl_drc[sz_idx]+=AAA_SPC93[trm_idx]*std::exp(-fff_SPC93[trm_idx]*std::pow(fct_log,PRC_CMP(2.0))); // [# m-2 s-1 m-1] Direct sea-salt production SPC93
      } // end loop over 
      // fxm: GMS02 drops this term
      flx_nbr_vrt_ssl_drc[sz_idx]*=0.506*std::pow(rds_prt_mcr,PRC_CMP(-0.024)); // [# m-2 s-1 m-1] Direct sea-salt production SPC93
    }else if(rds_prt_mcr<80.0){
      flx_nbr_vrt_ssl_drc[sz_idx]=3.5*cst_C1_And98*wnd_10m*0.977/rds_prt_mcr; // [# m-2 s-1 m-1] Direct sea-salt production And98
    }else if(rds_prt_mcr<220.0){
      flx_nbr_vrt_ssl_drc[sz_idx]=3.5*cst_C2_And98*wnd_10m*3.19*std::pow(rds_prt_mcr,PRC_CMP(-2.757)); // [# m-2 s-1 m-1] Direct sea-salt production And98
    }else if(rds_prt_mcr<562.0){
      flx_nbr_vrt_ssl_drc[sz_idx]=3.5*cst_C3_And98*wnd_10m*97.62*std::pow(rds_prt_mcr,PRC_CMP(-7.832)); // [# m-2 s-1 m-1] Direct sea-salt production And98
    }else{
      flx_nbr_vrt_ssl_drc[sz_idx]=0.0; // [# m-2 s-1 m-1] Direct sea-salt production And98
      err_prn(sbr_nm,"No production information for this size of sea-salt");
    } // endif
    // Sum direct + indirect to get total production
    flx_nbr_vrt_ssl[sz_idx]=flx_nbr_vrt_ssl_drc[sz_idx]+flx_nbr_vrt_ssl_ndr[sz_idx]; // [# m-2 s-1 m-1] Vertical mass flux of sea-salt

    // Convert number flux to mass flux
    // fxm: This assumes aerosol is same density as sodium chloride, must correct density for deliquescence
    dns_prt_ssl=dns_NaCl_std*1.0; // [kg m-3] Density of sea salt aerosol
    vlm_prt_mcr=(4.0/3.0)*cst_M_PIl*std::pow(rds_prt_mcr,PRC_CMP(3.0)); // [um] Volume at bin center
    flx_nbr_vrt_ssl_drc[sz_idx]=flx_nbr_vrt_ssl_drc[sz_idx]*vlm_prt_mcr*dns_prt_ssl; // [# m-2 s-1 m-1] Direct sea-salt production SPC93
    flx_mss_vrt_ssl_ndr[sz_idx]=flx_nbr_vrt_ssl_ndr[sz_idx]*vlm_prt_mcr*dns_prt_ssl; // [kg m-2 s-1 m-1] Indirect sea-salt production MSD86 GMS02 (2)

    // Diagnostic mass fluxes
    flx_mss_vrt_ssl[sz_idx]=flx_mss_vrt_ssl_drc[sz_idx]+flx_mss_vrt_ssl_ndr[sz_idx]; // [# m-2 s-1 m-1] Vertical mass flux of sea-salt
    ssl_flx_ndr_drc_rat[sz_idx]=flx_mss_vrt_ssl_ndr[sz_idx]/flx_mss_vrt_ssl[sz_idx]; // O [m-1] Ratio of direct to indirect sea-salt mass flux

    // Convert spectral mass to absolute mass flux
    flx_mss_vrt_ssl[sz_idx]*=sz_dlt[sz_idx]; // [kg m-2 s-1 m-1] Vertical mass flux of sea-salt
    flx_mss_vrt_ssl_drc[sz_idx]*=sz_dlt[sz_idx]; // [kg m-2 s-1 m-1] Direct sea-salt production And98
    flx_mss_vrt_ssl_ndr[sz_idx]*=sz_dlt[sz_idx]; // [kg m-2 s-1 m-1] Indirect sea-salt production MSD86

    // Size-integrated fluxes
    flx_mss_vrt_ssl_ttl+=flx_mss_vrt_ssl[sz_idx]; // [kg m-2 s-1 m-1] Vertical mass flux of sea-salt
    flx_mss_vrt_ssl_drc_ttl+=flx_mss_vrt_ssl_drc[sz_idx]; // [kg m-2 s-1 m-1] Direct sea-salt production And98
    flx_mss_vrt_ssl_ndr_ttl+=flx_mss_vrt_ssl_ndr[sz_idx]; // [kg m-2 s-1 m-1] Indirect sea-salt production MSD86
  }  // end loop over sz
  
  return rcd;
} // end flx_mss_vrt_ssl_get()

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of sea salt
flx_vrt_nbr_ssl_VDB01_get // [fnc] Sea salt generation function of VDB01
(prc_cmp rds_prt, // I [m] Particle radius at ocean surface = 100% RH
 prc_cmp wnd_10m) // I [m s-1] Wind speed at 10 m reference height
{
  /* Purpose: Return vertical number flux of sea salt according to VDB01
     Input is particle radius at 80% relative humidity and wind speed at 10 m */

  // Output
  //  int rcd(0); // O [rcd] Return success code
  // Local
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  const std::string sbr_nm("flx_vrt_nbr_ssl_VDB01_get"); // Subroutine name

  // Set invariant constants used in VDB01
  //  const prc_cmp rds_nma_mcr[3]={0.2,2.0,12.0}; // [um] Number median radius analytic 
  //  const prc_cmp gsd[3]={1.9,2.0,3.0}; // [frc] Geometric standard deviation

  // Most sea salt parameterizations, including VDB01, want radius in microns
  const prc_cmp rds_prt_mcr(rds_prt*1.0e6); // [um] Size at bin center
  // Note: VDB01 parameterized in terms of r(RH = 80%)
  // And98 suggest r(RH = 80%) = 0.518*r^0.976
  //  const prc_cmp rds_RH80pct_mcr(0.518*std::pow(rds_prt_mcr,PRC_CMP(0.976))); // [um] Particle radius at 80% RH And98 p. 2180 (3.4)

  prc_cmp flx_nbr_vrt_ssl_VDB01(0.0); // [# m-2 s-1 m-1] Direct sea-salt production VDB01
  prc_cmp cnc_ttl[3]; // [# cm-3] Number concentration in each mode
  cnc_ttl[0]=std::pow(PRC_CMP(10.0),PRC_CMP(0.095)*wnd_10m+PRC_CMP(0.283)); // [# cm-3]
  cnc_ttl[1]=std::pow(PRC_CMP(10.0),PRC_CMP(0.0422)*wnd_10m+PRC_CMP(0.288)); // [# cm-3]
  cnc_ttl[2]=std::pow(PRC_CMP(10.0),PRC_CMP(0.069)*wnd_10m-PRC_CMP(3.5)); // [# cm-3]
  for(unsigned short trm_idx=0;trm_idx<=2;trm_idx++){
    flx_nbr_vrt_ssl_VDB01+=cnc_ttl[trm_idx]; // [# m-2 s-1 m-1] Direct sea-salt production VDB01
  } // end loop over trm

  return flx_nbr_vrt_ssl_VDB01; // [# m-2 s-1 m-1] Direct sea-salt production VDB01
} // end flx_vrt_nbr_ssl_VDB01_get()

