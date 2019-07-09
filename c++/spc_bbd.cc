// $Id$ 

// Implementation (declaration) of blackbody spectra classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <spc_bbd.hh> // Blackbody spectra

// spc_bbd class

int spc_bbd_cls::nst_nbr=0; // [nbr] Number of instantiated class members
const long spc_bbd_cls::bnd_nbr_dfl(100); // [nbr] Default number of bands to discretize Planck function

spc_bbd_cls::spc_bbd_cls(const prc_cmp tpt_arg) // Default constructor
{
  // Purpose: Constructor for spc_bbd objects
  rcd=0; // [enm] Return success code

  bnd_nbr=bnd_nbr_dfl; // [nbr] Number of bands to discretize Planck function

  rcd+=tpt_set(tpt_arg); // [K] Blackbody temperature

  rcd+=recompute(); // [fnc] Recompute properties of object

  nst_nbr++; // [nbr] Number of instantiated class members
} // end spc_bbd_cls::spc_bbd_cls()

spc_bbd_cls::~spc_bbd_cls(){ // Destructor
  nst_nbr--; // [nbr] Number of instantiated class members
} // end spc_bbd_cls destructor

prc_cmp spc_bbd_cls::flx_ttl()const{return flx_hms;} // [W m-2] Hemispheric blackbody irradiance

prc_cmp spc_bbd_cls::tpt_get()const{return static_cast<prc_cmp>(tpt);} // [K] Blackbody temperature
int // [enm] Return success code
spc_bbd_cls::tpt_set(const prc_cmp &tpt_arg) // [K] Blackbody temperature
{
  // Purpose: Set temperature
  std::string sbr_nm("spc_bbd_cls::tpt_set"); // [sng] Subroutine name
  const prc_cmp tpt_max(1400.0); // [K] Maximum allowed blackbody temperature

  assert(tpt_arg >= 0.0);
  if(tpt_arg > tpt_max) err_prn(sbr_nm,"tpt_arg too large for natural source, exiting");
  tpt=tpt_arg; // [W m-2] Solar constant
  rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end spc_bbd_cls::tpt_set()

int // [enm] Return success code
spc_bbd_cls::recompute(){ // [fnc] Recompute properties of object
  // Set private members which lack set() functions
  std::string sbr_nm("spc_bbd_cls::recompute()"); // [sng] Subroutine name

  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using phc::cst_Stefan_Boltzmann; // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462
  // Blackbodies emit isotropically so hemispheric irradiance is pi times radiance
  flx_hms=cst_Stefan_Boltzmann*std::pow(tpt,4.0); // [W m-2] Hemispheric blackbody irradiance = cst_stf_blt * T^4
  ntn_ntg=flx_hms/mth::cst_M_PIl; // [W m-2 sr-1] Intensity of blackbody radiation = (cst_stf_blt*T^4)/pi

  return rcd; // [enm] Return success code
} // end spc_bbd_cls::recompute()

prc_cmp // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
spc_bbd_cls::eval(const prc_cmp &wvl){ // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  // Compute specific intensity of radiation at wavelength wvl emitted by blackbody of temperature tpt
  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant
  using phc::hc; // (1.986488377e-25) [J m] Planck constant times speed of light = hc
  using phc::hc2; // (5.9553531e-17) [J m2 s-1] Planck constant times speed of light squared = hc2
  xpn=hc/(wvl*cst_Boltzmann*tpt); // [frc] Exponent in Planck function
  // Guard against overflows
#ifdef PRC_FLT
  const prc_cmp xpn_bnd(std::log(FLT_MAX)); // [frc] Maximum exponent in Planck function
#else // !PRC_FLT
  const prc_cmp xpn_bnd(std::log(DBL_MAX)); // [frc] Maximum exponent in Planck function
#endif // !PRC_FLT
  if(xpn < xpn_bnd){
    dnm=std::exp(xpn)-1.0; // [frc] Denominator of Planck function
    ntn_bbd_wvl=2.0*hc2/(std::pow(wvl,PRC_CMP(5.0))*dnm); // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  }else{
    std::cerr << "spc_bbd_cls::eval(): tpt = " << tpt << " K, wvl = " << wvl*1.0e6 << " um, xpn = " << xpn << std::endl;
    wrn_prn("spc_bbd_cls::eval","Out of range, setting blackbody intensity to 0.0 W m-2 m-1 sr-1");
    ntn_bbd_wvl=0.0; // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  } // end else
  return ntn_bbd_wvl; // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
} // end spc_bbd_cls::eval()

int // O [enm] Return success code
spc_bbd_cls::flx_frc_get // [fnc] Fraction of blackbody emission in given spectral region
(const prc_cmp *wvl_min, // I [m] Minimum wavelength
 const prc_cmp *wvl_max, // I [m] Maximum wavelength
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 prc_cmp *flx_IR_frc, // O [frc] Fraction of infrared flux in band
 const long &bnd_nbr_arg) // I [nbr] Number of bands to discretize Planck function
{
  /* Purpose: Compute fraction of blackbody emission in given spectral region
     Routine uses following nomentclature (consistent with mie()):
     wavelengths (wvl): user-specified input/output grid
     bands (bnd): internal sub-grid to improve accuracy in wide wavelength bins
     The optional number bnd_nbr determines the number of sub-bands used
     to discretize the Planck function.
     bnd_nbr needs to be fairly large to accurately integrate the Planck
     function over wide spectral bands.
     fxm: find fast routine that optimizes quadrature of Planck function */
     
  std::string sbr_nm("spc_bbd_cls::flx_frc_get"); // [sng] Subroutine name
  long bnd_idx; // [idx] Counting index for bnd
  long wvl_idx; // [idx] Counting index for wvl
  double ntn_ttl; // [W m-2 sr-1] Radiance in specified spectral region

  // In case user specified a new discretization resolution...
  bnd_nbr=bnd_nbr_arg; // [nbr] Number of bands to discretize Planck function

  /* Instantiate wavelength grid over whole domain with only bnd_nbr elements
     This particular grid will be re-focused before it is used */
  prc_cmp wvl_mnm_arg(wvl_min[0]); // [m] Minimum wavelength
  prc_cmp wvl_mxm_arg(wvl_max[wvl_nbr-1]); // [m] Maximum wavelength
  wvl_grd_cls wvlgrd(static_cast<std::string>("regular"),wvl_mnm_arg,wvl_mxm_arg,bnd_nbr);

  /* Point to required center and width components
     Pointers will not change, although values will */
  prc_cmp * const bnd_ctr=wvlgrd.wvl_ctr_get(); // [m] Wavelength at band center
  prc_cmp * const bnd_dlt=wvlgrd.wvl_dlt_get(); // [m] Bandwidth

  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    ntn_ttl=0.0; // [W m-2 sr-1] Radiance in specified spectral region
    // Wavelengths must increase in each region to prevent negative flux fractions
    if(wvl_min[wvl_idx] > wvl_max[wvl_idx]) err_prn(sbr_nm,"wvl_min > wvl_max");
    // Re-grid bands to current wavelength
    wvlgrd.wvl_grd_min_max_set(wvl_min[wvl_idx],wvl_max[wvl_idx]);
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      ntn_ttl+=eval(bnd_ctr[bnd_idx])*bnd_dlt[bnd_idx]; // [W m-2 sr-1] Radiance in specified spectral region
    } // end loop over bnd
    // Convert [W m-2 sr-1] -> [frc] for fractional blackbody emission
    flx_IR_frc[wvl_idx]=ntn_ttl/ntn_ntg; // [frc] Fraction of infrared flux in band
  } // end loop over wvl
  
  return rcd; // [enm] Return success code
} // end spc_bbd_cls::flx_frc_get()

prc_cmp // [frc] Fraction of blackbody emission in given spectral region
spc_bbd_cls::flx_frc_get // [fnc] Fraction of blackbody emission in given spectral region
(const prc_cmp &wvl_min, // I [m] Minimum wavelength
 const prc_cmp &wvl_max) // I [m] Maximum wavelength
{
  /* Purpose: Compute fraction of blackbody emission in given spectral region
     Routine is simply a wrapper for array routine
     Error checking done in array routine */

  prc_cmp flx_IR_frc; // [frc] Fraction of IR flux in band
  long wvl_nbr(1L); // [nbr] Number of wavelength bands
  std::string sbr_nm("spc_bbd_cls::flx_frc_get"); // [sng] Subroutine name

  rcd+=flx_frc_get // [fnc] Fraction of blackbody emission in given spectral region
    (&wvl_min, // I [m] Minimum wavelength
     &wvl_max, // I [m] Maximum wavelength
     wvl_nbr, // I [nbr] Number of wavelength bands
     &flx_IR_frc); // O [frc] Fraction of IR flux in band

  return flx_IR_frc; // [frc] Fraction of IR flux in band
} // end spc_bbd_cls::flx_frc_get()

