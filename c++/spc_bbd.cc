// $Id$ 

// Implementation (declaration) of blackbody spectra classes

/* Copyright (C) 1997--present Charlie Zender
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
} // !spc_bbd_cls::spc_bbd_cls()

spc_bbd_cls::~spc_bbd_cls(){ // Destructor
  nst_nbr--; // [nbr] Number of instantiated class members
} // !spc_bbd_cls destructor

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
  tpt=tpt_arg; // [K] Temperature
  rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !spc_bbd_cls::tpt_set()

int // [enm] Return success code
spc_bbd_cls::recompute(){ // [fnc] Recompute properties of object
  // Set private members which lack set() functions
  std::string sbr_nm("spc_bbd_cls::recompute()"); // [sng] Subroutine name

  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using phc::cst_Stefan_Boltzmann; // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462
  // Blackbodies emit isotropically so hemispheric irradiance is pi times radiance
  flx_hms=cst_Stefan_Boltzmann*std::pow(tpt,4.0); // [W m-2] Hemispheric blackbody irradiance = cst_stf_blt * T^4
  ntn_bbd_ntg=flx_hms/mth::cst_M_PIl; // [W m-2 sr-1] Intensity of blackbody radiation = (cst_stf_blt*T^4)/pi

  return rcd; // [enm] Return success code
} // !spc_bbd_cls::recompute()

int // [enm] Return success code
spc_bbd_cls::plk_ntg_evl // [fnc] Compute integral of Planck function between two wavenumbers
(prc_cmp plk_ntg){ // [W m-2 sr-1] Integrated blackbody emission between wavenumbers

  // Set private members which lack set() functions
  std::string sbr_nm("spc_bbd_cls::planck_xpn()"); // [sng] Subroutine name

  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using phc::cst_Stefan_Boltzmann; // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462

  /* Original DISORT documentation and new variable names:
     A1,2,... : Power series coefficients
     C2       : h * c / k, in units cm*K (h = Plancks constant, c = speed of light, k = Boltzmann constant)
     D(I)     : Exponential series expansion of integral of Planck function from WNUMLO (i=1) or WNUMHI (i=2) to infinity
     EPSIL    : Smallest number such that 1+EPSIL .GT. 1 on computer
     EX       : EXP( - V(I) )
     EXM      : EX**M
     MMAX     : No. of terms to take in exponential series
     MV       : Multiples of V(I)
     P(I)     : Power series expansion of integral of Planck function from zero to WNUMLO (I=1) or WNUMHI (I=2)
     PI       : 3.14159...
     SIGMA    : Stefan-Boltzmann constant (W/m**2/K**4)
     SIGDPI   : SIGMA / PI
     SMALLV   : Number of times the power series is used (0,1,2)
     V(I)     : C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
     VCUT     : Power-series cutoff point
     VCP      : Exponential series cutoff points
     VMAX     : Largest allowable argument of EXP function */

#if 0
  const prc_cmp A1(1.0/3.0); // [] Power series coefficient
  const prc_cmp A2(-1.0/8.0); // [] Power series coefficient
  const prc_cmp A3(1.0/60.0); // [] Power series coefficient
  const prc_cmp A4(-1.0/5040.0); // [] Power series coefficient
  const prc_cmp A5(1.0/272160.0); // [] Power series coefficient
  const prc_cmp A6(-1.0/13305600.0); // [] Power series coefficient

  int i; // []
  int k; // []
  int m; // []
  int MMAX; // [] No. of terms to take in exponential series
  int SMALLV; // [] Number of times the power series is used (0,1,2)
  
  const prc_cmp C2(1.438786); // [] h * c / k, in units cm*K (h = Plancks constant, c = speed of light, k = Boltzmann constant)
  const prc_cmp VCUT(1.5); // [] Power-series cutoff point
  const prc_cmp VCP[7]={10.25,5.7,3.9,2.9,2.3,1.9,0.0}; // [] Exponential series cutoff points

  const prc_cmp CONC(15.0/(cst_M_PIl*cst_M_PIl*cst_M_PIl*cst_M_PIl)); // [] 
  const prc_cmp SIGDPI(cst_Stefan_Boltzmann/cst_M_PIl); // [] SIGMA / PI
  prc_cmp DEL; // [] 
  prc_cmp EPSIL; // [] Smallest number such that 1+EPSIL .GT. 1 on computer
  prc_cmp EX; // [] EXP( - V(I) )
  prc_cmp EXM; // [] EX**M
  prc_cmp HH; // [] 
  prc_cmp MV; // [] Multiples of V(I)
  prc_cmp VAL; // [] 
  prc_cmp VAL0; // [] 
  prc_cmp VMAX; // [] Largest allowable argument of EXP function */
  
  prc_cmp D[2]; // [] Exponential series expansion of integral of Planck function from WNUMLO (i=0) or WNUMHI (i=1) to infinity
  prc_cmp P[2]; // [] Power series expansion of integral of Planck function from zero to WNUMLO (I=0) or WNUMHI (I=1)
  prc_cmp V[2]; // [] C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
#endif // !0
  
  if(tpt < 1.0e-4){
    plk_ntg=0.0; // [W m-2 sr-1]
    return rcd;
  } // !tpt

  return rcd; // [enm] Return success code
} // !spc_bbd_cls::plk_ntg_evl()

prc_cmp // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
spc_bbd_cls::eval(const prc_cmp &wvl){ // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  // Evaluate (compute) specific intensity of radiation at wavelength wvl emitted by blackbody of temperature tpt
  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant
  using phc::hc; // (1.986488377e-25) [J m] Planck constant times speed of light = hc
  using phc::hc2; // (5.9553531e-17) [J m2 s-1] Planck constant times speed of light squared = hc2
  xpn=hc/(wvl*cst_Boltzmann*tpt); // [frc] Exponent in Planck function
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
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
    // 20190709: If warning is too annoying, adopt DISORT polynomial Planck function evaluation?
    if(dbg_lvl > dbg_scl) std::cerr << "spc_bbd_cls::eval(): tpt = " << tpt << " K, wvl = " << wvl*1.0e6 << " um, xpn = " << xpn << ", approximate blackbody emission as 0.0 W m-2 m-1 sr-1" << std::endl;
    ntn_bbd_wvl=0.0; // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  } // !xpn
  return ntn_bbd_wvl; // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
} // !spc_bbd_cls::eval()

int // O [enm] Return success code
spc_bbd_cls::flx_bbd_frc_get // [fnc] Fraction of blackbody emission in given spectral region
(const prc_cmp *wvl_min, // I [m] Minimum wavelength
 const prc_cmp *wvl_max, // I [m] Maximum wavelength
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 prc_cmp *flx_bbd_frc, // O [frc] Fraction of blackbody flux in band
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
     fxm: Implement fast routine that optimizes quadrature of Planck function */
     
  std::string sbr_nm("spc_bbd_cls::flx_bbd_frc_get"); // [sng] Subroutine name
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
    } // !bnd_idx
    // Convert [W m-2 sr-1] -> [frc] for fractional blackbody emission
    flx_bbd_frc[wvl_idx]=ntn_ttl/ntn_bbd_ntg; // [frc] Fraction of blackbody flux in band
  } // !wvl_idx
  
  return rcd; // [enm] Return success code
} // !spc_bbd_cls::flx_bbd_frc_get()

prc_cmp // [frc] Fraction of blackbody emission in given spectral region
spc_bbd_cls::flx_bbd_frc_get // [fnc] Fraction of blackbody emission in given spectral region
(const prc_cmp &wvl_min, // I [m] Minimum wavelength
 const prc_cmp &wvl_max) // I [m] Maximum wavelength
{
  /* Purpose: Compute fraction of blackbody emission in given spectral region
     Routine is simply a wrapper for array routine
     Error checking done in array routine */

  prc_cmp flx_bbd_frc; // [frc] Fraction of blackbody flux in band
  long wvl_nbr(1L); // [nbr] Number of wavelength bands
  std::string sbr_nm("spc_bbd_cls::flx_bbd_frc_get"); // [sng] Subroutine name

  rcd+=flx_bbd_frc_get // [fnc] Fraction of blackbody emission in given spectral region
    (&wvl_min, // I [m] Minimum wavelength
     &wvl_max, // I [m] Maximum wavelength
     wvl_nbr, // I [nbr] Number of wavelength bands
     &flx_bbd_frc); // O [frc] Fraction of blackbody flux in band

  return flx_bbd_frc; // [frc] Fraction of blackbody flux in band
} // !spc_bbd_cls::flx_bbd_frc_get()

int // O [enm] Return success code
flx_bbd_frc_get_WiW76 // [fnc] Fraction of blackbody emission in given spectral region
(const prc_cmp *wvn_grd, // I [cm-1] Wavenumber at band interfaces
 const long &wvn_nbr, // I [nbr] Number of wavenumber bands (interfaces minus one)
 const prc_cmp &tpt, // I [K] Temperature
 prc_cmp *flx_bbd_frc) // O [frc] Fraction of blackbody flux in band
{
  /* Purpose: Compute fraction of blackbody emission in each band of wavenumber grid wvn_grd
     Routine uses following nomentclature (consistent with mie()):
     wvn_grd[wvn_nbr+1]: User-specified input/output grid in CGS wavenumber [cm-1]
     flx_bbd_frc[wvn_nbr]: [frc] Fraction of blackbody flux in band
     
     Routine supports two normalization options:
     flx_frc_nrm: Re-normalize absolute fractions to sum to unity. This normalization ensures all thermal emission is represented in wvl_grd. The sum of flx_bbd_frc is exactly one.
     flx_frc_abs: Return absolute flux fractions. The fractions will sum to unity iff the grid spans all thermal wavelengths
     
     Assumptions:
     Wavenumber grid is contiguous, monotonic, and non-overlapping
     This allows planck_integral_WiW76() for each interior wavenumber interface to be re-used
     Wavenumber grid and temperature, together, are in Earth's parameter space */
     
  // These constants are only used diagnostically in this routine
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using phc::cst_Stefan_Boltzmann; // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462
  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant (2018 SI NIST)
  using phc::cst_Planck; // (6.62606876e-34) [J s] Planck's constant (CODATA, 2018 SI NIST) (exact)
  using phc::speed_of_light; // (2.99792458e+08) [m s-1] Speed of light in vacuo (CODATA, 2018 SI NIST)

  std::string sbr_nm("flx_bbd_frc_get_WiW76"); // [sng] Subroutine name
  std::string nrm_typ_sng("flx_frc_abs"); // [sng] Normalization type
  int rcd(0); // [enm] Return code
  
  double *ntn_bbd_blr=new double[wvn_nbr+1]; // [W m-2 sr-1] Integrated radiance bluer than corresponding wavenumber

  double ntn_bbd_ntg; // [W m-2 sr-1] Total integrated intensity of blackbody radiation = (cst_stf_blt*T^4)/pi
  double ntn_bbd_rsl; // [W m-2 sr-1] Integrated radiance resolved in single band
  double ntn_bbd_rsl_ttl; // [W m-2 sr-1] Total integrated radiance resolved in wavenumber grid
  double flx_bbd_frc_mss_blr; // [frc] Fractional "missing" radiance bluer (shorter wavelengths/larger wavenumbers) than wavenumber grid
  double flx_bbd_frc_mss_rdr; // [frc] Fractional "missing" radiance redder (longer wavelengths/smaller wavenumbers) than wavenumber grid
  double wvn_lo(DBL_MIN); // [cm-1] Lowest wavenumber in current bin

  long wvn_idx; // [idx] Counting index for wvn

  ntn_bbd_ntg=cst_Stefan_Boltzmann*std::pow(tpt,4.0)/cst_M_PIl; // [W m-2 sr-1] Intensity of blackbody radiation = (cst_stf_blt*T^4)/pi

  for(wvn_idx=0;wvn_idx<=wvn_nbr;wvn_idx++){ // NB: wvn_nbr+1 iterations

    // Wavenumbers must increase monotonically to prevent negative flux fractions
    if(wvn_idx < wvn_nbr){
            //      if(wvn_grd[wvn_idx+1] <= wvn_grd[wvn_idx]) err_prn(sbr_nm,"wvn_min > wvn_max");
      if(wvn_grd[wvn_idx] < wvn_grd[wvn_idx+1]){
	wvn_lo=wvn_grd[wvn_idx];
      }else{
	wvn_lo=wvn_grd[wvn_idx+1];
      } // !wvn_grd
    } // !wvn_idx
    // Compute integral of Planck function from wvn_lo to infinity
    ntn_bbd_blr[wvn_idx]=planck_integral_WiW76
      (wvn_lo, // [cm-1] Lower limit of Planck integral in wavenumbers
       tpt); // [K] Temperature

  } // !wvn_idx

  // Total integrated radiance resolved in wavenumber grid
  ntn_bbd_rsl_ttl=ntn_bbd_blr[0]-ntn_bbd_blr[wvn_nbr];

  // Fractional "missing" radiance bluer (shorter wavelengths/larger wavenumbers) than wavenumber grid
  flx_bbd_frc_mss_blr=ntn_bbd_blr[wvn_nbr]/ntn_bbd_ntg;

  // Fractional "missing" radiance redder (longer wavelengths/smaller wavenumbers) than wavenumber grid
  flx_bbd_frc_mss_rdr=(ntn_bbd_ntg-ntn_bbd_blr[0])/ntn_bbd_ntg;

  for(wvn_idx=0;wvn_idx<wvn_nbr;wvn_idx++){

    // Integrated radiance resolved in band is difference of consecutive half-open radiances to infinity
    ntn_bbd_rsl=ntn_bbd_blr[wvn_idx]-ntn_bbd_blr[wvn_idx+1];

    // Normalize [W m-2 sr-1] -> [frc] for fractional blackbody radiance/emission/flux
    if(nrm_typ_sng == "flx_frc_nrm") flx_bbd_frc[wvn_idx]=ntn_bbd_rsl/ntn_bbd_rsl_ttl; 
    else if(nrm_typ_sng == "flx_frc_abs") flx_bbd_frc[wvn_idx]=ntn_bbd_rsl/ntn_bbd_ntg;
    else err_prn(sbr_nm,"unknown normalization type "+nrm_typ_sng);

  } // !wvn_idx

  if(dbg_lvl_get() >= dbg_off){
    std::cout << "Diagnostics from " << sbr_nm << std::endl;
    std::cout << "Temperature = " << tpt << " K" << std::endl;
    std::cout << "Blackbody hemispheric irradiance = cst_stf_blt*T^4 = " << cst_Stefan_Boltzmann*std::pow(tpt,4.0) << " W m-2" << std::endl;
    std::cout << "Intensity of blackbody radiation = (cst_stf_blt*T^4)/pi = " << ntn_bbd_ntg << " W m-2 sr-1" << std::endl;
    std::cout << "Fractional missing radiance bluer (shorter wavelengths/larger wavenumbers) than wavenumber grid = " << flx_bbd_frc_mss_blr << std::endl;
    std::cout << "Fractional missing radiance redder (longer wavelengths/smaller wavenumbers) than wavenumber grid = " << flx_bbd_frc_mss_rdr << std::endl;
  
    std::cout << "idx\twvn_grd\t x_abc \tbbd_blr\tbbd_rdr" << std::endl;
    std::cout << "   \t cm-1  \t  frc  \tW/m2/sr\tW/m2/sr" << std::endl;
    for(wvn_idx=0;wvn_idx<=wvn_nbr;wvn_idx++) std::cout << wvn_idx << "\t" << wvn_grd[wvn_idx] << "\t" << cst_Planck*speed_of_light*wvn_grd[wvn_idx]*100/(cst_Boltzmann*tpt) << "\t" << ntn_bbd_blr[wvn_idx] << "\t" << ntn_bbd_ntg-ntn_bbd_blr[wvn_idx] << std::endl;

    std::cout << "idx\twvn_min\twvn_max\tflx_bbd_frc" << std::endl;
    std::cout << "   \t cm-1  \t cm-1  \t    frc    " << std::endl;
    for(wvn_idx=0;wvn_idx<wvn_nbr;wvn_idx++) std::cout << wvn_idx << "\t" << wvn_grd[wvn_idx] << "\t" << wvn_grd[wvn_idx+1] << "\t" << flx_bbd_frc[wvn_idx] << std::endl;
  } // !dbg

  delete []ntn_bbd_blr; // [W m-2 sr-1] Integrated radiance bluer than corresponding wavenumber

  return rcd; // [enm] Return success code
} // !spc_bbd_cls::flx_bbd_frc_get_WiW76()

double // [W m-2 sr-1] Integrated Planck function radiance from wvn_lo to infinity
planck_integral_WiW76 // [fnc] Compute integral of Planck function from wvn_lo to infinity
(double wvn_lo, // [cm-1] Lower limit of Planck integral in CGS wavenumbers
 double tpt){ // [K] Temperature
  /* Compute integral of Planck spectral radiance from wvn_lo [cm-1] to infinity
     Result in [W m-2 sr-1] is valid for [10 < wvn_lo < 10000 cm-1] at Earth's temperatures
     Result in [W m-2 sr-1] is valid for [0.2 < wvl_lo < 500 um] at Earth's temperatures
     Theory from Widger, W. K. and Woodall, M. P., Integration of the Planck blackbody radiation function, Bulletin of the Am. Meteorological Society, 57, 10, 1217-1219, Oct. 1976
     Based on C++ implementation from
     https://www.spectralcalc.com/blackbody/inband_radiance.html
     NB: For consistency with original sources, this function returns spectrally-integrated radiance
     Multiply result by pi to convert returned integrated radiance into spectrally-integrated hemispheric irradiance */

  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant (2018 SI NIST)
  using phc::cst_Planck; // (6.62606876e-34) [J s] Planck's constant (CODATA, 2018 SI NIST) (exact)
  using phc::speed_of_light; // (2.99792458e+08) [m s-1] Speed of light in vacuo (CODATA, 2018 SI NIST)

  const double spd_lgt_sqr(speed_of_light*speed_of_light); // [m2 s-2] Speed of light squared
  const double plk_spd_rcp_blt(cst_Planck*speed_of_light/cst_Boltzmann); // Radiation constant C1 from WiW76 (generally called c2 in other work)

  double x_abc=plk_spd_rcp_blt*100*wvn_lo/tpt; // [frc] Dimensionless abscissa (spectral coordinate) of integral
  // Compute re-used powers of x_abc, the dimensionless spectral coordinate
  double x_abc_sqr=x_abc*x_abc;
  double x_abc_cbd=x_abc*x_abc_sqr;
  
  /* How many terms of sum are needed?
     WiW76 shows that fewer than 200 terms suffice for NSD=10 in Earth's LW region [180 < T < 350 K], [3 < wvl < 100 um]
     Following WiW76 Figure 1, impose arbitrary maximum of 512 terms 
     20211123: 
     wvn=10cm-1 at 300 K = x_abc = 0.048 produces physically realistic (positive) bbd_rdr for itr_nbr=8, 16, 32 terms
     wvn=10cm-1 at 300 K = x_abc = 0.048 produces physically unrealistic (negative) bbd_rdr for itr_nbr=64, 128, 256, 512 terms (example shown below)
     Hence more terms are not a panacaea
     DISORT uses ~six terms in "power series" method, so set itr_nbr_max to 8
     Moreover DISORT switches from "power series" method to "exponential series" method for x (aka VCUT) < 1.5
     fxm: Implementing "exponential series" method for x_abc < ~1.5 seems most appropriate
     idx	wvn_grd	 x_abc 	bbd_blr	bbd_rdr
   	        cm-1  	  frc  	W/m2/sr	W/m2/sr
     0          10	0.048	146.199	-0.00068208 */
  const int itr_nbr_max(8); // [nbr] Maximum number of terms in power series method
  //std::string mth_sng("mth_WiW76_eqn6"); // [sng] Method of WiW76 Equation (6) (only)
  std::string mth_sng("mth_WiW76_appn"); // [sng] Method of WiW76 Appendix

  double itr_rcp_dbl; // [frc] Reciprocal of iteration index in double precision

  double itr_nbr_dbl=2.0+20.0/x_abc;
  itr_nbr_dbl=(itr_nbr_dbl < itr_nbr_max) ? itr_nbr_dbl : itr_nbr_max;
  int itr_nbr=int(itr_nbr_dbl);

  // Initialize series sum
  double srs_sum=0;

  if(mth_sng == "mth_WiW76_eqn6"){
    // Method of Equation 6 evaluates all four terms as series expansions
    for(int itr_idx=1;itr_idx<itr_nbr;itr_idx++){
      itr_rcp_dbl=1.0/itr_idx;
      srs_sum+=exp(-itr_idx*x_abc)*(x_abc_cbd+(3.0*x_abc_sqr+6.0*(x_abc+itr_rcp_dbl)*itr_rcp_dbl)*itr_rcp_dbl)*itr_rcp_dbl;
    } // !itr_idx
  }else if(mth_sng == "mth_WiW76_appn"){
    // Method in Appendix evaluates first term in Equantion (6) analytically, and remaining three terms as series expansions
    // Converges to 10 significant digits with 10-30% fewer terms than mth_WiW76_eqn6
    srs_sum+=-x_abc_cbd*log(1.0-exp(-x_abc));
    for(int itr_idx=1;itr_idx<itr_nbr;itr_idx++){
      itr_rcp_dbl=1.0/itr_idx;
      srs_sum+=exp(-itr_idx*x_abc)*(3.0*x_abc_sqr+6.0*(x_abc+itr_rcp_dbl)*itr_rcp_dbl)*itr_rcp_dbl*itr_rcp_dbl;
    } // !itr_idx
  } // !mth_sng
  
  // Return result in units of [W m2 sr-1]
  const double two_plk_spd_sqr(2.0*cst_Planck*spd_lgt_sqr); // [] Radiation constant C2 in WiW76 (often called c1L in other work)
  return two_plk_spd_sqr*pow(tpt/plk_spd_rcp_blt,4)*srs_sum;
} // !planck_integral_WiW76()
