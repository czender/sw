// $Id$ 

// Purpose: Boundary layer meteorology

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <blm.hh> // Boundary layer meteorology

#ifndef BLM_HH // Contents have not yet been inserted in current source file  
#define BLM_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <valarray> // STL valarray class template
#include <vector> // STL vector class template

// Standard C headers 
#include <cstdio> // stderr, FILE, NULL, etc.
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <a2d.hh> // Two dimensional arrays
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <phys_cst.hh> // Physical constants NB: Live code, not a header
#include <mth.hh> // Mathematical utilities, constants
#include <lsm.hh> // Land Surface Model
#include <tdy.hh> // Atmospheric thermodynamics
#include <cln.hh> // Calendar algorithms

namespace blm{ // [nms] Boundary layer meteorology namespace
  const prc_cmp rgh_mmn_ice_lak(0.04); // (0.04) [m] Roughness length over frozen lakes Bon96 p. 59
  const prc_cmp rgh_mmn_ice_lnd(0.05); // (0.05) [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
  const prc_cmp rgh_mmn_ice_ocn(0.0005); // (0.0005) [m] Roughness length over sea ice BKL97 p. F-3 (updated)
  //   const prc_cmp rgh_mmn_ice_ocn(0.04); // (0.04) [m] Roughness length over sea ice BKL97 p. F-3 CCM:dom/parpbl.h fxm
  //  const prc_cmp rgh_mmn_ice_ocn(0.05); // (0.05) [m] Roughness length over sea ice BKL97 p. F-3 (original) fxm
  const prc_cmp rgh_mmn_lak_wrm(0.001); // (0.001) [m] Roughness length over unfrozen lakes Bon96 p. 59
  const prc_cmp rgh_mmn_snw(0.04); // (0.04) [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
  const prc_cmp wnd_min_dps(1.0); // (1.0) [m s-1] Minimum windspeed used for deposition Bon96 p. 51
} // end Boundary layer meteorology namespace blm

 // Declare functions that have C++ linkages
std::string oro_sng_get(const prc_cmp oro);
std::string soi_typ_sng_get(const long soi_typ_idx);
std::string pnt_typ_sng_get(const long pnt_typ_idx);
std::string pft_sng_get(const long pft_idx);
std::string sfc_typ_sng_get(const long sfc_typ_idx);
std::string sfc_typ_dsc_get(const long sfc_typ_idx);

int // O [rcd] Return success code
lnd_frc_mbl_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp doy, // I [day] Day of year [1.0..367.0)
 const prc_cmp lat_rdn, // I [rdn] Latitude
 const prc_cmp *lnd_frc_dry, // I [frc] Dry land fraction
 prc_cmp *lnd_frc_mbl, // O [frc] Bare ground fraction
 const prc_cmp *oro, // I [frc] Orography
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc); // I [frc] Fraction of surface covered by snow

int // O [rcd] Return success code
wnd_rfr_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_mdp, // I [m] Midpoint height above surface
 const prc_cmp hgt_rfr, // I [m] Reference height for dust mobilization processes
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_mdp, // I [m s-1] Surface layer mean wind speed
 prc_cmp *wnd_rfr); // O [m s-1] Wind speed at reference height

int // O [rcd] Return success code
snw_frc_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *snw_frc, // O [frc] Fraction of surface covered by snow
 prc_cmp *snw_hgt, // O [m] Geometric bulk thickness of snow
 const prc_cmp *snw_hgt_lqd); // I [m] Equivalent liquid water snow depth

int // O [rcd] Return success code
rgh_zpd_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *oro,  // I [frc] Orography
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *wnd_10m, // I [m s-1] Reference height (i.e., 10 m) wind speed used for dust deposition
 prc_cmp *hgt_zpd); // O [m] Zero plane displacement

int // O [rcd] Return success code
mmn_dlt_evl
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_lwr, // I [m] Height of lower level
 const prc_cmp *hgt_upr, // I [m] Height of upper level
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 prc_cmp *mmn_dlt, // O [m s-1] Upper level - lower level windspeed
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *wnd_frc); // I [m s-1] Surface friction velocity

int // O [rcd] Return success code
rgh_mmn_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *oro,  // I [frc] Orography
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *wnd_10m); // I [m s-1] Reference height (i.e., 10 m) wind speed used for dust deposition

int // O [rcd] Return success code
rss_aer_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement height
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 prc_cmp *rss_aer, // O [s m-1] Aerodynamic resistance
 const prc_cmp *wnd_frc); // I [m s-1] Surface friction velocity

int // O [rcd] Return success code
zpd_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *hgt_zpd, // O [m] Zero plane displacement
 const prc_cmp *oro,  // I [frc] Orography
 const long *sfc_typ); // I [idx] LSM surface type (0..28)

int // O [rcd] Return success code
lnd_is_vgt
(const long lon_nbr, // I [nbr] Size of arrays
 const long *pnt_typ_idx, // I [idx] Plant type index 
 const prc_cmp *snw_hgt, // I [m] Geometric bulk thickness of snow
 bool *vgt); // O [flg] "Vegetated" flag

int // O [rcd] Return success code
trn_fsh_vpr_soi_atm_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp tpt_soi_frz, // I [K] Temperature of frozen soil
 prc_cmp *trn_fsh_vpr_soi_atm, // O [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
 const prc_cmp *vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
 const prc_cmp *vwc_sfc); // I [m3 m-3] Volumetric water content

// Define inline'd functions in header so source is visible to calling files
inline bool // O [flg] True if > 50% sea ice, false otherwise
oro_is_ice(const prc_cmp oro)
{
  // Purpose: Return true if > 50% sea ice, false otherwise
  // NB: Use floor() since C nint() is not standard (e.g., <sunmath.h>)
  return (std::floor(oro+0.5) == 2 ? true : false);
} // end oro_is_ice()

inline bool // O [flg] True if > 50% land, false otherwise
oro_is_lnd(const prc_cmp oro)
{
  // Purpose: Return true if > 50% land, false otherwise
  // NB: Use floor() since C nint() is not standard (e.g., <sunmath.h>)
  return (std::floor(oro+0.5) == 1 ? true : false);
} // end oro_is_lnd()

inline bool // O [flg] True if > 50% ocn, false otherwise
oro_is_ocn(const prc_cmp oro)
{
  // Purpose: Return true if > 50% ocean, false otherwise
  // NB: Use floor() since C nint() is not standard (e.g., <sunmath.h>)
  return (std::floor(oro+0.5) == 0 ? true : false);
} // end oro_is_ocn()

inline prc_cmp // O [frc] Neutral 10 m drag coefficient over ocean
xch_cff_mmn_ocn_ntr_get(prc_cmp wnd_10m_ntr)
{
  /* Purpose: Return neutral 10 m drag coefficient over ocean
     Input is neutral wind speed at 10 m from measurements or from model
     fxm: What is result of inputting ambient (non-neutral) 10 m wind speed here?
     LaP81 parameterization is 1.2e-3 for 4 < U10 < 11, 0.00049+0.000065*U10 for 11 < U10 < 25 m s-1
     LaP82 parameterization is 1.14e-3 for 3 < U10 < 10, 0.00049 + 0.000065*U10 for 10 < U10 < 25 m s-1
     LaP81 show CDN is independent of U10 for U10 < 10 m s-1, but significant for U10 > 10 m s-1
     W. Large uses following parameterization (presumably an updated version of LaP8X) in CSM */
  return 0.0027/wnd_10m_ntr+0.000142+0.0000764*wnd_10m_ntr; // [frc] CCM:dom/flxoce(), NOS97 p. I2
} // end xch_cff_mmn_ocn_ntr_get()

inline prc_cmp // O [frc] Stability correction for momentum, stable case
mno_stb_crc_mmn_stb_get(const prc_cmp mno_stb_prm){
  /* Purpose: Given the Monin-Obukhov stability parameter (usually called zeta) for 
     a stable atmosphere, return the stability correction factor for momentum, 
     usually called psi
     References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
     Lap82 p. 466 show better fits using 7.0 instead of 5.0
     Most authors assume the stability correction factors for heat and momentum are equal in stable (but not unstable) conditions
     Currently this function is BFB with CCM:dom/flxoce() */
  return -5.0*mno_stb_prm; // [frc]
} // end mno_stb_crc_mmn_stb_get()

inline prc_cmp // O [frc] Stability correction for momentum, unstable case
mno_stb_crc_mmn_uns_get(const prc_cmp sml_fnc_mmn_uns_rcp){
  /* Purpose: Given the reciprocal of the Monin-Obukhov similarity function 
     (usually called phi) for momentum in an unstable atmosphere, return the 
     stability correction factor for momentum, usually called psi 
     References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
     Currently this function is BFB with CCM:dom/flxoce() */
  return std::log((1.0+sml_fnc_mmn_uns_rcp*(2.0+sml_fnc_mmn_uns_rcp))*(1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/8.0)-2.0*std::atan(sml_fnc_mmn_uns_rcp)+1.571; // [frc]
} // end mno_stb_crc_mmn_uns_get()

inline prc_cmp // O [frc] Stability correction for heat, unstable case
mno_stb_crc_heat_uns_get(const prc_cmp sml_fnc_mmn_uns_rcp){
  /* Purpose: Given the reciprocal of the Monin-Obukhov similarity function 
     (usually called phi) for momentum in an unstable atmosphere, return the 
     stability correction factor for heat, usually called psi
     References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
     Currently this function is BFB with CCM:dom/flxoce() */
  return 2.0*std::log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0); // [frc]
} // end mno_stb_crc_heat_uns_get()

inline prc_cmp // O [frc] Stability correction for momentum
mno_stb_crc_mmn_get(const prc_cmp mno_stb_prm){
  /* Purpose: Given the Monin-Obukhov stability parameter z/L (usually called zeta),
     return the stability correction factor for momentum, usually called psi   
     References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
     Currently this function is BFB with CCM:dom/flxoce() */
  
  // Local
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  prc_cmp stb_val; // [flg] 1.0 if stable, 0.0 if unstable
  
  // Main Code
  stb_val=PRC_CMP(0.5)+sign_cpv(PRC_CMP(0.5),mno_stb_prm); // [flg] 1.0 if stable, 0.0 if unstable
  sml_fnc_mmn_uns_rcp=max_cpv(std::sqrt(PRC_CMP_ABS(PRC_CMP(1.0)-16.0*mno_stb_prm)),PRC_CMP(1.0)); // [frc] BKL97 p. F1
  sml_fnc_mmn_uns_rcp=std::sqrt(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
  
  return -PRC_CMP(5.0)*mno_stb_prm*stb_val+(PRC_CMP(1.0)-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp); // [frc] BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
  
} // end mno_stb_crc_mmn_uns_get()

#endif // BLM_HH  
