// $Id$ 

// Purpose: Description (definition) of line-by-line utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <lbl.hh> // Line-by-line utilities

#ifndef LBL_HH // Contents have not yet been inserted in current source file
#define LBL_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <htrn_c++.hh> // HITRAN line database definitions
#include <mth.hh> // Mathematical utilities, constants
#include <wvl_grd.hh> // Wavelength grids
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Forward declarations

// Typedefs
typedef struct{ // [sct] ln_sct Structure for transition lines
  short ln_iso; // [enm] HITRAN isotope number (1..9)
  short ln_mlc; // [enm] HITRAN molecule number (1..37)
} ln_sct; // [sct] Structure for transition lines

// Define lbl_cls class

// Prototype global functions with C++ linkages
int // O [enm] Return success code
rt_lbl // [fnc] Single line radiative transfer
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &prs_dlt, // I [Pa] Pressure thickness
 const prc_cmp &prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp &slr_zen_ngl_cos, // I [frc] Cosine solar zenith angle
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp &vmr_gas, // I [mlc mlc-1] Volume mixing ratio of gas
 const std::string &lbl_tst, // I [sng] Name of line-by-line test
 const wvl_grd_cls &wvlgrd); // I [m] Wavelength grid
// end rt_lbl() prototype

int // O [enm] Return success code
lnshp_lrn_vct // [fnc] Lorentzian line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_lrn_val); // O [Hz-1,cm] Lorentzian line shape profile
// end lnshp_lrn_vct() prototype

int // O [enm] Return success code
lnshp_dpp_vct // [fnc] Doppler line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWEM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_dpp_val); // O [Hz-1,cm] Doppler line shape profile
// end lnshp_dpp_vct() prototype

int // O [enm] Return success code
lnshp_vgt_vct // [fnc] Doppler line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWEM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const double &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_vgt_val); // O [Hz-1,cm] Doppler line shape profile
// end lnshp_vgt_vct() prototype

// Define inline'd functions in header so source is visible to calling files

template<class prc_T>prc_T // O [frc] Complex error function
err_fnc_cpx_Hum82 // [fnc] Complex error function of Humlicek (1982)
(prc_T zzz); // [frc] Argument to complex error function
// end err_fnc_cpx_Hum82() prototype

template<class prc_T>prc_T // O [frc] Complex error function
err_fnc_cpx_Hum82 // [fnc] Complex error function of Humlicek (1982)
(prc_T zzz) // [frc] Argument to complex error function
{
  /* Purpose: Evaluate complex probability function = complex error function 
     using method of Humlicek (1982) Hum82.
     Voigt function is real component of the complex error function
     See HAW78, Hum82, Kun95 have further details
     See Lio92 p. 30 (2.2.11) derives arguments connecting Voigt function to w(z)
     In Hum82 notation, complex error function is w(z)=u(x,y)+iv(x,y) 
     Voigt function is u(x,y)= Hum82 = K(x,y) = Lio92 = 
     \frac{\yyy}{\pi} \int_{-\infty}^{+\infty} 
     \frac{\me^{\ttt^2}}{\yyy^{2}+(\xxx-ttt)^{2}} \,\dft\ttt

     The complex error function is also related to, but should not be confused with,
     the complementary error function erfc(z), through the simple formula
     erfc(-iz) = exp(t^2) w(z)
     erfc(z) = exp(-t^2) w(iz)

     When I received this routine (from Martin Kuntz?) it had this header,
     which quantifies the accuracy of the routine:
     "Computes the complex probability function w(z)=exp(-z^2)*erfc(-i*z)
     in the upper half-plane z=x+iy (i.e., for y >= 0).
     Maximum relative error of both real and imaginary parts is < 1.e-4" */

  // Local
  const std::string sbr_nm("err_fnc_cpx_Hum82"); // [sng] Name of subroutine
  const double x(zzz.real()); // [frc]
  const double y(zzz.imag()); // [frc]
  const prc_T ttt(y,-x); // [frc]
  const double s(std::fabs(x)+y); // [frc]
  prc_T u; // [frc]
  prc_T vgt_Hum82; // [frc]

  // Main code
  // Region 1
  if(s<15.0) goto rgn_2;
  vgt_Hum82=ttt*0.5641896/(0.5+ttt*ttt);
  return vgt_Hum82; // O [frc] Complex error function

 rgn_2: // Region 2
  if(s<5.5) goto rgn_3;
  u=ttt*ttt;
  vgt_Hum82=ttt*(1.410474+u*0.5641896)/(0.75+u*(3.0+u));
  return vgt_Hum82; // O [frc] Complex error function

 rgn_3:  // Region 3
  if(y<0.195*std::fabs(x)-0.176) goto rgn_4;
  vgt_Hum82=(16.4955+ttt*(20.20933+ttt*(11.96482+ttt*(3.778987+ttt*0.5642236))))/
    (16.4955+ttt*(38.82363+ttt*(39.27121+ttt*(21.69274+ttt*(6.699398+ttt)))));
  return vgt_Hum82; // O [frc] Complex error function
				 
 rgn_4: // Region 4
  u=ttt*ttt;
  vgt_Hum82=std::exp(u)-ttt*
    (36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419))))))/
    (32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))));
  return vgt_Hum82;
} // end err_fnc_cpx_Hum82()

template<class prc_T>prc_T // O [Hz-1,cm] Lorentzian line shape profile
lnshp_lrn_fst_scl // [fnc] Lorentzian line shape profile
(const prc_T &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const prc_T &frq_dlt) // I [Hz,cm-1] Distance from line center
{
  /* Purpose: "Fast scalar" Lorentzian line shape profile
     In order to maximize chances of compiler inlining this function,
     it avoids error checking.
     Input must be HWHM and frequency offset in Hz or in wavenumber
     Output is line shape profile in same units as input */
  prc_T lnshp_lrn; // O [Hz-1,cm] Lorentzian line shape profile
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  lnshp_lrn=HWHM_lrn/(mth::cst_M_PIl*(frq_dlt*frq_dlt+HWHM_lrn*HWHM_lrn)); // O [Hz-1,cm] Lorentzian line shape profile
  // Main code
  return lnshp_lrn; // O [Hz-1,cm] Lorentzian line shape profile
} // end lnshp_lrn_fst_scl()

template<class prc_T>prc_T // O [Hz-1,cm] Doppler line shape profile
lnshp_dpp_fst_scl // [fnc] Doppler line shape profile
(const prc_T &HWEM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const prc_T &frq_dlt) // I [Hz,cm-1] Distance from line center
{
  /* Purpose: "Fast scalar" Doppler line shape profile
     In order to maximize chances of compiler inlining this function,
     it avoids error checking.
     Input must be Doppler width (HWEM not HWHM) and frequency offset in Hz or in wavenumber
     Output is line shape profile in same units as input */
  prc_T lnshp_dpp; // O [Hz-1,cm] Doppler line shape profile
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  lnshp_dpp=(1.0/(HWEM_dpp*std::sqrt(mth::cst_M_PIl)))*std::exp(-frq_dlt*frq_dlt/(HWEM_dpp*HWEM_dpp)); // O [Hz-1,cm] Doppler line shape profile
  // Main code
  return lnshp_dpp; // O [Hz-1,cm] Doppler line shape profile
} // end lnshp_dpp_fst_scl()

template<class prc_T>prc_T // O [Hz-1,cm] Voigt line shape profile
lnshp_vgt_fst_scl // [fnc] Voigt line shape profile
(const prc_T &HWEM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const prc_T &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const prc_T &frq_dlt) // I [Hz,cm-1] Distance from line center
{
  /* Purpose: "Fast scalar" Voigt line shape profile
     In order to maximize chances of compiler inlining this function,
     it avoids error checking.
     Input is Doppler width (HWEM not HWHM), Lorentz HWHM and frequency offset in Hz or in wavenumber
     Output is line shape profile in same units as input */

  // Numerical integration by method of Humlicek (1982)
  // fxm: this "fast scalar" code is buggy, off by factor near 2, don't know why, array version works
  const double sqrt_ln_2(std::sqrt(std::log(2.0))); // [Hz,cm-1] Doppler half-width at half-maximum
  const double HWHM_dpp(sqrt_ln_2*HWEM_dpp); // [Hz,cm-1] Doppler half-width at half-maximum
  const double vgt_pre_fct(1.0/HWHM_dpp); // [Hz-1,cm] Voigt function pre-factor
  const double arg_rl_fct(sqrt_ln_2/HWHM_dpp); // [Hz-1,cm] Factor in real part of complex error function argument
  const double err_fnc_cpx_arg_img(sqrt_ln_2*HWHM_lrn/HWHM_dpp); // [frc] Complex error function imaginary argument
  const double err_fnc_cpx_arg_rl(arg_rl_fct*frq_dlt); // [frc] Complex error function real argument
  const std::complex<double> err_fnc_cpx_arg=std::complex<double>(err_fnc_cpx_arg_rl,err_fnc_cpx_arg_img); // [frc] Complex error function argument
  const std::complex<double> err_fnc_cpx_val(err_fnc_cpx_Hum82(err_fnc_cpx_arg)); // [frc] Complex error function
  const prc_T lnshp_vgt_val=vgt_pre_fct*err_fnc_cpx_val.real(); // O [Hz-1,cm] Voigt line shape profile

  // Analytic approximation of fxm: who as reported in Lio92 p. 31 (2.2.17)
  /*  const prc_T HWHM_vgt= // [Hz,cm-1] "Half-width" of Voigt profile Lio92 p. 31 (2.2.17)
    0.5*(HWHM_lrn+std::sqrt(HWHM_lrn*HWHM_lrn+4.0*std::log(2.0)*HWEM_dpp*HWEM_dpp))+
    0.05*HWHM_lrn*(1.0-2.0*HWHM_lrn/std::sqrt(HWHM_lrn*HWHM_lrn+4.0*std::log(2.0)*HWEM_dpp*HWEM_dpp));
  const prc_T eta(frq_dlt/HWHM_vgt); // [frc]
  const prc_T xi(HWHM_lrn/HWHM_vgt); // [frc]
  const prc_T lnshp_vgt_val= // O [Hz-1,cm] Voigt line shape profile Lio92 p. 31 (2.2.16)
    std::sqrt(std::log(2.0)/mth::cst_M_PIl)*(1.0/HWHM_vgt)*(1.0-xi)*std::exp(-std::log(2.0)*eta*eta)+
    (1.0/(mth::cst_M_PIl*HWHM_vgt))*xi*(1.0/(1.0+eta*eta))-
    (1.0/(mth::cst_M_PIl*HWHM_vgt)*xi*(1.0-xi)*(xi+1.0+1.5/std::log(2.0)))*
    (0.066*std::exp(-0.4*eta*eta)-1.0/(40.0-5.5*eta*eta+eta*eta*eta*eta));  */

  return lnshp_vgt_val; // O [Hz-1,cm] Voigt line shape profile
} // end lnshp_vgt_fst_scl()

#endif // LBL_HH
