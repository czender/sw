// $Id$ 

// Purpose: Description (definition) of blackbody spectrum classes

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <spc_bbd.hh> // Blackbody spectra

#ifndef SPC_BBD_HH // Contents have not yet been inserted in current source file
#define SPC_BBD_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers
#include <cfloat> // Floating point representation, FLT_MAX, DBL_EPSILON...
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc. (required by spc_bbd_prn())

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <wvl_grd.hh> // Wavelength grids
#include <mth_cst.hh> // Mathematical constants, cst_M_PIl...
#include <phys_cst.hh> // Physical constants NB: Live code, not a header

// Typedefs

// Forward declarations
class spc_bbd_cls; // Blackbody spectra

class spc_bbd_cls{ // Blackbody spectra
public:
  // Static public members

  // Public member functions
  spc_bbd_cls // [fnc] Default constructor
  (const prc_cmp tpt_arg); // [K] Blackbody temperature

  ~spc_bbd_cls(); // [fnc] Destructor

  prc_cmp eval(const prc_cmp &wvl); // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation

  prc_cmp // [frc] Fraction of blackbody emission in given spectral region
  flx_frc_get // [fnc] Fraction of blackbody emission in given spectral region
  (const prc_cmp &wvl_min, // [m] Minimum wavelength
   const prc_cmp &wvl_max); // [m] Maximum wavelength

  int // O [enm] Return success code
  flx_frc_get // [fnc] Fraction of blackbody emission in given spectral region
  (const prc_cmp *wvl_min, // I [m] Minimum wavelength
   const prc_cmp *wvl_max, // I [m] Maximum wavelength
   const long &wvl_nbr, // I [nbr] Number of wavelength bands
   prc_cmp *flx_IR_frc, // O [frc] Fraction of infrared flux in band
   const long &bnd_nbr_arg=bnd_nbr_dfl); // I [nbr] Number of bands to discretize Planck function

  prc_cmp flx_ttl()const; // [W m-2] Hemispheric blackbody irradiance

  prc_cmp tpt_get()const; // [K] Blackbody temperature
  int // [enm] Return success code
  tpt_set(const prc_cmp &tpt_arg); // [K] Blackbody temperature

private:
  // Static private members
  static int nst_nbr; // [nbr] Number of instantiated class members
  static const long bnd_nbr_dfl; // [nbr] Default number of bands to discretize Planck function

  // Private members
  int rcd; // [enm] Return success code
  double tpt; // [K] Blackbody temperature
  double xpn; // [frc] Exponent in Planck function
  double dnm; // [frc] Denominator of Planck function
  double ntn_bbd_wvl; // [W m-2 m-1 sr-1] Specific intensity of blackbody radiation
  double ntn_ntg; // [W m-2 sr-1] Intensity of blackbody radiation = (cst_stf_blt*T^4)/pi
  double flx_hms; // [W m-2] Hemispheric blackbody irradiance = cst_stf_blt*T^4
  long bnd_nbr; // [nbr] Number of bands to discretize Planck function

  // Private member functions
  int recompute(); // [fnc] Recompute properties of object

}; // end class spc_bbd_cls

// Declare functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

#endif // SPC_BBD_HH  






