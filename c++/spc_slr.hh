// $Id$ 

// Purpose: Description (definition) of solar spectra classes

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <spc_slr.hh> // Solar spectra

#ifndef SPC_SLR_HH // Contents have not yet been inserted in current source file
#define SPC_SLR_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <sys/stat.h> // stat()
#include <cerrno> // perror()
#include <cassert> // Assertions

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
#include <dbg.hh> // Debugging constants
#include <fio.hh> // File Input/Output Layer
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <xtr.hh> // Extrapolation class
#include <vec.hh> // Vector functions ntp_vec(), rbn_vec()
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Forward declarations

// Typedefs
typedef prc_cmp (*flx_slr_frc_fnc_ptr_typ) // [fnc] Pointer to function returning fractional solar fluxes
(const prc_cmp &wvl_min, // [m] Minimum wavelength in band
 const prc_cmp &wvl_max); // [m] Maximum wavelength in band

// Solar flux sources
typedef struct{ // [sct] spc_slr_sct Solar flux source structure
  std::string abb; // [sng] Solar flux source abbreviation
  std::string dsc; // [sng] Solar flux source description
  std::string fl_slr_spc; // [sng] File containing solar spectrum
  flx_slr_frc_fnc_ptr_typ flx_slr_frc_fnc_ptr; // [fnc] Pointer to function returning fractional solar fluxes
} spc_slr_sct; // [sct] Solar flux source structure

typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map
typedef std::map<std::string,spc_slr_sct,std::less<std::string> > sng2spc_slr_sct_map; // String to spc_slr map

// Define spc_slr_cls class
class spc_slr_cls{

public:

  // Friends
  friend std::ostream & // [srm] Reference to output stream for cascading
  operator<< // [fnc] Stream insertion operator
  (std::ostream &srm_out, // [srm] Output stream
   const spc_slr_cls &spc_slr_obj); // [obj] Object to insert in stream

  // Static class members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  static std::string opt2abb(const std::string &opt_sng); // [sng] Convert option to abbreviation 
  static std::string abb2dsc(const std::string &abb_sng); // [sng] Convert abbreviation to description
  static std::string abb2fl_slr_spc(const std::string &abb_sng); // [sng] Convert abbreviation to file

  static flx_slr_frc_fnc_ptr_typ // [fnc] Pointer to function returning fractional solar fluxes
  abb2fnc // [fnc] Abbreviation to function mapper
  (const std::string &abb_sng); // [sng] Solar flux source abbreviation

  static int // O [enm] Return success code
  tst // [fnc] Self-test of spc_slr_cls class
  (long obj_nbr=100); // [nbr] Number of objects to create/destroy

  // Public member functions

  spc_slr_cls // [fnc] Constructor
  (const std::string &spc_slr_arg=spc_slr_typ_dfl, // [sng] Solar flux source abbreviation
   const prc_cmp &slr_cst_arg=slr_cst_dfl, // [W m-2] Solar constant
   const std::string &fl_slr_spc_arg=""); // [sng] File containing solar spectrum

  ~spc_slr_cls(); // [fnc] Destructor

  int // [enm] Return success code
  typ_set(const std::string &sng); // [sng] Set solar flux source type

  int // [enm] Return success code
  slr_cst_set(const prc_cmp &slr_cst_arg); // [W m-2] Solar constant

  int // [enm] Return success code
  fl_slr_spc_set(const std::string &fl_slr_spc_arg=""); // [sng] File containing solar spectrum

  int // [enm] Return success code
  flx_slr_frc_fnc_ptr_set(const flx_slr_frc_fnc_ptr_typ &flx_slr_frc_fnc_ptr_arg); // [fnc] Function to compute solar spectrum

  std::string dsc_get()const; // [sng] Solar flux description
  std::string typ_get()const; // [sng] Solar flux source type
  std::string fl_slr_spc_get()const; // [sng] File containing solar spectrum
  prc_cmp slr_cst_get()const; // [W m-2] Solar constant
  flx_slr_frc_fnc_ptr_typ flx_slr_frc_fnc_ptr_get()const; // [fnc] Function to compute solar spectrum

  prc_cmp // O [frc] Fraction of solar flux in band
  flx_frc_get // [fnc] Fraction of solar spectrum in given spectral region
  (const prc_cmp &wvl_min, // I [m] Minimum wavelength
   const prc_cmp &wvl_max)const; // I [m] Maximum wavelength
  // end flx_frc_get() prototype

  int // O [enm] Return success code
  flx_frc_get // [fnc] Fraction of solar spectrum in given spectral region
  (const prc_cmp *wvl_min, // I [m] Minimum wavelength
   const prc_cmp *wvl_max, // I [m] Maximum wavelength
   const long &wvl_nbr, // I [nbr] Number of wavelength bands
   prc_cmp *flx_slr_frc)const; // O [frc] Fraction of solar flux in band
  // end flx_frc_get() prototype

private:

  // Static private class members
  static int nst_nbr; // [nbr] Number of instantiated class members
  static prc_cmp slr_cst_dfl; // [W m-2] Solar constant, default
  static std::string spc_slr_typ_dfl; // [sng] Solar flux source abbreviation, default
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create option to abbreviation map
  static sng2sng_map opt2abb_map; // [map] Option to abbreviation map
  static sng2spc_slr_sct_map spc_slr_map_mk(); // [fnc] Create solar flux source map
  static sng2spc_slr_sct_map spc_slr_map; // [map] Solar flux source map

  // Private members

  /* Making rcd a class member is convenient but also problematic because 
     often rcd is only structure member a function modifies and if rcd were not 
     a member then entire function could be declared const.
     Hence the great compromise convention:
     Public get() functions always modify and return local rcd_lcl (not private member rcd)
     Private get() and Public and Private set() functions always modify and return private member rcd */

  bool flg_raw_spc_in_dyn_mmr; // [flg] Raw spectrum from file is in dynamic memory
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  flx_slr_frc_fnc_ptr_typ flx_slr_frc_fnc_ptr; // [fnc] Pointer to function returning fractional solar fluxes
  int rcd; // [enm] Return success code
  long wvl_nbr_in; // [nbr] Number of wavelength bands
  prc_cmp *flx_frc_blr_in; // [frc] Fraction of solar flux at shorter wavelengths
  prc_cmp *flx_slr_frc_in; // [frc] Fraction of solar flux in band
  prc_cmp *wvl_max_in; // [m] Maximum wavelength
  prc_cmp *wvl_min_in; // [m] Minimum wavelength
  prc_cmp slr_cst; // [W m-2] Solar constant
  std::string abb; // [sng] Abbreviation
  std::string dsc; // [sng] Description
  std::string fl_slr_spc; // [sng] File containing solar spectrum

  // Private member functions
  int allocate(); // [fnc] Allocate dynamic memory for object
  int deallocate(); // [fnc] Free dynamic memory for object
  int recompute(); // [fnc] Recompute properties of object

  int // O [enm] Return success code
  ntp_slr_flx_mnt // [fnc]  Interpolate raw solar spectrum to monotonic grid
  (const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
   const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
   const long wvl_nbr_out, // I [nbr] Number of wavelength bands
   prc_cmp *flx_frc_out, // O [frc] Fraction of solar flux in band
   const bool wrn_ntp_flg=false)const; // I [flg] Print WARNINGs from ntp_vec()
  // end ntp_slr_flx_mnt() prototype

  int // O [enm] Return success code
  ntp_slr_flx_nmn // [fnc] Interpolate raw solar spectrum to non-monotonic grid
  (const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
   const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
   const long wvl_nbr_out, // I [nbr] Number of wavelength bands
   prc_cmp *flx_frc_out)const; // O [frc] Fraction of solar flux in band
  // end ntp_slr_flx_nmn() prototype

  int // O [enm] Return success code
  ntp_slr_flx_nmn_CAM_SW // [fnc] Interpolate raw solar spectrum to CAM_SW grid
  (const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
   const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
   const long wvl_nbr_out, // I [nbr] Number of wavelength bands
   prc_cmp *flx_frc_out)const; // O [frc] Fraction of solar flux in band
  // end ntp_slr_flx_nmn_CAM_SW() prototype

  int // O [enm] Return success code
  spc_slr_raw_fl_get // [fnc] Load raw solar spectrum data from file
  (const std::string &fl_in, // I [sng] File containing solar spectrum
   prc_cmp *&wvl_min, // O [m] Minimum wavelength
   prc_cmp *&wvl_max, // O [m] Maximum wavelength
   long &wvl_nbr, // O [nbr] Number of wavelength bands
   prc_cmp *&flx_frc_blr, // O [frc] Fraction of solar flux at shorter wavelengths
   prc_cmp *&flx_slr_frc); // O [frc] Fraction of solar flux in band
  // end spc_slr_raw_fl_get() prototype

  void prn()const; // [fnc] Print object contents

}; // end class spc_slr_cls

// Prototype global functions with C++ linkages

prc_cmp // O [frc] Fraction of solar flux in band
flx_slr_frc_lsr // [fnc] Solar flux of laser
(const prc_cmp &wvl_min, // I [m] Minimum wavelength in band
 const prc_cmp &wvl_max); // I [m] Maximum wavelength in band
// end flx_slr_frc_lsr() prototype

int // O [enm] Return success code
flx_slr_frc_lsr // [fnc] Solar flux of laser
(const prc_cmp *wvl_min, // I [m] Minimum wavelength in band
 const prc_cmp *wvl_max, // I [m] Maximum wavelength in band
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 prc_cmp *flx_slr_frc); // O [frc] Fraction of solar flux in band
// end flx_slr_frc_lsr() prototype

// Define inline'd functions in header so source is visible to calling files

#endif // SPC_SLR_HH  






