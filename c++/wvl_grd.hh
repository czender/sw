// $Id$ 

// Purpose: Description (definition) of Wavelength grid classes

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <wvl_grd.hh> // Wavelength grids

#ifndef WVL_GRD_HH // Contents have not yet been inserted in current source file
#define WVL_GRD_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <phys_cst.hh> // Physical constants

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Define wvl_grd_cls class
class wvl_grd_cls{

public:

  // Friends
  friend std::ostream & // [srm] Reference to output stream for cascading
  operator<< // [fnc] Stream insertion operator
  (std::ostream &srm_out, // [srm] Output stream
   const wvl_grd_cls &wvl_obj); // [obj] Object to insert in stream
  
  // Static public members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  static std::string opt2abb(const std::string opt_sng); // Convert option to abbreviation 
  static const prc_cmp wvl_dbg_dfl; // [m] Debugging wavelength, default

  // Public member functions
  wvl_grd_cls // Default constructor
  (const std::string grd_sng_arg, // I [sng] Type of wavelength grid
   prc_cmp &wvl_mnm_arg, // I/O [m] Minimum wavelength
   prc_cmp &wvl_mxm_arg, // I/O [m] Maximum wavelength
   long &wvl_nbr_arg, // I/O [nbr] Number of wavelengths
   const bool flg_wvl_grd_arg=true, // I [flg] Compute wavelength grid
   const bool flg_wvn_grd_arg=true, // I [flg] Compute wavenumber grid
   const bool flg_frq_grd_arg=false); // I [flg] Compute frequency grid

  ~wvl_grd_cls(); // Destructor

  // fxm: implement typ() functions
  std::string typ()const; // [fnc] Type of wavelength grid
  void typ(const std::string sng); // [fnc] Type of wavelength grid

  long wvl_nbr_get()const; // [nbr] Number of wavelength bands
  int // [enm] Return success code
  wvl_nbr_set(const long &wvl_nbr_arg); // [nbr] Number of wavelength bands

  int // [enm] Return success code
  wvl_grd_min_max_set // [fnc] Reset wavelength grid boundaries
  (const prc_cmp &wvl_mnm_arg, // [m] Minimum wavelength
   const prc_cmp &wvl_mxm_arg); // [m] Maximum wavelength

  prc_cmp *wvl_ctr_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvl_dlt_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvl_grd_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvl_max_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvl_min_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvn_ctr_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvn_dlt_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvn_grd_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvn_max_get()const; // [fnc] Get pointer to const data
  prc_cmp *wvn_min_get()const; // [fnc] Get pointer to const data
  prc_cmp *frq_ctr_get()const; // [fnc] Get pointer to const data
  prc_cmp *frq_dlt_get()const; // [fnc] Get pointer to const data
  prc_cmp *frq_grd_get()const; // [fnc] Get pointer to const data
  prc_cmp *frq_max_get()const; // [fnc] Get pointer to const data
  prc_cmp *frq_min_get()const; // [fnc] Get pointer to const data

private:

  // Static private class members
  static int nst_nbr; // [nbr] Number of instantiated class members
  static const long wvl_nbr_max; // [nbr] Maximum number of wavelengths
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create abbreviation map
  static sng2sng_map opt2abb_map; // [fnc] Option to abbreviation map

  // Private members
  bool flg_dyn_mmr; // [flg] Dynamic memory arrays have been allocated
  bool flg_wvl_grd; // [flg] Compute wavelength grid
  bool flg_wvn_grd; // [flg] Compute wavenumber grid
  bool flg_frq_grd; // [flg] Compute frequency grid
  int rcd; // [enm] Return success code
  std::string grd_sng;
  long wvl_idx_dbg; // [idx] Debugging wavelength bin
  long wvl_nbr; // [nbr] Number of wavelengths
  prc_cmp wvl_mnm;
  prc_cmp wvl_mxm;
  // Other private members
  prc_cmp *wvl_ctr;
  prc_cmp *wvl_dlt;
  prc_cmp *wvl_grd;
  prc_cmp *wvl_max;
  prc_cmp *wvl_min;
  prc_cmp *wvn_ctr;
  prc_cmp *wvn_dlt;
  prc_cmp *wvn_grd;
  prc_cmp *wvn_max;
  prc_cmp *wvn_min;
  prc_cmp *frq_ctr;
  prc_cmp *frq_dlt;
  prc_cmp *frq_grd;
  prc_cmp *frq_max;
  prc_cmp *frq_min;
  prc_cmp wvn_mnm;
  prc_cmp wvn_mxm;
  long wvn_nbr;
  prc_cmp frq_mnm;
  prc_cmp frq_mxm;
  long frq_nbr;

  // Private member functions
  int allocate();
  int deallocate();
  int recompute();

}; // end class wvl_grd_cls

// Declare functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

#endif // WVL_GRD_HH  






