// $Id$ 

// Purpose: Description (definition) of particle size distribution classes

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <psd.hh> // Particle size distributions

#ifndef PSD_HH // Contents have not yet been inserted in current source file
#define PSD_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc. (required by psd_prn())

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <sz_grd.hh> // Size grids
#include <mth_cst.hh> // Mathematical constants, cst_M_PIl...

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Forward declarations

// Define psd_cls class
class psd_cls{ // Particle size distributions
public:
  // Static public members
  static std::string opt2abb(const std::string &opt_sng); // [fnc] Create option to abbreviation map
  static int nst_nbr_get(); // [nbr] Number of instantiated class members

  // Public member functions
  psd_cls // Default constructor
  (const prc_cmp &rds_nma_arg=1.0e-6, // [m] Number median radius analytic
   const prc_cmp &gsd_anl_arg=2.0, // [frc] Geometric standard deviation
   const prc_cmp &cnc_nbr_anl_arg=1.0, // [# m-3] Number concentration analytic
   const prc_cmp &dns_prt_arg=1000.0, // [kg m-3] Particle density
   const std::string &dst_sng_arg="lognormal", // [sng] Distribution string
   sz_grd_cls *sz_grd_arg=0); // [sct] Size grid

  ~psd_cls(); // [fnc] Destructor

  int // [enm] Return success code
  mss_frc_anl_set // [fnc] Set mass fraction analytic for each mode of ensemble
  (psd_cls *psd_lst, // [sct] List of aerosol size distributions
   const long &psd_nbr); // [nbr] Number of particle modes

  int // [enm] Return success code
  mss_frc_anl_set(const prc_cmp &mss_frc_anl_arg); // [frc] Mass fraction analytic

  int // [enm] Return success code
  psd_prn // [fnc] Print formatted list of properties of particle size distribution(s)
  (const psd_cls *psd_lst, // [sct] List of particle size distributions
   const long &psd_nbr); // [nbr] Number of particle modes

  int // [enm] Return success code
  prs_psd_sng // [fnc] Set size distribution parameters from string
  (const std::string &psd_sng); // [sng] Triplet specifying size distribution

  void typ_set(const std::string &sng); // [sng] Type of size distribution
  std::string typ_get()const; // [sng] Type of size distribution

  int // [enm] Return success code
  sz_grd_set(sz_grd_cls *grd); // [sct] Size grid

  int // [enm] Return success code
  rds_nma_set(const prc_cmp &rds_nma_arg); // [m] Number median radius analytic

  int // [enm] Return success code
  dmt_nma_set(const prc_cmp &dmt_nma_arg); // [m] Number median diameter analytic

  int // [enm] Return success code
  gsd_anl_set(const prc_cmp &gsd_anl_arg); // [frc] Geometric standard deviation

  int // [enm] Return success code
  mss_anl_set(const prc_cmp &mss_anl_arg); // [kg m-3] Mass concentration analytic

  int // [enm] Return success code
  dns_prt_set(const prc_cmp &dns_prt_arg); // [kg m-3] Particle density

  int // [enm] Return success code
  cnc_nbr_anl_set(const prc_cmp &cnc_nbr_anl); // [# m-3] Number concentration analytic

  prc_cmp cnc_nbr_anl_get()const; // [# m-3] Number concentration analytic
  prc_cmp rds_nma_get()const; // [m] Number median radius analytic
  prc_cmp gsd_anl_get()const; // [frc] Geometric standard deviation
  prc_cmp dmt_nma_get()const; // [m] Number median diameter analytic

  prc_cmp mss_anl_get()const; // [kg m-3] Mass concentration analytic
  prc_cmp dns_prt_get()const; // [kg m-3] Particle density
  prc_cmp mss_frc_anl_get()const; // [frc] Mass fraction analytic
private:
  // Private members
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  bool usr_mss_flg; // [flg] Mass flag
  bool mss_frc_anl_set_flg; // [flg] Mass fraction analytic has been set
  prc_cmp mss_frc_anl; // [frc] Mass fraction analytic
  int rcd; // [enm] Return success code
  std::string dst_sng;
  prc_cmp rds_nma; // [m] Number median radius analytic
  prc_cmp dmt_nma; // [m] Number median diameter analytic
  prc_cmp gsd_anl; // [frc] Geometric standard deviation 
  prc_cmp mss_anl; // [kg m-3] Mass concentration analytic
  prc_cmp dns_prt; // [kg m-3] Particle density
  prc_cmp cnc_nbr_anl; // [# m-3] Number concentration analytic
  long sz_nbr; // [nbr] Number of elements of size grid
  sz_grd_cls *sz_grd; // Size grid

  //Resolved # between r and r+dr per unit r
  prc_cmp *dst; // [# m-3 m-1] Number distribution 
  prc_cmp *cnc; // [# m-3] Number concentration
  prc_cmp *xsa; // [m2 m-3] Cross-sectional area concentration
  prc_cmp *mss; // [kg m-3] Mass concentration
  prc_cmp *sfc; // [m2 m-3] Surface area concentration
  prc_cmp *vlm; // [m3 m-3] Volume concentration
  // Private member functions
  int allocate(); // [fnc] Allocate dynamic memory for object
  int deallocate(); // [fnc] Free dynamic memory for object
  int recompute(); // [fnc] Recompute properties of object
  // Static private members
  static int nst_nbr; // Number of instantiated class members
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create abbreviation map
  static sng2sng_map opt2abb_map; // [fnc] Option to abbreviation map
}; // end class psd_cls

// Declare functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

#endif // PSD_HH  






