// $Id$ 

// Purpose: Description (definition) of radiative transfer classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <rt.hh> // Radiative transfer

#ifndef RT_HH // Contents have not yet been inserted in current source file
#define RT_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cerrno> /* perror() */

// 3rd party vendors
#include <netcdf.h> // netCDF C interface
#include <gsl/gsl_errno.h> // GNU Scientific Library error handling
#include <gsl/gsl_sf_expint.h> // GNU Scientific Library special functions Exponential integrals

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth_cst.hh> // Mathematical constants, cst_M_PIl...
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <flp.hh> // Floating point utilities, constants
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Forward declarations

// Typedefs
typedef 
int // O [enm] Return success code
(*rt_sln_fnc_ptr_typ) // [fnc] Pointer to function returning layer optical properties
  (const long &wvl_nbr, // I [nbr] Number of wavelength bands
   const prc_cmp &plrmu_avg, // I [frc] Mean zenith angle
   const prc_cmp *asm_prm, // I [frc] Asymmetry parameter
   const prc_cmp *ss_alb, // I [frc] Single scattering albedo
   const prc_cmp *tau_ext, // I [frc] Extinction optical depth
   prc_cmp *rfl_flx, // O [frc] Flux reflectance
   prc_cmp *trn_flx, // O [frc] Flux transmittance
   prc_cmp *abs_flx); // O [frc] Flux absorptance
// end rt_sln_fnc_ptr_typ() prototype

// Radiative transfer methods
typedef struct{ // [sct] rt_sct Radiative transfer method structure
  std::string abb; // [sng] Radiative transfer method abbreviation
  std::string dsc; // [sng] Radiative transfer method description
  rt_sln_fnc_ptr_typ rt_sln_fnc_ptr; // [fnc] Pointer to function returning layer optical properties
} rt_sct; // [sct] Radiative transfer method structure

typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map
typedef std::map<std::string,rt_sct,std::less<std::string> > sng2rt_sct_map; // String-to-rt_sct map
typedef std::map<std::string,var_mtd_sct,std::less<std::string> > sng2var_mtd_sct_map; // String-to-var_mtd_sct map

// Define rt_cls class
class rt_cls{
  
public:
  
  // Friends
  friend std::ostream & // [srm] Reference to output stream for cascading
  operator<< // [fnc] Stream insertion operator
  (std::ostream &srm_out, // [srm] Output stream
   const rt_cls &rt_obj); // [obj] Object to insert in stream
  
  static int wrn_nbr; // [nbr] Number of accumulated warnings

  // Static class members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  static std::string opt2abb(const std::string &opt_sng); // [sng] Convert option to abbreviation 
  static std::string abb2dsc(const std::string &abb_sng); // [sng] Convert abbreviation to description
  
  static rt_sln_fnc_ptr_typ // [fnc] Pointer to function returning layer optical properties
  abb2fnc // [fnc] Abbreviation to function mapper
  (const std::string &abb_sng); // [sng] Radiative transfer method abbreviation
  
  static int // O [enm] Return success code
  tst // [fnc] Self-test of rt_cls class
  (long obj_nbr=100); // [nbr] Number of objects to create/destroy
  
  // Public member functions
  
  rt_cls // [fnc] Constructor
  (const std::string &rt_arg=rt_typ_dfl, // I [sng] Radiative transfer method abbreviation
   const long &wvl_nbr_arg=wvl_nbr_dfl, // I [nbr] Number of wavelength bands
   const prc_cmp *asm_prm_arg=(prc_cmp *)NULL, // I [frc] Asymmetry parameter
   const prc_cmp *ss_alb_arg=(prc_cmp *)NULL, // I [frc] Single scattering albedo
   const prc_cmp *tau_ext_arg=(prc_cmp *)NULL); // I [frc] Extinction optical depth
  
  ~rt_cls(); // [fnc] Destructor
  
  int // [enm] Return success code
  typ_set(const std::string &sng); // [sng] Set radiative transfer method
  
  int // [enm] Return success code
  wvl_nbr_set(const long &wvl_nbr_arg); // [nbr] Number of wavelength bands
  
  int // [enm] Return success code
  plrmu_avg_set(const prc_cmp &plrmu_avg_arg); // [frc] Mean zenith angle
  
  int // [enm] Return success code
  nc_out_set(const std::string &fl_out); // [fl] netCDF output file
  
  int // [enm] Return success code
  nc_out_set(const int &nc_out); // [fl] netCDF output file

  int // [enm] Return success code
  rt_sln_fnc_ptr_set(const rt_sln_fnc_ptr_typ &rt_sln_fnc_ptr_arg); // [fnc] Function to compute layer optical properties
  
  const std::string dsc_get()const; // [sng] Radiative transfer method description
  const std::string typ_get()const; // [sng] Radiative transfer method abbreviation
  prc_cmp plrmu_avg_get()const; // [frc] Mean zenith angle
  rt_sln_fnc_ptr_typ rt_sln_fnc_ptr_get()const; // [fnc] Function to compute layer optical properties
  
  long wvl_nbr_get()const; // [nbr] Number of wavelength bands
  const prc_cmp *rfl_flx_get()const; // [frc] Flux reflectance
  const prc_cmp *tau_ext_get()const; // [frc] Extinction optical depth
  const prc_cmp *asm_prm_get()const; // [frc] Asymmetry parameter
  const prc_cmp *ss_alb_get()const; // [frc] Single scattering albedo
  const prc_cmp *abs_flx_get()const; // [frc] Flux absorptance
  const prc_cmp *trn_flx_get()const; // [frc] Flux transmittance
  
  int // O [enm] Return success code
  twg_set // [fnc] Set slab single scattering properties
  (const long &wvl_nbr_arg, // I [nbr] Number of wavelength bands
   const prc_cmp *asm_prm_arg, // I [frc] Asymmetry parameter
   const prc_cmp *ss_alb_arg, // I [frc] Single scattering albedo
   const prc_cmp *tau_ext_arg); // I [frc] Extinction optical depth
  // end twg_set() prototype
  
  int // O [enm] Return success code
  var_put // [fnc] Write variable to output file
  (const int &nc_out, // I [fl] netCDF file for output 
   const std::string var_nm); // I [sng] Name of variable to output

private:
  
  // Static private class members
  static int nst_nbr; // [nbr] Number of instantiated class members
  static prc_cmp plrmu_avg_dfl; // [frc] Mean zenith angle, default
  static std::string rt_typ_dfl; // [sng] Radiative transfer method abbreviation, default
  static long wvl_nbr_dfl; // [nbr] Number of wavelength bands, default
  static const long wvl_nbr_max; // [nbr] Maximum number of wavelengths
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create option to abbreviation map
  static sng2sng_map opt2abb_map; // [map] Option to abbreviation map
  static const sng2rt_sct_map rt_map_mk(); // [fnc] Create radiative transfer method map
  static const sng2rt_sct_map rt_map; // [map] Radiative transfer method map
  static sng2var_mtd_sct_map var_mtd_map; // [map] Variable metadata map

  static sng2var_mtd_sct_map 
  var_mtd_map_mk // [fnc] Create variable metadata map
  (const int nc_out); // I [fl] netCDF file for output 
  
  // Private members
  bool flg_rta_in_dyn_mmr; // [flg] Reflectance/transmittance/absorptance are in dynamic memory
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  rt_sln_fnc_ptr_typ rt_sln_fnc_ptr; // [fnc] Pointer to function returning layer optical properties
  prc_cmp plrmu_avg; // [frc] Mean zenith angle
  int rcd; // [enm] Return success code
  long wvl_nbr; // [nbr] Number of wavelength bands
  std::string abb; // [sng] Abbreviation
  std::string dsc; // [sng] Description

  // Physical properties derived from inputs
  prc_cmp *rfl_flx; // [frc] Flux reflectance
  prc_cmp *trn_flx; // [frc] Flux transmittance
  prc_cmp *abs_flx; // [frc] Flux absorptance
  
  // Local dynamic copies of inputs
  prc_cmp *asm_prm; // I [frc] Asymmetry parameter
  prc_cmp *ss_alb; // I [frc] Single scattering albedo
  prc_cmp *tau_ext; // I [frc] Extinction optical depth
  
  // Private member functions
  int allocate(); // [fnc] Allocate dynamic memory for object
  int deallocate(); // [fnc] Free dynamic memory for object
  int recompute(); // [fnc] Recompute properties of object
  
  void prn()const; // [fnc] Print object contents
  
  var_mtd_sct // [sct] O Variable metadata structure
  var2mtd // [fnc] Variable name to metadata mapper
  (const std::string &var_nm); // I [sng] Variable name

}; // end class rt_cls

// Prototype global functions with C++ linkages

int // O [enm] Return success code
two_srm_iso_sct // [fnc] Two stream approximation for isotropically scattering slab
(const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const prc_cmp &plrmu_avg, // I [frc] Mean zenith angle
 const prc_cmp *asm_prm, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext, // I [frc] Extinction optical depth
 prc_cmp *rfl_flx, // O [frc] Flux reflectance
 prc_cmp *trn_flx, // O [frc] Flux transmittance
 prc_cmp *abs_flx); // O [frc] Flux absorptance
// end two_srm_iso_sct() prototype

int // O [enm] Return success code
two_srm_asm_sct // [fnc] Two stream approximation for anisotropically scattering slab
(const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const prc_cmp &plrmu_avg, // I [frc] Mean zenith angle
 const prc_cmp *asm_prm, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext, // I [frc] Extinction optical depth
 prc_cmp *rfl_flx, // O [frc] Flux reflectance
 prc_cmp *trn_flx, // O [frc] Flux transmittance
 prc_cmp *abs_flx); // O [frc] Flux absorptance
// end two_srm_asm_sct() prototype

int // O [enm] Return success code
rt_lop // [fnc] Radiative transfer layer optical properties
(const int &nc_out, // I [fl] netCDF file for output 
 const rt_cls &rt_obj, // I [rt] Radiative transfer object
 const prc_cmp &slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle
// end rt_lop() prototype

int // O [enm] Return success code
rt_rfl_snw // [fnc] Radiative transfer for snow reflectance
(const int &nc_out, // I [fl] netCDF file for output 
 const prc_cmp *flx_frc_drc_sfc, // I [frc] Surface insolation fraction in direct beam
 const prc_cmp &rfl_gnd_dff, // I [frc] Diffuse reflectance of ground (beneath snow)
 const rt_cls &rt_obj, // I [rt] Radiative transfer object
 const prc_cmp &slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle
// end rt_rfl_snw() prototype

// Define inline'd functions in header so source is visible to calling files

#endif // RT_HH






