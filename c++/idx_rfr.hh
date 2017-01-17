// $Id$ 

// Purpose: Description (definition) of refractive index classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <idx_rfr.hh> // Refractive indices

#ifndef IDX_RFR_HH // Contents have not yet been inserted in current source file
#define IDX_RFR_HH

// Standard C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <string> // Standard C++ string class
#include <complex> // Standard C++ complex class

// Standard C headers
#include <cstring> // strcmp... 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdlib> // atof, atoi, malloc, getopt, strtod
#include <unistd.h> // All sorts of POSIX stuff  

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants eqn_qdr(), eqn_cbc(), sng2cpx()
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <vec.hh> // Vector functions ntp_vec(), rbn_vec()
#include <xtr.hh> // Extrapolation class

// aer_cls::abb2dsc() is called in idx_rfr_cls::recompute()
#include <aer.hh> // Aerosol physics
#include <idx_rfr_H2O.hh> // Refractive indices of water
#include <phys_cst.hh> // Physical constants NB: Live code, not a header

#include <libnco_c++.hh> // C++ interface to netCDF C library

// Namespaces

// Forward declarations:

class aer_cls; // [cls] Aerosol class
// idx_rfr_cls is used in typedefs
class idx_rfr_cls; // [cls] Refractive index class

// Typedefs
typedef
int // O [enm] Return success code
idx_rfr_fnc_typ // [fnc] Function returning refractive indices
  (const std::string cmp_sng_prt, // [sng] Composition of particle
   const prc_cmp idx_rfr_rl, // [frc] Real refractive index
   const prc_cmp idx_rfr_img, // [frc] Imaginary refractive index
   const prc_cmp wvl_ctr, // [m] Wavelength at band center
   const long &wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_fnc_typ() prototype

typedef
idx_rfr_fnc_typ // [fnc] Pointer to function returning refractive indices
*idx_rfr_fnc_ptr_typ; // [fnc] Pointer to function returning refractive indices
// end idx_rfr_fnc_ptr_typ() prototype

typedef std::map<std::string,idx_rfr_cls *,std::less<std::string> > sng2idx_rfr_map; // Stands for "string to idx_rfr map"

// Enums
enum ffc_mdm_typ{ // [enm] Effective medium type
  ffc_mdm_brg, // [enm] Effective medium approximation: Bruggeman
  ffc_mdm_mmg, // [enm] Effective medium approximation: Multiple Maxwell Garnett
  ffc_mdm_mxg, // [enm] Effective medium approximation: Maxwell Garnett
  ffc_mdm_pmr, // [enm] Effective medium approximation: Partial molar refraction
  ffc_mdm_vlw, // [enm] Effective medium approximation: Volume-weighted
  ffc_mdm_nil}; // [enm] Effective medium approximation: None (undefined, not-used)
// end ffc_mdm_enm for Effective medium type 

enum mca_typ{ // [enm] Multi-component aerosol type
  mca_cor, // [enm] Multi-component aerosol type: Core
  mca_mdm, // [enm] Multi-component aerosol type: Medium
  mca_mnt, // [enm] Multi-component aerosol type: Mantle
  mca_mtx, // [enm] Multi-component aerosol type: Matrix
  mca_ncl, // [enm] Multi-component aerosol type: Inclusion
  mca_prt, // [enm] Multi-component aerosol type: Particle
  mca_nil}; // [enm] Multi-component aerosol type: None (undefined, not-used)
// end mca_typ for Multi-component aerosol type

// Define idx_rfr_cls class
class idx_rfr_cls{
  
public:
  
  // Friends
  friend std::ostream & // [srm] Reference to output stream for cascading
  operator<< // [fnc] Stream insertion operator
  (std::ostream &srm_out, // [srm] Output stream
   const idx_rfr_cls &idx_rfr_obj); // [obj] Object to insert in stream
  
  // Static class members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  static std::string abb2fl_idx_rfr(const std::string &abb_sng); // [sng] Convert abbreviation to file
  
  static idx_rfr_fnc_ptr_typ // [fnc] Pointer to function returning refractive indices
  abb2fnc // [fnc] Abbreviation to function mapper
  (const std::string &abb_sng); // [sng] Refractive index abbreviation
  
  static idx_rfr_cls * // [ptr] Pointer to refractive index object
  abb2obj // [fnc] Abbreviation to object mapper
  (const std::string &abb_sng); // [fnc] Refractive index abbreviation
  
  static int // O [enm] Return success code
  tst // [fnc] Self-test idx_rfr_cls class
  (const std::string &abb_sng, // [fnc] Refractive index abbreviation
   const long obj_nbr); // [nbr] Number of objects to create/destroy
  
  // Public member functions
  
  idx_rfr_cls // [fnc] Constructor
  (const std::string &abb_sng_arg=idx_rfr_typ_dfl, // [sng] Refractive index abbreviation
   const std::string &fl_idx_rfr_arg="", // [sng] File containing refractive indices
   const prc_cmp &tpt_arg=tpt_dfl, // [K] Temperature
   const bool idx_rfr_usr_flg_arg=false, // [flg] Refractive index is user-specified
   const std::complex<prc_cmp> idx_rfr_usr_arg=idx_rfr_usr_dfl); // [frc] Refractive index, user-specified
  // end idx_rfr_cls::idx_rfr_cls() prototype
  
  ~idx_rfr_cls(); // [fnc] Destructor
  
  int // [enm] Return success code
  abb_set(const std::string &abb_sng_arg); // [sng] Refractive index type
  
  int // [enm] Return success code
  tpt_set(const prc_cmp &tpt_arg); // [K] Temperature
  
  int // [enm] Return success code
  idx_rfr_usr_set(const std::complex<prc_cmp> &idx_rfr_usr_arg); // [frc] Refractive index, user-specified
  
  int // [enm] Return success code
  idx_rfr_fnc_ptr_set(const idx_rfr_fnc_ptr_typ &idx_rfr_fnc_ptr_arg); // [fnc] Function to compute refractive indices
  
  int // [enm] Return success code
  fl_idx_rfr_set(const std::string &fl_idx_rfr_arg=""); // [sng] File containing refractive indices
  
  std::string abb_get()const; // [sng] Refractive index abbreviation
  std::string dsc_get()const; // [sng] Refractive index description
  std::string fl_idx_rfr_get()const; // [sng] File containing refractive indices
  std::complex<prc_cmp>idx_rfr_usr_get()const; // [frc] Refractive index, user-specified
  prc_cmp tpt_get()const; // [K] Temperature
  idx_rfr_fnc_ptr_typ idx_rfr_fnc_ptr_get()const; // [fnc] Function to compute refractive indices

  std::complex<prc_cmp> // O [frc] Refractive index
  idx_rfr_get // [fnc] Return refractive index at specified wavelength
  (const prc_cmp wvl_ctr)const; // I [m] Wavelength at band center
  // end idx_rfr_cls::idx_rfr_get() prototype
  
  int // O [enm] Return success code
  idx_rfr_get // [fnc] Return array of refractive indices
  (const prc_cmp *wvl_ctr, // I [m] Wavelength at band center
   std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index
   const long &wvl_nbr=1L, // I [nbr] Number of wavelength bands
   const bool wrn_ntp_flg=true)const; // I [flg] Print WARNINGs from ntp_vec()
  // end idx_rfr_cls::idx_rfr_get() prototype
  
  void prn()const; // [fnc] Print object contents
  
private:
  
  // Static private class members
  static const prc_cmp tpt_dfl; // [K] Temperature, default
  static const std::complex<prc_cmp> idx_rfr_usr_dfl; // [frc] Default refractive index, user-specified
  static const std::string idx_rfr_typ_dfl; // [sng] Refractive index abbreviation, default
  static int nst_nbr; // [nbr] Number of instantiated class members
  static sng2idx_rfr_map abb2obj_map; // [map] Abbreviation to object map
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create option to abbreviation map
  
  // Private members

  /* Making rcd a class member is convenient but also problematic because 
     often rcd is only structure member a function modifies and if rcd were not 
     a member then entire function could be declared const.
     Hence the great compromise convention:
     Public get() functions always modify and return local rcd_lcl (not private member rcd)
     Private get() and Public and Private set() functions always modify and return private member rcd */

  bool flg_raw_idx_rfr_in_dyn_mmr; // [flg] Raw refractive indices from file are in dynamic memory
  bool idx_rfr_usr_flg; // [flg] Refractive index is user-specified
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  idx_rfr_fnc_ptr_typ idx_rfr_fnc_ptr; // [fnc] Pointer to function returning refractive indices
  int rcd; // [enm] Return success code
  long wvl_nbr_tbl; // [nbr] Number of wavelength bands, tabulated
  prc_cmp *wvl_ctr_tbl; // [m] Wavelength at band center, tabulated
  prc_cmp tpt; // [K] Temperature
  std::complex<prc_cmp> *idx_rfr_tbl; // [frc] Refractive index, tabulated
  std::complex<prc_cmp> idx_rfr_usr; // [frc] Refractive index, user-specified
  std::string abb; // [sng] Abbreviation
  std::string dsc; // [sng] Description
  std::string fl_idx_rfr; // [sng] File containing refractive indices
  
  // Private member functions
  int allocate(); // [fnc] Allocate dynamic memory for object
  int deallocate(); // [fnc] Free dynamic memory for object
  int recompute(); // [fnc] Recompute properties of object
  
  int // [enm] Return success code
  abb2obj_map_add // [fnc] Register new refractive index object
  (idx_rfr_cls *idx_rfr_ptr); // [obj] Refractive index object
  
  int // [enm] Return success code
  idx_rfr_prn // [fnc] Print refractive indices in pretty table
  (prc_cmp *wvl_ctr, // I [m] Wavelength at band center
   long wvl_nbr, // I [nbr] Number of wavelength bands
   std::complex<prc_cmp> *idx_rfr)const; // I [frc] Refractive index
  // end idx_rfr_prn() prototype

  int // [enm] Return success code
  idx_rfr_raw_fl_get // [fnc] Load raw refractive index data from file
  (const std::string &cmp_sng_prt, // I [sng] Composition of particle
   const std::string &fl_in, // I [sng] File containing refractive indices
   prc_cmp *&wvl_ctr, // O [m] Wavelength at band center, tabulated
   long &wvl_nbr, // O [nbr] Number of wavelength bands, tabulated
   std::complex<prc_cmp> *&idx_rfr); // O [frc] Refractive index, tabulated
  // end idx_rfr_raw_fl_get() prototype

  int // O [enm] Return success code
  ntp_idx_rfr // [fnc] Interpolate raw refractive indices to requested grid
  (const prc_cmp *wvl_ctr_out, // I [m] Wavelength at band center
   const long wvl_nbr_out, // I [nbr] Number of wavelength bands
   prc_cmp *idx_rfr_out); // O [frc] Refractive index
  // end ntp_idx_rfr() prototype
  
}; // end class idx_rfr_cls

// Prototype functions with C++ linkages

int // [enm] Return success code
idx_rfr_usr_set // [fnc] Set refractive index to user-specified value
(const std::complex<prc_cmp> &idx_rfr_usr, // [frc] User-specified refractive index
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index
 const long wvl_nbr); // I [nbr] Number of output wavelength bands
// end idx_rfr_usr_set() prototype

int // [rcd] Return success code
idx_rfr_usr_set // [fnc] Set refractive index to user-specified value
(const prc_cmp &idx_rfr_rl_usr, // I [frc] User-specified real refractive index
 const prc_cmp &idx_rfr_img_usr, // I [frc] User-specified imaginary refractive index
 prc_cmp *idx_rfr_rl, // O [frc] Real refractive index
 prc_cmp *idx_rfr_img, // O [frc] Imaginary refractive index
 const long wvl_nbr); // I [nbr] Number of output wavelength bands
// end idx_rfr_usr_set() prototype

int // O [enm] Return success code
idx_rfr_air_get // [fnc] Compute refractive indices for air
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of air
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr); // I [nbr] Number of output wavelength bands
// end idx_rfr_air_get() prototype

int // O [enm] Return success code
idx_rfr_ChC90_get // [fnc] Compute refractive indices for soot aerosol of ChC90
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr); // I [nbr] Number of output wavelength bands
// end idx_rfr_ChC90_get() prototype

int // [enm] Return success code
idx_rfr_DKS91_get // [fnc] Compute refractive indices for soot aerosol of DKS91
(const std::string cmp_sng_prt, // [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // [m] Wavelength
 const long wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_DKS91_get() prototype

int // [rcd] Return success code
idx_rfr_LQB93_get // [fnc] Compute refractive indice for calcite or gypsum aerosol of LQB93
(const std::string cmp_sng_prt, // [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of calcite or gypsum
 const prc_cmp *wvl_ctr, // [m] Wavelength at band center
 const long wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_LQB93_get() prototype

int // [enm] Return success code
idx_rfr_WaW80_get // [fnc] Compute refractive indices for soot aerosol of WaW80
(const std::string cmp_sng_prt, // [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // [m] Wavelength
 const long wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_WaW80_get() prototype

int // [rcd] Return success code
idx_rfr_YZS00_get // [fnc] Compute refractive indices for soot aerosol of YZS00
(const std::string cmp_sng_prt, // [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl_ctr, // [m] Wavelength at band center
 const long wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_YZS00_get() prototype

int // [enm] Return success code
idx_rfr_ffc_get // [fnc] Compute effective optical constants of composites
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const int ffc_mdm_typ, // I [enm] Effective medium type
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction(s) in inclusion(s)
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> *idx_rfr_prt, // I [frc] Refractive index of particle
 std::complex<prc_cmp> *idx_rfr_ffc_brg, // O [frc] Effective refractive index, Bruggeman approximation
 std::complex<prc_cmp> *idx_rfr_ffc_mxg, // O [frc] Effective refractive index, Maxwell Garnett approximation
 std::complex<prc_cmp> *idx_rfr_ffc_pmr, // O [frc] Effective refractive index, partial molar refraction approximation
 std::complex<prc_cmp> *idx_rfr_ffc_vlw, // O [frc] Effective refractive index, volume-weighted approximation
 std::complex<prc_cmp> *idx_rfr_ffc); // O [frc] Effective refractive index of particle
// end idx_rfr_ffc_get() prototype

int // [enm] Return success code
idx_rfr_ffc_get // [fnc] Return specified Effective Medium Approximation
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const int ffc_mdm_typ, // I [enm] Effective medium type
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction(s) in inclusion(s)
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> *idx_rfr_prt, // I [frc] Refractive index of particle
 std::complex<prc_cmp> *idx_rfr_ffc); // O [frc] Effective refractive index of particle
// end idx_rfr_ffc_get() prototype

int // [enm] Return success code
idx_rfr_ffc_brg_get // [fnc] Effective Medium Approximation: Bruggeman
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_brg); // O [frc] Effective refractive index, Bruggeman approximation
// end idx_rfr_ffc_brg_get() prototype

int // [enm] Return success code
idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_mxg); // O [frc] Effective refractive index, Maxwell Garnett approximation
// end idx_rfr_ffc_mxg_get() prototype

int // [enm] Return success code
idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett, multiple inclusions
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> **idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_mxg); // O [frc] Effective refractive index, Maxwell Garnett approximation
// end idx_rfr_ffc_mxg_get() prototype

int // [enm] Return success code
idx_rfr_ffc_pmr_get // [fnc] Compute refractive indices using partial molar refraction approximation
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp dns_cor, // [kg m-3] Density of core
 const prc_cmp dns_mnt, // [kg m-3] Density of mantle
 const prc_cmp dns_ttl, // [kg m-3] Density of mixture
 const prc_cmp mlc_wgt_cor, // [kg mol-1] Molecular weight of core
 const prc_cmp mlc_wgt_mnt, // [kg mol-1] Molecular weight of mantle
 const prc_cmp mlc_wgt_ttl, // [kg mol-1] Molecular weight of mixture
 const prc_cmp vlm_frc_cor, // [frc] Volume fraction in core
 const prc_cmp vlm_frc_mnt, // I [frc] Volume fraction in mantle
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 std::complex<prc_cmp> *idx_rfr_ffc_pmr); // O [frc] Effective refractive index, partial molar refraction approximation
// end idx_rfr_ffc_pmr_get() prototype

int // [enm] Return success code
idx_rfr_ffc_pmr_lst_get // [fnc] Compute refractive indices using partial molar refraction approximation
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const aer_cls *aer_lst, // I [aer] List of aerosol components
 const long &aer_nbr, // I [frc] Number of aerosol components
 std::complex<prc_cmp> *idx_rfr_ffc_pmr); // O [frc] Effective refractive index, partial molar refraction approximation
// end idx_rfr_ffc_pmr_lst_get() prototype

// Define inline'd functions in header so source is visible to calling files

#endif // IDX_RFR_HH  
