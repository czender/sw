// $Id$ 

// Purpose: Aerosol physics

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <aer.hh> // Aerosol physics

#ifndef AER_HH // Contents have not yet been inserted in current source file  
#define AER_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cstdio> // stderr, FILE, NULL, etc.
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// 3rd party vendors
#include <gsl/gsl_sf_ellint.h> // GNU Scientific Library special functions Elliptic integrals

// Personal headers
#include <dbg.hh> // Debugging constants
#include <fio.hh> // File Input/Output Layer
#include <mth.hh> // Mathematical utilities, constants eqn_qdr(), eqn_cbc(), sng2cpx()
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

#include <phys_cst.hh> // Physical constants NB: Live code, not a header

#include <idx_rfr.hh> // Refractive indices

// Namespaces

// Forward declarations

/* Forward declaration of idx_rfr_cls breaks recursion caused by
   including idx_rfr_cls in aer_cls definition while including 
   aer.hh in idx_rfr.hh */
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

/*
typedef
int // O [enm] Return success code
(*idx_rfr_fnc_ptr_typ) // [fnc] Pointer to function returning refractive indices
  (const std::string cmp_sng_prt, // [sng] Composition of particle
   const prc_cmp idx_rfr_rl, // [frc] Real refractive index
   const prc_cmp idx_rfr_img, // [frc] Imaginary refractive index
   const prc_cmp wvl_ctr, // [m] Wavelength at band center
   const long &wvl_nbr); // [nbr] Number of output wavelength bands
// end idx_rfr_fnc_ptr_typ() prototype
*/

// Declare refractive index functions here
// Alternative is to #include <idx_rfr.hh> directly in aer.hh
extern idx_rfr_fnc_typ idx_rfr_YZS00_get; // [fnc] Pointer to function returning refractive indices

// Mineralogy structure
typedef struct{ // [sct] mnr_sct
  std::string abb; // [sng] Particle abbreviation
  std::string dsc; // [sng] Particle description
  //  std::string dsc_lng; // [sng] Long description
  std::string mlc; // [sng] Molecule
  std::string fl_idx_rfr; // [sng] File containing refractive indices  
  idx_rfr_fnc_ptr_typ idx_rfr_fnc_ptr; // [fnc] Pointer to function returning refractive indices
  prc_cmp mmw; // [kg mol-1] Mean molecular weight
  prc_cmp spc_heat; // [J kg-1 K-1] Specific heat capacity
  prc_cmp dns; // [kg m-3] Bulk density
} mnr_sct; // [sct] Mineralogy structure

typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map
typedef std::map<std::string,mnr_sct,std::less<std::string> > sng2mnr_sct_map; // Std::String to mnr_sct map

class aer_cls{

public:

  // Friends
  friend std::ostream & // [srm] Reference to output stream for cascading
  operator<< // [fnc] Stream insertion operator
  (std::ostream &srm_out, // [srm] Output stream
   const aer_cls &aer_obj); // [obj] Object to insert in stream

  // Static class members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  static std::string opt2abb(const std::string &opt_sng); // [sng] Convert option to abbreviation 
  static std::string abb2dsc(const std::string &abb_sng); // [sng] Convert abbreviation to description
  static std::string abb2fl_idx_rfr(const std::string &abb_sng); // [sng] Convert abbreviation to file
  static prc_cmp abb2dns(const std::string &abb_sng); // [kg m-3] Abbreviation to density mapper
  static prc_cmp abb2mmw(const std::string &abb_sng); // [kg mol-1] Abbreviation to molecular weight
  static prc_cmp abb2spc_heat(const std::string &abb_sng); // [J kg-1 K-1] Abbreviation to specific heat

  static idx_rfr_fnc_ptr_typ // [fnc] Pointer to function returning refractive indices
  abb2fnc // [fnc] Abbreviation to function mapper
  (const std::string &abb_sng); // [sng] Aerosol abbreviation

  int // O [enm] Return success code
  tst // [fnc] Self-test aer_cls class
  (const std::string &abb_sng, // [fnc] Aerosol abbreviation
   const long obj_nbr); // [nbr] Number of objects to create/destroy

  // Public member functions
  aer_cls // [fnc] Constructor
  (const std::string &abb_sng_arg=abb_sng_dfl, // [sng] Aerosol abbreviation
   const prc_cmp &dns_prt_arg=0.0, // [kg m-3] Density of particle
   const prc_cmp &tpt_arg=tpt_dfl, // [K] Temperature
   const std::string &fl_idx_rfr_arg="", // [sng] Aerosol input file
   const bool idx_rfr_usr_flg_arg=false, // [flg] Refractive index is user-specified
   const std::complex<prc_cmp> idx_rfr_usr_arg=static_cast<std::complex<prc_cmp> >(0.0)); // [frc] Refractive index, user-specified

  ~aer_cls(); // Destructor

  int // [enm] Return success code
  abb_set(const std::string &abb_sng_arg); // [sng] Aerosol abbreviation

  int // [enm] Return success code
  dns_set(const prc_cmp &dns_prt_arg); // [kg m-3] Density of particle

  int // [enm] Return success code
  fl_idx_rfr_set(const std::string &fl_idx_rfr_arg); // [sng] File containing refractive indices

  prc_cmp dns_get()const; // [kg m-3] Density of particle
  prc_cmp mmw_get()const; // [kg mol-1] Mean molecular weight
  prc_cmp spc_heat_get()const; // [J kg-1 K-1] Specific heat capacity
  std::string abb_get()const; // [sng] Aerosol abbreviation
  std::string dsc_get()const; // [sng] Aerosol description
  std::string fl_idx_rfr_get()const; // [sng] File containing refractive indices
  std::complex<prc_cmp>idx_rfr_usr_get()const; // [frc] Refractive index, user-specified

  std::complex<prc_cmp> // O [frc] Refractive index
  idx_rfr_get // [fnc] Return refractive index at specified wavelength
  (const prc_cmp wvl_ctr)const; // I [m] Wavelength at band center
  // end idx_rfr_get() prototype
  
  int // O [enm] Return success code
  idx_rfr_get // [fnc] Return array of refractive indices
  (const prc_cmp *wvl_ctr, // I [m] Wavelength at band center
   std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index
   const long &wvl_nbr=1L, // I [nbr] Number of wavelength bands
   const bool wrn_ntp_flg=true)const; // I [flg] Print WARNINGs from ntp_vec()
  // end idx_rfr_get() prototype

  void prn()const; // [fnc] Print object contents
  
private:

  // Static private class members
  static const prc_cmp tpt_dfl; // [K] Temperature, default
  static const sng2sng_map opt2abb_map; // [fnc] Option to abbreviation map
  static const sng2mnr_sct_map mnr_map; // [fnc] Mineral map
  static const std::complex<prc_cmp> idx_rfr_usr_dfl; // [frc] Default refractive index, user-specified
  static const std::string abb_sng_dfl; // [sng] Aerosol abbreviation, default
  static int nst_nbr; // [nbr] Number of instantiated class members
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create option to abbreviation map
  static sng2mnr_sct_map mnr_map_mk(); // [fnc] Create mineral map

  // Private members
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  idx_rfr_cls *idx_rfr; // [ptr] Refractive index object
  int rcd; // [enm] Return success code
  prc_cmp dns; // [kg m-3] Density
  prc_cmp mmw; // [kg mol-1] Mean molecular weight
  prc_cmp spc_heat; // [J kg-1 K-1] Specific heat capacity
  prc_cmp tpt; // [K] Temperature
  std::string abb; // [sng] Abbreviation
  std::string dsc; // [sng] Description

  // Private member functions
  int recompute();

}; // end class aer_cls

// Prototype functions with C++ linkages
int // O [rcd] Return success code
cff_drg_get
(const long in_nbr, // I [nbr] Size of arrays
 const prc_cmp *ryn_nbr, // I [frc] Reynolds number
 prc_cmp *cff_drg); // O [frc] Drag coefficient

int // O [rcd] Return success code
dmt_aer_get // [fnc] Compute aerodynamic diameter of particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_stk, // I [m] Stokes diameter
 const prc_cmp *slp_crc, // I [frc] Slip correction factor
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 prc_cmp *dmt_aer); // O [m] Aerodynamic diameter

int // O [rcd] Return success code
dmt_stk_get // [fnc] Compute Stokes diameter of particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 const prc_cmp *vlc_grv, // I [m s-1] Settling velocity
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *dmt_stk); // O [m] Stokes diameter

int // O [rcd] Return success code
vlc_grv_get // [fnc] Terminal fall speed of spherical particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp *slp_crc, // I [frc] Slip correction factor SeP97 p. 464
 const prc_cmp *vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
 prc_cmp *ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
 prc_cmp *stk_crc, // O [frc] Correction to Stokes settling velocity
 prc_cmp *vlc_grv); // O [m s-1] Settling velocity

int // O [rcd] Return success code
vlc_grv_get_Gin03 // [fnc] Terminal fall speed of aspherical particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *asp_rat_lps, // I [frc] Ellipsoidal aspect ratio
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp *slp_crc, // I [frc] Slip correction factor SeP97 p. 464
 const prc_cmp *vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
 prc_cmp *ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
 prc_cmp *stk_crc, // O [frc] Correction to Stokes settling velocity
 prc_cmp *vlc_grv); // O [m s-1] Settling velocity

int // O [rcd] Return success code
asp_rat_lps_get // [fnc] Get ellipsoidal aspect ratio
(const long sz_nbr, // I [nbr] Size of arrays
 const std::string cmp_sng_prt, // I [sng] Composition of particle
 const prc_cmp asp_rat_lps_dfl, // I [frc] Ellipsoidal aspect ratio
 prc_cmp *asp_rat_lps); // O [frc] Ellipsoidal aspect ratio

int // O [rcd] Return success code
vlc_grv_Gin03_get // [fnc] Terminal fall speed of aspherical particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *asp_rat_lps, // I [frc] Ellipsoidal aspect ratio
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp *slp_crc, // I [frc] Slip correction factor SeP97 p. 464
 const prc_cmp *vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
 prc_cmp *ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
 prc_cmp *stk_crc, // O [frc] Correction to Stokes settling velocity
 prc_cmp *vlc_grv); // O [m s-1] Settling velocity

int // O [rcd] Return success code
vlc_stk_get
(const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 const prc_cmp vsc_dyn_atm, // I [kg m-1 s-1] Dynamic viscosity of atmosphere
 prc_cmp *slp_crc, // O [frc] Slip correction factor
 prc_cmp *vlc_stk); // O [m s-1] Stokes' settling velocity (Re < 0.1)

int // O [rcd] Return success code
vnt_mss_get // [fnc] Mass ventilation coefficient
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *shm_nbr, // I [frc] Schmidt number
 const prc_cmp *ryn_nbr, // I [frc] Reynolds number
 prc_cmp *vnt_mss); // O [fnc] Mass ventilation coefficient

// Define inline'd functions in header so source is visible to calling files
inline prc_cmp // O [frc] Drag coefficient
cff_drg_fst_scl // [fnc] Drag coefficient
(const prc_cmp ryn_nbr) // I [frc] Reynolds number
{
  /* Purpose: "Fast scalar" drag coefficient for given Reynolds number
     Taken from Seinfeld and Pandis (1997) as reported in SeP97 p. 463 (8.32)
     Range of validity is < Re < 2.0e5
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics
     
     SeP97 expressions are not smoothly matched at boundaries
     This causes finite C_D jumps between continuously varying particle sizes
     For dust, jump occurs near D=80 um because Re changes from < 2 to > 2
     This jump nearly doubles sedimentation speed and is very unrealistic
     One approach to solve this is to blend solutions over limited range
     
     Pruppacher and Klett offer alternative formulations for intermediate C_D
     PrK78 p. 294 (10-50)--(10.52) and PrK98 p. 373 (10-51)--(10-53) have ln(Re/2)
     but SeP97 p. 463 (8.32) has ln(Re*2).
     Thus it appears one or the other has a typo

     GinO3 implies that Re < 2 for D < 100 um but really Re < 2	for D < 80 um
     so his figures are misleading. */
  using mth::cst_M_LN2l; // (0.6931471805599453094172321214581766L) [frc] log_e(2.0)
  using mth::cst_M_EULERl; // (0.577215664901532860606512090082L) [frc] Euler's constant
  const prc_cmp ryn_nbr_ntr_min(2.0); // [frc] Minimum Reynolds number for intermediate regime
  const prc_cmp ryn_nbr_ntr_max(5.0); // [frc] Maximum Reynolds number for intermediate regime
  prc_cmp wgt_mpr; // [frc] Weight of empirical parameterization regime
  prc_cmp cff_drg; // [frc] Drag coefficient
  if(ryn_nbr < 0.1){
    // Stokes regime
    cff_drg=24.0/ryn_nbr; // Stokes' law SeP97 p. 463 (8.32)
  }else if(ryn_nbr < ryn_nbr_ntr_min){ // SeP97 changes domains at Re=5.0
    // Asymptotic expansion regime
    // PrP57: Proudman and Pearson, J. Fluid Mech., 2, p. 237--262 (1957) 
    // cff_drg=(24.0/ryn_nbr)*(1.0+3.0*ryn_nbr/16.0+9.0*ryn_nbr*ryn_nbr*std::log(0.5*ryn_nbr)/160.0); // PrK78 p. 294 (10-51) PrK98 p. 373 (10-52)
    // Seinfeld and Pandis:
    // cff_drg=(24.0/ryn_nbr)*(1.0+3.0*ryn_nbr/16.0+9.0*ryn_nbr*ryn_nbr*std::log(2.0*ryn_nbr)/160.0); // SeP97 p. 463 (8.32)
    // ChB69: Chester and Breach, J. Fluid Mech., 37, p. 751--760 (1969) extension of PrP57:
    cff_drg=(24.0/ryn_nbr)*(1.0+3.0*ryn_nbr/16.0+9.0*ryn_nbr*ryn_nbr*(std::log(0.5*ryn_nbr)+mth::cst_M_EULERl+5.0*mth::cst_M_LN2l/3.0-323.0/360.0)/160.0+27.0*ryn_nbr*ryn_nbr*ryn_nbr*std::log(0.5*ryn_nbr)/640.0); // PrK78 p. 294 (10-52) PrK98 p. 373 (10-53)
  }else if(ryn_nbr < ryn_nbr_ntr_max){ // SeP97 changes domains at Re=5.0
    // Intermediate regime: Blend asymptotic expansion with pure empirical parameterization
    wgt_mpr=(ryn_nbr-ryn_nbr_ntr_min)/(ryn_nbr_ntr_max-ryn_nbr_ntr_min); // [frc] Weight of empirical parameterization regime
    cff_drg=(1.0-wgt_mpr)*(24.0/ryn_nbr)*(1.0+3.0*ryn_nbr/16.0+9.0*ryn_nbr*ryn_nbr*(std::log(0.5*ryn_nbr)+mth::cst_M_EULERl+5.0*mth::cst_M_LN2l/3.0-323.0/360.0)/160.0+27.0*ryn_nbr*ryn_nbr*ryn_nbr*std::log(0.5*ryn_nbr)/640.0)+wgt_mpr*(24.0/ryn_nbr)*(1.0+0.15*std::pow(ryn_nbr,PRC_CMP(0.687)));
  }else if(ryn_nbr < 500.0){ // Boundary of 500.0 seems arbitrary
    // Turbulent regime, pure empirical parameterization
    cff_drg=(24.0/ryn_nbr)*(1.0+0.15*std::pow(ryn_nbr,PRC_CMP(0.687))); // SeP97 p. 463 (8.32)
  }else if(ryn_nbr < 1.0e5){
    // fxm: Change to logarithmic relationship based on http://aerodyn.org/Drag/speed-drag.html
    cff_drg=0.44; // SeP97 p. 463 (8.32)
  }else if(ryn_nbr < 1.0e6){
    /* Turbulence transition regime
       Drag coeffients for Re > 20000 are elusive
       http://aerodyn.org/Drag/speed-drag.html shows data for different shapes
       Unclear what assumption are in these data
       Fgr. 4 appears to show that cff_dgr=0.44 decreases to cff_dgr ~ 0.1 
       as Re increases from ~ 1.0e5 to 1.0e6
       Decrease in drag with increase in Re at large Re is due to turbulence */
    cff_drg=0.44-0.34*(std::log10(ryn_nbr)-PRC_CMP(5.0)); // SeP97 p. 463 (8.32)
  }else{
    cff_drg=1.0e36; // Out of range error
  } // end else
  return cff_drg;
} // end cff_drg_fst_scl()

inline prc_cmp // O [frc] Drag coefficient for spheroids
cff_drg_Boo71_fst_scl // [fnc] Drag coefficient for spheroids
(const prc_cmp ryn_nbr, // I [frc] Reynolds number
 const prc_cmp sph_fct) // I [frc] Sphericity factor
{
  /* Purpose: "Fast scalar" drag coefficient for spheroids at given Reynolds number
     Analytic form for Re <= 2 taken from Boothroyd (1971) as described in Gin03
     Analytic form for Re >= 2 taken from Seinfeld and Pandis (1997) as reported in SeP97 p. 463 (8.32)
     fxm: must ensure continuity at Re == 2.0
     Range of physical validity is < Re < 2.0e5, but no aspherical data for Re > 2.0
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  prc_cmp cff_drg; // [frc] Drag coefficient

  // Determine aspherical drag as perturbation to spherical drag to avoid division by zero
  cff_drg=cff_drg_fst_scl(ryn_nbr);
  if(sph_fct != 1.0) cff_drg+=(24.0/ryn_nbr)*10.0*(1.0-sph_fct)*std::pow(ryn_nbr,PRC_CMP(0.35))/sph_fct; // Gin03 p. 2 (4)
  return cff_drg;
} // end cff_drg_Boo71_fst_scl()

inline prc_cmp // O [frc] Drag coefficient for hexagonal plates
cff_drg_plt_hxg_Wan02_fst_scl // [fnc] Drag coefficient for hexagonal plates
(const prc_cmp ryn_nbr) // I [frc] Reynolds number
{
  /* Purpose: "Fast scalar" drag coefficient for hexagonal plates at given Reynolds number
     Analytic form for Re < 0.22 and for Re > 150 from Seinfeld and Pandis (1997) as reported in SeP97 p. 463 (8.32)
     Analytic form for 0.2 < Re < 150 from Wang (2002)
     fxm: ensure continuity at regime interfaces
     Range of physical validity is < Re < 2.0e5, but no aspherical data for Re > 2.0
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  prc_cmp cff_drg; // [frc] Drag coefficient

  // Determine 
  if(ryn_nbr > 0.2 && ryn_nbr < 150.0){
    cff_drg=64.0*(1.0+0.078*pow(ryn_nbr,PRC_CMP(0.945)))/(mth::cst_M_PIl*ryn_nbr); // Wan02 p. 101 (3.19)
  }else{
    cff_drg=cff_drg_fst_scl(ryn_nbr);
  } // endif
  return cff_drg;
} // end cff_drg_plt_hxg_Wan02_fst_scl()

inline prc_cmp // O [frc] Equivalent diameter divided by semi-minor axis
psi_lps_fst_scl // [fnc] Equivalent diameter divided by semi-minor axis
(const prc_cmp asp_rat_lps) // I [frc] Ellipsoidal aspect ratio
{
  /* Purpose: "Fast scalar" equivalent diameter divided by semi-minor axis for given ellipsoidal aspect ratio
     Equivalent diameter is diameter of sphere with same surface area as given ellipsoid
     Psi is non-dimensional function defined such that Psi*b=dmt_eqv_sfc, so that
     Psi is surface area equivalent diameter divided by semi-minor axis
     Routine computes psi_lps from aspect ratio asp_rat_lps = a/b
     \lim_{\aspratlps \rightarrow 1} \Psi(\aspratlps) = 2
     Ginoux (2003) Gin03 p. 2 (10) has a typo/error in Psi definition
     Denominator under radical is "asp_rat_lps^2-1" not "1-asp_rat_lps^2"
     Former works fine for asp_rat >= 1 but latter is imaginary!
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  const prc_cmp asp_rat_lps_sqr(asp_rat_lps*asp_rat_lps); // [frc] Ellipsoidal aspect ratio squared
  return std::sqrt(PRC_CMP(2.0)+PRC_CMP(2.0)*asp_rat_lps_sqr/std::sqrt(asp_rat_lps_sqr-PRC_CMP(1.0))*std::asin(std::sqrt(PRC_CMP(1.0)-PRC_CMP(1.0)/asp_rat_lps_sqr))); // [frc] Equivalent diameter divided by semi-minor axis Gin03 p. 2 (10)
} // end psi_lps_fst_scl()

inline prc_cmp // O [m2] Surface area of prolate ellipsoid
lps_sfc_fst_scl // [fnc] Surface area of prolate ellipsoid
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" surface area of prolate ellipsoid
     Routine computes surface area from semi-major and minor axes a and b
     Ellipsoid is assumed to be prolate (rotated about semi-major axis) (a >= b)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  // Handle spheres first
  if(rds_lps_mjr == rds_lps_mnr) return PRC_CMP(4.0)*mth::cst_M_PIl*rds_lps_mjr*rds_lps_mjr;
  // Ellipsoidal area formula crashes when eccentricity is unity or degenerate
  assert(rds_lps_mjr*rds_lps_mnr != 0.0);
  const prc_cmp rds_lps_mjr_sqr(rds_lps_mjr*rds_lps_mjr); // [m2] 
  const prc_cmp rds_lps_mnr_sqr(rds_lps_mnr*rds_lps_mnr); // [m2] 
  const prc_cmp xcn_lps(std::sqrt(rds_lps_mjr_sqr-rds_lps_mnr_sqr)/rds_lps_mjr); // [frc] Eccentricity of ellipsoid
  return PRC_CMP(2.0)*mth::cst_M_PIl*(rds_lps_mjr_sqr+rds_lps_mnr_sqr*std::log((PRC_CMP(1.0)+xcn_lps)/(PRC_CMP(1.0)-xcn_lps))/(PRC_CMP(2.0)*xcn_lps)); // [m2] Surface area of prolate ellipsoid
} // end lps_sfc_fst_scl()

inline prc_cmp // O [m3] Volume of prolate ellipsoid
lps_vlm_fst_scl // [fnc] Volume of prolate ellipsoid
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" volume of prolate ellipsoid
     Routine computes volume from semi-major and minor axes a and b
     Ellipsoid is assumed to be prolate (rotated about semi-major axis) (a >= b)
     Routine easy to modify to work for any ellipse since volume formula is symmetric in a, b, c
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return PRC_CMP(4.0)*mth::cst_M_PIl*rds_lps_mjr*rds_lps_mnr*rds_lps_mnr/PRC_CMP(3.0); // [m3] Volume of prolate ellipsoid
} // end lps_vlm_fst_scl()

inline prc_cmp // O [m] Diameter of sphere with same surface area as ellipsoid
dmt_eqv_sfc_lps_fst_scl // [fnc] Diameter of sphere with same surface area as ellipsoid
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" equivalent diameter for given ellipsoid axes
     Equivalent diameter is diameter of sphere with same surface area as given ellipsoid
     Routine computes dmt_eqv_sfc_lps from semi-major and minor axes a and b
     Ellipsoid is assumed to be prolate (rotated about semi-major axis) (a >= b)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(lps_sfc_fst_scl(rds_lps_mjr,rds_lps_mnr)/mth::cst_M_PIl); // [m] Diameter of sphere with same surface area as ellipsoid
} // end dmt_eqv_sfc_lps_fst_scl()

inline prc_cmp // O [m] Diameter of sphere with same surface area as ellipsoid
dmt_eqv_sfc_lps_Gin03_fst_scl // [fnc] Diameter of sphere with same surface area as ellipsoid
(const prc_cmp asp_rat_lps, // I [frc] Ellipsoidal aspect ratio
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" equivalent diameter for given ellipsoidal aspect ratio
     Equivalent diameter is diameter of sphere with same surface area as given ellipsoid
     Routine computes dmt_eqv_sfc_lps from semi-minor axis b and aspect ratio asp_rat_lps = a/b
     Taken from Ginoux (2003) Gin03 p. 2 (10)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return rds_lps_mnr*psi_lps_fst_scl(asp_rat_lps); // [m] Diameter of sphere with same surface area as ellipsoid Gin03 p. 2 (10)
} // end dmt_eqv_sfc_lps_Gin03_fst_scl()

inline prc_cmp // O [m] Diameter of sphere with same volume as ellipsoid
dmt_eqv_vlm_lps_fst_scl // [fnc] Diameter of sphere with same volume as ellipsoid
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" equivalent diameter for given ellipsoid axes
     Equivalent diameter is diameter of sphere with same volume as given ellipsoid
     Routine computes dmt_eqv_vlm_lps from semi-major and minor axes a and b
     Ellipsoid is assumed to be prolate (rotated about semi-major axis) (a >= b)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::pow(PRC_CMP(6.0)*lps_vlm_fst_scl(rds_lps_mjr,rds_lps_mnr)/static_cast<prc_cmp>(mth::cst_M_PIl),PRC_CMP(1.0)/PRC_CMP(3.0)); // [frc] Diameter of sphere with same volume as ellipsoid
} // end dmt_eqv_vlm_lps_fst_scl()

inline prc_cmp // O [frc] Ellipticity factor
lpt_fct_fst_scl // [fnc] Ellipticity factor of given ellipsoidal aspect ratio
(const prc_cmp asp_rat_lps) // I [frc] Ellipsoidal aspect ratio
{
  /* Purpose: "Fast scalar" ellipticity factor for a given ellipsoidal aspect ratio
     Ellipticity factor lpt_fct of a prolate spheroid is 
     Routine computes lpt_fct from aspect ratio asp_rat_lps = a/b
     Taken from Ginoux (2003) Gin03 p. 2 (7), http://mathworld.wolfram.com/ProlateSpheroid.html
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(asp_rat_lps*asp_rat_lps-PRC_CMP(1.0))/asp_rat_lps; // [frc] Ellipticity factor Gin03 p. 2 (8)
} // end lpt_fct_fst_scl()

inline prc_cmp // O [frc] Eccentricity factor
xcn_lps_prl_fst_scl // [fnc] Eccentricity factor of given ellipsoidal aspect ratio
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" eccentricity factor for a given ellipsoidal aspect ratio
     Eccentricity factor xcn_lps_prl of prolate spheroid is sqrt(a^2-b^2)/a
     Routine computes xcn_lps_prl from semi-major and semi-minor axes a, b (a > b)
     Taken from Ginoux (2003) Gin03 p. 2 (7), http://mathworld.wolfram.com/ProlateSpheroid.html
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(rds_lps_mjr*rds_lps_mjr-rds_lps_mnr*rds_lps_mnr)/rds_lps_mjr; // [frc] Eccentricity factor of ellipsoid 
} // end xcn_lps_prl_fst_scl()

inline prc_cmp // O [frc] Eccentricity factor
xcn_lps_obl_fst_scl // [fnc] Eccentricity factor of given ellipsoidal aspect ratio
(const prc_cmp rds_lps_mjr, // I [m] Semi-major axis (a) of ellipsoid
 const prc_cmp rds_lps_mnr) // I [m] Semi-minor axis (b) of ellipsoid
{
  /* Purpose: "Fast scalar" eccentricity factor for a given ellipsoidal aspect ratio
     Eccentricity factor xcn_lps_obl of oblate spheroid is sqrt(a^2-b^2)/a
     Routine computes xcn_lps_obl from semi-major and semi-minor axes a, b (a > b)
     Taken from Ginoux (2003) Gin03 p. 2 (7), http://mathworld.wolfram.com/OblateSpheroid.html
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(rds_lps_mjr*rds_lps_mjr-rds_lps_mnr*rds_lps_mnr)/rds_lps_mjr; // [frc] Eccentricity factor of ellipsoid 
} // end xcn_lps_obl_fst_scl()

inline prc_cmp // O [frc] Eccentricity factor
xcn_lps_prl_fst_scl // [fnc] Eccentricity factor of given ellipsoidal aspect ratio
(const prc_cmp asp_rat_lps) // I [frc] Ellipsoidal aspect ratio
{
  /* Purpose: "Fast scalar" eccentricity factor for a given ellipsoidal aspect ratio
     Eccentricity factor xcn_lps_prl of prolate spheroid is sqrt(a^2-b^2)/a
     Routine computes xcn_lps_prl from aspect ratio asp_rat_lps = a/b
     Taken from Ginoux (2003) Gin03 p. 2 (7), http://mathworld.wolfram.com/ProlateSpheroid.html
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(asp_rat_lps*asp_rat_lps-PRC_CMP(1.0))/asp_rat_lps; // [frc] Eccentricity factor Gin03 p. 2 (7)
} // end xcn_lps_prl_fst_scl()

inline prc_cmp // O [frc] Eccentricity factor
xcn_lps_obl_fst_scl // [fnc] Eccentricity factor of given ellipsoidal aspect ratio
(const prc_cmp asp_rat_lps) // I [frc] Ellipsoidal aspect ratio
{
  /* Purpose: "Fast scalar" eccentricity factor for a given ellipsoidal aspect ratio
     Eccentricity factor xcn_lps_obl of oblate spheroid is sqrt(a^2-b^2)/a
     Routine computes xcn_lps_obl from aspect ratio asp_rat_lps = a/b
     Taken from Ginoux (2003) Gin03 p. 2 (7), http://mathworld.wolfram.com/OblateSpheroid.html
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  return std::sqrt(asp_rat_lps*asp_rat_lps-PRC_CMP(1.0))/asp_rat_lps; // [frc] Eccentricity factor Gin03 p. 2 (7)
} // end xcn_lps_obl_fst_scl()

inline prc_cmp // O [frc] Slip correction factor
slp_crc_fst_scl // [fnc] 
(const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 const prc_cmp dmt_prt) // I [m] Particle diameter
{
  /* Purpose: "Fast scalar" slip correction for a given particle diameter
     Taken from Seinfeld and Pandis (1997) SeP97 p. 464 (8.34)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  prc_cmp slp_crc; // [frc] Slip correction factor
  slp_crc=1.0+2.0*mfp_atm*(1.257+0.4*std::exp(-1.1*dmt_prt/(2.0*mfp_atm)))/dmt_prt; // [frc] Slip correction factor SeP97 p. 464 (8.34)
  return slp_crc; // [frc] Slip correction factor
} // end slp_crc_fst_scl()

inline prc_cmp // O [frc] Sphericity factor
sph_fct_fst_scl // [fnc] Sphericity factor of given ellipsoidal aspect ratio
(const prc_cmp asp_rat_lps) // I [frc] Ellipsoidal aspect ratio
{
  /* Purpose: "Fast scalar" sphericity factor for a given ellipsoidal aspect ratio
     Boo71 defines sphericity factor Phi as surface area of spherically equivalent volume
     divided by actual particle surface area
     Taken from Ginoux (2003) Gin03 p. 2 (8)
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  if(asp_rat_lps == PRC_CMP(1.0)) return PRC_CMP(1.0);
  const prc_cmp asp_rat_lps_sqr(asp_rat_lps*asp_rat_lps); // [frc] Ellipsoidal aspect ratio squared
  prc_cmp sph_fct; // [frc] Sphericity factor
  sph_fct=2.0*std::pow(asp_rat_lps,PRC_CMP(2.0)/PRC_CMP(3.0)); // [frc] Sphericity factor Gin03 p. 2 (8)
  sph_fct/=1.0+(asp_rat_lps_sqr/std::sqrt(asp_rat_lps_sqr-PRC_CMP(1.0))*std::asin(std::sqrt(PRC_CMP(1.0)-PRC_CMP(1.0)/asp_rat_lps_sqr))); // [frc] Sphericity factor Gin03 p. 2 (8)
  return sph_fct; // [frc] Sphericity factor
} // end sph_fct_fst_scl()

inline prc_cmp // O [frc] Mass or heat ventilation coefficient
vnt_cff_fst_scl // [fnc] Mass or heat ventilation coefficient
(const prc_cmp shm_nbr, // I [frc] Schmidt number
 const prc_cmp ryn_nbr) // I [frc] Reynolds number
{
  /* Purpose: "Fast scalar" mass ventilation coefficient for a given fluid numbers
     Taken from PrK98 p. 541 (13-60)
     Ventilation coefficients for mass and thermal diffusion are identical in form
     For mass ventilation coeffecient, feed this routine the Schmidt number for vapor diffusion = vsc_knm_air / dfs_cff_vpr
     For heat ventilation coeffecient, feed this routine the Schmidt number for thermal diffusion = vsc_knm_air / cnd_trm_air
     In order to maximize chances of compiler inlining this function,
     it avoids error checking and diagnostics */
  const prc_cmp vnt_fct(std::pow(shm_nbr,PRC_CMP(1.0)/PRC_CMP(3.0))*std::sqrt(ryn_nbr)); // [frc] Factor in ventilation coefficient PrK98 p. 541 (13-60)
  // Ensure domain is valid
  assert(vnt_fct < 51.4);
  return (vnt_fct < 1.4) ? 1.0+0.108*vnt_fct*vnt_fct : 0.78+0.308*vnt_fct; // [frc] Mass or heat ventilation coefficient PrK98 p. 541 (13-60)
} // end vnt_cff_fst_scl()

#endif // AER_HH  






