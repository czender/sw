// $Id$ 

// Purpose: Description (definition) of Mie scattering solutions and utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie_sln.hh> // Mie scattering solutions and utilities

#ifndef MIE_SLN_HH // Contents have not yet been inserted in current source file
#define MIE_SLN_HH

// Standard C++ headers 
#include <complex> // Standard C++ complex class
#include <valarray> // STL valarray class template

// Standard C headers 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc.

// 3rd party vendors
#include <gsl/gsl_sf_gamma.h> // GNU Scientific Library special functions gamma functions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants */
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <mie_MaS99.hh> // Mie scattering solutions from MaS99
#include <mie_Wis79.hh> // Mie scattering solutions from Wis79

// Namespaces
// fxm: TODO 75 this does not work
/*
namespace mie_ctl{ // [nms] Control parameters for Mie solutions
  long wvl_idx_dbg_mie; // [idx] Debugging wavelength bin
} // end namespace mie_ctl
*/

// Typedefs

// Forward declarations

// Prototype functions with C++ linkages

int // O [enm] Return success code
mie_prc // [fnc] Mie processor
(const bool abs_ncl_wk_mdm_flg, // I [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
 const bool coat_flg, // I [flg] Assume coated spheres
 const bool slf_tst_flg, // I [frg] Self-test flag
 const long bnd_idx, // I [idx] Counting index for band
 const long lgn_nbr, // I [nbr] Order of phase function Legendre expansion
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const long sz_idx, // I [idx] Counting index for size
 const long wvl_idx, // I [idx] Counting index for wavelength
 const long wvl_idx_dbg, // I [idx] Debugging wavelength bin
 const prc_cmp bnd_ctr, // I [m] Wavelength at band center
 const prc_cmp dmn_frc, // I [frc] Fractal dimensionality of inclusions
 const prc_cmp rds_cor, // I [m] Radius of core
 const prc_cmp rds_mnt, // I [m] Radius of mantle
 const prc_cmp sz_ctr_sph, // I [m] Equivalent sphere size at bin center
 const prc_cmp sz_prm_rsn_usr_spc, // I [m m-1] Size parameter resolution, user specified
 const std::complex<prc_cmp> idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> idx_rfr_ffc, // I [frc] Effective refractive index of particle
 const std::complex<prc_cmp> idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> idx_rfr_prt, // I [frc] Refractive index of particle
 const std::string slv_sng, // I [sng] Mie solver to use (BoH83 or Wis79)
 double &abs_fct_MaS99, // O [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms, // O [frc] Hemispheric backscatter
 double &q_abs, // O [frc] Absorption efficiency
 double &q_bck, // O [frc] Backscattering efficiency
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 double &spk_val, // O [frc] Mie coefficient of spike
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 prc_cmp * const phz_fnc, // O [sr-1] Phase function
 prc_cmp * const plz); // O [frc] Polarization
// !mie_prc()

int // O [enm] Return success code
mie_sph_BoH83 // [fnc] Mie solution for homogeneous spheres
(const double sz_prm, // I [m m-1] Size parameter
 const std::complex<double> &idx_rfr_rlt, // I [frc] Refractive index of particle (aerosol) relative to medium (air)
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms, // O [frc] Hemispheric backscatter
 double &q_bck, // O [frc] Backscattering efficiency
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 std::complex<double> * const s1, // O [frc] Scalar amplitude scattering matrix element S1
 std::complex<double> * const s2); // O [frc] Scalar amplitude scattering matrix element S2
// !mie_sph_BoH83()

int // O [enm] Return success code
mie_sph_coat_BoH83 // [fnc] Mie solution for coated spheres
(const double sz_prm_cor, // I [m m-1] Size parameter of core
 const double sz_prm_mnt, // I [m m-1] Size parameter of mantle
 const std::complex<double> idx_rfr_rlt_cor, // I [frc] Core index of refraction relative to medium
 const std::complex<double> idx_rfr_rlt_mnt, // I [frc] Mantle index of refraction relative to medium
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 double &q_bck, // O [frc] Backscattering efficiency
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms); // O [frc] Hemispheric backscatter
// !mie_sph_coat_BoH83()

int // O [rcd] Return success code
mie_bck_hms_Chy73 // [fnc] Mie solution for hemispheric backscatter
(const double sz_prm, // I [m m-1] Size parameter
 const long trm_nbr, // I [nbr] Number of terms
 const std::complex<double> * const an, // I [frc] an coefficients
 const std::complex<double> * const bn, // I [frc] bn coefficients
 double &bck_hms); // O [frc] Hemispheric backscatter
// !mie_bck_hms_Chy73()

int // [rcd] Return success code
mie_ngl_BoH83_csz // [fnc] Mie phase function for homogeneous spheres
(const double q_sct, // I [frc] Scattering efficiency
 const double sz_prm, // I [m m-1] Size parameter
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const std::complex<double> * const s1, // I [frc] Scalar amplitude scattering matrix element S1
 const std::complex<double> * const s2, // I [frc] Scalar amplitude scattering matrix element S2
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 prc_cmp *phz_fnc, // O [sr-1] Phase function
 prc_cmp *plz); // O [frc] Polarization
// !mie_ngl_BoH83_csz()

int // O [enm] Return success code
adt_apx // [fnc] Anomalous Diffraction Theory approximation
(const double &sz_prm, // I [m m-1] Size parameter
 const std::complex<prc_cmp> &idx_rfr_ffc, // I [frc] Effective refractive index of particle
 const std::complex<prc_cmp> &idx_rfr_mdm_img, // I [frc] Refractive index of medium
 double &q_abs, // I/O [frc] Absorption efficiency
 double &q_ext, // I/O [frc] Extinction efficiency
 double &q_sct); // I/O [frc] Scattering efficiency
// !adt_apx()

#endif // MIE_SLN_HH  
