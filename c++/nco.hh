// $Id$ 

// Purpose: Description (definition) of netCDF utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <nco.hh> // netCDF utilities

#ifndef NCO_HH // Contents have not yet been inserted in current source file
#define NCO_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers
#include <netcdf.h> // netCDF C interface

// 3rd party vendors
#include <netcdf.h> // netCDF C interface
#include <gsl/gsl_complex.h> // GNU Scientific Library complex types
#include <gsl/gsl_complex_math.h> // GNU Scientific Library complex math functions

// Personal headers
#include <a2d.hh> // Two dimensional arrays
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <aer.hh> // Aerosol physics
#include <spc_bbd.hh> // Blackbody spectra
#include <spc_slr.hh> // Solar spectra
#include <psd.hh> // Particle size distributions
#include <tdy.hh> // Atmospheric thermodynamics
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Typedefs

// Define nco_cls class

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const gsl_complex &gsl_cpx); // [obj] Object to insert in stream

// Prototype global functions with C++ linkages

int // O [enm] Return success code
aer_htg // [fnc] Determine aerosol heating characteristics
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp &prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp &wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 const prc_cmp &dmt_dtc, // I [m] Diameter of detector
 const spc_slr_cls &flx_slr_src, // I [obj] Solar spectrum
 const aer_cls &aer, // I [obj] Aerosol
 const psd_cls *psd_lst, // I [obj] Particle size distribution
 const prc_cmp *abs_fsh, // I [frc] Absorption efficiency
 const prc_cmp *cnc, // I [# m-3] Number concentration 
 const prc_cmp *mss, // I [kg] Mass 
 const prc_cmp *rds_ctr, // I [m] Radius at bin center
 const prc_cmp *xsa, // I [m2] Cross-sectional area
 const prc_cmp *ss_co_alb_fsh, // I [frc] Single scattering co-albedo
 const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp &abs_cff_mss, // I [m2 kg-1] Mass absorption coefficient
 const prc_cmp &mss_rsl, // I [kg m-3] Mass concentration resolved
 const prc_cmp &ss_co_alb); // I [frc] Single scattering co-albedo
// !aer_htg()

int // O [enm] Return success code
rnd_chm // [fnc] Raindrop chemistry
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *cnc, // I [# m-3] Number concentration 
 const prc_cmp *mss, // I [kg] Mass 
 const prc_cmp *rds_ctr, // I [m] Radius at bin center
 const prc_cmp &vmr_CO2); // [mlc mlc-1] Volume mixing ratio of CO2
// !rnd_chm()

int // O [enm] Return success code
plk_tbl_mk // [fnc] Build lookup-table for Planck function
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &tpt_min, // I [K] Minimum temperature in Planck-weight table
 const prc_cmp &tpt_max, // I [K] Maximum temperature in Planck-weight table
 const prc_cmp *wvn_grd, // I [cm-1] Wavenumber at band interfaces
 const long &wvn_nbr); // I [nbr] Number of wavenumber bands (interfaces minus one)
// !plk_tbl_mk()

int // O [enm] Return success code
rfl_frs // [fnc] Fresnel reflectance
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const long wvl_nbr, // I [nbr] Number of wavelength bins
 const std::complex<prc_cmp> *idx_rfr_1, // I [frc] Refractive index of transmitted medium
 const std::complex<prc_cmp> *idx_rfr_2, // I [frc] Refractive index of incident medium
 const prc_cmp &slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle
// !rfl_frs()

// Define inline'd functions in header so source is visible to calling files

#endif // NCO_HH  






