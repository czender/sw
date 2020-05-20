// $Id$ 

// Purpose: Prototypes, typedefs, and global variables for libcsm_c++  

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* libcsm_c++.hh headers depend on libcsz_c++.hh and libnco_c++ headers
   as well as dbg.hh, htrn_c++.hh, and phys_cst.hh */

// Usage:
// #include <libcsm_c++.hh> // Climate systems model library

#ifndef LIBCSM_CCC_HH // Contents have not yet been inserted in current source file  
#define LIBCSM_CCC_HH

// Personal headers
#include <blm.hh> // Boundary layer meteorology
#include <chm.hh> // Chemical utilities
#include <flx_sfc.hh> // Surface flux physics
#include <idx_rfr.hh> // Refractive indices
#include <lbl.hh> // Line-by-line utilities
#include <mie_sln.hh> // Mie scattering solutions and utilities
#include <mie_MaS99.hh> // Mie scattering solutions from MaS99
#include <mie_Sor01.hh> // Mie scattering solutions from Sor01
#include <mie_Wis79.hh> // Mie scattering solutions from Wis79
#include <tdy.hh> // Atmospheric thermodynamics
#include <aer.hh> // Aerosol physics
#include <flp.hh> // Floating point utilities, constants
#include <mnr_dst.hh> // Mineral dust physics
#include <nco.hh> // netCDF utilities
#include <pdf.hh> // Probability density functions
#include <phz_fnc.hh> // Phase functions
#include <psd.hh> // Particle size distributions
#include <rt.hh> // Radiative transfer
#include <spc_bbd.hh> // Blackbody spectra
#include <spc_slr.hh> // Solar spectra
#include <ssl.hh> // Sea salt physics
#include <sz_grd.hh> // Size grids
#include <wvl_grd.hh> // Wavelength grids

#endif // LIBCSM_CCC_HH  


