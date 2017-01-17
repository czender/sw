// $Id$ 

// Purpose:

/* Copyright (C) 2005--2017 Charlie Zender, Jorge Talamantes
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie_Sor01.hh> // Mie scattering solutions from MaS99

#ifndef MIE_SOR01_HH // Contents have not yet been inserted in current source file
#define MIE_SOR01_HH

// Standard C++ headers 
#include <complex> // Standard C++ complex class
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr

// Standard C headers 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdlib> // strtod, strtol, malloc, getopt  

// 3rd party vendors

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Namespaces

// Typedefs

// Forward declarations

// Prototype functions with C++ linkages

prc_cmp // O [m2] Aggregate cross-section for absorption (sigma_abs^agg)
xsx_abs_grg_Sor01 // [fnc] Compute aggregate cross-section for absorption (sigma_abs^agg)
(int mnm_nbr, // [nbr] Number of monomers in aggregate (N)
 prc_cmp rds_mnm, // [m] Monomer radius (a)
 prc_cmp wvn, // [m-1] Radiation wavenumber (k)
 std::complex<prc_cmp> idx_rfr); // [frc] Refractive index (m)
// end xsx_abs_grg_Sor01() prototype

prc_cmp // [m2] Aggregate cross-section for scattering (sigma_sca^agg)
xsx_sct_grg_Sor01 // [fnc] Compute aggregate cross-section for scattering (sigma_sca^agg)
(int mnm_nbr, // [nbr] Number of monomers in aggregate (N)
 prc_cmp dmn_frc, // [frc] Fractal dimension of aggregate (D)
 prc_cmp rds_gyr, // [m] Radius of gyration of aggregate (R_g)
 prc_cmp rds_mnm, // [m] Monomer radius (a)
 prc_cmp wvn, // [m-1] Radiation wavenumber (k)
 std::complex<prc_cmp> idx_rfr); // [frc] Refractive index (m)
// end xsx_sct_grg_Sor01() prototype

int // O [enm] Return success code
mie_Sor01_tst // [fnc] Test function for Sor01 physics
(void);
// end mie_Sor01_tst() prototype

#endif // MIE_SOR01_HH  
  
