// $Id$ 

// Purpose: Description (definition) of Markel and Shalaev (1999) Mie scattering solutions

/* Copyright (C) 2005--present Charlie Zender, Jorge Talamantes, Vadim Markel, and Vladimir Shalaev
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie_MaS99.hh> // Mie scattering solutions from MaS99

#ifndef MIE_MAS99_HH // Contents have not yet been inserted in current source file
#define MIE_MAS99_HH

// Standard C++ headers 
#include <complex> // Standard C++ complex class
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr

// Standard C headers 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdlib> // strtod, strtol, malloc, getopt  

// 3rd party vendors
#include <gsl/gsl_sf_bessel.h> // GNU Scientific Library special functions Bessel functions
#include <gsl/gsl_integration.h> //  GNU Scientific Library numerical integration

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Namespaces

// Typedefs

// Forward declarations

// Prototype functions with C++ linkages
double // O [fnc] Absorption enhancement for inclusions in weakly-absorbing spheres
mie_sph_abs_fct_MaS99 // [fnc] Absorption enhancement for inclusions in weakly-absorbing spheres
(const double sz_prm, // I [m m-1] Size parameter
 const double idx_rfr_mdm_rl, // I [frc] Refractive index of weakly-absorbing sphere, real component
 const double dmn_frc); // I [frc] Fractal dimensionality of inclusions
// end mie_sph_abs_fct_MaS99() prototype

int // O [enm] Return success code
mie_sph_abs_fct_MaS99_tst // [fnc] Test mie_sph_abs_fct_MaS99()
(void);
// end mie_sph_abs_fct_MaS99_tst() prototype
  
#endif // MIE_MAS99_HH  
  
