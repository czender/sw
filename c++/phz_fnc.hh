// $Id$ 

// Purpose: Description (definition) of phase function utilities

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <phz_fnc.hh> // Phase functions

#ifndef PHZ_FNC_HH // Contents have not yet been inserted in current source file
#define PHZ_FNC_HH

// Standard C++ header files 
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <valarray> // STL valarray class template

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 

// 3rd party vendors
#include <netcdf.h> // netCDF C interface
#include <gsl/gsl_sf_bessel.h> // GNU Scientific Library special functions Bessel functions
#include <gsl/gsl_sf_legendre.h> // GNU Scientific Library special functions Legendre functions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <a2d.hh> // Two dimensional arrays
#include <mth.hh> // Mathematical utilities, constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Forward declarations

// Typedefs

// Prototype global functions with C++ linkages
prc_cmp // O [frc] Legendre polynomial of order n+1
lgn_Pnp1 // [fnc] Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= -1)
 prc_cmp abc); // [frc] Abscissa [-1,1]
// end lgn_Pnp1() prototype

prc_cmp // O [frc] First derivative of Legendre polynomial of order n+1
lgn_frs_drv_Pnp1 // [fnc] First derivative of Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= -1)
 prc_cmp abc); // [frc] Abscissa [-1,1]
// end lgn_frs_drv_Pnp1() prototype

prc_cmp // O [frc] Second derivative of Legendre polynomial of order n+1
lgn_scn_drv_Pnp1 // [fnc] Second derivative of Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= 0)
 prc_cmp abc); // [frc] Abscissa [-1,1]
// end lgn_scn_drv_Pnp1() prototype

int // O [enm] Return success code
ngl_grd_get // [fnc] Create angular grid
(const long ngl_nbr, // I [nbr] Number of angles
 const std::string ngl_sng, // I [sng] Angle grid type
 prc_cmp * const ngl, // O [rdn] Angle
 prc_cmp * const ngl_dgr, // O [dgr] Angle degrees
 prc_cmp * const ngl_dlt, // O [rdn] Width of angle bin
 prc_cmp * const ngl_wgt); // O [rdn] Weight of angle bin
// end ngl_grd_get() prototype

int // O [enm] Return success code
phz_fnc_mdl // [fnc] Phase function diagnostics module
(const int dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file 
 const int nc_out, // I [fl] netCDF file for output 
 const long lgn_nbr, // I [nbr] Order of phase function Legendre expansion
 const long ngl_nbr, // I [nbr] Number of polar angles in one hemisphere
 const long wvl_idx_dbg, // I [idx] Debugging wavelength bin
 const long wvl_nbr, // I [nbr] Number of output wavelength bands
 const prc_cmp ngl_dbg_dgr, // I [dgr] Debugging angle
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dgr, // I [dgr] Angle degrees
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 const prc_cmp * const ngl_wgt, // I [rdn] Weight of angle bin
 const a2d_cls<prc_cmp> &phz_fnc_ffc, // I [sr-1] Effective phase function (weighted over size distribution and sub-bands)
 const prc_cmp * const phz_fnc_dgn, // I [sr-1] Phase function at diagnostic size, wavelength
 const prc_cmp * const plz_dgn); // I [frc] Degree of linear polarization at diagnostic size, wavelength
// end phz_fnc_mdl() template

// Define inline'd functions in header so source is visible to calling files

/* http://www.storage-b.com/c/18

   legendre.h — C++ functions to evaluate Legendre polynomials
 
   Copyright (C) 2005 by James A. Chappell
 
   Permission is hereby granted, free of charge, to any person
   obtaining a copy of this software and associated documentation
   files (the "Software"), to deal in the Software without
   restriction, including without limitation the rights to use,
   copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following
   condition:
 
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS IN THE SOFTWARE. */

namespace Legendre{
  inline double P0(double x){return 1.0;} // n=0
  inline double P1(double x){return x;} // n=1
  inline double P2(double x){return ((3.0*x*x)-1.0)*0.5;} // n=2
  inline double Pn(unsigned int n, double x){ // Pn(x)
    if(n == 0){return P0(x);}
    else if(n == 1){return P1(x);}
    else if(n == 2){return P2(x);}
    
    if (x == 1.0){return 1.0;}

    if (x == -1.0){return ((n % 2 == 0) ? 1.0 : -1.0);}

    if ((x == 0.0) && (n % 2)){return 0.0;}

    /* We could simply do this:
       return (double(((2 * n) - 1)) * x * frs_drv_Pn(n - 1, x) -
       (double(n - 1)) * frs_drv_Pn(n - 2, x)) / (double)n ;
       but it could be slow for large n */
  
    double pnm1(P2(x)) ;
    double pnm2(P1(x)) ;
    double pn(pnm1) ;
    
    for(unsigned int l=3;l<=n;l++){ 
      pn=(((2.0*(double)l)-1.0)*x*pnm1-(((double)l-1.0)*pnm2))/(double)l;
      pnm2=pnm1;
      pnm1=pn;
    } // end loop over l
    
    return pn;
  } // end Pn()

  // First derivatives
  inline double frs_drv_P0(double x){return 0.0;} // n=0
  inline double frs_drv_P1(double x){return 1;} // n=1
  inline double frs_drv_P2(double x){return 3.0*x;} // n=2
  inline double frs_drv_Pn(unsigned int n, double x){ // Pn(x)
    if(n == 0){return frs_drv_P0(x);}
    else if(n == 1){return frs_drv_P1(x);}
    else if(n == 2){return frs_drv_P2(x);}
    
    /* We could simply do this:
       return (double(((2 * n) - 1)) * x * frs_drv_Pn(n - 1, x) -
       (double(n) * Pn(n - 2, x)) / ((double)n-1) ;
       but it could be slow for large n */

    double frs_drv_pnm1(frs_drv_P2(x)) ;
    double frs_drv_pnm2(frs_drv_P1(x)) ;
    double frs_drv_pn(frs_drv_pnm1) ;
    
    for(unsigned int l=3;l<=n;l++){ 
      frs_drv_pn=((2.0*(double)l-1.0)*x*frs_drv_pnm1-(double)l*frs_drv_pnm2)/((double)l-1.0);
      frs_drv_pnm2=frs_drv_pnm1;
      frs_drv_pnm1=frs_drv_pn;
    } // end loop over l
    
    return frs_drv_pn;
  } // end frs_drv_Pn()

  // Second derivatives
  inline double scn_drv_P0(double x){return 0.0;} // n=0
  inline double scn_drv_P1(double x){return 0.0;} // n=1
  inline double scn_drv_P2(double x){return 3.0;} // n=2
  inline double scn_drv_Pn(unsigned int n, double x){ // Pn(x)
    if(n == 0){return scn_drv_P0(x);}
    else if(n == 1){return scn_drv_P1(x);}
    else if(n == 2){return scn_drv_P2(x);}
    
    /* We could simply do this:
       return (double(((2 * n) - 1)) * x * scn_drv_Pn(n - 1, x) -
       (double(n+1) * Pn(n - 2, x)) / ((double)n-2) ;
       but it could be slow for large n */

    double scn_drv_pnm1(scn_drv_P2(x)) ;
    double scn_drv_pnm2(scn_drv_P1(x)) ;
    double scn_drv_pn(scn_drv_pnm1) ;
    
    for(unsigned int l=3;l<=n;l++){ 
      scn_drv_pn=((2.0*(double)l-1.0)*x*scn_drv_pnm1-((double)l+1.0)*scn_drv_pnm2)/((double)l-2.0);
      scn_drv_pnm2=scn_drv_pnm1;
      scn_drv_pnm1=scn_drv_pn;
    } // end loop over l
    
    return scn_drv_pn;
  } // end scn_drv_Pn()

} // end namespace Legendre

#endif // PHZ_FNC_HH  
