// $Id$ 

// Purpose: Sea salt physics

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <ssl.hh> // Sea salt physics

#ifndef SSL_HH // Contents have not yet been inserted in current source file  
#define SSL_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <valarray> // STL valarray class template

// Standard C headers 
#include <cstdio> // stderr, FILE, NULL, etc.
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <mth_cst.hh> // Mathematical constants, cst_M_PIl...
#include <phys_cst.hh> // Physical constants NB: Live code, not a header
#include <blm.hh> // Boundary layer meteorology

// Define inline'd functions in header so source is visible to calling files

int // O [rcd] Return success code
cst_C1C2C3_And98_get // [fnc] 
(const prc_cmp wnd_10m, // I [m s-1] Wind speed at 10 m reference height
 prc_cmp &cst_C1_And98, // O [frc] Factor in direct production scheme And98
 prc_cmp &cst_C2_And98, // O [frc] Factor in direct production scheme And98
 prc_cmp &cst_C3_And98); // O [frc] Factor in direct production scheme And98
// end cst_C1C2C3_And98_get() prototype

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of 10 micron sea salt particles at 80% RH
flx_vrt_nbr_ssl_10mcr_SPC93_get // [fnc] Evaluate number flux of 10 micron sea salt particles using SPC93
(prc_cmp wnd_10m); // I [m s-1] Wind speed at 10 m reference height
// end flx_vrt_nbr_ssl_10mcr_SPC93_get() prototype

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of sea salt
flx_vrt_nbr_ssl_SPC93_get // [fnc] Sea salt generation function of SPC93
(prc_cmp rds_prt, // I [m] Particle radius at ocean surface = 100% RH
 prc_cmp wnd_10m); // I [m s-1] Wind speed at 10 m reference height
// end flx_vrt_nbr_ssl_SPC93_get() prototype

prc_cmp // O [# m-2 s-1 m-1] Spectral vertical number flux of sea salt
flx_vrt_nbr_ssl_VDB01_get // [fnc] Sea salt generation function of VDB01
(prc_cmp rds_prt, // I [m] Particle radius at ocean surface = 100% RH
 prc_cmp wnd_10m); // I [m s-1] Wind speed at 10 m reference height
// end flx_vrt_nbr_ssl_VDB01_get() prototype

int flx_mss_vrt_ssl_get
(const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *sz_ctr, // I [m] Size at bin center
 const prc_cmp *sz_dlt, // I [m] Width of size bin
 const prc_cmp &wnd_rfr, // I [m s-1] Wind speed at reference height
 prc_cmp *flx_mss_vrt_ssl_drc, // O [kg m-2 s-1 m-1] Direct sea-salt production And98
 prc_cmp *flx_mss_vrt_ssl_ndr, // O [kg m-2 s-1 m-1] Indirect sea-salt production MSD86
 prc_cmp *flx_mss_vrt_ssl, // O [kg m-2 s-1 m-1] Vertical mass flux of sea-salt
 prc_cmp *ssl_flx_ndr_drc_rat); // O [m-1] Ratio of direct to indirect sea-salt mass flux
// end flx_mss_vrt_ssl_get() prototype

#endif // SSL_HH  
