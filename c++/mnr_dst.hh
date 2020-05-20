// $Id$ 

// Purpose: Mineral dust physics

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mnr_dst.hh> // Mineral dust physics

#ifndef MNR_DST_HH // Contents have not yet been inserted in current source file  
#define MNR_DST_HH

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

int // O [rcd] Return success code
flx_mss_hrz_slt_ttl_Whi79_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 prc_cmp *flx_mss_hrz_slt_ttl, // O [kg m-1 s-1] Horizontal mass flux of saltators
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_frc_thr_slt); // I [m s-1] Threshold friction speed for saltation

int // O [rcd] Return success code
flx_mss_vrt_dst_ttl_MaB95_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *dst_slt_flx_rat_ttl, // O [frc] Ratio of vertical dust flux to horizontal saltator flux
 const prc_cmp *flx_mss_hrz_slt_ttl, // I [kg m-1 s-1] Horizontal mass flux of saltators
 prc_cmp *flx_mss_vrt_dst_ttl, // O [kg m-2 s-1] Total vertical mass flux of dust
 const prc_cmp *mss_frc_cly); // I [frc] Mass fraction clay 

int // O [rcd] Return success code
wnd_frc_thr_get
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_prt, // I [kg m-3] Density of particle
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *ryn_nbr_frc_thr, // O [frc] Threshold friction Reynolds number
 prc_cmp *ryn_nbr_frc_thr_prx, // O [frc] Threshold friction Reynolds number approximation
 prc_cmp *wnd_frc_thr, // O [m s-1] Threshold friction speed
 prc_cmp *wnd_frc_thr_prx); // O [m s-1] Threshold friction speed approximation

int // O [rcd] Return success code
frc_thr_ncr_drg_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *frc_thr_ncr_drg, // O [frc] Factor by which surface roughness increases threshold friction velocity
 prc_cmp *wnd_frc_fsh_frc, // O [frc] Efficient fraction of wind friction
 const prc_cmp rgh_mmn_mbl, // I [m] Roughness length momentum for erodible surfaces
 const prc_cmp rgh_mmn_smt); // I [m] Smooth roughness length

int // O [rcd] Return success code
frc_thr_ncr_wtr_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *frc_thr_ncr_wtr, // O [frc] Factor by which moisture increases threshold friction velocity
 prc_cmp *vwc_thr, // O [m3 m-3] Threshold volumetric water content to affect mobilization
 const prc_cmp *mss_frc_cly, // I [frc] Mass fraction of clay
 const prc_cmp *vwc_sfc); // I [frc] Volumetric water content

int // O [rcd] Return success code
wnd_frc_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *wnd_frc_mbl, // I [m s-1] Surface friction velocity
 prc_cmp *wnd_frc_slt, // O [m s-1] Saltating friction velocity
 prc_cmp *wnd_frc_slt_dlt, // O [m s-1] Friction velocity increase from saltation
 const prc_cmp *wnd_rfr_mbl, // I [m s-1] Wind speed at reference height
 const prc_cmp *wnd_rfr_thr_slt); // I [m s-1] Threshold 10 m wind speed for saltation

int // O [rcd] Return success code
wnd_frc_thr_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // [kg m-3] Midlayer density
 prc_cmp *wnd_frc_thr_slt); // O [m s-1] Threshold friction speed for saltation

int // O [rcd] Return success code
wnd_rfr_thr_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_frc_thr_slt, // I [m s-1] Threshold friction speed for saltation
 const prc_cmp *wnd_mdp, // I [m s-1] Surface layer mean wind speed
 const prc_cmp *wnd_rfr, // I [m s-1] Wind speed at reference height
 prc_cmp *wnd_rfr_thr_slt); // O [m s-1] Threshold 10 m wind speed for saltation

#endif // MNR_DST_HH  
