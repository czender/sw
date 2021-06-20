// $Id$ 

// Purpose: Surface flux physics

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <flx_sfc.hh> // Surface flux physics

#ifndef FLX_SFC_HH // Contents have not yet been inserted in current source file  
#define FLX_SFC_HH

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
#include <a2d.hh> // Two dimensional arrays
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <phys_cst.hh> // Physical constants NB: Live code, not a header
#include <lsm.hh> // Land Surface Model
#include <tdy.hh> // Atmospheric thermodynamics
#include <blm.hh> // Boundary layer meteorology

// Declare functions that have C++ linkages
int // O [rcd] Return success code
cnd_trm_soi_get // Thermal conductivity of soil
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnd_trm_soi, // O [W m-1 K-1] Soil thermal conductivity
 prc_cmp *lvl_dlt_snw, // O [m] Soil layer thickness including snow
 const prc_cmp *mss_frc_cly, // I [frc] Mass fraction clay 
 const prc_cmp *mss_frc_snd, // I [frc] Mass fraction sand
 const prc_cmp *snw_hgt, // I [m] Geometric bulk thickness of snow
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
 const prc_cmp *vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
 const prc_cmp *vwc_sat, // I [m3 m-3] Saturated volumetric water content
 const prc_cmp *vwc_sfc); // I [m3 m-3] Volumetric water content
// end cnd_trm_soi_get() prototype

int // O [rcd] Return success code
flx_sfc_lnd // Surface fluxes over land
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *vgt, // I [flg] "Vegetated" flag
 const prc_cmp *cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
 const prc_cmp *flx_SW_net_gnd, // I [W m-2] Solar flux absorbed by ground
 const prc_cmp *flx_SW_net_vgt, // I [W m-2] Solar flux absorbed by vegetation
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *lvl_dlt_snw, // I [m] Soil layer thickness including snow
 const prc_cmp *msv_gnd, // I [frc] Bare ground emissivity
 const prc_cmp *msv_snw, // I [frc] Snow emissivity
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 const prc_cmp tm_dlt, // I [s] Timestep
 const long *pnt_typ_idx, // I [idx] Plant type index (1..14)
 const long *soi_typ, // I [idx] LSM soil type (1..5)
 prc_cmp *cff_xch_heat, // O [frc] Exchange coefficient for heat transfer
 prc_cmp *cff_xch_mmn, // O [frc] Exchange coefficient for momentum transfer
 prc_cmp *cff_xch_mmn_ntr, // O [frc] Neutral drag coefficient hgt_mdp to z0m+zpd 
 prc_cmp *cff_xch_vpr, // O [frc] Exchange coefficient for vapor transfer
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm_ttl, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *flx_sns_gnd_ttl, // O [W m-2] Sensible heat flux to soil
 prc_cmp *flx_snw_mlt, // O [W m-2] Snow melt heat flux
 prc_cmp *mno_lng, // I/O [m] Monin-Obukhov length
 prc_cmp *msv_sfc, // O [frc] Surface (ground+snow) emissivity
 prc_cmp *rss_aer_heat_sfc, // O [s m-1] Aerodynamic resistance to heat transfer
 prc_cmp *rss_aer_mmn_sfc, // O [s m-1] Aerodynamic resistance to momentum transfer
 prc_cmp *rss_aer_vpr_sfc, // O [s m-1] Aerodynamic resistance to vapor transfer
 prc_cmp *tpt_aer, // O [K] "Aerodynamic" temperature at z=zpd+rgh_mmn
 prc_cmp *tpt_ash, // O [K] "Surface" temperature at z=zpd+rgh_heat
 prc_cmp *tpt_ash_p2m, // O [K] "Screen" temperature at z=zpd+rgh_heat+2m
 prc_cmp *tpt_gnd, // I/O [K] Ground temperature
 prc_cmp *tpt_msv, // O [K] Radiative emission temperature
 prc_cmp *tpt_vgt, // I/O [K] Vegetation temperature
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl); // O [kg m-1 s-2] Zonal wind stress
// end flx_sfc_lnd() prototype

int // O [rcd] Return success code
blm_mbl // Boundary layer meteorology over dust surfaces
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
 const prc_cmp *flx_SW_net, // I [W m-2] Solar flux absorbed by ground
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *lvl_dlt, // I [m] Soil layer thickness
 const prc_cmp *msv_gnd, // I [frc] Bare ground emissivity
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp *trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *mno_lng, // I/O [m] Monin-Obukhov length
 prc_cmp *tpt_gnd, // I/O [K] Ground temperature
 prc_cmp *wnd_frc); // O [m s-1] Surface friction velocity
// end blm_mbl() prototype

int // O [rcd] Return success code
blm_glb // Boundary layer meteorology for any surface type
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *oro, // I [frc] Orography
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling flux at surface
 prc_cmp *flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *hgt_zpd_dps, // O [m] Zero plane displacement
 prc_cmp *mno_lng_dps, // O [m] Monin-Obukhov length
 prc_cmp *rgh_mmn_dps, // O [m] Roughness length momentum
 prc_cmp *wnd_frc_dps, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_rfr_dps, // O [m s-1] Wind speed at reference height
 prc_cmp *wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl_dps); // O [kg m-1 s-2] Zonal wind stress
// end blm_glb() prototype

int // O [rcd] Return success code
blm_ice // Boundary layer meteorology over sea ice
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_ice, // I [flg] Sea ice flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl); // O [kg m-1 s-2] Zonal wind stress
// end blm_ice() prototype

int // O [rcd] Return success code
blm_lnd
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_lnd, // I [flg] Land flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Potential temperature
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl); // O [kg m-1 s-2] Zonal wind stress
// end blm_lnd() prototype

int // O [rcd] Return success code
blm_ocn // Boundary layer meteorology over ocean
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *flg_ocn, // I [flg] Ocean flag
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *prs_mdp, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp *tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
 const prc_cmp *tpt_sfc, // I [K] Surface temperature
 const prc_cmp *wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
 const prc_cmp *wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 prc_cmp *flx_LW_upw_sfc, // O [W m-2] Longwave upwelling flux at surface
 prc_cmp *flx_ltn, // O [W m-2] Latent heat flux to atmosphere
 prc_cmp *flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
 prc_cmp *flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
 prc_cmp *mno_lng, // O [m] Monin-Obukhov length
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 prc_cmp *wnd_frc, // O [m s-1] Surface friction velocity
 prc_cmp *wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
 prc_cmp *wnd_str_znl); // O [kg m-1 s-2] Zonal wind stress
// end blm_ocn() prototype

#endif // FLX_SFC_HH  
