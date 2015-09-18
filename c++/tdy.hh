// $Id$ 

// Purpose: Atmospheric thermodynamics

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <tdy.hh> // Atmospheric thermodynamics

#ifndef TDY_HH // Contents have not yet been inserted in current source file  
#define TDY_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <phys_cst.hh> // Physical constants NB: Live code, not a header

// Namespaces
using phc::tpt_frz_pnt; // [K] Kelvin--Celsius scale offset, see, e.g., Bol80
using phc::eps_H2O; // (0.622) [frc] molec wgt vapor/molec wgt dry air
using phc::one_mns_eps_H2O; // (0.378) [frc] 1 - eps_H2O Constant for saturation specific humidity
using phc::gas_cst_H2O; // (461.65) [J kg-1 K-1]
using phc::prs_STP; // [Pa] Standard pressure

// Declare functions with C++ linkages
int // O [rcd] Return success code
dsvddt_H2O_PrK78 // [fnc] Derivative of saturation vapor density over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvddt_H2O); // O [Pa K-1] Derivative of saturation vapor density over planar water

int // O [rcd] Return success code
dsvpdt_H2O_PrK78 // [fnc] Derivative of saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O_lqd); // O [Pa K-1] Derivative of saturation vapor pressure over planar water

int // O [rcd] Return success code
dsvpdt_H2O_lqd_PrK78 // [fnc] Derivative of saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O_lqd); // O [Pa K-1] Derivative of saturation vapor pressure over planar liquid water

int // O [rcd] Return success code
dsvpdt_H2O_ice_PrK78 // [fnc] Derivative of saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O_ice); // O [Pa K-1] Derivative of saturation vapor pressure over planar ice water

int // O [rcd] Return success code
svp_H2O_PrK78 // [fnc] Saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O); // O [Pa] Saturation vapor pressure over planar water

int // O [rcd] Return success code
svp_H2O_PrK98 // [fnc] Saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O); // O [Pa] Saturation vapor pressure over planar water

int // O [rcd] Return success code
svp_H2O_lqd_PrK78 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd); // O [Pa] Saturation vapor pressure over planar liquid water

int // O [rcd] Return success code
svp_H2O_ice_PrK78 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice); // O [Pa] Saturation vapor pressure over planar ice water

int // O [rcd] Return success code
svp_H2O_ice_PrK98 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice); // O [Pa] Saturation vapor pressure over planar ice water

int // O [rcd] Return success code
svp_H2O_lqd_Her87 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd); // O [Pa] Saturation vapor pressure over planar liquid water

int // O [rcd] Return success code
svp_H2O_ice_Her87 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice); // O [Pa] Saturation vapor pressure over planar ice water

int // O [rcd] Return success code
sfc_tns_wtr_lqd_PrK78 // [fnc] Surface tension of liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *sfc_tns_wtr_lqd); // [J m-2], [N m-1] Surface tension of liquid water

int // O [rcd] Return success code
svp_H2O_lqd_Bol80 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd); // O [Pa] Saturation vapor pressure over planar liquid water

int // O [rcd] Return success code
act_cff_TaM94 // [fnc] Water activity
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *mss_frc_slt, // I [frc] Mass fraction of solute
 prc_cmp *act_cff); // [frc] Water activity

int // O [rcd] Return success code
ppr_H2O_get // [fnc] Partial pressure of H2O
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *ppr_H2O, // O [Pa] Partial pressure of H2O
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr); // I [kg kg-1] Specific humidity

int // O [rcd] Return success code
tpt_ptn_1000_get // [fnc] Potential temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp *tpt, // I [K] Midlayer temperature
 prc_cmp *tpt_ptn); // O [K] Potential temperature

int // O [rcd] Return success code
tpt_ptn_rfr_get // [fnc] Potential temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp *prs_rfr, // I [Pa] Reference pressure
 const prc_cmp *tpt, // I [K] Midlayer temperature
 prc_cmp *tpt_ptn); // O [K] Potential temperature

int // O [rcd] Return success code
tpt_vrt_get // [fnc] Virtual temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *tpt, // I [K] Midlayer temperature
 prc_cmp *tpt_vrt); // O [K] Virtual temperature

int // O [rcd] Return success code
q_H2O_vpr_get // [fnc] Specific humidity
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *ppr_H2O, // I [Pa] Partial pressure of H2O
 const prc_cmp *prs, // I [Pa] Pressure
 prc_cmp *q_H2O_vpr); // O [kg kg-1] Specific humidity

// Define inline'd functions in header so source is visible to calling files
// Order of declaration of inline'd functions is significant
// Inline'd functions cannot "see" inline'd functions declared later

inline prc_cmp // O [C] Bounded temperature celsius
tpt_bnd_cls_get(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: Convert temperature in Kelvin to temperature in Celsius, 
     suitable for use in fast PrK78 formulae for saturation.
     Since the range of validity of the Lowe and Ficke (1974) formulae is
     -50 C < T < 50 C, this function bounds the output temperature to that range. */
  return max_cpv(min_cpv(tpt-tpt_frz_pnt,50.0),-50.0); // [K] Temperature
} // end tpt_bnd_cls_get()

inline prc_cmp // O [Pa K-1] Derivative of saturation vapor pressure over planar liquid water
dsvpdt_H2O_lqd_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" derivative of saturation vapor pressure over planar liquid water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is derivative of saturation vapor pressure over planar liquid water in Pa K-1
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  const double cff[]={
    4.438099984e-01,
    2.857002636e-02,
    7.938054040e-04,
    1.215215065e-05,
    1.036561403e-07,
    3.532421810e-10,
    -7.090244804e-13
  };
  // Main code
  return 100.0*(cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls)))))); // [Pa K-1]
} // end dsvpdt_H2O_lqd_PrK78_fst_scl()

inline prc_cmp // O [Pa K-1] Derivative of saturation vapor pressure over planar ice water
dsvpdt_H2O_ice_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" derivative of saturation vapor pressure over planar ice water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is derivative of saturation vapor pressure over planar ice water in Pa K-1
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  const double cff[]={
    5.030305237e-01,
    3.773255020e-02,
    1.267995369e-03,
    2.477563108e-05,
    3.005693132e-07,
    2.158542548e-09,
    7.131097725e-12
  };
  // Main code
  return 100.0*(cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls)))))); // [Pa K-1]
} // end dsvpdt_H2O_ice_PrK78_fst_scl()

inline prc_cmp // O [Pa K-1] Derivative of saturation vapor pressure over planar liquid water
dsvpdt_H2O_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" derivative of saturation vapor pressure over planar water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is derivative of saturation vapor pressure over planar water in Pa K-1
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  // Main code
  return (tpt_cls > 0.0) ? dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_cls) : dsvpdt_H2O_ice_PrK78_fst_scl(tpt_cls); // [Pa]
} // end dsvpdt_H2O_PrK78_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar liquid water
svp_H2O_lqd_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar liquid water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is saturation vapor pressure over planar liquid water in Pa
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  const double cff[]={
    6.107799961,
    4.436518521e-1,
    1.428945805e-2,
    2.650648471e-4,
    3.031240396e-6,
    2.034080948e-8,
    6.136820929e-11
  };
  // Main code
  return 100.0*(cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls)))))); // [Pa]
} // end svp_H2O_lqd_PrK78_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar ice water
svp_H2O_ice_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar ice water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is saturation vapor pressure over planar ice water in Pa
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  const double cff[]={
    6.109177956,
    5.034698970e-1,
    1.886013408e-2,
    4.176223716e-4,
    5.824720280e-6,
    4.838803174e-8,
    1.838826904e-10
  };
  // Main code
  return 100.0*(cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls)))))); // [Pa]
} // end svp_H2O_ice_PrK78_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar ice water
svp_H2O_ice_PrK98_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar ice water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is saturation vapor pressure over planar ice water in Pa
     Taken from Lowe and Ficke (1974) and Lowe (1981) as reported in PrK98 p. 854
     PrK98 values for ice differ from PrK78, values for liquid are identical */
  const double cff[]={
    6.10690449,
    5.02660639e-1,
    1.87743264e-2,
    4.13476180e-4,
    5.72333773e-6,
    4.71651246e-8,
    1.78086695e-10
  };
  // Main code
  return 100.0*(cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls)))))); // [Pa]
} // end svp_H2O_ice_PrK98_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar liquid water
svp_H2O_lqd_Her87_fst_scl(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar liquid water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 243.15 K < T < 308.15 K
     Output is saturation vapor pressure over planar liquid water in Pa
     Method is based on the Magnus equation, with coefficients from Herbert (1987) as reported in PrK98 p. 854 */
  // Main code
  return 100.0*6.1070*std::exp(17.15*(tpt-tpt_frz_pnt)/(tpt-38.25)); // [Pa] Her87
} // end svp_H2O_lqd_Her87_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar ice water
svp_H2O_ice_Her87_fst_scl(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar ice water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 243.15 K < T < 308.15 K
     Output is saturation vapor pressure over planar ice water in Pa
     Method is based on the Magnus equation, with coefficients from Herbert (1987) as reported in PrK98 p. 854 */
  // Main code
  return 100.0*6.1064*std::exp(21.88*(tpt-tpt_frz_pnt)/(tpt-7.65)); // [Pa] Her87
} // end svp_H2O_ice_Her87_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar liquid water
svp_H2O_Her87_fst_scl(const prc_cmp tpt) // I [K] Temperature celsius
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 243.15 K < T < 308.15 K
     Output is saturation vapor pressure over planar water in Pa
     Method is based on the Magnus equation, with coefficients from Herbert (1987) as reported in PrK98 p. 854 */

  // Main code
  return (tpt > tpt_frz_pnt) ? svp_H2O_lqd_Her87_fst_scl(tpt) : svp_H2O_ice_Her87_fst_scl(tpt); // [Pa]
} // end svp_H2O_Her87_fst_scl()

inline prc_cmp // O [Pa] Saturation vapor pressure over planar liquid water
svp_H2O_PrK78_fst_scl(const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" saturation vapor pressure over planar water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Celsius in range -50 C < T < 50 C
     Output is saturation vapor pressure over planar water in Pa
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  // Main code
  return (tpt_cls > 0.0) ? svp_H2O_lqd_PrK78_fst_scl(tpt_cls) : svp_H2O_ice_PrK78_fst_scl(tpt_cls); // [Pa]
} // end svp_H2O_PrK78_fst_scl()

inline prc_cmp // O [m s-1] Speed of sound in water
spd_snd_wtr_fst_scl
(const prc_cmp dpt, // I [m] Depth
 const prc_cmp sln_ppt, // I [ppt] Salinity
 const prc_cmp tpt_cls) // I [C] Temperature celsius
{
  /* Purpose: "Fast scalar" speed of sound in salt water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, salinity conversion, and included files.
     http://en.wikipedia.org/wiki/Speed_of_sound#Seawater
     c(T, S, z) = a1 + a2T + a3T^2 + a4T^3 + a5(S - 35) + a6z + a7z^2 + a8T(S - 35) + a9Tz^3 */

  const double cff[]={
    0.0,
    1448.96,
    4.591,
    -5.304e-2,
    2.374e-4,
    1.340,
    1.630e-2,
    1.675e-7,
    -1.025e-2,
    -7.139e-13 
  };
  // Main code
  return cff[1]+cff[2]*tpt_cls+cff[3]*tpt_cls*tpt_cls+cff[4]*tpt_cls*tpt_cls*tpt_cls+cff[5]*(sln_ppt-35.0)+cff[6]*dpt+cff[7]*dpt*dpt+cff[8]*tpt_cls*(sln_ppt-35.0)+cff[9]*tpt_cls*dpt*dpt*dpt; // [m s-1]
} // end spd_snd_wtr_fst_scl()

inline prc_cmp // O [kg kg-1] Specific humidity
q_H2O_vpr_fst_scl_get
(const prc_cmp ppr_H2O, // I [Pa] Partial pressure of H2O
 const prc_cmp prs) // I [Pa] Pressure
{
  /* Purpose: Given pressure and moisture values, return vapor mixing ratio
     Requires: <phys_cst.hh>
     Output
     Local
     Main code */
  return eps_H2O*ppr_H2O/(prs-one_mns_eps_H2O*ppr_H2O); // [kg kg] Specific humidity
} // end q_H2O_vpr_get()

inline prc_cmp // O [kg m-3 K-1] Derivative of saturation vapor density over planar water
dsvddt_H2O_PrK78_fst_scl // [fnc] Derivative of saturation vapor density over planar water
(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: "Fast scalar" derivative of saturation vapor density over planar water
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 223.15 K < T < 323.15 K
     Output is derivative of saturation vapor density over planar water in kg m-3 K-1
     Requires: <phys_cst.hh> for tpt_frz_pnt, gas_cst_H2O */
  // Main code
  const prc_cmp tpt_cls(tpt-tpt_frz_pnt); // [K] -> [C] Temperature
  return -gas_cst_H2O*svp_H2O_PrK78_fst_scl(tpt_cls)+(1.0/(gas_cst_H2O*tpt))*dsvpdt_H2O_PrK78_fst_scl(tpt_cls); // [kg m-3 K-1] Derivative of saturation vapor density over planar water
} // end dsvddt_H2O_PrK78_fst_scl()

inline prc_cmp // O [m2 s-1] Diffusivity of water vapor in air
dff_H2O_air_fst_scl // [fnc] Diffusivity of water vapor in air
(const prc_cmp tpt, // I [K] Temperature
 const prc_cmp prs) // I [Pa] Pressure
{
  /* Purpose: "Fast scalar" diffusivity of water vapor in air
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 233.15 K < T < 313.15 K
     Output is diffusivity of water vapor in air in m2 s-1
     Requires: <phys_cst.hh> for tpt_frz_pnt, prs_STP */
  // Main code
  assert(tpt > 233.15 && tpt < 313.15);
  return 2.11e-5*(prs_STP/prs)*std::pow(tpt/tpt_frz_pnt,1.94); // [m2 s-1] Diffusivity of water vapor in air PrK98 p. 503 (13-3)
} // end dff_H2O_air_fst_scl()

inline prc_cmp // O [kg m-1 s-1] Dynamic viscosity of atmosphere
vsc_dyn_atm_fst_scl // [fnc] Dynamic viscosity of atmosphere
(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: "Fast scalar" diffusivity of water vapor in air
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin in range 233.15 K < T < 313.15 K
     Output is dynamic viscosity of atmosphere in kg m-1 s-1 */
  // Main code
  // Approximations are given by RoY94 p. 102, PrK78 p. 323 (10-107), PrK98 p. 417 (10-141)
  /*  const prc_cmp tpt_cls(tpt-tpt_frz_pnt); // [K] -> [C] Temperature
      return (tpt_cls > 0.0) ? (1.718+0.0049*tpt_cls)*1.0e-5 : (1.718+0.0049*tpt_cls-1.2e-5*tpt_cls*tpt_cls)*1.0e-5; // [kg m-1 s-1] Dynamic viscosity of atmosphere PrK98 p. 417 (10-141) 
      return vsc_dyn_atm=1.832e-4*std::pow(tpt/296.16,1.5)*(296.16+120.0)/(tpt+120.0); // [g cm-1 s-1] CGS (Source: Probably FTV89 or PrK78) */
  return 1.72e-5*std::pow(tpt/273.0,1.5)*393.0/(tpt+120.0); // [kg m-1 s-1] Dynamic viscosity of atmosphere RoY94 p. 102
} // end vsc_dyn_atm_fst_scl()

inline prc_cmp // O [W m-1 K-1] Thermal conductivity of dry air
cnd_trm_dry_air_fst_scl // [fnc] Thermal conductivity of dry air
(const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: "Fast scalar" thermal conductivity of dry air
     In order to maximize chances of compiler inlining this function,
     it avoids error checking, temperature conversion, and included files
     Input must be temperature in Kelvin
     Output is thermal conductivity of dry air in W m-1 K-1 */
  // Main code
  // PrK98 p. 508 discus weighting cnd_trm_dry_air and cnd_trm_vpr to obtain cnd_trm_air
  //  return 100.0*joules_per_calorie*(5.69+0.017*(tpt-tpt_frz_pnt))*1.0e-5; // [W m-1 K-1] Thermal conductivity of dry air PrK98 p. 508 (13-18a)
  return 0.0043794+7.12e-5*tpt; // [W m-1 K-1] Thermal conductivity of dry air PrK98 p. 508 (13-18a)
} // end cnd_trm_dry_air_fst_scl()

#endif // TDY_HH  






