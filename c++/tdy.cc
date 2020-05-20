// $Id$ 

// Purpose: Thermodynamic utilities for C++ programs

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <tdy.hh> // Atmospheric thermodynamics

using phc::eps_H2O_rcp_m1; // (0.60777) [frc] Constant for virtual temperature
using phc::kappa_dry_air; // (0.286 = 2/7) [frc] Constant in potential temperature IrG81 p. 25, Tre922 p. 72 
using phc::prs_1000; // (100000.0) [Pa] Reference pressure for potential temperature
using phc::tpt_frz_pnt; // (273.15) [K] Kelvin--Celsius scale offset, see, e.g., Bol80

int // O [rcd] Return success code
svp_H2O_PrK78 // [fnc] Saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O) // O [Pa] Saturation vapor pressure over planar water
{
  /* Purpose: Compute saturation vapor pressure over planar water
     Calls appropriate routine for liquid or ice depending on temperature
     Requires: <phys_cst.hh> for tpt_frz_pnt */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // [idx] Counting index
  double tpt_cls; // [C] Temperature
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    if(tpt_cls > 0.0){
      svp_H2O[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    }else{
      svp_H2O[lon_idx]=svp_H2O_ice_PrK78_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    } // endif
  } // end loop over lon
  return rcd;
} // end svp_H2O_PrK78()

int // O [rcd] Return success code
svp_H2O_PrK98 // [fnc] Saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O) // O [Pa] Saturation vapor pressure over planar water
{
  /* Purpose: Compute saturation vapor pressure over planar water
     Calls appropriate routine for liquid or ice depending on temperature
     Requires: <phys_cst.hh> for tpt_frz_pnt */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    if(tpt_cls > 0.0){
      // PrK98 = PrK78 for liquid
      svp_H2O[lon_idx]=svp_H2O_lqd_PrK78_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    }else{
      svp_H2O[lon_idx]=svp_H2O_ice_PrK98_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    } // endif
  } // end loop over lon
  return rcd;
} // end svp_H2O_PrK98()

int // O [rcd] Return success code
svp_H2O_Her87 // [fnc] Saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O) // O [Pa] Saturation vapor pressure over planar water
{
  /* Purpose: Compute saturation vapor pressure over planar water
     Calls appropriate routine for liquid or ice depending on temperature
     Requires: <phys_cst.hh> for tpt_frz_pnt */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(tpt[lon_idx] > tpt_frz_pnt){
      svp_H2O[lon_idx]=svp_H2O_lqd_Her87_fst_scl(tpt[lon_idx]); // [Pa] Saturation vapor pressure over planar water
    }else{
      svp_H2O[lon_idx]=svp_H2O_ice_Her87_fst_scl(tpt[lon_idx]); // [Pa] Saturation vapor pressure over planar water
    } // endif
  } // end loop over lon
  return rcd;
} // end svp_H2O_Her87()

int // O [rcd] Return success code
svp_H2O_lqd_PrK78 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd) // O [Pa] Saturation vapor pressure over planar liquid water
{
  // Purpose: Compute saturation vapor pressure over planar liquid water
  // Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625
  // Range of validity is -50 C < T < 50 C
  // Requires: <phys_cst.hh> for tpt_frz_pnt
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double svp_H2O_lqd_dbl; // [Pa]
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
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    if(tpt_cls < -50.0 || tpt_cls > 50.0) rcd+=1; // Out of range error
    if(tpt_cls > -50.0){
      svp_H2O_lqd_dbl=cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls))))); // [mb]
      svp_H2O_lqd[lon_idx]=svp_H2O_lqd_dbl*100.0; // [mb] --> [Pa]
    }else{
      rcd+=svp_H2O_lqd_Bol80(lon_nbr,tpt,svp_H2O_lqd);
    } // endif
  } // end loop over lon
  return rcd;
} // end svp_H2O_lqd_PrK78()

int // O [rcd] Return success code
svp_H2O_ice_PrK78 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice) // O [Pa] Saturation vapor pressure over planar ice water
{
  // Purpose: Compute saturation vapor pressure over planar ice water
  // Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625
  // Range of validity is -50 C < T < 50 C
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double svp_H2O_ice_dbl; // [Pa]
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
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    svp_H2O_ice_dbl=cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls))))); // [mb]
    svp_H2O_ice[lon_idx]=svp_H2O_ice_dbl*100.0; // [mb] --> [Pa]
    if(tpt_cls < -50.0 || tpt_cls > 50.0) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end svp_H2O_ice_PrK78()

int // O [rcd] Return success code
svp_H2O_ice_PrK98 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice) // O [Pa] Saturation vapor pressure over planar ice water
{
  /* Purpose: Compute saturation vapor pressure over planar ice water
     Taken from Lowe and Ficke (1974) and Lowe (1981) as reported in PrK98 p. 854
     Range of validity is -50 C < T < 50 C 
     PrK98 values for ice differ from PrK78, values for liquid are identical */
   
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double svp_H2O_ice_dbl; // [Pa]
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
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    svp_H2O_ice_dbl=cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls))))); // [mb]
    svp_H2O_ice[lon_idx]=svp_H2O_ice_dbl*100.0; // [mb] --> [Pa]
    if(tpt_cls < -50.0 || tpt_cls > 50.0) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end svp_H2O_ice_PrK98()

int // O [rcd] Return success code
sfc_tns_wtr_lqd_PrK78 // [fnc] Surface tension of liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *sfc_tns_wtr_lqd) // [J m-2], [N m-1] Surface tension of liquid water
{
  /* Purpose: Compute surface tension of liquid water
     Input temperature in degrees kelvin
     Surface tension returned in [J m-2] = [N m-1]
     Taken from PrK78 p. 104 (5-12)
     Range of validity is 243 < T < 313 K */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  std::string sbr_nm("sfc_tns_wtr_lqd_PrK78"); // [sng Subroutine name
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    sfc_tns_wtr_lqd[lon_idx]=-1.55e-5*(tpt[lon_idx]-273.15)+7.610e-3; // [J m-2], [N m-1]
    if(tpt[lon_idx] < 243.0 || tpt[lon_idx] > 313.0) // Out of range error
      wrn_prn(sbr_nm,"Outside range of validity of 243 < T < 313 K for surface tension of liquid water");
  }  // end loop over lon
  return rcd;
} // end sfc_tns_wtr_lqd_PrK78()

int // O [rcd] Return success code
svp_H2O_lqd_Bol80 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd) // O [Pa] Saturation vapor pressure over planar liquid water
{
  /* Purpose: Compute saturation vapor pressure over planar liquid water
     Taken from Bol80
     This Bol80 parameterization has an accuracy of 0.1% for -30 < tpt_cls < 35 C */
  // Requires: <phys_cst.hh> for tpt_frz_pnt
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double svp_H2O_lqd_dbl; // [Pa]
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    svp_H2O_lqd_dbl=6.112*std::exp(17.67*tpt_cls/(tpt_cls+243.5)); // [mb] Bol80
    svp_H2O_lqd[lon_idx]=svp_H2O_lqd_dbl*100.0; // [mb] --> [Pa]
    if(tpt_cls < -30.0 || tpt_cls > 35.0) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end svp_H2O_lqd_Bol80()

int // O [rcd] Return success code
svp_H2O_lqd_Her87 // [fnc] Saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_lqd) // O [Pa] Saturation vapor pressure over planar liquid water
{
  /* Purpose: Compute saturation vapor pressure over planar liquid water
     Method is based on the Magnus equation, with coefficients from Herbert (1987) as reported in PrK98 p. 854 */
  // Requires: <phys_cst.hh> for tpt_frz_pnt
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double svp_H2O_lqd_dbl; // [Pa]
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    svp_H2O_lqd_dbl=6.1070*std::exp(17.15*(tpt[lon_idx]-tpt_frz_pnt)/(tpt[lon_idx]-38.25)); // [mb] Her87
    svp_H2O_lqd[lon_idx]=svp_H2O_lqd_dbl*100.0; // [mb] --> [Pa]
    // Range of validity of Her87 is unknown
    // Use same range as Bol80 since both are based on Magnus equation
    if(tpt[lon_idx] < 243.15 || tpt[lon_idx] > 308.15) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end svp_H2O_lqd_Her87()

int // O [rcd] Return success code
svp_H2O_ice_Her87 // [fnc] Saturation vapor pressure over planar ice water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *svp_H2O_ice) // O [Pa] Saturation vapor pressure over planar ice water
{
  /* Purpose: Compute saturation vapor pressure over planar ice water
     Method is based on the Magnus equation, with coefficients from Herbert (1987) as reported in PrK98 p. 854 */
  // Requires: <phys_cst.hh> for tpt_frz_pnt
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double svp_H2O_ice_dbl; // [Pa]
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    svp_H2O_ice_dbl=6.1064*std::exp(21.88*(tpt[lon_idx]-tpt_frz_pnt)/(tpt[lon_idx]-7.65)); // [mb] Her87
    svp_H2O_ice[lon_idx]=svp_H2O_ice_dbl*100.0; // [mb] --> [Pa]
    // Range of validity of Her87 is unknown
    // Use same range as Bol80 since both are based on Magnus equation
    if(tpt[lon_idx] < 243.15 || tpt[lon_idx] > 308.15) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end svp_H2O_ice_Her87()

int // O [rcd] Return success code
act_cff_TaM94 // [fnc] Water activity
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *mss_frc_slt, // I [frc] Mass fraction of solute
 prc_cmp *act_cff) // [frc] Water activity
{
  // Purpose: Compute water activity of solution
  // Input mass fraction of solute
  // Water activity returned in [frc]
  // Taken from TaM94 p. 18805 (2)
  // Range of validity is 0 < mss_frc_slt < 78%
// Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  prc_cmp mss_frc_slt_pct; // [pct] Mass fraction of solute in percent
  long lon_idx; // Counting index
  const prc_cmp cff_NH4_2_SO4[]={
    0.0,-2.715e-3,3.113e-5,-2.336e-6,1.412e-8
  };
  // Main code
  std::string sbr_nm("act_cff_TaM94");
  // Assign generic array to particular solute requested
  const prc_cmp *cff(cff_NH4_2_SO4); // Generic coefficients
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    mss_frc_slt_pct=mss_frc_slt[lon_idx]*100.0; // [pct] Mass fraction of solute in percent
    if(mss_frc_slt_pct < 78.0){
      const prc_cmp mss_frc_slt_pct_sqr(mss_frc_slt_pct*mss_frc_slt_pct); // [pct] Mass fraction of solute in percent, squared
      act_cff[lon_idx]=1.0+cff[1]*mss_frc_slt_pct+cff[2]*mss_frc_slt_pct_sqr+cff[3]*mss_frc_slt_pct*mss_frc_slt_pct_sqr+cff[4]*mss_frc_slt_pct_sqr*mss_frc_slt_pct_sqr; // [frc]
    }else{ // Out of range of validity
      rcd+=1; // Out of range error
      err_prn(sbr_nm,"Outside range of validity of molality");
    } // endif
  }  // end loop over lon
  return rcd;
} // end act_cff_TaM94()

int // O [rcd] Return success code
tpt_ptn_1000_get // [fnc] Potential temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *tpt_ptn) // O [K] Potential temperature
{
  // Purpose: Compute potential temperature relative to 1000 mb
  // Potential temperature is the temperature a parcel would have if moved adiabatically from prs to 1000 mb
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_ptn[lon_idx]=tpt[lon_idx]*std::pow(prs_1000/prs[lon_idx],kappa_dry_air); // [K] 
  }  // end loop over lon
  return rcd;
} // end tpt_ptn_1000_get()

int // O [rcd] Return success code
tpt_ptn_rfr_get // [fnc] Potential temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp *prs_rfr, // I [Pa] Reference pressure
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *tpt_ptn) // O [K] Potential temperature
{
  // Purpose: Compute potential temperature relative to supplied reference pressure
  // Potential temperature is the temperature a parcel would have if moved adiabatically from prs_mdp to prs_rfr
  // Normally, this routine is called with prs_rfr = prs_sfc 
  // The potential temperature returned is relative to prs_rfr, not to 1000 mb
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_ptn[lon_idx]=tpt[lon_idx]*std::pow(prs_rfr[lon_idx]/prs_mdp[lon_idx],static_cast<prc_cmp>(kappa_dry_air)); // [K] 
  }  // end loop over lon
  return rcd;
} // end tpt_ptn_rfr_get()

int // O [rcd] Return success code
tpt_vrt_get // [fnc] Virtual temperature
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *q_H2O_vpr, // I [kg kg-1] Specific humidity
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *tpt_vrt) // O [K] Virtual temperature
{
  // Purpose: Given temperature and moisture arrays, return virtual temperature
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_vrt[lon_idx]=tpt[lon_idx]*(1.0+eps_H2O_rcp_m1*q_H2O_vpr[lon_idx]); // [K] Virtual temperature
  }  // end loop over lon
  return rcd;
} // end tpt_vrt_get()

int // O [rcd] Return success code
ppr_H2O_get // [fnc] Partial pressure of H2O
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *ppr_H2O, // O [Pa] Partial pressure of H2O
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp *q_H2O_vpr) // I [kg kg-1] Specific humidity
{
  // Purpose: Given pressure and moisture arrays, return vapor pressure
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    ppr_H2O[lon_idx]=q_H2O_vpr[lon_idx]*prs[lon_idx]/(eps_H2O+one_mns_eps_H2O*q_H2O_vpr[lon_idx]); // [Pa] Partial pressure of H2O
  }  // end loop over lon
  return rcd;
} // end ppr_H2O_get()

int // O [rcd] Return success code
q_H2O_vpr_get // [fnc] Specific humidity
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *ppr_H2O, // I [Pa] Partial pressure of H2O
 const prc_cmp *prs, // I [Pa] Pressure
 prc_cmp *q_H2O_vpr) // O [kg kg-1] Specific humidity
{
  // Purpose: Given pressure and moisture arrays, return vapor mixing ratio
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    q_H2O_vpr[lon_idx]=eps_H2O*ppr_H2O[lon_idx]/(prs[lon_idx]-one_mns_eps_H2O*ppr_H2O[lon_idx]); // [kg kg] Specific humidity
  }  // end loop over lon
  return rcd;
} // end q_H2O_vpr_get()

int // O [rcd] Return success code
dsvpdt_H2O_lqd_PrK78 // [fnc] Derivative of saturation vapor pressure over planar liquid water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O_lqd) // O [Pa K-1] Derivative of saturation vapor pressure over planar liquid water
{
  /* Purpose: Compute derivative of saturation vapor pressure over planar liquid water
     Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625
     Range of validity is -50 C < T < 50 C
     d(svp)/dT is close to, but not exactly, the algebraic derivative of svp_H2O_lqd_PrK78
     Not sure why not, perhaps svp and d(svp)/dT were measured independently?
     Requires: <phys_cst.hh> for tpt_frz_pnt */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double dsvpdt_H2O_lqd_dbl; // [Pa K-1]
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
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    dsvpdt_H2O_lqd_dbl=cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls))))); // [mb K-1]
    dsvpdt_H2O_lqd[lon_idx]=dsvpdt_H2O_lqd_dbl*100.0; // [mb K-1] --> [Pa K-1]
    if(tpt_cls < -50.0 || tpt_cls > 50.0) rcd+=1; // Out of range error
  } // end loop over lon
  return rcd;
} // end dsvpdt_H2O_lqd_PrK78()

int // O [rcd] Return success code
dsvpdt_H2O_ice_PrK78
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O_ice) // O [Pa K-1] Derivative of saturation vapor pressure over planar ice water
{
  // Purpose: Compute derivative of saturation vapor pressure over planar ice water
  // Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625
  // Range of validity is -50 C < T < 50 C
  // Requires: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  double dsvpdt_H2O_ice_dbl; // [Pa K-1]
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
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    dsvpdt_H2O_ice_dbl=cff[0]+tpt_cls*(cff[1]+tpt_cls*(cff[2]+tpt_cls*(cff[3]+tpt_cls*(cff[4]+tpt_cls*(cff[5]+cff[6]*tpt_cls))))); // [mb K-1]
    dsvpdt_H2O_ice[lon_idx]=dsvpdt_H2O_ice_dbl*100.0; // [mb K-1] --> [Pa K-1]
    if(tpt_cls < -50.0 || tpt_cls > 50.0) rcd+=1; // Out of range error
  }  // end loop over lon
  return rcd;
} // end dsvpdt_H2O_ice_PrK78()

int // O [rcd] Return success code
dsvpdt_H2O_PrK78 // [fnc] Derivative of saturation vapor pressure over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvpdt_H2O) // O [Pa K-1] Derivative of saturation vapor pressure over planar water
{
  /* Purpose: Compute derivative of saturation vapor pressure over planar water
     Calls appropriate routine for liquid or ice depending on temperature
     Requires: <phys_cst.hh> for tpt_frz_pnt */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  double tpt_cls; // [C] Temperature
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    tpt_cls=tpt[lon_idx]-tpt_frz_pnt; // [K] -> [C] Temperature
    if(tpt_cls > 0.0){
      dsvpdt_H2O[lon_idx]=dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    }else{
      dsvpdt_H2O[lon_idx]=dsvpdt_H2O_ice_PrK78_fst_scl(tpt_cls); // [Pa] Saturation vapor pressure over planar water
    } // endif
  } // end loop over lon
  return rcd;
} // end dsvpdt_H2O_PrK78()

int // O [rcd] Return success code
dsvddt_H2O_PrK78 // [fnc] Derivative of saturation vapor density over planar water
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt, // I [K] Temperature
 prc_cmp *dsvddt_H2O) // O [kg m-3 K-1] Derivative of saturation vapor density over planar water
{
  /* Purpose: Compute derivative of saturation vapor density over planar water
     Wrapper for inline'd fast scalar function */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    dsvddt_H2O[lon_idx]=dsvddt_H2O_PrK78_fst_scl(tpt[lon_idx]); // [kg m-3 K-1] Derivative of saturation vapor density over planar water
  } // end loop over lon
  return rcd;
} // end dsvddt_H2O_PrK78()
