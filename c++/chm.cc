// $Id$ 

// Purpose: Chemical utilities for C++ programs

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <chm.hh> // Chemical utilities

int // O [rcd] Return success code
cnc_nbr_get // [fnc] Convert species mass mixing ratio and air density to species concentration
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnc_nbr, // O [mlc m-3] Number concentration of species
 const prc_cmp *dns, // I [kg m-3] Midlayer density
 const prc_cmp mmw, // I [kg mol-1] Mean molecular weight of species
 const prc_cmp *qmmr) // I [kg kg-1] Mass mixing ratio of species
{
  // Purpose: Given mass mixing ratio and density arrays, return species concentration
  // Requires: <phys_cst.hh>
  using phc::cst_Avagadro; // (6.022045e+23) [mlc mol-1] Avagadro's number
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    cnc_nbr[lon_idx]=qmmr[lon_idx]*dns[lon_idx]*cst_Avagadro/mmw; // [mlc m-3] Number concentration of species
  }  // end loop over lon
  return rcd;
} // end cnc_nbr_get()

int // O [rcd] Return success code
cnc_nbr_get // [fnc] Convert species mass mixing ratio and air pressure and temperature to species concentration
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnc_nbr, // O [mlc m-3] Number concentration of species
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp mmw, // I [kg mol-1] Mean molecular weight of species
 const prc_cmp *tpt, // I [K] Temperature
 const prc_cmp *qmmr) // I [kg kg-1] Mass mixing ratio of species
{
  // Purpose: Given mass mixing ratio and density arrays, return species concentration
  // Requires: <phys_cst.hh>
  using phc::cst_Avagadro; // (6.022045e+23) [mlc mol-1] Avagadro's number
  using phc::gas_cst_dry_air; // (287.05) [J kg-1 K-1] Gas constant of dry_air IrG81 p. 25, p. 245
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    cnc_nbr[lon_idx]=qmmr[lon_idx]*prs[lon_idx]*cst_Avagadro/(mmw*gas_cst_dry_air*tpt[lon_idx]); // [mlc m-3] Number concentration of species
  }  // end loop over lon
  return rcd;
} // end cnc_nbr_get()

////////////////////////////////////////////////////////////////////////

// Chm class

// Static members 
int Chm::nst_nbr=0; 
// Static member functions 
sng2sng_map Chm::opt2abb_map=Chm::opt2abb_map_mk();
sng2sng_map Chm::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;

  map_tmp.insert(sng2sng_map::value_type("HNO3","HNO3"));
  map_tmp.insert(sng2sng_map::value_type("nitric","HNO3"));
  map_tmp.insert(sng2sng_map::value_type("nitric acid","HNO3"));
  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end Chm::opt2abb_map_mk()
std::string Chm::opt2abb(const std::string opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("Chm::opt2abb()",opt_sng+" is unknown");
  return itr->second;
} // end Chm::opt2abb()

const sng2chm_sct_map Chm::chm_map=Chm::chm_map_mk();
const sng2chm_sct_map Chm::chm_map_mk(){ // Create chemical map
  using phc::mmw_HNO3; // (6.2995644e-02) [kg mol-1] Mean molecular weight of HNO3 HITRAN96
  using phc::cff_hnr_HNO3_H2O_298K; // (2.1e5) [mol ltr-1 atm-1] Henry's Law coefficient of HNO3 in liquid water at 298K SeP97 p. 341 Table 6.2
  const chm_sct chm[]={
    // {const std::string abb,const std::string dsc,const double mmw,const prc_cmp cff_hnr_H2O_298K,const prc_cmp rxn_ntp_H2O_298K}
    {"HNO3", // [sng] Standard abbreviation
     "Nitric acid", // [sng] Description
     mmw_HNO3, // [kg mol-1] Mean molecular weight
     cff_hnr_HNO3_H2O_298K, // [mol ltr-1 atm-1] Henry's Law coefficient in liquid water at 298 K
     72335.367, // [J mol-1] Reaction enthalpy of HNO3 at 298K
     0.0, // [frc] Normalized reactivity  of HNO3 for dry deposition SeP97 p. 973 (19.25)
     1.0e14 // [mol ltr-1 atm-1] Effective Henry's Law coefficient  of HNO3 in liquid water at 298 K at pH 6.5 SeP97 p. 975 Tbl. 19.2
    } // end HNO3
  }; // end chm
  int chm_nbr=sizeof(chm)/sizeof(chm_sct);
  long idx;
  sng2chm_sct_map chm_map_tmp;
  for(idx=0;idx<chm_nbr;idx++){
    /* NB: Define chemicals before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map. */
    chm_map_tmp.insert(sng2chm_sct_map::value_type(chm[idx].abb,chm[idx]));
  } // end loop over itr
  // NB: Return a value to initialize a static class member
  return chm_map_tmp;
} // end Chm::chm_map_mk()
prc_cmp Chm::abb2mmw(const std::string opt_sng){ // Abbreviation to mmw mapper
  // NB: Return a value to initialize a static class member
  sng2chm_sct_map::const_iterator chm_itr;
  chm_itr=chm_map.find(opt_sng);
  if(chm_itr == chm_map.end()) err_prn("Chm::abb2mmw()",opt_sng+" is unknown");
  return chm_itr->second.mmw;
} // end Chm::abb2mmw()
std::string Chm::abb2dsc(const std::string opt_sng){ // Abbreviation to description mapper
  // NB: Return a value to initialize a static class member
  sng2chm_sct_map::const_iterator chm_itr;
  chm_itr=chm_map.find(opt_sng);
  if(chm_itr == chm_map.end()) err_prn("Chm::abb2dsc()",opt_sng+" is unknown");
  return chm_itr->second.dsc;
} // end Chm::abb2dsc()

// Public member functions
// NB: Need to add overload constructors which allow command line specification of 
// chemical density, input file, and English description
Chm::Chm(const std::string chm_sng){ // Default constructor
  nst_nbr++; // Increment current # instantiated
  abb=chm_sng; // [sng] Chemical abbreviation
  recompute();
} // end Chm constructor
Chm::Chm(const std::string chm_sng,const prc_cmp mmw_chm){ // Constructor: specified mean molecular weight
  nst_nbr++; // Increment current # instantiated
  abb=chm_sng; // [sng] Chemical abbreviation
  if(mmw_chm > 0.0) mmw=mmw_chm; else mmw=abb2mmw(abb); // [kg mol-1] Mean molecular weight
  dsc=abb2dsc(abb); // [sng] Chemical description
} // end Chm constructor
Chm::~Chm(){ // Destructor
  nst_nbr--; // Decrement current # instantiated
} // end Chm destructor
void Chm::typ(const std::string sng){ // Set chemical type
  abb=sng; // [sng] Chemical abbreviation
  recompute();
} // end Chm::abb()
prc_cmp Chm::mmw_get()const{return mmw;} // [kg mol-1] Mean molecular weight
std::string Chm::abb_get()const{return abb;} // [sng] Chemical abbreviation
std::string Chm::dsc_get()const{return dsc;} // [sng] Chemical description

// Private member functions
void Chm::recompute(){
  mmw=abb2mmw(abb); // [kg mol-1] Mean molecular weight
  dsc=abb2dsc(abb); // [sng] Chemical description
} // end Chm::recompute()




