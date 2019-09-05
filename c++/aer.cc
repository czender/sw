// $Id$ 

// Purpose: Aerosol physics utilities for C++ programs

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <aer.hh> // Aerosol physics

// aer_cls class

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const aer_cls &aer_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for idx_rfr_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     cout << aer_obj; */
  srm_out << "Public contents of aer_obj: " << std::endl;
  srm_out << "aer_obj.abb_get() = " << aer_obj.abb_get() << std::endl;
  srm_out << "aer_obj.dns_get() = " << aer_obj.dns_get() << std::endl;
  srm_out << "aer_obj.dsc_get() = " << aer_obj.dsc_get() << std::endl;
  srm_out << "aer_obj.fl_idx_rfr_get() = " << aer_obj.fl_idx_rfr_get() << std::endl;
  srm_out << "aer_obj.mmw_get() = " << aer_obj.mmw_get() << std::endl;
  srm_out << "aer_obj.nst_nbr_get() = " << aer_obj.nst_nbr_get() << std::endl;
  srm_out << "aer_obj.spc_heat_get() = " << aer_obj.spc_heat_get() << std::endl;
  srm_out << "Private contents of aer_obj: " << std::endl;
  srm_out << "  Contents of idx_rfr: " << std::endl;
  // Rely on idx_rfr_cls::prn() to print itself
  // Not sure why this does not work, possiby mixing cout and srm_out?
  // aer_obj.idx_rfr->prn();
  srm_out << *(aer_obj.idx_rfr);
    
  return srm_out; // [srm] Reference to output stream for cascading
} // end aer_cls::operator<<()

// Friendly functions end
// Static members begin

const prc_cmp aer_cls::tpt_dfl=300.0; // [K] Default temperature
const std::string aer_cls::abb_sng_dfl="saharan_dust"; // [sng] Aerosol abbreviation, default
const std::complex<prc_cmp> aer_cls::idx_rfr_usr_dfl=sng2cpx("1.0+0.0i"); // [frc] Default refractive index, user-specified
int aer_cls::nst_nbr=0; // [nbr] Number of instantiated class members

// Static members end 
// Static member functions begin

int aer_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members

int // O [enm] Return success code
aer_cls::tst // [fnc] Self-test aer_cls class
(const std::string &abb_sng, // [fnc] Aerosol abbreviation
 const long obj_nbr) // [nbr] Number of objects to create/destroy
{
  // Purpose: Perform self test of aer_cls class
  int rcd_lcl(0); // [enm] Return success code
  std::string sbr_nm("aer_cls::tst"); // [sng] Subroutine name

  // Test for memory leaks
  std::cout << sbr_nm+"()"+" is testing for memory leaks by creating and destroying " << obj_nbr << " aer_cls objects..." << std::endl;
  long idx; // [idx] Counting index
  aer_cls *tst_obj; // [sct] Test object
  const prc_cmp wvl_ctr(0.5e-6); // [m] Wavelength at band center
  std::complex<prc_cmp> idx_rfr; // [frc] Refractive index
  for(idx=0;idx<obj_nbr;idx++){
    tst_obj=new aer_cls(abb_sng); // [sct] Test object
    idx_rfr=tst_obj->idx_rfr_get(wvl_ctr); // [frc] Refractive index
    std::cout << "idx = " << idx << ", nst_nbr = " << nst_nbr << std::endl;
    delete tst_obj; // [sct] Test object
  } // end loop over obj

  return rcd_lcl; // [enm] Return success code
} // end aer_cls::tst()

const sng2sng_map aer_cls::opt2abb_map=aer_cls::opt2abb_map_mk();
sng2sng_map aer_cls::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;

  /* NB: 20060621 Making attempt to standardize names to contain references
     
     "Volatile" indicates value corresponding to key may change
     Volatile keys are assigned to best known refractive indices
     As better data become available, key will be re-directed to new values
     Non-volatile keys are those which contain specific references
     References never change */

  map_tmp.insert(sng2sng_map::value_type("aeronet_Bhr","aeronet_Bhr"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_SdA","aeronet_SdA"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_CpV","aeronet_CpV"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_Bnz","aeronet_Bnz"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_Ogd","aeronet_Ogd"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_Brb","aeronet_Brb"));
  map_tmp.insert(sng2sng_map::value_type("aeronet_Mng","aeronet_Mng"));

  map_tmp.insert(sng2sng_map::value_type("afghan_dust","afghan_dust"));
  map_tmp.insert(sng2sng_map::value_type("afghan","afghan_dust"));
  map_tmp.insert(sng2sng_map::value_type("Afghan","afghan_dust"));
  map_tmp.insert(sng2sng_map::value_type("SAJ93","afghan_dust"));

  map_tmp.insert(sng2sng_map::value_type("air","air"));
  map_tmp.insert(sng2sng_map::value_type("Air","air"));

  map_tmp.insert(sng2sng_map::value_type("biomass","biomass"));
  map_tmp.insert(sng2sng_map::value_type("Biomass","biomass"));
  map_tmp.insert(sng2sng_map::value_type("bmb","biomass"));
  map_tmp.insert(sng2sng_map::value_type("BMB","biomass"));
  map_tmp.insert(sng2sng_map::value_type("Bio-mass","biomass"));
  map_tmp.insert(sng2sng_map::value_type("bio-mass","biomass"));

  map_tmp.insert(sng2sng_map::value_type("calcite_pellet_LQB93","calcite_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("calcite_LQB93","calcite_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("calcite","calcite_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("CaCO3","calcite_pellet_LQB93"));

  map_tmp.insert(sng2sng_map::value_type("calcite_eray_LQB93","calcite_eray_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("calcite_eray","calcite_eray_LQB93"));


  map_tmp.insert(sng2sng_map::value_type("calcite_oray_LQB93","calcite_oray_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("calcite_oray","calcite_oray_LQB93"));

  map_tmp.insert(sng2sng_map::value_type("limestone_dsp_QOL78","limestone_dsp_QOL78"));
  map_tmp.insert(sng2sng_map::value_type("limestone_dsp","limestone_dsp_QOL78"));

  map_tmp.insert(sng2sng_map::value_type("limestone_krk_QOL78","limestone_krk_QOL78"));
  map_tmp.insert(sng2sng_map::value_type("limestone_krk","limestone_krk_QOL78"));

  map_tmp.insert(sng2sng_map::value_type("dust_like","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("Dst","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("dst","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("Dust","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("dust","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("Dust-like","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("dust-like","dust_like"));
  map_tmp.insert(sng2sng_map::value_type("Dust_like","dust_like"));

  map_tmp.insert(sng2sng_map::value_type("Fe2O3_doccd","Fe2O3_doccd"));
  map_tmp.insert(sng2sng_map::value_type("Fe2O3","Fe2O3_doccd")); // Volatile
  map_tmp.insert(sng2sng_map::value_type("hematite","Fe2O3_doccd")); // Volatile
  map_tmp.insert(sng2sng_map::value_type("Hematite","Fe2O3_doccd")); // Volatile

  map_tmp.insert(sng2sng_map::value_type("Fe2O3_avg_hitran96","Fe2O3_avg_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("Fe2O3_oray_hitran96","Fe2O3_oray_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("Fe2O3_eray_hitran96","Fe2O3_eray_hitran96"));

  map_tmp.insert(sng2sng_map::value_type("Fe2O3_avg_roush","Fe2O3_avg_roush"));
  map_tmp.insert(sng2sng_map::value_type("Fe2O3_oray_roush","Fe2O3_oray_roush"));
  map_tmp.insert(sng2sng_map::value_type("Fe2O3_eray_roush","Fe2O3_eray_roush"));

  map_tmp.insert(sng2sng_map::value_type("gypsum","gypsum_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("CaSO4_2H2O","gypsum_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("gypsum_pellet_LQB93","gypsum_pellet_LQB93"));
  map_tmp.insert(sng2sng_map::value_type("gypsum_LQB93","gypsum_pellet_LQB93"));

  map_tmp.insert(sng2sng_map::value_type("gypsum_EX_LQB93","gypsum_EX_LQB93"));

  map_tmp.insert(sng2sng_map::value_type("h2o_lqd","h2o_lqd"));
  map_tmp.insert(sng2sng_map::value_type("H2O liquid","h2o_lqd"));
  map_tmp.insert(sng2sng_map::value_type("liquid","h2o_lqd"));
  map_tmp.insert(sng2sng_map::value_type("liq","h2o_lqd"));

  map_tmp.insert(sng2sng_map::value_type("h2o_ice","h2o_ice"));
  map_tmp.insert(sng2sng_map::value_type("H2O ice","h2o_ice"));
  map_tmp.insert(sng2sng_map::value_type("ice","h2o_ice"));
  map_tmp.insert(sng2sng_map::value_type("ice","h2o_ice"));

  map_tmp.insert(sng2sng_map::value_type("H2SO4_300K_PaW75","H2SO4_300K_PaW75"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_300K_PaW75","H2SO4_300K_PaW75"));
  map_tmp.insert(sng2sng_map::value_type("PaW75","H2SO4_300K_PaW75"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_PaW75","H2SO4_300K_PaW75"));

  map_tmp.insert(sng2sng_map::value_type("H2SO4_210K_61_NNM98","H2SO4_210K_61_NNM98"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_210K_61_NNM98","H2SO4_210K_61_NNM98"));

  map_tmp.insert(sng2sng_map::value_type("H2SO4_220K_72_NNM98","H2SO4_220K_72_NNM98"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_220K_72_NNM98","H2SO4_220K_72_NNM98"));
  map_tmp.insert(sng2sng_map::value_type("NNM98","H2SO4_220K_72_NNM98"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_NNM98","H2SO4_220K_72_NNM98"));

  map_tmp.insert(sng2sng_map::value_type("h2so4_300K","h2so4_300K"));
  map_tmp.insert(sng2sng_map::value_type("H2SO4","h2so4_300K"));
  map_tmp.insert(sng2sng_map::value_type("h2so4","h2so4_300K"));
  map_tmp.insert(sng2sng_map::value_type("H2SO4 at 300K","h2so4_300K"));

  map_tmp.insert(sng2sng_map::value_type("H2SO4 at 215K","h2so4_215K"));
  map_tmp.insert(sng2sng_map::value_type("h2so4_215K","h2so4_215K"));

  map_tmp.insert(sng2sng_map::value_type("illite","illite"));
  map_tmp.insert(sng2sng_map::value_type("Illite","illite"));

  map_tmp.insert(sng2sng_map::value_type("kaolinite","kaolinite"));
  map_tmp.insert(sng2sng_map::value_type("Kaolinite","kaolinite"));

  map_tmp.insert(sng2sng_map::value_type("soot","lac_ChC90"));
  map_tmp.insert(sng2sng_map::value_type("Soot","lac_ChC90"));
  map_tmp.insert(sng2sng_map::value_type("Black Carbon","lac_ChC90"));
  map_tmp.insert(sng2sng_map::value_type("bc","lac_ChC90"));
  map_tmp.insert(sng2sng_map::value_type("carbon","lac_ChC90"));

  map_tmp.insert(sng2sng_map::value_type("lac_ChC90","lac_ChC90"));
  map_tmp.insert(sng2sng_map::value_type("lac_DKS91","lac_DKS91"));
  map_tmp.insert(sng2sng_map::value_type("lac_HKS98","lac_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("lac_WaW80","lac_WaW80"));
  map_tmp.insert(sng2sng_map::value_type("lac_YZS00","lac_YZS00"));

  map_tmp.insert(sng2sng_map::value_type("meteoric_dust","meteoric_dust"));
  map_tmp.insert(sng2sng_map::value_type("Meteoric dust","meteoric_dust"));
  map_tmp.insert(sng2sng_map::value_type("meteoric","meteoric_dust"));
  map_tmp.insert(sng2sng_map::value_type("meteor","meteoric_dust"));
  map_tmp.insert(sng2sng_map::value_type("Meteor","meteoric_dust"));

  map_tmp.insert(sng2sng_map::value_type("MgSO4","MgSO4"));

  map_tmp.insert(sng2sng_map::value_type("mineral_dust","mineral_dust"));
  map_tmp.insert(sng2sng_map::value_type("Mineral dust","mineral_dust"));
  map_tmp.insert(sng2sng_map::value_type("mineral dust","mineral_dust"));
  map_tmp.insert(sng2sng_map::value_type("mnr_dst","mineral_dust"));
  map_tmp.insert(sng2sng_map::value_type("dst_mnr","mineral_dust"));

  map_tmp.insert(sng2sng_map::value_type("montmorillonite","montmorillonite"));
  map_tmp.insert(sng2sng_map::value_type("Montmorillonite","montmorillonite"));

  map_tmp.insert(sng2sng_map::value_type("NaCl","NaCl"));
  map_tmp.insert(sng2sng_map::value_type("NaCl_Fla04","NaCl_Fla04"));

  map_tmp.insert(sng2sng_map::value_type("oceanic","oceanic"));
  map_tmp.insert(sng2sng_map::value_type("Oceanic","oceanic"));
  map_tmp.insert(sng2sng_map::value_type("ocean","oceanic"));
  map_tmp.insert(sng2sng_map::value_type("Ocean","oceanic"));

  map_tmp.insert(sng2sng_map::value_type("saharan_dust","saharan_dust"));
  map_tmp.insert(sng2sng_map::value_type("saharan","saharan_dust"));
  map_tmp.insert(sng2sng_map::value_type("Saharan","saharan_dust"));
  map_tmp.insert(sng2sng_map::value_type("Vol73","saharan_dust"));
  map_tmp.insert(sng2sng_map::value_type("vol73","saharan_dust"));

  map_tmp.insert(sng2sng_map::value_type("sea_salt_Fla04","sea_salt_Fla04"));
  map_tmp.insert(sng2sng_map::value_type("seasalt_Fla04","sea_salt_Fla04"));

  map_tmp.insert(sng2sng_map::value_type("sea_salt","sea_salt_GADS"));
  map_tmp.insert(sng2sng_map::value_type("sea","sea_salt_GADS"));
  map_tmp.insert(sng2sng_map::value_type("sea_salt","sea_salt_GADS"));
  map_tmp.insert(sng2sng_map::value_type("sea_salt_GADS","sea_salt_GADS"));

  map_tmp.insert(sng2sng_map::value_type("SiO2","SiO2_avg_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("quartz","SiO2_avg_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("Quartz","SiO2_avg_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("SiO2_avg_hitran96","SiO2_avg_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("SiO2_oray_hitran96","SiO2_oray_hitran96"));
  map_tmp.insert(sng2sng_map::value_type("SiO2_eray_hitran96","SiO2_eray_hitran96"));

  map_tmp.insert(sng2sng_map::value_type("SiO2_avg_roush","SiO2_avg_roush"));
  map_tmp.insert(sng2sng_map::value_type("SiO2_oray_roush","SiO2_oray_roush"));
  map_tmp.insert(sng2sng_map::value_type("SiO2_eray_roush","SiO2_eray_roush"));

  map_tmp.insert(sng2sng_map::value_type("Sulfate","sulfate"));
  map_tmp.insert(sng2sng_map::value_type("sulfate","sulfate"));

  map_tmp.insert(sng2sng_map::value_type("toms_dust_STD03","toms_dust_STD03"));

  map_tmp.insert(sng2sng_map::value_type("volcanic_dust","volcanic_dust"));
  map_tmp.insert(sng2sng_map::value_type("Volcanic dust","volcanic_dust"));
  map_tmp.insert(sng2sng_map::value_type("volcanic dust","volcanic_dust"));
  map_tmp.insert(sng2sng_map::value_type("volcano","volcanic_dust"));
  map_tmp.insert(sng2sng_map::value_type("dst_vlc","volcanic_dust"));

  map_tmp.insert(sng2sng_map::value_type("wsoc80_HKS98","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("oc","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("OC","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("organic_carbon","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("waso","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("waso80","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc80_HKS98","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc80","wsoc80_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc50_HKS98","wsoc50_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc50","wsoc50_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc70_HKS98","wsoc70_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc70","wsoc70_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc90_HKS98","wsoc90_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc90","wsoc90_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc95_HKS98","wsoc95_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc95","wsoc95_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc98_HKS98","wsoc98_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc98","wsoc98_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc99_HKS98","wsoc99_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc99","wsoc99_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc00_HKS98","wsoc00_HKS98"));
  map_tmp.insert(sng2sng_map::value_type("wsoc00","wsoc00_HKS98"));

  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end aer_cls::opt2abb_map_mk()
std::string aer_cls::opt2abb(const std::string &opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("aer_cls::opt2abb",opt_sng+" is unknown");
  return itr->second;
} // end aer_cls::opt2abb()

const sng2mnr_sct_map aer_cls::mnr_map=aer_cls::mnr_map_mk();
sng2mnr_sct_map aer_cls::mnr_map_mk(){ // Create abbreviation map
  // Purpose: Initialize database of physical/optical/mineralogical properties

  using phc::dns_CaCO3; // (2710.0) [kg m-3] Density of calcite (calcium carbonate) (http://webmineral.com/data/Calcite.shtml)
  using phc::dns_CaSO4_2H2O; // (2300.0) [kg m-3] Density of gypsum (hydrous calcium sulfate) (http://webmineral.com/data/Gypsum.shtml)
  using phc::dns_Fe2O3; // (5260.0) [kg m-3] Density of Fe2O3
  using phc::dns_H2O_ice_std; // (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
  using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
  using phc::dns_H2SO4_61pct_273K_std; // (1526.16) [kg m-3] Density of 61% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
  using phc::dns_H2SO4_72pct_273K_std; // (1652.95) [kg m-3] Density of 72% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
  using phc::dns_H2SO4_75pct_273K_std; // (1688.83) [kg m-3] Density of 75% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
  using phc::dns_H2SO4_75pct_303K_std; // (1659.65) [kg m-3] Density of 75% H2SO4 solution (by weight) at 303 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
  using phc::dns_H2SO4_std; // (1851.69) [kg m-3] Density of 100% H2SO4 at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
  using phc::dns_MgSO4_std; // (2660.0) [kg m-3] Density of MgSO4 Yu's Thesis Ch. 4
  using phc::dns_NaCl_std; // (2170.0) [kg m-3] Density of NaCl http://www.crystran.co.uk/nacldata.htm
  using phc::dns_SiO2; // (2620.0) [kg m-3] Density of SiO2
  using phc::dns_dst_DKS91; // (1600.0) [kg m-3] Density of dust DKS91 p. 118
  using phc::dns_dst_PaG77; // (2500.0) [kg m-3] Density of dust PaG77 p. 2076
  using phc::dns_dst_Vol73; // (2650.0) [kg m-3] Density of dust Vol73
  using phc::dns_dst_mtr_DKS91; // (2500.0) [kg m-3] Density of meteoric dust DKS91 p. 118
  using phc::dns_dst_std; // (2500.0) [kg m-3] Standard density of dust
  using phc::dns_illite; // (2750.0) [kg m-3] Density of illite http://webmineral.com/data/Illite.shtml
  using phc::dns_kaolinite; // (2600.0) [kg m-3] Density of kaolinite http://webmineral.com/data/Kaolinite.shtml
  using phc::dns_montmorillonite; // (2350.0) [kg m-3] Density of montmorillonite http://webmineral.com/data/Montmorillonite.shtml
  using phc::dns_lac_DKS91; // (2300.0) [kg m-3] Standard density of soot DKS91 p. 118
  using phc::dns_lac_HKS98; // (1000.0) [kg m-3] Density of soot HKS98 p. 836
  using phc::dns_lac_WaW80; // (2050.0) [kg m-3] Density of soot WaW80 p. 
  using phc::dns_lac_YZS00; // (1860.0) [kg m-3] Density of sooty sulfate YZS00
  using phc::dns_lac_std; // (1800.0) [kg m-3] Standard density of soot
  using phc::dns_wsoc_HKS98; // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
  using phc::mmw_C; // (12.011e-03) [kg mol-1] Mean molecular weight of C IUPAC
  using phc::mmw_CaCO3; // (100.087e-03) [kg mol-1] Mean molecular weight of CaCO3
  using phc::mmw_CaSO4_2H2O; // (172.17e-03) [kg mol-1] Mean molecular weight of CaSO4_2H2O
  using phc::mmw_Fe2O3; // (159.692e-03) [kg mol-1] Mean molecular weight of Fe2O3
  using phc::mmw_H2O; // (1.8015259e-02) [kg mol-1] Mean molecular weight of H2O HITRAN96
  using phc::mmw_H2SO4; // [kg mol-1] Mean molecular weight of H2SO4
  using phc::mmw_MgSO4; // [kg mol-1] Mean molecular weight of MgSO4
  using phc::mmw_NH4NH4SO4; // (132.141e-03) [kg mol-1] Mean molecular weight of NH4NH4SO4
  using phc::mmw_NaCl; // (58.4425e-03) [kg mol-1] Mean molecular weight of NaCl
  using phc::mmw_SiO2; // (60.0843e-03) [kg mol-1] Mean molecular weight of SiO2
  using phc::mmw_dry_air; // (28.9644e-3) [kg mol-1] (Source: radcsw.F in CCM2/3)
  using phc::mmw_illite; // (389.34e-03) [kg mol-1] Mean molecular weight of Illite http://webmineral.com/data/Illite.shtml
  using phc::mmw_kaolinite; // (258.16e-03) [kg mol-1] Mean molecular weight of Kaolinite http://webmineral.com/data/Kaolinite.shtml
  using phc::mmw_montmorillonite; // (549.07e-03) [kg mol-1] Mean molecular weight of Montmorillonite http://webmineral.com/data/Montmorillonite.shtml
  using phc::spc_heat_CaCO3_sld; // () [J kg-1 K-1] Specific heat capacity of calcite (calcium carbonate) CRC95 p. 5-24
  using phc::spc_heat_CaSO4_sld; // () [J kg-1 K-1] Specific heat capacity of gypsum (hydrous calcium sulfate) CRC95 p. 5-9
  using phc::spc_heat_Fe2O3_sld; // () [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
  using phc::spc_heat_H2O_ice; // (2108.0) [J kg-1 K-1] Specific heat capacity of ice water (quora.com)
  using phc::spc_heat_H2O_lqd; // (4187.0) [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15 (fxm: better to derive this from density and volumetric spec. heat?)
  using phc::spc_heat_SiO2_sld; // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz CRC95 p. 5-21
  using phc::spc_heat_dry_air; // (1004.697) [J kg-1 K-1] IrG81 p. 25
  
  // fxm: Specific heat is either for water or for quartz
  // fxm: Distinguish functions from files as in spc_slr structure
  const mnr_sct mnr[]={
    {"H2SO4_210K_61_NNM98", // [sng] Particle abbreviation
     "61% H2SO4 at 210K (NNM98)", // [sng] Particle description
     "H2SO4_H2O", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_NNM98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_61pct_273K_std}, // [kg m-3] Density of 61% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
    {"H2SO4_220K_72_NNM98", // [sng] Particle abbreviation
     "72% H2SO4 at 220K (NNM98)", // [sng] Particle description
     "H2SO4-H2O", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_NNM98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_72pct_273K_std}, // [kg m-3] Density of 72% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
    {"H2SO4_300K_PaW75", // [sng] Particle abbreviation
     "75% H2SO4 at 300K (PaW75)", // [sng] Particle description
     "H2SO4-H2O", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_PaW75.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_75pct_303K_std}, // [kg m-3] Density of 75% H2SO4 solution (by weight) at 303 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
    {"aeronet_Bhr", // [sng] Particle abbreviation
     "Aeronet Bahrain", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_Bnz", // [sng] Particle abbreviation
     "Aeronet Banizoumbou", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_Brb", // [sng] Particle abbreviation
     "Aeronet Barbados", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_CpV", // [sng] Particle abbreviation
     "Aeronet Cape Verde", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_Mng", // [sng] Particle abbreviation
     "Aeronet Mongolia", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_Ogd", // [sng] Particle abbreviation
     "Aeronet Ougadougou", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"aeronet_SdA", // [sng] Particle abbreviation
     "Aeronet Saudi Arabia", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DHE02.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"afghan_dust", // [sng] Particle abbreviation
     "Afghan dust", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_SAJ93.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"air", // [sng] Particle abbreviation
     "Dry air", // [sng] Particle description
     "N2O2", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_air_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_dry_air, // [kg mol-1] Mean molecular weight
     spc_heat_dry_air, // [J kg-1 K-1] Specific heat capacity of dry air
     1.0}, // [kg m-3] Bulk density
    {"biomass", // [sng] Particle abbreviation
     "Biomass", // [sng] Particle description
     "C", // [sng] Molecular composition)
     fio::data_file_path_get("idx_rfr_DKS91.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     1000.0}, // [kg m-3] Bulk density Biomass fxm: get BMB density
    {"calcite_pellet_LQB93", // [sng] Particle abbreviation
     "Calcite (calcium carbonate)", // [sng] Particle description
     "CaCO3", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaCO3, // [kg mol-1] Mean molecular weight
     spc_heat_CaCO3_sld, // [J kg-1 K-1] Specific heat capacity of calcite
     dns_CaCO3}, // (2710.0) [kg m-3] Density of calcite
    {"calcite_oray_LQB93", // [sng] Particle abbreviation
     "Calcite (calcium carbonate)", // [sng] Particle description
     "CaCO3", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaCO3, // [kg mol-1] Mean molecular weight
     spc_heat_CaCO3_sld, // [J kg-1 K-1] Specific heat capacity of calcite
     dns_CaCO3}, // (2710.0) [kg m-3] Density of calcite
    {"calcite_eray_LQB93", // [sng] Particle abbreviation
     "Calcite (calcium carbonate)", // [sng] Particle description
     "CaCO3", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaCO3, // [kg mol-1] Mean molecular weight
     spc_heat_CaCO3_sld, // [J kg-1 K-1] Specific heat capacity of calcite
     dns_CaCO3}, // (2710.0) [kg m-3] Density of calcite
    {"gypsum_EX_LQB93", // [sng] Particle abbreviation
     "Gypsum (hydrous calcium sulfate)", // [sng] Particle description
     "CaSO4_2H2O", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaSO4_2H2O, // [kg mol-1] Mean molecular weight
     spc_heat_CaSO4_sld, // [J kg-1 K-1] Specific heat capacity of gypsum
     dns_CaSO4_2H2O}, // (2300.0) [kg m-3] Density of gypsum
    {"gypsum_pellet_LQB93", // [sng] Particle abbreviation
     "Gypsum (hydrous calcium sulfate)", // [sng] Particle description
     "CaSO4_2H2O", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaSO4_2H2O, // [kg mol-1] Mean molecular weight
     spc_heat_CaSO4_sld, // [J kg-1 K-1] Specific heat capacity of gypsum
     dns_CaSO4_2H2O}, // (2300.0) [kg m-3] Density of gypsum
    {"dust_like", // [sng] Particle abbreviation
     "Dust-like", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DKS91.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_DKS91}, // (1600.0) [kg m-3] Density of dust DKS91 p. 118
    {"Fe2O3_avg_hitran96", // [sng] Particle abbreviation
     "Hematite O/E-ray mean (Shettle/HITRAN96)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_oray_hitran96", // [sng] Particle abbreviation
     "Hematite O-ray (Shettle/HITRAN96)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_eray_hitran96", // [sng] Particle abbreviation
     "Hematite E-ray (Shettle/HITRAN96)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_avg_roush", // [sng] Particle abbreviation
     "Hematite O/E-ray mean (Querry/Roush)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Fe2O3.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_oray_roush", // [sng] Particle abbreviation
     "Hematite O-ray (Querry/Roush)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Fe2O3.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_eray_roush", // [sng] Particle abbreviation
     "Hematite E-ray (Querry/Roush)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Fe2O3.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"Fe2O3_doccd", // [sng] Particle abbreviation
     "Hematite (DOCCD)", // [sng] Particle description
     "Fe2O3", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_triaud_Fe2O3.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_Fe2O3, // (159.692e-03) [kg mol-1] Mean molecular weight
     spc_heat_Fe2O3_sld, // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
     dns_Fe2O3}, // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
    {"h2o_War84", // [sng] Particle abbreviation
     "H2O ice (War84)", // [sng] Particle description
     "H2O", // [sng] Molecular composition
     "Function ~/idx_rfr_H2O:idx_rfr_H2O_ice_get_War84()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2O, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_ice, // [J kg-1 K-1] Specific heat capacity of ice water 
     dns_H2O_ice_std}, // (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
    {"h2o_ice", // [sng] Particle abbreviation
     "H2O ice", // [sng] Particle description
     "H2O", // [sng] Molecular composition
     "Function ~/idx_rfr_H2O:idx_rfr_H2O_ice_get_WaB08()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2O, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_ice, // [J kg-1 K-1] Specific heat capacity of ice water 
     dns_H2O_ice_std}, // (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
    {"h2o_lqd", // [sng] Particle abbreviation
     "H2O liquid", // [sng] Particle description
     "H2O", // [sng] Molecular composition
     "Function ~/idx_rfr_H2O/idx_rfr_H2O_lqd_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2O, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2O_lqd_std}, // (1000.0) [kg m-3] Density of liquid water
    {"h2so4_215K", // [sng] Particle abbreviation
     "75% H2SO4 at 215K (HSL88)", // [sng] Particle description
     "H2SO4-H2O", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_75pct_273K_std}, // [kg m-3] Density of 75% H2SO4 solution (by weight) at 273 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
    {"h2so4_300K", // [sng] Particle abbreviation
     "75% H2SO4 at 300K (HSL88)", // [sng] Particle description
     "H2SO4-H2O", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_75pct_303K_std}, // [kg m-3] Density of 75% H2SO4 solution (by weight) at 303 K (Timmermans in ${HOME}/idx_rfr/hitran/timmerma.dat)
    {"illite", // [sng] Particle abbreviation
     "Illite (EgH79, Que87)", // [sng] Particle description
     "(K,H3O)(Al,Mg,Fe)2(Si,Al)4O10[(OH)2,(H2O)]", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Ill_Kao_Mnt.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_illite, // (389.34e-03) [kg mol-1] Mean molecular weight of Illite http://webmineral.com/data/Illite.shtml
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_illite}, // (2750.0) [kg m-3] Density of illite
    {"kaolinite", // [sng] Particle abbreviation
     "Kaolinite (EgH79, RPO91)", // [sng] Particle description
     "Al2Si2O5(OH)4", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Ill_Kao_Mnt.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_kaolinite, // (258.16e-03) [kg mol-1] Mean molecular weight of Kaolinite http://webmineral.com/data/Kaolinite.shtml
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_kaolinite}, // (2600.0) [kg m-3] Density of kaolinite
    {"lac_ChC90", // [sng] Particle abbreviation
     "Soot from Chang and Charalampopoulos (1990)", // [sng] Particle description
     "C", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_ChC90_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_lac_std}, // (1800.0) [kg m-3] Standard density of soot
    {"lac_DKS91", // [sng] Particle abbreviation
     "Soot", // [sng] Particle description
     "C", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_DKS91_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_lac_DKS91}, // (2300.0) [kg m-3] Density of soot DKS91 p. 118
    {"lac_HKS98", // [sng] Particle abbreviation
     "Soot (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_lac_HKS98}, // (1000.0) [kg m-3] Density of soot HKS98 p. 836
    {"lac_WaW80", // [sng] Particle abbreviation
     "Soot from Warren & Wiscombe (1980)", // [sng] Particle description
     "C", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_WaW80_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_lac_WaW80}, // (2050.0) [kg m-3] Bulk density soot 
    {"lac_YZS00", // [sng] Particle abbreviation
     "Sooty sulfate (YZS00)", // [sng] Particle description
     "H2SO4C", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_YZS00_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_C, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_lac_YZS00}, // (1860.0) [kg m-3] Bulk density Sooty sulfate of YZS00
    {"limestone_dsp_QOL78", // [sng] Particle abbreviation
     "Limestone (QOL78 DA)", // [sng] Particle description
     "CaCO3 (amorphous calcite = calcium carbonate)", // [sng] Molecular composition
     "Function ~/c++/idx_rfr.cc:idx_rfr_LQB93_get()", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaCO3, // [kg mol-1] Mean molecular weight
     spc_heat_CaCO3_sld, // [J kg-1 K-1] Specific heat capacity of calcite
     dns_CaCO3}, // (2710.0) [kg m-3] Density of calcite
    {"limestone_krk_QOL78", // [sng] Particle abbreviation
     "Limestone (QOL78 K-K)", // [sng] Particle description
     "CaCO3 (amorphous calcite = calcium carbonate)", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_limestone.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_CaCO3, // [kg mol-1] Mean molecular weight
     spc_heat_CaCO3_sld, // [J kg-1 K-1] Specific heat capacity of calcite
     dns_CaCO3}, // (2710.0) [kg m-3] Density of calcite
    {"montmorillonite", // [sng] Particle abbreviation
     "Montmorillonite (EgH79, RPO91)", // [sng] Particle description
     "(Na,Ca)0,3(Al,Mg)2Si4O10(OH)2.n(H2O)", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Ill_Kao_Mnt.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_montmorillonite, // (549.07e-03) [kg mol-1] Mean molecular weight of Montmorillonite http://webmineral.com/data/Montmorillonite.shtml
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_montmorillonite}, // (2350.0) [kg m-3] Density of montmorillonite
    {"meteoric_dust", // [sng] Particle abbreviation
     "Meteoric dust (Shettle/HITRAN96)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_mtr_DKS91}, // [kg m-3] Bulk density Meteoric dust DKS91 p. 118
    {"mineral_dust", // [sng] Particle abbreviation
     "Mineral dust (DKS91)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DKS91.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"NaCl", // [sng] Particle abbreviation
     "Sodium chloride", // [sng] Particle description
     "NaCl", // [sng] Molecular composition
     "", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NaCl, // [kg mol-1] Mean molecular weight % fxm 
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water
     dns_NaCl_std}, // [kg m-3] Bulk density % http://www.crystran.co.uk/nacldata.htm 
    {"NaCl_Fla04", // [sng] Particle abbreviation
     "Sodium chloride", // [sng] Particle description
     "NaCl", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Fla04.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NaCl, // [kg mol-1] Mean molecular weight % fxm 
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water
     dns_NaCl_std}, // [kg m-3] Bulk density % http://www.crystran.co.uk/nacldata.htm 
    {"MgSO4", // [sng] Particle abbreviation
     "Magnesium sulfate", // [sng] Particle description
     "MgSO4", // [sng] Molecular composition
     "", // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_MgSO4, // [kg mol-1] Mean molecular weight % fxm 
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water
     dns_MgSO4_std}, // [kg m-3] Bulk density from Yu's thesis
    {"oceanic", // [sng] Particle abbreviation
     "Oceanic (DKS91)", // [sng] Particle description
     "NaCl MgSO4", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DKS91.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NaCl, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     1000.0}, // [kg m-3] Bulk density % fxm
    {"saharan_dust", // [sng] Particle abbreviation
     "Saharan dust (Vol73)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Vol73.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_Vol73}, // (2650.0) [kg m-3] Density of dust Vol73
    {"sea_salt_GADS", // [sng] Particle abbreviation
     "Sea salt", // [sng] Particle description
     "NaCl", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_GADS.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NaCl, // [kg mol-1] Mean molecular weight % fxm 
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_NaCl_std}, // [kg m-3] Bulk density % http://www.crystran.co.uk/nacldata.htm 
    {"sea_salt_Fla04", // [sng] Particle abbreviation
     "Sea salt", // [sng] Particle description
     "NaCl", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_Fla04.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NaCl, // [kg mol-1] Mean molecular weight % fxm 
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_NaCl_std}, // [kg m-3] Bulk density % http://www.crystran.co.uk/nacldata.htm 
    {"SiO2_avg_hitran96", // [sng] Particle abbreviation
     "Quartz O/E-ray mean (Shettle/HITRAN96)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"SiO2_oray_hitran96", // [sng] Particle abbreviation
     "Quartz O-ray (Shettle/HITRAN96)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"SiO2_eray_hitran96", // [sng] Particle abbreviation
     "Quartz E-ray (Shettle/HITRAN96)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle_Fe2O3_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"SiO2_avg_roush", // [sng] Particle abbreviation
     "Quartz O/E-ray mean (Querry/Roush)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"SiO2_oray_roush", // [sng] Particle abbreviation
     "Quartz O-ray (Querry/Roush)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"SiO2_eray_roush", // [sng] Particle abbreviation
     "Quartz E-ray (Querry/Roush)", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_SiO2.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // (60.0843e-03) [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz
     dns_SiO2}, // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
    {"sulfate", // [sng] Particle abbreviation
     "Sulfate", // [sng] Particle description
     "H2SO4", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_DKS91.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_H2SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_H2SO4_std}, // [kg m-3] Bulk density Sulfate (same density as H2SO4 for now)
    {"toms_dust_STD03", // [sng] Particle abbreviation
     "TOMS dust retrievals in Tor02/STD03/Tor03", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_STD03.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     dns_dst_std}, // (2500.0) [kg m-3] Standard density of dust
    {"volcanic_dust", // [sng] Particle abbreviation
     "Volcanic dust", // [sng] Particle description
     "SiO2", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_shettle.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_SiO2, // [kg mol-1] Mean molecular weight
     spc_heat_SiO2_sld, // [J kg-1 K-1] Specific heat capacity of quartz
     2760.0}, // [kg m-3] Bulk density Volcanic dust DKS91 p. 118 (Volz 1973)
    {"wsoc00_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.00 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc50_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.50 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc70_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.70 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc80_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.80 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc90_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.90 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc95_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.95 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc98_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.98 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc99_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.99 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98}, // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
    {"wsoc00_HKS98", // [sng] Particle abbreviation
     "Water soluble organic/inorganic mixture at RH=0.0 (HKS98)", // [sng] Particle description
     "C", // [sng] Molecular composition
     fio::data_file_path_get("idx_rfr_HKS98.nc"), // [sng] File containing refractive indices
     (idx_rfr_fnc_ptr_typ)CEWI_NULL, // [fnc] Function to compute refractive indices
     mmw_NH4NH4SO4, // [kg mol-1] Mean molecular weight
     spc_heat_H2O_lqd, // [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15
     dns_wsoc_HKS98} // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
  }; // end mnr_sct mnr[]
  long idx; // [idx] Counting index
  int mnr_nbr=sizeof(mnr)/sizeof(mnr_sct); // [nbr] Number mineral types
  sng2mnr_sct_map mnr_map_tmp; // [sct] Map with key=abbreviation, value=mineral structure
  for(idx=0;idx<mnr_nbr;idx++){
    /* fxm: Define variables before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map. */
    mnr_map_tmp.insert(sng2mnr_sct_map::value_type(mnr[idx].abb,mnr[idx])); // [sct] Map with key=abbreviation, value=mineral structure
  } // end loop over itr
  // Return a value to initialize a static class member
  return mnr_map_tmp;
} // end aer_cls::mnr_map_mk()

prc_cmp // [kg m-3] Bulk density
aer_cls::abb2dns(const std::string &abb_sng){ // [fnc] Abbreviation to density mapper
  // Purpose: Return density associated with mineral abbreviation string
  return mnr_map.find(abb_sng)->second.dns; // [kg m-3] Bulk density
} // end aer_cls::abb2dns()

prc_cmp // [kg mol-1] Mean molecular weight
aer_cls::abb2mmw(const std::string &abb_sng){ // [fnc] Abbreviation to molecular weight mapper
  // Purpose: Return molecular weight associated with mineral abbreviation string
  return mnr_map.find(abb_sng)->second.mmw; // [kg mol-1] Mean molecular weight
} // end aer_cls::abb2mmw()

prc_cmp // [J kg-1 K-1] Specific heat capacity
aer_cls::abb2spc_heat(const std::string &abb_sng){ // [fnc] Abbreviation to specific heat mapper
  // Purpose: Return specific heat associated with mineral abbreviation string
  return mnr_map.find(abb_sng)->second.spc_heat; // [J kg-1 K-1] Specific heat capacity
} // end aer_cls::abb2spc_heat()

std::string // [sng] Description string
aer_cls::abb2dsc(const std::string &abb_sng){ // [fnc] Abbreviation to description mapper
  // Purpose: Return description of mineral
  return mnr_map.find(abb_sng)->second.dsc; // [sng] Description string
} // end aer_cls::abb2dsc()

idx_rfr_fnc_ptr_typ // [fnc] Pointer to function returning refractive indices
aer_cls::abb2fnc // [fnc] Abbreviation to function mapper
(const std::string &abb_sng){ // [sng] Abbreviation string
  // Purpose: Return pointer to function which computes refractive indices
  return mnr_map.find(abb_sng)->second.idx_rfr_fnc_ptr; // [fnc] Pointer to function returning refractive indices
} // end aer_cls::abb2fnc()

std::string // [sng] File containing refractive indices
aer_cls::abb2fl_idx_rfr(const std::string &abb_sng){ // [fnc] Abbreviation to file mapper
  // Purpose: Return file containing mineral refractive indices
  return mnr_map.find(abb_sng)->second.fl_idx_rfr; // [sng] File containing refractive indices
} // end aer_cls::abb2fl_idx_rfr()

// Static member functions end
// Public member functions begin

aer_cls::aer_cls // [fnc] Constructor
(const std::string &abb_sng_arg, // [sng] Aerosol abbreviation
 const prc_cmp &dns_prt_arg, // [kg m-3] Density of particle
 const prc_cmp &tpt_arg, // [K] Temperature
 const std::string &fl_idx_rfr_arg, // [sng] Aerosol input file
 const bool idx_rfr_usr_flg_arg, // [flg] Refractive index is user-specified
 const std::complex<prc_cmp> idx_rfr_usr_arg) // [frc] Refractive index, user-specified
{
  // Purpose: Constructor for aerosol class

  rcd=0; // [enm] Return success code

  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls

  rcd+=abb_set(abb_sng_arg); // [sng] Aerosol abbreviation
  rcd+=dns_set(dns_prt_arg); // [kg m-3] Density of particle

  idx_rfr=new idx_rfr_cls(abb,fl_idx_rfr_arg,tpt_arg,idx_rfr_usr_flg_arg,idx_rfr_usr_arg); // [ptr] Refractive index object

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  nst_nbr++; // [nbr] Number of instantiated class members
} // end aer_cls constructor

aer_cls::~aer_cls(){ // Destructor
  // fxm: does delete idx_rfr call destructor automatically?
  delete idx_rfr; // [ptr] Refractive index object
  nst_nbr--; // [nbr] Number of instantiated class members
} // end aer_cls destructor

int // [enm] Return success code
aer_cls::abb_set(const std::string &abb_sng_arg){ // [fnc] Set aerosol type
  // Purpose: Set particle type
  abb= (abb_sng_arg == "") ? abb_sng_dfl : abb_sng_arg; // [sng] Aerosol abbreviation
  /* If constructor is called with abbreviation string next line is redundant
     However, this allows constructor to be called with option string, as well */
  abb=opt2abb(abb); // [sng] Aerosol abbreviation
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end aer_cls::abb_set()

int // [enm] Return success code
aer_cls::dns_set(const prc_cmp &dns_prt_arg){ // [fnc] Set aerosol density
  // Purpose: Set particle density
  std::string sbr_nm("aer_cls::dns_set"); // [sng] Subroutine name
  prc_cmp dns_prt_max(1.0e6); // [kg m-3] Maximum allowed particle density

  if(dns_prt_arg < 0.0 || dns_prt_arg > dns_prt_max) err_prn(sbr_nm,"dns_prt_arg is out of range, proceeding anyway");
  dns= (dns_prt_arg > 0.0 && dns_prt_arg < dns_prt_max) ? dns_prt_arg : abb2dns(abb); // [kg m-3] Density of particle
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end aer_cls::dns_set()

int // [enm] Return success code
aer_cls::fl_idx_rfr_set(const std::string &fl_idx_rfr_arg){ // [fnc] File containing refractive indices
  // Purpose: Set file containing refractive indices
  rcd+=idx_rfr->fl_idx_rfr_set(fl_idx_rfr_arg);
  return rcd; // [enm] Return success code
} // end aer_cls::fl_idx_rfr_set()

prc_cmp aer_cls::dns_get()const{return dns;} // [kg m-3] Density of particle
prc_cmp aer_cls::mmw_get()const{return mmw;} // [kg mol-1] Molecular weight
prc_cmp aer_cls::spc_heat_get()const{return spc_heat;} // [J kg-1 K-1] Specific heat capacity
std::complex<prc_cmp> aer_cls::idx_rfr_usr_get()const{return idx_rfr->idx_rfr_usr_get();} // [frc] Refractive index, user-specified
std::string aer_cls::abb_get()const{return abb;} // [sng] Aerosol abbreviation
std::string aer_cls::dsc_get()const{return dsc;} // [sng] Aerosol description
std::string aer_cls::fl_idx_rfr_get()const{return idx_rfr->fl_idx_rfr_get();} // [sng] Aerosol input file
void aer_cls::prn()const{std::cout << this;} // [fnc] Print object contents

std::complex<prc_cmp> // O [frc] Refractive index
aer_cls::idx_rfr_get // [fnc] Return refractive index at specified wavelength
(const prc_cmp wvl_ctr) // I [m] Wavelength at band center
  const{
  // Purpose: Wrapper for idx_rfr_cls::idx_rfr_get()
  return idx_rfr->idx_rfr_get(wvl_ctr); // [fnc] Refractive index at band center
} // end idx_rfr_get()

int // O [enm] Return success code
aer_cls::idx_rfr_get // [fnc] Return array of refractive indices
(const prc_cmp *wvl_ctr, // I [m] Wavelength at band center
 std::complex<prc_cmp> *idx_rfr_out, // O [frc] Refractive index
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const bool wrn_ntp_flg) // I [flg] Print WARNINGs from ntp_vec()
  const{
  // Purpose: Wrapper for idx_rfr_cls::idx_rfr_get()
  // NB: nomenclature idx_rfr_out avoids name clash with private member idx_rfr
  int rcd_lcl(0); // [enm] Return success code
  
  rcd_lcl+=idx_rfr->idx_rfr_get
    (wvl_ctr, // I [m] Wavelength at band center
     idx_rfr_out, // O [frc] Refractive index
     wvl_nbr, // I [nbr] Number of wavelength bands
     wrn_ntp_flg); // I [flg] Print WARNINGs from ntp_vec()

  return rcd_lcl; // [enm] Return success code
} // end idx_rfr_get() prototype

// Public member functions end
// Private member functions begin

int // [enm] Return success code
aer_cls::recompute(){ // [fnc] Recompute properties of object
  // Purpose: Recompute properties of object

  // Set private members which lack set() functions
  dsc=abb2dsc(abb); // [sng] Aerosol description
  mmw=abb2mmw(abb); // [kg mol-1] Molecular weight
  spc_heat=abb2spc_heat(abb); // [J kg-1 K-1] Specific heat capacity

  // fxm: add idx_rfr::recompute() here?

  return rcd; // [enm] Return success code
} // end aer_cls::recompute()

// Private member functions end
// Global functions with C++ linkages begin

int // O [rcd] Return success code
cff_drg_get
(const long in_nbr, // I [nbr] Size of arrays
 const prc_cmp *ryn_nbr, // I [frc] Reynolds number
 prc_cmp *cff_drg) // O [frc] Drag coefficient
{
  // Purpose: Compute drag coefficient for given Reynolds number
  // Range of validity is Re < 2.0e5 as reported in SeP97 p. 463 (8.32)
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string sbr_nm("cff_drg_get"); // [sng] Name of subroutine
  long idx; // Counting index
  // Main code
  for(idx=0;idx<in_nbr;idx++){
    if(ryn_nbr[idx] < 1.0e6){
      cff_drg[idx]=cff_drg_fst_scl(ryn_nbr[idx]);
    }else{
      std::cerr << prg_nm_get() << ": WARNING Reynolds number " << ryn_nbr[idx] << " exceeds 1.0e6" << std::endl;
      rcd+=1; // Out of range error
      err_prn(sbr_nm,"Reynolds number exceeds bounds of drag coefficient parameterization");
    } // end else
  } // end loop over lon
  return rcd;
} // end cff_drg_get()

int // O [rcd] Return success code
vlc_stk_get
(const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 const prc_cmp vsc_dyn_atm, // I [kg m-1 s-1] Dynamic viscosity of atmosphere
 prc_cmp *slp_crc, // O [frc] Slip correction factor
 prc_cmp *vlc_stk) // O [m s-1] Stokes' settling velocity (Re < 0.1)
{
  /* Purpose: Compute slip correction factor and Stokes velocity of particles
     Solution for vlc_stk is analytic and "exact" but only realistic when Re < 0.1 
     and when there are no slipstream (near-field) effects.
     Near field effects may occur when particles of like sizes collect eachother
     Because they travel at nearly same speed and thus have ample time to perturb
     the Stokes flow before collecting, validity of Stokes flow is reduced.
     See discussion in PrK98 for details */
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string sbr_nm("vlc_stk_get"); // [sng] Name of subroutine
  long idx; // Counting index
  // Main code
  for(idx=0;idx<sz_nbr;idx++){
    slp_crc[idx]=1.0+2.0*mfp_atm*(1.257+0.4*std::exp(-1.1*dmt_ctr[idx]/(2.0*mfp_atm)))/dmt_ctr[idx]; // [frc] Slip correction factor SeP97 p. 464 (8.34)
    // NB: Approximate (dns_prt - dns_mdp) as dns_prt SeP97 p. 466 (8.42)
    vlc_stk[idx]=(1.0/18.0)*dmt_ctr[idx]*dmt_ctr[idx]*dns_prt*grv_sfc_mean*slp_crc[idx]/vsc_dyn_atm; // [m s-1] Stokes' settling velocity (Re < 0.1) SeP97 p. 466 (8.42)
  } // end loop over sz
  return rcd;
} // end vlc_stk_get()

int // O [rcd] Return success code
vlc_grv_get // [fnc] Terminal fall speed of spherical particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp *slp_crc, // I [frc] Slip correction factor
 const prc_cmp *vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
 prc_cmp *ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
 prc_cmp *stk_crc, // O [frc] Correction to Stokes settling velocity
 prc_cmp *vlc_grv) // O [m s-1] Settling velocity
{
  /* Purpose: Compute terminal fall speed of particles and raindrops using 
     iterative numerical solution to full non-linear equation of motion
     For Reynolds number flows Re < 0.1 Stokes' velocity is valid for vlc_grv SeP97 p. 466 (8.42)
     For larger Re, inertial effects become important and empirical drag coefficients must be employed
     Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
     Using Stokes' velocity rather than iterative solution with empirical drag coefficient causes 60% errors for D = 200 um SeP97 p. 468 */
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp vlc_grv_old; // [m s-1] Previous gravitational settling velocity
  long itr_idx; // [idx] Counting index

  // Main code
  const std::string sbr_nm("vlc_grv_get"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Sanity check
  assert(dns_mdp < phc::dns_mdp_max); // [kg m-3] Make sure aerosol and environmental densities are not switched
  assert(dns_prt > phc::dns_prt_min); // [kg m-3] Make sure aerosol and environmental densities are not switched

  // Iterative solution for drag coefficient, Reynolds number, and terminal velocity
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for vlc_grv is exact for Re < 0.1
    vlc_grv[idx]=vlc_stk[idx]; // [m s-1] Settling velocity
    if(dbg_lvl > dbg_crr){
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s %8s %8s\n",
		    "idx ","itr ","dmt_ctr","  C_c  ","  Re   ","  C_d  "," v_Stk "," v_grv ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s %8s %8s\n",
		    "    ","    ","  um   ","  frc  ","  frc  ","  frc  "," m s-1 "," m s-1 ","  frc  ");
    } // end if dbg
    while(eps_crr > eps_max){
      // Save terminal velocity for convergence test
      vlc_grv_old=vlc_grv[idx]; // [m s-1] 
      ryn_nbr_grv[idx]=vlc_grv[idx]*dmt_ctr[idx]/vsc_knm_atm; // [frc] SeP97 p. 460
      if(ryn_nbr_grv[idx] > 1.0e6){ // [frc] Limit set in aer.hh:cff_drg_fst_scl()
	std::cerr << prg_nm << ": WARNING Reynolds number of particle diameter " << dmt_ctr[idx]*1.0e6 << " um is " << ryn_nbr_grv[idx] << std::endl;
	err_prn(prg_nm,sbr_nm,"Reynolds number exceeds bounds of drag coefficient parameterization");
      } // end else
      // Update drag coefficient based on new Reynolds number
      cff_drg_grv[idx]=cff_drg_fst_scl(ryn_nbr_grv[idx]); // [frc] Drag coefficient at terminal velocity
      // Update terminal velocity based on new Reynolds number and drag coefficient
      vlc_grv[idx]=std::sqrt(4.0*grv_sfc_mean*dmt_ctr[idx]*slp_crc[idx]*(dns_prt-dns_mdp)/(3.0*cff_drg_grv[idx]*dns_mdp)); // [m s-1] Terminal velocity SeP97 p. 467 (8.44)
      eps_crr=PRC_CMP_ABS((vlc_grv[idx]-vlc_grv_old)/vlc_grv[idx]); // Relative convergence
      if(dbg_lvl > dbg_crr){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %7.5f %8.6f %8.6f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_ctr[idx]*1.0e6,slp_crc[idx],ryn_nbr_grv[idx],cff_drg_grv[idx],vlc_stk[idx],vlc_grv[idx],eps_crr);
      } // end if dbg
      if(itr_idx >= 12){
	// Numerical ping-pong may occur when Re = 0.1, 2.0, or 500.0 due to discontinuities in derivative of drag coefficient
	wrn_prn(prg_nm,sbr_nm,"Averaging old and new vlc_grv to attempt to avoid ping pong in iteration "+nbr2sng(itr_idx));
	vlc_grv[idx]=0.5*(vlc_grv[idx]+vlc_grv_old); // [m s-1]
      } // endif
      if(itr_idx > 20){
	std::cerr << "Final: ryn_nbr_grv = " << ryn_nbr_grv[idx] << ", vlc_grv_old = " << vlc_grv_old << " m s-1, vlc_grv = " << vlc_grv[idx] << " m s-1, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Terminal velocity not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
    stk_crc[idx]=vlc_grv[idx]/vlc_stk[idx]; // [frc] Correction to Stokes settling velocity
  } // end loop over sz
  return rcd;
} // end vlc_grv_get()

int // O [rcd] Return success code
vlc_grv_get_Gin03 // [fnc] Terminal fall speed of aspherical particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *asp_rat_lps, // I [frc] Ellipsoidal aspect ratio
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp *slp_crc, // I [frc] Slip correction factor
 const prc_cmp *vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
 prc_cmp *ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
 prc_cmp *stk_crc, // O [frc] Correction to Stokes settling velocity
 prc_cmp *vlc_grv) // O [m s-1] Settling velocity
{
  /* Purpose: Compute terminal fall speed of aspherical particles using 
     iterative numerical solution to full non-linear equation of motion
     For Reynolds number flows Re < 0.1 Stokes' velocity is valid for vlc_grv SeP97 p. 466 (8.42)
     For larger Re, inertial effects become important and empirical drag coefficients must be employed
     Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
     Using Stokes' velocity rather than iterative solution with empirical drag coefficient causes 60% errors for D = 200 um SeP97 p. 468
     mie -dbg --no_mie --asp_rat=2.0 --sz_mnm=0.0 --sz_mxm=50.0 --sz_nbr=100 --sz_grd=lnr
     ncks -C -F -d sz,10.0e-6 -v asp_rat_lps,dmt_ctr,dmt_mjr,dmt_mnr,dmt_aer,dmt_eqv_sfc,dmt_eqv_vlm,dmt_stk,stk_crc,vlc_grv,vlc_stk,vlc_grv_nwr,vlc_grv_vwr ${DATA}/mie/mie.nc */
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp vlc_grv_old; // [m s-1] Previous gravitational settling velocity
  long itr_idx; // [idx] Counting index

  // Variable required by aspherical considerations
  gsl_mode_t gsl_mode(GSL_PREC_DOUBLE); // [enm] GSL precision mode
  gsl_sf_result nsw_dbl; // [frc] GSL result structure
  prc_cmp fct_asp; // [frc] Factor by which spherical vlc_grv^2 differs from spherical
  prc_cmp xcn_lps_prl; // [frc] Eccentricity of ellipsoid
  prc_cmp psi_lps; // [frc] Equivalent diameter divided by semi-minor axis
  prc_cmp sph_fct; // [frc] Sphericity factor

  // Main code
  const std::string sbr_nm("vlc_grv_Gin03_get"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Sanity check
  assert(dns_mdp < phc::dns_mdp_max); // [kg m-3] Make sure aerosol and environmental densities are not switched
  assert(dns_prt > phc::dns_prt_min); // [kg m-3] Make sure aerosol and environmental densities are not switched

  // Iterative solution for drag coefficient, Reynolds number, and terminal velocity
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for vlc_grv is exact for Re < 0.1
    vlc_grv[idx]=vlc_stk[idx]; // [m s-1] Settling velocity
    xcn_lps_prl=xcn_lps_prl_fst_scl(asp_rat_lps[idx]); // [frc] Eccentricity of ellipsoid
    assert(xcn_lps_prl < 1.0);
    psi_lps=psi_lps_fst_scl(asp_rat_lps[idx]); // [frc] Equivalent diameter divided by semi-minor axis Gin03 p. 2 (10)
    sph_fct=sph_fct_fst_scl(asp_rat_lps[idx]); // [frc] Sphericity factor
    if(dbg_lvl > dbg_crr){
      (void)std::fprintf(stderr,"%4s %4s %7s %8s %8s %8s %8s %8s %7s %8s %8s %8s %8s %8s\n",
		    "idx ","itr ","dmt_ctr","asp_rat","xcn_lps","sph_fct","psi_lps"," CEISK ","  C_c  ","  Re   ","  C_d  "," v_Stk "," v_grv ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %8s %8s %8s %8s %8s %7s %8s %8s %8s %8s %8s\n",
		    "    ","    ","  um   ","  frc  ","  frc  ","  frc  ","  frc  ","  frc  ","  frc  ","  frc  ","  frc  "," m s-1 "," m s-1 ","  frc  ");
    } // end if dbg
    while(eps_crr > eps_max){
      // Save terminal velocity for convergence test
      vlc_grv_old=vlc_grv[idx]; // [m s-1] 
      ryn_nbr_grv[idx]=vlc_grv[idx]*dmt_ctr[idx]/vsc_knm_atm; // [frc] SeP97 p. 460
      if(ryn_nbr_grv[idx] > 1.0e6){ // [frc] Limit set in aer.hh:cff_drg_fst_scl()
	std::cerr << prg_nm << ": WARNING Reynolds number of particle diameter " << dmt_ctr[idx]*1.0e6 << " um is " << ryn_nbr_grv[idx] << std::endl;
	err_prn(prg_nm,sbr_nm,"Reynolds number exceeds bounds of drag coefficient parameterization");
      } // end else
      // Update drag coefficient based on new Reynolds number
      cff_drg_grv[idx]=cff_drg_Boo71_fst_scl(ryn_nbr_grv[idx],sph_fct); // [frc] Drag coefficient at terminal velocity
      rcd+=gsl_sf_ellint_Ecomp_e((double)xcn_lps_prl,gsl_mode,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      /* Aspherical vlc_grv^2 differs from spherical vlc_grv^2 by factor of fct_asp
	 \lim_{asp_rat_lps \rightarrow 1} fct_asp = 1 since
	 \lim_{asp_rat_lps \rightarrow 1} CEISK(xcn_lps_prl) = pi/2
	 where CEISK = Complete Elliptic Integral of the Second Kind and 
	 \lim_{asp_rat_lps \rightarrow 1} psi_lps = 2 */
      fct_asp=mth::cst_M_PIl/(nsw_dbl.val*psi_lps); // [frc] Factor by which spherical vlc_grv^2 differs from spherical
      // Update terminal velocity based on new Reynolds number and drag coefficient
      // fxm: change both soln's to (dns_prt-dns_mdp) in numerator
      vlc_grv[idx]=std::sqrt(4.0*grv_sfc_mean*dmt_ctr[idx]*slp_crc[idx]*(dns_prt-dns_mdp)*fct_asp/(3.0*cff_drg_grv[idx]*dns_mdp)); // [m s-1] Terminal velocity SeP97 p. 467 (8.44)
      eps_crr=PRC_CMP_ABS((vlc_grv[idx]-vlc_grv_old)/vlc_grv[idx]); // Relative convergence
      if(dbg_lvl > dbg_crr){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %8.6f %8.6f %8.6f %8.6f %8.6f %7.5f %8.6f %8.6f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_ctr[idx]*1.0e6,asp_rat_lps[idx],xcn_lps_prl,sph_fct,psi_lps,nsw_dbl.val,slp_crc[idx],ryn_nbr_grv[idx],cff_drg_grv[idx],vlc_stk[idx],vlc_grv[idx],eps_crr);
      } // end if dbg
      if(itr_idx >= 12){
	// Numerical ping-pong may occur when Re = 0.1, 2.0, or 500.0 due to discontinuities in derivative of drag coefficient
	wrn_prn(prg_nm,sbr_nm,"Averaging old and new vlc_grv to attempt to avoid ping pong in iteration "+nbr2sng(itr_idx));
	vlc_grv[idx]=0.5*(vlc_grv[idx]+vlc_grv_old); // [m s-1]
      } // endif
      if(itr_idx > 20){
	std::cerr << "Final: ryn_nbr_grv = " << ryn_nbr_grv[idx] << ", vlc_grv_old = " << vlc_grv_old << " m s-1, vlc_grv = " << vlc_grv[idx] << " m s-1, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Terminal velocity not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
    stk_crc[idx]=vlc_grv[idx]/vlc_stk[idx]; // [frc] Correction to Stokes settling velocity
  } // end loop over sz
  return rcd;
} // end vlc_grv_Gin03_get()

int // O [rcd] Return success code
asp_rat_lps_get // [fnc] Determine particle aspect ratio
(const long sz_nbr, // I [nbr] Size of arrays
 const std::string cmp_sng_prt, // I [sng] Composition of particle
 const prc_cmp asp_rat_lps_dfl, // I [frc] Ellipsoidal aspect ratio
 prc_cmp *asp_rat_lps) // O [frc] Ellipsoidal aspect ratio
{
  // Purpose: Determine particle aspect ratio
  int rcd(0); // O [rcd] Return success code
  /* Gin03 shows that aspect ratio may change with size for realistic minerals
     Ignore this for now and prescribe size-independent aspect ratio */
  for(long sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    asp_rat_lps[sz_idx]=asp_rat_lps_dfl; // [frc] Ellipsoidal aspect ratio
  } // end loop over idx
  return rcd;
}// end asp_rat_lps_get()

int // O [rcd] Return success code
dmt_aer_get // [fnc] Compute aerodynamic diameter of particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_stk, // I [m] Stokes diameter
 const prc_cmp *slp_crc, // I [frc] Slip correction factor
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 prc_cmp *dmt_aer) // O [m] Aerodynamic diameter
{
  /* Purpose: Compute aerodynamic diameter of particles using iterative numerical solution to full implicit equation
     Aerodynamic diameter is diameter of a unit density sphere with same terminal velocity as particle
     Input requires particle diameter which, for generality, we call the Stokes diameter
     Stokes diameter is diameter of sphere with same terminal settling velocity and density as particle
     Stokes diameter equals particle diameter when particle is spherical
     Input slip correction factor is for Stokes diameter grid and does not change
     Slip correction factor for aerodynamic diameter is computed interactively */

  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  const prc_cmp dns_rat(dns_prt/1000.0); // [frc] Density ratio
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp dmt_aer_old; // [m] Previous aerodynamic diameter
  prc_cmp slp_crc_aer; // [frc] Slip correction factor for aerodynamic size/density
  long itr_idx; // [idx] Counting index

  // Main code
  const std::string sbr_nm("dmt_aer_get"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Sanity check
  assert(mfp_atm < 1.0); // [m] Mean free path of atmosphere
  assert(dns_prt > phc::dns_prt_min); // [kg m-3] Make sure aerosol and environmental densities are not switched

  // Iterative solution for aerodynamic diameter
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for dmt_aer is Stokes diameter
    dmt_aer[idx]=dmt_stk[idx]; // [m] Aerodynamic diameter
    if(dbg_lvl == dbg_old){
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "idx ","itr ","dmt_stk","  C_c  ","dmt_aer","  C_c  ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "    ","    ","  um   ","  frc  ","  um   ","  frc  ","  frc  ");
    } // end if dbg
    while(eps_crr > eps_max){
      // Save aerodynamic diameter for convergence test
      dmt_aer_old=dmt_aer[idx]; // [m] Aerodynamic diameter
      // Update slip correction based on new aerodynamic diameter
      slp_crc_aer=slp_crc_fst_scl(mfp_atm,dmt_aer[idx]); // [frc] Slip correction factor
      // Update aerodynamic diameter based on new slip correction factor
      // dns_rat is ratio of particle density to unity density in CGS ( = 1000 kg m-3)
      dmt_aer[idx]=dmt_stk[idx]*std::sqrt(dns_rat*slp_crc[idx]/slp_crc_aer); // [m] Aerodynamic diameter SeP97 p. 488 (8.106)
      eps_crr=PRC_CMP_ABS((dmt_aer[idx]-dmt_aer_old)/dmt_aer[idx]); // Relative convergence
      if(dbg_lvl == dbg_old){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %7.5f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_stk[idx]*1.0e6,slp_crc[idx],dmt_aer[idx]*1.0e6,slp_crc_aer,eps_crr);
      } // end if dbg
      if(itr_idx == 15){
	// Routine converges slowly, 12 iterations is common for small, heavy particles
	wrn_prn(prg_nm,sbr_nm,"Averaging old and new dmt_aer in a vain attempt to avoid ping pong...");
	dmt_aer[idx]=0.5*(dmt_aer[idx]+dmt_aer_old); // [m s-1]
      } // endif
      if(itr_idx > 20){
	std::cerr << "Final: dmt_aer_old = " << dmt_aer_old*1.0e6 << " um, dmt_aer = " << dmt_aer[idx]*1.0e6 << " um, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Aerodynamic diameter not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
  } // end loop over sz
  return rcd;
} // end dmt_aer_get()

int // O [rcd] Return success code
dmt_stk_get // [fnc] Compute Stokes diameter of particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 const prc_cmp *vlc_grv, // I [m s-1] Settling velocity
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *dmt_stk) // O [m] Stokes diameter
{
  /* Purpose: Compute Stokes diameter of particles using iterative numerical solution to full implicit equation
     Stokes diameter is diameter of sphere with same terminal settling velocity and density as particle
     Stokes diameter equals particle diameter when particle is spherical
     Input slip correction factor is for Stokes diameter grid and does not change
     Slip correction factor for aerodynamic diameter is computed interactively */
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp dmt_stk_old; // [m] Previous Stokes diameter
  prc_cmp cff_drg_grv; // [frc] Drag coefficient at terminal velocity
  prc_cmp ryn_nbr_grv; // [frc] Reynolds number at terminal velocity
  prc_cmp slp_crc; // [frc] Slip correction factor
  long itr_idx; // [idx] Counting index

  // Main code
  const std::string sbr_nm("dmt_stk_get"); // [sng] Subroutine name
  const std::string prg_nm(prg_nm_get()); // [sng] Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // [enm] Debugging level

  // Sanity check
  assert(mfp_atm < 1.0); // [m] Mean free path of atmosphere
  assert(dns_prt > phc::dns_prt_min); // [kg m-3] Make sure aerosol and environmental densities are not switched

  // Iterative solution for aerodynamic diameter
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for dmt_stk is bin center diameter
    dmt_stk[idx]=dmt_ctr[idx]; // [m] Aerodynamic diameter
    if(dbg_lvl > dbg_crr){
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "idx ","itr ","dmt_stk","  C_c  ","  Re   ","  C_D  ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "    ","    ","  um   ","  frc  ","  frc  ","  frc  ","  frc  ");
    } // end if dbg
    while(eps_crr > eps_max){
      // Save settling speed for convergence test
      dmt_stk_old=dmt_stk[idx]; // [m s-1] Previous gravitational settling velocity
      // Update slip correction based on new Stokes diameter
      slp_crc=slp_crc_fst_scl(mfp_atm,dmt_stk[idx]); // [frc] Slip correction factor
      ryn_nbr_grv=vlc_grv[idx]*dmt_stk[idx]/vsc_knm_atm; // [frc] SeP97 p. 460
      // Update drag coefficient based on new Reynolds number
      cff_drg_grv=cff_drg_fst_scl(ryn_nbr_grv); // [frc] Drag coefficient at terminal velocity
      // Update Stokes diameter based on new settling velocity
      dmt_stk[idx]=3.0*cff_drg_grv*dns_mdp*vlc_grv[idx]*vlc_grv[idx]/(4.0*grv_sfc_mean*slp_crc*(dns_prt-dns_mdp)); // [m] Stokes diameter SeP97 p. 488 (8.106)
      eps_crr=PRC_CMP_ABS((dmt_stk[idx]-dmt_stk_old)/dmt_stk[idx]); // Relative convergence
      if(dbg_lvl > dbg_crr){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %7.5f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_stk[idx]*1.0e6,slp_crc,ryn_nbr_grv,cff_drg_grv,eps_crr);
      } // end if dbg
      if(itr_idx == 15){
	// Routine converges slowly, 12 iterations is common for small, heavy particles
	wrn_prn(prg_nm,sbr_nm,"Averaging old and new dmt_stk in a vain attempt to avoid ping pong...");
	dmt_stk[idx]=0.5*(dmt_stk[idx]+dmt_stk_old); // [m s-1]
      } // endif
      if(itr_idx > 20){
	std::cerr << "Final: dmt_stk_old = " << dmt_stk_old*1.0e6 << " um, dmt_stk = " << dmt_stk[idx]*1.0e6 << " um, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Stokes diameter not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
  } // end loop over sz
  return rcd;
} // end dmt_stk_get()

int // O [rcd] Return success code
dlc_fnc_ffc_brg_get // [fnc] Compute aerodynamic diameter of particles
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_stk, // I [m] Stokes diameter
 const prc_cmp *slp_crc, // I [frc] Slip correction factor
 const prc_cmp dns_prt, // I [kg m-3] Particle density
 const prc_cmp mfp_atm, // I [m] Mean free path of atmosphere
 prc_cmp *dmt_aer) // O [m] Aerodynamic diameter
{
  /* Purpose: Compute aerodynamic diameter of particles using iterative numerical solution to full implicit equation
     Aerodynamic diameter is diameter of a unit density sphere with same terminal velocity as particle
     Input requires particle diameter which, for generality, we call the Stokes diameter
     Stokes diameter is diameter of a sphere with same terminal settling velocity and density as particle
     Stokes diameter equals particle diameter when particle is spherical
     Input slip correction factor is for Stokes diameter grid and does not change
     Slip correction factor for aerodynamic diameter is computed interactively */

  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  const prc_cmp dns_rat(dns_prt/1000.0); // [frc] Density ratio
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp dmt_aer_old; // [m s-1] Previous gravitational settling velocity
  prc_cmp slp_crc_aer; // [frc] Slip correction factor for aerodynamic size/density
  long itr_idx; // [idx] Counting index

  // Main code
  const std::string sbr_nm("dlc_fnc_ffc_brg_get"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Sanity check
  assert(mfp_atm < 1.0); // [m] Mean free path of atmosphere
  assert(dns_prt > phc::dns_prt_min); // [kg m-3] Make sure aerosol and environmental densities are not switched

  // Iterative solution for aerodynamic diameter
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for dmt_aer is Stokes diameter
    dmt_aer[idx]=dmt_stk[idx]; // [m] Aerodynamic diameter
    if(dbg_lvl > dbg_crr){
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "idx ","itr ","dmt_stk","  C_c  ","dmt_aer","  C_c  ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s\n",
		    "    ","    ","  um   ","  frc  ","  um   ","  frc  ","  frc  ");
    } // end if dbg
    while(eps_crr > eps_max){
      // Save aerodynamic diameter for convergence test
      dmt_aer_old=dmt_aer[idx]; // [m] Aerodynamic diameter
      // Update slip correction based on new aerodynamic diameter
      slp_crc_aer=slp_crc_fst_scl(mfp_atm,dmt_aer[idx]); // [frc] Slip correction factor
      // Update aerodynamic diameter based on new slip correction factor
      // dns_rat is ratio of particle density to unity density in CGS ( = 1000 kg m-3)
      dmt_aer[idx]=dmt_stk[idx]*std::sqrt(dns_rat*slp_crc[idx]/slp_crc_aer); // [m] Aerodynamic diameter SeP97 p. 488 (8.106)
      eps_crr=PRC_CMP_ABS((dmt_aer[idx]-dmt_aer_old)/dmt_aer[idx]); // Relative convergence
      if(dbg_lvl > dbg_crr){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %7.5f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_stk[idx]*1.0e6,slp_crc[idx],dmt_aer[idx]*1.0e6,slp_crc_aer,eps_crr);
      } // end if dbg
      if(itr_idx == 15){
	// Routine converges slowly, 12 iterations is common for small, heavy particles
	wrn_prn(prg_nm,sbr_nm,"Averaging old and new dmt_aer in a vain attempt to avoid ping pong...");
	dmt_aer[idx]=0.5*(dmt_aer[idx]+dmt_aer_old); // [m s-1]
      } // endif
      if(itr_idx > 20){
	std::cerr << "Final: dmt_aer_old = " << dmt_aer_old*1.0e6 << " um, dmt_aer = " << dmt_aer[idx]*1.0e6 << " um, eps_crr = " << eps_crr << std::endl;
	wrn_prn(prg_nm,sbr_nm,"Aerodynamic diameter not converging, breaking loop...");
	break;
      } // endif
      itr_idx++;
    } // end loop over itr
  } // end loop over sz
  return rcd;
} // end dlc_fnc_ffc_brg_get()

int // O [rcd] Return success code
vnt_mss_get // [fnc] Mass ventilation coefficient
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *shm_nbr, // I [frc] Schmidt number
 const prc_cmp *ryn_nbr, // I [frc] Reynolds number
 prc_cmp *vnt_mss) // O [fnc] Mass ventilation coefficient
{
  /* Purpose: Compute mass ventilation coefficient for given fluid numbers
     Taken from PrK98 p. 541 (13-60) */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long sz_idx; // [idx] Counting index
  // Main code
  for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    vnt_mss[sz_idx]= // O [frc] Mass ventilation coefficient
      vnt_cff_fst_scl // [fnc] Mass or heat ventilation coefficient
      (shm_nbr[sz_idx], // I [frc] Schmidt number
       ryn_nbr[sz_idx]); // I [frc] Reynolds number
  } // end loop over sz
  return rcd;
} // end vnt_mss_get()

