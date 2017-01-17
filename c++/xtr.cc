// $Id$ 

// Implementation (declaration) of extrapolation class

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <xtr.hh> // Class-specific definitions

// Xtr_cls class

// Non-static members
Xtr_cls::Xtr_cls(const std::string &opt_sng) // Overload constructor
{
  // Save descriptive string
  xtr_sng=opt_sng;
  // Set all flags
  flg_set(xtr_sng);
} // end Default constructor

void Xtr_cls::flg_fls(){ // Set flags to false
  // Initialize all flags to false
  xtr_err=false; // [flg] No extrapolation is allowed
  xtr_fll=false; // [flg] Full extrapolation is allowed (implies xtr_prt)
  xtr_fll_lnr=false; // [flg] Perform linear extrapolation using two nearest valid neighbors
  xtr_fll_ngh=false; // [flg] Set extrapolated value to value of nearest valid neighbor
  xtr_fll_nil=false; // [flg] Set extrapolated value to 0.0
  xtr_prt=false; // [flg] Partial extrapolation is allowed
  xtr_prt_frc=false; // [flg] Use average value of overlap region
  xtr_prt_lnr=false; // [flg] Perform linear extrapolation using two nearest valid neighbors
  xtr_prt_ngh=false; // [flg] Set extrapolated value to value of nearest valid neighbor
  xtr_prt_nil=false; // [flg] Set extrapolated value to 0.0
  xtr_prt_wgt=false; // [flg] Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc)
  xtr_vrb=false; // [flg] Print verbose warning when extrapolation is performed
} // end Xtr::flg_fls()

void Xtr_cls::flg_set(const std::string &xtr_sng){ // Set flags
  unsigned long flg_crr_end;
  long flg_crr_len;
  sng2sng_map::const_iterator itr;
  std::string flg_crr_sng;
  std::string opt_lst_crr;
  std::string sbr_nm("Xtr_cls::flg_set");

  // Initialize all flags to false
  flg_fls();

  // Loop over flags
  unsigned long flg_crr_srt=xtr_sng.find("xtr");
  while(flg_crr_srt < std::string::npos){
    // Isolate current flag
    flg_crr_end=xtr_sng.find_first_of(" , + \t \n ;",flg_crr_srt+1); // These are allowed separators
    if(flg_crr_end < std::string::npos) flg_crr_len=flg_crr_end-flg_crr_srt; else flg_crr_len=xtr_sng.size()-flg_crr_srt;
    flg_crr_sng=xtr_sng.substr(flg_crr_srt,flg_crr_len);

    // Search flag abbreviations for current flag string
    itr=abb2dsc_map.find(flg_crr_sng);
    if(itr == abb2dsc_map.end()) err_prn(sbr_nm,flg_crr_sng+" is unknown");

    if(dbg_lvl_get() == dbg_io) err_prn(sbr_nm,"Setting flag: "+flg_crr_sng+", Description: "+itr->second);

    // Set the appropriate flag
    if(flg_crr_sng == "xtr_err") xtr_err=true; // [flg] No extrapolation is allowed
    if(flg_crr_sng == "xtr_fll") xtr_fll=true; // [flg] Full extrapolation is allowed (implies xtr_prt)
    if(flg_crr_sng == "xtr_fll_lnr") xtr_fll_lnr=true; // [flg] Perform linear extrapolation using two nearest valid neighbors
    if(flg_crr_sng == "xtr_fll_ngh") xtr_fll_ngh=true; // [flg] Set extrapolated value to value of nearest valid neighbor
    if(flg_crr_sng == "xtr_fll_nil") xtr_fll_nil=true; // [flg] Set extrapolated value to 0.0
    if(flg_crr_sng == "xtr_prt") xtr_prt=true; // [flg] Partial extrapolation is allowed
    if(flg_crr_sng == "xtr_prt_frc") xtr_prt_frc=true; // [flg] Use average value of overlap region
    if(flg_crr_sng == "xtr_prt_lnr") xtr_prt_lnr=true; // [flg] Perform linear extrapolation using two nearest valid neighbors
    if(flg_crr_sng == "xtr_prt_ngh") xtr_prt_ngh=true; // [flg] Set extrapolated value to value of nearest valid neighbor
    if(flg_crr_sng == "xtr_prt_nil") xtr_prt_nil=true; // [flg] Set extrapolated value to 0.0
    if(flg_crr_sng == "xtr_prt_wgt") xtr_prt_wgt=true; // [flg] Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc)
    if(flg_crr_sng == "xtr_vrb") xtr_vrb=true; // [flg] Print verbose warning when extrapolation is performed

    // Is there another flag to process?
    flg_crr_srt=xtr_sng.find("xtr",flg_crr_srt+1);
  } // end loop over flags

  flg_vld(); // Validate flags

} // end Xtr_cls::flg_set()

void Xtr_cls::flg_vld(){ // Validate flags
  // Purpose: Perform a consistency check on the flags
  // NB: nil, ngh, lnr, are mutually exclusive, frc implies prt, and wgt implies frc
  if(xtr_fll_nil || xtr_fll_ngh || xtr_fll_lnr) xtr_fll=true; else xtr_fll=false;
  if(xtr_prt_nil || xtr_prt_ngh || xtr_prt_lnr || xtr_prt_frc || xtr_prt_wgt) xtr_prt=true; else xtr_prt=false;
  if(xtr_prt_wgt) xtr_prt_frc=true;
} // end Xtr_cls::flg_vld()

// Initialize static members
sng2sng_map Xtr_cls::abb2dsc_map=Xtr_cls::abb2dsc_map_mk();
sng2sng_map Xtr_cls::abb2dsc_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;
  map_tmp.insert(sng2sng_map::value_type("xtr_err","No extrapolation is allowed"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt","Partial extrapolation is allowed"));
  map_tmp.insert(sng2sng_map::value_type("xtr_fll","Full extrapolation is allowed (implies xtr_prt)"));
  map_tmp.insert(sng2sng_map::value_type("xtr_fll_nil","Set extrapolated value to 0.0"));
  map_tmp.insert(sng2sng_map::value_type("xtr_fll_ngh","Set extrapolated value to value of nearest valid neighbor"));
  map_tmp.insert(sng2sng_map::value_type("xtr_fll_lnr","Perform linear extrapolation using two nearest valid neighbors"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt_nil","Set extrapolated value to 0.0"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt_ngh","Set extrapolated value to value of nearest valid neighbor"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt_lnr","Perform linear extrapolation using two nearest valid neighbors"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt_frc","Use average value of overlap region"));
  map_tmp.insert(sng2sng_map::value_type("xtr_prt_wgt","Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc)"));
  map_tmp.insert(sng2sng_map::value_type("xtr_vrb","Print verbose warnings whenever extrapolation is performed"));
  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end Xtr::abb2dsc_map_mk()

std::string Xtr_cls::abb2dsc(const std::string &flg_sng){ // Abbreviation to description mapper
  sng2sng_map::const_iterator itr;
  itr=abb2dsc_map.find(flg_sng);
  std::string sbr_nm("Xtr_cls::abb2dsc");
  if(itr == abb2dsc_map.end()) err_prn(sbr_nm,flg_sng+" is unknown");
  // NB: Return a value to initialize a static class member
  return itr->second;
} // end Xtr_cls::abb2dsc()
