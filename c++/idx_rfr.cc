// $Id$ 

// Purpose: Implementation (declaration) of refractive index classes

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <idx_rfr.hh> // Refractive indices

// idx_rfr_cls class

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const idx_rfr_cls &idx_rfr_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for idx_rfr_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     cout << idx_rfr_obj; */
  srm_out << "Public contents of idx_rfr_obj: " << std::endl;
  srm_out << "idx_rfr_obj.nst_nbr_get() = " << idx_rfr_obj.nst_nbr_get() << std::endl;
  srm_out << "idx_rfr_obj.abb_get() = " << idx_rfr_obj.abb_get() << std::endl;
  srm_out << "idx_rfr_obj.dsc_get() = " << idx_rfr_obj.dsc_get() << std::endl;
  srm_out << "idx_rfr_obj.fl_idx_rfr_get() = " << idx_rfr_obj.fl_idx_rfr_get() << std::endl;
  srm_out << "idx_rfr_obj.idx_rfr_fnc_ptr_get() = " << idx_rfr_obj.idx_rfr_fnc_ptr_get() << std::endl;
  srm_out << "idx_rfr_obj.tpt_get() = " << idx_rfr_obj.tpt_get() << std::endl;
  srm_out << "Private contents of idx_rfr_obj: " << std::endl;
  srm_out << "idx_rfr_usr_flg = " << idx_rfr_obj.idx_rfr_usr_flg << std::endl;
  if(idx_rfr_obj.idx_rfr_usr_flg){
    srm_out << "idx_rfr_usr = " << idx_rfr_obj.idx_rfr_usr << std::endl;
  } // endif user-specified refractive indices
  if(idx_rfr_obj.flg_raw_idx_rfr_in_dyn_mmr){
    srm_out << "fl_idx_rfr = " << idx_rfr_obj.fl_idx_rfr << std::endl;
    srm_out << "wvl_nbr_tbl = " << idx_rfr_obj.wvl_nbr_tbl << std::endl;
    srm_out << "wvl_ctr_tbl[0] = " << idx_rfr_obj.wvl_ctr_tbl[0] << std::endl;
    srm_out << "wvl_ctr_tbl[wvl_nbr_tbl-1] = " << idx_rfr_obj.wvl_ctr_tbl[idx_rfr_obj.wvl_nbr_tbl-1] << std::endl;
    idx_rfr_obj.idx_rfr_prn(idx_rfr_obj.wvl_ctr_tbl,idx_rfr_obj.wvl_nbr_tbl,idx_rfr_obj.idx_rfr_tbl);
    // if(dbg_lvl_get() >= dbg_vrb) idx_rfr_prn(wvl_ctr,wvl_nbr,idx_rfr);
    // srm_out << "idx_rfr_tbl[0] = " << idx_rfr_obj.idx_rfr_tbl[0] << std::endl;
    // srm_out << "idx_rfr_tbl[wvl_nbr_tbl-1] = " << idx_rfr_obj.idx_rfr_tbl[idx_rfr_obj.wvl_nbr_tbl-1] << std::endl;
  }else{
    srm_out << "No raw refractive indices in dynamic memory" << std::endl;
  } // endif
    
  return srm_out; // [srm] Reference to output stream for cascading
} // end idx_rfr_cls::operator<<()

// Friendly functions end
// Static members begin

const prc_cmp idx_rfr_cls::tpt_dfl=300.0; // [K] Default temperature
const std::string idx_rfr_cls::idx_rfr_typ_dfl="lac_ChC90"; // [sng] Default refractive index abbreviation
const std::complex<prc_cmp> idx_rfr_cls::idx_rfr_usr_dfl=sng2cpx("1.0+0.0i"); // [frc] Default refractive index, user-specified
int idx_rfr_cls::nst_nbr=0; // [nbr] Number of instantiated class members
sng2idx_rfr_map idx_rfr_cls::abb2obj_map; // [map] Abbreviation to object map

idx_rfr_fnc_ptr_typ // [fnc] Pointer to function returning refractive indices
idx_rfr_cls::abb2fnc // [fnc] Abbreviation to function mapper
(const std::string &abb_sng){ // [fnc] Abbreviation to function mapper
  // Purpose: Return pointer to function which computes refractive indices
  return abb2obj_map.find(abb_sng)->second->idx_rfr_fnc_ptr; // [fnc] Pointer to function returning refractive indices
} // end idx_rfr_cls::abb2fnc()

// Static members end
// Static member functions begin

int idx_rfr_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members

int // O [enm] Return success code
idx_rfr_cls::tst // [fnc] Self-test idx_rfr_cls class
(const std::string &abb_sng, // [fnc] Refractive index abbreviation
 const long obj_nbr) // [nbr] Number of objects to create/destroy
{
  // Purpose: Perform self test of idx_rfr_cls class
  int rcd_lcl(0); // [enm] Return success code
  std::string sbr_nm("idx_rfr_cls::tst"); // [sng] Subroutine name

  // Test for memory leaks
  std::cout << sbr_nm+"()"+" is testing for memory leaks by creating and destroying " << obj_nbr << " idx_rfr_cls objects..." << std::endl;
  long idx; // [idx] Counting index
  idx_rfr_cls *tst_obj; // [sct] Test object
  const prc_cmp wvl_ctr(0.5e-6); // [m] Wavelength at band center
  std::complex<prc_cmp> idx_rfr; // [frc] Refractive index
  for(idx=0;idx<obj_nbr;idx++){
    tst_obj=new idx_rfr_cls(abb_sng); // [sct] Test object
    idx_rfr=tst_obj->idx_rfr_get(wvl_ctr); // [frc] Refractive index
    std::cout << "idx = " << idx << ", nst_nbr = " << nst_nbr << std::endl;
    delete tst_obj; // [sct] Test object
  } // end loop over obj

  // fxm Print all refractive indices
  tst_obj=new idx_rfr_cls("lac_ChC90"); // [sct] Test object
  idx_rfr=tst_obj->idx_rfr_get(wvl_ctr); // [frc] Refractive index
  std::cout << "abb = " << tst_obj->abb << ", wvl_ctr = " << wvl_ctr << ", idx_rfr = " << idx_rfr << std::endl;
  delete tst_obj; // [sct] Test object

  return rcd_lcl; // [enm] Return success code
} // end idx_rfr_cls::tst()

idx_rfr_cls * // [ptr] Pointer to refractive index object
idx_rfr_cls::abb2obj // [fnc] Abbreviation to object mapper
(const std::string &abb_sng){ // [sng] Abbreviation string
  // Purpose: Return pointer to refractive index object
  return abb2obj_map.find(abb_sng)->second; // [fnc] Pointer to refractive index object
} // end idx_rfr_cls::abb2obj()

int // [enm] Return success code
idx_rfr_cls::abb2obj_map_add // [fnc] Register new refractive index object
(idx_rfr_cls *idx_rfr_ptr) // [obj] Refractive index object
{
  // Purpose: Register refractive index object with list of existing objects
  int rcd_lcl(0); // [enm] Return success code
  assert(idx_rfr_ptr != NULL);
  std::string sbr_nm("idx_rfr_cls::abb2obj_map_add"); // [sng] Subroutine name
  if(dbg_lvl_get() == dbg_crr) err_prn(sbr_nm,"Registering "+idx_rfr_ptr->abb_get());
  abb2obj_map.insert(sng2idx_rfr_map::value_type(idx_rfr_ptr->abb_get(),idx_rfr_ptr)); // [map] Abbreviation to object map
  return rcd_lcl; // [enm] Return success code
} // end idx_rfr_cls::abb2obj_map_add()

// Static member functions end
// Public member functions begin

idx_rfr_cls::idx_rfr_cls // [fnc] Constructor
(const std::string &abb_sng_arg, // [sng] Refractive index abbreviation
 const std::string &fl_idx_rfr_arg, // [sng] Aerosol input file
 const prc_cmp &tpt_arg, // [K] Temperature
 const bool idx_rfr_usr_flg_arg, // [flg] Refractive index is user-specified
 const std::complex<prc_cmp> idx_rfr_usr_arg) // [frc] Refractive index, user-specified
{
  // Purpose: Constructor
  std::string sbr_nm("idx_rfr_cls::idx_rfr_cls"); // [sng] Subroutine name

  rcd=0; // [enm] Return success code

  flg_raw_idx_rfr_in_dyn_mmr=false; // [flg] Raw refractive indices from file are in dynamic memory

  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls

  // Set type before file
  rcd+=abb_set(abb_sng_arg); // [sng] Refractive index type

  // Constructor defaults to reasonable value to allow for calls that omit tpt
  rcd+=tpt_set(tpt_arg); // [K] Temperature

  idx_rfr_usr_flg=idx_rfr_usr_flg_arg; // [flg] Refractive index is user-specified

  if(idx_rfr_usr_flg){
    idx_rfr_usr=idx_rfr_usr_arg; // [frc] Refractive index, user-specified
  }else{
    // Load raw index of refraction data into memory
    rcd+=fl_idx_rfr_set(fl_idx_rfr_arg); // [sng] File containing refractive indices
  } // not idx_rfr_usr_flg

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  if(rcd != 0) err_prn(sbr_nm,"Failed constructor"); // [fnc] Print uniform error message and exit

  nst_nbr++; // [nbr] Number of instantiated class members
} // end idx_rfr_cls constructor

idx_rfr_cls::~idx_rfr_cls(){ // [fnc] Destructor
  rcd+=deallocate(); // [fnc] Free dynamic memory for object
  nst_nbr--; // [nbr] Number of instantiated class members
} // end idx_rfr_cls destructor

idx_rfr_fnc_ptr_typ idx_rfr_cls::idx_rfr_fnc_ptr_get()const{return idx_rfr_fnc_ptr;} // [fnc] Function to compute refractive indices
prc_cmp idx_rfr_cls::tpt_get()const{return tpt;} // [K] Temperature 
std::string idx_rfr_cls::abb_get()const{return abb;} // [sng] Refractive index abbreviation
std::string idx_rfr_cls::dsc_get()const{return dsc;} // [sng] Refractive index description
std::string idx_rfr_cls::fl_idx_rfr_get()const{return fl_idx_rfr;} // [sng] File containing refractive indices
std::complex<prc_cmp> idx_rfr_cls::idx_rfr_usr_get()const{return idx_rfr_usr;} // [frc] Refractive index, user-specified
void idx_rfr_cls::prn()const{std::cout << this;} // [fnc] Print object contents

int // [enm] Return success code
idx_rfr_cls::idx_rfr_fnc_ptr_set(const idx_rfr_fnc_ptr_typ &idx_rfr_fnc_ptr_arg){ // [fnc] Function to compute refractive indices
  std::string sbr_nm("idx_rfr_cls::idx_rfr_fnc_ptr_set"); // [sng] Subroutine name
  idx_rfr_fnc_ptr=idx_rfr_fnc_ptr_arg; // [fnc] Function to compute refractive indices
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::idx_rfr_fnc_ptr_set()

int // [enm] Return success code
idx_rfr_cls::abb_set(const std::string &abb_sng_arg){ // [fnc] Set refractive index type
  std::string sbr_nm("idx_rfr_cls::abb_set"); // [sng] Subroutine name
  abb= (abb_sng_arg == "") ? idx_rfr_typ_dfl : abb_sng_arg; // [sng] Refractive index abbreviation
  /* If constructor is called with abbreviation string next line is redundant
     However, this allows constructor to be called with option string, as well */
  abb=aer_cls::opt2abb(abb); // [sng] Refractive index abbreviation
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::abb_set()

int // [enm] Return success code
idx_rfr_cls::tpt_set(const prc_cmp &tpt_arg) // [K] Temperature
{
  // Purpose: Set particle temperature
  std::string sbr_nm("idx_rfr_cls::tpt_set"); // [sng] Subroutine name

  assert(tpt_arg > 0.0 && tpt_arg <= 400.0); // [frc] Allowed temperature range
  tpt=tpt_arg; // [K] Temperature
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::tpt_set()

int // [enm] Return success code
idx_rfr_cls::idx_rfr_usr_set(const std::complex<prc_cmp> &idx_rfr_usr_arg) // [frc] Refractive index, user-specified
{
  // Purpose: Set user-specified refractive index
  std::string sbr_nm("idx_rfr_cls::idx_rfr_usr_set"); // [sng] Subroutine name

  assert(idx_rfr_usr_arg.imag() > 0.0 && idx_rfr_usr_arg.real() > 0.0); // [frc] Allowed refractive index range
  idx_rfr_usr=idx_rfr_usr_arg; // [frc] Refractive index, user-specified
  idx_rfr_usr_flg=true; // [flg] Refractive index is user-specified

  // Remove tabular refractive index data, if any, from memory
  rcd+=deallocate(); // [fnc] Free dynamic memory for object

  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::tpt_set()

int // [enm] Return success code
idx_rfr_cls::fl_idx_rfr_set(const std::string &fl_idx_rfr_arg){ // [fnc] File containing refractive indices
  std::string sbr_nm("idx_rfr_cls::fl_idx_rfr_set"); // [sng] Subroutine name

  fl_idx_rfr= (fl_idx_rfr_arg == "") ? fio::data_file_path_get(aer_cls::abb2fl_idx_rfr(abb)) : fl_idx_rfr_arg; // [fnc] File containing refractive indices

  // In case object already has raw refractive index data in dynamic memory
  rcd+=deallocate(); // [fnc] Free dynamic memory for object

  // fxm: Searching for Function is cheesy---should add enumerated type to mnr_sct which indicates for file or function (or other?)
  if(fl_idx_rfr.find("Function") == std::string::npos){
    rcd+=idx_rfr_raw_fl_get // [fnc] Load raw refractive index data from file
      (abb, // I [sng] Composition of particle
       fl_idx_rfr, // I [sng] File containing refractive indices
       wvl_ctr_tbl, // O [m] Wavelength at band center, tabulated
       wvl_nbr_tbl, // O [nbr] Number of wavelength bands, tabulated
       idx_rfr_tbl); // O [frc] Refractive index, tabulated
    
    flg_raw_idx_rfr_in_dyn_mmr=true; // [flg] Raw spectrum from file is in dynamic memory
  } // endif
  
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object

  if(rcd != 0) err_prn(sbr_nm,"Failed raw file get or recompute()"); // [fnc] Print uniform error message and exit

  return rcd; // [enm] Return success code
} // end idx_rfr_cls::fl_idx_rfr_set()

std::complex<prc_cmp> // O [frc] Refractive index
idx_rfr_cls::idx_rfr_get // [fnc] Return refractive index at specified wavelength
(const prc_cmp wvl_ctr) // I [m] Wavelength at band center
  const{
  /* Purpose: Return refractive index at specified wavelength
     Function is overloaded */
  int rcd_lcl(0); // [enm] Return success code
  std::complex<prc_cmp> idx_rfr; // [frc] Refractive index
  const long wvl_nbr_tmp(1L); // [nbr] Number of wavelength bands
  const prc_cmp *wvl_ctr_tmp(&wvl_ctr); // [m] Wavelength at band center

  rcd_lcl+=idx_rfr_get // [fnc] Return array of refractive indices
    (wvl_ctr_tmp, // I [m] Wavelength at band center
     &idx_rfr, // O [frc] Refractive index
     wvl_nbr_tmp); // I [nbr] Number of wavelength bands

  return idx_rfr; // [frc] Refractive index
} // end idx_rfr_get()

int // O [enm] Return success code
idx_rfr_cls::idx_rfr_get // [fnc] Return array of refractive indices
(const prc_cmp *wvl_ctr, // I [m] Wavelength at band center
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const bool wrn_ntp_flg) // I [flg] Print WARNINGs from ntp_vec()
  const{
  /* Purpose: Return refractive indices at specified wavelengths
     Routine works with monotonic and non-monotonic arrays */
  int rcd_lcl(0); // [enm] Return success code

  std::string sbr_nm("idx_rfr_cls::idx_rfr_get"); // [sng] Subroutine name

  if(idx_rfr_usr_flg){
    for(long wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      idx_rfr[wvl_idx]=idx_rfr_usr; // [frc] Refractive index
    } // end loop over wvl 
    return rcd_lcl; // [enm] Return success code
  } // not idx_rfr_usr_flg

  // If idx_rfr is parameterized then get values and return
  // fxm: systematize this to work with pointers to functions
  if(abb == "air"){
    rcd_lcl+=idx_rfr_air_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of air
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif air

  if(abb.find("LQB93") != std::string::npos || abb.find("limestone_dsp_QOL78") != std::string::npos){
    // [fnc] Compute refractive indice for calcite or gypsum aerosol of LQB93 or limestone of QOL78
    rcd_lcl+=idx_rfr_LQB93_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of calcite or gypsum
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif idx_rfr_LQB93

  if(abb == "h2o_lqd" || abb == "h2o_ice"){
    // fxm: overload these functions to have array interface
    for(long wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      if(abb == "h2o_lqd") rcd_lcl+=idx_rfr_H2O_lqd_get(wvl_ctr[wvl_idx],tpt,idx_rfr+wvl_idx);
      if(abb == "h2o_ice") rcd_lcl+=idx_rfr_H2O_ice_get(wvl_ctr[wvl_idx],tpt,idx_rfr+wvl_idx);
    } // end loop over wvl 
    return rcd_lcl; // [enm] Return success code
  } // endif h2o_ice

  if(abb == "lac_ChC90"){
    rcd_lcl+=idx_rfr_ChC90_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of soot
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif lac_ChC90
  if(abb == "lac_DKS91"){
    rcd_lcl+=idx_rfr_DKS91_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of soot
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif lac_DKS91
  if(abb == "lac_WaW80"){
    rcd_lcl+=idx_rfr_WaW80_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of soot
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif lac_WaW80
  if(abb == "lac_YZS00"){
    rcd_lcl+=idx_rfr_YZS00_get
      (abb, // [sng] Composition of particle
       idx_rfr, // [frc] Refractive index of soot
       wvl_ctr, // [m] Wavelength at band center
       wvl_nbr); // [nbr] Number of output wavelength bands
    if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
    return rcd_lcl; // [enm] Return success code
  } // endif lac_YZS00

  /* If requested refractive indices are predicted by a function, they have 
     been returned and function has been exited.
     If function reaches here, assume refractive indices reside in disk file 
     and were already read into memory in fl_idx_rfr_set() */

  // Interpolate disk values to output coordinates
  std::string xtr_sng("xtr_fll_ngh"); // Set extrapolated value to value of nearest valid neighbor
  if(wrn_ntp_flg) xtr_sng+="+xtr_vrb"; // Add verbose diagnostics if requested
  Xtr_cls xtr(xtr_sng); // [enm] Extrapolation flags

  // Interpolate refractive indices from input to output grid
  rcd_lcl+=ntp_vec // [fnc] Interpolate array to requested grid
    (wvl_nbr_tbl, // I [nbr] Input coordinate size
     wvl_ctr_tbl, // I [crd] Input coordinate
     idx_rfr_tbl, // I [frc] Input data (on coordinate)
     wvl_nbr, // I [nbr] Output coordinate size
     wvl_ctr, // I [crd] Output coordinate
     idx_rfr, // O [frc] Output (interpolated) data (on coordinate)
     xtr, // I [enm] LHS extrapolation flags
     xtr); // I [enm] RHS extrapolation flags

  if(dbg_lvl_get() == dbg_vrb){
    std::cout << "Interpolated refractive indices for " << abb << std::endl;
    std::cout << "idx\twvl_ctr\treal\timag" << std::endl;
    std::cout << "\tum\tfrc\tfrc" << std::endl;
    std::cout.setf(std::ios::fixed);
    std::cout.setf(std::ios::showpoint);
    std::cout.setf(std::ios::showpos);
    std::cout.setf(std::ios::internal,std::ios::adjustfield);
    std::cout.setf(std::ios::left);
    std::cout.precision(3);
    for(long wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      // std::cout.setf(std::ios::scientific);
      std::cout << std::setw(8) << " " << wvl_idx 
		<< std::setw(8) << " " << wvl_ctr[wvl_idx]*1.0e6
		<< std::setw(8) << " " << idx_rfr[wvl_idx].real()
		<< std::setw(8) << " " << idx_rfr[wvl_idx].imag()
		<< std::setw(8) << " " << std::endl;
    } // end loop over wvl_ctr
    // Restore defaults
    std::cout.precision(8);
    std::cout.unsetf(std::ios::floatfield | std::ios::showpoint | std::ios::showpos);
  } // end dbg
  
  return rcd_lcl; // O [enm] Return success code

} // end idx_rfr_cls::idx_rfr_get()

// Public member functions end
// Private member functions begin

int // [enm] Return success code
idx_rfr_cls::allocate(){ // [fnc] Allocate dynamic memory for object
  // Purpose: Allocate dynamic memory for object
  // Currently a no-op
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::allocate()

int // [enm] Return success code
idx_rfr_cls::deallocate(){ // [fnc] Free dynamic memory for object
  // Purpose: Free dynamic memory for object

  // Delete any arrays allocated by private member functions
  if(flg_raw_idx_rfr_in_dyn_mmr){
    // netCDF library in idx_rfr_raw_fl_get() creates these arrays
    delete []wvl_ctr_tbl; // [m] Wavelength at band center, tabulated
    delete []idx_rfr_tbl; // [frc] Refractive index, tabulated

    // Reset flag
    flg_raw_idx_rfr_in_dyn_mmr=false; // [flg] Raw refractive indices from file are in dynamic memory
  } // endif

  return rcd; // [enm] Return success code
} // end idx_rfr_cls::deallocate()

int // [enm] Return success code
idx_rfr_cls::recompute(){ // [fnc] Recompute properties of object
  std::string sbr_nm("idx_rfr_cls::recompute()"); // [sng] Subroutine name

  // Set private members which lack set() functions
  dsc=aer_cls::abb2dsc(abb); // [sng] Refractive index description
  idx_rfr_fnc_ptr=aer_cls::abb2fnc(abb); // [fnc] Pointer to function returning refractive indices
  
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::recompute()

int // [enm] Return success code
idx_rfr_cls::idx_rfr_prn // [fnc] Print refractive indices in pretty table
(prc_cmp *wvl_ctr, // I [m] Wavelength at band center
 long wvl_nbr, // I [nbr] Number of wavelength bands
 std::complex<prc_cmp> *idx_rfr) // I [frc] Refractive index
  const{
  // Purpose: Print refractive indices in pretty table
  int rcd_lcl(0); // [enm] Return success code
  long wvl_idx; // [idx] Counting index for wavelength
  
  std::cout << "idx\twvl_ctr\treal\timag" << std::endl;
  // std::cout.setf(ios::scientific);
  std::cout.setf(std::ios::fixed);
  std::cout.setf(std::ios::showpoint);
  std::cout.setf(std::ios::showpos);
  std::cout.setf(std::ios::internal,std::ios::adjustfield);
  std::cout.setf(std::ios::left);
  std::cout.precision(3);
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    std::cout << std::setw(8) << " " << wvl_idx 
	      << std::setw(8) << " " << wvl_ctr[wvl_idx]*1.0e6
	      << std::setw(8) << " " << idx_rfr[wvl_idx].real()
	      << std::setw(8) << " " << idx_rfr[wvl_idx].imag()
	      << std::setw(8) << " " << std::endl;
  } // end loop over wvl_ctr 
    // Restore defaults
  std::cout.precision(8);
  std::cout.unsetf(std::ios::floatfield | std::ios::showpoint | std::ios::showpos);

  return rcd_lcl; // [enm] Return success code
} // end idx_rfr_prn()

int // [enm] Return success code
idx_rfr_cls::idx_rfr_raw_fl_get // [fnc] Load raw refractive index data from file
(const std::string &cmp_sng_prt, // I [sng] Composition of particle
 const std::string &fl_in, // I [sng] File containing refractive indices
 prc_cmp *&wvl_ctr, // O [m] Wavelength at band center, tabulated
 long &wvl_nbr, // O [nbr] Number of wavelength bands, tabulated
 std::complex<prc_cmp> *&idx_rfr) // O [frc] Refractive index, tabulated
{
  // Purpose: Load raw refractive index data from file
  int rcd(0); // [enm] Return success code
  int rcd_opt(NC_ENOTVAR); // [enm] Acceptable return code
  std::string sbr_nm("idx_rfr_raw_fl_get"); // [sng] Subroutine name

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  long wvl_idx; // [idx] Counting index for wavelength
  prc_cmp *idx_rfr_rl_tbl; // [frc] Real refractive index, tabulated
  prc_cmp *idx_rfr_img_tbl; // [frc] Imaginary refractive index, tabulated
  const std::string wvl_nm_gnr("bnd"); // [sng] Name of wavelength coordinate in refractive index file
  std::string wvl_nm; // [sng] Name of wavelength coordinate in refractive index file
  
  // Check for file existance
  struct stat stat_sct;
  rcd+=stat(fl_in.c_str(),&stat_sct);
  if(rcd == -1) std::perror("ERROR std::stat()");

  // Open input file
  int nc_id=nco_open(fl_in,NC_NOWRITE); // [fnc] Open netCDF file
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Opening "+fl_in);
  
  // Name for wavelenth is bnd_cmp_sng for newer files, bnd for older files
  wvl_nm=wvl_nm_gnr+"_"+cmp_sng_prt;
  int var_id; // [id] Variable ID
  rcd=nco_inq_varid(nc_id,wvl_nm,var_id,rcd_opt);
  /* If wavelength variable was not found with mineral-specific name, 
     then use wavelength variable with generic name */
  if(rcd == rcd_opt) wvl_nm=wvl_nm_gnr;

  wvl_nbr=nco_inq_dimlen(nc_id,wvl_nm); // [nbr] Number of wavelength bands
  
  const std::string idx_rfr_rl_nm("idx_rfr_"+cmp_sng_prt+"_rl");
  const std::string idx_rfr_img_nm("idx_rfr_"+cmp_sng_prt+"_img");
  
  // User must free this memory when no longer needed
  // For objects in idx_rfr_cls, this is done in fl_idx_rfr_set() which calls deallocate()
  rcd=nco_get_var(nc_id,wvl_nm,wvl_ctr); // [m] Band nominal wavelength, tabulated
  rcd=nco_get_var(nc_id,idx_rfr_rl_nm,idx_rfr_rl_tbl); // [frc] Real refractive index, tabulated
  rcd=nco_get_var(nc_id,idx_rfr_img_nm,idx_rfr_img_tbl); // [frc] Imaginary refractive index, tabulated
  
  // Check for monotonicity of input grid
  if(!mnt_chk(wvl_ctr,wvl_nbr)) err_prn(sbr_nm,"wvl_ctr not monotonic");

  // Assemble complex values from real and imaginary components
  idx_rfr=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index, tabulated
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(idx_rfr_rl_tbl[wvl_idx],idx_rfr_img_tbl[wvl_idx]); // [frc] Refractive index, tabulated
  } // end for
  delete []idx_rfr_rl_tbl; // [frc] Real refractive index, tabulated
  delete []idx_rfr_img_tbl; // [frc] Imaginary refractive index, tabulated

  std::string wvl_units_sng; // [sng] Value of "units" attribute
  rcd=nco_get_att(nc_id,wvl_nm,static_cast<std::string>("units"),wvl_units_sng);
  if(wvl_units_sng.find("micron") != std::string::npos){
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) wvl_ctr[wvl_idx]/=1.0e6;
  } // end if units are microns
  
  if(dbg_lvl_get() == dbg_vrb) idx_rfr_prn(wvl_ctr,wvl_nbr,idx_rfr);

  // Close input file
  rcd=nco_close(nc_id); // [fnc] Close netCDF file
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Closing "+fl_in);
  
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // end idx_rfr_cls::idx_rfr_raw_fl_get()

// Private member functions end
// Global functions with C++ linkages begin

int // O [enm] Return success code
idx_rfr_air_get // [fnc] Compute refractive indices for air
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of air
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices of air
  // fxm: Default to prs_STP, tpt_STP, but allow user to optionally over-ride
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  std::string sbr_nm("idx_rfr_air_get"); // [sng] Subroutine name
  prc_cmp idx_rfr_rl; // [frc] Real refractive index of air
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt != "air") err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);

  long wvl_idx;
  using phc::tpt_STP; // (273.16) [K] Standard temperature
  using phc::prs_STP; // (101325.0) [Pa] Standard pressure
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    /* Index of refraction of dry air at STP is described by 
       Len93 p. 155, BrS84 p. 107, GoY89 p. 297, Ste94 p. 161, ThS99 p. 73, Lio92 p. 166
       Method of Len93 depends on temperature and pressure, which make sense
       Method of Lio92 has no temperature or pressure dependence
       This compares some methods:
       wvl       Lio92     Len93    
       0.3 um  
       0.5 um  1.00028   1.00029 */

    // fxm: idx_rfr_air depends on local air density, not on STP air density
    idx_rfr_rl= // [frc] Real refractive index of air Len93 p. 155
      1.0+1.0e-6*(77.46+0.459/(1.0e12*wvl[wvl_idx]*wvl[wvl_idx]))*prs_STP*0.01/tpt_STP;
    /*    idx_rfr_rl[wvl_idx]= // [frc] Real refractive index of air Lio92 p. 166 (3.8.34)
      1.0+1.0e-8*(6.43e3+
		  2.95e6/(146.0-1.0/(1.0e12*wvl[wvl_idx]*wvl[wvl_idx]))+
		  2.55e4/(41.0-1.0/(1.0e12*wvl[wvl_idx]*wvl[wvl_idx]))); */
    // Neglect gaseous absorption so that imaginary refractive index of air = 0.0
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(idx_rfr_rl,0.0); // [frc] 
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_air_get()

int // O [enm] Return success code
idx_rfr_ChC90_get // [fnc] Compute refractive indices for soot aerosol of ChC90
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices for soot aerosol of ChC90
  /* ChC90: Flatau recommends this source
     http://atol.ucsd.edu/~pflatau/refrtab/soot
     and provides a spreadsheet of values
     http://atol.ucsd.edu/~pflatau/refrtab/soot/chang1990.xls */
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  prc_cmp wvl_mcr; // [um] Wavelength
  prc_cmp ln1_wvl_mcr; // [um] Log wavelength
  prc_cmp ln2_wvl_mcr; // [um] Log squared wavelength
  prc_cmp ln3_wvl_mcr; // [um] Log cubed wavelength
  prc_cmp idx_rfr_rl; // [frc] Refractive index, real component
  prc_cmp idx_rfr_img; // [frc] Refractive index, imaginary component
  std::string sbr_nm("idx_rfr_ChC90_get"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt != "lac_ChC90") err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);


  long wvl_idx;
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_mcr=1.0e6*wvl[wvl_idx]; // [um] Wavelength
    if(wvl_mcr < 0.19 || wvl_mcr > 30.0) wrn_prn(sbr_nm,"Called "+sbr_nm+" with wavelength = "+nbr2sng(wvl_mcr)+" um. Original range of data for parameterization is 0.2 < lambda < 30 um."+cmp_sng_prt);
    ln1_wvl_mcr=std::log(wvl_mcr); // [um] Log wavelength
    ln2_wvl_mcr=ln1_wvl_mcr*ln1_wvl_mcr; // [um] Log squared wavelength
    ln3_wvl_mcr=ln1_wvl_mcr*ln2_wvl_mcr; // [um] Log cubed wavelength
    idx_rfr_rl=1.811+0.1263*ln1_wvl_mcr+0.027*ln2_wvl_mcr+0.017*ln3_wvl_mcr; // [frc]
    idx_rfr_img=0.5821+0.1213*ln1_wvl_mcr+0.2309*ln2_wvl_mcr-0.01*ln3_wvl_mcr; // [frc]
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(idx_rfr_rl,idx_rfr_img); // [frc]
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ChC90_get()

int // O [enm] Return success code
idx_rfr_DKS91_get // [fnc] Compute refractive indices for soot aerosol of DKS91
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices for soot aerosol of DKS91
  /* DKS91: Assume index of refraction is constant with wavelength 
     Assumption is used by CLV96, who refer to DKS91 as the source */
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  prc_cmp wvl_mcr; // [um] Wavelength
  std::string sbr_nm("idx_rfr_DKS91_get"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt != "lac_DKS91") err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);

  long wvl_idx;
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_mcr=1.0e6*wvl[wvl_idx]; // [um] Wavelength
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(1.75,0.44); // [frc] CLV96 p. 23366
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_DKS91_get()

int // O [enm] Return success code
idx_rfr_LQB93_get // [fnc] Compute refractive indice for calcite or gypsum aerosol of LQB93
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of calcite or gypsum
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices for calcite or gypsum aerosol of LQB93
  /* LQB93 determined refractive indices using surface reflectance techniques
     Normal incident reflectance was measured then fit with a finite series of Lorentz lines.
     A line strength ln_str, damping constant cst_dmp, and line center frequency rsn_ctr
     characterize each line.
     LQB93 tabulates the "modified line strength" ln_str_mod = ln_str/rsn_ctr_sqr
     A high frequency dielectric constant cst_lct sets the non-resonant properties.
     Measurements were made from 30 < nubar < 4000.0 cm-1 by combining multiple instruments.

     Results presented in:
     Long, L. L., M. R. Querry, R. J. Bell, and R. W. Alexander (1993), Optical properties of calcite and gypsum in crystalline and powdered form in the infrared and far-infrared, Infrared Physics 34(2), 191-201.

     Results give no indication that behavior would change drastically outside this range
     I arbitrarily extend the range of validity to 20 < nubar < 5000 cm-1.
     However, this still leaves out the entire visible spectrum.
     Refractive indices must be known to 0.2 um (50000 cm-1) for full solar radiative transfer.
     What to do in this range?
     Sokolik and Toon seem to either extrapolate to the visible, or to use
     amorphous calcite (limestone) properties (from QOL78) in the visible

     20060620: Compared calcite and gypsum indices generated by this routine to LQB93.
     The fit appears perfect and verifies that no bugs are in routines which
     turn LQB93 dispersive analysis parameters into refractive indices.

     For convenience, a single calcite type from QOL78 is stored in this routine
     This compound is called limestone_dsp_QOL78 ("dsp" is for "dispersive")
     Limestone is primarily randomly oriented microcrystals of calcite, with 
     dolomite and organic fossil material as secondary constituents.
     Hence limestone is also known as amorphous calcite.
     QOL78 measurements by:
     Professor Emeritus Marvin Querry, Physics Dept., UMKC, retired in 2003. 
     Querry is pronounced "Query". 
     Querry et al. measured Bethany Falls limestone and QOL78 reports table from 0.2-32 um:
     Querry, M.R, G.C. Osbourn, K. Lies, R. Jordon, R. Coveney, Jr. (1978), Complex refractive index of limestone in the visible and infrared, Appl. Opt., 17(3), 353-356.
     Separate measurements were made to facilitate two analysis techniques, 
     Kramers-Kronig (KK) and Lorentz Dispersion Analysis (DA).
     The refractive indices resulting show "reasonable agreement" from 0.2--25 um.
     QOL78 found DA real and imaginary indices higher than KK in 25--32.8 um.
     Probably due to assumptions in KK analysis.
     DA technique is definitely preferred for lambda > 25 um.

     QOL78 report different dispersion parameters than LQB93:
     QOL78 report line strength rho_j = A_j/4*pi = S_j/(4*pi*nubar*nubar)
     QOL78 report damping parameter gamma_QOL78 = gamma_LQB93/nubar
     LQB93 report line strength A_j = S_j/(nubar*nubar) = 4*pi*rho_j
     LQB93 report damping parameter gamma_LQB93 = nubar*gamma_QOL78
     This routine stores A_j for LQB93 and rho_j for QOL78
     Routine converts both to S_j before applying LQB93 algorithm

     20060620: Compared limestone generated from DA by this routine to QOL78.
     The fit of DA to QOL78 figures is not expected to appears perfect 
     because QOL78 plot only the KK analysis results.
     To the extent intercomparison is possible, the DA reconstruction by this
     routine appears to be accurate.
     The major discrepancy between DA and KK analyses is that DA leads to much
     lower refractive indices (~10^-5 instead of ~10^-2) everywhere except 14 um.
     Also the amplitude of the resonance at 14 um is much larger in DA than KK

     Limestone results may be reported in NOAA 76022603:
     Holland, W.E., M.R. Querry, and R.M. Coveney, Jr. (1975), Measurements of spectral reflectance and optical constants of selected rock samples for application to remote sensing of soil moisture, NOAA Report 76022603: 78 pp.
     which may be another name for
     Holland, W.E., M.R. Querry, and R.M. Coveney, Jr., Measurements of spectral reflectance and optical constants of selected rock samples for application to remote sensing of soil moisture, Completion Report, U.S. Dept. of Commerce grant 04-4-158-27, available from National Technical Information Service, Springfield, VA 22151
     There is also
     Jordan, Ray E., Jr., Querry, Marvin R., Holland, Wayne E., Osbourn, Gordon C., and Coveney, Raymond M., Jr. (1975), Reflectance and complex refractive index of polycrystalline limestone, J. Opt. Soc. Am. A, vol. 65, page 1170. */

  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  const std::complex<prc_cmp> i_img(0,1); // [nbr] Imaginary unit
  const int rsn_prm_nbr(3); // [nbr] Number of parameters per resonance

  const int ln_str_mod_idx(0); // [idx] Lorentz line strength column in resonance table
  const int dmp_idx(1); // [idx] Lorentz damping constant column in resonance table
  const int rsn_ctr_idx(2); // [idx] Resonance center column in resonance table

  const int rsn_nbr_limestone_dsp(4); // [nbr] Number of resonances
  const int rsn_nbr_calcite_oray(5); // [nbr] Number of resonances
  const int rsn_nbr_calcite_eray(3); // [nbr] Number of resonances
  const int rsn_nbr_calcite_pellet(7); // [nbr] Number of resonances
  const int rsn_nbr_gypsum_EX(7); // [nbr] Number of resonances
  const int rsn_nbr_gypsum_pellet(13); // [nbr] Number of resonances
  const int rsn_nbr_max(13); // [nbr] Maximum number of resonances in any mineral

  const int mnr_idx_limestone_dsp(0); // [idx] Position of mineral in resonance table
  const int mnr_idx_calcite_oray(1); // [idx] Position of mineral in resonance table
  const int mnr_idx_calcite_eray(2); // [idx] Position of mineral in resonance table
  const int mnr_idx_calcite_pellet(3); // [idx] Position of mineral in resonance table
  const int mnr_idx_gypsum_EX(4); // [idx] Position of mineral in resonance table
  const int mnr_idx_gypsum_pellet(5); // [idx] Position of mineral in resonance table
  const int mnr_nbr(6); // [nbr] Number of minerals in table

  const prc_cmp lctdlc_hgh_limestone_dsp(2.38); // [F m-1] Dielectric constant high frequency
  const prc_cmp lctdlc_hgh_calcite_oray(2.625); // [F m-1] Dielectric constant high frequency
  const prc_cmp lctdlc_hgh_calcite_eray(2.170); // [F m-1] Dielectric constant high frequency
  const prc_cmp lctdlc_hgh_calcite_pellet(2.169); // [F m-1] Dielectric constant high frequency
  const prc_cmp lctdlc_hgh_gypsum_EX(2.625); // [F m-1] Dielectric constant high frequency
  const prc_cmp lctdlc_hgh_gypsum_pellet(1.975); // [F m-1] Dielectric constant high frequency

  const prc_cmp wvn_max_vld(5000.0); // [cm-1] Maximum wavenumber of applicability
  const prc_cmp wvn_min_vld(20.0); // [cm-1] Minimum wavenumber of applicability

  const prc_cmp rsn_tbl[mnr_nbr][rsn_nbr_max][rsn_prm_nbr]={ // [sct] Resonance parameters
    // ln_str cst_dmp rsn_ctr
    {{0.0320, 0.0521, 1449.0}, // Limestone QOL78 p. 355 Tbl. 1
     {0.0055, 0.0188, 0877.0}, // Band 2
     {0.0054, 0.0070, 0714.0}, // Band 3
     {0.1400, 0.0620, 0303.0}, // Band 4
     {0.0000, 00.000, 0000.0}, // NB: ln_str_QOL78 and cst_dmp_QOL78 
     {0.0000, 00.000, 0000.0}, // differ from corresponding LQB93 definitions
     {0.0000, 00.000, 0000.0}, // Read aer.pdf for details
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0}}, // end Limestone
    // ln_str cst_dmp rsn_ctr
    {{2.7751, 06.037, 0104.1}, // Calcite O-ray LQB93 p. 193 Tbl. 1
     {0.9400, 12.051, 0223.3}, // Band 2
     {1.6657, 13.507, 0295.6}, // Band 3
     {0.0175, 07.517, 0713.2}, // Band 4
     {0.5507, 08.519, 1406.8}, // Band 5
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0}}, // end Calcite O-ray
    // ln_str cst_dmp rsn_ctr
    {{4.6412, 06.797, 0094.6}, // Calcite E-ray LQB93 p. 193 Tbl. 2
     {1.3350, 10.451, 0304.9}, // Band 2
     {0.0817, 01.560, 0870.8}, // Band 3
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0}}, // end Calcite E-ray
    // ln_str cst_dmp rsn_ctr
    {{2.5628, 13.994, 0101.0}, // Calcite pellet LQB93 p. 194 Tbl. 3
     {0.5175, 21.046, 0225.8}, // Band 2
     {1.2617, 27.014, 0302.0}, // Band 3
     {0.0123, 01.600, 0711.0}, // Band 4
     {0.0277, 12.846, 0877.6}, // Band 5
     {0.1739, 41.793, 1420.1}, // Band 6
     {0.1419, 67.328, 1466.2}, // Band 7
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0}}, // end Calcite pellet
    // ln_str cst_dmp rsn_ctr
    {{2.9672, 08.064, 0117.6}, // Gypsum E // X axis LQB93 p. 197 Tbl. 4
     {1.1337, 17.840, 0207.9}, // Band 2
     {2.6027, 97.859, 0210.4}, // Band 3
     {1.0642,128.976, 0333.7}, // Band 4
     {0.4265, 91.616, 0447.6}, // Band 5
     {0.0956, 08.844, 0666.5}, // Band 6
     {0.2421, 18.170, 1135.1}, // Band 7
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0},
     {0.0000, 00.000, 0000.0}}, // end Gypsum E // X axis
    // ln_str cst_dmp rsn_ctr
    {{0.7219, 16.200, 0119.8}, // Gypsum pellet LQB93 p. 197 Tbl. 7
     {0.3389, 20.051, 0174.8}, // Band 2
     {1.4537, 51.256, 0216.8}, // Band 3
     {2.6548, 59.066, 0303.4}, // Band 4
     {0.6069,204.392, 0454.8}, // Band 5
     {0.0566, 19.981, 0599.3}, // Band 6
     {0.0402, 13.672, 0670.5}, // Band 7
     {0.0254,180.041, 0939.0}, // Band 8
     {0.2264, 32.075, 1124.0}, // Band 9
     {0.0061, 13.477, 1621.2}, // Band 10
     {0.0720,751.938, 1707.4}, // Band 11
     {0.0148, 64.124, 3410.2}, // Band 12
     {0.0031, 82.332, 3538.8}} // end Gypsum pellet
  }; // end rsn_tbl

  // Output
  int rcd(0); // [enm] Return success code
  // Local
  bool flg_QOL78(false); // [flg] Convert tabular parameters from QOL78->LQB93 definitions
  int mnr_idx(CEWI_int); // [idx] Position of mineral in resonance table
  int rsn_nbr(CEWI_int); // [nbr] Number of resonances
  prc_cmp lctdlc_hgh(CEWI_cpv); // [F m-1] Dielectric constant high frequency
  prc_cmp ln_str; // [] Lorentz line strength
  prc_cmp cst_dmp; // [wvn] Damping constant
  prc_cmp rsn_ctr_sqr; // [cm-2] Resonance center location, squared
  prc_cmp wvl_mcr; // [um] Wavelength
  prc_cmp wvn; // [cm-1] Wavenumber
  prc_cmp wvn_sqr; // [cm-1] Wavenumber squared
  // 20060605: Following statement breaks xlC and icc 9.0 "error: unsupported underlying vla type (non-POD (Plain Old Data) class type)"
  // std::complex<prc_cmp> idx_rfr_sqr[wvl_nbr]; // [frc] Refractive index of calcite or gypsum
  std::complex<prc_cmp> idx_rfr_sqr; // [frc] Refractive index of calcite or gypsum
  std::string sbr_nm("idx_rfr_LQB93_get"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt.find("LQB93") == std::string::npos && cmp_sng_prt.find("QOL78") == std::string::npos) err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);

  if(cmp_sng_prt.find("calcite_oray") != std::string::npos){
    mnr_idx=mnr_idx_calcite_oray; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_calcite_oray; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_calcite_oray; // [F m-1] Dielectric constant high frequency
  }else if(cmp_sng_prt.find("calcite_eray") != std::string::npos){
    mnr_idx=mnr_idx_calcite_eray; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_calcite_eray; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_calcite_eray; // [F m-1] Dielectric constant high frequency
  }else if(cmp_sng_prt.find("calcite_pellet") != std::string::npos){
    mnr_idx=mnr_idx_calcite_pellet; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_calcite_pellet; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_calcite_pellet; // [F m-1] Dielectric constant high frequency
  }else if(cmp_sng_prt.find("gypsum_EX") != std::string::npos){
    mnr_idx=mnr_idx_gypsum_EX; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_gypsum_EX; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_gypsum_EX; // [F m-1] Dielectric constant high frequency
  }else if(cmp_sng_prt.find("gypsum_pellet") != std::string::npos){
    mnr_idx=mnr_idx_gypsum_pellet; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_gypsum_pellet; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_gypsum_pellet; // [F m-1] Dielectric constant high frequency
  }else if(cmp_sng_prt.find("limestone_dsp") != std::string::npos){
    mnr_idx=mnr_idx_limestone_dsp; // [idx] Position of mineral in resonance table
    rsn_nbr=rsn_nbr_limestone_dsp; // [nbr] Number of resonances
    lctdlc_hgh=lctdlc_hgh_limestone_dsp; // [F m-1] Dielectric constant high frequency
    flg_QOL78=true; // [flg] Convert tabular parameters from QOL78->LQB93 definitions
  } // end else

  for(long wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_mcr=1.0e6*wvl[wvl_idx]; // [um] Wavelength
    wvn=0.01/wvl[wvl_idx]; // [cm-1] Wavenumber
    if(wvn < wvn_min_vld || wvn > wvn_max_vld) 
      if(dbg_lvl_get() > 1)
	wrn_prn(sbr_nm,"Called "+sbr_nm+" with wavelength = "+nbr2sng(wvl_mcr)+" um = "+nbr2sng(wvn)+" cm-1 which is significantly outside bounds of LQB93 measurements, 30 < nubar < 4000.0 cm-1. Will apply Lorentz resonance formula outside the measured domain based on comment of SoT99 p. 9428 that calcite and gypsum have almost no absorption in UV and visible.");
    wvn_sqr=wvn*wvn; // [cm-1] Wavenumber squared
    // Initialize sum that accumulates
    idx_rfr_sqr=lctdlc_hgh; // [F m-1] Electric permittivity
    for(long rsn_idx=0;rsn_idx<rsn_nbr;rsn_idx++){
      rsn_ctr_sqr=rsn_tbl[mnr_idx][rsn_idx][rsn_ctr_idx]*rsn_tbl[mnr_idx][rsn_idx][rsn_ctr_idx];
      if(flg_QOL78){
	// Convert from QOL78 rho_j stored in table to LQB93 S_j used in formula
	ln_str=4.0*mth::cst_M_PIl*rsn_tbl[mnr_idx][rsn_idx][ln_str_mod_idx]*rsn_ctr_sqr; // LQB93 p. 193
	cst_dmp=rsn_tbl[mnr_idx][rsn_idx][rsn_ctr_idx]*rsn_tbl[mnr_idx][rsn_idx][dmp_idx]; // [wvn] Damping constant
      }else{
	// Convert from A_j stored in table to S_j used in formula
	ln_str=rsn_tbl[mnr_idx][rsn_idx][ln_str_mod_idx]*rsn_ctr_sqr; // LQB93 p. 193
	cst_dmp=rsn_tbl[mnr_idx][rsn_idx][dmp_idx]; // [wvn] Damping constant
      } // !flg_QOL78
      idx_rfr_sqr+= // [F m-1] Electric permittivity LQB93 p. 192
	ln_str/(rsn_ctr_sqr-wvn_sqr-i_img*wvn*cst_dmp);
    } // end loop over rsn
    idx_rfr[wvl_idx]=std::sqrt(idx_rfr_sqr); // [frc] Refractive index
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_LQB93_get()

int // O [enm] Return success code
idx_rfr_WaW80_get // [fnc] Compute refractive indices for soot aerosol of WaW80
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices for soot aerosol of WaW80
  /* WaW80 : Personal communication from Mark Flanner
     WaW80 cite Twitty and Weinmann (1971) as source of refractive index = (1.8,0.5)
     Assume for now that index of refraction is constant with wavelength */
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  prc_cmp wvl_mcr; // [um] Wavelength
  std::string sbr_nm("idx_rfr_WaW80_get"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt != "lac_WaW80") err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);

  long wvl_idx;
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_mcr=1.0e6*wvl[wvl_idx]; // [um] Wavelength
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(1.8,0.5); // [frc]
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_WaW80_get()

int // O [enm] Return success code
idx_rfr_YZS00_get // [fnc] Compute refractive indices for soot aerosol of YZS00
(const std::string cmp_sng_prt, // I [sng] Composition of particle
 std::complex<prc_cmp> *idx_rfr, // O [frc] Refractive index of soot
 const prc_cmp *wvl, // I [m] Wavelength
 const long wvl_nbr) // I [nbr] Number of output wavelength bands
{
  // Purpose: Compute and return refractive indices for soot aerosol of YZS00
  /* 
     YZS00 : Personal communication from Shaocai Yu 
     Real and imaginary parts refractive index were inferred at 0.5 um
     or 0.598 um from field measurements made by Saxena's group at NCSU.
     The technique is documented in numerous papers.
     Shaocai then applied an extrapolation formula to obtain the spectral 
     variation of the refractive index throughout the solar spectrum.
     The method of Kent et al. (1983) was used for the real part.
     
     n(wvl)=n(wvl=0.598 um)-0.03*(wvl-0.598)
     
     The real part was calculated by Kent's method: (assume real(0.5)=1.50)
     real(wvl)=1.4971-0.03*(wvl-0.598)
     
     The wavelength dependence of imaginary part was parameterized from
     BC data from DKS91:
     
     img=(2.0e-05*wvl^6-0.0006*wvl^5+0.0085*wvl^4-0.0553*wvl^3+0.1748*wvl^2-0.2004*wvl+0.5104)*0.113966
     
     where wvl=wavelength in microns, and it is assumed that img=0.051 at 0.50 um.
     Thus the aerosol has the same absorption characteristics as black carbon,
     which is a reasonable assumption given the location of the measurements.
  */
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  prc_cmp wvl_mcr; // [um] Wavelength
  prc_cmp idx_rfr_rl; // [frc]
  prc_cmp idx_rfr_img; // [frc]

  std::string sbr_nm("idx_rfr_YZS00_get"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(cmp_sng_prt != "lac_YZS00") err_prn(sbr_nm,"Called "+sbr_nm+" with invalid aerosol type = "+cmp_sng_prt);

  long wvl_idx;
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_mcr=1.0e6*wvl[wvl_idx]; // [um] Wavelength
    idx_rfr_rl=1.4971-0.03*(wvl_mcr-0.598); // [frc]
    idx_rfr_img=(2.0e-05*std::pow(wvl_mcr,PRC_CMP(6.0))-0.0006*std::pow(wvl_mcr,PRC_CMP(5.0))+0.0085*std::pow(wvl_mcr,PRC_CMP(4.0))-0.0553*std::pow(wvl_mcr,PRC_CMP(3.0))+0.1748*std::pow(wvl_mcr,PRC_CMP(2.0))-0.2004*wvl_mcr+0.5104)*0.113966; // [frc]
    idx_rfr[wvl_idx]=std::complex<prc_cmp>(idx_rfr_rl,idx_rfr_img); // [frc]
  } // end loop over wvl
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_YZS00_get()

int // [enm] Return success code
idx_rfr_ffc_get // [fnc] Compute effective optical constants of composites
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const int ffc_mdm_typ, // I [enm] Effective medium type
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction(s) in inclusion(s)
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> *idx_rfr_prt, // I [frc] Refractive index of particle
 std::complex<prc_cmp> *idx_rfr_ffc_brg, // O [frc] Effective refractive index, Bruggeman approximation
 std::complex<prc_cmp> *idx_rfr_ffc_mxg, // O [frc] Effective refractive index, Maxwell Garnett approximation
 std::complex<prc_cmp> *idx_rfr_ffc_pmr, // O [frc] Effective refractive index, partial molar refraction approximation
 std::complex<prc_cmp> *idx_rfr_ffc_vlw, // O [frc] Effective refractive index, volume-weighted approximation
 std::complex<prc_cmp> *idx_rfr_ffc) // O [frc] Effective refractive index of particle
{
  /* Purpose: Compute effective refractive index based on various assumptions
     to provide complete idx_rfr diagnostic output
     Fill idx_rfr_ffc array with refractive indices selected by user
     Routine is used for implementing and testing new effective medium approximations 
     Routine should probably be deprecated in favor of overloaded routine which 
     only computes user-selected effective medium approximation */
  /* Notation: 
     Homogeneous aerosol composition: (Default case)
     Binary aerosol composition:
     1. Case: Coated Sphere
     Core refers to inner sphere, Mantle to outer shell
     2. Case: Effective medium approximations
     2a. Volume-weighted approximation
     Core refers to either substance
     Mantle refers to either substance
     Distinction is not significant since approximation is invariant under exchange of materials (for equal volume fractions)
     2b. Partial molar refraction approximation
     Core refers to non-water substance (e.g., salt)
     Mantle refers to H2O
     Distinction is significant since approximation is not invariant under exchange of materials (not sure, not fully implemented yet)
     2c. Maxwell Garnett approximation
     Matrix refers to matrix composition (usually bulk of particle)
     Inclusion refers to inclusion composition (usually less than bulk)
     Distinction is significant since approximation is not invariant under exchange of materials
     2d. Bruggeman approximation
     Matrix refers to matrix composition (usually bulk of particle)
     Inclusion refers to inclusion composition (usually less than bulk)
     Distinction is not significant since approximation is invariant under exchange of materials */

  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_get"); // [sng] Subroutine name

  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  if(ncl_nbr == 1){
    // Bruggeman approximation
    rcd+=idx_rfr_ffc_brg_get // [fnc] Effective Medium Approximation: Bruggeman
      (wvl_nbr, // I [nbr] Number of wavelengths
       vlm_frc_ncl[0], // I [frc] Volume fraction in inclusion
       idx_rfr_mtx, // I [frc] Refractive index of matrix
       idx_rfr_ncl, // I [frc] Refractive index of inclusion
       idx_rfr_ffc_brg); // O [frc] Effective refractive index, Bruggeman approximation
    
    // Maxwell Garnett approximation
    rcd+=idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett
      (wvl_nbr, // I [nbr] Number of wavelengths
       vlm_frc_ncl[0], // I [frc] Volume fraction in inclusion
       idx_rfr_mtx, // I [frc] Refractive index of matrix
       idx_rfr_ncl, // I [frc] Refractive index of inclusion
       idx_rfr_ffc_mxg); // O [frc] Effective refractive index, Maxwell Garnett approximation

  }else{ // ncl_nbr > 1
    
    // [fnc] Effective Medium Approximation: Maxwell Garnett, multiple inclusions
    const std::complex<prc_cmp> *idx_rfr_ncl_lst[]={idx_rfr_cor,idx_rfr_mnt,idx_rfr_ncl,idx_rfr_prt}; // [frc] Refractive index of inclusion

    rcd+=idx_rfr_ffc_mxg_get
      (ncl_nbr, // I [nbr] Number of inclusions
       wvl_nbr, // I [nbr] Number of wavelengths
       vlm_frc_ncl, // I [frc] Volume fraction in inclusion
       idx_rfr_mtx, // I [frc] Refractive index of matrix
       idx_rfr_ncl_lst, // I [frc] Refractive index of inclusion
       idx_rfr_ffc_mxg); // O [frc] Effective refractive index, Maxwell Garnett approximation

  } // ncl_nbr > 1
  
  // Partial molar refraction approximation
  // fxm: Calling function must provide all CEWI's here
  prc_cmp dns_mtx(CEWI_cpv); // [kg m-3] Density of matrix
  prc_cmp dns_ncl(CEWI_cpv); // [kg m-3] Density of inclusion
  prc_cmp dns_ttl(CEWI_cpv); // [kg m-3] Density of mixture
  prc_cmp mlc_wgt_mtx(CEWI_cpv); // [kg mol-1] Molecular weight of matrix
  prc_cmp mlc_wgt_ncl(CEWI_cpv); // [kg mol-1] Molecular weight of inclusion
  prc_cmp mlc_wgt_ttl(CEWI_cpv); // [kg mol-1] Molecular weight of mixture
  prc_cmp mlr_frc_mtx(CEWI_cpv); // [mol mol-1] Molar fraction of matrix
  prc_cmp mlr_frc_ncl(CEWI_cpv); // [mol mol-1] Molar fraction of inclusion
  prc_cmp mlr_vlm_mtx; // [m3 mol-1] Molar volume of matrix
  prc_cmp mlr_vlm_ncl; // [m3 mol-1] Molar volume of inclusion
  prc_cmp mlr_vlm_ttl; // [m3 mol-1] Molar volume of mixture
  prc_cmp vlm_frc_mtx(CEWI_cpv); // [frc] Volume fraction of matrix
  std::complex<prc_cmp> mlr_rfr_mtx; // [frc] Molar refraction of matrix
  std::complex<prc_cmp> mlr_rfr_ncl; // [frc] Molar refraction of inclusion
  std::complex<prc_cmp> mlr_rfr_ttl; // [frc] Molar refraction of mixture
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    mlr_vlm_mtx=mlc_wgt_mtx/dns_mtx; // [m3 mol-1] Molar volume of matrix
    mlr_vlm_ncl=mlc_wgt_ncl/dns_ncl; // [m3 mol-1] Molar volume of inclusion
    mlr_vlm_ttl=mlr_frc_mtx*mlc_wgt_mtx/dns_mtx+mlr_frc_ncl*mlc_wgt_ncl/dns_ncl; // [m3 mol-1] Molar volume of mixture Ste90
    mlr_rfr_mtx=mlr_vlm_mtx*(idx_rfr_mtx[wvl_idx]*idx_rfr_mtx[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_mtx[wvl_idx]*idx_rfr_mtx[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of matrix Ste90
    mlr_rfr_ncl=mlr_vlm_ncl*(idx_rfr_ncl[wvl_idx]*idx_rfr_ncl[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_ncl[wvl_idx]*idx_rfr_ncl[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of inclusion Ste90
    mlr_frc_mtx=vlm_frc_mtx*(mlc_wgt_ttl/mlc_wgt_mtx)*(dns_mtx/dns_ttl); // [mol mol-1] Molar fraction of matrix
    mlr_frc_ncl=vlm_frc_ncl[0]*(mlc_wgt_ttl/mlc_wgt_ncl)*(dns_ncl/dns_ttl); // [mol mol-1] Molar fraction of inclusion
    mlr_rfr_ttl=mlr_frc_mtx*mlr_rfr_mtx+mlr_frc_ncl*mlr_rfr_ncl; // [frc] Molar refraction of mixture Ste90
    idx_rfr_ffc_pmr[wvl_idx]=std::sqrt((PRC_CMP(1.0)+PRC_CMP(2.0)*mlr_rfr_ttl/mlr_vlm_ttl)/(PRC_CMP(1.0)-mlr_rfr_ttl/mlr_vlm_ttl)); // [frc] Effective refractive index, partial molar refraction approximation Ste90 p. 1677 (5)
  } // end loop over wvl

  // Volume-weighted approximation
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    idx_rfr_ffc_vlw[wvl_idx]=(PRC_CMP(1.0)-vlm_frc_ncl[0])*idx_rfr_mtx[wvl_idx]+vlm_frc_ncl[0]*idx_rfr_ncl[wvl_idx]; // [frc] Effective refractive index, volume-weighted approximation
  } // end loop over wvl

  // Set effective refractive index used in generation of optical properties to user-selected type
  switch(ffc_mdm_typ){
    // Set effective refractive index of composite to appropriate type
  case ffc_mdm_brg:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_ffc_brg[wvl_idx];
    break;
  case ffc_mdm_mxg:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_ffc_mxg[wvl_idx];
    break;
  case ffc_mdm_nil:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_prt[wvl_idx];
    break;
  case ffc_mdm_pmr:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_ffc_pmr[wvl_idx];
    break;
  case ffc_mdm_vlw:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_ffc_vlw[wvl_idx];
    break;
  default:
    err_prn(prg_nm_get(),sbr_nm,"Unknown ffc_mdm_typ");
    break;
  } // end switch 
  
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_get()

int // [enm] Return success code
idx_rfr_ffc_get // [fnc] Return specified Effective Medium Approximation
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const int ffc_mdm_typ, // I [enm] Effective medium type
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction(s) in inclusion(s)
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> *idx_rfr_prt, // I [frc] Refractive index of particle
 std::complex<prc_cmp> *idx_rfr_ffc) // O [frc] Effective refractive index of particle
{
  // Purpose: Compute effective refractive index based on specified Effective Medium Approximation
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_get"); // [sng] Subroutine name
  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  switch(ffc_mdm_typ){
  case ffc_mdm_brg:
    rcd+=idx_rfr_ffc_brg_get // [fnc] Effective Medium Approximation: Bruggeman
      (wvl_nbr, // I [nbr] Number of wavelengths
       vlm_frc_ncl[0], // I [frc] Volume fraction in inclusion
       idx_rfr_mtx, // I [frc] Refractive index of matrix
       idx_rfr_ncl, // I [frc] Refractive index of inclusion
       idx_rfr_ffc); // O [frc] Effective refractive index, Bruggeman approximation
    break;
  case ffc_mdm_mxg:
    if(ncl_nbr == 1){
      
      rcd+=idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett
	(wvl_nbr, // I [nbr] Number of wavelengths
	 vlm_frc_ncl[0], // I [frc] Volume fraction in inclusion
	 idx_rfr_mtx, // I [frc] Refractive index of matrix
	 idx_rfr_ncl, // I [frc] Refractive index of inclusion
	 idx_rfr_ffc); // O [frc] Effective refractive index, Maxwell Garnett approximation
    }else{ // ncl_nbr != 1
      
      // [fnc] Effective Medium Approximation: Maxwell Garnett, multiple inclusions
      const std::complex<prc_cmp> *idx_rfr_ncl_lst[]={idx_rfr_cor,idx_rfr_mnt,idx_rfr_ncl,idx_rfr_prt}; // [frc] Refractive index of inclusion
      
      rcd+=idx_rfr_ffc_mxg_get
	(ncl_nbr, // I [nbr] Number of inclusions
	 wvl_nbr, // I [nbr] Number of wavelengths
	 vlm_frc_ncl, // I [frc] Volume fraction in inclusion
	 idx_rfr_mtx, // I [frc] Refractive index of matrix
	 idx_rfr_ncl_lst, // I [frc] Refractive index of inclusion
	 idx_rfr_ffc); // O [frc] Effective refractive index, Maxwell Garnett approximation

    } // ncl_nbr != 1
    
    break;
  case ffc_mdm_nil:
    // When no effective medium is selected, use "particle" values un-modified
      for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc[wvl_idx]=idx_rfr_prt[wvl_idx];
      break;
  case ffc_mdm_pmr:
    err_prn(sbr_nm,"Effective medium type ffc_mdm_pmr not implemented in "+sbr_nm+" yet");
    break;
  case ffc_mdm_vlw:
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      idx_rfr_ffc[wvl_idx]=(PRC_CMP(1.0)-vlm_frc_ncl[0])*idx_rfr_mtx[wvl_idx]+vlm_frc_ncl[0]*idx_rfr_ncl[wvl_idx]; // [frc] Effective refractive index, volume-weighted approximation
    } // end loop over wvl
    break;
  default:
    err_prn(sbr_nm,"Unknown ffc_mdm_typ");
    break;
  } // end switch 

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_get()

int // [enm] Return success code
idx_rfr_ffc_brg_get // [fnc] Effective Medium Approximation: Bruggeman
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_brg) // O [frc] Effective refractive index, Bruggeman approximation
{
  // Purpose: Compute effective refractive index based on Bruggeman approximation
  /* Distinction between matrix and inclusion is not significant---approximation is invariant under exchange of materials
     Solution has correct limits implemented for vlm_frc_ncl=0.0,1.0
     However, solution jumps significantly as vlm_frc_ncl->0.0
     In other words, behavior for vlm_frc_ncl near 0.0 not so hot
     fxm: investigate vlm_frc_ncl behavrio near 0.0 and improve if possible */
  
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  bool flg_fll_wvl_sln_rqr(true); // [flg] Full wavelength-by-wavelength solution required
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_brg_get"); // [sng] Subroutine name
  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  // Solve special cases first
  if(vlm_frc_ncl == 0.0){
    // Degenerate case occurs when no inclusion
    flg_fll_wvl_sln_rqr=false; // [flg] Full wavelength-by-wavelength solution required
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc_brg[wvl_idx]=idx_rfr_mtx[wvl_idx]; // [frc] Effective dielectric function, Bruggeman approximation
  }else if(vlm_frc_ncl == 1.0){
    // Degenerate case occurs when all inclusion
    flg_fll_wvl_sln_rqr=false; // [flg] Full wavelength-by-wavelength solution required
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++) idx_rfr_ffc_brg[wvl_idx]=idx_rfr_ncl[wvl_idx]; // [frc] Effective dielectric function, Bruggeman approximation
  } // endif

  if(flg_fll_wvl_sln_rqr){
    // Bruggeman approximation
    const std::complex<prc_cmp> prm_a(2.0); // [frc] Parameter a of polynomial equation
    std::complex<prc_cmp> idx_rfr_mtx_sqr;
    std::complex<prc_cmp> idx_rfr_ncl_sqr;
    std::complex<prc_cmp> idx_rfr_ffc_brg_sqr;
    std::complex<prc_cmp> prm_b; // [frc] Parameter b of polynomial equation
    std::complex<prc_cmp> prm_c; // [frc] Parameter c of polynomial equation
    std::complex<prc_cmp> sln_1; // [frc] First root of complex equation
    std::complex<prc_cmp> sln_2; // [frc] Second root of complex equation
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      // Solve wavelength-specific special case first
      if(idx_rfr_mtx[wvl_idx] == idx_rfr_ncl[wvl_idx]){
	// Degenerate case occurs when matrix = inclusion
	idx_rfr_ffc_brg[wvl_idx]=idx_rfr_mtx[wvl_idx]; // [frc] Effective dielectric function, Bruggeman approximation BoH83 p. 217 (8.51)
      }else{ // endif matrix != inclusion
	idx_rfr_mtx_sqr=idx_rfr_mtx[wvl_idx]*idx_rfr_mtx[wvl_idx];
	idx_rfr_ncl_sqr=idx_rfr_ncl[wvl_idx]*idx_rfr_ncl[wvl_idx];
	prm_b=(PRC_CMP(1.0)-PRC_CMP(3.0)*vlm_frc_ncl)*idx_rfr_ncl_sqr+(PRC_CMP(3.0)*vlm_frc_ncl-PRC_CMP(2.0))*idx_rfr_mtx_sqr; // [frc] Parameter b of polynomial equation
	prm_c=-idx_rfr_ncl_sqr*idx_rfr_mtx_sqr; // [frc] Parameter c of polynomial equation
	// Quadratic equation defines effective dielectric function
	eqn_qdr_slvr(prm_a,prm_b,prm_c,&sln_1,&sln_2);
	// Verify root does not violate Bruggeman definition
	assert(sln_1 != PRC_CMP(-0.5)*idx_rfr_mtx_sqr && sln_1 != PRC_CMP(-0.5)*idx_rfr_ncl_sqr);
	assert(sln_2 != PRC_CMP(-0.5)*idx_rfr_mtx_sqr && sln_2 != PRC_CMP(-0.5)*idx_rfr_ncl_sqr);
	idx_rfr_ffc_brg_sqr=(sln_1.imag() > 0.0 ? sln_1 : sln_2); // [frc] Effective dielectric function, Bruggeman approximation BoH83 p. 217 (8.51)
	idx_rfr_ffc_brg[wvl_idx]=std::sqrt(idx_rfr_ffc_brg_sqr); // [frc] Effective refractive index, Bruggeman approximation BoH83 p. 217 (8.51)
      } // endelse
    } // end loop over wvl
  } // endif flg_fll_wvl_sln_rqr
  
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_brg_get()

int // [enm] Return success code
idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> *idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_mxg) // O [frc] Effective refractive index, Maxwell Garnett approximation
{
  // Purpose: Compute effective refractive index based on Maxwell Garnett approximation
  // Distinction between matrix and inclusion is significant since approximation is not invariant under exchange of materials
  
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_mxg_get"); // [sng] Subroutine name
  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Maxwell Garnett approximation
  std::complex<prc_cmp> idx_rfr_mtx_sqr; // [F m-1] Electric permittivity of matrix
  std::complex<prc_cmp> idx_rfr_ncl_sqr; // [F m-1] Electric permittivity of inclusion
  std::complex<prc_cmp> idx_rfr_ffc_sqr; // [F m-1] Effective electric permittivity
  std::complex<prc_cmp> mxg_fct_sum; // [F m-1]
  std::complex<prc_cmp> mxg_fct_dff; // [F m-1]
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    // Factors in Maxwell Garnett approximation ViC98 p. 3 (5), BoH83 p. 217 (8.50)
    idx_rfr_mtx_sqr=idx_rfr_mtx[wvl_idx]*idx_rfr_mtx[wvl_idx]; // [F m-1] Electric permittivity of matrix
    idx_rfr_ncl_sqr=idx_rfr_ncl[wvl_idx]*idx_rfr_ncl[wvl_idx]; // [F m-1] Electric permittivity of inclusion
    mxg_fct_dff=idx_rfr_ncl_sqr-idx_rfr_mtx_sqr; // [F m-1]
    mxg_fct_sum=idx_rfr_ncl_sqr+PRC_CMP(2.0)*idx_rfr_mtx_sqr; // [F m-1]
    idx_rfr_ffc_sqr=idx_rfr_mtx_sqr*(mxg_fct_sum+PRC_CMP(2.0)*vlm_frc_ncl*mxg_fct_dff)/(mxg_fct_sum-vlm_frc_ncl*mxg_fct_dff);
    idx_rfr_ffc_mxg[wvl_idx]=std::sqrt(idx_rfr_ffc_sqr); // [frc] Effective refractive index, Maxwell Garnett approximation
  } // end loop over wvl

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_mxg_get()

int // [enm] Return success code
idx_rfr_ffc_mxg_get // [fnc] Effective Medium Approximation: Maxwell Garnett, multiple inclusions
(const long ncl_nbr, // I [nbr] Number of inclusions
 const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp * const vlm_frc_ncl, // I [frc] Volume fraction in inclusion
 const std::complex<prc_cmp> *idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> **idx_rfr_ncl, // I [frc] Refractive index of inclusion
 std::complex<prc_cmp> *idx_rfr_ffc_mxg) // O [frc] Effective refractive index, Maxwell Garnett approximation
{
  /* Purpose: Compute effective refractive index based on Maxwell Garnett approximation for multiple inclusions
     Distinction between matrix and inclusion is significant since approximation is not invariant under exchange of materials
     Multiple inclusions are treated as in BoH83 p. 216 */
  
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long ncl_idx; // [idx] Counting index for inclusion
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_mxg_get"); // [sng] Subroutine name
  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Maxwell Garnett approximation
  std::complex<prc_cmp> idx_rfr_mtx_sqr; // [F m-1] Electric permittivity of matrix
  std::complex<prc_cmp> idx_rfr_ncl_sqr; // [F m-1] Electric permittivity of inclusion
  std::complex<prc_cmp> idx_rfr_ffc_sqr; // [F m-1] Effective electric permittivity
  std::complex<prc_cmp> mxg_fct_sum; // [F m-1]
  std::complex<prc_cmp> mxg_fct_dff; // [F m-1]
  std::complex<prc_cmp> smm_dnm; // [frc] Denominator summation of f_j*beta_j
  std::complex<prc_cmp> smm_nmr; // [frc] Numerator summation of f_j*beta_j*epsilon_j
  std::complex<prc_cmp> beta_shape; // [frc] Shape factor beta BoH83 p. 216 (8.48)
  std::complex<prc_cmp> vlm_frc_beta_shape; // [frc] Product f_j*beta_j

  prc_cmp vlm_frc_ncl_ttl(0.0); // [frc] Total volume fraction of inclusions
  for(ncl_idx=0;ncl_idx<ncl_nbr;ncl_idx++){
    vlm_frc_ncl_ttl+=vlm_frc_ncl[ncl_idx]; // [frc] Total volume fraction of inclusions
  } // end loop over ncl
  const prc_cmp vlm_frc_mtx(1.0-vlm_frc_ncl_ttl); // [frc] Volume fraction of matrix
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    // Initialize incremented variables
    smm_dnm=0.0; // [frc] Denominator summation of f_j*beta_j
    smm_nmr=0.0; // [frc] Numerator summation of f_j*beta_j*epsilon_j
    for(ncl_idx=0;ncl_idx<ncl_nbr;ncl_idx++){
      // Factors in Maxwell Garnett approximation ViC98 p. 3 (5), BoH83 p. 217 (8.50)
      idx_rfr_mtx_sqr=idx_rfr_mtx[wvl_idx]*idx_rfr_mtx[wvl_idx]; // [F m-1] Electric permittivity of matrix
      idx_rfr_ncl_sqr=idx_rfr_ncl[ncl_idx][wvl_idx]*idx_rfr_ncl[ncl_idx][wvl_idx]; // [F m-1] Electric permittivity of inclusion
      mxg_fct_dff=idx_rfr_ncl_sqr-idx_rfr_mtx_sqr; // [F m-1]
      mxg_fct_sum=idx_rfr_ncl_sqr+PRC_CMP(2.0)*idx_rfr_mtx_sqr; // [F m-1]
      beta_shape=PRC_CMP(3.0)*idx_rfr_mtx_sqr/mxg_fct_sum; // [frc] Shape factor beta BoH83 p. 216 (8.48)
      vlm_frc_beta_shape=vlm_frc_ncl[ncl_idx]*beta_shape; // [frc] Product f_j*beta_j
      smm_dnm+=vlm_frc_beta_shape; // [frc] Denominator summation of f_j*beta_j
      smm_nmr+=vlm_frc_beta_shape*idx_rfr_ncl_sqr; // [frc] Numerator summation of f_j*beta_j*epsilon_j
    } // end loop over ncl
    idx_rfr_ffc_sqr=(vlm_frc_mtx*idx_rfr_mtx_sqr+smm_nmr)/(vlm_frc_mtx+smm_dnm);
    idx_rfr_ffc_mxg[wvl_idx]=std::sqrt(idx_rfr_ffc_sqr); // [frc] Effective refractive index, Maxwell Garnett approximation
  } // end loop over wvl

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_mxg_get()

int // [enm] Return success code
idx_rfr_ffc_pmr_get // [fnc] Compute refractive indices using partial molar refraction approximation
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const prc_cmp dns_cor, // [kg m-3] Density of core
 const prc_cmp dns_mnt, // [kg m-3] Density of mantle
 const prc_cmp dns_ttl, // [kg m-3] Density of mixture
 const prc_cmp mlc_wgt_cor, // [kg mol-1] Molecular weight of core
 const prc_cmp mlc_wgt_mnt, // [kg mol-1] Molecular weight of mantle
 const prc_cmp mlc_wgt_ttl, // [kg mol-1] Molecular weight of mixture
 const prc_cmp vlm_frc_cor, // [frc] Volume fraction in core
 const prc_cmp vlm_frc_mnt, // I [frc] Volume fraction in mantle
 const std::complex<prc_cmp> *idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> *idx_rfr_mnt, // I [frc] Refractive index of mantle
 std::complex<prc_cmp> *idx_rfr_ffc_pmr) // O [frc] Effective refractive index, partial molar refraction approximation
{
  // Purpose: Compute effective refractive index based on partial molar refraction approximation
  /* Notation: 
     2b. Partial molar refraction approximation
     Core refers to non-water substance (e.g., salt)
     Mantle refers to H2O
     Distinction is not significant since approximation is invariant under exchange of materials
  */

  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wavelength
  std::string sbr_nm("idx_rfr_ffc_pmr_get"); // [sng] Subroutine name
  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Partial molar refraction approximation
  // fxm: Calling function must provide all CEWI's here
  prc_cmp mlr_frc_cor(CEWI_cpv); // [mol mol-1] Molar fraction of core
  prc_cmp mlr_frc_mnt(CEWI_cpv); // [mol mol-1] Molar fraction of mantle
  prc_cmp mlr_vlm_cor; // [m3 mol-1] Molar volume of core
  prc_cmp mlr_vlm_mnt; // [m3 mol-1] Molar volume of mantle
  prc_cmp mlr_vlm_ttl; // [m3 mol-1] Molar volume of mixture
  std::complex<prc_cmp> mlr_rfr_cor; // [frc] Molar refraction of core
  std::complex<prc_cmp> mlr_rfr_mnt; // [frc] Molar refraction of mantle
  std::complex<prc_cmp> mlr_rfr_ttl; // [frc] Molar refraction of mixture
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    mlr_vlm_cor=mlc_wgt_cor/dns_cor; // [m3 mol-1] Molar volume of core
    mlr_vlm_mnt=mlc_wgt_mnt/dns_mnt; // [m3 mol-1] Molar volume of mantle
    mlr_vlm_ttl=mlr_frc_cor*mlc_wgt_cor/dns_cor+mlr_frc_mnt*mlc_wgt_mnt/dns_mnt; // [m3 mol-1] Molar volume of mixture Ste90
    mlr_rfr_cor=mlr_vlm_cor*(idx_rfr_cor[wvl_idx]*idx_rfr_cor[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_cor[wvl_idx]*idx_rfr_cor[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of core Ste90
    mlr_rfr_mnt=mlr_vlm_mnt*(idx_rfr_mnt[wvl_idx]*idx_rfr_mnt[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_mnt[wvl_idx]*idx_rfr_mnt[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of mantle Ste90
    mlr_frc_cor=vlm_frc_cor*(mlc_wgt_ttl/mlc_wgt_cor)*(dns_cor/dns_ttl); // [mol mol-1] Molar fraction of core
    mlr_frc_mnt=vlm_frc_mnt*(mlc_wgt_ttl/mlc_wgt_mnt)*(dns_mnt/dns_ttl); // [mol mol-1] Molar fraction of mantle
    mlr_rfr_ttl=mlr_frc_cor*mlr_rfr_cor+mlr_frc_mnt*mlr_rfr_mnt; // [frc] Molar refraction of mixture Ste90
    idx_rfr_ffc_pmr[wvl_idx]=std::sqrt((PRC_CMP(1.0)+PRC_CMP(2.0)*mlr_rfr_ttl/mlr_vlm_ttl)/(PRC_CMP(1.0)-mlr_rfr_ttl/mlr_vlm_ttl)); // [frc] Effective refractive index, partial molar refraction approximation Ste90 p. 1677 (5)
  } // end loop over wvl

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end idx_rfr_ffc_pmr_get()

int // [enm] Return success code
idx_rfr_ffc_pmr_lst_get // [fnc] Compute refractive indices using partial molar refraction approximation
(const long wvl_nbr, // I [nbr] Number of wavelengths
 const aer_cls *aer_lst, // I [aer] List of aerosol components
 const long &aer_nbr, // I [frc] Number of aerosol components
 std::complex<prc_cmp> *idx_rfr_ffc_pmr) // O [frc] Effective refractive index, partial molar refraction approximation
{
  /* Purpose: Compute effective refractive index based on partial molar refraction approximation
     fxm: finish variables routine requires so that it accepts arbitrary lists
     of aerosols and computes partial molar refraction approximation on result */
  // Output
  int rcd(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wavelength
  long aer_idx; // [idx] Counting index for aerosol mode
  const prc_cmp vlm_frc_mnt(0.0); // I [frc] Volume fraction in mantle
  std::string sbr_nm("idx_rfr_ffc_pmr_lst_get"); // [sng] Subroutine name

  // Main code
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Sanity check
  if(aer_nbr <= 0) err_prn("aer_cls::mss_frc_anl_set","aer_nbr <= 0");

  // Partial molar refraction approximation
  prc_cmp dns_cor(CEWI_cpv); // [kg m-3] Density of core
  prc_cmp dns_mnt(CEWI_cpv); // [kg m-3] Density of mantle
  prc_cmp dns_ttl(CEWI_cpv); // [kg m-3] Density of mixture
  prc_cmp mlc_wgt_cor(CEWI_cpv); // [kg mol-1] Molecular weight of core
  prc_cmp mlc_wgt_mnt(CEWI_cpv); // [kg mol-1] Molecular weight of mantle
  prc_cmp mlc_wgt_ttl(CEWI_cpv); // [kg mol-1] Molecular weight of mixture
  prc_cmp mlr_frc_cor(CEWI_cpv); // [mol mol-1] Molar fraction of core
  prc_cmp mlr_frc_mnt(CEWI_cpv); // [mol mol-1] Molar fraction of mantle
  prc_cmp mlr_vlm_cor; // [m3 mol-1] Molar volume of core
  prc_cmp mlr_vlm_mnt; // [m3 mol-1] Molar volume of mantle
  prc_cmp mlr_vlm_ttl; // [m3 mol-1] Molar volume of mixture
  prc_cmp vlm_frc_cor(CEWI_cpv); // [frc] Volume fraction in core
  std::complex<prc_cmp> mlr_rfr_cor; // [frc] Molar refraction of core
  std::complex<prc_cmp> mlr_rfr_mnt; // [frc] Molar refraction of mantle
  std::complex<prc_cmp> mlr_rfr_ttl; // [frc] Molar refraction of mixture

  // fxm: initialize these before using them
  std::complex<prc_cmp> *idx_rfr_cor=new std::complex<prc_cmp>[wvl_nbr]; // I [frc] Refractive index of core
  std::complex<prc_cmp> *idx_rfr_mnt=new std::complex<prc_cmp>[wvl_nbr]; // I [frc] Refractive index of mantle
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    idx_rfr_cor[wvl_idx]=std::complex<prc_cmp>(0.0,0.0); // [frc] Refractive index of core
    idx_rfr_mnt[wvl_idx]=std::complex<prc_cmp>(0.0,0.0); // [frc] Refractive index of mantle
 } // end loop over wvl

  for(aer_idx=0;aer_idx<aer_nbr;aer_idx++){
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      mlr_vlm_cor=mlc_wgt_cor/dns_cor; // [m3 mol-1] Molar volume of core
      mlr_vlm_mnt=mlc_wgt_mnt/dns_mnt; // [m3 mol-1] Molar volume of mantle
      mlr_vlm_ttl=mlr_frc_cor*mlc_wgt_cor/dns_cor+mlr_frc_mnt*mlc_wgt_mnt/dns_mnt; // [m3 mol-1] Molar volume of mixture Ste90
      mlr_rfr_cor=mlr_vlm_cor*(idx_rfr_cor[wvl_idx]*idx_rfr_cor[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_cor[wvl_idx]*idx_rfr_cor[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of core Ste90
      mlr_rfr_mnt=mlr_vlm_mnt*(idx_rfr_mnt[wvl_idx]*idx_rfr_mnt[wvl_idx]-PRC_CMP(1.0))/(idx_rfr_mnt[wvl_idx]*idx_rfr_mnt[wvl_idx]+PRC_CMP(2.0)); // [frc] Molar refraction of mantle Ste90
      mlr_frc_cor=vlm_frc_cor*(mlc_wgt_ttl/mlc_wgt_cor)*(dns_cor/dns_ttl); // [mol mol-1] Molar fraction of core
      mlr_frc_mnt=vlm_frc_mnt*(mlc_wgt_ttl/mlc_wgt_mnt)*(dns_mnt/dns_ttl); // [mol mol-1] Molar fraction of mantle
      mlr_rfr_ttl=mlr_frc_cor*mlr_rfr_cor+mlr_frc_mnt*mlr_rfr_mnt; // [frc] Molar refraction of mixture Ste90
      idx_rfr_ffc_pmr[wvl_idx]=std::sqrt((PRC_CMP(1.0)+PRC_CMP(2.0)*mlr_rfr_ttl/mlr_vlm_ttl)/(PRC_CMP(1.0)-mlr_rfr_ttl/mlr_vlm_ttl)); // [frc] Effective refractive index, partial molar refraction approximation Ste90 p. 1677 (5)
    } // end loop over wvl
  } // end loop over aer
  
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");

  if(&aer_lst == &aer_lst){;} // CEWU Compiler Error Warning Usage

  return rcd;
} // end idx_rfr_ffc_pmr_lst_get()

// Global functions with C++ linkages end
