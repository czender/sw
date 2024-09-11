// $Id$ 

// Implementation (declaration) of solar spectra classes

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <spc_slr.hh> // Solar spectra

// spc_slr_cls class

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const spc_slr_cls &spc_slr_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for spc_slr_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     std::cout << spc_slr_obj; */
  srm_out << "Public contents of spc_slr_obj: " << std::endl;
  srm_out << "spc_slr_obj.nst_nbr_get() = " << spc_slr_obj.nst_nbr_get() << std::endl;
  srm_out << "spc_slr_obj.typ_get() = " << spc_slr_obj.typ_get() << std::endl;
  srm_out << "spc_slr_obj.dsc_get() = " << spc_slr_obj.dsc_get() << std::endl;
  srm_out << "spc_slr_obj.fl_slr_spc_get() = " << spc_slr_obj.fl_slr_spc_get() << std::endl;
  srm_out << "spc_slr_obj.flx_slr_frc_fnc_ptr_get() = " << spc_slr_obj.flx_slr_frc_fnc_ptr_get() << std::endl;
  srm_out << "spc_slr_obj.slr_cst_get() = " << spc_slr_obj.slr_cst_get() << std::endl;
  srm_out << "Private contents of spc_slr_obj: " << std::endl;
  srm_out << "wvl_nbr_in = " << spc_slr_obj.wvl_nbr_in << std::endl;
  if(spc_slr_obj.flg_raw_spc_in_dyn_mmr){
    srm_out << "wvl_min_in[0] = " << spc_slr_obj.wvl_min_in[0] << std::endl;
    srm_out << "wvl_max_in[wvl_nbr_in-1] = " << spc_slr_obj.wvl_max_in[spc_slr_obj.wvl_nbr_in-1] << std::endl;
    srm_out << "flx_slr_frc_in[wvl_nbr_in-1] = " << spc_slr_obj.flx_slr_frc_in[spc_slr_obj.wvl_nbr_in-1] << std::endl;
    srm_out << "flx_frc_blr_in[wvl_nbr_in-1] = " << spc_slr_obj.flx_frc_blr_in[spc_slr_obj.wvl_nbr_in-1] << std::endl;
  }else{
    srm_out << "No raw spectrum in dynamic memory" << std::endl;
  } // endif
    
  return srm_out; // [srm] Reference to output stream for cascading
} // !operator<<()

// Friendly functions end
// Static members begin

int spc_slr_cls::nst_nbr=0; // [nbr] Number of instantiated class members
// 20240911 Change default solar spectrum, constant from LaN68, 1367 W/m2 to FDE24, 1361.353 W/m2 from CMIP7
//prc_cmp spc_slr_cls::slr_cst_dfl=1367.0; // [W m-2] Default solar constant (CCM3)
//std::string spc_slr_cls::spc_slr_typ_dfl="LaN68"; // [sng] Default solar flux source abbreviation
prc_cmp spc_slr_cls::slr_cst_dfl=1361.353; // [W m-2] Default solar constant (CCM3)
std::string spc_slr_cls::spc_slr_typ_dfl="FDE24"; // [sng] Default solar flux source abbreviation

// Static members end
// Static member functions begin

int // O [enm] Return success code
spc_slr_cls::tst(long obj_nbr){ // [fnc] Self-test of spc_slr_cls class
  // Purpose: Perform self test of spc_slr_cls class
  int rcd_lcl(0); // [enm] Return success code
  std::string sbr_nm("spc_slr_cls::tst"); // [sng] Subroutine name

  // Test for memory leaks
  std::cout << sbr_nm+"()"+" is testing for memory leaks by creating and destroying " << obj_nbr << " spc_slr_cls objects..." << std::endl;
  long idx; // [idx] Counting index
  spc_slr_cls *tst_obj; // [sct] Test object
  prc_cmp flx_frc; // [fnc] Fraction of solar spectrum in a single spectral region
  for(idx=0;idx<obj_nbr;idx++){
    tst_obj=new spc_slr_cls; // [sct] Test object
    flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
    std::cout << "idx = " << idx << ", nst_nbr = " << nst_nbr << ", flx_frc = " << flx_frc << std::endl;
    delete tst_obj; // [sct] Test object
  } // !loop over obj

  // Test all solar spectra
  tst_obj=new spc_slr_cls("LaN68"); // [sct] Test object
  flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  delete tst_obj; // [sct] Test object

  tst_obj=new spc_slr_cls("ThD71"); // [sct] Test object
  flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  delete tst_obj; // [sct] Test object

  tst_obj=new spc_slr_cls("Kur95_20wvn"); // [sct] Test object
  flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  delete tst_obj; // [sct] Test object

  tst_obj=new spc_slr_cls("Kur95_01wvn"); // [sct] Test object
  flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  delete tst_obj; // [sct] Test object

  tst_obj=new spc_slr_cls("FDE24"); // [sct] Test object
  flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  delete tst_obj; // [sct] Test object

  // 20240911: Testing lsr source breaks because there is no openable file with its spectrum
  //tst_obj=new spc_slr_cls("lsr"); // [sct] Test object
  //flx_frc=tst_obj->flx_frc_get(0.2e-6,1.0e-5); // [fnc] Fraction of solar spectrum in a single spectral region
  //std::cout << "abb = " << tst_obj->abb << ", flx_frc = " << flx_frc << std::endl;
  //delete tst_obj; // [sct] Test object

  return rcd_lcl; // [enm] Return success code
} // !spc_slr_cls::tst()

int spc_slr_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members
sng2sng_map spc_slr_cls::opt2abb_map=spc_slr_cls::opt2abb_map_mk(); // [map] Option to abbreviation map
sng2sng_map spc_slr_cls::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp; // [map] Option to abbreviation map

  map_tmp.insert(sng2sng_map::value_type("ThD71","ThD71"));
  map_tmp.insert(sng2sng_map::value_type("Thekeakara","ThD71"));
  map_tmp.insert(sng2sng_map::value_type("thekeakara","ThD71"));

  map_tmp.insert(sng2sng_map::value_type("LaN68","LaN68"));
  map_tmp.insert(sng2sng_map::value_type("Labs","LaN68"));
  map_tmp.insert(sng2sng_map::value_type("labs","LaN68"));

  map_tmp.insert(sng2sng_map::value_type("NeL84","NeL84"));
  map_tmp.insert(sng2sng_map::value_type("Neckel","NeL84"));
  map_tmp.insert(sng2sng_map::value_type("neckel","NeL84"));

  map_tmp.insert(sng2sng_map::value_type("Kur95_01wvn","Kur95_01wvn"));
  map_tmp.insert(sng2sng_map::value_type("Kur95","Kur95_01wvn"));
  map_tmp.insert(sng2sng_map::value_type("kurucz","Kur95_01wvn"));
  map_tmp.insert(sng2sng_map::value_type("Kurucz","Kur95_01wvn"));

  map_tmp.insert(sng2sng_map::value_type("Kur95_20wvn","Kur95_20wvn"));

  map_tmp.insert(sng2sng_map::value_type("JHC21","JHC21"));
  map_tmp.insert(sng2sng_map::value_type("Jing","JHC21"));
  map_tmp.insert(sng2sng_map::value_type("Huang","JHC21"));

  map_tmp.insert(sng2sng_map::value_type("FDE24","FDE24"));
  map_tmp.insert(sng2sng_map::value_type("Funke","FDE24"));
  map_tmp.insert(sng2sng_map::value_type("CMIP7","FDE24"));

  map_tmp.insert(sng2sng_map::value_type("lsr","lsr"));
  map_tmp.insert(sng2sng_map::value_type("laser","lsr"));
  map_tmp.insert(sng2sng_map::value_type("Laser","lsr"));

  return map_tmp;
} // !spc_slr_cls::opt2abb_map_mk()

std::string // O [sng] Solar flux source abbreviation
spc_slr_cls::opt2abb // [fnc] Option to abbreviation mapper
(const std::string &opt_sng){ // [sng] Shorthand option
  // Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("spc_slr_cls::opt2abb",opt_sng+" is unknown");
  return itr->second; // O [sng] Solar flux source abbreviation
} // !spc_slr_cls::opt2abb()

sng2spc_slr_sct_map spc_slr_cls::spc_slr_map=spc_slr_cls::spc_slr_map_mk(); // [map] Solar flux source map
sng2spc_slr_sct_map spc_slr_cls::spc_slr_map_mk(){ // [fnc] Create solar flux source map
  const spc_slr_sct spc_slr[]={
    {"ThD71", // [sng] Solar flux source abbreviation
     "Thekeakara and Drummond (1971)", // [sng] Solar flux source description
     fio::data_file_path_get("spc_ThD71.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"LaN68", // [sng] Solar flux source abbreviation
     "Labs and Neckel (1968)", // [sng] Solar flux source description
     fio::data_file_path_get("spc_LaN68.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"NeL84", // [sng] Solar flux source abbreviation
     "Neckel and Labs (1984)", // [sng] Solar flux source description
     fio::data_file_path_get("spc_NeL84.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"Kur95_20wvn", // [sng] Solar flux source abbreviation
     "Kurucz (1995) 20 cm-1", // [sng] Solar flux source description
     fio::data_file_path_get("spc_Kur95_20wvn.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"Kur95_01wvn", // [sng] Solar flux source abbreviation
     "Kurucz (1995) 01 cm-1", // [sng] Solar flux source description
     fio::data_file_path_get("spc_Kur95_01wvn.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"JHC21", // [sng] Solar flux source abbreviation
     "Jing, Huang, Chen, et al. (2021)", // [sng] Solar flux source description
     fio::data_file_path_get("spc_JHC21.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"FDE24", // [sng] Solar flux source abbreviation
     "Funke, Dudok de Wit, Ermolli, et al. (2024) CMIP7", // [sng] Solar flux source description
     fio::data_file_path_get("spc_FDE24.nc"), // [sng] File containing solar spectrum
     (flx_slr_frc_fnc_ptr_typ)CEWI_NULL}, // [fnc] Function to compute solar spectrum
    {"lsr", // [sng] Solar flux source abbreviation
     "Laser (delta function)", // [sng] Solar flux source description
     fio::data_file_path_get("spc_lsr_foo.nc"), // [sng] File containing solar spectrum
     flx_slr_frc_lsr} // [fnc] Function to compute solar spectrum
  }; // !spc_slr_sct spc_slr[]
  long idx; // [idx] Counting index
  int spc_slr_nbr=sizeof(spc_slr)/sizeof(spc_slr_sct); // [nbr] Number of solar flux source structures
  sng2spc_slr_sct_map spc_slr_map_tmp; // [sct] Map with key=abbreviation, value=solar flux source structure
  for(idx=0;idx<spc_slr_nbr;idx++){
    /* fxm: Define variables before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map */
    spc_slr_map_tmp.insert(sng2spc_slr_sct_map::value_type(spc_slr[idx].abb,spc_slr[idx])); // [sct] Map with key=abbreviation, value=solar flux source structure
  } // !idx
  // Return a value to initialize a static class member
  return spc_slr_map_tmp;
} // !spc_slr_cls::spc_slr_map_mk()

std::string // [sng] Description string
spc_slr_cls::abb2dsc(const std::string &abb_sng){ // [fnc] Abbreviation to description mapper
  // Purpose: Return description of solar flux source
  return spc_slr_map.find(abb_sng)->second.dsc; // [sng] Description string
} // !spc_slr_cls::abb2dsc()

flx_slr_frc_fnc_ptr_typ // [fnc] Pointer to function returning fractional solar fluxes
spc_slr_cls::abb2fnc // [fnc] Abbreviation to function mapper
(const std::string &abb_sng){ // [sng] Abbreviation string
  // Purpose: Return pointer to function which computes solar flux source
  return spc_slr_map.find(abb_sng)->second.flx_slr_frc_fnc_ptr; // [fnc] Pointer to function returning fractional solar fluxes
} // !spc_slr_cls::abb2fnc()

std::string // [sng] File containing solar spectrum
spc_slr_cls::abb2fl_slr_spc(const std::string &abb_sng){ // [fnc] Abbreviation to file mapper
  // Purpose: Return file containing solar spectrum
  return spc_slr_map.find(abb_sng)->second.fl_slr_spc; // [sng] File containing solar spectrum
} // !spc_slr_cls::abb2fl_slr_spc()

// Static member functions end
// Public member functions begin

spc_slr_cls::spc_slr_cls // [fnc] Constructor
(const std::string &spc_slr_arg, // [sng] Solar spectrum string
 const prc_cmp &slr_cst_arg, // [W m-2] Solar constant
 const std::string &fl_slr_spc_arg) // [sng] File containing solar spectrum
{
  // Purpose: Constructor
  rcd=0; // [enm] Return success code
  flg_raw_spc_in_dyn_mmr=false; // [flg] Raw spectrum is in dynamic memory

  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls

  // Set type before solar constant to prevent warnings for powerful lasers
  rcd+=typ_set(spc_slr_arg); // [sng] Set solar flux source

  // Constructor defaults to reasonable value to allow for calls that omit slr_cst
  rcd+=slr_cst_set(slr_cst_arg); // [W m-2] Solar constant

  // Load raw solar spectrum data into memory
  rcd+=fl_slr_spc_set(fl_slr_spc_arg); // [sng] File containing solar spectrum

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  nst_nbr++; // [nbr] Number of instantiated class members
} // !spc_slr_cls constructor

spc_slr_cls::~spc_slr_cls(){ // Destructor
  rcd+=deallocate(); // [fnc] Free dynamic memory for object
  nst_nbr--; // [nbr] Number of instantiated class members
} // !spc_slr_cls destructor

void spc_slr_cls::prn()const{std::cout << this;} // [fnc] Print object contents

std::string spc_slr_cls::dsc_get()const{return dsc;} // [sng] Solar flux source description
std::string spc_slr_cls::fl_slr_spc_get()const{return fl_slr_spc;} // [sng] File containing solar spectrum
flx_slr_frc_fnc_ptr_typ spc_slr_cls::flx_slr_frc_fnc_ptr_get()const{return flx_slr_frc_fnc_ptr;} // [fnc] Function to compute solar spectrum
std::string spc_slr_cls::typ_get()const{return abb;} // [sng] Solar flux source abbreviation
prc_cmp spc_slr_cls::slr_cst_get()const{return slr_cst;} // [W m-2] Solar constant

int // [enm] Return success code
spc_slr_cls::fl_slr_spc_set(const std::string &fl_slr_spc_arg){ // [fnc] File containing solar spectrum
  fl_slr_spc= (fl_slr_spc_arg == "") ? abb2fl_slr_spc(abb) : fl_slr_spc_arg; // [fnc] File containing solar spectrum
  fl_slr_spc=drc_pfx(fio::data_path_get(),fl_slr_spc); // [fnc] File containing solar spectrum

  // In case object already has raw spectrum data in dynamic memory
  rcd+=deallocate(); // [fnc] Free dynamic memory for object

  if(fl_slr_spc != ""){
    rcd+=spc_slr_raw_fl_get // [fnc] Load raw solar spectrum data from file
      (fl_slr_spc, // I [sng] File containing solar spectrum
       wvl_min_in, // O [m] Minimum wavelength
       wvl_max_in, // O [m] Maximum wavelength
       wvl_nbr_in, // O [nbr] Number of wavelength bands
       flx_frc_blr_in, // O [frc] Fraction of solar flux at shorter wavelengths
       flx_slr_frc_in); // O [frc] Fraction of solar flux in band
    
    flg_raw_spc_in_dyn_mmr=true; // [flg] Raw spectrum from file is in dynamic memory
  } // endif
  
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !spc_slr_cls::fl_slr_spc_set()

int // [enm] Return success code
spc_slr_cls::flx_slr_frc_fnc_ptr_set(const flx_slr_frc_fnc_ptr_typ &flx_slr_frc_fnc_ptr_arg){ // [fnc] Function to compute solar spectrum
  std::string sbr_nm("spc_slr_cls::flx_slr_frc_fnc_ptr_set"); // [sng] Subroutine name
  flx_slr_frc_fnc_ptr=flx_slr_frc_fnc_ptr_arg; // [fnc] Function to compute solar spectrum
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !spc_slr_cls::flx_slr_frc_fnc_ptr_set()

int // [enm] Return success code
spc_slr_cls::typ_set(const std::string &spc_slr_arg){ // [fnc] Set solar flux source
  std::string sbr_nm("spc_slr_cls::typ_set"); // [sng] Subroutine name
  abb= (spc_slr_arg == "") ? spc_slr_typ_dfl : spc_slr_arg; // [sng] Solar flux source abbreviation
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !spc_slr_cls::typ_set()

int // [enm] Return success code
spc_slr_cls::slr_cst_set(const prc_cmp &slr_cst_arg) // [W m-2] Solar constant
{
  // Purpose: Set solar constant
  std::string sbr_nm("spc_slr_cls::slr_cst_set"); // [sng] Subroutine name
  prc_cmp slr_cst_max(1400.0); // [W m-2] Maximum allowed solar constant for natural sources

  assert(slr_cst_arg >= 0.0);
  if(slr_cst_arg > slr_cst_max && abb != "lsr") err_prn(sbr_nm,"slr_cst_arg too large for natural source, proceeding anyway");
  slr_cst=slr_cst_arg; // [W m-2] Solar constant
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !spc_slr_cls::slr_cst_set()

prc_cmp // O [frc] Fraction of solar flux in band
spc_slr_cls::flx_frc_get // [fnc] Fraction of solar spectrum in a single spectral region
(const prc_cmp &wvl_min, // I [m] Minimum wavelength
 const prc_cmp &wvl_max) // I [m] Maximum wavelength
  const{
  /* Purpose: Compute fraction of solar spectrum in a single spectral region
     Routine is called with scalar input values and returns a scalar output value
     This convenience routine is an overloaded wrapper for detailed flx_frc_get() */
  int rcd_lcl(0); // O [rcd] Return success code
  prc_cmp flx_slr_frc; // [frc] Fraction of solar flux in band
  long wvl_nbr(1L); // [nbr] Number of wavelength bands
  std::string sbr_nm("spc_slr_cls::flx_frc_get"); // [sng] Subroutine name

  // Check input
  if(wvl_min > wvl_max) err_prn(sbr_nm,"wvl_min > wvl_max");
  
  rcd_lcl+=flx_frc_get // [fnc] Fraction of solar spectrum in given spectral region
    (&wvl_min, // I [m] Minimum wavelength
     &wvl_max, // I [m] Maximum wavelength
     wvl_nbr, // I [nbr] Number of wavelength bands
     &flx_slr_frc); // O [frc] Fraction of solar flux in band
  
  return flx_slr_frc; // [frc] Fraction of solar flux in band
} // !spc_slr_cls::flx_frc_get()

int // O [enm] Return success code
spc_slr_cls::flx_frc_get // [fnc] Compute fraction of solar spectrum in given spectral region
(const prc_cmp *wvl_min, // I [m] Minimum wavelength
 const prc_cmp *wvl_max, // I [m] Maximum wavelength
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 prc_cmp *flx_slr_frc) // O [frc] Fraction of solar flux in band
  const{
  /* Purpose: Compute fraction of solar spectrum in given spectral region
     Routine is called with two lists of boundaries
     The corresponding elements of each list are the minimum and maximum wavelengths
     Single element lists do not require specification of wvl_nbr which defaults to 1
     Multiple band lists must be monotonic in wavelength
     Non-monotonic lists, e.g., CAM_SW, are more difficult
     Ideally, a special go-between function, e.g., slr_spc_get_CAM_SW(), handles this */
  int rcd_lcl(0); // [enm] Return success code
  long wvl_idx; // [idx] Counting index for wvl
  std::string sbr_nm("spc_slr_cls::flx_frc_get"); // [sng] Subroutine name

  if(flx_slr_frc_fnc_ptr != (flx_slr_frc_fnc_ptr_typ)CEWI_NULL){

    // Solar spectra is available as PDF(lambda)
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      flx_slr_frc[wvl_idx]=(*flx_slr_frc_fnc_ptr)(wvl_min[wvl_idx],wvl_max[wvl_idx]); // [frc] Fraction of solar flux in band
    } // !loop over wvl_idx

  }else{ // !flx_slr_frc_fnc_ptr != CEWI_NULL

    // Solar spectra is only available in file form
    // Monotonic wavelength grids may be passed directly to ntp_slr_flx_mnt()
    // Massage non-monotonic grids, e.g., CAM_SW, RRTMG_SW, before calling ntp_slr_flx_mnt()
    if(!mnt_chk(wvl_min,wvl_nbr) || !mnt_chk(wvl_max,wvl_nbr)){
    
      rcd_lcl+=ntp_slr_flx_nmn // [fnc] Interpolate raw solar spectrum to non-monotonic grid
	(wvl_min, // I [m] Minimum wavelength in band
	 wvl_max, // I [m] Maximum wavelength in band
	 wvl_nbr, // I [nbr] Number of wavelength bands
	 flx_slr_frc); // O [frc] Fraction of solar flux in band

    }else{ // !mnt_chk
      
      rcd_lcl+=ntp_slr_flx_mnt // [fnc]  Interpolate raw solar spectrum to monotonic grid
	(wvl_min, // I [m] Minimum wavelength in band
	 wvl_max, // I [m] Maximum wavelength in band
	 wvl_nbr, // I [nbr] Number of wavelength bands
	 flx_slr_frc); // O [frc] Fraction of solar flux in band

    } // !mnt_chk

  } // endelse flx_slr_frc_fnc_ptr != CEWI_NULL
  
  return rcd_lcl; // [enm] Return success code
} // !spc_slr_cls::flx_frc_get()

// Public member functions end
// Private member functions begin

int // [enm] Return success code
spc_slr_cls::allocate(){ // [fnc] Allocate dynamic memory for object
  // Purpose: Allocate dynamic memory for object
  return rcd; // [enm] Return success code
} // !spc_slr_cls_cls::allocate()

int // [enm] Return success code
spc_slr_cls::deallocate(){ // [fnc] Free dynamic memory for object
  // Purpose: Free dynamic memory for object

  // Delete any arrays allocated by private member functions
  if(flg_raw_spc_in_dyn_mmr){
    // netCDF library in spc_slr_raw_fl_get() creates these arrays
    delete []wvl_min_in; // [m] Minimum wavelength
    delete []wvl_max_in; // [m] Maximum wavelength
    delete []flx_frc_blr_in; // [frc] Fraction of solar flux at shorter wavelengths
    delete []flx_slr_frc_in; // [frc] Fraction of solar flux in band
    // Reset flag
    flg_raw_spc_in_dyn_mmr=false; // [flg] Raw spectrum from file is in dynamic memory
  } // endif

  return rcd; // [enm] Return success code
} // !spc_slr_cls::deallocate()

int // [enm] Return success code
spc_slr_cls::recompute(){ // [fnc] Recompute properties of object
  std::string sbr_nm("spc_slr_cls::recompute()"); // [sng] Subroutine name

  // Set private members which lack set() functions
  dsc=abb2dsc(abb); // [sng] Solar flux source description
  flx_slr_frc_fnc_ptr=abb2fnc(abb); // [fnc] Pointer to function returning fractional solar fluxes
  
  /* Most spectral solar flux sources (ThD71 LaN68 Kur95 FDE24) are stored in files
     These sources are automatically loaded whenever the external file is changed in fl_spc_src_set() 
     Other sources use functions and here is where we specify the function
     This could be done by expanding the slr_spc struct to contain a function pointer
     However this method is more flexible in that it provides a hook for user-specified functions */
  if(abb == "NeL84") err_prn(sbr_nm,"Neckel & Labs option unavailable at this time."); 

  return rcd; // [enm] Return success code
} // !spc_slr_cls::recompute()

int // O [enm] Return success code
spc_slr_cls::ntp_slr_flx_mnt // [fnc]  Interpolate raw solar spectrum to monotonic grid
(const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
 const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
 const long wvl_nbr_out, // I [nbr] Number of wavelength bands
 prc_cmp *flx_frc_out, // O [frc] Fraction of solar flux in band
 const bool wrn_ntp_flg) // I [flg] Print WARNINGs from ntp_vec()
  const{
  /* Purpose:  Interpolate raw solar spectrum to monotonic grid
     Suffix _in refers to arrays stored in file (raw solar spectrum data)
     Suffix _out refers to arrays passed in to and out of this subroutine 
     Routine is C++ version of ${HOME}/f/csz_F77.F:slr_spc_get()

     Usage:
     Output wavelength grid couplets need not be contiguous in wavelength space
     wvl_min_out and wvl_max_out must be monotonic and increase in same direction

     Method:
     flx_frc_out is determined by interpolating flx_frc_blr to bin boundaries and differencing
     This ensures that total flx_frc_out=1.0 when raw spectrum data is re-gridded to any contiguous wavelength grid spanning entire solar spectrum
     Interpolating flx_slr_frc itself would not guarantee conservation of energy */

  int rcd_lcl(0); // [enm] Return success code
  long wvl_idx; // [idx] Counting index 
  std::string sbr_nm("ntp_slr_flx_mnt"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  prc_cmp *flx_frc_blr_max_out=new prc_cmp[wvl_nbr_out]; // [frc] Fraction of solar flux at shorter wavelengths
  prc_cmp *flx_frc_blr_min_out=new prc_cmp[wvl_nbr_out]; // [frc] Fraction of solar flux at shorter wavelengths

  // Print requested output wavelength grid
  if(dbg_lvl_get() >= dbg_io){ 
    for(wvl_idx=0;wvl_idx<wvl_nbr_out;wvl_idx++){
      std::cout << "wvl_min_out[" << wvl_idx << "] = " << wvl_min_out[wvl_idx] << ", wvl_max_out[" << wvl_idx << "] = " << wvl_max_out[wvl_idx] << std::endl;
    } // !wvl_idx
  } // !dbg
  
  // Check for monotonicity of output grid
  if(!mnt_chk(wvl_min_out,wvl_nbr_out)) err_prn(sbr_nm,"wvl_min_out not monotonic");
  if(!mnt_chk(wvl_max_out,wvl_nbr_out)) err_prn(sbr_nm,"wvl_max_out not monotonic");
  const bool mnt_ncr_flg=mnt_ncr(wvl_min_out,wvl_nbr_out);

  // Set extrapolation flags
  std::string xtr_sng("xtr_fll_ngh"); // [sng] Set extrapolated value to value of nearest valid neighbor
  // 20040831 fxm: c++ TODO 110 Implement way to pass xtr_vrb in constructor. Until 20040831, I set xtr_vrb by default. However, this generates ntp_vec() warnings whenever solar wavelength range is outside of file range, and these messages were overwhelming size parameter mie runs
  if(wrn_ntp_flg) xtr_sng+="+xtr_vrb"; // Add verbose diagnostics if requested
  Xtr_cls xtr_typ(xtr_sng); // [sct] Extrapolation type

  // flx_frc_blr[wvl_idx] is fractional flux bluer than wvl_max[wvl_idx], NOT wvl_ctr[wvl_idx] or wvl_min[wvl_idx]

  // Interpolate flx_frc_blr_in to wvl_min_out
  rcd_lcl+=ntp_vec // [fnc] Interpolate array to requested grid
    (wvl_nbr_in, // I [nbr] Input coordinate size
     wvl_max_in, // I [crd] Input coordinate
     flx_frc_blr_in, // I [frc] Input data (on coordinate)
     wvl_nbr_out, // I [nbr] Output coordinate size
     wvl_min_out, // I [crd] Output coordinate
     flx_frc_blr_min_out, // O [frc] Output (interpolated) data (on coordinate)
     xtr_typ, // I [enm] LHS extrapolation flags
     xtr_typ); // I [enm] RHS extrapolation flags

  // Interpolate flx_frc_blr_in to wvl_max_out
  rcd_lcl+=ntp_vec // [fnc] Interpolate array to requested grid
    (wvl_nbr_in, // I [nbr] Input coordinate size
     wvl_max_in, // I [crd] Input coordinate
     flx_frc_blr_in, // I [frc] Input data (on coordinate)
     wvl_nbr_out, // I [nbr] Output coordinate size
     wvl_max_out, // I [crd] Output coordinate
     flx_frc_blr_max_out, // O [frc] Output (interpolated) data (on coordinate)
     xtr_typ, // I [enm] LHS extrapolation flags
     xtr_typ); // I [enm] RHS extrapolation flags

  /* TOMS is an example of a monotonic wavelength grid that is not contiguous
     The TOMS grid increases monotonically
     The separation between each channel means the grid is not contiguous */
     
  // Check if wavelength grid is contiguous 
  bool ctg_flg(true); // [flg] Wavelength grid is contiguous 
  if(mnt_ncr_flg){
    for(wvl_idx=0;wvl_idx<wvl_nbr_out-1;wvl_idx++){ // Loop ends at wvl_nbr_out-2
      if(wvl_min_out[wvl_idx+1] != wvl_max_out[wvl_idx]){
	if(dbg_lvl_get() > 2) wrn_prn(sbr_nm,"Monotonically increasing grid is not contiguous wvl_min_out["+nbr2sng(wvl_idx)+"+1] != wvl_max_out["+nbr2sng(wvl_idx)+"]");
	ctg_flg=false;
	break;
      } // !wvl_min_out
    } // !wvl_idx
  }else{ // !mnt_ncr_flg
    for(wvl_idx=0;wvl_idx<wvl_nbr_out-1;wvl_idx++){ // Loop ends at wvl_nbr_out-2
      if(wvl_min_out[wvl_idx] != wvl_max_out[wvl_idx+1]){
	if(dbg_lvl_get() > 2) wrn_prn(sbr_nm,"Monotonically decreasing grid is not contiguous wvl_min_out["+nbr2sng(wvl_idx)+"] != wvl_max_out["+nbr2sng(wvl_idx)+"+1]");
	ctg_flg=false;
	break;
      } // !wvl_min_out
    } // !wvl_idx
  } // !mnt_ncr_flg

  // Check that bins abut if wavelength grid is contiguous 
  if(ctg_flg){
    if(mnt_ncr_flg){
      for(wvl_idx=0;wvl_idx<wvl_nbr_out-1;wvl_idx++) // Loop ends at wvl_nbr_out-2
	if(flx_frc_blr_min_out[wvl_idx+1] != flx_frc_blr_max_out[wvl_idx]) err_prn(sbr_nm,"Supposedly monotonically increasing grid has flx_frc_blr_min_out["+nbr2sng(wvl_idx)+"+1] != flx_frc_blr_max_out["+nbr2sng(wvl_idx)+"]");
    }else{ // !mnt_ncr_flg
      for(wvl_idx=0;wvl_idx<wvl_nbr_out-1;wvl_idx++) // Loop ends at wvl_nbr_out-2
	if(flx_frc_blr_min_out[wvl_idx] != flx_frc_blr_max_out[wvl_idx+1]) err_prn(sbr_nm,"Supposedly monotonically decreasing grid has flx_frc_blr_min_out["+nbr2sng(wvl_idx)+"] != flx_frc_blr_max_out["+nbr2sng(wvl_idx)+"+1]");
    } // !mnt_ncr_flg
  } // !ctg_flg
  
  // Convert flx_frc_blr_[min,max]_out to flx_frc_out
  for(wvl_idx=0;wvl_idx<wvl_nbr_out;wvl_idx++) flx_frc_out[wvl_idx]=flx_frc_blr_max_out[wvl_idx]-flx_frc_blr_min_out[wvl_idx];
  
  // Delete dynamically-allocated memory
  delete[] flx_frc_blr_max_out; // [frc] Fraction of solar flux at shorter wavelengths
  delete[] flx_frc_blr_min_out; // [frc] Fraction of solar flux at shorter wavelengths

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd_lcl; // [enm] Return success code
} // !spc_slr_cls::ntp_slr_flx_mnt()

int // O [enm] Return success code
spc_slr_cls::ntp_slr_flx_nmn // [fnc] Interpolate raw solar spectrum to non-monotonic grid
(const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
 const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
 const long wvl_nbr_out, // I [nbr] Number of wavelength bands
 prc_cmp *flx_frc_out) // O [frc] Fraction of solar flux in band
  const{
  /* Purpose: Wrapper routine for ntp_vec() which handles generic non-monotonic grids
     Input coordinate vector must be monotonic but output vector need not be
     Routine sorts grids and feeds them to ntp_vec() one monotonic grid at a time, 
     then reassembles them into original order
     Routine is useful for re-gridding monotonic grids to non-monotonic grids, e.g., CAM_SW
     A note on CAM_SW grid in particular is in order
     The CAM_SW grid is constructed for a k-correlation radiative transfer treatment of the near infrared (NIR)
     The mapping of these k-bands to wavelength space is made possible by the nearly monotonically increasing
     strength of the water vapor vibration-rotation bands in the NIR.
     Weighting and mapping optical properties on k-band is necessarily approximate and inexact
     Finally, CAM_SW wvl_min and wvl_max arrays to sort equivalently into same chunks of (16,1,2) monotonic bands
     More general non-montonic grids are imaginable in which wvl_min and wvl_max arrays to not sort into same chunks
     For this reason we assert() that wvl_max grid has same number of monotonic elements as wvl_min grid */ 

  // Output
  int rcd_lcl(0); // [enm] Return success code
  // Local
  long wvl_idx; // [idx] Counting index for wvl
  // Main Code
  std::string sbr_nm("ntp_slr_flx_nmn"); // [sng] Subroutine name
  unsigned short dbg_lvl(dbg_lvl_get()); // [enm] Debugging level
  long nmn_cnt(0); // [cnt] Number of monotonic chunks used in interpolation

  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  /* g++ does not allow following ntp_slr_flx_mnt() arguments to be passed by reference,
     even though they are passed by reference when ntp_slr_flx_mnt() is called directly in flx_frc_get() 
     Re-instantiate wvl_nbr_tmp as const each time because const_cast<> only works on pointers */

  // Loop through input arrays finding boundaries of monotonic chunks
  // Send each monotonic chunk to ntp_slr_flx_mnt()
  long wvl_idx_srt(0); // [idx] Starting index of current output block
  long wvl_nbr_rmn(wvl_nbr_out); // [nbr] Number of remaining wavelengths to interpolate
  long wvl_nbr_tmp; // [nbr] Number of consecutive monotonic elements of wvl_min array
  long wvl_nbr_tmp2; // [nbr] Number of consecutive monotonic elements of wvl_max array
  while(wvl_nbr_rmn > 0){

    wvl_nbr_tmp= // O [nbr] Number of consecutive monotonic elements of wvl_min array
      mnt_nbr // [fnc] Bracket monotonic portion of array
      (wvl_min_out, // I [frc] Array
       wvl_nbr_rmn, // I [nbr] Number of elements
       wvl_idx_srt); // I [idx] Index of starting element

    wvl_nbr_tmp2= // O [nbr] Number of consecutive monotonic elements of wvl_max array
      mnt_nbr // [fnc] Bracket monotonic portion of array
      (wvl_max_out, // I [frc] Array
       wvl_nbr_rmn, // I [nbr] Number of elements
       wvl_idx_srt); // I [idx] Index of starting element

    if(wvl_nbr_tmp != wvl_nbr_tmp2){
      /* CAM_LW grid triggers this problem */
      wrn_prn(sbr_nm,"Solar fluxes requested on psychotic grid type where wvl_nbr_tmp = "+nbr2sng(wvl_nbr_tmp)+" != "+nbr2sng(wvl_nbr_tmp2)+" = wvl_nbr_tmp2. Setting wvl_nbr_tmp = 1.");
      wvl_nbr_tmp=1;
    } /* !wvl_nbr_tmp */

    if(dbg_lvl_get() >= dbg_crr){
      dbg_prn(sbr_nm,"Calling ntp_slr_flx_mnt() for "+nbr2sng(wvl_nbr_tmp)+" bands beginning with band "+nbr2sng(wvl_idx_srt));
      std::cout << "idx\twvl_min\twvl_max\tflx_frc_out" << std::endl;
      std::cout << "   \t  um   \t  um   \t    frc    " << std::endl;
      for(wvl_idx=wvl_idx_srt;wvl_idx<wvl_idx_srt+wvl_nbr_tmp;wvl_idx++) std::cout << wvl_idx << "\t" << wvl_min_out[wvl_idx]*1.0e6 << "\t" << wvl_max_out[wvl_idx]*1.0e6 << "\t" << flx_frc_out[wvl_idx] << std::endl;
    } // !dbg

    rcd_lcl+=ntp_slr_flx_mnt // [fnc]  Interpolate raw solar spectrum to monotonic grid
      (const_cast<const prc_cmp *>(wvl_min_out+wvl_idx_srt), // I [m] Minimum wavelength in band
       const_cast<const prc_cmp *>(wvl_max_out+wvl_idx_srt), // I [m] Maximum wavelength in band
       wvl_nbr_tmp, // I [nbr] Number of wavelength bands
       flx_frc_out+wvl_idx_srt); // O [frc] Fraction of solar flux in band

    // Increment starting position for next iteration
    wvl_idx_srt+=wvl_nbr_tmp; // [idx] Index of starting element
    wvl_nbr_rmn-=wvl_nbr_tmp; // [nbr] Number of remaining wavelengths to interpolate
    // Increment counter
    nmn_cnt++; // [cnt] Number of monotonic chunks used in interpolation
  } // !wvl_nbr_rmn

  if(dbg_lvl == dbg_crr) std::cout << prg_nm_get() << ": INFO "+sbr_nm+"() reports calling ntp_slr_flx_mnt() " << nmn_cnt << " times" << std::endl;
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd_lcl;
} // !spc_slr_cls::ntp_slr_flx_nmn()

int // O [enm] Return success code
spc_slr_cls::ntp_slr_flx_nmn_CAM_SW // [fnc] Interpolate raw solar spectrum to CAM_SW grid
(const prc_cmp *wvl_min_out, // I [m] Minimum wavelength in band
 const prc_cmp *wvl_max_out, // I [m] Maximum wavelength in band
 const long wvl_nbr_out, // I [nbr] Number of wavelength bands
 prc_cmp *flx_frc_out) // O [frc] Fraction of solar flux in band
  const{
  /* Purpose: Wrapper routine for ntp_vec() which specifically handles CAM_SW grid
     Input coordinate vector must CAM_SW grid
     Routine sorts grids and feeds them to ntp_vec() one monotonic grid at a time, 
     then reassembles them into original order
     Routine is useful for re-gridding to non-monotonic grids, e.g., CAM_SW
     Routine is similar to csz_f77.F:rbn_vec_CCM() */ 
  // Output
  int rcd_lcl(0); // [enm] Return success code
  // Local
  long wvl_idx_srt; // [idx] Starting index of current output block
  // Main Code
  std::string sbr_nm("ntp_slr_flx_nmn_CAM_SW");
  unsigned short dbg_lvl=dbg_lvl_get();
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  if(wvl_nbr_out != 19) err_prn(sbr_nm,"Should be using generic non-monotonic solar flux interpolation ntp_slr_flx_nmn()");

  // Assume 19 bins signifies CAM_SW grid
  if(wvl_nbr_out == 19){
    // fxm: Write a more elegant solution!
    /* g++ will not allow following ntp_slr_flx_mnt() arguments to be passed by reference,
       even though they are passed by reference when ntp_slr_flx_mnt() is called directly in flx_frc_get()
       Re-instantiate wvl_nbr_tmp as const each time because const_cast<> only works on pointers */

    // Get first 16 bands
    wvl_idx_srt=0; // [idx] Starting index of current output block
    { // New scope to hide re-declaration of const wvl_nbr_tmp
      const long wvl_nbr_tmp(16); // [nbr] Number of output bins in current output block
      rcd_lcl+=ntp_slr_flx_mnt(const_cast<const prc_cmp *>(wvl_min_out+wvl_idx_srt),const_cast<const prc_cmp *>(wvl_max_out+wvl_idx_srt),wvl_nbr_tmp,flx_frc_out+wvl_idx_srt);
    } // !scope

    // Get band 17
    wvl_idx_srt=16; // [idx] Starting index of current output block
    { // New scope to hide re-declaration of const wvl_nbr_tmp
      const long wvl_nbr_tmp(1); // [nbr] Number of output bins in current output block
      rcd_lcl+=ntp_slr_flx_mnt(const_cast<const prc_cmp *>(wvl_min_out+wvl_idx_srt),const_cast<const prc_cmp *>(wvl_max_out+wvl_idx_srt),wvl_nbr_tmp,flx_frc_out+wvl_idx_srt);
    } // !scope

    // Get bands 18 and 19
    wvl_idx_srt=17; // [idx] Starting index of current output block
    { // New scope to hide re-declaration of const wvl_nbr_tmp
      const long wvl_nbr_tmp(2); // [nbr] Number of output bins in current output block
      rcd_lcl+=ntp_slr_flx_mnt(const_cast<const prc_cmp *>(wvl_min_out+wvl_idx_srt),const_cast<const prc_cmp *>(wvl_max_out+wvl_idx_srt),wvl_nbr_tmp,flx_frc_out+wvl_idx_srt);
    } // !scope
  } // endif CAM_SW

  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd_lcl;
} // !spc_slr_cls::ntp_slr_flx_nmn_CAM_SW()

int // [enm] Return success code
spc_slr_cls::spc_slr_raw_fl_get // [fnc] Load raw solar spectrum data from file
(const std::string &fl_in, // I [sng] File containing solar spectrum
 prc_cmp *&wvl_min, // O [m] Minimum wavelength
 prc_cmp *&wvl_max, // O [m] Maximum wavelength
 long &wvl_nbr, // O [nbr] Number of wavelength bands
 prc_cmp *&flx_frc_blr, // O [frc] Fraction of solar flux at shorter wavelengths
 prc_cmp *&flx_slr_frc) // O [frc] Fraction of solar flux in band
{
  // Purpose: Load raw solar spectrum data from file
  // fxm: probably should make this a private (or at least public) class function

  /* Ensure smallest wavelength bin has no flux so that rdr and blr equals [1,0] are exact!
     This is why I inserted first line of NeL84 file and last line of Kur95 files
     flx_frc_blr for bin "idx" is fraction of flux occuring at wavelengths less than wvl_max[idx] NOT wvl_ctr[idx] or wvl_min[idx]
     This convention must be supported by input text files
     This convention is similar to BPB's implementation of Lan68 and ThD71 */
  
  int rcd(0); // [enm] Return success code
  std::string sbr_nm("spc_slr_raw_fl_get"); // [sng] Subroutine name

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Check for file existance
  struct stat stat_sct;
  rcd+=stat(fl_in.c_str(),&stat_sct);
  // 20020126: fxm: icc 5.0 does not recognize std::stat
  //  rcd+=std::stat(fl_in.c_str(),&stat_sct);
  if(rcd == -1) std::perror("ERROR std::stat()");

  // Open input file
  int nc_id=nco_open(fl_in,NC_NOWRITE); // [fnc] Open netCDF file
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Opening "+fl_in);
  
  wvl_nbr=nco_inq_dimlen(nc_id,static_cast<std::string>("wvl")); // [nbr] Number of wavelengths

  // User must free this memory when no longer needed
  // For objects in spc_slr_cls, this is done in fl_slr_spc_set() which calls deallocate()
  rcd=nco_get_var(nc_id,static_cast<std::string>("wvl_min"),wvl_min); // [m] Minimum wavelength in band
  rcd=nco_get_var(nc_id,static_cast<std::string>("wvl_max"),wvl_max); // [m] Maximum wavelength in band
  rcd=nco_get_var(nc_id,static_cast<std::string>("flx_frc_blr"),flx_frc_blr); // [frc] Fraction of solar flux at shorter wavelengths
  // fxm: a better name might be flx_frc since solar flux is not required
  rcd=nco_get_var(nc_id,static_cast<std::string>("flx_frc"),flx_slr_frc); // [frc] Fraction of solar flux in band
  
  // Check for monotonicity of input grid
  if(!mnt_chk(wvl_min,wvl_nbr)) err_prn(sbr_nm,"wvl_min not monotonic");
  
  /* fxm: 20040705 valgrind reports problems here
     ==7924== Invalid read of size 1
     ==7924==    at 0x3C020CC1: strlen (mac_replace_strmem.c:189)
     ==7924==    by 0x3C0C7219: nco_get_att(int const&, int const&, std::string const&, std::string&) (char_traits.h:143)
     ==7924==    by 0x80B6F16: int nco_get_att<std::string>(int const&, std::string const&, std::string const&, std::string&) (nco_att.hh:324)
     ==7924==    by 0x80FE6E2: spc_slr_raw_fl_get(std::string const&, double*&, double*&, long&, double*&, double*&) (stl_alloc.h:652)
     ==7924==  Address 0x3C6A9801 is 0 bytes after a block of size 5 alloc'd
     ==7924==    at 0x3C0217AF: operator new[](unsigned) (vg_replace_malloc.c:113)
     ==7924==    by 0x3C0C71DA: nco_get_att(int const&, int const&, std::string const&, std::string&) (nco_att.cc:446)
     ==7924==    by 0x80B6F16: int nco_get_att<std::string>(int const&, std::string const&, std::string const&, std::string&) (nco_att.hh:324)
     ==7924==    by 0x80FE6E2: spc_slr_raw_fl_get(std::string const&, double*&, double*&, long&, double*&, double*&) (stl_alloc.h:652)
     ==7924== 
  */
  std::string wvl_units_sng;
  rcd=nco_get_att(nc_id,static_cast<std::string>("wvl"),static_cast<std::string>("units"),wvl_units_sng);
  long idx; // [idx] Counting index for wvl
  if(wvl_units_sng.find("micron") != std::string::npos){
    for(idx=0;idx<wvl_nbr;idx++){
      wvl_min[idx]/=1.0e6; // [um] -> [m]
      wvl_max[idx]/=1.0e6; // [um] -> [m]
    } // !idx
  } // !wvl_units_sng
  
  // Close input file
  rcd=nco_close(nc_id); // [fnc] Close netCDF file
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Closing "+fl_in);
  
  if(dbg_lvl_get() > dbg_old){
    std::cout << "idx\twvl_min\twvl_max\tflx_frc_blr" << std::endl;
    std::cout << "   \t  um   \t  um   \t    frc    " << std::endl;
    for(idx=0;idx<wvl_nbr;idx++) std::cout << idx << "\t" << wvl_min[idx]*1.0e6 << "\t" << wvl_max[idx]*1.0e6 << "\t" << flx_frc_blr[idx] << std::endl;
  } // !dbg

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !spc_slr_cls::spc_slr_raw_fl_get()

// Private member functions end
// Global functions with C++ linkages begin

prc_cmp // [frc] Fraction of solar flux in band
flx_slr_frc_lsr // [fnc] Solar flux of laser
(const prc_cmp &wvl_min, // [m] Minimum wavelength in band
 const prc_cmp &wvl_max) // [m] Maximum wavelength in band
{
  /* Purpose: Return fractional solar flux in band
     Currently function is dummy function to simulate laser source
     Function returns unity for any input */

  // Output
  //  int rcd(0); // [enm] Return success code
  std::string sbr_nm("flx_slr_frc_lsr"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  if(wvl_min == wvl_max){;} // CEWU Compiler Error Warning Usage

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return 1.0;
} // !flx_slr_frc_lsr()

int // O [enm] Return success code
flx_slr_frc_lsr // [fnc] Solar flux of laser
(const prc_cmp *wvl_min, // I [m] Minimum wavelength in band
 const prc_cmp *wvl_max, // I [m] Maximum wavelength in band
 const long &wvl_nbr, // I [nbr] Number of wavelength bands
 prc_cmp *flx_slr_frc) // O [frc] Fraction of solar flux in band
{
  /* Purpose: Return fractional solar flux in band
     Currently function is dummy function to simulate laser source
     Function homogeneous square wave normalized to unity for any input */

  // Output
  int rcd(0); // [enm] Return success code
  long wvl_idx; // [idx] Counting index for wvl
  std::string sbr_nm("flx_slr_frc_lsr"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  prc_cmp flx_tmp(1.0/wvl_nbr); // [frc] Fraction of solar flux in band
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    flx_slr_frc[wvl_idx]=flx_tmp; // [frc] Fraction of solar flux in band
  } // !wvl_idx

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");

  if(wvl_min[0] == wvl_max[0]){;} // CEWU Compiler Error Warning Usage

  return rcd;
} // !flx_slr_frc_lsr()

// fxm: This function does nothing useful
prc_cmp // [frc] Fraction of solar flux in band
flx_slr_frc_ThD71 // [fnc] Fraction of solar flux in given spectral region of ThD71
(const prc_cmp &wvl_min_arg, // [m] Minimum wavelength
 const prc_cmp &wvl_max_arg) // [m] Maximum wavelength
{
  // Purpose: Compute fraction of solar flux in given spectral region of ThD71
  std::string sbr_nm("flx_slr_frc_ThD71"); // [sng] Subroutine name
  long idx; // [idx] Counting idx
  const long wvl_nbr_ThD71(68); // [nbr] Number of wavelength bins

  if(wvl_min_arg > wvl_max_arg) err_prn(sbr_nm,"wvl_min_arg > wvl_max_arg");
  const long wvl_nbr(wvl_nbr_ThD71); // [nbr] Number of wavelength bands

  // Create standard arrays
  prc_cmp *wvl=new prc_cmp[wvl_nbr]; // [m] Nominal wavelength
  prc_cmp *wvl_dlt=new prc_cmp[wvl_nbr]; // [m] Bandwidth
  prc_cmp *wvl_min=new prc_cmp[wvl_nbr]; // [m] Minimum wavelength in band
  prc_cmp *wvl_max=new prc_cmp[wvl_nbr]; // [m] Maximum wavelength in band
  prc_cmp *wvl_ctr=new prc_cmp[wvl_nbr]; // [m] Wavelength at band center
  prc_cmp *flx_slr_frc=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux in band
  prc_cmp foo(0.0); // [frc] fxm

  /* Solar flux spectrum of Thaekeakara and Drummond (1971)
     CSZ added (wvl_mcr,flx_frc_blr) = (10.0,1.0)
     As input by BPB, last element was (5.0,0.9951) */
  static prc_cmp wvl_max_tmp[wvl_nbr_ThD71]={
    .15 ,
    .20 , .22 , .23 , .24 , .25 , .26  ,  .27 ,
    .28 ,
    .29 , .30 , .31 , .32 , .33 , .34  ,  .35 ,
    .36 , .37 , .38 , .39 , .40 , .41  ,  .42 ,
    .43 , .44 , .45 , .46 , .47 , .48  ,  .49 ,
    .50 , .51 , .52 , .53 , .54 , .55  ,  .56 ,
    .57 ,
    .58 , .59 , .60 , .62 , .64 , .66  ,  .68 ,
    .70 , .72 , .75 , .80 , .90 ,1.00  , 1.20 ,
    1.40 ,1.60 ,1.80 ,2.00 ,2.20 ,2.40  , 2.60 ,
    2.80 ,3.00 ,3.20 ,3.40 ,3.60 ,3.80  , 4.00 ,
    5.00,
    10.00
  };
  
  prc_cmp flx_frc_blr_tmp[wvl_nbr_ThD71]={
    .0000,
    .0001,.0005,.0010,.0014,.0019,.0027 ,.0041 ,
    .0056,
    .0081,.0121,.0165,.0222,.0293,.0372 ,.0452 ,
    .0532,.0615,.0700,.0782,.0873,.0992 ,.1122 ,
    .1247,.1373,.1514,.1665,.1817,.1968 ,.2115 ,
    .2260,.2401,.2538,.2674,.2808,.2938 ,.3065 ,
    .3191,
    .3318,.3444,.3568,.3810,.4042,.4266 ,.4481 ,
    .4688,.4886,.5169,.5602,.6336,.6946 ,.7839 ,
    .8434,.8861,.9159,.9349,.9483,.9589 ,.9667 ,
    .9731,.9783,.9822,.9850,.9872,.9891 ,.9906 ,
    .9951,
    1.
  };
  
  // Convert tabular data to SI units
  for(idx=0;idx<wvl_nbr;idx++) wvl_max_tmp[idx]*=1.0e-6; // [um]->[m]
  if(dbg_lvl_get() > 1){
    std::cout << "idx\twvl_max\tflx_frc_blr" << std::endl;
    std::cout << "#\tum\t" << std::endl;
    for(idx=0;idx<wvl_nbr;idx++) std::cout << idx << "\t" << wvl_max_tmp[idx]*1.0e6 << "\t" << flx_frc_blr_tmp[idx] << std::endl;
  } // !dbg
  
  // Create standard fields
  for(idx=0;idx<wvl_nbr;idx++){
    wvl_max[idx]=wvl_max_tmp[idx]; // [m] Maximum wavelength in band
  } // !idx
  
  wvl_min[0]=wvl_max[0]-(wvl_max[1]-wvl_max[0]);
  for(idx=1;idx<wvl_nbr;idx++){
    wvl_min[idx]=wvl_max[idx-1]; // [m] Minimum wavelength in band
  } // !idx
  
  for(idx=0;idx<wvl_nbr;idx++){
    wvl[idx]=wvl_ctr[idx]=0.5*(wvl_max[idx]+wvl_min[idx]); // [m] Wavelength at band center
    wvl_dlt[idx]=wvl_max[idx]-wvl_min[idx]; // [m] Bandwidth
  } // !idx
  
  flx_slr_frc[0]=flx_frc_blr_tmp[0]; // [W m-2] Solar flux in band
  for(idx=1;idx<wvl_nbr;idx++){
    flx_slr_frc[idx]=flx_frc_blr_tmp[idx]-flx_frc_blr_tmp[idx-1]; // [frc] Fraction of solar flux in band
  } // !idx
  
  // Free dynamic memory
  delete[] wvl; // [m] Nominal wavelength
  delete[] wvl_dlt; // [m] Bandwidth
  delete[] wvl_min; // [m] Minimum wavelength in band
  delete[] wvl_max; // [m] Maximum wavelength in band
  delete[] wvl_ctr; // [m] Wavelength at band center
  delete[] flx_slr_frc; // [frc] Fraction of solar flux in band

  // fxm: return whole array
  return foo; // [frc]
} // !flx_slr_frc_ThD71()

// Global functions with C++ linkages end
