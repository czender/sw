// $Id$ 

// Implementation (declaration) of wavelength grid class

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <wvl_grd.hh> // Wavelength grids

// wvl_grd_cls class

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const wvl_grd_cls &wvl_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for wvl_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     std::cout << wvl_obj; */
  const long wvl_nbr(wvl_obj.wvl_nbr_get()); // [nbr] Number of wavelengths
  const long wvl_idx_dbg(wvl_obj.wvl_idx_dbg); // [idx] Debugging wavelength bin

  srm_out << "Public contents of wvl_obj: " << std::endl;
  srm_out << "wvl_obj.nst_nbr_get() = " << wvl_obj.nst_nbr_get() << std::endl;
  srm_out << "wvl_obj.wvl_nbr_get() = " << wvl_obj.wvl_nbr_get() << std::endl;
  srm_out << "Private contents of wvl_obj: " << std::endl;
  srm_out << "wvl_obj.grd_sng = " << wvl_obj.grd_sng << std::endl;
  srm_out << "wvl_idx_dbg = " << wvl_obj.wvl_idx_dbg << std::endl;
  if(wvl_obj.flg_dyn_mmr){
    if(wvl_obj.flg_wvl_grd){
      srm_out << "wvl_min = " << wvl_obj.wvl_min[0]*1.0e6 << " um" << std::endl;
      srm_out << "wvl_max = " << wvl_obj.wvl_max[wvl_nbr-1]*1.0e6 << " um" << std::endl;
      srm_out << "wvl_dbg = " << wvl_obj.wvl_ctr[wvl_idx_dbg]*1.0e6 << " um" << std::endl;
    } // endif flg_wvl_grd
    if(wvl_obj.flg_wvn_grd){
      srm_out << "wvn_min = " << wvl_obj.wvn_min[0] << " cm-1" << std::endl;
      srm_out << "wvn_max = " << wvl_obj.wvn_max[wvl_nbr-1] << " cm-1" << std::endl;
      srm_out << "wvn_dbg = " << wvl_obj.wvn_ctr[wvl_idx_dbg] << " cm-1" << std::endl;
    } // endif flg_wvn_grd
    if(wvl_obj.flg_frq_grd){
      srm_out << "frq_min = " << wvl_obj.frq_min[0]/1.0e12 << " THz" << std::endl;
      srm_out << "frq_max = " << wvl_obj.frq_max[wvl_nbr-1]/1.0e12 << " THz" << std::endl;
      srm_out << "frq_dbg = " << wvl_obj.frq_ctr[wvl_idx_dbg]/1.0e12 << " THz"  << std::endl;
    } // endif flg_frq_grd
  }else{
    srm_out << "No wavelength grids in dynamic memory" << std::endl;
  } // endif
  srm_out << "sizeof(wvl_obj) = " << sizeof(wvl_obj) << " B" << std::endl;
    
  return srm_out; // [srm] Reference to output stream for cascading
} // end operator<<()

// Friendly functions end
// Static members begin

int wvl_grd_cls::nst_nbr=0; // [nbr] Number of instantiated class members
const long wvl_grd_cls::wvl_nbr_max(100000); // [nbr] Maximum number of wavelengths
const prc_cmp wvl_grd_cls::wvl_dbg_dfl(0.5e-6); // [m] Debugging wavelength, default

// Static members end
// Static member functions begin
int wvl_grd_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members
sng2sng_map wvl_grd_cls::opt2abb_map=wvl_grd_cls::opt2abb_map_mk();
sng2sng_map wvl_grd_cls::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;
  map_tmp.insert(sng2sng_map::value_type("log","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("lgr","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("logarithmic","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("lin","regular"));
  map_tmp.insert(sng2sng_map::value_type("lnr","regular"));
  map_tmp.insert(sng2sng_map::value_type("regular","regular"));
  map_tmp.insert(sng2sng_map::value_type("linear","regular"));
  map_tmp.insert(sng2sng_map::value_type("rgl","regular"));
  map_tmp.insert(sng2sng_map::value_type("RGL","regular"));
  map_tmp.insert(sng2sng_map::value_type("REG","regular"));
  map_tmp.insert(sng2sng_map::value_type("reg","regular"));
  map_tmp.insert(sng2sng_map::value_type("wvn_rgl","wvn_rgl"));
  map_tmp.insert(sng2sng_map::value_type("sns","sensor"));
  map_tmp.insert(sng2sng_map::value_type("sensor","sensor"));
  map_tmp.insert(sng2sng_map::value_type("sensors","sensor"));
  map_tmp.insert(sng2sng_map::value_type("rnt","AERONET"));
  map_tmp.insert(sng2sng_map::value_type("aeronet","AERONET"));
  map_tmp.insert(sng2sng_map::value_type("AERONET","AERONET"));
  map_tmp.insert(sng2sng_map::value_type("MODIS","MODIS"));
  map_tmp.insert(sng2sng_map::value_type("MDS","MODIS"));
  map_tmp.insert(sng2sng_map::value_type("mds","MODIS"));
  map_tmp.insert(sng2sng_map::value_type("modis","MODIS"));
  map_tmp.insert(sng2sng_map::value_type("toms","TOMS"));
  map_tmp.insert(sng2sng_map::value_type("tms","TOMS"));
  map_tmp.insert(sng2sng_map::value_type("TOMS","TOMS"));
  map_tmp.insert(sng2sng_map::value_type("CAM_SW","CAM_SW"));
  map_tmp.insert(sng2sng_map::value_type("CCM_SW","CAM_SW"));
  map_tmp.insert(sng2sng_map::value_type("CCSM_SW","CAM_SW"));
  map_tmp.insert(sng2sng_map::value_type("SW","CAM_SW"));
  map_tmp.insert(sng2sng_map::value_type("sw","CAM_SW"));
  map_tmp.insert(sng2sng_map::value_type("fastj_aer","fastj_aer"));
  map_tmp.insert(sng2sng_map::value_type("FASTJ_aer","fastj_aer"));
  map_tmp.insert(sng2sng_map::value_type("FastJ_aer","fastj_aer"));
  map_tmp.insert(sng2sng_map::value_type("fastj_gas","fastj_gas"));
  map_tmp.insert(sng2sng_map::value_type("FASTJ_gas","fastj_gas"));
  map_tmp.insert(sng2sng_map::value_type("FastJ_gas","fastj_gas"));
  map_tmp.insert(sng2sng_map::value_type("SWNB","SWNB"));
  map_tmp.insert(sng2sng_map::value_type("swnb","SWNB"));
  map_tmp.insert(sng2sng_map::value_type("WVL_DBG","WVL_DBG"));
  map_tmp.insert(sng2sng_map::value_type("dbg","WVL_DBG"));
  map_tmp.insert(sng2sng_map::value_type("DBG","WVL_DBG"));
  map_tmp.insert(sng2sng_map::value_type("CAM_LW","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("CCM_LW","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("CCSM_LW","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("RRTM_LW","RRTM_LW"));
  map_tmp.insert(sng2sng_map::value_type("RRTM_SW","RRTM_SW"));
  map_tmp.insert(sng2sng_map::value_type("LW","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("lw","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("CCM LW","CAM_LW"));
  map_tmp.insert(sng2sng_map::value_type("GSFC_LW","GSFC_LW"));
  map_tmp.insert(sng2sng_map::value_type("GSFC","GSFC_LW"));
  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end wvl_grd_cls::opt2abb_map_mk()
std::string wvl_grd_cls::opt2abb(const std::string opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("wvl_grd_cls::opt2abb",opt_sng+" is unknown");
  return itr->second;
} // end wvl_grd_cls::opt2abb()

// Static member functions end
// Public member functions begin

wvl_grd_cls::wvl_grd_cls // [fnc] Default constructor
(const std::string grd_sng_arg, // I [sng] Type of wavelength grid
 prc_cmp &wvl_mnm_arg, // I/O [m] Minimum wavelength
 prc_cmp &wvl_mxm_arg, // I/O [m] Maximum wavelength
 long &wvl_nbr_arg, // I/O [nbr] Number of wavelengths
 const bool flg_wvl_grd_arg, // I [flg] Compute wavelength grid
 const bool flg_wvn_grd_arg, // I [flg] Compute wavenumber grid
 const bool flg_frq_grd_arg) // I [flg] Compute frequency grid
{
  /* Purpose: Default wavelength grid constructor
     Note that specifying a pre-defined wavelength grid sets private values 
     of wvl_nbr, wvl_mnm, wvl_mxm. 
     Routine automatically alters calling values to values imposed by grid choice
     Doing this from the calling program using the *_get() functions is error-prone */
  std::string sbr_nm("wvl_grd_cls::wvl_grd_cls"); // [sng] Subroutine name

  rcd=0; // [enm] Return success code
  flg_dyn_mmr=false; // [flg] Dynamic memory arrays have been allocated
  wvl_mnm=wvl_mnm_arg; // [m] Minimum wavelength
  wvl_mxm=wvl_mxm_arg; // [m] Maximum wavelength
  wvl_nbr=wvl_nbr_arg; // [nbr] Number of wavelengths
  grd_sng=grd_sng_arg; // [sng] Type of wavelength grid
  // For now, always compute wavelength grid
  assert(flg_wvl_grd_arg == true);
  flg_wvl_grd=flg_wvl_grd_arg; // [flg] Compute wavelength grid
  flg_wvn_grd=flg_wvn_grd_arg; // [flg] Compute wavenumber grid
  flg_frq_grd=flg_frq_grd_arg; // [flg] Compute frequency grid
  
  /* Allocate space here for user-specified grids only
     recompute() calls allocate() for all pre-defined grids
     Rationale is explained in recompute() */
  if(grd_sng == "regular" || grd_sng == "logarithmic" || grd_sng == "wvn_rgl") rcd+=allocate(); // [fnc] Allocate dynamic memory for object
  rcd+=recompute(); // [fnc] Recompute properties of object
  
  /* Ensure calling program values reflect pre-defined grid, if any */
  wvl_mnm_arg=wvl_mnm; // [m] Minimum wavelength
  wvl_mxm_arg=wvl_mxm; // [m] Maximum wavelength
  wvl_nbr_arg=wvl_nbr; // [nbr] Number of wavelengths

  if(rcd != 0) err_prn(sbr_nm,"Failed constructor"); // [fnc] Print uniform error message and exit

  nst_nbr++; // [nbr] Number of instantiated class members
} // end wvl_grd_cls constructor

wvl_grd_cls::~wvl_grd_cls(){ // Destructor
  rcd+=deallocate(); // [fnc] Free dynamic memory of object
  nst_nbr--; // [nbr] Number of instantiated class members
} // end wvl_grd_cls destructor

int // [enm] Return success code
wvl_grd_cls::wvl_grd_min_max_set // [fnc] Reset wavelength grid boundaries
(const prc_cmp &wvl_mnm_arg, // [m] Minimum wavelength
 const prc_cmp &wvl_mxm_arg) // [m] Maximum wavelength
{
  rcd=0; // [enm] Return success code
  wvl_mnm=wvl_mnm_arg; // [m] Minimum wavelength
  wvl_mxm=wvl_mxm_arg; // [m] Maximum wavelength
  rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end wvl_grd_cls::wvl_grd_min_max_set()

std::string wvl_grd_cls::typ()const{return grd_sng;} // [fnc] Get grid string
void wvl_grd_cls::typ(const std::string sng){ // [fnc] Set grid type
  grd_sng=sng;
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end wvl_grd_cls::typ()

long wvl_grd_cls::wvl_nbr_get()const{return wvl_nbr;} // [nbr] Number of wavelength bands
int // [enm] Return success code
wvl_grd_cls::wvl_nbr_set(const long &wvl_nbr_arg) // [nbr] Number of wavelength bands
{
  // Purpose: Set number of wavelength bands
  // Use routine to change resolution for given spectral region
  std::string sbr_nm("wvl_grd_cls::wvl_nbr_set"); // [sng] Subroutine name

  // Sanity check
  assert(wvl_nbr_arg > 0 && wvl_nbr_arg <= wvl_nbr_max); // [nbr] Allowed number of wavelength bands
  wvl_nbr=wvl_nbr_arg; // [nbr] Number of wavelength bands

  if(grd_sng != "regular" && grd_sng != "logarithmic") err_prn(sbr_nm,"Attempt to set wavelength number for pre-defined grid. Re-setting wavelength number would violate grid definition.");
  rcd+=recompute(); // [fnc] Recompute properties of object

  return rcd; // [enm] Return success code
} // end wvl_grd_cls::wvl_nbr_set()

prc_cmp *wvl_grd_cls::wvl_ctr_get()const{assert(flg_wvl_grd == true); return wvl_ctr;} // [m] Wavelength at band center
prc_cmp *wvl_grd_cls::wvl_dlt_get()const{assert(flg_wvl_grd == true); return wvl_dlt;} // [m] Bandwidth
prc_cmp *wvl_grd_cls::wvl_grd_get()const{assert(flg_wvl_grd == true); return wvl_grd;} // [m] Wavelength grid
prc_cmp *wvl_grd_cls::wvl_max_get()const{assert(flg_wvl_grd == true); return wvl_max;} // [m] Maximum wavelength in band
prc_cmp *wvl_grd_cls::wvl_min_get()const{assert(flg_wvl_grd == true); return wvl_min;} // [m] Minimum wavelength in band

prc_cmp *wvl_grd_cls::wvn_ctr_get()const{assert(flg_wvn_grd == true); return wvn_ctr;} // [cm-1] Wavenumber at band center
prc_cmp *wvl_grd_cls::wvn_dlt_get()const{assert(flg_wvn_grd == true); return wvn_dlt;} // [cm-1] Bandwidth
prc_cmp *wvl_grd_cls::wvn_grd_get()const{assert(flg_wvn_grd == true); return wvn_grd;} // [cm-1] Wavenumber grid
prc_cmp *wvl_grd_cls::wvn_max_get()const{assert(flg_wvn_grd == true); return wvn_max;} // [cm-1] Maximum wavenumber in band
prc_cmp *wvl_grd_cls::wvn_min_get()const{assert(flg_wvn_grd == true); return wvn_min;} // [cm-1] Minimum wavenumber in band

prc_cmp *wvl_grd_cls::frq_ctr_get()const{assert(flg_frq_grd == true); return frq_ctr;} // [Hz] Frequency at band center
prc_cmp *wvl_grd_cls::frq_dlt_get()const{assert(flg_frq_grd == true); return frq_dlt;} // [Hz] Bandwidth
prc_cmp *wvl_grd_cls::frq_grd_get()const{assert(flg_frq_grd == true); return frq_grd;} // [Hz] Frequency grid
prc_cmp *wvl_grd_cls::frq_max_get()const{assert(flg_frq_grd == true); return frq_max;} // [Hz] Maximum frequency in band
prc_cmp *wvl_grd_cls::frq_min_get()const{assert(flg_frq_grd == true); return frq_min;} // [Hz] Minimum frequency in band

// Public member functions end
// Private member functions begin

int // [enm] Return success code
wvl_grd_cls::allocate(){ // [fnc] Allocate dynamic memory for object
  if(flg_wvl_grd){
    wvl_ctr=new prc_cmp[wvl_nbr]; // [m] Wavelength at band center
    wvl_dlt=new prc_cmp[wvl_nbr]; // [m] Bandwidth
    wvl_grd=new prc_cmp[wvl_nbr+1]; // [m] Wavelength grid
    wvl_max=new prc_cmp[wvl_nbr]; // [m] Maximum wavelength in band
    wvl_min=new prc_cmp[wvl_nbr]; // [m] Minimum wavelength in band
  } // endif flg_wvl_grd
  if(flg_wvn_grd){
    wvn_ctr=new prc_cmp[wvl_nbr]; // [cm-1] Wavenumber at band center
    wvn_dlt=new prc_cmp[wvl_nbr]; // [cm-1] Bandwidth
    wvn_grd=new prc_cmp[wvl_nbr+1]; // [cm-1] Wavenumber grid
    wvn_max=new prc_cmp[wvl_nbr]; // [cm-1] Maximum wavenumber in band
    wvn_min=new prc_cmp[wvl_nbr]; // [cm-1] Minimum wavenumber in band
  } // endif flg_wvn_grd
  if(flg_frq_grd){
    frq_ctr=new prc_cmp[wvl_nbr]; // [Hz] Frequency at band center
    frq_dlt=new prc_cmp[wvl_nbr]; // [Hz] Bandwidth
    frq_grd=new prc_cmp[wvl_nbr+1]; // [Hz] Frequency grid
    frq_max=new prc_cmp[wvl_nbr]; // [Hz] Maximum frequency in band
    frq_min=new prc_cmp[wvl_nbr]; // [Hz] Minimum frequency in band
  } // endif flg_frq_grd
  flg_dyn_mmr=true; // [flg] Dynamic memory arrays have been allocated
  return rcd; // [enm] Return success code
} // end wvl_grd_cls::allocate()

int // [enm] Return success code
wvl_grd_cls::deallocate(){ // [fnc] Free dynamic memory of object
  if(flg_dyn_mmr){ // [flg] Dynamic memory arrays have been allocated
    if(flg_wvl_grd){
      delete []wvl_ctr; // [m] Wavelength at band center
      delete []wvl_dlt; // [m] Bandwidth
      delete []wvl_grd; // [m] Wavelength grid
      delete []wvl_max; // [m] Maximum wavelength in band
      delete []wvl_min; // [m] Minimum wavelength in band
    } // endif flg_wvl_grd
    if(flg_wvn_grd){
      delete []wvn_ctr; // [cm-1] Wavenumber at band center
      delete []wvn_dlt; // [cm-1] Bandwidth
      delete []wvn_grd; // [cm-1] Wavenumber grid
      delete []wvn_max; // [cm-1] Maximum wavenumber in band
      delete []wvn_min; // [cm-1] Minimum wavenumber in band
    } // endif flg_wvn_grd
    if(flg_frq_grd){
      delete []frq_ctr; // [Hz] Frequency at band center
      delete []frq_dlt; // [Hz] Bandwidth
      delete []frq_grd; // [Hz] Frequency grid
      delete []frq_max; // [Hz] Maximum frequency in band
      delete []frq_min; // [Hz] Minimum frequency in band
    } // endif flg_frq_grd
    flg_dyn_mmr=false; // [flg] Dynamic memory arrays have been allocated
  } // [flg] Dynamic memory arrays have been allocated
  return rcd; // [enm] Return success code
} // end wvl_grd_cls::deallocate()

int // [enm] Return success code
wvl_grd_cls::recompute(){ // [fnc] Recompute properties of object
  /* Purpose: Generate a wavelength grid
     Grid-specific portion of routine creates wvl_min and wvl_max arrays
     Common code then fills in rest of structure based on wvl_min and wvl_max
     Routine also calls allocate() for all wavelength grids except for completely 
     user-specified grids, i.e., regular and logarithmic grids
     Constructor or wvl_nbr_set() allocates memory for user-specified grids
     This ensures user-specified grids are not automatically allocated/de-allocated
     with each re-compute */
  std::string sbr_nm("wvl_grd_cls::recompute"); // [sng] Subroutine name
  long wvl_idx; // [idx] Counting index for wvl
  if(grd_sng == "CAM_SW" || grd_sng == "CAM_LW" || grd_sng == "GSFC_LW" || grd_sng == "RRTM_SW"){
    std::cerr << prg_nm_get() << ": HINT: Monotonicity WARNINGs may be safely ignored for this gridtype " << grd_sng << std::endl;
  } // endif
  // Create wvl_min and wvl_max
  if(grd_sng == "logarithmic"){ // Logarithmically spaced wavelength arrays
    prc_cmp max_min_ratio=wvl_mxm/wvl_mnm;
    prc_cmp series_ratio=std::pow(max_min_ratio,PRC_CMP(1.0)/wvl_nbr);
    if(wvl_mnm==0.0) err_prn(sbr_nm,"wvl_mnm = 0.0 for wavelength grid type = "+grd_sng);
    wvl_min[0]=wvl_mnm;
    for(wvl_idx=1;wvl_idx<wvl_nbr;wvl_idx++){ // NB: Loop starts from 1
      wvl_min[wvl_idx]=wvl_min[wvl_idx-1]*series_ratio;
      wvl_max[wvl_idx-1]=wvl_min[wvl_idx];
    } // end loop over wvl_idx
    wvl_max[wvl_nbr-1]=wvl_mxm; // [m] Ensure final wvl_max is not affected by roundoff
  }else if(grd_sng == "regular"){ // Regularly spaced wavelength arrays
    // Do not deallocate()/allocate() memory for regular grids here, see function header for why
    prc_cmp wvl_ncr=(wvl_mxm-wvl_mnm)/wvl_nbr; // [m]
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_mnm+wvl_idx*wvl_ncr; // [m]
      wvl_max[wvl_idx]=wvl_mnm+(wvl_idx+1)*wvl_ncr; // [m]
    } // end loop over wvl_idx
    wvl_max[wvl_nbr-1]=wvl_mxm; // [m] Ensure final wvl_max is not affected by roundoff
  }else if(grd_sng == "wvn_rgl"){ // Regularly spaced wavenumber arrays
    // Do not deallocate()/allocate() memory for regular grids here, see function header for why
    wvn_mnm=0.01/wvl_mxm; // [cm-1] Minimum wavenumber
    wvn_mxm=0.01/wvl_mnm; // [cm-1] Maximum wavenumber
    wvn_nbr=wvl_nbr;
    prc_cmp wvn_ncr=(wvn_mxm-wvn_mnm)/wvn_nbr; // [cm-1]
    long wvn_idx; // [idx] Counting index for wvn
    for(wvn_idx=0;wvn_idx<wvn_nbr;wvn_idx++){
      wvn_min[wvn_idx]=wvn_mnm+wvn_idx*wvn_ncr; // [cm-1]
      wvn_max[wvn_idx]=wvn_mnm+(wvn_idx+1)*wvn_ncr; // [cm-1]
    } // end loop over wvn_idx
    wvn_max[wvn_nbr-1]=wvn_mxm; // [m] Ensure final wvn_max is not affected by roundoff
    // Reverse monotonicity so output increases with wavelength not wavenumber
    //rvr_vec(wvn_min,wvn_nbr);
    //rvr_vec(wvn_max,wvn_nbr);
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=0.01/wvn_max[wvl_idx]; // [cm-1] -> [m]
      wvl_max[wvl_idx]=0.01/wvn_min[wvl_idx]; // [cm-1] -> [m]
    } // end loop over wvl_idx
  }else if(grd_sng == "WVL_DBG"){ // Wavelength grid for debugging
    wvl_mnm=0.18e-6; // [m]
    wvl_mxm=0.22e-6; // [m]
    const long wvl_nbr_WVL_DBG(10); // [nbr]
    if(wvl_nbr != wvl_nbr_WVL_DBG){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_WVL_DBG; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    prc_cmp wvl_ncr=(wvl_mxm-wvl_mnm)/wvl_nbr; // [m]
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_mnm+wvl_idx*wvl_ncr; // [m]
      wvl_max[wvl_idx]=wvl_mnm+(wvl_idx+1)*wvl_ncr; // [m]
    } // end loop over wvl_idx
    wvl_max[wvl_nbr-1]=wvl_mxm; // [m] Ensure final wvl_max is not affected by roundoff
  }else if(grd_sng == "CAM_SW"){ // 
    /* JGW supplied wavelength grid for CCM NIR used in WKB98
       "The band borders were calculated using those "p-weights" that
       Bruce Briegleb lists in his paper on the delta-Eddington"
       19980715: Added Nimbus 7 NIR band to agree with CCM3.6 */
    /* 01      0.2000  0.2450
       02      0.2450  0.2650
       03      0.2650  0.2750
       04      0.2750  0.2850
       05      0.2850  0.2950
       06      0.2950  0.3050
       07      0.3050  0.3500
       08      0.3500  0.6400
       09      0.6400  0.7000
       10      0.7000  1.0890
       11      1.0890  1.4080
       12      1.4080  1.7490
       13      1.7490  2.1520
       14      2.1520  2.7140
       15      2.7140  3.4570
       16      3.4570  5.0000
       17      2.6300  2.8600
       18      4.1600  4.4750
       19      4.4750  4.5500 */
    const prc_cmp wvl_min_CAM_SW_mcr[]={0.2000, 0.2450, 0.2650, 0.2750, 0.2850,
					0.2950, 0.3050, 0.3500, 0.6400, 0.7000,  
					1.0890, 1.4080, 1.7490, 2.1520, 2.7140, 
					3.4570, 2.6300, 4.1600, 4.4750};
    const prc_cmp wvl_max_CAM_SW_mcr[]={0.2450, 0.2650, 0.2750, 0.2850, 0.2950,
					0.3050, 0.3500, 0.6400, 0.7000, 1.0890, 
					1.4080, 1.7490, 2.1520, 2.7140, 3.4570, 
					5.0000, 2.8600, 4.4750, 4.5500};
    const long wvl_nbr_CAM_SW(sizeof(wvl_min_CAM_SW_mcr)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_CAM_SW){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_CAM_SW; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_CAM_SW_mcr[wvl_idx]*1.0e-6; // [um] -> [m]
      wvl_max[wvl_idx]=wvl_max_CAM_SW_mcr[wvl_idx]*1.0e-6; // [um] -> [m]
    } // end for
    // end CAM_SW-specific initialization
  }else if(grd_sng == "fastj_aer"){ // 
    /* Fast J aerosol wavelength grid
       MJP's group list these bins as having effective (solar-irradiance weighted) centers at:
       355.0, 500.0, 800.0 nm, respectively */
    const prc_cmp wvl_min_fastj_aer_nm[]={300.0,400.0,600.0};
    const prc_cmp wvl_max_fastj_aer_nm[]={400.0,600.0,999.0};
    const long wvl_nbr_fastj_aer_nm(sizeof(wvl_min_fastj_aer_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_fastj_aer_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_fastj_aer_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_fastj_aer_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_fastj_aer_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end fastj_aer-specific initialization
  }else if(grd_sng == "fastj_gas"){ // 
    /* Fast J gas wavelength grid
       MJP's group list these bins as having effective (solar-irradiance weighted) centers at:
       294.0, 303.0, 310.0, 316.0, 333.0, 380.0, 574.0 nm, respectively */
    const prc_cmp wvl_min_fastj_gas_nm[]={289.00,298.25,307.45,312.45,320.30,345.00,412.45};
    const prc_cmp wvl_max_fastj_gas_nm[]={298.25,307.45,312.45,320.30,345.00,412.45,850.00};
    const long wvl_nbr_fastj_gas_nm(sizeof(wvl_min_fastj_gas_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_fastj_gas_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_fastj_gas_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_fastj_gas_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_fastj_gas_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end fastj_gas-specific initialization
  }else if(grd_sng == "AERONET"){ // 
    /* AERONET sensor retrievals described/used in DHE02, STD03
       AERONET inversion algorithm provides optical properties at four wavelengths
       (Somewhere I recall reading that the AERONET sensor bandpass is 10 nm) */
    const prc_cmp wvl_min_AERONET_nm[]={435.0,665.0,865.0,1015.0};
    const prc_cmp wvl_max_AERONET_nm[]={445.0,675.0,875.0,1025.0};
    const long wvl_nbr_AERONET_nm(sizeof(wvl_min_AERONET_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_AERONET_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_AERONET_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_AERONET_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_AERONET_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end AERONET-specific initialization
  }else if(grd_sng == "sensor"){ // 
    /* "sensor" wavelengths are all AERONET, MODIS, TOMS, and WMO channels */
    const prc_cmp wvl_min_sensor_nm[]={330.7,339.16,359.90,379.45,435.0,459,495.0,545,620,665.0,840,865.0,1015.0,1230,1628,2105};
    const prc_cmp wvl_max_sensor_nm[]={331.7,340.16,360.90,380.45,445.0,479,505.0,565,670,675.0,880,875.0,1025.0,1250,1654,2155};
    const long wvl_nbr_sensor_nm(sizeof(wvl_min_sensor_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_sensor_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_sensor_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_sensor_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_sensor_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end sensor-specific initialization
  }else if(grd_sng == "MODIS"){ // 
    /* MODIS sensor wavelength grid
       MODIS Level 1B land bands 1, 2, 3, 4, 5, 6, and 7 (centered at 648 nm,
       858 nm, 470 nm, 555 nm, 1240 nm, 1640 nm, and 2130 nm, respectively). 
       MODIS band ranges in um (half power bandwidth) from
       http://stratus.ssec.wisc.edu/streamer_web/userman/bandweights.html
       Spectral response functions themselves are at:
       http://daac.gsfc.nasa.gov/MODIS/FAQ/A_sci_spectral_response.shtml
       1 0.62  - 0.67
       2 0.84  - 0.88
       3 0.459 - 0.479
       4 0.545 - 0.565
       5 1.230 - 1.250
       6 1.628 - 1.654
       7 2.105 - 2.155 */
    const prc_cmp wvl_min_MODIS_nm[]={620,840,459,545,1230,1628,2105};
    const prc_cmp wvl_max_MODIS_nm[]={670,880,479,565,1250,1654,2155};
    const long wvl_nbr_MODIS_nm(sizeof(wvl_min_MODIS_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_MODIS_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_MODIS_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_MODIS_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_MODIS_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end MODIS-specific initialization
  }else if(grd_sng == "TOMS"){ // 
    /* TOMS sensor wavelength grid
       According to Masaru Yoshioka, TOMS sensors are centered at
       0.33120 0.33966 0.36040 0.37995 microns with 1 nm bandpass. 
       (Somewhere I recall reading that the TOMS sensor bandpass is 10 nm, oh well)
       0.55 and 0.63 are diagnostic with 10 nm bandpass */
    const prc_cmp wvl_min_TOMS_nm[]={330.7,339.16,359.90,379.45,545.0,625.0};
    const prc_cmp wvl_max_TOMS_nm[]={331.7,340.16,360.90,380.45,555.0,635.0};
    const long wvl_nbr_TOMS_nm(sizeof(wvl_min_TOMS_nm)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_TOMS_nm){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_TOMS_nm; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=wvl_min_TOMS_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
      wvl_max[wvl_idx]=wvl_max_TOMS_nm[wvl_idx]*1.0e-9; // [nm] -> [m]
    } // end for
    // end TOMS-specific initialization
  }else if(grd_sng == "SWNB"){ // 
    /* Define wavelength grid approximately equal to original SWNB grid
       1590 10 cm-1 bands from 2000--17900 cm-1
       100 bands from 176--558.65 nm
       Algorithm based on rt_utl.F:wvl_grd_mk() */
    // Generate wavenumber grid first (as in SWNB) then reverse
    const prc_cmp wvl_grd_min(175.4e-09); // [m]
    const prc_cmp wvl_grd_max(5.0e-06); // [m]
    const long nbr_pure_O3_bnd(100); // [nbr]
    prc_cmp wvl_grd_rsn; // [m]
    prc_cmp wvn_max_lcl; // [cm-1]
    const long wvl_nbr_SWNB(1690); // [nbr]
    if(wvl_nbr != wvl_nbr_SWNB){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_SWNB; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=1;wvl_idx<=wvl_nbr-nbr_pure_O3_bnd;wvl_idx++){
      wvn_max_lcl=2000.0+wvl_idx*10; // [cm-1]
      wvl_min[wvl_idx-1]=1.0/(100.0*wvn_max_lcl); // [m]
    } // end loop over wvl
    wvl_grd_rsn=(wvl_min[wvl_nbr-nbr_pure_O3_bnd-1]-wvl_grd_min)/nbr_pure_O3_bnd;
    for(wvl_idx=wvl_nbr;wvl_idx>=wvl_nbr-(nbr_pure_O3_bnd-1);wvl_idx--){
      wvl_min[wvl_idx-1]=wvl_grd_min+wvl_grd_rsn*(wvl_nbr-wvl_idx); // [m]
    } // end loop over wvl
    // Reverse monotonicity so output increases with wavelength not wavenumber
    rvr_vec(wvl_min,wvl_nbr);
    for(wvl_idx=0;wvl_idx<wvl_nbr-1;wvl_idx++){ // Loop ends at wvl_nbr-1
      wvl_max[wvl_idx]=wvl_min[wvl_idx+1]; // [m]
    } // end loop over wvl
    wvl_max[wvl_nbr-1]=wvl_grd_max; // [m]
    // end SWNB-specific initialization
  }else if(grd_sng == "RRTM_LW"){ // 
    /* RRTM spectral grid from Rachel Scanza 20120323      
       Default dust input files for CAM4 are on GLADE:
       /glade/proj3/cseg/inputdata/atm/cam/physprops/dust[1-4]_camrt_c080918.nc
       Default dust input file  for CAM5 is  on GLADE:
       /glade/proj3/cseg/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc */
    const prc_cmp wvn_min_RRTM_LW[]={  10.0, 350.0, 500.0, 630.0, 700.0, 820.0, 980.0,1080.0,
				     1180.0,1390.0,1480.0,1800.0,2080.0,2250.0,2390.0,2600.0}; // [cm-1]
    const prc_cmp wvn_max_RRTM_LW[]={ 350.0, 500.0, 630.0, 700.0, 820.0, 980.0,1080.0,1180.0,
				     1390.0,1480.0,1800.0,2080.0,2250.0,2390.0,2600.0,3250.0}; // [cm-1]
    const long wvl_nbr_RRTM_LW(sizeof(wvn_min_RRTM_LW)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_RRTM_LW){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_RRTM_LW; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=0.01/wvn_max_RRTM_LW[wvl_idx]; // [cm-1] -> [m]
      wvl_max[wvl_idx]=0.01/wvn_min_RRTM_LW[wvl_idx]; // [cm-1] -> [m]
    } // end loop over wvl_idx
    // end RRTM_LW-specific initialization
  }else if(grd_sng == "RRTM_SW"){ // 
    /* RRTM spectral grid from Rachel Scanza 20120323
       NB: RRTM_SW increases monotonically with wavenumber except for last bin */
    const prc_cmp wvn_min_RRTM_SW[]=
      { 2600.0, 3250.0, 4000.0, 4650.0, 5150.0,6150.0,7700.0, 8050.0,
       12850.0,16000.0,22650.0,29000.0,38000.0, 820.0}; // [cm-1]
    const prc_cmp wvn_max_RRTM_SW[]=
      { 3250.0, 4000.0, 4650.0, 5150.0, 6150.0,7700.0,8050.0,12850.0,
       16000.0,22650.0,29000.0,38000.0,50000.0,2600.0}; // [cm-1]
    const long wvl_nbr_RRTM_SW(sizeof(wvn_min_RRTM_SW)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_RRTM_SW){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_RRTM_SW; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=0.01/wvn_max_RRTM_SW[wvl_idx]; // [cm-1] -> [m]
      wvl_max[wvl_idx]=0.01/wvn_min_RRTM_SW[wvl_idx]; // [cm-1] -> [m]
    } // end loop over wvl_idx
    // end RRTM_SW-specific initialization
  }else if(grd_sng == "CAM_LW"){ // 
    /* Following bins are seven trace gas overlap bands from KHB96:
       const prc_cmp wvn_min_CAM_LW[]={500.0,750.0,820.0,880.0,900.0,1000.0,1120.0,1170.0};
       const prc_cmp wvn_max_CAM_LW[]={750.0,820.0,880.0,900.0,1000.0,1120.0,1170.0,1500.0};
       Following bins are the six dominant H2O bands (500--800 cm-1 H2O rotation-H2O continuum overlap and 800--1200 H2O window are split in half as per RaD86)
       Non-standard modifications:
       To capture all LW absorption, may wish to extend CAM_LW bounds to 200, 2200 cm-1:
       Minimum wavenumber of 200 cm-1 (50 um), upper bound of index of refraction data
       Maximum wavenumber to 2200 cm-1 (4.5 um), lower bound of index of refraction data */
    /* 19990724: Change CAM_LW grid band 1 from 250--500 cm-1 to 200--800 cm-1
       const prc_cmp wvn_max_CAM_LW[]={500.0,650.0,800.0,1000.0,1200.0,2000.0}; // [cm-1]
       20050130: Change CAM_LW grid band 1 from 200--800 cm-1 to 250--500 cm-1 
       Change back because CAM_LW is too psychotic to code generically for
       const prc_cmp wvn_min_CAM_LW[]={200.0,500.0,650.0,800.0,1000.0,1200.0}; // [cm-1]
       const prc_cmp wvn_max_CAM_LW[]={800.0,650.0,800.0,1000.0,1200.0,2200.0}; // [cm-1]

       20050320: Change from CCM_LW to CAM_LW grid based on CCM3 H2SO4 code
       CAM_LW bands are documented in CRB04 p. 120 Table 4.2 and in 
       CAM3:models/atm/cam/src/physics/cam1/volcrad.F90 
       CAM3:models/atm/cam/src/physics/cam1/aerosol_radiation_interface.F90 
       CAM3 code uses first two bands for H2O non-window and window 
       CAM3 band 1 is 0--800+1200--2200 cm-1
       CAM3 band 2 is 800--1200 cm-1
       Five last CAM3 bands (3,4,5,6,7) correspond to bands 3-7 below */
    const prc_cmp wvn_min_CAM_LW[]={200.0, 800.0,500.0,650.0, 800.0,1000.0,1200.0}; // [cm-1]
    const prc_cmp wvn_max_CAM_LW[]={800.0,1200.0,650.0,800.0,1000.0,1200.0,2000.0}; // [cm-1]
    const long wvl_nbr_CAM_LW(sizeof(wvn_min_CAM_LW)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_CAM_LW){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_CAM_LW; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=0.01/wvn_max_CAM_LW[wvl_idx]; // [cm-1] -> [m]
      wvl_max[wvl_idx]=0.01/wvn_min_CAM_LW[wvl_idx]; // [cm-1] -> [m]
    } // end loop over wvl_idx
    // end CAM_LW-specific initialization
  }else if(grd_sng == "GSFC_LW"){ // 
    /* Following 10 bands are from ChS94:
       ChS94 subdivide band 3 (540--800 cm-1) into three subbands: 540--620, 620--720, 720-800 cm-1
       [fnc] Set minimum wavenumber to 200 cm-1 (50 um), upper bound of index of refraction data
       [fnc] Set maximum wavenumber to 3000 cm-1 (3.3 um), lower bound of index of refraction data */
    const prc_cmp wvn_min_GSFC_LW[]={200.0,340.0,540.0,800.0, 980.0,1100.0,1215.0,1380.0,1900.0,540.0}; // [cm-1]
    const prc_cmp wvn_max_GSFC_LW[]={340.0,540.0,800.0,980.0,1100.0,1215.0,1380.0,1900.0,3000.0,620.0}; // [cm-1]
    const long wvl_nbr_GSFC_LW(sizeof(wvn_min_GSFC_LW)/sizeof(prc_cmp));
    if(wvl_nbr != wvl_nbr_GSFC_LW){
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      wvl_nbr=wvl_nbr_GSFC_LW; // [nbr]
    } // endif
    rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      wvl_min[wvl_idx]=0.01/wvn_max_GSFC_LW[wvl_idx]; // [cm-1] -> [m]
      wvl_max[wvl_idx]=0.01/wvn_min_GSFC_LW[wvl_idx]; // [cm-1] -> [m]
    } // end loop over wvl_idx
    // end GSFC_LW-specific initialization
  }else{
    err_prn(sbr_nm,grd_sng+" is not a valid grid type");
  } // end else

  // wvl_min and wvl_max arrays are defined, derive other properties from these 
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    wvl_dlt[wvl_idx]=wvl_max[wvl_idx]-wvl_min[wvl_idx]; // [m]
    wvl_ctr[wvl_idx]=0.5*(wvl_max[wvl_idx]+wvl_min[wvl_idx]); // [m]
  } // end loop over wvl_idx
  bool ncr_wvl; // [flg] Wavelength array increases
  if(grd_sng == "CAM_LW"){
    // fxm: CAM_LW decreases (sort of) but not monotonically so this would cause mnt_ncr() to fail
    ncr_wvl=false;
  }else if(grd_sng == "CAM_SW"){
    // fxm: CAM_SW increases (sort of) but not monotonically so this would cause mnt_ncr() to fail
    ncr_wvl=true;
  }else if(grd_sng == "RRTM_SW"){
    // fxm: RRTM_SW decreases monotonically except for last bin
    ncr_wvl=false;
  }else ncr_wvl=mnt_ncr(wvl_min,wvl_nbr);
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    if(ncr_wvl) wvl_grd[wvl_idx]=wvl_min[wvl_idx]; else wvl_grd[wvl_idx]=wvl_max[wvl_idx]; // [m]
  } // end loop over wvl_idx
  if(ncr_wvl) wvl_grd[wvl_nbr]=wvl_max[wvl_nbr-1]; else wvl_grd[wvl_nbr]=wvl_min[wvl_nbr-1]; // [m]
  wvl_mnm=vec_min(wvl_min,wvl_nbr); // [m]
  wvl_mxm=vec_max(wvl_max,wvl_nbr); // [m]

  // At this point, wavelength grid is defined for all grid-types
  // If requested, derive wavenumber and/or frequency grids too

  // fxm: move to flg_wvn_grd_set()
  if(flg_wvn_grd){
    // Recompute wavenumber grid
    long wvn_idx; // [idx] Counting index for wvn
    wvn_nbr=wvl_nbr;
    for(wvn_idx=0;wvn_idx<wvn_nbr;wvn_idx++){
      wvn_min[wvn_idx]=1.0/(100.0*wvl_max[wvn_idx]); // [m] -> [cm-1]
      wvn_max[wvn_idx]=1.0/(100.0*wvl_min[wvn_idx]); // [m] -> [cm-1]
      wvn_ctr[wvn_idx]=0.5*(wvn_min[wvn_idx]+wvn_max[wvn_idx]); // [cm-1]
      wvn_dlt[wvn_idx]=wvn_max[wvn_idx]-wvn_min[wvn_idx]; // [cm-1]
      wvn_grd[wvn_idx]=1.0/(100.0*wvl_grd[wvn_idx]); // [m] -> [cm-1]
    } // end loop over wvn 
    wvn_grd[wvn_nbr]=1.0/(100.0*wvl_grd[wvl_nbr]); // [m] -> [cm-1]
    wvn_mnm=vec_min(wvn_min,wvn_nbr); // [cm-1]
    wvn_mxm=vec_max(wvn_max,wvn_nbr); // [cm-1]
  } // endif flg_wvn_grd
  
  // fxm: move to flg_frq_grd_set()
  if(flg_frq_grd){
    // Recompute frequency grid
    using phc::speed_of_light; // (2.99793e+08) [m s-1] Speed of light in vacuo
    long frq_idx; // [idx] Counting index for frq
    frq_nbr=wvl_nbr;
    for(frq_idx=0;frq_idx<frq_nbr;frq_idx++){
      frq_min[frq_idx]=speed_of_light/wvl_max[frq_idx]; // [m] -> [Hz]
      frq_max[frq_idx]=speed_of_light/wvl_min[frq_idx]; // [m] -> [Hz]
      frq_ctr[frq_idx]=0.5*(frq_min[frq_idx]+frq_max[frq_idx]); // [Hz]
      frq_dlt[frq_idx]=frq_max[frq_idx]-frq_min[frq_idx]; // [Hz]
      frq_grd[frq_idx]=speed_of_light/wvl_grd[frq_idx]; // [m] -> [Hz]
    } // end loop over frq 
    frq_grd[frq_nbr]=speed_of_light/wvl_grd[wvl_nbr]; // [m] -> [Hz]
    frq_mnm=vec_min(frq_min,frq_nbr); // [Hz]
    frq_mxm=vec_max(frq_max,frq_nbr); // [Hz]
  } // endif flg_frq_grd

  // Now that grid is know, interpolate to find debuggin wavelength
  wvl_idx_dbg=vec_val2idx(wvl_ctr,wvl_nbr,wvl_dbg_dfl); // [idx] Wavelength bin for debugging

  return rcd; // [enm] Return success code
} // end wvl_grd_cls::recompute()

// Private member functions end
// Global functions with C++ linkages begin

// Global functions with C++ linkages end
