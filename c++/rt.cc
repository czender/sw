// $Id$ 

// Purpose: Implementation (declaration) of radiative transfer classes

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <rt.hh> // Radiative transfer

// rt_cls class

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const rt_cls &rt_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for rt_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     std::cout << rt_obj; */
  size_t dyn_sz(0); // [sz] Total size of dynamic buffers owned by object

  srm_out << "Public contents of rt_obj: " << std::endl;
  srm_out << "rt_obj.nst_nbr_get() = " << rt_obj.nst_nbr_get() << std::endl;
  srm_out << "rt_obj.typ_get() = " << rt_obj.typ_get() << std::endl;
  srm_out << "rt_obj.dsc_get() = " << rt_obj.dsc_get() << std::endl;
  srm_out << "rt_obj.plrmu_avg_get() = " << rt_obj.plrmu_avg_get() << std::endl;
  srm_out << "rt_obj.rt_sln_fnc_ptr_get() = " << rt_obj.rt_sln_fnc_ptr_get() << std::endl;
  srm_out << "rt_obj.wvl_nbr_get() = " << rt_obj.wvl_nbr_get() << std::endl;
  srm_out << "Private contents of rt_obj: " << std::endl;
  srm_out << "plrmu_avg = " << rt_obj.plrmu_avg << std::endl;
  if(rt_obj.flg_rta_in_dyn_mmr){
    srm_out << "rfl_flx[0] = " << rt_obj.rfl_flx[0] << std::endl;
    srm_out << "trn_flx[0] = " << rt_obj.trn_flx[0] << std::endl;
    srm_out << "abs_flx[0] = " << rt_obj.abs_flx[0] << std::endl;
    srm_out << "tau_ext[0] = " << rt_obj.tau_ext[0] << std::endl;
    srm_out << "ss_alb[0] = " << rt_obj.ss_alb[0] << std::endl;
    srm_out << "asm_prm[0] = " << rt_obj.asm_prm[0] << std::endl;
    // fxm: is this totally bogus!! is sizeof just returning size of pointer not array?
    dyn_sz+=sizeof(rt_obj.rfl_flx)+sizeof(rt_obj.trn_flx)+sizeof(rt_obj.abs_flx)+sizeof(rt_obj.tau_ext)+sizeof(rt_obj.ss_alb)+sizeof(rt_obj.asm_prm); // [sz] Total size of dynamic buffers owned by object
  }else{
    srm_out << "No layer optical properties in dynamic memory" << std::endl;
  } // !rt_obj.flg_rta_in_dyn_mmr
  srm_out << "sizeof(rt_obj) = " << sizeof(rt_obj) << " B" << std::endl;
  srm_out << "Size of dynamic arrays = " << dyn_sz << " B" << std::endl;
  srm_out << "Total memory owned by object = " << sizeof(rt_obj)+dyn_sz << " B" << std::endl;
    
  return srm_out; // [srm] Reference to output stream for cascading
} // !operator<<()

// Friendly functions end
// Static members begin

int rt_cls::nst_nbr(0); // [nbr] Number of instantiated class members
int rt_cls::wrn_nbr(0); // [nbr] Number of accumulated warnings
prc_cmp rt_cls::plrmu_avg_dfl(0.5); // [frc] Mean zenith angle, default
std::string rt_cls::rt_typ_dfl("two_srm_iso_sct"); // [sng] Radiative transfer method, default
long rt_cls::wvl_nbr_dfl(-1); // [nbr] Number of wavelength bands, default
const long rt_cls::wvl_nbr_max(100000); // [nbr] Maximum number of wavelengths

// Static members end
// Static member functions begin

int // O [enm] Return success code
rt_cls::tst(long obj_nbr){ // [fnc] Self-test of rt_cls class
  // Purpose: Perform self-test of rt_cls class
  int rcd(0); // [enm] Return success code
  std::string sbr_nm("rt_cls::tst"); // [sng] Subroutine name
  long wvl_nbr_foo(10000); // [nbr] Number of wavelength bands

  // Test for memory leaks
  std::cout << sbr_nm+"()"+" is testing for memory leaks by creating and destroying " << obj_nbr << " rt_cls objects..." << std::endl;
  long idx; // [idx] Counting index
  rt_cls *tst_obj; // [sct] Test object
  prc_cmp plrmu_avg_tmp; // [frc] Mean zenith angle
  for(idx=0;idx<obj_nbr;idx++){
    tst_obj=new rt_cls; // [sct] Test object
    tst_obj->wvl_nbr_set(wvl_nbr_foo); // [nbr] Number of wavelength bands
    plrmu_avg_tmp=tst_obj->plrmu_avg_get(); // [frc] Mean zenith angle
    std::cout << "idx = " << idx << ", nst_nbr = " << nst_nbr << ", plrmu_avg = " << plrmu_avg_tmp << std::endl;
    delete tst_obj; // [sct] Test object
  } // !idx

  // Test all radiative transfer methods
  tst_obj=new rt_cls; // [sct] Test object
  std::cout << sbr_nm+"()"+" sizeof(rt_cls) object is " << sizeof tst_obj << " B and sizeof(rt_cls) type is " << sizeof(rt_cls) << " B" << std::endl;
  plrmu_avg_tmp=tst_obj->plrmu_avg_get(); // [frc] Mean zenith angle
  std::cout << "abb = " << tst_obj->abb << ", plrmu_avg_tmp = " << plrmu_avg_tmp << std::endl;
  delete tst_obj; // [sct] Test object

  return rcd; // [enm] Return success code
} // !rt_cls::tst()

int rt_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members
sng2sng_map rt_cls::opt2abb_map=rt_cls::opt2abb_map_mk(); // [map] Option to abbreviation map
sng2sng_map rt_cls::opt2abb_map_mk(){ // [fnc] Create option to abbreviation map
  sng2sng_map map_tmp; 

  map_tmp.insert(sng2sng_map::value_type("two_srm_iso_sct","two_srm_iso_sct"));

  map_tmp.insert(sng2sng_map::value_type("two_srm_asm_sct","two_srm_asm_sct"));

  return map_tmp; // [map] Option to abbreviation map
} // !rt_cls::opt2abb_map_mk()

std::string // O [sng] Radiative transfer method abbreviation
rt_cls::opt2abb // [fnc] Option to abbreviation mapper
(const std::string &opt_sng){ // [sng] Shorthand option
  // Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("rt_cls::opt2abb",opt_sng+" is unknown");
  return itr->second; // [sng] Radiative transfer method abbreviation
} // !rt_cls::opt2abb()

const sng2rt_sct_map rt_cls::rt_map=rt_cls::rt_map_mk(); // [map] Radiative transfer method map
const sng2rt_sct_map rt_cls::rt_map_mk(){ // [fnc] Create radiative transfer method map
  const rt_sct rt[]={
    {"two_srm_iso_sct", // [sng] Radiative transfer method abbreviation
     "Two stream isotropically scattering", // [sng] Radiative transfer method description
     two_srm_iso_sct}, // [fnc] Function to compute layer optical properties
    {"two_srm_asm_sct", // [sng] Radiative transfer method abbreviation
     "Two stream anisotropically scattering", // [sng] Radiative transfer method description
     two_srm_asm_sct} // [fnc] Function to compute layer optical properties
  }; // !rt_sct rt[]
  long idx; // [idx] Counting index
  int rt_nbr=sizeof(rt)/sizeof(rt_sct); // [nbr] Number of radiative transfer methods
  sng2rt_sct_map rt_map_tmp; // [sct] Map with key=abbreviation, value=radiative transfer method structure
  for(idx=0;idx<rt_nbr;idx++){
    /* fxm: Define variables before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map. */
    rt_map_tmp.insert(sng2rt_sct_map::value_type(rt[idx].abb,rt[idx])); // [sct] Map with key=abbreviation, value=radiative transfer method structure
  } // !idx
  // Return a value to initialize a static class member
  return rt_map_tmp;
} // !rt_cls::rt_map_mk()

sng2var_mtd_sct_map rt_cls::var_mtd_map; // [map] Variable metadata map
sng2var_mtd_sct_map rt_cls::var_mtd_map_mk // [fnc] Create variable metadata map
(const int nc_out) // I [fl] netCDF file for output 
{
  /* Purpose: Create variable metadata map
     Map is used by class members wishing to output their contents
     Routine must not be called until dimensions are defined in output file
     Usage:
     Once nc_out is defined, call
     rt_cls::var_mtd_map=var_mtd_map_mk(nc_out); // [map] Variable metadata map
  */

  const std::string sbr_nm("rt_cls::var_mtd_map_mk"); // [sng] Subroutine name
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // Wavelength dimension
  const int *dmn_wvl(&wvl_dmn); // Pointer to wavelength dimension
  //  const int *dmn_scl(&wvl_dmn); // Dummy argument, not used

  var_mtd_sct var_mtd[]={
    {0,"abs_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux absorptance","units","fraction"},
    {0,"rfl_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux reflectance","units","fraction"},
    {0,"trn_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux transmittance","units","fraction"}
  }; // !var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  long idx; // [idx] Counting index
  sng2var_mtd_sct_map var_mtd_map_tmp; // [sct] Map with key=abbreviation, value=variable metadata structure
  for(idx=0;idx<var_mtd_nbr;idx++){
    var_mtd_map_tmp.insert(sng2var_mtd_sct_map::value_type(var_mtd[idx].nm,var_mtd[idx])); // [sct] Map with key=abbreviation, value=variable metadata structure
  } // !idx
  // Return value to initialize static class member
  return var_mtd_map_tmp;
} // !rt_cls::var_mtd_map_mk()

var_mtd_sct // [sct] Variable metadata structure
rt_cls::var2mtd // [fnc] Variable name to metadata mapper
(const std::string &var_nm){ // [sng] Variable name
  // Purpose: Return metadata structure for requested variable
  return var_mtd_map.find(var_nm)->second; // [sct] Variable metadata structure
} // !rt_cls::var2mtd()

std::string // [sng] Description string
rt_cls::abb2dsc(const std::string &abb_sng){ // [fnc] Abbreviation to description mapper
  // Purpose: Return description of radiative transfer method
  return rt_map.find(abb_sng)->second.dsc; // [sng] Description string
} // !rt_cls::abb2dsc()

rt_sln_fnc_ptr_typ // [fnc] Pointer to function returning layer optical properties
rt_cls::abb2fnc // [fnc] Abbreviation to function mapper
(const std::string &abb_sng){ // [sng] Radiative transfer method abbreviation
  // Purpose: Return pointer to function which computes layer optical properties
  return rt_map.find(abb_sng)->second.rt_sln_fnc_ptr; // [fnc] Pointer to function returning layer optical properties
} // !rt_cls::abb2fnc()

// Static member functions end
// Public member functions begin

rt_cls::rt_cls // [fnc] Constructor
(const std::string &rt_arg, // I [sng] Radiative transfer method abbreviation
 const long &wvl_nbr_arg, // I [nbr] Number of wavelength bands
 const prc_cmp *asm_prm_arg, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb_arg, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext_arg) // I [frc] Extinction optical depth
{
  /* Purpose: Constructor for rt_cls
     Usage: 
     rt_cls rt("two_srm_iso_sct"); // [rt] Radiative transfer object
  */

  rcd=0; // [enm] Return success code
  flg_rta_in_dyn_mmr=false; // [flg] Reflectance/transmittance/absorptance are in dynamic memory

  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls

  // Set type before mean zenith angle
  rcd+=typ_set(rt_arg); // [sng] Set radiative transfer method

  // Constructor defaults to reasonable value to allow for calls that omit plrmu_avg
  rcd+=plrmu_avg_set(plrmu_avg_dfl); // [frc] Mean zenith angle

  // rt_cls object is worthless without single scattering properties
  if(wvl_nbr_arg > wvl_nbr_dfl){
    // Set single scattering properties with arguments to constructor
    rcd+=twg_set // [fnc] Set slab single scattering properties
      (wvl_nbr_arg, // I [nbr] Number of wavelength bands
       asm_prm_arg, // I [frc] Asymmetry parameter
       ss_alb_arg, // I [frc] Single scattering albedo
       tau_ext_arg); // I [frc] Extinction optical depth
  }else{
    // Set single scattering properties with useful defaults
    // Must define temporary "_tmp" variables because "_arg" variables are const
    const long wvl_nbr_tmp(1); // [nbr] Number of wavelength bands
    const prc_cmp asm_prm_tmp(0.8); // [frc] Asymmetry parameter
    const prc_cmp ss_alb_tmp(0.8); // [frc] Single scattering albedo
    const prc_cmp tau_ext_tmp(1.0); // [frc] Extinction optical depth
    rcd+=twg_set // [fnc] Set slab single scattering properties
      (wvl_nbr_tmp, // I [nbr] Number of wavelength bands
       &asm_prm_tmp, // I [frc] Asymmetry parameter
       &ss_alb_tmp, // I [frc] Single scattering albedo
       &tau_ext_tmp); // I [frc] Extinction optical depth
  } // endif 

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  nst_nbr++; // [nbr] Number of instantiated class members
} // !rt_cls constructor

rt_cls::~rt_cls(){ // [fnc] Destructor
  rcd+=deallocate(); // [fnc] Free dynamic memory for object
  nst_nbr--; // [nbr] Number of instantiated class members
} // !rt_cls destructor

void rt_cls::prn()const{std::cout << this;} // [fnc] Print object contents

const std::string rt_cls::dsc_get()const{return dsc;} // [sng] Radiative transfer method description
const std::string rt_cls::typ_get()const{return abb;} // [sng] Radiative transfer method abbreviation
prc_cmp rt_cls::plrmu_avg_get()const{return plrmu_avg;} // [frc] Mean zenith angle
rt_sln_fnc_ptr_typ rt_cls::rt_sln_fnc_ptr_get()const{return rt_sln_fnc_ptr;} // [fnc] Function to compute layer optical properties

long rt_cls::wvl_nbr_get()const{return wvl_nbr;} // [nbr] Number of wavelength bins
const prc_cmp *rt_cls::rfl_flx_get()const{return rfl_flx;} // [frc] Flux reflectance
const prc_cmp *rt_cls::abs_flx_get()const{return abs_flx;} // [frc] Flux absorptance
const prc_cmp *rt_cls::trn_flx_get()const{return trn_flx;} // [frc] Flux transmittance
const prc_cmp *rt_cls::tau_ext_get()const{return tau_ext;} // [frc] Extinction optical depth
const prc_cmp *rt_cls::asm_prm_get()const{return asm_prm;} // [frc] Asymmetry parameter
const prc_cmp *rt_cls::ss_alb_get()const{return ss_alb;} // [frc] Single scattering albedo

int // [enm] Return success code
rt_cls::rt_sln_fnc_ptr_set(const rt_sln_fnc_ptr_typ &rt_sln_fnc_ptr_arg){ // [fnc] Function to compute layer optical properties
  std::string sbr_nm("rt_cls::rt_sln_fnc_ptr_set"); // [sng] Subroutine name
  rt_sln_fnc_ptr=rt_sln_fnc_ptr_arg; // [fnc] Function to compute layer optical properties
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !rt_cls::rt_sln_fnc_ptr_set()

int // [enm] Return success code
rt_cls::typ_set(const std::string &rt_arg){ // [fnc] Set radiative transfer method
  std::string sbr_nm("rt_cls::typ_set"); // [sng] Subroutine name
  abb= (rt_arg == "") ? rt_typ_dfl : rt_arg; // [sng] Radiative transfer method abbreviation
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !rt_cls::typ_set()

int // [enm] Return success code
rt_cls::wvl_nbr_set(const long &wvl_nbr_arg) // [nbr] Number of wavelength bands
{
  // Purpose: Set number of wavelength bands
  // Changing number of wavelength bands alone is no reason to recompute
  std::string sbr_nm("rt_cls::wvl_nbr_set"); // [sng] Subroutine name

  // Sanity check
  assert(wvl_nbr_arg > 0 && wvl_nbr_arg <= wvl_nbr_max); // [nbr] Allowed number of wavelength bands
  wvl_nbr=wvl_nbr_arg; // [nbr] Number of wavelength bands

  // Prepare space to store layer optical properties
  rcd+=allocate(); // [fnc] Allocate dynamic memory for object

  return rcd; // [enm] Return success code
} // !rt_cls::wvl_nbr_set()

int // [enm] Return success code
rt_cls::nc_out_set // [fl] Set netCDF output file
(const int &nc_out_arg) // [id] netCDF output file
{
  // Purpose: Set netCDF output file
  std::string sbr_nm("rt_cls::nc_out_set"); // [sng] Subroutine name

  //  nc_out=nc_out_arg; // [fl] netCDF file ID for output

  return rcd; // [enm] Return success code
} // !rt_cls::nc_out_set()
  
int // [enm] Return success code
rt_cls::nc_out_set // [fl] Set netCDF output file
(const std::string &fl_out) // [fl] netCDF output file
{
  /* Purpose: Set netCDF output file
     Routine is overloaded */
  std::string sbr_nm("rt_cls::nc_out_set"); // [sng] Subroutine name

  //  rcd=nco_create(fl_out,NC_CLOBBER,nc_out); 

  return rcd; // [enm] Return success code
} // !rt_cls::nc_out_set()
  
int // [enm] Return success code
rt_cls::plrmu_avg_set(const prc_cmp &plrmu_avg_arg) // [frc] Mean zenith angle
{
  // Purpose: Set hemispheric mean zenith angle for two stream approximation
  std::string sbr_nm("rt_cls::plrmu_avg_set"); // [sng] Subroutine name

  assert(plrmu_avg_arg > 0.0 && plrmu_avg_arg < 1.0); // [frc] Allowed zenith angle range
  plrmu_avg=plrmu_avg_arg; // [frc] Mean zenith angle
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !rt_cls::plrmu_avg_set()

int // O [enm] Return success code
rt_cls::twg_set // [fnc] Set slab single scattering properties
(const long &wvl_nbr_arg, // I [nbr] Number of wavelength bands
 const prc_cmp *asm_prm_arg, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb_arg, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext_arg) // I [frc] Extinction optical depth
{
  /* Purpose: Set slab single scattering properties
     asm_prm, ss_alb, and tau_ext are expected to be valid
     Note that rt_cls() expects calling routine to perform memory management for asm_prm,ss_alb,tau_ext */

  std::string sbr_nm("rt_cls::twg_set"); // [sng] Subroutine name

  // Free existing dynamic memory before changing dimension sizes
  if(flg_rta_in_dyn_mmr && wvl_nbr_arg != wvl_nbr) 
    rcd+=deallocate(); // [fnc] Free dynamic memory for object

  // Memory allocation depends on wvl_nbr
  rcd+=wvl_nbr_set(wvl_nbr_arg); // [nbr] Number of wavelength bands

  // Free this memory in rt_cls::deallocate()
  asm_prm=vec_cpy(asm_prm_arg,wvl_nbr_arg); // [frc] Asymmetry parameter
  ss_alb=vec_cpy(ss_alb_arg,wvl_nbr_arg); // [frc] Single scattering albedo
  tau_ext=vec_cpy(tau_ext_arg,wvl_nbr_arg); // [frc] Extinction optical depth

  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // !rt_cls::twg_set()

int // O [enm] Return success code
rt_cls::var_put // [fnc] Write variable to output file
(const int &nc_out, // I [fl] netCDF file for output 
 const std::string var_nm) // I [sng] Name of variable to output
{
  /* Purpose: Write specified variable to output file 
     Usage: rcd+=rt_obj.var_put("rfl_flx");
     Outputting a variable changes the file so file may not be declared const */
  // fxm: Not sure why but this is broken, and segfaults downstream in nco_var_dfn() because the dimension does not seem to be valid

  std::string sbr_nm("rt_cls::var_put"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  if(dbg_lvl_get() == dbg_crr){
    std::cout << "Size of wvl dimension on disk = " << nco_inq_dimlen(nc_out,static_cast<std::string>("wvl")) << ", Size in structure = " << wvl_nbr << std::endl;
  } // endif dbg

  // If metadata structure has not yet been defined, create it
  if(var_mtd_map.size() <= 0) var_mtd_map=var_mtd_map_mk(nc_out); // [map] Variable metadata map

  // Find metadata structure for this particular variable
  var_mtd_sct var_mtd(var2mtd(var_nm)); // [sct] Variable metadata structure
  const int var_mtd_nbr(1); // [nbr] Number of variables to output
  const int dmn_nbr_max(1); // [nbr] Maximum number of dimensions allowed in single variable in output file

  rcd+=nco_var_dfn // [fnc] Define variables in output netCDF file
    (nc_out, // I [fl] netCDF file for output 
     &var_mtd, // I/O [sct] Array of structures containing variable metadata
     var_mtd_nbr, // I [nbr] Number of variables in array
     dmn_nbr_max); // I [nbr] Maximum number of dimensions allowed in single variable in output file

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !rt_cls::var_put()

// Public member functions end
// Private member functions begin

int // [enm] Return success code
rt_cls::allocate(){ // [fnc] Allocate dynamic memory for object
  // Purpose: Allocate dynamic memory for object
  // Routine should only be called from twg_set() since wvl_nbr must be known
  
  // Make sure wvl_nbr has been initialized correctly
  assert(wvl_nbr > 0 && wvl_nbr < wvl_nbr_max);

  if(!flg_rta_in_dyn_mmr){
    // Reflectance/transmittance/absorptance de-allocated in rt_cls::deallocate()
    rfl_flx=new prc_cmp[wvl_nbr]; // [frc] Flux reflectance
    trn_flx=new prc_cmp[wvl_nbr]; // [frc] Flux transmittance
    abs_flx=new prc_cmp[wvl_nbr]; // [frc] Flux absorptance

    // Reset flag
    flg_rta_in_dyn_mmr=true; // [flg] Reflectance/transmittance/absorptance are in dynamic memory
  } // endif

  return rcd; // [enm] Return success code
} // !rt_cls_cls::allocate()

int // [enm] Return success code
rt_cls::deallocate(){ // [fnc] Free dynamic memory for object
  // Purpose: Free dynamic memory for object

  // Delete any arrays allocated by private member functions
  if(flg_rta_in_dyn_mmr){
    // Reflectance/transmittance/absorptance allocated in rt_cls::allocate()
    delete []rfl_flx; // [frc] Flux reflectance
    delete []trn_flx; // [frc] Flux transmittance
    delete []abs_flx; // [frc] Flux absorptance

    // Tau, omega, g allocated by vec_cpy() in rt_cls::twg_set()
    delete []asm_prm; // [frc] Asymmetry parameter
    delete []ss_alb; // [frc] Single scattering albedo
    delete []tau_ext; // [frc] Extinction optical depth

    // Reset flag
    flg_rta_in_dyn_mmr=false; // [flg] Reflectance/transmittance/absorptance are in dynamic memory
  } // endif

  return rcd; // [enm] Return success code
} // !rt_cls::deallocate()

int // [enm] Return success code
rt_cls::recompute(){ // [fnc] Recompute properties of object
  std::string sbr_nm("rt_cls::recompute"); // [sng] Subroutine name

  // Set private members which lack set() functions
  dsc=abb2dsc(abb); // [sng] Radiative transfer method description
  rt_sln_fnc_ptr=abb2fnc(abb); // [fnc] Pointer to function returning layer optical properties
  
  if(abb == "two_srn_asm_sct") err_prn(sbr_nm,"Anisotropic scattering option unavailable at this time."); 

  // Space for answers must have been allocated
  if(!flg_rta_in_dyn_mmr) err_prn(sbr_nm," attempt to recompute() before optical properties set and allocate()d"); 

  rcd+=(*rt_sln_fnc_ptr) // [fnc] Function to compute layer optical properties
    (wvl_nbr, // I [nbr] Number of wavelength bands
     plrmu_avg, // I [frc] Mean zenith angle
     asm_prm, // I [frc] Asymmetry parameter
     ss_alb, // I [frc] Single scattering albedo
     tau_ext, // I [frc] Extinction optical depth
     rfl_flx, // O [frc] Flux reflectance
     trn_flx, // O [frc] Flux transmittance
     abs_flx); // O [frc] Flux absorptance

  return rcd; // [enm] Return success code
} // !rt_cls::recompute()

// Private member functions end
// Global functions with C++ linkages begin

int // O [enm] Return success code
two_srm_iso_sct // [fnc] Two stream approximation for isotropically scattering slab
(const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const prc_cmp &plrmu_avg, // I [frc] Mean zenith angle
 const prc_cmp *asm_prm, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext, // I [frc] Extinction optical depth
 prc_cmp *rfl_flx, // O [frc] Flux reflectance
 prc_cmp *trn_flx, // O [frc] Flux transmittance
 prc_cmp *abs_flx) // O [frc] Flux absorptance
{
  /* Purpose: Two stream approximation for isotropically scattering slab 
     Handles conservative and non-conservative scattering cases correctly
     Is robust for large optical depths */

  // Output
  const long wrn_nbr_max(100); // [nbr] Maximum number of warning messages to print
  int rcd(0); // [enm] Return success code
  std::string sbr_nm("two_srm_iso_sct"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  long wvl_idx; // [idx] Counting index for wvl
  prc_cmp *rho_inf=new prc_cmp[wvl_nbr]; // [frc] Reflectance at infinite optical depth
  // Local variables
  double d_str; // [frc] D star
  double gmm; // [frc] Gamma
  double gmm_sqr; // [frc] Gamma squared
  double exp_mns_2gmm_tau; // [frc] exp(-2*gmm*tau_ext)
  double exp_mns_gmm_tau; // [frc] exp(-gmm*tau_ext)
  prc_cmp ss_coalb_sqrt; // [frc] Square root of co-albedo
#ifdef PRC_FLT
  const prc_cmp ss_alb_bnd(1.0-1.0e-6); // [frc] Maximum single scattering albedo treated as non-conservative
#else // !PRC_FLT
  const prc_cmp ss_alb_bnd(1.0-1.0e-12); // [frc] Maximum single scattering albedo treated as non-conservative
#endif // !PRC_FLT

  assert(plrmu_avg != 0.0);
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    assert(ss_alb[wvl_idx] <= 1.0);
    if(ss_alb[wvl_idx] > ss_alb_bnd){
      if(ss_alb[wvl_idx] != 1.0){
	if(rt_cls::wrn_nbr < wrn_nbr_max) std::cout << "WARNING: "+sbr_nm+"() using conservative scattering approximation for ss_alb = " << ss_alb[wvl_idx] << " at wvl_idx = " << wvl_idx << std::endl;
	if(rt_cls::wrn_nbr == wrn_nbr_max) std::cout << "WARNING: "+sbr_nm+"() Reached " << wrn_nbr_max << " warning messages. No more warnings will be printed" << std::endl;
	// fxm: modifying a static class member by a non-member function is NG
	rt_cls::wrn_nbr++; // [nbr] Number of accumulated warnings
      } // endelse
      rfl_flx[wvl_idx]=tau_ext[wvl_idx]/(1.0+tau_ext[wvl_idx]); // [frc] Flux reflectance ThS99 p. 233 (7.53)
      trn_flx[wvl_idx]=1.0/(1.0+tau_ext[wvl_idx]); // [frc] Flux transmittance ThS99 p. 233 (7.53)
    }else{
      gmm=std::sqrt(1.0-ss_alb[wvl_idx])/plrmu_avg; // [frc] Gamma ThS99 p. 228 (7.30)
      gmm_sqr=gmm*gmm; // [frc] Gamma squared
      ss_coalb_sqrt=std::sqrt(1.0-ss_alb[wvl_idx]); // [frc] Square root of co-albedo
      rho_inf[wvl_idx]=(1.0-ss_coalb_sqrt)/(1.0+ss_coalb_sqrt); // [frc] Reflectance at infinite optical depth ThS99 p. 229 (7.33)
      exp_mns_2gmm_tau=std::exp(-2.0*gmm*tau_ext[wvl_idx]); // [frc] exp(-2*gmm*tau_ext)
      exp_mns_gmm_tau=std::exp(-gmm*tau_ext[wvl_idx]); // [frc] exp(-gmm*tau_ext)
      d_str=1.0-rho_inf[wvl_idx]*rho_inf[wvl_idx]*exp_mns_2gmm_tau; // [frc] D star ThS99 p. 229 (7.38)
      rfl_flx[wvl_idx]=rho_inf[wvl_idx]*(1.0-exp_mns_2gmm_tau)/d_str; // [frc] Flux reflectance ThS99 p. 230 (7.42)
      trn_flx[wvl_idx]=(1.0-rho_inf[wvl_idx]*rho_inf[wvl_idx])*exp_mns_gmm_tau/d_str; // [frc] Flux transmittance ThS99 p. 230 (7.43)
    } // !ss_alb
    abs_flx[wvl_idx]=1.0-rfl_flx[wvl_idx]-trn_flx[wvl_idx]; // [frc] Flux absorptance ThS99 p. 230 (7.44)
  } // !wvl_idx
  
  // Perhaps rho_inf is worth saving?
  delete []rho_inf; // [frc] Reflectance at infinite optical depth

  if(rt_cls::wrn_nbr >= wrn_nbr_max) std::cout << "WARNING: "+sbr_nm+"() is exited with " << rt_cls::wrn_nbr << " accumulated warnings, " << rt_cls::wrn_nbr-wrn_nbr_max << " of which were not printed" << std::endl;

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !two_srm_iso_sct()

int // O [enm] Return success code
two_srm_asm_sct // [fnc] Two stream approximation for anisotropically scattering slab
(const long &wvl_nbr, // I [nbr] Number of wavelength bands
 const prc_cmp &plrmu_avg, // I [frc] Mean zenith angle
 const prc_cmp *asm_prm, // I [frc] Asymmetry parameter
 const prc_cmp *ss_alb, // I [frc] Single scattering albedo
 const prc_cmp *tau_ext, // I [frc] Extinction optical depth
 prc_cmp *rfl_flx, // O [frc] Flux reflectance
 prc_cmp *trn_flx, // O [frc] Flux transmittance
 prc_cmp *abs_flx) // O [frc] Flux absorptance
{
  /* Purpose: Two stream approximation for anisotropically scattering slab */

  // Output
  int rcd(0); // [enm] Return success code
  std::string sbr_nm("two_srm_asm_sct"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  long wvl_idx; // [idx] Counting index for wvl
  prc_cmp *rho_inf=new prc_cmp[wvl_nbr]; // [frc] Reflectance at infinite optical depth
  // Local variables
  double d_str; // [frc] D star
  double gmm; // [frc] Gamma
  double gmm_sqr; // [frc] Gamma squared
  double exp_gmm_tau; // [frc] exp(gmm*tau_ext)
  double exp_mns_gmm_tau; // [frc] exp(-gmm*tau_ext)
  prc_cmp ss_coalb_sqrt; // [frc] Square root of co-albedo

  assert(plrmu_avg != 0.0);
  err_prn(sbr_nm,"Solution is same as isotropically scattering slab...");
  // fxm: exponentials are ill-conditioned for large optical depth, divide solution forms by exp(tau_ext) and rewrite
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    assert(ss_alb[wvl_idx] != 1.0);
    gmm=std::sqrt(1.0-ss_alb[wvl_idx])/plrmu_avg; // [frc] Gamma ThS99 p. 228 (7.30)
    gmm_sqr=gmm*gmm; // [frc] Gamma squared
    ss_coalb_sqrt=std::sqrt(1.0-ss_alb[wvl_idx]); // [frc] Square root of co-albedo
    rho_inf[wvl_idx]=(1.0-ss_coalb_sqrt)/(1.0+ss_coalb_sqrt); // [frc] Reflectance at infinite optical depth ThS99 p. 229 (7.33)
    exp_gmm_tau=std::exp(gmm*tau_ext[wvl_idx]); // [frc] exp(gmm*tau_ext)
    exp_mns_gmm_tau=std::exp(-gmm*tau_ext[wvl_idx]); // [frc] exp(-gmm*tau_ext)
    d_str=exp_gmm_tau-rho_inf[wvl_idx]*rho_inf[wvl_idx]*exp_mns_gmm_tau; // [frc] D star ThS99 p. 229 (7.38)
    rfl_flx[wvl_idx]=rho_inf[wvl_idx]*(exp_gmm_tau-exp_mns_gmm_tau)/d_str; // [frc] Flux reflectance ThS99 p. 230 (7.42)
    trn_flx[wvl_idx]=(1.0-rho_inf[wvl_idx]*rho_inf[wvl_idx])/d_str; // [frc] Flux transmittance ThS99 p. 230 (7.43)
    abs_flx[wvl_idx]=1.0-rfl_flx[wvl_idx]-trn_flx[wvl_idx]; // [frc] Flux absorptance ThS99 p. 230 (7.44)
  } // !wvl_idx

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !two_srm_asm_sct()

int // O [enm] Return success code
rt_lop // [fnc] Radiative transfer layer optical properties
(const int &nc_out, // I [fl] netCDF file for output 
 const rt_cls &rt_obj, // I [rt] Radiative transfer object
 const prc_cmp &slr_zen_ngl_cos) // I [frc] Cosine solar zenith angle
{
  /* Purpose: Determine layer optical properties
     nc_out is modified as variables are defined
     fxm: This definition causes problems with g++ on SUNMP and SGIMP64 architectures
     ncks -C -F -q -d wvl,0.5e-6 -v tau_ext,ss_alb,asm_prm,rfl_flx,trn_flx,abs_flx ${DATA}/mie/out.nc
  */

  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("rt_lop"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  const prc_cmp *tau_ext(rt_obj.tau_ext_get()); // [frc] Extinction optical depth
  const prc_cmp *asm_prm(rt_obj.asm_prm_get()); // [frc] Asymmetry parameter
  const prc_cmp *ss_alb(rt_obj.ss_alb_get()); // [frc] Single scattering albedo

  const prc_cmp *rfl_flx=rt_obj.rfl_flx_get(); // [frc] Flux reflectance
  const prc_cmp *trn_flx=rt_obj.trn_flx_get(); // [frc] Flux transmittance
  const prc_cmp *abs_flx=rt_obj.abs_flx_get(); // [frc] Flux absorptance

  const long wvl_nbr(rt_obj.wvl_nbr_get()); // [nbr] Number of wavelength bins

  for(idx=0;idx<wvl_nbr;idx++){
    ;
  } // !idx

  // Sanity check
  // Flux reflectance, transmittance and absorptance should sum to unity

  if(true){
    std::cout << "Radiative Transfer:" << std::endl;
    std::cout << "  Aerosol: Optical depth = " << tau_ext[0] << ", Single scattering albedo = " << ss_alb[0] << ", Asymmetry parameter = " << asm_prm[0] << std::endl;
    std::cout << "  Layer optical properties: Reflectance = " << rfl_flx[0] << ", trn_flx = " << trn_flx[0] << ", Absorptance = " << abs_flx[0] << std::endl;
    std::cout << "  Energy conservation: rfl_flx + trn_flx + abs_flx = " << rfl_flx[0]+trn_flx[0]+abs_flx[0] << std::endl;
  } // endif true

  // Sanity check
  bool apx_eql; // [flg] Arguments are indistinguishable
  std::string nfo_msg("Vetting radiative energy conservation");
  for(idx=0;idx<wvl_nbr;idx++){
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       PRC_CMP(1.0), // I [frc] Target argument
       rfl_flx[idx]+trn_flx[idx]+abs_flx[idx], // I [frc] Approximation to target argument
       PRC_CMP(1.0e-4), // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
  } // !idx

  // Delete obsolete arrays

  // Output results of module
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // Wavelength dimension
  const int *dmn_wvl(&wvl_dmn); // Pointer to wavelength dimension
  //  const int *dmn_scl(&wvl_dmn); // Dummy argument, not used

  var_mtd_sct var_mtd[]={
    {0,"abs_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux absorptance","units","fraction"},
    {0,"rfl_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux reflectance","units","fraction"},
    {0,"trn_flx",NC_FLOAT,1,dmn_wvl,"long_name","Flux transmittance","units","fraction"}
  }; // !var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  const int dmn_nbr_max(1); // [nbr] Maximum number of dimensions allowed in single variable in output file
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  // After writing, do NOT delete arrays because deallocate() does that and [rfl,trn,abs]_flx are const
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_flx"),rfl_flx); // delete []rfl_flx;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_flx"),trn_flx); // delete []trn_flx;
  rcd=nco_put_var(nc_out,static_cast<std::string>("abs_flx"),abs_flx); // delete []abs_flx;

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !rt_lop()

int // O [enm] Return success code
rt_rfl_snw // [fnc] Radiative transfer for snow reflectance
(const int &nc_out, // I [fl] netCDF file for output 
 const prc_cmp *flx_frc_drc_sfc, // I [frc] Surface insolation fraction in direct beam
 const prc_cmp &rfl_gnd_dff, // I [frc] Diffuse reflectance of ground (beneath snow)
 const rt_cls &rt_obj, // I [rt] Radiative transfer object
 const prc_cmp &slr_zen_ngl_cos) // I [frc] Cosine solar zenith angle
{
  /* Purpose: Determine snow reflectance using method of Wiscombe and Warren (1980) (WiW80)
     WiW80 use delta-Eddington approximation
     nc_out is modified as variables are defined
     ncks -C -F -q -d wvl,0.5e-6 -v tau_ext,ss_alb,asm_prm,'rfl_spc_snw.?' ${DATA}/mie/mie.nc
     LGGE tests:
     
  */

  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("rt_rfl_snw"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Obtain snow layer optical properties
  const prc_cmp *tau_ext(rt_obj.tau_ext_get()); // [frc] Extinction optical depth
  const prc_cmp *asm_prm(rt_obj.asm_prm_get()); // [frc] Asymmetry parameter
  const prc_cmp *ss_alb(rt_obj.ss_alb_get()); // [frc] Single scattering albedo

  const long wvl_nbr(rt_obj.wvl_nbr_get()); // [nbr] Number of wavelength bins

  // Output arrays
  prc_cmp *rfl_spc_snw_dff=new prc_cmp[wvl_nbr]; // [frc] Snow spectral flux reflectance for diffuse radiation
  prc_cmp *rfl_spc_snw_dff_apx=new prc_cmp[wvl_nbr]; // [frc] Snow spectral flux reflectance for diffuse radiation, semi-inifinite approximation
  prc_cmp *rfl_spc_snw_drc=new prc_cmp[wvl_nbr]; // [frc] Snow spectral flux reflectance for direct radiation
  prc_cmp *rfl_spc_snw_drc_apx=new prc_cmp[wvl_nbr]; // [frc] Snow spectral flux reflectance for direct radiation
  prc_cmp *rfl_spc_snw=new prc_cmp[wvl_nbr]; // [frc] Snow spectral flux reflectance

  // Delta-scaled optical properties
  prc_cmp tau_ext_scl; // [frc] Extinction optical depth, scaled
  prc_cmp asm_prm_scl; // [frc] Asymmetry parameter, scaled
  prc_cmp ss_alb_scl; // [frc] Single scattering albedo, scaled

  // Terms defined in Equation 3
  prc_cmp P; // [frc] Parameter in WiW80 Equation 3
  prc_cmp Q; // [frc] Parameter in WiW80 Equation 3
  prc_cmp Q_mns; // [frc] Parameter in WiW80 Equation 3
  prc_cmp Q_pls; // [frc] Parameter in WiW80 Equation 3
  prc_cmp a_str; // [frc] Parameter in WiW80 Equation 3
  prc_cmp b_str; // [frc] Parameter in WiW80 Equation 3
  prc_cmp gamma; // [frc] Parameter in WiW80 Equation 3
  prc_cmp xi; // [frc] Parameter in WiW80 Equation 3

  gsl_sf_result nsw_dbl; // [frc] GSL result structure
  gsl_sf_result nsw_dbl2; // [frc] GSL result structure
  prc_cmp asm_prm_sqr; // [frc] Asymmetry parameter squared
  prc_cmp eqn3_RHS; // [frc] RHS of Equation 3 
  prc_cmp eqn3_RHS_trm1; // [frc] Term 1 on RHS of Equation 3 
  prc_cmp eqn3_RHS_trm2; // [frc] Term 2 on RHS of Equation 3 
  prc_cmp eqn3_RHS_trm3; // [frc] Term 3 on RHS of Equation 3 
  prc_cmp eqn3_RHS_trm4; // [frc] Term 4 on RHS of Equation 3 
  prc_cmp eqn4_RHS_trm1; // [frc] Term 1 on RHS of Equation 4 
  prc_cmp eqn4_RHS_trm2; // [frc] Term 2 on RHS of Equation 4 
  prc_cmp eqn6_RHS; // [frc] RHS of Equation 6 
  prc_cmp eqn6_RHS_trm1; // [frc] Term 1 on RHS of Equation 6 
  prc_cmp eqn6_RHS_trm2; // [frc] Term 2 on RHS of Equation 6 
  prc_cmp eqn6_RHS_trm3; // [frc] Term 3 on RHS of Equation 6 
  prc_cmp eqn6_RHS_trm4; // [frc] Term 4 on RHS of Equation 6 
  prc_cmp eqn7_RHS_trm1; // [frc] Term 1 on RHS of Equation 7 
  prc_cmp eqn7_RHS_trm2; // [frc] Term 2 on RHS of Equation 7 
  prc_cmp tmp1; // [frc] Temporary variable
  prc_cmp tmp2; // [frc] Temporary variable

  // Perform delta-scaling of layer optical properties
  for(idx=0;idx<wvl_nbr;idx++){
    assert(tau_ext[idx] != 0.0);
    asm_prm_sqr=asm_prm[idx]*asm_prm[idx]; // [frc] Asymmetry parameter squared
    tmp1=ss_alb[idx]*asm_prm_sqr;
    assert(tmp1 != 1.0);
    tau_ext_scl=(1.0-tmp1)*tau_ext[idx]; // [frc] Extinction optical depth, scaled
    assert(asm_prm[idx] != -1.0);
    asm_prm_scl=asm_prm[idx]/(1.0+asm_prm[idx]); // [frc] Asymmetry parameter, scaled
    ss_alb_scl=(1.0-asm_prm_sqr)*ss_alb[idx]/(1.0-asm_prm_sqr*ss_alb[idx]); // [frc] Single scattering albedo, scaled

    a_str=1.0-ss_alb_scl*asm_prm_scl;
    b_str=asm_prm_scl/a_str;
    xi=std::sqrt(3.0*a_str*(1.0-ss_alb_scl));
    P=2.0*xi/(3.0*a_str);
    gamma=(1.0-rfl_gnd_dff)/(1.0+rfl_gnd_dff);

    // Avoid potential numerical problems by first computing semi-infinite solutions
    // Equation 4 is semi-infinite direct albedo
    eqn4_RHS_trm1=ss_alb_scl/(1.0+P); // [frc] Term 1 on RHS of Equation 4, semi-infinite direct albedo
    eqn4_RHS_trm2=(1.0-b_str*xi*slr_zen_ngl_cos)/(1.0+xi*slr_zen_ngl_cos); // [frc] Term 2 on RHS of Equation 4, semi-infinite direct albedo
    rfl_spc_snw_drc_apx[idx]=eqn4_RHS_trm1*eqn4_RHS_trm2; // [frc] Snow spectral flux reflectance for direct radiation, semi-infinite approximation

    // Equation 7 is semi-infinite diffuse albedo
    eqn7_RHS_trm1=2.0*ss_alb_scl/(1.0+P);
    eqn7_RHS_trm2=(1.0+b_str)*(xi-std::log(1.0+xi))/(xi*xi)-b_str/2.0;
    rfl_spc_snw_dff_apx[idx]=eqn7_RHS_trm1*eqn7_RHS_trm2; // [frc] Snow spectral flux reflectance for diffuse radiation, semi-inifinite approximation

    if(tau_ext_scl < 1000.0){
      
      // Equation 3 is analytic delta-Eddington solution to direct reflectance
      Q_pls=(gamma+P)*std::exp(+xi*tau_ext_scl); // fxm: potential overflow?
      Q_mns=(gamma-P)*std::exp(-xi*tau_ext_scl);
      Q=(1.0+P)*Q_pls-(1.0-P)*Q_mns; // fxm: potential underflow?
      
      tmp2=xi*slr_zen_ngl_cos;
      eqn3_RHS_trm1=P*(1.0-gamma+ss_alb_scl*b_str)+ss_alb_scl*(1.0+b_str)*(gamma*xi*slr_zen_ngl_cos-P)/(1.0-tmp2*tmp2);
      eqn3_RHS_trm2=std::exp(-tau_ext_scl/slr_zen_ngl_cos);
      eqn3_RHS_trm3=ss_alb_scl*b_str*(Q_pls-Q_mns);
      eqn3_RHS_trm4=ss_alb_scl*(1.0+b_str)*(Q_pls/(1.0+tmp2)-Q_mns/(1-tmp2));
      eqn3_RHS=2.0*eqn3_RHS_trm1*eqn3_RHS_trm2-eqn3_RHS_trm3+eqn3_RHS_trm4;
      
      rfl_spc_snw_drc[idx]=eqn3_RHS/Q; // [frc] Snow spectral flux reflectance for direct radiation
      
      // Equation 6 is analytic delta-Eddington solution to diffuse reflectance
      eqn6_RHS_trm1=(1.0-gamma+ss_alb_scl*b_str)*(1.0-tau_ext_scl)-gamma*ss_alb_scl*(1.0+b_str)/(1.0-ss_alb_scl);
      eqn6_RHS_trm1*=2.0*P*std::exp(-tau_ext_scl);
      eqn6_RHS_trm2=ss_alb_scl*(1.0+b_str)*(2.0/(xi*xi)+gamma*tau_ext_scl/(1.0-ss_alb_scl))+(1.0-gamma+ss_alb_scl*b_str)*tau_ext_scl*tau_ext_scl;
      
      gsl_error_handler_t *old_gsl_error_handler; // [fnc] GSL error handler
      old_gsl_error_handler=gsl_set_error_handler_off(); // [fnc] GSL error handler
      
      // gsl_sf_expint_Ei_e() computes Exponential integral
      rcd+=gsl_sf_expint_Ei_e(-tau_ext_scl,&nsw_dbl); // fxm: potential underflow?
      eqn6_RHS_trm2*=2.0*P*nsw_dbl.val;
      eqn6_RHS_trm3=2.0*ss_alb_scl*(1.0+b_str)/(xi*xi);
      rcd+=gsl_sf_expint_Ei_e(-(1.0+xi)*tau_ext_scl,&nsw_dbl); // fxm: potential underflow?
      rcd+=gsl_sf_expint_Ei_e(-(1.0-xi)*tau_ext_scl,&nsw_dbl2); // fxm: potential underflow?
      eqn6_RHS_trm3*=Q_pls*(nsw_dbl.val+xi-std::log(1.0+xi))-Q_mns*(nsw_dbl2.val-xi-std::log(1.0-xi)); // fxm: potential underflow?
      eqn6_RHS_trm4=ss_alb_scl*b_str*(Q_pls-Q_mns); // fxm: potential underflow?
      eqn6_RHS=eqn6_RHS_trm1-eqn6_RHS_trm2+eqn6_RHS_trm3-eqn6_RHS_trm4;
      rfl_spc_snw_dff[idx]=eqn6_RHS/Q; // [frc] Snow spectral flux reflectance for diffuse radiation
      
      // Restore original handler
      old_gsl_error_handler=gsl_set_error_handler(NULL); // [fnc] GSL error handler
      
    }else{ // !tau_ext_scl
      if(dbg_lvl_get() > dbg_sbr) std::cout << "INFO: Using semi-infinite snowpack approximation" << std::endl;
      rfl_spc_snw_drc[idx]=rfl_spc_snw_drc_apx[idx]; // [frc] Snow spectral flux reflectance for direct radiation
      rfl_spc_snw_dff[idx]=rfl_spc_snw_dff_apx[idx]; // [frc] Snow spectral flux reflectance for diffuse radiation
    } // !tau_ext_scl
    
    rfl_spc_snw[idx]=flx_frc_drc_sfc[idx]*rfl_spc_snw_drc[idx]+(1.0-flx_frc_drc_sfc[idx])*rfl_spc_snw_dff[idx]; // [frc] Snow spectral flux reflectance

  } // !idx
  
    // Sanity check
    // Flux reflectance, transmittance and absorptance should sum to unity
  if(true){
    std::cout << "Snowpack optics:" << std::endl;
    std::cout << "  Snow particles: Optical depth = " << tau_ext[0] << ", Single scattering albedo = " << ss_alb[0] << ", Asymmetry parameter = " << asm_prm[0] << std::endl;
    std::cout << "  Solar zenith angle: Cosine = " << slr_zen_ngl_cos << ", Degrees = " << std::acos(slr_zen_ngl_cos)*180.0/mth::cst_M_PIl << std::endl;
    std::cout << "  Diffuse reflectance of underlying ground = " << rfl_gnd_dff << ", Fraction of incident flux in direct beam = " << flx_frc_drc_sfc[0] << std::endl;
    std::cout << "  Reflectance: Direct = " << rfl_spc_snw_drc[0] << " (direct semi-infinite = " << rfl_spc_snw_drc_apx[0] << "), Diffuse = " << rfl_spc_snw_dff[0] << " (diffuse semi-infinite = " << rfl_spc_snw_dff_apx[0] << "), Total = " << rfl_spc_snw[0] << std::endl;
  } // !true
  
  // Delete obsolete arrays
  
  // Output module results
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // Wavelength dimension
  const int *dmn_wvl(&wvl_dmn); // Pointer to wavelength dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars CLIP

  /* Check that current dimension size in output file equals array length
     (writing record variables fails silently until record dimension is written) */
  const long wvl_nbr_crr(nco_inq_dimlen(nc_out,wvl_dmn)); // [nbr] Current number of wavelength bins in output file
  assert(wvl_nbr == wvl_nbr_crr);

  var_mtd_sct var_mtd[]={
    {0,"flx_frc_drc_sfc",NC_FLOAT,1,dmn_wvl,"long_name","Surface insolation fraction in direct beam","units","fraction"},
    {0,"rfl_spc_snw",NC_FLOAT,1,dmn_wvl,"long_name","Snow spectral flux reflectance","units","fraction"},
    {0,"rfl_spc_snw_dff",NC_FLOAT,1,dmn_wvl,"long_name","Snow spectral flux reflectance for diffuse radiation","units","fraction"},
    {0,"rfl_spc_snw_dff_apx",NC_FLOAT,1,dmn_wvl,"long_name","Snow spectral flux reflectance for diffuse radiation, semi-inifinite approximation","units","fraction"},
    {0,"rfl_spc_snw_drc",NC_FLOAT,1,dmn_wvl,"long_name","Snow spectral flux reflectance for direct radiation","units","fraction"},
    {0,"rfl_spc_snw_drc_apx",NC_FLOAT,1,dmn_wvl,"long_name","Snow spectral flux reflectance for direct radiation, semi-infinite approximation","units","fraction"},
    {0,"rfl_gnd_dff",NC_FLOAT,0,dmn_scl,"long_name","Diffuse reflectance of ground (beneath snow)","units","fraction"}
  }; // !var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  const int dmn_nbr_max(1); // [nbr] Maximum number of dimensions allowed in single variable in output file
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_drc_sfc"),flx_frc_drc_sfc); delete []flx_frc_drc_sfc;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_spc_snw_drc"),rfl_spc_snw_drc); delete []rfl_spc_snw_drc;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_spc_snw_drc_apx"),rfl_spc_snw_drc_apx); delete []rfl_spc_snw_drc_apx;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_spc_snw_dff"),rfl_spc_snw_dff); delete []rfl_spc_snw_dff;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_spc_snw_dff_apx"),rfl_spc_snw_dff_apx); delete []rfl_spc_snw_dff_apx;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_spc_snw"),rfl_spc_snw); delete []rfl_spc_snw;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_gnd_dff"),rfl_gnd_dff);

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // !rt_rfl_snw()

// Global functions with C++ linkages end
