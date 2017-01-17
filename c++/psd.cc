// $Id$ 

// Implementation (declaration) of particle size distribution class

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <psd.hh> // Particle size distribution

// psd_cls class

// Static members 
int psd_cls::nst_nbr=0; // [nbr] Number of instantiated class members

// Initialize static member functions 
int psd_cls::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members
sng2sng_map psd_cls::opt2abb_map=psd_cls::opt2abb_map_mk();
sng2sng_map psd_cls::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;
  map_tmp.insert(sng2sng_map::value_type("log","lognormal"));
  map_tmp.insert(sng2sng_map::value_type("lgn","lognormal"));
  map_tmp.insert(sng2sng_map::value_type("lognormal","lognormal"));
  map_tmp.insert(sng2sng_map::value_type("gam","gamma"));
  map_tmp.insert(sng2sng_map::value_type("gmm","gamma"));
  map_tmp.insert(sng2sng_map::value_type("gamma","gamma"));
  map_tmp.insert(sng2sng_map::value_type("Gamma","gamma"));
  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end psd_cls::opt2abb_map_mk()
std::string psd_cls::opt2abb(const std::string &opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("psd_cls::opt2abb",opt_sng+" is unknown");
  return itr->second;
} // end psd_cls::opt2abb()

// Public member functions
psd_cls::psd_cls // [fnc] Default constructor
(const prc_cmp &rds_nma_arg, // [m] Number median radius analytic
 const prc_cmp &gsd_anl_arg, // [frc] Geometric standard deviation
 const prc_cmp &cnc_nbr_anl_arg, // [# m-3] Number concentration analytic
 const prc_cmp &dns_prt_arg, // [kg m-3] Particle density
 const std::string &dst_sng_arg, // [sng] Distribution string
 sz_grd_cls *sz_grd_arg) // [sct] Size grid
{
  // Purpose: Default constructor for size distribution objects

  std::string sbr_nm("psd_cls::psd_cls"); // [sng] Name of subroutine

  // Error status of this object should remain zero
  rcd=0; // [enm] Return success code
  
  // Set rcm_flg to false until essential members have been set
  // set() functions will still perform range-checking and other diagnostics
  // but overhead of numeric computations will be avoided
  rcm_flg=false; // [flg] Invoke recompute() on set() calls
  mss_frc_anl_set_flg=false; // [flg] Mass fraction analytic has been set

  // fxm: valgrind flags this unusual allocation
  sz_grd_cls *sz_grd_dfl=new sz_grd_cls[1]; // Default size grid

  // Set size grid first so correct size is known for allocate()
  if(sz_grd_arg == 0){
    rcd+=sz_grd_set(sz_grd_dfl); // [sct] Size grid
  }else{
    delete []sz_grd_dfl; // Default size grid
    rcd+=sz_grd_set(sz_grd_arg); // [sct] Size grid
  } // endif

  rcd+=allocate(); // [fnc] Allocate dynamic memory for object

  // Copy arguments
  dst_sng=dst_sng_arg; // [sng] Distribution string

  // Use set() functions so range-checking and other diagnostics are performed
  rcd+=rds_nma_set(rds_nma_arg); // [m] Number median radius analytic
  rcd+=gsd_anl_set(gsd_anl_arg); // [frc] Geometric standard deviation
  rcd+=dns_prt_set(dns_prt_arg); // [kg m-3] Particle density
  rcd+=cnc_nbr_anl_set(cnc_nbr_anl_arg); // [# m-3] Number concentration analytic

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  if(rcd != 0) err_prn(sbr_nm,"Error constructing object");    

  nst_nbr++; // [nbr] Number of instantiated class members
} // end psd_cls::psd_cls()

psd_cls::~psd_cls(){ // Destructor
  std::string sbr_nm("psd_cls::~psd_cls"); // [sng] Name of subroutine
  rcd+=deallocate(); // [fnc] Free dynamic memory for object
  if(rcd != 0) err_prn(sbr_nm,"Error destroying object");
  nst_nbr--; // [nbr] Number of instantiated class members
} // end psd_cls::~psd_cls()

std::string psd_cls::typ_get()const{return dst_sng;} // [fnc] Get distribution string
void psd_cls::typ_set(const std::string &sng){ // [fnc] Set distribution type
  dst_sng=sng;
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end psd_cls::typ_set()

prc_cmp psd_cls::rds_nma_get()const{return rds_nma;} // [m] Number median radius analytic
prc_cmp psd_cls::dmt_nma_get()const{return dmt_nma;} // [m] Number median diameter analytic
prc_cmp psd_cls::gsd_anl_get()const{return gsd_anl;} // [frc] Geometric standard deviation
prc_cmp psd_cls::mss_anl_get()const{return mss_anl;} // [kg m-3] Mass concentration analytic
prc_cmp psd_cls::dns_prt_get()const{return dns_prt;} // [kg m-3] Particle density
prc_cmp psd_cls::cnc_nbr_anl_get()const{return cnc_nbr_anl;} // [# m-3] Number concentration analytic

prc_cmp // [kg m-3] Mass fraction analytic
psd_cls::mss_frc_anl_get()const{
  if(!mss_frc_anl_set_flg) err_prn("psd_cls::mss_frc_anl_get","Attempt to get undefined mss_frc_anl");
  return mss_frc_anl;
} // end psd_cls::mss_frc_anl_get()

int // [enm] Return success code
psd_cls::mss_frc_anl_set // [fnc] Set mass fraction analytic for each mode of ensemble
(psd_cls *psd_lst, // [sct] List of particle size distributions
 const long &psd_nbr) // [nbr] Number of particle modes
{
  // Purpose: Compute and store mass fraction for each mode in a list of modes

  long psd_idx; // [idx] Counting index for aerosol mode
  prc_cmp mss_anl_ttl(0.0); // [kg m-3] Mass concentration analytic total

  // Sanity check
  if(psd_nbr <= 0) err_prn("psd_cls::mss_frc_anl_set","psd_nbr <= 0");
  // Sum up mass concentration in all modes
  for(psd_idx=0;psd_idx<psd_nbr;psd_idx++) mss_anl_ttl+=psd_lst[psd_idx].mss_anl_get(); // [kg m-3] Mass concentration analytic total

  // Store mass fraction in each mode
  for(psd_idx=0;psd_idx<psd_nbr;psd_idx++) psd_lst[psd_idx].mss_frc_anl_set(psd_lst[psd_idx].mss_anl_get()/mss_anl_ttl); // [frc] Mass fraction analytic

  mss_frc_anl_set_flg=true; // [flg] Mass fraction analytic has been set

  return rcd; // [enm] Return success code
} // end psd_cls::mss_frc_anl_set()

int // [enm] Return success code
psd_cls::mss_frc_anl_set(const prc_cmp &mss_frc_anl_arg) // [frc] Mass fraction analytic
{
  // Purpose: Set size distribution mass fraction analytic
  mss_frc_anl=mss_frc_anl_arg; // [frc] Mass fraction analytic
  if(mss_frc_anl < 0.0 || mss_frc_anl > 1.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::mss_frc_anl_set() reports out of range error with mss_frc_anl = " << mss_frc_anl << std::endl;
  mss_frc_anl_set_flg=true; // [flg] Mass fraction analytic has been set
  return rcd; // [enm] Return success code
} // end psd_cls::mss_frc_anl_set()

int // [enm] Return success code
psd_cls::psd_prn // [fnc] Print formatted list of properties of particle size distribution(s)
(const psd_cls *psd_lst, // [sct] List of particle size distributions
 const long &psd_nbr) // [nbr] Number of particle modes
{
  // Purpose: Print formatted list of properties of particle size distribution(s)

  long psd_idx; // [idx] Counting index for particle mode
  std::string sbr_nm("psd_prn"); // [sng] Name of subroutine

  // Sanity check
  if(psd_nbr <= 0) err_prn("psd_cls::psd_prn","psd_nbr <= 0");

  if(psd_lst[0].typ_get() == "lognormal"){
    (void)std::fprintf(stdout,"  %4s %9s %9s %9s %9s\n","Mode","dmt_nma","gsd_anl","mss_frc","cnc_nbr");
    (void)std::fprintf(stdout,"  %4s %9s %9s %9s %9s\n","    ","  um   ","  frc  ","  frc  "," # m-3 ");
    for(psd_idx=0;psd_idx<psd_nbr;psd_idx++){
      (void)std::fprintf(stdout,"  %4ld %9.5f %9.5f %9.5f %9.3e\n",psd_idx,psd_lst[psd_idx].dmt_nma_get()*1.0e6,psd_lst[psd_idx].gsd_anl_get(),psd_lst[psd_idx].mss_frc_anl_get(),psd_lst[psd_idx].cnc_nbr_anl_get());
    } /* end loop over psd_idx */
  }else if(psd_lst[0].typ_get() == "gamma"){
    // fxm: Gamma distributions needs to be re-validated, they are broken now
    err_prn(sbr_nm,"Unable to print gamma distributions");
  }else{ 
    err_prn(sbr_nm,"Unknown psd_lst[0].typ_get()");
  } // end else

  return rcd; // [enm] Return success code
} // end psd_cls::psd_prn()

int // [enm] Return success code
psd_cls::prs_psd_sng  // [fnc] Set size distribution parameters from string
(const std::string &psd_sng) // [sng] Triplet specifying size distribution
{
  // Purpose: Set size distribution parameters from a command line string argument

  prc_cmp cnc_nbr_anl_lcl; // [# m-3] Number concentration analytic
  prc_cmp rds_nma_lcl; // [m] Number median radius analytic
  prc_cmp rds_nma_mcr; // [um] Number median radius analytic
  prc_cmp gsd_anl_lcl; // [frc] Geometric standard deviation
  prc_cmp mss_anl_lcl; // [kg m-3] Mass concentration analytic

  std::string sbr_nm("psd_cls::prs_psd_sng"); // [sng] Name of subroutine
  std::string dlm_sng(","); // [sng] Delimiter string
  unsigned long dlm_1_lcn; // [idx] Location of first comma
  unsigned long dlm_2_lcn; // [idx] Location of second comma
  
  // Find positions of commas and number of characters between (non-inclusive) them
  dlm_1_lcn=psd_sng.find(dlm_sng); // [idx] Location of first comma
  if(dlm_1_lcn == std::string::npos) err_prn(sbr_nm,"Invalid psd specification, dlm_1_lcn");
  
  dlm_2_lcn=psd_sng.find(dlm_sng,dlm_1_lcn+1); // [idx] Location of second comma
  if(dlm_2_lcn == std::string::npos) err_prn(sbr_nm,"Invalid psd specification, dlm_2_lcn");
  
  // Exit if any argument is length zero
  if(dlm_1_lcn == 0 || dlm_2_lcn == dlm_1_lcn+1 || dlm_2_lcn == psd_sng.length()-1) err_prn(sbr_nm,"Invalid psd specification, zero length");
    
  bool mss_anl_flg(false); // [enm] User-specified Mass concentration
  if(psd_sng.find("mass") != std::string::npos || (psd_sng.find("Mass") != std::string::npos) || (psd_sng.find("mss") != std::string::npos)){
    mss_anl_flg=true; // [enm] User-specified Mass concentration
  } // endif
  
  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls

  rds_nma_mcr=static_cast<prc_cmp>(std::strtod(psd_sng.substr(0,dlm_1_lcn).c_str(),(char **)NULL)); // [um] Number median radius analytic
  rds_nma_lcl=rds_nma_mcr*1.0e-6; // [m] Number median radius analytic
  rcd+=rds_nma_set(rds_nma_lcl); // [m] Number median radius analytic
  
  gsd_anl_lcl=static_cast<prc_cmp>(std::strtod(psd_sng.substr(dlm_1_lcn+1,(dlm_2_lcn-dlm_1_lcn)-1).c_str(),(char **)NULL)); // [frc] Geometric standard deviation
  rcd+=gsd_anl_set(gsd_anl_lcl); // [frc] Geometric standard deviation
  
  // User is expected to enter mass fractions as simple mass concentrations
  if(mss_anl_flg){ // [enm] User-specified mass concentration
    mss_anl_lcl=static_cast<prc_cmp>(std::strtod(psd_sng.substr(dlm_2_lcn+1).c_str(),(char **)NULL)); // [kg m-3] Mass concentration analytic
    rcd+=mss_anl_set(mss_anl_lcl); // [kg m-3] Mass concentration analytic
  }else{
    cnc_nbr_anl_lcl=static_cast<prc_cmp>(std::strtod(psd_sng.substr(dlm_2_lcn+1).c_str(),(char **)NULL)); // [# m-3] Total concentration
    rcd+=cnc_nbr_anl_set(cnc_nbr_anl_lcl); // [# m-3] Total concentration
  } // endif

  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls

  rcd+=recompute(); // [fnc] Recompute properties of object

  if(rcd != 0) err_prn(sbr_nm,"Error parsing size distribution structures");    

  return rcd; // [enm] Return success code
} // end psd_cls::prs_psd_sng()

int // [enm] Return success code
psd_cls::dmt_nma_set(const prc_cmp &dmt_nma_arg) // [m] Number median diameter analytic
{
  // Purpose: Set size distribution number median diameter
  dmt_nma=dmt_nma_arg; // [m] Number median diameter analytic
  if(dmt_nma < 0.0 || dmt_nma > 1.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::dmt_nma_set() reports out of range error with dmt_nma = " << dmt_nma << std::endl;
  rds_nma=dmt_nma/2.0; // [m] Number median radius analytic
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::dmt_nma_set()

int // [enm] Return success code
psd_cls::rds_nma_set(const prc_cmp &rds_nma_arg) // [m] Number median radius analytic
{
  // Purpose: Set size distribution number median radius
  rds_nma=rds_nma_arg; // [m] Number median radius analytic
  if(rds_nma < 0.0 || rds_nma > 1.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::rds_nma_set() reports out of range error with rds_nma = " << rds_nma << std::endl;
  dmt_nma=rds_nma*2.0; // [m] Number median diameter analytic
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::rds_nma_set()

int // [enm] Return success code
psd_cls::gsd_anl_set(const prc_cmp &gsd_anl_arg) // [frc] Geometric standard deviation
{
  // Purpose: Set size distribution geometric standard deviation
  gsd_anl=gsd_anl_arg; // [frc] Geometric standard deviation
  if(gsd_anl < 0.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::gsd_anl_set() reports out of range error with gsd_anl = " << gsd_anl << std::endl;
  if(gsd_anl > 3.0) std::cerr << prg_nm_get() << ": WARNING psd_cls::gsd_anl_set() reports unusually large gsd_anl = " << gsd_anl << std::endl;
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::gsd_anl_set()

int // [enm] Return success code
psd_cls::mss_anl_set(const prc_cmp &mss_anl_arg) // [kg m-3] Mass concentration analytic
{
  /* Purpose: Set size distribution mass concentration
     Mass concentrations may be set by using the "mass" keyword to an aerosol mode
     Mass fractions may be implemented very simply using this function by 
     setting the analytic mass concentration numerically equal to the 
     desired mass fraction.
     Thus if the user desires a mass fraction of 0.5, then setting the aerosol
     analytic mass concentration to 0.5 kg m-3 will, in conjunction with 
     setting the other modes appropriately, work. */
  prc_cmp cnc_nbr_anl_lcl; // [# m-3] Number concentration analytic
  // Set mss_anl now, then set cnc_nbr_anl to be consistent
  mss_anl=mss_anl_arg; // [kg m-3] Mass concentration analytic
  if(mss_anl < 0.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::mss_anl_set() reports out of range error with mss_anl = " << mss_anl << std::endl;
  if(mss_anl > 1.0) std::cerr << prg_nm_get() << ": WARNING psd_cls::mss_anl_set() reports unusually large mss_anl = " << mss_anl << std::endl;
  // Convert specified mass concentration to number concentration
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  cnc_nbr_anl_lcl=mss_anl/(dns_prt*(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_nma,PRC_CMP(3.0))*std::exp(4.5*std::log(gsd_anl)*std::log(gsd_anl))); // [# m-3] Number concentration analytic
  rcd+=cnc_nbr_anl_set(cnc_nbr_anl_lcl); // [# m-3] Number concentration analytic
  return rcd; // [enm] Return success code
} // end psd_cls::mss_anl_set()

int // [enm] Return success code
psd_cls::dns_prt_set(const prc_cmp &dns_prt_arg) // [kg m-3] Particle density
{
  // Purpose: Set size distribution particle density
  dns_prt=dns_prt_arg; // [kg m-3] Particle density
  if(dns_prt < 0.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::dns_prt_set() reports out of range error with dns_prt = " << dns_prt << std::endl;
  if(dns_prt > 1.0e6) std::cerr << prg_nm_get() << ": WARNING psd_cls::dns_prt_set() reports unusually large dns_prt = " << dns_prt << std::endl;
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::dns_prt_set()

int // [enm] Return success code
psd_cls::cnc_nbr_anl_set(const prc_cmp &cnc_nbr_anl_arg) // [# m-3] Number concentration analytic
{
  // Purpose: Set size distribution total number concentration
  cnc_nbr_anl=cnc_nbr_anl_arg; // [# m-3] Number concentration analytic
  if(cnc_nbr_anl < 0.0) std::cerr << prg_nm_get() << ": ERROR psd_cls::cnc_nbr_anl_set() reports out of range error with cnc_nbr_anl = " << cnc_nbr_anl << std::endl;
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::cnc_nbr_anl_set()

int // [enm] Return success code
psd_cls::sz_grd_set(sz_grd_cls *sz_grd_arg) // [sct] Size grid
{
  // Purpose: Set size grid for size distribution
  sz_grd=sz_grd_arg; // [sct] Size grid
  sz_nbr=sz_grd->nbr_get(); // [nbr] Number of elements of size grid
  if(sz_nbr <= 0) std::cerr << prg_nm_get() << ": ERROR psd_cls::sz_grd_set() reports out of range error with sz_nbr = " << sz_nbr << std::endl;
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end psd_cls::sz_grd_set()

// Private member functions
int // [enm] Return success code
psd_cls::allocate(){ // [fnc] Allocate dynamic memory for object
  // Purpose: Allocate dynamic memory for object
  dst=new prc_cmp[sz_nbr]; // [# m-3 m-1] Number distribution
  cnc=new prc_cmp[sz_nbr]; // [# m-3] Number concentration
  xsa=new prc_cmp[sz_nbr]; // [m2] Cross-sectional area
  sfc=new prc_cmp[sz_nbr]; // [m2] Surface area
  vlm=new prc_cmp[sz_nbr]; // [m3] Volume
  mss=new prc_cmp[sz_nbr]; // [kg] Mass
  return rcd; // [enm] Return success code
} // end psd_cls::allocate()

int // [enm] Return success code
psd_cls::deallocate(){ // [fnc] Free dynamic memory for object
  // Purpose: Free dynamic memory for object
  delete []dst; // [# m-3 m-1] Number distribution
  delete []cnc; // [# m-3] Number concentration
  delete []xsa; // [m2] Cross-sectional area
  delete []sfc; // [m2] Surface area
  delete []vlm; // [m3] Volume
  delete []mss; // [kg] Mass
  return rcd; // [enm] Return success code
} // end psd_cls::deallocate()

int // [enm] Return success code
psd_cls::recompute(){ // [fnc] Recompute properties of object
  // Purpose: Update size distribution due to changes in any member
  // Ensure all dependencies are correctly handled
  const prc_cmp *sz_ctr=sz_grd->sz_ctr_get();
  const prc_cmp *sz_dlt=sz_grd->sz_dlt_get();
  prc_cmp cnc_nbr_rsl(0.0); // [# m-3] Number concentration resolved
  prc_cmp mss_rsl(0.0); // [kg m-3] Mass concentration resolved
  prc_cmp sfc_rsl(0.0); // [m2 m-3] Surface area concentration resolved
  prc_cmp vlm_rsl(0.0); // [m3 m-3] Volume concentration resolved
  prc_cmp xsa_rsl(0.0); // [m2 m-3] Cross-sectional area concentration resolved
  long idx;
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  for(idx=0;idx<sz_nbr;idx++){
    /* fxm: c++ TODO #107: Must set dst based on real PDF before this works
       Right now, the continuous arrays in the psd class are bogus
       Only the analytic numbers work */
    dst[idx]=0.0; // [# m-3] Number concentration
    cnc[idx]=dst[idx]*sz_dlt[idx]; // [# m-3] Number concentration
    cnc_nbr_rsl=cnc_nbr_rsl+cnc[idx]; // [# m-3] Number concentration resolved

    xsa[idx]=mth::cst_M_PIl*sz_ctr[idx]*sz_ctr[idx]; // [m2] Cross-sectional area of sphere of given size
    xsa_rsl+=xsa[idx]*cnc[idx]; // [m2 m-3] Cross-sectional area concentration resolved

    sfc[idx]=4.0*xsa[idx]; // [m2] Surface area of sphere of given size
    sfc_rsl+=sfc[idx]*cnc[idx]; // [m2 m-3] Surface area concentration resolved

    vlm[idx]=(4.0/3.0)*mth::cst_M_PIl*std::pow(sz_ctr[idx],PRC_CMP(3.0)); // [m3] Volume of sphere of given size
    vlm_rsl+=vlm[idx]*cnc[idx]; // [m3 m-3] Volume concentration resolved

    mss[idx]=vlm[idx]*dns_prt; // [kg] Mass of sphere of given size
    mss_rsl+=mss[idx]*cnc[idx]; // [kg m-3] Mass concentration resolved
  } // end loop over sz
  // fxm  mss_anl=vlm_anl*dns_prt; // [kg m-3] Mass concentration analytic

  mss_anl=dns_prt*(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_nma,PRC_CMP(3.0))*cnc_nbr_anl*std::exp(4.5*std::log(gsd_anl)*std::log(gsd_anl)); // [kg m-3] Mass concentration analytic

  return rcd; // [enm] Return success code
} // end psd_cls::recompute()

