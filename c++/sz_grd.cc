// $Id$ 

// Implementation (declaration) of size grid class

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <sz_grd.hh> // Size grids

// sz_grd_cls class

// Initialize static members 
sng2sng_map sz_grd_cls::opt2abb_map=sz_grd_cls::opt2abb_map_mk();
sng2sng_map sz_grd_cls::opt2abb_map_mk(){ // Create abbreviation map
  sng2sng_map map_tmp;
  map_tmp.insert(sng2sng_map::value_type("log","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("lgr","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("logarithmic","logarithmic"));
  map_tmp.insert(sng2sng_map::value_type("lin","linear"));
  map_tmp.insert(sng2sng_map::value_type("lnr","linear"));
  map_tmp.insert(sng2sng_map::value_type("linear","linear"));
  map_tmp.insert(sng2sng_map::value_type("regular","linear"));
  map_tmp.insert(sng2sng_map::value_type("rgl","linear"));
  map_tmp.insert(sng2sng_map::value_type("RGL","linear"));
  map_tmp.insert(sng2sng_map::value_type("REG","linear"));
  // NB: Return a value to initialize a static class member
  return map_tmp;
} // end sz_grd_cls::opt2abb_map_mk()
std::string sz_grd_cls::opt2abb(const std::string opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  itr=opt2abb_map.find(opt_sng);
  if(itr == opt2abb_map.end()) err_prn("sz_grd_cls::opt2abb",opt_sng+" is unknown");
  return itr->second;
} // end sz_grd_cls::opt2abb()

// Public member functions
sz_grd_cls::sz_grd_cls // Default constructor
(const prc_cmp sz_mnm_arg, // [m] Minimum size in distribution
 const prc_cmp sz_mxm_arg, // [m] Maximum size in distribution
 const long sz_nbr_arg, // [nbr] Number of size bins
 const std::string grd_sng_arg) // [sng] Type of size grid
{
  rcd=0; // [enm] Return success code
  sz_mnm=sz_mnm_arg; // [m] Minimum size in distribution
  sz_mxm=sz_mxm_arg; // [m] Maximum size in distribution
  assert(sz_mnm != sz_mxm);
  sz_nbr=sz_nbr_arg; // [nbr] Number of size bins 
  grd_sng=grd_sng_arg; // [sng] Type of size grid
  rcd+=allocate(); // [fnc] Allocate dynamic memory for object
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end sz_grd_cls constructor
sz_grd_cls::~sz_grd_cls(){ // Destructor
  rcd+=deallocate(); // [fnc] Free dynamic memory of object
} // end sz_grd_cls destructor
sz_grd_cls::sz_grd_cls(const sz_grd_cls &orig):sz_nbr(orig.sz_nbr) // Copy constructor
{
  sz_mnm=orig.sz_mnm;
  sz_mxm=orig.sz_mxm;
  rcd+=allocate(); // [fnc] Allocate dynamic memory for object
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end sz_grd_cls copy constructor
const sz_grd_cls &sz_grd_cls::operator=(const sz_grd_cls &RHS) // Assignment operator
{
  if(&RHS != this){
    if(sz_nbr != RHS.sz_nbr){
      sz_nbr=RHS.sz_nbr;
      rcd+=deallocate(); // [fnc] Free dynamic memory of object
      rcd+=allocate(); // [fnc] Allocate dynamic memory for object
    } // endif
    sz_mnm=RHS.sz_mnm;
    sz_mxm=RHS.sz_mxm;
    rcd+=recompute(); // [fnc] Recompute properties of object
  } // endif
  return *this;
} // end sz_grd_cls assignment operator

std::string sz_grd_cls::typ_get()const{return grd_sng;} // [fnc] Get grid string
void sz_grd_cls::typ_set(const std::string sng){ // [fnc] Set grid type
  grd_sng=sng;
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end sz_grd_cls::typ_set()

const prc_cmp *sz_grd_cls::sz_ctr_get()const{return sz_ctr;} // [fnc] Get pointer to const data
const prc_cmp *sz_grd_cls::sz_dlt_get()const{return sz_dlt;} // [fnc] Get pointer to const data
const prc_cmp *sz_grd_cls::sz_grd_get()const{return sz_grd;} // [fnc] Get pointer to const data
const prc_cmp *sz_grd_cls::sz_max_get()const{return sz_max;} // [fnc] Get pointer to const data
const prc_cmp *sz_grd_cls::sz_min_get()const{return sz_min;} // [fnc] Get pointer to const data
long sz_grd_cls::nbr_get()const{return sz_nbr;} // [nbr] Number of elements of size grid
void sz_grd_cls::nbr_set // [fnc] Set number of elements of size grid
(const long sz_nbr_arg) // [nbr] Number of elements of size grid
{
  // Purpose: Set number of elements of size grid  
  if(sz_nbr_arg != 1){;} // CEWU Compiler Error Warning Usage
  err_prn("sz_grd_cls::nbr not fully implemented yet");
} // end sz_grd_cls::nbr_set()

// Private member functions
int // [enm] Return success code
sz_grd_cls::allocate(){ // [fnc] Allocate dynamic memory
  sz_ctr=new prc_cmp[sz_nbr];
  sz_dlt=new prc_cmp[sz_nbr];
  sz_grd=new prc_cmp[sz_nbr+1];
  sz_max=new prc_cmp[sz_nbr];
  sz_min=new prc_cmp[sz_nbr];
  return rcd; // [enm] Return success code
} // end sz_grd_cls::allocate()

int // [enm] Return success code
sz_grd_cls::deallocate(){ // [fnc] Free dynamic memory of object
  delete []sz_ctr;
  delete []sz_dlt;
  delete []sz_grd;
  delete []sz_max;
  delete []sz_min;
  return rcd; // [enm] Return success code
} // end sz_grd_cls::deallocate()

int // [enm] Return success code
sz_grd_cls::recompute(){ // [fnc] Recompute properties of object
  int idx; // [idx] counting index
  const std::string sbr_nm("sz_grd_cls::recompute"); // [sng] Name of subroutine
  if(grd_sng == "logarithmic"){ // Space size bins logarithmically
    prc_cmp max_min_ratio=sz_mxm/sz_mnm;
    prc_cmp series_ratio=std::pow(max_min_ratio,PRC_CMP(1.0)/sz_nbr);
    if(sz_mnm==0.0) err_prn(sbr_nm,"sz_mnm = 0.0 for size grid type = "+grd_sng);
    sz_min[0]=sz_mnm;
    for(idx=1;idx<sz_nbr;idx++){ // NB: Loop starts from 1
      sz_min[idx]=sz_min[idx-1]*series_ratio;
      sz_max[idx-1]=sz_min[idx];
    } // end loop over idx
  }else if(grd_sng == "linear"){ // Space size bins linearly
    prc_cmp sz_ncr=(sz_mxm-sz_mnm)/sz_nbr;
    for(idx=0;idx<sz_nbr;idx++){
      sz_min[idx]=sz_mnm+idx*sz_ncr;
      sz_max[idx]=sz_mnm+(idx+1)*sz_ncr;
    } // end loop over idx
  }else{
    err_prn(sbr_nm,grd_sng+" is not a valid grid type");
  } // end else
  sz_max[sz_nbr-1]=sz_mxm; // Ensure final max is not affected by roundoff
  for(idx=0;idx<sz_nbr;idx++) sz_dlt[idx]=sz_max[idx]-sz_min[idx];
  if(grd_sng == "logarithmic"){ // Space size bins logarithmically
    prc_cmp log_min; 
    prc_cmp log_max;
    prc_cmp log_ctr;
    for(idx=0;idx<sz_nbr;idx++){
      // sz_ctr=logarithmic mean for logarithmic grids
      log_min=std::log(sz_min[idx]);
      log_max=std::log(sz_max[idx]);
      log_ctr=0.5*(log_max+log_min);
      sz_ctr[idx]=std::exp(log_ctr);
    } // end loop over idx
  }else if(grd_sng == "linear"){ // Space size bins linearly
    for(idx=0;idx<sz_nbr;idx++){
      // sz_ctr=linear mean for linear grids
      sz_ctr[idx]=0.5*(sz_max[idx]+sz_min[idx]);
    } // end loop over idx
  } // end else
  for(idx=0;idx<sz_nbr;idx++) sz_grd[idx]=sz_min[idx];
  sz_grd[sz_nbr]=sz_max[sz_nbr-1];
  return rcd; // [enm] Return success code
} // end sz_grd_cls::recompute()

