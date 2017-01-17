// $Id$ 

// Purpose: Utilities for libcsz_c++  

/* Copyright (C) 1997--2017 Charlie Zender
   This software is distributed under the terms of the General Public License
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Source: Routines are subset of CSZ's ~/c++/utl.[cc/hh]

// Usage:
// #include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

#ifndef UTL_HH // Contents have not yet been inserted in current source file  
#define UTL_HH

// C++ headers
#include <complex> // Standard C++ complex class
#include <iostream> // Standard C++ I/O streams cout, cin
#include <sstream> // Standard C++ string stream processing
#include <string> // Standard C++ string class

// Standard C headers
#include <cassert> // Assertions
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc. (needed for f77_prm_dat_prn())
#include <cstdlib> // abort, exit, getopt, malloc, strtod, strtol
#include <ctime> // Machine time

// 3rd party vendors

// Personal headers
#include <dbg.hh> // Debugging constants

// Namespaces

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif // EXIT_SUCCESS is not defined in SUN4

// Prototype functions that have C linkages
extern "C" {
  char *prg_nm_get(void);
  unsigned short dbg_lvl_get(void);

  void f77_blk_dat_prn
  (const int dat_nbr_per_ln, // [nbr] Number of datum per line
   const char *dat_sng, // [sng] Fortran data dimension statement
   const long dat_nbr, // [nbr] Number of data
   const prc_cmp *dat); // [frc] Data
  // end f77_blk_dat_prn() prototype

  void 
  f90_prm_dat_prn
  (const bool ftn_fxd, // [flg] Fixed format
   const int dat_nbr_per_ln, // [nbr] Number of datum per line
   const char *dat_sng, // [sng] Fortran data dimension statement
   const long dat_nbr, // [nbr] Number of data
   const prc_cmp *dat); // [frc] Data
  // end f90_prm_dat_prn() prototype

} // end extern C

// Define inline'd functions in header so source is visible to calling files
inline prc_cmp sign_cpv(const prc_cmp a,const prc_cmp b){return (b > 0 ? a : -a);}
inline prc_cmp min_cpv(const prc_cmp a,const prc_cmp b){return (a < b ? a : b);}
inline prc_cmp max_cpv(const prc_cmp a,const prc_cmp b){return (a > b ? a : b);}
inline long min_lng(const long a,const long b){return (a < b ? a : b);}
inline long max_lng(const long a,const long b){return (a > b ? a : b);}

// Prototype functions with C++ linkages
void vec_set // [fnc] Broadcast scalar value to vector
(prc_cmp *vec, // O [frc] Vector to fill with value
 const long vec_nbr, // I [nbr] Vector size
 const prc_cmp val); // I [frc] Value to fill vector with
// end vec_set() prototype

int // O [enm] Return success code
vec_set // [fnc] Fill complex vector with real, imaginary components
(const long vec_nbr, // I [nbr] Vector size
 const prc_cmp *val_rl, // I [frc] Real array
 const prc_cmp *val_img, // I [frc] Imaginary array
 std::complex<prc_cmp> *vec); // O [frc] Complex array
// end vec_set() prototype

bool // [flg] Array is monotonic
mnt_chk // [fnc] Determine monotonicity of array
(const prc_cmp * const grd, // [frc] Array
 const long &grd_nbr); // [nbr] Number of elements
// end mnt_chk() prototype

long // O [nbr] Number of consecutive monotonic elements
mnt_nbr // [fnc] Bracket monotonic portion of array
(const prc_cmp * const grd, // I [frc] Array
 const long &grd_nbr, // I [nbr] Number of elements
 const long &idx_srt=0); // I [idx] Index of starting element (optional)
// end mnt_nbr() prototype

bool mnt_dcr(const prc_cmp * const grd,const long &grd_nbr);
bool mnt_ncr(const prc_cmp * const grd,const long &grd_nbr);
prc_cmp vec_max(const prc_cmp * const crd,const long &crd_nbr);
prc_cmp vec_min(const prc_cmp * const crd,const long &crd_nbr);
long vec_val2idx(const prc_cmp * const grd,const long &grd_nbr,const prc_cmp &val);
void whenfvlt(const long crd_nbr,const prc_cmp *crd,const prc_cmp *val,long *vld_idx,long &vld_nbr);
void whenfgt(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr);
void whenflt(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr);
void whenfeq(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr);

std::string CVS_Date_2_sng(std::string CVS_Date);
void dbg_lvl_tst(void);
void Exit_gracefully(void);
void dbg_prn(std::string msg);
void dbg_prn(std::string sbr_nm,std::string msg);
void dbg_prn(std::string prg_nm,std::string sbr_nm,std::string msg);
void err_prn(std::string msg);
void err_prn(std::string sbr_nm,std::string msg);
void err_prn(std::string prg_nm,std::string sbr_nm,std::string msg);
void wrn_prn(std::string sbr_nm,std::string msg);
void wrn_prn(std::string prg_nm,std::string sbr_nm,std::string msg);

std::string // O [sng] Parsed command line
cmd_ln_sng // [fnc] Re-construct command line from arguments
(const int argc, // I [nbr] Argument count
 const char * const * const argv); // I [sng] Command line argument values
// end cmd_ln_sng() prototype

std::string // O [sng] Full file name = drc/fl_nm
drc_pfx // [fnc] Intelligently return drc/fl_nm
(const std::string drc, // I [nbr] Directory to prepend
 const std::string fl_nm); // I [sng] File name
// end drc_pfx() prototype

// Templates
template<class prc_T>bool // O [flg] Arguments are indistinguishable
apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
(const std::string &sbr_nm_cll, // I [sng] Subroutine name of calling routine
 const bool &flg_vrb, // I [flg] Verbose output
 const prc_T &arg_trg, // I [frc] Target argument
 const prc_T &arg_apx, // I [frc] Approximation to target argument
 const prc_T &rlt_prc, // I [frc] Relative precision
 const std::string &msg); // I [sng] Descriptive message of context
// end apx_eql_chk() prototype
template<class prc_T>bool // O [flg] Arguments are indistinguishable
apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
(const std::string &sbr_nm_cll, // I [sng] Subroutine name of calling routine
 const bool &flg_vrb, // I [flg] Verbose output
 const prc_T &arg_trg, // I [frc] Target argument
 const prc_T &arg_apx, // I [frc] Approximation to target argument
 const prc_T &rlt_prc, // I [frc] Relative precision
 const std::string &msg) // I [sng] Descriptive message of context
{
  /* Purpose: Uniform procedure to determine whether arguments are indistinguishable
     Relative precision to which "approximation" must match "target argument" is user-specified
     Arguments are indistinguishable iff fabs((arg_trg-arg_apx)/arg_trg) <= rlt_prc */
  const std::string sbr_nm("apx_eql_chk"); // [sng] Subroutine name
  bool apx_eql_chk_flg(false); // [flg] Alice and Bob are indistinguishable

  // Handle degenerate case where all values are zero
  if(arg_trg == 0.0 && arg_apx == 0.0) return true;

  if(arg_trg != 0.0)
    if(std::fabs((arg_trg-arg_apx)/arg_trg) <= rlt_prc) 
      return true;

  if(arg_trg == 0.0)
    if(std::fabs(arg_trg-arg_apx) <= rlt_prc) 
      return true;

  if(flg_vrb && !apx_eql_chk_flg){
    if(arg_trg == 0.0){
      std::cout << prg_nm_get() << ": "+sbr_nm+"() called from "+sbr_nm_cll+"(): "+msg+" std::fabs((" << arg_trg << " - " << arg_apx << ")" << " = " << std::fabs(arg_trg-arg_apx) << " > " << rlt_prc << std::endl;
    }else{
      std::cout << prg_nm_get() << ": "+sbr_nm+"() called from "+sbr_nm_cll+"(): "+msg+" std::fabs((" << arg_trg << " - " << arg_apx << ")/" << arg_trg << " = " << std::fabs((arg_trg-arg_apx)/arg_trg) << " > " << rlt_prc << std::endl;
    } // endif
  } // endif

  return false; // [flg] Arguments are indistinguishable
} // end apx_eql_chk()

template<class val_T>std::string nbr2sng(const val_T nbr); // O [sng] Number stored as string
template<class val_T> // [obj] Object type
std::string // O [sng] Number stored as string
  nbr2sng // [fnc] Convert number to string
  (const val_T val) // I [frc] Number to convert to string
{
  /* Purpose: Convert number to string
     Method taken from CUED C++ FAQ
     NB: String streams require GCC >= 3.x */
  const std::string sbr_nm("nbr2sng"); // [sng] Subroutine name
  std::ostringstream sng_srm_out; // [srm] Output string stream
  if(sng_srm_out << val) return sng_srm_out.str(); 
  // If control reaches this point then error was encountered
  err_prn(sbr_nm,"Unable to convert number to string");
  // Need return value here to avoid compiler warnings
  return sng_srm_out.str(); // 
} // end nbr2sng()

template<class val_T>std::string nbr2sng(const val_T nbr,int dcm_plc_prc); // O [sng] Number stored as string
template<class val_T> // [obj] Object type
std::string // [sng] Number stored as string
  nbr2sng // [fnc] Convert number to string
  (const val_T val, // I [frc] Number to convert to string
   const int dcm_plc_prc) // I [nbr] Decimal places of precision
{
  /* Purpose: Convert number to string
     Method taken from CUED C++ FAQ
     NB: String streams require GCC 3.x */
  const std::string sbr_nm("nbr2sng"); // [sng] Subroutine name
  std::ostringstream sng_srm_out; // [srm] Output string stream
  sng_srm_out.precision(dcm_plc_prc);
  if(sng_srm_out << val) return sng_srm_out.str(); 
  // If control reaches this point then error was encountered
  err_prn(sbr_nm,"Unable to convert number to string");
  // Need return value here to avoid compiler warnings
  return sng_srm_out.str(); // 
} // end nbr2sng()

template<class val_T>void sng2nbr(const std::string sng,val_T *nbr); // O [sng] String stored as number
template<class val_T> // [obj] Object type
void
sng2nbr // [fnc] Convert string to number
(const std::string sng, // I [frc] String to convert to number
 val_T *nbr) // O [frc] Number that was converted from string
{
  /* Purpose: Convert string to number
     Method inverse of nbr2sng() (see below)
     This method suggested by Martin York in
     comp.lang.c++.moderated on 20050522 in response to my thread on 
     "Cross-platform strtoll() functionality" */
  const std::string sbr_nm("sng2nbr"); // [sng] Subroutine name
  std::stringstream sng_srm_in(sng); // [srm] Input string stream
  if(!(sng_srm_in >> *nbr)) err_prn(sbr_nm,"Unable to convert string \""+sng+"\" to number");
} // end sng2nbr()

template<class val_T>val_T sng2nbr(const std::string sng,const val_T nbr); // O [sng] String stored as number
template<class val_T> // [obj] Object type
val_T // O [sng] String stored as number
sng2nbr // [fnc] Convert string to number
(const std::string sng, // I [frc] String to convert to number
 const val_T nbr) // I [frc] Number of type to convert to (not touched)
{
  /* Purpose: Convert string to number
     Method inverse of nbr2sng() (see below)
     This method suggested by Martin York in
     comp.lang.c++.moderated on 20050522 in response to my thread on 
     "Cross-platform strtoll() functionality" */
  val_T val_out; // O [nbr] String stored as number
  const std::string sbr_nm("sng2nbr"); // [sng] Subroutine name
  std::stringstream sng_srm_in(sng); // [srm] Input string stream
  if(sng_srm_in >> val_out) return val_out; else err_prn(sbr_nm,"Unable to convert string \""+sng+"\" to number");
  // Redundant return value here avoids two compiler warnings
  return nbr; // O [nbr] String stored as number
} // end sng2nbr()

template<class prc_T>void
rvr_vec // [ptr] Reverse an array in place
(prc_T *dat_out, // [frc] Array to be reversed
 const long &out_nbr) // [nbr] Number of element in array
{
  // Purpose: Reverse an array in place
  // Usage: rvr_vec(crd,nbr)
  std::string sbr_nm("rvr_vec"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  long idx; // [idx] Counting index
  prc_T *dat_swp=new prc_T[out_nbr]; // Swap array

  // Copy input array into local swap array
  for(idx=0;idx<out_nbr;idx++){
    dat_swp[idx]=dat_out[idx]; // Swap array
  } // end loop over idx
  // Reverse local swap array into output array
  for(idx=0;idx<out_nbr;idx++){
    dat_out[idx]=dat_swp[out_nbr-idx-1]; // Array to be reversed
  } // end loop over idx

  delete []dat_swp; // Swap array
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
} // end rvr_vec()

template<class prc_T>prc_T * // O [ptr] Copy of vec_in
vec_cpy // [fnc] Copy input vector
(const prc_T * const vec_in, // I [frc] Vector
 const long &in_nbr); // I [frc] Number of elements
// end vec_cpy() prototype
template<class prc_T>prc_T * // O [ptr] Copy of vec_in
vec_cpy // [fnc] Copy input vector
(const prc_T * const vec_in, // I [frc] Vector
 const long &in_nbr) // I [frc] Number of elements
{
  /* Purpose: Return copy of input array
     Warning: Calling routine is responsible for freeing this memory
     Usage: crd_cpy=vec_cpy(crd,nbr) */
  std::string sbr_nm("vec_cpy"); // [sng] Subroutine name
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  long idx; // [idx] Counting index
  prc_T *vec_out=new prc_T[in_nbr]; // Swap array

  // Copy input array into duplicate array
  for(idx=0;idx<in_nbr;idx++){
    vec_out[idx]=vec_in[idx]; // Copy of input array
  } // end loop over idx

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return vec_out; // Copy of input array
} // end vec_cpy()

#endif // UTL_HH  
