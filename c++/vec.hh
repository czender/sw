// $Id$ 

// Purpose: Vector utilities for C++ programs

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <vec.hh> // Vector functions ntp_vec(), rbn_vec()

#ifndef VEC_HH // Contents have not yet been inserted in current source file  
#define VEC_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers

// 3rd party vendors

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <xtr.hh> // Extrapolation class

// Prototype functions with C++ linkages

void 
vec_tst // [fnc] Test vector utilities
(Xtr_cls xtr_LHS, // I [enm] LHS extrapolation flags
 Xtr_cls xtr_RHS); // I [enm] RHS extrapolation flags

int // [enm] Return success code
rbn_vec // [fnc] Rebin a vector from one grid onto another grid
(const long in_nbr, // (I) [nbr] Input grid size
 prc_cmp *grd_in, // (I) Input grid
 prc_cmp *dat_in, // (I) Input data
 const long out_nbr, // (I) [nbr] Output grid size
 prc_cmp *grd_out, // (I) Output grid
 prc_cmp *dat_out, // (O) Output (re-binned) data
 const Xtr_cls &xtr_LHS, // (I) [flg] LHS extrapolation flags
 const Xtr_cls &xtr_RHS); // (I) [flg] RHS extrapolation flags
// end rbn_vec() prototype

void 
trp_area_vec // [fnc] Compute area of array of trapezoids
(const long in_nbr, // (I) Input coordinate size
 const prc_cmp *grd_in, // (I) Input grid
 const prc_cmp *dat_in, // (I) Input data
 prc_cmp *trp_area); // (O) Area of each trapezoid defined by input data
// end trp_area_vec() prototype

// Templates
template<class prc_T>int // O [enm] Return success code
ntp_vec // [fnc] Project a vector onto another vector
(const long in_nbr, // I [nbr] Input coordinate size
 const prc_cmp *crd_in, // I [crd] Input coordinate
 const prc_T *dat_in, // I [frc] Input data (on coordinate)
 const long out_nbr, // I [nbr] Output coordinate size
 const prc_cmp *crd_out, // I [crd] Output coordinate
 prc_T *dat_out, // O [frc] Output (interpolated) data (on coordinate)
 const Xtr_cls &xtr_LHS, // I [enm] LHS extrapolation flags
 const Xtr_cls &xtr_RHS); // I [enm] RHS extrapolation flags
// end ntp_vec() prototype
template<class prc_T>int // O [enm] Return success code
ntp_vec // [fnc] Project a vector onto another vector
(const long in_nbr, // I [nbr] Input coordinate size
 const prc_cmp *crd_in, // I [crd] Input coordinate
 const prc_T *dat_in, // I [frc] Input data (on coordinate)
 const long out_nbr, // I [nbr] Output coordinate size
 const prc_cmp *crd_out, // I [crd] Output coordinate
 prc_T *dat_out, // O [frc] Output (interpolated) data (on coordinate)
 const Xtr_cls &xtr_LHS, // I [enm] LHS extrapolation flags
 const Xtr_cls &xtr_RHS) // I [enm] RHS extrapolation flags
{
  /* Purpose: Project a vector onto another vector
     Input and output coordinate vectors should be monotonic, either increasing or decreasing
     Internally, this routine only works on monotonically increasing arrays
     If input array(s) is/are decreasing, routine will automatically reverse array(s) before arithmetic, then reverse again before passing answers back */

  /* Nomenclature:
     crd has a slightly different meaning in ntp_vec() than in rbn_vec()
     Following is the mathematical representation for ntp_vec():
     crd_in defines a coordinate grid associated with values dat_in
     crd_out is a new (user) coordinate grid onto which dat_in is linearly interpolated as dat_out

     ^
     |     
     |                     dat_in[2]-> x
     |     
     |                            x <-dat_out[1]
     |     
     |           dat_in[1]-> x 
     |                                         dat_in[N]-> x 
     |                  x <-dat_out[0]
     |     
     | dat_in[0]-> x                                            x <-dat_out[N]
     |                            
     --------------^----^----^----^----^---------|---------^----^--->
     |             |    *    |    *    |         |         |    *
     |             |    *    |    *    |         |         |    *
     Input:    crd_in[0]*crd_in[1]*crd_in[2]  ......   crd_in[N]* 
     |                  *         *                             *
     Output:        crd_out[0]crd_out[1]      ......        crd_out[K] 

     Part 1:
     First, crd_in and crd_out are checked for required monotonicity
     Copies of crd_in and crd_out are made and reversed if necessary so both are monotonically increasing
     dat_in undergoes the same procedure as crd_in to be consistent
     Arrays based on monotonically increasing coordinates are given an _mnt suffix: crd_in_mnt, crd_out_mnt, dat_in_mnt

     Part 2:
     Once monotonically increasing coordinate arrays are ready, main loop is entered
     Interpolation algorithm is a one-pass, non-vectorized loop through out_idx=[0..out_nbr-1]
     Loop begins with identifying brk_lft_idx and brk_rgt_idx for crd_out[out_idx]

     brk_lft_idx is the C-index of the point in crd_in[] immediately to the left of crd_out[out_idx]
     Legal values of brk_lft_idx are, therefore, [0...in_nbr-1]
     brk_rgt_idx is the C-index of the point in crd_in[] immediately to the right of crd_out[out_idx]
     Legal values of brk_rgt_idx are, therefore, [0...in_nbr-1]
     "Legal values" of brk_lft_idx are values which may be accessed in crd_in_mnt,dat_in_mnt
     Not all possible values of brk_lft_idx are "legal"

     brk_lft_idx == -1 when crd_out_mnt[out_idx] < crd_in_mnt[0]
     This case requires LHS extrapolation and crd_in[brk_lft_idx] should never be accessed
     In most cases, 0 <= brk_lft_idx <= in_nbr-2 
     In these cases, brk_rgt_idx=brk_lft_idx+1 and interpolation is straightforward
     brk_lft_idx == in_nbr-1 when crd_out_mnt[out_idx] > crd_in_mnt[in_nbr-1] so brk_rgt_idx DNE
     This case requires RHS extrapolation and crd_in[brk_rgt_idx] should never be accessed

     Loop ends when crd_out[out_idx] has been set
     For efficiency, next iteration of loop starts with brk_lft_idx from previous iteration

     Part 3:
     crd_in_mnt, crd_out_mnt, dat_in_mnt are deleted since they are (possibly reversed) copies of input
     dat_out is reversed if crd_out is monotonically decreasing
  */

  // Output
  int rcd(0); // [enm] Return success code
  // Local Workspace
  long brk_rgt_idx; // [idx] Counting index for RHS bracket
  long out_idx; // [idx] Counting index for output grid

  // Main Code
  std::string sbr_nm("ntp_vec");
  unsigned short dbg_lvl=dbg_lvl_get();
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Vet input
  if(in_nbr <= 1) err_prn(sbr_nm,"in_nbr <= 1");
  if(out_nbr < 1) err_prn(sbr_nm,"out_nbr < 1");

  // Find direction of monotonicity
  bool in_ncr=mnt_ncr(crd_in,in_nbr);
  bool in_dcr=!in_ncr;
  bool out_ncr=mnt_ncr(crd_out,out_nbr);
  bool out_dcr=!out_ncr;
      
  // Convert input and output arrays to monotonic increasing arrays
  // fxm: Is it possible to perform this copy operation only when necessary rather than always?
  prc_cmp *crd_in_mnt=vec_cpy(crd_in,in_nbr); // [frc] Input coordinate monotonically increasing
  prc_cmp *crd_out_mnt=vec_cpy(crd_out,out_nbr); // [frc] Output coordinate monotonically increasing
  prc_T *dat_in_mnt=vec_cpy(dat_in,in_nbr); // [frc] Input data on monotonically increasing coordinate
  if(in_dcr){
    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Reversing input grid");
    rvr_vec(crd_in_mnt,in_nbr);
    rvr_vec(dat_in_mnt,in_nbr);
  } // endif in_dcr
  if(out_dcr){
    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Reversing output grid");
    rvr_vec(crd_out_mnt,out_nbr);
  } // endif out_dcr

  // Initialize bracketing index
  long brk_lft_idx(0);
  // Loop over desired output coordinates
  for(out_idx=0;out_idx<out_nbr;out_idx++){
    // Order of conditions is important since second condition is illegal if brk_lft_idx >= in_nbr
    while((brk_lft_idx < in_nbr) && (crd_in_mnt[brk_lft_idx] < crd_out_mnt[out_idx])){
      brk_lft_idx++;
    } // end while
    brk_lft_idx--;
    // Handle identity interpolation separately to preserve symmetry in extrapolation code 
    if(brk_lft_idx != in_nbr-1){
      if(crd_in_mnt[brk_lft_idx+1] == crd_out_mnt[out_idx]){
	dat_out[out_idx]=dat_in_mnt[brk_lft_idx+1];
	if(brk_lft_idx == -1) brk_lft_idx=0; // Reset brk_lft_idx to 0 so next while loop works
	continue; // Jump to next iteration
      } // endif
    } // endif
    if(brk_lft_idx == -1){
      // LHS Extrapolation required
      // Degenerate case: crd_out_mnt[out_idx] < crd_in_mnt[0]
      brk_lft_idx=0; // Reset brk_lft_idx to 0 so next while loop works
      if(xtr_LHS.xtr_vrb) std::cerr << "WARNING: ntp_vec() output value dat_out[" << out_idx << "] at coordinate crd_out_mnt[" << out_idx << "] = " << crd_out_mnt[out_idx] << " requires LHS extrapolation beyond leftmost valid coordinate at crd_in_mnt[" << brk_lft_idx << "] = " << crd_in_mnt[brk_lft_idx] << ". Nearest valid datum is dat_in_mnt[" << brk_lft_idx << "] = " << dat_in_mnt[brk_lft_idx] << std::endl;
      // Extrapolation options are presented in decreasing order of preference
      if(!xtr_LHS.xtr_fll){
	err_prn(sbr_nm,"Full LHS extrapolation required but not permitted");
      }else if(xtr_LHS.xtr_fll_nil){
	dat_out[out_idx]=0.0;
      }else if(xtr_LHS.xtr_fll_ngh){
	dat_out[out_idx]=dat_in_mnt[0];
      }else if(xtr_LHS.xtr_fll_lnr){
	dat_out[out_idx]=dat_in_mnt[0]-
	  (crd_in_mnt[0]-crd_out_mnt[out_idx])*
	  (dat_in_mnt[1]-dat_in_mnt[0])/(crd_in_mnt[1]-crd_in_mnt[0]);
      }else{
	err_prn(sbr_nm,"Unknown xtr_LHS");
      } // endif xtr_LHS
      if(xtr_LHS.xtr_vrb) std::cerr << xtr_LHS.xtr_sng << " yields dat_out[" << out_idx << "] = " << dat_out[out_idx] << std::endl;
    }else if(brk_lft_idx < in_nbr-1){
      // Normal case: crd_out_mnt is interpolable
      brk_rgt_idx=brk_lft_idx+1; 
      // NB: brk_rgt_idx is ALWAYS greater than brk_lft_idx
      // This simulaneously meets two criteria:
      // 1. Divide by zero errors are impossible in the next step
      // 2. The identity interpolation is satisfied since crd_dlt == 0.0: 
      // i.e., If crd_out_mnt[idx] == crd_in_mnt[brk_lft_idx] then dat_out[out_idx] := dat_in_mnt[brk_lft_idx]
      // Linearly interpolate
      dat_out[out_idx]=
	dat_in_mnt[brk_lft_idx]+
	(crd_out_mnt[out_idx]-crd_in_mnt[brk_lft_idx])*
	(dat_in_mnt[brk_rgt_idx]-dat_in_mnt[brk_lft_idx])/
	(crd_in_mnt[brk_rgt_idx]-crd_in_mnt[brk_lft_idx]);
    }else if(brk_lft_idx == in_nbr-1){
      // RHS Extrapolation required
      // Degenerate case: brk_lft_idx is last element of crd_in_mnt 
      brk_rgt_idx=brk_lft_idx;
      if(xtr_RHS.xtr_vrb) std::cerr << "WARNING: ntp_vec() output value dat_out[" << out_idx << "] at coordinate crd_out_mnt[" << out_idx << "] = " << crd_out_mnt[out_idx] << " requires RHS extrapolation beyond rightmost valid coordinate at crd_in_mnt[" << brk_rgt_idx << "] = " << crd_in_mnt[brk_rgt_idx] << ". Nearest valid datum is dat_in_mnt[" << brk_rgt_idx << "] = " << dat_in_mnt[brk_rgt_idx] << std::endl;
      // Extrapolation options are presented in decreasing order of preference
      if(!xtr_RHS.xtr_fll){
	err_prn(sbr_nm,"Full RHS extrapolation required but not permitted");
      }else if(xtr_RHS.xtr_fll_nil){
	dat_out[out_idx]=0.0;
      }else if(xtr_RHS.xtr_fll_ngh){
	dat_out[out_idx]=dat_in_mnt[in_nbr-1];
      }else if(xtr_RHS.xtr_fll_lnr){
	dat_out[out_idx]=dat_in_mnt[in_nbr-1]+
	  (crd_out_mnt[out_idx]-crd_in_mnt[in_nbr-1])*
	  (dat_in_mnt[in_nbr-1]-dat_in_mnt[in_nbr-2])/
	  (crd_in_mnt[in_nbr-1]-crd_in_mnt[in_nbr-2]);
      }else{
	err_prn(sbr_nm,"Unknown xtr_RHS");
      } // endif xtr_RHS
      if(xtr_RHS.xtr_vrb) std::cerr << xtr_RHS.xtr_sng << " yields dat_out[" << out_idx << "] = " << dat_out[out_idx] << std::endl; 
    }else{
      err_prn(sbr_nm,"Unforeseen value of brk_lft_idx");
    } // end else RHS
  } // end loop over output coordinates
      
  // Convert input and output arrays to original direction of monotonicity
  delete []crd_in_mnt;
  delete []dat_in_mnt;
  delete []crd_out_mnt;
  //  if(in_dcr){
  //    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Un-reversing input grid");
    //    rvr_vec(crd_in,in_nbr);
    //    rvr_vec(dat_in,in_nbr);
  //  } // endif in_dcr
  if(out_dcr){
    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Un-reversing output grid");
    // rvr_vec(crd_out,out_nbr);
    rvr_vec(dat_out,out_nbr);
  } // endif in_dcr
  
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end ntp_vec()

template<class prc_T>prc_T // O [frc] Arrival point
ntp_vec_one // [fnc] Perform one point interpolation call to ntp_vec()
(const long crd_nbr, // I [nbr] Number of elements in input arrays
 const prc_cmp *crd_in, // I [crd] Abscissa values of input
 const prc_T *dat_in, // I [frc] Ordinate values of input
 prc_cmp crd_out); // O [frc] Abscissa value of output
// end ntp_vec_one() prototype
template<class prc_T>prc_T // O [frc] Arrival point
ntp_vec_one // [fnc] Perform one point interpolation call to ntp_vec()
(const long crd_nbr, // I [nbr] Number of elements in input arrays
 const prc_cmp *crd_in, // I [crd] Abscissa values of input
 const prc_T *dat_in, // I [frc] Ordinate values of input
 prc_cmp crd_out) // O [frc] Abscissa value of output
{
  /* Purpose: Perform one point interpolation call to ntp_vec()
     This function is fairly useful, e.g., at inverting functions
     The function takes care of the messy overhead needed to call ntp_vec()
     fxm: crd_in and dat_in are actually const but may be doubly-reversed in ntp_vec(), so find a way to change these back to const */
  int rcd(0); // [enm] Return success code
  const long out_nbr(1);
  prc_cmp *crd_out_ptr=&crd_out;
  prc_T dat_out;
  Xtr_cls xtr_fll_lnr("xtr_fll_lnr");
  rcd+=ntp_vec(crd_nbr,crd_in,dat_in,out_nbr,crd_out_ptr,&dat_out,xtr_fll_lnr,xtr_fll_lnr);
  return dat_out;
} // end ntp_vec_one()

template<class prc_T>int // O [enm] Return success code
ntp_vec_wrp // [fnc] Wrapper for vector interpolation
(const long in_nbr, // I [nbr] Input coordinate size
 const prc_cmp *crd_in, // I [crd] Input coordinate
 const prc_T *dat_in, // I [frc] Input data (on coordinate)
 const long &out_nbr, // I [nbr] Output coordinate size
 const prc_cmp *crd_out, // I [crd] Output coordinate
 prc_T *dat_out, // O [frc] Output (interpolated) data (on coordinate)
 const Xtr_cls &xtr_LHS, // I [enm] LHS extrapolation flags
 const Xtr_cls &xtr_RHS); // I [enm] RHS extrapolation flags
// end ntp_vec_wrp() prototype
template<class prc_T>int // O [enm] Return success code
ntp_vec_wrp // [fnc] Wrapper for vector interpolation
(const long in_nbr, // I [nbr] Input coordinate size
 const prc_cmp *crd_in, // I [crd] Input coordinate
 const prc_T *dat_in, // I [frc] Input data (on coordinate)
 const long &out_nbr, // I [nbr] Output coordinate size
 const prc_cmp *crd_out, // I [crd] Output coordinate
 prc_T *dat_out, // O [frc] Output (interpolated) data (on coordinate)
 const Xtr_cls &xtr_LHS, // I [enm] LHS extrapolation flags
 const Xtr_cls &xtr_RHS) // I [enm] RHS extrapolation flags
{
  /* Purpose: Wrapper routine for interpolating vector
     Routine decides which specific interpolation best fits input data:
     1. Call ntp_vec() once if output coordinate is monotonic
     2. Call ntp_vec() as many times as necessary for non-monotonic output coordinates */

  // Output
  int rcd(0); // [enm] Return success code
  // Main Code
  std::string sbr_nm("ntp_vec_wrp"); // [sng] Subroutine name
  unsigned short dbg_lvl(dbg_lvl_get()); // [enm] Debugging level
  long nmn_cnt(0); // [cnt] Number of monotonic chunks used in interpolation
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  if(mnt_chk(crd_out,out_nbr)){
     
    rcd+=ntp_vec // [fnc] Project a vector onto another vector
      (in_nbr, // I [nbr] Input coordinate size
       crd_in, // I [crd] Input coordinate
       dat_in, // I [frc] Input data (on coordinate)
       out_nbr, // I [nbr] Output coordinate size
       crd_out, // I [crd] Output coordinate
       dat_out, // O [frc] Output (interpolated) data (on coordinate)
       xtr_LHS, // I [enm] LHS extrapolation flags
       xtr_RHS); // I [enm] RHS extrapolation flags
    // end ntp_vec() arguments
    
  }else{ // not monotonic

    // Loop through input arrays finding boundaries of monotonic chunks
    // Send each monotonic chunk separately to ntp_vec()
    long out_idx_srt(0); // [idx] Starting index of current output block
    long out_nbr_rmn(out_nbr); // [nbr] Number of remaining elements to interpolate
    while(out_nbr_rmn > 0){

      const long out_nbr_tmp= // O [nbr] Number of consecutive monotonic elements of crd_out array
	mnt_nbr // [fnc] Bracket monotonic portion of array
	(crd_out, // I [frc] Array
	 out_nbr, // I [nbr] Number of elements
	 out_idx_srt); // I [idx] Index of starting element
      
      rcd+=ntp_vec // [fnc] Project a vector onto another vector
	(in_nbr, // I [nbr] Input coordinate size
	 crd_in, // I [crd] Input coordinate
	 dat_in, // I [frc] Input data (on coordinate)
	 out_nbr_tmp, // I [nbr] Output coordinate size
	 const_cast<const prc_cmp *>(crd_out+out_idx_srt), // I [crd] Output coordinate
	 dat_out+out_idx_srt, // O [frc] Output (interpolated) data (on coordinate)
	 xtr_LHS, // I [enm] LHS extrapolation flags
	 xtr_RHS); // I [enm] RHS extrapolation flags
      // end ntp_vec() arguments
      
      // Increment starting position for next iteration
      out_idx_srt+=out_nbr_tmp; // [idx] Index of starting element
      out_nbr_rmn-=out_nbr_tmp; // [nbr] Number of remaining elements to interpolate
      // Increment counter
      nmn_cnt++; // [cnt] Number of monotonic chunks used in interpolation

      // Sanity check
      if(nmn_cnt > 100) std::cout << prg_nm_get() << ": WARNING "+sbr_nm+"() reports calling ntp_vec() " << nmn_cnt << " times" << std::endl;

    } // end loop over out

    if(dbg_lvl == dbg_crr) std::cout << prg_nm_get() << ": INFO "+sbr_nm+"() reports calling ntp_vec() " << nmn_cnt << " times" << std::endl;

  } // endelse not monotonic

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // end ntp_vec_wrp()

#endif // VEC_HH  
