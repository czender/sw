// $Id$ 

// Purpose: Standalone utilities for C++ programs

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// fxm: many functions defined here, e.g., rvr_vec, mnt*, should be rewritten as generic templates

#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Global functions with C linkages begin

void 
f77_blk_dat_prn
(const int dat_nbr_per_ln, // [nbr] Number of datum per line
 const char *dat_sng, // [sng] Fortran data dimension statement
 const long dat_nbr, // [nbr] Number of data
 const prc_cmp *dat) // [frc] Data
{
  // Purpose: Print array data to stderr in fixed format for Fortran77 block data

  long idx; // [idx] Counting index

  (void)std::fprintf(stderr,"%s",dat_sng);
  for(idx=0;idx<dat_nbr;idx++){
    if(!(idx % dat_nbr_per_ln)) (void)std::fprintf(stderr,"\n     $     ");
    (void)std::fprintf(stderr,"%.7e",dat[idx]);
    if(idx != dat_nbr-1) (void)std::fprintf(stderr,", ");
  } // end loop over wvl
  (void)std::fprintf(stderr," /\n\n");
  
} // end f77_blk_dat_prn()

void 
f90_prm_dat_prn
(const bool ftn_fxd, // [flg] Fixed format
 const int dat_nbr_per_ln, // [nbr] Number of datum per line
 const char *dat_sng, // [sng] Fortran data dimension statement
 const long dat_nbr, // [nbr] Number of data
 const prc_cmp *dat) // [frc] Data
{
  // Purpose: Print array data to stderr in Fortran90 parameter format

  long idx; // [idx] Counting index

  (void)std::fprintf(stderr,"%s",dat_sng);
  for(idx=0;idx<dat_nbr;idx++){
    if(!(idx % dat_nbr_per_ln)){
      if(ftn_fxd) (void)std::fprintf(stderr,"\n     $     "); else (void)std::fprintf(stderr,"&\n           ");
    } // endif
    (void)std::fprintf(stderr,"%.7e",dat[idx]);
    if(idx != dat_nbr-1) (void)std::fprintf(stderr,", ");
  } // end loop over wvl
  (void)std::fprintf(stderr," /)\n\n");
  
} // end f90_prm_dat_prn()

// Declare global functions with C++ linkages

void 
dbg_lvl_tst(void) // [fnc] Return value of dbg_lvl
{
  // Purpose: Return the value of dbg_lvl seen within libcsz_c++
  std::cerr << "libcsz_c++ reports dbg_lvl_get() returns " << dbg_lvl_get() << std::endl;
} // end dbg_lvl_tst()

std::string 
CVS_Date_2_sng(std::string CVS_Date)
{
  // Purpose: Given CVS Date string, return the date string component of it
  std::string sng;

  // If dollarsign Date is present, then keyword or value or both are present
  if(CVS_Date.find("$Date") == std::string::npos){;}
  return sng;
} // end CVS_Date_2_sng()

void 
err_prn(std::string prg_nm,std::string sbr_nm,std::string msg) // [fnc] Print uniform error message and exit 
{
  // Purpose: Print a uniform error message and exit 
  std::cerr << prg_nm << ": ERROR " << sbr_nm << "(): "+msg << std::endl;
#ifdef ABORT_ON_ERROR
  // abort() produces a core dump and traceback information useful to debuggers
  std::abort(); // [fnc] Exit with core dump
#else // !ABORT_ON_ERROR
  // exit() produces no core dump or useful debugger information
  std::exit(EXIT_FAILURE); // [fnc] Exit nicely
#endif // !ABORT_ON_ERROR
} // end err_prn()

void 
err_prn(std::string sbr_nm,std::string msg) // [fnc] Print uniform error message and exit 
{
  // Purpose: Print a uniform error message and exit 
  std::cerr << prg_nm_get() << ": ERROR " << sbr_nm << "(): "+msg << std::endl;
#ifdef ABORT_ON_ERROR
  std::abort(); // [fnc] Produce core dump
#else
  std::exit(EXIT_FAILURE); // [fnc] Exit nicely
#endif // !ABORT_ON_ERROR
} // end err_prn()

void 
err_prn(std::string msg) // [fnc] Print uniform error message and exit
{
  // Purpose: Print uniform error message and exit 
  std::cerr << msg << std::endl;
#ifdef ABORT_ON_ERROR
  std::abort(); // [fnc] Produce core dump
#else
  std::exit(EXIT_FAILURE); // [fnc] Exit nicely
#endif // !ABORT_ON_ERROR
} // end err_prn()

void 
wrn_prn(std::string prg_nm,std::string sbr_nm,std::string msg) // [fnc] Print uniform warning message and exit
{
  // Purpose: Print a uniform warning message
  std::cerr << prg_nm << ": WARNING " << sbr_nm << "(): "+msg << std::endl;
} // end wrn_prn()

void 
wrn_prn(std::string sbr_nm,std::string msg) // [fnc] Print uniform warning message and exit
{
  // Purpose: Print a uniform warning message
  std::cerr << prg_nm_get() << ": WARNING " << sbr_nm << "(): "+msg << std::endl;
} // end wrn_prn()

void 
dbg_prn(std::string prg_nm,std::string sbr_nm,std::string msg) // [fnc] Print uniform debugging message
{
  // Purpose: Print a uniform debugging message
  std::cerr << prg_nm << ": DEBUG " << sbr_nm << "(): "+msg << std::endl;
} // end dbg_prn()

void 
dbg_prn(std::string sbr_nm,std::string msg) // [fnc] Print uniform debugging message
{
  // Purpose: Print a uniform debugging message
  std::cerr << prg_nm_get() << ": DEBUG " << sbr_nm << "(): "+msg << std::endl;
} // end dbg_prn()

void 
dbg_prn(std::string msg) // [fnc] Print uniform debugging message
{
  // Purpose: Print a uniform debugging message
  std::cerr << "DEBUG: "+msg << std::endl;
} // end dbg_prn()

void 
Exit_gracefully(void) // [fnc] Stop clock and exit program gracefully
{
  // Purpose: Stop clock and exit program gracefully
  const std::time_t time_crr_time_t=std::time((std::time_t *)NULL); // [tm] Current date and time
  const std::string time_bfr_end=std::ctime(&time_crr_time_t); // [sng] Current date and time
  std::cerr << "\tEnd = " << time_bfr_end;

  std::exit(EXIT_SUCCESS);
} // end Exit_gracefully() 

bool
mnt_ncr(const prc_cmp * const grd,const long &grd_nbr)
{
  /* Purpose: Determine whether input array is monotonically increasing
     Returns true when array is monotonically increasing, false otherwise
     Routine prints an error when input array is not monotonic */
  std::string sbr_nm("mnt_ncr"); // [sng] Subroutine name
  if(grd_nbr <= 1) return true;
  const bool mnt(mnt_chk(grd,grd_nbr));
  bool ncr(false);
  if(mnt){
    if(grd[0] < grd[1]) ncr=true; 
  }else{
    // fxm: 20030718 use wrn_prn() rather than err_prn() so CAM_SW/CAM_LW grids may be processed
    wrn_prn(sbr_nm,"Array not monotonic");
  }// endif
  return ncr;
} // end mnt_ncr()

bool
mnt_dcr(const prc_cmp * const grd,const long &grd_nbr)
{
  /* Purpose: Determine whether input array is monotonically decreasing
     Returns true when array is monotonically decreasing, false otherwise
     Routine prints an error when input array is not monotonic */
  std::string sbr_nm("mnt_dcr"); // [sng] Subroutine name
  if(grd_nbr <= 1) return true;
  const bool mnt(mnt_chk(grd,grd_nbr));
  bool dcr(false);
  if(mnt){
    if(grd[0] > grd[1]) dcr=true;
  }else{
    // fxm: 20030718 use wrn_prn() rather than err_prn() so CAM_SW/CAM_LW grids may be processed
    wrn_prn(sbr_nm,"Array not monotonic");
  }// endif
  return dcr;
} // end mnt_dcr()

prc_cmp
vec_min(const prc_cmp * const crd,const long &crd_nbr)
{
  // Purpose: Return minimum value of array
  prc_cmp min=crd[0];
  long idx; // [idx] Counting index
  for(idx=1;idx<crd_nbr;idx++){ // NB: Loop starts at 1
    if(crd[idx] < min) min=crd[idx];
  } // end loop over idx
  return min;
} // end vec_min()

prc_cmp
vec_max(const prc_cmp * const crd,const long &crd_nbr)
{
  // Purpose: Return maximum value of array
  prc_cmp max=crd[0];
  long idx; // [idx] Counting index
  for(idx=1;idx<crd_nbr;idx++){ // NB: Loop starts at 1
    if(crd[idx] > max) max=crd[idx];
  } // end loop over idx
  return max;
} // end vec_max()

long 
vec_val2idx(const prc_cmp * const crd,const long &crd_nbr,const prc_cmp &val)
{
  // Purpose: Locate index of array member closest to specified value
  // Vet input:
  if(crd_nbr < 1) err_prn("vec_val2idx()","crd_nbr < 1");
  long idx; // [idx] Counting index
  long idx_val(0);
  prc_cmp dst_new;
  prc_cmp dst_old=PRC_CMP_ABS(crd[0]-val);
  for(idx=0;idx<crd_nbr;idx++){
    if((dst_new=PRC_CMP_ABS(crd[idx]-val)) < dst_old){
      idx_val=idx;
      dst_old=dst_new;
    } // endif
  } // end loop over idx
  return idx_val;
} // end vec_val2idx()

bool // [flg] Array is monotonic
mnt_chk // [fnc] Determine monotonicity of array
(const prc_cmp * const grd, // [frc] Array
 const long &grd_nbr) // [nbr] Number of elements
{
  // Purpose: Determine whether input array is monotonic
  // Returns true when array is monotonic, false otherwise

  bool ncr(true); // [bln] Array monotonically increases
  bool dcr; // [bln] Array monotonically decreases
  bool mnt; // [bln] Array is monotonic
  long idx(0); // [idx] Counting index
  
  if(grd_nbr <= 1) return true;

  if(grd[1]-grd[0] < 0.0) ncr=false; // [bln] Array monotonically increases
  dcr=!ncr; // [bln] Array monotonically decreases

  // First check is re-done for completeness
  if(ncr){
    for(idx=1;idx<grd_nbr;idx++){
      if(grd[idx]-grd[idx-1] < 0.0) break;
    } // end loop over idx
  }else{ // not ncr
    for(idx=1;idx<grd_nbr;idx++){
      if(grd[idx]-grd[idx-1] > 0.0) break;
    } // end loop over idx
  } // not ncr

  if(idx == grd_nbr){
    mnt=true; // [bln] Array is monotonic
  }else{
    mnt=false; // [bln] Array is monotonic
    if(dbg_lvl_get() > dbg_io){
      // Print non-monotonic points
      std::cerr << "mnt_chk() reports non-monotonic array:\nBased on the first two points, this array should " << (ncr ? "increase" : "decrease") << " between the following two points:\n" << "grd[" << idx-1 << "] = " << grd[idx-1] << ", grd[" << idx << "] = " << grd[idx] << std::endl;
    } // endif dbg
  } // endif not monotonic

  return mnt; // [bln] Array is monotonic
} // end mnt_chk()

long // O [nbr] Number of consecutive monotonic elements
mnt_nbr // [fnc] Bracket monotonic portion of array
(const prc_cmp * const grd, // I [frc] Array
 const long &grd_nbr, // I [nbr] Number of elements
 const long &idx_srt) // I [idx] Index of starting element (optional)
{
  /* Purpose: Determine number of consecutive monotonic elements
     Count begins at idx_srt=0 by default, though any idx_srt < grd_nbr is valid
     Returns greatest number of elements which can be treated as monotonic array
     where elements are counted from element idx_srt
     Since one element array is always monotonic, minimal return value is 1 
     Minimal return value is 2 for input arrays of size 2 or greater
     Algorithm adapted from mnt_chk() to return number rather than flag */

  bool ncr(true); // [bln] Array monotonically increases
  long idx(idx_srt); // [idx] Counting index
  
  assert(idx_srt >= 0);
  if(grd_nbr-idx_srt == 1) return 1;

  if(grd[idx_srt+1]-grd[idx_srt] < 0.0) ncr=false; // [bln] Array monotonically increases
  // Re-do first check for completeness
  if(ncr){
    for(idx=idx_srt+1;idx<grd_nbr;idx++){
      if(grd[idx]-grd[idx-1] < 0.0) break;
    } // end loop over idx
  }else{ // not ncr
    for(idx=idx_srt+1;idx<grd_nbr;idx++){
      if(grd[idx]-grd[idx-1] > 0.0) break;
    } // end loop over idx
  } // not ncr

  // On exit, idx contains number (not index) of consecutive monotonic elements
  return idx-idx_srt; // [nbr] Number of consecutive monotonic elements
} // end mnt_nbr()

void
vec_set // [fnc] Broadcast scalar value to vector
(prc_cmp *vec, // O [frc] Vector to fill with value
 const long vec_nbr, // I [nbr] Vector size
 const prc_cmp val) // I [frc] Value to fill vector with
{
  // Purpose: Set a vector vec of vec_nbr real's to the value val
  // NB: Routine has same interface as CCM:control/resetr()
  // Local
  long idx; // [idx] Counting index
  for(idx=0;idx<vec_nbr;idx++){
    vec[idx]=val;
  } // end loop over idx
} // end vec_set()

int // O [enm] Return success code
vec_set // [fnc] Fill complex vector with real, imaginary components
(const long vec_nbr, // I [nbr] Vector size
 const prc_cmp *val_rl, // I [frc] Real array
 const prc_cmp *val_img, // I [frc] Imaginary array
 std::complex<prc_cmp> *vec) // O [frc] Complex array
{
  /* Purpose: Fill complex vector with real, imaginary components
     Function is overloaded */
  // Local
  int rcd(0); // O [enm] Return success code
  long idx; // [idx] Counting index
  for(idx=0;idx<vec_nbr;idx++) vec[idx]=std::complex<prc_cmp>(val_rl[idx],val_img[idx]);
  return rcd; // [enm] Return success code
} // end vec_set()

void
whenfgt(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr)
{
  // Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) > val
  // Function arguments are in same order as in CCM:srchutil/whenfgt()
  long idx; // [idx] Counting index
  // Initialize counter
  vld_nbr=0;
  for(idx=0;idx<crd_nbr;idx++){
    if(crd[idx] > val) vld_idx[vld_nbr++]=idx;
  } // end loop over idx
} // end whenfgt()

void
whenflt(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr)
{
  // Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) < val
  // Function arguments are in same order as in CCM:srchutil/whenflt()
  long idx; // [idx] Counting index
  // Initialize counter
  vld_nbr=0;
  for(idx=0;idx<crd_nbr;idx++){
    if(crd[idx] < val) vld_idx[vld_nbr++]=idx;
  } // end loop over idx
} // end whenflt()

void
whenfvlt(const long crd_nbr,const prc_cmp *crd,const prc_cmp *val,long *vld_idx,long &vld_nbr)
{
  // Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) < val(vld_idx)
  long idx; // [idx] Counting index
  // Initialize counter
  vld_nbr=0;
  for(idx=0;idx<crd_nbr;idx++){
    if(crd[idx] < val[idx]) vld_idx[vld_nbr++]=idx;
  } // end loop over idx
} // end whenfvlt()

void
whenfeq(const long crd_nbr,const prc_cmp *crd,const prc_cmp val,long *vld_idx,long &vld_nbr)
{
  // Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) == val
  // Function arguments are in same order as in CCM:srchutil/whenfeq()
  long idx; // [idx] Counting index
  // Initialize counter
  vld_nbr=0;
  for(idx=0;idx<crd_nbr;idx++){
    if(crd[idx] == val) vld_idx[vld_nbr++]=idx;
  } // end loop over idx
} // end whenfeq()

std::string // O [sng] Parsed command line
cmd_ln_sng // [fnc] Re-construct command line from arguments
(const int argc, // I [nbr] Argument count
 const char * const * const argv) // I [sng] Command line argument values
{
  /* Purpose: Re-construct command line from argument list and number
     Routine is identical to C version in NCO:csz.c
     Having C++ version allows C++ programs to be independent of NCO */

  std::string cmd_ln; // [sng] Parsed command line
  std::string arg_sng; // [sng] String representation of current optarg, if any
  int idx; // [idx] Index

  if(argc > 0){
    // Initialize with program name (presumably)
    cmd_ln=argv[0]; // [sng] Parsed command line
    for(idx=1;idx<argc;idx++){
      if(argv[idx]) arg_sng=argv[idx]; // Change C string into C++ string
      // Add single spaces between arguments
      cmd_ln+=" "+arg_sng;
    } // end loop over args
  } // end if

  return cmd_ln; // [sng] Parsed command line
} // end cmd_ln_sng()

std::string // O [sng] Full file name = drc/fl_nm
drc_pfx // [fnc] Intelligently return drc/fl_nm
(const std::string drc, // I [nbr] Directory to prepend
 const std::string fl_nm) // I [sng] File name
{
  /* Purpose: fl_nm := drc/fl_nm
     Routine is similar Fortran version in sng_f77.F:drcpfx()
     1. Store result in second string not first string
     2. Trim strings before concatenation
     3. Do not alter filenames which already contain slashes
     Usage:
     fl_nm=drc_pfx(drc,fl_nm) // [sbr] fl_nm := drc/fl_nm */
  const std::string sbr_nm("drc_pfx"); // [sng] Subroutine name

  // If there is no directory to prepend then return filename unaltered
  if(drc.length() <= 0L) return fl_nm;

  // If filename already contains slash then return filename unaltered
  if(fl_nm.find("/") != std::string::npos) return fl_nm;

  // Otherwise prepend directory name to filename
  // Initialize output to filename
  if(fl_nm.length() == 0) err_prn(sbr_nm,"fl_nm has zero length");
  std::string drc_fl_nm(fl_nm); // O [sng] Full file name = drc/fl_nm
  if(drc.substr(drc.length()-1L,1L) != "/") drc_fl_nm="/"+drc_fl_nm; // O [sng] Full file name = drc/fl_nm
  drc_fl_nm=drc+drc_fl_nm; // O [sng] Full file name = drc/fl_nm

  return drc_fl_nm; // O [sng] Full file name = drc/fl_nm
} // end drc_pfx()
