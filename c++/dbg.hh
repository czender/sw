// $Id$ 

/* Purpose: Systematic nomenclature for debugging/verbosity levels
   These debugging variables are used nearly globally
   Therefore complete dbg namespace is declared accessible by all 
   compilation units (files) which #include this dbg.hh header 
   dbg.hh is used to define global typedefs such as prc_cmp
   Thus dbg.hh must be first file in "personal header" section */

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Usage:
#include <dbg.hh> // Debugging constants
using namespace dbg; // [nms] Debugging namespace
*/

#ifndef DBG_HH // Contents have not yet been inserted in current source file  
#define DBG_HH

/* Pointer initialization: C++ does not define NULL, which is a C feature
   Hence C++ pointers are often initialized to 0 instead 
   CEWI_NULL replaces these zeros with a named token 
   As its name implies, CEWI_NULL is used to quiet compiler warnings
   To keep this CEWI useful and traceable, we use a token not a variable
   To make it a useful variable, there would need to be one of every type, e.g.,
   CEWI_flt_ptr, CEWI_dbl_ptr, etc. */
#define CEWI_NULL 0 

/* Precision: Either define PRC_FLT (default) (reals are usually 32 bits)...
   ...or define PRC_DBL (reals are usually 64 bits)...
   ...or define PRC_LDB (reals are 64, 96, or 128 bits)...
   _pcs stands for "precision computational suffix"
   Corresponding PRC_CMP_DCM_PLC (computational precision decimal places) are 7, 15, 22
   cout.setprecision(PRC_CMP_DCM_PLC) then prints numbers with all significan decimal places 
   for computational precision */
#if (! defined PRC_FLT) && (! defined PRC_DBL) && (! defined PRC_LDB)
#define PRC_FLT
#endif // PRC_FLT, PRC_DBL, PRC_LDB
#ifdef PRC_FLT
typedef float prc_cmp; // [prc] Computational precision is float (usually 4 bytes)
typedef double prc_cmp_min_dbl; // [prc] Computational precision at least double
/* Double ## sign is standard pre-processor syntax to join tokens
   Append "f" to literal constant to signify float */
#define PRC_CMP(x) x##f
#define PRC_CMP_MIN_DBL(x) x
//#define PRC_CMP_ABS std::fabsf
#define PRC_CMP_ABS std::fabs
#define PRC_CMP_DCM_PLC 7
#endif // !PRC_FLT
#ifdef PRC_DBL
typedef double prc_cmp; // [prc] Computational precision is double (usually 8 bytes)
typedef double prc_cmp_min_dbl; // [prc] Computational precision at least double
// Literal constants (naked constants) without specifiers default to double
#define PRC_CMP(x) x
#define PRC_CMP_MIN_DBL(x) x
#define PRC_CMP_ABS std::fabs
#define PRC_CMP_DCM_PLC 14
#endif // !PRC_DBL
#ifdef PRC_LDB
typedef long double prc_cmp; // [prc] Computational precision is long double (usually 8, 12, or 16 bytes)
typedef long double prc_cmp_min_dbl; // [prc] Computational precision at least double
/* Double ## sign is standard pre-processor syntax to join tokens
   Append "l" to literal constant to signify long double */
#define PRC_CMP(x) x##l
#define PRC_CMP_MIN_DBL(x) x##l
#define PRC_CMP_ABS std::fabsl
#define PRC_CMP_DCM_PLC 21
#endif // !PRC_DBL

// NB: dbg.hh (and dbg.cc) is natural place to define dbg_lvl_get() commands
//extern "C" {
//  unsigned short dbg_lvl=0; // [enm] Debugging level
//  unsigned short dbg_lvl_get(void){return dbg_lvl;} // end dbg_lvl_get()
//} // end extern C
  
//std::string dbg_sng(""); // [sng] Debugging string
//inline std::string dbg_sng_get(void){return dbg_sng;} // [sng] Debugging string
//inline void dbg_sng_set(std::string dbg_sng_arg){dbg_sng=dbg_sng_arg;} // [sng] Debugging string

namespace dbg{ // [nms] Debugging namespace
  // Initialization constants to be overridden before usage
  const prc_cmp cmd_ln_dfl(1.0e36); // Default for command line override values

  // Purpose of CLIP_int is to be replaced by C99-compliant compound literal, e.g.,
  // const int *dmn_scl=(const int[]){0}; // CLIP Compound literal initialization placeholder

  const char CEWI_schr(-127); // Compiler Error Warning Initializer for signed char
  const double CEWI_dbl(9.9692099683868690e+36); // Compiler Error Warning Initializer for double
  const float CEWI_flt(9.9692099683868690e+36f); // Compiler Error Warning Initializer for float 
  const int CEWI_int(-2147483647); // Compiler Error Warning Initializer for int
  const long CEWI_lng(-2147483647L); // Compiler Error Warning Initializer for long
#ifdef HAVE_LONG_LONG
  const long long CEWI_lng_lng(-2147483647LL); // Compiler Error Warning Initializer for long long
#endif // !HAVE_LONG_LONG
  const prc_cmp CEWI_cpv=PRC_CMP(9.9692099683868690e+36); // Compiler Error Warning Initializer for computational precision
  const short CEWI_sht(-32767); // Compiler Error Warning Initializer for short
  const size_t CEWI_szt(0); // Compiler Error Warning Initializer for size_t
  const unsigned char CEWI_uchr(0); // Compiler Error Warning Initializer for unsigned char
  // Defining CEWIs for complex requires #include'ing <complex>
  //  const std::complex<prc_cmp> CEWI_cpx_cpv=static_cast<std::complex<prc_cmp> >(9.9692099683868690e+36,9.9692099683868690e+36);
  
  // fxm: Declare dbg_* as signed chars to reduce storage to 1-byte each?
  // fxm: Enumerate these instead?
  // Debugging levels:
  const short dbg_nbr(9); // Number of different debugging levels
  
  const short dbg_off(0); // Production mode. Debugging is turned off.
  const short dbg_fl(1); // Filenames
  const short dbg_scl(2); // Scalars
  const short dbg_crr(3); // Current task
  const short dbg_sbr(4); // Subroutine names on entry and exit
  const short dbg_io(5); // Subroutine I/O
  const short dbg_vec(6); // Entire vectors
  const short dbg_vrb(7); // Everything
  const short dbg_old(8); // Old debugging blocks not used anymore
  
} // end debugging namespace dbg

/* Namespace dbg is required by most compilation units (files)
   Implementing "using" command here, immediately after namespace definition,
   avoids need to place command separately at top of each source file */
using namespace dbg; // [nms] Debugging namespace

#endif // DBG_HH  
