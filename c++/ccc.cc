// $Id$

// Purpose: (1) C++ template (2) System testing (3) Library driver

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text
   The original author of this software, Charlie Zender, seeks to improve
   it with your suggestions, contributions, bug-reports, and patches.
   Charlie Zender <zender at uci dot edu>
   Department of Earth System Science
   University of California, Irvine
   Irvine, CA 92697-3100 */

/* Compilation: 
   NB: libnco_c++.a must be built with same compiler
   NB: CPP token controls OpenMP activation so ccache may use non-OpenMP objects
   make -W ccc.cc OPTS=D VRS_SNG=2.6.0 ccc
   make OPTS=D ccc
   // With netCDF3
   cd ${HOME}/sw/c++;make -W ccc.cc OPTS=D ccc;cd -
   cd ${HOME}/nco/src/nco_c++;make -f Makefile.old lib_cln cln;make -f Makefile.old OPTS=D lib;cd -
   // With netCDF4
   spectral.ess.uci.edu (as of 20240808):
   cd ~/sw/c++;make CPPFLAGS="-DABORT_ON_ERROR -DHAVE_LONG_LONG -DPRC_DBL -DLINUX -I ${HOME}/include -I/opt/netcdf/include -I/opt/homebrew/include" CFLAGS="-g" LDFLAGS="-L/opt/netcdf/lib -L/opt/homebrew/lib -L${HOME}/lib -lcsz_c++ -lcsm_c++ -lnco_c++ -lnetcdf -lomp";cd - # Use -std=c++17 library
   cd ${HOME}/sw/c++;make -W ccc.cc OPTS=D OMP=N NETCDF4=Y NETCDF_INC=/opt/netcdf/include NETCDF_LIB=/opt/netcdf/lib ccc;cd -
   cd ${HOME}/sw/c++;make -W ccc.cc OPTS=D OMP=N NETCDF4=Y ccc;cd -
   cd ${HOME}/nco/src/nco_c++;make -f Makefile.old OMP=Y NETCDF4=Y lib_cln cln;make -f Makefile.old OMP=Y OPTS=D NETCDF4=Y lib;cd -
   cd ${HOME}/nco/src/nco_c++;make -f Makefile.old OMP=N NETCDF4=Y lib_cln cln;make -f Makefile.old OMP=N OPTS=D NETCDF4=Y lib;cd -
   cd ${HOME}/sw/c++;make lib_cln cln;make --jobs=1 OPTS=D ccc;cd -
   cd ${HOME}/mie;make lib_cln cln;make OPTS=D mie;cd -
   scp ~/sw/c++/ccc.cc sand.ess.uci.edu:sw/c++
   scp ~/sw/c++/ccc.cc esmf.ess.uci.edu:sw/c++
   scp ~/sw/c++/ccc.cc goldhill.cgd.ucar.edu:sw/c++ */

// etags ~/sw/c++/*.cc ~/sw/c++/*.hh ~/sw/mie/*.cc ~/sw/mie/*.hh ~/sw/slr_spc/*.cc ~/sw/slr_spc/*.hh ~/sw/ck/htrn.c

/* Usage:
   ccc --dbg=3 --function=sin --flt=1.5
   ccc --fl_in=${HOME}/nco/data/in.nc --fl_out=${HOME}/sw/c++/foo.nc
   ccc --drc_in=${HOME}/nco/data --fl_in=in.nc --drc_out=${HOME}/sw/c++ --fl_out=foo.nc
   ccc --drc_dat=${DATA}/aca --drc_in=${HOME}/nco/data --fl_in=in.nc --drc_out=${HOME}/sw/c++ --fl_out=foo.nc
   ccc --dbg=3 --xtr_LHS="xtr_prt_wgt+xtr_fll_nil+xtr_vrb" */

// Standard C++ headers
#include <complex> // Standard C++ complex class
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr
#include <limits> // STL limits numeric_limits
#include <map> // STL multimap and map classes
#include <new> // Standard C++ new handler set_new_handler()
#include <sstream> // Standard C++ string stream processing
#include <string> // Standard C++ string class
#include <thread> // Standard C++-11 thread class
#include <typeinfo> // Standard C++ header for typeid, type_info
#include <valarray> // STL valarray class template

// Internationalization i18n
// libintl.h header is required in every file with gettext() functions
#if ( !defined ALPHA ) && ( !defined MACOS ) && ( !defined SGI6 ) && ( !defined SGI64 ) && ( !defined SGIMP64 )
# include <libintl.h> // Internationalization i18n
#endif // !OS has libintl.h
// Linux and SGI libintl.h load locale.h themselves, but Solaris libintl.h does not so do it manually
#include <locale.h> // Locale setlocale()
#ifndef _LIBINTL_H
// Define stub for gettext to allow compiling when libintl.h is not available
# define gettext(foo) foo
#endif // !_LIBINTL_H

// Standard C headers
#include <cfloat> // Floating point representation, FLT_MAX, DBL_EPSILON...
#include <climits> // Integer representation, INT_MIN, INT_MAX...
#include <cmath> // sin cos cos sin 3.14159
//#include <cstdio> // stderr, EOF, FILE, NULL, etc.
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv
#include <cstring> // strcmp...
#include <ctime> // Machine time

#if !(defined __xlC__) && !(defined SGIMP64) // C++ compilers that do not allow stdint.h
# include <stdint.h> // Required by g++ for LLONG_MAX, ULLONG_MAX, by icpc for int64_t
#endif // C++ compilers that do not allow stdint.h
#include <unistd.h> // All sorts of POSIX stuff 
// #include <sys/time.h> 
// #include <sys/types.h> // 
#include <sys/resource.h> // Resource usage and limits

// 3rd party vendors
#ifdef ENABLE_MPI
# include <mpi.h> // MPI definitions
#endif // !ENABLE_MPI
#ifdef _OPENMP
# include <omp.h> // OpenMP pragmas
#endif // !_OPENMP
#include <gsl/gsl_errno.h> // GNU Scientific Library error handling
#include <gsl/gsl_sf_airy.h> // GNU Scientific Library special functions Airy functions
#include <gsl/gsl_sf_bessel.h> // GNU Scientific Library special functions Bessel functions
#include <gsl/gsl_sf_ellint.h> // GNU Scientific Library special functions Elliptic integrals
#include <gsl/gsl_sf_expint.h> // GNU Scientific Library special functions Exponential integrals
#include <gsl/gsl_sf_erf.h> // GNU Scientific Library special functions error functions
#include <gsl/gsl_sf_gamma.h> // GNU Scientific Library special functions gamma functions
#include <gsl/gsl_sf_legendre.h> // GNU Scientific Library special functions Legendre functions
#include <gsl/gsl_ieee_utils.h> // GNU Scientific Library IEEE utilities
#include <netcdf.h> // netCDF C interface

// Personal headers
// #define MAIN_PROGRAM_FILE MUST precede #include ccc.hh
#define MAIN_PROGRAM_FILE
#include "ccc.hh" // Program-specific definitions
#include "cls.hh" // Program-specific class definitions
#include <libcsz_c++.hh> // Personal C++ library
#include <libcsm_c++.hh> // Climate systems modeling library
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Namespaces
//using namespace std; // [nms] std is namespace for all standard C++ headers
namespace nmspc1{
  const prc_cmp pi(3.0); // [frc] 3.0
} // !namespace nmspc1
namespace nmspc2{
  const prc_cmp pi(4.0); // [frc] 4.0
} // !namespace nmspc2

// Important to #define personal definitions after system #include's
#define NBR_LMN(array) (sizeof(array)/sizeof(array[0]))

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map
typedef std::map<std::string,var_mtd_sct,std::less<std::string> > sng2var_mtd_map; // String-to-var_mtd_sct map

// Create mineralogy structure
typedef struct{ // [sct] mnr_tst_sct
  std::string abb; // [sng] Abbreviation
  std::string mlc; // [sng] Molecule
  prc_cmp mmw; // [kg mol-1] Mean molecular weight
  prc_cmp dns; // [kg m-3] Bulk density
} mnr_tst_sct;
typedef std::map<std::string,mnr_tst_sct,std::less<std::string> > sng2mnr_tst_map; // String-to-mnr_tst_sct map

// Prototypes
extern "C" {
#include "getopt.h" // GNU getopt() functionality from BSD my_getopt()
} // !extern

int main(int argc,char **argv)
{
  // Prototypes  
  void usg_prn(const char *opt_sng); // [fnc] Print correct usage of program
  void my_new_handler(void); // [fnc] Handler for bad_alloc()'s thrown by new[]

  // Environment
  const std::string nvr_DATA((std::getenv("DATA")) ? std::getenv("DATA") : ""); // [sng] Environment variable DATA
  const std::string nvr_DATA_RT((std::getenv("DATA_RT")) ? std::getenv("DATA_RT") : ""); // [sng] Environment variable DATA_RT
  const std::string nvr_HOME((std::getenv("HOME")) ? std::getenv("HOME") : ""); // [sng] Environment variable HOME
  const std::string nvr_HOSTNAME((std::getenv("HOSTNAME")) ? std::getenv("HOSTNAME") : ""); // [sng] Environment variable HOSTNAME
  const std::string nvr_USER((std::getenv("USER")) ? std::getenv("USER") : ""); // [sng] Environment variable USER

  // Locals
  long idx; // [idx] Counting index
  long bnd_idx; // [idx] Counting index for bnd
  long sz_idx; // [idx] Counting index for sz
  long lat_idx; // [idx] Counting index for lat
  long lon_idx; // [idx] Counting index for lon
  const int prc_cmp_dcm_plc(PRC_CMP_DCM_PLC); // [nbr] Decimal places to right of decimal to represent prc_cmp
  int rcd(0); // [rcd] Return success code
  int thr_nbr(0); // [nbr] Number of threads

  const std::string CVS_Date("$Date$"); // [sng] CVS date string
  const std::string CVS_Header("$Id$"); // [sng] CVS header string
  const std::string CVS_Id("$Id$"); // [sng] CVS identification string
  const std::string CVS_Revision("$Revision$"); // [sng] CVS revision string
  const std::string date_cvs(CVS_Date.length() > 7 ? CVS_Date.substr(7,19) : "Unknown"); // [sng] Date from CVS
  const std::string sbr_nm("main"); // [sng] Subroutine name
  const std::string vrs_cvs(CVS_Revision.length() > 10 ? CVS_Revision.substr(10,4) : "Unknown"); // [sng] Version from CVS

  /* C pre-processor macros for instantiating variable values with string tokens
     Macros for token pasting described at http://www.parashift.com/c++-faq-lite
     Layer of indirection is required, use public macro to call private macro */
#define TKN2SNG_PRV(x) #x
#define TKN2SNG(x) TKN2SNG_PRV(x)
  const std::string date_cpp(__DATE__); // [sng] Date from C pre-processor
  const std::string time_cpp(__TIME__); // [sng] Time from C pre-processor
  const std::string bufsiz_cpp(TKN2SNG(BUFSIZ)); // [sng] Version from C pre-processor
  const std::string vrs_cpp(TKN2SNG(VERSION)); // [sng] Version from C pre-processor
  const std::string hst_cpp(TKN2SNG(HOSTNAME)); // [sng] Hostname from C pre-processor
  const std::string usr_cpp(TKN2SNG(USER)); // [sng] Hostname from C pre-processor

  // Start clock and save command line
  const std::time_t time_crr_time_t(std::time((std::time_t *)NULL)); // [tm] Current date and time
  const std::string time_bfr_srt(std::ctime(&time_crr_time_t)); // [sng] Current date and time
  std::cerr << "\tStart = " << time_bfr_srt;
  prg_nm=std::strrchr(argv[0],'/'); // [sng] Program name
  if(prg_nm == NULL) prg_nm=argv[0]; else ++prg_nm; // [sng] Program name
  if(vrs_cvs == "Unknown") std::cerr << prg_nm << " version " << vrs_cpp << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  if(vrs_cvs != "Unknown") std::cerr << prg_nm << " version " << vrs_cvs << " last modified " << date_cvs << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  const std::string cmd_ln(cmd_ln_sng(argc,argv)); // [sng] Parsed command line
  std::cout << prg_nm << ": INFO Command line = " << cmd_ln << std::endl;
  
  // Set defaults for command line options 
  // Option name is variable name, e.g., --lng_foo=3, unless otherwise indicated
  Xtr_cls xtr_LHS("xtr_prt_wgt+xtr_fll_nil+xtr_vrb"); // [msk] Extrapolation flags for LHS
  Xtr_cls xtr_RHS("xtr_prt_wgt+xtr_fll_nil+xtr_vrb"); // [msk] Extrapolation flags for RHS
  bool flg_flg(false); // [flg] Temporary flag variable
  double dbl_foo(0.0); // [frc] Intrinsic double temporary variable
  double gsl_a(0.0); // [frc] Parameter a of GSL evaluations
  double gsl_b(0.0); // [frc] Parameter b of GSL evaluations
  unsigned int gsl_uint(1); // [nbr] Unsigned int for GSL evaluations
  float flt_foo(0.0f); // [frc] Intrinsic float temporary variable
  prc_cmp cpv_foo(PRC_CMP(0.0)); // [frc] Intrinsic computational precision temporary variable
  std::complex<prc_cmp> prm_a(PRC_CMP(1.0),PRC_CMP(0.0)); // [frc] Parameter a of equation
  std::complex<prc_cmp> prm_b(PRC_CMP(0.0),PRC_CMP(0.0)); // [frc] Parameter b of equation
  std::complex<prc_cmp> prm_c(PRC_CMP(-1.0),PRC_CMP(0.0)); // [frc] Parameter c of equation
  int arr_nbr; // [nbr] Number of modes
  int int_foo(0); // [nbr] Intrinsic int temporary variable
  int fl_out_fmt(NCO_FORMAT_UNDEFINED); // [enm] Output file format
  long bnd_nbr(1L); // [nbr] Number of bands
  long lat_nbr(1L); // [nbr] Number of latitudes
  long lng_foo(0L); // [nbr] Intrinsic long temporary variable
  unsigned char ubyte_foo(0); // [nbr] Intrinsic unsigned char temporary variable
#ifdef HAVE_LONG_LONG
  long long lng_lng_foo(0LL); // [nbr] Intrinsic long long temporary variable
  unsigned long long ulng_lng_foo(0ULL); // [nbr] Intrinsic long long temporary variable
#endif // !HAVE_LONG_LONG
  long lon_nbr(1L); // [nbr] Number of longitudes
  long sz_nbr(1L); // [nbr] Number of sizes
  std::string dmn_rcd(""); // [sng] Record dimension name
  std::string eqn_sng("lnstr_rat"); // [sng] Name of equation to solve
  std::string fl_err("foo.stderr"); // [sng] Error file
  std::string fl_in("in.nc"); // [sng] Input file
  std::string fl_out("foo.nc"); // [sng] Output file
  std::string fnc_sng("sin"); // [sng] Function name
  std::string lat_nm("lat"); // [sng] Name of latitude-like dimension
  std::string lon_nm("lon"); // [sng] Name of longitude-like dimension
  std::string sng_foo(""); // [sng] Intrinsic string temporary variable
  std::string slt_sng("$1$charlie"); // [sng] Salt for UNIX crypt function
  std::string tst_sng(""); // [sng] Name of test to perform

  // Derived fields
  std::string drc_dat((nvr_DATA_RT.length() > 0) ? nvr_DATA_RT : "/data/zender/aca"); // [sng] Data directory
  std::string drc_in((nvr_HOME.length() > 0) ? nvr_HOME+"/nco/data" : "/home/zender/nco/data"); // [sng] Input directory
  std::string drc_out((nvr_HOME.length() > 0) ? nvr_HOME+"/sw/c++" : ""); // [sng] Output directory

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}
       If flag is non-zero, getopt_long() returns zero and flag is set to val
       If flag is zero, getopt_long() returns contents of val */
    // Long options with no argument, no short option counterpart
    {"flg_flg",no_argument,0,0}, // [flg] Flag flag
    // Long options with argument, no short option counterpart
    {"arr_val",required_argument,0,0}, // [frc] Mode parameters
    {"bnd_nbr",required_argument,0,0}, // [nbr] Number of bands
    {"data_path",required_argument,0,0}, // [sng] Data directory
    {"dbl_foo",required_argument,0,0}, // [nbr] Intrinsic double temporary variable 
    {"dmn_rcd",required_argument,0,0}, // [sng] Record dimension name
    {"drc_dat",required_argument,0,0}, // [sng] Data directory
    {"drc_in",required_argument,0,0}, // [sng] Input directory
    {"drc_out",required_argument,0,0}, // [sng] Output directory
    {"eqn_sng",required_argument,0,0}, // [sng] Name of equation to solve
    {"fl_fmt",required_argument,0,0}, // [enm] Output file format
    {"file_format",required_argument,0,0}, // [enm] Output file format
    {"function",required_argument,0,0}, // [sng] Function name
    {"gsl_a",required_argument,0,0}, // [frc] Parameter a of GSL evaluations
    {"gsl_b",required_argument,0,0}, // [frc] Parameter b of GSL evaluations
    {"gsl_uint",required_argument,0,0}, // [nbr] Unsigned int for GSL evaluations
    {"int_foo",required_argument,0,0}, // [nbr] Intrinsic int temporary variable
    {"lat_nbr",required_argument,0,0}, // [nbr] Number of latitudes
    {"lat_nm",required_argument,0,0}, // [sng] Name of latitude-like dimension
    {"lng_foo",required_argument,0,0}, // [nbr] Intrinsic long temporary variable 
#ifdef HAVE_LONG_LONG
    {"lng_lng_foo",required_argument,0,0}, // [nbr] Intrinsic long long temporary variable 
    {"ulng_lng_foo",required_argument,0,0}, // [nbr] Intrinsic unsigned long long temporary variable 
#endif // !HAVE_LONG_LONG
    {"lon_nbr",required_argument,0,0}, // [nbr] Number of longitudes
    {"lon_nm",required_argument,0,0}, // [sng] Name of longitude-like dimension
    {"prm_a",required_argument,0,0}, // [frc] Parameter a of equation
    {"prm_b",required_argument,0,0}, // [frc] Parameter b of equation
    {"prm_c",required_argument,0,0}, // [frc] Parameter c of equation
    {"slt_sng",required_argument,0,0}, // [sng] Salt for UNIX crypt function
    {"sng_foo",required_argument,0,0}, // [sng] Intrinsic string temporary variable
    {"sz_nbr",required_argument,0,0}, // [nbr] Number of sizes
    {"tst_sng",required_argument,0,0}, // [sng] Name of test to perform
    {"ubyte_foo",required_argument,0,0}, // [nbr] Intrinsic unsigned char temporary variable 
    {"xtr_LHS",required_argument,0,0}, // [msk] Extrapolation flags for LHS
    {"xtr_RHS",required_argument,0,0}, // [msk] Extrapolation flags for RHS
    // Long options with optional argument, no short option counterpart
    // Long options with short counterparts
    {"4",no_argument,0,'4'}, // [enm] Output file format
    {"64bit",no_argument,0,'4'}, // [enm] Output file format
    {"netcdf4",no_argument,0,'4'}, // [enm] Output file format
    {"dbg_lvl",required_argument,0,'D'}, // [enm] Debugging level
    {"error",required_argument,0,'e'}, // [sng] Error file
    {"fl_in",required_argument,0,'i'}, // [sng] Input file
    {"fl_out",required_argument,0,'o'}, // [sng] Output file
    {"flt_foo",required_argument,0,'f'}, // [frc] Intrinsic float temporary variable
    {"cpv_foo",required_argument,0,'r'}, // [frc] Intrinsic computational precision temporary variable
    {"help",no_argument,0,'h'},
    {"input",required_argument,0,'i'}, // [sng] Input file
    {"output",required_argument,0,'o'}, // [sng] Input file
    {"version",no_argument,0,'v'},
    // Last option named "0" signals getopt_long() to stop processing  
    {0,0,0,0}
  }; // !opt_lng
  
  // Short options: no colon = no arg, one colon = required arg, two colons = optional arg
  const char * const opt_sht_lst("34D::f:hi:o:r:v"); // [sng] List of single-letter (C-style) option abbreviations
  extern char *optarg; // [sng] char * representation of current optarg, if any (this memory is owned by system)
  extern int optind; // [idx] extern enumerating cardinal of current option
  int opt; // [enm] Value is zero if current argument is long type, else value contains single letter version of command line argument
  int opt_idx(0); // [idx] Index of current long option into opt_lng array
  std::string opt_crr; // [sng] String representation of current long-option name
  std::string opt_sng; // [sng] String representation of current optarg, if any
 
  // Parse command line arguments 
  while(1){
    // getopt_long_only() allows a single dash '-' to prefix long options as well
    opt=getopt_long_only(argc,argv,opt_sht_lst,opt_lng,&opt_idx);
    // NB: access to opt_crr is only valid when long_opt was detected
    opt_crr=opt_lng[opt_idx].name;  
    if(optarg) opt_sng=optarg; // Change C string into C++ string
    if(opt == EOF) break; // Parse positional arguments once getopt_long_only() returns EOF
    // Process long options without short option counterparts
    if(opt == 0){
      if(dbg_lvl >= dbg_io) std::cerr << "Long option name: " << opt_crr << (optarg ? ",  Argument: "+opt_sng : ", No Argument") << std::endl;
      if(opt_crr == "flg_flg"){flg_flg=true;std::cout << "INFO: User explicitly set flg_flg to " << flg_flg << std::endl;}
      if(opt_crr == "bnd_nbr") bnd_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10);
      if(opt_crr == "sz_nbr") sz_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10);
      if(opt_crr == "arr_val"){
	char **lst_prs(char *sng_in,const char *dlm_sng,int *nbr_lst);
	// char *strdup(const char *); /* strdup() is not strict ANSI */
	// lst_prs() parses argument to list of strings, then convert to floats
	// fxm: Improve this method, should at least free memory in arg_lst_sng
	char *opt_sng_tmp=new char[opt_sng.length()+1];
	opt_sng_tmp=std::strcpy(opt_sng_tmp,opt_sng.c_str());
	/* std::string::c_str() is always const, so copy opt_sng.c_str() into C string,
	   but use std::strcpy() manually because std::strdup() is non-ANSI */
	char **arg_lst_sng=lst_prs(opt_sng_tmp,",",&arr_nbr);
	prc_cmp *arr_val=new prc_cmp[static_cast<unsigned int>(arr_nbr)]; // [frc] Mode parameters
	for(idx=0;idx<arr_nbr;idx++){
	  arr_val[idx]=static_cast<prc_cmp>(std::strtod(arg_lst_sng[idx],(char **)NULL)); // [frc]
	  std::cout << "arr_val[" << idx << "] = " << arr_val[idx] << std::endl;
	} // !idx
	// Free dynamic memory
	delete []arr_val; // [frc] Mode parameters
	delete []arg_lst_sng;
	delete []opt_sng_tmp;
      } // endif

      if(opt_crr == "dbl_foo") dbl_foo=std::strtod(opt_sng.c_str(),(char **)NULL); // [nbr] Intrinsic double temporary variable
      if(opt_crr == "dmn_rcd") dmn_rcd=opt_sng; // [sng] Record dimension name
      if(opt_crr == "drc_dat") drc_dat=opt_sng; // [sng] Data directory
      if(opt_crr == "drc_in") drc_in=opt_sng; // [sng] Input directory
      if(opt_crr == "drc_out") drc_out=opt_sng; // [sng] Output directory
      if(opt_crr == "eqn_sng") eqn_sng=opt_sng; // [sng] Name of equation to solve
      if(opt_crr == "fl_fmt" || opt_crr == "file_format") rcd=nco_create_mode_prs(opt_sng,fl_out_fmt); // [enm] Output file format
      if(opt_crr == "function") fnc_sng=SzDstFnc::opt2abb(opt_sng); // [sng] Function name 
      if(opt_crr == "gsl_a") gsl_a=std::strtod(opt_sng.c_str(),(char **)NULL); // [frc] Parameter a of GSL evaluation
      if(opt_crr == "gsl_b") gsl_b=std::strtod(opt_sng.c_str(),(char **)NULL); // [frc] Parameter b of GSL evaluation
      if(opt_crr == "gsl_uint") gsl_uint=static_cast<unsigned int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Unsigned int for GSL evaluations
      if(opt_crr == "int_foo") int_foo=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Intrinsic int temporary variable
      if(opt_crr == "lat_nbr") lat_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of latitudes
      if(opt_crr == "lat_nm") lat_nm=opt_sng; // [sng] Name of latitude-like dimension
      if(opt_crr == "lng_foo") lng_foo=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Intrinsic long temporary variable
#ifdef HAVE_LONG_LONG
      // Future C++ standard will support std::strtoll() and std::strtoull()
      if(opt_crr == "lng_lng_foo") sng2nbr(opt_sng,&lng_lng_foo); // [nbr] Intrinsic long long temporary variable
      if(opt_crr == "ulng_lng_foo") sng2nbr(opt_sng,&ulng_lng_foo); // [nbr] Intrinsic unsigned long long temporary variable
#endif // !HAVE_LONG_LONG
      if(opt_crr == "lon_nbr") lon_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of longitudes
      if(opt_crr == "lon_nm") lon_nm=opt_sng; // [sng] Name of longitude-like dimension
      if(opt_crr == "prm_a") prm_a=sng2cpx(opt_sng); // [frc] Parameter a of equation
      if(opt_crr == "prm_b") prm_b=sng2cpx(opt_sng); // [frc] Parameter b of equation
      if(opt_crr == "prm_c") prm_c=sng2cpx(opt_sng); // [frc] Parameter c of equation
      if(opt_crr == "slt_sng") slt_sng=opt_sng; // [sng] Salt for UNIX crypt function
      if(opt_crr == "sng_foo") sng_foo=opt_sng; // [sng] Intrinsic string temporary variable
      if(opt_crr == "tst_sng") tst_sng=opt_sng; // [sng] Name of test to perform
      if(opt_crr == "ubyte_foo") sng2nbr(opt_sng,&ubyte_foo); // [nbr] Intrinsic unsigned char temporary variable
      //if(opt_crr == "ubyte_foo") static_cast<unsigned char>(std::strtoul(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Intrinsic unsigned char temporary variable
      if(opt_crr == "xtr_LHS") xtr_LHS.flg_set(opt_sng); // [msk] Extrapolation options for LHS
      if(opt_crr == "xtr_RHS") xtr_RHS.flg_set(opt_sng); // [msk] Extrapolation options for RHS
    } // opt != 0
    switch(opt){
    case 0: // Long options have already been processed, return
      break;
    case '3': // Prescribe netCDF3 output storage format
      fl_out_fmt=NC_FORMAT_CLASSIC;
      break;
    case '4': // Catch-all to prescribe output storage format
      if(opt_sng == "64bit") fl_out_fmt=NC_FORMAT_64BIT; else fl_out_fmt=NC_FORMAT_NETCDF4;
      break;
    case 'D': // Debugging level (default is 0) 
      if(optarg) dbg_lvl=static_cast<unsigned short int>(std::strtoul(opt_sng.c_str(),(char **)NULL,10)); else dbg_lvl=dbg_fl;
      break;
    case 'f': // Set generic tuning parameter (default is 0.0)
      flt_foo=static_cast<float>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Intrinsic float temporary variable
      break;
    case 'r': // Set generic tuning parameter (default is 0.0)
      cpv_foo=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Intrinsic computational precision temporary variable
      break;
    case 'i': 
      fl_in=opt_sng; // [sng] Input file
      break;
    case 'h':
      (void)usg_prn(opt_sht_lst); // [fnc] Print proper usage 
      Exit_gracefully();
      break;
    case 'o':
      fl_out=opt_sng; // [sng] Output file
      break;
    case 'v': // Print CVS program info 
      std::cerr << CVS_Header << std::endl;
      Exit_gracefully();
      break;
    default:
      (void)usg_prn(opt_sht_lst); // Print proper usage 
      return EXIT_FAILURE;
    } // !opt
  } // !while(1)
  
  // Process positional arguments  
  if(optind < argc){
    int psn_arg_nbr=argc-optind; // [nbr] Positional argument number
    if(psn_arg_nbr > 2){
      err_prn(prg_nm,sbr_nm,"Too many positional arguments");
    }else if(psn_arg_nbr == 1){
      fl_out=argv[optind++]; // [sng] Output file
    }else if(psn_arg_nbr == 2){
      fl_in=argv[optind++]; // [sng] Input file
      fl_out=argv[optind++]; // [sng] Output file
    } // !psn_arb_nbr
  } // !optind

  // Compute quantities that might depend on command line input
  (void)fio::data_path_set(drc_dat); // [sng] Data directory
  // Prepend user-specified path, if any, to input data file names
  if(fl_in.length() > 0) fl_in=drc_pfx(drc_in,fl_in); // [sng] Input file
  // Prepend user-specified path, if any, to output data file names
  if(fl_out.length() > 0) fl_out=drc_pfx(drc_out,fl_out); // [sng] Output file
  // Diagnose properties of input file and filesystem
  (void)fio::fl_stat_dgn(fl_in);

  // Open input file
  int nc_id=nco_open(fl_in,NC_NOWRITE); // [fnc] Open netCDF file
  // Input required data
  const long wvl_nbr(nco_inq_dimlen(nc_id,static_cast<std::string>("wvl"))); // [nbr] Number of wavelengths
  if(dbg_lvl > dbg_off) std::cerr << "Number of wavelengths in "+fl_in+" is " << wvl_nbr << std::endl;
  /* netCDF C++ interface automatically allocates memory required by get_var()
     User must explicitly free() this memory when no longer needed
     Currently do this after writing input array to output file */
  prc_cmp *wvl; // [m] Wavelength
  rcd=nco_get_var(nc_id,static_cast<std::string>("wvl"),wvl); // [m] Wavelength
  if(dbg_lvl > dbg_off) std::cerr << "Value of wvl[0] in "+fl_in+" is " << wvl[0] << std::endl;
  rcd=nco_close(nc_id); // [fnc] Close netCDF file
  std::cerr << prg_nm_get() << ": INFO Ingested " << fl_in << std::endl;

  // Main body of code 
  if(dbg_lvl >= dbg_sbr) dbg_lvl_tst();
  //  prc_cmp tpt_d2d[bnd_nbr][sz_nbr]; // [K] Temperature (dynamic 2-dimensional)
  //  std::vector< std::vector<prc_cmp> > tpt(bnd_nbr,std::vector<prc_cmp>(sz_nbr)); // [K] Temperature (vector 2-dimensional)
  std::vector<prc_cmp> tpt_v1d(bnd_nbr*sz_nbr); // [K] Temperature (vector 1-dimensional)
  // Generate dummy data
  prc_cmp *sz=new prc_cmp[sz_nbr]; // [m] Size at bin center
  for(idx=0;idx<sz_nbr;idx++) sz[idx]=0.5e-6+idx*0.1e-6; // [m] Size at bin center

  if(dbg_lvl == dbg_off && tst_sng == "") std::cout << "HINT: Invoke with " << prg_nm << " --tst=tst_sng to see specific test or --dbg=3 to see current test" << std::endl;

  a2d_cls<prc_cmp> tpt_a2d(bnd_nbr,sz_nbr); // [K] Temperature (array two-dimensional)
  if(dbg_lvl == dbg_old || tst_sng == "a2d"){
    std::cout << "Testing a2d_cls behavior..." << std::endl;
    std::cout << "ccc --tst=a2d --bnd_nbr=5 --sz_nbr=3" << std::endl;
    std::cout << "ncks -C -v tpt_v1d,tpt_a2d ~/c++/foo.nc" << std::endl;
    tpt_a2d=0.0; // [K] Temperature (array two-dimensional)
    std::cout << "Testing ostream overload operator..." << std::endl;
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
	tpt_a2d(bnd_idx,sz_idx)=bnd_idx*sz_nbr+sz_idx; // [K] Temperature (array two-dimensional)
	std::cout << "tpt_a2d(" << bnd_idx << "," << sz_idx << ") = " << tpt_a2d(bnd_idx,sz_idx) << std::endl;
      } // !sz_idx
    } // !bnd_idx
    // Print contents using overloaded stream insertion operator
    std::cout << tpt_a2d << std::endl; // [fnc]
    // std::cout << "Access a2d_cls out of range: " << tpt_a2d(bnd_nbr*sz_nbr) << std::endl; // [fnc]
  } // !dbg || tst_sng == "a2d"

  if(dbg_lvl == dbg_old || tst_sng == "aer"){
    std::cout << "Testing aer_cls..." << std::endl;
    std::cout << "ccc --tst=aer --lng_foo=1 --sng_foo='saharan_dust'" << std::endl;
    aer_cls aer(sng_foo);
    std::cout << aer << std::endl;
    aer.tst(sng_foo,lng_foo);
  } // !dbg || tst_sng == "aer"

  if(dbg_lvl == dbg_old || tst_sng == "arr"){
    // ccc --tst=arr --arr_val=-1.0,1.0e1,1.0e36
    std::cout << "Testing arrays in command line arguments..." << std::endl;
  } // !dbg || tst_sng == "arr"

  if(dbg_lvl == dbg_old || tst_sng == "bnr"){
    /* ccc --dbg=3 --int_foo=0x01
       const int int_foo(0x01); // Hexadecimal notation */
    std::cout << "Testing internal binary representation of all types..." << std::endl;
    std::cout << "HINT: use --tst=gsl --gsl_a=3.141 to see GSL's intrinsic functions for binary representation of floating point types..." << std::endl;
    std::cout << "Seed values with --int_foo, e.g., ccc --tst=bnr --int_foo=-4 --flt_foo=255" << std::endl;
    const bool bln_bnr(true);
    const char chr_bnr(int_foo);
    const double dbl_bnr(dbl_foo);
    const float flt_bnr(flt_foo);
    const int int_bnr(int_foo);
    const long lng_bnr(int_foo);
    const short sht_bnr(int_foo);
    const signed char schr_bnr(int_foo);
    const unsigned char uchr_bnr(std::abs(int_foo));
    const unsigned int uint_bnr(std::abs(int_foo));
    const unsigned long ulng_bnr(std::abs(int_foo));
    const unsigned short usht_bnr(std::abs(int_foo));
#ifdef HAVE_LONG_LONG
    // long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
    const long long lng_lng_bnr(int_foo);
    const unsigned long long ulng_lng_bnr(std::abs(int_foo));
#endif // !HAVE_LONG_LONG
    std::cout << "Binary of float " << flt_bnr << " is " << bit_sng_tpl(flt_bnr) << std::endl;
    flt_foo=flt_bnr; // Copy memory so bitshifting exponent to byte boundary does not corrupt original
    unsigned int u32_val=*(unsigned int *)(&flt_foo); // Reinterpret as integer (reinterpret_cast also works)
    u32_val<<=1; // Bit-shift so exponent is on byte boundaries
    unsigned char u8_val=*(unsigned char *)(&u32_val); // This should be the 8-bit exponent
    int int32_val=(int)u8_val; // For printing
    std::cout << "u32_val = " << u32_val << std::endl;
    std::cout << "u8_val = " << u8_val << std::endl;
    std::cout << "int32_val = " << int32_val << std::endl;
    std::cout << "Binary of u8_val " << uchr_bnr << " is " << bit_sng_tpl(u8_val) << std::endl;
    std::cout << "Binary of double " << dbl_bnr << " is " << bit_sng_tpl(dbl_bnr) << std::endl;
    std::cout << "Binary of bool " << bln_bnr << " is " << bit_sng_tpl(bln_bnr) << std::endl;
    std::cout << "Binary of char " << chr_bnr << " is " << bit_sng_tpl(chr_bnr) << std::endl;
    std::cout << "Binary of int " << int_bnr << " is " << bit_sng_tpl(int_bnr) << std::endl;
    std::cout << "Binary of long " << lng_bnr << " is " << bit_sng_tpl(lng_bnr) << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "Binary of long long " << lng_lng_bnr << " is " << bit_sng_tpl(lng_lng_bnr) << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "Binary of short " << sht_bnr << " is " << bit_sng_tpl(sht_bnr) << std::endl;
    std::cout << "Binary of signed char " << schr_bnr << " is " << bit_sng_tpl(schr_bnr) << std::endl;
    std::cout << "Binary of unsigned char " << uchr_bnr << " is " << bit_sng_tpl(uchr_bnr) << std::endl;
    std::cout << "Binary of unsigned int " << uint_bnr << " is " << bit_sng_tpl(uint_bnr) << std::endl;
    std::cout << "Binary of unsigned long " << ulng_bnr << " is " << bit_sng_tpl(ulng_bnr) << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "Binary of unsigned long long " << ulng_lng_bnr << " is " << bit_sng_tpl(ulng_lng_bnr) << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "Binary of unsigned short " << usht_bnr << " is " << bit_sng_tpl(usht_bnr) << std::endl;
  } // !dbg || tst_sng == "bnr"
  
  if(dbg_lvl == dbg_old || tst_sng == "uint"){
    /* ccc --dbg=3 --int_foo=0x01
       const int int_foo(0x01); // Hexadecimal notation */
    std::cout << "Testing representation of all types as uint32, uint64..." << std::endl;
    std::cout << "HINT: This shows how to encode typed values as uint32/uint64, e.g., for passing as HDF5 filter parameters..." << std::endl;
    std::cout << "Seed values with --int_foo, e.g., ccc --tst=uint --int_foo=4" << std::endl;
    const bool bln_bnr(true);
    const char chr_bnr(int_foo);
    const double dbl_bnr(dbl_foo);
    const float flt_bnr(flt_foo);
    const int int_bnr(int_foo);
    const long lng_bnr(int_foo);
    const short sht_bnr(int_foo);
    const signed char schr_bnr(int_foo);
    const unsigned char uchr_bnr(std::abs(int_foo));
    const unsigned int uint_bnr(std::abs(int_foo));
    const unsigned long ulng_bnr(std::abs(int_foo));
    const unsigned short usht_bnr(std::abs(int_foo));
#ifdef HAVE_LONG_LONG
    // long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
    const long long lng_lng_bnr(int_foo);
    const unsigned long long ulng_lng_bnr(std::abs(int_foo));
#endif // !HAVE_LONG_LONG
    unsigned int ui32_val; // Target representation is four byte unsigned int
    unsigned long long int ui64_val; // Target representation is eight byte unsigned long long int
    ui32_val=*(unsigned int *)(&flt_bnr); // Reinterpret as ui32 (reinterpret_cast also works)
    std::cout << "ui32 of float " << flt_bnr << "f is " << ui32_val << "u" << std::endl;
    ui64_val=*(unsigned long long int *)(&dbl_bnr); // Reinterpret as ui64 (reinterpret_cast also works)
    std::cout << "ui64 of double " << dbl_bnr << "d is " << ui64_val << "u" << std::endl;
    ui32_val=*(unsigned int *)(&int_bnr);
    std::cout << "ui32 of int " << int_bnr << " is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of bool " << bln_bnr << "fxm is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of char " << chr_bnr << "b is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of long " << lng_bnr << " is " << ui32_val << "u" << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "ui64 of long long " << lng_lng_bnr << "l is " << ui64_val << "u" << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "ui32 of short " << sht_bnr << "s is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of signed char " << schr_bnr << "b is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of unsigned char " << uchr_bnr << "ub is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of unsigned int " << uint_bnr << "u is " << ui32_val << "u" << std::endl;
    std::cout << "ui32 of unsigned long " << ulng_bnr << "u is " << ui32_val << "u" << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "ui64 of unsigned long long " << ulng_lng_bnr << "ul is " << ui64_val << "u" << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "ui32 of unsigned short " << usht_bnr << "us is " << ui32_val << "u" << std::endl;
  } // !dbg || tst_sng == "uint"

  if(dbg_lvl == dbg_old || tst_sng == "cmp"){
    // ccc --dbg=3 --tst=cmp
    std::cout << "Testing compiler behavior..." << std::endl;
    // Purpose: Test compiler behavior
    std::string cmp_nm("Searched for known compiler tokens in main(), none found, compiler is unknown"); // [sng] Compiler name
    // Begin tokens usually found on AIX
#ifdef _AIX
    std::cout << "Token _AIX defined in main(), probably compiled on AIX" << std::endl;
#endif // !_AIX
#ifdef _AIX32
    std::cout << "Token _AIX32 defined in main(), compiled on AIX version >= 3.2" << std::endl;
#endif // !_AIX32
#ifdef _AIX41
    std::cout << "Token _AIX41 defined in main(), compiled on AIX version >= 4.1" << std::endl;
#endif // !_AIX41
#ifdef _AIX43
    std::cout << "Token _AIX43 defined in main(), compiled on AIX version >= 4.3" << std::endl;
#endif // !_AIX43
#ifdef _AIX51
    std::cout << "Token _AIX51 defined in main(), compiled on AIX version >= 5.1" << std::endl;
#endif // !_AIX51
#ifdef _AIX52
    std::cout << "Token _AIX52 defined in main(), compiled on AIX version >= 5.2" << std::endl;
#endif // !_AIX52
#ifdef _AIX53
    std::cout << "Token _AIX53 defined in main(), compiled on AIX version >= 5.3" << std::endl;
#endif // !_AIX53
#ifdef _AIX61
    std::cout << "Token _AIX61 defined in main(), compiled on AIX version >= 6.1" << std::endl;
#endif // !_AIX61
#ifdef _BIG_ENDIAN
    std::cout << "Token _BIG_ENDIAN defined in main(), CPU is big-endian" << std::endl;
#endif // !_BIG_ENDIAN
#ifdef __xlC__
    // Obtain xlC compiler version with: xlc -qversion or xlC -qversion
    cmp_nm="Token __xlC__ defined in main(), probably compiled with AIX xlC_r or xlC"; // [sng] Compiler name
#endif // !__xlC__
    // End tokens usually found on AIX
    // Begin clang tokens
#ifdef __clang
    const std::string cmp_vrs(TKN2SNG(__clang)); // [sng] Version from C pre-processor
    cmp_nm="Token __clang defined in main(), probably compiled with LLVM clang"; // [sng] Compiler name
    std::cout << "Clang compiler version from __clang is " << cmp_vrs << std::endl;
#endif // !__clang
    // End clang tokens
    // Begin GCC tokens
#if defined(__GNUC__) && !defined(__clang) && !defined(__INTEL_COMPILER) && !defined(__PATHCC__)
    cmp_nm="Token __GNUC__ defined in main(), probably compiled with GNU g++"; // [sng] Compiler name
    const std::string cmp_vrs(TKN2SNG(__VERSION__)); // [sng] Compiler version
    const std::string cmp_vrs_mjr(TKN2SNG(__GNUC__)); // [sng] Compiler major version
    const std::string cmp_vrs_mnr(TKN2SNG(__GNUC_MINOR__)); // [sng] Compiler minor version
    const std::string cmp_vrs_pch(TKN2SNG(__GNUC_PATCHLEVEL__)); // [sng] Compiler patch version
    std::cout << "GCC major version is " << cmp_vrs_mjr << std::endl;
    std::cout << "GCC minor version is " << cmp_vrs_mnr << std::endl;
    std::cout << "GCC patch version is " << cmp_vrs_pch << std::endl;
    std::cout << "GCC version is " << cmp_vrs << std::endl;
#endif // !__GNUC__
    // End GCC tokens
#ifdef __INTEL_COMPILER
    const std::string cmp_vrs(TKN2SNG(__INTEL_COMPILER)); // [sng] Version from C pre-processor
    cmp_nm="Token __INTEL_COMPILER defined in main(), probably compiled with Intel icpc"; // [sng] Compiler name
    std::cout << "Intel compiler version from __INTEL_COMPILER is " << cmp_vrs << std::endl;
#endif // !__INTEL_COMPILER
#ifdef __PATHCC__
    // Pathscale pre-defined macros in User's Guide p. 57
    cmp_nm="Token __PATHCC__ defined in main(), probably compiled with Pathscale (Qlogic) pathCC"; // [sng] Compiler name
    const std::string cmp_vrs(TKN2SNG(__VERSION__)); // [sng] Compiler version
    const std::string cmp_vrs_mjr(TKN2SNG(__PATHCC__)); // [sng] Compiler major version
    const std::string cmp_vrs_mnr(TKN2SNG(__PATHCC_MINOR__)); // [sng] Compiler minor version
    const std::string cmp_vrs_pch(TKN2SNG(__PATHCC_PATCHLEVEL__)); // [sng] Compiler minor version
    std::cout << "pathCC major version is " << cmp_vrs_mjr << std::endl;
    std::cout << "pathCC minor version is " << cmp_vrs_mnr << std::endl;
    std::cout << "pathCC patch version is " << cmp_vrs_pch << std::endl;
    std::cout << "pathCC version is " << cmp_vrs << std::endl;
#endif // !__PATHCC__
#ifdef PGI_CXX
    cmp_nm="Token PGI_CXX defined in main(), probably compiled with PGI pgCC"; // [sng] Compiler name
#endif // !PGI_CXX
    std::cout << cmp_nm << std::endl;
  } // !dbg || tst_sng == "cmp"

  if(dbg_lvl == dbg_old || tst_sng == "cpp"){
    std::cout << "Testing C pre-proprocessor behavior..." << std::endl;
    // Purpose: Test pre-processor behavior
    std::cerr << "date_cpp = " << date_cpp << std::endl;
    std::cerr << "time_cpp = " << time_cpp << std::endl;
    std::cerr << "vrs_cpp = " << vrs_cpp << std::endl;
    std::cout << "BUFSIZ during compilation = " << bufsiz_cpp << std::endl;

    std::string pvm_arch_cmp(TKN2SNG(PVM_ARCH)); // [sng] Architecture name
    if(pvm_arch_cmp == "") pvm_arch_cmp="Undefined"; // [sng] Architecture name
    if(pvm_arch_cmp == "Undefined") err_prn(sbr_nm," PVM_ARCH = "+pvm_arch_cmp+" during compilation");
    std::cout << "PVM_ARCH during compilation = " << pvm_arch_cmp << std::endl;
    std::string nvr_PVM_ARCH((std::getenv("PVM_ARCH")) ? std::getenv("PVM_ARCH") : ""); // [sng] Environment variable PVM_ARCH
    std::cout << "PVM_ARCH current environment = " << nvr_PVM_ARCH << std::endl;
    if(pvm_arch_cmp != nvr_PVM_ARCH) err_prn(sbr_nm," PVM_ARCH during compilation != PVM_ARCH in current environment");
  } // !dbg || tst_sng == "cpp"

  if(dbg_lvl == dbg_old || tst_sng == "cpx"){
    std::cout << "Testing complex number arithmetic..." << std::endl;
    std::cout << "ccc --tst=cpx --prm_a=0.70711+0.70711i" << std::endl;
    std::string sng_cpx; // [sng] String representing complex number
    sng_cpx="2.0-3.2i"; // [sng] String representing complex number
    std::cout << "sng2cpx("+sng_cpx+") = " << sng2cpx(sng_cpx) << std::endl;
    std::complex<float> fcx_foo(std::sqrt(2.0),-std::sqrt(2.0));
    std::complex<double> dcx_foo(0.0,1.2345678901234567);
    std::complex<prc_cmp> two_cpx(2.0,0.0);
    std::cout << "prm_a = " << prm_a << std::endl;
    std::cout << "two_cpx = " << two_cpx << std::endl;
    std::cout << "fcx_foo = " << fcx_foo << std::endl;
    std::cout << "dcx_foo = " << dcx_foo << std::endl;
    std::cout << "std::sqrt(fcx_foo) = " << std::sqrt(fcx_foo) << std::endl;
    std::cout << "std::sqrt(dcx_foo) = " << std::sqrt(dcx_foo) << std::endl;
    std::cout << "norm(fcx_foo) = " << norm(fcx_foo) << std::endl;
    std::cout << "norm(dcx_foo) = " << norm(dcx_foo) << std::endl;
    std::cout << "prm_a*prm_a = " << prm_a*prm_a << std::endl;
    // std::cout << "std::pow(prm_a,2.0) = " << std::pow(prm_a,2.0) << std::endl; // NB: GCC emits warning here
    std::cout << "std::pow(prm_a,two_cpx) = " << std::pow(prm_a,two_cpx) << std::endl; // NB: GCC does not emit warning here, therefore GCC wants both arguments of pow to be the same type
    std::cout << "std::exp(prm_a) = " << std::exp(prm_a) << std::endl;
    std::cout << "std::sin(prm_a) = " << std::sin(prm_a) << std::endl;
    dcx_foo=static_cast<std::complex<double> >(fcx_foo);
    std::cout << "after cast dcx_foo = " << dcx_foo << std::endl;
    // std::string sng_foo("(1.5,-0.3)");
    // fcx_foo=static_cast<std::complex<float> >(sng_foo); // ERROR: casting string to complex is not supported
    // std::cout << "after cast fcx_foo = " << dcx_foo << std::endl;
    dcx_foo=std::complex<double>(0.0,1.0);
    std::cout << "after assignment dcx_foo=std::complex<double>(0.0,1.0), dcx_foo = " << dcx_foo << std::endl;
    dcx_foo=std::complex<double>(fcx_foo);
    std::cout << "after assignment dcx_foo=std::complex<double>(fcx_foo), dcx_foo = " << dcx_foo << std::endl;
  } // !dbg || tst_sng == "cpx"
  
  if(dbg_lvl == dbg_old || tst_sng == "crp"){
    std::cout << "Testing cryptography..." << std::endl;
    std::cout << "ccc --dbg=3 --tst=crp --sng='cracker'" << std::endl;
    // Salting with two character string tells crypt() to use old UNIX-standard DES hash (insecure)
    std::cout << "ccc --dbg=3 --tst=crp --slt='cz' --sng='cracker'" << std::endl;
    // Salting with $1$ (followed by up to eight additional characters) tells crypt() to use new MD5 hash
    std::cout << "ccc --dbg=3 --tst=crp --slt='$1$charlie' --sng='cracker'" << std::endl;
    const string pwd_ucr(sng_foo); // [sng] Un-encrypted password
    string pwd_ncr; // [sng] Encrypted password
#ifndef PGI_CXX
    pwd_ncr=crypt(pwd_ucr.c_str(),slt_sng.c_str()); // [sng] Encrypted password
#else // PGI_CXX
    std::cout << "PGI compiler pgCC does not support crypt() function, skipping test..." << std::endl;
#endif // PGI_CXX
    std::cout << "pwd_ucr = " << pwd_ucr << ", salt = " << slt_sng << ", pwd_ncr = " << pwd_ncr << std::endl;
  } // !dbg || tst_sng == "crp"

  if(dbg_lvl == dbg_old || tst_sng == "ddra"){
    std::cout << "Testing DDRA model..." << std::endl;
    std::cout << "Try, e.g.,\nccc --tst=ddra" << std::endl;
    // Compute mean dimension sizes
    int dmnsz_stl_lat=2160;
    int dmnsz_stl_lon=4320;
    float dmnszavg_stl=std::sqrt(static_cast<float>(dmnsz_stl_lat*dmnsz_stl_lon));
    long varsz_stl_1D=dmnsz_stl_lat;
    long varsz_stl_2D=dmnsz_stl_lat*dmnsz_stl_lon;
    int dmnsz_gcm_scl=1;
    int dmnsz_gcm_time=8;
    int dmnsz_gcm_lev=32;
    int dmnsz_gcm_lat=128;
    int dmnsz_gcm_lon=256;
    int dmnsz_gcm_cst=1;
    int varnbr_gcm_0D=8;
    int varnbr_gcm_1D=8;
    int varnbr_gcm_2D=16;
    int varnbr_gcm_3D=64;
    int varnbr_gcm_4D=32;
    int varnbr_gcm_5D=0;
    long varsz_gcm_0D=dmnsz_gcm_scl;
    long varsz_gcm_1D=dmnsz_gcm_time;
    long varsz_gcm_2D=dmnsz_gcm_lat*dmnsz_gcm_lon;
    long varsz_gcm_3D=dmnsz_gcm_time*dmnsz_gcm_lat*dmnsz_gcm_lon;
    long varsz_gcm_4D=dmnsz_gcm_time*dmnsz_gcm_lev*dmnsz_gcm_lat*dmnsz_gcm_lon;
    long varsz_gcm_5D=dmnsz_gcm_time*dmnsz_gcm_lev*dmnsz_gcm_lat*dmnsz_gcm_lon*dmnsz_gcm_cst;
    int varnbr_gcm_ttl=varnbr_gcm_1D+varnbr_gcm_2D+varnbr_gcm_3D+varnbr_gcm_4D+varnbr_gcm_5D;
    long long lmnnbr_gcm_0D_ttl=varsz_gcm_0D*varnbr_gcm_0D;
    long long lmnnbr_gcm_1D_ttl=varsz_gcm_1D*varnbr_gcm_1D;
    long long lmnnbr_gcm_2D_ttl=varsz_gcm_2D*varnbr_gcm_2D;
    long long lmnnbr_gcm_3D_ttl=varsz_gcm_3D*varnbr_gcm_3D;
    long long lmnnbr_gcm_4D_ttl=varsz_gcm_4D*varnbr_gcm_4D;
    long long lmnnbr_gcm_5D_ttl=varsz_gcm_5D*varnbr_gcm_5D;
    long long lmnnbr_gcm_ttl=lmnnbr_gcm_0D_ttl+lmnnbr_gcm_1D_ttl+lmnnbr_gcm_2D_ttl+lmnnbr_gcm_3D_ttl+lmnnbr_gcm_4D_ttl+lmnnbr_gcm_5D_ttl;
    float dmnszavg_gcm_0D=1.0;
    float dmnszavg_gcm_1D=std::pow(varsz_gcm_1D,1.0/1.0);
    float dmnszavg_gcm_2D=std::pow(varsz_gcm_2D,1.0/2.0);
    float dmnszavg_gcm_3D=std::pow(varsz_gcm_3D,1.0/3.0);
    float dmnszavg_gcm_4D=std::pow(varsz_gcm_4D,1.0/4.0);
    float dmnszavg_gcm_5D=std::pow(varsz_gcm_4D,1.0/5.0);
    float dmnszavg_gcm_avg_wgt_varnbr=(varnbr_gcm_0D*dmnszavg_gcm_0D+varnbr_gcm_1D*dmnszavg_gcm_1D+varnbr_gcm_2D*dmnszavg_gcm_2D+varnbr_gcm_3D*dmnszavg_gcm_3D+varnbr_gcm_4D*dmnszavg_gcm_4D+varnbr_gcm_5D*dmnszavg_gcm_5D)/varnbr_gcm_ttl;
    float dmnszavg_gcm_avg_wgt_lmnnbr=(lmnnbr_gcm_0D_ttl*dmnszavg_gcm_0D+lmnnbr_gcm_1D_ttl*dmnszavg_gcm_1D+lmnnbr_gcm_2D_ttl*dmnszavg_gcm_2D+lmnnbr_gcm_3D_ttl*dmnszavg_gcm_3D+lmnnbr_gcm_4D_ttl*dmnszavg_gcm_4D+lmnnbr_gcm_5D_ttl*dmnszavg_gcm_5D)/lmnnbr_gcm_ttl;
    
    // Define constants normally in headers, e.g., enums
    short nco_op_typ_sbt=0; /* [enm] Operation type */
    short nco_op_typ_avg=1; /* [enm] Operation type */
    short fl_typ_gcm=0;
    short fl_typ_stl=1;
    
    /* Cumulative file costs */
    long long lmn_nbr_ttl=0LL; /* I/O [nbr] Cumulative variable size */
    long long ntg_nbr_ttl=0LL; /* I/O [nbr] Cumulative integer operations */
    long long flp_nbr_ttl=0LL; /* I/O [nbr] Cumulative floating point operations */
    
    /* Cumulative times */
    float tm_ntg_ttl=0.0f; /* I/O [s] Cumulative integer time */
    float tm_flp_ttl=0.0f; /* I/O [s] Cumulative floating point time */
    float tm_rd_ttl=0.0f; /* I/O [s] Cumulative read time */
    float tm_wrt_ttl=0.0f; /* I/O [s] Cumulative write time */
    float tm_io_ttl=0.0f; /* [s] I/O time */
    float tm_ttl=0.0f; /* I/O [s] Cumulative time */
    
    /* Current variable costs */
    float tm_ntg; /* [s] Integer time */
    float tm_flp; /* [s] Floating point time */
    float tm_rd; /* [s] Read time */
    float tm_wrt; /* [s] Write time */
    float tm_io; /* [s] I/O time */
    float tm_crr; /* [s] Time for this variable */
    long long ntg_nbr; /* [nbr] Integer operations */
    long long flp_nbr; /* [nbr] Floating point operations */
    long long rd_nbr_byt=CEWI_int; /* [B] Bytes read */
    long long wrt_nbr_byt=CEWI_int; /* [B] Bytes written */
    
    /* Summary statistics */
    float tm_frc_flp_ttl; /* [frc] Floating point time fraction */
    float tm_frc_io_ttl; /* [frc] I/O time fraction */
    float tm_frc_ntg_ttl; /* [frc] Integer time fraction */
    float tm_frc_rd_ttl; /* [frc] Read time fraction */
    float tm_frc_wrt_ttl; /* [frc] Write time fraction */
    
    /* Default algorithm costs if invoked */
    long long ntg_nbr_byt_swp_dfl; /* [nbr] Integer operations for byte-swap */
    long long ntg_nbr_brd_dfl; /* [nbr] Integer operations for broadcasting */
    long long ntg_nbr_clc_dfl; /* [nbr] Integer operations for collection */
    long long ntg_nbr_rdc_dfl; /* [nbr] Integer operations for reduction */
    long long ntg_nbr_nrm_dfl; /* [nbr] Integer operations for normalization */
    long long flp_nbr_bnr_dfl; /* [nbr] Floating point operations for binary arithmetic */
    long long flp_nbr_rdc_dfl; /* [nbr] Floating point operations for reduction */
    long long flp_nbr_nrm_dfl; /* [nbr] Floating point operations for normalization */
    
    /* Initialize all algorithm counts for this variable to zero then increment */
    long long ntg_nbr_byt_swp=0LL; /* [nbr] Integer operations for byte-swap */
    long long ntg_nbr_brd=0LL; /* [nbr] Integer operations for broadcasting */
    long long ntg_nbr_clc=0LL; /* [nbr] Integer operations for collection */
    long long ntg_nbr_rdc=0LL; /* [nbr] Integer operations for reduction */
    long long ntg_nbr_nrm=0LL; /* [nbr] Integer operations for normalization */
    long long flp_nbr_bnr=0LL; /* [nbr] Floating point operations for binary arithmetic */
    long long flp_nbr_rdc=0LL; /* [nbr] Floating point operations for reduction */
    long long flp_nbr_nrm=0LL; /* [nbr] Floating point operations for normalization */
    
    // From nco_ddra() in nco_ctl.c
    const float ntg_nbr_brd_fdg_fct=1.8; /* [frc] Empirical correction to broadcasting */
    const float spd_flp_ncwa=153e6; /* [# s-1] Floating point operation speed */
    const float spd_ntg_ncwa=200e6; /* [# s-1] Integer operation speed */
    const float spd_flp_ncbo=353.2e6; /* [# s-1] Floating point operation speed */
    const float spd_ntg_ncbo=1386.54e6; /* [# s-1] Integer operation speed */
    const float spd_rd=63.375e6; /* [B s-1] Disk read bandwidth */
    const float spd_wrt=57.865e6; /* [B s-1] Disk write bandwidth */
    
    float spd_flp=CEWI_flt; /* [# s-1] Floating point operation speed */
    float spd_ntg=CEWI_flt; /* [# s-1] Integer operation speed */
    int rcd=NC_NOERR; /* [rcd] Return code */
    
    // Initialize variables that could be passed into nco_ddra()
    short MRV_flg=1; /* [flg] Avergaging dimensions are MRV dimensions */
    int rnk_avg=CEWI_int; /* [nbr] Rank of averaging space */
    int rnk_var=CEWI_int; /* [nbr] Variable rank (in input file) */
    int rnk_wgt=1; /* [nbr] Rank of weight */
    short wgt_flg=0; /* [flt] Weight name */
    short wgt_brd_flg=0; /* [flg] Broadcast weight for this variable */
    int wrd_sz=4; /* [B] Bytes per element */
    long long lmn_nbr=0LL; /* [nbr] Variable size */
    long long lmn_nbr_avg=0LL; /* [nbr] Averaging block size */
    long long lmn_nbr_wgt=0LL; /* [nbr] Weight size */
    
    // Locals
    long long lmn_nbr_out=CEWI_lng_lng; /* [nbr] Output elements */
    int var_nbr_apx=0;
    int var_idx=0; 
    
    // Choose operation type
    short nco_op_typ=nco_op_typ_avg; /* [enm] Operation type */
    short fl_typ=fl_typ_gcm;
    
    /* Notation:
       flg_acm: Indicates nco_ctl.c code accumulates with += 
       Usually ddra.nco does not because it computes all variables at once
       flg_typ_lng_lng: Indicates type changed to float to hold long longs */
    
    if(fl_typ==fl_typ_gcm){
      var_nbr_apx=32;
      lmn_nbr=var_nbr_apx*varsz_gcm_4D; /* [nbr] Variable size */ /* flg_typ_lng_lng */
      if(nco_op_typ==nco_op_typ_avg){
	lmn_nbr_avg=var_nbr_apx*varsz_gcm_4D; /* [nbr] Averaging block size */ /* flg_typ_lng_lng */
	lmn_nbr_wgt=dmnsz_gcm_lat; /* [nbr] Weight size */
	rnk_var=rnk_avg=1; /* [nbr] Rank of averaging space */
      } // !nco_op_typ_avg
    }else if(fl_typ==fl_typ_stl){
      var_nbr_apx=8;
      lmn_nbr=var_nbr_apx*varsz_stl_2D; /* [nbr] Variable size */ /* flg_typ_lng_lng */
      if(nco_op_typ==nco_op_typ_avg){
	lmn_nbr_avg=var_nbr_apx*varsz_stl_2D; /* [nbr] Averaging block size */ /* flg_typ_lng_lng */
	lmn_nbr_wgt=dmnsz_stl_lat; /* [nbr] Weight size */
	rnk_var=rnk_avg=2; /* [nbr] Rank of averaging space */
      } // !nco_op_typ_avg
    } // !fl_typ
    
    if(nco_op_typ==nco_op_typ_sbt){
      spd_flp=spd_flp_ncbo; /* [# s-1] Floating point operation speed */
      spd_ntg=spd_ntg_ncbo; /* [# s-1] Integer operation speed */
      lmn_nbr_out=lmn_nbr; /* [nbr] Output elements */
    }else if(nco_op_typ==nco_op_typ_avg){
      spd_flp=spd_flp_ncwa; /* [# s-1] Floating point operation speed */
      spd_ntg=spd_ntg_ncwa; /* [# s-1] Integer operation speed */
      lmn_nbr_out=lmn_nbr/lmn_nbr_avg; /* [nbr] Output elements */
    } // !nco_op_typ_avg
    
    flp_nbr_bnr_dfl=lmn_nbr; /* [nbr] Floating point operations for binary arithmetic */
    flp_nbr_nrm_dfl=lmn_nbr_out; /* [nbr] Floating point operations for normalization */
    flp_nbr_rdc_dfl=lmn_nbr; /* [nbr] Floating point operations for reduction */
    
    /* Integer operations for broadcasting weight */
    ntg_nbr_brd_dfl=(long long)(ntg_nbr_brd_fdg_fct*lmn_nbr*(6*rnk_var+8*rnk_wgt+2)); /* [nbr] N(6R+8R_w+2) */
    
    /* Byte-swap integer operations per element */
    ntg_nbr_byt_swp_dfl=wrd_sz+2; /* [nbr nbr-1] W+2 */
    
    /* Integer operations for collection */
    ntg_nbr_clc_dfl=lmn_nbr*(14*rnk_var+4); /* [nbr] N(14R+4) */
    
    /* Integer operations for normalization */
    ntg_nbr_nrm_dfl=4*lmn_nbr_out; /* [nbr] 4N/N_A = 4N_O */
    
    /* Integer operations for reduction */
    ntg_nbr_rdc_dfl=lmn_nbr*6+lmn_nbr_out; /* [nbr] N(6+N/N_A) */
    
    if(nco_op_typ==nco_op_typ_sbt){
      /* Subtraction computation assumes variables are same size
	 fxm: Account for broadcasting */
      /* One floating point (add/subtract/multiply/divide) per element */
      flp_nbr_bnr=flp_nbr_bnr_dfl;
      /* Byte-swap elements from two input files and one output file */
      ntg_nbr_byt_swp=3*lmn_nbr*ntg_nbr_byt_swp_dfl; /* 3N(W+2) */
      rd_nbr_byt=2*lmn_nbr*wrd_sz; /* [B] Bytes read */
      wrt_nbr_byt=lmn_nbr_out*wrd_sz; /* [B] Bytes written */
    }else if(nco_op_typ==nco_op_typ_avg){
      rd_nbr_byt=lmn_nbr*wrd_sz; /* [B] Bytes read */
      wrt_nbr_byt=lmn_nbr_out*wrd_sz; /* [B] Bytes written */
      /* One floating point add per input element to sum numerator */
      flp_nbr_rdc=lmn_nbr;
      /* One floating point divide per output element to normalize numerator by denominatro (tally) */
      flp_nbr_nrm=lmn_nbr_out;
      /* Byte-swap elements from one input file and one (rank-reduced) output file */
      ntg_nbr_byt_swp=(lmn_nbr+lmn_nbr_out)*ntg_nbr_byt_swp_dfl;
      if(!MRV_flg){
	/* Collection required for numerator */
	ntg_nbr_clc=ntg_nbr_clc_dfl; /* flg_acm */
      } /* !MRV_flg */
      if(wgt_flg){
	if(var_idx == 0){
	  /* Set cost = 0 after first variable since only read weight once */
	  rd_nbr_byt=lmn_nbr_wgt*wrd_sz; /* [B] Bytes read */ /* flg_acm */
	  /* Byte-swap cost for first weight input is usually negligible */
	  ntg_nbr_byt_swp=ntg_nbr_byt_swp+lmn_nbr_wgt*ntg_nbr_byt_swp_dfl; /* OK flg_acm */
	} /* var_idx != 0 */
	/* One floating point multiply per input element for weight*value in numerator,
	   and one floating point add per input element to sum weight in denominator */
	flp_nbr_rdc=flp_nbr_rdc+2*lmn_nbr; /* flg_acm */
	/* One floating point divide per output element to normalize denominator by tally */
	flp_nbr_nrm=flp_nbr_nrm+lmn_nbr_out; /* real flg_acm */
	if(wgt_brd_flg){
	  /* fxm: Charge for broadcasting weight at least once */
	  /* Broadcasting cost for weight */
	  ntg_nbr_brd=ntg_nbr_brd_dfl;
	} /* !wgt_brd_flg */
	if(!MRV_flg){
	  /* Collection required for denominator */
	  ntg_nbr_clc=ntg_nbr_clc+ntg_nbr_clc_dfl; /* real acm flg */
	} /* !MRV_flg */
      } /* !wgt_flg */
    } // !nco_op_typ_avg
    
    flp_nbr= /* [nbr] Floating point operations */
      flp_nbr_bnr+ /* [nbr] Floating point operations for binary arithmetic */
      flp_nbr_rdc+ /* [nbr] Floating point operations for reduction */
      flp_nbr_nrm; /* [nbr] Floating point operations for normalization */
    
    ntg_nbr= /* [nbr] Integer operations */
      ntg_nbr_byt_swp+ /* [nbr] Integer operations for byte-swap */
      ntg_nbr_brd+ /* [nbr] Integer operations for broadcasting */
      ntg_nbr_clc+ /* [nbr] Integer operations for collection */
      ntg_nbr_rdc+ /* [nbr] Integer operations for reduction */
      ntg_nbr_nrm; /* [nbr] Integer operations for normalization */
    
    tm_ntg=ntg_nbr/spd_ntg; /* [s] Integer time */
    tm_flp=flp_nbr/spd_flp; /* [s] Floating point time */
    tm_rd=rd_nbr_byt/spd_rd; /* [s] Read time */
    tm_wrt=wrt_nbr_byt/spd_wrt; /* [s] Write time */
    
    tm_io=tm_rd+tm_wrt; /* [s] I/O time */
    tm_crr=tm_ntg+tm_flp+tm_rd+tm_wrt; /* [s] Time for this variable */
    
    tm_ntg_ttl=tm_ntg_ttl+tm_ntg; /* [s] Cumulative integer time */
    tm_flp_ttl=tm_flp_ttl+tm_flp; /* [s] Cumulative floating point time */
    tm_rd_ttl=tm_rd_ttl+tm_rd; /* [s] Cumulative read time */
    tm_wrt_ttl=tm_wrt_ttl+tm_wrt; /* [s] Cumulative write time */
    tm_io_ttl=tm_io_ttl+tm_io; /* [s] Cumulative I/O time */
    tm_ttl=tm_ttl+tm_crr; /* [s] Cumulative time */
    
    tm_frc_flp_ttl=tm_flp_ttl/tm_ttl; /* [frc] Floating point time fraction */
    tm_frc_io_ttl=tm_io_ttl/tm_ttl; /* [frc] I/O time fraction */
    tm_frc_ntg_ttl=tm_ntg_ttl/tm_ttl; /* [frc] Integer time fraction */
    tm_frc_rd_ttl=tm_rd_ttl/tm_ttl; /* [frc] Read time fraction */
    tm_frc_wrt_ttl=tm_wrt_ttl/tm_ttl; /* [frc] Write time fraction */
    
    // Remove unused variable warnings
    dmnszavg_stl+=0.0*(varsz_stl_1D+dmnszavg_gcm_avg_wgt_varnbr+dmnszavg_gcm_avg_wgt_lmnnbr+lmn_nbr_ttl+ntg_nbr_ttl+flp_nbr_ttl+rcd+rnk_avg);

    std::cout << "DDRA results: " << std::endl;
    std::cout << "dmnszavg_stl = " << dmnszavg_stl << std::endl;
    std::cout << "tm_frc_flp_ttl = " << tm_frc_flp_ttl << std::endl;
    std::cout << "tm_frc_io_ttl = " << tm_frc_io_ttl << std::endl;
    std::cout << "tm_frc_ntg_ttl = " << tm_frc_ntg_ttl << std::endl;
    std::cout << "tm_frc_rd_ttl = " << tm_frc_rd_ttl << std::endl;
    std::cout << "tm_frc_wrt_ttl = " << tm_frc_wrt_ttl << std::endl;
    std::cout << "ntg_nbr_rdc_dfl = " << ntg_nbr_rdc_dfl << std::endl;
    std::cout << "ntg_nbr_nrm_dfl = " << ntg_nbr_nrm_dfl << std::endl;
    std::cout << "flp_nbr_rdc_dfl = " << flp_nbr_rdc_dfl << std::endl;
    std::cout << "flp_nbr_nrm_dfl = " << flp_nbr_nrm_dfl << std::endl;

  } // !dbg || tst_sng == "ddra"
  
  if(dbg_lvl == dbg_old || tst_sng == "dyn"){
    std::cout << "Testing dynamic allocation of multidimensional fixed float arrays..." << std::endl;
    std::cout << "Try, e.g.,\nccc --bnd_nbr=3 --sz_nbr=2 --tst=dyn" << std::endl;
    std::cout << "Defined tpt[bnd_nbr][sz_nbr] = tpt_d2d[" << bnd_nbr << "][" << sz_nbr << "]" << std::endl;
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
	//	tpt_d2d[bnd_idx][sz_idx]=bnd_idx*sz_nbr+sz_idx; // [K] Temperature (dynamic 2-dimensional)
	//	std::cout << "tpt_d2d[" << bnd_idx << "][" << sz_idx << "] = " << tpt_d2d[bnd_idx][sz_idx] << std::endl;
	tpt_v1d[bnd_idx*sz_nbr+sz_idx]=bnd_idx*sz_nbr+sz_idx; // [K] Temperature (vector 1-dimensional)
	std::cout << "tpt_v1d[" << bnd_idx*sz_nbr+sz_idx << "] = " << tpt_v1d[bnd_idx*sz_nbr+sz_idx] << std::endl;
      } // !sz_idx
    } // !bnd_idx
  } // !dbg || tst_sng == "dyn"
  
  if(dbg_lvl == dbg_old || tst_sng == "lmt" || tst_sng == "flt"){
    std::cout << "Testing floating point representation using <cfloat>..." << std::endl;
    std::cout << "FLT_DIG = " << FLT_DIG << " = Precision in decimal digits" << std::endl;
    std::cout << "FLT_EPSILON = " << FLT_EPSILON << " = Smallest X such that 1.0 + X != 1.0" << std::endl;
    std::cout << "FLT_MANT_DIG = " << FLT_MANT_DIG << " = Number of mantissa digits, base FLT_RADIX" << std::endl;
    std::cout << "FLT_MAX = " << FLT_MAX << " = Largest finite representable value" << std::endl;
    std::cout << "FLT_MAX_10_EXP = " << FLT_MAX_10_EXP << " = Maximum integer X, such that 10^X is finite representable value" << std::endl;
    std::cout << "FLT_MAX_EXP = " << FLT_MAX_EXP << " = Maximum integer X, such that FLT_RADIX^(X - 1) is finite representable value" << std::endl;
    std::cout << "FLT_MIN = " << FLT_MIN << " = Smallest normalized, finite representable value" << std::endl;
    std::cout << "FLT_MIN_10_EXP = " << FLT_MIN_10_EXP << " = Minimum integer X such that 10^X is normalized, finite representable value" << std::endl;
    std::cout << "FLT_MIN_EXP = " << FLT_MIN_EXP << " = minimum integer X such that FLT_RADIX^(X - 1) is a normalized, finite representable value" << std::endl;
    std::cout << "FLT_RADIX = " << FLT_RADIX << " = Radix of all floating-point representations" << std::endl;
    std::cout << "FLT_ROUNDS = " << FLT_ROUNDS << " = Rounding algorithm" << std::endl;
    std::cout << "DBL_DIG = " << DBL_DIG << " = Precision in decimal digits" << std::endl;
    std::cout << "DBL_EPSILON = " << DBL_EPSILON << " = Smallest X such that 1.0 + X != 1.0" << std::endl;
    std::cout << "DBL_MANT_DIG = " << DBL_MANT_DIG << " = Number of mantissa digits, base FLT_RADIX" << std::endl;
    std::cout << "DBL_MAX = " << DBL_MAX << " = Largest finite representable value" << std::endl;
    std::cout << "DBL_MAX_10_EXP = " << DBL_MAX_10_EXP << " = Maximum integer X, such that 10^X is finite representable value" << std::endl;
    std::cout << "DBL_MAX_EXP = " << DBL_MAX_EXP << " = Maximum integer X, such that FLT_RADIX^(X - 1) is finite representable value" << std::endl;
    std::cout << "DBL_MIN = " << DBL_MIN << " = Smallest normalized, finite representable value" << std::endl;
    std::cout << "DBL_MIN_10_EXP = " << DBL_MIN_10_EXP << " = Minimum integer X such that 10^X is normalized, finite representable value" << std::endl;
    std::cout << "DBL_MIN_EXP = " << DBL_MIN_EXP << " = Minimum integer X such that FLT_RADIX^(X - 1) is normalized, finite representable value" << std::endl;
#ifdef HAVE_LONG_LONG
    // fxm: Comeau has problems here, although it claims long double support
    std::cout << "LDBL_DIG = " << LDBL_DIG << " = Precision in decimal digits" << std::endl;
    std::cout << "LDBL_EPSILON = " << LDBL_EPSILON << " = Smallest X such that 1.0 + X != 1.0" << std::endl;
    std::cout << "LDBL_MANT_DIG = " << LDBL_MANT_DIG << " = Number of mantissa digits, base FLT_RADIX" << std::endl;
    std::cout << "LDBL_MAX = " << LDBL_MAX << " = Largest finite representable value" << std::endl;
    std::cout << "LDBL_MAX_10_EXP = " << LDBL_MAX_10_EXP << " = Maximum integer X, such that 10^X is finite representable value" << std::endl;
    std::cout << "LDBL_MAX_EXP = " << LDBL_MAX_EXP << " = Maximum integer X, such that FLT_RADIX^(X - 1) is finite representable value" << std::endl;
    std::cout << "LDBL_MIN = " << LDBL_MIN << " = Smallest normalized, finite representable value" << std::endl;
    std::cout << "LDBL_MIN_10_EXP = " << LDBL_MIN_10_EXP << " = Minimum integer X such that 10^X is normalized, finite representable value" << std::endl;
    std::cout << "LDBL_MIN_EXP = " << LDBL_MIN_EXP << " = minimum integer X such that FLT_RADIX^(X - 1) is a normalized, finite representable value" << std::endl;
#endif // !HAVE_LONG_LONG
#ifdef HAVE_NUMERIC_LIMITS
    // fxm: Get this working
    // numeric_limits<T> templates can replace all of <cfloat> and <climits>
    // Discussion on Jos99 p. 64
    std::cout << "Testing floating point representation using numeric_limits<T>..." << std::endl;
    std::numeric_limits<float>::max() << std::endl;
#endif // !HAVE_NUMERIC_LIMITS
  } // !dbg || tst_sng == "flt"
  
  if(dbg_lvl == dbg_old || tst_sng == "fnc"){
    // ccc --dbg=3 --function=sine --flt=0.5
    std::cout << "Testing execution time object manipulation..." << std::endl;
    std::cout << "fnc_sng = " << fnc_sng << std::endl;
    std::cout << "cpv_foo = " << cpv_foo << std::endl;
    //  std::cout << "szdstfnc(cpv_foo) = " << szdstfnc(cpv_foo) << std::endl;
    Sine sin_fnc(fnc_sng);
    Cosine cos_fnc(fnc_sng);
    Base base(&sin_fnc);
    if(fnc_sng == "sin"){
      base.set_fnc(&sin_fnc);
    }else if(fnc_sng == "cos"){
      base.set_fnc(&cos_fnc);
    } // endif
    std::cout << "sin_fnc(cpv_foo) = " << sin_fnc(cpv_foo) << std::endl;
    std::cout << "cos_fnc(cpv_foo) = " << cos_fnc(cpv_foo) << std::endl;
    std::cout << "base.fnc->eval(cpv_foo) = " << base.fnc->eval(cpv_foo) << std::endl;
    //  std::cout << "base.eval(cpv_foo) = " << base.eval(cpv_foo) << std::endl;
    //  virtualViaReference(&sin_fnc,cpv_foo);
  } // !dbg || tst_sng == "fnc"
  
  if(dbg_lvl == dbg_old || tst_sng == "gmt"){
    std::cout << "Testing GMT time conversion utilities..." << std::endl;
    const double time_unix_dbl(dbl_foo); // [s] Seconds since 1969
    const long time_unix_lng(static_cast<long>(std::floor(time_unix_dbl))); /* [s] UNIX time stored as intrinsic long */
    const std::time_t time_unix_time_t(time_unix_lng); // [s] Seconds since 1969
    const double tm_sec_frc(time_unix_dbl-time_unix_lng); /* [s] Fractional seconds */
    const struct tm *gmt_tm=gmtime(&time_unix_time_t); // [tm] Time structure
    std::cout << "gmtime info for UNIX time " << time_unix_time_t << " s since 1969:" << std::endl;
    std::cout << "sec = gmtime(" << time_unix_time_t << ")->sec = " << gmt_tm->tm_sec << std::endl;
    std::cout << "mnt = gmtime(" << time_unix_time_t << ")->mnt = " << gmt_tm->tm_min << std::endl;
    std::cout << "hr = gmtime(" << time_unix_time_t << ")->hr = " << gmt_tm->tm_hour << std::endl;
    std::cout << "mdy = gmtime(" << time_unix_time_t << ")->mday = " << gmt_tm->tm_mday << std::endl;
    std::cout << "mth = gmtime(" << time_unix_time_t << ")->mon+1 = " << gmt_tm->tm_mon+1 << std::endl;
    std::cout << "ydy = gmtime(" << time_unix_time_t << ")->yday+1 = " << gmt_tm->tm_yday+1 << std::endl;
    std::cout << "yr = gmtime(" << time_unix_time_t << ")->year+1900 = " << gmt_tm->tm_year+1900 << std::endl;
    int tm_mth(gmt_tm->tm_mon+1); // [mth] 1-based month of year [1..12]
    int tm_day(gmt_tm->tm_mday); // [day] 1-based day of month [1..31]
    int tm_ydy(gmt_tm->tm_yday+1); // [day] 1-based day of year [1..366]
    int tm_yr(gmt_tm->tm_year+1900); // [yr] Christian Year
    int tm_hr(gmt_tm->tm_hour); // [hr] 0-based, hour of day [0..23]
    int tm_mnt(gmt_tm->tm_min); // [mnt] 0-based minute of second [0..59]
    int tm_sec(gmt_tm->tm_sec); // [s] 0-based second of minute [0..59]
    (void)std::fprintf(stdout,"%4.4d/%2.2d/%2.2d YDY %d GMT %2.2d:%2.2d:%f\n",tm_yr,tm_mth,tm_day,tm_ydy,tm_hr,tm_mnt,tm_sec+tm_sec_frc);
  } // !dbg || tst_sng == "gmt"

  if(dbg_lvl == dbg_old || tst_sng == "gsl"){
    std::cout << "Testing GSL behavior..." << std::endl;
    std::cout << "Seed values with --gsl_[a,b,c], e.g., ccc --tst=gsl --gsl_a=1.0 --gsl_b=1.0" << std::endl;
    std::string nfo_msg; // [sng] Descriptive message of context

    /* GSL special function library takes mode argument which can be 
       GSL_PREC_DOUBLE, GSL_PREC_SINGLE, GSL_PREC_APPROX (GDT01 p. 27)
       Set mode in shell before calling executable, e.g.,
       GSL_PREC=GSL_PREC_DOUBLE ccc --tst=gsl --gsl_a=1.0 */
    std::string nvr_GSL_PREC((std::getenv("GSL_PREC")) ? std::getenv("GSL_PREC") : ""); // [sng] Environment variable GSL_PREC
    std::cout << "nvr_GSL_PREC = " << nvr_GSL_PREC << std::endl;
    gsl_mode_t gsl_mode(GSL_PREC_DOUBLE); // [enm] GSL precision mode

    gsl_sf_result nsw_dbl; // [frc] GSL result structure
    bool apx_eql; // [flg] Arguments are indistinguishable
    
    /* Test error handling GDT07 p. 13 */
    gsl_error_handler_t *old_gsl_error_handler; // [fnc] GSL error handler

    // Turn off error handler
    std::cout << "Turning off default GSL error handler..." << std::endl;
    old_gsl_error_handler=gsl_set_error_handler_off(); // [fnc] GSL error handler
    std::cout << "Testing bad numerics with GSL error handler turned off..." << std::endl;

    prc_cmp tau_ext_scl(1000.0); // [frc] Optical depth
    
    std::cout << "Calling gsl_sf_expint_Ei_e(-" << tau_ext_scl << ",&nsw_dbl) to trigger underflow..." << std::endl;
    rcd=gsl_sf_expint_Ei_e(-tau_ext_scl,&nsw_dbl); // [fnc] Exponential integral

    // Save original handler, install new handler
    std::cout << "Installing custom handler..." << std::endl;
    old_gsl_error_handler=gsl_set_error_handler(gsl_error_handler_csz); // [fnc] GSL error handler

    std::cout << "Calling gsl_sf_expint_Ei_e(-" << tau_ext_scl << ",&nsw_dbl) to trigger underflow..." << std::endl;
    rcd=gsl_sf_expint_Ei_e(-tau_ext_scl,&nsw_dbl); // [fnc] Exponential integral
    if(rcd == 15) std::cout << "Presence of WARNING: rcd = 15 in above line indicates custom error handler works." << std::endl; else std::cout << "WARNING: previous call did not produce rcd = 15, custom error handler may not work." << std::endl;

    // Restore original handler
    std::cout << "Restoring default handler..." << std::endl;
    old_gsl_error_handler=gsl_set_error_handler(old_gsl_error_handler); // [fnc] GSL error handler

    if(dbg_lvl == dbg_old){
      // Test non-GSL numeric problems
      std::cout << "Test non-GSL numeric problems..." << std::endl;
      std::cout << "exp(-1000.0) = " << std::exp(-1000.0) << std::endl;
      int cst_one(1); // [nbr] One
      int cst_zero(0); // [nbr] Zero
      std::cout << "One divided by zero (1/0) = " << cst_one/cst_zero << std::endl;
    } // !dbg || tst_sng == "idx_rfr"

    // Complete elliptic integral of the second kind
    std::cout << "Complete elliptic integral of the second kind gsl_sf_ellint_Ecomp(0.0)=pi/2...";
    // ccc --tst=gsl --gsl_a=0.5
    // GDT01 p. 43, AbS64 p. 589
    // http://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html
    rcd=gsl_sf_ellint_Ecomp_e(0.0,gsl_mode,&nsw_dbl); // [fnc] Complete elliptic integral of the second kind
    nfo_msg="Vetting Complete elliptic integral of the second kind gsl_sf_ellint_Ecomp_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       mth::cst_M_PIl/2.0, // I [frc] Target argument: E(pi/2,0.0)=pi/2
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Complete elliptic integral of the second kind gsl_sf_ellint_Ecomp() error");
    if(gsl_a > -1.0 && gsl_a < 1.0){ // ellint_Ecomp is undefined
      rcd=gsl_sf_ellint_Ecomp_e(gsl_a,gsl_mode,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_ellint_Ecomp(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    } // endif
    // !Complete elliptic integral of the second kind test

    // Error function
    std::cout << "Error function gsl_sf_erf(1.0)=0.842701...";
    rcd=gsl_sf_erf_e(1.0,&nsw_dbl); // [fnc] Error function
    nfo_msg="Vetting error function gsl_sf_erf_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       0.842701, // I [frc] Target argument
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-6, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Error function gsl_sf_erf() error");
    rcd=gsl_sf_erf_e(gsl_a,&nsw_dbl);
    if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
    std::cout << "gsl_sf_erf(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    // !Error function test

    // Exponential integral Ei
    // AbS64 p. 233 (5.3) Example 2
    std::cout << "Exponential integral gsl_sf_expint_Ei(8.0)=440.38...";
    rcd=gsl_sf_expint_Ei_e(8.0,&nsw_dbl); // [fnc] Exponential integral
    nfo_msg="Vetting exponential integral gsl_sf_expint_Ei_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       440.38, // I [frc] Target argument
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-6, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Exponential integral gsl_sf_expint_Ei() error");
    if(gsl_a != 0.0){ // Testing Ei(0.0) causes domain error
      rcd=gsl_sf_expint_Ei_e(gsl_a,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_expint_Ei(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    } // endif
    // !Exponential integral Ei test

    // Exponential integral E1
    // NB: E1 and Ei are both called Exponential Integrals, yet are completely different
    std::cout << "Exponential integral gsl_sf_expint_E1(1.0)=0.219384...";
    rcd=gsl_sf_expint_E1_e(1.0,&nsw_dbl); // [fnc] Exponential integral
    nfo_msg="Vetting exponential integral gsl_sf_expint_E1_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       0.219384, // I [frc] Target argument
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-6, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Exponential integral gsl_sf_expint_E1() error");
    if(gsl_a != 0.0){ // Testing E1(0.0) causes domain error
      rcd=gsl_sf_expint_E1_e(gsl_a,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_expint_E1(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    } // endif
    // !Exponential integral E1 test

    // Gamma function
    std::cout << "Gamma function gsl_sf_gamma(0.5)=sqrt(pi)...";
    rcd=gsl_sf_gamma_e(0.5,&nsw_dbl); // [fnc] Gamma function
    using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
    nfo_msg="Vetting gamma function gsl_sf_gamma_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       std::sqrt(mth::cst_M_PIl), // I [frc] Target argument: gamma(0.5)=sqrt(pi)
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Gamma function gsl_sf_gamma() error");
    if(gsl_a != 0.0){ // Testing gamma(0.0) will cause core dump
      rcd=gsl_sf_gamma_e(gsl_a,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_gamma(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    } // !gsl_a
    // !Gamma function test

    // Legendre function
    std::cout << "Legendre function gsl_sf_legendre(2,1.0/sqrt(3))=0.0...";
    // ccc --tst=gsl --gsl_a=2.0 --gsl_b=0.5773
    // GDT01 p. 57
    rcd=gsl_sf_legendre_Pl_e(2,1.0/std::sqrt(3.0),&nsw_dbl); // [fnc] Legendre function
    nfo_msg="Vetting Legendre function gsl_sf_legendre_Pl_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       0.0, // I [frc] Target argument: P2(1.0/sqrt(3))=0.0
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Legendre function gsl_sf_legendre() error");
    int legendre_rdr(static_cast<int>(gsl_a)); // [nbr] Order of Legendre polynomial
    if(legendre_rdr >= 0 && std::fabs(gsl_b) <= 1.0){ // legendre(n < 0, |x| > 1.0) is undefined
      rcd=gsl_sf_legendre_Pl_e(legendre_rdr,gsl_b,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_legendre(" << gsl_a << ", " << gsl_b << ") = " << nsw_dbl.val << std::endl;
    } // endif
    // !Legendre function test

    // Log double factorial function
    std::cout << "Log double factorial function gsl_sf_lndoublefact(10)=log(3840.0)...";
    // http://functions.wolfram.com/06.02.03.0019
    rcd=gsl_sf_lndoublefact_e(10u,&nsw_dbl); // [fnc] Log double factorial ln(x!!)
    nfo_msg="Vetting log double factorial function gsl_sf_lndoublefact_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       std::log(3840.0), // I [frc] Target argument: lndoublefact(10)=log(3840.0)
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR log double factorial function gsl_sf_lndoublefact() error");
    rcd=gsl_sf_lndoublefact_e(static_cast<unsigned int>(gsl_a),&nsw_dbl);
    if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
    std::cout << "gsl_sf_lndoublefact(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    // !Log double factorial function test

    // Incomplete gamma function
    std::cout << "Incomplete gamma function gsl_sf_gamma_inc_Q(1.0,1.0)=e^{-1}...";
    rcd=gsl_sf_gamma_inc_Q_e(1.0,1.0,&nsw_dbl); // [fnc] Normalized incomplete gamma function
    using mth::cst_M_El; // (2.7182818284590452353602874713526625L) [frc] exp(1.0)
    nfo_msg="Vetting incomplete gamma function"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       1.0/cst_M_El, // I [frc] Target argument: gamma_inc_Q(1.0,1.0)=e^{-1}
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-5, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Gamma function gsl_sf_gamma_inc_Q() error");
    if(!(gsl_a == 0.0 && gsl_b == 0.0) && gsl_a >= 0.0){ // Testing gamma_inc(0.0,0.0) will cause core dump
      rcd=gsl_sf_gamma_inc_Q_e(gsl_a,gsl_b,&nsw_dbl);
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      std::cout << "gsl_sf_gamma_inc_Q(" << gsl_a << ", " << gsl_b << ") = " << nsw_dbl.val << std::endl;
    } // endif
    // !Incomplete gamma function test

    // Airy function
    std::cout << "Airy function gsl_sf_airy_Ai_e(0)=1.0/(pow(3,2/3)*gamma(2/3))...";
    // http://functions.wolfram.com/BesselAiryStruve/AiryAi/
    rcd=gsl_sf_airy_Ai_e(0.0,gsl_mode,&nsw_dbl); // [fnc] Airy function
    nfo_msg="Vetting Airy function gsl_sf_airy_Ai_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       1.0/(std::pow(3.0,2.0/3.0)*gsl_sf_gamma(2.0/3.0)), // I [frc] Target argument: Ai(0)=1.0/(pow(3,2/3)*gamma(2/3))
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Airy function gsl_sf_airy_Ai() error");
    rcd=gsl_sf_airy_Ai_e(gsl_a,gsl_mode,&nsw_dbl);
    if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
    std::cout << "gsl_sf_airy_Ai(" << gsl_a << ") = " << nsw_dbl.val << std::endl;
    // !Airy function test

    // Regular cylindrical Bessel function of fractional order
    std::cout << "Bessel function gsl_sf_bessel_Jnu_e(0)...";
    rcd=gsl_sf_bessel_Jnu_e(0.5,0.0,&nsw_dbl); // [fnc] Regular cylindrical Bessel function of fractional order
    nfo_msg="Vetting regular cylindrical Bessel function of fractional order gsl_sf_bessel_Jnu_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       0.0, // I [frc] Target argument: J_{1/2}(0)=
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Bessel function gsl_sf_bessel_Jnu() error");
    if(gsl_a >= 0.0){ // Testing bessel_Jnu(0.5,0.0) will cause core dump
      rcd=gsl_sf_bessel_Jnu_e(0.5,gsl_a,&nsw_dbl);
    } // endif
    if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
    std::cout << "gsl_sf_bessel_Jnu(0.5," << gsl_a << ") = " << nsw_dbl.val << std::endl;
    // !Regular cylindrical Bessel function of fractional order test

    // Zeros of regular Bessel function of order 1
    // ccc --tst=gsl --gsl_uint=5
    std::cout << "Zeros of regular Bessel function of order 1 gsl_sf_bessel_zero_J1_e(5)...";
    /* http://mathworld.wolfram.com/BesselFunctionZeros.html has j_1..j_5 of J_1
       Higher precision and zeros for derivatives of Bessel functions at
       http://webcomputing.bio.bas.bg/webMathematica/webComputing/BesselZeros.jsp */
    rcd=gsl_sf_bessel_zero_J1_e(5,&nsw_dbl); // [fnc] Zeros of regular Bessel function of order 1
    nfo_msg="Vetting zeros of Bessel function of order 1 gsl_sf_bessel_zero_J1_e()"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       16.47063005087763, // I [frc] Target argument: J_{1}(x_5)=0.0
       nsw_dbl.val, // I [frc] Approximation to target argument
       1.0e-7, // I [frc] Relative precision
       nfo_msg); // I [sng] Descriptive message of context
    if(apx_eql) std::cout << "OK" << std::endl; else err_prn(prg_nm,sbr_nm,"ERROR Bessel function gsl_sf_bessel_zero_J1() error");
    if(gsl_uint > 0){ // Testing bessel_zero_J1(0) will cause core dump
      rcd=gsl_sf_bessel_zero_J1_e(gsl_uint,&nsw_dbl);
    } // endif
    if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
    std::cout << "gsl_sf_bessel_zero_J1_e(" << gsl_uint << ") = " << nsw_dbl.val << std::endl;
    // !Zeros of regular Bessel function of order 1 test

    /* Binary representation of floating point variables
       Excellent discussion in GDT01 p. 364 */
    float gsl_a_flt=static_cast<float>(gsl_a);
    double gsl_a_flt_dbl=static_cast<double>(gsl_a_flt);
    std::cout << "         gsl_a = " << gsl_a << std::endl;
    std::cout << "   float gsl_a = "; gsl_ieee_fprintf_float(stdout,&gsl_a_flt);std::cout << std::endl;
    std::cout << "promoted gsl_a = "; gsl_ieee_fprintf_double(stdout,&gsl_a_flt_dbl);std::cout << std::endl;
    std::cout << "  double gsl_a = "; gsl_ieee_fprintf_double(stdout,&gsl_a);std::cout << std::endl;
    
  } // !dbg || tst_sng == "gsl"
  
  if(dbg_lvl == dbg_old || tst_sng == "idx_rfr"){
    std::cout << "Testing idx_rfr_cls..." << std::endl;
    std::cout << "ccc --tst=idx_rfr --lng_foo=1 --sng_foo='saharan_dust'" << std::endl;
    idx_rfr_cls idx_rfr(sng_foo);
    std::cout << idx_rfr << std::endl;
    idx_rfr.tst(sng_foo,lng_foo);
    idx_rfr_H2O_tst();
  } // !dbg || tst_sng == "idx_rfr"

  if(dbg_lvl == dbg_old || tst_sng == "ld"){
#ifdef AIX
    std::cout << "Testing linker data..." << std::endl;
    /* External symbols defined by linker, see `man ld`
       These symbols are used by taking their addresses 
       Do NOT redefine them */
    extern char _text[]; // [] First location of program
    extern char _etext[]; // [] First location after program
    extern char _data[]; // [] First location of data
    extern char _edata[]; // [] First location after initialized data
    extern char _end[]; // [] First location after all data
#else // !AIX
    std::cout << "Linker test only available with AIX..." << std::endl;
#endif // !AIX
  } // !dbg || tst_sng == "ld"

  if(dbg_lvl == dbg_old || tst_sng == "lmt" || tst_sng == "ntg" ){
    std::cout << "Testing integer representation using <climits>..." << std::endl;
    std::cout << "SCHAR_MIN = " << SCHAR_MIN << " = Minimum value a `signed char' can hold" << std::endl;
    std::cout << "SCHAR_MAX = " << SCHAR_MAX << " = Maximum value a `signed char' can hold" << std::endl;
    std::cout << "UCHAR_MAX = " << UCHAR_MAX << " = Maximum value an `unsigned char' can hold" << std::endl;
    std::cout << "SHRT_MIN = " << SHRT_MIN << " = Minimum value a `signed short int' can hold" << std::endl;
    std::cout << "SHRT_MAX = " << SHRT_MAX << " = Maximum value a `signed short int' can hold" << std::endl;
    std::cout << "USHRT_MAX = " << USHRT_MAX << " = Maximum value an `unsigned short int' can hold" << std::endl;
    std::cout << "INT_MIN = " << INT_MIN << " = Minimum value a `signed int' can hold" << std::endl;
    std::cout << "INT_MAX = " << INT_MAX << " = Maximum value a `signed int' can hold" << std::endl;
    std::cout << "UINT_MAX = " << UINT_MAX << " = Maximum value an `unsigned int' can hold" << std::endl;
    std::cout << "LONG_MIN = " << LONG_MIN << " = Minimum value a `signed long int' can hold" << std::endl;
    std::cout << "LONG_MAX = " << LONG_MAX << " = Maximum value a `signed long int' can hold" << std::endl;
    std::cout << "ULONG_MAX = " << ULONG_MAX << " = Maximum value an `unsigned long int' can hold" << std::endl;
    // 20070518: g++ 4.0.3 on sand seems to have LLONG_MIN but not LLONG_MAX nor ULLONG_MAX
    // 20070530: g++ 4.1.2 finds LLONG_MAX, ULLONG_MAX iff #include stdint.h 
#if (defined LLONG_MIN)
    std::cout << "LLONG_MIN = " << LLONG_MIN << " = Minimum value a `signed long long int' can hold" << std::endl;
#endif // !LLONG_MIN
#if (defined LLONG_MAX)
    std::cout << "LLONG_MAX = " << LLONG_MAX << " = Maximum value a `signed long long int' can hold" << std::endl;
#endif // !LLONG_MAX
#if (defined ULLONG_MAX)
    std::cout << "ULLONG_MAX = " << ULLONG_MAX << " = Maximum value an `unsigned long long int' can hold" << std::endl;
#endif // !ULLONG_MAX
  } // !dbg || tst_sng == "lmt" || tst_sng == "ntg"

  if(dbg_lvl == dbg_old || tst_sng == "map"){
    std::cout << "Testing STL map containers" << std::endl;
    using phc::mmw_C; // (12.011e-03) [kg mol-1] Mean molecular weight of C IUPAC
    using phc::mmw_SiO2; // [kg mol-1] Mean molecular weight of SiO2
    std::string cmp_sng_prt="biomass"; // [sng] Aerosol abbreviation
    prc_cmp dns_prt; // [kg m-3] Density of particle
    // Test simple one-to-one maps
    typedef std::map<std::string,prc_cmp,std::less<std::string> > sng2cpv_map; // String-to-computational precision map
    sng2cpv_map aer_abb2dns; // Key is string abbreviation, value is aerosol dry density
    aer_abb2dns.insert(sng2cpv_map::value_type("biomass",1000.0)); // kg m-3 Biomass ???
    aer_abb2dns.insert(sng2cpv_map::value_type("dust_like",2500.0)); // kg m-3 Dust-like BPB ARESE report 2/96, Table 1, PaG77 p. 2076, (DKS91 p. 118 has 1600, Volz 1972)
    dns_prt=aer_abb2dns[cmp_sng_prt]; // [kg m-3] Density of particle
    std::cerr << "dns_prt = " << dns_prt << " kg m-3" << std::endl;
    std::cerr << "aer_abb2dns[\"dust_like\"] = " << aer_abb2dns["dust_like"] << std::endl;
    // Test more complicated one-to-many maps
    const mnr_tst_sct mnr[]={
      {"biomass", // [sng] Particle abbreviation
       "C", // [sng] Molecular composition
       (prc_cmp)mmw_C, // [kg mol-1] Mean molecular weight
       1000.0}, // [kg m-3] Bulk density fxm: get Biomass density
      {"dust_like", // [sng] Particle abbreviation
       "SiO2", // [sng] Molecular composition
       (prc_cmp)mmw_SiO2, // [kg mol-1] Mean molecular weight
       2500.0} // [kg m-3] Bulk density BPB ARESE report 2/96, Table 1, PaG77 p. 2076, (DKS91 p. 118 has 1600, Volz 1972)
      /*
      {"", // [sng] Particle abbreviation
       "", // [sng] Molecular composition
       , // [kg mol-1] Mean molecular weight
       }, // [kg m-3] Bulk density fxm: get Biomass density
      */
    }; // !mnr_tst_sct mnr[]
    int mnr_nbr=sizeof(mnr)/sizeof(mnr_tst_sct);
    sng2mnr_tst_map mnr_map; // [sct] Map with key=abbreviation, value=mineral structure
    for(idx=0;idx<mnr_nbr;idx++){
      /* fxm: Define variables before inserting into map, because map values 
	 seem to be unwritable (read-only) once they are in map. */
      mnr_map.insert(sng2mnr_tst_map::value_type(mnr[idx].abb,mnr[idx])); // [sct] Map with key=abbreviation, value=mineral structure
    } // !idx
    // Inserting and instantiating a structure at the same time does not work
    /*    mnr_map.insert(sng2cpv_map::value_type("afghan_dust",
    {"afghan_dust", // [sng] Particle abbreviation
       "SiO2", // [sng] Molecular composition
       (prc_cmp)mmw_SiO2, // [kg mol-1] Mean molecular weight
       2500.0})); // [kg m-3] Bulk density SAJ93 (same density as Dust-like, for convenience, for now) */
    dns_prt=mnr_map.find(cmp_sng_prt)->second.dns; // [kg m-3] Bulk density
    prc_cmp mmw_aer=mnr_map.find(cmp_sng_prt)->second.mmw; // [kg mol-1] Mean molecular weight
    std::cerr << "dns_prt = " << dns_prt << " kg m-3" << std::endl;
    std::cerr << "mmw_aer = " << mmw_aer << " kg mol-1" << std::endl;
    // fxm: indexing the map by the key (like hashing in Perl) does not work
    //    std::cerr << "mnr_map[\"dust_like\"] = " << mnr_map["dust_like"] << std::endl;
  } // !dbg || tst_sng == "map"
  
  if(dbg_lvl == dbg_old || tst_sng == "mmr"){
    std::cout << "Testing memory usage statistics gleaned by getrusage()..." << std::endl;
    int sz_pg(getpagesize()); /* [B] Page size in Bytes */
    std::cout << "getpagesize() reports page size = " << sz_pg << "B" << std::endl;
#ifdef AIX
    std::cout << "System type is AIX so rusage uses kilobytes [kB] for size and seconds [s] for time" << std::endl;
#endif // !AIX
#ifdef LINUX
    std::cout << "System type is LINUX so, as of kernel 2.6.9, rusage only maintains the fields ru_utime, ru_stime, ru_minflt, ru_majflt, and ru_nswap. rusage does not yet implement ru_maxrss, ru_ixrss, ru_idrss, and ru_idrss." << std::endl;
#endif /* !LINUX */
#ifdef SUNMP
    std::cout << "System type is SUNMP so rusage uses pages [pg] for size and ticks [tck] for time." << std::endl;
#endif // !SUNMP

    struct rusage usg;
    rcd+=getrusage(RUSAGE_SELF,&usg);
    std::cout << "rusage.ru_utime.tv_sec = user time used, seconds = " << usg.ru_utime.tv_sec << std::endl;
    std::cout << "rusage.ru_utime.tv_usec = user time used, microseconds = " << usg.ru_utime.tv_usec << std::endl;
    std::cout << "rusage.ru_stime.tv_sec = system time used, seconds = " << usg.ru_stime.tv_sec << std::endl;
    std::cout << "rusage.ru_stime.tv_usec = system time used, microseconds = " << usg.ru_stime.tv_usec << std::endl;
    std::cout << "rusage.ru_maxrss = maximum resident set size = " << usg.ru_maxrss << std::endl;
    std::cout << "rusage.ru_ixrss = integral shared memory size = " << usg.ru_ixrss << std::endl;
    std::cout << "rusage.ru_idrss = integral unshared data size = " << usg.ru_idrss << std::endl;
    std::cout << "rusage.ru_isrss = integral unshared stack size = " << usg.ru_isrss << std::endl;
    std::cout << "rusage.ru_minflt = page reclaims = " << usg.ru_minflt << std::endl;
    std::cout << "rusage.ru_majflt = page faults = " << usg.ru_majflt << std::endl;
    std::cout << "rusage.ru_nswap = swaps = " << usg.ru_nswap << std::endl;

    std::cout << "Testing std::malloc() routines..." << std::endl;
    std::cout << "HINT: run with MALLOC_CHECK_=1 ccc --tst=mmr --lng=4" << std::endl;
    std::cout << prg_nm_get() << ": INFO std::malloc()'ing " << lng_foo << " bytes..." << std::endl;
    char *chr_ptr=(char *)std::malloc(lng_foo); // [sng] Pointer to buffer
    //    char *chr_usr_ptr=(char *)std::malloc_usable_space(lng_foo); // [sng] Memory allocated for user manipulation Linux Journal #87 (July, 2001, p.~82)
    std::cout << prg_nm_get() << ": INFO std::strcpy()'ing to overflow allocated memory..." << std::endl;
    (void)std::strcpy(chr_ptr,"This string overflows allocated memory");
    std::cout << prg_nm_get() << ": INFO std::fprintf()'ing string with overflow..." << std::endl;
    (void)std::fprintf(stdout,"%s\n",chr_ptr);
    std::cout << prg_nm_get() << ": INFO std::free()'ing allocated memory..." << std::endl;
    std::free(chr_ptr); // [sng] Pointer to buffer
    std::cout << prg_nm_get() << ": INFO std::strcpy()'ing to std::free()'d memory..." << std::endl;
    (void)std::strcpy(chr_ptr,"This string written to already std::free()'d memory");
    std::cout << prg_nm_get() << ": INFO std::fprintf()'ing string written to already std::free()'d memory..." << std::endl;
    (void)std::fprintf(stdout,"%s\n",chr_ptr);
    std::cout << prg_nm_get() << ": INFO std::free()'ing already std::free()'d memory..." << std::endl;
    std::free(chr_ptr); // [sng] Pointer to buffer
    std::cout << prg_nm_get() << ": INFO std::realloc()'ing already std::free()'d memory..." << std::endl;
    chr_ptr=(char *)std::realloc(chr_ptr,lng_foo); // [sng] Pointer to buffer
    std::cout << prg_nm_get() << ": INFO std::strcpy()'ing to erroneous memory..." << std::endl;
    (void)std::strcpy(chr_ptr,"This string std::strcpy()'d to erroneous memory");
    std::cout << prg_nm_get() << ": INFO std::fprintf()'ing string std::strcpy()'d to erroneous memory..." << std::endl;
    (void)std::fprintf(stdout,"%s\n",chr_ptr);
  } // !dbg || tst_sng == "mmr"

  if(dbg_lvl == dbg_old || tst_sng == "mpi"){
    /* mpirun -np 2 ccc --dbg=0 --tst=mpi
       mpiexec -n 2 ccc --dbg=0 --tst=mpi 
       mpiexec -gdb -n 2 ccc --dbg=0 --tst=mpi */
    /* Message Passing Interface (MPI) parallelization for clusters
       NB: Modern MPI implementations have an OO C++ interface
       Not sure if/how standard C++ interface is
       Using C interface for now
       Hello World example at
       http://beige.ucs.indiana.edu/B673/node100.html */
    std::cout << "Testing MPI implementation..." << std::endl;
#ifdef ENABLE_MPI
    std::cout << "MPI Activation token ENABLE_MPI is true" << std::endl; 
#ifdef HOST_NAME_MAX
#define NCO_HOST_NAME_MAX HOST_NAME_MAX
#else // !HOST_NAME_MAX
#define NCO_HOST_NAME_MAX 256
#endif // !HOST_NAME_MAX
    char mpi_hst_nm[NCO_HOST_NAME_MAX];
    int mpi_hst_nm_lng;
    int mpi_prc_id; /* [id] Process ID */
    int mpi_prc_nbr(0); /* [nbr] MPI process number */
    // MPI Initialization
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_prc_nbr);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_prc_id);
    MPI_Get_processor_name(mpi_hst_nm,&mpi_hst_nm_lng);
    if(mpi_prc_id == 0){ /* MPI manager code */
      std::cout << "MPI: MPI process number mpi_prc_nbr = " << mpi_prc_nbr << std::endl;
    } // mpi_prc_id != 0
    std::cout << "MPI: MPI process ID mpi_prc_id = " << mpi_prc_id << " executing on host " << mpi_hst_nm << std::endl;
    MPI_Finalize();
    const int mpi_prc_nbr_max_fsh(4); // [nbr] Maximum number of processes program can use efficiently
    std::valarray<int> prc_idx; // [idx] Array to hold process indices
    if(mpi_prc_nbr > mpi_prc_nbr_max_fsh){
      std::cout << prg_nm_get() << ": INFO Reducing number of processes from " << mpi_prc_nbr << " to " << mpi_prc_nbr_max_fsh << " since " << prg_nm_get() << " hits I/O bottleneck above " << mpi_prc_nbr_max_fsh << " processes" << std::endl;
    } // endif
#else // !ENABLE_MPI
    std::cout << "MPI Activation token ENABLE_MPI is false" << std::endl; 
    std::cout << prg_nm_get() << ": INFO Not attempting MPI processes" << std::endl;
#endif // !ENABLE_MPI
  } // !dbg || tst_sng == "mpi"

  if(dbg_lvl == dbg_old || tst_sng == "mth"){
    std::cout << "Testing mathematical utilities, constants..." << std::endl;
    std::cout << "ccc --tst_sng=mth --prm_a=1.0 --prm_b=0.0 --prm_c=-1.0" << std::endl;
    std::cout << "sgn(PRC_CMP(-1.0)) = " << sgn(PRC_CMP(-1.0)) << std::endl;
    std::cout << "sgn(PRC_CMP(-1.0)) = " << sgn(PRC_CMP(-1.0)) << std::endl;
    using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
    using mth::cst_M_LN2l; // (0.6931471805599453094172321214581766L) [frc] log_e(2.0)
    using mth::cst_M_El; // (2.7182818284590452353602874713526625L) [frc] exp(1.0)
    using mth::cst_M_SQRT2l; // (1.4142135623730950488016887242096981L) [frc] sqrt(2.0)
    std::cout << "mth::cst_M_PIl = " << mth::cst_M_PIl << std::endl;
    std::cout << "cst_M_LN2l = " << cst_M_LN2l << std::endl;
    std::cout << "cst_M_El = " << cst_M_El << std::endl;
    std::cout << "cst_M_SQRT2l = " << cst_M_SQRT2l << std::endl;
    std::complex<prc_cmp> sln_1; // [frc] First root of complex equation
    std::complex<prc_cmp> sln_2; // [frc] Second root of complex equation
    std::complex<prc_cmp> sln_3; // [frc] Third root of complex equation
    // Solve quadratic equation with real coefficients and complex solutions
    rcd+=eqn_qdr_slvr(prm_a.real(),prm_b.real(),prm_c.real(),&sln_1,&sln_2);
    std::cout << "Quadratic with real coefficients a = " << prm_a.real() << ", b = " << prm_b.real() << ", c = " << prm_c.real() << " has solutions: sln_1 = " << sln_1 << ", sln_2 = " << sln_2 << std::endl;
    // Solve quadratic equation with complex coefficients and complex solutions
    rcd+=eqn_qdr_slvr(prm_a,prm_b,prm_c,&sln_1,&sln_2);
    std::cout << "Quadratic with complex coefficients a = " << prm_a << ", b = " << prm_b << ", c = " << prm_c << " has solutions: sln_1 = " << sln_1 << ", sln_2 = " << sln_2 << std::endl;
    // Solve cubic equation with real coefficients and complex solutions
    rcd+=eqn_cbc_slvr(prm_a.real(),prm_b.real(),prm_c.real(),&sln_1,&sln_2,&sln_3);
    std::cout << "Cubic with real coefficients a = " << prm_a.real() << ", b = " << prm_b.real() << ", c = " << prm_c.real() << " has solutions sln_1 = " << sln_1 << ", sln_2 = " << sln_2 << ", sln_3 = " << sln_3 << std::endl;
    // Solve cubic equation with complex coefficients and complex solutions
    rcd+=eqn_cbc_slvr(prm_a,prm_b,prm_c,&sln_1,&sln_2,&sln_3);
    std::cout << "Cubic with complex coefficients a = " << prm_a << ", b = " << prm_b << ", c = " << prm_c << " has solutions: sln_1 = " << sln_1 << ", sln_2 = " << sln_2 << ", sln_3 = " << sln_3 << std::endl;
  } // !dbg || tst_sng == "mth"

  if(dbg_lvl == dbg_old || tst_sng == "nms"){
    // ccc --dbg=3 --tst=nms
    std::cout << "Testing namespaces..." << std::endl;
    std::cout << "nmspc1::pi = " << nmspc1::pi << std::endl;
    std::cout << "nmspc2::pi = " << nmspc2::pi << std::endl;
  } // !dbg || tst_sng == "nms"
  
  if(dbg_lvl == dbg_old || tst_sng == "ntp" || tst_sng == "rbn" || tst_sng == "vec"){
    // ccc --dbg=3 --xtr_LHS="xtr_prt_wgt+xtr_fll_nil+xtr_vrb" --xtr_RHS="xtr_prt_wgt+xtr_fll_nil+xtr_vrb"
    // ccc --dbg=3 --xtr_LHS="xtr_prt_wgt+xtr_fll_lnr+xtr_vrb" --xtr_RHS="xtr_prt_wgt+xtr_fll_lnr+xtr_vrb"
    // ccc --dbg=3 --xtr_LHS="xtr_prt_wgt+xtr_fll_ngh+xtr_vrb" --xtr_RHS="xtr_prt_wgt+xtr_fll_ngh+xtr_vrb"
    std::cout << "Testing ntp_vec() and rbn_vec()..." << std::endl;
    vec_tst(xtr_LHS,xtr_RHS);
  } // !dbg || tst_sng == "ntp", "rbn", "vec"
  
  if(dbg_lvl == dbg_old || tst_sng == "nvr"){
    std::cout << "Testing environment variables..." << std::endl;
    std::string nvr_TZ((std::getenv("TZ")) ? std::getenv("TZ") : ""); // [sng] Environment variable TZ
    std::string nvr_DATA_tst((std::getenv("DATA")) ? std::getenv("DATA") : ""); // [sng] Environment variable DATA
    std::string nvr_HOME_tst((std::getenv("HOME")) ? std::getenv("HOME") : ""); // [sng] Environment variable HOME
    std::string nvr_USER_tst((std::getenv("USER")) ? std::getenv("USER") : ""); // [sng] Environment variable USER
    std::string nvr_BUFSIZ((std::getenv("BUFSIZ")) ? std::getenv("BUFSIZ") : ""); // [sng] Environment variable BUFSIZ
    std::string nvr_QUARK((std::getenv("QUARK")) ? std::getenv("QUARK") : ""); // [sng] Environment variable QUARK
    std::cout << "nvr_TZ = " << nvr_TZ << ", nvr_HOME_tst = " << nvr_HOME_tst << ", nvr_USER_tst = " << nvr_USER_tst << ", nvr_DATA_tst = " << nvr_DATA_tst << ", nvr_BUFSIZ = " << nvr_BUFSIZ << ", nvr_QUARK = " << nvr_QUARK << std::endl;

    std::cout << "Complete list of environment variables..." << std::endl;
    extern char **environ; // [sng] System-owned pointer to list of pointers to environment variables
    char **nvr=environ; // [sng] Copy of system-owned pointer to list of pointers to environment variables
    char *nvr_crr; // [sng] System-owned pointer to current environment variable
    // Step through list of pointers until NULL pointer termination is reached
    // Inner parentheses allow nvr_crr to be updated, outer parens for boolean evaluation of loop
    // Technique from UNIX programming FAQ
    while((nvr_crr=*nvr++)) std::cout << nvr_crr << std::endl;
  } // !dbg || tst_sng == "nvr"

  if(dbg_lvl == dbg_old || tst_sng == "nwt_rph"){
    std::cout << "Testing Newton-Raphson root finder..." << std::endl;
    /* prm_a is bnd_LHS
       prm_b is bnd_RHS
       prm_c is cnv_thr
       Usage:
       ccc --dbg=5 --tst_sng=nwt_rph --eqn_sng=tst --prm_a=-10.0 --prm_b=20.0 --prm_c=0.0001
       ccc --dbg=5 --tst_sng=nwt_rph --eqn_sng=lnstr_rat --prm_a=1.0 --prm_b=1.0e6 --prm_c=0.0001
       ccc --dbg=5 --tst_sng=nwt_rph --eqn_sng=wien --prm_a=1.0 --prm_b=10.0 --prm_c=0.0001
       ccc --dbg=5 --tst_sng=nwt_rph --eqn_sng=sci_prg --prm_a=0.0 --prm_b=10.0 --prm_c=1.0e-5
    */
    prc_cmp sln; // [frc] Solution to equation
    rcd+=nwt_rph_slvr(eqn_sng,prm_a.real(),prm_b.real(),prm_c.real(),&sln);
    std::cout << "Solving equation " << eqn_sng << std::endl;
    std::cout << "Solution is x = " << sln << std::endl;
  } // !dbg || tst_sng == "nwt_rph"

  if(dbg_lvl == dbg_old || tst_sng == "ocn"){
    std::cout << "Testing exchange coefficient over ocean..." << std::endl;
    using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant
    long wnd_nbr=100;
    long wnd_idx;
    prc_cmp wnd_10m_ntr;
    prc_cmp rgh_mmn_ocn;
    prc_cmp xch_cff_mmn_ocn_ntr;
    for(wnd_idx=0;wnd_idx<wnd_nbr;wnd_idx++){
      wnd_10m_ntr=0.1+50.0*wnd_idx/static_cast<prc_cmp>(wnd_nbr); // [m s-1]
      xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_ntr); // [frc]
      rgh_mmn_ocn=10.0*std::exp(-cst_von_krm/std::sqrt(xch_cff_mmn_ocn_ntr)); // [m]
      std::cout << "idx " << wnd_idx << ", wnd_10m_ntr = " << wnd_10m_ntr << " m s-1, xch_cff_mmn_ocn_ntr = " << xch_cff_mmn_ocn_ntr << ", rgh_mmn_ocn = " << rgh_mmn_ocn << " m" << std::endl;
    } // end loop over wnd
    /* IDL code:
       sz_nbr=100
       cst_von_krm=0.4
       wnd_10m_ntr=0.1+50.0*findgen(sz_nbr)/sz_nbr
       xch_cff_mmn_ocn_ntr=0.0027/wnd_10m_ntr+0.000142+0.0000764*wnd_10m_ntr
       rgh_mmn_ocn=10.0*std::exp(-cst_von_krm/std::sqrt(xch_cff_mmn_ocn_ntr))
       plot,wnd_10m_ntr,xch_cff_mmn_ocn_ntr
    */
  } // !dbg || tst_sng == "ocn"
  
  if(dbg_lvl == dbg_old || tst_sng == "cxx11"){
    // ccc --dbg=1 --tst=cxx11
      /* C++-11 threads
	 http://baptiste-wicht.com/posts/2012/03/cpp11-concurrency-part1-start-threads.html
      */
    std::cout << "Testing C++-11 features..." << std::endl; 
    std::vector<std::thread> threads;
    for(int i=0;i<5;++i) threads.push_back(std::thread(cxx11_hello));
    for(auto& thread: threads) thread.join();
  } // !dbg || tst_sng == "cxx11"

  if(dbg_lvl == dbg_old || tst_sng == "omp"){
    // ccc --dbg=1 --tst=omp
      /* OpenMP parallelization for SMP systems
       OpenMP API:
       http://www.openmp.org/specs
       file:/usr/local/pgi/doc/pgiws_ug/pgi31u.htm
       Sample OpenMP programs:
       Fortran: babyblue.ucar.edu:/usr/local/examples/openmp
       C: utefe.ucar.edu:~rosinski/timing */
    std::cout << "Testing OpenMP implementation..." << std::endl; 
#ifdef _OPENMP // OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification
    std::cout << "OpenMP Activation token _OPENMP is true" << std::endl; 
    int dyn_thr(1); // [flg] Allow system to dynamically set number of threads
    int ntg_OMP_NUM_THREADS=CEWI_int; // [nbr] OMP_NUM_THREADS environment variable
    const int thr_nbr_max_fsh(4); // [nbr] Maximum number of threads program can use efficiently
    const int thr_nbr_max(omp_get_max_threads()); // [nbr] Maximum number of threads system allows
    const int prc_nbr_max(omp_get_num_procs()); // [nbr] Number of processors available
    const std::string nvr_OMP_NUM_THREADS((std::getenv("OMP_NUM_THREADS")) ? std::getenv("OMP_NUM_THREADS") : ""); // [sng] Environment variable OMP_NUM_THREADS

    std::valarray<int> thr_idx; // [idx] Array to hold thread indices
    if(dbg_lvl > 0){
      std::cout << prg_nm_get() << ": INFO Value of _OPENMP token is " << _OPENMP << std::endl;
      if(nvr_OMP_NUM_THREADS != "") ntg_OMP_NUM_THREADS=static_cast<int>(std::strtol(nvr_OMP_NUM_THREADS.c_str(),(char **)NULL,10));
      std::cout << prg_nm_get() << ": INFO Environment variable OMP_NUM_THREADS ";
      if(ntg_OMP_NUM_THREADS > 0) std::cout << "= " << ntg_OMP_NUM_THREADS; else std::cout << "does not exist";
      std::cout << std::endl;
      std::cout << prg_nm_get() << ": INFO Maximum number of threads system allows = omp_get_max_threads() = " << thr_nbr_max << std::endl;
      std::cout << prg_nm_get() << ": INFO Number of processors available = omp_get_num_procs() = " << prc_nbr_max << std::endl;
    } // !dbg_lvl
    if(thr_nbr_max > thr_nbr_max_fsh){
      std::cout << prg_nm_get() << ": INFO Reducing number of threads from " << thr_nbr_max << " to " << thr_nbr_max_fsh << " since " << prg_nm_get() << " hits I/O bottleneck above " << thr_nbr_max_fsh << " threads" << std::endl;
      (void)omp_set_num_threads(thr_nbr_max_fsh); // [nbr] Maximum number of threads system is allowed
    } // endif
    (void)omp_set_dynamic(dyn_thr); // [flg] Allow system to dynamically set number of threads
    if(dbg_lvl > 0) std::cout << prg_nm_get() << ": INFO " << (dyn_thr ? "Allowing" : "Not Allowing") << " OS to utilize dynamic threading" << std::endl;
    dyn_thr=omp_get_dynamic(); // [flg] Allow system to dynamically set number of threads
    if(dbg_lvl > 0) std::cout << prg_nm_get() << ": INFO System will " << (dyn_thr ? "" : "not ") << "utilize dynamic threading" << std::endl;
#pragma omp parallel
    { // begin OpenMP parallel
#pragma omp single nowait
      { // begin OpenMP single
	std::cout << prg_nm_get() << ": INFO OpenMP threading with omp_get_num_threads() = " << omp_get_num_threads() << " threads" << std::endl;
      } // end OpenMP single
    } // end OpenMP parallel
#else // !_OPENMP
    std::cout << "OpenMP Activation token _OPENMP is false" << std::endl; 
    std::cout << prg_nm_get() << ": INFO Not attempting OpenMP threading" << std::endl;
#endif // !_OPENMP
#ifdef _OPENMP // OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification
#pragma omp parallel default(none) shared(thr_idx,thr_nbr)
    { // begin omp parallel
#pragma omp single
      { // begin OpenMP single
	thr_nbr=omp_get_num_threads(); // [nbr] Number of threads
	thr_idx.resize(static_cast<size_t>(thr_nbr),-1); // [idx] Array to hold thread indices
      } // end OpenMP single
      thr_idx[static_cast<size_t>(omp_get_thread_num())]=omp_get_thread_num();
    } // end omp parallel
    std::cout << "Values != -1 prove that parallel loop was executed:" << std::endl;
    for(idx=0;idx<thr_nbr;idx++){
      std::cout << "thr_idx[" << idx << "]=" << thr_idx[idx] << std::endl;
    } // end loop over idx
#endif // endif _OPENMP
  } // !dbg || tst_sng == "omp"

  if(dbg_lvl == dbg_old || tst_sng == "oom"){
    std::cout << "Testing out of memory (OOM) behavior and my_new_handler()..." << std::endl;
    std::cout << "WARNING: This test may cause OOM condition and huge core file..." << std::endl;
    /* Instead of throwing bad_alloc(), call homebrew handler when new[] fails
       my_new_handler() defined in libcsz_c++:utl.cc */
    std::set_new_handler(my_new_handler); // [fnc] Handle failures thrown by new[]
    const long mmr_blk_nbr(1000000); // [nbr] Number of blocks to attempt to allocate
    const long mmr_blk_sz(100000000); // [B] Size of each block
    const float mmr_blk_ttl_max(100.0e9); // [B] Maximum amount to try to allocate
    char *mmr_blk_ptr[mmr_blk_nbr]; // [ptr] Pointer for each block
    long mmr_blk_cnt(0L); // [nbr] Count of blocks allocated
#ifdef HAVE_LONG_LONG
    // long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
    long long mmr_blk_ttl(0LL); // Cumulative allocated memory
#else // !HAVE_LONG_LONG
    long mmr_blk_ttl(0L); // [B] Cumulative allocated memory
#endif // !HAVE_LONG_LONG
    for(idx=0;idx<mmr_blk_nbr;idx++){
      mmr_blk_ptr[idx]=new char[mmr_blk_sz];
      mmr_blk_cnt++;
      mmr_blk_ttl+=mmr_blk_sz;
      std::cout << idx+1 << " blocks of size " << mmr_blk_sz/1.0e6 << " MB = " << mmr_blk_ttl/1.0e6 << " MB total" << std::endl;
      /* AIX allocates at least 7000 GB without complaining
	 Stop at 100 GB rather than bringing down the system
	 Allocating more than 100 GB is currently pointless */
      if(mmr_blk_ttl > mmr_blk_ttl_max){
	std::cout << "Successfully allocated more than " << mmr_blk_ttl_max/1.0e6 << " MB, ending OOM test" << std::endl;
	break;
      } // !mmr_blk_ptr
    } // !idx
    for(idx=0;idx<mmr_blk_cnt;idx++){
      delete []mmr_blk_ptr[idx];
    } // !idx
  } // !dbg || tst_sng == "oom"
  
  if(dbg_lvl == dbg_old || tst_sng == "phys_cst"){
    std::cout << "Testing physical constants..." << std::endl;
    using phc::dns_dst_std; // (2500.0) [kg m-3] Standard density of dust
    using phc::dns_Fe2O3; // (5260.0) [kg m-3] Density of Fe2O3
    using phc::dns_SiO2; // (2620.0) [kg m-3] Density of SiO2
    using phc::mmw_C; // (12.011e-03) [kg mol-1] Mean molecular weight of C IUPAC
    using phc::mmw_CaCO3; // (100.087e-03) [kg mol-1] Mean molecular weight of CaCO3
    using phc::mmw_Fe2O3; // () [kg mol-1] Mean molecular weight of Fe2O3
    using phc::mmw_HO2; // (33.0067e-03) [kg mol-1] Mean molecular weight of HO2
    using phc::mmw_MgSO4; // (120.369e-03) [kg mol-1] Mean molecular weight of MgSO4
    using phc::mmw_N2O5; // (108.01e-03) [kg mol-1] Mean molecular weight of N2O5
    using phc::mmw_NH4; // (18.0385e-03) [kg mol-1] Mean molecular weight of NH4
    using phc::mmw_NO3; // (62.0049e-03) [kg mol-1] Mean molecular weight of NO3
    using phc::mmw_NaCl; // (58.4425e-03) [kg mol-1] Mean molecular weight of NaCl
    using phc::mmw_SO4; // (96.0636e-03) [kg mol-1] Mean molecular weight of SO4
    using phc::mmw_SiO2; // (60.0843e-03) [kg mol-1] Mean molecular weight of SiO2
    using phc::spc_heat_SiO2_sld; // () [J kg-1 K-1] Specific heat capacity of quartz CRC95 p. 5-21
    using phc::spc_heat_Fe2O3_sld; // () [J kg-1 K-1] Specific heat capacity of hematite CRC95 p. 5-21
    std::cout 	 
      << "dns_dst_std = " << dns_dst_std << " kg m-3, "
      << "dns_Fe2O3 = " << dns_Fe2O3 << " kg m-3, "
      << "dns_SiO2 = " << dns_SiO2 << " kg m-3, "
      << "mmw_CaCO3 = " << mmw_CaCO3 << " kg mol-1, " 
      << "mmw_Fe2O3 = " << mmw_Fe2O3 << " kg mol-1, "
      << "mmw_HO2 = " << mmw_HO2 << " kg mol-1, "
      << "mmw_MgSO4 = " << mmw_MgSO4 << " kg mol-1, "
      << "mmw_N2O5 = " << mmw_N2O5 << " kg mol-1, "
      << "mmw_NH4 = " << mmw_NH4 << " kg mol-1, "
      << "mmw_NO3 = " << mmw_NO3 << " kg mol-1, "
      << "mmw_NaCL = " << mmw_NaCl << " kg mol-1, "
      << "mmw_SO4 = " << mmw_SO4 << " kg mol-1, "
      << "mmw_SiO2 = " << mmw_SiO2 << " kg mol-1, "
      << "spc_heat_Fe2O3_sld = " << spc_heat_Fe2O3_sld << " J kg-1 K-1, "
      << "spc_heat_SiO2_sld = " << spc_heat_SiO2_sld << " J kg-1 K-1, "
      << std::endl;
  } // !dbg || tst_sng == "phys_cst"

  if(dbg_lvl == dbg_old || tst_sng == "prn"){
    // Test counter printing
    std::cout << "Testing printing gymnastics..." << std::endl;
    std::cout << "Set loop counter with --int_foo, e.g., ccc --tst=prn --int_foo=9" << std::endl;
    std::cout << "Wavelength index current: ";
    std::string bck_sng("\b\b\b\b\b\b");
    for(idx=0;idx<int_foo;idx++){
      if(idx > 0) std::cout << bck_sng;
      //      std::cout << std::setw(6) << idx;
      std::cout << "\rMie computation: " << std::setw(6) << idx;
    } // end loop over idx
    std::cout << std::endl;
  } // !dbg || tst_sng == "prn"

  if(dbg_lvl == dbg_old || tst_sng == "ptr"){
    std::cout << "Testing pointer nomenclature..." << std::endl;
    std::cout << "static const int * const ( * ( * const p [10] ) (const int *) ) ( void ): p is a static array of ten constant pointers to functions that pass one argument of type const int * and return a pointer to a function that passes no arguments and returns a constant pointer to a constant int. This demonstrates the \"Right-Left\" rule: 1. Locate the identifier in the expression. It is not a keyword, operator or anything extraneous (not required). 2) Look to the Right, and interpret (and keep going Right), unless you encounter a \')\' or are done, then go Left. Step 3) When going Left, interpret (and keep going Left), unless you encounter a \'(\', then go Right." << std::endl;
    std::cout << "int * const foo: foo is a constant pointer to an integer. foo cannot change, but the integer(s) can, default for arrays." << std::endl;
    std::cout << "int &foo: foo is a reference to an integer. If foo changes, the integer changes." << std::endl;
    std::cout << "float *&foo: foo is a reference to a float pointer. Change foo to change the original float pointer." << std::endl;
    std::cout << "const int * foo: foo is a pointer to an integer constant. foo cannot be used to change the integer" << std::endl;
    std::cout << "const int * const foo: foo is a constant pointer to an integer constant. Neither foo nor the integer may be changed." << std::endl;
    std::cout << "float (*foo)(const float *bar): foo is a pointer to a function that takes a single float pointer argument and returns a float" << std::endl;
    std::cout << "float *(*foo)(const float &bar): foo is a pointer to a function that takes a float argument (by reference) and returns a pointer to a float" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "int * restrict foo: foo is a restricted pointer to an integer. If the integer changes, no other pointer will access the integer" << std::endl;
    std::cout << "int * restrict const foo: foo is a constant restricted pointer to an integer. If the integer changes, no other pointer will access the integer" << std::endl;
  } // !dbg || tst_sng == "ptr"

  if(dbg_lvl == dbg_old || tst_sng == "rnd" || tst_sng == "round" ){
    std::cout << "Testing round() functionality..." << std::endl;
    std::cout << "Testing rounding ..." << std::endl;
    std::cout << "round(-1.0) = " << mth_fnc::round(-1.0) << std::endl;
    std::cout << "round(-0.5) = " << mth_fnc::round(-0.5) << std::endl;
    std::cout << "round(-0.5+1.0e-17) = " << mth_fnc::round(-0.5+1.0e-17) << std::endl;
    std::cout << "round(-0.5+1.0e-16) = " << mth_fnc::round(-0.5+1.0e-16) << std::endl;
    std::cout << "round(0.0) = " << mth_fnc::round(0.0) << std::endl;
    std::cout << "round(0.5-1.0e-16) = " << mth_fnc::round(0.5-1.0e-16) << std::endl;
    std::cout << "round(0.5-1.0e-17) = " << mth_fnc::round(0.5-1.0e-17) << std::endl;
    std::cout << "round(0.5) = " << mth_fnc::round(0.5) << std::endl;
    std::cout << "round(dbl_foo) = round(" << dbl_foo << ") = " << mth_fnc::round(dbl_foo) << std::endl;
  } // !dbg || tst_sng == "rnd"

  if(dbg_lvl == dbg_old || tst_sng == "rtti"){
    std::cout << "Testing run-time-type-instantiation (RTTI)..." << std::endl;
    std::cout << "Homebrew typeid2sng() function reporting typeid.name() of all types..." << std::endl;
    /* Suffixes for integer and floating point constants DeD01 p. 925
       Double precision is assumed default for real literal ("naked") constants
       u,l,f qualifiers are case-insensitive 
       u following integer specifies unsigned
       l following integer specifies long int
       ul following integer specifies unsigned long
       l following decimal point sets long double
       ll following integer specifies long long
       ull following integer specifies unsigned long long */
    std::cout << "typeid(3.14159).name() = " << typeid(3.14159).name() << std::endl;
    std::cout << "typeid(3.14159f).name() = " << typeid(3.14159f).name() << std::endl;
    std::cout << "typeid(3.14159l).name() = " << typeid(3.14159l).name() << std::endl;
    std::cout << "typeid(3).name() = " << typeid(3).name() << std::endl;
    std::cout << "typeid(3u).name() = " << typeid(3u).name() << std::endl;
    std::cout << "typeid(3l).name() = " << typeid(3l).name() << std::endl;
    std::cout << "typeid(3ul).name() = " << typeid(3ul).name() << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "typeid(3ll).name() = " << typeid(3ll).name() << std::endl;
    std::cout << "typeid(3ull).name() = " << typeid(3ull).name() << std::endl;
#endif // !HAVE_LONG_LONG
    // Define all other types in terms of command line arguments flt_foo and int_foo
    const bool bln_foo(true);
    const char chr_foo(int_foo);
    const double dbl_foo(flt_foo);
    //    const float flt_foo(flt_foo);
    //    const int int_foo(int_foo);
    const long double ldbl_foo(flt_foo);
    const long lng_foo(int_foo);
#ifdef HAVE_LONG_LONG
    const long long lng_lng_foo(int_foo);
#endif // !HAVE_LONG_LONG
    const short sht_foo(int_foo);
    const signed char schr_foo(int_foo);
    const std::complex<double> dcx_foo(0.0,0.0);
    const std::complex<float> fcx_foo(0.0,0.0);
    const unsigned char uchr_foo(std::abs(int_foo));
    const unsigned int uint_foo(std::abs(int_foo));
#ifdef HAVE_LONG_LONG
    const unsigned long long ulng_lng_foo(std::abs(int_foo));
#endif // !HAVE_LONG_LONG
    const unsigned long ulng_foo(std::abs(int_foo));
    const unsigned short usht_foo(int_foo);
    std::cout << "typeid(bool).name() = " << typeid(bln_foo).name() << std::endl;
    std::cout << "typeid(char).name() = " << typeid(chr_foo).name() << std::endl;
    std::cout << "typeid(double).name() = " << typeid(dbl_foo).name() << std::endl;
    std::cout << "typeid(float).name() = " << typeid(flt_foo).name() << std::endl;
    std::cout << "typeid(int).name() = " << typeid(int_foo).name() << std::endl;
    std::cout << "typeid(long double).name() = " << typeid(ldbl_foo).name() << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "typeid(long long).name() = " << typeid(lng_lng_foo).name() << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "typeid(long).name() = " << typeid(lng_foo).name() << std::endl;
    std::cout << "typeid(short).name() = " << typeid(sht_foo).name() << std::endl;
    std::cout << "typeid(signed char).name() = " << typeid(schr_foo).name() << std::endl;
    std::cout << "typeid(std::complex<double>).name() = " << typeid(dcx_foo).name() << std::endl;
    std::cout << "typeid(std::complex<float>).name() = " << typeid(fcx_foo).name() << std::endl;
    std::cout << "typeid(unsigned char).name() = " << typeid(uchr_foo).name() << std::endl;
    std::cout << "typeid(unsigned int).name() = " << typeid(uint_foo).name() << std::endl;
#ifdef HAVE_LONG_LONG
    std::cout << "typeid(unsigned long long).name() = " << typeid(ulng_lng_foo).name() << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "typeid(unsigned long).name() = " << typeid(ulng_foo).name() << std::endl;
    std::cout << "typeid(unsigned short).name() = " << typeid(usht_foo).name() << std::endl;

/* RTTI results from various compilers
   NB: Certain algebraic solutions in mth.hh depend on RTTI for complex types
   
   AIX xlC RTTI:
   typeid(bool).name() = bool
   typeid(char).name() = char
   typeid(double).name() = double
   typeid(float).name() = float
   typeid(int).name() = int
   typeid(long double).name() = long double
   typeid(long long).name() = long long
   typeid(long).name() = long
   typeid(short).name() = short
   typeid(signed char).name() = signed char
   typeid(std::complex<double>).name() = std::complex<double>
   typeid(std::complex<float>).name() = std::complex<float>
   typeid(unsigned char).name() = unsigned char
   typeid(unsigned int).name() = unsigned int
   typeid(unsigned long long).name() = unsigned long long
   typeid(unsigned long).name() = unsigned long
   typeid(unsigned short).name() = unsigned short
   
   GNU g++ 3.4.2 RTTI (g++ changes RTTI all the time):
   typeid(bool).name() = b
   typeid(char).name() = c
   typeid(double).name() = d
   typeid(float).name() = f
   typeid(int).name() = i
   typeid(long double).name() = e
   typeid(long long).name() = x
   typeid(long).name() = l
   typeid(short).name() = s
   typeid(signed char).name() = a
   typeid(std::complex<double>).name() = St7complexIdE
   typeid(std::complex<float>).name() = St7complexIfE
   typeid(unsigned char).name() = h
   typeid(unsigned int).name() = j
   typeid(unsigned long long).name() = y
   typeid(unsigned long).name() = m
   typeid(unsigned short).name() = t */
  } // !dbg || tst_sng == "rtti"

  if(dbg_lvl == dbg_old || tst_sng == "simd"){
    // ccc --dbg=1 --tst=simd
    std::cout << "Testing OpenMP SIMD implementation..." << std::endl; 
#ifdef _OPENMP // OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification
#else // !_OPENMP
    std::cout << "OpenMP Activation token _OPENMP is false" << std::endl; 
    std::cout << prg_nm_get() << ": INFO Not attempting OpenMP threading" << std::endl;
#endif // !_OPENMP
#ifdef _OPENMP // OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification
#endif // endif _OPENMP
  } // !dbg || tst_sng == "simd"

  if(dbg_lvl == dbg_old || tst_sng == "sizeof"){
    std::cout << "Testing dynamic sizing ability..." << std::endl;
    std::cout << "sizeof is compile time operator, not runtime function DeD01 p. 324" << std::endl;
    const std::string grd_typ_sng[]={
      "regular",
      "CAM_SW",
      "CAM_LW"};
    const int grd_typ_nbr((sizeof grd_typ_sng)/(sizeof(std::string)));
    std::cerr << "sizeof grd_typ_sng = " << sizeof grd_typ_sng << " B" << std::endl;
    std::cerr << "sizeof(std::string) = " << sizeof(std::string) << " B" << std::endl;
    std::cerr << "grd_typ_nbr = " << grd_typ_nbr << std::endl;
  } // !dbg || tst_sng == "sizeof"
  
  if(dbg_lvl == dbg_old || tst_sng == "spc_slr"){
    std::cout << "Testing spc_slr class..." << std::endl;
    // Instantiate solar flux source
    spc_slr_cls flx_slr_src; // [sct] Solar flux source
    // Print solar flux source object
    std::cout << flx_slr_src; // [sct] Solar flux source
    // Test solar flux source class for memory leaks by initializing lng_foo members
    rcd+=spc_slr_cls::tst(lng_foo); // [sct] Solar flux source
  } // !dbg || tst_sng == "spc_slr"

  if(dbg_lvl == dbg_old || tst_sng == "sng" || tst_sng == "srm"){
    /* Test string and stream handling DeD97 p. 915, DeD01 p. 960
       ccc --tst=srm --sz_nbr=15 --dbl_foo=73.0e-73 --flt_foo=37.0e-37 --sng_foo=9223372036854775807
       Following two lines should work in ANSI C++ when istringstreams are supported
       istringstream idx_rfr_srm(idx_rfr_sng); // See DeD97 p. 915, DeD01 p. 960
       idx_rfr_srm >> idx_rfr_rl_usr >> idx_rfr_img_usr; */
    std::cout << "Testing string and stream handling..." << std::endl;
    std::ostringstream sng_srm_out; // [srm] Output string stream
    std::stringstream sng_srm_in(sng_foo); // [srm] Input string stream
    std::string sng_val; // [sng] Value expressed as string 
    sng_val=nbr2sng(dbl_foo); // [sng] Value expressed as string 
#ifdef HAVE_LONG_LONG
    // long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
    if(sng_foo.size()==0L) std::cout << "HINT: Set sng_foo to a number or conversion functions will fail" << std::endl;
    if(sng_srm_in >> lng_lng_foo) std::cout << "String \"" << sng_foo << "\" converted to atomic type long long: " << lng_lng_foo << std::endl; else err_prn(prg_nm,sbr_nm,"Unable to convert string \""+sng_foo+"\" to atomic type long long");
    std::cout << "String \"" << sng_foo << "\" converted to atomic type unsigned long long: " << sng2nbr(sng_foo,ulng_lng_foo) << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "Variable ubyte_foo has value = " << ubyte_foo << " after command line processing" << std::endl;
    std::cout << "String \"" << sng_foo << "\" converted to atomic type unsigned char: " << sng2nbr(sng_foo,ubyte_foo) << std::endl;
    if(sng_srm_out << dbl_foo) sng_val=sng_srm_out.str(); else err_prn(prg_nm,sbr_nm,"Unable to convert float to string");
    std::cout << "Value " << dbl_foo << " expressed as string: \"" << sng_val << "\"" << std::endl;
    std::cout << "Value " << dbl_foo << " expressed with nbr2sng(dbl_foo): \"" << nbr2sng(dbl_foo) << "\"" << std::endl;
    std::cout << "Value " << flt_foo << " expressed with nbr2sng(flt_foo,prc_cmp_dcm_plc): \"" << nbr2sng(flt_foo,prc_cmp_dcm_plc) << "\"" << std::endl;
    std::cout << "Value " << dbl_foo << " expressed with nbr2sng(dbl_foo,prc_cmp_dcm_plc): \"" << nbr2sng(dbl_foo,prc_cmp_dcm_plc) << "\"" << std::endl;

    // Demonstrate stream manipulation by printing a table
    std::cout << "Printing a table..." << std::endl;

    // floatfield is mask for fixed, scientific, None (default)
    // std::cout.unsetf(std::ios::floatfield); // [msk] Mask: fixed, scientific
    std::cout.setf(std::ios::fixed); // [flg] Fixed number of digits to right of decimal
    // std::cout.setf(ios::scientific); // [flg] Print in scientific notation
    std::cout.setf(std::ios::showpoint); // [flg] Always output floats with decimal point

    // adjustfield is mask for left, right, internal, None (default)
    // std::cout.unsetf(std::ios::adjustfield); // [msk] Mask: None, left, right, internal
    // std::cout.setf(std::ios::left); // Left justify output
    std::cout.setf(std::ios::right); // [flg] Right justify output
    // std::cout.setf(std::ios::internal); // Left-justify sign, right-justify magnitude

    std::cout.setf(std::ios::uppercase); // [flg] Uppercase letters
    std::cout.setf(std::ios::showpos); // [flg] Write + sign for positive numbers
    std::cout.precision(3); // Number of digits to right of decimal
    std::cout << "0         1         2         3         4         5         6         7         " << std::endl;
    std::cout << "01234567890123456789012345678901234567890123456789012345678901234567890123456789" << std::endl;
    /* Table style: Set variable with all table-wide flags first
       Format exceptions in place, then reset to defaults using variable
       Streams should have defaults set on exit from table
       Use setw() once per desired element (it is not "sticky")
       Use setiosflags() and resetiosflags() to embed flags interactively */
    std::ios::fmtflags ios_flg_dfl(std::cout.flags()); // [flg] I/O stream defaults
    std::ios::fmtflags ios_flg_tbl(std::ios::fixed | std::ios::showpoint | std::ios::right); // [flg] I/O stream table format
    int tbl_clm_wdt(8); // [nbr] Width of table columns
    // Initialize table format
    std::cout.setf(ios_flg_tbl); // [srm] Output stream
    std::cout << std::setw(3) << "idx" << std::setw(tbl_clm_wdt) << "sz" << std::setw(tbl_clm_wdt) << "sz" << std::setw(tbl_clm_wdt) << "sz" << std::setw(tbl_clm_wdt) << "sz" << std::setw(tbl_clm_wdt) << "sz" << std::setw(tbl_clm_wdt) << "sz" << std::endl;
    std::cout << std::setw(3) << "#" << std::setw(tbl_clm_wdt) << "um" << std::setw(tbl_clm_wdt) << "um" << std::setw(tbl_clm_wdt) << "um" << std::setw(tbl_clm_wdt) << "um" << std::setw(tbl_clm_wdt) << "um" << std::setw(tbl_clm_wdt) << "um" << std::endl;
    for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
      std::cout << std::setiosflags(ios_flg_tbl) << std::setprecision(3) << std::setw(3) << sz_idx << std::setw(tbl_clm_wdt) << sz[sz_idx]*1.e6 << std::setiosflags(std::ios::showpos) << std::setw(tbl_clm_wdt) << sz[sz_idx]*1.0e6+10.0 << std::resetiosflags(std::ios::showpos) << std::setw(tbl_clm_wdt) << -sz_nbr+2*sz_idx << std::setprecision(2) << std::setw(tbl_clm_wdt) << sz[sz_idx]*1.0e6*pow(2.0,static_cast<double>(sz_idx)) << std::setprecision(3) << std::setw(tbl_clm_wdt) << sz[sz_idx]*1.0e6 << std::setiosflags(std::ios::scientific) << std::setprecision(3) << std::setw(tbl_clm_wdt) << sz[sz_idx]*1.0e6 << std::setiosflags(ios_flg_tbl) << std::endl;
    } // end loop over sz
    // Restore defaults
    std::cout.precision(8); // Number of digits to right of decimal
    std::cout.flags(ios_flg_dfl); // [srm] Output stream
  } // !dbg || tst_sng == "sng" || tst_sng == "srm"

  if(dbg_lvl == dbg_old || tst_sng == "sizeof" || tst_sng == "sz"){
    std::cout << "Testing native type sizes..." << std::endl;
    std::cout << "sizeof(bool) = " << sizeof(bool) << std::endl;
    std::cout << "sizeof(char) = " << sizeof(char) << std::endl;
    std::cout << "sizeof(double) = " << sizeof(double) << std::endl;
    std::cout << "sizeof(float) = " << sizeof(float) << std::endl;
    std::cout << "sizeof(int) = " << sizeof(int) << std::endl;
    std::cout << "sizeof(long double) = " << sizeof(long double) << std::endl;
#ifdef HAVE_LONG_LONG
    // long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
    std::cout << "sizeof(long long) = " << sizeof(long long) << std::endl;
#endif // !HAVE_LONG_LONG
    std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
    std::cout << "sizeof(short) = " << sizeof(short) << std::endl;
    std::cout << "sizeof(size_t) = " << sizeof(size_t) << std::endl;
  } // !dbg || tst_sng == "sizeof"
  
  if(dbg_lvl == dbg_old || tst_sng == "c99"){
    /* Test C99 syntax
       Try to keep these tests consistent with c.c
       Many tests from http://www.comeaucomputing.com/features.html */
    std::cout << "Testing C99 syntax..." << std::endl;

    // Designated initializers
#if 0 // HAVE_C99
    std::cout << "Testing designated initializers..." << std::endl;
    struct sct_typ{int lmn_int;float lmn_flt;int lmn_int_arr[3];};
    std::cout << "Initialize all structure members in forward order..." << std::endl;
#if defined(_AIX) // Compiler is xlc
    // xlc does not allow specifying range of elements to initialize
    struct sct_typ sct_1={.lmn_int=3,.lmn_flt=3.123,.lmn_int_arr[0]=-3};
#else // Other compiler (GCC...)
    // GCC allows this non-ISO construct
    // c.c:222: warning: ISO C forbids specifying range of elements to initialize
    struct sct_typ sct_1={.lmn_int=3,.lmn_flt=3.123,.lmn_int_arr[0 ... 2]=-3};
#endif // !_AIX
    std::cout << "sct_1.lmn_int = " << sct_1.lmn_int <<  "sct_1.lmn_flt = " << sct_1.lmn_flt << "sct_1.lmn_int_arr[0] = " << sct_1.lmn_int_arr[0] << std::endl;
    std::cout << "Initialize some structure members in reverse order..." << std::endl;
    struct sct_typ sct_2={.lmn_flt=3.123,.lmn_int=3};
    std::cout << "Initialize individual array elements within structure..." << std::endl;
    struct sct_typ sct_3={.lmn_int_arr[2]=2};
    std::cout << "Initialize individual array elements of plain array..." << std::endl;
    int int_arr_dsg_ntl[3]={[0]=1,[2]=3};
    struct sct_nst_typ{int lmn_int;struct sct_typ sct_typ_lmn;};
    std::cout << "Nested initialization of elements of structure and sub-structure..." << std::endl;
    struct sct_nst_typ sct_nst={.lmn_int=11,.sct_typ_lmn.lmn_int_arr={0,1,2}};
    sct_1=sct_1; // CEWI
    sct_2=sct_2; // CEWI
    sct_3=sct_3; // CEWI
    //int_arr_dsg_ntl=int_arr_dsg_ntl; // CEWI
    sct_nst=sct_nst; // CEWI
#else // !HAVE_C99
    std::cout << "INFO: Non-C99 compiler, not testing designated initializers..." << std::endl;
#endif // !HAVE_C99

    // Compound literals    
#if 0 // HAVE_C99
    std::cout << "Testing compound literals..." << std::endl;
    std::cout << "Define un-named int array..." << std::endl;
    int *int_arr=(int[]){1,2}; // Compound literal
    std::cout << "Define const un-named int array..." << std::endl;
    const int *cst_int_arr=(const int[]){1,2}; // Compound literal
    std::cout << "Define pointer to int value 1..." << std::endl;
    int *int_ptr=&(int){1}; // Compound literal
    std::cout << "Define const pointer to int value 1..." << std::endl;
    const int *cst_int_ptr=&(const int){1}; // Compound literal
    std::cout << "Define const string pointer to value \"Compound literal\"..." << std::endl;
    const char *cst_sng=(const char[]){"Compound literal"}; // Compound literal
#else // !HAVE_C99
    std::cout << "INFO: Non-GCC compiler, not testing compound literals..." << std::endl;
#endif // !HAVE_C99

    // Restrict type qualifier
#if 1 // HAVE_C99
    /* Restrict type qualifier. See, e.g., IBM C/C++ Guide p. 69.
       restrict is a type qualifier.
       "restrict" promises that if memory pointed to by "restrict"-qualified pointer is modified,
       no other pointer will access that same memory.
       This is a promise the programmer must keep.
       It is an assumption the compiler makes, but cannot enforce
       Memory that is not modified may be accessed (read) by multiple restricted pointers
       
       Understanding restrict is difficult
       Some useful discussions are at
       http://developers.sun.com/tools/cc/articles/cc_restrict.html 
       is the most useful restrict illustration found so far.
       It has good examples of all flavors of restrict: prototypes, structures, etc. 

       Other sources are:
       http://www.lysator.liu.se/c/restrict.html
       http://www.cbau.freeserve.co.uk/Compiler/RestrictPointers.html

       IBM xlC consumes and ignores "restrict" for compatibility with C99 
       
       GCC g++ uses `__restrict__', or `__restrict'
       GCC requires restrict in function definitions, but not in function prototypes */
    std::cout << "Testing restrict keyword..." << std::endl;
    // Substitute whitespace for __restrict in all compilers except g++
#ifndef __GNUG__
    // Define __restrict as empty string on non-g++ compilers. May be un-safe. Generally dangerous to fiddle with compiler __variables.
#define __restrict
    std::cout << "INFO: Compiler is not g++, substituting whitespace for __restrict" << std::endl;
#endif // !__GNUG__
    int * __restrict int_ptr_rst_1=new int[bnd_nbr];
    int * __restrict int_ptr_rst_2;
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      int_ptr_rst_1[bnd_idx]=bnd_idx;
    } // end loop over bnd
    // Assign one restricted pointer to another. Allowed
    int_ptr_rst_2=int_ptr_rst_1;
    // Access already modified value through second restricted pointer. Dis-allowed
    int_ptr_rst_2[0]=73; // Dangerous. Should cause compiler warning?
    delete []int_ptr_rst_1;
#else // !HAVE_C99
    std::cout << "INFO: Compiler is not xlC or g++, not testing restrict type qualifier..." << std::endl;
#endif // !HAVE_C99
  } // !dbg || tst_sng == "c99"

  if(dbg_lvl == dbg_old || tst_sng == "stat"){
    // ccc --tst=stat --fl_in=~/foo.nc
    // Test stat() functions
    // NB: stat() is system-specific and structure element names may need to be ifdef'd to work on non-Linux system
    int rcd_sys;
    struct stat stat_sct;
    rcd_sys=stat(fl_in.c_str(),&stat_sct);
    std::cout << "Testing stat() functions on " << fl_in << ":" << std::endl;
    std::cout << "Return code = rcd_sys = " << rcd_sys << std::endl;
    std::cout << "Inode number = st_ino = " << stat_sct.st_ino << std::endl;
    std::cout << "Protection mode = st_mode = " << stat_sct.st_mode << std::endl;
    std::cout << "Number of hard links = st_nlink = " << stat_sct.st_nlink << std::endl;
    std::cout << "User ID of owner = st_uid = " << stat_sct.st_uid << std::endl;
    std::cout << "Group ID of owner = st_gid = " << stat_sct.st_gid << std::endl;
    std::cout << "Total size = st_size = " << stat_sct.st_size << " B" << std::endl;
    std::cout << "File system I/O blocksize = st_blksize = " << stat_sct.st_blksize << " B" << std::endl;
    std::cout << "Number of 512B blocks allocated = st_blocks = " << stat_sct.st_blocks << std::endl;
    std::cout << "Time of last access = st_atime = " << stat_sct.st_atime << std::endl;
    std::cout << "Time of last modification = st_mtime = " << stat_sct.st_mtime << std::endl;
    std::cout << "Time of last ctime = st_ctime = " << stat_sct.st_ctime << std::endl;
  } // !dbg || tst_sng == "stat"

  if(dbg_lvl == dbg_old || tst_sng == "sys"){
    // ccc --tst=sys
    // Test systems functions
    std::cout << "Testing systems functions..." << std::endl;
    std::cout << "User ID = getuid() = " << getuid() << std::endl;
    std::cout << "Effective user ID = geteuid() = " << geteuid() << std::endl;
    std::cout << "Process ID = getpid() = " << getpid() << std::endl;
    std::cout << "Parent process ID = getppid() = " << getppid() << std::endl;
  } // !dbg || tst_sng == "sys"

  if(dbg_lvl == dbg_old || tst_sng == "snd"){
    // ccc --tst=snd
    // Test speed of sound algorithm
    // http://en.wikipedia.org/wiki/Speed_of_sound#Seawater
    // c(T, S, z) = a1 + a2T + a3T2 + a4T3 + a5(S - 35) + a6z + a7z2 + a8T(S - 35) + a9Tz3 
    /* cat > ~/spd_snd_wtr.nco << EOF
       tpt=303.16;tpt_cls=tpt-273.16;sln=0.035;sln_ppt=sln*1000;dpt=1000;a1=1448.96;a2=4.591;a3=-5.304e-2;a4=2.374e-4;a5=1.340;a6=1.630e-2;a7=1.675e-7;a8=-1.025e-2;a9=-7.139e-13;c_wtr=a1+a2*tpt_cls+a3*tpt_cls^2+a4*tpt_cls^3+a5*(sln_ppt-35.0)+a6*dpt+a7*dpt^2+a8*tpt_cls*(sln_ppt-35.0)+a9*tpt_cls*dpt^3;
EOF
       ncap2 -O -v -S ~/spd_snd_wtr.nco ~/nco/data/in.nc ~/foo.nc */
    std::cout << "Testing speed of sound algorithm..." << std::endl;
    prc_cmp dpt=1000.0; // [m] Depth
    prc_cmp tpt=303.16; // [K] Temperature
    prc_cmp tpt_cls=tpt-tpt_frz_pnt; // [C] Temperature celsius
    prc_cmp sln=0.035; // [ppt] Salinity
    prc_cmp sln_ppt=sln*1000; // [ppt] Salinity
    std::cout << "Temperature = " << tpt_cls << " C, Salinity = " << sln_ppt << " ppt, Depth = " << dpt << " m, Speed of sound = " << spd_snd_wtr_fst_scl(dpt,sln_ppt,tpt_cls) << " m s-1" << std::endl;
  } // !dbg || tst_sng == "snd"

  if(dbg_lvl == dbg_old || tst_sng == "try"){
#ifndef PGI_CXX
    std::cout << "Testing try/throw/catch exceptions..." << std::endl;
    try{
      if(dbg_lvl == 0){
	throw 1;
      } // endif
    }catch(...){ // Catch all exceptions
      std::cout << "Caught exception generated by testing dbg_lvl == 0" << std::endl;
    } // end try
#else // PGI_CXX
    std::cout << "PGI compiler pgCC does not support try/throw/catch exceptions, skipping test..." << std::endl;
#endif // PGI_CXX
  } // !dbg || tst_sng == "try"

  if(dbg_lvl == dbg_old || tst_sng == "tst"){
    Test test;
    std::cout << "Testing Test class..." << std::endl;
    std::cout << "Default constructor no argument test.tst_nbr_get() = " << test.tst_nbr_get() << std::endl;
    Test test2(4);
    std::cout << "Default constructor one argument test2.tst_nbr_get() = " << test2.tst_nbr_get() << std::endl;
  } // !dbg || tst_sng == "tst"

  // Allocate dynamic arrays
  prc_cmp *lat=new prc_cmp[lat_nbr]; // [dgr] Latitude
  prc_cmp *lon=new prc_cmp[lon_nbr]; // [dgr] Longitude

  /* Create regular grid for appending to SeaWiFS data files
     SeaWiFs level 3 monthly data comes on 2048x4096 grid
     Grid spacing is same for both Lat and Lon = 180.0/2048 = 0.087890625
     Latitude runs from South to North
     Longitudes runs Eastward starting from date line
     SW point is centered at -89.9560546875 S, -179.9560546875 E

     Some flexibility has been included for renaming dimensions and changing dimension sizes
     Thus the code below is easy to change for arbitrary, even irregular, grids
     SeaWiFS grid is obtained by specifying following command line arguments
     cd ~/c++;ccc --lat_nm="lat" --lon_nm="lon" --lat_nbr=2048 --lon_nbr=4096;cd -
     Append SeaWiFS grid to SeaWiFS HDF file using
     ncks -A -v lat,lon ~/c++/foo.nc ~/nco/data/tst.hdf
  */
  for(lat_idx=0;lat_idx<lat_nbr;lat_idx++){
    lat[lat_idx]=-89.9560546875+(static_cast<prc_cmp>(lat_idx)/static_cast<prc_cmp>(lat_nbr))*180.0; // [dgr] Latitude (SeaWiFS)
    lat[lat_idx]=-89.958333333333333333+(static_cast<prc_cmp>(lat_idx)/static_cast<prc_cmp>(lat_nbr))*180.0; // [dgr] Latitude (GFDL TerrainBase)
  } // end loop over lat
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    lon[lon_idx]=-179.9560546875+(static_cast<prc_cmp>(lon_idx)/static_cast<prc_cmp>(lon_nbr))*360.0; // [dgr] Longitude (SeaWiFS)
    lon[lon_idx]=0.0416666666666666666667+(static_cast<prc_cmp>(lon_idx)/static_cast<prc_cmp>(lon_nbr))*360.0; // [dgr] Longitude (UCSB & GFDL TerrainBase)
  } // end loop over lon

  // Map grid
  std::string edg_est_nm("EDGEE"); // [sng] East edge of grid variable name
  std::string edg_nrt_nm("EDGEN"); // [sng] North edge of grid variable name
  std::string edg_sth_nm("EDGES"); // [sng] South edge of grid variable name
  std::string edg_wst_nm("EDGEW"); // [sng] West edge of grid variable name
  std::string lat_ctr_2d_nm("LATIXY"); // [sng] Latitude at gridcell center variable name
  std::string lon_ctr_2d_nm("LONGXY"); // [sng] Longitude at gridcell center variable name
  prc_cmp edg_nrt(90.0); // [dgr] North edge of grid
  prc_cmp edg_est(360.0); // [dgr] East edge of grid
  prc_cmp edg_sth(-90.0); // [dgr] South edge of grid
  prc_cmp edg_wst(0.0); // [dgr] West edge of grid
  prc_cmp edg_grd[4]={edg_nrt,edg_est,edg_sth,edg_wst}; // [dgr] Grid edges (north, east, south, west)
  a2d_cls<prc_cmp> lat_ctr_2d(lat_nbr,lon_nbr); // [dgr] Latitude at gridcell center
  a2d_cls<prc_cmp> lon_ctr_2d(lat_nbr,lon_nbr); // [dgr] Longitude at gridcell center
  // Fill in 2D coordinate grid
  for(lat_idx=0;lat_idx<lat_nbr;lat_idx++){
    for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
      lat_ctr_2d(lat_idx,lon_idx)=lat[lat_idx]; // [dgr] Latitude at gridcell center
      lon_ctr_2d(lat_idx,lon_idx)=lon[lon_idx]; // [dgr] Longitude at gridcell center
    } // end loop over lon
  } // end loop over lat

  // Open output file
  int nccreate_mode(NC_CLOBBER); // [enm] Mode flag for nco_create() call
#ifdef ENABLE_NETCDF4
  if(fl_out_fmt==NCO_FORMAT_UNDEFINED) fl_out_fmt=NC_FORMAT_NETCDF4;
  if(fl_out_fmt == NC_FORMAT_64BIT){
    nccreate_mode|=NC_64BIT_OFFSET;
  }else if(fl_out_fmt == NC_FORMAT_NETCDF4){
    nccreate_mode|=NC_NETCDF4;
  }else if(fl_out_fmt == NC_FORMAT_NETCDF4_CLASSIC){
    nccreate_mode|=NC_NETCDF4|NC_CLASSIC_MODEL;
  } /* end else fl_out_fmt */
#else // !ENABLE_NETCDF4
  if(fl_out_fmt==NCO_FORMAT_UNDEFINED) fl_out_fmt=NC_FORMAT_CLASSIC;
  if(fl_out_fmt == NC_FORMAT_CLASSIC) nccreate_mode+=0; // CEWI
#endif // !ENABLE_NETCDF4
  const int nc_out(nco_create(fl_out,nccreate_mode)); 
  const nc_type nco_xtyp(nco_get_xtype(PRC_CMP(1.0))); // [enm] External netCDF type
  if(dbg_lvl > dbg_off || tst_sng == "nc"){
    std::cout << prg_nm << ": INFO External netCDF type of prc_cmp variables will be " << nco_typ_sng(nco_xtyp) << std::endl;
    std::cout << prg_nm << ": INFO Record dimension will be " << (dmn_rcd == "" ? "non-existent" : dmn_rcd) << std::endl;
  } // !dbg
 
  // Create dimensions
  const int wvl_dmn(nco_def_dim(nc_out,static_cast<std::string>("wvl"),dmn_rcd == "wvl" ? NC_UNLIMITED : wvl_nbr)); // [dmn] Wavelength dimension
  const int sz_dmn(nco_def_dim(nc_out,static_cast<std::string>("sz"),dmn_rcd == "sz" ? NC_UNLIMITED : sz_nbr)); // [dmn] Size dimension
  const int bnd_dmn(nco_def_dim(nc_out,static_cast<std::string>("bnd"),dmn_rcd == "bnd" ? NC_UNLIMITED : bnd_nbr)); // [dmn] Band dimension
  const int lat_dmn(nco_def_dim(nc_out,lat_nm,dmn_rcd == lat_nm ? NC_UNLIMITED : lat_nbr)); // [dmn] Latitude dimension
  const int lon_dmn(nco_def_dim(nc_out,lon_nm,dmn_rcd == lon_nm ? NC_UNLIMITED : lon_nbr)); // [dmn] Longitude dimension
  // Derived dimensions
  const int dmn_bnd_sz[2]={bnd_dmn,sz_dmn}; // [dmn] 
  const int dmn_lat_lon[2]={lat_dmn,lon_dmn}; // [dmn]
  const int *dmn_lat(&lat_dmn); // [dmn] Pointer to latitude dimension
  const int *dmn_lon(&lon_dmn); // [dmn] Pointer to longitude dimension
  const int *dmn_sz(&sz_dmn); // [dmn] Pointer to size dimension
  const int *dmn_wvl(&wvl_dmn); // [dmn] Pointer to wavelength dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars CLIP

  if(dbg_lvl == dbg_old || tst_sng == "nc"){
    rcd=nco_inq_ndims(nc_out,int_foo);
    std::cerr << "Currently there are " << int_foo << " dimensions defined" << std::endl;
  } // !dbg
  // Global attributes
  rcd=nco_put_att(nc_out,NC_GLOBAL,"history",time_bfr_srt.substr(0,time_bfr_srt.size()-1)+": "+cmd_ln);
  rcd=nco_put_att(nc_out,NC_GLOBAL,"CVS_Id",CVS_Id);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("omp_thread_number"),thr_nbr);

  var_mtd_sct var_mtd[]={
    {0,"tpt_v1d",nco_xtyp,2,dmn_bnd_sz,"long_name","Temperature (vector 1-dimensional)","units","kelvin"},
    {0,"tpt_a2d",nco_xtyp,2,dmn_bnd_sz,"long_name","Temperature (array 2-dimensional)","units","kelvin"},
    {0,"wvl",nco_xtyp,1,dmn_wvl,"long_name","Wavelength at band center","units","meter"},
    {0,"sz",nco_xtyp,1,dmn_sz,"long_name","Size at bin center","units","meter"},
    {0,"sz_nbr",NC_INT,0,dmn_scl,"long_name","Number of sizes","units","number"},
    {0,"float_var",NC_FLOAT,0,dmn_scl,"long_name","float test variable","units","number"},
    {0,"double_var",NC_DOUBLE,0,dmn_scl,"long_name","double test variable","units","number"},
    {0,"int_var",NC_INT,0,dmn_scl,"long_name","int test variable","units","number"},
    {0,"short_var",NC_SHORT,0,dmn_scl,"long_name","short test variable","units","number"},
    {0,"char_var",NC_CHAR,0,dmn_scl,"long_name","char test variable","units","number"},
    {0,"byte_var",NC_BYTE,0,dmn_scl,"long_name","byte test variable","units","number"},
#ifdef ENABLE_NETCDF4
    {0,"ubyte_var",NC_UBYTE,0,dmn_scl,"long_name","ubyte test variable","units","number"},
    {0,"ushort_var",NC_USHORT,0,dmn_scl,"long_name","ushort test variable","units","number"},
    {0,"uint_var",NC_UINT,0,dmn_scl,"long_name","uint test variable","units","number"},
    {0,"int64_var",NC_INT64,0,dmn_scl,"long_name","int64 test variable","units","number"},
    {0,"uint64_var",NC_UINT64,0,dmn_scl,"long_name","uint64 test variable","units","number"},
    //    {0,"string_var",NC_STRING,0,dmn_scl,"long_name","string test variable","units","number"},
#endif // !ENABLE_NETCDF4
    // Map grid variables
    {0,edg_est_nm.c_str(),nco_xtyp,0,dmn_scl,"long_name","eastern edge of surface grid","units","degrees east"},
    {0,edg_nrt_nm.c_str(),nco_xtyp,0,dmn_scl,"long_name","northern edge of surface grid","units","degrees north"},
    {0,edg_sth_nm.c_str(),nco_xtyp,0,dmn_scl,"long_name","southern edge of surface grid","units","degrees north"},
    {0,edg_wst_nm.c_str(),nco_xtyp,0,dmn_scl,"long_name","western edge of surface grid","units","degrees east"},
    {0,lat_ctr_2d_nm.c_str(),nco_xtyp,2,dmn_lat_lon,"long_name","latitude-2d","units","degrees north"},
    {0,lat_nm.c_str(),nco_xtyp,1,dmn_lat,"long_name","Latitude","units","degrees north"},
    {0,lon_ctr_2d_nm.c_str(),nco_xtyp,2,dmn_lat_lon,"long_name","longitude-2d","units","degrees east"},
    {0,lon_nm.c_str(),nco_xtyp,1,dmn_lon,"long_name","Longitude","units","degrees east"}
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct));

  sng2var_mtd_map var_mtd_map;
  for(idx=0;idx<var_mtd_nbr;idx++){
    /* fxm: Define variables before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map. */
    rcd=nco_def_var(nc_out,var_mtd[idx].nm,var_mtd[idx].type,var_mtd[idx].dmn_nbr,var_mtd[idx].dmn_id,var_mtd[idx].var_id);
    var_mtd_map.insert(sng2var_mtd_map::value_type(var_mtd[idx].nm,var_mtd[idx]));
  } // end loop over itr

  sng2var_mtd_map::const_iterator var_mtd_itr;
  for(var_mtd_itr=var_mtd_map.begin();var_mtd_itr!=var_mtd_map.end();++var_mtd_itr){
    // Write first attribute (long_name)
    rcd=nco_put_att(nc_out,var_mtd_itr->second.var_id,var_mtd_itr->second.att_1_nm,var_mtd_itr->second.att_1_val);
    // Write second attribute (units)
    rcd=nco_put_att(nc_out,var_mtd_itr->second.var_id,var_mtd_itr->second.att_2_nm,var_mtd_itr->second.att_2_val);
    std::cout << "Defined " << var_mtd_itr->first << " with long_name = " << var_mtd_itr->second.att_1_val << " and units = " << var_mtd_itr->second.att_2_val << std::endl;
  } // end loop over var_mtd_itr

  std::cerr << var_mtd_map.find("wvl")->second.att_1_val << std::endl;

  // Leave define mode
  rcd=nco_enddef(nc_out,NC_ENOTINDEFINE); // [fnc] Leave define mode

  /* Write at least one variable containing record dimension, if any,
     with nco_put_vara() routine to make output file aware of record dimension size.
     Or else calls to nco_put_var() routines assume record dimension is size zero 
     Only generic way to do this is to write all coordinates with this formalism */
  rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("sz"),sz_nbr,sz); delete []sz;
  rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("wvl"),wvl_nbr,wvl); delete []wvl;

  // Write data and delete dynamic arrays
  // Syntax valid for statically allocated arrays foo[dim1][dim2]
  // rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_d2d"),&tpt_d2d[0]);
  // Syntax valid for vector< vector<prc_cmp> > foo(dim1,vector<prc_cmp>(dim2))
  // rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_v2d"),&tpt_v2d[0][0]);
  // Syntax valid for vector<prc_cmp> foo(dim1*dim2)
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_v1d"),&tpt_v1d[0]);
  // Syntax valid for a2d_cls<prc_cmp> foo(dim1,dim2)
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_a2d"),&tpt_a2d(0,0));
  rcd=nco_put_var(nc_out,static_cast<std::string>("sz_nbr"),sz_nbr);
  // netCDF3 types
  float float_var(FLT_MAX);
  double double_var(DBL_MAX);
  int int_var(INT_MAX);
  short short_var(SHRT_MAX);
  char char_var(CHAR_MAX);
  signed char byte_var(CHAR_MAX);
  // netCDF4 types
#ifdef ENABLE_NETCDF4
  unsigned char ubyte_var(UCHAR_MAX);
  unsigned short ushort_var(USHRT_MAX);
  unsigned int uint_var(UINT_MAX);
#if (defined LLONG_MAX)
  long long int64_var(LLONG_MAX);
#else // !LLONG_MAX
  long long int64_var(0LL);
#endif // !LLONG_MAX
#if (defined ULLONG_MAX)
  unsigned long long uint64_var(ULLONG_MAX);
#else // !ULLONG_MAX
  unsigned long long uint64_var(0ULL);
#endif // !ULLONG_MAX
  std::string string_var("flamenco uber alles");
#endif // !ENABLE_NETCDF4
  rcd=nco_put_var(nc_out,static_cast<std::string>("float_var"),float_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("double_var"),double_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("int_var"),int_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("short_var"),short_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("char_var"),char_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("byte_var"),byte_var);
#ifdef ENABLE_NETCDF4
  rcd=nco_put_var(nc_out,static_cast<std::string>("ubyte_var"),ubyte_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ushort_var"),ushort_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("uint_var"),uint_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("int64_var"),int64_var);
  rcd=nco_put_var(nc_out,static_cast<std::string>("uint64_var"),uint64_var);
  // fxm: Outputting a single string may not work until netCDF4-beta supports more get/put var?_string() functions
  //  rcd=nco_put_var(nc_out,static_cast<std::string>("string_var"),string_var.c_str());
#endif // !ENABLE_NETCDF4

  // Map grid
  rcd=nco_put_vara_crd(nc_out,lat_nm,lat_nbr,lat); delete []lat;
  rcd=nco_put_vara_crd(nc_out,lon_nm,lon_nbr,lon); delete []lon;
  rcd=nco_put_var(nc_out,edg_nrt_nm,edg_grd[0]);
  rcd=nco_put_var(nc_out,edg_est_nm,edg_grd[1]);
  rcd=nco_put_var(nc_out,edg_sth_nm,edg_grd[2]);
  rcd=nco_put_var(nc_out,edg_wst_nm,edg_grd[3]);
  rcd=nco_put_var(nc_out,lat_ctr_2d_nm,&lat_ctr_2d(0,0));
  rcd=nco_put_var(nc_out,lon_ctr_2d_nm,&lon_ctr_2d(0,0));

  if(dbg_lvl == dbg_old || tst_sng == "rt"){
    std::cout << "Testing rt_cls..." << std::endl;
    // Instantiate radiative transfer object
    rt_cls rt_obj(static_cast<std::string>("two_srm_iso_sct")); // [rt] Radiative transfer object
    // Print radiative transfer object
    std::cout << rt_obj; // [rt] Radiative transfer object
    // Test radiative transfer object for memory leaks by initializing lng_foo members
    rcd+=rt_cls::tst(lng_foo); // [rt] Radiative transfer 
    // Test output functions
    //    rcd+=rt_obj.var_put(nc_out,"rfl_flx"); // [frc] Flux reflectance
    const prc_cmp tmp_foo(0.5); // [frc] Cosine solar zenith angle
    rcd+=rt_lop(nc_out,rt_obj,tmp_foo); // [frc] Flux reflectance
  } // !dbg || "rt"

  if(dbg_lvl == dbg_old || tst_sng == "wvl"){
    // ccc --tst=wvl
    std::cout << "Testing wvl_grd_cls..." << std::endl;
    // Instantiate wavelength grid object
    prc_cmp wvl_mnm_arg(-1); // I/O [m] Minimum wavelength
    prc_cmp wvl_mxm_arg(-1); // I/O [m] Maximum wavelength
    long wvl_nbr_arg(-1); // I/O [nbr] Number of wavelengths
    wvl_grd_cls wvl_obj(static_cast<std::string>("SWNB"),wvl_mnm_arg,wvl_mxm_arg,wvl_nbr_arg); // [obj] Wavelength grid
    // Print wavelength grid object
    std::cout << wvl_obj; // [obj] Wavelength grid object
    // Test wavelength grid object for memory leaks by initializing lng_foo members
    //rcd+=wvl_grd_cls::tst(lng_foo); // [obj] Wavelength grid 
  } // !dbg || "wvl"

  // Close output file
  rcd=nco_close(nc_out); // [fnc] Close netCDF file
  std::cerr << "Wrote results to " << fl_out << std::endl;
  std::cerr << "ncks: ncks -C -F -d wvl,0.5e-6 -v wvl " << fl_out << std::endl;
  std::cerr << "ncks: ncks --no_crd --prn --ftn --mtd --units --dmn wvl,0.5e-6 --var wvl " << fl_out << std::endl;

  Exit_gracefully();
  return EXIT_SUCCESS;
} // !main()

void usg_prn(const char *opt_sng)
{
  // Purpose: Print correct usage of program
  std::cerr << "\nusage: " << prg_nm << " [-options] where options are one or more of:\n\n";
  std::cerr << "--dbl_foo dbl_foo Set generic double (default is 0)" << std::endl;
  std::cerr << "--dbg,-D dbg_lvl Set debugging level (default is 0)" << std::endl;
  std::cerr << "--flt_foo,-f flt_foo Set generic float (default is 0)" << std::endl;
  std::cerr << "--cpv_foo,-r cpv_foo Set generic computational precision (default is 0)" << std::endl;
  std::cerr << "-i fl_in Set input file name. (default is in.nc)" << std::endl;
  std::cerr << "-o fl_out Set output file name. Default is out.nc" << std::endl;
  std::cerr << "-v print the CVS program version" << std::endl;
  std::cerr << "--tst_sng tst_typ Perform system/compiler/software for test type:" << std::endl;
  std::cerr << "\ta2d\tArrays of two dimensions" << std::endl;
  std::cerr << "\taer\tAerosol class" << std::endl;
  std::cerr << "\tbnr\tBinary representation" << std::endl;
  std::cerr << "\tc99\tC99 syntax features" << std::endl;
  std::cerr << "\tcmp\tCompiler" << std::endl;
  std::cerr << "\tcpp\tC preprocessor" << std::endl;
  std::cerr << "\tcpx\tComplex numbers" << std::endl;
  std::cerr << "\tcrp\tCryptography" << std::endl;
  std::cerr << "\tdyn\tDynamic allocation" << std::endl;
  std::cerr << "\tflt\tFloating point representation <cfloat>" << std::endl;
  std::cerr << "\tfnc\tFunction instantiation" << std::endl;
  std::cerr << "\tgmt\tTime" << std::endl;
  std::cerr << "\tgsl\tGSL (GNU Scientific Library) functions" << std::endl;
  std::cerr << "\tidx_rfr\tRefractive index class" << std::endl;
  std::cerr << "\tld\tLinker data" << std::endl;
  std::cerr << "\tlmt, ntg\tLimits, integer representation <climits>" << std::endl;
  std::cerr << "\tmap\tMap class" << std::endl;
  std::cerr << "\tmmr\tMemory and malloc() handling" << std::endl;
  std::cerr << "\tmth\tMathematical utilities, constants" << std::endl;
  std::cerr << "\tnc\tnetCDF" << std::endl;
  std::cerr << "\tnms\tNamespaces" << std::endl;
  std::cerr << "\tntg, lmt\tInteger representation, limits <climits>" << std::endl;
  std::cerr << "\tntp\tVector interpolation and rebinning" << std::endl;
  std::cerr << "\tnvr\tEnvironment variables" << std::endl;
  std::cerr << "\tnwt_rph\tNewton-Raphson root finding" << std::endl;
  std::cerr << "\tocn\tDrag coefficient over ocean" << std::endl;
  std::cerr << "\tomp\tOpenMP" << std::endl;
  std::cerr << "\toom\tOut of Memory behavior" << std::endl;
  std::cerr << "\tphys_cst\tPhysical Constants" << std::endl;
  std::cerr << "\tptr\tPointers" << std::endl;
  std::cerr << "\trbn\tVector interpolation and rebinning" << std::endl;
  std::cerr << "\trnd\tround() function" << std::endl;
  std::cerr << "\trt\tRadiative transfer class" << std::endl;
  std::cerr << "\trtti\tRun time type instantiation" << std::endl;
  std::cerr << "\tsizeof,sz\tsizeof() operations and results" << std::endl;
  std::cerr << "\tsng,srm\tStrings and streams" << std::endl;
  std::cerr << "\tspc_slr\tSolar spectrum class" << std::endl;
  std::cerr << "\tsys\tSystem functions" << std::endl;
  std::cerr << "\ttry\tTry/throw/catch exceptions" << std::endl;
  std::cerr << "\ttst\tTest" << std::endl;
  std::cerr << "\tvec\tVector interpolation and rebinning" << std::endl;
  std::cerr << "\twvl\tWavelength grid class" << std::endl;
  std::cerr << std::endl;

  if(&opt_sng == &opt_sng){;} // CEWU Compiler Error Warning Usage
} // end usg_prn() 

void virtualViaReference(const SzDstFnc &fnc,const prc_cmp rds)
{
  //  fnc(rds);
  if(&fnc == &fnc){;} // CEWU Compiler Error Warning Usage
  if(&rds == &rds){;} // CEWU Compiler Error Warning Usage
} // end virtualViaReference

char ** /* O [sng] Array of list elements */
lst_prs /* [fnc] Create list of strings from given string and arbitrary delimiter */
(char *sng_in, /* I/O [sng] Delimited argument list (delimiters are changed to NULL on output */
 const char *dlm_sng, /* I [sng] delimiter string */
 int *nbr_lst) /* O [nbr] number of elements in list */
{
  /* Purpose: Create list of strings from given string and arbitrary delimiter
     This routine is often called with system memory, e.g., with strings from
     command line arguments whose memory was allocated by the shell or getopt().
     A conservative policy would be, therefore, to never modify the input string
     However, we are safe if any modifications do not extend the input string
     Thus this routine is allowed to replace delimiter strings by NULs */

  /* Number of list members is always one more than number of delimiters, e.g.,
     foo,,3, has 4 arguments: "foo", "", "3" and "".
     A delimiter without an argument is valid syntax to indicate default argument
     Therefore a storage convention is necessary to indicate default argument was selected
     Either NULL or '\0' can be used without requiring additional flag
     NULL is not printable, but is useful as a logical flag since its value is False
     On the other hand, '\0', the empty string, can be printed but is not as useful as a flag
     Currently, NCO implements the former convention, where default selections are set to NULL */
  
  char *sng_in_ptr;

  int idx;
  size_t dlm_lng;

  /* Delimiter must be NUL-terminated (a string) so we may find its length */
  dlm_lng=strlen(dlm_sng); 

  /* Do not increment actual sng_in pointer while searching for delimiters---increment a dummy pointer instead. */
  sng_in_ptr=sng_in; 

  /* First element does not require a delimiter in front of it */
  *nbr_lst=1;

  /* Count list members */
  while((sng_in_ptr=strstr(sng_in_ptr,dlm_sng))){
    sng_in_ptr+=dlm_lng;
    (*nbr_lst)++;
  } /* end while */

  // lst=(char **)std::malloc(*nbr_lst*sizeof(char *));
  // 20030928 Hard to believe this funky syntax is legal
  char **lst=new char *[static_cast<size_t>(*nbr_lst)];

  sng_in_ptr=sng_in; 
  lst[0]=sng_in;
  idx=0;
  while((sng_in_ptr=strstr(sng_in_ptr,dlm_sng))){
    /* NUL-terminate previous arg */
    *sng_in_ptr='\0';
    sng_in_ptr+=dlm_lng;
    lst[++idx]=sng_in_ptr;
  } /* end while */

  /* A default list member is assumed whenever two delimiters are adjacent to eachother, such that
     the length of the string between them is 0. If the list ends with a delimiter, then the last
     element of the list is also assumed to be a default list member. */
  /* This loop sets default list members to NULL */
  for(idx=0;idx<*nbr_lst;idx++)
    if(strlen(lst[idx]) == 0) lst[idx]=(char *)NULL;

  if(dbg_lvl == 5){
    (void)std::fprintf(stderr,"%d elements in list delimited by \"%s\"\n",*nbr_lst,dlm_sng);
    for(idx=0;idx<*nbr_lst;idx++) 
      (void)std::fprintf(stderr,"lst[%d] = %s\n",idx,(lst[idx] == NULL) ? "NULL" : lst[idx]);
    (void)std::fprintf(stderr,"\n");
    (void)fflush(stderr);
  } /* end debug */

  return lst;
} /* end lst_prs() */

void 
my_new_handler(void) // [fnc] Handler for bad_alloc()'s thrown by new[]
{
  /* Purpose: Generic handler for bad_alloc()'s thrown by new[]
     Usage: Register in main() with set_new_handler()
     DeD01 p. 744 
     NB: AIX xlC cannot optimize code with custom new_handler()
     Hence, my_new_handler() is in main.c rather than library */
  const std::string sbr_nm("my_new_handler"); // [sng] Subroutine name
  std::cerr << prg_nm_get() << ": ERROR Memory allocation with new[] failed" << std::endl;
  std::cerr << prg_nm_get() << ": Homebrew new[] handler "+sbr_nm+"() has been called" << std::endl;
  std::cerr << prg_nm_get() << ": HINT Possible Out-of-Memory (OOM) condition, try killing non-vital processes and running again" << std::endl;
  std::cerr << prg_nm_get() << ": Calling std::abort() (i.e., dumping core)" << std::endl;
  std::abort(); // [fnc] Dump core file
  // std::exit(); // [fnc] Exit without core dump
} // end my_new_handler()
