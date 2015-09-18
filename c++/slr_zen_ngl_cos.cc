// $Id$

// Purpose: Compute and return solar zenith angle for given coordinates

/* Compilation:
   cd ~/c++;make slr_zen_ngl_cos;cd -
   g++ -O2 -g -Wall -Werror -o ${MY_BIN_DIR}/slr_zen_ngl_cos ~/c++/slr_zen_ngl_cos.cc
   xlC_r -O -q64 -o ${MY_BIN_DIR}/slr_zen_ngl_cos ${MY_OBJ_DIR}/getopt_bsd.o ~/c++/slr_zen_ngl_cos.cc */

/* Distribution:
   scp ~/c++/slr_zen_ngl_cos.cc dust.ess.uci.edu:c++
   scp ~/c++/slr_zen_ngl_cos.cc esmf.ess.uci.edu:c++ */

/* Usage: 
   slr_zen_ngl_cos
   slr_zen_ngl_cos --lat_dgr=+23.44 --doy=172.5
   slr_zen_ngl_cos --lat_dgr=-23.44 --doy=355.55

   Production: */

/* History: 
   20080901 Initial version */

// Standard C++ headers
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr
#include <string> // Standard C++ string class

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159
#include <cstdio> // stderr, EOF, FILE, NULL, etc.
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv

// Prototypes
extern "C" {
#include "getopt.h" // GNU getopt() functionality
} // end extern

int main(int argc,char **argv)
{
  // Functions
  int // O [enm] Return success code
    slr_crd_Bri92 // [fnc] Compute solar zenith angle and eccentricity
    (const double lat, // I [rdn] Latitude
     const double lcl_yr_day, // I [day] Local year day
     double &slr_zen_ngl_cos, // O [frc] Solar zenith angle cosine
     double &xnt_fac); // O [frc] Eccentricity factor
  // end slr_crd_Bri92() prototype

  // Parameters
  // Physical and numerical constants
  const double CEWI_dbl(9.9692099683868690e+36); // Compiler Error Warning Initializer for double
  const double cst_M_PIl(3.1415926535897932384626433832795029L); // [frc] 3

  // Set defaults for command line options 
  // Option name is variable name, e.g., --lng_foo=3, unless otherwise indicated
  bool flg_flg(false); // [flg] Temporary flag variable
  bool flg_txt(false); // [flg] Print day/night answer in text
  double dbl_foo(0.0); // [frc] Intrinsic double temporary variable
  double lcl_yr_day(172.5); // [day] Local year day
  double lat_dgr(23.44); // [dgr] Latitude
  double lat(CEWI_dbl); // [rdn] Latitude
  double lat_rdn(CEWI_dbl); // [rdn] Latitude
  float flt_foo(0.0f); // [frc] Intrinsic float temporary variable
  int int_foo(0); // [nbr] Intrinsic int temporary variable
  long lng_foo(0L); // [nbr] Intrinsic long temporary variable
  std::string sng_foo(""); // [sng] Intrinsic string temporary variable
  unsigned short dbg_lvl(0); // [enm] Debugging level

  // Derived constants

  // Doubly-derived constants

  // Trebly-derived variables

  // Locals requiring initialization 
  int rcd(0); // [enm] Return success code

  // Locals
  double slr_zen_ngl_cos; // O [frc] Solar zenith angle cosine
  double xnt_fac; // O [frc] Eccentricity factor

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}
       If flag is non-zero, getopt_long() returns zero and flag is set to val
       If flag is zero, getopt_long() returns contents of val */
    // Long options with no argument, no short option counterpart
    {"flg_flg",no_argument,0,0}, // [flg] Flag flag
    {"flg_txt",no_argument,0,0}, // [flg] Print day/night answer in text
    // Long options with argument, no short option counterpart
    {"dbl_foo",required_argument,0,0}, // [nbr] Intrinsic double temporary variable 
    {"doy",required_argument,0,0}, // [day] Local year day
    {"lcl_yr_day",required_argument,0,0}, // [day] Local year day
    {"lat_rdn",required_argument,0,0}, // [rdn] Latitude
    {"lat_dgr",required_argument,0,0}, // [dgr] Latitude
    {"int_foo",required_argument,0,0}, // [nbr] Intrinsic int temporary variable
    {"lng_foo",required_argument,0,0}, // [nbr] Intrinsic long temporary variable 
    {"sng_foo",required_argument,0,0}, // [sng] Intrinsic string temporary variable
    // Long options with optional argument, no short option counterpart
    // Long options with short counterparts
    {"dbg_lvl",optional_argument,0,'D'}, // [enm] Debugging level
    {"flt_foo",required_argument,0,'f'}, // [frc] Intrinsic float temporary variable
    {"help",no_argument,0,'h'},
    {"version",no_argument,0,'v'},
    // Last option named "0" signals getopt_long() to stop processing  
    {0,0,0,0}
  }; // end opt_lng
  
  // Short options: no colon = no arg, one colon = required arg, two colons = optional arg
  const char * const opt_sht_lst("D:f:"); // [sng] List of single-letter (C-style) option abbreviations
  extern char *optarg; // [sng] char * representation of current optarg, if any (this memory is owned by system)
  // extern int optind; // [idx] extern enumerating cardinal of current option
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
      if(dbg_lvl >= 5) std::cerr << "Long option name: " << opt_crr << (optarg ? ",  Argument: "+opt_sng : ", No Argument") << std::endl;
      if(opt_crr == "dbl_foo") dbl_foo=std::strtod(opt_sng.c_str(),(char **)NULL); // [nbr] Intrinsic double temporary variable
      if(opt_crr == "doy" || opt_crr == "lcl_yr_day") lcl_yr_day=std::strtod(opt_sng.c_str(),(char **)NULL); // [day] Local year day
      if(opt_crr == "lat_dgr"){
	lat_dgr=std::strtod(opt_sng.c_str(),(char **)NULL); // [dgr] Latitude
	lat=lat_rdn=cst_M_PIl*lat_dgr/180.0; // [rdn] Latitude
      } // endif lat_dgr
      if(opt_crr == "lat_rdn"){
	lat=lat_rdn=std::strtod(opt_sng.c_str(),(char **)NULL); // [rdn] Latitude
	lat_dgr=180.0*lat_rdn/cst_M_PIl; // [dgr] Latitude
      } // endif lat_rdn
      if(opt_crr == "flg_flg") flg_flg=true;
      if(opt_crr == "flg_txt") flg_txt=true; // [flg] Print day/night answer in text
      if(opt_crr == "int_foo") int_foo=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Intrinsic int temporary variable
      if(opt_crr == "lng_foo") lng_foo=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Intrinsic long temporary variable
      if(opt_crr == "sng_foo") sng_foo=opt_sng; // [sng] Intrinsic string temporary variable
      // Multi-line initializations break alphabetical-order-by-line rule
    } // opt != 0
    switch(opt){
    case 0: // Long options have already been processed, return
      break;
    case 'D': // Debugging level (default is 0) 
      if(optarg) dbg_lvl=static_cast<unsigned short int>(std::strtoul(opt_sng.c_str(),(char **)NULL,10)); else dbg_lvl=1;
      break;
    case 'f': // Set generic tuning parameter (default is 0.0)
      flt_foo=static_cast<float>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Intrinsic float temporary variable
      break;
    default:
      return EXIT_FAILURE;
    } // end opt switch 
  } // end while loop 

  if(argc <= 1){
    std::cout << "ERROR: need latitude and/or day-of-year argument" << std::endl;
    std::cout << "Usage: slr_zen_ngl_cos [--lat_dgr] [--doy]" << std::endl;
    std::cout << "e.g., slr_zen_ngl_cos --lat_dgr=+23.44 --doy=172.5" << std::endl;
    std::cout << "Exiting...try again" << std::endl;
    return EXIT_FAILURE;
  } // endif argc
  
  // Derived variables that depend on command-line input
  if(lat==CEWI_dbl) lat=lat_rdn=cst_M_PIl*lat_dgr/180.0; // [rdn] Latitude

  // Doubly-derived variables that depend on command-line input

  // Finished with inputs and derived variables
  // Validate geometric configuration

  if(dbg_lvl == 1){
    std::cout << "Initialization State:" << std::endl;
    std::cout << "lcl_yr_day = " << lcl_yr_day << std::endl;
    std::cout << "lat = " << lat << std::endl;
  } // endif dbg
  
  // Compute solar zenith angle and eccentricity
  rcd=slr_crd_Bri92
    (lat, // I [rdn] Latitude
     lcl_yr_day, // I [day] Local year day
     slr_zen_ngl_cos, // O [frc] Solar zenith angle cosine
     xnt_fac); // O [frc] Eccentricity factor

  if(flg_txt){
    std::cout << (slr_zen_ngl_cos > 0.0 ? "day" : "night") << std::endl;
  }else{
    std::cout << slr_zen_ngl_cos << std::endl;
  } // !flg_txt  

} // end snw

int // O [enm] Return success code
slr_crd_Bri92 // [fnc] Compute solar zenith angle and eccentricity
(const double lat, // I [rdn] Latitude
 const double lcl_yr_day, // I [day] Local year day
 double &slr_zen_ngl_cos, // O [frc] Solar zenith angle cosine
 double &xnt_fac) // O [frc] Eccentricity factor
{
  /* Purpose: Compute solar geometry
     Reference: B. P. Briegleb's routine used in CCM2, CCM3 */

  int rcd(0); // [enm] Return success code
  const double days_per_year(365.0); // [day]
  double delta; // [rdn] Solar declination
  double phi; // [rdn] Local phase angle (0 is midnight)
  double theta; // [rdn] Solar polar coordinate angle
  const double cst_M_PIl(3.1415926535897932384626433832795029L); // [frc] 3
  
  // Compute eccentricity factor (Sun-Earth distance factor)
  theta=2.0*cst_M_PIl*lcl_yr_day/days_per_year;
  xnt_fac=1.000110+0.034221*cos(theta)+0.001280*sin(theta)+ 
    0.000719*cos(2.0*theta)+0.000077*sin(2.0*theta);
    
  // Solar declination in radians
  delta=0.006918-0.399912*cos(theta)+0.070257*sin(theta)- 
    0.006758*cos(2.0*theta)+0.000907*sin(2.0*theta)- 
    0.002697*cos(3.0*theta)+0.001480*sin(3.0*theta);
    
  phi=2.0*cst_M_PIl*lcl_yr_day;
  slr_zen_ngl_cos=
    sin(lat)*sin(delta)-
    cos(lat)*cos(delta)*cos(phi);

  return rcd; // [enm] Return success code
} // end slr_crd_Bri92()
