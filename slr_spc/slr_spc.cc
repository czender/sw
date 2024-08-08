// $Id$ 

// Purpose: Compute and store TOA solar spectrum properties

/* ThD71 binned their data every 0.01 microns from 0.2 up to 1.0 micron 
   then by 0.2 microns up to 5.0 microns for a total of 120 intervals. 
   Valid wavelengths are between 0.15 and 4.95 microns inclusive.
  
   LaN68 binned their data every 0.01 microns from 0.2 up to 0.6 microns 
   then by 0.1 microns up to 5 microns for a total of 120 intervals. 
   Valid wavelengths are between 0.2 and 5.0 microns inclusive.
  
   Users of various solar spectral flux datasets:
   LaN68: Bri92, CCM2, CCM3, CCM4
   ThD71: TaL89, Ste78, Sli89
   NeL84: EbC92 
   Kur95: ZBP97
   JHC21: JHC21 */

/* Strategy: 
   1. Front end section: acquire data from measurement source in 
   whatever format and units measurements are supplied.
   Convert input to standard fields.
   Standard fields from which all else is computed are:
   wvl // [m] Nominal wavelength
   wvl_dlt // [m] Bandwidth
   wvl_min // [m] Minimum wavelength in band
   wvl_max // [m] Maximum wavelength in band
   wvl_ctr // [m] Wavelength at band center
   wvl_nbr // [nbr] Number of wavelength bands
   flx_bnd // [W m-2] Solar flux in band
   2. Create standard derived fields
   3. Back end: Write out standard derived fields */

/* Example usage:
   
   Build (as of 20240808) on spectral.ess.uci.edu:
   make CPPFLAGS="-DLINUX -I ${HOME}/include -I/opt/netcdf/include -I/opt/homebrew/include" LDFLAGS="-L/opt/netcdf/lib -L/opt/homebrew/lib -L${HOME}/lib -lcsz_c++ -lnco_c++ -lnetcdf"

   Debugging:
   slr_spc --dbg=1 --wvl_nbr=2497
   slr_spc --dbg=1 --wvl_nbr=2497 ${DATA}/aca/spc_Kur95_20wvn.txt ${DATA}/aca/spc_Kur95_20wvn.nc
   slr_spc --dbg=1 --wvl_lmt=5.0 --wvl_nbr=2497 ${DATA}/aca/spc_Kur95_20wvn.txt ${DATA}/aca/spc_Kur95_20wvn.nc

   Production:
   slr_spc --dbg=1 --ncr_wvl ${DATA}/aca/tsis_ssi_L3_c24h_20240807.nc ${DATA}/aca/spc_JHC21.nc
   slr_spc --dbg=1 --ncr_wvl --wvl_nbr=49934 ${DATA}/aca/spc_Kur95_01wvn.txt ${DATA}/aca/spc_Kur95_01wvn.nc
   slr_spc --dbg=1 --ncr_wvl --wvl_nbr=2497 ${DATA}/aca/spc_Kur95_20wvn.txt ${DATA}/aca/spc_Kur95_20wvn.nc
   slr_spc --dbg=1 --flx_slr_src=NeL84 ${DATA}/aca/spc_NeL84.txt ${DATA}/aca/spc_NeL84.nc
   slr_spc --dbg=1 --flx_slr_src=LaN68 -o ${DATA}/aca/spc_LaN68.nc
   slr_spc --dbg=1 --flx_slr_src=ThD71 -o ${DATA}/aca/spc_ThD71.nc

   Validation:

   */

// Standard C++ header files 
#include <iostream> // Standard C++ I/O streams cout, cin
#include <iomanip> // setw()
#include <string> // Standard C++ string class
 
// Standard C header files 
#include <cstdio> // stderr, FILE, NULL, etc.
#include <ctime> // Machine time
#include <cstring> // strcmp()...
#include <cmath> // sin cos cos sin 3.14159
#include <cstdlib> // strtod, strtol, malloc, getopt 

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// #define MAIN_PROGRAM_FILE MUST precede #include slr_spc.hh
#define MAIN_PROGRAM_FILE
#include "slr_spc.hh" // Program-specific definitions
#include <libcsz_c++.hh> // Personal C++ library
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Namespace
//using namespace std; // std is default namespace 

// Prototypes
extern "C" {
#include "getopt.h" // GNU getopt() functionality from BSD my_getopt()
} // end extern

int main(const int argc,char **argv)
{
  // Prototypes  
  void usg_prn(const char *opt_sng);

  // Enums
  enum flx_slr_src{LaN68,ThD71,NeL84,Kur95,JHC21};
  const std::string flx_slr_sng[5]={
    "Labs and Neckel (1968) (LaN68)",
    "Thekeakara and Drummond (1971) (ThD71)",
    "Neckel and Labs (via WDC) (1984) (NeL84)",
    "Kurucz (1995) (Kur95)",
    "Jing et al. (2021)"
  };
  int flx_slr_src=JHC21; // default

  // Locals
  FILE *fp_in; // instead of cin, for now
  char bff[1000]; // Space for scanf() input
  float (*flx_frc_fnc)(const float *wvl_min_mcr,const float *wvl_max_mcr);
  float flt_foo2;
  float flt_foo;
  prc_cmp wvl_ctr_mcr;
  int int_foo;
  int rcd(0);
  long bnd_idx;
  long idx;
  long idx_max_wvl; // [idx] Index of maximum wavelenth
  long idx_min_wvl; // [idx] Index of minimum wavelenth
  long ncr_ncr_wvl; // [nbr] Increment to increasing wavelength
  long wvl_idx;
  const std::string CVS_Date("$Date$"); // [sng] CVS date string
  const std::string CVS_Header("$Id$"); // [sng] CVS header string
  const std::string CVS_Id("$Id$"); // [sng] CVS identification string
  const std::string CVS_Revision("$Revision$"); // [sng] CVS revision string
  const std::string date_cvs(CVS_Date.length() > 7 ? CVS_Date.substr(7,19) : "Unknown"); // [sng] Date from CVS
  const std::string sbr_nm("main"); // [sng] Subroutine name
  const std::string vrs_cvs(CVS_Revision.length() > 10 ? CVS_Revision.substr(10,4) : "Unknown"); // [sng] Version from CVS

  // C pre-processor macros for instantiating variable values with string tokens
#define XTKN2SNG(x) #x
#define TKN2SNG(x) XTKN2SNG(x)
  const std::string date_cpp(__DATE__); // [sng] Date from C pre-processor
  const std::string time_cpp(__TIME__); // [sng] Time from C pre-processor
  const std::string vrs_cpp(TKN2SNG(VERSION)); // [sng] Version from C pre-processor
  const std::string hst_cpp(TKN2SNG(HOSTNAME)); // [sng] Hostname from C pre-processor
  const std::string usr_cpp(TKN2SNG(USER)); // [sng] Hostname from C pre-processor

  // Start clock and save command line  
  const std::time_t time_crr_time_t(std::time((std::time_t *)NULL)); // [tm] Current date and time
  const std::string time_bfr_srt(std::ctime(&time_crr_time_t)); // [sng] Current date and time
  std::cerr << "\tStart = " << time_bfr_srt;
  prg_nm=((prg_nm=std::strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++; // [sng] Name of program
  if(vrs_cvs == "Unknown") std::cerr << prg_nm << " version " << vrs_cpp << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  if(vrs_cvs != "Unknown") std::cerr << prg_nm << " version " << vrs_cvs << " last modified " << date_cvs << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  const std::string cmd_ln(cmd_ln_sng(argc,argv)); // [sng] Parsed command line
  std::cerr << "Command Line: " << cmd_ln << std::endl;

  // Set defaults for command line options
  bool flg_lmt(true); // [flg] Option lmt, no_lmt
  bool flg_slr_spc(true); // [flg] Option slr, no_slr
  bool flg_ncr_wvl(true); // [flg] Option ncr_wvl, dcr_wvl
  float flt_foo_input(0.0); // Option f  
  prc_cmp slr_cst(1367.0); // [W m-2] Option slr_cst
  prc_cmp wvl_lmt_mcr(5.0); // [um] Option wvl_lmt_mcr
  prc_cmp wvl_mxm_mcr(1.1); // [um] Option wvl_mxm_mcr
  prc_cmp wvl_mnm_mcr(1.0); // [um] Option wvl_mnm_mcr
  long bnd_nbr(10); // [nbr] Option bnd_nbr
  long wvl_nbr(10); // [nbr] Option wvl_nbr
  std::string fl_err("cerr"); // Option e
  std::string fl_in("/data/zender/aca/spc_Kur95_20wvn.txt"); // Option i
  std::string fl_out("/data/zender/aca/out.nc"); // Option o

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}, 
       If flag is non-zero, getopt_long() returns zero and flag is set to val.
       If flag is zero, getopt_long() returns contents of val.
       */
    // Long options with no argument
    {"dcr_wvl",0,0,0},
    {"ncr_wvl",0,0,0},
    {"slr_spc",0,0,0},
    {"no_slr_spc",0,0,0},
    // Long options with argument
    {"bnd_nbr",1,0,0},
    {"flx_slr_src",1,0,0},
    {"wvl_mxm_mcr",1,0,0},
    {"wvl_mnm_mcr",1,0,0},
    {"wvl_lmt_mcr",1,0,0},
    {"wvl_nbr",1,0,0},
    // Long options with optional argument
    // Long options with short counterparts
    {"slr_cst",1,0,'s'},
    {"dbg",2,0,'D'},
    {"error",1,0,'e'},
    {"help",0,0,'h'},
    {"input",1,0,'i'},
    {"output",1,0,'o'},
    {"version",0,0,'v'},
    // Last option named "0" to signal getopt_long() to stop processing  
    {0,0,0,0}
  }; // end opt_lng
  
  // Short options: no colon = no arg, one colon = required arg, two colons = optional arg
  const char *opt_short_lst="D::e:f:hi:o:s:t:v"; // List of single-letter (C-style) option abbreviations
  extern char *optarg; // char * representation of current optarg, if any (this memory is owned by system)
  extern int optind; // extern enumerating cardinal of current option
  int opt; // Value is zero if current argument is long type, else value contains single letter version of command line argument
  int opt_idx=0; // Index of current long option into opt_lng array
  std::string opt_crr; // String representation of current long-option name
  std::string opt_sng; // String representation of current optarg, if any
 
  // Parse command line arguments 
  while(1){
    // getopt_long_only() allows a single dash '-' to prefix long options as well
    opt=getopt_long_only(argc,argv,opt_short_lst,opt_lng,&opt_idx);
    // NB: access to opt_crr is only valid when long_opt was detected
    opt_crr=opt_lng[opt_idx].name;  
    if(optarg) opt_sng=optarg; // Change C string into C++ string
    if(opt == EOF) break; // Parse positional arguments once getopt_long_only() returns EOF
    // Process long options without short option counterparts
    if(opt == 0){
      if(dbg_lvl >= dbg_io) std::cerr << "Long option name: " << opt_crr << (optarg ? ",  Argument: "+opt_sng : ", No Argument") << std::endl;
      if(opt_crr == "bnd_nbr") bnd_nbr=strtol(optarg,(char **)NULL,10);
      if(opt_crr == "ncr_wvl") flg_ncr_wvl=true;
      if(opt_crr == "dcr_wvl") flg_ncr_wvl=false;
      if(opt_crr == "slr_spc") flg_slr_spc=true;
      if(opt_crr == "no_slr_spc") flg_slr_spc=false;
      if(opt_crr == "lmt") flg_lmt=true;
      if(opt_crr == "no_lmt") flg_lmt=false;
      if(opt_crr == "wvl_nbr") wvl_nbr=strtol(optarg,(char **)NULL,10);
      if(opt_crr == "wvl_mnm_mcr") wvl_mnm_mcr=static_cast<prc_cmp>(strtod(optarg,(char **)NULL));
      if(opt_crr == "wvl_lmt_mcr") wvl_lmt_mcr=static_cast<prc_cmp>(strtod(optarg,(char **)NULL));
      if(opt_crr == "wvl_mxm_mcr") wvl_mxm_mcr=static_cast<prc_cmp>(strtod(optarg,(char **)NULL));
      if(opt_crr == "flx_slr_src"){
	// Get author of solar flux source. Default is ThD71
	if((flx_slr_sng[ThD71].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("ThD71") != std::string::npos)) 
	  flx_slr_src=ThD71;
	else if((flx_slr_sng[LaN68].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("LaN68") != std::string::npos)) 
	  flx_slr_src=LaN68;
	else if((flx_slr_sng[NeL84].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("NeL84") != std::string::npos)) 
	  flx_slr_src=NeL84;
	else if((flx_slr_sng[Kur95].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("Kur95") != std::string::npos)) 
	  flx_slr_src=Kur95;
	else err_prn(prg_nm,sbr_nm,"Unknown flx_slr_src");
      } // end if "flx_slr_src"
    } // opt != 0
    
    switch(opt){
    case 0: // Long options have already been processed, return
      break;
    case 'D': // The debugging level.  Default is 0. 
      if(optarg) dbg_lvl=(unsigned short int)strtol(optarg,(char **)NULL,10); else dbg_lvl=1;
      break;
    case 'f': // Set the generic tuning parameter.  Default is 0. 
      flt_foo_input=static_cast<float>(strtod(optarg,(char **)NULL));
      break;
    case 'i': // Get the input file name. Default is cout 
      fl_in=optarg;
      break;
    case 's': // Set the solar constant.  Default is 1367 W m-2
      slr_cst=static_cast<prc_cmp>(strtod(optarg,(char **)NULL));
      break;
    case 'h': // Print proper usage 
      //      (void)usg_prn(opt_sng);
      return EXIT_SUCCESS;
    case 'o': // Get the output file name. Default is cout 
      fl_out=optarg;
      break;
    case 'v': // Print the CVS program info 
      std::cerr << CVS_Header << std::endl;
      return EXIT_SUCCESS;
      break;
    default: // Print proper usage 
      //      (void)usg_prn(opt_sng);
      return EXIT_FAILURE;
    } // end switch 
  } // end while loop 
  
  // Process the positional arguments  
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      err_prn(prg_nm,sbr_nm,"Too many positional arguments");
    }else if(psn_arg_nbr == 1){
      fl_in=argv[optind++];
    }else if(psn_arg_nbr == 2){
      fl_in=argv[optind++];
      fl_out=argv[optind++];
    } // end else
  } // end if  

  // Convert input data to SI units
  //  prc_cmp wvl_lmt(wvl_lmt_mcr*1.0e-6); // um -> m
  //  prc_cmp wvl_min=wvl_mnm_mcr*1.0e-6; // um -> m
  //  prc_cmp wvl_max=wvl_mxm_mcr*1.0e-6; // um -> m
  //  prc_cmp wvl_ctr=0.5*(wvl_min+wvl_max);
  //  prc_cmp wvl_ncr=(wvl_max-wvl_min)/wvl_nbr;

  // Main body of code  
  dbg_lvl_tst();

  if(flx_slr_src == LaN68){
    wvl_nbr=122;
    // 20240808 Change to new Fortran function names
    flx_frc_fnc=FORTRAN_slffln; // [fnc] Fractional solar flux function
  }else if(flx_slr_src == ThD71){
    wvl_nbr=68;
    // 20240808 Change to new Fortran function names
    flx_frc_fnc=FORTRAN_slfftd; // [fnc] Fractional solar flux function
  }else if(flx_slr_src == NeL84){
    wvl_nbr=921;
    flx_frc_fnc=FORTRAN_slffln; // Bogus, not used
  }else if(flx_slr_src == Kur95){
    flx_frc_fnc=FORTRAN_slffln; // Bogus, not used
  } // end else
  if(dbg_lvl > 1) std::cerr << "flx_slr_src = " << flx_slr_sng[flx_slr_src] << std::endl;

  // Create standard arrays
  prc_cmp *wvl=new prc_cmp[wvl_nbr]; // [m] Nominal wavelength
  prc_cmp *wvl_dlt=new prc_cmp[wvl_nbr]; // [m] Bandwidth
  prc_cmp *wvl_min=new prc_cmp[wvl_nbr]; // [m] Minimum wavelength in band
  prc_cmp *wvl_max=new prc_cmp[wvl_nbr]; // [m] Maximum wavelength in band
  prc_cmp *wvl_ctr=new prc_cmp[wvl_nbr]; // [m] Wavelength at band center
  prc_cmp *flx_bnd=new prc_cmp[wvl_nbr]; // [W m-2] Solar flux in band

  if(flx_slr_src == LaN68){ // LaN68
    fl_in="RAM-based";
    prc_cmp wvl_max_tmp[122]={
      0.200,0.205,0.215,0.225,0.235,0.245,0.255,0.265,0.275,
      0.285,0.295,
      0.305,0.315,0.325,0.335,0.345,0.355,0.365,0.375,0.385,0.395,
      0.405,0.415,0.425,0.435,0.445,0.455,0.465,0.475,0.485,0.495,
      0.505,0.515,0.525,0.535,0.545,0.555,0.565,0.575,0.585,0.595,
      0.605,0.615,0.625,0.635,0.645,0.655,0.665,0.675,0.685,0.695,
      0.705,0.715,0.725,0.735,0.745,0.755,0.765,0.775,0.785,0.795,
      0.805,0.815,0.825,0.835,0.845,0.855,0.865,0.875,0.885,0.895,
      0.905,0.915,0.925,0.935,0.945,0.955,0.965,0.975,0.985,0.995,
      1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,
      2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,
      3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,
      4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,
      5.00
    };
    
    // Flux shortwave (bluer) than corresponding wavelength
    prc_cmp flx_frc_blr_tmp[122]={ // [frc] Fraction of solar flux at shorter wavelengths
      0.0000,
      0.0001,0.0003,0.0006,0.0010,0.0015,0.0020,0.0029,0.0042,0.0059,0.0088,
      0.0127,0.0172,0.0227,0.0291,0.0358,0.0427,0.0502,0.0580,0.0654,0.0732,
      0.0835,0.0959,0.1085,0.1209,0.1343,0.1489,0.1638,0.1786,0.1930,0.2073,
      0.2217,0.2356,0.2493,0.2633,0.2774,0.2911,0.3047,0.3183,0.3318,0.3451,
      0.3581,0.3709,0.3834,0.3956,0.4076,0.4191,0.4304,0.4417,0.4528,0.4636,
      0.4741,0.4844,0.4945,0.5043,0.5139,0.5232,0.5324,0.5414,0.5502,0.5588,
      0.5672,0.5755,0.5835,0.5913,0.5989,0.6062,0.6134,0.6204,0.6273,0.6340,
      0.6407,0.6472,0.6536,0.6598,0.6660,0.6719,0.6778,0.6835,0.6892,0.6947,
      0.7230,0.7671,0.8034,0.8336,0.8590,0.8808,0.8992,0.9144,0.9269,0.9372,
      0.9457,0.9529,0.9590,0.9642,0.9686,0.9724,0.9757,0.9786,0.9811,0.9834,
      0.9853,0.9870,0.9885,0.9899,0.9911,0.9922,0.9931,0.9940,0.9948,0.9955,
      0.9962,0.9968,0.9973,0.9978,0.9982,0.9987,0.9990,0.9994,0.9997,0.9999,
      1.0
    };

    // Convert input data to SI units
    for(idx=0;idx<wvl_nbr;idx++) wvl_max_tmp[idx]*=1.0e-6;

    if(dbg_lvl > 1){
      std::cerr << "idx\twvl_max\tflx_frc_blr" << std::endl;
      std::cerr << "#\tum\t" << std::endl;
      for(idx=0;idx<wvl_nbr;idx++) std::cerr << idx << "\t" << wvl_max_tmp[idx]*1.0e6 << "\t" << flx_frc_blr_tmp[idx] << std::endl;
    } // end if dbg

    // Create standard fields
    for(idx=0;idx<wvl_nbr;idx++){
      wvl_max[idx]=wvl_max_tmp[idx];
    } // end loop over wvl 

    wvl_min[0]=wvl_max[0]-(wvl_max[1]-wvl_max[0]);
    for(idx=1;idx<wvl_nbr;idx++){
      wvl_min[idx]=wvl_max[idx-1];
    } // end loop over wvl 

    for(idx=0;idx<wvl_nbr;idx++){
      wvl[idx]=wvl_ctr[idx]=0.5*(wvl_max[idx]+wvl_min[idx]);
      wvl_dlt[idx]=wvl_max[idx]-wvl_min[idx];
    } // end loop over wvl 

    flx_bnd[0]=flx_frc_blr_tmp[0]; // [W m-2] Solar flux in band
    for(idx=1;idx<wvl_nbr;idx++){
      flx_bnd[idx]=flx_frc_blr_tmp[idx]-flx_frc_blr_tmp[idx-1]; // [W m-2] Solar flux in band
    } // end loop over wvl 

  } // end LaN68

  if(flx_slr_src == ThD71){ // ThD71
    fl_in="RAM-based";

    // CSZ added (wvl_mcr,flx_frc_blr) = (10.0,1.0)
    // As input by BPB, last element was (5.0,0.9951)

    const long wvl_nbr_ThD71(68); // [nbr] Number of wavelength bins
    wvl_nbr=wvl_nbr_ThD71; // [nbr] Number of wavelength bins
    prc_cmp wvl_max_tmp[wvl_nbr_ThD71]={
      0.15,
      0.20,0.22,0.23,0.24,0.25,0.26,0.27,
      0.28,
      0.29,0.30,0.31,0.32,0.33,0.34,0.35,
      0.36,0.37,0.38,0.39,0.40,0.41,0.42,
      0.43,0.44,0.45,0.46,0.47,0.48,0.49,
      0.50,0.51,0.52,0.53,0.54,0.55,0.56,
      0.57,
      0.58,0.59,0.60,0.62,0.64,0.66,0.68,
      0.70,0.72,0.75,0.80,0.90,1.00,1.20,
      1.40,1.60,1.80,2.00,2.20,2.40,2.60,
      2.80,3.00,3.20,3.40,3.60,3.80,4.00,
      5.00,
      10.00
    };
    
    prc_cmp flx_frc_blr_tmp[wvl_nbr_ThD71]={ // [frc] Fraction of solar flux at shorter wavelengths
      0.0000,
      0.0001,0.0005,0.0010,0.0014,0.0019,0.0027,0.0041,
      0.0056,
      0.0081,0.0121,0.0165,0.0222,0.0293,0.0372,0.0452,
      0.0532,0.0615,0.0700,0.0782,0.0873,0.0992,0.1122,
      0.1247,0.1373,0.1514,0.1665,0.1817,0.1968,0.2115,
      0.2260,0.2401,0.2538,0.2674,0.2808,0.2938,0.3065,
      0.3191,
      0.3318,0.3444,0.3568,0.3810,0.4042,0.4266,0.4481,
      0.4688,0.4886,0.5169,0.5602,0.6336,0.6946,0.7839,
      0.8434,0.8861,0.9159,0.9349,0.9483,0.9589,0.9667,
      0.9731,0.9783,0.9822,0.9850,0.9872,0.9891,0.9906,
      0.9951,
      1.0
    };

    // Convert input data to SI units
    for(idx=0;idx<wvl_nbr;idx++) wvl_max_tmp[idx]*=1.0e-6;

    if(dbg_lvl > 1){
      std::cerr << "idx\twvl_max\tflx_frc_blr" << std::endl;
      std::cerr << "#\tum\t" << std::endl;
      for(idx=0;idx<wvl_nbr;idx++) std::cerr << idx << "\t" << wvl_max_tmp[idx]*1.0e6 << "\t" << flx_frc_blr_tmp[idx] << std::endl;
    } // end if dbg

    // Create standard fields
    for(idx=0;idx<wvl_nbr;idx++){
      wvl_max[idx]=wvl_max_tmp[idx]; // [m] Maximum wavelength in band
    } // end loop over wvl 

    wvl_min[0]=wvl_max[0]-(wvl_max[1]-wvl_max[0]);
    for(idx=1;idx<wvl_nbr;idx++){
      wvl_min[idx]=wvl_max[idx-1]; // [m] Minimum wavelength in band
    } // end loop over wvl 

    for(idx=0;idx<wvl_nbr;idx++){
      wvl[idx]=wvl_ctr[idx]=0.5*(wvl_max[idx]+wvl_min[idx]); // [m] Wavelength at band center
      wvl_dlt[idx]=wvl_max[idx]-wvl_min[idx]; // [m] Bandwidth
    } // end loop over wvl 

    flx_bnd[0]=flx_frc_blr_tmp[0]; // [W m-2] Solar flux in band
    for(idx=1;idx<wvl_nbr;idx++){
      flx_bnd[idx]=flx_frc_blr_tmp[idx]-flx_frc_blr_tmp[idx-1]; // [W m-2] Solar flux in band
    } // end loop over wvl 

  } // end ThD71

  if(flx_slr_src == NeL84){ // NeL84
    /* NB: 19970507 

       I carefully examined the NeL84 input file from WDC because
       there are numerical noise spikes apparent when plotting the
       derived solar spectrum.  Most of the spikes in the resulting
       solar spectrum seemed to be caused by numerical artifacts in
       the data, in which only 4-significant digits of flux are
       retained and differenced and spread onto 4--5 significant digits
       of wavelength.  Combining this with the high wavelength
       resolution results in significant jumps in the derived spectral
       flux.  So in the ultraviolet the precision is very high,
       roughly 0.001 W m-2 nm-1, while in the NIR the precision is
       degraded to 1.000 W m-2 (10 nm)-1.

       Performing the operations in double precision would not solve
       this problem.  More siginificant figures are required in the
       flux values relative to the wavelength values.

       The problem is artifical in one way: Most of the solar spectrum
       is not very numerically noisy, and it oscillates around the
       correct values.  Thus Radiative transfer routines run at lower
       wavelength resolution than the input data will not see this
       noise because the spectral flux they see is interpolated from
       the accumulated flux between the wavelength bounds supplied by
       the RT code.

       The spike at 330 nm can be seen comparing the output of the following:
       ncks -C -H -d wvl,163,165 -v flx_bnd,flx_frc_blr,wvl_min,wvl_dlt,wvl_max,flx_spc_wvl ${DATA}/aca/spc_NeL84.nc | m
       to the corresponding input data:
       162	330.0 .3586e+02
       163	330.4 .3597e+02
       164	330.5 .3698e+02
       165      331.5 .3795e+02 
       This spike is exactly where spectrum of Neckel and Labs is joined to spectrum of Arveson
       This suggests this particular spike, at least, is due to differences in measurements

       Spike at 630 nm appears due to changing wavelength resolution from 1 nm to 2 nm
       Spike at 1 um appears due to changing wavelength resolution from 2 nm to 5 nm
       Remaining spikes are numerical noise due to loss of precision in flux

       */

    prc_cmp *flx_blr=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux at shorter wavelengths

    fp_in=fopen(fl_in.c_str(),"r"); 
    if(fp_in == NULL) err_prn(prg_nm,sbr_nm,"Unable to open "+fl_in);
    for(idx=0;idx<22;idx++){
      rcd=fscanf(fp_in," %[^\n] ",bff);
      if(dbg_lvl > 1) (void)fprintf(stderr,"%s\n",bff);
    } // end loop over idx
    
    // Read until EOF or until allocated memory runs out
    idx=0L;
    while(rcd != EOF && idx < wvl_nbr){
      rcd=fscanf(fp_in,"%f %f ",&flt_foo,&flt_foo2);
      if(rcd > 0){
	// Convert input data to SI units
	wvl_max[idx]=flt_foo*1.0e-9; // nm --> m
	flx_blr[idx]=flt_foo2; // [frc] Fraction of solar flux at shorter wavelengths
	idx++;
      } // end if
    } // end while
    std::cerr << "Ingested data from " << fl_in << std::endl;
    (void)fclose(fp_in); 
    
    // Set size of wavelength and wavenumber arrays to number in input file
    wvl_nbr=idx;

    if(dbg_lvl > 1){
      std::cerr << "idx\twvl_max\tflx_blr" << std::endl;
      std::cerr << "#\tnm\tW m-2" << std::endl;
      for(idx=0;idx<wvl_nbr;idx++) std::cerr << idx << "\t" << wvl_max[idx]*1.0e9 << "\t" << flx_blr[idx] << std::endl;
    } // end if dbg

    // Create standard fields
    wvl_min[0]=wvl_max[0]-(wvl_max[1]-wvl_max[0]); // [m] Minimum wavelength in band
    for(idx=1;idx<wvl_nbr;idx++){
      wvl_min[idx]=wvl_max[idx-1]; // [m] Minimum wavelength in band
    } // end loop over wvl 

    for(idx=0;idx<wvl_nbr;idx++){
      wvl[idx]=wvl_ctr[idx]=0.5*(wvl_max[idx]+wvl_min[idx]); // [m] Wavelength at band center
      wvl_dlt[idx]=(double)wvl_max[idx]-(double)wvl_min[idx]; // [m] Bandwidth
    } // end loop over wvl 

    flx_bnd[0]=flx_blr[0]; // [W m-2] Solar flux in band
    for(idx=1;idx<wvl_nbr;idx++){
      flx_bnd[idx]=(double)flx_blr[idx]-(double)flx_blr[idx-1]; // [W m-2] Solar flux in band
    } // end loop over wvl 

    delete []flx_blr; // [frc] Fraction of solar flux at shorter wavelengths
  } // end NeL84

  if(flx_slr_src == Kur95){ // Kur95
    prc_cmp *wvn_ctr=new prc_cmp[wvl_nbr]; // [cm-1] Wavenumber at band center
    prc_cmp *flx_spc_wvn=new prc_cmp[wvl_nbr];

    fp_in=fopen(fl_in.c_str(),"r"); 

    for(idx=0;idx<7;idx++){
      rcd=fscanf(fp_in," %[^\n] ",bff);
      if(dbg_lvl > 1) (void)fprintf(stderr,"%s\n",bff);
    } // end loop over idx
    
    // Read until EOF or until allocated memory runs out
    idx=0L;
    while(rcd != EOF && idx < wvl_nbr){
      rcd=fscanf(fp_in,"%f %f ",&flt_foo,&flt_foo2);
      //    if(rcd == EOF) wrn_prn("Reached EOF in input file");
      if(rcd > 0){
	wvn_ctr[idx]=flt_foo; // [cm-1] Wavenumber at band center
	flx_spc_wvn[idx]=flt_foo2;
	idx++;
      } // end if
    } // end while
    std::cerr << "Ingested data from " << fl_in << std::endl;
    (void)fclose(fp_in); 

    // Set size of wavelength and wavenumber arrays to number in input file
    wvl_nbr=idx; // [nbr] Number of wavelength bands
  
    if(dbg_lvl > 1){
      std::cerr << "idx\twvn_ctr\tflx_spc_wvn" << std::endl;
      std::cerr << "#\tcm-1\tW m-2/cm-1" << std::endl;
      for(idx=0;idx<wvl_nbr;idx++) std::cerr << idx << "\t" << wvn_ctr[idx] << "\t" << flx_spc_wvn[idx] << std::endl;
    } // end if dbg

    prc_cmp *wvn_dlt=new prc_cmp[wvl_nbr]; // [cm-1] Bandwidth
    prc_cmp *wvn_min=new prc_cmp[wvl_nbr]; // [cm-1] Minimum wavenumber in band
    prc_cmp *wvn_max=new prc_cmp[wvl_nbr]; // [cm-1] Maximum wavenumber in band

    prc_cmp wvn_ncr=wvn_ctr[1]-wvn_ctr[0];
    // Create standard fields
    for(idx=0;idx<wvl_nbr;idx++){
      wvn_dlt[idx]=wvn_ncr; // [cm-1] Bandwidth
      wvn_min[idx]=wvn_ctr[idx]-0.5*wvn_ncr; // [cm-1] Minimum wavenumber in band
      wvn_max[idx]=wvn_ctr[idx]+0.5*wvn_ncr; // [cm-1] Maximum wavenumber in band
      wvl[idx]=1.0/(100.*wvn_ctr[idx]); // [m] Nominal wavelength
      wvl_ctr[idx]=wvl[idx]; // [m] Wavelength at band center
      wvl_min[idx]=1.0/(100.*wvn_max[idx]); // [m] Minimum wavelength in band
      wvl_max[idx]=1.0/(100.*wvn_min[idx]); // [m] Maximum wavelength in band
      wvl_dlt[idx]=wvl_max[idx]-wvl_min[idx]; // [m] Bandwidth
      flx_bnd[idx]=wvn_dlt[idx]*flx_spc_wvn[idx]; // [W m-2] Solar flux in band
    } // end loop over wvl 

    delete []wvn_ctr; // [m] Wavelength at band center
    delete []flx_spc_wvn;
  } // end if Kur95
  
  // Make sure standard fields are in requested order
  if(flg_ncr_wvl && flx_slr_src == Kur95){
    std::cerr << "Reversing input data to be in order of increasing wavelength" << std::endl;
    (void)rvr_vec(wvl,wvl_nbr);
    (void)rvr_vec(wvl_dlt,wvl_nbr);
    (void)rvr_vec(wvl_min,wvl_nbr);
    (void)rvr_vec(wvl_max,wvl_nbr);
    (void)rvr_vec(wvl_ctr,wvl_nbr);
    (void)rvr_vec(flx_bnd,wvl_nbr);
  } // end if flg_ncr_wvl

  if(flg_ncr_wvl){
    idx_max_wvl=wvl_nbr-1L;
    idx_min_wvl=0L;
    ncr_ncr_wvl=1L;
  }else{
    idx_max_wvl=0L;
    idx_min_wvl=wvl_nbr-1L;
    ncr_ncr_wvl=-1L;
  } // endif flg_ncr_wvl

  // Derive fields from standard fields

  prc_cmp *wvl_mcr=new prc_cmp[wvl_nbr];
  for(idx=0;idx<wvl_nbr;idx++){
    wvl_mcr[idx]=wvl[idx]*1.0e6;
  } // end loop over wvl 

  // Create wavenumber grid
  long wvn_nbr=wvl_nbr;
  prc_cmp *wvn=new prc_cmp[wvn_nbr]; // [cm-1] Nominal wavenumber
  prc_cmp *wvn_ctr=new prc_cmp[wvn_nbr]; // [cm-1] Wavenumber at band center
  prc_cmp *wvn_dlt=new prc_cmp[wvn_nbr]; // [cm-1] Bandwidth
  prc_cmp *wvn_min=new prc_cmp[wvn_nbr]; // [cm-1] Minimum wavenumber in band
  prc_cmp *wvn_max=new prc_cmp[wvn_nbr]; // [cm-1] Maximum wavenumber in band
  for(idx=0;idx<wvn_nbr;idx++){
    wvn_min[idx]=1.0/(100.0*wvl_max[idx]); // [cm-1] Minimum wavenumber in band
    wvn_max[idx]=1.0/(100.0*wvl_min[idx]); // [cm-1] Maximum wavenumber in band
    wvn_dlt[idx]=wvn_max[idx]-wvn_min[idx]; // [cm-1] Bandwidth
    wvn[idx]=wvn_ctr[idx]=0.5*(wvn_max[idx]+wvn_min[idx]); // [cm-1] Wavenumber at band center
  } // end loop over wvn 

  // Derive useful flux information
  prc_cmp flx_ttl=0.0; // [W m-2]
  prc_cmp flx_ttl_foo=0.0; // [W m-2]
  prc_cmp *flx_frc=new prc_cmp[wvl_nbr]; // [W m-2 m-1]
  prc_cmp *flx_spc_wvn=new prc_cmp[wvn_nbr]; // [W m-2 (cm-1)-1]
  prc_cmp *flx_spc_wvl=new prc_cmp[wvl_nbr]; // [W m-2 m-1]
  for(idx=0;idx<wvl_nbr;idx++){
    flx_spc_wvl[idx]=flx_bnd[idx]/wvl_dlt[idx]; // [W m-2 m-1]
    flx_spc_wvn[idx]=flx_bnd[idx]/wvn_dlt[idx]; // [W m-2 (cm-1)-1]
    flx_ttl+=flx_bnd[idx]; // [W m-2]
  } // end loop over wvl 
  for(idx=0;idx<wvl_nbr;idx++){
    /* Normalize input by total of input so that sum of flx_frc = 1.0
       Until this is done flx_bnd represents different quantities for different sources:
       LaN68 is already normalized so flx_bnd is flx_frc already
       ThD71 is already normalized so flx_bnd is flx_frc already
       Kur95 is absolute so flx_bnd is absolute flux */
    flx_frc[idx]=flx_bnd[idx]/flx_ttl; // [W m-2 m-1]
    flx_ttl_foo+=flx_ttl*flx_frc[idx]; // [W m-2]
  } // end loop over wvl 

  // Compute partial sums
  prc_cmp *flx_frc_rdr=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux at longer wavelengths
  prc_cmp *flx_frc_blr=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux at shorter wavelengths
  /* Ensure smallest wavelength bin has no flux so that rdr and blr equals [1,0] are exact!
     This is why I inserted first line of NeL84 file and last line of Kur95 files
     flx_frc_blr for bin "idx" is fraction of flux occuring at wavelengths less than wvl_max[idx] NOT wvl_ctr[idx] or wvl_min[idx]
     This convention must be supported by input text files
     This convention is similar to BPB's implementation of Lan68 and ThD71 */
  if(flx_frc[idx_min_wvl] != 0.0){
    wrn_prn(prg_nm,sbr_nm,"flx_frc[idx_min_wvl] != 0.0");
    std::cerr << "flx_frc[" << idx_min_wvl << "] = " << flx_frc[idx_min_wvl] << std::endl;
  } // endif
  if(flg_ncr_wvl){
    flx_frc_rdr[idx_min_wvl]=1.0; // [frc] Fraction of solar flux at longer wavelengths
    flx_frc_blr[idx_min_wvl]=0.0; // [frc] Fraction of solar flux at shorter wavelengths
    for(idx=1;idx<wvl_nbr;idx++){
      flx_frc_rdr[idx]=flx_frc_rdr[idx-1]-flx_frc[idx]; // [frc] Fraction of solar flux at longer wavelengths
      flx_frc_blr[idx]=flx_frc_blr[idx-1]+flx_frc[idx]; // [frc] Fraction of solar flux at shorter wavelengths
    } // end loop over wvl 
    flx_frc_rdr[idx_max_wvl]=0.0; // [frc] Fraction of solar flux at longer wavelengths
    flx_frc_blr[idx_max_wvl]=1.0; // [frc] Fraction of solar flux at shorter wavelengths
  }else{
    flx_frc_rdr[idx_max_wvl]=0.0; // [frc] Fraction of solar flux at longer wavelengths
    flx_frc_blr[idx_max_wvl]=1.0; // [frc] Fraction of solar flux at shorter wavelengths
    for(idx=1;idx<wvl_nbr;idx++){
      flx_frc_rdr[idx]=flx_frc_rdr[idx-1]+flx_frc[idx]; // [frc] Fraction of solar flux at longer wavelengths
      flx_frc_blr[idx]=flx_frc_blr[idx-1]-flx_frc[idx]; // [frc] Fraction of solar flux at shorter wavelengths
    } // end loop over wvl 
    flx_frc_rdr[idx_min_wvl]=1.0; // [frc] Fraction of solar flux at longer wavelengths
    flx_frc_blr[idx_min_wvl]=0.0; // [frc] Fraction of solar flux at shorter wavelengths
  } // end else

  // Normalize input absolute fluxes to standard solar constant
  prc_cmp nrm_fac=slr_cst/flx_ttl; // [frc]
  prc_cmp flx_ttl_foo2=0; // [W m-2]
  for(idx=0;idx<wvl_nbr;idx++){
    // flx_bnd currently contains fractional solar flux, normalize it to absolute flux
    flx_bnd[idx]*=nrm_fac; // [W m-2] Solar flux in band
    flx_spc_wvn[idx]*=nrm_fac; // [W m-2 (cm-1)-1]
    flx_spc_wvl[idx]*=nrm_fac; // [W m-2 m-1]
    flx_ttl_foo2+=flx_bnd[idx]; // [W m-2]
  } // end loop over wvl 

  if(true){
    std::cerr << "Initialization state:" << std::endl;
    std::cerr << "Input file: " << fl_in << std::endl;
    std::cerr << "Spectral solar flux source: "+flx_slr_sng[flx_slr_src] << std::endl;
    std::cerr << "Number of input spectral bins = " << wvn_nbr << std::endl;
    std::cerr << "Smallest wavenumber = " << wvn_min[idx_max_wvl] << " cm-1" << std::endl;
    std::cerr << "Largest wavenumber = " << wvn_max[idx_min_wvl] << " cm-1" << std::endl;
    std::cerr << "Shortest wavelength = " << wvl_min[idx_min_wvl]*1.0e6 << " um" << std::endl;
    std::cerr << "Longest wavelength = " << wvl_max[idx_max_wvl]*1.0e6 << " um" << std::endl;
    std::cerr << "Input solar constant = " << flx_ttl << " W m-2" << std::endl;
    std::cerr << "Sum of flx_frc = " << flx_ttl_foo << " W m-2" << std::endl;
    std::cerr << "User solar constant = " << slr_cst << " W m-2" << std::endl;
    std::cerr << "Output solar constant = " << flx_ttl_foo2 << " W m-2" << std::endl;
    std::cerr << "Number of output spectral bins = " << wvl_nbr << std::endl;
    std::cerr << "Number of sub-bins in each spectral bin = " << bnd_nbr << std::endl;
    std::cerr << "Output file: " << fl_out << std::endl;
    std::cerr << std::endl;
  } /* end if */

  prc_cmp *flx_frc_LaN68=new prc_cmp[wvl_nbr];
  prc_cmp *flx_frc_ThD71=new prc_cmp[wvl_nbr];
  prc_cmp *flx_frc_NeL84=new prc_cmp[wvl_nbr];
  prc_cmp *flx_frc_Kur95=new prc_cmp[wvl_nbr];
  prc_cmp *flx_slr=new prc_cmp[wvl_nbr];
  prc_cmp *flx_spc_slr=new prc_cmp[wvl_nbr];
  float wvl_min_mcr;
  float wvl_max_mcr;
  for(idx=0;idx<wvl_nbr;idx++){
    wvl_min_mcr=1.0e6*wvl_min[idx];
    wvl_max_mcr=1.0e6*wvl_max[idx];
    // NB: arguments in microns 
    flx_frc_LaN68[idx]=static_cast<prc_cmp>(FORTRAN_slffln(&wvl_min_mcr,&wvl_max_mcr));
    flx_frc_ThD71[idx]=static_cast<prc_cmp>(FORTRAN_slfftd(&wvl_min_mcr,&wvl_max_mcr));
    flx_frc_NeL84[idx]=0.0;
    flx_frc_Kur95[idx]=flx_frc[idx];
    flx_slr[idx]=slr_cst*flx_frc[idx];
    flx_spc_slr[idx]=flx_slr[idx]/wvl_dlt[idx];
  } // end loop over wvl 

  prc_cmp *wvl_wgt=new prc_cmp[wvl_nbr];
  for(idx=0;idx<wvl_nbr;idx++){
    // Compute non-weighted indices of refraction
    wvl_ctr_mcr=wvl[idx]*1.0e6;
    wvl_wgt[idx]=0.0;
  } // end loop over wvl 

  // Create arrays to contain sub-bands of each interval
  prc_cmp *bnd=new prc_cmp[bnd_nbr];
  prc_cmp *bnd_sz_mcr=new prc_cmp[bnd_nbr];
  prc_cmp *bnd_min_mcr=new prc_cmp[bnd_nbr];
  prc_cmp *bnd_max_mcr=new prc_cmp[bnd_nbr];
  prc_cmp *bnd_ctr_mcr=new prc_cmp[bnd_nbr];
  prc_cmp *flx_frc_bnd=new prc_cmp[bnd_nbr];
  prc_cmp *flx_slr_bnd=new prc_cmp[bnd_nbr];

  // For each spectral interval in the final output array
  if(dbg_lvl > 0) std::cerr << "Printing one `.' before each wvl loop:" << std::endl;
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){

    // Divide each output spectral interval into a number of subintervals called bands
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      bnd_sz_mcr[bnd_idx]=(wvl_max[wvl_idx]-wvl_min[wvl_idx])*1.0e6/bnd_nbr;
      bnd_min_mcr[bnd_idx]=wvl_min[wvl_idx]*1.0e6+bnd_idx*bnd_sz_mcr[bnd_idx];
      bnd_max_mcr[bnd_idx]=wvl_min[wvl_idx]*1.0e6+(bnd_idx+1)*bnd_sz_mcr[bnd_idx];
      bnd_ctr_mcr[bnd_idx]=0.5*(bnd_max_mcr[bnd_idx]+bnd_min_mcr[bnd_idx]);
      bnd[bnd_idx]=1.0e-6*bnd_ctr_mcr[bnd_idx];
    } // end loop over bnd
    
    // Routines which only take scalars
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      /* NB: Of course these numbers are only valid for LaN68, ThD71 or 
	 other sources which have spectra in fractional form */
      flx_frc_bnd[bnd_idx]=(*flx_frc_fnc)(static_cast<float *>(bnd_min_mcr+bnd_idx),static_cast<float *>(bnd_max_mcr+bnd_idx));
    } // end loop over bnd
      
    // Weight of current spectral interval
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      flx_slr_bnd[bnd_idx]=slr_cst*flx_frc_bnd[bnd_idx];
      wvl_wgt[wvl_idx]+=bnd_ctr_mcr[bnd_idx]*flx_frc_bnd[bnd_idx];
    } // end loop over bnd

    // Normalize by solar fraction to get solar spectral weighted average
    wvl_wgt[wvl_idx]/=flx_frc[wvl_idx];
    
  } // end loop over wvl 
  if(dbg_lvl > 0) std::cerr << std::endl;

  // Destroy arrays that contain sub-bands of each interval
  delete []bnd;
  delete []bnd_sz_mcr;
  delete []bnd_min_mcr;
  delete []bnd_max_mcr;
  delete []bnd_ctr_mcr;
  delete []flx_frc_bnd;
  delete []flx_slr_bnd;

  // Open output file
  const int nc_out(nco_create(fl_out,NC_CLOBBER)); 
  const nc_type nco_xtyp(nco_get_xtype(static_cast<prc_cmp>(1.0))); // [enm] External netCDF type
  if(dbg_lvl > dbg_off) std::cout << prg_nm << ": INFO External netCDF type of prc_cmp variables will be " << nco_typ_sng(nco_xtyp) << std::endl;

  // Create dimensions
  const int wvn_dmn=nco_def_dim(nc_out,static_cast<std::string>("wvn"),wvn_nbr);
  const int wvl_dmn=nco_def_dim(nc_out,static_cast<std::string>("wvl"),wvl_nbr);
  const int *dmn_wvl=&wvl_dmn;
  const int *dmn_wvn=&wvn_dmn;
  if(dbg_lvl > 3){
    rcd=nco_inq_ndims(nc_out,int_foo);
    std::cerr << "Currently there are " << int_foo << " dimensions defined" << std::endl;
  } // end if dbg
  const int *dmn_scl=&wvl_dmn; // dummy argument, not used

  // Global attributes
  rcd=nco_put_att(nc_out,NC_GLOBAL,"history",time_bfr_srt.substr(0,time_bfr_srt.size()-1)+": "+cmd_ln);
  rcd=nco_put_att(nc_out,NC_GLOBAL,"CVS_Id",CVS_Id);
  rcd=nco_put_att(nc_out,NC_GLOBAL,"solar_flux_source",flx_slr_sng[flx_slr_src]);

  var_mtd_sct var_mtd[]={
    {0,"flx_bnd",nco_xtyp,1,dmn_wvl,"long_name","Absolute solar flux in band","units","watt meter-2"},
    {0,"flx_frc",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux in band","units","fraction"},
    {0,"flx_frc_Kur95",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Kurucz (1995)","units","fraction"},
    {0,"flx_frc_LaN68",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Labs & Neckel (1968)","units","fraction"},
    {0,"flx_frc_NeL84",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Neckel & Labs (1984)","units","fraction"},
    {0,"flx_frc_ThD71",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Thekeakara & Drummond (1971)","units","fraction"},
    {0,"flx_frc_blr",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux at shorter wavelengths","units","fraction"},
    {0,"flx_frc_rdr",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux at longer wavelengths","units","fraction"},
    {0,"flx_slr",nco_xtyp,1,dmn_wvl,"long_name","Absolute solar flux in band","units","watt meter-2"},
    {0,"flx_spc_wvl",nco_xtyp,1,dmn_wvl,"long_name","Solar spectral flux per unit wavelength","units","watt meter-2 meter-1"},
    {0,"flx_spc_wvn",nco_xtyp,1,dmn_wvl,"long_name","Solar spectral flux per unit wavenumber","units","watt meter-2 meter-1"},
    {0,"slr_cst",nco_xtyp,0,dmn_scl,"long_name","Solar constant","units","watt meter-2"},
    {0,"wvl",nco_xtyp,1,dmn_wvl,"long_name","Wavelength at band center","units","meter"},
    {0,"wvl_dlt",nco_xtyp,1,dmn_wvl,"long_name","Bandwidth","units","meter"},
    {0,"wvl_max",nco_xtyp,1,dmn_wvl,"long_name","Maximum wavelength in band","units","meter"},
    {0,"wvl_min",nco_xtyp,1,dmn_wvl,"long_name","Minimum wavelength in band","units","meter"},
    {0,"wvl_wgt",nco_xtyp,1,dmn_wvl,"long_name","Solar flux-weighted wavelength","units","meter"},
    {0,"wvn",nco_xtyp,1,dmn_wvn,"long_name","Wavenumber at band center","units","centimeter-1"},
    {0,"wvn_dlt",nco_xtyp,1,dmn_wvn,"long_name","Bandwidth","units","centimeter-1"},
    {0,"wvn_max",nco_xtyp,1,dmn_wvn,"long_name","Maximum wavenumber in band","units","centimeter-1"},
    {0,"wvn_min",nco_xtyp,1,dmn_wvn,"long_name","Minimum wavenumber in band","units","centimeter-1"}
    //    0,"foo",nco_xtyp,1,dmn_wvl,"long_name","foo","units","foo",
  }; // end var_mtd
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct));
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr); // [fnc] Define variables in output netCDF file

  // Leave define mode, allow file to already be in define mode
  rcd=nco_enddef(nc_out,NC_ENOTINDEFINE); // [fnc] Leave define mode

  // Write data and delete dynamic arrays
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_bnd"),flx_bnd); delete []flx_bnd;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc"),flx_frc); delete []flx_frc;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_Kur95"),flx_frc_Kur95); delete []flx_frc_Kur95;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_LaN68"),flx_frc_LaN68); delete []flx_frc_LaN68;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_NeL84"),flx_frc_NeL84); delete []flx_frc_NeL84;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_ThD71"),flx_frc_ThD71); delete []flx_frc_ThD71;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_blr"),flx_frc_blr); delete []flx_frc_blr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_frc_rdr"),flx_frc_rdr); delete []flx_frc_rdr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr"),flx_slr); delete []flx_slr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_spc_wvl"),flx_spc_wvl); delete []flx_spc_wvl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("flx_spc_wvn"),flx_spc_wvn); delete []flx_spc_wvn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("slr_cst"),slr_cst);
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl"),wvl); delete []wvl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_dlt"),wvl_dlt); delete []wvl_dlt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_max"),wvl_max); delete []wvl_max;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_min"),wvl_min); delete []wvl_min;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_wgt"),wvl_wgt); delete []wvl_wgt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvn"),wvn); delete []wvn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_dlt"),wvn_dlt); delete []wvn_dlt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_max"),wvn_max); delete []wvn_max;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_min"),wvn_min); delete []wvn_min;

  // Close output file
  rcd=nco_close(nc_out); // [fnc] Close netCDF file
  std::cerr << "Wrote results to " << fl_out << std::endl;
  std::cerr << "ncks: ncks -C -H -F -m -u -d wvl,0.5e-6 -v wvl " << fl_out << std::endl;
    
  Exit_gracefully();
} // end main() 

void usg_prn(const char *opt_sng)
{
  std::cerr << "\nusage: " << prg_nm << " [-options] where options are one or more of:\n\n";
  std::cerr << opt_sng << "\n\n:";
  std::cerr << "-D dbg_lvl Set debugging level (default is 0)\n";
  std::cerr << "-f flt_foo Set generic float (default is 0)\n";
  std::cerr << "-i fl_in Set input file name. (default is in.nc)\n";
  std::cerr << "-o fl_out Set output file name. Default is out.nc\n";
  std::cerr << "-v print the CVS program version\n";
  std::cerr << std::endl;
} // end usg_prn() 
