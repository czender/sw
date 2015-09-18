// $Id$

// Purpose:  Template for using GSL, the GNU Scientific Library

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Compilation
   make -W gsl.cc OPTS=D gsl
   make OPTS=D gsl
   cd $HOME/c;make -W OPTS=D gsl.cc gsl;cd -
   cd $HOME/c;make OPTS=D gsl;cd - */

/* Usage:
   gsl --dbg=3 --dbl=1.0 */

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <complex> // Standard C++ complex class
#include <map> // STL multimap and map classes
//#include <typeinfo> // Standard C++ header for typeid, type_info

// Standard C headers
// #include <cstdio> // stderr, FILE, NULL, etc.
#include <ctime> // Machine time
#include <cstring> // strcmp()...
#include <cmath> // sin cos cos sin 3.14159
#include <cstdlib> // atof, atoi, strtol, strtod, malloc, getopt 
// #include <unistd.h> // All sorts of POSIX stuff 

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
// #define MAIN_PROGRAM_FILE MUST precede #include ccc.hh
#define MAIN_PROGRAM_FILE
#include <gsl.hh> // Program-specific definitions
#include <libcsz_c++.hh> // Personal C++ library
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Namespace
using namespace std; // std is default namespace 

// NB: Important to #define my stuff after system #include's
#define NUM(array) (sizeof(array)/sizeof(array[0]))

// Typedefs
typedef map<std::string,std::string,less<std::string> > sng2sng_map; // String-to-string map
typedef map<std::string,var_mtd_sct,less<std::string> > sng2var_mtd_map; // String-to-var_mtd_sct map

// Prototypes
extern "C" {
#include <gsl/gsl_sf_erf.h> // GNU Scientific Library special functions error functions
#include <gsl/gsl_sf_gamma.h> // GNU Scientific Library special functions gamma functions
#if (defined LINUX) || (defined LINUXALPHA) || (defined SGIMP64) || (defined SGI6)
#include <getopt.h> // GNU getopt() is standard on Linux
#else // not LINUX
#include "getopt.h" // GNU getopt()
#endif // not LINUX
} // end extern

int main(int argc,char **argv)
{
  // Prototypes  
  void usg_prn(const char *opt_sng);

  // Locals
  int idx;
  int int_foo;
  std::string CVS_Date="$Date$";
  std::string CVS_Header="$Id$";
  std::string CVS_Id="$Id$";
  std::string CVS_Revision="$Revision$";
  std::string date_sng=(CVS_Date.length() > 7 ? CVS_Date.substr(7,19) : "Unknown");
  std::string sbr_nm("main");
  std::string vrs_sng=(CVS_Revision.length() > 10 ? CVS_Revision.substr(10,4) : "Unknown");

  // Start the clock and save the command line  
  const std::time_t time_crr_time_t=std::time((std::time_t *)NULL);
  const std::string time_bfr_srt=std::ctime(&time_crr_time_t);
  std::cerr << "\tStart = " << time_bfr_srt;
  prg_nm=((prg_nm=std::strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++;
  std::cerr << "This is " << prg_nm << " version " << vrs_sng << " dated " << date_sng << std::endl;
  const std::string cmd_ln=cmd_ln_sng(argc,argv);
  std::cerr << "Command Line: " << cmd_ln << std::endl;

  // Set defaults for command line options
  bool flg_flg(false); // Option --flg_flg
  float flt_foo(0.0); // Option -f --flt_foo 
  double dbl_foo(0.0); // Option --dbl_foo 
  long bnd_nbr(1); // Option --bnd_nbr
  std::string fl_err="stderr"; // Option -e --error
  std::string fl_in="~/nc/nco/data/in.nc"; // Option -i --input
  std::string fl_out="foo.nc"; // Option -o --output
  Xtr_cls xtr_LHS("xtr_prt_wgt+xtr_fll_nil+xtr_vrb"); // Option xtr_LHS
  Xtr_cls xtr_RHS("xtr_prt_wgt+xtr_fll_nil+xtr_vrb"); // Option xtr_RHS

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}, 
       If flag is non-zero, getopt_long() returns zero and flag is set to val.
       If flag is zero, getopt_long() returns contents of val.
       */
    // Long options with no argument
    {"flg_flg",0,0,0},
    // Long options with argument
    {"bnd_nbr",1,0,0},
    {"function",1,0,0},
    {"xtr_LHS",1,0,0},
    {"xtr_RHS",1,0,0},
    {"dbl_foo",1,0,0},
    // Long options with optional argument
    // Long options with short counterparts
    {"dbg",2,0,'D'},
    {"flt_foo",1,0,'f'},
    {"error",1,0,'e'},
    {"help",0,0,'h'},
    {"input",1,0,'i'},
    {"output",1,0,'o'},
    {"version",0,0,'v'},
    // Last option named "0" to signal getopt_long() to stop processing  
    {0,0,0,0}
  }; // end opt_lng
  
  // Short options: no colon = no arg, one colon = required arg, two colons = optional arg
  const char *opt_short_lst="D::f:hi:o:v"; // List of single-letter (C-style) option abbreviations
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
      if(opt_crr == "flg_flg") flg_flg=true;
      if(opt_crr == "bnd_nbr") bnd_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10);
      if(opt_crr == "xtr_LHS") xtr_LHS.flg_set(opt_sng);
      if(opt_crr == "xtr_RHS") xtr_RHS.flg_set(opt_sng);
      if(opt_crr == "dbl_foo") dbl_foo=std::strtod(opt_sng.c_str(),(char **)NULL);
    } // opt != 0
    switch(opt){
    case 0: // Long options have already been processed, return
      break;
    case 'D': // Debugging level (default is 0) 
      if(optarg) dbg_lvl=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); else dbg_lvl=dbg_fl;
      break;
    case 'f': // Set generic tuning parameter (default is 0.0)
      flt_foo=static_cast<float>(std::strtod(opt_sng.c_str(),(char **)NULL));
      break;
    case 'i': // Get input file name (default is in.nc)
      fl_in=opt_sng;
      break;
    case 'h': // Print proper usage 
      (void)usg_prn(opt_short_lst);
      Exit_gracefully();
      break;
    case 'o': // Get output file name (default is foo.nc)
      fl_out=opt_sng;
      break;
    case 'v': // Print CVS program info 
      std::cerr << CVS_Header << std::endl;
      Exit_gracefully();
      break;
    default: // Print proper usage 
      (void)usg_prn(opt_short_lst);
      return EXIT_FAILURE;
    } // end switch 
  } // end while loop 
  
  // Process positional arguments  
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      err_prn(prg_nm,sbr_nm,"Too many positional arguments");
    }else if(psn_arg_nbr == 1){
      fl_out=argv[optind++];
    }else if(psn_arg_nbr == 2){
      fl_in=argv[optind++];
      fl_out=argv[optind++];
    } // end else
  } // end if  

  // Main body of code 
  int rcd(0); // Return code

  // Special function library takes mode agrument which can be GSL_PREC_DOUBLE, GSL_PREC_SINGLE, GSL_PREC_APPROX
  std::string nvr_GSL_PREC((std::getenv("GSL_PREC")) ? std::getenv("GSL_PREC") : ""); // [sng] Environment variable GSL_PREC
  gsl_sf_result answer; // GSL result structure
  bool apx_eql; // [flg] Arguments are indistinguishable
  double dbl_one(1.0); // [frc] One

  rcd+=gsl_sf_erf_e(dbl_one,&answer); // [fnc] Erro function
  if(std::fabs((0.8427-answer.val)/0.8427) > 0.001){
    std::cout << "erf(" << dbl_one << ") = " << answer.val << std::endl;
    std::cout << "gsl_sf_erf.val = " << answer.val << std::endl;
    std::cout << "gsl_sf_erf.err = " << answer.err << std::endl;
    std::cout << "rcd = " << rcd << std::endl;
    err_prn(prg_nm,sbr_nm,"ERROR Error function erf() error");
  } // endif err

  std::string tst_sng;
  if(dbg_lvl == dbg_crr || tst_sng == "erf"){
    std::cout << "Vetting ... error function" << std::endl;
    std::cout << "Seed values with --int foo, e.g., gsl --tst=erf --dbl_foo=1.0" << std::endl;
    std::string nfo_msg; // [sng] Descriptive message of context

    rcd+=gsl_sf_erf_e(dbl_one,&answer); // [fnc] Erro function
    nfo_msg="Vetting error function"; // [sng] Descriptive message of context
    apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
      (sbr_nm, // I [sng] Subroutine name of calling routine
       true, // I [flg] Verbose output
       0.8427, // I [frc] Target argument
       answer.val, // I [frc] Approximation to target argument
       1.0e-4, // I [frc] Relative precision
       "Vetting error function"); // I [sng] Descriptive message of context

    rcd+=gsl_sf_erf_e(dbl_foo,&answer);
    std::cout << "erf(" << dbl_foo << ") = " << answer.val << std::endl;
  } // end if dbg 
  
  rcd+=gsl_sf_gamma_e(dbl_one,&answer); // [fnc] Gamma function
  rcd+=gsl_sf_gamma_inc_Q_e(dbl_one,dbl_one,&answer); // [fnc] Normalized incomplete gamma function

  // Dummy data for netCDF output
  long wvl_nbr=1;
  prc_cmp *wvl=new prc_cmp[wvl_nbr];
  wvl[0]=0.5e-6;

  // Open output file
  const int nc_out(nco_create(fl_out,NC_CLOBBER)); 
 
  // Create dimensions
  const int wvl_dim=nco_def_dim(nc_out,static_cast<std::string>("wvl"),wvl_nbr);
  const int *dmn_wvl=&wvl_dim;
  if(dbg_lvl == dbg_old){
    rcd=nco_inq_ndims(nc_out,int_foo);
    std::cerr << "Currently there are " << int_foo << " dimensions defined" << std::endl;
  } // end if dbg
  const int *dmn_scl=&wvl_dim; // Dummy argument, not used
  // Global attributes
  rcd=nco_put_att(nc_out,NC_GLOBAL,"history",time_bfr_srt.substr(0,time_bfr_srt.size()-1)+": "+cmd_ln);
  rcd=nco_put_att(nc_out,NC_GLOBAL,"CVS_Id",CVS_Id);

  var_mtd_sct var_mtd[]={
    {0,"wvl",NC_FLOAT,1,dmn_wvl,"long_name","Wavelength at band center","units","meter"},
    {0,"wvl_nbr",NC_FLOAT,0,dmn_scl,"long_name","Number of wavelengths","units","number"}
  };
  const int var_mtd_nbr=sizeof(var_mtd)/sizeof(var_mtd_sct);
  
  sng2var_mtd_map var_map;
  for(idx=0;idx<var_mtd_nbr;idx++){
    /* fxm: Define variables before inserting into map, because map values 
       seem to be unwritable (read-only) once they are in map. */
    rcd=nco_def_var(nc_out,var_mtd[idx].nm,var_mtd[idx].type,var_mtd[idx].dmn_nbr,var_mtd[idx].dmn_id,var_mtd[idx].var_id);
    var_map.insert(sng2var_mtd_map::value_type(var_mtd[idx].nm,var_mtd[idx]));
  } // end loop over itr

  sng2var_mtd_map::const_iterator var_mtd_itr;
  for(var_mtd_itr=var_map.begin();var_mtd_itr!=var_map.end();++var_mtd_itr){
    // Write first attribute (long_name)
    rcd=nco_put_att(nc_out,var_mtd_itr->second.var_id,var_mtd_itr->second.att_1_nm,var_mtd_itr->second.att_1_val);
    // Write second attribute (units)
    rcd=nco_put_att(nc_out,var_mtd_itr->second.var_id,var_mtd_itr->second.att_2_nm,var_mtd_itr->second.att_2_val);
    if(dbg_lvl == dbg_old) std::cout << "Defined " << var_mtd_itr->first << " with long_name " << var_mtd_itr->second.att_1_val << std::endl;
  } // end loop over var_mtd_itr

  if(dbg_lvl == dbg_old) std::cerr << var_map.find("wvl")->second.att_1_val << std::endl;
  if(dbg_lvl == dbg_old) std::cerr << var_map.find("wvl_min")->second.att_1_val << std::endl;

  // Leave define mode
  rcd=nco_enddef(nc_out,NC_ENOTINDEFINE); // [fnc] Leave define mode

  // After writing, delete arrays 
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl"),wvl); delete []wvl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_nbr"),wvl_nbr);

  // Close output file
  rcd=nco_close(nc_out); // [fnc] Close netCDF file
  std::cerr << "Wrote results to " << fl_out << std::endl;
  std::cerr << "ncks: ncks -C -F -d wvl,0.5e-6 -v wvl " << fl_out << std::endl;

  Exit_gracefully();
  return EXIT_SUCCESS;
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

