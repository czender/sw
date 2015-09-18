// $Id$

// Purpose: Simulate GDZ09 integrating sphere geometry

/* Compilation:
   cd ~/c++;make snw;cd -
   g++ -O2 -g -Wall -Werror -o ${MY_BIN_DIR}/snw ~/c++/snw.cc
   xlC_r -O -q64 -o ${MY_BIN_DIR}/snw ${MY_OBJ_DIR}/getopt_bsd.o ~/c++/snw.cc */

/* Distribution:
   scp ~/c++/snw.cc dust.ess.uci.edu:c++
   scp ~/c++/snw.cc esmf.ess.uci.edu:c++ */

/* Usage: 
   snw
   snw --llm=drc
   snw --apx=bff --llm=drc
   snw --apx=bff --llm=drc --wvl=1310
   snw --apx=bff --llm=drc --wvl=1310 --rfl_srt=0.01 --rfl_dlt=0.1 --drc_dff_fct=0.90
   snw --azi_nbr=1 --dpt_nbr=1 --rds_nbr=1 --dbg_lvl=1
   snw --azi_nbr=10 --dpt_usr=0.0 --rds_usr=0.0 --dbg_lvl=1

   Production:
   drc_dff_fct obtained as ratio alb_drc/alb_dff from swnb2
   rfl_usr     obtained as       alb_dff         from swnb2
   snw --wvl=0635 --drc_dff_fct=0.97 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.89 --dpt_typ=GDZ09_13mm | m
   snw --wvl=0635 --drc_dff_fct=0.97 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.89 --dpt_typ=GDZ09_25mm | m
   Ascending order of BC concentration:
   snw --wvl=0635 --drc_dff_fct=0.97 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.89 --dpt_typ=GDZ09_25mm | m # ppb_000.0
   snw --wvl=0635 --drc_dff_fct=0.96 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.88 --dpt_typ=GDZ09_25mm | m # ppb_250.0
   snw --wvl=0635 --drc_dff_fct=0.96 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.86 --dpt_typ=GDZ09_25mm | m # ppb_750.0
   snw --wvl=0635 --drc_dff_fct=0.95 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.84 --dpt_typ=GDZ09_25mm | m # ppm_001.5
   snw --wvl=0635 --drc_dff_fct=0.91 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.70 --dpt_typ=GDZ09_25mm | m # ppm_007.5
   snw --wvl=0635 --drc_dff_fct=0.86 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.60 --dpt_typ=GDZ09_25mm | m # ppm_015.9
   snw --wvl=0635 --drc_dff_fct=0.57 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.18 --dpt_typ=GDZ09_25mm | m # ppm_221.0

   snw --wvl=1310 --frc_msn_drc=0.95 --drc_dff_fct=0.90 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.67 --dpt_typ=GDZ09_13mm | m
   snw --wvl=1310 --frc_msn_drc=0.95 --drc_dff_fct=0.90 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.67 --dpt_typ=GDZ09_25mm | m
   snw --wvl=1550 --frc_msn_drc=0.95 --drc_dff_fct=0.57 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.18 --dpt_typ=GDZ09_13mm | m
   snw --wvl=1550 --frc_msn_drc=0.95 --drc_dff_fct=0.57 --azi_nbr=10 --rds_nbr=10 --rfl_usr=0.18 --dpt_typ=GDZ09_25mm | m */

/* Evaluation:
   // Sample edge at base, zero azimuth, has small planar angle (6.3 dgr):
   snw --llm=dff --azi_usr=0 --dpt_usr=2.5 --rds_usr=1.9 --dbg_lvl=4
   // Sample edge at base, zero azimuth, has large planar angle (56.8 dgr):
   snw --llm=dff --azi_usr=90 --dpt_usr=2.5 --rds_usr=1.9 --dbg_lvl=4
   // Perfect agreement beween analytic and numerical integration:
   snw --llm=drc --azi_usr=90 --dpt_usr=2.5 --rds_nbr=100 --dbg_lvl=4

   fxm: azi_nbr > 2 gives sames psi_bar as azi_nbr=2, why?, e.g., change azi_nbr on
   snw --llm=drc --azi_nbr=4 --dpt_usr=2.5 --rds_usr=0.4 --dbg_lvl=4
   same problem occurs with llm=drc and dff
   snw --llm=dff --azi_nbr=4 --dpt_usr=2.5 --rds_usr=0.4 --dbg_lvl=4

   fxm: omega_bar is negative? I thought everything was positive-definite
   snw --llm=drc --azi_nbr=1 --dpt_usr=2.5 --rds_usr=0.4 --dbg_lvl=4

   From Ghislain Picard 20080525 diffusefraction2.dat:
   rfl  ryt       flt         flt+bff
   0.01 0.00465129 0.00472593 0.00472701 
   0.03 0.0140894 0.0143131 0.0143164 
   0.05 0.0237127 0.0240849 0.0240906 
   0.07 0.0335268 0.0340469 0.034055 
   0.09 0.0435372 0.0442046 0.0442152 
   0.11 0.05375 0.0545639 0.054577 
   0.21 0.108076 0.109604 0.109632 
   0.31 0.16851 0.170706 0.170752 
   0.41 0.236145 0.23893 0.238997 
   0.51 0.312348 0.315597 0.315692 
   0.61 0.398859 0.402378 0.402507 
   0.71 0.497919 0.501417 0.501589 
   0.81 0.612475 0.615508 0.615735 */

/* History: 
   20080114 Received HRL06 starter program from Ghislain Picard
   20080407 Trimmed to enable/use only HRL06 method
   20080420--20080424 Implemented HRL06 flat sample and baffle variants
   20080515 Allow command-line inputs
   20080516 Restore aperture size as parameter in indefinite integral and
            properly incorporate beam width as limits of integration
   20080527 Zero known bugs on planar angles
   20080529 Integrate areal-mean solid angle within radius loop */

/* Original message:
   En fait, ce n'est pas si negligable que ca. Voici un graphe (dont je ne suis pas totalement sur). En gros pour une reflectance de 50%, la sample recoit 30% de diffus (en plus de 100% de directe). Ce graphe est valable qq soit la longueur d'onde. 
   Ghislain 

   In fact, it's not so negligable as that. Here is a graph (of which I am not totally sure). Basically for a 50% reflectance, the sample receives 30% of diffuse (in addition to 100% of direct). This graph is valid qq is the wavelength.
   Ghislain */

// Standard C++ headers
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr
#include <string> // Standard C++ string class

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159
#include <cstdio> // stderr, EOF, FILE, NULL, etc.
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv

// Prototypes
extern "C" {
#include <getopt.h> // GNU getopt() functionality
} // end extern

int main(int argc,char **argv)
{
  // Functions
  double fnc_ngl_ntg_ndf_azi_90 // [fnc] Indefinite integral for psi_1+psi_2
    (const double a, // a = z_s = Depth in snow sample plus lip thickness
     const double b, // b = R_a = Radius of aperature to sample holder
     const double x); // x = r  = Radius from center of sample

  int // O [enm] Return success code
    fnc_ngl_psi_one_psi_two // [fnc] Azimuthal angles psi_1 and psi_2
    (const double rds, // I [cm] Radius from center of sample
     const double azi, // I [rdn] Azimuthal angle
     const double dpt, // I [cm] Depth in snow sample plus lip thickness
     const double rds_smp_apr, // I [cm] Radius of aperature to sample holder
     double &ngl_apr_psi_1, // O [rdn] Supplementary planar angle one subtended by point
     double &ngl_apr_psi_2); // O [rdn] Supplementary planar angle two subtended by point
  // end fnc_ngl_psi_one_psi_two() prototype

  // Enumerated types
  enum apx_typ{apx_smp,apx_flt,apx_bff}; // [enm] Approximation type
  const std::string apx_sng[]={"smp = Simplest sphere in HRL06","flt = Sphere with flat sample","bff = Sphere with flat sample and baffle to specular reflection"};
  int apx_typ(apx_bff); // [enm] Approximation type

  enum dpt_typ{dpt_rgl,dpt_sgl,dpt_GDZ09_13mm,dpt_GDZ09_25mm}; // [enm] Depth-grid type
  const std::string dpt_sng[]={"rgl = Regular (uniformly spaced) grid","usr = User-specified single depth","GDZ09_13mm = Stretched 32-layer, 13mm grid from GDZ09","GDZ09_25mm = Stretched 32-layer, 25mm grid from GDZ09"};
  int dpt_typ(dpt_rgl); // [enm] Depth-grid type

  enum llm_typ{llm_drc,llm_dff}; // [enm] Illumination type
  const std::string llm_sng[]={"drc = direct = collimated beam","dff = diffuse = isotropic incidence"};
  int llm_typ(llm_dff); // [enm] Illumination type

  enum wvl_typ{wvl_0635,wvl_1310,wvl_1550}; // [enm] Wavelength type
  const std::string wvl_sng[]={"0635 = 635 nm = 0.635 um","1310 = 1310 nm = 1.310 um","1550 = 1550 nm = 1.550 um"};
  int wvl_typ(wvl_1310); // [enm] Wavelength type

  // Parameters
  // Physical and numerical constants
  const double CEWI_dbl(9.9692099683868690e+36); // Compiler Error Warning Initializer for double
  const double cst_M_PIl(3.1415926535897932384626433832795029L); // [frc] 3
  const double cmxinch(2.54); // [cm] Centimeters per inch
  
  // Pre-set depth-grid choices
  // GDZ09 stretched depth-grids
  const double dpt_grd_GDZ09_13mm[]={ // [cm] Depth (centered)
    0.0005,0.0015,0.0025,0.0035,0.0045,
    0.00625,0.00875,0.01125,0.01375,0.0175,
    0.025,0.035,0.045,0.055,0.065,
    0.075,0.09,0.11,0.13,0.16,
    0.205,0.28,0.38,0.48,0.58,
    0.68,0.78,0.88,0.98,1.08,
    1.18,1.28};
  const double dpt_grd_GDZ09_25mm[]={ // [cm] Depth (centered)
    0.0005,0.0015,0.0025,0.0035,0.0045,
    0.00625,0.00875,0.01125,0.01375,0.0175,
    0.025,0.04,0.075,0.125,0.175,
    0.225,0.275,0.325,0.375,0.425,
    0.475,0.55,0.65,0.75,0.85,
    0.95,1.125,1.375,1.625,1.875,
    2.125,2.375};

  // Integrating sphere dimensions, properties from JC Gallet 20080424, 20080915
  // Blueprint in ${DATA}/icr/snw_lgge_ntg_sph.pdf
  const double dpt_smp_hld_13mm(1.33) ; // [cm] Depth of 13mm-deep "small" sample holder (it is not exactly 13mm deep)
  const double dpt_smp_hld_25mm(2.50) ; // [cm] Depth of 25mm-deep "large" sample holder
  const double hgt_lip(0.0); // [cm] Thickness of lip from sphere to sample top (z_s)
  const double rds_cll_0635nm(1.0/2); // [cm] Collimated (direct) beam radius (R_b)
  const double rds_cll_1310nm(1.0/2); // [cm] Collimated (direct) beam radius (R_b)
  const double rds_cll_1550nm(0.4/2); // [cm] Collimated (direct) beam radius (R_b)
  const double rds_dtc(0.1835); // [cm] Detector (photodiode) radius
  const double rds_smp_apr(1.5*cmxinch/2); // 1.905 [cm] Aperture (exposed sample) radius (R_a)
  const double rds_smp_hld(6.33/2); // 3.165 [cm] Radius of sample holder (R_E)
  const double rds_sph_ntg(6.0*cmxinch/2); // 7.62 [cm] Integrating sphere radius
  const double rds_src(0.58*cmxinch/2); // 0.736 [cm] Source (laser diode) radius
  const double rfl_wll_0635nm(0.9875); // [frc] Wall reflectance at 0635 nm
  const double rfl_wll_1310nm(0.972); // [frc] Wall reflectance at 1310 nm
  const double rfl_wll_1550nm(0.965); // [frc] Wall reflectance at 1550 nm

  // Set defaults for command line options 
  // Option name is variable name, e.g., --lng_foo=3, unless otherwise indicated
  bool flg_azi_usr(false); // [flg] Azimuthal angle user-specified
  bool flg_dpt_usr(false); // [flg] Depth user-specified
  bool flg_flg(false); // [flg] Temporary flag variable
  bool flg_rds_usr(false); // [flg] Radius user-specified
  bool flg_rfl_usr(false); // [flg] Reflectance user-specified
  double azi_dgr_usr(CEWI_dbl); // [dgr] Azimuthal angle user-specified
  double dbl_foo(0.0); // [frc] Intrinsic double temporary variable
  double dpt_usr(CEWI_dbl); // [cm] Snow depth user-specified
  double drc_dff_fct(0.9); // [frc] Direct reflectance as fraction of diffuse
  double frc_msn_drc(1.0); // [frc] Fraction radiation emitted as direct, collimated light
  double rds_usr(CEWI_dbl); // [cm] Radius user-specified
  double rfl_end(1.0); // [frc] Sample reflectance, diffuse incidence, end
  double rfl_srt(0.0); // [frc] Sample reflectance, diffuse incidence, start
  double rfl_usr(CEWI_dbl); // [frc] Reflectance user-specified
  float flt_foo(0.0f); // [frc] Intrinsic float temporary variable
  int int_foo(0); // [nbr] Intrinsic int temporary variable
  long azi_nbr(1L); // [nbr] Number of azimuthal angles
  long dpt_nbr(10L); // [nbr] Number of layers
  long lng_foo(0L); // [nbr] Intrinsic long temporary variable
  long rds_nbr(1L); // [nbr] Number of radii
  long rfl_nbr(10L); // [nbr] Number of reflectances
  std::string sng_foo(""); // [sng] Intrinsic string temporary variable
  std::string tst_sng(""); // [sng] Name of test to perform
  unsigned short dbg_lvl(1); // [enm] Debugging level

  // Derived constants
  const double azi_ntv(cst_M_PIl/2.0); // [frc] Azimuth integration interval
  const double azi_srt(0.0); // [frc] Azimuth integration lower bound
  const double dpt_srt(0.0); // [cm] Depth integration lower bound
  const double hgt_smp(rds_sph_ntg-sqrt(rds_sph_ntg*rds_sph_ntg-rds_smp_apr*rds_smp_apr)); // [cm] Height of sample surface above tangent plane
  const double rds_srt(0.0); // [frc] Radius integration lower bound
  const double sfc_sph(4.0*cst_M_PIl*rds_sph_ntg*rds_sph_ntg); // [cm2] Surface area of sphere

  // Doubly-derived constants
  const double sfc_dtc_nrm(cst_M_PIl*rds_dtc*rds_dtc/sfc_sph); // [frc] Normalized detector (photodiode) area
  const double sfc_src_nrm(cst_M_PIl*rds_src*rds_src/sfc_sph); // [frc] Normalized source (laser diode) area
  // Sample-area accounts for "spherical cap" effect
  const double sfc_smp_nrm(cst_M_PIl*(rds_smp_apr*rds_smp_apr+hgt_smp*hgt_smp)/sfc_sph); // [frc] Normalized sample area

  // Trebly-derived variables
  const double sfc_wll_nrm(1.0-sfc_smp_nrm-sfc_dtc_nrm-sfc_src_nrm); // [frc] Normalized wall area excluding holes HRL06 p. 5250

  // Locals requiring initialization 
  int rcd(0); // [enm] Return success code

  // Locals
  double azi; // [rdn] Azimuthal angle
  double dpt; // [cm] Snow depth
  double dpt_dlt(CEWI_dbl); // [cm] Depth increment
  double dpt_snw_pls_lip; // [cm] Depth in snow sample plus lip thickness
  double dpt_smp_hld(CEWI_dbl) ; // [cm] Depth of 13mm-deep "small" sample holder
  double fct_mtx_nvr(CEWI_dbl); // [frc] Matrix diagonalization factor
  double fct_tmp; // [frc] Temporary factor in matrix elements
  double flx_frc_drc; // [frc] Insolation fraction in direct beam (for swnb2)
  double frc_bck_RS; // [frc] Fraction radiation backscattered into state RS
  double frc_bck_RW; // [frc] Fraction radiation backscattered into state RW
  double frc_dff_dwn; // [frc] Fraction diffuse downwelling radiation (dff/drc)
  double ngl_apr_pln; // [rdn] Planar angle subtended by point
  double ngl_apr_psi_1; // [rdn] Supplementary planar angle one subtended by point
  double ngl_apr_psi_2; // [rdn] Supplementary planar angle two subtended by point
  double ngl_apr_sld; // [sr] Solid angle subtended from point to aperature
  double ngl_pln_avg(CEWI_dbl); // [rdn] Mean planar angle subtended by point
  double ngl_pln_dbl_avg(CEWI_dbl); // [rdn] Mean planar angle subtended in sphere by layer
  double ngl_pln_dbl_avg_azi_90; // [rdn] Mean planar angle subtended in sphere by layer for azimuth=90 degrees
  double ngl_sld_avg(CEWI_dbl); // [sr] Mean solid angle subtended in sphere by layer
  double ngl_sld_avg_azi_90; // [sr] Mean solid angle subtended in sphere by layer for azimuth=90 degrees
  double ntg_dfn_val; // [frc] Definite integral value
  double prb_RS_RW; // [frc] P(5,6) = Probability of RS->RW transition: photon reflected from sample then reflected by wall
  double prb_lng_RS_AD; // [frc] QR(1,1) = Probability of RS->AD transition: photon reflected first from sample enters detector
  double prb_lng_RW_AD; // [frc] QR(2,1) = Probability of RW->AD transition: photon reflected first from wall   enters detector
  double prb_lng_RW_AS; // [frc] QR(2,2) = Probability of RW->AS transition: photon reflected from wall then absorbed by sample
  double prb_lng_RW_RS; // [frc] QT(2,1) = Probability of RW->RS transition: photon reflected from wall then reflected by sample
  double prb_lng_RW_RS_or_RW_AS; // [frc] QT(2,1)+QR(2,2) = Union of RW->RS and RW->AS: photon reflected from wall incident on sample
  double rds; // [cm] Radius
  double rds_cll; // [cm] Collimated (direct) beam radius (R_b)
  double rds_llm; // [cm] Illuminated sample radius (R=R_b or R_a)
  double rds_ntv; // [frc] Radius integration interval
  double rfl_smp(CEWI_dbl); // [frc] Sample reflectance, diffuse incidence
  double rfl_smp_drc(CEWI_dbl); // [frc] Sample reflectance, direct beam
  double rfl_wll; // [frc] Wall reflectance
  double sgn_msr_RS_AD; // [frc] Measured signal from collimated illumination of sample
  double sgn_msr_RW_AD; // [frc] Measured signal from diffuse    illumination of sample
  long azi_idx; // [idx] Counting index for azimuth
  long dpt_idx; // [idx] Counting index for depth
  long rds_idx; // [idx] Counting index for radius
  long rfl_idx; // [idx] Counting index for reflectance

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}
       If flag is non-zero, getopt_long() returns zero and flag is set to val
       If flag is zero, getopt_long() returns contents of val */
    // Long options with no argument, no short option counterpart
    {"flg_flg",no_argument,0,0}, // [flg] Flag flag
    // Long options with argument, no short option counterpart
    {"apx_typ",required_argument,0,0}, // [enm] Approximation type
    {"azi_nbr",required_argument,0,0}, // [nbr] Number of azimuthal angles
    {"azi_usr_flg",required_argument,0,0}, // [flg] Azimuthal angle user-specified
    {"dbl_foo",required_argument,0,0}, // [nbr] Intrinsic double temporary variable 
    {"dpt_nbr",required_argument,0,0}, // [nbr] Number of layers
    {"dpt_typ",required_argument,0,0}, // [enm] Depth-grid type
    {"dpt_usr_flg",required_argument,0,0}, // [flg] Depth user-specified
    {"drc_dff_fct",required_argument,0,0}, // [frc] Direct reflectance as fraction of diffuse
    {"frc_msn_drc",required_argument,0,0}, // [frc] Fraction radiation emitted as direct, collimated light
    {"int_foo",required_argument,0,0}, // [nbr] Intrinsic int temporary variable
    {"llm_typ",required_argument,0,0}, // [enm] Illumination type
    {"lng_foo",required_argument,0,0}, // [nbr] Intrinsic long temporary variable 
    {"rds_nbr",required_argument,0,0}, // [nbr] Number of radii
    {"rds_usr_flg",required_argument,0,0}, // [flg] Radius user-specified
    {"rfl_end",required_argument,0,0}, // [frc] Sample reflectance, diffuse incidence, end
    {"rfl_nbr",required_argument,0,0}, // [frc] Number of reflectances
    {"rfl_srt",required_argument,0,0}, // [frc] Sample reflectance, diffuse incidence, start
    {"rfl_usr_flg",required_argument,0,0}, // [flg] Reflectance user-specified
    {"sng_foo",required_argument,0,0}, // [sng] Intrinsic string temporary variable
    {"tst_sng",required_argument,0,0}, // [sng] Name of test to perform
    {"wvl_typ",required_argument,0,0}, // [enm] Wavelength type
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
      if(opt_crr == "azi_nbr") azi_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of azimuthal angles
      if(opt_crr == "azi_usr_flg"){azi_nbr=1L;flg_azi_usr=true;azi_dgr_usr=std::strtod(opt_sng.c_str(),(char **)NULL);} // [dgr] Azimuthal angle user-specified
      if(opt_crr == "dbl_foo") dbl_foo=std::strtod(opt_sng.c_str(),(char **)NULL); // [nbr] Intrinsic double temporary variable
      if(opt_crr == "dpt_nbr") dpt_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of layers
      if(opt_crr == "dpt_usr_flg"){dpt_nbr=1L;dpt_typ=dpt_sgl;flg_dpt_usr=true;dpt_usr=std::strtod(opt_sng.c_str(),(char **)NULL);} // [cm] Depth user-specified
      if(opt_crr == "drc_dff_fct") drc_dff_fct=std::strtod(opt_sng.c_str(),(char **)NULL); // [frc] Direct reflectance as fraction of diffuse
      if(opt_crr == "frc_msn_drc") frc_msn_drc=std::strtod(opt_sng.c_str(),(char **)NULL); // [frc] Fraction radiation emitted as direct, collimated light
      if(opt_crr == "flg_flg") flg_flg=true;
      if(opt_crr == "int_foo") int_foo=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Intrinsic int temporary variable
      if(opt_crr == "lng_foo") lng_foo=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Intrinsic long temporary variable
      if(opt_crr == "rds_nbr") rds_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of radii
      if(opt_crr == "rds_usr_flg"){rds_nbr=1L;flg_rds_usr=true;rds_usr=std::strtod(opt_sng.c_str(),(char **)NULL);} // [cm] Radius user-specified
      if(opt_crr == "rfl_end") rfl_end=std::strtod(opt_sng.c_str(),(char **)NULL); // Sample reflectance, diffuse incidence, ending value
      if(opt_crr == "rfl_nbr") rfl_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of reflectances
      if(opt_crr == "rfl_srt") rfl_srt=std::strtod(opt_sng.c_str(),(char **)NULL); // Sample reflectance, diffuse incidence, starting value
      if(opt_crr == "rfl_usr_flg"){rfl_nbr=1L;flg_rfl_usr=true;rfl_usr=std::strtod(opt_sng.c_str(),(char **)NULL);} // [frc] Reflectance user-specified
      if(opt_crr == "sng_foo") sng_foo=opt_sng; // [sng] Intrinsic string temporary variable
      if(opt_crr == "tst_sng") tst_sng=opt_sng; // [sng] Name of test to perform
      // Multi-line initializations break alphabetical-order-by-line rule
      if(opt_crr == "apx_typ"){
	// Decode approximation type
	if((apx_sng[apx_smp].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("apx_smp") != std::string::npos)){ 
	  apx_typ=apx_smp; // [enm] Approximation type
	}else if((apx_sng[apx_flt].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("apx_flt") != std::string::npos)){
	  apx_typ=apx_flt; // [enm] Approximation type
	}else if((apx_sng[apx_bff].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("apx_bff") != std::string::npos)){
	  apx_typ=apx_bff; // [enm] Approximation type
	}else{
	  std::cerr << "Unknown apx_typ" << std::endl;
	  return EXIT_FAILURE;
	} // end else
      } // end if "apx_typ"
      if(opt_crr == "dpt_typ"){
	// Decode depth type
	if((dpt_sng[dpt_rgl].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("dpt_rgl") != std::string::npos)){ 
	  dpt_typ=dpt_rgl; // [enm] Depth-grid type
	}else if((dpt_sng[dpt_sgl].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("dpt_sgl") != std::string::npos)){
	  dpt_typ=dpt_sgl; // [enm] Depth-grid type
	}else if((dpt_sng[dpt_GDZ09_13mm].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("dpt_GDZ09_13mm") != std::string::npos)){
	  dpt_typ=dpt_GDZ09_13mm; // [enm] Depth-grid type
	}else if((dpt_sng[dpt_GDZ09_25mm].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("dpt_GDZ09_25mm") != std::string::npos)){
	  dpt_typ=dpt_GDZ09_25mm; // [enm] Depth-grid type
	}else{
	  std::cerr << "Unknown dpt_typ" << std::endl;
	  return EXIT_FAILURE;
	} // end else
      } // end if "dpt_typ"
      if(opt_crr == "llm_typ"){
	// Decode illumination type
	if((llm_sng[llm_drc].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("llm_drc") != std::string::npos)){ 
	  llm_typ=llm_drc; // [enm] Illumination type
	}else if((llm_sng[llm_dff].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("llm_dff") != std::string::npos)){
	  llm_typ=llm_dff; // [enm] Illumination type
	}else{
	  std::cerr << "Unknown llm_typ" << std::endl;
	  return EXIT_FAILURE;
	} // end else
      } // end if "llm_typ"
      if(opt_crr == "wvl_typ"){
	// Decode wavelength type
	if((wvl_sng[wvl_0635].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("635") != std::string::npos)){ 
	  wvl_typ=wvl_0635; // [enm] Wavelength type
	}else if((wvl_sng[wvl_1310].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("1310") != std::string::npos)){
	  wvl_typ=wvl_1310; // [enm] Wavelength type
	}else if((wvl_sng[wvl_1550].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("1550") != std::string::npos)){
	  wvl_typ=wvl_1550; // [enm] Wavelength type
	}else{
	  std::cerr << "Unknown wvl_typ" << std::endl;
	  return EXIT_FAILURE;
	} // end else
      } // end if "wvl_typ"
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

  // Derived variables that depend on command-line input
  switch(wvl_typ){ // [enm] Wavelength type
  case wvl_0635: // [flg] Wavelength is 0.635 um
    rds_cll=rds_cll_0635nm; // [cm] Collimated (direct) beam radius (R_b)
    rfl_wll=rfl_wll_0635nm; // [frc] Wall reflectance
    break;
  case wvl_1310: // [flg] Wavelength is 1.310 um
    rds_cll=rds_cll_1310nm; // [cm] Collimated (direct) beam radius (R_b)
    rfl_wll=rfl_wll_1310nm; // [frc] Wall reflectance
    break;
  case wvl_1550: // [flg] Wavelength is 1.550 um
    rds_cll=rds_cll_1550nm; // [cm] Collimated (direct) beam radius (R_b)
    rfl_wll=rfl_wll_1550nm; // [frc] Wall reflectance
    break;
  default:
    std::cerr << "ERROR: Unknown wvl_typ" << std::endl;
    return EXIT_FAILURE;
  } // end wvl_typ switch 

  const double azi_dlt(azi_ntv/azi_nbr); // [frc] Azimuthal angle increment
  const double frc_msn_dff(1.0-frc_msn_drc); // [frc] Fraction radiation emitted as diffuse "stray light"
  const double rfl_dlt((rfl_end-rfl_srt)/rfl_nbr); // [frc] Reflectance increment

  if(dpt_typ == dpt_GDZ09_13mm) dpt_smp_hld=dpt_smp_hld_13mm; else dpt_smp_hld=dpt_smp_hld_25mm; // [cm] Depth integration interval
  // Scalar (uniform) depth increment needed only for regular grids
  const double dpt_ntv(dpt_smp_hld); // [cm] Depth integration interval
  switch(dpt_typ){ // [enm] Depth-grid type
  case dpt_rgl: // [flg] Regular (uniformly spaced) grid
    dpt_dlt=dpt_ntv/dpt_nbr; // [cm] Depth increment
    break;
  case dpt_sgl: // [flg] User-specified single depth
    // No-op
    break;
  case dpt_GDZ09_13mm: // [flg] Stretched 32-layer, 13mm grid from GDZ09
    dpt_nbr=sizeof(dpt_grd_GDZ09_13mm)/sizeof(double); // 32 [nbr] Number of layers
    break;
  case dpt_GDZ09_25mm: // [flg] Stretched 32-layer, 25mm grid from GDZ09
    dpt_nbr=sizeof(dpt_grd_GDZ09_25mm)/sizeof(double); // 32 [nbr] Number of layers
    break;
  default:
    std::cerr << "ERROR: Unknown dpt_typ" << std::endl;
    return EXIT_FAILURE;
  } // end dpt_typ switch 
  double *ngl_sld_avg_hms=new double[dpt_nbr]; // [hms] Mean solid angle subtended in sphere by layer

  // Radial integration limit is beam width and aperture radius for collimated 
  // (direct beam) and diffuse illumination, respectively
  switch(llm_typ){ // [enm] Illumination type
  case llm_drc: // [flg] Direct, collimated beam illumination
    rds_ntv=rds_llm=rds_cll; // [cm] Illuminated sample radius (R=R_b or R_a)
    break;
  case llm_dff: // [flg] Diffuse, isotropic illumination
    rds_ntv=rds_llm=rds_smp_apr; // [cm] Illuminated sample radius (R=R_b or R_a)
    break;
  default:
    std::cerr << "ERROR: Unknown llm_typ" << std::endl;
    return EXIT_FAILURE;
  } // end llm_typ switch 

  // Doubly-derived variables that depend on command-line input
  const double rds_dlt(rds_ntv/rds_nbr); // [cm] Radius increment

  // Finished with inputs and derived variables
  // Validate geometric configuration
  if(flg_rds_usr && rds_usr>rds_smp_apr){
    std::cout << "ERROR: rds_usr = " << rds_usr << " cm > rds_smp_apr = " << rds_smp_apr << " cm" << std::endl;
    return EXIT_FAILURE;
  } // endif err
  if(flg_rds_usr && rds_usr>rds_llm){
    std::cout << "ERROR: rds_usr = " << rds_usr << " cm > rds_llm = " << rds_llm << " cm" << std::endl;
    return EXIT_FAILURE;
  } // endif err
  if(flg_dpt_usr && dpt_usr>dpt_smp_hld){
    std::cout << "ERROR: dpt_usr = " << dpt_usr << " cm > dpt_smp_hld = " << dpt_smp_hld << " cm" << std::endl;
    return EXIT_FAILURE;
  } // endif err

  /* Nomenclature for alpha, d, h, m, r, s, w, follows HRL06:
     double alpha=sfc_wll_nrm; // [frc] Normalized wall area excluding holes HRL06 p. 5250
     double d=sfc_dtc_nrm; // [frc] Normalized detector (photodiode) area
     double h=sfc_src_nrm; // [frc] Normalized source (laser diode) area
     double m=sgn_msr_RS_AD; // [frc] Measured signal, direct incidence
     double r=rfl_smp; // [frc] Sample reflectance
     double s=sfc_smp_nrm; // [frc] Normalized sample area
     double w=rfl_wll; // [frc] Wall reflectance */

  /* Nomenclature for long-term probabilities (LTPs):
     prb_lng_RS_AD is LTP that photon reflected first from sample enters detector
     prb_lng_RW_AD is LTP that photon reflected first from wall   enters detector
     prb_lng_RW_AS is LTP that photon reflected first from wall is absorbed by sample
     prb_lng_RW_RS is LTP that photon reflected first from wall is reflected from sample
     Ghislain calls prb_lng_RW_RS the "back-intensity" and "re-illumination" */
  
  // Loop over snow reflectance
  for(rfl_idx=0;rfl_idx<rfl_nbr;rfl_idx++){
    if(flg_rfl_usr) rfl_smp=rfl_usr; else rfl_smp=rfl_srt+(0.5+rfl_idx)*rfl_dlt; // [frc] Sample reflectance, diffuse incidence

    /* Direct reflectance is factor times diffuse reflectance
       These are within ~5% of eachother for visible light in real snowpack
       Can be ~40% different for NIR light in real snowpack */
    rfl_smp_drc=drc_dff_fct*rfl_smp; // [frc] Sample reflectance, direct beam

    switch(apx_typ){ // [enm] Approximation type
    case apx_smp: // [flg] Integrating sphere as HRL06 simplest case
      // [frc] Diagonalization factor 1/(ad-bc) in Q matrix HRL06 p. 5250 (5)
      fct_mtx_nvr=1.0/(1.0-rfl_smp*sfc_smp_nrm-rfl_wll*sfc_wll_nrm);
      prb_RS_RW=rfl_wll*sfc_wll_nrm; // [frc] P(5,6) = Probability of RS->RW transition
      prb_lng_RS_AD=fct_mtx_nvr*sfc_dtc_nrm; // [frc] QR(1,1) = Probability of RS->AD transition
      prb_lng_RW_AD=fct_mtx_nvr*sfc_dtc_nrm; // [frc] QR(2,1) = Probability of RW->AD transition
      prb_lng_RW_AS=fct_mtx_nvr*sfc_smp_nrm*(1.0-rfl_smp); // [frc] QR(2,2) = Probability of RW->AS transition
      prb_lng_RW_RS=rfl_smp*sfc_smp_nrm*fct_mtx_nvr; // [frc] QT(2,1) = Probability of RW->RS transition
      break;
    case apx_flt: // [flg] Integrating sphere with flat sample
      fct_tmp=1.0-sfc_smp_nrm*(1.0-rfl_smp); // [frc] 1-s*(1-r)
      // [frc] Diagonalization factor 1/(ad-bc) in Q matrix HRL06 p. 5250 (5)
      fct_mtx_nvr=(1.0-sfc_smp_nrm)/(1.0-sfc_smp_nrm-rfl_wll*sfc_wll_nrm*fct_tmp);
      prb_RS_RW=rfl_wll*sfc_wll_nrm/(1.0-sfc_smp_nrm); // [frc] P(5,6) = Probability of RS->RW transition
      prb_lng_RS_AD=fct_mtx_nvr*sfc_dtc_nrm/(1.0-sfc_smp_nrm); // [frc] QR(1,1) = Probability of RS->AD transition
      prb_lng_RW_AD=fct_mtx_nvr*sfc_dtc_nrm*fct_tmp/(1.0-sfc_smp_nrm); // [frc] QR(2,1) = Probability of RW->AD transition
      prb_lng_RW_AS=fct_mtx_nvr*sfc_smp_nrm*(1.0-rfl_smp); // [frc] QR(2,2) = Probability of RW->AS transition
      prb_lng_RW_RS=fct_mtx_nvr*rfl_smp*sfc_smp_nrm; // [frc] QT(2,1) = Probability of RW->RS transition
      break;
    case apx_bff: // [flg] Integrating sphere with flat sample and baffle
      fct_tmp=1.0-(sfc_smp_nrm+sfc_dtc_nrm); // [frc] 1-(s+d)
      // [frc] Diagonalization factor 1/(ad-bc) in Q matrix HRL06 p. 5250 (5)
      fct_mtx_nvr=fct_tmp/(fct_tmp-rfl_wll*sfc_wll_nrm*(fct_tmp+rfl_smp*sfc_smp_nrm));
      prb_RS_RW=rfl_wll*sfc_wll_nrm/fct_tmp; // [frc] P(5,6) = Probability of RS->RW transition
      prb_lng_RS_AD=fct_mtx_nvr*sfc_dtc_nrm*rfl_wll*sfc_wll_nrm/fct_tmp; // [frc] QR(1,1) = Probability of RS->AD transition
      prb_lng_RW_AD=fct_mtx_nvr*sfc_dtc_nrm; // [frc] QR(2,1) = Probability of RW->AD transition
      prb_lng_RW_AS=fct_mtx_nvr*sfc_smp_nrm*(1.0-rfl_smp); // [frc] QR(2,2) = Probability of RW->AS transition
      prb_lng_RW_RS=fct_mtx_nvr*rfl_smp*sfc_smp_nrm; // [frc] QT(2,1) = Probability of RW->RS transition
      break;
    default:
      std::cerr << "ERROR: Unknown apx_typ" << std::endl;
      return EXIT_FAILURE;
    } // end apx_typ switch 
    
    prb_lng_RW_RS_or_RW_AS=prb_lng_RW_RS+prb_lng_RW_AS; // [frc] QT(2,1)+QR(2,2) = Union of RW->RS and RW->AS
    
    /* HRL06 calls   sgn_msr_RW_AD "measured spectrum, m" just after equation (6)
       HRL06 defines sgn_msr_RS_AD                        just after equation (7)
       The measured spectrum is useful for predicting measurements
       JC and FD do this calibration with the integrating sphere
       I usually try to predict the _corrections_ to the measurements
       HRL06 show how to invert sgn_msr for rfl_smp of unknown samples */
    sgn_msr_RS_AD=rfl_smp*prb_lng_RS_AD; // [frc] Measured signal from collimated illumination of sample
    sgn_msr_RW_AD=rfl_wll*prb_lng_RW_AD; // [frc] Measured signal from diffuse    illumination of sample

    /* Following implements ppr_ZGD09 equations eqn:frc_dff_dwn_RW_dfn, 
       eqn:frc_dff_dwn_flt_RW_dfn, and eqn:frc_dff_dwn_bff_RW_dfn 
       frc_bck_RS and frc_bck_RW are slightly ambiguous in manuscript
       Ambiguity arises from pre-factors that convert relative to absolute fractions */
    frc_bck_RS=frc_msn_drc*rfl_smp_drc*prb_RS_RW; // [frc] Fraction radiation backscattered into state RS
    frc_bck_RW=frc_msn_dff*rfl_wll; // [frc] Fraction radiation backscattered into state RW
    frc_dff_dwn=(frc_bck_RS+frc_bck_RW)*prb_lng_RW_RS_or_RW_AS; // [frc] Fraction diffuse downwelling radiation (f_dff/f_drc)   

    flx_frc_drc=1.0/(1.0+frc_dff_dwn); // [frc] Insolation fraction in direct beam ( = f_drc/(f_drc+f_dff) for swnb2)
    
    std::cout << "HRL06 rfl=" << rfl_smp << ": RW->AS = " << prb_lng_RW_AS << ", RW->RS = " << prb_lng_RW_RS << ", RW->(RS+AS) = " << prb_lng_RW_RS_or_RW_AS << ", frc_dff_dwn = " << frc_dff_dwn << ", flx_frc_drc = " << flx_frc_drc << std::endl;
    // Place output used by GDZ09 script on own line for simplicity
    std::cout << "flx_frc_drc_GDZ09 = " << flx_frc_drc << std::endl;

  } // end loop over sample reflectance

  /* Mathematica gives analytic solution to indefinite integrals of
     psi_1="x ArcTan[a/(b-x)]" and psi_2="x ArcTan[a/(b+x)]" at
     http://integrals.wolfram.com/index.jsp
     Scribblings in my LGGE notebook simplify their sum = psi_1+psi_2 
     psi_1="ArcTan[a/(Sqrt[b-Cos[x]Cos[x]]-Sin[x])" has no solution
     Notation is that 
     x = r     = Radius from center of sample
     a = z+z_s = Depth in snow sample (plus any lip thickness)
     b = R_a   = Radius of aperature to sample holder */

  // Loop over depth
  for(dpt_idx=0;dpt_idx<dpt_nbr;dpt_idx++){
    switch(dpt_typ){ // [enm] Depth-grid type
    case dpt_rgl: // [flg] Regular (uniformly spaced) grid
      dpt=dpt_srt+(0.5+dpt_idx)*dpt_dlt; // [cm] Depth (centered)
      break;
    case dpt_sgl: // [flg] User-specified single depth
      dpt=dpt_usr; // [cm] Depth (centered)
      break;
    case dpt_GDZ09_13mm: // [flg] Stretched 32-layer, 13mm grid from GDZ09
      dpt=dpt_grd_GDZ09_13mm[dpt_idx]; // [cm] Depth (centered)
      break;
    case dpt_GDZ09_25mm: // [flg] Stretched 32-layer, 25mm grid from GDZ09
      dpt=dpt_grd_GDZ09_25mm[dpt_idx]; // [cm] Depth (centered)
      break;
    default:
      std::cerr << "ERROR: Unknown dpt_typ" << std::endl;
      return EXIT_FAILURE;
    } // end dpt_typ switch 
    
    // Initialize integrands
    ngl_pln_dbl_avg=0.0; // [rdn] Mean planar angle subtended in sphere by layer
    ngl_sld_avg=0.0; // [sr] Mean solid angle subtended in sphere by layer

    dpt_snw_pls_lip=dpt+hgt_lip; // [cm] Depth in snow sample plus lip thickness
    ntg_dfn_val=fnc_ngl_ntg_ndf_azi_90(dpt_snw_pls_lip,rds_smp_apr,rds_llm)-fnc_ngl_ntg_ndf_azi_90(dpt_snw_pls_lip,rds_smp_apr,0.0);
    ngl_pln_dbl_avg_azi_90=cst_M_PIl-2.0*ntg_dfn_val/(rds_llm*rds_llm); // [rdn] Mean planar angle subtended in sphere by layer for azimuth=90 degrees
    ngl_sld_avg_azi_90=2.0*cst_M_PIl*(1.0-cos(ngl_pln_dbl_avg_azi_90/2.0)); // [sr] Mean solid angle subtended in sphere by layer for azimuth=90 degrees
    
    // Loop over azimuthal radius
    for(rds_idx=0;rds_idx<rds_nbr;rds_idx++){
      if(flg_rds_usr) rds=rds_usr; else rds=rds_srt+(0.5+rds_idx)*rds_dlt; // [cm] Radius (centered)

      // Initialize integrand
      ngl_pln_avg=0.0; // [rdn] Mean planar angle subtended by point
      
      // Loop over azimuthal angle
      for(azi_idx=0;azi_idx<azi_nbr;azi_idx++){
	if(flg_azi_usr) azi=azi_dgr_usr*cst_M_PIl/180.0; else azi=azi_srt+(0.5+azi_idx)*azi_dlt; // [rdn] Azimuthal angle (centered)
	
	rcd+=fnc_ngl_psi_one_psi_two // [fnc] Azimuthal angles psi_1 and psi_2
	  (rds, // I [cm] Radius from center of sample
	   azi, // I [rdn] Azimuthal angle
	   dpt_snw_pls_lip, // I [cm] Depth in snow sample plus lip thickness
	   rds_smp_apr, // I [cm] Radius of aperature to sample holder
	   ngl_apr_psi_1, // O [rdn] Supplementary planar angle one subtended by point
	   ngl_apr_psi_2); // O [rdn] Supplementary planar angle two subtended by point
	
	ngl_apr_pln=cst_M_PIl-(ngl_apr_psi_1+ngl_apr_psi_2); // [rdn] Planar angle subtended by point
	ngl_pln_avg+=azi_dlt*ngl_apr_pln; // [rdn] Mean planar angle subtended by point
	if(dbg_lvl == 4){
	  std::cout << "dpt[" << dpt_idx << "]=" << dpt << " cm, rds[" << rds_idx << "]=" << rds << " cm, azi[" << azi_idx << "]=" << azi*180.0/cst_M_PIl << " dgr: psi = " << ngl_apr_pln*180.0/cst_M_PIl << " dgr, psi_1 = " << ngl_apr_psi_1*180.0/cst_M_PIl << " dgr, psi_2 = " << ngl_apr_psi_2*180.0/cst_M_PIl << " dgr" << std::endl;
	} // endif dbg

      } // end loop over azimuth

      // Normalize integrand by integration interval
      ngl_pln_avg/=azi_ntv; // [rdn] Mean planar angle subtended by point
      ngl_apr_sld=2.0*cst_M_PIl*(1.0-cos(ngl_pln_avg/2.0)); // [sr] Solid angle subtended from point to aperature
      
      if(dbg_lvl >= 3){
	std::ios::fmtflags ios_flg_dfl(std::cout.flags()); // [flg] I/O stream defaults
	std::cout.setf(std::ios::fixed); // [flg] Fixed number of digits to right of decimal
	std::cout << "dpt[" << dpt_idx << "]=" << dpt << " cm, rds[" << rds_idx << "]=" << rds << " cm: psi_bar=" << std::setprecision(8) << ngl_pln_avg*180.0/cst_M_PIl << " dgr, omega=" << ngl_apr_sld/(2.0*cst_M_PIl) << " hms" << std::endl;
	// Restore defaults
	std::cout.precision(3); // Number of digits to right of decimal
	std::cout.flags(ios_flg_dfl); // [srm] Output stream
      } // endif dbg

      ngl_pln_dbl_avg+=rds*ngl_pln_avg*rds_dlt; // [rdn] Mean planar angle subtended in sphere by layer
    
      ngl_sld_avg+=rds*cos(ngl_pln_avg/2.0)*rds_dlt; // [sr] Mean solid angle subtended in sphere by layer

    } // end loop over radius
    
    // Normalize integrand (azimuthal integration factor azi_ntv already applied)
    ngl_pln_dbl_avg=ngl_pln_dbl_avg*2.0/(rds_ntv*rds_ntv); // [rdn] Mean planar angle subtended in sphere by layer

    ngl_sld_avg=2.0*cst_M_PIl*(1.0-2.0*ngl_sld_avg/(rds_ntv*rds_ntv)); // [sr] Mean solid angle subtended in sphere by layer

    // Archive for easier output diagnostics
    ngl_sld_avg_hms[dpt_idx]=ngl_sld_avg/(2.0*cst_M_PIl); // [hms] Mean solid angle subtended in sphere by layer
	
    std::cout << "dpt[" << dpt_idx << "]=" << dpt << " cm: psi_dbl_bar(azi=90) = " << ngl_pln_dbl_avg_azi_90*180.0/cst_M_PIl << " dgr, omega_bar(azi=90) = " << ngl_sld_avg_azi_90/(2.0*cst_M_PIl) << " hms" << std::endl;
    
    std::cout << "dpt[" << dpt_idx << "]=" << dpt << " cm: psi_dbl_bar         = " << ngl_pln_dbl_avg       *180.0/cst_M_PIl << " dgr, omega_bar         = " << ngl_sld_avg       /(2.0*cst_M_PIl) << " hms" << std::endl;
    
  } // end loop over snow depth

  if(flg_rds_usr && rds_usr==0.0) std::cout << "WARNING: Setting --rds_usr=0.0 also zeros area averages because of \"r\" weight. psi_dbl_bar will be zero, though psi and psi_bar are valid." << std::endl;
  if(rds_nbr==1) std::cout << "WARNING: rds_nbr=1 so do not trust areal averages" << std::endl;
  if(azi_nbr==1) std::cout << "WARNING: azi_nbr=1 so do not trust azimuthal averages" << std::endl;

  if(dbg_lvl == 1){
    std::cout << "Initialization State:" << std::endl;
    std::cout << "Approximation type: " << apx_sng[apx_typ] << std::endl;
    std::cout << "Depth-grid    type: " << dpt_sng[dpt_typ] << std::endl;
    std::cout << "Illumination  type: " << llm_sng[llm_typ] << std::endl;
    std::cout << "Wavelength        : " << wvl_sng[wvl_typ] << std::endl;
    std::cout << "dpt_smp_hld=" << dpt_smp_hld << " [cm] Depth of sample holder" << std::endl;
    std::cout << "drc_dff_fct=" << drc_dff_fct << " [frc] Direct reflectance as fraction of diffuse" << std::endl;
    std::cout << "fct_mtx_nvr=" << fct_mtx_nvr << " [frc] Matrix diagonalization factor" << std::endl;
    std::cout << "frc_msn_dff=" << frc_msn_dff << " [frc] Fraction radiation emitted as diffuse \"stray light\"" << std::endl;
    std::cout << "frc_msn_drc=" << frc_msn_drc << " [frc] Fraction radiation emitted as direct, collimated light" << std::endl;
    std::cout << "hgt_lip=" << hgt_lip << " [cm] Thickness of lip from sphere to sample top (z_s)" << std::endl;
    std::cout << "hgt_smp=" << hgt_smp << " [cm] Height of sample surface above tangent plane" << std::endl;
    std::cout << "rds_cll=" << rds_cll << " [cm] Collimated (direct) beam radius (R_b)" << std::endl;
    std::cout << "rds_dtc=" << rds_dtc << " [cm] Detector (photodiode) radius (R_d)" << std::endl;
    std::cout << "rds_llm=" << rds_llm << " [cm] Illuminated sample radius (R=R_b or R_a)" << std::endl;
    std::cout << "rds_smp_apr=" << rds_smp_apr << " [cm] Aperture (exposed sample) radius (R_a)" << std::endl;
    std::cout << "rds_smp_hld=" << rds_smp_hld << " [cm] Radius of sample holder (R_E)" << std::endl;
    std::cout << "rds_sph_ntg=" << rds_sph_ntg << " [cm] Integrating sphere radius (R_IS)" << std::endl;
    std::cout << "rds_src=" << rds_src << " [cm] Source (laser diode) radius (R_h)" << std::endl;
    std::cout << "rfl_smp_dff=" << rfl_smp << " [frc] Sample diffuse reflectance \"r\" (loop exit value)" << std::endl;
    std::cout << "rfl_smp_drc=" << rfl_smp_drc << " [frc] Sample direct reflectance \"r^.\" (loop exit value)" << std::endl;
    std::cout << "rfl_wll=" << rfl_wll << " [frc] Wall reflectance \"w\"" << std::endl;
    std::cout << "rfl_wll_0635nm=" << rfl_wll_0635nm << " [frc] Wall reflectance at 0635 nm" << std::endl;
    std::cout << "rfl_wll_1310nm=" << rfl_wll_1310nm << " [frc] Wall reflectance at 1310 nm" << std::endl;
    std::cout << "rfl_wll_1550nm=" << rfl_wll_1550nm << " [frc] Wall reflectance at 1550 nm" << std::endl;
    std::cout << "sfc_dtc_nrm=" << sfc_dtc_nrm << " [frc] Normalized detector (photodiode) area \"d\"" << std::endl;
    std::cout << "sfc_smp_nrm=" << sfc_smp_nrm << " [frc] Normalized sample area \"s\"" << std::endl;
    std::cout << "sfc_src_nrm=" << sfc_src_nrm << " [frc] Normalized source (laser diode) area \"h\"" << std::endl;
    std::cout << "sfc_wll_nrm=" << sfc_wll_nrm << " [frc] Normalized wall area excluding holes \"alpha\"" << std::endl;

    // Diagnostic FOV in easy ncap2 format
    std::cout << "ngl_sld_avg_hms[levp_snw]={";
    for(dpt_idx=0;dpt_idx<dpt_nbr;dpt_idx++) std::cout << ngl_sld_avg_hms[dpt_idx] << ",";
    // Repeat last layer value for bottom reflectance
    std::cout << ngl_sld_avg_hms[dpt_nbr-1] << "};" << std::endl;

  } // endif dbg
  
  // Memory clean-up
  delete []ngl_sld_avg_hms; // [hms] Mean solid angle subtended in sphere by layer

} // end snw

int // O [enm] Return success code
fnc_ngl_psi_one_psi_two // [fnc] Azimuthal angles psi_1 and psi_2
(const double rds, // I [cm] Radius from center of sample
 const double azi, // I [rdn] Azimuthal angle
 const double dpt, // I [cm] Depth in snow sample plus lip thickness
 const double rds_smp_apr, // I [cm] Radius of aperature to sample holder
 double &ngl_apr_psi_1, // O [rdn] Supplementary planar angle one subtended by point
 double &ngl_apr_psi_2) // O [rdn] Supplementary planar angle two subtended by point
{
  // Purpose: Compute psi_1
  int rcd(0); // [enm] Return success code
  double trm_1;
  double trm_2;
  double trm_3;
  trm_1=rds*cos(azi);
  trm_2=rds*sin(azi);
  trm_3=sqrt(rds_smp_apr*rds_smp_apr-trm_1*trm_1);
  ngl_apr_psi_1=atan(dpt/(trm_3-trm_2));
  ngl_apr_psi_2=atan(dpt/(trm_3+trm_2));
  return rcd; // [enm] Return success code
} // end fnc_ngl_psi_one()

double 
fnc_ngl_ntg_ndf_azi_90 // [fnc] Indefinite integral for psi_1+psi_2 for azi=90
(const double a, // a = z+z_s = Depth in snow sample plus lip thickness
 const double b, // b = R_a   = Radius of aperature to sample holder
 const double x) // x = r     = Radius from center of sample
{
  // Purpose: Indefinite integral for psi_1+psi_2
  const double cst_M_PIl(3.1415926535897932384626433832795029L); // [frc] 3
  double trm_1;
  double trm_2;
  double trm_3;
  trm_1=(x*x/2.0)*(atan(a/(b+x))+atan(a/(b-x)));
  /* Avoid divide-by-zero's in exact expression
     trm_2=((b*b-a*a)/2.0)*(atan((x+b)/a)-atan((x-b)/a)); */
  trm_2=(b*b-a*a)/2.0;
  if(a==0.0) trm_2*=cst_M_PIl; else trm_2*=atan((x+b)/a)-atan((x-b)/a);
  trm_3=(a*b/2.0)*log((a*a+b*b+x*x+2*b*x)*(a*a+b*b+x*x-2*b*x));
  return trm_1+trm_2-trm_3;
} // end fnc_ngl_ntg_ndf_azi_90()
