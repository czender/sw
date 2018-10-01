/* $Id$ */

/* htrn -- Parse and compute HITRAN statistics */

/* Copyright (C) 1994--2018 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text
   The original author of this software, Charlie Zender, seeks to improve
   it with your suggestions, contributions, bug-reports, and patches.
   Charlie Zender <zender at uci dot edu>
   Department of Earth System Science
   University of California, Irvine
   Irvine, CA 92697-3100 */

/* fxm: distinguish between isotope, isotopomer, and isotopologue:
   14C is an isotope of carbon, as is 13C
   Thus carbon has 3 abundant isotopes, 12C, 13C, 14C
   All molecules are isotopomers
   There are about 90 isotopomers in HITRAN, but no single species has more than about 10 isotopes
   Should switch to itp = isotope, ipm = isotopomer 
   RRG98 use iso = isotope 
   fxm: remove datestamp in output files? */

/* Usage:
   Compilation:
   cd ~/ck;make htrn

   Debugging:
   htrn --dbg=2 --pth_out=${HOME}/sw/ck --pth_out_F=${HOME}/sw/ck --pth_out_F90=${HOME}/sw/ck --pth_out_cc=${HOME}/sw/ck # borken

   Production:
   htrn
   htrn --dbg=1
   htrn --dbg=2

   By default, htrn overwrites following files:
   ${HOME}/sw/ck/HITRAN.pm
   ${HOME}/sw/f/htrn_mdl.F90
   ${HOME}/sw/aca/hitran.F
   ${HOME}/sw/aca/hitran.com
   ${HOME}/sw/ck/hitran.h
   ${HOME}/sw/c++/htrn_c++.cc
   ${HOME}/sw/c++/htrn_c++.hh
   To change this, specify different output directories, e.g., home directory with "sw/ck/" prefix as above */

/* Purpose: Parse HITRAN "molparam.txt" molecule parameter file and auotmatically
   generate header files for C, Fortran, and Perl programs.

   Contents of "molparam.txt" are expected to look like:
Molecule # Iso Abundance     Q(296K)      gj    Molar Mass(g)
   H2O (1)
         161  .997317E+00   .174626E+03    1     18.010565
         181  .199983E-02   .176141E+03    1     20.014811 */

/* History:
   20010212: Filled in Q and gj columns of molparam.txt for HITRAN2000
   H2O isotopes 5 and 6 (which were blank and thus crashing htrn)
   20090512: Renamed input molparamYY.txt. Tried with HITRAN2008.
   20130611: Updated to HITRAN2012
   20181001: Updated to HITRAN2016 */

/* Standard C headers */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <sys/stat.h> /* stat() */
#include <time.h> /* machine time */
#include <unistd.h> /* all sorts of POSIX stuff */

#ifndef LINUX
#include "getopt.h" /* GNU getopt() */
#else /* LINUX */
#include <getopt.h> /* GNU getopt() is standard on Linux */
#endif /* LINUX */

/* 3rd party vendors */

/* Personal headers */
#include <dbg.h> /* Debugging constants */

/* #define MAIN_PROGRAM_FILE MUST precede #include nc.h */
#define MAIN_PROGRAM_FILE

#define MLC_NBR_MAX_HTRN 50 /* Maximum number of molecules in HITRAN */
#define ISO_PER_MLC_NBR_MAX_HTRN 20 /* Maximum number of isotopes per molecule in HITRAN */
#define ISO_NBR_MAX_HTRN 150 /* Maximum number of isotopomers in HITRAN */

enum mlc_grp{
  unknown_structure, /* Molecule has not been assigned a group */
  asymmetric_rotor, /* Group 1. Asymmetric Rotors */
  diatomic_and_linear_molecule_with_integer_J, /* Group 2. Diatomic and Linear Molecules with Integer J */
  spherical_rotor, /* Group 3. Spherical Rotors */
  symmetric_rotor, /* Group 4. Symmetric Rotors */
  triplet_ground_electronic_state, /* Group 5. Triplet Ground Electronic States */
  doublet_ground_electronic_state_with_half_integer_J /* Group 6. Doublet Ground Electronic States (Half Integer J) */
};
char *mlc_grp_sng[7]={
  "Unknown group",
  "Asymmetric Rotors",
  "Diatomic and Linear Molecules with Integer J",
  "Spherical Rotors",
  "Symmetric Rotors",
  "Triplet Ground Electronic States ",
  "Doublet Ground Electronic States (Half Integer J)"
};

typedef struct{
  char *nm; /* HITRAN name. Currently isotopic abbreviation */
  char *nm_sf; /* Name used in C and Fortran programs (HITRAN name stripped of illegal characters like '+') */
  char *iso_abb; /* HITRAN isotopic abbreviation. */
  double TIPS; /* Total internal partition sum (i.e., partition function) at 296K */
  double cnc_frc; /* Fractional abundance of isotope in Earth's atmosphere */
  double mmw; /* Mean molecular weight of isotope */
  short gj; /* Spin statistical weight */
  short iso_idx; /* HITRAN isotope index. 1 = most abundant in atmosphere, 2 = next most abundant... */
  short iso_cnsc_idx; /* HITRAN consecutive isotope index. 1 = 1H2_16O, ... 5 = 12C_16O2... */
  short mlc_idx; /* HITRAN molecule index. 1 = H2O, ... */
  short nm_len; /* Length of nm_sf */
} iso_sct;

typedef struct{
  char *nm; /* HITRAN name */
  char *nm_sf; /* Name used in C and Fortran programs (HITRAN name stripped of illegal characters like '+') */
  double mmw; /* Mean mean molecular weight of gas obtained by weighting isotopic mmw by isotopic concentration */
  iso_sct **iso; /* Pointer to list of isotopes */
  float rtt_fnc_tpt_xpn; /* Rotational partition function temperature dependent exponent */
  short iso_nbr; /* Number of isotopes */
  short mlc_idx; /* HITRAN molecule index. 1 = H2O, ... */
  short nm_len; /* Length of nm_sf */
  short mlc_grp; /* Structural class */
} mlc_sct;

/* Global variables declared here */
char *prg_nm;

unsigned short int dbg_lvl=0; /* Option D */
unsigned short int dbg_val=0; /* Option d */

int main(int argc,char **argv)
{
  char *nm_mk_sf(char *);
  char *nm_mk_iso(char *,short);
  char *cmd_ln_sng(int,char **);
  void mlc_grp_set(mlc_sct *);
  void rtt_fnc_tpt_xpn_set(mlc_sct *);
  void Exit_gracefully(void);
  void usg_prn(char *);
  
  char mmw_mlc_get_fnc_prt[1000];
  char mmw_iso_get_fnc_prt[1000];
  char mlc_sng_get_fnc_prt[1000];
  char iso_sng_get_fnc_prt[1000];
  char rtl_fnc_tpt_xpn_get_fnc_prt[1000];
  char iso_idx_map_fnc_prt[1000];

  FILE *fp_in;
  FILE *fp_htrn_F90;
  FILE *fp_htrn_F;
  FILE *fp_htrn_com;
  FILE *fp_htrn_h;
  FILE *fp_htrn_cc;
  FILE *fp_htrn_hh;
  FILE *fp_htrn_pm;
  
  char *cmd_ln;
  char *fl_in=NULL;
  char *fl_out=NULL;
  char *fl_htrn_F90=NULL;
  char *fl_htrn_F=NULL;
  char *fl_htrn_com=NULL;
  char *fl_htrn_hh=NULL;
  char *fl_htrn_cc=NULL;
  char *fl_htrn_h=NULL;
  char *fl_htrn_pm=NULL;
  char mlc_sng[8];
  char iso_sng[8];
  char *time_bfr_srt;
  char *pth_in=NULL; /* [sng] Input directory */
  char *pth_out=NULL; /* [sng] Output directory */
  char *pth_out_F=NULL; /* [sng] Output directory for Fortran files */
  char *pth_out_F90=NULL; /* [sng] Output directory for Fortran90 files */
  char *pth_out_cc=NULL; /* [sng] Output directory for C++ files */
  char *home_nvr=NULL;
  char bfr[1000]; /* Space for scanf() input */
  char CVS_Id[]="$Id$\n";
  char CVS_Revision[]="$Revision$\n";
  char sbr_nm[]="main()";
  
  double cnc_frc_ttl;
  double cnc_frc;
  double TIPS;
  double mmw;
  short gj;
  
  double dbl_foo=0.0; /* Option dbl_foo */
  float flt_foo=0.0; /* Option f or flt_foo */
  int rcd=0; /* [rcd] Return code */
  long pth_out_sng_lng; /* [nbr] Length of pth_out string */
  short sht_foo=73; /* Option sht_foo */
  
  iso_sct **iso;
  mlc_sct **mlc;
  
  short idx;
  short iso_idx=0; /* CEWI */
  short iso_cnsc_idx;
  short iso_nbr;
  short max_nbr_iso_per_mlc;
  short mlc_idx;
  short mlc_idx_htrn;
  short mlc_nbr;
  
  /* Variables for option processing */
  char *opt_crr; /* String representation of current long-option name */
  char *opt_short_lst; /* String representation of current optarg, if any */
  char *opt_sng=NULL; /* String representation of current optarg, if any */ /* CEWI */
  extern char *optarg; /* char * representation of current optarg, if any (this memory is owned by system) */
  extern int optind; /* extern enumerating cardinal of current option */
  int opt; /* Value is zero if current argument is long type, else value contains single letter version of command line argument */
  int opt_idx=0; /* Index of current long option into opt_lng array */
  static struct option opt_lng[] =
  {
    /* Option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is compared to enum _argtype{no_argument,required_argument,optional_argument}, 
       flag points to a variable that gets set to val whenever the name option is set.
       For long options that have zero flag field, getopt() returns contents of val
    */
    /* Long options with no argument */
    /* Long options with argument */
    {"dbl_foo",1,0,0}, /* [frc] Intrinsic double temporary variable */
    {"sht_foo",1,0,0}, /* [frc] Intrinsic short temporary variable */
    {"pth_in",1,0,0}, /* [sng] Input directory */
    {"pth_out",1,0,0}, /* [sng] Output directory */
    {"pth_out_F",1,0,0}, /* [sng] Output directory for Fortran files*/
    {"pth_out_F90",1,0,0}, /* [sng] Output directory for Fortran90 files*/
    {"pth_out_cc",1,0,0}, /* [sng] Output directory for C++ files*/
    /* Long options with optional argument */
    /* Long options with short counterparts */
    {"dbg_lvl",1,0,'D'}, /* [enm] Debugging level */
    {"debug",1,0,'D'},
    {"input",1,0,'i'},
    {"output",1,0,'o'},
    {"version",0,0,'v'},
    {"verbose",0,0,'D'},
    /* Last option named "0" to signal getopt_long() to stop processing */
    {0,0,0,0}
  }; /* end opt_lng */
  
  time_t clock;
  
  /* Set defaults */
  
  /* Default file names */
  home_nvr=getenv("HOME");
  pth_in=(char *)strdup(home_nvr); /* [sng] Input directory */
  if(pth_in != NULL){
    /* Add one for NUL-terminator */
    char pth_sfx_in[]="/sw/ck/"; /* [sng] */
    pth_in=realloc(pth_in,strlen(pth_sfx_in)+strlen(pth_in)+1L);
    pth_in=strcat(pth_in,(char *)strdup(pth_sfx_in));
  } /* end if */

  pth_out=(char *)strdup(home_nvr); /* [sng] Output directory */
  pth_out_sng_lng=strlen(pth_out); /* [nbr] Length of pth_out string */
  if(pth_out != NULL){
    /* Add one for NUL-terminator */
    char pth_sfx_out[]="/sw/ck/"; /* [sng] */
    pth_out=realloc(pth_out,strlen(pth_sfx_out)+strlen(pth_out)+1L);
    pth_out=strcat(pth_out,(char *)strdup(pth_sfx_out));
  } /* end if */

  pth_out_F=(char *)strdup(home_nvr); /* [sng] Output directory for Fortran files */
  if(pth_out_F != NULL){
    /* Add one for NUL-terminator */
    char pth_sfx_F[]="/sw/aca/"; /* [sng] */
    pth_out_F=realloc(pth_out_F,strlen(pth_sfx_F)+strlen(pth_out_F)+1L);
    pth_out_F=strcat(pth_out_F,(char *)strdup(pth_sfx_F));
  } /* end if */

  pth_out_F90=(char *)strdup(home_nvr); /* [sng] Output directory for Fortran90 files */
  if(pth_out_F90 != NULL){
    /* Add one for NUL-terminator */
    char pth_sfx_F90[]="/sw/f/"; /* [sng] */
    pth_out_F90=realloc(pth_out_F90,strlen(pth_sfx_F90)+strlen(pth_out_F90)+1L);
    pth_out_F90=strcat(pth_out_F90,(char *)strdup(pth_sfx_F90));
  } /* end if */

  pth_out_cc=(char *)strdup(home_nvr); /* [sng] Output directory for C++ files */
  if(pth_out_cc != NULL){
    /* Add one for NUL-terminator */
    char pth_sfx_cc[]="/sw/c++/"; /* [sng] */
    pth_out_cc=realloc(pth_out_cc,strlen(pth_sfx_cc)+strlen(pth_out_cc)+1L);
    pth_out_cc=strcat(pth_out_cc,(char *)strdup(pth_sfx_cc));
  } /* end if */

  /* Start the clock and save the command line */
  cmd_ln=cmd_ln_sng(argc,argv);
  clock=time((time_t *)NULL);
  time_bfr_srt=ctime(&clock);
  (void)fprintf(stderr,"\tStart = %s",time_bfr_srt);
  prg_nm=((prg_nm=strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++;
  
  /* Short options: no colon = no arg, one colon = required arg, two colons = optional arg */
  opt_short_lst="D:f:i:o:v"; /* List of single-letter (C-style) option abbreviations */
  /* Parse command line arguments */
  while(1){
    /* getopt_long_only() allows a single dash '-' to prefix long options as well */
    opt=getopt_long_only(argc,argv,opt_short_lst,opt_lng,&opt_idx);
    /* NB: access to opt_crr is only valid when long_opt was detected */
    opt_crr=(char *)strdup(opt_lng[opt_idx].name);
    if(optarg) opt_sng=optarg; /* Copy system memory into local string for safer operations */
    if(opt == EOF) break; /* Parse positional arguments once getopt_long_only() returns EOF */
    /* Process long options without short option counterparts */
    if(opt == 0){
      if(dbg_lvl >= dbg_io){
	(void)fprintf(stdout,"Long option name: %s, ",opt_lng[opt_idx].name);
	if(optarg) (void)fprintf(stdout," Argument: %s\n",opt_sng); else (void)fprintf(stdout," No argument\n");
      } /* end if dbg */
      if(!strcmp(opt_crr,"pth_in")) pth_in=(char *)strdup(opt_sng); /* [sng] Input directory */
      if(!strcmp(opt_crr,"pth_out")){
	pth_out=(char *)strdup(opt_sng); /* [sng] Output directory */
	if(pth_out[pth_out_sng_lng-1L] != '/'){
	  pth_out=(char *)realloc(pth_out,(pth_out_sng_lng+2L)*sizeof(char));
	  pth_out[pth_out_sng_lng]='/';
	  pth_out[pth_out_sng_lng+1L]='\0';
	} /* endif */
      } /* endif */
      if(!strcmp(opt_crr,"pth_out_F")){
	pth_out_F=(char *)strdup(opt_sng); /* [sng] Output directory for Fortran files */
	if(pth_out_F[pth_out_sng_lng-1L] != '/'){
	  pth_out_F=(char *)realloc(pth_out_F,(pth_out_sng_lng+2L)*sizeof(char));
	  pth_out_F[pth_out_sng_lng]='/';
	  pth_out_F[pth_out_sng_lng+1L]='\0';
	} /* endif */
      } /* endif */
      if(!strcmp(opt_crr,"pth_out_F90")){
	pth_out_F90=(char *)strdup(opt_sng); /* [sng] Output directory for Fortran files */
	if(pth_out_F90[pth_out_sng_lng-1L] != '/'){
	  pth_out_F90=(char *)realloc(pth_out_F90,(pth_out_sng_lng+2L)*sizeof(char));
	  pth_out_F90[pth_out_sng_lng]='/';
	  pth_out_F90[pth_out_sng_lng+1L]='\0';
	} /* endif */
      } /* endif */
      if(!strcmp(opt_crr,"pth_out_cc")){
	pth_out_cc=(char *)strdup(opt_sng); /* [sng] Output directory for C++ files */
	if(pth_out_cc[pth_out_sng_lng-1L] != '/'){
	  pth_out_cc=(char *)realloc(pth_out_cc,(pth_out_sng_lng+2L)*sizeof(char));
	  pth_out_cc[pth_out_sng_lng]='/';
	  pth_out_cc[pth_out_sng_lng+1L]='\0';
	} /* endif */
      } /* endif */
      if(!strcmp(opt_crr,"dbl_foo")) dbl_foo=strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"flt_foo")) flt_foo=(float)strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"sht_foo")) sht_foo=(short)strtol(opt_sng,(char **)NULL,10);
    } /* opt != 0 */
    switch(opt){
    case 0: /* Long options have already been processed, return */
      break;
    case 'D': /* Debugging level.  Default is 0. */
      if(optarg) dbg_lvl=(unsigned short int)strtol(optarg,(char **)NULL,10); else dbg_lvl=dbg_fl;
      break;
    case 'f': /* Generic tuning parameter.  Default is 0. */
      flt_foo=(float)strtod(optarg,(char **)NULL);
      break;
    case 'i': /* Input file name. Default is stdin */
      fl_in=opt_sng;
      break;
    case 'o': /* Output file name. Default is stdout */
      fl_out=opt_sng;
      break;
    case 'v': /* CVS program info */
      (void)fprintf(stderr,"%s %s\n",CVS_Revision,CVS_Id);
      exit(EXIT_SUCCESS);
      break;
    default: /* Print proper usage */
      (void)usg_prn(opt_short_lst);
      exit(EXIT_FAILURE);
    } /* end switch */
  } /* end while loop */
  
  /* Process positional arguments */
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      (void)fprintf(stdout,"%s: ERROR %s reports too many positional arguments\n",prg_nm,sbr_nm);
      exit(EXIT_FAILURE);
    }else if(psn_arg_nbr == 1){
      fl_out=argv[optind++];
    }else if(psn_arg_nbr == 2){
      fl_in=argv[optind++];
      fl_out=argv[optind++];
    } /* end else */
  } /* end if */

  fl_in=(char *)strdup(pth_in);
  fl_in=(char *)realloc(fl_in,(strlen(pth_in)+strlen("molparam12.txt")+1L)*sizeof(char));
  (void)strcat(fl_in,"molparam16.txt");
  fl_htrn_com=(char *)strdup(pth_out_F);
  fl_htrn_com=(char *)realloc(fl_htrn_com,(strlen(pth_out_F)+strlen("hitran.com")+1L)*sizeof(char));
  (void)strcat(fl_htrn_com,"hitran.com");
  fl_htrn_F=(char *)strdup(pth_out_F);
  fl_htrn_F=(char *)realloc(fl_htrn_F,(strlen(pth_out_F)+strlen("hitran.F")+1L)*sizeof(char));
  (void)strcat(fl_htrn_F,"hitran.F");
  fl_htrn_F90=(char *)strdup(pth_out_F90);
  fl_htrn_F90=(char *)realloc(fl_htrn_F90,(strlen(pth_out_F90)+strlen("htrn_mdl.F90")+1L)*sizeof(char));
  (void)strcat(fl_htrn_F90,"htrn_mdl.F90");
  fl_htrn_hh=(char *)strdup(pth_out_cc);
  fl_htrn_hh=(char *)realloc(fl_htrn_hh,(strlen(pth_out_cc)+strlen("htrn_c++.hh")+1L)*sizeof(char));
  (void)strcat(fl_htrn_hh,"htrn_c++.hh");
  fl_htrn_cc=(char *)strdup(pth_out_cc);
  fl_htrn_cc=(char *)realloc(fl_htrn_cc,(strlen(pth_out_cc)+strlen("htrn_c++.cc")+1L)*sizeof(char));
  (void)strcat(fl_htrn_cc,"htrn_c++.cc");
  fl_htrn_h=(char *)strdup(pth_out);
  fl_htrn_h=(char *)realloc(fl_htrn_h,(strlen(pth_out)+strlen("hitran.h")+1L)*sizeof(char));
  (void)strcat(fl_htrn_h,"hitran.h");
  fl_htrn_pm=(char *)strdup(pth_out);
  fl_htrn_pm=(char *)realloc(fl_htrn_pm,(strlen(pth_out)+strlen("HITRAN.pm")+1L)*sizeof(char));
  (void)strcat(fl_htrn_pm,"HITRAN.pm");
  
  /* Main body of code */
  fp_in=fopen(fl_in,"r"); 
  if(fp_in == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open input file %s\n",prg_nm,fl_in);
      exit(EXIT_FAILURE);
  } /* endif */
  for(idx=0;idx<1;idx++){ /* Skip one line */
    rcd=fscanf(fp_in," %[^\n] ",bfr);
    if(dbg_lvl > 1) (void)fprintf(stderr,"%s\n",bfr);
  } /* end loop over idx */
  
  /* Allocate space for the molecule list */
  mlc=(mlc_sct **)malloc(MLC_NBR_MAX_HTRN*sizeof(mlc_sct *));

  /* Allocate space for the isotope list */
  iso=(iso_sct **)malloc(ISO_NBR_MAX_HTRN*sizeof(iso_sct *));
  
  /* Reset molecule counter */
  idx=-1;
  iso_cnsc_idx=0;
  /* Read until EOF or until allocated memory runs out */
  while(rcd != EOF && idx < MLC_NBR_MAX_HTRN){
    /* Take next line from input file, store it in a string buffer */
    rcd=fscanf(fp_in," %[^\n] ",bfr);
    /* 20001009: As of GCC 2.96, fscanf() began returning 0 rather than EOF here at end of file */
    if(rcd == EOF || rcd == 0) break; /* Exit while loop */
    if(dbg_lvl > dbg_scl) (void)fprintf(stderr,"%s: DEBUG new line contents: %s\n",prg_nm,bfr);
    
    /* If string contains (*) then it starts a new molecule */
    if((strstr(bfr,"(") != NULL) && (strstr(bfr,")") != NULL)) rcd=1; else rcd=0;

    /* Prepare new molecule structure */
    if(rcd){
      /* Ingest molecule information */
      /* fxm: Solaris core dumps here */
      rcd=sscanf(bfr,"%s (%hd)",mlc_sng,&mlc_idx_htrn);
      if(rcd != 2){
	(void)fprintf(stderr,"%s: ERROR Could not find molecule header, exiting...\n",prg_nm);
	exit(EXIT_FAILURE);
      } /* endif */

      /* Reset isotope counter */
      iso_idx=0;

      /* Increment molecule counter */
      idx++;

      if(dbg_lvl > 0) (void)fprintf(stderr,"idx = %d, mlc_sng = %s, mlc_idx_htrn = %d\n",idx,mlc_sng,mlc_idx_htrn);

      /* Allocate space for molecule */
      mlc[idx]=(mlc_sct *)malloc(sizeof(mlc_sct));
      
      /* Fill in molecule structure */
      mlc[idx]->nm=(char *)strdup(mlc_sng);
      mlc[idx]->nm_sf=nm_mk_sf(mlc_sng);
      mlc[idx]->nm_len=strlen(mlc[idx]->nm_sf);
      (void)mlc_grp_set(mlc[idx]);
      (void)rtt_fnc_tpt_xpn_set(mlc[idx]);

      /* Allocate space for isotope list */
      mlc[idx]->iso=&iso[iso_cnsc_idx];
      mlc[idx]->mlc_idx=mlc_idx_htrn;
      mlc[idx]->mmw=0.;

    }else{ /* If string does not contain (*) then it starts a new isotope for current molecule */
      /* Ingest isotope information */
      rcd=sscanf(bfr,"%s %lf %lf %hd %lf",iso_sng,&cnc_frc,&TIPS,&gj,&mmw);
      if(rcd != 5){
	(void)fprintf(stderr,"%s: ERROR Could not parse isotope info\n",prg_nm);
	exit(EXIT_FAILURE);
      } /* endif */
      
      if(dbg_lvl > 1) (void)fprintf(stderr,"iso_idx = %d, iso_sng = %s, cnc_frc = %f, TIPS = %f, gj = %d, mmw = %f\n",iso_idx,iso_sng,cnc_frc,TIPS,gj,mmw);
      
      /* Allocate space for isotope */
      iso[iso_cnsc_idx]=mlc[idx]->iso[iso_idx]=(iso_sct *)malloc(sizeof(iso_sct));
      
      /* Convert input data to SI units */
      mmw*=1.0e-3; /* [mol g-1] --> [mol kg-1] */
      
      /* Fill in isotope structure */
      mlc[idx]->iso[iso_idx]->nm=nm_mk_iso(mlc[idx]->nm,iso_idx+1);
      mlc[idx]->iso[iso_idx]->nm_sf=nm_mk_sf(mlc[idx]->iso[iso_idx]->nm);
      mlc[idx]->iso[iso_idx]->nm_len=strlen(mlc[idx]->iso[iso_idx]->nm_sf);
      mlc[idx]->iso[iso_idx]->iso_abb=(char *)strdup(iso_sng);
      mlc[idx]->iso[iso_idx]->mlc_idx=idx;
      mlc[idx]->iso[iso_idx]->iso_idx=iso_idx+1; /* 0-based --> 1-based */
      mlc[idx]->iso[iso_idx]->iso_cnsc_idx=iso_cnsc_idx+1; /* 0-based --> 1-based */
      mlc[idx]->iso[iso_idx]->mmw=mmw;
      mlc[idx]->iso[iso_idx]->cnc_frc=cnc_frc;
      mlc[idx]->iso[iso_idx]->gj=gj;
      mlc[idx]->iso[iso_idx]->TIPS=TIPS;
      
      /* Increment isotope counters */
      iso_idx++;
      iso_cnsc_idx++;

      /* This is an inelegant way of setting final count */
      mlc[idx]->iso_nbr=iso_idx;
    
    } /* end else isotope line */
    
  } /* end while loop over molecules */

  mlc_nbr=idx+1; /* Number of molecules processed */
  iso_nbr=iso_cnsc_idx; /* Total number of isotopomers */

  /* Compute mean, mean molecular weight of molecule for standard isotopic distribution in Earth's atmosphere */
  for(idx=0;idx<mlc_nbr;idx++){
    mlc[idx]->mmw=0.;
    cnc_frc_ttl=0.;
    for(iso_idx=0;iso_idx<mlc[idx]->iso_nbr;iso_idx++){
      mlc[idx]->mmw+=mlc[idx]->iso[iso_idx]->mmw*mlc[idx]->iso[iso_idx]->cnc_frc;
      cnc_frc_ttl+=mlc[idx]->iso[iso_idx]->cnc_frc;
    } /* end loop over iso */
    mlc[idx]->mmw/=cnc_frc_ttl;
  } /* end loop over mlc */

  /* Search for greatest number of isotopes per molecule */
  max_nbr_iso_per_mlc=0;
  for(idx=0;idx<mlc_nbr;idx++)
    if(mlc[idx]->iso_nbr > max_nbr_iso_per_mlc) max_nbr_iso_per_mlc=mlc[idx]->iso_nbr;
  
  (void)fprintf(stderr,"Ingested data from %s\n",fl_in);
  if(dbg_lvl > 0) (void)fprintf(stderr,"Processed %d molecules, %d isotopomers\n",mlc_nbr,iso_nbr);
  (void)fclose(fp_in); 

  /* Now we have isotopomer data, print it out in desired format */

  /* Create HITRAN definitions for Perl programs */
  fp_htrn_pm=fopen(fl_htrn_pm,"w"); 
  if(fp_htrn_pm == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_pm);
      exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_pm,"# $Id$\n\n");
  (void)fprintf(fp_htrn_pm,"# Purpose: Perl module for use by HITRAN programs.\n\n");
  (void)fprintf(fp_htrn_pm,"# Usage: %s automatically generated by %s on %s\n\n",fl_htrn_pm,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_pm,"# Command: %s\n\n",cmd_ln);

  /* Set up module */
  (void)fprintf(fp_htrn_pm,"package HITRAN;\n");
  (void)fprintf(fp_htrn_pm,"require Exporter;\n");
  (void)fprintf(fp_htrn_pm,"@ISA=qw(Exporter);\n");
  (void)fprintf(fp_htrn_pm,"@EXPORT=qw($mlc_nbr_max_htrn $iso_nbr_max_htrn $iso_per_mlc_nbr_max_htrn %%mlc_sng %%iso_sng %%mlc_iso ");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_pm,"%%%s ",mlc[idx]->nm_sf);
  (void)fprintf(fp_htrn_pm,"); # Export by default\n");

  /* Set dimensions of molecular database arrays */
  (void)fprintf(fp_htrn_pm,"$mlc_nbr_max_htrn=%i; # Number of gases in HITRAN database\n",mlc_nbr);
  (void)fprintf(fp_htrn_pm,"$iso_nbr_max_htrn=%i; # Number of isotopomers in HITRAN database\n",iso_nbr);
  (void)fprintf(fp_htrn_pm,"$iso_per_mlc_nbr_max_htrn=%i; # Maximum number of isotopes of a molecule in HITRAN database\n\n",max_nbr_iso_per_mlc);

  /* Hash of molecule strings */
  (void)fprintf(fp_htrn_pm,"%%mlc_sng=(\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_pm,"\t%d => \'%s\',\n",idx+1,mlc[idx]->nm_sf);
  (void)fprintf(fp_htrn_pm,"\t); # end %%mlc_sng\n\n");

  /* Hash of isotope strings */
  (void)fprintf(fp_htrn_pm,"%%iso_sng=(\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_pm,"\t%d => \'%s\',\n",idx+1,iso[idx]->nm_sf);
  (void)fprintf(fp_htrn_pm,"\t); # end %%iso_sng\n\n");

  /* Hash of lists of isotope numbers */
  (void)fprintf(fp_htrn_pm,"%%mlc_iso=(\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_pm,"\t%s => [",mlc[mlc_idx]->nm_sf);
    for(iso_idx=0;iso_idx<mlc[mlc_idx]->iso_nbr;iso_idx++){
      (void)fprintf(fp_htrn_pm,"%d, ",mlc[mlc_idx]->iso[iso_idx]->iso_cnsc_idx);
    } /* end loop over iso */
    (void)fprintf(fp_htrn_pm,"],\n");
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_pm,"\t); # end %%mlc_iso\n\n");

  /* Set of hashes of isotope strings */
  idx=0;
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_pm,"%%%s=( # HITRAN molecule %d\n",mlc[mlc_idx]->nm_sf,mlc_idx+1);
    for(iso_idx=0;iso_idx<mlc[mlc_idx]->iso_nbr;iso_idx++,idx++) (void)fprintf(fp_htrn_pm,"\t%d => \'%s\',\n",iso_idx+1,iso[idx]->nm_sf);
    (void)fprintf(fp_htrn_pm,"\t); # end %%%s\n\n",mlc[mlc_idx]->nm_sf);
  } /* end loop over molecules */

  (void)fprintf(stderr,"Wrote HITRAN Perl subroutines to %s\n",fl_htrn_pm);
  (void)fclose(fp_htrn_pm); 

  /* Create HITRAN definitions for C programs */
  fp_htrn_h=fopen(fl_htrn_h,"w"); 
  if(fp_htrn_h == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_h);
      exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_h,"/* $Id$ */\n\n");
  (void)fprintf(fp_htrn_h,"/* Purpose: HITRAN definitions used by C programs */\n\n");
  (void)fprintf(fp_htrn_h,"/* Usage: %s automatically generated by %s on %s */\n\n",fl_htrn_h,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_h,"/* Command: %s */\n\n",cmd_ln);
  (void)fprintf(fp_htrn_h,"/* Defining mean molecular weight in SI units instead of g mol-1 keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0 */\n\n");
  for(idx=0;idx<mlc_nbr;idx++){
    (void)fprintf(fp_htrn_h,"const double mmw_%s=%.7e; /* [kg mol-1] mmw of %s */\n",mlc[idx]->nm_sf,mlc[idx]->mmw,mlc[idx]->nm);
  } /* end loop over molecules */
  for(idx=0;idx<iso_nbr;idx++){
    (void)fprintf(fp_htrn_h,"const double mmw_%s=%.7e; /* [kg mol-1] mmw of %s isotope %1d, %s */\n",iso[idx]->nm_sf,iso[idx]->mmw,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  } /* end loop over molecules */

  (void)fprintf(stderr,"Wrote HITRAN C definitions to %s\n",fl_htrn_h);
  (void)fclose(fp_htrn_h); 

  /* Create HITRAN functions for C++ programs */
  fp_htrn_cc=fopen(fl_htrn_cc,"w"); 
  if(fp_htrn_cc == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_cc);
      exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_cc,"// $Id$\n\n");
  (void)fprintf(fp_htrn_cc,"// Purpose: Implementation (declaration) of HITRAN classes\n\n");
  (void)fprintf(fp_htrn_cc,"/* Copyright (C) 1997--2018 Charlie Zender\n   This software is distributed under the terms of the GNU General Public License (GPL) Version 3\n   See http://www.gnu.org/copyleft/gpl.html for full license text */\n\n");
  (void)fprintf(fp_htrn_cc,"/* Usage: %s automatically generated by %s on %s\n",fl_htrn_cc,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_cc,"   Command: %s\n\n",cmd_ln);
  (void)fprintf(fp_htrn_cc,"   Compilation:\n"); 
  (void)fprintf(fp_htrn_cc,"   g++ -Wall -c -O -D$PVM_ARCH -I${HOME}/include ${HOME}/sw/c++/htrn_c++.cc -o ${MY_OBJ_DIR}/htrn_c++.o\n");
  (void)fprintf(fp_htrn_cc,"*/\n");
  (void)fprintf(fp_htrn_cc,"#include <htrn_c++.hh> // HITRAN definitions\n\n");
  (void)fprintf(fp_htrn_cc,"// Namespaces\n");
  (void)fprintf(fp_htrn_cc,"using namespace htrn; // [nms] HITRAN namespace\n\n");
  (void)fprintf(fp_htrn_cc,"// htrn_cls class\n\n");
  (void)fprintf(fp_htrn_cc,"// Friendly functions begin\n\n");

  /* Function to fill in array with mmw of HITRAN molecules */
  (void)sprintf(mmw_mlc_get_fnc_prt,"int // O [enm] Return success code\nmmw_mlc_get // [fnc] Mean molecular weight of HITRAN molecules\n(double *mmw_mlc) // O [kg mol-1] Mean molecular weight of HITRAN molecules\n");
  (void)fprintf(fp_htrn_cc,"%s{\n",mmw_mlc_get_fnc_prt);
  (void)fprintf(fp_htrn_cc,"  // Purpose: Provide mean molecular weight of HITRAN molecules\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  mmw_mlc[0]=1.0e36; // Zeroth element enables array addressing by mlc_id\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_cc,"  mmw_mlc[idx_%s]=%.7e; // [kg mol-1] Mean molecular weight of %s\n",mlc[mlc_idx]->nm_sf,mlc[mlc_idx]->mmw,mlc[mlc_idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"} // end mmw_mlc_get()\n\n");

  /* Function to fill in array with mmw of HITRAN isotopomers */
  (void)sprintf(mmw_iso_get_fnc_prt,"int // O [enm] Return success code\nmmw_iso_get // [fnc] Mean molecular weight of HITRAN isotopomers\n(double *mmw_iso) // O [kg mol-1] Mean molecular weight of HITRAN isotopomers\n");
  (void)fprintf(fp_htrn_cc,"%s{\n",mmw_iso_get_fnc_prt);
  (void)fprintf(fp_htrn_cc,"  // Purpose: Provide mean molecular weight of HITRAN isotopomers\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  mmw_iso[0]=1.0e36; // Zeroth element enables array addressing by iso_id\n");
  for(iso_idx=0;iso_idx<iso_nbr;iso_idx++){
    (void)fprintf(fp_htrn_cc,"  mmw_iso[idx_%s]=%.7e; // [kg mol-1] Mean molecular weight of %s isotope %1d, %s\n",iso[iso_idx]->nm_sf,iso[iso_idx]->mmw,mlc[iso[iso_idx]->mlc_idx]->nm,iso[iso_idx]->iso_idx,iso[iso_idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"} // end mmw_iso_get()\n\n");

  /* Function to fill in array with molecule names */
  (void)sprintf(mlc_sng_get_fnc_prt,"int // O [enm] Return success code\nmlc_sng_get // [fnc] HITRAN molecule names\n(std::string *mlc_sng) // O [sng] HITRAN molecule names\n"); 
  (void)fprintf(fp_htrn_cc,"%s{\n",mlc_sng_get_fnc_prt);
  (void)fprintf(fp_htrn_cc,"  // Purpose: Get array of molecule strings\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  mlc_sng[0]=\"Undefined\"; // Zeroth element enables array addressing by mlc_id\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_cc,"  mlc_sng[idx_%s]=\"%s\";\n",mlc[mlc_idx]->nm_sf,mlc[mlc_idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"} // end mlc_sng_get()\n\n");

  /* Function to fill in array with isotopomer names */
  (void)sprintf(iso_sng_get_fnc_prt,"int // O [enm] Return success code\niso_sng_get // [fnc] HITRAN isotopomer names\n(std::string *iso_sng) // O [sng] HITRAN isotopomer names\n"); 
  (void)fprintf(fp_htrn_cc,"%s{\n",iso_sng_get_fnc_prt);
  (void)fprintf(fp_htrn_cc,"  // Purpose: Get array of isotopomer strings\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  iso_sng[0]=\"Undefined\"; // Zeroth element enables array addressing by iso_id\n");
  for(iso_idx=0;iso_idx<iso_nbr;iso_idx++){
    (void)fprintf(fp_htrn_cc,"  iso_sng[idx_%s]=\"%s\";\n",iso[iso_idx]->nm_sf,iso[iso_idx]->nm);
  } /* end loop over isotopomers */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"} // end iso_sng_get()\n\n");

  /* Function to fill in array with rotational partition function temperature dependent exponent */
  (void)sprintf(rtl_fnc_tpt_xpn_get_fnc_prt,"int // O [enm] Return success code\nrtl_fnc_tpt_xpn_get // [fnc] Exponent defining temperature dependence of rotational partition function\n(double *xpn) // O [frc] Exponent defining temperature dependence of rotational partition function\n");
  (void)fprintf(fp_htrn_cc,"%s{\n",rtl_fnc_tpt_xpn_get_fnc_prt);
  (void)fprintf(fp_htrn_cc,"  /* Purpose: Fills an array containing the exponent which defines\n");
  (void)fprintf(fp_htrn_cc,"     the temperature dependence of the rotational partition function\n");
  (void)fprintf(fp_htrn_cc,"     for each HITRAN molecule.\n");
  (void)fprintf(fp_htrn_cc,"     Exponent is 1 for linear molecules (e.g., CO2, NO, N2O, CO), 3/2 for nonlinear molecules (e.g., H2O, O3, CH4)\n");
  (void)fprintf(fp_htrn_cc,"     See Lio92 p. 33 (2.2.22a) */\n\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  xpn[0]=1.0e36; // Zeroth element enables array addressing by mlc_id\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_cc,"  xpn[idx_%s]=%3.1f; // Group %d, %s\n",mlc[mlc_idx]->nm_sf,mlc[mlc_idx]->rtt_fnc_tpt_xpn,mlc[mlc_idx]->mlc_grp,mlc_grp_sng[mlc[mlc_idx]->mlc_grp]);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"} // end rtl_fnc_tpt_xpn_get()\n\n");

  /* Subroutine to return map between isotopomer index and given molecular and isotopic indices */
  (void)sprintf(iso_idx_map_fnc_prt,"int // O [enm] Return success code\niso_idx_map_get // [fnc] Map [mlc_id,iso_id]->istpmr_id\n(short map[htrn::mlc_nbr_max_htrn+1][htrn::iso_per_mlc_nbr_max_htrn+1]) // O [map] Map [mlc_id,iso_id]->istpmr_id\n");
  (void)fprintf(fp_htrn_cc,"%s{\n",iso_idx_map_fnc_prt);
  (void)fprintf(fp_htrn_cc,"/* Purpose: Fill array which maps molecular and isotopic indices to isotopomer index */\n\n");
  (void)fprintf(fp_htrn_cc,"  int rcd(0); // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"  int iso_idx; // [idx] HITRAN isotope index\n");
  (void)fprintf(fp_htrn_cc,"  int mlc_idx; // [idx] HITRAN molecule index\n\n");
  (void)fprintf(fp_htrn_cc,"  // Initialize the map\n");
  (void)fprintf(fp_htrn_cc,"  for(mlc_idx=0;mlc_idx<=mlc_nbr_max_htrn;mlc_idx++){\n");
  (void)fprintf(fp_htrn_cc,"    for(iso_idx=0;iso_idx<=iso_per_mlc_nbr_max_htrn;iso_idx++){\n");
  (void)fprintf(fp_htrn_cc,"       map[mlc_idx][iso_idx]=0;\n");
  (void)fprintf(fp_htrn_cc,"    } // end loop over iso\n");
  (void)fprintf(fp_htrn_cc,"  } // end loop over mlc\n\n");
  (void)fprintf(fp_htrn_cc,"  // Fill in map for all valid isotopomers\n");
  idx=0;
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_cc,"\n  //  Molecule %d, %s\n",mlc_idx+1,mlc[mlc_idx]->nm);
    for(iso_idx=0;iso_idx<mlc[mlc_idx]->iso_nbr;iso_idx++,idx++){
      (void)fprintf(fp_htrn_cc,"  map[%d][%d]=%d; // [idx] %s\n",mlc_idx+1,iso_idx+1,idx+1,iso[idx]->nm);
    } /* end loop over iso */
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_cc,"  return rcd; // [enm] Return success code\n");
  (void)fprintf(fp_htrn_cc,"\n} // end iso_idx_map_get()\n\n");

  (void)fprintf(stderr,"Wrote HITRAN C++ functions to %s\n",fl_htrn_cc);
  (void)fclose(fp_htrn_cc); 

  /* Create HITRAN definitions for C++ programs */
  fp_htrn_hh=fopen(fl_htrn_hh,"w"); 
  if(fp_htrn_hh == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_hh);
      exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_hh,"// $Id$\n\n");
  (void)fprintf(fp_htrn_hh,"// Purpose: HITRAN definitions used by C++ programs\n\n");
  (void)fprintf(fp_htrn_hh,"/* Copyright (C) 1997--2018 Charlie Zender\n   This software is distributed under the terms of the GNU General Public License (GPL) Version 3\n   See http://www.gnu.org/copyleft/gpl.html for full license text */\n\n");
  (void)fprintf(fp_htrn_hh,"/* Usage: %s automatically generated by %s on %s",fl_htrn_hh,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_hh,"   Command: %s */\n\n",cmd_ln);
  (void)fprintf(fp_htrn_hh,"// #include <htrn_c++.hh> // HITRAN line database definitions\n\n");
  (void)fprintf(fp_htrn_hh,"#ifndef HTRN_CCC_HH // Contents have not yet been inserted in current source file\n");
  (void)fprintf(fp_htrn_hh,"#define HTRN_CCC_HH\n\n");
  (void)fprintf(fp_htrn_hh,"// C++ headers\n");
  (void)fprintf(fp_htrn_hh,"#include <string> // Standard C++ string class\n\n");
  (void)fprintf(fp_htrn_hh,"// Standard C headers\n\n");
  (void)fprintf(fp_htrn_hh,"// Personal headers\n\n");
  (void)fprintf(fp_htrn_hh,"// Forward declarations\n\n");
  (void)fprintf(fp_htrn_hh,"// Namespaces\n\n");
  (void)fprintf(fp_htrn_hh,"// Typedefs\n\n");
  (void)fprintf(fp_htrn_hh,"namespace htrn{ // [nms] HITRAN namespace\n");
  (void)fprintf(fp_htrn_hh,"\n  // HITRAN dimensions\n");
  (void)fprintf(fp_htrn_hh,"  const int mlc_nbr_max_htrn(%d); // [nbr] Number of gases in HITRAN database\n",mlc_nbr);
  (void)fprintf(fp_htrn_hh,"  const int iso_nbr_max_htrn(%d); // [nbr] Number of isotopomers in HITRAN database\n",iso_nbr);
  (void)fprintf(fp_htrn_hh,"  const int iso_per_mlc_nbr_max_htrn(%d); // [nbr] Maximum number of isotopes of a molecule in HITRAN database\n",max_nbr_iso_per_mlc);
  (void)fprintf(fp_htrn_hh,"\n  // Defining mean molecular weight in SI units instead of g mol-1 keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0\n");
  for(idx=0;idx<mlc_nbr;idx++){
    (void)fprintf(fp_htrn_hh,"  const double mmw_%s(%.7e); // [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->mmw,mlc[idx]->nm);
  } /* end loop over molecules */
  for(idx=0;idx<iso_nbr;idx++){
    (void)fprintf(fp_htrn_hh,"  const double mmw_%s(%.7e); // [kg mol-1] mmw of %s isotope %1d, %s\n",iso[idx]->nm_sf,iso[idx]->mmw,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_hh,"\n  // Integer Fortran (1-based) indices for all HITRAN molecules\n");
  for(idx=0;idx<mlc_nbr;idx++){
    (void)fprintf(fp_htrn_hh,"  const int idx_%s(%d); // [enm] HITRAN molecule number for %s\n",mlc[idx]->nm_sf,mlc[idx]->mlc_idx,mlc[idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_hh,"\n  // Integer Fortran (1-based) indices for all HITRAN isotopomers\n");
  for(idx=0;idx<iso_nbr;idx++){
    (void)fprintf(fp_htrn_hh,"  const int idx_%s(%d); // [enm] HITRAN isotopomer number for %s isotope %d, %s\n",iso[idx]->nm_sf,iso[idx]->iso_cnsc_idx,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  } /* end loop over molecules */
  (void)fprintf(fp_htrn_hh,"} // end HITRAN namespace htrn\n\n");

  (void)fprintf(fp_htrn_hh,"// Prototype global functions with C++ linkages\n");
  (void)fprintf(fp_htrn_hh,"%s; // end mmw_mlc_get_fnc() prototype\n\n",mmw_mlc_get_fnc_prt);
  (void)fprintf(fp_htrn_hh,"%s; // end mmw_iso_get_fnc() prototype\n\n",mmw_iso_get_fnc_prt);
  (void)fprintf(fp_htrn_hh,"%s; // end mlc_sng_get_fnc() prototype\n\n",mlc_sng_get_fnc_prt);
  (void)fprintf(fp_htrn_hh,"%s; // end iso_sng_get_fnc() prototype\n\n",iso_sng_get_fnc_prt);
  (void)fprintf(fp_htrn_hh,"%s; // end rtl_fnc_tpt_xpn_get_fnc() prototype\n\n",rtl_fnc_tpt_xpn_get_fnc_prt);
  (void)fprintf(fp_htrn_hh,"%s; // end iso_idx_map_fnc() prototype\n\n",iso_idx_map_fnc_prt);

  (void)fprintf(fp_htrn_hh,"// Define inline'd functions in header so source is visible to calling files\n");
  (void)fprintf(fp_htrn_hh,"\n#endif // HTRN_CCC_HH\n");
  (void)fprintf(stderr,"Wrote HITRAN C++ definitions to %s\n",fl_htrn_hh);
  (void)fclose(fp_htrn_hh); 

  /* Create subroutines for Fortran77 programs */
  fp_htrn_F=fopen(fl_htrn_F,"w"); 
  if(fp_htrn_F == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_F);
    exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_F,"c     $Id$\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: HITRAN subroutines used by swnb, nbm, ck, and lbl\n\n");
  (void)fprintf(fp_htrn_F,"c     Compilation: \n");
  (void)fprintf(fp_htrn_F,"c     pgf90 -c -fast -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign -DLINUX -I. -I${HOME}/include -Di386 -I/usr/local/include -o ${HOME}/obj/LINUX/hitran.o %s\n",fl_htrn_F);
  (void)fprintf(fp_htrn_F,"c     g77 -c -O -ffixed-line-length-132 -DLINUX -I. -I${HOME}/include -Di386 -I/usr/local/include -o ${HOME}/obj/LINUX/hitran.o %s\n",fl_htrn_F);
  (void)fprintf(fp_htrn_F,"c     f90 -c -xs -stackvar -e  -fast -DSUNMP -I. -I${HOME}/include -I/contrib/include  -o ${HOME}/SUNMP/hitran.o %s\n",fl_htrn_F);
  (void)fprintf(fp_htrn_F,"c     f90 -cpp -c -64 -mips4 -extend_source -mp -mpio -O2 -DSGI64 -I. -I${HOME}/include -I/usr/local/include -o ${HOME}/obj/SGIMP64/hitran.o %s\n\n",fl_htrn_F);
  (void)fprintf(fp_htrn_F,"c     Usage: %s automatically generated by %s on %s\n",fl_htrn_F,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_F,"c     Command: %s\n\n",cmd_ln);

  /* Define tokens to take care of autodoubling */
  (void)fprintf(fp_htrn_F,"#ifdef DOUBLE\n");
  (void)fprintf(fp_htrn_F,"#define COMPUTATIONAL_PRECISION double precision\n");
  (void)fprintf(fp_htrn_F,"#else\n");
  (void)fprintf(fp_htrn_F,"#define COMPUTATIONAL_PRECISION real\n");
  (void)fprintf(fp_htrn_F,"#endif\n\n");

  /* Subroutine to fill in array with rotational partition function temperature dependent exponent */
  (void)fprintf(fp_htrn_F,"      subroutine rtl_fnc_tpt_xpn_get(xpn)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills an array containing the exponent which defines\n");
  (void)fprintf(fp_htrn_F,"c     the temperature dependence of the rotational partition function\n");
  (void)fprintf(fp_htrn_F,"c     for each HITRAN molecule\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      real xpn(mlc_nbr_max_htrn)  ! Exponent defining temperature dependence of rotational partition function\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F,"      xpn(idx_%s)=%3.1f             ! Group %d, %s\n",mlc[idx]->nm_sf,mlc[idx]->rtt_fnc_tpt_xpn,mlc[idx]->mlc_grp,mlc_grp_sng[mlc[idx]->mlc_grp]);
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end rtl_fnc_tpt_xpn_get()\n\n");

  /* Subroutine to fill in array with molecular mmw */
  (void)fprintf(fp_htrn_F,"      subroutine mmw_mlc_get(mmw)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills an array containing the mean molecular weight\n");
  (void)fprintf(fp_htrn_F,"c     for each HITRAN molecule\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      real mmw(mlc_nbr_max_htrn)  ! [kg mol-1] mean molecular weight of all HITRAN molecules\n\n");
  (void)fprintf(fp_htrn_F,"c     Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_F,"c     This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F,"      mmw(idx_%s)=%.7e            ! [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->mmw,mlc[idx]->nm);
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end mmw_mlc_get()\n\n");

  /* Subroutine to fill in array with isotopic mmw */
  (void)fprintf(fp_htrn_F,"      subroutine mmw_iso_get(mmw)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills an array containing the mean molecular weight\n");
  (void)fprintf(fp_htrn_F,"c     for each HITRAN isotope\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      real mmw(iso_nbr_max_htrn)  ! [kg mol-1] mean molecular weight of all HITRAN isotopomers\n\n");
  (void)fprintf(fp_htrn_F,"c     Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_F,"c     This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_F,"      mmw(idx_%s)=%.7e            ! [kg mol-1] mmw of %s\n",iso[idx]->nm_sf,iso[idx]->mmw,iso[idx]->nm);
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end mmw_iso_get()\n\n");

  /* Subroutine to return map between isotopomer index and given molecular and isotopic indices */
  (void)fprintf(fp_htrn_F,"      subroutine iso_idx_map_get(map)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills in an array which maps molecular and isotopic indices to isotopomer index\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      integer map(iso_per_mlc_nbr_max_htrn,mlc_nbr_max_htrn) ! map of mlc,iso indices to isotopomer index\n\n");
  (void)fprintf(fp_htrn_F,"c     Local variables\n");
  (void)fprintf(fp_htrn_F,"      integer iso_idx             ! HITRAN isotope index\n");
  (void)fprintf(fp_htrn_F,"      integer mlc_idx             ! HITRAN molecule index\n\n");
  (void)fprintf(fp_htrn_F,"c     Initialize the map\n");
  (void)fprintf(fp_htrn_F,"      do mlc_idx=1,mlc_nbr_max_htrn\n");
  (void)fprintf(fp_htrn_F,"         do iso_idx=1,iso_per_mlc_nbr_max_htrn\n");
  (void)fprintf(fp_htrn_F,"            map(iso_idx,mlc_idx)=0\n");
  (void)fprintf(fp_htrn_F,"         enddo                  ! end loop over iso\n");
  (void)fprintf(fp_htrn_F,"      enddo                     ! end loop over mlc\n\n");
  (void)fprintf(fp_htrn_F,"c     Fill in the map for all valid isotopomers\n");
  idx=0;
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_F,"\nc     Molecule %d, %s\n",mlc_idx+1,mlc[mlc_idx]->nm);
    for(iso_idx=0;iso_idx<mlc[mlc_idx]->iso_nbr;iso_idx++,idx++){
      (void)fprintf(fp_htrn_F,"      map(%d,%d)=%d           ! %s\n",iso_idx+1,mlc_idx+1,idx+1,iso[idx]->nm);
    } /* end loop over iso */
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end iso_idx_map_get()\n\n");

  /* Subroutine to fill in array with molecule strings */
  (void)fprintf(fp_htrn_F,"      subroutine mlc_sng_get(mlc_sng)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills in an array containing all the HITRAN molecule names\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      character(10) mlc_sng(mlc_nbr_max_htrn) ! Contains all HITRAN molecule names\n\n");
  (void)fprintf(fp_htrn_F,"c     Fill in the array for all molecules\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_F,"      mlc_sng(idx_%s)=\'%s\'\n",mlc[mlc_idx]->nm_sf,mlc[mlc_idx]->nm);
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end mlc_sng_get()\n\n");

  /* Subroutine to fill in array with isotopomer strings */
  (void)fprintf(fp_htrn_F,"      subroutine iso_sng_get(iso_sng)\n\n");
  (void)fprintf(fp_htrn_F,"c     Purpose: Fills in an array containing all the HITRAN isotopomer names\n\n");
  (void)fprintf(fp_htrn_F,"      implicit none\n\n");
  (void)fprintf(fp_htrn_F,"#include <hitran.com>\n\n");
  (void)fprintf(fp_htrn_F,"c     Input/Output args\n");
  (void)fprintf(fp_htrn_F,"      character(20) iso_sng(iso_nbr_max_htrn) ! Contains all HITRAN isotopomer names\n\n");
  (void)fprintf(fp_htrn_F,"c     Fill in the array for all isotopomers\n");
  for(iso_idx=0;iso_idx<iso_nbr;iso_idx++){
    (void)fprintf(fp_htrn_F,"      iso_sng(idx_%s)=\'%s\'\n",iso[iso_idx]->nm_sf,iso[iso_idx]->nm);
  } /* end loop over iso */
  (void)fprintf(fp_htrn_F,"\n      return\n");
  (void)fprintf(fp_htrn_F,"      end                       ! end iso_sng_get()\n\n");

  (void)fprintf(stderr,"Wrote HITRAN Fortran subroutines to %s\n",fl_htrn_F);
  (void)fclose(fp_htrn_F); 

  /* Create HITRAN common block for Fortran programs */
  fp_htrn_com=fopen(fl_htrn_com,"w"); 
  if(fp_htrn_com == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_com);
      exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_com,"c     $Id$\n\n");
  (void)fprintf(fp_htrn_com,"c     Purpose: HITRAN common blocks used by swnb, ck, and lbl\n\n");
  (void)fprintf(fp_htrn_com,"c     Usage: %s automatically generated by %s on %s\n",fl_htrn_com,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_com,"c     Command: %s\n\n",cmd_ln);

  /* Set dimensions of molecular database arrays */
  (void)fprintf(fp_htrn_com,"      integer mlc_nbr_max_htrn        ! Number of gases in HITRAN database\n");
  (void)fprintf(fp_htrn_com,"      parameter(mlc_nbr_max_htrn=%i) ! Number of gases in HITRAN database\n\n",mlc_nbr);
  (void)fprintf(fp_htrn_com,"      integer iso_nbr_max_htrn        ! Number of isotopomers in HITRAN database\n");
  (void)fprintf(fp_htrn_com,"      parameter(iso_nbr_max_htrn=%i) ! Number of isotopomers in HITRAN database\n\n",iso_nbr);
  (void)fprintf(fp_htrn_com,"      integer iso_per_mlc_nbr_max_htrn ! Maximum number of isotopes of a molecule in HITRAN database\n");
  (void)fprintf(fp_htrn_com,"      parameter(iso_per_mlc_nbr_max_htrn=%i) ! Maximum number of isotopes of a molecule in HITRAN database\n\n",max_nbr_iso_per_mlc);

  /* Define mmw for each molecule */
  (void)fprintf(fp_htrn_com,"c     Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_com,"c     This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_com,"      real mmw_%s             ! [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->nm);
  (void)fprintf(fp_htrn_com,"      parameter(\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_com,"     $     mmw_%s=%.7e%s       ! [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->mmw,(idx != mlc_nbr-1) ? "," : "",mlc[idx]->nm);
  (void)fprintf(fp_htrn_com,"     $)\n\n");

  /* Define mmw for each isotopomer */
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_com,"      real mmw_%s             ! [kg mol-1] mmw of %s isotope %1d, %s\n",iso[idx]->nm_sf,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_com,"      parameter(\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_com,"     $     mmw_%s=%.7e%s       ! [kg mol-1] mmw of %s isotope %1d, %s\n",iso[idx]->nm_sf,iso[idx]->mmw,(idx != iso_nbr-1) ? "," : "",mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_com,"     $)\n\n");

  /* Write integer Fortran indices for all molecules */
  (void)fprintf(fp_htrn_com,"c     Integer Fortran (1-based) indices for all HITRAN molecules\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_com,"      integer idx_%s             ! HITRAN molecule number for %s\n",mlc[idx]->nm_sf,mlc[idx]->nm);
  (void)fprintf(fp_htrn_com,"      parameter(\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_com,"     $     idx_%s=%d%s       ! HITRAN molecule number for %s\n",mlc[idx]->nm_sf,mlc[idx]->mlc_idx,(idx != mlc_nbr-1) ? "," : "",mlc[idx]->nm);
  (void)fprintf(fp_htrn_com,"     $)\n\n");

  /* Write integer Fortran indices for all isotopomers */
  (void)fprintf(fp_htrn_com,"c     Integer Fortran (1-based) indices for all HITRAN isotopomers\n\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_com,"      integer idx_%s             ! HITRAN isotopomer number for %s isotope %d, %s\n",iso[idx]->nm_sf,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_com,"      parameter(\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_com,"     $     idx_%s=%d%s       ! HITRAN isotopomer number for %s isotope %d, %s\n",iso[idx]->nm_sf,iso[idx]->iso_cnsc_idx,(idx != iso_nbr-1) ? "," : "",mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_com,"     $)\n\n");

  (void)fprintf(stderr,"Wrote HITRAN Fortran common blocks to %s\n",fl_htrn_com);
  (void)fclose(fp_htrn_com); 

  /* Define public parameters for Fortran90 programs */
  fp_htrn_F90=fopen(fl_htrn_F90,"w"); 
  if(fp_htrn_F90 == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_htrn_F90);
    exit(EXIT_FAILURE);
  } /* endif */
  (void)fprintf(fp_htrn_F90,"! $Id$\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: HITRAN constants and subroutines used by radiative transfer programs swnb, nbm, ck, lbl\n\n");
  (void)fprintf(fp_htrn_F90,"! Copyright (C) 1994--2018 Charlie Zender\n! This software is distributed under the terms of the GNU General Public License (GPL) Version 3\n! See http://www.gnu.org/copyleft/gpl.html for full license text\n\n");
  (void)fprintf(fp_htrn_F90,"! Compilation: \n");
  (void)fprintf(fp_htrn_F90,"! pgf90 -c -fast -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign -DLINUX -I. -I${HOME}/include -Di386 -I/usr/local/include -o ${HOME}/obj/LINUX/hitran.o %s\n",fl_htrn_F90);
  (void)fprintf(fp_htrn_F90,"! f90 -c -xs -stackvar -e  -fast -DSUNMP -I. -I${HOME}/include -I/contrib/include  -o ${HOME}/SUNMP/hitran.o %s\n",fl_htrn_F90);
  (void)fprintf(fp_htrn_F90,"! f90 -cpp -c -64 -mips4 -extend_source -mp -mpio -O2 -DSGIMP64 -I. -I${HOME}/include -I/usr/local/include -o ${HOME}/obj/SGIMP64/hitran.o %s\n\n",fl_htrn_F90);
  (void)fprintf(fp_htrn_F90,"! Generation: %s automatically generated by %s on %s\n",fl_htrn_F90,prg_nm,time_bfr_srt);
  (void)fprintf(fp_htrn_F90,"! Command: %s\n\n",cmd_ln);
  (void)fprintf(fp_htrn_F90,"! Usage:\n! use htrn_mdl ! [mdl] HITRAN constants, subroutines\n\n");
  (void)fprintf(fp_htrn_F90,"module htrn_mdl\n\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n\n");

  /* Set dimensions of molecular database arrays */
  (void)fprintf(fp_htrn_F90,"  ! Dimensions of molecular database arrays\n");
  (void)fprintf(fp_htrn_F90,"  integer,parameter::mlc_nbr_max_htrn=%i ! Number of gases in HITRAN database\n",mlc_nbr);
  (void)fprintf(fp_htrn_F90,"  integer,parameter::iso_nbr_max_htrn=%i ! Number of isotopomers in HITRAN database\n",iso_nbr);
  (void)fprintf(fp_htrn_F90,"  integer,parameter::iso_per_mlc_nbr_max_htrn=%i ! Maximum number of isotopes of a molecule in HITRAN database\n",max_nbr_iso_per_mlc);

  /* Define mmw for each molecule */
  (void)fprintf(fp_htrn_F90,"! Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_F90,"! This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F90,"  real,parameter::mmw_%s=%.7e ! [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->mmw,mlc[idx]->nm);
  (void)fprintf(fp_htrn_F90,"\n");

  /* Define mmw for each isotopomer */
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_F90,"  real,parameter::mmw_%s=%.7e ! [kg mol-1] mmw of %s isotope %1d, %s\n",iso[idx]->nm_sf,iso[idx]->mmw,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_F90,"\n");

  /* Write integer Fortran indices for all molecules */
  (void)fprintf(fp_htrn_F90,"! Integer Fortran (1-based) indices for all HITRAN molecules\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F90,"  integer,parameter::idx_%s=%d ! [idx] HITRAN molecule number for %s\n",mlc[idx]->nm_sf,mlc[idx]->mlc_idx,mlc[idx]->nm);
  (void)fprintf(fp_htrn_F90,"\n");

  /* Write integer Fortran indices for all isotopomers */
  (void)fprintf(fp_htrn_F90,"! Integer Fortran (1-based) indices for all HITRAN isotopomers\n\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_F90,"  integer,parameter::idx_%s=%d ! [idx] HITRAN isotopomer number for %s isotope %d, %s\n",iso[idx]->nm_sf,iso[idx]->iso_cnsc_idx,mlc[iso[idx]->mlc_idx]->nm,iso[idx]->iso_idx,iso[idx]->nm);
  (void)fprintf(fp_htrn_F90,"\n");

  /* Define module subroutines for Fortran90 programs */
  /* Subroutine to fill in array with rotational partition function temperature dependent exponent */
  (void)fprintf(fp_htrn_F90,"contains\n\n");
  (void)fprintf(fp_htrn_F90,"subroutine rtl_fnc_tpt_xpn_get(xpn)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills array with exponent which defines\n");
  (void)fprintf(fp_htrn_F90,"! temperature dependence of rotational partition function\n");
  (void)fprintf(fp_htrn_F90,"! of each HITRAN molecule\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:mlc_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output\n");
  (void)fprintf(fp_htrn_F90,"  real,intent(out)::xpn(mlc_nbr_max_htrn) ! [frc] Exponent defining temperature dependence of rotational partition function\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F90,"  xpn(idx_%s)=%3.1f ! [frc] Group %d, %s\n",mlc[idx]->nm_sf,mlc[idx]->rtt_fnc_tpt_xpn,mlc[idx]->mlc_grp,mlc_grp_sng[mlc[idx]->mlc_grp]);
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine rtl_fnc_tpt_xpn_get\n\n");

  /* Subroutine to fill in array with molecular mmw */
  (void)fprintf(fp_htrn_F90,"subroutine mmw_mlc_get(mmw)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills an array containing the mean molecular weight\n");
  (void)fprintf(fp_htrn_F90,"! for each HITRAN molecule\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:mlc_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output\n");
  (void)fprintf(fp_htrn_F90,"  real,intent(out)::mmw(mlc_nbr_max_htrn) ! [kg mol-1] mean molecular weight of all HITRAN molecules\n\n");
  (void)fprintf(fp_htrn_F90,"! Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_F90,"! This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<mlc_nbr;idx++) (void)fprintf(fp_htrn_F90,"  mmw(idx_%s)=%.7e ! [kg mol-1] mmw of %s\n",mlc[idx]->nm_sf,mlc[idx]->mmw,mlc[idx]->nm);
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine mmw_mlc_get\n\n");

  /* Subroutine to fill in array with isotopic mmw */
  (void)fprintf(fp_htrn_F90,"subroutine mmw_iso_get(mmw)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills array with mean molecular weight\n");
  (void)fprintf(fp_htrn_F90,"! of each HITRAN isotope\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:iso_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output\n");
  (void)fprintf(fp_htrn_F90,"  real,intent(out)::mmw(iso_nbr_max_htrn) ! [kg mol-1] mean molecular weight of all HITRAN isotopomers\n\n");
  (void)fprintf(fp_htrn_F90,"! Define mean molecular weight in SI units instead of g mol-1\n");
  (void)fprintf(fp_htrn_F90,"! This keeps Avagadro's number fixed and gets rid of many factors of 1000.0\n\n");
  for(idx=0;idx<iso_nbr;idx++) (void)fprintf(fp_htrn_F90,"  mmw(idx_%s)=%.7e ! [kg mol-1] mmw of %s\n",iso[idx]->nm_sf,iso[idx]->mmw,iso[idx]->nm);
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine mmw_iso_get\n\n");

  /* Subroutine to return map between isotopomer index and given molecular and isotopic indices */
  (void)fprintf(fp_htrn_F90,"subroutine iso_idx_map_get(map)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills array which maps molecular and isotopic indices to isotopomer index\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:iso_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output args\n");
  (void)fprintf(fp_htrn_F90,"  integer,intent(out)::map(iso_per_mlc_nbr_max_htrn,mlc_nbr_max_htrn) ! [idx] Map of mlc,iso indices to isotopomer index\n\n");
  (void)fprintf(fp_htrn_F90,"! Initialize map\n");
  (void)fprintf(fp_htrn_F90,"  map(:,:)=0\n");
  (void)fprintf(fp_htrn_F90,"! Fill in map for all valid isotopomers\n");
  idx=0;
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_F90,"\n! Molecule %d, %s\n",mlc_idx+1,mlc[mlc_idx]->nm);
    for(iso_idx=0;iso_idx<mlc[mlc_idx]->iso_nbr;iso_idx++,idx++){
      (void)fprintf(fp_htrn_F90,"  map(%d,%d)=%d ! [idx] %s\n",iso_idx+1,mlc_idx+1,idx+1,iso[idx]->nm);
    } /* end loop over iso */
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine iso_idx_map_get\n\n");

  /* Subroutine to fill in array with molecule strings */
  (void)fprintf(fp_htrn_F90,"subroutine mlc_sng_get(mlc_sng)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills in an array containing all the HITRAN molecule names\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:mlc_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output\n");
  (void)fprintf(fp_htrn_F90,"  character(10),intent(out)::mlc_sng(mlc_nbr_max_htrn) ! [sng] Contains all HITRAN molecule names\n\n");
  (void)fprintf(fp_htrn_F90,"! Fill in array for all molecules\n");
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    (void)fprintf(fp_htrn_F90,"  mlc_sng(idx_%s)=\"%s\"\n",mlc[mlc_idx]->nm_sf,mlc[mlc_idx]->nm);
  } /* end loop over mlc */
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine mlc_sng_get\n\n");

  /* Subroutine to fill in array with isotopomer strings */
  (void)fprintf(fp_htrn_F90,"subroutine iso_sng_get(iso_sng)\n\n");
  (void)fprintf(fp_htrn_F90,"! Purpose: Fills array containing all HITRAN isotopomer names\n\n");
  (void)fprintf(fp_htrn_F90,"! use htrn_cst,only:iso_nbr_max_htrn ! [mdl] HITRAN constants, subroutines\n");
  (void)fprintf(fp_htrn_F90,"  implicit none\n");
  (void)fprintf(fp_htrn_F90,"! Output\n");
  (void)fprintf(fp_htrn_F90,"  character(20),intent(out)::iso_sng(iso_nbr_max_htrn) ! [sng] Contains all HITRAN isotopomer names\n\n");
  (void)fprintf(fp_htrn_F90,"! Fill array for all isotopomers\n");
  for(iso_idx=0;iso_idx<iso_nbr;iso_idx++){
    (void)fprintf(fp_htrn_F90,"  iso_sng(idx_%s)=\"%s\"\n",iso[iso_idx]->nm_sf,iso[iso_idx]->nm);
  } /* end loop over iso */
  (void)fprintf(fp_htrn_F90,"  return\n");
  (void)fprintf(fp_htrn_F90,"end subroutine iso_sng_get\n\n");
  (void)fprintf(fp_htrn_F90,"end module htrn_mdl\n\n");

  (void)fprintf(stderr,"Wrote HITRAN Fortran90 module to %s\n",fl_htrn_F90);
  (void)fclose(fp_htrn_F90); 
  
  Exit_gracefully();
  return EXIT_SUCCESS;
} /* end main() */

void 
usg_prn(char *opt_sng)
{
  (void)fprintf(stderr,"\nusage: htrn [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",opt_sng);
  (void)fprintf(stderr,"-D dbg_lvl debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-d dbg_val second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-f flt_foo Set generic float. Default is 0.\n");
  (void)fprintf(stderr,"-v print CVS program version\n");
  (void)fprintf(stderr,"\n");
} /* end usg_prn() */

void 
Exit_gracefully(void)
{
  char *time_bfr_end;
  time_t clock;

  /* end clock */ 
  
  clock=time((time_t *)NULL);
  time_bfr_end=ctime(&clock);
  (void)fprintf(stderr,"\tfinish = %s\n",time_bfr_end);

  (void)fclose(stderr);
  (void)fclose(stdin);
  (void)fclose(stdout);

  exit(0);
} /* end Exit_gracefully() */

char *
cmd_ln_sng(int argc,char **argv)
/* 
   int argc: input argument count
   char **argv: input argument list
   char *cmd_ln_sng(): output command line
*/
{
  char *cmd_ln;
  
  int cmd_ln_sz=0;
  int idx;

  for(idx=0;idx<argc;idx++){
    cmd_ln_sz+=(int)strlen(argv[idx])+1;
  } /* end loop over args */
  cmd_ln=(char *)malloc(cmd_ln_sz*sizeof(char));
  if(argc <= 0){
    cmd_ln=(char *)malloc(sizeof(char));
    cmd_ln[0]='\0';
  }else{
    (void)strcpy(cmd_ln,argv[0]);
    for(idx=1;idx<argc;idx++){
      (void)strcat(cmd_ln," ");
      (void)strcat(cmd_ln,argv[idx]);
    } /* end loop over args */
  } /* end else */

  return cmd_ln;
} /* end cmd_ln_sng() */

char *
nm_mk_sf(char *sng_in)
/*  
 char *sng_in: input string to make namespace-safe
 char *nm_mk_sf(): output namespace-safe string
*/
{
  /* Routine creates a copy of sng_in and, if necessary, modifies it to be safe to use in C and Fortran programs */  

  char *sng_out;
  char *sng_ptr;

  if(dbg_lvl == dbg_crr) (void)fprintf(stderr,"%s: DEBUG nm_mk_sf() reports sng_in = %s\n",prg_nm,sng_in);
  if(!sng_in) (void)fprintf(stderr,"%s: ERROR nm_mk_sf() reports sng_in is NULL\n",prg_nm);
  sng_out=(char *)strdup(sng_in);

  /* Currently only unsafe HITRAN name is "NO+", because it contains a '+' */      
  if(strchr(sng_in,'+') != NULL){
    sng_out=(char *)strdup(sng_in);
    sng_out=realloc(sng_out,(strlen(sng_in)+3+1)*sizeof(char));
    sng_ptr=strchr(sng_out,'+');
    (void)strcpy(sng_ptr,"_ion"); /* Replace '+' with "_ion" */
  } /* end if */
  
  return sng_out;
} /* end nm_mk_sf() */

void
rtt_fnc_tpt_xpn_set(mlc_sct *mlc)
/*  
    mlc_sct mlc: pointer to input molecule structure
    void rtt_fnc_tpt_xpn_set():
*/
{
  /* Routine sets exponent defining temperature dependence of molecule's rotational partition function */

  /* Purpose: Return exponent which defines temperature dependence of rotational partition function.
    Exponent is generally taken to be 1.5 for non-linear molecules, and 1.0 for linear molecules.
    See discussion in Lio92 p. 33, RGT92 Table 5 p. 494 (updated version of RGG87 Table III p. 4061),  
    which segregates all molecules into one of six classes:
    1. Asymmetric Rotors
    2. Diatomic and Linear Molecules with Integer J
    3. Spherical Rotors
    4. Symmetric Rotors
    5. Triplet Ground Electronic States
    6. Doublet Ground Electronic States (Half Integer J)
    My system for assigning temperature exponent of rotational partition function is as follows:
    Molecules in HITRAN groups 2,5,6 are assigned linear molecule exponent, i.e., n=1.0
    Molecules in HITRAN groups 1,3,4 are assigned nonlinear molecule exponent, i.e., n=1.5 */

  float tpt_xpn;

  short mlc_grp;

  mlc_grp=mlc->mlc_grp;

  if(mlc_grp < 0 || mlc_grp > 6) (void)fprintf(stderr,"%s: ERROR molecule %s has unknown molecular structure group %d in rtt_fnc_tpt_xpn_set()\n",prg_nm,mlc->nm,mlc_grp);
  if(mlc_grp == unknown_structure){
    (void)fprintf(stderr,"%s: WARNING Setting tpt_xpn=1.5 for molecule %s in rtt_fnc_tpt_xpn_set(). \n",prg_nm,mlc->nm);
    mlc_grp=asymmetric_rotor;
  } /* endif */

  switch(mlc_grp){
  case diatomic_and_linear_molecule_with_integer_J:
  case triplet_ground_electronic_state:
  case doublet_ground_electronic_state_with_half_integer_J:
    tpt_xpn=1.0;
    break;
  case asymmetric_rotor:
  case spherical_rotor:
  case symmetric_rotor:
    tpt_xpn=1.5;
    break;
  default:
    tpt_xpn=1.5;
    (void)fprintf(stderr,"%s: WARNING Setting tpt_xpn = %f for mlc_grp %d in rtt_fnc_tpt_xpn_set()\n",prg_nm,tpt_xpn,mlc_grp);
    break;
  } /* end switch */

  mlc->rtt_fnc_tpt_xpn=tpt_xpn;

} /* end rtt_fnc_tpt_xpn_set() */

void
mlc_grp_set(mlc_sct *mlc)
/*  
 mlc_sct mlc: pointer to input molecule structure
 void mlc_grp_set(): 
*/
{
  /* Routine returns index corresponding to HITRAN molecular structure group 
     Group members are defined in RRG98 p. 670 Tbl. 4 */
  short mlc_grp;
  char *mlc_sng;

  mlc_grp=-1;
  mlc_sng=mlc->nm;

  if(!strcmp(mlc_sng,"H2O")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"CO2")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"O3")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"N2O")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"CO")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"CH4")) mlc_grp=spherical_rotor;
  if(!strcmp(mlc_sng,"O2")) mlc_grp=triplet_ground_electronic_state;
  if(!strcmp(mlc_sng,"NO")) mlc_grp=doublet_ground_electronic_state_with_half_integer_J;
  if(!strcmp(mlc_sng,"SO2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"NO2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"NH3")) mlc_grp=symmetric_rotor;
  if(!strcmp(mlc_sng,"HNO3")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"OH")) mlc_grp=doublet_ground_electronic_state_with_half_integer_J;
  if(!strcmp(mlc_sng,"HF")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"HCl")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"HBr")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"HI")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"ClO")) mlc_grp=doublet_ground_electronic_state_with_half_integer_J;
  if(!strcmp(mlc_sng,"OCS")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"H2CO")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"HOCl")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"N2")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"HCN")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"CH3Cl")) mlc_grp=symmetric_rotor;
  if(!strcmp(mlc_sng,"H2O2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"C2H2")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"C2H6")) mlc_grp=symmetric_rotor;
  if(!strcmp(mlc_sng,"PH3")) mlc_grp=symmetric_rotor;
  if(!strcmp(mlc_sng,"COF2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"SF6")) mlc_grp=symmetric_rotor;
  if(!strcmp(mlc_sng,"H2S")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"HCOOH")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"HO2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"O")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"ClONO2")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"NO+")) mlc_grp=diatomic_and_linear_molecule_with_integer_J;
  if(!strcmp(mlc_sng,"HOBr")) mlc_grp=asymmetric_rotor;
  if(!strcmp(mlc_sng,"C2H4")) mlc_grp=unknown_structure; /* fxm: Here and below are HITRAN 2008, 2012---need their groups! */
  if(!strcmp(mlc_sng,"CH3OH")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"CH3Br")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"CH3CN")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"CF4")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"C4H2")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"HC3N")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"H2")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"CS")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"SO3")) mlc_grp=unknown_structure;
  if(!strcmp(mlc_sng,"C2N2")) mlc_grp=unknown_structure; /* fxm: Here and below are HITRAN 2016---need their groups! */
  if(!strcmp(mlc_sng,"COCl2")) mlc_grp=unknown_structure;

  if(mlc_grp == unknown_structure) (void)fprintf(stderr,"%s: WARNING molecule %s has unknown molecular structure group in mlc_grp_set()\n",prg_nm,mlc_sng);
  if(mlc_grp == -1){
    (void)fprintf(stderr,"%s: ERROR unknown molecule %s in mlc_grp_set()\n",prg_nm,mlc_sng);
    exit(EXIT_FAILURE);
  } /* end if */
  
  mlc->mlc_grp=mlc_grp;

} /* end mlc_grp_set() */

char *
nm_mk_iso(char *mlc_sng,short iso_idx)
/*  
 char *mlc_sng: input molecule string
 short iso_idx: input isotope index
 char *nm_mk_iso(): output isotope string
*/
{
  /* Routine creates a copy of sng_in and, if necessary, modifies it to be safe to use in C and Fortran programs */  

  char *iso_sng=NULL;

  if(!strcmp(mlc_sng,"H2O"))
    switch(iso_idx){
    case 1: iso_sng="1H2_16O"; break;
    case 2: iso_sng="1H2_18O"; break;
    case 3: iso_sng="1H2_17O"; break;
    case 4: iso_sng="1H_2H_16O"; break;
    case 5: iso_sng="1H_2H_18O"; break;
    case 6: iso_sng="1H_2H_17O"; break;
    case 7: iso_sng="2H_2H_16O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CO2"))
    switch(iso_idx){
    case 1: iso_sng="12C_16O2"; break;
    case 2: iso_sng="13C_16O2"; break;
    case 3: iso_sng="16O_12C_18O"; break;
    case 4: iso_sng="16O_12C_17O"; break;
    case 5: iso_sng="16O_13C_18O"; break;
    case 6: iso_sng="16O_13C_17O"; break;
    case 7: iso_sng="12C_18O2"; break;
    case 8: iso_sng="17O_12C_18O"; break;
    case 9: iso_sng="12C_17O2"; break;
    case 10: iso_sng="13C_18O2"; break;
    case 11: iso_sng="18O_13C_17O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"O3"))
    switch(iso_idx){
    case 1: iso_sng="16O3"; break;
    case 2: iso_sng="16O_16O_18O"; break;
    case 3: iso_sng="16O_18O_16O"; break;
    case 4: iso_sng="16O_16O_17O"; break;
    case 5: iso_sng="16O_17O_16O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */

  if(!strcmp(mlc_sng,"N2O"))
    switch(iso_idx){
    case 1: iso_sng="14N2_16O"; break;
    case 2: iso_sng="14N_15N_16O"; break;
    case 3: iso_sng="15N_14N_16O"; break;
    case 4: iso_sng="14N2_18O"; break;
    case 5: iso_sng="14N2_17O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CO"))
    switch(iso_idx){
    case 1: iso_sng="12C_16O"; break;
    case 2: iso_sng="13C_16O"; break;
    case 3: iso_sng="12C_18O"; break;
    case 4: iso_sng="12C_17O"; break;
    case 5: iso_sng="13C_18O"; break;
    case 6: iso_sng="13C_17O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CH4"))
    switch(iso_idx){
    case 1: iso_sng="12C_1H4"; break;
    case 2: iso_sng="13C_1H4"; break;
    case 3: iso_sng="12C_1H3_2H"; break;
    case 4: iso_sng="13C_1H3_2H"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"O2"))
    switch(iso_idx){
    case 1: iso_sng="16O2"; break;
    case 2: iso_sng="16O_18O"; break;
    case 3: iso_sng="16O_17O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"NO"))
    switch(iso_idx){
    case 1: iso_sng="14N_16O"; break;
    case 2: iso_sng="15N_16O"; break;
    case 3: iso_sng="14N_18O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"SO2"))
    switch(iso_idx){
    case 1: iso_sng="32S_16O2"; break;
    case 2: iso_sng="34S_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"NO2"))
    switch(iso_idx){
    case 1: iso_sng="14N_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"NH3"))
    switch(iso_idx){
    case 1: iso_sng="14N_1H3"; break;
    case 2: iso_sng="15N_1H3"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HNO3"))
    switch(iso_idx){
    case 1: iso_sng="1H_14N_16O3"; break;
    case 2: iso_sng="1H_15N_16O3"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"OH"))
    switch(iso_idx){
    case 1: iso_sng="16O_1H"; break;
    case 2: iso_sng="18O_1H"; break;
    case 3: iso_sng="16O_2H"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HF"))
    switch(iso_idx){
    case 1: iso_sng="1H_19F"; break;
    case 2: iso_sng="2H_19F"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HCl"))
    switch(iso_idx){
    case 1: iso_sng="1H_35Cl"; break;
    case 2: iso_sng="1H_37Cl"; break;
    case 3: iso_sng="2H_35Cl"; break;
    case 4: iso_sng="2H_37Cl"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HBr"))
    switch(iso_idx){
    case 1: iso_sng="1H_79Br"; break;
    case 2: iso_sng="1H_81Br"; break;
    case 3: iso_sng="2H_79Br"; break;
    case 4: iso_sng="2H_81Br"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HI"))
    switch(iso_idx){
    case 1: iso_sng="1H_127I"; break;
    case 2: iso_sng="2H_127I"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"ClO"))
    switch(iso_idx){
    case 1: iso_sng="35Cl_16O"; break;
    case 2: iso_sng="37Cl_16O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"OCS"))
    switch(iso_idx){
    case 1: iso_sng="16O_12C_32S"; break;
    case 2: iso_sng="16O_12C_34S"; break;
    case 3: iso_sng="16O_13C_32S"; break;
    case 4: iso_sng="16O_12C_33S"; break;
    case 5: iso_sng="18O_12C_32S"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"H2CO"))
    switch(iso_idx){
    case 1: iso_sng="1H2_12C_16O"; break;
    case 2: iso_sng="1H2_13C_16O"; break;
    case 3: iso_sng="1H2_12C_18O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HOCl"))
    switch(iso_idx){
    case 1: iso_sng="1H_16O_35Cl"; break;
    case 2: iso_sng="1H_16O_37Cl"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"N2"))
    switch(iso_idx){
    case 1: iso_sng="14N2"; break;
    case 2: iso_sng="14N_15N"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HCN"))
    switch(iso_idx){
    case 1: iso_sng="1H_12C_14N"; break;
    case 2: iso_sng="1H_13C_14N"; break;
    case 3: iso_sng="1H_12C_15N"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CH3Cl"))
    switch(iso_idx){
    case 1: iso_sng="12C_1H_35Cl"; break;
    case 2: iso_sng="12C_1H_37Cl"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"H2O2"))
    switch(iso_idx){
    case 1: iso_sng="1H2_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"C2H2"))
    switch(iso_idx){
    case 1: iso_sng="12C2_1H2"; break;
    case 2: iso_sng="1H_12C_13C_1H"; break;
    case 3: iso_sng="1H_12C_12C_2H"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"C2H6"))
    switch(iso_idx){
    case 1: iso_sng="12C2_1H6"; break;
    case 2: iso_sng="12C_13C_1H6"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"PH3"))
    switch(iso_idx){
    case 1: iso_sng="31P_1H3"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"COF2"))
    switch(iso_idx){
    case 1: iso_sng="12C_16O_19F2"; break;
    case 2: iso_sng="13C_16O_19F2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"SF6"))
    switch(iso_idx){
    case 1: iso_sng="32S_19F6"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"H2S"))
    switch(iso_idx){
    case 1: iso_sng="1H2_32S"; break;
    case 2: iso_sng="1H_34S_1H"; break;
    case 3: iso_sng="1H_33S_1H"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HCOOH"))
    switch(iso_idx){
    case 1: iso_sng="1H2_12C_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HO2"))
    switch(iso_idx){
    case 1: iso_sng="1H_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"O"))
    switch(iso_idx){
    case 1: iso_sng="16O"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"ClONO2"))
    switch(iso_idx){
    case 1: iso_sng="35Cl_16O_14N_16O2"; break;
    case 2: iso_sng="37Cl_16O_14N_16O2"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"NO+"))
    switch(iso_idx){
    case 1: iso_sng="14N_16O+"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HOBr"))
    switch(iso_idx){
    case 1: iso_sng="1H_16O_79Br"; break;
    case 2: iso_sng="1H_16O_81Br"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"C2H4"))
    switch(iso_idx){
    case 1: iso_sng="12C2_1H4"; break;
    case 2: iso_sng="12C_13C_1H4"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CH3OH")) /* Methanol */
    switch(iso_idx){
    case 1: iso_sng="12C2_1H3_16O_1H"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CH3Br")) /* Methyl Bromide */
    switch(iso_idx){
    case 1: iso_sng="12C_1H3_79Br"; break;
    case 2: iso_sng="12C_1H3_81Br"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CH3CN")) /* Methyl cyanide (acetonitrile) */
    switch(iso_idx){
    case 1: iso_sng="12C2_1H3_14N"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"CF4")) /* Tetrafluorocarbon (CFC-14) */
    switch(iso_idx){
    case 1: iso_sng="12C_F4"; break;
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"C4H2")) /* Diacetylene (aka butadiyne) */
    switch(iso_idx){
    case 1: iso_sng="12C4_1H2"; break; /* fxm: 20130611 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"HC3N")) /* Foo */
    switch(iso_idx){
    case 1: iso_sng="1H_12C3_14N"; break; /* fxm: 20130611 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"H2")) /* Hydrogen (dihydrogen) */
    switch(iso_idx){
    case 1: iso_sng="1H"; break; /* fxm: 20130611 isopotomer not yet verified */
    case 2: iso_sng="2H"; break; /* fxm: 20130611 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */

  if(!strcmp(mlc_sng,"CS")) /* Carbon monosulfide fxm */
    switch(iso_idx){
    case 1: iso_sng="12C_32S"; break; /* fxm: 20130611 isopotomer not yet verified */
    case 2: iso_sng="13C_32S"; break; /* fxm: 20130611 isopotomer not yet verified */
    case 3: iso_sng="12C_34S"; break; /* fxm: 20130611 isopotomer not yet verified */
    case 4: iso_sng="13C_34S"; break; /* fxm: 20130611 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"SO3")) /* Sulfur trioxide */
    switch(iso_idx){
    case 1: iso_sng="32S_16O3"; break; /* fxm: 20130611 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"C2N2")) /* Cyanogen */
    switch(iso_idx){
    case 1: iso_sng="12C2_14N2"; break; /* fxm: 20181001 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(!strcmp(mlc_sng,"COCl2")) /* Phosgene */
    switch(iso_idx){
    case 1: iso_sng="12C_16O_35Cl2"; break; /* fxm: 20181001 isopotomer not yet verified */
    case 2: iso_sng="12C_16O_35Cl_37Cl"; break; /* fxm: 20181001 isopotomer not yet verified */
    default: (void)fprintf(stdout,"%s: ERROR unknown %s isotope number %d\n",prg_nm,mlc_sng,iso_idx); break;
    } /* end switch */
  
  if(iso_sng == NULL) (void)fprintf(stderr,"%s: ERROR unknown molecule %s in nm_mk_iso()\n",prg_nm,mlc_sng);

  return iso_sng;
} /* end nm_mk_iso() */
