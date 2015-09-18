/* $Id$ */

/* Purpose: Template for C/C99 programs
   All constructs used are ISO C99-compliant
   No C++ features are used */

/* Usage:
   ./c --dbg=1 --dbl=1.0
   scp ~/c/c.c esmf.ess.uci.edu:c

   gcc -std=c99 -DHAVE_C99 -I${MY_INC_DIR} ${MY_OBJ_DIR}/getopt_bsd.o -o c c.c -lgsl -lgslcblas -lm
   gcc -DHAVE_C99 -maix -Wconversion -o c c.c -lm # GCC on AIX
   pathcc -std=c99 -g -mp -DHAVE_C99 -I${MY_INC_DIR} ${MY_OBJ_DIR}/getopt_bsd.o -o c c.c -lgsl -lgslcblas -lm
   pgcc -c9x -g -mp -DHAVE_C99 -DPGI_CC -I${MY_INC_DIR} ${MY_OBJ_DIR}/getopt_bsd.o -o c c.c -lgsl -lgslcblas -lm
   xlc_r -o c c.c -lm
   xlc_r -qlanglvl=extc99 -DHAVE_C99 -I${MY_INC_DIR} -I/usr/local/include -qlanglvl=extended ${MY_OBJ_DIR}/getopt_bsd.o -o c c.c -L/usr/local/lib -lgsl -lgslcblas -lm
*/

#define bool int

/* Standard C headers */
#include <float.h> /* DBL_MAX, FLT_MAX, LDBL_MAX */
#include <limits.h> /* CHAR_MAX, SHRT_MAX, INT_MAX */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <regex.h> /* POSIX regular expressions library */
#include <stdarg.h> /* va_start, va_arg, va_end */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* strtod, strtol, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* POSIX compliance, exit() */
#include "getopt.h" /* GNU getopt() functionality from BSD my_getopt() */

/* 3rd party vendors */
#include <gsl/gsl_sf_erf.h> /* GNU Scientific Library special functions error functions */
#include <gsl/gsl_sf_gamma.h> /* GNU Scientific Library special functions gamma functions */

#ifdef _OPENMP
#include <omp.h> /* OpenMP pragmas */
#endif /* not _OPENMP */

/* Personal headers */
#include <dbg.h> /* Debugging constants */

/* Global variables declared here */
unsigned short int dbg_lvl=0; /* Debugging level */
unsigned short int dbg_val=0; /* Debugging value */

char *prg_nm; /* [sng] Name of program */
char *prg_nm_get(void); /* [sng] Program name */

#define TKN2SNG_PRV(x) #x
#define TKN2SNG(x) TKN2SNG_PRV(x)

char *tkn_sng_get(void); /* [sng] Missing value attribute name */

#ifndef TKN_SNG
  /*# define TKN_SNG missing_value */
# define TKN_SNG _FillValue
#endif /* TKN_SNG */

/*char nco_mss_val_sng[]=TKN2SNG(TKN_SNG);*/ /* [sng] Missing value attribute name */
char nco_mss_val_sng[]="_FillValue"; /* [sng] Missing value attribute name */
char *tkn_sng_get(void){return nco_mss_val_sng;} /* [sng] Missing value attribute name */

int main(int argc,char **argv)
{

  void * /* O [ptr] Pointer to allocated memory */
    ncz_malloc /* [fnc] Wrapper for malloc() */
    (const size_t sz, /* I [nbr] Number of bytes to allocate */
     int arg_nbr, /* [nbr] Number of arguments */
     ...); /* I [llp] Ellipsis defined in stdarg.h */
  void * void_ptr;

  char *cmd_ln_sng(int argc,char **argv); /* [fnc] Parse command line */
  void Exit_gracefully(void); /* [fnc] Exit program */
  void usg_prn(const char * const opt_sht_lst); /* [fnc] Print correct usage of program */
  char *cmd_ln; /* [sng] Parsed command line */
  char *fl_in; /* [sng] Input file */
  char *fl_out; /* [sng] Output file */
  char *time_bfr_srt; /* [sng] Current date and time */

  char CVS_Id[]="$Id$"; /* [sng] CVS identification string */
  char CVS_Revision[]="$Revision$"; /* [sng] CVS revision string */
  char sbr_nm[]="main()"; /* [sng] Subroutine name */

  double dbl_foo=0.0; /* [frc] Intrinsic double temporary variable */
  double val_dbl; /* Holds short cast as double */

  float flt_foo=0.0; /* [frc] Intrinsic float temporary variable */
  gsl_sf_result nsw_dbl; /* [gsl] GSL result structure */

  int rcd=0; /* [rcd] Return code */

  regmatch_t *result;
  regex_t *rx;

  short sht_foo=73; /* [frc] Intrinsic short temporary variable */
  short val_sht; /* Holds double cast as short */

  time_t time_crr_time_t; /* [tm] Current date and time */

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
  
  /* Set defaults */

  /* Start clock and save command line */
  cmd_ln=cmd_ln_sng(argc,argv); /* [sng] Parsed command line */
  time_crr_time_t=time((time_t *)NULL); /* [tm] Current date and time */
  time_bfr_srt=ctime(&time_crr_time_t); /* [tm] Current date and time */
  (void)fprintf(stderr,"\tStart = %s",time_bfr_srt);
  prg_nm=((prg_nm=strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++; /* [sng] Name of program */

  (void)fprintf(stdout,"%s: tkn_sng_get() returns token TKN_SNG converted with TKN2SNG(x) as char * = \"%s\"\n",prg_nm_get(),tkn_sng_get());

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
      if(dbg_lvl >= dbg_io){
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
      if(optarg) dbg_lvl=(unsigned short int)strtoul(optarg,(char **)NULL,10); else dbg_lvl=dbg_fl;
      break;
    case 'f': /* Generic tuning parameter.  Default is 0. */
      flt_foo=(float)strtod(optarg,(char **)NULL);
      break;
    case 'o': /* Get the output file name. Default is stdout */
      fl_out=opt_sng;
      break;
    case 'v': /* CVS program info */
      (void)fprintf(stderr,"%s %s\n",CVS_Revision,CVS_Id);
      exit(EXIT_SUCCESS);
      break;
    case '-': /* Long options are not allowed */
      (void)fprintf(stderr,"Long options are not available in this build. Use single letter options instead.\n");
      exit(EXIT_FAILURE);
      break;
    default: /* Print proper usage */
      (void)usg_prn(opt_sht_lst);
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

  /* Main body of code */

  FILE * const fp_stdout=stdout; /* [fl] stdout filehandle CEWI */
  /* Using naked stdin/stdout/stderr in parallel region generates xlc/xlC warning
     Copy appropriate filehandle to variable scoped shared in parallel clause */
#ifdef _OPENMP /* OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification */
  int thr_nbr_max; /* [nbr] Maximum number of threads system/user allow program to use */  
#endif /* not _OPENMP */
#ifdef _OPENMP /* OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification */
  thr_nbr_max=omp_get_max_threads(); /* [nbr] Maximum number of threads system/user allow program to use */
#pragma omp parallel default(none) shared(dbg_lvl)
  { /* begin OpenMP parallel */
#pragma omp single nowait
    { /* begin OpenMP single */
      if(dbg_lvl >= 1) (void)fprintf(fp_stdout,"%s: INFO OpenMP threading with %d threads\n",prg_nm_get(),omp_get_num_threads());
    } /* end OpenMP single */
  } /* end OpenMP parallel */
#else /* not _OPENMP */
  (void)fprintf(fp_stdout,"%s: INFO Not attempting OpenMP threading\n",prg_nm_get());
#endif /* not _OPENMP */

  /* Test GSL math routines */
  rcd+=gsl_sf_erf_e(1.0,&nsw_dbl);
  if(fabs(0.8427-nsw_dbl.val)/0.8427 > 0.001){
    (void)fprintf(stdout,"%s: ERROR Error function error\n",prg_nm);
    exit(EXIT_FAILURE);
  } /* endif err */
  rcd+=gsl_sf_erf_e(dbl_foo,&nsw_dbl);
  (void)fprintf(stdout,"gsl_sf_erf(%f) = %f\n",dbl_foo,nsw_dbl.val);

  /* Test intrinsic math routines */
#ifndef AIX
  nsw_dbl.val=tgamma(0.5);
  (void)fprintf(stdout,"tgamma(0.5) = %f\n",nsw_dbl.val);
  if(fabs(1.77245385-nsw_dbl.val)/1.77245385 > 0.001){
    (void)fprintf(stdout,"%s: ERROR Gamma function tgamma() error\n",prg_nm);
    exit(EXIT_FAILURE);
  } /* endif err */
  nsw_dbl.val=tgamma(dbl_foo);
  (void)fprintf(stdout,"tgamma(%f) = %f\n",dbl_foo,nsw_dbl.val);
#endif /* !AIX */

  /* Test variadic malloc() routine */
  size_t mmr_sz=73;
  if((void_ptr=(void *)ncz_malloc(mmr_sz,1)) == NULL){
    (void)fprintf(stdout,"%s: ERROR Unable to malloc() %zu bytes in ncz_var_get()\n",prg_nm,mmr_sz);
    exit(EXIT_FAILURE);
  } /* end if */

  /* Prevent "variable set but not used" warnings */
  cmd_ln=cmd_ln;
  fl_in=fl_in;
  fl_out=fl_out;

  val_dbl=(double)sht_foo; /* Cast short to double */
  val_sht=(short)val_dbl; /* Cast double back to short */
  (void)fprintf(stdout,"%s: sht_foo = %d, val_dbl=%f, val_sht=%d\n",prg_nm,sht_foo,val_dbl,val_sht);
  val_sht=(short)flt_foo; /* Cast arbitrary float to short */
  (void)fprintf(stdout,"%s: flt_foo = %f, val_sht=%d\n",prg_nm,flt_foo,val_sht);

  /* Using C99 features (// comments, mixed declarations, etc.) is OK only in C99 sections
     Many tests from http://www.comeaucomputing.com/features.html */
  (void)fprintf(stdout,"Attempting to test C99 syntax...\n");

  /* Designated initializers */
#ifdef HAVE_C99
  (void)fprintf(stdout,"Testing designated initializers...\n");
  struct sct_typ{int lmn_int;float lmn_flt;int lmn_int_arr[3];};
  (void)fprintf(stdout,"Initialize all structure members in forward order...\n");
#if defined(_AIX) || defined(__INTEL_COMPILER) || defined(__PATHCC__)
  // Normal (non-GCC) C99-compliant compilers
  struct sct_typ sct_1={.lmn_int=3,.lmn_flt=3.123,.lmn_int_arr[0]=-3};
#else // !_AIX, etc.
# if defined(PGI_CC) // Compiler is pgcc
  // Broken, C99-non-compliant compilers
  struct sct_typ sct_1;
  sct_1.lmn_int=3;
  sct_1.lmn_flt=3.123;
  sct_1.lmn_int_arr[0]=-3;
# endif // !PGI_CC
# if defined(__GNUC__) // Compiler is gcc
  // GCC allows this non-ISO construct, i.e., specifying range of elements to initialize
  // c.c:222: warning: ISO C forbids specifying range of elements to initialize
  struct sct_typ sct_1={.lmn_int=3,.lmn_flt=3.123,.lmn_int_arr[0 ... 2]=-3};
# endif // !__GNUC__
#endif // !_AIX, etc.
  (void)fprintf(stdout,"sct_1.lmn_int=%d, sct_1.lmn_flt=%f, sct_1.lmn_int_arr[0]=%d\n",sct_1.lmn_int,sct_1.lmn_flt,sct_1.lmn_int_arr[0]);
  (void)fprintf(stdout,"Initialize some structure members in reverse order...\n");
  struct sct_typ sct_2={.lmn_flt=3.123,.lmn_int=3};
  (void)fprintf(stdout,"Initialize individual array elements within structure...\n");
  struct sct_typ sct_3={.lmn_int_arr[2]=2};
  (void)fprintf(stdout,"Initialize individual array elements of plain array...\n");
  int int_arr_dsg_ntl[3]={[0]=1,[2]=3};
  struct sct_nst_typ{int lmn_int;struct sct_typ sct_typ_lmn;};
  (void)fprintf(stdout,"Nested initialization of elements of structure and sub-structure...\n");
  struct sct_nst_typ sct_nst={.lmn_int=11,.sct_typ_lmn.lmn_int_arr={0,1,2}};
  sct_1=sct_1; // CEWI
  sct_2=sct_2; // CEWI
  sct_3=sct_3; // CEWI
  //int_arr_dsg_ntl=int_arr_dsg_ntl; // CEWI
  sct_nst=sct_nst; // CEWI
#else /* !HAVE_C99 */
  (void)fprintf(stdout,"INFO: Assuming non-C99 compiler, not testing designated initializers...\n");
#endif /* !HAVE_C99 */
  
  /* Compound literals */
#ifdef HAVE_C99
  (void)fprintf(stdout,"Testing compound literals...\n");
  (void)fprintf(stdout,"Define un-named int array...\n");
  int *int_arr=(int[]){1,2}; // Compound literal
  (void)fprintf(stdout,"Define const un-named int array...\n");
  const int *cst_int_arr=(const int[]){1,2}; // Compound literal
  (void)fprintf(stdout,"Define pointer to int value 1...\n");
  int *int_ptr=&(int){1}; // Compound literal
  (void)fprintf(stdout,"Define const pointer to int value 1...\n");
  const int *cst_int_ptr=&(const int){1}; // Compound literal
  (void)fprintf(stdout,"Define const string pointer to value \"Compound literal\"...\n");
  const char *cst_sng=(const char[]){"Compound literal"}; // Compound literal
  cst_sng=cst_sng; // CEWI
  int_arr=int_arr; // CEWI
  cst_int_arr=cst_int_arr; // CEWI
  int_ptr=int_ptr; // CEWI
  cst_int_ptr=cst_int_ptr; // CEWI
#else /* !HAVE_C99 */
  (void)fprintf(stdout,"INFO: Non-GCC compiler, not testing compound literals...\n");
#endif /* !HAVE_C99 */

  /* Restrict type qualifier */
#ifdef HAVE_C99 
  /* Restrict keyword. See, e.g., IBM C/C++ Guide p. 69.
     restrict is a type qualifier.
     "restrict" promises that if memory pointed to by "restrict"-qualified pointer is modified,
     no other pointer will access that same memory.
     This is a promise the programmer must keep.
     It is an assumption the compiler makes, but cannot enforce
     Memory that is not modified may be accessed (read) by multiple restricted pointers

     IBM xlC consumes and ignores "restrict" for compatibility with C99 
     Hence, "restict" has no effect on code compiled by IBM xlC

     "restrict" generates errors when compiled by GNU g++
     Instead, GCC g++ implements (a superset of) "restrict" and requires that the
     type qualifier by `__restrict__' or `__restrict'
     GCC requires restrict in function definitions, but not in function prototypes

     Torvalds seems to think GCC treatment of restrict is broken until GCC 3.3 */
  (void)fprintf(stdout,"Testing restrict keyword...\n");
  int * restrict int_ptr_rst_1;
  int * restrict int_ptr_rst_2;
#ifndef AIX
  /* AIX xlc does not allow assignments between restricted pointers */
  int_ptr_rst_1=int_ptr_rst_1; // CEWI
  int_ptr_rst_2=int_ptr_rst_2; // CEWI
#endif /* !AIX */
#else /* !HAVE_C99 */
  (void)fprintf(stdout,"INFO: Non-GCC compiler, not testing restrict type qualifier...\n");
#endif /* !HAVE_C99 */

  Exit_gracefully();
  return EXIT_SUCCESS;
} /* end main() */

void usg_prn(const char * const opt_sng)
{
  (void)fprintf(stderr,"\nusage: c [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",opt_sng);
  (void)fprintf(stderr,"-D dbg_lvl The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-d dbg_val The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-f flt_foo Set the generic float. Default is 0.0\n");
  (void)fprintf(stderr,"-v print the CVS program version\n");
  (void)fprintf(stderr,"\n");
} /* end usg_prn() */

void Exit_gracefully(void)
{
  char *time_bfr_end;
  time_t time_crr_time_t;

  /* end the clock */
  
  time_crr_time_t=time((time_t *)NULL);
  time_bfr_end=ctime(&time_crr_time_t);
  (void)fprintf(stderr,"\tFinish = %s",time_bfr_end);

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

void * /* O [ptr] Pointer to allocated memory */
ncz_malloc /* [fnc] Wrapper for malloc() */
(const size_t sz, /* I [nbr] Number of bytes to allocate */
 int arg_nbr, /* [nbr] Number of arguments */
 ...) /* I [llp] Ellipsis defined in stdarg.h */
{
  /* Purpose: Custom wrapper for malloc()
     Routine prints error when malloc() returns a NULL pointer 
     Routine does not call malloc() when sz == 0 */
  va_list arg_lst; /* [] Variable argument list */
  /* int arg_idx; *//* [idx] Argument index */

  void *ptr; /* [ptr] Pointer to new buffer */
  
  /* malloc(0) is ANSI-legal, albeit unnecessary
     NCO sometimes employs this degenerate case behavior of malloc() to simplify code
     Some debugging tools like Electric Fence consider any NULL returned by malloc() to be an error
     So circumvent malloc() calls when sz == 0 */
  if(sz == 0) return NULL;
  
  ptr=malloc(sz); /* [ptr] Pointer to new buffer */
  if(ptr == NULL){

    /* Begin variable argument list access */
    va_start(arg_lst,arg_nbr);
  
    /* End variable argument list access */
    va_end(arg_lst);

    (void)fprintf(stdout,"%s: ERROR ncz_malloc() unable to allocate %li bytes\n",prg_nm_get(),(long)sz);
    /* fxm: Should be exit(8) on ENOMEM errors? */
    exit(EXIT_FAILURE);
  } /* endif */
  return ptr; /* [ptr] Pointer to new buffer */
} /* ncz_malloc() */

char *prg_nm_get(void){return prg_nm;} /* [sng] Program name */
