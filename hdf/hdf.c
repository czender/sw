/* $Id$ */

/* Purpose: Template for HDF programs
   All constructs used are ANSI C compliant
   No C++ features are used */

/* Usage:
   hdf --dbl=1.0
   cd ~/hdf;hdf;h5dump $HOME/hdf/foo.h5;cd -
   h5dump $HOME/hdf/foo.h5
   h5ls $HOME/hdf/foo.h5
   h5toh4 $HOME/hdf/foo.h5
   h5debug $HOME/hdf/foo.h5
   hdfls $HOME/hdf/foo.hdf
   h4toh5 $HOME/hdf/foo.hdf
 */

#define bool int

/* Standard C headers */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* all sorts of POSIX stuff */

#if (defined LINUX) || (defined LINUXALPHA)
#include <getopt.h> /* GNU getopt() is standard on Linux */
#else /* not LINUX */
#include "getopt.h" /* GNU getopt() */
#endif /* not LINUX */

/* 3rd party vendors */
#include <hdf5.h> /* HDF5 C Interface */

/* Personal headers */
#include <dbg.h> /* Debugging constants */

/* Global variables declared here */
unsigned short int dbg_lvl=0; /* Debugging level */
unsigned short int dbg_val=0; /* Debugging value */

int main(int argc,char **argv)
{
  char *cmd_ln_sng(int,char **); /* [fnc] Parse command line */
  void Exit_gracefully(void); /* [fnc] Exit program */
  void usg_prn(char *); /* [fnc] Print correct usage of program */
  char *cmd_ln; /* [sng] Parsed command line */
  char *fl_in=NULL; /* [sng] Input file */
  char *fl_out=NULL; /* [sng] Output file */
  char *prg_nm; /* [sng] Name of program */
  char *time_bfr_srt; /* [sng] Current date and time */
  char *pth_in=NULL; /* [drc] Input path */
  char *pth_out=NULL; /* [drc] Output path */
  char *home_nvr=NULL; /* [sng] HOME environment variable */

  char CVS_Id[]="$Id$"; /* [sng] CVS identification string */
  char CVS_Revision[]="$Revision$"; /* [sng] CVS revision string */
  char sbr_nm[]="main()"; /* [sng] Subroutine name */

  /* HDF variables */
  char var_nm[]="IntArray"; /* [sng] Variable name */
  herr_t h5_rcd; /* [enm] HDF return code */
  hid_t dat_set; /* [hnd] HDF dataset handle */
  hid_t dat_spc; /* [hnd] HDF dataspace handle */
  hid_t hdf_out; /* [hnd] HDF file handle */
  hid_t var_typ; /* [hnd] HDF datatype handle */
  long dmn_nbr=2; /* [nbr] Number of dimensions */
  long lat_nbr=5; /* [nbr] Number of latitudes */
  long lon_nbr=6; /* [nbr] Number of longitudes */
  register long lat_idx; /* [idx] Counting index for lat */
  register long lon_idx; /* [idx] Counting index for lon */

  /* Derived variables */
  hsize_t dmn_vct[2]; /* Dataset dimensions */
  int var_val[lon_nbr][lat_nbr]; /* Data to write */

  double dbl_foo=0.0; /* Option dbl_foo */
  float flt_foo=0.0; /* Option f or flt_foo */
  int rcd=0; /* [rcd] Return code */
  short sht_foo=73; /* Option sht_foo */
  double val_dbl; /* Holds short cast as double */
  short val_sht; /* Holds double cast as short */

  time_t clock; /* [tm] Current date and time */

  /* Variables for option processing */
  char *opt_crr; /* String representation of current long-option name */
  char *opt_short_lst; /* String representation of current optarg, if any */
  char *opt_sng=NULL; /* String representation of current optarg, if any */ /* CEWI */
  extern char *optarg; /* char * representation of current optarg, if any (this memory is owned by system) */
  extern int optind; /* extern enumerating cardinal of current option */
  int opt; /* Value is zero if current argument is long type, else value contains single letter version of command line argument */
  int opt_idx=0; /* Index of current long option into opt_lng array */
  static struct option opt_lng[]=
  {
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is compared to enum _argtype{no_argument,required_argument,optional_argument}, 
       flag points to a variable that gets set to val whenever the name option is set.
       For long options that have a zero flag field, getopt() returns the contents of val.
    */
    /* Long options with no argument */
    /* Long options with argument */
    {"dbl_foo",1,0,0}, /* [frc] Intrinsic double temporary variable */
    {"sht_foo",1,0,0}, /* [frc] Intrinsic short temporary variable */
    /* Long options with optional argument */
    /* Long options with short counterparts */
    {"debug",1,0,'D'}, /* [enm] Debugging level */
    {"flt_foo",1,0,'f'}, /* [frc] Intrinsic float temporary variable */
    {"output",1,0,'o'},
    {"version",0,0,'v'},
    {"verbose",0,0,'D'},
    /* Last option named "0" to signal getopt_long() to stop processing */
    {0,0,0,0}
  }; /* end opt_lng */
  
  /* Set defaults */

  /* Default file names */
  home_nvr=getenv("HOME"); /* [sng] HOME environment variable */
  pth_in=(char *)strdup(home_nvr); /* [drc] Input path */
  if(pth_in != NULL){
    /* 5 is size of /ck/, then add one for NUL-terminator */
    pth_in=realloc(pth_in,(5+strlen(pth_in)+1)); /* [drc] Input path */
    pth_in=strcat(pth_in,"/hdf/"); /* [drc] Input path */
  } /* end if */
  fl_in=(char *)strdup(pth_in); /* [sng] Input file name */
  fl_in=(char *)realloc(fl_in,(strlen(pth_in)+strlen("foo.txt")+1)*sizeof(char)); /* [sng] Input file name */
  (void)strcat(fl_in,"foo.txt"); /* [sng] Input file name */
  fl_out=(char *)strdup(pth_in);
  fl_out=(char *)realloc(fl_out,(strlen(pth_in)+strlen("foo.h5")+1)*sizeof(char));
  (void)strcat(fl_out,"foo.h5");

  /* Start clock and save command line */
  cmd_ln=cmd_ln_sng(argc,argv); /* [sng] Parsed command line */
  clock=time((time_t *)NULL); /* [tm] Current date and time */
  time_bfr_srt=ctime(&clock); /* [tm] Current date and time */
  (void)fprintf(stderr,"\tStart = %s",time_bfr_srt);
  prg_nm=((prg_nm=strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++; /* [sng] Name of program */

  /* Short options: no colon = no arg, one colon = required arg, two colons = optional arg */
  opt_short_lst="D:f:o:v"; /* List of single-letter (C-style) option abbreviations */
  /* Parse command line arguments */
  while(1){
    /* getopt_long_only() allows a single dash '-' to prefix long options as well */
    opt=getopt_long_only(argc,argv,opt_short_lst,opt_lng,&opt_idx);
    /* NB: access to opt_crr is only valid when long_opt was detected */
    opt_crr=(char *)strdup(opt_lng[opt_idx].name);
    if(optarg) opt_sng=optarg; /*  Copy system memory into local string for safer operations */
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
      if(optarg) dbg_lvl=(unsigned short int)strtol(optarg,(char **)NULL,10); else dbg_lvl=dbg_fl;
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

  /* Prevent "variable set but not used" warnings */
  cmd_ln=cmd_ln;
  fl_in=fl_in;
  fl_out=fl_out;

  /* Main body of code */

  /* Initialize data */
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++)
    for(lat_idx=0;lat_idx<lat_nbr;lat_idx++)
      var_val[lon_idx][lat_idx]=lat_idx+lon_idx;
  /*
   * 0 1 2 3 4 5 
   * 1 2 3 4 5 6
   * 2 3 4 5 6 7
   * 3 4 5 6 7 8
   * 4 5 6 7 8 9
   */
  
  /* Create new file using H5F_ACC_TRUNC access,
   * default file creation properties, and default file
   * access properties. */
  hdf_out=H5Fcreate(fl_out,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  
  /* Describe size of array and create data space for fixed size dataset */
  dmn_vct[0]=lon_nbr;
  dmn_vct[1]=lat_nbr;
  dat_spc=H5Screate_simple(dmn_nbr,dmn_vct,NULL);
  
  /* Define datatype for data in file
     Store as little endian INT numbers */
  var_typ=H5Tcopy(H5T_NATIVE_INT);
  h5_rcd=H5Tset_order(var_typ,H5T_ORDER_LE);
  
  /* Create new dataset within file using defined dat_spc and var_typ and default dataset creation properties */
  dat_set=H5Dcreate(hdf_out,var_nm,var_typ,dat_spc,H5P_DEFAULT);
  
  /* Write data to dataset using default transfer properties */
  h5_rcd=H5Dwrite(dat_set,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var_val);
  
  /* Close/release resources */
  H5Sclose(dat_spc); /* Close dataspace */
  H5Tclose(var_typ); /* Close datatype */
  H5Dclose(dat_set); /* Close dataset */
  H5Fclose(hdf_out); /* Close output file */
  
  Exit_gracefully();
  return EXIT_SUCCESS;
} /* end main() */

void usg_prn(char *opt_sng)
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
  char *time_bfr_finish;
  time_t clock;

  /* end the clock */
  
  clock=time((time_t *)NULL);
  time_bfr_finish=ctime(&clock);
  (void)fprintf(stderr,"\tFinish = %s",time_bfr_finish);

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





