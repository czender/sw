static char rcs_Id[] = "$Id$\n";
static char rcs_Revision[] = "$Revision$";

/* $Author: zender $
 * $Date$
 * $Locker:  $
 * $RCSfile: io_diff.c,v $
 * $Source: /home/zender/cvs/cld/io_diff.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose: Prints out the difference between two netCDF cloud files. */ 

/* Example usage:
This will produce plots of certain fields in foo.nc minus the corresponding
field in foo.X.nc:
io_diff -D 12 -i foo.nc -j foo.X.nc
*/

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:40  zender
/* Imported sources
/*
 * Revision 4.12  1993/08/24  16:47:08  zender
 * about to branch to a cray compatible version which works
 * with default ncar graphics version 3.2 installation. this
 * mainly means changing to the GKS c-structures and changes
 * in the library location in the makefiles.
 *
 * Revision 4.10  1993/06/11  04:40:29  zender
 * changed CRAY2 to CRAY, i think io_diff is way different, not sure!
 *
 * Revision 4.5  1993/05/25  02:11:55  zender
 * Adapted to read in full descriptions of two files and
 * print out certain plots of both and their difference.
 *
 *  */

#ifdef CRAY
#define ROUNDER Cray_rounder
#else
/* The Cray doesn't know about this IEEE rounding function */ 
#define ROUNDER irint 
#endif

#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include "/usr/local/include/netcdf.h" /* netCDF routines */ 

/* my header files */
#include "defs.h"               /* YES, NO, FILESIZE, and more */
#include "globals.h"            /* externs announced in main() */

#define NUM_IWC_DIFF_CONTOURS 42

static int iasf[13]={
  1,1,1,1,1,1,1,1,1,1,1,1,1};

static char *IWC_llbs[]={
  "0.",
  ".1",".5","1.","2.","5.",
  "10.","15.","20.","25.","30.",
  "40.","50.","60.","70.","80.",
  "90.","100.","125.","150.","175.",
  "200."};

#define NUM_IWC_CONTOURS 21
static float IWC_contour_level[NUM_IWC_CONTOURS+1]={
  0.,
  .1,.5,1.,2.,5.,
  10.,15.,20.,25.,30.,
  40.,50.,60.,70.,80.,
  90.,100.,125.,150.,175.,
  200.};

static char *dist_llbs[]={
  "0.",
  ".03125",".0625",".125",
  ".25",".5","1","2","4",
  "8","16","32","64","128",
  "256","512","1k","2k","4k",
  "8k","16k","32k","64k","128k",
  "256k"};

#define NUM_DIST_CONTOURS 24
static float dist_contour_level[NUM_DIST_CONTOURS+1]={
  0.,
  .03125,.0625,.125,
  .25,.5,1.,2.,4.,
  8.,16.,32.,64.,128.,
  256.,512.,1024.,2048.,4096.,
  8192.,16384.,32768.,65536.,131072.,
  262144.};

static char *dist_diff_llbs[]={
/*  "16k","8k",*/
/*  "4k","2k","1k","512","256",*/
/*  "128","64","32","16","8",*/
/*  "4","2","1",".5",".25",*/
/*  ".125",".0625",".03125",*/
/*  "0",*/
/*  ".03125",".0625",".125",*/
/*  ".25",".5","1","2","4",*/
/*  "8","16","32","64","128",*/
/*  "256","512","1k","2k","4k",*/
/*  "8k","16k"*/
/*  };*/
  "16k","8k",
  "4k","2k","1k","512","256",
  "128","64","32","16","8",
  "4","2","1",".5",".25",
  "0",
  ".25",".5","1","2","4",
  "8","16","32","64","128",
  "256","512","1k","2k","4k",
  "8k","16k"
  };

#define NUM_DIST_DIFF_CONTOURS 34
static float dist_diff_contour_level[NUM_DIST_DIFF_CONTOURS+1]={
/*  0.,*/
/*  -16384.,-8192.,*/
/*  -4096.,-2148.,-1024.,-512.,-256.,*/
/*  -128.,-64.,-32.,-16.,-8.,*/
/*  -4.,-2.,-1.,-.5,-.25,*/
/*  -.125,-.0625,-.03125,*/
/*  0.,*/
/*  .03125,.0625,.125,*/
/*  .25,.5,1.,2.,4.,*/
/*  8.,16.,32.,64.,128.,*/
/*  256.,512.,1024.,2048.,4096.,*/
/*  8192.*/
/*  };*/
  0.,
  -16384.,-8192.,
  -4096.,-2148.,-1024.,-512.,-256.,
  -128.,-64.,-32.,-16.,-8.,
  -4.,-2.,-1.,-.5,-.25,
  0.,
  .25,.5,1.,2.,4.,
  8.,16.,32.,64.,128.,
  256.,512.,1024.,2048.,4096.,
  8192.
  };

float *abscissa;
float *ordinate;

int debug=0; /* Option D */
int debug_value = 0; /* Option d */

FILE *fp_out,*fp_err,*fp_in;

main(int argc,char **argv)
{
  /* Note: the NCAR graphics routines are too numerous to bother declaring,
   they all return type int, i hope */ 
  float equiv_rad_to_bullet_length(float /* length */);
  float max_VEC();
  float min_VEC();
  float total_VEC();

  int Cray_rounder();
  int colorscale();
  int greyscale();
  int prescale();
  int red_blue_scale();
  int zscale3();
  /* unlikely as it may seems, these routines work on the cray without 
   having to be called as COLRAM and SHADAM because the names are passed by me */
  int dist_colram_();
  int IWC_colram_();
  int shadam_();

  void Exit_gracefully(void);
  void print_usage();

  Boolean EXIT_FLAG = False;
  Boolean NETCDF_FORMAT;
  Boolean PAUSE;
  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *char_ptr_foo;
  char *time_buf_start;
  char cmdline[CMDLINE_SIZE];
  char err_file[FILESIZE];
  char in_file_1[FILESIZE];
  char in_file_2[FILESIZE];
  char out_file[FILESIZE];
  char line_label[LABEL_SIZE];
  
  extern char *optarg;
  extern int optind;
  
  float **IWC_snapshot;
  float **number_snapshot;
  float **layer_frame;
  float **layer_size;
  float **concentration;
  float **distribution;
  float **crystal_heat_LW;
  float **crystal_heat_SW;

  float * IWC_snapshot_ptr;
  float * number_snapshot_ptr;
  float * layer_frame_ptr;
  float * layer_size_ptr;
  float * concentration_ptr;
  float * distribution_ptr;
  float * crystal_heat_LW_ptr;
  float * crystal_heat_SW_ptr;

  float *rwrk;

  float *IWC_1;
  float *IWP_1;
  float *IWP_of_time_1;
  float *LW_cloud_forcing_1;
  float *SW_cloud_forcing_1;
  float *advective_heating_1;
  float *albedo_of_time_1;
  float *altitude_1;
  float *cloud_base_of_time_1;
  float *cloud_effective_radius_1;
  float *cloud_optical_depth_1;
  float *cloud_top_of_time_1;
  float *crystal_length_1;
  float *crystal_mass_1;
  float *delta_length_1;
  float *delta_mass_1;
  float *effective_length_1;
  float *effective_radius_1;
  float *emissivity_of_time_1;
  float *env_density_1;
  float *env_pressure_1;
  float *env_temperature_1;
  float *flux_down_LW_1;
  float *flux_down_SW_1;
  float *flux_up_LW_1;
  float *flux_up_SW_1;
  float *heating_rate_LW_1;
  float *heating_rate_SW_1;
  float *heating_rate_net_1;
  float *latent_heating_1;
  float *max_length_of_time_1;
  float *mmr_vapor_1;
  float *optical_depth_1;
  float *orig_env_temp_1;
  float *orig_mmr_vapor_1;
  float *potential_temp_1;
  float *pp_vapor_1;
  float *saturation_ice_of_time_1;
  float *saturation_ice_1;
  float *surface_area_of_time_1;
  float *time_array_1;
  float *time_snapshot_1;
  float *vapor_density_1;
  float *vapor_path_1;
  float *vapor_path_of_time_1;
  float *water_path_of_time_1;
  float *wind_speed_1;
  
  float *IWC_2;
  float *IWP_2;
  float *IWP_of_time_2;
  float *LW_cloud_forcing_2;
  float *SW_cloud_forcing_2;
  float *advective_heating_2;
  float *albedo_of_time_2;
  float *altitude_2;
  float *cloud_base_of_time_2;
  float *cloud_effective_radius_2;
  float *cloud_optical_depth_2;
  float *cloud_top_of_time_2;
  float *crystal_length_2;
  float *crystal_mass_2;
  float *delta_length_2;
  float *delta_mass_2;
  float *effective_length_2;
  float *effective_radius_2;
  float *emissivity_of_time_2;
  float *env_density_2;
  float *env_pressure_2;
  float *env_temperature_2;
  float *flux_down_LW_2;
  float *flux_down_SW_2;
  float *flux_up_LW_2;
  float *flux_up_SW_2;
  float *heating_rate_LW_2;
  float *heating_rate_SW_2;
  float *heating_rate_net_2;
  float *latent_heating_2;
  float *max_length_of_time_2;
  float *mmr_vapor_2;
  float *optical_depth_2;
  float *orig_env_temp_2;
  float *orig_mmr_vapor_2;
  float *potential_temp_2;
  float *pp_vapor_2;
  float *saturation_ice_of_time_2;
  float *saturation_ice_2;
  float *surface_area_of_time_2;
  float *time_array_2;
  float *time_snapshot_2;
  float *vapor_density_2;
  float *vapor_path_2;
  float *vapor_path_of_time_2;
  float *water_path_of_time_2;
  float *wind_speed_2;
  
  float abscissa_max;
  float abscissa_min;
  float dt;
  float dz;
  float cloud_effective_length;
  float float_foo2;
  float float_foo;
  float ord_height;
  float ordinate_max;
  float ordinate_min;
  float vpb;
  float vpl;
  float vpr;
  float vpt;
  float xaxis_max;
  float xaxis_min;
  float yaxis_max;
  float yaxis_min;
  
  int *IWC_crosshatch_lind;
  int *IWC_diff_grey_lind;
  int *dist_crosshatch_lind;
  int *dist_grey_lind;
  int *dist_diff_grey_lind;
  
  int *iam;
  int *iai;
  int *iag;
  float *xra;
  float *yra;
  int *iwrk;

  int agg_step;
  int rad_step;
  int plot_step;
  int iftp;
  int lam;
  int nra;
  int mai;
  int bufa;
  int conid;
  int frame;
  int ier;
  int int_foo;
  int jcrt,jsize;
  int kiwk;
  int krwk;
  int layer;
  int level;
  int num_band;
  int num_ccm2_level;
  int num_cloudy_layers;
  int num_contour_level;
  int num_frame;
  int num_layer;
  int num_rad_steps;
  int num_size;
  int num_step;
  int num_ticks;
  int opt;
  int size;
  int step;
  int tick;
  int wkid;
  int wtype;
  int xdim;
  int ydim;

  time_t clock;
  
  /* netCDF declarations */ 
  int cdfid_1,cdfid_2;  

  /* dimension sizes (which are longs) and array to hold them */ 
  long count[2];
  long num_layerp2,num_bandp2,num_framep1,num_sizep2,num_stepp1;

  /* indicial starting offsets (which are longs) */
  long start[2];
  
  /* scalar variable ids */
  int dt_id, dz_id, num_band_id, num_ccm2_level_id, num_cloudy_layers_id, num_frame_id, num_layer_id, num_size_id, num_step_id, agg_step_id, rad_step_id, plot_step_id, PAUSE_id;

  /* two dimensional variable ids */
  int IWC_snapshot_id, number_snapshot_id, concentration_id, distribution_id, crystal_heat_SW_id, crystal_heat_LW_id;

  /* one dimensional variable ids */
  int IWC_id, IWP_id, IWP_of_time_id, LW_cloud_forcing_id, SW_cloud_forcing_id, advective_heating_id, albedo_of_time_id, altitude_id, cloud_base_of_time_id, cloud_effective_radius_id, cloud_top_of_time_id, crystal_length_id, crystal_mass_id, delta_length_id, delta_mass_id, effective_radius_id, emissivity_of_time_id, env_density_id, env_pressure_id, env_temperature_id, flux_down_LW_id, flux_down_SW_id, flux_up_LW_id, flux_up_SW_id, heating_rate_LW_id, heating_rate_SW_id, heating_rate_net_id, latent_heating_id, max_length_of_time_id, mmr_vapor_id, cloud_optical_depth_id, orig_env_temp_id, orig_mmr_vapor_id, potential_temp_id, pp_vapor_id, saturation_ice_of_time_id, saturation_ice_id, surface_area_of_time_id, time_array_id, time_snapshot_id, vapor_density_id, vapor_path_id, vapor_path_of_time_id, optical_depth_id, water_path_of_time_id, wind_speed_id;

  /* for deciding when runs have different # sizes */ 
  Boolean PAUSE_2;
  char cmdline_2[CMDLINE_SIZE];
  float dt_2;
  float dz_2;
  int num_band_2;
  int num_ccm2_level_2;
  int num_cloudy_layers_2;
  int num_frame_2;
  int num_layer_2;
  int num_size_2;
  int num_step_2;
  int agg_step_2;
  int rad_step_2;
  int plot_step_2;
  
  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = False; /* Option I */
  STDOUT = False; /* Option O */
  NETCDF_FORMAT = True;
  VERBOSE = True; /* print out WARNINGS? */ /* Option V */

  (void)strcpy(err_file,"stderr"); /* Option e */
  (void)strcpy(in_file_1,"cloud.nc"); /* Option i */
  (void)strcpy(in_file_2,"cloud.nc"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  
  /* Set GKS variables */
  /* note that 300000 and 30000 are enough */ 
  lam = 500000; /* size of NCAR iam array */ 
  nra = 50000; /* size of NCAR xra and yra arrays */ 
  mai = 10; /* size of NCAR iai and iag arrays */ 

  conid = 2;
  ier = 6;
  jcrt = 10;
  jsize = 12;
  kiwk = 1000; /* size of NCAR CONPACK int workspace */
  krwk = 5000; /* size of NCAR CONPACK real workspace */
  wkid = 1;
  wtype = 1;

  /* parse command line arguments */
  while((opt = getopt(argc, argv,"D:d:Ee:Ii:j:Oo:Vv"))!= EOF){
    switch(opt){
    case 'D':
      /* The debugging level.  Default is 0 */
      debug = atoi(optarg);
      break;
    case 'd':
      /* the spare int input for debugging. Default is 0 */
      debug_value = atoi(optarg);
      break;
    case 'E':
      /* Toggle the error file stream. Default is True */
      STDERR = !STDERR;
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_file,optarg);
      break;
    case 'I':
      /* Toggle the input file stream. Default is False */
      STDIN = !STDIN;
      break;
    case 'i':
      /* get the input file name of the first. Default is cloud.nc */
      (void)strcpy(in_file_1,optarg);
      break;
    case 'j':
      /* get the input file name of the second file. Default is cloud.nc */
      (void)strcpy(in_file_2,optarg);
      break;
    case 'O':
      /* Toggle the output file stream. Default is True */
      STDOUT = !STDOUT;
      break;
    case 'o':
      /* get the output file name. Default is stdout */
      (void)strcpy(out_file,optarg);
      break;
    case 'v':
      /* print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: io_diff.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/cld/io_diff.c,v $\n");
      (void)fprintf(stderr,"$Id$\n");
      (void)fprintf(stderr,"$State: Exp $\n");
      exit(0);
      break;
    case 'V':
      /* toggle verbose printing out of WARNINGS. Default is True */ 
      VERBOSE=!VERBOSE;
      break;
    case '?':
      /* print proper usage */
      (void)print_usage();
      exit(1);
    } /* end switch */
  } /* end while loop */
  
  /* start the clock */  
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stderr,"\tstart = %s",time_buf_start);
  
  if(STDOUT){
    fp_out = stdout;
  }else{
    if( (fp_out = fopen( out_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening output file %s\n",out_file);
      exit(1);
    } /* end if */
  } /* end if */ 
  
  if(STDERR){
    fp_err = stderr;
  }else{
    if( (fp_err = fopen( err_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening error file %s\n",err_file);
      exit(1);
    } /* end if */
  } /* end if */
  
  /* allocate storage for all the non-input arrays, i.e. GKS stuff */
  if(
     ((IWC_crosshatch_lind=(int *)malloc((NUM_IWC_CONTOURS)*sizeof(int))) == NULL ) ||
     ((IWC_diff_grey_lind=(int *)malloc((NUM_IWC_DIFF_CONTOURS)*sizeof(int))) == NULL ) || 
     ((dist_crosshatch_lind=(int *)malloc((NUM_DIST_CONTOURS)*sizeof(int))) == NULL ) || 
     ((dist_grey_lind=(int *)malloc((NUM_DIST_CONTOURS)*sizeof(int))) == NULL ) || 
     ((dist_diff_grey_lind=(int *)malloc((NUM_DIST_DIFF_CONTOURS)*sizeof(int))) == NULL ) ||
     ((iam=(int *)malloc((lam)*sizeof(int))) == NULL ) ||
     ((iai=(int *)malloc((mai)*sizeof(int))) == NULL ) ||
     ((iag=(int *)malloc((mai)*sizeof(int))) == NULL ) ||
     ((xra=(float *)malloc((nra)*sizeof(float))) == NULL ) ||
     ((yra=(float *)malloc((nra)*sizeof(float))) == NULL ) ||
     ((iwrk=(int *)malloc((kiwk)*sizeof(int))) == NULL ) ||
     ((rwrk=(float *)malloc((krwk)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  /* set the contour index arrays */ 
  for(level=0;level<NUM_IWC_CONTOURS;level++){
    IWC_crosshatch_lind[level] = level+1;
  } /* end loop over levels */
  for(level=0;level<NUM_IWC_DIFF_CONTOURS;level++){
    IWC_diff_grey_lind[level] = level+2;
  } /* end loop over levels */
  for(level=0;level<NUM_DIST_CONTOURS;level++){
    dist_crosshatch_lind[level] = level+1;
  } /* end loop over levels */
  for(level=0;level<NUM_DIST_CONTOURS;level++){
    dist_grey_lind[level] = level+2;
  } /* end loop over levels */
  for(level=0;level<NUM_DIST_DIFF_CONTOURS;level++){
    dist_diff_grey_lind[level] = level+2;
  } /* end loop over levels */

  if(NETCDF_FORMAT){
    /* first get the scalar values and dimensions for malloc'ing */ 
    start[0]=0L;
    cdfid_1 = ncopen(in_file_1, NC_NOWRITE);

    ncattget(cdfid_1,NC_GLOBAL,"cmdline",(void *)cmdline);
    dt_id=ncvarid(cdfid_1,"dt");
    ncvarget1(cdfid_1,dt_id,start,(void *)&dt);
    dz_id=ncvarid(cdfid_1,"dz");
    ncvarget1(cdfid_1,dz_id,start,(void *)&dz);
    num_band_id=ncvarid(cdfid_1,"num_band");
    ncvarget1(cdfid_1,num_band_id,start,(void *)&num_band);
    num_ccm2_level_id=ncvarid(cdfid_1,"num_ccm2_level");
    ncvarget1(cdfid_1,num_ccm2_level_id,start,(void *)&num_ccm2_level);
    num_cloudy_layers_id=ncvarid(cdfid_1,"num_cloudy_layers");
    ncvarget1(cdfid_1,num_cloudy_layers_id,start,(void *)&num_cloudy_layers);
    num_frame_id=ncvarid(cdfid_1,"num_frame");
    ncvarget1(cdfid_1,num_frame_id,start,(void *)&num_frame);
    num_layer_id=ncvarid(cdfid_1,"num_layer");
    ncvarget1(cdfid_1,num_layer_id,start,(void *)&num_layer);
    num_size_id=ncvarid(cdfid_1,"num_size");
    ncvarget1(cdfid_1,num_size_id,start,(void *)&num_size);
    num_step_id=ncvarid(cdfid_1,"num_step");
    ncvarget1(cdfid_1,num_step_id,start,(void *)&num_step);
    agg_step_id=ncvarid(cdfid_1,"agg_step");
    ncvarget1(cdfid_1,agg_step_id,start,(void *)&agg_step);
    rad_step_id=ncvarid(cdfid_1,"rad_step");
    ncvarget1(cdfid_1,rad_step_id,start,(void *)&rad_step);
    plot_step_id=ncvarid(cdfid_1,"plot_step");
    ncvarget1(cdfid_1,plot_step_id,start,(void *)&plot_step);
    PAUSE_id=ncvarid(cdfid_1,"PAUSE");
    ncvarget1(cdfid_1,PAUSE_id,start,(void *)&PAUSE);
  } /* end if NETCDF_FORMAT */ 
  
  if(debug == 13){
    (void)fprintf(fp_err,"\n");
    (void)fprintf(fp_err,"Command Line = %s\n",cmdline);
    (void)fprintf(fp_err,"dt = %g\n",dt);
    (void)fprintf(fp_err,"dz = %g\n",dz);
    (void)fprintf(fp_err,"num_band = %i\n",num_band);
    (void)fprintf(fp_err,"num_ccm2_level = %i\n",num_ccm2_level);
    (void)fprintf(fp_err,"num_cloudy_layers = %i\n",num_cloudy_layers);
    (void)fprintf(fp_err,"num_frame = %i\n",num_frame);
    (void)fprintf(fp_err,"num_layer = %i\n",num_layer);
    (void)fprintf(fp_err,"num_size = %i\n",num_size);
    (void)fprintf(fp_err,"num_step = %i\n",num_step);
    (void)fprintf(fp_err,"agg_step = %i\n",agg_step);
    (void)fprintf(fp_err,"rad_step = %i\n",rad_step);
    (void)fprintf(fp_err,"plot_step = %i\n",plot_step);
    (void)fprintf(fp_err,"PAUSE = %i\n",PAUSE);
    (void)fprintf(fp_err,"\n");
  } /* end debug */

  /* allocate storage for all the dynamic input arrays from the first file */
  if(
     ((IWC_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((LW_cloud_forcing_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((SW_cloud_forcing_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((advective_heating_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((albedo_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((altitude_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_base_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_effective_radius_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_top_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((crystal_length_1=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_mass_1=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_length_1=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_mass_1=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((effective_length_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((effective_radius_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((emissivity_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((env_density_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_temperature_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_LW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_SW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_LW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_SW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_LW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_SW_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_net_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((latent_heating_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((max_length_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((mmr_vapor_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_optical_depth_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((orig_env_temp_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((orig_mmr_vapor_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((potential_temp_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((pp_vapor_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((saturation_ice_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((saturation_ice_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((surface_area_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_array_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_snapshot_1=(float *)malloc((num_frame+1)*sizeof(float))) == NULL ) ||
     ((vapor_density_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((optical_depth_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((water_path_of_time_1=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((wind_speed_1=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  /* allocate storage for all the dynamic input arrays from the second file */
  if(
     ((IWC_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((LW_cloud_forcing_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((SW_cloud_forcing_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((advective_heating_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((albedo_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((altitude_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_base_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_effective_radius_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_top_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((crystal_length_2=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_mass_2=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_length_2=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_mass_2=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((effective_length_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((effective_radius_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((emissivity_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((env_density_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_temperature_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_LW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_SW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_LW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_SW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_LW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_SW_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_net_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((latent_heating_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((max_length_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((mmr_vapor_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_optical_depth_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((orig_env_temp_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((orig_mmr_vapor_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((potential_temp_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((pp_vapor_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((saturation_ice_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((saturation_ice_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((surface_area_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_array_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_snapshot_2=(float *)malloc((num_frame+1)*sizeof(float))) == NULL ) ||
     ((vapor_density_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((optical_depth_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((water_path_of_time_2=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((wind_speed_2=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  int_foo=(num_layer > num_size) ? num_layer : num_size;
  int_foo=(int_foo > num_step+1) ? int_foo : num_step+1;
  if(
     ((abscissa=(float *)malloc((int_foo)*sizeof(float))) == NULL ) ||
     ((ordinate=(float *)malloc((int_foo)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  if(
     ((IWC_snapshot_ptr=
       (float *)malloc((num_frame+1)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((number_snapshot_ptr=
       (float *)malloc((num_frame+1)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((layer_frame_ptr=
       (float *)malloc((num_layer)*(num_frame+1)*sizeof(float))) == NULL ) ||
     ((layer_size_ptr=
       (float *)malloc((num_layer)*(num_size)*sizeof(float))) == NULL ) ||
     ((concentration_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((distribution_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_heat_LW_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_heat_SW_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  if(
     ((IWC_snapshot=
       (float **)malloc((num_frame+1)*sizeof(float *))) == NULL ) ||
     ((number_snapshot=
       (float **)malloc((num_frame+1)*sizeof(float *))) == NULL ) ||
     ((layer_frame=
       (float **)malloc((num_layer)*sizeof(float *))) == NULL ) ||
     ((layer_size=
       (float **)malloc((num_layer)*sizeof(float *))) == NULL ) ||
     ((concentration=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((distribution=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((crystal_heat_LW=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((crystal_heat_SW=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  /* REMEMBER: these matrices are indexed as [layer][size] */
  for(layer=0;layer<=num_layer+1;layer++){
    concentration[layer]=concentration_ptr+(num_size+2)*layer;    
    distribution[layer]=distribution_ptr+(num_size+2)*layer;    
    crystal_heat_LW[layer]=crystal_heat_LW_ptr+(num_size+2)*layer;    
    crystal_heat_SW[layer]=crystal_heat_SW_ptr+(num_size+2)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these are the only snug matrices indexed as [layer][size] */
  for(layer=0;layer<=num_layer-1;layer++){
    layer_size[layer]=layer_size_ptr+(num_size)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these are the only snug matrices indexed as [layer][frame] */
  for(layer=0;layer<=num_layer-1;layer++){
    layer_frame[layer]=layer_frame_ptr+(num_frame+1)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these are the only matrices indexed as [frame][layer] */
  for(frame=0;frame<=num_frame;frame++){
    IWC_snapshot[frame]=IWC_snapshot_ptr+(num_layer+2)*frame;    
    number_snapshot[frame]=number_snapshot_ptr+(num_layer+2)*frame;    
  } /* end loop over frames */
  
  /* Now that appropriate space has been allocated, read in the arrays 
   from the first file */
  if(NETCDF_FORMAT){
    /* retrieve the one-dimensional arays */ 
    start[0]=0L;

    /* define sizes of dimensions */
    num_layerp2 = (long)(num_layer+2);
    num_bandp2 = (long)(num_band+2);
    num_framep1 = (long)(num_frame+1);
    num_sizep2 = (long)(num_size+2);
    num_stepp1 = (long)(num_step+1);

    IWC_id=ncvarid(cdfid_1,"IWC");
    ncvarget(cdfid_1,IWC_id,start,&num_layerp2,(void *)IWC_1);
    IWP_id=ncvarid(cdfid_1,"IWP");
    ncvarget(cdfid_1,IWP_id,start,&num_layerp2,(void *)IWP_1);
    IWP_of_time_id=ncvarid(cdfid_1,"IWP_of_time");
    ncvarget(cdfid_1,IWP_of_time_id,start,&num_stepp1,(void *)IWP_of_time_1);
    LW_cloud_forcing_id=ncvarid(cdfid_1,"LW_cloud_forcing");
    ncvarget(cdfid_1,LW_cloud_forcing_id,start,&num_stepp1,(void *)LW_cloud_forcing_1);
    SW_cloud_forcing_id=ncvarid(cdfid_1,"SW_cloud_forcing");
    ncvarget(cdfid_1,SW_cloud_forcing_id,start,&num_stepp1,(void *)SW_cloud_forcing_1);
    advective_heating_id=ncvarid(cdfid_1,"advective_heating");
    ncvarget(cdfid_1,advective_heating_id,start,&num_layerp2,(void *)advective_heating_1);
    albedo_of_time_id=ncvarid(cdfid_1,"albedo_of_time");
    ncvarget(cdfid_1,albedo_of_time_id,start,&num_stepp1,(void *)albedo_of_time_1);
    altitude_id=ncvarid(cdfid_1,"altitude");
    ncvarget(cdfid_1,altitude_id,start,&num_layerp2,(void *)altitude_1);
    cloud_base_of_time_id=ncvarid(cdfid_1,"cloud_base_of_time");
    ncvarget(cdfid_1,cloud_base_of_time_id,start,&num_stepp1,(void *)cloud_base_of_time_1);
    cloud_effective_radius_id=ncvarid(cdfid_1,"cloud_effective_radius");
    ncvarget(cdfid_1,cloud_effective_radius_id,start,&num_stepp1,(void *)cloud_effective_radius_1);
    cloud_top_of_time_id=ncvarid(cdfid_1,"cloud_top_of_time");
    ncvarget(cdfid_1,cloud_top_of_time_id,start,&num_stepp1,(void *)cloud_top_of_time_1);
    crystal_length_id=ncvarid(cdfid_1,"crystal_length");
    ncvarget(cdfid_1,crystal_length_id,start,&num_sizep2,(void *)crystal_length_1);
    crystal_mass_id=ncvarid(cdfid_1,"crystal_mass");
    ncvarget(cdfid_1,crystal_mass_id,start,&num_sizep2,(void *)crystal_mass_1);
    delta_length_id=ncvarid(cdfid_1,"delta_length");
    ncvarget(cdfid_1,delta_length_id,start,&num_sizep2,(void *)delta_length_1);
    delta_mass_id=ncvarid(cdfid_1,"delta_mass");
    ncvarget(cdfid_1,delta_mass_id,start,&num_sizep2,(void *)delta_mass_1);
    effective_radius_id=ncvarid(cdfid_1,"effective_radius");
    ncvarget(cdfid_1,effective_radius_id,start,&num_layerp2,(void *)effective_radius_1);
    emissivity_of_time_id=ncvarid(cdfid_1,"emissivity_of_time");
    ncvarget(cdfid_1,emissivity_of_time_id,start,&num_stepp1,(void *)emissivity_of_time_1);
    env_density_id=ncvarid(cdfid_1,"env_density");
    ncvarget(cdfid_1,env_density_id,start,&num_layerp2,(void *)env_density_1);
    env_pressure_id=ncvarid(cdfid_1,"env_pressure");
    ncvarget(cdfid_1,env_pressure_id,start,&num_layerp2,(void *)env_pressure_1);
    env_temperature_id=ncvarid(cdfid_1,"env_temperature");
    ncvarget(cdfid_1,env_temperature_id,start,&num_layerp2,(void *)env_temperature_1);
    flux_down_LW_id=ncvarid(cdfid_1,"flux_down_LW");
    ncvarget(cdfid_1,flux_down_LW_id,start,&num_layerp2,(void *)flux_down_LW_1);
    flux_down_SW_id=ncvarid(cdfid_1,"flux_down_SW");
    ncvarget(cdfid_1,flux_down_SW_id,start,&num_layerp2,(void *)flux_down_SW_1);
    flux_up_LW_id=ncvarid(cdfid_1,"flux_up_LW");
    ncvarget(cdfid_1,flux_up_LW_id,start,&num_layerp2,(void *)flux_up_LW_1);
    flux_up_SW_id=ncvarid(cdfid_1,"flux_up_SW");
    ncvarget(cdfid_1,flux_up_SW_id,start,&num_layerp2,(void *)flux_up_SW_1);
    heating_rate_LW_id=ncvarid(cdfid_1,"heating_rate_LW");
    ncvarget(cdfid_1,heating_rate_LW_id,start,&num_layerp2,(void *)heating_rate_LW_1);
    heating_rate_SW_id=ncvarid(cdfid_1,"heating_rate_SW");
    ncvarget(cdfid_1,heating_rate_SW_id,start,&num_layerp2,(void *)heating_rate_SW_1);
    heating_rate_net_id=ncvarid(cdfid_1,"heating_rate_net");
    ncvarget(cdfid_1,heating_rate_net_id,start,&num_layerp2,(void *)heating_rate_net_1);
    latent_heating_id=ncvarid(cdfid_1,"latent_heating");
    ncvarget(cdfid_1,latent_heating_id,start,&num_layerp2,(void *)latent_heating_1);
    max_length_of_time_id=ncvarid(cdfid_1,"max_length_of_time");
    ncvarget(cdfid_1,max_length_of_time_id,start,&num_stepp1,(void *)max_length_of_time_1);
    mmr_vapor_id=ncvarid(cdfid_1,"mmr_vapor");
    ncvarget(cdfid_1,mmr_vapor_id,start,&num_layerp2,(void *)mmr_vapor_1);
    cloud_optical_depth_id=ncvarid(cdfid_1,"cloud_optical_depth");
    ncvarget(cdfid_1,cloud_optical_depth_id,start,&num_stepp1,(void *)cloud_optical_depth_1);
    orig_env_temp_id=ncvarid(cdfid_1,"orig_env_temp");
    ncvarget(cdfid_1,orig_env_temp_id,start,&num_layerp2,(void *)orig_env_temp_1);
    orig_mmr_vapor_id=ncvarid(cdfid_1,"orig_mmr_vapor");
    ncvarget(cdfid_1,orig_mmr_vapor_id,start,&num_layerp2,(void *)orig_mmr_vapor_1);
    potential_temp_id=ncvarid(cdfid_1,"potential_temp");
    ncvarget(cdfid_1,potential_temp_id,start,&num_layerp2,(void *)potential_temp_1);
    pp_vapor_id=ncvarid(cdfid_1,"pp_vapor");
    ncvarget(cdfid_1,pp_vapor_id,start,&num_layerp2,(void *)pp_vapor_1);
    saturation_ice_of_time_id=ncvarid(cdfid_1,"saturation_ice_of_time");
    ncvarget(cdfid_1,saturation_ice_of_time_id,start,&num_stepp1,(void *)saturation_ice_of_time_1);
    saturation_ice_id=ncvarid(cdfid_1,"saturation_ice");
    ncvarget(cdfid_1,saturation_ice_id,start,&num_layerp2,(void *)saturation_ice_1);
    surface_area_of_time_id=ncvarid(cdfid_1,"surface_area_of_time");
    ncvarget(cdfid_1,surface_area_of_time_id,start,&num_stepp1,(void *)surface_area_of_time_1);
    time_array_id=ncvarid(cdfid_1,"time_array");
    ncvarget(cdfid_1,time_array_id,start,&num_stepp1,(void *)time_array_1);
    time_snapshot_id=ncvarid(cdfid_1,"time_snapshot");
    ncvarget(cdfid_1,time_snapshot_id,start,&num_framep1,(void *)time_snapshot_1);
    vapor_density_id=ncvarid(cdfid_1,"vapor_density");
    ncvarget(cdfid_1,vapor_density_id,start,&num_layerp2,(void *)vapor_density_1);
    vapor_path_id=ncvarid(cdfid_1,"vapor_path");
    ncvarget(cdfid_1,vapor_path_id,start,&num_layerp2,(void *)vapor_path_1);
    vapor_path_of_time_id=ncvarid(cdfid_1,"vapor_path_of_time");
    ncvarget(cdfid_1,vapor_path_of_time_id,start,&num_stepp1,(void *)vapor_path_of_time_1);
    optical_depth_id=ncvarid(cdfid_1,"optical_depth");
    ncvarget(cdfid_1,optical_depth_id,start,&num_layerp2,(void *)optical_depth_1);
    water_path_of_time_id=ncvarid(cdfid_1,"water_path_of_time");
    ncvarget(cdfid_1,water_path_of_time_id,start,&num_stepp1,(void *)water_path_of_time_1);
    wind_speed_id=ncvarid(cdfid_1,"wind_speed");
    ncvarget(cdfid_1,wind_speed_id,start,&num_layerp2,(void *)wind_speed_1);

    /* now get the two-dimensional arays */ 
    start[0]=0L;
    start[1]=0L;

    count[0] = num_framep1;
    count[1] = num_layerp2;
    IWC_snapshot_id=ncvarid(cdfid_1,"IWC_snapshot");
    ncvarget(cdfid_1,IWC_snapshot_id,start,count,(void *)IWC_snapshot_ptr);
    count[0] = num_framep1;
    count[1] = num_layerp2;
    number_snapshot_id=ncvarid(cdfid_1,"number_snapshot");
    ncvarget(cdfid_1,number_snapshot_id,start,count,(void *)number_snapshot_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid_1,"concentration");
    ncvarget(cdfid_1,concentration_id,start,count,(void *)concentration_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    distribution_id=ncvarid(cdfid_1,"distribution");
    ncvarget(cdfid_1,distribution_id,start,count,(void *)distribution_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    crystal_heat_SW_id=ncvarid(cdfid_1,"crystal_heat_SW");
    ncvarget(cdfid_1,crystal_heat_SW_id,start,count,(void *)crystal_heat_SW_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    crystal_heat_LW_id=ncvarid(cdfid_1,"crystal_heat_LW");
    ncvarget(cdfid_1,crystal_heat_LW_id,start,count,(void *)crystal_heat_LW_ptr);

    (void)fprintf(fp_err,"Reading first file's netCDF format cloud data from %s\n",in_file_1);
  } /* end if NETCDF_FORMAT */
  
  /* now read in and get the arrays from the second file */ 
  if(NETCDF_FORMAT){
    /* first get the scalar values and dimensions for malloc'ing */ 
    start[0]=0L;
    cdfid_2=ncopen(in_file_2, NC_NOWRITE);

    ncattget(cdfid_2,NC_GLOBAL,"cmdline",(void *)cmdline_2);
    dt_id=ncvarid(cdfid_2,"dt");
    ncvarget1(cdfid_2,dt_id,start,(void *)&dt_2);
    dz_id=ncvarid(cdfid_2,"dz");
    ncvarget1(cdfid_2,dz_id,start,(void *)&dz_2);
    num_band_id=ncvarid(cdfid_2,"num_band");
    ncvarget1(cdfid_2,num_band_id,start,(void *)&num_band_2);
    num_ccm2_level_id=ncvarid(cdfid_2,"num_ccm2_level");
    ncvarget1(cdfid_2,num_ccm2_level_id,start,(void *)&num_ccm2_level_2);
    num_cloudy_layers_id=ncvarid(cdfid_2,"num_cloudy_layers");
    ncvarget1(cdfid_2,num_cloudy_layers_id,start,(void *)&num_cloudy_layers_2);
    num_frame_id=ncvarid(cdfid_2,"num_frame");
    ncvarget1(cdfid_2,num_frame_id,start,(void *)&num_frame_2);
    num_layer_id=ncvarid(cdfid_2,"num_layer");
    ncvarget1(cdfid_2,num_layer_id,start,(void *)&num_layer_2);
    num_size_id=ncvarid(cdfid_2,"num_size");
    ncvarget1(cdfid_2,num_size_id,start,(void *)&num_size_2);
    num_step_id=ncvarid(cdfid_2,"num_step");
    ncvarget1(cdfid_2,num_step_id,start,(void *)&num_step_2);
    agg_step_id=ncvarid(cdfid_2,"agg_step");
    ncvarget1(cdfid_2,agg_step_id,start,(void *)&agg_step_2);
    rad_step_id=ncvarid(cdfid_2,"rad_step");
    ncvarget1(cdfid_2,rad_step_id,start,(void *)&rad_step_2);
    plot_step_id=ncvarid(cdfid_2,"plot_step");
    ncvarget1(cdfid_2,plot_step_id,start,(void *)&plot_step_2);
    PAUSE_id=ncvarid(cdfid_2,"PAUSE");
    ncvarget1(cdfid_2,PAUSE_id,start,(void *)&PAUSE_2);
  } /* end if NETCDF_FORMAT */ 
  
  if(debug == 13){
    (void)fprintf(fp_err,"\n");
    (void)fprintf(fp_err,"Command Line = %s\n",cmdline);
    (void)fprintf(fp_err,"dt = %g\n",dt);
    (void)fprintf(fp_err,"dz = %g\n",dz);
    (void)fprintf(fp_err,"num_band = %i\n",num_band);
    (void)fprintf(fp_err,"num_ccm2_level = %i\n",num_ccm2_level);
    (void)fprintf(fp_err,"num_cloudy_layers = %i\n",num_cloudy_layers);
    (void)fprintf(fp_err,"num_frame = %i\n",num_frame);
    (void)fprintf(fp_err,"num_layer = %i\n",num_layer);
    (void)fprintf(fp_err,"num_size = %i\n",num_size);
    (void)fprintf(fp_err,"num_step = %i\n",num_step);
    (void)fprintf(fp_err,"agg_step = %i\n",agg_step);
    (void)fprintf(fp_err,"rad_step = %i\n",rad_step);
    (void)fprintf(fp_err,"plot_step = %i\n",plot_step);
    (void)fprintf(fp_err,"PAUSE = %i\n",PAUSE);
    (void)fprintf(fp_err,"\n");
  } /* end debug */

  /* Now that appropriate space has been allocated, read in the arrays */
  if(NETCDF_FORMAT){
    /* retrieve the one-dimensional arrays */ 
    start[0]=0L;

    /* define sizes of dimensions */
    num_layerp2 = (long)(num_layer+2);
    num_bandp2 = (long)(num_band+2);
    num_framep1 = (long)(num_frame+1);
    num_sizep2 = (long)(num_size+2);
    num_stepp1 = (long)(num_step+1);

    IWC_id=ncvarid(cdfid_2,"IWC");
    ncvarget(cdfid_2,IWC_id,start,&num_layerp2,(void *)IWC_2);
    IWP_id=ncvarid(cdfid_2,"IWP");
    ncvarget(cdfid_2,IWP_id,start,&num_layerp2,(void *)IWP_2);
    IWP_of_time_id=ncvarid(cdfid_2,"IWP_of_time");
    ncvarget(cdfid_2,IWP_of_time_id,start,&num_stepp1,(void *)IWP_of_time_2);
    LW_cloud_forcing_id=ncvarid(cdfid_2,"LW_cloud_forcing");
    ncvarget(cdfid_2,LW_cloud_forcing_id,start,&num_stepp1,(void *)LW_cloud_forcing_2);
    SW_cloud_forcing_id=ncvarid(cdfid_2,"SW_cloud_forcing");
    ncvarget(cdfid_2,SW_cloud_forcing_id,start,&num_stepp1,(void *)SW_cloud_forcing_2);
    advective_heating_id=ncvarid(cdfid_2,"advective_heating");
    ncvarget(cdfid_2,advective_heating_id,start,&num_layerp2,(void *)advective_heating_2);
    albedo_of_time_id=ncvarid(cdfid_2,"albedo_of_time");
    ncvarget(cdfid_2,albedo_of_time_id,start,&num_stepp1,(void *)albedo_of_time_2);
    altitude_id=ncvarid(cdfid_2,"altitude");
    ncvarget(cdfid_2,altitude_id,start,&num_layerp2,(void *)altitude_2);
    cloud_base_of_time_id=ncvarid(cdfid_2,"cloud_base_of_time");
    ncvarget(cdfid_2,cloud_base_of_time_id,start,&num_stepp1,(void *)cloud_base_of_time_2);
    cloud_effective_radius_id=ncvarid(cdfid_2,"cloud_effective_radius");
    ncvarget(cdfid_2,cloud_effective_radius_id,start,&num_stepp1,(void *)cloud_effective_radius_2);
    cloud_top_of_time_id=ncvarid(cdfid_2,"cloud_top_of_time");
    ncvarget(cdfid_2,cloud_top_of_time_id,start,&num_stepp1,(void *)cloud_top_of_time_2);
    crystal_length_id=ncvarid(cdfid_2,"crystal_length");
    ncvarget(cdfid_2,crystal_length_id,start,&num_sizep2,(void *)crystal_length_2);
    crystal_mass_id=ncvarid(cdfid_2,"crystal_mass");
    ncvarget(cdfid_2,crystal_mass_id,start,&num_sizep2,(void *)crystal_mass_2);
    delta_length_id=ncvarid(cdfid_2,"delta_length");
    ncvarget(cdfid_2,delta_length_id,start,&num_sizep2,(void *)delta_length_2);
    delta_mass_id=ncvarid(cdfid_2,"delta_mass");
    ncvarget(cdfid_2,delta_mass_id,start,&num_sizep2,(void *)delta_mass_2);
    effective_radius_id=ncvarid(cdfid_2,"effective_radius");
    ncvarget(cdfid_2,effective_radius_id,start,&num_layerp2,(void *)effective_radius_2);
    emissivity_of_time_id=ncvarid(cdfid_2,"emissivity_of_time");
    ncvarget(cdfid_2,emissivity_of_time_id,start,&num_stepp1,(void *)emissivity_of_time_2);
    env_density_id=ncvarid(cdfid_2,"env_density");
    ncvarget(cdfid_2,env_density_id,start,&num_layerp2,(void *)env_density_2);
    env_pressure_id=ncvarid(cdfid_2,"env_pressure");
    ncvarget(cdfid_2,env_pressure_id,start,&num_layerp2,(void *)env_pressure_2);
    env_temperature_id=ncvarid(cdfid_2,"env_temperature");
    ncvarget(cdfid_2,env_temperature_id,start,&num_layerp2,(void *)env_temperature_2);
    flux_down_LW_id=ncvarid(cdfid_2,"flux_down_LW");
    ncvarget(cdfid_2,flux_down_LW_id,start,&num_layerp2,(void *)flux_down_LW_2);
    flux_down_SW_id=ncvarid(cdfid_2,"flux_down_SW");
    ncvarget(cdfid_2,flux_down_SW_id,start,&num_layerp2,(void *)flux_down_SW_2);
    flux_up_LW_id=ncvarid(cdfid_2,"flux_up_LW");
    ncvarget(cdfid_2,flux_up_LW_id,start,&num_layerp2,(void *)flux_up_LW_2);
    flux_up_SW_id=ncvarid(cdfid_2,"flux_up_SW");
    ncvarget(cdfid_2,flux_up_SW_id,start,&num_layerp2,(void *)flux_up_SW_2);
    heating_rate_LW_id=ncvarid(cdfid_2,"heating_rate_LW");
    ncvarget(cdfid_2,heating_rate_LW_id,start,&num_layerp2,(void *)heating_rate_LW_2);
    heating_rate_SW_id=ncvarid(cdfid_2,"heating_rate_SW");
    ncvarget(cdfid_2,heating_rate_SW_id,start,&num_layerp2,(void *)heating_rate_SW_2);
    heating_rate_net_id=ncvarid(cdfid_2,"heating_rate_net");
    ncvarget(cdfid_2,heating_rate_net_id,start,&num_layerp2,(void *)heating_rate_net_2);
    latent_heating_id=ncvarid(cdfid_2,"latent_heating");
    ncvarget(cdfid_2,latent_heating_id,start,&num_layerp2,(void *)latent_heating_2);
    max_length_of_time_id=ncvarid(cdfid_2,"max_length_of_time");
    ncvarget(cdfid_2,max_length_of_time_id,start,&num_stepp1,(void *)max_length_of_time_2);
    mmr_vapor_id=ncvarid(cdfid_2,"mmr_vapor");
    ncvarget(cdfid_2,mmr_vapor_id,start,&num_layerp2,(void *)mmr_vapor_2);
    cloud_optical_depth_id=ncvarid(cdfid_2,"cloud_optical_depth");
    ncvarget(cdfid_2,cloud_optical_depth_id,start,&num_stepp1,(void *)cloud_optical_depth_2);
    orig_env_temp_id=ncvarid(cdfid_2,"orig_env_temp");
    ncvarget(cdfid_2,orig_env_temp_id,start,&num_layerp2,(void *)orig_env_temp_2);
    orig_mmr_vapor_id=ncvarid(cdfid_2,"orig_mmr_vapor");
    ncvarget(cdfid_2,orig_mmr_vapor_id,start,&num_layerp2,(void *)orig_mmr_vapor_2);
    potential_temp_id=ncvarid(cdfid_2,"potential_temp");
    ncvarget(cdfid_2,potential_temp_id,start,&num_layerp2,(void *)potential_temp_2);
    pp_vapor_id=ncvarid(cdfid_2,"pp_vapor");
    ncvarget(cdfid_2,pp_vapor_id,start,&num_layerp2,(void *)pp_vapor_2);
    saturation_ice_of_time_id=ncvarid(cdfid_2,"saturation_ice_of_time");
    ncvarget(cdfid_2,saturation_ice_of_time_id,start,&num_stepp1,(void *)saturation_ice_of_time_2);
    saturation_ice_id=ncvarid(cdfid_2,"saturation_ice");
    ncvarget(cdfid_2,saturation_ice_id,start,&num_layerp2,(void *)saturation_ice_2);
    surface_area_of_time_id=ncvarid(cdfid_2,"surface_area_of_time");
    ncvarget(cdfid_2,surface_area_of_time_id,start,&num_stepp1,(void *)surface_area_of_time_2);
    time_array_id=ncvarid(cdfid_2,"time_array");
    ncvarget(cdfid_2,time_array_id,start,&num_stepp1,(void *)time_array_2);
    time_snapshot_id=ncvarid(cdfid_2,"time_snapshot");
    ncvarget(cdfid_2,time_snapshot_id,start,&num_framep1,(void *)time_snapshot_2);
    vapor_density_id=ncvarid(cdfid_2,"vapor_density");
    ncvarget(cdfid_2,vapor_density_id,start,&num_layerp2,(void *)vapor_density_2);
    vapor_path_id=ncvarid(cdfid_2,"vapor_path");
    ncvarget(cdfid_2,vapor_path_id,start,&num_layerp2,(void *)vapor_path_2);
    vapor_path_of_time_id=ncvarid(cdfid_2,"vapor_path_of_time");
    ncvarget(cdfid_2,vapor_path_of_time_id,start,&num_stepp1,(void *)vapor_path_of_time_2);
    optical_depth_id=ncvarid(cdfid_2,"optical_depth");
    ncvarget(cdfid_2,optical_depth_id,start,&num_layerp2,(void *)optical_depth_2);
    water_path_of_time_id=ncvarid(cdfid_2,"water_path_of_time");
    ncvarget(cdfid_2,water_path_of_time_id,start,&num_stepp1,(void *)water_path_of_time_2);
    wind_speed_id=ncvarid(cdfid_2,"wind_speed");
    ncvarget(cdfid_2,wind_speed_id,start,&num_layerp2,(void *)wind_speed_2);

    (void)fprintf(fp_err,"Reading second file's netCDF format cloud data from %s\n",in_file_2);
  } /* end if NETCDF_FORMAT */
  
  c_gopks(ier,bufa);
  c_gopwk(wkid,conid,wtype);
  c_gacwk(wkid); 
  
  c_cpgetr("VPL",&vpl);
  c_cpgetr("VPR",&vpr);
  c_cpgetr("VPB",&vpb);
  c_cpgetr("VPT",&vpt);
  if(debug == 33){
    (void)fprintf(fp_err,"vpl = %g, vpr = %g, vpb = %g, vpt = %g\n",
		  vpl,vpr,vpb,vpt);
  }

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,0.,0.,0.);

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* draw the Radiative Time Evolution plots using conpack routines */
    /* Plot the shortwave cloud forcing */
    for(step=0;step<=num_step-1;step++){
      abscissa[step]=time_array_1[step]*MINUTES_PER_SECOND; /* minutes */
      ordinate[step]=SW_cloud_forcing_1[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    abscissa_min=min_VEC(abscissa,num_step);
    abscissa_max=max_VEC(abscissa,num_step);
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Radiative Evolution of Cirrus Cloud");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step-1]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Time (min.)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Cloud Forcing (W/m:S:2:N:)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F5.2)","(F5.0)",0,0,10,10,0,15,0); 
    /* set the best domain for multiple plots of same dimensions */ 
    yaxis_min=ordinate_min;
    yaxis_max=ordinate_max;
    yaxis_min=(yaxis_min < (float_foo2=min_VEC(LW_cloud_forcing_1,num_step))) ?
      yaxis_min : float_foo2;
    yaxis_min=(yaxis_min < (float_foo2=min_VEC(LW_cloud_forcing_2,num_step))) ?
      yaxis_min : float_foo2;
    yaxis_min=(yaxis_min < (float_foo2=min_VEC(SW_cloud_forcing_2,num_step))) ?
      yaxis_min : float_foo2;
    yaxis_max=(yaxis_max > (float_foo2=max_VEC(LW_cloud_forcing_1,num_step))) ?
      yaxis_max : float_foo2;
    yaxis_max=(yaxis_max > (float_foo2=max_VEC(LW_cloud_forcing_2,num_step))) ?
      yaxis_max : float_foo2;
    yaxis_max=(yaxis_max > (float_foo2=max_VEC(SW_cloud_forcing_2,num_step))) ?
      yaxis_max : float_foo2;
    for(step=0;step<=num_step-1;step++){
      yaxis_min=(yaxis_min < 
		 (float_foo2=LW_cloud_forcing_1[step]+SW_cloud_forcing_1[step])) ?
		   yaxis_min : float_foo2;
      yaxis_max=(yaxis_max > float_foo2) ? yaxis_max : float_foo2;
      yaxis_min=(yaxis_min < 
		 (float_foo2=LW_cloud_forcing_2[step]+SW_cloud_forcing_2[step])) ?
		   yaxis_min : float_foo2;
      yaxis_max=(yaxis_max > float_foo2) ? yaxis_max : float_foo2;
    } /* end loop over steps */

    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,yaxis_min,yaxis_max,1);
    /* ensure the gridal routine accesses the highest quality characters */ 
/*    c_gaseti("LTY", 1);*/
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'$'$'$'$'$'$'$'$'$'$'$'$'$SWCF1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* print the header */
    float_foo=.5;
    (void)sprintf(line_label,"Ordinate [Minimum, Maximum, Ending] (units)");
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* actually draw in the abscissa = 0 line */
    c_line(abscissa_min,0.,abscissa_max,0.);
    
    /* print the normalization factor */
    (void)sprintf(line_label,"Shortwave Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;

    /* Overplot shortwave cloud forcing from the second file */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=SW_cloud_forcing_2[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* print the normalization factor */
    (void)sprintf(line_label,"Shortwave Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'$'$'$'$'$'$'$'$'$'$'$'$'$SWCF2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* Overplot longwave cloud forcing */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=LW_cloud_forcing_1[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* print the normalization factor */
    (void)sprintf(line_label,"Longwave Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$$'$$$'$$'$$$'$$LWCF1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* Overplot longwave cloud forcing from the second file*/
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=LW_cloud_forcing_2[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* print the normalization factor */
    (void)sprintf(line_label,"Longwave Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$$'$$$'$$'$$$'$$LWCF2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* Overplot net cloud forcing */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=LW_cloud_forcing_1[step]+SW_cloud_forcing_1[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* print the normalization factor */
    (void)sprintf(line_label,"Net Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$$$$NETCF1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* Overplot net cloud forcing from the second file */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=LW_cloud_forcing_2[step]+SW_cloud_forcing_2[step]; /* W/m^2 */ 
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* print the normalization factor */
    (void)sprintf(line_label,"Net Cloud Forcing [%.2f, %.2f, %.2f] W/m:S:2",
		  ordinate_min,ordinate_max,ordinate[num_step-1]);
    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),
	     line_label,-.7,0.,-1.);
    float_foo-=.05;
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$$$$NETCF2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* Overplot the visible optical depth */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=cloud_optical_depth_1[step];
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);
    /* set the best domain for multiple plots of same dimensions */ 
    yaxis_min=ordinate_min;
    yaxis_max=ordinate_max;
    yaxis_min=(yaxis_min < (float_foo2=min_VEC(cloud_optical_depth_2,num_step))) ?
      yaxis_min : float_foo2;
    yaxis_max=(yaxis_max > (float_foo2=max_VEC(cloud_optical_depth_2,num_step))) ?
      yaxis_max : float_foo2;

    /* normalize coordinates around the right-hand y-axis region */
    c_set(.9,1.,.1,.9,0.,1.,0.,1.,1);
    c_plchhq(.7,.5,"Visible Optical Depth :PGL1:S:R:",-.7,90.,0.);
    /* format the axis labels, putting the Y-axis labels on the right */
    c_labmod("(F6.2)","(F5.2)",0,0,10,10,1,0,0); 
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,0.,yaxis_max,1);
    /* arrange for the tick spacing and actually draw the axes, no x ticks */
    c_gridal(0,0,9,4,-1,1,5,0,0); 

    /* print the normalization factor */
/*    (void)sprintf(line_label,"Visible Optical Depth (.55 :PGL1:L:R:) :PGL1:S:R: [%.2f, %.2f, %.2f]",*/
/*		  ordinate_min,ordinate_max,ordinate[num_step-1]);*/
/*    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),*/
/*	     line_label,-.7,0.,-1.);*/
    float_foo-=.05;

    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$'$'$$$$$'$'$$$$$TAU1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);

    /* Overplot the visible optical depth from the second file */
    for(step=0;step<=num_step-1;step++){
      ordinate[step]=cloud_optical_depth_2[step];
    } /* end loop over time steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step);
    ordinate_min=min_VEC(ordinate,num_step);

    /* print the normalization factor */
/*    (void)sprintf(line_label,"Visible Optical Depth (.55 :PGL1:L:R:) :PGL1:S:R: [%.2f, %.2f, %.2f]",*/
/*		  ordinate_min,ordinate_max,ordinate[num_step-1]);*/
/*    c_plchhq(.1*abscissa_max,yaxis_min+float_foo*(yaxis_max-yaxis_min),*/
/*	     line_label,-.7,0.,-1.);*/
    float_foo-=.05;

    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$'$'$$$$$'$'$$$$$TAU2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa,ordinate,num_step);
    
    /* advance the frame */
    c_frame();
  } /* end if */ 
  
  if(True){
    /* draw the emissivity vs albedo scatterplot using conpack routines */
    num_rad_steps=num_step/rad_step;
    for(step=0;step<=num_rad_steps;step++){
      abscissa[step]=emissivity_of_time_1[step*rad_step];
      ordinate[step]=albedo_of_time_1[step*rad_step];
    } /* end loop over steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_rad_steps+1);
    ordinate_min=min_VEC(ordinate,num_rad_steps+1);
    abscissa_max=max_VEC(abscissa,num_rad_steps+1);
    abscissa_min=min_VEC(abscissa,num_rad_steps+1);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Emissivity vs. Albedo");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Emissivity :PGL1:E:R:",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Albedo A",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F5.2)","(F5.2)",0,0,10,10,0,15,0); 
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,
	  0.,1.,0.,.75,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,10,4,1,1,0,0,0); 
    /* plot the points */
    c_points(abscissa,ordinate,num_rad_steps+1,-2,0);
    
    /* Overplot the Emissivity vs. Albedo for the second file */ 
    for(step=0;step<=num_rad_steps;step++){
      abscissa[step]=emissivity_of_time_2[step*rad_step];
      ordinate[step]=albedo_of_time_2[step*rad_step];
    } /* end loop over steps */
    /* plot the points */
    c_points(abscissa,ordinate,num_rad_steps+1,-3,0);
    
    /* advance the frame */
    c_frame();
  } /* end if */
  
  if(True){
    /* draw the Effective Radius plot of both files using conpack routines */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=effective_radius_1[layer]*MICRONS_PER_METER; /* microns */
      ordinate[layer]=altitude_1[layer]*KILOMETERS_PER_METER; /* km */
      if(altitude_1[layer] < cloud_base_of_time_1[num_step] ||
	 altitude_1[layer] > cloud_top_of_time_1[num_step])
	abscissa[layer]=0.;
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    for(layer=1;layer<=num_layer;layer++){
      /*    abscissa[layer]/=abscissa_max;*/
    } /* end loop over layers */
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Effective Radius");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Effective Radius r:B3:eff    (:PGL1:L:R:)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F6.2)","(F5.1)",0,0,10,10,0,15,0); 
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,
	  abscissa_min,abscissa_max,ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
    /* print the header */
    float_foo=ordinate_min+.4*(ordinate_max-ordinate_min);
    (void)sprintf(line_label,"Abscissa [Minimum, Maximum] (units)");
    c_plchhq(abscissa_min+.2*(abscissa_max-abscissa_min),
	     float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* print the normalization factor */
    (void)sprintf(line_label,"Effective radius r:B:eff:N: [%.1f, %.1f] :PGL1:L:R:",
		  abscissa_min,abscissa_max);
    c_plchhq(abscissa_min+.2*(abscissa_max-abscissa_min),
	     float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$'$'$$$$$'$'$$$$$REFF1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the second file effective radius */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=effective_radius_2[layer]*MICRONS_PER_METER; /* microns */
      if(altitude_1[layer] < cloud_base_of_time_2[num_step] ||
	 altitude_1[layer] > cloud_top_of_time_2[num_step])
	abscissa[layer]=0.;
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Effective radius r:B:eff:N: [%.1f, %.1f] :PGL1:L:R:",
		  abscissa_min,abscissa_max);
    c_plchhq(abscissa_min+.2*(abscissa_max-abscissa_min),
	     float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$'$'$$$$$'$'$$$$$REFF2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* advance the frame */
    c_frame();
  } /* end if */
  
  if(True){
    /* draw the Effective Radius of the first files minus the second file */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=(effective_radius_1[layer]-effective_radius_2[layer])*
	MICRONS_PER_METER; /* microns */
      ordinate[layer]=altitude_1[layer]*KILOMETERS_PER_METER; /* km */
      if(altitude_1[layer] < cloud_base_of_time_1[num_step] ||
	 altitude_1[layer] > cloud_top_of_time_1[num_step])
	abscissa[layer]=0.;
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    for(layer=1;layer<=num_layer;layer++){
      /*    abscissa[layer]/=abscissa_max;*/
    } /* end loop over layers */
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Difference in Effective Radius");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Effective Radius r:B3:eff    (:PGL1:L:R:)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F6.2)","(F5.1)",0,0,10,10,0,15,0); 
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,
	  abscissa_min,abscissa_max,ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
    /* print the header */
    float_foo=ordinate_min+.4*(ordinate_max-ordinate_min);
    (void)sprintf(line_label,"Abscissa [Minimum, Maximum] (units)");
    c_plchhq(abscissa_min+.2*(abscissa_max-abscissa_min),
	     float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* print the normalization factor */
    (void)sprintf(line_label,"Effective radius r:B:eff:N: [%.1f, %.1f] :PGL1:L:R:",
		  abscissa_min,abscissa_max);
    c_plchhq(abscissa_min+.2*(abscissa_max-abscissa_min),
	     float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$'$'$$$$$'$'$$$$$REFF");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    /* advance the frame */
    c_frame();
  } /* end if */
  
  if(num_size == num_size_2){
    /* contour plot the final distribution of the first array minus the second */ 
    /* Mark the contours as a percent change from the first population */ 
    /* get the first two-dimensional array */ 
    start[0]=0L;
    start[1]=0L;
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid_1,"concentration");
    ncvarget(cdfid_1,concentration_id,start,count,(void *)concentration_ptr);
    
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert #/m^3/m --> #/m^3/micron and transpose, FORTRAN indexing */ 
	layer_size[layer-1][size-1]=concentration[layer][size]/
	  (delta_length_1[size]*MICRONS_PER_METER);
      } /* end loop over sizes */
    } /* end loop over layers */
    
    /* get the second two-dimensional array */ 
    start[0]=0L;
    start[1]=0L;
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid_2,"concentration");
    ncvarget(cdfid_2,concentration_id,start,count,(void *)concentration_ptr);
    
    /* subtract the second from the first */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert #/m^3/m --> #/m^3/micron and transpose, FORTRAN indexing */ 
	float_foo=concentration[layer][size]/
	  (delta_length_2[size]*MICRONS_PER_METER);
	if(layer_size[layer-1][size-1] > .03125){
	  layer_size[layer-1][size-1]=
	    100.*(layer_size[layer-1][size-1]-float_foo)/
	      layer_size[layer-1][size-1];
	}else{ 
	  layer_size[layer-1][size-1]=0.;
	} /* end else */
      } /* end loop over sizes */
    } /* end loop over layers */
    
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length_1[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude_1[layer]/1000.; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"%% Difference in Size Distributions");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]/60.); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Crystal Length (microns)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F5.0)","(F5.1)",0,0,10,10,0,15,0); 
    /* only plot the interesting regions */ 
    abscissa_min=0.;
    abscissa_max=600.;
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    /* actually draw in the abscissa = reff line */
    cloud_effective_length=equiv_rad_to_bullet_length(cloud_effective_radius_1[num_step]);
    if(debug == 16){
      (void)fprintf(fp_err,"cloud_effective_radius = %f, cloud_effective_length = %f\n",
		    cloud_effective_radius_1[num_step]*MICRONS_PER_METER,
		    cloud_effective_length*MICRONS_PER_METER);
    } /* end debug */
    c_line(cloud_effective_length*MICRONS_PER_METER,ordinate_min,
	   cloud_effective_length*MICRONS_PER_METER,ordinate_max);
    
    /* Overplot the effective length in the layer */
    for(layer=1;layer<=num_layer;layer++){
      effective_length_1[layer]=equiv_rad_to_bullet_length(effective_radius_1[layer])*
	MICRONS_PER_METER; /* microns */
      if(effective_length_1[layer] > abscissa_max){
	;
      } /* end if */
    } /* end loop over layers */
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(effective_length_1+1,ordinate+1,num_layer);
    
    xdim=num_size;
    ydim=num_layer;
    /* turn on the clipping indicator */
    c_gsclip(1);
    /* force plotchar to use characters of the highest quality */
    c_pcseti("QU - QUALITY FLAG",0);
    /* ensure that conpack doesn't call set itself */ 
    c_cpseti("SET",0);
    /* ensure that calls to cppkcl don't do the level picking */ 
    c_cpseti("CLS -- contour level selection",0);
    /* set the type of coordinate transforation used in cpmpxy */ 
    c_cpseti("MAP",4); 
    /* associate X and Y coordinates with their physical values */ 
    c_cpsetr("XC1",1);
    c_cpsetr("XCM",num_size);
    c_cpsetr("YC1",ordinate_min);
    c_cpsetr("YCN",ordinate_max);
    /* world coordinates around the plot region for the non-linear mapping */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=41;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      float_foo2=-40.+(level-1.)*2.;
/*      float_foo2=-100.+(level-1.)*5.;*/
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f #/m^3/u\n",level,float_foo2);
      } /* end debug */
      if(float_foo2 >= 1.){
	(void)sprintf(line_label,"%i",ROUNDER(float_foo2));
      }else if(float_foo2 >= .5){
	(void)sprintf(line_label,"%.1f",float_foo2);
      }else if(float_foo2 >= .1){
	(void)sprintf(line_label,"%.2f",float_foo2);
      }else{
	(void)sprintf(line_label,"%.3f",float_foo2);
      } /* end else */ 
      c_cpseti("PAI -- parameter array index",level);
      c_cpsetr("CLV -- contour level value",float_foo2);
      c_cpseti("CLU -- contour level use",3);
      c_cpsetr("CLL -- contour level line width",1.);
      c_cpseti("AIA -- area identifier above line",level+1);
      c_cpseti("AIB -- area identifier below line",level);
      if(float_foo2 < 0.)
	c_cpsetc("CLD -- contour line dash pattern","$'$'$'$'$'$'$'$'");
      else
	c_cpsetc("CLD -- contour line dash pattern","$$$$$$$$$$$$$$$$");
      c_cpsetc("LLT -- contour label text",line_label);
      c_cpseti("CLC -- contour line color",-1);
      c_cpseti("LLC -- contour line label color",-1);
    } /* end loop over levels */
    /* initialize the drawing of the contour plot */
    c_cprect(layer_size_ptr,
	     xdim,xdim,ydim,
	     rwrk,krwk,iwrk,kiwk);
    /* draw contour lines and labels */
    c_cpcldr(layer_size_ptr,rwrk,iwrk);
    /* add the informational label and the high/low labels */
    (void)sprintf(line_label,"Values in %% from $ZMN$ to $ZMX$ by steps of 2%%");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */
  
  if(True){
    /* contour plot the final distribution of the first array minus the second */ 
    /* get the first two-dimensional array */ 
    start[0]=0L;
    start[1]=0L;
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid_1,"concentration");
    ncvarget(cdfid_1,concentration_id,start,count,(void *)concentration_ptr);
    
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert #/m^3/m --> #/m^3/micron and transpose, FORTRAN indexing */ 
	layer_size[layer-1][size-1]=concentration[layer][size]/
	  (delta_length_1[size]*MICRONS_PER_METER);
      } /* end loop over sizes */
    } /* end loop over layers */
    
    /* get the second two-dimensional array */ 
    start[0]=0L;
    start[1]=0L;
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid_2,"concentration");
    ncvarget(cdfid_2,concentration_id,start,count,(void *)concentration_ptr);
    
    /* subtract the second from the first */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert #/m^3/m --> #/m^3/micron and transpose, FORTRAN indexing */ 
	layer_size[layer-1][size-1]-=concentration[layer][size]/
	  (delta_length_2[size]*MICRONS_PER_METER);
      } /* end loop over sizes */
    } /* end loop over layers */
    
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length_1[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude_1[layer]/1000.; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Difference in Size Distributions");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]/60.); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Crystal Length (microns)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F5.0)","(F5.1)",0,0,10,10,0,15,0); 
    /* only plot the interesting regions */ 
    abscissa_min=0.;
    abscissa_max=600.;
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    /* actually draw in the abscissa = reff line */
    cloud_effective_length=equiv_rad_to_bullet_length(cloud_effective_radius_1[num_step]);
    if(debug == 16){
      (void)fprintf(fp_err,"cloud_effective_radius = %f, cloud_effective_length = %f\n",
		    cloud_effective_radius_1[num_step]*MICRONS_PER_METER,
		    cloud_effective_length*MICRONS_PER_METER);
    } /* end debug */
    c_line(cloud_effective_length*MICRONS_PER_METER,ordinate_min,
	   cloud_effective_length*MICRONS_PER_METER,ordinate_max);
    
    /* Overplot the effective length in the layer */
    for(layer=1;layer<=num_layer;layer++){
      effective_length_1[layer]=equiv_rad_to_bullet_length(effective_radius_1[layer])*
	MICRONS_PER_METER; /* microns */
      if(effective_length_1[layer] > abscissa_max){
	;
      } /* end if */
    } /* end loop over layers */
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(effective_length_1+1,ordinate+1,num_layer);
    
    xdim=num_size;
    ydim=num_layer;
    /* turn on the clipping indicator */
    c_gsclip(1);
    /* force plotchar to use characters of the highest quality */
    c_pcseti("QU - QUALITY FLAG",0);
    /* ensure that conpack doesn't call set itself */ 
    c_cpseti("SET",0);
    /* ensure that calls to cppkcl do nothing, so i can pick 'em myself */ 
    c_cpseti("CLS -- contour level selection",0);
    /* set the type of coordinate transforation used in cpmpxy */ 
    c_cpseti("MAP",4); 
    /* associate X and Y coordinates with their physical values */ 
    c_cpsetr("XC1",1);
    c_cpsetr("XCM",num_size);
    c_cpsetr("YC1",ordinate_min);
    c_cpsetr("YCN",ordinate_max);
    /* world coordinates around the plot region for the non-linear mapping */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=NUM_DIST_DIFF_CONTOURS;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      float_foo2=dist_diff_contour_level[level];
      /*    float_foo2=pow(2.,(float)(level-6));*/
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f #/m^3/u\n",level,float_foo2);
      } /* end debug */
      if(float_foo2 >= 1.){
	(void)sprintf(line_label,"%i",ROUNDER(float_foo2));
      }else if(float_foo2 >= .5){
	(void)sprintf(line_label,"%.1f",float_foo2);
      }else if(float_foo2 >= .1){
	(void)sprintf(line_label,"%.2f",float_foo2);
      }else{
	(void)sprintf(line_label,"%.3f",float_foo2);
      } /* end else */ 
      c_cpseti("PAI -- parameter array index",level);
      c_cpsetr("CLV -- contour level value",float_foo2);
      if(float_foo2 != 0.)
	c_cpseti("CLU -- contour level use",3);
      else
	c_cpseti("CLU -- contour level use",0);
      c_cpsetr("CLL -- contour level line width",1.);
      c_cpseti("AIA -- area identifier above line",level+1);
      c_cpseti("AIB -- area identifier below line",level);
      if(float_foo2 < 0.)
	c_cpsetc("CLD -- contour line dash pattern","$'$'$'$'$'$'$'$'");
      else
	c_cpsetc("CLD -- contour line dash pattern","$$$$$$$$$$$$$$$$");
      c_cpsetc("LLT -- contour label text",line_label);
      c_cpseti("CLC -- contour line color",-1);
      c_cpseti("LLC -- contour line label color",-1);
    } /* end loop over levels */
    /* initialize the drawing of the contour plot */
    c_cprect(layer_size_ptr,
	     xdim,xdim,ydim,
	     rwrk,krwk,iwrk,kiwk);
    /* draw contour lines and labels */
    c_cpcldr(layer_size_ptr,rwrk,iwrk);
    /* add the informational label and the high/low labels */
    (void)sprintf(line_label,"Values (in #/m:S:3:N:/:PGL1:L:R:) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */
  
  if(debug == 12){
    /* contour plot the final distribution with colors */ 
    red_blue_scale(NUM_DIST_DIFF_CONTOURS);
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length_1[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude_1[layer]/1000.; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Difference in Size Distributions (#-m:S:-3:N:-:PGL1:L:R:m:S:-1:N:) Columns - Spheres");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]/60.); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
/*    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); */
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Crystal Length (:PGL1:L:R:m)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* only plot the interesting regions */ 
    abscissa_min=3.;
    abscissa_max=300.;
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    
    xdim=num_size;
    ydim=num_layer;
    /* turn on the clipping indicator */
    c_gsclip(1);
    /* force plotchar to use characters of the highest quality */
    c_pcseti("QU - QUALITY FLAG",0);
    /* ensure that conpack doesn't call set itself */ 
    c_cpseti("SET",0);
    /* ensure that calls to cppkcl do nothing, so i can pick 'em myself */ 
    c_cpseti("CLS -- contour level selection",0);
    /* set the type of coordinate transforation used in cpmpxy */ 
    c_cpseti("MAP",4); 
    /* associate X and Y coordinates with their physical values */ 
    c_cpsetr("XC1",1);
    c_cpsetr("XCM",num_size);
    c_cpsetr("YC1",ordinate_min);
    c_cpsetr("YCN",ordinate_max);
    /* world coordinates around the plot region for the non-linear mapping */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=NUM_DIST_DIFF_CONTOURS;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      float_foo2=dist_diff_contour_level[level];
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f #/m^3/u\n",level,float_foo2);
      } /* end debug */
      if(float_foo2 >= 1.){
	(void)sprintf(line_label,"%i",ROUNDER(float_foo2));
      }else if(float_foo2 >= .5){
	(void)sprintf(line_label,"%.1f",float_foo2);
      }else if(float_foo2 >= .1){
	(void)sprintf(line_label,"%.2f",float_foo2);
      }else{
	(void)sprintf(line_label,"%.3f",float_foo2);
      } /* end else */ 
      c_cpseti("PAI -- parameter array index",level);
      c_cpsetr("CLV -- contour level value",float_foo2);
      if(float_foo2 != 0.)
	c_cpseti("CLU -- contour level use",1);
      else
	c_cpseti("CLU -- contour level use",0);
      c_cpsetr("CLL -- contour level line width",1.);
      c_cpseti("AIA -- area identifier above line",level);
      c_cpseti("AIB -- area identifier below line",level-1);
      if(float_foo2 < 0.)
	c_cpsetc("CLD -- contour line dash pattern","$'$'$'$'$'$'$'$'");
      else
	c_cpsetc("CLD -- contour line dash pattern","$$$$$$$$$$$$$$$$");
      c_cpsetc("LLT -- contour label text",line_label);
      c_cpseti("CLC -- contour line color",-1);
      c_cpseti("LLC -- contour line label color",-1);
    } /* end loop over levels */
    /* set all the GKS aspect source flags to "individual" */ 
    c_gsasf(iasf);
    /* initialize the software fill package to do the desired type of fill */ 
    c_sfseti("TYPE OF FILL",0);
    /* force solid fill */ 
    c_gsfais(1);
    /* initialize the drawing of the contour plot */
    c_cprect(layer_size_ptr,
	     xdim,xdim,ydim,
	     rwrk,krwk,iwrk,kiwk);
    /* initialize area map */ 
    c_arinam(iam,lam);
    /* add the contour lines to the area map */ 
    c_cpclam(layer_size_ptr,rwrk,iwrk,iam);
    /* color in the area map */ 
    c_arscam(iam,xra,yra,nra,
	     iai,iag,mai,dist_colram_);
    /* thicken the lines of the label bar */ 
    c_lbsetr("WBL -- width of box lines",2.);
    /* Put a Vertical label bar to the right of the contour plot */ 
    c_lblbar(
	     1, /* IHOV */ 
	     .91, /* XLEB */ 
	     1., /* XREB */ 
	     .1, /* YBEB */ 
	     .9, /* YTEB */ 
	     num_contour_level, /* NBOX */ 
	     .5, /* WSFB */ 
	     1., /* HSFB */ 
	     dist_diff_grey_lind, /* LFIN */ 
	     0, /* IFTP */ 
	     dist_diff_llbs, /* LLBS */ 
	     num_contour_level+1, /* NLBS */ 
	     1 /* LBAB */ 
	     );
    /* draw contour lines and labels */
    c_cpcldr(layer_size_ptr,rwrk,iwrk);
    /* format the informational label and the high/low labels */
    (void)sprintf(line_label,"Values (in #-m:S:-3:N:-:PGL1:L:R:m:S:-1:N:) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetc("HIT -- high label text string","");
    c_cpsetc("LOT -- low label text string","");
    c_cpsetr("ILY -- y coord of info label",.05);
    /* Draw Informational, High, and Low labels. */
/*    c_cplbdr(layer_size_ptr,rwrk,iwrk);*/
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* format the axis labels */
    c_labmod("(F5.0)","(F5.1)",0,0,10,10,0,15,0); 
    /* set the tick-mark length and direction */ 
    c_tick4(-12,-8,-12,-8);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    /* actually draw in the abscissa = reff line */
/*    cloud_effective_length=equiv_rad_to_bullet_length(cloud_effective_radius_1[num_step]);*/
    if(debug == 16){
      (void)fprintf(fp_err,"cloud_effective_radius = %f, cloud_effective_length = %f\n",
		    cloud_effective_radius_1[num_step]*MICRONS_PER_METER,
		    cloud_effective_length*MICRONS_PER_METER);
    } /* end debug */
/*    c_line(cloud_effective_length*MICRONS_PER_METER,ordinate_min,*/
/*	   cloud_effective_length*MICRONS_PER_METER,ordinate_max);*/
    
    /* Overplot the effective length in the layer */
    for(layer=1;layer<=num_layer;layer++){
      effective_length_1[layer]=equiv_rad_to_bullet_length(effective_radius_1[layer])*
	MICRONS_PER_METER; /* microns */
      if(effective_length_1[layer] > abscissa_max){
	;
      } /* end if */
    } /* end loop over layers */
    /* set the dash characteristics of the line */
/*    (void)sprintf(line_label,"$'");*/
/*    c_dashdc(line_label,jcrt,jsize);*/
    /* draw the curve */
/*    c_curved(effective_length_1+1,ordinate+1,num_layer);*/
    /* advance the frame */
    c_frame();
  } /* end debug */
  
  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* draw the shortwave, longwave, and net heating rates from both files */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_net_1[layer]*SECONDS_PER_HOUR; /* K/hour */
      ordinate[layer]=altitude_1[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Heating Rates");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Heating Rate (K/hour)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F6.2)","(F5.1)",0,0,10,10,0,15,0); 
    /* set the best domain for multiple plots of same dimensions */ 
    xaxis_min=abscissa_min;
    xaxis_max=abscissa_max;
    xaxis_min=(xaxis_min < 
	       (float_foo2=
		SECONDS_PER_HOUR*min_VEC(heating_rate_SW_1+1,num_layer))) ? 
		  abscissa_min : float_foo2;
    xaxis_min=(xaxis_min <
	       (float_foo2=
		SECONDS_PER_HOUR*min_VEC(heating_rate_LW_1+1,num_layer))) ? 
		  xaxis_min : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=
		SECONDS_PER_HOUR*max_VEC(heating_rate_SW_1+1,num_layer))) ? 
		  abscissa_max : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=
		SECONDS_PER_HOUR*max_VEC(heating_rate_LW_1+1,num_layer))) ? 
		  xaxis_max : float_foo2;


    xaxis_min=(xaxis_min < 
	       (float_foo2=
		SECONDS_PER_HOUR*min_VEC(heating_rate_SW_2+1,num_layer))) ? 
		  abscissa_min : float_foo2;
    xaxis_min=(xaxis_min < 
	       (float_foo2=
		SECONDS_PER_HOUR*min_VEC(heating_rate_SW_2+1,num_layer))) ? 
		  abscissa_min : float_foo2;
    xaxis_min=(xaxis_min <
	       (float_foo2=
		SECONDS_PER_HOUR*min_VEC(heating_rate_LW_2+1,num_layer))) ? 
		  xaxis_min : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=
		SECONDS_PER_HOUR*max_VEC(heating_rate_SW_2+1,num_layer))) ? 
		  abscissa_max : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=
		SECONDS_PER_HOUR*max_VEC(heating_rate_SW_2+1,num_layer))) ? 
		  abscissa_max : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=
		SECONDS_PER_HOUR*max_VEC(heating_rate_LW_2+1,num_layer))) ? 
		  xaxis_max : float_foo2;

    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,xaxis_min,xaxis_max,
	  ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
    /* print the header */
    float_foo=ordinate_min+.3*(ordinate_max-ordinate_min);
    (void)sprintf(line_label,"Abscissa [Minimum, Maximum] (units)");
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* actually draw in the abscissa = 0 line */
    c_line(0.,ordinate_min,0.,ordinate_max);
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$Net1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    (void)sprintf(line_label,"Net Radiative Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* Overplot the second file net heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_net_2[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Net Radiative Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$Net2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the first file shortwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_SW_1[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Shortwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'$'$'$'$'$'$'$'$'$'$'SW1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the second file shortwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_SW_2[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Shortwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'$'$'$'$'$'$'$'$'$'$'SW2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the first file longwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_LW_1[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Longwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$LW1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the second file longwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_LW_2[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Longwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$LW2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the first file latent heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=latent_heating_1[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Latent Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$'$$$'$'$$$'$'$$$'$'$$$'$'$$$LAT1");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the second file latent heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=latent_heating_2[layer]*SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Latent Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$'$$$'$'$$$'$'$$$'$'$$$'$'$$$LAT2");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    if(debug == 10){
      /* Overplot the first file advective heating rate */
      for(layer=1;layer<=num_layer;layer++){
	abscissa[layer]=advective_heating_1[layer]*SECONDS_PER_HOUR; /* K/hour */
      } /* end loop over layers */
      /* Normalize the abscissas */
      abscissa_max=max_VEC(abscissa+1,num_layer);
      abscissa_min=min_VEC(abscissa+1,num_layer);
      /* print the normalization factor */
      (void)sprintf(line_label,"Advective Heating Rate [%.1f, %.1f] K/hour",
		    abscissa_min,abscissa_max);
      c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
      float_foo-=.04*(ordinate_max-ordinate_min);
      /* set the dash characteristics of the line */
      (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$ADV1");
      c_dashdc(line_label,jcrt,jsize);
      /* draw the curve */
      c_curved(abscissa+1,ordinate+1,num_layer);

      /* Overplot the second file advective heating rate */
      for(layer=1;layer<=num_layer;layer++){
	abscissa[layer]=advective_heating_2[layer]*SECONDS_PER_HOUR; /* K/hour */
      } /* end loop over layers */
      /* Normalize the abscissas */
      abscissa_max=max_VEC(abscissa+1,num_layer);
      abscissa_min=min_VEC(abscissa+1,num_layer);
      /* print the normalization factor */
      (void)sprintf(line_label,"Advective Heating Rate [%.1f, %.1f] K/hour",
		    abscissa_min,abscissa_max);
      c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
      float_foo-=.04*(ordinate_max-ordinate_min);
      /* set the dash characteristics of the line */
      (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$ADV2");
      c_dashdc(line_label,jcrt,jsize);
      /* draw the curve */
      c_curved(abscissa+1,ordinate+1,num_layer);
    } /* end if */
    
    /* advance the frame */
    c_frame();
  } /* end if */ 

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* draw the differences in the shortwave, longwave, and net heating rates */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=(heating_rate_net_1[layer]-heating_rate_net_2[layer])*SECONDS_PER_HOUR; /* K/hour */
      ordinate[layer]=altitude_1[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Difference in Heating Rates");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array_1[num_step]*MINUTES_PER_SECOND); 
    c_plchhq(.5,.2,line_label,-.7,0.,0.); 
    /* plotchar the command line */
    c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
    /* normalize coordinates around the x-axis region */
    c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
    c_plchhq(.5,.5,"Heating Rate (K/hour)",-.7,0.,0.);
    /* normalize coordinates around the y-axis region */
    c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
    c_plchhq(.3,.5,"Altitude z (km)",-.7,90.,0.);
    /* format the axis labels */
    c_labmod("(F6.2)","(F5.1)",0,0,10,10,0,15,0); 
    /* set the best domain for multiple plots of same dimensions */ 
    xaxis_min=abscissa_min;
    xaxis_max=abscissa_max;
    for(layer=1;layer<=num_layer;layer++){
    xaxis_min=(xaxis_min < 
	       (float_foo2=(heating_rate_SW_1[layer]-heating_rate_SW_2[layer])*
		SECONDS_PER_HOUR)) ? xaxis_min : float_foo2;
    xaxis_min=(xaxis_min < 
	       (float_foo2=(heating_rate_LW_1[layer]-heating_rate_LW_2[layer])*
		SECONDS_PER_HOUR)) ? xaxis_min : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=(heating_rate_SW_1[layer]-heating_rate_SW_2[layer])*
		SECONDS_PER_HOUR)) ? xaxis_max : float_foo2;
    xaxis_max=(xaxis_max > 
	       (float_foo2=(heating_rate_LW_1[layer]-heating_rate_LW_2[layer])*
		SECONDS_PER_HOUR)) ? xaxis_max : float_foo2;
    } /* end loop over layers */
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,xaxis_min,xaxis_max,
	  ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
    /* print the header */
    float_foo=ordinate_min+.3*(ordinate_max-ordinate_min);
    (void)sprintf(line_label,"Abscissa [Minimum, Maximum] (units)");
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* actually draw in the abscissa = 0 line */
    c_line(0.,ordinate_min,0.,ordinate_max);
    
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$Net");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    (void)sprintf(line_label,"Net Radiative Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* Overplot the shortwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=(heating_rate_SW_1[layer]-heating_rate_SW_2[layer])*
	SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Shortwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$'$'$'$'$'$'$'$'$'$'$'SW");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the longwave heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=(heating_rate_LW_1[layer]-heating_rate_LW_2[layer])*
	SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Longwave Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$LW");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    /* Overplot the latent heating rate */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=(latent_heating_1[layer]-latent_heating_2[layer])*
	SECONDS_PER_HOUR; /* K/hour */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    /* print the normalization factor */
    (void)sprintf(line_label,"Latent Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* set the dash characteristics of the line */
    (void)sprintf(line_label,"$$$'$'$$$'$'$$$'$'$$$'$'$$$'$'$$$LAT");
    c_dashdc(line_label,jcrt,jsize);
    /* draw the curve */
    c_curved(abscissa+1,ordinate+1,num_layer);
    
    if(debug == 10){
      /* Overplot the advective heating rate */
      for(layer=1;layer<=num_layer;layer++){
	abscissa[layer]=(advective_heating_1[layer]-advective_heating_2[layer])*
	  SECONDS_PER_HOUR; /* K/hour */
      } /* end loop over layers */
      /* Normalize the abscissas */
      abscissa_max=max_VEC(abscissa+1,num_layer);
      abscissa_min=min_VEC(abscissa+1,num_layer);
      /* print the normalization factor */
      (void)sprintf(line_label,"Advective Heating Rate [%.1f, %.1f] K/hour",
		    abscissa_min,abscissa_max);
      c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
      float_foo-=.04*(ordinate_max-ordinate_min);
      /* set the dash characteristics of the line */
      (void)sprintf(line_label,"$$$'$$$'$$$'$$$'$$$'$$$ADV");
      c_dashdc(line_label,jcrt,jsize);
      /* draw the curve */
      c_curved(abscissa+1,ordinate+1,num_layer);
    } /* end if */
    
    /* advance the frame */
    c_frame();
  } /* end if */ 

  if(NETCDF_FORMAT){
    ncclose(cdfid_1);
    ncclose(cdfid_2);
  } /* end if */

  c_gdawk(wkid);
  c_gclwk(wkid);
  c_gclks(); 

  Exit_gracefully();
  exit(0);
} /* end of main() */
 
#ifdef CRAY
void CPMPXY(imap,xinp,yinp,xotp,yotp)
#else
void cpmpxy_(imap,xinp,yinp,xotp,yotp)
#endif
     int *imap;
     float *xinp,*yinp,*xotp,*yotp;
{
  /* Recall that this routine is called internally by FORTRAN programs,
     so it's kind of amazing that it works at all, clearly all variables
     are passed by reference as all I/O is to FORTRAN */
  
#ifdef CRAY
  void MAPTRN();
#else
  void maptrn_();
#endif

  int size;

  if(*imap == 1){
#ifdef CRAY
  MAPTRN(yinp,xinp,xotp,yotp);
#else
  maptrn_(yinp,xinp,xotp,yotp);
#endif
  }else if(*imap == 2){
    *xotp=*xinp*cos(.017453292519943**yinp);
    *yotp=*xinp*cos(.017453292519943**yinp);
  }else if(*imap == 3){
    *xotp=*xinp;
    *yotp=*yinp;
  }else if(*imap == 4){
    /* REMEMBER to fix this to be exact some day !!!  i.e. stop sending
     (xdim+2)*(ydim+2) arrays with meaningless edges */ 
    size=(int)(*xinp);
    *xotp=abscissa[size]+(*xinp-(float)size)*
      (abscissa[size+1]-abscissa[size]);
/*    *xotp=((float)size+1.-*xinp)*abscissa[size-1]+*/
/*      (*xinp-(float)size)*abscissa[size];*/
    *yotp=*yinp;
  }else{
    *xotp=*xinp;
    *yotp=*yinp;
  } /* end else */
  if(debug == 58){
    (void)fprintf(fp_err,"*imap = %i, *xinp = %g, *yinp = %g, *xotp = %g, *yotp = %g\n",*imap,*xinp,*yinp,*xotp,*yotp); 
  } /* end debug */
  
} /* end cpmpxy_() or CPMPXY() */ 

int Cray_rounder(float_val)
     float float_val;
{
  int nearest_int;

  nearest_int = ( float_val >= 0.) ? floor(float_val+.5) : floor(float_val-.5);

  return nearest_int;
} /* end Cray_rounder() */ 

void Exit_gracefully(void)
{
  char *time_buf_finish;
  time_t clock;

  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(fp_err,"\tfinish = %s\n",time_buf_finish);

  (void)fclose(fp_err);
  (void)fclose(fp_in);
  (void)fclose(fp_out);
}

void print_usage()
{
  (void)fprintf(stderr,
		"\nusage: iographs [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"D:d:Ee:Ii:j:Oo:vV\n\n");
  (void)fprintf(stderr,"-D debug Set the debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is False\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"-d debug_value Default is 0\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-i in_file_1 get the 1st file name. Default is cloud.nc\n");
  (void)fprintf(stderr,"-j in_file_2 get the 2nd file name. Default is cloud.nc\n");
  (void)fprintf(stderr,"-o outfile get the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"\n");
}

int shadam_(xra,yra,nra,iai,iag,nai)
     float *xra,*yra;
     int *nra,*iai,*iag,*nai;
{
  /* Internal routine supplied to arscam() */ 

  /* This version of SHADAM shades the area-map polygon whose edge is 
   defined by the points ((XRA(I),YRA(I),I=1,NRA) if and only, relative
   to edge group 3, its area identifier is between 1 and 20.  The
   density of the shading increases with the area identifier. */ 
  
  /* Declare and allocate the workspaces required by SFSGFA */ 
  float *dst;

  int *ind;

  int area_identifier;
  int num_area_identifiers;
  int dst_size;
  int ind_size;
  int ish;

  dst_size = 2**nra;
  ind_size = 3**nra;

  if(
     ((ind=(int *)malloc((ind_size)*sizeof(int))) == NULL ) ||
     ((dst=(float *)malloc((dst_size)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in shadam_\n");
    exit(1);
  } /* end if */

  /* Turn off shading */ 
  ish=0;

  /* If the area identifier for group 3 is in the right range, turn on
   shading. */ 
  num_area_identifiers=*nai;
  for(area_identifier=0;area_identifier<=num_area_identifiers-1;area_identifier++){
    if(iag[area_identifier] == 3 && 
       iai[area_identifier] >= .1 && 
       iai[area_identifier] <= 200.){
      ish=iai[area_identifier];
    } /* end if */
  } /* end loop over area_identifiers */

  /* If shading is turned on, shade the area.  The last point of the 
   edge is redundant and may be omitted. */ 
  if(ish != 0){
    c_sfsgfa(
	     xra,
	     yra,
	     *nra,
	     dst,
	     dst_size,
	     ind,
	     ind_size,
	     ish-1
	     );
  } /* end if */

  /* free the unneeded memory */ 
/*  free(ind);*/
/*  free(dst);*/

} /* end shadam_() */ 

int IWC_colram_(xra,yra,nra,iai,iag,nai)
     float *xra,*yra;
     int *nra,*iai,*iag,*nai;
{
  /* Internal routine supplied to arscam() */ 

  /* The arrays XRA and YRA, for indices 1 to NRA, contain the X and Y
   coordinates of points defining a polygon.  The area identifiers in
   the array IAI, each with an associated group indentifier in the array
   IAG, tell us whether the polygon is to be color-filled or not. */ 

  int area_identifier;
  int num_area_identifiers;
  int ifll;

  /* assume the polygon will be filled until we find otherwise */ 
  ifll=1;

  /* if any of the area identifiers is negative, don't fill the polygon */ 
  num_area_identifiers=*nai;
  for(area_identifier=0;area_identifier<=num_area_identifiers-1;area_identifier++){
    if(iai[area_identifier] < 0){
      ifll=0;
    } /* end if */
  } /* end loop over area_identifiers */

  /* Otherwise, fill the polygon in the color implied by its area
     identifier relative to edge group 3 (the contour-line group). */ 
  if(ifll != 0){
    ifll=0;
    for(area_identifier=0;area_identifier<=num_area_identifiers-1;area_identifier++){
      if(iag[area_identifier] == 3){
	 ifll=iai[area_identifier];
      } /* end if */
    }/* end loop over area_identifiers */
    if(ifll > 0 && ifll < NUM_IWC_CONTOURS){
      c_gsfaci(ifll+1);
      c_gfa(*nra-1,xra,yra);
    } /* end if */
  } /* end if */

} /* end IWC_colram_() */ 

int dist_colram_(xra,yra,nra,iai,iag,nai)
     float *xra,*yra;
     int *nra,*iai,*iag,*nai;
{
  /* Internal routine supplied to arscam() */ 

  /* The arrays XRA and YRA, for indices 1 to NRA, contain the X and Y
   coordinates of points defining a polygon.  The area identifiers in
   the array IAI, each with an associated group indentifier in the array
   IAG, tell us whether the polygon is to be color-filled or not. */ 

  int area_identifier;
  int num_area_identifiers;
  int ifll;

  /* assume the polygon will be filled until we find otherwise */ 
  ifll=1;

  /* if any of the area identifiers is negative, don't fill the polygon */ 
  num_area_identifiers=*nai;
  for(area_identifier=0;area_identifier<=num_area_identifiers-1;area_identifier++){
    if(iai[area_identifier] < 0){
      ifll=0;
    } /* end if */
  } /* end loop over area_identifiers */

  /* Otherwise, fill the polygon in the color implied by its area
     identifier relative to edge group 3 (the contour-line group). */ 
  if(ifll != 0){
    ifll=0;
    for(area_identifier=0;area_identifier<=num_area_identifiers-1;area_identifier++){
      if(iag[area_identifier] == 3){
	 ifll=iai[area_identifier];
      } /* end if */
    }/* end loop over area_identifiers */
    if(ifll > 0 && ifll < NUM_DIST_DIFF_CONTOURS){
      c_gsfaci(ifll+1);
      c_gfa(*nra-1,xra,yra);
    } /* end if */
  } /* end if */

} /* end dist_colram_() */ 

int greyscale(num_colors)
     int num_colors;
{
/* define a set of grey scale color triples for colors 2 through 22 */ 

  int color;

  float red;
  float green;
  float blue;

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,0.,0.,0.);
  
  for(color=1;color<=num_colors;color++){
    red=green=blue=1.-((float)color-1.)/((float)num_colors-1.);
    c_gscr(1,color+1,red,green,blue);
  } /* end loop over colors */

} /* end greyscale() */ 

int colorscale(num_colors)
     int num_colors;
{
/* define a set of grey scale color triples for colors 2 through 25 */ 

  int color;

  float red;
  float green;
  float blue;
  float cmin;
  float cmax;
  float hue;
  float hue1;

  cmin=-20.;
  cmax=250.;

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,0.,0.,0.);
  /* this makes color 2 white */ 
  c_gscr(1,2,1.,1.,1.);
  
  for(color=1;color<=num_colors;color++){
    hue1 =(color-1.)/(num_colors-1.);
    hue = hue1 + .05*sin(4*3.1415*hue1);
    hue = cmin + hue*(cmax - cmin);
    c_hsvrgb(hue,1.,1.,&red,&green,&blue);
    c_gscr(1,color+2,red,green,blue);
  } /* end loop over colors */

} /* end colorscale() */ 

int red_blue_scale(num_colors)
     int num_colors;
{
  /* define a set of color triples for colors 2 through num_colors + 1 */ 

  int color;

  float red;
  float green;
  float blue;
  float hue;
  float saturation;
  float value;

  value=1.;

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,0.,0.,0.);
  /* this makes color 2 white */ 
/*  c_gscr(1,2,1.,1.,1.);*/
  
  /* make the first half of the colors in a reddish scale, and
     make the second half in a bluish scale */ 
  hue=215.;
  for(color=1;color<=num_colors/2;color++){
    saturation=1.-(color-1)/((num_colors/2.)-1);
    c_hsvrgb(hue,saturation,value,&red,&green,&blue);
    c_gscr(1,color+1,red,green,blue);
    if(debug == 12){
/*      (void)fprintf(fp_err,"color+1=%i,r=%f,g=%f,b=%f,hue=%f,sat=%f,val=%f\n",*/
/*		    color+1,red,green,blue,hue,saturation,value);*/
    } /* end debug */
  } /* end loop over colors */

  hue=0.;
  for(color=1+num_colors/2;color<=num_colors+1;color++){
    saturation=(color-1.-num_colors/2.)/(num_colors/2.);
    c_hsvrgb(hue,saturation,value,&red,&green,&blue);
    c_gscr(1,color+1,red,green,blue);
    if(debug == 12){
      (void)fprintf(fp_err,"color+1=%i,r=%f,g=%f,b=%f,hue=%f,sat=%f,val=%f\n",
		    color+1,red,green,blue,hue,saturation,value);
    } /* end debug */
  } /* end loop over colors */

  /* set the two neutral colors to gray */ 
  c_gscr(1,num_colors/2+1,.8,.8,.8);
  c_gscr(1,num_colors/2+2,.8,.8,.8);

} /* end red_blue_scale() */ 

int zscale3(num_colors)
     int num_colors;
{
/* define a set of color triples for colors 2 through num_colors+1 */ 

  int color;

  float red;
  float green;
  float blue;

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,0.,0.,0.);
  /* this makes color 2 white */ 
  c_gscr(1,2,1.,1.,1.);
  
  for(color=3;color<=num_colors/3;color++){
    blue=1.;
    green=1.-((color-3.)/(num_colors/3.-3.));
    red=(color-3.)/(num_colors/3.-3.);
    c_gscr(1,color,red,green,blue);
  } /* end loop over colors */

  for(color=1;color<=(num_colors+3)/3;color++){
    red=1.;
    blue=1.-((color)/((num_colors+3)/3.));
    green=(color)/((num_colors+3)/3.);
    c_gscr(1,color+num_colors/3,red,green,blue);
  } /* end loop over colors */

  for(color=1;color<=(num_colors+3)/3;color++){
    green=1.;
    red=1.-((color)/((num_colors+3)/3.));
    blue=(color)/((num_colors+3)/3.);
    c_gscr(1,color+2*num_colors/3,red,green,blue);
  } /* end loop over colors */

} /* end zscale3() */ 

int prescale(num_colors)
     int num_colors;
{
/* define a set of grey scale color triples for colors 2 through 25 */ 

  static float triplet[NUM_DIST_DIFF_CONTOURS+2][3]={
    {1.,1.,1.},

    {.0,.0,.0},
    {1.,1.,1.},
    {1.,.75,1.}, /* color 3 */ 
    {1.,.38,1.},
    {1.,1.,.12},

    {1.,0.,0.},
    {1.,.63,.88},
    {.75,.38,.12},
    {.75,.88,.25},
    {.75,.88,1.},

    {.5,.63,1.},
    {.5,0.,1.},
    {.5,.38,.63},
    {.5,.5,0.},
    {.5,0.,0.},

    {.5,1.,0.},
    {.25,.63,.0},
    {.25,.5,.5},
    {.25,.12,1.},
    {.25,0.,.25},

    {0.,.75,1.},
    {0.,1.,.5},
    {0.,.63,.0},
    {0.,.12,.5},
    {0.,0.,0.}
  };

  int color;

  float red;
  float green;
  float blue;

  /* this makes color 0 white */ 
  c_gscr(1,0,1.,1.,1.);
  /* this makes color 1 black */ 
  c_gscr(1,1,.0,.0,.0);
  /* this makes color 2 white */ 
  c_gscr(1,2,1.,1.,1.);
  
  for(color=3;color<=num_colors+1;color++){
    red=triplet[color][0];
    green=triplet[color][1];
    blue=triplet[color][2];
    c_gscr(1,color,red,green,blue);
  } /* end loop over colors */

} /* end prescale() */ 
