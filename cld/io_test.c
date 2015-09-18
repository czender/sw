static char rcsid[] = "$Id$\n";

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:42  zender
/* Imported sources
/*
 * Revision 4.12  1993/08/24  16:47:08  zender
 * about to branch to a cray compatible version which works
 * with default ncar graphics version 3.2 installation. this
 * mainly means changing to the GKS c-structures and changes
 * in the library location in the makefiles.
 *
 * Revision 4.10  1993/06/11  04:41:41  zender
 * changed CRAY2 to CRAY, i think io_diff is way different, not sure!
 *
 * Revision 1.1  1993/06/11  04:39:11  zender
 * Initial revision
 *
 * Revision 4.2  1993/03/08  00:57:47  zender
 * this version contains working working data i/o in 3 formats:
 * ASCII, CZ_BINARY (my vanilla binary format), and NETCDF.
 * since netCDF is transparent across platforms (i.e., cray and sun),
 * the next version will only allow netCDF i/o.
 *
 * Revision 4.1  1993/02/27  00:05:18  zender
 * incorporated ascii writing automatically on the CRAY in main,
 * merged ioinput and iooutput into iographs. added albedo vs.
 * emissivity graph. alphabetized output. might work on netCDF now.
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

#define NUM_IWC_CONTOURS 21
#define NUM_DIST_CONTOURS 24

static int iasf[13]={
  1,1,1,1,1,1,1,1,1,1,1,1,1};

static char *IWC_llbs[]={
  "0.",
  ".1",".5","1.","2.","5.",
  "10.","15.","20.","25.","30.",
  "40.","50.","60.","70.","80.",
  "90.","100.","125.","150.","175.",
  "200."};

static float IWC_contour_level[NUM_IWC_CONTOURS+1]={
  0.,
  .1,.5,1.,2.,5.,
  10.,15.,20.,25.,30.,
  40.,50.,60.,70.,80.,
  90.,100.,125.,150.,175.,
  200.};

static int IWC_crosshatch_lind[NUM_IWC_CONTOURS]={
  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};

static int IWC_grey_lind[NUM_IWC_CONTOURS]={
  2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};

static char *dist_llbs[]={
  "0.",
  ".03125",".0625",".125",
  ".25",".5","1","2","4",
  "8","16","32","64","128",
  "256","512","1k","2k","4k",
  "8k","16k","32k","64k","128k",
  "256k"};

static float dist_contour_level[NUM_DIST_CONTOURS+1]={
  0.,
  .03125,.0625,.125,
  .25,.5,1.,2.,4.,
  8.,16.,32.,64.,128.,
  256.,512.,1024.,2048.,4096.,
  8192.,16384.,32768.,65536.,131072.,
  262144.};

static int dist_crosshatch_lind[NUM_DIST_CONTOURS]={
  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};

static int dist_grey_lind[NUM_DIST_CONTOURS]={
  2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

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
  char in_file[FILESIZE];
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
  float **daltitude_dtime;

  float * IWC_snapshot_ptr;
  float * number_snapshot_ptr;
  float * layer_frame_ptr;
  float * layer_size_ptr;
  float * concentration_ptr;
  float * distribution_ptr;
  float * crystal_heat_LW_ptr;
  float * crystal_heat_SW_ptr;
  float * daltitude_dtime_ptr;

  float *SW_cloud_forcing;
  float *LW_cloud_forcing;
  float *cloud_optical_depth;
  float *IWP_of_time;
  float *altitude;
  float *cloud_base_of_time;
  float *cloud_effective_radius;
  float *cloud_top_of_time;
  float *crystal_length;
  float *crystal_mass;
  float *delta_length;
  float *delta_mass;
  float *effective_length;
  float *effective_radius;
  float *env_density;
  float *env_pressure;
  float *env_temperature;
  float *flux_up_SW;
  float *flux_down_SW;
  float *flux_up_LW;
  float *flux_down_LW;
  float *heating_rate_LW;
  float *heating_rate_SW;
  float *heating_rate_net;
  float *albedo_of_time;
  float *emissivity_of_time;
  float *IWC;
  float *IWP;
  float *latent_heating;
  float *max_length_of_time;
  float *mmr_vapor;
  float *orig_env_temp;
  float *orig_mmr_vapor;
  float *potential_temp;
  float *pp_vapor;
  float *rwrk;
  float *saturation_ice_of_time;
  float *saturation_ice;
  float *advective_heating;
  float *surface_area_of_time;
  float *time_array;
  float *time_snapshot;
  float *vapor_density;
  float *vapor_path;
  float *vapor_path_of_time;
  float *optical_depth;
  float *water_path_of_time;
  float *wind_speed;
  
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
  int int_foo2;
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
  int cdfid;  

  /* dimension sizes (which are longs) and array to hold them */ 
  long count[2];
  long num_layerp2,num_bandp2,num_framep1,num_sizep2,num_stepp1;

  /* indicial starting offsets (which are longs) */
  long start[2];
  
  /* scalar variable ids */
  int dt_id, dz_id, num_band_id, num_ccm2_level_id, num_cloudy_layers_id, num_frame_id, num_layer_id, num_size_id, num_step_id, agg_step_id, rad_step_id, plot_step_id, PAUSE_id;

  /* two dimensional variable ids */
  int IWC_snapshot_id, number_snapshot_id, concentration_id, distribution_id, crystal_heat_SW_id, crystal_heat_LW_id, daltitude_dtime_id;

  /* one dimensional variable ids */
  int IWC_id, IWP_id, IWP_of_time_id, LW_cloud_forcing_id, SW_cloud_forcing_id, advective_heating_id, albedo_of_time_id, altitude_id, cloud_base_of_time_id, cloud_effective_radius_id, cloud_top_of_time_id, crystal_length_id, crystal_mass_id, delta_length_id, delta_mass_id, effective_radius_id, emissivity_of_time_id, env_density_id, env_pressure_id, env_temperature_id, flux_down_LW_id, flux_down_SW_id, flux_up_LW_id, flux_up_SW_id, heating_rate_LW_id, heating_rate_SW_id, heating_rate_net_id, latent_heating_id, max_length_of_time_id, mmr_vapor_id, cloud_optical_depth_id, orig_env_temp_id, orig_mmr_vapor_id, potential_temp_id, pp_vapor_id, saturation_ice_of_time_id, saturation_ice_id, surface_area_of_time_id, time_array_id, time_snapshot_id, vapor_density_id, vapor_path_id, vapor_path_of_time_id, optical_depth_id, water_path_of_time_id, wind_speed_id;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = False; /* Option I */
  STDOUT = False; /* Option O */
  NETCDF_FORMAT = True;
  VERBOSE = True; /* print out WARNINGS? */ /* Option V */

  (void)strcpy(err_file,"stderr"); /* Option e */
  (void)strcpy(in_file,"cloud.nc"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  
  /* Set GKS variables */
  lam = 200000; /* size of NCAR iam array */ 
  nra = 20000; /* size of NCAR xra and yra arrays */ 
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
  while((opt = getopt(argc, argv,"D:d:Ee:Ii:Oo:Vv"))!= EOF){
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
      /* get the input file name. Default is cloud.nc */
      (void)strcpy(in_file,optarg);
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
      /* print the RCS program version */
      (void)fprintf(stderr,rcsid);
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
  
  if(STDIN){
    fp_in = stdin;
  }else{
    if( (fp_in = fopen( in_file, "r")) == NULL) {
      (void)fprintf(stderr,"\nError in opening input file %s\n",in_file);
      exit(1);
    } /* end if */
  } /* end if */
  
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

  if(NETCDF_FORMAT){
    /* first get the scalar values and dimensions for malloc'ing */ 
    start[0]=0L;
    cdfid = ncopen(in_file, NC_NOWRITE);

    ncattget(cdfid,NC_GLOBAL,"cmdline",(void *)cmdline);
    dt_id=ncvarid(cdfid,"dt");
    ncvarget1(cdfid,dt_id,start,(void *)&dt);
    dz_id=ncvarid(cdfid,"dz");
    ncvarget1(cdfid,dz_id,start,(void *)&dz);
    num_band_id=ncvarid(cdfid,"num_band");
    ncvarget1(cdfid,num_band_id,start,(void *)&num_band);
    num_ccm2_level_id=ncvarid(cdfid,"num_ccm2_level");
    ncvarget1(cdfid,num_ccm2_level_id,start,(void *)&num_ccm2_level);
    num_cloudy_layers_id=ncvarid(cdfid,"num_cloudy_layers");
    ncvarget1(cdfid,num_cloudy_layers_id,start,(void *)&num_cloudy_layers);
    num_frame_id=ncvarid(cdfid,"num_frame");
    ncvarget1(cdfid,num_frame_id,start,(void *)&num_frame);
    num_layer_id=ncvarid(cdfid,"num_layer");
    ncvarget1(cdfid,num_layer_id,start,(void *)&num_layer);
    num_size_id=ncvarid(cdfid,"num_size");
    ncvarget1(cdfid,num_size_id,start,(void *)&num_size);
    num_step_id=ncvarid(cdfid,"num_step");
    ncvarget1(cdfid,num_step_id,start,(void *)&num_step);
    agg_step_id=ncvarid(cdfid,"agg_step");
    ncvarget1(cdfid,agg_step_id,start,(void *)&agg_step);
    rad_step_id=ncvarid(cdfid,"rad_step");
    ncvarget1(cdfid,rad_step_id,start,(void *)&rad_step);
    plot_step_id=ncvarid(cdfid,"plot_step");
    ncvarget1(cdfid,plot_step_id,start,(void *)&plot_step);
    PAUSE_id=ncvarid(cdfid,"PAUSE");
    ncvarget1(cdfid,PAUSE_id,start,(void *)&PAUSE);
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

  /* allocate storage for all the dynamic input arrays */
  if(
     ((IWC=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWP_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((LW_cloud_forcing=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((SW_cloud_forcing=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((advective_heating=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((albedo_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((altitude=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_base_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_effective_radius=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_top_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((crystal_length=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_mass=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_length=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_mass=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((effective_length=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((effective_radius=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((emissivity_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((env_density=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_temperature=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_net=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((latent_heating=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((max_length_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((mmr_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_optical_depth=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((orig_env_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((orig_mmr_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((potential_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((pp_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((saturation_ice_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((saturation_ice=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((surface_area_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_array=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_snapshot=(float *)malloc((num_frame+1)*sizeof(float))) == NULL ) ||
     ((vapor_density=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((optical_depth=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((water_path_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((wind_speed=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
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
     ((daltitude_dtime_ptr=
       (float *)malloc((num_size+2)*(num_layer+2)*sizeof(float))) == NULL ) ||
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
     ((daltitude_dtime=
       (float **)malloc((num_size+2)*sizeof(float *))) == NULL ) ||
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

  /* REMEMBER: these are the only matrices indexed as [size][layer] */
  for(size=0;size<num_size+2;size++){
    daltitude_dtime[size]=daltitude_dtime_ptr+(num_layer+2)*size;
  } /* end loop over sizes */
  
  /* REMEMBER: these are the only snug matrices indexed as [layer][frame] */
  for(layer=0;layer<=num_layer-1;layer++){
    layer_frame[layer]=layer_frame_ptr+(num_frame+1)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these are the only matrices indexed as [frame][layer] */
  for(frame=0;frame<=num_frame;frame++){
    IWC_snapshot[frame]=IWC_snapshot_ptr+(num_layer+2)*frame;    
    number_snapshot[frame]=number_snapshot_ptr+(num_layer+2)*frame;    
  } /* end loop over frames */
  
  /* Now that appropriate space has been allocated, read in the arrays */
  if(NETCDF_FORMAT){
    /* retrieve the one-dimensional arays */ 
    start[0]=0L;

    /* define sizes of dimensions */
    num_layerp2 = (long)(num_layer+2);
    num_bandp2 = (long)(num_band+2);
    num_framep1 = (long)(num_frame+1);
    num_sizep2 = (long)(num_size+2);
    num_stepp1 = (long)(num_step+1);

    IWC_id=ncvarid(cdfid,"IWC");
    ncvarget(cdfid,IWC_id,start,&num_layerp2,(void *)IWC);
    IWP_id=ncvarid(cdfid,"IWP");
    ncvarget(cdfid,IWP_id,start,&num_layerp2,(void *)IWP);
    IWP_of_time_id=ncvarid(cdfid,"IWP_of_time");
    ncvarget(cdfid,IWP_of_time_id,start,&num_stepp1,(void *)IWP_of_time);
    LW_cloud_forcing_id=ncvarid(cdfid,"LW_cloud_forcing");
    ncvarget(cdfid,LW_cloud_forcing_id,start,&num_stepp1,(void *)LW_cloud_forcing);
    SW_cloud_forcing_id=ncvarid(cdfid,"SW_cloud_forcing");
    ncvarget(cdfid,SW_cloud_forcing_id,start,&num_stepp1,(void *)SW_cloud_forcing);
    advective_heating_id=ncvarid(cdfid,"advective_heating");
    ncvarget(cdfid,advective_heating_id,start,&num_layerp2,(void *)advective_heating);
    albedo_of_time_id=ncvarid(cdfid,"albedo_of_time");
    ncvarget(cdfid,albedo_of_time_id,start,&num_stepp1,(void *)albedo_of_time);
    altitude_id=ncvarid(cdfid,"altitude");
    ncvarget(cdfid,altitude_id,start,&num_layerp2,(void *)altitude);
    cloud_base_of_time_id=ncvarid(cdfid,"cloud_base_of_time");
    ncvarget(cdfid,cloud_base_of_time_id,start,&num_stepp1,(void *)cloud_base_of_time);
    cloud_effective_radius_id=ncvarid(cdfid,"cloud_effective_radius");
    ncvarget(cdfid,cloud_effective_radius_id,start,&num_stepp1,(void *)cloud_effective_radius);
    cloud_top_of_time_id=ncvarid(cdfid,"cloud_top_of_time");
    ncvarget(cdfid,cloud_top_of_time_id,start,&num_stepp1,(void *)cloud_top_of_time);
    crystal_length_id=ncvarid(cdfid,"crystal_length");
    ncvarget(cdfid,crystal_length_id,start,&num_sizep2,(void *)crystal_length);
    crystal_mass_id=ncvarid(cdfid,"crystal_mass");
    ncvarget(cdfid,crystal_mass_id,start,&num_sizep2,(void *)crystal_mass);
    delta_length_id=ncvarid(cdfid,"delta_length");
    ncvarget(cdfid,delta_length_id,start,&num_sizep2,(void *)delta_length);
    delta_mass_id=ncvarid(cdfid,"delta_mass");
    ncvarget(cdfid,delta_mass_id,start,&num_sizep2,(void *)delta_mass);
    effective_radius_id=ncvarid(cdfid,"effective_radius");
    ncvarget(cdfid,effective_radius_id,start,&num_layerp2,(void *)effective_radius);
    emissivity_of_time_id=ncvarid(cdfid,"emissivity_of_time");
    ncvarget(cdfid,emissivity_of_time_id,start,&num_stepp1,(void *)emissivity_of_time);
    env_density_id=ncvarid(cdfid,"env_density");
    ncvarget(cdfid,env_density_id,start,&num_layerp2,(void *)env_density);
    env_pressure_id=ncvarid(cdfid,"env_pressure");
    ncvarget(cdfid,env_pressure_id,start,&num_layerp2,(void *)env_pressure);
    env_temperature_id=ncvarid(cdfid,"env_temperature");
    ncvarget(cdfid,env_temperature_id,start,&num_layerp2,(void *)env_temperature);
    flux_down_LW_id=ncvarid(cdfid,"flux_down_LW");
    ncvarget(cdfid,flux_down_LW_id,start,&num_layerp2,(void *)flux_down_LW);
    flux_down_SW_id=ncvarid(cdfid,"flux_down_SW");
    ncvarget(cdfid,flux_down_SW_id,start,&num_layerp2,(void *)flux_down_SW);
    flux_up_LW_id=ncvarid(cdfid,"flux_up_LW");
    ncvarget(cdfid,flux_up_LW_id,start,&num_layerp2,(void *)flux_up_LW);
    flux_up_SW_id=ncvarid(cdfid,"flux_up_SW");
    ncvarget(cdfid,flux_up_SW_id,start,&num_layerp2,(void *)flux_up_SW);
    heating_rate_LW_id=ncvarid(cdfid,"heating_rate_LW");
    ncvarget(cdfid,heating_rate_LW_id,start,&num_layerp2,(void *)heating_rate_LW);
    heating_rate_SW_id=ncvarid(cdfid,"heating_rate_SW");
    ncvarget(cdfid,heating_rate_SW_id,start,&num_layerp2,(void *)heating_rate_SW);
    heating_rate_net_id=ncvarid(cdfid,"heating_rate_net");
    ncvarget(cdfid,heating_rate_net_id,start,&num_layerp2,(void *)heating_rate_net);
    latent_heating_id=ncvarid(cdfid,"latent_heating");
    ncvarget(cdfid,latent_heating_id,start,&num_layerp2,(void *)latent_heating);
    max_length_of_time_id=ncvarid(cdfid,"max_length_of_time");
    ncvarget(cdfid,max_length_of_time_id,start,&num_stepp1,(void *)max_length_of_time);
    mmr_vapor_id=ncvarid(cdfid,"mmr_vapor");
    ncvarget(cdfid,mmr_vapor_id,start,&num_layerp2,(void *)mmr_vapor);
    cloud_optical_depth_id=ncvarid(cdfid,"cloud_optical_depth");
    ncvarget(cdfid,cloud_optical_depth_id,start,&num_stepp1,(void *)cloud_optical_depth);
    orig_env_temp_id=ncvarid(cdfid,"orig_env_temp");
    ncvarget(cdfid,orig_env_temp_id,start,&num_layerp2,(void *)orig_env_temp);
    orig_mmr_vapor_id=ncvarid(cdfid,"orig_mmr_vapor");
    ncvarget(cdfid,orig_mmr_vapor_id,start,&num_layerp2,(void *)orig_mmr_vapor);
    potential_temp_id=ncvarid(cdfid,"potential_temp");
    ncvarget(cdfid,potential_temp_id,start,&num_layerp2,(void *)potential_temp);
    pp_vapor_id=ncvarid(cdfid,"pp_vapor");
    ncvarget(cdfid,pp_vapor_id,start,&num_layerp2,(void *)pp_vapor);
    saturation_ice_of_time_id=ncvarid(cdfid,"saturation_ice_of_time");
    ncvarget(cdfid,saturation_ice_of_time_id,start,&num_stepp1,(void *)saturation_ice_of_time);
    saturation_ice_id=ncvarid(cdfid,"saturation_ice");
    ncvarget(cdfid,saturation_ice_id,start,&num_layerp2,(void *)saturation_ice);
    surface_area_of_time_id=ncvarid(cdfid,"surface_area_of_time");
    ncvarget(cdfid,surface_area_of_time_id,start,&num_stepp1,(void *)surface_area_of_time);
    time_array_id=ncvarid(cdfid,"time_array");
    ncvarget(cdfid,time_array_id,start,&num_stepp1,(void *)time_array);
    time_snapshot_id=ncvarid(cdfid,"time_snapshot");
    ncvarget(cdfid,time_snapshot_id,start,&num_framep1,(void *)time_snapshot);
    vapor_density_id=ncvarid(cdfid,"vapor_density");
    ncvarget(cdfid,vapor_density_id,start,&num_layerp2,(void *)vapor_density);
    vapor_path_id=ncvarid(cdfid,"vapor_path");
    ncvarget(cdfid,vapor_path_id,start,&num_layerp2,(void *)vapor_path);
    vapor_path_of_time_id=ncvarid(cdfid,"vapor_path_of_time");
    ncvarget(cdfid,vapor_path_of_time_id,start,&num_stepp1,(void *)vapor_path_of_time);
    optical_depth_id=ncvarid(cdfid,"optical_depth");
    ncvarget(cdfid,optical_depth_id,start,&num_layerp2,(void *)optical_depth);
    water_path_of_time_id=ncvarid(cdfid,"water_path_of_time");
    ncvarget(cdfid,water_path_of_time_id,start,&num_stepp1,(void *)water_path_of_time);
    wind_speed_id=ncvarid(cdfid,"wind_speed");
    ncvarget(cdfid,wind_speed_id,start,&num_layerp2,(void *)wind_speed);

    /* now get the two-dimensional arays */ 
    start[0]=0L;
    start[1]=0L;

    count[0] = num_framep1;
    count[1] = num_layerp2;
    IWC_snapshot_id=ncvarid(cdfid,"IWC_snapshot");
    ncvarget(cdfid,IWC_snapshot_id,start,count,(void *)IWC_snapshot_ptr);
    count[0] = num_framep1;
    count[1] = num_layerp2;
    number_snapshot_id=ncvarid(cdfid,"number_snapshot");
    ncvarget(cdfid,number_snapshot_id,start,count,(void *)number_snapshot_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    concentration_id=ncvarid(cdfid,"concentration");
    ncvarget(cdfid,concentration_id,start,count,(void *)concentration_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    distribution_id=ncvarid(cdfid,"distribution");
    ncvarget(cdfid,distribution_id,start,count,(void *)distribution_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    crystal_heat_SW_id=ncvarid(cdfid,"crystal_heat_SW");
    ncvarget(cdfid,crystal_heat_SW_id,start,count,(void *)crystal_heat_SW_ptr);
    count[0] = num_layerp2;
    count[1] = num_sizep2;
    crystal_heat_LW_id=ncvarid(cdfid,"crystal_heat_LW");
    ncvarget(cdfid,crystal_heat_LW_id,start,count,(void *)crystal_heat_LW_ptr);
    count[0] = num_sizep2;
    count[1] = num_layerp2;
    daltitude_dtime_id=ncvarid(cdfid,"daltitude_dtime");
    ncvarget(cdfid,daltitude_dtime_id,start,count,(void *)daltitude_dtime_ptr);

    ncclose(cdfid);
    (void)fprintf(fp_err,"Read netCDF format cloud data from %s\n",in_file);
  } /* end if NETCDF_FORMAT */
  
  if(debug == 55){
    for(frame=0;frame<=num_frame;frame++){
      for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(fp_err,"IWC_snapshot[%i][%i] = %g mg/m^3\n",
		      frame,layer,IWC_snapshot[frame][layer]*MILLIGRAMS_PER_KILOGRAM);
      } /* end loop over layers */
    } /* end loop over frames */
  } /* end debug */
  
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
    /* draw the shortwave, longwave, and net heating rates */
    for(layer=1;layer<=num_layer;layer++){
      abscissa[layer]=heating_rate_net[layer]*SECONDS_PER_HOUR; /* K/hour */
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_layer);
    abscissa_min=min_VEC(abscissa+1,num_layer);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    c_plchhq(.5,.5,"Heating Rates",-.7,0.,0.);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array[num_step]*MINUTES_PER_SECOND); 
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
/*    c_labmod("(F6.2)","(F5.1)",0,0,10,10,0,15,0); */
    /* set the best domain for multiple plots of same dimensions */ 
/*    xaxis_min=abscissa_min;*/
/*    xaxis_max=abscissa_max;*/
    /* world coordinates around the plot region */
/*    c_set(.1,.9,.1,.9,xaxis_min,xaxis_max,*/
/*	  ordinate_min,ordinate_max,1);*/

    /* arrange for the tick spacing and actually draw the axes */
/*    c_gridal(10,4,9,4,1,1,5,0,0); */
    /* print the header */
    float_foo=ordinate_min+.3*(ordinate_max-ordinate_min);
    (void)sprintf(line_label,"Abscissa [Minimum, Maximum] (units)");
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    /* actually draw in the abscissa = 0 line */
    c_line(0.,ordinate_min,0.,ordinate_max);
    
    /* set the dash characteristics of the line */
/*    (void)sprintf(line_label,"$$$$$$$$$$$$$$$$$Net");*/
/*    c_dashdc(line_label,jcrt,jsize);*/
    /* draw the curve */
/*    c_curved(abscissa+1,ordinate+1,num_layer);*/
    (void)sprintf(line_label,"Net Radiative Heating Rate [%.1f, %.1f] K/hour",
		  abscissa_min,abscissa_max);
    c_plchhq(xaxis_min+.3*(xaxis_max-xaxis_min),float_foo,line_label,-.7,0.,-1.);
    float_foo-=.04*(ordinate_max-ordinate_min);
    
    /* advance the frame */
    c_frame();
  } /* end if */ 

  if(debug == 15){
    /* draw the IWP vs Reff scatterplot using conpack routines */
    for(step=0;step<=num_step;step++){
      abscissa[step]=IWP_of_time[step]*1000.; /* g/m^2 */
      ordinate[step]=cloud_effective_radius[step]*MICRONS_PER_METER; /* microns */
    } /* end loop over layers */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step+1);
    ordinate_min=min_VEC(ordinate,num_step+1);
    abscissa_max=max_VEC(abscissa,num_step+1);
    abscissa_min=min_VEC(abscissa,num_step+1);
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,
	  abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    if(debug_value == 1){
      c_set(.1,.9,.1,.9,
	    0.,abscissa_max,
	    20.,55.,1);
    } /* end if */
    if(debug_value == 2){
      c_set(.1,.9,.1,.9,
	    0.,75.,
	    20.,42.,1);
    } /* end if */
    /*	abscissa_min,abscissa_max,*/
    /*	ordinate_min,ordinate_max,1);*/
    /* plot the points */
    c_points(abscissa,ordinate,num_step+1,-2,0);
    /* advance the frame */
    c_frame();
    
    c_gdawk(wkid);
    c_gclwk(wkid);
    c_gclks(); 
    
    Exit_gracefully();
    exit(0);
  } /* end debug */
  
  if(True){
    /* draw the emissivity vs albedo scatterplot using conpack routines */
    num_rad_steps=num_step/rad_step;
    for(step=0;step<=num_step;step++){
      abscissa[step]=emissivity_of_time[step];
      ordinate[step]=albedo_of_time[step];
    } /* end loop over steps */
    /* Normalize the ordinates */
    ordinate_max=max_VEC(ordinate,num_step+1);
    ordinate_min=min_VEC(ordinate,num_step+1);
    abscissa_max=max_VEC(abscissa,num_step+1);
    abscissa_min=min_VEC(abscissa,num_step+1);

/*    num_rad_steps=num_step/rad_step;*/
/*    for(step=0;step<=num_rad_steps;step++){*/
/*      abscissa[step]=emissivity_of_time[step*rad_step];*/
/*      ordinate[step]=albedo_of_time[step*rad_step];*/
/*    } *//* end loop over steps */
    /* Normalize the ordinates */
/*    ordinate_max=max_VEC(ordinate,num_rad_steps+1);*/
/*    ordinate_min=min_VEC(ordinate,num_rad_steps+1);*/
/*    abscissa_max=max_VEC(abscissa,num_rad_steps+1);*/
/*    abscissa_min=min_VEC(abscissa,num_rad_steps+1);*/
    if(debug != 73){
      /* normalize coordinates around the title region */
      c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
      /* plotchar the title of the graph */
      (void)sprintf(line_label,"Emissivity vs. Albedo");
      c_plchhq(.5,.5,line_label,-.7,0.,0.);
      (void)fprintf(fp_err,"Plotting...%s\n",line_label);
      /* plotchar the time the cloud has evolved */
      (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		    time_array[num_step]*MINUTES_PER_SECOND); 
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
    } /* end debug */
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,
	  0.,1.,0.,.75,1);
    if(debug != 73){
      /* arrange for the tick spacing and actually draw the axes */
      c_gridal(10,4,10,4,1,1,0,0,0); 
    } /* end debug */
    /* change the size of the polymarkers */ 
/*    c_gsmksc(3.);*/
    /* plot the points */
    debug_value = (debug_value == 0) ? 1 : debug_value;
/*    (void)fprintf(fp_err,"print debug_value = %i\n",debug_value);*/
    c_points(abscissa,ordinate,num_step+1,-debug_value,0);
/*    c_points(abscissa,ordinate,num_rad_steps+1,-2,0);*/
    /* reset the size of the polymarkers */ 
/*    c_gsmksc(1.);*/
    
    /* advance the frame */
    c_frame();
  } /* end if */

  /* draw the IWP vs Reff scatterplot using conpack routines */
  for(step=0;step<=num_step;step++){
    abscissa[step]=IWP_of_time[step]*1000.; /* g/m^2 */
    ordinate[step]=cloud_effective_radius[step]*MICRONS_PER_METER; /* microns */
  } /* end loop over layers */
  /* Normalize the ordinates */
  ordinate_max=max_VEC(ordinate,num_step+1);
  ordinate_min=min_VEC(ordinate,num_step+1);
  abscissa_max=max_VEC(abscissa,num_step+1);
  abscissa_min=min_VEC(abscissa,num_step+1);
  /* normalize coordinates around the title region */
  c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
  /* plotchar the title of the graph */
  c_plchhq(.5,.5,"Ice Water Path IWP vs. Effective Radius r:B:eff:N:",-.7,0.,0.);
  /* plotchar the time the cloud has evolved */
  (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		time_array[num_step]/60.); 
  c_plchhq(.5,.2,line_label,-.7,0.,0.); 
  /* plotchar the command line */
  c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
  /* normalize coordinates around the x-axis region */
  c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
  c_plchhq(.5,.5,"IWP (g/m:S:2:N:)",-.7,0.,0.);
  /* normalize coordinates around the y-axis region */
  c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
  c_plchhq(.3,.5,"Effective Radius r:B3:eff    (:PGL1:L:R:)",-.7,90.,0.);
  /* format the axis labels */
  c_labmod("(F5.1)","(F5.1)",0,0,10,10,0,15,0); 
  /* world coordinates around the plot region */
  c_set(.1,.9,.1,.9,
	abscissa_min,abscissa_max,
	ordinate_min,ordinate_max,1);
  if(debug_value == 1){
    c_set(.1,.9,.1,.9,
	  0.,abscissa_max,
	  20.,55.,1);
  } /* end debug */
  if(debug_value == 2){
    c_set(.1,.9,.1,.9,
	  0.,75.,
	  20.,42.,1);
  } /* end debug */
  /* arrange for the tick spacing and actually draw the axes */
  c_gridal(10,4,9,4,1,1,5,0,0); 
  /* plot the points */
  c_points(abscissa,ordinate,num_step+1,-2,0);
  /* advance the frame */
  c_frame();

  /* draw the IWC vs. Temperature plot using conpack routines */
  for(layer=1;layer<=num_layer;layer++){
    abscissa[layer]=env_temperature[layer];
    ordinate[layer]=IWC[layer]*MILLIGRAMS_PER_KILOGRAM; /* mg/m^3 */
  } /* end loop over layers */
  /* Normalize the abscissas */
  abscissa_max=max_VEC(abscissa+1,num_layer);
  abscissa_min=min_VEC(abscissa+1,num_layer);
  ordinate_max=max_VEC(ordinate+1,num_layer);
  ordinate_min=min_VEC(ordinate+1,num_layer);
  /* normalize coordinates around the title region */
  c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
  /* plotchar the title of the graph */
  c_plchhq(.5,.5,"Ice Water Content IWC vs. Temperature T",-.7,0.,0.);
  /* plotchar the time the cloud has evolved */
  (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		time_array[num_step]*MINUTES_PER_SECOND); 
  c_plchhq(.5,.2,line_label,-.7,0.,0.); 
  /* plotchar the command line */
  c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
  /* normalize coordinates around the x-axis region */
  c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
  c_plchhq(.5,.5,"Temperature T (K)",-.7,0.,0.);
  /* normalize coordinates around the y-axis region */
  c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
  c_plchhq(.3,.5,"IWC (mg/m:S:3:N:)",-.7,90.,0.);
  /* world coordinates around the plot region */
  c_set(.1,.9,.1,.9,
	abscissa_min,abscissa_max,ordinate_min,ordinate_max,1);
  /* format the axis labels */
  c_labmod("(F5.1)","(F5.1)",0,0,10,10,0,15,0); 
  /* arrange for the tick spacing and actually draw the axes */
  c_gridal(10,4,9,4,1,1,5,0,0); 
  /* plot the points */
  c_points(abscissa,ordinate,num_layer,-2,0);
  /* advance the frame */
  c_frame();

  /* draw the IWC_snapshot vs. Temperature plot using conpack routines */
  for(layer=1;layer<=num_layer;layer++){
    abscissa[layer]=env_temperature[layer];
    ordinate[layer]=IWC[layer]*MILLIGRAMS_PER_KILOGRAM; /* mg/m^3 */
  } /* end loop over layers */
  /* Normalize the abscissas */
  abscissa_max=max_VEC(abscissa+1,num_layer);
  abscissa_min=min_VEC(abscissa+1,num_layer);
  ordinate_max=max_VEC(ordinate+1,num_layer);
  ordinate_min=min_VEC(ordinate+1,num_layer);
  /* normalize coordinates around the title region */
  c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
  /* plotchar the title of the graph */
  c_plchhq(.5,.5,"Time Evolution of Ice Water Content IWC vs. Temperature T",
	   -.7,0.,0.);
  /* plotchar the time the cloud has evolved */
  (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		time_array[num_step]*MINUTES_PER_SECOND); 
  c_plchhq(.5,.2,line_label,-.7,0.,0.); 
  /* plotchar the command line */
  c_plchhq(.1,.05,cmdline,-.4,0.,-1.); 
  /* normalize coordinates around the x-axis region */
  c_set(0.,1.,0.,.1,0.,1.,0.,1.,1);
  c_plchhq(.5,.5,"Temperature T (K)",-.7,0.,0.);
  /* normalize coordinates around the y-axis region */
  c_set(0.,.1,0.,1.,0.,1.,0.,1.,1);
  c_plchhq(.3,.5,"IWC (mg/m:S:3:N:)",-.7,90.,0.);
  /* world coordinates around the plot region */
  c_set(.1,.9,.1,.9,
	abscissa_min,abscissa_max,ordinate_min,ordinate_max,1);
  /* format the axis labels */
  c_labmod("(F5.1)","(F5.1)",0,0,10,10,0,15,0); 
  /* arrange for the tick spacing and actually draw the axes */
  c_gridal(10,4,9,4,1,1,5,0,0); 
  /* plot the points */
  c_points(abscissa,ordinate,num_layer,-2,0);

  /* Overplot the same points from other times in the evolution */
  int_foo=0;
  for(frame=1;frame<num_frame;frame+=10){
    for(layer=1;layer<=num_layer;layer++){
      /* assume the temperature in each layer was static */ 
      ordinate[layer]=IWC_snapshot[frame][layer]*MILLIGRAMS_PER_KILOGRAM; /* mg/m^3 */
    } /* end loop over layers */
    /* plot the points */
    c_points(abscissa,ordinate,num_layer,(--int_foo)%6,0);
  } /* end loop over frames */
  /* advance the frame */
  c_frame();

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* contour plot the crystal_heat_LW */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert J/s to nJ/s and transpose */ 
	layer_size[layer-1][size-1]=NANOWATTS_PER_WATT*
	  crystal_heat_LW[layer][size];
      } /* end loop over sizes */
    } /* end loop over layers */
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    c_plchhq(.5,.5,"Longwave Heat Budget of an Individual Crystal",-.7,0.,0.);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array[num_step]*MINUTES_PER_SECOND); 
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
    c_set(.1,.9,.1,.9,0.,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=32;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      int_foo=-num_contour_level/2.+level;
      if(int_foo < 0){
	float_foo2=-1.*pow(2.,(float)(-int_foo));
      }else{
	float_foo2=pow(2.,(float)int_foo);
      } /* end else */
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f nJ/s\n",level,float_foo2);
      } /* end debug */
      (void)sprintf(line_label,"%i",ROUNDER(float_foo2));
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
    (void)sprintf(line_label,"Values (in nano-Watts) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* contour plot the crystal_heat_SW */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert J/s to nJ/s and transpose */ 
	layer_size[layer-1][size-1]=NANOWATTS_PER_WATT*
	  crystal_heat_SW[layer][size];
      } /* end loop over sizes */
    } /* end loop over layers */
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    c_plchhq(.5,.5,"Shortwave Heat Budget of an Individual Crystal",-.7,0.,0.);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array[num_step]*MINUTES_PER_SECOND); 
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
    c_set(.1,.9,.1,.9,0.,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=32;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      int_foo=-num_contour_level/2.+level;
      if(int_foo < 0){
	float_foo2=-1.*pow(2.,(float)(-int_foo));
      }else{
	float_foo2=pow(2.,(float)int_foo);
      } /* end else */
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f nJ/s\n",level,float_foo2);
      } /* end debug */
      (void)sprintf(line_label,"%i",ROUNDER(float_foo2));
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
    (void)sprintf(line_label,"Values (in nano-Watts) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* contour plot the crystal_heat_net */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert J/s to nJ/s and transpose */ 
	layer_size[layer-1][size-1]=NANOWATTS_PER_WATT*
	  (crystal_heat_LW[layer][size]+crystal_heat_SW[layer][size]);
      } /* end loop over sizes */
    } /* end loop over layers */
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    c_plchhq(.5,.5,"Net Heat Budget of an Individual Crystal",-.7,0.,0.);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array[num_step]*MINUTES_PER_SECOND); 
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
    c_set(.1,.9,.1,.9,0.,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=32;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      int_foo=-num_contour_level/2.+level;
      if(int_foo < 0){
	float_foo2=-1.*pow(2.,(float)(-int_foo));
      }else{
	float_foo2=pow(2.,(float)int_foo);
      } /* end else */
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f nJ/s\n",level,float_foo2);
      } /* end debug */
      (void)sprintf(line_label,"%i",ROUNDER(float_foo2));
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
    (void)sprintf(line_label,"Values (in nano-Watts) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */

  if((char_ptr_foo=strstr(cmdline,"-R")) == NULL){
    /* contour plot the crystal_heat_net*concentration */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert J/s to nJ/s and transpose */ 
	layer_size[layer-1][size-1]=MICROWATTS_PER_WATT*
	  (crystal_heat_LW[layer][size]+crystal_heat_SW[layer][size])*
	    concentration[layer][size];
      } /* end loop over sizes */
    } /* end loop over layers */
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    c_plchhq(.5,.5,"Integrated Heat Budget of Crystal Distribution",-.7,0.,0.);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Elapsed time :PGU1:D:RU:t = %.1f min.", 
		  time_array[num_step]*MINUTES_PER_SECOND); 
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
    c_set(.1,.9,.1,.9,0.,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=32;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      int_foo=-num_contour_level/2.+level;
      if(int_foo < 0){
	float_foo2=-1.*pow(2.,(float)(-int_foo));
      }else{
	float_foo2=pow(2.,(float)int_foo);
      } /* end else */
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f nJ/s\n",level,float_foo2);
      } /* end debug */
      (void)sprintf(line_label,"%i",ROUNDER(float_foo2));
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
    (void)sprintf(line_label,"Values (in :PGL1:L:R:Watts/m:S:3:N:) from $ZMN$ to $ZMX$ contoured by factors of 2.");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
  } /* end if */

  if(True){
    /* contour plot the fall speed daltitude_dtime */ 
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* convert m/s to cm/s and plot define falling to be positive */ 
	layer_size[layer-1][size-1]=
	  -CENTIMETERS_PER_METER*daltitude_dtime[size][layer];
      } /* end loop over sizes */
    } /* end loop over layers */
    for(size=1;size<=num_size;size++){
      abscissa[size]=crystal_length[size]*MICRONS_PER_METER;
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      ordinate[layer]=altitude[layer]*KILOMETERS_PER_METER; /* km */
    } /* end loop over layers */
    /* Normalize the abscissas */
    abscissa_max=max_VEC(abscissa+1,num_size);
    abscissa_min=min_VEC(abscissa+1,num_size);
    ordinate_max=max_VEC(ordinate+1,num_layer);
    ordinate_min=min_VEC(ordinate+1,num_layer);
    /* normalize coordinates around the title region */
    c_set(0.,1.,.9,1.,0.,1.,0.,1.,1);
    /* plotchar the title of the graph */
    (void)sprintf(line_label,"Net Fall Speed of Crystal Distribution");
    c_plchhq(.5,.5,line_label,-.7,0.,0.);
    (void)fprintf(fp_err,"Plotting...%s\n",line_label);
    /* plotchar the time the cloud has evolved */
    (void)sprintf(line_label,"Maximum Speed in Updraft Profile = %.1f cm/s", 
		  max_VEC(wind_speed+1,num_layer)*CENTIMETERS_PER_METER); 
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
    abscissa_max=1000.;
    /* world coordinates around the plot region */
    c_set(.1,.9,.1,.9,abscissa_min,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* arrange for the tick spacing and actually draw the axes */
    c_gridal(10,4,9,4,1,1,5,0,0); 
    
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
    c_set(.1,.9,.1,.9,0.,abscissa_max,
	  ordinate_min,ordinate_max,1);
    /* set the number of contour levels to be drawn */ 
    num_contour_level=50;
    c_cpseti("NCL -- number contour levels",num_contour_level);
    /* reset the contour levels */ 
    for(level=1;level<=num_contour_level;level++){
      int_foo=-num_contour_level/2.+level;
      int_foo2=(int_foo > 0) ? int_foo : -int_foo;
      if(int_foo2 < 6){
	float_foo2=int_foo2*2.;
      }else{
	float_foo2=10.+(int_foo2-5)*5.; /* contour every 5 cm/s */ 
      } /* end else */
      float_foo2=(int_foo > 0) ? float_foo2 : -float_foo2;
      if(debug == 59){
	(void)fprintf(fp_err,"contour level #%i at %.1f nJ/s\n",level,float_foo2);
      } /* end debug */
      (void)sprintf(line_label,"%i",ROUNDER(float_foo2));
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
    (void)sprintf(line_label,"Values (in cm/s) from $ZMN$ to $ZMX$ contoured by levels of 2 and then 5");
    c_cpsetc("ILT -- info label text string",line_label);
    c_cpsetr("ILY -- y coord of info label",.05);
    c_cplbdr(layer_size_ptr,rwrk,iwrk);
    /* advance the frame */
    c_frame();
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

void Exit_gracefully()
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
  (void)fprintf(stderr,"D:d:Ee:Ii:Oo:vV\n\n");
  (void)fprintf(stderr,"-D debug Set the debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is False\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"-d debug_value Default is 0\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-i in_file get the input file name. Default is stdin\n");
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
    if(ifll > 0 && ifll < NUM_DIST_CONTOURS){
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

  static float triplet[NUM_DIST_CONTOURS+2][3]={
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
