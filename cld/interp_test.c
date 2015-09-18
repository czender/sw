#ifndef lint
  static char rcsid[] = "$Id$";
#endif

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:45  zender
/* Imported sources
/* */

#define Boolean int
#define True 1
#define False 0
#define YES 1
#define NO 0
#define FILESIZE 80

#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */

int debug=0; /* Option D */
int debug_value=0; /* Option d */

#define NUM_STD_ATM_LEVELS 38
#define NUM_CCM2_LEVELS 35

main(argc,argv)
     int argc;
     char * argv[];
{
  void spline_();

  void Exit_gracefully();
  void ratint();
  void polint();
  void print_usage();

  Boolean VERBOSE;

  char *time_buf_start;
  char out_file[FILESIZE];
  char in_file[FILESIZE];
  char err_file[FILESIZE];
  
  extern char *optarg;
  extern int optind;
  
  float *altitude;
  float *env_pressure;
  float *env_temperature;
  float *std_atm_alt;
  float *std_atm_pres;
  float *std_atm_temp;

  float dz;
  float float_foo;

  int int_foo;
  int layer;
  int num_layer;
  int num_std_atm_level;
  int opt;
  
  time_t clock;

  static float altitude_temperature_pressure[NUM_STD_ATM_LEVELS+2][3]={
    /* 45 degrees N July profile USSA 1966 p 115 */ 
    /* format is pressure (mb), temperature (K), altitude (Pa) */ 
    0.,0.,0., /* element [0] */ 

    4000.,273.57,6.280e4,
    4250.,272.04,6.087e4,
    4500.,270.51,5.899e4,
    4750.,268.98,5.715e4,

    5000.,267.45,5.536e4,
    5250.,265.92,5.362e4,
    5500.,264.39,5.192e4,
    5750.,262.86,5.027e4,

    6000.,261.33,4.866e4,
    6250.,259.70,4.709e4,
    6500.,258.07,4.557e4,
    6750.,256.44,4.408e4,
    
    7000.,254.81,4.264e4,
    7250.,253.17,4.123e4,
    7500.,251.54,3.986e4,
    7750.,249.91,3.853e4,

    8000.,248.28,3.724e4,
    8250.,246.65,3.598e4,
    8500.,245.03,3.475e4,
    8750.,243.40,3.356e4,

    9000.,241.77,3.240e4,
    9250.,240.15,3.128e4,
    9500.,238.52,3.018e4,
    9750.,236.90,2.912e4,

    10000.,235.27,2.809e4,
    10250.,233.65,2.709e4,
    10500.,232.02,2.612e4,
    10750.,230.40,2.517e4,

    11000.,228.77,2.426e4,
    11500.,225.53,2.250e4,

    12000.,222.30,2.086e4,
    12500.,219.06,1.931e4,

    13000.,215.82,1.786e4,
    13500.,215.65,1.650e4,
    
    14000.,215.65,1.525e4,
    14500.,215.65,1.409e4,

    15000.,215.65,1.525e4,
    15500.,215.65,1.409e4,
    
    0.,0.,0. /* element [NUM_STD_ATM_LEVELS+1] */
    };

  /* set defaults */
  num_layer=15.;
  num_std_atm_level=NUM_STD_ATM_LEVELS;

  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */
  
  /* parse command line arguments */
  while((opt = getopt(argc, argv, "D:e:i:l:o:Vv")) != EOF){
    switch(opt){
    case 'D':
      /* The debugging level.  Default is 0. */
      debug = (unsigned short int)atoi(optarg);
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_file,optarg);
      break;
    case 'i':
      /* get the input file name. Default is stdin */
      (void)strcpy(in_file,optarg);
      break;
    case 'l':
      /* Default is 15 */
      num_layer = atoi(optarg); 
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
  
  if(
     ((altitude=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_temperature=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((std_atm_alt=(float *)malloc((NUM_STD_ATM_LEVELS+2)*sizeof(float))) == NULL ) ||
     ((std_atm_pres=(float *)malloc((NUM_STD_ATM_LEVELS+2)*sizeof(float))) == NULL ) ||
     ((std_atm_temp=(float *)malloc((NUM_STD_ATM_LEVELS+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  for(layer=1;layer<=NUM_STD_ATM_LEVELS;layer++){
    std_atm_alt[layer]=altitude_temperature_pressure[layer][0];
    std_atm_temp[layer]=altitude_temperature_pressure[layer][1];
    std_atm_pres[layer]=altitude_temperature_pressure[layer][2];
  } /* end loop over layers */

  (void)fprintf(stdout,"Initial Thermodynamic Grid:\n");
  (void)fprintf(stdout,"layer\talt.\tT\tp\n");
  (void)fprintf(stdout,"\t km\t K\t mbar\n");
  for(layer=1;layer<=NUM_STD_ATM_LEVELS;layer++){
  (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,std_atm_alt[layer]/1000.,
		std_atm_temp[layer],std_atm_pres[layer]/100.);
  } /* end loop over layers */

  for(layer=1;layer<=num_layer;layer++){
    altitude[layer]=6000.+(layer-1)*6000./num_layer+3000./num_layer;
  } /* end loop over layers */

  int_foo=num_std_atm_level-1;
  /* set debug_value = 2 for spline debugging */ 
  spline_(std_atm_alt+1,std_atm_temp+1,&int_foo,altitude+1,
	  env_temperature+1,&num_layer,&debug_value);
  spline_(std_atm_alt+1,std_atm_pres+1,&int_foo,altitude+1,
	  env_pressure+1,&num_layer,&debug_value);

  (void)fprintf(stdout,"\n");
  (void)fprintf(stdout,"Spline fitting results spline.f:\n");
  (void)fprintf(stdout,"layer\talt.\tT\tp\n");
  (void)fprintf(stdout,"\t km\t K\t mbar\n");
  for(layer=1;layer<=num_layer;layer++){
    (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,altitude[layer]/1000.,
		  env_temperature[layer],env_pressure[layer]/100.);
  } /* end loop over layers */

  for(layer=1;layer<=num_layer;layer++){
    polint(std_atm_alt,std_atm_temp,num_std_atm_level,altitude[layer],
	   env_temperature+layer,&float_foo);
    polint(std_atm_alt,std_atm_pres,num_std_atm_level,altitude[layer],
	   env_pressure+layer,&float_foo);
  } /* end loop over layers */

  (void)fprintf(stdout,"\n");
  (void)fprintf(stdout,"Polynomial Interpolation polint.c results\n");
  (void)fprintf(stdout,"layer\talt.\tT\tp\n");
  (void)fprintf(stdout,"\t km\t K\t mbar\n");
  for(layer=1;layer<=num_layer;layer++){
    (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,altitude[layer]/1000.,
		  env_temperature[layer],env_pressure[layer]/100.);
  } /* end loop over layers */

  for(layer=1;layer<=num_layer;layer++){
    ratint(std_atm_alt,std_atm_temp,num_std_atm_level,altitude[layer],
	   env_temperature+layer,&float_foo);
    ratint(std_atm_alt,std_atm_pres,num_std_atm_level,altitude[layer],
	   env_pressure+layer,&float_foo);
  } /* end loop over layers */

  (void)fprintf(stdout,"\n");
  (void)fprintf(stdout,"Rational function interpolation ratint.c:\n");
  (void)fprintf(stdout,"layer\talt.\tT\tp\n");
  (void)fprintf(stdout,"\t km\t K\t mbar\n");
  for(layer=1;layer<=num_layer;layer++){
    (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,altitude[layer]/1000.,
		  env_temperature[layer],env_pressure[layer]/100.);
  } /* end loop over layers */

  Exit_gracefully();
} /* end main() */

void print_usage()
{
  (void)fprintf(stderr,"\nusage: cloud [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"D:e:i:o:vV\n\n");
  (void)fprintf(stderr,"-D debug The debugging level.  Default is 0.\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-i in_file get the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-l num_layer Default is 20\n");
  (void)fprintf(stderr,"-o out_file get the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"\n");
}

void Exit_gracefully()
{
  char *time_buf_finish;
  time_t clock;

  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stdout,"\tfinish = %s\n",time_buf_finish);

  exit(1);
}



