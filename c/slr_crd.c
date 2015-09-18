#define Boolean int
#define True 1
#define False 0
#define YES 1
#define NO 0
#define FILESIZE 80
#define CMDLINE_SIZE 200

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

/* Personal headers */

FILE *fp_out,*fp_err,*fp_in;

/* Global variables declared here */
int dbg=0; /* Option D */
int dbg_value=0; /* Option d */

int main(int argc,char **argv)
{
  /* call with something like . . .
     ./sun -n 1000 -l 70. -s 150. -f 240. -D 1
     */

  void Exit_gracefully(void);
  void print_usage(void);
  void string_cmdline(int,char *[],char *,int);

  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *time_buf_start;
  char cmdline[CMDLINE_SIZE];
  char out_file[FILESIZE];
  char in_file[FILESIZE];
  char err_file[FILESIZE];
  char CVS_Id[] = "$Id$\n";
  char CVS_Revision[] = "$Revision$\n";
  
  extern char *optarg;
  
  int opt;
  
  time_t clock;

  void solar_geometry(float,float,int,float *,float *,float *);
  float d2r(float);
  float r2d(float);

  float *cosSZA;
  float *local_time;

  float calendar_day_of_year;
  float dt_days;
  float eccentricity_factor;
  float finish_cal_day;
  float latitude_deg;
  float latitude_rad;
  float mean_SZA;
  float mean_cosSZA;
  float start_cal_day;
  float daylight_SZA_total;
  float daylight_cosSZA_total;

  int longitude;
  int num_longitudes;
  int num_steps;
  int step;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = True; /* Option I */
  STDOUT = True; /* Option O */
  
  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */
  
  latitude_deg=0.; /* Option l */
  num_longitudes=1;
  num_steps=2000; /* Option n */
  start_cal_day=80.; /* Option s */
  finish_cal_day=90.; /* Option f */

  /* parse command line arguments */
  while((opt = getopt(argc, argv, "D:d:Ee:f:Ii:l:n:Oo:s:Vv")) != EOF){
    switch(opt){
    case 'f':
      /* The finish calendar day.  Default is 90. */
      finish_cal_day = (float)atof(optarg);
      break;
    case 's':
      /* The start calendar day.  Default is 80. */
      start_cal_day = (float)atof(optarg);
      break;
    case 'n':
      /* The number of quadrature points.  Default is 100 */
      num_steps = (int)atoi(optarg);
      break;
    case 'l':
      /* The latitude.  Default is 0. */
      latitude_deg = (float)atof(optarg);
      break;
    case 'D':
      /* The debugging level.  Default is 0. */
      dbg = (unsigned short int)atoi(optarg);
      break;
    case 'd':
      /* The second debugging level.  Default is 0. */
      dbg_value = (unsigned short int)atoi(optarg);
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
      /* Toggle the input file stream. Default is True */
      STDIN = !STDIN;
      break;
    case 'i':
      /* get the input file name. Default is stdin */
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
      /* print the CVS program version */
      (void)fprintf(stderr,"%s %s\n",CVS_Revision,CVS_Id);
      exit(EXIT_SUCCESS);
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
  
  if(STDERR){
    fp_err = stderr;
  }else{
    if( (fp_err = fopen( err_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening error file %s\n",err_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  if(STDIN){
    fp_in = stdin;
  }else{
    if( (fp_in = fopen( in_file, "r")) == NULL) {
      (void)fprintf(stderr,"\nError in opening input file %s\n",in_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  if(STDOUT){
    fp_out = stdout;
  }else{
    if( (fp_out = fopen( out_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening output file %s\n",out_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  /* start the clock and save the command line */
  string_cmdline( argc, argv, cmdline, CMDLINE_SIZE );
  (void)fprintf(fp_err,"Command Line: %s\n",cmdline);
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stderr,"\tstart = %s",time_buf_start);


  latitude_rad=d2r(latitude_deg);
  dt_days=(finish_cal_day-start_cal_day)/num_steps;

  if(
     ((local_time=(float *)malloc((num_longitudes+2)*sizeof(float))) == NULL ) ||
     ((cosSZA=(float *)malloc((num_longitudes+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  daylight_SZA_total=0.;
  daylight_cosSZA_total=0.;
  calendar_day_of_year=start_cal_day;
  for(step=1;step<=num_steps;step++){
    solar_geometry
      (latitude_rad,
       calendar_day_of_year,
       num_longitudes,
       local_time,
       cosSZA,
       &eccentricity_factor);

    /* make sure you're not counting darkness */
    if(cosSZA[1] >= 0.){
      daylight_SZA_total+=acos(cosSZA[1]);
      daylight_cosSZA_total+=cosSZA[1];
    } /* end if */

    if(dbg == 1){
      for(longitude=1;longitude<=num_longitudes;longitude++){
	(void)fprintf(fp_err,"calendar_day_of_year = %g\n",calendar_day_of_year);
	(void)fprintf(fp_err,"local_time = %g\n",local_time[longitude]);
	(void)fprintf(fp_err,"cosSZA = %g, i.e., %g degrees\n",cosSZA[longitude],r2d(acos(cosSZA[longitude])));
	(void)fprintf(fp_err,"mean SZA so far is = %g degrees\n",r2d(daylight_SZA_total)/step);
	(void)fprintf(fp_err,"mean cosSZA so far is = %g, i.e., %g degrees\n",daylight_cosSZA_total/step,r2d(acos(daylight_cosSZA_total/step)));
	(void)fprintf(fp_err,"\n");
      } /* end loop over longitudes */
    } /* end dbg */
    
    calendar_day_of_year+=dt_days;
  } /* end loop over steps */
  mean_SZA=daylight_SZA_total/num_steps;
  mean_cosSZA=daylight_cosSZA_total/num_steps;

  (void)fprintf(fp_err,"num_longitudes = %i\n",num_longitudes);
  (void)fprintf(fp_err,"latitude_rad = %g radians\n",latitude_rad);
  (void)fprintf(fp_err,"latitude_deg = %g degrees\n",latitude_deg);
  (void)fprintf(fp_err,"\n");

  (void)fprintf(fp_err,"num_steps = %i\n",num_steps);
  (void)fprintf(fp_err,"start_cal_day = %g\n",start_cal_day);
  (void)fprintf(fp_err,"finish_cal_day = %g\n",finish_cal_day);
  (void)fprintf(fp_err,"dt_days = %g\n",dt_days);
  (void)fprintf(fp_err,"\n");

  (void)fprintf(fp_err,"eccentricity_factor = %g\n",eccentricity_factor);
  (void)fprintf(fp_err,"mean_SZA = %g degrees\n",r2d(mean_SZA));
  (void)fprintf(fp_err,"mean_cosSZA = %g, i.e., %g degrees\n",mean_cosSZA,
		r2d(acos(mean_cosSZA)));

  Exit_gracefully();
  return EXIT_SUCCESS;
} /* end main() */

void print_usage(void)
{
  (void)fprintf(stderr,"\nusage: clouds [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"D:d:Ee:Ii:Oo:vV\n\n");
  (void)fprintf(stderr,"-D dbg The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-d dbg_value The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-i in_file get the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-o out_file get the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"\n");
} /* end print_usage() */

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

  exit(0);
} /* end Exit_gracefully() */

/* Module:	cmdparse.c (Command Parse)
 * Purpose:	Set options from command line
 * Subroutine:	parse_cmdline()			returns: int
 * Subroutine:	string_cmdline()		returns: void
 * Subroutine:	usage()				returns: int
 * Xlib calls:	none
 * Copyright:	1989 Smithsonian Astrophysical Observatory
 *		You may do anything you like with this file except remove
 *		this copyright.  The Smithsonian Astrophysical Observatory
 *		makes no representations about the suitability of this
 *		software for any purpose.  It is provided "as is" without
 *		express or implied warranty.
 * Modified:	{0} Michael VanHilst	initial version	       5 January 1989
 *              {1} MVH BSDonly strings.h compatability           19 Feb 1990
 *		{n} <who> -- <does what> -- <when>
 */

void string_cmdline(int argc,char *argv[],char *cmdline,int linemax )
{
  int i;

  if( argc <= 0 ) {
    cmdline[0] = '\0';
  } else {
    (void)strcpy(cmdline, argv[0]);
    for( i=1; i<argc; i++ ) {
      (void)strncat(cmdline, " ", linemax);
      (void)strncat(cmdline, argv[i], linemax);
    }
  }
} /* end string_cmdline() */

#define DAYS_PER_YEAR 365.
void solar_geometry
  (float latitude_rad,
   float calendar_day_of_year,
   int num_longitudes,
   float *local_time,
   float *cosSZA,
   float *eccentricity_factor)
/* compute solar zenith angle and other diagnositics from 
   the time of year and latitude: this is a hack of cmpsol.f
   from the CCM2 code by Kiehl and Briegleb. */
{
  float cos_lat;
  float cos_delta;
  float cphase;
  float delta;
  float phi;
  float sin_lat;
  float sin_delta;
  float theta;
  
  int longitude;

  /* compute eccentricity factor (sun-earth distance factor) */
  theta=2.*M_PI*calendar_day_of_year/DAYS_PER_YEAR;
  *eccentricity_factor=1.000110+.034221*cos(theta)+.001280*sin(theta)+
    .000719*cos(2.*theta)+.000077*sin(2.*theta);

  /* solar declination in radians: */
  delta=.006918-.399912*cos(theta)+.070257*sin(theta)-
    .006758*cos(2.*theta)+.000907*sin(2.*theta)-
      .002697*cos(3.*theta)+.001480*sin(3.*theta);

  /* compute local cosine solar zenith angle: */
  sin_lat=sin(latitude_rad);
  sin_delta=sin(delta);
  cos_lat=cos(latitude_rad);
  cos_delta=cos(delta);

  /* calendar_day_of_year is the calender day for greenwich, including fraction
     of day; the fraction of the day represents a local time at
     greenwich; to adjust this to produce a true instantaneous time
     for other longitudes, we must correct for the local time change: */

  for(longitude=1;longitude<=num_longitudes;longitude++){
    phi=calendar_day_of_year+((float)(longitude-1)/(float)(num_longitudes));
    cphase=cos(2.*M_PI*phi);
    *(cosSZA+longitude)=sin_lat*sin_delta-cos_lat*cos_delta*cphase;
    *(local_time+longitude)=12.*acos(cphase)/M_PI;
  } /* end loop over longitudes */

} /* end solar_geometry() */

float d2r(float theta_degrees)
{
  return theta_degrees*M_PI/180.;
} /* end d2r */

float r2d(float theta_radians)
{
  return theta_radians*180./M_PI;
} /* end r2d */

