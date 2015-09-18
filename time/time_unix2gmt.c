/* $Id$ */

/* Purpose: Convert UNIX time to GMT string */ 

/* Usage:
   time_unix2gmt 815079107.6
   1995/10/30 DOY 303 GMT 18:51:47.600000
 */ 

/* Standard C headers */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* all sorts of POSIX stuff */

int main(int argc,char **argv)
{
  double sec_frc; /* [s] Fractional seconds */
  double time_unix_dbl; /* [s] UNIX time stored as double */
  double hr_dcm; /* [frc] Decimal hours */

  int mth; /* [mth] 1-based month of year [1..12] */
  int day; /* [day] 1-based day of month [1..31] */
  int ydy; /* [day] 1-based day of year [1..366] */
  int yr; /* [yr] Christian Year */
  int hr; /* [hr] 0-based, hour of day [0..23] */
  int mnt; /* [mnt] 0-based minute of second [0..59] */
  int sec; /* [s] 0-based second of minute [0..59] */

  long time_unix_lng; /* [s] UNIX time stored as intrinsic long */

  struct tm *gmt_tm_ptr; /* [tm] UNIX time structure */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  if(argc < 2){
     fprintf(stderr,"\nUsage:  %s <unix_time> \n",argv[0]);
     exit(1);
  } /* endif */ 

  time_unix_dbl=atof(argv[1]); /* [s] UNIX time stored as double */
  time_unix_lng=(long)floor(time_unix_dbl); /* [s] UNIX time stored as intrinsic long */
  sec_frc=time_unix_dbl-time_unix_lng; /* [s] Fractional seconds */
  time_unix_time_t=(time_t)time_unix_lng; /* [s] UNIX time stored as time_t */

  gmt_tm_ptr=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */

  mth=gmt_tm_ptr->tm_mon+1; /* [mth] 1-based month of year [1..12] */
  day=gmt_tm_ptr->tm_mday; /* [day] 1-based day of month [1..31] */
  ydy=gmt_tm_ptr->tm_yday+1; /* [day] 1-based day of year [1..366] */
  yr=gmt_tm_ptr->tm_year+1900; /* [yr] Christian Year */
  hr=gmt_tm_ptr->tm_hour; /* [hr] 0-based, hour of day [0..23] */
  mnt=gmt_tm_ptr->tm_min; /* [mnt] 0-based minute of second [0..59] */
  sec=gmt_tm_ptr->tm_sec; /* [s] 0-based second of minute [0..59] */

  hr_dcm=gmt_tm_ptr->tm_hour+gmt_tm_ptr->tm_min/60.0+gmt_tm_ptr->tm_sec/3600.0; /* [frc] Decimal hours */
  hr_dcm+=0; /* [frc] Decimal hours */

  (void)fprintf(stdout,"%4.4d/%2.2d/%2.2d YDY %d GMT %2.2d:%2.2d:%f\n",yr,mth,day,ydy,hr,mnt,sec+sec_frc);

  return EXIT_SUCCESS;
} /* end main() */ 




