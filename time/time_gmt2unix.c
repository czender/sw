/* $Id$ */

/* Purpose:  Convert GMT time string to UNIX time */ 

/* Usage:
   GMT time should be specified in this format:
   YYYYMMDD HH:MM:SS
   time_gmt2unix 19700101 00:00:00
   time_gmt2unix 19950101 00:00:00
   time_gmt2unix 19500101 00:00:00
   time_gmt2unix 19640312 12:09:00

   If YYYYMMDD HH:MM:SS is not specified on command line, then built in default is used:

   time_gmt2unix 
   1995/10/30 18:51:47 GMT is 815079107 seconds since 00:00:00.00 Jan 1, 1970
   1995/01/01 00:00:00 GMT is 788918400 seconds since 00:00:00.00 Jan 1, 1970
*/ 

/* Standard C headers */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */
#include <unistd.h>             /* all sorts of POSIX stuff */

int main(int argc,char **argv)
{
  char *tz_sys=NULL; /* [sng] Pointer to TZ variable in environment */
  char *tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */
  char *tz_set_cmd=NULL; /* [sng] Command for putenv() to set TZ */
  char sbr_nm[]="time_gmt2unix.c"; /* [sng] Subroutine name */
  char *yyyymmdd_sng;
  char *hhmmss_sng;

#ifdef AIX
  extern long timezone; /* [mnt] Offset of current TZ from GMT */
#else /* not AIX */ 
  extern time_t timezone; /* [mnt] Offset of current TZ from GMT */
#endif /* not AIX */ 
#ifndef LINUX
  extern time_t altzone; /* [mnt] Offset of alternate TZ from GMT */
#endif /* not LINUX */ 
  extern int daylight; /* [flg] Daylight savings time */
  extern char *tzname[2]; /* [sng] UNIX code of main and alternate TZs */

  int mth; /* [mth] 1-based month of year [1..12] */
  int day; /* [day] 1-based day of month [1..31] */
  int yr; /* [yr] Christian Year */
  int hr; /* [hr] 0-based, hour of day [0..23] */
  int mnt; /* [mnt] 0-based minute of second [0..59] */
  int sec; /* [s] 0-based second of minute [0..59] */
  int rcd=0; /* [rcd] Return code */
  
  long yyyymmdd_lng;
  long time_unix_long; /* [s] UNIX time stored as intrinsic long */

  struct tm *gmt_tm; /* [tm] GMT time structure for desired input date */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  /* Set defaults */
  yr=1995;
  mth=1; /* January is 1 */
  day=1; /* The first is 1 */
  hr=0; /* 0-based, 24 hour */
  mnt=0;
  sec=0;

  /* Get command line time if any */
  if(argc == 3){
    yyyymmdd_sng=strdup(argv[1]);
    hhmmss_sng=strdup(argv[2]);
    /* Parse these strings */
    yyyymmdd_lng=atol(yyyymmdd_sng);
    day=yyyymmdd_lng%100;
    mth=((yyyymmdd_lng%10000)-day)/100;
    yr=(int)(yyyymmdd_lng/10000);

    /* Overwrite colons so that hh, mm, and ss components of hhmmss_sng are all NUL-terminated */
    hhmmss_sng[2]='\0';
    hhmmss_sng[5]='\0';
    hr=(int)atol(hhmmss_sng);
    mnt=(int)atol(hhmmss_sng+3);
    sec=(int)atol(hhmmss_sng+6);
    
    if(hhmmss_sng != NULL){free(hhmmss_sng); hhmmss_sng=NULL;}
    if(yyyymmdd_sng != NULL){free(yyyymmdd_sng); yyyymmdd_sng=NULL;}
  } /* endif */ 

  /* Save original timezone environment */ 
  tz_sys=getenv("TZ"); /* [sng] Pointer to TZ variable in environment */
  if(tz_sys != NULL) tz_lcl=(char *)strdup(tz_sys); else (void)fprintf(stderr,"WARNING: %s unable to find environment variable TZ\n",sbr_nm); /* [sng] Local copy of TZ environment variable */

  /* Change to GMT */ 
  if(tz_sys != NULL){
    char tz_gmt_set_cmd[]="TZ=GMT"; /* [sng] Command to copy to use in putenv() */
    /* Allocate space for command and terminating NUL */
    tz_set_cmd=(char *)malloc((strlen(tz_gmt_set_cmd)+1)*sizeof(char)); /* [sng] Command for putenv() to set TZ */
    (void)sprintf(tz_set_cmd,"%s",tz_gmt_set_cmd);
    rcd+=putenv(tz_set_cmd);
    (void)tzset();
    /* Do not free tz_set_cmd as as per caveats for Linux putenv() */
    /*    if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;*/
  } /* endif TZ exists */

  /* Initialize time structure */ 
  time_unix_time_t=(time_t)0; /* [s] UNIX time stored as time_t */
  gmt_tm=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */
  
  gmt_tm->tm_sec=sec; /* [s] 0-based second of minute [0..59] */
  gmt_tm->tm_min=mnt; /* [mnt] 0-based minute of second [0..59] */
  gmt_tm->tm_hour=hr; /* [hr] 0-based, hour of day [0..23] */
  gmt_tm->tm_mday=day; /* [day] 1-based day of month [1..31] */
  gmt_tm->tm_mon=mth-1; /* [mth] 0-based month [0..11] */
  gmt_tm->tm_year=yr-1900; /* [yr] Christian Year - 1900 */
  
  /* mktime() ignores ydy going in this direction */ 
  time_unix_long=(long)mktime(gmt_tm); /* [s] UNIX time stored as intrinsic long */

  if(1){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",(long)timezone);
#ifndef LINUX
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",(long)altzone);
#endif /* not LINUX */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   

  /* Restore original timezone environment */ 
  if(tz_sys != NULL){
    /* Allocate space for "TZ=" + tz_lcl + terminating NUL */
    tz_set_cmd=(char *)malloc((strlen(tz_lcl)+4)*sizeof(char)); /* [sng] Command for putenv() to set TZ */
    (void)sprintf(tz_set_cmd,"TZ=%s",tz_lcl); /* [sng] Command for putenv() to set TZ */
    rcd+=putenv(tz_set_cmd);
    (void)tzset();
    /* Do not free tz_set_cmd as per caveats on Linux putenv() man page */
    /*    if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;*/
  } /* endif TZ exists */

  if(tz_lcl != NULL) (void)free(tz_lcl); tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */

  (void)fprintf(stdout,"%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d GMT is %ld seconds since 00:00:00.00 Jan 1, 1970\n",yr,mth,day,hr,mnt,sec,time_unix_long);
  /*  (void)fprintf(stdout,"%d seconds since 00:00:00.00 Jan 1, 1970\n",sec);*/

  return EXIT_SUCCESS;
} /* end main() */ 
