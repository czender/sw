/* $Id$ */

/* Purpose:  Test conversions between GMT time and UNIX time */ 

/* Example Usage (place mousable command lines here):
   time_tst 
   1995/10/30 18:51:47 GMT is 815079107 seconds since 00:00:00.00 Jan 1, 1970
   1995/01/01 00:00:00 GMT is 788918400 seconds since 00:00:00.00 Jan 1, 1970
*/ 

/* Standard C headers */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */

int main(int argc,char **argv)
{
  char *tz_lcl=NULL;
  char *tz_set_cmd=NULL;

  extern time_t timezone;
#ifndef LINUX
  extern time_t altzone;
#endif /* not LINUX */ 
  extern int daylight;
  extern char *tzname[2];

  int mth;
  int day;
  int yr;
  int hr;
  int mnt;
  int sec;
  
  long time_unix_long;

  struct tm *gmt_tm;

  time_t time_unix_time_t;

  /* Initialize time structure */ 
  time_unix_time_t=(time_t)0;
  gmt_tm=gmtime(&time_unix_time_t);
  
  yr=1995;
  mth=1;
  day=1;
  hr=0;
  mnt=0;
  sec=0;

  gmt_tm->tm_sec=sec;
  gmt_tm->tm_min=mnt;
  gmt_tm->tm_hour=hr;
  gmt_tm->tm_mday=day;
  gmt_tm->tm_mon=mth-1;
  gmt_tm->tm_year=yr-1900;
  
  /* NB: These time externs are set automatically in mktime().
     Calling timegm() does not set them. 
     The user can set the externs with tzset(). */ 
  /* Quoth the Solaris man page for ctime(): 
     "ctime(), localtime(), mktime(), and strftime() will also update these
     external variables as if they had called tzset() at the time
     specified by the time_t or struct tm  value  that  they  are
     converting." */

  /* "The external time_t variable altzone  contains  the  differ-
     ence, in seconds, between Coordinated Universal Time and the
     alternate  time  zone.   The  external   variable   timezone
     contains  the  difference, in seconds, between UTC and local
     standard time.  The  external  variable  daylight  indicates
     whether  time  should  reflect  daylight savings time.  Both
     timezone and altzone default to 0 (UTC).  The external vari-
     able  daylight is non-zero if an alternate time zone exists." */

  /* localtime() and gmtime() return pointers  to  tm  structures
     (see  below).   localtime()  corrects for the main time zone
     and possible alternate  (``daylight  savings'')  time  zone;
     gmtime()  converts  directly  to  Coordinated Universal Time
     (UTC), which is what the UNIX system uses internally. */

  (void)fprintf(stdout,"Original timezone environment:\n");
  tz_lcl=(char *)strdup(getenv("TZ"));
  (void)fprintf(stderr,"Local Time Zone = $TZ = tz_lcl = %s\n",tz_lcl);
  (void)tzset();
  if(1){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",timezone);
#ifndef LINUX
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",altzone);
#endif /* not LINUX */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   

#if ( defined LINUX )
  time_unix_long=(long)timegm(gmt_tm);
  (void)fprintf(stdout,"LINUX calling timegm():\n");
  (void)fprintf(stdout,"%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d is %li seconds since 00:00:00.00 Jan 1, 1970\n",yr,mth,day,hr,mnt,sec,time_unix_long);
  (void)putenv("TZ=GMT");
  (void)tzset();
  if(1){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",timezone);
#ifndef LINUX
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",altzone);
#endif /* not LINUX */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   
  time_unix_long=(long)mktime(gmt_tm);
  (void)fprintf(stdout,"LINUX calling mktime():\n");
  (void)fprintf(stdout,"%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d is %li seconds since 00:00:00.00 Jan 1, 1970\n",yr,mth,day,hr,mnt,sec,time_unix_long);
#else
  (void)putenv("TZ=GMT");
  (void)tzset();
  time_unix_long=(long)mktime(gmt_tm);
  (void)fprintf(stdout,"SUNMP calling mktime():\n");
  (void)fprintf(stdout,"%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d is %d seconds since 00:00:00.00 Jan 1, 1970\n",yr,mth,day,hr,mnt,sec,time_unix_long);
#endif

  /* Restore original timezone environment */ 
  (void)fprintf(stdout,"Restoring original timezone environment...\n");
  tz_set_cmd=(char *)malloc((strlen(tz_lcl)+4)*sizeof(char));
  (void)sprintf(tz_set_cmd,"TZ=%s",tz_lcl);
  (void)putenv(tz_set_cmd);
  (void)tzset();
  if(1){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",timezone);
#ifndef LINUX
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",altzone);
#endif /* not LINUX */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   
  if(tz_lcl != NULL) (void)free(tz_lcl); tz_lcl=NULL;
  if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;

  /*  (void)fprintf(stdout,"%d seconds since 00:00:00.00 Jan 1, 1970\n",sec);*/

  exit(0);
} /* end main() */ 


