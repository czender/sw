/* $Id$ */

/* Purpose: Fortran interface to C time routines */

/* Copyright (C) 1994--2009 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* NB: Fortran calling C naming convention for LINUX is bizarre:
   Subroutines passing a string need two underscores, without a string only one */ 

/* NB: Compile this with C not C++ compiler since symbols must be read by Fortran */

/* Source: Routine date_time_() rewritten from datetime.c by Brian Eaton */ 

/* Usage:
   make ${MY_OBJ_DIR}/date_time.o 2>&1 | more

   Compilation:
   gcc -Wall -c -O -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o
   gcc -Wall -c -g -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o
   pgcc -c -g -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o

   pathcc -c -g -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o

   pgcc -c -g -DPGI_CC -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o

   xlc_r -c -O -q64 -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o

   cc -64 -c -O -D${PVM_ARCH} ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o
   gcc -DSGIMP64 -mabi=64 -O2 -Wall -c ${HOME}/sw/c/date_time.c -o ${MY_OBJ_DIR}/date_time.o */ 

/* Standard C headers */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt, putenv */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* all sorts of POSIX stuff */

#undef UNDERSCORE_ZERO
#undef UNDERSCORE_ONE
#undef UNDERSCORE_TWO

#if defined(AIX)
#define UNDERSCORE_ZERO
#else /* not AIX */
#define UNDERSCORE_ONE
#endif /* not AIX */

#ifdef UNDERSCORE_ZERO
#define FNC_C2F(x) x
#endif /* !UNDERSCORE_ZERO */
#ifdef UNDERSCORE_ONE
#define FNC_C2F(x) x##_
#endif /* !UNDERSCORE_ONE */
#ifdef UNDERSCORE_TWO
#define FNC_C2F(x) x##__
#endif /* !UNDERSCORE_TWO */

void
FNC_C2F(date_time)
     (char *bfr, /* O [sng] Date/time string */
      int bfr_lng) /* I [nbr] Length of bfr */
{
  /* Purpose: Packs a fortran-supplied input buffer with date/time in format 
     Fri Sep 16 13:11:08 1994 GMT
     Fortran calling conventions are assumed, arguments are passed as pointers, strings have buffer lengths tacked onto end of argument list
     Usage:
     Declare in Fortran with, e.g.,
     character*26 lcl_date_time      ! Store current time
     external date_time            ! Wrapper for C ctime() routine
     
     Call from Fortran with, e.g.,
     call date_time(lcl_date_time)
     write (0,'(a13,a26)') 'Start time = ',lcl_date_time
  */
  
  char *chr_ptr; /* [ptr] Pointer to current character */
  
  time_t time_time_t; /* [s] Time stored as time_t */

  /* Guard against overflow */
  if(bfr_lng < 26){
    (void)fprintf(stderr,"WARNING time buffer too small in date_time()\n");
    return;
  } /* endif */ 
  
  time_time_t=time((time_t *)NULL); /* [s] Time stored as time_t */
  (void)strcpy(bfr,ctime(&time_time_t)); /* [sng] Date/time string */
  
  /* Remove newline and pack remainder of bfr with blanks */
  chr_ptr=&bfr[24]; /* [ptr] Pointer to current character */
  while(chr_ptr-bfr < bfr_lng) 
    *chr_ptr++=' '; /* [ptr] Pointer to current character */

  return;
} /* end date_time() */ 

void 
/* [fnc] Compute UNIX time (seconds since 1969) of given GMT date-time */
FNC_C2F(gmt2unix)
     (int *yr, /* I [yr] Christian Year */
      int *mth, /* I [mth] 1-based month of year [1..12] */
      int *day, /* I [day] 1-based day of month [1..31] */
      int *hr, /* I [hr] 0-based, hour of day [0..23] */
      int *mnt, /* I [mnt] 0-based minute of second [0..59] */
      int *sec, /* I [s] 0-based second of minute [0..59] */
      int *time_unix_long) /* O [s] Pointer to UNIX time */
{
  /* Purpose: Compute UNIX time corresponding to given broken down GMT date-time structure
     Fortran calling conventions are assumed, arguments are passed as pointers, strings have buffer lengths tacked onto end of argument list
     Usage:
     c     Get UNIX time offset of beginning of year
     gmt_yr=1995
     int_foo=0
     int_foo2=1                ! Recall gmt2unix uses 0-based months and days of month
     call gmt2unix(gmt_yr,int_foo2,int_foo2,int_foo,int_foo,int_foo,time_unix_yr_srt)
  */

  char *tz_sys=NULL; /* [sng] Pointer to TZ variable in environment */
  char *tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */
  char *tz_set_cmd=NULL; /* [sng] Command to set time-zone */
  char sbr_nm[]="gmt2unix_"; /* [sng] Subroutine name */

#if defined(AIX) || defined(ALPHA)
  extern long timezone; /* [mnt] Offset of current TZ from GMT */
#else /* not AIX */ 
  extern time_t timezone; /* [mnt] Offset of current TZ from GMT */
#endif /* not AIX */ 
#if !defined(LINUX) && !defined(LINUXAMD64) && !defined(ALPHA)
  extern time_t altzone; /* [mnt] Offset of alternate TZ from GMT */
#endif /* not LINUX or ALPHA */ 
  extern int daylight; /* [flg] Daylight savings time */
  extern char *tzname[2]; /* [sng] UNIX code of main and alternate TZs */

  int rcd=0; /* [rcd] Return code */

  struct tm *gmt_tm; /* [tm] GMT time structure for desired input date */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  /* Save original timezone environment */ 
  tz_sys=getenv("TZ"); /* [sng] Pointer to TZ variable in environment */
  if(tz_sys != NULL) tz_lcl=(char *)strdup(tz_sys); else (void)fprintf(stderr,"WARNING: %s unable to find environment variable TZ\n",sbr_nm); /* [sng] Local copy of TZ environment variable */

  /* Change to GMT */ 
  if(tz_sys != NULL){
    char tz_gmt_set_cmd[]="TZ=GMT"; /* [sng] Command to copy to use in putenv() */
    /* Allocate space for command and terminating NUL */
    tz_set_cmd=(char *)malloc((strlen(tz_gmt_set_cmd)+1)*sizeof(char));
    (void)sprintf(tz_set_cmd,"%s",tz_gmt_set_cmd);
    /* 20060913: CentOS 4.1 needs __USE_SVID in stdlib.h to prototype putenv()
       Implicit declaration is fine since function returns int
       Better long-term workaround may be to switch to BSD setenv() */
    rcd+=putenv(tz_set_cmd);
    (void)tzset();
    /* Do not free tz_set_cmd as per caveats on Linux putenv() man page */
    /*    if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;*/
  } /* endif TZ exists */

  /* Initialize time structure */ 
  time_unix_time_t=(time_t)0; /* [s] UNIX time stored as time_t */
  gmt_tm=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */
  
  gmt_tm->tm_sec=*sec; /* [s] 0-based second of minute [0..59] */
  gmt_tm->tm_min=*mnt; /* [mnt] 0-based minute of second [0..59] */
  gmt_tm->tm_hour=*hr; /* [hr] 0-based, hour of day [0..23] */
  gmt_tm->tm_mday=*day; /* [day] 1-based day of month [1..31] */
  gmt_tm->tm_mon=*mth-1; /* [mth] 0-based month [0..11] */
  gmt_tm->tm_year=*yr-1900; /* [yr] Christian Year - 1900 */

  /* mktime() ignores ydy going in this direction */ 
  *time_unix_long=(long)mktime(gmt_tm); /* [s] UNIX time stored as intrinsic long */

  if(0){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",(long)timezone);
#if !defined(LINUX) && !defined(LINUXAMD64) && !defined(ALPHA)
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",(long)altzone);
#endif /* not LINUX or ALPHA */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   

  /* Warning: Behavior of putenv under Linux is very subtle---read man page for latest on 
     re-entrancy, memory leaks, use of automatic variable and BSD compatibility of this routine 
     Safest policy is to malloc() memory here, not to free it, and to allow memory leak */

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

  if(tz_lcl != NULL) (void)free(tz_lcl);
  tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */

  return;
} /* end gmt2unix() */ 

void 
FNC_C2F(unix2gmt)
     (double *time_unix, /* I [s] Pointer to UNIX time stored in double precision */
      int *yr, /* O [yr] Christian Year */
      int *ydy, /* O [day] 1-based day of year [1..366] */
      int *mth, /* O [mth] 1-based month of year [1..12] */
      int *day, /* O [day] 1-based day of month [1..31] */
      int *hr, /* O [hr] 0-based, hour of day [0..23] */
      int *mnt, /* O [mnt] 0-based minute of hour [0..59] */
      int *sec) /* O [s] 0-based second of minute [0..59] */
{
  /* Purpose: 
     Fortran calling conventions are assumed, arguments are passed as pointers, strings have buffer lengths tacked onto end of argument list
     Usage: 
     integer gmt_day
     integer gmt_hr
     integer gmt_mnt
     integer gmt_mth
     integer gmt_sec
     integer gmt_ydy
     integer gmt_yr
     double precision time_unix
     character*4 tz_sng
     tz_sng='GMT'
     call unix2gmt(time_unix,gmt_yr,gmt_ydy,gmt_mth,gmt_day,gmt_hr,gmt_mnt,gmt_sec)
     call unix2gmt_sng(time_unix,gmt_sng,tz_sng)
  */
  long time_unix_long; /* [s] UNIX time stored as intrinsic long */

  struct tm *gmt_tm_ptr; /* [tm] UNIX time structure */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  /* Initialize time structure */ 
  time_unix_long=(long)floor(*time_unix); /* [s] UNIX time stored as intrinsic long */
  time_unix_time_t=(time_t)time_unix_long; /* [s] UNIX time stored as time_t */

  /* Calling putenv() in separate subroutine with automatic string argument 
     may cause gmtime() to fail here, depending on putenv() implementation */
  gmt_tm_ptr=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */
  
  *sec=gmt_tm_ptr->tm_sec; /* [s] 0-based second of minute [0..59] */
  *mnt=gmt_tm_ptr->tm_min; /* [mnt] 0-based minute of second [0..59] */
  *hr=gmt_tm_ptr->tm_hour; /* [hr] 0-based, hour of day [0..23] */
  *day=gmt_tm_ptr->tm_mday; /* [day] 1-based day of month [1..31] */
  *mth=gmt_tm_ptr->tm_mon+1; /* [mth] 1-based month of year [1..12] */
  *ydy=gmt_tm_ptr->tm_yday+1; /* [day] 1-based day of year [1..366] */
  *yr=gmt_tm_ptr->tm_year+1900; /* [yr] Christian Year */

  return;
} /* end unix2gmt() */ 

void 
FNC_C2F(gmtime_tst)
     (double *time_unix) /* I [s] Pointer to UNIX time stored in double precision */
{
  /* [fnc] Test C library gmtime() routine from Fortran */

  long time_unix_long; /* [s] UNIX time stored as intrinsic long */

  struct tm *gmt_tm_ptr; /* [tm] UNIX time structure */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  /* Initialize time structure */ 
  time_unix_long=(long)floor(*time_unix); /* [s] UNIX time stored as intrinsic long */
  time_unix_time_t=(time_t)time_unix_long; /* [s] UNIX time stored as time_t */

  /* Calling putenv() in separate subroutine with automatic string argument 
     may cause gmtime() to fail here, depending on putenv() implementation */
  gmt_tm_ptr=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */

  gmt_tm_ptr=gmt_tm_ptr; /* CEWI: Avoids set but not used message */
  
  return;
} /* end gmtime_tst() */ 

void 
 /* [fnc] Create standard date-time string for given decimal time and time-zone */
FNC_C2F(unix2gmt_sng)
     (double *time_unix, /* I [s] Pointer to UNIX time_t */
      char *bfr, /* O [sng] Character buffer to hold output date-time string */
      char *tz_nm, /* I [sng] String to use for time-zone in output date-time string */
      int bfr_lng, /* I [nbr] Length of bfr string (passed automatically by Fortran-C convention) */
      int tz_lng) /* I [nbr] Length of tz_nm string (passed automatically by Fortran-C convention) */
{
  /* Determine standard date-time string corresponding to given decimal time and time-zone
     e.g., Fri Sep 16 13:11:08 1994 GMT 
     Fortran calling conventions are assumed, arguments are passed as pointers, strings have buffer lengths tacked onto end of argument list
     Usage:
     character*4 tz_sng
     character*32 gmt_sng
     tz_sng='GMT'
     double precision time_unix
     time_unix=time_unix_yr_srt+(gmt_doy-1.0)*86400.0
     call unix2gmt_sng(time_unix,gmt_sng,tz_sng)
  */

  char *tz_sys=NULL; /* [sng] Pointer to TZ variable in environment */
  char *tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */
  char *tz_set_cmd=NULL; /* [sng] Command for putenv() to set TZ */

#if (defined AIX) || (defined ALPHA)
  extern long timezone; /* [mnt] Offset of current TZ from GMT */
#else /* not AIX */ 
  extern time_t timezone; /* [mnt] Offset of current TZ from GMT */
#endif /* not AIX */ 
#if !defined(LINUX) && !defined(LINUXAMD64) && !defined(ALPHA)
  extern time_t altzone; /* [mnt] Offset of alternate TZ from GMT */
#endif /* not LINUX or ALPHA */ 
  extern int daylight; /* [flg] Daylight savings time */
  extern char *tzname[2]; /* [sng] UNIX code of main and alternate TZs */

  double time_unix_dbl; /* [s] UNIX time stored as double */
  double sec_frc; /* [s] Fractional seconds */

  int rcd=0; /* [rcd] Return code */

  long time_unix_long; /* [s] UNIX time stored as intrinsic long */
  long tz_lng_max=4; /* [nbr] Maximum length of TZ string (including trailing NUL) */

  static char *day_nm[]={"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};

  static char *mth_nm[]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

  struct tm *gmt_tm_ptr; /* [tm] UNIX time structure */

  time_t time_unix_time_t; /* [s] UNIX time stored as time_t */

  /* Guard against overflow */
  if(bfr_lng < 32){
    (void)fprintf(stderr,"WARNING time buffer too small in unix2gmt_sng()\n");
    return;
  } /* endif */ 

  if(tz_lng > tz_lng_max){
    (void)fprintf(stderr,"WARNING tz_lng = %d in unix2gmt_sng()\n",tz_lng);
    return;
  } /* endif */ 

  if(tz_lng < tz_lng_max){
    (void)fprintf(stderr,"WARNING tz_lng = %d in unix2gmt_sng()\n",tz_lng);
    return;
  } /* endif */ 

  /* NUL-terminate time-zone string */ 
  tz_nm[3]='\0';

  /* Initialize time structure */ 
  time_unix_dbl=*time_unix; /* [s] UNIX time stored as double */
  time_unix_long=(long)floor(time_unix_dbl); /* [s] UNIX time stored as intrinsic long */
  sec_frc=time_unix_dbl-time_unix_long; /* [s] Fractional seconds */
  time_unix_time_t=(time_t)time_unix_long; /* [s] UNIX time stored as time_t */

  /* NB: LINUX gmtime() assumes argument is UNIX time and returns the appropriate time sct */ 
  /* NB: SUNMP gmtime() assumes argument is local time and returns time sct accordingly
     To circumvent this feature, we reset timezone to GMT internally before calling gmtime() */ 
  /* Save original timezone environment */ 
  tz_lcl=(char *)strdup(getenv("TZ")); /* [sng] Local copy of TZ environment variable */

  /* Change to GMT */ 
  if(tz_sys != NULL){
    char tz_gmt_set_cmd[]="TZ=GMT"; /* [sng] Command to copy to use in putenv() */
    /* Allocate space for command and terminating NUL */
    tz_set_cmd=(char *)malloc((strlen(tz_gmt_set_cmd)+1)*sizeof(char)); /* [sng] Command for putenv() to set TZ */
    (void)sprintf(tz_set_cmd,"%s",tz_gmt_set_cmd); /* [sng] Command for putenv() to set TZ */
    rcd+=putenv(tz_set_cmd);
    (void)tzset();
    /* Do not free tz_set_cmd as as per caveats for Linux putenv() */
    /*    if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;*/
  } /* endif TZ exists */

  /* Calling putenv() in separate subroutine with automatic string argument 
     may cause gmtime() to fail here, depending on putenv() implementation */
  gmt_tm_ptr=gmtime(&time_unix_time_t); /* [tm] UNIX time structure */

  /* Restore original timezone environment */ 
  if(tz_sys != NULL){
    /* Allocate space for "TZ=" + tz_lcl + terminating NUL */
    tz_set_cmd=(char *)malloc((strlen(tz_lcl)+4)*sizeof(char)); /* [sng] Command for putenv() to set TZ */
    (void)sprintf(tz_set_cmd,"TZ=%s",tz_lcl); /* [sng] Command for putenv() to set TZ */
    rcd+=putenv(tz_set_cmd);
    (void)tzset();
    /* Do not free tz_set_cmd as as per caveats for Linux putenv() */
    /*    if(tz_set_cmd != NULL) (void)free(tz_set_cmd); tz_set_cmd=NULL;*/
  } /* endif TZ exists */

  if(tz_lcl != NULL) (void)free(tz_lcl);
  tz_lcl=NULL; /* [sng] Local copy of TZ environment variable */

  /* Standard mnemonics are not used so comment out 
  int mth; // [mth] 1-based month of year [1..12]
  int day; // [day] 1-based day of month [1..31]
  int ydy; // [day] 1-based day of year [1..366]
  int yr; // [yr] Christian Year
  int hr; // [hr] 0-based, hour of day [0..23]
  int mnt; // [mnt] 0-based minute of second [0..59]
  int sec; // [s] 0-based second of minute [0..59]
  mth=gmt_tm_ptr->tm_mon+1; // [mth] 1-based month of year [1..12]
  day=gmt_tm_ptr->tm_mday; // [day] 1-based day of month [1..31]
  ydy=gmt_tm_ptr->tm_yday+1; // [day] 1-based day of year [1..366]
  yr=gmt_tm_ptr->tm_year+1900; // [yr] Christian Year
  hr=gmt_tm_ptr->tm_hour; // [hr] 0-based, hour of day [0..23]
  mnt=gmt_tm_ptr->tm_min; // [mnt] 0-based minute of second [0..59]
  sec=gmt_tm_ptr->tm_sec; // [s] 0-based second of minute [0..59]
 */

  /*  (void)sprintf(bfr,"%4.4d/%2.2d/%2.2d YDY %d GMT %2.2d:%2.2d:%f\n",yr,mth,day,ydy,hr,mnt,sec+sec_frc);*/
  /* Currently 31 characters + NULL terminator*/
  (void)sprintf(bfr,"%3s %3s %2.2d %2.2d:%2.2d:%5.2f %4.4d %3s",
		day_nm[gmt_tm_ptr->tm_wday],
		mth_nm[gmt_tm_ptr->tm_mon],
		gmt_tm_ptr->tm_mday,
		gmt_tm_ptr->tm_hour,
		gmt_tm_ptr->tm_min,
		gmt_tm_ptr->tm_sec+sec_frc,
		gmt_tm_ptr->tm_year+1900,
		tz_nm
		);
  
  if(0){
    (void)fprintf(stderr,"Main Time Zone = MTZ = tzname[0] = %s\n",tzname[0]);
    (void)fprintf(stderr,"Alternate Time Zone = ATZ = tzname[1] = %s\n",tzname[1]);
    (void)fprintf(stderr,"UTC - MTZ = timezone = %li\n",(long)timezone);
#if !defined(LINUX) && !defined(LINUXAMD64) && !defined(ALPHA)
    (void)fprintf(stderr,"UTC - ATZ = altzone = %li\n",(long)altzone);
#endif /* not LINUX or ALPHA */ 
    (void)fprintf(stderr,"daylight = %d\n",daylight);
  } /* endif */   

  return;
} /* end unix2gmt_sng() */ 
