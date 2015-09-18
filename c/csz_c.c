/* $Id$ */

/* Purpose: Standalone utilities for C programs (no netCDF required) */

/* Standard header files */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */
#include <string.h>             /* strcmp. . . */
#include <sys/stat.h>           /* stat() */
#include <time.h>               /* machine time */
#include <unistd.h>             /* all sorts of POSIX stuff */
/* #include <errno.h> */        /* errno */
/* #include <malloc.h>    */    /* malloc() stuff */
/* #include <assert.h> */       /* assert() debugging macro */

#include <sys/types.h>          /* needed for _res */
#include <netinet/in.h>         /* needed for _res */
#include <pwd.h>                /* password structures for getpwuid() */
#include <arpa/nameser.h>       /* needed for _res */
#include <resolv.h>             /* Internet structures for _res */

/* I'm only keeping these netCDF include files around because i'm worried that the
   function prototypes in nc.h are needed here. Eventually the prototypes for these
   routines should be broken into separate files, like csz.h... */
/*#include <netcdf.h>*/        /* netCDF def'ns */
/*#include "nc.h"*/            /* netCDF operator universal def'ns */

#ifndef bool
#define bool int
#endif /* bool */
#ifndef True
#define True 1
#endif /* True */
#ifndef False
#define False 0
#endif /* False */

#ifdef MAIN_PROGRAM_FILE /* Current file contains main() */

/* Global variables and variables with scope limited to main.c allocated here */

int prg; // [enm] Program ID
int prg_get(void){return prg;} // [enm] Program ID

char *prg_nm; // [sng] Program name
char *prg_nm_get(void){return prg_nm;} // [sng] Program name

unsigned short dbg_lvl=0; // [enm] Debugging level
unsigned short dbg_lvl_get(void){return dbg_lvl;} // [enm] Debugging level

#else /* MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main() */

/* External references to global variables are declared as extern here
   Variables with local file scope in all files except the main.c file are allocated here */

#endif /* MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main() */

void 
#ifdef CRAY
DAT_PRN_F77
#else /* not CRAY */
dat_prn_f77_
#endif /* not CRAY */
(int *dat_nbr,float *dat,int *nbr_dat_per_ln,char *dat_sng,int sng_len)
{
  /* Routine prints input data out to stderr in Fortran 77 block data format */

  long idx;

  (void)fprintf(stderr,"      data %s /",dat_sng);
  for(idx=0;idx<*dat_nbr;idx++){
    if(!(idx % *nbr_dat_per_ln)) (void)fprintf(stderr,"\n     $     ");
    (void)fprintf(stderr,"%.7e",dat[idx]);
    if(idx != *dat_nbr-1) (void)fprintf(stderr,", ");
  } /* end loop over wvl */
  (void)fprintf(stderr," /\n\n");
  
} /* end dat_prn_f77() */
