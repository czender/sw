static char rcs_Id[] = "$Id$\n";
static char rcs_Revision[] = "$Revision$\n";

/* $Id$ */

/* Purpose:  Compute indices of refraction weighted over a given solar wavelength interval */  
/* Inputs: the beginning and ending of the wavelength interval in microns */ 
/* Outputs: the fractional flux, absolute flux, and indice of refraction */

/* Example Usage (place mousable command lines here):
   here are the default settings:
   idx_rfr -a 0 -l .5 -h 1.1 -s 1367. -t 230. -n 1 -D 0
   idx_rfr -a 0 -l .695 -h 5. -s 1367. -t 230. -n 1 -D 0
   idx_rfr -a 0 -l .695 -h 2.7 -s 1367. -t 230. -n 1 -D 0
 */ 

#define bool int
#define True 1
#define False 0
#define FL_SZ 80
#define CMD_LN_SIZE 200
#define METERS_PER_MICRON 1.e-6

/* First author of the solar flux data to employ */ 
enum slr_flx_author{
  THEKEAKARA, /* 0 = ThD71 */ 
  LABS, /* 1 = LaN68 */ 
  NECKEL}; /* 2 = NeL84 */ 

/* Standard header files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */ 
#include <unistd.h>             /* all sorts of POSIX stuff */ 
#if ( defined SUN4SOL2 ) || ( defined SUNMP ) 
#include "/opt/SUNWspro/SC4.0/include/cc/sunmath.h"            /* IEEE signal handling */
#include <siginfo.h>            /* IEEE stuff, je ne sais pas */ 
#include <ucontext.h>           /* IEEE stuff, je ne sais pas */ 
#endif
#if defined (SUN4)
#include <signal.h>            /* IEEE stuff, je ne sais pas */ 
#endif

#ifdef CRAY
#define FORTRAN_refice REFICE
#define FORTRAN_slfftd SLFFTD
#define FORTRAN_slffln SLFFLN
#endif
#if ( defined SUN4 ) || ( defined SUN4SOL2 ) || ( defined SUNMP ) || ( defined SGI5 )
#define FORTRAN_refice refice_
#define FORTRAN_slfftd slfftd_
#define FORTRAN_slffln slffln_
#endif
#ifdef RS6K
#define FORTRAN_refice refice
#define FORTRAN_slfftd slfftd
#define FORTRAN_slffln slffln
#endif

/* Global variables declared here */ 
int debug=0; /* Option D */
int debug_value=0; /* Option d */

int main(int argc,char **argv)
{
  float FORTRAN_slfftd(float */* LMIN */, float */* LMAX */ );
  float FORTRAN_slffln(float */* LMIN */, float */* LMAX */ );

  void FORTRAN_refice();

  extern char *ctime(const time_t *);  
  extern time_t time(time_t *);
  extern void exit(int);

  void Exit_gracefully(void);
  void prn_usg(char *);
  void sng_cmd_ln(int,char *[],char *,int);

  bool BATCH_JOB;
  bool WVN_INPUT;
  bool VERBOSE;

  char *time_buf_srt;
  char *opt_sng;

  char *slr_flx_sng[3]={
    "Thekeakara and Drummond (1971)",
    "Labs and Neckel (1968)",
    "Neckel and Labs (1984)"};

  char cmd_ln[CMD_LN_SIZE];
  char out_fl[FL_SZ];
  char in_fl[FL_SZ];
  char err_fl[FL_SZ];
  
  extern char *optarg;
  extern int optind;
  
  float abs_slr_flx;
  float bnd_ctr_mcr;
  float bnd_max_mcr;
  float bnd_min_mcr;
  float bnd_sz_mcr;
  float wvn_max;
  float wvn_min;
  float wvn_sz;
  float float_foo_input;
  float frc_slr_flx;
  float ice_temperature;
  float idx_rfr_img;
  float idx_rfr_rl;
  float slr_cst;
  float sub_bnd_ctr_mcr;
  float sub_bnd_max_mcr;
  float sub_bnd_min_mcr;
  float sub_bnd_sz_mcr;
  float ttl_frc_slr_flx;
  float ttl_idx_rfr_img;
  float ttl_idx_rfr_rl;
  float ttl_wvl_mcr;
  float wvl_mcr;

  int int_foo;
  int sub_bnd;
  int nbr_sub_bnd;
  int opt;
  int slr_flx_author;

  time_t clock;

  /* set defaults */
  BATCH_JOB = False; /* Option B */ 
  WVN_INPUT = False; /* Option w */ 

  float_foo_input = 0.; /* Option f */ 
  ice_temperature = 230.; /* Option t */ 
  nbr_sub_bnd = 1; /* Option n */ 
  slr_cst = 1367.; /* Option s */ 
  slr_flx_author = THEKEAKARA; /* Option a */ 
  bnd_max_mcr = 1.1; /* Option h */ 
  bnd_min_mcr = 1.0; /* Option l */ 

  (void)strcpy(in_fl,"stdin"); /* Option i */
  (void)strcpy(out_fl,"stdout"); /* Option o */
  (void)strcpy(err_fl,"stderr"); /* Option e */
  
  /* parse command line arguments */
  opt_sng="a:BD:d:e:f:h:i:l:n:o:s:t:Vvw";
  while((opt = getopt(argc,argv,opt_sng)) != EOF){
    switch(opt){
    case 'a':
      /* The solar flux study to use.  Default is THEKEAKARA. */
      slr_flx_author = atoi(optarg);
      break;
    case 'B':
      /* toggle printing to stdout of batch output. Default is False */ 
      BATCH_JOB=!BATCH_JOB;
      break;
    case 'D':
      /* The debugging level.  Default is 0. */
      debug = atoi(optarg);
      break;
    case 'd':
      /* The second debugging level.  Default is 0. */
      debug_value = atoi(optarg);
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_fl,optarg);
      break;
    case 'f':
      /* Set the generic tuning parameter.  Default is 0. */
      float_foo_input = atof(optarg);
      break;
    case 'h':
      /* maximum wavelength in interval, in microns. Default is 1.1 */
      bnd_max_mcr = atof(optarg);
      break;
    case 'i':
      /* get the input file name. Default is stdin */
      (void)strcpy(in_fl,optarg);
      break;
    case 'l':
      /* minimum wavelength in interval, in microns. Default is 1.0 */
      bnd_min_mcr = atof(optarg);
      break;
    case 'n':
      /* number of sub-intervals to divide band into.  Default is 1 */
      nbr_sub_bnd = atoi(optarg);
      break;
    case 'o':
      /* get the output file name. Default is stdout */
      (void)strcpy(out_fl,optarg);
      break;
    case 's':
      /* integrated TOA solar flux in W/m^2. Default is 1367. */
      slr_cst = atof(optarg);
      break;
    case 't':
      /* Set the ice temperature.  Default is 230. */
      ice_temperature = atof(optarg);
      break;
    case 'v':
      /* print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: idx_rfr.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/cld/idx_rfr.c,v $\n");
      (void)fprintf(stderr,"$Id$\n");
      (void)fprintf(stderr,"$State: Exp $\n");
      exit(0);
      break;
    case 'V':
      /* toggle verbose printing out of WARNINGS. Default is True */ 
      VERBOSE=!VERBOSE;
      break;
    case 'w':
      /* wavenumber input */ 
      WVN_INPUT=!WVN_INPUT;
      break;
    case '?':
      /* print proper usage */
      (void)prn_usg(opt_sng);
      exit(1);
    } /* end switch */
  } /* end while loop */
  
  /* start the clock and save the command line */  
  sng_cmd_ln(argc,argv,cmd_ln,CMD_LN_SIZE);
  (void)fprintf(stdout,"Command Line: %s\n",cmd_ln);
  clock=time((time_t *)NULL);
  time_buf_srt=ctime(&clock);
  (void)fprintf(stdout,"\tstart = %s",time_buf_srt);

  if((bnd_min_mcr == bnd_max_mcr) && (nbr_sub_bnd >= 2)){
    (void)fprintf(stdout,"Only use one band if you're not weighting the index of refraction over a finite interval.\nProceeding anyway...\n\n");
  } /* end if */

  if(WVN_INPUT){
    /* in this case -h specified the maximum wavenumber, and -l specified the minimum wavenumber */ 
    wvn_max=bnd_max_mcr;
    wvn_min=bnd_min_mcr;
    wvn_sz=wvn_max-wvn_min;

    bnd_max_mcr=1.e6/(wvn_min*100.);
    bnd_min_mcr=1.e6/(wvn_max*100.);
  }else{
    wvn_max=1.e6/(bnd_min_mcr*100.);
    wvn_min=1.e6/(bnd_max_mcr*100.);
    wvn_sz=wvn_max-wvn_min;
  } /* end if */

  bnd_sz_mcr=bnd_max_mcr-bnd_min_mcr;
  bnd_ctr_mcr=.5*(bnd_max_mcr+bnd_min_mcr);
  sub_bnd_sz_mcr=(bnd_max_mcr-bnd_min_mcr)/nbr_sub_bnd;

  if(True){
    (void)fprintf
      (stdout,"Initialization state:\n");
    (void)fprintf(stdout,"source of spectral solar flux = %s\n",
		  slr_flx_sng[slr_flx_author]);
    (void)fprintf(stdout,"solar constant = %10.1f W/m^2\n",slr_cst);
    (void)fprintf(stdout,"ice temperature = %10.1f K\n",ice_temperature);
    (void)fprintf(stdout,"minimum wavelength = %10.2f microns = %g cm-1\n",
		  bnd_min_mcr,wvn_max);
    (void)fprintf(stdout,"maximum wavelength = %10.2f microns = %g cm-1\n",
		  bnd_max_mcr,wvn_min);
    (void)fprintf(stdout,"midpoint wavelength = %10.2f microns\n",
		  bnd_ctr_mcr);
    (void)fprintf(stdout,"number of sub-bands = %i\n",nbr_sub_bnd);
    (void)fprintf(stdout,"width of each sub-band = %10.3f microns\n",
		  sub_bnd_sz_mcr);
    (void)fprintf(stdout,"\n");
  } /* end if */

  /* initialize the accumulating variables */ 
  ttl_idx_rfr_img=0.;
  ttl_idx_rfr_rl=0.;
  ttl_wvl_mcr=0.;
  ttl_frc_slr_flx=0.;
  sub_bnd_min_mcr=bnd_min_mcr;
  sub_bnd_max_mcr=bnd_min_mcr+sub_bnd_sz_mcr;
  sub_bnd_ctr_mcr=.5*(sub_bnd_max_mcr+sub_bnd_min_mcr);

  /* title the current sub-interval's properties */ 
  if(debug == 1){
    (void)fprintf
      (stdout,"Properties of the %i sub-intervals:\n",nbr_sub_bnd);
    (void)fprintf
      (stdout,"%10s%10s%10s%10s%10s%10s%10s%10s\n",
       "sub-band","lambd min","lambd max","lambd ctr","real indx",
       "imag indx","frac flux","abs flux");
    (void)fprintf
      (stdout,"%10s%10s%10s%10s%10s%10s%10s%10s\n",
       "index","microns","microns","microns","#",
       "#","%","W/m^2");
  } /* end debug */
  
  for(sub_bnd=1;sub_bnd<=nbr_sub_bnd;sub_bnd++){
    if(slr_flx_author == THEKEAKARA){
      frc_slr_flx=
	FORTRAN_slfftd
	  (&sub_bnd_min_mcr,  /* LMIN */
	   &sub_bnd_max_mcr);  /* LMAX */
    }else if(slr_flx_author == LABS){
      frc_slr_flx=
	FORTRAN_slffln
	  (&sub_bnd_min_mcr,  /* LMIN */
	   &sub_bnd_max_mcr);  /* LMAX */
    }else if(slr_flx_author == NECKEL){
      (void)fprintf(stdout,"Neckel & Labs option unavailable at this time.\n");
      Exit_gracefully();
      exit(1);
    } /* end else */
    
    FORTRAN_refice
      (&sub_bnd_ctr_mcr,&ice_temperature,
       &idx_rfr_rl,&idx_rfr_img);
    
    abs_slr_flx=slr_cst*frc_slr_flx;

    /* report the current sub-interval's properties */ 
    if(debug == 1){
      (void)fprintf(stdout,"%10i%10.3f%10.3f%10.3f%10.7f%10.3e%10.2f%10.4f\n",
		    sub_bnd,
		    sub_bnd_min_mcr,
		    sub_bnd_max_mcr,
		    sub_bnd_ctr_mcr,
		    idx_rfr_rl,
		    idx_rfr_img,
		    frc_slr_flx*100.,
		    abs_slr_flx);
    } /* end debug */

    /* add the weighted index of refraction to the accumulative totals */ 
    ttl_idx_rfr_img+=idx_rfr_img*frc_slr_flx;
    ttl_idx_rfr_rl+=idx_rfr_rl*frc_slr_flx;
    ttl_wvl_mcr+=sub_bnd_ctr_mcr*frc_slr_flx;
    ttl_frc_slr_flx+=frc_slr_flx;
    sub_bnd_max_mcr+=sub_bnd_sz_mcr;
    sub_bnd_ctr_mcr+=sub_bnd_sz_mcr;
    sub_bnd_min_mcr+=sub_bnd_sz_mcr;
  } /* end loop over sub_bnd */

  /* normalize by the total flux to find the weighted average */ 
  /* in the non-averaging case, just report the index at this wavelength,
     which is already contained in idx_rfr_img, idx_rfr_rl. the fluxes should
     be 0. */ 
  if(bnd_min_mcr != bnd_max_mcr){
    /* proceed normally */ 
    idx_rfr_img=ttl_idx_rfr_img/ttl_frc_slr_flx;
    idx_rfr_rl=ttl_idx_rfr_rl/ttl_frc_slr_flx;
    wvl_mcr=ttl_wvl_mcr/ttl_frc_slr_flx;
  }else{
    wvl_mcr=bnd_ctr_mcr;
  } /* endelse */ 
  abs_slr_flx=slr_cst*ttl_frc_slr_flx;

  if(True){
    (void)fprintf(stdout,"\n");
    (void)fprintf(stdout,"Computed properties of the entire band:\n");
    (void)fprintf(stdout,"flux weighted real index = %10.7f\n",idx_rfr_rl);
    (void)fprintf(stdout,"flux weighted imag index = %10.7f = %7.3e\n",
		  idx_rfr_img,idx_rfr_img);
    (void)fprintf(stdout,"flux weighted wavelength = %10.3f microns\n",
		  wvl_mcr);
    (void)fprintf(stdout,"fraction TOA solar flux in interval = %10.2f %%\n",
		  ttl_frc_slr_flx*100.);
    (void)fprintf(stdout,"absolute TOA solar flux in interval = %10.2f W/m^2\n",
		  abs_slr_flx);
    (void)fprintf(stdout,"TOA solar spectral flux in interval = %g W/m^2/micron\n",
		  abs_slr_flx/bnd_sz_mcr);
    (void)fprintf(stdout,"TOA solar spectral flux in interval = %g W/m^2/(cm-1)\n",
		  abs_slr_flx/wvn_sz);
  } /* end if */

  if(BATCH_JOB){
    (void)fprintf(stdout,"flux weighted real index = %10.7f = %7.3e\n",
		  idx_rfr_rl,idx_rfr_rl);
    (void)fprintf(stdout,"flux weighted imag index = %10.7f = %7.3e\n",
		  idx_rfr_img,idx_rfr_img);
    (void)fprintf(stdout,"flux weighted wavelength = %10.3f = %7.3e microns\n",
		  wvl_mcr,wvl_mcr);
  } /* end if */

  Exit_gracefully();
} /* end main() */

void prn_usg(char *opt_sng)
{
  (void)fprintf(stderr,"\nusage: cloud [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",opt_sng);
  (void)fprintf(stderr,"-B toggle batch output to stdout. Default is False\n");
  (void)fprintf(stderr,"-D debug The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"-a slr_flx_author  Default is 0 (THEKEAKARA)\n");
  (void)fprintf(stderr,"-d debug_value The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-e err_fl Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-f float_foo_input Set the generic float. Default is 0.\n");
  (void)fprintf(stderr,"-h bnd_max_mcr Default is 1.1\n");
  (void)fprintf(stderr,"-i in_fl get the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-l bnd_min_mcr  Default is 1.0\n");
  (void)fprintf(stderr,"-n nbr_sub_bnd # of sub-intervals Default is 1\n");
  (void)fprintf(stderr,"-o out_fl get the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-s slr_cst Default is 1367.\n");
  (void)fprintf(stderr,"-t ice_temperature Default is 230.\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"\n");
} /* end prn_usg() */ 

void Exit_gracefully(void)
{
  extern char *ctime(const time_t *);  
  extern time_t time(time_t *);

  char *time_buf_finish;
  time_t clock;

  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stdout,"\tfinish = %s\n",time_buf_finish);

  exit(0);
} /* end Exit_gracefully() */ 

/* Module:	cmdparse.c (Command Parse)
 * Purpose:	Set options from command line
 * Subroutine:	parse_cmd_ln()			returns: int
 * Subroutine:	sng_cmd_ln()		returns: void
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
 *              {2} CZ ANSI C compatability                       29 May 1993
 *		{n} <who> -- <does what> -- <when>
 */

void sng_cmd_ln(int argc,char *argv[],char *cmd_ln,int linemax)
{
  int current_arg;

  if(argc <= 0){
    cmd_ln[0]='\0';
  }else{
    (void)strcpy(cmd_ln, argv[0]);
    for(current_arg=1;current_arg<argc;current_arg++){
      (void)strncat(cmd_ln," ",linemax);
      (void)strncat(cmd_ln,argv[current_arg],linemax);
    } /* end loop over args */ 
  } /* end else */ 
} /* end sng_cmd_ln() */ 
