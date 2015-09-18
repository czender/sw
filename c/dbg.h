/* $Id$ */

/* Purpose: Systematic nomenclature for debugging/verbosity levels */

/* Usage: */
/* #include <dbg.h> *//* Debugging constants */

#ifndef DBG_H /* This include file has not yet been defined in the current source file */
#define DBG_H

/*
  NB: dbg.h (and dbg.c) is the natural place to define the dbg_lvl_get functions
  extern "C" {
  unsigned short dbg_lvl=0; // Option D 
  unsigned short dbg_lvl_get(void){return dbg_lvl;} // end dbg_lvl_get()
  } 
*/
 
  /* Debugging levels are as follows: */
const int dbg_nbr=9; /* Number of different debugging levels */

const int dbg_off=0; /* Production mode. Debugging is turned off. */
const int dbg_fl=1; /* Filenames */
const int dbg_scl=2; /* Scalars */
const int dbg_crr=3; /* Current task */
const int dbg_sbr=4; /* Subroutine names on entry and exit */
const int dbg_io=5; /* Subroutine I/O */
const int dbg_vec=6; /* Entire vectors */
const int dbg_vrb=7; /* Everything */
const int dbg_old=8; /* Old debugging blocks not used anymore */

#endif /* DBG_H */
