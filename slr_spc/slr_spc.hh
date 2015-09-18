// $Id$ 

// Typedefs and global variables for C++  

#ifndef _SLR_SPC_H // This include file has not yet been defined in the current source file  
#define _SLR_SPC_H

#include <string> // Standard C++ string class

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif // SUN4

#ifdef MAIN_PROGRAM_FILE // The current file contains main() 

// Global variables and variables with scope limited to the main.c file allocated here  

extern "C" {
int prg;
int prg_get(void){return prg;}

char *prg_nm;
char *prg_nm_get(void){return prg_nm;}

unsigned short dbg_lvl=0; // Option D 
unsigned short dbg_lvl_get(void){return dbg_lvl;}
} // end extern C

#else // MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main()  

// External references to global variables are declared as extern here.
// Variables with local file scope in all files except the main.c file are allocated here.  

#endif // MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main() 

#ifdef CRAY
#define FORTRAN_slfftd SLFFTD
#define FORTRAN_slffln SLFFLN
#endif
#if ( defined SUN4 ) || ( defined SUN4SOL2 ) || ( defined SUNMP ) || ( defined SGI5 ) || ( defined SGI64 ) || ( defined SGIMP64 )
#define FORTRAN_slfftd slfftd_
#define FORTRAN_slffln slffln_
#endif
// NB: g77 subroutines have two underscores by default, but g77 functions have only one!
// Use switches to change
#ifdef LINUX
#define FORTRAN_slfftd slfftd_
#define FORTRAN_slffln slffln_
#endif
#ifdef RS6K
#define FORTRAN_slfftd slfftd
#define FORTRAN_slffln slffln
#endif

// Declare functions that have FORTRAN linkages (g++ does not recognize extern "FORTRAN" yet)
extern "C" float FORTRAN_slfftd(const float *wvl_min_mcr,const float *wvl_max_mcr);
extern "C" float FORTRAN_slffln(const float *wvl_min_mcr,const float *wvl_max_mcr);

// Declare functions that have C linkages
extern "C" {
extern char *prg_nm_get(void);
extern int prg_get(void);
extern unsigned short dbg_lvl_get(void);
} // end extern C

// Declare functions that have C++ linkages
void Exit_gracefully(void);
void usg_prn(const char *opt_sng);

#endif // _SLR_SPC_H  






