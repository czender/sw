// $Id$ 

// Purpose: Typedefs and global variables for gsl

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#ifndef GSL_HH // Contents have not yet been inserted in current source file  
#define GSL_HH

// Standard C++ headers
#include <string> // Standard C++ string class
#include <complex> // Standard C++ complex class

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif // EXIT_SUCCESS is not defined in SUN4

#ifdef MAIN_PROGRAM_FILE // The current file contains main() 

// Global variables and variables with scope limited to the main.c file allocated here  

//std::string prg_nm;
//std::string prg_nm_get(void){return prg_nm;}
extern "C" {
  char *prg_nm;
  char *prg_nm_get(void){return prg_nm;}

  unsigned short dbg_lvl=0; // Option D 
  unsigned short dbg_lvl_get(void){return dbg_lvl;}
} // end extern C

#else // MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main()  

// External references to global variables are declared as extern here.
// Variables with local file scope in all files except the main.c file are allocated here.  

#endif // MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main() 

// Define global structures, enums, and functions
// Declare functions that have FORTRAN linkages (g++ does not recognize extern "FORTRAN" yet)
// Declare functions that have C linkages
extern "C" {
  extern char *prg_nm_get(void);
  extern int prg_get(void);
  extern unsigned short dbg_lvl_get(void);
  extern void *nc_lib_vrs_prn();
} // end extern C
//std::string prg_nm_get(void);

// Declare functions that have C++ linkages
//extern std::string prg_nm_get(void);

#endif // GSL_HH  






