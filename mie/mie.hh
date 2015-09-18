// $Id$ 

// Purpose: Typedefs and global variables for mie program

/* Copyright (C) 1997--2014 Charlie Zender
   You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie.hh> // Global variables and functions

#ifndef MIE_HH // Contents have not yet been inserted in current source file 
#define MIE_HH

// EXIT_SUCCESS is not defined in SUN4
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif // EXIT_SUCCESS

#ifdef MAIN_PROGRAM_FILE // Current file contains main() 

// Global variables and variables with scope limited to main.c allocated here

// Implementation (declaration) functions that have C linkages
extern "C" {
char *prg_nm; // [sng] Program name
char *prg_nm_get(void){return prg_nm;} // [sng] Program name

unsigned short int dbg_lvl=0; // [enm] Debugging level
unsigned short int dbg_lvl_get(void){return dbg_lvl;} // [enm] Debugging level
} // end extern C
  
#else // MAIN_PROGRAM_FILE is NOT defined, i.e., the current file does not contain main()  

// External references to global variables are declared as extern here
// Variables with local file scope in all files except the main.c file are allocated here 

#endif // MAIN_PROGRAM_FILE is NOT defined, i.e., current file does not contain main() 

// Typedefs

// Prototype functions that have C linkages
extern "C" {
  extern char *prg_nm_get(void); // [sng] Program name
  extern int prg_get(void); // [enm] Program ID
  extern unsigned short dbg_lvl_get(void); // [enm] Debugging level
} // end extern C

// Prototype functions that have C++ linkages
//void usg_prn(const char *opt_sng);

#endif // MIE_HH  




