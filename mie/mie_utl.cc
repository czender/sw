// $Id$ 

// Purpose: Implementation (declaration) of new utilities for Mie program

/* Copyright (C) 1997--2014 Charlie Zender
   You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Standard C++ header files 
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class
#include <complex> // Standard C++ complex class

// Standard C header files 
#include <cstring> // strcmp. . . 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdlib> // atof, atoi, malloc, getopt, strtod
#include <unistd.h> // All sorts of POSIX stuff  

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
#include <libcsz_c++.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <libnco_c++.hh> // C++ interface to netCDF C library

int foo_fnc(){return 1;} // [fnc] Dummy function to keep file from being empty
