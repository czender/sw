// $Id$ 

// Description (definition) of new classes for Mie program

/* Copyright (C) 1997--2014 Charlie Zender
   You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie_cls.hh> // Program-specific class definitions

#ifndef MIE_CLS_HH // Contents have not yet been inserted in current source file
#define MIE_CLS_HH

// Standard C++ header files 
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 

// 3rd party vendors
#include <netcdf.h> // netCDF C interface

// Personal headers
#include <libcsz_c++.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Typedefs

// Define Test class
class Test{
public:
  // Static public members
  // Public member functions
  Test(long tst_nbr=1); // [fnc] Default constructor
  long tst_nbr_get()const; // [nbr] Number of tst
private:
  // Static private members
  // Private members
  long tst_nbr; // [nbr] Number of tst
}; // end class Test

// Prototype functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

#endif // MIE_CLS_HH  
