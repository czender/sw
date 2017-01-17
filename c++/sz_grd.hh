// $Id$ 

// Purpose: Description (definition) of Size grid classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <sz_grd.hh> // Size grids

#ifndef SZ_GRD_HH // Contents have not yet been inserted in current source file
#define SZ_GRD_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cassert> // Assertions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Forward declarations
class sz_grd_cls; // Size grids

// Define sz_grd_cls class
class sz_grd_cls{
public:
  // Static public members
  static std::string opt2abb(const std::string opt_sng); // [fnc] Create option to abbreviation map

  // Public member functions
  sz_grd_cls // Default constructor
  (const prc_cmp sz_mnm_arg=0.9e-6, // [m] Minimum size in distribution
   const prc_cmp sz_mxm_arg=1.1e-6, // [m] Maximum size in distribution
   const long sz_nbr_arg=1L, // [nbr] Number of size bins
   const std::string grd_sng_arg="logarithmic"); // [sng] Type of size grid
  sz_grd_cls // Copy constructor
  (const sz_grd_cls &orig);
  ~sz_grd_cls(); // Destructor
  const sz_grd_cls &operator= // Assignment operator
  (const sz_grd_cls &RHS);

  void typ_set(const std::string sng); // [fnc] Type of size grid
  std::string typ_get()const; // [fnc] Type of size grid
  void nbr_set(const long sz_nbr_arg); // [fnc] Set number of elements in size grid
  long nbr_get()const; // [fnc] Get number of elements in size grid
  const prc_cmp *sz_ctr_get()const; // [fnc] Get pointer to const data
  const prc_cmp *sz_grd_get()const; // [fnc] Get pointer to const data
  const prc_cmp *sz_max_get()const; // [fnc] Get pointer to const data
  const prc_cmp *sz_min_get()const; // [fnc] Get pointer to const data
  const prc_cmp *sz_dlt_get()const; // [fnc] Get pointer to const data
private:
  // Private members
  int rcd; // [enm] Return success code
  std::string grd_sng;
  prc_cmp sz_mnm;
  prc_cmp sz_mxm;
  long sz_nbr; // if grd_nbr were const then it could not be dynamically resized
  prc_cmp *sz_ctr;
  prc_cmp *sz_min;
  prc_cmp *sz_max;
  prc_cmp *sz_dlt;
  prc_cmp *sz_grd;
  // Private member functions
  int allocate();
  int deallocate();
  int recompute();
  // Static private members
  static sng2sng_map opt2abb_map_mk(); // [fnc] Create abbreviation map
  static sng2sng_map opt2abb_map; // [fnc] Option to abbreviation map
}; // end class sz_grd_cls

// Declare functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

#endif // SZ_GRD_HH  






