// $Id$

// Purpose: Description (definition) of extrapolation classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage: 
// Variables in this header are set with routine xtr_ini() in libcsz_c++
// #include <xtr.hh> // Extrapolation class

#ifndef XTR_HH // Contents have not yet been inserted in current source file
#define XTR_HH

// Standard C++ headers
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Namespaces

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Forward declarations
class Xtr_cls;

// Define Xtr_cls class
class Xtr_cls{
public:
  // Static public members
  static std::string abb2dsc(const std::string &abb_sng); // Given abbreviation, return description
  // Public member functions
  Xtr_cls(const std::string &opt_sng); // Default constructor
  void flg_set(const std::string &xtr_sng); // Set flags

  std::string xtr_sng; // Concise string description of current configuation
  bool xtr_err; // No extrapolation is allowed
  bool xtr_fll; // Full extrapolation is allowed (implies xtr_prt)
  bool xtr_fll_lnr; // Perform linear extrapolation using two nearest valid neighbors
  bool xtr_fll_ngh; // Set extrapolated value to value of nearest valid neighbor
  bool xtr_fll_nil; // Set extrapolated value to 0.0
  bool xtr_prt; // Partial extrapolation is allowed
  bool xtr_prt_frc; // Use average value of overlap region
  bool xtr_prt_lnr; // Perform linear extrapolation using two nearest valid neighbors
  bool xtr_prt_ngh; // Set extrapolated value to value of nearest valid neighbor
  bool xtr_prt_nil; // Set extrapolated value to 0.0
  bool xtr_prt_wgt; // Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc)
  bool xtr_vrb; // Print verbose warning when extrapolation is performed
private:
  // Static private members
  static sng2sng_map abb2dsc_map_mk(); // Create abbreviation to description map
  static sng2sng_map abb2dsc_map; // Abbreviation to descritption map
  // Private members
  void flg_vld(); // Validate flags
  void flg_fls(); // Initialize all flags to false
}; // end class Xtr_cls

#endif // XTR_HH  
