// $Id$ 

// Purpose: Chemical utilities

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <chm.hh> // Chemical utilities

#ifndef CHM_HH // Contents have not yet been inserted in current source file  
#define CHM_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <phys_cst.hh> // Physical constants NB: Live code, not a header

// Chemical structure
typedef struct{
  std::string abb; // [sng] Standard abbreviation
  std::string dsc; // [sng] Description
  double mmw; // [kg mol-1] Mean molecular weight
  prc_cmp cff_hnr_H2O_298K; // [mol ltr-1 atm-1] Henry's Law coefficient in liquid water at 298 K
  prc_cmp rxn_ntp_H2O_298K; // [J mol-1] Reaction enthalpy at 298K
  prc_cmp rct_nrm; // [frc] Normalized reactivity for dry deposition SeP97 p. 973 (19.25)
  prc_cmp cff_hnr_H2O_298K_ffc_dps; // [mol ltr-1 atm-1] Effective Henry's Law coefficient in liquid water at 298 K at pH 6.5
} chm_sct;

// Typedefs
typedef std::map<std::string,chm_sct,std::less<std::string> > sng2chm_sct_map; // String-to-chm_sct map
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Forward class declarations
class Chm;

////////////////////////////////////////////////////////////////////////

class Chm{
public:
  // Static class members
  static int nst_nbr_get(); // Get current number of instantiated Chm objects
  static std::string opt2abb(const std::string opt_sng); // Convert option to abbreviation 
  static std::string abb2dsc(const std::string opt_sng); // Convert abbreviation to description
  static prc_cmp abb2mmw(const std::string opt_sng); // Convert abbreviation to mmw

  // Public member functions
  Chm(const std::string chm_sng="HNO3"); // Default constructor
  Chm(const std::string chm_sng,const prc_cmp mmw_chm); // Overloaded constructor
  ~Chm(); // Destructor
  prc_cmp mmw_get()const; // Get const data
  void typ(const std::string sng); // Set chemical abbreviation
  std::string abb_get()const; // Get chemical abbreviation
  std::string dsc_get()const; // Get chemical description
private:
  // Private members
  prc_cmp mmw; // Mean molecular weight
  std::string abb; // Abbreviation
  std::string dsc; // Description
  static int nst_nbr; // Current number of instantiated Chm objects
  // Private member functions
  void recompute();
  // Static private class members
  static sng2sng_map opt2abb_map_mk(); // Create option to abbreviation map
  static sng2sng_map opt2abb_map; // Option to abbreviation map
  static const sng2chm_sct_map chm_map_mk(); // Create chemical map
  static const sng2chm_sct_map chm_map; // Chemical map
}; // end class Chm

////////////////////////////////////////////////////////////////////////

// Declare functions with C++ linkages
int // O [rcd] Return success code
cnc_nbr_get // [fnc] Convert species mass mixing ratio and air density to species concentration
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnc_nbr, // O [mlc m-3] Number concentration of species
 const prc_cmp *dns, // I [kg m-3] Midlayer density
 const prc_cmp mmw, // I [kg mol-1] Mean molecular weight of species
 const prc_cmp *qmmr); // I [kg kg-1] Mass mixing ratio of species

int // O [rcd] Return success code
cnc_nbr_get // [fnc] Convert species mass mixing ratio and air pressure and temperature to species concentration
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *cnc_nbr, // O [mlc m-3] Number concentration of species
 const prc_cmp *prs, // I [Pa] Pressure
 const prc_cmp mmw, // I [kg mol-1] Mean molecular weight of species
 const prc_cmp *tpt, // I [K] Temperature
 const prc_cmp *qmmr); // I [kg kg-1] Mass mixing ratio of species

// Define inline'd functions in header so source is visible to calling files

#endif // CHM_HH  






