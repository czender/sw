// $Id$ 

// Purpose: Description (definition) of classes

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <cls.hh> // Program-specific class definitions

#ifndef CLS_HH // This include file has not yet been defined in current source file
#define CLS_HH

// Standard C++ headers
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 

// Personal headers
#include <ccc.hh> // Program-specific definitions
#include <libcsz_c++.hh> // Personal C++ library

// Namespaces

// Typedefs
typedef std::map<std::string,prc_cmp,std::less<std::string> > sng2cpv; // String-to-computational precision map
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Forward class declarations
class Base;
class SzDstFnc;
class Sine;
class Cosine;

////////////////////////////////////////////////////////////////////////

// Define Base class
class Base{
public:
  // Static public members
  // Public member functions
  Base(SzDstFnc *szdstfnc); // Default constructor
  void set_fnc(SzDstFnc *szdstfnc); // Set function type
  SzDstFnc *fnc;
  //  prc_cmp eval()const; // Evaluate fnc
private:
  // Static private members
  // Private members
  //prc_cmp foo;
}; // end class Base

// Define Test class
class Test{
public:
  // Static public members
  // Public member functions
  // Note that initializing arguments to default constructor must use "=" and may not use instantiation ()
  Test(long tst_nbr=1); // Default constructor
  //  Test(); // Default constructor (ambiguous now that other constructor has complete set of default arguments, making it the default constructor)
  long tst_nbr_get()const; // Get tst_nbr
  //  prc_cmp eval()const; // Evaluate fnc
private:
  // Static private members
  // Private members
  long tst_nbr;
}; // end class Test

// Define SzDstFnc class
class SzDstFnc{
public:
  // Static public members
  static std::string opt2abb(const std::string opt_sng); // Create option to abbreviation map
  // Public member functions
  SzDstFnc(const std::string sng); // [fnc] Default constructor
  virtual ~SzDstFnc(); // [fnc] Destructor
  virtual prc_cmp operator()(const prc_cmp &rds)const=0; // Get # between r and r+dr per unit r
  virtual prc_cmp eval(const prc_cmp &rds)const=0; // Get # between r and r+dr per unit r
private:
  // Static private members
  static sng2sng_map opt_map_mk(); // Create abbreviation map
  static sng2sng_map opt_map; // Option to abbreviation map
  // Private members
  std::string fnc_sng;
  //prc_cmp foo;
}; // end class SzDstFnc

// Define Sine class
class Sine : public SzDstFnc{
public:
  // Static public members
  // Public member functions
  Sine(const std::string sng); // Default constructor
  prc_cmp operator()(const prc_cmp &ngl)const;
  prc_cmp eval(const prc_cmp &ngl)const;
private:
  // Static private members
  // Private members
  prc_cmp foo;
}; // end class Sine

// Define Cosine class
class Cosine : public SzDstFnc{
public:
  // Static public members
  // Public member functions
  Cosine(const std::string sng); // Default constructor
  prc_cmp operator()(const prc_cmp &ngl)const;
  prc_cmp eval(const prc_cmp &ngl)const;
private:
  // Static private members
  // Private members
  prc_cmp foo;
}; // end class Cosine

#endif // CLS_HH  








