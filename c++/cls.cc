// $Id$ 

// Purpose: Implementation (declaration) of classes

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include "cls.hh" // Class-specific definitions
////////////////////////////////////////////////////////////////////////

// Base class

// Non-static members
Base::Base(SzDstFnc *szdstfnc) // Overload constructor
{
  fnc=szdstfnc;
} // end Default constructor
void Base::set_fnc(SzDstFnc *szdstfnc){fnc=szdstfnc;} // [fnc] Set function type

// Initialize static members

////////////////////////////////////////////////////////////////////////

// Test class

// Non-static members
// Test::Test:tst_nbr(1){;} // Overload constructor
Test::Test(long tst_nbr_arg) // Overload constructor
{
  tst_nbr=tst_nbr_arg;
} // end Default constructor
long Test::tst_nbr_get()const{return tst_nbr;} // Get tst_nbr

// Initialize static members

////////////////////////////////////////////////////////////////////////

// SzDstFnc class

// Non-static members
SzDstFnc::SzDstFnc(const std::string sng) // Default constructor
{
  fnc_sng=sng;
} // end Default constructor
SzDstFnc::~SzDstFnc(){ // 
  // Purpose: Destructor for SzDstFnc class
  // Since SzDstFnc class has virtual functions, it should have virtual destructor DeD98 p. 580
  ;
} // end SzDstFnc destructor

// Initialize static members
sng2sng_map SzDstFnc::opt_map=SzDstFnc::opt_map_mk();
sng2sng_map SzDstFnc::opt_map_mk(){ // Create abbreviation map
  sng2sng_map opt_map_tmp;
  opt_map_tmp.insert(sng2sng_map::value_type("sin","sin"));
  opt_map_tmp.insert(sng2sng_map::value_type("sine","sin"));
  opt_map_tmp.insert(sng2sng_map::value_type("cos","cos"));
  opt_map_tmp.insert(sng2sng_map::value_type("cosine","cos"));
  // NB: Return a value to initialize a static class member
  return opt_map_tmp;
} // end SzDstFnc::opt_map_mk()
std::string SzDstFnc::opt2abb(const std::string opt_sng){ // Option to abbreviation mapper
  // NB: Return a value to initialize a static class member
  sng2sng_map::const_iterator itr;
  if(dbg_lvl_get() == dbg_crr){
    std::cout << "There are currently " << opt_map.size() << " option abbreviations" << std::endl;
    std::cout << "Option\t" << "Abbreviation" << std::endl;
    for(itr=opt_map.begin();itr!=opt_map.end();++itr)
      std::cout << itr->first << '\t' << itr->second << std::endl;
  } // endif dbg
  itr=opt_map.find(opt_sng);
  if(itr == opt_map.end()) err_prn("SzDstFnc::opt2abb()",opt_sng+" is unknown");
  return itr->second;
} // end SzDstFnc::opt2abb()
//prc_cmp SzDstFnc::operator()(const prc_cmp &r)const{ // Evaluate function
//  std::cerr << "ERROR calling SzDstFnc::operator()"" << std::endl;
//} // end Sine::operator

////////////////////////////////////////////////////////////////////////

// Sine class

Sine::Sine(const std::string sng) : // Default constructor
  SzDstFnc(sng) // Call base class constructor
{
  foo=1.0;
} // end Default constructor
prc_cmp Sine::operator()(const prc_cmp &ngl)const{ // Evaluate function
  return std::sin(ngl);
} // end Sine::operator
prc_cmp Sine::eval(const prc_cmp &ngl)const{return std::sin(ngl);} // Evaluate function

////////////////////////////////////////////////////////////////////////

// Cosine class

Cosine::Cosine(const std::string sng) : // Default constructor
  SzDstFnc(sng) // Call base class constructor
{
  foo=1.0;
} // end Default constructor
prc_cmp Cosine::operator()(const prc_cmp &ngl)const{ // Evaluate function
  return std::cos(ngl);
} // end Cosine::operator
prc_cmp Cosine::eval(const prc_cmp &ngl)const{return std::cos(ngl);} // Evaluate function

////////////////////////////////////////////////////////////////////////

/* From g++FAQ 2.7.2: 
g++ does not automatically instantiate templates defined in other
files.  Because of this, code written for cfront will often produce
undefined symbol errors when compiled with g++.  You need to tell g++
which template instances you want, by explicitly instantiating them in
the file where they are defined.  For instance, given the files */
// 1998/05/27 The above caveat seems to hold for gcc-2.8.x as well
//template class LNF<prc_cmp>;

