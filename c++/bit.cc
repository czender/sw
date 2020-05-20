// $Id$ 

// Purpose: Bit manipulation utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <bit.hh> // Bit manipulation utilities

// Global functions with C linkages begin

// Declare global functions with C++ linkages

void 
bit_prn_uint // [fnc] Display bits
(const unsigned int val) // [msk] Value to print in binary
{
  // Purpose: Print bits of input variable
  // NB: This routine is deprecated in favor of bit_sng_tpl()

  long idx; // [idx] Counting index
  const int byt_nbr(sizeof(unsigned int)); // [nbr] Number of bytes in input value
  const short bit_per_byt(8); // [nbr] Bits per byte
  const int bit_nbr(byt_nbr*bit_per_byt); // [nbr] Number of bits in input value

  // Initialize mask to 1 in leftmost bit, 0 elsewhere by leftshifting binary 1 by bit_nbr-1 bits
  // fxm: next line generates "bit.cc", line 24: warning #68: integer conversion resulted in a change of sign
  unsigned int msk(1 << (bit_nbr-1)); // [msk] Mask

  // Technique modified from DeD98 p. 808 
  std::cout << "Internal representation of " << val << " is ";
  for(idx=0;idx<bit_nbr;idx++){
    std::cout << (msk & val ? '1' : '0');
    // Rightshift mask to examine next bit of value in next iteration
    msk >>= 1; // [msk] Mask
  } // end loop over idx
  std::cout << std::endl;

} // end bit_prn_uint()

std::string // O [sng] String describing typeid
typeid2sng // [fnc] Convert typeid to descriptive string
(std::string typeid_nm_sng) // I [sng] typeid.name()
{
  // Purpose: Convert RTTI type to string
  sng2sng_map typeid2sng_map; // Key is GCC typeid.name(), value is C++ type name
    
  const bool bln_foo('1'); typeid2sng_map.insert(sng2sng_map::value_type(typeid(bln_foo).name(),"bool")); // [sng]
  const char chr_foo('1'); typeid2sng_map.insert(sng2sng_map::value_type(typeid(chr_foo).name(),"char")); // [sng]
  const double dbl_foo(1.0); typeid2sng_map.insert(sng2sng_map::value_type(typeid(dbl_foo).name(),"double")); // [sng]
  const std::complex<double> dcx_foo(1.0,0.0); typeid2sng_map.insert(sng2sng_map::value_type(typeid(dcx_foo).name(),"std::complex<double>")); // [sng]
  const float flt_foo(1.0); typeid2sng_map.insert(sng2sng_map::value_type(typeid(flt_foo).name(),"float")); // [sng]
  const std::complex<float> fcx_foo(1.0,0.0); typeid2sng_map.insert(sng2sng_map::value_type(typeid(fcx_foo).name(),"std::complex<float>")); // [sng]
  const int int_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(int_foo).name(),"int")); // [sng]
  const long double lng_dbl_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(lng_dbl_foo).name(),"long double")); // [sng]
  const long lng_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(lng_foo).name(),"long")); // [sng]
  const long long lng_lng_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(lng_lng_foo).name(),"long long")); // [sng]
  const short sht_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(sht_foo).name(),"short")); // [sng]
  const signed char schr_foo('1'); typeid2sng_map.insert(sng2sng_map::value_type(typeid(schr_foo).name(),"signed char")); // [sng]
  const std::string sng_foo("1"); typeid2sng_map.insert(sng2sng_map::value_type(typeid(sng_foo).name(),"std::string")); // [sng]
  const unsigned char uchr_foo('1'); typeid2sng_map.insert(sng2sng_map::value_type(typeid(uchr_foo).name(),"unsigned char")); // [sng]
  const unsigned int uint_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(uint_foo).name(),"unsigned int")); // [sng]
  const signed int sint_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(sint_foo).name(),"signed int")); // [sng]
  const unsigned long long ulng_lng_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(ulng_lng_foo).name(),"unsigned long long")); // [sng]
  const unsigned long ulng_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(ulng_foo).name(),"unsigned long")); // [sng]
  const unsigned short usht_foo(1); typeid2sng_map.insert(sng2sng_map::value_type(typeid(usht_foo).name(),"unsigned short")); // [sng]
  
  if(dbg_lvl_get() == dbg_crr){
    // fxm: A more useful method is to make this map a static class member callable anytime
    // so map is only computed once
    sng2sng_map::const_iterator typeid2sng_itr;
    for(typeid2sng_itr=typeid2sng_map.begin();typeid2sng_itr!=typeid2sng_map.end();++typeid2sng_itr)
      std::cout << "C++ type " << typeid2sng_itr->second << " has OS-specific RTTI typeid.name() = " << typeid2sng_itr->first << std::endl;
  } // end if dbg

  return typeid2sng_map[typeid_nm_sng];
} // end typeid2sng()
