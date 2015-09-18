// $Id$ 

// Purpose: Description (definition) bit manipulation utilities

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <bit.hh> // Bit manipulation utilities

#ifndef BIT_HH // Contents have not yet been inserted in current source file
#define BIT_HH

// C++ headers
#include <cstring> // strcmp...
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL multimap and map classes
#include <typeinfo> // Standard C++ header for typeid, type_info

// Standard C headers

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Namespaces

// 3rd party vendors

// Typedefs
typedef std::map<std::string,std::string,std::less<std::string> > sng2sng_map; // String-to-string map

// Prototype global functions with C++ linkages

void 
bit_prn_uint // [fnc] Display bits
(const unsigned int val); // [msk] Value to print in binary
// end bit_prn_uint() prototype

std::string // O [sng] String describing typeid
typeid2sng // [fnc] Convert typeid to descriptive string
(std::string typeid_nm_sng); // I [sng] typeid.name()
// end typeid2sng() prototype

/* Template definitions must be included in every file (translation unit)
   unless the compiler supports the "export" keyword (gcc 2.96 does not). 
   Thus actual function definitions must appear in .hh file being #include'd */

template<class prc_T>std::string // O [sng] String containing binary representation of value
bit_sng_tpl // [fnc] Construct string with binary representation of value
(const prc_T val); // [msk] Value to print in binary
// end bit_sng_tpl() prototype

template<class prc_T>std::string // O [sng] String containing binary representation of value
bit_sng_tpl // [fnc] Construct string with binary representation of value
(const prc_T val) // [msk] Value to print in binary
{
  /* Purpose: Construct string with binary representation of value
     Routine is a good example of usefulness of C++ templates
     Routine works for [signed/unsigned] char, short, int, long, long long
     20150123: Trying to get routine to work for float, double
     http://math.stackexchange.com/questions/144659/an-algorithm-to-convert-float-number-to-binary-representation
     Wikipedia IEEE-754:
     http://en.wikipedia.org/wiki/Single-precision_floating-point_format
     Online Calculator:
     http://www.h-schmidt.net/FloatConverter/IEEE754.html

     Techniques for roundig by masking and shifting
     http://blog.frama-c.com/index.php?post/2013/05/03/nearbyintf2
     assert (sizeof(unsigned int) == sizeof(float));
     unsigned int u;
     memcpy(&u,&f,sizeof(float));
     Exponent:
     int exp = ((u>>23) & 255) - 127;
     Explicit significand: u & 0x7fffff */
  long idx; // [idx] Counting index
  const int byt_nbr(sizeof(prc_T)); // [nbr] Number of bytes in input value
  const short bit_per_byt(8); // [nbr] Bits per byte
  const int bit_nbr(byt_nbr*bit_per_byt); // [nbr] Number of bits in input value
  std::string sng(""); // [sng] String containing binary representation of value
  const std::string sbr_nm("bit_sng_tpl"); // [sng] Subroutine name

  unsigned char msk_chr(1);
  unsigned short msk_sht(1);
  unsigned int msk_int(1);
  unsigned long msk_lng(1);
#ifdef HAVE_LONG_LONG
  unsigned long long msk_lng_lng(1);
#endif // !HAVE_LONG_LONG

  if(dbg_lvl_get() == dbg_io){
    std::cout << "byt_nbr = " << byt_nbr << ", bit_nbr = " << bit_nbr << std::endl;
    std::cout << sbr_nm+"() reports sizeof() storage class prc_T is " << sizeof(prc_T) << " bytes = " << sizeof(prc_T)*bit_per_byt << " bits" << std::endl;
  } // endif dbg

  // Sanity check
  /* How should this algorithm should work for negative integers?
     Represenation of negative integers and two's complements in DeD98 p. 1075 is good
     Definition of binary representation of negative integer is such that negative integer
     plus positive integer must equal zero in all bits
     NyL97 p. A70 give very nice discussion of internal representation of floats and integers */

  /* Copy mask into integer variable with same number of bytes as val
     Initialize mask to 1 in leftmost bit, 0 elsewhere by leftshifting binary 1 by bit_nbr-1 bits
     Technique modified from DeD98 p. 808 */
  switch(byt_nbr){
  case 1:
    unsigned char val_chr;
    std::memcpy(&val_chr,&val,sizeof(prc_T));
    /* NB: for some reason type of arguments to bit-shifting must not exceed shifted type
       "msk_sht<<=bit_nbr_msk-1" produces warning but msk_sht<<=15" does not
       Hard-coded naked constants do not emit warnings in clang or g++
       Arithmetic with constants emits warnings g++, not clang
       Casted arithmetic also emits warnings g++, not clang
       Hence we simply hardcode the shift by type to keep g++ happy */
    //msk_chr<<=bit_nbr_msk-1;
    //msk_chr<<=(unsigned char)(bit_nbr_msk-1);
    //msk_chr<<=8*sizeof(prc_T)-1;
    //msk_chr<<=(unsigned char)(8*sizeof(prc_T)-1);
    msk_chr<<=7;
    for(idx=0;idx<bit_nbr;idx++){
      sng+=(msk_chr & val_chr ? "1" : "0");
      msk_chr >>= 1;
    } // end loop over idx
    break;
  case 2:
    unsigned short val_sht;
    std::memcpy(&val_sht,&val,sizeof(prc_T));
    msk_sht<<=15;
    //    msk_sht<<=bit_nbr_msk-1;
    for(idx=0;idx<bit_nbr;idx++){
      sng+=(msk_sht & val_sht ? "1" : "0");
      msk_sht >>= 1;
    } // end loop over idx
    break;
  case 4:
    unsigned int val_int;
    std::memcpy(&val_int,&val,sizeof(prc_T));
    msk_int<<=31;
    for(idx=0;idx<bit_nbr;idx++){
      sng+=(msk_int & val_int ? "1" : "0");
      msk_int >>= 1;
    } // end loop over idx
    break;
  case 8:
    unsigned long val_lng;
    std::memcpy(&val_lng,&val,sizeof(prc_T));
    msk_lng<<=8*sizeof(prc_T)-1;
    for(idx=0;idx<bit_nbr;idx++){
      sng+=(msk_lng & val_lng ? "1" : "0");
      msk_lng >>= 1; // [msk] Mask
    } // end loop over idx
    break;
#ifdef HAVE_LONG_LONG
  case 16:
    unsigned long long val_lng_lng;
    std::memcpy(&val_lng_lng,&val,sizeof(prc_T));
    msk_lng_lng<<=8*sizeof(prc_T)-1;
    for(idx=0;idx<bit_nbr;idx++){
      sng+=(msk_lng_lng & val_lng_lng ? "1" : "0");
      msk_lng_lng >>= 1; // [msk] Mask
    } // end loop over idx
    break;
#endif // !HAVE_LONG_LONG
  default:
    std::cout << "ERROR: missed case statement in bit.hh" << std::endl;
    return sng;
  } // end switch 

  return sng; // [sng] String containing binary representation of value
} // end bit_sng_tpl()

// Define inline'd functions in header so source is visible to calling files

#endif // BIT_HH  






