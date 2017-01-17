// $Id$ 

// Purpose: Calendar algorithms

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <cln.hh> // Calendar algorithms

#ifndef CLN_HH
#define CLN_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers 
#include <cassert> // Assertions
#include <cstdio> // stderr, FILE, NULL, etc.
#include <cstdlib> // abort, exit, getopt, malloc, strtod, strtol

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <cln_cst.hh> // Calendar constants: Live code, not a header

// Typedefs
typedef struct {
  short yyyy; // [yyyy] Year
  short mm; // [mth] Month [1..12]
  short dd; // [day] Day [1..31]
} yyyymmdd_sct;

// Declare functions that have C++ linkages
yyyymmdd_sct yyyymmdd_prs(const long yyyymmdd);
bool is_leap_yr(const short yyyy);
short dpy_get(const short yyyy);
short dpm_get(const short mm);

short // [day] Days in month
dpmyr
(const short mm, // [mth] Month [1..12]
 const short yyyy); // [yyyy] Year

long // O [yyyymmdd] Date at end
yyyymmdd_ncr
(const long yyyymmdd_srt, // I [yyyymmdd] Date at start
 const long day_nbr); // I [day] Days to increment

long // O [yyyymmdd] Date at end
yrdoy2yyyymmdd
(const short yyyy, // I [yyyy] Year 
 const prc_cmp doy); // I [day] Day of year

prc_cmp // O [day] Day of year
yyyymmdd2doy(const long yyyymmdd); // I [yyyymmdd] Date   

long // O [yyyymmdd]
yyyymmdd2lng
(const short yyyy, // I [yyyy] Year in yyyy format
 const short mm, // I [mth] Month [1..12]
 const short dd); // I [day] Day [1..31]

long // O [yyyymmdd]
yyyymmdd2lng(const yyyymmdd_sct yyyymmdd); // I [yyyymmdd] Date

std::string // O [sng] Date in yyyymmdd format
date2vrs // [fnc] Convert date string into yyyymmdd
(const std::string &date_sng); // [sng] POSIX date

long // O [yyyymmdd] Date
date2yyyymmdd // [fnc] Convert date string into yyyymmdd
(const std::string &date_sng); // [sng] POSIX date

#endif // CLN_HH  




