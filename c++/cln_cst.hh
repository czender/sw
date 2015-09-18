// $Id$ 

// Purpose: Calendar constants

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage: This header should only be included by cln.hh
// #include <cln_cst.hh> // Calendar constants: Live code, not a header

#ifndef CLN_CST_HH // Contents have not yet been inserted in current source file  
#define CLN_CST_HH

// Variables with file scope

namespace cln{ // [nms] Calendar namespace
  
  // Following variables describe 365 day model year
  const long mpy(12); // (12) [nbr] Model months per year
  const long dpy_365(365); // (365) [nbr] Model days per year
  /* const long dpm[mpy+1]= // [day] Days per month
     {0, // Zeroth element enables array addressing by month number
     31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // end dpm[] */
  const long fdom_doy_365[mpy+1]= // [day] Day of year of first day of month 
  {0, // Zeroth element enables array addressing by month number
   01, 32, 60, 91,121,152,182,213,244,274,305,335}; // end fdom_doy[]
  const long ldom_doy_365[mpy+1]= // [day] Day of year of last day of month
  {0, // Zeroth element enables array addressing by month number
   31, 59, 90,120,151,181,212,243,273,304,334,365}; // end ldom_doy[]
  const std::string day_nm[7]={"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};
  const std::string mth_nm[mpy]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
  
} // end calendar namespace cln

#endif // CLN_CST_HH  




