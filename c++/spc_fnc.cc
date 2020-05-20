// $Id$ 

// Implementation (declaration) of special functions

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <spc_fnc.hh> // Special functions

void nr_hook()
{
  // Purpose: hook for Numerical Recipes Routines
  std::string sbr_nm("nr_hook");
  err_prn(sbr_nm,"Routine should not be called");
} // end nr_hook()
