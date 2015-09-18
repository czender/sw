// $Id$ 

// Purpose: Implementation (declaration) of new classes for Mie program

/* Copyright (C) 1997--2014 Charlie Zender
   You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <mie_cls.hh> // Class-specific definitions

// Test class

// Non-static members
Test::Test(long tst_nbr_arg) // [fnc] Constructor
{
  tst_nbr=tst_nbr_arg; // [nbr] Number of tst
} // end Default constructor
long Test::tst_nbr_get()const{return tst_nbr;} // [nbr] Number of tst

// Initialize static members


