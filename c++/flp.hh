// $Id$ 

// Purpose: Floating point utilities, constants

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Header includes mathematical constants in mth_cst.hh
   If routine requires utility functions AND constants, then #include flp.hh
   If routine requires only constants, then #include mth_cst.hh
   Usage:
   #include <flp.hh> // Floating point utilities, constants */

#ifndef FLP_HH // Contents have not yet been inserted in current source file  
#define FLP_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers
#include <cstdlib> // abort, exit, getopt, malloc, strtod, strtol

// 3rd party vendors
#include <gsl/gsl_errno.h> // GNU Scientific Library error handling

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Forward declarations

// Typedefs

// Declare functions with C++ linkages

void gsl_error_handler_csz // Custom GSL error handler
(const char * msg, // I [sng] Reason error occurred
 const char * fl_err, // I [sng] File where error occurred
 int ln_err, // I [idx] Line where error occurred
 int gsl_errno); // I [enm] GSL error code in gsl_errno.h

#endif // FLP_HH  






