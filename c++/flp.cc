// $Id$ 

// Purpose: Floating point utilities, constants for C++ programs

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <flp.hh> // Floating point utilities, constants

void gsl_error_handler_csz // Custom GSL error handler
(const char * msg, // I [sng] Reason error occurred
 const char * fl_err, // I [sng] File where error occurred
 int ln_err, // I [idx] Line where error occurred
 int gsl_errno) // I [enm] GSL error code in gsl_errno.h
{
  /* Purpose: Custom GSL error handler
     Usage: Register in main() with gsl_set_error_handler(gsl_crror_handler_csz)
     See gsl.pdf Chapter 3 */
  const std::string sbr_nm("gsl_error_handler_csz"); // [sng] Subroutine name
  std::cerr << prg_nm_get() << ": Homebrew GSL error handler "+sbr_nm+"() called from File = " << fl_err << ", line = " << ln_err << std::endl;
  std::cerr << prg_nm_get() << ": gsl_errno = " << gsl_errno << ", message = " << msg << std::endl;
  switch(gsl_errno){
  case GSL_SUCCESS: // [enm] No problemo
    break;
  case GSL_EUNDRFLW: // [enm] Underflow
    std::cerr << prg_nm_get() << ": " << sbr_nm << " handler allows underflow" << std::endl;
    break;
  default:
    std::cerr << prg_nm_get() << ": Calling std::abort() (i.e., dumping core)" << std::endl;
    std::abort(); // [fnc] Dump core-file
    // std::exit(); // [fnc] Exit without core dump
  } // end switch
} // end gsl_error_handler_csz()
