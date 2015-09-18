// $Header: 
// Purpose: OpenMP (OMP) shared memory parallelism (SMP) utilities

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <omp.hh> // OpenMP utilities

#ifndef OMP_HH
#define OMP_HH

// Standard C++ headers
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class
#include <thread> // Standard C++-11 thread class

// Standard C headers
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv
#include <cstring> // strcmp...

// 3rd party vendors
#ifdef _OPENMP
# include <omp.h> /* OpenMP pragmas */
#endif /* not _OPENMP */

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Forward declarations

// Typedefs

// Declare functions with C++ linkages

void
cxx11_hello(); // [fnc] Test CXX-11 multi-threading environment

int // O [nbr] Thread number in parallel regions
cxx11_ini // [fnc] Set up CXX-11 multi-threading environment
(const int thr_nbr); // I [nbr] User-requested thread number

int // O [nbr] Thread number in parallel regions
openmp_ini // [fnc] Set up OpenMP multi-threading environment
(const int thr_nbr); // I [nbr] User-requested thread number

#endif // OMP_HH
