// $Header: 
// Purpose: Description (definition) of memory management utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mmr.hh> // Memory management

#ifndef MMR_HH // Contents have not yet been inserted in current source file
#define MMR_HH

// Standard C++ headers
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates
#include <time.h> // Standard Time library

// Standard C headers

// 3rd party vendors

// Personal headers

// To handle file and line number, we include the necessary preproc macro (Stroustrup Section 24.3.7)

#define TRACK (__FILE__, __LINE__) // Track individual files and lines

// Prototype global functions with C++ linkages

/* User-Defined new and delete Cannot Be Declared in a Namespace in gcc > 2.94
   Thus, we now use function calls instead of operators */
namespace mmr{ // [nmr] Memory Management Namespace
  // Include file and line info, so calls to new() are delinated as:
  // char *p = new(char[100],TRACK);
  void *newVar(size_t size,const char *file,const int line);
  void *newArr(size_t size,const char *file,const int line);
  void deleteVar(void *ptr,const char *file,const int line);
  void deleteArr(void *ptr,const char *file,const int line);
} // end of nms mmr

// Define inline'd functions in header so source is visible to calling files

#endif // MMR_HH
