// $Id$ 

// Purpose: Provides a File Input/Output Layer

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <fio.hh> // File Input/Output Layer

#ifndef FIO_HH
#define FIO_HH

// Standard C++ headers
#include <iostream> // std::cout, cin, std::cerr
#include <string> // Standard C++ string class

// Standard C headers
#include <sys/stat.h> // stat()

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// fxm integrate input and output paths into this design
namespace fio{ // [nms] File Input Output Layer Namespace

  static std::string data_path=""; // [sng] Data Directory

  std::string data_path_get(void); // [sng] Get data file path

  std::string // [fnc] File name prepended by data directory
  data_file_path_get // [fnc] Expand file path to include data directory
  (const std::string data_file); // [sng] File name
  // end fio::data_file_path_get() prototype

  int // [enm] Return success code
  data_path_set // [fnc] Set data directory
  (const std::string data_path_arg); // [sng] Data Directory
  // end fio::data_path_set() prototype

  // Diagnose file and filesystem properties
  int // [enm] Return success code
  fl_stat_dgn // [fnc] Diagnose file and filesystem properties
  (const std::string fl_nm); // [sng] Filename

  int // [enm] Return success code
  tst // [fnc] Test IO layer
  (const std::string tst_sng); // [sng] Test string
  // end fio::tst() prototype

} // end namespace fio

#endif // FIO_HH

