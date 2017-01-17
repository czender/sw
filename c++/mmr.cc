// $Id$ 

// Purpose: Implementation of memory management namspace functions

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Eschew storing memory statistics now, to improve performance
   It is more efficent to analyze usage using perl than to handle usage stats during execution */

#include <mmr.hh> // Memory Management

 void *mmr::newVar(size_t size,const char *file,const int line)
 {
   void *ptr=new char[size];
   std::cerr << "new: Allocating |" << size << "| bytes in file |" << file  << "|, line |" << line << "|, address is |" << ptr << std::endl;
   return ptr;
 }

 void *mmr::newArr(size_t size,const char *file,const int line)
 {
   void *ptr = new char[size];
   std::cerr << "new[]: Allocating |" << size << "| bytes in file |" << file << "|, line |" << line << "|, address is |" << ptr << std::endl;
   return ptr;
 }

 void mmr::deleteVar(void *ptr,const char *file,const int line)
 {
   std::cerr << "delete: Freeing memory allocated at file |" << file << "|, line |" << line  << "|, address is |" << ptr << std::endl;
   delete[] (char *) ptr;
 }

 void mmr::deleteArr(void *ptr,const char *file,const int line)
 {
   std::cerr << "delete[]: Freeing memory allocated at file |" << file << "|, line |" << line << "|, address is |" << ptr << std::endl;
   delete[] (char *) ptr;
 }

