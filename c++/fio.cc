// $Id$ 

// Implementation of File Input/Output Layer

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <fio.hh> // File IO Layer

std::string fio::data_path_get(void){return data_path;} // [sng] Data path 

std::string // [fnc] File name prepended by data directory
fio::data_file_path_get // [fnc] Expand file path to include data directory
(const std::string data_file) // [sng] File name
{
  // Purpose: Expand file path to include data directory
  std::string fl_nm_fll; // [sng] Full file name
  // fl_nm_fll=drc_pfx(data_path,data_file);
  fl_nm_fll=drc_pfx("/data/zender/aca",data_file);
  // if(!data_path.empty()) fl_nm_fll=drc_pfx(data_path,data_file);
  // if(data_path.empty()) fl_nm_fll=data_file;
  /* 20120409: ccc segfaults here for unknown reason
     20130119: seems to be because data_path is undefined (NULL)
     This even though ccc calls data_path_set() early
     Possibly, though, after mnr structure is initialized
     So possibly not early enough 
     How to set data_path before structures initialized? */
  //  return drc_pfx(data_path,data_file);
  return fl_nm_fll;
} // end fio::data_path_get()

int // [enm] Return success code
fio::data_path_set // [fnc] Set data directory
(const std::string data_path_arg) // [sng] Data Directory
{
  // Purpose: Set data directorypp
  const std::string sbr_nm("fio::data_path_set"); // [sng] Subroutine name
  if(data_path_arg.empty()) return false; else data_path=data_path_arg;
  
  // Check for directory existance
  int rcd(0); // [enm] Return success code
  struct stat stat_sct;
  rcd+=stat(data_path.c_str(),&stat_sct);
  if(rcd == -1) err_prn(prg_nm_get(),sbr_nm,"Unable to open "+data_path);
  return rcd; // [enm] Return success code
} // end fio::data_path_set()

// Diagnose file and filesystem properties
int // [enm] Return success code
fio::fl_stat_dgn // [fnc] Diagnose file and filesystem properties
(const std::string fl_nm) // [sng] Filename
{
  // Purpose: Diagnose file and filesystem properties
  std::string sbr_nm("fio:fl_stat_dgn"); /* [sng] Subroutine name */
  
  //blksize_t fl_sys_blk_sz; /* [nbr] File system blocksize for I/O */

  int rcd_sys;
  
  mode_t fl_md;
  mode_t fl_usr_md;
  mode_t fl_usr_wrt_md;
  
  struct stat stat_sct;

  /* Output file now guaranteed to exist. Perform stat() to check its permissions. */
  rcd_sys=stat(fl_nm.c_str(),&stat_sct);

  /* 20120228 Ensure output file is writable even when input file is not 
     stat structure includes st_mode field which includes following flags:
     mode_t st_mode
     S_IRWXU    00700     mask for file owner permissions
     S_IWUSR    00200     owner has write permission
     Method of checking: 
     First  bit-wise "and" (& S_IRWXU) uses mask to strips full, multibyte, file mode flag of all but user/owner byte 
     Second bit-wise "and" (& S_IWUSR) is only "true" (non-zero) is owner write permission is set */
  fl_md=stat_sct.st_mode;
  /* Blocksize information in stat structure:
     blksize_t st_blksize blocksize for file system I/O */
  // fl_sys_blk_sz=stat_sct.st_blksize;
  fl_usr_md=fl_md & S_IRWXU;
  fl_usr_wrt_md=fl_usr_md & S_IWUSR;
  if(dbg_lvl_get() >= dbg_scl) std::cerr << prg_nm_get() << ": " << sbr_nm << " reports permissions for file " << fl_nm << " are (octal) " << (unsigned long)fl_md << std::endl;
  // if(dbg_lvl_get() >= nco_dbg_scl) (void)fprintf(stderr,"%s: %s reports preferred filesystem I/O block size: %ld bytes\n",prg_nm_get(),fnc_nm,(long)fl_sys_blk_sz);

  return rcd_sys; // [enm] Return success code
} // end fio::fl_stat_dgn()

// fxm write test code for fio.cc
int // [enm] Return success code
fio::tst // [fnc] Test IO layer
(const std::string tst_sng) // [sng] Test string
{
  // Purpose: Test fio namespace
  int rcd(0); // [enm] Return success code
  return rcd; // [enm] Return success code
} // end fio::tst()

// fio class
