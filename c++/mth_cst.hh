// $Id$

// Purpose: Mathematical constants

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* This header is automatically included in mth.hh but must be
   manually included in headers that do not use mth.hh */
// Usage:
// #include <mth_cst.hh> // Mathematical constants, cst_M_PIl...

#ifndef MTH_CST_HH // Contents have not yet been inserted in current source file  
#define MTH_CST_HH

// Variables with file scope

namespace mth{ // [nms] Mathematical constant namespace
  /* Math constants taken from GNU /usr/include/math.h
     Constants are suitable for use with IEEE 128 bit quadruple precision */
  const double cst_M_EULERl(0.577215664901532860606512090082L); // (0.577215664901532860606512090082L) [frc] Euler's constant
  const double cst_M_PIl(3.1415926535897932384626433832795029L); // (3.1415926535897932384626433832795029L) [frc] 3
  const double cst_M_LN2l(0.6931471805599453094172321214581766L); // (0.6931471805599453094172321214581766L) [frc] log_e(2.0)
  const double cst_M_El(2.7182818284590452353602874713526625L); // (2.7182818284590452353602874713526625L) [frc] exp(1.0)
  const double cst_M_SQRT2l(1.4142135623730950488016887242096981L); // (1.4142135623730950488016887242096981L) [frc] sqrt(2.0)
} // end Mathematical constant namespace mth

#endif // MTH_CST_HH  




