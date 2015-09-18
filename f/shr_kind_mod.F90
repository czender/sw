! $Id$ -*-f90-*-

! Purpose: Portable precision definitions for floating point and integer variables

! Copyright (C) 1994--2014 Charlie Zender
! License: GNU General Public License (GPL) Version 3
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...

! Compilation
! cd ~/f;f90 -c -I$HOME/include -o $MY_OBJ_DIR/shr_kind_mod.o shr_kind_mod.F90
! cd ~/f;pgf90 -c -I$HOME/include -o $MY_OBJ_DIR/shr_kind_mod.o shr_kind_mod.F90
! cd ~/f;pgf90 -c -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign -Ktrap=fp -fast -DLINUX -I. -I$HOME/include -I/usr/local/include -o $MY_OBJ_DIR/shr_kind_mod.o shr_kind_mod.F90

! Standard precisions may be selected with predefined kind numbers
! Real kinds:
! 4 = Single precision values with approximately  7 significant digits (4 bytes)
! 8 = Double precision values with approximately 14 significant digits (8 bytes)
! IEEE and Cray storage both return 8 byte reals for precision=12 (but not p=15)
! Integer kinds:
! 1 =  8-bit integers [-2^7,2^7-1]
! 2 = 16-bit integers [-2^15,2^15-1]
! 3 = 32-bit integers [-2^31,2^31-1]
! 4 = 64-bit integers [-2^63,2^63-1]

module shr_kind_mod

  ! Definitions are consistent with CCM and MATCH
  ! Real numbers
  integer,parameter::shr_kind_r16=selected_real_kind(p=24) ! r16: 16B (C long double)
  integer,parameter::shr_kind_r4=selected_real_kind(p=6) ! r4: 4B (C float) default, 8B (C double) possible
#if defined(PRC_DBL) || !defined(PRC_FLT)
  integer,parameter::shr_kind_r8=selected_real_kind(p=12) ! r8: 8B (C double) default, 4B (C float) possible
#endif /* !PRC_DBL */
#ifdef PRC_FLT
  integer,parameter::shr_kind_r8=selected_real_kind(p=6) ! r8: 8B (C double) default, 4B (C float) possible
#endif /* !PRC_FLT */
  integer,parameter::r_ntv=kind(1.0) ! r_ntv: Native real kind
  ! Michael Metcalf's recommendation in comp.lang.fortran 20040824
  integer, parameter::d_ntv=selected_real_kind(2*precision(1.0_r_ntv)) ! d_ntv: double native real kind precision, i.e., double precision
  integer,parameter::DBLKIND=selected_real_kind(p=12) ! DBLKIND: Fixed 8B (double precision) always

  ! Integers
  integer,parameter::shr_kind_i8=selected_int_kind(13) ! i8: 8B (C long long) default, 4B (C int) possible
  integer,parameter::shr_kind_i4=selected_int_kind(6) ! i4: 4B (C int) default, 8B (C long long) possible
  integer,parameter::i_ntv=selected_int_kind(1) ! i_ntv: Native integer kind

  ! Characters
  integer,parameter::chr_lng=256 ! chr_lng: long char
  integer,parameter::chr_sht=80  ! chr_sht: short char

end module shr_kind_mod

