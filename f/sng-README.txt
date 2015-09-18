$Id$

Last Updated: 18 March 2009

sng: A Portable Fortran9X/2K getopt_long()-like Command-Line Parser 
     and String Manipulation Library.

Welcome to the only known documentation for sng!
In about 1997 I wrote a package of Fortran utilities to emulate a
subset of the GNU-POSIX getopt() command line argument parsing.
It has run trouble-free since about 2001, so I'm giving it away.
Perhaps others will make it fully GNU getopt() compatible.
The package is hereby known as "sng", which sounds like "string" 
if you lisp and mumble it under a strong running shower.
sng was developed mainly as an alternative to the namelist mechanism, 
to allow altering program inputs efficiently without altering 
a namelist file or recompiling. The sng tarball is available from

http://dust.ess.uci.edu/f

The sng code is well-commented but there is no other documentation. 
This package works with all known Fortran9X compilers. 
(Although PGI pgf90 seems to have undocumented behavior).
It is standard-compliant free-format Fortran-9X upwardly
compatible with F2K (see below). sng contains a nearly complete
collection of Fortran analogues to the C string library:

ftn_strlen() 
ftn_strstr() 
ftn_strcat()
ftn_strcpy()
ftn_strcmp()

And some functions that ease treating Fortran strings like C strings: 
ftn_strini()
ftn_strnul()
ftn_strlsc()
ftn_strfic()
ftn_strprn()
ftn_strpfx()
ftn_drcpfx()
ftn_strcpylsc()

SNG uses these routines internally to construct getopt()-like
routines for parsing command-lines with a simple, efficient API.
SNG parses arguments of any type (float, double, logical, character)
with the same (overloaded) routine, ftn_arg_get().
Yes, the library includes a routine (cmd_ln_get()) which returns the
entire command line for further manipulation and storage.

The library emulates the most important part of GNU-POSIX extended
command line parsing: allowing multi-character option descriptions.
The API supports single-dash and double-dash command line arguments,
e.g., sng --dbg=2 or sng -D 2. Probably the most significant GNU
features lacking in sng are the absence of support for resolution of 
unambiguous abbreviations for options, and the lack of any support
for argument re-ordering after parsing. The front-end tester, sng.F90,
demonstrates a work-around for the former, but there is none for the
latter. I am very interested in accepting patches which improve the
behavior of sng and bring it towards fully GNU-POSIX getopt()
compatibility. 

sng uses the proposed F2K standard routines command_argument_count()
and get_command_argument() to access the command line.
Since few if any F2K compilers exist, sng defines these routines
as wrappers for the de-facto standard (on UNIX, anyway) routines
iargc() and getarg(). Thus sng is a stand-alone self-contained
package. Code written to use sng will not change when F2K compilers
are used. The compiler-native command_argument_count() and
get_command_argument() routines will simply supercede the wrappers
used internally by sng. This will be transparent to the users.

For those compilers which do not define iargc() and getarg(),
sng relies on the f2kcli package

http://www.winteracter.com/f2kcli/downld.htm

f2kcli provides command_argument_count() and get_command_argument()
routines for basically all known Fortran compilers.
Simply compile sng with the -DF2KCLI token and it will "use f2kcli"
for these routines rather than defining wrappers internally.

The parts of the sng package are:
sng-README.txt	Documentation (this file)
sng.F90		Driver program which demonstrates usage of sng
sng_mdl.F90	Module which contains all string manipulation utilities
dbg_mdl.F90	Debugging module used by some sng routines

Note that sng_mdl.F90 uses a free, non-restrictive, non-copyleft
X11-style license, whereas sng.F90 and dbg_mdl.F90 are public domain.
Instructions for compiling and testing sng are contained at the top
of the sng.F90 driver program. They are short and simple so I repeat
them here, with the example output:

zender@lanina:~/f$ pgf90 -o dbg_mdl.o -c dbg_mdl.F90
zender@lanina:~/f$ pgf90 -o sng_mdl.o -c sng_mdl.F90
zender@lanina:~/f$ pgf90 -o sng sng_mdl.o dbg_mdl.o sng.F90
zender@lanina:~/f$ sng --dbg=0 --dbl=1.234e234 --flt=-5.67e-37 --int=12345678 --lgc=T --sng="GNU's Not UNIX" --drc_out=${HOME}/output/directory
 Before command line parsing and string manipulation: 
 main() reports cmd_ln = 
 main() reports dbl_foo =    0.0000000000000000E+000
 main() reports flt_foo =    0.0000000E+00
 main() reports int_foo =             1
 main() reports lgc_foo =   F
 main() reports sng_foo = Default value
 main() reports fl_in = in.nc
 main() reports fl_out = foo.nc
 main() reports drc_in = /home/user/input/dir
 main() reports drc_out = 
 After command line parsing and string manipulation: 
 main() reports cmd_ln = 
 sng --dbg=0 --dbl=1.234e234 --flt=-5.67e-37 --int=12345678 --lgc=T --sng=GNU's Not UNIX --drc_out=/home/zender/output/directory
 main() reports dbl_foo =    1.2339999999999977E+234
 main() reports flt_foo =   -5.6700002E-37
 main() reports int_foo =      12345678
 main() reports lgc_foo =   T
 main() reports sng_foo = GNU's Not UNIX
 main() reports fl_in = /home/user/input/dir/in.nc
 main() reports fl_out = /home/zender/output/directory/foo.nc
 main() reports drc_in = /home/user/input/dir
 main() reports drc_out = /home/zender/output/directory

Enjoy!
Charlie
