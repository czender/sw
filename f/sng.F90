program sng

  ! Purpose: Demonstrate command-line argument-handling and string manipulation
  ! This program is intended to be fully F90-standard compliant except it relies
  ! on the presence of external iargc() and getarg() routines normally supplied
  ! by Fortran compilers. 

  ! Copyright (C) 1997--present Charlie Zender
  ! License: This file is in the Public Domain

  ! The original author of this software, Charlie Zender, wants to receive
  ! your suggestions, thanks, bug-reports, and patches to improve it.
  ! Charlie Zender <zender at uci dot edu>
  ! Department of Earth System Science
  ! University of California, Irvine
  ! Irvine, CA 92697-3100
  
  ! Compilation: 
  ! g95:
  ! g95 -o dbg_mdl.o -c dbg_mdl.F90
  ! g95 -o sng_mdl.o -c sng_mdl.F90
  ! g95 -o sng sng_mdl.o dbg_mdl.o sng.F90
  ! gfortran:
  ! gfortran -o dbg_mdl.o -c dbg_mdl.F90
  ! gfortran -o sng_mdl.o -c sng_mdl.F90
  ! gfortran -o sng sng_mdl.o dbg_mdl.o sng.F90
  ! Intel ifc compiler:
  ! ifc -o dbg_mdl.o -c dbg_mdl.F90
  ! ifc -o sng_mdl.o -c sng_mdl.F90
  ! ifc -o sng sng_mdl.o dbg_mdl.o sng.F90 -lPEPCF90
  ! Portland Group pgf90 compiler:
  ! pgf90 -o dbg_mdl.o -c dbg_mdl.F90
  ! pgf90 -o sng_mdl.o -c sng_mdl.F90
  ! pgf90 -o sng sng_mdl.o dbg_mdl.o sng.F90

  ! Usage:
  ! sng --dbg=0 --dbl=1.234e234 --flt=-5.67e-37 --int=12345678 --lgc=T --sng="GNU's Not UNIX" --drc_out=${HOME}/output/directory
  ! sng -D 5 --dbl 1.234e234

  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use sng_mdl ! [mdl] String manipulation

  implicit none
  ! Parameters
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$HeadURL$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)

  ! Command-line parsing
  character(80)::arg_val ! [sng] command-line argument value
  character(200)::cmd_ln ! [sng] command-line
  character(80)::opt_sng ! [sng] Option string
  character(2)::dsh_key ! [sng] command-line dash and switch
  character(200)::prg_ID ! [sng] Program ID

  integer::arg_idx ! [idx] Counting index
  integer::arg_nbr ! [nbr] Number of command-line arguments
  integer::opt_lng ! [nbr] Length of option

  ! Set defaults for command-line options 
  character(80)::drc_in="/home/user/input/dir"//nlc ! [sng] Input directory
  character(80)::drc_out="" ! [sng] Output directory
  character(80)::fl_in="in.nc"//nlc ! [sng] Input file
  character(80)::fl_out="foo.nc"//nlc ! [sng] Output file
  character(80)::sng_foo="Default value"//nlc ! [sng] String
  integer::int_foo=1 ! [nbr] Integer
  logical::lgc_foo=.false. ! [flg] Logical
  real(selected_real_kind(p=6))::flt_foo=0.0 ! [frc] Float
  real(selected_real_kind(p=12))::dbl_foo=0.0 ! [frc] Double

  ! Main code
  dbg_lvl=0 ! [idx] Causes DDD source window to display this file
  ! CEWI for lf95: Initialize cmd_ln before taking length 
  call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

  if (.true.) then
     write (6,*) "Before command-line parsing and string manipulation: "
     write (6,*) "main() reports cmd_ln = ",cmd_ln(1:ftn_strlen(cmd_ln))
     write (6,*) "main() reports dbl_foo = ",dbl_foo
     write (6,*) "main() reports flt_foo = ",flt_foo
     write (6,*) "main() reports int_foo = ",int_foo
     write (6,*) "main() reports lgc_foo = ",lgc_foo
     write (6,*) "main() reports sng_foo = ",sng_foo(1:ftn_strlen(sng_foo))
     write (6,*) "main() reports fl_in = ",fl_in(1:ftn_strlen(fl_in))
     write (6,*) "main() reports fl_out = ",fl_out(1:ftn_strlen(fl_out))
     write (6,*) "main() reports drc_in = ",drc_in(1:ftn_strlen(drc_in))
     write (6,*) "main() reports drc_out = ",drc_out(1:ftn_strlen(drc_out))
  endif ! endif .true.

  call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

  arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments
  arg_idx=1 ! [idx] Counting index
  loop_while_options: do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
     dsh_key=arg_val(1:2) ! [sng] First two characters of option
     if (dsh_key == "--") then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) error stop "Long option has no name"
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (dbg_lvl >= dbg_io) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
             ": DEBUG Double hyphen indicates multi-character option: ", &
             "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
        ! fxm: Change if else if construct to select case but how to handle fall-through cases elegantly?
        if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == "dbl" .or. opt_sng == "dbl_foo" ) then
           call ftn_arg_get(arg_idx,arg_val,dbl_foo) ! [frc] Double
        else if (opt_sng == "drc_in") then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == "drc_out") then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == "fl_in") then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
        else if (opt_sng == "fl_out") then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
        else if (opt_sng == "flt" .or. opt_sng == "flt_foo" ) then
           call ftn_arg_get(arg_idx,arg_val,flt_foo) ! [frc] Float
        else if (opt_sng == "int" .or. opt_sng == "int_foo" ) then
           call ftn_arg_get(arg_idx,arg_val,int_foo) ! [nbr] Integer
        else if (opt_sng == "lgc" .or. opt_sng == "lgc_foo" ) then
           call ftn_arg_get(arg_idx,arg_val,lgc_foo) ! [lgc] Logical
        else if (opt_sng == "sng" .or. opt_sng == "sng_foo") then
           call ftn_arg_get(arg_idx,arg_val,sng_foo) ! [sng] String
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     else if (dsh_key(1:1) == '-') then ! endif long option
        ! Handle short options
        if (dsh_key == "-D") then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (dsh_key == "-f") then
           call ftn_arg_get(arg_idx,arg_val,flt_foo) ! [frc] Float
        else if (dsh_key == "-i") then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
        else if (dsh_key == "-l") then
           call ftn_arg_get(arg_idx,arg_val,int_foo)
        else if (dsh_key == "-p") then
           lgc_foo=.not.lgc_foo
        else if (dsh_key == "-o") then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
        else if (dsh_key == "-s") then
           call ftn_arg_get(arg_idx,arg_val,sng_foo) ! [sng] String
        else if (dsh_key == "-v") then
           goto 1000 ! Goto exit with error status
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif arg_val
     else ! endif short option
        ! Last argument(s) may be positional
        if (arg_idx == arg_nbr) then
           call ftn_strcpylsc(fl_in,arg_val) ! [sng] Input file
        else if (arg_idx == arg_nbr+1) then
           call ftn_strcpylsc(fl_out,arg_val) ! [sng] Output file
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! arg_idx != arg_nbr
     end if ! end position arguments
  end do loop_while_options ! end while (arg_idx <= arg_nbr)
  
  ! Compute any quantities that might depend on command-line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file

  if (.true.) then
     write (6,*) "After command-line parsing and string manipulation: "
     write (6,*) "main() reports cmd_ln = ",cmd_ln(1:ftn_strlen(cmd_ln))
     write (6,*) "main() reports dbl_foo = ",dbl_foo
     write (6,*) "main() reports flt_foo = ",flt_foo
     write (6,*) "main() reports int_foo = ",int_foo
     write (6,*) "main() reports lgc_foo = ",lgc_foo
     write (6,*) "main() reports sng_foo = ",sng_foo(1:ftn_strlen(sng_foo))
     write (6,*) "main() reports fl_in = ",fl_in(1:ftn_strlen(fl_in))
     write (6,*) "main() reports fl_out = ",fl_out(1:ftn_strlen(fl_out))
     write (6,*) "main() reports drc_in = ",drc_in(1:ftn_strlen(drc_in))
     write (6,*) "main() reports drc_out = ",drc_out(1:ftn_strlen(drc_out))
  endif ! endif .true.

1000 continue ! Jumping point for quick exit

end program sng
