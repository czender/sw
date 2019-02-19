! $Id$ -*-f90-*-

program fff

  ! Purpose: Provide a standard template for developing Fortran programs

  ! Copyright (C) 1994--present Charlie Zender
  ! This software is distributed under the terms of the GNU General Public License
  ! See http://www.gnu.org/copyleft/gpl.html for full license text
  ! The original author of this software, Charlie Zender, seeks to improve
  ! it with your suggestions, contributions, bug-reports, and patches.
  ! Charlie Zender <zender at uci dot edu>
  ! Department of Earth System Science
  ! University of California, Irvine
  ! Irvine, CA 92697-3100
  
  ! History:
  ! 1994--1998: "f" is originally a Fortran77 template program
  ! 1999: "f" incorporates some Fortran90 constructs, but remains fixed-form
  ! 2000: "f" evolves into "fff.F90", a free format Fortran90 template
  ! 2001: "f" is deprecated, Fortran77 code turned into modules
  ! 2002: "fff" is completely modular

  ! NB: Preprocessed Fortran code that contains C comment characters causes warnings
#if ( !defined AIX ) && ( !defined CRAY ) && ( !defined LINUX )
  ! Tag all Fortran code for modification
  ! etags ~/f/*.F ~/f/*.F90 ~/aer/*.F90 ~/aer/*.h ~/map/*.F90 ~/aca/*.F ~/aca/*.F90
  ! */
#endif /* not CRAY */

  ! Compiler-defined (built-in) preprocessor tokens (macros):
  ! cd ${HOME}/f;touch foo.f90;gfortran -cpp -dM foo.f90
  ! cd ${HOME}/f;touch foo.f90;g95 -cpp -dM foo.f90

  ! Compilation:
  ! export NETCDF_ROOT='/usr';export NETCDF4_ROOT='/usr'
  ! cd ${HOME}/f; make -W fff.F90 OPTS=D NETCDF4=Y fff; cd -
  ! cd ${HOME}/f; make -W fff.F90 OPTS=D fff; cd -
  ! cd ${HOME}/f; make -W fff.F90 fff; cd -
  ! cd ${HOME}/f; make OPTS=D fff; cd -

  ! Compilation with F90:
  ! f90 -c -g -e -DSUNMP -I./ -I${HOME}/include -I/contrib/include -o ${HOME}/obj/SUNMP/fff.o fff.F90
  ! f90 -o ${HOME}/bin/SUNMP/fff ${HOME}/obj/SUNMP/fff.o -lcsz -L/contrib/lib -lnetcdf -lsunmath -lnsl

  ! Test g95 and gfortran compilers:
  ! -P omits CPP line # information which confuses Fortran compilers
  ! -C preserves C/C++ comments and thus preserves // = Fortran concatenate operator
#if 0
  FC='gfortran'
  FC='g95'
  FC='ifc'
  source ${HOME}/bin/sh/cmp_chg.sh gcc ${FC};sudo ${HOME}/bin/sh/cmp_chg.sh gcc ${FC}
  ln -s -f ${DATA}/g95/g95-install /tmp/g95
  src_lst='dbg_mdl sng_mdl sng'
  src_lst='shr_kind_mod dbg_mdl sng_mdl utl_mdl xtr_mdl vec_mdl erf_mdl gmm_mdl mmr_mdl nf90_utl fff_mdl fff' # NB: Order is important
  cd ~/f
  for fl_stb in ${src_lst}; do
     /bin/rm -f ${fl_stb}.mod ${fl_stb}.o 
  done # end loop over fl_stb
  /bin/rm -f netcdf.* typesizes.*
  ln -s -f /usr/local/include/netcdf.mod.${FC} ~/f/netcdf.mod # G95 has trouble with module files in other directories
  ln -s -f /usr/local/include/typesizes.mod.${FC} ~/f/typesizes.mod
  for fl_stb in ${src_lst}; do
     printf "Compiling ${fl_stb}.F90...\n"
     ${FC} -Wp,"-P,-C,-DPRC_FLT,-DLINUX" -I ./ -I/usr/local/include -c -o ${fl_stb}.o -fno-second-underscore ${fl_stb}.F90
  done # end loop over fl_stb
  ${FC} -o sng dbg_mdl.o sng_mdl.o sng.o
  ${FC} -I ./ -I/usr/local/include -o fff ${MY_OBJ_DIR}/date_time.o dbg_mdl.o erf_mdl.o fff.o fff_mdl.o gmm_mdl.o mmr_mdl.o nf90_utl.o shr_kind_mod.o sng_mdl.o utl_mdl.o vec_mdl.o xtr_mdl.o -L/usr/local/lib -lnetcdf

  fl_stb='foo'
  cat dbg_mdl.F90 sng_mdl.F90 sng.F90 > ${fl_stb}.F90
  cat dbg_mdl.F90 shr_kind_mod.F90 sng_mdl.F90 utl_mdl.F90 vec_mdl.F90 xtr_mdl.F90 fff_mdl.F90 fff.F90 > ${fl_stb}.F90
  cpp -P -C -DPRC_FLT -DLINUX -I. -I/usr/local/include ${fl_stb}.F90 > ${fl_stb}.f90
  ${FC} ${fl_stb}.f90
#endif /* endif False */

  ! Usage:
  ! fff -f 1.0
  ! fff -Wl,-T -f 1.0 # Convert to/from big-endian with lf95 executables
  ! fff -f 1.0 -i /fs/cgd/home0/zender/tmp/in.nc
  ! fff -f 1.0 -i /gs/zender/nco/data/in.nc
  ! fff -f 1.0 --drc_in=${HOME}/nco/data --tst=gsl
  ! scp ~/f/fff.F90 dust.ess.uci.edu:f
  ! scp ~/f/fff.F90 givre.ess.uci.edu:f

  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use fff_mdl ! [mdl] Fortran features and regression tests
  use mmr_mdl,only:wrp_allocate,wrp_deallocate ! [mdl] Memory management
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use shr_kind_mod,only:r8=>shr_kind_r8,r4=>shr_kind_r4,DBLKIND ! [mdl] Precision r8, i8, ...
  use sng_mdl ! [mdl] String manipulation
  ! use typesizes ! [mdl] netCDF internal bytesizeok()
  use xtr_mdl ! [mdl] Extrapolation/interpolation handling
  use utl_mdl,only:date_time_get ! [mdl] Utility functions (date_time_get,mnt_chk...)

  implicit none

#if 0
  ! Hacks to try to solve PGI compiler problem
  interface
     integer function ftn_strcmp(sng1,sng2) ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
       character(len=*),intent(in)::sng1
       character(len=*),intent(in)::sng2
     end function ftn_strcmp
  end interface
#endif /* !0 */

  ! Parameters
  character(len=*),parameter::CVS_Date='$Date$' ! [sng] Date string
  character(len=*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(len=*),parameter::CVS_Name='$HeadURL$' ! [sng] File name string
  character(len=*),parameter::CVS_Revision='$Revision$' ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  character(len=*),parameter::sbr_nm='fff' ! [sng] Subroutine name
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer,parameter::bnd_nbr_max=100 ! [nbr] Maximum number of bands
  integer,parameter::dmn_nbr_max=3 ! [nbr] Maximum rank of variable
  integer,parameter::fl_in_unit=73 ! [idx] Input file unit
  integer,parameter::fl_out_unit=74 ! [idx] Output file unit
  integer,parameter::lat_nbr_max=64 ! [nbr] Maximum number of latitudes
  integer,parameter::len_sng_max=100 ! [nbr] Maximum length of string
  integer,parameter::lev_nbr_max=19 ! [nbr] Maximum number of levels
  integer,parameter::lon_nbr_max=128 ! [nbr] Maximum number of longitudes

  ! Set defaults for command-line options 
  character(80)::drc_in='/home/zender/nco/data'//nlc ! [sng] Input directory
  character(80)::drc_out='' ! [sng] Output directory
  character(80)::fl_in='in.nc'//nlc ! [sng] Input file
  character(80)::fl_out='foo.nc'//nlc ! [sng] Output file
  character(80)::sng_foo='Default value'//nlc ! [sng] String
  character(80)::tst_sng='No test specified'//nlc ! [sng] Name of test to perform
  integer::int_foo=1 ! [nbr] Integer
  logical::lgc_foo=.false. ! [flg] Logical
  real(r8)::cmp_prc_foo=0.0_r8 ! [frc] Computational precision
  real(r4)::flt_foo=0.0_r4 ! [frc] Float
  real(DBLKIND)::dbl_foo=0.0_DBLKIND ! [frc] Double

  ! Derived fields

  ! Locals with simple initialization and no command-line override
  ! integer::exit_status=0 ! [enm] Program exit status (non-standard Fortran)
  integer::rcd=nf90_noerr ! [rcd] Return success code
  integer::xtr_typ_LHS=xtr_prt_wgt+xtr_fll_nil+xtr_vrb ! [enm] Extrapolation type
  integer::xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil+xtr_vrb ! [enm] Extrapolation type

  character(2)::dsh_key ! [sng] command-line dash and switch
  character(200)::cmd_ln ! [sng] command-line
  character(200)::prg_ID ! [sng] Program ID
  character(26)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(80)::arg_val ! [sng] command-line argument value
  character(80)::opt_sng ! [sng] Option string

  ! Command-line parsing
  integer::arg_idx ! [idx] Counting index
  integer::arg_nbr ! [nbr] Number of command-line arguments
  integer::opt_lng ! [nbr] Length of option

  ! Model grid sizes
  integer::bnd_nbr ! [nbr] Number of bands
  integer::lat_idx ! [idx] Counting index for lat
  integer::lat_nbr ! [nbr] Number of latitudes
  integer::lev_idx ! [idx] Counting index for lev
  integer::lev_nbr ! [nbr] Number of levels
  integer::lon_idx ! [idx] Counting index for lon
  integer::lon_nbr ! [nbr] Number of longitudes

  ! netCDF bookkeeping
  integer::bnd_dmn_id ! [id] Dimension ID for bnd
  integer::cnt_lon_lev_lat(3) ! [nbr] Dimension sizes
  integer::dmn_lon_lev_lat(3) ! [id] Dimension IDs
  integer::lat_dmn_id ! [id] Dimension ID for lat
  integer::lev_dmn_id ! [id] Dimension ID for lev
  integer::lon_dmn_id ! [id] Dimension ID for lon
  integer::nc_id ! [id] netCDF file handle
  integer::lat_id ! [id] Coordinate ID
  integer::lev_id ! [id] Coordinate ID
  integer::lon_id ! [id] Coordinate ID
  integer::one_dmn_var_id ! [id] Variable ID
  integer::scalar_var_id ! [id] Variable ID
  integer::three_dmn_var_id ! [id] Variable ID
  integer::nf90_r8 ! [enm] External netCDF type for r8 kind
  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call

  ! Allocatable variables
  real(r8),pointer::lat(:) ! Coordinate variable
  real(r8),dimension(:),allocatable::lev ! Coordinate variable
  real(r8),dimension(:),allocatable::lon ! Coordinate variable
  real(r8),dimension(:),allocatable::one_dmn_var
  real(r8),dimension(:,:,:),allocatable::three_dmn_var
  real(r8),dimension(:,:,:),allocatable::three_dmn_var_crd
  real(r8),dimension(:),pointer::one_dmn_var_ptr

  ! Static variables
  real(r8) pi ! [frc] 3
  real(r8) scalar_var

  integer::mcdate(8)
  data mcdate /000101,011231,851201,850101,000312,000911,010312,010228/

  ! Main code
  dbg_lvl=0                 ! [idx] Causes DDD source window to display this file
  ! fxm: Get value of ${DATA} from environment
  ! drc_in=getenv('DATA')     ! [sng] Input directory
  ! drc_in=pxfgetenv('DATA')     ! [sng] Input directory
  ! write (6,'(a)') drc_in(1:ftn_strlen(drc_in))

  ! Locals requiring complex initialization expressions
  pi=4.0_r8*atan(1.0_r8) ! [frc] 3
  pi=0.0_r8+pi ! [frc] 3 CEWI

  ! Retrieve command-line arguments
  call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID
  call date_time_get(lcl_date_time) ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY  
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments
  arg_idx=1 ! [idx] Counting index
  loop_while_options: do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg(), increment arg_idx
     dsh_key=arg_val(1:2) ! [sng] First two characters of option
     if (dsh_key == '--') then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop 'Long option has no name'
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (dbg_lvl >= dbg_io) write (6,'(5a,i3)') prg_nm(1:ftn_strlen(prg_nm)), &
             ': DEBUG Double hyphen indicates multi-character option: ', &
             'opt_sng = ',opt_sng(1:ftn_strlen(opt_sng)),', opt_lng = ',opt_lng
        ! fxm: Change if else if construct to select case but how to handle fall-through cases elegantly?
        if (opt_sng == 'dbg' .or. opt_sng == 'dbg_lvl' ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == 'dbl' .or. opt_sng == 'dbl_foo' ) then
           call ftn_arg_get(arg_idx,arg_val,dbl_foo) ! [frc] Double
        else if (opt_sng == 'dfl' .or. opt_sng == 'deflate' ) then
           call ftn_arg_get(arg_idx,arg_val,dfl_lvl) ! [enm] Deflate level
        else if (opt_sng == 'drc_in') then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == 'drc_out') then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == 'fl_in') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
        else if (opt_sng == 'fl_out') then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
        else if (opt_sng == 'cmp_prc' .or. opt_sng == 'cmp_prc_foo' ) then
           call ftn_arg_get(arg_idx,arg_val,cmp_prc_foo) ! [frc] Computational precision
        else if (opt_sng == 'flt' .or. opt_sng == 'flt_foo' ) then
           call ftn_arg_get(arg_idx,arg_val,flt_foo) ! [frc] Float
        else if (opt_sng == 'int' .or. opt_sng == 'int_foo' ) then
           call ftn_arg_get(arg_idx,arg_val,int_foo) ! [nbr] Integer
        else if (opt_sng == 'lgc' .or. opt_sng == 'lgc_foo' ) then
           call ftn_arg_get(arg_idx,arg_val,lgc_foo) ! [lgc] Logical
        else if (opt_sng == 'tst' .or. opt_sng == 'tst_sng') then
           call ftn_arg_get(arg_idx,arg_val,tst_sng) ! [sng] Name of test to perform
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     else if (dsh_key(1:1) == '-') then ! endif long option
        ! Handle short options
        if (dsh_key == '-3') then
           fl_out_fmt=nf90_format_classic ! [enm] Output file format
        else if (dsh_key == '-4') then
           fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
        else if (dsh_key == '-D') then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (dsh_key == '-f') then
           call ftn_arg_get(arg_idx,arg_val,flt_foo) ! [frc] Float
        else if (dsh_key == '-i') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
        else if (dsh_key == '-l') then
           call ftn_arg_get(arg_idx,arg_val,int_foo)
        else if (dsh_key == '-L') then
           call ftn_arg_get(arg_idx,arg_val,dfl_lvl) ! [enm] Deflate level
        else if (dsh_key == '-p') then
           lgc_foo=.not.lgc_foo
        else if (dsh_key == '-o') then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
        else if (dsh_key == '-r') then
           call ftn_arg_get(arg_idx,arg_val,xtr_typ_RHS) ! [enm] Extrapolation type
        else if (dsh_key == '-s') then
           call ftn_arg_get(arg_idx,arg_val,sng_foo) ! [sng] String
        else if (dsh_key == '-v') then
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

  ! Check for netCDF Fortran90 interface conformance
  !if(.not.byteSizesOK()) then ! Requires use typesizes.mod explicitly
  !   stop 'Compiler does not support required kinds of variables'
  !end if ! endif bytesizesok()

  ! Read in netCDF data
  rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
  if(rcd /= nf90_noerr) then
     write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': HINT Invoke with:'
     write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)), &
          ' --dbg=1 --drc_in=${HOME}/nco/data --tst=foo'
     call nf90_err_exit(rcd,'Unable to open '//fl_in(1:ftn_strlen(fl_in))//'in '//__FILE__)
  endif ! endif err
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lat',lat_dmn_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dmn_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dmn_id)
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr),sbr_nm//': inquire_dim bnd')
  if(bnd_nbr > bnd_nbr_max) stop 'bnd_nbr > bnd_nbr_max'
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr),sbr_nm//': inquire_dim lat')
  if(lat_nbr > lat_nbr_max) stop 'lat_nbr > lat_nbr_max'
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_dmn_id,len=lev_nbr),sbr_nm//': inquire_dim lev')
  if(lev_nbr > lev_nbr_max) stop 'lev_nbr > lev_nbr_max'
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr),sbr_nm//': inquire_dim lon')
  if(lon_nbr > lon_nbr_max) stop 'lon_nbr > lon_nbr_max'
  if(rcd /= nf90_noerr) stop 'Error retrieving dimension sizes'

  ! Allocate space for dynamic arrays
  ! fxm: TODO 13
  ! call wrp_allocate(lat_nbr,lat,'lat') ! Coordinate variable
  allocate(lat(lat_nbr),stat=rcd) ! Coordinate variable
  if(rcd /= 0) stop 'allocate() failed for lat'
  allocate(lev(lev_nbr),stat=rcd) ! Coordinate variable
  if(rcd /= 0) stop 'allocate() failed for lev'
  allocate(lon(lon_nbr),stat=rcd) ! Coordinate variable
  if(rcd /= 0) stop 'allocate() failed for lon'
  allocate(one_dmn_var(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop 'allocate() failed for one_dmn_var'
  allocate(three_dmn_var(lon_nbr,lev_nbr,lat_nbr),stat=rcd)
  if(rcd /= 0) stop 'allocate() failed for three_dmn_var'
  allocate(three_dmn_var_crd(lon_nbr,lat_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop 'allocate() failed for three_dmn_var_crd'
  ! Assemble count vectors for each multidimensional combination of dimensions
  cnt_lon_lev_lat=(/lon_nbr,lev_nbr,lat_nbr/)
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'lat',lat_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lev',lev_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lon',lon_id)
  rcd=nf90_wrp_inq_varid(nc_id,'scalar_var',scalar_var_id)
  rcd=nf90_wrp_inq_varid(nc_id,'one_dmn_var',one_dmn_var_id)
  rcd=nf90_wrp_inq_varid(nc_id,'three_dmn_var',three_dmn_var_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,lat_id,lat),sbr_nm//': get_var lat')
  rcd=nf90_wrp(nf90_get_var(nc_id,lev_id,lev),sbr_nm//': get_var lev')
  rcd=nf90_wrp(nf90_get_var(nc_id,lon_id,lon),sbr_nm//': get_var lon')
  rcd=nf90_wrp(nf90_get_var(nc_id,scalar_var_id,scalar_var),sbr_nm//': get_var scalar_var')
  rcd=nf90_wrp(nf90_get_var(nc_id,one_dmn_var_id,one_dmn_var),sbr_nm//': get_var one_dmn_var')
  rcd=nf90_wrp(nf90_get_var(nc_id,three_dmn_var_id,three_dmn_var),sbr_nm//': get_var three_dmn_var')
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_in,'Ingested') ! [fnc] Close file

  ! Main code goes here
  if (dbg_lvl == dbg_off) write(6,'(3a)') 'HINT: Invoke with ', &
       prg_nm(1:ftn_strlen(prg_nm)), &
       ' --tst=tst_sng to see specific test or --dbg=3 to see current test'

  if (.false.) then
     write (6,'(a,i2)') 'ftn_strlen(tst_sng) = ',ftn_strlen(tst_sng)
     write (6,'(a,i2)') 'ftn_strlen(''prn'') = ',ftn_strlen('prn')
     write (6,'(2a)') 'tst_sng(1:ftn_strlen(tst_sng)) = ',tst_sng(1:ftn_strlen(tst_sng))
     write (6,'(2a)') 'tst_sng(1:len_trim(tst_sng)) = ',tst_sng(1:len_trim(tst_sng))
     call ftn_strprn(tst_sng) ! [fnc] Print character values of string
     call ftn_strprn('prn') ! [fnc] Print character values of string
     ! Avoid recursive I/O by storing then printing result of routines that may themselves print
     int_foo=ftn_strcmp(tst_sng,'prn')
     write (6,'(a,i2)') 'ftn_strcmp(tst_sng,''prn'') = ',int_foo
  endif ! endif .false.

  if (.true.) then
     write (6,*) 'main() reports cmp_prc_foo = ',cmp_prc_foo
     write (6,*) 'main() reports dbl_foo = ',dbl_foo
     write (6,*) 'main() reports flt_foo = ',flt_foo
     write (6,*) 'main() reports int_foo = ',int_foo
     write (6,*) 'main() reports lgc_foo = ',lgc_foo
     write (6,*) 'main() reports sng_foo = ',sng_foo
  endif ! endif .true.

  ! Test date and time routines
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'gmt')==0) call date_time_tst(lcl_date_time)

  ! Test C-interface time routines
  ! NB: argument MUST be equivalent to double precision C (8 bytes usually)
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'gmt')==0) call gmtime_tst(dbl_foo)

  ! Test Scientific library functions
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'gsl')==0) call gsl_tst(cmp_prc_foo)

  ! Test IEEE error trapping
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'ieee')==0) call ieee_tst(cmp_prc_foo)

  ! Test rbn_vec() subroutine and subsidiaries
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'ntp')==0) call rbn_tst(xtr_typ_LHS,xtr_typ_RHS)

  ! Test OpenMP routines
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'omp')==0) call omp_tst(cmp_prc_foo,int_foo)

  ! Test precision
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'prc')==0) call prc_tst(cmp_prc_foo,dbl_foo,flt_foo,int_foo)

  ! Test formatted print statements
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'prn')==0) call prn_tst(cmp_prc_foo)

  ! Pointers
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'ptr')==0) call ptr_tst(one_dmn_var,one_dmn_var_ptr)

  ! Test Fortran syntax
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'syn')==0) call syn_tst(cmp_prc_foo)

  ! Test netCDF I/O
  if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'nc')==0 .or. ftn_strcmp(tst_sng,'nco')==0) then
     write (6,'(a)') 'Printing contents of input file...'
     do lat_idx=1,lat_nbr
        write (6,'(a4,i3,a2,f10.5)') 'lat(',lat_idx,')=',lat(lat_idx)
     enddo ! end loop over lat
     do lev_idx=1,lev_nbr
        write (6,'(a4,i3,a2,f10.5)') 'lev(',lev_idx,')=',lev(lev_idx)
     enddo ! end loop over lev
     do lon_idx=1,lon_nbr
        write (6,'(a4,i3,a2,f10.5)') 'lon(',lon_idx,')=',lon(lon_idx)
     enddo ! end loop over lon
     write (6,'(a)') 'three_dmn_var in CCM order: (lon,lev,lat)'
     do lat_idx=1,lat_nbr
        do lev_idx=1,lev_nbr
           do lon_idx=1,lon_nbr
              write (6,'(a18,3(i3,a5),a1,f9.5)') &
                   'three_dmn_var(lon=',lon_idx,',lev=',lev_idx,',lat=',lat_idx,')','=',three_dmn_var(lon_idx,lev_idx,lat_idx)
           enddo ! end loop over lat
        enddo ! end loop over lev
     enddo ! end loop over lon
     write (6,'(a)') 'three_dmn_var in COORDS order: (lon,lat,lev)'
     do lev_idx=1,lev_nbr
        do lat_idx=1,lat_nbr
           do lon_idx=1,lon_nbr
              write (6,'(a18,3(i3,a5),a1,f9.5)') &
                   'three_dmn_var(lon=',lon_idx,',lev=',lev_idx,',lat=',lat_idx,')','=',three_dmn_var(lon_idx,lev_idx,lat_idx)
           enddo ! end loop over lon
        enddo ! end loop over lat
     enddo ! end loop over lev
     ! ncks -C -F -H -v three_dmn_var,three_dmn_var_crd ${HOME}/nco/data/in.nc
     write (6,'(a)') 'three_dmn_var_crd=reshape(three_dmn_var,shape=(/lon_nbr,lat_nbr,lev_nbr/),order=(/1,3,2/))'
     three_dmn_var_crd=reshape(three_dmn_var,shape=(/lon_nbr,lat_nbr,lev_nbr/),order=(/1,3,2/))
     do lev_idx=1,lev_nbr
        do lat_idx=1,lat_nbr
           do lon_idx=1,lon_nbr
              write (6,'(a22,3(i3,a5),a1,f9.5)') &
                   'three_dmn_var_crd(lon=',lon_idx,',lat=',lat_idx,',lev=',lev_idx,')','=', &
                   three_dmn_var_crd(lon_idx,lat_idx,lev_idx)
           enddo ! end loop over lon
        enddo ! end loop over lat
     enddo ! end loop over lev
     
     ! Test array manipulation
     ! if (dbg_lvl==dbg_old .or. ftn_strcmp(tst_sng,'arr')==0) call arr_tst(three_dim_var)

  endif ! endif tst_sng==nc
  
  ! Begin netCDF output routines
#ifdef ENABLE_NETCDF4
  if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
  if (fl_out_fmt == nf90_format_64bit) then
     nf90_create_mode=nf90_create_mode+nf90_64bit_offset
  else if (fl_out_fmt == nf90_format_netcdf4) then
     nf90_create_mode=nf90_create_mode+nf90_netcdf4
  else if (fl_out_fmt == nf90_format_netcdf4_classic) then
     nf90_create_mode=nf90_create_mode+(nf90_classic_model+nf90_netcdf4)
  end if ! end else fl_out_fmt
#else /* !ENABLE_NETCDF4 */
  if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_classic ! [enm] Output file format
  if(fl_out_fmt == nf90_format_classic) nf90_create_mode=nf90_create_mode+0 ! CEWI
  dfl_lvl=dfl_lvl+0 ! CEWI
  flg_dfl=flg_dfl+0 ! CEWI
  flg_shf=flg_shf+0 ! CEWI
#endif /* !ENABLE_NETCDF4 */
  rcd=nf90_wrp_create(fl_out,nf90_create_mode,nc_id,sbr_nm=sbr_nm)
  ! Define dimension IDs
  rcd=nf90_wrp(nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dmn_id),sbr_nm//': def_dim bnd in '//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'lat',lat_nbr,lat_dmn_id),sbr_nm//': def_dim lat in '//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'lev',lev_nbr,lev_dmn_id),sbr_nm//': def_dim lev in '//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'lon',lon_nbr,lon_dmn_id),sbr_nm//': def_dim lon in '//__FILE__)
  ! Assemble ID and count vectors for each multidimensional combination of dimensions
  dmn_lon_lev_lat=(/lon_dmn_id,lev_dmn_id,lat_dmn_id/)
  cnt_lon_lev_lat=(/lon_nbr,lev_nbr,lat_nbr/)
  cnt_lon_lev_lat=0+cnt_lon_lev_lat ! CEWI
  ! Variable definitions
  nf90_r8=nf90_xtype_r8_get() ! [enm] External netCDF type for r8 kind
  ! Wrapped
  rcd=nf90_wrp(nf90_def_var(nc_id,'three_dmn_var',nf90_r8,dmn_lon_lev_lat,three_dmn_var_id), &
       sbr_nm//': def_var three_dmn_var in '//__FILE__)
  ! Set HDF Lempel-Ziv compression level, if requested
#ifdef ENABLE_NETCDF4
  if (dfl_lvl > 0) rcd=nf90_def_var_deflate(nc_id,three_dmn_var_id,flg_shf,flg_dfl,dfl_lvl)
#endif /* !ENABLE_NETCDF4 */

  ! Not wrapped
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat',nf90_r8,lat_dmn_id,lat_id),sbr_nm//': def_var lat in '//__FILE__)
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev',nf90_r8,lev_dmn_id,lev_id),sbr_nm//': def_var lev in '//__FILE__)
  rcd=nf90_wrp(nf90_def_var(nc_id,'lon',nf90_r8,lon_dmn_id,lon_id),sbr_nm//': def_var lon in '//__FILE__)
  rcd=nf90_wrp(nf90_def_var(nc_id,'one_dmn_var',nf90_r8,bnd_dmn_id,one_dmn_var_id),sbr_nm//': def_var one_dmn_var in '//__FILE__)
  rcd=nf90_wrp(nf90_def_var(nc_id,'scalar_var',nf90_r8,scalar_var_id),sbr_nm//': def_var scalar_var in '//__FILE__)
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': put_att CVS_Id in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm//': put_att creation_date in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'history',lcl_date_time//': '// &
       cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm//': put_att history in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm//': put_att prg_ID in '//__FILE__)
  ! Add English text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'long_name','Latitude'),sbr_nm//': put_att lat in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'long_name','Layer pressure'),sbr_nm//': put_att lev in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'long_name','Longitude'),sbr_nm//': put_att lon in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,one_dmn_var_id,'long_name','Description'),sbr_nm//': put_att one_dmn_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,scalar_var_id,'long_name','Description'),sbr_nm//': put_att scalar_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,three_dmn_var_id,'long_name','Description'),sbr_nm//': put_att three_dmn_var in '//__FILE__)
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'units','degree'),sbr_nm//': put_att lat in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'units','pascal'),sbr_nm//': put_att lev in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'units','degree'),sbr_nm//': put_att lon in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,one_dmn_var_id,'units','unknown'),sbr_nm//': put_att one_dmn_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,scalar_var_id,'units','unknown'),sbr_nm//': put_att scalar_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,three_dmn_var_id,'units','unknown'),sbr_nm//': put_att three_dmn_var in '//__FILE__)
  ! Now that all dimensions, variables, and attributes have been defined, make call to end define mode
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_id,lat),sbr_nm//': put_var lat in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_id,lev),sbr_nm//': put_var lev in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lon_id,lon),sbr_nm//': put_var lon in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,one_dmn_var_id,one_dmn_var),sbr_nm//': put_var one_dmn_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,scalar_var_id,scalar_var),sbr_nm//': put_var scalar_var in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,three_dmn_var_id,three_dmn_var),sbr_nm//': put_var three_dmn_var in '//__FILE__)
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to')

  ! De-allocate dynamic variables
  call wrp_deallocate(lat,rcd,'lat in '//__FILE__) ! Coordinate variable
  if (allocated(lev)) deallocate(lev,stat=rcd) ! Coordinate variable
  if(rcd /= 0) stop 'deallocate() failed for lev'
  if (allocated(lon)) deallocate(lon,stat=rcd) ! Coordinate variable
  if(rcd /= 0) stop 'deallocate() failed for lon'
  if (allocated(one_dmn_var)) deallocate(one_dmn_var,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for one_dmn_var'
  if (allocated(three_dmn_var)) deallocate(three_dmn_var,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for three_dmn_var'

#ifdef CPP_FOO
  write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': CPP token CPP_FOO is defined'
#else /* !CPP_FOO */
  write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': CPP token CPP_FOO is not defined'
#endif /* !CPP_FOO */

1000 continue ! Jumping point for quick exit

  ! call exit(exit_status)    ! [enm] Exit with current exit status (non-standard Fortran)
end program fff

