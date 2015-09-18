# $Id$

# Purpose: Perl command line scripts to do global text substitutions on all files in specified directories

# -pi conjoins two separate Perl switches, i.e., same as -p -i
# -p causes Perl script to loop over filename arguments, and to print each input line
# -i causes Perl to edit the files in-place
# -e specifies that following argument is Perl command

# Change CVS keywords to SVN keywords
perl -pi -e 's/\$Header\$/\$Id\$/g;' ~/aeroce/README

# Change UNIX LF to Windows CR/LF
perl -pi -e 's/\012/\012\015/g;' ~/job/cv.txt

# Change location of CaCO3 data
perl -pi -e 's/\/CaCO3/rsmas/g;' `/bin/ls *`

# Change NCO code
perl -pi -e 's///g;' `/bin/ls *.c`

# Change paper title containing single quote
perl -pi -e 's/1990.s/1990s/g;' `/bin/ls *.tex`

# Convert physics course to boundary layer course
perl -pi -e 's/phz/bnd/g;' `/bin/ls *`

# Convert #endif /* __cplusplus */ to #endif /* !__cplusplus */
perl -pi -e 's/#endif \/\* __cplusplus/#endif \/* !__cplusplus/g;' `/bin/ls *.h`

# Convert <INSTALLDIR> to /opt/intel_cc_80
perl -pi -e 's/<INSTALLDIR>/\/opt\/intel_cc_80/g;' i[cf]c* i[cf]cvars* icpc* 
# Convert <INSTALLDIR> to /opt/intel_fc_80
perl -pi -e 's/<INSTALLDIR>/\/opt\/intel_fc_80/g;' i[cf]c* i[cf]cvars* icpc* 
# Convert <INSTALLDIR> to /opt/intel_idb_80
perl -pi -e 's/<INSTALLDIR>/\/opt\/intel_idb_80/g;' i[cf]c* i[cf]cvars* icpc* 

# Convert machine domains
perl -pi -e 's/.ps.uci.edu/.ess.uci.edu/g;' `/bin/ls *.tex`

# Convert UCI class years
perl -pi -e 's/02s/03s/g;s/s02/s03/g;s/S02/S03/g;' `/bin/ls *.tex`

# Convert $HOME to ${HOME}
perl -pi -e 's/\$HOME/\$\{HOME\}/g;s/\$DATA/\$\{DATA\}/g;' `/bin/ls *.tex`

# Fortran77 -> Fortran90 conversion
perl -pi -e 's/enddo/end do/g' `/bin/ls *.F *.F90`
perl -pi -e 's/character*8 /character(8)/g;s/character fl_out*80 /character(80) fl_out/g;s/character*80 /character\(80\)/g;s/character\*(\*)/character(len=*)/g;' `/bin/ls *.F *.F90 *.h *.com`

# Generate netCDF wrappers
perl -pi -e 's/rcd=rcd+nf90_def_dim(nc_id,\'lev\',lev_nbr,lev_dim_id)//g' `/bin/ls *.F90`

# Miscellaneous rules
perl -pi -e 's/dbg\.h/dbg.com/g;s/xtr\.h/xtr.com/g;s/ini\.F/ini.com/g;s/srt\.F/srt.com/g;s/end\.F/end.com/g;s/\.inc/.com/g;s/netcdf\.com/netcdf.inc/g;' `/bin/ls *.F *.com`
perl -pi -e 's/rcd\+nf90_inq_dimid/nf90_wrp_inq_dimid/g;s/rcd\+nf90_inq_varid/nf90_wrp_inq_varid/g;' `/bin/ls *.F90`
perl -pi -e 's/flamenco/heinlein/g' `/bin/ls *`
perl -pi -e 's/ZBN02/ZBN03/g;s/dst_mdl/ppr_ZBN03/g' `/bin/ls *`
perl -pi -e 's/nccat/ncrcat/g' `/bin/ls *`
perl -pi -e 's/max_nbr_mlc_htrn/mlc_nbr_max_htrn/g;s/max_nbr_iso_htrn/iso_nbr_max_htrn/g;s/max_nbr_iso_per_mlc_htrn/iso_per_mlc_nbr_max_htrn/g;' `/bin/ls *`
perl -pi -e 's/\:LT/':LT/g' `/bin/ls *.sh`
perl -pi -e 's/1.E36...LT/1.E35,':LT'/g' `/bin/ls *.sh`
perl -pi -e 's/erbe_b/erbe_a/g' `/bin/ls *.sh`
perl -pi -e 's/\.jul/_jul/g' `/bin/ls *.sh`
perl -pi -e 's/\.out/.txt/g' `/bin/ls *.sh`
perl -pi -e 's/adv_psd_01/psd_01_adv/g;s/adv_prg_anl_01/prg_anl_01_adv/g;s/adv_sdn_eos/sdn_adv_eos/g;s/~\/ess/~\/hire/' `/bin/ls *`
perl -pi -e 's/hmwskpshrt/hmwskpsht/g' `/bin/ls *.tex`
perl -pi -e 's/strstr/ftn_strstr/g;s/strlen/ftn_strlen/g;s/strcpy/ftn_strcpy/g;s/strnulc/ftn_strnulc/g;s/strlsc/ftn_strlsc/g;s/strfic/ftn_strfic/g;s/strcat/ftn_strcat/g;s/strprn/ftn_strprn/g;s/ strcpylsc/ ftn_strcpylsc/g;s/strpfx/ftn_strpfx/g;s/strini/ftn_strini/g;s/drcpfx/ftn_drcpfx/g;s/ strnul/ ftn_strnul/g;s/prg_ID_mk/ftn_prg_ID_mk/g;s/cmd_ln_sng/ftn_cmd_ln_sng/g;s/date2sng/ftn_date2sng/g;s/sec2sng/ftn_sec2sng/g;s/getarg_wrp/ftn_getarg_wrp/g;s/getarg_err/ftn_getarg_err/g;s/sng_arg_get/ftn_sng_arg_get/g;s/flt_arg_get/ftn_flt_arg_get/g;s/int_arg_get/ftn_int_arg_get/g;' `/bin/ls *.F *.F90`
perl -pi -e 's/\$\(MY_OBJ_DIR\)/\$\{MY_OBJ_DIR\}/g;s/\$\(MY_BIN_DIR\)/\$\{MY_BIN_DIR\}/g;s/\$\(HOME\)/\$\{HOME\}/g;s/\$\(MY_SHR_DIR\)/\$\{MY_SHR_DIR\}/g;s/\$\(GSL_LIB\)/\$\{GSL_LIB\}/g;s/\$\(MY_DPN_DIR\)/\$\{MY_DPN_DIR\}/g;s/\$\(CC\)/\$\{CC\}/g;s/\$\(CPPFLAGS\)/\$\{CPPFLAGS\}/g;s/\$\(CFLAGS\)/\$\{CFLAGS\}/g;s/\$\(FC\)/\$\{FC\}/g;s/\$\(FFLAGS\)/\$\{FFLAGS\}/g;s/\$\(MY_INC_DIR\)/\$\{MY_INC_DIR\}/g;s/\$\(CPP\)/\$\{CPP\}/g;s/\$\(NETCDF_INC\)/\$\{NETCDF_INC\}/g;s/\$\(NETCDF_LIB\)/\$\{NETCDF_LIB\}/g;s/\$\(PVM_ARCH\)/\$\{PVM_ARCH\}/g;s/\$\(C\+\+FLAGS\)/\$\{C\+\+FLAGS\}/g;s/\$\(OPTS\)/\$\{OPTS\}/g;s/\$\(MY_LIB_DIR\)/\$\{MY_LIB_DIR\}/g;s/\$\(SGI_ABI_FLG\)/\$\{SGI_ABI_FLG\}/g;s/\$\(GCC_ABI_FLG\)/\$\{GCC_ABI_FLG\}/g;s/\$\(CPP_PTH\)/\$\{CPP_PTH\}/g;s/\$\(LDFLAGS\)/\$\{LDFLAGS\}/g;s/\$\(UNAMES\)/\$\{UNAMES\}/g;s/\$\(MAKE\)/\$\{MAKE\}/g;s/\$\(C\+\+\)/\$\{C\+\+\}/g;s/\$\(DATA\)/\$\{DATA\}/g;s/\$\(MY_BLD_DIR\)/\$\{MY_BLD_DIR\}/g;s/\$\(OMP\)/\$\{OMP\}/g;s/\$\(MY_DAT_DIR\)/\$\{MY_DAT_DIR\}/g;s/\$\(MDL_INC\)/\$\{MDL_INC\}/g;s/\$\(MDL_DPN\)/\$\{MDL_DPN\}/g;s/\$\(MDL_PTH\)/\$\{MDL_PTH\}/g;s/\$\(MDL_BIN_TRG\)/\$\{MDL_BIN_TRG\}/g;s/\$\(MDL_BIN_SYM_LNK\)/\$\{MDL_BIN_SYM_LNK\}/g;' \
perl -pi -e 's/\$\(MDL_BIN_STB\)/\$\{MDL_BIN_STB\}/g;s/\$\(MDL_INC_TRG\)/\$\{MDL_INC_TRG\}/g;s/\$\(MDL_INC_SYM_LNK\)/\$\{MDL_INC_SYM_LNK\}/g;s/\$\(DPN_GNR\)/\$\{DPN_GNR\}/g;s/\$\(null\)/\$\{null\}/g;s/\$\(slash\)/\$\{slash\}/g;s/\$\(slash_rx\)/\$\{slash_rx\}/g;s/\$\(MDL_SRC\)/\$\{MDL_SRC\}/g;s/\$\(dir\)/\$\{dir\}/g;s/\$\(MY_OBJ_DIR_RX\)/\$\{MY_OBJ_DIR_RX\}/g;s/\$\(MY_DPN_DIR_RX\)/\$\{MY_DPN_DIR_RX\}/g;s/\$\(dollar\)/\$\{dollar\}/g;s/\$\(MDL_OBJ\)/\$\{MDL_OBJ\}/g;s/\$\(HOSTNAME\)/\$\{HOSTNAME\}/g;s/\$\(HOST\)/\$\{HOST\}/g;s/\$\(MY_DOC_DIR\)/\$\{MY_DOC_DIR\}/g;s/\$\(MY_ES_DIR\)/\$\{MY_ES_DIR\}/g;s/\$\(VPATH\)/\$\{VPATH\}/g;s/\$\(FIXEDFLAGS\)/\$\{FIXEDFLAGS\}/g;s/\$\(FREEFLAGS\)/\$\{FREEFLAGS\}/g;s/\$\(MAKECMDGOALS\)/\$\{MAKECMDGOALS\}/g;s/\$\(GOALS_WHICH_DELETE_DEPENDENCY_FILES\)/\$\{GOALS_WHICH_DELETE_DEPENDENCY_FILES\}/g;s/\$\(INCLUDE_DPN\)/\$\{INCLUDE_DPN\}/g;s/\$\(MDL_BIN\)/\$\{MDL_BIN\}/g;s/\$\(space\)/\$\{space\}/g;s/\$\(MDL_INC_STB\)/\$\{MDL_INC_STB\}/g;s/\$\(F90FLAGS\)/\$\{F90FLAGS\}/g;s/\$\(CCM_TRG_FLG\)/\$\{CCM_TRG_FLG\}/g;s/\$\(CCM_RSN\)/\$\{CCM_RSN\}/g;s/\$\(DMN_SPC_TKN\)/\$\{DMN_SPC_TKN\}/g;s/\$\(MY_SRC_DIR\)/\$\{MY_SRC_DIR\}/g;s/\$\(FOO\)/\$\{FOO\}/g;' \
perl -pi -e 's/\$\(CCM_TRG_FLG\)/\$\{CCM_TRG_FLG\}/g;s/\$\(CCM_RSN\)/\$\{CCM_RSN\}/g;s/\$\(DMN_SPC_TKN\)/\$\{DMN_SPC_TKN\}/g;s/\$\(MY_SRC_DIR\)/\$\{MY_SRC_DIR\}/g;s/\$\(FOO\)/\$\{FOO\}/g;' \
${HOME}/aca/Makefile \
${HOME}/arese/Makefile \
${HOME}/avhrr/Makefile \
${HOME}/bin/sh/Makefile \
${HOME}/c++/Makefile \
${HOME}/c/Makefile \
${HOME}/ccm/Makefile \
${HOME}/ccm_o2x/Makefile \
${HOME}/ck/Makefile \
${HOME}/crm/crm/bld/Makefile \
${HOME}/crr/Makefile \
${HOME}/dot/Makefile \
${HOME}/dst/Makefile \
${HOME}/erbe/Makefile \
${HOME}/ess/Makefile \
${HOME}/ess_atm/Makefile \
${HOME}/ess_rdn/Makefile \
${HOME}/f/Makefile \
${HOME}/frc/Makefile \
${HOME}/fsf/Makefile \
${HOME}/hdf/Makefile \
${HOME}/hire/Makefile \
${HOME}/idl/Makefile \
${HOME}/idx_rfr/Makefile \
${HOME}/igpp/Makefile \
${HOME}/job/Makefile \
${HOME}/jrn/Makefile \
${HOME}/linux/Makefile \
${HOME}/ltr/Makefile \
${HOME}/map/Makefile \
${HOME}/match/Makefile \
${HOME}/match_dst/dst/Makefile \
${HOME}/mie/Makefile \
${HOME}/mny/Makefile \
${HOME}/nco/bld/Makefile \
${HOME}/poem/Makefile \
${HOME}/prp/Makefile \
${HOME}/rvw/Makefile \
${HOME}/seawifs/Makefile \
${HOME}/slr_spc/Makefile \
${HOME}/tex/Makefile \
${HOME}/time/Makefile \
${HOME}/www/Makefile
