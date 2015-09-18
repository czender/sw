# $Id$

# Purpose: Makefile for shell scripts module sh

# Shell scripts are installed here
MY_BIN_DIR := ${HOME}/bin/sh

# Bin targets
ifneq (,$(findstring .uci.edu,${HOST}))
 MDL_BIN_TRG := msrcp msread mswrite rdsp
else
 MDL_BIN_TRG := 
endif
#MDL_BIN_TRG := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_TRG})
MDL_BIN_SYM_LNK := $(basename ${MDL_BIN_TRG}) # Symbolic links
MDL_BIN_STB := ${MDL_BIN_SYM_LNK} # All files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Rule for installing shell script file
# % : %.sh 
# 	ln -s $(notdir $<) $@ 
% : %.pl
	ln -s $(notdir $<) $@ 
msread : msrcp.pl
	ln -s $(notdir $<) $@ 
mswrite : msrcp.pl
	ln -s $(notdir $<) $@ 
rdsp : msrcp.pl
	ln -s $(notdir $<) $@ 

# Default targets
all: bin
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln dbg 
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .sh .pl .csh
# Delete targets which were not successfully made
.DELETE_ON_ERROR:

# Targets in bin
bin: ${MDL_BIN_TRG}
bin_cln :
	rm -f ${MDL_BIN}
# Targets in bld
buildclean: bld_cln
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	printf "MDL_BIN = ${MDL_BIN}\n"
	printf "MDL_BIN_TRG = ${MDL_BIN_TRG}\n"
	printf "MDL_BIN_SYM_LNK = ${MDL_BIN_SYM_LNK}\n"
	printf "MDL_BIN_STB = ${MDL_BIN_STB}\n"
	printf "FOO = ${FOO}\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : bin 
tags : 
	etags * Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln
# Targets in dat
data:	dat
dat : 
dat_cln :
# Targets in dpn
depend: dpn
dpn : 
dpn_cln :
# Targets in inc
include: inc
inc : ${MDL_INC}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
library: lib
lib : 
libclean: lib_cln
lib_cln :
# Targets in obj
object: obj
obj : 
objclean: obj_cln
obj_cln :

