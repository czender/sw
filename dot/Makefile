# $Id$

# Purpose: Makefile for system files module dot

# System files are installed here
MY_DOT_DIR := ${HOME}
MY_SSH_DIR := ${HOME}/.ssh

# Variables having to do with header files created by module
MDL_DOT_TRG := ChangeLog cshrc mailrc README fvwm2rc rhosts login xinitrc.old signature crontab emacs mwmrc twmrc fvwmrc hluresfile idl idle latex2html-init lightning lightningmask logout mailcap whale xmodmaprc profile.old fvwm95rc plan gnus bashrc bash_profile bash_logout dbxrc gdbinit pgdbgrc cvDebugger Xdefaults cvsignore inputrc forward netrc mpd.conf dupload.conf rpmmacros pbuilderrc wgetrc dput.cf gitconfig soft condarc
# 20161222 deprecate daprc and dodsrc until working
# Dot targets
MDL_DOT_SYM_LNK := # Symbolic links
MDL_DOT_STB := ${MDL_DOT_TRG} ${MDL_DOT_SYM_LNK} # All header files in MY_DOT_DIR
MDL_DOT := $(addprefix ${MY_DOT_DIR}/.,${MDL_DOT_STB}) # dst_cln removes these files

# Rule for installing system dot file is prepend a dot and copy to home directory
${MY_DOT_DIR}/.% : %
	cp -f -p $(notdir $<) $@
${MY_SSH_DIR}/% : %
	cp -f -p $(notdir $<) $@
${MY_DOT_DIR}/.Xdefaults : Xdefaults
	cp -f -p $(notdir $<) $@
# Running xrdb to update the X resource database requires X-server access
# This may not be possible when a remote machine has no X connection
# Thus error handling should be implemented in this rule
# One method would be to check whether xhost returns an error
# If not then go ahead and run xrdb
#	xrdb -merge $@

# Default targets
all: dot ssh

# Targets in bin
bin : 
bin_cln :
# Targets in bld
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	printf "MDL_DOT = ${MDL_DOT}\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : dot 
tags : 
	etags * Makefile
uninstall : dst_cln bin_cln dot_cln lib_cln
# Targets in dat
dat : 
dat_cln :
# Targets in dpn
dpn : 
dpn_cln :
# Targets in dot
dot : ${MDL_DOT}
dot_cln :
	rm -f ${MDL_DOT}
fs : 
#	cp -f -p $(notdir $<) /fs/cgd/home0/zender/.$<
	for fl in ${MDL_DOT_TRG}; do cp -f -p $$fl /fs/cgd/home0/zender/.$$fl; done
tempest : ${MDL_DOT}
	for fl in ${MDL_DOT_TRG}; do rcp -p $$fl tempest.ucar.edu:/home/tempest/zender/.$$fl; done
ouray : ${MDL_DOT}
	for fl in ${MDL_DOT_TRG}; do rcp -p $$fl ouray.ucar.edu:/home0/zender/.$$fl; done
acd : ${MDL_DOT}
	rsh garcia.acd.ucar.edu "mv .fvwm2rc .fvwm2rc.cgd; mv .xinitrc .xinitrc.cgd"
	rsh acd.ucar.edu "mv .fvwm2rc .fvwm2rc.cgd; mv .xinitrc .xinitrc.cgd"
# Targets in lib
lib : 
lib_cln :
# Targets in obj
obj : 
obj_cln :
# Targets in ssh
# $^ is name of all dependencies = foo.c bar.c
ssh : ${MY_SSH_DIR}/authorized_keys
	-cp -f $(notdir $^) ${MY_SSH_DIR}
	-scp -p $^ zender@spectral.ess.uci.edu:.ssh
	-scp -p $^ zender@e3sm.ess.uci.edu:.ssh
	-scp -p $^ zender@hpc3.rcic.uci.edu:.ssh
	-scp -p $^ zender@imua.ess.uci.edu:.ssh
	-scp -p $^ zender@maluhia.ess.uci.edu:.ssh
	-scp -p $^ zender@sastrugi.ess.uci.edu:.ssh
	-scp -p $^ zender@home.ps.uci.edu:.ssh
	-scp -p $^ ac.zender@blues.lcrc.anl.gov:.ssh
#	-scp -p $^ zender@compy.pnl.gov:.ssh
#	-scp -p $^ zender@katabatic.ess.uci.edu:.ssh
#	-scp -p $^ zender@grele.ess.uci.edu:.ssh
#	-scp -p $^ zender@gplogin2.ps.uci.edu:.ssh
#	-scp -p $^ zender@babyblue.ucar.edu:.ssh
#	-scp -p $^ zender@biogenic.ess.uci.edu:.ssh
#	-scp -p $^ zender@blackforest.ucar.edu:.ssh
#	-scp -p $^ zender@bluefire.ucar.edu:.ssh
#	-scp -p $^ zender@bluesky.ucar.edu:.ssh
#	-scp -p $^ zender@chinookfe.ucar.edu:.ssh
#	-scp -p $^ zender@edison.nersc.gov:.ssh
#	-scp -p $^ zender@elnino.ess.uci.edu:.ssh
#	-scp -p $^ zender@firn.ess.uci.edu:.ssh
#	-scp -p $^ zender@frazil.ess.uci.edu:.ssh
#	-scp -p $^ zender@givre.ess.uci.edu:.ssh
#	-scp -p $^ zender@glace.ess.uci.edu:.ssh
#	-scp -p $^ zender@goldhill.cgd.ucar.edu:.ssh
#	-scp -p $^ zender@haze.ess.uci.edu:.ssh
#	-scp -p $^ zender@hpc.oit.uci.edu:.ssh
#	-scp -p $^ zender@krein.math.uci.edu:.ssh
#	-scp -p $^ zender@lanina.ess.uci.edu:.ssh
#	-scp -p $^ zender@longs.ucar.edu:.ssh
#	-scp -p $^ zender@mirage0.ucar.edu:.ssh
#	-scp -p $^ zender@nco.sourceforge.net:.ssh
#	-scp -p $^ zender@neige.ess.uci.edu:.ssh
#	-scp -p $^ zender@pbs.ess.uci.edu:.ssh
#	-scp -p $^ zender@roulee.ess.uci.edu:.ssh
#	-scp -p $^ zender@sand.ess.uci.edu:.ssh
#	-scp -p $^ zender@seasalt.ess.uci.edu:.ssh
#	-scp -p $^ zender@skyglow.ess.uci.edu:.ssh
#	-scp -p $^ zender@snow.engin.umich.edu:.ssh
#	-scp -p $^ zender@soot.ess.uci.edu:.ssh
#	-scp -p $^ zender@tempest.ucar.edu:.ssh
#	-scp -p $^ zender@tephra.ess.uci.edu:.ssh
#	-scp -p $^ zender@virga.ess.uci.edu:.ssh
ssh_cln :
	rm -f ${MY_SSH_DIR}/$@ 


