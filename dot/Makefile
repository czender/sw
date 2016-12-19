# $Id$

# Purpose: Makefile for system files module dot

# System files are installed here
MY_DOT_DIR := ${HOME}
MY_SSH_DIR := ${HOME}/.ssh

# Variables having to do with header files created by module
MDL_DOT_TRG := ChangeLog cshrc mailrc README fvwm2rc rhosts login xinitrc.old signature crontab emacs mwmrc twmrc fvwmrc hluresfile idl idle latex2html-init lightning lightningmask logout mailcap whale xmodmaprc profile.old fvwm95rc plan gnus bashrc bash_profile bash_logout dbxrc gdbinit pgdbgrc cvDebugger Xdefaults cvsignore inputrc forward netrc mpd.conf dupload.conf rpmmacros pbuilderrc wgetrc dput.cf daprc gitconfig soft
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
	-scp -p $^ dust.ess.uci.edu:.ssh
	-scp -p $^ frazil.ess.uci.edu:.ssh
	-scp -p $^ skyglow.ess.uci.edu:.ssh
	-scp -p $^ glace.ess.uci.edu:.ssh
	-scp -p $^ grele.ess.uci.edu:.ssh
	-scp -p $^ gplogin2.ps.uci.edu:.ssh
	-scp -p $^ home.ps.uci.edu:.ssh
	-scp -p $^ givre.ess.uci.edu:.ssh
	-scp -p $^ snow.engin.umich.edu:.ssh
	-scp -p $^ pbs.ess.uci.edu:.ssh
	-scp -p $^ neige.ess.uci.edu:.ssh
	-scp -p $^ roulee.ess.uci.edu:.ssh
#	-scp -p $^ hpc.oit.uci.edu:.ssh
#	-scp -p $^ tephra.ess.uci.edu:.ssh
#	-scp -p $^ virga.ess.uci.edu:.ssh
#	-scp -p $^ seasalt.ess.uci.edu:.ssh
#	-scp -p $^ lanina.ess.uci.edu:.ssh
#	-scp -p $^ sand.ess.uci.edu:.ssh
#	-scp -p $^ biogenic.ess.uci.edu:.ssh
#	-scp -p $^ soot.ess.uci.edu:.ssh
#	-scp -p $^ goldhill.cgd.ucar.edu:.ssh
#	-scp -p $^ haze.ess.uci.edu:.ssh
#	-scp -p $^ elnino.ess.uci.edu:.ssh
#	-scp -p $^ krein.math.uci.edu:.ssh
#	-scp -p $^ nco.sourceforge.net:.ssh
#	-scp -p $^ longs.ucar.edu:.ssh
#	-scp -p $^ tempest.ucar.edu:.ssh
#	-scp -p $^ chinookfe.ucar.edu:.ssh
#	-scp -p $^ blackforest.ucar.edu:.ssh
#	-scp -p $^ bluesky.ucar.edu:.ssh
#	-scp -p $^ bluefire.ucar.edu:.ssh
#	-scp -p $^ babyblue.ucar.edu:.ssh
#	-scp -p $^ mirage0.ucar.edu:.ssh
ssh_cln :
	rm -f ${MY_SSH_DIR}/$@ 


