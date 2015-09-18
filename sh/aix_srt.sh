#!/bin/sh 
# $Id$

# Purpose: Executed once per X session

xmodmap -e "pointer = 3 2 1" &
if [ -f ~/.xmodmaprc ] ; then 
    xmodmap ~/.xmodmaprc & 
fi      # endif
if [ -f $DATA/pix/ovum_1280x1024.xbm ] ; then 
    echo
    xsetroot -bitmap $DATA/pix/ovum_1280x1024.xbm &
fi      # endif

if [ -f $HOME/.Xdefaults ] ; then 
    xrdb -merge $HOME/.Xdefaults	
fi	# endif Xdefaults
if [ -f $HOME/idl/ibp_Xdefaults.txt ] ; then 
    xrdb -merge $HOME/idl/ibp_Xdefaults.txt
fi	# endif Xdefaults

(xhost \
	acd.ucar.edu \
	antero.ucar.edu \
	antero2.ucar.edu \
	aztec.ucar.edu \
	aztec2.ucar.edu \
	bearmtn.cgd.ucar.edu \
	chipeta.ucar.edu \
	chipeta2.ucar.edu \
	dsl.ucar.edu \
	garcia.acd.ucar.edu \
	goldhill.cgd.ucar.edu \
	greenmtn.cgd.ucar.edu \
	dust.acd.ucar.edu \
	flagstaf.cgd.ucar.edu \
	neit.cgd.ucar.edu \
	meeker.ucar.edu \
	odin.cgd.ucar.edu \
	ouray.ucar.edu \
	ouray2.ucar.edu \
	paiute.ucar.edu \
	paiute2.ucar.edu \
	sanitas.cgd.ucar.edu \
	st-elmo.ucar.edu \
	st-elmo2.ucar.edu \
	dataproc.ucar.edu \
	dataproc2.ucar.edu \
	dataproc.cgd.ucar.edu \
	ute.ucar.edu \
	ute2.ucar.edu \
	z.ppp.ucar.edu \
	)

