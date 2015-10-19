# $Id$

# Purpose: Makefile for web box model files module bxm

# Usage:
# cd ~/sw/bxm;sudo make bxm;cd - # Install box model scripts on web server
# cd ~/sw/bxm;sudo make bxm_dead;cd - # Install current private executable on web server

# System files are installed here
MDL_NM := bxm

MY_BXM_DIR := ${HOME}/sw/${MDL_NM}
MY_DOC_DIR := .
MY_FTP_DIR := /var/ftp/pub/zender
MY_HTTP_DIR := Sites/dead
MY_HTTP_MCH := dust.ess.uci.edu
MY_TMP_DIR := /tmp/${MDL_NM}
ifndef VERSION
 VERSION := $(shell cat ${MY_DOC_DIR}/VERSION)
endif

# Files to install on website
DST_FL := bxm.ncl bxm_cnf.html bxm_dwn.html bxm_home.html bxm_mch.pl bxm_run.py bxm_run.sh bxm_sz.ncl bxm_vzn.html # Distribution files
DST_MFL := ChangeLog Makefile README TODO VERSION # Distribution metafiles
#DST_MFL += FAQ index.shtml INSTALL MANIFEST

# Install box model scripts on web server
bxm : ${MY_BXM_DIR}/bxm_home.html
	-./bxm_nst.sh

# Install current private executable on web server
bxm_dead : ${MY_BIN_DIR}/dead
	-/bin/cp -f $^ ${MY_HTTP_DIR}
	-chmod 744 ${MY_HTTP_DIR}/dead
	-chown apache ${MY_HTTP_DIR}/dead
	-chgrp apache ${MY_HTTP_DIR}/dead

# Make distribution
dst : ${DST_FL} ${DST_MFL}
	tar cvzf ${DATA}/tmp/${MDL_NM}-${VERSION}.tar.gz ${DST_FL} ${DST_MFL}
	scp -p ${DST_FL} ${DATA}/tmp/${MDL_NM}-${VERSION}.tar.gz ${MY_HTTP_MCH}:${MY_HTTP_DIR}
	-ssh ${MY_HTTP_MCH} "cd ${MY_HTTP_DIR};ln -s -f ${MDL_NM}-${VERSION}.tar.gz ${MDL_NM}.tar.gz"

tags : 
	etags *.html *.ncl *.pl *.py *.sh README Makefile
