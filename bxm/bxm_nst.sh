#!/bin/sh -x

# $Id$

# Purpose: Install on-line box model on website
# Automate permissions and directories, which are a nightmare, with this script
# bxm_nst.sh is called by Makefile: sudo make bxm

# Usage:
# sudo make bxm
# scp ~/bxm/bxm_nst.sh dust.ess.uci.edu:bxm

MDL_NM=bxm

MY_BXM_DIR=${HOME}/${MDL_NM}
MY_FTP_DIR=/var/ftp/pub/zender
MY_HTTP_DIR=/var/www/html/dead
MY_HTTP_MCH=dust.ess.uci.edu
MY_TMP_DIR=/tmp/${MDL_NM}

if [ 1 = 1 ]; then
# Install on local machine
    mkdir -p ${MY_TMP_DIR}
    /bin/cp -p bxm_run.sh ${MY_TMP_DIR}
    /bin/cp -p bxm.ncl ${MY_HTTP_DIR}
    /bin/cp -p bxm_sz.ncl ${MY_HTTP_DIR}
    /bin/cp -p bxm_cnf.html ${MY_HTTP_DIR}
    /bin/cp -p bxm_dwn.html ${MY_HTTP_DIR}
    /bin/cp -p bxm_home.html ${MY_HTTP_DIR}
    /bin/cp -p bxm_mch.pl ${MY_TMP_DIR}
    /bin/cp -p bxm_vzn.html ${MY_HTTP_DIR}
    touch ${MY_TMP_DIR}/bxm_txt.html
    ln -s -f ${MY_TMP_DIR}/bxm_txt.html ${MY_HTTP_DIR}/bxm_txt.html
    touch ${MY_TMP_DIR}/bxm.jpg
    touch ${MY_HTTP_DIR}/bxm.jpg
    chmod 744 ${MY_HTTP_DIR}/bxm.jpg
    chown -R apache ${MY_HTTP_DIR}/bxm.jpg
    chgrp -R apache ${MY_HTTP_DIR}/bxm.jpg
#    ln -s -f ${MY_TMP_DIR}/bxm.jpg ${MY_HTTP_DIR}/bxm.jpg
    /bin/cp -p bxm_run.py /var/www/cgi-bin
    chown -R apache ${MY_TMP_DIR}
    chgrp -R apache ${MY_TMP_DIR}
    chmod 755 ${MY_TMP_DIR}
    chown apache /var/www/cgi-bin/bxm_run.py ${MY_HTTP_DIR}/dead
    chgrp apache /var/www/cgi-bin/bxm_run.py ${MY_HTTP_DIR}/dead
    mkdir -p /var/ftp/dead
    chown apache /var/ftp/dead
    chgrp apache /var/ftp/dead
else
# Install on ${MY_HTTP_MCH}
# OpenSSH requires logging in as root each time ssh or scp is used! PITA!
    ssh ${MY_HTTP_MCH} mkdir -p ${MY_TMP_DIR}
    scp -p bxm_run.sh ${MY_HTTP_MCH}:${MY_TMP_DIR}
    scp -p bxm.ncl ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    scp -p bxm_sz.ncl ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    scp -p bxm_cnf.html ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    scp -p bxm_dwn.html ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    scp -p bxm_home.html ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    scp -p bxm_mch.pl ${MY_HTTP_MCH}:${MY_TMP_DIR}
    scp -p bxm_vzn.html ${MY_HTTP_MCH}:${MY_HTTP_DIR}
    ssh ${MY_HTTP_MCH} touch ${MY_TMP_DIR}/bxm_txt.html
    ssh ${MY_HTTP_MCH} touch ${MY_TMP_DIR}/bxm.jpg
    ssh ${MY_HTTP_MCH} ln -s -f ${MY_TMP_DIR}/bxm_txt.html ${MY_HTTP_DIR}/bxm_txt.html
#    ssh ${MY_HTTP_MCH} ln -s -f ${MY_TMP_DIR}/bxm.jpg ${MY_HTTP_DIR}/bxm.jpg
    scp -p bxm_run.py ${MY_HTTP_MCH}:/var/www/cgi-bin
    ssh ${MY_HTTP_MCH} 'chown -R apache ${MY_TMP_DIR};chgrp -R apache ${MY_TMP_DIR}'
    ssh ${MY_HTTP_MCH} 'chown apache /var/www/cgi-bin/bxm_run.py ${MY_HTTP_DIR}/dead;chgrp apache /var/www/cgi-bin/bxm_run.py ${MY_HTTP_DIR}/dead'
fi # endif
