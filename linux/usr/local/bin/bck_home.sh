#!/bin/bash

# $Header: /home/zender/cvs/linux/usr/local/bin/bck_home.sh,v 1.3 2004-05-10 15:26:56 zender Exp $

# Purpose: Cron script for root to backup user home directories

# NB: Place in priveleged directory to satisfy cron security measures

# Usage: 

# bck_home.sh <usr_nm>
# sudo cp ~/linux/usr/local/bin/bck_home.sh /usr/local/bin
# sudo scp ~/linux/usr/local/bin/bck_home.sh dust.ess.uci.edu:/usr/local/bin

if [ $# -ne 1 ]; then
    echo "usage: $0 user"
    exit 1
fi

yyyymmdd=`date +"%Y%m%d"`
yyyymm=`date +"%Y%m"`
hst_nm=`hostname`

usr_nm=$1
tar_fl=/tmp/bck_${hst_nm}_home_${usr_nm}.tar.gz
tar cvfz ${tar_fl} /home/${usr_nm}
#cp ${tar_fl} /dust/bck
cp ${tar_fl} /biogenic/bck
rm -rf ${tar_fl}
