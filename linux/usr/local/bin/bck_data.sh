#!/bin/bash

# $Header: /home/zender/cvs/linux/usr/local/bin/bck_data.sh,v 1.5 2004-05-11 18:38:46 zender Exp $

# Purpose: Cron script for root to backup important data

# Usage: 
# bck_home.sh <user>
# sudo cp ~/linux/usr/local/bin/bck_data.sh /usr/local/bin
# sudo scp ~/linux/usr/local/bin/bck_data.sh dust.ess.uci.edu:/usr/local/bin

if [ $# -ne 1 ]; then
    echo "usage: $0 suffix"
    echo "where suffix describes backup frequency"
    echo "e.g.,: $0 dly"
    echo "e.g.,: $0 wk"
    echo "e.g.,: $0 mth"
    echo "e.g.,: $0 yr"
    exit 1
fi

yyyymmdd=`date +"%Y%m%d"`
yyyymm=`date +"%Y%m"`
hst_nm=`hostname`

sfx=$1
tar_fl=/tmp/bck_${hst_nm}_data_${sfx}.tar.gz
tar cvfz ${tar_fl} /data/zender/dst/data /data/zender/fgr /data/zender/ppr
#cp ${tar_fl} /dust/bck
cp ${tar_fl} /biogenic/bck
rm -rf ${tar_fl}
