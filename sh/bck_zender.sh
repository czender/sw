#!/bin/bash

# $Id$

# Purpose: Cron script for zender to backup important data

# NB: Place in personal directory

# Usage: 
# bck_zender.sh dly
# ~/bin/sh/bck_zender.sh dly

if [ $# -ne 1 ]; then
    echo "usage: $0 suffix"
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
tar_fl="/tmp/bck_${hst_nm}_zender_${sfx}.tar.gz"
# User-priveleged tar commands may not start with '/'
# Do not use -v or --verbose since logs are ~1 MB
xcl='--exclude="data/*"' # Exclude files in /home/zender/data directory
# Use --ignore-failed-read or else command exits on machines without files in exclude pattern
# --wildcards and --wildcards-match-slash are both on by default
# Make them explicit since we use them in exclude list
cd /;tar --create --file ${tar_fl} --gzip --ignore-failed-read --preserve --wildcards --wildcards-match-slash ${xcl} home/zender
scp -p -C ${tar_fl} goldhill.cgd.ucar.edu:/home/zender/data/bck
rm -rf ${tar_fl}
