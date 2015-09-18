#!/bin/bash

# Purpose: quick stats of some variable in netCDF file
#
# Examples:
#   qstats2.sh map.nc hgt_sfc
#   qstats2.sh dstmch14_clm_0112.nc DSTODXC

#---------------------------------------------------------------
usage="usage: $0 file var"
if [ $# -ne 2 ]; then
    echo $usage
    exit 1
fi
file=$1
var=$2
tmpfile=/usr/tmp/tmp.nc
#---------------------------------------------------------------

echo "avg"
ncwa -O -v $var -y avg $file $tmpfile &> /dev/null
ncks -H -v $var              $tmpfile

echo "max"
ncwa -O -v $var -y max $file $tmpfile &> /dev/null
ncks -H -v $var              $tmpfile

echo "min"
ncwa -O -v $var -y min $file $tmpfile &> /dev/null
ncks -H -v $var              $tmpfile

rm -f $tmpfile
