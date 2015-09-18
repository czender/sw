#!/bin/sh

# $Id$

# Purpose: Process data from LSM dust module

# Usage: 
# $HOME/lsm/dstlsm.sh

for caseid in dstlsm10 ; do

wrk_dir="/data/zender/${caseid}"
yr='8589'
cd ${wrk_dir}

time_avg="mth"
msk_ocn="-m ORO -o lt -M 0.5 "
msk_lnd="-m ORO -o gt -M 0.5 "

case "$caseid" in
    dstlsm10 ) fld="DSTODXC,ORO,DSTQ,date,time";time_avg="dly"
    ;;
    * )	echo "rgn.sh: ERROR Unable to match caseid"
	exit
    ;;	
esac # endcase $caseid

# Do not use ocean mask for land points!!!

if [ "${time_avg}" = "dly" ]; then
    date_sng="19980101_19981231"
elif [ "${time_avg}" = "mth" ]; then
    date_sng="${yr}"
    date_nsm_sng="${yr}_0112"
else 
    printf "ERROR: \$time_avg is unknown\n"
fi # endif time_avg

for mth in 01 02 03 04 05 06 07 08 09 10 11 12 ; do 
#    mv ${wrk_dir}/lsmh_${mth}.nc ${wrk_dir}/${caseid}_${mth}.nc
    ncwa -O -a lon -C -v LATIXY ${wrk_dir}/${caseid}_${mth}.nc foo.nc
    ncwa -A -a lat -C -v LONGXY ${wrk_dir}/${caseid}_${mth}.nc foo.nc
    ncwa -A -a lat,lon -C -v ZSOI ${wrk_dir}/${caseid}_${mth}.nc foo.nc
    ncrename -v LONGXY,lon -v LATIXY,lat -v ZSOI,lev foo.nc
    ncks -A -C -v lat,lon,lev foo.nc ${wrk_dir}/${caseid}_${mth}.nc 
#    ncks -A -C -v lsmlev ${HOME}/nco/data/in.nc ${wrk_dir}/${caseid}_${mth}.nc 
#    ncrename -v lsmlev,lev ${wrk_dir}/${caseid}_${mth}.nc 
done # end loop over foo

done # end loop over caseid
