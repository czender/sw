#!/bin/csh -f -v

# Purpose: Retrieve and process LSM data
# Currently only surface fields are retained

# Usage: /home/zender/lsm/lsm.sh
# NB: u3527a is the 15 year CCM 3.6 AMIP control run

setenv caseid lsm
setenv src_dir /ROSINSKI/csm/u3527a/lnd/hist
#setenv wrk_dir /fs/cgd/data0/zender/${caseid}
#setenv wrk_dir /data/zender/${caseid}
setenv wrk_dir /usr/tmp/zender/${caseid}
mkdir $wrk_dir
cd $wrk_dir

foreach yr (78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93)
foreach mth (01 02 03 04 05 06 07 08 09 10 11 12)
echo "Retrieving ${src_dir}/lsmh_19${yr}-${mth}.nc"
msrcp mss:${src_dir}/lsmh_19${yr}-${mth}.nc $wrk_dir
end
end

foreach mth (01 02 03 04 05 06 07 08 09 10 11 12)
ncra -O -D 3 -d soil_layer_index,0 lsmh_19??-${mth}.nc lsm_7993_${mth}.nc
# LSM uses coordinate names "latitude" and "longitude" instead of "lat" and "lon"
ncrename -v latitude,lat -v longitude,lon lsm_7993_${mth}.nc
# Get rid of soil level dimension
ncwa -O -D 3 -a soil_layer_index lsm_7993_${mth}.nc lsm_7993_${mth}.nc
end

ncrcat -O -D 3 lsm_7993_??.nc lsm_7993_0112.nc
ncra -O -D 3 lsm_7993_0112.nc lsm_7993.nc

rcp $wrk_dir/${caseid}_*.nc odin.cgd.ucar.edu:/data/zender/${caseid}
