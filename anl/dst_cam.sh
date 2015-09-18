# $Id$

# Purpose: Script to generate data and figures for CAM/SNICAR work

# Testing:
# /bin/cp -r ${DATA}/sncpd02 ${DATA}/sncpd02_bck
# /bin/rm -r ${DATA}/sncpd02;/bin/cp -r ${DATA}/sncpd02_bck ${DATA}/sncpd02

# Production:
# cd ~/anl;dst_cam.sh
# cd ~/anl;time dst_cam.sh > foo 2>&1 &

set -o xtrace

export caseid='sncpd02'
export yr_end='0020' # End year
export yr_srt='0001' # Start year
export trg_yr_sng='clm' # Target year label (what ensembles will be called)

# Derived environment variables
export drc_in="${DATA}/${caseid}"
export drc_out="${DATA}/${caseid}"
export wrk_dir="${DATA}/${caseid}"

# Go to working directory
cd ${wrk_dir}

# Remove clutter from output directory
/bin/rm ${drc_out}/${caseid}_clm??.nc # Remove climatological monthly averages

# Generate derived fields with ncap
# Following fields are non-linear so averaging is meaningless:
# RELHUM, ALBEDO, DPSWETFRC, WND_SPD
for mth in {01..12}; do
  mm=`printf "%02d" $mth`
  fl_lst_in=''
  for yr in `seq ${yr_srt} ${yr_end}`; do
    yyyy=`printf "%04d" ${yr}`
    fl_in=${caseid}.cam2.h0.${yyyy}-${mm}.nc
    fl_lst_in="${fl_lst_in} ${caseid}.cam2.h0.${yyyy}-${mm}.nc"
    printf "Adding to input file list: ${fl_in}\n" ${fl_in}
  done # end loop over yr
  ncra -O -D 1 -o ${drc_out}/${caseid}_clm${mm}.nc -p ${drc_in} ${fl_lst_in}
done # end loop over mth

# Create ensemble monthly means and their zonal and global means
for mth in {01..12}; do
    mm=`printf "%02d" $mth`
    ncra -O -D 3 -p ${wrk_dir} ${caseid}_000[4-6]${mm}.nc ${caseid}_${trg_yr_sng}${mm}.nc
    ncwa -O -D 3 -a lon -p ${wrk_dir} ${caseid}_${trg_yr_sng}${mm}.nc ${caseid}_${trg_yr_sng}${mm}_x.nc
    ncwa -O -D 3 -a time,lat,lon -w gw -p ${wrk_dir} ${caseid}_${trg_yr_sng}${mm}.nc ${caseid}_${trg_yr_sng}${mm}_xy.nc
done # end loop over mth

# Create and globally average ensemble annual mean files
for yr in {${yr_srt}..${yr_end}}; do
    yyyy=`printf "%04d" $yr`
    ncra -O -D 3 -p ${wrk_dir} ${caseid}_${yyyy}[0-1][0-9].nc ${caseid}_${trg_yr_sng}.nc
    ncwa -O -a time,lat,lon -w gw -p ${wrk_dir} ${caseid}_${yyyy}.nc ${caseid}_${trg_yr_sng}_xy.nc
    ncwa -O -a time,lon -p ${wrk_dir} ${caseid}_${yyyy}.nc ${caseid}_${trg_yr_sng}_x.nc
done # end loop over yr

# Make seasonal means
ncra -O -n 3,2,1 -p ${wrk_dir} ${caseid}_${trg_yr_sng}03.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0305.nc
ncra -O -n 3,2,1 -p ${wrk_dir} ${caseid}_${trg_yr_sng}06.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0608.nc
ncra -O -n 3,2,1 -p ${wrk_dir} ${caseid}_${trg_yr_sng}09.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0911.nc
ncra -O -n 3,2,1,12 -p ${wrk_dir} ${caseid}_${trg_yr_sng}12.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_1202.nc
# Create zonal means of seasonal means
ncwa -O -a lon -p ${wrk_dir} ${caseid}_${trg_yr_sng}_0305.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0305_x.nc
ncwa -O -a lon -p ${wrk_dir} ${caseid}_${trg_yr_sng}_0608.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0608_x.nc
ncwa -O -a lon -p ${wrk_dir} ${caseid}_${trg_yr_sng}_0911.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_0911_x.nc
ncwa -O -a lon -p ${wrk_dir} ${caseid}_${trg_yr_sng}_1202.nc ${wrk_dir}/${caseid}_${trg_yr_sng}_1202_x.nc


if [ '1' = '1' ]; then # dbg

cd ${DATA}
for caseid in 'snclgm05 snclgm06 cssncpi03b cssncpi04b sncpd05 sncpd06 cssnc2050_02b cssnc2050_03b' ; do
  ncks -H -v 
done # end loop over caseid
fi # !dbg
