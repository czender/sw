# Generate the precipitation timeseries
ncecat -O -D 3 -d lon,120.,130. -d lat,8.,22. -v PRECT,PRECC,ORO omega.??.nc gcm_0114.nc
ncrename -O -D 3 -d record,day gcm_0114.nc gcm_0114.nc
ncecat -O -D 3 gcm_0114.nc gcm_0114.nc
ncwa -O -D 3 -m ORO -M .1 -o le -a record gcm_0114.nc gcm_ocn_0114.nc

ncwa -O -D 3 -o gt -m PRECC -M 1.736e-7 -a lon,lat gcm_ocn_0114.nc gcm_ocn_xyavg_0114.nc
ncwa -O -D 3 -a lon,lat gcm_ocn_0114.nc gcm_ocn_xyavg_0114.nc
ncwa -O -D 3 -a day,time gcm_ocn_xyavg_0114.nc gcm_ocn_txyavg_0114.nc
ncks -H -C -v PRECT gcm_ocn_txyavg_0114.nc

ncwa -O -D 3 -o gt -m PRECC -M 1.736e-7 -a lon,lat,time,record,day gcm_0114.nc gcm_txyavg_0114.nc
ncwa -O -D 3 -a lon,lat,time,record,day gcm_0114.nc gcm_txyavg_0114.nc
ncks -H -C -v PRECT gcm_txyavg_0114.nc

# Get the average CMFMC for anv_cem Mc plot
ncea -O -D 3 -d lon,120.,130. -d lat,8.,22. -v CMFMC,ORO omega.??.nc omega_0114.nc
ncwa -O -D 3 -o le -m ORO -M .1 -a time,lat,lon omega_0114.nc omega_txyavg_0114.nc
