# Create ECMWF files
ncrename -v ZA,Z2TEST /data2/zender/ecmwf/ecmwf_pres_9095_01.nc
ncrename -v ZA,Z2TEST /data2/zender/ecmwf/ecmwf_pres_8994_07.nc

# Meridionally average the pressure files
ncwa -O -D 3 -a lat -w gw -v U,T,OMEGA,Q,MPSI,Z2TEST -d lat,0.,20. /data2/zender/ecmwf/ecmwf_pres_9095_01.nc /data2/zender/ecmwf/ecmwf_pres_yavg_00N20N_9095_01.nc
ncwa -O -D 3 -a lat -w gw -v U,T,OMEGA,Q,MPSI,Z2TEST -d lat,0.,20. /data2/zender/ecmwf/ecmwf_pres_8994_07.nc /data2/zender/ecmwf/ecmwf_pres_yavg_00N20N_8994_07.nc

# Zonal average the pressure files
ncwa -O -D 3 -a lon -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/ecmwf/ecmwf_pres_9095_01.nc /data2/zender/ecmwf/ecmwf_pres_xavg_9095_01.nc
ncwa -O -D 3 -a lon -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/ecmwf/ecmwf_pres_8994_07.nc /data2/zender/ecmwf/ecmwf_pres_xavg_8994_07.nc

# Difference the pressure files
ncdiff -O -D 3 -v U,Q,Z2TEST,CHI,PSI /data2/zender/amip5/amip5_pres_8589_01.nc /data2/zender/ecmwf/ecmwf_pres_9095_01.nc /data2/zender/amip5/amip5_8589_ecmwf_9095_pres_01.nc
ncdiff -O -D 3 -v U,Q,Z2TEST,CHI,PSI /data2/zender/amip5/amip5_pres_8589_07.nc /data2/zender/ecmwf/ecmwf_pres_8994_07.nc /data2/zender/amip5/amip5_8589_ecmwf_8994_pres_07.nc
ncdiff -O -D 3 -v U,Q,Z2TEST,CHI,PSI /data2/zender/spcp_85/spcp_85_pres_8589_01.nc /data2/zender/ecmwf/ecmwf_pres_9095_01.nc /data2/zender/spcp_85/spcp_85_8589_ecmwf_9095_pres_01.nc
ncdiff -O -D 3 -v U,Q,Z2TEST,CHI,PSI /data2/zender/spcp_85/spcp_85_pres_8589_07.nc /data2/zender/ecmwf/ecmwf_pres_8994_07.nc /data2/zender/spcp_85/spcp_85_8589_ecmwf_8994_pres_07.nc

# Difference the zonal average pressure files
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,PSI /data2/zender/amip5/amip5_pres_xavg_8589_01.nc /data2/zender/ecmwf/ecmwf_pres_xavg_9095_01.nc /data2/zender/amip5/amip5_8589_ecmwf_9095_pres_xavg_01.nc
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,PSI /data2/zender/amip5/amip5_pres_xavg_8589_07.nc /data2/zender/ecmwf/ecmwf_pres_xavg_8994_07.nc /data2/zender/amip5/amip5_8589_ecmwf_8994_pres_xavg_07.nc

ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,PSI /data2/zender/spcp_85/spcp_85_pres_xavg_8589_01.nc /data2/zender/ecmwf/ecmwf_pres_xavg_9095_01.nc /data2/zender/spcp_85/spcp_85_8589_ecmwf_9095_pres_xavg_01.nc
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,PSI /data2/zender/spcp_85/spcp_85_pres_xavg_8589_07.nc /data2/zender/ecmwf/ecmwf_pres_xavg_8994_07.nc /data2/zender/spcp_85/spcp_85_8589_ecmwf_8994_pres_xavg_07.nc

