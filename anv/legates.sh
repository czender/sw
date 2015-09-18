# Create Legates files
ncrename -v TS,TS1 /data2/zender/legates/legates_2080_01.nc
ncrename -v TS,TS1 /data2/zender/legates/legates_2080_01.nc

# Difference the surface files
ncdiff -O -D 3 -v TS1,PRECT,gw /data2/zender/spcp_85/spcp_85_8589_01.nc /data2/zender/legates/legates_2080_01.nc /data2/zender/spcp_85/spcp_85_8589_legates_2080_01.nc
ncdiff -O -D 3 -v TS1,PRECT,gw /data2/zender/spcp_85/spcp_85_8589_07.nc /data2/zender/legates/legates_2080_07.nc /data2/zender/spcp_85/spcp_85_8589_legates_2080_07.nc

ncdiff -O -D 3 -v TS1,PRECT,gw /data2/zender/amip5/amip5_8589_01.nc /data2/zender/legates/legates_2080_01.nc /data2/zender/amip5/amip5_8589_legates_2080_01.nc
ncdiff -O -D 3 -v TS1,PRECT,gw /data2/zender/amip5/amip5_8589_07.nc /data2/zender/legates/legates_2080_07.nc /data2/zender/amip5/amip5_8589_legates_2080_07.nc
