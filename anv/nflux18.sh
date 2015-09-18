#! /bin/tcsh -f
foreach fl (`ls nflux18_sho??.nc`)
echo "renaming $fl"
ncrename -v oro,ORO $fl
end

-v gw,LWCF,SWCF,TS1,FLNT,FSNT,CLDLOW,CLDMED,CLDHGH,CLDTOT,TOTCWP

# Create files of the interannual anomalies
ncrcat -D 3 -O -v gw,LWCF,SWCF,TS1 /data2/zender/spcp_85/spcp_85_8589_0112.nc /data2/zender/spcp_85/spcp_85_8589_0112.nc /data2/zender/spcp_85/spcp_85_8589_0112.nc /data2/zender/spcp_85/spcp_85_8589_0112.nc /data2/zender/spcp_85/spcp_85_8589_0112.nc /data2/zender/spcp_85/spcp_85_8589_0160_ac_foo.nc  

ncdiff -D 3 -O -v gw,LWCF,SWCF,TS1 /data2/zender/spcp_85/spcp_85_8589_0160.nc /data2/zender/spcp_85/spcp_85_8589_0160_ac_foo.nc /data2/zender/spcp_85/spcp_85_ia_anom_8589_0160.nc

ncwa -A -D 3 -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/spcp_85/spcp_85_ia_anom_8589_0160.nc

ncwa -O -D 3 -a lat -w gw -d lat,-10.,10. /data2/zender/spcp_85/spcp_85_ia_anom_8589_0160.nc /data2/zender/spcp_85/spcp_85_ia_anom_yavg_10S10N_8589_0160.nc

