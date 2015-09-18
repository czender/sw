ncwa -O -D 3 -v CLDLOW,CLDMED,CLDHGH,CLDTOT -d lat,-20.,20. -a lat,lon -w gw -m ORO -o lt -M .5 isccp_8388.nc foo.nc
ncwa -O -D 3 -v CLDLOW,CLDMED,CLDHGH,CLDTOT -d lat,-20.,20. -a lat,lon -w gw -m ORO -o gt -M .5 isccp_8388.nc foo.nc
ncwa -O -D 3 -v CLDLOW,CLDMED,CLDHGH,CLDTOT -d lat,-20.,20. -a lat,lon -w gw isccp_8388.nc foo.nc
ncks -C -H -v CLDLOW,CLDMED,CLDHGH,CLDTOT foo.nc

ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/isccp/isccp_8388_0112.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-15.,5. -d lon,60.,80. /data2/zender/isccp/isccp_8388_0112.nc /data2/zender/isccp/isccp_xyavg_reg_Indian_Central_8388_0112.nc
