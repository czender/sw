#!/bin/csh

# Multi-month (ensemble) averages produced by this script are only valid for fluxes
# The albedos are erroneous since they are averaged albedos and not albedos computed
# from averaged fluxes

# Get ERBE files to conform to CCM names
foreach fl (`ls /data/zender/erbe_b/erbe_b_??_??.nc`)
echo "renaming $fl"
ncrename -O -D 3 -v .MLWD,FLNT -v .MLWCSD,FLNTC -v .ALBD,ALBEDO -v .TSOLRDCS,SOLIN -v .SSTK,TS1 $fl
end

# Fill in the missing ERBE months to complete 8589_0160 files
/bin/rm /data/zender/erbe_b/erbe_b_85_01.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_85_01.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_06.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_06.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_07.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_07.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_08.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_08.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_09.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_09.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_10.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_10.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_11.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_11.nc
/bin/rm /data/zender/erbe_b/erbe_b_89_12.nc; ln -s /data/zender/erbe_b/erbe_b_85_02.nc /data/zender/erbe_b/erbe_b_89_12.nc

ncea -O -D 3 erbe_b_86_01.nc erbe_b_87_01.nc erbe_b_88_01.nc erbe_b_89_01.nc erbe_b_8689_01.nc
ncea -O -D 3 erbe_b_85_02.nc erbe_b_86_02.nc erbe_b_87_02.nc erbe_b_88_02.nc erbe_b_89_02.nc erbe_b_8589_02.nc
ncea -O -D 3 erbe_b_85_03.nc erbe_b_86_03.nc erbe_b_87_03.nc erbe_b_88_03.nc erbe_b_89_03.nc erbe_b_8589_03.nc
ncea -O -D 3 erbe_b_85_04.nc erbe_b_86_04.nc erbe_b_87_04.nc erbe_b_88_04.nc erbe_b_89_04.nc erbe_b_8589_04.nc
ncea -O -D 3 erbe_b_85_05.nc erbe_b_86_05.nc erbe_b_87_05.nc erbe_b_88_05.nc erbe_b_89_05.nc erbe_b_8589_05.nc
ncea -O -D 3 erbe_b_85_06.nc erbe_b_86_06.nc erbe_b_87_06.nc erbe_b_88_06.nc erbe_b_8588_06.nc
ncea -O -D 3 erbe_b_85_07.nc erbe_b_86_07.nc erbe_b_87_07.nc erbe_b_88_07.nc erbe_b_8588_07.nc
ncea -O -D 3 erbe_b_85_08.nc erbe_b_86_08.nc erbe_b_87_08.nc erbe_b_88_08.nc erbe_b_8588_08.nc
ncea -O -D 3 erbe_b_85_09.nc erbe_b_86_09.nc erbe_b_87_09.nc erbe_b_88_09.nc erbe_b_8588_09.nc
ncea -O -D 3 erbe_b_85_10.nc erbe_b_86_10.nc erbe_b_87_10.nc erbe_b_88_10.nc erbe_b_8588_10.nc
ncea -O -D 3 erbe_b_85_11.nc erbe_b_86_11.nc erbe_b_87_11.nc erbe_b_88_11.nc erbe_b_8588_11.nc
ncea -O -D 3 erbe_b_85_12.nc erbe_b_86_12.nc erbe_b_87_12.nc erbe_b_88_12.nc erbe_b_8588_12.nc

# Fill in the missing ERBE months to complete 8589 files
ln -s /data/zender/erbe_b/erbe_b_8689_01.nc /data/zender/erbe_b/erbe_b_8589_01.nc
ln -s /data/zender/erbe_b/erbe_b_8588_06.nc /data/zender/erbe_b/erbe_b_8589_06.nc
ln -s /data/zender/erbe_b/erbe_b_8588_07.nc /data/zender/erbe_b/erbe_b_8589_07.nc
ln -s /data/zender/erbe_b/erbe_b_8588_08.nc /data/zender/erbe_b/erbe_b_8589_08.nc
ln -s /data/zender/erbe_b/erbe_b_8588_09.nc /data/zender/erbe_b/erbe_b_8589_09.nc
ln -s /data/zender/erbe_b/erbe_b_8588_10.nc /data/zender/erbe_b/erbe_b_8589_10.nc
ln -s /data/zender/erbe_b/erbe_b_8588_11.nc /data/zender/erbe_b/erbe_b_8589_11.nc
ln -s /data/zender/erbe_b/erbe_b_8588_12.nc /data/zender/erbe_b/erbe_b_8589_12.nc

ncra -O -D 3 erbe_b_8689_01.nc erbe_b_8589_02.nc erbe_b_8589_03.nc erbe_b_8589_04.nc erbe_b_8589_05.nc erbe_b_8588_06.nc erbe_b_8588_07.nc erbe_b_8588_08.nc erbe_b_8588_09.nc erbe_b_8588_10.nc erbe_b_8588_11.nc erbe_b_8588_12.nc erbe_b_8589.nc

ncrcat -O -D 3 erbe_b_8689_01.nc erbe_b_8589_02.nc erbe_b_8589_03.nc erbe_b_8589_04.nc erbe_b_8589_05.nc erbe_b_8588_06.nc erbe_b_8588_07.nc erbe_b_8588_08.nc erbe_b_8588_09.nc erbe_b_8588_10.nc erbe_b_8588_11.nc erbe_b_8588_12.nc erbe_b_8589_0112.nc

ncea -O -D 3 erbe_b_86_01.nc erbe_b_88_01.nc erbe_b_89_01.nc erbe_b_8689x87_01.nc
ncea -O -D 3 erbe_b_85_02.nc erbe_b_86_02.nc erbe_b_88_02.nc erbe_b_89_02.nc erbe_b_8589x87_02.nc
ncea -O -D 3 erbe_b_85_03.nc erbe_b_86_03.nc erbe_b_88_03.nc erbe_b_89_03.nc erbe_b_8589x87_03.nc
ncea -O -D 3 erbe_b_85_04.nc erbe_b_86_04.nc erbe_b_88_04.nc erbe_b_89_04.nc erbe_b_8589x87_04.nc
ncea -O -D 3 erbe_b_85_05.nc erbe_b_86_05.nc erbe_b_88_05.nc erbe_b_89_05.nc erbe_b_8589x87_05.nc
ncea -O -D 3 erbe_b_85_06.nc erbe_b_86_06.nc erbe_b_88_06.nc erbe_b_8588x87_06.nc
ncea -O -D 3 erbe_b_85_07.nc erbe_b_86_07.nc erbe_b_88_07.nc erbe_b_8588x87_07.nc
ncea -O -D 3 erbe_b_85_08.nc erbe_b_86_08.nc erbe_b_88_08.nc erbe_b_8588x87_08.nc
ncea -O -D 3 erbe_b_85_09.nc erbe_b_86_09.nc erbe_b_88_09.nc erbe_b_8588x87_09.nc
ncea -O -D 3 erbe_b_85_10.nc erbe_b_86_10.nc erbe_b_88_10.nc erbe_b_8588x87_10.nc
ncea -O -D 3 erbe_b_85_11.nc erbe_b_86_11.nc erbe_b_88_11.nc erbe_b_8588x87_11.nc
ncea -O -D 3 erbe_b_85_12.nc erbe_b_86_12.nc erbe_b_88_12.nc erbe_b_8588x87_12.nc

ncea -O -D 3 erbe_b_8689x87_01.nc erbe_b_8589x87_02.nc erbe_b_8589x87_03.nc erbe_b_8589x87_04.nc erbe_b_8589x87_05.nc erbe_b_8588x87_06.nc erbe_b_8588x87_07.nc erbe_b_8588x87_08.nc erbe_b_8588x87_09.nc erbe_b_8588x87_10.nc erbe_b_8588x87_11.nc erbe_b_8588x87_12.nc erbe_b_8589x87.nc

ncrcat -O -D 3 erbe_b_8689x87_01.nc erbe_b_8589x87_02.nc erbe_b_8589x87_03.nc erbe_b_8589x87_04.nc erbe_b_8589x87_05.nc erbe_b_8588x87_06.nc erbe_b_8588x87_07.nc erbe_b_8588x87_08.nc erbe_b_8588x87_09.nc erbe_b_8588x87_10.nc erbe_b_8588x87_11.nc erbe_b_8588x87_12.nc erbe_b_8589x87_0112.nc

# Create the gridpoint interannual file
#ncrcat -D 3 -O /data/zender/erbe_b/erbe_b_??_??.nc erbe_b_8589_0160.nc
ncrcat -O -R -D 3 -l ./ -p /ZENDER/proc/erbe_b erbe_b_85_01.nc erbe_b_85_02.nc erbe_b_85_03.nc erbe_b_85_04.nc erbe_b_85_05.nc erbe_b_85_06.nc erbe_b_85_07.nc erbe_b_85_08.nc erbe_b_85_09.nc erbe_b_85_10.nc erbe_b_85_11.nc erbe_b_85_12.nc erbe_b_86_01.nc erbe_b_86_02.nc erbe_b_86_03.nc erbe_b_86_04.nc erbe_b_86_05.nc erbe_b_86_06.nc erbe_b_86_07.nc erbe_b_86_08.nc erbe_b_86_09.nc erbe_b_86_10.nc erbe_b_86_11.nc erbe_b_86_12.nc erbe_b_87_01.nc erbe_b_87_02.nc erbe_b_87_03.nc erbe_b_87_04.nc erbe_b_87_05.nc erbe_b_87_06.nc erbe_b_87_07.nc erbe_b_87_08.nc erbe_b_87_09.nc erbe_b_87_10.nc erbe_b_87_11.nc erbe_b_87_12.nc erbe_b_88_01.nc erbe_b_88_02.nc erbe_b_88_03.nc erbe_b_88_04.nc erbe_b_88_05.nc erbe_b_88_06.nc erbe_b_88_07.nc erbe_b_88_08.nc erbe_b_88_09.nc erbe_b_88_10.nc erbe_b_88_11.nc erbe_b_88_12.nc erbe_b_89_01.nc erbe_b_89_02.nc erbe_b_89_03.nc erbe_b_89_04.nc erbe_b_89_05.nc erbe_b_89_06.nc erbe_b_89_07.nc erbe_b_89_08.nc erbe_b_89_09.nc erbe_b_89_10.nc erbe_b_89_11.nc erbe_b_89_12.nc erbe_b_8589_0160.nc

# Create the gridpoint interannual anomaly file from the monthly data and one of three means:
# The first mean averages the seasonal cycle, so the ensemble average months are equally weighted.
# This is comparable to GCM data, where all 60 months exist.
ncwa -O -D 3 -a time erbe_b_8589.nc erbe_b_8589_foo.nc
# The second mean averages only the valid ERBE data, 52 months of it. 
# This ensures the 52-month anomalies should add up to 0., but should not be used in production.
ncwa -O -D 3 -a time -d time,1,52 erbe_b_8589_0160.nc erbe_b_8589_foo.nc
# Subtract the mean from the monthly to get the anomaly
ncdiff -D 3 -O erbe_b_8589_0160.nc erbe_b_8589_foo.nc erbe_b_anom_8589_0160.nc
/bin/rm -f erbe_b_8589_foo.nc

# Add the ORO field from SEP1
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8589.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8589x87.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8589x87_0112.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_anom_8589x87_0112.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8589_0112.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8589_0160.nc
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_anom_8589_0160.nc

# Difference the hybrid files
ncdiff -O -D 3 -v gw,ORO,LWCF,SWCF,TS1,GCLD,GCLR,FSNT,FLNT /data/zender/amip5/amip5_8589_01.nc /data/zender/erbe_b/erbe_b_8589_01.nc /data/zender/amip5/amip5_8589_erbe_b_8589_01.nc
ncdiff -O -D 3 -v gw,ORO,LWCF,SWCF,TS1,GCLD,GCLR,FSNT,FLNT /data/zender/amip5/amip5_8589_07.nc /data/zender/erbe_b/erbe_b_8589_07.nc /data/zender/amip5/amip5_8589_erbe_b_8589_07.nc
ncdiff -O -D 3 -v gw,ORO,LWCF,SWCF,TS1,GCLD,GCLR,FSNT,FLNT /data/zender/spcp_85/spcp_85_8589_01.nc /data/zender/erbe_b/erbe_b_8589_01.nc /data/zender/spcp_85/spcp_85_8589_erbe_b_8589_01.nc
ncdiff -O -D 3 -v gw,ORO,LWCF,SWCF,TS1,GCLD,GCLR,FSNT,FLNT /data/zender/spcp_85/spcp_85_8589_07.nc /data/zender/erbe_b/erbe_b_8589_07.nc /data/zender/spcp_85/spcp_85_8589_erbe_b_8589_07.nc

# Meridionally average (and avoid leaving ORO in a source ensemble file)
cp /data/zender/erbe_b/erbe_b_8689_01.nc /data/zender/erbe_b/erbe_b_8689_01.nc_foo 
cp /data/zender/erbe_b/erbe_b_8588_07.nc  /data/zender/erbe_b/erbe_b_8588_07.nc_foo
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8689_01.nc_foo
ncwa -A -C -a time -v ORO /data/zender/ccm/SEP1.nc /data/zender/erbe_b/erbe_b_8588_07.nc_foo
ncwa -O -a lat -w gw -m ORO -M 0. -o eq -d lat,0.,20. /data/zender/erbe_b/erbe_b_8689_01.nc_foo /data/zender/erbe_b/erbe_b_ocn_yavg_00N20N_8589_01.nc
ncwa -O -a lat -w gw -m ORO -M 0. -o eq -d lat,0.,20. /data/zender/erbe_b/erbe_b_8588_07.nc_foo /data/zender/erbe_b/erbe_b_ocn_yavg_00N20N_8589_07.nc
/bin/rm /data/zender/erbe_b/*foo*

# Create the interannual files of meridionally averaged anomalies
ncwa -O -D 3 -a lat -w gw -d lat,-10.,10. erbe_b_anom_8589_0160.nc erbe_b_anom_yavg_10S10N_8589_0160.nc

# Create the interannual files of regional averages for DNS paper
# Absolute files
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,140.,270. /data/zender/erbe_b/erbe_b_8589_0160.nc /data/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,310.,360. /data/zender/erbe_b/erbe_b_8589_0160.nc /data/zender/erbe_b/erbe_b_xyavg_reg_Atlantic_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,50.,100. /data/zender/erbe_b/erbe_b_8589_0160.nc /data/zender/erbe_b/erbe_b_xyavg_reg_Indian_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,0.,360. /data/zender/erbe_b/erbe_b_8589_0160.nc /data/zender/erbe_b/erbe_b_xyavg_reg_Ocean_Tropical_8589_0160.nc

# Anomaly files
# Does it make sense to spatially average the global anomaly file? 
# Yes and No: the temporal and spatial averaging operators are only linear as long as there is no missing data
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,140.,270. /data/zender/erbe_b/erbe_b_anom_8589_0160.nc /data/zender/erbe_b/erbe_b_anom_xyavg_reg_Pacific_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,310.,360. /data/zender/erbe_b/erbe_b_anom_8589_0160.nc /data/zender/erbe_b/erbe_b_anom_xyavg_reg_Atlantic_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,50.,100. /data/zender/erbe_b/erbe_b_anom_8589_0160.nc /data/zender/erbe_b/erbe_b_anom_xyavg_reg_Indian_Tropical_8589_0160.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,0.,360. /data/zender/erbe_b/erbe_b_anom_8589_0160.nc /data/zender/erbe_b/erbe_b_anom_xyavg_reg_Ocean_Tropical_8589_0160.nc

# Quality control: true anomalies sum to 0; temporal and spatial averaging is linear (order is unimportant)
# This must work to machine precision for fields without any missing values, see CZP IV p. 16
# FLNT has no missing values, SWCF does
# First create new, consistent anomaly files
ncwa -O -D 3 -a time -v ORO,FLNT,SWCF,gw -p /data/zender/erbe_b erbe_b_8589_0160.nc foo_8589.nc
ncdiff -O -D 3 -v ORO,FLNT,SWCF,gw -p /data/zender/erbe_b erbe_b_8589_0160.nc foo_8589.nc foo_anom_8589_0160.nc

# Temporally average the spatially averaged anomalies 
ncwa -O -D 3 -a lon,lat -w gw -v FLNT,SWCF -d lat,-20.,20. -d lon,140.,270. -p /data/zender/erbe_b foo_anom_8589_0160.nc foo_foo_anom_8589_0160.nc
ncwa -O -D 3 -a time -v FLNT,SWCF -p /data/zender/erbe_b foo_foo_anom_8589_0160.nc foo_foo_foo_anom_8589_0160.nc

# Spatially average the temporally averaged anomalies 
# NB: In this order of operations, the spatial average operates on points that are all zero to machine precision 
# Hence, this order is numerically unstable. To verify this, sum the spatial region without weights.
ncwa -O -D 3 -a time -v FLNT,SWCF,gw -p /data/zender/erbe_b foo_anom_8589_0160.nc foo_foo_anom_8589_0160.nc
ncwa -O -D 3 -a lon,lat -w gw -v FLNT,SWCF -d lat,-20.,20. -d lon,140.,270. -p /data/zender/erbe_b foo_foo_anom_8589_0160.nc foo_foo_foo_anom_8589_0160.nc
ncwa -O -D 3 -n -a lon,lat -v FLNT,SWCF -d lat,-20.,20. -d lon,140.,270. -p /data/zender/erbe_b foo_foo_anom_8589_0160.nc foo_foo_foo_anom_8589_0160.nc

# Spatially and temporally average simultaneously
ncwa -O -D 3 -a lon,lat,time -w gw -v FLNT,SWCF,gw -d lat,-20.,20. -d lon,140.,270. -p /data/zender/erbe_b foo_anom_8589_0160.nc foo_foo_foo_anom_8589_0160.nc

ncks -H -C -v FLNT,SWCF -d lat,0 -d lon,0 -p /data/zender/erbe_b foo_8589.nc
ncks -H -C -v FLNT,SWCF -d lat,0 -d lon,0 -p /data/zender/erbe_b erbe_b_8589_0160.nc
ncks -H -C -v FLNT,SWCF -d lat,0 -d lon,0 -p /data/zender/erbe_b foo_anom_8589_0160.nc
ncks -H -C -v FLNT,SWCF -d lat,0. -d lon,0. -p /data/zender/erbe_b foo_foo_anom_8589_0160.nc
ncks -H -C -v FLNT,SWCF -p /data/zender/erbe_b foo_foo_anom_8589_0160.nc
ncks -H -C -v FLNT,SWCF -p /data/zender/erbe_b foo_foo_foo_anom_8589_0160.nc
#/bin/rm -f /data/zender/erbe_b/*foo*

