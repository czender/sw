#!/bin/csh -v

# Purpose: Create seasonal cycle files of regionally averaged anomalies

# Land
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,10.,30. -d lon,90.,110. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_India_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,10.,25. -d lon,250.,280. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_America_Central_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-60.,-25. -d lon,280.,310. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_America_South_South_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-90.,-65. -d lon,0.,360. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Antarctica_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-10.,10. -d lon,90.,150. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indonesia_Land_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,10.,30. -d lon,90.,120. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indochina_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-10.,5. -d lon,10.,30. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Congo_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,55.,70. -d lon,5.,60. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Europe_Northern_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,30.,50. -d lon,270.,290. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_US_Eastern_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,30.,50. -d lon,250.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_US_Central_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,30.,50. -d lon,230.,250. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_US_Western_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-10.,0. -d lon,290.,310. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Amazon_Basin_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,60.,90. -d lon,300.,340. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Greenland_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,50.,70. -d lon,60.,90. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Siberia_Western_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,50.,70. -d lon,90.,140. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Siberia_Eastern_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,30.,40. -d lon,80.,100. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Tibetan_Plateau_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-40.,-10. -d lon,110.,160. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Australia_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,-35.,-10. -d lon,10.,40. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Africa_South_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,50.,70. -d lon,190.,260. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Alaska_NW_Canada_8589_0112.nc

# Ocean
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,140.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Equatorial_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,140.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Tropical_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-60.,-30. -d lon,150.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_South_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,30.,50. -d lon,150.,210. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_North_8589_0112.nc

ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,140.,180. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Tropical_Western_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,180.,230. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Tropical_Central_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,20. -d lon,230.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Tropical_Eastern_8589_0112.nc

ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,140.,180. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Equatorial_Western_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,180.,230. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Equatorial_Central_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,230.,270. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Equatorial_Eastern_8589_0112.nc

ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,30.,60. -d lon,300.,330. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Atlantic_North_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-20.,10. -d lon,50.,100. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Tropical_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,90.,150. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indonesia_Ocean_8589_0112.nc

ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,140.,170. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_Western_Warm_Pool_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,0.,15. -d lon,140.,200. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_ITCZ_West_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-10.,10. -d lon,50.,100. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Equatorial_8589_0112.nc
ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-15.,5. -d lon,60.,80. -v LWCF,SWCF /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Central_8589_0112.nc

# Both (land and ocean)

ncwa -O -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,10.,25. -d lon,50.,70. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Arabian_Sea_8589_0112.nc

ncks -O -d lat,15.,35. -d lon,340.,60. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/foo.nc
ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,15.,35. -d lon,340.,60. /data2/zender/erbe_b/foo.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Sahara_Arabia_8589_0112.nc

ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,36.6048 -d lon,262.52. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_US_SGP_CART_8589_0112.nc

