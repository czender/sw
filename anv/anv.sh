#! /bin/tcsh -f
# $Id$

# Purpose: Generate useful hyperslabs from monthly CCM data

# Create ensemble averages
ncra -O -D 3 sld012d_??_??.nc sld012d_8589.nc
ncra -O -D 3 sld012d_85_??.nc sld012d_86_??.nc sld012d_87_??.nc sld012d_88_??.nc sld012d_89_??.nc sld012d_8589.nc
ncra -O -D 3 sld012d_8589_??.nc sld012d_8589.nc
mswrite -t 365 sld012d_8589.nc /ZENDER/proc/sld012d/sld012d_8589.nc &
rcp sld012d_8589.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &

ncra -O -D 3 sld012d_85_01.nc sld012d_86_01.nc sld012d_87_01.nc sld012d_88_01.nc sld012d_89_01.nc sld012d_8589_01.nc
ncra -O -D 3 sld012d_85_02.nc sld012d_86_02.nc sld012d_87_02.nc sld012d_88_02.nc sld012d_89_02.nc sld012d_8589_02.nc
ncra -O -D 3 sld012d_85_03.nc sld012d_86_03.nc sld012d_87_03.nc sld012d_88_03.nc sld012d_89_03.nc sld012d_8589_03.nc
ncra -O -D 3 sld012d_85_04.nc sld012d_86_04.nc sld012d_87_04.nc sld012d_88_04.nc sld012d_89_04.nc sld012d_8589_04.nc
ncra -O -D 3 sld012d_85_05.nc sld012d_86_05.nc sld012d_87_05.nc sld012d_88_05.nc sld012d_89_05.nc sld012d_8589_05.nc
ncra -O -D 3 sld012d_85_06.nc sld012d_86_06.nc sld012d_87_06.nc sld012d_88_06.nc sld012d_89_06.nc sld012d_8589_06.nc
ncra -O -D 3 sld012d_85_07.nc sld012d_86_07.nc sld012d_87_07.nc sld012d_88_07.nc sld012d_89_07.nc sld012d_8589_07.nc
ncra -O -D 3 sld012d_85_08.nc sld012d_86_08.nc sld012d_87_08.nc sld012d_88_08.nc sld012d_89_08.nc sld012d_8589_08.nc
ncra -O -D 3 sld012d_85_09.nc sld012d_86_09.nc sld012d_87_09.nc sld012d_88_09.nc sld012d_89_09.nc sld012d_8589_09.nc
ncra -O -D 3 sld012d_85_10.nc sld012d_86_10.nc sld012d_87_10.nc sld012d_88_10.nc sld012d_89_10.nc sld012d_8589_10.nc
ncra -O -D 3 sld012d_85_11.nc sld012d_86_11.nc sld012d_87_11.nc sld012d_88_11.nc sld012d_89_11.nc sld012d_8589_11.nc
ncra -O -D 3 sld012d_85_12.nc sld012d_86_12.nc sld012d_87_12.nc sld012d_88_12.nc sld012d_89_12.nc sld012d_8589_12.nc

rcp sld012d_8589_01.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
rcp sld012d_8589_07.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
foreach fl (`ls sld012d_8589_??.nc`)
echo "mswriting $fl"
mswrite -t 365 $fl /ZENDER/proc/sld012d/$fl &
end

ncra -O -D 3 sld012d_8589_12.nc sld012d_8589_01.nc sld012d_8589_02.nc sld012d_8589_1202.nc
ncra -O -D 3 sld012d_8589_03.nc sld012d_8589_04.nc sld012d_8589_05.nc sld012d_8589_0305.nc
ncra -O -D 3 sld012d_8589_06.nc sld012d_8589_07.nc sld012d_8589_08.nc sld012d_8589_0608.nc
ncra -O -D 3 sld012d_8589_09.nc sld012d_8589_10.nc sld012d_8589_11.nc sld012d_8589_0911.nc
rcp sld012d_8589_1202.nc sld012d_8589_0305.nc sld012d_8589_0608.nc sld012d_8589_0911.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
foreach fl (`ls sld012d_8589_????.nc`)
echo "mswriting $fl"
mswrite -t 365 $fl /ZENDER/proc/sld012d/$fl &
end

ncra -O -D 3 sld012d_85_??.nc sld012d_86_??.nc sld012d_88_??.nc sld012d_89_??.nc sld012d_8589x87.nc
mswrite -t 365 sld012d_8589x87.nc /ZENDER/proc/sld012d/sld012d_8589x87.nc &
rcp sld012d_8589x87.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &

ncra -O -D 3 sld012d_85_01.nc sld012d_86_01.nc sld012d_88_01.nc sld012d_89_01.nc sld012d_8589x87_01.nc
ncra -O -D 3 sld012d_85_02.nc sld012d_86_02.nc sld012d_88_02.nc sld012d_89_02.nc sld012d_8589x87_02.nc
ncra -O -D 3 sld012d_85_03.nc sld012d_86_03.nc sld012d_88_03.nc sld012d_89_03.nc sld012d_8589x87_03.nc
ncra -O -D 3 sld012d_85_04.nc sld012d_86_04.nc sld012d_88_04.nc sld012d_89_04.nc sld012d_8589x87_04.nc
ncra -O -D 3 sld012d_85_05.nc sld012d_86_05.nc sld012d_88_05.nc sld012d_89_05.nc sld012d_8589x87_05.nc
ncra -O -D 3 sld012d_85_06.nc sld012d_86_06.nc sld012d_88_06.nc sld012d_89_06.nc sld012d_8589x87_06.nc
ncra -O -D 3 sld012d_85_07.nc sld012d_86_07.nc sld012d_88_07.nc sld012d_89_07.nc sld012d_8589x87_07.nc
ncra -O -D 3 sld012d_85_08.nc sld012d_86_08.nc sld012d_88_08.nc sld012d_89_08.nc sld012d_8589x87_08.nc
ncra -O -D 3 sld012d_85_09.nc sld012d_86_09.nc sld012d_88_09.nc sld012d_89_09.nc sld012d_8589x87_09.nc
ncra -O -D 3 sld012d_85_10.nc sld012d_86_10.nc sld012d_88_10.nc sld012d_89_10.nc sld012d_8589x87_10.nc
ncra -O -D 3 sld012d_85_11.nc sld012d_86_11.nc sld012d_88_11.nc sld012d_89_11.nc sld012d_8589x87_11.nc
ncra -O -D 3 sld012d_85_12.nc sld012d_86_12.nc sld012d_88_12.nc sld012d_89_12.nc sld012d_8589x87_12.nc

rcp sld012d_8589x87_01.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
rcp sld012d_8589x87_07.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
foreach fl (`ls sld012d_8589x87_??.nc`)
echo "mswriting $fl"
mswrite -t 365 $fl /ZENDER/proc/sld012d/$fl &
end

# Difference the hybrid files
foreach mth (01 02 03 04 05 06 07 08 09 10 11 12)
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA sld012d_8589_$mth.nc sld012d_8589_$mth.nc sld012d_8589_sld012d_8589_$mth.nc
end
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA /data2/zender/sld012d/sld012d_8589_04.nc /data2/zender/sld012d/sld012d_8589_04.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_04.nc
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA /data2/zender/sld012d/sld012d_8589_10.nc /data2/zender/sld012d/sld012d_8589_10.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_10.nc

# Zonal average the hybrid files
foreach mth (01 02 03 04 05 06 07 08 09 10 11 12)
ncwa -O -D 3 -a lon /data2/zender/sld012d/sld012d_8589_$mth.nc /data2/zender/sld012d/sld012d_xavg_8589_$mth.nc
end
ncwa -O -D 3 -a lon /data2/zender/sld012d/sld012d_8589_01.nc /data2/zender/sld012d/sld012d_xavg_8589_01.nc
ncwa -O -D 3 -a lon /data2/zender/sld012d/sld012d_8589_07.nc /data2/zender/sld012d/sld012d_xavg_8589_07.nc

# Difference the zonal average hybrid files
foreach mth (01 02 03 04 05 06 07 08 09 10 11 12)
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA /data2/zender/sld012d/sld012d_xavg_8589_$mth.nc /data2/zender/sld012d/sld012d_xavg_8589_$mth.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_xavg_$mth.nc
end
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA /data2/zender/sld012d/sld012d_xavg_8589_01.nc /data2/zender/sld012d/sld012d_xavg_8589_01.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_xavg_01.nc
ncdiff -O -D 3 -x -v FICE,REI,FSDS,GBIWC,FSDSTOA /data2/zender/sld012d/sld012d_xavg_8589_07.nc /data2/zender/sld012d/sld012d_xavg_8589_07.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_xavg_07.nc

# Meridionally average the hybrid files
ncwa -O -D 3 -a lat -w gw -d lat,0.,30. /data2/zender/sld012d/sld012d_8589_01.nc /data2/zender/sld012d/sld012d_yavg_00N30N_8589_01.nc
ncwa -O -D 3 -a lat -w gw -d lat,0.,30. /data2/zender/sld012d/sld012d_8589_07.nc /data2/zender/sld012d/sld012d_yavg_00N30N_8589_07.nc
ncdiff -O -D 3 -v CMFMC,CLOUD,QRS,QRL,RADD,QICE,QC,Q,HGS,CMFDT,QDIABAT,T,VMAGSFC /data2/zender/sld012d/sld012d_yavg_00N30N_8589_01.nc /data2/zender/sld012d/sld012d_yavg_00N30N_8589_01.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_yavg_00N30N_01.nc
ncdiff -O -D 3 -v CMFMC,CLOUD,QRS,QRL,RADD,QICE,QC,Q,HGS,CMFDT,QDIABAT,T,VMAGSFC /data2/zender/sld012d/sld012d_yavg_00N30N_8589_07.nc /data2/zender/sld012d/sld012d_yavg_00N30N_8589_07.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_yavg_00N30N_07.nc
ncwa -O -D 3 -a lat -w gw -m ORO -M 0. -o eq -d lat,0.,30. /data2/zender/sld012d/sld012d_8589_01.nc /data2/zender/sld012d/sld012d_ocn_yavg_00N30N_8589_01.nc
ncwa -O -D 3 -a lat -w gw -m ORO -M 0. -o eq -d lat,0.,30. /data2/zender/sld012d/sld012d_8589_07.nc /data2/zender/sld012d/sld012d_ocn_yavg_00N30N_8589_07.nc

# Difference the pressure files
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,DIV,PSI /data2/zender/sld012d/sld012d_pres_8589_01.nc /data2/zender/sld012d/sld012d_pres_8589_01.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_pres_01.nc
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST,CHI,DIV,PSI /data2/zender/sld012d/sld012d_pres_8589_07.nc /data2/zender/sld012d/sld012d_pres_8589_07.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_pres_07.nc

# Create and difference yearly pressure files
ncrcat -O -D 3 -v Z2TEST /data2/zender/sld012d/sld012d_pres_??_01.nc /data2/zender/sld012d/sld012d_pres_8589_0101.nc
ncrcat -O -D 3 -v Z2TEST /data2/zender/sld012d/sld012d_pres_??_01.nc /data2/zender/sld012d/sld012d_pres_8589_0101.nc
ncdiff -O -D 3 /data2/zender/sld012d/sld012d_pres_8589_0101.nc /data2/zender/sld012d/sld012d_pres_8589_0101.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_pres_0101.nc

# Meridionally average the pressure files
ncwa -O -D 3 -a lat -w gw -v U,T,OMEGA,Q,MPSI,Z2TEST -d lat,0.,30. /data2/zender/sld012d/sld012d_pres_8589_01.nc /data2/zender/sld012d/sld012d_pres_yavg_00N30N_8589_01.nc
ncwa -O -D 3 -a lat -w gw -v U,T,OMEGA,Q,MPSI,Z2TEST -d lat,0.,30. /data2/zender/sld012d/sld012d_pres_8589_07.nc /data2/zender/sld012d/sld012d_pres_yavg_00N30N_8589_07.nc

# Zonal average the pressure files
ncwa -O -D 3 -a lon -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/sld012d/sld012d_pres_8589_01.nc /data2/zender/sld012d/sld012d_pres_xavg_8589_01.nc
ncwa -O -D 3 -a lon -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/sld012d/sld012d_pres_8589_07.nc /data2/zender/sld012d/sld012d_pres_xavg_8589_07.nc

# Difference the zonal average pressure files
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/sld012d/sld012d_pres_xavg_8589_01.nc /data2/zender/sld012d/sld012d_pres_xavg_8589_01.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_pres_xavg_01.nc
ncdiff -O -D 3 -v U,T,OMEGA,Q,MPSI,Z2TEST /data2/zender/sld012d/sld012d_pres_xavg_8589_07.nc /data2/zender/sld012d/sld012d_pres_xavg_8589_07.nc /data2/zender/sld012d/sld012d_8589_sld012d_8589_pres_xavg_07.nc

# Add the ORO field from SEP1
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_8589.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_8589x87.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_8589x87_0112.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_anom_8589x87_0112.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_8589_0160.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_anom_8589_0160.nc

# Create the SEB files
ncrcat -O -R -D 3 -v gw,LWCF,SWCF,TS1,FSDS,FSNTC,FSNT,ALBEDO,FLNT,FSNS,FSNSC,FLNS,VMAGSFC,SHFLX,LHFLX,E_P,NET,QRS,QRL,RADD,QDIABAT,HGS,HDFF,CMFDT,CMFMC,Q,QC,CLDLOW,CLDMED,CLDHGH,CLDTOT,TOTCWP -l ./ -p /ZENDER/proc/sld012d -n 12,2,1 sld012d_8589_01.nc sld012d_8589_0112.nc
#ncrcat -O -R -D 3 -v gw,LWCF,SWCF,TS1,FLNT,FSNT,CLDLOW,CLDMED,CLDHGH,CLDTOT,TOTCWP -l ./ -p /ZENDER/proc/sld012d -n 12,2,1 sld012d_8589_01.nc sld012d_8589_0112.nc
mswrite -t 365 sld012d_8589_0112.nc /ZENDER/proc/sld012d/sld012d_8589_0112.nc &
rcp sld012d_8589_0112.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_8589_0112.nc

ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,0.,15. -d lon,140.,200. /data2/zender/erbe_bd/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_ITCZ_West_8589_0112.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,-15.,5. -d lon,60.,80. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Central_8589_0112.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 0. -w gw -d lat,10.,25. -d lon,50.,70. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Arabian_Sea_8589_0112.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,15.,35. -d lon,0.,60. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Sahara_Arabia_8589_0112.nc
ncwa -O -D 3 -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,36.6048 -d lon,262.52. /data2/zender/erbe_b/erbe_b_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_US_SGP_CART_8589_0112.nc
ncdiff -O /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_ITCZ_West_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Pacific_ITCZ_West_8589_0112.nc /data2/zender/erbe_b/erbe_b_8589_erbe_b_8589_xyavg_reg_Pacific_ITCZ_West_0112.nc
ncdiff -O /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Central_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Indian_Central_8589_0112.nc /data2/zender/erbe_b/erbe_b_8589_erbe_b_8589_xyavg_reg_Indian_Central_0112.nc
ncdiff -O /data2/zender/erbe_b/erbe_b_xyavg_reg_Arabian_Sea_8589_0112.nc /data2/zender/erbe_b/erbe_b_xyavg_reg_Arabian_Sea_8589_0112.nc /data2/zender/erbe_b/erbe_b_8589_erbe_b_8589_xyavg_reg_Arabian_Sea_0112.nc

# Create the gridpoint seasonal cycle anomaly file
ncwa -O -a time -v FLNT,LWCF,SWCF,TS1,gw /data2/zender/sld012d/sld012d_8589.nc /data2/zender/sld012d/sld012d_8589_foo.nc
ncdiff -D 3 -O -v FLNT,LWCF,SWCF,TS1,gw /data2/zender/sld012d/sld012d_8589_0112.nc /data2/zender/sld012d/sld012d_8589_foo.nc /data2/zender/sld012d/sld012d_anom_8589_0112.nc
ncwa -A -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_anom_8589_0112.nc
/bin/rm -f /data2/zender/sld012d/sld012d_8589_foo.nc

# Create the gridpoint interannual anomaly file
#ncrcat -O -D 3 -v LWCF,SWCF,FLNT,FLNTC,ALBEDO,TS1,FLUS,GCLD,GCLR,gw /usr/tmp/zender/sld012d_??_??.nc sld012d_8589_0160.nc
#ncrcat -O -D 3 -v gw,LWCF,SWCF,TS1,FLNT,FSNT,CLDLOW,CLDMED,CLDHGH,CLDTOT,TOTCWP sld012d_85_??.nc sld012d_86_??.nc sld012d_87_??.nc sld012d_88_??.nc sld012d_89_??.nc sld012d_8589_0160.nc
ncrcat -O -D 3 -v SHFLX,LHFLX,FSNT,PRECT,PRECC,TOTCWP,CLDTOT,CLDHGH,CLDLOW,CLDMED,LWCF,SWCF,FLNT,FLNTC,ALBEDO,TS1,FLUS,GCLD,GCLR,gw /usr/tmp/zender/sld012d_??_??.nc sld012d_8589_0160.nc
ncrcat -O -R -D 3 -v LWCF,SWCF,FLNT,FLNTC,ALBEDO,TS1,FLUS,GCLD,GCLR,gw -l ./ -p /ZENDER/proc/sld012d sld012d_85_01.nc sld012d_85_02.nc sld012d_85_03.nc sld012d_85_04.nc sld012d_85_05.nc sld012d_85_06.nc sld012d_85_07.nc sld012d_85_08.nc sld012d_85_09.nc sld012d_85_10.nc sld012d_85_11.nc sld012d_85_12.nc sld012d_86_01.nc sld012d_86_02.nc sld012d_86_03.nc sld012d_86_04.nc sld012d_86_05.nc sld012d_86_06.nc sld012d_86_07.nc sld012d_86_08.nc sld012d_86_09.nc sld012d_86_10.nc sld012d_86_11.nc sld012d_86_12.nc sld012d_87_01.nc sld012d_87_02.nc sld012d_87_03.nc sld012d_87_04.nc sld012d_87_05.nc sld012d_87_06.nc sld012d_87_07.nc sld012d_87_08.nc sld012d_87_09.nc sld012d_87_10.nc sld012d_87_11.nc sld012d_87_12.nc sld012d_88_01.nc sld012d_88_02.nc sld012d_88_03.nc sld012d_88_04.nc sld012d_88_05.nc sld012d_88_06.nc sld012d_88_07.nc sld012d_88_08.nc sld012d_88_09.nc sld012d_88_10.nc sld012d_88_11.nc sld012d_88_12.nc sld012d_89_01.nc sld012d_89_02.nc sld012d_89_03.nc sld012d_89_04.nc sld012d_89_05.nc sld012d_89_06.nc sld012d_89_07.nc sld012d_89_08.nc sld012d_89_09.nc sld012d_89_10.nc sld012d_89_11.nc sld012d_89_12.nc sld012d_8589_0160.nc
mswrite -t 365 sld012d_8589_0160.nc /ZENDER/proc/sld012d/sld012d_8589_0160.nc &
rcp sld012d_8589_0160.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
ncwa -O -D 3 -a time -v SHFLX,LHFLX,FSNT,PRECT,PRECC,TOTCWP,CLDTOT,CLDHGH,CLDLOW,CLDMED,LWCF,SWCF,FLNT,FLNTC,ALBEDO,TS1,FLUS,GCLD,GCLR,gw /data2/zender/sld012d/sld012d_8589.nc /data2/zender/sld012d/sld012d_8589_foo.nc
ncdiff -O -D 3 -v SHFLX,LHFLX,FSNT,PRECT,PRECC,TOTCWP,CLDTOT,CLDHGH,CLDLOW,CLDMED,LWCF,SWCF,FLNT,FLNTC,ALBEDO,TS1,FLUS,GCLD,GCLR,gw /data2/zender/sld012d/sld012d_8589_0160.nc /data2/zender/sld012d/sld012d_8589_foo.nc /data2/zender/sld012d/sld012d_anom_8589_0160.nc
ncwa -A -D 3 -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/sld012d_anom_8589_0160.nc
/bin/rm -f /data2/zender/sld012d/sld012d_8589_foo.nc

# Create seasonal cycle files of zonal average anomalies
ncwa -O -D 3 -a lon -v LWCF,SWCF,gw /data2/zender/sld012d/sld012d_anom_8589_0112.nc /data2/zender/sld012d/sld012d_anom_xavg_8589_0112.nc

# Create interannual files of meridionally averaged anomalies
ncwa -O -D 3 -a lat -w gw -d lat,-10.,10. /data2/zender/sld012d/sld012d_anom_8589_0160.nc /data2/zender/sld012d/sld012d_anom_yavg_10S10N_8589_0160.nc
ncwa -O -D 3 -a lat -w gw -d lat,-05.,05. /data2/zender/sld012d/sld012d_anom_8589_0160.nc /data2/zender/sld012d/sld012d_anom_yavg_05S05N_8589_0160.nc

# Create interannual files of meridional averages
#foreach xpt (sld012d sld012d 422 spcp_85 erbe_b)
foreach xpt (sld012d 422 spcp_85 erbe_b)
echo "case $xpt"
ncwa -O -D 3 -a lat -w gw -d lat,-10.,10. /data2/zender/${xpt}/${xpt}_8589_0160.nc /data2/zender/${xpt}/${xpt}_yavg_10S10N_8589_0160.nc
ncwa -O -D 3 -a lat -w gw -d lat,-05.,05. /data2/zender/${xpt}/${xpt}_8589_0160.nc /data2/zender/${xpt}/${xpt}_yavg_05S05N_8589_0160.nc
end
ncwa -O -D 3 -a lat -w gw -d lat,-10.,10. /data2/zender/sld012d/sld012d_8589_0160.nc /data2/zender/sld012d/sld012d_yavg_10S10N_8589_0160.nc
ncwa -O -D 3 -a lat -w gw -d lat,-05.,05. /data2/zender/sld012d/sld012d_8589_0160.nc /data2/zender/sld012d/sld012d_yavg_05S05N_8589_0160.nc

# Create the difference files for 1987-1985 (all 12-months) then hack on those locally
ncrcat -O -R -D 3 -v LWCF,SWCF,TS1,TMQ,CMFMC,QC,CLDTOT,CLDHGH,gw -l ./ -p /ZENDER/proc/sld012d -n 12,2,1 sld012d_85_01.nc sld012d_85_0112.nc
ncrcat -O -R -D 3 -v LWCF,SWCF,TS1,TMQ,CMFMC,QC,CLDTOT,CLDHGH,gw -l ./ -p /ZENDER/proc/sld012d -n 12,2,1 sld012d_87_01.nc sld012d_87_0112.nc
mswrite -t 365 sld012d_85_0112.nc /ZENDER/proc/sld012d/sld012d_85_0112.nc &
mswrite -t 365 sld012d_87_0112.nc /ZENDER/proc/sld012d/sld012d_87_0112.nc &
rcp sld012d_??_0112.nc heinlein.cgd.ucar.edu:/data2/zender/sld012d &
# NB: ncdiff has a bug when hyperslabbing while differencing, so difference then hyperslab for now.
ncdiff -O -D 3 -l ./ -p /ZENDER/proc/sld012d sld012d_87_0112.nc sld012d_85_0112.nc sld012d_87m85_0112.nc
# Make seasonal cycle files
ncks -O -D 3 -d ilev,500. sld012d_87m85_0112.nc foo.nc
ncwa -O -D 3 -a ilev /data2/zender/sld012d/foo.nc /data2/zender/sld012d/foo.nc
ncks -O -D 3 -C -x -v ilev /data2/zender/sld012d/foo.nc /data2/zender/sld012d/foo.nc
ncwa -A -D 3 -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/foo.nc 
ncecat -O -D 3 /data2/zender/sld012d/foo.nc /data2/zender/sld012d/foo_ocn_87m85_0112.nc
ncwa -O -D 3 -a record -m ORO -M 0. -o eq /data2/zender/sld012d/foo_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc

# Make special averaged files
/bin/cp /data2/zender/sld012d/sld012d_87m85_0112.nc /data2/zender/sld012d/foo.nc
ncwa -A -D 3 -C -a time -v ORO /data2/zender/ccm/SEP1.nc /data2/zender/sld012d/foo.nc
#ncwa -O -D 3 -m ORO -M 0. -o eq -a time,lat -w gw -d time,2,4 -d lat,-10.,10. /data2/zender/sld012d/foo.nc /data2/zender/sld012d/sld012d_yavg_10S10N_87m85_MAM.nc
ncwa -O -D 3 -a time,lat -w gw -d time,2,4 -d lat,-10.,10. /data2/zender/sld012d/foo.nc /data2/zender/sld012d/sld012d_yavg_10S10N_87m85_MAM.nc
/bin/rm *foo*

ncks -O -d lat,-5.,5. -d lon,175.,185. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_tst_87m85_0112.nc

ncks -O -d lat,-20.,20. -d lon,140.,180. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Tropical_Western_87m85_0112.nc
ncks -O -d lat,-20.,20. -d lon,180.,230. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Tropical_Central_87m85_0112.nc
ncks -O -d lat,-20.,20. -d lon,230.,270. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Tropical_Eastern_87m85_0112.nc
ncks -O -d lat,-20.,20. -d lon,140.,270. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Tropical_87m85_0112.nc

ncks -O -d lat,-10.,10. -d lon,140.,180. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Equatorial_Western_87m85_0112.nc
ncks -O -d lat,-10.,10. -d lon,180.,230. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Equatorial_Central_87m85_0112.nc
ncks -O -d lat,-10.,10. -d lon,230.,270. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Equatorial_Eastern_87m85_0112.nc
ncks -O -d lat,-10.,10. -d lon,140.,270. /data2/zender/sld012d/sld012d_ocn_87m85_0112.nc /data2/zender/sld012d/sld012d_reg_Pacific_Equatorial_87m85_0112.nc




