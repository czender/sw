# Purpose: Process dust map/erodibility data

# Usage: cd ~/map;chmod a+x bds.sh;${HOME}/map/bds.sh > ~/map.txt 2>&1 &
# scp ~/map/bds.sh dust.ess.uci.edu:map

# URL: http://dust.ess.uci.edu/dead/data

caseid='map'
#drc_bds="${DATA}/data" # Directory for bds output
drc_bds='/var/www/html/dead/data' # Directory for bds output when flg_cst_tpr=T
#drc_bds='/var/www/html/dead/data/nocoast' # Directory for bds output when flg_cst_tpr=F
#drc_bds='/var/www/html/dead/data/norgnxcl' # Directory for bds output when flg_rgn_xcl=F

cd ${HOME}/map # Must be in 'map' directory to access local region mask files

CCSM_PST_PRC_FLG='1' # [flg] Post-process CCSM Gaussian and FV grids
CCSM_GNR_FLG='1' # [flg] Generate CCSM Gaussian and FV grids
FALSE='0' # [flg] Flag that is never true

if [ ${CCSM_GNR_FLG} = '1' ]; then
# Create dust source maps for CAM3/CLM Gaussian and FV grids
# Note that NCAR uses "lat by lon" = "latxlon" = "Y by X" nomenclature
# On their grids, this is the same as "coarser dimension by finer dimension"
# FV grids have map_typ = 33

# 2.5x3.3 fxm
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/ | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 108 -y 72 -m 33 -t 19801980 -o ${drc_bds}/dst_2.5x3.3.nc

# 2x2.5 /fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-09-01_2x2.5_L26_c040615.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/cami_0000-09-01_2x2.5_L26_c020430.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 144 -y 91 -m 33 -t 19801980 -o ${drc_bds}/dst_2x2.5.nc

# 1.9x2.5 /fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-09-01_1.9x2.5_L26_c040809.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/ | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 144 -y 96 -m 33 -t 19801980 -o ${drc_bds}/dst_1.9x2.5.nc

# 1x1.25 /fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-09-01_1x1.25_L26_c030918.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/cami_0000-09-01_1x1.25_L26_c030918.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 288 -y 181 -m 33 -t 19801980 -o ${drc_bds}/dst_1x1.25.nc

# 1x1 fxm
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 360 -y 181 -m 33 -t 19801980 -o ${drc_bds}/dst_1x1.nc

# 0.9x1.25 dst_0.9x1.25_c080724.nc
# ncks -H -v lon,lat,slat,slon -m /home/zender/Desktop/dst_0.9x1.25_c080724.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 288 -y 192 -m 33 -t 19801980 -o ${drc_bds}/dst_0.9x1.25.nc

# 0.47x0.63 dst_0.47x0.63_c081217.nc
# ncks -H -v lon,lat,slat,slon -m /home/zender/Desktop/dst_0.47x0.63_c081217.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 576 -y 384 -m 33 -t 19801980 -o ${drc_bds}/dst_0.47x0.63.nc

# 0.23x0.31 dst_0.23x0.31_c090117.nc
# ncks -H -v lon,lat,slat,slon -m /home/zender/Desktop/dst_0.23x0.31_c090117.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 1152 -y 768 -m 33 -t 19801980 -o ${drc_bds}/dst_0.23x0.31.nc

# 10x15 /fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-09-01_10x15_L26_c020430.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/cami_0000-09-01_10x15_L26_c020430.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 24 -y 19 -m 33 -t 19801980 -o ${drc_bds}/dst_10x15.nc

# 4x5 /fs/cgd/csm/inputdata/atm/cam/inic/fv/cami_0000-09-01_4x5_L26_c031217.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/fv/cami_0000-09-01_4x5_L26_c031217.nc | m
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 72 -y 46 -m 33 -t 19801980 -o ${drc_bds}/dst_4x5.nc

# Spectral grids: map_typ = 23
# 32x64 /fs/cgd/csm/inputdata/atm/cam/inic/gaus/cami_0000-09-01_32x64_L26_c030918.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/gaus/cami_0000-01-01_32x64_L26_c030228.nc
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 64 -y 32 -m 23 -t 19801980 -o ${drc_bds}/dst_32x64.nc

# 48x96 (i.e., T31) /fs/cgd/csm/inputdata/atm/cam/inic/gaus/cami_0000-09-01_48x96_L26_c040420.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/gaus/cami_0000-09-01_48x96_L26_c040420.nc
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 96 -y 48 -m 23 -t 19801980 -o ${drc_bds}/dst_48x96.nc

# 64x128 (i.e., T42) /fs/cgd/csm/inputdata/atm/cam/inic/gaus/
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/gaus/cami_0000-09-01_64x128_L26_c030918.nc
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 128 -y 64 -m 23 -t 19801980 -o ${drc_bds}/dst_64x128.nc

# 94x192 (i.e., T62)
# ncks -H -v lon,lat,slat,slon -m 
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 192 -y 94 -m 23 -t 19801980 -o ${drc_bds}/dst_94x192.nc

# 128x256  (i.e., T85) /fs/cgd/csm/inputdata/atm/cam/inic/gaus/cami_0000-01-01_128x256_L26_c030918.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/gaus/cami_0000-01-01_128x256_L26_c030918.nc
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 256 -y 128 -m 23 -t 19801980 -o ${drc_bds}/dst_128x256.nc

# 8x16 /fs/cgd/csm/inputdata/atm/cam/inic/gaus/cami_0000-01-01_8x16_L26_c030228.nc
# ncks -H -v lon,lat,slat,slon -m /datashare/inputdata/csm/atm/cam2/inic/gaus/cami_0000-01-01_8x16_L26_c030228.nc
bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 16 -y 8 -m 23 -t 19801980 -o ${drc_bds}/dst_8x16.nc
fi # !CCSM_GNR_FLG

if [ ${CCSM_PST_PRC_FLG} = '1' ]; then
# Strip boundary data sets of everything but erodibility factor
#for rsn in 2.5x3.3 2x2.5 1.9x2.5 1x1.25 1x1 0.9x1.25 0.47x0.63 0.23x0.31 10x15 4x5 32x64 48x96 64x128 94x192 128x256 8x16; do
for rsn in 1x1; do
# Following factors are tuned for T62
    ncap2 -O -s "rdb_fct_unity=0.0*mbl_bsn_fct+0.1614" -s "rdb_fct_tpg=mbl_bsn_fct" -s "rdb_fct_gmr=5.707*sfc_acm_fct" -s "rdb_fct_hyd=22.429*flw_acm_fct" -s "rdb_fct_rfl_mds_lnr=rfl_sfc_lnd_nrm_mds_lnr" -s "rdb_fct_rfl_mds_sqr=rfl_sfc_lnd_nrm_mds_sqr" ${drc_bds}/dst_${rsn}.nc ${drc_bds}/dst_${rsn}.nc
    ncks -O -v area,bsn_enm,bsn_sz,lat_grd,lon_grd,lat,lat_fv,lon,lon_sz,lat_sz,oro,lnd_frc,lnd_msk,mss_frc_cly,mss_frc_slt,mss_frc_snd,rdb_fct_gmr,rdb_fct_unity,rdb_fct_tpg,rdb_fct_hyd,rdb_fct_rfl_mds_lnr,rdb_fct_rfl_mds_sqr,slat,slon,vai_ttl_clm ${drc_bds}/dst_${rsn}.nc ${drc_bds}/dst_${rsn}.nc
    ncap2 -O -s "mbl_bsn_fct=rdb_fct_gmr" ${drc_bds}/dst_${rsn}.nc ${drc_bds}/dst_${rsn}.nc
#   ncks -H -v lon,lat,slat,slon -m ${drc_bds}/dst_${rsn}.nc | m
#    scp ${drc_bds}/dst_${rsn}.nc dust.ess.uci.edu:/var/www/html/dead/data/nocoast
#    scp ${drc_bds}/dst_${rsn}.nc dust.ess.uci.edu:/var/www/html/dead/data
#    cp ${drc_bds}/dst_${rsn}.nc /var/www/html/dead/data/nocoast
done
fi # !CCSM_PST_PRC

if [ ${FALSE} = '1' ]; then
# Generate most frequently used dust boundary data sets
bds --bsn_fct_hrz=F -x 180 -y 90 -m 13 -t 19801980 -o ${drc_bds}/dst_2x2.nc
bds --bsn_fct_hrz=F -x 16 -y 8 -m 23 -t 19801980 -o ${drc_bds}/dst_T5.nc
bds --bsn_fct_hrz=F -x 96 -y 48 -m 23 -t 19801980 -o ${drc_bds}/dst_T31.nc
bds --bsn_fct_hrz=F -x 128 -y 64 -m 23 -t 19801980 -o ${drc_bds}/dst_T42.nc
bds --bsn_fct_hrz=F -x 192 -y 94 -m 23 -t 19801980 -o ${drc_bds}/dst_T62.nc
bds --bsn_fct_hrz=F -x 320 -y 160 -m 23 -t 19801980 -o ${drc_bds}/dst_T106.nc
bds --bsn_fct_hrz=F -x 360 -y 180 -m 13 -t 19801980 -o ${drc_bds}/dst_1x1.nc
bds --bsn_fct_hrz=F -x 384 -y 190 -m 23 -t 19801980 -o ${drc_bds}/dst_T126.nc
bds --bsn_fct_hrz=F -x 512 -y 256 -m 23 -t 19801980 -o ${drc_bds}/dst_T170.nc

# Make basin factors on high resolution grid
bds --bsn_fct_hrz=T -x 96 -y 48 -m 23 -t 19801980 -o ${drc_bds}/dst_T31_hrz.nc
bds --bsn_fct_hrz=T -x 128 -y 64 -m 23 -t 19801980 -o ${drc_bds}/dst_T42_hrz.nc
bds --bsn_fct_hrz=T -x 192 -y 94 -m 23 -t 19801980 -o ${drc_bds}/dst_T62_hrz.nc
bds --bsn_fct_hrz=T -x 360 -y 180 -m 13 -t 19801980 -o ${drc_bds}/dst_1x1_hrz.nc

# Strip boundary data sets of unnecessary fields to conserve size
for rsn in 1x1 2x2 T5 T31 T42 T62 T106 T126 T170 T31_hrz T42_hrz T62_hrz 1x1_hrz ; do
# Leave 2-D field that alphabetically precedes gw so IBP does not try to load gw
    ncks -O -v gw,fsh_fct,mss_frc_snd,mss_frc_cly,lnd_frc_dry,sfc_typ,time,mbl_bsn_fct,erd_fct_unity,erd_fct_tpg,erd_fct_gmr,erd_fct_hyd,lai_lsm,vai_lsm,vai_ttl_clm,mss_frc_CaCO3,sfc_frc_bln ${drc_bds}/dst_${rsn}.nc ${drc_bds}/dst_${rsn}.nc
    scp ${drc_bds}/dst_${rsn}.nc dust.ess.uci.edu:/var/www/html/tmp
done

# Split bds output into months (helps with IBP viewing)
for mth in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
    echo "Processing caseid ${caseid} month ${mth}"
    ncks -O -F -d time,${mth} ${wrk_drc}/${caseid}_8589_0112.nc ${wrk_drc}/${caseid}_8589_$mth.nc
done
fi # !FALSE
