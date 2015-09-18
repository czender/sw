# $Id$

# Purpose: Generate data for CAM/SNICAR work

# Testing:
# /bin/cp -r ${DATA}/sncpd02 ${DATA}/sncpd02_bck
# /bin/rm -r ${DATA}/sncpd02;/bin/cp -r ${DATA}/sncpd02_bck ${DATA}/sncpd02

# Production:
# cd ~/anl;anl_cam_snc.sh PD
# cd ~/anl;chmod a+x anl_cam_snc.sh;anl_cam_snc.sh
# cd ~/anl;time anl_cam_snc.sh > foo 2>&1 &

# Uncomment next line to trace program execution
#set -o xtrace # dbg

# 0. Command-line options
#clm_id='LGM' # [sng] Last Glacial Maximum
#clm_id='PI'  # [sng] Pre-industrial
clm_id='PD'  # [sng] Present day
#clm_id='2050A2'  # [sng] IPCC A2 Scenario 2050

if [ -n "${1}" ]; then
    clm_id=${1}
fi # !$1

if [ ${clm_id} = 'LGM' ]; then
    clm_sng='lgm';caseid_xpt='snclgm05';caseid_ctl='snclgm06' # LGM
elif [ ${clm_id} = 'PI' ]; then
    clm_sng='pi';caseid_xpt='cssncpi03b';caseid_ctl='cssncpi04b' # PI
elif [ ${clm_id} = 'PD' ]; then
    clm_sng='pd';caseid_xpt='sncpd05';caseid_ctl='sncpd06' # PD
elif [ ${clm_id} = '2050A2' ]; then
    clm_sng='2050a2';caseid_xpt='cssnc2050_02b';caseid_ctl='cssnc2050_03b' # 2050A2
fi # end if

# Print initial state
if [ '1' = '1' ]; then # dbg
    printf "${0}: Initialization State:\n"
    printf "${0}: clm_id=${clm_id}\n"
    printf "${0}: clm_sng=${clm_sng}\n"
    printf "${0}: caseid_xpt=${caseid_xpt}\n"
    printf "${0}: caseid_ctl=${caseid_ctl}\n"
fi # !dbg

# 1. Primary variables

# 2. Derived variables

fld_lev='VLWCSFC+,WATSAT,lev+,TSOI+,TLAKE+,SOILLIQ,SOILICE,H2OSOI,QCHANR,QCHOCNR,latrot,lonrof'
# ncks -x -v ${fld_lev} -H ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_clm_xyt.nc | /bin/more
for caseid in ${caseid_ctl} ${caseid_xpt}; do
    ncwa -O ${DATA}/anl_${caseid}/${caseid}_clm_xyt.nc ${DATA}/anl_${caseid}/${caseid}_clm_txyz.nc
    ncap -O -s "FTNS=FSNS-FLNS" ${DATA}/anl_${caseid}/${caseid}_clm_txyz.nc ${DATA}/anl_${caseid}/${caseid}_clm_txyz.nc
done # end ctl/xpt loop
ncap -O -s "SNOBCFRC2L=0.0" -s "SNODSTFRC2L=0.0" -s "SNOOCFRC2L=0.0" -s "SNOBCFRCI=0.0" -s "SNODSTFRCI=0.0" -s "SNOBCFRCL=0.0" -s "SNODSTFRCL=0.0" -s "SNOOCFRCL=0.0" ${DATA}/anl_${caseid_ctl}/${caseid_ctl}_clm_txyz.nc ${DATA}/anl_${caseid_ctl}/${caseid_ctl}_clm_txyz.nc
ncdiff -O ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_clm_txyz.nc ${DATA}/anl_${caseid_ctl}/${caseid_ctl}_clm_txyz.nc ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_${caseid_ctl}_clm_txyz.nc
ncdiff -O -v TREFHT,SNOAERFRC,SNOBCFRC,SNODSTFRC,FTNT ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_clm_txyz.nc ${DATA}/anl_${caseid_ctl}/${caseid_ctl}_clm_txyz.nc ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_${caseid_ctl}_clm_txyz.nc
ncks -H -v TREFHT,FTNT ${DATA}/anl_${caseid_xpt}/${caseid_xpt}_${caseid_ctl}_clm_txyz.nc
# Diagnose Efficacy
ncap -O -v \
-s "SFCFRC2TRPFRC=0.91; // [frc] FZR06 ZBP97" \
-s "CO2SNS=0.69; // [K (W m-2)-1] KSH06" \
-s "SNOAERFRC_TOA=SFCFRC2TRPFRC*SNOAERFRC" \
-s "SNOAERSNS=TREFHT/SNOAERFRC_TOA" \
-s "SNOAERFFC=SNOAERSNS/CO2SNS" \
-s "SNOBCSNS=0.15/(0.054*SFCFRC2TRPFRC); // FZR07 BC-Only sensitivity" \
-s "SNOBCFFC=SNOBCSNS/CO2SNS" \
-s "SNODSTFRC_TOA=SFCFRC2TRPFRC*SNODSTFRC" \
-s "SNOBCFRC_TOA=SFCFRC2TRPFRC*SNOBCFRC" \
-s "SNODSTSNS=(SNOAERFRC_TOA*SNOAERSNS-SNOBCFRC_TOA*SNOBCSNS)/SNODSTFRC_TOA; // Climate sensitivity to dust computed as residual of total and FZR07 BC-only changes" \
-s "SNODSTFFC=SNODSTSNS/CO2SNS" \
${DATA}/anl_${caseid_xpt}/${caseid_xpt}_${caseid_ctl}_clm_txyz.nc ~/foo_${clm_sng}.nc
# Temperature change due to both aerosols is 0.154 K, 2.75 K (PD, LGM) (snc???05/06)
# Temperature change due to both aerosols is 0.223 K, 2.49 K (PD, LGM) (snc???02/03)
# Temperature change due to both aerosols is 0.163 K, 0.952 K (PD, LGM) (cspdfXX)
# Kiehl et al. (2006) have CO2 doubling sensitivity as 0.69 K (W m-2)-1 based on
# Regressions of 20 years of CAMSOM for doubled (355->710) CO2 forcing of 3.58 W m-2
# causing temperature change of 2.47 K (equilibrium response) 
# (1.48 K is transient response centered on year 70)
ncks -H ~/foo_${clm_sng}.nc

if [ '0' = '1' ]; then # dbg
# Diagnose surface deposition fluxes and mass paths from dry fluxes and wet tendencies
    ncrcat -O -D 1 -v 'DSTX0[1-4]PP,DSTX0[1-4]DD,DSTX0[1-4]SF,^DSTX0[1-4]$,PDELDRY,gw' ${DATA}/${caseid_xpt}/${caseid_xpt}.cam2.h0.????-[0-1][0-9].nc ${DATA}/tmp/${caseid_xpt}_sfc_flx_4D_foo.nc
    ncap -O -D 1 -s "DSTX01DW=DSTX01PP*PDELDRY/9.8" -s "DSTX02DW=DSTX02PP*PDELDRY/9.8" -s "DSTX03DW=DSTX03PP*PDELDRY/9.8" -s "DSTX04DW=DSTX04PP*PDELDRY/9.8" -s "MPDSTX01=DSTX01*PDELDRY/9.8" -s "MPDSTX02=DSTX02*PDELDRY/9.8" -s "MPDSTX03=DSTX03*PDELDRY/9.8" -s "MPDSTX04=DSTX04*PDELDRY/9.8" ${DATA}/tmp/${caseid_xpt}_sfc_flx_4D_foo.nc ${DATA}/tmp/${caseid_xpt}_sfc_flx_4D_prc_foo.nc
    ncwa -O -D 1 -N -a lev ${DATA}/tmp/${caseid_xpt}_sfc_flx_4D_prc_foo.nc ${DATA}/tmp/${caseid_xpt}_sfc_flx_3D_foo.nc
    ncra -O -D 1 ${DATA}/tmp/${caseid_xpt}_sfc_flx_3D_foo.nc ${DATA}/tmp/${caseid_xpt}_sfc_flx_2D_foo.nc
    ncap -O -D 1 -s "DSTSFWET=DSTX01DW+DSTX02DW+DSTX03DW+DSTX04DW" -s "DSTSFDRY=DSTX01DD+DSTX02DD+DSTX03DD+DSTX04DD" -s "DSTSFMBL=DSTX01SF+DSTX02SF+DSTX03SF+DSTX04SF" -s "DSTSFDPS01=DSTX01DW+DSTX01DD" -s "DSTSFDPS02=DSTX02DW+DSTX02DD" -s "DSTSFDPS03=DSTX03DW+DSTX03DD" -s "DSTSFDPS04=DSTX04DW+DSTX04DD" -s "DSTSFDPS=DSTSFDRY+DSTSFWET" -s "DSTMPC=MPDSTX01+MPDSTX02+MPDSTX03+MPDSTX04" ${DATA}/tmp/${caseid_xpt}_sfc_flx_2D_foo.nc ${DATA}/tmp/${caseid_xpt}_sfc_flx.nc
    ncwa -O -D 1 -a lat,lon -w gw ${DATA}/tmp/${caseid_xpt}_sfc_flx.nc ${DATA}/tmp/${caseid_xpt}_sfc_flx_xyt.nc
    scp esmf.ess.uci.edu:/data/zender/tmp/sncpd05_sfc_flx.nc esmf.ess.uci.edu:/data/zender/tmp/sncpd05_sfc_flx_xyt.nc /data/zender/tmp
fi # endif dbg

if [ '0' = '1' ]; then # dbg
# Generate annual mean files of whole experiments:
    xpt='snclgm05 snclgm06 cssncpi03b cssncpi04b sncpd05 sncpd06 cssnc2050_02b cssnc2050_03b'
    for caseid in ${xpt}; do
	cd ${DATA}/${caseid}
	printf "Processing climatology for caseid = ${caseid}...\n"
	ncra -O -D 1 ${caseid}.cam2.h0.000[6-9]-??.nc ${caseid}.cam2.h0.00[1-2]?-??.nc ${DATA}/tmp/${caseid}.cam2_clm.nc
	ncra -O -D 1 ${caseid}.clm2.h0.000[6-9]-??.nc ${caseid}.clm2.h0.00[1-2]?-??.nc ${DATA}/tmp/${caseid}.clm2_clm.nc
    done # end loop over caseid
    cd ${DATA}/tmp
    rsync 'esmf.ess.uci.edu:/data/zender/tmp/snc*_clm.nc' ${DATA}/tmp
fi # endif dbg

if [ '0' = '1' ]; then # dbg
# Collect annual mean files:
    xpt='snclgm05 snclgm06 cssncpi03b cssncpi04b sncpd05 sncpd06 cssnc2050_02b cssnc2050_03b'
    for caseid in ${xpt}; do
	printf "Grabbing climatology for caseid = ${caseid}...\n"
	mkdir -p ${DATA}/anl_${caseid}
	cd ${DATA}/anl_${caseid}
	rsync sand.ess.uci.edu:/data/zender/anl_${caseid} ${DATA}
#	rsync sand.ess.uci.edu:/data/zender/anl_${caseid}/${caseid}_clm.nc ${DATA}/anl_${caseid}
#	rsync sand.ess.uci.edu:/data/zender/anl_${caseid}/${caseid}_clm_xyt.nc ${DATA}/anl_${caseid}
    done # end loop over caseid
fi # endif dbg

if [ '0' = '1' ]; then # dbg
# Difference control from experiment
    xpt='snclgm05'
    ctl='snclgm06'
    ncdiff -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
    xpt='sncpd05'
    ctl='sncpd06'
    ncdiff -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
    xpt='cssncpi03b'
    ctl='cssncpi04b'
    ncdiff -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
    xpt='cssnc2050_02b'
    ctl='cssnc2050_03b'
    ncdiff -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 

# Difference present climate from other climates
    xpt='snclgm06'
    ctl='sncpd06'
    var_xcl=''
    ncdiff ${var_xcl} -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
    xpt='cssncpi04b'
    ctl='sncpd06'
    var_pd_not_pi='ALBNS_SNO,ALBS_SNO,ALBVS_SNO,FSDSN_SNO,FSDSV_SNO,FSDS_SNO,SABGN_SNO,SABGV_SNO,SABG_SNO,SABG_SNO0,SABG_SNO1,SABG_SNO2,SABG_SNO3,SABG_SNO4,SNOAGE,SNOAGEI,SNOAGEL,SNOLIQI,TSNO0,TSNO1,TSNO2,TSNO3,TSNO4,'
    var_pi_not_pd='SNOFSDSND,SNOFSDSNI,SNOFSDSVD,SNOFSDSVI,SNOFSRND,SNOFSRNI,SNOFSRVD,SNOFSRVI,SNOTTOPL'
    var_xcl="-x -v ${var_pd_not_pi}${var_pi_not_pd}"
    ncdiff ${var_xcl} -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
    xpt='cssnc2050_03b'
    ctl='sncpd06'
    ncdiff ${var_xcl} -O ${DATA}/anl_${xpt}/${xpt}_clm.nc ${DATA}/anl_${ctl}/${ctl}_clm.nc ${DATA}/anl_${xpt}/${xpt}_${ctl}_clm.nc 
fi # endif dbg

if [ '0' = '1' ]; then # dbg
# Convert double to single precision and merge ice+land fields fxm: do in anl_* batch
    xpt='snclgm05 snclgm06 cssncpi03b cssncpi04b sncpd05 sncpd06 cssnc2050_02b cssnc2050_03b'
    for caseid in ${xpt}; do
	printf "Converting double to single precision and merging ice+land fields for caseid = ${caseid}...\n"
	cd ${DATA}/anl_${caseid}
	ncap -O -s 'FSNOI=float(FSNOI)' ${DATA}/anl_${caseid}/${caseid}_clm.nc ${DATA}/anl_${caseid}/${caseid}_clm.nc
	ncatted -O -a missing_value,FSNOI,o,f,1.0e36 -a _FillValue,FSNOI,o,f,1.0e36 ${DATA}/anl_${caseid}/${caseid}_clm.nc
#	ncks -m -v FSNOI ${DATA}/anl_${caseid}/${caseid}_clm.nc
	ncap2 -O -s 'FSNOL=FSNO;FSNO=landfrac*FSNOL+ICEFRAC*FSNOI' ${DATA}/anl_${caseid}/${caseid}_clm.nc ${DATA}/anl_${caseid}/${caseid}_clm.nc
    done # end loop over caseid
fi # endif dbg
