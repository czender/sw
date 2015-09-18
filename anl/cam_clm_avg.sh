# $Id$

# Purpose: Average and hyperslab monthly data into climatological and regional data

# Testing:
# /bin/cp -r ${DATA}/sncpd02 ${DATA}/sncpd02_bck
# /bin/rm -r ${DATA}/sncpd02;/bin/cp -r ${DATA}/sncpd02_bck ${DATA}/sncpd02

# Distribution:
# scp ~/anl/cam_clm_avg.sh dust.ess.uci.edu:anl

# Production:
# cd ~/anl;cam_clm_avg.sh
# cd ~/anl;time ./cam_clm_avg.sh > foo.avg 2>&1 &
# cd ~/anl;time ./cam_clm_avg.sh sncpd02 0000 0020

# set -o xtrace # dbg

# Analysis script for CAM/CLM output.
# Produces:
#   -climate mean
#   -monthly means 
#   -seasonal means 
#   -ensemble annual mean
#   -ensemble seasonal means
#   -global means of all of the above
#   -select regional means

# User sets following environmental variables:
#  caseid:  caseid ( original data MUST be in directory ${DATA}/${caseid} )
#  drc_out: Output directory for processed, reduced data. 
#           Script does not modify original data. 
#           User needs write permission to ${DATA}/${caseid}
#  CC:      Two-digit century, i.e., first two digits of year in caseid. 
#           NOTE: Change convention if analyzing more than 100 years of data
#  yyyy_srt: Year of first December to analyze.
#           Convention is to include Dec of yyyy_srt in seasonal 
#           and monthly analysis, and exclude Dec of yyyy_end. 
#           All other months in yyyy_srt are excluded.
#           Annual means and climate means exclude all months of yyyy_srt,
#           and include all months of yyyy_end.  
#           Seasonal averages are not 100% consistent with annual averages.  
#           But to have continuity in the DJF analysis, there is no alternative.
#  yyyy_end: Last year of analysis. 
#           Dec of yyyy_end is excluded from the seasonal and monthly analysis.  
#           For a proper analysis, yyyy_end must be a complete year of output.
#  MDL:     This tag determines whether CLM or CAM or (other?) output is analyzed
#           It refers to the character sequence in the history files:
#           caseid.MDL.h0.YR-MTH.nc (valid examples are 'clm2' and 'cam2')        
#           NOTE: right now 'h0' is fixed in the file name, but an environmental 
#           variable could easily replace this to analyze other history tapes.

function trim_leading_zeros {
# Purpose: Trim leading zeros from string representing an integer
# Bash treats zero-padded integers as octal
# This is surprisingly hard to workaround
# My way around this is to remove leading zeros prior to arithmetic
# Usage: trim_leading zeros ${sng}
    sng_trm=${1} # [sng] Trimmed string
# Remove up to three leading zeros, one at a time, using bash 2.0 pattern matching
    sng_trm=${sng_trm##0} # NeR98 p. 99
    sng_trm=${sng_trm##0} # NeR98 p. 99
    sng_trm=${sng_trm##0} # NeR98 p. 99
# If all zeros removed, replace with single zero
    if [ ${sng_trm} = '' ]; then 
	sng_trm='0'
    fi # endif
} # end trim_leading_zeros

# 0. Command-line options
caseid='sncpd02'
yyyy_srt='0000' # Start year (NB: Include year 0 for now)
yyyy_end='0020' # End year

if [ -n "${1}" ]; then
    caseid=${1}
fi # !$1
if [ -n "${2}" ]; then
    yyyy_srt=${2}
fi # !$2
if [ -n "${3}" ]; then
    yyyy_end=${3}
fi # !$3

# 1. Primary variables
MDL='clm2'

# 2. Derived variables
drc_in=${DATA}/${caseid}
drc_out=${DATA}/anl_${caseid}

# 3. Doubly-derived variables
wrk_dir=${drc_out}

# Determine century and first full year
trim_leading_zeros ${yyyy_srt}
yr_srt_rth=${sng_trm}
let CC_tmp=${yr_srt_rth}/100
CC=`printf "%02d\n" ${CC_tmp}`
let yr_one=${yr_srt_rth}+1
#let yr_one=`printf "%02d\n" ${yr_srt_rth}`+1
printf "${0}: yyyy_srt = ${yyyy_srt}\n"
printf "${0}: yr_srt_rth = ${yr_srt_rth}\n"
printf "${0}: CC_tmp = ${CC_tmp}\n"
printf "${0}: CC = ${CC}\n"
printf "${0}: yr_one = ${yr_one}\n"

# Create output directory
mkdir -p ${drc_out}
mkdir -p ${drc_in}/tmp

# Go to working directory
cd ${wrk_dir}

# Move any existing ${drc_in}/tmp files to original directories (precautionary)
mv ${drc_in}/tmp/* ${drc_in}

# Move data from year 0 out of the way
mv ${drc_in}/${caseid}.${MDL}.h0.${yyyy_srt}-[0-1][0-9].nc ${drc_in}/tmp

# STEP 1: Create ensemble annual and seasonal means
# Bring in December of year 0 for seasonal mean
mv ${drc_in}/tmp/${caseid}.${MDL}.h0.${yyyy_srt}-12.nc ${drc_in}
for yr in `seq $yr_one $yyyy_end`; do
    YYYY=`printf "%04d" ${yr}`
    let yrm1=yr-1
    YYYYm1=`printf "%04d" ${yrm1}`

    # Package all months from each year into a timeseries and apply dpm weights:
    printf "Creating ensemble annual and seasonal means for year ${YYYY}...\n"
    ncrcat -O ${drc_in}/${caseid}.${MDL}.h0.${YYYY}-[0-1][0-9].nc ${drc_out}/foo.nc
    ncwa -O -x -v time_written,date_written -a time -b -w dpm ${drc_out}/foo.nc ${drc_out}/${caseid}_${YYYY}.nc

    # Package all months from the year, plus December from previous year,
    # and average over appropriate hyper-slab, applying dpm weights
    ncrcat -O ${drc_in}/${caseid}.${MDL}.h0.${YYYYm1}-12.nc ${drc_in}/${caseid}.${MDL}.h0.${YYYY}-[0-1][0-9].nc ${drc_out}/foo2.nc
    ncwa -O -x -v time_written,date_written -a time -b -w dpm -d time,0,2 ${drc_out}/foo2.nc ${drc_out}/${caseid}_${YYYY}_DJF.nc
    ncwa -O -x -v time_written,date_written -a time -b -w dpm -d time,3,5 ${drc_out}/foo2.nc ${drc_out}/${caseid}_${YYYY}_MAM.nc
    ncwa -O -x -v time_written,date_written -a time -b -w dpm -d time,6,8 ${drc_out}/foo2.nc ${drc_out}/${caseid}_${YYYY}_JJA.nc
    ncwa -O -x -v time_written,date_written -a time -b -w dpm -d time,9,11 ${drc_out}/foo2.nc ${drc_out}/${caseid}_${YYYY}_SON.nc
done
rm -f ${drc_out}/foo.nc
rm -f ${drc_out}/foo2.nc

# Package annual means into time series:
printf "Packaging annual means into time series...\n"
ncrcat -O ${drc_out}/${caseid}_${CC}[0-9][0-9].nc ${drc_out}/${caseid}_ts_ANN.nc
# Remove originals:
rm -f ${drc_out}/${caseid}_${CC}[0-9][0-9].nc

# Package seasonal means into time series:
printf "Packaging seasonal means into time series...\n"
ncrcat -O ${drc_out}/${caseid}_${CC}[0-9][0-9]_DJF.nc ${drc_out}/${caseid}_ts_DJF.nc
ncrcat -O ${drc_out}/${caseid}_${CC}[0-9][0-9]_MAM.nc ${drc_out}/${caseid}_ts_MAM.nc
ncrcat -O ${drc_out}/${caseid}_${CC}[0-9][0-9]_JJA.nc ${drc_out}/${caseid}_ts_JJA.nc
ncrcat -O ${drc_out}/${caseid}_${CC}[0-9][0-9]_SON.nc ${drc_out}/${caseid}_ts_SON.nc
# Remove originals
rm -f ${drc_out}/${caseid}_${CC}[0-9][0-9]_DJF.nc
rm -f ${drc_out}/${caseid}_${CC}[0-9][0-9]_MAM.nc
rm -f ${drc_out}/${caseid}_${CC}[0-9][0-9]_JJA.nc
rm -f ${drc_out}/${caseid}_${CC}[0-9][0-9]_SON.nc

# STEP 2: Create monthly means
# Move Dec of last year out (Dec of year 0 is already in)
# so ensemble seasonal means are consistent with monthly means:
mv ${drc_in}/${caseid}.${MDL}.h0.${yyyy_end}-12.nc ${drc_in}/tmp/
for mth in {01..12}; do
    MM=`printf "%02d" $mth`
    printf "Creating monthly means for month ${MM}...\n"
    ncra -O -x -v time_written,date_written ${drc_in}/${caseid}.${MDL}.h0.${CC}[0-9][0-9]-${MM}.nc ${drc_out}/${caseid}_clm_${MM}.nc
done 

# Package monthly means into a time series:
printf "Packaging monthly means into time series...\n"
ncrcat -O ${drc_out}/${caseid}_clm_[0-1][0-9].nc ${drc_out}/${caseid}_clm_0112.nc
# Create special timeseries for DJF analysis:
ncrcat -O ${drc_out}/${caseid}_clm_12.nc ${drc_out}/${caseid}_clm_01.nc ${drc_out}/${caseid}_clm_02.nc ${drc_out}/DJF_foo.nc

# Remove originals
rm -f ${drc_out}/${caseid}_clm_[0-1][0-9].nc

# STEP 3: Create seasonal means
printf "Creating seasonal means...\n"
ncwa -O -a time -b -w dpm -d time,2,4 ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_MAM.nc
ncwa -O -a time -b -w dpm -d time,5,7 ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_JJA.nc
ncwa -O -a time -b -w dpm -d time,8,10 ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_SON.nc
ncwa -O -a time -b -w dpm -d time,0,2 ${drc_out}/DJF_foo.nc ${drc_out}/${caseid}_clm_DJF.nc
rm -f ${drc_out}/DJF_foo.nc

# STEP 4: Create climate-mean average
# Method: Average annual-mean ensembles
# (alternatively, we could average over the monthly means, applying dpm weights)
printf "Creating climate mean...\n"
ncra -O ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_clm.nc

# STEP 5: Create globally-averaged means:
printf "Creating global means...\n"
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_xyt.nc

ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_ts_DJF.nc ${drc_out}/${caseid}_ts_DJF_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_ts_MAM.nc ${drc_out}/${caseid}_ts_MAM_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_ts_JJA.nc ${drc_out}/${caseid}_ts_JJA_xyt.nc
ncwa -O -a lat,lon -w area ${drc_out}/${caseid}_ts_SON.nc ${drc_out}/${caseid}_ts_SON_xyt.nc

# STEP 6: Create regional averages (OPTIONAL):
printf "Creating regional means...\n"

# North Polar region:
ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP.nc
#ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP.nc
#ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP.nc
#ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP.nc
#ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP.nc
ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP.nc
ncwa -O -a lat,lon -d lat,80.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP.nc

# North Polar (expanded) region2:
ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP2.nc
#ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP2.nc
#ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP2.nc
#ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP2.nc
#ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP2.nc
ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP2.nc
ncwa -O -a lat,lon -d lat,70.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP2.nc

# North Polar (expanded) region3:
ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP3.nc
#ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP3.nc
#ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP3.nc
#ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP3.nc
#ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP3.nc
ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP3.nc
ncwa -O -a lat,lon -d lat,60.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP3.nc

# North Polar (expanded) region4:
ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP4.nc
#ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP4.nc
#ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP4.nc
#ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP4.nc
#ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP4.nc
ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP4.nc
ncwa -O -a lat,lon -d lat,50.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP4.nc

# North Polar (expanded) region4:
ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP5.nc
#ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP5.nc
#ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP5.nc
#ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP5.nc
#ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP5.nc
ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP5.nc
ncwa -O -a lat,lon -d lat,40.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP5.nc

# North Polar (expanded) region5 land-only for comparison to HaQ06:
# landmask variable seems to create a problem with NP6 region
ncwa -O -a lat,lon -d lat,30.,90. -x -v landmask -B 'landfrac=1' -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NP6.nc
#ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NP6.nc
#ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NP6.nc
#ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NP6.nc
#ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NP6.nc
ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NP6.nc
ncwa -O -a lat,lon -d lat,30.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NP6.nc

# Greenland:
ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_GR.nc
#ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_GR.nc
#ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_GR.nc
#ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_GR.nc
#ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_GR.nc
ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_GR.nc
ncwa -O -a lat,lon -d lat,60.,90. -d lon,300.,340. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_GR.nc

# Northern Hemisphere:
ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_NH.nc
#ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_NH.nc
#ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_NH.nc
#ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_NH.nc
#ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_NH.nc
ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_NH.nc
ncwa -O -a lat,lon -d lat,0.,90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_NH.nc

# Tibet:
ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_TP.nc
#ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_TP.nc
#ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_TP.nc
#ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_TP.nc
#ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_TP.nc
ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_TP.nc
ncwa -O -a lat,lon -d lat,30.,40. -d lon,80.,100. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_TP.nc

# South Pole:
ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_SP.nc
#ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_SP.nc
#ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_SP.nc
#ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_SP.nc
#ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_SP.nc
ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_SP.nc
ncwa -O -a lat,lon -d lat,-90.,-90. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_SP.nc

# Antarctica:
ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm.nc ${drc_out}/${caseid}_clm_AN.nc
#ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm_DJF.nc ${drc_out}/${caseid}_clm_DJF_AN.nc
#ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm_MAM.nc ${drc_out}/${caseid}_clm_MAM_AN.nc
#ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm_JJA.nc ${drc_out}/${caseid}_clm_JJA_AN.nc
#ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm_SON.nc ${drc_out}/${caseid}_clm_SON_AN.nc
ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_ts_ANN.nc ${drc_out}/${caseid}_ts_ANN_AN.nc
ncwa -O -a lat,lon -d lat,-90.,-60. -w area ${drc_out}/${caseid}_clm_0112.nc ${drc_out}/${caseid}_clm_0112_AN.nc

# Move /tmp files back into original directories:
mv ${drc_in}/tmp/* ${drc_in}
rm -rf ${drc_in}/tmp
