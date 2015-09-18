#!/bin/sh

# Purpose: Climatology script for ACME CAM/CLM output

# Distribution:
# scp ~/anl/acme_clm.sh rhea.ccs.ornl.gov:anl

# Usage:
# chmod a+x ~/anl/acme_clm.sh
# cd ~/anl;./acme_clm.sh $CASEID $yyyy_srt $yyyy_end
# cd ~/anl;./acme_clm.sh ne30_gx1.B1850c5d 0005 0007

# Production:
# cd ~/anl;acme_clm.sh
# cd ~/anl;time ./acme_clm.sh > foo.avg 2>&1 &
# cd ~/anl;time ./acme_clm.sh ne30_gx1.B1850c5d 0005 0007
 
# set -o xtrace # dbg

# ACME Climatology Guidelines:
# https://acme-climate.atlassian.net/wiki/display/ATM/Climo+Files+-+v0.3+AMIP+runs
# Produces:
# ensemble monthly means 
# ensemble seasonal means
# ensemble annual mean

# User sets following environmental variables:
# caseid:  Simulation name (filenames must start with ${caseid})
# drc_out: Output directory for processed, reduced data
#          Script does not modify original data
#          User needs write permission to ${DATA}/${caseid}
# yyyy_srt: Year of first December to analyze
#          Convention is to include Dec of yyyy_srt in seasonal 
#          and monthly analysis, and exclude Dec of yyyy_end
#          All other months in yyyy_srt are excluded
#          Annual means and climate means exclude all months of yyyy_srt,
#          and include all months of yyyy_end.  
#          Seasonal averages are not 100% consistent with annual averages
#          To have continuity in the DJF analysis, there is no alternative
# yyyy_end: Last year of analysis
#          Dec of yyyy_end is excluded from the seasonal and monthly analysis
#          For a proper annual analysis, yyyy_end must be a complete year of output
# mdl:     This tag determines whether CLM or CAM or (other?) output is analyzed
#          It refers to the character sequence in the history files:
#          caseid.mdl.h0.YYYY-MM.nc (valid examples are 'clm2' and 'cam2')        
#          NB: Currently 'h0' is fixed in the file name, but an environmental 
#          variable could easily replace this to analyze other history tapes.

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
caseid='ne30_gx1.B1850c5d'
yyyy_srt='0005' # [yyyy] Start year (NB: fine to include year 0000)
yyyy_end='0007' # [yyyy] End year

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
mdl='cam'

# 2. Derived variables
drc_in="${MY_PROJWORK}/ne30/raw" # [sng] Directory for input files
drc_out="${MY_PROJWORK}/ne30/clm" # [sng] Directory for output files

# 3. Doubly-derived variables
wrk_dir=${drc_out}

# Determine century and first full year
trim_leading_zeros ${yyyy_srt}
yr_srt_rth=${sng_trm}
let yr_srtp1=${yr_srt_rth}+1
trim_leading_zeros ${yyyy_end}
yr_end_rth=${sng_trm}
let yr_endm1=${yr_end_rth}-1

# Print initial state
if [ '0' = '1' ]; then # dbg
    printf "${0}: yyyy_srt = ${yyyy_srt}\n"
    printf "${0}: yr_srt_rth = ${yr_srt_rth}\n"
    printf "${0}: yr_srtp1 = ${yr_srtp1}\n"
    printf "${0}: yr_endm1 = ${yr_endm1}\n"
fi # !dbg

# Create output directory, go to working directory
mkdir -p ${drc_out}
cd ${wrk_dir}

# Step 0: Implement variable subsetting (optional yet fast so helps with debugging)
var_opt='-v FSNT,SOLIN,AODVIS'

# Step 1: Create monthly means
for MM in {01..12}; do
    printf "Creating climate-mean month ${MM} ...\n"
    yr_fl=''
    for yr in `seq $yr_srtp1 $yyyy_end`; do
	YYYY=`printf "%04d" $yr`
	yr_fl="${yr_fl} ${caseid}.${mdl}.h0.${YYYY}-${MM}.nc"
    done 
    cmd="ncra -O ${var_opt} -p ${drc_in} ${yr_fl} ${drc_out}/${caseid}_${MM}_climo.nc"
#    echo ${cmd}
    eval ${cmd}
done

# Step 2: Create special December monthly mean for continual DJF seasonal means
MM='12'
printf "Creating monthly mean for December for seasonally continuous means...\n"
yr_fl=''
for yr in `seq $yr_srt_rth $yr_endm1`; do
    YYYY=`printf "%04d" $yr`
    yr_fl="${yr_fl} ${caseid}.${mdl}.h0.${YYYY}-${MM}.nc"
done 
cmd="ncra -O ${var_opt} -p ${drc_in} ${yr_fl} ${drc_out}/${caseid}_${MM}_ym1_climo.nc"
#echo ${cmd}
eval ${cmd}

# Step 3: Create seasonal means, where default DJF label incorporates seasonally continual months
printf "Creating seasonal means...\n"
ncra -O -w 31,30,31 ${drc_out}/${caseid}_0[3-5]_climo.nc         ${drc_out}/${caseid}_MAM_climo.nc
ncra -O -w 30,31,31 ${drc_out}/${caseid}_0[6-8]_climo.nc         ${drc_out}/${caseid}_JJA_climo.nc
ncra -O -w 30,31,30 ${drc_out}/${caseid}_{09,10,11}_climo.nc     ${drc_out}/${caseid}_SON_climo.nc
ncra -O -w 31,31,28 ${drc_out}/${caseid}_{12,01,02}_climo.nc     ${drc_out}/${caseid}_DJF_ann_climo.nc
ncra -O -w 31,31,28 ${drc_out}/${caseid}_{12_ym1,01,02}_climo.nc ${drc_out}/${caseid}_DJF_climo.nc

# Step 4: Create climate-mean where DJF input does not cross year boundaries
printf "Creating climate mean...\n"
ncra -O ${drc_out}/${caseid}_{MAM,JJA,SON,DJF_ann}_climo.nc ${drc_out}/${caseid}_ANN_climo.nc
