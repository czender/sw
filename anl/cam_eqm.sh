#!/bin/sh

# $Id$ 

# Purpose: Diagnose CAM simulations for equilibrium conditions
# Assumes monthly mean data are in ${DATA}/${caseid}
# Input files are not altered in any way

# Conventions:
# yr_* variables may be in yyyy format, but need not be
# yr_*_rth variables have no leading zeros (automatically trimmed for arithmetic)
# yyyy_* variables are in yyyy format

# Usage:
# chmod a+x ~/anl/cam_eqm.sh
# cd ~/anl;./cam_eqm.sh
# cd ~/anl;./cam_eqm.sh sncpd02 0000 0020
# cd ~/anl;./cam_eqm.sh snclgm03 0000 0030
# cd ~/anl;./cam_eqm.sh bc_seta 0001 0011
# cd ~/anl;./cam_eqm.sh bc_setb 0001 0008
# cd ~/anl;./cam_eqm.sh bc_setb_aus 0001 0015
# cd ~/ess_gcm;./cam_eqm.sh essgcm12 0001 0010

# Production:
# cd ~/anl;./cam_eqm.sh
# cd ~/anl;./cam_eqm.sh > ~/foo.eqm 2>&1
# cd ~/anl;./cam_eqm.sh ${caseid} ${yr_srt} ${yr_end}

function trim_leading_zeros {
# Purpose: Trim leading zeros from strings that represent integers
# Bash treats zero-padded integers as octal
# This is surprisingly hard to work around
# My approach is to remove leading zeros prior to arithmetic
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

#set -o xtrace # dbg

# 0. Command line options
caseid='sncpd03' # [sng] Case ID
yr_srt='0' # [nbr] Start year
yr_end='20' # [nbr] End year

if [ -n "${1}" ]; then
    caseid=${1}
fi # !$1
if [ -n "${2}" ]; then
    yr_srt=${2}
fi # !$2
if [ -n "${3}" ]; then
    yr_end=${3}
fi # !$3

# 1. Primary variables
drc_in="${DATA}/${caseid}" # [sng] Directory for input files
drc_out="${DATA}/anl" # [sng] Directory for output files

mkdir -p ${drc_out}
cd ${drc_out}

#fld_lst='FSNT,FLNT,SNODSTCI,SNODSTMCI,SNOBCCI,SNOBCMCI,TREFHT,gw' # SNICAR runs
fld_lst='FSNT,FLNT,LANDFRAC,OCNFRAC,PRECC,PRECL,TREFHT,gw' # Standard CAM runs
#fld_lst='FSNT,FLNT,LANDFRAC,OCNFRAC,PRECC,PRECL,PRECSL,PRECSC,TREFHT,BCDEPDRY,BCDEPWET,BCPHISF,BCPHOSF,gw' # BC runs
for yr in `seq $yr_srt $yr_end` ; do
    yyyy=`printf "%04d" ${yr}`
    printf "Averaging ${caseid} data for year ${yyyy}...\n"
# NB: Weighting each month by days_per_month would be more accurate
    if ! ncra -O -v ${fld_lst} -n 12,2,1 -p ${drc_in} ${caseid}.cam2.h0.${yyyy}-01.nc ${drc_out}/${caseid}_${yyyy}.nc ; then
	echo "${0}: ncra-averaging failed for year ${yr}, aborting..." && exit 1
    fi # endif
    if ! ncks -A -C -v area -p ${drc_in} ${caseid}.clm2.h0.${yyyy}-01.nc ${drc_out}/${caseid}_${yyyy}.nc ; then
	echo "${0}: Adding area field failed for year ${yr}, aborting..." && exit 1
    fi # endif
    if ! ncap2 -O -s 'FTNT=FSNT-FLNT;PRECT=PRECC+PRECL' ${drc_out}/${caseid}_${yyyy}.nc ${drc_out}/${caseid}_${yyyy}.nc ; then
#    if ! ncap2 -O -s 'FTNT=FSNT-FLNT;PRECT=PRECC+PRECL;PRECS=PRECSL+PRECSC;SOTDEP=-(BCDEPDRY+BCDEPWET)' ${drc_out}/${caseid}_${yyyy}.nc ${drc_out}/${caseid}_${yyyy}.nc ; then
    echo "${0}: ncap2 processing failed, aborting..." && exit 1
    fi # endif
done # end loop over yr

# NB: Bash treats zero-padded integers as octal (UGGH)
# Remove leading zeros prior to arithmetic:
trim_leading_zeros ${yr_srt}
yr_srt_rth=${sng_trm}
trim_leading_zeros ${yr_end}
yr_end_rth=${sng_trm}
let yr_nbr="${yr_end_rth}-${yr_srt_rth}+1"

# Convert yr_srt_rth, yr_end_rth to yyyy-format for output filenames
yyyy_srt=`printf "%04d" ${yr_srt_rth}`
yyyy_end=`printf "%04d" ${yr_end_rth}`
if ! ncrcat -O -n ${yr_nbr},4,1 ${drc_out}/${caseid}_${yyyy_srt}.nc ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}.nc ; then
    echo "${0}: ncrcat concatenation failed, aborting..." && exit 1
fi # endif
if ! ncra -O ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}.nc ${drc_out}/${caseid}_clm.nc ; then
    echo "${0}: ncra temporal averaging failed, aborting..." && exit 1
fi # endif
if ! ncwa -O -a lat,lon -w gw ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}.nc ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}_txy.nc ; then
    echo "${0}: ncwa spatial averaging failed, aborting..." && exit 1
fi # endif
if ! ncwa -O -a time ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}_txy.nc ${drc_out}/${caseid}_clm_txy.nc ; then
    echo "${0}: ncwa time averaging failed, aborting..." && exit 1
fi # endif

printf "Annual timeseries of global-mean TREFHT, F^N_R(TOA) complete.\n"
printf "Computed ${yr_nbr} annual means from year ${yr_srt_rth} to year ${yr_end_rth}.\n"
printf "View timeseries and spatial data, respectively, with:\n"
printf "ncview ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}_txy.nc &\n"
printf "ncview ${drc_out}/${caseid}_${yyyy_srt}_${yyyy_end}.nc &\n"
