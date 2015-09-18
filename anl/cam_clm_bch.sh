# $Id$

# Purpose: Batch driver for processing CAM/CLM/SNICAR simulations

# Conventions:
# yr_* variables may be in yyyy format, but need not be
# yr_*_rth variables have no leading zeros (automatically trimmed for arithmetic)
# yyyy_* variables are in yyyy format

# Testing:
# /bin/cp -p -r ${DATA}/sncpd10 ${DATA}/sncpd10_bck
# /bin/rm -r -f ${DATA}/sncpd10 ${DATA}/sncpd10;/bin/cp -p -r ${DATA}/sncpd10_bck ${DATA}/sncpd10

# Production:
# cd ~/anl;chmod a+x *.sh
# cd ~/anl;./cam_clm_bch.sh
# cd ~/anl;time ./cam_clm_bch.sh > ~/foo.bch 2>&1 &
# cd ~/anl;time ./cam_clm_bch.sh sncpd05 0005 0025 1 1 1 > ~/foo.bch 2>&1 &

# Latest runs:
# cd /data/zender;mkdir sncpd10 sncpd03 snclgm02 snclgm03
# cd /data/zender/sncpd05;rsync *000[5-9]* *001[0-9]* *002[0-5]* sand.ess.uci.edu:/data/zender/sncpd05
# cd /data/zender/sncpd06;rsync *000[5-9]* *001[0-9]* *002[0-5]* sand.ess.uci.edu:/data/zender/sncpd06
# cd /data/zender/snclgm05;rsync *000[5-9]* *001[0-9]* *002[0-5]* sand.ess.uci.edu:/data/zender/snclgm05
# cd /data/zender/snclgm06;rsync *000[5-9]* *001[0-9]* *002[0-5]* sand.ess.uci.edu:/data/zender/snclgm06
# cd ~/anl;time ./cam_clm_bch.sh sncpd05 0005 0025 1 1 1 > ~/foo.bch 2>&1 &
# cd ~/anl;time ./cam_clm_bch.sh sncpd06 0005 0025 1 1 1 > ~/foo.bch 2>&1 &
# cd ~/anl;time ./cam_clm_bch.sh snclgm05 0005 0025 1 1 1 > ~/foo.bch 2>&1 &
# cd ~/anl;time ./cam_clm_bch.sh snclgm06 0005 0025 1 1 1 > ~/foo.bch 2>&1 &

# 0. Command-line options
caseid='sncpd10'
# flg_dst: Similar to flg_sot, but not all options implemented
# 0: Exclude all dust fields
# 1: Include all (FRC+FDB) dust atmosphere and snow fields
flg_dst='1' # [flg] Files contain dust
# flg_sot:
# 0: Exclude all soot fields
# 1: Include all (FRC+FDB) soot atmosphere and snow fields
# 2: Include just atmosphere (FRC+FDB) soot fields
# 3: Include SNW+ATM FDB soot fields
# 4: Include ATM FDB soot fields
# 5: Include SOTDEP fields (passive tracer experiments)
flg_sot='1' # [enm] Files contain soot
flg_snc='1' # [flg] Files contain SNICAR fields
yr_srt='0000' # Start year (NB: Include year 0 for now)
yr_end='0020' # End year

if [ -n "${1}" ]; then
    caseid=${1}
fi # !$1
if [ -n "${2}" ]; then
    yr_srt=${2}
fi # !$2
if [ -n "${3}" ]; then
    yr_end=${3}
fi # !$3
if [ -n "${4}" ]; then
    flg_dst=${4}
fi # !$4
if [ -n "${5}" ]; then
    flg_sot=${5}
fi # !$5
if [ -n "${6}" ]; then
    flg_snc=${6}
fi # !$6

# Print initial state
if [ '1' = '1' ]; then # dbg
    printf "${0}: Initialization State:\n"
    printf "${0}: caseid=${caseid}\n"
    printf "${0}: yr_srt=${yr_srt}\n"
    printf "${0}: yr_end=${yr_end}\n"
    printf "${0}: flg_dst=${flg_dst}\n"
    printf "${0}: flg_sot=${flg_sot}\n"
    printf "${0}: flg_snc=${flg_snc}\n"
fi # !dbg

# 1. Primary variables
flg_CAM='1' # [enm] 1 for CAM, 2 for CLM

# 2. Derived variables

if [ ${flg_CAM} = '1' ]; then

#    swm_ssn_key='zender01'; # Name of persistent SWAMP session
#    python ~/swamp/src/swamp_client.py -d ${swm_ssn_key}
#    python ~/swamp/src/swamp_client.py -n ${swm_ssn_key}
#    python ~/swamp/src/swamp_client.py -k ${swm_ssn_key} ${HOME}/anl/cam_clm_mrg.sh ${caseid} ${yr_srt} ${yr_end} ${flg_dst} ${flg_sot} ${flg_snc}
    ${HOME}/anl/cam_clm_mrg.sh ${caseid} ${yr_srt} ${yr_end} ${flg_dst} ${flg_sot} ${flg_snc}
    ${HOME}/anl/cam_clm_avg.sh ${caseid} ${yr_srt} ${yr_end}
    ${HOME}/anl/cam_clm_drv.sh ${caseid} ${flg_dst} ${flg_sot}

elif [ ${flg_CAM} = '2' ]; then

# csz scripts do not have clm-only counterparts yet
    caseid=clm24
    yr_srt=1979
    yr_end=1999
    
    ${HOME}/anl/cam_clm_mrg.sh ${caseid} ${yr_srt} ${yr_end}
    ${HOME}/anl/cam_clm_avg.sh ${caseid} ${yr_srt} ${yr_end}
    ${HOME}/anl/cam_clm_drv.sh ${caseid} ${flg_dst} ${flg_sot}
fi # !flg_CAM
