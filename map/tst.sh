#!/bin/sh

# $Id$

# Purpose: Run Dust boundary dataset generator in NQS environment
# Source: NQS setup follows MATCH and CCM templates

# Usage:
# NCAR AIX: llsubmit tst.sh
# NCAR ALPHA: qsub -q reg -o ~/map/map_tst.txt tst.sh
# NCAR SGI: qsub -q ded_4 -o ~/map/map_tst.txt tst.sh
# UCI Linux: tst.sh >map_tst.txt 2>&1 &
# UCI SGI: qsub -q q4 -o ~/map/map_tst.txt tst.sh
# scp ~/map/tst.sh dataproc.ucar.edu:/fs/cgd/home0/zender/map/tst.sh
# scp ~/map/tst.sh krein.math.uci.edu:map/tst.sh

# AIX LoadLeveler batch system
# @ job_name    = map_tst
# @ output      = $(job_name).txt
# @ error       = $(job_name).txt
# @ node_usage  = not_shared
# @ class       = share
# @ queue

# NQS tokens
# -A nbr        Project number (CMS=03010063, AMWG=93300075)
# -J m          Append NQS output to stdout
# -eo           Glue together stderr and stdout
# -l mpp_p=nbr  Maximum # of CPUs = OMP_NUM_THREADS+4
# -lm sz -lm sz Maximum memory footprint
# -lt tm -lT tm Maximum run time in queue
# -me           Send mail when job ends
# -o fl_out     Output file
# -q queue_name Queue name
# -s shell      Shell to use (default is /bin/sh)

# UCI and NCAR both
#QSUB -me
#QSUB -J m

# UCI only
#QSUB -q q4
##QSUB -q genq1_4
#QSUB -l mpp_p=8
#QSUB -eo -o /home/ess/zender/zender/map/map_tst.txt

# NCAR only
##QSUB -q ded_4
##QSUB -eo -o /fs/cgd/home0/zender/map/map_tst.txt
##QSUB -lt 420000 -lT 420000
##QSUB -lm 650Mb -lM 650Mb

# Main code
set timestamp
set echo
date

# OS-specific
case "${PVM_ARCH}" in 
    AIX* ) 
	data_drc='/ptmp/zender/map'
	export OMP_NUM_THREADS=2
    ;; # endif AIX*
    LINUX* ) 
	data_drc='/usr/tmp/zender/map'
	export OMP_NUM_THREADS=2
    ;; # endif LIN*
    SGI* )
	data_drc='/usr/tmp/zender/map'
	export OMP_NUM_THREADS=4
    ;; # endif SGI
    * )
	data_drc='/usr/tmp/zender/map'
	export OMP_NUM_THREADS=1
    ;; # endif default
esac # endcase ${PVM_ARCH}

# Job-specific
CASE=map_tst
# Run program
mkdir -p ${data_drc}
cd ${data_drc}
PVM_ARCH=`${HOME}/bin/sh/pvmgetarch`
MY_BIN_DIR=${HOME}/bin/${PVM_ARCH}
EXE_CMD="${MY_BIN_DIR}/bds -D 4 -x 90 -y 90 -m 13 -t 19801980 -o ${data_drc}/dst_1x1.nc"

printf "Running job ${CASE}\n"
if [ "${DBG}" = '1' ]; then
    printf "WARNING: Debugging activated in run script\n"
fi # endif DBG
printf "OMP_NUM_THREADS = ${OMP_NUM_THREADS}\n"
printf "Executing command: ${EXE_CMD}\n"
case ${PVM_ARCH} in
    AIX* ) ${EXE_CMD} ; ;; # endif AIX*
    LINUX* ) nice -n 19 ${EXE_CMD} ; ;; # endif LINUX*
    * ) ${EXE_CMD} ; ;; # end default
esac # end case

if [ "${status}" != '0' ]; then
 goto err
fi # endif err

# Successful exit branch
date;pwd;ls -la
exit 0

# Error exit branch
err:
date;pwd;ls -la
exit 1
