#! /bin/sh

# $Id$ -*-shell-script-*-

# Purpose: Run mie optical property generator in batch environment
# Source: NQS setup follows MATCH and CCM templates

# Usage:
# NCAR AIX: llsubmit mie.sh (monitor with llq, kill with llcancel)
# NCAR Cray: qsub mie.sh 
# NCAR SGI: npri -h 250 mie.sh >mie.txt 2>&1 &
# UCI Linux: mie.sh >mie.txt 2>&1 &
# UCI SGI: qsub -q q4 -o ~/mie/mie.txt mie.sh
# scp ~/mie/mie.sh ~/mie/mie.pl esmf.ess.uci.edu:mie
# scp ~/mie/mie.sh ~/mie/mie.pl goldhill.cgd.ucar.edu:/fs/cgd/home0/zender/mie

# AIX LoadLeveler batch system
# class: Queue name. "llclass" lists available queues.
# node: Number of nodes. 
# tasks_per_node: MPI processes per node, set to 1 for Hybrid OpenMP/MPI codes 
# output: Script output for STDOUT
# error: Script output for STDERR
# job_type: parallel declares that multiple nodes will be used
# network.MPI: Network connection type between nodes. Leave this be.
# node_usage: not_shared acquires dedicated access to nodes
# queue: Tells Loadleveler to submit job

#@ job_name       = mie20
##@ class          = com_rg8
##@ class          = com_rg32
#@ class          = com_node03
#@ node           = 1
#@ tasks_per_node = 1
#@ output         = $(job_name).txt
#@ error          = $(job_name).txt
#@ job_type       = parallel
#@ network.MPI    = csss,shared,us
#@ node_usage     = not_shared
##@ account_no     = 36271012
##@ wall_clock_limit = 3800
#@ queue

# NQS tokens
# -A nbr        Project number (CMS=03010063, AMWG=93300075, ESMF=36271012)
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
#QSUB -eo -o /home/ess/zender/zender/mie/mie.txt

# NCAR only
##QSUB -q ded_4
##QSUB -eo -o /fs/cgd/home0/zender/mie/mie.txt
##QSUB -lt 420000 -lT 420000
##QSUB -lm 650Mb -lM 650Mb

# Main code
set timestamp
set echo
date
unalias -a # Turn off aliases

# OS-generic
export NTASKS=1 # MPI
PVM_ARCH=`~zender/bin/sh/pvmgetarch`
HOST=`hostname`
# OS-specific
case "${PVM_ARCH}" in 
    AIX* ) 
	export DATA=/ptmp/${USER}
	case "${HOST}" in 
	    esmf* ) # UCI
		case "${HOST}" in 
		    esmf0[1-7]* ) # UCI
		    export NTHREADS=8 # OpenMP
		    ;; # endif UCI
		    esmf0[8]* ) # UCI
		    export NTHREADS=32 # OpenMP
		    ;; # endif UCI
		esac # endcase ${HOST}
		;; # endif UCI
	    b[bfs]*en ) # babyblue, blackforest, bluesky
	    export NTHREADS=8 # OpenMP
	    ;; # endif NCAR
	esac # endcase ${HOST}
	;; # endif AIX*
    LINUX* ) 
	export DATA=/data/${USER}
	export NTHREADS=2
# Attempt to get sufficient stack memory
	ulimit -s unlimited
    ;; # endif LIN*
    SGI* )
	export DATA=/data/${USER}
	export NTHREADS=4
    ;; # endif SGI
    * )
	echo "ERROR: $PVM_ARCH is unsupported operating system"
	exit 1
    ;; # endif default
esac # endcase ${PVM_ARCH}

# Job-specific
PRG_NM='mie' # [sng] Program name, semantic identifier
FL_NM_SH='mie.sh' # [sng] Shell batch script name
FL_NM_PL='mie.pl' # [sng] Perl batch script name
CASEID='mie21' # [sng] Case ID
# NB: Whitespace in this string currently breaks mie.pl restart capability
XPT_DSC='Size_parameter_experiment_with_liquid_cloud_droplets' # [sng] Experiment description
# Derive data paths
export DATA_RT=${DATA}/aca
export DATA_OUT=${DATA}/${PRG_NM}/${CASEID}
# Run program
if [ ! -d ${DATA_OUT} ]; then
    mkdir -p ${DATA_OUT}
fi # endif
cd ${DATA_OUT}
export LID="`date +%Y%m%d-%H%M%S`"
EXE=${DATA_OUT}/${PRG_NM}
FL_PL=${DATA_OUT}/${FL_NM_PL}
/bin/cp -f -p ~/bin/${PVM_ARCH}/${PRG_NM} ${DATA_OUT} || exit 1
/bin/cp -f -p ~/mie/${FL_NM_PL} ${DATA_OUT} || exit 1
/bin/cp -f -p ~/mie/${FL_NM_SH} ${DATA_OUT} || exit 1
# Copy from Production lines in mie.pl:
# Sulfate
#CMD_LN="${FL_PL} --dbg=1 --xpt_dsc=${XPT_DSC} --fl_lbl=${CASEID} --fl_nbr=480 --cmp_prt=sulfate --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=1 --bnd_nbr_dfl=100000 --sz_nbr_dfl=30 --sz_mnm_dfl=0.001 --sz_mxm_dfl=1.0 --rds_nma_dfl=0.08 --gsd_anl_dfl=1.2 --drc_dat=${DATA_RT} --drc_out=${DATA_OUT}"
# Sample all size parameters
CMD_LN="${FL_PL} --dbg=1 --thr_nbr=1 --idx_rfr_prt_dfl='1.33+1.0e-6i' --dmn_nbr_max=1 --xpt_dsc='${XPT_DSC}' --fl_lbl=${CASEID} --fl_nbr=480 --cmp_prt=h2o_lqd --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=100000 --bnd_nbr_dfl=1 --sz_nbr_dfl=1 --sz_mnm_dfl=9.999 --sz_mxm_dfl=10.001 --rds_swa_dfl=10.0 --drc_dat=${DATA_RT} --drc_out=${DATA_OUT}"
# Liquid Cloud single thread
#CMD_LN="${FL_PL} --dbg=1 --thr_nbr=1 --dmn_nbr_max=1 --xpt_dsc='${XPT_DSC}' --fl_lbl=${CASEID} --fl_nbr=480 --cmp_prt=h2o_lqd --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=100000 --bnd_nbr_dfl=1 --sz_nbr_dfl=30 --sz_mnm_dfl=0.1 --sz_mxm_dfl=30.0 --rds_swa_dfl=10.0 --gsd_anl_dfl=1.6 --drc_dat=${DATA_RT} --drc_out=${DATA_OUT}"
# Liquid Cloud
#CMD_LN="${FL_PL} --dbg=1 --xpt_dsc='${XPT_DSC}' --fl_lbl=${CASEID} --fl_nbr=480 --cmp_prt=h2o_lqd --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=1 --bnd_nbr_dfl=1000000 --sz_nbr_dfl=30 --sz_mnm_dfl=0.1 --sz_mxm_dfl=30.0 --rds_swa_dfl=10.0 --gsd_anl_dfl=1.6 --drc_dat=${DATA_RT} --drc_out=${DATA_OUT}"
FL_STDOUT="${PRG_NM}.log.${LID}"

echo "Timestamp ${LID}"
echo "Batch shell script $0 running CASEID = ${CASEID} on machine ${HOST}"
echo "Invoking executable with ${CMD_LN}"
echo "STDOUT/STDERR re-directed to file:"
echo "/bin/more ${DATA_OUT}/${FL_STDOUT}"
case "${PVM_ARCH}" in 
    AIX* ) 
# Set POE environment for interactive jobs  
# LoadLeveler batch jobs ignore these settings
# MP_NODES is node number
# XLSMPOPTS thread stack size
	export MP_EUILIB='us'
	export MP_NODES="${NTASKS}"
	export MP_TASKS_PER_NODE='1'
	export MP_RMPOOL='1'
	export XLSMPOPTS='stack=86000000'
	if [ ${NTASKS} -gt 1 ]; then
	    poe ${CMD_LN} > ${FL_STDOUT} 2>&1
	else
	    env OMP_NUM_THREADS="${NTHREADS}" PATH=${PATH}\:/usr/local/bin\:${DATA_OUT} ${CMD_LN} > ${FL_STDOUT} 2>&1
	fi # end else OpenMP
	;; # endif AIX*
    LINUX* ) 
	if [ ${NTASKS} -gt 1 ]; then
	    mpirun -np ${NTASKS} ${CMD_LN} > ${FL_STDOUT} 2>&1
	else     
	    env OMP_NUM_THREADS="${NTHREADS}" MPSTKZ="128M" ${CMD_LN} > ${FL_STDOUT} 2>&1
	fi # end else OpenMP
	;; # endif LIN*
    SGI* )
	export TRAP_FPE='UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE; INVALID=ABORT,TRACE'
	export OMP_DYNAMIC='FALSE'
	export _DSM_PLACEMENT='ROUND_ROBIN'
	export _DSM_WAIT='SPIN'
	export MPC_GANG='OFF'
# MP_SLAVE_STACKSIZE sets size of thread stack
	export MP_SLAVE_STACKSIZE='40000000'
# Run pure SPMD or pure OpenMP 
	if [ ${NTASKS} -gt 1 ]; then
	    mpirun -np ${NTASKS} ${CMD_LN} > ${FL_STDOUT} 2>&1
	else
	    env MP_SET_NUMTHREADS="${NTHREADS}" ${CMD_LN} > ${FL_STDOUT} 2>&1
	fi # end else OpenMP
	;; # endif SGI
    * )
	echo "ERROR: ${PVM_ARCH} is unsupported operating system"
	exit 1
	;; # endif default
esac # endcase ${PVM_ARCH}

exit_status=$?
if [ "${exit_status}" ]; then
# Successful exit branch
    echo "SUCCESS execution of ${FL_NM_PL} completed successfully"
    date;pwd;ls -la
    echo "$0 batch script completed successfully at `date +%Y%m%d-%H%M%S`"
    exit 0
else
# Error exit branch
    echo "ERROR execution of ${FL_NM_PL} failed"
    date;pwd;ls -la
    echo "$0 batch script error exit at `date +%Y%m%d-%H%M%S`"
    exit 1
fi # endif err

