#! /bin/csh
# -f = fast start, ignore .cshrc and .login

# $Id$

# Purpose: Run MATCH-Dust model on all machines
# Source: http://www.cgd.ucar.edu/csm/working_groups/Atmosphere/scriptccm3.html

# Usage:
# NCAR AIX: llsubmit match.sh (monitor with llq, kill with llcancel)
# NCAR SGI: npri -h 250 match.sh >dstmch87.txt 2>&1 &
# NCAR Cray: qsub match.sh 
# UCI Linux: match.sh >dstmch87.txt 2>&1 &
# UCI SGI: qsub match.sh
# scp ~/match/match.sh dataproc.ucar.edu:/fs/cgd/home0/zender/match/match.sh
# scp ~/match/match.sh krein.math.uci.edu:match/match.sh

# Debugging: 
# Use totalview for interactive debugger
# Use debugview for diagnosing batch failures
# totalview -c core match -L -include_files
# totalview -c core match

# 20000622: match-4_0_beta2-dst-1_1_1 T62 pcnst=5 PRC=S requires 175 MB prc-1 on Linux
# 20000622: match-3_3_23-dst-0_9_1 T62 pcnst=5 requires 62 hrs wallclock per year with no debugging (except qneg,fixmas) on krein
# 20000620: match_dst T62 pcnst=5 requires 81 hrs wallclock per year with full debugging (qneg,DST_DBG,DST_MSS_BDG,fixmas) on krein
# 19990601: match_dst T42 pcnst=5 requires 216 MB, 320+ MB with -g
# 19991210: match_dst T62 pcnst=5 requires 284 MB prc-1 on SGI
# AIX LoadLeveler batch system
#@ job_name    = dstmch87
#@ class       = share
#@ node_usage  = not_shared
#@ node        = 1
#@ job_type    = parallel
#@ network.MPI = css0,shared,us
#@ total_tasks = 1
#@ error       = $(job_name).txt
#@ output      = $(job_name).txt
#@ account_no  = 93300075
#@ queue

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
#QSUB -J m
#QSUB -me
#QSUB -mu zender@uci.edu
#QSUB -s /bin/csh

# UCI only
#QSUB -q q4 -l mpp_p=8
##QSUB -q genq1_4 -l mpp_p=8
##QSUB -q pub_8 -l mpp_p=12
##QSUB -q pub_16 -l mpp_p=20
#QSUB -eo -o /home/ess/zender/zender/match/dstmch87.txt

# NCAR only
##QSUB -q fesb -l mpp_p=8
###QSUB -q share_16 -l mpp_p=12
###QSUB -q ded_16 -l mpp_p=20
###QSUB -q ded_16
##QSUB -eo -o /fs/cgd/home0/zender/match/dstmch87.txt
###QSUB -lt 10000 -lT 10000
####QSUB -lt 70000 -lT 70000
####QSUB -lt 420000 -lT 420000
###QSUB -lm 280Mb -lM 280Mb

# Main code
set timestamp
set echo
date

# Job-specific
setenv CASE dstmch87
if (${HOST} =~ *.uci.edu) then
    set wrk_dir = /data/zender/match/${CASE}
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute* ) then
    set wrk_dir = /ptmp/${USER}/match/${CASE}
else
    printf "ERROR: Unknown HOST in run script\n"
    exit
endif # end domain-specific switches
setenv EXE ${wrk_dir}/match
if ( ! -d ${wrk_dir} ) mkdir -p ${wrk_dir}
cd ${wrk_dir} || exit 1

# Domain-specific
if (${HOST} =~ *.uci.edu) then
    set MY_BIN_DIR = ${HOME}/match/bin/${PVM_ARCH} # Executable
    set MY_DAT_DIR = /data/zender/data # Boundary data
    set MY_DGN_DIR = ${DATA}/dgn # Diagnostic output data
    set MY_RUN_DIR = ${HOME}/match # Run script
    set MY_SRC_DIR = ${HOME}/match_dst # Source code
    /bin/ln -s -f ${HOME}/bin/sh/msread msread || exit 1
    /bin/ln -s -f ${HOME}/bin/sh/mswrite mswrite || exit 1
    /bin/ln -s -f ${HOME}/bin/sh/rdsp rdsp || exit 1
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute* ) then
    set MY_BIN_DIR = /fs/cgd/home0/${USER}/match/bin/${PVM_ARCH} # Executable
    set MY_DAT_DIR = /fs/cgd/data0/zender/data # Boundary data
    set MY_DGN_DIR = /fs/cgd/data0/${USER}/dgn # Diagnostic output data
    set MY_RUN_DIR = /fs/cgd/home0/${USER}/match # Run script
    set MY_SRC_DIR = /fs/cgd/data0/${USER}/match_dst # Source code
    /bin/ln -s -f ${MY_SRC_DIR}/utils/rdsp rdsp || exit 1
else
    printf "ERROR: Unknown HOST in run script\n"
    exit
endif # end domain-specific switches

set DBG = 0
limit stacksize unlimited
limit datasize  unlimited
limit memoryuse unlimited
# OS-specific
if ("${PVM_ARCH}" == "AIX") then
    setenv MP_EUILIB us
    setenv MP_NODES 1
    setenv MP_RMPOOL 1
    setenv MP_TASKS_PER_NODE 1
    setenv OMP_NUM_THREADS 4
    setenv XLSMPOPTS "stack=255000000:parthds=${OMP_NUM_THREADS}"
else if ("${PVM_ARCH}" == "CRAY") then # endif AIX
    setenv NCPUS 8
else if ("${PVM_ARCH}" == "LINUX") then # endif CRAY
    if (${HOST} =~ dust* || ${HOST} =~ biogenic* ) then
	setenv OMP_NUM_THREADS 2
    else
	setenv OMP_NUM_THREADS 1
    endif # end host-specific switches
else if ("${PVM_ARCH}" == "SGIMP64") then # endif LINUX
    setenv _DSM_PLACEMENT ROUND_ROBIN # ROUND_ROBIN ensures efficient memory allocation
    #setenv _DSM_VERBOSE # Print diagnostic info about system-level software
    setenv MPC_GANG OFF # Gang scheduling does not work well so always turn off
    setenv MP_SLAVE_STACKSIZE 32000000 # ???
    setenv OMP_NUM_THREADS 4
    if (${DBG} == '1') then
	setenv TRAP_FPE "UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE" # Ensure abort IEEE errors, requires using LDFLAGS+= -lfpe on SGI
	setenv F90_BOUNDS_CHECK_ABORT YES
    endif # endif DBG
    # Environment-specific
    if ($?ENVIRONMENT) then
	if ( "${ENVIRONMENT}" == "BATCH" ) then
	    setenv _DSM_WAIT SPIN # Set to SPIN iff have dedicated access to processors
	    setenv OMP_DYNAMIC FALSE # TRUE = give up processors when machine is busy, use FALSE in dedicated queues
	endif # endif NQS batch job
    endif # endif ENVIRONMENT exists
else # endif SGIMP64
    printf "ERROR: Unknown PVM_ARCH in run script\n"
endif # endif PVM_ARCH
if (${DBG} == '1') then
    setenv OMP_NUM_THREADS 1
endif # endif DBG

# Cleanup previous run
/bin/rm -f core
/bin/rm -f perf.data
setenv RESTART False
if ("${PVM_ARCH}" == "AIX") then
    if ( -f ${MY_RUN_DIR}/${CASE}.rsb ) then
        setenv RESTART True
    endif
endif
if ("$RESTART" != "True") then
    /bin/rm -f r0???
    /bin/rm -f h0???
    /bin/rm -f h0???.nc
    /bin/rm -f namelist
    /bin/rm -f rstrt
    /bin/rm -f match
    /bin/rm -f r??-??*
    /bin/rm -f ??-??
    /bin/rm -f ???r_????????_?????*
    /bin/rm -f dst_mss_bdg.nc
# Remove symbolic links in case source directory changed names
    /bin/rm -f T62n.T42n.avg.nc
    /bin/rm -f dst_T42.nc
    /bin/rm -f dst_T62.nc
    /bin/rm -f grb2d.list
    /bin/rm -f grbsanl.list
    /bin/rm -f grib.table2.rean
# Assemble new run
    /bin/cp -f ${MY_BIN_DIR}/match ${EXE} || exit 1
endif # endif !RESTART

# Copy boundary files to local run directory in case boundary files change during run
# Select boundary data with appropriate resolution and erodibility factor
#/bin/cp -f ${MY_DAT_DIR}/dst_T42.nc dst_bnd.nc || exit 1
#/bin/cp -f ${MY_DAT_DIR}/dst_T62.nc dst_bnd.nc || exit 1
#/bin/cp -f ${MY_DAT_DIR}/dst_T62_notaper.nc dst_bnd.nc || exit 1
/bin/cp -f ${MY_DAT_DIR}/dst_T62_geo.nc dst_bnd.nc || exit 1

if (${HOST} =~ *.uci.edu) then
    set rpthdyn = '/DSS' # NCEP dynamics input remote directory
    set lpthdyn = '/data/zender/DSS' # NCEP dynamics input local directory
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute* ) then
    set rpthdyn = '/DSS' # NCEP dynamics input remote directory
    set lpthdyn = "/ptmp/${USER}/NCEP" # NCEP dynamics input local directory
    #set rpthdyn  = "/ZENDER/csm/dst01/ccm3/hist/h0001" # CCM dynamics input remote directory
    #set lpthdyn = "/ptmp/${USER}/CCM" # CCM dynamics input local directory
endif # end domain-specific switches
if ( ! -d ${lpthdyn} ) mkdir -p ${lpthdyn}
# NCEP dynamics need local GRIB lookup tables
set grb_dir = ${MY_SRC_DIR}/readers/ncep/data
if ( ! -e grb2d.list ) /bin/ln -s ${grb_dir}/grb2d.list grb2d.list
if ( ! -e grbsanl.list ) /bin/ln -s ${grb_dir}/grbsanl.list grbsanl.list
if ( ! -e grib.table2.rean ) /bin/ln -s ${grb_dir}/grib.table2.rean grib.table2.rean
# NCEP regrids to T42 using these map overlap weights
if ( ! -e T62n.T42n.avg.nc ) /bin/ln -s ${grb_dir}/T62n.T42n.avg.nc T62n.T42n.avg.nc

# OS-specific run parameters
if ("${PVM_ARCH}" == "CRAY") then
    set ndens = 4 # ndens > 1 only works in double precision (PRC=D)
else if ("${PVM_ARCH}" == "AIX" || "${PVM_ARCH}" == "LINUX" || "${PVM_ARCH}" == "SGIMP64") then
# UCI SGI (krein) does not pack tapes correctly
    set ndens = 1 # ndens > 1 only works in double precision (PRC=D)
else
    printf "ERROR: Unknown PVM_ARCH in run script\n"
endif # endif PVM_ARCH

# Testing:
# nestep   = 4
# nhtfrq   = 1
# Production:
# nestep   = -396
# nhtfrq   = -24
# Warning: MATCH does not have nelapse, remember to change nestep on restarts!

# Initialize values, make sure 
set rstflg = 0 
set nestep_max = 4 # Total number of steps in simulation
set nestep_ncr = 2  # Number of steps to compute in each queue residence
set nestep = ${nestep_ncr} # Number of steps since beginning of run
# nhtfrq and mfilt must be commensurate with nestep_ncr
set nhtfrq = 2 # Averaging (write) frequency
set mfilt  = 2 # Averaging periods per tape
# Update values for re-submitted runs
if ("${PVM_ARCH}" == "AIX") then
    if ( -f ${MY_RUN_DIR}/${CASE}.rsb ) then
	set rstflg = 1
	set nestep = `head -n 1 ${MY_RUN_DIR}/${CASE}.rsb`
	printf "Ingested rstflg = ${rstflg} and nestep = ${nestep} from existing ${CASE}.rsb file\n"
    endif # endif ${CASE}.rsb
endif # endif AIX

cat >! namelist << EOF || exit 1
 &NLIST
 tracnam  = 'DSTQ01','DSTQ02','DSTQ03','DSTQ04'
 title    = '${CASE} MATCH 4.0beta2 NCEP T62 dst-1.1.5: fdg=0.7e-3 vai=0.3 area nocoast nopole'
 fixmas   = .true.
 qrelax   = 0.0
 delt     = 2400
 wpasswd  = 'passwd'
 rstflg   = ${rstflg}
 rmout    = .true.
 rmrst    = .true.
 nestep   = ${nestep}
 nhtfrq   = ${nhtfrq}
 mfilt    = ${mfilt}
 ndens    = ${ndens}
 irt      = 1000
 rstrrt   = 365
 icdate   = 20000701
 icdatesec    = 0
 rpthdyn  = '${rpthdyn}'
 lpthdyn  = '${lpthdyn}'
 dynoffset= 0
 outinst  = 'ORO','PHIS'
 outtimav = 'U','V','PS','T','Q','RELHUM','CLOUD'
	    ,'PRECC','PRECL','ZPRECC','CWAT','PRECT','PRECT2'
	    ,'BSN_FCT','FRC_WET','TPT_GND','VWC_SFC','GWC_SFC','LND_MBL','VAI_DST'
	    ,'WND_FRC','WND_FRCS','WND_FRCT','WND_RFR','WND_RFRT'
	    ,'DSTQ','DSTQ01','DSTQ02','DSTQ03','DSTQ04'
	    ,'DSTMPC','DSTMPC01','DSTMPC02','DSTMPC03','DSTMPC04'
	    ,'DSTODXC','DSTODX01','DSTODX02','DSTODX03','DSTODX04'
	    ,'DSTSFMBL','DSTSFM01','DSTSFM02','DSTSFM03','DSTSFM04'
	    ,'DSTSFPCP','DSTSFP01','DSTSFP02','DSTSFP03','DSTSFP04'
	    ,'DSTSFTRB','DSTSFT01','DSTSFT02','DSTSFT03','DSTSFT04'
	    ,'DSTSFGRV','DSTSFG01','DSTSFG02','DSTSFG03','DSTSFG04'
            ,'DSTSSPCP','DSTSSEVP','DSTSSDRY'
 /
EOF

pwd;ls -la

if ("${PVM_ARCH}" == "CRAY") then
    ja
endif
printf "Running job ${CASE}\n"
if (${DBG} == '1') then
    printf "WARNING: Debugging activated in run script\n"
endif # endif DBG
printf "OMP_NUM_THREADS = ${OMP_NUM_THREADS}\n"
limit -h
if ("${PVM_ARCH}" == "AIX") then
#    poe ${EXE} < namelist -procs 1 -pgmmodel spmd
    ${EXE} < namelist
else if ("${PVM_ARCH}" == "LINUX") then
#    nice --adjustment=19 ${EXE} < namelist
    nice -19 ${EXE} < namelist
else
    ${EXE} < namelist
endif
if ( $status != 0 ) goto err

# Resubmit until finished
if ("${PVM_ARCH}" == "AIX") then
    cd ${MY_RUN_DIR} || exit 1
    # Set rstflg and increment nestep
    @ nestep += ${nestep_ncr}
    if ( ${nestep} < 0 ) then # nestep in days
        if ( ${nestep} >= ${nestep_max} ) then
	    echo ${nestep} >! ${CASE}.rsb
	    /bin/mv -f ${CASE}.txt ${CASE}.txt.${nestep}
	    llsubmit ${0}
        endif
    endif
    if ( ${nestep} > 0 ) then # nestep in timesteps
        if ( ${nestep} <= ${nestep_max} ) then
	    echo ${nestep} >! ${CASE}.rsb
	    /bin/mv -f ${CASE}.txt ${CASE}.txt.${nestep}
	    llsubmit ${0}
        endif
    endif
    cd ${wrk_dir} || exit 1
endif # endif AIX

# Successful exit branch
if ( -e dst_mss_bdg.nc ) /bin/cp -f dst_mss_bdg.nc ${MY_DGN_DIR}/${CASE}_bdg.nc
date;pwd;ls -la
if ("${PVM_ARCH}" == "CRAY") then
    ja -sclhft
endif
exit 0

# Error exit branch
err:
date;pwd;ls -la
exit 1
