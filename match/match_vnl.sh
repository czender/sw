#! /bin/csh -f
# -f = fast start, ignore .cshrc and .login

# $Id$

# Purpose: Run vanilla MATCH model on Cray and SGI machines
# Source: http://www.cgd.ucar.edu/csm/working_groups/Atmosphere/scriptccm3.html

# Usage:
# NCAR AIX: llsubmit match_vnl.sh
# NCAR SGI: npri -h 250 match_vnl.sh >dstmch73.txt 2>&1 &
# NCAR Cray: qsub match_vnl.sh 
# UCI Linux: match_vnl.sh >dstmch73.txt 2>&1 &
# UCI SGI: qsub match_vnl.sh

# Debugging: 
# Use totalview for interactive debugger
# Use debugview for diagnosing batch failures
# totalview -c core match -L -include_files
# totalview -c core match

# 20000622: match-3_3_23-dst-0_9_1 T62 pcnst=5 requires 62 hrs wallclock per year with no debugging (except qneg,fixmas) on krein
# 20000620: match_dst T62 pcnst=5 requires 81 hrs wallclock per year with full debugging (qneg,DST_DBG,DST_MSS_BDG,fixmas) on krein
# 19990601: match_dst T42 pcnst=5 requires 216 Mb, 320+ Mb with -g
# 19991210: match_dst T62 pcnst=5 requires 284 Mb prc-1 on SGI

# AIX LoadLeveler batch system
# @ job_name    = dstmch73
# @ output      = $(job_name).txt
# @ error       = $(job_name).txt
# @ job_type    = parallel
# @ network.MPI = css0,shared,us
# @ total_tasks = 1
# @ node        = 1
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
#QSUB -s /bin/csh
#QSUB -J m

# UCI only
#QSUB -q q4
##QSUB -q genq1_4
##QSUB -q pub_8
##QSUB -q pub_16
#QSUB -l mpp_p=8
#QSUB -eo -o /home/ess/zender/zender/match/dstmch73.txt

# NCAR only
##QSUB -q reg
##QSUB -eo -o /fs/cgd/home0/zender/match/dstmch73.txt
##QSUB -lt 10000 -lT 10000
##QSUB -lt 70000 -lT 70000
##QSUB -lt 420000 -lT 420000
##QSUB -lm 280Mb -lM 280Mb

# Main code
set timestamp
set echo
date

# Job-specific
setenv CASE dstmch73
if (${HOST} =~ *.uci.edu) then
    set wrk_dir = /data/zender/match/${CASE}
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute) then
    set wrk_dir = /ptmp/${USER}/${CASE}
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
    set MY_DAT_DIR = ${DATA}/data # Boundary data
    set MY_DGN_DIR = ${DATA}/dgn # Diagnostic output data
    set MY_SRC_DIR = ${HOME}/match_dst # Source code
    /bin/ln -s -f ${HOME}/bin/sh/msread msread || exit 1
    /bin/ln -s -f ${HOME}/bin/sh/mswrite mswrite || exit 1
    /bin/ln -s -f ${HOME}/bin/sh/rdsp rdsp || exit 1
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute) then
    set MY_BIN_DIR = /fs/cgd/home0/${USER}/match/bin/${PVM_ARCH} # Executable
    set MY_DAT_DIR = /fs/cgd/data0/zender/data # Boundary data
    set MY_DGN_DIR = /fs/cgd/data0/${USER}/dgn # Diagnostic output data
    set MY_SRC_DIR = /fs/cgd/data0/${USER}/match_dst # Source code
else
    printf "ERROR: Unknown HOST in run script\n"
    exit
endif # end domain-specific switches

# OS-specific
if ("${PVM_ARCH}" == "AIX") then
    limit stacksize unlimited
    limit datasize  unlimited
    limit memoryuse unlimited
    limit
    setenv OMP_NUM_THREADS 4
else if ("${PVM_ARCH}" == "CRAY") then # endif AIX
    setenv NCPUS 8
else if ("${PVM_ARCH}" == "LINUX") then # endif CRAY
    if (${HOST} =~ dust || ${HOST} =~ biogenic) then
	setenv OMP_NUM_THREADS 2
    else
	setenv OMP_NUM_THREADS 1
    endif # end host-specific switches
else if ("${PVM_ARCH}" == "SGIMP64") then # endif LINUX
    setenv OMP_NUM_THREADS 4
    setenv _DSM_PLACEMENT ROUND_ROBIN # ROUND_ROBIN ensures efficient memory allocation
    #setenv _DSM_VERBOSE # Print diagnostic info about system-level software
    setenv MPC_GANG OFF # Gang scheduling does not work well so always turn off
    setenv TRAP_FPE "UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE" # Ensure abort IEEE errors, requires using LDFLAGS+= -lfpe on SGI
    # Environment-specific
    if ($?ENVIRONMENT) then
	if ( "$ENVIRONMENT" == "BATCH" ) then
	    setenv _DSM_WAIT SPIN # Set to SPIN iff have dedicated access to processors
	    setenv OMP_DYNAMIC FALSE # TRUE = give up processors when machine is busy, use FALSE in dedicated queues
	endif # endif NQS batch job
    endif # endif ENVIRONMENT exists
else # endif SGIMP64
    printf "ERROR: Unknown PVM_ARCH in run script\n"
endif # endif PVM_ARCH

# Cleanup previous run
/bin/rm -f core
/bin/rm -f perf.data
setenv RESTART False
if ("$RESTART" != "True") then
    /bin/rm -f r0???
    /bin/rm -f h0???
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
endif # endif !RESTART

# Assemble new run
/bin/cp -f ${MY_BIN_DIR}/match ${EXE} || exit 1
# SGI: link permanent files to local run directory
# Resolution-dependent part
#/bin/ln -s -f ${MY_DAT_DIR}/dst_T42.nc dst_bnd.nc || exit 1
#/bin/ln -s -f ${MY_DAT_DIR}/dst_T62.nc dst_bnd.nc || exit 1
# fxm: Make this default dataset soon
/bin/ln -s -f ${MY_DAT_DIR}/dst_T62.nc dst_bnd.nc || exit 1

if (${HOST} =~ *.uci.edu) then
    set rpthdyn  = '/DSS' # NCEP dynamics input remote directory
    set lpthdyn = "/data/zender/DSS" # NCEP dynamics input local directory
else if (${HOST} =~ b[bf]*en || ${HOST} =~ dataproc || ${HOST} =~ ute) then
    set rpthdyn  = '/DSS' # NCEP dynamics input remote directory
    set lpthdyn = "/ptmp/${USER}/NCEP" # NCEP dynamics input local directory
    #set rpthdyn  = "/ZENDER/csm/dst01/ccm3/hist/h0001" # CCM dynamics input remote directory
    #set lpthdyn = "/ptmp/${USER}/CCM" # CCM dynamics input local directory
endif # end domain-specific switches
if ( ! -d $lpthdyn ) mkdir -p $lpthdyn
# NCEP dynamics need local GRIB lookup tables
set grb_dir = ${MY_SRC_DIR}/readers/ncep/data
if ( ! -e grb2d.list ) /bin/ln -s ${grb_dir}/grb2d.list grb2d.list
if ( ! -e grbsanl.list ) /bin/ln -s ${grb_dir}/grbsanl.list grbsanl.list
if ( ! -e grib.table2.rean ) /bin/ln -s ${grb_dir}/grib.table2.rean grib.table2.rean
# NCEP at T42 needs these map overlap weights
if ( ! -e T62n.T42n.avg.nc ) /bin/ln -s ${grb_dir}/T62n.T42n.avg.nc T62n.T42n.avg.nc
# src/drydep.F:inidrydep() needs surface type of each grid point
# 1999/09/13: landuse.nc is no longer used by default MATCH
#if ( ! -e landuse.nc ) /bin/ln -s /fs/cgd/csm/people/eaton/data/T42/landuse.nc landuse.nc

# OS-specific run parameters
if ("${PVM_ARCH}" == "AIX" || "${PVM_ARCH}" == "CRAY") then
    set ndens = 4
else if ("${PVM_ARCH}" == "LINUX" || "${PVM_ARCH}" == "SGIMP64") then
# UCI SGI (krein) does not pack tapes correctly
    set ndens = 1
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

cat >! namelist << EOF || exit 1
 &NLIST
 title    = '${CASE} MATCH 4.0beta2 NCEP T62 Vanilla: PRC=S'
 fixmas   = .true.
 qrelax   = 0.0
 delt     = 2400
 wpasswd  = 'passwd'
 rstflg   = 0
 rmout    = .true.
 rmrst    = .true.
 nestep   = 10
 nhtfrq   = 10
 mfilt    = 10
 ndens    = ${ndens}
 irt      = 1000
 rstrrt   = 365
 icdate   = 19980101
 icdatesec    = 0
 rpthdyn  = '${rpthdyn}'
 lpthdyn  = '${lpthdyn}'
 dynoffset= 0
 outinst  = 'ORO','PHIS'
 outtimav = 'U','V','PS','T','Q','RELHUM','CLOUD'
 /
EOF

pwd;ls -la

if ("${PVM_ARCH}" == "CRAY") then
    ja
endif
if ("${PVM_ARCH}" == "AIX") then
    poe ${EXE} < namelist -procs 1 -pgmmodel spmd
else
    ${EXE} < namelist
endif
if ( $status != 0 ) goto err

# Successful exit branch
if ( -e dst_mss_bdg.nc ) /bin/cp -f dst_mss_bdg.nc $MY_DGN_DIR/${CASE}_bdg.nc
date;pwd;ls -la
if ("${PVM_ARCH}" == "CRAY") then
    ja -sclhft
endif
exit 0

# Error exit branch
err:
# CRAY
debugview -c core -B ${EXE}
# SGI
#dbx -I ${MY_SRC_DIR}/src ${EXE}
date;pwd;ls -la
exit 1
