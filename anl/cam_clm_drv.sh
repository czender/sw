# $Id$

# Purpose: Derive fields for all climatological and regional data files

# Usage:
# cd ~/anl;cam_clm_drv.sh
# cd ~/anl;time ./cam_clm_drv.sh > foo.drv 2>&1 &
# cd ~/anl;time ./cam_clm_drv.sh sncpd02 1 1

# This script performs the following operations on all files 
# in the CAM/CLM analysis directory, defined below:

# 1. Derive new fields, as defined in ${HOME}/anl/cam_clm_drv.nco 
# 2. Exclude some variables (FLD_XCL) to conserve memory (optional)

# 0. Command-line options
caseid='sncpd02'
# flg_dst: Similar to flg_sot, but not all options implemented
# 0: Exclude all dust fields
# 1: Include all (FRC+FDB) dust atmosphere and snow fields
flg_dst='1' # [flg] Files contain dust
# flg_sot: 
# 0: Exclude all soot fields
# 1: Include all (FRC+FDB) soot atmosphere and snow fields
flg_sot='1' # [flg] Files contain soot

if [ -n "${1}" ]; then
    caseid=${1}
fi # !$1
if [ -n "${2}" ]; then
    flg_dst=${2}
fi # !$2
if [ -n "${3}" ]; then
    flg_sot=${3}
fi # !$3

# Print initial state
if [ '1' = '1' ]; then # dbg
    printf "${0}: Initialization State:\n"
    printf "${0}: caseid=${caseid}\n"
    printf "${0}: flg_dst=${flg_dst}\n"
    printf "${0}: flg_sot=${flg_sot}\n"
fi # !dbg

# 1. Primary variables

# 2. Derived variables
drc_in=${DATA}/anl_${caseid}
drc_out=${DATA}/anl_${caseid}

# 3. Doubly-derived variables
wrk_dir=${drc_out}

# Create output directory
mkdir -p ${drc_out}

# Go to working directory
cd ${wrk_dir}

##################  3.  ###################
# Exclude variables to conserve memory:
# Note: this may screw up CLM analysis script.
# For NCAR analysis script, use original CLM files
FLD_XCL="ZSOI,DZSOI,SUCSAT,BSW,BTRAN,ERRH2O,ERRSOL,ERRSOI"
FLD_XCL="${FLD_XCL},FCEV,FCTR,FPSN,FSH_V,RSSHA,RSSUN"
FLD_XCL="${FLD_XCL},TV,ZBOT"
#FLD_XCL="${FLD_XCL},QOVER,QDRAI,QRGWL"
FLD_XCL="${FLD_XCL},ELAI,ESAI,FSA,FSH,TSA"
#FLD_XCL="${FLD_XCL},FSDSND,FSDSVD,FSDSNI,FSDSVI"
#FLD_XCL="${FLD_XCL},PRECSL,PRECSC"

# OPTIONAL: REMOVE seasonal timeseries, as these take a long time to process
# rm -f ${drc_in}/${caseid}_ts_DJF.nc
# rm -f ${drc_in}/${caseid}_ts_MAM.nc
# rm -f ${drc_in}/${caseid}_ts_JJA.nc
# rm -f ${drc_in}/${caseid}_ts_SON.nc

/bin/rm -f ${drc_in}/*.tmp
for fl in `ls ${drc_in}`; do
    printf "Deriving new fields for file: ${fl}...\n"

    ncap2 -O -S ${HOME}/anl/cam_clm_drv.nco ${drc_in}/${fl} ${drc_out}/${fl}

    if [ ${flg_dst} = '1' -o ${flg_sot} = '1' ]; then
	ncatted -D 0 -h -O -a missing_value,H2OSNO_TOP,m,f,0.0 -a missing_value,H2OSNO_TOP,d,, ${drc_out}/${fl}
    fi # !flg_dst && !flg_sot

    if [ ${flg_dst} = '1' ]; then
	# Compute surface layer dust concentration
        # Sea-ice SNICAR has one layer so SNODSTMCI=SNODSTMSI
	# Scale soot by 1.0e9 and dust by 1.0e6
	ncatted -D 0 -h -O -a missing_value,SNODSTMSL,m,f,0.0 -a missing_value,SNODSTMSL,d,, ${drc_out}/${fl}
	ncap2 -O -s "SNODSTC = 1.0e6*(SNODSTMSL*landfrac + SNODSTMCI*ICEFRAC)/(H2OSNO_TOP*landfrac + SNOWHICE*1000*ICEFRAC + 1.0e-10f)" ${drc_out}/${fl} ${drc_out}/${fl}
    fi # !flg_dst

    if [ ${flg_sot} = '1' ]; then
	# Computing surface layer BC concentration
        # Sea-ice SNICAR has one layer so SNOBCMCI=SNOBCMSI
	# Scale soot by 1.0e9 and dust by 1.0e6
	ncatted -D 0 -h -O -a missing_value,SNOBCMSL,m,f,0.0 -a missing_value,SNOBCMSL,d,, ${drc_out}/${fl}
	ncap2 -O -s "SNOBCC = 1.0e9*(SNOBCMSL*landfrac + SNOBCMCI*ICEFRAC)/(H2OSNO_TOP*landfrac + SNOWHICE*1000*ICEFRAC + 1.0e-10f)" ${drc_out}/${fl} ${drc_out}/${fl}
    fi # !flg_sot

    # Exclude fields
    ncks -O -x -v ${FLD_XCL} ${drc_out}/${fl} ${drc_out}/${fl}

done # end loop over fl
