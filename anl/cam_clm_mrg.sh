# $Id$

# Purpose: Merge subsets CAM and CLM variables into common monthly data

# Testing:
# /bin/cp -p -r ${DATA}/sncpd10 ${DATA}/sncpd10_bck
# /bin/rm -r ${DATA}/sncpd10;/bin/cp -p -r ${DATA}/sncpd10_bck ${DATA}/sncpd10
# /bin/rm -r ${DATA}/sncpd10;mkdir ${DATA}/sncpd10;/bin/cp -p ${DATA}/sncpd10_bck/sncpd10*.h0.002[01]*.nc ${DATA}/sncpd10

# Distribution:
# scp ~/anl/cam_clm_mrg.sh dust.ess.uci.edu:anl

# Production:
# cd ~/anl;./cam_clm_mrg.sh
# cd ~/anl;time ./cam_clm_mrg.sh > foo.mrg 2>&1 &
# cd ~/anl;time ./cam_clm_mrg.sh sncpd02 0000 0020 1 1 1
# cd ~/anl;time ./cam_clm_mrg.sh sncpd02 0000 0020 1 1 1 > ~/foo.mrg 2>&1 &

# set -o xtrace # dbg

# 0. Command-line options
caseid='sncpd02'
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
yr_srt='0000' # [yr] Start year (NB: Include year 0 for now)
yr_end='0020' # [yr] End year

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
# NB: Turn this on to process control runs which lack FLD_FK variables
# Currently their are no variables in experiments which are not in controls
# so this is a no-op
flg_fk='1' # [flg] Fake variable SNODSTFRC=SNOBCFRC=0.0

# 2. Derived variables
drc_in="${DATA}/${caseid}"
drc_out="${DATA}/${caseid}"
wrk_dir="${DATA}/${caseid}"

# Go to working directory
cd ${wrk_dir}

# Define CAM fields to add to CLM output
FLD_APN="FLNT,FLUT,FSNT,FSNTOA,FSNTOAC,FLNS,LHFLX,SHFLX,SRFRAD"
#FLD_APN="${FLD_APN},gw"
FLD_APN="${FLD_APN},FSDS,FSNS"
FLD_APN="${FLD_APN},SOLL,SOLLD,SOLS,SOLSD"
FLD_APN="${FLD_APN},PRECC,PRECL,PRECSC,PRECSL,SNOWHLND,SNOWHICE"
FLD_APN="${FLD_APN},ICEFRAC,OCNFRAC,TREFHT,TS,PS"
FLD_APN="${FLD_APN},FSNSC,FSNTC,FSDSC,FSNSOI,SOLIN"
FLD_APN="${FLD_APN},CLDLOW,CLDMED,CLDHGH,CLDTOT"
FLD_APN="${FLD_APN},FLNSC,FLNTC,LWCF"
# 20060829: New aerosol diagnostics
#FLD_APN="${FLD_APN},fxm"

# CAM/CSIM fields for SNICAR
if [ ${flg_snc} = '1' ]; then
    FLD_APN="${FLD_APN},FSNSVD,FSNSVI,FSNSND,FSNSNI"
    FLD_APN="${FLD_APN},SNORDSI,SNOdTdzI,SNOLIQI"
fi
# CAM/CSIM fields for Dust
if [ ${flg_dst} = '1' ]; then
# Soot-like fields first, then dust-specific fields, outfld() routine in brackets
    FLD_APN="${FLD_APN},DUSTOD_v" # Total visible optical depth [radiation.F90:radiation_tend()]
    FLD_APN="${FLD_APN},SNODSTMCI" # Mass Path Column of dust in snow on sea-ice
    FLD_APN="${FLD_APN},SNODSTFRCI" # Forcing by dust in snow on sea-ice
    FLD_APN="${FLD_APN},DSTODXC" # Optical depth (zero)
#    FLD_APN="${FLD_APN},DSTDEP" # Dust deposition (dry+wet) from atmosphere [dust_intr.F90:dust_drydep_intr(), dust_emis_intr(), dust_wet_intr()]
    FLD_APN="${FLD_APN},DSTSFDRY,DSTSFMBL,DSTSFWET" # Total Dry, Mobilization, Wet fluxes, [dust_intr.F90:dust_drydep_intr(), dust_emis_intr(), dust_wet_intr()]
    FLD_APN="${FLD_APN},DSTX01DD,DSTX02DD,DSTX03DD,DSTX04DD" # Dry Deposition (Gravitational+Turbulent) [dust_intr.F90:dust_drydep_intr()]
    FLD_APN="${FLD_APN},DSTX01WD,DSTX02WD,DSTX03WD,DSTX04WD" # Wet Deposition [dust_intr.F90:dust_wet_intr()]
    FLD_APN="${FLD_APN},DSTX01MP,DSTX02MP,DSTX03MP,DSTX04MP" # Mass Path [aerosol_radiation_interface.F90:aerosol_diagnostics()]
#    FLD_APN="${FLD_APN},DSTX01,DSTX02,DSTX03,DSTX04" # Mixing ratio NB: 4-D []
#    FLD_APN="${FLD_APN},DSTX01DT,DSTX02DT,DSTX03DT,DSTX04DT" # Deposition tendency NB: 4-D [dust_intr.F90:dust_drydep_intr()]
#    FLD_APN="${FLD_APN},DSTX01DV,DSTX02DV,DSTX03DV,DSTX04DV" # Deposition velocity NB: 4-D [dust_intr.F90:dust_drydep_intr()]
#    FLD_APN="${FLD_APN},DSTX01PP,DSTX02PP,DSTX03PP,DSTX04PP" # Precipitation tendency NB: 4-D [dust_intr.F90:dust_wet_intr()]
    FLD_APN="${FLD_APN},DSTX01GV,DSTX02GV,DSTX03GV,DSTX04GV" # Gravitational setting [dust_intr.F90:dust_drydep_intr()]
    FLD_APN="${FLD_APN},DSTX01OD,DSTX02OD,DSTX03OD,DSTX04OD" # Optical depth [radiation.F90:radiation_tend()]
    FLD_APN="${FLD_APN},DSTX01SF,DSTX02SF,DSTX03SF,DSTX04SF" # Emissions [dust_intr.F90:dust_emis_intr()]
    FLD_APN="${FLD_APN},DSTX01TB,DSTX02TB,DSTX03TB,DSTX04TB" # Turbulent depostion [dust_intr.F90:dust_drydep_intr()]
    FLD_APN="${FLD_APN},SFDSTX01,SFDSTX02,SFDSTX03,SFDSTX04" # Total Deposition (Dry+Wet) (zero)
fi # !flg_dst
# CAM/CSIM fields for Soot
if [ ${flg_sot} = '1' ]; then
    FLD_APN="${FLD_APN},BCOD_v,OCOD_v,AEROD_v,AERSSA_v,AERASM_v"
    FLD_APN="${FLD_APN},BCDEPDRY,BCDEPWET,SNOBCMCI"
    FLD_APN="${FLD_APN},FSNS_RF,FSNT_RF,FSNSC_RF,FSNTC_RF"
    FLD_APN="${FLD_APN},SNOBCFRCI,SNOAERFRCI"
elif [ ${flg_sot} = '2' ]; then
    FLD_APN="${FLD_APN},BCOD_v,OCOD_v,AEROD_v,AERSSA_v,AERASM_v"
    FLD_APN="${FLD_APN},FSNS_RF,FSNT_RF,FSNSC_RF,FSNTC_RF"
    FLD_APN="${FLD_APN},BCDEPDRY,BCDEPWET"
elif [ ${flg_sot} = '3' ]; then
    FLD_APN="${FLD_APN},BCOD_v,OCOD_v,AEROD_v,AERSSA_v,AERASM_v"
    FLD_APN="${FLD_APN},BCDEPDRY,BCDEPWET,SNOBCMCI"
elif [ ${flg_sot} = '4' ]; then
    FLD_APN="${FLD_APN},BCOD_v,OCOD_v,AEROD_v,AERSSA_v,AERASM_v"
elif [ ${flg_sot} = '5' ]; then
    FLD_APN="BCDEPDRY,BCDEPWET"
fi # !flg_sot

##################  2.  ###################
# Variables to merge from both CLM and CAM
# Variable name must conform to: VAR_I in CAM and VAR_L in CLM 
# FLD_MRG1: re-establish zero as missing value after merging
# FLD_MRG2: do not re-establish missing value
# SNICAR fields
if [ ${flg_snc} = '1' ]; then
    FLD_MRG1="SNORDS SNOdTdz"
fi # !flg_snc
# Aerosol fields:
FLD_MRG2=""
if [ ${flg_sot} = '1' ]; then
    FLD_MRG2="${FLD_MRG2} SNOBCFRC SNOAERFRC"
fi # !flg_sot
if [ ${flg_dst} = '1' ]; then
    FLD_MRG2="${FLD_MRG2} SNODSTFRC"
fi # !flg_dst

##################  3.  ###################
# Change radiative flux variables to double precision to prevent underflow
# (necessary for single fire emission events)
FLD_DBL=""
if [ ${flg_dst} = '1' -o ${flg_sot} = '1' -o ${flg_sot} = '2' ]; then
    FLD_DBL="FSNT FSNT_RF FSNTC FSNTC_RF FSNS FSNS_RF FSNSC FSNSC_RF"
fi # !flg_dst && !flg_sot
if [ ${flg_dst} = '1' ]; then
    FLD_DBL="${FLD_DBL} SNODSTFRCL SNODSTFRCI SNODSTFRC"
fi # !flg_dst
if [ ${flg_sot} = '1' ]; then
    FLD_DBL="${FLD_DBL} SNOBCFRCL SNOBCFRCI SNOBCFRC SNOAERFRCL SNOAERFRCI SNOAERFRC"
fi # !flg_sot

unset dpm # Days per month
declare -a dpm
dpm=(0 31 28.25 31 30 31 30 31 31 30 31 30 31) # Allows 1-based indexing

for yr in `seq $yr_srt $yr_end`; do
    YYYY=`printf "%04d" ${yr}`
    for mth in {01..12}; do
	printf "Appending new variables to CLM output files for year ${YYYY}, month ${mth}...\n"
	MM=`printf "%02d" ${mth}`
	# Append days_per_month
	ncap2 -O -s "dpm=${dpm[${mth}]}" ${drc_in}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	# Appropriate place to add "fake" variables that make datasets subtractable
	if [ ${flg_fk} = '1' ]; then
	    echo 'Hello world' > /dev/null;
	fi # !flg_fk
	# Append CAM variables to CLM files
	ncks -A -C -v ${FLD_APN} ${drc_in}/${caseid}.cam2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	# Multiply aerosol surface forcing on land by landfrac
        # (forcing on sea-ice is interactively multiplied by icefrac, and operates over snow-covered and bare sea-ice)
        # (forcing on land is interactively weighted by snow cover fraction)
	if [ ${flg_dst} = '1' ]; then
	    ncap2 -O -s 'SNODSTFRCL=SNODSTFRCL*landfrac' ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	fi
	if [ ${flg_sot} = '1' ]; then
	    ncap2 -O -s 'SNOBCFRCL=SNOBCFRCL*landfrac' -s 'SNOAERFRCL=SNOAERFRCL*landfrac' ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	fi
	# Derive sea-ice snow cover fraction:
	ncap2 -O -s 'FSNOI=SNOWHICE*3.0303f/(SNOWHICE*3.0303f+0.02f)' ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc

	# Merge land and ice fields into single field:
	for fld in ${FLD_MRG1}; do
	    # Zero missing values for land and sea-ice, remove missing_value attribute
	    # NB: Order matters for these -a att_edits
	    ncatted -D 0 -h -O -a missing_value,${fld}I,m,f,0.0 -a missing_value,${fld}L,m,f,0.0 -a missing_value,${fld}I,d,, -a missing_value,${fld}L,d,, ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	
	    # Combine land with sea-ice fields, weighting by land_frac and ice_frac
	    # Do not weight by FSO and FSNOI, because then the answer would be sensitive to snow cover fraction
	    #ncap2 -O -s "${fld}=(${fld}L*landfrac*FSNO+${fld}I*ICEFRAC*FSNOI)/(landfrac*FSNO+ICEFRAC*FSNOI+1.0e-10f)" ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	    ncap2 -O -s "${fld}=(${fld}L*landfrac+${fld}I*ICEFRAC)/(landfrac+ICEFRAC +1.0e-10f)" ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	
            # Re-establish missing_values for the original land and sea-ice fields
	    # NB: Order matters for these -a att_edits
	    ncatted -D 0 -h -O -a missing_value,${fld}I,c,f,0.0 -a missing_value,${fld}L,c,f,0.0 -a missing_value,${fld}I,m,f,1.0e36 -a missing_value,${fld}L,m,f,1.0e36 ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc

	    # Establish missing_value for merged land/ice fields
	    # NB: Order matters for these -a att_edits
	    ncatted -D 0 -h -O -a missing_value,${fld},c,d,0.0 -a missing_value,${fld},m,d,1.0e36 ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	done # end loop over fld

	# Aerosol fields
	for fld in ${FLD_MRG2}; do
	    # Zero all missing values for land and sea-ice, remove missing_value attribute
	    # NB: Order matters for these -a att_edits
	    ncatted -D 0 -h -O -a missing_value,${fld}L,m,f,0.0 -a missing_value,${fld}L,d,, ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	
	    # Combine land sea-ice fields
	    # Fields have already been multiplied by respective snow fractions,
	    # so they are simply added together
	    ncap2 -O -s "${fld}=${fld}L+${fld}I" ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	done # end loop over fld

	# Set certain variables to double precision
	# Be careful to change missing value precision consistently
	for fld in ${FLD_DBL}; do
	    ncap2 -O -s "${fld}=double(${fld})" ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc ${drc_out}/${caseid}.clm2.h0.${YYYY}-${MM}.nc
	done
    done # end loop over mth
done # end loop over yr
