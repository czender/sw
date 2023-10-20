#!/usr/bin/env bash
# Purpose: Simulate standard profiles with standard surface reflectances
# Usage: swn_prf.sh [no options yet]
# Required Environment:
# DATA_RT is SWNB2 input data directory (typically ~zender/data/aca)
# DATA is grandparent of SWNB2 output directory (typically ~zender/data/rrtmgp/swnb2)

# Define I/O directories
drc_in=${DATA_RT}
drc_out=${DATA}/rrtmgp/swnb2
spt_src="${BASH_SOURCE[0]}"
spt_nm=$(basename ${spt_src}) # [sng] Script name (unlike $0, ${BASH_SOURCE[0]} works well with 'source <script>')

# Main loop over profiles
#for prf in mls mlw sas saw std tro; do
for prf in tro; do
    mkdir -p ${drc_out}

    # Select surface reflectance appropriate for profile
    case ${prf} in
	tro* ) rfl_sfc='--rfl=rfl_spc_sfc_ocn_chl0.3mgxm3_csza0.5.nc' ; ;;
	mlw* | sas* | saw* ) rfl_sfc='--rfl=rfl_spc_sfc_snc_100um_csza0.5.nc' ; ;;
	mls* | std* ) rfl_sfc='--alb_sfc_vsb=0.1 --alb_sfc_NIR=0.1' ; ;;
	* ) printf "${spt_nm}: ERROR profile = prf = ${prf} has no associated surface reflectance\n" ; ;;
    esac # !prf

    # Simulate atmosphere
    cmd_swn="swnb2 --drc_in=${drc_in} --lqd=aer_h2o_lqd_rds_swa_10.nc --prf=${prf}_icrccm_50lvl.nc ${rfl_sfc} --out=${drc_out}/swn_${prf}_clr.nc > ${drc_out}/swn_${prf}_clr.txt 2>&1"
    eval ${cmd_swn}
    if [ "$?" -ne 0 ]; then
	printf "${spt_nm}: ERROR Failed to simulate atmosphere. Debug this:\n${cmd_swn}\n"
	exit 1
    fi # !err

    # Reverse spectral dimension so it increases (rather than decreases) with wavelength
    cmd_pdq="ncpdq -O --rdr=-bnd ${drc_out}/swn_${prf}_clr.nc ${drc_out}/swn_${prf}_clr.nc >> ${drc_out}/swn_${prf}_clr.txt 2>&1"
    eval ${cmd_pdq}
    if [ "$?" -ne 0 ]; then
	printf "${spt_nm}: ERROR Failed to reverse spectral dimension. Debug this:\n${cmd_pdq}\n"
	exit 1
    fi # !err

    # Analyze SWNB2 simulation
    cmd_anl="ncap2 -O -v --fl_spt ${HOME}/sw/aca/swn_rrtmgp_bnd5.nco ${drc_out}/swn_${prf}_clr.nc ${drc_out}/swn_${prf}_bnd5_clr.nc >> ${drc_out}/swn_${prf}_clr.txt 2>&1"
    # Allow antlr_ret=133 for NCO < 5.1.9
    if [ "$?" -ne 0 ] && [ "$?" -ne 133 ] ; then
	printf "${spt_nm}: ERROR Failed to analyze SWNB2 simulation. Debug this:\n${cmd_anl}\n"
	exit 1
    fi # !err

done # !prf

exit 0
