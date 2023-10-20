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
#for prf in mls mlw sas saw std tro; do
for prf in tro; do
    mkdir -p ${drc_out}
    if [ ${prf} = 'tro' ]; then
       rfl_sfc=rfl_spc_sfc_snc_100um_csza0.5.nc
    else
       rfl_sfc=rfl_spc_sfc_ocn_chl0.0_csza0.5.nc
    fi # !prf
    # Simulate atmosphere
    cmd_swn="swnb2 --drc_in=${drc_in} --lqd=aer_h2o_lqd_rds_swa_10.nc --prf=${prf}_icrccm_50lvl.nc --rfl=${rfl_sfc} --out=${drc_out}/swn_${prf}_clr.nc > ${drc_out}/swnb2_${prf}_clr.txt 2>&1"
    eval ${cmd_swn}
    if [ "$?" -ne 0 ]; then
printf "${spt_nm}: ERROR Failed to simulate atmosphere. Debug this:\n${cmd_swn}\n"
exit 1
    fi # !err
    # Re-order wavelengths (optional)
    cmd_pdq="ncpdq -O --rdr=-bnd ${drc_out}/swn_${prf}_clr.nc ${drc_out}/swn_${prf}_clr.nc"
    eval ${cmd_pdq}
    if [ "$?" -ne 0 ]; then
printf "${spt_nm}: ERROR Failed to reverse dimensions. Debug this:\n${cmd_pdq}\n"
exit 1
    fi # !err
done
