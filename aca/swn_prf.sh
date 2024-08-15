#!/usr/bin/env bash
# Purpose: Simulate standard profiles with standard surface reflectances
# Usage: swn_prf.sh [no options yet]
# Required Environment:
# DATA_RT is SWNB2 input data directory (typically ~zender/data/aca)
# DATA is grandparent of SWNB2 output directory (typically ~zender/data/rrtmgp/slr_Kur95)

# Define spectral indices
idx_vsb_min=1296 # Bluest SWNB2 band fully in RRTMGP (1-based) band 5
idx_vsb_max=1466 # Bluest SWNB2 band fully in EAM VIS
idx_nir_min=1467 # Reddest SWNB2 band fully in EAM NIR
idx_nir_max=1610 # Reddest SWNB2 band fully in RRTMGP (1-based) band 5

# Specify solar spectrum
#spc=Kur95
spc=FDE24

# Define I/O directories
drc_in=${DATA_RT}
drc_out=${DATA}/rrtmgp/spc_${spc}
spt_src="${BASH_SOURCE[0]}"
spt_nm=$(basename ${spt_src}) # [sng] Script name (unlike $0, ${BASH_SOURCE[0]} works well with 'source <script>')

# Main loop over profiles
for prf in mls mlw sas saw std tro; do
#for prf in std; do

    # Loop over clear/cloudy columns
    for sky in clr cld; do

	mkdir -p ${drc_out}
	
	# Select surface reflectance appropriate for profile
	case ${prf} in
	    tro* ) rfl_sfc='--rfl=rfl_spc_sfc_ocn_chl0.3mgxm3_csza0.5.nc' ; ;;
	    mlw* | sas* | saw* ) rfl_sfc='--rfl=rfl_spc_sfc_snc_100um_csza0.5.nc' ; ;;
	    mls* | std* ) rfl_sfc='--alb_sfc_vsb=0.1 --alb_sfc_NIR=0.1' ; ;;
	    * ) printf "${spt_nm}: ERROR profile = prf = ${prf} has no associated surface reflectance\n" ; ;;
	esac # !prf
	
	# Select cloud thickness appropriate for sky-type
	case ${sky} in
	    clr* ) mpc_CWP=0.0 ; ;;
	    cld* ) mpc_CWP=0.1 ; ;;
	    * ) printf "${spt_nm}: ERROR sky = ${sky} has no associated cloud thickness\n" ; ;;
	esac # !sky
	
	# Select appropriate spectral flux
	case ${spc} in
	    Kur95* ) spc_slr="--spc_slr=spc_${spc}.nc" ; ;;
	    FDE24* ) spc_slr="--spc_slr=spc_${spc}.nc" ; ;;
	    * ) printf "${spt_nm}: ERROR spc = ${spc} has no associated spectral file\n" ; ;;
	esac # !sky
	
	printf "\n${spt_nm}: Simulating/analysing atmosphere with profile = ${prf}, sky = ${sky}, spc = ${spc}\n"

	# Simulate atmospheric RT
	if true; then
	    cmd_swn="swnb2 --drc_in=${drc_in} --lqd=aer_h2o_lqd_rds_swa_10.nc --prf=${prf}_icrccm_50${sky}.nc --mpc_CWP=${mpc_CWP} ${rfl_sfc} ${spc_slr} --out=${drc_out}/swn_${prf}_${sky}.nc > ${drc_out}/swn_${prf}_${sky}.txt 2>&1"
	    echo ${cmd_swn}
	    eval ${cmd_swn}
	    if [ "$?" -ne 0 ]; then
		printf "${spt_nm}: ERROR Failed to simulate atmosphere. Debug this:\n${cmd_swn}\n"
		exit 1
	    fi # !err
	fi # !false
	
	# Reverse spectral dimension so it increases (rather than decreases) with wavelength
	if true; then
	    cmd_pdq="ncpdq -O --rdr=-bnd ${drc_out}/swn_${prf}_${sky}.nc ${drc_out}/swn_${prf}_${sky}.nc"
	    echo ${cmd_pdq}
	    eval ${cmd_pdq}
	    if [ "$?" -ne 0 ]; then
		printf "${spt_nm}: ERROR Failed to reverse spectral dimension. Debug this:\n${cmd_pdq}\n"
		exit 1
	    fi # !err
	fi # !false
	
	# Hyperslab file to retain only RRTMGP Band 5 wavelengths, and select scalars
	if true; then
	    cmd_hyp="ncks -O -d bnd,${idx_vsb_min},${idx_nir_max} -v flx_spc_.?,wvl_.?,wvn_.? ${drc_out}/swn_${prf}_${sky}.nc ${drc_out}/swn_${prf}_${sky}_bnd5.nc"
	    echo ${cmd_hyp}
	    eval ${cmd_hyp}
	    if [ "$?" -ne 0 ]; then
		printf "${spt_nm}: ERROR Failed to hyperslab spectral dimension. Debug this:\n${cmd_hyp}\n"
		exit 1
	    fi # !err
	fi # !false
	
	# Analyze SWNB2 simulation
	# Currently borken on MacOS not Linux, not sure why
	if true; then
	    cmd_anl="ncap2 -O -v --fl_spt=${HOME}/sw/aca/swn_rrtmgp_anl_bnd5.nco ${drc_out}/swn_${prf}_${sky}.nc ${drc_out}/swn_${prf}_anl_bnd5_${sky}.nc"
	    if [ "${PVM_ARCH}" != 'MACOS' ]; then
		echo ${cmd_anl}
		eval ${cmd_anl}
		if [ "$?" -ne 0 ]; then
		    printf "${spt_nm}: ERROR Failed to analyze SWNB2 simulation. Debug this:\n${cmd_anl}\n"
		    exit 1
		fi # !err
	    else
		# ...instead, on MacOS, mouse this into terminal to process all data
		printf "${spt_nm}: WARNING MacOS requires manual analysis...mouse this into terminal window at end:\n"
		printf "for prf in mls mlw sas saw std tro; do\n\tfor sky in clr cld; do\n\t\tncap2 -O -v --fl_spt=${HOME}/sw/aca/swn_rrtmgp_anl_bnd5.nco ${drc_out}/swn_\${prf}_\${sky}.nc ${drc_out}/swn_\${prf}_anl_bnd5_\${sky}.nc\n\t\tncks -v flx_frc_d.._vsb ${drc_out}/swn_\${prf}_anl_bnd5_\${sky}.nc\n\tdone\ndone\n"
	    fi # !PVM_ARCH
	fi # !false
	
    done # !sky

done # !prf
    
exit 0
    
