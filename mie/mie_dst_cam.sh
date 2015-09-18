# $Id$

# Purpose: Modify default CAM SW aerosol optics file (aeroptics) to contain new 
# dust optical properties, and make CAM LW array easy to edit

# NB: Script intended to be used where SW and LW optical properties are generated
# Usually this is a fast machine like tephra or esmf
# Installation:
# scp ~/mie/mie_dst_cam.sh esmf.ess.uci.edu:mie
# scp ~/mie/mie_dst_cam.sh tephra.ess.uci.edu:mie

# Usage:
# rsync 'tephra.ess.uci.edu:${DATA}/dst/mie/*dst_bln*.nc' ${DATA}/dst/mie
# ~/mie/mie_dst_cam.sh
# chmod a+x ~/mie/mie_dst_cam.sh;~/mie/mie_dst_cam.sh
# scp ${DATA}/inputdata_cam/extra/${fl_cam_xpt} 'esmf.ess.uci.edu:/ptmp/${USER}/inputdata_cam/extra'

aer_cmp_sng='dst_bln_20060904' # [sng] Aerosol composition string
drc_mie="${DATA}/dst/mie" # [drc] Mie output location
drc_cam="${DATA}/inputdata_cam/extra" # [drc] CAM input location
#fl_cam_ctl_stb="AerosolsOptics_cam3019" # [fl] CAM control optical properties
fl_cam_ctl_stb="AerosolsOptics_extended3_BCc" # [fl] CAM control optical properties
fl_cam_ctl="${fl_cam_ctl_stb}.nc" # [fl] CAM control optical properties
fl_cam_xpt="${fl_cam_ctl_stb}_${aer_cmp_sng}.nc" # [fl] CAM experiment optical properties

# Generate new optical properties with mie
# NB: This takes a long time---only re-run when truly necessary
# ${HOME}/mie/psd.pl --dbg=1 --CAM_SW --EMA=20060621 --psd_nbr=4 --dns_prt=2500.0 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_SW.txt 2>&1 &

for wvl_grd in CAM_SW CAM_LW ; do 
    for psd_idx in 01 02 03 04 ; do
# Put optical properties into files in CAM input data directory
	ncks -O -v abs_cff_mss,asm_prm,ext_cff_mss,ss_alb -p ${drc_mie} aer_${aer_cmp_sng}_${psd_idx}_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${psd_idx}_${wvl_grd}.nc
    done # end loop over psd_idx
    
# Concatenate size-bins into single file
    cd ${drc_cam} # Change directory so wildcarding works in next step
    ncecat -O -p ${drc_cam} aer_${aer_cmp_sng}_??_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${wvl_grd}.nc
# Convert mass extinction from [m2 kg-1] to CAM units [m2 g-1] where necessary
# 20060622: aeroptics file for SW code definitely uses [m2 g-1] 
# SW extinctions are converted to [m2 kg-1] for multiplication by aer_mass [kg m-2]
# by "conversion" factor in aer_optics.F90:aer_optics_initialize()
# LW code in aerosol_radiation_interface appears to use m2 kg-1 entire time
# Hence convert ext_cff_mss to [m2 g-1] and leave abs_cff_mss as [m2 kg-1]
# Finally, ensure cooordinate variable wvl_cam is in [um] not [m]
    ncap2 -O -s "wvl=wvl*1.0e6" -s "abs_cff_mss=abs_cff_mss*1.0" -s "ext_cff_mss=ext_cff_mss/1000.0" -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${wvl_grd}.nc
    ncatted -O -a units,ext_cff_mss,o,c,"meter2 gram-1" -a units,abs_cff_mss,o,c,"meter2 kilogram-1" -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc
# Rename variables and dimensions to CAM standards
    ncrename -O -d wvl,wvl_cam -d record,dst -v abs_cff_mss,abs_cam_dust_mod -v asm_prm,asm_cam_dust_mod -v ext_cff_mss,ext_cam_dust_mod -v ss_alb,ssa_cam_dust_mod -v wvl,wvl_cam -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc
# Re-order dimensions
    ncpdq -O -a wvl_cam,dst -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${wvl_grd}.nc
# Remove record dimension?
    ncecat -O -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${wvl_grd}.nc 
    ncwa -O -a record -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc ${drc_cam}/aer_${aer_cmp_sng}_${wvl_grd}.nc 
# Append new attribute to remind people that properties are not CAM-default
    ncatted -O -a source,,o,c,"Optical properties in this variable are not CAM-default. They are developmental and were created by C. Zender using mie." -p ${drc_cam} aer_${aer_cmp_sng}_${wvl_grd}.nc
done # end loop over wvl_grd

# Copy all CAM_SW control properties into CAM_SW experiment file
/bin/cp ${drc_cam}/${fl_cam_ctl} ${drc_cam}/${fl_cam_xpt}
# Place new dust properties in CAM_SW experiment
ncks -A -v asm_cam_dust_mod,ext_cam_dust_mod,ssa_cam_dust_mod -p ${drc_cam} aer_${aer_cmp_sng}_CAM_SW.nc ${drc_cam}/${fl_cam_xpt}

# Diagnostics and sanity checks
if [ "0" = "1" ]; then
    ncks -H -m -v asm_cam_dust_mod,ext_cam_dust_mod,ssa_cam_dust_mod -p ${drc_cam} aer_${aer_cmp_sng}_CAM_SW.nc
    ncks -H -m -v asm_cam_dust_mod,ext_cam_dust_mod,ssa_cam_dust_mod -p ${drc_cam} ${fl_cam_xpt}
    ncbo -O -v asm_cam_dust_mod,ext_cam_dust_mod,ssa_cam_dust_mod -p ${drc_cam} aer_${aer_cmp_sng}_CAM_SW.nc ${fl_cam_xpt} ${drc_cam}/foo.nc
    ncks ${drc_cam}/foo.nc
fi # endif

if [ "1" = "1" ]; then
    # Special treatment for CAM longwave properties on CAM_LW grid
    # CAM3 band 1 is, approximately, H2O non-window, i.e., 0--800+1200--2200 cm-1
    # Compute this as average of 200--800 cm-1 and 1200--2000 properties
    # For mie CAM_LW grid, this is average of bands 1 and 7 (1-based enumeration)
    # Remember that ncap uses C-based indexing and CAM uses Fortran-based
    # NB: Script must always be run from scratch so command does not accumulate!
    # NB: Must hand-edit results into array in aerosol_radiation_interface.F90
    ncap2 -O -s 'abs_cam_dust_mod(0,:)=0.5*(abs_cam_dust_mod(0,:)+abs_cam_dust_mod(6,:))' -p ${drc_cam} aer_${aer_cmp_sng}_CAM_LW.nc aer_${aer_cmp_sng}_CAM_LW.nc

    cdl_sng='netcdf foo_lw {
    dimensions:
    wvl_cam=7,naer_ttl=12;
    variables:
    float abs_cff_mss_aer(wvl_cam,naer_ttl);
    float abs_cff_mss_vlc(wvl_cam);
    data:
    abs_cff_mss_vlc=70.257384,285.282943,1.0273851e+02,6.3073303e+01,1.2039569e+02,3.6343643e+02,2.7138528e+02;
    }' # end cdl_sng
    echo ${cdl_sng} > ~/foo_lw.cdl
    ncgen -b -o ~/foo_lw.nc ~/foo_lw.cdl
    ncks -A -C -v abs_cff_mss_aer,abs_cff_mss_vlc ~/foo_lw.nc ${drc_cam}/aer_${aer_cmp_sng}_CAM_LW.nc
    ncap2 -O -s 'abs_cff_mss_aer(:,:)=0.0;abs_cff_mss_aer(:,2:5)=abs_cam_dust_mod;abs_cff_mss_aer(:,11)=abs_cff_mss_vlc' -p ${drc_cam} aer_${aer_cmp_sng}_CAM_LW.nc ${drc_cam}/aer_${aer_cmp_sng}_CAM_LW.nc
    /bin/rm -f ~/foo.txt
    printf "!     Dust Longwave Optical properties for %s\n" ${aer_cmp_sng} > ~/foo.txt
    date_sng=`date`
    printf "!     Generated by ${USER} with mie, mie_dst_cam.sh on ${date_sng}\n" >> ~/foo.txt
    printf '      real(r8), public, parameter, dimension (naer_all,bnd_nbr_LW)::  &\n' >> ~/foo.txt
    printf '       abs_cff_mss_aer =  reshape ( (/ &\n' >> ~/foo.txt
    for wvl_idx in 1 2 3 4 5 6 7 ; do
	printf '                   ' >> ~/foo.txt
	ncks -H -F -C -s '%15.9e, ' -v abs_cff_mss_aer -d wvl_cam,${wvl_idx} -p ${drc_cam} aer_${aer_cmp_sng}_CAM_LW.nc >> ~/foo.txt
	printf '&\n' >> ~/foo.txt
    done
    printf '                   /),(/naer_all,bnd_nbr_LW/) )' >> ~/foo.txt

    ncks -F -v abs_cam_dust_mod,abs_cff_mss_aer -p ${drc_cam} aer_${aer_cmp_sng}_CAM_LW.nc > ~/foo
fi # endif

# Copy new CAM-conformant aeroptics file to ESMF and IPCC
scp ${drc_cam}/${fl_cam_xpt} 'esmf.ess.uci.edu:/ptmp/${USER}/inputdata_cam/extra'
scp ${drc_cam}/${fl_cam_xpt} 'ipcc.ess.uci.edu:/ptmp/${USER}/inputdata_cam/extra'

# Diagnostics and sanity checks
if [ "0" = "1" ]; then
    ncks -F -d wvl_cam,0.5 -C -H -v ssa_cam_dust,ssa_cam_dust_mod ${DATA}/inputdata_cam/extra/AerosolsOptics_cam3019.nc
    ncks -F -d wvl_cam,0.5 -C -H -v ssa_cam_dust,ssa_cam_dust_mod ${DATA}/inputdata_cam/extra/AerosolsOptics_extended3_BCc_dst_bln_20060621.nc
    ncks -F -d wvl_cam,0.5 -C -H -v ssa_cam_dust,ssa_cam_dust_mod ${DATA}/inputdata_cam/extra/AerosolsOptics_extended3_BCc_dst_bln_20060904.nc
    ncks -F -d wvl_cam,0.5 -C -H -v ssa_cam_dust,ssa_cam_dust_mod ${DATA}/inputdata_cam/extra/AerosolsOptics_extended3_BCc_dst_bln_20060907.nc
    for psd_idx in 01 02 03 04 ; do
	ncks -F -d wvl,0.5e-6 -C -H -v ss_alb ${DATA}/dst/mie/aer_dst_bln_20060907_${psd_idx}_CAM_SW.nc
    done # end loop over psd_idx
    for psd_idx in 01 02 03 04 ; do
	ncks -F -u -d wvl,10.0e-6 -C -H -v abs_cff_mss ${DATA}/dst/mie/aer_dst_bln_20060907_${psd_idx}_CAM_LW.nc
    done # end loop over psd_idx
fi # endif False
