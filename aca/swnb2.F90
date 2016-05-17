! $Id$

program swnb2
  
  ! Purpose: Narrow band model employing Malkmus random line approximation
  ! Effective monochromatic optical depth from gaseous line absorption is combined  
  ! with monochromatic continuum absorption and scattering cross sections.
  ! DISORT computes discrete ordinates radiative transfer solution
  ! I/O uses netCDF interface

  ! Copyright (C) 1994--2016 Charlie Zender
  ! License: GNU General Public License (GPL) Version 3
  ! See http://www.gnu.org/copyleft/gpl.html for full license text
  ! The original author of this software, Charlie Zender, seeks to improve
  ! it with your suggestions, contributions, bug-reports, and patches.
  ! Charlie Zender <zender at uci dot edu>
  ! Department of Earth System Science
  ! University of California, Irvine
  ! Irvine, CA 92697-3100
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W swnb2.F90 OMP=N OPTS=D NETCDF4=Y swnb2; cd - # grele, neige, virga
  ! cd ~/aca;NETCDF_ROOT='/usr' NETCDF_INC='/usr/include' NETCDF_LIB='/usr/lib' make NETCDF4=N NETCDFF=Y; cd - # netCDF3 givre, neige
  ! cd ${HOME}/aca; make -W swnb2.F90 MPI=Y OMP=N OPTS=D NETCDF4=Y swnb2; cd - # greenplanet
  ! cd ${HOME}/aca; make -W swnb2.F90 OMP=N OPTS=D NETCDF4=N swnb2; cd - # sand
  ! Debugging compilation:
  ! cd ${HOME}/f  ;make cln;make OMP=N OPTS=X NETCDF4=Y
  ! cd ${HOME}/aca;make cln;make OMP=N OPTS=X NETCDF4=Y
  ! Production compilation:
  ! cd ${HOME}/f  ;make cln;make OMP=Y OPTS=O
  ! cd ${HOME}/aca;make cln;make OMP=N OPTS=O
  ! scp ~/aca/swnb2.F90 esmf.ess.uci.edu:aca
  ! scp esmf.ess.uci.edu:aca/swnb2.F90 ~/aca

  ! NB: Compile and run swnb2 in double precision
  ! swnb2 is intimately linked to fortran library libcsz_f90
  ! Compile these two pieces with same intrinsic precision
  ! AIX/xlf underflows on attempts to read small O(10^-40) abs_xsx_O2O2
  ! (produced by O2X) into single-precision variable
  ! Long term solution is generic interface to SP and DP versions of rbn_vec(), ntp_vec()
  ! Then always read and re-grid O2X files in double precision, independently
  ! of rest of program
  ! Short term solution is do not trap underflows in single precision

  ! Graphics: 
  ! IDL procedure ~/idl/mie.pro:odx_htg_gph() plots heating profiles
  ! NCL procedure ~/ncl/odxc.ncl: plots spectral optical depths
  ! IDL procedure ~/idl/mie.pro:trn_abs_gph() plots transmission vs. absorption
  ! Execution scripts:
  ! ~/dst/swnb.sh: Run swnb2 for comparison to crm aerosol radiative forcing
  ! ~/aca/aca.pl: Run swnb2 forced with ARESE and IOP data
  
  ! Nomenclature:
  ! Introducing snow causes some semantic difficulties
  ! Formerly, TOA and sfc were interfaces bounding atm
  ! Now, atm means different things in different places:
  ! _atm suffix denotes the column above the snowpack when 
  ! applied to interface levels, e.g., lev_atm and levp_atm, 
  ! _atm suffix denotes entire column between TOA and sfc when
  ! applied to archived fluxes etc., e.g., flx_abs_atm.
  ! _snw suffix denotes entire snowpack when
  ! applied to archived fluxes etc., e.g., flx_abs_snw.
  ! _snw suffix denotes properties of pure snow when
  ! applied to archived species-specific optical properties, e.g., odxc_obs_snw
  ! This confusing state of affairs arose because I did not want
  ! to change the physical definitions of atm when implementing snw.
  ! Guessing a variable's definition from its name is now ambiguous!

  ! For now, we have:
  ! TOA      <- interface
  ! 	atm  <- layer
  ! snw      <- interface
  ! 	snw  <- layer
  ! sfc      <- interface
  ! so that:
  ! TOA and sfc are generally considered interfaces
  ! atm is generally considered a layer, either to the snow or to the surface
  ! snw is considered both a layer and an interface

  ! In the future, we may want to be more precise, e.g.:
  ! 1. Change sfc to mean snow top
  ! 2. Use gnd to mean ground
  ! 3. Use snp (as distinct from snw) to refer to vertically extensive snowpack
  ! 4. Use snw to refer top snow interface when present, ground when not
  ! In thinking through future nomenclature changes, priorities are:
  ! 1. Meaning of _atm does not extend into snowpack. _atm refers to above-snowpack.
  ! 2. Meaning of _atm does not extend into snowpack. _atm refers to above-snowpack.

  ! Command line options:
  ! -a aerosol file
  ! -A turn off aerosol
  ! -b background aerosol file
  ! -B turn off background aerosol
  ! -c, --CO2 fl_CO2: CO2 file
  ! -C turn off CO2
  ! -d fl_out ! [sng] Output file
  ! --drc_in drc_in ! [sng] Input directory
  ! --drc_out drc_out ! [sng] Output directory
  ! -D debug level
  ! -e band index for single band computation
  ! -E single band computation
  ! -f force liquid phase
  ! -F force ice phase
  ! -g, --CH4 fl_CH4: CH4 file
  ! -G turn off CH4
  ! -h, --H2O fl_H2O: H2O file
  ! -H turn off H2O
  ! -i --ice fl_ice: Ice file
  ! -I turn off Ice
  ! -j use non-Lambertian surface
  ! -J turn off thermal source
  ! -k --O2O2 fl_O2O2: O2O2 file
  ! -K turn off O2O2
  ! -l, --lqd fl_lqd: Liquid file
  ! -L turn off Liquid
  ! -m cloud water path
  ! -M aerosol optical depth
  ! -n instrument filter file
  ! -N turn off instrument filter
  ! -o, --O2O2: file
  ! -O turn off O2
  ! -p column file
  ! -P Henyey-Greenstein test case tst_case_HG
  ! -q, --H2OH2O: H2OH2O file
  ! -Q turn _on_ (!!!) H2OH2O
  ! -r surface albedo
  ! -R turn off Rayleigh scattering
  ! -s number of streams
  ! -S solar constant W m-2
  ! -t Rayleigh test case tst_case_Rayleigh
  ! -T solar spectrum file
  ! -u number of user polar angles
  ! -U turn off O2N2
  ! -v surface properties file
  ! -V BRDF type
  ! -w O3, --O3: file
  ! -W turn off O3
  ! -x, --NO2: NO2 file
  ! -X turn off NO2
  ! -y, --OH: OH file
  ! -Y turn off OH
  ! -z solar zenith angle cosine
  ! -Z number of user azimuthal angles
  
  ! Key band indices:
  ! bnd_idx = 800, wavelength = 1.0005 microns
  ! bnd_idx = 1593, wavelength = 0.55 microns
  ! bnd_idx = 1603, wavelength = 0.5 microns
  ! bnd_idx = 1229, wavelength = 0.700035 microns
  ! bnd_idx = 1643, wavelength = 0.30075 microns
  ! bnd_idx = 744, wavelength = 1.05988 microns (middle of O2-O2 NIR band)
  
  ! TDDR channels (1995): not done yet
  ! TDDR channels (1994):
  ! channel center = 0.380, bnd_idx = 1627, wavelength = 0.380 microns
  ! channel center = 0.412, bnd_idx = 1621, wavelength = 0.410 microns
  ! channel center = 0.500, bnd_idx = 1603, wavelength = 0.500 microns
  ! channel center = 0.675, bnd_idx = 1282, wavelength = 0.674992 microns
  ! channel center = 0.862, bnd_idx = 961, wavelength = 0.861698 microns
  ! channel center = 1.064, bnd_idx = 740, wavelength = 1.0644 microns
  ! channel center = 1.640, bnd_idx = 410, wavelength = 1.64069 microns
  
  ! Debugging usage:
  ! Perform full setup for all bands but only call DISORT() for #1603:
  ! swnb2 --drc_in=${DATA}/aca -D 1 -E -e 1603 -d ~/foo.nc 
  ! Intercompare model on different machines:
  ! swnb2 --drc_in=${DATA}/aca -d ~/foo.nc;ncks -H -v flx_bb_dwn_sfc ~/foo.nc
  
  ! Production usage:
  ! ncks command to print out diagnostic fields, for many test cases, is
  ! ncks -H -C -v flx_bb_abs_atm,flx_bb_dwn_sfc ${DATA}/aca/swnb.nc
  ! swnb2 --drc_in=${DATA}/aca --drc_out=${DATA}/aca -D 1 -d ~/foo.nc 
  ! Snow:
  ! swnb2 -p ${DATA}/aca/mls_snw_JRI.nc -d ${DATA}/aca/swnb_snw.nc
  ! swnb2 --thermal=false -s 4 -p ${DATA}/aca/mls_snw_JRI.nc -d ~/foo.nc > ~/foo_snw 2>&1 &
  ! swnb2 --thermal=false -s 4 -p ${DATA}/aca/mls_snw_JRI.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ~/foo.nc > ~/foo_snw 2>&1 &
  ! ncks -F -H -C -v '.snw' ~/foo.nc | m

  ! Zenith angles:
  ! swnb2 ${crd_arg} -d ~/foo.nc
  ! crd_arg='--lat=+23.44 --doy=172.5' # Summer solstice NH
  ! crd_arg='--lat=-23.44 --doy=355.5' # Summer solstice SH
  ! crd_arg='--lat=+00.00 --doy=083.5' # Spring equinox  NH
  ! crd_arg='--lat=+00.00 --doy=002.5' # Perihelion
  ! crd_arg='--lat=+00.00 --doy=185.5' # Aphelion
  ! ncks -Q -F -H -C -d bnd,0.5e-6 -v slr_zen_ngl_cos,lcl_yr_day,lat_dgr,xnt_fac ~/foo.nc

  ! Stefan Kinne's test cases:
  ! NB: Thermal emission (flg_Planck) does not work with SK test cases
  ! This is because AFGL temperature gradients are too large at TOA
  ! All gases, no aerosol, cloud, or thermal emission, in MLS, TRP atm with cos theta=0.8660,0.2588 alb=0.2,
  ! swnb2 -z 0.8660 -r 0.2 -s 4 -p ${DATA}/aca/mls_afgl_73lvl.nc -d ${DATA}/tmp/sk_01.nc > ~/foo_1 2>&1 &
  ! swnb2 -z 0.2588 -r 0.2 -s 4 -p ${DATA}/aca/mls_afgl_73lvl.nc -d ${DATA}/tmp/sk_02.nc > ~/foo_2 2>&1 &
  ! swnb2 -z 0.8660 -r 0.2 -s 4 -p ${DATA}/aca/trp_afgl_73lvl.nc -d ${DATA}/tmp/sk_03.nc > ~/foo_3 2>&1 &
  ! swnb2 -z 0.2588 -r 0.2 -s 4 -p ${DATA}/aca/trp_afgl_73lvl.nc -d ${DATA}/tmp/sk_04.nc > ~/foo_4 2>&1 &
  ! swnb2 -z 0.8660 -r 0.2 -s 4 -p ${DATA}/aca/sas_afgl_73lvl.nc -d ${DATA}/tmp/sk_05.nc > ~/foo_5 2>&1 &
  ! ncks -H -C -d bnd,0.7e-6 -v flx_bb_dwn_sfc ${DATA}/tmp/sk_01.nc
  
  ! Test cases: Set CO2 vmr to 300 ppm for most of these!
  ! NB: BPB used only most common isotope of each species for these
  
  ! O2 absorption only:
  ! swnb2 -A -B -C -G -H -I -J -K -L -R -U -W -X -Y -D 1 -r 0.1 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -o ${DATA}/aca/swnb_16O2.nc
  
  ! CO2 absorption only:
  ! swnb2 -A -B -G -H -I -J -K -L -O -R -U -W -X -Y -D 1 -r 0.1 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -c ${DATA}/aca/swnb_12C_16O2.nc
  
  ! H2O absorption only:
  ! swnb2 -A -B -C -G -I -J -K -L -O -R -U -W -X -Y -D 1 -r 0.0 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -h ${DATA}/aca/swnb_1H2_16O.nc
  ! swnb2 -A -B -C -G -I -J -K -L -O -R -U -W -X -Y -D 1 -r 0.0 -z 0.8660254 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -h ${DATA}/aca/swnb_1H2_16O.nc
  ! swnb2 -A -B -C -G -I -J -K -L -O -R -U -W -X -Y -D 1 -r 0.0 -z 0.8660254 -S 1370.0 -p ${DATA}/aca/mls_icrccm_92lvl.nc -h ${DATA}/aca/swnb_H2O.nc
  
  ! O3 absorption only:
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -X -Y -D 1 -r 0.0 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -d ~/swnb.nc 
  
  ! Rayleigh scattering only:
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -U -W -X -Y -r 0.0 -z 0.5 -d ~/swnb.nc
  
  ! Mie scattering/absorption only:
  ! swnb2 -C -G -H -J -K -O -R -U -W -X -Y -r 0.0 -z 0.5 -m 0.1 -d ~/swnb.nc
  
  ! Full physics: 18 layer MLS, cos theta=.5, A=.2, CWP=100
  ! swnb2 -z 0.5 -r 0.2 -m 0.1 -s 4 -d ~/swnb.nc
  ! ncks -C -H -d bnd,1.26e-6 -v flx_bb_dwn_sfc,odxc_spc_O2O2,odxc_spc_O2N2 ${DATA}/aca/swnb.nc
  
  ! Sensitivity studies: Everything except given process, then that process alone
  ! Control:
  ! swnb2 -r 0.1 -z 0.5 -m 0.0 -d ~/swnb.nc
  ! ncks -C -H -v flx_bb_dwn_sfc,flx_bb_abs_atm ~/swnb.nc
  
  ! H2OH2O absorption:
  ! swnb2 -Q -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_H2OH2O.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_H2OH2O.nc
  
  ! O2-O2 absorption:
  ! swnb2 -K -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2O2.nc
  ! swnb2 -A -B -C -G -H -I -J -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2O2.nc
  
  ! O2-N2 absorption:
  ! swnb2 -U -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2N2.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2N2.nc
  
  ! O2 absorption:
  ! swnb2 -O -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O2.nc
  
  ! O3 absorption:
  ! swnb2 -W -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O3.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_O3.nc
  
  ! OH absorption:
  ! swnb2 -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_OH.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -W -X -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_OH.nc
  
  ! CH4 absorption:
  ! swnb2 -G -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_CH4.nc
  ! swnb2 -A -B -C -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_CH4.nc
  
  ! NO2 absorption:
  ! swnb2 -X -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_NO2.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -W -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_NO2.nc
  
  ! CO2 absorption:
  ! swnb2 -C -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_CO2.nc
  ! swnb2 -A -B -G -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_CO2.nc
  
  ! H2O absorption:
  ! swnb2 -H -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_H2O.nc
  ! swnb2 -A -B -C -G -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_H2O.nc
  
  ! Rayleigh scattering:
  ! swnb2 -R -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_Ray.nc
  ! swnb2 -A -B -C -G -H -I -J -K -L -O -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_Ray.nc
  
  ! Thermal absorption/emission:
  ! swnb2 -J -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_Plk.nc
  ! swnb2 -A -B -C -G -H -I -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_Plk.nc
  
  ! Aerosol scattering/absorption:
  ! swnb2 -A -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_aer.nc
  ! swnb2 -B -C -G -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_aer.nc
  
  ! Background aerosol scattering/absorption:
  ! swnb2 -B -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_bga.nc
  ! swnb2 -A -C -G -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ~/swnb_bga.nc
  
  ! Liquid scattering/absorption:
  ! swnb2 -L -r 0.1 -z 0.5 -m 0.01 -f -d ~/swnb_lqd.nc
  ! swnb2 -A -B -C -G -H -I -J -K -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.01 -f -d ~/swnb_lqd.nc
  
  ! Ice scattering/absorption:
  ! swnb2 -I -r 0.1 -z 0.5 -m 0.01 -F -d ~/swnb_ice.nc
  ! swnb2 -A -B -C -G -H -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.01 -F -d ~/swnb_ice.nc
  
  ! ARESE simulations: 
  ! swnb2 -p ${DATA}/arese/clm/951011_1200_arese_clm.nc -d ~/foo.nc
  ! ncks -C -H -v flx_bb_dwn_sfc,flx_bb_abs_atm ~/foo.nc
  ! ncks -C -H -d bnd,0.5e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ~/foo.nc
  ! ncks -C -H -d bnd,0.413e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ~/foo.nc
  ! ncks -C -H -d bnd,0.860e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ~/foo.nc
  
  ! Dust simulations
  ! j_NO2
  ! swnb2 -A -z 0.866 -p ${DATA}/swnb2/arese_951011_1200_clm_clr_aer_dst.nc -d ${DATA}/swnb2/arese_951011_1200_mdl_clr_cln.nc > ~/foo.cln 2>&1 &
  ! swnb2 -z 0.866 -a ${DATA}/aca/aer_mineral_dust.nc -p ${DATA}/swnb2/arese_951011_1200_clm_clr_aer_dst.nc -d ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_dst.nc > ~/foo.dst 2>&1 &
  ! swnb2 -z 0.866 -a ${DATA}/aca/aer_h2so4_300K.nc -p ${DATA}/swnb2/arese_951011_1200_clm_clr_aer_slf.nc -d ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_slf.nc > ~/foo.slf 2>&1 &
  ! ncdiff -v j_NO2,lev,z,bnd,wvl_ctr,flx_spc_act_pht_TOA,flx_spc_act_pht_sfc ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_dst.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_cln.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_dstmcln.nc
  ! ncks -O -v j_NO2,lev,z,bnd,wvl_ctr,flx_spc_act_pht_TOA,flx_spc_act_pht_sfc -u ${DATA}/swnb2/arese_951011_1200_mdl_clr_cln.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_cln.nc
  ! ncks -O -v j_NO2,lev,z,bnd,wvl_ctr,flx_spc_act_pht_TOA,flx_spc_act_pht_sfc -u ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_dst.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_dst.nc 
  ! ncks -O -v j_NO2,lev,z,bnd,wvl_ctr,flx_spc_act_pht_TOA,flx_spc_act_pht_sfc -u ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_slf.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_aer_slf.nc
  ! ncks -O -v j_NO2,lev,z,bnd,wvl_ctr,flx_spc_act_pht_TOA,flx_spc_act_pht_sfc -u ${DATA}/swnb2/arese_951011_1200_mdl_clr_dstmcln.nc ${DATA}/swnb2/arese_951011_1200_mdl_clr_dstmcln.nc
  
  ! New Options for DISORT2
  ! -j use non-Lambertian surface
  ! -v surface spectral properties file
  ! -V BRDF type
  
  ! Code makes heavy use of acronyms and abbreviations in variable names
  ! These are defined in ~/sw/crr/abb.tex and http://dust.ess.uci.edu/doc/abb/abb.pdf
  
  ! Code proceeds as follows:
  ! Section 1: Initialization
  ! Section 2: Main computation loop over spectral bands
  ! Section 3: Post-processing
  
  ! Begin section 1: Initialization
  ! initialization of program defaults
  ! processing of command line switches
  ! netCDF input of atmospheric thermodynamic and cloud profile
  ! netCDF input of gaseous absorption and aerosol scattering and absorption characteristics
  ! override netCDF input data according to command line switches
  ! band-independent DISORT() initialization
  ! band loops for figuring out where to splice input data together
  ! initialization of  band quantities that are level independent
  ! initialization of level quantities that are band independent
  ! a band-level loop for initializing band-level quantities
  
  use clm_mdl,only:slr_crd_Bri92 ! [mdl] Column (CLM) processing
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use flx_slr_mdl ! [mdl] Solar spectral fluxes
  use mmr_mdl,only:wrp_allocate,wrp_deallocate ! [mdl] Memory management
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use phys_cst_mdl ! [mdl] Fundamental and derived physical constants
  use rt_mdl ! [mdl] Radiative transfer utilities
  use sng_mdl ! [mdl] String manipulation
  use tdy_mdl,only:svp_H2O_lqd_PrK78,svp_H2O_ice_PrK78 ! [mdl] Thermodynamics
  use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
  use wvl_mdl ! [mdl] Wavelength grid parameters
  use xtr_mdl ! [mdl] Extrapolation/interpolation handling
  implicit none
  ! Parameters
  ! NB: DISORT uses fortran unit 99 for its own purposes.
  character(*),parameter::sbr_nm='swnb2' ! [sng] Subroutine name
  character(*),parameter::CVS_Date='$Date$' ! [sng] Date string
  character(*),parameter::CVS_Header='$Id$' ! [sng] Full CVS Header
  character(*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(*),parameter::CVS_Name='$HeadURL$' ! [sng] File name string
  character(*),parameter::CVS_Revision='$Revision$' ! [sng] File revision string
  character(*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer,parameter::azi_nbr_max=16 ! arbitrary
  integer,parameter::bnd_nbr_CO2_max=3180 ! 5 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_H2OH2O_max=8192 ! Old 3 cm-1 #s from Chy98
  integer,parameter::bnd_nbr_H2O_max=1590 ! 10 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_NO2_max=750 ! 10 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_O2O2_max=4086 ! Greenblatt et al. (1992) resolution + SPS98 1.27 um band
  integer,parameter::bnd_nbr_O2_max=1590 ! 10 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_O3_max=158 ! WMO85 resolution
  integer,parameter::bnd_nbr_OH_max=1590 ! 10 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_CH4_max=3180 ! 5 cm-1 resolution from 0.56--5.0 um
  integer,parameter::bnd_nbr_aer_max=1690 ! WMO-spliced
  integer,parameter::bnd_nbr_bga_max=480 ! 0.01 um resolution from 0.2--5.0 um
  integer,parameter::bnd_nbr_ice_max=480 ! 0.01 um resolution from 0.2--5.0 um
  integer,parameter::bnd_nbr_lqd_max=480 ! 0.01 um resolution from 0.2--5.0 um
  integer,parameter::bnd_nbr_lmn_max=81 ! CIE photopic luminosity function 5 nm resolution from 380-780 nm
  integer,parameter::bnd_nbr_nst_max=1000 ! new FSBR resolution
  integer,parameter::bnd_nbr_rfl_max=1690 ! WMO-spliced
  integer,parameter::bnd_nbr_mpr_max=1690 ! WMO-spliced
  integer,parameter::bnd_nbr_snw_max=1690 ! WMO-spliced
  integer,parameter::bnd_nbr_max=1690 ! WMO-spliced onto 10 cm-1 resolution from 0.56--5.0 um
  integer,parameter::chn_nbr_max=12 ! [nbr] Arbitrary
  integer,parameter::lev_snw_nbr_max=5 ! CLM snow model resolution
  integer,parameter::lev_nbr_max=110 ! roughly 10 mb resolution from 1010 mb to TOA
  integer,parameter::igbp_nbr=16 ! [nbr] number of IGBP surface types
  integer,parameter::plr_nbr_max=16 ! # user polar angles
  integer,parameter::str_nbr_max=16 ! # computational polar angles
  integer,parameter::sng_lng_dfl_fl=80 ! [nbr] Default filename string length
  integer,parameter::sng_lng_dfl_stt=200 ! [nbr] Default statement string length
  real,parameter::tpt_Malkmus_fit=250.0 ! reference temperature for Malkmus parameters
  real,parameter::mss_val=1.0e36 ! Missing value = missing_value and/or _FillValue
  real,parameter::real_tiny=1.0e-10 ! Tiny value to avoid divide-by-zero errors
  ! integer,parameter::bnd_nbr_nst_max=235 ! FSBR resolution
  
  ! Derived parameters
  integer,parameter::asm_prm_mmn_idx=1 ! [idx] Asymmetry parameter Legendre index
  integer,parameter::mmn_nbr_max=str_nbr_max ! # phase function moments
  integer,parameter::bnd_nbr_brdf_max=bnd_nbr_max ! maximum number of BRDF bands
  integer,parameter::levp_nbr_max=lev_nbr_max+1 ! Interface levels
  integer,parameter::tau_nbr_max=lev_nbr_max+1 ! # user optical depths NB: should equal levp_nbr_max
  
  ! Input Arguments
  ! Input/Output Arguments
  ! Output Arguments
  ! Local workspace
  ! 20100108: Replace declarations like "character src_rfr_sng*200" with initialization
  ! statements so that valgrind thinks strings are initialized when used
  character(sng_lng_dfl_fl)::arg_val=nlc      ! [sng] Command line argument value
  character cmd_ln*500      ! [sng] Command line
  character dsh_key*2       ! [sng] Command line dash and switch
  character(sng_lng_dfl_fl)::drc_in=nlc       ! [sng] Input directory
  character(sng_lng_dfl_fl)::drc_out=nlc      ! [sng] Output directory
  character(sng_lng_dfl_fl)::opt_sng=nlc      ! [sng] Option string
  character aer_sng*100
  character mpr_sng*100
  character snw_sng*100
  character azi_sng*100
  character bga_sng*100
  character char_foo*100
  character(sng_lng_dfl_fl)::fl_CO2=nlc  
  character(sng_lng_dfl_fl)::fl_H2OH2O=nlc
  character(sng_lng_dfl_fl)::fl_H2O=nlc
  character(sng_lng_dfl_fl)::fl_NO2=nlc
  character(sng_lng_dfl_fl)::fl_O2=nlc
  character(sng_lng_dfl_fl)::fl_O3=nlc
  character(sng_lng_dfl_fl)::fl_O2O2=nlc
  character(sng_lng_dfl_fl)::fl_OH=nlc
  character(sng_lng_dfl_fl)::fl_CH4=nlc
  character(sng_lng_dfl_fl)::fl_aer=nlc
  character(sng_lng_dfl_fl)::fl_bga=nlc
  character(sng_lng_dfl_fl)::fl_chn=nlc
  character(sng_lng_dfl_fl)::fl_clm=nlc
  character(sng_lng_dfl_fl)::fl_ice=nlc
  character(sng_lng_dfl_fl)::fl_lqd=nlc
  character(sng_lng_dfl_fl)::fl_lmn=nlc
  character(sng_lng_dfl_fl)::fl_nst=nlc
  character(sng_lng_dfl_fl)::fl_out=nlc
  character(sng_lng_dfl_fl)::fl_rfl=nlc
  character(sng_lng_dfl_fl)::fl_mpr=nlc
  character(sng_lng_dfl_fl)::fl_snw=nlc
  character(sng_lng_dfl_fl)::fl_slr=nlc
  character(sng_lng_dfl_fl)::fl_brdf=nlc
  character*26 lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character opt_dep_sng*100
  character plr_sng*100
  character prf_sng*100
  character(sng_lng_dfl_stt)::prg_ID=nlc
  character(sng_lng_dfl_stt)::src_rfr_sng=nlc 
  character(sng_lng_dfl_stt)::str_sng=nlc
  character(sng_lng_dfl_stt)::stt_CO2=nlc
  character(sng_lng_dfl_stt)::stt_H2O=nlc
  character(sng_lng_dfl_stt)::stt_H2OH2O=nlc
  character(sng_lng_dfl_stt)::stt_Herzberg=nlc
  character(sng_lng_dfl_stt)::stt_NO2=nlc
  character(sng_lng_dfl_stt)::stt_O2=nlc
  character(sng_lng_dfl_stt)::stt_O3=nlc
  character(sng_lng_dfl_stt)::stt_O2O2=nlc
  character(sng_lng_dfl_stt)::stt_O2N2=nlc
  character(sng_lng_dfl_stt)::stt_OH=nlc
  character(sng_lng_dfl_stt)::stt_CH4=nlc
  character(sng_lng_dfl_stt)::stt_Planck=nlc
  character(sng_lng_dfl_stt)::stt_Rayleigh=nlc
  character(sng_lng_dfl_stt)::stt_aer=nlc
  character(sng_lng_dfl_stt)::stt_bga=nlc
  character(sng_lng_dfl_stt)::stt_ice=nlc
  character(sng_lng_dfl_stt)::stt_lqd=nlc
  character(sng_lng_dfl_stt)::stt_slr=nlc
  character(sng_lng_dfl_stt)::stt_flt_lmn=nlc
  character(sng_lng_dfl_stt)::stt_flt_nst=nlc
  character(sng_lng_dfl_stt)::stt_cld_sat=nlc
  character(sng_lng_dfl_stt)::stt_sat_cld=nlc
  character(sng_lng_dfl_stt)::stt_sct_lqd=nlc
  character(sng_lng_dfl_stt)::stt_rfl=nlc
  character(sng_lng_dfl_stt)::stt_mpr=nlc
  character(sng_lng_dfl_stt)::stt_msm=nlc
  character(sng_lng_dfl_stt)::stt_mie=nlc
  character(sng_lng_dfl_stt)::stt_snw=nlc
  character(sng_lng_dfl_stt)::stt_vpr_H2O_abs_cld=nlc
  character(sng_lng_dfl_stt)::stt_top_lvl=nlc
  
  integer arg_idx           ! [idx] Counting index
  integer arg_nbr           ! [nbr] Number of command line arguments
  integer bnd_dbg           ! [idx] Debugging band
  integer exit_status       ! [enm] Program exit status
  integer int_foo           ! [nbr] Integer
  integer opt_lng           ! [nbr] Length of option
  integer rcd               ! [rcd] Return success code
  
  logical cmd_ln_alb
  logical cmd_ln_alb_sfc_NIR_dff
  logical cmd_ln_alb_sfc_NIR_drc
  logical cmd_ln_alb_sfc_vsb_dff
  logical cmd_ln_alb_sfc_vsb_drc
  logical cmd_ln_dns_snw
  logical cmd_ln_dpt_snw
  logical cmd_ln_lat_dgr
  logical cmd_ln_lcl_yr_day
  logical cmd_ln_mmr_mpr_snw
  logical cmd_ln_mpc_CWP
  logical cmd_ln_odxc_obs_aer
  logical cmd_ln_odxc_obs_snw
  logical cmd_ln_slr_cst
  logical cmd_ln_slr_zen_ngl_cos
  logical flg_CO2
  logical flg_H2OH2O
  logical flg_H2O
  logical flg_Herzberg
  logical flg_NO2
  logical flg_O2
  logical flg_O3
  logical flg_O2O2
  logical flg_O2N2
  logical flg_OH
  logical flg_CH4
  logical flg_Planck
  logical flg_Rayleigh
  logical flg_aer ! [flg] Aerosol is optically active
  logical flg_bga ! [flg] Background aerosol is optically active
  logical flg_ice ! [flg] Ice water clouds are optically active
  logical flg_lqd ! [flg] Liquid water clouds are optically active
  logical flg_bfb ! [flg] Produce bit-for-bit output
  logical flg_rfl ! [flg] Incorporate surface reflectance dataset
  logical flg_mie ! [flg] Use Legendre polynomial phase function expansion not HG
  logical flg_mie_aer ! [flg] Require phase function expansion for aerosol
  logical flg_mie_bga ! [flg] Require phase function expansion for background aerosol
  logical flg_mie_ice ! [flg] Require phase function expansion for ice
  logical flg_mie_lqd ! [flg] Require phase function expansion for liquid
  logical flg_mie_mpr ! [flg] Require phase function expansion for snow impurities
  logical flg_mie_snw ! [flg] Require phase function expansion for snow
  logical flg_mpr ! [flg] Snow grain impurities are optically active
  logical flg_msm ! [flg] Use multi-layer snow model if snow structure present
  logical flg_snw ! [flg] Snow grains are optically active
  logical flg_cld_sat ! [flg] Force cloudy and snowy layers to be saturated
  logical flg_sat_cld ! [flg] Force saturated layers to be cloudy
  logical flg_sct_lqd ! [flg] Liquid cloud droplets are pure scatterers
  logical flg_xtr_aer_snw ! [flg] Extrapolate surface aerosol into snowpack
  logical flg_vpr_H2O_abs_cld ! [flg] H2O vapor absorption in cloud allowed
  logical flt_lmn ! [flg] Apply luminosity filter
  logical flt_nst
  logical force_ice_phz
  logical force_lqd_phz
  logical single_bnd_computation
  logical sv_cmp_tau
  logical sv_cmp_plr_ngl
  logical sv_ntn
  logical top_lvl
  logical tst_case_HG
  logical tst_case_Rayleigh
  logical mode_std
  
  integer bnd_nbr           ! dimension size
  integer bnd_nbr_CO2       ! dimension size
  integer bnd_nbr_H2O       ! dimension size
  integer bnd_nbr_H2OH2O    ! dimension size
  integer bnd_nbr_NO2       ! dimension size
  integer bnd_nbr_O2        ! dimension size
  integer bnd_nbr_O2O2      ! dimension size
  integer bnd_nbr_O3        ! dimension size
  integer bnd_nbr_OH        ! dimension size
  integer bnd_nbr_CH4       ! dimension size
  integer bnd_nbr_aer       ! dimension size
  integer bnd_nbr_brdf      ! dimension size BRDF
  integer bnd_nbr_bga       ! dimension size
  integer bnd_nbr_ice       ! dimension size
  integer bnd_nbr_lqd       ! dimension size
  integer bnd_nbr_rfl       ! dimension size
  integer bnd_nbr_mpr       ! dimension size
  integer bnd_nbr_snw       ! dimension size
  integer bnd_nbr_non_O3    ! dimension size
  integer bnd_nbr_lmn       ! dimension size
  integer bnd_nbr_nst       ! dimension size
  integer lev_nbr           ! dimension size
  integer lev_atm_nbr       ! dimension size
  integer lev_snw_nbr       ! dimension size
  integer lev_bnd_mpr_nbr   ! dimension size
  integer lev_bnd_snw_nbr   ! dimension size
  integer levp_snw_nbr      ! dimension size
  integer levp_atm_nbr      ! dimension size
  integer levp_nbr          ! dimension size
  integer azi_nbr           ! dimension size
  integer chn_nbr           ! dimension size
  integer mmn_nbr           ! dimension size
  integer mmn_spc_nbr       ! dimension size
  integer lgn_dmn_sz        ! dimension size
  integer plr_nbr           ! dimension size
  integer bnd_nbr_pure_O3   ! dimension size
  integer str_nbr           ! dimension size
  integer tau_nbr           ! dimension size
#ifdef _OPENMP
  integer omp_get_num_threads ! [nbr] OpenMP number of threads
  integer omp_get_thread_num ! [idx] OpenMP thread index
  integer thr_nbr           ! [nbr] OpenMP number of threads
#endif /* not _OPENMP */
  
  !  integer bnd_idx_brdf	! counting index BRDF
  !integer bnd_idx_O2        ! counting index
  !integer bnd_idx_OH        ! counting index
  integer azi_idx           ! counting index
  integer bnd_idx           ! counting index
  integer bnd_idx_CH4       ! counting index
  integer bnd_idx_CO2       ! counting index
  integer bnd_idx_H2O       ! counting index
  integer bnd_idx_O3        ! counting index
  integer bnd_idx_nst       ! counting index
  integer bnd_idx_tmp_O3    ! counting index
  integer chn_idx           ! counting index
  integer i21               ! counting index
  integer lev_TOA           ! counting index
  integer levp_TOA          ! counting index
  integer lev_bnd_snw_idx   ! counting index
  integer lev_idx           ! counting index
  integer lev_sfc           ! counting index
  integer levp_sfc          ! counting index
  integer lev_snw_idx       ! counting index
  integer levp_idx          ! counting index
  integer levp_snw_idx      ! counting index
  integer mmn_idx           ! counting index
  integer plr_idx           ! counting index
  integer plr_ndr           ! counting index
  integer plr_zen           ! counting index
  integer tau_idx           ! counting index
  
  integer azi_dmn_id        ! dimension ID for azi
  integer bnd_dmn_id        ! dimension ID for bnd
  integer chn_dmn_id        ! dimension ID for chn
  integer grd_dmn_id        ! dimension ID for grid
  integer dmn_id_vec_foo(4)
  integer srt_one(1)
  integer srt_one_one(2)
  integer srt_one_one_one(3)
  integer cnt_bnd(1)
  integer cnt_bndp(1)
  integer cnt_chn(1)
  integer cnt_lev(1)
  integer cnt_levp(1)
  integer cnt_bnd_lev(2)
  integer cnt_mmn_bnd(2)
  integer cnt_mmn_bnd_lev(3)
  integer dim_azi_plr_chn(3)
  integer dim_azi_plr_bnd_levp(4)
  integer dim_azi_plr_levp(3)
  integer dim_bnd_lev(2)
  integer dim_bnd_levp(2)
  integer dim_bnd_lev_snw(2)
  integer dim_bnd_levp_snw(2)
  integer dim_plr_bnd(2)
  integer dim_mmn_bnd_lev(3)
  integer dim_plr_bnd_levp(3)
  integer dim_plr_levp(2)
  integer lev_dmn_id        ! dimension ID for lev
  integer levp_dmn_id       ! dimension ID for levp
  integer lev_bnd_mpr_dmn_id ! dimension ID for lev in impurity file
  integer lev_bnd_snw_dmn_id ! dimension ID for lev in snow file
  integer lev_snw_dmn_id    ! dimension ID for lev_snw
  integer levp_snw_dmn_id   ! dimension ID for levp_snw
  integer mmn_spc_dmn_id    ! dimension ID for mmn
  integer nc_id             ! file handle
  integer mmn_dmn_id        ! dimension ID for mmn
  integer plr_dmn_id        ! dimension ID for plr
  integer tau_dmn_id        ! dimension ID for tau
  integer nf90_r8 ! [enm] External netCDF type for r8 kind
  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call
  
  ! netCDF output variables not contained in input files
  integer abs_bb_SAS_id
  integer abs_bb_atm_id
  integer abs_bb_sfc_id
  integer abs_nst_SAS_id
  integer abs_nst_atm_id
  integer abs_nst_sfc_id
  integer abs_spc_SAS_id
  integer abs_spc_atm_id
  integer abs_spc_sfc_id
  integer alb_sfc_id
  integer azi_dgr_id          
  integer azi_id            ! coordinate ID
  integer bnd_id            ! coordinate ID
  integer flx_abs_atm_rdr_id
  integer flx_bb_abs_atm_id
  integer flx_bb_abs_id
  integer flx_bb_abs_sfc_id
  integer flx_bb_abs_ttl_id
  integer flx_bb_dwn_TOA_id
  integer flx_bb_upw_TOA_id
  integer flx_bb_dwn_dff_id
  integer flx_bb_dwn_drc_id
  integer flx_bb_dwn_dff_sfc_id
  integer flx_bb_dwn_drc_sfc_id
  integer flx_bb_dwn_id
  integer flx_bb_dwn_sfc_id
  integer flx_bb_net_id
  integer flx_bb_upw_id
  integer flx_chn_dwn_TOA_id
  integer flx_chn_upw_TOA_id
  integer flx_frc_dwn_sfc_blr_id
  integer flx_frc_dwn_sfc_id
  integer flx_nst_abs_atm_id
  integer flx_nst_abs_id
  integer flx_nst_abs_sfc_id
  integer flx_nst_abs_ttl_id
  integer flx_nst_dwn_TOA_id
  integer flx_nst_dwn_id
  integer flx_nst_dwn_sfc_id
  integer flx_nst_net_id
  integer flx_nst_upw_id
  integer flx_spc_abs_SAS_id
  integer flx_spc_abs_atm_id
  integer flx_spc_abs_id
  integer flx_spc_abs_sfc_id
  integer flx_spc_act_pht_TOA_id
  integer flx_spc_act_pht_sfc_id
  integer flx_spc_dwn_TOA_id
  integer flx_spc_dwn_dff_id
  integer flx_spc_dwn_drc_id
  integer flx_spc_dwn_dff_sfc_id
  integer flx_spc_dwn_drc_sfc_id
  integer flx_spc_dwn_id
  integer flx_spc_dwn_sfc_id
  integer flx_spc_pht_dwn_sfc_id
  integer flx_spc_upw_id
  integer flx_slr_frc_id
  integer htg_rate_bb_id
  integer ilm_dwn_TOA_id
  integer ilm_dwn_id
  integer ilm_dwn_sfc_id
  integer ilm_upw_id
  integer j_NO2_id
  integer j_spc_NO2_sfc_id
  integer lmn_bb_aa_id
  integer lmn_spc_aa_ndr_id
  integer lmn_spc_aa_ndr_TOA_id
  integer lmn_spc_aa_ndr_sfc_id
  integer lmn_spc_aa_sfc_id
  integer ntn_bb_aa_id
  integer ntn_bb_mean_id
  integer rfl_chn_TOA_id
  integer ntn_spc_aa_ndr_id
  integer ntn_spc_aa_ndr_TOA_id
  integer ntn_spc_aa_ndr_sfc_id
  integer ntn_spc_aa_sfc_id
  integer ntn_spc_aa_zen_id
  integer ntn_spc_aa_zen_sfc_id
  integer ntn_spc_chn_id
  integer ntn_spc_mean_id
  integer odxc_spc_CO2_id
  integer odxc_spc_H2OH2O_id
  integer odxc_spc_H2O_id
  integer odxc_spc_NO2_id
  integer odxc_spc_O2_id
  integer odxc_spc_O3_id
  integer odxc_spc_O2O2_id
  integer odxc_spc_O2N2_id
  integer odxc_spc_OH_id
  integer odxc_spc_CH4_id
  integer odxc_spc_Ray_id
  integer odac_spc_ice_id
  integer odac_spc_lqd_id
  integer odac_spc_aer_id
  integer odac_spc_mpr_id
  integer odac_spc_snw_id
  integer odac_spc_bga_id
  integer odxc_spc_aer_id
  integer odxc_spc_mpr_id
  integer odxc_spc_snw_id
  integer odxc_spc_bga_id
  integer odxc_spc_ice_id
  integer odxc_spc_lqd_id
  integer odxc_spc_ttl_id
  integer nrg_pht_id
  integer plr_cos_id          
  integer plr_dgr_id          
  integer plr_id            ! coordinate ID
  integer rfl_bb_SAS_id
  integer rfl_bb_sfc_id
  integer rfl_nst_SAS_id
  integer rfl_nst_sfc_id
  integer rfl_spc_SAS_id
  integer rfl_ddm_spc_snw_id
  integer trn_ddm_spc_snw_id
  integer rfl_dff_spc_snw_id
  integer trn_dff_spc_snw_id
  integer rfl_drc_spc_snw_id
  integer trn_ttl_drc_spc_snw_id
  integer alb_spc_snw_id
  integer rfl_dff_spc_snw_cnt_id
  integer trn_dff_spc_snw_cnt_id
  integer rfl_drc_upw_snp_id
  integer rfl_dff_upw_snp_id
  integer rfl_dff_dwn_snp_id
  integer trn_ttl_drc_snp_id
  integer trn_drc_drc_snp_id
  integer alb_dff_spc_snw_dea_id
  integer alb_drc_spc_snw_dea_id
  integer flx_spc_abs_snw_id
  integer flx_spc_dwn_snw_id
  integer flx_spc_upw_snw_id
  integer flx_bb_abs_snw_id
  integer flx_bb_dwn_snw_id
  integer rfl_bb_snw_id
  integer abs_bb_snw_id
  integer trn_bb_snw_id
  integer rfl_spc_sfc_id
  integer tau_id            ! coordinate ID
  integer tau_prs_id          
  integer trn_bb_atm_id
  integer trn_nst_atm_id
  integer trn_spc_atm_CO2_id
  integer trn_spc_atm_H2OH2O_id
  integer trn_spc_atm_H2O_id
  integer trn_spc_atm_NO2_id
  integer trn_spc_atm_O2_id
  integer trn_spc_atm_O3_id
  integer trn_spc_atm_O2O2_id
  integer trn_spc_atm_O2N2_id
  integer trn_spc_atm_OH_id
  integer trn_spc_atm_CH4_id
  integer trn_spc_atm_Ray_id
  integer trn_spc_atm_aer_id
  integer trn_spc_atm_mpr_id
  integer trn_spc_atm_snw_id
  integer trn_spc_atm_bga_id
  integer trn_spc_atm_ice_id
  integer trn_spc_atm_lqd_id
  integer trn_spc_atm_ttl_id
  integer wvl_id
  integer wvl_ctr_id
  integer wvl_grd_id
  integer wvl_max_id
  integer wvl_min_id
  integer wvl_dlt_id
  integer wvn_ctr_id
  integer wvn_max_id
  integer wvn_min_id
  integer wvn_dlt_id
  integer alt_id
  integer alt_ntf_id
  
  ! WMO input variables
  integer abs_xsx_O2_id
  integer abs_xsx_O3_id
  integer wvl_max_O3_id
  integer wvl_min_O3_id
  integer wvl_ctr_O3_id
  
  ! O2-O2 input variables
  integer abs_xsx_O2O2_id
  integer wvl_grd_O2O2_id
  
  ! NO2 input variables
  integer abs_xsx_NO2_id
  integer qnt_yld_NO2_id
  integer wvl_grd_NO2_id
  
  ! Reflectance input variables
  integer wvl_grd_rfl_id

  ! H2OH2O input variables
  integer abs_xsx_H2OH2O_id
  integer wvl_grd_H2OH2O_id
  
  ! Narrow band H2O input variables
  integer A_phi_H2O_id
  integer A_psi_H2O_id
  integer B_phi_H2O_id
  integer B_psi_H2O_id
  integer S_d_abs_cff_mss_H2O_id
  integer S_p_abs_cff_mss_H2O_id
  integer wvl_max_H2O_id
  integer wvl_min_H2O_id
  integer wvl_ctr_H2O_id
  
  ! Narrow band CO2 input variables
  integer A_phi_CO2_id
  integer A_psi_CO2_id
  integer B_phi_CO2_id
  integer B_psi_CO2_id
  integer S_d_abs_cff_mss_CO2_id
  integer S_p_abs_cff_mss_CO2_id
  integer wvl_max_CO2_id
  integer wvl_min_CO2_id
  integer wvl_ctr_CO2_id
  
  ! Narrow band OH input variables
  integer A_phi_OH_id
  integer A_psi_OH_id
  integer B_phi_OH_id
  integer B_psi_OH_id
  integer S_d_abs_cff_mss_OH_id
  integer S_p_abs_cff_mss_OH_id
  integer wvl_max_OH_id
  integer wvl_min_OH_id
  integer wvl_ctr_OH_id
  
  ! Narrow band CH4 input variables
  integer A_phi_CH4_id
  integer A_psi_CH4_id
  integer B_phi_CH4_id
  integer B_psi_CH4_id
  integer S_d_abs_cff_mss_CH4_id
  integer S_p_abs_cff_mss_CH4_id
  integer wvl_max_CH4_id
  integer wvl_min_CH4_id
  integer wvl_ctr_CH4_id
  
  ! Narrow band O2 input variables
  integer A_phi_O2_id
  integer A_psi_O2_id
  integer B_phi_O2_id
  integer B_psi_O2_id
  integer S_d_abs_cff_mss_O2_id
  integer S_p_abs_cff_mss_O2_id
  integer wvl_max_O2_id
  integer wvl_min_O2_id
  integer wvl_ctr_O2_id
  
  ! Aerosol input variables
  integer lgn_xpn_cff_aer_id
  integer asm_prm_aer_id
  integer ext_cff_mss_aer_id
  integer abs_cff_mss_aer_id
  integer sca_cff_mss_aer_id
  integer wvl_grd_aer_id
  
  ! Background aerosol input variables
  integer lgn_xpn_cff_bga_id
  integer asm_prm_bga_id
  integer abs_cff_mss_bga_id
  integer sca_cff_mss_bga_id
  integer wvl_grd_bga_id
  
  ! Ice water input variables
  integer lgn_xpn_cff_ice_id
  integer asm_prm_ice_id
  integer abs_cff_mss_ice_id
  integer sca_cff_mss_ice_id
  integer wvl_grd_ice_id
  
  ! Liquid water input variables
  integer lgn_xpn_cff_lqd_id
  integer asm_prm_lqd_id
  integer abs_cff_mss_lqd_id
  integer sca_cff_mss_lqd_id
  integer wvl_grd_lqd_id
  
  ! Impurity input variables
  integer lgn_xpn_cff_mpr_id
  integer asm_prm_mpr_id
  integer abs_cff_mss_mpr_id
  integer sca_cff_mss_mpr_id
  integer wvl_grd_mpr_id

  ! Snow input variables
  integer lgn_xpn_cff_snw_id
  integer asm_prm_snw_id
  integer abs_cff_mss_snw_id
  integer sca_cff_mss_snw_id
  integer wvl_grd_snw_id

  ! Luminosity input variables
  integer lmn_SRF_id
  integer wvl_grd_lmn_id
  
  ! Instrument input variables
  integer chn_SRF_id
  integer nst_SRF_id
  integer wvl_ctr_nst_id
  integer wvl_max_nst_id
  integer wvl_min_nst_id
  
  ! CLM input variables
  integer lev_id            ! coordinate ID
  integer levp_id           ! coordinate ID
  integer lev_snw_id        ! coordinate ID
  integer levp_snw_id       ! coordinate ID
  integer dns_snw_id
  integer dpt_dlt_snw_id
  integer dpt_ntf_snw_id
  integer dpt_snw_id
  integer mmr_mpr_snw_id
  integer rds_ffc_snw_id
  integer tpt_snw_id
  integer tpt_ntf_snw_id
  integer foo_snw_id
  integer lgn_xpn_cff_Mie_ttl_id
  integer RH_lqd_id
  integer alb_sfc_NIR_dff_id
  integer alb_sfc_NIR_drc_id
  integer alb_sfc_vsb_dff_id
  integer alb_sfc_vsb_drc_id
  integer slr_zen_ngl_cos_id
  integer prs_dlt_id
  integer frc_ice_id
  integer frc_ice_ttl_id
  integer grv_id
  integer lat_dgr_id
  integer lcl_time_hr_id
  integer lcl_yr_day_id
  integer lon_dgr_id
  integer mmw_mst_air_id
  integer mpc_CWP_id
  integer mpl_CO2_id
  integer mpl_CWP_id
  integer mpl_H2O_id
  integer mpl_IWP_id
  integer mpl_LWP_id
  integer mpl_O2_id
  integer mpl_OH_id
  integer mpl_CH4_id
  integer mpl_aer_id
  integer mpl_bga_id
  integer mpl_mst_air_id
  integer npl_NO2_id
  integer npl_O2_id
  integer npl_O3_id
  integer npl_O2O2_id
  integer npl_H2OH2O_id
  integer odxc_obs_aer_id
  integer odxc_obs_bga_id
  integer odal_obs_aer_id
  integer odal_obs_bga_id
  integer odsl_obs_aer_id
  integer odsl_obs_bga_id
  integer odxl_obs_aer_id
  integer odxl_obs_bga_id
  integer odxc_obs_mpr_id
  integer odal_obs_mpr_id
  integer odsl_obs_mpr_id
  integer odxl_obs_mpr_id
  integer odxc_obs_snw_id
  integer odal_obs_snw_id
  integer odsl_obs_snw_id
  integer odxl_obs_snw_id
  integer prs_id
  integer prs_ntf_id
  integer q_CO2_id
  integer q_H2O_id
  integer q_O2_id
  integer q_OH_id
  integer q_CH4_id
  integer spc_heat_mst_air_id
  integer tpt_id
  integer tpt_ntf_id
  integer tpt_skn_id
  integer xnt_fac_id
  integer wvl_obs_aer_id
  integer wvl_obs_mpr_id
  integer wvl_obs_snw_id
  integer wvl_obs_bga_id
  integer alt_cld_btm_id
  integer alt_cld_thick_id
  
  ! BRDF input variables
  integer brdf_typ_id
  integer wvl_ctr_brdf_id
  integer wvl_max_brdf_id
  integer wvl_min_brdf_id
  integer f_iso_spc_id
  integer f_vol_spc_id
  integer f_geo_spc_id
  
  integer ocn_msk_id
  integer wnd_spd_id
  real ocn_msk
  
  ! netCDF output variables
  ! Scalars
  real abs_bb_SAS
  real abs_bb_atm
  real abs_bb_sfc
  real abs_nst_SAS
  real abs_nst_atm
  real abs_nst_sfc
  real alb_sfc
  real flx_bb_abs_snw
  real flx_bb_dwn_snw
  real rfl_bb_snw
  real abs_bb_snw
  real trn_bb_snw
  real flx_bb_abs_atm
  real flx_bb_abs_sfc
  real flx_bb_abs_ttl
  real flx_bb_dwn_TOA
  real flx_bb_upw_TOA
  real flx_bb_dwn_sfc
  real flx_bb_dwn_dff_sfc
  real flx_bb_dwn_drc_sfc
  real flx_nst_abs_atm
  real flx_nst_abs_sfc
  real flx_nst_abs_ttl
  real flx_nst_dwn_TOA
  real flx_nst_dwn_sfc
  real ilm_dwn_TOA
  real ilm_dwn_sfc
  real rfl_bb_SAS
  real rfl_bb_sfc
  real rfl_nst_SAS
  real rfl_nst_sfc
  real trn_bb_atm
  real trn_nst_atm

  ! Array dimensions: azi
  real,dimension(:),allocatable::azi !
  real,dimension(:),allocatable::azi_dgr !

  ! Array dimensions: bnd,lev_snw
  real,dimension(:,:),allocatable::rfl_ddm_spc_snw ! [frc] Snow layer spectral flux reflectance (direct+diffuse mean)
  real,dimension(:,:),allocatable::trn_ddm_spc_snw ! [frc] Snow layer spectral flux transmittance (direct+diffuse mean)
  real,dimension(:,:),allocatable::rfl_dff_spc_snw ! [frc] Snow layer spectral flux reflectance, isotropic incidence
  real,dimension(:,:),allocatable::trn_dff_spc_snw ! [frc] Snow layer spectral flux transmittance, isotropic incidence
  real,dimension(:,:),allocatable::rfl_drc_spc_snw ! [frc] Snow layer spectral flux reflectance, direct beam
  real,dimension(:,:),allocatable::trn_ttl_drc_spc_snw ! [frc] Snow layer spectral flux total transmittance to direct beam
  real,dimension(:,:),allocatable::trn_drc_drc ! [frc] Snow layer spectral flux direct beam transmittance

  ! Array dimensions: bnd,lev_bnd_snw
  real,dimension(:,:),allocatable::abs_cff_mss_snw
  real,dimension(:,:),allocatable::sca_cff_mss_snw
  real,dimension(:,:),allocatable::asm_prm_snw
  ! Array dimensions: mmn,bnd,lev_snw
  real,dimension(:,:,:),allocatable::lgn_xpn_cff_snw

  ! Array dimensions: mmn,bnd,lev
  real,dimension(:,:),allocatable::lgn_xpn_cff_aer
  real,dimension(:,:),allocatable::lgn_xpn_cff_bga
  real,dimension(:,:),allocatable::lgn_xpn_cff_ice
  real,dimension(:,:),allocatable::lgn_xpn_cff_lqd
  real,dimension(:,:),allocatable::lgn_xpn_cff_mpr

  ! Array dimensions: bnd
  real,dimension(:),allocatable::ext_cff_mss_aer
  real,dimension(:),allocatable::abs_cff_mss_aer
  real,dimension(:),allocatable::asm_prm_aer
  real,dimension(:),allocatable::sca_cff_mss_aer
  real,dimension(:),allocatable::abs_cff_mss_bga
  real,dimension(:),allocatable::asm_prm_bga
  real,dimension(:),allocatable::sca_cff_mss_bga
  real,dimension(:),allocatable::abs_cff_mss_ice
  real,dimension(:),allocatable::asm_prm_ice
  real,dimension(:),allocatable::sca_cff_mss_ice
  real,dimension(:),allocatable::abs_cff_mss_lqd
  real,dimension(:),allocatable::asm_prm_lqd
  real,dimension(:),allocatable::sca_cff_mss_lqd
  real,dimension(:),allocatable::abs_cff_mss_mpr
  real,dimension(:),allocatable::asm_prm_mpr
  real,dimension(:),allocatable::sca_cff_mss_mpr
  real,dimension(:),allocatable::abs_spc_SAS
  real,dimension(:),allocatable::abs_spc_atm
  real,dimension(:),allocatable::abs_spc_sfc
  real,dimension(:),allocatable::bnd
  real,dimension(:),allocatable::flx_abs_atm_rdr
  real,dimension(:),allocatable::flx_frc_dwn_sfc_blr
  real,dimension(:),allocatable::flx_frc_dwn_sfc
  real,dimension(:),allocatable::flx_slr_frc
  real,dimension(:),allocatable::flx_spc_abs_SAS
  real,dimension(:),allocatable::flx_spc_abs_atm
  real,dimension(:),allocatable::flx_spc_abs_sfc
  real,dimension(:),allocatable::flx_spc_act_pht_TOA
  real,dimension(:),allocatable::flx_spc_act_pht_sfc
  real,dimension(:),allocatable::flx_spc_dwn_TOA
  real,dimension(:),allocatable::flx_spc_dwn_sfc
  real,dimension(:),allocatable::flx_spc_dwn_dff_sfc
  real,dimension(:),allocatable::flx_spc_dwn_drc_sfc
  real,dimension(:),allocatable::flx_spc_pht_dwn_sfc
  real,dimension(:),allocatable::j_spc_NO2_sfc
  real,dimension(:),allocatable::lmn_SRF
  real,dimension(:),allocatable::lmn_spc_aa_ndr_TOA
  real,dimension(:),allocatable::lmn_spc_aa_ndr_sfc
  real,dimension(:),allocatable::nrg_pht
  real,dimension(:),allocatable::ntn_spc_aa_ndr_TOA
  real,dimension(:),allocatable::ntn_spc_aa_ndr_sfc
  real,dimension(:),allocatable::ntn_spc_aa_zen_sfc
  real,dimension(:),allocatable::odac_spc_aer
  real,dimension(:),allocatable::odac_spc_mpr
  real,dimension(:),allocatable::odac_spc_snw
  real,dimension(:),allocatable::odac_spc_bga
  real,dimension(:),allocatable::odac_spc_ice
  real,dimension(:),allocatable::odac_spc_lqd
  real,dimension(:),allocatable::odxc_spc_CO2
  real,dimension(:),allocatable::odxc_spc_H2O
  real,dimension(:),allocatable::odxc_spc_H2OH2O
  real,dimension(:),allocatable::odxc_spc_NO2
  real,dimension(:),allocatable::odxc_spc_O2
  real,dimension(:),allocatable::odxc_spc_O2N2
  real,dimension(:),allocatable::odxc_spc_O2O2
  real,dimension(:),allocatable::odxc_spc_O3
  real,dimension(:),allocatable::odxc_spc_OH
  real,dimension(:),allocatable::odxc_spc_CH4
  real,dimension(:),allocatable::odxc_spc_Ray
  real,dimension(:),allocatable::odxc_spc_aer
  real,dimension(:),allocatable::odxc_spc_mpr
  real,dimension(:),allocatable::odxc_spc_snw
  real,dimension(:),allocatable::odxc_spc_bga
  real,dimension(:),allocatable::odxc_spc_ice
  real,dimension(:),allocatable::odxc_spc_lqd
  real,dimension(:),allocatable::odxc_spc_ttl
  real,dimension(:),allocatable::rfl_spc_SAS
  real,dimension(:),allocatable::alb_spc_snw ! [frc] Snowpack spectral flux reflectance
  real,dimension(:),allocatable::flx_spc_abs_snw ! [frc] Spectral flux absorbed by snowpack
  real,dimension(:),allocatable::flx_spc_dwn_snw ! [frc] Spectral insolation at snowpack
  real,dimension(:),allocatable::flx_spc_upw_snw ! [frc] Spectral upwelling flux at snowpack
  real,dimension(:),allocatable::rfl_spc_sfc
  real,dimension(:),allocatable::trn_spc_atm_CO2
  real,dimension(:),allocatable::trn_spc_atm_H2O
  real,dimension(:),allocatable::trn_spc_atm_H2OH2O
  real,dimension(:),allocatable::trn_spc_atm_NO2
  real,dimension(:),allocatable::trn_spc_atm_O2
  real,dimension(:),allocatable::trn_spc_atm_O2N2
  real,dimension(:),allocatable::trn_spc_atm_O2O2
  real,dimension(:),allocatable::trn_spc_atm_O3
  real,dimension(:),allocatable::trn_spc_atm_OH
  real,dimension(:),allocatable::trn_spc_atm_CH4
  real,dimension(:),allocatable::trn_spc_atm_Ray
  real,dimension(:),allocatable::trn_spc_atm_aer
  real,dimension(:),allocatable::trn_spc_atm_mpr
  real,dimension(:),allocatable::trn_spc_atm_snw
  real,dimension(:),allocatable::trn_spc_atm_bga
  real,dimension(:),allocatable::trn_spc_atm_ice
  real,dimension(:),allocatable::trn_spc_atm_lqd
  real,dimension(:),allocatable::trn_spc_atm_ttl
  real,dimension(:),allocatable::wvl
  real,dimension(:),allocatable::wvl_ctr
  real,dimension(:),allocatable::wvl_dlt
  real,dimension(:),allocatable::wvl_max
  real,dimension(:),allocatable::wvl_min
  real,dimension(:),allocatable::wvn
  real,dimension(:),allocatable::wvn_ctr
  real,dimension(:),allocatable::wvn_dlt
  real,dimension(:),allocatable::wvn_max
  real,dimension(:),allocatable::wvn_min
  
  ! Array dimensions: chn
  real,dimension(:),allocatable::flx_chn_dwn_TOA
  real,dimension(:),allocatable::flx_chn_upw_TOA

  ! Array dimensions: grd
  real,dimension(:),allocatable::wvl_grd

  ! Array dimensions: lev
  real,dimension(:),allocatable::flx_bb_abs !
  real,dimension(:),allocatable::flx_nst_abs !
  real,dimension(:),allocatable::j_NO2 !
  real,dimension(:),allocatable::htg_rate_bb !
  real,dimension(:),allocatable::ntn_bb_mean !

  ! Array dimensions: levp
  real,dimension(:),allocatable::flx_bb_dwn !
  real,dimension(:),allocatable::flx_bb_dwn_dff !
  real,dimension(:),allocatable::flx_bb_dwn_drc !
  real,dimension(:),allocatable::flx_bb_net !
  real,dimension(:),allocatable::flx_bb_upw !
  real,dimension(:),allocatable::flx_nst_dwn !
  real,dimension(:),allocatable::flx_nst_net !
  real,dimension(:),allocatable::flx_nst_upw !
  real,dimension(:),allocatable::ilm_dwn !
  real,dimension(:),allocatable::ilm_upw !

  ! Array dimensions: plr
  real,dimension(:),allocatable::plr
  real,dimension(:),allocatable::plr_cos
  real,dimension(:),allocatable::plr_dgr

  ! Array dimensions: tau
  real,dimension(:),allocatable::tau
  real,dimension(:),allocatable::tau_prs

  ! Array dimensions: mmn,bnd,lev
  real,dimension(:,:,:),allocatable::lgn_xpn_cff_Mie_ttl

  ! Array dimensions: bnd,lev
  real,dimension(:,:),allocatable::asm_prm_HG_ttl
  real,dimension(:,:),allocatable::flx_spc_abs
  real,dimension(:,:),allocatable::ntn_spc_mean
  real,dimension(:,:),allocatable::odxl_spc_ttl
  real,dimension(:,:),allocatable::ss_alb_fct

  ! Array dimensions: bnd,levp
  real,dimension(:,:),allocatable::flx_spc_dwn
  real,dimension(:,:),allocatable::flx_spc_dwn_dff
  real,dimension(:,:),allocatable::flx_spc_dwn_drc
  real,dimension(:,:),allocatable::flx_spc_upw
  real,dimension(:,:),allocatable::lmn_spc_aa_ndr
  real,dimension(:,:),allocatable::ntn_spc_aa_ndr
  real,dimension(:,:),allocatable::ntn_spc_aa_zen

  ! Array dimensions: plr,bnd
  real,dimension(:,:),allocatable::lmn_spc_aa_sfc
  real,dimension(:,:),allocatable::ntn_spc_aa_sfc

  ! Array dimensions: plr,levp
  real,dimension(:,:),allocatable::lmn_bb_aa
  real,dimension(:,:),allocatable::ntn_bb_aa

  ! Array dimensions: azi,plr,bnd
  real,dimension(:,:,:),allocatable::ntn_spc_TOA

  ! Array dimensions: azi,plr,chn
  real,dimension(:,:,:),allocatable::ntn_chn_TOA
  real,dimension(:,:,:),allocatable::rfl_chn_TOA

  ! Array dimensions: azi,plr,levp
  real,dimension(:,:,:),allocatable::ntn_spc_chn

  ! WMO input variables
  real abs_xsx_O2(bnd_nbr_O3_max)
  real abs_xsx_O3(bnd_nbr_O3_max)
  real wvl_max_O3(bnd_nbr_O3_max)
  real wvl_min_O3(bnd_nbr_O3_max)
  real wvl_ctr_O3(bnd_nbr_O3_max)
  
  ! O2-O2 input variables
  real abs_xsx_O2O2_dsk(bnd_nbr_O2O2_max)
  real abs_xsx_O2O2(bnd_nbr_max)
  real wvl_grd_O2O2(bnd_nbr_O2O2_max+1)
  
  ! NO2 input variables
  real abs_xsx_NO2_dsk(bnd_nbr_NO2_max)
  real qnt_yld_NO2_dsk(bnd_nbr_NO2_max)
  real abs_xsx_NO2(bnd_nbr_max)
  real qnt_yld_NO2(bnd_nbr_max)
  real wvl_grd_NO2(bnd_nbr_NO2_max+1)
  
  ! H2OH2O input variables
  real abs_xsx_H2OH2O_dsk(bnd_nbr_H2OH2O_max)
  real abs_xsx_H2OH2O(bnd_nbr_max)
  real wvl_grd_H2OH2O(bnd_nbr_H2OH2O_max+1)
  
  ! Narrow band H2O input variables
  real A_phi_H2O(bnd_nbr_H2O_max)
  real A_psi_H2O(bnd_nbr_H2O_max)
  real B_phi_H2O(bnd_nbr_H2O_max)
  real B_psi_H2O(bnd_nbr_H2O_max)
  real S_d_abs_cff_mss_H2O(bnd_nbr_H2O_max)
  real S_p_abs_cff_mss_H2O(bnd_nbr_H2O_max)
  real wvl_max_H2O(bnd_nbr_H2O_max)
  real wvl_min_H2O(bnd_nbr_H2O_max)
  real wvl_ctr_H2O(bnd_nbr_H2O_max)
  
  ! Narrow band CO2 input variables
  real A_phi_CO2(bnd_nbr_CO2_max)
  real A_psi_CO2(bnd_nbr_CO2_max)
  real B_phi_CO2(bnd_nbr_CO2_max)
  real B_psi_CO2(bnd_nbr_CO2_max)
  real S_d_abs_cff_mss_CO2(bnd_nbr_CO2_max)
  real S_p_abs_cff_mss_CO2(bnd_nbr_CO2_max)
  real wvl_max_CO2(bnd_nbr_CO2_max)
  real wvl_min_CO2(bnd_nbr_CO2_max)
  real wvl_ctr_CO2(bnd_nbr_CO2_max)
  
  ! Narrow band OH input variables
  real A_phi_OH(bnd_nbr_OH_max)
  real A_psi_OH(bnd_nbr_OH_max)
  real B_phi_OH(bnd_nbr_OH_max)
  real B_psi_OH(bnd_nbr_OH_max)
  real S_d_abs_cff_mss_OH(bnd_nbr_OH_max)
  real S_p_abs_cff_mss_OH(bnd_nbr_OH_max)
  real wvl_max_OH(bnd_nbr_OH_max)
  real wvl_min_OH(bnd_nbr_OH_max)
  real wvl_ctr_OH(bnd_nbr_OH_max)
  
  ! Narrow band CH4 input variables
  real A_phi_CH4(bnd_nbr_CH4_max)
  real A_psi_CH4(bnd_nbr_CH4_max)
  real B_phi_CH4(bnd_nbr_CH4_max)
  real B_psi_CH4(bnd_nbr_CH4_max)
  real S_d_abs_cff_mss_CH4(bnd_nbr_CH4_max)
  real S_p_abs_cff_mss_CH4(bnd_nbr_CH4_max)
  real wvl_max_CH4(bnd_nbr_CH4_max)
  real wvl_min_CH4(bnd_nbr_CH4_max)
  real wvl_ctr_CH4(bnd_nbr_CH4_max)
  
  ! Narrow band O2 input variables
  real A_phi_O2(bnd_nbr_O2_max)
  real A_psi_O2(bnd_nbr_O2_max)
  real B_phi_O2(bnd_nbr_O2_max)
  real B_psi_O2(bnd_nbr_O2_max)
  real S_d_abs_cff_mss_O2(bnd_nbr_O2_max)
  real S_p_abs_cff_mss_O2(bnd_nbr_O2_max)
  real wvl_max_O2(bnd_nbr_O2_max)
  real wvl_min_O2(bnd_nbr_O2_max)
  real wvl_ctr_O2(bnd_nbr_O2_max)
  
  ! Aerosol input variables
  real,dimension(:),allocatable::abs_cff_mss_aer_dsk
  real,dimension(:),allocatable::asm_prm_aer_dsk
  real,dimension(:),allocatable::ext_cff_mss_aer_dsk
  real,dimension(:),allocatable::sca_cff_mss_aer_dsk
  real,dimension(:),allocatable::wvl_grd_aer
  real,dimension(:,:),allocatable::lgn_xpn_cff_aer_dsk
  
  ! Background aerosol input variables
  real,dimension(:),allocatable::abs_cff_mss_bga_dsk
  real,dimension(:),allocatable::asm_prm_bga_dsk
  real,dimension(:),allocatable::sca_cff_mss_bga_dsk
  real,dimension(:),allocatable::wvl_grd_bga
  real,dimension(:,:),allocatable::lgn_xpn_cff_bga_dsk
  
  ! Ice water input variables
  real,dimension(:),allocatable::abs_cff_mss_ice_dsk
  real,dimension(:),allocatable::asm_prm_ice_dsk
  real,dimension(:),allocatable::sca_cff_mss_ice_dsk
  real,dimension(:),allocatable::wvl_grd_ice
  real,dimension(:,:),allocatable::lgn_xpn_cff_ice_dsk
  
  ! Liquid water input variables
  real,dimension(:),allocatable::abs_cff_mss_lqd_dsk
  real,dimension(:),allocatable::asm_prm_lqd_dsk
  real,dimension(:),allocatable::sca_cff_mss_lqd_dsk
  real,dimension(:),allocatable::wvl_grd_lqd
  real,dimension(:,:),allocatable::lgn_xpn_cff_lqd_dsk
  
  ! Impurity input variables
  real,dimension(:),allocatable::abs_cff_mss_mpr_dsk
  real,dimension(:),allocatable::asm_prm_mpr_dsk
  real,dimension(:),allocatable::sca_cff_mss_mpr_dsk
  real,dimension(:),allocatable::wvl_grd_mpr
  real,dimension(:,:),allocatable::lgn_xpn_cff_mpr_dsk
  
  ! Snow input variables
  real,dimension(:,:),allocatable::abs_cff_mss_snw_dsk
  real,dimension(:,:),allocatable::asm_prm_snw_dsk
  real,dimension(:,:),allocatable::sca_cff_mss_snw_dsk
  real,dimension(:),allocatable::wvl_grd_snw
  real,dimension(:,:,:),allocatable::lgn_xpn_cff_snw_dsk
  
  ! Reflectance input variables
  real,dimension(:),allocatable::rfl_spc_sfc_dsk
  real,dimension(:),allocatable::wvl_grd_rfl
  
  ! Luminosity input variables
  real,dimension(:),allocatable::lmn_SRF_dsk
  real,dimension(:),allocatable::wvl_grd_lmn
  
  ! Instrument input variables
  integer chn_SRF_msk(bnd_nbr_max)
  real chn_SRF(bnd_nbr_max,chn_nbr_max)
  real nst_SRF(bnd_nbr_nst_max)
  real wvl_ctr_nst(bnd_nbr_nst_max)
  real wvl_max_nst(bnd_nbr_nst_max)
  real wvl_min_nst(bnd_nbr_nst_max)
  
  ! CLM input variables
  real(selected_real_kind(p=12))::lat      ! [rdn] Latitude
  real(selected_real_kind(p=12))::lat_dgr
  real(selected_real_kind(p=12))::lcl_time_hr
  real(selected_real_kind(p=12))::lcl_yr_day
  real(selected_real_kind(p=12))::lon_dgr
  real(selected_real_kind(p=12))::slr_zen_ngl_cos
  
  ! Array dimensions: lev
  real,dimension(:),allocatable::RH_lqd ! [frc] Relative humidity w/r/t liquid
  real,dimension(:),allocatable::lev
  real,dimension(:),allocatable::alt
  real,dimension(:),allocatable::odal_obs_aer
  real,dimension(:),allocatable::odal_obs_mpr
  real,dimension(:),allocatable::odal_obs_snw
  real,dimension(:),allocatable::odal_obs_bga
  real,dimension(:),allocatable::odsl_obs_aer
  real,dimension(:),allocatable::odsl_obs_mpr
  real,dimension(:),allocatable::odsl_obs_snw
  real,dimension(:),allocatable::odsl_obs_bga
  real,dimension(:),allocatable::odxl_obs_aer
  real,dimension(:),allocatable::odxl_obs_mpr
  real,dimension(:),allocatable::odxl_obs_snw
  real,dimension(:),allocatable::odxl_obs_bga
  real,dimension(:),allocatable::tpt
  ! Array dimensions: levp
  real,dimension(:),allocatable::alt_ntf
  real,dimension(:),allocatable::levp
  real,dimension(:),allocatable::tpt_ntf
  ! Array dimensions: lev_snw
  real,dimension(:),allocatable::lev_snw ! [m] Snow depth
  real,dimension(:),allocatable::dns_snw ! [kg m-3] Snow density
  real,dimension(:),allocatable::dpt_dlt_snw ! [m] Snow layer thickness
  real,dimension(:),allocatable::dpt_snw ! [m] Snow depth
  real,dimension(:),allocatable::mmr_mpr_snw ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow
  real,dimension(:),allocatable::rds_ffc_snw ! [m] Snow effective radius
  real,dimension(:),allocatable::tpt_snw ! [K] Snow temperature
  real,dimension(:),allocatable::tpt_ntf_snw ! [K] Snow temperature interfaces
  real,dimension(:),allocatable::foo_snw ! [K] Snow foo
  ! Array dimensions: lev_snw
  real,dimension(:),allocatable::levp_snw ! [m] Snow depth interfaces
  real,dimension(:),allocatable::dpt_ntf_snw ! [m] Snow depth interfaces

  real alb_sfc_NIR_dff
  real alb_sfc_NIR_drc
  real alb_sfc_vsb_dff
  real alb_sfc_vsb_drc
  real alt_cld_btm
  real alt_cld_thick
  real azi_nbr_wvl_dlt
  real frc_ice(lev_nbr_max)
  real frc_ice_ttl
  real grv(lev_nbr_max)
  real mmw_mst_air(lev_nbr_max)
  real mpc_CWP
  real mpl_CO2(lev_nbr_max)
  real mpl_CWP(lev_nbr_max)
  real mpl_H2O(lev_nbr_max)
  real mpl_IWP(lev_nbr_max)
  real mpl_LWP(lev_nbr_max)
  real mpl_O2(lev_nbr_max)
  real mpl_OH(lev_nbr_max)
  real mpl_CH4(lev_nbr_max)
  real mpl_aer(lev_nbr_max)
  real mpl_bga(lev_nbr_max)
  real mpl_mst_air(lev_nbr_max)
  real npl_H2OH2O(lev_nbr_max)
  real npl_NO2(lev_nbr_max)
  real npl_O2(lev_nbr_max)
  real npl_O2O2(lev_nbr_max)
  real npl_O3(lev_nbr_max)
  real odxc_obs_aer
  real odxc_obs_mpr
  real odxc_obs_snw
  real odxc_obs_bga
  real prs(lev_nbr_max)
  real prs_dlt(lev_nbr_max)
  real prs_ntf(levp_nbr_max)
  real q_CO2(lev_nbr_max)
  real q_H2O(lev_nbr_max)
  real q_O2(lev_nbr_max)
  real q_OH(lev_nbr_max)
  real q_CH4(lev_nbr_max)
  real spc_heat_mst_air(lev_nbr_max)
  real tpt_skn
  real wvl_obs_aer
  real wvl_obs_mpr
  real wvl_obs_snw
  real wvl_obs_bga
  real xnt_fac
  
  ! BRDF input variables
  real wvl_max_brdf(bnd_nbr_brdf_max)
  real wvl_min_brdf(bnd_nbr_brdf_max)
  real wvl_ctr_brdf(bnd_nbr_brdf_max)
  real f_iso_spc(bnd_nbr_brdf_max,igbp_nbr)
  real f_vol_spc(bnd_nbr_brdf_max,igbp_nbr)
  real f_geo_spc(bnd_nbr_brdf_max,igbp_nbr)
  ! Cox and Munk spectral parameters
  !real idx_rfr_sfc_spc(bnd_nbr_brdf_max) ! surface index of refraction
  !real nrm_cff_CM_spc(bnd_nbr_brdf_max) ! normalization coefficient
  
  ! Surface bidirectional reflectivity BRDF parameters
  ! Surface type
  common /brdf_com/ brdf_typ
  integer brdf_typ
  ! EOS MODIS BRDF/Albedo Product parameters
  common /brdf_EOS/ f_iso, f_vol, f_geo
  real f_iso                ! linear coefficient of isotropic kernel
  real f_vol                ! linear coefficient of RossThick kernel
  real f_geo                ! linear coefficient of LiSparse kernel
  ! Cox and Munk parameters
  common /brdf_CoxMunk_prm/ nrm_cff_CM, wnd_spd, idx_rfr_sfc
  real nrm_cff_CM           ! empirical normalization constant
  real wnd_spd              ! suface windspeed m / s
  real idx_rfr_sfc          ! surface index of refraction
  ! Hapke parameters
  common /brdf_Hapke_prm/ b0, hh, w
  real b0                   ! empirical factor to account for the finite size of particles
  real hh                   ! angular width parameter of opposition effect
  real w                    ! single scattering albedo
  ! Minnaert parameters
  common /brdf_Minnaert_prm/ nrm_rfl_M, k_cff_M
  real nrm_rfl_M            ! normal reflectance
  real k_cff_M              ! coefficient for power dependence on mu, mup
  ! RahmanPinty parameters
  common /brdf_RahmanPinty_prm/ nrm_cff_RP, k_cff_RP, g_phs
  real nrm_cff_RP           ! reflectance coefficient
  real k_cff_RP             ! coefficient for power dependence on mu, mup
  real g_phs                ! asymmetry parameter in Henyey Greenstein function
  ! LommelSeeliger parameters
  common /brdf_LommelSeeliger_prm/ nrm_rfl_ls
  real nrm_rfl_LS           ! normal reflectance
  ! LiSparse parameters
  common /brdf_LiSparse_prm/ hb, br
  real hb                   ! height to vertical radius ratio
  real br                   ! vertical radius to horizontal radius ratio
  
  ! Local arrays
  real(selected_real_kind(p=12))::lat_dgr_cmd_ln
  real(selected_real_kind(p=12))::lcl_yr_day_cmd_ln
  real(selected_real_kind(p=12))::slr_zen_ngl_cos_cmd_ln
  real(selected_real_kind(p=12))::pi

  ! Allocatable local arrays
  real,dimension(:),allocatable::mpl_mpr ! [kg m-2] Mass path of impurity in layer
  real,dimension(:),allocatable::mpl_snw ! [kg m-2] Mass path of snow in layer
  real,dimension(:),allocatable::frc_sfc_snw ! [frc] Fraction of lowest layer air mass within snow layer
  real,dimension(:),allocatable::odal_mpr ! [frc] Layer impurity absorption optical depth
  real,dimension(:),allocatable::odsl_mpr ! [frc] Layer impurity scattering optical depth
  real,dimension(:),allocatable::odal_snw ! [frc] Layer snow absorption optical depth
  real,dimension(:),allocatable::odsl_snw ! [frc] Layer snow scattering optical depth
  real,dimension(:),allocatable::odal_aer ! [frc] Layer aerosol absorption optical depth
  real,dimension(:),allocatable::odsl_aer ! [frc] Layer aerosol scattering optical depth
  real,dimension(:),allocatable::odal_bga ! [frc] Layer background aerosol absorption optical depth
  real,dimension(:),allocatable::odsl_bga ! [frc] Layer background aerosol scattering optical depth
  real,dimension(:),allocatable::odal_ice ! [frc] Layer ice absorption optical depth
  real,dimension(:),allocatable::odsl_ice ! [frc] Layer ice scattering optical depth
  real,dimension(:),allocatable::odal_lqd ! [frc] Layer liquid absorption optical depth
  real,dimension(:),allocatable::odsl_lqd ! [frc] Layer liquid scattering optical depth
  ! Diagnostic snow reflectance contributions from delta-Eddington/adding method
  real,dimension(:),allocatable::alb_dff_spc_snw_dea ! [frc] Snowpack spectral albedo to isotropic illumination from delta-Eddington/adding
  real,dimension(:),allocatable::alb_drc_spc_snw_dea ! [frc] Snowpack spectral albedo to direct beam from delta-Eddington/adding
  real,dimension(:,:),allocatable::rfl_dff_spc_snw_cnt ! [frc] Layer contribution to snowpack spectral albedo
  real,dimension(:,:),allocatable::trn_dff_spc_snw_cnt ! [frc] Layer contribution to snowpack spectral transmittance
  real,dimension(:,:),allocatable::rfl_drc_upw_snp ! [frc] Spectral reflectance of snowpack beneath interface to downwelling direct beam
  real,dimension(:,:),allocatable::rfl_dff_upw_snp ! [frc] Spectral reflectance of snowpack beneath interface to downwelling diffuse radiation
  real,dimension(:,:),allocatable::rfl_dff_dwn_snp ! [frc] Spectral reflectance of snowpack above interface to upwelling diffuse radiation
  real,dimension(:,:),allocatable::trn_ttl_drc_snp ! [frc] Total transmittance of snowpack above interface to downwelling direct beam
  real,dimension(:,:),allocatable::trn_drc_drc_snp ! [frc] Direct beam (scaled) transmittance of snowpack above interface to downwelling direct beam

  integer aer_lvl_nbr
  integer bnd_obs_aer
  integer bnd_obs_mpr
  integer bnd_obs_snw
  integer bnd_obs_bga
  integer cld_lvl_nbr
  integer mpr_lvl_nbr
  integer slr_spc_xtr_typ
  integer xtr_typ_LHS
  integer xtr_typ_RHS
  
  real xpn_arg ! [frc] Bri92 limit
  real xpn_val ! [frc] Bri92 limit
  real mu_CCY83 ! [frc] Cosine solar zenith angle D/E CCY83
  real alpha_CCY83 ! [frc] Intermediate factor D/E CCY83
  real gamma_CCY83 ! [frc] Intermediate factor D/E CCY83
  real tau_ext ! [frc] Extinction optical depth D/E CCY83
  real sng_sct ! [frc] Single scattering albedo D/E CCY83
  real asm_prm ! [frc] Asymmetry parameter D/E CCY83
  real fwd_frc ! [frc] Forward scattering fraction D/E CCY83
  real tau_ext_scl ! [frc] Scaled extinction optical depth D/E CCY83
  real sng_sct_scl ! [frc] Scaled single scattering albedo D/E CCY83
  real asm_prm_scl ! [frc] Scaled asymmetry parameter D/E CCY83
  real lambda_CCY83 ! [frc] Intermediate factor D/E CCY83
  real u_CCY83 ! [frc] Intermediate factor D/E CCY83
  real N_CCY83 ! [frc] Intermediate factor D/E CCY83
  real tmp_fct_mtx_nvr ! [frc] Temporary factor: Matrix inversion (1-R2*R1)^{-1}
  real rfl_dff_spc_snw_crr ! [frc] Snowpack reflectance through current layer
  real rfl_dff_spc_snw_nxt ! [frc] Snowpack reflectance through next layer
  real trn_dff_spc_snw_crr ! [frc] Snowpack transmittance through current layer
  real trn_dff_spc_snw_nxt ! [frc] Snowpack transmittance through next layer
  real rfl_drc_spc_snw_crr ! [frc] Snowpack reflectance through current layer
  real rfl_drc_spc_snw_nxt ! [frc] Snowpack reflectance through next layer
  real trn_ttl_drc_spc_snw_crr ! [frc] Snowpack transmittance through current layer
  real trn_ttl_drc_spc_snw_nxt ! [frc] Snowpack transmittance through next layer
  real trn_drc_drc_crr ! [frc] Direct beam transmittance through current layer
  real trn_drc_drc_nxt ! [frc] Direct beam transmittance through next layer
  real dpt_snw_crc_fct ! [frc] Snow depth correction factor
  real tpt_vrt_sfc ! [K] Virtual temperature at surface
  real dns_mst_air_sfc ! [kg m-3] Density of moist air at surface
  real sat_vpr_lqd ! [Pa] Saturation vapor pressure w/r/t liquid H2O
  real sat_vpr_ice ! [Pa] Saturation vapor pressure w/r/t ice H2O
  real qst_H2O_lqd ! [kg kg-1] Saturation mixing ratio w/r/t liquid H2O
  real qst_H2O_ice ! [kg kg-1] Saturation mixing ratio w/r/t ice H2O
  real odxl_tmp ! [frc] Temporary optical depth
  real mpl_tmp ! [kg m-2] Temporary optical depth
  real flx_frc_drc_TOA ! [frc] TOA insolation fraction in direct beam
  real alb_cmd_ln
  real bnd_wgt(bnd_nbr_max)
  real bnd_wgt_lmn(bnd_nbr_max)
  real float_foo
  real flx_spc_act
  real flx_spc_act_pht
  real flx_spc_net(bnd_nbr_max,levp_nbr_max)
  real idx_rfr_air_STP(bnd_nbr_max)
  real j_spc_NO2
  real mpc_CWP_cmd_ln
  real dns_snw_cmd_ln ! [kg m-3] Snow density
  real dpt_snw_cmd_ln ! [m] Snowpack thickness
  real mmr_mpr_snw_cmd_ln
  real mpc_IWP
  real odal_CO2(lev_nbr_max)
  real odal_H2O(lev_nbr_max)
  real odal_H2OH2O(lev_nbr_max)
  real odal_NO2(lev_nbr_max)
  real odal_O2(lev_nbr_max)
  real odal_O2N2(lev_nbr_max)
  real odal_O2O2(lev_nbr_max)
  real odal_O3(lev_nbr_max)
  real odal_OH(lev_nbr_max)
  real odal_CH4(lev_nbr_max)
  real odal_spc_ttl(bnd_nbr_max,lev_nbr_max)
  real odsl_tmp
  real odsl_Ray(lev_nbr_max)
  real odsl_spc_ttl(bnd_nbr_max,lev_nbr_max)
  real odsl_HG(lev_nbr_max)
  real odsl_Mie(lev_nbr_max)
  real odxc_obs_aer_cmd_ln
  real odxc_obs_snw_cmd_ln
  real opt_dep_ITOD_CO2_hires(levp_nbr_max,2)
  real opt_dep_ITOD_CH4_hires(levp_nbr_max,2)
  real opt_dep_ITOD_H2O(levp_nbr_max)
  real opt_dep_ITOD_O2(levp_nbr_max)
  real opt_dep_ITOD_OH(levp_nbr_max)
  real opt_dep_LTOD_CO2_hires
  real opt_dep_LTOD_CH4_hires
  real phi_wgt(lev_nbr_max)
  real prs_bar(levp_nbr_max)
  real psi_wgt(lev_nbr_max)
  real sca_cff_mss_Ray(lev_nbr_max)
  real sca_frc_Ray(lev_nbr_max) ! [frc] Scattering fraction treated with Rayleigh phase function
  real sca_frc_HG(lev_nbr_max) ! [frc] Scattering fraction treated with HG phase function
  real sca_frc_Mie(lev_nbr_max) ! [frc] Scattering fraction treated with Mie phase function
  real slr_cst
  real slr_cst_cmd_ln
  real slr_cst_xnt_fac
  real tpt_dlt_Mlk(lev_nbr_max)
  real tpt_dlt_Mlk_sqr(lev_nbr_max)
  real trn_ALT_CO2
  real trn_LT_CO2_hires(lev_nbr_max,2)
  real trn_ALT_CH4
  real trn_LT_CH4_hires(lev_nbr_max,2)
  real u_bar(levp_nbr_max)
  real wvl_Planck
  
  ! The following line is from the DISORT() subroutine: 
  !  PARAMETER ( MXCLY = 92, MXULV = 93, MXCMU = 16, MXUMU = 16, &
  !           MXPHI = 16
  ! parameter ( mxcly = 85, mxulv = 85, mxcmu = 49, mxumu = 16, &
  !      maxpphi = 3)
  ! Parameters declared below should match the values from DISORT():
  integer maxcly            ! [nbr] Maximum number of computational layers
  integer maxcmu            ! [nbr] Maximum number of computational polar angles
  integer maxphi            ! [nbr] Maximum number of output azimuthal angles
  integer maxulv            ! [nbr] Maximum number of output layers
  integer maxumu            ! [nbr] Maximum number of output polar angles
  integer maxmom            ! [nbr] Maximum number of moments? for DISORT2
  parameter(maxcly=tau_nbr_max, &
       maxcmu = str_nbr_max, &
       maxphi = azi_nbr_max, &
       maxulv = tau_nbr_max, &
       maxumu = plr_nbr_max, &
       maxmom = str_nbr_max) ! for DISORT2
  
  ! DISORT() input variables:
  character  header*127
  logical  lamber, plank, onlyfl, prnt(5), usrang, usrtau
  integer  ibcnd, nlyr, numu, nstr, nphi, ntau
  integer nmom              ! for DISORT2
  real accur, albedo, btemp, dtauc( maxcly ), fbeam, fisot, &
       ! swnb2 does not use hl()
       !       hl( 0:maxcmu ), phi( maxphi ), pmom( 0:maxmom, maxcly ), & ! DISORT2
       phi( maxphi ), pmom( 0:maxmom, maxcly ), & ! DISORT2
       phi0, ssalb( maxcly ), temper( 0:maxcly ), temis, ttemp, &
       wvnmlo, wvnmhi, umu( maxumu ), umu0, utau( maxulv )
  
  ! DISORT() output variables:
  real     rfldir( maxulv ), rfldn( maxulv ), flup( maxulv ), &
       dfdt( maxulv ), uavg( maxulv ), u0u( maxumu, maxulv ), &
       uu( maxumu, maxulv, maxphi ), albmed( maxumu ), &
       trnmed( maxumu )
  
  ! Main code
  dbg_lvl=0                 ! Causes DDD source window to display this file
  
  ! Initialize default values
  drc_in='/data/zender/aca'//nlc ! [sng] Input directory
  drc_out=''                ! [sng] Output directory
  fl_CO2='swnb_CO2.nc'//nlc
  fl_H2OH2O='abs_xsx_H2OH2O.nc'//nlc
  fl_H2O='swnb_H2O.nc'//nlc
  fl_OH='swnb_OH.nc'//nlc
  fl_CH4='mlk_CH4.nc'//nlc
  fl_O2='swnb_O2.nc'//nlc
  fl_O3='abs_xsx_O3.nc'//nlc
  fl_O2O2='abs_xsx_O2O2.nc'//nlc
  fl_NO2='abs_xsx_NO2.nc'//nlc
  fl_clm='mls_clr.nc'//nlc
  fl_ice='aer_h2o_ice_rds_swa_20.nc'//nlc
  fl_lqd='aer_h2o_lqd_rds_swa_10.nc'//nlc
  fl_lmn='lmn_CIE.nc'//nlc
  fl_nst='nst_FSBR.nc'//nlc
  fl_chn='epic.nc'//nlc
  fl_aer='aer_sulfate.nc'//nlc
  fl_bga='aer_sulfate.nc'//nlc
  fl_rfl='rfl_spc_sfc.nc'//nlc
  fl_mpr='aer_lac_phb_ctr_FZR07.nc'//nlc
  fl_snw='aer_snw_rds_ffc_100um.nc'//nlc
  fl_out='swnb.nc'//nlc
  fl_slr='spc_Kur95_01wvn.nc'//nlc
  fl_brdf='brdf.nc'//nlc  ! BRDF file
  
  azi_nbr=azi_nbr_max
  bnd_dbg=1603
  brdf_typ=1      ! Default BRDF
  cmd_ln_alb=.false.
  cmd_ln_alb_sfc_NIR_dff=.false.
  cmd_ln_alb_sfc_NIR_drc=.false.
  cmd_ln_alb_sfc_vsb_dff=.false.
  cmd_ln_alb_sfc_vsb_drc=.false.
  cmd_ln_dns_snw=.false.
  cmd_ln_dpt_snw=.false.
  cmd_ln_lat_dgr=.false.
  cmd_ln_lcl_yr_day=.false.
  cmd_ln_mmr_mpr_snw=.false.
  cmd_ln_mpc_CWP=.false.
  cmd_ln_odxc_obs_aer=.false.
  cmd_ln_odxc_obs_snw=.false.
  cmd_ln_slr_cst=.false.
  cmd_ln_slr_zen_ngl_cos=.false.
  exit_status=0             ! [enm] Program exit status
  flg_CH4=.true.
  flg_CO2=.true.
  flg_H2O=.true.
  flg_H2OH2O=.false. ! [flg] H2OH2O is the only gas turned off by default
  flg_Herzberg=.true.
  flg_NO2=.true.
  flg_O2=.true.
  flg_O2N2=.true.
  flg_O2O2=.true.
  flg_O3=.true.
  flg_OH=.true.
  flg_Planck=.true.
  flg_Rayleigh=.true.
  flg_aer=.true.
  flg_bfb=.false. ! [flg] Produce bit-for-bit output
  flg_bga=.true.
  flg_cld_sat=.false. ! [flg] Force cloudy and snowy layers to be saturated
  flg_ice=.true.
  flg_lqd=.true.
  flg_mie=.false. ! [flg] Use phase function expansion where available
  flg_mie_aer=.false. ! [flg] Require phase function expansion for aerosol
  flg_mie_bga=.false. ! [flg] Require phase function expansion for background
  flg_mie_ice=.false. ! [flg] Require phase function expansion for ice
  flg_mie_lqd=.false. ! [flg] Require phase function expansion for liquid
  flg_mie_mpr=.false. ! [flg] Require phase function expansion for snow impurities
  flg_mie_snw=.false. ! [flg] Require phase function expansion for snow
  flg_mpr=.true. ! [flg] Incorporate snowpack impurities (assumes flg_msm)
  flg_msm=.true. ! [flg] Use multi-layer snow model if snow structure present
  flg_rfl=.false. ! [flg] Incorporate surface reflectance dataset
  flg_sat_cld=.false. ! [flg] Force saturated layers to be cloudy
  flg_sct_lqd=.false. ! [flg] Liquid cloud droplets are pure scatterers
  flg_snw=.true. ! [flg] Snow grains are optically active in snow model
  ! flg_toa_isot=.false. ! [flg] TOA isotropic flux as sky brightness in nL?
  flg_vpr_H2O_abs_cld=.true. ! [flg] H2O vapor absorption in cloud allowed
  flg_xtr_aer_snw=.false. ! [flg] Extrapolate surface aerosol into snowpack
  float_foo=0.0
  flt_lmn=.true.
  flt_nst=.true.
  flx_frc_drc_TOA=1.0 ! [frc] TOA insolation fraction in direct beam
  force_ice_phz=.false.
  force_lqd_phz=.false.
  lamber=.true.
  lat_dgr_cmd_ln=mss_val  ! [dgr] Latitude
  lcl_yr_day_cmd_ln=mss_val ! [day] Local year day
  mode_std=.true. ! Standard run mode (non-Fillmore, i.e., no ocean mask or wind speeds or channels)
  odxc_obs_mpr=0.0 ! [frc] Column impurity extinction optical depth 
  odxc_obs_snw=0.0 ! [frc] Column snow extinction optical depth 
  pi=4.0*atan(1.0)
  plr_nbr=2
  rcd=nf90_noerr              ! nf90_noerr == 0
  single_bnd_computation=.false.
  slr_cst=slr_cst_CCM
  slr_zen_ngl_cos_cmd_ln=mss_val ! [frc] Cosine solar zenith angle
  str_nbr=4
  sv_cmp_plr_ngl=.true.
  sv_cmp_tau=.true.
  sv_ntn=.false.
  top_lvl=.false.
  tst_case_HG=.false.
  tst_case_Rayleigh=.false.
  wvl_Planck=2.0e-6
  wvl_obs_mpr=0.5e-6
  wvl_obs_snw=0.5e-6
  
  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command line into single string
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count()           ! [nbr] Number of command line arguments
  arg_idx=1                 ! [idx] Counting index
  loop_while_options: do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
     dsh_key=arg_val(1:2)   ! [sng] First two characters of option
     if_dbl_dsh: if (dsh_key == '--') then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop 'Long option has no name'
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (opt_sng == 'dbg' .or. opt_sng == 'dbg_lvl' ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == 'fl_aer' .or. opt_sng == 'aer') then
           call ftn_arg_get(arg_idx,arg_val,fl_aer) ! [sng] Aerosol file
        else if (opt_sng == 'alb' .or. opt_sng == 'alb_sfc') then
           cmd_ln_alb=.not.cmd_ln_alb
           call ftn_arg_get(arg_idx,arg_val,alb_cmd_ln) ! [frc] Surface albedo
        else if (opt_sng == 'alb_sfc_NIR') then ! [frc] Direct+diffuse NIR albedo
           cmd_ln_alb_sfc_NIR_dff=.not.cmd_ln_alb_sfc_NIR_dff
           cmd_ln_alb_sfc_NIR_drc=.not.cmd_ln_alb_sfc_NIR_drc
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_NIR_dff) ! [frc] 
           alb_sfc_NIR_drc=alb_sfc_NIR_dff ! [frc] 
        else if (opt_sng == 'alb_sfc_NIR_dff') then ! [frc] Diffuse NIR albedo
           cmd_ln_alb_sfc_NIR_dff=.not.cmd_ln_alb_sfc_NIR_dff
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_NIR_dff) ! [frc] 
        else if (opt_sng == 'alb_sfc_NIR_drc') then ! [frc] Direct NIR albedo
           cmd_ln_alb_sfc_NIR_drc=.not.cmd_ln_alb_sfc_NIR_drc
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_NIR_drc) ! [frc] 
        else if (opt_sng == 'alb_sfc_vsb') then ! [frc] Direct+diffuse visible albedo
           cmd_ln_alb_sfc_vsb_dff=.not.cmd_ln_alb_sfc_vsb_dff
           cmd_ln_alb_sfc_vsb_drc=.not.cmd_ln_alb_sfc_vsb_drc
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_vsb_dff) ! [frc] 
           alb_sfc_vsb_drc=alb_sfc_vsb_dff ! [frc] 
        else if (opt_sng == 'alb_sfc_vsb_dff') then ! [frc] Diffuse visible albedo
           cmd_ln_alb_sfc_vsb_dff=.not.cmd_ln_alb_sfc_vsb_dff
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_vsb_dff) ! [frc] 
        else if (opt_sng == 'alb_sfc_vsb_drc') then ! [frc] Direct visible albedo
           cmd_ln_alb_sfc_vsb_drc=.not.cmd_ln_alb_sfc_vsb_drc
           call ftn_arg_get(arg_idx,arg_val,alb_sfc_vsb_drc) ! [frc] 
        else if (opt_sng == 'dns_snw' .or. opt_sng == 'snw_dns') then
           cmd_ln_dns_snw=.not.cmd_ln_dns_snw
           call ftn_arg_get(arg_idx,arg_val,dns_snw_cmd_ln) ! [kg m-3] Snow density
        else if (opt_sng == 'dpt_snw' .or. opt_sng == 'snw_dpt') then
           cmd_ln_dpt_snw=.not.cmd_ln_dpt_snw
           call ftn_arg_get(arg_idx,arg_val,dpt_snw_cmd_ln) ! [m] Snowpack thickness
        else if (opt_sng == 'drc_in') then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == 'drc_out') then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == 'fl_clm' .or. opt_sng == 'input' .or. opt_sng == 'input') then
           call ftn_arg_get(arg_idx,arg_val,fl_clm) ! [sng] Column profile
        else if (opt_sng == 'fl_lqd' .or. opt_sng == 'lqd') then
           call ftn_arg_get(arg_idx,arg_val,fl_lqd) ! [sng] H2O liquid file
        else if (opt_sng == 'fl_out' .or. opt_sng == 'output') then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output profile
        else if (opt_sng == 'fl_mpr' .or. opt_sng == 'mpr') then
           call ftn_arg_get(arg_idx,arg_val,fl_mpr) ! [sng] Impurity file
        else if (opt_sng == 'fl_rfl' .or. opt_sng == 'rfl') then
           flg_rfl=.not.flg_rfl
           call ftn_arg_get(arg_idx,arg_val,fl_rfl) ! [sng] Spectral reflectance file
        else if (opt_sng == 'fl_snw' .or. opt_sng == 'snw') then
           call ftn_arg_get(arg_idx,arg_val,fl_snw) ! [sng] Snow file
        else if (opt_sng == 'flg_bfb') then
           flg_bfb=.not.flg_bfb ! [flg] Produce bit-for-bit output
        else if (opt_sng == 'flt_lmn') then
           flt_lmn=.not.flt_lmn ! [flg] Use luminosity filter
        else if (opt_sng == 'flg_mie') then
           flg_mie=.not.flg_mie ! [flg] Use phase function expansion where available
        else if (opt_sng == 'flg_mpr') then
           flg_mpr=.not.flg_mpr ! [flg] Snow grain impurities are optically active
        else if (opt_sng == 'flg_msm') then
           flg_msm=.not.flg_msm ! [flg] Use multi-layer snow model if snow structure present
        else if (opt_sng == 'flg_snw') then
           flg_snw=.not.flg_snw ! [flg] Snow grains are optically active within snow model
        else if (opt_sng == 'flg_xtr_aer_snw') then
           flg_xtr_aer_snw=.true. ! [flg] Extrapolate surface aerosol into snowpack
        else if (opt_sng == 'force_lqd_phz' .or. opt_sng == 'frc_lqd') then
           force_lqd_phz=.not.force_lqd_phz ! [flg] Force liquid phase
        else if (opt_sng == 'lat_dgr' .or. opt_sng == 'lat') then ! [dgr] Latitude
           cmd_ln_lat_dgr=.not.cmd_ln_lat_dgr
           call ftn_arg_get(arg_idx,arg_val,lat_dgr_cmd_ln)
        else if (opt_sng == 'lcl_yr_day' .or. opt_sng == 'doy') then ! [day] Local year day
           cmd_ln_lcl_yr_day=.not.cmd_ln_lcl_yr_day
           call ftn_arg_get(arg_idx,arg_val,lcl_yr_day_cmd_ln)
        else if (opt_sng == 'mode_std' .or. opt_sng == 'std') then
           mode_std=.not.mode_std
        else if (opt_sng == 'mpc_CWP' .or. opt_sng == 'CWP') then
           cmd_ln_mpc_CWP=.not.cmd_ln_mpc_CWP
           call ftn_arg_get(arg_idx,arg_val,mpc_CWP_cmd_ln) ! [kg m-2] Condensed Water Path
        else if (opt_sng == 'mmr_mpr_snw' .or. opt_sng == 'mmr_mpr') then
           cmd_ln_mmr_mpr_snw=.not.cmd_ln_mmr_mpr_snw
           call ftn_arg_get(arg_idx,arg_val,mmr_mpr_snw_cmd_ln) ! [kg kg-1] Impurity mass mixing ratio in snow
        else if (opt_sng == 'odxc_snw') then
           cmd_ln_odxc_obs_snw=.not.cmd_ln_odxc_obs_snw
           call ftn_arg_get(arg_idx,arg_val,odxc_obs_snw_cmd_ln)
        else if (opt_sng == 'plk' .or. opt_sng == 'planck' .or. opt_sng == 'thermal' ) then
           call ftn_arg_get(arg_idx,arg_val,flg_Planck) ! [flg] Include thermal emission of atmosphere
        else if (opt_sng == 'cld_sat') then
           flg_cld_sat=.true. ! [flg] Force cloudy and snowy layers to be saturated
        else if (opt_sng == 'flx_frc_drc_TOA' .or. opt_sng == 'flx_frc_drc') then
           call ftn_arg_get(arg_idx,arg_val,flx_frc_drc_TOA)
        else if (opt_sng == 'sat_cld') then
           flg_sat_cld=.true. ! [flg] Force saturated layers to be cloudy
        else if (opt_sng == 'sct_lqd') then
           call ftn_arg_get(arg_idx,arg_val,flg_sct_lqd) ! [flg] Liquid cloud droplets are pure scatterers
        else if (opt_sng == 'slr_cst') then
           cmd_ln_slr_cst=.not.cmd_ln_slr_cst
           call ftn_arg_get(arg_idx,arg_val,slr_cst_cmd_ln)
        else if (opt_sng == 'slr_zen_ngl_cos') then
           cmd_ln_slr_zen_ngl_cos=.not.cmd_ln_slr_zen_ngl_cos
           call ftn_arg_get(arg_idx,arg_val,slr_zen_ngl_cos_cmd_ln) ! [frc] Cosine solar zenith angle
        else if (opt_sng == 'srm' .or. opt_sng == 'streams') then
           call ftn_arg_get(arg_idx,arg_val,str_nbr)
        else if (opt_sng == 'vpr_H2O_abs_cld') then
           call ftn_arg_get(arg_idx,arg_val,flg_vpr_H2O_abs_cld) ! [flg] H2O vapor absorption in cloud allowed
        else if (opt_sng == 'zen') then ! Shortcut for --slr_zen_ngl_cos=1.0
           cmd_ln_slr_zen_ngl_cos=.true.
           slr_zen_ngl_cos_cmd_ln=1.0 ! [frc] Cosine solar zenith angle
        else                ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif               ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     endif if_dbl_dsh            ! endif long option
     ! Handle short options
     if_sgl_dsh: if (dsh_key == '-3') then
        fl_out_fmt=nf90_format_classic ! [enm] Output file format
     else if (dsh_key == '-4') then
        fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
     else if (dsh_key == '-A') then
        flg_aer=.not.flg_aer
     else if (dsh_key == '-a') then
        call ftn_arg_get(arg_idx,arg_val,fl_aer)
     else if (dsh_key == '-B') then
        flg_bga=.not.flg_bga
     else if (dsh_key == '-b') then
        call ftn_arg_get(arg_idx,arg_val,fl_bga)
     else if (dsh_key == '-C') then
        flg_CO2=.not.flg_CO2
     else if (dsh_key == '-c') then
        call ftn_arg_get(arg_idx,arg_val,fl_CO2)
     else if (dsh_key == '-D') then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
     else if (dsh_key == '-d') then
        call ftn_arg_get(arg_idx,arg_val,fl_out)
     else if (dsh_key == '-E') then
        single_bnd_computation=.not.single_bnd_computation
     else if (dsh_key == '-e') then
        call ftn_getarg_wrp(arg_idx,arg_val)
        read (arg_val,'(i4)') bnd_dbg
     else if (dsh_key == '-F') then
        force_ice_phz=.not.force_ice_phz
     else if (dsh_key == '-f') then
        force_lqd_phz=.not.force_lqd_phz
     else if (dsh_key == '-G') then
        flg_CH4=.not.flg_CH4
     else if (dsh_key == '-g') then
        call ftn_arg_get(arg_idx,arg_val,fl_CH4)
     else if (dsh_key == '-H') then
        flg_H2O=.not.flg_H2O
     else if (dsh_key == '-h') then
        call ftn_arg_get(arg_idx,arg_val,fl_H2O)
     else if (dsh_key == '-I') then
        flg_ice=.not.flg_ice
     else if (dsh_key == '-i') then
        call ftn_arg_get(arg_idx,arg_val,fl_ice)
     else if (dsh_key == '-J') then
        flg_Planck=.not.flg_Planck
     else if (dsh_key == '-j') then
        lamber=.false.      ! DISORT --> DISORT2
     else if (dsh_key == '-K') then
        flg_O2O2=.not.flg_O2O2
     else if (dsh_key == '-k') then
        call ftn_arg_get(arg_idx,arg_val,fl_O2O2)
     else if (dsh_key == '-L') then
        flg_lqd=.not.flg_lqd
     else if (dsh_key == '-l') then
        call ftn_arg_get(arg_idx,arg_val,fl_lqd)
     else if (dsh_key == '-M') then
        cmd_ln_odxc_obs_aer=.not.cmd_ln_odxc_obs_aer
        call ftn_arg_get(arg_idx,arg_val,odxc_obs_aer_cmd_ln)
     else if (dsh_key == '-m') then
        cmd_ln_mpc_CWP=.not.cmd_ln_mpc_CWP
        call ftn_arg_get(arg_idx,arg_val,mpc_CWP_cmd_ln)
     else if (dsh_key == '-N') then
        flt_nst=.not.flt_nst
     else if (dsh_key == '-n') then
        call ftn_arg_get(arg_idx,arg_val,fl_nst)
     else if (dsh_key == '-O') then
        flg_O2=.not.flg_O2
     else if (dsh_key == '-o') then
        call ftn_arg_get(arg_idx,arg_val,fl_O2)
     else if (dsh_key == '-P') then
        tst_case_HG=.not.tst_case_HG
     else if (dsh_key == '-p') then
        call ftn_arg_get(arg_idx,arg_val,fl_clm)
     else if (dsh_key == '-Q') then
        flg_H2OH2O=.not.flg_H2OH2O
        ! else if (dsh_key == '-q') then
     else if (dsh_key == '-R') then
        flg_Rayleigh=.not.flg_Rayleigh
     else if (dsh_key == '-r') then
        cmd_ln_alb=.not.cmd_ln_alb
        call ftn_arg_get(arg_idx,arg_val,alb_cmd_ln)
     else if (dsh_key == '-S') then
        cmd_ln_slr_cst=.not.cmd_ln_slr_cst
        call ftn_arg_get(arg_idx,arg_val,slr_cst_cmd_ln)
     else if (dsh_key == '-s') then
        call ftn_arg_get(arg_idx,arg_val,str_nbr)
     else if (dsh_key == '-T') then
        call ftn_arg_get(arg_idx,arg_val,fl_slr)
     else if (dsh_key == '-t') then
        tst_case_Rayleigh=.not.tst_case_Rayleigh
     else if (dsh_key == '-U') then
        flg_O2N2=.not.flg_O2N2
     else if (dsh_key == '-u') then
        call ftn_arg_get(arg_idx,arg_val,plr_nbr)
        sv_cmp_plr_ngl=.false.
     else if (dsh_key == '-V') then
        call ftn_arg_get(arg_idx,arg_val,brdf_typ) ! BRDF type
     else if (dsh_key == '-v') then
        call ftn_arg_get(arg_idx,arg_val,fl_brdf) ! BRDF file
     else if (dsh_key == '-W') then
        flg_O3=.not.flg_O3
     else if (dsh_key == '-w') then
        call ftn_arg_get(arg_idx,arg_val,fl_O3)
     else if (dsh_key == '-X') then
        flg_NO2=.not.flg_NO2
     else if (dsh_key == '-x') then
        call ftn_arg_get(arg_idx,arg_val,fl_NO2)
     else if (dsh_key == '-Y') then
        flg_OH=.not.flg_OH
     else if (dsh_key == '-y') then
        call ftn_arg_get(arg_idx,arg_val,fl_OH)
     else if (dsh_key == '-Z') then
        call ftn_arg_get(arg_idx,arg_val,azi_nbr)
     else if (dsh_key == '-z') then
        cmd_ln_slr_zen_ngl_cos=.not.cmd_ln_slr_zen_ngl_cos
        call ftn_arg_get(arg_idx,arg_val,slr_zen_ngl_cos_cmd_ln)
     else                   ! Option not recognized
        arg_idx=arg_idx-1   ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif if_sgl_dsh       ! endif arg_val
  end do loop_while_options ! end while (arg_idx <= arg_nbr)
  
  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     call ftn_drcpfx(drc_in,fl_CO2) ! [sng] CO2 file
     call ftn_drcpfx(drc_in,fl_H2O) ! [sng] H20 file
     call ftn_drcpfx(drc_in,fl_H2OH2O) ! [sng] H2OH2O file
     call ftn_drcpfx(drc_in,fl_NO2) ! [sng] NO2 file
     call ftn_drcpfx(drc_in,fl_O2) ! [sng] O2 file
     call ftn_drcpfx(drc_in,fl_O2O2) ! [sng] O2O2 file
     call ftn_drcpfx(drc_in,fl_O3) ! [sng] O3 file
     call ftn_drcpfx(drc_in,fl_OH) ! [sng] OH file
     call ftn_drcpfx(drc_in,fl_CH4) ! [sng] CH4 file
     call ftn_drcpfx(drc_in,fl_aer) ! [sng] Aerosol file
     call ftn_drcpfx(drc_in,fl_bga) ! [sng] Background aerosol file
     call ftn_drcpfx(drc_in,fl_clm) ! [sng] Column file
     call ftn_drcpfx(drc_in,fl_ice) ! [sng] H2O ice file
     call ftn_drcpfx(drc_in,fl_lqd) ! [sng] H2O liquid file
     call ftn_drcpfx(drc_in,fl_lmn) ! [sng] Luminosity file
     call ftn_drcpfx(drc_in,fl_nst) ! [sng] Instrument file
     call ftn_drcpfx(drc_in,fl_rfl) ! [sng] Spectral reflectance file
     call ftn_drcpfx(drc_in,fl_mpr) ! [sng] Impurity file
     call ftn_drcpfx(drc_in,fl_snw) ! [sng] Snow file
     call ftn_drcpfx(drc_in,fl_slr) ! [sng] Solar spectrum file
     call ftn_drcpfx(drc_in,fl_brdf) ! [sng] BRDF file
  endif                     ! endif drc_in
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file
  
  flg_Herzberg=flg_O3
  mmn_nbr=str_nbr           ! # moments always equals # streams
  if (sv_cmp_plr_ngl) then
     plr_nbr=str_nbr
  endif
  ! Set permanent indices to avoid barenaked constants like '1'
  plr_ndr=1
  plr_zen=plr_nbr
  
  ! Write diagnostic strings
  call ftn_strcpylsc(src_rfr_sng,'Model reference is Zender et al. (1997) (ZBP97)'//nlc)
  if (top_lvl) then
     call ftn_strcpylsc(stt_top_lvl,'Pure absorbing layer to space at model top not yet implemented')
  else
     call ftn_strcpylsc(stt_top_lvl,'Pure absorbing layer to space at model top not yet implemented')
  endif
  if (flg_CO2) then
     call ftn_strcpylsc(stt_CO2,'CO2 line absorption: Malkmus narrow band parameters from '//fl_CO2) 
  else
     call ftn_strcpylsc(stt_CO2,'CO2 line absorption: Off')
  endif
  if (flg_H2O) then
     call ftn_strcpylsc(stt_H2O,'H2O line absorption: Malkmus narrow band parameters from '//fl_H2O)
  else
     call ftn_strcpylsc(stt_H2O,'H2O line absorption: Off')
  endif
  if (flg_H2OH2O) then
     call ftn_strcpylsc(stt_H2OH2O,'H2O-H2O dimer absorption: Continuum absorption cross sections from '//fl_H2OH2O)
  else
     call ftn_strcpylsc(stt_H2OH2O,'H2O-H2O dimer absorption: Off')
  endif
  if (flg_OH) then
     call ftn_strcpylsc(stt_OH,'OH line absorption: Malkmus narrow band parameters from '//fl_OH)
  else
     call ftn_strcpylsc(stt_OH,'OH line absorption: Off')
  endif
  if (flg_CH4) then
     call ftn_strcpylsc(stt_CH4,'CH4 line absorption: Malkmus narrow band parameters from '//fl_CH4)
  else
     call ftn_strcpylsc(stt_CH4,'CH4 line absorption: Off')
  endif
  if (flg_O2) then
     call ftn_strcpylsc(stt_O2,'O2 line absorption: Malkmus narrow band parameters from '//fl_O2)
  else
     call ftn_strcpylsc(stt_O2,'O2 line absorption: Off')
  endif
  if (flg_Herzberg) then
     call ftn_strcpylsc(stt_Herzberg,'O2 Herzberg bands: Continuum absorption cross sections from '//fl_O3)
  else
     call ftn_strcpylsc(stt_Herzberg,'O2 Herzberg bands: Off')
  endif
  if (flg_O3) then
     call ftn_strcpylsc(stt_O3,'O3 Hartley, Huggins, and Chappuis bands: Continuum absorption cross sections from '//fl_O3)
  else
     call ftn_strcpylsc(stt_O3,'O3 Hartley, Huggins, and Chappuis bands: Off')
  endif
  if (flg_O2O2) then
     call ftn_strcpylsc(stt_O2O2,'O2-O2 collision-induced absorption: Continuum absorption cross sections from '//fl_O2O2)
  else
     call ftn_strcpylsc(stt_O2O2,'O2-O2 collision-induced absorption: Off')
  endif
  if (flg_O2N2) then
     call ftn_strcpylsc(stt_O2N2,'O2-N2 collision-induced absorption: Occurs in 1.26 micron with 0.2 efficiency of O2-O2')
  else
     call ftn_strcpylsc(stt_O2N2,'O2-N2 collision-induced absorption: Off')
  endif
  if (flg_NO2) then
     call ftn_strcpylsc(stt_NO2,'NO2 absorption: Continuum absorption cross sections from '//fl_NO2)
  else
     call ftn_strcpylsc(stt_NO2,'NO2 absorption: Off')
  endif
  if (flg_ice) then
     call ftn_strcpylsc(stt_ice,'Ice water crystal scattering and absorption: Mie theory from '//fl_ice)
  else
     call ftn_strcpylsc(stt_ice,'Ice water crystal scattering and absorption: Off')
  endif
  if (flg_lqd) then
     call ftn_strcpylsc(stt_lqd,'Liquid water droplet scattering and absorption: Mie theory from '//fl_lqd)
  else
     call ftn_strcpylsc(stt_lqd,'Liquid water droplet scattering and absorption: Off')
  endif
  if (flg_cld_sat) then
     call ftn_strcpylsc(stt_cld_sat,'Saturation in clouds: Forcing cloudy levels to be saturated')
  else                      ! not flg_cld_sat
     call ftn_strcpylsc(stt_cld_sat,'Saturation in clouds: Allowing cloudy levels to be unsaturated')
  endif                     ! not flg_cld_sat
  if (flg_sat_cld) then
     call ftn_strcpylsc(stt_sat_cld,'Clouds in saturation: Forcing saturated levels to be cloudy')
  else                      ! not flg_sat_cld
     call ftn_strcpylsc(stt_sat_cld,'Clouds in saturation: Allowing saturated levels to be clear')
  endif                     ! not flg_sat_cld
  if (flg_sct_lqd) then
     call ftn_strcpylsc(stt_sct_lqd,'Absorption by cloud droplets: Not allowed')
  else                      ! not flg_sct_lqd
     call ftn_strcpylsc(stt_sct_lqd,'Absorption by cloud droplets: Allowed')
  endif                     ! not flg_sct_lqd
  if (flg_vpr_H2O_abs_cld) then
     call ftn_strcpylsc(stt_vpr_H2O_abs_cld,'Water vapor absorption in clouds: Allowed')
  else                      ! not flg_vpr_H2O_abs_cld
     call ftn_strcpylsc(stt_vpr_H2O_abs_cld,'Water vapor absorption in clouds: Not allowed')
  endif                     ! not flg_vpr_H2O_abs_cld
  if (flg_aer) then
     call ftn_strcpylsc(stt_aer,'Aerosol scattering and absorption: Mie theory from '//fl_aer)
  else
     call ftn_strcpylsc(stt_aer,'Aerosol scattering and absorption: Off')
  endif
  if (flg_bga) then
     call ftn_strcpylsc(stt_bga,'Background aerosol scattering and absorption: Mie theory from '//fl_bga)
  else
     call ftn_strcpylsc(stt_bga,'Background aerosol scattering and absorption: Off')
  endif
  if (flg_mpr) then
     call ftn_strcpylsc(stt_mpr,'Snow impurities: Impurity optical properties from '//fl_mpr)
  else
     call ftn_strcpylsc(stt_mpr,'Snow impurities: Off')
  endif
  if (flg_snw) then
     call ftn_strcpylsc(stt_snw,'Snow grains: Snow optical properties from '//fl_snw)
  else
     call ftn_strcpylsc(stt_snw,'Snow grains: Off')
  endif
  if (flg_mie) then
     call ftn_strcpylsc(stt_mie,'Phase function: Legendre expansion of Mie phase function')
  else
     call ftn_strcpylsc(stt_mie,'Phase function: Henyey-Greenstein approximation')
  endif
  if (flg_Planck) then
     call ftn_strcpylsc(stt_Planck,'Thermal emission of atmosphere: Included for lambda > 2 microns')
  else
     call ftn_strcpylsc(stt_Planck,'Thermal emission of atmosphere: Off')
  endif
  if (flg_Rayleigh) then
     call ftn_strcpylsc(stt_Rayleigh,'Rayleigh scattering: method of Lenoble (1993)')
  else
     call ftn_strcpylsc(stt_Rayleigh,'Rayleigh scattering: Off')
  endif
  if (flt_lmn) then
     call ftn_strcpylsc(stt_flt_lmn,'Luminosity spectral response function: Curve from '//fl_lmn)
  else
     call ftn_strcpylsc(stt_flt_lmn,'Luminosity spectral spectral response function: Off')
  endif
  if (flt_nst) then
     call ftn_strcpylsc(stt_flt_nst,'Instrument filter spectral response function: instrument from '//fl_nst)
  else
     call ftn_strcpylsc(stt_flt_nst,'Instrument filter spectral response function: Off')
  endif
  if (flg_rfl) then
     call ftn_strcpylsc(stt_rfl,'Surface reflectance: Spectral surface reflectance from '//fl_rfl)
  else
     call ftn_strcpylsc(stt_rfl,'Surface reflectance: From broadband parameters')
  endif
  call ftn_strcpylsc(stt_slr,'TOA solar spectrum from '//fl_slr)
  
  ! Ingest fl_clm
  rcd=nf90_wrp_open(fl_clm,nf90_nowrite,nc_id)
  ! Get global attributes
  rcd=rcd+nf90_get_att(nc_id,nf90_global,'prf_sng',prf_sng)
  
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dmn_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'levp',levp_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_dmn_id,len=lev_atm_nbr),sbr_nm//": inquire_dim lev")
  if (lev_atm_nbr>lev_nbr_max) stop 'lev_atm_nbr>lev_nbr_max'
  srt_one(:)=1
  srt_one_one(:)=1
  srt_one_one_one(:)=1
  cnt_lev(1)=lev_atm_nbr
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,levp_dmn_id,len=levp_atm_nbr),sbr_nm//": inquire_dim levp")
  if (levp_atm_nbr>levp_nbr_max) stop 'levp_atm_nbr>levp_nbr_max'
  cnt_levp(1)=levp_atm_nbr

  ! Get potentially missing dimension IDs and sizes
  ! First set layer structure assuming no snow layers
  lev_nbr=lev_atm_nbr
  ! NB: lev_snw_nbr must be non-zero to be a dimension in netCDF3 output files
  ! Otherwise netCDF3 dim_def() thinks size=0 indicates an unlimited dimension
  ! netCDF3 support for only one unlimited dimension then causes more problems...
  lev_snw_nbr=0
  levp_snw_nbr=0
  ! Modify layer structures if snow requested (with flg_msm) and present
  if (flg_msm) then
     rcd=nf90_wrp_inq_dimid(nc_id,'lev_snw',lev_snw_dmn_id,nf90_ebaddim)
     if(rcd == nf90_noerr) then
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_snw_dmn_id,len=lev_snw_nbr),sbr_nm//": inquire_dim lev_snw")
        if (lev_snw_nbr>0) then
           write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': Found and will use snow'
           if (lev_snw_nbr>lev_snw_nbr_max) stop 'lev_snw_nbr>lev_snw_nbr_max'
           if (lev_snw_nbr+lev_atm_nbr>lev_nbr_max) stop 'lev_snw_nbr+lev_atm_nbr>lev_nbr_max'
           rcd=nf90_wrp_inq_dimid(nc_id,'levp_snw',levp_snw_dmn_id)
           ! Account for snow in vertical grid
           lev_nbr=lev_atm_nbr+lev_snw_nbr
           levp_snw_nbr=lev_snw_nbr+1
        endif ! lev_snw_nbr>0
     endif ! lev_snw dimension exists
     if (lev_snw_nbr<0) stop 'lev_snw_nbr<0'
     if (lev_snw_nbr==0) then
        write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': Snow not found, not used'
        ! Set flg_msm false internally to skip unnecessary work
        flg_msm=.false.
        flg_snw=.false.
     endif
     ! NB: flg_msm meaning changes hereafter from
     ! [flg] Use multi-layer snow model if snow structure present
     ! to
     ! [flg] Snow model is requested, present, and active
  endif ! !flg_msm
  if (flg_msm) then
     call ftn_strcpylsc(stt_msm,'Multi-layer snow model: On')
  else
     call ftn_strcpylsc(stt_msm,'Multi-layer snow model: Off')
  endif
  levp_nbr=lev_nbr+1

  ! Set permanent indices to avoid barenaked constants like '1'
  lev_sfc=lev_nbr
  levp_sfc=levp_nbr
  lev_TOA=1
  levp_TOA=1
  
  ! Array dimensions: lev
  allocate(RH_lqd(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for RH_lqd"
  allocate(alt(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for alt"
  allocate(lev(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lev"
  allocate(odal_obs_aer(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_obs_aer"
  allocate(odal_obs_bga(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_obs_bga"
  allocate(odsl_obs_aer(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_obs_aer"
  allocate(odsl_obs_bga(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_obs_bga"
  allocate(odxl_obs_aer(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxl_obs_aer"
  allocate(odxl_obs_bga(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxl_obs_bga"
  allocate(odal_obs_mpr(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_obs_mpr"
  allocate(odsl_obs_mpr(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_obs_mpr"
  allocate(odxl_obs_mpr(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxl_obs_mpr"
  allocate(odal_obs_snw(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_obs_snw"
  allocate(odsl_obs_snw(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_obs_snw"
  allocate(odxl_obs_snw(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxl_obs_snw"
  allocate(tpt(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tpt"
  
  ! Array dimensions: levp
  allocate(alt_ntf(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for alt_ntf"
  allocate(levp(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for levp"
  allocate(tpt_ntf(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tpt_ntf"

  ! Array dimensions: lev_snw
  allocate(lev_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lev_snw"
  allocate(dns_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for dns_snw"
  allocate(dpt_dlt_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for dpt_dlt_snw"
  allocate(dpt_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for dpt_snw"
  allocate(mmr_mpr_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for mmr_mpr_snw"
  allocate(rds_ffc_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for rds_ffc_snw"
  allocate(tpt_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tpt_snw"
  allocate(foo_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for foo_snw"
  
  ! Array dimensions: levp_snw
  allocate(levp_snw(levp_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for levp_snw"
  allocate(dpt_ntf_snw(levp_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for dpt_ntf_snw"
  allocate(tpt_ntf_snw(levp_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tpt_ntf_snw"
  
  ! Derived variables with array dimensions: lev
  allocate(odal_mpr(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_mpr"
  allocate(odsl_mpr(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_mpr"
  allocate(odal_snw(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_snw"
  allocate(odsl_snw(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_snw"
  allocate(odal_aer(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_aer"
  allocate(odsl_aer(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_aer"
  allocate(odal_bga(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_bga"
  allocate(odsl_bga(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_bga"
  allocate(odal_ice(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_ice"
  allocate(odsl_ice(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_ice"
  allocate(odal_lqd(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odal_lqd"
  allocate(odsl_lqd(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odsl_lqd"

  ! Derived variables with array dimensions: lev_snw
  allocate(mpl_mpr(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for mpl_mpr"
  allocate(mpl_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for mpl_snw"
  allocate(frc_sfc_snw(lev_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for frc_sfc_snw"

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'RH_lqd',RH_lqd_id)
  if (.not.cmd_ln_alb_sfc_NIR_dff) rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_NIR_dff',alb_sfc_NIR_dff_id)
  if (.not.cmd_ln_alb_sfc_NIR_drc) rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_NIR_drc',alb_sfc_NIR_drc_id)
  if (.not.cmd_ln_alb_sfc_vsb_dff) rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_vsb_dff',alb_sfc_vsb_dff_id)
  if (.not.cmd_ln_alb_sfc_vsb_drc) rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_vsb_drc',alb_sfc_vsb_drc_id)
  rcd=nf90_wrp_inq_varid(nc_id,'alt',alt_id)
  rcd=nf90_wrp_inq_varid(nc_id,'alt_cld_btm',alt_cld_btm_id)
  rcd=nf90_wrp_inq_varid(nc_id,'alt_cld_thick',alt_cld_thick_id)
  rcd=nf90_wrp_inq_varid(nc_id,'alt_ntf',alt_ntf_id)
  rcd=nf90_wrp_inq_varid(nc_id,'frc_ice',frc_ice_id)
  rcd=nf90_wrp_inq_varid(nc_id,'frc_ice_ttl',frc_ice_ttl_id)
  rcd=nf90_wrp_inq_varid(nc_id,'grv',grv_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lat_dgr',lat_dgr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lcl_time_hr',lcl_time_hr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lcl_yr_day',lcl_yr_day_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lev',lev_id)
  rcd=nf90_wrp_inq_varid(nc_id,'levp',levp_id)
  rcd=nf90_wrp_inq_varid(nc_id,'lon_dgr',lon_dgr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mmw_mst_air',mmw_mst_air_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpc_CWP',mpc_CWP_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_CO2',mpl_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_CWP',mpl_CWP_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_H2O',mpl_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_IWP',mpl_IWP_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_LWP',mpl_LWP_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_O2',mpl_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_OH',mpl_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_CH4',mpl_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_aer',mpl_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_bga',mpl_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'mpl_mst_air',mpl_mst_air_id)
  rcd=nf90_wrp_inq_varid(nc_id,'npl_H2OH2O',npl_H2OH2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'npl_NO2',npl_NO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'npl_O2',npl_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'npl_O2O2',npl_O2O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'npl_O3',npl_O3_id)
  rcd=nf90_wrp_inq_varid(nc_id,'odxc_obs_aer',odxc_obs_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'odxc_obs_bga',odxc_obs_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'odxl_obs_aer',odxl_obs_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'odxl_obs_bga',odxl_obs_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'prs',prs_id)
  rcd=nf90_wrp_inq_varid(nc_id,'prs_dlt',prs_dlt_id)
  rcd=nf90_wrp_inq_varid(nc_id,'prs_ntf',prs_ntf_id)
  rcd=nf90_wrp_inq_varid(nc_id,'q_CO2',q_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'q_H2O',q_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'q_O2',q_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'q_OH',q_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'q_CH4',q_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'slr_zen_ngl_cos',slr_zen_ngl_cos_id)
  rcd=nf90_wrp_inq_varid(nc_id,'spc_heat_mst_air',spc_heat_mst_air_id)
  rcd=nf90_wrp_inq_varid(nc_id,'tpt',tpt_id)
  rcd=nf90_wrp_inq_varid(nc_id,'tpt_ntf',tpt_ntf_id)
  rcd=nf90_wrp_inq_varid(nc_id,'tpt_skn',tpt_skn_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_obs_aer',wvl_obs_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_obs_bga',wvl_obs_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'xnt_fac',xnt_fac_id)
  if (flg_msm) then
     rcd=nf90_wrp_inq_varid(nc_id,'lev_snw',lev_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'levp_snw',levp_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'dns_snw',dns_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'dpt_dlt_snw',dpt_dlt_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'dpt_ntf_snw',dpt_ntf_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'dpt_snw',dpt_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'mmr_mpr_snw',mmr_mpr_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'rds_ffc_snw',rds_ffc_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'tpt_snw',tpt_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'tpt_ntf_snw',tpt_ntf_snw_id)
     rcd=nf90_wrp_inq_varid(nc_id,'foo_snw',foo_snw_id)
  endif ! !flg_msm
  
  if (.not.mode_std) then
     rcd=nf90_wrp_inq_varid(nc_id,'sea_mask',ocn_msk_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wnd_spd',wnd_spd_id)
  endif
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,lat_dgr_id,lat_dgr),"gv lat_dgr")
  rcd=nf90_wrp(nf90_get_var(nc_id,lcl_time_hr_id,lcl_time_hr),"gv lcl_time_hr")
  rcd=nf90_wrp(nf90_get_var(nc_id,lcl_yr_day_id,lcl_yr_day),"gv lcl_yr_day")
  rcd=nf90_wrp(nf90_get_var(nc_id,lon_dgr_id,lon_dgr),"gv lon_dgr")
  rcd=nf90_wrp(nf90_get_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos),"gv slr_zen_ngl_cos")

  rcd=nf90_wrp(nf90_get_var(nc_id,RH_lqd_id,RH_lqd,srt_one,cnt_lev),"gv RH_lqd")
  if (.not.cmd_ln_alb_sfc_NIR_dff) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_NIR_dff_id,alb_sfc_NIR_dff),"gv alb_sfc_NIR_dff")
  if (.not.cmd_ln_alb_sfc_NIR_drc) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_NIR_drc_id,alb_sfc_NIR_drc),"gv alb_sfc_NIR_drc")
  if (.not.cmd_ln_alb_sfc_vsb_dff) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_vsb_dff_id,alb_sfc_vsb_dff),"gv alb_sfc_vsb_dff")
  if (.not.cmd_ln_alb_sfc_vsb_drc) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_vsb_drc_id,alb_sfc_vsb_drc),"gv alb_sfc_vsb_drc")
  rcd=nf90_wrp(nf90_get_var(nc_id,alt_cld_btm_id,alt_cld_btm),"gv alt_cld_btm")
  rcd=nf90_wrp(nf90_get_var(nc_id,alt_cld_thick_id,alt_cld_thick),"gv alt_cld_thick")
  rcd=nf90_wrp(nf90_get_var(nc_id,alt_id,alt,srt_one,cnt_lev),"gv alt")
  rcd=nf90_wrp(nf90_get_var(nc_id,alt_ntf_id,alt_ntf,srt_one,cnt_levp),"gv alt_ntf")
  rcd=nf90_wrp(nf90_get_var(nc_id,frc_ice_id,frc_ice,srt_one,cnt_lev),"gv frc_ice")
  rcd=nf90_wrp(nf90_get_var(nc_id,frc_ice_ttl_id,frc_ice_ttl),"gv frc_ice_ttl")
  rcd=nf90_wrp(nf90_get_var(nc_id,grv_id,grv,srt_one,cnt_lev),"gv grv")
  rcd=nf90_wrp(nf90_get_var(nc_id,lev_id,lev,srt_one,cnt_lev),"gv lev")
  rcd=nf90_wrp(nf90_get_var(nc_id,levp_id,levp,srt_one,cnt_levp),"gv levp")
  rcd=nf90_wrp(nf90_get_var(nc_id,mmw_mst_air_id,mmw_mst_air,srt_one,cnt_lev),"gv mmw_mst_air")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpc_CWP_id,mpc_CWP),"gv mpc_CWP")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_CO2_id,mpl_CO2,srt_one,cnt_lev),"gv mpl_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_CWP_id,mpl_CWP,srt_one,cnt_lev),"gv mpl_CWP")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_H2O_id,mpl_H2O,srt_one,cnt_lev),"gv mpl_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_IWP_id,mpl_IWP,srt_one,cnt_lev),"gv mpl_IWP")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_LWP_id,mpl_LWP,srt_one,cnt_lev),"gv mpl_LWP")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_O2_id,mpl_O2,srt_one,cnt_lev),"gv mpl_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_OH_id,mpl_OH,srt_one,cnt_lev),"gv mpl_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_CH4_id,mpl_CH4,srt_one,cnt_lev),"gv mpl_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_aer_id,mpl_aer,srt_one,cnt_lev),"gv mpl_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_bga_id,mpl_bga,srt_one,cnt_lev),"gv mpl_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,mpl_mst_air_id,mpl_mst_air,srt_one,cnt_lev),"gv mpl_mst_air")
  rcd=nf90_wrp(nf90_get_var(nc_id,npl_H2OH2O_id,npl_H2OH2O,srt_one,cnt_lev),"gv npl_H2OH2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,npl_NO2_id,npl_NO2,srt_one,cnt_lev),"gv npl_NO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,npl_O2O2_id,npl_O2O2,srt_one,cnt_lev),"gv npl_O2O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,npl_O2_id,npl_O2,srt_one,cnt_lev),"gv npl_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,npl_O3_id,npl_O3,srt_one,cnt_lev),"gv npl_O3")
  rcd=nf90_wrp(nf90_get_var(nc_id,odxc_obs_aer_id,odxc_obs_aer),"gv odxc_obs_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,odxc_obs_bga_id,odxc_obs_bga),"gv odxc_obs_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,odxl_obs_aer_id,odxl_obs_aer,srt_one,cnt_lev),"gv odxl_obs_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,odxl_obs_bga_id,odxl_obs_bga,srt_one,cnt_lev),"gv odxl_obs_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,prs_dlt_id,prs_dlt,srt_one,cnt_lev),"gv prs_dlt")
  rcd=nf90_wrp(nf90_get_var(nc_id,prs_id,prs,srt_one,cnt_lev),"gv prs")
  rcd=nf90_wrp(nf90_get_var(nc_id,prs_ntf_id,prs_ntf,srt_one,cnt_levp),"gv prs_ntf")
  rcd=nf90_wrp(nf90_get_var(nc_id,q_CO2_id,q_CO2,srt_one,cnt_lev),"gv q_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,q_H2O_id,q_H2O,srt_one,cnt_lev),"gv q_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,q_O2_id,q_O2,srt_one,cnt_lev),"gv q_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,q_OH_id,q_OH,srt_one,cnt_lev),"gv q_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,q_CH4_id,q_CH4,srt_one,cnt_lev),"gv q_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,spc_heat_mst_air_id,spc_heat_mst_air,srt_one,cnt_lev),"gv spc_heat_mst_air")
  rcd=nf90_wrp(nf90_get_var(nc_id,tpt_id,tpt,srt_one,cnt_lev),"gv tpt")
  rcd=nf90_wrp(nf90_get_var(nc_id,tpt_ntf_id,tpt_ntf,srt_one,cnt_levp),"gv tpt_ntf")
  rcd=nf90_wrp(nf90_get_var(nc_id,tpt_skn_id,tpt_skn),"gv tpt_skn")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_obs_aer_id,wvl_obs_aer),"gv wvl_obs_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_obs_bga_id,wvl_obs_bga),"gv wvl_obs_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,xnt_fac_id,xnt_fac),"gv xnt_fac")
  if (flg_msm) then
     rcd=nf90_wrp(nf90_get_var(nc_id,lev_snw_id,lev_snw),"gv lev_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,levp_snw_id,levp_snw),"gv levp_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,dns_snw_id,dns_snw),"gv dns_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,dpt_dlt_snw_id,dpt_dlt_snw),"gv dpt_dlt_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,dpt_ntf_snw_id,dpt_ntf_snw),"gv dpt_ntf_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,dpt_snw_id,dpt_snw),"gv dpt_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,mmr_mpr_snw_id,mmr_mpr_snw),"gv mmr_mpr_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,rds_ffc_snw_id,rds_ffc_snw),"gv rds_ffc_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,tpt_snw_id,tpt_snw),"gv tpt_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,tpt_ntf_snw_id,tpt_ntf_snw),"gv tpt_ntf_snw")
     rcd=nf90_wrp(nf90_get_var(nc_id,foo_snw_id,foo_snw),"gv foo_snw")
     if (cmd_ln_dns_snw) then
        do lev_snw_idx=1,lev_snw_nbr
           dns_snw(lev_snw_idx)=dns_snw_cmd_ln ! [kg m-3] Snow density
        enddo ! end loop over lev_snw
     endif ! !cmd_ln_dpt_snw)
     if (cmd_ln_dpt_snw) then
        ! Uniformly stretch snowpack to user-specified depth
        dpt_snw_crc_fct=dpt_snw_cmd_ln/dpt_ntf_snw(levp_snw_nbr)
        do lev_snw_idx=1,lev_snw_nbr
           dpt_dlt_snw(lev_snw_idx)=dpt_dlt_snw(lev_snw_idx)*dpt_snw_crc_fct
           dpt_snw(lev_snw_idx)=dpt_snw(lev_snw_idx)*dpt_snw_crc_fct
           lev_snw(lev_snw_idx)=lev_snw(lev_snw_idx)*dpt_snw_crc_fct
        enddo ! end loop over lev_snw
        do levp_snw_idx=1,levp_snw_nbr
           dpt_ntf_snw(levp_snw_idx)=dpt_ntf_snw(levp_snw_idx)*dpt_snw_crc_fct
           levp_snw(levp_snw_idx)=levp_snw(levp_snw_idx)*dpt_snw_crc_fct
        enddo ! end loop over levp_snw
     endif ! !cmd_ln_dpt_snw)
  endif ! !flg_msm
  if (.not.mode_std) then
     rcd=nf90_wrp(nf90_get_var(nc_id,ocn_msk_id,ocn_msk),"gv ocn_msk")
     rcd=nf90_wrp(nf90_get_var(nc_id,wnd_spd_id,wnd_spd),"gv wnd_spd")
  else
     ocn_msk=0.0 ! CEWI lf95
     wnd_spd=0.0 ! [m s-1]
  endif

  ! Fill-in atmospheric profile arrays inside snowpack
  if (flg_msm) then
     ! Override impurities with command-line arguments
     if (cmd_ln_mmr_mpr_snw) then
        mpr_lvl_nbr=0
        do lev_snw_idx=1,lev_snw_nbr
           if (mmr_mpr_snw(lev_snw_idx)>0.0) mpr_lvl_nbr=mpr_lvl_nbr+1
        enddo                  ! end loop over lev
        if (mpr_lvl_nbr==0) stop 'mpr_lvl_nbr==0'
        do lev_snw_idx=1,lev_snw_nbr
           if (mmr_mpr_snw(lev_snw_idx)>0.0) mmr_mpr_snw(lev_snw_idx)=mmr_mpr_snw_cmd_ln
        enddo                  ! end loop over lev
     endif                     ! end if overriding CLM profile mmr_mpr_snw
     ! Bolt-on coordinate grid
     tpt_vrt_sfc=(1.0+eps_H2O_rcp_m1*q_H2O(lev_atm_nbr))*tpt_ntf(levp_atm_nbr)
     dns_mst_air_sfc=prs_ntf(levp_atm_nbr)/(gas_cst_dry_air*tpt_vrt_sfc)
     do lev_snw_idx=1,lev_snw_nbr
        lev_idx=lev_snw_idx+lev_atm_nbr
        tpt(lev_idx)=tpt_snw(lev_snw_idx)
        ! NB: Add snow mass to pressure beneath pore close-off?
        prs_dlt(lev_idx)=dpt_dlt_snw(lev_snw_idx)*dns_mst_air_sfc*grv(lev_atm_nbr)
     enddo ! end loop over lev_snw
     do levp_idx=levp_atm_nbr+1,levp_nbr
        lev_idx=levp_idx-1
        prs_ntf(levp_idx)=prs_ntf(levp_idx-1)+prs_dlt(lev_idx)
        levp(levp_idx)=prs_ntf(levp_idx)
     enddo ! end loop over levp
     do lev_idx=lev_atm_nbr+1,lev_nbr
        levp_idx=lev_idx
        prs(lev_idx)=0.5*(prs_ntf(levp_idx)+prs_ntf(levp_idx+1))
        lev(lev_idx)=prs(lev_idx)
     enddo ! end loop over lev
     do levp_snw_idx=1,levp_snw_nbr
        levp_idx=levp_atm_nbr+levp_snw_idx-1
        tpt_ntf(levp_idx)=tpt_ntf_snw(levp_snw_idx)
        alt_ntf(levp_idx)=-dpt_ntf_snw(levp_snw_idx)
     enddo ! end loop over levp_snw
     ! Fraction of lowest atmospheric layer air in each snowpack layer
     do lev_snw_idx=1,lev_snw_nbr
        lev_idx=lev_snw_idx+lev_atm_nbr
        ! Vertical fraction is straightforward
        frc_sfc_snw(lev_snw_idx)=dpt_dlt_snw(lev_snw_idx)/(alt_ntf(levp_atm_nbr-1)-alt_ntf(levp_atm_nbr))
        ! Remove volume taken by non-pores (i.e., by snow)
        frc_sfc_snw(lev_snw_idx)=frc_sfc_snw(lev_snw_idx)*dns_snw(lev_snw_idx)/dns_H2O_ice_std
        alt(lev_idx)=-dpt_snw(lev_snw_idx)
     enddo ! end loop over lev_snw
     do lev_snw_idx=1,lev_snw_nbr
        lev_idx=lev_atm_nbr+lev_snw_idx
        grv(lev_idx)=grv(lev_atm_nbr)
        spc_heat_mst_air(lev_idx)=spc_heat_mst_air(lev_atm_nbr)
        mmw_mst_air(lev_idx)=mmw_mst_air(lev_atm_nbr)
        ! Following are mass paths that SWNB2 uses for RT
        ! Copy appropriate fraction of lowest atmospheric level into snow pores
        ! Extrapolating surface-layer gases into snowpack seems always appropriate
        mpl_CO2(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_CO2(lev_atm_nbr)
        mpl_H2O(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_H2O(lev_atm_nbr)
        mpl_O2(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_O2(lev_atm_nbr)
        mpl_OH(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_OH(lev_atm_nbr)
        mpl_CH4(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_CH4(lev_atm_nbr)
        mpl_mst_air(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_mst_air(lev_atm_nbr)
        npl_NO2(lev_idx)=frc_sfc_snw(lev_snw_idx)*npl_NO2(lev_atm_nbr)
        npl_O2(lev_idx)=frc_sfc_snw(lev_snw_idx)*npl_O2(lev_atm_nbr)
        npl_O3(lev_idx)=frc_sfc_snw(lev_snw_idx)*npl_O3(lev_atm_nbr)
        npl_O2O2(lev_idx)=frc_sfc_snw(lev_snw_idx)*npl_O2O2(lev_atm_nbr)
        npl_H2OH2O(lev_idx)=frc_sfc_snw(lev_snw_idx)*npl_H2OH2O(lev_atm_nbr)
        ! No reason to put "clouds" in snowpack unless/until treatment improves
        ! LWP could be used to represent liquid in melting snowpack...
        ! However, current liquid properties assume spheres, not interstitial glop
        mpl_CWP(lev_idx)=0.0
        mpl_IWP(lev_idx)=0.0
        mpl_LWP(lev_idx)=0.0
        if (flg_xtr_aer_snw) then
           ! Extrapolate surface-layer particulates into snowpack only by request
           ! This avoids assumption that aerosol dry-dep is efficient 
           mpl_aer(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_aer(lev_atm_nbr)
           mpl_bga(lev_idx)=frc_sfc_snw(lev_snw_idx)*mpl_bga(lev_atm_nbr)
        else
           mpl_aer(lev_idx)=0.0
           mpl_bga(lev_idx)=0.0
        end if ! endif false
     enddo ! end loop over lev_snw
  endif ! !flg_msm

  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_clm,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_H2O
  rcd=nf90_wrp_open(fl_H2O,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_H2O),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_H2O>bnd_nbr_H2O_max) stop 'bnd_nbr_H2O>bnd_nbr_H2O_max'
  cnt_bnd(1)=bnd_nbr_H2O
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_H2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_H2O_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_H2O_id,A_phi_H2O,srt_one,cnt_bnd),"gv A_phi_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_H2O_id,A_psi_H2O,srt_one,cnt_bnd),"gv A_psi_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_H2O_id,B_phi_H2O,srt_one,cnt_bnd),"gv B_phi_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_H2O_id,B_psi_H2O,srt_one,cnt_bnd),"gv B_psi_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_H2O_id,S_d_abs_cff_mss_H2O,srt_one,cnt_bnd),"gv S_d_abs_cff_mss_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_H2O_id,S_p_abs_cff_mss_H2O,srt_one,cnt_bnd),"gv S_p_abs_cff_mss_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_H2O_id,wvl_max_H2O,srt_one,cnt_bnd),"gv wvl_max_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_H2O_id,wvl_min_H2O,srt_one,cnt_bnd),"gv wvl_min_H2O")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_H2O_id,wvl_ctr_H2O,srt_one,cnt_bnd),"gv wvl_ctr_H2O")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_H2O,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_CO2
  rcd=nf90_wrp_open(fl_CO2,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_CO2),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_CO2>bnd_nbr_CO2_max) stop 'bnd_nbr_CO2>bnd_nbr_CO2_max'
  cnt_bnd(1)=bnd_nbr_CO2
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_CO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_CO2_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_CO2_id,A_phi_CO2,srt_one,cnt_bnd),"gv A_phi_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_CO2_id,A_psi_CO2,srt_one,cnt_bnd),"gv A_psi_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_CO2_id,B_phi_CO2,srt_one,cnt_bnd),"gv B_phi_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_CO2_id,B_psi_CO2,srt_one,cnt_bnd),"gv B_psi_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_CO2_id,S_d_abs_cff_mss_CO2,srt_one,cnt_bnd),"gv S_d_abs_cff_mss_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_CO2_id,S_p_abs_cff_mss_CO2,srt_one,cnt_bnd),"gv S_p_abs_cff_mss_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_CO2_id,wvl_max_CO2,srt_one,cnt_bnd),"gv wvl_max_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_CO2_id,wvl_min_CO2,srt_one,cnt_bnd),"gv wvl_min_CO2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_CO2_id,wvl_ctr_CO2,srt_one,cnt_bnd),"gv wvl_ctr_CO2")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_CO2,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_OH
  rcd=nf90_wrp_open(fl_OH,nf90_nowrite,nc_id)
  
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_OH),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_OH>bnd_nbr_OH_max) stop 'bnd_nbr_OH>bnd_nbr_OH_max'
  cnt_bnd(1)=bnd_nbr_OH

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_OH_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_OH_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_OH_id,A_phi_OH,srt_one,cnt_bnd),"gv A_phi_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_OH_id,A_psi_OH,srt_one,cnt_bnd),"gv A_psi_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_OH_id,B_phi_OH,srt_one,cnt_bnd),"gv B_phi_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_OH_id,B_psi_OH,srt_one,cnt_bnd),"gv B_psi_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_OH_id,S_d_abs_cff_mss_OH,srt_one,cnt_bnd),"gv S_d_abs_cff_mss_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_OH_id,S_p_abs_cff_mss_OH,srt_one,cnt_bnd),"gv S_p_abs_cff_mss_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_OH_id,wvl_max_OH,srt_one,cnt_bnd),"gv wvl_max_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_OH_id,wvl_min_OH,srt_one,cnt_bnd),"gv wvl_min_OH")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_OH_id,wvl_ctr_OH,srt_one,cnt_bnd),"gv wvl_ctr_OH")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_OH,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_OH
  rcd=nf90_wrp_open(fl_CH4,nf90_nowrite,nc_id)
  
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_CH4),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_CH4>bnd_nbr_CH4_max) stop 'bnd_nbr_CH4>bnd_nbr_CH4_max'
  cnt_bnd(1)=bnd_nbr_CH4

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_CH4_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_CH4_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_CH4_id,A_phi_CH4,srt_one,cnt_bnd),"gv A_phi_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_CH4_id,A_psi_CH4,srt_one,cnt_bnd),"gv A_psi_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_CH4_id,B_phi_CH4,srt_one,cnt_bnd),"gv B_phi_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_CH4_id,B_psi_CH4,srt_one,cnt_bnd),"gv B_psi_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_CH4_id,S_d_abs_cff_mss_CH4,srt_one,cnt_bnd),"gv S_d_abs_cff_mss_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_CH4_id,S_p_abs_cff_mss_CH4,srt_one,cnt_bnd),"gv S_p_abs_cff_mss_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_CH4_id,wvl_max_CH4,srt_one,cnt_bnd),"gv wvl_max_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_CH4_id,wvl_min_CH4,srt_one,cnt_bnd),"gv wvl_min_CH4")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_CH4_id,wvl_ctr_CH4,srt_one,cnt_bnd),"gv wvl_ctr_CH4")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_CH4,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_O2
  rcd=nf90_wrp_open(fl_O2,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_O2),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_O2>bnd_nbr_O2_max) stop 'bnd_nbr_O2>bnd_nbr_O2_max'
  cnt_bnd(1)=bnd_nbr_O2
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_O2_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_O2_id,A_phi_O2,srt_one,cnt_bnd),"gv A_phi_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_O2_id,A_psi_O2,srt_one,cnt_bnd),"gv A_psi_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_O2_id,B_phi_O2,srt_one,cnt_bnd),"gv B_phi_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_O2_id,B_psi_O2,srt_one,cnt_bnd),"gv B_psi_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_O2_id,S_d_abs_cff_mss_O2,srt_one,cnt_bnd),"gv S_d_abs_cff_mss_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_O2_id,S_p_abs_cff_mss_O2,srt_one,cnt_bnd),"gv S_p_abs_cff_mss_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_O2_id,wvl_max_O2,srt_one,cnt_bnd),"gv wvl_max_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_O2_id,wvl_min_O2,srt_one,cnt_bnd),"gv wvl_min_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_O2_id,wvl_ctr_O2,srt_one,cnt_bnd),"gv wvl_ctr_O2")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_O2,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_O3
  rcd=nf90_wrp_open(fl_O3,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_O3),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_O3>bnd_nbr_O3_max) stop 'bnd_nbr_O3>bnd_nbr_O3_max'
  cnt_bnd(1)=bnd_nbr_O3
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_O2',abs_xsx_O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_O3_cold',abs_xsx_O3_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_O3_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_O3_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_O3_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_O2_id,abs_xsx_O2,srt_one,cnt_bnd),"gv abs_xsx_O2")
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_O3_id,abs_xsx_O3,srt_one,cnt_bnd),"gv abs_xsx_O3")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_O3_id,wvl_max_O3,srt_one,cnt_bnd),"gv wvl_max_O3")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_O3_id,wvl_min_O3,srt_one,cnt_bnd),"gv wvl_min_O3")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_O3_id,wvl_ctr_O3,srt_one,cnt_bnd),"gv wvl_ctr_O3")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_O3,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_O2O2
  rcd=nf90_wrp_open(fl_O2O2,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_O2O2),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_O2O2>bnd_nbr_O2O2_max) stop 'bnd_nbr_O2O2>bnd_nbr_O2O2_max'
  cnt_bnd(1)=bnd_nbr_O2O2
  cnt_bndp(1)=bnd_nbr_O2O2+1
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_O2O2',abs_xsx_O2O2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_O2O2_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2_dsk,srt_one,cnt_bnd),"gv abs_xsx_O2O2_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_O2O2_id,wvl_grd_O2O2,srt_one,cnt_bndp),"gv wvl_grd_O2O2")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_O2O2,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_NO2
  rcd=nf90_wrp_open(fl_NO2,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_NO2),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_NO2>bnd_nbr_NO2_max) stop 'bnd_nbr_NO2>bnd_nbr_NO2_max'
  cnt_bnd(1)=bnd_nbr_NO2
  cnt_bndp(1)=bnd_nbr_NO2+1
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_NO2',abs_xsx_NO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'qnt_yld_NO2',qnt_yld_NO2_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_NO2_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_NO2_id,abs_xsx_NO2_dsk,srt_one,cnt_bnd),"gv abs_xsx_NO2_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,qnt_yld_NO2_id,qnt_yld_NO2_dsk,srt_one,cnt_bnd),"gv qnt_yld_NO2_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_NO2_id,wvl_grd_NO2,srt_one,cnt_bndp),"gv wvl_grd_NO2")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_NO2,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_H2OH2O
  rcd=nf90_wrp_open(fl_H2OH2O,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_H2OH2O),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_H2OH2O>bnd_nbr_H2OH2O_max) stop 'bnd_nbr_H2OH2O>bnd_nbr_H2OH2O_max'
  cnt_bnd(1)=bnd_nbr_H2OH2O
  cnt_bndp(1)=bnd_nbr_H2OH2O+1
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_H2OH2O',abs_xsx_H2OH2O_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_H2OH2O_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_H2OH2O_id,abs_xsx_H2OH2O_dsk,srt_one,cnt_bnd),"gv abs_xsx_H2OH2O_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_H2OH2O_id,wvl_grd_H2OH2O,srt_one,cnt_bndp),"gv wvl_grd_H2OH2O")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_H2OH2O,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_aer
  rcd=nf90_wrp_open(fl_aer,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_aer),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_aer>bnd_nbr_aer_max) stop 'bnd_nbr_aer>bnd_nbr_aer_max'
  cnt_bnd(1)=bnd_nbr_aer
  cnt_bndp(1)=bnd_nbr_aer+1
  cnt_mmn_bnd(2)=bnd_nbr_aer
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_aer=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd(1)=mmn_nbr+1
  ! Array dimensions: bnd_nbr_aer
  allocate(abs_cff_mss_aer_dsk(bnd_nbr_aer),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_aer_dsk"
  allocate(ext_cff_mss_aer_dsk(bnd_nbr_aer),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ext_cff_mss_aer_dsk"
  allocate(sca_cff_mss_aer_dsk(bnd_nbr_aer),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_aer_dsk"
  allocate(asm_prm_aer_dsk(bnd_nbr_aer),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_aer_dsk"
  ! Array dimensions: bnd_nbr_aer+1
  allocate(wvl_grd_aer(bnd_nbr_aer+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_aer"
  ! Array dimensions: mmn_nbr,bnd_nbr_aer
  if (flg_mie_aer) then
     allocate(lgn_xpn_cff_aer_dsk(0:mmn_nbr,bnd_nbr_aer),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_aer_dsk"
  endif ! !flg_mie_aer
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'ext_cff_mss',ext_cff_mss_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_aer_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_aer_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_aer_id,abs_cff_mss_aer_dsk,srt_one,cnt_bnd),"gv abs_cff_mss_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_aer_id,asm_prm_aer_dsk,srt_one,cnt_bnd),"gv asm_prm_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,ext_cff_mss_aer_id,ext_cff_mss_aer_dsk,srt_one,cnt_bnd),"gv ext_cff_mss_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_aer_id,sca_cff_mss_aer_dsk,srt_one,cnt_bnd),"gv sca_cff_mss_aer")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_aer_id,wvl_grd_aer,srt_one,cnt_bndp),"gv wvl_grd_aer")
  ! Get optional data
  if (flg_mie_aer) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_aer_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_aer_id,lgn_xpn_cff_aer_dsk,srt_one_one,cnt_mmn_bnd),"gv lgn_xpn_cff_aer")
     do bnd_idx=1,cnt_bnd(1)
        ! Sanity check on phase function expansion
        if (abs(lgn_xpn_cff_aer_dsk(1,bnd_idx)-asm_prm_aer_dsk(bnd_idx)) > 0.40) then
           write (6,'(3a,/,2(a,i4,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR inaccurate mie phase function expansion in ',fl_aer, &
                'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx,') = ',lgn_xpn_cff_aer_dsk(1,bnd_idx), &
                ' != asm_prm(bnd_idx=',bnd_idx,') = ',asm_prm_aer_dsk(bnd_idx), &
                '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
           call abort
        endif ! endif err
     enddo ! end loop over bnd
  end if ! !flg_mie_aer
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_aer,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_bga
  rcd=nf90_wrp_open(fl_bga,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_bga),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_bga>bnd_nbr_bga_max) stop 'bnd_nbr_bga>bnd_nbr_bga_max'
  cnt_bnd(1)=bnd_nbr_bga
  cnt_bndp(1)=bnd_nbr_bga+1
  cnt_mmn_bnd(2)=bnd_nbr_bga
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_bga=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd(1)=mmn_nbr+1
  ! Array dimensions: bnd_nbr_bga
  allocate(abs_cff_mss_bga_dsk(bnd_nbr_bga),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_bga_dsk"
  allocate(sca_cff_mss_bga_dsk(bnd_nbr_bga),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_bga_dsk"
  allocate(asm_prm_bga_dsk(bnd_nbr_bga),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_bga_dsk"
  ! Array dimensions: bnd_nbr_bga+1
  allocate(wvl_grd_bga(bnd_nbr_bga+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_bga"
  ! Array dimensions: mmn_nbr,bnd_nbr_bga
  if (flg_mie_bga) then
     allocate(lgn_xpn_cff_bga_dsk(0:mmn_nbr,bnd_nbr_bga),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_bga_dsk"
  endif ! !flg_mie_bga
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_bga_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_bga_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_bga_id,abs_cff_mss_bga_dsk,srt_one,cnt_bnd),"gv abs_cff_mss_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_bga_id,asm_prm_bga_dsk,srt_one,cnt_bnd),"gv asm_prm_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_bga_id,sca_cff_mss_bga_dsk,srt_one,cnt_bnd),"gv sca_cff_mss_bga")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_bga_id,wvl_grd_bga,srt_one,cnt_bndp),"gv wvl_grd_bga")
  ! Get optional data
  if (flg_mie_bga) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_bga_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_bga_id,lgn_xpn_cff_bga_dsk,srt_one_one,cnt_mmn_bnd),"gv lgn_xpn_cff_bga")
     do bnd_idx=1,cnt_bnd(1)
        ! Sanity check on phase function expansion
        if (abs(lgn_xpn_cff_bga_dsk(1,bnd_idx)-asm_prm_bga_dsk(bnd_idx)) > 0.40) then
           write (6,'(3a,/,2(a,i4,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR inaccurate mie phase function expansion in ',fl_bga, &
                'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx,') = ',lgn_xpn_cff_bga_dsk(1,bnd_idx), &
                ' != asm_prm(bnd_idx=',bnd_idx,') = ',asm_prm_bga_dsk(bnd_idx), &
                '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
           call abort
        endif ! endif err
     enddo ! end loop over bnd
  end if ! !flg_mie_bga
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_bga,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_ice
  rcd=nf90_wrp_open(fl_ice,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_ice),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_ice>bnd_nbr_ice_max) stop 'bnd_nbr_ice>bnd_nbr_ice_max'
  cnt_bnd(1)=bnd_nbr_ice
  cnt_bndp(1)=bnd_nbr_ice+1
  cnt_mmn_bnd(2)=bnd_nbr_ice
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_ice=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd(1)=mmn_nbr+1
  ! Array dimensions: bnd_nbr_ice
  allocate(abs_cff_mss_ice_dsk(bnd_nbr_ice),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_ice_dsk"
  allocate(sca_cff_mss_ice_dsk(bnd_nbr_ice),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_ice_dsk"
  allocate(asm_prm_ice_dsk(bnd_nbr_ice),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_ice_dsk"
  ! Array dimensions: bnd_nbr_ice+1
  allocate(wvl_grd_ice(bnd_nbr_ice+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_ice"
  ! Array dimensions: mmn_nbr,bnd_nbr_ice
  if (flg_mie_ice) then
     allocate(lgn_xpn_cff_ice_dsk(0:mmn_nbr,bnd_nbr_ice),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_ice_dsk"
  endif ! !flg_mie_ice
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_ice_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_ice_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_ice_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_ice_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_ice_id,abs_cff_mss_ice_dsk,srt_one,cnt_bnd),"gv abs_cff_mss_ice")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_ice_id,asm_prm_ice_dsk,srt_one,cnt_bnd),"gv asm_prm_ice")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_ice_id,sca_cff_mss_ice_dsk,srt_one,cnt_bnd),"gv sca_cff_mss_ice")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_ice_id,wvl_grd_ice,srt_one,cnt_bndp),"gv wvl_grd_ice")
  ! Get optional data
  if (flg_mie_ice) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_ice_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_ice_id,lgn_xpn_cff_ice_dsk,srt_one_one,cnt_mmn_bnd),"gv lgn_xpn_cff_ice")
     do bnd_idx=1,cnt_bnd(1)
        ! Sanity check on phase function expansion
        if (abs(lgn_xpn_cff_ice_dsk(1,bnd_idx)-asm_prm_ice_dsk(bnd_idx)) > 0.40) then
           write (6,'(3a,/,2(a,i4,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR inaccurate mie phase function expansion in ',fl_ice, &
                'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx,') = ',lgn_xpn_cff_ice_dsk(1,bnd_idx), &
                ' != asm_prm(bnd_idx=',bnd_idx,') = ',asm_prm_ice_dsk(bnd_idx), &
                '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
           call abort
        endif ! endif err
     enddo ! end loop over bnd
  end if ! !flg_mie_ice
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_ice,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_lqd
  rcd=nf90_wrp_open(fl_lqd,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_lqd),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_lqd>bnd_nbr_lqd_max) stop 'bnd_nbr_lqd>bnd_nbr_lqd_max'
  cnt_bnd(1)=bnd_nbr_lqd
  cnt_bndp(1)=bnd_nbr_lqd+1
  cnt_mmn_bnd(2)=bnd_nbr_lqd
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_lqd=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd(1)=mmn_nbr+1
  ! Array dimensions: bnd_nbr_lqd
  allocate(abs_cff_mss_lqd_dsk(bnd_nbr_lqd),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_lqd_dsk"
  allocate(sca_cff_mss_lqd_dsk(bnd_nbr_lqd),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_lqd_dsk"
  allocate(asm_prm_lqd_dsk(bnd_nbr_lqd),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_lqd_dsk"
  ! Array dimensions: bnd_nbr_lqd+1
  allocate(wvl_grd_lqd(bnd_nbr_lqd+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_lqd"
  ! Array dimensions: mmn_nbr,bnd_nbr_lqd
  if (flg_mie_lqd) then
     allocate(lgn_xpn_cff_lqd_dsk(0:mmn_nbr,bnd_nbr_lqd),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_lqd_dsk"
  endif ! !flg_mie_lqd
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_lqd_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_lqd_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_lqd_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_lqd_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_lqd_id,abs_cff_mss_lqd_dsk,srt_one,cnt_bnd),"gv abs_cff_mss_lqd")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_lqd_id,asm_prm_lqd_dsk,srt_one,cnt_bnd),"gv asm_prm_lqd")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_lqd_id,sca_cff_mss_lqd_dsk,srt_one,cnt_bnd),"gv sca_cff_mss_lqd")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_lqd_id,wvl_grd_lqd,srt_one,cnt_bndp),"gv wvl_grd_lqd")
  ! Get optional data
  if (flg_mie_lqd) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_lqd_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_lqd_id,lgn_xpn_cff_lqd_dsk,srt_one_one,cnt_mmn_bnd),"gv lgn_xpn_cff_lqd")
     do bnd_idx=1,cnt_bnd(1)
        ! Sanity check on phase function expansion
        if (abs(lgn_xpn_cff_lqd_dsk(1,bnd_idx)-asm_prm_lqd_dsk(bnd_idx)) > 0.40) then
           write (6,'(3a,/,2(a,i4,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR inaccurate mie phase function expansion in ',fl_lqd, &
                'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx,') = ',lgn_xpn_cff_lqd_dsk(1,bnd_idx), &
                ' != asm_prm(bnd_idx=',bnd_idx,') = ',asm_prm_lqd_dsk(bnd_idx), &
                '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
           call abort
        endif ! endif err
     enddo ! end loop over bnd
  end if ! !flg_mie_lqd
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_lqd,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_mpr
  rcd=nf90_wrp_open(fl_mpr,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_mpr),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_mpr>bnd_nbr_mpr_max) stop 'bnd_nbr_mpr>bnd_nbr_mpr_max'
  cnt_bnd(1)=bnd_nbr_mpr
  cnt_bndp(1)=bnd_nbr_mpr+1
  cnt_mmn_bnd(2)=bnd_nbr_mpr
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_mpr=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd(1)=mmn_nbr+1
  ! Get potentially missing dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_bnd_mpr_dmn_id,nf90_ebaddim)
  if (rcd == nf90_noerr) then
     ! Input file has multiple levels of optical properties
     rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_bnd_mpr_dmn_id,len=lev_bnd_mpr_nbr),sbr_nm//": inquire_dim lev")
     if (lev_bnd_mpr_nbr /= lev_snw_nbr) then
        stop 'Multi-layer impurity optical properties must agree with snow vertical grid'
     endif ! incommensurate vertical grids
     stop 'Input files with multi-layer impurity optical properties not yet supported'
  else
     ! Input file has single level of optical properties
     lev_bnd_mpr_nbr=1
  endif ! rcd != nf90_noerr
  ! Array dimensions: bnd_nbr_mpr
  allocate(abs_cff_mss_mpr_dsk(bnd_nbr_mpr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_mpr_dsk"
  allocate(sca_cff_mss_mpr_dsk(bnd_nbr_mpr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_mpr_dsk"
  allocate(asm_prm_mpr_dsk(bnd_nbr_mpr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_mpr_dsk"
  ! Array dimensions: bnd_nbr_mpr+1
  allocate(wvl_grd_mpr(bnd_nbr_mpr+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_mpr"
  ! Array dimensions: mmn_nbr,bnd_nbr_mpr
  if (flg_mie_mpr) then
     allocate(lgn_xpn_cff_mpr_dsk(0:mmn_nbr,bnd_nbr_mpr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_mpr_dsk"
  endif ! !flg_mie_mpr
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_mpr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_mpr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_mpr_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_mpr_id)
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_mpr_id,abs_cff_mss_mpr_dsk,srt_one,cnt_bnd),"gv abs_cff_mss_mpr")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_mpr_id,asm_prm_mpr_dsk,srt_one,cnt_bnd),"gv asm_prm_mpr")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_mpr_id,sca_cff_mss_mpr_dsk,srt_one,cnt_bnd),"gv sca_cff_mss_mpr")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_mpr_id,wvl_grd_mpr,srt_one,cnt_bndp),"gv wvl_grd_mpr")
  ! Get optional data
  if (flg_mie_mpr) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_mpr_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_mpr_id,lgn_xpn_cff_mpr_dsk,srt_one_one,cnt_mmn_bnd),"gv lgn_xpn_cff_mpr")
     do bnd_idx=1,cnt_bnd(1)
        ! Sanity check on phase function expansion
        if (abs(lgn_xpn_cff_mpr_dsk(1,bnd_idx)-asm_prm_mpr_dsk(bnd_idx)) > 0.40) then
           write (6,'(3a,/,2(a,i4,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR inaccurate mie phase function expansion in ',fl_mpr, &
                'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx,') = ',lgn_xpn_cff_mpr_dsk(1,bnd_idx), &
                ' != asm_prm(bnd_idx=',bnd_idx,') = ',asm_prm_mpr_dsk(bnd_idx), &
                '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
           call abort
        endif ! endif err
     enddo ! end loop over bnd
  end if ! !flg_mie_mpr
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_mpr,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_snw
  rcd=nf90_wrp_open(fl_snw,nf90_nowrite,nc_id)
  ! Get required dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_snw),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_snw>bnd_nbr_snw_max) stop 'bnd_nbr_snw>bnd_nbr_snw_max'
  cnt_bnd_lev(1)=bnd_nbr_snw
  cnt_bndp(1)=bnd_nbr_snw+1
  cnt_mmn_bnd_lev(2)=bnd_nbr_snw
  ! Get potentially missing dimension IDs and sizes
  if (flg_mie) then
     flg_mie_snw=.true. ! [flg] Require phase function expansion
     rcd=nf90_wrp_inq_dimid(nc_id,'lgn',mmn_spc_dmn_id,nf90_ebaddim)
     if (rcd == nf90_noerr) then
        ! Input file has phase function moments
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,mmn_spc_dmn_id,len=lgn_dmn_sz),sbr_nm//": inquire_dim lgn")
        ! Dimension size on disk is one more than number of Legendre moments
        mmn_spc_nbr=lgn_dmn_sz-1
        if (mmn_spc_nbr < mmn_nbr) then
           stop 'Species has fewer phase function moments than RT streams'
        endif ! Incommensurate expansions
     else
        ! Input file lacks phase function moments
        stop 'Species optical properties file lacks phase function moments'
     endif ! rcd != nf90_noerr
  endif ! !flg_mie
  cnt_mmn_bnd_lev(1)=mmn_nbr+1
  
  ! Get potentially missing dimension IDs and sizes
  rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_bnd_snw_dmn_id,nf90_ebaddim)
  if (rcd == nf90_noerr) then
     ! Input file has multiple levels of optical properties
     rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_bnd_snw_dmn_id,len=lev_bnd_snw_nbr),sbr_nm//": inquire_dim lev")
     if (lev_bnd_snw_nbr /= lev_snw_nbr) then
        stop 'Multi-layer snow optical properties must agree with snow vertical grid'
     endif ! Incommensurate vertical grids
  else
     ! Input file has single level of optical properties
     lev_bnd_snw_nbr=1
  endif ! rcd != nf90_noerr
  cnt_bnd_lev(2)=lev_bnd_snw_nbr
  cnt_mmn_bnd_lev(3)=lev_bnd_snw_nbr
  
  ! Array dimensions: bnd_nbr_snw,lev_bnd_snw_nbr
  allocate(abs_cff_mss_snw_dsk(bnd_nbr_snw,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_snw_dsk"
  allocate(sca_cff_mss_snw_dsk(bnd_nbr_snw,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_snw_dsk"
  allocate(asm_prm_snw_dsk(bnd_nbr_snw,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_snw_dsk"
  ! Array dimensions: bnd_nbr_snw+1
  allocate(wvl_grd_snw(bnd_nbr_snw+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_snw"
  ! Array dimensions: mmn_nbr,bnd_nbr_snw,lev_bnd_snw_nbr
  if (flg_mie_snw) then
     allocate(lgn_xpn_cff_snw_dsk(0:mmn_nbr,bnd_nbr_snw,lev_bnd_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_snw_dsk"
  endif ! !flg_mie_snw

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_snw_id)
  rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_snw_id)
  rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_snw_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_snw_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_snw_id,abs_cff_mss_snw_dsk,srt_one_one,cnt_bnd_lev),"gv abs_cff_mss_snw")
  rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_snw_id,sca_cff_mss_snw_dsk,srt_one_one,cnt_bnd_lev),"gv sca_cff_mss_snw")
  rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_snw_id,asm_prm_snw_dsk,srt_one_one,cnt_bnd_lev),"gv asm_prm_snw")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_snw_id,wvl_grd_snw,srt_one,cnt_bndp),"gv wvl_grd_snw")
  ! Get optional data
  if (flg_mie_snw) then
     rcd=nf90_wrp_inq_varid(nc_id,'lgn_xpn_cff',lgn_xpn_cff_snw_id)
     rcd=nf90_wrp(nf90_get_var(nc_id,lgn_xpn_cff_snw_id,lgn_xpn_cff_snw_dsk,srt_one_one_one,cnt_mmn_bnd_lev),"gv lgn_xpn_cff_snw")
     do bnd_idx=1,cnt_mmn_bnd_lev(2)
        do lev_idx=1,cnt_mmn_bnd_lev(3)
           ! Sanity check on phase function expansion
           if (abs(lgn_xpn_cff_snw_dsk(1,bnd_idx,lev_idx)-asm_prm_snw_dsk(bnd_idx,lev_idx)) > 0.40) then
              write (6,'(3a,/,2(a,i4,a,i2,a,e9.2),/,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                   ': ERROR inaccurate mie phase function expansion in ',fl_snw, &
                   'lgn_xpn_cff(mmn_idx=1,bnd_idx=',bnd_idx, &
                   ',lev_idx=',lev_idx,') = ', &
                   lgn_xpn_cff_snw_dsk(1,bnd_idx,lev_idx), &
                   ' != asm_prm(bnd_idx=',bnd_idx, &
                   ',lev_idx=',lev_idx,') = ', &
                   asm_prm_snw_dsk(bnd_idx,lev_idx), &
                   '. HINT: Re-generate file or use Henyey-Greenstein approximation.'
              call abort
           endif ! endif err
        enddo ! end loop over lev
     enddo ! end loop over bnd
  end if ! !flg_mie_snw
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_snw,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_rfl
  rcd=nf90_wrp_open(fl_rfl,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_rfl),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_rfl>bnd_nbr_rfl_max) stop 'bnd_nbr_rfl>bnd_nbr_rfl_max'
  cnt_bnd(1)=bnd_nbr_rfl
  cnt_bndp(1)=bnd_nbr_rfl+1
  
  ! Array dimensions: bnd_nbr_rfl
  allocate(rfl_spc_sfc_dsk(bnd_nbr_rfl),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for rfl_spc_sfc_dsk"
  allocate(wvl_grd_rfl(bnd_nbr_rfl+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_rfl"

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'rfl_spc_sfc',rfl_spc_sfc_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_rfl_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,rfl_spc_sfc_id,rfl_spc_sfc_dsk,srt_one,cnt_bnd),"gv rfl_spc_sfc_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_rfl_id,wvl_grd_rfl,srt_one,cnt_bndp),"gv wvl_grd")

  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_rfl,'Ingested') ! [fnc] Close file

  ! Ingest fl_lmn
  rcd=nf90_wrp_open(fl_lmn,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_lmn),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_lmn>bnd_nbr_lmn_max) stop 'bnd_nbr_lmn>bnd_nbr_lmn_max'
  cnt_bnd(1)=bnd_nbr_lmn
  cnt_bndp(1)=bnd_nbr_lmn+1
  
  ! Array dimensions: bnd_nbr_lmn
  allocate(lmn_SRF_dsk(bnd_nbr_lmn),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_SRF_dsk"
  allocate(wvl_grd_lmn(bnd_nbr_lmn+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd_lmn"

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'lmn_SRF',lmn_SRF_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_lmn_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,lmn_SRF_id,lmn_SRF_dsk,srt_one,cnt_bnd),"gv lmn_SRF_dsk")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_lmn_id,wvl_grd_lmn,srt_one,cnt_bndp),"gv wvl_grd")

  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_lmn,'Ingested') ! [fnc] Close file
  
  ! Ingest fl_nst
  rcd=nf90_wrp_open(fl_nst,nf90_nowrite,nc_id)
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dmn_id)
  
  ! Get dimension sizes
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_nst),sbr_nm//": inquire_dim bnd")
  if (bnd_nbr_nst>bnd_nbr_nst_max) stop 'bnd_nbr_nst>bnd_nbr_nst_max'
  cnt_bnd(1)=bnd_nbr_nst
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,'nst_SRF',nst_SRF_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_nst_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_nst_id)
  rcd=nf90_wrp_inq_varid(nc_id,'wvl',wvl_ctr_nst_id)
  
  ! Get data
  rcd=nf90_wrp(nf90_get_var(nc_id,nst_SRF_id,nst_SRF,srt_one,cnt_bnd),"gv nst_SRF")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_nst_id,wvl_max_nst,srt_one,cnt_bnd),"gv wvl_max_nst")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_nst_id,wvl_min_nst,srt_one,cnt_bnd),"gv wvl_min_nst")
  rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_nst_id,wvl_ctr_nst,srt_one,cnt_bnd),"gv wvl_ctr_nst")
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_nst,'Ingested') ! [fnc] Close file
  
  if (mode_std) then
     ! Do not read in ocean mask and wind speeds
     ! Set channel dimension size to one for convenience/safety
     chn_nbr=1
  else 
     ! Ingest fl_chn
     rcd=nf90_wrp_open(fl_chn,nf90_nowrite,nc_id)
     ! Get dimension IDs
     ! rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
     rcd=nf90_wrp_inq_dimid(nc_id,'chn',chn_dmn_id)
     
     ! Get dimension sizes
     ! rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_nst),sbr_nm//": inquire_dim bnd")
     ! if (bnd_nbr_nst>bnd_nbr_nst_max) stop 'bnd_nbr_nst>bnd_nbr_nst_max'
     rcd=nf90_wrp(nf90_inquire_dimension(nc_id,chn_dmn_id,len=chn_nbr),sbr_nm//": inquire_dim chn")
     if (chn_nbr>chn_nbr_max) stop 'chn_nbr>chn_nbr_max'
     cnt_chn(1)=chn_nbr
     
     ! Get variable IDs
     rcd=nf90_wrp_inq_varid(nc_id,'chn_SRF',chn_SRF_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_nst_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_nst_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_nst_id)
     
     ! Get data
     rcd=nf90_wrp(nf90_get_var(nc_id,chn_SRF_id,chn_SRF,srt_one,cnt_chn),"gv chn_SRF")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_nst_id,wvl_max_nst,srt_one,cnt_chn),"gv wvl_max_nst")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_nst_id,wvl_min_nst,srt_one,cnt_chn),"gv wvl_min_nst")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_nst_id,wvl_ctr_nst,srt_one,cnt_chn),"gv wvl_ctr_nst")
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_chn,'Ingested') ! [fnc] Close file
  endif ! !mode_std
  
  ! if (.not.lamber) then
  if (.false.) then         ! Do not read brdf file for now
     ! Ingest fl_brdf ! BRDF
     rcd=nf90_wrp_open(fl_brdf,nf90_nowrite,nc_id)
     ! Get dimension IDs
     rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
     ! Get dimension sizes
     rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr_brdf),sbr_nm//": inquire_dim bnd")
     if (bnd_nbr_brdf>bnd_nbr_brdf_max) stop 'bnd_nbr_brdf>bnd_nbr_brdf_max'
     cnt_bnd(1)=bnd_nbr_brdf
     ! Get variable IDs
     rcd=nf90_wrp_inq_varid(nc_id,'brdf_type',brdf_typ_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_brdf_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_brdf_id)
     rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_brdf_id)
     rcd=nf90_wrp_inq_varid(nc_id,'f_iso_spc',f_iso_spc_id)
     rcd=nf90_wrp_inq_varid(nc_id,'f_vol_spc',f_vol_spc_id)
     rcd=nf90_wrp_inq_varid(nc_id,'f_geo_spc',f_geo_spc_id)
     ! Get data
     rcd=nf90_wrp(nf90_get_var(nc_id,brdf_typ_id,brdf_typ),"gv brdf_typ")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_brdf_id,wvl_max_brdf,srt_one,cnt_bnd),"gv wvl_max_brdf")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_brdf_id,wvl_min_brdf,srt_one,cnt_bnd),"gv wvl_min_brdf")
     rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_brdf_id,wvl_ctr_brdf,srt_one,cnt_bnd),"gv wvl_ctr_brdf")
     rcd=nf90_wrp(nf90_get_var(nc_id,f_iso_spc_id,f_iso_spc,srt_one,cnt_bnd),"gv f_iso_spc")
     rcd=nf90_wrp(nf90_get_var(nc_id,f_vol_spc_id,f_vol_spc,srt_one,cnt_bnd),"gv f_vol_spc")
     rcd=nf90_wrp(nf90_get_var(nc_id,f_geo_spc_id,f_geo_spc,srt_one,cnt_bnd),"gv f_geo_spc")
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_brdf,'Ingested') ! [fnc] Close file
  endif                     ! end if not lamber
  
  if (.not.mode_std) then
     if (ocn_msk < 1.0) then
        print *, 'ocn_msk < 1.0', ocn_msk
        stop
     endif
  endif
  
  ! All necessary input data have been read

  ! Allocate arrays whose dimension sizes depend on input data and/or command line switches
  ! Array dimensions: azi
  allocate(azi(azi_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for azi"
  allocate(azi_dgr(azi_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for azi_dgr"
  ! Array dimensions: bnd (still unknown, see below)
  ! Array dimensions: chn
  allocate(flx_chn_dwn_TOA(chn_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_chn_dwn_TOA"
  allocate(flx_chn_upw_TOA(chn_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_chn_upw_TOA"
  ! Array dimensions: grd (still unknown, see below)
  ! Array dimensions: lev
  allocate(flx_bb_abs(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_abs"
  allocate(flx_nst_abs(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_nst_abs"
  allocate(htg_rate_bb(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for htg_rate_bb"
  allocate(j_NO2(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for j_NO2"
  allocate(ntn_bb_mean(lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_bb_mean"
  ! Array dimensions: levp
  allocate(flx_bb_dwn(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_dwn"
  allocate(flx_bb_dwn_dff(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_dwn_dff"
  allocate(flx_bb_dwn_drc(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_dwn_drc"
  allocate(flx_bb_net(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_net"
  allocate(flx_bb_upw(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_bb_upw"
  allocate(flx_nst_dwn(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_nst_dwn"
  allocate(flx_nst_net(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_nst_net"
  allocate(flx_nst_upw(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_nst_upw"
  allocate(ilm_dwn(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ilm_dwn"
  allocate(ilm_upw(levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ilm_upw"
  ! Array dimensions: plr
  allocate(plr(plr_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for plr"
  allocate(plr_cos(plr_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for plr_cos"
  allocate(plr_dgr(plr_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for plr_dgr"
  ! Array dimensions: tau (still unknown, see below)
  ! Array dimensions: bnd,lev (still unknown, see below)
  ! Array dimensions: bnd,levp (still unknown, see below)
  ! Array dimensions: plr,bnd (still unknown, see below)
  ! Array dimensions: plr,levp
  allocate(lmn_bb_aa(plr_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_bb_aa"
  allocate(ntn_bb_aa(plr_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_bb_aa"
  ! Array dimensions: azi,plr,bnd (still unknown, see below)
  ! Array dimensions: azi,plr,chn
  allocate(ntn_chn_TOA(azi_nbr,plr_nbr,chn_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_chn_TOA"
  allocate(rfl_chn_TOA(azi_nbr,plr_nbr,chn_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for rfl_chn_TOA"
  ! Array dimensions: azi,plr,levp
  allocate(ntn_spc_chn(azi_nbr,plr_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_chn"

  ! Initialize input data that may overridden by command line switches
  if (cmd_ln_odxc_obs_aer) then
     aer_lvl_nbr=0
     do lev_idx=1,lev_nbr
        if (odxl_obs_aer(lev_idx) > 0.0) aer_lvl_nbr=aer_lvl_nbr+1
     enddo                  ! end loop over lev
     if (aer_lvl_nbr == 0) stop 'aer_lvl_nbr==0'
     do lev_idx=1,lev_nbr
        if (odxl_obs_aer(lev_idx) > 0.0) odxl_obs_aer(lev_idx)=odxc_obs_aer_cmd_ln/real(aer_lvl_nbr)
     enddo                  ! end loop over lev
     ! Recompute column path
     if (odxc_obs_aer_cmd_ln /= odxc_obs_aer) then
        write (6,'(a,a,f9.6,a,f9.6,a,f9.6,a)') prg_nm(1:ftn_strlen(prg_nm)), &
             ': INFO User-specified optical depth of ',odxc_obs_aer_cmd_ln, &
             ' differs from optical depth in input profile of ',odxc_obs_aer, &
             ' (both taken to be at ',wvl_obs_aer*1.0e6,' um)'
        write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
             ': INFO Expect WARNING messages as aerosol mass paths are changed'
     endif                  ! endif
     odxc_obs_aer=0.0
     do lev_idx=1,lev_nbr
        odxc_obs_aer=odxc_obs_aer+odxl_obs_aer(lev_idx)
     enddo                  ! end loop over lev
  endif                     ! end if overriding CLM profile odxc_obs_aer
  
  if (cmd_ln_mpc_CWP) then
     cld_lvl_nbr=0
     do lev_idx=1,lev_nbr
        if (mpl_CWP(lev_idx)>0.0) cld_lvl_nbr=cld_lvl_nbr+1
     enddo                  ! end loop over lev
     if (cld_lvl_nbr==0 .and. mpc_CWP_cmd_ln/=0.0) stop 'cld_lvl_nbr==0'
     do lev_idx=1,lev_nbr
        if (mpl_CWP(lev_idx)>0.0) mpl_CWP(lev_idx)=mpc_CWP_cmd_ln/real(cld_lvl_nbr)
     enddo                  ! end loop over lev
  endif                     ! end if overriding CLM profile mpl_CWP
  if (force_ice_phz) then
     do lev_idx=1,lev_nbr
        frc_ice(lev_idx)=1.0
     enddo                  ! end loop over lev
  endif                     ! end if forcing ice phase mpl_CWP
  if (force_lqd_phz) then
     do lev_idx=1,lev_nbr
        frc_ice(lev_idx)=0.0
     enddo                  ! end loop over lev
  endif                     ! end if forcing liquid phase mpl_CWP
  if (cmd_ln_mpc_CWP.or.force_ice_phz.or.force_lqd_phz) then
     mpc_CWP=0.0
     mpc_IWP=0.0
     do lev_idx=1,lev_nbr
        mpl_IWP(lev_idx)=frc_ice(lev_idx)*mpl_CWP(lev_idx)
        mpl_LWP(lev_idx)=max(0.0,mpl_CWP(lev_idx)-mpl_IWP(lev_idx))
        mpc_CWP=mpc_CWP+mpl_CWP(lev_idx)
        mpc_IWP=mpc_IWP+mpl_IWP(lev_idx)
     enddo                  ! end loop over lev
     if (mpc_CWP /= 0.0) then
        frc_ice_ttl=mpc_IWP/mpc_CWP
     else                   ! end if column is cloudy
        frc_ice_ttl=0.0
        alt_cld_btm=0.0
        alt_cld_thick=0.0
     endif                  ! end if column is clear
  endif                     ! end if recomputing mpl_CWP, mpl_IWP, mpl_LWP
  ! flg_cld_sat and flg_sat_cld must not both be true
  if (flg_cld_sat) flg_sat_cld=.false.
  if (flg_cld_sat) then
     ! Force cloudy and snowy layers to be saturated
     do lev_idx=1,lev_nbr
        if (mpl_CWP(lev_idx)>0.0) then
           ! fxm: Propogate this correction to q_H2OH2O as well
           sat_vpr_lqd=svp_H2O_lqd_PrK78(tpt(lev_idx))
           qst_H2O_lqd=eps_H2O*sat_vpr_lqd/(prs(lev_idx)-one_mns_eps_H2O*sat_vpr_lqd)
           q_H2O(lev_idx)=qst_H2O_lqd
        endif ! endif cloudy layer
     enddo                  ! end loop over lev
     ! Force snowy layers to be saturated with respect to ice
     ! fxm: not yet implemented because program architecture makes it hard
     ! swnb2 does not use q_H2O. swnb2 relies on/extrapolates from mpl_H2O.
     ! clm does not yet keep track of q_H2O in snow grid
     if (flg_cld_sat) then
        do lev_snw_idx=1,lev_snw_nbr
           sat_vpr_ice=svp_H2O_ice_PrK78(tpt_snw(lev_snw_idx))
           qst_H2O_ice=eps_H2O*sat_vpr_ice/(prs(lev_atm_nbr+lev_snw_idx)-one_mns_eps_H2O*sat_vpr_ice)
           
           ! q_H2O(lev_snw_idx)=qst_H2O_ice
        enddo                  ! end loop over lev_snw
     endif
  endif ! endif flg_cld_sat
  if (flg_sat_cld) then
     ! Force saturated layers to be cloudy (currently implemented by clm not swnb2)
  endif ! endif flg_sat_cld
  if (cmd_ln_lcl_yr_day) then
     ! Derive solar geometry based on command-line input year of day
     ! GMT Year-day [1.0..366.0) (365-day GCM calendar) is processed as in CCM3/Bri92
     ! True local sidereal noon in this framework is julian day fraction 0.5
     ! Solar coordinates are "climatological" late 20th century
     ! The "exact" solution, like Mic88 used in CLM, requires UNIX time, etc.
     ! Latitude is read from command-line or from CLM input file
     ! Longitude is implicitly zero (because year-day is GMT for simplicity)
     lcl_yr_day=lcl_yr_day_cmd_ln
     lcl_time_hr=(lcl_yr_day-real(int(lcl_yr_day)))*24.0 ! [hr] Local time hour
     lon_dgr=0.0            ! [dgr]
     ! Command-line lat_dgr only has effect if lcl_yr_day also command-line specified
     if (cmd_ln_lat_dgr) then
        lat_dgr=lat_dgr_cmd_ln
     endif ! end if overriding CLM profile latitude
     lat=pi*lat_dgr/180.0   ! [dgr] -> [rdn]
     call slr_crd_Bri92( &
          lat,                 & ! I [rdn] Latitude
          lcl_yr_day,          & ! I [day] Local year day
          slr_zen_ngl_cos,     & ! O [frc] Solar zenith angle cosine
          xnt_fac)             ! O [frc] Eccentricity factor
     if (dbg_lvl>=dbg_scl) then
        write (6,'(a,f15.12)') 'xnt_fac from Bri92 = ',xnt_fac
     endif ! endif dbg
     ! Set eccentricity factor to one for consistency with climatological assumptions
     xnt_fac=1.0
  endif ! end if overriding CLM profile local year day
  if (cmd_ln_slr_zen_ngl_cos) then
     slr_zen_ngl_cos=slr_zen_ngl_cos_cmd_ln
     ! Set eccentricity factor to one to ease intercomparison with other RT models
     xnt_fac=1.0
  endif                     ! end if overriding CLM profile zen ang
  if (dbg_lvl>=dbg_scl) then
     write (6,'(a,f15.12)') 'slr_zen_ngl_cos = ',slr_zen_ngl_cos
  endif ! endif dbg
  if (cmd_ln_alb) then
     alb_sfc_NIR_dff=alb_cmd_ln
     alb_sfc_NIR_drc=alb_cmd_ln
     alb_sfc_vsb_dff=alb_cmd_ln
     alb_sfc_vsb_drc=alb_cmd_ln
     alb_sfc=alb_cmd_ln
  else
     alb_sfc=0.5*(alb_sfc_vsb_drc+alb_sfc_NIR_drc)
  endif                     ! end if overriding CLM profile albedo
  ! Sanity check
  if (alb_sfc>1.0.or.alb_sfc<0.0) stop 'alb_sfc>1.0.or.alb_sfc<0.0 in swnb2()'
  if (cmd_ln_slr_cst) then
     slr_cst=slr_cst_cmd_ln
  endif                     ! end if overriding solar constant
  
  ! There is no easier way to turn off Herzberg continuum
  ! but keep O2 line absorption on than to zero absorption cross-sections here.
  if (.not.flg_Herzberg) then
     do bnd_idx_O3=1,bnd_nbr_O3
        abs_xsx_O2(bnd_idx_O3)=0.0
     enddo                  ! end loop over O3 bands
  endif                     ! end if turning off Herzberg continuum
  
  ! Many DISORT arguments only need to be set once
  ! Set them here, outside band loop
  nstr=str_nbr
  nlyr=lev_nbr
  nmom=mmn_nbr              ! for DISORT2
  
  ! Temperatures are specified at interfaces
  do levp_idx=0,levp_nbr-1
     temper(levp_idx)=tpt_ntf(levp_idx+1)
  end do
  
  ! NB: Setting usrtau to .true. has profound consequences on 
  ! array lengths returned by DISORT.
  ! It is "safest" (easiest) to keep usrtau = .false. 
  ! Then DISORT returns all radiant quantities at computational layer interfaces
  usrtau=.not.sv_cmp_tau
  if (sv_cmp_tau) then 
     tau_nbr=levp_nbr
  else
     tau_nbr=levp_nbr
  endif                     ! end if no user defined levels
  ! Array dimensions: tau
  allocate(tau(tau_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tau"
  allocate(tau_prs(tau_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for tau_prs"
  ntau=tau_nbr
  do tau_idx=1,tau_nbr
     tau_prs(tau_idx)=prs_ntf(tau_idx)
  end do                    ! end loop over tau
  
  ! Set angles at which to report intensities
  ! Array must be in increasing order of cosine polar angle
  usrang=.not.sv_cmp_plr_ngl
  numu=plr_nbr
  if (usrang) then
     plr(1)=pi
     do plr_idx=2,plr_nbr-1
        plr(plr_idx)=pi-real((plr_idx-1))*pi/real((plr_nbr-1))
     enddo                  ! end loop over plr
     plr(plr_nbr)=0.0
     
     do plr_idx=1,plr_nbr
        plr_dgr(plr_idx)=180.0*plr(plr_idx)/pi
     enddo                  ! end loop over plr
     
     plr_cos(1)=-1.0        ! 180 degrees, i.e., towards nadir
     do plr_idx=2,plr_nbr-1
        plr_cos(plr_idx)=cos(plr(plr_idx))
     enddo                  ! end loop over plr
     plr_cos(plr_nbr)=1.0   ! 0 degrees, i.e., towards zenith
     
     do plr_idx=1,plr_nbr
        umu(plr_idx)=plr_cos(plr_idx)
     enddo                  ! end loop over plr
  endif                     ! endif if looking at user angles
  
  ! Set azimuthal angles at which to report intensities.
  ! nphi=0 is valid only when onlyfl=.true.
  nphi=azi_nbr
  do azi_idx=1,azi_nbr
     azi(azi_idx)=real((azi_idx-1))*2.0*pi/real(azi_nbr)
     azi_dgr(azi_idx)=real((azi_idx-1))*360.0/real(azi_nbr)
     phi(azi_idx)=azi_dgr(azi_idx)
  enddo                     ! end loop over azimuthal angles
  
  ! Set boundary conditions
  ! 0: General case, includes beam illumination from top
  ibcnd=0
  
  ! Set solar zenith angle cosine
  ! Beware of using 0.5 when str_nbr/2 is odd---
  ! it is a quadrature point and can crash DISORT.
  umu0=slr_zen_ngl_cos
  
  ! Set azimuth angle of incident sunlight
  phi0=0.0
  
  ! Specify bottom boundary temperature btemp
  ! Btemp for atmosphere is skin temperature in CCM/CAM parlance
  ! Bottom boundary emissivity is derived from albedo above
  ! Top boundary temperature and emissivity must also be specified
  ! Of course these are only used when plank is .true.
  if (flg_msm) then 
     btemp=tpt_ntf_snw(levp_snw_nbr)
  else ! !flg_msm
     btemp=tpt_skn 
  endif ! !flg_msm
  ttemp=tpt_ntf(1)
  temis=0.0 ! 20160513: Emissivity of upper boundary is usually 0.0
  
  ! Set control flags:
  ! Set deltam=.true. unless looking at radiances within 10 degrees of forward peak
  ! When deltam=.true., the returned downwelling direct and diffuse fluxes are the "true" fluxes which have been recovered from scaled fluxes
  ! Upwelling flux is (theoretically) not affected by delta scaling
  ! Thus the direct/diffuse ratio computed with deltam=.true. approximately equals the ratio computed with deltam=.false.
  ! Only forward beam radiances change significantly when deltam=.true.
  !deltam=.true.
  
  ! Sanity check
  if (slr_zen_ngl_cos<=0.0.and..not.flg_Planck) then 
     write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': WARNING Sun beneath horizon for pure solar calculation'
  endif
  
  ! Set onlyfl=.false. when looking at user angle radiances
  ! When onlyfl=.true., azimuthally averaged radiances at 
  ! computational angles--NOT user angles--will be reported 
  ! in u0u whenever maxumu >= nstr.
  onlyfl=.false.
  
  ! accur is maximum relative error in last three terms 
  ! in azimuthal series, and it determines convergence
  ! accur value does not seem to contribute to any single precision problems 
  ! Set 0.0 < accur < 0.01
  accur=0.0
  
  ! Set printing flags for DISORT:
  prnt(1)=.false.           ! for DISORT2
  prnt(2)=.false.
  prnt(3)=.false.
  prnt(4)=.false.
  prnt(5)=.false.
  
  ! Band-independent DISORT() initialization is now complete
  ! Remaining DISORT arguments need to be set inside main band loop
  
  ! Compute how many total bands will be used in entire band
  ! computation by adding the number of narrow band H2O
  ! data to the O3-O2 continua. Number of pure O2-O3 bands
  ! will be number of bands from the O2-O3 continua data
  ! truncated at shortest wavelength narrow band data.
  do bnd_idx_O3=1,bnd_nbr_O3
     if (wvl_max_O3(bnd_idx_O3)>wvl_min_H2O(bnd_nbr_H2O)) then
        bnd_nbr_pure_O3=bnd_idx_O3
        goto 100
     endif
  enddo                     ! end loop over O3 bands
100 continue
  do bnd_idx_H2O=1,bnd_nbr_H2O
     if (wvl_min_H2O(bnd_idx_H2O)<wvl_max_O3(bnd_nbr_O3)) then
        bnd_nbr_non_O3=bnd_idx_H2O-1
        goto 110
     endif
  enddo                     ! end loop over H2O bands
110 continue
  
  bnd_nbr=bnd_nbr_pure_O3+bnd_nbr_H2O

  ! Array dimensions: bnd
  allocate(ext_cff_mss_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ext_cff_mss_aer"
  allocate(abs_cff_mss_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_aer"
  allocate(sca_cff_mss_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_aer"
  allocate(asm_prm_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_aer"
  allocate(abs_cff_mss_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_bga"
  allocate(sca_cff_mss_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_bga"
  allocate(asm_prm_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_bga"
  allocate(abs_cff_mss_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_ice"
  allocate(sca_cff_mss_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_ice"
  allocate(asm_prm_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_ice"
  allocate(abs_cff_mss_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_lqd"
  allocate(sca_cff_mss_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_lqd"
  allocate(asm_prm_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_lqd"
  allocate(abs_cff_mss_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_mpr"
  allocate(sca_cff_mss_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_mpr"
  allocate(asm_prm_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_mpr"
  allocate(abs_spc_SAS(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_spc_SAS"
  allocate(abs_spc_atm(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_spc_atm"
  allocate(abs_spc_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_spc_sfc"
  allocate(bnd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for bnd"
  allocate(flx_abs_atm_rdr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_abs_atm_rdr"
  allocate(flx_frc_dwn_sfc_blr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_frc_dwn_sfc_blr"
  allocate(flx_frc_dwn_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_frc_dwn_sfc"
  allocate(flx_slr_frc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_slr_frc"
  allocate(flx_spc_abs_SAS(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_abs_SAS"
  allocate(flx_spc_abs_atm(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_abs_atm"
  allocate(flx_spc_abs_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_abs_sfc"
  allocate(flx_spc_act_pht_TOA(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_act_pht_TOA"
  allocate(flx_spc_act_pht_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_act_pht_sfc"
  allocate(flx_spc_dwn_TOA(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_TOA"
  allocate(flx_spc_dwn_dff_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_dff_sfc"
  allocate(flx_spc_dwn_drc_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_drc_sfc"
  allocate(flx_spc_dwn_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_sfc"
  allocate(flx_spc_pht_dwn_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_pht_dwn_sfc"
  allocate(j_spc_NO2_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for j_spc_NO2_sfc"
  allocate(lmn_SRF(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_SRF"
  allocate(nrg_pht(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for nrg_pht"
  allocate(lmn_spc_aa_ndr_TOA(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_spc_aa_ndr_TOA"
  allocate(lmn_spc_aa_ndr_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_spc_aa_ndr_sfc"
  allocate(ntn_spc_aa_ndr_TOA(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_ndr_TOA"
  allocate(ntn_spc_aa_ndr_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_ndr_sfc"
  allocate(ntn_spc_aa_zen_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_zen_sfc"
  allocate(odac_spc_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_aer"
  allocate(odac_spc_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_mpr"
  allocate(odac_spc_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_snw"
  allocate(odac_spc_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_bga"
  allocate(odac_spc_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_ice"
  allocate(odac_spc_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odac_spc_lqd"
  allocate(odxc_spc_CO2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_CO2"
  allocate(odxc_spc_H2O(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_H2O"
  allocate(odxc_spc_H2OH2O(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_H2OH2O"
  allocate(odxc_spc_NO2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_NO2"
  allocate(odxc_spc_O2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_O2"
  allocate(odxc_spc_O2N2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_O2N2"
  allocate(odxc_spc_O2O2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_O2O2"
  allocate(odxc_spc_O3(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_O3"
  allocate(odxc_spc_OH(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_OH"
  allocate(odxc_spc_CH4(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_CH4"
  allocate(odxc_spc_Ray(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_Ray"
  allocate(odxc_spc_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_aer"
  allocate(odxc_spc_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_mpr"
  allocate(odxc_spc_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_snw"
  allocate(odxc_spc_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_bga"
  allocate(odxc_spc_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_ice"
  allocate(odxc_spc_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_lqd"
  allocate(odxc_spc_ttl(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxc_spc_ttl"
  allocate(rfl_spc_SAS(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for rfl_spc_SAS"
  allocate(alb_spc_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for alb_spc_snw"
  allocate(flx_spc_abs_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_abs_snw"
  allocate(flx_spc_dwn_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_snw"
  allocate(flx_spc_upw_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_upw_snw"
  allocate(rfl_spc_sfc(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for rfl_spc_sfc"
  allocate(trn_spc_atm_CO2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_CO2"
  allocate(trn_spc_atm_H2O(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_H2O"
  allocate(trn_spc_atm_H2OH2O(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_H2OH2O"
  allocate(trn_spc_atm_NO2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_NO2"
  allocate(trn_spc_atm_O2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_O2"
  allocate(trn_spc_atm_O2N2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_O2N2"
  allocate(trn_spc_atm_O2O2(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_O2O2"
  allocate(trn_spc_atm_O3(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_O3"
  allocate(trn_spc_atm_OH(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_OH"
  allocate(trn_spc_atm_CH4(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_CH4"
  allocate(trn_spc_atm_Ray(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_Ray"
  allocate(trn_spc_atm_aer(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_aer"
  allocate(trn_spc_atm_mpr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_mpr"
  allocate(trn_spc_atm_snw(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_snw"
  allocate(trn_spc_atm_bga(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_bga"
  allocate(trn_spc_atm_ice(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_ice"
  allocate(trn_spc_atm_lqd(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_lqd"
  allocate(trn_spc_atm_ttl(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for trn_spc_atm_ttl"
  allocate(wvl(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl"
  allocate(wvl_ctr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_ctr"
  allocate(wvl_dlt(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_dlt"
  allocate(wvl_max(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_max"
  allocate(wvl_min(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_min"
  allocate(wvn(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvn"
  allocate(wvn_ctr(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvn_ctr"
  allocate(wvn_dlt(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvn_dlt"
  allocate(wvn_max(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvn_max"
  allocate(wvn_min(bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvn_min"
  
  ! Array dimensions: grd
  allocate(wvl_grd(bnd_nbr+1),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for wvl_grd"

  ! Array dimensions: bnd,lev_bnd_snw
  allocate(abs_cff_mss_snw(bnd_nbr,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for abs_cff_mss_snw"
  allocate(sca_cff_mss_snw(bnd_nbr,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for sca_cff_mss_snw"
  allocate(asm_prm_snw(bnd_nbr,lev_bnd_snw_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_snw"

  if (flg_mie) then
     ! Array dimensions: mmn,bnd
     allocate(lgn_xpn_cff_aer(mmn_nbr,bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_aer"
     allocate(lgn_xpn_cff_bga(mmn_nbr,bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_bga"
     allocate(lgn_xpn_cff_ice(mmn_nbr,bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_ice"
     allocate(lgn_xpn_cff_lqd(mmn_nbr,bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_lqd"
     allocate(lgn_xpn_cff_mpr(mmn_nbr,bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_mpr"
  end if ! !flg_mie
  if (flg_mie_snw) then
     ! Array dimensions: mmn,bnd,lev_bnd_snw
     allocate(lgn_xpn_cff_snw(mmn_nbr,bnd_nbr,lev_bnd_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_snw"
  end if ! !flg_mie_snw

  if (flg_mie) then
     ! Array dimensions: mmn,bnd,lev
     allocate(lgn_xpn_cff_Mie_ttl(mmn_nbr,bnd_nbr,lev_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for lgn_xpn_cff_Mie_ttl"
  end if ! !flg_mie

  ! Array dimensions: bnd,lev
  allocate(asm_prm_HG_ttl(bnd_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for asm_prm_HG_ttl"
  allocate(flx_spc_abs(bnd_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_abs"
  allocate(ntn_spc_mean(bnd_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_mean"
  allocate(odxl_spc_ttl(bnd_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for odxl_spc_ttl"
  allocate(ss_alb_fct(bnd_nbr,lev_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ss_alb_fct"

  ! Array dimensions: bnd,levp
  allocate(flx_spc_dwn(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn"
  allocate(flx_spc_dwn_dff(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_dff"
  allocate(flx_spc_dwn_drc(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_dwn_drc"
  allocate(flx_spc_upw(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for flx_spc_upw"
  allocate(lmn_spc_aa_ndr(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_spc_aa_ndr"
  allocate(ntn_spc_aa_ndr(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_ndr"
  allocate(ntn_spc_aa_zen(bnd_nbr,levp_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_zen"

  ! Array dimensions: plr,bnd
  allocate(lmn_spc_aa_sfc(plr_nbr,bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for lmn_spc_aa_sfc"
  allocate(ntn_spc_aa_sfc(plr_nbr,bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_aa_sfc"

  ! Array dimensions: azi,plr,bnd
  allocate(ntn_spc_TOA(azi_nbr,plr_nbr,bnd_nbr),stat=rcd)
  if(rcd /= 0) stop "allocate() failed for ntn_spc_TOA"

  ! Initialize level-independent arrays that depend on bnd_nbr
  do bnd_idx=1,bnd_nbr_H2O
     wvl_max(bnd_idx)=wvl_max_H2O(bnd_idx)
     wvl_min(bnd_idx)=wvl_min_H2O(bnd_idx)
  enddo                     ! end loop over bnd
  wvl_max(bnd_nbr_H2O+1)=wvl_min_H2O(bnd_nbr_H2O)
  wvl_min(bnd_nbr_H2O+1)=wvl_min_O3(bnd_nbr_pure_O3)
  do bnd_idx=bnd_nbr_H2O+2,bnd_nbr
     bnd_idx_O3=bnd_nbr_pure_O3-(bnd_idx-bnd_nbr_H2O-1)
     wvl_max(bnd_idx)=wvl_max_O3(bnd_idx_O3)
     wvl_min(bnd_idx)=wvl_min_O3(bnd_idx_O3)
  enddo                     ! end loop over bnd
  do bnd_idx=1,bnd_nbr
     wvl_ctr(bnd_idx)= &
          0.5*(wvl_max(bnd_idx)+wvl_min(bnd_idx))
     bnd(bnd_idx)=wvl_ctr(bnd_idx)
     wvl(bnd_idx)=wvl_ctr(bnd_idx)
     wvl_grd(bnd_idx)=wvl_max(bnd_idx)
     wvl_dlt(bnd_idx)= &
          wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvn_max(bnd_idx)=1.0/(100.0*wvl_min(bnd_idx))
     wvn_min(bnd_idx)=1.0/(100.0*wvl_max(bnd_idx))
     wvn_ctr(bnd_idx)= &
          0.5*(wvn_min(bnd_idx)+wvn_max(bnd_idx))
     wvn(bnd_idx)=wvn_ctr(bnd_idx)
     wvn_dlt(bnd_idx)= &
          wvn_max(bnd_idx)-wvn_min(bnd_idx)
  enddo                     ! end loop over bnd
  wvl_grd(bnd_nbr+1)=wvl_min(bnd_nbr)
  do bnd_idx=1,bnd_nbr
     if ((wvl_min(bnd_idx)<=wvl_obs_aer).and. &
          (wvl_max(bnd_idx)>wvl_obs_aer)) bnd_obs_aer=bnd_idx
     if ((wvl_min(bnd_idx)<=wvl_obs_mpr).and. &
          (wvl_max(bnd_idx)>wvl_obs_mpr)) bnd_obs_mpr=bnd_idx
     if ((wvl_min(bnd_idx)<=wvl_obs_snw).and. &
          (wvl_max(bnd_idx)>wvl_obs_snw)) bnd_obs_snw=bnd_idx
     if ((wvl_min(bnd_idx)<=wvl_obs_bga).and. &
          (wvl_max(bnd_idx)>wvl_obs_bga)) bnd_obs_bga=bnd_idx
  enddo                     ! end loop over bnd
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh ! Use xtr_fll_ngh on solar spectra
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Let user know where wavelength chips have fallen
  write (str_sng,'(a45,i2,a9)') 'Discrete ordinate computation performed with ',str_nbr,' streams'//nlc
  if (sv_cmp_plr_ngl) then
     write (plr_sng,'(a24,i2,a28)') 'Intensities reported at ',plr_nbr,' computational polar angles'//nlc
  else
     write (plr_sng,'(a24,i2,a33)') 'Intensities reported at ',plr_nbr,' evenly spaced user polar angles'//nlc
  endif
  write (azi_sng,'(a24,i2,a53,i4,a3,f7.5,a4)')  &
       'Intensities reported at ',azi_nbr, &
       ' evenly spaced azimuthal angles for (1-based) band = ', &
       bnd_dbg,' = ',wvl_ctr(bnd_dbg)*1.0e6,' um'//nlc
  write (aer_sng,'(a56,i4,a3,f7.5,a4)')  &
       'Layer aerosol optical depths saved for (1-based) band = ', &
       bnd_obs_aer,' = ',wvl_ctr(bnd_obs_aer)*1.0e6,' um'//nlc
  write (mpr_sng,'(a62,i4,a3,f7.5,a4)')  &
       'Layer snow impurity optical depths saved for (1-based) band = ', &
       bnd_obs_mpr,' = ',wvl_ctr(bnd_obs_mpr)*1.0e6,' um'//nlc
  write (snw_sng,'(a53,i4,a3,f7.5,a4)')  &
       'Layer snow optical depths saved for (1-based) band = ', &
       bnd_obs_snw,' = ',wvl_ctr(bnd_obs_snw)*1.0e6,' um'//nlc
  write (bga_sng,'(a67,i4,a3,f7.5,a4)')  &
       'Layer background aerosol optical depths saved for (1-based) band = ', &
       bnd_obs_bga,' = ',wvl_ctr(bnd_obs_bga)*1.0e6,' um'//nlc
  if (sv_cmp_tau) then
     write (opt_dep_sng,'(a33,i3,a27)') 'Radiative quantities reported at ',tau_nbr, &
          ' interface pressure levels'//nlc
  else
     write (opt_dep_sng,'(a33,i3,a29)') 'Radiative quantities reported at ',tau_nbr, &
          ' user-defined optical depths'//nlc
  endif
  ! End block of formatting diagnostic strings
  
  if (dbg_lvl>dbg_off) then
     write (6,'(35(a,/))')        &
          str_sng(1:ftn_strlen(str_sng)), &
          prf_sng(1:ftn_strlen(prf_sng)), &
          plr_sng(1:ftn_strlen(plr_sng)), &
          opt_dep_sng(1:ftn_strlen(opt_dep_sng)), &
          azi_sng(1:ftn_strlen(azi_sng)), &
          aer_sng(1:ftn_strlen(aer_sng)), &
          snw_sng(1:ftn_strlen(snw_sng)), &
          bga_sng(1:ftn_strlen(bga_sng)), &
          stt_CO2(1:ftn_strlen(stt_CO2)), &
          stt_H2O(1:ftn_strlen(stt_H2O)), &
          stt_H2OH2O(1:ftn_strlen(stt_H2OH2O)), &
          stt_Herzberg(1:ftn_strlen(stt_Herzberg)), &
          stt_ice(1:ftn_strlen(stt_ice)), &
          stt_lqd(1:ftn_strlen(stt_lqd)), &
          stt_cld_sat(1:ftn_strlen(stt_cld_sat)), &
          stt_sat_cld(1:ftn_strlen(stt_sat_cld)), &
          stt_sct_lqd(1:ftn_strlen(stt_sct_lqd)), &
          stt_vpr_H2O_abs_cld(1:ftn_strlen(stt_vpr_H2O_abs_cld)), &
          stt_OH(1:ftn_strlen(stt_OH)), &
          stt_CH4(1:ftn_strlen(stt_CH4)), &
          stt_O2(1:ftn_strlen(stt_O2)), &
          stt_O3(1:ftn_strlen(stt_O3)), &
          stt_O2O2(1:ftn_strlen(stt_O2O2)), &
          stt_O2N2(1:ftn_strlen(stt_O2N2)), &
          stt_NO2(1:ftn_strlen(stt_NO2)), &
          stt_Rayleigh(1:ftn_strlen(stt_Rayleigh)), &
          stt_aer(1:ftn_strlen(stt_aer)), &
          stt_bga(1:ftn_strlen(stt_bga)), &
          stt_rfl(1:ftn_strlen(stt_rfl)), &
          stt_mpr(1:ftn_strlen(stt_mpr)), &
          stt_snw(1:ftn_strlen(stt_snw)), &
          stt_msm(1:ftn_strlen(stt_msm)), &
          stt_flt_nst(1:ftn_strlen(stt_flt_nst)), &
          stt_slr(1:ftn_strlen(stt_slr)), &
          stt_Planck(1:ftn_strlen(stt_Planck)), &
          stt_top_lvl(1:ftn_strlen(stt_top_lvl))
  endif                     ! endif dbg
  
  if (dbg_lvl>dbg_off) then
     write (6,'(20(a,i4,/))')        &
          '# input atmosphere levels lev_atm_nbr = ',lev_atm_nbr, &
          '# input snow levels lev_snw_nbr = ',lev_snw_nbr, &
          '# total levels lev_nbr = lev_atm_nbr+lev_snw_nbr = ',lev_nbr, &
          '# input O3 bands bnd_nbr_O3 = ',bnd_nbr_O3, &
          '# input H2OH2O bands bnd_nbr_H2OH2O = ',bnd_nbr_H2OH2O, &
          '# input H2O narrow bands bnd_nbr_H2O = ',bnd_nbr_H2O, &
          '# input CO2 narrow bands bnd_nbr_CO2 = ',bnd_nbr_CO2, &
          '# input OH narrow bands bnd_nbr_OH = ',bnd_nbr_OH, &
          '# input CH4 narrow bands bnd_nbr_CH4 = ',bnd_nbr_CH4, &
          '# input O2 narrow bands bnd_nbr_O2 = ',bnd_nbr_O2, &
          '# input O2-O2 bands bnd_nbr_O2O2 = ',bnd_nbr_O2O2, &
          '# input NO2 bands bnd_nbr_NO2 = ',bnd_nbr_NO2, &
          '# input ice bands bnd_nbr_ice = ',bnd_nbr_ice, &
          '# input liq bands bnd_nbr_lqd = ',bnd_nbr_lqd, &
          '# input aer bands bnd_nbr_aer = ',bnd_nbr_aer, &
          '# input snw bands bnd_nbr_snw = ',bnd_nbr_snw, &
          '# input bga bands bnd_nbr_bga = ',bnd_nbr_bga, &
          '# input O3 bands outside H2O bands (no H2O overlap) bnd_nbr_pure_O3 = ',bnd_nbr_pure_O3, &
          'idx last pure H2O band (no O2-O3 overlap) bnd_nbr_non_O3 = ',bnd_nbr_non_O3, &
          '# total bands bnd_nbr = bnd_nbr_pure_O3+bnd_nbr_H2O = ',bnd_nbr
  endif                     ! end if dbg
  
  slr_cst_xnt_fac=slr_cst*xnt_fac
  do bnd_idx=1,bnd_nbr
     odxc_spc_OH(bnd_idx)=0.0
     odxc_spc_CH4(bnd_idx)=0.0
     odxc_spc_O2(bnd_idx)=0.0
     odxc_spc_CO2(bnd_idx)=0.0
     odxc_spc_H2OH2O(bnd_idx)=0.0
     odxc_spc_H2O(bnd_idx)=0.0
     odxc_spc_O3(bnd_idx)=0.0
     odxc_spc_O2O2(bnd_idx)=0.0
     odxc_spc_O2N2(bnd_idx)=0.0
     odxc_spc_NO2(bnd_idx)=0.0
     odxc_spc_ice(bnd_idx)=0.0
     odxc_spc_lqd(bnd_idx)=0.0
     odac_spc_aer(bnd_idx)=0.0
     odac_spc_snw(bnd_idx)=0.0
     odac_spc_bga(bnd_idx)=0.0
     odac_spc_ice(bnd_idx)=0.0
     odac_spc_lqd(bnd_idx)=0.0
     odxc_spc_aer(bnd_idx)=0.0
     odxc_spc_snw(bnd_idx)=0.0
     odxc_spc_bga(bnd_idx)=0.0
     odxc_spc_Ray(bnd_idx)=0.0
     odxc_spc_ttl(bnd_idx)=0.0
  enddo                     ! end loop over bnd
  
  ! Initialize band-independent arrays that depend on level
  do lev_idx=1,lev_nbr
     tpt_dlt_Mlk(lev_idx)=tpt(lev_idx)-tpt_Malkmus_fit
     tpt_dlt_Mlk_sqr(lev_idx)=tpt_dlt_Mlk(lev_idx)*tpt_dlt_Mlk(lev_idx)
  enddo                     ! end loop over lev
  do lev_snw_idx=1,lev_snw_nbr
     mpl_snw(lev_snw_idx)=dns_snw(lev_snw_idx)*dpt_dlt_snw(lev_snw_idx) ! [kg m-2] Mass path of snow in layer
     mpl_mpr(lev_snw_idx)=mmr_mpr_snw(lev_snw_idx)*mpl_snw(lev_snw_idx) ! [kg m-2] Mass path of impurity in layer
  enddo                  ! end loop over lev_snw
  
  ! Zero broad-band arrays that will accumulate spectral fluxes
  do lev_idx=1,lev_nbr
     j_NO2(lev_idx)=0.0
     ntn_bb_mean(lev_idx)=0.0
  enddo                     ! end loop over lev
  do lev_idx=1,levp_nbr
     do plr_idx=1,plr_nbr
        lmn_bb_aa(plr_idx,lev_idx)=0.0
        ntn_bb_aa(plr_idx,lev_idx)=0.0
     enddo                  ! end loop over plr
     flx_bb_dwn_drc(lev_idx)=0.0
     flx_bb_dwn_dff(lev_idx)=0.0
     flx_bb_upw(lev_idx)=0.0
  enddo                     ! end loop over lev
  
  ! Initialize counters either incremented or decremented in main loop over bands
  
  ! Rebin continuum cross sections and quantum yields onto RT wavelength grid
  ! Quantum yields should be one at high energy and zero at low energy
  xtr_typ_LHS=xtr_prt_ngh+xtr_fll_ngh
  xtr_typ_RHS=xtr_prt_ngh+xtr_fll_ngh
  call rbn_vec(bnd_nbr_NO2,wvl_grd_NO2,qnt_yld_NO2_dsk, &
       bnd_nbr,wvl_grd,qnt_yld_NO2, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Re-bin spectral surface reflectance using nearest neighbor extrapolation
  if (flg_rfl) then
     call rbn_vec(bnd_nbr_rfl,wvl_grd_rfl,rfl_spc_sfc_dsk, &
          bnd_nbr,wvl_grd,rfl_spc_sfc, &
          xtr_typ_LHS,xtr_typ_RHS)
     if (allocated(rfl_spc_sfc_dsk)) deallocate(rfl_spc_sfc_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for rfl_spc_sfc_dsk"
     if (allocated(wvl_grd_rfl)) deallocate(wvl_grd_rfl,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_rfl"
  endif ! !flg_rfl
  ! Re-bin optical properties using nearest neighbor extrapolation
  if (flg_aer) then
     call rbn_vec(bnd_nbr_aer,wvl_grd_aer,abs_cff_mss_aer_dsk, &
          bnd_nbr,wvl_grd,abs_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_aer,wvl_grd_aer,ext_cff_mss_aer_dsk, &
          bnd_nbr,wvl_grd,ext_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_aer,wvl_grd_aer,sca_cff_mss_aer_dsk, &
          bnd_nbr,wvl_grd,sca_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_aer,wvl_grd_aer,asm_prm_aer_dsk, &
          bnd_nbr,wvl_grd,asm_prm_aer,xtr_typ_LHS,xtr_typ_RHS)
     if (flg_mie_aer) then
        do mmn_idx=1,mmn_nbr
           call rbn_vec(bnd_nbr_aer,wvl_grd_aer,lgn_xpn_cff_aer_dsk(mmn_idx,:), &
                bnd_nbr,wvl_grd,lgn_xpn_cff_aer(mmn_idx,:),xtr_typ_LHS,xtr_typ_RHS)
        end do ! end loop over mmn
     end if ! !flg_mie_aer
     if (allocated(abs_cff_mss_aer_dsk)) deallocate(abs_cff_mss_aer_dsk,stat=rcd) !
     if (allocated(ext_cff_mss_aer_dsk)) deallocate(ext_cff_mss_aer_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for ext_cff_mss_aer_dsk"
     if (allocated(sca_cff_mss_aer_dsk)) deallocate(sca_cff_mss_aer_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_aer_dsk"
     if (allocated(asm_prm_aer_dsk)) deallocate(asm_prm_aer_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_aer_dsk"
     if (allocated(wvl_grd_aer)) deallocate(wvl_grd_aer,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_aer"
     if (allocated(lgn_xpn_cff_aer_dsk)) deallocate(lgn_xpn_cff_aer_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_aer_dsk"
  endif ! !flg_aer
  if (flg_bga) then
     call rbn_vec(bnd_nbr_bga,wvl_grd_bga,abs_cff_mss_bga_dsk, &
          bnd_nbr,wvl_grd,abs_cff_mss_bga,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_bga,wvl_grd_bga,sca_cff_mss_bga_dsk, &
          bnd_nbr,wvl_grd,sca_cff_mss_bga,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_bga,wvl_grd_bga,asm_prm_bga_dsk, &
          bnd_nbr,wvl_grd,asm_prm_bga,xtr_typ_LHS,xtr_typ_RHS)
     if (flg_mie_bga) then
        do mmn_idx=1,mmn_nbr
           call rbn_vec(bnd_nbr_bga,wvl_grd_bga,lgn_xpn_cff_bga_dsk(mmn_idx,:), &
                bnd_nbr,wvl_grd,lgn_xpn_cff_bga(mmn_idx,:),xtr_typ_LHS,xtr_typ_RHS)
        end do ! end loop over mmn
     end if ! !flg_mie_bga
     if (allocated(abs_cff_mss_bga_dsk)) deallocate(abs_cff_mss_bga_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for abs_cff_mss_bga_dsk"
     if (allocated(sca_cff_mss_bga_dsk)) deallocate(sca_cff_mss_bga_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_bga_dsk"
     if (allocated(asm_prm_bga_dsk)) deallocate(asm_prm_bga_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_bga_dsk"
     if (allocated(wvl_grd_bga)) deallocate(wvl_grd_bga,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_bga"
     if (allocated(lgn_xpn_cff_bga_dsk)) deallocate(lgn_xpn_cff_bga_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_bga_dsk"
  endif ! !flg_bga
  if (flg_ice) then
     call rbn_vec(bnd_nbr_ice,wvl_grd_ice,abs_cff_mss_ice_dsk, &
          bnd_nbr,wvl_grd,abs_cff_mss_ice,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_ice,wvl_grd_ice,sca_cff_mss_ice_dsk, &
          bnd_nbr,wvl_grd,sca_cff_mss_ice,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_ice,wvl_grd_ice,asm_prm_ice_dsk, &
          bnd_nbr,wvl_grd,asm_prm_ice,xtr_typ_LHS,xtr_typ_RHS)
     if (flg_mie_ice) then
        do mmn_idx=1,mmn_nbr
           call rbn_vec(bnd_nbr_ice,wvl_grd_ice,lgn_xpn_cff_ice_dsk(mmn_idx,:), &
                bnd_nbr,wvl_grd,lgn_xpn_cff_ice(mmn_idx,:),xtr_typ_LHS,xtr_typ_RHS)
        end do ! end loop over mmn
     end if ! !flg_mie_ice
     if (allocated(abs_cff_mss_ice_dsk)) deallocate(abs_cff_mss_ice_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for abs_cff_mss_ice_dsk"
     if (allocated(sca_cff_mss_ice_dsk)) deallocate(sca_cff_mss_ice_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_ice_dsk"
     if (allocated(asm_prm_ice_dsk)) deallocate(asm_prm_ice_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_ice_dsk"
     if (allocated(wvl_grd_ice)) deallocate(wvl_grd_ice,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_ice"
     if (allocated(lgn_xpn_cff_ice_dsk)) deallocate(lgn_xpn_cff_ice_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_ice_dsk"
  endif ! !flg_ice
  if (flg_lqd) then
     call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,abs_cff_mss_lqd_dsk, &
          bnd_nbr,wvl_grd,abs_cff_mss_lqd,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,sca_cff_mss_lqd_dsk, &
          bnd_nbr,wvl_grd,sca_cff_mss_lqd,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,asm_prm_lqd_dsk, &
          bnd_nbr,wvl_grd,asm_prm_lqd,xtr_typ_LHS,xtr_typ_RHS)
     if (flg_mie_lqd) then
        do mmn_idx=1,mmn_nbr
           call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,lgn_xpn_cff_lqd_dsk(mmn_idx,:), &
                bnd_nbr,wvl_grd,lgn_xpn_cff_lqd(mmn_idx,:),xtr_typ_LHS,xtr_typ_RHS)
        end do ! end loop over mmn
     end if ! !flg_mie_lqd
     if (allocated(abs_cff_mss_lqd_dsk)) deallocate(abs_cff_mss_lqd_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for abs_cff_mss_lqd_dsk"
     if (allocated(sca_cff_mss_lqd_dsk)) deallocate(sca_cff_mss_lqd_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_lqd_dsk"
     if (allocated(asm_prm_lqd_dsk)) deallocate(asm_prm_lqd_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_lqd_dsk"
     if (allocated(wvl_grd_lqd)) deallocate(wvl_grd_lqd,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_lqd"
     if (allocated(lgn_xpn_cff_lqd_dsk)) deallocate(lgn_xpn_cff_lqd_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_lqd_dsk"
  endif ! !flg_lqd
  if (flg_mpr) then
     call rbn_vec(bnd_nbr_mpr,wvl_grd_mpr,abs_cff_mss_mpr_dsk, &
          bnd_nbr,wvl_grd,abs_cff_mss_mpr,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_mpr,wvl_grd_mpr,sca_cff_mss_mpr_dsk, &
          bnd_nbr,wvl_grd,sca_cff_mss_mpr,xtr_typ_LHS,xtr_typ_RHS)
     call rbn_vec(bnd_nbr_mpr,wvl_grd_mpr,asm_prm_mpr_dsk, &
          bnd_nbr,wvl_grd,asm_prm_mpr,xtr_typ_LHS,xtr_typ_RHS)
     if (flg_mie_mpr) then
        do mmn_idx=1,mmn_nbr
           call rbn_vec(bnd_nbr_mpr,wvl_grd_mpr,lgn_xpn_cff_mpr_dsk(mmn_idx,:), &
                bnd_nbr,wvl_grd,lgn_xpn_cff_mpr(mmn_idx,:),xtr_typ_LHS,xtr_typ_RHS)
        end do ! end loop over mmn
     end if ! !flg_mie_mpr
     if (allocated(abs_cff_mss_mpr_dsk)) deallocate(abs_cff_mss_mpr_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for abs_cff_mss_mpr_dsk"
     if (allocated(sca_cff_mss_mpr_dsk)) deallocate(sca_cff_mss_mpr_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_mpr_dsk"
     if (allocated(asm_prm_mpr_dsk)) deallocate(asm_prm_mpr_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_mpr_dsk"
     if (allocated(wvl_grd_mpr)) deallocate(wvl_grd_mpr,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_mpr"
     if (allocated(lgn_xpn_cff_mpr_dsk)) deallocate(lgn_xpn_cff_mpr_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_mpr_dsk"
  endif ! !flg_mpr
  if (flg_snw) then
     do lev_bnd_snw_idx=1,lev_bnd_snw_nbr
        call rbn_vec(bnd_nbr_snw,wvl_grd_snw,abs_cff_mss_snw_dsk(:,lev_bnd_snw_idx), &
             bnd_nbr,wvl_grd,abs_cff_mss_snw(:,lev_bnd_snw_idx),xtr_typ_LHS,xtr_typ_RHS)
        call rbn_vec(bnd_nbr_snw,wvl_grd_snw,sca_cff_mss_snw_dsk(:,lev_bnd_snw_idx), &
             bnd_nbr,wvl_grd,sca_cff_mss_snw(:,lev_bnd_snw_idx),xtr_typ_LHS,xtr_typ_RHS)
        call rbn_vec(bnd_nbr_snw,wvl_grd_snw,asm_prm_snw_dsk(:,lev_bnd_snw_idx), &
             bnd_nbr,wvl_grd,asm_prm_snw(:,lev_bnd_snw_idx),xtr_typ_LHS,xtr_typ_RHS)

        if (flg_mie_snw) then
           do mmn_idx=1,mmn_nbr
              call rbn_vec(bnd_nbr_snw,wvl_grd_snw,lgn_xpn_cff_snw_dsk(mmn_idx,:,lev_bnd_snw_idx), &
                   bnd_nbr,wvl_grd,lgn_xpn_cff_snw(mmn_idx,:,lev_bnd_snw_idx),xtr_typ_LHS,xtr_typ_RHS)
           end do ! end loop over mmn_idx
        end if ! !flg_mie_snw
     end do ! end loop over lev_bnd_snw_idx
     if (flg_mie) write (6,'(a,a,i2,a)') prg_nm(1:ftn_strlen(prg_nm)),': Utilizing ',mmn_nbr,' phase function moments'
     if (allocated(abs_cff_mss_snw_dsk)) deallocate(abs_cff_mss_snw_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for abs_cff_mss_snw_dsk"
     if (allocated(sca_cff_mss_snw_dsk)) deallocate(sca_cff_mss_snw_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for sca_cff_mss_snw_dsk"
     if (allocated(asm_prm_snw_dsk)) deallocate(asm_prm_snw_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for asm_prm_snw_dsk"
     if (allocated(lgn_xpn_cff_snw_dsk)) deallocate(lgn_xpn_cff_snw_dsk,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for lgn_xpn_cff_snw_dsk"
     if (allocated(wvl_grd_snw)) deallocate(wvl_grd_snw,stat=rcd) !
     if(rcd /= 0) stop "deallocate() failed for wvl_grd_snw"
  endif ! !flg_snw
  if (flg_msm) then
     ! Space for diagnostic CCY83 (delta-Eddington/adding) computation
     ! Array dimensions: bnd
     allocate(alb_dff_spc_snw_dea(bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for alb_dff_spc_snw_dea"
     allocate(alb_drc_spc_snw_dea(bnd_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for alb_drc_spc_snw_dea"
     ! Array dimensions: bnd,levp_snw
     allocate(rfl_dff_spc_snw_cnt(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_dff_spc_snw_cnt"
     allocate(trn_dff_spc_snw_cnt(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_dff_spc_snw_cnt"
     allocate(rfl_drc_upw_snp(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_drc_upw_snp"
     allocate(rfl_dff_upw_snp(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_dff_upw_snp"
     allocate(rfl_dff_dwn_snp(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_dff_dwn_snp"
     allocate(trn_ttl_drc_snp(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_ttl_drc_snp"
     allocate(trn_drc_drc_snp(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_drc_drc_snp"
     ! Array dimensions: bnd,lev_snw
     allocate(rfl_ddm_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_ddm_spc_snw"
     allocate(trn_ddm_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_ddm_spc_snw"
     allocate(rfl_dff_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_dff_spc_snw"
     allocate(trn_dff_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_dff_spc_snw"
     allocate(rfl_drc_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for rfl_drc_spc_snw"
     allocate(trn_ttl_drc_spc_snw(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_ttl_drc_spc_snw"
     allocate(trn_drc_drc(bnd_nbr,levp_snw_nbr),stat=rcd)
     if(rcd /= 0) stop "allocate() failed for trn_drc_drc"
  endif ! !flg_msm
  ! Absorption cross-sections should be 0 outside data range
  xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
  xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
  call rbn_vec(bnd_nbr_NO2,wvl_grd_NO2,abs_xsx_NO2_dsk, &
       bnd_nbr,wvl_grd,abs_xsx_NO2, &
       xtr_typ_LHS,xtr_typ_RHS)
  call rbn_vec(bnd_nbr_O2O2,wvl_grd_O2O2,abs_xsx_O2O2_dsk, &
       bnd_nbr,wvl_grd,abs_xsx_O2O2, &
       xtr_typ_LHS,xtr_typ_RHS)
  call rbn_vec(bnd_nbr_H2OH2O,wvl_grd_H2OH2O,abs_xsx_H2OH2O_dsk, &
       bnd_nbr,wvl_grd,abs_xsx_H2OH2O, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Luminosity response should be weighted 0 outside data range
  xtr_typ_LHS=xtr_prt_wgt+xtr_fll_nil
  xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil
  call rbn_vec(bnd_nbr_lmn,wvl_grd_lmn,lmn_SRF_dsk, &
       bnd_nbr,wvl_grd,lmn_SRF, &
       xtr_typ_LHS,xtr_typ_RHS)
  if (allocated(lmn_SRF_dsk)) deallocate(lmn_SRF_dsk,stat=rcd) !
  if(rcd /= 0) stop "deallocate() failed for lmn_SRF_dsk"
  if (allocated(wvl_grd_lmn)) deallocate(wvl_grd_lmn,stat=rcd) !
  if(rcd /= 0) stop "deallocate() failed for wvl_grd_lmn"
  do bnd_idx=1,bnd_nbr
     bnd_wgt_lmn(bnd_idx)=lmn_SRF(bnd_idx)
  enddo                  ! end loop over bnd
  
  chn_SRF_msk(:)=0 ! CEWI lf95
  if (.not.mode_std) then
     do bnd_idx=1,bnd_nbr
        do chn_idx=1,chn_nbr
           if (chn_SRF(bnd_idx,chn_idx) > 0.0) then
              chn_SRF_msk(bnd_idx)=1
           endif
        enddo
     enddo
  endif ! endif mode_std
  
  ! If necessary, make user-specified extinction optical depth
  ! consistent with particle optical properties
  if (cmd_ln_odxc_obs_aer) then
     do lev_idx=1,lev_nbr
        odxl_tmp=mpl_aer(lev_idx)*ext_cff_mss_aer(bnd_obs_aer)
        mpl_tmp=mpl_aer(lev_idx)
        if (odxl_tmp /= odxl_obs_aer(lev_idx)) then
           if (ext_cff_mss_aer(bnd_obs_aer) == 0.0) then
              write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR ext_cff_mss_aer(bnd_obs_aer) = 0.0 in denominator'
              call abort
           endif            ! endif
           mpl_aer(lev_idx)=odxl_obs_aer(lev_idx)/ext_cff_mss_aer(bnd_obs_aer)
           write (6,'(a,a,i4,a,e9.2,a,e9.2,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': WARNING changing mpl_aer(',lev_idx,') from = ',mpl_tmp,' to ',mpl_aer(lev_idx),' kg m-2'
        endif               ! endif adjusting mass path
     enddo                  ! end loop over lev
  endif                     ! end if overriding CLM profile odxc_obs_aer

  ! End section 1: Initialization
  ! Begin section 2: Main computation loop over all bands
  ! begin main loop over bands
  ! assignment of counting indices into input data arrays
  ! a level loop for zeroing band quantities that depend on level
  ! computation and binning of process-specific optical properties
  ! summation of individual optical properties into bulk layer props
  ! assignment of bulk properties into actual DISORT() input arrays
  ! call DISORT()
  ! assignment of DISORT() output to final output storage arrays
  ! end main loop over bands
  
  ! OpenMP threading:
  ! Distinguish shared from private values
  ! Private variables which depend on previous work should be firstprivate
  ! LHS variables with a bnd_idx dependence are shared
  ! RHS variables (not written to in parallel region) are shared to save memory
  ! First seven lines of shared() and private() lists follow disort2() calling list
  ! Normal continuation line restrictions apply, no commas between clauses, no parameters in clauses
  !$omp parallel default(none) 
  !$omp$private(dtauc,ssalb,pmom,wvnmlo)
  !$omp$private(wvnmhi,utau)
  !$omp$private(fbeam)
  !$omp$private(albedo)
  !$omp$private(plank,header)
  !$omp$private(rfldir,rfldn)
  !$omp$private(flup,dfdt,uavg,uu,albmed,trnmed,u0u)
  !$omp$private(bnd_idx,lev_idx,bnd_idx_O3,bnd_idx_CH4,bnd_idx_H2O,bnd_idx_CO2,bnd_idx_tmp_O3,i21)
  !$omp$private(mmn_idx,tau_idx,plr_idx,azi_idx)
  !$omp$private(phi_wgt,psi_wgt,u_bar,prs_bar,float_foo,tau)
  !$omp$private(opt_dep_ITOD_H2O,opt_dep_ITOD_O2,opt_dep_ITOD_OH,opt_dep_ITOD_CH4)
  !$omp$private(opt_dep_ITOD_CO2_hires,opt_dep_LTOD_CO2_hires,trn_LT_CO2_hires,trn_ALT_CO2)
  !$omp$private(odsl_Ray,odsl_ice,odsl_lqd,odsl_aer,odsl_bga,odal_ice,odal_lqd,odal_aer,odal_bga,odal_OH,odal_CH4)
  !$omp$private(odal_O2,odal_O3,odal_O2O2,odal_O2N2,odal_NO2,odal_CO2,odal_H2OH2O,odal_H2O)
  !$omp$private(sca_cff_mss_Ray,sca_frc_HG,sca_frc_Mie)
  !$omp$shared(nlyr,nmom,temper)
  !$omp$shared(usrtau,ntau,nstr,usrang,numu)
  !$omp$shared(umu,nphi,phi,ibcnd,umu0,phi0)
  !$omp$shared(fisot,lamber,btemp,ttemp,temis)
  !$omp$shared(onlyfl,accur,prnt)
  !$omp$shared(dbg_lvl,thr_nbr,brdf_typ)
  !$omp$shared(idx_rfr_air_STP,prg_nm,sv_ntn,sv_cmp_tau,single_bnd_computation)
  !$omp$shared(hb,br,f_iso,f_vol,f_geo,nrm_cff_CM,idx_rfr_sfc,b0,hh,w,nrm_rfl_M,k_cff_M,nrm_cff_RP,k_cff_RP,g_phs,nrm_rfl_LS)
  !$omp$shared(wvl_ctr,wvl_dlt,wvl_min,wvl_max,wvn_min,wvn_max,bnd_dbg,tst_case_Rayleigh,tst_case_HG)
  !$omp$shared(flg_Rayleigh,flg_ice,flg_lqd,flg_aer,flg_bga,flg_H2O,flg_H2OH2O,flg_OH,flg_CH4,flg_O2,flg_CO2)
  !$omp$shared(flg_O3,flg_O2O2,flg_O2N2,flg_NO2,bnd_obs_aer,bnd_obs_bga,flg_Planck,wvl_Planck,mode_std)
  !$omp$shared(prs,prs_ntf,tpt,mmw_mst_air,mpl_mst_air,grv,pi,ocn_msk)
  !$omp$shared(odal_obs_aer,odsl_obs_aer,odxl_obs_aer,odal_obs_bga,odsl_obs_bga,odxl_obs_bga)
  !$omp$shared(slr_zen_ngl_cos,alb_sfc_vsb_drc,alb_sfc_vsb_dff,alb_sfc_NIR_drc,alb_sfc_NIR_dff)
  !$omp$shared(slr_cst_xnt_fac,flx_slr_frc,chn_SRF_msk)
  !$omp$shared(lev_nbr,levp_nbr,plr_nbr,mmn_nbr,azi_nbr,tau_nbr,bnd_nbr_H2O,bnd_nbr_pure_O3,bnd_nbr_non_O3,bnd_nbr_O3)
  !$omp$shared(abs_xsx_O3,npl_O3,wvl_min_O3,wvl_max_O3,abs_xsx_O2,npl_O2,bnd_nbr,abs_xsx_H2OH2O,npl_H2OH2O)
  !$omp$shared(abs_xsx_O2O2,npl_O2O2,abs_xsx_NO2,npl_NO2)
  !$omp$shared(abs_cff_mss_lqd,mpl_LWP,sca_cff_mss_lqd,bnd_nbr_lqd,wvl_min_lqd,wvl_max_lqd,asm_prm_lqd)
  !$omp$shared(abs_cff_mss_aer,mpl_aer,sca_cff_mss_aer,bnd_nbr_aer,wvl_min_aer,wvl_max_aer,asm_prm_aer)
  !$omp$shared(abs_cff_mss_bga,mpl_bga,sca_cff_mss_bga,bnd_nbr_bga,wvl_min_bga,wvl_max_bga,asm_prm_bga)
  !$omp$shared(abs_cff_mss_ice,mpl_IWP,sca_cff_mss_ice,bnd_nbr_ice,wvl_min_ice,wvl_max_ice,asm_prm_ice)
  !$omp$shared(lgn_xpn_cff_Mie_ttl,asm_prm_HG_ttl,asm_prm_lqd,asm_prm_ice,asm_prm_aer,asm_prm_bga)
  !$omp$shared(flx_spc_dwn,flx_spc_dwn_dff,flx_spc_dwn_drc,flx_spc_net,flx_spc_upw)
  !$omp$shared(ntn_bb_aa,ntn_spc_aa_ndr,ntn_spc_aa_ndr_TOA,ntn_spc_aa_ndr_sfc,ntn_spc_TOA)
  !$omp$shared(ntn_spc_aa_sfc,ntn_spc_aa_zen,ntn_spc_aa_zen_sfc,ntn_spc_chn,ntn_spc_mean)
  !$omp$shared(odac_spc_aer,odac_spc_snw,odac_spc_mpr,odac_spc_bga,odac_spc_ice,odac_spc_lqd,odal_spc_ttl)
  !$omp$shared(odsl_spc_ttl,odxc_spc_CO2,odxc_spc_H2O,odxc_spc_H2OH2O,odxc_spc_NO2)
  !$omp$shared(odxc_spc_O2,odxc_spc_O2N2,odxc_spc_O2O2,odxc_spc_O3,odxc_spc_OH,odxc_spc_CH4)
  !$omp$shared(odxc_spc_Ray,odxc_spc_aer,odxc_spc_aer,odxc_spc_bga,odxc_spc_ice,odxc_spc_lqd)
  !$omp$shared(odxl_spc_ttl,odxc_spc_ttl)
  !$omp$shared(S_d_abs_cff_mss_H2O,S_p_abs_cff_mss_H2O,mpl_H2O,A_phi_H2O,B_phi_H2O,A_psi_H2O,B_psi_H2O,q_H2O)
  !$omp$shared(S_d_abs_cff_mss_OH,S_p_abs_cff_mss_OH,mpl_OH,A_phi_OH,B_phi_OH,A_psi_OH,B_psi_OH,q_OH)
  !$omp$shared(S_d_abs_cff_mss_CH4,S_p_abs_cff_mss_CH4,mpl_CH4,A_phi_CH4,B_phi_CH4,A_psi_CH4,B_psi_CH4,q_CH4)
  !$omp$shared(S_d_abs_cff_mss_O2,S_p_abs_cff_mss_O2,mpl_O2,A_phi_O2,B_phi_O2,A_psi_O2,B_psi_O2,q_O2)
  !$omp$shared(S_d_abs_cff_mss_CO2,S_p_abs_cff_mss_CO2,mpl_CO2,A_phi_CO2,B_phi_CO2,A_psi_CO2,B_psi_CO2,q_CO2)
  !$omp$shared(ss_alb_fct,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr)
#ifdef _OPENMP
  !$omp single
  thr_nbr=omp_get_num_threads() ! [nbr] OpenMP number of threads
  write (6,'(a,a36,i3,a8)') prg_nm(1:ftn_strlen(prg_nm)),': INFO OpenMP multi-threading using ',thr_nbr,' threads'
  !$omp end single
#else /* not _OPENMP */
  write (6,'(a,a45)') prg_nm(1:ftn_strlen(prg_nm)),': INFO Not attempting OpenMP multi-threading'
#endif /* not _OPENMP */
  !$omp do
  do bnd_idx=1,bnd_nbr
     
     ! Skip calculation if band is outside of all instrument channels
     if (.not.mode_std.and.chn_SRF_msk(bnd_idx) /= 1) then
        goto 999
     endif
     
     !bnd_idx_brdf=bnd_idx
     ! Surface BRDF parametrization
     ! EOS
     if (brdf_typ == 0) then
        hb=2.0
        br=1.0
        f_iso=0.1
        f_vol=0.1
        f_geo=0.1
     endif
     ! Cox and Munk
     if (brdf_typ == 1) then
        ! nrm_cff_CM=nrm_cff_CM_spc(bnd_idx_brdf)
        ! idx_rfr_sfc=idx_rfr_sfc_spc(bnd_idx_brdf)
        ! wnd_spd=4.0 ! wind speed is read in from fl_clm
        nrm_cff_CM=pi / 4.0
        idx_rfr_sfc=1.4
     endif
     ! Hapke
     if (brdf_typ == 2) then
        b0=1.0
        hh=0.06
        w=0.6
     endif
     ! Minnaert
     if (brdf_typ == 3) then
        nrm_rfl_M=0.01
        k_cff_M=0.9
     endif
     ! RahmanPinty
     if (brdf_typ == 4) then
        nrm_cff_RP=0.1
        k_cff_RP=0.9
        g_phs=0.9
     endif
     ! LommelSeeliger
     if (brdf_typ == 5) then
        nrm_rfl_LS=0.1
     endif
     
     if (dbg_lvl>dbg_off) write (6,'(a1)',advance="no") '.'
     if (dbg_lvl==dbg_old) write (6,'(i4)',advance="no") bnd_idx
     
     ! O3 data is reverse indexed with an offset relative to 
     ! combined band data, so index juggling is necessary.
     bnd_idx_O3=bnd_nbr_pure_O3-(bnd_idx-bnd_nbr_H2O-1)
     
     ! OH data is indexed identically to the combined band data
     !bnd_idx_OH=bnd_idx
     
     ! CH4 data is indexed identically to the combined band data
     bnd_idx_CH4=bnd_idx
     
     ! O2 data is indexed identically to the combined band data
     !bnd_idx_O2=bnd_idx
     
     ! H2O data is indexed identically to the combined band data
     bnd_idx_H2O=bnd_idx
     
     ! Before anything is actually computed for this band, 
     ! zero out all optical depth arrays. This is important
     ! both as a failsafe programming technique and for the case
     ! when arguments are given to ignore specified radiative
     ! processes, e.g., the pure Rayleigh scattering case.
     ! Also, at least one optical property, odal_O2, is
     ! accumulated between two processes (Herzberg continuum absorption
     ! and line absorption near 0.69, 0.77, and 1.28 um in the 
     ! case of O2).
     do lev_idx=1,lev_nbr
        odsl_Ray(lev_idx)=0.0
        odsl_ice(lev_idx)=0.0
        odsl_lqd(lev_idx)=0.0
        odsl_aer(lev_idx)=0.0
        odsl_mpr(lev_idx)=0.0
        odsl_snw(lev_idx)=0.0
        odsl_bga(lev_idx)=0.0
        odal_ice(lev_idx)=0.0
        odal_lqd(lev_idx)=0.0
        odal_aer(lev_idx)=0.0
        odal_mpr(lev_idx)=0.0
        odal_snw(lev_idx)=0.0
        odal_bga(lev_idx)=0.0
        odal_OH(lev_idx)=0.0
        odal_CH4(lev_idx)=0.0
        odal_O2(lev_idx)=0.0
        odal_O3(lev_idx)=0.0
        odal_O2O2(lev_idx)=0.0
        odal_O2N2(lev_idx)=0.0
        odal_NO2(lev_idx)=0.0
        odal_CO2(lev_idx)=0.0
        odal_H2OH2O(lev_idx)=0.0
        odal_H2O(lev_idx)=0.0
     enddo                  ! end loop over lev
     
     ! Malkmus narrow band parameters need to be evaluated for H2O, OH, CH4, O2 and CO2 bands
     ! H2O band data is indexed exactly as the combined band data so no index juggling is necessary
     ! O2 data is presumably indexed the same as H2O
     ! OH data is presumably indexed the same as H2O
     ! CH4 data covers the same range at twice the resolution
     ! CO2 data covers the same range at twice the resolution
     ! The same physics applies to all three gases
     if (bnd_idx<=bnd_nbr_H2O) then
        
        ! Line absorption computation #1: H2O
        do lev_idx=1,lev_nbr
           ! Phi and psi weights are mid-layer quantities
           phi_wgt(lev_idx)= &
                exp(A_phi_H2O(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_phi_H2O(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
           psi_wgt(lev_idx)= &
                exp(A_psi_H2O(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_psi_H2O(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
        enddo               ! end loop over lev
        
        ! U_bar and prs_bar are layer-interface quantities
        ! Compute top interface
        u_bar(1)=q_H2O(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
        prs_bar(1)=q_H2O(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
        
        ! Sum-up integrands
        ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g., q_H2O*dp/g = mpl_H2O ([kg m-2])
        ! Using stored values saves lots of floating point operations
        do lev_idx=2,levp_nbr
           u_bar(lev_idx)=u_bar(lev_idx-1)+ &
                mpl_H2O(lev_idx-1)*phi_wgt(lev_idx-1)
           prs_bar(lev_idx)=prs_bar(lev_idx-1)+ &
                mpl_H2O(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
        enddo               ! end loop over lev
        
        ! Normalize pressure-weighted mass path
        do lev_idx=1,levp_nbr
           prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
        enddo               ! end loop over lev
        
        ! Compute absorption optical depth
        ! NB: The transmission being computed below is an interface quantity
        ! We could also have computed the transmission of a layer as a mid-layer 
        ! quantity directly, but BPB advises against this for subtle reasons
        do lev_idx=1,levp_nbr
           float_foo= &
                sqrt( &
                1.0+4.0*S_p_abs_cff_mss_H2O(bnd_idx)*u_bar(lev_idx)/ &
                (prs_bar(lev_idx)*slr_zen_ngl_cos) &
                )
           opt_dep_ITOD_H2O(lev_idx)=(float_foo-1.0)* &
                0.5*S_d_abs_cff_mss_H2O(bnd_idx)*prs_bar(lev_idx)/ &
                S_p_abs_cff_mss_H2O(bnd_idx)
        enddo               ! end loop over lev
        
        ! Compute layer "monochromatic" absorption optical depth 
        ! directly from the difference of the column "monochromatic"
        ! optical depths rather than from the log of the difference of
        ! the column transmissions in order to avoid numerical 
        ! difficulties with saturated lines when the transmission is 0.0
        ! This log transmission method fails in single precision
        do lev_idx=1,lev_nbr
           odal_H2O(lev_idx)=slr_zen_ngl_cos* &
                (opt_dep_ITOD_H2O(lev_idx+1)-opt_dep_ITOD_H2O(lev_idx))
           if (.not.flg_vpr_H2O_abs_cld) then
              ! If H2O vapor absorption in cloud is not allowed...
              if (mpl_CWP(lev_idx) > 0.0) odal_H2O(lev_idx)=0.0
           endif ! endif flg_vpr_H2O_abs_cld
        enddo               ! end loop over lev
        
        ! Line absorption computation #2: OH
        do lev_idx=1,lev_nbr
           
           ! Phi and psi weights are mid-layer quantities.
           phi_wgt(lev_idx)= &
                exp(A_phi_OH(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_phi_OH(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
           psi_wgt(lev_idx)= &
                exp(A_psi_OH(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_psi_OH(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
        enddo               ! end loop over lev
        
        ! U_bar and prs_bar are layer-interface quantities
        ! Compute top interface
        u_bar(1)=q_OH(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
        prs_bar(1)=q_OH(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
        
        ! Sum-up integrands
        ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g.,
        ! q_OH*dp/g = mpl_OH ([kg m-2]). Using stored values
        ! saves lots of floating point operations.
        do lev_idx=2,levp_nbr
           u_bar(lev_idx)=u_bar(lev_idx-1)+ &
                mpl_OH(lev_idx-1)*phi_wgt(lev_idx-1)
           prs_bar(lev_idx)=prs_bar(lev_idx-1)+ &
                mpl_OH(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
        enddo               ! end loop over lev
        
        ! Normalize pressure-weighted mass path
        do lev_idx=1,levp_nbr
           if (u_bar(lev_idx)>0.0) then
              prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
           else
              ! NB: Allowing OH concentration to be zero requires a more bulletproof narrow band code than otherwise
              ! If u_bar == 0.0 then prs_bar := 1 to avoid overflow in optical depth later
              prs_bar(lev_idx)=1.0
           endif
        enddo               ! end loop over lev
        
        ! Compute absorption optical depth
        ! NB: The transmission being computed below is an interface
        ! quantity. We could also have computed the transmission of
        ! a layer as a mid-layer quantity directly, but BPB advises
        ! against this for subtle reasons.
        do lev_idx=1,levp_nbr
           float_foo= &
                sqrt( &
                1.0+4.0*S_p_abs_cff_mss_OH(bnd_idx)*u_bar(lev_idx)/ &
                (prs_bar(lev_idx)*slr_zen_ngl_cos) &
                )
           ! NB: OH line strengths are often underflows in single precision,
           ! so S_p_abs_cff_mss_OH must either be checked for this, or else
           ! set to a reasonable value in htrn2nb.
           opt_dep_ITOD_OH(lev_idx)=(float_foo-1.0)* &
                0.5*S_d_abs_cff_mss_OH(bnd_idx)*prs_bar(lev_idx)/ &
                S_p_abs_cff_mss_OH(bnd_idx)
        enddo               ! end loop over lev
        
        ! Compute layer "monochromatic" absorption optical depth 
        ! directly from the difference of the column "monochromatic"
        ! optical depths rather than from the log of the difference of
        ! the column transmissions in order to avoid numerical 
        ! difficulties with saturated lines when the transmission is 0.0
        ! This log transmission method fails in single precision.
        do lev_idx=1,lev_nbr
           odal_OH(lev_idx)=odal_OH(lev_idx)+ &
                slr_zen_ngl_cos* &
                (opt_dep_ITOD_OH(lev_idx+1)-opt_dep_ITOD_OH(lev_idx))
        enddo               ! end loop over lev
        
        ! Line absorption computation #4: O2
        do lev_idx=1,lev_nbr
           
           ! Phi and psi weights are mid-layer quantities.
           phi_wgt(lev_idx)= &
                exp(A_phi_O2(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_phi_O2(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
           psi_wgt(lev_idx)= &
                exp(A_psi_O2(bnd_idx)*tpt_dlt_Mlk(lev_idx)+ &
                B_psi_O2(bnd_idx)*tpt_dlt_Mlk_sqr(lev_idx))
        enddo               ! end loop over lev
        
        ! U_bar and prs_bar are layer-interface quantities
        ! Compute top interface
        u_bar(1)=q_O2(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
        prs_bar(1)=q_O2(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
        
        ! Sum-up integrands
        ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g.,
        ! q_O2*dp/g = mpl_O2 ([kg m-2]). Using stored values
        ! saves lots of floating point operations.
        do lev_idx=2,levp_nbr
           u_bar(lev_idx)=u_bar(lev_idx-1)+ &
                mpl_O2(lev_idx-1)*phi_wgt(lev_idx-1)
           prs_bar(lev_idx)=prs_bar(lev_idx-1)+ &
                mpl_O2(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
        enddo               ! end loop over lev
        
        ! Normalize pressure-weighted mass path
        do lev_idx=1,levp_nbr
           prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
        enddo               ! end loop over lev
        
        ! Compute absorption optical depth
        ! NB: The transmission being computed below is an interface
        ! quantity. We could also have computed the transmission of
        ! a layer as a mid-layer quantity directly, but BPB advises
        ! against this for subtle reasons.
        do lev_idx=1,levp_nbr
           float_foo= &
                sqrt( &
                1.0+4.0*S_p_abs_cff_mss_O2(bnd_idx)*u_bar(lev_idx)/ &
                (prs_bar(lev_idx)*slr_zen_ngl_cos) &
                )
           opt_dep_ITOD_O2(lev_idx)=(float_foo-1.0)* &
                0.5*S_d_abs_cff_mss_O2(bnd_idx)*prs_bar(lev_idx)/ &
                S_p_abs_cff_mss_O2(bnd_idx)
        enddo               ! end loop over lev
        
        ! Compute layer "monochromatic" absorption optical depth 
        ! directly from the difference of the column "monochromatic"
        ! optical depths rather than from the log of the difference of
        ! the column transmissions in order to avoid numerical 
        ! difficulties with saturated lines when the transmission is 0.0
        ! This log transmission method fails in single precision.
        do lev_idx=1,lev_nbr
           odal_O2(lev_idx)=odal_O2(lev_idx)+ &
                slr_zen_ngl_cos* &
                (opt_dep_ITOD_O2(lev_idx+1)-opt_dep_ITOD_O2(lev_idx))
        enddo               ! end loop over lev
        
        ! Line absorption computation #5: CO2
        ! bnd_idx_CO2 is set in sequence to the the two indices of CO2
        ! data spanning the same wavenumber interval as the combined data.
        ! Since the narrowband data is in order of increasing wavenumber,
        ! decreasing wavelength, the CO2 data to be considered for
        ! this bnd_idx are 2*bnd_idx-1 and 2*bnd_idx.
        do bnd_idx_CO2=2*bnd_idx-1,2*bnd_idx
           
           ! i21 toggles between 1 and 2 (really between 2 and 1).
           ! i21 is 2: lower wavenumber, lower index, higher wavelength
           ! i21 is 1: higher wavenumber, higher index, lower wavelength
           i21=mod(bnd_idx_CO2,2)+1
           do lev_idx=1,lev_nbr
              
              ! Phi and psi weights are mid-layer quantities.
              phi_wgt(lev_idx)= &
                   exp(A_phi_CO2(bnd_idx_CO2)*tpt_dlt_Mlk(lev_idx)+ &
                   B_phi_CO2(bnd_idx_CO2)*tpt_dlt_Mlk_sqr(lev_idx))
              psi_wgt(lev_idx)= &
                   exp(A_psi_CO2(bnd_idx_CO2)*tpt_dlt_Mlk(lev_idx)+ &
                   B_psi_CO2(bnd_idx_CO2)*tpt_dlt_Mlk_sqr(lev_idx))
           enddo            ! end loop over lev
           
           ! U_bar and prs_bar are layer-interface quantities
           ! Compute top interface
           u_bar(1)=q_CO2(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
           prs_bar(1)=q_CO2(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
           
           ! Sum-up integrands
           ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g.,
           ! q_CO2*dp/g = mpl_CO2 ([kg m-2]). Using the stored values
           ! saves lots of floating point operations.
           do lev_idx=2,levp_nbr
              u_bar(lev_idx)=u_bar(lev_idx-1)+ &
                   mpl_CO2(lev_idx-1)*phi_wgt(lev_idx-1)
              prs_bar(lev_idx)=prs_bar(lev_idx-1)+ &
                   mpl_CO2(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
           enddo            ! end loop over lev
           
           ! Normalize pressure-weighted mass path
           do lev_idx=1,levp_nbr
              prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
           enddo            ! end loop over lev
           
           ! Compute absorption optical depth
           ! NB: The transmission being computed below is an interface
           ! quantity. We could also have computed the transmission of
           ! a layer as a mid-layer quantity directly, but BPB advises
           ! against this for subtle reasons.
           do lev_idx=1,levp_nbr
              float_foo= &
                   sqrt( &
                   1.0+4.0*S_p_abs_cff_mss_CO2(bnd_idx_CO2)*u_bar(lev_idx)/ &
                   (prs_bar(lev_idx)*slr_zen_ngl_cos) &
                   )
              opt_dep_ITOD_CO2_hires(lev_idx,i21)=(float_foo-1.0)* &
                   0.5*S_d_abs_cff_mss_CO2(bnd_idx_CO2)*prs_bar(lev_idx)/ &
                   S_p_abs_cff_mss_CO2(bnd_idx_CO2)
           enddo            ! end loop over lev
        enddo               ! end loop over both hires CO2 bands
        
        ! What I call "interface transmission optical depth" (ITOD) is 
        ! the argument in the exponential that determines transmission.
        ! ITOD has a one over slr_zen_ngl_cos dependence built into it by
        ! the Malkmus formula, so ITOD represents full slant path 
        ! optical depth, NOT zenith optical depth. Perhaps ITOD
        ! should be given a diffusivity factor correction for high
        ! scattering optical depths because the diffuse radiation does
        ! not follow the full slant path. ITOD is an interface quantity.
        ! "interface transmission" (IT) is the transmission associated with
        ! ITOD.
        
        ! What I call "layer transmission optical depth" (LTOD) is the 
        ! difference of two adjacent ITODs. Thus, LTOD has slant path
        ! dependence built in and does not need to be divided by the
        ! cosine of the solar zenith angle. BPB defines the LMOD of a 
        ! layer as the cosine of the solar zenith angle times the LTOD.
        ! Thus, LTOD is the MOD except for a factor of slr_zen_ngl_cos.
        ! "layer transmission" (LT) is the transmission associated with 
        ! LTOD.
        
        ! What I call "layer monochromatic optical depth" (LMOD) is the
        ! optical depth which defines the "layer monochromatic transmission"
        ! (LMT) in the standard RT form of transmission being the 
        ! exponential of minus the zenith optical depth over the cosine
        ! of the solar zenith angle. LMOD is a full layer quantity which
        ! BPB has defined to be the product of a LTOD with 
        ! the cosine of the solar zenith angle.
        ! Note that LMT and LT are identical, they are just
        ! defined distinctly to identify their different methods of 
        ! computation: LMT is the exponential of the LMOD over mu, 
        ! while LT is the exponential of LTOD.
        
        ! What I call "average layer monochromatic optical depth" 
        ! (ALMOD) is the result of combining LMODs from adjacent bands. 
        ! One (easy) way to compute the ALMOD is to simply average the two 
        ! adjacent LMODs. This gives a different answer from computing 
        ! ALMOD as BPB recommends: average the LT of the bands to be 
        ! combined, forming an "average layer transmission" (ALT); compute 
        ! the "average layer transmission optical depth" (ALTOD) by taking
        ! the logarithm of the ALT; compute ALMOD by multiplying 
        ! ALT by the cosine of the solar zenith angle.
        
        ! ALMOD only need be computed for CO2 and CH4
        do i21=2,1,-1
           do lev_idx=1,lev_nbr
              opt_dep_LTOD_CO2_hires= &
                   opt_dep_ITOD_CO2_hires(lev_idx+1,i21)- &
                   opt_dep_ITOD_CO2_hires(lev_idx,i21)
              trn_LT_CO2_hires(lev_idx,i21)= &
                   exp(-opt_dep_LTOD_CO2_hires)
           enddo            ! end loop over lev
        enddo               ! end loop over both hires CO2 bands
        
        do lev_idx=1,lev_nbr
           trn_ALT_CO2= &
                0.5*(trn_LT_CO2_hires(lev_idx,1)+ &
                trn_LT_CO2_hires(lev_idx,2))
           odal_CO2(lev_idx)= &
                -slr_zen_ngl_cos*log(trn_ALT_CO2)
        enddo               ! end loop over lev
        
        ! Line absorption computation #6: CH4
        ! bnd_idx_CH4 is set in sequence to the the two indices of CH4
        ! data spanning the same wavenumber interval as the combined data.
        ! Since the narrowband data is in order of increasing wavenumber,
        ! decreasing wavelength, the CH4 data to be considered for
        ! this bnd_idx are 2*bnd_idx-1 and 2*bnd_idx.
        do bnd_idx_CH4=2*bnd_idx-1,2*bnd_idx
           
           ! i21 toggles between 1 and 2 (really between 2 and 1)
           ! i21 is 2: lower wavenumber, lower index, higher wavelength
           ! i21 is 1: higher wavenumber, higher index, lower wavelength
           i21=mod(bnd_idx_CH4,2)+1
           do lev_idx=1,lev_nbr
              
              ! Phi and psi weights are mid-layer quantities.
              phi_wgt(lev_idx)= &
                   exp(A_phi_CH4(bnd_idx_CH4)*tpt_dlt_Mlk(lev_idx)+ &
                   B_phi_CH4(bnd_idx_CH4)*tpt_dlt_Mlk_sqr(lev_idx))
              psi_wgt(lev_idx)= &
                   exp(A_psi_CH4(bnd_idx_CH4)*tpt_dlt_Mlk(lev_idx)+ &
                   B_psi_CH4(bnd_idx_CH4)*tpt_dlt_Mlk_sqr(lev_idx))
           enddo            ! end loop over lev
           
           ! U_bar and prs_bar are layer-interface quantities
           ! Compute top interface
           u_bar(1)=q_CH4(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
           prs_bar(1)=q_CH4(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
           
           ! Sumuup integrands
           ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g.,
           ! q_CH4*dp/g = mpl_CH4 ([kg m-2]). Using the stored values
           ! saves lots of floating point operations.
           do lev_idx=2,levp_nbr
              u_bar(lev_idx)=u_bar(lev_idx-1)+ &
                   mpl_CH4(lev_idx-1)*phi_wgt(lev_idx-1)
              prs_bar(lev_idx)=prs_bar(lev_idx-1)+ &
                   mpl_CH4(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
           enddo            ! end loop over lev
           
           ! Normalize pressure-weighted mass path
           do lev_idx=1,levp_nbr
              prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
           enddo            ! end loop over lev
           
           ! Compute absorption optical depth
           ! NB: The transmission being computed below is an interface
           ! quantity. We could also have computed the transmission of
           ! a layer as a mid-layer quantity directly, but BPB advises
           ! against this for subtle reasons.
           do lev_idx=1,levp_nbr
              float_foo= &
                   sqrt( &
                   1.0+4.0*S_p_abs_cff_mss_CH4(bnd_idx_CH4)*u_bar(lev_idx)/ &
                   (prs_bar(lev_idx)*slr_zen_ngl_cos) &
                   )
              opt_dep_ITOD_CH4_hires(lev_idx,i21)=(float_foo-1.0)* &
                   0.5*S_d_abs_cff_mss_CH4(bnd_idx_CH4)*prs_bar(lev_idx)/ &
                   S_p_abs_cff_mss_CH4(bnd_idx_CH4)
           enddo            ! end loop over lev
        enddo               ! end loop over both hires CH4 bands
        
        ! What i call "interface transmission optical depth" (ITOD) is 
        ! the argument in the exponential that determines transmission.
        ! ITOD has a one over slr_zen_ngl_cos dependence built into it by
        ! the Malkmus formula, so ITOD represents full slant path 
        ! optical depth, NOT zenith optical depth. Perhaps ITOD
        ! should be given a diffusivity factor correction for high
        ! scattering optical depths because the diffuse radiation does
        ! not follow the full slant path. ITOD is an interface quantity.
        ! "interface transmission" (IT) is the transmission associated with
        ! ITOD.
        
        ! What I call "layer transmission optical depth" (LTOD) is the 
        ! difference of two adjacent ITODs. Thus, LTOD has slant path
        ! dependence built in and does not need to be divided by the
        ! cosine of the solar zenith angle. BPB defines the LMOD of a 
        ! layer as the cosine of the solar zenith angle times the LTOD.
        ! Thus, LTOD is the MOD except for a factor of slr_zen_ngl_cos.
        ! "layer transmission" (LT) is the transmission associated with 
        ! LTOD.
        
        ! What I call "layer monochromatic optical depth" (LMOD) is the
        ! optical depth which defines the "layer monochromatic transmission"
        ! (LMT) in the standard RT form of transmission being the 
        ! exponential of minus the zenith optical depth over the cosine
        ! of the solar zenith angle. LMOD is a full layer quantity which
        ! BPB has defined to be the product of a LTOD with 
        ! the cosine of the solar zenith angle.
        ! Note that LMT and LT are identical, they are just
        ! defined distinctly to identify their different methods of 
        ! computation: LMT is the exponential of the LMOD over mu, 
        ! while LT is the exponential of LTOD.
        
        ! What I call "average layer monochromatic optical depth" 
        ! (ALMOD) is the result of combining LMODs from adjacent bands. 
        ! One (easy) way to compute the ALMOD is to simply average the two 
        ! adjacent LMODs. This gives a different answer from computing 
        ! ALMOD as BPB recommends: average the LT of the bands to be 
        ! combined, forming an "average layer transmission" (ALT); compute 
        ! the "average layer transmission optical depth" (ALTOD) by taking
        ! the logarithm of the ALT; compute ALMOD by multiplying 
        ! ALT by the cosine of the solar zenith angle.
        
        ! ALMOD only need be computed for CO2 and CH4.
        do i21=2,1,-1
           do lev_idx=1,lev_nbr
              opt_dep_LTOD_CH4_hires= &
                   opt_dep_ITOD_CH4_hires(lev_idx+1,i21)- &
                   opt_dep_ITOD_CH4_hires(lev_idx,i21)
              trn_LT_CH4_hires(lev_idx,i21)= &
                   exp(-opt_dep_LTOD_CH4_hires)
           enddo            ! end loop over lev
        enddo               ! end loop over both hires CH4 bands
        
        do lev_idx=1,lev_nbr
           trn_ALT_CH4= &
                0.5*(trn_LT_CH4_hires(lev_idx,1)+ &
                trn_LT_CH4_hires(lev_idx,2))
           odal_CH4(lev_idx)= &
                -slr_zen_ngl_cos*log(trn_ALT_CH4)
        enddo               ! end loop over lev
        
     endif                  ! end if band is in line data region
     
     bnd_idx_tmp_O3=bnd_nbr_O3
     if ((bnd_idx_O3>bnd_nbr_pure_O3).and. &
          (bnd_idx>bnd_nbr_non_O3)) then
        ! We are in overlap region where there is both H2O and continuum data
        ! We need to interpolate O3,O2 continuum absorption data from
        ! O3 bins onto the H2O bins. This is messy.
        ! For now we just search for the O3 band which contains the H2O
        ! band center, and assign that O3 cross-section to that H2O band.
120     if ((wvl_ctr(bnd_idx)<wvl_min_O3(bnd_idx_tmp_O3)).or. &
             (wvl_ctr(bnd_idx)>wvl_max_O3(bnd_idx_tmp_O3))) then
           bnd_idx_tmp_O3=bnd_idx_tmp_O3-1
           goto 120
        endif
        
        do lev_idx=1,lev_nbr
           odal_O3(lev_idx)= &
                abs_xsx_O3(bnd_idx_tmp_O3)*npl_O3(lev_idx)
           odal_O2(lev_idx)=odal_O2(lev_idx)+ &
                abs_xsx_O2(bnd_idx_tmp_O3)*npl_O2(lev_idx)
        enddo               ! end loop over lev
        
     endif                  ! endif in overlap region
     
     if (bnd_idx_O3==bnd_nbr_pure_O3) then
        ! Take special care with this band as it is the band we truncated
        ! O3,O2 cross-sections here should be interpolated
        do lev_idx=1,lev_nbr
           odal_O3(lev_idx)= &
                abs_xsx_O3(bnd_idx_O3)*npl_O3(lev_idx)
           odal_O2(lev_idx)=odal_O2(lev_idx)+ &
                abs_xsx_O2(bnd_idx_O3)*npl_O2(lev_idx)
        enddo               ! end loop over lev
     endif                  ! endif band is the splice band
     
     if ((bnd_idx_O3<bnd_nbr_pure_O3).and. &
          (bnd_idx<=bnd_nbr)) then
        ! We are in pure continuum region
        ! Absorption will be due solely to O3, O2
        do lev_idx=1,lev_nbr
           odal_O3(lev_idx)= &
                abs_xsx_O3(bnd_idx_O3)*npl_O3(lev_idx)
           odal_O2(lev_idx)=odal_O2(lev_idx)+ &
                abs_xsx_O2(bnd_idx_O3)*npl_O2(lev_idx)
        enddo               ! end loop over lev
     endif                  ! end if outside H2O data
     
     ! H2OH2O absorption
     ! call odal_H2OH2O_Chy97(wvn_ctr(bnd_idx),lev_nbr,t,mpl_H2O,RH_lqd,odal_H2OH2O,dbg_lvl)
     ! abs_xsx_H2OH2O(bnd_idx)=abs_xsx_H2OH2O_CFT99(wvn_ctr(bnd_idx))
     ! abs_xsx_H2OH2O(bnd_idx)=abs_xsx_H2OH2O_Chy97(wvn_ctr(bnd_idx))
     do lev_idx=1,lev_nbr
        odal_H2OH2O(lev_idx)= &
             abs_xsx_H2OH2O(bnd_idx)*npl_H2OH2O(lev_idx)
     enddo                  ! end loop over lev
     
     ! O2-O2 continuum absorption
     do lev_idx=1,lev_nbr
        odal_O2O2(lev_idx)=abs_xsx_O2O2(bnd_idx)*npl_O2O2(lev_idx)
     enddo                  ! end loop over lev
     
     ! O2-N2 continuum absorption
     ! Compute 02-N2 absorption optical depth based on O2-O2 optical depth
     ! 1.26 micron band O2-O2 has measured FWHM = 182.5 cm-1 and intensity 1.60 times O2-O2 1.06 micron band. 
     ! Restrict scaling of O2-N2 absorption in 1.26 micron band to between 1.2 and 1.35 um.
     ! SPS98 Figure 1a shows O2-N2 absorption cross sections in units of m5/mlc2
     ! Dimension is absorption cross section of O2-N2 per unit concentration of O2
     ! O2:N2 ratio in lower atmosphere is fixed, so simply scale absorption optical
     ! depth of O2-O2 to O2-N2, rather than cross sections.
     if ((wvl_ctr(bnd_idx)>=1.2e-6).and.(wvl_ctr(bnd_idx)<=1.35e-6)) then
        do lev_idx=1,lev_nbr
           ! fxm: Change efficiency factor (0.2) to parameter or cmd_ln_arg
           ! 20080616: find measurements of O2-N2 cross-sections instead of 0.2
           odal_O2N2(lev_idx)=odal_O2O2(lev_idx)*0.2*N2_per_O2 ! O2-N2 absorption efficiency is 20% of O2-O2 SPS98 p. 12
        enddo               ! end loop over lev
     endif                  ! endif 1.26 micron band
     
     ! NO2 continuum absorption
     do lev_idx=1,lev_nbr
        odal_NO2(lev_idx)=abs_xsx_NO2(bnd_idx)*npl_NO2(lev_idx)
     enddo                  ! end loop over lev
     
     ! Compute aerosol scattering/absorption optical depths
     do lev_idx=1,lev_nbr
        odal_aer(lev_idx)= &
             abs_cff_mss_aer(bnd_idx)*mpl_aer(lev_idx)
        odsl_aer(lev_idx)= &
             sca_cff_mss_aer(bnd_idx)*mpl_aer(lev_idx)
     enddo                  ! end loop over lev

     ! Compute background aerosol scattering/absorption optical depths
     do lev_idx=1,lev_nbr
        odal_bga(lev_idx)= &
             abs_cff_mss_bga(bnd_idx)*mpl_bga(lev_idx)
        odsl_bga(lev_idx)= &
             sca_cff_mss_bga(bnd_idx)*mpl_bga(lev_idx)
     enddo                  ! end loop over lev
     
     ! Compute ice scattering/absorption optical depths
     do lev_idx=1,lev_nbr
        odal_ice(lev_idx)= &
             abs_cff_mss_ice(bnd_idx)*mpl_IWP(lev_idx)
        odsl_ice(lev_idx)= &
             sca_cff_mss_ice(bnd_idx)*mpl_IWP(lev_idx)
     enddo                  ! end loop over lev
     
     ! Compute liquid scattering/absorption optical depths
     do lev_idx=1,lev_nbr
        odal_lqd(lev_idx)= &
             abs_cff_mss_lqd(bnd_idx)*mpl_LWP(lev_idx)
        ! [flg] Liquid cloud droplets are pure scatterers
        if (flg_sct_lqd) odal_lqd(lev_idx)=0.0
        odsl_lqd(lev_idx)= &
             sca_cff_mss_lqd(bnd_idx)*mpl_LWP(lev_idx)
     enddo                  ! end loop over lev
     
     ! Compute snow-impurity scattering/absorption optical depths
     do lev_snw_idx=1,lev_snw_nbr
        lev_idx=lev_snw_idx+lev_atm_nbr
        odal_mpr(lev_idx)= &
             abs_cff_mss_mpr(bnd_idx)*mpl_mpr(lev_snw_idx)
        odsl_mpr(lev_idx)= &
             sca_cff_mss_mpr(bnd_idx)*mpl_mpr(lev_snw_idx)
     enddo                  ! end loop over lev_snw
     
     ! Compute snow scattering/absorption optical depths
     do lev_snw_idx=1,lev_snw_nbr
        lev_bnd_snw_idx=1
        if (lev_bnd_snw_nbr > 1) lev_bnd_snw_idx=lev_snw_idx
        lev_idx=lev_snw_idx+lev_atm_nbr
        odal_snw(lev_idx)= &
             abs_cff_mss_snw(bnd_idx,lev_bnd_snw_idx)*mpl_snw(lev_snw_idx)
        odsl_snw(lev_idx)= &
             sca_cff_mss_snw(bnd_idx,lev_bnd_snw_idx)*mpl_snw(lev_snw_idx)
     enddo                  ! end loop over lev_snw
     
     ! Rayleigh scattering optical depth
     ! See Len93 p. 154 for details. See also BrS84 p. 107, GoY89 p. 297
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+.459/(1.0e12*wvl_ctr(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
     float_foo=32.0*(idx_rfr_air_STP(bnd_idx)-1.0)**2/3.0
     do lev_idx=1,lev_nbr
        ! Following line causes overflow in single precision if not handled properly
        ! Moving factor of N_STP around avoids this problem
        sca_cff_mss_Ray(lev_idx)= &
             pi**3*float_foo*(Avagadro/N_STP)/ &
             (mmw_mst_air(lev_idx)*N_STP*wvl_ctr(bnd_idx)**4)
        odsl_Ray(lev_idx)= &
             mpl_mst_air(lev_idx)*sca_cff_mss_Ray(lev_idx)
     enddo                  ! end loop over lev
     ! End continuum processes
     
     ! All individual optical properties for each radiative process
     ! have been computed. Look to see if user has specified any
     ! processes to be ignored, i.e., compute clear-sky case, pure
     ! Rayleigh scattering case, etc. Zero out appropriate 
     ! optical constants then proceed normally to calculation 
     ! of bulk layer properties. 
     ! NB: Mie asymmetry parameter never needs to be zeroed 
     ! (e.g., in clear sky case) because it is always weighted 
     ! by Mie optical depth before being used in a computation.
     if (.not.flg_Rayleigh) then
        do lev_idx=1,lev_nbr
           odsl_Ray(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no Rayleigh processes
     if (.not.flg_ice) then
        do lev_idx=1,lev_nbr
           odal_ice(lev_idx)=0.0
           odsl_ice(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no ice processes
     if (.not.flg_lqd) then
        do lev_idx=1,lev_nbr
           odal_lqd(lev_idx)=0.0
           odsl_lqd(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no liq processes
     if (.not.flg_aer) then
        do lev_idx=1,lev_nbr
           odal_aer(lev_idx)=0.0
           odsl_aer(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no aer processes
     if (.not.flg_mpr) then
        do lev_idx=1,lev_nbr
           odal_mpr(lev_idx)=0.0
           odsl_mpr(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no mpr processes
     if (.not.flg_snw) then
        do lev_idx=1,lev_nbr
           odal_snw(lev_idx)=0.0
           odsl_snw(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no snw processes
     if (.not.flg_bga) then
        do lev_idx=1,lev_nbr
           odal_bga(lev_idx)=0.0
           odsl_bga(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no bga processes
     if (.not.flg_H2O) then
        do lev_idx=1,lev_nbr
           odal_H2O(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no H2O processes
     if (.not.flg_H2OH2O) then
        do lev_idx=1,lev_nbr
           odal_H2OH2O(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no H2OH2O processes
     if (.not.flg_CO2) then
        do lev_idx=1,lev_nbr
           odal_CO2(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no CO2 processes
     if (.not.flg_OH) then
        do lev_idx=1,lev_nbr
           odal_OH(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no OH processes
     if (.not.flg_CH4) then
        do lev_idx=1,lev_nbr
           odal_CH4(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no CH4 processes
     if (.not.flg_O2) then
        do lev_idx=1,lev_nbr
           odal_O2(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no O2 processes
     if (.not.flg_O3) then
        do lev_idx=1,lev_nbr
           odal_O3(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no O3 processes
     if (.not.flg_O2O2) then
        do lev_idx=1,lev_nbr
           odal_O2O2(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no O2-O2 processes
     if (.not.flg_O2N2) then
        do lev_idx=1,lev_nbr
           odal_O2N2(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no O2N2 processes
     if (.not.flg_NO2) then
        do lev_idx=1,lev_nbr
           odal_NO2(lev_idx)=0.0
        enddo               ! end loop over lev
     endif                  ! end if no NO2 processes
     
     ! Now that processes have been turned on/off, save diagnostic values
     if (bnd_idx==bnd_obs_aer) then
        do lev_idx=1,lev_nbr
           odal_obs_aer(lev_idx)=odal_aer(lev_idx)
           odsl_obs_aer(lev_idx)=odsl_aer(lev_idx)
           odxl_obs_aer(lev_idx)=odal_aer(lev_idx)+odsl_aer(lev_idx)
        enddo               ! end loop over lev
     endif                  ! end if bnd_obs_aer
     if (bnd_idx==bnd_obs_mpr) then
        odxc_obs_mpr=0.0
        do lev_idx=1,lev_nbr
           odal_obs_mpr(lev_idx)=odal_mpr(lev_idx)
           odsl_obs_mpr(lev_idx)=odsl_mpr(lev_idx)
           odxl_obs_mpr(lev_idx)=odal_mpr(lev_idx)+odsl_mpr(lev_idx)
           odxc_obs_mpr=odxc_obs_mpr+odxl_obs_mpr(lev_idx)
        enddo               ! end loop over lev
     endif                  ! end if bnd_obs_mpr
     if (bnd_idx==bnd_obs_snw) then
        odxc_obs_snw=0.0
        do lev_idx=1,lev_nbr
           odal_obs_snw(lev_idx)=odal_snw(lev_idx)
           odsl_obs_snw(lev_idx)=odsl_snw(lev_idx)
           odxl_obs_snw(lev_idx)=odal_snw(lev_idx)+odsl_snw(lev_idx)
           odxc_obs_snw=odxc_obs_snw+odxl_obs_snw(lev_idx)
        enddo               ! end loop over lev
     endif                  ! end if bnd_obs_snw
     if (bnd_idx==bnd_obs_bga) then
        do lev_idx=1,lev_nbr
           odal_obs_bga(lev_idx)=odal_bga(lev_idx)
           odsl_obs_bga(lev_idx)=odsl_bga(lev_idx)
           odxl_obs_bga(lev_idx)=odal_bga(lev_idx)+odsl_bga(lev_idx)
        enddo               ! end loop over lev
     endif                  ! end if bnd_obs_bga
     
     ! Place scattering into correct phase function
     do lev_idx=1,lev_nbr
        odsl_HG(lev_idx)=0.0 ! [frc] Scattering approximated with HG
        odsl_Mie(lev_idx)=0.0 ! [frc] Scattering treated with Mie phase function
        odsl_tmp= & ! [frc] Scattering optical depth
             odsl_aer(lev_idx)+odsl_bga(lev_idx)+ &
             odsl_lqd(lev_idx)+odsl_ice(lev_idx)+ &
             odsl_mpr(lev_idx)+odsl_snw(lev_idx)
        ! Treat scattering with which phase function? Mie or HG
        if (flg_mie) then
           odsl_Mie(lev_idx)=odsl_tmp 
        else 
           odsl_HG(lev_idx)=odsl_tmp
        endif ! !flg_mie
     enddo ! end loop over lev

     ! Weight optical parameters as per CCY83
     do lev_idx=1,lev_nbr
        odsl_spc_ttl(bnd_idx,lev_idx)= &
             odsl_Ray(lev_idx)+ &
             odsl_aer(lev_idx)+odsl_bga(lev_idx)+ &
             odsl_lqd(lev_idx)+odsl_ice(lev_idx)+ &
             odsl_mpr(lev_idx)+odsl_snw(lev_idx)
        odal_spc_ttl(bnd_idx,lev_idx)= &
             odal_H2O(lev_idx)+odal_CO2(lev_idx)+ &
             odal_O2(lev_idx)+odal_O3(lev_idx)+ &
             odal_NO2(lev_idx)+odal_OH(lev_idx)+odal_CH4(lev_idx)+ &
             odal_O2O2(lev_idx)+odal_O2N2(lev_idx)+ &
             odal_H2OH2O(lev_idx)+ &
             odal_aer(lev_idx)+odal_bga(lev_idx)+ &
             odal_lqd(lev_idx)+odal_ice(lev_idx)+ &
             odal_mpr(lev_idx)+odal_snw(lev_idx)
        odxl_spc_ttl(bnd_idx,lev_idx)= &
             odsl_spc_ttl(bnd_idx,lev_idx)+odal_spc_ttl(bnd_idx,lev_idx)
     enddo                  ! end loop over lev
     
     do lev_idx=1,lev_nbr
        ! Initialize some properties to zero
        ! These are over-written later in the loop if they are non-zero
        asm_prm_HG_ttl(bnd_idx,lev_idx)=0.0
        if (flg_mie) then
           do mmn_idx=1,mmn_nbr
              lgn_xpn_cff_Mie_ttl(mmn_idx,bnd_idx,lev_idx)=0.0
           end do ! end loop over mmn
        end if ! !flg_mie

        ! Avoid divide-by-zero conditions in rare, diagnostic cases where ...
        if (odsl_spc_ttl(bnd_idx,lev_idx)<=0.0) then
           ! ... All scattering is turned off ... 
           sca_frc_Ray(lev_idx)=0.0
           sca_frc_HG(lev_idx)=0.0
           sca_frc_Mie(lev_idx)=0.0
           ss_alb_fct(bnd_idx,lev_idx)=0.0
        else if (odsl_spc_ttl(bnd_idx,lev_idx)==odsl_Ray(lev_idx)) then
           ! ... All scattering is Rayleigh scattering ...
           sca_frc_Ray(lev_idx)=1.0
           sca_frc_HG(lev_idx)=0.0
           sca_frc_Mie(lev_idx)=0.0
           ss_alb_fct(bnd_idx,lev_idx)=1.0
        else ! endif no scattering whatsoever
           ! ... Some scattering is particle scattering ...
           ! Single scattering albedo is ill-conditioned when 
           ! scattering optical depth is zero and there is no absorption.
           ! This may only be a problem in single precision---I am not sure
           ss_alb_fct(bnd_idx,lev_idx)= &
                odsl_spc_ttl(bnd_idx,lev_idx)/ &
                odxl_spc_ttl(bnd_idx,lev_idx)
           
           ! Weighted asymmetry parameters and scattering fraction due to 
           ! each process (Rayleigh, HG, and Mie) must be saved for all 
           ! bands and levels until DISORT() is called.
           ! They are used to weight moments of total phase function 
           ! between Rayleigh, Henyey-Greenstein, and Mie components.
           sca_frc_Ray(lev_idx)=odsl_Ray(lev_idx)/odsl_spc_ttl(bnd_idx,lev_idx)
           sca_frc_HG(lev_idx)=odsl_HG(lev_idx)/odsl_spc_ttl(bnd_idx,lev_idx)
           sca_frc_Mie(lev_idx)=odsl_Mie(lev_idx)/odsl_spc_ttl(bnd_idx,lev_idx)
              
              ! See CZP III p. #115 for discussion of effective asymmetry parameter
           lev_bnd_snw_idx=lev_idx-lev_atm_nbr
           if (lev_bnd_snw_nbr==1 .or. lev_bnd_snw_idx < 1) lev_bnd_snw_idx=1

           ! Always compute the weighted asymmetry parameter 
           ! mie computes asm_prm independently of lgn_xpn_cff(1)
           ! lgn_xpn_cff(1) depends on the angular discretization
           ! mie computes asm_prm directly from the a + b mie expansion terms
           ! Hence asm_prm is more accurate than lgn_xpn_cff(1)
           ! asm_prm can therefore be used in place of the first Legendre moment
           ! asm_prm can also be compared to the diagnosed first Legendre moment
           asm_prm_HG_ttl(bnd_idx,lev_idx)= &
                (asm_prm_aer(bnd_idx)*odsl_aer(lev_idx)+ &
                asm_prm_bga(bnd_idx)*odsl_bga(lev_idx)+ &
                asm_prm_lqd(bnd_idx)*odsl_lqd(lev_idx)+ &
                asm_prm_ice(bnd_idx)*odsl_ice(lev_idx)+ &
                asm_prm_mpr(bnd_idx)*odsl_mpr(lev_idx)+ &
                asm_prm_snw(bnd_idx,lev_bnd_snw_idx)*odsl_snw(lev_idx))
           if (flg_mie) then
              asm_prm_HG_ttl(bnd_idx,lev_idx)= &
                   asm_prm_HG_ttl(bnd_idx,lev_idx)/odsl_Mie(lev_idx)
              ! Line above normalizes asm_prm_HG_ttl by odsl_Mie
              ! This ensures asm_prm_HG_ttl is defined when flg_mie is true
           else
              asm_prm_HG_ttl(bnd_idx,lev_idx)= &
                   asm_prm_HG_ttl(bnd_idx,lev_idx)/odsl_HG(lev_idx)
           endif ! !flg_mie

           if (flg_mie) then
              ! Use Legendre expansion for (currently) all particles
              do mmn_idx=1,mmn_nbr
                 lgn_xpn_cff_Mie_ttl(mmn_idx,bnd_idx,lev_idx)= &
                      (lgn_xpn_cff_aer(mmn_idx,bnd_idx)*odsl_aer(lev_idx)+ &
                      lgn_xpn_cff_bga(mmn_idx,bnd_idx)*odsl_bga(lev_idx)+ &
                      lgn_xpn_cff_ice(mmn_idx,bnd_idx)*odsl_ice(lev_idx)+ &
                      lgn_xpn_cff_lqd(mmn_idx,bnd_idx)*odsl_lqd(lev_idx)+ &
                      lgn_xpn_cff_mpr(mmn_idx,bnd_idx)*odsl_mpr(lev_idx)+ &
                      lgn_xpn_cff_snw(mmn_idx,bnd_idx,lev_bnd_snw_idx)*odsl_snw(lev_idx))/ &
                      odsl_Mie(lev_idx)
              end do ! end loop over mmn
           end if ! !flg_mie_snw
        endif ! endif there is particulate scattering
     enddo ! end loop over lev
     
     ! Sanity check for unphysical single scattering albedos
     ! omega's as large as 1.0000656 can occur (at least under LINUX) for very
     ! small Rayleigh scattering optical depths (p <~ 1 Pa) in single precision
     do lev_idx=1,lev_nbr
        if (ss_alb_fct(bnd_idx,lev_idx)>1.0) then
           write (6,'(a,a,i4,a,i3,a,f10.7)')  &
                prg_nm(1:ftn_strlen(prg_nm)),': WARNING ss_alb_fct(',bnd_idx,',',lev_idx,') = ',ss_alb_fct(bnd_idx,lev_idx)
           ss_alb_fct(bnd_idx,lev_idx)=1.0
        endif
     enddo                  ! end loop over lev
     
     ! Compute diagnostics
     do lev_idx=1,lev_nbr
        ! Gaseous absorption
        odxc_spc_CH4(bnd_idx)=odxc_spc_CH4(bnd_idx)+odal_CH4(lev_idx)
        odxc_spc_CO2(bnd_idx)=odxc_spc_CO2(bnd_idx)+odal_CO2(lev_idx)
        odxc_spc_H2O(bnd_idx)=odxc_spc_H2O(bnd_idx)+odal_H2O(lev_idx)
        odxc_spc_H2OH2O(bnd_idx)=odxc_spc_H2OH2O(bnd_idx)+odal_H2OH2O(lev_idx)
        odxc_spc_NO2(bnd_idx)=odxc_spc_NO2(bnd_idx)+odal_NO2(lev_idx)
        odxc_spc_O2(bnd_idx)=odxc_spc_O2(bnd_idx)+odal_O2(lev_idx)
        odxc_spc_O2N2(bnd_idx)=odxc_spc_O2N2(bnd_idx)+odal_O2N2(lev_idx)
        odxc_spc_O2O2(bnd_idx)=odxc_spc_O2O2(bnd_idx)+odal_O2O2(lev_idx)
        odxc_spc_O3(bnd_idx)=odxc_spc_O3(bnd_idx)+odal_O3(lev_idx)
        odxc_spc_OH(bnd_idx)=odxc_spc_OH(bnd_idx)+odal_OH(lev_idx)
        ! Gaseous scattering
        odxc_spc_Ray(bnd_idx)=odxc_spc_Ray(bnd_idx)+odsl_Ray(lev_idx)
        ! Particle absorption
        odac_spc_aer(bnd_idx)=odac_spc_aer(bnd_idx)+odal_aer(lev_idx)
        odac_spc_bga(bnd_idx)=odac_spc_bga(bnd_idx)+odal_bga(lev_idx)
        odac_spc_ice(bnd_idx)=odac_spc_ice(bnd_idx)+odal_ice(lev_idx)
        odac_spc_lqd(bnd_idx)=odac_spc_lqd(bnd_idx)+odal_lqd(lev_idx)
        odac_spc_mpr(bnd_idx)=odac_spc_mpr(bnd_idx)+odal_mpr(lev_idx)
        odac_spc_snw(bnd_idx)=odac_spc_snw(bnd_idx)+odal_snw(lev_idx)
        ! Particle extinction
        odxc_spc_aer(bnd_idx)=odxc_spc_aer(bnd_idx)+odal_aer(lev_idx)+odsl_aer(lev_idx)
        odxc_spc_bga(bnd_idx)=odxc_spc_bga(bnd_idx)+odal_bga(lev_idx)+odsl_bga(lev_idx)
        odxc_spc_ice(bnd_idx)=odxc_spc_ice(bnd_idx)+odal_ice(lev_idx)+odsl_ice(lev_idx)
        odxc_spc_lqd(bnd_idx)=odxc_spc_lqd(bnd_idx)+odal_lqd(lev_idx)+odsl_lqd(lev_idx)
        odxc_spc_mpr(bnd_idx)=odxc_spc_mpr(bnd_idx)+odal_mpr(lev_idx)+odsl_mpr(lev_idx)
        odxc_spc_snw(bnd_idx)=odxc_spc_snw(bnd_idx)+odal_snw(lev_idx)+odsl_snw(lev_idx)
     enddo                  ! end loop over lev
     odxc_spc_ttl(bnd_idx)= &
          ! Gaseous absorption
          odxc_spc_CH4(bnd_idx)+ &
          odxc_spc_CO2(bnd_idx)+ &
          odxc_spc_H2O(bnd_idx)+ &
          odxc_spc_H2OH2O(bnd_idx)+ &
          odxc_spc_NO2(bnd_idx)+ &
          odxc_spc_O2(bnd_idx)+ &
          odxc_spc_O2N2(bnd_idx)+ &
          odxc_spc_O2O2(bnd_idx)+ &
          odxc_spc_O3(bnd_idx)+ &
          odxc_spc_OH(bnd_idx)+ &
          ! Gaseous scattering
          odxc_spc_Ray(bnd_idx)+ &
          ! Particle extinction
          odxc_spc_aer(bnd_idx)+ &
          odxc_spc_bga(bnd_idx)+ &
          odxc_spc_ice(bnd_idx)+ &
          odxc_spc_lqd(bnd_idx)+ &
          odxc_spc_mpr(bnd_idx)+ &
          odxc_spc_snw(bnd_idx)
     
     ! Spectral information needed for DISORT and analysis has now been saved
     ! Massage data for input to DISORT
     ! DISORT uses a top down numbering scheme like CCM2, CCM3, CAM
     ! (Reverse layer ordering from that of HPMM cirrus cloud model)
     
     ! Assign optical depth and single scattering albedos
     do lev_idx=1,lev_nbr
        dtauc(lev_idx)=odxl_spc_ttl(bnd_idx,lev_idx)
        ssalb(lev_idx)=ss_alb_fct(bnd_idx,lev_idx)
     end do ! end loop over lev
     
     ! Phase-function information
     do lev_idx=1,lev_nbr
        pmom(0,lev_idx)=1.0
        ! fxm: Until mie computes PMOM correctly use g for PMOM(1)
        pmom(1,lev_idx)=asm_prm_HG_ttl(bnd_idx,lev_idx)
     end do ! end loop over lev
     
     ! All moments of Rayleigh phase-function except second are zero
     ! To blend Rayleigh phase-function with HG and Mie scattering phase-functions,
     ! weight coefficient contribution by scattering fraction due to its process
     do mmn_idx=2,mmn_nbr
        do lev_idx=1,lev_nbr
           ! Combine HG with "exact" (for spheres) Mie phase functions
           if (flg_mie) then
              pmom(mmn_idx,lev_idx)= &
                   sca_frc_Mie(lev_idx)*lgn_xpn_cff_Mie_ttl(mmn_idx,bnd_idx,lev_idx)
           else ! !flg_mie
              pmom(mmn_idx,lev_idx)= &
                   sca_frc_HG(lev_idx)*asm_prm_HG_ttl(bnd_idx,lev_idx)**mmn_idx
           end if ! !flg_mie
        end do ! end loop over lev
     end do ! end loop over mmn
     
     ! Add Rayleigh-scattered fraction times Rayleigh expansion
     do lev_idx=1,lev_nbr
        ! Rayleigh scattering contributes only to second moment
        pmom(2,lev_idx)=pmom(2,lev_idx)+sca_frc_Ray(lev_idx)*0.1
     end do

     if (.not.sv_cmp_tau) then 
        
        ! If usrtau==true. then place levels at top and bottom of atmosphere 
        ! and at computational levels in between. 
        ! This is same as setting usrtau=.false.
        ! There is not much point (yet) setting user levels to something else, e.g., 
        ! evenly spaced in optical depth, or at cloud top and bottom of cloud. 
        ! However, a reason may arise someday 
        ! In that case all radiant quantities output to netCDF must be 
        ! dimensioned by tau_nbr instead of levp_nbr.
        tau(1)=0.0
        do tau_idx=2,tau_nbr
           tau(tau_idx)=tau(tau_idx-1)+odxl_spc_ttl(bnd_idx,tau_idx-1)
        end do              ! end loop over tau
        
        ntau=tau_nbr
        do tau_idx=1,tau_nbr
           utau(tau_idx)=tau(tau_idx)
        end do              ! end loop over tau
     else
        do tau_idx=1,tau_nbr
           tau(tau_idx)=utau(tau_idx)
        end do              ! end loop over tau
     endif                  ! end if setting user defined levels
     
     ! Set wavenumbers for this spectral interval.
     wvnmlo=wvn_min(bnd_idx)
     wvnmhi=wvn_max(bnd_idx)
     
     ! Set plank=.true. whenever considering thermal emission
     ! Computing thermal source function generates underflows in single precision when lambda < XXX um
     if (wvl_ctr(bnd_idx)>wvl_Planck) then
        plank=flg_Planck
     else 
        plank=.false.
     endif
     
     ! Intensity of incident collimated beam at top boundary
     ! Units are arbitrary though must match fisot in solar case
     ! Must be in W m-2 if there is any thermal emission (plank=.true.)
     ! fbeam is flux normal to earth-sun path, not normal to ground!
     ! Spectral integral of fbeam at TOA is "solar constant", ~1367 W m-2
     ! DISORT interally adjusts incoming flux for zenith angle
     ! DO NOT pre-multiply fbeam by solar zenith angle cosine
     fbeam=slr_cst_xnt_fac*flx_slr_frc(bnd_idx)*flx_frc_drc_TOA
     
     ! Incident isotropic intensity at top boundary
     ! DISORT does not (and should not) apply zenith angle corrections to fisot
     ! Hence changing problem from (non-zenith) direct to diffuse
     ! using flx_frc_drc_TOA causes different total downwelling flux
     ! Units of FISOT and factor of pi are subtle---from disort2.doc:
     ! "FISOT: Intensity of top-boundary isotropic illumination.
     ! (same units as PLKAVG (default W/sq m) if thermal sources active, 
     ! otherwise arbitrary units).
     ! Corresponding incident flux is pi (3.14159...) (steradians) times FISOT."
     ! In other words, DISORT wants fisot to be input in radiance units
     ! Confusingly, DISORT documentation treats steradians as non-dimensional
     ! i.e., states units of PLKAVG as W m-2 instead of W m-2 sr-1
     ! Two clues that DISORT does in fact want fisot input as radiance: 
     ! 1. Units are described as "intensity"
     ! 2. Instruction to multiply by pi to convert isotropic fisot to horizontal flux
     ! Next line converts horizontal diffuse flux to isotropic radiance with 1/pi
     fisot=slr_cst_xnt_fac*flx_slr_frc(bnd_idx)*(1.0-flx_frc_drc_TOA)/pi
     ! test for JGG 20080925
     ! fisot=0.2*fbeam
  
     ! Give details about lower boundary reflectance
     ! When lamber is true, albedo specifies isotropic reflectance
     ! Otherwise hl array must be specified to give BRDF of bottom boundary
     ! lamber=.true. DISORT2
     if (flg_rfl) then
        ! Apply spectral reflectance from boundary condition file
        albedo=rfl_spc_sfc(bnd_idx)
        if(albedo<0.0.or.albedo>1.0) stop "rfl_spc_sfc out of bounds"
     else ! flg_rfl
        ! Determine spectral reflectance from broadband parameters
        if (wvl_ctr(bnd_idx)<0.7e-6) then
           if (slr_zen_ngl_cos>0.5) then
              albedo=alb_sfc_vsb_drc
           else
              albedo=alb_sfc_vsb_dff
           endif
        else
           if (slr_zen_ngl_cos>0.5) then
              albedo=alb_sfc_NIR_drc
           else
              albedo=alb_sfc_NIR_dff
           endif
        endif
     endif ! !flg_rfl
     ! Ocean albedo
     if (.not. mode_std .and. ocn_msk == 1.0) then
        albedo = 0.026 / (slr_zen_ngl_cos**1.7 + 0.065) + &
             0.15 * (slr_zen_ngl_cos - 0.1) * (slr_zen_ngl_cos - 0.5) &
             * (slr_zen_ngl_cos - 1.0)
     endif
     ! Archive spectral reflectance actually used
     rfl_spc_sfc(bnd_idx)=albedo
     
     if (single_bnd_computation) then
        
        ! Set header which DISORT will use printing results
        ! If header has length greater than zero, DISORT() prints
        ! annoying header message on each call.
        write (header,'(a,i4,a,f6.4,a,es8.1,a,es8.1)') &
             'swnb: i = ',bnd_idx,', lambda = ', &
             wvl_ctr(bnd_idx)*1.0e6, &
             ', flx_slr_frc = ',flx_slr_frc(bnd_idx),', bandwith = ', &
             wvl_dlt(bnd_idx)
        
        if (bnd_idx==bnd_dbg) then
           
           ! Hardcode numbers to test radiative code against benchmark values
           if (tst_case_Rayleigh) then
              ! Pure Rayleigh scattering
              ! swnb2 -A -B -C -G -H -I -J -K -L -O -U -W -X -Y -t -s 16 -D 1 -r 0.0 -E -e 1 -p ${DATA}/aca/mls_clr.nc -d ~/swnb.nc
              ! Answer: Albedo = 0.796920
              fbeam=pi
              fisot=0.0
              albedo=0.0
              umu0=0.5
              do lev_idx=1,lev_nbr
                 dtauc(lev_idx)=4.284/real(lev_nbr)
                 odsl_Ray(lev_idx)=dtauc(lev_idx)
                 odsl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                 odxl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
              enddo         ! end loop over lev
           endif            ! end if tst_case_Rayleigh
           if (tst_case_HG) then
              ! Pure Henyey-Greenstein scattering
              ! Answers: Albedo = 0.123420, Atmospheric absorptance = 0.360522
              ! swnb2 -A -B -C -G -H -I -J -K -L -O -R -U -W -X -Y -P -s 16 -D 1 -r 0.0 -E -e 1 -p ${DATA}/aca/mls_clr.nc -d ~/swnb.nc
              ! Check answers with:
              ! ncks -F -C -d bnd,1 -H -v abs_spc_atm,rfl_spc_SAS,rfl_spc_sfc,trn_spc_atm_ttl ~/swnb.nc
              fbeam=pi
              fisot=0.0
              albedo=0.0
              umu0=0.5
              do lev_idx=1,lev_nbr
                 dtauc(lev_idx)=1.0/real(lev_nbr)
                 ssalb(lev_idx)=0.8
                 odsl_HG(lev_idx)=dtauc(lev_idx)
                 odsl_spc_ttl(bnd_idx,lev_idx)=odsl_HG(lev_idx)
                 odxl_spc_ttl(bnd_idx,lev_idx)=odsl_HG(lev_idx)
                 ss_alb_fct(bnd_idx,lev_idx)=ssalb(lev_idx)
                 sca_frc_HG(lev_idx)=1.0
                 sca_frc_Mie(lev_idx)=0.0
                 asm_prm_lqd(bnd_idx)=0.75
                 asm_prm_HG_ttl(bnd_idx,lev_idx)=asm_prm_lqd(bnd_idx)
              enddo         ! end loop over lev
              do lev_idx=1,lev_nbr
                 pmom(0,lev_idx)=1.0
              end do
              do mmn_idx=1,mmn_nbr
                 do lev_idx=1,lev_nbr
                    pmom(mmn_idx,lev_idx)=0.75**mmn_idx
                 end do
              end do
           endif            ! end if tst_case_HG
           if (dbg_lvl>dbg_off) then
              write (6,'(/,a,i4,a4,f9.6,a2,f9.6,a11,f9.3,a2,f9.3,a5)') &
                   'Detailed debugging information for bnd(',bnd_idx,') = ', &
                   wvl_min(bnd_idx)*1.0e6,'--',wvl_max(bnd_idx)*1.0e6,' um = ', &
                   wvnmlo,'--',wvnmhi,' cm-1'
              write (6,'(a4,20(1x,a8))') &
                   'lev_idx', &
                   'prs_mdp', &
                   'tpt_mdp', &
                   'tau_sca', &
                   'tau_abs', &
                   'tau_ext', &
                   'ss_alb', &
                   'frc_HG', &
                   'frc_Mie', &
                   'sct_Ray', &
                   'sct_HG', &
                   'sct_Mie', &
                   'asm_HG', &
                   'pth_wet', &
                   'dtauc', &
                   'ssalb', &
                   'pmom(1)', &
                   'pmom(2)', &
                   'pmom(3)', &
                   'pmom(4)'
              do lev_idx=1,lev_nbr
                 write (6,'(i4,20(1x,es8.1))') &
                      lev_idx, &
                      prs(lev_idx), &
                      tpt(lev_idx), &
                      odsl_spc_ttl(bnd_idx,lev_idx), &
                      odal_spc_ttl(bnd_idx,lev_idx), &
                      odxl_spc_ttl(bnd_idx,lev_idx), &
                      ss_alb_fct(bnd_idx,lev_idx), &
                      sca_frc_HG(lev_idx), &
                      sca_frc_Mie(lev_idx), &
                      odsl_Ray(lev_idx), &
                      odsl_HG(lev_idx), &
                      odsl_Mie(lev_idx), &
                      asm_prm_HG_ttl(bnd_idx,lev_idx), &
                      mpl_mst_air(lev_idx), &
                      dtauc(lev_idx), &
                      ssalb(lev_idx), &
                      pmom(1,lev_idx), &
                      pmom(2,lev_idx), &
                      pmom(3,lev_idx), &
                      pmom(4,lev_idx)
              enddo         ! end loop over lev
           endif            ! endif dbg
        else                ! end if dbg band
           goto 999
        endif               ! end if not dbg band
     endif                  ! end if single_bnd_computation
     
     ! Snow pre-processing
     if_flg_msm: if (flg_msm) then

        ! CCY83 Appendix B
        do lev_snw_idx=1,lev_snw_nbr
           lev_idx=lev_snw_idx+lev_atm_nbr
           ! Copy base parameters from DISORT input
           tau_ext=odxl_spc_ttl(bnd_idx,lev_idx)
           sng_sct=ss_alb_fct(bnd_idx,lev_idx)
           asm_prm=pmom(1,lev_idx)/3.0 ! CCY83 (B9)
           ! Fraction of energy in delta-forward-scattered peak
           fwd_frc=2.0*pmom(2,lev_idx)/5.0 ! CCY83 (B10)
           ! Delta-scaled single-scattering parameters
           tau_ext_scl=(1.0-sng_sct*fwd_frc)*tau_ext ! CCY83 (B6), CAM3 (4.208)
           sng_sct_scl=(1.0-fwd_frc)*sng_sct/(1.0-sng_sct*fwd_frc) ! CCY83 (B7), CAM3 (4.209)
           asm_prm_scl=(asm_prm-fwd_frc)/(1.0-fwd_frc) ! CCY83 (B8), CAM3 (4.210)
           ! Parameters for delta-Eddington diffuse R/T
           lambda_CCY83=sqrt(3.0*(1.0-sng_sct_scl)*(1.0-sng_sct_scl*asm_prm_scl)) ! CCY83 (B5), CAM3 (4.219)
           ! Parameters for delta-Eddington direct R/T
           mu_CCY83=slr_zen_ngl_cos
           alpha_CCY83=0.75*sng_sct_scl*mu_CCY83* & ! CCY83 (B3) CAM3 (4.215)
                (1.0+asm_prm_scl*(1.0-sng_sct_scl))/ &
                (1.0-lambda_CCY83*lambda_CCY83*mu_CCY83*mu_CCY83)
           gamma_CCY83=0.5*sng_sct_scl* & ! CCY83 (B4) CAM3 (4.216)
                (3.0*asm_prm_scl*(1.0-sng_sct_scl)*mu_CCY83*mu_CCY83+1.0)/ &
                (1.0-lambda_CCY83*lambda_CCY83*mu_CCY83*mu_CCY83)
           ! Layer isotropic reflectance/transmittance r/t
           if (sng_sct<1.0) then
              u_CCY83=1.5*(1.0-sng_sct_scl*asm_prm_scl)/lambda_CCY83 ! CCY83 (B15), CAM3 (4.218)
              ! NB: eliminate positive exponentials as in TMA89?
              xpn_arg=min(25.0,lambda_CCY83*tau_ext_scl) ! Bri92 limit
              xpn_val=exp(-xpn_arg) ! e^{-lambda*tau_star}
              N_CCY83= & ! CCY83 (B14), CAM3 (4.217)
                   (u_CCY83+1.0)*(u_CCY83+1.0)/xpn_val- &
                   (u_CCY83-1.0)*(u_CCY83-1.0)*xpn_val
              ! [frc] Snow layer spectral flux reflectance
              rfl_dff_spc_snw(bnd_idx,lev_snw_idx)= & ! CCY83 (B12), CAM3 (4.213)
                   (u_CCY83+1.0)*(u_CCY83-1.0)* &
                   (1.0/xpn_val-xpn_val)/N_CCY83
              if(rfl_dff_spc_snw(bnd_idx,lev_snw_idx)<0.0) then 
                 write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                      ': INFO negative isotropic reflectance in D/E approximation'
                 rfl_dff_spc_snw(bnd_idx,lev_snw_idx)=0.0
              endif ! endif error
              ! [frc] Snow layer spectral flux transmittance
              trn_dff_spc_snw(bnd_idx,lev_snw_idx)=4.0*u_CCY83/N_CCY83 ! CCY83 (B13), CAM3 (4.214)
           else ! sng_sct=1.0
              ! [frc] Snow layer spectral flux reflectance
              rfl_dff_spc_snw(bnd_idx,lev_snw_idx)= & ! CCY83 (B16)
                   (1.0-asm_prm)*tau_ext/ &
                   ((1.0-asm_prm)*tau_ext+4.0/3.0)
              trn_dff_spc_snw(bnd_idx,lev_snw_idx)= & ! CCY83 (B17)
                   1.0-rfl_dff_spc_snw(bnd_idx,lev_snw_idx)
           endif ! sng_sct=1.0
           xpn_arg=min(25.0,tau_ext_scl/mu_CCY83) ! Bri92 limit
           xpn_val=exp(-xpn_arg) ! e^{-tau_star/mu}
           ! R from CCY83/CAM3 is reflectance to direct-incident radiation
           ! T from CCY83/CAM3 is "total transmittance" to direct-incident radiation
           ! NB: CCY83 p. 135 (B1) differs from CAM3 p. 114 (4.211)
           ! Assume CAM3, and CAM3 implementation, are correct and that
           ! CCY83 (B1) should multiply R_diffuse by (alpha+gamma) not (alpha-gamma)
           rfl_drc_spc_snw(bnd_idx,lev_snw_idx)= & ! CCY83 (B1), CAM3 p. 114 (4.211)
                (alpha_CCY83-gamma_CCY83)* &
                (trn_dff_spc_snw(bnd_idx,lev_snw_idx)*xpn_val-1.0)+ &
                (alpha_CCY83+gamma_CCY83)* &
                rfl_dff_spc_snw(bnd_idx,lev_snw_idx)
           ! NB: Transmission to direct beam is total (direct+diffuse) downwelling 
           ! flux exiting layer divided by downwelling direct flux entering layer.
           ! This "total transmission" differs from direct transmission of the
           ! direct beam which is simply Beer's Law.
           trn_ttl_drc_spc_snw(bnd_idx,lev_snw_idx)= & ! CCY83 (B2), CAM3 p. 114 (4.212)
                (alpha_CCY83+gamma_CCY83)*trn_dff_spc_snw(bnd_idx,lev_snw_idx)+ &
                xpn_val* &
                ((alpha_CCY83-gamma_CCY83)*rfl_dff_spc_snw(bnd_idx,lev_snw_idx)- &
                (alpha_CCY83+gamma_CCY83-1.0))
           ! Direct-beam transmittance
           trn_drc_drc(bnd_idx,lev_snw_idx)=xpn_val
        enddo ! end loop over lev_snw

        ! Mathematically treat bottom boundary condition as extra snow layer
        ! Reflectance properties are of surface (ice, ground...) underlying snow
        ! This allows "vectorization" of layer combination loops
        rfl_dff_spc_snw(bnd_idx,levp_snw_nbr)=rfl_spc_sfc(bnd_idx)
        rfl_drc_spc_snw(bnd_idx,levp_snw_nbr)=rfl_spc_sfc(bnd_idx)
        trn_dff_spc_snw(bnd_idx,levp_snw_nbr)=0.0
        trn_ttl_drc_spc_snw(bnd_idx,levp_snw_nbr)=0.0
        trn_drc_drc(bnd_idx,levp_snw_nbr)=0.0

        ! Diagnostic layer properties
        do levp_snw_idx=1,levp_snw_nbr
           rfl_ddm_spc_snw(bnd_idx,levp_snw_idx)= &
                flx_frc_drc_TOA*rfl_drc_spc_snw(bnd_idx,levp_snw_idx)+ &
                (1.0-flx_frc_drc_TOA)*rfl_dff_spc_snw(bnd_idx,levp_snw_idx)
           trn_ddm_spc_snw(bnd_idx,levp_snw_idx)=mss_val ! fxm
        enddo ! end loop over levp_snw

        ! delta-Eddington layer optical properties have now been computed
        ! Initialize inhomogenous (layer combination) properties for adding
        rfl_dff_spc_snw_crr=rfl_dff_spc_snw(bnd_idx,1)
        trn_dff_spc_snw_crr=trn_dff_spc_snw(bnd_idx,1)
        rfl_drc_spc_snw_crr=rfl_drc_spc_snw(bnd_idx,1)
        trn_ttl_drc_spc_snw_crr=trn_ttl_drc_spc_snw(bnd_idx,1)
        trn_drc_drc_crr=trn_drc_drc(bnd_idx,1)
        ! Initialize diagnostic outputs
        rfl_dff_spc_snw_cnt(bnd_idx,1)=rfl_dff_spc_snw_crr
        trn_dff_spc_snw_cnt(bnd_idx,1)=trn_dff_spc_snw_crr

        ! delta-Eddington snowpack (_snp) combined optical properties on interfaces
        ! These are combined optical properties of partial inhomogeneous snowpack 
        ! "Partial snowpack" includes all snow layers (down from)/(up to) interface
        ! Snowpack properties and boundary fluxes determine all interior fluxes
        ! Nomenclature "upw" and "dwn" for snowpack properties follows Bri92/CAM3:
        ! upw reflectance properties of partial snowpack apply to downwelling light
        ! dwn reflectance properties of partial snowpack apply to upwelling   light
        !   transmittance properties of partial snowpack apply to downwelling light

        ! fxm: initialize and compute these five diagnostics
        ! Compute the two "upw" properties in upward loop through layers
        rfl_drc_upw_snp(bnd_idx,1)=mss_val ! CAM3 p. 115 R_dwn(munot)
        rfl_dff_upw_snp(bnd_idx,1)=mss_val ! CAM3 p. 115 R_bar_dwn
        ! Compute the three "dwn" properties in downward loop through layers
        rfl_dff_dwn_snp(bnd_idx,1)=0.0 ! CAM3 p. 115 R_bar_upw
        trn_ttl_drc_snp(bnd_idx,1)=1.0 ! CAM3 p. 115 T_dwn(munot)
        trn_drc_drc_snp(bnd_idx,1)=1.0 ! CAM3 p. 115 T_drc_drc(munot)
        do lev_snw_idx=2,levp_snw_nbr
           ! Fake loop to fill snowpack properties fxm split into upw and dwn loops
           rfl_drc_upw_snp(bnd_idx,lev_snw_idx)=mss_val ! CAM3 p. 115 R_dwn(munot)
           rfl_dff_upw_snp(bnd_idx,lev_snw_idx)=mss_val ! CAM3 p. 115 R_bar_dwn
           rfl_dff_dwn_snp(bnd_idx,lev_snw_idx)=mss_val ! CAM3 p. 115 R_bar_upw
           trn_ttl_drc_snp(bnd_idx,lev_snw_idx)=mss_val ! CAM3 p. 115 T_dwn(munot)
           trn_drc_drc_snp(bnd_idx,lev_snw_idx)=mss_val ! CAM3 p. 115 T_drc_drc(munot)
        enddo                  ! end loop over lev_snw

        ! Proceed top-down through layers computing inhomogeneous R/T to interfaces
        do levp_snw_idx=2,levp_snw_nbr
           tmp_fct_mtx_nvr=1.0/(1.0-rfl_dff_spc_snw(bnd_idx,levp_snw_idx)*rfl_dff_spc_snw_crr)
           ! Diffuse flux addition formulae
           ! CCY83 p. 136 (B21), CAM3 p. 115 (4.222)
           rfl_dff_spc_snw_nxt=rfl_dff_spc_snw_crr+ &
                trn_dff_spc_snw_crr*tmp_fct_mtx_nvr* &
                rfl_dff_spc_snw(bnd_idx,levp_snw_idx)*trn_dff_spc_snw_crr
           ! CCY83 p. 136 (B22), CAM3 p. 115 (4.223)
           trn_dff_spc_snw_nxt=trn_dff_spc_snw_crr*tmp_fct_mtx_nvr*trn_dff_spc_snw(bnd_idx,levp_snw_idx)
           
           ! Direct transmittance of direct beam
           trn_drc_drc_nxt=trn_drc_drc_crr*trn_drc_drc(bnd_idx,levp_snw_idx)
           ! NB: 
           ! CCY83 (1) and (2) are delta-Eddington apx. Figure B1 method 1 
           ! CCY83 (B19) and (B20) are delta-Eddington apx. Figure B1 method 2 
           ! CCY83 (B21) and (B22) are used in two-stream apx. Figure B1
           ! delta-Eddington is superior, so ignore two-stream discussion
           
           ! CCY83 p. 136 (B19), CAM3 p. 115 (4.220)
           rfl_drc_spc_snw_nxt=rfl_drc_spc_snw_crr+trn_dff_spc_snw_crr* &
                ((trn_ttl_drc_spc_snw_crr-trn_drc_drc_crr)*rfl_dff_spc_snw(bnd_idx,levp_snw_idx)+ &
                trn_drc_drc_crr*rfl_drc_spc_snw(bnd_idx,levp_snw_idx))*tmp_fct_mtx_nvr
           ! CCY83 p. 136 (B20), CAM3 p. 115 (4.221)
           trn_ttl_drc_spc_snw_nxt=trn_drc_drc_crr*trn_ttl_drc_spc_snw(bnd_idx,levp_snw_idx)+ &
                trn_dff_spc_snw(bnd_idx,levp_snw_idx)*tmp_fct_mtx_nvr* &
                ((trn_ttl_drc_spc_snw_crr-trn_drc_drc_crr)+ &
                trn_drc_drc_crr*rfl_drc_spc_snw(bnd_idx,levp_snw_idx)*rfl_dff_spc_snw_crr)
           
           rfl_dff_spc_snw_cnt(bnd_idx,levp_snw_idx)=rfl_dff_spc_snw_nxt-rfl_dff_spc_snw_crr
           trn_dff_spc_snw_cnt(bnd_idx,levp_snw_idx)=trn_dff_spc_snw_crr-trn_dff_spc_snw_nxt
           if (bnd_idx==bnd_dbg.and.dbg_lvl>dbg_off) then
              write (6,'(a,i4,4(a25,f9.6))') &
                   'levp_snw = ',levp_snw_idx, &
                   ', rfl_dff_spc_snw_crr = ',rfl_dff_spc_snw_crr, &
                   ', trn_dff_spc_snw_crr = ',trn_dff_spc_snw_crr, &
                   ', rfl_dff_spc_snw_nxt = ',rfl_dff_spc_snw_nxt, &
                   ', trn_dff_spc_snw_nxt = ',trn_dff_spc_snw_nxt
           endif                ! end if dbg band
           rfl_dff_spc_snw_crr=rfl_dff_spc_snw_nxt
           trn_dff_spc_snw_crr=trn_dff_spc_snw_nxt
           rfl_drc_spc_snw_crr=rfl_drc_spc_snw_nxt
           trn_ttl_drc_spc_snw_crr=trn_ttl_drc_spc_snw_nxt
           trn_drc_drc_crr=trn_drc_drc_nxt
        enddo                  ! end loop over levp_snw

        ! Diagnostic surface albedos
        alb_dff_spc_snw_dea(bnd_idx)=rfl_dff_spc_snw_nxt
        alb_drc_spc_snw_dea(bnd_idx)=rfl_drc_spc_snw_nxt

     end if if_flg_msm ! !flg_msm pre-processing
     
     ! call DISORT( nlyr, dtauc, ssalb, pmom, temper, wvnmlo, &
     !         wvnmhi, usrtau, ntau, utau, nstr, usrang, &
     !         numu, umu, nphi, phi, ibcnd, fbeam, umu0, &
     !         phi0, fisot, lamber, albedo, hl, btemp, ttemp, &
     !         temis, deltam, plank, onlyfl, accur, prnt, &
     !         header, maxcly, maxulv, maxumu, maxcmu, &
     !         maxphi, rfldir, rfldn, flup, dfdt, uavg, &
     !         uu, u0u, albmed, trnmed )
     
     ! Call DISORT2
     ! new arguments: nmom, maxmom
     ! obsolete arguments: hl, deltam, maxcmu, u0u
     u0u(:,:)=0.0 ! disort1 u0u held azimuthal average intensity at user angles
     
     ! disort1 (output) u0u held azimuthal average intensity at user angles
     ! disort2 (output) uavg holds mean intensity (azimuthal-and-polar-averaged) including the direct beam
     ! disort2 (internal) u0c: azimuthal average intensity at computational angles
     ! disort2 (internal) u0u: azimuthal average diffuse intensity (ONLYFL=false)
     ! disort2 (output) rfldir: Direct-beam flux (unscaled) 
     ! disort2 (output) rfldn: Diffuse down-flux (total minus direct, unscaled) 
     ! disort2 (output) flup: Diffuse up-flux
     
     ! call t_startf('disort')
     
     if (slr_zen_ngl_cos > 0.0) then ! daytime
        call DISORT( nlyr, dtauc, ssalb, nmom, pmom, temper, wvnmlo, &
             wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, &
             umu, nphi, phi, ibcnd, fbeam, umu0, phi0, &
             fisot, lamber, albedo, btemp, ttemp, temis, &
             plank, onlyfl, accur, prnt, header, maxcly, &
             maxulv, maxumu, maxphi, maxmom, rfldir, rfldn, &
             flup, dfdt, uavg, uu, albmed, trnmed )
     else ! night-time
        ! Initialize outputs with netCDF-writable numbers

        ! DISORT outputs
        albmed(:)=0.0
        dfdt(:)=0.0
        rfldir(:)=0.0
        rfldn(:)=0.0
        trnmed(:)=0.0
        u0u(:,:)=0.0
        uavg(:)=0.0
        uu(:,:,:)=0.0

        ! Gaseous absorption
        odxc_spc_CH4(:)=0.0
        odxc_spc_CO2(:)=0.0
        odxc_spc_H2O(:)=0.0
        odxc_spc_H2OH2O(:)=0.0
        odxc_spc_NO2(:)=0.0
        odxc_spc_O2(:)=0.0
        odxc_spc_O2N2(:)=0.0
        odxc_spc_O2O2(:)=0.0
        odxc_spc_O3(:)=0.0
        odxc_spc_OH(:)=0.0
        ! Gaseous scattering
        odxc_spc_Ray(:)=0.0
        ! Particle absorption
        odac_spc_aer(:)=0.0
        odac_spc_bga(:)=0.0
        odac_spc_ice(:)=0.0
        odac_spc_lqd(:)=0.0
        odac_spc_mpr(:)=0.0
        odac_spc_snw(:)=0.0
        ! Particle extinction
        odxc_spc_aer(:)=0.0
        odxc_spc_bga(:)=0.0
        odxc_spc_ice(:)=0.0
        odxc_spc_lqd(:)=0.0
        odxc_spc_mpr(:)=0.0
        odxc_spc_snw(:)=0.0

        ! Diagnostic two stream properties
        trn_ttl_drc_spc_snw(:,:)=0.0
        rfl_ddm_spc_snw(:,:)=0.0
        rfl_drc_spc_snw(:,:)=0.0
        rfl_dff_spc_snw(:,:)=0.0
     endif ! night-time

     ! call t_stopf('disort')
     
     ! Store spectral fluxes returned from DISORT
     ! DISORT uses a funny indexing system: 
     ! Returned flux arrays, which are on layer interfaces, 
     ! are always indexed starting from 1 at TOA even though same levels 
     ! start with index 0 on input to DISORT, e.g., temper().
     ! Know what you are doing before changing these indices
     do lev_idx=1,lev_nbr
        ntn_spc_mean(bnd_idx,lev_idx)= &
             0.5*(uavg(lev_idx)+uavg(lev_idx+1))/ &
             wvl_dlt(bnd_idx)
     enddo                  ! end loop over lev
     do lev_idx=1,levp_nbr
        flx_spc_dwn_drc(bnd_idx,lev_idx)=rfldir(lev_idx)/ &
             wvl_dlt(bnd_idx)
        flx_spc_dwn_dff(bnd_idx,lev_idx)=rfldn(lev_idx)/ &
             wvl_dlt(bnd_idx)
        flx_spc_dwn(bnd_idx,lev_idx)= &
             flx_spc_dwn_drc(bnd_idx,lev_idx)+ &
             flx_spc_dwn_dff(bnd_idx,lev_idx)
        flx_spc_upw(bnd_idx,lev_idx)=flup(lev_idx)/ &
             wvl_dlt(bnd_idx)
        flx_spc_net(bnd_idx,lev_idx)= &
             flx_spc_dwn(bnd_idx,lev_idx)-flx_spc_upw(bnd_idx,lev_idx)
     enddo                  ! end loop over lev
     
     ! In DISORT1 u0u returned azimuthally averaged intensities at user polar angles
     ! DISORT2 does not by default return any azimuthally averaged intensities
     ! uu returns exact intensities at each user azimuthal angle at each user polar angle.
     ! User polar angles are ordered by increasing cosine of polar angle, 
     ! i.e., from downwelling radiance to upwelling radiance. 
     
     ! NB: Possible bug in current DISORT() documentation
     ! Intensities are returned in arrays in order of descending cosine of polar angle
     ! ndr = nadir  = travelling towards nadir
     ! zen = zenith = travelling towards zenith
     ! 20160515: Recover azimuthally averaged intensites from DISORT2 quantities
     azi_nbr_wvl_dlt=azi_nbr*wvl_dlt(bnd_idx)
     do lev_idx=1,levp_nbr
        do plr_idx=1,plr_nbr
           do azi_idx=1,azi_nbr
              lmn_bb_aa(plr_idx,lev_idx)=lmn_bb_aa(plr_idx,lev_idx)+ &
                   uu(plr_idx,lev_idx,azi_idx)*bnd_wgt_lmn(bnd_idx)
              ntn_bb_aa(plr_idx,lev_idx)=ntn_bb_aa(plr_idx,lev_idx)+ &
                   uu(plr_idx,lev_idx,azi_idx)
           enddo            ! end loop over azi
           ! Normalize _aa_ intensities by azi_nbr
           ! Do NOT normalize _bb_ intensities by wvl_dlt
           lmn_bb_aa(plr_idx,lev_idx)=lumens_per_Watt_555nm* &
                lmn_bb_aa(plr_idx,lev_idx)/azi_nbr
           ntn_bb_aa(plr_idx,lev_idx)= &
                ntn_bb_aa(plr_idx,lev_idx)/azi_nbr
        enddo               ! end loop over plr
        
        do azi_idx=1,azi_nbr
           lmn_spc_aa_ndr(bnd_idx,lev_idx)= &
                lmn_spc_aa_ndr(bnd_idx,lev_idx)+ &
                uu(plr_ndr,lev_idx,azi_idx)*bnd_wgt_lmn(bnd_idx)
           ntn_spc_aa_ndr(bnd_idx,lev_idx)= &
                ntn_spc_aa_ndr(bnd_idx,lev_idx)+uu(plr_ndr,lev_idx,azi_idx)
           ntn_spc_aa_zen(bnd_idx,lev_idx)= &
                ntn_spc_aa_zen(bnd_idx,lev_idx)+ &
                uu(plr_zen,lev_idx,azi_idx)
        enddo            ! end loop over azi
        ! Normalize _spc_aa_ intensities by wvl_dlt and azi_nbr
        lmn_spc_aa_ndr(bnd_idx,lev_idx)=lumens_per_Watt_555nm* &
             lmn_spc_aa_ndr(bnd_idx,lev_idx)/azi_nbr_wvl_dlt
        ntn_spc_aa_ndr(bnd_idx,lev_idx)= &
             ntn_spc_aa_ndr(bnd_idx,lev_idx)/azi_nbr_wvl_dlt
        ntn_spc_aa_zen(bnd_idx,lev_idx)= &
             ntn_spc_aa_zen(bnd_idx,lev_idx)/azi_nbr_wvl_dlt
     enddo                  ! end loop over levp
     do plr_idx=1,plr_nbr
        do azi_idx=1,azi_nbr
           lmn_spc_aa_sfc(plr_idx,bnd_idx)= &
                lmn_spc_aa_sfc(plr_idx,bnd_idx)+ &
                uu(plr_idx,levp_sfc,azi_idx)*bnd_wgt_lmn(bnd_idx)
           ntn_spc_aa_sfc(plr_idx,bnd_idx)= &
                ntn_spc_aa_sfc(plr_idx,bnd_idx)+ &
                uu(plr_idx,levp_sfc,azi_idx)
        enddo            ! end loop over azi
        lmn_spc_aa_sfc(plr_idx,bnd_idx)=lumens_per_Watt_555nm* &
             lmn_spc_aa_sfc(plr_idx,bnd_idx)/azi_nbr_wvl_dlt
        ntn_spc_aa_sfc(plr_idx,bnd_idx)= &
          ntn_spc_aa_sfc(plr_idx,bnd_idx)/azi_nbr_wvl_dlt
     enddo                  ! end loop over plr
     lmn_spc_aa_ndr_TOA(bnd_idx)=lmn_spc_aa_ndr(bnd_idx,levp_TOA)
     lmn_spc_aa_ndr_sfc(bnd_idx)=lmn_spc_aa_sfc(plr_ndr,bnd_idx)
     ntn_spc_aa_ndr_TOA(bnd_idx)=ntn_spc_aa_ndr(bnd_idx,levp_TOA)
     ntn_spc_aa_ndr_sfc(bnd_idx)=ntn_spc_aa_sfc(plr_ndr,bnd_idx)
     ntn_spc_aa_zen_sfc(bnd_idx)=ntn_spc_aa_sfc(plr_zen,bnd_idx)
     if (.false.) then
        ! 20160515: old DISORT1 method
        ntn_spc_aa_ndr_TOA(bnd_idx)=u0u(plr_ndr,1)/ &
             wvl_dlt(bnd_idx)
        ntn_spc_aa_ndr_sfc(bnd_idx)=u0u(plr_ndr,levp_sfc)/ &
             wvl_dlt(bnd_idx)
        ntn_spc_aa_zen_sfc(bnd_idx)=u0u(plr_zen,levp_sfc)/ &
             wvl_dlt(bnd_idx)

        do lev_idx=1,levp_nbr
           do plr_idx=1,plr_nbr
              ntn_bb_aa(plr_idx,lev_idx)=ntn_bb_aa(plr_idx,lev_idx)+ &
                   u0u(plr_idx,lev_idx)
           enddo               ! end loop over plr
           ntn_spc_aa_ndr(bnd_idx,lev_idx)= &
                u0u(plr_ndr,lev_idx)/wvl_dlt(bnd_idx)
           ntn_spc_aa_zen(bnd_idx,lev_idx)= &
                u0u(plr_zen,lev_idx)/wvl_dlt(bnd_idx)
        enddo                  ! end loop over lev
        do plr_idx=1,plr_nbr
           ntn_spc_aa_sfc(plr_idx,bnd_idx)= &
                u0u(plr_idx,levp_sfc)/wvl_dlt(bnd_idx)
        enddo                  ! end loop over lev
     endif ! endif false
        
     if (bnd_idx==bnd_dbg) then
        do lev_idx=1,levp_nbr
           do plr_idx=1,plr_nbr
              do azi_idx=1,azi_nbr
                 ntn_spc_chn(azi_idx,plr_idx,lev_idx)= &
                      uu(plr_idx,lev_idx,azi_idx)/wvl_dlt(bnd_idx)
              enddo   ! end loop over azi
           enddo      ! end loop over plr
        enddo         ! end loop over lev
     endif            ! end if dbg band
     
     do azi_idx=1,azi_nbr
        do plr_idx=1,plr_nbr
           ntn_spc_TOA(azi_idx,plr_idx,bnd_idx)= &
                uu(plr_idx,1,azi_idx)/wvl_dlt(bnd_idx)
        enddo               ! end loop over plr
     enddo                  ! end loop over azi
     
!     if (sv_ntn) then
!        do azi_idx=1,azi_nbr
!           do plr_idx=1,plr_nbr
!             do lev_idx=1,levp_nbr
!                 ntn_spc_aa(plr_idx,bnd_idx,lev_idx)= &
!                      u0u(plr_idx,lev_idx)/wvl_dlt(bnd_idx)
!              enddo         ! end loop over lev
!           enddo            ! end loop over plr
!        enddo               ! end loop over azi
!     endif                  ! end if saving full intensity arrays
     
999  continue
     
  enddo                     ! end loop over bnd
  
  ! call t_prf(0)
  
  !$omp end do
  !$omp end parallel
  
  ! End section 2: Main computation loop over all bands
  ! Begin section 3: Post-processing
  ! computation of diagnostic arrays from output storage arrays
  ! summation of spectral diagnostic arrays to integrated arrays
  ! reduction of integrated arrays to diagnostic scalars
  ! netCDF output
  
  if (dbg_lvl>dbg_off) write (6,'(/)')
  
  ! If we are storing values at computational angles then
  ! they should have been returned in umu array.
  if (sv_cmp_plr_ngl) then
     do plr_idx=1,plr_nbr
        plr_cos(plr_idx)=umu(plr_idx)
        plr(plr_idx)=acos(plr_cos(plr_idx))
        plr_dgr(plr_idx)=180.0*plr(plr_idx)/pi
     enddo                  ! end loop over plr
  endif                     ! endif looking at computational angles
  
  do bnd_idx=1,bnd_nbr
     trn_spc_atm_CO2(bnd_idx)=exp(-min(odxc_spc_CO2(bnd_idx),25.0))
     trn_spc_atm_H2O(bnd_idx)=exp(-min(odxc_spc_H2O(bnd_idx),25.0))
     trn_spc_atm_H2OH2O(bnd_idx)=exp(-min(odxc_spc_H2OH2O(bnd_idx),25.0))
     trn_spc_atm_ice(bnd_idx)=exp(-min(odxc_spc_ice(bnd_idx),25.0))
     trn_spc_atm_lqd(bnd_idx)=exp(-min(odxc_spc_lqd(bnd_idx),25.0))
     trn_spc_atm_aer(bnd_idx)=exp(-min(odxc_spc_aer(bnd_idx),25.0))
     trn_spc_atm_mpr(bnd_idx)=exp(-min(odxc_spc_mpr(bnd_idx),25.0))
     trn_spc_atm_snw(bnd_idx)=exp(-min(odxc_spc_snw(bnd_idx),25.0))
     trn_spc_atm_bga(bnd_idx)=exp(-min(odxc_spc_bga(bnd_idx),25.0))
     trn_spc_atm_OH(bnd_idx)=exp(-min(odxc_spc_OH(bnd_idx),25.0))
     trn_spc_atm_CH4(bnd_idx)=exp(-min(odxc_spc_CH4(bnd_idx),25.0))
     trn_spc_atm_O2(bnd_idx)=exp(-min(odxc_spc_O2(bnd_idx),25.0))
     trn_spc_atm_O3(bnd_idx)=exp(-min(odxc_spc_O3(bnd_idx),25.0))
     trn_spc_atm_O2O2(bnd_idx)=exp(-min(odxc_spc_O2O2(bnd_idx),25.0))
     trn_spc_atm_O2N2(bnd_idx)=exp(-min(odxc_spc_O2N2(bnd_idx),25.0))
     trn_spc_atm_NO2(bnd_idx)=exp(-min(odxc_spc_NO2(bnd_idx),25.0))
     trn_spc_atm_Ray(bnd_idx)=exp(-min(odxc_spc_Ray(bnd_idx),25.0))
  enddo                     ! end loop over bnd
  do bnd_idx=1,bnd_nbr
     flx_spc_dwn_TOA(bnd_idx)=flx_spc_dwn(bnd_idx,1)
     flx_spc_dwn_sfc(bnd_idx)=flx_spc_dwn(bnd_idx,levp_sfc)
     flx_spc_dwn_dff_sfc(bnd_idx)=flx_spc_dwn_dff(bnd_idx,levp_sfc)
     flx_spc_dwn_drc_sfc(bnd_idx)=flx_spc_dwn_drc(bnd_idx,levp_sfc)
     flx_spc_dwn_snw(bnd_idx)=flx_spc_dwn(bnd_idx,levp_atm_nbr)
     flx_spc_upw_snw(bnd_idx)=flx_spc_upw(bnd_idx,levp_atm_nbr)
     ! Compute absorbed spectral fluxes
     flx_spc_abs_SAS(bnd_idx)=flx_spc_net(bnd_idx,1)
     flx_spc_abs_sfc(bnd_idx)=flx_spc_net(bnd_idx,levp_sfc)
     flx_spc_abs_atm(bnd_idx)= &
          flx_spc_abs_SAS(bnd_idx)-flx_spc_abs_sfc(bnd_idx)
     flx_spc_abs_snw(bnd_idx)= &
          flx_spc_net(bnd_idx,levp_atm_nbr)-flx_spc_abs_sfc(bnd_idx)
     
     do lev_idx=1,lev_nbr
        flx_spc_abs(bnd_idx,lev_idx)= &
             max(0.0, &
             flx_spc_net(bnd_idx,lev_idx)-flx_spc_net(bnd_idx,lev_idx+1))
     enddo                  ! end loop over lev

  enddo                     ! end loop over bnd
  
  ! Compute running sum of atmospheric absorption
  ! NB: flx_abs_atm_rdr holds total atmospheric absorption
  ! occuring red-ward (at longer wavelengths) of a given band. 
  ! This quantity, at certain wavelengths, can be backed out of 
  ! fractional spectral instruments, like Valero's FSBR.
  flx_abs_atm_rdr(1)=flx_spc_abs_atm(1)*wvl_dlt(1)
  do bnd_idx=2,bnd_nbr
     flx_abs_atm_rdr(bnd_idx)=flx_abs_atm_rdr(bnd_idx-1)+ &
          flx_spc_abs_atm(bnd_idx)*wvl_dlt(bnd_idx)
  enddo                     ! end loop over bnd
  
  ! Instrument fluxes are a paradigm useful for any fractional band measurement
  ! Therefore we compute all instrument diagnostics in one location
  ! The computation of instrument diagnostics is perfectly analogous to 
  ! broadband (bb) diagnostics; reuse code by changing bb to nst.
  ! All that need be changed to adapt to a new instrument are the
  ! lower and upper limits of the (rectangular window) bandpass. 
  ! Initialize basic instrument quantities which will be incremented
  if (flt_lmn) then
     ! Treat luminosity computation similarly to instrument bandpass
     do levp_idx=1,levp_nbr
        ilm_dwn(levp_idx)=0.0
        ilm_upw(levp_idx)=0.0
     enddo                     ! end loop over levp
     do bnd_idx=1,bnd_nbr
        do levp_idx=1,levp_nbr
           ilm_dwn(levp_idx)=ilm_dwn(levp_idx)+ &
                lumens_per_Watt_555nm* &
                flx_spc_dwn(bnd_idx,levp_idx)* &
                wvl_dlt(bnd_idx)*bnd_wgt_lmn(bnd_idx)
           ilm_upw(levp_idx)=ilm_upw(levp_idx)+ &
                lumens_per_Watt_555nm* &
                flx_spc_upw(bnd_idx,levp_idx)* &
                wvl_dlt(bnd_idx)*bnd_wgt_lmn(bnd_idx)
        enddo                  ! end loop over levp
     enddo                     ! end loop over bnd
     ilm_dwn_TOA=ilm_dwn(1)
     ilm_dwn_sfc=ilm_dwn(levp_sfc)
  endif ! endif flt_lmn

  do levp_idx=1,levp_nbr
     flx_nst_dwn(levp_idx)=0.0
     flx_nst_upw(levp_idx)=0.0
  enddo                     ! end loop over levp
  ! Initialize default instrument spectral response function
  do bnd_idx=1,bnd_nbr
     bnd_wgt(bnd_idx)=1.0
  enddo                     ! end loop over bnd
  if (flt_nst) then
     ! Instrument filter response affects all spectral bands
     bnd_idx_nst=bnd_nbr_nst
     do bnd_idx=1,bnd_nbr
170     continue
        if (wvl_ctr(bnd_idx)<wvl_min_nst(1)) then
           bnd_idx_nst=1
        else if (wvl_ctr(bnd_idx)>wvl_max_nst(bnd_nbr_nst)) then
           bnd_idx_nst=bnd_nbr_nst
        else if ((wvl_ctr(bnd_idx)>=wvl_min_nst(bnd_idx_nst)).and. &
             (wvl_ctr(bnd_idx)<=wvl_max_nst(bnd_idx_nst))) then
           bnd_idx_nst=bnd_idx_nst
        else
           bnd_idx_nst=bnd_idx_nst-1
           goto 170
        endif
        bnd_wgt(bnd_idx)=nst_SRF(bnd_idx_nst)
     enddo                  ! end loop over bnd
  endif                     ! endif flt_nst
  do bnd_idx=1,bnd_nbr
     do levp_idx=1,levp_nbr
        flx_nst_dwn(levp_idx)=flx_nst_dwn(levp_idx)+ &
             flx_spc_dwn(bnd_idx,levp_idx)*wvl_dlt(bnd_idx)*bnd_wgt(bnd_idx)
        flx_nst_upw(levp_idx)=flx_nst_upw(levp_idx)+ &
             flx_spc_upw(bnd_idx,levp_idx)*wvl_dlt(bnd_idx)*bnd_wgt(bnd_idx)
     enddo                  ! end loop over levp
  enddo                     ! end loop over bnd
  do levp_idx=1,levp_nbr
     flx_nst_net(levp_idx)=flx_nst_dwn(levp_idx)-flx_nst_upw(levp_idx)
  enddo                     ! end loop over levp
  do lev_idx=1,lev_nbr
     flx_nst_abs(lev_idx)=flx_nst_net(lev_idx)-flx_nst_net(lev_idx+1)
  enddo                     ! end loop over lev
  if (flx_nst_dwn(1)>0.0) then 
     abs_nst_SAS=flx_nst_net(1)/flx_nst_dwn(1)
     abs_nst_atm=(flx_nst_net(1)-flx_nst_net(levp_sfc))/flx_nst_dwn(1)
     abs_nst_sfc=flx_nst_net(levp_sfc)/flx_nst_dwn(1)
     rfl_nst_SAS=flx_nst_upw(1)/flx_nst_dwn(1)
     rfl_nst_sfc=flx_nst_upw(levp_sfc)/flx_nst_dwn(levp_sfc)
     trn_nst_atm=flx_nst_dwn(levp_sfc)/flx_nst_dwn(1)
  else
     abs_nst_SAS=0.0
     abs_nst_atm=0.0
     abs_nst_sfc=0.0
     rfl_nst_SAS=0.0
     rfl_nst_sfc=0.0
     trn_nst_atm=0.0
  endif
  flx_nst_abs_atm=flx_nst_net(1)-flx_nst_net(levp_sfc)
  flx_nst_abs_sfc=flx_nst_net(levp_sfc)
  flx_nst_abs_ttl=flx_nst_net(1)
  flx_nst_dwn_TOA=flx_nst_dwn(1)
  flx_nst_dwn_sfc=flx_nst_dwn(levp_sfc)
  ! End instrument computations
  ! Multi-channel instrument computations
  if (.not.mode_std) then
     do chn_idx=1,chn_nbr
        ! do bnd_idx=1,bnd_nbr
        ! bnd_idx_nst=bnd_nbr_nst
        ! 770          continue
        ! if (wvl_ctr(bnd_idx)<wvl_min_nst(1)) then
        ! bnd_idx_nst=1
        ! else if (wvl_ctr(bnd_idx)>wvl_max_nst(bnd_nbr_nst)) then
        ! bnd_idx_nst=bnd_nbr_nst
        ! else if ((wvl_ctr(bnd_idx)>=wvl_min_nst(bnd_idx_nst)).and. &
        !                 (wvl_ctr(bnd_idx)<=wvl_max_nst(bnd_idx_nst))) then
        ! bnd_idx_nst=bnd_idx_nst
        ! else
        ! bnd_idx_nst=bnd_idx_nst-1
        ! goto 770
        ! endif
        ! bnd_chn_wgt(bnd_idx,chn_idx)=chn_SRF(bnd_idx_nst,chn_idx)
        ! write (11, *) 'chn ', chn_idx, 'bnd ', bnd_idx, 'bnd_nst ', bnd_idx_nst, &
        !           bnd_chn_wgt(bnd_idx,chn_idx)
        ! enddo                  ! end loop over bnd
        do azi_idx=1,azi_nbr
           do plr_idx=1,plr_nbr
              ntn_chn_TOA(azi_idx,plr_idx,chn_idx)=0.0
              do bnd_idx=1,bnd_nbr
                 ntn_chn_TOA(azi_idx,plr_idx,chn_idx)= &
                      ntn_chn_TOA(azi_idx,plr_idx,chn_idx)+ &
                      ntn_spc_TOA(azi_idx,plr_idx,bnd_idx)* &
                      wvl_dlt(bnd_idx)*chn_SRF(bnd_idx,chn_idx)
              enddo         ! end loop over bnd
           enddo            ! end loop over plr
        enddo               ! end loop over azi
        flx_chn_upw_TOA(chn_idx)=0.0
        flx_chn_dwn_TOA(chn_idx)=0.0
        do bnd_idx=1,bnd_nbr
           flx_chn_upw_TOA(chn_idx)=flx_chn_upw_TOA(chn_idx)+ &
                flx_spc_upw(bnd_idx,1)* &
                wvl_dlt(bnd_idx)*chn_SRF(bnd_idx,chn_idx)
           flx_chn_dwn_TOA(chn_idx)=flx_chn_dwn_TOA(chn_idx)+ &
                flx_spc_dwn(bnd_idx,1)* &
                wvl_dlt(bnd_idx)*chn_SRF(bnd_idx,chn_idx)
        enddo               ! end loop over bnd
        do azi_idx=1,azi_nbr
           do plr_idx=1,plr_nbr
              if (flx_chn_dwn_TOA(chn_idx) /= 0.0) then
                 rfl_chn_TOA(azi_idx,plr_idx,chn_idx)= &
                      pi*ntn_chn_TOA(azi_idx,plr_idx,chn_idx)/ &
                      flx_chn_dwn_TOA(chn_idx)
              else
                 rfl_chn_TOA(azi_idx,plr_idx,chn_idx)=0.0
              endif
           enddo            ! end loop over plr
        enddo               ! end loop over azi
     enddo                  ! end loop over chn
  endif
  ! End Multi-channel instrument computations
  ! Define system transmittance, reflectance, and absorptance.
  ! These definitions must be made in a conditional clause because
  ! they are all normalized by insolation, which may be zero.
  do bnd_idx=1,bnd_nbr
     if (flx_spc_dwn_TOA(bnd_idx)>0.0) then
        trn_spc_atm_ttl(bnd_idx)= &
             flx_spc_dwn_sfc(bnd_idx)/ &
             flx_spc_dwn_TOA(bnd_idx)
        rfl_spc_SAS(bnd_idx)= &
             flx_spc_upw(bnd_idx,1)/ &
             flx_spc_dwn_TOA(bnd_idx)
        ! Layer absorptance is absorbed flux normalized by flux entering layer
        ! Define absorptance so surface + atmospheric absorptances sum to total SAS absorptance,
        ! i.e., as fraction of insolation absorbed by atmosphere, surface, and SAS, respectively.
        abs_spc_SAS(bnd_idx)= &
             flx_spc_net(bnd_idx,1)/flx_spc_dwn_TOA(bnd_idx)
        abs_spc_sfc(bnd_idx)= &
             flx_spc_net(bnd_idx,levp_sfc)/flx_spc_dwn_TOA(bnd_idx)
     else ! flx_spc_dwn_TOA<=0
        trn_spc_atm_ttl(bnd_idx)=0.0
        rfl_spc_SAS(bnd_idx)=0.0
        abs_spc_SAS(bnd_idx)=0.0
        abs_spc_sfc(bnd_idx)=0.0
     endif ! flx_spc_dwn_TOA<=0
     if (flx_spc_dwn_snw(bnd_idx)>0.0) then
        alb_spc_snw(bnd_idx)= & ! [frc] Snowpack spectral flux reflectance
             flx_spc_upw(bnd_idx,levp_atm_nbr)/ &
             flx_spc_dwn(bnd_idx,levp_atm_nbr)
     else ! flx_spc_dwn_snw<=0
        alb_spc_snw(bnd_idx)=0.0 ! [frc] Snowpack spectral flux reflectance
     endif ! flx_spc_dwn_snw<=0
  enddo                     ! end loop over bnd
  do bnd_idx=1,bnd_nbr
     abs_spc_atm(bnd_idx)=abs_spc_SAS(bnd_idx)-abs_spc_sfc(bnd_idx)
  enddo                     ! end loop over bnd
  
  do bnd_idx=1,bnd_nbr
     ! Compute actinic fluxes
     nrg_pht(bnd_idx)=Planck*speed_of_light/wvl(bnd_idx) ! J pht-1
     flx_spc_pht_dwn_sfc(bnd_idx)=flx_spc_dwn_sfc(bnd_idx)/nrg_pht(bnd_idx) ! pht m-2 s-1 m-1
     do lev_idx=1,lev_nbr
        flx_spc_act=4.0*pi*ntn_spc_mean(bnd_idx,lev_idx) ! W m-2 m-1 sr-1 --> W m-2 m-1
        flx_spc_act_pht=flx_spc_act/nrg_pht(bnd_idx) ! W m-2 m-1 --> pht m-2 s-1 m-1
        j_spc_NO2=          & ! s-1 m-1
             flx_spc_act_pht*abs_xsx_NO2(bnd_idx)*qnt_yld_NO2(bnd_idx)
        j_NO2(lev_idx)=     & ! s-1 m-1 --> s-1 
             j_NO2(lev_idx)+j_spc_NO2*wvl_dlt(bnd_idx)
        ! j_NO2 refers to photolysis of NO2 into O(3P) + NO
     enddo                  ! end loop over lev
     flx_spc_act_pht_TOA(bnd_idx)=4.0*pi*ntn_spc_mean(bnd_idx,lev_TOA)/nrg_pht(bnd_idx)
     flx_spc_act_pht_sfc(bnd_idx)=flx_spc_act_pht
     j_spc_NO2_sfc(bnd_idx)=j_spc_NO2
     
     ! Accumulate spectral fluxes into broadband arrays
     do lev_idx=1,lev_nbr
        ntn_bb_mean(lev_idx)=ntn_bb_mean(lev_idx)+ &
             ntn_spc_mean(bnd_idx,lev_idx)* &
             wvl_dlt(bnd_idx)
     enddo                  ! end loop over lev
     do levp_idx=1,levp_nbr
        flx_bb_dwn_drc(levp_idx)=flx_bb_dwn_drc(levp_idx)+ &
             flx_spc_dwn_drc(bnd_idx,levp_idx)* &
             wvl_dlt(bnd_idx)
        flx_bb_dwn_dff(levp_idx)=flx_bb_dwn_dff(levp_idx)+ &
             flx_spc_dwn_dff(bnd_idx,levp_idx)* &
             wvl_dlt(bnd_idx)
        flx_bb_upw(levp_idx)=flx_bb_upw(levp_idx)+ &
             flx_spc_upw(bnd_idx,levp_idx)* &
             wvl_dlt(bnd_idx)
     enddo                  ! end loop over levp
  enddo                     ! end loop over bnd
  
  ! Process broadband fluxes
  do levp_idx=1,levp_nbr
     flx_bb_dwn(levp_idx)=flx_bb_dwn_drc(levp_idx)+flx_bb_dwn_dff(levp_idx)
     flx_bb_net(levp_idx)=flx_bb_dwn(levp_idx)-flx_bb_upw(levp_idx)
  enddo                     ! end loop over levp
  
  ! Compute scalar diagnostics
  ! NB: Somewhat strange definitions here---
  ! Atmospheric quantities normalized by TOA insolation
  ! Snowpack quantities normalized by snowpack insolation
  abs_bb_SAS=flx_bb_net(1)/max(flx_bb_dwn(1),real_tiny)
  abs_bb_atm=(flx_bb_net(1)-flx_bb_net(levp_sfc))/max(flx_bb_dwn(1),real_tiny)
  abs_bb_sfc=flx_bb_net(levp_sfc)/max(flx_bb_dwn(1),real_tiny)
  abs_bb_snw=flx_bb_net(levp_atm_nbr)/max(flx_bb_dwn(levp_atm_nbr),real_tiny)
  rfl_bb_SAS=flx_bb_upw(1)/max(flx_bb_dwn(1),real_tiny)
  rfl_bb_sfc=flx_bb_upw(levp_sfc)/max(flx_bb_dwn(levp_sfc),real_tiny)
  rfl_bb_snw=flx_bb_upw(levp_atm_nbr)/max(flx_bb_dwn(levp_atm_nbr),real_tiny)
  trn_bb_atm=flx_bb_dwn(levp_sfc)/max(flx_bb_dwn(1),real_tiny)
  trn_bb_snw=flx_bb_dwn(levp_sfc)/max(flx_bb_dwn(levp_atm_nbr),real_tiny)
  flx_bb_abs_ttl=flx_bb_net(1)
  flx_bb_abs_sfc=flx_bb_net(levp_sfc)
  flx_bb_abs_atm=flx_bb_net(1)-flx_bb_net(levp_sfc)
  flx_bb_abs_snw=flx_bb_net(levp_atm_nbr)-flx_bb_net(levp_sfc)
  flx_bb_dwn_TOA=flx_bb_dwn(1)
  flx_bb_upw_TOA=flx_bb_upw(1)
  flx_bb_dwn_sfc=flx_bb_dwn(levp_sfc)
  flx_bb_dwn_dff_sfc=flx_bb_dwn_dff(levp_sfc)
  flx_bb_dwn_drc_sfc=flx_bb_dwn_drc(levp_sfc)
  flx_bb_dwn_snw=flx_bb_dwn(levp_atm_nbr)
  
  do lev_idx=1,lev_nbr
     flx_bb_abs(lev_idx)=flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1)
     htg_rate_bb(lev_idx)= &
          (flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1))*grv(lev_idx)/ &
          (spc_heat_mst_air(lev_idx)*prs_dlt(lev_idx))
!     write (6,'(a,i2,a,es8.1)') ' flurt flx_bb_abs(',lev_idx,') =  ',flx_bb_abs(lev_idx)
  enddo                     ! end loop over lev
  ! write (6,'(i4,a,es8.1)') bnd_idx,' flurt flx_spc_dwn(1690,levp_sfc) =  ',flx_spc_dwn(1690,levp_sfc)
  
  ! Compute spectral diagnostics that depend on total fluxes
  do bnd_idx=1,bnd_nbr
     flx_frc_dwn_sfc(bnd_idx)=flx_spc_dwn_sfc(bnd_idx)*wvl_dlt(bnd_idx)/max(flx_bb_dwn_sfc,real_tiny)
  enddo                     ! end loop over bnd
  flx_frc_dwn_sfc_blr(1)=1.0
  do bnd_idx=2,bnd_nbr
     flx_frc_dwn_sfc_blr(bnd_idx)=flx_frc_dwn_sfc_blr(bnd_idx-1)-flx_frc_dwn_sfc(bnd_idx)
  enddo                     ! end loop over bnd
  flx_frc_dwn_sfc_blr(bnd_nbr)=0.0

  if (mode_std) then
     ! Begin netCDF output routines
#ifdef ENABLE_NETCDF4
     if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_classic ! [enm] Output file format
     if (fl_out_fmt == nf90_format_64bit) then
        nf90_create_mode=nf90_create_mode+nf90_64bit_offset
     else if (fl_out_fmt == nf90_format_netcdf4) then
        nf90_create_mode=nf90_create_mode+nf90_netcdf4
     else if (fl_out_fmt == nf90_format_netcdf4_classic) then
        nf90_create_mode=nf90_create_mode+(nf90_classic_model+nf90_netcdf4)
     end if ! end else fl_out_fmt
#else /* !ENABLE_NETCDF4 */
     if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_classic ! [enm] Output file format
     if(fl_out_fmt == nf90_format_classic) nf90_create_mode=nf90_create_mode+0 ! CEWI
#endif /* !ENABLE_NETCDF4 */
     dfl_lvl=dfl_lvl+0 ! CEWI
     flg_dfl=flg_dfl+0 ! CEWI
     flg_shf=flg_shf+0 ! CEWI
     rcd=nf90_wrp_create(fl_out,nf90_create_mode,nc_id,sbr_nm=sbr_nm)
     
     ! Define dimension IDs
     rcd=nf90_wrp(nf90_def_dim(nc_id,'azi',azi_nbr,azi_dmn_id),sbr_nm//': def_dim azi in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dmn_id),sbr_nm//': def_dim bnd in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'grd',bnd_nbr+1,grd_dmn_id),sbr_nm//': def_dim grd in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'lev',lev_nbr,lev_dmn_id),sbr_nm//': def_dim lev in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'levp',levp_nbr,levp_dmn_id),sbr_nm//': def_dim levp in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'plr',plr_nbr,plr_dmn_id),sbr_nm//': def_dim plr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'tau',tau_nbr,tau_dmn_id),sbr_nm//': def_dim tau in '//__FILE__)
     if (lev_snw_nbr > 0) then
        rcd=nf90_wrp(nf90_def_dim(nc_id,'lev_snw',lev_snw_nbr,lev_snw_dmn_id),sbr_nm//': def_dim lev_snw in '//__FILE__)
        rcd=nf90_wrp(nf90_def_dim(nc_id,'levp_snw',levp_snw_nbr,levp_snw_dmn_id),sbr_nm//': def_dim levp_snw in '//__FILE__)
        dim_bnd_lev_snw=(/bnd_dmn_id,lev_snw_dmn_id/)
        dim_bnd_levp_snw=(/bnd_dmn_id,levp_snw_dmn_id/)
     endif ! lev_snw_nbr==0
     if (flg_mie) then
        rcd=nf90_wrp(nf90_def_dim(nc_id,'mmn',mmn_nbr+1,mmn_dmn_id),sbr_nm//': def_dim mmn in '//__FILE__)
     endif ! !flg_mie

     ! Assemble ID and count vectors for each multidimensional combination of dimensions
     dim_bnd_lev=(/bnd_dmn_id,lev_dmn_id/)
     dim_bnd_levp=(/bnd_dmn_id,levp_dmn_id/)
     dim_plr_levp=(/plr_dmn_id,levp_dmn_id/)
     dim_plr_bnd=(/plr_dmn_id,bnd_dmn_id/)
     dim_plr_bnd_levp=(/plr_dmn_id,bnd_dmn_id,levp_dmn_id/)
     dim_mmn_bnd_lev=(/mmn_dmn_id,bnd_dmn_id,lev_dmn_id/)
     dim_azi_plr_levp=(/azi_dmn_id,plr_dmn_id,levp_dmn_id/)
     dim_azi_plr_bnd_levp=(/azi_dmn_id,plr_dmn_id,bnd_dmn_id,levp_dmn_id/)
     ! Variable definitions
     nf90_r8=nf90_xtype_r8_get() ! [enm] External netCDF type for r8 kind
     
     ! Variable definitions
     if (flg_mie) then
        rcd=nf90_wrp(nf90_def_var(nc_id,'lgn_xpn_cff_Mie_ttl',nf90_float,dim_mmn_bnd_lev,lgn_xpn_cff_Mie_ttl_id), &
             sbr_nm//': dv lgn_xpn_cff_Mie_ttl')
     endif ! !flg_mie
     if (lev_snw_nbr > 0) then
        rcd=nf90_wrp(nf90_def_var(nc_id,'lev_snw',nf90_float,lev_snw_dmn_id,lev_snw_id),sbr_nm//': dv lev_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'levp_snw',nf90_float,levp_snw_dmn_id,levp_snw_id),sbr_nm//': dv levp_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'dns_snw',nf90_float,lev_snw_dmn_id,dns_snw_id),sbr_nm//': dv dns_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_dlt_snw',nf90_float,lev_snw_dmn_id,dpt_dlt_snw_id),sbr_nm//': dv dpt_dlt_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_ntf_snw',nf90_float,levp_snw_dmn_id,dpt_ntf_snw_id),sbr_nm//': dv dpt_ntf_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_snw',nf90_float,lev_snw_dmn_id,dpt_snw_id),sbr_nm//': dv dpt_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'mmr_mpr_snw',nf90_float,lev_snw_dmn_id,mmr_mpr_snw_id),sbr_nm//': dv mmr_mpr_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rds_ffc_snw',nf90_float,lev_snw_dmn_id,rds_ffc_snw_id),sbr_nm//': dv rds_ffc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_snw',nf90_float,lev_snw_dmn_id,tpt_snw_id),sbr_nm//': dv tpt_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ntf_snw',nf90_float,levp_snw_dmn_id,tpt_ntf_snw_id),sbr_nm//': dv tpt_ntf_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'foo_snw',nf90_float,lev_snw_dmn_id,foo_snw_id),sbr_nm//': dv foo_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_ddm_spc_snw',nf90_float,dim_bnd_levp_snw,rfl_ddm_spc_snw_id), &
             sbr_nm//': dv rfl_ddm_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_ddm_spc_snw',nf90_float,dim_bnd_levp_snw,trn_ddm_spc_snw_id), &
             sbr_nm//': dv trn_ddm_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'alb_dff_spc_snw_dea',nf90_float,bnd_dmn_id,alb_dff_spc_snw_dea_id), &
             sbr_nm//': dv alb_dff_spc_snw_dea')
        rcd=nf90_wrp(nf90_def_var(nc_id,'alb_drc_spc_snw_dea',nf90_float,bnd_dmn_id,alb_drc_spc_snw_dea_id), &
             sbr_nm//': dv alb_drc_spc_snw_dea')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_dff_spc_snw',nf90_float,dim_bnd_levp_snw,rfl_dff_spc_snw_id), &
             sbr_nm//': dv rfl_dff_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_dff_spc_snw',nf90_float,dim_bnd_levp_snw,trn_dff_spc_snw_id), &
             sbr_nm//': dv trn_dff_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_drc_spc_snw',nf90_float,dim_bnd_levp_snw,rfl_drc_spc_snw_id), &
             sbr_nm//': dv rfl_drc_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_ttl_drc_spc_snw',nf90_float,dim_bnd_levp_snw,trn_ttl_drc_spc_snw_id), &
             sbr_nm//': dv trn_ttl_drc_spc_snw')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_dff_spc_snw_cnt',nf90_float,dim_bnd_levp_snw,rfl_dff_spc_snw_cnt_id), &
             sbr_nm//': dv rfl_dff_spc_snw_cnt')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_dff_spc_snw_cnt',nf90_float,dim_bnd_levp_snw,trn_dff_spc_snw_cnt_id), &
             sbr_nm//': dv trn_dff_spc_snw_cnt')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_drc_upw_snp',nf90_float,dim_bnd_levp_snw,rfl_drc_upw_snp_id), &
             sbr_nm//': dv rfl_drc_upw_snp')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_dff_upw_snp',nf90_float,dim_bnd_levp_snw,rfl_dff_upw_snp_id), &
             sbr_nm//': dv rfl_dff_upw_snp')
        rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_dff_dwn_snp',nf90_float,dim_bnd_levp_snw,rfl_dff_dwn_snp_id), &
             sbr_nm//': dv rfl_dff_dwn_snp')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_ttl_drc_snp',nf90_float,dim_bnd_levp_snw,trn_ttl_drc_snp_id), &
             sbr_nm//': dv trn_ttl_drc_snp')
        rcd=nf90_wrp(nf90_def_var(nc_id,'trn_drc_drc_snp',nf90_float,dim_bnd_levp_snw,trn_drc_drc_snp_id), &
             sbr_nm//': dv trn_drc_drc_snp')
     endif ! lev_snw_nbr==0
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_spc_SAS',nf90_float,bnd_dmn_id,abs_spc_SAS_id),sbr_nm//': dv abs_spc_SAS')
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_spc_atm',nf90_float,bnd_dmn_id,abs_spc_atm_id),sbr_nm//': dv abs_spc_atm')
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_spc_sfc',nf90_float,bnd_dmn_id,abs_spc_sfc_id),sbr_nm//': dv abs_spc_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'alt_ntf',nf90_float,levp_dmn_id,alt_ntf_id),sbr_nm//': dv alt_ntf')
     rcd=nf90_wrp(nf90_def_var(nc_id,'azi',nf90_float,azi_dmn_id,azi_id),sbr_nm//': dv azi')
     rcd=nf90_wrp(nf90_def_var(nc_id,'azi_dgr',nf90_float,azi_dmn_id,azi_dgr_id),sbr_nm//': dv azi_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'bnd',nf90_float,bnd_dmn_id,bnd_id),sbr_nm//': dv bnd')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_abs',nf90_float,lev_dmn_id,flx_bb_abs_id),sbr_nm//': dv flx_bb_abs')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn',nf90_float,levp_dmn_id,flx_bb_dwn_id),sbr_nm//': dv flx_bb_dwn')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_dff',nf90_float,levp_dmn_id,flx_bb_dwn_dff_id),sbr_nm//': dv flx_bb_dwn_dff')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_drc',nf90_float,levp_dmn_id,flx_bb_dwn_drc_id),sbr_nm//': dv flx_bb_dwn_drc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_net',nf90_float,levp_dmn_id,flx_bb_net_id),sbr_nm//': dv flx_bb_net')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_upw',nf90_float,levp_dmn_id,flx_bb_upw_id),sbr_nm//': dv flx_bb_upw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_abs',nf90_float,lev_dmn_id,flx_nst_abs_id),sbr_nm//': dv flx_nst_abs')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_dwn',nf90_float,levp_dmn_id,flx_nst_dwn_id),sbr_nm//': dv flx_nst_dwn')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_net',nf90_float,levp_dmn_id,flx_nst_net_id),sbr_nm//': dv flx_nst_net')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_upw',nf90_float,levp_dmn_id,flx_nst_upw_id),sbr_nm//': dv flx_nst_upw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dmn_id,flx_slr_frc_id),sbr_nm//': dv flx_slr_frc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_abs',nf90_float,dim_bnd_lev,flx_spc_abs_id),sbr_nm//': dv flx_spc_abs')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn',nf90_float,dim_bnd_levp,flx_spc_dwn_id),sbr_nm//': dv flx_spc_dwn')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_upw',nf90_float,dim_bnd_levp,flx_spc_upw_id),sbr_nm//': dv flx_spc_upw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'htg_rate_bb',nf90_float,lev_dmn_id,htg_rate_bb_id),sbr_nm//': dv htg_rate_bb')
     rcd=nf90_wrp(nf90_def_var(nc_id,'j_NO2',nf90_float,lev_dmn_id,j_NO2_id),sbr_nm//': dv j_NO2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'j_spc_NO2_sfc',nf90_float,bnd_dmn_id,j_spc_NO2_sfc_id),sbr_nm//': dv j_spc_NO2_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lev',nf90_float,lev_dmn_id,lev_id),sbr_nm//': dv lev')
     rcd=nf90_wrp(nf90_def_var(nc_id,'levp',nf90_float,levp_dmn_id,levp_id),sbr_nm//': dv levp')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_SRF',nf90_float,bnd_dmn_id,lmn_SRF_id),sbr_nm//': dv lmn_SRF')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ilm_dwn',nf90_float,levp_dmn_id,ilm_dwn_id),sbr_nm//': dv ilm_dwn')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ilm_upw',nf90_float,levp_dmn_id,ilm_upw_id),sbr_nm//': dv ilm_upw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_bb_aa',nf90_double,dim_plr_levp,lmn_bb_aa_id),sbr_nm//': dv lmn_bb_aa')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_spc_aa_ndr',nf90_double,dim_bnd_levp,lmn_spc_aa_ndr_id),sbr_nm//': dv lmn_spc_aa_ndr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_spc_aa_sfc',nf90_double,dim_plr_bnd,lmn_spc_aa_sfc_id),sbr_nm//': dv lmn_spc_aa_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dmn_id,nrg_pht_id),sbr_nm//': dv nrg_pht')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_bb_aa',nf90_float,dim_plr_levp,ntn_bb_aa_id),sbr_nm//': dv ntn_bb_aa')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_bb_mean',nf90_float,lev_dmn_id,ntn_bb_mean_id),sbr_nm//': dv ntn_bb_mean')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_ndr',nf90_float,dim_bnd_levp,ntn_spc_aa_ndr_id),sbr_nm//': dv ntn_spc_aa_ndr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_sfc',nf90_float,dim_plr_bnd,ntn_spc_aa_sfc_id),sbr_nm//': dv ntn_spc_aa_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_zen',nf90_float,dim_bnd_levp,ntn_spc_aa_zen_id),sbr_nm//': dv ntn_spc_aa_zen')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_chn',nf90_float,dim_azi_plr_levp,ntn_spc_chn_id),sbr_nm//': dv ntn_spc_chn')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_mean',nf90_float,dim_bnd_lev,ntn_spc_mean_id),sbr_nm//': dv ntn_spc_mean')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_aer',nf90_float,bnd_dmn_id,odac_spc_aer_id),sbr_nm//': dv odac_spc_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_mpr',nf90_float,bnd_dmn_id,odac_spc_mpr_id),sbr_nm//': dv odac_spc_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_snw',nf90_float,bnd_dmn_id,odac_spc_snw_id),sbr_nm//': dv odac_spc_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_bga',nf90_float,bnd_dmn_id,odac_spc_bga_id),sbr_nm//': dv odac_spc_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_ice',nf90_float,bnd_dmn_id,odac_spc_ice_id),sbr_nm//': dv odac_spc_ice')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odac_spc_lqd',nf90_float,bnd_dmn_id,odac_spc_lqd_id),sbr_nm//': dv odac_spc_lqd')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odal_obs_aer',nf90_float,lev_dmn_id,odal_obs_aer_id),sbr_nm//': dv odal_obs_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odal_obs_mpr',nf90_float,lev_dmn_id,odal_obs_mpr_id),sbr_nm//': dv odal_obs_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odal_obs_snw',nf90_float,lev_dmn_id,odal_obs_snw_id),sbr_nm//': dv odal_obs_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odal_obs_bga',nf90_float,lev_dmn_id,odal_obs_bga_id),sbr_nm//': dv odal_obs_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odsl_obs_aer',nf90_float,lev_dmn_id,odsl_obs_aer_id),sbr_nm//': dv odsl_obs_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odsl_obs_mpr',nf90_float,lev_dmn_id,odsl_obs_mpr_id),sbr_nm//': dv odsl_obs_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odsl_obs_snw',nf90_float,lev_dmn_id,odsl_obs_snw_id),sbr_nm//': dv odsl_obs_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odsl_obs_bga',nf90_float,lev_dmn_id,odsl_obs_bga_id),sbr_nm//': dv odsl_obs_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_CO2',nf90_float,bnd_dmn_id,odxc_spc_CO2_id),sbr_nm//': dv odxc_spc_CO2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_H2O',nf90_float,bnd_dmn_id,odxc_spc_H2O_id),sbr_nm//': dv odxc_spc_H2O')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_H2OH2O',nf90_float,bnd_dmn_id,odxc_spc_H2OH2O_id),sbr_nm//': dv odxc_spc_H2OH2O')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_NO2',nf90_float,bnd_dmn_id,odxc_spc_NO2_id),sbr_nm//': dv odxc_spc_NO2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_O2',nf90_float,bnd_dmn_id,odxc_spc_O2_id),sbr_nm//': dv odxc_spc_O2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_O2N2',nf90_float,bnd_dmn_id,odxc_spc_O2N2_id),sbr_nm//': dv odxc_spc_O2N2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_O2O2',nf90_float,bnd_dmn_id,odxc_spc_O2O2_id),sbr_nm//': dv odxc_spc_O2O2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_O3',nf90_float,bnd_dmn_id,odxc_spc_O3_id),sbr_nm//': dv odxc_spc_O3')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_OH',nf90_float,bnd_dmn_id,odxc_spc_OH_id),sbr_nm//': dv odxc_spc_OH')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_CH4',nf90_float,bnd_dmn_id,odxc_spc_CH4_id),sbr_nm//': dv odxc_spc_CH4')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_Ray',nf90_float,bnd_dmn_id,odxc_spc_Ray_id),sbr_nm//': dv odxc_spc_Ray')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_aer',nf90_float,bnd_dmn_id,odxc_spc_aer_id),sbr_nm//': dv odxc_spc_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_mpr',nf90_float,bnd_dmn_id,odxc_spc_mpr_id),sbr_nm//': dv odxc_spc_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_snw',nf90_float,bnd_dmn_id,odxc_spc_snw_id),sbr_nm//': dv odxc_spc_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_bga',nf90_float,bnd_dmn_id,odxc_spc_bga_id),sbr_nm//': dv odxc_spc_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_ice',nf90_float,bnd_dmn_id,odxc_spc_ice_id),sbr_nm//': dv odxc_spc_ice')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_lqd',nf90_float,bnd_dmn_id,odxc_spc_lqd_id),sbr_nm//': dv odxc_spc_lqd')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_spc_ttl',nf90_float,bnd_dmn_id,odxc_spc_ttl_id),sbr_nm//': dv odxc_spc_ttl')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_aer',nf90_float,lev_dmn_id,odxl_obs_aer_id),sbr_nm//': dv odxl_obs_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_mpr',nf90_float,lev_dmn_id,odxl_obs_mpr_id),sbr_nm//': dv odxl_obs_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_snw',nf90_float,lev_dmn_id,odxl_obs_snw_id),sbr_nm//': dv odxl_obs_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_bga',nf90_float,lev_dmn_id,odxl_obs_bga_id),sbr_nm//': dv odxl_obs_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr',nf90_float,plr_dmn_id,plr_id),sbr_nm//': dv plr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr_cos',nf90_float,plr_dmn_id,plr_cos_id),sbr_nm//': dv plr_cos')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr_dgr',nf90_float,plr_dmn_id,plr_dgr_id),sbr_nm//': dv plr_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_spc_SAS',nf90_float,bnd_dmn_id,rfl_spc_SAS_id),sbr_nm//': dv rfl_spc_SAS')
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_spc_snw',nf90_float,bnd_dmn_id,alb_spc_snw_id),sbr_nm//': dv alb_spc_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_abs_snw',nf90_float,bnd_dmn_id,flx_spc_abs_snw_id),sbr_nm//': dv flx_spc_abs_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_snw',nf90_float,bnd_dmn_id,flx_spc_dwn_snw_id),sbr_nm//': dv flx_spc_dwn_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_upw_snw',nf90_float,bnd_dmn_id,flx_spc_upw_snw_id),sbr_nm//': dv flx_spc_upw_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_abs_snw',nf90_float,flx_bb_abs_snw_id),sbr_nm//': dv flx_bb_abs_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_snw',nf90_float,flx_bb_dwn_snw_id),sbr_nm//': dv flx_bb_dwn_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_bb_snw',nf90_float,rfl_bb_snw_id),sbr_nm//': dv rfl_bb_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_bb_snw',nf90_float,abs_bb_snw_id),sbr_nm//': dv abs_bb_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_bb_snw',nf90_float,trn_bb_snw_id),sbr_nm//': dv trn_bb_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_spc_sfc',nf90_float,bnd_dmn_id,rfl_spc_sfc_id),sbr_nm//': dv rfl_spc_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'tau',nf90_float,tau_dmn_id,tau_id),sbr_nm//': dv tau')
     rcd=nf90_wrp(nf90_def_var(nc_id,'tau_prs',nf90_float,tau_dmn_id,tau_prs_id),sbr_nm//': dv tau_prs')
     rcd=nf90_wrp(nf90_def_var(nc_id,'tpt',nf90_float,lev_dmn_id,tpt_id),sbr_nm//': dv tpt')
     rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ntf',nf90_float,levp_dmn_id,tpt_ntf_id),sbr_nm//': dv tpt_ntf')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_CO2',nf90_float,bnd_dmn_id,trn_spc_atm_CO2_id),sbr_nm//': dv trn_spc_atm_CO2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_H2O',nf90_float,bnd_dmn_id,trn_spc_atm_H2O_id),sbr_nm//': dv trn_spc_atm_H2O')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_NO2',nf90_float,bnd_dmn_id,trn_spc_atm_NO2_id),sbr_nm//': dv trn_spc_atm_NO2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_O2',nf90_float,bnd_dmn_id,trn_spc_atm_O2_id),sbr_nm//': dv trn_spc_atm_O2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_O3',nf90_float,bnd_dmn_id,trn_spc_atm_O3_id),sbr_nm//': dv trn_spc_atm_O3')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_OH',nf90_float,bnd_dmn_id,trn_spc_atm_OH_id),sbr_nm//': dv trn_spc_atm_OH')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_CH4',nf90_float,bnd_dmn_id,trn_spc_atm_CH4_id),sbr_nm//': dv trn_spc_atm_CH4')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_Ray',nf90_float,bnd_dmn_id,trn_spc_atm_Ray_id),sbr_nm//': dv trn_spc_atm_Ray')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_aer',nf90_float,bnd_dmn_id,trn_spc_atm_aer_id),sbr_nm//': dv trn_spc_atm_aer')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_mpr',nf90_float,bnd_dmn_id,trn_spc_atm_mpr_id),sbr_nm//': dv trn_spc_atm_mpr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_snw',nf90_float,bnd_dmn_id,trn_spc_atm_snw_id),sbr_nm//': dv trn_spc_atm_snw')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_bga',nf90_float,bnd_dmn_id,trn_spc_atm_bga_id),sbr_nm//': dv trn_spc_atm_bga')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_ice',nf90_float,bnd_dmn_id,trn_spc_atm_ice_id),sbr_nm//': dv trn_spc_atm_ice')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_lqd',nf90_float,bnd_dmn_id,trn_spc_atm_lqd_id),sbr_nm//': dv trn_spc_atm_lqd')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_ttl',nf90_float,bnd_dmn_id,trn_spc_atm_ttl_id),sbr_nm//': dv trn_spc_atm_ttl')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl',nf90_float,bnd_dmn_id,wvl_id),sbr_nm//': dv wvl')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dmn_id,wvl_ctr_id),sbr_nm//': dv wvl_ctr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dmn_id,wvl_dlt_id),sbr_nm//': dv wvl_dlt')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dmn_id,wvl_grd_id),sbr_nm//': dv wvl_grd')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dmn_id,wvl_max_id),sbr_nm//': dv wvl_max')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dmn_id,wvl_min_id),sbr_nm//': dv wvl_min')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_ctr',nf90_float,bnd_dmn_id,wvn_ctr_id),sbr_nm//': dv wvn_ctr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_dlt',nf90_float,bnd_dmn_id,wvn_dlt_id),sbr_nm//': dv wvn_dlt')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_max',nf90_float,bnd_dmn_id,wvn_max_id),sbr_nm//': dv wvn_max')
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_min',nf90_float,bnd_dmn_id,wvn_min_id),sbr_nm//': dv wvn_min')
     rcd=nf90_wrp(nf90_def_var(nc_id,'alt',nf90_float,lev_dmn_id,alt_id),sbr_nm//': dv alt')
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_bb_SAS',nf90_float,abs_bb_SAS_id),sbr_nm//': dv abs_bb_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_bb_atm',nf90_float,abs_bb_atm_id),sbr_nm//': dv abs_bb_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_bb_sfc',nf90_float,abs_bb_sfc_id),sbr_nm//': dv abs_bb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_nst_SAS',nf90_float,abs_nst_SAS_id),sbr_nm//': dv abs_nst_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_nst_atm',nf90_float,abs_nst_atm_id),sbr_nm//': dv abs_nst_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'abs_nst_sfc',nf90_float,abs_nst_sfc_id),sbr_nm//': dv abs_nst_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc',nf90_float,alb_sfc_id),sbr_nm//': dv alb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_btm',nf90_float,alt_cld_btm_id),sbr_nm//': dv alt_cld_btm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_thick',nf90_float,alt_cld_thick_id),sbr_nm//': dv alt_cld_thick in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_abs_atm',nf90_float,flx_bb_abs_atm_id),sbr_nm//': dv flx_bb_abs_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_abs_sfc',nf90_float,flx_bb_abs_sfc_id),sbr_nm//': dv flx_bb_abs_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_abs_ttl',nf90_float,flx_bb_abs_ttl_id),sbr_nm//': dv flx_bb_abs_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_TOA',nf90_float,flx_bb_dwn_TOA_id),sbr_nm//': dv flx_bb_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_upw_TOA',nf90_float,flx_bb_upw_TOA_id),sbr_nm//': dv flx_bb_upw_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'frc_ice_ttl',nf90_float,frc_ice_ttl_id),sbr_nm//': dv frc_ice_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'lat_dgr',nf90_double,lat_dgr_id),sbr_nm//': dv lat_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'lcl_time_hr',nf90_double,lcl_time_hr_id),sbr_nm//': dv lcl_time_hr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'lcl_yr_day',nf90_double,lcl_yr_day_id),sbr_nm//': dv lcl_yr_day in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CWP',nf90_float,mpc_CWP_id),sbr_nm//': dv mpc_CWP in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_aer',nf90_float,odxc_obs_aer_id),sbr_nm//': dv odxc_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_mpr',nf90_float,odxc_obs_mpr_id),sbr_nm//': dv odxc_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_snw',nf90_float,odxc_obs_snw_id),sbr_nm//': dv odxc_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_bga',nf90_float,odxc_obs_bga_id),sbr_nm//': dv odxc_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_bb_SAS',nf90_float,rfl_bb_SAS_id),sbr_nm//': dv rfl_bb_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_bb_sfc',nf90_float,rfl_bb_sfc_id),sbr_nm//': dv rfl_bb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_nst_SAS',nf90_float,rfl_nst_SAS_id),sbr_nm//': dv rfl_nst_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_nst_sfc',nf90_float,rfl_nst_sfc_id),sbr_nm//': dv rfl_nst_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_bb_atm',nf90_float,trn_bb_atm_id),sbr_nm//': dv trn_bb_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_nst_atm',nf90_float,trn_nst_atm_id),sbr_nm//': dv trn_nst_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_aer',nf90_float,wvl_obs_aer_id),sbr_nm//': dv wvl_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_mpr',nf90_float,wvl_obs_mpr_id),sbr_nm//': dv wvl_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_snw',nf90_float,wvl_obs_snw_id),sbr_nm//': dv wvl_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_bga',nf90_float,wvl_obs_bga_id),sbr_nm//': dv wvl_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_abs_SAS',nf90_float,bnd_dmn_id,flx_spc_abs_SAS_id),sbr_nm//': dv flx_spc_abs_SAS')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_abs_atm',nf90_float,bnd_dmn_id,flx_spc_abs_atm_id),sbr_nm//': dv flx_spc_abs_atm')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_abs_sfc',nf90_float,bnd_dmn_id,flx_spc_abs_sfc_id),sbr_nm//': dv flx_spc_abs_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dmn_id,flx_spc_dwn_TOA_id),sbr_nm//': dv flx_spc_dwn_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_sfc',nf90_float,bnd_dmn_id,flx_spc_dwn_sfc_id),sbr_nm//': dv flx_spc_dwn_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_frc_dwn_sfc',nf90_float,bnd_dmn_id,flx_frc_dwn_sfc_id),sbr_nm//': dv flx_frc_dwn_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_O2N2',nf90_float,bnd_dmn_id,trn_spc_atm_O2N2_id),sbr_nm//': dv trn_spc_atm_O2N2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_O2O2',nf90_float,bnd_dmn_id,trn_spc_atm_O2O2_id),sbr_nm//': dv trn_spc_atm_O2O2')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_dff',nf90_float,dim_bnd_levp,flx_spc_dwn_dff_id),sbr_nm//': dv flx_spc_dwn_dff')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_drc',nf90_float,dim_bnd_levp,flx_spc_dwn_drc_id),sbr_nm//': dv flx_spc_dwn_drc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_abs_atm_rdr',nf90_float,bnd_dmn_id,flx_abs_atm_rdr_id),sbr_nm//': dv flx_abs_atm_rdr')
     ! Wrap
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_NIR_dff',nf90_float,alb_sfc_NIR_dff_id), &
          sbr_nm//': dv alb_sfc_NIR_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_NIR_drc',nf90_float,alb_sfc_NIR_drc_id), &
          sbr_nm//': dv alb_sfc_NIR_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_vsb_dff',nf90_float,alb_sfc_vsb_dff_id), &
          sbr_nm//': dv alb_sfc_vsb_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_vsb_drc',nf90_float,alb_sfc_vsb_drc_id), &
          sbr_nm//': dv alb_sfc_vsb_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'trn_spc_atm_H2OH2O',nf90_float,bnd_dmn_id,trn_spc_atm_H2OH2O_id), &
          sbr_nm//': dv trn_spc_atm_H2OH2O')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_act_pht_TOA',nf90_float,bnd_dmn_id,flx_spc_act_pht_TOA_id), &
          sbr_nm//': dv flx_spc_act_pht_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_act_pht_sfc',nf90_float,bnd_dmn_id,flx_spc_act_pht_sfc_id), &
          sbr_nm//': dv flx_spc_act_pht_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_dff_sfc',nf90_float,bnd_dmn_id,flx_spc_dwn_dff_sfc_id), &
          sbr_nm//': dv flx_spc_dwn_dff_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_drc_sfc',nf90_float,bnd_dmn_id,flx_spc_dwn_drc_sfc_id), &
          sbr_nm//': dv flx_spc_dwn_drc_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_frc_dwn_sfc_blr',nf90_float,bnd_dmn_id,flx_frc_dwn_sfc_blr_id), &
          sbr_nm//': dv flx_frc_dwn_sfc_blr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_pht_dwn_sfc',nf90_float,bnd_dmn_id,flx_spc_pht_dwn_sfc_id), &
          sbr_nm//': dv flx_spc_pht_dwn_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_spc_aa_ndr_TOA',nf90_float,bnd_dmn_id,lmn_spc_aa_ndr_TOA_id), &
          sbr_nm//': dv lmn_spc_aa_ndr_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lmn_spc_aa_ndr_sfc',nf90_float,bnd_dmn_id,lmn_spc_aa_ndr_sfc_id), &
          sbr_nm//': dv lmn_spc_aa_ndr_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_ndr_TOA',nf90_float,bnd_dmn_id,ntn_spc_aa_ndr_TOA_id), &
          sbr_nm//': dv ntn_spc_aa_ndr_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_ndr_sfc',nf90_float,bnd_dmn_id,ntn_spc_aa_ndr_sfc_id), &
          sbr_nm//': dv ntn_spc_aa_ndr_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa_zen_sfc',nf90_float,bnd_dmn_id,ntn_spc_aa_zen_sfc_id), &
          sbr_nm//': dv ntn_spc_aa_zen_sfc')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_sfc',nf90_float,flx_bb_dwn_sfc_id), &
          sbr_nm//': dv flx_bb_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_dff_sfc',nf90_float,flx_bb_dwn_dff_sfc_id), &
          sbr_nm//': dv flx_bb_dwn_dff_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bb_dwn_drc_sfc',nf90_float,flx_bb_dwn_drc_sfc_id), &
          sbr_nm//': dv flx_bb_dwn_drc_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_abs_atm',nf90_float,flx_nst_abs_atm_id), &
          sbr_nm//': dv flx_nst_abs_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_abs_sfc',nf90_float,flx_nst_abs_sfc_id), &
          sbr_nm//': dv flx_nst_abs_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_abs_ttl',nf90_float,flx_nst_abs_ttl_id), &
          sbr_nm//': dv flx_nst_abs_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_dwn_TOA',nf90_float,flx_nst_dwn_TOA_id), &
          sbr_nm//': dv flx_nst_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_nst_dwn_sfc',nf90_float,flx_nst_dwn_sfc_id), &
          sbr_nm//': dv flx_nst_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'ilm_dwn_TOA',nf90_float,ilm_dwn_TOA_id), &
          sbr_nm//': dv ilm_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'ilm_dwn_sfc',nf90_float,ilm_dwn_sfc_id), &
          sbr_nm//': dv ilm_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_def_var(nc_id,'slr_zen_ngl_cos',nf90_double,slr_zen_ngl_cos_id), &
          sbr_nm//': dv slr_zen_ngl_cos in '//__FILE__)

     ! Add global attributes
     ! If BFB files are important then do not archive dates and command lines
     if (.not.flg_bfb) then
        rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date', &
             lcl_date_time),sbr_nm//': pa creation_date in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln', & 
             cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm//': pa cmd_ln in '//__FILE__)
     endif ! !flg_bfb
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': pa CVS_Id in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm//': pa prg_ID in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'aer_sng',aer_sng(1:ftn_strlen(aer_sng))),sbr_nm//': pa aer_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'mpr_sng',mpr_sng(1:ftn_strlen(mpr_sng))),sbr_nm//': pa mpr_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'snw_sng',snw_sng(1:ftn_strlen(snw_sng))),sbr_nm//': pa snw_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'azi_sng',azi_sng(1:ftn_strlen(azi_sng))),sbr_nm//': pa azi_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'bga_sng',bga_sng(1:ftn_strlen(bga_sng))),sbr_nm//': pa bga_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'plr_sng',plr_sng(1:ftn_strlen(plr_sng))),sbr_nm//': pa plr_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prf_sng',prf_sng(1:ftn_strlen(prf_sng))),sbr_nm//': pa prf_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'str_sng',str_sng(1:ftn_strlen(str_sng))),sbr_nm//': pa str_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_CO2',stt_CO2(1:ftn_strlen(stt_CO2))),sbr_nm//': pa stt_CO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_H2O',stt_H2O(1:ftn_strlen(stt_H2O))),sbr_nm//': pa stt_H2O in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_NO2',stt_NO2(1:ftn_strlen(stt_NO2))),sbr_nm//': pa stt_NO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_O2',stt_O2(1:ftn_strlen(stt_O2))),sbr_nm//': pa stt_O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_O3',stt_O3(1:ftn_strlen(stt_O3))),sbr_nm//': pa stt_O3 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_OH',stt_OH(1:ftn_strlen(stt_OH))),sbr_nm//': pa stt_OH in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_CH4',stt_CH4(1:ftn_strlen(stt_CH4))),sbr_nm//': pa stt_CH4 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_aer',stt_aer(1:ftn_strlen(stt_aer))),sbr_nm//': pa stt_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_bga',stt_bga(1:ftn_strlen(stt_bga))),sbr_nm//': pa stt_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_ice',stt_ice(1:ftn_strlen(stt_ice))),sbr_nm//': pa stt_ice in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_lqd',stt_lqd(1:ftn_strlen(stt_lqd))),sbr_nm//': pa stt_lqd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_mpr',stt_mpr(1:ftn_strlen(stt_mpr))),sbr_nm//': pa stt_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_snw',stt_snw(1:ftn_strlen(stt_snw))),sbr_nm//': pa stt_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_msm',stt_msm(1:ftn_strlen(stt_msm))),sbr_nm//': pa stt_msm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_mie',stt_mie(1:ftn_strlen(stt_mie))),sbr_nm//': pa stt_mie in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_rfl',stt_rfl(1:ftn_strlen(stt_rfl))),sbr_nm//': pa stt_rfl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_slr',stt_slr(1:ftn_strlen(stt_slr))),sbr_nm//': pa stt_slr in '//__FILE__)
     
     ! Wrap
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_Planck',stt_Planck(1:ftn_strlen(stt_Planck))), &
          sbr_nm//': pa stt_Planck in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_Rayleigh',stt_Rayleigh(1:ftn_strlen(stt_Rayleigh))), &
          sbr_nm//': pa stt_Rayleigh in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_flt_nst',stt_flt_nst(1:ftn_strlen(stt_flt_nst))), &
          sbr_nm//': pa stt_flt_nst in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_cld_sat',stt_cld_sat(1:ftn_strlen(stt_cld_sat))), &
          sbr_nm//': pa stt_cld_sat in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_sat_cld',stt_sat_cld(1:ftn_strlen(stt_sat_cld))), &
          sbr_nm//': pa stt_sat_cld in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_sct_lqd',stt_sct_lqd(1:ftn_strlen(stt_sct_lqd))), &
          sbr_nm//': pa stt_sct_lqd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng))), &
          sbr_nm//': pa src_rfr_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'opt_dep_sng',opt_dep_sng(1:ftn_strlen(opt_dep_sng))), &
          sbr_nm//': pa opt_dep_sng in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_H2OH2O',stt_H2OH2O(1:ftn_strlen(stt_H2OH2O))), &
          sbr_nm//': pa stt_H2OH2O in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_Herzberg',stt_Herzberg(1:ftn_strlen(stt_Herzberg))), &
          sbr_nm//': pa stt_Herzberg in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_O2N2',stt_O2N2(1:ftn_strlen(stt_O2N2))), &
          sbr_nm//': pa stt_O2N2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_O2O2',stt_O2O2(1:ftn_strlen(stt_O2O2))), &
          sbr_nm//': pa stt_O2O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_top_lvl',stt_top_lvl(1:ftn_strlen(stt_top_lvl))), &
          sbr_nm//': pa stt_top_lvl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'stt_vpr_H2O_abs_cld',stt_vpr_H2O_abs_cld(1:ftn_strlen(stt_vpr_H2O_abs_cld))), &
          sbr_nm//': pa stt_vpr_H2O_abs_cld in '//__FILE__)

     ! Add english text descriptions
     if (flg_mie) then
        rcd=nf90_wrp(nf90_put_att(nc_id,lgn_xpn_cff_Mie_ttl_id, &
             'long_name','Phase function Legendre polynomial expansion coefficients'), &
             sbr_nm//': pa long_name in '//__FILE__)
     endif ! !flg_mie
     if (lev_snw_nbr > 0) then
        rcd=nf90_wrp(nf90_put_att(nc_id,lev_snw_id,'long_name','Snow depth'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,levp_snw_id,'long_name','Snow depth interfaces'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dns_snw_id,'long_name','Snow density'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_dlt_snw_id,'long_name','Snow layer thickness'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_ntf_snw_id,'long_name','Snow depth interfaces'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_snw_id,'long_name','Snow depth'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,mmr_mpr_snw_id,'long_name','Mass mixing ratio of impurities in snow'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rds_ffc_snw_id,'long_name','Snow effective radius'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,tpt_snw_id,'long_name','Snow temperature'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_snw_id,'long_name','Snow temperature interfaces'), &
          sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,foo_snw_id,'long_name','foo snow'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_ddm_spc_snw_id,'long_name','Snow spectral flux reflectance (direct+diffuse mean)'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ddm_spc_snw_id,'long_name','Snow spectral flux transmittance (direct+diffuse mean)'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_spc_snw_id,'long_name','Snow spectral flux reflectance, diffuse incidence'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_dff_spc_snw_id,'long_name','Snow spectral flux transmittance, diffuse incidence'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_drc_spc_snw_id,'long_name','Snow spectral flux reflectance, direct incidence'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ttl_drc_spc_snw_id,'long_name','Layer total transmittance to direct beam'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_spc_snw_cnt_id,'long_name', &
             'Layer contribution to snowpack spectral albedo to diffuse radiation'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_dff_spc_snw_cnt_id,'long_name', &
             'Layer contribution to snowpack spectral transmittance to diffuse radiation'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_drc_upw_snp_id,'long_name', &
             'hello'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_upw_snp_id,'long_name', &
             'hello'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_dwn_snp_id,'long_name', &
             'hello'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ttl_drc_snp_id,'long_name', &
             'hello'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_drc_drc_snp_id,'long_name', &
             'hello'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,alb_dff_spc_snw_dea_id,'long_name', &
             'Snowpack spectral flux reflectance to isotropic illumination, delta-Eddington/adding method'), &
             sbr_nm//': pa long_name in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,alb_drc_spc_snw_dea_id,'long_name', &
             'Snowpack spectral flux reflectance to direct beam, delta-Eddington/adding method'), &
             sbr_nm//': pa long_name in '//__FILE__)
     endif ! lev_snw_nbr==0
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_SAS_id,'long_name','Broadband absorptance of surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_atm_id,'long_name','Broadband absorptance of surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_sfc_id,'long_name','Broadband absorptance of atmosphere'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_SAS_id,'long_name','FSBR absorptance of surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_atm_id,'long_name','FSBR absorptance of surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_sfc_id,'long_name','FSBR absorptance of atmosphere'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_SAS_id,'long_name','Spectral absorptance of surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_atm_id,'long_name','Spectral absorptance of atmosphere'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_sfc_id,'long_name','Spectral absorptance of surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_id,'long_name','Broadband Lambertian surface albedo'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_dff_id,'long_name','NIR reflectance large zenith angles'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_drc_id,'long_name','NIR reflectance small zenith angles'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_dff_id,'long_name','Visible reflectance at large zenith angles'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_drc_id,'long_name','Visible reflectance at small zenith angles'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_btm_id,'long_name','Highest interface beneath all clouds in column'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_thick_id,'long_name','Thickness of region containing all clouds'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'long_name','Altitude'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'long_name','Interface altitude'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,azi_dgr_id,'long_name','Azimuthal angle (degrees)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,azi_id,'long_name','Azimuthal angle (radians)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'long_name','Midpoint wavelength'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_abs_atm_rdr_id,'long_name','Flux absorbed in atmosphere at longer wavelengths'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_atm_id,'long_name','Broadband flux absorbed by atmospheric column only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_id,'long_name','Broadband flux absorbed by layer'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_sfc_id,'long_name','Broadband flux absorbed by surface only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_ttl_id,'long_name','Broadband flux absorbed by surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_TOA_id,'long_name','Broadband incoming flux at TOA (total insolation)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_upw_TOA_id,'long_name','Broadband upwelling flux at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_dff_id,'long_name','Diffuse downwelling broadband flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_drc_id,'long_name','Direct downwelling broadband flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_id,'long_name','Total downwelling broadband flux (direct + diffuse)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_sfc_id,'long_name','Broadband downwelling flux at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_dff_sfc_id,'long_name','Broadband downwelling diffuse flux at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_drc_sfc_id,'long_name','Broadband downwelling direct flux at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_net_id,'long_name','Net broadband flux (downwelling - upwelling)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_upw_id,'long_name','Upwelling broadband flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_frc_dwn_sfc_blr_id,'long_name','Fraction of insolation at shorter wavelengths'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_frc_dwn_sfc_id,'long_name','Fraction of insolation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_atm_id,'long_name','FSBR flux absorbed by atmospheric column only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_id,'long_name','FSBR flux absorbed by layer'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_sfc_id,'long_name','FSBR flux absorbed by surface only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_ttl_id,'long_name','FSBR flux absorbed by surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_TOA_id,'long_name','FSBR incoming flux at TOA (total insolation)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_id,'long_name','Total downwelling FSBR flux (direct + diffuse)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_sfc_id,'long_name','FSBR downwelling flux at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_net_id,'long_name','Net FSBR flux (downwelling - upwelling)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_upw_id,'long_name','Upwelling FSBR flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_SAS_id,'long_name','Spectral flux absorbed by surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_atm_id,'long_name','Spectral flux absorbed by atmospheric column only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_id,'long_name','Spectral flux absorbed by layer'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_sfc_id,'long_name','Spectral flux absorbed by surface only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_act_pht_TOA_id,'long_name','Spectral actinic photon flux at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_act_pht_sfc_id,'long_name','Spectral actinic photon flux at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral insolation at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_dff_id,'long_name','Spectral diffuse downwelling flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_drc_id,'long_name','Spectral direct downwelling flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_id,'long_name','Spectral downwelling flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_dff_sfc_id,'long_name','Spectral insolation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_drc_sfc_id,'long_name','Spectral insolation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_sfc_id,'long_name','Spectral insolation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_pht_dwn_sfc_id,'long_name','Spectral photon flux downwelling at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_upw_id,'long_name','Spectral upwelling flux'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_ttl_id,'long_name','Fraction of column condensate that is ice'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,htg_rate_bb_id,'long_name','Broadband heating rate'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,j_NO2_id,'long_name','Photolysis rate for NO2 + hv --> O(3P) + NO'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,j_spc_NO2_sfc_id,'long_name','Spectral photolysis rate at sfc for NO2+hv --> O(3P)+NO'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'long_name','Latitude (degrees)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lcl_time_hr_id,'long_name','Local day hour'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lcl_yr_day_id,'long_name','Day of year in local time'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'long_name','Layer pressure'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'long_name','Interface pressure'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_TOA_id,'long_name','Incoming illuminance at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_id,'long_name','Total downwelling illuminance (direct + diffuse)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_sfc_id,'long_name','Downwelling illuminance at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_upw_id,'long_name','Upwelling illuminance'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_SRF_id,'long_name','Luminosity spectral response function'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_id,'long_name','Spectral luminance of nadir radiation'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_bb_aa_id,'long_name','Broadband azimuthally averaged luminance'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_TOA_id,'long_name','Spectral luminance of nadir radiation at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_sfc_id,'long_name','Spectral luminance of nadir radiation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_sfc_id,'long_name','Spectral luminance of radiation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CWP_id,'long_name','Total column Condensed Water Path'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_bb_aa_id,'long_name','Broadband azimuthally averaged intensity'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_bb_mean_id,'long_name','Broadband mean intensity'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_id,'long_name','Spectral intensity of nadir radiation'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_TOA_id,'long_name','Spectral intensity of nadir radiation at TOA'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_sfc_id,'long_name','Spectral intensity of nadir radiation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_sfc_id,'long_name','Spectral intensity of radiation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_zen_id,'long_name','Spectral intensity of zenith radiation'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_zen_sfc_id,'long_name','Spectral intensity of zenith radiation at surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_chn_id,'long_name','Full spectral intensity of particular band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_mean_id,'long_name','Spectral mean intensity'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_aer_id,'long_name','Aerosol absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_mpr_id,'long_name','Impurity absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_snw_id,'long_name','Snow absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_bga_id,'long_name','Background aerosol absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_ice_id,'long_name','Liquid water absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_lqd_id,'long_name','Ice water absorption optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_aer_id,'long_name','Layer aerosol absorption optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_mpr_id,'long_name','Layer impurity absorption optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_snw_id,'long_name','Layer snow absorption optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_bga_id,'long_name','Layer background aerosol absorption optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_aer_id,'long_name','Layer aerosol scattering optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_mpr_id,'long_name','Layer impurity scattering optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_snw_id,'long_name','Layer snow scattering optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_bga_id,'long_name','Layer background aerosol scattering optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_aer_id,'long_name','Column aerosol extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_mpr_id,'long_name','Column impurity extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_snw_id,'long_name','Column snow extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_bga_id,'long_name','Column background aerosol extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_CO2_id,'long_name','CO2 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_H2OH2O_id,'long_name','H2O dimer optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_H2O_id,'long_name','H2O optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_NO2_id,'long_name','NO2 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2N2_id,'long_name','O2N2 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2O2_id,'long_name','O2O2 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2_id,'long_name','O2 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O3_id,'long_name','O3 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_OH_id,'long_name','OH optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_CH4_id,'long_name','CH4 optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_Ray_id,'long_name','Rayleigh scattering optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_aer_id,'long_name','Aerosol extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_mpr_id,'long_name','Snow impurity extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_snw_id,'long_name','Snow extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_bga_id,'long_name','Background aerosol extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_ice_id,'long_name','Ice water extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_lqd_id,'long_name','Liquid water extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_ttl_id,'long_name','Total extinction optical depth to surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_aer_id,'long_name','Layer aerosol extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_mpr_id,'long_name','Layer snow impurity extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_snw_id,'long_name','Layer snow extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_bga_id,'long_name','Layer background aerosol extinction optical depth'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_cos_id,'long_name','Cosine polar angle (degrees)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_dgr_id,'long_name','Polar angle (degrees)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_id,'long_name','Polar angle (radians)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_SAS_id,'long_name','Broadband albedo of entire surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_sfc_id,'long_name','Broadband albedo of surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_nst_SAS_id,'long_name','FSBR albedo of entire surface-atmosphere system'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_nst_sfc_id,'long_name','FSBR albedo of surface'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_spc_SAS_id,'long_name','Spectral planetary flux reflectance'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_spc_snw_id,'long_name','Snowpack spectral flux reflectance'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_snw_id,'long_name','Spectral flux absorbed by snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_snw_id,'long_name','Spectral insolation at snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_upw_snw_id,'long_name','Spectral upwelling flux at snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_snw_id,'long_name','Broadband flux absorbed by snowpack only'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_snw_id,'long_name','Broadband downwelling flux at snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_snw_id,'long_name','Broadband albedo of snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_snw_id,'long_name','Broadband absorptance of snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_bb_snw_id,'long_name','Broadband transmittance of snowpack'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_spc_sfc_id,'long_name','Spectral surface flux reflectance'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_cos_id,'long_name','Cosine solar zenith angle'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_id,'long_name','Optical level (optical depth)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_prs_id,'long_name','Optical level (pressure)'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'long_name','Layer Temperature'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_id,'long_name','Interface temperature'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_bb_atm_id,'long_name','Broadband transmittance of atmospheric column'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_nst_atm_id,'long_name','FSBR transmittancea of atmospheric column'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_CO2_id,'long_name','Column transmission due to CO2 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_H2OH2O_id,'long_name','Column transmission due to H2O dimer absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_H2O_id,'long_name','Column transmission due to H2O absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_NO2_id,'long_name','Column transmission due to NO2 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2N2_id,'long_name','Column transmission due to O2-N2 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2O2_id,'long_name','Column transmission due to O2-O2 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2_id,'long_name','Column transmission due to O2 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O3_id,'long_name','Column transmission due to O3 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_OH_id,'long_name','Column transmission due to OH absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_CH4_id,'long_name','Column transmission due to CH4 absorption'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_Ray_id,'long_name','Column transmission due to Rayleigh scattering'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_aer_id,'long_name','Column transmission due to aerosol extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_mpr_id,'long_name','Column transmission due to impurity extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_snw_id,'long_name','Column transmission due to snow extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_bga_id,'long_name','Column transmission due to background aerosol extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_ice_id,'long_name','Column transmission due to ice extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_lqd_id,'long_name','Column transmission due to liquid extinction'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_ttl_id,'long_name','Spectral flux transmission of entire column'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_id,'long_name','Nominal wavelength in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'long_name','Midpoint wavelength in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'long_name','Width of band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'long_name','Maximum wavelength in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'long_name','Minimum wavelength in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_aer_id,'long_name','Wavelength of aerosol optical depth specification'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_mpr_id,'long_name','Wavelength of snow impurity optical depth specification'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_snw_id,'long_name','Wavelength of snow optical depth specification'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_bga_id,'long_name','Wavelength of background aerosol optical depth specification'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'long_name','Midpoint wavenumber in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth in wavenumbers'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'long_name','Maximum wavenumber in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'long_name','Minimum wavenumber in band'), &
          sbr_nm//': pa long_name in '//__FILE__)
     
     ! Add units
     if (flg_mie) then
        rcd=nf90_wrp(nf90_put_att(nc_id,lgn_xpn_cff_Mie_ttl_id,'units','fraction'), &
             sbr_nm//': pa units in '//__FILE__)
     endif ! !flg_mie
     if (lev_snw_nbr > 0) then
        rcd=nf90_wrp(nf90_put_att(nc_id,lev_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,levp_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dns_snw_id,'units','kilogram meter-3'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_dlt_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_ntf_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,dpt_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,mmr_mpr_snw_id,'units','kilogram kilogram-1'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rds_ffc_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,tpt_snw_id,'units','kelvin'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_snw_id,'units','kelvin'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,foo_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_ddm_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ddm_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_dff_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_drc_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ttl_drc_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_spc_snw_cnt_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_dff_spc_snw_cnt_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_drc_upw_snp_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_upw_snp_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,rfl_dff_dwn_snp_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_ttl_drc_snp_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,trn_drc_drc_snp_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,alb_dff_spc_snw_dea_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
        rcd=nf90_wrp(nf90_put_att(nc_id,alb_drc_spc_snw_dea_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     endif ! lev_snw_nbr==0
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_atm_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_atm_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_nst_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_atm_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_spc_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_dff_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_drc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_dff_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_drc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_btm_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_thick_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,azi_dgr_id,'units','degree'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,azi_id,'units','radian'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_abs_atm_rdr_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_atm_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_ttl_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_TOA_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_upw_TOA_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_dff_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_drc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_dff_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_drc_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_net_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_upw_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_frc_dwn_sfc_blr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_frc_dwn_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_atm_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_abs_ttl_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_TOA_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_dwn_sfc_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_net_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_nst_upw_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_SAS_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_atm_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_sfc_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_dff_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_drc_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_dff_sfc_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_drc_sfc_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_sfc_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_upw_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_ttl_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,htg_rate_bb_id,'units','kelvin second-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_TOA_id,'units','lumen meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_id,'units','lumen meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_dwn_sfc_id,'units','lumen meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ilm_upw_id,'units','lumen meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,j_NO2_id,'units','second-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,j_spc_NO2_sfc_id,'units','second-1 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'units','degree'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lcl_time_hr_id,'units','hour'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lcl_yr_day_id,'units','day'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'units','pascal'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'units','pascal'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_SRF_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_bb_aa_id,'units','lumen meter-2 sterradian-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CWP_id,'units','kilogram meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_bb_aa_id,'units','watt meter-2 sterradian-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_bb_mean_id,'units','watt meter-2 sterradian-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_ice_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odac_spc_lqd_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odal_obs_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odsl_obs_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_CO2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_H2OH2O_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_H2O_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_NO2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2N2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2O2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_O3_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_OH_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_CH4_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_Ray_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_ice_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_lqd_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxc_spc_ttl_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_cos_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_dgr_id,'units','degree'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,plr_id,'units','radian'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_nst_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_nst_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_spc_SAS_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alb_spc_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_abs_snw_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_snw_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_upw_snw_id,'units','watt meter-2 meter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_abs_snw_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_bb_dwn_snw_id,'units','watt meter-2'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_bb_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,abs_bb_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_bb_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rfl_spc_sfc_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_cos_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_prs_id,'units','pascal'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'units','kelvin'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_id,'units','kelvin'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_bb_atm_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_nst_atm_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_CO2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_H2OH2O_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_H2O_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_NO2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2N2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2O2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O2_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_O3_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_OH_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_CH4_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_Ray_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_aer_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_mpr_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_snw_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_bga_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_ice_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_lqd_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,trn_spc_atm_ttl_id,'units','fraction'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_aer_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_mpr_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_snw_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_bga_id,'units','meter'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1'),sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1'),sbr_nm//': pa units in '//__FILE__)
     ! Wrap
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_act_pht_TOA_id,'units','photon meter-2 second-1 meter-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_act_pht_sfc_id,'units','photon meter-2 second-1 meter-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_pht_dwn_sfc_id,'units','photon meter-2 second-1 meter-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_id,'units','lumen meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_sfc_id,'units','lumen meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_TOA_id,'units','lumen meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lmn_spc_aa_ndr_sfc_id,'units','lumen meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_TOA_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_sfc_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_ndr_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_sfc_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_zen_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_aa_zen_sfc_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_chn_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ntn_spc_mean_id,'units','watt meter-2 meter-1 sterradian-1'), &
          sbr_nm//': pa units in '//__FILE__)

     ! Axis attributes
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'axis','Z'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'axis','Z'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'axis','Z'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'axis','Z'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'positive','up'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'positive','up'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'positive','down'),sbr_nm//': pa long_name in '//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'positive','down'),sbr_nm//': pa long_name in '//__FILE__)

     ! All dimensions, variables, and attributes have been defined
     rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
     
     ! Write data
     if (flg_mie) then
        rcd=nf90_wrp(nf90_put_var(nc_id,lgn_xpn_cff_Mie_ttl_id,lgn_xpn_cff_Mie_ttl),sbr_nm//': pv lgn_xpn_cff_Mie_ttl'//__FILE__)
     endif ! !flg_mie
     if (lev_snw_nbr > 0) then
        
        rcd=nf90_wrp(nf90_put_var(nc_id,lev_snw_id,lev_snw),sbr_nm//': pv lev_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,levp_snw_id,levp_snw),sbr_nm//': pv levp_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,dns_snw_id,dns_snw),sbr_nm//': pv dns_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,dpt_dlt_snw_id,dpt_dlt_snw),sbr_nm//': pv dpt_dlt_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,dpt_ntf_snw_id,dpt_ntf_snw),sbr_nm//': pv dpt_ntf_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,dpt_snw_id,dpt_snw),sbr_nm//': pv dpt_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,mmr_mpr_snw_id,mmr_mpr_snw),sbr_nm//': pv mmr_mpr_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rds_ffc_snw_id,rds_ffc_snw),sbr_nm//': pv rds_ffc_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,tpt_snw_id,tpt_snw),sbr_nm//': pv tpt_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ntf_snw_id,tpt_ntf_snw),sbr_nm//': pv tpt_ntf_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,foo_snw_id,foo_snw),sbr_nm//': pv foo_snw'//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_ddm_spc_snw_id,rfl_ddm_spc_snw),sbr_nm//': pv rfl_ddm_spc_snw in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_ddm_spc_snw_id,trn_ddm_spc_snw),sbr_nm//': pv trn_ddm_spc_snw in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_dff_spc_snw_id,rfl_dff_spc_snw),sbr_nm//': pv rfl_dff_spc_snw in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_dff_spc_snw_id,trn_dff_spc_snw),sbr_nm//': pv trn_dff_spc_snw in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_drc_spc_snw_id,rfl_drc_spc_snw),sbr_nm//': pv rfl_drc_spc_snw in '//__FILE__)
        ! Wrap
        rcd=nf90_wrp(nf90_put_var(nc_id,alb_dff_spc_snw_dea_id,alb_dff_spc_snw_dea), &
             sbr_nm//': pv alb_dff_spc_snw_dea in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,alb_drc_spc_snw_dea_id,alb_drc_spc_snw_dea), &
             sbr_nm//': pv alb_drc_spc_snw_dea in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_dff_spc_snw_cnt_id,rfl_dff_spc_snw_cnt), &
             sbr_nm//': pv rfl_dff_spc_snw_cnt in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_dff_spc_snw_cnt_id,trn_dff_spc_snw_cnt), &
             sbr_nm//': pv trn_dff_spc_snw_cnt in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_drc_upw_snp_id,rfl_drc_upw_snp), &
             sbr_nm//': pv rfl_drc_upw_snp in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_dff_upw_snp_id,rfl_dff_upw_snp), &
             sbr_nm//': pv rfl_dff_upw_snp in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,rfl_dff_dwn_snp_id,rfl_dff_dwn_snp), &
             sbr_nm//': pv rfl_dff_dwn_snp in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_ttl_drc_snp_id,trn_ttl_drc_snp), &
             sbr_nm//': pv trn_ttl_drc_snp in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_drc_drc_snp_id,trn_drc_drc_snp), &
             sbr_nm//': pv trn_drc_drc_snp in '//__FILE__)
        rcd=nf90_wrp(nf90_put_var(nc_id,trn_ttl_drc_spc_snw_id,trn_ttl_drc_spc_snw), &
             sbr_nm//': pv trn_ttl_drc_spc_snw in '//__FILE__)
     endif ! lev_snw_nbr==0
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_bb_SAS_id,abs_bb_SAS),sbr_nm//': pv abs_bb_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_bb_atm_id,abs_bb_atm),sbr_nm//': pv abs_bb_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_bb_sfc_id,abs_bb_sfc),sbr_nm//': pv abs_bb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_nst_SAS_id,abs_nst_SAS),sbr_nm//': pv abs_nst_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_nst_atm_id,abs_nst_atm),sbr_nm//': pv abs_nst_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_nst_sfc_id,abs_nst_sfc),sbr_nm//': pv abs_nst_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_spc_SAS_id,abs_spc_SAS),sbr_nm//': pv abs_spc_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_spc_atm_id,abs_spc_atm),sbr_nm//': pv abs_spc_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_spc_sfc_id,abs_spc_sfc),sbr_nm//': pv abs_spc_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_id,alb_sfc),sbr_nm//': pv alb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_NIR_dff_id,alb_sfc_NIR_dff),sbr_nm//': pv alb_sfc_NIR_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_NIR_drc_id,alb_sfc_NIR_drc),sbr_nm//': pv alb_sfc_NIR_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_vsb_dff_id,alb_sfc_vsb_dff),sbr_nm//': pv alb_sfc_vsb_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_vsb_drc_id,alb_sfc_vsb_drc),sbr_nm//': pv alb_sfc_vsb_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_btm_id,alt_cld_btm),sbr_nm//': pv alt_cld_btm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_thick_id,alt_cld_thick),sbr_nm//': pv alt_cld_thick in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alt_id,alt),sbr_nm//': pv alt in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alt_ntf_id,alt_ntf),sbr_nm//': pv alt_ntf in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,azi_dgr_id,azi_dgr),sbr_nm//': pv azi_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,azi_id,azi),sbr_nm//': pv azi in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,bnd_id,bnd),sbr_nm//': pv bnd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_abs_atm_rdr_id,flx_abs_atm_rdr),sbr_nm//': pv flx_abs_atm_rdr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_abs_atm_id,flx_bb_abs_atm),sbr_nm//': pv flx_bb_abs_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_abs_id,flx_bb_abs),sbr_nm//': pv flx_bb_abs in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_abs_sfc_id,flx_bb_abs_sfc),sbr_nm//': pv flx_bb_abs_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_abs_ttl_id,flx_bb_abs_ttl),sbr_nm//': pv flx_bb_abs_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_TOA_id,flx_bb_dwn_TOA),sbr_nm//': pv flx_bb_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_upw_TOA_id,flx_bb_upw_TOA),sbr_nm//': pv flx_bb_upw_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_dff_id,flx_bb_dwn_dff),sbr_nm//': pv flx_bb_dwn_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_drc_id,flx_bb_dwn_drc),sbr_nm//': pv flx_bb_dwn_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_id,flx_bb_dwn),sbr_nm//': pv flx_bb_dwn in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_sfc_id,flx_bb_dwn_sfc),sbr_nm//': pv flx_bb_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_dff_sfc_id,flx_bb_dwn_dff_sfc),sbr_nm//': pv flx_bb_dwn_dff_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_drc_sfc_id,flx_bb_dwn_drc_sfc),sbr_nm//': pv flx_bb_dwn_drc_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_net_id,flx_bb_net),sbr_nm//': pv flx_bb_net in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_upw_id,flx_bb_upw),sbr_nm//': pv flx_bb_upw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_frc_dwn_sfc_id,flx_frc_dwn_sfc),sbr_nm//': pv flx_frc_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_abs_atm_id,flx_nst_abs_atm),sbr_nm//': pv flx_nst_abs_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_abs_id,flx_nst_abs),sbr_nm//': pv flx_nst_abs in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_abs_sfc_id,flx_nst_abs_sfc),sbr_nm//': pv flx_nst_abs_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_abs_ttl_id,flx_nst_abs_ttl),sbr_nm//': pv flx_nst_abs_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_dwn_TOA_id,flx_nst_dwn_TOA),sbr_nm//': pv flx_nst_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_dwn_id,flx_nst_dwn),sbr_nm//': pv flx_nst_dwn in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_dwn_sfc_id,flx_nst_dwn_sfc),sbr_nm//': pv flx_nst_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_net_id,flx_nst_net),sbr_nm//': pv flx_nst_net in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_nst_upw_id,flx_nst_upw),sbr_nm//': pv flx_nst_upw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc),sbr_nm//': pv flx_slr_frc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_abs_SAS_id,flx_spc_abs_SAS),sbr_nm//': pv flx_spc_abs_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_abs_atm_id,flx_spc_abs_atm),sbr_nm//': pv flx_spc_abs_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_abs_id,flx_spc_abs),sbr_nm//': pv flx_spc_abs in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_abs_sfc_id,flx_spc_abs_sfc),sbr_nm//': pv flx_spc_abs_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA),sbr_nm//': pv flx_spc_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_dff_id,flx_spc_dwn_dff),sbr_nm//': pv flx_spc_dwn_dff in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_drc_id,flx_spc_dwn_drc),sbr_nm//': pv flx_spc_dwn_drc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_id,flx_spc_dwn),sbr_nm//': pv flx_spc_dwn in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_sfc_id,flx_spc_dwn_sfc),sbr_nm//': pv flx_spc_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_upw_id,flx_spc_upw),sbr_nm//': pv flx_spc_upw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,frc_ice_ttl_id,frc_ice_ttl),sbr_nm//': pv frc_ice_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,htg_rate_bb_id,htg_rate_bb),sbr_nm//': pv htg_rate_bb in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ilm_dwn_TOA_id,ilm_dwn_TOA),sbr_nm//': pv ilm_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ilm_dwn_id,ilm_dwn),sbr_nm//': pv ilm_dwn in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ilm_dwn_sfc_id,ilm_dwn_sfc),sbr_nm//': pv ilm_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ilm_upw_id,ilm_upw),sbr_nm//': pv ilm_upw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,j_NO2_id,j_NO2),sbr_nm//': pv j_NO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,j_spc_NO2_sfc_id,j_spc_NO2_sfc),sbr_nm//': pv j_spc_NO2_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lat_dgr_id,lat_dgr),sbr_nm//': pv lat_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lcl_time_hr_id,lcl_time_hr),sbr_nm//': pv lcl_time_hr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lcl_yr_day_id,lcl_yr_day),sbr_nm//': pv lcl_yr_day in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lev_id,lev),sbr_nm//': pv lev in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,levp_id,levp),sbr_nm//': pv levp in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_SRF_id,lmn_SRF),sbr_nm//': pv lmn_SRF in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_bb_aa_id,lmn_bb_aa),sbr_nm//': pv lmn_bb_aa in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_spc_aa_ndr_id,lmn_spc_aa_ndr),sbr_nm//': pv lmn_spc_aa_ndr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_spc_aa_sfc_id,lmn_spc_aa_sfc),sbr_nm//': pv lmn_spc_aa_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CWP_id,mpc_CWP),sbr_nm//': pv mpc_CWP in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,nrg_pht_id,nrg_pht),sbr_nm//': pv nrg_pht in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_bb_aa_id,ntn_bb_aa),sbr_nm//': pv ntn_bb_aa in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_bb_mean_id,ntn_bb_mean),sbr_nm//': pv ntn_bb_mean in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_ndr_id,ntn_spc_aa_ndr),sbr_nm//': pv ntn_spc_aa_ndr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_sfc_id,ntn_spc_aa_sfc),sbr_nm//': pv ntn_spc_aa_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_zen_id,ntn_spc_aa_zen),sbr_nm//': pv ntn_spc_aa_zen in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_chn_id,ntn_spc_chn),sbr_nm//': pv ntn_spc_chn in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_mean_id,ntn_spc_mean),sbr_nm//': pv ntn_spc_mean in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_aer_id,odac_spc_aer),sbr_nm//': pv odac_spc_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_mpr_id,odac_spc_mpr),sbr_nm//': pv odac_spc_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_snw_id,odac_spc_snw),sbr_nm//': pv odac_spc_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_bga_id,odac_spc_bga),sbr_nm//': pv odac_spc_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_ice_id,odac_spc_ice),sbr_nm//': pv odac_spc_ice in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odac_spc_lqd_id,odac_spc_lqd),sbr_nm//': pv odac_spc_lqd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odal_obs_aer_id,odal_obs_aer),sbr_nm//': pv odal_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odal_obs_mpr_id,odal_obs_mpr),sbr_nm//': pv odal_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odal_obs_snw_id,odal_obs_snw),sbr_nm//': pv odal_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odal_obs_bga_id,odal_obs_bga),sbr_nm//': pv odal_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odsl_obs_aer_id,odsl_obs_aer),sbr_nm//': pv odsl_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odsl_obs_mpr_id,odsl_obs_mpr),sbr_nm//': pv odsl_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odsl_obs_snw_id,odsl_obs_snw),sbr_nm//': pv odsl_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odsl_obs_bga_id,odsl_obs_bga),sbr_nm//': pv odsl_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_aer_id,odxc_obs_aer),sbr_nm//': pv odxc_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_mpr_id,odxc_obs_mpr),sbr_nm//': pv odxc_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_snw_id,odxc_obs_snw),sbr_nm//': pv odxc_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_bga_id,odxc_obs_bga),sbr_nm//': pv odxc_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_CO2_id,odxc_spc_CO2),sbr_nm//': pv odxc_spc_CO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_H2OH2O_id,odxc_spc_H2OH2O),sbr_nm//': pv odxc_spc_H2OH2O in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_H2O_id,odxc_spc_H2O),sbr_nm//': pv odxc_spc_H2O in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_NO2_id,odxc_spc_NO2),sbr_nm//': pv odxc_spc_NO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_O2N2_id,odxc_spc_O2N2),sbr_nm//': pv odxc_spc_O2N2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_O2O2_id,odxc_spc_O2O2),sbr_nm//': pv odxc_spc_O2O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_O2_id,odxc_spc_O2),sbr_nm//': pv odxc_spc_O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_O3_id,odxc_spc_O3),sbr_nm//': pv odxc_spc_O3 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_OH_id,odxc_spc_OH),sbr_nm//': pv odxc_spc_OH in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_CH4_id,odxc_spc_CH4),sbr_nm//': pv odxc_spc_CH4 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_Ray_id,odxc_spc_Ray),sbr_nm//': pv odxc_spc_Ray in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_aer_id,odxc_spc_aer),sbr_nm//': pv odxc_spc_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_mpr_id,odxc_spc_mpr),sbr_nm//': pv odxc_spc_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_snw_id,odxc_spc_snw),sbr_nm//': pv odxc_spc_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_bga_id,odxc_spc_bga),sbr_nm//': pv odxc_spc_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_ice_id,odxc_spc_ice),sbr_nm//': pv odxc_spc_ice in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_lqd_id,odxc_spc_lqd),sbr_nm//': pv odxc_spc_lqd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxc_spc_ttl_id,odxc_spc_ttl),sbr_nm//': pv odxc_spc_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_aer_id,odxl_obs_aer),sbr_nm//': pv odxl_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_mpr_id,odxl_obs_mpr),sbr_nm//': pv odxl_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_snw_id,odxl_obs_snw),sbr_nm//': pv odxl_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_bga_id,odxl_obs_bga),sbr_nm//': pv odxl_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_cos_id,plr_cos),sbr_nm//': pv plr_cos in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_dgr_id,plr_dgr),sbr_nm//': pv plr_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_id,plr),sbr_nm//': pv plr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_bb_SAS_id,rfl_bb_SAS),sbr_nm//': pv rfl_bb_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_bb_sfc_id,rfl_bb_sfc),sbr_nm//': pv rfl_bb_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_nst_SAS_id,rfl_nst_SAS),sbr_nm//': pv rfl_nst_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_nst_sfc_id,rfl_nst_sfc),sbr_nm//': pv rfl_nst_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_spc_SAS_id,rfl_spc_SAS),sbr_nm//': pv rfl_spc_SAS in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,alb_spc_snw_id,alb_spc_snw),sbr_nm//': pv alb_spc_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_abs_snw_id,flx_spc_abs_snw),sbr_nm//': pv flx_spc_abs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_snw_id,flx_spc_dwn_snw),sbr_nm//': pv flx_spc_dwn_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_upw_snw_id,flx_spc_upw_snw),sbr_nm//': pv flx_spc_upw_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_abs_snw_id,flx_bb_abs_snw),sbr_nm//': pv flx_bb_abs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_bb_dwn_snw_id,flx_bb_dwn_snw),sbr_nm//': pv flx_bb_dwn_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_bb_snw_id,rfl_bb_snw),sbr_nm//': pv rfl_bb_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,abs_bb_snw_id,abs_bb_snw),sbr_nm//': pv abs_bb_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_bb_snw_id,trn_bb_snw),sbr_nm//': pv trn_bb_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_spc_sfc_id,rfl_spc_sfc),sbr_nm//': pv rfl_spc_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos),sbr_nm//': pv slr_zen_ngl_cos in '//__FILE__)
     if (dbg_lvl == 3) write(6,*) 'fxm gfortran tau = ',tau ! fxm gfortran
     rcd=nf90_wrp(nf90_put_var(nc_id,tau_id,tau),sbr_nm//': pv tau in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tau_prs_id,tau_prs),sbr_nm//': pv tau_prs in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tpt_id,tpt),sbr_nm//': pv tpt in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ntf_id,tpt_ntf),sbr_nm//': pv tpt_ntf in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_bb_atm_id,trn_bb_atm),sbr_nm//': pv trn_bb_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_nst_atm_id,trn_nst_atm),sbr_nm//': pv trn_nst_atm in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_CO2_id,trn_spc_atm_CO2),sbr_nm//': pv trn_spc_atm_CO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_H2O_id,trn_spc_atm_H2O),sbr_nm//': pv trn_spc_atm_H2O in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_NO2_id,trn_spc_atm_NO2),sbr_nm//': pv trn_spc_atm_NO2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_O2N2_id,trn_spc_atm_O2N2),sbr_nm//': pv trn_spc_atm_O2N2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_O2O2_id,trn_spc_atm_O2O2),sbr_nm//': pv trn_spc_atm_O2O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_O2_id,trn_spc_atm_O2),sbr_nm//': pv trn_spc_atm_O2 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_O3_id,trn_spc_atm_O3),sbr_nm//': pv trn_spc_atm_O3 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_OH_id,trn_spc_atm_OH),sbr_nm//': pv trn_spc_atm_OH in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_CH4_id,trn_spc_atm_CH4),sbr_nm//': pv trn_spc_atm_CH4 in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_Ray_id,trn_spc_atm_Ray),sbr_nm//': pv trn_spc_atm_Ray in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_aer_id,trn_spc_atm_aer),sbr_nm//': pv trn_spc_atm_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_mpr_id,trn_spc_atm_mpr),sbr_nm//': pv trn_spc_atm_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_snw_id,trn_spc_atm_snw),sbr_nm//': pv trn_spc_atm_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_bga_id,trn_spc_atm_bga),sbr_nm//': pv trn_spc_atm_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_ice_id,trn_spc_atm_ice),sbr_nm//': pv trn_spc_atm_ice in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_lqd_id,trn_spc_atm_lqd),sbr_nm//': pv trn_spc_atm_lqd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_ttl_id,trn_spc_atm_ttl),sbr_nm//': pv trn_spc_atm_ttl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_id,wvl),sbr_nm//': pv wvl in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr),sbr_nm//': pv wvl_ctr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt),sbr_nm//': pv wvl_dlt in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_grd_id,wvl_grd),sbr_nm//': pv wvl_grd in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_max_id,wvl_max),sbr_nm//': pv wvl_max in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_min_id,wvl_min),sbr_nm//': pv wvl_min in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_aer_id,wvl_obs_aer),sbr_nm//': pv wvl_obs_aer in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_mpr_id,wvl_obs_mpr),sbr_nm//': pv wvl_obs_mpr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_snw_id,wvl_obs_snw),sbr_nm//': pv wvl_obs_snw in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_bga_id,wvl_obs_bga),sbr_nm//': pv wvl_obs_bga in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr),sbr_nm//': pv wvn_ctr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt),sbr_nm//': pv wvn_dlt in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvn_max_id,wvn_max),sbr_nm//': pv wvn_max in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wvn_min_id,wvn_min),sbr_nm//': pv wvn_min in '//__FILE__)
     ! Wrap
     ! fxm g95 test character and line counting
#if 0
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_frc_dwn_sfc_blr_id,flx_frc_dwn_sfc_blr),sbr_nm//': pv flx_frc_dwn_sfc_blr in this big old hanging line')
#endif /* !0 */
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_dff_sfc_id,flx_spc_dwn_dff_sfc), &
          sbr_nm//': pv flx_spc_dwn_dff_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_drc_sfc_id,flx_spc_dwn_drc_sfc), &
          sbr_nm//': pv flx_spc_dwn_drc_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_frc_dwn_sfc_blr_id,flx_frc_dwn_sfc_blr), &
          sbr_nm//': pv flx_frc_dwn_sfc_blr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_act_pht_TOA_id,flx_spc_act_pht_TOA), &
          sbr_nm//': pv flx_spc_act_pht_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_act_pht_sfc_id,flx_spc_act_pht_sfc), &
          sbr_nm//': pv flx_spc_act_pht_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_pht_dwn_sfc_id,flx_spc_pht_dwn_sfc), &
          sbr_nm//': pv flx_spc_pht_dwn_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_spc_aa_ndr_TOA_id,lmn_spc_aa_ndr_TOA), &
          sbr_nm//': pv lmn_spc_aa_ndr_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lmn_spc_aa_ndr_sfc_id,lmn_spc_aa_ndr_sfc), &
          sbr_nm//': pv lmn_spc_aa_ndr_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_ndr_TOA_id,ntn_spc_aa_ndr_TOA), &
          sbr_nm//': pv ntn_spc_aa_ndr_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_ndr_sfc_id,ntn_spc_aa_ndr_sfc), &
          sbr_nm//': pv ntn_spc_aa_ndr_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ntn_spc_aa_zen_sfc_id,ntn_spc_aa_zen_sfc), &
          sbr_nm//': pv ntn_spc_aa_zen_sfc in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,trn_spc_atm_H2OH2O_id,trn_spc_atm_H2OH2O), &
          sbr_nm//': pv trn_spc_atm_H2OH2O in '//__FILE__)
     if (sv_ntn) then
        ! Comment-out for now because this consumes too much disk space
        ! ntn_spc_aa_id,
        ! ntn_spc_aa(plr_nbr_max,bnd_nbr_max,levp_nbr_max),
        ! rcd=rcd+nf90_redef(nc_id)
        ! rcd=nf90_wrp(nf90_def_var(nc_id,'ntn_spc_aa',nf90_float,dim_plr_bnd_levp,ntn_spc_aa_id),sbr_nm//': dv ntn_spc_aa')
        ! rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_id,'long_name','Spectral intensity of radiation')
        ! rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_id,'units','watt meter-2 meter-1 sterradian-1')
        ! rcd=rcd+nf90_enddef(nc_id)
        ! rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_id,ntn_spc_aa)
     endif                  ! end if saving full intensity arrays
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
     
  endif                     ! end if mode_std
  
  if (.not. mode_std) then
     rcd=nf90_wrp_create(fl_out,nf90_clobber,nc_id,sbr_nm=sbr_nm)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'azi',azi_nbr,azi_dmn_id),sbr_nm//': def_dim azi in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'plr',plr_nbr,plr_dmn_id),sbr_nm//': def_dim plr in '//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'chn',chn_nbr,chn_dmn_id),sbr_nm//': def_dim chn in '//__FILE__)
     dim_azi_plr_chn=(/azi_dmn_id,plr_dmn_id,chn_dmn_id/)
     rcd=nf90_wrp(nf90_def_var(nc_id,'slr_zen_ngl_cos',nf90_double,slr_zen_ngl_cos_id),sbr_nm//': dv slr_zen_ngl_cos')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lon_dgr',nf90_double,lon_dgr_id),sbr_nm//': dv lon_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'lat_dgr',nf90_double,lat_dgr_id),sbr_nm//': dv lat_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'azi',nf90_float,azi_dmn_id,azi_id),sbr_nm//': dv azi')
     rcd=nf90_wrp(nf90_def_var(nc_id,'azi_dgr',nf90_float,azi_dmn_id,azi_dgr_id),sbr_nm//': dv azi_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr',nf90_float,plr_dmn_id,plr_id),sbr_nm//': dv plr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr_cos',nf90_float,plr_dmn_id,plr_cos_id),sbr_nm//': dv plr_cos')
     rcd=nf90_wrp(nf90_def_var(nc_id,'plr_dgr',nf90_float,plr_dmn_id,plr_dgr_id),sbr_nm//': dv plr_dgr')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_chn_dwn_TOA',nf90_float,chn_dmn_id,flx_chn_dwn_TOA_id),sbr_nm//': dv flx_chn_dwn_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'flx_chn_upw_TOA',nf90_float,chn_dmn_id,flx_chn_upw_TOA_id),sbr_nm//': dv flx_chn_upw_TOA')
     rcd=nf90_wrp(nf90_def_var(nc_id,'rfl_chn_TOA',nf90_float,dim_azi_plr_chn,rfl_chn_TOA_id),sbr_nm//': dv rfl_chn_TOA')
     rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
     rcd=nf90_wrp(nf90_put_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos),sbr_nm//': pv slr_zen_ngl_cos in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lon_dgr_id,lon_dgr),sbr_nm//': pv lon_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lat_dgr_id,lat_dgr),sbr_nm//': pv lat_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,azi_id,azi),sbr_nm//': pv azi in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,azi_dgr_id,azi_dgr),sbr_nm//': pv azi_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_id,plr),sbr_nm//': pv plr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_cos_id,plr_cos),sbr_nm//': pv plr_cos in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,plr_dgr_id,plr_dgr),sbr_nm//': pv plr_dgr in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_chn_dwn_TOA_id,flx_chn_dwn_TOA),sbr_nm//': pv flx_chn_dwn_TOA in '//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,flx_chn_upw_TOA_id,flx_chn_upw_TOA),sbr_nm//': pv flx_chn_upw_TOA in '//__FILE__)
     ! do azi_idx = 1, azi_nbr
     ! do plr_idx = 1, plr_nbr
     ! do chn_idx = 1, chn_nbr
     ! print *, rfl_chn_TOA(azi_idx,plr_idx,chn_idx)
     ! enddo
     ! enddo
     ! enddo
     rcd=nf90_wrp(nf90_put_var(nc_id,rfl_chn_TOA_id,rfl_chn_TOA),sbr_nm//': pv rfl_chn_TOA in '//__FILE__)
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to')
  endif                     ! end if not mode_std
  
  if (rcd /= nf90_noerr) write (6,'(a,a,i4,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR rcd = ',rcd,' on exit'
  
  if (dbg_lvl==dbg_crr) then
     ! fxm:
     ! ncks -v ntn_spc_chn -F -C -m -H ${DATA}/aca/swnb.nc | m
     ! rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to')
     write (6,'(a,3i2)') 'azi_nbr, plr_nbr, levp_nbr = ',azi_nbr,plr_nbr,levp_nbr
     rcd=nf90_wrp_open(fl_out,nf90_nowrite,nc_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ntn_spc_chn',ntn_spc_chn_id)
     rcd=nf90_wrp(nf90_inquire_variable(nc_id,ntn_spc_chn_id,ndims=int_foo),sbr_nm//': inquire_var ntn_spc_chn')
     write (6,'(a,i2)') 'nf90_inquire_variable(ntn_spc_chn) ndims = ',int_foo
     rcd=nf90_wrp(nf90_inquire_variable(nc_id,ntn_spc_chn_id,dimids=dmn_id_vec_foo),sbr_nm//': inquire_var ntn_spc_chn')
     char_foo='   '
     call ftn_strini(char_foo)
     do lev_idx=1,int_foo
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,dmn_id_vec_foo(lev_idx),name=char_foo),sbr_nm//': inquire_dim char_foo')
        rcd=nf90_wrp(nf90_inquire_dimension(nc_id,dmn_id_vec_foo(lev_idx),len=int_foo),sbr_nm//': inquire_dim int_foo')
        write (6,'(3a,i3,a,i3)') 'nf90_inquire_dimension(', &
             char_foo(1:ftn_strlen(char_foo)),'), len= ',int_foo, &
             'id = ',dmn_id_vec_foo(lev_idx)
     enddo                  ! end loop over lev
     write (6,'(a)') 'Zeroing ntn_spc_chn...'
     call vec_set(ntn_spc_chn,azi_nbr*plr_nbr*levp_nbr,0.0)
     write (6,'(a)') 'Reading ntn_spc_chn...'
     rcd=nf90_wrp(nf90_get_var(nc_id,ntn_spc_chn_id,ntn_spc_chn),'gv ntn_spc_chn')
     write (6,'(a,es10.2)') 'ntn_spc_chn(1,4,2) = ',ntn_spc_chn(1,4,2) ! fxm
     rcd=nf90_wrp_close(nc_id,fl_out,'Read something from')
     goto 1000              ! Goto exit with error status
  endif                     ! endif dbg
  
  if (rcd /= nf90_noerr) write (6,'(a,a,i4,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR rcd = ',rcd,' on exit'
  
1000 continue                  ! Jumping point
  
  ! De-allocate dynamic variables
  ! Array dimensions: azi
  if (allocated(azi)) deallocate(azi,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for azi'
  if (allocated(azi_dgr)) deallocate(azi_dgr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for azi_dgr'
  ! Array dimensions: bnd
  if (allocated(abs_spc_SAS)) deallocate(abs_spc_SAS,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for abs_spc_SAS'
  if (allocated(abs_spc_atm)) deallocate(abs_spc_atm,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for abs_spc_atm'
  if (allocated(abs_spc_sfc)) deallocate(abs_spc_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for abs_spc_sfc'
  if (allocated(bnd)) deallocate(bnd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for bnd'
  if (allocated(flx_abs_atm_rdr)) deallocate(flx_abs_atm_rdr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_abs_atm_rdr'
  if (allocated(flx_frc_dwn_sfc_blr)) deallocate(flx_frc_dwn_sfc_blr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_frc_dwn_sfc_blr'
  if (allocated(flx_frc_dwn_sfc)) deallocate(flx_frc_dwn_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_frc_dwn_sfc'
  if (allocated(flx_slr_frc)) deallocate(flx_slr_frc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_slr_frc'
  if (allocated(flx_spc_abs_SAS)) deallocate(flx_spc_abs_SAS,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_abs_SAS'
  if (allocated(flx_spc_abs_atm)) deallocate(flx_spc_abs_atm,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_abs_atm'
  if (allocated(flx_spc_abs_sfc)) deallocate(flx_spc_abs_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_abs_sfc'
  if (allocated(flx_spc_act_pht_TOA)) deallocate(flx_spc_act_pht_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_act_pht_TOA'
  if (allocated(flx_spc_act_pht_sfc)) deallocate(flx_spc_act_pht_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_act_pht_sfc'
  if (allocated(flx_spc_dwn_TOA)) deallocate(flx_spc_dwn_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_TOA'
  if (allocated(flx_spc_dwn_dff_sfc)) deallocate(flx_spc_dwn_dff_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_dff_sfc'
  if (allocated(flx_spc_dwn_drc_sfc)) deallocate(flx_spc_dwn_drc_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_drc_sfc'
  if (allocated(flx_spc_dwn_sfc)) deallocate(flx_spc_dwn_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_sfc'
  if (allocated(flx_spc_pht_dwn_sfc)) deallocate(flx_spc_pht_dwn_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_pht_dwn_sfc'
  if (allocated(j_spc_NO2_sfc)) deallocate(j_spc_NO2_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for j_spc_NO2_sfc'
  if (allocated(lmn_SRF)) deallocate(lmn_SRF,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lmn_SRF'
  if (allocated(lmn_spc_aa_ndr_TOA)) deallocate(lmn_spc_aa_ndr_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lmn_spc_aa_ndr_TOA'
  if (allocated(nrg_pht)) deallocate(nrg_pht,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for nrg_pht'
  if (allocated(ntn_spc_aa_ndr_TOA)) deallocate(ntn_spc_aa_ndr_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_ndr_TOA'
  if (allocated(ntn_spc_aa_ndr_sfc)) deallocate(ntn_spc_aa_ndr_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_ndr_sfc'
  if (allocated(ntn_spc_aa_zen_sfc)) deallocate(ntn_spc_aa_zen_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_zen_sfc'
  if (allocated(odac_spc_aer)) deallocate(odac_spc_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_aer'
  if (allocated(odac_spc_mpr)) deallocate(odac_spc_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_mpr'
  if (allocated(odac_spc_snw)) deallocate(odac_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_snw'
  if (allocated(odac_spc_bga)) deallocate(odac_spc_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_bga'
  if (allocated(odac_spc_ice)) deallocate(odac_spc_ice,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_ice'
  if (allocated(odac_spc_lqd)) deallocate(odac_spc_lqd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odac_spc_lqd'
  if (allocated(odxc_spc_CO2)) deallocate(odxc_spc_CO2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_CO2'
  if (allocated(odxc_spc_H2O)) deallocate(odxc_spc_H2O,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_H2O'
  if (allocated(odxc_spc_H2OH2O)) deallocate(odxc_spc_H2OH2O,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_H2OH2O'
  if (allocated(odxc_spc_NO2)) deallocate(odxc_spc_NO2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_NO2'
  if (allocated(odxc_spc_O2)) deallocate(odxc_spc_O2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_O2'
  if (allocated(odxc_spc_O2N2)) deallocate(odxc_spc_O2N2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_O2N2'
  if (allocated(odxc_spc_O2O2)) deallocate(odxc_spc_O2O2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_O2O2'
  if (allocated(odxc_spc_O3)) deallocate(odxc_spc_O3,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_O3'
  if (allocated(odxc_spc_OH)) deallocate(odxc_spc_OH,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_OH'
  if (allocated(odxc_spc_CH4)) deallocate(odxc_spc_CH4,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_CH4'
  if (allocated(odxc_spc_Ray)) deallocate(odxc_spc_Ray,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_Ray'
  if (allocated(odxc_spc_aer)) deallocate(odxc_spc_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_aer'
  if (allocated(odxc_spc_mpr)) deallocate(odxc_spc_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_mpr'
  if (allocated(odxc_spc_snw)) deallocate(odxc_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_snw'
  if (allocated(odxc_spc_bga)) deallocate(odxc_spc_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_bga'
  if (allocated(odxc_spc_ice)) deallocate(odxc_spc_ice,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_ice'
  if (allocated(odxc_spc_lqd)) deallocate(odxc_spc_lqd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_lqd'
  if (allocated(odxc_spc_ttl)) deallocate(odxc_spc_ttl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxc_spc_ttl'
  if (allocated(rfl_spc_SAS)) deallocate(rfl_spc_SAS,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_spc_SAS'
  if (allocated(alb_spc_snw)) deallocate(alb_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for alb_spc_snw'
  if (allocated(flx_spc_abs_snw)) deallocate(flx_spc_abs_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_abs_snw'
  if (allocated(flx_spc_dwn_snw)) deallocate(flx_spc_dwn_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_snw'
  if (allocated(flx_spc_upw_snw)) deallocate(flx_spc_upw_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_upw_snw'
  if (allocated(rfl_spc_sfc)) deallocate(rfl_spc_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_spc_sfc'
  if (allocated(trn_spc_atm_CO2)) deallocate(trn_spc_atm_CO2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_CO2'
  if (allocated(trn_spc_atm_H2O)) deallocate(trn_spc_atm_H2O,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_H2O'
  if (allocated(trn_spc_atm_H2OH2O)) deallocate(trn_spc_atm_H2OH2O,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_H2OH2O'
  if (allocated(trn_spc_atm_NO2)) deallocate(trn_spc_atm_NO2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_NO2'
  if (allocated(trn_spc_atm_O2)) deallocate(trn_spc_atm_O2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_O2'
  if (allocated(trn_spc_atm_O2N2)) deallocate(trn_spc_atm_O2N2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_O2N2'
  if (allocated(trn_spc_atm_O2O2)) deallocate(trn_spc_atm_O2O2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_O2O2'
  if (allocated(trn_spc_atm_O3)) deallocate(trn_spc_atm_O3,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_O3'
  if (allocated(trn_spc_atm_OH)) deallocate(trn_spc_atm_OH,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_OH'
  if (allocated(trn_spc_atm_CH4)) deallocate(trn_spc_atm_CH4,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_CH4'
  if (allocated(trn_spc_atm_Ray)) deallocate(trn_spc_atm_Ray,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_Ray'
  if (allocated(trn_spc_atm_aer)) deallocate(trn_spc_atm_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_aer'
  if (allocated(trn_spc_atm_mpr)) deallocate(trn_spc_atm_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_mpr'
  if (allocated(trn_spc_atm_bga)) deallocate(trn_spc_atm_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_bga'
  if (allocated(trn_spc_atm_ice)) deallocate(trn_spc_atm_ice,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_ice'
  if (allocated(trn_spc_atm_lqd)) deallocate(trn_spc_atm_lqd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_lqd'
  if (allocated(trn_spc_atm_ttl)) deallocate(trn_spc_atm_ttl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_spc_atm_ttl'
  if (allocated(wvl)) deallocate(wvl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl'
  if (allocated(wvl_ctr)) deallocate(wvl_ctr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl_ctr'
  if (allocated(wvl_dlt)) deallocate(wvl_dlt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl_dlt'
  if (allocated(wvl_max)) deallocate(wvl_max,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl_max'
  if (allocated(wvl_min)) deallocate(wvl_min,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl_min'
  if (allocated(wvn)) deallocate(wvn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvn'
  if (allocated(wvn_ctr)) deallocate(wvn_ctr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvn_ctr'
  if (allocated(wvn_dlt)) deallocate(wvn_dlt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvn_dlt'
  if (allocated(wvn_max)) deallocate(wvn_max,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvn_max'
  if (allocated(wvn_min)) deallocate(wvn_min,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvn_min'
  ! Array dimensions: chn
  if (allocated(flx_chn_dwn_TOA)) deallocate(flx_chn_dwn_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_chn_dwn_TOA'
  if (allocated(flx_chn_upw_TOA)) deallocate(flx_chn_upw_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_chn_upw_TOA'
  ! Array dimensions: grd
  if (allocated(wvl_grd)) deallocate(wvl_grd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for wvl_grd'
  ! Array dimensions: lev
  if (allocated(RH_lqd)) deallocate(RH_lqd,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for RH_lqd'
  if (allocated(alt)) deallocate(alt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for alt'
  if (allocated(lev)) deallocate(lev,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lev'
  if (allocated(flx_bb_abs)) deallocate(flx_bb_abs,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_abs'
  if (allocated(flx_nst_abs)) deallocate(flx_nst_abs,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_abs'
  if (allocated(htg_rate_bb)) deallocate(htg_rate_bb,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for htg_rate_bb'
  if (allocated(j_NO2)) deallocate(j_NO2,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for j_NO2'
  if (allocated(ntn_bb_mean)) deallocate(ntn_bb_mean,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_bb_mean'
  if (allocated(odal_obs_aer)) deallocate(odal_obs_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odal_obs_aer'
  if (allocated(odal_obs_mpr)) deallocate(odal_obs_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odal_obs_mpr'
  if (allocated(odal_obs_snw)) deallocate(odal_obs_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odal_obs_snw'
  if (allocated(odal_obs_bga)) deallocate(odal_obs_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odal_obs_bga'
  if (allocated(odsl_obs_aer)) deallocate(odsl_obs_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odsl_obs_aer'
  if (allocated(odsl_obs_mpr)) deallocate(odsl_obs_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odsl_obs_mpr'
  if (allocated(odsl_obs_snw)) deallocate(odsl_obs_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odsl_obs_snw'
  if (allocated(odsl_obs_bga)) deallocate(odsl_obs_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odsl_obs_bga'
  if (allocated(odxl_obs_aer)) deallocate(odxl_obs_aer,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxl_obs_aer'
  if (allocated(odxl_obs_mpr)) deallocate(odxl_obs_mpr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxl_obs_mpr'
  if (allocated(odxl_obs_snw)) deallocate(odxl_obs_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxl_obs_snw'
  if (allocated(odxl_obs_bga)) deallocate(odxl_obs_bga,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxl_obs_bga'
  if (allocated(tpt)) deallocate(tpt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tpt'
  ! Array dimensions: levp
  if (allocated(alt_ntf)) deallocate(alt_ntf,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for alt_ntf'
  if (allocated(flx_bb_dwn_dff)) deallocate(flx_bb_dwn_dff,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_dwn_dff'
  if (allocated(flx_bb_dwn_drc)) deallocate(flx_bb_dwn_drc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_dwn_drc'
  if (allocated(flx_bb_dwn)) deallocate(flx_bb_dwn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_dwn'
  if (allocated(flx_bb_net)) deallocate(flx_bb_net,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_net'
  if (allocated(flx_bb_upw)) deallocate(flx_bb_upw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_bb_upw'
  if (allocated(flx_nst_dwn)) deallocate(flx_nst_dwn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_dwn'
  if (allocated(flx_nst_net)) deallocate(flx_nst_net,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_net'
  if (allocated(flx_nst_upw)) deallocate(flx_nst_upw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_upw'
  if (allocated(flx_nst_dwn)) deallocate(flx_nst_dwn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_dwn'
  if (allocated(flx_nst_upw)) deallocate(flx_nst_upw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_nst_upw'
  if (allocated(levp)) deallocate(levp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for levp'
  if (allocated(tpt_ntf)) deallocate(tpt_ntf,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tpt_ntf'
  ! Array dimensions: plr
  if (allocated(plr)) deallocate(plr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for plr'
  if (allocated(plr_cos)) deallocate(plr_cos,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for plr_cos'
  if (allocated(plr_dgr)) deallocate(plr_dgr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for plr_dgr'
  ! Array dimensions: lev_snw, levp_snw
  if (allocated(lev_snw)) deallocate(lev_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lev_snw'
  if (allocated(levp_snw)) deallocate(levp_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for levp_snw'
  if (allocated(dns_snw)) deallocate(dns_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for dns_snw'
  if (allocated(dpt_dlt_snw)) deallocate(dpt_dlt_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for dpt_dlt_snw'
  if (allocated(dpt_ntf_snw)) deallocate(dpt_ntf_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for dpt_ntf_snw'
  if (allocated(dpt_snw)) deallocate(dpt_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for dpt_snw'
  if (allocated(mmr_mpr_snw)) deallocate(mmr_mpr_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for mmr_mpr_snw'
  if (allocated(rds_ffc_snw)) deallocate(rds_ffc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rds_ffc_snw'
  if (allocated(tpt_snw)) deallocate(tpt_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tpt_snw'
  if (allocated(tpt_ntf_snw)) deallocate(tpt_ntf_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tpt_ntf_snw'
  if (allocated(foo_snw)) deallocate(foo_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for foo_snw'
  if (allocated(alb_dff_spc_snw_dea)) deallocate(alb_dff_spc_snw_dea,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for alb_dff_spc_snw_dea'
  if (allocated(alb_drc_spc_snw_dea)) deallocate(alb_drc_spc_snw_dea,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for alb_drc_spc_snw_dea'
  ! Array dimensions: bnd,lev_snw
  if (allocated(rfl_ddm_spc_snw)) deallocate(rfl_ddm_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_ddm_spc_snw'
  if (allocated(trn_ddm_spc_snw)) deallocate(trn_ddm_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_ddm_spc_snw'
  if (allocated(rfl_dff_spc_snw)) deallocate(rfl_dff_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_dff_spc_snw'
  if (allocated(trn_dff_spc_snw)) deallocate(trn_dff_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_dff_spc_snw'
  if (allocated(rfl_drc_spc_snw)) deallocate(rfl_drc_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_drc_spc_snw'
  if (allocated(trn_ttl_drc_spc_snw)) deallocate(trn_ttl_drc_spc_snw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_ttl_drc_spc_snw'
  if (allocated(trn_drc_drc)) deallocate(trn_drc_drc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_drc_drc'
  ! Array dimensions: bnd,levp_snw
  if (allocated(rfl_dff_spc_snw_cnt)) deallocate(rfl_dff_spc_snw_cnt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_dff_spc_snw_cnt'
  if (allocated(trn_dff_spc_snw_cnt)) deallocate(trn_dff_spc_snw_cnt,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_dff_spc_snw_cnt'
  if (allocated(rfl_drc_upw_snp)) deallocate(rfl_drc_upw_snp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_drc_upw_snp'
  if (allocated(rfl_dff_upw_snp)) deallocate(rfl_dff_upw_snp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_dff_upw_snp'
  if (allocated(rfl_dff_dwn_snp)) deallocate(rfl_dff_dwn_snp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_dff_dwn_snp'
  if (allocated(trn_ttl_drc_snp)) deallocate(trn_ttl_drc_snp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_ttl_drc_snp'
  if (allocated(trn_drc_drc_snp)) deallocate(trn_drc_drc_snp,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for trn_drc_drc_snp'
  ! Array dimensions: tau
  if (allocated(tau)) deallocate(tau,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tau'
  if (allocated(tau_prs)) deallocate(tau_prs,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for tau_prs'
  ! Array dimensions: bnd,lev
  if (allocated(lgn_xpn_cff_Mie_ttl)) deallocate(lgn_xpn_cff_Mie_ttl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lgn_xpn_cff_Mie_ttl'
  if (allocated(asm_prm_HG_ttl)) deallocate(asm_prm_HG_ttl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for asm_prm_HG_ttl'
  if (allocated(flx_spc_abs)) deallocate(flx_spc_abs,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_abs'
  if (allocated(ntn_spc_mean)) deallocate(ntn_spc_mean,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_mean'
  if (allocated(odxl_spc_ttl)) deallocate(odxl_spc_ttl,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for odxl_spc_ttl'
  if (allocated(ss_alb_fct)) deallocate(ss_alb_fct,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ss_alb_fct'
  ! Array dimensions: bnd,levp
  if (allocated(flx_spc_dwn)) deallocate(flx_spc_dwn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn'
  if (allocated(flx_spc_dwn_dff)) deallocate(flx_spc_dwn_dff,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_dff'
  if (allocated(flx_spc_dwn_drc)) deallocate(flx_spc_dwn_drc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_dwn_drc'
  if (allocated(flx_spc_upw)) deallocate(flx_spc_upw,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for flx_spc_upw'
  if (allocated(lmn_spc_aa_ndr)) deallocate(lmn_spc_aa_ndr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lmn_spc_aa_ndr'
  if (allocated(ntn_spc_aa_ndr)) deallocate(ntn_spc_aa_ndr,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_ndr'
  if (allocated(ntn_spc_aa_zen)) deallocate(ntn_spc_aa_zen,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_zen'
  ! Array dimensions: plr,bnd
  if (allocated(lmn_spc_aa_sfc)) deallocate(lmn_spc_aa_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lmn_spc_aa_sfc'
  if (allocated(ntn_spc_aa_sfc)) deallocate(ntn_spc_aa_sfc,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_aa_sfc'
  ! Array dimensions: plr,levp
  if (allocated(lmn_bb_aa)) deallocate(lmn_bb_aa,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for lmn_bb_aa'
  if (allocated(ntn_bb_aa)) deallocate(ntn_bb_aa,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_bb_aa'
  ! Array dimensions: azi,plr,bnd
  if (allocated(ntn_spc_TOA)) deallocate(ntn_spc_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_TOA'
  ! Array dimensions: azi,plr,chn
  if (allocated(ntn_chn_TOA)) deallocate(ntn_chn_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_chn_TOA'
  if (allocated(rfl_chn_TOA)) deallocate(rfl_chn_TOA,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for rfl_chn_TOA'
  ! Array dimensions: azi,plr,levp
  if (allocated(ntn_spc_chn)) deallocate(ntn_spc_chn,stat=rcd)
  if(rcd /= 0) stop 'deallocate() failed for ntn_spc_chn'
  if (flg_mie) then
     ! Array dimensions: mmn,bnd,lev
     if (allocated(lgn_xpn_cff_Mie_ttl)) deallocate(lgn_xpn_cff_Mie_ttl,stat=rcd)
     if(rcd /= 0) stop 'deallocate() failed for lgn_xpn_cff_Mie_ttl'
  end if ! !flg_mie

  call exit(exit_status)
end program swnb2 ! end swnb2()
