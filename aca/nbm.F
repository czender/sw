c     $Id$

      program nbm

c     Purpose: Narrow band model employing Malkmus random line approximation
c     Effective monochromatic optical depth from gaseous line absorption is combined
c     with continuum absorption and scattering cross sections and fed to DISORT. 
c     All I/O uses netCDF interface

c     Compilation: 
c     cd $HOME/aca; make -W nbm.F OPTS=D nbm; cd -
c     cd $HOME/aca; make -W nbm.F nbm; cd -
c     cd $HOME/aca; make OPTS=D nbm; cd -

c     Debugging usage:
c     Perform full setup for all bands but only call DISORT() for #1603:
c     nbm -D 1 -E -e 0.5e-6 -d foo.nc 

c     Longwave usage:
c     NB: LW runs should never set -J, which turns off thermal emission
c     Switches -v -V -q and -j are provisional
c     nbm -D 1 -E -e 0.5e-6 -v 0.176e-6 -V 5.0e-6 -q 5.0e-9 
c     nbm -D 3 -v 5.0e-6 -V 50.0e-6 -q 1.0e-6 
c     nbm -D 3 -j 3
c     H2O vapor + Thermal emission only:
c     nbm -v 5.0e-6 -V 50.0e-6 -q 1.0e-6 -A -B -C -I -K -L -O -R -Q -U -W -X -Y -r 0.1 -z 0.5 -m 0.0
c     ncks -H -u -C -F -v flx_bb_dwn_sfc ${DATA}/aca/nbm.nc

c     Production usage:

c     Stefan Kinne's test cases:
c     NB: Thermal emission (flg_Planck) does not work with SK test cases 
c     This is because the AFGL temperature steps are too large near TOA
c     All gases, no aerosol, cloud, or thermal emission, in MLS, TRP atm with cos theta=0.8660,0.2588 alb=0.2,
c     nbm -z 0.8660 -r 0.2 -s 4 -p ${DATA}/aca/mls_afgl_73lvl.nc -d ${DATA}/tmp/sk_01.nc >&! foo_1 &
c     nbm -z 0.2588 -r 0.2 -s 4 -p ${DATA}/aca/mls_afgl_73lvl.nc -d ${DATA}/tmp/sk_02.nc >&! foo_2 &
c     nbm -z 0.8660 -r 0.2 -s 4 -p ${DATA}/aca/trp_afgl_73lvl.nc -d ${DATA}/tmp/sk_03.nc >&! foo_3 &
c     nbm -z 0.2588 -r 0.2 -s 4 -p ${DATA}/aca/trp_afgl_73lvl.nc -d ${DATA}/tmp/sk_04.nc >&! foo_4 &
c     ncks -H -C -d bnd,0.7e-6 -v flx_bb_dwn_sfc ${DATA}/tmp/sk_01.nc

c     Test cases: Remember to set CO2 vmr = 300 ppm for most of these!
c     NB: BPB used only the most common isotope of each species for these
c     The ncks command to print out diagnostic fields, for many test cases, is
c     ncks -u -H -C -v flx_bb_abs_atm,flx_bb_dwn_sfc ${DATA}/aca/nbm.nc

c     Rayleigh scattering phase function: tst_case_Rayleigh
c     nbm -A -B -C -H -I -J -K -L -O -Q -U -W -X -Y -t -D 1 -r 0.0 -E -e 0.5e-6 -p ${DATA}/aca/mls_clr.nc

c     Henyey-Greenstein phase function: tst_case_HG
c     nbm -A -B -C -H -I -J -K -L -O -Q -U -W -X -Y -P -D 1 -r 0.0 -E -e 0.5e-6 -p ${DATA}/aca/mls_clr.nc

c     O2 absorption only:
c     nbm -A -B -C -H -I -J -K -L -Q -R -U -W -X -Y -D 1 -r 0.1 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -o ${DATA}/aca/mlk_16O2.nc

c     CO2 absorption only:
c     nbm -A -B -H -I -J -K -L -O -Q -R -U -W -X -Y -D 1 -r 0.1 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -c ${DATA}/aca/mlk_12C_16O2.nc

c     H2O absorption only:
c     nbm -A -B -C -I -J -K -L -O -Q -R -U -W -X -Y -D 1 -r 0.0 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -h ${DATA}/aca/mlk_1H2_16O.nc
c     nbm -A -B -C -I -J -K -L -O -Q -R -U -W -X -Y -D 1 -r 0.0 -z 0.8660254 -S 1370.0 -p ${DATA}/aca/mls_clr.nc -h ${DATA}/aca/mlk_1H2_16O.nc
c     nbm -A -B -C -I -J -K -L -O -Q -R -U -W -X -Y -D 1 -r 0.0 -z 0.8660254 -S 1370.0 -p ${DATA}/aca/mls_icrccm_92lvl.nc -h ${DATA}/aca/mlk_H2O.nc

c     O3 absorption only:
c     nbm -A -B -C -H -I -J -K -L -O -Q -R -U -X -Y -D 1 -r 0.0 -z 0.5 -S 1370.0 -p ${DATA}/aca/mls_clr.nc 

c     Rayleigh scattering only:
c     nbm -A -B -C -H -I -J -K -L -O -Q -U -W -X -Y -r 0.0 -z 0.5

c     Mie scattering/absorption only:
c     nbm -C -H -J -K -O -Q -R -U -W -X -Y -r 0.0 -z 0.5 -m 0.1

c     Full physics: 18 layer MLS, cos theta=0.5, A=0.2, CWP=100
c     nbm -z 0.5 -r 0.2 -m 0.1 -s 4 
c     ncks -C -H -d bnd,1.27e-6 -v flx_bb_dwn_sfc,odxc_spc_O2O2,odxc_spc_O2N2 ${DATA}/aca/nbm.nc

c     Sensitivity studies: Everything except given process, then that process alone
c     Control:
c     nbm -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm.nc
c     ncks -C -H -v flx_bb_dwn_sfc,flx_bb_abs_atm ${DATA}/tmp/nbm.nc

c     H2OH2O absorption:
c     nbm -Q -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_H2OH2O.nc
c     nbm -A -B -C -H -I -J -K -L -O -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_H2OH2O.nc

c     O2-O2 absorption:
c     nbm -K -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2O2.nc
c     nbm -A -B -C -H -I -J -L -O -R -Q -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2O2.nc

c     O2-N2 absorption:
c     nbm -U -r .1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2N2.nc
c     nbm -A -B -C -H -I -J -K -L -O -R -Q -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2N2.nc

c     O2 absorption:
c     nbm -O -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2.nc
c     nbm -A -B -C -H -I -J -K -L -Q -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O2.nc

c     O3 absorption:
c     nbm -W -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O3.nc
c     nbm -A -B -C -H -I -J -K -L -O -Q -R -U -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_O3.nc

c     OH absorption:
c     nbm -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_OH.nc
c     nbm -A -B -C -H -I -J -K -L -O -Q -R -U -W -X -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_OH.nc

c     NO2 absorption:
c     nbm -X -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_NO2.nc
c     nbm -A -B -C -H -I -J -K -L -O -Q -R -U -W -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_NO2.nc

c     CO2 absorption:
c     nbm -C -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_CO2.nc
c     nbm -A -B -H -I -J -K -L -O -R -Q -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_CO2.nc

c     H2O absorption:
c     nbm -H -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_H2O.nc
c     nbm -A -B -C -I -J -K -L -O -R -Q -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_H2O.nc

c     Rayleigh scattering:
c     nbm -R -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_Ray.nc
c     nbm -A -B -C -H -I -J -K -L -O -Q -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_Ray.nc

c     Thermal absorption/emission:
c     nbm -J -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_Plk.nc
c     nbm -A -B -C -H -I -K -L -O -Q -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_Plk.nc

c     Aerosol scattering/absorption:
c     nbm -A -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_aer.nc
c     nbm -B -C -H -I -J -K -L -O -Q -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_aer.nc

c     Background aerosol scattering/absorption:
c     nbm -B -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_bga.nc
c     nbm -A -C -H -I -J -K -L -O -Q -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.0 -d ${DATA}/tmp/nbm_bga.nc

c     Liquid scattering/absorption:
c     nbm -L -r .1 -z .5 -m .01 -f -d ${DATA}/tmp/nbm_lqd.nc
c     nbm -A -B -C -H -I -J -K -O -Q -R -U -W -X -Y -r .1 -z .5 -m .01 -f -d ${DATA}/tmp/nbm_lqd.nc

c     Ice scattering/absorption:
c     nbm -I -r 0.1 -z 0.5 -m 0.01 -F -d ${DATA}/tmp/nbm_ice.nc
c     nbm -A -B -C -H -J -K -L -O -Q -R -U -W -X -Y -r 0.1 -z 0.5 -m 0.01 -F -d ${DATA}/tmp/nbm_ice.nc

c     ARESE simulations: 
c     nbm -Q -p ${DATA}/arese/clm/951011_1200_arese_clm.nc -d ${DATA}/tmp/951011_1200_arese_mdl_clr_aer.nc
c     nbm -Q -p ${DATA}/arese/clm/951015_1200_arese_clm.nc -d ${DATA}/tmp/951015_1200_arese_mdl_clr_aer.nc
c     ncks -C -H -v flx_bb_dwn_sfc,flx_bb_abs_atm ${DATA}/tmp/foo.nc
c     ncks -C -H -d bnd,0.5e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ${DATA}/tmp/foo.nc
c     ncks -C -H -d bnd,0.413e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ${DATA}/tmp/foo.nc
c     ncks -C -H -d bnd,0.860e-6 -v wvl_obs_aer,wvl_obs_bga,odxc_obs_aer,odxc_obs_bga,odxc_spc_aer,odxc_spc_bga ${DATA}/tmp/foo.nc

c     Dust simulations
c     j_NO2
c     nbm -A -z 0.866 -p ${DATA}/dst/arese_951011_1200_clm_clr_aer_dst.nc -d ${DATA}/dst/arese_951011_1200_mdl_clr_cln.nc >&! foo.cln &
c     nbm -z 0.866 -a ${DATA}/aca/aer_mineral_dust.nc -p ${DATA}/dst/arese_951011_1200_clm_clr_aer_dst.nc -d ${DATA}/dst/arese_951011_1200_mdl_clr_aer_dst.nc >&! foo.dst &
c     nbm -z 0.866 -a ${DATA}/aca/aer_h2so4_300K.nc -p ${DATA}/dst/arese_951011_1200_clm_clr_aer_slf.nc -d ${DATA}/dst/arese_951011_1200_mdl_clr_aer_slf.nc >&! foo.slf &
c     ncdiff -v j_NO2 ${DATA}/dst/arese_951011_1200_mdl_clr_aer_dst.nc ${DATA}/dst/arese_951011_1200_mdl_clr_cln.nc ${DATA}/dst/arese_951011_1200_mdl_clr_dstmcln.nc
c     ncks -C -H -v j_NO2 -u ${DATA}/dst/arese_951011_1200_mdl_clr_cln.nc
c     ncks -C -H -v j_NO2 -u ${DATA}/dst/arese_951011_1200_mdl_clr_aer_dst.nc
c     ncks -C -H -v j_NO2 -u ${DATA}/dst/arese_951011_1200_mdl_clr_aer_slf.nc
c     ncks -C -H -v j_NO2 -u ${DATA}/dst/arese_951011_1200_mdl_clr_dstmcln.nc

c     The code makes heavy use of acronyms and abbreviations in variable names
c     The master abbreviation list is at ${HOME}/crr/abb.tex 
c     The abbreviation list URL is http://www.cgd.ucar.edu/~zender/abb

c     Code proceeds as follows:
c     Section 1: Initialization
c     Section 2: Main computation loop over spectral bands
c     Section 3: Postprocessing

c     Begin section 1: Initialization
c     initialization of program defaults
c     processing of command line switches
c     netCDF input of atmospheric thermodynamic and cloud profile
c     netCDF input of gaseous absorption and aerosol scattering and absorption characteristics
c     override netCDF input data according to command line switches
c     band-independent DISORT() initialization
c     band loops for figuring out where to splice input data together
c     initialization of  band quantities that are level independent
c     initialization of level quantities that are band independent
c     a band-level loop for initializing band-level quantities

      use netcdf ! [mdl] netCDF interface
      implicit none
c     Parameters
      character*(*) CVS_Id

c     Add to lbl
      integer wvl_grd_typ       ! Wavelength grid type
      integer wvl_idx           ! counting index
      integer wvl_nbr           ! # wavelength bins
      real wvl_grd_min          ! [m] Minimum wavelength
      real wvl_grd_max          ! [m] Maximum wavelength
      real wvl_grd_rsn          ! [m] Resolution

c     NB: DISORT uses fortran unit 99 for its own purposes
      integer azi_nbr_max       ! # user azimuthal angles
      integer bnd_nbr_CO2_max
      integer bnd_nbr_H2OH2O_max
      integer bnd_nbr_H2O_max
      integer bnd_nbr_Mlk_max
      integer bnd_nbr_NO2_max
      integer bnd_nbr_O2O2_max
      integer bnd_nbr_O2_max
      integer bnd_nbr_O3_max
      integer bnd_nbr_OH_max
      integer bnd_nbr_CH4_max
      integer bnd_nbr_N2O_max
      integer bnd_nbr_aer_max
      integer bnd_nbr_bga_max
      integer bnd_nbr_ice_max
      integer bnd_nbr_lqd_max
      integer bnd_nbr_max
      integer bnd_nbr_nst_max
      integer lev_nbr_max
      integer levp_nbr_max
      integer plr_nbr_max       ! # user polar angles
      integer str_nbr_max       ! # computational polar angles
      integer tau_nbr_max       ! # user optical depths
      real tpt_Malkmus_fit   ! reference temperature for Malkmus parameters

      parameter(CVS_Id='$Id$',
     $     azi_nbr_max=16,      ! Arbitrary (must match DISORT parameter)
     $     bnd_nbr_CO2_max=5000, ! 5 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_H2OH2O_max=8192, ! Broadband averages from TGC98
     $     bnd_nbr_H2O_max=2500, ! 10 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_Mlk_max=10000, ! max(bnd_nbr_CO2_max,bnd_nbr_H2O_max,bnd_nbr_O2_max=1590,bnd_nbr_OH_max, ...)
     $     bnd_nbr_NO2_max=750, ! 
     $     bnd_nbr_O2O2_max=4086, ! Greenblatt et al. (1992) resolution + SPS98 1.27 um band
     $     bnd_nbr_O2_max=5000, ! 5 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_O3_max=158,  ! WMO85 resolution
     $     bnd_nbr_OH_max=5000, ! 5 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_CH4_max=5000, ! 5 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_N2O_max=5000, ! 5 cm-1 resolution from 0-25000 cm-1
     $     bnd_nbr_aer_max=200,  ! From 0.2--5.0 microns
     $     bnd_nbr_bga_max=200,  ! From 0.2--5.0 microns
     $     bnd_nbr_ice_max=194, ! War84 resolution
     $     bnd_nbr_lqd_max=74,  ! ? resolution
c     $     bnd_nbr_max=1690, ! 1 nm resolution from 176 nm to 5.0 um
     $     bnd_nbr_max=4824,    ! 1 nm resolution from 176 nm to 5.0 um
     $     bnd_nbr_nst_max=235, ! new FSBR resolution
     $     lev_nbr_max=110,     ! Roughly 10 mb resolution from 1010 mb to TOA
     $     levp_nbr_max=lev_nbr_max+1,
     $     plr_nbr_max=16,      ! Arbitrary (must match DISORT parameter)
     $     str_nbr_max=16,      ! Arbitrary (must match DISORT parameter)
     $     tau_nbr_max=lev_nbr_max+1,
     $     tpt_Malkmus_fit=250.0) ! [K]

c     Commons
#include <parameter.com> /* Fundamental and derived physical constants */ 
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation/interpolation handling */
#include <wvl_grd.com> /* Wavelength grid constants */
c     Input Arguments
c     Input/Output Arguments
c     Output Arguments
c     Local workspace
      character argv*80
      character aer_sng*100
      character azi_sng*100
      character bga_sng*100
      character char_foo*100
      character cmd_ln*200
      character fl_CO2*80
      character fl_H2OH2O*80
      character fl_H2O*80
      character fl_NO2*80
      character fl_O2*80
      character fl_O3*80
      character fl_O2O2*80
      character fl_OH*80
      character fl_CH4*80
      character fl_N2O*80
      character fl_aer*80
      character fl_bga*80
      character fl_clm*80
      character fl_ice*80
      character fl_lqd*80
      character fl_nst*80
      character fl_out*80
      character fl_slr*80
      character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
      character opt_dep_sng*100
      character plr_sng*100
      character prf_sng*100
      character prg_ID*200
      character CVS_Date*28
      character CVS_Revision*80
      character src_rfr_sng*200
      character str_sng*200
      character stt_CO2*200
      character stt_H2O*200
      character stt_H2OH2O*200
      character stt_Herzberg*200
      character stt_NO2*200
      character stt_O2*200
      character stt_O3*200
      character stt_O2O2*200
      character stt_O2N2*200
      character stt_OH*200
      character stt_CH4*200
      character stt_N2O*200
      character stt_Planck*200
      character stt_Rayleigh*200
      character stt_aer*200
      character stt_bga*200
      character stt_ice*200
      character stt_lqd*200
      character stt_slr*200
      character stt_flt_nst*200
      character stt_top_lvl*200

      integer arg
      integer exit_status       ! Program exit status
      integer arg_nbr
      integer int_foo
      integer rcd               ! Return success code

      logical cmd_ln_mpc_CWP
      logical cmd_ln_odxc_obs_aer
      logical cmd_ln_alb
      logical cmd_ln_slr_zen_ngl_cos
      logical cmd_ln_slr_cst
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
      logical flg_N2O
      logical flg_Planck
      logical flg_Rayleigh
      logical flg_aer
      logical flg_bga
      logical flg_ice
      logical flg_lqd
      logical flt_nst
      logical force_ice_phz
      logical force_lqd_phz
      logical logical_foo
      logical mnt_ncr           ! Monotonic and increasing flag
      logical single_bnd_computation
      logical sv_cmp_tau
      logical sv_cmp_plr_ngl
      logical sv_ntn
      logical top_lvl
      logical tst_case_HG
      logical tst_case_Rayleigh

      integer bnd_nbr_CO2       ! Dimension size
      integer bnd_nbr_H2OH2O    ! Dimension size
      integer bnd_nbr_H2O       ! Dimension size
      integer bnd_nbr_NO2       ! Dimension size
      integer bnd_nbr_O2        ! Dimension size
      integer bnd_nbr_O3        ! Dimension size
      integer bnd_nbr_O2O2      ! Dimension size
      integer bnd_nbr_OH        ! Dimension size
      integer bnd_nbr_CH4        ! Dimension size
      integer bnd_nbr_N2O        ! Dimension size
      integer bnd_nbr_aer       ! Dimension size
      integer bnd_nbr_bga       ! Dimension size
      integer azi_nbr           ! Dimension size
      integer bnd_nbr           ! Dimension size
      integer bnd_nbr_ice       ! Dimension size
      integer lev_nbr           ! Dimension size
      integer levp_nbr          ! Dimension size
      integer bnd_nbr_lqd       ! Dimension size
      integer nbr_mom           ! Dimension size
      integer bnd_nbr_nst       ! Dimension size
      integer plr_nbr           ! Dimension size
      integer str_nbr           ! Dimension size
      integer tau_nbr           ! Dimension size
      
      integer azi_idx           ! Counting index
      integer bnd_idx           ! Counting index
      integer lev_idx           ! Counting index
      integer mom_idx           ! Counting index
      integer plr_idx           ! Counting index
      integer tau_idx           ! Counting index
      
      integer azi_dim_id        ! Dimension ID for azi
      integer bnd_dim_id        ! Dimension ID for bnd
      integer fll_mode_old      ! Old fill mode
      integer grd_dim_id        ! Dimension ID for grid
      integer dim_id_vec_foo(4)
      integer dim_azi_plr_bnd_levp(4)
      integer dim_azi_plr_levp(3)
      integer dim_bnd_lev(2)
      integer dim_bnd_levp(2)
      integer dim_plr_bnd(2)
      integer dim_plr_bnd_levp(3)
      integer dim_plr_levp(2)
      integer lev_dim_id        ! Dimension ID for lev
      integer levp_dim_id       ! Dimension ID for levp
      integer nc_id             ! File handle
      integer plr_dim_id        ! Dimension ID for plr
      integer tau_dim_id        ! Dimension ID for tau
      
c     netCDF output variables not contained in input files
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
      integer bnd_var_id
      integer flx_abs_atm_rdr_id
      integer flx_bb_abs_atm_id
      integer flx_bb_abs_id
      integer flx_bb_abs_sfc_id
      integer flx_bb_abs_ttl_id
      integer flx_bb_dwn_TOA_id
      integer flx_bb_dwn_dff_id
      integer flx_bb_dwn_drc_id
      integer flx_bb_dwn_id
      integer flx_bb_dwn_sfc_id
      integer flx_bb_net_id
      integer flx_bb_upw_id
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
      integer flx_spc_dwn_id
      integer flx_spc_dwn_sfc_id
      integer flx_spc_pht_dwn_sfc_id
      integer flx_spc_upw_id
      integer flx_slr_frc_id
      integer htg_rate_bb_id
      integer j_NO2_id
      integer j_spc_NO2_sfc_id
      integer lev_id            ! coordinate ID
      integer levp_id           ! coordinate ID
      integer ntn_bb_aa_id
      integer ntn_bb_mean_id
      integer ntn_spc_aa_ndr_id
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
      integer odxc_spc_N2O_id
      integer odxc_spc_Ray_id
      integer odac_spc_aer_id
      integer odac_spc_bga_id
      integer odac_spc_ice_id
      integer odac_spc_lqd_id
      integer odxc_spc_aer_id
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
      integer trn_spc_atm_N2O_id
      integer trn_spc_atm_Ray_id
      integer trn_spc_atm_aer_id
      integer trn_spc_atm_bga_id
      integer trn_spc_atm_ice_id
      integer trn_spc_atm_lqd_id
      integer trn_spc_atm_ttl_id
      integer wvl_ctr_id
      integer wvl_grd_id
      integer wvl_max_id
      integer wvl_min_id
      integer wvl_dlt_id
      integer wvn_ctr_id
      integer wvn_grd_id
      integer wvn_max_id
      integer wvn_min_id
      integer wvn_dlt_id
      integer alt_id
      integer alt_ntf_id

c     WMO input variables
      integer abs_xsx_O2_id
      integer abs_xsx_O3_id
      integer wvl_grd_O3_id

c     O2-O2 input variables
      integer abs_xsx_O2O2_id
      integer wvl_grd_O2O2_id

c     NO2 input variables
      integer abs_xsx_NO2_id
      integer qnt_yld_NO2_id
      integer wvl_grd_NO2_id

c     H2OH2O input variables
      integer abs_xsx_H2OH2O_id
      integer wvl_grd_H2OH2O_id

c     Narrow band H2O input variables
      integer A_phi_H2O_id
      integer A_psi_H2O_id
      integer B_phi_H2O_id
      integer B_psi_H2O_id
      integer S_d_abs_cff_mss_H2O_id
      integer S_p_abs_cff_mss_H2O_id
      integer wvl_grd_H2O_id

c     Narrow band CO2 input variables
      integer A_phi_CO2_id
      integer A_psi_CO2_id
      integer B_phi_CO2_id
      integer B_psi_CO2_id
      integer S_d_abs_cff_mss_CO2_id
      integer S_p_abs_cff_mss_CO2_id
      integer wvl_grd_CO2_id

c     Narrow band OH input variables
      integer A_phi_OH_id
      integer A_psi_OH_id
      integer B_phi_OH_id
      integer B_psi_OH_id
      integer S_d_abs_cff_mss_OH_id
      integer S_p_abs_cff_mss_OH_id
      integer wvl_grd_OH_id

c     Narrow band OH input variables
      integer A_phi_CH4_id
      integer A_psi_CH4_id
      integer B_phi_CH4_id
      integer B_psi_CH4_id
      integer S_d_abs_cff_mss_CH4_id
      integer S_p_abs_cff_mss_CH4_id
      integer wvl_grd_CH4_id

c     Narrow band OH input variables
      integer A_phi_N2O_id
      integer A_psi_N2O_id
      integer B_phi_N2O_id
      integer B_psi_N2O_id
      integer S_d_abs_cff_mss_N2O_id
      integer S_p_abs_cff_mss_N2O_id
      integer wvl_grd_N2O_id

c     Narrow band O2 input variables
      integer A_phi_O2_id
      integer A_psi_O2_id
      integer B_phi_O2_id
      integer B_psi_O2_id
      integer S_d_abs_cff_mss_O2_id
      integer S_p_abs_cff_mss_O2_id
      integer wvl_grd_O2_id

c     Ice water input variables
      integer asm_prm_ice_id
      integer abs_cff_mss_ice_id
      integer sca_cff_mss_ice_id
      integer wvl_grd_ice_id

c     Liquid water input variables
      integer asm_prm_lqd_id
      integer abs_cff_mss_lqd_id
      integer sca_cff_mss_lqd_id
      integer wvl_grd_lqd_id

c     Instrument input variables
      integer nst_SRF_id
      integer wvl_grd_nst_id

c     Aerosol input variables
      integer asm_prm_aer_id
      integer abs_cff_mss_aer_id
      integer ext_cff_mss_aer_id
      integer sca_cff_mss_aer_id
      integer wvl_grd_aer_id

c     Background aerosol input variables
      integer asm_prm_bga_id
      integer abs_cff_mss_bga_id
      integer sca_cff_mss_bga_id
      integer wvl_grd_bga_id

c     CLM input variables
      integer RH_lqd_id
      integer alb_sfc_NIR_drc_id
      integer alb_sfc_NIR_dff_id
      integer alb_sfc_vsb_drc_id
      integer alb_sfc_vsb_dff_id
      integer slr_zen_ngl_cos_id
      integer prs_dlt_id
      integer frc_ice_id
      integer frc_ice_ttl_id
      integer grv_id
      integer lat_dgr_id
      integer lcl_time_hr_id
      integer lcl_yr_day_id
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
      integer mpl_N2O_id
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
      integer prs_id
      integer prs_ntf_id
      integer q_CO2_id
      integer q_H2O_id
      integer q_O2_id
      integer q_OH_id
      integer q_CH4_id
      integer q_N2O_id
      integer spc_heat_mst_air_id
      integer tpt_id
      integer tpt_ntf_id
      integer tpt_skn_id
      integer xnt_fac_id
      integer wvl_obs_aer_id
      integer wvl_obs_bga_id
      integer alt_cld_btm_id
      integer alt_cld_thick_id

c     netCDF output variables
      real nrg_pht(bnd_nbr_max)
      real flx_spc_pht_dwn_sfc(bnd_nbr_max)
      real abs_bb_SAS
      real abs_bb_atm
      real abs_bb_sfc
      real abs_nst_SAS
      real abs_nst_atm
      real abs_nst_sfc
      real abs_spc_SAS(bnd_nbr_max)
      real abs_spc_atm(bnd_nbr_max)
      real abs_spc_sfc(bnd_nbr_max)
      real alb_sfc
      real asm_prm_ttl(bnd_nbr_max,lev_nbr_max)
      real azi(azi_nbr_max)     ! Coordinate variable
      real azi_dgr(azi_nbr_max)
      real bnd(bnd_nbr_max)     ! Coordinate variable
      real bnd_var(bnd_nbr_max) ! Extra variable
      real flx_abs_atm_rdr(bnd_nbr_max)
      real flx_bb_abs(lev_nbr_max)
      real flx_bb_abs_atm
      real flx_bb_abs_sfc
      real flx_bb_abs_ttl
      real flx_bb_dwn(levp_nbr_max)
      real flx_bb_dwn_TOA
      real flx_bb_dwn_dff(levp_nbr_max)
      real flx_bb_dwn_drc(levp_nbr_max)
      real flx_bb_dwn_sfc
      real flx_bb_net(levp_nbr_max)
      real flx_bb_up(levp_nbr_max)
      real flx_nst_abs(lev_nbr_max)
      real flx_nst_abs_atm
      real flx_nst_abs_sfc
      real flx_nst_abs_ttl
      real flx_nst_dwn(levp_nbr_max)
      real flx_nst_dwn_TOA
      real flx_nst_dwn_sfc
      real flx_nst_net(levp_nbr_max)
      real flx_nst_up(levp_nbr_max)
      real flx_spc_abs(bnd_nbr_max,lev_nbr_max)
      real flx_spc_abs_SAS(bnd_nbr_max)
      real flx_spc_abs_atm(bnd_nbr_max)
      real flx_spc_abs_sfc(bnd_nbr_max)
      real flx_spc_act_pht_TOA(bnd_nbr_max)
      real flx_spc_act_pht_sfc(bnd_nbr_max)
      real flx_spc_dwn(bnd_nbr_max,levp_nbr_max)
      real flx_spc_dwn_TOA(bnd_nbr_max)
      real flx_spc_dwn_dff(bnd_nbr_max,levp_nbr_max)
      real flx_spc_dwn_drc(bnd_nbr_max,levp_nbr_max)
      real flx_spc_dwn_sfc(bnd_nbr_max)
      real flx_spc_up(bnd_nbr_max,levp_nbr_max)
      real flx_slr_frc(bnd_nbr_max)
      real htg_rate_bb(lev_nbr_max)
      real j_NO2(lev_nbr_max)
      real j_spc_NO2_sfc(bnd_nbr_max)
      real lev(lev_nbr_max)     ! coordinate variable
      real levp(levp_nbr_max)   ! coordinate variable
      real ntn_bb_aa(plr_nbr_max,levp_nbr_max)
      real ntn_bb_mean(lev_nbr_max)
      real ntn_spc_aa_ndr(bnd_nbr_max,levp_nbr_max)
      real ntn_spc_aa_ndr_sfc(bnd_nbr_max)
      real ntn_spc_aa_sfc(plr_nbr_max,bnd_nbr_max)
      real ntn_spc_aa_zen(bnd_nbr_max,levp_nbr_max)
      real ntn_spc_aa_zen_sfc(bnd_nbr_max)
      real ntn_spc_chn(azi_nbr_max,plr_nbr_max,levp_nbr_max)
      real ntn_spc_mean(bnd_nbr_max,lev_nbr_max)
      real odxc_spc_CO2(bnd_nbr_max)
      real odxc_spc_H2OH2O(bnd_nbr_max)
      real odxc_spc_H2O(bnd_nbr_max)
      real odxc_spc_NO2(bnd_nbr_max)
      real odxc_spc_O2(bnd_nbr_max)
      real odxc_spc_O3(bnd_nbr_max)
      real odxc_spc_O2O2(bnd_nbr_max)
      real odxc_spc_O2N2(bnd_nbr_max)
      real odxc_spc_OH(bnd_nbr_max)
      real odxc_spc_CH4(bnd_nbr_max)
      real odxc_spc_N2O(bnd_nbr_max)
      real odxc_spc_Ray(bnd_nbr_max)
      real odac_spc_aer(bnd_nbr_max)
      real odac_spc_bga(bnd_nbr_max)
      real odac_spc_ice(bnd_nbr_max)
      real odac_spc_lqd(bnd_nbr_max)
      real odxc_spc_aer(bnd_nbr_max)
      real odxc_spc_bga(bnd_nbr_max)
      real odxc_spc_ice(bnd_nbr_max)
      real odxc_spc_lqd(bnd_nbr_max)
      real odxc_spc_ttl(bnd_nbr_max)
      real odxl_spc_ttl(bnd_nbr_max,lev_nbr_max)
      real plr(plr_nbr_max)     ! coordinate variable
      real plr_cos(plr_nbr_max)
      real plr_dgr(plr_nbr_max)
      real rfl_bb_SAS
      real rfl_bb_sfc
      real rfl_nst_SAS
      real rfl_nst_sfc
      real rfl_spc_SAS(bnd_nbr_max)
      real ss_alb_fct(bnd_nbr_max,lev_nbr_max)
      real tau(tau_nbr_max)     ! coordinate variable
      real tau_prs(tau_nbr_max)
      real trn_bb_atm
      real trn_nst_atm
      real trn_spc_atm_CO2(bnd_nbr_max)
      real trn_spc_atm_H2OH2O(bnd_nbr_max)
      real trn_spc_atm_H2O(bnd_nbr_max)
      real trn_spc_atm_NO2(bnd_nbr_max)
      real trn_spc_atm_O2(bnd_nbr_max)
      real trn_spc_atm_O3(bnd_nbr_max)
      real trn_spc_atm_O2O2(bnd_nbr_max)
      real trn_spc_atm_O2N2(bnd_nbr_max)
      real trn_spc_atm_OH(bnd_nbr_max)
      real trn_spc_atm_CH4(bnd_nbr_max)
      real trn_spc_atm_N2O(bnd_nbr_max)
      real trn_spc_atm_Ray(bnd_nbr_max)
      real trn_spc_atm_aer(bnd_nbr_max)
      real trn_spc_atm_bga(bnd_nbr_max)
      real trn_spc_atm_ice(bnd_nbr_max)
      real trn_spc_atm_lqd(bnd_nbr_max)
      real trn_spc_atm_ttl(bnd_nbr_max)
      real wvl(bnd_nbr_max)
      real wvl_ctr(bnd_nbr_max)
      real wvl_grd(bnd_nbr_max+1)
      real wvl_max(bnd_nbr_max)
      real wvl_min(bnd_nbr_max)
      real wvl_dlt(bnd_nbr_max)
      real wvn(bnd_nbr_max)
      real wvn_ctr(bnd_nbr_max)
      real wvn_grd(bnd_nbr_max+1)
      real wvn_max(bnd_nbr_max)
      real wvn_min(bnd_nbr_max)
      real wvn_dlt(bnd_nbr_max)

c     WMO input variables
      real abs_xsx_O2_obs(bnd_nbr_O3_max)
      real abs_xsx_O3_obs(bnd_nbr_O3_max)
      real abs_xsx_O2(bnd_nbr_max)
      real abs_xsx_O3(bnd_nbr_max)
      real wvl_grd_O3(bnd_nbr_O3_max+1)

c     O2-O2 input variables
      real abs_xsx_O2O2_obs(bnd_nbr_O2O2_max)
      real abs_xsx_O2O2(bnd_nbr_max)
      real wvl_grd_O2O2(bnd_nbr_O2O2_max+1)

c     NO2 input variables
      real abs_xsx_NO2_obs(bnd_nbr_NO2_max)
      real qnt_yld_NO2_obs(bnd_nbr_NO2_max)
      real abs_xsx_NO2(bnd_nbr_max)
      real qnt_yld_NO2(bnd_nbr_max)
      real wvl_grd_NO2(bnd_nbr_NO2_max+1)

c     H2OH2O input variables
      real abs_xsx_H2OH2O_obs(bnd_nbr_H2OH2O_max)
      real abs_xsx_H2OH2O(bnd_nbr_max)
      real wvl_grd_H2OH2O(bnd_nbr_H2OH2O_max+1)

c     Narrow band H2O input variables
      real A_phi_H2O(bnd_nbr_max)
      real A_psi_H2O(bnd_nbr_max)
      real B_phi_H2O(bnd_nbr_max)
      real B_psi_H2O(bnd_nbr_max)
      real S_d_abs_cff_mss_H2O(bnd_nbr_max)
      real S_p_abs_cff_mss_H2O(bnd_nbr_max)
      real wvl_grd_H2O(bnd_nbr_H2O_max+1)

c     Narrow band CO2 input variables
      real A_phi_CO2(bnd_nbr_max)
      real A_psi_CO2(bnd_nbr_max)
      real B_phi_CO2(bnd_nbr_max)
      real B_psi_CO2(bnd_nbr_max)
      real S_d_abs_cff_mss_CO2(bnd_nbr_max)
      real S_p_abs_cff_mss_CO2(bnd_nbr_max)
      real wvl_grd_CO2(bnd_nbr_CO2_max+1)

c     Narrow band OH input variables
      real A_phi_OH(bnd_nbr_max)
      real A_psi_OH(bnd_nbr_max)
      real B_phi_OH(bnd_nbr_max)
      real B_psi_OH(bnd_nbr_max)
      real S_d_abs_cff_mss_OH(bnd_nbr_max)
      real S_p_abs_cff_mss_OH(bnd_nbr_max)
      real wvl_grd_OH(bnd_nbr_OH_max+1)

c     Narrow band OH input variables
      real A_phi_CH4(bnd_nbr_max)
      real A_psi_CH4(bnd_nbr_max)
      real B_phi_CH4(bnd_nbr_max)
      real B_psi_CH4(bnd_nbr_max)
      real S_d_abs_cff_mss_CH4(bnd_nbr_max)
      real S_p_abs_cff_mss_CH4(bnd_nbr_max)
      real wvl_grd_CH4(bnd_nbr_CH4_max+1)

c     Narrow band OH input variables
      real A_phi_N2O(bnd_nbr_max)
      real A_psi_N2O(bnd_nbr_max)
      real B_phi_N2O(bnd_nbr_max)
      real B_psi_N2O(bnd_nbr_max)
      real S_d_abs_cff_mss_N2O(bnd_nbr_max)
      real S_p_abs_cff_mss_N2O(bnd_nbr_max)
      real wvl_grd_N2O(bnd_nbr_N2O_max+1)

c     Narrow band O2 input variables
      real A_phi_O2(bnd_nbr_max)
      real A_psi_O2(bnd_nbr_max)
      real B_phi_O2(bnd_nbr_max)
      real B_psi_O2(bnd_nbr_max)
      real S_d_abs_cff_mss_O2(bnd_nbr_max)
      real S_p_abs_cff_mss_O2(bnd_nbr_max)
      real wvl_grd_O2(bnd_nbr_O2_max+1)

c     Ice water input variables
      real asm_prm_ice_obs(bnd_nbr_ice_max)
      real abs_cff_mss_ice_obs(bnd_nbr_ice_max)
      real sca_cff_mss_ice_obs(bnd_nbr_ice_max)
      real asm_prm_ice(bnd_nbr_max)
      real abs_cff_mss_ice(bnd_nbr_max)
      real sca_cff_mss_ice(bnd_nbr_max)
      real wvl_grd_ice(bnd_nbr_ice_max+1)

c     Liquid water input variables
      real asm_prm_lqd_obs(bnd_nbr_lqd_max)
      real abs_cff_mss_lqd_obs(bnd_nbr_lqd_max)
      real sca_cff_mss_lqd_obs(bnd_nbr_lqd_max)
      real asm_prm_lqd(bnd_nbr_max)
      real abs_cff_mss_lqd(bnd_nbr_max)
      real sca_cff_mss_lqd(bnd_nbr_max)
      real wvl_grd_lqd(bnd_nbr_lqd_max+1)

c     Instrument input variables
      real nst_SRF_obs(bnd_nbr_nst_max)
      real nst_SRF(bnd_nbr_max)
      real wvl_grd_nst(bnd_nbr_nst_max+1)

c     Aerosol input variables
      real asm_prm_aer_obs(bnd_nbr_aer_max)
      real abs_cff_mss_aer_obs(bnd_nbr_aer_max)
      real sca_cff_mss_aer_obs(bnd_nbr_aer_max)
      real ext_cff_mss_aer_obs(bnd_nbr_aer_max)
      real asm_prm_aer(bnd_nbr_max)
      real abs_cff_mss_aer(bnd_nbr_max)
      real ext_cff_mss_aer(bnd_nbr_max)
      real sca_cff_mss_aer(bnd_nbr_max)
      real wvl_grd_aer(bnd_nbr_aer_max+1)

c     Background aerosol input variables
      real asm_prm_bga_obs(bnd_nbr_bga_max)
      real abs_cff_mss_bga_obs(bnd_nbr_bga_max)
      real sca_cff_mss_bga_obs(bnd_nbr_bga_max)
      real asm_prm_bga(bnd_nbr_max)
      real abs_cff_mss_bga(bnd_nbr_max)
      real sca_cff_mss_bga(bnd_nbr_max)
      real wvl_grd_bga(bnd_nbr_bga_max+1)

c     CLM input variables
      double precision lat_dgr
      double precision lcl_time_hr
      double precision lcl_yr_day
      double precision slr_zen_ngl_cos

      real RH_lqd(lev_nbr_max)
      real alb_sfc_NIR_drc
      real alb_sfc_NIR_dff
      real alb_sfc_vsb_drc
      real alb_sfc_vsb_dff
      real prs_dlt(lev_nbr_max)
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
      real mpl_N2O(lev_nbr_max)
      real mpl_aer(lev_nbr_max)
      real mpl_bga(lev_nbr_max)
      real mpl_mst_air(lev_nbr_max)
      real npl_NO2(lev_nbr_max)
      real npl_O2(lev_nbr_max)
      real npl_O3(lev_nbr_max)
      real npl_O2O2(lev_nbr_max)
      real npl_H2OH2O(lev_nbr_max)
      real odxc_obs_aer
      real odxc_obs_bga
      real odal_obs_aer(lev_nbr_max)
      real odal_obs_bga(lev_nbr_max)
      real odsl_obs_aer(lev_nbr_max)
      real odsl_obs_bga(lev_nbr_max)
      real odxl_obs_aer(lev_nbr_max)
      real odxl_obs_bga(lev_nbr_max)
      real prs(lev_nbr_max)
      real prs_ntf(levp_nbr_max)
      real q_CO2(lev_nbr_max)
      real q_H2O(lev_nbr_max)
      real q_O2(lev_nbr_max)
      real q_OH(lev_nbr_max)
      real q_CH4(lev_nbr_max)
      real q_N2O(lev_nbr_max)
      real spc_heat_mst_air(lev_nbr_max)
      real tpt(lev_nbr_max)
      real tpt_ntf(levp_nbr_max)
      real tpt_skn
      real wvl_obs_aer
      real wvl_obs_bga
      real xnt_fac
      real alt(lev_nbr_max)
      real alt_cld_btm
      real alt_cld_thick
      real alt_ntf(levp_nbr_max)

c     Local arrays
      double precision slr_zen_ngl_cos_cmd_ln
      double precision pi

      integer aer_lvl_nbr
      integer bnd_idx_aer
      integer bnd_idx_bga
      integer bnd_idx_dbg
      integer cld_lvl_nbr
      integer slr_spc_xtr_typ
      integer xtr_typ_LHS
      integer xtr_typ_RHS

      real odxl_tmp             ! [frc] Temporary optical depth
      real mpl_tmp              ! [kg m-2] Temporary optical depth
      real alb_cmd_ln
      real wvl_dbg
      real float_foo
      real flx_spc_act
      real flx_spc_act_pht
      real flx_spc_net(bnd_nbr_max,levp_nbr_max)
      real idx_rfr_air_STP(bnd_nbr_max)
      real j_spc_NO2
      real mpc_CWP_cmd_ln
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
      real odal_N2O(lev_nbr_max)
      real odal_aer(lev_nbr_max)
      real odal_bga(lev_nbr_max)
      real odal_ice(lev_nbr_max)
      real odal_lqd(lev_nbr_max)
      real odal_spc_ttl(bnd_nbr_max,lev_nbr_max)
      real odsl_Ray(lev_nbr_max)
      real odsl_aer(lev_nbr_max)
      real odsl_bga(lev_nbr_max)
      real odsl_ice(lev_nbr_max)
      real odsl_lqd(lev_nbr_max)
      real odsl_spc_ttl(bnd_nbr_max,lev_nbr_max)
      real odxc_obs_aer_cmd_ln
      real sca_cff_mss_Ray(lev_nbr_max)
      real sca_frc_ttl(lev_nbr_max)
      real slr_cst
      real slr_cst_cmd_ln
      real slr_cst_xnt_fac
      real spc_npt_gnr(bnd_nbr_Mlk_max)
      real tpt_dlt_Mlk(lev_nbr_max)
      real tpt_dlt_Mlk_sqr(lev_nbr_max)
      real wvl_Planck

c     The following line is from the DISORT() subroutine. 
c     parameter ( mxcly = 85, mxulv = 85, mxcmu = 49, mxumu = 16,
c     $     maxpphi = 3)
c     The parameters declared below should match the values from
c     the DISORT() subroutine.
      integer maxcly
      integer maxcmu
      integer maxphi
      integer maxulv
      integer maxumu
      parameter(maxcly=tau_nbr_max,
     $     maxcmu = str_nbr_max,
     $     maxphi = azi_nbr_max,
     $     maxulv = tau_nbr_max,
     $     maxumu = plr_nbr_max)

c     DISORT() input variables:
      character  header*127
      logical  deltam, lamber, plank, onlyfl, prnt(7), usrang, usrtau
      integer  ibcnd, nlyr, numu, nstr, nphi, ntau
      real     accur, albedo, btemp, dtauc( maxcly ), fbeam, fisot,
     $     hl( 0:maxcmu ), phi( maxphi ), pmom( 0:maxcmu, maxcly ),
     $     phi0, ssalb( maxcly ), temper( 0:maxcly ), temis, ttemp,
     $     wvnmlo, wvnmhi, umu( maxumu ), umu0, utau( maxulv )

c     DISORT() output variables:
      real     rfldir( maxulv ), rfldn( maxulv ), flup( maxulv ),
     $     dfdt( maxulv ), uavg( maxulv ), u0u( maxumu, maxulv ),
     $     uu( maxumu, maxulv, maxphi ), albmed( maxumu ),
     $     trnmed( maxumu )

c     Externals
      integer iargc
      integer ftn_strlen            ! Return length of string (libcsz_f77)
      integer nf90_wrp_close      ! Close netCDF file (libcsz_f77)
      integer nf90_wrp_open       ! Open existing netCDF file (libcsz_f77)
      integer vec_val2idx       ! Locate index of array member closest to specified value (libcsz_f77)
      logical mnt_ncr_chk       ! Returns true for monotonic increasing arrays (libcsz_f77)
c      real abs_xsx_H2OH2O_Chy97 ! m2 mlc-1

c     Main code
      dbg_lvl=0                 ! Causes DDD source window to display this file

c     Initialize default values
#if ( defined RS6K )
      fl_CO2='/fs/cgd/data0/zender/aca/mlk_CO2.nc'
      fl_H2OH2O='/fs/cgd/data0/zender/aca/abs_xsx_H2OH2O.nc'
      fl_H2O='/fs/cgd/data0/zender/aca/mlk_H2O.nc'
      fl_OH='/fs/cgd/data0/zender/aca/mlk_OH.nc'
      fl_CH4='/fs/cgd/data0/zender/aca/mlk_CH4.nc'
      fl_N2O='/fs/cgd/data0/zender/aca/mlk_N2O.nc'
      fl_O2='/fs/cgd/data0/zender/aca/mlk_O2.nc'
      fl_O3='/fs/cgd/data0/zender/aca/abs_xsx_O3.nc'
      fl_O2O2='/fs/cgd/data0/zender/aca/abs_xsx_O2O2.nc'
      fl_NO2='/fs/cgd/data0/zender/aca/abs_xsx_NO2.nc'
      fl_clm='/fs/cgd/data0/zender/aca/mls_cld.nc'
      fl_ice='/fs/cgd/data0/zender/aca/ice_20.nc'
      fl_lqd='/fs/cgd/data0/zender/aca/lqd_10.nc'
      fl_nst='/fs/cgd/data0/zender/aca/nst_FSBR.nc'
      fl_aer='/fs/cgd/data0/zender/aca/aer_mineral_dust.nc'
      fl_bga='/fs/cgd/data0/zender/aca/aer_h2so4_215K.nc'
      fl_out='/fs/cgd/data0/zender/aca/nbm.nc'
      fl_slr='/fs/cgd/data0/zender/aca/spc_Kur95_01wvn.nc'
#else /* not RS6K */ 
      fl_CO2='/data/zender/aca/mlk_CO2.nc'
      fl_H2OH2O='/data/zender/aca/abs_xsx_H2OH2O.nc'
      fl_H2O='/data/zender/aca/mlk_H2O.nc'
      fl_OH='/data/zender/aca/mlk_OH.nc'
      fl_CH4='/data/zender/aca/mlk_CH4.nc'
      fl_N2O='/data/zender/aca/mlk_N2O.nc'
      fl_O2='/data/zender/aca/mlk_O2.nc'
      fl_O3='/data/zender/aca/abs_xsx_O3.nc'
      fl_O2O2='/data/zender/aca/abs_xsx_O2O2.nc'
      fl_NO2='/data/zender/aca/abs_xsx_NO2.nc'
      fl_clm='/data/zender/aca/mls_cld.nc'
      fl_ice='/data/zender/aca/ice_20.nc'
      fl_lqd='/data/zender/aca/lqd_10.nc'
      fl_nst='/data/zender/aca/nst_FSBR.nc'
      fl_aer='/data/zender/aca/aer_mineral_dust.nc'
      fl_bga='/data/zender/aca/aer_h2so4_215K.nc'
      fl_out='/data/zender/aca/nbm.nc'
      fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
#endif /* not RS6K */ 
      CVS_Date='$Date$'
      CVS_Revision='$Revision$'
      azi_nbr=azi_nbr_max
      wvl_dbg=0.5e-6            ! [um] Wavelength of debugging band
      cmd_ln_alb=.false.
      cmd_ln_mpc_CWP=.false.
      cmd_ln_odxc_obs_aer=.false.
      cmd_ln_slr_cst=.false.
      cmd_ln_slr_zen_ngl_cos=.false.
      dbg_lvl=0
      exit_status=0
      flg_CO2=.true.
      flg_H2O=.true.
      flg_H2OH2O=.false.
      flg_Herzberg=.true.
      flg_NO2=.true.
      flg_O2=.true.
      flg_O2N2=.false.
      flg_O2O2=.false.
      flg_O3=.true.
      flg_OH=.true.
      flg_CH4=.true.
      flg_N2O=.true.
      flg_Planck=.true.
      flg_Rayleigh=.true.
      flg_aer=.true.
      flg_bga=.true.
      flg_ice=.true.
      flg_lqd=.true.
      float_foo=0.
      flt_nst=.true.
      force_ice_phz=.false.
      force_lqd_phz=.false.
      logical_foo=.true.
      pi=4.0*atan(1.0)          ! 3
      plr_nbr=2
      rcd=nf90_noerr              ! nf90_noerr == 0
      wvl_grd_rsn=5.0e-9        ! [m] 
      single_bnd_computation=.false.
      slr_cst=slr_cst_CCM
      str_nbr=4
      sv_cmp_tau=.true.
      sv_cmp_plr_ngl=.true.
      sv_ntn=.false.
      top_lvl=.false.
      tst_case_HG=.false.
      tst_case_Rayleigh=.false.
      wvl_Planck=2.0e-6
      wvl_grd_max=5.0e-6        ! [m] H2O Malkmus data ends at 5.0 um
      wvl_grd_min=176.0e-9      ! [m] WMO O2 & O3 data starts at 175.4 nm
      wvl_grd_typ=wvl_grd_dfl   ! Wavelength grid type

c     Retrieve command line arguments
      call date_time_get(lcl_date_time)
      call ftn_cmd_ln_sng(cmd_ln)
      call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
      write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
#ifdef LINUX
         write (6,'(5a)') prg_nm(1:ftn_strlen(prg_nm)),
     $        ': WARNING 32-bit Intel Linux may cause SIGFPE and dump core when compiled with -Ktrap=fp with pgf90.',
     $        '  64 bit machines do not seem to suffer from this problem.',
     $        '  Trigger for SIGFPE is invalid operation, divide by zero, or overflow',
     $        '  See file:/usr/local/pgi/doc/pgiws_ug/pgi31u02.htm#Heading6'
#endif /* endif LINUX */
      arg_nbr=iargc()
      do arg=1,arg_nbr
         call getarg(arg,argv)
         if (argv(1:2) == '-A') then
            flg_aer=.not.flg_aer
         endif
         if (argv(1:2) == '-a') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_aer
         endif
         if (argv(1:2) == '-B') then
            flg_bga=.not.flg_bga
         endif
         if (argv(1:2) == '-b') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_bga
         endif
         if (argv(1:2) == '-C') then
            flg_CO2=.not.flg_CO2
         endif
         if (argv(1:2) == '-c') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_CO2
         endif
         if (argv(1:2) == '-D') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') dbg_lvl
         endif
         if (argv(1:2) == '-d') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_out
         endif
         if (argv(1:2) == '-E') then
            single_bnd_computation=.not.single_bnd_computation
         endif
         if (argv(1:2) == '-e') then
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') wvl_dbg ! [m]
         endif
         if (argv(1:2) == '-F') then
            force_ice_phz=.not.force_ice_phz
         endif
         if (argv(1:2) == '-f') then
            force_lqd_phz=.not.force_lqd_phz
         endif
         if (argv(1:2) == '-H') then
            flg_H2O=.not.flg_H2O
         endif
         if (argv(1:2) == '-h') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_H2O
         endif
         if (argv(1:2) == '-I') then
            flg_ice=.not.flg_ice
         endif
         if (argv(1:2) == '-i') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_ice
         endif
         if (argv(1:2) == '-J') then
            flg_Planck=.not.flg_Planck
         endif
         if (argv(1:2) == '-j') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') wvl_grd_typ
         endif
         if (argv(1:2) == '-K') then
            flg_O2O2=.not.flg_O2O2
         endif
         if (argv(1:2) == '-k') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O2O2
         endif
         if (argv(1:2) == '-L') then
            flg_lqd=.not.flg_lqd
         endif
         if (argv(1:2) == '-l') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_lqd
         endif
         if (argv(1:2) == '-M') then
            cmd_ln_odxc_obs_aer=.not.cmd_ln_odxc_obs_aer
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') odxc_obs_aer_cmd_ln
         endif
         if (argv(1:2) == '-m') then
            cmd_ln_mpc_CWP=.not.cmd_ln_mpc_CWP
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') mpc_CWP_cmd_ln
         endif
         if (argv(1:2) == '-N') then
            flt_nst=.not.flt_nst
         endif
         if (argv(1:2) == '-n') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_nst
         endif
         if (argv(1:2) == '-O') then
            flg_O2=.not.flg_O2
         endif
         if (argv(1:2) == '-o') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O2
         endif
         if (argv(1:2) == '-P') then
            tst_case_HG=.not.tst_case_HG
         endif
         if (argv(1:2) == '-p') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_clm
         endif
         if (argv(1:2) == '-Q') then
            flg_H2OH2O=.not.flg_H2OH2O
         endif
         if (argv(1:2) == '-q') then
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') wvl_grd_rsn ! [m]
         endif
         if (argv(1:2) == '-R') then
            flg_Rayleigh=.not.flg_Rayleigh
         endif
         if (argv(1:2) == '-r') then
            cmd_ln_alb=.not.cmd_ln_alb
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') alb_cmd_ln
         endif
         if (argv(1:2) == '-S') then
            cmd_ln_slr_cst=.not.cmd_ln_slr_cst
            call getarg(arg+1,argv)
            read (argv,'(f8.3)') slr_cst_cmd_ln
         endif
         if (argv(1:2) == '-s') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') str_nbr
         endif
         if (argv(1:2) == '-T') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_slr
         endif
         if (argv(1:2) == '-t') then
            tst_case_Rayleigh=.not.tst_case_Rayleigh
         endif
         if (argv(1:2) == '-U') then
            flg_O2N2=.not.flg_O2N2
         endif
         if (argv(1:2) == '-u') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') plr_nbr
         endif
         if (argv(1:2) == '-V') then
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') wvl_grd_max ! [m]
         endif
         if (argv(1:2) == '-v') then
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') wvl_grd_min ! [m]
         endif
         if (argv(1:2) == '-W') then
            flg_O3=.not.flg_O3
         endif
         if (argv(1:2) == '-w') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_O3
         endif
         if (argv(1:2) == '-X') then
            flg_NO2=.not.flg_NO2
         endif
         if (argv(1:2) == '-x') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_NO2
         endif
         if (argv(1:2) == '-Y') then
            flg_OH=.not.flg_OH
            flg_CH4=.not.flg_CH4
            flg_N2O=.not.flg_N2O
         endif
         if (argv(1:2) == '-y') then
            call getarg(arg+1,argv)
            read (argv,'(a)') fl_OH
         endif
         if (argv(1:2) == '-Z') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') azi_nbr
         endif
         if (argv(1:2) == '-z') then
            cmd_ln_slr_zen_ngl_cos=.not.cmd_ln_slr_zen_ngl_cos
            call getarg(arg+1,argv)
            read (argv,'(f15.8)') slr_zen_ngl_cos_cmd_ln
         endif
      end do

c     Compute any quantities that might depend on command line input
      call ftn_strnul(aer_sng)
      call ftn_strnul(azi_sng)
      call ftn_strnul(bga_sng)
      call ftn_strnul(fl_CO2)
      call ftn_strnul(fl_H2O)
      call ftn_strnul(fl_H2OH2O)
      call ftn_strnul(fl_NO2)
      call ftn_strnul(fl_O2)
      call ftn_strnul(fl_O2O2)
      call ftn_strnul(fl_O3)
      call ftn_strnul(fl_OH)
      call ftn_strnul(fl_CH4)
      call ftn_strnul(fl_N2O)
      call ftn_strnul(fl_aer)
      call ftn_strnul(fl_bga)
      call ftn_strnul(fl_clm)
      call ftn_strnul(fl_ice)
      call ftn_strnul(fl_lqd)
      call ftn_strnul(fl_nst)
      call ftn_strnul(fl_out)
      call ftn_strnul(fl_slr)
      call ftn_strnul(opt_dep_sng)
      call ftn_strnul(plr_sng)
      call ftn_strnul(prf_sng)
      call ftn_strnul(str_sng)
      call ftn_strnul(src_rfr_sng)
      call ftn_strcpy(src_rfr_sng,'Model reference is Zender et al. (1997) (ZBP97)')

      flg_Herzberg=flg_O3
      nbr_mom=str_nbr           ! # moments always equals # streams
      if (sv_cmp_plr_ngl) then
         plr_nbr=str_nbr
      endif                     ! endif

c     Define wavelength grid
      call wvl_grd_mk(wvl_grd_typ, ! Input
     $     wvl_grd_min,wvl_grd_max,wvl_grd_rsn, ! Input/Output
     $     wvl_ctr,wvl_dlt,wvl_min,wvl_max,wvl_grd,wvl_nbr) ! Output
      bnd_nbr=wvl_nbr
c     Enough memory?
      if (bnd_nbr > bnd_nbr_max) stop 'bnd_nbr > bnd_nbr_max in nbm()'
c     Update fields which depend on wavelength grid
      if (wvl_dbg < wvl_grd_min.or.wvl_dbg > wvl_grd_max) wvl_dbg=0.5*(wvl_grd_min+wvl_grd_max)
      bnd_idx_dbg=vec_val2idx(wvl_ctr,wvl_nbr,wvl_dbg) ! Locate debugging wavelength bin

c     Derived spectral grid quantities
      do bnd_idx=1,bnd_nbr
         bnd(bnd_idx)=wvl_ctr(bnd_idx) ! [m]
         wvl(bnd_idx)=wvl_ctr(bnd_idx) ! [m]
         wvn_max(bnd_idx)=1.0/(100.0*wvl_min(bnd_idx)) ! [cm-1]
         wvn_min(bnd_idx)=1.0/(100.0*wvl_max(bnd_idx)) ! [cm-1]
         wvn_ctr(bnd_idx)=0.5*(wvn_min(bnd_idx)+wvn_max(bnd_idx)) ! [cm-1]
         wvn(bnd_idx)=wvn_ctr(bnd_idx) ! [cm-1]
         wvn_dlt(bnd_idx)=wvn_max(bnd_idx)-wvn_min(bnd_idx) ! [cm-1]
      enddo                     ! end loop over bnd

c     Wavenumber grid depends on direction of monotonicity of wavelength grid
      mnt_ncr=mnt_ncr_chk(wvl_grd,wvl_nbr+1)
      if (mnt_ncr) then
         do bnd_idx=1,bnd_nbr
            wvn_grd(bnd_idx)=wvn_max(bnd_idx)
         enddo                  ! end loop over bnd
         wvn_grd(bnd_nbr+1)=wvn_min(bnd_nbr)
      else                      ! not mnt_ncr
         do bnd_idx=1,bnd_nbr
            wvn_grd(bnd_idx)=wvn_min(bnd_idx)
         enddo                  ! end loop over bnd
         wvn_grd(bnd_nbr+1)=wvn_max(bnd_nbr)
      endif                     ! not mnt_ncr

c     Define initialization state
      if (flg_CO2) then
         call ftn_strcpy(stt_CO2,'CO2 line absorption: Malkmus narrow band parameters from '//fl_CO2) 
      else
         call ftn_strcpy(stt_CO2,'CO2 line absorption: Off')
      endif                     ! not
      if (flg_H2O) then
         call ftn_strcpy(stt_H2O,'H2O line absorption: Malkmus narrow band parameters from '//fl_H2O)
      else
         call ftn_strcpy(stt_H2O,'H2O line absorption: Off')
      endif                     ! not
      if (flg_H2OH2O) then
         call ftn_strcpy(stt_H2OH2O,'H2O-H2O dimer absorption: Continuum absorption cross sections from '//fl_H2OH2O)
      else
         call ftn_strcpy(stt_H2OH2O,'H2O-H2O dimer absorption: Off')
      endif
      if (flg_OH) then
         call ftn_strcpy(stt_OH,'OH line absorption: Malkmus narrow band parameters from '//fl_OH)
      else
         call ftn_strcpy(stt_OH,'OH line absorption: Off')
      endif
      if (flg_CH4) then
         call ftn_strcpy(stt_CH4,'CH4 line absorption: Malkmus narrow band parameters from '//fl_CH4)
      else
         call ftn_strcpy(stt_CH4,'CH4 line absorption: Off')
      endif
      if (flg_N2O) then
         call ftn_strcpy(stt_N2O,'N2O line absorption: Malkmus narrow band parameters from '//fl_N2O)
      else
         call ftn_strcpy(stt_N2O,'N2O line absorption: Off')
      endif
      if (flg_O2) then
         call ftn_strcpy(stt_O2,'O2 line absorption: Malkmus narrow band parameters from '//fl_O2)
      else
         call ftn_strcpy(stt_O2,'O2 line absorption: Off')
      endif
      if (flg_Herzberg) then
         call ftn_strcpy(stt_Herzberg,'O2 Herzberg bands: Continuum absorption cross sections from '//fl_O3)
      else
         call ftn_strcpy(stt_Herzberg,'O2 Herzberg bands: Off')
      endif
      if (flg_O3) then
         call ftn_strcpy(stt_O3,'O3 Hartley, Huggins, and Chappuis bands: Continuum absorption cross sections from '//fl_O3)
      else
         call ftn_strcpy(stt_O3,'O3 Hartley, Huggins, and Chappuis bands: Off')
      endif
      if (flg_O2O2) then
         call ftn_strcpy(stt_O2O2,'O2-O2 collision-induced absorption: Continuum absorption cross sections from '//fl_O2O2)
      else
         call ftn_strcpy(stt_O2O2,'O2-O2 collision-induced absorption: Off')
      endif
      if (flg_O2N2) then
         call ftn_strcpy(stt_O2N2,'O2-N2 collision-induced absorption: Occurs in 1.27 um band with 0.2 efficiency of O2-O2')
      else
         call ftn_strcpy(stt_O2N2,'O2-N2 collision-induced absorption: Off')
      endif
      if (flg_NO2) then
         call ftn_strcpy(stt_NO2,'NO2 absorption: Continuum absorption cross sections from '//fl_NO2)
      else
         call ftn_strcpy(stt_NO2,'NO2 absorption: Off')
      endif
      if (flg_ice) then
         call ftn_strcpy(stt_ice,'Ice water crystal scattering and absorption: Mie theory from '//fl_ice)
      else
         call ftn_strcpy(stt_ice,'Ice water crystal scattering and absorption: Off')
      endif
      if (flg_lqd) then
         call ftn_strcpy(stt_lqd,'Liquid water droplet scattering and absorption: Mie theory from '//fl_lqd)
      else
         call ftn_strcpy(stt_lqd,'Liquid water droplet scattering and absorption: Off')
      endif
      if (flg_Planck) then
         call ftn_strcpy(stt_Planck,'Thermal emission of atmosphere: Included for lambda > 2 um')
      else
         call ftn_strcpy(stt_Planck,'Thermal emission of atmosphere: Off')
      endif
      if (flg_Rayleigh) then
         call ftn_strcpy(stt_Rayleigh,'Rayleigh scattering: method of Lenoble (1993)')
      else
         call ftn_strcpy(stt_Rayleigh,'Rayleigh scattering: Off')
      endif
      if (flg_aer) then
         call ftn_strcpy(stt_aer,'Aerosol scattering and absorption: Mie theory from '//fl_aer)
      else
         call ftn_strcpy(stt_aer,'Aerosol scattering and absorption: Off')
      endif
      if (flg_bga) then
         call ftn_strcpy(stt_bga,'Background aerosol scattering and absorption: Mie theory from '//fl_bga)
      else
         call ftn_strcpy(stt_bga,'Background aerosol scattering and absorption: Off')
      endif
      if (flt_nst) then
         call ftn_strcpy(stt_flt_nst,'Instrument filter spectral response function: instrument from '//fl_nst)
      else
         call ftn_strcpy(stt_flt_nst,'Instrument filter spectral response function: Off')
      endif
      call ftn_strcpy(stt_slr,'TOA solar spectrum from '//fl_slr)
      if (top_lvl) then
         call ftn_strcpy(stt_top_lvl,'Pure absorbing layer to space at model top not yet implemented')
      else
         call ftn_strcpy(stt_top_lvl,'Pure absorbing layer to space at model top not yet implemented')
      endif                     ! not

c     Ingest fl_clm
      rcd=rcd+nf90_wrp_open(fl_clm,nf90_nowrite,nc_id)
c     Get global attributes
      rcd=rcd+nf90_get_att(nc_id,nf90_global,'prf_sng',prf_sng)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'lev',lev_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'levp',levp_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,lev_dim_id,len=lev_nbr)
      if (lev_nbr > lev_nbr_max) stop 'lev_nbr > lev_nbr_max'
      rcd=rcd+nf90_inquire_dimension(nc_id,levp_dim_id,len=levp_nbr)
      if (levp_nbr > levp_nbr_max) stop 'levp_nbr > levp_nbr_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'RH_lqd',RH_lqd_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alb_sfc_NIR_drc',alb_sfc_NIR_drc_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alb_sfc_NIR_dff',alb_sfc_NIR_dff_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alb_sfc_vsb_drc',alb_sfc_vsb_drc_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alb_sfc_vsb_dff',alb_sfc_vsb_dff_id)
      rcd=rcd+nf90_inq_varid(nc_id,'slr_zen_ngl_cos',slr_zen_ngl_cos_id)
      rcd=rcd+nf90_inq_varid(nc_id,'prs_dlt',prs_dlt_id)
      rcd=rcd+nf90_inq_varid(nc_id,'frc_ice',frc_ice_id)
      rcd=rcd+nf90_inq_varid(nc_id,'frc_ice_ttl',frc_ice_ttl_id)
      rcd=rcd+nf90_inq_varid(nc_id,'grv',grv_id)
      rcd=rcd+nf90_inq_varid(nc_id,'lat_dgr',lat_dgr_id)
      rcd=rcd+nf90_inq_varid(nc_id,'lcl_time_hr',lcl_time_hr_id)
      rcd=rcd+nf90_inq_varid(nc_id,'lcl_yr_day',lcl_yr_day_id)
      rcd=rcd+nf90_inq_varid(nc_id,'lev',lev_id)
      rcd=rcd+nf90_inq_varid(nc_id,'levp',levp_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mmw_mst_air',mmw_mst_air_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpc_CWP',mpc_CWP_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_CO2',mpl_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_CWP',mpl_CWP_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_H2O',mpl_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_IWP',mpl_IWP_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_LWP',mpl_LWP_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_O2',mpl_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_OH',mpl_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_CH4',mpl_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_N2O',mpl_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_aer',mpl_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_bga',mpl_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'mpl_mst_air',mpl_mst_air_id)
      rcd=rcd+nf90_inq_varid(nc_id,'npl_NO2',npl_NO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'npl_O2',npl_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'npl_O3',npl_O3_id)
      rcd=rcd+nf90_inq_varid(nc_id,'npl_O2O2',npl_O2O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'npl_H2OH2O',npl_H2OH2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'odxc_obs_aer',odxc_obs_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'odxc_obs_bga',odxc_obs_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'odxl_obs_aer',odxl_obs_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'odxl_obs_bga',odxl_obs_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'prs',prs_id)
      rcd=rcd+nf90_inq_varid(nc_id,'prs_ntf',prs_ntf_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_CO2',q_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_H2O',q_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_O2',q_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_OH',q_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_CH4',q_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'q_N2O',q_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'spc_heat_mst_air',spc_heat_mst_air_id)
      rcd=rcd+nf90_inq_varid(nc_id,'tpt',tpt_id)
      rcd=rcd+nf90_inq_varid(nc_id,'tpt_ntf',tpt_ntf_id)
      rcd=rcd+nf90_inq_varid(nc_id,'tpt_skn',tpt_skn_id)
      rcd=rcd+nf90_inq_varid(nc_id,'xnt_fac',xnt_fac_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_obs_aer',wvl_obs_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_obs_bga',wvl_obs_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alt_cld_btm',alt_cld_btm_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alt_cld_thick',alt_cld_thick_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alt',alt_id)
      rcd=rcd+nf90_inq_varid(nc_id,'alt_ntf',alt_ntf_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,lat_dgr_id,lat_dgr)
      rcd=rcd+nf90_get_var(nc_id,lcl_time_hr_id,lcl_time_hr)
      rcd=rcd+nf90_get_var(nc_id,lcl_yr_day_id,lcl_yr_day)
      rcd=rcd+nf90_get_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos)
      rcd=rcd+nf90_get_var(nc_id,RH_lqd_id,RH_lqd)
      rcd=rcd+nf90_get_var(nc_id,alb_sfc_NIR_drc_id,alb_sfc_NIR_drc)
      rcd=rcd+nf90_get_var(nc_id,alb_sfc_NIR_dff_id,alb_sfc_NIR_dff)
      rcd=rcd+nf90_get_var(nc_id,alb_sfc_vsb_drc_id,alb_sfc_vsb_drc)
      rcd=rcd+nf90_get_var(nc_id,alb_sfc_vsb_dff_id,alb_sfc_vsb_dff)
      rcd=rcd+nf90_get_var(nc_id,prs_dlt_id,prs_dlt)
      rcd=rcd+nf90_get_var(nc_id,frc_ice_id,frc_ice)
      rcd=rcd+nf90_get_var(nc_id,frc_ice_ttl_id,frc_ice_ttl)
      rcd=rcd+nf90_get_var(nc_id,grv_id,grv)
      rcd=rcd+nf90_get_var(nc_id,lev_id,lev)
      rcd=rcd+nf90_get_var(nc_id,levp_id,levp)
      rcd=rcd+nf90_get_var(nc_id,mmw_mst_air_id,mmw_mst_air)
      rcd=rcd+nf90_get_var(nc_id,mpc_CWP_id,mpc_CWP)
      rcd=rcd+nf90_get_var(nc_id,mpl_CO2_id,mpl_CO2)
      rcd=rcd+nf90_get_var(nc_id,mpl_CWP_id,mpl_CWP)
      rcd=rcd+nf90_get_var(nc_id,mpl_H2O_id,mpl_H2O)
      rcd=rcd+nf90_get_var(nc_id,mpl_IWP_id,mpl_IWP)
      rcd=rcd+nf90_get_var(nc_id,mpl_LWP_id,mpl_LWP)
      rcd=rcd+nf90_get_var(nc_id,mpl_O2_id,mpl_O2)
      rcd=rcd+nf90_get_var(nc_id,mpl_OH_id,mpl_OH)
      rcd=rcd+nf90_get_var(nc_id,mpl_CH4_id,mpl_CH4)
      rcd=rcd+nf90_get_var(nc_id,mpl_N2O_id,mpl_N2O)
      rcd=rcd+nf90_get_var(nc_id,mpl_aer_id,mpl_aer)
      rcd=rcd+nf90_get_var(nc_id,mpl_bga_id,mpl_bga)
      rcd=rcd+nf90_get_var(nc_id,mpl_mst_air_id,mpl_mst_air)
      rcd=rcd+nf90_get_var(nc_id,npl_NO2_id,npl_NO2)
      rcd=rcd+nf90_get_var(nc_id,npl_O2O2_id,npl_O2O2)
      rcd=rcd+nf90_get_var(nc_id,npl_O2_id,npl_O2)
      rcd=rcd+nf90_get_var(nc_id,npl_O3_id,npl_O3)
      rcd=rcd+nf90_get_var(nc_id,npl_H2OH2O_id,npl_H2OH2O)
      rcd=rcd+nf90_get_var(nc_id,odxc_obs_aer_id,odxc_obs_aer)
      rcd=rcd+nf90_get_var(nc_id,odxc_obs_bga_id,odxc_obs_bga)
      rcd=rcd+nf90_get_var(nc_id,odxl_obs_aer_id,odxl_obs_aer)
      rcd=rcd+nf90_get_var(nc_id,odxl_obs_bga_id,odxl_obs_bga)
      rcd=rcd+nf90_get_var(nc_id,prs_id,prs)
      rcd=rcd+nf90_get_var(nc_id,prs_ntf_id,prs_ntf)
      rcd=rcd+nf90_get_var(nc_id,q_CO2_id,q_CO2)
      rcd=rcd+nf90_get_var(nc_id,q_H2O_id,q_H2O)
      rcd=rcd+nf90_get_var(nc_id,q_O2_id,q_O2)
      rcd=rcd+nf90_get_var(nc_id,q_OH_id,q_OH)
      rcd=rcd+nf90_get_var(nc_id,q_CH4_id,q_CH4)
      rcd=rcd+nf90_get_var(nc_id,q_N2O_id,q_N2O)
      rcd=rcd+nf90_get_var(nc_id,spc_heat_mst_air_id,spc_heat_mst_air)
      rcd=rcd+nf90_get_var(nc_id,tpt_id,tpt)
      rcd=rcd+nf90_get_var(nc_id,tpt_ntf_id,tpt_ntf)
      rcd=rcd+nf90_get_var(nc_id,tpt_skn_id,tpt_skn)
      rcd=rcd+nf90_get_var(nc_id,wvl_obs_aer_id,wvl_obs_aer)
      rcd=rcd+nf90_get_var(nc_id,wvl_obs_bga_id,wvl_obs_bga)
      rcd=rcd+nf90_get_var(nc_id,xnt_fac_id,xnt_fac)
      rcd=rcd+nf90_get_var(nc_id,alt_cld_btm_id,alt_cld_btm)
      rcd=rcd+nf90_get_var(nc_id,alt_cld_thick_id,alt_cld_thick)
      rcd=rcd+nf90_get_var(nc_id,alt_id,alt)
      rcd=rcd+nf90_get_var(nc_id,alt_ntf_id,alt_ntf)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_clm,'Ingested')

c     Ingest fl_H2O
      rcd=rcd+nf90_wrp_open(fl_H2O,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_H2O)
      if (bnd_nbr_H2O > bnd_nbr_H2O_max) stop 'bnd_nbr_H2O > bnd_nbr_H2O_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_H2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_H2O_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_H2O_id,wvl_grd_H2O)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_H2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_H2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_H2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_H2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_H2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_H2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_H2O,wvl_grd_H2O,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_H2O,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_H2O,'Ingested')

c     Ingest fl_OH
      rcd=rcd+nf90_wrp_open(fl_OH,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_OH)
      if (bnd_nbr_OH > bnd_nbr_OH_max) stop 'bnd_nbr_OH > bnd_nbr_OH_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_OH_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_OH_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_OH_id,wvl_grd_OH)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_OH,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_OH,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_OH,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_OH,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_OH,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_OH_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_OH,wvl_grd_OH,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_OH,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_OH,'Ingested')

c     Ingest fl_CH4
      rcd=rcd+nf90_wrp_open(fl_CH4,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_CH4)
      if (bnd_nbr_CH4 > bnd_nbr_CH4_max) stop 'bnd_nbr_CH4 > bnd_nbr_CH4_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_CH4_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_CH4_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_CH4_id,wvl_grd_CH4)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_CH4,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_CH4,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_CH4,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_CH4,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_CH4,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_CH4_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CH4,wvl_grd_CH4,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_CH4,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_CH4,'Ingested')

c     Ingest fl_N2O
      rcd=rcd+nf90_wrp_open(fl_N2O,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_N2O)
      if (bnd_nbr_N2O > bnd_nbr_N2O_max) stop 'bnd_nbr_N2O > bnd_nbr_N2O_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_N2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_N2O_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_N2O_id,wvl_grd_N2O)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_N2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_N2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_N2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_N2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_N2O,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_N2O_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_N2O,wvl_grd_N2O,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_N2O,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_N2O,'Ingested')

c     Ingest fl_CO2
      rcd=rcd+nf90_wrp_open(fl_CO2,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_CO2)
      if (bnd_nbr_CO2 > bnd_nbr_CO2_max) stop 'bnd_nbr_CO2 > bnd_nbr_CO2_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_CO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_CO2_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_CO2_id,wvl_grd_CO2)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_CO2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_CO2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_CO2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_CO2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_CO2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_CO2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_CO2,wvl_grd_CO2,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_CO2,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_CO2,'Ingested')

c     Ingest fl_O2
      rcd=rcd+nf90_wrp_open(fl_O2,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_O2)
      if (bnd_nbr_O2 > bnd_nbr_O2_max) stop 'bnd_nbr_O2 > bnd_nbr_O2_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'A_phi',A_phi_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'A_psi',A_psi_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_phi',B_phi_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'B_psi',B_psi_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_O2_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_O2_id,wvl_grd_O2)
c     Malkmus parameters should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      rcd=rcd+nf90_get_var(nc_id,A_phi_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,A_phi_O2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,A_psi_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,A_psi_O2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_phi_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,B_phi_O2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,B_psi_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,B_psi_O2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_d_abs_cff_mss_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,S_d_abs_cff_mss_O2,xtr_typ_LHS,xtr_typ_RHS)
      rcd=rcd+nf90_get_var(nc_id,S_p_abs_cff_mss_O2_id,spc_npt_gnr)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,spc_npt_gnr,bnd_nbr,wvl_grd,S_p_abs_cff_mss_O2,xtr_typ_LHS,xtr_typ_RHS)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_O2,'Ingested')

c     Ingest fl_O3
      rcd=rcd+nf90_wrp_open(fl_O3,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_O3)
      if (bnd_nbr_O3 > bnd_nbr_O3_max) stop 'bnd_nbr_O3 > bnd_nbr_O3_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_xsx_O2',abs_xsx_O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'abs_xsx_O3_cold',abs_xsx_O3_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_O3_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_xsx_O2_id,abs_xsx_O2_obs)
      rcd=rcd+nf90_get_var(nc_id,abs_xsx_O3_id,abs_xsx_O3_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_O3_id,wvl_grd_O3)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_O3,'Ingested')
c     Absorption cross-sections should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      call rbn_vec(bnd_nbr_O3,wvl_grd_O3,abs_xsx_O3_obs,bnd_nbr,wvl_grd,abs_xsx_O3,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_O2,wvl_grd_O2,abs_xsx_O2_obs,bnd_nbr,wvl_grd,abs_xsx_O2,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_O2O2
      rcd=rcd+nf90_wrp_open(fl_O2O2,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_O2O2)
      if (bnd_nbr_O2O2 > bnd_nbr_O2O2_max) stop 'bnd_nbr_O2O2 > bnd_nbr_O2O2_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_xsx_O2O2',abs_xsx_O2O2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_O2O2_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_O2O2_id,wvl_grd_O2O2)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_O2O2,'Ingested')
c     Absorption cross-sections should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      call rbn_vec(bnd_nbr_O2O2,wvl_grd_O2O2,abs_xsx_O2O2_obs,bnd_nbr,wvl_grd,abs_xsx_O2O2,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_NO2
      rcd=rcd+nf90_wrp_open(fl_NO2,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_NO2)
      if (bnd_nbr_NO2 > bnd_nbr_NO2_max) stop 'bnd_nbr_NO2 > bnd_nbr_NO2_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_xsx_NO2',abs_xsx_NO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'qnt_yld_NO2',qnt_yld_NO2_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_NO2_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_xsx_NO2_id,abs_xsx_NO2_obs)
      rcd=rcd+nf90_get_var(nc_id,qnt_yld_NO2_id,qnt_yld_NO2_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_NO2_id,wvl_grd_NO2)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_NO2,'Ingested')
c     Quantum yields should be 1 at high energy and 0 at low energy
      xtr_typ_LHS=xtr_prt_ngh+xtr_fll_ngh
      xtr_typ_RHS=xtr_prt_ngh+xtr_fll_ngh
      call rbn_vec(bnd_nbr_NO2,wvl_grd_NO2,qnt_yld_NO2_obs,bnd_nbr,wvl_grd,qnt_yld_NO2,xtr_typ_LHS,xtr_typ_RHS)
c     Absorption cross-sections should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      call rbn_vec(bnd_nbr_NO2,wvl_grd_NO2,abs_xsx_NO2_obs,bnd_nbr,wvl_grd,abs_xsx_NO2,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_H2OH2O
      rcd=rcd+nf90_wrp_open(fl_H2OH2O,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_H2OH2O)
      if (bnd_nbr_H2OH2O > bnd_nbr_H2OH2O_max) stop 'bnd_nbr_H2OH2O > bnd_nbr_H2OH2O_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_xsx_H2OH2O',abs_xsx_H2OH2O_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_H2OH2O_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_xsx_H2OH2O_id,abs_xsx_H2OH2O_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_H2OH2O_id,wvl_grd_H2OH2O)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_H2OH2O,'Ingested')
c     Absorption cross-sections should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
      call rbn_vec(bnd_nbr_H2OH2O,wvl_grd_H2OH2O,abs_xsx_H2OH2O_obs,bnd_nbr,wvl_grd,abs_xsx_H2OH2O,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_ice
      rcd=rcd+nf90_wrp_open(fl_ice,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_ice)
      if (bnd_nbr_ice > bnd_nbr_ice_max) stop 'bnd_nbr_ice > bnd_nbr_ice_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_ice_id)
      rcd=rcd+nf90_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_ice_id)
      rcd=rcd+nf90_inq_varid(nc_id,'asm_prm',asm_prm_ice_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_ice_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_cff_mss_ice_id,abs_cff_mss_ice_obs)
      rcd=rcd+nf90_get_var(nc_id,sca_cff_mss_ice_id,sca_cff_mss_ice_obs)
      rcd=rcd+nf90_get_var(nc_id,asm_prm_ice_id,asm_prm_ice_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_ice_id,wvl_grd_ice)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_ice,'Ingested')
c     Absorption and scattering coefficients should be constant outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_ngh
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_ngh
      call rbn_vec(bnd_nbr_ice,wvl_grd_ice,abs_cff_mss_ice_obs,bnd_nbr,wvl_grd,abs_cff_mss_ice,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_ice,wvl_grd_ice,sca_cff_mss_ice_obs,bnd_nbr,wvl_grd,sca_cff_mss_ice,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_ice,wvl_grd_ice,asm_prm_ice_obs,bnd_nbr,wvl_grd,asm_prm_ice,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_lqd
      rcd=rcd+nf90_wrp_open(fl_lqd,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'bnd',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_lqd)
      if (bnd_nbr_lqd > bnd_nbr_lqd_max) stop 'bnd_nbr_lqd > bnd_nbr_lqd_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_lqd_id)
      rcd=rcd+nf90_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_lqd_id)
      rcd=rcd+nf90_inq_varid(nc_id,'asm_prm',asm_prm_lqd_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_lqd_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_cff_mss_lqd_id,abs_cff_mss_lqd_obs)
      rcd=rcd+nf90_get_var(nc_id,sca_cff_mss_lqd_id,sca_cff_mss_lqd_obs)
      rcd=rcd+nf90_get_var(nc_id,asm_prm_lqd_id,asm_prm_lqd_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_lqd_id,wvl_grd_lqd)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_lqd,'Ingested')
c     Absorption and scattering coefficients should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_ngh
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_ngh
      call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,abs_cff_mss_lqd_obs,bnd_nbr,wvl_grd,abs_cff_mss_lqd,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,sca_cff_mss_lqd_obs,bnd_nbr,wvl_grd,sca_cff_mss_lqd,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_lqd,wvl_grd_lqd,asm_prm_lqd_obs,bnd_nbr,wvl_grd,asm_prm_lqd,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_aer
      rcd=rcd+nf90_wrp_open(fl_aer,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'wvl',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_aer)
      if (bnd_nbr_aer > bnd_nbr_aer_max) stop 'bnd_nbr_aer > bnd_nbr_aer_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'ext_cff_mss',ext_cff_mss_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'asm_prm',asm_prm_aer_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_aer_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_cff_mss_aer_id,abs_cff_mss_aer_obs)
      rcd=rcd+nf90_get_var(nc_id,ext_cff_mss_aer_id,ext_cff_mss_aer_obs)
      rcd=rcd+nf90_get_var(nc_id,sca_cff_mss_aer_id,sca_cff_mss_aer_obs)
      rcd=rcd+nf90_get_var(nc_id,asm_prm_aer_id,asm_prm_aer_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_aer_id,wvl_grd_aer)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_aer,'Ingested')
c     Absorption and scattering coefficients should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_ngh
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_ngh
      call rbn_vec(bnd_nbr_aer,wvl_grd_aer,abs_cff_mss_aer_obs,bnd_nbr,wvl_grd,abs_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_aer,wvl_grd_aer,ext_cff_mss_aer_obs,bnd_nbr,wvl_grd,ext_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_aer,wvl_grd_aer,sca_cff_mss_aer_obs,bnd_nbr,wvl_grd,sca_cff_mss_aer,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_aer,wvl_grd_aer,asm_prm_aer_obs,bnd_nbr,wvl_grd,asm_prm_aer,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_bga
      rcd=rcd+nf90_wrp_open(fl_bga,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'wvl',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_bga)
      if (bnd_nbr_bga > bnd_nbr_bga_max) stop 'bnd_nbr_bga > bnd_nbr_bga_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'asm_prm',asm_prm_bga_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_bga_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,abs_cff_mss_bga_id,abs_cff_mss_bga_obs)
      rcd=rcd+nf90_get_var(nc_id,sca_cff_mss_bga_id,sca_cff_mss_bga_obs)
      rcd=rcd+nf90_get_var(nc_id,asm_prm_bga_id,asm_prm_bga_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_bga_id,wvl_grd_bga)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_bga,'Ingested')
c     Absorption and scattering coefficients should be 0 outside data range
      xtr_typ_LHS=xtr_prt_nil+xtr_fll_ngh
      xtr_typ_RHS=xtr_prt_nil+xtr_fll_ngh
      call rbn_vec(bnd_nbr_bga,wvl_grd_bga,abs_cff_mss_bga_obs,bnd_nbr,wvl_grd,abs_cff_mss_bga,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_bga,wvl_grd_bga,sca_cff_mss_bga_obs,bnd_nbr,wvl_grd,sca_cff_mss_bga,xtr_typ_LHS,xtr_typ_RHS)
      call rbn_vec(bnd_nbr_bga,wvl_grd_bga,asm_prm_bga_obs,bnd_nbr,wvl_grd,asm_prm_bga,xtr_typ_LHS,xtr_typ_RHS)

c     Ingest fl_nst
      rcd=rcd+nf90_wrp_open(fl_nst,nf90_nowrite,nc_id)
c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'wvl',bnd_dim_id)
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,bnd_dim_id,len=bnd_nbr_nst)
      if (bnd_nbr_nst > bnd_nbr_nst_max) stop 'bnd_nbr_nst > bnd_nbr_nst_max'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'nst_SRF',nst_SRF_id)
      rcd=rcd+nf90_inq_varid(nc_id,'wvl_grd',wvl_grd_nst_id)
c     Get data
      rcd=rcd+nf90_get_var(nc_id,nst_SRF_id,nst_SRF_obs)
      rcd=rcd+nf90_get_var(nc_id,wvl_grd_nst_id,wvl_grd_nst)
c     Close file
      rcd=rcd+nf90_wrp_close(nc_id,fl_nst,'Ingested')
      if (flt_nst) then
c     Instrument filter response should be 0 outside data range
         xtr_typ_LHS=xtr_prt_nil+xtr_fll_nil
         xtr_typ_RHS=xtr_prt_nil+xtr_fll_nil
         call rbn_vec(bnd_nbr_nst,wvl_grd_nst,nst_SRF_obs,bnd_nbr,wvl_grd,nst_SRF,xtr_typ_LHS,xtr_typ_RHS)
      else                      ! not flt_nst
c     Initialize a default instrument spectral response function
         call vec_set(nst_SRF,bnd_nbr,1.0)
      endif                     ! not flt_nst

c     All necessary input data has been read
c     Set wavelength bins at which vertical profile of aerosol optical depths are saved
c     fxm: This will not work with LW runs when wvl_obs ~ 0.5 um
      bnd_idx_aer=vec_val2idx(wvl_ctr,wvl_nbr,wvl_obs_aer) ! Locate aerosol bin to archive
      bnd_idx_bga=vec_val2idx(wvl_ctr,wvl_nbr,wvl_obs_bga) ! Locate background aerosol bin to archive
      
c     If necessary, alter mass path of aerosol to make user-specified extinction optical depth
c     consistent with particle optical properties
      if (cmd_ln_odxc_obs_aer) then
         do lev_idx=1,lev_nbr
            odxl_tmp=mpl_aer(lev_idx)*ext_cff_mss_aer(bnd_idx_aer)
            mpl_tmp=mpl_aer(lev_idx)
            if (odxl_tmp /= odxl_obs_aer(lev_idx)) then
               if (ext_cff_mss_aer(bnd_idx_aer) == 0.0) then
                  write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR ext_cff_mss_aer(bnd_idx_aer) = 0.0 in denominator'
                  call abort
               endif            ! endif
               mpl_aer(lev_idx)=odxl_obs_aer(lev_idx)/ext_cff_mss_aer(bnd_idx_aer)
               write (6,'(a,a,i4,a,e9.2,a,e9.2,a)') prg_nm(1:ftn_strlen(prg_nm)),
     $              ': WARNING changing mpl_aer(',lev_idx,') from = ',mpl_tmp,' to ',mpl_aer(lev_idx),' kg m-2'
            endif               ! endif adjusting mass path
         enddo                  ! end loop over lev
      endif                     ! end if overriding CLM profile odxc_obs_aer

c     Get TOA solar spectrum
      slr_spc_xtr_typ=xtr_fll_ngh ! Use xtr_fll_ngh on solar spectra
      call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
      
c     Initialize input data that may overridden by command line switches
      if (cmd_ln_odxc_obs_aer) then
         aer_lvl_nbr=0
         do lev_idx=1,lev_nbr
            if (odxl_obs_aer(lev_idx) > 0.) aer_lvl_nbr=aer_lvl_nbr+1
         enddo                  ! end loop over lev
         if (aer_lvl_nbr == 0) stop 'aer_lvl_nbr == 0'
         do lev_idx=1,lev_nbr
            if (odxl_obs_aer(lev_idx) > 0.) odxl_obs_aer(lev_idx)=odxc_obs_aer_cmd_ln/aer_lvl_nbr
         enddo                  ! end loop over lev
c     Recompute column path
         if (odxc_obs_aer_cmd_ln /= odxc_obs_aer) then
            write (6,'(a,a,f9.6,a,f9.6,a,f9.6,a)') prg_nm(1:ftn_strlen(prg_nm)),
     $           ': INFO User-specified optical depth of ',odxc_obs_aer_cmd_ln,
     $           ' differs from optical depth in input profile of ',odxc_obs_aer,
     $           ' (both taken to be at ',wvl_obs_aer*1.0e6,' um)'
            write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),
     $           ': INFO Expect WARNING messages as aerosol mass paths are changed'
         endif                  ! endif
         odxc_obs_aer=0.0
         do lev_idx=1,lev_nbr
               odxc_obs_aer=odxc_obs_aer+odxl_obs_aer(lev_idx)
         enddo                  ! end loop over lev
      endif                     ! end if overriding CLM profile odxc_obs_aer

      if (cmd_ln_mpc_CWP) then
         cld_lvl_nbr=0
         do lev_idx=1,lev_nbr
            if (mpl_CWP(lev_idx) > 0.) cld_lvl_nbr=cld_lvl_nbr+1
         enddo                  ! end loop over lev
         if (cld_lvl_nbr == 0) stop 'cld_lvl_nbr == 0'
         do lev_idx=1,lev_nbr
            if (mpl_CWP(lev_idx) > 0.) mpl_CWP(lev_idx)=mpc_CWP_cmd_ln/cld_lvl_nbr
         enddo                  ! end loop over lev
      endif                     ! end if overriding CLM profile mpl_CWP
      if (force_ice_phz) then
         call vec_set(frc_ice,lev_nbr,1.0) ! [frc]
      endif                     ! end if forcing ice phase mpl_CWP
      if (force_lqd_phz) then
         call vec_set(frc_ice,lev_nbr,0.0) ! [frc]
      endif                     ! end if forcing liquid phase mpl_CWP
      if (cmd_ln_mpc_CWP.or.force_ice_phz.or.force_lqd_phz) then
         mpc_CWP=0.0
         do lev_idx=1,lev_nbr
               mpl_IWP(lev_idx)=frc_ice(lev_idx)*mpl_CWP(lev_idx)
               mpl_LWP(lev_idx)=max(0.,mpl_CWP(lev_idx)-mpl_IWP(lev_idx))
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
      if (cmd_ln_slr_zen_ngl_cos) then
         slr_zen_ngl_cos=slr_zen_ngl_cos_cmd_ln
         xnt_fac=1.0
      endif                     ! end if overriding CLM profile zen ang
      if (cmd_ln_alb) then
         alb_sfc_NIR_drc=alb_cmd_ln
         alb_sfc_NIR_dff=alb_cmd_ln
         alb_sfc_vsb_drc=alb_cmd_ln
         alb_sfc_vsb_dff=alb_cmd_ln
         alb_sfc=alb_cmd_ln
      else
         alb_sfc=0.5*(alb_sfc_vsb_drc+alb_sfc_NIR_drc)
      endif                     ! end if overriding CLM profile albedo
c     Sanity check
      if (alb_sfc > 1.0.or.alb_sfc < 0.0) stop 'alb_sfc > 1.0.or.alb_sfc < 0.0 in nbm()'
      if (cmd_ln_slr_cst) then
         slr_cst=slr_cst_cmd_ln
      endif                     ! end if overriding solar constant

c     There is no easier way to keep O2 line absorption but turn off Herzberg continuum
c     than to zero Herzberg absorption cross-sections here.
      if (.not.flg_Herzberg) call vec_set(abs_xsx_O2,bnd_nbr,0.0)

c     Many DISORT arguments only need to be set once. 
c     Set them here, outside band loop.
      nstr=str_nbr
      nlyr=lev_nbr

c     Temperatures are specified at interfaces.
      do lev_idx=0,lev_nbr
         temper(lev_idx)=tpt_ntf(lev_idx+1)
      end do

c     NB: Setting usrtau to .true. has grave consequences on 
c     lengths of arrays returned by DISORT. In particular,
c     it is safer to keep usrtau = .false. because then all 
c     radiant quantities are returned at boundaries to 
c     computational layers.
      usrtau=.not.sv_cmp_tau
      if (sv_cmp_tau) then 
         tau_nbr=levp_nbr
      else
         tau_nbr=levp_nbr
      endif                     ! end if no user defined levels
      ntau=tau_nbr
      do tau_idx=1,tau_nbr
         tau_prs(tau_idx)=prs_ntf(tau_idx)
      end do                    ! end loop over tau

c     Set angles at which to report intensities
c     Array must be in order of increasing cosine polar angle, e.g., [-180.0,0.0] dgr
      usrang=.not.sv_cmp_plr_ngl
      numu=plr_nbr
      if (usrang) then
         plr(1)=pi
         do plr_idx=2,plr_nbr-1
            plr(plr_idx)=pi-(plr_idx-1)*pi/(plr_nbr-1)
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

c     Set azimuthal angles at which to report intensities
c     nphi=0 is valid only when onlyfl=.true.
      nphi=azi_nbr
      do azi_idx=1,azi_nbr
         azi(azi_idx)=(azi_idx-1)*2.0*pi/azi_nbr
         azi_dgr(azi_idx)=(azi_idx-1)*360.0/azi_nbr
         phi(azi_idx)=azi_dgr(azi_idx)
      enddo                     ! end loop over azimuthal angles

c     Set boundary conditions
c     0 includes case of beam illumination from top
      ibcnd=0

c     Set cosine of the solar zenith angle 
c     Beware of using 0.5 when nbr_streams/2 is odd---it's a quadrature point and can crash DISORT
      umu0=slr_zen_ngl_cos

c     Set azimuth angle of incident sunlight
      phi0=0.

c     Set incident isotropic illumination at top boundary
      fisot=0.

c     Specify bottom boundary temperature btemp 
c     Btemp is skin temperature in CCM parlance
c     Bottom boundary emissivity is derived from albedo above 
c     Top boundary temperature and emissivity must also be specified
c     Of course these are only used when plank is .true.
      btemp=tpt_skn
      ttemp=tpt_ntf(1)
      temis=0.0                  ! I have no idea what this should be

c     Set control flags:
c     Set deltam=.true. unless looking at radiances within 10 degrees of forward peak
c     When deltam=.true., the returned downwelling direct and diffuse fluxes are the "true" fluxes which have been recovered from scaled fluxes
c     Upwelling flux is (theoretically) not affected by delta scaling
c     Thus the direct/diffuse ratio computed with deltam=.true. approximately equals the ratio computed with deltam=.false. 
c     Only forward beam radiances change significantly when deltam=.true.
      deltam=.true.

c     Sanity check
      if (slr_zen_ngl_cos < 0..and..not.flg_Planck) then 
         write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR Sun beneath horizon for pure solar calculation'
         call abort
      endif

c     Set onlyfl=.false. when looking at user angle radiances.
c     When onlyfl=.true., azimuthally averaged radiances at 
c     computational angles--NOT user angles--will be reported 
c     in u0u whenever maxumu  >=  nstr.
      onlyfl=.false.

c     accur is maximum relative error in last three terms 
c     in azimuthal series, and it determines convergence.
c     Value of accur does not seem to contribute to any problems
c     with single precision.
c     Set 0.0 < accur < 0.01
      accur=0.0

c     Set printing flags for DISORT:
c     prnt(1)=.true.: Print input variables (except pmom)
      prnt(1)=.false.

c     prnt(2)=.true.: Print fluxes
      prnt(2)=.false.

c     prnt(3)=.true.: Print azimuthally-averaged intensities at user levels and computational polar angles
      prnt(3)=.false.

c     prnt(4)=.true.: Print azimuthally-averaged intensities at user levels and user polar angles
      prnt(4)=.false.

c     prnt(5)=.true.: Print intensities at user levels and angles
      prnt(5)=.false.

c     prnt(6)=.true.: Print planar transmissivity and planar albedo 
c     as a function of solar zenith angle (for ibcnd=1)
      prnt(6)=.false.

c     prnt(7)=.true.: Print phase function moments 'pmom' for each
c     layer (only if prnt(1)=.true., and only for layers with scattering)
      prnt(7)=.false.

c     Band-independent DISORT() initialization is now complete
c     Remaining DISORT arguments need to be set inside main band loop

c     Let user know where the wavelength chips have fallen
      if (dbg_lvl > dbg_off) then
         write (str_sng,'(a45,i2,a8)') 'Discrete ordinate computation performed with ',str_nbr,' streams'
         if (sv_cmp_plr_ngl) then
            write (plr_sng,'(a24,i2,a28)') 'Intensities reported at ',plr_nbr,' computational polar angles.'
         else
            write (plr_sng,'(a24,i2,a33)') 'Intensities reported at ',plr_nbr,' evenly spaced user polar angles.'
         endif
         write (azi_sng,'(a24,i2,a48,i4,a3,f9.5,a3)') 
     $        'Intensities reported at ',azi_nbr,
     $        ' evenly spaced user azimuthal angles for band = ',
     $        bnd_idx_dbg,' = ',wvl_ctr(bnd_idx_dbg)*1.0e6,' um'
         write (aer_sng,'(a46,i4,a3,f7.5,a3)') 
     $        'Layer aerosol optical depths saved for band = ',
     $        bnd_idx_aer,' = ',wvl_ctr(bnd_idx_aer)*1.0e6,' um'
         write (bga_sng,'(a57,i4,a3,f7.5,a3)') 
     $        'Layer background aerosol optical depths saved for band = ',
     $        bnd_idx_bga,' = ',wvl_ctr(bnd_idx_bga)*1.0e6,' um'
         if (sv_cmp_tau) then
            write (opt_dep_sng,'(a33,i3,a26)') 'Radiative quantities reported at ',tau_nbr,' interface pressure levels'
         else
            write (opt_dep_sng,'(a33,i3,a28)') 'Radiative quantities reported at ',tau_nbr,' user defined optical depths'
         endif
      endif                     ! endif dbg
      
      if (dbg_lvl > dbg_off) then
         write (6,'(28(a,/))')       
     $        str_sng(1:ftn_strlen(str_sng)),
     $        prf_sng(1:ftn_strlen(prf_sng)),
     $        plr_sng(1:ftn_strlen(plr_sng)),
     $        opt_dep_sng(1:ftn_strlen(opt_dep_sng)),
     $        azi_sng(1:ftn_strlen(azi_sng)),
     $        aer_sng(1:ftn_strlen(aer_sng)),
     $        bga_sng(1:ftn_strlen(bga_sng)),
     $        stt_CH4(1:ftn_strlen(stt_CH4)),
     $        stt_CO2(1:ftn_strlen(stt_CO2)),
     $        stt_H2O(1:ftn_strlen(stt_H2O)),
     $        stt_H2OH2O(1:ftn_strlen(stt_H2OH2O)),
     $        stt_Herzberg(1:ftn_strlen(stt_Herzberg)),
     $        stt_N2O(1:ftn_strlen(stt_N2O)),
     $        stt_NO2(1:ftn_strlen(stt_NO2)),
     $        stt_O2(1:ftn_strlen(stt_O2)),
     $        stt_O2N2(1:ftn_strlen(stt_O2N2)),
     $        stt_O2O2(1:ftn_strlen(stt_O2O2)),
     $        stt_O3(1:ftn_strlen(stt_O3)),
     $        stt_OH(1:ftn_strlen(stt_OH)),
     $        stt_Planck(1:ftn_strlen(stt_Planck)),
     $        stt_Rayleigh(1:ftn_strlen(stt_Rayleigh)),
     $        stt_aer(1:ftn_strlen(stt_aer)),
     $        stt_bga(1:ftn_strlen(stt_bga)),
     $        stt_flt_nst(1:ftn_strlen(stt_flt_nst)),
     $        stt_ice(1:ftn_strlen(stt_ice)),
     $        stt_lqd(1:ftn_strlen(stt_lqd)),
     $        stt_slr(1:ftn_strlen(stt_slr)),
     $        stt_top_lvl(1:ftn_strlen(stt_top_lvl))
      endif                     ! endif dbg
      
      if (dbg_lvl > dbg_off) then
         write (6,'(16(a,i4,/))')       
     $        '# input O3 bands bnd_nbr_O3 = ',bnd_nbr_O3,
     $        '# input H2OH2O bands bnd_nbr_H2OH2O = ',bnd_nbr_H2OH2O,
     $        '# input H2O narrow bands bnd_nbr_H2O = ',bnd_nbr_H2O,
     $        '# input CO2 narrow bands bnd_nbr_CO2 = ',bnd_nbr_CO2,
     $        '# input OH narrow bands bnd_nbr_OH = ',bnd_nbr_OH,
     $        '# input CH4 narrow bands bnd_nbr_CH4 = ',bnd_nbr_CH4,
     $        '# input N2O narrow bands bnd_nbr_N2O = ',bnd_nbr_N2O,
     $        '# input O2 narrow bands bnd_nbr_O2 = ',bnd_nbr_O2,
     $        '# input O2-O2 bands bnd_nbr_O2O2 = ',bnd_nbr_O2O2,
     $        '# input NO2 bands bnd_nbr_NO2 = ',bnd_nbr_NO2,
     $        '# input ice bands bnd_nbr_ice = ',bnd_nbr_ice,
     $        '# input liq bands bnd_nbr_lqd = ',bnd_nbr_lqd,
     $        '# input aer bands bnd_nbr_aer = ',bnd_nbr_aer,
     $        '# input bga bands bnd_nbr_bga = ',bnd_nbr_bga,
     $        '# input nst bands bnd_nbr_nst = ',bnd_nbr_nst,
     $        '# total bands bnd_nbr = ',bnd_nbr
      endif                     ! end if dbg

      slr_cst_xnt_fac=slr_cst*xnt_fac
      do bnd_idx=1,bnd_nbr
         odac_spc_aer(bnd_idx)=0.0
         odac_spc_bga(bnd_idx)=0.0
         odac_spc_ice(bnd_idx)=0.0
         odac_spc_lqd(bnd_idx)=0.0
         odxc_spc_CH4(bnd_idx)=0.0
         odxc_spc_CO2(bnd_idx)=0.0
         odxc_spc_H2O(bnd_idx)=0.0
         odxc_spc_H2OH2O(bnd_idx)=0.0
         odxc_spc_N2O(bnd_idx)=0.0
         odxc_spc_NO2(bnd_idx)=0.0
         odxc_spc_O2(bnd_idx)=0.0
         odxc_spc_O2N2(bnd_idx)=0.0
         odxc_spc_O2O2(bnd_idx)=0.0
         odxc_spc_O3(bnd_idx)=0.0
         odxc_spc_OH(bnd_idx)=0.0
         odxc_spc_Ray(bnd_idx)=0.0
         odxc_spc_aer(bnd_idx)=0.0
         odxc_spc_bga(bnd_idx)=0.0
         odxc_spc_ice(bnd_idx)=0.0
         odxc_spc_lqd(bnd_idx)=0.0
         odxc_spc_ttl(bnd_idx)=0.0
      enddo                     ! end loop over bnd

c     The code may run in wavelength regions outside where Malkmus band parameters (from htrn2nb) are available
c     In these regions the S_d and S_p parameters are extrapolated to be 0.0
c     An S_p of 0.0, however, causes an overflow so we reset it to unity here
      do bnd_idx=1,bnd_nbr
         if (S_d_abs_cff_mss_H2O(bnd_idx) == 0.0.and.S_p_abs_cff_mss_H2O(bnd_idx) == 0.0) S_p_abs_cff_mss_H2O(bnd_idx)=1.0
         if (S_d_abs_cff_mss_OH(bnd_idx) == 0.0.and.S_p_abs_cff_mss_OH(bnd_idx) == 0.0) S_p_abs_cff_mss_OH(bnd_idx)=1.0
         if (S_d_abs_cff_mss_CH4(bnd_idx) == 0.0.and.S_p_abs_cff_mss_CH4(bnd_idx) == 0.0) S_p_abs_cff_mss_CH4(bnd_idx)=1.0
         if (S_d_abs_cff_mss_N2O(bnd_idx) == 0.0.and.S_p_abs_cff_mss_N2O(bnd_idx) == 0.0) S_p_abs_cff_mss_N2O(bnd_idx)=1.0
         if (S_d_abs_cff_mss_CO2(bnd_idx) == 0.0.and.S_p_abs_cff_mss_CO2(bnd_idx) == 0.0) S_p_abs_cff_mss_CO2(bnd_idx)=1.0
         if (S_d_abs_cff_mss_O2(bnd_idx) == 0.0.and.S_p_abs_cff_mss_O2(bnd_idx) == 0.0) S_p_abs_cff_mss_O2(bnd_idx)=1.0
         if (single_bnd_computation.and.bnd_idx == bnd_idx_dbg) then
            write (6,'(a10,7(a3,7x))') 'Parameter','H2O','CO2','O2','OH','OH','OH'
            write (6,'(a10,7(es10.3,1x))') 'S_d_abs_cff_mss',
     $           S_d_abs_cff_mss_H2O(bnd_idx),S_d_abs_cff_mss_CO2(bnd_idx),S_d_abs_cff_mss_O2(bnd_idx),
     $           S_d_abs_cff_mss_OH(bnd_idx),S_d_abs_cff_mss_CH4(bnd_idx),S_d_abs_cff_mss_N2O(bnd_idx)
            write (6,'(1p,a10,7(es10.3,1x))') 'S_p_abs_cff_mss',
     $           S_p_abs_cff_mss_H2O(bnd_idx),S_p_abs_cff_mss_CO2(bnd_idx),S_p_abs_cff_mss_O2(bnd_idx),
     $           S_p_abs_cff_mss_OH(bnd_idx),S_p_abs_cff_mss_CH4(bnd_idx),S_p_abs_cff_mss_N2O(bnd_idx)
         endif                  ! endif dbg
         bnd_var(bnd_idx)=S_p_abs_cff_mss_OH(bnd_idx) ! Generic spectral variable
      enddo                     ! end loop over bnd

c     Initialize band-independent arrays that depend on level
      do lev_idx=1,lev_nbr
         tpt_dlt_Mlk(lev_idx)=tpt(lev_idx)-tpt_Malkmus_fit
         tpt_dlt_Mlk_sqr(lev_idx)=tpt_dlt_Mlk(lev_idx)*tpt_dlt_Mlk(lev_idx)
      enddo                     ! end loop over lev

c     Zero broad-band arrays which will accumulate spectral fluxes
      call vec_set(j_NO2,lev_nbr,0.0) ! [s-1]
      call vec_set(ntn_bb_mean,lev_nbr,0.0) ! [W m-2 sr-1]
      call vec_set(ntn_bb_aa,plr_nbr*levp_nbr,0.0) ! [W m-2 sr-1]
      call vec_set(flx_bb_dwn_drc,levp_nbr,0.0) ! [W m-2]
      call vec_set(flx_bb_dwn_dff,levp_nbr,0.0) ! [W m-2]
      call vec_set(flx_bb_up,levp_nbr,0.0) ! [W m-2]

c     End section 1: Initialization
c     Begin section 2: Main computation loop over all bands
c     begin main loop over bands
c     assignment of counting indices into input data arrays
c     a level loop for zeroing band quantities that depend on level
c     computation and binning of process-specific optical properties
c     summation of individual optical properties into bulk layer props
c     assignment of bulk properties into actual DISORT() input arrays
c     call DISORT()
c     assignment of DISORT() output to final output storage arrays
c     end main loop over bands

      do bnd_idx=1,bnd_nbr

         if (dbg_lvl > dbg_off) write (6,'(a1)',advance="no") '.'
         if (dbg_lvl == dbg_old) write (6,'(i4)',advance="no") bnd_idx

c     Zero out all optical depth arrays before computing anything for this band 
c     This is important as defensive programming and when arguments are given to ignore specified radiative processes, e.g., the pure Rayleigh scattering case
c     Also, O2, O3, (and, eventually, H2O) cause both line and continuum absorption
         do lev_idx=1,lev_nbr
            odal_CH4(lev_idx)=0.0
            odal_CO2(lev_idx)=0.0
            odal_H2O(lev_idx)=0.0
            odal_H2OH2O(lev_idx)=0.0
            odal_N2O(lev_idx)=0.0
            odal_NO2(lev_idx)=0.0
            odal_O2(lev_idx)=0.0
            odal_O2N2(lev_idx)=0.0
            odal_O2O2(lev_idx)=0.0
            odal_O3(lev_idx)=0.0
            odal_OH(lev_idx)=0.0
            odal_aer(lev_idx)=0.0
            odal_bga(lev_idx)=0.0
            odal_ice(lev_idx)=0.0
            odal_lqd(lev_idx)=0.0
            odsl_Ray(lev_idx)=0.0
            odsl_aer(lev_idx)=0.0
            odsl_bga(lev_idx)=0.0
            odsl_ice(lev_idx)=0.0
            odsl_lqd(lev_idx)=0.0
         enddo                  ! end loop over lev

c     Begin line absorption
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_H2O(bnd_idx),A_psi_H2O(bnd_idx),B_phi_H2O(bnd_idx),B_psi_H2O(bnd_idx), ! Input
     $     S_d_abs_cff_mss_H2O(bnd_idx),S_p_abs_cff_mss_H2O(bnd_idx), ! Input
     $        grv,q_H2O,mpl_H2O,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_H2O)         ! Output
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_OH(bnd_idx),A_psi_OH(bnd_idx),B_phi_OH(bnd_idx),B_psi_OH(bnd_idx), ! Input
     $        S_d_abs_cff_mss_OH(bnd_idx),S_p_abs_cff_mss_OH(bnd_idx), ! Input
     $        grv,q_OH,mpl_OH,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_OH)          ! Output
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_O2(bnd_idx),A_psi_O2(bnd_idx),B_phi_O2(bnd_idx),B_psi_O2(bnd_idx), ! Input
     $        S_d_abs_cff_mss_O2(bnd_idx),S_p_abs_cff_mss_O2(bnd_idx), ! Input
     $        grv,q_O2,mpl_O2,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_O2)          ! Output
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_CO2(bnd_idx),A_psi_CO2(bnd_idx),B_phi_CO2(bnd_idx),B_psi_CO2(bnd_idx), ! Input
     $        S_d_abs_cff_mss_CO2(bnd_idx),S_p_abs_cff_mss_CO2(bnd_idx), ! Input
     $        grv,q_CO2,mpl_CO2,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_CO2)         ! Output
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_CH4(bnd_idx),A_psi_CH4(bnd_idx),B_phi_CH4(bnd_idx),B_psi_CH4(bnd_idx), ! Input
     $        S_d_abs_cff_mss_CH4(bnd_idx),S_p_abs_cff_mss_CH4(bnd_idx), ! Input
     $        grv,q_CH4,mpl_CH4,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_CH4)         ! Output
         call mlk_abs(
     $        lev_nbr,levp_nbr,prs_HITRAN,slr_zen_ngl_cos, ! Input
     $        A_phi_N2O(bnd_idx),A_psi_N2O(bnd_idx),B_phi_N2O(bnd_idx),B_psi_N2O(bnd_idx), ! Input
     $        S_d_abs_cff_mss_N2O(bnd_idx),S_p_abs_cff_mss_N2O(bnd_idx), ! Input
     $        grv,q_N2O,mpl_N2O,prs,prs_ntf,tpt_dlt_Mlk,tpt_dlt_Mlk_sqr, ! Input
     $        odal_N2O)         ! Output
c     End line absorption
         
         if (.false.) then
            if (dbg_lvl>dbg_off) then
               write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': DEBUG O3 cross-section'
               write (6,'(a,i4,a,es8.1)') 'abs_xsx_O3(',bnd_idx,') = ',abs_xsx_O3(bnd_idx)
               write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': DEBUG npl_O3'
               do lev_idx=1,lev_nbr
                  write (6,'(a,i4,a,es8.1)') 'npl_O3(',lev_idx,') = ',npl_O3(lev_idx)
               enddo            ! end loop over bnd
            endif               ! endif dbg
         endif                  ! endif
         
c     Begin continuum processes
c     O3 continuum absorption
         do lev_idx=1,lev_nbr
            odal_O3(lev_idx)=abs_xsx_O3(bnd_idx)*npl_O3(lev_idx)
            odal_O2(lev_idx)=odal_O2(lev_idx)+abs_xsx_O2(bnd_idx)*npl_O2(lev_idx)
         enddo                  ! end loop over lev

c     H2OH2O continuum absorption
c     call odal_H2OH2O_Chy97(wvn_ctr(bnd_idx),lev_nbr,t,mpl_H2O,RH_lqd,odal_H2OH2O,dbg_lvl)
c     abs_xsx_H2OH2O(bnd_idx)=abs_xsx_H2OH2O_CFT99(wvn_ctr(bnd_idx))
c     abs_xsx_H2OH2O(bnd_idx)=abs_xsx_H2OH2O_Chy97(wvn_ctr(bnd_idx))
         do lev_idx=1,lev_nbr
            odal_H2OH2O(lev_idx)=abs_xsx_H2OH2O(bnd_idx)*npl_H2OH2O(lev_idx)
         enddo                  ! end loop over lev
         
c     O2-O2 continuum absorption
         do lev_idx=1,lev_nbr
            odal_O2O2(lev_idx)=abs_xsx_O2O2(bnd_idx)*npl_O2O2(lev_idx)
         enddo                  ! end loop over lev
         
c     O2-N2 continuum absorption
c     Compute 02-N2 absorption optical depth based on O2-O2 optical depth
c     1.27 um band O2-O2 has measured FWHM = 182.5 cm-1 and intensity 1.60 times O2-O2 1.06 um band 
c     Restrict scaling of O2-N2 absorption in 1.27 um band to between 1.2 and 1.35 um
c     SPS98 Figure 1a shows O2-N2 absorption cross sections in units of m5 mlc-2
c     These funny units represent the absorption cross section of O2-N2 per unit concentration of O2
c     Since the ratio of O2 to N2 in the atmosphere is fixed, we simply scale the absorption optical
c     depth of O2-O2 to O2-N2, rather than scaling the cross sections themselves.
         if ((wvl_ctr(bnd_idx) >= 1.2e-6).and.(wvl_ctr(bnd_idx) <= 1.35e-6)) then
            do lev_idx=1,lev_nbr
c     fxm: Should make efficiency factor (0.2) a parameter or cmd ln arg
               odal_O2N2(lev_idx)=odal_O2O2(lev_idx)*0.2*N2_per_O2 ! O2-N2 absorption efficiency is 20% of O2-O2 SPS98 p. 12
            enddo               ! end loop over lev
         endif                  ! endif 1.27 um band
         
c     NO2 continuum absorption
         do lev_idx=1,lev_nbr
            odal_NO2(lev_idx)=abs_xsx_NO2(bnd_idx)*npl_NO2(lev_idx)
         enddo                  ! end loop over lev
         
c     Compute liquid Mie scattering/absorption optical depths
         do lev_idx=1,lev_nbr
            odal_lqd(lev_idx)=abs_cff_mss_lqd(bnd_idx)*mpl_LWP(lev_idx)
            odsl_lqd(lev_idx)=sca_cff_mss_lqd(bnd_idx)*mpl_LWP(lev_idx)
         enddo                  ! end loop over lev

c     Compute ice Mie scattering/absorption optical depths
         do lev_idx=1,lev_nbr
            odal_ice(lev_idx)=abs_cff_mss_ice(bnd_idx)*mpl_IWP(lev_idx)
            odsl_ice(lev_idx)=sca_cff_mss_ice(bnd_idx)*mpl_IWP(lev_idx)
         enddo                  ! end loop over lev

c     Compute aerosol Mie scattering/absorption optical depths
         do lev_idx=1,lev_nbr
            odal_aer(lev_idx)=abs_cff_mss_aer(bnd_idx)*mpl_aer(lev_idx)
            odsl_aer(lev_idx)=sca_cff_mss_aer(bnd_idx)*mpl_aer(lev_idx)
         enddo                  ! end loop over lev

c     Compute background aerosol Mie scattering/absorption optical depths
         do lev_idx=1,lev_nbr
            odal_bga(lev_idx)=abs_cff_mss_bga(bnd_idx)*mpl_bga(lev_idx)
            odsl_bga(lev_idx)=sca_cff_mss_bga(bnd_idx)*mpl_bga(lev_idx)
         enddo                  ! end loop over lev

c     Figure out Rayleigh scattering optical depth
         do lev_idx=1,lev_nbr

c     See Len93 p. 154 for details. See also BrS84 p. 107, GoY89 p. 297
            idx_rfr_air_STP(bnd_idx)=
     $           1.0+
     $           1.0e-6*(77.46+0.459/(1.0e12*wvl_ctr(bnd_idx)**2))*
     $           prs_STP*0.01/tpt_STP
            float_foo=32.0*(idx_rfr_air_STP(bnd_idx)-1.0)**2/3.0

c     Place factor of N_STP so as to avoid overflow in single precision
            sca_cff_mss_Ray(lev_idx)=
     $           pi**3*float_foo*(Avagadro/N_STP)/
     $           (mmw_mst_air(lev_idx)*N_STP*wvl_ctr(bnd_idx)**4)
            odsl_Ray(lev_idx)=
     $           mpl_mst_air(lev_idx)*sca_cff_mss_Ray(lev_idx)
         enddo                  ! end loop over lev
c     End continuum processes

c     Individual optical properties for each radiative process have been computed 
c     See if user specified any processes to intentionally neglect
c     Zero out appropriate optical constants calculate bulk layer properties normally 
c     NB: Mie asymmetry parameter never needs to be zeroed (e.g., in clear sky case)
c     because it is always weighted by Mie optical depth before being used in a computation.
         if (.not.flg_Rayleigh) then
            call vec_set(odsl_Ray,lev_nbr,0.0) ! [frc]
         endif                  ! end if no Rayleigh processes
         if (.not.flg_ice) then
            call vec_set(odal_ice,lev_nbr,0.0) ! [frc]
            call vec_set(odsl_ice,lev_nbr,0.0) ! [frc]
         endif                  ! end if no ice processes
         if (.not.flg_lqd) then
            call vec_set(odal_lqd,lev_nbr,0.0) ! [frc]
            call vec_set(odsl_lqd,lev_nbr,0.0) ! [frc]
         endif                  ! end if no liq processes
         if (.not.flg_aer) then
            call vec_set(odal_aer,lev_nbr,0.0) ! [frc]
            call vec_set(odsl_aer,lev_nbr,0.0) ! [frc]
         endif                  ! end if no aer processes
         if (.not.flg_bga) then
            call vec_set(odal_bga,lev_nbr,0.0) ! [frc]
            call vec_set(odsl_bga,lev_nbr,0.0) ! [frc]
         endif                  ! end if no bga processes
         if (.not.flg_H2O) then
            call vec_set(odal_H2O,lev_nbr,0.0) ! [frc]
         endif                  ! end if no H2O processes
         if (.not.flg_H2OH2O) then
            call vec_set(odal_H2OH2O,lev_nbr,0.0) ! [frc]
         endif                  ! end if no H2OH2O processes
         if (.not.flg_CO2) then
            call vec_set(odal_CO2,lev_nbr,0.0) ! [frc]
         endif                  ! end if no CO2 processes
         if (.not.flg_OH) then
            call vec_set(odal_OH,lev_nbr,0.0) ! [frc]
         endif                  ! end if no OH processes
         if (.not.flg_CH4) then
            call vec_set(odal_CH4,lev_nbr,0.0) ! [frc]
         endif                  ! end if no OH processes
         if (.not.flg_N2O) then
            call vec_set(odal_N2O,lev_nbr,0.0) ! [frc]
         endif                  ! end if no OH processes
         if (.not.flg_O2) then
            call vec_set(odal_O2,lev_nbr,0.0) ! [frc]
         endif                  ! end if no O2 processes
         if (.not.flg_O3) then
            call vec_set(odal_O3,lev_nbr,0.0) ! [frc]
         endif                  ! end if no O3 processes
         if (.not.flg_O2O2) then
            call vec_set(odal_O2O2,lev_nbr,0.0) ! [frc]
         endif                  ! end if no O2-O2 processes
         if (.not.flg_O2N2) then
            call vec_set(odal_O2N2,lev_nbr,0.0) ! [frc]
         endif                  ! end if no O2N2 processes
         if (.not.flg_NO2) then
            call vec_set(odal_NO2,lev_nbr,0.0) ! [frc]
         endif                  ! end if no NO2 processes

c     Now that processes have been turned on/off, save diagnostic values
c     Note that, lamely, bnd_idx_aer in nbm is bnd_obs_aer in swnb
         if (bnd_idx == bnd_idx_aer) then
            do lev_idx=1,lev_nbr
               odal_obs_aer(lev_idx)=odal_aer(lev_idx)
               odsl_obs_aer(lev_idx)=odsl_aer(lev_idx)
               odxl_obs_aer(lev_idx)=odal_aer(lev_idx)+odsl_aer(lev_idx)
            enddo               ! end loop over lev
         endif                  ! end if bnd_idx_aer
         if (bnd_idx == bnd_idx_bga) then
            do lev_idx=1,lev_nbr
               odal_obs_bga(lev_idx)=odal_bga(lev_idx)
               odsl_obs_bga(lev_idx)=odsl_bga(lev_idx)
               odxl_obs_bga(lev_idx)=odal_bga(lev_idx)+odsl_bga(lev_idx)
            enddo               ! end loop over lev
         endif                  ! end if bnd_idx_bga

c     Weight optical parameters as per CCY83
         do lev_idx=1,lev_nbr
            odsl_spc_ttl(bnd_idx,lev_idx)=
     $           odsl_Ray(lev_idx)+
     $           odsl_lqd(lev_idx)+odsl_ice(lev_idx)+
     $           odsl_aer(lev_idx)+odsl_bga(lev_idx)
            odal_spc_ttl(bnd_idx,lev_idx)=
     $           odal_H2O(lev_idx)+odal_CO2(lev_idx)+
     $           odal_O2(lev_idx)+odal_O3(lev_idx)+
     $           odal_NO2(lev_idx)+odal_OH(lev_idx)+
     $           odal_CH4(lev_idx)+odal_N2O(lev_idx)+
     $           odal_O2O2(lev_idx)+odal_O2N2(lev_idx)+
     $           odal_H2OH2O(lev_idx)+
     $           odal_lqd(lev_idx)+odal_ice(lev_idx)+
     $           odal_aer(lev_idx)+odal_bga(lev_idx)
            odxl_spc_ttl(bnd_idx,lev_idx)=
     $           odsl_spc_ttl(bnd_idx,lev_idx)+odal_spc_ttl(bnd_idx,lev_idx)
         enddo                  ! end loop over lev

c     In rare, diagnostic cases where all scattering is turned off,
c     the following lines have a zero over zero condition and
c     will crash under many compilers unless protected by conditions.
         do lev_idx=1,lev_nbr
            if (odsl_spc_ttl(bnd_idx,lev_idx) <= 0.0) then
               ss_alb_fct(bnd_idx,lev_idx)=0.0
               sca_frc_ttl(lev_idx)=0.0
               asm_prm_ttl(bnd_idx,lev_idx)=0.0
            else                ! endif no scattering whatsoever

c     Single scattering albedo is ill-conditioned when 
c     scattering optical depth is zero and there is no absorption.
c     This may only be a problem in single precision---I'm not sure
               ss_alb_fct(bnd_idx,lev_idx)=
     $              odsl_spc_ttl(bnd_idx,lev_idx)/
     $              odxl_spc_ttl(bnd_idx,lev_idx)

c     Mie scattering asymmetry parameter and fraction of
c     scattering due to Mie scattering must be saved for all 
c     bands and levels until DISORT() is called. 
c     They are needed to correctly weight the moments of the total phase function 
c     between Rayleigh and Henyey-Greenstein components.
               sca_frc_ttl(lev_idx)=
     $              (odsl_lqd(lev_idx)+odsl_ice(lev_idx)+
     $              odsl_aer(lev_idx)+odsl_bga(lev_idx))/
     $              odsl_spc_ttl(bnd_idx,lev_idx)

c     CZP III p. #115 has discussion of effective asymmetry parameter
               asm_prm_ttl(bnd_idx,lev_idx)=
     $              (asm_prm_lqd(bnd_idx)*odsl_lqd(lev_idx)+
     $              asm_prm_ice(bnd_idx)*odsl_ice(lev_idx)+
     $              asm_prm_aer(bnd_idx)*odsl_aer(lev_idx)+
     $              asm_prm_bga(bnd_idx)*odsl_bga(lev_idx))/
     $              odsl_spc_ttl(bnd_idx,lev_idx)
            endif               ! endif there is scattering
         enddo                  ! end loop over lev

c     Sanity check for unphysical single scattering albedos
c     omegas as large as 1.0000656 can occur (at least under LINUX) for very
c     small Rayleigh scattering optical depths (p <~ 1 Pa) in single precision
         do lev_idx=1,lev_nbr
            if (ss_alb_fct(bnd_idx,lev_idx) > 1.0) then
               write (6,'(a,a,i4,a,i3,a,f10.7)') 
     $              prg_nm(1:ftn_strlen(prg_nm)),': WARNING ss_alb_fct(',bnd_idx,',',lev_idx,') = ',ss_alb_fct(bnd_idx,lev_idx)
               ss_alb_fct(bnd_idx,lev_idx)=1.0
            endif
         enddo                  ! end loop over lev
            
c     Compute diagnostics
         do lev_idx=1,lev_nbr
            odac_spc_aer(bnd_idx)=odac_spc_aer(bnd_idx)+odal_aer(lev_idx)
            odac_spc_bga(bnd_idx)=odac_spc_bga(bnd_idx)+odal_bga(lev_idx)
            odac_spc_ice(bnd_idx)=odac_spc_ice(bnd_idx)+odal_ice(lev_idx)
            odac_spc_lqd(bnd_idx)=odac_spc_lqd(bnd_idx)+odal_lqd(lev_idx)
            odxc_spc_CO2(bnd_idx)=odxc_spc_CO2(bnd_idx)+odal_CO2(lev_idx)
            odxc_spc_H2O(bnd_idx)=odxc_spc_H2O(bnd_idx)+odal_H2O(lev_idx)
            odxc_spc_H2OH2O(bnd_idx)=odxc_spc_H2OH2O(bnd_idx)+odal_H2OH2O(lev_idx)
            odxc_spc_NO2(bnd_idx)=odxc_spc_NO2(bnd_idx)+odal_NO2(lev_idx)
            odxc_spc_O2(bnd_idx)=odxc_spc_O2(bnd_idx)+odal_O2(lev_idx)
            odxc_spc_O2N2(bnd_idx)=odxc_spc_O2N2(bnd_idx)+odal_O2N2(lev_idx)
            odxc_spc_O2O2(bnd_idx)=odxc_spc_O2O2(bnd_idx)+odal_O2O2(lev_idx)
            odxc_spc_O3(bnd_idx)=odxc_spc_O3(bnd_idx)+odal_O3(lev_idx)
            odxc_spc_OH(bnd_idx)=odxc_spc_OH(bnd_idx)+odal_OH(lev_idx)
            odxc_spc_CH4(bnd_idx)=odxc_spc_CH4(bnd_idx)+odal_CH4(lev_idx)
            odxc_spc_N2O(bnd_idx)=odxc_spc_N2O(bnd_idx)+odal_N2O(lev_idx)
            odxc_spc_Ray(bnd_idx)=odxc_spc_Ray(bnd_idx)+odsl_Ray(lev_idx)
            odxc_spc_aer(bnd_idx)=odxc_spc_aer(bnd_idx)+odal_aer(lev_idx)+odsl_aer(lev_idx)
            odxc_spc_bga(bnd_idx)=odxc_spc_bga(bnd_idx)+odal_bga(lev_idx)+odsl_bga(lev_idx)
            odxc_spc_ice(bnd_idx)=odxc_spc_ice(bnd_idx)+odal_ice(lev_idx)+odsl_ice(lev_idx)
            odxc_spc_lqd(bnd_idx)=odxc_spc_lqd(bnd_idx)+odal_lqd(lev_idx)+odsl_lqd(lev_idx)
         enddo                  ! end loop over lev
         odxc_spc_ttl(bnd_idx)=
     $        odxc_spc_CO2(bnd_idx)+
     $        odxc_spc_H2O(bnd_idx)+
     $        odxc_spc_H2OH2O(bnd_idx)+
     $        odxc_spc_NO2(bnd_idx)+
     $        odxc_spc_O2(bnd_idx)+
     $        odxc_spc_O2N2(bnd_idx)+
     $        odxc_spc_O2O2(bnd_idx)+
     $        odxc_spc_O3(bnd_idx)+
     $        odxc_spc_OH(bnd_idx)+
     $        odxc_spc_CH4(bnd_idx)+
     $        odxc_spc_N2O(bnd_idx)+
     $        odxc_spc_Ray(bnd_idx)+
     $        odxc_spc_aer(bnd_idx)+
     $        odxc_spc_bga(bnd_idx)+
     $        odxc_spc_ice(bnd_idx)+
     $        odxc_spc_lqd(bnd_idx)

c     All spectral information needed for DISORT and data reduction has now been saved
c     Massage data for input to DISORT 
c     DISORT uses a top down numbering scheme like CCM
c     Layer indices from HPMM cirrus cloud model must be reversed

c     Assign optical depth and single scattering albedos
         do lev_idx=1,lev_nbr
            dtauc(lev_idx)=odxl_spc_ttl(bnd_idx,lev_idx)
            ssalb(lev_idx)=ss_alb_fct(bnd_idx,lev_idx)
         end do

c     Compute phase function moments (currently a weighted mixture of Henyey-Greenstein and Rayleigh phase functions)
         do lev_idx=1,lev_nbr
            pmom(0,lev_idx)=1.0
         end do

c     All moments of Rayleigh phase function except second are zero
c     To blend this with HG phase function from Mie
c     scattering we weight HG coefficients by fraction of
c     scattering due to Mie processes.
         do mom_idx=1,nbr_mom
            do lev_idx=1,lev_nbr
               pmom(mom_idx,lev_idx)=
     $              sca_frc_ttl(lev_idx)*
     $              asm_prm_ttl(bnd_idx,lev_idx)**mom_idx
            end do
         end do

c     For the second moment, we form a weighted average of 
c     HG scattering fraction times the HG moment (g^2), 
c     plus Rayleigh scattered fraction times Rayleigh moment (0.1).
         do lev_idx=1,lev_nbr
            pmom(2,lev_idx)=
     $           (1.0-sca_frc_ttl(lev_idx))*0.1+
     $           sca_frc_ttl(lev_idx)*
     $           asm_prm_ttl(bnd_idx,lev_idx)**2
         end do
         if (.not.sv_cmp_tau) then 

c     If usrtau == true. then place levels at top and 
c     bottom of atmosphere and at computational levels in between. 
c     This should be same as setting usrtau=.false.
c     There is not much point in setting user levels to anything
c     else, i.e., evenly spaced in optical depth, or at top and
c     bottom of cloud, yet. But there may be soon. In that case
c     all radiant quantities stored in netCDF will need to be 
c     made dimensioned by tau_nbr instead of levp_nbr.
            tau(1)=0.0
            float_foo=0.0
            do tau_idx=1,tau_nbr
               float_foo=float_foo+odxl_spc_ttl(bnd_idx,tau_idx)
               tau(tau_idx+1)=float_foo
            end do              ! end loop over tau

            ntau=tau_nbr
            do tau_idx=1,tau_nbr
               utau(tau_idx)=tau(tau_idx)
            end do              ! end loop over tau
         endif                  ! end if setting user defined levels

c     Set wavenumbers for this spectral interval.
         wvnmlo=wvn_min(bnd_idx)
         wvnmhi=wvn_max(bnd_idx)

c     Set plank=.true. whenever considering thermal emission.
c     Computing thermal source function generates underflows in single precision when lambda < XXX um.
         if (wvl_ctr(bnd_idx) > wvl_Planck) then
            plank=flg_Planck
         else 
            plank=.false.
         endif

c     Set intensity of incident parallel beam at top boundary.
c     Units are arbitrary but must match fisot in solar case.
c     Must be in W/m^2 if there is any thermal emission (plank=.true.)
c     Recall that fbeam is flux normal to earth-sun path,
c     not normal to ground, i.e., zenith angle correction
c     to incoming flux is made within DISORT() so DO NOT
c     pre-multiply fbeam by cosine of the solar zenith angle.
         fbeam=slr_cst_xnt_fac*flx_slr_frc(bnd_idx)

c     Give details about lower boundary reflectance.
c     When lamber is true, albedo specifies isotropic reflectance.
c     Otherwise hl array must be specified to give BDRF of bottom boundary.
         lamber=.true.
         if (wvl_ctr(bnd_idx) < .7e-6) then
            if (slr_zen_ngl_cos > .5) then
               albedo=alb_sfc_vsb_drc
            else
               albedo=alb_sfc_vsb_dff
            endif
         else
            if (slr_zen_ngl_cos > .5) then
               albedo=alb_sfc_NIR_drc
            else
               albedo=alb_sfc_NIR_dff
            endif
         endif
      
         if (single_bnd_computation) then

c     Set header which DISORT will use printing results
c     If header has length greater than zero, DISORT() prints annoying header message on each call
            write (header,'(a,i4,a,f6.4,a,es8.1,a,es8.1)')
     $           'nbm: i = ',bnd_idx,', lambda = ',
     $           wvl_ctr(bnd_idx)*1.0e6,
     $           ', flx_slr_frc = ',flx_slr_frc(bnd_idx),', bandwith = ',
     $           wvl_dlt(bnd_idx)

            if (bnd_idx == bnd_idx_dbg) then

c     Hardcode numbers to test radiative code against benchmark values
               if (tst_case_Rayleigh) then

c     Pure Rayleigh scattering. Answer: Albedo = 0.796920
                  fbeam=pi
                  albedo=0.0
                  umu0=0.5
                  do lev_idx=1,lev_nbr
                     dtauc(lev_idx)=4.284/lev_nbr
                     odsl_Ray(lev_idx)=dtauc(lev_idx)
                     odsl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     odxl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                  enddo         ! end loop over lev
               endif            ! end if tst_case_Rayleigh
               if (tst_case_HG) then

c     Pure Henyey-Greenstein scattering. Answer: Albedo = .123420
c     Atmospheric absorptance = 0.360522
                  fbeam=pi
                  albedo=0.0
                  umu0=0.5
                  do lev_idx=1,lev_nbr
                     dtauc(lev_idx)=1.0/lev_nbr
                     ssalb(lev_idx)=0.8
                     odsl_Ray(lev_idx)=dtauc(lev_idx)
                     odsl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     odxl_spc_ttl(bnd_idx,lev_idx)=dtauc(lev_idx)
                     ss_alb_fct(bnd_idx,lev_idx)=ssalb(lev_idx)
                     sca_frc_ttl(lev_idx)=1.0
                     asm_prm_lqd(bnd_idx)=0.75
                     asm_prm_ttl(bnd_idx,lev_idx)=0.75
                     odsl_Ray(lev_idx)=0.
                  enddo         ! end loop over lev
                  do lev_idx=1,lev_nbr
                     pmom(0,lev_idx)=1.0
                  end do
                  do mom_idx=1,nbr_mom
                     do lev_idx=1,lev_nbr
                        pmom(mom_idx,lev_idx)=0.75**mom_idx
                     end do
                  end do
               endif            ! end if tst_case_HG
               if (dbg_lvl > dbg_off) then
                  write (6,'(/,a,i4,a4,f9.6,a2,f9.6,a11,f9.3,a2,f9.3,a5)')
     $                 'Detailed debugging information for bnd(',bnd_idx,') = ',
     $                 wvl_min(bnd_idx)*1.0e6,'--',wvl_max(bnd_idx)*1.0e6,' um = ',
     $                 wvnmlo,'--',wvnmhi,' cm-1'
                  write (6,'(a4,18(1x,a8))')
     $                 'lev_idx',
     $                 'prs_mdp',
     $                 'tpt_mdp',
     $                 'tau_sca',
     $                 'tau_abs',
     $                 'tau_ext',
     $                 'ss_alb',
     $                 'Mie_frc',
     $                 'lqd_asym',
     $                 'tot_asym',
     $                 'Ray_msig',
     $                 'pth_wet',
     $                 'Ray_tau',
     $                 'dtauc',
     $                 'ssalb',
     $                 'pmom(1)',
     $                 'pmom(2)',
     $                 'pmom(3)',
     $                 'pmom(4)'
                  do lev_idx=1,lev_nbr
                     write (6,'(1p,i4,18(1x,es8.1))')
     $                    lev_idx,
     $                    prs(lev_idx),
     $                    tpt(lev_idx),
     $                    odsl_spc_ttl(bnd_idx,lev_idx),
     $                    odal_spc_ttl(bnd_idx,lev_idx),
     $                    odxl_spc_ttl(bnd_idx,lev_idx),
     $                    ss_alb_fct(bnd_idx,lev_idx),
     $                    sca_frc_ttl(lev_idx),
     $                    asm_prm_lqd(bnd_idx),
     $                    asm_prm_ttl(bnd_idx,lev_idx),
     $                    sca_cff_mss_Ray(lev_idx),
     $                    mpl_mst_air(lev_idx),
     $                    odsl_Ray(lev_idx),
     $                    dtauc(lev_idx),
     $                    ssalb(lev_idx),
     $                    pmom(1,lev_idx),
     $                    pmom(2,lev_idx),
     $                    pmom(3,lev_idx),
     $                    pmom(4,lev_idx)
                  enddo         ! end loop over lev
               endif            ! endif dbg
            else                ! end if dbg band
               goto 999
            endif               ! end if not dbg band
         endif                  ! end if single_bnd_computation
         
         call DISORT( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,
     $        wvnmhi, usrtau, ntau, utau, nstr, usrang,
     $        numu, umu, nphi, phi, ibcnd, fbeam, umu0,
     $        phi0, fisot, lamber, albedo, hl, btemp, ttemp,
     $        temis, deltam, plank, onlyfl, accur, prnt,
     $        header, maxcly, maxulv, maxumu, maxcmu,
     $        maxphi, rfldir, rfldn, flup, dfdt, uavg,
     $        uu, u0u, albmed, trnmed )

c     Store spectral fluxes returned from DISORT
c     DISORT uses a funny indexing system: 
c     Returned flux arrays, which are on layer interfaces, 
c     are always indexed starting from 1 at TOA even though same levels 
c     start with index 0 on input to DISORT, i.e., temper().
c     Make sure you know what you're doing before changing these indices
         do lev_idx=1,lev_nbr
            ntn_spc_mean(bnd_idx,lev_idx)=
     $           0.5*(uavg(lev_idx)+uavg(lev_idx+1))/
     $           wvl_dlt(bnd_idx)
         enddo                  ! end loop over lev
         do lev_idx=1,levp_nbr
            flx_spc_dwn_drc(bnd_idx,lev_idx)=rfldir(lev_idx)/wvl_dlt(bnd_idx)
            flx_spc_dwn_dff(bnd_idx,lev_idx)=rfldn(lev_idx)/wvl_dlt(bnd_idx)
            flx_spc_dwn(bnd_idx,lev_idx)=
     $           flx_spc_dwn_drc(bnd_idx,lev_idx)+
     $           flx_spc_dwn_dff(bnd_idx,lev_idx)
            flx_spc_up(bnd_idx,lev_idx)=flup(lev_idx)/wvl_dlt(bnd_idx)
            flx_spc_net(bnd_idx,lev_idx)=flx_spc_dwn(bnd_idx,lev_idx)-flx_spc_up(bnd_idx,lev_idx)
         enddo                  ! end loop over lev

c     u0u returns azimuthally averaged intensities at user polar angles 
c     uu returns exact intensities at each user azimuthal angle at each user polar angle
c     User polar angles are ordered by increasing cosine of polar angle, 
c     i.e., from downwelling radiance to upwelling radiance. 

c     NB: Possible bug in current DISORT() documentation 
c     Intensities are returned in arrays in order of descending cosine of polar angle
         ntn_spc_aa_ndr_sfc(bnd_idx)=u0u(1,levp_nbr)/wvl_dlt(bnd_idx)
         ntn_spc_aa_zen_sfc(bnd_idx)=u0u(plr_nbr,levp_nbr)/wvl_dlt(bnd_idx)

         do lev_idx=1,levp_nbr
            do plr_idx=1,plr_nbr
               ntn_bb_aa(plr_idx,lev_idx)=ntn_bb_aa(plr_idx,lev_idx)+
     $              u0u(plr_idx,lev_idx)
            enddo               ! end loop over plr
            ntn_spc_aa_ndr(bnd_idx,lev_idx)=
     $           u0u(1,lev_idx)/wvl_dlt(bnd_idx)
            ntn_spc_aa_zen(bnd_idx,lev_idx)=
     $           u0u(plr_nbr,lev_idx)/wvl_dlt(bnd_idx)
         enddo                  ! end loop over lev

         do plr_idx=1,plr_nbr
            ntn_spc_aa_sfc(plr_idx,bnd_idx)=
     $           u0u(plr_idx,levp_nbr)/wvl_dlt(bnd_idx)
         enddo                  ! end loop over lev

         if (bnd_idx == bnd_idx_dbg) then
            do azi_idx=1,azi_nbr
               do plr_idx=1,plr_nbr
                  do lev_idx=1,levp_nbr
                     ntn_spc_chn(azi_idx,plr_idx,lev_idx)=
     $                    uu(plr_idx,lev_idx,azi_idx)/wvl_dlt(bnd_idx)
                  enddo         ! end loop over lev
               enddo            ! end loop over plr
            enddo               ! end loop over azimuthal angles
         endif                  ! end if dbg band

         if (sv_ntn) then
c            do azi_idx=1,azi_nbr
c               do plr_idx=1,plr_nbr
c                  do lev_idx=1,levp_nbr
c                     ntn_spc_aa(plr_idx,bnd_idx,lev_idx)=
c     $                    u0u(plr_idx,lev_idx)/wvl_dlt(bnd_idx)
c                  enddo         ! end loop over lev
c               enddo            ! end loop over plr
c            enddo               ! end loop over azimuthal angles
         endif                  ! end if saving full intensity arrays

 999     continue

      enddo                     ! end loop over bnd

c     End section 2: Main computation loop over all bands
c     Begin section 3: Postprocessing
c     computation of diagnostic arrays from output storage arrays
c     summation of spectral diagnostic arrays to integrated arrays
c     reduction of integrated arrays to diagnostic scalars
c     netCDF output

      if (dbg_lvl > dbg_off) write (6,'(/)')

      if (sv_cmp_plr_ngl) then
c     Computational angles are returned in umu array
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
         trn_spc_atm_NO2(bnd_idx)=exp(-min(odxc_spc_NO2(bnd_idx),25.0))
         trn_spc_atm_O2(bnd_idx)=exp(-min(odxc_spc_O2(bnd_idx),25.0))
         trn_spc_atm_O2N2(bnd_idx)=exp(-min(odxc_spc_O2N2(bnd_idx),25.0))
         trn_spc_atm_O2O2(bnd_idx)=exp(-min(odxc_spc_O2O2(bnd_idx),25.0))
         trn_spc_atm_O3(bnd_idx)=exp(-min(odxc_spc_O3(bnd_idx),25.0))
         trn_spc_atm_OH(bnd_idx)=exp(-min(odxc_spc_OH(bnd_idx),25.0))
         trn_spc_atm_CH4(bnd_idx)=exp(-min(odxc_spc_CH4(bnd_idx),25.0))
         trn_spc_atm_N2O(bnd_idx)=exp(-min(odxc_spc_N2O(bnd_idx),25.0))
         trn_spc_atm_Ray(bnd_idx)=exp(-min(odxc_spc_Ray(bnd_idx),25.0))
         trn_spc_atm_aer(bnd_idx)=exp(-min(odxc_spc_aer(bnd_idx),25.0))
         trn_spc_atm_bga(bnd_idx)=exp(-min(odxc_spc_bga(bnd_idx),25.0))
         trn_spc_atm_ice(bnd_idx)=exp(-min(odxc_spc_ice(bnd_idx),25.0))
         trn_spc_atm_lqd(bnd_idx)=exp(-min(odxc_spc_lqd(bnd_idx),25.0))
      enddo                     ! end loop over bnd
      do bnd_idx=1,bnd_nbr
         flx_spc_dwn_TOA(bnd_idx)=flx_spc_dwn(bnd_idx,1)
         flx_spc_dwn_sfc(bnd_idx)=flx_spc_dwn(bnd_idx,levp_nbr)

c     Compute absorbed spectral fluxes
         flx_spc_abs_SAS(bnd_idx)=flx_spc_net(bnd_idx,1)
         flx_spc_abs_sfc(bnd_idx)=flx_spc_net(bnd_idx,levp_nbr)
         flx_spc_abs_atm(bnd_idx)=
     $        flx_spc_abs_SAS(bnd_idx)-flx_spc_abs_sfc(bnd_idx)

         do lev_idx=1,lev_nbr
            flx_spc_abs(bnd_idx,lev_idx)=
     $           max(0.,
     $           flx_spc_net(bnd_idx,lev_idx)-flx_spc_net(bnd_idx,lev_idx+1))
         enddo                  ! end loop over lev
      enddo                     ! end loop over bnd

c     Compute running sum of atmospheric absorption
c     NB: flx_abs_atm_rdr holds total atmospheric absorption
c     occuring redward (at longer wavelengths) of a given band. 
c     This quantity, at certain wavelengths, can be backed out of 
c     fractional spectral instruments, like Valero's FSBR.
      flx_abs_atm_rdr(1)=flx_spc_abs_atm(1)*wvl_dlt(1)
      do bnd_idx=2,bnd_nbr
         flx_abs_atm_rdr(bnd_idx)=flx_abs_atm_rdr(bnd_idx-1)+
     $        flx_spc_abs_atm(bnd_idx)*wvl_dlt(bnd_idx)
      enddo                     ! end loop over bnd

c     Instrument fluxes are a paradigm useful for any fractional band measurement 
c     Therefore we compute all instrument diagnostics in one location
c     The computation of instrument diagnostics is perfectly analogous to 
c     broadband (bb) diagnostics; reuse code by changing bb to nst.
c     All that need be changed to adapt to a new instrument are the
c     lower and upper limits of the (rectangular window) bandpass. 
c     Initialize basic instrument quantities which will be incremented
      call vec_set(flx_nst_dwn,levp_nbr,0.0) ! [W m-2]
      call vec_set(flx_nst_up,levp_nbr,0.0) ! [W m-2]
      do bnd_idx=1,bnd_nbr
         do lev_idx=1,levp_nbr
            flx_nst_dwn(lev_idx)=flx_nst_dwn(lev_idx)+
     $           flx_spc_dwn(bnd_idx,lev_idx)*wvl_dlt(bnd_idx)*nst_SRF(bnd_idx)
            flx_nst_up(lev_idx)=flx_nst_up(lev_idx)+
     $           flx_spc_up(bnd_idx,lev_idx)*wvl_dlt(bnd_idx)*nst_SRF(bnd_idx)
         enddo                  ! end loop over lev
      enddo                     ! end loop over bnd
      do lev_idx=1,levp_nbr
         flx_nst_net(lev_idx)=flx_nst_dwn(lev_idx)-flx_nst_up(lev_idx)
      enddo                     ! end loop over lev
      do lev_idx=1,lev_nbr
         flx_nst_abs(lev_idx)=flx_nst_net(lev_idx)-flx_nst_net(lev_idx+1)
      enddo                  ! end loop over lev
      if (flx_nst_dwn(1) > 0.0) then 
         abs_nst_SAS=flx_nst_net(1)/flx_nst_dwn(1)
         abs_nst_atm=(flx_nst_net(1)-flx_nst_net(levp_nbr))/flx_nst_dwn(1)
         abs_nst_sfc=flx_nst_net(levp_nbr)/flx_nst_dwn(1)
         rfl_nst_SAS=flx_nst_up(1)/flx_nst_dwn(1)
         rfl_nst_sfc=flx_nst_up(levp_nbr)/flx_nst_dwn(levp_nbr)
         trn_nst_atm=flx_nst_dwn(levp_nbr)/flx_nst_dwn(1)
      else
         abs_nst_SAS=0.0
         abs_nst_atm=0.0
         abs_nst_sfc=0.0
         rfl_nst_SAS=0.0
         rfl_nst_sfc=0.0
         trn_nst_atm=0.0
      endif
      flx_nst_abs_atm=flx_nst_net(1)-flx_nst_net(levp_nbr)
      flx_nst_abs_sfc=flx_nst_net(levp_nbr)
      flx_nst_abs_ttl=flx_nst_net(1)
      flx_nst_dwn_TOA=flx_nst_dwn(1)
      flx_nst_dwn_sfc=flx_nst_dwn(levp_nbr)
c     End instrument computations

c     Define system transmittance, reflectance, and absorptance.
c     Because these quantities are all normalized by the insolation,
c     which may be zero, these definitions must be made in a 
c     conditional clause.
      do bnd_idx=1,bnd_nbr
         if (flx_spc_dwn_TOA(bnd_idx) > 0.0) then
            trn_spc_atm_ttl(bnd_idx)=
     $           flx_spc_dwn_sfc(bnd_idx)/
     $           flx_spc_dwn_TOA(bnd_idx)
            rfl_spc_SAS(bnd_idx)=
     $           flx_spc_up(bnd_idx,1)/
     $           flx_spc_dwn_TOA(bnd_idx)

c     Layer absorptance is absorbed flux in a layer normalized by total flux entering layer. 
c     Define absorptance so surface + atmospheric absorptances sum to total SAS absorptance,
c     i.e., as fraction of solar insolation absorbed by atmosphere, surface, and SAS, respectively.
            abs_spc_SAS(bnd_idx)=
     $           flx_spc_net(bnd_idx,1)/flx_spc_dwn_TOA(bnd_idx)
            abs_spc_sfc(bnd_idx)=
     $           flx_spc_net(bnd_idx,levp_nbr)/flx_spc_dwn_TOA(bnd_idx)
         else
            trn_spc_atm_ttl(bnd_idx)=0.0
            rfl_spc_SAS(bnd_idx)=0.0
            abs_spc_SAS(bnd_idx)=0.0
            abs_spc_sfc(bnd_idx)=0.0
         endif
      enddo                     ! end loop over bnd
      do bnd_idx=1,bnd_nbr
         abs_spc_atm(bnd_idx)=abs_spc_SAS(bnd_idx)-abs_spc_sfc(bnd_idx)
      enddo                     ! end loop over bnd

      do bnd_idx=1,bnd_nbr
c     Compute actinic fluxes
         nrg_pht(bnd_idx)=Planck*speed_of_light/wvl(bnd_idx) ! [J pht-1]
         flx_spc_pht_dwn_sfc(bnd_idx)=flx_spc_dwn_sfc(bnd_idx)/nrg_pht(bnd_idx) ! [pht m-2 s-1 m-1]
         do lev_idx=1,lev_nbr
            flx_spc_act=4.0*pi*ntn_spc_mean(bnd_idx,lev_idx) ! [W m-2 m-1 sr-1] --> [W m-2 m-1]
            flx_spc_act_pht=flx_spc_act/nrg_pht(bnd_idx) ! [W m-2 m-1] --> [pht m-2 s-1 m-1]
            j_spc_NO2=          ! [s-1 m-1]
     $           flx_spc_act_pht*abs_xsx_NO2(bnd_idx)*qnt_yld_NO2(bnd_idx)
            j_NO2(lev_idx)=     ! [s-1 m-1] --> [s-1] 
     $           j_NO2(lev_idx)+j_spc_NO2*wvl_dlt(bnd_idx)
c     j_NO2 refers to photolysis of NO2 into O(3P) + NO
         enddo                  ! end loop over lev
         flx_spc_act_pht_TOA(bnd_idx)=4.0*pi*ntn_spc_mean(bnd_idx,1)/nrg_pht(bnd_idx) ! [pht m-2 s-1 m-1]
         flx_spc_act_pht_sfc(bnd_idx)=flx_spc_act_pht ! [pht m-2 s-1 m-1]
         j_spc_NO2_sfc(bnd_idx)=j_spc_NO2 ! [s-1 m-1]

c     Accumulate spectral fluxes into broadband arrays
         do lev_idx=1,lev_nbr
            ntn_bb_mean(lev_idx)=ntn_bb_mean(lev_idx)+
     $           ntn_spc_mean(bnd_idx,lev_idx)*
     $           wvl_dlt(bnd_idx)
         enddo                  ! end loop over lev
         do lev_idx=1,levp_nbr
            flx_bb_dwn_drc(lev_idx)=flx_bb_dwn_drc(lev_idx)+
     $           flx_spc_dwn_drc(bnd_idx,lev_idx)*
     $           wvl_dlt(bnd_idx)
            flx_bb_dwn_dff(lev_idx)=flx_bb_dwn_dff(lev_idx)+
     $           flx_spc_dwn_dff(bnd_idx,lev_idx)*
     $           wvl_dlt(bnd_idx)
            flx_bb_up(lev_idx)=flx_bb_up(lev_idx)+
     $           flx_spc_up(bnd_idx,lev_idx)*
     $           wvl_dlt(bnd_idx)
         enddo                  ! end loop over lev
      enddo                  ! end loop over bnd

c     Process broadband fluxes
      do lev_idx=1,levp_nbr
         flx_bb_dwn(lev_idx)=flx_bb_dwn_drc(lev_idx)+flx_bb_dwn_dff(lev_idx)
         flx_bb_net(lev_idx)=flx_bb_dwn(lev_idx)-flx_bb_up(lev_idx)
      enddo                  ! end loop over lev

c     Compute scalar diagnostics
      abs_bb_SAS=flx_bb_net(1)/flx_bb_dwn(1)
      abs_bb_atm=(flx_bb_net(1)-flx_bb_net(levp_nbr))/flx_bb_dwn(1)
      abs_bb_sfc=flx_bb_net(levp_nbr)/flx_bb_dwn(1)
      rfl_bb_SAS=flx_bb_up(1)/flx_bb_dwn(1)
      rfl_bb_sfc=flx_bb_up(levp_nbr)/flx_bb_dwn(levp_nbr)
      trn_bb_atm=flx_bb_dwn(levp_nbr)/flx_bb_dwn(1)
      flx_bb_abs_ttl=flx_bb_net(1)
      flx_bb_abs_sfc=flx_bb_net(levp_nbr)
      flx_bb_abs_atm=flx_bb_net(1)-flx_bb_net(levp_nbr)
      flx_bb_dwn_TOA=flx_bb_dwn(1)
      flx_bb_dwn_sfc=flx_bb_dwn(levp_nbr)

      do lev_idx=1,lev_nbr
         flx_bb_abs(lev_idx)=flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1)
         htg_rate_bb(lev_idx)=
     $        (flx_bb_net(lev_idx)-flx_bb_net(lev_idx+1))*grv(lev_idx)/
     $        (spc_heat_mst_air(lev_idx)*prs_dlt(lev_idx))
      enddo                  ! end loop over lev
c      write (6,'(i4,a,es8.1)') bnd_idx,' flurt flx_spc_dwn(1690,levp_nbr) =  ',flx_spc_dwn(1690,levp_nbr)

c     Begin netCDF output routines
      rcd=rcd+nf90_create(fl_out,nf90_clobber,nc_id)
      rcd=rcd+nf90_set_fill(nc_id,nf90_nofill,fll_mode_old)

c     Define dimension IDs
      rcd=rcd+nf90_def_dim(nc_id,'azi',azi_nbr,azi_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'grd',bnd_nbr+1,grd_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'lev',lev_nbr,lev_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'levp',levp_nbr,levp_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'plr',plr_nbr,plr_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'tau',tau_nbr,tau_dim_id)

c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_bnd_lev(1)=bnd_dim_id
      dim_bnd_lev(2)=lev_dim_id

      dim_bnd_levp(1)=bnd_dim_id
      dim_bnd_levp(2)=levp_dim_id

      dim_plr_levp(1)=plr_dim_id
      dim_plr_levp(2)=levp_dim_id

      dim_plr_bnd(1)=plr_dim_id
      dim_plr_bnd(2)=bnd_dim_id

      dim_plr_bnd_levp(1)=plr_dim_id
      dim_plr_bnd_levp(2)=bnd_dim_id
      dim_plr_bnd_levp(3)=levp_dim_id

      dim_azi_plr_levp(1)=azi_dim_id
      dim_azi_plr_levp(2)=plr_dim_id
      dim_azi_plr_levp(3)=levp_dim_id

      dim_azi_plr_bnd_levp(1)=azi_dim_id
      dim_azi_plr_bnd_levp(2)=plr_dim_id
      dim_azi_plr_bnd_levp(3)=bnd_dim_id
      dim_azi_plr_bnd_levp(4)=levp_dim_id

c     Variable definitions
      rcd=rcd+nf90_def_var(nc_id,'abs_bb_SAS',nf90_float,abs_bb_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_bb_atm',nf90_float,abs_bb_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_bb_sfc',nf90_float,abs_bb_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_nst_SAS',nf90_float,abs_nst_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_nst_atm',nf90_float,abs_nst_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_nst_sfc',nf90_float,abs_nst_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_spc_SAS',nf90_float,bnd_dim_id,abs_spc_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_spc_atm',nf90_float,bnd_dim_id,abs_spc_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_spc_sfc',nf90_float,bnd_dim_id,abs_spc_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'alb_sfc',nf90_float,alb_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'azi',nf90_float,azi_dim_id,azi_id)
      rcd=rcd+nf90_def_var(nc_id,'azi_dgr',nf90_float,azi_dim_id,azi_dgr_id)
      rcd=rcd+nf90_def_var(nc_id,'bnd',nf90_float,bnd_dim_id,bnd_id)
      rcd=rcd+nf90_def_var(nc_id,'bnd_var',nf90_float,bnd_dim_id,bnd_var_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_abs_atm_rdr',nf90_float,bnd_dim_id,flx_abs_atm_rdr_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_abs',nf90_float,lev_dim_id,flx_bb_abs_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_abs_atm',nf90_float,flx_bb_abs_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_abs_sfc',nf90_float,flx_bb_abs_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_abs_ttl',nf90_float,flx_bb_abs_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_dwn',nf90_float,levp_dim_id,flx_bb_dwn_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_dwn_TOA',nf90_float,flx_bb_dwn_TOA_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_dwn_dff',nf90_float,levp_dim_id,flx_bb_dwn_dff_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_dwn_drc',nf90_float,levp_dim_id,flx_bb_dwn_drc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_dwn_sfc',nf90_float,flx_bb_dwn_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_net',nf90_float,levp_dim_id,flx_bb_net_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_bb_up',nf90_float,levp_dim_id,flx_bb_upw_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_abs',nf90_float,lev_dim_id,flx_nst_abs_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_abs_atm',nf90_float,flx_nst_abs_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_abs_sfc',nf90_float,flx_nst_abs_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_abs_ttl',nf90_float,flx_nst_abs_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_dwn',nf90_float,levp_dim_id,flx_nst_dwn_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_dwn_TOA',nf90_float,flx_nst_dwn_TOA_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_dwn_sfc',nf90_float,flx_nst_dwn_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_net',nf90_float,levp_dim_id,flx_nst_net_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_nst_up',nf90_float,levp_dim_id,flx_nst_upw_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dim_id,flx_slr_frc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_abs',nf90_float,dim_bnd_lev,flx_spc_abs_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_abs_SAS',nf90_float,bnd_dim_id,flx_spc_abs_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_abs_atm',nf90_float,bnd_dim_id,flx_spc_abs_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_abs_sfc',nf90_float,bnd_dim_id,flx_spc_abs_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn',nf90_float,dim_bnd_levp,flx_spc_dwn_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_dff',nf90_float,dim_bnd_levp,flx_spc_dwn_dff_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_drc',nf90_float,dim_bnd_levp,flx_spc_dwn_drc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_sfc',nf90_float,bnd_dim_id,flx_spc_dwn_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_up',nf90_float,dim_bnd_levp,flx_spc_upw_id)
      rcd=rcd+nf90_def_var(nc_id,'frc_ice_ttl',nf90_float,frc_ice_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'htg_rate_bb',nf90_float,lev_dim_id,htg_rate_bb_id)
      rcd=rcd+nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dim_id,nrg_pht_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_act_pht_TOA',nf90_float,bnd_dim_id,flx_spc_act_pht_TOA_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_act_pht_sfc',nf90_float,bnd_dim_id,flx_spc_act_pht_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_spc_pht_dwn_sfc',nf90_float,bnd_dim_id,flx_spc_pht_dwn_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'j_NO2',nf90_float,lev_dim_id,j_NO2_id)
      rcd=rcd+nf90_def_var(nc_id,'j_spc_NO2_sfc',nf90_float,bnd_dim_id,j_spc_NO2_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'lat_dgr',nf90_double,lat_dgr_id)
      rcd=rcd+nf90_def_var(nc_id,'lcl_time_hr',nf90_double,lcl_time_hr_id)
      rcd=rcd+nf90_def_var(nc_id,'lcl_yr_day',nf90_double,lcl_yr_day_id)
      rcd=rcd+nf90_def_var(nc_id,'lev',nf90_float,lev_dim_id,lev_id)
      rcd=rcd+nf90_def_var(nc_id,'levp',nf90_float,levp_dim_id,levp_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_CWP',nf90_float,mpc_CWP_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_bb_aa',nf90_float,dim_plr_levp,ntn_bb_aa_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_bb_mean',nf90_float,lev_dim_id,ntn_bb_mean_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa_ndr',nf90_float,dim_bnd_levp,ntn_spc_aa_ndr_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa_ndr_sfc',nf90_float,bnd_dim_id,ntn_spc_aa_ndr_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa_sfc',nf90_float,dim_plr_bnd,ntn_spc_aa_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa_zen',nf90_float,dim_bnd_levp,ntn_spc_aa_zen_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa_zen_sfc',nf90_float,bnd_dim_id,ntn_spc_aa_zen_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_chn',nf90_float,dim_azi_plr_levp,ntn_spc_chn_id)
      rcd=rcd+nf90_def_var(nc_id,'ntn_spc_mean',nf90_float,dim_bnd_lev,ntn_spc_mean_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_obs_aer',nf90_float,odxc_obs_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_obs_bga',nf90_float,odxc_obs_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_CO2',nf90_float,bnd_dim_id,odxc_spc_CO2_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_H2O',nf90_float,bnd_dim_id,odxc_spc_H2O_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_NO2',nf90_float,bnd_dim_id,odxc_spc_NO2_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_O2',nf90_float,bnd_dim_id,odxc_spc_O2_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_O2N2',nf90_float,bnd_dim_id,odxc_spc_O2N2_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_O2O2',nf90_float,bnd_dim_id,odxc_spc_O2O2_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_O3',nf90_float,bnd_dim_id,odxc_spc_O3_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_OH',nf90_float,bnd_dim_id,odxc_spc_OH_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_CH4',nf90_float,bnd_dim_id,odxc_spc_CH4_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_N2O',nf90_float,bnd_dim_id,odxc_spc_N2O_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_Ray',nf90_float,bnd_dim_id,odxc_spc_Ray_id)
      rcd=rcd+nf90_def_var(nc_id,'odac_spc_aer',nf90_float,bnd_dim_id,odac_spc_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odac_spc_bga',nf90_float,bnd_dim_id,odac_spc_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'odac_spc_ice',nf90_float,bnd_dim_id,odac_spc_ice_id)
      rcd=rcd+nf90_def_var(nc_id,'odac_spc_lqd',nf90_float,bnd_dim_id,odac_spc_lqd_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_aer',nf90_float,bnd_dim_id,odxc_spc_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_bga',nf90_float,bnd_dim_id,odxc_spc_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_H2OH2O',nf90_float,bnd_dim_id,odxc_spc_H2OH2O_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_ice',nf90_float,bnd_dim_id,odxc_spc_ice_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_lqd',nf90_float,bnd_dim_id,odxc_spc_lqd_id)
      rcd=rcd+nf90_def_var(nc_id,'odxc_spc_ttl',nf90_float,bnd_dim_id,odxc_spc_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'odal_obs_aer',nf90_float,lev_dim_id,odal_obs_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odal_obs_bga',nf90_float,lev_dim_id,odal_obs_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'odsl_obs_aer',nf90_float,lev_dim_id,odsl_obs_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odsl_obs_bga',nf90_float,lev_dim_id,odsl_obs_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'odxl_obs_aer',nf90_float,lev_dim_id,odxl_obs_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'odxl_obs_bga',nf90_float,lev_dim_id,odxl_obs_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'plr',nf90_float,plr_dim_id,plr_id)
      rcd=rcd+nf90_def_var(nc_id,'plr_cos',nf90_float,plr_dim_id,plr_cos_id)
      rcd=rcd+nf90_def_var(nc_id,'plr_dgr',nf90_float,plr_dim_id,plr_dgr_id)
      rcd=rcd+nf90_def_var(nc_id,'rfl_bb_SAS',nf90_float,rfl_bb_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'rfl_bb_sfc',nf90_float,rfl_bb_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'rfl_nst_SAS',nf90_float,rfl_nst_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'rfl_nst_sfc',nf90_float,rfl_nst_sfc_id)
      rcd=rcd+nf90_def_var(nc_id,'rfl_spc_SAS',nf90_float,bnd_dim_id,rfl_spc_SAS_id)
      rcd=rcd+nf90_def_var(nc_id,'slr_zen_ngl_cos',nf90_double,slr_zen_ngl_cos_id)
      rcd=rcd+nf90_def_var(nc_id,'tpt',nf90_float,lev_dim_id,tpt_id)
      rcd=rcd+nf90_def_var(nc_id,'tpt_ntf',nf90_float,levp_dim_id,tpt_ntf_id)
      rcd=rcd+nf90_def_var(nc_id,'tau',nf90_float,tau_dim_id,tau_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_prs',nf90_float,tau_dim_id,tau_prs_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_bb_atm',nf90_float,trn_bb_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_nst_atm',nf90_float,trn_nst_atm_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_CO2',nf90_float,bnd_dim_id,trn_spc_atm_CO2_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_H2O',nf90_float,bnd_dim_id,trn_spc_atm_H2O_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_NO2',nf90_float,bnd_dim_id,trn_spc_atm_NO2_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_O2',nf90_float,bnd_dim_id,trn_spc_atm_O2_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_O2N2',nf90_float,bnd_dim_id,trn_spc_atm_O2N2_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_O2O2',nf90_float,bnd_dim_id,trn_spc_atm_O2O2_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_O3',nf90_float,bnd_dim_id,trn_spc_atm_O3_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_OH',nf90_float,bnd_dim_id,trn_spc_atm_OH_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_CH4',nf90_float,bnd_dim_id,trn_spc_atm_CH4_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_N2O',nf90_float,bnd_dim_id,trn_spc_atm_N2O_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_Ray',nf90_float,bnd_dim_id,trn_spc_atm_Ray_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_aer',nf90_float,bnd_dim_id,trn_spc_atm_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_bga',nf90_float,bnd_dim_id,trn_spc_atm_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_H2OH2O',nf90_float,bnd_dim_id,trn_spc_atm_H2OH2O_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_ice',nf90_float,bnd_dim_id,trn_spc_atm_ice_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_lqd',nf90_float,bnd_dim_id,trn_spc_atm_lqd_id)
      rcd=rcd+nf90_def_var(nc_id,'trn_spc_atm_ttl',nf90_float,bnd_dim_id,trn_spc_atm_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dim_id,wvl_ctr_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dim_id,wvl_grd_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dim_id,wvl_max_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dim_id,wvl_min_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_obs_aer',nf90_float,wvl_obs_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_obs_bga',nf90_float,wvl_obs_bga_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dim_id,wvl_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'wvn_ctr',nf90_float,bnd_dim_id,wvn_ctr_id)
      rcd=rcd+nf90_def_var(nc_id,'wvn_grd',nf90_float,grd_dim_id,wvn_grd_id)
      rcd=rcd+nf90_def_var(nc_id,'wvn_max',nf90_float,bnd_dim_id,wvn_max_id)
      rcd=rcd+nf90_def_var(nc_id,'wvn_min',nf90_float,bnd_dim_id,wvn_min_id)
      rcd=rcd+nf90_def_var(nc_id,'wvn_dlt',nf90_float,bnd_dim_id,wvn_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'z',nf90_float,lev_dim_id,alt_id)
      rcd=rcd+nf90_def_var(nc_id,'alt_cld_btm',nf90_float,alt_cld_btm_id)
      rcd=rcd+nf90_def_var(nc_id,'alt_cld_thick',nf90_float,alt_cld_thick_id)
      rcd=rcd+nf90_def_var(nc_id,'alt_ntf',nf90_float,levp_dim_id,alt_ntf_id)

c     Add global attributes
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng)))

      rcd=rcd+nf90_put_att(nc_id,nf90_global,'aer_sng',aer_sng(1:ftn_strlen(aer_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'azi_sng',azi_sng(1:ftn_strlen(azi_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'bga_sng',bga_sng(1:ftn_strlen(bga_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'opt_dep_sng',opt_dep_sng(1:ftn_strlen(opt_dep_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'plr_sng',plr_sng(1:ftn_strlen(plr_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'prf_sng',prf_sng(1:ftn_strlen(prf_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'str_sng',str_sng(1:ftn_strlen(str_sng)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_CO2',stt_CO2(1:ftn_strlen(stt_CO2)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_H2O',stt_H2O(1:ftn_strlen(stt_H2O)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_H2OH2O',stt_H2OH2O(1:ftn_strlen(stt_H2OH2O)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_Herzberg',stt_Herzberg(1:ftn_strlen(stt_Herzberg)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_NO2',stt_NO2(1:ftn_strlen(stt_NO2)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_O2',stt_O2(1:ftn_strlen(stt_O2)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_O2N2',stt_O2N2(1:ftn_strlen(stt_O2N2)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_O2O2',stt_O2O2(1:ftn_strlen(stt_O2O2)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_O3',stt_O3(1:ftn_strlen(stt_O3)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_OH',stt_OH(1:ftn_strlen(stt_OH)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_CH4',stt_CH4(1:ftn_strlen(stt_CH4)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_N2O',stt_N2O(1:ftn_strlen(stt_N2O)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_Planck',stt_Planck(1:ftn_strlen(stt_Planck)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_Rayleigh',stt_Rayleigh(1:ftn_strlen(stt_Rayleigh)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_aer',stt_aer(1:ftn_strlen(stt_aer)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_bga',stt_bga(1:ftn_strlen(stt_bga)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_flt_nst',stt_flt_nst(1:ftn_strlen(stt_flt_nst)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_ice',stt_ice(1:ftn_strlen(stt_ice)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_lqd',stt_lqd(1:ftn_strlen(stt_lqd)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_slr',stt_slr(1:ftn_strlen(stt_slr)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'stt_top_lvl',stt_top_lvl(1:ftn_strlen(stt_top_lvl)))

c     Add english text descriptions
c     First do any description commands which exceed 132 characters
      rcd=rcd+nf90_put_att(nc_id,abs_bb_SAS_id,'long_name','Broadband absorptance of surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,abs_bb_atm_id,'long_name','Broadband absorptance of surface')
      rcd=rcd+nf90_put_att(nc_id,abs_bb_sfc_id,'long_name','Broadband absorptance of atmosphere')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_SAS_id,'long_name','FSBR absorptance of surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_atm_id,'long_name','FSBR absorptance of surface')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_sfc_id,'long_name','FSBR absorptance of atmosphere')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_SAS_id,'long_name','Spectral absorptance of surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_atm_id,'long_name','Spectral absorptance of surface')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_sfc_id,'long_name','Spectral absorptance of atmosphere')
      rcd=rcd+nf90_put_att(nc_id,alb_sfc_id,'long_name','Specified Lambertian surface albedo')
      rcd=rcd+nf90_put_att(nc_id,alt_cld_btm_id,'long_name','Highest interface beneath all clouds in column')
      rcd=rcd+nf90_put_att(nc_id,alt_cld_thick_id,'long_name','Thickness of region containing all clouds')
      rcd=rcd+nf90_put_att(nc_id,alt_id,'long_name','Altitude')
      rcd=rcd+nf90_put_att(nc_id,alt_ntf_id,'long_name','Interface altitude')
      rcd=rcd+nf90_put_att(nc_id,azi_dgr_id,'long_name','Azimuthal angle (degrees)')
      rcd=rcd+nf90_put_att(nc_id,azi_id,'long_name','Azimuthal angle (radians)')
      rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Midpoint wavelength')
      rcd=rcd+nf90_put_att(nc_id,bnd_var_id,'long_name','Generic spectral variable')
      rcd=rcd+nf90_put_att(nc_id,flx_abs_atm_rdr_id,'long_name','Flux absorbed in atmosphere at longer wavelengths')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_atm_id,'long_name','Broadband flux absorbed by atmospheric column only')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_id,'long_name','Broadband flux absorbed by layer')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_sfc_id,'long_name','Broadband flux absorbed by surface only')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_ttl_id,'long_name','Broadband flux absorbed by surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_TOA_id,'long_name','Broadband incoming flux at TOA (total insolation)')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_dff_id,'long_name','Diffuse downwelling broadband flux')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_drc_id,'long_name','Direct downwelling broadband flux')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_id,'long_name','Total downwelling broadband flux (direct + diffuse)')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_sfc_id,'long_name','Broadband downwelling flux at surface')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_net_id,'long_name','Net broadband flux (downwelling - upwelling)')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_upw_id,'long_name','Upwelling broadband flux')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_atm_id,'long_name','FSBR flux absorbed by atmospheric column only')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_id,'long_name','FSBR flux absorbed by layer')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_sfc_id,'long_name','FSBR flux absorbed by surface only')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_ttl_id,'long_name','FSBR flux absorbed by surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_TOA_id,'long_name','FSBR incoming flux at TOA (total insolation)')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_id,'long_name','Total downwelling FSBR flux (direct + diffuse)')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_sfc_id,'long_name','FSBR downwelling flux at surface')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_net_id,'long_name','Net FSBR flux (downwelling - upwelling)')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_upw_id,'long_name','Upwelling FSBR flux')
      rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_SAS_id,'long_name','Spectral flux absorbed by surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_atm_id,'long_name','Spectral flux absorbed by atmospheric column only')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_id,'long_name','Spectral flux absorbed by layer')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_sfc_id,'long_name','Spectral flux absorbed by surface only')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_act_pht_TOA_id,'long_name','Spectral actinic photon flux at TOA')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_act_pht_sfc_id,'long_name','Spectral actinic photon flux at surface')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_dff_id,'long_name','Spectral diffuse downwelling flux')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_drc_id,'long_name','Spectral direct downwelling flux')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_id,'long_name','Spectral downwelling flux')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_sfc_id,'long_name','Spectral solar insolation at surface')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_pht_dwn_sfc_id,'long_name','Spectral photon flux downwelling at surface')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_upw_id,'long_name','Spectral upwelling flux')
      rcd=rcd+nf90_put_att(nc_id,frc_ice_ttl_id,'long_name','Fraction of column condensate that is ice')
      rcd=rcd+nf90_put_att(nc_id,htg_rate_bb_id,'long_name','Broadband heating rate')
      rcd=rcd+nf90_put_att(nc_id,j_NO2_id,'long_name','Photolysis rate for NO2 + hv --> O(3P) + NO')
      rcd=rcd+nf90_put_att(nc_id,j_spc_NO2_sfc_id,'long_name','Spectral photolysis rate at surface for NO2+hv --> O(3P)+NO')
      rcd=rcd+nf90_put_att(nc_id,lat_dgr_id,'long_name','Latitude (degrees)')
      rcd=rcd+nf90_put_att(nc_id,lcl_time_hr_id,'long_name','Local day hour')
      rcd=rcd+nf90_put_att(nc_id,lcl_yr_day_id,'long_name','Day of year in local time')
      rcd=rcd+nf90_put_att(nc_id,lev_id,'long_name','Layer pressure')
      rcd=rcd+nf90_put_att(nc_id,levp_id,'long_name','Interface pressure')
      rcd=rcd+nf90_put_att(nc_id,mpc_CWP_id,'long_name','Total column Condensed Water Path')
      rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center')
      rcd=rcd+nf90_put_att(nc_id,ntn_bb_aa_id,'long_name','Broadband azimuthally averaged intensity')
      rcd=rcd+nf90_put_att(nc_id,ntn_bb_mean_id,'long_name','Broadband mean intensity')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_ndr_id,'long_name','Spectral intensity of nadir radiation')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_ndr_sfc_id,'long_name','Spectral intensity of nadir radiation at surface')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_sfc_id,'long_name','Spectral intensity of radiation at surface')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_zen_id,'long_name','Spectral intensity of zenith radiation')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_zen_sfc_id,'long_name','Spectral intensity of zenith radiation at surface')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_chn_id,'long_name','Full spectral intensity of particular band')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_mean_id,'long_name','Spectral mean intensity')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_aer_id,'long_name','Aerosol absorption optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_bga_id,'long_name','Background aerosol absorption optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_ice_id,'long_name','Liquid water absorption optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_lqd_id,'long_name','Ice water absorption optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odal_obs_aer_id,'long_name','Layer aerosol absorption optical depth')
      rcd=rcd+nf90_put_att(nc_id,odal_obs_bga_id,'long_name','Layer background aerosol absorption optical depth')
      rcd=rcd+nf90_put_att(nc_id,odsl_obs_aer_id,'long_name','Layer aerosol scattering optical depth')
      rcd=rcd+nf90_put_att(nc_id,odsl_obs_bga_id,'long_name','Layer background aerosol scattering optical depth')
      rcd=rcd+nf90_put_att(nc_id,odxc_obs_aer_id,'long_name','Column aerosol extinction optical depth')
      rcd=rcd+nf90_put_att(nc_id,odxc_obs_bga_id,'long_name','Column background aerosol extinction optical depth')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_CH4_id,'long_name','CH4 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_CO2_id,'long_name','CO2 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_H2OH2O_id,'long_name','H2O dimer optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_H2O_id,'long_name','H2O optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_N2O_id,'long_name','N2O optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_NO2_id,'long_name','NO2 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2N2_id,'long_name','O2N2 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2O2_id,'long_name','O2O2 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2_id,'long_name','O2 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O3_id,'long_name','O3 optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_OH_id,'long_name','OH optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_Ray_id,'long_name','Rayleigh scattering optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_aer_id,'long_name','Aerosol extinction optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_bga_id,'long_name','Background aerosol extinction optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_ice_id,'long_name','Ice water extinction optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_lqd_id,'long_name','Liquid water extinction optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_ttl_id,'long_name','Total extinction optical depth to surface')
      rcd=rcd+nf90_put_att(nc_id,odxl_obs_aer_id,'long_name','Layer aerosol extinction optical depth')
      rcd=rcd+nf90_put_att(nc_id,odxl_obs_bga_id,'long_name','Layer background aerosol extinction optical depth')
      rcd=rcd+nf90_put_att(nc_id,plr_cos_id,'long_name','Cosine polar angle (degrees)')
      rcd=rcd+nf90_put_att(nc_id,plr_dgr_id,'long_name','Polar angle (degrees)')
      rcd=rcd+nf90_put_att(nc_id,plr_id,'long_name','Polar angle (radians)')
      rcd=rcd+nf90_put_att(nc_id,rfl_bb_SAS_id,'long_name','Broadband albedo of entire surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,rfl_bb_sfc_id,'long_name','Broadband albedo of surface')
      rcd=rcd+nf90_put_att(nc_id,rfl_nst_SAS_id,'long_name','FSBR albedo of entire surface-atmosphere system')
      rcd=rcd+nf90_put_att(nc_id,rfl_nst_sfc_id,'long_name','FSBR albedo of surface')
      rcd=rcd+nf90_put_att(nc_id,rfl_spc_SAS_id,'long_name','Spectral planetary flux reflectance')
      rcd=rcd+nf90_put_att(nc_id,slr_zen_ngl_cos_id,'long_name','Cosine solar zenith angle')
      rcd=rcd+nf90_put_att(nc_id,tau_id,'long_name','Optical level (optical depth)')
      rcd=rcd+nf90_put_att(nc_id,tau_prs_id,'long_name','Optical level (pressure)')
      rcd=rcd+nf90_put_att(nc_id,tpt_id,'long_name','Layer Temperature')
      rcd=rcd+nf90_put_att(nc_id,tpt_ntf_id,'long_name','Interface temperature')
      rcd=rcd+nf90_put_att(nc_id,trn_bb_atm_id,'long_name','Broadband transmission of atmospheric column')
      rcd=rcd+nf90_put_att(nc_id,trn_nst_atm_id,'long_name','FSBR transmission of atmospheric column')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_CH4_id,'long_name','Column transmission due to CH4 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_CO2_id,'long_name','Column transmission due to CO2 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_H2OH2O_id,'long_name','Column transmission due to H2O dimer absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_H2O_id,'long_name','Column transmission due to H2O absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_N2O_id,'long_name','Column transmission due to N2O absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_NO2_id,'long_name','Column transmission due to NO2 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2N2_id,'long_name','Column transmission due to O2-N2 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2O2_id,'long_name','Column transmission due to O2-O2 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2_id,'long_name','Column transmission due to O2 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O3_id,'long_name','Column transmission due to O3 absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_OH_id,'long_name','Column transmission due to OH absorption')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_Ray_id,'long_name','Column transmission due to Rayleigh scattering')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_aer_id,'long_name','Column transmission due to aerosol extinction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_bga_id,'long_name','Column transmission due to background aerosol extinction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_ice_id,'long_name','Column transmission due to ice extinction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_lqd_id,'long_name','Column transmission due to liquid extinction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_ttl_id,'long_name','Spectral flux transmission of entire column')
      rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'long_name','Midpoint wavelength in band')
      rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'long_name','Width of band')
      rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid')
      rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'long_name','Maximum wavelength in band')
      rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'long_name','Minimum wavelength in band')
      rcd=rcd+nf90_put_att(nc_id,wvl_obs_aer_id,'long_name','Wavelength of aerosol optical depth specification')
      rcd=rcd+nf90_put_att(nc_id,wvl_obs_bga_id,'long_name','Wavelength of background aerosol optical depth specification')
      rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'long_name','Midpoint wavenumber in band')
      rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth in wavenumbers')
      rcd=rcd+nf90_put_att(nc_id,wvn_grd_id,'long_name','Wavenumber grid')
      rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'long_name','Maximum wavenumber in band')
      rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'long_name','Minimum wavenumber in band')

c     Add units
      rcd=rcd+nf90_put_att(nc_id,abs_bb_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_bb_atm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_bb_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_atm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_nst_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_atm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_spc_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,alb_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,azi_dgr_id,'units','degree')
      rcd=rcd+nf90_put_att(nc_id,azi_id,'units','radian')
      rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,bnd_var_id,'units','unknown')
      rcd=rcd+nf90_put_att(nc_id,slr_zen_ngl_cos_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,flx_abs_atm_rdr_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_atm_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_sfc_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_abs_ttl_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_TOA_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_dff_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_drc_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_dwn_sfc_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_net_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_bb_upw_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_atm_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_sfc_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_abs_ttl_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_TOA_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_dwn_sfc_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_net_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_nst_upw_id,'units','watt meter-2')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_SAS_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_atm_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_abs_sfc_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_dff_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_drc_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_sfc_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_upw_id,'units','watt meter-2 meter-1')
      rcd=rcd+nf90_put_att(nc_id,frc_ice_ttl_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,htg_rate_bb_id,'units','kelvin second-1')
      rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_act_pht_TOA_id,'units','photon meter-2 second-1 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_act_pht_sfc_id,'units','photon meter-2 second-1 meter-1')
      rcd=rcd+nf90_put_att(nc_id,flx_spc_pht_dwn_sfc_id,'units','photon meter-2 second-1 meter-1')
      rcd=rcd+nf90_put_att(nc_id,j_NO2_id,'units','second-1')
      rcd=rcd+nf90_put_att(nc_id,j_spc_NO2_sfc_id,'units','second-1 meter-1')
      rcd=rcd+nf90_put_att(nc_id,lat_dgr_id,'units','degree')
      rcd=rcd+nf90_put_att(nc_id,lcl_time_hr_id,'units','hour')
      rcd=rcd+nf90_put_att(nc_id,lcl_yr_day_id,'units','day')
      rcd=rcd+nf90_put_att(nc_id,lev_id,'units','pascal')
      rcd=rcd+nf90_put_att(nc_id,levp_id,'units','pascal')
      rcd=rcd+nf90_put_att(nc_id,mpc_CWP_id,'units','kilogram meter-2')
      rcd=rcd+nf90_put_att(nc_id,ntn_bb_aa_id,'units','watt meter-2 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_bb_mean_id,'units','watt meter-2 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_ndr_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_ndr_sfc_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_sfc_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_zen_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_zen_sfc_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_chn_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,ntn_spc_mean_id,'units','watt meter-2 meter-1 sterradian-1')
      rcd=rcd+nf90_put_att(nc_id,odxc_obs_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_obs_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odal_obs_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odal_obs_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odsl_obs_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odsl_obs_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxl_obs_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxl_obs_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_CO2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_H2OH2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_H2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_NO2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O3_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2O2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_O2N2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_OH_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_CH4_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_N2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_Ray_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_ice_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odac_spc_lqd_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_ice_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_lqd_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,odxc_spc_ttl_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,plr_cos_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,plr_dgr_id,'units','degree')
      rcd=rcd+nf90_put_att(nc_id,plr_id,'units','radian')
      rcd=rcd+nf90_put_att(nc_id,rfl_bb_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,rfl_bb_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,rfl_nst_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,rfl_nst_sfc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,rfl_spc_SAS_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,tpt_id,'units','kelvin')
      rcd=rcd+nf90_put_att(nc_id,tpt_ntf_id,'units','kelvin')
      rcd=rcd+nf90_put_att(nc_id,tau_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,tau_prs_id,'units','pascal')
      rcd=rcd+nf90_put_att(nc_id,trn_bb_atm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_nst_atm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_CO2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_H2OH2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_H2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_NO2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O3_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2O2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_O2N2_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_OH_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_CH4_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_N2O_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_Ray_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_aer_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_bga_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_ice_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_lqd_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,trn_spc_atm_ttl_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_obs_aer_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_obs_bga_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1')
      rcd=rcd+nf90_put_att(nc_id,wvn_grd_id,'units','centimeter-1')
      rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1')
      rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1')
      rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1')
      rcd=rcd+nf90_put_att(nc_id,alt_cld_btm_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,alt_cld_thick_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,alt_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,alt_ntf_id,'units','meter')

c     All dimensions, variables, and attributes have been defined
      rcd=rcd+nf90_enddef(nc_id)

c     Write data
      rcd=rcd+nf90_put_var(nc_id,lat_dgr_id,lat_dgr)
      rcd=rcd+nf90_put_var(nc_id,lcl_time_hr_id,lcl_time_hr)
      rcd=rcd+nf90_put_var(nc_id,lcl_yr_day_id,lcl_yr_day)
      rcd=rcd+nf90_put_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos)
      rcd=rcd+nf90_put_var(nc_id,abs_bb_SAS_id,abs_bb_SAS)
      rcd=rcd+nf90_put_var(nc_id,abs_bb_atm_id,abs_bb_atm)
      rcd=rcd+nf90_put_var(nc_id,abs_bb_sfc_id,abs_bb_sfc)
      rcd=rcd+nf90_put_var(nc_id,abs_nst_SAS_id,abs_nst_SAS)
      rcd=rcd+nf90_put_var(nc_id,abs_nst_atm_id,abs_nst_atm)
      rcd=rcd+nf90_put_var(nc_id,abs_nst_sfc_id,abs_nst_sfc)
      rcd=rcd+nf90_put_var(nc_id,abs_spc_SAS_id,abs_spc_SAS)
      rcd=rcd+nf90_put_var(nc_id,abs_spc_atm_id,abs_spc_atm)
      rcd=rcd+nf90_put_var(nc_id,abs_spc_sfc_id,abs_spc_sfc)
      rcd=rcd+nf90_put_var(nc_id,alb_sfc_id,alb_sfc)
      rcd=rcd+nf90_put_var(nc_id,azi_dgr_id,azi_dgr)
      rcd=rcd+nf90_put_var(nc_id,azi_id,azi)
      rcd=rcd+nf90_put_var(nc_id,bnd_id,bnd)
      rcd=rcd+nf90_put_var(nc_id,bnd_var_id,bnd_var)
      rcd=rcd+nf90_put_var(nc_id,flx_abs_atm_rdr_id,flx_abs_atm_rdr)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_abs_atm_id,flx_bb_abs_atm)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_abs_id,flx_bb_abs)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_abs_sfc_id,flx_bb_abs_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_abs_ttl_id,flx_bb_abs_ttl)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_dwn_TOA_id,flx_bb_dwn_TOA)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_dwn_dff_id,flx_bb_dwn_dff)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_dwn_drc_id,flx_bb_dwn_drc)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_dwn_id,flx_bb_dwn)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_dwn_sfc_id,flx_bb_dwn_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_net_id,flx_bb_net)
      rcd=rcd+nf90_put_var(nc_id,flx_bb_upw_id,flx_bb_up)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_abs_atm_id,flx_nst_abs_atm)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_abs_id,flx_nst_abs)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_abs_sfc_id,flx_nst_abs_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_abs_ttl_id,flx_nst_abs_ttl)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_dwn_TOA_id,flx_nst_dwn_TOA)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_dwn_id,flx_nst_dwn)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_dwn_sfc_id,flx_nst_dwn_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_net_id,flx_nst_net)
      rcd=rcd+nf90_put_var(nc_id,flx_nst_upw_id,flx_nst_up)
      rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_abs_SAS_id,flx_spc_abs_SAS)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_abs_atm_id,flx_spc_abs_atm)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_abs_id,flx_spc_abs)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_abs_sfc_id,flx_spc_abs_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_dff_id,flx_spc_dwn_dff)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_drc_id,flx_spc_dwn_drc)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_id,flx_spc_dwn)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_sfc_id,flx_spc_dwn_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_upw_id,flx_spc_up)
      rcd=rcd+nf90_put_var(nc_id,frc_ice_ttl_id,frc_ice_ttl)
      rcd=rcd+nf90_put_var(nc_id,htg_rate_bb_id,htg_rate_bb)
      rcd=rcd+nf90_put_var(nc_id,nrg_pht_id,nrg_pht)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_act_pht_TOA_id,flx_spc_act_pht_TOA)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_act_pht_sfc_id,flx_spc_act_pht_sfc)
      rcd=rcd+nf90_put_var(nc_id,flx_spc_pht_dwn_sfc_id,flx_spc_pht_dwn_sfc)
      rcd=rcd+nf90_put_var(nc_id,j_NO2_id,j_NO2)
      rcd=rcd+nf90_put_var(nc_id,j_spc_NO2_sfc_id,j_spc_NO2_sfc)
      rcd=rcd+nf90_put_var(nc_id,lxxxxxev_id,lev)
      rcd=rcd+nf90_put_var(nc_id,levp_id,levp)
      rcd=rcd+nf90_put_var(nc_id,mpc_CWP_id,mpc_CWP)
      rcd=rcd+nf90_put_var(nc_id,ntn_bb_aa_id,ntn_bb_aa)
      rcd=rcd+nf90_put_var(nc_id,ntn_bb_mean_id,ntn_bb_mean)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_ndr_id,ntn_spc_aa_ndr)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_ndr_sfc_id,ntn_spc_aa_ndr_sfc)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_sfc_id,ntn_spc_aa_sfc)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_zen_id,ntn_spc_aa_zen)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_zen_sfc_id,ntn_spc_aa_zen_sfc)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_chn_id,ntn_spc_chn)
      rcd=rcd+nf90_put_var(nc_id,ntn_spc_mean_id,ntn_spc_mean)
      rcd=rcd+nf90_put_var(nc_id,odxc_obs_aer_id,odxc_obs_aer)
      rcd=rcd+nf90_put_var(nc_id,odxc_obs_bga_id,odxc_obs_bga)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_CO2_id,odxc_spc_CO2)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_H2O_id,odxc_spc_H2O)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_NO2_id,odxc_spc_NO2)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_O2N2_id,odxc_spc_O2N2)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_O2O2_id,odxc_spc_O2O2)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_O2_id,odxc_spc_O2)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_O3_id,odxc_spc_O3)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_OH_id,odxc_spc_OH)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_CH4_id,odxc_spc_CH4)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_N2O_id,odxc_spc_N2O)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_Ray_id,odxc_spc_Ray)
      rcd=rcd+nf90_put_var(nc_id,odac_spc_aer_id,odac_spc_aer)
      rcd=rcd+nf90_put_var(nc_id,odac_spc_bga_id,odac_spc_bga)
      rcd=rcd+nf90_put_var(nc_id,odac_spc_ice_id,odac_spc_ice)
      rcd=rcd+nf90_put_var(nc_id,odac_spc_lqd_id,odac_spc_lqd)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_aer_id,odxc_spc_aer)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_bga_id,odxc_spc_bga)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_H2OH2O_id,odxc_spc_H2OH2O)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_ice_id,odxc_spc_ice)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_lqd_id,odxc_spc_lqd)
      rcd=rcd+nf90_put_var(nc_id,odxc_spc_ttl_id,odxc_spc_ttl)
      rcd=rcd+nf90_put_var(nc_id,odal_obs_aer_id,odal_obs_aer)
      rcd=rcd+nf90_put_var(nc_id,odal_obs_bga_id,odal_obs_bga)
      rcd=rcd+nf90_put_var(nc_id,odsl_obs_aer_id,odsl_obs_aer)
      rcd=rcd+nf90_put_var(nc_id,odsl_obs_bga_id,odsl_obs_bga)
      rcd=rcd+nf90_put_var(nc_id,odxl_obs_aer_id,odxl_obs_aer)
      rcd=rcd+nf90_put_var(nc_id,odxl_obs_bga_id,odxl_obs_bga)
      rcd=rcd+nf90_put_var(nc_id,plr_cos_id,plr_cos)
      rcd=rcd+nf90_put_var(nc_id,plr_dgr_id,plr_dgr)
      rcd=rcd+nf90_put_var(nc_id,plr_id,plr)
      rcd=rcd+nf90_put_var(nc_id,rfl_bb_SAS_id,rfl_bb_SAS)
      rcd=rcd+nf90_put_var(nc_id,rfl_bb_sfc_id,rfl_bb_sfc)
      rcd=rcd+nf90_put_var(nc_id,rfl_nst_SAS_id,rfl_nst_SAS)
      rcd=rcd+nf90_put_var(nc_id,rfl_nst_sfc_id,rfl_nst_sfc)
      rcd=rcd+nf90_put_var(nc_id,rfl_spc_SAS_id,rfl_spc_SAS)
      rcd=rcd+nf90_put_var(nc_id,tpt_id,tpt)
      rcd=rcd+nf90_put_var(nc_id,tpt_ntf_id,tpt_ntf)
      rcd=rcd+nf90_put_var(nc_id,tau_id,tau)
      rcd=rcd+nf90_put_var(nc_id,tau_prs_id,tau_prs)
      rcd=rcd+nf90_put_var(nc_id,trn_bb_atm_id,trn_bb_atm)
      rcd=rcd+nf90_put_var(nc_id,trn_nst_atm_id,trn_nst_atm)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_CO2_id,trn_spc_atm_CO2)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_H2O_id,trn_spc_atm_H2O)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_NO2_id,trn_spc_atm_NO2)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_O2N2_id,trn_spc_atm_O2N2)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_O2O2_id,trn_spc_atm_O2O2)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_O2_id,trn_spc_atm_O2)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_O3_id,trn_spc_atm_O3)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_OH_id,trn_spc_atm_OH)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_CH4_id,trn_spc_atm_CH4)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_N2O_id,trn_spc_atm_N2O)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_Ray_id,trn_spc_atm_Ray)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_aer_id,trn_spc_atm_aer)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_bga_id,trn_spc_atm_bga)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_H2OH2O_id,trn_spc_atm_H2OH2O)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_ice_id,trn_spc_atm_ice)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_lqd_id,trn_spc_atm_lqd)
      rcd=rcd+nf90_put_var(nc_id,trn_spc_atm_ttl_id,trn_spc_atm_ttl)
      rcd=rcd+nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr)
      rcd=rcd+nf90_put_var(nc_id,wvl_grd_id,wvl_grd)
      rcd=rcd+nf90_put_var(nc_id,wvl_max_id,wvl_max)
      rcd=rcd+nf90_put_var(nc_id,wvl_min_id,wvl_min)
      rcd=rcd+nf90_put_var(nc_id,wvl_obs_aer_id,wvl_obs_aer)
      rcd=rcd+nf90_put_var(nc_id,wvl_obs_bga_id,wvl_obs_bga)
      rcd=rcd+nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt)
      rcd=rcd+nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr)
      rcd=rcd+nf90_put_var(nc_id,wvn_grd_id,wvn_grd)
      rcd=rcd+nf90_put_var(nc_id,wvn_max_id,wvn_max)
      rcd=rcd+nf90_put_var(nc_id,wvn_min_id,wvn_min)
      rcd=rcd+nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt)
      rcd=rcd+nf90_put_var(nc_id,alt_cld_btm_id,alt_cld_btm)
      rcd=rcd+nf90_put_var(nc_id,alt_cld_thick_id,alt_cld_thick)
      rcd=rcd+nf90_put_var(nc_id,alt_id,alt)
      rcd=rcd+nf90_put_var(nc_id,alt_ntf_id,alt_ntf)

      if (sv_ntn) then
c     $     ntn_spc_aa_id,
c     $     ntn_spc_aa(plr_nbr_max,bnd_nbr_max,levp_nbr_max),
c     rcd=rcd+nf90_redef(nc_id)
c     rcd=rcd+nf90_def_var(nc_id,'ntn_spc_aa',nf90_float,dim_plr_bnd_levp,    ntn_spc_aa_id)
c     rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_id,'long_name','Spectral intensity of radiation')
c     rcd=rcd+nf90_put_att(nc_id,ntn_spc_aa_id,'units','watt meter-2 meter-1 sterradian-1')
c     rcd=rcd+nf90_enddef(nc_id)
c     rcd=rcd+nf90_put_var(nc_id,ntn_spc_aa_id,ntn_spc_aa)
      endif                     ! end if saving full intensity arrays

      rcd=rcd+nf90_wrp_close(nc_id,fl_out,'Wrote output to ')

      if (dbg_lvl == dbg_crr) then
c     fxm:
c     ncks -v ntn_spc_chn -F -C -m -H ${DATA}/aca/nbm.nc | m
c     rcd=rcd+nf90_close(nc_id)
c     write (6,'(a16,1x,a)') 'Wrote results to',fl_out(1:ftn_strlen(fl_out))
         write (6,'(a,3i2)') 'azi_nbr, plr_nbr, levp_nbr = ',azi_nbr,plr_nbr,levp_nbr
         rcd=rcd+nf90_wrp_open(fl_out,nf90_nowrite,nc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'ntn_spc_chn',ntn_spc_chn_id)
         rcd=rcd+nf90_inquire_variable(nc_id,ntn_spc_chn_id,ndims=int_foo)
         write (6,'(a,i3)') 'nf90_inquire_variable(ntn_spc_chn) ndims = ',int_foo
         rcd=rcd+nf90_inquire_variable(nc_id,ntn_spc_chn_id,dimids=dim_id_vec_foo)
         char_foo='   '
         call ftn_strnul(char_foo)
         do lev_idx=1,int_foo
            rcd=rcd+nf90_inquire_dimension(nc_id,dim_id_vec_foo(lev_idx),name=char_foo)
            rcd=rcd+nf90_inquire_dimension(nc_id,dim_id_vec_foo(lev_idx),len=int_foo)
            write (6,'(3a,i3,a,i3)') 'nf90_inquire_dimension(',
     $           char_foo(1:ftn_strlen(char_foo)),'), len= ',int_foo,
     $           'id = ',dim_id_vec_foo(lev_idx)
         enddo                  ! end loop over lev
         write (6,'(a)') 'Zeroing ntn_spc_chn...'
         call vec_set(ntn_spc_chn,azi_nbr*plr_nbr*levp_nbr,0.0)
         write (6,'(a)') 'Reading ntn_spc_chn...'
         rcd=rcd+nf90_get_var(nc_id,ntn_spc_chn_id,ntn_spc_chn)
         write (6,'(a,es10.2)') 'ntn_spc_chn(1,4,2) = ',ntn_spc_chn(1,4,2) ! fxm
         rcd=rcd+nf90_wrp_close(nc_id,fl_out,'Ingested')
         goto 1000
      endif                     ! endif dbg
      
 1000 continue
      
      call exit(exit_status)
      end                       ! end nbm()
