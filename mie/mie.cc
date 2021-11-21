// $Id$ 

/* Purpose: Compute and store Mie scattering properties for various aerosols
   Weight properties by radiant flux, either TOA solar spectrum or Planck function
   High resolution template for developing coarse resolution aerosol parameterizations */

/* Copyright (C) 1997--2014 Charlie Zender
   You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text
   The original author of this software, Charlie Zender, wants to improve it
   with the help of your suggestions, improvements, bug-reports, and patches.
   Charlie Zender <zender at uci dot edu>
   Department of Earth System Science
   University of California, Irvine
   Irvine, CA 92697-3100 */

// Editing, Compilation:
// etags ~/sw/aer/*.F ~/sw/aer/*.h ~/sw/c++/*.cc ~/sw/c++/*.hh ~/sw/mie/*.cc ~/sw/mie/*.hh ~/sw/slr_spc/*.cc ~/sw/slr_spc/*.hh
/* cd ~/sw/mie;make OPTS=D VRS_SNG=3.4.2 mie;cd -
   cd ~/sw/mie;make OPTS=D NETCDF4=Y mie;cd -
   cd ~/sw/c++;make OMP=N OPTS=D NETCDF4=Y UDUNITS2=Y;cd ~/sw/mie;make OMP=N OPTS=D NETCDF4=Y UDUNITS2=Y
   cd ~/sw/c++;make OMP=N OPTS=D NETCDF4=N UDUNITS2=N;cd ~/sw/mie;make OMP=N OPTS=D NETCDF4=N UDUNITS2=N
   cd ~/sw/c++;make dst_cln;cd ~/nco/src/nco_c++;make -f Makefile.old dst_cln;cd ~/sw/mie;make dst_cln;make OPTS=X mie
   unset LM_LICENSE_FILE;${HOME}/bin/sh/cmp_chg.sh icc;cd ~/nco/src/nco_c++;make -f Makefile.old dst_cln all;cd ~/sw/c++;make dst_cln;cd ~/sw/mie;make dst_cln;make --jobs=1 OPTS=X OMP=Y mie 
   cd ~/sw/c++;make dst_cln;cd ~/nco/src/nco_c++;make -f Makefile.old dst_cln;cd ~/sw/mie;make dst_cln;make --jobs=1 OPTS=D mie */

/* Scripts and Documentation:
   Documentation of command line switches: ${HOME}/crr/psd.tex
   Scripts to compute optical and microphysical properties:
   Default aerosols and snow, aerosols for ARESE: ${HOME}/mie/mie_arese.sh
   Mineral dust aerosols for CAM/CCM/GSFC/MATCH dust model: ${HOME}/mie/psd.pl
   Liquid cloud optics and resonance studies: ${HOME}/mie/mie_rsn.sh
   Multi-component aerosols, snow: ${HOME}/mie/mie_mca.sh
   Scripts to plot output data:
   NCL: Optical properties: ${HOME}/anl/mie_xv.ncl
   IDL: Optical properties, Angstrom exponent: ${HOME}/idl/mie.pro:aer_gph()
   IDL: Particle size distributions, Settling velocity: ${HOME}/idl/mie.pro:psd_bch() 
   IDL: Scavenging coefficients: ${HOME}/idl/mie.pro:scv_bch()
   IDL: Phase function: ${HOME}/idl/mie.pro:phz_fnc(),aer_gph() */

/* Usage:
   mie -dbg --wvl_mnm=1.0 --wvl_mxm=2.0 --wvl_nbr=10 2>&1 &
   mie -dbg --wvl_grd=CAM_SW
   mie -dbg --wvl_mnm=1.0 --wvl_mxm=2.0 --wvl_nbr=100
   mie -dbg --wvl_mnm=0.3 --wvl_mxm=4.0 --wvl_nbr=100 --bnd_nbr=10
   mie -dbg --cmp_prt=dust_like --wvl_mnm=2.0 --wvl_mxm=3.0 --wvl_nbr=1 --bnd_nbr=1
   mie -dbg --sz_mnm=0.0 --sz_mxm=10.0 --sz_nbr=10 --bnd_nbr=1 --wvl_nbr=1 --sz_grd=lin
   mie -dbg --wvl_nbr=1 --bnd_nbr=1 --sz_nbr=1 --ngl_nbr=91 foo.nc
   mie -dbg --wvl_nbr=1 --bnd_nbr=1 --sz_nbr=1 foo.nc
   scp ~/sw/mie/mie.cc esmf.ess.uci.edu:mie

// Validations:
// BoH83 p. 482 test case for mie_sph_BoH83()
   mie -dbg --slf_tst=BoH83 --mie_slv=BoH83
// BoH83 p. 482 test case for mie_sph_coat_BoH83()
   mie -dbg --slf_tst=BoH83 --coat
// ViC98 p. 5 Tbl. 3 test case for EMAs
   mie -dbg --slf_tst=ViC98
// Planck function
   mie --dbg=3 --cpv_foo=1.0e-6
   mie --dbg=3 --tpt_bbd_wgt=300.0
   mie --wvl_mnm=3.0 --wvl_mxm=40.0 --wvl_nbr=1

// Nussenzweig: Use --tst=nsz to omit memory-intensive non-Mie routines
// BoH83 p. 300 Figure 11.7, 11.8 from CKK78a m = 1.33+1.0e-8i near X = 50.337
   mie --sz_mnm=4.999 --sz_mxm=5.001 --sz_nbr=1 --wvl_nbr=100 --wvl_mnm=0.06240872 --wvl_mxm=0.06241244 --cmp_prt=h2o_lqd --idx_rfr_prt="1.33+1.0e-8" ${DATA}/mie/mie.nc 2>&1 > ${DATA}/mie/foo.txt &

// SeaWiFS:
   mie --no_wrn_ntp --dbg=1 --spc_idx_sng=01 --cmp_prt=saharan_dust --psd_typ=lognormal --wvl_grd=reg --wvl_nbr=1 --wvl_mnm=0.860 --wvl_mxm=0.870 --bnd_nbr=2 --sz_grd=log --sz_mnm=0.05 --sz_mxm=0.5 --sz_nbr=200 --rds_nma=0.275 --gsd_anl=2.2 ${DATA}/dst/mie/aer_saharan_dust_01_086.nc

// Aerosol deposition and mobilization properties
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.01 --sz_mxm=1.0 --sz_nbr=200 --gsd_anl=2.2
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.25 --sz_mxm=100.0 --sz_nbr=200 --gsd_anl=2.0
// ZSK94 p. 816 Tbl 1 Bin 3 validation   
   ncks -C -F -d sz,38.29e-6 -v q_HNO3_gas,cnc_HNO3_gas,vmr_HNO3_gas,mfp_HNO3_air,dff_HNO3_air,vlc_mwb_HNO3,knd_nbr_HNO3_air,vnt_crc_aer,cnt_rgm_bll_crc_HNO3_aer,dff_HNO3_aer,mss_acm_cff_HNO3_aer,mss_upt_cff_HNO3_aer ${DATA}/mie/mie.nc
   ncks -C -F -v tpt_mdp,q_H2O_vpr,RH_lqd,RH_ice,qst_H2O_lqd,qst_H2O_ice,svp_H2O_lqd,svp_H2O_ice ${DATA}/mie/mie.nc
   ncks -C -F -v tpt_gnd,tpt_vgt,mno_lng_mbl,wnd_frc_mbl,flx_LW_upw_sfc,flx_ltn,flx_q_H2O,wnd_str_mrd,wnd_str_znl ${DATA}/mie/mie.nc
   ncks -C -F -d sz,38.29e-6 -v mss_frc_cly,flx_mss_hrz_slt_ttl,flx_mss_vrt_dst_ttl,dst_slt_flx_rat_ttl ${DATA}/mie/mie.nc
   ncks -C -F -d sz,38.29e-6 -v dns_mdp,mss_frc_cly,mss_frc_snd,frc_thr_ncr_wtr,frc_thr_ncr_drg,snw_hgt_lqd,snw_frc,vwc_sfc,gwc_sfc,dns_blk_sfc,dns_blk_dry,vwc_rel,vwc_thr,wnd_mdp,wnd_rfr_mbl,wnd_frc_mbl,wnd_frc_thr,wnd_frc_thr_slt,wnd_frc_thr_opt,wnd_rfr_thr_slt,rgh_mmn_mbl,hgt_zpd_mbl,flx_mss_hrz_slt_ttl,flx_mss_vrt_dst_ttl,dst_slt_flx_rat_ttl ${DATA}/mie/mie.nc
   ncks -C -F -d sz,38.29e-6 -v mss_frc_cly,mss_frc_snd,frc_thr_ncr_wtr,smp_sat,smp_sfc,smp_xpn_b,snw_hgt_lqd,vwc_dry,vwc_opt,vwc_rel,vwc_sat,vwc_sfc,gwc_sfc,dns_blk_sfc,dns_blk_dry,dns_prt_sfc,vwc_thr,cnd_trm_soi,trn_fsh_vpr_soi_atm ${DATA}/mie/mie.nc
   ncks -C -F -d sz,10.0e-6 -v mlmic,mss_wtr_rat,svp_fct_klv,svp_fct_slt,svp_fct_ttl,act_cff,mss_frc_solute ${DATA}/mie/mie.nc
   ncks -C -F -d sz,10.0e-6 -v ryn_nbr_grv,vlc_grv,vlc_stk,stk_crc,slp_crc ${DATA}/mie/mie.nc
   ncks -C -F -d sz,10.0e-6 -v dmt_pcp_nma,ryn_nbr_pcp,vlc_pcp,vlc_stk_pcp,tau_rlx,stk_nbr_rlt,cll_fsh ${DATA}/mie/mie.nc
   ncks -C -F -d sz,10.0e-6 -v vlc_grv,vlc_trb,vlc_dry,rss_lmn,rss_aer_mmn ${DATA}/mie/mie.nc
   ncks -C -F -d sz,38.29e-6 -v dns_mdp,dmt_ctr,dmt_slt_opt,vsc_knm_atm,vsc_dyn_atm,dns_prt,wnd_frc_thr,ryn_nbr,ryn_nbr_frc,ryn_nbr_frc_thr,ryn_nbr_frc_thr_prx,wnd_frc_thr_prx,wnd_frc_thr_slt ${DATA}/mie/mie.nc
   ncks -C -F -d sz,38.29e-6 -v dmt_ctr,wbl_scl,wbl_shp,rss_aer,str_shr,cff_xch_mmn,wnd_rfr_mbl,wnd_frc_mbl,pugtut,wnd_rfr_thr ${DATA}/mie/mie.nc
   ncks -C -F -v dns_mdp,mno_lng_mbl,prs_mdp,tpt_mdp,tpt_vrt,wnd_frc_mbl,hgt_mdp,wnd_mdp,prs_ntf,q_H2O_vpr,sfc_typ,rgh_mmn_mbl,hgt_zpd_mbl,wnd_rfr_mbl,tpt_ptn_mdp,tpt_ptn_vrt_mdp ${DATA}/mie/mie.nc
   ncks -C -F -v rss_lmn_vwr,shm_nbr_vwr,stk_nbr_vwr,vlc_grv_vwr ${DATA}/mie/mie.nc
   ncks -C -F -v cnc_spc_rsl,xsa_spc_rsl,vlm_spc_rsl ${DATA}/mie/mie.nc
   ncks -C -F -v dmt_nwr,dmt_swr,dmt_vwr,dmt_nmr,dmt_smr,dmt_vmr ${DATA}/mie/mie.nc

// Aerosol multimodal capability
// Verify one mode yields same answers as two modes with half the concentration
   mie --dbg=1 --psd=0.2986,2.0,1.0 --sz_mnm=0.05 --sz_mxm=5.0 --sz_nbr=400 ${DATA}/mie/mie.nc
   mie --dbg=1 --psd=0.2986,2.0,0.5 --psd=0.2986,2.0,0.5 --sz_mnm=0.05 --sz_mxm=5.0 --sz_nbr=400 ${DATA}/mie/mie2.nc
   mie --dbg=1 --psd=0.2986,2.0,0.5,mss --psd=0.2986,2.0,0.5,mss --sz_mnm=0.05 --sz_mxm=5.0 --sz_nbr=400 ${DATA}/mie/mie2.nc
// Construct bi-modal AERONET Bahrain distribution from DHE02 p. 596 Tbl. 1: Vcoarse/Vfine ~ 10 (15--5), Rcoarse = 2.54 (2.6--1.9) um, Rfine ~ 0.15 (0.12--0.16) um, sigmacoarse ~ 0.61, sigmafine ~ 0.42, Ccoarse=
   ncap2 -O -v -s 'rds_vma_crs=2.54;ln_gsd_crs=0.61;tau_1020=0.5;C_crs=-0.02+tau_1020*0.92;gsd_anl_crs=exp(ln_gsd_crs);rds_nma_crs=rds_vma_crs*exp(-3.0*(ln_gsd_crs^2));C_fn=0.02+tau_1020*0.10;rds_vma_fn=0.15;ln_gsd_fn=0.42;gsd_anl_fn=exp(ln_gsd_fn);rds_nma_fn=rds_vma_fn*exp(-3.0*(ln_gsd_fn^2));' ~/nco/data/in.nc ~/foo.nc;ncks ~/foo.nc
   mie --dbg=1 --psd=0.0884,1.52,0.5,mss --psd=0.8318,1.84,0.5,mss --sz_mnm=0.05 --sz_mxm=5.0 --sz_nbr=400 ~/foo_Bhr.nc

// Verify precipitation size distribution is identical to multimodal particle size distribution in concentration and fall speed
   mie --no_wrn_ntp -no_mie --cmp_prt=h2o_lqd --sz_grd=log --dsd_mnm=0.0018 --dsd_mxm=2000.0 --dsd_nbr=200 --gsd_pcp_anl=2.0 --dmt_pcp_nma=2.0 --sz_mnm=0.0009 --sz_mxm=1000.0 --sz_nbr=200 --psd=1.0,2.0,0.5 --psd=1.0,2.0,0.5 --gsd_anl=2.0
   ncks -C -F -d sz,1.0e-6 -d dsd_sz,2.0e-6 -v vlc_grv,vlc_pcp,dmt_ctr,dmt_pcp,cnc,cnc_pcp ${DATA}/mie/mie.nc
// Wet deposition validation against SeP97 p. 1020 Fig. 20.10, 20.11
   mie -dbg --no_wrn_ntp -no_mie --dmt_pcp_nma=200.0 --flx_vlm_pcp=0.277e-6 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.0009 --sz_mxm=1000.0 --sz_nbr=200
   ncks -C -F -d sz,1.0e-6 -v cll_fsh,cll_fsh_brn_dff,cll_fsh_ntc,cll_fsh_mpc,stc_fsh,clc_fsh,scv_cff,scv_cff_pcp_nrm,scv_cff_mss_avg,scv_cff_mss_avg_pcp_nrm,scv_cff_nbr_avg ${DATA}/mie/mie.nc
// Wet deposition validation against SeP97 p. 1023 Fig. 20.12
   mie -dbg --no_wrn_ntp -no_mie --dmt_pcp_nma=400.0 --flx_vlm_pcp=1.0e-6 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.0009 --sz_mxm=1000.0 --sz_nbr=200 --rds_nma=1.0 --gsd_anl=2.0
// Wet deposition for scavenging figures
   mie -dbg --no_wrn_ntp -no_mie --dmt_pcp_nma=1000.0 --flx_vlm_pcp=1.0e-6 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.0009 --sz_mxm=1000.0 --sz_nbr=200 --rds_nma=1.0 --gsd_anl=2.0
// Wet deposition for mineral dust aerosol
// Convective
   mie --no_wrn_ntp -no_mie --hrz --dbg=1 --spc_idx_sng=01 --dsd_mnm=1.0 --dsd_mxm=4000.0 --dsd_nbr=200 --gsd_pcp_anl=1.86 --dmt_pcp_nma=1000.0 --flx_vlm_pcp=2.77e-6 --cmp_prt=saharan_dust --psd_typ=lognormal --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 ${DATA}/mie/scv_cnv.nc 2>&1 | m
   mie --no_wrn_ntp -no_mie --hrz --dbg=1 --spc_idx_sng=01 --dsd_mnm=999.0 --dsd_mxm=1000.0 --dsd_nbr=1 --gsd_pcp_anl=1.01 --dmt_pcp_nma=1000.0 --flx_vlm_pcp=2.77e-6 --cmp_prt=saharan_dust --psd_typ=lognormal --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 ${DATA}/mie/scv_cnv_mono.nc 2>&1 | m
// Stratiform
   mie --no_wrn_ntp -no_mie --hrz --dbg=1 --spc_idx_sng=01 --dsd_mnm=1.0 --dsd_mxm=4000.0 --dsd_nbr=200 --gsd_pcp_anl=1.86 --dmt_pcp_nma=400.0 --flx_vlm_pcp=0.277e-6 --cmp_prt=saharan_dust --psd_typ=lognormal --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 ${DATA}/mie/scv_str.nc 2>&1 | m
   mie --no_wrn_ntp -no_mie --hrz --dbg=1 --spc_idx_sng=01 --dsd_mnm=399.0 --dsd_mxm=401.0 --dsd_nbr=1 --gsd_pcp_anl=1.01 --dmt_pcp_nma=400.0 --flx_vlm_pcp=0.277e-6 --cmp_prt=saharan_dust --psd_typ=lognormal --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 ${DATA}/mie/scv_str_mono.nc 2>&1 | m

// Dry deposition validation against SeP97 p. 971 Fig. 19.4
   mie -dbg -no_mie --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.25 --sz_mxm=100.0 --sz_nbr=200 --gsd_anl=2.0 --rgh_mmn_mbl=0.001 --wnd_znl_mdp=1.0 --flx_SW_net_gnd=0.0 ${DATA}/mie/dps.nc
   ncks -C -F -d sz,0.5e-6 -v vlc_grv,vlc_trb,vlc_dry,rss_lmn,rss_aer_mmn,rgh_mmn_mbl,mno_lng_mbl ${DATA}/mie/mie.nc
// Dry deposition validation against Seh84 p. 553 Fig. 12.3
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.005 --sz_mxm=100.0 --sz_nbr=200 --rgh_mmn_dps=0.1 --wnd_frc_dps=0.30 --wnd_znl_mdp=10.0 --mno_lng_dps=100.0 --hgt_mdp=100.0 ${DATA}/mie/dps_01.nc
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.005 --sz_mxm=100.0 --sz_nbr=200 --rgh_mmn_dps=0.01 --wnd_frc_dps=0.30 --wnd_znl_mdp=10.0 --mno_lng_dps=100.0 --hgt_mdp=100.0 ${DATA}/mie/dps_02.nc
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.005 --sz_mxm=100.0 --sz_nbr=200 --rgh_mmn_dps=0.001 --wnd_frc_dps=0.30 --wnd_znl_mdp=10.0 --mno_lng_dps=100.0 --hgt_mdp=100.0 ${DATA}/mie/dps_03.nc
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.005 --sz_mxm=100.0 --sz_nbr=200 --rgh_mmn_dps=0.0001 --wnd_frc_dps=0.30 --wnd_znl_mdp=10.0 --mno_lng_dps=100.0 --hgt_mdp=100.0 ${DATA}/mie/dps_04.nc
   mie -dbg -no_mie --cmp_prt=saharan_dust --sz_grd=log --sz_mnm=0.005 --sz_mxm=100.0 --sz_nbr=200 --rgh_mmn_dps=0.1 --wnd_frc_dps=0.30 --wnd_znl_mdp=10.0 --mno_lng_dps=100.0 --hgt_mdp=100.0 ${DATA}/mie/dps_tst.nc

// Clouds
   mie -dbg --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=200 --rds_nma=5.75 --gsd_anl=1.6 --wvl_mnm=0.25 --wvl_mxm=2.05 --wvl_nbr=18
   mie -dbg --spc_idx_sng=01 --cmp_prt=h2o_lqd --wvl_grd=CAM_SW --bnd_nbr=3 --sz_grd=log --sz_mnm=1.0 --sz_mxm=30.0 --sz_nbr=25 --rds_nma=5.75 --gsd_anl=1.6 ${DATA}/mie/aer_H2O_lqd_SW.nc
   mie -dbg --spc_idx_sng=01 --cmp_prt=h2o_lqd --wvl_grd=CAM_LW --bnd_SW_LW=4.5e-6 --bnd_nbr=10 --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=200 --rds_nma=5.75 --gsd_anl=1.6 ${DATA}/mie/aer_H2O_lqd_LW.nc
   mie -dbg --spc_idx_sng=01 --cmp_prt=h2o_lqd --wvl_grd=GSFC_LW --bnd_SW_LW=4.5e-6 --bnd_nbr=10 --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=200 --rds_nma=5.75 --gsd_anl=1.6 ${DATA}/mie/aer_H2O_lqd_GSFC_LW.nc
   ncks -C -F -v rds_swr,rds_swa,rds_nmr,rds_vmr,ext_cff_mss,bck_cff_mss,asm_prm,ss_alb,abs_cff_mss_bb_LW ${DATA}/mie/aer_H2O_lqd_LW.nc | m

// LW properties
   mie -dbg --spc_idx_sng=01 --wvl_mnm=9.0 --wvl_mxm=11.0 --wvl_nbr=1 --bnd_nbr=100 --sz_grd=log --sz_mnm=0.01 --sz_mxm=10.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 --idx_rfr_prt="1.55+0.9i" ${DATA}/mie/aer_saharan_dust_LW.nc
   mie -dbg --spc_idx_sng=01 --cmp_prt=saharan_dust --wvl_grd=CAM_LW --bnd_SW_LW=4.5e-6 --bnd_nbr=10 --sz_grd=log --sz_mnm=0.01 --sz_mxm=10.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/aer_saharan_dust_01_LW.nc &
   ncks -C -F -v abs_cff_mss,asm_prm,msv_1km ${DATA}/mie/aer_saharan_dust_01_LW.nc

// Shaocai Yu: rg=0.11 um, sigma(g)=1.71, N=150 cm(-3), m=1.50+0.051i), dry aerosol
   mie -dbg --wvl_grd=CAM_SW --bnd_nbr=3 --sz_grd=log --sz_mnm=0.001 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.11 --gsd_anl=1.71 --idx_rfr_prt="1.50+0.051i" --dns_prt=1860.0 ${DATA}/mie/aer_ncsu.nc
   mie -dbg --wvl_grd=CAM_SW --cmp_prt=lac_YZS00 --bnd_nbr=3 --sz_grd=log --sz_mnm=0.001 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.11 --gsd_anl=1.71 ${DATA}/mie/aer_YZS00.nc

// Refractive indices mineral dust aerosol DKS91 
   mie --dbg=3 -no_mie --cmp_prt=mineral_dust --sz_nbr=1 --wvl_nbr=1 --bnd_nbr=1

// Size distribution WKB96 Fig.2
   mie -dbg -no_mie --psd_typ=lognormal --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.0 --sz_nbr=100 --rds_nma=0.05 --gsd_anl=2.01 --wvl_nbr=1 --bnd_nbr=1

// Size distribution HaT74 Fig.7, p. 550
   mie -dbg -no_mie --psd_typ=gamma --sz_grd=log --sz_mnm=0.1 --sz_mxm=10.0 --sz_nbr=100 --rds_ffc_gmm=1.0 --var_ffc_gmm=0.01 --wvl_nbr=1 --bnd_nbr=1

// Size distribution HaT74 Fig.7, p. 550
   mie -dbg -no_mie --psd_typ=gamma --sz_grd=log --sz_mnm=0.1 --sz_mxm=10.0 --sz_nbr=100 --rds_ffc_gmm=10.0 --var_ffc_gmm=0.25 --wvl_nbr=1 --bnd_nbr=1

// Scattering efficiency HaT74 Fig.8, p. 551
   mie -dbg --psd_typ=gamma --sz_grd=log --sz_mnm=.01 --sz_mxm=10 --sz_nbr=1000 --rds_ffc_gmm=10.0 --var_ffc_gmm=0.01 --wvl_nbr=1 --bnd_nbr=1 --idx_rfr_prt="1.33"

// Scattering efficiency HaT74 Fig.9, p. 552
   mie -dbg --psd_typ=gamma --sz_grd=log --sz_mnm=0.01 --sz_mxm=10 --sz_nbr=1000 --rds_ffc_gmm=10 --var_ffc_gmm=0.01 --wvl_nbr=1 --bnd_nbr=1 --idx_rfr_prt="1.33+0.01i"

// TeL96 Size bin 1 (of 8), fully resolved, compare to Table 1, p. 19240
   mie -dbg --psd_typ=gamma --sz_grd=log --sz_mnm=0.01 --sz_mxm=1.0 --sz_nbr=100 --rds_ffc_gmm=0.15 --var_ffc_gmm=0.20 --wvl_nbr=1 --wvl_mnm=0.525 --wvl_mxm=0.575 --bnd_nbr=10
// TeL96 Size bin 5 (of 8), fully resolved, compare to Table 1, p. 19240
   mie -dbg --psd_typ=gamma --sz_grd=log --sz_mnm=0.01 --sz_mxm=10.0 --sz_nbr=100 --rds_ffc_gmm=1.5 --var_ffc_gmm=0.20 --wvl_nbr=1 --wvl_mnm=0.525 --wvl_mxm=0.575 --bnd_nbr=10
   ncks -C -F -v rds_swr,rds_nwr,ext_fsh_ffc ${DATA}/mie/mie.nc

// Desert dust q_ext, omega, g LaM95 Plate 2.1, p. 17 (light blue line)
   mie -dbg --cmp_prt=dust_like --psd_typ=gamma --wvl_mnm=0.3 --wvl_mxm=3.0 --wvl_nbr=100 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-2 --sz_mxm=5.0 --sz_nbr=10 --rds_ffc_gmm=0.5 --var_ffc_gmm=0.2  

// CCM liquid water clouds
mie --dbg=1 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=10 --rds_swa=7.0 --gsd_anl=1.6 --wvl_mnm=0.25 --wvl_mxm=5.0 --wvl_nbr=475 

// Optical properties dust_like aerosol Bri96 ARESE Table 1 (based on PaG77,DKS91)
   mie -dbg --slr_spc=Labs --cmp_prt=dust_like --psd_typ=lognormal --wvl_mnm=0.32799 --wvl_mxm=0.32801 --wvl_nbr=1 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 --idx_rfr_prt="1.53+0.008i"
ncks -C -v ext_cff_mss,ss_alb,asm_prm ${DATA}/mie/mie.nc

// Ice water properties from BPB ${DATA}/aca/ice_20.dat
   mie -dbg --cmp_prt=ice --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=100 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1000.0 --sz_nbr=100 --rds_nma=20.0 --gsd_anl=1.6
   mie -dbg -no_mie --psd_typ=lognormal --sz_grd=log --sz_mnm=1.0 --sz_mxm=100.0 --sz_nbr=100 --rds_nma=20.0 --gsd_anl=1.6 --wvl_nbr=1 --bnd_nbr=1

// Generate deposition velocities: compare dff_aer to SeP97 p. 474 
   mie -dbg --no_mie --psd_typ=lognormal --sz_grd=log --sz_mnm=0.01 --sz_mxm=1.0 --sz_nbr=100 --rds_nma=0.7 --gsd_anl=2.2 

// Generate aerosols with specified single scattering albedos
   mie -dbg --ss_alb=0.50 --slr_spc=Labs --cmp_prt=mineral_dust --psd_typ=lognormal --wvl_mnm=0.3 --wvl_mxm=5.0 --wvl_nbr=94 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/tmp/aer_ssalb_050.nc

   Production usage:
   mie --cmp_prt=dust_like --wvl_mnm=0.3 --wvl_mxm=4.0 --wvl_nbr=100 --bnd_nbr=10
   end debugging */

// Global variables and variables with scope limited to the main.c file allocated here  

// Standard C++ header files 
#include <complex> // Standard C++ complex class
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr
#include <new> // Standard C++ new handler set_new_handler()
#include <string> // Standard C++ string class

// Standard C header files 
#include <cassert> // Assertions
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, EOF, FILE, NULL, etc.
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv
#include <cstring> // strcmp... 
#include <ctime> // Machine time 
#include <unistd.h> // All sorts of POSIX stuff  

// 3rd party vendors
#ifdef _OPENMP
# include <omp.h> // OpenMP pragmas
#endif // not _OPENMP
#include <gsl/gsl_sf_erf.h> // GNU Scientific Library error function
#include <gsl/gsl_sf_gamma.h> // GNU Scientific Library special functions gamma functions
#include <netcdf.h> // netCDF C interface
extern "C" {
#include "getopt.h" // GNU getopt() functionality from BSD my_getopt()
} // end extern

// Personal headers
// #define MAIN_PROGRAM_FILE MUST precede #include mie.hh
#define MAIN_PROGRAM_FILE
#include <mie.hh> // Program-specific definitions
#include <mie_cls.hh> // Program-specific class definitions
#include <libcsz_c++.hh> // Personal C++ library
#include <libcsm_c++.hh> // Climate systems modeling library
#include <libnco_c++.hh> // C++ interface to netCDF C library

// Namespaces

int main(int argc,char **argv)
{
  /* Instead of throwing bad_alloc(), call homebrew handler when new[] fails
     my_new_handler() defined in libcsz_c++:utl.cc
     20040309: AIX optimization has problems with routines with custom handlers */
  // std::set_new_handler(my_new_handler); // [fnc] Handle failures thrown by new[]

  // Enums
  enum slf_tst_typ{BoH83,ViC98}; // [enm] Self-test type
  const std::string slf_tst_sng[]={"Bohren & Huffman (1983) p. 482","Videen & Chylek (1998) p. 5 Tbl. 3"};
  int slf_tst_typ(BoH83); // [enm] Self-test type

  const std::string ffc_mdm_sng[]={"Bruggeman","Multiple Maxwell Garnett","Maxwell Garnett approximation","Partial molar refraction approximation","None, use particle refractive index directly","Volume-weighted approximation"};
  int ffc_mdm_typ(ffc_mdm_nil); // [enm] Effective medium type

  enum cnt_rgm_bll_crc{SeP97,FuS70,ZSK94,ZhC99};
  //  const std::string flx_cnt_rgm_bll_crc_sng[4]={"Seinfeld and Pandis (1997)","Fuchs and Sutugin (1970)","Fuchs and Sutugin (1970)","Fuchs and Sutugin (1970)"};
  int cnt_rgm_bll_crc(FuS70); // default
  cnt_rgm_bll_crc=1.0*cnt_rgm_bll_crc; // default fxm: CEWI Not completely implemented yet

  // Environment
  std::string nvr_DATA((std::getenv("DATA")) ? std::getenv("DATA") : ""); // [sng] Environment variable DATA
  std::string nvr_DATA_RT((std::getenv("DATA_RT")) ? std::getenv("DATA_RT") : ""); // [sng] Environment variable DATA_RT
  std::string nvr_HOME((std::getenv("HOME")) ? std::getenv("HOME") : ""); // [sng] Environment variable HOME
  std::string nvr_HOSTNAME((std::getenv("HOSTNAME")) ? std::getenv("HOSTNAME") : ""); // [sng] Environment variable HOSTNAME
  std::string nvr_USER((std::getenv("USER")) ? std::getenv("USER") : ""); // [sng] Environment variable USER

  // Locals
  bool ffc_mdm_flg(false); // [flg] Employ effective medium approximation

  const prc_cmp ss_alb_bnd_CAM_SW(0.999999); // [frc] Maximum single scattering albedo printed for CAM_SW bands
  const prc_cmp wnd_min_mbl(1.0); // [m s-1] Minimum windspeed used for mobilization
  const long lon_nbr(1); // [nbr] Number of longitudes
  const int prc_cmp_dcm_plc(PRC_CMP_DCM_PLC); // [nbr] Decimal places to right of decimal to represent prc_cmp

  int rcd(0); // [enm] Return code

  // Volume fractions for mantle, matrix, and particle kludged for diagnostics
  prc_cmp vlm_frc_mnt(CEWI_cpv); // [frc] Volume fraction in mantle
  prc_cmp vlm_frc_mtx(CEWI_cpv); // [frc] Volume fraction in matrix
  prc_cmp vlm_frc_prt(CEWI_cpv); // [frc] Volume fraction in particle

  long psd_idx; // [idx] Counting index for aerosol mode
#ifdef PGI_CXX
  // fxm: 20060817 pgCC requires int, rather than long, parallel loop variables
  int bnd_idx; // [idx] Counting index for band
#else // !PGI_CXX
  long bnd_idx; // [idx] Counting index for band
#endif // !PGI_CXX
  long idx; // [idx] Counting index 
  long ncl_idx; // [idx] Counting index for inclusion
  long ngl_idx; // [idx] Counting index for angle
  long sz_idx; // [idx] Counting index for size
  long dsd_idx; // [idx] Counting index for raindrop size
  long wvl_idx; // [idx] Counting index for wavelength
  long wvl_bnd_sz_idx(0); // [idx] Counting index for wavelength/band/size loop

  const std::string CVS_Date("$Date$"); // [sng] CVS date string
  const std::string CVS_Header("$Id$"); // [sng] CVS header string
  const std::string CVS_Id("$Id$"); // [sng] CVS identification string
  const std::string CVS_Name("$HeadURL$"); // [sng] CVS name string
  const std::string CVS_Revision("$Revision$"); // [sng] CVS revision string
  const std::string date_cvs(CVS_Date.length() > 7 ? CVS_Date.substr(7,19) : "Unknown"); // [sng] Date from CVS
  const std::string sbr_nm("main"); // [sng] Subroutine name
  const std::string vrs_cvs(CVS_Revision.length() > 10 ? CVS_Revision.substr(10,4) : "Unknown"); // [sng] Version 

#define TKN2SNG_PRV(x) #x
#define TKN2SNG(x) TKN2SNG_PRV(x)
  const std::string date_cpp(__DATE__); // [sng] Date from C pre-processor
  const std::string time_cpp(__TIME__); // [sng] Time from C pre-processor
  const std::string vrs_cpp(TKN2SNG(VERSION)); // [sng] Version from C pre-processor
  const std::string hst_cpp(TKN2SNG(HOSTNAME)); // [sng] Hostname from C pre-processor
  const std::string usr_cpp(TKN2SNG(USER)); // [sng] Hostname from C pre-processor

  // Start clock and save command line   
  const std::time_t time_crr_time_t(std::time((std::time_t *)NULL)); // [tm] Current date and time
  const std::string time_bfr_srt(std::ctime(&time_crr_time_t)); // [sng] Current date and time
  std::cerr << "\tStart = " << time_bfr_srt;
  prg_nm=((prg_nm=std::strrchr(argv[0],'/')) == NULL) ? argv[0] : ++prg_nm; // [sng] Name of program
  if(vrs_cvs == "Unknown") std::cout << prg_nm << " version " << vrs_cpp << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  if(vrs_cvs != "Unknown") std::cout << prg_nm << " version " << vrs_cvs << " last modified " << date_cvs << " built " << date_cpp << " on " << hst_cpp << " by " << usr_cpp << std::endl;
  const std::string cmd_ln(cmd_ln_sng(argc,argv)); // [sng] Parsed command line
  std::cout << prg_nm << ": INFO Command line = " << cmd_ln << std::endl;

  // Aerosol modes
  const long psd_nbr_max(3); // [nbr] Maximum number of aerosol modes
  long psd_nbr(0); // [nbr] Number of aerosol modes
  bool dsd_usr_flg(false); // [flg] User-specified raindrop mode
  std::string psd_arg[psd_nbr_max]; // [sng] Aerosol mode argument string
  std::string dsd_arg; // [sng] Raindrop mode argument string

  // Aerosol inclusions
  const long ncl_nbr_max(10); // [nbr] Maximum number of inclusions
  long ncl_nbr(0); // [nbr] Number of inclusions
  std::string ncl_arg[ncl_nbr_max]; // [sng] Inclusion properties argument string

  // Set defaults for command line options
  // Control generic program behavior
  // Option string for pure Boolean flags is same as Boolean flag variable name
  bool abc_flg(true); // [flg] Alphabetize output with ncks
  bool abs_ncl_wk_mdm_flg(false); // [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
  bool bch_flg(false); // [flg] Batch behavior (do not print progress to screen)
  bool coat_flg(false); // [flg] Assume coated spheres
  bool coat_nrm_by_core_flg(false); // [flg] Normalize intensive optical properties of coated spheres to core properties
  bool drv_rds_nma_flg(false); // [flg] Derive rds_nma from bin boundaries
  bool fdg_flg(false); // [flg] Tune extinction of particular band
  bool ftn_fxd_flg(false); // [flg] Fortran fixed format
  bool hrz_flg(false); // [flg] Print size-resolved optical properties at debug wavelength
  bool hxg_flg(true); // [flg] Aspherical particles are hexagonal prisms
  bool vts_flg(false); // [flg] Apply equal-V/S approximation for aspherical optical properties
  bool idx_rfr_cor_usr_flg(false); // [flg] Refractive index of core is user-specified
  bool idx_rfr_mdm_usr_flg(false); // [flg] Refractive index of medium is user-specified
  bool idx_rfr_mnt_usr_flg(false); // [flg] Refractive index of mantle is user-specified
  bool idx_rfr_mtx_usr_flg(false); // [flg] Refractive index of matrix is user-specified
  bool idx_rfr_ncl_usr_flg(false); // [flg] Refractive index of inclusion is user-specified
  bool idx_rfr_prt_usr_flg(false); // [flg] Refractive index of particle is user-specified
  bool mie_flg(true); // [flg] Perform mie scattering calculation
  bool phz_flg(false); // [flg] Perform phase function diagnostics
  bool slf_tst_flg(false); // [flg] Perform self-test
  bool ss_alb_flg(false); // [flg] Manually set single scattering albedo
  bool wrn_ntp_flg(true); // [flg] Print WARNINGs from ntp_vec()
  float flt_foo(0.0); // [frc] Intrinsic float temporary variable
  int dmn_nbr_max(2); // [nbr] Maximum number of dimensions allowed in single variable in output file 
  int fl_out_fmt(NCO_FORMAT_UNDEFINED); // [enm] Output file format
  int thr_nbr(0); // [nbr] Thread number
  long bnd_nbr(1); // [nbr] Number of sub-bands per output band
  long dsd_nbr(1); // [nbr] Number of raindrop size bins
  long fdg_idx(0); // [idx] Band to tune by fdg_val
  long lgn_nbr(8); // [nbr] Order of phase function Legendre expansion
  long lng_foo(0); // [nbr] Intrinsic long temporary variable
  long ngl_nbr(11); // [nbr] Number of polar angles in one hemisphere NB: ngl_nbr affects scattering matrix (and hence phase function) but not q_ext, q_sca, q_back
  long pnt_typ_idx(14); // [idx] Plant type index 
  long sfc_typ(2); // [idx] LSM surface type (0..28)
  long soi_typ(1); // [idx] LSM soil type (1..5)
  long sz_nbr(1); // [nbr] Number of particle size bins
  long wvl_nbr(1); // [nbr] Number of output wavelength bands
  long wvn_nbr(1); // [nbr] Number of output wavenumber bands
  prc_cmp RH_lqd(0.8); // [frc] Relative humidity w/r/t liquid water
  prc_cmp asp_rat_hxg_dfl(1.0); // [frc] Hexagonal prism aspect ratio
  prc_cmp asp_rat_lps_dfl(1.0); // [frc] Ellipsoidal aspect ratio
  prc_cmp bnd_SW_LW(5.0e-6); // [m] Boundary between SW and LW weighting
  prc_cmp cnc_nbr_anl_dfl(1.0); // [# m-3] Number concentration analytic, default
  prc_cmp cnc_nbr_pcp_anl(1.0); // [# m-3] Number concentration analytic, raindrop
  prc_cmp cpv_foo(0.0); // [frc] Intrinsic computational precision temporary variable
  prc_cmp dmn_frc(3.0); // I [frc] Fractal dimensionality of inclusions
  prc_cmp dmt_dtc(0.001); // [m] Diameter of detector
  prc_cmp dmt_nma_mcr(cmd_ln_dfl); // [um] Number median diameter analytic, microns
  prc_cmp dmt_pcp_nma_mcr(1000.0); // [um] Diameter number median analytic, raindrop, microns NGD94 p. 2337 Table 1
  prc_cmp dmt_swa_mcr(cmd_ln_dfl); // [um] Surface area weighted mean diameter analytic, microns
  prc_cmp dmt_vma_mcr(cmd_ln_dfl); // [um] Volume median diameter analytic, microns
  prc_cmp dns_cor(0.0); // [kg m-3] Density of core
  prc_cmp dns_mdm(0.0); // [kg m-3] Density of medium
  prc_cmp dns_mnt(0.0); // [kg m-3] Density of mantle
  prc_cmp dns_mtx(0.0); // [kg m-3] Density of matrix
  prc_cmp dns_ncl(0.0); // [kg m-3] Density of inclusion
  prc_cmp dns_prt(0.0); // [kg m-3] Density of particle
  prc_cmp doy(135.0); // [day] Day of year [1.0..367.0)
  prc_cmp dsd_dbg_mcr(1000.0); // [um] Debugging size for raindrops
  prc_cmp dsd_mnm_mcr(999.0); // [um] Minimum diameter in raindrop distribution
  prc_cmp dsd_mxm_mcr(1001.0); // [um] Maximum diameter in raindrop distribution
  prc_cmp fdg_val(1.0); // [frc] Tuning factor for all bands
  prc_cmp flx_LW_dwn_sfc(350.0); // [W m-2] Longwave downwelling flux at surface
  prc_cmp flx_SW_net_gnd(450.0); // [W m-2] Solar flux absorbed by ground
  prc_cmp flx_SW_net_vgt(0.0); // [W m-2] Solar flux absorbed by vegetation
  prc_cmp flx_frc_drc_sfc_cmd_ln(0.85); // [frc] Surface insolation fraction in direct beam
  prc_cmp flx_vlm_pcp_rsl(-1.0); // [m3 m-2 s-1]=[m s-1] Precipitation volume flux, resolved
  prc_cmp gsd_anl_dfl(2.0); // [frc] Geometric standard deviation, default
  prc_cmp gsd_pcp_anl(1.86); // [frc] Geometric standard deviation, raindrop NGD94 p. 2337 Table 1
  prc_cmp hgt_mdp(95.0); // [m] Midlayer height above surface
  prc_cmp hgt_rfr(10.0); // [m] Reference height (i.e., 10 m) at which surface winds are evaluated for dust mobilization
  prc_cmp hgt_zpd_dps_cmd_ln(cmd_ln_dfl); // [m] Zero plane displacement height
  prc_cmp hgt_zpd_mbl(0.0); // [m] Zero plane displacement height for erodible surfaces
  prc_cmp lat_dgr(40.0); // [dgr] Latitude
  prc_cmp lnd_frc_dry(1.0); // [frc] Dry land fraction
  prc_cmp mie_cnv_eps(1.0e-6); // [frc] Mie convergence precision
  prc_cmp mmw_prt(0.0); // [kg mol-1] Mean molecular weight of particle
  prc_cmp mno_lng_dps_cmd_ln(cmd_ln_dfl); // [m] Monin-Obukhov length
  prc_cmp mss_frc_cly(0.19); // [frc] Mass fraction clay
  prc_cmp mss_frc_snd(0.777); // [frc] Mass fraction sand
  prc_cmp mss_frc_cor(cmd_ln_dfl); // [frc] Mass fraction in core
  prc_cmp mss_frc_ncl(cmd_ln_dfl); // [frc] Mass fraction in inclusion
  prc_cmp msv_gnd(cmd_ln_dfl); // [frc] Bare ground emissivity
  prc_cmp msv_snw(cmd_ln_dfl); // [frc] Emissivity of snow
  prc_cmp ngl_dbg_dgr(0.0); // [dgr] Debugging angle
  prc_cmp oro(1.0); // [frc] Orography: ocean=0.0, land=1.0, sea ice=2.0
  prc_cmp prs_mdp(100825.0); // [Pa] Environmental pressure
  prc_cmp prs_ntf(prs_STP); // [Pa] Environmental surface pressure
  prc_cmp q_H2O_vpr(cmd_ln_dfl); // [kg kg-1] Specific humidity
  prc_cmp rds_ffc_gmm_mcr(50.0); // [um] Effective radius of Gamma distribution
  prc_cmp rds_frc_cor(cmd_ln_dfl); // [frc] Radius fraction in core
  prc_cmp rds_frc_ncl(cmd_ln_dfl); // [frc] Radius fraction in inclusion
  prc_cmp rds_nma_mcr(0.2986); // [um] Number median radius analytic, microns
  prc_cmp sfc_spc_cm2xg(cmd_ln_dfl); // [cm2 g-1] Specific surface area, CGS
  prc_cmp rds_swa_mcr(cmd_ln_dfl); // [um] Surface area weighted mean radius analytic, microns
  prc_cmp rds_vma_mcr(cmd_ln_dfl); // [um] Volume median radius analytic, microns
  prc_cmp rfl_gnd_dff(0.2); // [frc] Diffuse reflectance of ground (beneath snow)
  prc_cmp rgh_mmn_dps_cmd_ln(cmd_ln_dfl); // [m] Roughness length momentum
  prc_cmp rgh_mmn_ice_std(0.0005); // [m] Roughness length over sea ice BKL97 p. F-3 (updated)
  prc_cmp rgh_mmn_mbl(100.0e-6); // [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
  prc_cmp rgh_mmn_smt(10.0e-6); // [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207 fxm: aer uses 30.0e-6
  prc_cmp slr_cst(1367.0); // [W m-2] Solar constant
  prc_cmp slr_zen_ngl_cos(1.0); // [frc] Cosine solar zenith angle
  prc_cmp snw_hgt_lqd(0.0); // [m] Equivalent liquid water snow depth
  prc_cmp spc_heat_prt(0.0); // [J kg-1 K-1] Specific heat capacity of particle
  prc_cmp ss_alb_cmd_ln(1.0); // [frc] Single scattering albedo
  prc_cmp sz_dbg_mcr(1.0); // [um] Debugging size
  prc_cmp sz_mnm_mcr(0.9); // [um] Minimum size in distribution
  prc_cmp sz_mxm_mcr(1.1); // [um] Maximum size in distribution
  prc_cmp sz_prm_rsn_usr_spc(0.1); // [m m-1] Size parameter resolution, user specified
  prc_cmp tm_dlt(1200.0); // [s] Timestep
  prc_cmp tpt_bbd_wgt(273.15); // [K] Blackbody temperature of radiation
  prc_cmp tpt_gnd(300.0); // [K] Ground temperature
  prc_cmp tpt_ice(tpt_frz_pnt); // [K] Ice temperature
  prc_cmp tpt_mdp(300.0); // [K] Environmental temperature
  prc_cmp tpt_prt(273.15); // [K] Particle temperature
  prc_cmp tpt_soi(297.0); // [K] Soil temperature
  prc_cmp tpt_sst(300.0); // [K] Sea surface temperature
  prc_cmp tpt_vgt(300.0); // [K] Vegetation temperature
  prc_cmp var_ffc_gmm(1.0); // [frc] Effective variance of Gamma distribution
  prc_cmp vlm_frc_cor(0.5); // [frc] Volume fraction in core
  prc_cmp vlm_frc_ncl(0.5); // [frc] Volume fraction in inclusion
  prc_cmp vmr_CO2(355e-6); // [mlc mlc-1] Volume mixing ratio of CO2
  prc_cmp vmr_HNO3_gas(0.05e-9); // [mlc mlc-1] Volume mixing ratio of gaseous HNO3 ZhC99 p. 358 Tbl. 3
  prc_cmp vwc_sfc(0.03); // [m3 m-3] Volumetric water content
  prc_cmp wbl_shp(2.4); // [frc] Weibull distribution shape parameter
  prc_cmp wnd_frc_dps_cmd_ln(cmd_ln_dfl); // [m s-1] Friction speed
  prc_cmp wnd_mrd_mdp(0.0); // [m s-1] Surface layer meridional wind speed
  prc_cmp wnd_znl_mdp(10.0); // [m s-1] Surface layer zonal wind speed
  prc_cmp wvl_dbg_mcr(0.50); // [um] Debugging wavelength
  prc_cmp wvl_dlt_mcr(cmd_ln_dfl); // [um] Bandwidth
  prc_cmp wvl_mdp_mcr(cmd_ln_dfl); // [um] Midpoint wavelength
  prc_cmp wvl_mnm_mcr(0.45); // [um] Minimum wavelength
  prc_cmp wvl_mxm_mcr(0.55); // [um] Maximum wavelength
  prc_cmp wvn_dlt_xcm(cmd_ln_dfl); // [cm-1] Bandwidth
  prc_cmp wvn_mdp_xcm(cmd_ln_dfl); // [cm-1] Midpoint wavenumber
  prc_cmp wvn_mnm_xcm(cmd_ln_dfl); // [cm-1] Minimum wavenumber
  prc_cmp wvn_mxm_xcm(cmd_ln_dfl); // [cm-1] Maximum wavenumber
  std::complex<prc_cmp> idx_rfr_cor_usr(1.0,0.0); // [frc] Refractive index of core
  std::complex<prc_cmp> idx_rfr_mdm_usr(1.0,0.0); // [frc] Refractive index of medium
  std::complex<prc_cmp> idx_rfr_mnt_usr(1.33,0.0); // [frc] Refractive index of mantle
  std::complex<prc_cmp> idx_rfr_mtx_usr(1.0,0.0); // [frc] Refractive index of matrix
  std::complex<prc_cmp> idx_rfr_ncl_usr(1.0,0.0); // [frc] Refractive index of inclusion
  std::complex<prc_cmp> idx_rfr_prt_usr(1.33,0.0); // [frc] Refractive index of particle
  std::string dmn_rcd("wvl"); // [sng] Record dimension name
  std::string fl_err("foo.stderr"); // [sng] File for error messages
  std::string fl_idx_rfr_cor(""); // [sng] File or function for refractive indices of core
  std::string fl_idx_rfr_mdm(""); // [sng] File or function for refractive indices of medium
  std::string fl_idx_rfr_mnt(""); // [sng] File or function for refractive indices of mantle
  std::string fl_idx_rfr_mtx(""); // [sng] File or function for refractive indices of matrix
  std::string fl_idx_rfr_ncl(""); // [sng] File or function for refractive indices of inclusion
  std::string fl_idx_rfr_prt(""); // [sng] File or function for refractive indices of particle
  std::string fl_out("mie.nc"); // [sng] Output file
  std::string fl_slr_spc(""); // [sng] File or function for solar spectrum
  std::string lbl_tst("CO2"); // [sng] Name of line-by-line test
  std::string cmp_sng_cor("saharan_dust"); // [sng] Composition of core
  std::string cmp_sng_mdm("air"); // [sng] Composition of medium
  std::string cmp_sng_mnt("h2o_lqd"); // [sng] Composition of mantle
  std::string cmp_sng_mtx("SiO2"); // [sng] Composition of matrix
  std::string cmp_sng_ncl("Fe2O3"); // [sng] Composition of inclusion
  std::string cmp_sng_prt("saharan_dust"); // [sng] Composition of particle
  std::string psd_typ("lognormal"); // [sng] Particle size distribution type
  std::string slr_spc_key("LaN68"); // [sng] Solar spectrum string
  std::string slv_sng("Wis79"); // [sng] Mie solver to use (BoH83 or Wis79)
  std::string spc_idx_sng("foo"); // [sng] Label for FORTRAN block data
  std::string spc_wgt_sng("Default"); // [sng] Spectral weight
  std::string sz_grd_sng("logarithmic"); // [sng] Type of size grid
  std::string ngl_sng("lobatto"); // [sng] Angle grid type
  std::string tst_sng(""); // [sng] Name of test to perform
  std::string wvl_grd_sng("regular"); // [sng] Type of wavelength grid
  std::string xpt_dsc(""); // [sng] Experiment description

  // Derived fields
  std::string drc_dat((nvr_DATA_RT.length() > 0) ? nvr_DATA_RT : "/data/zender/aca"); // [sng] Data directory
  std::string drc_in((nvr_HOME.length() > 0) ? nvr_HOME+"/nco" : "/home/zender/nco/data"); // [sng] Input directory
  std::string drc_out((nvr_DATA.length() > 0) ? nvr_DATA+"/mie" : ""); // [sng] Output directory
  std::string spc_abb_sng(cmp_sng_prt); // [sng] Species abbreviation for Fortran data

  static struct option opt_lng[]={
    /* The option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is enum _argtype{no_argument,required_argument,optional_argument}
       If flag is non-zero, getopt_long() returns zero and flag is set to val
       If flag is zero, getopt_long() returns contents of val */
    // Long options with no argument, no short option counterpart
    {"abc_flg",no_argument,0,0}, // [flg] Alphabetize output with ncks
    {"abs_ncl_wk_mdm_flg",no_argument,0,0}, // [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
    {"bch_flg",no_argument,0,0}, // [flg] Batch behavior (do not print progress to screen)
    {"coat_flg",no_argument,0,0}, // [flg] Assume coated spheres
    {"coat_nrm_by_core_flg",no_argument,0,0}, // [flg] Normalize intensive optical properties of coated spheres to core properties
    {"ftn_fxd_flg",no_argument,0,0}, // [flg] Fortran fixed format
    {"drv_rds_nma_flg",no_argument,0,0}, // [flg] Derive rds_nma from bin boundaries
    {"hrz_flg",no_argument,0,0}, // [flg] Print size-resolved optical properties at debug wavelength
    {"hxg_flg",no_argument,0,0}, // [flg] Aspherical particles are hexagonal prisms
    {"vts_flg",no_argument,0,0}, // [flg] Apply equal-V/S approximation for aspherical optical properties
    {"mie_flg",no_argument,0,0}, // [flg] Perform mie scattering calculation
    {"no_abc_flg",no_argument,0,0}, // [flg] Do not alphabetize output with ncks
    {"no_bch_flg",no_argument,0,0}, // [flg] Interactive/console behavior
    {"no_hrz_flg",no_argument,0,0}, // [flg] Do not print size-resolved optical properties at debug wavelength
    {"no_mie_flg",no_argument,0,0}, // [flg] Do not perform mie scattering calculation
    {"no_wrn_ntp_flg",no_argument,0,0}, // [flg] Do not print WARNINGs from ntp_vec()
    // Long options with optional argument, no short option counterpart
    // Long options with argument, no short option counterpart
    {"RH_lqd",required_argument,0,0}, // [frc] Relative humidity w/r/t liquid water
    {"asp_rat_hxg_dfl",required_argument,0,0}, // [frc] Hexagonal prism aspect ratio
    {"asp_rat_lps_dfl",required_argument,0,0}, // [frc] Ellipsoidal aspect ratio
    {"bnd_SW_LW",required_argument,0,0}, // [m] Boundary between SW and LW weighting
    {"bnd_nbr",required_argument,0,0}, // [nbr] Number of sub-bands per output band
    {"cnc_nbr_anl_dfl",required_argument,0,0}, // [# m-3] Number concentration analytic, default
    {"cpv_foo",required_argument,0,0}, // [frc] Intrinsic computational precision temporary variable
    {"distribution",required_argument,0,0}, // [sng] Type of size distribution
    {"dmn_frc",required_argument,0,0}, // I [frc] Fractal dimensionality of inclusions
    {"dmn_nbr_max",required_argument,0,0}, // [nbr] Maximum number of dimensions allowed in single variable in output file
    {"dmn_rcd",required_argument,0,0}, // [sng] Record dimension name
    {"dmt_dtc",required_argument,0,0}, // [m] Diameter of detector
    {"dmt_nma_mcr",required_argument,0,0}, // [um] Number median diameter analytic, microns
    {"dmt_pcp_nma_mcr",required_argument,0,0}, // [um] Diameter number median analytic, raindrop, microns
    {"dmt_swa_mcr",required_argument,0,0}, // [um] Surface area weighted mean diameter analytic, microns
    {"dmt_vma_mcr",required_argument,0,0}, // [um] Volume median diameter analytic, microns
    {"dns_prt",required_argument,0,0}, // [kg m-3] Density of particle
    {"dns_cor",required_argument,0,0}, // [kg m-3] Density of core
    {"dns_mdm",required_argument,0,0}, // [kg m-3] Density of medium
    {"dns_mnt",required_argument,0,0}, // [kg m-3] Density of mantle
    {"dns_mtx",required_argument,0,0}, // [kg m-3] Density of matrix
    {"dns_ncl",required_argument,0,0}, // [kg m-3] Density of inclusion
    {"drc_dat",required_argument,0,0}, // [sng] Data directory
    {"drc_in",required_argument,0,0}, // [sng] Input directory
    {"drc_out",required_argument,0,0}, // [sng] Output directory
    {"dsd",required_argument,0,0}, // [sng] Raindrop mode argument string
    {"dsd_dbg_mcr",required_argument,0,0}, // [um] Debugging size for raindrops
    {"dsd_mnm_mcr",required_argument,0,0}, // [um] Minimum diameter in raindrop distribution
    {"dsd_mxm_mcr",required_argument,0,0}, // [um] Maximum diameter in raindrop distribution
    {"dsd_nbr",required_argument,0,0}, // [nbr] Number of raindrop size bins
    {"fdg_idx",required_argument,0,0}, // [idx] Band to tune by fdg_val
    {"fdg_val",required_argument,0,0}, // [frc] Tuning factor for all bands
    {"ffc_mdm_typ",required_argument,0,0}, // [enm] Effective medium type
    {"fl_idx_rfr_cor",required_argument,0,0}, // [sng] File or function for refractive indices of core
    {"fl_idx_rfr_mdm",required_argument,0,0}, // [sng] File or function for refractive indices of medium
    {"fl_idx_rfr_mnt",required_argument,0,0}, // [sng] File or function for refractive indices of mantle
    {"fl_idx_rfr_mtx",required_argument,0,0}, // [sng] File or function for refractive indices of matrix
    {"fl_idx_rfr_ncl",required_argument,0,0}, // [sng] File or function for refractive indices of inclusion
    {"fl_fmt",required_argument,0,0}, // [enm] Output file format
    {"file_format",required_argument,0,0}, // [enm] Output file format
    {"fl_slr_spc",required_argument,0,0}, // [sng] File containing solar spectrum
    {"flt_foo",required_argument,0,0}, // [frc] Intrinsic float temporary variable
    {"flx_LW_dwn_sfc",required_argument,0,0}, // [W m-2] Longwave downwelling flux at surface
    {"flx_SW_net_gnd",required_argument,0,0}, // [W m-2] Solar flux absorbed by ground
    {"flx_SW_net_vgt",required_argument,0,0}, // [W m-2] Solar flux absorbed by vegetation
    {"flx_frc_drc_sfc_cmd_ln",required_argument,0,0}, // [frc] Surface insolation fraction in direct beam
    {"flx_vlm_pcp_rsl",required_argument,0,0}, // [m3 m-2 s-1]=[m s-1] Precipitation volume flux, resolved
    {"gsd_anl_dfl",required_argument,0,0}, // [frc] Geometric standard deviation, default
    {"gsd_pcp_anl",required_argument,0,0}, // [frc] Geometric standard deviation, raindrop
    {"hgt_mdp",required_argument,0,0}, // [m] Midlayer height above surface
    {"hgt_rfr",required_argument,0,0}, // [m] Reference height (i.e., 10 m) at which surface winds are evaluated for dust mobilization
    {"hgt_zpd_dps_cmd_ln",required_argument,0,0}, // [m] Zero plane displacement height
    {"hgt_zpd_mbl",required_argument,0,0}, // [m] Zero plane displacement height for erodible surfaces
    {"idx_rfr_cor_usr",required_argument,0,0}, // [frc] Refractive index of core
    {"idx_rfr_mdm_usr",required_argument,0,0}, // [frc] Refractive index of medium
    {"idx_rfr_mtx_usr",required_argument,0,0}, // [frc] Refractive index of matrix
    {"idx_rfr_mnt_usr",required_argument,0,0}, // [flg] Refractive index of mantle
    {"idx_rfr_ncl_usr",required_argument,0,0}, // [frc] Refractive index of inclusion
    {"idx_rfr_prt_usr",required_argument,0,0}, // [frc] Refractive index of particle
    {"lbl_tst",required_argument,0,0}, // [sng] Name of line-by-line test
    {"lng_foo",required_argument,0,0}, // [nbr] Intrinsic long temporary variable 
    {"lnd_frc_dry",required_argument,0,0}, // [frc] Dry land fraction
    {"cmp_prt",required_argument,0,0}, // [sng] Composition of particle
    {"cmp_cor",required_argument,0,0}, // [sng] Composition of core
    {"cmp_mdm",required_argument,0,0}, // [sng] Composition of medium
    {"cmp_mnt",required_argument,0,0}, // [sng] Composition of mantle
    {"cmp_mtx",required_argument,0,0}, // [sng] Composition of matrix
    {"cmp_ncl",required_argument,0,0}, // [sng] Composition of inclusion
    {"mie_cnv_eps",required_argument,0,0}, // [frc] Mie convergence precision
    {"mie_slv",required_argument,0,0}, // [sng] Mie solver to use (BoH83 or Wis79)
    {"mno_lng_dps_cmd_ln",required_argument,0,0}, // [m] Monin-Obukhov length
    {"mss_frc_cly",required_argument,0,0}, // [frc] Mass fraction clay 
    {"mss_frc_snd",required_argument,0,0}, // [frc] Mass fraction sand
    {"mss_frc_cor",required_argument,0,0}, // [frc] Mass fraction in core
    {"mss_frc_ncl",required_argument,0,0}, // [frc] Mass fraction in inclusion
    {"msv_gnd",required_argument,0,0}, // [frc] Bare ground emissivity
    {"msv_snw",required_argument,0,0}, // [frc] Snow emissivity
    {"lgn_nbr",required_argument,0,0}, // [nbr] Order of phase function Legendre expansion
    {"ngl_dbg_dgr",required_argument,0,0}, // [dgr] Debugging angle
    {"ngl_nbr",required_argument,0,0}, // [nbr] Number of polar angles in one hemisphere
    {"ngl_sng",required_argument,0,0}, // [sng] Angle grid type
    {"oro",required_argument,0,0}, // [frc] Orography: ocean=0.0, land=1.0, sea ice=2.0
    {"pnt_typ_idx",required_argument,0,0}, // [idx] Plant type index 
    {"prs_mdp",required_argument,0,0}, // [Pa] Environmental pressure
    {"prs_ntf",required_argument,0,0}, // [Pa] Environmental surface pressure
    {"ncl",required_argument,0,0}, // [sng] Inclusion properties
    {"psd",required_argument,0,0}, // [sng] Aerosol mode argument string
    {"psd_typ",required_argument,0,0}, // [sng] Particle size distribution type
    {"q_H2O_vpr",required_argument,0,0}, // [kg kg-1] Specific humidity
    {"rds_ffc_gmm_mcr",required_argument,0,0}, // [um] Effective radius of Gamma distribution
    {"rds_nma_mcr",required_argument,0,0}, // [um] Number median radius analytic, microns
    {"sfc_spc_cm2xg",required_argument,0,0}, // [cm2 g-1] Specific surface area, CGS
    {"rds_swa_mcr",required_argument,0,0}, // [um] Surface area weighted mean radius analytic, microns
    {"rds_vma_mcr",required_argument,0,0}, // [um] Volume median radius analytic, microns
    {"rgh_mmn_dps_cmd_ln",required_argument,0,0}, // [m] Roughness length momentum
    {"rgh_mmn_ice_std",required_argument,0,0}, // [m] Roughness length over sea ice
    {"rgh_mmn_mbl",required_argument,0,0}, // [m] Roughness length momentum for erodible surfaces
    {"rgh_mmn_smt",required_argument,0,0}, // [m] Smooth roughness length
    {"rfl_gnd_dff",required_argument,0,0}, // [frc] Diffuse reflectance of ground (beneath snow)
    {"sfc_typ",required_argument,0,0}, // [idx] LSM surface type (0..28)
    {"slr_spc_key",required_argument,0,0}, // [sng] Solar spectrum string
    {"slr_zen_ngl_cos",required_argument,0,0}, // [frc] Cosine solar zenith angle
    {"snw_hgt_lqd",required_argument,0,0}, // [m] Equivalent liquid water snow depth
    {"soi_typ",required_argument,0,0}, // [idx] LSM soil type (1..5)
    {"spc_abb_sng",required_argument,0,0}, // [sng] Species abbreviation for Fortran data
    {"spc_idx_sng",required_argument,0,0}, // [sng] Species index for Fortran data
    {"ss_alb",required_argument,0,0}, // [frc] Single scattering albedo
    {"sz_dbg_mcr",required_argument,0,0}, // [um] Debugging size
    {"sz_grd_typ",required_argument,0,0}, // [sng] Type of size grid
    {"sz_mnm_mcr",required_argument,0,0}, // [um] Minimum size in distribution
    {"sz_mxm_mcr",required_argument,0,0}, // [um] Maximum size in distribution
    {"sz_prm_rsn_usr_spc",required_argument,0,0}, // [m m-1] Size parameter resolution, user specified
    {"sz_nbr",required_argument,0,0}, // [nbr] Number of particle size bins
    {"thr_nbr",required_argument,0,0}, // [nbr] Thread number
    {"tm_dlt",required_argument,0,0}, // [s] Timestep
    {"tpt_bbd_wgt",required_argument,0,0}, // [K] Blackbody temperature of radiation
    {"tpt_gnd",required_argument,0,0}, // [K] Ground temperature
    {"tpt_ice",required_argument,0,0}, // [K] Ice temperature
    {"tpt_mdp",required_argument,0,0}, // [K] Environmental temperature
    {"tpt_soi",required_argument,0,0}, // [K] Soil temperature
    {"tpt_sst",required_argument,0,0}, // [K] Sea surface temperature
    {"tpt_vgt",required_argument,0,0}, // [K] Vegetation temperature
    {"slf_tst_typ",required_argument,0,0}, // [enm] Self-test type
    {"tst_sng",required_argument,0,0}, // [sng] Name of test to perform
    {"var_ffc_gmm",required_argument,0,0}, // [frc] Effective variance of Gamma distribution
    {"rds_frc_cor",required_argument,0,0}, // [frc] Radius fraction in core
    {"rds_frc_ncl",required_argument,0,0}, // [frc] Radius fraction in inclusion
    {"vlm_frc_cor",required_argument,0,0}, // [frc] Volume fraction in core
    {"vlm_frc_ncl",required_argument,0,0}, // [frc] Volume fraction in inclusion
    {"vmr_CO2",required_argument,0,0}, // [mlc mlc-1] Volume mixing ratio of CO2
    {"vmr_HNO3_gas",required_argument,0,0}, // [mlc mlc-1] Volume mixing ratio of gaseous HNO3
    {"vwc_sfc",required_argument,0,0}, // [m3 m-3] Volumetric water content
    {"wbl_shp",required_argument,0,0}, // [frc] Weibull distribution shape parameter
    {"wnd_frc_dps_cmd_ln",required_argument,0,0}, // [m s-1] Friction speed
    {"wnd_mrd_mdp",required_argument,0,0}, // [m s-1] Surface layer meridional wind speed
    {"wnd_znl_mdp",required_argument,0,0}, // [m s-1] Surface layer zonal wind speed
    {"wvl_dbg_mcr",required_argument,0,0}, // [um] Debugging wavelength
    {"wvl_dlt_mcr",required_argument,0,0}, // [um] Bandwidth
    {"wvl_grd_typ",required_argument,0,0}, // [enm] Wavelength grid type
    {"wvl_mdp_mcr",required_argument,0,0}, // [um] Midpoint wavelength
    {"wvl_mnm_mcr",required_argument,0,0}, // [um] Minimum wavelength
    {"wvl_mxm_mcr",required_argument,0,0}, // [um] Maximum wavelength
    {"wvl_nbr",required_argument,0,0}, // [nbr] Number of output wavelength bands
    {"wvn_dlt_xcm",required_argument,0,0}, // [cm-1] Bandwidth
    {"wvn_mdp_xcm",required_argument,0,0}, // [cm-1] Midpoint wavenumber
    {"wvn_mnm_xcm",required_argument,0,0}, // [cm-1] Minimum wavenumber
    {"wvn_mxm_xcm",required_argument,0,0}, // [cm-1] Maximum wavenumber
    {"wvn_nbr",required_argument,0,0}, // [nbr] Number of output wavenumber bands
    {"xpt_dsc",required_argument,0,0}, // [sng] Experiment description
    // Long options with short counterparts
    {"4",no_argument,0,'4'}, // [enm] Output file format
    {"64bit",no_argument,0,'4'}, // [enm] Output file format
    {"netcdf4",no_argument,0,'4'}, // [enm] Output file format
    {"dbg_lvl",optional_argument,0,'D'}, // [enm] Debugging level
    {"fl_idx_rfr_prt",required_argument,0,'i'}, // [sng] File containing refractive indices of particle
    {"help",no_argument,0,'h'}, // [flg] Print proper usage 
    {"fl_in",required_argument,0,'i'}, // [sng] File containing refractive indices of particle
    {"fl_out",required_argument,0,'o'}, // [sng] Output file
    {"slr_cst",required_argument,0,'s'}, // [W m-2] Solar constant
    {"tpt_prt",required_argument,0,'t'}, // [K] Particle temperature
    {"version",no_argument,0,'v'}, // [flg] Print CVS information
    // Last option named "0" signals getopt_long() to stop processing  
    {0,0,0,0}
  }; // end opt_lng
  
  // Short options: no colon = no arg, one colon = required arg, two colons = optional arg
  const char * const opt_sht_lst="34D::e:hi:o:s:t:v"; // [sng] List of single-letter (C-style) option abbreviations
  extern char *optarg; // [sng] char * representation of current optarg, if any (this memory is owned by system)
  extern int optind; // [idx] extern enumerating cardinal of current option
  int opt; // [idx] Value is zero if current argument is long type, else value contains single letter version of command line argument
  int opt_idx=0; // [idx] Index of current long option into opt_lng array
  std::string opt_crr; // [sng] String representation of current long-option name
  std::string opt_sng; // [sng] String representation of current optarg, if any
  const std::string opt_ndc_sng("--"); // [sng] Command line option indicator string
 
  // Parse command line arguments 
  while(1){
    // getopt_long_only() allows a single dash '-' to prefix long options as well
    opt=getopt_long_only(argc,argv,opt_sht_lst,opt_lng,&opt_idx);
    // NB: access to opt_crr is only valid when long_opt was detected
    opt_crr=opt_lng[opt_idx].name; // [sng] String representation of current long-option name
    if(optarg) opt_sng=optarg; // Copy system memory into C++ string for safer operations 
    // Attempt to spot mal-formed command lines
    if(opt_sng.find(opt_ndc_sng) != std::string::npos) err_prn(sbr_nm,"Option indicator \""+opt_ndc_sng+"\" present inside option value \""+opt_sng+"\". Command line syntax error? HINT: Make sure no options are inadvertently conjoined with other options.");
    if(opt == EOF) break; // Parse positional arguments once getopt_long_only() returns EOF
    // Process long options without short option counterparts
    if(opt == 0){
      // fxm: TODO #12 Optimize option parsing to avoid unecessary if statements, perhaps can generate a hash key?
      // fxm: TODO #13 Need mechanism to warn/exit when user option is found in option database structure opt_lng but not acted on later
      if(dbg_lvl >= dbg_io) std::cerr << "Long option name: " << opt_crr << (optarg ? ",  Argument: "+opt_sng : ", No Argument") << std::endl;
      if(opt_crr == "RH_lqd") RH_lqd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc]
      if(opt_crr == "dsd"){
	dsd_usr_flg=true; // [flg] User specified raindrop mode
	dsd_arg=opt_sng; // [sng] Raindrop mode argument string
      } // endif "dsd"
      if(opt_crr == "ncl"){
	ncl_arg[ncl_nbr]=opt_sng; // [frc] Inclusion properties argument string
	ncl_nbr++; // [nbr] Number of inclusions
	assert(ncl_nbr <= ncl_nbr_max);
      } // endif "ncl"
      if(opt_crr == "psd"){
	psd_arg[psd_nbr]=opt_sng; // [sng] Aerosol mode argument string
	psd_nbr++; // [nbr] Number of aerosol modes
	assert(psd_nbr <= psd_nbr_max);
      } // endif "psd"
      if(opt_crr == "abc_flg") abc_flg=true; // [flg] Alphabetize output with ncks
      if(opt_crr == "abs_ncl_wk_mdm_flg") abs_ncl_wk_mdm_flg=true; // [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
      if(opt_crr == "cmp_prt"){cmp_sng_prt=aer_cls::opt2abb(opt_sng); spc_abb_sng=cmp_sng_prt;} // [sng] Composition of particle
      if(opt_crr == "cmp_cor") cmp_sng_cor=aer_cls::opt2abb(opt_sng); // [sng] Composition of core
      if(opt_crr == "cmp_mdm") cmp_sng_mdm=aer_cls::opt2abb(opt_sng); // [sng] Composition of medium
      if(opt_crr == "cmp_mnt") cmp_sng_mnt=aer_cls::opt2abb(opt_sng); // [sng] Composition of mantle
      if(opt_crr == "cmp_mtx") cmp_sng_mtx=aer_cls::opt2abb(opt_sng); // [sng] Composition of matrix
      if(opt_crr == "cmp_ncl") cmp_sng_ncl=aer_cls::opt2abb(opt_sng); // [sng] Composition of inclusion
      if(opt_crr == "asp_rat_hxg_dfl") asp_rat_hxg_dfl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Hexagonal prism aspect ratio
      if(opt_crr == "asp_rat_lps_dfl") asp_rat_lps_dfl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Ellipsoidal aspect ratio
      if(opt_crr == "bch_flg") bch_flg=true; // [flg] Batch behavior (do not print progress to screen)
      if(opt_crr == "bnd_SW_LW") bnd_SW_LW=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Boundary between SW and LW weighting
      if(opt_crr == "bnd_nbr") bnd_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of sub-bands per output band
      if(opt_crr == "cnc_nbr_anl_dfl") cnc_nbr_anl_dfl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL));
      if(opt_crr == "coat_flg") coat_flg=true; // [flg] Assume coated spheres
      if(opt_crr == "coat_nrm_by_core_flg") coat_nrm_by_core_flg=true; // [flg] Normalize intensive optical properties of coated spheres to core properties
      if(opt_crr == "cpv_foo") cpv_foo=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Intrinsic computational precision temporary variable
      if(opt_crr == "dmn_frc") dmn_frc=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // I [frc] Fractal dimensionality of inclusions
      if(opt_crr == "dmn_nbr_max") dmn_nbr_max=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Maximum number of dimensions allowed in single variable in output file
      if(opt_crr == "dmn_rcd") dmn_rcd=opt_sng; // [sng] Record dimension name
      if(opt_crr == "dmt_dtc") dmt_dtc=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Diameter of detector
      if(opt_crr == "dmt_nma_mcr") dmt_nma_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Number median diameter analytic, microns
      if(opt_crr == "dmt_pcp_nma_mcr") dmt_pcp_nma_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Diameter number median analytic, raindrop, microns
      if(opt_crr == "dmt_swa_mcr") dmt_swa_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Surface area weighted mean diameter analytic, microns
      if(opt_crr == "dmt_vma_mcr") dmt_vma_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Volume median diameter analytic, microns
      if(opt_crr == "dns_prt") dns_prt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of particle
      if(opt_crr == "dns_cor") dns_cor=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of core
      if(opt_crr == "dns_mdm") dns_mdm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of medium
      if(opt_crr == "dns_mnt") dns_mnt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of mantle
      if(opt_crr == "dns_mtx") dns_mtx=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of matrix
      if(opt_crr == "dns_ncl") dns_ncl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg m-3] Density of inclusion
      if(opt_crr == "doy") doy=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [day] Day of year
      if(opt_crr == "drc_dat") drc_dat=opt_sng; // [sng] Data directory
      if(opt_crr == "drc_in") drc_in=opt_sng; // [sng] Input directory
      if(opt_crr == "drc_out") drc_out=opt_sng; // [sng] Output directory
      if(opt_crr == "drv_rds_nma_flg") drv_rds_nma_flg=true; // [flg] Derive rds_nma from bin boundaries
      if(opt_crr == "dsd_dbg_mcr") dsd_dbg_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Debugging size for raindrops
      if(opt_crr == "dsd_mnm_mcr") dsd_mnm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Minimum diameter in raindrop distribution
      if(opt_crr == "dsd_mxm_mcr") dsd_mxm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Maximum diameter in raindrop distribution
      if(opt_crr == "dsd_nbr") dsd_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of raindrop size bins
      if(opt_crr == "spc_abb_sng") spc_abb_sng=opt_sng; // [sng] Species abbreviation for Fortran data
      if(opt_crr == "spc_idx_sng") spc_idx_sng=opt_sng; // [sng] Species index for Fortran data
      if(opt_crr == "fdg_idx") fdg_idx=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [idx] Band to tune by fdg_val
      if(opt_crr == "fl_idx_rfr_cor") fl_idx_rfr_cor=opt_sng; // [sng] File or function for refractive indices of core
      if(opt_crr == "fl_idx_rfr_mdm") fl_idx_rfr_mdm=opt_sng; // [sng] File or function for refractive indices of medium
      if(opt_crr == "fl_idx_rfr_mnt") fl_idx_rfr_mnt=opt_sng; // [sng] File or function for refractive indices of mantle
      if(opt_crr == "fl_idx_rfr_mtx") fl_idx_rfr_mtx=opt_sng; // [sng] File or function for refractive indices of matrix
      if(opt_crr == "fl_idx_rfr_ncl") fl_idx_rfr_ncl=opt_sng; // [sng] File or function for refractive indices of inclusion
      if(opt_crr == "fl_fmt" || opt_crr == "file_format") rcd=nco_create_mode_prs(opt_sng,fl_out_fmt); // [enm] Output file format
      if(opt_crr == "fl_slr_spc") fl_slr_spc=opt_sng; // [sng] File or function for solar spectrum
      if(opt_crr == "flt_foo") flt_foo=static_cast<float>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Intrinsic float temporary variable
      if(opt_crr == "flx_LW_dwn_sfc") flx_LW_dwn_sfc=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [W m-2] Longwave downwelling flux at surface
      if(opt_crr == "flx_SW_net_gnd") flx_SW_net_gnd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [W m-2] Solar flux absorbed by ground
      if(opt_crr == "flx_SW_net_vgt") flx_SW_net_vgt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [W m-2] Solar flux absorbed by vegetation
      if(opt_crr == "flx_frc_drc_sfc_cmd_ln") flx_frc_drc_sfc_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Surface insolation fraction in direct beam
      if(opt_crr == "flx_vlm_pcp_rsl") flx_vlm_pcp_rsl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m3 m-2 s-1]=[m s-1] Precipitation volume flux, resolved
      if(opt_crr == "ftn_fxd_flg") ftn_fxd_flg=true; // [flg] Fortran fixed format
      if(opt_crr == "gsd_anl_dfl") gsd_anl_dfl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Geometric standard deviation, default
      if(opt_crr == "gsd_pcp_anl") gsd_pcp_anl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Geometric standard deviation, raindrop
      if(opt_crr == "hgt_mdp") hgt_mdp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Midlayer height above surface
      if(opt_crr == "hgt_rfr") hgt_rfr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Reference height (i.e., 10 m) at which surface winds are evaluated for dust mobilization
      if(opt_crr == "hgt_zpd_dps_cmd_ln") hgt_zpd_dps_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Zero plane displacement height
      if(opt_crr == "hgt_zpd_mbl") hgt_zpd_mbl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Zero plane displacement height
      if(opt_crr == "hrz_flg") hrz_flg=true; // [flg] Print size-resolved optical properties at debug wavelength
      if(opt_crr == "hxg_flg") hxg_flg=true; // [flg] Aspherical particles are hexagonal prisms
      if(opt_crr == "vts_flg") vts_flg=true; // [flg] Apply equal-V/S approximation for aspherical optical properties
      if(opt_crr == "lat_dgr") lat_dgr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [dgr] Latitude
      if(opt_crr == "lbl_tst") lbl_tst=opt_sng; // [sng] Name of line-by-line test
      if(opt_crr == "lgn_nbr") lgn_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Order of phase function Legendre expansion
      if(opt_crr == "lnd_frc_dry") lnd_frc_dry=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Dry land fraction
      if(opt_crr == "lng_foo") lng_foo=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Intrinsic long temporary variable
      if(opt_crr == "mie_cnv_eps") mie_cnv_eps=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Mie convergence precision
      if(opt_crr == "mie_flg") mie_flg=true; // [flg] Perform mie scattering calculation
      if(opt_crr == "mno_lng_dps_cmd_ln") mno_lng_dps_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Monin-Obukhov length
      if(opt_crr == "mss_frc_cly") mss_frc_cly=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Mass fraction clay
      if(opt_crr == "mss_frc_snd") mss_frc_snd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Mass fraction sand
      if(opt_crr == "mss_frc_cor") mss_frc_cor=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Mass fraction in core
      if(opt_crr == "mss_frc_ncl") mss_frc_ncl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Mass fraction in inclusion
      if(opt_crr == "msv_gnd") msv_gnd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Bare ground emissivity
      if(opt_crr == "msv_snw") msv_snw=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Snow emissivity
      if(opt_crr == "ngl_dbg_dgr") ngl_dbg_dgr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [dgr] Debugging angle
      if(opt_crr == "ngl_nbr") ngl_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of angles in Mie computation
      if(opt_crr == "no_abc_flg") abc_flg=false; // [flg] Alphabetize output with ncks 
      if(opt_crr == "no_bch_flg") bch_flg=false; // [flg] Batch behavior (do not print progress to screen)
      if(opt_crr == "no_hrz_flg") hrz_flg=false; // [flg] Print size-resolved optical properties at debug wavelength
      if(opt_crr == "no_mie_flg") mie_flg=false; // [flg] Perform mie scattering calculation
      if(opt_crr == "no_phz_flg") phz_flg=false; // [flg] Perform phase function diagnostics
      if(opt_crr == "no_wrn_ntp_flg") wrn_ntp_flg=false; // [flg] Print WARNINGs from ntp_vec()
      if(opt_crr == "oro") oro=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Orography: ocean=0.0, land=1.0, sea ice=2.0
      if(opt_crr == "phz_flg") phz_flg=true; // [flg] Perform phase function diagnostics
      if(opt_crr == "pnt_typ_idx") pnt_typ_idx=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [idx] Plant type index 
      if(opt_crr == "prs_mdp") prs_mdp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [Pa] Environmental pressure
      if(opt_crr == "prs_ntf") prs_ntf=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [Pa] Environmental surface pressure
      if(opt_crr == "psd_typ") psd_typ=psd_cls::opt2abb(opt_sng); // [sng] Particle size distribution type
      if(opt_crr == "q_H2O_vpr") q_H2O_vpr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [kg kg-1] Specific humidity
      if(opt_crr == "rds_ffc_gmm_mcr") rds_ffc_gmm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Effective radius of Gamma distribution
      if(opt_crr == "rds_nma_mcr") rds_nma_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Number median radius analytic, microns
      if(opt_crr == "sfc_spc_cm2xg") sfc_spc_cm2xg=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [cm2 g-1] Specific surface area, CGS
      if(opt_crr == "rds_swa_mcr") rds_swa_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Surface area weighted mean radius analytic, microns
      if(opt_crr == "rds_vma_mcr") rds_vma_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Volume median radius analytic, microns
      if(opt_crr == "rgh_mmn_dps_cmd_ln") rgh_mmn_dps_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Roughness length momentum
      if(opt_crr == "rgh_mmn_ice_std") rgh_mmn_ice_std=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Roughness length over sea ice
      if(opt_crr == "rgh_mmn_mbl") rgh_mmn_mbl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Roughness length momentum for erodible surfaces
      if(opt_crr == "rgh_mmn_smt") rgh_mmn_smt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Smooth roughness length
      if(opt_crr == "rfl_gnd_dff") rfl_gnd_dff=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Diffuse reflectance of ground (beneath snow)
      if(opt_crr == "sfc_typ") sfc_typ=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [idx] LSM surface type (0..28)
      if(opt_crr == "slr_spc_key") slr_spc_key=spc_slr_cls::opt2abb(opt_sng); // [sng] Solar spectrum string
      if(opt_crr == "slr_zen_ngl_cos") slr_zen_ngl_cos=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Cosine solar zenith angle
      if(opt_crr == "snw_hgt_lqd") snw_hgt_lqd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m] Equivalent liquid water snow depth
      if(opt_crr == "soi_typ") soi_typ=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [idx] LSM soil type (1..5)
      if(opt_crr == "spc_wgt_sng") spc_wgt_sng=opt_sng; // [sng] Spectral weight
      if(opt_crr == "sz_dbg_mcr") sz_dbg_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Debugging size
      if(opt_crr == "sz_grd_typ") sz_grd_sng=sz_grd_cls::opt2abb(opt_sng); // Type of size grid
      if(opt_crr == "ngl_sng") ngl_sng=opt_sng; // [sng] Angle grid type
      if(opt_crr == "sz_mnm_mcr") sz_mnm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Minimum size in distribution
      if(opt_crr == "sz_mxm_mcr") sz_mxm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Maximum size in distribution
      if(opt_crr == "sz_nbr") sz_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of particle size bins
      if(opt_crr == "sz_prm_rsn_usr_spc") sz_prm_rsn_usr_spc=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m m-1] Size parameter resolution, user specified
      if(opt_crr == "thr_nbr") thr_nbr=static_cast<int>(std::strtol(opt_sng.c_str(),(char **)NULL,10)); // [nbr] Thread number
      if(opt_crr == "tm_dlt") tm_dlt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [s] Timestep
      if(opt_crr == "tpt_bbd_wgt") tpt_bbd_wgt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Blackbody temperature of radiation
      if(opt_crr == "tpt_gnd") tpt_gnd=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Ground temperature
      if(opt_crr == "tpt_ice") tpt_ice=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Ice temperature
      if(opt_crr == "tpt_mdp") tpt_mdp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Environmental temperature
      if(opt_crr == "tpt_soi") tpt_soi=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Soil temperature
      if(opt_crr == "tpt_sst") tpt_sst=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Sea surface temperature
      if(opt_crr == "tpt_vgt") tpt_vgt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K] Vegetation temperature
      if(opt_crr == "tst_sng") tst_sng=opt_sng; // [sng] Name of test to perform
      if(opt_crr == "var_ffc_gmm") var_ffc_gmm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Effective variance of Gamma distribution
      if(opt_crr == "rds_frc_cor") rds_frc_cor=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Radius fraction in core
      if(opt_crr == "rds_frc_ncl") rds_frc_ncl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Radius fraction in inclusion
      if(opt_crr == "vlm_frc_cor") vlm_frc_cor=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Volume fraction in core
      if(opt_crr == "vlm_frc_ncl") vlm_frc_ncl=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Volume fraction in inclusion
      if(opt_crr == "vmr_CO2") vmr_CO2=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [mlc mlc-1] Volume mixing ratio of CO2
      if(opt_crr == "vmr_HNO3_gas") vmr_HNO3_gas=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [mlc mlc-1] Volume mixing ratio of gaseous HNO3
      if(opt_crr == "vwc_sfc") vwc_sfc=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m3 m-3] Volumetric water content
      if(opt_crr == "wbl_shp") wbl_shp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Weibull distribution shape parameter
      if(opt_crr == "wnd_frc_dps_cmd_ln") wnd_frc_dps_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m s-1]
      if(opt_crr == "wnd_mrd_mdp") wnd_mrd_mdp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m s-1] Surface layer meridional wind speed
      if(opt_crr == "wnd_znl_mdp") wnd_znl_mdp=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [m s-1] Surface layer zonal wind speed
      if(opt_crr == "wrn_ntp_flg") wrn_ntp_flg=true; // [flg] Print WARNINGs from ntp_vec()
      if(opt_crr == "wvl_dbg_mcr") wvl_dbg_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Debugging wavelength
      if(opt_crr == "wvl_dlt_mcr") wvl_dlt_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Bandwidth
      if(opt_crr == "wvl_grd_typ") wvl_grd_sng=wvl_grd_cls::opt2abb(opt_sng); // [sng] Wavelength grid type
      if(opt_crr == "wvl_mdp_mcr") wvl_mdp_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Midpoint wavelength
      if(opt_crr == "wvl_mnm_mcr") wvl_mnm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Minimum wavelength
      if(opt_crr == "wvl_mxm_mcr") wvl_mxm_mcr=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [um] Maximum wavelength
      if(opt_crr == "wvl_nbr") wvn_nbr=wvl_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of output wavelength bands
      if(opt_crr == "wvn_dlt_xcm") wvn_dlt_xcm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [cm-1] Bandwidth
      if(opt_crr == "wvn_mdp_xcm") wvn_mdp_xcm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [cm-1] Midpoint wavenumber
      if(opt_crr == "wvn_mnm_xcm") wvn_mnm_xcm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [cm-1] Minimum wavenumber
      if(opt_crr == "wvn_mxm_xcm") wvn_mxm_xcm=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [cm-1] Maximum wavenumber
      if(opt_crr == "wvn_nbr") wvn_nbr=wvl_nbr=std::strtol(opt_sng.c_str(),(char **)NULL,10); // [nbr] Number of output wavenumber bands
      if(opt_crr == "xpt_dsc") xpt_dsc=opt_sng; // [sng] Experiment description
      if(opt_crr == "fdg_val"){
	fdg_flg=true; // [flg] Tune extinction of particular band
	fdg_val=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Tuning factor for all bands
	assert(fdg_val >= 0.0);
      } // endif fdg_flg
      if(opt_crr == "ss_alb"){
	ss_alb_flg=true; // [flg] Manually set single scattering albedo
	ss_alb_cmd_ln=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [frc] Single scattering albedo
	assert(ss_alb_cmd_ln >= 0.0 && ss_alb_cmd_ln <= 1.0);
      } // endif ss_alb_flg
      if(opt_crr == "idx_rfr_cor_usr"){
	if(optarg) idx_rfr_cor_usr=sng2cpx(opt_sng); // [frc] Refractive index of core
	idx_rfr_cor_usr_flg=true; // [flg] Refractive index of core is user-specified
      } // end if "idx_rfr_cor_usr"
      if(opt_crr == "idx_rfr_mdm_usr"){
	if(optarg) idx_rfr_mdm_usr=sng2cpx(opt_sng); // [frc] Refractive index of medium
	idx_rfr_mdm_usr_flg=true; // [flg] Refractive index of medium is user-specified
      } // end if "idx_rfr_mdm_usr"
      if(opt_crr == "idx_rfr_mnt_usr"){
	if(optarg) idx_rfr_mnt_usr=sng2cpx(opt_sng); // [frc] Refractive index of mantle
	idx_rfr_mnt_usr_flg=true; // [flg] Refractive index of mantle is user-specified
      } // end if "idx_rfr_mnt_usr"
      if(opt_crr == "idx_rfr_mtx_usr"){
	if(optarg) idx_rfr_mtx_usr=sng2cpx(opt_sng); // [frc] Refractive index of medium
	idx_rfr_mtx_usr_flg=true; // [flg] Refractive index of medium is user-specified
      } // end if "idx_rfr_mtx_usr"
      if(opt_crr == "idx_rfr_ncl_usr"){
	if(optarg) idx_rfr_ncl_usr=sng2cpx(opt_sng); // [frc] Refractive index of inclusion
	idx_rfr_ncl_usr_flg=true; // [flg] Refractive index of inclusion is user-specified
      } // end if "idx_rfr_ncl_usr"
      if(opt_crr == "idx_rfr_prt_usr"){
	if(optarg) idx_rfr_prt_usr=sng2cpx(opt_sng); // [frc] Refractive index of particle
	idx_rfr_prt_usr_flg=true; // [flg] Refractive index of particle is user-specified
      } // end if "idx_rfr_prt_usr"
      if(opt_crr == "mie_slv"){
	// [sng] Mie solver to use (BoH83 or Wis79)
	if((opt_sng == "BoH83") || (opt_sng == "Wis79")) slv_sng=opt_sng; else err_prn(prg_nm,sbr_nm,opt_sng+ " is not a optics valid solver");
      } // endif "mie_slv"
      if(opt_crr == "slf_tst_typ"){
	// Decode self-test type
	if((slf_tst_sng[BoH83].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("BoH83") != std::string::npos)){ 
	  slf_tst_typ=BoH83; // [enm] Self-test type
	  //ngl_sng="regular"; // [sng] Angle grid type
	}else if((slf_tst_sng[ViC98].find(opt_sng) != std::string::npos) || 
		 (opt_sng.find("ViC98") != std::string::npos)){
	  slf_tst_typ=ViC98; // [enm] Self-test type
	}else{
	  err_prn(prg_nm,sbr_nm,"Unknown slf_tst_typ");
	} // end else
	slf_tst_flg=true; // [flg] Perform self-test
      } // end if "slf_tst_typ"
      if(opt_crr == "ffc_mdm_typ"){
	// Get type of effective medium approximation
	ffc_mdm_flg=true; // [flg] Employ effective medium approximation
	if((ffc_mdm_sng[ffc_mdm_brg].find(opt_sng) != std::string::npos) || 
	   (opt_sng.find("brg") != std::string::npos)) 
	  ffc_mdm_typ=ffc_mdm_brg; // [enm] Effective medium type
	else if((ffc_mdm_sng[ffc_mdm_mxg].find(opt_sng) != std::string::npos) || 
		(opt_sng.find("mxg") != std::string::npos)) 
	  ffc_mdm_typ=ffc_mdm_mxg; // [enm] Effective medium type
	else if((ffc_mdm_sng[ffc_mdm_pmr].find(opt_sng) != std::string::npos) || 
		(opt_sng.find("pmr") != std::string::npos)) 
	  ffc_mdm_typ=ffc_mdm_pmr; // [enm] Effective medium type
	else if((ffc_mdm_sng[ffc_mdm_vlw].find(opt_sng) != std::string::npos) || 
		(opt_sng.find("vlw") != std::string::npos)) 
	  ffc_mdm_typ=ffc_mdm_vlw; // [enm] Effective medium type
	else if((ffc_mdm_sng[ffc_mdm_nil].find(opt_sng) != std::string::npos) || 
		(opt_sng.find("nil") != std::string::npos)){
	  ffc_mdm_typ=ffc_mdm_nil; // [enm] Effective medium type
	  ffc_mdm_flg=false;} // [flg] Employ effective medium approximation
	else err_prn(prg_nm,sbr_nm,"Unknown ffc_mdm_typ");
      } // end if "ffc_mdm_typ"
    } // opt != 0
    switch(opt){
    case 0: // Long options have already been processed, return
      break;
    case '3': // Prescribe netCDF3 output storage format
      fl_out_fmt=NC_FORMAT_CLASSIC;
      break;
    case '4': // Catch-all to prescribe output storage format */
      if(opt_sng == "64bit") fl_out_fmt=NC_FORMAT_64BIT; else fl_out_fmt=NC_FORMAT_NETCDF4;
      break;
    case 'D': // The debugging level. Default is 0 
      if(optarg) dbg_lvl=static_cast<unsigned short int>(std::strtoul(opt_sng.c_str(),(char **)NULL,10)); else dbg_lvl=dbg_fl;
      break;
    case 'i': // Input file name. Default uses file in aerosol database
      fl_idx_rfr_prt=opt_sng; // [sng] File containing refractive indices of particle
      break;
    case 's': // Solar constant. Default is 1367 W m-2
      slr_cst=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [W m-2]
      break;
    case 't': // Particle temperature. Default is 273.15 K 
      tpt_prt=static_cast<prc_cmp>(std::strtod(opt_sng.c_str(),(char **)NULL)); // [K]
      break;
    case 'h': // Print proper usage 
      //      (void)usg_prn(opt_sht_lst);
      return EXIT_SUCCESS;
    case 'o': // Get output file name. Default is std::cout 
      fl_out=opt_sng; // [sng] Output file
      break;
    case 'v': // Print CVS information
      std::cerr << CVS_Header << std::endl;
      return EXIT_SUCCESS;
    default: // Print proper usage 
      // (void)usg_prn(opt_sht_lst);
      return EXIT_FAILURE;
    } // end switch 
  } // end while loop 
  
  // Process positional arguments  
  if(optind < argc){
    int psn_arg_nbr=argc-optind; // [nbr] Positional argument number
    if(psn_arg_nbr > 2){
      err_prn(prg_nm,sbr_nm,"Too many positional arguments");
    }else if(psn_arg_nbr == 1){
      // fxm: TODO#17 conflicts with -o if present
      fl_out=argv[optind++]; // [sng] Output file
    }else if(psn_arg_nbr == 2){
      fl_idx_rfr_prt=argv[optind++]; // [sng] File containing refractive indices of particle
      // fxm: TODO#17 conflicts with -o if present
      fl_out=argv[optind++]; // [sng] Output file
    } // end else
  } // end if  

  // Compute quantities that might depend on command line input
  (void)fio::data_path_set(drc_dat); // [sng] Data directory
  // Prepend user-specified path, if any, to input data file names
  if(fl_idx_rfr_cor.length() > 0) fl_idx_rfr_cor=drc_pfx(fio::data_path_get(),fl_idx_rfr_cor); // [sng] File or function for refractive indices of core
  if(fl_idx_rfr_mdm.length() > 0) fl_idx_rfr_mdm=drc_pfx(fio::data_path_get(),fl_idx_rfr_mdm); // [sng] File or function for refractive indices of medium
  if(fl_idx_rfr_mnt.length() > 0) fl_idx_rfr_mnt=drc_pfx(fio::data_path_get(),fl_idx_rfr_mnt); // [sng] File or function for refractive indices of mantle
  if(fl_idx_rfr_mtx.length() > 0) fl_idx_rfr_mtx=drc_pfx(fio::data_path_get(),fl_idx_rfr_mtx); // [sng] File or function for refractive indices of matrix
  if(fl_idx_rfr_ncl.length() > 0) fl_idx_rfr_ncl=drc_pfx(fio::data_path_get(),fl_idx_rfr_ncl); // [sng] File or function for refractive indices of inclusion
  if(fl_idx_rfr_prt.length() > 0) fl_idx_rfr_prt=drc_pfx(fio::data_path_get(),fl_idx_rfr_prt); // [sng] File or function for refractive indices of particle
  // Prepend user-specified path, if any, to output data file names
  if(fl_out.length() > 0) fl_out=drc_pfx(drc_out,fl_out); // [sng] Output file
  const prc_cmp one_plus_mie_cnv_eps(1.0+mie_cnv_eps); // [frc] Mie convergence precision factor

  // Main code  
  if(dbg_lvl > dbg_sbr) dbg_lvl_tst();
  if(!wrn_ntp_flg) wrn_prn(sbr_nm,"Turning off verbose error reports from ntp_vec()");

  // Initialize thread information
  thr_nbr=openmp_ini(thr_nbr);

  /* Spectral grid: Program works internally in wavelength space
     User may specify grid using either wavelength in microns or wavenumber in cm-1 
     Grid boundaries specified as midpoint + full range or as mininum and maximum
     Grid resolution is set with wvl_nbr or wvn_nbr, as appropriate
     --wvn_dlt_xcm [cm-1] Bandwidth
     --wvn_mdp_xcm [cm-1] Midpoint wavenumber
     --wvn_mnm_xcm [cm-1] Minimum wavenumber
     --wvn_mxm_xcm [cm-1] Maximum wavenumber
     --wvl_dlt_mcr [um] Bandwidth
     --wvl_mdp_mcr [um] Midpoint wavelength
     --wvl_mnm_mcr [um] Minimum wavelength
     --wvl_mxm_mcr [um] Maximum wavelength
     Translate input to wvl_mxm_mcr,wvl_mxm,wvl_mnm_mcr,wvl_mnm before creating grid */
  assert(wvn_mnm_xcm != 0.0 && wvn_mxm_xcm != 0.0);
  if(wvn_mdp_xcm != cmd_ln_dfl){
    if(wvn_dlt_xcm == cmd_ln_dfl) wvn_dlt_xcm=1.0; // [cm-1] Bandwidth
    // Place wvn_mdp_xcm at bin midpoint if wvn_nbr is odd, at bin interface if wvn_nbr is even
    wvn_mnm_xcm=wvn_mdp_xcm-0.5*wvn_dlt_xcm; // [cm-1] Minimum wavenumber
    wvn_mxm_xcm=wvn_mdp_xcm+0.5*wvn_dlt_xcm; // [cm-1] Maximum wavenumber
  } // endif
  if(wvn_mnm_xcm != cmd_ln_dfl) wvl_mxm_mcr=1.0e6*0.01/wvn_mnm_xcm; // [um] Maximum wavelength
  if(wvn_mxm_xcm != cmd_ln_dfl) wvl_mnm_mcr=1.0e6*0.01/wvn_mxm_xcm; // [um] Minimum wavelength
  if(wvl_mdp_mcr != cmd_ln_dfl){
    if(wvl_dlt_mcr == cmd_ln_dfl) wvl_dlt_mcr=0.1; // [um] Bandwidth
    // Place wvl_mdp_mcr at bin midpoint if wvl_nbr is odd, at bin interface if wvl_nbr is even
    wvl_mnm_mcr=wvl_mdp_mcr-0.5*wvl_dlt_mcr; // [um] Minimum wavelength
    wvl_mxm_mcr=wvl_mdp_mcr+0.5*wvl_dlt_mcr; // [um] Maximum wavelength
  } // endif

  // Perform error checking/correction on inputs
  if(fdg_flg && ss_alb_flg) err_prn(prg_nm,sbr_nm," fdg_flg and ss_alb_flg are both specified");

  if(sz_dbg_mcr < sz_mnm_mcr || sz_dbg_mcr > sz_mxm_mcr) sz_dbg_mcr=0.5*(sz_mnm_mcr+sz_mxm_mcr); // [um] Debugging size
  if(wvl_dbg_mcr < wvl_mnm_mcr || wvl_dbg_mcr > wvl_mxm_mcr) wvl_dbg_mcr=0.5*(wvl_mnm_mcr+wvl_mxm_mcr); // [um] Debugging wavelength
  if(dsd_dbg_mcr < dsd_mnm_mcr || dsd_dbg_mcr > dsd_mxm_mcr) dsd_dbg_mcr=0.5*(dsd_mnm_mcr+dsd_mxm_mcr); // [um] Debugging size for raindrop

  // Set parameters for test cases
  if(slf_tst_flg){
    wvl_nbr=1; // [nbr] Number of output wavelength bands
    bnd_nbr=1; // [nbr] Number of sub-bands per output band
    sz_nbr=1; // [nbr] Number of particle size bins
    ngl_nbr=11; // [nbr] Number of polar angles in one hemisphere
    if(slf_tst_typ == BoH83){ // Bohren & Huffman (1983) p. 482
      /* NB: mie --slf_tst=BoH83 should yield same results as command line implementation of same test:
	 mie --dbg=3 --wvl_mnm=0.6327 --wvl_mxm=0.6329 --wvl_nbr=1 --bnd_nbr=1 --sz_grd=linear --sz_mnm=0.524999 --sz_mxm=0.525001 --sz_nbr=1 --ngl_nbr=11 --idx_rfr_mdm=1.0 --idx_rfr_prt=1.55 ${DATA}/mie/mie.nc; ncks -C -v ext_fsh,sca_fsh,abs_fsh,bck_fsh ${DATA}/mie/mie.nc */
      idx_rfr_mdm_usr_flg=true; // [flg] Refractive index of medium is user-specified
      idx_rfr_mdm_usr=std::complex<prc_cmp>(1.0,0.0); // [frc] Refractive index of medium BoH83 p. 482
      idx_rfr_prt_usr_flg=true; // [flg] Refractive index of particle is user-specified
      idx_rfr_prt_usr=std::complex<prc_cmp>(1.55,0.0); // [frc] Refractive index of particle BoH83 p. 482
      sz_mnm_mcr=0.524999; // [um] Minimum size BoH83 p. 482
      sz_mxm_mcr=0.525001; // [um] Maximum size BoH83 p. 482
      wvl_mnm_mcr=0.6328; // [um] Minimum wavelength BoH83 p. 482
      wvl_mxm_mcr=0.6328; // [um] Maximum wavelength BoH83 p. 482
      if(coat_flg){ // Bohren & Huffman (1983) p. 489
      // NB: mie --slf_tst=BoH83 --coat should yield same results as command line implementation of same test:
      // mie --dbg=3 --coat --dbg=3 --rds_frc_cor=0.02729449321628092577 --wvl_mnm=2.999 --wvl_mxm=3.001 --wvl_nbr=1 --bnd_nbr=1 --sz_grd=linear --sz_mnm=0.6249 --sz_mxm=0.6251 --sz_nbr=1 --ngl_nbr=11 --idx_rfr_cor=1.59+0.66i --idx_rfr_mdm=1.0 --idx_rfr_mnt=1.409+0.1747i ${DATA}/mie/mie.nc; ncks -C -H -v ext_fsh,sca_fsh,abs_fsh,bck_fsh ${DATA}/mie/mie.nc
	idx_rfr_cor_usr_flg=true; // [flg] Refractive index of core is user-specified
	idx_rfr_cor_usr=std::complex<prc_cmp>(1.59,0.66); // [frc] Refractive index of core BoH83 p. 489
	idx_rfr_mdm_usr_flg=true; // [flg] Refractive index of medium is user-specified
	idx_rfr_mdm_usr=std::complex<prc_cmp>(1.0,0.0); // [frc] Refractive index of medium
	idx_rfr_mnt_usr_flg=true; // [flg] Mantle refractive index is user-specified
	idx_rfr_mnt_usr=std::complex<prc_cmp>(1.409,0.1747); // [frc] Refractive index of mantle BoH83 p. 489
	// BoH83 coat test must use rds_mnt[sz_idx]=6.265e-6; // [m] Radius of mantle
	rds_frc_cor=0.02729449321628092577; // [frc] Radius fraction in core
	sz_mnm_mcr=6.2649; // [um] Minimum size BoH83 p. 489
	sz_mxm_mcr=6.2651; // [um] Maximum size BoH83 p. 489
	wvl_mnm_mcr=2.999; // [um] Minimum wavelength BoH83 p. 489
	wvl_mxm_mcr=3.001; // [um] Maximum wavelength BoH83 p. 489
      } // endif coat_flg
    } // not BoH83
    if(slf_tst_typ == ViC98){ // Videen & Chylek (1998) p. 1
      // NB: mie --slf_tst=ViC98 should yield same results as command line implementation of same test:
      // mie --dbg=3 --ffc_mdm=mxg --vlm_frc_ncl=0.1 --wvl_mnm=0.999 --wvl_mxm=1.001 --wvl_nbr=1 --bnd_nbr=1 --sz_grd=linear --sz_mnm=0.999 --sz_mxm=1.001 --sz_nbr=1 --ngl_nbr=11 --idx_rfr_mdm=1.0  --idx_rfr_mtx=1.335 --idx_rfr_ncl=1.94+0.66i ${DATA}/mie/mie.nc; ncks -C -m -v ext_fsh,sca_fsh,abs_fsh,bck_fsh,idx_rfr_ffc_mxg_rl,idx_rfr_ffc_mxg_img,idx_rfr_ffc_brg_rl,idx_rfr_ffc_brg_img ${DATA}/mie/mie.nc
      ffc_mdm_typ=ffc_mdm_mxg; // I [enm] Effective medium type
      idx_rfr_mdm_usr_flg=true; // [flg] Refractive index of medium is user-specified
      idx_rfr_mdm_usr=std::complex<prc_cmp>(1.0,0.0); // [frc] Refractive index of medium
      idx_rfr_mtx_usr_flg=true; // [flg] Refractive index of matrix is user-specified
      idx_rfr_mtx_usr=std::complex<prc_cmp>(1.335,0.0); // [frc] Refractive index of matrix ViC98 p. 5 Tbl. 3
      idx_rfr_ncl_usr_flg=true; // [flg] Inclusion refractive index is user-specified
      idx_rfr_ncl_usr=std::complex<prc_cmp>(1.94,0.66); // [frc] Refractive index of inclusion ViC98 p. 5 Tbl. 3
      sz_mnm_mcr=0.999; // [um] Minimum size ViC98 p. 5 Tbl. 3
      sz_mxm_mcr=1.001; // [um] Maximum size ViC98 p. 5 Tbl. 3
      wvl_mnm_mcr=0.999; // [um] Minimum wavelength ViC98 p. 5 Tbl. 3
      wvl_mxm_mcr=1.001; // [um] Maximum wavelength ViC98 p. 5 Tbl. 3
      vlm_frc_ncl=0.1; // [frc] Volume fraction in core ViC98 p. 5 Tbl. 3
      // vlm_frc_ncl=0.01; // [frc] Volume fraction in core ViC98 p. 5 Tbl. 3
      // vlm_frc_ncl=0.001; // [frc] Volume fraction in core ViC98 p. 5 Tbl. 3
      // vlm_frc_ncl=0.0001; // [frc] Volume fraction in core ViC98 p. 4 Fgr. 3
      // rds_frc_ncl=0.464; // [frc] Radius fraction in core ViC98 p. 3 Tbl. 2
      // rds_frc_ncl=0.215; // [frc] Radius fraction in core ViC98 p. 3 Tbl. 1
    } // not ViC98
  } // not slf_tst_flg
  
  // Convert input data to SI units
  prc_cmp wvl_mnm=wvl_mnm_mcr*1.0e-6; // [um] -> [m] Minimum wavelength 
  prc_cmp wvl_mxm=wvl_mxm_mcr*1.0e-6; // [um] -> [m] Maximum wavelength
  prc_cmp wvl_dbg=wvl_dbg_mcr*1.0e-6; // [um] -> [m] Debugging wavelength
  prc_cmp sz_mnm=sz_mnm_mcr*1.0e-6; // [um] -> [m] Minimum size
  prc_cmp sz_mxm=sz_mxm_mcr*1.0e-6; // [um] -> [m] Maximum size
  prc_cmp sz_dbg=sz_dbg_mcr*1.0e-6; // [um] -> [m] Debugging size

  // Error checking/correction on input
  if(wvl_mnm_mcr == wvl_mxm_mcr && wvl_nbr*bnd_nbr != 1) wrn_prn(prg_nm,sbr_nm,"wvl_mnm_mcr == wvl_mxm_mcr, may cause problems with flx_wgt_frc");
  // Setting minimum size equals maximum size could be made to work but would open up many pointless ambiguities such as how is number concentration determined (since dst_nbr*sz_dlt = 0.0)? Should number concentrations of 0.0 work? 
  if(sz_mnm_mcr == sz_mxm_mcr) err_prn(prg_nm,sbr_nm,"sz_mnm_mcr == sz_mxm_mcr, degenerate size distribution HINT: to consider mono-disperse size distributions, set sz_mnm_mcr = sz_mxm_mcr + epsilon or use narrow standard deviation in size distribution, e.g., --gsl_anl = 1.01");

  // fxm: Attach aer_cls object to psd_cls object so density and file source is always defined?
  // Instantiate aerosol components
  aer_cls cmp_cor(cmp_sng_cor,dns_cor,tpt_prt,fl_idx_rfr_cor,idx_rfr_cor_usr_flg,idx_rfr_cor_usr); // [aer] Core component
  aer_cls cmp_mdm(cmp_sng_mdm,dns_mdm,tpt_mdp,fl_idx_rfr_mdm,idx_rfr_mdm_usr_flg,idx_rfr_mdm_usr); // [aer] Medium component
  aer_cls cmp_mnt(cmp_sng_mnt,dns_mnt,tpt_prt,fl_idx_rfr_mnt,idx_rfr_mnt_usr_flg,idx_rfr_mnt_usr); // [aer] Mantle component
  aer_cls cmp_mtx(cmp_sng_mtx,dns_mtx,tpt_prt,fl_idx_rfr_mtx,idx_rfr_mtx_usr_flg,idx_rfr_mtx_usr); // [aer] Matrix component
  aer_cls cmp_ncl(cmp_sng_ncl,dns_ncl,tpt_prt,fl_idx_rfr_ncl,idx_rfr_ncl_usr_flg,idx_rfr_ncl_usr); // [aer] Inclusion component
  fl_idx_rfr_cor=cmp_cor.fl_idx_rfr_get(); // [sng] File containing refractive indices of core
  fl_idx_rfr_mdm=cmp_mdm.fl_idx_rfr_get(); // [sng] File containing refractive indices of medium
  fl_idx_rfr_mnt=cmp_mnt.fl_idx_rfr_get(); // [sng] File containing refractive indices of mantle
  fl_idx_rfr_mtx=cmp_mtx.fl_idx_rfr_get(); // [sng] File containing refractive indices of matrix
  fl_idx_rfr_ncl=cmp_ncl.fl_idx_rfr_get(); // [sng] File containing refractive indices of inclusion
  dns_cor=cmp_cor.dns_get(); // [kg m-3] Density of core
  dns_mdm=cmp_mdm.dns_get(); // [kg m-3] Density of medium
  dns_mnt=cmp_mnt.dns_get(); // [kg m-3] Density of mantle
  dns_mtx=cmp_mtx.dns_get(); // [kg m-3] Density of matrix
  dns_ncl=cmp_ncl.dns_get(); // [kg m-3] Density of inclusion

  // Determine relative weights of aerosol components
  /* Generally, optics models rely on volume fractions not mass fractions
     However, usually more convenient to specify mass or radius fraction
     Convert mss_frc or rds_frc when specified, to vlm_frc */
  if(mss_frc_cor != cmd_ln_dfl) vlm_frc_cor=mss_frc_cor*cmp_mnt.dns_get()/cmp_cor.dns_get(); // [frc] Volume fraction in core
  if(mss_frc_ncl != cmd_ln_dfl) vlm_frc_ncl=mss_frc_ncl*cmp_mtx.dns_get()/cmp_ncl.dns_get(); // [frc] Volume fraction in inclusion
  if(rds_frc_cor != cmd_ln_dfl) vlm_frc_cor=std::pow(rds_frc_cor,PRC_CMP(3.0)); // [frc] Volume fraction in core
  if(rds_frc_ncl != cmd_ln_dfl) vlm_frc_ncl=std::pow(rds_frc_ncl,PRC_CMP(3.0)); // [frc] Volume fraction in inclusion

  // (re-)Derive diagnostic mss_frc and rds_frc from vlm_frc
  mss_frc_cor=vlm_frc_cor*cmp_cor.dns_get()/(vlm_frc_cor*cmp_cor.dns_get()+(1.0-vlm_frc_cor)*cmp_mnt.dns_get()); // [frc] Mass fraction in cor
  mss_frc_ncl=vlm_frc_ncl*cmp_ncl.dns_get()/(vlm_frc_ncl*cmp_ncl.dns_get()+(1.0-vlm_frc_ncl)*cmp_mtx.dns_get()); // [frc] Mass fraction in inclusion
  rds_frc_cor=std::pow(vlm_frc_cor,PRC_CMP(1.0)/PRC_CMP(3.0)); // [frc] Radius fraction in core
  rds_frc_ncl=std::pow(vlm_frc_ncl,PRC_CMP(1.0)/PRC_CMP(3.0)); // [frc] Radius fraction in inclusion

  /* Generic properties are taken to be those of particle component
     Particle component properties may be derived from other components
     For instance, coated aerosols have an exact Mie solution,
     and specific optical properties (e.g., mass absorption coefficient)
     must be normalized to density of core/mantle system.
     Effective medium approximation (EMA) specific optical properties
     should be normalized to the density of the inclusion/matrix system. */

  // Allow user-specified particle density to override derived mixture
  if(coat_flg && dns_prt == 0.0){
    const prc_cmp spc_heat_cor(cmp_cor.spc_heat_get()); // [J kg-1 K-1] Specific heat capacity of core
    const prc_cmp mmw_cor(cmp_cor.mmw_get()); // [kg mol-1] Mean molecular weight of core
    const prc_cmp spc_heat_mnt(cmp_mnt.spc_heat_get()); // [J kg-1 K-1] Specific heat capacity of mantle
    const prc_cmp mmw_mnt(cmp_mnt.mmw_get()); // [kg mol-1] Mean molecular weight of mantle

    wrn_prn(sbr_nm,"Computing particle physical properties (e.g., density) from specified blend of core with mantle. This automatically makes output specific (per-unit mass) optical properties consistent with multi-component aerosol (MCA).");

    // Density of particle component equals volume mean density of core/mantle
    dns_prt=vlm_frc_cor*dns_cor+(1.0-vlm_frc_cor)*dns_mnt; // [kg m-3] Density of particle
    // Mean molecular weight of particle component equals mass-mean ...
    mmw_prt=mss_frc_cor*mmw_cor+(1.0-mss_frc_cor)*mmw_mnt; // [J kg-1 K-1] Specific heat capacity of particle
    // Specific heat capacity of particle component equals mass-mean ...
    spc_heat_prt=mss_frc_cor*spc_heat_cor+(1.0-mss_frc_cor)*spc_heat_mnt; // [J kg-1 K-1] Specific heat capacity of particle
  } // !coat_flg || dns_prt != 0.0

  // Allow user-specified particle density to override derived mixture
  if(ffc_mdm_flg && dns_prt == 0.0){
    const prc_cmp spc_heat_ncl(cmp_ncl.spc_heat_get()); // [J kg-1 K-1] Specific heat capacity of inclusion
    const prc_cmp mmw_ncl(cmp_ncl.mmw_get()); // [kg mol-1] Mean molecular weight of inclusion
    const prc_cmp spc_heat_mtx(cmp_mtx.spc_heat_get()); // [J kg-1 K-1] Specific heat capacity of matrix
    const prc_cmp mmw_mtx(cmp_mtx.mmw_get()); // [kg mol-1] Mean molecular weight of matrix

    wrn_prn(sbr_nm,"Computing particle physical properties (e.g., density) from specified blend of inclusion with matrix. This automatically makes output specific (per-unit mass) optical properties consistent with multi-component aerosol (MCA).");

    // Density of particle component equals volume mean density of inclusion/matrix
    dns_prt=vlm_frc_ncl*dns_ncl+(1.0-vlm_frc_ncl)*dns_mtx; // [kg m-3] Density of particle
    // Mean molecular weight of particle component equals mass-mean ...
    mmw_prt=mss_frc_ncl*mmw_ncl+(1.0-mss_frc_ncl)*mmw_mtx; // [J kg-1 K-1] Specific heat capacity of particle
    // Specific heat capacity of particle component equals mass-mean ...
    spc_heat_prt=mss_frc_ncl*spc_heat_ncl+(1.0-mss_frc_ncl)*spc_heat_mtx; // [J kg-1 K-1] Specific heat capacity of particle
  } // !ffc_mdm_flg || dns_prt != 0.0

  aer_cls cmp_prt(cmp_sng_prt,dns_prt,tpt_prt,fl_idx_rfr_prt,idx_rfr_prt_usr_flg,idx_rfr_prt_usr); // [aer] Particle component
  fl_idx_rfr_prt=cmp_prt.fl_idx_rfr_get(); // [sng] File containing refractive indices of particle
  // Density in aerosol component object already accounts for coating
  dns_prt=cmp_prt.dns_get(); // [kg m-3] Density of particle
  // Do no overwrite MCA specific heat or mean molecular weight, if any
  if(mmw_prt == 0.0) mmw_prt=cmp_prt.mmw_get(); // [kg mol-1] Mean molecular weight of particle
  if(spc_heat_prt == 0.0) spc_heat_prt=cmp_prt.spc_heat_get(); // [J kg-1 K-1] Specific heat capacity of particle

  // Instantiate size distributions
  prc_cmp rds_ffc_gmm; // [m] Effective radius for gamma distribution
  prc_cmp rds_nma_dfl; // [m] Radius number median analytic, default
  rds_ffc_gmm=rds_ffc_gmm_mcr*1.0e-6; // [m] Effective radius for gamma distribution
  /* Mie program works internally with radius 
     Lognormal distributions, e.g., are specified with rds_nma
     For monomodal distributions, user may instead specify 
     --rds_swa_mcr Surface area weighted mean radius analytic, microns (i.e, effective radius)
     --rds_vma_mcr Volume median radius analytic, microns (mass-median radius) 
     --sfc_spc_cm2xg Specific surface area, CGS [cm2 g-1] 
     --dmt_nma_mcr Number median diameter analytic, microns
     --dmt_swa_mcr Surface area weighted mean diameter analytic, microns (i.e., effective diameter)
     --dmt_vma_mcr Volume median diameter analytic, microns (mass-median diameter) 
     Convert option to equivalent rds_nma before instantiating distribution(s) */
  if(sfc_spc_cm2xg != cmd_ln_dfl){
    prc_cmp sfc_spc_m2xkg(sfc_spc_cm2xg*0.1); // [cm2 g-1 -> m2 kg-1]
    prc_cmp rds_swa_tmp(3.0/(dns_prt*sfc_spc_m2xkg)); // [m]
    rds_nma_mcr=1.0e6*swa2nma(rds_swa_tmp,gsd_anl_dfl); // [m->um] Number median radius analytic, microns
    //    std::cout << "quark foo = " << foo << std::endl;
  } // endif sfc_spc
  if(rds_swa_mcr != cmd_ln_dfl) rds_nma_mcr=swa2nma(rds_swa_mcr,gsd_anl_dfl); // [um] Number median radius analytic, microns
  if(rds_vma_mcr != cmd_ln_dfl) rds_nma_mcr=vma2nma(rds_vma_mcr,gsd_anl_dfl); // [um] Number median radius analytic, microns
  if(dmt_nma_mcr != cmd_ln_dfl) rds_nma_mcr=0.5*dmt_nma_mcr; // [um] Number median radius analytic, microns
  if(dmt_swa_mcr != cmd_ln_dfl) rds_nma_mcr=0.5*swa2nma(dmt_swa_mcr,gsd_anl_dfl); // [um] Number median radius analytic, microns
  if(dmt_vma_mcr != cmd_ln_dfl) rds_nma_mcr=0.5*vma2nma(dmt_vma_mcr,gsd_anl_dfl); // [um] Number median radius analytic, microns
  // Older versions of aerosol model specified bin midpoint as rds_nma
  // Setting drv_rds_nma_flg invokes this older method
  if(drv_rds_nma_flg) rds_nma_dfl=0.5*(sz_mnm+sz_mxm); else rds_nma_dfl=rds_nma_mcr*1.0e-6; // [m] Radius number median analytic
  psd_cls psd_dfl(rds_nma_dfl,gsd_anl_dfl,cnc_nbr_anl_dfl); // [sct] Default size distribution
  psd_cls *psd_lst=new psd_cls[max_lng(psd_nbr,1L)]; // [sct] List of aerosol size distributions

  // Parse aerosol modes
  if(dbg_lvl > dbg_crr) std::cout << "Current number of size distributions is " << psd_lst[0].nst_nbr_get() << std::endl;
  if(psd_nbr > 0){
    // Aerosol modes were specified by user
    for(psd_idx=0;psd_idx<psd_nbr;psd_idx++) rcd+=psd_lst[psd_idx].prs_psd_sng(psd_arg[psd_idx]);
  }else{
    // No modes were specified by user, use default mode
    psd_nbr=1; // [nbr] Number of aerosol modes
    psd_lst=&psd_dfl; // [sct] List of aerosol modes
  } // endif psd_nbr > 0
  // Set size distribution to correct density
  for(psd_idx=0;psd_idx<psd_nbr;psd_idx++) rcd+=psd_lst[psd_idx].dns_prt_set(dns_prt); // [kg m-3] Particle density
  // Set size distribution mass fractions
  rcd+=psd_lst[0].mss_frc_anl_set(psd_lst,psd_nbr); // [frc] Mass fraction analytic

  // Convenient arrays which group like mode parameters together
  prc_cmp *psd_dmt_nma=new prc_cmp[psd_nbr]; // [m] Number median diameter analytic
  prc_cmp *psd_rds_nma=new prc_cmp[psd_nbr]; // [m] Number median radius analytic
  prc_cmp *psd_gsd_anl=new prc_cmp[psd_nbr]; // [frc] Geometric standard deviation
  prc_cmp *psd_cnc_nbr_anl=new prc_cmp[psd_nbr]; // [# m-3] Number concentration analytic
  prc_cmp *psd_mss_frc_anl=new prc_cmp[psd_nbr]; // [frc] Mass fraction analytic
  for(idx=0;idx<psd_nbr;idx++){
    psd_rds_nma[idx]=psd_lst[idx].rds_nma_get(); // [m] Number median radius analytic
    psd_dmt_nma[idx]=psd_lst[idx].dmt_nma_get(); // [m] Number median diameter analytic
    psd_gsd_anl[idx]=psd_lst[idx].gsd_anl_get(); // [frc] Geometric standard deviation
    psd_cnc_nbr_anl[idx]=psd_lst[idx].cnc_nbr_anl_get(); // [# m-3] Number concentration analytic
    psd_mss_frc_anl[idx]=psd_lst[idx].mss_frc_anl_get(); // [frc] Mass fraction analytic
  } // end loop over modes

  // Instantiate size grid
  assert(psd_nbr > 0L);
  sz_grd_cls psd_szgrd(sz_mnm,sz_mxm,sz_nbr,sz_grd_sng); // [m] Size grid object
  const prc_cmp *sz_min=psd_szgrd.sz_min_get(); // [m] Minimum size in bin
  const prc_cmp *sz_max=psd_szgrd.sz_max_get(); // [m] Maximum size in bin
  const prc_cmp *sz_ctr=psd_szgrd.sz_ctr_get(); // [m] Size at bin center
  const prc_cmp *sz_grd=psd_szgrd.sz_grd_get(); // [m] Size grid
  const prc_cmp *sz_dlt=psd_szgrd.sz_dlt_get(); // [m] Width of size bin
  /* 20210418: Conflict between GCC and Clang OpenMP implementations
     GCC disallows const scalars to be specificially enumerated in OpenMP shared ("error: 'sz_idx_dbg' is predetermined 'shared' for 'shared'")
     Clang requires const scalars to be specificially enumerated in OpenMP shared
     Workaround: declare sz_idx_dbg, wvl_bnd_sz_nbr, wvl_idx_dbg as non-const */
  //  const long sz_idx_dbg=vec_val2idx(sz_ctr,sz_nbr,sz_dbg); // [idx] Size bin for debugging
  long sz_idx_dbg=vec_val2idx(sz_ctr,sz_nbr,sz_dbg); // [idx] Size bin for debugging

  // Define radius grid
  prc_cmp *rds_ctr=new prc_cmp[sz_nbr]; // [m] Radius at bin center
  prc_cmp *rds_dlt=new prc_cmp[sz_nbr]; // [m] Width of radius bin
  prc_cmp *rds_grd=new prc_cmp[sz_nbr+1]; // [m] Radius grid
  prc_cmp *rds_max=new prc_cmp[sz_nbr]; // [m] Maximum radius in bin
  prc_cmp *rds_min=new prc_cmp[sz_nbr]; // [m] Minimum radius in bin
  for(idx=0;idx<sz_nbr;idx++){
    rds_ctr[idx]=sz_ctr[idx]; // [m] Radius at bin center
    rds_dlt[idx]=sz_dlt[idx]; // [m] Width of radius bin
    rds_grd[idx]=sz_grd[idx]; // [m] Radius grid
    rds_max[idx]=sz_max[idx]; // [m] Maximum radius in bin
    rds_min[idx]=sz_min[idx]; // [m] Minimum radius in bin
  } // end loop over sz
  rds_grd[sz_nbr]=sz_grd[sz_nbr]; // [m] Radius grid
  const prc_cmp rds_ctr_ctr(0.5*(rds_grd[0]+rds_grd[sz_nbr])); // [m] Mean grid radius
  const prc_cmp rds_min_min(rds_grd[0]); // [m] Minimum grid radius
  const prc_cmp rds_max_max(rds_grd[sz_nbr]); // [m] Maximum grid radius

  // Define diameter grid
  prc_cmp *dmt_ctr=new prc_cmp[sz_nbr]; // [m] Diameter at bin center
  prc_cmp *dmt_dlt=new prc_cmp[sz_nbr]; // [m] Width of diameter bin
  prc_cmp *dmt_grd=new prc_cmp[sz_nbr+1]; // [m] Diameter grid
  prc_cmp *dmt_max=new prc_cmp[sz_nbr]; // [m] Maximum diameter in bin
  prc_cmp *dmt_min=new prc_cmp[sz_nbr]; // [m] Minimum diameter in bin
  for(idx=0;idx<sz_nbr;idx++){
    dmt_ctr[idx]=2.0*rds_ctr[idx]; // [m] Diameter at bin center
    dmt_dlt[idx]=2.0*rds_dlt[idx]; // [m] Width of diameter bin
    dmt_grd[idx]=2.0*rds_grd[idx]; // [m] Diameter grid
    dmt_max[idx]=2.0*rds_max[idx]; // [m] Maximum diameter in bin
    dmt_min[idx]=2.0*rds_min[idx]; // [m] Minimum diameter in bin
  } // end loop over sz
  dmt_grd[sz_nbr]=2.0*rds_grd[sz_nbr]; // [m] Diameter grid
  const prc_cmp dmt_ctr_ctr(0.5*(dmt_grd[0]+dmt_grd[sz_nbr])); // [m] Mean grid diameter
  const prc_cmp dmt_min_min(dmt_grd[0]); // [m] Minimum grid diameter
  const prc_cmp dmt_max_max(dmt_grd[sz_nbr]); // [m] Maximum grid diameter

  /* Notes on proposed size distribution binning method:
     The following method is not yet implemented in mie(), but is the logical next step
     The method rests on the idea that the true or exact size distribution over which a function is to be integrated is known a priori
     An arbitrary function, e.g., fall speed or optical property, is to be discretized and stored to disk on a coarser resolution distribution
     We affix "_tsd" (for "true size distribution") to properties of the true distribution 
     We affix "_crs" (for "coarse") or nothing (since it is the default assumption) to properties of the coarse distribution 
     This true distribution is either computed (or read in from measurements) and kept in a first array
     Currently the true distribution is assumed to be a single lognormal (or gamma) distribution
     This is too inflexible but can be simply extended 
     The first step to handling more realistic distributions is to allow the true distribution to be constructed as a multi-modal superposition of distributions
     When transitioning to this new method many of the parameters currently labeled as "analytic" ("_anl" suffix) should perhaps be relabeled "_tsd"
     Since the true size distributions can be either multi-modal lognormal distributions or measurements, "_anl" does not really fit anymore
   */

  // Instantiate distribution functions
  // fxm: Make it unneccessary to carry all types of size distributions
  LogNormalFunction<prc_cmp> dst_lgn_aer;
  GammaFunction<prc_cmp> dst_gmm_aer;

  // Convolve size grid with particle distribution
  prc_cmp *dst=new prc_cmp[sz_nbr]; // [# m-3 m-1] Number distribution
  (void)vec_set(dst,sz_nbr,0.0);

  prc_cmp rds_nma(0.0); // [m] Number median radius analytic
  prc_cmp rds_nwa(0.0); // [m] Number weighted mean radius analytic
  prc_cmp rds_sma(0.0); // [m] Surface median radius analytic
  prc_cmp rds_swa(0.0); // [m] Surface area weighted mean radius analytic
  prc_cmp rds_vma(0.0); // [m] Volume median radius analytic
  prc_cmp rds_vwa(0.0); // [m] Volume weighted mean radius analytic
  // "Area" is too ambiguous to be useful in this program
  // Instead, we use xsa (cross-sectional area) or sfc (surface area)
  prc_cmp cnc_nbr_anl(0.0); // [# m-3] Number concentration analytic
  prc_cmp sfc_anl(0.0); // [m2 m-3] Surface area concentration analytic
  prc_cmp vlm_anl(0.0); // [m3 m-3] Volume concentration analytic
  prc_cmp xsa_anl(0.0); // [m2 m-3] Cross-sectional area concentration analytic
  for(psd_idx=0;psd_idx<psd_nbr;psd_idx++){
    dst_lgn_aer.prm_rst(psd_lst[psd_idx]); // Reset parameters of distribution function
    // fxm: Create better way to specify multi-modal gamma distributions
    dst_gmm_aer.rds_ffc_set(psd_lst[psd_idx].rds_nma_get()); 
    dst_gmm_aer.var_ffc_set(psd_lst[psd_idx].gsd_anl_get()); // 
    dst_gmm_aer.cnc_nbr_anl_set(psd_lst[psd_idx].cnc_nbr_anl_get()); // 
    if(psd_lst[psd_idx].typ_get() == "lognormal"){
      // Mass of aerosol per unit volume of air equals aerosol mass density times aerosol volume per unit volume of air
      cnc_nbr_anl+=dst_lgn_aer.cnc_nbr_anl_get(); // [# m-3] Number concentration analytic
      xsa_anl+=dst_lgn_aer.xsa_get(); // [m2 m-3] Cross-sectional area concentration analytic
      sfc_anl+=dst_lgn_aer.sfc_get(); // [m2 m-3] Surface area concentration analytic
      vlm_anl+=dst_lgn_aer.vlm_get(); // [m3 m-3] Volume concentration analytic

      // fxm: Most analytic weighted radii are still ill-defined
      rds_nwa+=dst_lgn_aer.cnc_nbr_anl_get()*dst_lgn_aer.rds_nwa_get(); // [m] Number weighted mean radius analytic
      rds_swa+=dst_lgn_aer.sfc_get()*dst_lgn_aer.rds_swa_get(); // [m] Surface area weighted mean radius analytic
      rds_vwa+=dst_lgn_aer.vlm_get()*dst_lgn_aer.rds_vwa_get(); // [m] Volume weighted mean radius analytic
      for(idx=0;idx<sz_nbr;idx++) dst[idx]+=dst_lgn_aer(sz_ctr[idx]); // [# m-3 m-1]
    }else if(psd_lst[psd_idx].typ_get() == "gamma"){
      cnc_nbr_anl+=dst_gmm_aer.cnc_nbr_anl_get(); // [# m-3] Number concentration analytic
      xsa_anl+=dst_gmm_aer.xsa_get(); // [m2 m-3] Cross-sectional area concentration analytic
      sfc_anl+=dst_gmm_aer.sfc_get(); // [m2 m-3] Surface area concentration analytic
      vlm_anl+=dst_gmm_aer.vlm_get(); // [m3 m-3] Volume concentration analytic
      for(idx=0;idx<sz_nbr;idx++) dst[idx]+=dst_gmm_aer(sz_ctr[idx]); // [# m-3 m-1]
    }else{
      err_prn(prg_nm,sbr_nm,psd_lst[psd_idx].typ_get()+" is not a valid psd_cls::Dst_typ");
    } // end if
  } // end loop over psd
  // Normalize by appropriate denominators
  rds_nwa/=cnc_nbr_anl; // [m] Number weighted mean radius analytic
  rds_swa/=sfc_anl; // [m] Surface area weighted mean radius analytic
  rds_vwa/=vlm_anl; // [m] Volume weighted mean radius analytic

  // For single mode distributions, analytic sizes are already known
  // For multimode distributions, many "analytic" sizes must be numerically interpolated
  if(psd_nbr == 1){
    rds_nma=dst_lgn_aer.rds_nma_get(); // [m] Number median radius analytic
    rds_sma=dst_lgn_aer.rds_sma_get(); // [m] Surface median radius analytic
    rds_vma=dst_lgn_aer.rds_vma_get(); // [m] Volume median radius analytic
  }else{
    // For a single mode distribution, analytic sizes are already known
    prc_cmp *nbr_prt_anl=new prc_cmp[sz_nbr]; // [# m-3] Number concentration of smaller particles analytic
    prc_cmp *sfc_prt_anl=new prc_cmp[sz_nbr]; // [m2 m-3] Surface area concentration of smaller particles analytic
    prc_cmp *vlm_prt_anl=new prc_cmp[sz_nbr]; // [m3 m-3] Volume concentration of smaller particles analytic
    prc_cmp erf_arg; // [frc] Argument to error function
    gsl_sf_result gsl_rsl; // [frc] GSL result structure
    nbr_prt_anl[0]=0.0; // [# m-3] Number concentration of smaller particles analytic
    sfc_prt_anl[0]=0.0; // [m2 m-3] Surface area concentration of smaller particles analytic
    vlm_prt_anl[0]=0.0; // [m3 m-3] Volume concentration of smaller particles analytic
    for(psd_idx=0;psd_idx<psd_nbr;psd_idx++){
      for(idx=0;idx<sz_nbr;idx++){
	erf_arg=std::log(rds_ctr[idx]/psd_rds_nma[psd_idx])/(std::sqrt(2.0)*std::log(psd_gsd_anl[psd_idx])); // [frc] Argument to error function
	rcd+=gsl_sf_erf_e(erf_arg,&gsl_rsl); // [frc] Error function
	nbr_prt_anl[idx]+=0.5*psd_cnc_nbr_anl[psd_idx]*(1.0+gsl_rsl.val); // [# m-3] Number concentration of smaller particles analytic
	// fxm: nbr_prt_anl is correct but sfc_prt_anl and vlm_prt_anl are bogus, must be derived analytically
	sfc_prt_anl[idx]+=0.5*psd_cnc_nbr_anl[psd_idx]*(1.0+gsl_rsl.val); // [m2 m-3] Surface area concentration of smaller particles analytic
	vlm_prt_anl[idx]+=0.5*psd_cnc_nbr_anl[psd_idx]*(1.0+gsl_rsl.val); // [m3 m-3] Volume concentration of smaller particles analytic
      } // end loop over size
    } // end loop over size distribution
    // Interpolate to find median size of total analytic distribution 
    // fxm: these ntp_vec_one() calls cause non-monotonic warnings because first few elements of *_prt_anl can be identically 0.0. Hmm. Nothing to worry about?
    wrn_prn(prg_nm,sbr_nm,"Possibly three non-monotonic warnings will follow, triggered by ntp_vec_one() calls to search for rds_*ma in multimodal distributions. These messages may safely be ignored.");
    if(sz_nbr > 1) rds_nma=ntp_vec_one(sz_nbr,nbr_prt_anl,rds_ctr,cnc_nbr_anl/2.0); else rds_nma=rds_ctr[0]; // [m] Number median radius analytic
    if(sz_nbr > 1) rds_sma=ntp_vec_one(sz_nbr,sfc_prt_anl,rds_ctr,sfc_anl/2.0); else rds_sma=rds_ctr[0]; // [m] Surface median radius analytic
    if(sz_nbr > 1) rds_vma=ntp_vec_one(sz_nbr,vlm_prt_anl,rds_ctr,vlm_anl/2.0); else rds_vma=rds_ctr[0]; // [m] Volume median radius analytic
    delete []nbr_prt_anl; // [# m-3] Number concentration of smaller particles analytic
    delete []sfc_prt_anl; // [m2 m-3] Surface area concentration of smaller particles analytic
    delete []vlm_prt_anl; // [m3 m-3] Volume concentration of smaller particles analytic
  } // endif psd_nbr > 1

  // Mass intensive properties
  const prc_cmp mss_anl(dns_prt*vlm_anl); // [kg m-3] Mass concentration analytic
  const prc_cmp cnc_spc_anl(cnc_nbr_anl/mss_anl); // [# kg-1] Specific concentration analytic
  const prc_cmp xsa_spc_anl(xsa_anl/mss_anl); // [m2 kg-1] Specific cross-sectional area analytic
  const prc_cmp sfc_spc_anl(sfc_anl/mss_anl); // [m2 kg-1] Specific Surface area analytic
  const prc_cmp vlm_spc_anl(vlm_anl/mss_anl); // [m3 kg-1] Specific analytic volume

  // Derive diameters from radii
  const prc_cmp dmt_nma(2.0*rds_nma); // [m] Number median diameter analytic
  const prc_cmp dmt_nwa(2.0*rds_nwa); // [m] Number weighted mean diameter analytic
  const prc_cmp dmt_sma(2.0*rds_sma); // [m] Surface median diameter analytic
  const prc_cmp dmt_swa(2.0*rds_swa); // [m] Surface area weighted mean diameter analytic
  const prc_cmp dmt_vma(2.0*rds_vma); // [m] Volume median diameter analytic
  const prc_cmp dmt_vwa(2.0*rds_vwa); // [m] Volume weighted mean diameter analytic

  prc_cmp *dst_mss=new prc_cmp[sz_nbr]; // [kg m-3 m-1] Mass distribution
  prc_cmp *dst_rds=new prc_cmp[sz_nbr]; // [m m-3 m-1] Radius distribution
  prc_cmp *dst_sfc=new prc_cmp[sz_nbr]; // [m2 m-3 m-1] Surface area distribution
  prc_cmp *dst_vlm=new prc_cmp[sz_nbr]; // [m3 m-3 m-1] Volume distribution
  prc_cmp *dst_xsa=new prc_cmp[sz_nbr]; // [m2 m-3 m-1] Cross-sectional area distribution
  prc_cmp *cnc=new prc_cmp[sz_nbr]; // [# m-3] Number concentration
  prc_cmp *mss=new prc_cmp[sz_nbr]; // [kg] Mass
  prc_cmp *sfc=new prc_cmp[sz_nbr]; // [m2] Surface area
  prc_cmp *vlm=new prc_cmp[sz_nbr]; // [m3] Volume
  prc_cmp *xsa=new prc_cmp[sz_nbr]; // [m2] Cross-sectional area

  prc_cmp *mss_prt_rsl=new prc_cmp[sz_nbr]; // [kg m-3] Mass concentration of smaller particles
  prc_cmp *mss_prt_rsl_frc=new prc_cmp[sz_nbr]; // [frc] Fraction of mass concentration from smaller particles
  prc_cmp *nbr_prt_rsl=new prc_cmp[sz_nbr];  // [# m-3] Number concentration of smaller particles
  prc_cmp *nbr_prt_rsl_frc=new prc_cmp[sz_nbr]; // [frc] Fraction of number concentration from smaller particles
  prc_cmp *sfc_prt_rsl=new prc_cmp[sz_nbr]; // [m2 m-3] Surface area concentration of smaller particles
  prc_cmp *sfc_prt_rsl_frc=new prc_cmp[sz_nbr]; // [frc] Fraction of surface area concentration from smaller particles
  prc_cmp *vlm_prt_rsl=new prc_cmp[sz_nbr]; // [m3 m-3] Volume concentration of smaller particles
  prc_cmp *vlm_prt_rsl_frc=new prc_cmp[sz_nbr]; // [frc] Fraction of volume concentration from smaller particles
  prc_cmp *xsa_prt_rsl=new prc_cmp[sz_nbr]; // [m2 m-3] Cross-sectional area concentration of smaller particles
  prc_cmp *xsa_prt_rsl_frc=new prc_cmp[sz_nbr]; // [frc] Fraction of cross-sectional area concentration from smaller particles
  prc_cmp cnc_nbr_rsl(0.0); // [# m-3] Number concentration resolved
  prc_cmp mss_rsl(0.0); // [kg m-3] Mass concentration resolved
  prc_cmp rds_nwr(0.0); // [m] Number weighted mean radius
  prc_cmp rds_swr(0.0); // [m] Surface area weighted mean radius
  prc_cmp rds_vwr(0.0); // [m] Mass weighted mean radius
  prc_cmp sfc_rsl(0.0); // [m2 m-3] Surface area concentration resolved
  prc_cmp vlm_rsl(0.0); // [m3 m-3] Volume concentration resolved
  prc_cmp xsa_rsl(0.0); // [m2 m-3] Cross-sectional area concentration resolved
  for(idx=0;idx<sz_nbr;idx++){
    // No need for loop over size modes here since dst[] includes all modes already
    cnc[idx]=dst[idx]*sz_dlt[idx]; // [# m-3] Number concentration
    cnc_nbr_rsl+=cnc[idx]; // [# m-3] Number concentration resolved
    nbr_prt_rsl[idx]=cnc_nbr_rsl; // [# m-3] Number concentration of smaller particles

    xsa[idx]=mth::cst_M_PIl*rds_ctr[idx]*rds_ctr[idx]; // [m2] Cross-sectional area of sphere of given size
    xsa_rsl+=xsa[idx]*cnc[idx]; // [m2 m-3] Cross-sectional area concentration resolved
    xsa_prt_rsl[idx]=xsa_rsl; // [m2 m-3] Cross-sectional area concentration of smaller particles

    sfc[idx]=4.0*xsa[idx]; // [m2] Surface area of sphere of given size
    sfc_rsl+=sfc[idx]*cnc[idx]; // [m2 m-3] Surface area concentration resolved
    sfc_prt_rsl[idx]=sfc_rsl; // [m2 m-3] Surface area concentration of smaller particles

    vlm[idx]=(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_ctr[idx],PRC_CMP(3.0)); // [m3] Volume of sphere of given size
    vlm_rsl+=vlm[idx]*cnc[idx]; // [m3 m-3] Volume concentration resolved
    vlm_prt_rsl[idx]=vlm_rsl; // [m3 m-3] Volume concentration of smaller particles

    mss[idx]=vlm[idx]*dns_prt; // [kg] Mass of sphere of given size
    mss_rsl+=mss[idx]*cnc[idx]; // [kg m-3] Mass concentration resolved
    mss_prt_rsl[idx]=mss_rsl; // [kg m-3] Mass concentration of smaller particles

    dst_rds[idx]=dst[idx]*sz_ctr[idx]; // [m m-3 m-1] Radius distribution
    dst_xsa[idx]=dst[idx]*xsa[idx]; // [m2 m-3 m-1] Cross-sectional area distribution
    dst_sfc[idx]=dst[idx]*sfc[idx]; // [m2 m-3 m-1] Surface area distribution
    dst_vlm[idx]=dst[idx]*vlm[idx]; // [m3 m-3 m-1] Volume distribution
    dst_mss[idx]=dst[idx]*mss[idx]; // [kg m-3 m-1] Mass distribution
    rds_nwr+=rds_ctr[idx]*cnc[idx]; // Arithmetic mean radius resolved
    rds_swr+=rds_ctr[idx]*sfc[idx]*cnc[idx]; // Surface area weighted radius resolved
    rds_vwr+=rds_ctr[idx]*vlm[idx]*cnc[idx]; // Mass weighted mean radius resolved
  } // end loop over sz

  // Normalize weighted integrals by resolved concentrations
  rds_nwr/=cnc_nbr_rsl; // [m] Arithmetic mean radius resolved
  rds_swr/=sfc_rsl; // [m] Surface area weighted mean radius resolved
  rds_vwr/=vlm_rsl; // [m] Mass weighted mean radius resolved

  // Mass-specific resolved quantities
  const prc_cmp cnc_spc_rsl(cnc_nbr_rsl/mss_rsl); // [# kg-1] Specific concentration resolved
  const prc_cmp xsa_spc_rsl(xsa_rsl/mss_rsl); // [m2 kg-1] Specific area resolved 
  const prc_cmp sfc_spc_rsl(sfc_rsl/mss_rsl); // [m2 kg-1] Specific Surface area resolved 
  const prc_cmp vlm_spc_rsl(vlm_rsl/mss_rsl); // [m3 kg-1] Specific volume resolved

  for(idx=0;idx<sz_nbr;idx++){
    nbr_prt_rsl_frc[idx]=nbr_prt_rsl[idx]/cnc_nbr_rsl; // [frc] Fraction of number concentration from smaller particles
    xsa_prt_rsl_frc[idx]=xsa_prt_rsl[idx]/xsa_rsl; // [frc] Fraction of cross-sectional area concentration from smaller particles
    sfc_prt_rsl_frc[idx]=sfc_prt_rsl[idx]/sfc_rsl; // [frc] Fraction of surface area concentration from smaller particles
    vlm_prt_rsl_frc[idx]=vlm_prt_rsl[idx]/vlm_rsl; // [frc] Fraction of volume concentration from smaller particles
    mss_prt_rsl_frc[idx]=mss_prt_rsl[idx]/mss_rsl; // [frc] Fraction of mass concentration from smaller particles
  } // end loop over sz

  // Find resolved median sizes by interpolating inverted size -> [number,surface,volume] relationships
  prc_cmp rds_nmr; // [m] Number median radius resolved
  prc_cmp rds_smr; // [m] Surface median radius resolved
  prc_cmp rds_vmr; // [m] Volume median radius resolved
  if(sz_nbr > 1) rds_nmr=ntp_vec_one(sz_nbr,nbr_prt_rsl_frc,rds_ctr,0.5); else rds_nmr=rds_ctr[0]; // [m] Number median radius resolved
  if(sz_nbr > 1) rds_smr=ntp_vec_one(sz_nbr,sfc_prt_rsl_frc,rds_ctr,0.5); else rds_smr=rds_ctr[0]; // [m] Surface median radius resolved
  if(sz_nbr > 1) rds_vmr=ntp_vec_one(sz_nbr,vlm_prt_rsl_frc,rds_ctr,0.5); else rds_vmr=rds_ctr[0]; // [m] Volume median radius resolved

  // Resolved fractions (i.e., discretization accuracy)
  const prc_cmp cnc_nbr_rsl_frc(cnc_nbr_rsl/cnc_nbr_anl); // [frc] Resolved fraction of number concentration
  const prc_cmp xsa_rsl_frc(xsa_rsl/xsa_anl); // [frc] Resolved fraction of cross-sectional area concentration
  const prc_cmp sfc_rsl_frc(sfc_rsl/sfc_anl); // [frc] Resolved fraction of surface area concentration
  const prc_cmp vlm_rsl_frc(vlm_rsl/vlm_anl); // [frc] Resolved fraction of volume concentration

  // Derive diameters from radii
  const prc_cmp dmt_nwr(2.0*rds_nwr); // [m] Number weighted mean diameter
  const prc_cmp dmt_swr(2.0*rds_swr); // [m] Surface area weighted mean diameter
  const prc_cmp dmt_vwr(2.0*rds_vwr); // [m] Mass weighted mean diameter
  const prc_cmp dmt_nmr(2.0*rds_nmr); // [m] Number median diameter resolved
  const prc_cmp dmt_smr(2.0*rds_smr); // [m] Surface median diameter resolved
  const prc_cmp dmt_vmr(2.0*rds_vmr); // [m] Volume median diameter resolved

  // Aspherical crystal effects for hexagonal prisms
  prc_cmp *asp_rat_hxg=new prc_cmp[sz_nbr]; // [frc] Hexagonal prism aspect ratio
  prc_cmp *cnc_vts=new prc_cmp[sz_nbr]; // [# m-3] Number concentration of equal V/S spheres
  prc_cmp *dmt_hxg=new prc_cmp[sz_nbr]; // [m] Length c of hexagonal prism
  prc_cmp *nbr_vts_per_hxg=new prc_cmp[sz_nbr]; // [nbr] Number equal V/S spheres per hexagonal prism
  prc_cmp *rds_ctr_vts=new prc_cmp[sz_nbr]; // [m] Radius at bin center of equal V/S spheres
  prc_cmp *rds_hxg=new prc_cmp[sz_nbr]; // [m] Half-width a of basal face of hexagonal prism
  prc_cmp *xsa_vts=new prc_cmp[sz_nbr]; // [m2] Equal-V/S sphere cross-sectional area
  prc_cmp *vlm_vts=new prc_cmp[sz_nbr]; // [m3] Equal-V/S sphere volume
  prc_cmp cnc_nbr_vts_rsl(0.0); // [frc] Equal-V/S sphere number concentration resolved
  prc_cmp hxg_fct_tmp; // [frc] Factor in computation of nbr_vts_per_hxg
  prc_cmp hxg_fct_tmp2; // [frc] Factor in computation of nbr_vts_per_hxg
  prc_cmp xsa_vts_rsl(0.0); // [m2 m-3] Equal-V/S sphere cross-sectional area concentration resolved
  prc_cmp mss_vts_rsl(0.0); // [kg m-3] Equal-V/S sphere mass concentration resolved
  for(idx=0;idx<sz_nbr;idx++){
    /* Conceptual and algorithmic framework for treating hexagonal prisms:
       Adopt all symbols and nomenclature from GrW99/NGW03
       Two general approaches are possible:
       1. Assume size distribution is prism half-width a, rds_ctr=rds_hxg 
       Combine rds_ctr with asp_rat_hxg to determine dmt_hxg
       Pros: Methods for specifying rds_[nma,swa,vma] may be re-used, although
             they now specify lognormal (usually) properties of the corresponding
	     aspherical dimension, e.g., hexagonal prism width
             Output of mie_slv() routine is weighted by cnc_vts rather than cnc
       Cons: Must evolve much code away from spherical assumption
             xsa, vlm arrays are for spheres, not hexagons
             Size argument to mie_slv() routine is not rds_ctr
             Breaks paradigm that dmt,rds,sfc,xsa,vlm are consistent
       2. Assume size distribution is spherical rds_hxg=rds_ctr*fxm
       Pros: dmt,sfc,xsa,vlm arrays are also spherical, hence self-consistent
       Cons: Size argument to mie_slv() routine is fxm 
       
       Approach #1 is the most intuitive and extensible approach because
       the simplest nomenclature (e.g., rds[idx]) specifies the real particles.
       Adopting approach #1 means eventually dropping the spherical assumption */ 

    // Assume aspect ratio is size invariant for now
    asp_rat_hxg[idx]=asp_rat_hxg_dfl; // [frc] Hexagonal prism aspect ratio NGW03 p. 3 (5)
    rds_hxg[idx]=rds_ctr[idx]; // [m] Half-width a of basal face of hexagonal prism NGW03 p. 3 (5)
    dmt_hxg[idx]=asp_rat_hxg[idx]*2.0*rds_hxg[idx]; // [m] Length c of hexagonal prism NGW03 p. 3 (5)
    hxg_fct_tmp=PRC_CMP(4.0)*asp_rat_hxg[idx]+std::sqrt(PRC_CMP(3.0)); // [frc] Factor in computation of nbr_vts_per_hxg NGW03 p. 3 (6b)
    hxg_fct_tmp2=std::pow(hxg_fct_tmp,PRC_CMP(3.0)); // [frc] Factor in computation of nbr_vts_per_hxg NGW03 p. 3 (6b)
    // rds_ctr_vts is radius at bin center of spheres with same V/S ratio as hexagonal prism of size (a,Gamma)=(rds_ctr,asp_rat_hxg)=(rds_hxg,asp_rat_hxg)
    rds_ctr_vts[idx]=3.0*std::sqrt(PRC_CMP(3.0))*rds_hxg[idx]*asp_rat_hxg[idx]/hxg_fct_tmp; // [m] Radius at bin center of equal V/S spheres NGW03 p.3 (6a)
    nbr_vts_per_hxg[idx]=hxg_fct_tmp2/(PRC_CMP(36.0)*mth::cst_M_PIl*asp_rat_hxg[idx]*asp_rat_hxg[idx]); // [nbr] Number equal V/S spheres per hexagonal prism NGW03 p. 3 (6b)
    // cnc_vts is number of equal-V/S-spheres with same volume (and area) as cnc hexagonal prism of size (a,Gamma)=(rds_ctr,asp_rat_hxg)=(rds_hxg,asp_rat_hxg)
    cnc_vts[idx]=nbr_vts_per_hxg[idx]*cnc[idx]; // [# m-3] Number concentration of equal V/S spheres
    cnc_nbr_vts_rsl+=cnc_vts[idx]; // [frc] Equal-V/S sphere number concentration resolved
    xsa_vts[idx]=mth::cst_M_PIl*rds_ctr_vts[idx]*rds_ctr_vts[idx]; // [m2] Equal-V/S sphere cross-sectional area
    xsa_vts_rsl+=xsa_vts[idx]*cnc_vts[idx]; // [m2 m-3] Equal-V/S sphere cross-sectional area concentration resolved
    vlm_vts[idx]=(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_ctr_vts[idx],PRC_CMP(3.0)); // [m3] Equal-V/S sphere volume
    mss_vts_rsl+=vlm_vts[idx]*cnc_vts[idx]*dns_prt; // [kg m-3] Equal-V/S sphere mass concentration resolved
  } // end loop over sz
  /* Analytic expressions exist for total number concentration, surface area, and volume of hexagonal prisms
     These expressions are similar to statistics of spheres except for pre-factors
     However, number concentration of equivalent V/S-spheres is complex function of 
     particle shape (see nbr_vts_per_hxg above).
     Analytic number concentration of V/S-spheres is integral of 
     nbr_vts_per_hxg times lognormal distribution function for hexagonal prisms
     Total concentration of equivalent V/S-spheres corresponding to
     known distribution of hexagonal prisms must be computed numerically unless 
     size dependence of aspect ratio is analytic.
     Simplest case is distribution of hexagonal prisms with constant aspect ratio
     Then the ratio \cncvts/\cnchxg is constant throughout the distribution
     Ergo the analytic number concentration of equivalent V/S-spheres
     is simply \cncvts/\cnchxg times the analytic number concentration of 
     hexagonal prisms which is presumably known directly from the lognormal
     size distribution parameters.

     Alternative is to compute cnc_nbr_vts_anl consistent with cnc_nbr_vts_rsl, e.g., with (temporary) wider size distribution at higher size resolution
     However this alternative method would set extremely undesirable precedent of looking outside user-specified size boundaries 
     For now we just assume the analytic and resolved are equal... */
  prc_cmp cnc_nbr_vts_anl=cnc_nbr_vts_rsl; // [# m-3] Equal-V/S sphere number concentration analytic
  cnc_nbr_vts_anl+=0.0; // CEWI [# m-3] Equal-V/S sphere number concentration analytic
  // For now, assume same as resolved fraction of spheres for first order correction
  const prc_cmp cnc_nbr_vts_rsl_frc(cnc_nbr_rsl_frc); // [frc] Equal-V/S sphere resolved fraction of number concentration

  // Equal-V/S optical properties
  if(vts_flg && !hxg_flg) err_prn(prg_nm,sbr_nm,"Attempted aspherical optical properties currently requires hexagonal prism assumption");
  prc_cmp *cnc_sph; // [# m-3] Equivalent sphere number concentration
  prc_cmp *sz_ctr_sph; // [m] Equivalent sphere size at bin center
  prc_cmp *xsa_sph; // [m2] Equivalent sphere cross-sectional area
  prc_cmp *vlm_sph; // [m3] Equivalent sphere volume
  prc_cmp mss_sph_rsl; // [kg m-3] Equivalent sphere mass concentration resolved
  prc_cmp vlm_sph_rsl; // [m3 m-3] Equivalent sphere volume concentration resolved
  prc_cmp xsa_sph_rsl; // [m2 m-3] Equivalent sphere cross-sectional area concentration resolved
  prc_cmp cnc_nbr_sph_rsl_frc; // [frc] Equivalent sphere resolved fraction of number concentration
  if(hxg_flg && vts_flg){
    cnc_nbr_sph_rsl_frc=cnc_nbr_vts_rsl_frc; // [frc] Equivalent sphere resolved fraction of number concentration
    cnc_sph=cnc_vts; // [# m-3] Equivalent sphere number concentration
    mss_sph_rsl=mss_vts_rsl; // [kg m-3] Equivalent sphere mass concentration resolved
    vlm_sph_rsl=mss_sph_rsl/dns_prt; // [m3 m-3] Equivalent sphere volume concentration resolved
    sz_ctr_sph=rds_ctr_vts; // [m] Equivalent sphere size at bin center
    xsa_sph=xsa_vts; // [m2] Equivalent sphere cross-sectional area
    xsa_sph_rsl=xsa_vts_rsl; // [m2 m-3] Equivalent sphere cross-sectional area concentration resolved
    vlm_sph=vlm_vts; // [m3] Equivalent sphere volume
  }else{ // not hexagonal prisms
    cnc_nbr_sph_rsl_frc=cnc_nbr_rsl_frc; // [frc] Equivalent sphere resolved fraction of number concentration
    cnc_sph=cnc; // [# m-3] Equivalent sphere number concentration
    mss_sph_rsl=mss_rsl; // [kg m-3] Equivalent sphere mass concentration resolved
    vlm_sph_rsl=vlm_rsl; // [m3 m-3] Equivalent sphere volume concentration resolved
    sz_ctr_sph=rds_ctr; // [m] Equivalent sphere size at bin center
    xsa_sph=xsa; // [m2] Equivalent sphere cross-sectional area
    xsa_sph_rsl=xsa_rsl; // [m2 m-3] Equivalent sphere cross-sectional area concentration resolved
    vlm_sph=vlm; // [m3] Equivalent sphere volume
  } // not equal-V/S optical properties

  // Determine RH and vapor properties
  prc_cmp RH_ice; // [frc] Relative humidity w/r/t ice water
  prc_cmp qst_H2O_lqd; // [kg kg-1] Saturation mixing ratio of H2O w/r/t liquid
  prc_cmp qst_H2O_ice; // [kg kg-1] Saturation mixing ratio of H2O w/r/t ice
  prc_cmp svp_H2O_lqd; // [Pa] Saturation vapor pressure of H2O w/r/t liquid
  prc_cmp svp_H2O_ice; // [Pa] Saturation vapor pressure of H2O w/r/t ice
  rcd+=svp_H2O_lqd_PrK78(1,&tpt_mdp,&svp_H2O_lqd);
  rcd+=svp_H2O_ice_PrK78(1,&tpt_mdp,&svp_H2O_ice);
  if(prs_mdp > one_mns_eps_H2O*svp_H2O_lqd){
    qst_H2O_lqd=eps_H2O*svp_H2O_lqd/(prs_mdp-one_mns_eps_H2O*svp_H2O_lqd); // [kg kg-1] Saturation mixing ratio of H2O w/r/t liquid
  }else{
    std::cerr << prg_nm << ": WARNING svp_H2O_lqd = " << svp_H2O_lqd/100.0 << " mb and (1 - eps_H2O)*svp_H2O_lqd = " << one_mns_eps_H2O*svp_H2O_lqd/100.0 << " mb at ambient pressure = " << prs_mdp/100.0 << " mb, and temperature = " << tpt_mdp-tpt_frz_pnt << " C\nSetting qst_H2O_lqd = RH_lqd = 0.0" << std::endl;
    qst_H2O_lqd=0.0; // [kg kg-1] Saturation mixing ratio of H2O w/r/t liquid
  } // endif 
  if(prs_mdp > one_mns_eps_H2O*svp_H2O_ice){
    qst_H2O_ice=eps_H2O*svp_H2O_ice/(prs_mdp-one_mns_eps_H2O*svp_H2O_ice); // [kg kg-1] Saturation mixing ratio of H2O w/r/t liquid
  }else{
    std::cerr << prg_nm << ": WARNING svp_H2O_ice = " << svp_H2O_ice/100.0 << " mb and (1 - eps_H2O)*svp_H2O_ice = " << one_mns_eps_H2O*svp_H2O_ice/100.0 << " mb at ambient pressure = " << prs_mdp/100.0 << " mb, and temperature = " << tpt_mdp-tpt_frz_pnt << " C\nSetting qst_H2O_ice = RH_ice = 0.0" << std::endl;
    qst_H2O_ice=0.0; // [kg kg-1] Saturation mixing ratio of H2O w/r/t liquid
  } // endif 
  if(q_H2O_vpr != cmd_ln_dfl){
    // Derive RH_lqd from specified q_H2O_vpr, T, p
    RH_lqd=q_H2O_vpr/qst_H2O_lqd; // [kg kg-1] Specific humidity
  }else{
    // Convert specified RH_lqd into q_H2O_vpr for given T, p
    q_H2O_vpr=qst_H2O_lqd*RH_lqd; // [kg kg-1] Specific humidity
  } // endif
  RH_ice=q_H2O_vpr/qst_H2O_ice; // [frc] Relative humidity w/r/t ice water

  // Stokes' gravitational settling velocity
  prc_cmp mfp_atm; // [m] Mean free path of atmosphere
  using phc::eps_H2O_rcp_m1; // (0.60777) [frc] Constant for virtual temperature
  prc_cmp tpt_vrt=tpt_mdp*(1.0+eps_H2O_rcp_m1*q_H2O_vpr); // [K] Virtual temperature
  using phc::gas_cst_dry_air; // (287.05) [J kg-1 K-1] Gas constant of dry_air IrG81 p. 25, p. 245
  prc_cmp dns_mdp=prs_mdp/(gas_cst_dry_air*tpt_vrt); // [kg m-3] Midlayer density
  prc_cmp mmr_rsl=mss_rsl/dns_mdp; // [kg kg-1] Aerosol mixing ratio resolved
  prc_cmp *slp_crc=new prc_cmp[sz_nbr]; // [frc] Slip correction factor
  prc_cmp *vlc_stk=new prc_cmp[sz_nbr]; // [m s-1] Stokes' settling velocity (Re < 0.1)
  const prc_cmp vsc_dyn_atm(vsc_dyn_atm_fst_scl(tpt_mdp)); // [kg m-1 s-1] Dynamic viscosity of atmosphere
  using phc::mmw_dry_air; // (28.9644e-3) [kg mol-1] (Source: radcsw.F in CCM2/3)
  using phc::gas_cst_unv; // (8.31441) [J mol-1 K-1] Universal gas constant
  mfp_atm=2.0*vsc_dyn_atm/(prs_mdp*std::sqrt(8.0*mmw_dry_air/(mth::cst_M_PIl*gas_cst_unv*tpt_mdp))); // [m] Mean free path of atmosphere SeP97 p. 455
  rcd+=vlc_stk_get
    (sz_nbr, // I [nbr] Number of size bins
     dmt_ctr, // I [m] Diameter at bin center
     dns_prt, // I [kg m-3] Density of particle
     mfp_atm, // I [m] Mean free path of atmosphere
     vsc_dyn_atm, // I [kg m-1 s-1] Dynamic viscosity of atmosphere
     slp_crc, // O [frc] Slip correction factor
     vlc_stk); // O [m s-1] Stokes' settling velocity (Re < 0.1)

  // Aerosol gravitational settling velocity
  // Using Stokes' velocity rather than iterative solution with empirical drag coefficient causes 60% errors for D = 200 um SeP97 p. 468
  prc_cmp *cff_drg_grv=new prc_cmp[sz_nbr]; // [frc] Drag coefficient at terminal velocity
  prc_cmp *ryn_nbr_grv=new prc_cmp[sz_nbr]; // [frc] Reynolds number at terminal velocity
  prc_cmp *stk_crc=new prc_cmp[sz_nbr]; // [frc] Correction to Stokes settling velocity
  prc_cmp *vlc_grv=new prc_cmp[sz_nbr]; // [m s-1] Settling velocity
  const prc_cmp vsc_knm_atm(vsc_dyn_atm/dns_mdp); // [m2 s-1] Kinematic viscosity of atmosphere 
  prc_cmp *asp_rat_lps=new prc_cmp[sz_nbr]; // [frc] Ellipsoidal aspect ratio
  rcd+=asp_rat_lps_get // [fnc] Determine particle aspect ratio
    (sz_nbr, // I [nbr] Size of arrays
     cmp_sng_prt, // I [sng] Composition of particle
     asp_rat_lps_dfl, // I [frc] Ellipsoidal aspect ratio
     asp_rat_lps); // O [frc] Ellipsoidal aspect ratio
  if(asp_rat_lps_dfl == PRC_CMP(1.0)) 
    rcd+=vlc_grv_get
      (sz_nbr, // I [nbr] Size of arrays
       dmt_ctr, // I [m] Diameter at bin center
       slp_crc, // I [frc] Slip correction factor SeP97 p. 464
       vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
       dns_mdp, // I [kg m-3] Midlayer density
       dns_prt, // I [kg m-3] Density of particle
       vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
       cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
       ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
       stk_crc, // O [frc] Correction to Stokes settling velocity
       vlc_grv); // O [m s-1] Settling velocity
  else 
    rcd+=vlc_grv_get_Gin03
      (sz_nbr, // I [nbr] Size of arrays
       asp_rat_lps, // I [frc] Ellipsoidal aspect ratio
       dmt_ctr, // I [m] Diameter at bin center
       slp_crc, // I [frc] Slip correction factor SeP97 p. 464
       vlc_stk, // I [m s-1] Stokes' settling velocity (Re < 0.1)
       dns_mdp, // I [kg m-3] Midlayer density
       dns_prt, // I [kg m-3] Density of particle
       vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
       cff_drg_grv, // O [frc] Drag coefficient at terminal velocity
       ryn_nbr_grv, // O [frc] Reynolds number at terminal velocity
       stk_crc, // O [frc] Correction to Stokes settling velocity
       vlc_grv); // O [m s-1] Settling velocity
    
  // Weight by number and by mass
  prc_cmp vlc_grv_nwr(0.0); // [m s-1] Number weighted terminal velocity
  prc_cmp vlc_grv_vwr(0.0); // [m s-1] Mass weighted terminal velocity
  for(idx=0;idx<sz_nbr;idx++){
    vlc_grv_nwr+=vlc_grv[idx]*cnc[idx]; // [m s-1] Number weighted terminal velocity
    vlc_grv_vwr+=vlc_grv[idx]*vlm[idx]*cnc[idx]; // [m s-1] Mass weighted terminal velocity
  } // end loop over sz
  vlc_grv_nwr/=cnc_nbr_rsl; // [m s-1] Number weighted terminal velocity
  vlc_grv_vwr/=vlm_rsl; // [m s-1] Mass weighted terminal velocity

  // Diagnostic diameters
  prc_cmp *dmt_mjr=new prc_cmp[sz_nbr]; // [m] Major axis of spheroid
  prc_cmp *dmt_mnr=new prc_cmp[sz_nbr]; // [m] Minor axis of spheroid
  prc_cmp *dmt_stk=new prc_cmp[sz_nbr]; // [m] Stokes diameter
  prc_cmp *dmt_aer=new prc_cmp[sz_nbr]; // [m] Aerodynamic diameter
  prc_cmp *dmt_eqv_sfc=new prc_cmp[sz_nbr]; // [m] Diameter of sphere with same surface area
  prc_cmp *dmt_eqv_vlm=new prc_cmp[sz_nbr]; // [m] Diameter of sphere with same volume
  for(idx=0;idx<sz_nbr;idx++){
    /* Particles are ellipsoidal when asp_rat_lps != 1.0
       In this case, diameter does not suffice to characterize particle dimension
       First, dmt_ctr remains the canonical particle size and so must be assigned a physically meaningful definition
       Options are to set dmt_ctr to one of 
       dmt_stk: Diameter of sphere with same terminal settling velocity and density
       dmt_aer: Diameter of sphere with same terminal settling velocity but unit density
       dmt_eqv_sfc: Diameter of sphere with same surface area
       dmt_eqv_vlm: Diameter of sphere with same surface volume
       dmt_mjr: Major axis of ellipsoid is diameter along "a" axis, i.e., 2*a
       dmt_mnr: Minor axis of ellipsoid is diameter along "b" axis, i.e., 2*b
       We set dmt_ctr equal to diameter of sphere with same surface area as in Gin03 p. 2 (10)
       This is sometimes called the surface equivalent diameter */
    if(asp_rat_lps[idx] == 1.0){
      dmt_mjr[idx]=dmt_ctr[idx]; // [m] Major axis of ellipsoid
      dmt_mnr[idx]=dmt_ctr[idx]; // [m] Minor axis of ellipsoid
    }else{
      // Minor axis of ellipsoid is diameter along "b" axis, i.e., 2*b
      dmt_mnr[idx]=2.0*dmt_ctr[idx]/psi_lps_fst_scl(asp_rat_lps[idx]); // [m] Minor axis of ellipsoid Gin03 p.2 (10)
      // Major axis of ellipsoid is diameter along "a" axis, i.e., 2*a
      dmt_mjr[idx]=dmt_mnr[idx]*asp_rat_lps[idx]; // [m] Major axis of ellipsoid
    } // endif ellipsoid
    dmt_eqv_sfc[idx]=dmt_eqv_sfc_lps_fst_scl(0.5*dmt_mjr[idx],0.5*dmt_mnr[idx]); // [m] Diameter of sphere with same surface area
    dmt_eqv_vlm[idx]=dmt_eqv_vlm_lps_fst_scl(0.5*dmt_mjr[idx],0.5*dmt_mnr[idx]); // [m] Diameter of sphere with same volume
  } // end loop over sz

  /* Stokes diameter is diameter of sphere of same density with same terminal settling velocity
     Stokes diameter equals particle diameter when particle is spherical */
  rcd+=dmt_stk_get // [fnc] Compute Stokes diameter of particles
    (sz_nbr, // I [nbr] Size of arrays
     dmt_ctr, // I [m] Diameter at bin center
     dns_mdp, // I [kg m-3] Midlayer density
     dns_prt, // I [kg m-3] Particle density
     mfp_atm, // I [m] Mean free path of atmosphere
     vlc_grv, // I [m s-1] Settling velocity
     vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
     dmt_stk); // O [m] Stokes diameter

  // Aerodynamic diameter is diameter of unit density sphere with same terminal velocity
  rcd+=dmt_aer_get // [fnc] Compute aerodynamic diameter of particles
    (sz_nbr, // I [nbr] Size of arrays
     dmt_stk, // I [m] Stokes diameter
     slp_crc, // I [frc] Slip correction factor
     dns_prt, // I [kg m-3] Particle density
     mfp_atm, // I [m] Mean free path of atmosphere
     dmt_aer); // O [m] Aerodynamic diameter

  /* 19990731: Imposing LSM-like defaults for non-land causes problems since
     many land routines are automatically called anyway, with, e.g., vwc_sfc=1.0e30
  // Make land surface properties consistent with orography
  if(oro_is_ocn(oro)){ // Ocean
    lnd_frc_dry=mss_frc_cly=mss_frc_snd=0.0; // [frc] 
    vwc_sfc=1.0e30; // [m3 m-3] Volumetric water content
  }else if(oro_is_ice(oro)){ // Sea ice
    lnd_frc_dry=mss_frc_cly=mss_frc_snd=0.0; // [frc] 
    vwc_sfc=1.0; // [m3 m-3] Volumetric water content
  } // end if non-land
  */
  prc_cmp tpt_sfc; // [K] Surface temperature
  if(oro_is_lnd(oro))
    tpt_sfc=tpt_gnd; // [K] Surface temperature
  else if(oro_is_ocn(oro)) 
    tpt_sfc=tpt_sst; // [K] Surface temperature
  else if(oro_is_ice(oro)) 
    tpt_sfc=tpt_ice; // [K] Surface temperature
  else 
    err_prn(prg_nm,sbr_nm,"Unknown orography type");

  /* Default reference height at 10 m works for atmospheric aerosol, not for snow
     For snow simulations we often want layer thicknesses much less than 10 m
     Assume that mid-layer thicknesses < 10 m indicates primary interest is snow 
     In that case, reset reference height to small value */
  prc_cmp snw_hgt; // [m] Geometric bulk thickness of snow
  if(hgt_rfr > hgt_mdp){
    std::cout << prg_nm_get() << ": INFO " << sbr_nm << ": Re-setting reference height from " << hgt_rfr << " m to " << hgt_mdp*2.0 << " m to allow realistic surface aerosol (snow) layer thicknesses. Diagnostics which relied on 10 m reference height (like interpolated winds, Owens' effect for dust) may not be trustworthy." << std::endl;
    hgt_rfr=hgt_mdp*2.0; // [m] Reference height
    using phc::dns_H2O_snw_gnd_std; // (100.0) [kg m-3] Standard bulk density of snow on ground WiW80 p. 2724, 2725, CCM:physics/tsinti()
    using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
    std::cout << prg_nm_get() << ": INFO " << sbr_nm << ": Setting snow height to " << hgt_mdp*2.0 << " m for consistency" << std::endl;
    snw_hgt=hgt_mdp*2.0; // [m] Geometric bulk thickness of snow
    std::cout << prg_nm_get() << ": INFO " << sbr_nm << ": Re-setting liquid water equivalent snow depth (snw_hgt_lqd) from " << snw_hgt_lqd;
    snw_hgt_lqd=snw_hgt*dns_H2O_snw_gnd_std/dns_H2O_lqd_std; // [m] Equivalent liquid water snow depth
    std::cout << " m to " << snw_hgt_lqd << " m for consistency" << std::endl;
  } // endif hgt_rfr > hgt_mdp

  // Fraction of surface covered by snow
  prc_cmp snw_frc; // [frc] Fraction of surface covered by snow
  rcd+=snw_frc_get
    (lon_nbr, // I [nbr] Size of arrays
     &snw_frc, // O [frc] Fraction of surface covered by snow
     &snw_hgt, // O [m] Geometric bulk thickness of snow
     &snw_hgt_lqd); // I [m] Equivalent liquid water snow depth

  // Begin mobilization routines
  prc_cmp wnd_mdp; // [m s-1] Surface layer mean wind speed
  wnd_mdp=std::sqrt(wnd_znl_mdp*wnd_znl_mdp+wnd_mrd_mdp*wnd_mrd_mdp); // [m s-1] Surface layer mean wind speed
  prc_cmp wnd_mdp_bnd; // [m s-1] Bounded surface layer mean wind speed
  wnd_mdp_bnd=max_cpv(wnd_mdp,wnd_min_mbl); // [m s-1] Bounded surface layer mean wind speed

  // Virtual temperature
  prc_cmp tpt_vrt_mdp; // [K] Virtual temperature
  rcd+=tpt_vrt_get
    (lon_nbr, // I [nbr] Size of arrays
     &q_H2O_vpr, // I [kg kg-1] Specific humidity
     &tpt_mdp, // I [K] Midlayer temperature
     &tpt_vrt_mdp); // O [K] Virtual temperature

  // Potential temperature
  prc_cmp tpt_ptn_mdp; // [K] Midlayer local potential temperature
  rcd+=tpt_ptn_rfr_get
    (lon_nbr, // I [nbr] Size of arrays
     &prs_mdp, // I [Pa] Midlayer pressure
     &prs_ntf, // I [Pa] Reference pressure
     &tpt_mdp, // I [K] Temperature
     &tpt_ptn_mdp); // O [K] Potential temperature

  // Virtual potential temperature
  prc_cmp tpt_ptn_vrt_mdp; // [K] Virtual potential temperature
  rcd+=tpt_vrt_get
    (lon_nbr, // I [nbr] Size of arrays
     &q_H2O_vpr, // I [kg kg-1] Specific humidity
     &tpt_ptn_mdp, // I [K] Potential temperature
     &tpt_ptn_vrt_mdp); // O [K] Virtual potential temperature

  // Hydrologic Properties
  prc_cmp mss_frc_slt(1.0-mss_frc_snd-mss_frc_cly); // [frc] Mass fraction silt
  prc_cmp smp_sat; // [mm H2O] Saturated soil matric potential (sand-dependent)
  prc_cmp smp_sfc; // [mm H2O] Soil matric potential
  prc_cmp smp_xpn_b; // [frc] Exponent "b" for smp (clay-dependent)
  prc_cmp vwc_dry; // [m3 m-3] Dry volumetric water content (no E-T)
  prc_cmp vwc_opt; // [m3 m-3] E-T optimal volumetric water content 
  prc_cmp vwc_rel; // [frc] Water content relative to saturation
  prc_cmp vwc_sat; // [m3 m-3] Saturated volumetric water content (sand-dependent)
 
  smp_xpn_b= // [frc] Exponent "b" for smp (clay-dependent)
    2.91+0.159*mss_frc_cly*100.0; // Bon96 p. 98
  // NB: Adopt convention that matric potential is positive definite
  smp_sat= // [mm H2O] Saturated soil matric potential (sand-dependent)
    10.0*std::pow(PRC_CMP(10.0),PRC_CMP(1.88)-PRC_CMP(0.0131)*mss_frc_snd*PRC_CMP(100.0)); // Bon96 p. 98
  vwc_sat= // [m3 m-3] Saturated volumetric water content (sand-dependent)
    0.489-0.00126*mss_frc_snd*100.0; // Bon96 p. 98
  assert(mss_frc_slt >= 0.0 && mss_frc_slt <= 1.0);

  if(oro_is_lnd(oro)) assert(vwc_sfc <= vwc_sat);
  vwc_dry= // [m3 m-3] Dry volumetric water content (no E-T)
    vwc_sat*std::pow(PRC_CMP(316230.0)/smp_sat,PRC_CMP(-1.0)/smp_xpn_b); // Bon96 p. 98
  vwc_opt= // [m3 m-3] E-T optimal volumetric water content
    vwc_sat*std::pow(PRC_CMP(158490.0)/smp_sat,PRC_CMP(-1.0)/smp_xpn_b); // Bon96 p. 98
  vwc_rel= // [frc] Water content relative to saturation
    max_cpv(vwc_sfc/vwc_sat,0.05); // Bon96 p. 97
  // NB: Adopt convention that matric potential is positive definite
  smp_sfc= // [mm H2O] Soil matric potential
    smp_sat*std::pow(vwc_rel,-smp_xpn_b); // Bon96 p. 97

  // Determine gravimetric water content from volumetric water content and soil density
  using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
  using phc::dns_dst_std; // (2500.0) [kg m-3] Standard density of dust
  const prc_cmp dns_prt_sfc(dns_dst_std); // [kg m-3] Dry density of soil particles (excluding pores)
  const prc_cmp dns_blk_dry(dns_prt_sfc*(1.0-vwc_sat)); // [kg m-3] Bulk density of dry surface soil
  const prc_cmp dns_blk_sfc(dns_blk_dry+vwc_sfc*dns_H2O_lqd_std); // [kg m-3] Bulk density of surface soil
  const prc_cmp gwc_sfc(vwc_sfc*dns_H2O_lqd_std/dns_blk_dry); // [kg kg-1] Gravimetric water content
  
  // Transfer efficiency of vapor from soil to atmosphere
  prc_cmp trn_fsh_vpr_soi_atm; // [frc] Transfer efficiency of vapor from soil to atmosphere
  rcd+=trn_fsh_vpr_soi_atm_get
    (lon_nbr, // I [nbr] Size of arrays
     &tpt_soi, // I [K] Soil temperature
     tpt_frz_pnt, // I [K] Temperature of frozen soil
     &trn_fsh_vpr_soi_atm, // O [frc] Transfer efficiency of vapor from soil to atmosphere
     &vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
     &vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
     &vwc_sfc); // I [m3 m-3] Volumetric water content

  // Soil thermal conductivity
  prc_cmp lvl_dlt_snw; // [m] Soil layer thickness including snow
  prc_cmp cnd_trm_soi; // [W m-1 K-1] Soil thermal conductivity
  rcd+=cnd_trm_soi_get
    (lon_nbr, // I [nbr] Size of arrays
     &cnd_trm_soi, // O [W m-1 K-1] Soil thermal conductivity
     &lvl_dlt_snw, // O [m] Soil layer thickness including snow
     &mss_frc_cly, // I [frc] Mass fraction clay 
     &mss_frc_snd, // I [frc] Mass fraction sand
     &snw_hgt, // I [m] Geometric bulk thickness of snow
     &tpt_soi, // I [K] Soil temperature
     &vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
     &vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
     &vwc_sat, // I [m3 m-3] Saturated volumetric water content
     &vwc_sfc); // I [m3 m-3] Volumetric water content
  
  // Boundary layer exchange over erodible surfaces
  prc_cmp wnd_frc_mbl; // [m s-1] Surface friction velocity for erodible surface
  prc_cmp mno_lng_mbl(-15.0); // [m] Monin-Obukhov length for erodible surface
  prc_cmp tpt_gnd_mbl(tpt_gnd); // [K] Ground temperature for erodible surface
  using lsm::sfc_ems; // [frc] Bare ground emissivity
  if(msv_gnd == cmd_ln_dfl) msv_gnd=sfc_ems[soi_typ]; // [frc] Bare ground emissivity
  rcd+=blm_mbl
    (lon_nbr, // I [nbr] Size of arrays
     &cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
     &dns_mdp, // I [kg m-3] Midlayer density
     &flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
     &flx_SW_net_gnd, // I [W m-2] Solar flux absorbed by ground
     &hgt_mdp, // I [m] Midlayer height above surface
     &hgt_zpd_mbl, // I [m] Zero plane displacement for erodible surface
     &lvl_dlt_snw, // I [m] Soil layer thickness including snow
     &msv_gnd, // I [frc] Bare ground emissivity
     &prs_mdp, // I [Pa] Pressure
     &q_H2O_vpr, // I [kg kg-1] Specific humidity
     &rgh_mmn_mbl, // I [m] Roughness length momentum for erodible surface
     &tpt_mdp, // I [K] Midlayer temperature
     &tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
     &tpt_soi, // I [K] Soil temperature
     &trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
     &wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     &wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     &mno_lng_mbl, // I/O [m] Monin-Obukhov length for erodible surface
     &tpt_gnd_mbl, // I/O [K] Ground temperature for erodible surface
     &wnd_frc_mbl); // O [m s-1] Surface friction velocity for erodible surface
    
  // Integrate to wind speed at 10 m from midlayer wind speed
  prc_cmp rgh_heat(rgh_mmn_mbl); // [m] Roughness length heat
  prc_cmp rgh_vpr(rgh_mmn_mbl); // [m] Roughness length vapor
  prc_cmp hgt_ash(rgh_heat+hgt_zpd_mbl); // [m] Apparent sink for heat
  prc_cmp hgt_asm(rgh_mmn_mbl+hgt_zpd_mbl); // [m] Apparent sink for momentum
  prc_cmp hgt_asv(rgh_vpr+hgt_zpd_mbl); // [m] Apparent sink for vapor
  prc_cmp mmn_dlt; // [m s-1] Upper level - lower level windspeed
  prc_cmp wnd_rfr_mbl; // [m s-1] Wind speed at reference height for erodible surface
  rcd+=mmn_dlt_evl
    (lon_nbr, // I [nbr] Size of arrays
     &hgt_rfr, // I [m] Reference height for dust mobilization processes
     &hgt_mdp, // I [m] Midlayer height above surface
     &hgt_zpd_mbl, // I [m] Zero plane displacement
     &mmn_dlt, // O [m s-1] Upper level - lower level windspeed
     &mno_lng_mbl, // I [m] Monin-Obukhov length
     &wnd_frc_mbl); // I [m s-1] Surface friction velocity
  wnd_rfr_mbl=wnd_mdp-mmn_dlt; // [m s-1] Wind speed at reference height
  if(wnd_rfr_mbl < wnd_min_mbl){
    if(dbg_lvl >= dbg_sbr) std::cerr << prg_nm_get() << ": INFO " << sbr_nm << ": Reset U(10 m) for mobilization from " << wnd_rfr_mbl << " m s-1 to " << wnd_min_mbl << " m s-1" << std::endl;
    wnd_rfr_mbl=wnd_min_mbl;
  } // endif
  assert(wnd_rfr_mbl >= 0.0);
  prc_cmp wnd_mrd_rfr; // [m s-1] Reference level meridional wind speed
  prc_cmp wnd_znl_rfr; // [m s-1] Reference level zonal wind speed
  wnd_mrd_rfr=wnd_mrd_mdp*wnd_rfr_mbl/wnd_mdp_bnd; // [m s-1] Reference level meridional wind speed
  wnd_znl_rfr=wnd_znl_mdp*wnd_rfr_mbl/wnd_mdp_bnd; // [m s-1] Reference level zonal wind speed

  if(true){
    /* Integrating to wind speed at 10 m downward from midlayer should be
       equivalent to interpolating upward from apparent sink for momentum.
       There may be conditions where this ideal test fails, however
       Stability functions are only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
       For highly unstable conditions, -1.0 < L < 0.0 and zeta << -100.0
       This only occurs for very small mean wind speeds, i.e., wnd_mdp < 2.0 m s-1
       Then, integrating from either direction exposes an instability
       fxm: Under exactly what conditions does integrating top-->btm != btm-->top fail?
       I suspect it is only when |zeta| >> 10.0
    */
    prc_cmp wnd_rfr_dwn(wnd_rfr_mbl); // [m s-1] Wind speed obtained integrating downward
    prc_cmp wnd_rfr_up; // [m s-1] Wind speed obtained integrating upward

    // Integrate to wind speed at 10 m from apparent sink for momentum
    rcd+=mmn_dlt_evl
      (lon_nbr, // I [nbr] Size of arrays
       &hgt_asm, // I [m] Apparent sink for momentum
       &hgt_rfr, // I [m] Reference height for dust mobilization processes
       &hgt_zpd_mbl, // I [m] Zero plane displacement
       &mmn_dlt, // O [m s-1] Upper level - lower level windspeed
       &mno_lng_mbl, // I [m] Monin-Obukhov length
       &wnd_frc_mbl); // I [m s-1] Surface friction velocity
    wnd_rfr_up=mmn_dlt; // [m s-1] Wind speed at reference height
    assert(wnd_rfr_up >= 0.0);
    wnd_rfr_up=max_cpv(wnd_rfr_up,wnd_min_mbl); // [m s-1]

    // Relative convergence
    if(std::fabs((wnd_rfr_up-wnd_rfr_dwn)/wnd_rfr_dwn) > 1.0e-2) std::cerr << prg_nm << ": WARNING " << sbr_nm << ": Upward and downward integrations to reference wind speed disagree: " << wnd_rfr_up << " m s-1 != " << wnd_rfr_dwn << " m s-1" << std::endl;
  } // endif true

  // Drag partition to account for non-erodible roughness elements 
  prc_cmp frc_thr_ncr_drg(1.0); // [frc] Factor by which surface roughness increases threshold friction velocity
  prc_cmp wnd_frc_fsh_frc(1.0); // [frc] Efficient fraction of wind friction
  rcd+=frc_thr_ncr_drg_get
    (lon_nbr, // I [nbr] Size of arrays
     &frc_thr_ncr_drg, // O [frc] Factor by which surface roughness increases threshold friction velocity
     &wnd_frc_fsh_frc, // O [frc] Efficient fraction of wind friction
     rgh_mmn_mbl, // I [m] Roughness length momentum for erodible surfaces
     rgh_mmn_smt); // I [m] Smooth roughness length

  // Account for drag partition to determine friction speed over erodible surfaces
  prc_cmp wnd_frc_smt; // [m s-1] Friction speed over erodible surface
  wnd_frc_smt=wnd_frc_mbl*wnd_frc_fsh_frc; // [m s-1] Friction speed over erodible surface MaB95 p. 16420, GMB98 p. 6207

  // Factor by which soil moisture increases threshold friction velocity  
  prc_cmp frc_thr_ncr_wtr(1.0); // [frc] Factor by which moisture increases threshold friction velocity
  prc_cmp vwc_thr; // [m3 m-3] Threshold volumetric water content to affect mobilization
  rcd+=frc_thr_ncr_wtr_get
    (lon_nbr, // I [nbr] Size of arrays
     &frc_thr_ncr_wtr, // O [frc] Factor by which moisture increases threshold friction velocity
     &vwc_thr, // O [m3 m-3] Threshold volumetric water content to affect mobilization
     &mss_frc_cly, // I [frc] Mass fraction of clay
     &vwc_sfc); // I [frc] Volumetric water content % fxm: should be gwc_sfc

  // Threshold friction speed for saltation over dry, bare, smooth ground
  prc_cmp wnd_frc_thr_slt; // [m s-1] Threshold friction speed for saltation
  rcd+=wnd_frc_thr_slt_get
    (lon_nbr, // I [nbr] Size of arrays
     &dns_mdp, // I [kg m-3] Midlayer density
     &wnd_frc_thr_slt); // O [m s-1] Threshold friction speed for saltation

  // Adjust threshold friction velocity to account for moisture and roughness
  const prc_cmp wnd_frc_thr_slt_dry_flt(wnd_frc_thr_slt); // [m s-1] Threshold friction velocity for dry, smooth ground
  wnd_frc_thr_slt*=frc_thr_ncr_wtr*frc_thr_ncr_drg; // [m s-1] Threshold friction speed for saltation

  // Threshold windspeeds
  prc_cmp wnd_rfr_thr_slt; // [m s-1] Threshold 10 m wind speed for saltation
  prc_cmp wnd_mdp_thr_slt; // [m s-1] Threshold midlayer wind speed for saltation
  wnd_rfr_thr_slt=wnd_rfr_mbl*wnd_frc_thr_slt/wnd_frc_mbl; // [m s-1] Threshold 10 m wind speed for saltation
  wnd_mdp_thr_slt=wnd_mdp_bnd*wnd_frc_thr_slt/wnd_frc_mbl; // [m s-1] Threshold midlayer wind speed for saltation

  // Saltating friction velocity equals friction velocity if no saltation
  prc_cmp wnd_frc_slt; // [m s-1] Saltating friction velocity
  prc_cmp wnd_frc_slt_dlt; // [m s-1] Friction velocity increase from saltation
  rcd+=wnd_frc_slt_get
    (lon_nbr, // I [nbr] Size of arrays
     &wnd_frc_mbl, // I [m s-1] Surface friction velocity
     &wnd_frc_slt, // O [m s-1] Saltating friction velocity
     &wnd_frc_slt_dlt, // O [m s-1] Friction velocity increase from saltation
     &wnd_rfr_mbl, // I [m s-1] Wind speed at reference height
     &wnd_rfr_thr_slt); // I [m s-1] Threshold 10 m wind speed for saltation
    
  // Horizontal streamwise mass flux
  prc_cmp flx_mss_hrz_slt_ttl; // [kg m-1 s-1] Vertically integrated streamwise mass flux
  rcd+=flx_mss_hrz_slt_ttl_Whi79_get
    (lon_nbr, // I [nbr] Size of arrays
     &dns_mdp, // I [kg m-3] Midlayer density
     &flx_mss_hrz_slt_ttl, // O [kg m-1 s-1] Vertically integrated streamwise mass flux
     &wnd_frc_slt, // I [m s-1] Saltating friction velocity
     &wnd_frc_thr_slt); // I [m s-1] Threshold friction speed for saltation

  prc_cmp lnd_frc_mbl; // [frc] Bare ground fraction
  prc_cmp lat(lat_dgr); // [dgr] Latitude
  prc_cmp lat_rdn(mth::cst_M_PIl*lat_dgr/180.0); // [rdn] Latitude
  rcd+=lnd_frc_mbl_get
    (lon_nbr, // I [nbr] Size of arrays
     doy, // I [day] Day of year [1.0..367.0)
     lat_rdn, // I [rdn] Latitude
     &lnd_frc_dry, // I [frc] Dry land fraction
     &lnd_frc_mbl, // O [frc] Bare ground fraction
     &oro, // I [frc] Orography
     &sfc_typ, // I [idx] LSM surface type (0..28)
     &snw_frc); // I [frc] Fraction of surface covered by snow

  // Apply land surface and vegetation limitations
  flx_mss_hrz_slt_ttl*=lnd_frc_mbl; // [kg m-1 s-1] Vertically integrated streamwise mass flux

  // Vertical dust mass flux
  // NB: Computing flx_mss_vrt_dst is done in full dust model but not here because of complexity
  //  prc_cmp *flx_mss_vrt_dst=new prc_cmp[sz_nbr]; // [kg m-2 s-1] Vertical mass flux of dust
  prc_cmp flx_mss_vrt_dst_ttl; // [kg m-2 s-1] Total vertical mass flux of dust
  prc_cmp dst_slt_flx_rat_ttl; // [m-1] Ratio of vertical dust flux to streamwise mass flux
  rcd+=flx_mss_vrt_dst_ttl_MaB95_get
    (lon_nbr, // I [nbr] Size of arrays
     &dst_slt_flx_rat_ttl, // O [m-1] Ratio of vertical dust flux to streamwise mass flux
     &flx_mss_hrz_slt_ttl, // I [kg m-1 s-1] Vertically integrated streamwise mass flux
     &flx_mss_vrt_dst_ttl, // O [kg m-2 s-1] Total vertical mass flux of dust
     &mss_frc_cly); // I [frc] Mass fraction clay 

  bool vgt; // [flg] "Vegetated" flag
  rcd+=lnd_is_vgt
    (lon_nbr, // I [nbr] Size of arrays
     &pnt_typ_idx, // I [idx] Plant type index 
     &snw_hgt, // I [m] Geometric bulk thickness of snow
     &vgt); // O [flg] "Vegetated" flag

  // Everything else    
  using phc::msv_snw_std; // [frc] Emissivity of snow CCM:lsm/snoconi
  if(msv_snw == cmd_ln_dfl) msv_snw=msv_snw_std; // [frc] Snow emissivity

  prc_cmp cff_xch_heat; // [frc] Exchange coefficient for heat transfer
  prc_cmp cff_xch_mmn; // [frc] Exchange coefficient for momentum transfer
  prc_cmp cff_xch_mmn_ntr; // [frc] Neutral drag coefficient hgt_mdp to z0m+zpd 
  prc_cmp cff_xch_vpr; // [frc] Exchange coefficient for vapor transfer
  prc_cmp flx_LW_upw_sfc; // [W m-2] Longwave upwelling flux at surface
  prc_cmp flx_ltn; // [W m-2] Latent heat flux to atmosphere
  prc_cmp flx_q_H2O; // [kg m-2 s-1] Moisture flux at surface
  prc_cmp flx_sns_atm; // [W m-2] Sensible heat flux to atmosphere
  prc_cmp flx_sns_gnd; // [W m-2] Sensible heat flux to soil
  prc_cmp flx_snw_mlt; // [W m-2] Snow melt heat flux
  prc_cmp msv_sfc; // [frc] Surface (bare ground+snow) emissivity
  prc_cmp rss_aer_heat; // [s m-1] Aerodynamic resistance to heat transfer
  prc_cmp rss_aer_mmn; // [s m-1] Aerodynamic resistance to momentum transfer
  // 20210418: fill-in rss_aer_mmn_ntr
  prc_cmp rss_aer_mmn_ntr(cmd_ln_dfl); // [s m-1] Neutral aerodynamic resistance to momentum transfer
  prc_cmp rss_aer_vpr; // [s m-1] Aerodynamic resistance to vapor transfer
  prc_cmp tpt_aer; // [K] "Aerodynamic" temperature at z=zpd+rgh_mmn
  prc_cmp tpt_ash; // [K] "Surface" temperature at z=zpd+rgh_heat
  prc_cmp tpt_ash_p2m; // [K] "Screen" temperature at z=zpd+rgh_heat+2m
  prc_cmp tpt_ffc; // [K] Radiative effective temperature
  prc_cmp wnd_str_mrd; // [kg m-1 s-2] Meridional wind stress
  prc_cmp wnd_str_znl; // [kg m-1 s-2] Zonal wind stress
  rcd+=flx_sfc_lnd
    (lon_nbr, // I [nbr] Size of arrays
     &vgt, // I [flg] "Vegetated" flag
     &cnd_trm_soi, // I [W m-1 K-1] Soil thermal conductivity
     &dns_mdp, // I [kg m-3] Midlayer density
     &flx_LW_dwn_sfc, // I [W m-2] Longwave downwelling flux at surface
     &flx_SW_net_gnd, // I [W m-2] Solar flux absorbed by ground
     &flx_SW_net_vgt, // I [W m-2] Solar flux absorbed by vegetation
     &hgt_mdp, // I [m] Midlayer height above surface
     &hgt_zpd_mbl, // I [m] Zero plane displacement
     &lvl_dlt_snw, // I [m] Soil layer thickness including snow
     &msv_gnd, // I [frc] Bare ground emissivity
     &msv_snw, // I [frc] Snow emissivity
     &prs_mdp, // I [Pa] Pressure
     &q_H2O_vpr, // I [kg kg-1] Specific humidity
     &rgh_mmn_mbl, // I [m] Roughness length momentum
     &snw_frc, // I [frc] Fraction of surface covered by snow
     &snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
     &tpt_mdp, // I [K] Midlayer temperature
     &tpt_ptn_mdp, // I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
     &tpt_soi, // I [K] Soil temperature
     &trn_fsh_vpr_soi_atm, // I [frc] Transfer efficiency of vapor from soil to atmosphere
     &wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     &wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     tm_dlt, // I [s] Timestep
     &pnt_typ_idx, // I [idx] Plant type index (1..14)
     &soi_typ, // I [idx] LSM soil type (1..5)
     &cff_xch_heat, // O [frc] Exchange coefficient for heat transfer
     &cff_xch_mmn, // O [frc] Exchange coefficient for momentum transfer
     &cff_xch_mmn_ntr, // O [frc] Neutral drag coefficient hgt_mdp to z0m+zpd 
     &cff_xch_vpr, // O [frc] Exchange coefficient for vapor transfer
     &flx_LW_upw_sfc, // O [W m-2] Longwave upwelling flux at surface
     &flx_ltn, // O [W m-2] Latent heat flux to atmosphere
     &flx_q_H2O, // O [kg m-2 s-1] Moisture flux to atmosphere
     &flx_sns_atm, // O [W m-2] Sensible heat flux to atmosphere
     &flx_sns_gnd, // O [W m-2] Sensible heat flux to soil
     &flx_snw_mlt, // O [W m-2] Snow melt heat flux
     &mno_lng_mbl, // I/O [m] Monin-Obukhov length
     &msv_sfc, // O [frc] Surface (bare ground+snow) emissivity
     &rss_aer_heat, // O [s m-1] Aerodynamic resistance to heat transfer
     &rss_aer_mmn, // O [s m-1] Aerodynamic resistance to momentum transfer
     &rss_aer_vpr, // O [s m-1] Aerodynamic resistance to vapor transfer
     &tpt_aer, // O [K] "Aerodynamic" temperature at z=zpd+rgh_mmn
     &tpt_ash, // O [K] "Surface" temperature at z=zpd+rgh_heat
     &tpt_ash_p2m, // O [K] "Screen" temperature at z=zpd+rgh_heat+2m
     &tpt_ffc, // O [K] Radiative effective temperature
     &tpt_gnd, // I/O [K] Ground temperature
     &tpt_vgt, // I/O [K] Vegetation temperature
     &wnd_frc_mbl, // O [m s-1] Surface friction velocity
     &wnd_str_mrd, // O [kg m-1 s-2] Meridional wind stress
     &wnd_str_znl); // O [kg m-1 s-2] Zonal wind stress

  // Boundary layer meteorology for deposition processes
  prc_cmp hgt_zpd_dps; // [m] Zero plane displacement
  prc_cmp mno_lng_dps; // [m] Monin-Obukhov length
  prc_cmp rgh_mmn_dps; // [m] Roughness length momentum
  prc_cmp wnd_frc_dps; // [m s-1] Surface friction velocity
  prc_cmp wnd_rfr_dps; // [m s-1] Wind speed at reference height
  // Surface energy fluxes in deposition region
  prc_cmp flx_LW_upw_sfc_dps; // [W m-2] Longwave upwelling flux at surface
  prc_cmp flx_ltn_dps; // [W m-2] Latent heat flux to atmosphere
  prc_cmp flx_q_H2O_dps; // [kg m-2 s-1] Moisture flux to atmosphere
  prc_cmp flx_sns_atm_dps; // [W m-2] Sensible heat flux to atmosphere
  prc_cmp wnd_str_mrd_dps; // [kg m-1 s-2] Meridional wind stress
  prc_cmp wnd_str_znl_dps; // [kg m-1 s-2] Zonal wind stress

  rcd+=blm_glb // Solve boundary layer meteorology on global scale
    (lon_nbr, // I [nbr] Size of arrays
     &dns_mdp, // I [kg m-3] Midlayer density
     &hgt_mdp, // I [m] Midlayer height above surface
     &oro, // I [frc] Orography
     &prs_mdp, // I [Pa] Pressure
     &q_H2O_vpr, // I [kg kg-1] Specific humidity
     &sfc_typ, // I [idx] LSM surface type (0..28)
     &snw_hgt_lqd, // I [m] Equivalent liquid water snow depth
     &tpt_mdp, // I [K] Midlayer temperature
     &tpt_ptn_mdp, // I [K] Potential temperature
     &tpt_sfc, // I [K] Surface temperature
     &wnd_mrd_mdp, // I [m s-1] Surface layer meridional wind speed
     &wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
     &flx_LW_upw_sfc_dps, // O [W m-2] Longwave upwelling flux at surface
     &flx_ltn_dps, // O [W m-2] Latent heat flux to atmosphere
     &flx_q_H2O_dps, // O [kg m-2 s-1] Moisture flux to atmosphere
     &flx_sns_atm_dps, // O [W m-2] Sensible heat flux to atmosphere
     &hgt_zpd_dps, // O [m] Zero plane displacement
     &mno_lng_dps, // O [m] Monin-Obukhov length
     &rgh_mmn_dps, // O [m] Roughness length momentum
     &wnd_frc_dps, // O [m s-1] Surface friction velocity
     &wnd_rfr_dps, // O [m s-1] Wind speed at reference height
     &wnd_str_mrd_dps, // O [kg m-1 s-2] Meridional wind stress
     &wnd_str_znl_dps); // O [kg m-1 s-2] Zonal wind stress
  assert(rgh_mmn_dps > 0.0 && rgh_mmn_dps < 10.0);
  assert(hgt_zpd_dps >= 0.0 && hgt_zpd_dps < 30.0 && hgt_zpd_dps < hgt_mdp);
  if(dbg_lvl >= dbg_sbr) std::cerr << prg_nm << ": INFO hgt_rfr = " << hgt_rfr << " m, wnd_mdp = " << wnd_mdp << " m s-1, wnd_rfr_dps = " << wnd_rfr_dps << " m s-1" << std::endl;
  assert(wnd_rfr_dps >= 0.0 && (wnd_rfr_dps < wnd_mdp || wnd_mdp <= 1.0));
  if(hgt_zpd_dps_cmd_ln != cmd_ln_dfl) hgt_zpd_dps=hgt_zpd_dps_cmd_ln; // [m]
  if(rgh_mmn_dps_cmd_ln != cmd_ln_dfl) rgh_mmn_dps=rgh_mmn_dps_cmd_ln; // [m]
  if(mno_lng_dps_cmd_ln != cmd_ln_dfl) mno_lng_dps=mno_lng_dps_cmd_ln; // [m]
  if(wnd_frc_dps_cmd_ln != cmd_ln_dfl) wnd_frc_dps=wnd_frc_dps_cmd_ln; // [m s-1]

  // Aerodynamic resistance
  cpv_foo=rss_aer_mmn; // [s m-1]
  rcd+=rss_aer_get
    (lon_nbr, // I [nbr] Size of arrays
     &hgt_mdp, // I [m] Midlayer height above surface
     &hgt_zpd_dps, // I [m] Zero plane displacement height
     &mno_lng_dps, // I [m] Monin-Obukhov length
     &rgh_mmn_dps, // I [m] Roughness length momentum
     &rss_aer_mmn, // O [s m-1] Aerodynamic resistance
     &wnd_frc_dps); // I [m s-1] Surface friction velocity
  // Check that SeP97 method agrees with Bon96 method
  //  assert(std::fabs((rss_aer_mmn-cpv_foo)/rss_aer_mmn) < 1.0e-4);

  // Use settling velocities to compute Stokes number
  prc_cmp *dff_aer=new prc_cmp[sz_nbr]; // [m2 s-1] Brownian diffusivity of particle
  prc_cmp *rss_lmn=new prc_cmp[sz_nbr]; // [s m-1] Quasi-laminar layer resistance
  prc_cmp *shm_nbr=new prc_cmp[sz_nbr]; // [frc] Schmidt number
  prc_cmp *pcl_nbr=new prc_cmp[sz_nbr]; // [frc] Peclet number (for completeness)
  prc_cmp *stk_nbr=new prc_cmp[sz_nbr]; // [frc] Stokes number
  prc_cmp dff_aer_vwr(0.0); // [m2 s-1] Mass weighted Brownian diffusivity of particle
  prc_cmp dff_aer_nwr(0.0); // [m2 s-1] Number weighted Brownian diffusivity of particle
  prc_cmp rss_lmn_vwr(0.0); // [s m-1] Mass weighted laminar resistance
  prc_cmp rss_lmn_nwr(0.0); // [s m-1] Number weighted laminar resistance
  prc_cmp shm_nbr_vwr(0.0); // [frc] Mass weighted Schmidt number
  prc_cmp shm_nbr_nwr(0.0); // [frc] Number weighted Schmidt number
  prc_cmp stk_nbr_vwr(0.0); // [frc] Mass weighted Stokes number
  prc_cmp stk_nbr_nwr(0.0); // [frc] Number weighted Stokes number
  prc_cmp shm_nbr_xpn; // [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt number
  // Use shm_nbr^(-2/3) over solid surfaces and shm_nbr^(-1/2) over liquid surfaces SlS80 p. 1014
  if(oro_is_ocn(oro)) shm_nbr_xpn=-0.5; else shm_nbr_xpn=-0.6666667; // [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt number
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant
  for(idx=0;idx<sz_nbr;idx++){
    stk_nbr[idx]=vlc_grv[idx]*wnd_frc_dps*wnd_frc_dps/(grv_sfc_mean*vsc_knm_atm); // [frc] Stokes number SeP97 p. 965
    // Without slip correction factor, following is known as Stokes-Einstein relation
    dff_aer[idx]=cst_Boltzmann*tpt_mdp*slp_crc[idx]/(3.0*mth::cst_M_PIl*vsc_dyn_atm*dmt_ctr[idx]); // [m2 s-1] Brownian diffusivity SeP97 p. 474 (8.73), PeE92 p. 2558 (21)
    shm_nbr[idx]=vsc_knm_atm/dff_aer[idx]; // [frc] Schmidt number SeP97 p. 972
    pcl_nbr[idx]=vlc_grv[idx]*rds_ctr[idx]/dff_aer[idx]; // [frc] Peclet number Sli82 p. 323
    rss_lmn[idx]=1.0/(wnd_frc_dps*(std::pow(shm_nbr[idx],shm_nbr_xpn)+std::pow(PRC_CMP(10.0),PRC_CMP(-3.0)/stk_nbr[idx]))); // [s m-1] Quasi-laminar layer resistance SeP97 p. 972 (19.18), 965 (19.18), SlS80 p. 1014 (5), SHH78 p. 2087 (C.6), p. 2070 (43)
  } // end loop over sz

  // Weight by number and by mass
  for(idx=0;idx<sz_nbr;idx++){
    dff_aer_vwr+=dff_aer[idx]*vlm[idx]*cnc[idx]; // Mass weighted Brownian diffusivity of particle
    dff_aer_nwr+=dff_aer[idx]*cnc[idx]; // Number weighted Brownian diffusivity of particle
    rss_lmn_vwr+=rss_lmn[idx]*vlm[idx]*cnc[idx]; // Mass weighted laminar resistance
    rss_lmn_nwr+=rss_lmn[idx]*cnc[idx]; // Number weighted laminar resistance
    shm_nbr_vwr+=shm_nbr[idx]*vlm[idx]*cnc[idx]; // Mass weighted Schmidt number
    shm_nbr_nwr+=shm_nbr[idx]*cnc[idx]; // Number weighted Schmidt number
    stk_nbr_vwr+=stk_nbr[idx]*vlm[idx]*cnc[idx]; // Mass weighted Stokes number
    stk_nbr_nwr+=stk_nbr[idx]*cnc[idx]; // Number weighted Stokes number
  } // end loop over sz
  dff_aer_vwr/=vlm_rsl; // [m2 s-1] Mass weighted Brownian diffusivity of particle
  dff_aer_nwr/=cnc_nbr_rsl; // [m2 s-1] Number weighted Brownian diffusivity of particle
  rss_lmn_vwr/=vlm_rsl; // [s m-1] Mass weighted laminar resistance
  rss_lmn_nwr/=cnc_nbr_rsl; // [s m-1] Number weighted laminar resistance
  shm_nbr_vwr/=vlm_rsl; // [frc] Mass weighted Schmidt number
  shm_nbr_nwr/=cnc_nbr_rsl; // [frc] Number weighted Schmidt number
  stk_nbr_vwr/=vlm_rsl; // [frc] Mass weighted Stokes number
  stk_nbr_nwr/=cnc_nbr_rsl; // [frc] Number weighted Stokes number

  // Combine resistances and settling velocity to determine total deposition velocity
  prc_cmp *rss_trb=new prc_cmp[sz_nbr]; // [s m-1] Resistance to turbulent deposition
  prc_cmp *vlc_dry=new prc_cmp[sz_nbr]; // [m s-1] Total dry deposition velocity
  prc_cmp *vlc_trb=new prc_cmp[sz_nbr]; // [m s-1] Turbulent deposition velocity
  for(idx=0;idx<sz_nbr;idx++){
    rss_trb[idx]=rss_aer_mmn+rss_lmn[idx]+rss_aer_mmn*rss_lmn[idx]*vlc_grv[idx]; // [s m-1] Resistance to turbulent deposition SeP97 p. 961 (19.7)
    vlc_trb[idx]=1.0/rss_trb[idx]; // [m s-1] Turbulent deposition velocity SeP97 p. 961 (19.7)
    vlc_dry[idx]=vlc_grv[idx]+vlc_trb[idx]; // [m s-1] Total dry deposition velocity SeP97 p. 961 (19.7)
  } // end loop over sz

  // Wet deposition
  // Instantiate raindrop size grid
  const prc_cmp dns_pcp(dns_H2O_lqd_std); // [kg m-3] Precipitation density
  prc_cmp dsd_mnm(dsd_mnm_mcr*1.0e-6); // [um]->[m] Minimum diameter in raindrop distribution
  prc_cmp dsd_mxm(dsd_mxm_mcr*1.0e-6); // [um]->[m] Maximum diameter in raindrop distribution
  prc_cmp dsd_dbg=dsd_dbg_mcr*1.0e-6; // [um]->[m] Debugging size for raindrops
  // Instantiate raindrop size grid
  sz_grd_cls dsd_szgrd(dsd_mnm,dsd_mxm,dsd_nbr,sz_grd_sng); // [m] Raindrop size grid object
  const prc_cmp *dmt_pcp=dsd_szgrd.sz_ctr_get(); // [m] Diameter at bin center, raindrop
  const prc_cmp *dsd_sz=dsd_szgrd.sz_ctr_get(); // [m] Diameter at bin center, raindrop
  const long dsd_idx_dbg=vec_val2idx(dmt_pcp,dsd_nbr,dsd_dbg); // [idx] Size bin for debugging, raindrop

  // Instantiate default raindrop size distribution
  prc_cmp dmt_pcp_nma(dmt_pcp_nma_mcr*1.0e-6); // [um]->[m] Diameter number median analytic, raindrop
  prc_cmp rds_pcp_nma(dmt_pcp_nma/2.0); // [m] Radius number median analytic, raindrop
  psd_cls psd_pcp(rds_pcp_nma,gsd_pcp_anl,cnc_nbr_pcp_anl); // [sct] Raindrop size distribution
  rcd+=psd_pcp.dns_prt_set(dns_pcp); // [kg m-3] Particle density
  rcd+=psd_pcp.mss_frc_anl_set(1.0); // [frc] Mass fraction analytic
  LogNormalFunction<prc_cmp> dst_lgn_pcp;
  // Set precipitation distribution to be per unit diameter
  rcd+=dst_lgn_pcp.abc_typ_set("dmt");  // [fnc] Set abc_typ to rds_typ or dmt_typ
  rcd+=dst_lgn_pcp.prm_rst(psd_pcp); // [fnc] Reset parameters of distribution function
  // Override defaults with user-specified mode, if any
  if(dsd_usr_flg) rcd+=psd_pcp.prs_psd_sng(dsd_arg);

  // Useful fundamental properties of raindrop size distribution
  const prc_cmp *dmt_dlt_pcp=dsd_szgrd.sz_dlt_get(); // [m] Width of size bin, raindrop
  prc_cmp *rds_pcp=new prc_cmp[dsd_nbr]; // [m] Radius at bin center, raindrop
  prc_cmp *sfc_pcp=new prc_cmp[dsd_nbr]; // [m2] Surface area of precipitation
  prc_cmp *xsa_pcp=new prc_cmp[dsd_nbr]; // [m2] Cross-sectional area of precipitation
  prc_cmp *mss_pcp=new prc_cmp[dsd_nbr]; // [kg] Precipitation mass
  prc_cmp *vlm_pcp=new prc_cmp[dsd_nbr]; // [m3] Precipitation volume
  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++){
    rds_pcp[dsd_idx]=0.5*dmt_pcp[dsd_idx]; // [m] Radius at bin center, raindrop
    sfc_pcp[dsd_idx]=mth::cst_M_PIl*dmt_pcp[dsd_idx]*dmt_pcp[dsd_idx]; // [m2] Surface area of precipitation
    xsa_pcp[dsd_idx]=mth::cst_M_PIl*dmt_pcp[dsd_idx]*dmt_pcp[dsd_idx]*0.25; // [m2] Cross-sectional area of precipitation
    vlm_pcp[dsd_idx]=mth::cst_M_PIl*std::pow(dmt_pcp[dsd_idx],PRC_CMP(3.0))/6.0; // [m3] Precipitation volume
    mss_pcp[dsd_idx]=mth::cst_M_PIl*std::pow(dmt_pcp[dsd_idx],PRC_CMP(3.0))*dns_pcp/6.0; // [kg] Precipitation mass
  } // end loop over dsd

  // Stokes' gravitational settling velocity, raindrops
  prc_cmp *slp_crc_pcp=new prc_cmp[dsd_nbr]; // [frc] Slip correction factor SeP97 of raindrop
  prc_cmp *vlc_stk_pcp=new prc_cmp[dsd_nbr]; // [m s-1] Stokes' settling velocity (Re < 0.1) of raindrop
  rcd+=vlc_stk_get
    (dsd_nbr, // I [nbr] Number of size bins
     dmt_pcp, // I [m] Raindrop diameter
     dns_pcp, // I [kg m-3] Raindrop density
     mfp_atm, // I [m] Mean free path of atmosphere
     vsc_dyn_atm, // I [kg m-1 s-1] Dynamic viscosity of atmosphere
     slp_crc_pcp, // O [frc] Slip correction factor of raindrop
     vlc_stk_pcp); // O [m s-1] Stokes' settling velocity (Re < 0.1) of raindrop

  // Raindrop gravitational settling velocity
  prc_cmp *cff_drg_pcp=new prc_cmp[dsd_nbr]; // [frc] Drag coefficient at terminal velocity of raindrop
  prc_cmp *ryn_nbr_pcp=new prc_cmp[dsd_nbr]; // [frc] Reynolds number at terminal velocity of raindrop
  prc_cmp *stk_crc_pcp=new prc_cmp[dsd_nbr]; // [frc] Correction to Stokes settling velocity of raindrop
  prc_cmp *vlc_pcp=new prc_cmp[dsd_nbr]; // [m s-1] Settling velocity of raindrop
  rcd+=vlc_grv_get
    (dsd_nbr, // I [nbr] Size of arrays
     dmt_pcp,  // I [m] Raindrop diameter
     slp_crc_pcp, // I [frc] Slip correction factor SeP97 of raindrop
     vlc_stk_pcp, // I [m s-1] Stokes' settling velocity (Re < 0.1) of raindrop
     dns_mdp, // I [kg m-3] Midlayer density
     dns_pcp, // I [kg m-3] Raindrop density
     vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
     cff_drg_pcp, // O [frc] Drag coefficient at terminal velocity of raindrop
     ryn_nbr_pcp, // O [frc] Reynolds number at terminal velocity of raindrop
     stk_crc_pcp, // O [frc] Correction to Stokes settling velocity of raindrop
     vlc_pcp); // O [m s-1] Settling velocity of raindrop
  // Free un-needed memory from recent calculations
  delete []cff_drg_pcp; // [frc] Drag coefficient at terminal velocity of raindrop
  delete []slp_crc_pcp; // [frc] Slip correction factor SeP97 of raindrop
  delete []stk_crc_pcp; // [frc] Correction to Stokes settling velocity of raindrop

  // Convolve size grid with particle distribution
  prc_cmp *dst_pcp=new prc_cmp[dsd_nbr]; // [# m-3 m-1] Number distribution, raindrops
  (void)vec_set(dst_pcp,dsd_nbr,0.0);
  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++) dst_pcp[dsd_idx]+=dst_lgn_pcp(dmt_pcp[dsd_idx]); // [# m-3 m-1] Raindrop number distribution

  // Compute precipitation intensity with provisional raindrop number distribution
  prc_cmp flx_vlm_pcp_rsl_prv(0.0); // [m3 m-2 s-1]=[m s-1] Raindrop volume flux, provisional
  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++) flx_vlm_pcp_rsl_prv+=vlm_pcp[dsd_idx]*vlc_pcp[dsd_idx]*dst_pcp[dsd_idx]*dmt_dlt_pcp[dsd_idx]; // [m3 m-2 s-1]=[m s-1] Raindrop volume flux, provisional

  // If precipitation intensity was specified...
  if(flx_vlm_pcp_rsl > 0.0){
    // Normalize provisional precipitation intensity to specified intensity by linearly adjusting distribution
    const prc_cmp pcp_rat(flx_vlm_pcp_rsl/flx_vlm_pcp_rsl_prv); // [frc] Ratio of provisional to specified precipitation intensity
    cnc_nbr_pcp_anl*=pcp_rat; // [# m-3] Raindrop number concentration, analytic
    rcd+=psd_pcp.cnc_nbr_anl_set(cnc_nbr_pcp_anl); // [# m-3] Raindrop number concentration, analytic
    rcd+=dst_lgn_pcp.cnc_nbr_anl_set(cnc_nbr_pcp_anl); // [# m-3] Raindrop number concentration, analytic
    for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++) dst_pcp[dsd_idx]*=pcp_rat; // [# m-3 m-1] Raindrop number distribution
  } // endif flx_vlm_pcp_rsl > 0.0
  
  // Precipitation diagnostics
  prc_cmp flx_mss_pcp_rsl(0.0); // [kg m-2 s-1] Raindrop mass flux, resolved
  prc_cmp flx_cnc_pcp_rsl(0.0); // [# m-2 s-1] Raindrop number flux, resolved
  prc_cmp cnc_nbr_pcp_rsl(0.0); // [# m-3] Raindrop number concentration, resolved
  prc_cmp sfc_pcp_rsl(0.0); // [m2 m-3] Raindrop surface area concentration resolved
  prc_cmp xsa_pcp_rsl(0.0); // [m2 m-3] Raindrop cross-sectional area concentration resolved
  prc_cmp vlm_pcp_rsl(0.0); // [m3 m-3] Raindrop volume concentration resolved
  prc_cmp mss_pcp_rsl(0.0); // [kg m-3] Raindrop mass concentration resolved
  prc_cmp *cnc_pcp=new prc_cmp[dsd_nbr]; // [# m-3] Raindrop number concentration
  prc_cmp *flx_cnc_spc_pcp=new prc_cmp[dsd_nbr]; // [# m-2 s-1 m-1] Raindrop spectral number flux
  prc_cmp *flx_mss_spc_pcp=new prc_cmp[dsd_nbr]; // [kg m-2 s-1 m-1] Raindrop spectral mass flux
  prc_cmp *flx_vlm_spc_pcp=new prc_cmp[dsd_nbr]; // [m3 m-2 s-1 m-1] Raindrop spectral volume flux
  flx_vlm_pcp_rsl=0.0; // [m3 m-2 s-1]=[m s-1] Raindrop volume flux, resolved
  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++){
    cnc_pcp[dsd_idx]=dst_pcp[dsd_idx]*dmt_dlt_pcp[dsd_idx]; // [# m-3] Raindrop number concentration

    flx_cnc_spc_pcp[dsd_idx]=vlc_pcp[dsd_idx]*dst_pcp[dsd_idx]; // [# m-2 s-1 m-1] Raindrop spectral number flux
    flx_vlm_spc_pcp[dsd_idx]=vlm_pcp[dsd_idx]*flx_cnc_spc_pcp[dsd_idx]; // [m3 m-2 s-1]=[m s-1] Raindrop volume flux, provisional
    flx_mss_spc_pcp[dsd_idx]=dns_pcp*flx_vlm_spc_pcp[dsd_idx]; // [kg m-2 s-1 m-1] Raindrop spectral mass flux
    cnc_nbr_pcp_rsl+=cnc_pcp[dsd_idx]; // [# m-3] Raindrop number concentration, resolved
    sfc_pcp_rsl+=sfc_pcp[dsd_idx]*cnc_pcp[dsd_idx]; // [m2 m-3] Raindrop surface area concentration resolved
    xsa_pcp_rsl+=xsa_pcp[dsd_idx]*cnc_pcp[dsd_idx]; // [m2 m-3] Raindrop cross-sectional area concentration resolved
    vlm_pcp_rsl+=vlm_pcp[dsd_idx]*cnc_pcp[dsd_idx]; // [m3 m-3] Raindrop volume concentration resolved
    flx_vlm_pcp_rsl+=flx_vlm_spc_pcp[dsd_idx]*dmt_dlt_pcp[dsd_idx]; // [m3 m-2 s-1]=[m s-1] Raindrop volume flux, resolved
    flx_cnc_pcp_rsl+=flx_cnc_spc_pcp[dsd_idx]*dmt_dlt_pcp[dsd_idx]; // [# m-2 s-1] Raindrop number flux, resolved
  } // end loop over dsd
  flx_mss_pcp_rsl=dns_pcp*flx_vlm_pcp_rsl; // [kg m-2 s-1] Raindrop mass flux, resolved
  mss_pcp_rsl=dns_pcp*vlm_pcp_rsl; // [kg m-3] Raindrop mass concentration resolved
  prc_cmp dmt_pcp_vwr_Mas71; // [m] Diameter volume weighted resolved, raindrop, Mas71 parameterzation
  dmt_pcp_vwr_Mas71=0.0305*std::pow(flx_vlm_pcp_rsl,PRC_CMP(0.25)); // [m] Diameter volume weighted resolved, raindrop Mas71, Sli84 p. 476 (11.33)

  // Resolved fractions (i.e., discretization efficiency)
  const prc_cmp cnc_nbr_rsl_frc_pcp(cnc_nbr_pcp_rsl/cnc_nbr_pcp_anl); // [frc] Resolved fraction of number concentration, raindrops
  const prc_cmp sfc_rsl_frc_pcp(sfc_pcp_rsl/dst_lgn_pcp.sfc_get()); // [frc] Resolved fraction of surface area concentration, raindrops
  const prc_cmp xsa_rsl_frc_pcp(xsa_pcp_rsl/dst_lgn_pcp.xsa_get()); // [frc] Resolved fraction of cross-sectional area concentration, raindrops
  const prc_cmp vlm_rsl_frc_pcp(vlm_pcp_rsl/dst_lgn_pcp.vlm_get()); // [frc] Resolved fraction of volume concentration, raindrops

  // Distribution diagnostics
  prc_cmp vlc_pcp_nwr(0.0); // [m s-1] Number weighted terminal velocity, raindrops
  prc_cmp vlc_pcp_vwr(0.0); // [m s-1] Mass weighted terminal velocity, raindrops
  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++){
    vlc_pcp_nwr+=vlc_pcp[dsd_idx]*cnc_pcp[dsd_idx]; // Number weighted terminal velocity, raindrops
    vlc_pcp_vwr+=vlc_pcp[dsd_idx]*vlm_pcp[dsd_idx]*cnc_pcp[dsd_idx]; // Mass weighted terminal velocity, raindrops
  } // end loop over dsd
  vlc_pcp_nwr/=cnc_nbr_pcp_rsl; // [m s-1] Number weighted terminal velocity, raindrops
  vlc_pcp_vwr/=vlm_pcp_rsl; // [m s-1] Mass weighted terminal velocity, raindrops

  // Calculation of collision efficiency E discussed on SeP97 p. 1019, DaH76, NGD94
  // All modern sources seem based on work of W. G. N. Slinn: SHH78 Sli82, Sli83, Sli84
  prc_cmp dmt_rat;  // [frc] Ratio of collectee to collector sizes
  a2d_cls<prc_cmp> cll_fsh(dsd_nbr,sz_nbr); // [frc] Collision efficiency
  a2d_cls<prc_cmp> clc_fsh(dsd_nbr,sz_nbr); // [frc] Collection efficiency
  prc_cmp *stc_fsh=new prc_cmp[sz_nbr]; // [frc] Sticking (coalescence) efficiency
  // The viscosity ratio used in PrK98 replaces the ratio of raindrop internal circulation to raindrop fall speed used in Sli82
  // Results for the latter are not given in Sli82, but are discussed extensively in PrK98 Section 10.3.1 pp. 386--393
  using phc::vsc_dyn_H2O; // (0.001) [kg m-1 s-1] Dynamic viscosity of liquid H2O PrK98 p. 387
  const prc_cmp vsc_rat(vsc_dyn_H2O/vsc_dyn_atm); // [frc] Ratio of liquid H2O to atmospheric dynamic viscosity
  assert(vsc_rat > 20.0 && vsc_rat < 100.0); // Based on PrK98 p. 387 discussion
  const prc_cmp dns_rat(dns_pcp/dns_prt); // [frc] Ratio of raindrop density to aerosol density
  // Parameterization for critical Stokes number is discussed in Sli82 p. 332
  // Formula comes from fit to data of Beard and Grover (1974)
  prc_cmp *stk_nbr_crt=new prc_cmp[dsd_nbr]; // [frc] Critical Stokes number
  a2d_cls<prc_cmp> stk_nbr_rlt(dsd_nbr,sz_nbr); // [frc] Stokes number of relative flow
  prc_cmp *tau_rlx=new prc_cmp[sz_nbr]; // [s] Relaxation timescale
  a2d_cls<prc_cmp> cll_fsh_brn_dff(dsd_nbr,sz_nbr); // [frc] Collision efficiency of Brownian diffusion
  a2d_cls<prc_cmp> cll_fsh_ntc(dsd_nbr,sz_nbr); // [frc] Collision efficiency of interception
  a2d_cls<prc_cmp> cll_fsh_mpc(dsd_nbr,sz_nbr); // [frc] Collision efficiency of impaction
  prc_cmp mpc_fct; // [frc] Factor in impaction contribution to collision efficiency
  // Relaxation timescale depends only on particle and atmospheric properties
  for(idx=0;idx<sz_nbr;idx++){
    // Relaxation timescale is approximate time (truly valid for Stokes flow only) for particle to reach terminal velocity from standing start
    tau_rlx[idx]=mss[idx]*slp_crc[idx]/(3.0*mth::cst_M_PIl*vsc_dyn_atm*dmt_ctr[idx]); // [s] Relaxation timescale Sep97 p. 465 (8.38) (validate against Table 8.4)
  } // end loop over sz

  prc_cmp *vlc_dry_SeH78=new prc_cmp[sz_nbr]; // [m s-1] Total dry deposition velocity based on SeH78
  prc_cmp *vlc_trb_SeH78=new prc_cmp[sz_nbr]; // [m s-1] Turbulent deposition velocity based on SeH78
  prc_cmp tau_rlx_scl; // [frc] Non-dimensional relaxation timescale
  prc_cmp int_3; // [frc] Resistance integral Int_3 from SeH78
  for(idx=0;idx<sz_nbr;idx++){
    // fxm: how to non-dimensionalize relaxation timescale?
    tau_rlx_scl=tau_rlx[idx]*dns_prt*wnd_frc_dps*wnd_frc_dps/vsc_knm_atm; // [frc] Non-dimensional relaxation timescale SeH782 p. 2.12
    int_3=-std::exp(-378.051+16.498*std::log(shm_nbr[idx])-12.804*std::log(dmt_ctr[idx])
	       +std::log(tau_rlx_scl)*(-11.818-0.2863*std::log(tau_rlx_scl)+0.3226*std::log(dmt_ctr[idx]/rgh_mmn_dps)
				   -0.3385*std::log(dff_aer[idx]/(rgh_mmn_dps*wnd_frc_dps)))); 
    vlc_trb_SeH78[idx]=int_3; // [m s-1] Turbulent deposition velocity based on SeH78
    vlc_dry_SeH78[idx]=vlc_trb_SeH78[idx]+vlc_grv[idx]; // [m s-1] Total dry deposition velocity based on SeH78
  } // end loop over sz

  for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++){
    // Critical Stokes number depends only on raindrop properties
    stk_nbr_crt[dsd_idx]=(1.2+(1.0/12.0)*std::log(1.0+ryn_nbr_pcp[dsd_idx]))/(1.0+std::log(1.0+ryn_nbr_pcp[dsd_idx])); // [frc] Critical Stokes number Sli82 p. 329 and p. 331, SeP97 p. 1019 (20.57), NGD94 p. 2336 (8)

    // Properties that depend on both particle and precipitation properties
    for(idx=0;idx<sz_nbr;idx++){

      // Relative Stokes number depends on both aerosol and raindrop properties:
      // Use relative velocity and normalize by raindrop diameter
      stk_nbr_rlt(dsd_idx,idx)=2.0*tau_rlx[idx]*std::fabs((vlc_pcp[dsd_idx]-vlc_grv[idx]))/dmt_pcp[dsd_idx]; // [frc] Stokes number of relative flow SeP97 p. 1019
      
      // Collection due to Brownian Diffusion accounts for particles "random walking" across streamlines
      // Brownian diffusion is primarily important for particles D < 2 um
      cll_fsh_brn_dff(dsd_idx,idx)=4.0/(ryn_nbr_pcp[dsd_idx]*shm_nbr[idx]); // (frc] Collision efficiency of Brownian diffusion SeP97 p. 1019 (20.56) Sli82 p. 324 (45), NGD94 p. 2336 (10)
      // The last term of the Brownian diffusion expression attempts to account for the internal circulation of raindrops
      cll_fsh_brn_dff(dsd_idx,idx)*=1.0+0.4*std::sqrt(ryn_nbr_pcp[dsd_idx])*std::pow(shm_nbr[idx],PRC_CMP(1.0)/PRC_CMP(3.0))+0.16*std::sqrt(ryn_nbr_pcp[dsd_idx]*shm_nbr[idx]); // [frc] Collision efficiency of Brownian diffusion SeP97 p. 1019 (20.56), Sli82 p. 324 (45) and p. 331 caption of Figure 11, NGD94 p. 2336 (10)
      
      // Collection due to impaction accounts for particles which are unable to follow the streamlines of fluid flow around an approaching raindrop
      // Due to their large size (Stokes number) and the sharpness of the streamlines, these particles are unable to move out of the way of the raindrops so collection occurs
      // Contribution from impaction should vanish unless St > St* Sli82 p. 332
      // However, most references show a continuous contribution from impaction, e.g., DaH76 p. 46 Fgr. 1 
      // Impaction is primarily important for particles D > 2 um
      cll_fsh_mpc(dsd_idx,idx)=0.0; // [frc] Collision efficiency of impaction SeP97 p. 1019 (20.56)
      if (stk_nbr_rlt(dsd_idx,idx) >= stk_nbr_crt[dsd_idx]){
	mpc_fct=(stk_nbr_rlt(dsd_idx,idx)-stk_nbr_crt[dsd_idx])/(stk_nbr_rlt(dsd_idx,idx)-stk_nbr_crt[dsd_idx]+2.0/3.0); // [frc] Factor in impaction contribution to collision efficiency SeP97 p. 1019 (20.56) Sli82 p. 331, NGD94 p. 2336 (9)
	cll_fsh_mpc(dsd_idx,idx)=std::sqrt(dns_rat)*std::pow(mpc_fct,PRC_CMP(1.5)); // [frc] Collision efficiency of impaction SeP97 p. 1019 (20.56), NGD94 p. 2336 (9)
      }else{
	;
	// Older formula for impaction which preserves continuous behavior
	// This formula is superceded by the critical Stokes number theory of Slinn
	//      mpc_fct=(stk_nbr_rlt(dsd_idx,idx)-1.0/12.0)/(stk_nbr_rlt(dsd_idx,idx)+7.0/12.0); // [frc] Factor in impaction contribution to collision efficiency DaH76 p. 46 Table 1
	//      cll_fsh_mpc(dsd_idx,idx)=std::sqrt(dns_rat)*std::pow(mpc_fct,PRC_CMP(1.5)); // [frc] Collision efficiency of impaction DaH76 p. 46 Table 1
      } // endif
      
      /* Collection due to interception accounts for particles which follow the streamlines of flow around an approaching raindrop
	 When the particle radius is larger than the distance of the streamline to the raindrop, interception occurs
	 Thus interception is strictly due to particle size, not mass
	 Interception is primarily important for particles 1 < D < 5 um
	 Ratio of sizes depends on relative size of aerosol and precipitation
	 Smaller particle ("collectee") must be in numerator */
      if(dmt_ctr[idx] < dmt_pcp[dsd_idx]) dmt_rat=dmt_ctr[idx]/dmt_pcp[dsd_idx]; else dmt_rat=dmt_pcp[dsd_idx]/dmt_ctr[idx]; // [frc] Ratio of collectee to collector sizes
      cll_fsh_ntc(dsd_idx,idx)=4.0*dmt_rat*((1.0/vsc_rat)+dmt_rat*(1.0+2.0*std::sqrt(ryn_nbr_pcp[dsd_idx]))); // [frc] Collision efficiency of interception SeP97 p. 1019 (20.56), NGD94 p. 2336 (8)
      // If interception efficiency is too large limit it to Fuchs (1964) result for flow about a sphere
      if(cll_fsh_ntc(dsd_idx,idx) > 1.0) cll_fsh_ntc(dsd_idx,idx)=3.0*dmt_rat; // [frc] Collision efficiency of interception Sli82 p. 325 (47), DaH76 p. 46 Tbl. 1
      // If interception efficiency is still too large then limit it to pure geometric result
      if(cll_fsh_ntc(dsd_idx,idx) > 1.0) cll_fsh_ntc(dsd_idx,idx)=dmt_rat*dmt_rat; // [frc] Collision efficiency of interception Sli82 p. 325 (46)
      
      // Sum diffusion, interception, and impaction contributions to total collision efficiency
      cll_fsh(dsd_idx,idx)=cll_fsh_brn_dff(dsd_idx,idx)+cll_fsh_ntc(dsd_idx,idx)+cll_fsh_mpc(dsd_idx,idx); // [frc] Collision efficiency SeP97 p. 1019 (20.56)
      
      // Adjust each contribution by same relative amount
      if(cll_fsh(dsd_idx,idx) > 1.0){
	// std::cerr << prg_nm << ": WARNING cll_fsh = " << cll_fsh(dsd_idx,idx) << " for dmt_ctr[" << idx << "] = " << dmt_ctr[idx] << ", reducing cll_fsh_brn_dff,cll_fsh_ntc,cll_fsh_mpc to compensate" << std::endl;
	// Interception efficiency seems the most logical to adjust
	cll_fsh_brn_dff(dsd_idx,idx)/=cll_fsh(dsd_idx,idx); // [frc] Collision efficiency of Brownian diffusion
	cll_fsh_ntc(dsd_idx,idx)/=cll_fsh(dsd_idx,idx); // [frc] Collision efficiency of interception
	cll_fsh_mpc(dsd_idx,idx)/=cll_fsh(dsd_idx,idx); // [frc] Collision efficiency of impaction
	cll_fsh(dsd_idx,idx)=1.0; // [frc] Collision efficiency
      } // endif
      //    std::cerr << "cll_fsh[" << idx << "] = " << cll_fsh(dsd_idx,idx) << std::endl;
      assert(cll_fsh(dsd_idx,idx) <= 1.0);
      
      // Sticking efficiency is probably unity for D < 20 um Sli82 p. 330
      stc_fsh[idx]=1.0; // [frc] Sticking (coalescence) efficiency
      clc_fsh(dsd_idx,idx)=stc_fsh[idx]*cll_fsh(dsd_idx,idx); // [frc] Collection efficiency
    } // end loop over sz
  } // end loop over dsd
  prc_cmp cll_vlm_xsa; // [m2] Cross-sectional area of geometric collision volume
  prc_cmp cll_vlm_vlc; // [m s-1] Velocity of geometric collision volume
  prc_cmp *scv_cff=new prc_cmp[sz_nbr]; // [s-1] Scavenging coefficient
  prc_cmp *scv_cff_pcp_nrm=new prc_cmp[sz_nbr]; // [m2 kg-1] Scavenging coefficient, precipitation-normalized
  prc_cmp scv_cff_mss_avg(0.0); // [s-1] Mass mean scavenging coefficient of aerosol
  prc_cmp scv_cff_mss_avg_pcp_nrm(0.0); // [m2 kg-1] Mass mean scavenging coefficient of aerosol, precipitation-normalized
  prc_cmp scv_cff_nbr_avg(0.0); // [s-1] Number mean scavenging coefficient of aerosol
  prc_cmp scv_cff_nbr_avg_pcp_nrm(0.0); // [m2 kg-1] Number mean scavenging coefficient of aerosol, precipitation-normalized
  /* Scavenging coefficient varies as square of raindrop diameter 
     Raindrop diameters are of course event dependent
     NGD94 p. 2337 Table 1 suggest dmt_nma=0.4 mm, gsd_anl=1.86 for drizzle and dmt_nma=1.0 mm, gsd_anl=1.86 for convective rain */
  for(idx=0;idx<sz_nbr;idx++){
    // Initalize integrals
    scv_cff[idx]=0.0; // [s-1] Scavenging coefficient
    for(dsd_idx=0;dsd_idx<dsd_nbr;dsd_idx++){
      cll_vlm_xsa=0.25*mth::cst_M_PIl*(dmt_pcp[dsd_idx]+dmt_ctr[idx])*(dmt_pcp[dsd_idx]+dmt_ctr[idx]); // [m2] Cross-sectional area of geometric collision volume
      cll_vlm_vlc=vlc_pcp[dsd_idx]-vlc_grv[idx]; // [m s-1] Velocity of geometric collision volume
      scv_cff[idx]+=cll_vlm_xsa*cll_vlm_vlc*clc_fsh(dsd_idx,idx)*cnc_pcp[dsd_idx]; // [s-1] Scavenging coefficient SeP97 p. 1021 (20.53, 20.58)
      // Standard approximation is (Dp + dp)^2 ~ Dp^2 and Vg-vg ~ Vg
      //      scv_cff[idx]+=xsa_pcp[dsd_idx]*vlc_pcp[dsd_idx]*clc_fsh[idx]*cnc_pcp[dsd_idx]; // [s-1] Scavenging coefficient SeP97 p. 1021 (20.53, 20.58)
    } // end loop over dsd
    // Normalize scavenging coefficient by precipitation mass flux
    scv_cff_pcp_nrm[idx]=scv_cff[idx]/flx_mss_pcp_rsl; // [m2 kg-1] Scavenging coefficient, precipitation-normalized SeP97 p. 1021 (20.53, 20.58)
    scv_cff_mss_avg+=scv_cff[idx]*cnc[idx]*mss[idx]; // [s-1] Mass mean scavenging coefficient of aerosol SeP97 p. 1022 (20.63)
    scv_cff_nbr_avg+=scv_cff[idx]*cnc[idx]; // [s-1] Number mean scavenging coefficient of aerosol SeP97 p. 1022 (20.68)
  } // end loop over sz
  scv_cff_mss_avg/=mss_rsl; // [s-1] Mass mean scavenging scavenging coefficient
  scv_cff_mss_avg_pcp_nrm=scv_cff_mss_avg/flx_mss_pcp_rsl; // [m2 kg-1] Mass mean scavenging scavenging coefficient, precipitation-normalized
  scv_cff_nbr_avg/=cnc_nbr_rsl; // [s-1] Number mean scavenging scavenging coefficient
  scv_cff_nbr_avg_pcp_nrm=scv_cff_nbr_avg/flx_mss_pcp_rsl; // [m2 kg-1] Number mean scavenging scavenging coefficient, precipitation-normalized

  /*
  // Approximate wet deposition
  prc_cmp scv_cff_mss_avg_apx/=mss_rsl; // [s-1] Approximate mass mean scavenging scavenging coefficient
  prc_cmp scv_cff_mss_avg_pcp_nrm_apx; // [m2 kg-1] Approximate mass mean scavenging scavenging coefficient, precipitation-normalized
  scv_cff_mss_avg_apx=scv_cff_apx[idx]*cnc[idx]*mss[idx]; // [s-1] Mass mean scavenging coefficient of aerosol SeP97 p. 1022 (20.63)
  scv_cff_mss_avg_pcp_nrm_apx=scv_cff_mss_avg_apx/flx_mss_pcp_rsl; // [m2 kg-1] Approximate mass mean scavenging scavenging coefficient, precipitation-normalized
  */

  // Partition vertical dust flux into size bins
  
  // Reynolds numbers
  prc_cmp *ryn_nbr=new prc_cmp[sz_nbr]; // [frc] Reynolds number ambient windspeed
  prc_cmp *ryn_nbr_frc=new prc_cmp[sz_nbr]; // [frc] Friction Reynolds number
  // Reynolds numbers for ambient windspeed and friction speed
  for(idx=0;idx<sz_nbr;idx++){
    ryn_nbr[idx]=wnd_mdp*dmt_ctr[idx]/vsc_knm_atm; // [frc] SeP97 p. 460
    ryn_nbr_frc[idx]=wnd_frc_mbl*dmt_ctr[idx]/vsc_knm_atm; // [frc] SeP97 p. 460
  } // end loop over sz

  // Threshold friction speed (iterative computation)
  prc_cmp *ryn_nbr_frc_thr=new prc_cmp[sz_nbr]; // [frc] Threshold friction Reynolds number
  prc_cmp *ryn_nbr_frc_thr_prx=new prc_cmp[sz_nbr]; // [frc] Threshold friction Reynolds number approximation
  prc_cmp *wnd_frc_thr=new prc_cmp[sz_nbr]; // [m s-1] Threshold friction speed
  prc_cmp *wnd_frc_thr_prx=new prc_cmp[sz_nbr]; // [m s-1] Threshold friction speed approximation
  rcd+=wnd_frc_thr_get
    (sz_nbr, // I [nbr] Size of arrays
     dmt_ctr, // I [m] Diameter at bin center
     dns_prt, // I [kg m-3] Density of particle
     dns_mdp, // I [kg m-3] Midlayer density
     vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
     ryn_nbr_frc_thr, // O [frc] Threshold friction Reynolds number
     ryn_nbr_frc_thr_prx, // O [frc] Threshold friction Reynolds number approximation
     wnd_frc_thr, // O [m s-1] Threshold friction speed
     wnd_frc_thr_prx); // O [m s-1] Threshold friction speed approximation

  // Optimal threshold friction speed
  prc_cmp wnd_frc_thr_opt; // [m s-1] Optimal threshold friction speed for saltation
  prc_cmp dmt_slt_opt; // [m] Optimal diameter for saltation
  wnd_frc_thr_opt=vec_min(wnd_frc_thr,sz_nbr); // [m s-1]
  dmt_slt_opt=dmt_ctr[vec_val2idx(wnd_frc_thr,sz_nbr,wnd_frc_thr_opt)]; // [m] Locate bin with minima first

  /* Heterogeneous chemistry
     Validation #1: HNO3 chemistry from ZSK94 PEM case
     ZSK94 p. 816 Table 1 size bin 3, in particular
     Diameter = 2.0 um, T = 283 K, alpha = 0.1, 
     z = 4 km --> p ~ 700 mb 
     HNO3 vmr = 0.1 pptv
     Total dust mass concentration = 100.0 ug m-3, bin 3 is 26% of mass = 26 ug m-3
     26.0e-9 kg m-3 x 1.11994e+14 # kg-1 = 2.91e6 # m-3
     These conditions are implemented with the following command:
     mie -dbg -no_mie --sz_mnm=1.0 --sz_mxm=1.01 --sz_dbg_mcr=1.0 --dns_prt=2100.0 --tpt_mdp=283.0 --prs_mdp=70000.0 --gsd_anl=1.5 --cnc_nbr_anl=2.91e6
     Validation #2: HNO3 chemistry from ZhC99 p. 359 PEM-West-A case
     HNO3 vmr = 0.05 ppbv
     Total dust mass concentration = 220.0 ug m-3 or 100.0 ug m-3
     220.0e-9 kg m-3 x 3.80637e+13 # kg-1 = 8.36e6 # m-3 or 3.80e6 # m-3
     These conditions are implemented with the following command:
     mie -dbg -no_mie --sz_mnm=0.1 --sz_mxm=40.0 --sz_nbr=100 --sz_dbg_mcr=0.5 --dns_prt=2600.0 --tpt_mdp=300.0 --prs_mdp=101300.0 --rds_nma=0.88 --gsd_anl=1.7 --cnc_nbr_anl=3.80e6
     19990811: Current model predictions of dff_HNO3_aer and rxrc_HNO3_aer are ~3x too large and too small, respectively, relative to ZhC99 p. 360 Fgr. 3. This is explained by: Binary diffusivity of HNO3 in air is factor of 1.8 larger than ZhC99 assume
     One problem may be ambiguous definitions of continuum regime diffusion correction, cnt_rgm_bll_crc_HNO3_aer */
  
  // Chemical constants (perhaps should be in chm.hh)
  using phc::mmw_HNO3; // (6.2995644e-02) [kg mol-1] Mean molecular weight of HNO3 HITRAN96
  const prc_cmp mss_rat_HNO3(mmw_HNO3/mmw_dry_air); // [frc] Mass ratio of HNO3 to dry air SeP97 p. 457 
  const prc_cmp mss_acm_cff_HNO3_aer(0.01); // [frc] Mass accomodation coefficient of HNO3 to aerosol ZhC99 p. 357 Tbl. 2, ZSK94 p. 817
  const prc_cmp mss_upt_cff_HNO3_aer(mss_acm_cff_HNO3_aer); // [frc] Mass uptake coefficient of HNO3 to aerosol SeP97 p. 633 (11.121)
 
  // Environmental chemical properties
  prc_cmp cff_hnr_HNO3_H2O; // [mol ltr-1 atm-1] Henry's Law coefficient for HNO3(g)<-->H2O(l) equilibrium
  prc_cmp cnc_dry_air(prs_mdp/(cst_Boltzmann*tpt_mdp)); // [mlc m-3] Concentration of dry air
  prc_cmp cnc_HNO3_gas(vmr_HNO3_gas*cnc_dry_air); // [mlc m-3] Concentration of HNO3
  prc_cmp dff_HNO3_air; // [mlc m2 s-1] Binary diffusivity of HNO3 in air
  prc_cmp dmt_cll_avg; // [m] Mean collision diameter
  prc_cmp mfp_HNO3_air; // [m] Mean free path of HNO3 in air
  prc_cmp vlc_mwb_HNO3; // [m s-1] Thermal speed of HNO3
  using phc::cst_Avagadro; // (6.022045e+23) [mlc mol-1] Avagadro's number
  prc_cmp q_HNO3_gas(cnc_HNO3_gas*mmw_HNO3/(dns_mdp*cst_Avagadro)); // [kg kg-1] Mixing ratio of HNO3
  // Temperature dependence of Henry's Law coefficient for HNO3
  using phc::cff_hnr_HNO3_H2O_298K; // (2.1e5) [mol ltr-1 atm-1] Henry's Law coefficient of HNO3 in liquid water at 298K SeP97 p. 341 Table 6.2
  using phc::rxn_ntp_HNO3_H2O_298K; // [J mol-1] Reaction enthalpy (heat of dissolution) at 298K inferred from E/R quoted by Sep97 p. 391 Tbl. 6.A.1
  using phc::tpt_298; // (298.00) [K] Standard temperature for Henry's Law
  cff_hnr_HNO3_H2O=cff_hnr_HNO3_H2O_298K*std::exp(rxn_ntp_HNO3_H2O_298K*(1.0/tpt_298-1.0/tpt_mdp)/gas_cst_unv); // [mol ltr-1 atm-1] SeP97 p. 342 (6.5)
  vlc_mwb_HNO3=std::sqrt(8.0*gas_cst_unv*tpt_mdp/(mth::cst_M_PIl*mmw_HNO3)); // [m s-1] SeP97 p. 453 (8.2)
  using phc::dmt_cll_HNO3; // (3.5e-10) [m] fxm: Pure guess
  using phc::dmt_cll_air; // (3.46e-10) [m] Mean collision diameter of air SeP97 p. 1292 Table A.7
  dmt_cll_avg=0.5*(dmt_cll_HNO3+dmt_cll_air); // [m] SeP97 p. 457 (8.10)
  mfp_HNO3_air=1.0/(mth::cst_M_PIl*std::sqrt(2.0)*cnc_HNO3_gas*dmt_cll_HNO3*dmt_cll_HNO3+mth::cst_M_PIl*std::sqrt(1.0+mss_rat_HNO3)*cnc_dry_air*dmt_cll_avg*dmt_cll_avg); // [m] SeP97 p. 457 (8.9)
  //  mfp_HNO3_air=1.0/(mth::cst_M_PIl*std::sqrt(1.0+mss_rat_HNO3)*cnc_dry_air*dmt_cll_avg*dmt_cll_avg); // [m] SeP97 p. 457 (8.11)
  //  mfp_HNO3_air=3.0*dff_HNO3_air/vlc_mwb_HNO3; // [m] SeP97 p. 604 Tbl. 11.1
  dff_HNO3_air=vlc_mwb_HNO3*3.0*mth::cst_M_PIl*(1.0+mss_rat_HNO3)*mfp_HNO3_air/32.0; // [m2 s-1] SeP97 p. 457 (8.13)

  // Size-dependent mass transfer properties
  prc_cmp *dff_HNO3_aer=new prc_cmp[sz_nbr]; // [m3 s-1 prt-1] Normalized diffusion rate of HNO3 to aerosol
  prc_cmp *dff_HNO3_aer_cnt=new prc_cmp[sz_nbr]; // [m3 s-1 prt-1] Normalized diffusion rate of HNO3 to aerosol
  prc_cmp *dff_HNO3_aer_knt=new prc_cmp[sz_nbr]; // [m3 s-1 prt-1] Normalized diffusion rate of HNO3 to aerosol
  prc_cmp *vnt_crc_aer=new prc_cmp[sz_nbr]; // [frc] Ventilation correction factor for gas phase diffusion to aerosol surface
  prc_cmp *cnt_rgm_bll_crc_HNO3_aer=new prc_cmp[sz_nbr]; // [frc] Correction to continuum regime diffusion approximation due to ballistic effects
  prc_cmp *knd_nbr_HNO3_air=new prc_cmp[sz_nbr]; // [frc] Knudsen number of HNO3 in air
  prc_cmp *knt_rgm_dff_crc_HNO3_aer=new prc_cmp[sz_nbr]; // [frc] Correction to kinetic regime ballistic approximation due to diffusion effects
  prc_cmp *cnt_knt_rgm_rat_HNO3_aer=new prc_cmp[sz_nbr]; // [frc] Ratio of continuum to kinetic regime solutions
  for(idx=0;idx<sz_nbr;idx++){
    // Knudsen number
    knd_nbr_HNO3_air[idx]=2.0*mfp_HNO3_air/dmt_ctr[idx]; // [frc] SeP97 p. 456 (8.8)

    // Ventilation correction is negligible for D < 50 um
    if(ryn_nbr_grv[idx] <= 2.5) vnt_crc_aer[idx]=1.0+0.09635*ryn_nbr_grv[idx]; else vnt_crc_aer[idx]=0.78+0.27477*std::sqrt(ryn_nbr_grv[idx]); // [frc] Ventilation correction factor for gas phase diffusion to aerosol surface RoY94 p. 116, PrK78 p. 443 (13-57), PrK98 p. 541 (13-60)
    /* The kinetic uptake correction factor is quite significant for D < 50 um: 
       cnt_rgm_bll_crc_HNO3_aer = 0.04 and 0.40 for D = 1 and 10 um, respectively
       Unfortunately, the correction takes many forms, all complex algebraic expressions
       SeP97 p. 604 notes that all forms give comparable results for D > 0.2 um
       However, it is important to use the expression for the mean free path that
       is consistent with the derivation of the diffusion correction factor.
       SeP97 p. 604 Table 11.1 lists all the correction factors and
       their associated definition of the mean free path.
       The expression of Fuchs and Sutugin (1971) (SeP97 p. 606 (11.43)), even
       though is was assumed the molecule was much lighter than air in deriving it.
    */
    cnt_rgm_bll_crc_HNO3_aer[idx]=0.75*mss_acm_cff_HNO3_aer*(1.0+knd_nbr_HNO3_air[idx])/(knd_nbr_HNO3_air[idx]*knd_nbr_HNO3_air[idx]+knd_nbr_HNO3_air[idx]+0.283*knd_nbr_HNO3_air[idx]*mss_acm_cff_HNO3_aer+0.75*mss_acm_cff_HNO3_aer); // [frc] SeP97 p. 606 (11.43)
    // cnt_rgm_bll_crc_HNO3_aer[idx]=1.0/(1.0+knd_nbr_HNO3_air[idx]*(mfp_HNO3_air+4.0*(1.0-mss_acm_cff_HNO3_aer)/(3.0*mss_acm_cff_HNO3_aer))); // [frc] ZSK94 p. 817 (3), ZhC99 p. 356 (9)
    /* The transport of HNO3 to the aerosol surface is 
       Jc = 4 * pi * R * Dab * fv * F SeP97 p. 636 (11.128), ZSK94 p. 817 (2)
       where the RHS is the solution to the diffusion equation (valid for Kn << 1)
       modified by the ventilation factor fv and the correction factor F to the 
       diffusion equation caused by ballistic effects (important at Kn ~ 1) and by
       uptake effects (i.e., alpha != 1)
    */

    /* Continuum regime approximation
       Jc = 4 * pi * radius * (uptake coefficient ?) * binary diffusivity
       19990813: fxm: Should Jc include alpha as Jk does? I think so, does SeP97?
       Doing so makes sense to me but brings results far away from ZhC99 */
    dff_HNO3_aer_cnt[idx]=4.0*mth::cst_M_PIl*sz_ctr[idx]*dff_HNO3_air; // [m3 s-1 prt-1] SeP97 p. 636 (11.128), ZSK94 p. 817 (2)

    // Kinetic regime approximation
    // Jk = ( Aerosol surface area * uptake coefficient * molecular velocity ) / 4.0
    dff_HNO3_aer_knt[idx]=4.0*mth::cst_M_PIl*sz_ctr[idx]*sz_ctr[idx]*mss_upt_cff_HNO3_aer*vlc_mwb_HNO3/4.0; // [m3 s-1 prt-1] SeP97 p. 600 (11.23) p. 618 (11.80) BOT99 p. 118 (3.40)

    cnt_knt_rgm_rat_HNO3_aer[idx]=dff_HNO3_aer_cnt[idx]/dff_HNO3_aer_knt[idx]; // [frc] Ratio of continuum to kinetic regime solutions

    /* Transition regime approximation obeys correct limits as Kn-->0 and Kn-->infty
       Implement full solution as correction to the continuum regime since 
       tropospheric heterogeneous chemistry is closer to the continuum regime. 
       Also apply ventilation correction here for simplicity. */
    dff_HNO3_aer[idx]=dff_HNO3_aer_cnt[idx]*vnt_crc_aer[idx]*cnt_rgm_bll_crc_HNO3_aer[idx]; // [m3 s-1 prt-1] SeP97 p. 636 (11.128), ZSK94 p. 817 (2)

    knt_rgm_dff_crc_HNO3_aer[idx]=dff_HNO3_aer[idx]/dff_HNO3_aer_knt[idx]; // [frc]

  } // end loop over sz
  // Free un-needed memory from recent calculations
  delete []cnt_knt_rgm_rat_HNO3_aer; // [frc] Ratio of continuum to kinetic regime solutions

  // Heterogeneous removal rates and mixing ratio adjustment
  prc_cmp *rxr_HNO3_gas_aer_vmr=new prc_cmp[sz_nbr]; // [mlc mlc-1 s-1] Mean rate of HNO3 removal by aerosol
  prc_cmp *rxr_nst_HNO3_gas_aer_vmr=new prc_cmp[sz_nbr]; // [mlc mlc-1 s-1] Instantaneous rate of HNO3 removal by aerosol
  prc_cmp *rxrc_HNO3_aer=new prc_cmp[sz_nbr]; // [s-1] Pseudo first order rate coefficient for HNO3 removal by aerosol
  prc_cmp rxr_HNO3_gas_aer_vmr_ttl(0.0); // [mlc mlc-1 s-1] Total mean rate of HNO3 removal by aerosol
  prc_cmp rxr_nst_HNO3_gas_aer_vmr_ttl(0.0); // [mlc mlc-1 s-1] Total instantaneous rate of HNO3 removal by aerosol
  prc_cmp rxrc_HNO3_aer_ttl(0.0); // [s-1] Total pseudo first order rate coefficient for HNO3 removal by aerosol
  prc_cmp vmr_HNO3_gas_dlt; // [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio
  prc_cmp *rxrc_HNO3_aer_knt=new prc_cmp[sz_nbr]; // [s-1] Pseudo first order rate coefficient for HNO3 removal by aerosol, kinetic approximation
  prc_cmp rxrc_HNO3_aer_knt_ttl(0.0); // [s-1] Total pseudo first order rate coefficient for HNO3 removal by aerosol, kinetic approximation
  for(idx=0;idx<sz_nbr;idx++){
    // Pseudo first order rate coefficients kinetic approximation
    rxrc_HNO3_aer_knt[idx]=cnc[idx]*dff_HNO3_aer_knt[idx]; // [s-1] ZSK94 p. 817 (1), SeP97 p. 636 (11.136) ZhC99 p. 356 (9)
    rxrc_HNO3_aer_knt_ttl+=rxrc_HNO3_aer_knt[idx]; // [s-1] ZSK94 p. 817 (1), SeP97 p. 636 (11.136) ZhC99 p. 356 (9)

    // Pseudo first order rate coefficients
    rxrc_HNO3_aer[idx]=cnc[idx]*dff_HNO3_aer[idx]; // [s-1] ZSK94 p. 817 (1), SeP97 p. 636 (11.136) ZhC99 p. 356 (9)
    rxrc_HNO3_aer_ttl+=rxrc_HNO3_aer[idx]; // [s-1] ZSK94 p. 817 (1), SeP97 p. 636 (11.136) ZhC99 p. 356 (9)

    // Sink
    vmr_HNO3_gas_dlt=vmr_HNO3_gas*(1.0-std::exp(-rxrc_HNO3_aer[idx]*tm_dlt)); // [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio

    // Mean reaction rate over a timestep in volume mixing ratio units
    rxr_HNO3_gas_aer_vmr[idx]=vmr_HNO3_gas_dlt/tm_dlt; // [mlc mlc-1 s-1]
    rxr_HNO3_gas_aer_vmr_ttl+=rxr_HNO3_gas_aer_vmr[idx]; // [mlc mlc-1 s-1]

    // Instantaneous reaction rates in volume mixing ratio units
    rxr_nst_HNO3_gas_aer_vmr[idx]=vmr_HNO3_gas*rxrc_HNO3_aer[idx]; // [mlc mlc-1 s-1]
    rxr_nst_HNO3_gas_aer_vmr_ttl+=rxr_nst_HNO3_gas_aer_vmr[idx]; // [mlc mlc-1 s-1]

    // Adjust concentration
    vmr_HNO3_gas=vmr_HNO3_gas-vmr_HNO3_gas_dlt; // [mlc mlc-1] Gaseous HNO3 volume mixing ratio
    cnc_HNO3_gas=vmr_HNO3_gas*cnc_dry_air; // [mlc m-3] Concentration of HNO3
    q_HNO3_gas=cnc_HNO3_gas*mmw_HNO3/(dns_mdp*cst_Avagadro); // [kg kg-1] Mixing ratio of HNO3
  } // end loop over s
  // Free un-needed memory from recent calculations
  delete []rxr_nst_HNO3_gas_aer_vmr; // [mlc mlc-1 s-1] Instantaneous rate of HNO3 removal by aerosol

  /* Timescales for aqueous chemistry
     Species-specific aqueous diffusion rates are hard to find
     PrK78 say 1.0e-9 < Dl < 2.0e-9 m2 s-1, and use 1.15e-9 m2 s-1 for S(IV) */
  const prc_cmp dff_HNO3_H2O_lqd(1.0e-9); // [m2 s-1] Diffusion coefficient of HNO3 in liquid H2O SeP97 p. 616 PrK98 p. 764
  prc_cmp tau_rxrc_HNO3_aer_ttl(1.0/rxrc_HNO3_aer_ttl); // [s] Total timescale for HNO3 uptake by aerosol
  prc_cmp tau_rxrc_HNO3_aer_knt_ttl(1.0/rxrc_HNO3_aer_knt_ttl); // [s] Total timescale for HNO3 uptake by aerosol, kinetic approximation
  prc_cmp *tau_rxrc_HNO3_aer=new prc_cmp[sz_nbr]; // [s] Timescale for HNO3 uptake by aerosol NB: variable must usually be stored as NC_DOUBLE
  prc_cmp *tau_rxrc_HNO3_aer_knt=new prc_cmp[sz_nbr]; // [s] Timescale for HNO3 uptake by aerosol, kinetic approximation
  prc_cmp *tau_gpd_HNO3_aer=new prc_cmp[sz_nbr]; // [s] Timescale for HNO3 gas-phase diffusion to aerosol
  prc_cmp *tau_ntf_eqm_HNO3_aer=new prc_cmp[sz_nbr]; // [s] Timescale for HNO3 to achieve phase equilibrium at interface
  prc_cmp *tau_aqs_dss_HNO3_aer=new prc_cmp[sz_nbr]; // [s] Timescale for aqueous dissociation of HNO3 to NO3-
  prc_cmp *tau_aqs_dff_HNO3_aer=new prc_cmp[sz_nbr]; // [s] Timescale for aqueous phase diffusion of HNO3 in H2O droplet
  for(idx=0;idx<sz_nbr;idx++){
    // e-folding time of HNO3 due to irreversible uptake by current size bin
    tau_rxrc_HNO3_aer[idx]=1.0/rxrc_HNO3_aer[idx]; // [s] Timescale for HNO3 uptake by aerosol
    tau_rxrc_HNO3_aer_knt[idx]=1.0/rxrc_HNO3_aer_knt[idx]; // [s] Timescale for HNO3 uptake by aerosol, kinetic approximation
    if(dbg_lvl == dbg_old) std::cout << "tau_rxrc_HNO3_aer[" << idx << "] = " << tau_rxrc_HNO3_aer[idx] << std::endl;
    // Timescale for gaseous diffusion-limited HNO3 uptake to bring a droplet to Henry's law equilibrium?
    tau_gpd_HNO3_aer[idx]=0.25*sz_ctr[idx]*sz_ctr[idx]/dff_HNO3_air; // [s] SeP97 p. 611 (11.49) PrK98 p. 759 
    // Timescale to reach balance between uptake and evaporation of HNO3 at surface of droplet
    // Highly soluble gases (like HNO3):
    tau_ntf_eqm_HNO3_aer[idx]=sz_ctr[idx]*cff_hnr_HNO3_H2O*std::sqrt(2.0*mth::cst_M_PIl*mmw_HNO3*gas_cst_unv*tpt_mdp)/(3.0*mss_acm_cff_HNO3_aer); // [s] SeP97 p. 613 (11.61)
    // Insoluble gases (like O3):
    // tau_ntf_eqm_HNO3_aer[idx]=sz_ctr[idx]*sz_ctr[idx]/(mth::cst_M_PIl*mth::cst_M_PIl*dff_HNO3_H2O_lqd); // [s] SeP97 p. 613 (11.62)
    // Timescale for hydrolyzed HNO3 to dissociate to NO3- + H+
    tau_aqs_dss_HNO3_aer[idx]=-1.0; // [s] fxm: Timescale for aqueous dissociation of HNO3 to NO3- + H+ SeP97 p. 615 (11.71) 
    // Timescale for hydrolyzed HNO3 to diffuse from surface to center of droplet
    tau_aqs_dff_HNO3_aer[idx]=sz_ctr[idx]*sz_ctr[idx]/(mth::cst_M_PIl*mth::cst_M_PIl*dff_HNO3_H2O_lqd); // [s] Timescale for aqueous phase diffusion of HNO3 in H2O droplet SeP97 p. 616 (11.75)
  } // end loop over sz

  // Hygroscopic growth measurements (increasing RH) from Han76 p. 115, Table IV
  const prc_cmp hyg_dst_RH[]={
    0.204, 0.349, 0.457, 0.590, 0.648,
    0.700, 0.751, 0.789, 0.843, 0.896,
    0.900, 0.955, 0.971, 0.976, 0.986,
    0.990, 0.997 
  };
  const prc_cmp hyg_dst_mlmic[]={
    0.032, 0.024, 0.021, 0.024, 0.033,
    0.088, 0.130, 0.133, 0.120, 0.107,
    0.105, 0.075, 0.072, 0.071, 0.067, 
    0.065, 0.065
  };
  /*  const prc_cmp hyg_dst_mwrmd[]={
    0.008, 0.013, 0.018, 0.034, 0.060, 
    0.206, 0.392, 0.497, 0.644, 0.922,
    0.945, 1.604, 2.411, 2.887, 4.719,
    6.435, 20.90
    }; */
  long RH_nbr=(sizeof hyg_dst_RH)/sizeof(prc_cmp); // Length of hygroscopic data vectors
  prc_cmp mlmic; // [frc] Mean linear mass increase coefficient
  mlmic=ntp_vec_one(RH_nbr,hyg_dst_RH,hyg_dst_mlmic,RH_lqd); // [frc]
  prc_cmp *mss_wtr=new prc_cmp [sz_nbr]; // [kg] Mass of water on aerosol
  prc_cmp *svp_fct_klv=new prc_cmp [sz_nbr]; // [frc] Saturation vapor pressure Kelvin factor
  prc_cmp *svp_fct_slt=new prc_cmp [sz_nbr]; // [frc] Saturation vapor pressure solute factor
  prc_cmp *svp_fct_ttl=new prc_cmp [sz_nbr]; // [frc] Saturation vapor pressure total factor
  prc_cmp *mss_wtr_rat=new prc_cmp [sz_nbr]; // [frc] Ratio of water mass to dry aerosol mass
  prc_cmp *RH_lqd_sfc=new prc_cmp [sz_nbr]; // [frc] Relative humidity w/r/t liquid water at particle surface
  prc_cmp *act_cff=new prc_cmp [sz_nbr]; // [frc] Activity coefficient of solution
  prc_cmp *mss_frc_solute=new prc_cmp [sz_nbr]; // [frc] Mass fraction of solute
  for(idx=0;idx<sz_nbr;idx++){
    mss_frc_solute[idx]=0.5; // [frc] Mass fraction of solute
  } // end loop over sz
  // [frc] Activity coefficient
  rcd+=act_cff_TaM94(sz_nbr,mss_frc_solute,act_cff);
  prc_cmp sfc_tns_wtr_lqd; // [J m-2] Surface tension of liquid water
  rcd+=sfc_tns_wtr_lqd_PrK78(1,&tpt_mdp,&sfc_tns_wtr_lqd);
  using phc::mmw_H2O; // (1.8015259e-02) [kg mol-1] Mean molecular weight of H2O HITRAN96
  prc_cmp svp_fct_klv_xpn; // [frc] Exponent in Kelvin effect IrG81 p. 89 (5) 
  for(idx=0;idx<sz_nbr;idx++){
    // Currently we neglect solute effect but account for curvature effect
    svp_fct_klv_xpn=2.0*mmw_H2O*sfc_tns_wtr_lqd/(dns_H2O_lqd_std*gas_cst_unv*tpt_mdp*sz_ctr[idx]); // [frc] Kelvin effect IrG81 p. 89 (5) 
    svp_fct_klv[idx]=std::exp(svp_fct_klv_xpn); // [frc] Kelvin effect IrG81 p. 89 (5) 
    if(svp_fct_klv_xpn > 1.0) wrn_prn(sbr_nm,"Kelvin (curvature) effect problem: svp_fct_klv_xpn = "+nbr2sng(svp_fct_klv_xpn)+" > 1.0. Reducing Kelvin effect to unity. (HINT: aerosol may be too small to apply classical Kelvin theory?).");
    svp_fct_klv[idx]=min(1.0,svp_fct_klv[idx]); // [frc] Kelvin effect IrG81 p. 89 (5) 
    svp_fct_slt[idx]=1.0; // [frc] Solute effect
    svp_fct_ttl[idx]=svp_fct_klv[idx]*svp_fct_slt[idx]; // [frc] Total effect
    RH_lqd_sfc[idx]=RH_lqd*svp_fct_ttl[idx]; // [frc]
    assert(RH_lqd_sfc[idx] <= 1.0);
    mss_wtr[idx]=vlm[idx]*dns_prt*mlmic*RH_lqd_sfc[idx]/(1.0-RH_lqd_sfc[idx]); // [kg]
    mss_wtr_rat[idx]=mss_wtr[idx]/mss[idx]; // [frc]
  } // end loop over sz
  // Free un-needed memory from recent calculations
  delete []mss_wtr; // [kg] Mass of water on aerosol
  delete []RH_lqd_sfc; // [frc] Relative humidity w/r/t liquid water at particle surface
  
  // Wind statistics
  prc_cmp *pugtut=new prc_cmp[sz_nbr]; // [frc] Probability windspeed exceeds threshold
  prc_cmp *wnd_rfr_thr=new prc_cmp[sz_nbr]; // [m s-1] Threshold 10 m wind speed
  prc_cmp str_shr=dns_mdp*wnd_rfr_mbl/rss_aer_mmn; // [kg m-1 s-2] Shear stress Bon96 p. 47
  bool JHM78_20(true); // Use JHM78 (20) approximation
  if(JHM78_20) wbl_shp=0.94*std::sqrt(wnd_rfr_mbl); // [frc] JHM78 (20)
  if(std::fabs(std::sqrt(mth::cst_M_PIl)-std::exp(gsl_sf_lngamma(0.5)))/std::sqrt(mth::cst_M_PIl) > 1.0e-5) err_prn(prg_nm,sbr_nm,"Gamma function error"); // Sanity check
  prc_cmp wbl_scl=wnd_rfr_mbl/std::exp(gsl_sf_lngamma((wbl_shp+1.0)/wbl_shp)); // [m s-1] JHM78 (16)
  for(idx=0;idx<sz_nbr;idx++){
    wnd_rfr_thr[idx]=wnd_frc_thr[idx]*wnd_frc_thr[idx]*rss_aer_mmn/dns_mdp; // [m s-1] 
    pugtut[idx]=std::exp(-std::pow((wnd_rfr_thr[idx]/wbl_scl),wbl_shp)); // [frc] JHM78 (2), SRL96
  } // end loop over sz
  
  // Instantiate wavelength grid
  wvl_grd_cls wvlgrd(wvl_grd_sng,wvl_mnm,wvl_mxm,wvl_nbr); // [obj] Wavelength grid
  //  if(dbg_lvl >= dbg_fl) std::cout << "quark1 wvl_mnm = " << wvl_mnm << ", wvl_mxm = " << wvl_mxm << std::endl;
  //if(dbg_lvl >= dbg_fl) std::cout << "quark2 Wavelength grid = " << wvlgrd << std::endl;
  // Update fields which may be overridden by constructors
  wvl_nbr=wvlgrd.wvl_nbr_get(); // [m] Number of wavelength bands
  wvl_dlt_mcr=1.0e6*(wvl_mxm-wvl_mnm); // [um] Bandwidth
  wvl_mdp_mcr=0.5*1.0e6*(wvl_mnm+wvl_mxm); // [um] Midpoint wavelength
  wvl_mnm_mcr=1.0e6*wvl_mnm; // [um] Minimum wavelength
  wvl_mxm_mcr=1.0e6*wvl_mxm; // [um] Maximum wavelength
  wvn_mnm_xcm=0.01/wvl_mxm; // [cm-1] Minimum wavenumber
  wvn_mxm_xcm=0.01/wvl_mnm; // [cm-1] Maximum wavenumber
  // Gather public pointers to data
  const prc_cmp *wvl_ctr=wvlgrd.wvl_ctr_get(); // [m] Wavelength at band center
  const prc_cmp *wvl=wvl_ctr; // [m] Nominal wavelength
  const prc_cmp *wvl_dlt=wvlgrd.wvl_dlt_get(); // [m] Bandwidth
  const prc_cmp *wvl_grd=wvlgrd.wvl_grd_get(); // [m] Wavelength grid
  const prc_cmp *wvl_max=wvlgrd.wvl_max_get(); // [m] Maximum wavelength in band
  const prc_cmp *wvl_min=wvlgrd.wvl_min_get(); // [m] Minimum wavelength in band
  const prc_cmp *wvn_ctr=wvlgrd.wvn_ctr_get(); // [cm-1] Wavenumber at band center
  const prc_cmp *wvn=wvn_ctr; // [cm-1] Nominal wavenumber
  const prc_cmp *wvn_dlt=wvlgrd.wvn_dlt_get(); // [cm-1] Bandwidth
  const prc_cmp *wvn_grd=wvlgrd.wvn_grd_get(); // [cm-1] Wavenumber grid
  const prc_cmp *wvn_max=wvlgrd.wvn_max_get(); // [cm-1] Maximum wavenumber in band
  const prc_cmp *wvn_min=wvlgrd.wvn_min_get(); // [cm-1] Minimum wavenumber in band
  // Derive debugging information
  // 20210418: Explanation for non-const wvl_idx_dbg is near sz_idx_dbg declaration
  // const long wvl_idx_dbg=vec_val2idx(wvl_ctr,wvl_nbr,wvl_dbg); // [idx] Debugging wavelength bin
  long wvl_idx_dbg=vec_val2idx(wvl_ctr,wvl_nbr,wvl_dbg); // [idx] Debugging wavelength bin

  // Instantiate spectral solar irradiance sources
  spc_slr_cls flx_slr_src(slr_spc_key,slr_cst,fl_slr_spc); // [sct] Solar flux source
  fl_slr_spc=flx_slr_src.fl_slr_spc_get(); // [sng] File containing solar spectrum
  spc_slr_cls *flx_slr_src_ThD71=new spc_slr_cls("ThD71"); // [sct] Solar flux source ThD71
  spc_slr_cls *flx_slr_src_LaN68=new spc_slr_cls("LaN68"); // [sct] Solar flux source LaN68
  prc_cmp *flx_slr_frc_ThD71=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux in band ThD71
  prc_cmp *flx_slr_frc_LaN68=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux in band LaN68
  prc_cmp *flx_slr_frc=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux in band
  rcd+=flx_slr_src_ThD71->flx_frc_get(wvl_min,wvl_max,wvl_nbr,flx_slr_frc_ThD71); // [frc] Fraction of solar flux in band ThD71
  rcd+=flx_slr_src_LaN68->flx_frc_get(wvl_min,wvl_max,wvl_nbr,flx_slr_frc_LaN68); // [frc] Fraction of solar flux in band LaN68
  rcd+=flx_slr_src.flx_frc_get(wvl_min,wvl_max,wvl_nbr,flx_slr_frc); // [frc] Fraction of solar flux in band
  // These two spc_slr_cls objects were needed only for diagnostic output fluxes
  delete flx_slr_src_LaN68; // [sct] Solar flux source LaN68
  delete flx_slr_src_ThD71; // [sct] Solar flux source ThD71

  // Infrared spectral irrandiance
  spc_bbd_cls spc_bbd(tpt_bbd_wgt); // [bbd] Blackbody spectrum object
  prc_cmp *flx_IR_frc=new prc_cmp[wvl_nbr]; // [frc] Fraction of infrared flux in band
  rcd+=spc_bbd.flx_frc_get(wvl_min,wvl_max,wvl_nbr,flx_IR_frc); // [frc] Fraction of infrared flux in band

  // Irradiance diagnostics for solar and infrared wavelengths
  prc_cmp *flx_slr=new prc_cmp[wvl_nbr]; // [W m-2] Solar flux in band
  prc_cmp *flx_slr_frc_blr=new prc_cmp[wvl_nbr]; // [frc] Fraction of solar flux at shorter wavelengths
  prc_cmp *flx_spc_slr=new prc_cmp[wvl_nbr]; // [W m-2 m-1] Solar spectral flux in band
  prc_cmp *flx_spc_slr_pht=new prc_cmp[wvl_nbr]; // [pht m-2 s-1 m-1] Solar spectral photon flux in band
  prc_cmp *nrg_pht=new prc_cmp[wvl_nbr]; // [J pht-1] Energy of photon at band center
  using phc::cst_Planck; // (6.62620e-34) [J s] Planck's constant
  using phc::speed_of_light; // (2.99793e+08) [m s-1] Speed of light in vacuo 
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    flx_slr[wvl_idx]=slr_cst*flx_slr_frc[wvl_idx]; // [W m-2] Solar flux in band
    flx_spc_slr[wvl_idx]=flx_slr[wvl_idx]/wvl_dlt[wvl_idx]; // [W m-2 m-1] Solar spectral flux in band
    nrg_pht[wvl_idx]=cst_Planck*speed_of_light/wvl_ctr[wvl_idx]; // [J pht-1] Energy of photon at band center
    flx_spc_slr_pht[wvl_idx]=flx_spc_slr[wvl_idx]/nrg_pht[wvl_idx]; // [pht m-2 s-1 m-1] Solar spectral photon flux in band
  } // end loop over wvl 
  if(wvl_nbr > 1 && mnt_chk(wvl,wvl_nbr)){
    // Generate series of partial sums series with fraction of flux bluer
    // NB: Technique is ill-defined for non-monotonic wavelength grids
    const bool flg_ncr_wvl(wvl[1] > wvl[0] ? true : false); // [flg] Wavelength grid increases
    long idx_max_wvl; // [idx] Index of maximum wavelenth
    long idx_min_wvl; // [idx] Index of minimum wavelenth
    if(flg_ncr_wvl){
      idx_max_wvl=wvl_nbr-1L;
      idx_min_wvl=0L;
    }else{
      idx_max_wvl=0L;
      idx_min_wvl=wvl_nbr-1L;
    } // endif flg_ncr_wvl
    if(flg_ncr_wvl){
      flx_slr_frc_blr[idx_min_wvl]=0.0; // [frc] Fraction of solar flux at shorter wavelengths
      for(wvl_idx=1;wvl_idx<wvl_nbr;wvl_idx++){
	flx_slr_frc_blr[wvl_idx]=flx_slr_frc_blr[wvl_idx-1]+flx_slr_frc[wvl_idx]; // [frc] Fraction of solar flux at shorter wavelengths
      } // end loop over wvl 
      flx_slr_frc_blr[idx_max_wvl]=1.0; // [frc] Fraction of solar flux at shorter wavelengths
    }else{
      flx_slr_frc_blr[idx_max_wvl]=1.0; // [frc] Fraction of solar flux at shorter wavelengths
      for(wvl_idx=1;wvl_idx<wvl_nbr;wvl_idx++){
	flx_slr_frc_blr[wvl_idx]=flx_slr_frc_blr[wvl_idx-1]-flx_slr_frc[wvl_idx]; // [frc] Fraction of solar flux at shorter wavelengths
      } // end loop over wvl 
      flx_slr_frc_blr[idx_min_wvl]=0.0; // [frc] Fraction of solar flux at shorter wavelengths
    } // end else flg_ncr_wvl
  } // end else wvl_nbr > 1 and monotonic grid
  
  std::complex<prc_cmp> *idx_rfr_ffc_wgt=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Spectral flux-weighted effective refractive index
  prc_cmp *wvl_wgt=new prc_cmp[wvl_nbr]; // [m] Solar flux-weighted wavelength
  // Initialize weighted refractive indices
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    idx_rfr_ffc_wgt[wvl_idx]=std::complex<prc_cmp>(0.0,0.0); // [frc] Spectral flux-weighted effective refractive index
    wvl_wgt[wvl_idx]=0.0; // [m] Solar flux-weighted wavelength
  } // end loop over wvl 

  // Get tabulated refractive indices
  std::complex<prc_cmp> *idx_rfr_cor=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of core
  std::complex<prc_cmp> *idx_rfr_mdm=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of medium
  std::complex<prc_cmp> *idx_rfr_mnt=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of mantle
  std::complex<prc_cmp> *idx_rfr_mtx=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of matrix
  std::complex<prc_cmp> *idx_rfr_ncl=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of inclusion
  std::complex<prc_cmp> *idx_rfr_prt=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Refractive index of particle
  rcd+=cmp_cor.idx_rfr_get(wvl_ctr,idx_rfr_cor,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of core
  rcd+=cmp_mdm.idx_rfr_get(wvl_ctr,idx_rfr_mdm,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of medium
  rcd+=cmp_mnt.idx_rfr_get(wvl_ctr,idx_rfr_mnt,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of mantle
  rcd+=cmp_mtx.idx_rfr_get(wvl_ctr,idx_rfr_mtx,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of matrix
  rcd+=cmp_ncl.idx_rfr_get(wvl_ctr,idx_rfr_ncl,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of inclusion
  rcd+=cmp_prt.idx_rfr_get(wvl_ctr,idx_rfr_prt,wvl_nbr,wrn_ntp_flg); // [frc] Refractive index of particle

  // Always get and archive refractive indices of air for diagnostics
  std::complex<prc_cmp> *idx_rfr_air=new std::complex<prc_cmp>[wvl_nbr]; // O [frc] Refractive index of air
  rcd+=idx_rfr_air_get // [fnc] Compute refractive indices for air
    (static_cast<std::string>("air"), // I [sng] Composition of particle
     idx_rfr_air, // O [frc] Refractive index of air
     wvl_ctr, // I [m] Wavelength
     wvl_nbr); // I [nbr] Number of output wavelength bands

  prc_cmp *ncl_lst=new prc_cmp[max_lng(ncl_nbr,1L)]; // [sct] List of inclusion properties
  prc_cmp *vlm_frc_ncl_mlt=new prc_cmp[max_lng(ncl_nbr,1L)]; // [frc] Volume fraction(s) of inclusion(s)

  // Parse inclusions
  if(dbg_lvl > dbg_crr) std::cout << "Current number of inclusions is " << ncl_nbr << std::endl;
  if(ncl_nbr > 0){
    // Inclusion list specified by user
    // Currently inclusion properties are simply volume fractions of inclusions
    for(ncl_idx=0;ncl_idx<ncl_nbr;ncl_idx++) ncl_lst[ncl_idx]=static_cast<prc_cmp>(std::strtod(ncl_arg[ncl_idx].c_str(),(char **)NULL)); // [frc] Volume fraction(s) of inclusion(s)
  }else{
    // No inclusions were specified by user, use default inclusion
    ncl_nbr=1; // [nbr] Number of inclusions
    ncl_lst=&vlm_frc_ncl; // [sct] List of inclusions
  } // endif ncl_nbr > 0
  for(ncl_idx=0;ncl_idx<ncl_nbr;ncl_idx++) vlm_frc_ncl_mlt[ncl_idx]=ncl_lst[ncl_idx]; // [frc] Volume fraction(s) of inclusion(s)
  if(ncl_nbr == 4){ // Kludge to systematically assign fractions to components for blends
    vlm_frc_cor=vlm_frc_ncl_mlt[0];
    vlm_frc_mnt=vlm_frc_ncl_mlt[1];
    vlm_frc_ncl=vlm_frc_ncl_mlt[2];
    vlm_frc_prt=vlm_frc_ncl_mlt[3];
    // 20130119: fxm is this a bug or is there a reason vlm_frc_mtx > 1?
    vlm_frc_mtx=(1.0-(vlm_frc_cor+vlm_frc_mnt+vlm_frc_ncl+vlm_frc_prt));
  } // ncl_nbr != 4

  // Compute effective refractive indices
  std::complex<prc_cmp> *idx_rfr_ffc=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Effective refractive index of particle
  std::complex<prc_cmp> *idx_rfr_ffc_brg=new std::complex<prc_cmp>[wvl_nbr]; // O [frc] Effective refractive index, Bruggeman approximation
  std::complex<prc_cmp> *idx_rfr_ffc_mxg=new std::complex<prc_cmp>[wvl_nbr]; // O [frc] Effective refractive index, Maxwell Garnett approximation
  std::complex<prc_cmp> *idx_rfr_ffc_pmr=new std::complex<prc_cmp>[wvl_nbr]; // O [frc] Effective refractive index, partial molar refraction approximation
  std::complex<prc_cmp> *idx_rfr_ffc_vlw=new std::complex<prc_cmp>[wvl_nbr]; // O [frc] Effective refractive index, volume-weighted approximation
  // Compute all effective medium approximations on wavelength grid in CL1
  rcd+=idx_rfr_ffc_get // [fnc] Compute effective optical constants of composites
    (ncl_nbr, // I [nbr] Number of inclusions
     wvl_nbr, // I [nbr] Number of wavelengths
     ffc_mdm_typ, // I [enm] Effective medium type
     vlm_frc_ncl_mlt, // I [frc] Volume fraction(s) of inclusion(s)
     idx_rfr_cor, // I [frc] Refractive index of core
     idx_rfr_mdm, // I [frc] Refractive index of medium
     idx_rfr_mnt, // I [frc] Refractive index of mantle
     idx_rfr_mtx, // I [frc] Refractive index of matrix
     idx_rfr_ncl, // I [frc] Refractive index of inclusion
     idx_rfr_prt, // I [frc] Refractive index of particle
     idx_rfr_ffc_brg, // O [frc] Effective refractive index, Bruggeman approximation
     idx_rfr_ffc_mxg, // O [frc] Effective refractive index, Maxwell Garnett approximation
     idx_rfr_ffc_pmr, // O [frc] Effective refractive index, partial molar refraction approximation
     idx_rfr_ffc_vlw, // O [frc] Effective refractive index, volume-weighted approximation
     idx_rfr_ffc); // O [frc] Effective refractive index of particle

  // Initialize core and mantle sizes
  prc_cmp *rds_cor=new prc_cmp[sz_nbr]; // [m] Radius of core
  prc_cmp *rds_mnt=new prc_cmp[sz_nbr]; // [m] Radius of mantle
  for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    // Core radius measured from center of particle to core/mantle interface
    // Mantle radius measured from center of particle to outside surface
    rds_mnt[sz_idx]=sz_ctr_sph[sz_idx]; // [m] Radius of mantle
    // Derive core dimension from user-specified core radius fraction 
    rds_cor[sz_idx]=rds_mnt[sz_idx]*rds_frc_cor; // [m] Radius of core
    assert(rds_mnt[sz_idx] >= rds_cor[sz_idx]);
  } // end loop over sz

  // 20210418: Explanation for non-const wvl_bnd_sz_nbr is near sz_idx_dbg declaration
  //const long wvl_bnd_sz_nbr(wvl_nbr*bnd_nbr*sz_nbr); // [nbr] Number of wavelength/band/size loops
  long wvl_bnd_sz_nbr(wvl_nbr*bnd_nbr*sz_nbr); // [nbr] Number of wavelength/band/size loops
  const prc_cmp sz_prm_dbg(2.0*mth::cst_M_PIl*rds_ctr[sz_idx_dbg]/wvl_ctr[wvl_idx_dbg]); // [frc] Size parameter at debug size, wavelength
  const prc_cmp sz_prm_rds_min_wvl_max(2.0*mth::cst_M_PIl*rds_min[0]/wvl_max[wvl_nbr-1]); // [frc] Size parameter, smallest particle, longest wavelength
  const prc_cmp sz_prm_rds_max_wvl_min(2.0*mth::cst_M_PIl*rds_max[sz_nbr-1]/wvl_min[0]); // [frc] Size parameter, largest particle, shortest wavelength
  const prc_cmp sz_prm_rsn_avg((sz_prm_rds_max_wvl_min-sz_prm_rds_min_wvl_max)/wvl_bnd_sz_nbr); // [frc] Size parameter, mean computational resolution
  prc_cmp *sz_prm_swa=new prc_cmp[wvl_nbr]; // [frc] Size parameter at rds_swa
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    sz_prm_swa[wvl_idx]=2.0*mth::cst_M_PIl*rds_swa/wvl_ctr[wvl_idx]; // [frc] Size parameter at rds_swa
  } // end loop over wvl 
  
  if(true){
    std::cout << "Initialization state:" << std::endl;
    std::cout << "Wet Deposition:" << std::endl;
    //    std::cout << "  Raindrop: Number Median Diameter = " << dmt_nmr_pcp*1.0e6 << " um, Mass Median Diameter = " << dmt_vmr_pcp*1.0e6 << " um, Fall speed = " << vlc_pcp << " m s-1, Reynolds number = " << ryn_nbr_pcp << std::endl;
    std::cout << "  Raindrop: Diameter = " << dmt_pcp[dsd_idx_dbg]*1.0e6 << " um, Mass = " << mss_pcp[dsd_idx_dbg]*1.0e9 << " ug, Fall speed = " << vlc_pcp[dsd_idx_dbg] << " m s-1, Reynolds number = " << ryn_nbr_pcp[dsd_idx_dbg] << std::endl;
    std::cout << "  Precipitation: Concentration = " << cnc_nbr_pcp_rsl << " # m-3, Number flux = " << flx_cnc_pcp_rsl << " # m-2 s-1" << std::endl;
    std::cout << "  Analytic concentrations: number = " << cnc_nbr_pcp_anl << " # m-3, surface = " << dst_lgn_pcp.sfc_get() << " m2 m-3, volume =  " << dst_lgn_pcp.vlm_get() << " m3 m-3, mass = " << dns_pcp*dst_lgn_pcp.vlm_get() << " kg m-3" << std::endl;
    std::cout << "  Resolved concentrations: number = " << cnc_nbr_pcp_rsl << " # m-3, surface = " << sfc_pcp_rsl << " m2 m-3, volume = " << vlm_pcp_rsl << " m3 m-3, mass = " << mss_pcp_rsl << " kg m-3" << std::endl;
    std::cout << "  Fractions resolved: number = " << cnc_nbr_rsl_frc_pcp << ", surface = " << sfc_rsl_frc_pcp << ", volume & mass = " << vlm_rsl_frc_pcp << std::endl;
    std::cout << "  Precipitation: " << flx_vlm_pcp_rsl << " m3 m-2 s-1 = " << flx_vlm_pcp_rsl*1000.0*3600.0 << " mm hr-1 = " << flx_mss_pcp_rsl << " kg m-2 s-1" << std::endl;
    std::cout << "  Density ratio (H2O/aer) = " << dns_rat << ", Viscosity ratio (H2O/air) = " << vsc_rat << std::endl;
    std::cout << "  Efficiencies for (raindrop, aerosol) pair = (" << dmt_pcp[dsd_idx_dbg]*1.0e6 << "," << dmt_ctr[sz_idx_dbg]*1.0e6 << ") um: Diffusion = " << cll_fsh_brn_dff(dsd_idx_dbg,sz_idx_dbg) << ", Interception = " << cll_fsh_ntc(dsd_idx_dbg,sz_idx_dbg) << ", Impaction = " << cll_fsh_mpc(dsd_idx_dbg,sz_idx_dbg) << ", Total = " << cll_fsh(dsd_idx_dbg,sz_idx_dbg) << std::endl;
    std::cout << "  Collision = " << cll_fsh(dsd_idx_dbg,sz_idx_dbg) << ", Sticking = " << stc_fsh[sz_idx_dbg] << ", Collection = " << clc_fsh(dsd_idx_dbg,sz_idx_dbg) << std::endl;
    std::cout << "  Scavenging coefficient = " << scv_cff[sz_idx_dbg] << " s-1 = " << scv_cff[sz_idx_dbg]*3600.0 << " hr-1" << std::endl;
    std::cout << "  Mass mean scavenging coefficient = " << scv_cff_mss_avg << " s-1, normalized by precipitation volume flux = " << 0.001*scv_cff_mss_avg/flx_vlm_pcp_rsl << " mm-1" << std::endl;
    std::cout << "  Number mean scavenging coefficient = " << scv_cff_nbr_avg << " s-1, normalized by precipitation volume flux = " << 0.001*scv_cff_nbr_avg/flx_vlm_pcp_rsl << " mm-1" << std::endl;
    std::cout << "Environment:" << std::endl;
    std::cout << "  Solar constant = " << slr_cst << " W m-2" << std::endl;
    std::cout << "  Solar spectrum is " << flx_slr_src.dsc_get() << " from " << (flx_slr_src.fl_slr_spc_get() != "" ? flx_slr_src.fl_slr_spc_get() : " a function") << std::endl;
    std::cout << "  Fraction of solar flux in " << wvl_mnm*1.0e6 << "--" << wvl_mxm*1.0e6 << " um is " << flx_slr_src.flx_frc_get(wvl_min[wvl_idx_dbg],wvl_max[wvl_idx_dbg]) << std::endl;
    std::cout << "  Blackbody temperature of particles = " << tpt_prt << " K" << std::endl;
    std::cout << "  Blackbody temperature of radiation = " << tpt_bbd_wgt << " K, hemispheric blackbody emission = " << spc_bbd.flx_ttl() << " W m-2, fraction in " << wvl_mnm*1.0e6 << "--" << wvl_mxm*1.0e6 << " um is " << spc_bbd.flx_frc_get(wvl_mnm,wvl_mxm) << " = " << spc_bbd.flx_ttl()*spc_bbd.flx_frc_get(wvl_mnm,wvl_mxm) << " W m-2" << std::endl;
    std::cout << "  Pressure = " << prs_mdp/100.0 << " mb, Temperature = " << tpt_mdp << " K, Density = " << dns_mdp << " kg m-3" << std::endl;
    std::cout << "  Relative humidity w/r/t liquid H2O: = " << RH_lqd << ", ice H2O = " << RH_ice << ", q(H2O) = " << q_H2O_vpr << " kg kg-1" << std::endl;
    std::cout << "  Kinematic vsc. = " << vsc_knm_atm << " m2 s-1, Dynamic vsc. = " << vsc_dyn_atm << " kg m-1 s-1, Mean free path of air = " << mfp_atm << " m" << std::endl;
    std::cout << "Surface:" << std::endl;
    std::cout << "  Orography = " << oro << " = " << oro_sng_get(oro) << std::endl;
    std::cout << "  Dry land fraction = " << lnd_frc_dry << ", Bare ground fraction = " << lnd_frc_mbl << std::endl;
    std::cout << "  "+sfc_typ_dsc_get(sfc_typ) << std::endl;
    using lsm::rgh_lng; // [m] Roughness length momentum
    std::cout << "  Sub-grid \"soil type\" properties: Type " << soi_typ << " = " << soi_typ_sng_get(soi_typ) << ", Bare surface roughness = " << rgh_lng[soi_typ] << " m, Bare ground emissivity = " << msv_gnd << std::endl;
    std::cout << "  Soil texture: " << mss_frc_cly*100.0 << "% clay, " << mss_frc_slt*100.0 << "% silt, " << mss_frc_snd*100.0 << "% sand" << std::endl;
    std::cout << "  Snow height (liquid) = " << snw_hgt_lqd*100.0 << " cm, Snow height (geometric) = " << snw_hgt*100.0 << " cm, Snow fraction = " << snw_frc << std::endl;
    std::cout << "  Roughness length deposition = " << rgh_mmn_dps << " m, Zero-plane displacement height deposition = " << hgt_zpd_dps << " m" << std::endl; 
    std::cout << "  Soil temperature = " << tpt_soi << " K, Thermal conductivity = " << cnd_trm_soi << " W m-1 K-1" << std::endl;
    std::cout << "  Volumetric water content = " << vwc_sfc << " m3 m-3, Saturation fraction = " << vwc_rel*100.0 << "%" << std::endl;
    std::cout << "  Poreless soil particle density = " << dns_prt_sfc << " kg m-3, Bulk dry soil density = " << dns_blk_dry << " kg m-3, Bulk soil density = " << dns_blk_sfc << " kg m-3" << std::endl;
    std::cout << "  Gravimetric water content = " << gwc_sfc << " kg kg-1" << std::endl;
    std::cout << "  Dry volumetric water content (no E-T) = " << vwc_dry << " m3 m-3, E-T optimal volumetric water content = " << vwc_opt << " m3 m-3" << std::endl;
    std::cout << "  Transfer efficiency of vapor from soil to atmosphere = " << trn_fsh_vpr_soi_atm*100.0 << " %" << std::endl;
    std::cout << "Wind properties:" << std::endl;
    std::cout << "  Wind speeds at z = " << hgt_mdp << " m: U = " << wnd_znl_mdp << " m s-1, V = " << wnd_mrd_mdp << " m s-1, Total = " << wnd_mdp << " m s-1" << std::endl;
    std::cout << "  Wind speeds at z = " << hgt_rfr << " m: U = " << wnd_znl_rfr << " m s-1, V = " << wnd_mrd_rfr << " m s-1, Total = " << wnd_rfr_mbl << " m s-1" << std::endl;
    std::cout << "  Wind friction speeds: Non-saltating = " << wnd_frc_mbl << " m s-1, Owen's effect = " << wnd_frc_slt_dlt << " m s-1, Saltating = " << wnd_frc_slt << " m s-1" << std::endl;
    std::cout << "  Threshold friction velocity increase factor from moisture = " << frc_thr_ncr_wtr << std::endl;
    std::cout << "  Threshold friction velocity increase factor from roughness = " << frc_thr_ncr_drg << std::endl;
    std::cout << "  Threshold friction speeds for saltation: Dry, flat ground = " << wnd_frc_thr_slt_dry_flt << " m s-1, Including moisture = " << wnd_frc_thr_slt_dry_flt*frc_thr_ncr_wtr << " m s-1, Including roughness = " << wnd_frc_thr_slt_dry_flt*frc_thr_ncr_drg << " m s-1, Including moisture and roughness = " << wnd_frc_thr_slt << " m s-1" << std::endl;
    std::cout << "  Threshold wind speeds: At reference = " << hgt_rfr << " m = " << wnd_rfr_thr_slt << " m s-1, at midlayer = " << hgt_mdp << " m = " << wnd_mdp_thr_slt << " m s-1" << std::endl;
    std::cout << "  Weibull distribution: Scaling parameter = " << wbl_scl << " m s-1, Shape parameter = " << wbl_shp << std::endl;
    std::cout << "Dust mobilization:" << std::endl;
    std::cout << "  Vertically integrated streamwise mass flux = " << flx_mss_hrz_slt_ttl << " kg m-1 s-1" << std::endl; 
    std::cout << "  Ratio of vertical dust flux to vertically integrated streamwise mass flux = " << dst_slt_flx_rat_ttl << " m-1" << std::endl; 
    std::cout << "  Vertical dust flux = " << flx_mss_vrt_dst_ttl << " kg m-2 s-1" << std::endl; 
    std::cout << "Surface Energy Budget mobilization:" << std::endl;
    std::cout << "  Flux LW up = " << flx_LW_upw_sfc << " W m-2, Flux LW down = " << flx_LW_dwn_sfc << " W m-2, Flux LW net = " << flx_LW_dwn_sfc-flx_LW_upw_sfc << " W m-2" << std::endl;
    std::cout << "  Flux SW net ground = " << flx_SW_net_gnd << " W m-2, Flux SW net veg = " << flx_SW_net_vgt << " W m-2, Flux SW net total = " << flx_SW_net_gnd+flx_SW_net_vgt << " W m-2" << std::endl;
    std::cout << "  Latent heat = " << flx_ltn << " W m-2, Sensible heat to atmosphere = " << flx_sns_atm << " W m-2, Sensible heat to soil = " << flx_sns_gnd << " W m-2, Flux snow melt = " << flx_snw_mlt << " W m-2" << std::endl;
    std::cout << "  Water vapor flux to atmosphere = " << flx_q_H2O*1.0e6 << " mg m-2 s-1" << std::endl;
    std::cout << "  Temperatures at z = " << hgt_mdp << " m: T = " << tpt_mdp << " K, Tv = " << tpt_vrt_mdp << " K, Theta = " << tpt_ptn_mdp << " K, Thetav = " << tpt_ptn_vrt_mdp << " K" << std::endl;
    std::cout << "  Temperatures: GCM = " << tpt_mdp << " K, Screen (z0h+d+2m) = " << tpt_ash_p2m << " K, Aerodynamic (z0m+d)= " << tpt_aer << " K, Surface (z0h+d) = " << tpt_ash << " K, Vegetation = " << tpt_vgt << " K, Ground = " << tpt_gnd << " K, Soil = " << tpt_soi << " K, Effective = " << tpt_ffc << " K" << std::endl;
    std::cout << "Surface Exchange Properties:" << std::endl;
    std::cout << "  Roughness lengths: Mobilization = " << rgh_mmn_mbl << " m, Smooth = " << rgh_mmn_smt << " m, Deposition = " << rgh_mmn_dps << " m" << std::endl; 
    std::cout << "  Transfer properties between hgt_mdp = " << hgt_mdp << " m and z0m+zpd = " << hgt_asm << " m, and z0h+zpd = z0v+zpd = " << hgt_ash << " m:" << std::endl;
    std::cout << "  Dimensionless exchange coefficients: Momentum (neutral) = " << cff_xch_mmn_ntr << ", Momentum (actual) = " << cff_xch_mmn << ", Heat = " << cff_xch_heat << ", Vapor = " << cff_xch_vpr << std::endl;
    std::cout << "  Aerodynamic resistances to surface transfer of momentum:  Momentum = " << rss_aer_mmn << " s m-1, Heat = " << rss_aer_heat << " s m-1, Vapor = " << rss_aer_vpr << " s m-1" << std::endl;
    std::cout << "Deposition properties:" << std::endl;
    std::cout << "  Roughness length: " << rgh_mmn_dps << " m" << std::endl; 
    std::cout << "  Monin-Obukhov length: " << mno_lng_dps << " m" << std::endl; 
    std::cout << "  Wind friction speed: " << wnd_frc_dps << " m s-1" << std::endl; 
    std::cout << "Hygroscopic growth:" << std::endl;
    std::cout << "  Mean linear mass increase coefficient = " << mlmic << ", (mass water)/(mass particle) = " << mss_wtr_rat[sz_idx_dbg] << std::endl;
    std::cout << "  Curvature (Kelvin) factor = " << svp_fct_klv[sz_idx_dbg] << ", Solute factor = " << svp_fct_slt[sz_idx_dbg] << std::endl;
    std::cout << "  Activity coefficient = " << act_cff[sz_idx_dbg] << ", Mass fraction solute = " << mss_frc_solute[sz_idx_dbg] << std::endl;
    std::cout << "Aerosol chemistry with HNO3:" << std::endl;
    std::cout << "  Accomodation coefficient (alpha) = " << mss_acm_cff_HNO3_aer << ", Uptake coefficient (gamma) = " << mss_upt_cff_HNO3_aer << std::endl;
    std::cout << "  " << "HNO3 volume mixing ratio = " << vmr_HNO3_gas*1.0e9 << " ppbv, HNO3 number concentration = " << cnc_HNO3_gas << " mlc m-3, HNO3 mass mixing ratio = " << q_HNO3_gas << " kg kg-1" << std::endl;
    std::cout << "  " << "Binary diffusivity HNO3 in air = " << dff_HNO3_air << " m2 s-1, Mean free path of HNO3 in air = " << mfp_HNO3_air << " m, Thermal speed of HNO3 = " << vlc_mwb_HNO3 << " m s-1" << std::endl;
    std::cout << "  " << "Properties for HNO3 interaction with aerosol diameter = " << dmt_ctr[sz_idx_dbg]*1.0e6 << " um:" << std::endl; 
    std::cout << "\tKnudsen number = " << knd_nbr_HNO3_air[sz_idx_dbg] << std::endl;
    std::cout << "\tVentilation correction = " << vnt_crc_aer[sz_idx_dbg] << std::endl; 
    std::cout << "\tKinetic regime Brownian diffusion rate to aerosol (includes alpha) = " << dff_HNO3_aer_knt[sz_idx_dbg] << " m3 s-1 prt-1" << std::endl; 
    std::cout << "\tContinuum regime fluid diffusion rate to aerosol (excludes alpha) = " << dff_HNO3_aer_cnt[sz_idx_dbg] << " m3 s-1 prt-1" << std::endl; 
    std::cout << "\tCorrection to kinetic regime ballistic approximation due to diffusion effects = " << knt_rgm_dff_crc_HNO3_aer[sz_idx_dbg] << std::endl; 
    std::cout << "\tCorrection to continuum regime diffusion approximation due to ballistic effects = " << cnt_rgm_bll_crc_HNO3_aer[sz_idx_dbg] << std::endl; 
    std::cout << "\tTransition regime net diffusion rate to aerosol = " << dff_HNO3_aer[sz_idx_dbg] << " m3 s-1 prt-1" << std::endl; 
    std::cout << "\tPseudo first order rate coefficient = " << rxrc_HNO3_aer[sz_idx_dbg] << " s-1, Timescale = " << tau_rxrc_HNO3_aer[sz_idx_dbg] << " s" << std::endl; 
    std::cout << "\tPseudo first order rate coefficient, kinetic appoximation = " << rxrc_HNO3_aer_knt[sz_idx_dbg] << " s-1, Timescale = " << tau_rxrc_HNO3_aer_knt[sz_idx_dbg] << " s" << std::endl; 
    std::cout << "\tMean rate of removal by aerosol = " << rxr_HNO3_gas_aer_vmr[sz_idx_dbg] << " mlc mlc-1 s-1" << std::endl;
    std::cout << "  " << "Total pseudo first order rate coefficient for HNO3 removal by aerosol = " << rxrc_HNO3_aer_ttl << " s-1, Timescale = " << tau_rxrc_HNO3_aer_ttl << " s" << std::endl;
    std::cout << "  " << "Total pseudo first order rate coefficient for HNO3 removal by aerosol, kinetic approximation = " << rxrc_HNO3_aer_knt_ttl << " s-1, Timescale = " << tau_rxrc_HNO3_aer_knt_ttl << " s" << std::endl;
    std::cout << "  " << "Total rates of removal by aerosol: Mean = " << rxr_HNO3_gas_aer_vmr_ttl << " mlc mlc-1 s-1, Instantaneous = " << rxr_nst_HNO3_gas_aer_vmr_ttl << " mlc mlc-1 s-1" << std::endl;
    std::cout << "  " << "Timescales for aqueous chemistry of HNO3 with liquid H2O diameter = " << dmt_ctr[sz_idx_dbg]*1.0e6 << " um: \n\tPresumed e-folding time due to irreversible uptake = " << tau_rxrc_HNO3_aer[sz_idx_dbg] << " s, \n\tGas phase diffusion to droplet surface = " << tau_gpd_HNO3_aer[sz_idx_dbg] << " s, \n\tEquilibrium at gas-droplet interface = " << tau_ntf_eqm_HNO3_aer[sz_idx_dbg] << " s, \n\tAqueous dissociation HNO3(aq)<-->NO3-(aq) = " << tau_aqs_dss_HNO3_aer[sz_idx_dbg] << " s, \n\tAqueous phase diffusion = " << tau_aqs_dff_HNO3_aer[sz_idx_dbg] << " s" << std::endl;
    std::cout << "  " << "Properties of liquid phase: \n\tHenry's Law coefficient at T = 298 K is " << cff_hnr_HNO3_H2O_298K << " mol ltr-1 atm-1, at " <<  tpt_mdp << " K is " << cff_hnr_HNO3_H2O << " mol ltr-1 atm-1" << std::endl;
    std::cout << "Aerosol composition:" << std::endl;
    std::cout << "  " << cmp_prt.dsc_get() << " aerosol with density " << dns_prt << " kg m-3" << std::endl;
    std::cout << "  Density(aerosol)/density(atmosphere) = " << dns_prt/dns_mdp << std::endl;
    if(ss_alb_flg) std::cout << "  Single scattering properties scaled to user-specified value of omega = " << ss_alb_cmd_ln << " at wvl[" << wvl_idx_dbg << "] = " << wvl_ctr[wvl_idx_dbg]*1.0e6 << " um" << std::endl;
    std::cout << "Wavelength grid:" << std::endl;
    std::cout << "  Grid type is " << wvl_grd_sng << std::endl;
    std::cout << "  Shortest wavelength = " << wvl_mnm*1.0e6 << " um = " << 1.0e-2/wvl_mnm << " cm-1" << std::endl;
    std::cout << "  Longest wavelength = " << wvl_mxm*1.0e6 << " um = " << 1.0e-2/wvl_mxm << " cm-1" << std::endl;
    std::cout << "  Center wavelength = " << 0.5*(wvl_mnm+wvl_mxm)*1.0e6 << " um = " << 2.0e-2/(wvl_mxm+wvl_mnm) << " cm-1" << std::endl;
    std::cout << "  " << wvl_nbr << " output spectral bin(s) are each evaluated using " << bnd_nbr << " sub-bin(s)" << std::endl;
    std::cout << "Raindrop Mode:" << std::endl;
    psd_pcp.psd_prn(&psd_pcp,1L); // [fnc] Print formatted list of properties of particle size distribution(s)
    std::cout << "Aerosol Mode(s):" << std::endl;
    std::cout << "  Shape is " << psd_lst[0].typ_get() << std::endl;
    if(psd_lst[0].typ_get() == "lognormal"){
      psd_lst[0].psd_prn(psd_lst,psd_nbr); // [fnc] Print formatted list of properties of particle size distribution(s)
    }else if(psd_lst[0].typ_get() == "gamma"){
      // fxm: Gamma distributions needs to be re-validated, they are broken now
      std::cout << "rds_ffc = " << dst_gmm_aer.rds_ffc_get()*1.0e6 << " um, var_ffc = " << dst_gmm_aer.var_ffc_get() << ", cnc_nbr_anl = " << dst_gmm_aer.cnc_nbr_anl_get() << " m-3" << std::endl;
    }else{ 
      err_prn(prg_nm,sbr_nm,"Unknown psd_lst[0].typ_get()");
    } // end else
    std::cout << "  Size grid is " << psd_szgrd.typ_get() << " with " << sz_nbr << " bins spanning domain " << sz_mnm*1.0e6 << "--" << sz_mxm*1.0e6 << " um" << std::endl;
    std::cout << "  Weighted radii, resolved: cnc = " << rds_nwr*1.0e6 << " um, sfc = " << rds_swr*1.0e6 << " um, vlm = " << rds_vwr*1.0e6 << " um" << std::endl;
    std::cout << "  Weighted settling velocities: cnc = " << vlc_grv_nwr << " m s-1, vlm = " << vlc_grv_vwr << " m s-1" << std::endl;
    std::cout << "  Analytic concentrations: number = " << cnc_nbr_anl << " # m-3, surface = " << sfc_anl << " m2 m-3, volume =  " << vlm_anl << " m3 m-3, mass = " << mss_anl << " kg m-3" << std::endl;
    std::cout << "  Resolved concentrations: number = " << cnc_nbr_rsl << " # m-3, surface = " << sfc_rsl << " m2 m-3, volume = " << vlm_rsl << " m3 m-3, mass = " << mss_rsl << " kg m-3" << std::endl;
    std::cout << "  Fractions resolved: number = " << cnc_nbr_rsl_frc << ", surface = " << sfc_rsl_frc << ", volume & mass = " << vlm_rsl_frc << std::endl;
    std::cout << "  Specific concentration analytic: number = " << cnc_spc_anl << " # kg-1, surface = " << sfc_spc_anl << " m2 kg-1, volume = " << vlm_spc_anl << " m3 kg-1" << std::endl;
    std::cout << "  Specific concentration resolved: number = " << cnc_spc_rsl << " # kg-1, surface = " << sfc_spc_rsl << " m2 kg-1, volume = " << vlm_spc_rsl << " m3 kg-1" << std::endl;
    std::cout << "Aspherical corrections:" << std::endl;
    std::cout << "  Ellipsoidal aspect ratio (for aerodynamics): " << asp_rat_lps[sz_idx_dbg]  << std::endl;
    std::cout << "  Hexagonal prism aspect ratio (for optics): " << asp_rat_hxg[sz_idx_dbg] << std::endl;
    std::cout << "Mie options:" << std::endl;
    std::cout << "  Mie solver algorithm = " << slv_sng << std::endl;
    std::cout << "  Order of phase function Legendre expansion = " << lgn_nbr << std::endl;
    std::cout << "  Angles per hemisphere = " << ngl_nbr << std::endl;
    std::cout << "Diagnostics:" << std::endl;
    std::cout << "  Diagnostic wavelength is bin " << wvl_idx_dbg << " from " << wvl_min[wvl_idx_dbg]*1.0e6 << "--" << wvl_max[wvl_idx_dbg]*1.0e6 << " um = " << wvn_min[wvl_idx_dbg] << "--" << wvn_max[wvl_idx_dbg] << " cm-1" << std::endl;
    std::cout << "  Diagnostic size is bin " << sz_idx_dbg << " from " << sz_min[sz_idx_dbg]*1.0e6 << "--" << sz_max[sz_idx_dbg]*1.0e6 << " um" << std::endl;
    std::cout << "  Resolved concentrations for this size bin: number = " << cnc[sz_idx_dbg] << " # m-3, surface = " << sfc[sz_idx_dbg]*cnc[sz_idx_dbg] << " m2 m-3, volume = " << vlm[sz_idx_dbg]*cnc[sz_idx_dbg] << " m3 m-3, mass = " << mss[sz_idx_dbg]*cnc[sz_idx_dbg] << " kg m-3" << std::endl;
    std::cout << "  Size Parameter Minimum (2.0*pi*rds_min/wvl_max): " << sz_prm_rds_min_wvl_max << ", Maximum (2.0*pi*rds_max/wvl_min): " << sz_prm_rds_max_wvl_min << ", Number = " << wvl_bnd_sz_nbr << std::endl;
    std::cout << "  Size Parameter Mean resolution (chi_max-chi_min)/nbr: " << sz_prm_rsn_avg << std::endl;
    std::cout << "Composition/Coatings/Optics: " << std::endl;
    std::cout << "  Coatings: " << (coat_flg ? "Active" : "Inactive") << std::endl;
    std::cout << "    Radius of core = " << rds_cor[sz_idx_dbg]*1.0e6 << " um" << std::endl;
    std::cout << "    Radius of mantle = " << rds_mnt[sz_idx_dbg]*1.0e6 << " um" << std::endl;
    std::cout << "    Size parameter of core = " << 2.0*mth::cst_M_PIl*rds_cor[sz_idx_dbg]*real(idx_rfr_mdm[wvl_idx_dbg])/wvl_ctr[wvl_idx_dbg] << std::endl;
    std::cout << "    Size parameter of mantle = " << 2.0*mth::cst_M_PIl*rds_mnt[sz_idx_dbg]*real(idx_rfr_mdm[wvl_idx_dbg])/wvl_ctr[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive indices of core = " << cmp_cor.dsc_get() << " ";
    if(idx_rfr_cor_usr_flg) std::cout << "are user-specified as " << idx_rfr_cor_usr << std::endl; else std::cout << "from " << fl_idx_rfr_cor << std::endl;
    std::cout << "  Refractive indices of medium = " << cmp_mdm.dsc_get() << " ";
    if(idx_rfr_mdm_usr_flg) std::cout << "are user-specified as " << idx_rfr_mdm_usr << std::endl; else std::cout << "from " << fl_idx_rfr_mdm << std::endl;
    std::cout << "  Refractive indices of mantle = " << cmp_mnt.dsc_get() << " ";
    if(idx_rfr_mnt_usr_flg) std::cout << "are user-specified as " << idx_rfr_mnt_usr << std::endl; else std::cout << "from " << fl_idx_rfr_mnt << std::endl;
    std::cout << "  Refractive indices of matrix = " << cmp_mtx.dsc_get() << " ";
    if(idx_rfr_mtx_usr_flg) std::cout << "are user-specified as " << idx_rfr_mtx_usr << std::endl; else std::cout << "from " << fl_idx_rfr_mtx << std::endl;
    std::cout << "  Refractive indices of inclusion = " << cmp_ncl.dsc_get() << " ";
    if(idx_rfr_ncl_usr_flg) std::cout << "are user-specified as " << idx_rfr_ncl_usr << std::endl; else std::cout << "from " << fl_idx_rfr_ncl << std::endl;
    std::cout << "  Refractive indices of particle = " << cmp_prt.dsc_get() << " ";
    if(idx_rfr_prt_usr_flg) std::cout << "are user-specified as " << idx_rfr_prt_usr << std::endl; else std::cout << "from " << fl_idx_rfr_prt << std::endl;
    std::cout << "  Refractive index of core = " << idx_rfr_cor[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive index of medium = " << idx_rfr_mdm[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive index of mantle = " << idx_rfr_mnt[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive index of matrix = " << idx_rfr_mtx[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive index of inclusion = " << idx_rfr_ncl[wvl_idx_dbg] << std::endl;
    std::cout << "  Refractive index of particle = " << idx_rfr_prt[wvl_idx_dbg] << std::endl;
    std::cout << "  Effective refractive index of particle = " << idx_rfr_ffc[wvl_idx_dbg] << std::endl;
    std::cout << "  Effective medium approximation: " << (ffc_mdm_flg ? "Active---using type "+ffc_mdm_sng[ffc_mdm_typ] : "Inactive") << std::endl;
    std::cout << "    Diagnostic effective refractive indices from all EMAs: Bruggeman = " << idx_rfr_ffc_brg[wvl_idx_dbg] << ", Maxwell Garnett = " << idx_rfr_ffc_mxg[wvl_idx_dbg] << ", Partial molar refraction = " << idx_rfr_ffc_pmr[wvl_idx_dbg] << ", Volume-weighted = " << idx_rfr_ffc_vlw[wvl_idx_dbg] << std::endl;
    if(slf_tst_flg && slf_tst_typ == ViC98){
      std::cerr << "\nCompare these results to Videen & Chylek (1998), \n\"Scattering by a composite sphere with an absorbing inclusion and effective medium approximations\", Opt. Comm., 158, pp. 1--6\n" << std::endl;
      std::cerr << "p. 5 Tbl. 3 quotes values for host sphere of water n=(1.335,0.0), size equal to wavelength, containing a carbon inclusion n=(1.94,0.66).\nVolume fraction V=0.10:  n(EEMA) = (1.337,0.00367), n(Bruggeman) = (1.398,0.0537), n(Maxwell Garnett) = (1.399,0.0515)\n" << std::endl;
    } // endif self-test
    std::cout << "    Core fractions: rds_frc_cor = " << rds_frc_cor << ", vlm_frc_cor = " << vlm_frc_cor << ", mss_frc_cor = " << mss_frc_cor << std::endl;
    std::cout << "    Inclusion fractions: rds_frc_ncl = " << rds_frc_ncl << ", vlm_frc_ncl = " << vlm_frc_ncl << ", mss_frc_ncl = " << mss_frc_ncl << std::endl;
    std::cout << "  Enhanced absorption by inclusions in weakly-absorbing spheres: " << (abs_ncl_wk_mdm_flg ? "Active, treated with method of MaS99" : "Inactive") << std::endl;
    std::cout << std::endl;
  } /* end if true */

  if(dbg_lvl >= dbg_io){
    std::cerr << "Aerodynamic size distribution diagnostics from main():" << std::endl;
    if(psd_lst[0].typ_get() == "lognormal"){
      std::cerr << "rds_nma =  " << dst_lgn_aer.rds_nma_get() << ", gsd = " << dst_lgn_aer.gsd_anl_get() << ", cnc_nbr_anl =  " << dst_lgn_aer.cnc_nbr_anl_get() << std::endl;
    }else if(psd_lst[0].typ_get() == "gamma"){
      std::cerr << "rds_ffc =  " << dst_gmm_aer.rds_ffc_get() << ", var_ffc = " << dst_gmm_aer.var_ffc_get() << ", cnc_nbr_anl =  " << dst_gmm_aer.cnc_nbr_anl_get() << std::endl;
    }else{
      err_prn(prg_nm,sbr_nm,psd_lst[0].typ_get()+" is not a valid psd_cls::Dst_typ");
    } // end if
    std::cerr << "\nidx\tr\tdr\tn(r)\tN(r)\tN<r" << std::endl;
    for(idx=0;idx<sz_nbr;idx++){
      std::cerr << " " << idx << " " << sz_ctr[idx] << " " << sz_dlt[idx] << " " << dst[idx] << " " << cnc[idx] << " " << cnc_nbr_rsl << std::endl;
    } // end loop over sz
    std::cerr << "cnc_nbr_rsl = " << cnc_nbr_rsl << ", cnc_nbr_rsl_frc = " << cnc_nbr_rsl/dst_lgn_aer.cnc_nbr_anl_get() << " cnc_nbr_anl = " << dst_lgn_aer.cnc_nbr_anl_get() << std::endl;
  } // end if dbg

  if(dbg_lvl == dbg_io){
    dbg_prn("Size grid in main():");
    std::cerr << "sz_idx_dbg = " << sz_idx_dbg << std::endl;
    std::cerr << "idx\t" << "sz_min\t" << "sz_max\t" << "sz_grd\t" << std::endl;
    for(idx=0;idx<sz_nbr+1;idx++) std::cerr << idx << "\t" << (idx < sz_nbr ? sz_min[idx] : 0.0) << "\t" << (idx < sz_nbr ? sz_max[idx] : 0.0) << "\t" << sz_grd[idx] << std::endl;
    dbg_prn("Wavelength grid in main():");
    std::cerr << "wvl_idx_dbg = " << wvl_idx_dbg << std::endl;
    std::cerr << "idx\t" << "wvl_min\t" << "wvl_max\t" << "wvl_grd\t" << std::endl;
    for(idx=0;idx<wvl_nbr+1;idx++) std::cerr << idx << "\t" << (idx < wvl_nbr ? wvl_min[idx] : 0.0) << "\t" << (idx < wvl_nbr ? wvl_max[idx] : 0.0) << "\t" << wvl_grd[idx] << std::endl;

    dbg_prn("Wavenumber grid in main():");
    std::cerr << "idx\t" << "wvn_min\t" << "wvn_max\t" << "wvn_grd\t" << std::endl;
    for(idx=0;idx<wvl_nbr+1;idx++) std::cerr << idx << "\t" << (idx < wvl_nbr ? wvn_min[idx] : 0.0) << "\t" << (idx < wvl_nbr ? wvn_max[idx] : 0.0) << "\t" << wvn_grd[idx] << std::endl;
  } // endif dbg

  if(dbg_lvl == dbg_old){
    dbg_prn("Testing spc_bbd() from main()");
    std::cerr << "Temperature is " << spc_bbd.tpt_get() << " K" << std::endl;
    std::cerr << "spc_bbd.flx_ttl() returns " << spc_bbd.flx_ttl() << " W m-2" << std::endl;
    std::cerr << "Wavelength is " << cpv_foo*1.0e6 << " um" << std::endl;
    std::cerr << "spc_bbd.eval(" << cpv_foo << ") returns " << spc_bbd.eval(cpv_foo) << " W m-2 m-1 sr-1" << std::endl;
  } // endif dbg

  /* Create arrays to hold integrated optical properties
     These arrays arr archived to output file
     Hence, all these arrays are shared */
  prc_cmp *abs_fct_MaS99=new prc_cmp[wvl_nbr]; // [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
  prc_cmp *asm_prm=new prc_cmp[wvl_nbr]; // [frc] Asymmetry parameter
  prc_cmp *ang_xpn=new prc_cmp[wvl_nbr]; // [frc] Angstrom exponent
  prc_cmp *abs_cff_vlm=new prc_cmp[wvl_nbr]; // [m2 m-3] Volume absorption coefficient
  prc_cmp *sca_cff_vlm=new prc_cmp[wvl_nbr]; // [m2 m-3] Volume scattering coefficient
  prc_cmp *ext_cff_vlm=new prc_cmp[wvl_nbr]; // [m2 m-3] Volume extinction coefficient
  prc_cmp *bck_cff_vlm=new prc_cmp[wvl_nbr]; // [m2 m-3] Volume backscattering coefficient
  prc_cmp *abs_cff_mss=new prc_cmp[wvl_nbr]; // [m2 kg-1] Mass absorption coefficient
  prc_cmp *ext_cff_mss=new prc_cmp[wvl_nbr]; // [m2 kg-1] Mass extinction coefficient
  prc_cmp *bck_cff_mss=new prc_cmp[wvl_nbr]; // [m2 kg-1] Mass backscattering coefficient
  prc_cmp *sca_cff_mss=new prc_cmp[wvl_nbr]; // [m2 kg-1] Mass scattering coefficient
  prc_cmp *ss_alb=new prc_cmp[wvl_nbr]; // [frc] Single scattering albedo 
  prc_cmp *ss_co_alb=new prc_cmp[wvl_nbr]; // [frc] Single scattering co-albedo
//   prc_cmp *abs_xsx=new prc_cmp[wvl_nbr];
//   prc_cmp *ext_xsx=new prc_cmp[wvl_nbr];
//   prc_cmp *sca_xsx=new prc_cmp[wvl_nbr];
//   prc_cmp *flx_planck_frc=new prc_cmp[wvl_nbr];
  prc_cmp *msv_1km=new prc_cmp[wvl_nbr]; // [frc] Emissivity of 1 km column
  prc_cmp *ext_fsh_ffc=new prc_cmp[wvl_nbr]; // [frc] Effective extinction efficiency
  prc_cmp *bck_fsh_ffc=new prc_cmp[wvl_nbr]; // [frc] Effective backscattering efficiency
  prc_cmp *sca_fsh_ffc=new prc_cmp[wvl_nbr]; // [frc] Effective scattering efficiency
  prc_cmp *abs_fsh_ffc=new prc_cmp[wvl_nbr]; // [frc] Effective absorption efficiency

  /* Create arrays to contain sub-band properties for each interval
     These arrays are temporary---used within CL1 and then discarded
     For parallelization over CL1, make these arrays automatic and allocated within CL1 if possible
     These arrays should be private within CL1, public within CL2 and CL3 */
  prc_cmp *bnd_ctr=new prc_cmp[bnd_nbr]; // [m] Wavelength at band center
  prc_cmp *bnd_dlt=new prc_cmp[bnd_nbr]; // [m] Bandwidth
  prc_cmp *bnd_min_wvn=new prc_cmp[bnd_nbr]; // [cm-1] Minimum wavenumber in band
  prc_cmp *bnd_max_wvn=new prc_cmp[bnd_nbr]; // [cm-1] Maximum wavenumber in band
  prc_cmp *bnd_min=new prc_cmp[bnd_nbr]; // [m] Minimum wavelength in band
  prc_cmp *bnd_max=new prc_cmp[bnd_nbr]; // [m] Maximum wavelength in band
  prc_cmp *flx_slr_frc_bnd=new prc_cmp[bnd_nbr]; // [frc] Fraction of solar flux in band
  prc_cmp *flx_IR_frc_bnd=new prc_cmp[bnd_nbr]; // [frc] Fraction of infrared flux in band
  std::complex<prc_cmp> *idx_rfr_cor_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of core
  std::complex<prc_cmp> *idx_rfr_ffc_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Effective refractive index of particle
  std::complex<prc_cmp> *idx_rfr_mdm_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of medium
  std::complex<prc_cmp> *idx_rfr_mnt_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of mantle
  std::complex<prc_cmp> *idx_rfr_mtx_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of matrix
  std::complex<prc_cmp> *idx_rfr_ncl_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of inclusion
  std::complex<prc_cmp> *idx_rfr_prt_bnd=new std::complex<prc_cmp>[bnd_nbr]; // [frc] Refractive index of particle

  /* These sz_nbr arrays are shared diagnostics accumulated within CL3 
     They are only updated within critical regions */
  prc_cmp *abs_fsh=new prc_cmp[sz_nbr]; // [frc] Absorption efficiency
  prc_cmp *asm_prm_fsh=new prc_cmp[sz_nbr]; // [frc] Asymmetry parameter
  prc_cmp *ext_fsh=new prc_cmp[sz_nbr]; // [frc] Extinction efficiency
  prc_cmp *bck_fsh=new prc_cmp[sz_nbr]; // [frc] Backscattering efficiency
  prc_cmp *sca_fsh=new prc_cmp[sz_nbr]; // [frc] Scattering efficiency
  prc_cmp *ss_alb_fsh=new prc_cmp[sz_nbr]; // [frc] Single scattering albedo
  prc_cmp *ss_co_alb_fsh=new prc_cmp[sz_nbr]; // [frc] Single scattering co-albedo
  prc_cmp *ext_spc=new prc_cmp[sz_nbr]; // [m2 kg-1] Specific extinction
  prc_cmp *bck_spc=new prc_cmp[sz_nbr]; // [m2 kg-1] Specific backscattering
  prc_cmp *sca_spc=new prc_cmp[sz_nbr]; // [m2 kg-1] Specific scattering
  prc_cmp *abs_spc=new prc_cmp[sz_nbr]; // [m2 kg-1] Specific absorption

  /* There are 2*ngl_nbr-1 values of ngl 
     fxm: Change ngl_nbr to ngl_plr_hms_nbr and 2*ngl_nbr-1 to ngl_plr_sph_nbr?
     ngl_nbr is number of polar angles in one hemisphere [0 <= theta <= pi/2]
     2*ngl_nbr-1 is number of polar angles in sphere [0 <= theta <= pi]
     Store these values in ngl[0..2*ngl_nbr-2] = ngl[0..2*(ngl_nbr-1)]
     2*ngl_nbr is size of amplitude scattering matrix ("S Matrix") for BoH83
     (Recall that azimuthal angle spans  [0 <= phi <= 2*pi])
     To compute and print results every degree, set ngl_nbr=91 so 2*ngl_nbr-1=181
     Then ngl[0]=0, ngl[1]=1, ..., ngl[180]=180 [dgr]

     Best angle choice for Legendre expansion coefficients is not regular
     Wis771 shows that Gauss-Lobatto quadrature is most efficient for highly
     asymmetric phase functions, e.g., cloud droplet scattering.
     Thus, default choice of angles should be quadrature angles
     Use regularly spaced human-readable angles only for self-tests */
  prc_cmp *ngl=new prc_cmp[static_cast<unsigned int>(2*ngl_nbr-1)]; // [rdn] Scattering angle
  prc_cmp *ngl_dgr=new prc_cmp[static_cast<unsigned int>(2*ngl_nbr-1)]; // [dgr] Angle degrees
  prc_cmp *ngl_dlt=new prc_cmp[static_cast<unsigned int>(2*ngl_nbr-1)]; // [rdn] Width of angle bin
  prc_cmp *ngl_wgt=new prc_cmp[static_cast<unsigned int>(2*ngl_nbr-1)]; // [rdn] Weight of angle bin
  rcd+=ngl_grd_get // [fnc] Create angular grid
    (ngl_nbr, // I [nbr] Number of angles
     ngl_sng, // I [sng] Angle grid type
     ngl, // O [rdn] Angle
     ngl_dgr, // O [dgr] Angle degrees
     ngl_dlt, // O [rdn] Width of angle bin
     ngl_wgt); // O [rdn] Weight of angle bin
  
  prc_cmp *phz_fnc_crr; // [sr-1] Phase function at current size, wavelength
  prc_cmp *plz_crr; // [frc] Degree of linear polarization at current size, wavelength
  prc_cmp *phz_fnc_dgn=new prc_cmp[2*ngl_nbr-1]; // [sr-1] Phase function at diagnostic size, wavelength
  a2d_cls<prc_cmp> phz_fnc_ffc(wvl_nbr,2*ngl_nbr-1); // [sr-1] Effective phase function (weighted over size distribution and sub-bands)
  prc_cmp *plz_dgn=new prc_cmp[2*ngl_nbr-1]; // [frc] Degree of linear polarization at diagnostic size, wavelength
  // Scalar initialization
  phz_fnc_ffc=0.0; // [sr-1] Effective phase function (weighted over size distribution and sub-bands)

  // Initialize diagnostic integrands that dependent on size grid
  for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    ext_fsh[sz_idx]=0.0; // [frc] Extinction efficiency
    bck_fsh[sz_idx]=0.0; // [frc] Backscattering efficiency
    sca_fsh[sz_idx]=0.0; // [frc] Scattering efficiency
    abs_fsh[sz_idx]=0.0; // [frc] Absorption efficiency
    asm_prm_fsh[sz_idx]=0.0; // [frc] Asymmetry parameter
  } // end loop over sz
  
  // Answers from Mie routines are always in double precision
  double abs_fct_MaS99_scl(CEWI_dbl); // [frc] Absorption enhancement for inclusions in weakly-absorbing spheres, scalar
  double asm_prm_scl(CEWI_dbl); // [frc] Asymmetry parameter, scalar
  double bck_hms_scl(CEWI_dbl); // [frc] Hemispheric backscatter, scalar
  double q_abs(CEWI_dbl); // [frc] Absorption eficiency
  double q_bck(CEWI_dbl); // [frc] Backscattering eficiency
  double q_ext(CEWI_dbl); // [frc] Extinction efficiency
  double q_sct(CEWI_dbl); // [frc] Scattering efficiency
  double spk_val(CEWI_dbl); // [frc] Mie coefficient of spike
 // Work variables
  double wgt_xsa; // [m2 m-3] Weight for cross-sectional area and spectral flux
  double wgt_vlm; // [m3 m-3] Weight of volume and spectral flux
  const long prg_mtr_rsn(100); // [nbr] Number of loops between progress meter updates

  /* Computation loop (CL) structure:
     Shared-memory (OpenMP) imeplementataion should sit over CL2 or CL3
     Typically, wvl_nbr ~ O(100) sz_nbr ~ O(100), bnd_nbr ~ O(100000)

     CL1: Outermost loop wvl_idx=[0..wvl_nbr-1]
     Case 1: OpenMP Shared Memory creates new thread for each wavelength
     Typically, there are 1--100 wavelengths
     Parallel region would be large and threads last entire length of program
     May require substantially more memory usage?
     Resonance runs often set wvl_nbr = 1, which would destroy parallelism

     CL2: Middle loop bnd_idx=[0..bnd_nbr-1]
     Case 2: OpenMP Shared Memory creates new thread for each band
     Typically, there are 10--100000 bands 
     Since bnd_nbr >> wvl_nbr, parallel region is large and threads
     will last virtually the entire execution time
     
     CL3: Innermost loop sz_idx=[0..sz_nbr-1]
     Case 3: OpenMP Shared Memory creates new thread for each size 
     Typically, there are 10--100 size bins, so parallel region is small
     With so few sizes, the overhead of creating/destroying threads is high

     Initial OpenMP implementation was Case 3, parallelize over size loop CL3 
     This led to modest benchmark improvements with typical thirty sizes
     Second OpenMP implementation was Case 2, parallelize over band loop CL2
     This dramatically improved throughput by increasing thread lifetime

     Initial values of private variables in parallel regions is indeterminate
     To copy-construct private variables from first value, use firstprivate() not private() */
  
  // NB: CL1 Parallelization is NOT READY YET
#ifdef _OPENMP // OpenMP
#ifdef PARALLELIZE_OVER_CL1
#pragma omp parallel for default(none) firstprivate(rcd) private(abs_fct_MaS99_scl,asm_prm_scl,bck_hms_scl,bnd_idx,ngl_idx,phz_fnc_crr,plz_crr,q_abs,q_bck,q_ext,q_sct,spk_val,sz_idx,wgt_xsa,wvl_idx) shared(abs_cff_vlm,abs_fct_MaS99,abs_fsh,abs_fsh_ffc,abs_ncl_wk_mdm_flg,asm_prm,asm_prm_fsh,bch_flg,bck_cff_vlm,bck_fsh,bck_fsh_ffc,bnd_ctr,bnd_nbr,cnc_sph,coat_flg,dbg_io,dbg_lvl,dbg_off,dbg_old,dbg_scl,dmn_frc,dns_prt,ext_cff_vlm,ext_fsh,ext_fsh_ffc,ffc_mdm_typ,flx_wgt_frc,flx_slr_frc_bnd,flx_IR_frc_bnd,flx_wgt_frc_bnd,idx_rfr_cor_bnd,idx_rfr_ffc_bnd,idx_rfr_mdm_bnd,idx_rfr_mnt_bnd,idx_rfr_mtx_bnd,idx_rfr_ncl_bnd,idx_rfr_prt_bnd,lgn_nbr,lng_foo,mie_flg,ngl,ngl_dlt,ngl_nbr,phz_fnc_ffc,phz_fnc_dgn,plz_dgn,prg_mtr_rsn,rds_cor,rds_mnt,sbr_nm,sca_cff_vlm,sca_fsh,sca_fsh_ffc,slf_tst_flg,slf_tst_typ,slv_sng,std::cerr,sz_ctr_sph,sz_idx_dbg,sz_nbr,sz_prm_rsn_usr_spc,tst_sng,vlm_sph,vlm_frc_ncl,wgt_vlm,wvl_bnd_sz_idx,wvl_bnd_sz_nbr,wvl_idx_dbg,wvl_nbr,xsa_sph)
#endif // !PARALLELIZE_OVER_CL1
#endif // !_OpenMP
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){ // CL1: Start outermost loop
    
    // Divide output wavelength bins into bnd_nbr regular spectral subintervals called bands
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      bnd_dlt[bnd_idx]=(wvl_max[wvl_idx]-wvl_min[wvl_idx])/bnd_nbr; // [m] Positive definite
      bnd_min[bnd_idx]=wvl_min[wvl_idx]+bnd_idx*bnd_dlt[bnd_idx]; // [m] Minimum wavelength in band
      bnd_max[bnd_idx]=wvl_min[wvl_idx]+(bnd_idx+1)*bnd_dlt[bnd_idx]; // [m] Maximum wavelength in band
      bnd_ctr[bnd_idx]=0.5*(bnd_max[bnd_idx]+bnd_min[bnd_idx]); // [m] Wavelength at band center
      bnd_min_wvn[bnd_idx]=1.0/(100.0*bnd_max[bnd_idx]); // [m] -> [cm-1] Minimum wavenumber in band
      bnd_max_wvn[bnd_idx]=1.0/(100.0*bnd_min[bnd_idx]); // [m] -> [cm-1] Maximum wavenumber in band
    } // end loop over bnd
    
    rcd+=flx_slr_src.flx_frc_get(bnd_min,bnd_max,bnd_nbr,flx_slr_frc_bnd); // [frc] Fraction of solar flux in band
    rcd+=spc_bbd.flx_frc_get(bnd_min,bnd_max,bnd_nbr,flx_IR_frc_bnd); // [frc] Fraction of infrared flux in band

    // Assign refractive indices
    rcd+=cmp_cor.idx_rfr_get(bnd_ctr,idx_rfr_cor_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of core
    rcd+=cmp_mdm.idx_rfr_get(bnd_ctr,idx_rfr_mdm_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of medium
    rcd+=cmp_mnt.idx_rfr_get(bnd_ctr,idx_rfr_mnt_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of mantle
    rcd+=cmp_mtx.idx_rfr_get(bnd_ctr,idx_rfr_mtx_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of matrix
    rcd+=cmp_ncl.idx_rfr_get(bnd_ctr,idx_rfr_ncl_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of inclusion
    rcd+=cmp_prt.idx_rfr_get(bnd_ctr,idx_rfr_prt_bnd,bnd_nbr,wrn_ntp_flg); // [frc] Refractive index of particle
    
    // Compute only user-selected effective medium approximations on band grid in CL2
    rcd+=idx_rfr_ffc_get // [fnc] Compute effective optical constants of composites
      (ncl_nbr, // I [nbr] Number of inclusions
       bnd_nbr, // I [nbr] Number of wavelengths (or, in this case, bands)
       ffc_mdm_typ, // I [enm] Effective medium type
       vlm_frc_ncl_mlt, // I [frc] Volume fraction(s) of inclusion(s)
       idx_rfr_cor_bnd, // I [frc] Refractive index of core
       idx_rfr_mdm_bnd, // I [frc] Refractive index of medium
       idx_rfr_mnt_bnd, // I [frc] Refractive index of mantle
       idx_rfr_mtx_bnd, // I [frc] Refractive index of matrix
       idx_rfr_ncl_bnd, // I [frc] Refractive index of inclusion
       idx_rfr_prt_bnd, // I [frc] Refractive index of particle
       idx_rfr_ffc_bnd); // O [frc] Effective refractive index of particle
    
    // Assign appropriate type weight to current spectral interval
    const prc_cmp *flx_wgt_frc_bnd(CEWI_NULL); // [frc] Radiant flux weighting factor
    prc_cmp flx_wgt_frc(CEWI_cpv); // [frc] Radiant flux weighting factor
    if(spc_wgt_sng == "Default"){
      // Default weight is spectral flux (Planck or Rayleigh weighting)
      if(bnd_ctr[0] < bnd_SW_LW && bnd_ctr[bnd_nbr-1] > bnd_SW_LW){
	/* Spectral region straddles solar and infrared */
	flx_wgt_frc_bnd=flx_IR_frc_bnd; // [frc] Radiant flux weighting factor
	flx_wgt_frc=flx_IR_frc[wvl_idx]; // [frc] Radiant flux weighting factor
	wrn_prn(prg_nm,sbr_nm,"Using LW weights in partly solar wavelength band:");
	std::cerr << "bnd_ctr[0] = " << bnd_ctr[0]*1.0e6 << " um, " << "bnd_ctr[" << bnd_nbr-1 << "] = " << bnd_ctr[bnd_nbr-1]*1.0e6 << " um" << std::endl;
	std::cerr << "Solar--Longwave boundary is " << bnd_SW_LW*1.0e6 << " um" << std::endl;
	std::cerr << "Straddled by wvl[" << wvl_idx << "] = " << wvl_min[wvl_idx]*1.0e6 << "--" << wvl_max[wvl_idx]*1.0e6 << " um = " << wvn_min[wvl_idx] << "--" << wvn_max[wvl_idx] << " cm-1" << std::endl;
	std::cerr << "idx\tbnd_min\tbnd_max\twvn_min\twvn_max\tSW_wgt\tLW_wgt" << std::endl;
	std::cerr << "idx\tum\tum\tcm-1\tcm-1\t\t" << std::endl;
	for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
	  std::cerr << bnd_idx << "\t" << bnd_min[bnd_idx]*1.0e6 << "\t" << bnd_max[bnd_idx]*1.0e6 << "\t" << bnd_min_wvn[bnd_idx] << "\t" << bnd_max_wvn[bnd_idx] << "\t" << flx_slr_frc_bnd[bnd_idx] << "\t" << flx_IR_frc_bnd[bnd_idx] << std::endl;
	} // !bnd_idx
      }else if(wvl_ctr[wvl_idx] < bnd_SW_LW){
	/* Spectral region completely solar */
	flx_wgt_frc_bnd=flx_slr_frc_bnd; // [frc] Radiant flux weighting factor
	flx_wgt_frc=flx_slr_frc[wvl_idx]; // [frc] Radiant flux weighting factor
	if(flx_wgt_frc == 0.0){
	  /* flx_wgt_frc == 0.0 for bands blue-ward of input solar spectral flux data, e.g., < 0.2 um
	     Set equal band weights flx_wgt_frc_bnd to avoid NaN's caused by 0.0 in denominator */
	  for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++) flx_slr_frc_bnd[bnd_idx]=1.0; // [frc] Fraction of solar flux in band
	  // Weight over wavelength bin is sum of bnd_nbr band weights to give band-equal weighting
	  flx_wgt_frc=bnd_nbr*1.0; // [frc] Radiant flux weighting factor
	} // !wvl_ctr
      }else{
	/* Spectral region completely infrared */
	flx_wgt_frc_bnd=flx_IR_frc_bnd; // [frc] Radiant flux weighting factor
	flx_wgt_frc=flx_IR_frc[wvl_idx]; // [frc] Radiant flux weighting factor
      } // !wvl_ctr
    }else if(spc_wgt_sng == "Reflectance"){ // endif spc_wgt_sng == "Default"
      err_prn(prg_nm,sbr_nm,"spc_wgt_sng == Reflectance weight not implemented yet");
    }else{ // !spc_wgt_sng == "Reflectance"
      err_prn(prg_nm,sbr_nm,"Unknown spc_wgt_sng");
    } // !spc_wgt_sng == "Default"
    assert(flx_wgt_frc > 0.0 || wvl_nbr*bnd_nbr == 1);
    
    // Weight refractive indices over current spectral interval
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){
      if(idx_rfr_prt_bnd[bnd_idx].imag() < 0.0) err_prn(prg_nm,sbr_nm,"idx_rfr_prt_bnd.imag() < 0.0. HINT: Imaginary part of refractive index must be positive");
      idx_rfr_ffc_wgt[wvl_idx]+=idx_rfr_ffc_bnd[bnd_idx]*flx_wgt_frc_bnd[bnd_idx]; // [frc] Spectral flux-weighted effective refractive index
      wvl_wgt[wvl_idx]+=bnd_ctr[bnd_idx]*flx_wgt_frc_bnd[bnd_idx]; // [m] Solar flux-weighted wavelength
    } // end loop over bnd

    // Normalize by total fraction to get appropriate spectral weighted average
    idx_rfr_ffc_wgt[wvl_idx]/=flx_wgt_frc; // [frc] Spectral flux-weighted effective refractive index
    wvl_wgt[wvl_idx]/=flx_wgt_frc; // [m] Solar flux-weighted wavelength
    
    /* Rationale: Weighting and computation loops

       Weight optical properties by fraction of solar or terrestrial flux in wavelength bin 
       Rationale is intuitive: we want to weight optical properties of each bin 
       by variation of refractive indices and solar flux within the bin. 
       Two sensible methods accomplish this---one fast, one slow
       
       Fast method is to first weight refractive indices by solar flux over bin
       It takes only O(wvl_nbr*bnd_nbr) computations to obtain wavelength-weighted refractive indices and wavelength for each wavelength bin
       Next, compute Mie properties for each aerosol size in each wavelength bin using weighted refractive indices and weighted wavelength 
       Total number of mie_sph_BoH83() calls for this fast method is wvl_nbr*sz_nbr
       
       Slow method is to perform a Mie calculation for each sub-band within each wavelength bin for each size
       Total number of mie_sph_BoH83() calls for this slow method is wvl_nbr*bnd_nbr*sz_nbr

       Same rationale holds for longwave radiation, but optical properties are weighted by Planck function of local temperature rather than solar flux */

    // Initialize integrands
    abs_cff_vlm[wvl_idx]=0.0; // [m2 m-3] Volume absorption coefficient
    abs_fct_MaS99[wvl_idx]=0.0; // [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
    abs_fsh_ffc[wvl_idx]=0.0; // [frc] Effective absorption efficiency
    asm_prm[wvl_idx]=0.0; // [frc] Asymmetry parameter
    bck_cff_vlm[wvl_idx]=0.0; // [m2 m-3] Volume backscattering coefficient
    bck_fsh_ffc[wvl_idx]=0.0; // [frc] Effective backscattering efficiency
    ext_cff_vlm[wvl_idx]=0.0; // [m2 m-3] Volume extinction coefficient
    ext_fsh_ffc[wvl_idx]=0.0; // [frc] Effective extinction efficiency
    sca_cff_vlm[wvl_idx]=0.0; // [m2 m-3] Volume scattering coefficient
    sca_fsh_ffc[wvl_idx]=0.0; // [frc] Effective scattering efficiency
    wgt_vlm=0.0; // [m3 m-3] Weight of volume and spectral flux

    // Default behavior: Parallelize over CL2 which is faster for typical resonance run geometries
#ifndef PARALLELIZE_OVER_CL3
#define PARALLELIZE_OVER_CL2 1
#endif // PARALLELIZE_OVER_CL3
#ifdef _OPENMP // OpenMP
#ifdef PARALLELIZE_OVER_CL2
#pragma omp parallel for default(none) firstprivate(rcd) private(abs_fct_MaS99_scl,asm_prm_scl,bck_hms_scl,bnd_idx,ngl_idx,phz_fnc_crr,plz_crr,q_abs,q_bck,q_ext,q_sct,spk_val,sz_idx,wgt_xsa) shared(abs_cff_vlm,abs_fct_MaS99,abs_fsh,abs_fsh_ffc,abs_ncl_wk_mdm_flg,asm_prm,asm_prm_fsh,bch_flg,bck_cff_vlm,bck_fsh,bck_fsh_ffc,bnd_ctr,bnd_nbr,cnc_sph,coat_flg,dbg_lvl,dmn_frc,dns_prt,ext_cff_vlm,ext_fsh,ext_fsh_ffc,ffc_mdm_typ,flx_wgt_frc,flx_slr_frc_bnd,flx_IR_frc_bnd,flx_wgt_frc_bnd,idx_rfr_cor_bnd,idx_rfr_ffc_bnd,idx_rfr_mdm_bnd,idx_rfr_mnt_bnd,idx_rfr_mtx_bnd,idx_rfr_ncl_bnd,idx_rfr_prt_bnd,lgn_nbr,lng_foo,mie_flg,ngl,ngl_dlt,ngl_nbr,phz_fnc_ffc,phz_fnc_dgn,plz_dgn,rds_cor,rds_mnt,sca_cff_vlm,sca_fsh,sca_fsh_ffc,slf_tst_flg,slf_tst_typ,slv_sng,std::cerr,sz_ctr_sph,sz_idx_dbg,sz_nbr,sz_prm_rsn_usr_spc,tst_sng,vlm_sph,vlm_frc_ncl,wgt_vlm,wvl_bnd_sz_idx,wvl_bnd_sz_nbr,wvl_idx,wvl_idx_dbg,wvl_nbr,xsa_sph)
#endif // !PARALLELIZE_OVER_CL2
#endif // !_OpenMP
    for(bnd_idx=0;bnd_idx<bnd_nbr;bnd_idx++){ // CL2: Start middle loop

#ifdef _OPENMP // OpenMP
#ifdef PARALLELIZE_OVER_CL3
#pragma omp parallel for default(none) firstprivate(rcd) private(abs_fct_MaS99_scl,asm_prm_scl,bck_hms_scl,ngl_idx,phz_fnc_crr,plz_crr,q_abs,q_bck,q_ext,q_sct,spk_val,sz_idx,wgt_xsa) shared(abs_cff_vlm,abs_fct_MaS99,abs_fsh,abs_fsh_ffc,abs_ncl_wk_mdm_flg,asm_prm,asm_prm_fsh,bch_flg,bck_cff_vlm,bck_fsh,bck_fsh_ffc,bnd_ctr,bnd_idx,bnd_nbr,cnc_sph,coat_flg,dbg_io,dbg_lvl,dbg_off,dbg_old,dbg_scl,dmn_frc,dns_prt,ext_cff_vlm,ext_fsh,ext_fsh_ffc,ffc_mdm_typ,flx_IR_frc_bnd,flx_slr_frc_bnd,flx_wgt_frc,flx_wgt_frc_bnd,idx_rfr_cor_bnd,idx_rfr_ffc_bnd,idx_rfr_mdm_bnd,idx_rfr_mnt_bnd,idx_rfr_mtx_bnd,idx_rfr_ncl_bnd,idx_rfr_prt_bnd,lgn_nbr,lng_foo,mie_flg,ngl,ngl_dlt,ngl_nbr,phz_fnc_ffc,phz_fnc_dgn,plz_dgn,prg_mtr_rsn,rds_cor,rds_mnt,sbr_nm,sca_cff_vlm,sca_fsh,sca_fsh_ffc,slf_tst_flg,slf_tst_typ,slv_sng,std::cerr,sz_ctr_sph,sz_idx_dbg,sz_nbr,sz_prm_rsn_usr_spc,tst_sng,vlm_sph,vlm_frc_ncl,wgt_vlm,wvl_bnd_sz_idx,wvl_bnd_sz_nbr,wvl_idx,wvl_idx_dbg,wvl_nbr,xsa_sph)
#endif // !PARALLELIZE_OVER_CL3
#endif // endif OpenMP
      for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){ // CL3: Start innermost loop
	// Allocate workspace variables needed within main loop
	/* Threads need private phase function diagnostic phz_fnc_crr and plz_crr
	   OpenMP would only create private copies of array pointers, not buffers */
	phz_fnc_crr=new prc_cmp[2*ngl_nbr-1]; // [sr-1] Phase function at current size, band
	plz_crr=new prc_cmp[2*ngl_nbr-1]; // [frc] Degree of linear polarization at current size, band

	if(mie_flg){

	  rcd+=mie_prc // [fnc] Mie processor
	    (abs_ncl_wk_mdm_flg, // I [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
	     coat_flg, // I [flg] Assume coated spheres
	     slf_tst_flg, // I [flg] Self-test flag
	     bnd_idx, // I [idx] Counting index for band
	     lgn_nbr, // I [nbr] Order of phase function Legendre expansion
	     ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
	     sz_idx, // I [idx] Counting index for size
	     wvl_idx, // I [idx] Counting index for wavelength
	     wvl_idx_dbg, // I [idx] Debugging wavelength bin
	     bnd_ctr[bnd_idx], // I [m] Wavelength at band center
	     dmn_frc, // I [frc] Fractal dimensionality of inclusions
	     rds_cor[sz_idx], // I [m] Radius of core
	     rds_mnt[sz_idx], // I [m] Radius of mantle
	     sz_ctr_sph[sz_idx], // I [m] Size at bin center
	     sz_prm_rsn_usr_spc, // I [m m-1] Size parameter resolution, user specified
	     idx_rfr_cor_bnd[bnd_idx], // I [frc] Refractive index of core
	     idx_rfr_ffc_bnd[bnd_idx], // I [frc] Effective refractive index of particle
	     idx_rfr_mdm_bnd[bnd_idx], // I [frc] Refractive index of medium
	     idx_rfr_mnt_bnd[bnd_idx], // I [frc] Refractive index of mantle
	     idx_rfr_mtx_bnd[bnd_idx], // I [frc] Refractive index of matrix
	     idx_rfr_ncl_bnd[bnd_idx], // I [frc] Refractive index of inclusion
	     idx_rfr_prt_bnd[bnd_idx], // I [frc] Refractive index of particle
	     slv_sng, // I [str] Mie solver to use (BoH83 or Wis79)
	     // Output
	     abs_fct_MaS99_scl, // O [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
	     asm_prm_scl, // O [frc] Asymmetry parameter
	     bck_hms_scl, // O [frc] Hemispheric backscatter
	     q_abs, // O [frc] Absorption efficiency
	     q_bck, // O [frc] Backscattering efficiency
	     q_ext, // O [frc] Extinction efficiency
	     q_sct, // O [frc] Scattering efficiency
	     spk_val, // O [frc] Mie coefficient of spike
	     ngl, // I [rdn] Scattering angle
	     ngl_dlt, // I [rdn] Width of angle bin
	     phz_fnc_crr, // O [sr-1] Phase function at current size, band
	     plz_crr); // O [frc] Degree of linear polarization at current size, band
	} // end if mie_flg	

	/* Test integration over size distribution  */
	if(tst_sng == "psd_ntg_dgn"){
	  /* Aggregating optical properties over a size distribution is error-prone
	     To test size distribution integration, set optical efficiency and flux weight to one
	     Then carry on with normal integration of optical properties over size distribution
	     Algorithm then computes geometric moments such as rds_xsa
	     Ideally, degenerate extinctions exactly equal analytic size distribution moments
	     When done correctly, the following equivalences will occur:
	     Let A = xsa_sph_rsl, V = vlm_sph_rsl
	     ext_cff_mss=A/V exactly (NB: must set density dns_prt=1.0)
	     ext_cff_vlm=A (NB: must set flx_wgt_frc=1.0 and cnc_nbr_sph_rsl_frc=1.0)
	     ext_fsh_ffc=1 (NB: must set flx_wgt_frc=1.0)
	     ext_cff_vlm/ext_cff_mss=V
	     Normalization factors flx_wgt_frc and cnc_nbr_sph_rsl_frc are computed outside CL3
	     and so must be individually set to unity.
	     This may make some other diagnosed quantities internally inconsistent
	     Hence only pay attention to degenerate extinctions in psd_ntg_dgn output files
	     
	     For lognormal distributions with cnc_nbr_anl_ttl=1.0, 
	     A=(pi/4)*dmt_nma^2*exp(2*ln(gsd_anl)^2)=pi*rds_nma^2*exp(2*ln(gsd_anl)^2)
	     V=(pi/6)*dmt_nma^3*exp(4.5*ln(gsd_anl)^2)=(4*pi/3)*rds_nma^3*exp(4.5*ln(gsd_anl)^2)
	     For dmt_nma=1.0 and gsd=e=exp(1)=2.718,
	     A=(pi/4)*1.0^2*exp(2*ln(e)^2)=(pi/4)*exp(2)
	     V=(pi/6)*1.0^3*exp(4.5*ln(e)^2)=(pi/4)*exp(4.5)
	     
	     Check that these relations hold.
	     Multiple wavelengths and bands stress test the PSD integration algorithm
	     pi=3.1415926532d
	     dmt_nma=1.0
	     gsd_anl_dfl=2.71828182846
	     mie --dbg=3 --bch --tst_sng=psd_ntg_dgn --wvl_nbr=1 --bnd_nbr=10 --sz_nbr=100 --sz_mnm=0.01 --sz_mxm=50.0 --dmt_nma=${dmt_nma} --psd_typ=lognormal --sz_grd=log --gsd_anl=${gsd_anl_dfl} --cnc_nbr_anl_dfl=1.0 ${DATA}/mie/mie.nc > foo 2>&1
	     ncap -D 1 -O -s "ext_cff_vlm_over_ext_cff_mss=ext_cff_vlm/ext_cff_mss" -s "xsa_xpc=(${pi}/4.0d)*dmt_nma^2*exp(2.0d*(ln(1.0*gsd_anl_dfl)^2))" -s "vlm_xpc=(${pi}/6.0d)*dmt_nma^3*exp(4.5d*(ln(1.0*gsd_anl_dfl)^2))" ${DATA}/mie/mie.nc ${DATA}/mie/foo.nc
	     ncks -C -v ext_cff_mss,ext_cff_vlm,ext_fsh_ffc,ext_cff_vlm_over_ext_cff_mss,xsa_rsl,xsa_xpc,vlm_rsl,vlm_xpc,cnc_nbr_rsl_frc,dmt_nma,gsd_anl_dfl ${DATA}/mie/foo.nc
	     Results should show that:
	     ext_cff_vlm_over_ext_cff_mss =~ vlm_rsl  ( =~ vlm_xpc when vlm_rsl_frc ~1)
	     ext_cff_vlm =~ xsa_rsl ( =~ xsa_xpc when xsa_rsl_frc ~1)
	     ext_fsh_ffc =~ 1 */
	  q_abs=1.0; // O [frc] Absorption efficiency
	  q_bck=1.0; // O [frc] Backscattering efficiency
	  q_ext=1.0; // O [frc] Extinction efficiency
	  q_sct=1.0; // O [frc] Scattering efficiency
	  asm_prm_scl=1.0; // O [frc] Asymmetry parameter
	  //	     cnc_nbr_sph_rsl_frc=1.0; // [frc] Resolved fraction of number concentration
	  dns_prt=1.0; // [kg m-3] Density of particle
	  flx_wgt_frc=bnd_nbr; // [frc] Radiant flux weighting factor
	  // flx_wgt_frc_bnd points to flx_slr_frc_bnd or flx_IR_frc_bnd so change those
	  flx_slr_frc_bnd[bnd_idx]=1.0; // [frc] Fraction of solar flux in band
	  flx_IR_frc_bnd[bnd_idx]=1.0; // [frc] Fraction of infrared flux in band
	} /* endif tst_sng=="psd_ntg_dgn" */
	
	/* wgt_xsa weights size distribution by cross-sectional area and spectral flux
	   wgt_xsa appears in integrand in numerator (in CL3) and so does not accumulate */
       	wgt_xsa=xsa_sph[sz_idx]*cnc_sph[sz_idx]*flx_wgt_frc_bnd[bnd_idx]; // [m2 m-3] Weight of cross-sectional area and spectral flux
	
	/* Atomic Updates and Critical Regions:
	   Operation in parallel regions are safe iff shared variables are on RHS
	   Shared variables on LHS of equation must be protected if multiple threads
	   might have simultaneous access
	   When this occurs results are indeterminate
	   Accumulating integrands is the perfect example */
#ifdef _OPENMP // OpenMP
#pragma omp atomic
#endif // !_OpenMP
	wvl_bnd_sz_idx++; // [idx] Counting index for wavelength/band/size loop
	// Monitor progress to completion
	if(dbg_lvl > dbg_off && !bch_flg && wvl_bnd_sz_idx%prg_mtr_rsn == 0) std::cerr << std::setw(6) << "\rwvl = " << wvl_idx << "/" << wvl_nbr << ", bnd = " << bnd_idx << "/" << bnd_nbr << ", sz = " << sz_idx << "/" << sz_nbr << ", count = " << wvl_bnd_sz_idx << "/" << wvl_bnd_sz_nbr << ": " << 100.0*wvl_bnd_sz_idx/wvl_bnd_sz_nbr << "% done";

	/* Test handling of shared variables exposed to multiple simultaneous writes
	   lng_foo should be updated atomically (or in critical region)
	   On exit, lng_foo equals wvl_nbr*bnd_nbr*sz_nbr iff no write contentions
	   Deviation of lng_foo from wvl_bnd_sz_nbr measures "molecularity"
	   ncks -C -v lng_foo,wvl_bnd_sz_nbr ${DATA}/mie/mie.nc */
	lng_foo++; // [nbr] Intrinsic long temporary variable
	
#ifdef _OPENMP // OpenMP
#pragma omp critical (CL3_answers_update_wvl_arrays)
	{ // begin OpenMP critical
#endif // !_OpenMP
	  /* wgt_vlm weights size distribution by volume and spectral flux
	     wgt_vlm is used once denominator integral completes (in CL1 after CL2) */
	  wgt_vlm+=vlm_sph[sz_idx]*cnc_sph[sz_idx]*flx_wgt_frc_bnd[bnd_idx]; // [m3 m-3] Weight of volume and spectral flux
	  
	  /* Weight optical efficiencies by cross-sectional area and spectral flux
	     Aerosol Optical Nomenclature (AON) Remarks #1:
	     *_cff_vlm coefficients are initially computed per unit volume air, subsequently
	     converted to and ultimately output as per unit volume aerosol (not air).
	     Extinction efficiencies q_* convert aerosol cross-sectional area [m2 aer] to
	     extinguished light cross-sectional area of air [m2 air].
	     One could define the extinguished light in [m2 shd] (shd means "shade") for
	     pedagogic reasons, but introducing another dimension like shade is unnecessary.
	     Shade, or extinguished light, has same dimensions as air which it partitions
	     Dimensionless q_* efficiencies are really [m2 air m-2 aer]=[m2 shd m-2 aer]
	     Initially *_cff_vlm coefficients are 
	     q_???*wgt_xsa = [m2 air m-2 aer] * [m2 aer m-3 air] = [m-1 air]
	     After we normalize *_cff_vlm coefficients by wgt_vlm [m3 aer m-3 air] then
	     *_cff_vlm coefficients are in [m2 air m-3 aer] = [m2 air m-3 aer]. 
	     These are the final dimensions for the volumetric coefficients *_cff_vlm
	     Optical depth is *_cff_vlm times aerosol volume path [m3 aer m-2 air]
	     Dividing volumetric coefficients by aerosol density [kg aer m-3 aer] 
	     yields mass-specific coefficients *_cff_mss [m-2 air kg-1 aer]. */
	  ext_cff_vlm[wvl_idx]+=q_ext*wgt_xsa; // [m2 m-3] Volume extinction coefficient
	  bck_cff_vlm[wvl_idx]+=q_bck*wgt_xsa; // [m2 m-3] Volume backscattering coefficient
	  sca_cff_vlm[wvl_idx]+=q_sct*wgt_xsa; // [m2 m-3] Volume scattering coefficient
	  abs_cff_vlm[wvl_idx]+=q_abs*wgt_xsa; // [m2 m-3] Volume absorption coefficient
	  asm_prm[wvl_idx]+=asm_prm_scl*q_sct*wgt_xsa; // [frc] Asymmetry parameter
	  abs_fct_MaS99[wvl_idx]+=abs_fct_MaS99_scl*q_sct*wgt_xsa; // [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
	  // Effective efficiencies are size-distribution integrated efficiencies at a given wavelength
	  ext_fsh_ffc[wvl_idx]+=q_ext*wgt_xsa; // [frc] Effective extinction efficiency
	  bck_fsh_ffc[wvl_idx]+=q_bck*wgt_xsa; // [frc] Effective backscattering efficiency
	  sca_fsh_ffc[wvl_idx]+=q_sct*wgt_xsa; // [frc] Effective scattering efficiency
	  abs_fsh_ffc[wvl_idx]+=q_abs*wgt_xsa; // [frc] Effective absorption efficiency
	  
	  // Weight contributions to effective phase function by scattering efficiency, cross-sectional area of scatterers, and spectral flux
	  /* See integration notes in phz_fnc.cc. Summary is:
	     Integratate intensities not probabilities
	     I integrate phase functions (not moments) over the size distribution
	     This differs from integrating moments over the size distribution
	     I think my approach is valid, but I am not yet sure */
	  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++) phz_fnc_ffc(wvl_idx,ngl_idx)+=phz_fnc_crr[ngl_idx]*q_sct*wgt_xsa; // [sr-1] Effective phase function (weighted over size distribution and sub-bands)
	  
#ifdef _OPENMP // OpenMP
	} // end OpenMP critical CL3_answers_update_wvl_arrays
#endif // !_OpenMP
	
	// Accumulate diagnostic integrands
	if(wvl_idx == wvl_idx_dbg){
#ifdef _OPENMP // OpenMP
#pragma omp critical (CL3_answers_update_sz_arrays)
	  { // begin OpenMP critical
#endif // !_OpenMP
	    ext_fsh[sz_idx]+=q_ext*flx_wgt_frc_bnd[bnd_idx]; // [frc] Extinction efficiency
	    bck_fsh[sz_idx]+=q_bck*flx_wgt_frc_bnd[bnd_idx]; // [frc] Backscattering efficiency
	    sca_fsh[sz_idx]+=q_sct*flx_wgt_frc_bnd[bnd_idx]; // [frc] Scattering efficiency
	    abs_fsh[sz_idx]+=q_abs*flx_wgt_frc_bnd[bnd_idx]; // [frc] Absorption efficiency
	    asm_prm_fsh[sz_idx]+=asm_prm_scl*flx_wgt_frc_bnd[bnd_idx]; // [frc] Asymmetry parameter
#ifdef _OPENMP // OpenMP
	  } // end OpenMP critical CL3_answers_update_sz_arrays
#endif // !_OpenMP
	} // end if debug wavelength
	
	// Archive phase function for diagnostic size and first sub-band of diagnostic wavelength
	if(sz_idx == sz_idx_dbg && bnd_idx == 0){
	  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
	    phz_fnc_dgn[ngl_idx]=phz_fnc_crr[ngl_idx]; // [sr-1] Phase function at diagnostic size, wavelength
	    plz_dgn[ngl_idx]=plz_crr[ngl_idx]; // [frc] Degree of linear polarization at diagnostic size, wavelength
	  } // end loop over ngl
	} // endif sz_dbg
	
	// Delete workspace variables needed within main loop but no longer
	delete []phz_fnc_crr; // [sr-1] Phase function at current size, band
	delete []plz_crr; // [frc] Degree of linear polarization at current size, band
	
      } // CL3: End innermost loop sz_idx=[0..sz_nbr-1]
      
    } // CL2: End middle loop bnd_idx=[0..bnd_nbr-1]
    
    /* Compute single scattering albedo, asymmetry parameter, and mass-intensive optical properties before normalizing  
       This reduces error propogation 
       Mass extinction coefficient Psi is ratio of computed extinction efficiency to computed mass density of aerosol */
    // Check for roundoff error
    if(mie_flg && sca_cff_vlm[wvl_idx] > ext_cff_vlm[wvl_idx]){
      wrn_prn(sbr_nm,"Numerical roundoff error caused sca_cff_vlm > ext_cff_vlm so that ss_alb > 1.0 for wvl_idx = "+nbr2sng(wvl_idx)+". Exact numbers are sca_cff_vlm = "+nbr2sng(sca_cff_vlm[wvl_idx],prc_cmp_dcm_plc)+", ext_cff_vlm = "+nbr2sng(ext_cff_vlm[wvl_idx],prc_cmp_dcm_plc)+", and ss_alb = "+nbr2sng(sca_cff_vlm[wvl_idx]/ext_cff_vlm[wvl_idx],prc_cmp_dcm_plc)+".");
      if(sca_cff_vlm[wvl_idx] <= one_plus_mie_cnv_eps*ext_cff_vlm[wvl_idx]){
	sca_cff_vlm[wvl_idx]=ext_cff_vlm[wvl_idx];
	wrn_prn(sbr_nm,"Error is beneath roundoff precision threshhold mie_cnv_eps = "+nbr2sng(mie_cnv_eps)+". Re-setting sca_cff_vlm=ext_cff_vlm so ss_alb=1.0.");
      }else{ // endif roundoff error
	assert(sca_cff_vlm[wvl_idx] <= ext_cff_vlm[wvl_idx]);
      } // endif error
    } // endif roundoff error      
    ss_alb[wvl_idx]=sca_cff_vlm[wvl_idx]/ext_cff_vlm[wvl_idx]; // [frc] Single scattering albedo 
    ss_co_alb[wvl_idx]=1.0-ss_alb[wvl_idx]; // [frc] Single scattering co-albedo 
    asm_prm[wvl_idx]/=sca_cff_vlm[wvl_idx]; // [frc] Asymmetry parameter
    abs_fct_MaS99[wvl_idx]/=sca_cff_vlm[wvl_idx]; // [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
    
    // Extinction, scattering, absorption, and backscattering per unit volume aerosol
    // See Aerosol Optical Nomenclature Remarks #1, #2 above and below, respectively
    ext_cff_vlm[wvl_idx]/=wgt_vlm; // [m2 m-3] Volume extinction coefficient
    bck_cff_vlm[wvl_idx]/=wgt_vlm; // [m2 m-3] Volume backscattering coefficient
    sca_cff_vlm[wvl_idx]/=wgt_vlm; // [m2 m-3] Volume scattering coefficient
    abs_cff_vlm[wvl_idx]/=wgt_vlm; // [m2 m-3] Volume absorption coefficient

    // Specific extinction, scattering, absorption, and backscattering
    ext_cff_mss[wvl_idx]=ext_cff_vlm[wvl_idx]/dns_prt; // [m2 kg-1] Mass extinction coefficient
    bck_cff_mss[wvl_idx]=bck_cff_vlm[wvl_idx]/dns_prt; // [m2 kg-1] Mass backscattering coefficient
    sca_cff_mss[wvl_idx]=sca_cff_vlm[wvl_idx]/dns_prt; // [m2 kg-1] Mass scattering coefficient
    abs_cff_mss[wvl_idx]=abs_cff_vlm[wvl_idx]/dns_prt; // [m2 kg-1] Mass absorption coefficient
    msv_1km[wvl_idx]=1.0-std::exp(-1.66*abs_cff_mss[wvl_idx]*1000.0*mss_sph_rsl); // [frc] Emissivity of 1 km column

    // Normalize intensive optical properties of coated spheres to core properties
    if(coat_nrm_by_core_flg){
      /* Normalize *_cff_vlm and *_cff_mss to core (not mixture) volume and mass
	 Intensive optical properties for coated spheres have, by default, dimensions
	 abs_cff_vlm [m2 air m-3 mixture]
	 abs_cff_mss [m2 air kg-1 mixture] 
	 Setting this flag changes renormalizes these optical properties, and
	 only these optical properties, to dimensions of 
	 abs_cff_vlm [m2 air m-3 core]
	 abs_cff_mss [m2 air kg-1 core] */
      abs_cff_vlm[wvl_idx]=abs_cff_vlm[wvl_idx]/vlm_frc_cor; // [m2 air m-3 mixture] * [m3 mixture m-3 core] = [m2 air m-3 core]
      abs_cff_mss[wvl_idx]=abs_cff_vlm[wvl_idx]/dns_cor; // [m2 air m-3 core] * [m3 core kg-1 core] = [m2 air kg-1 core]
      sca_cff_vlm[wvl_idx]=sca_cff_vlm[wvl_idx]/vlm_frc_cor; // [m2 air m-3 mixture] * [m3 mixture m-3 core] = [m2 air m-3 core]
      sca_cff_mss[wvl_idx]=sca_cff_vlm[wvl_idx]/dns_cor; // [m2 air m-3 core] * [m3 core kg-1 core] = [m2 air kg-1 core]
      ext_cff_vlm[wvl_idx]=ext_cff_vlm[wvl_idx]/vlm_frc_cor; // [m2 air m-3 mixture] * [m3 mixture m-3 core] = [m2 air m-3 core]
      ext_cff_mss[wvl_idx]=ext_cff_vlm[wvl_idx]/dns_cor; // [m2 air m-3 core] * [m3 core kg-1 core] = [m2 air kg-1 core]
    } // !coat_nrm_by_core_flg
    
    /* Normalize by resolved fraction of total concentration:
       This ensures that incompletely resolving distribution does not produce
       systematic error in optical properties. For example, say specified
       total number concentration of distribution is 100 cm-3. If, due to poor
       choice of quadrature points, or inadaquate number of quadrature points,
       the definite integral of n*dr is only 90 cm-3, then must divide computed
       optical properties by 90/100=0.9 to avoid biasing them. 
       Recall dN=n(r)*dr weights numerator of effective efficiency computation.
       Thus normalizing by _fractional_ error keeps the extinction coefficient,
       say, proportional to total number concentration. 
    
       Prior to 20060926 normalized *_cff_vlm by cnc_nbr_sph_rsl_frc*flx_wgt_frc
       However, now output *_cff_vlm in [m2 air m-3 aer] not [m2 air m-3 air] 
       Now require that *_cff_vlm and *_cff_mss are related solely by aerosol density
       Apply rsl vs. anl normalization factor to *_cff_dst instead */
    ext_fsh_ffc[wvl_idx]/=(xsa_sph_rsl*flx_wgt_frc); // [frc] Effective extinction efficiency
    bck_fsh_ffc[wvl_idx]/=(xsa_sph_rsl*flx_wgt_frc); // [frc] Effective backscattering efficiency
    sca_fsh_ffc[wvl_idx]/=(xsa_sph_rsl*flx_wgt_frc); // [frc] Effective scattering efficiency
    abs_fsh_ffc[wvl_idx]/=(xsa_sph_rsl*flx_wgt_frc); // [frc] Effective absorption efficiency
    
    // Normalize diagnostic integrands
    if(wvl_idx == wvl_idx_dbg){
      for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
	ext_fsh[sz_idx]/=flx_wgt_frc; // [frc] Extinction efficiency
	bck_fsh[sz_idx]/=flx_wgt_frc; // [frc] Backscattering efficiency
	sca_fsh[sz_idx]/=flx_wgt_frc; // [frc] Scattering efficiency
	abs_fsh[sz_idx]/=flx_wgt_frc; // [frc] Absorption efficiency
	asm_prm_fsh[sz_idx]/=flx_wgt_frc; // [frc] Asymmetry parameter
	ss_alb_fsh[sz_idx]=sca_fsh[sz_idx]/ext_fsh[sz_idx]; // [frc] Single scattering albedo
	ss_co_alb_fsh[sz_idx]=1.0-ss_alb_fsh[sz_idx]; // [frc] Single scattering co-albedo
      } // end loop over sz
    } // end if debug wavelength
    
      // Normalize effective phase function by summed scattering efficiency, cross-sectional area of scatterers, and spectral flux
    for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++) phz_fnc_ffc(wvl_idx,ngl_idx)/=sca_fsh_ffc[wvl_idx]*xsa_sph_rsl*flx_wgt_frc; // [sr-1] Effective phase function (weighted over size distribution and sub-bands)
    
  } // CL1: End outer loop wvl_idx=[0..wvl_nbr-1]
  if(dbg_lvl > dbg_off) std::cerr << std::endl;
  /* Here ends main loop over wavelength
     Apply user-specified adjustments to computed optical properties
     Then it will be safe to derive diagnostic quantities from optical properties */

  // Destroy arrays that contain sub-bands of each interval
  delete []bnd_ctr; // [m] Wavelength at band center
  delete []bnd_dlt; // [m] Bandwidth
  delete []bnd_min; // [m] Minimum wavelength in band
  delete []bnd_max; // [m] Maximum wavelength in band
  delete []flx_slr_frc_bnd; // [frc] Fraction of solar flux in band
  delete []flx_IR_frc_bnd; // [frc] Fraction of infrared flux in band
  delete []idx_rfr_cor_bnd; // [frc] Refractive index of core
  delete []idx_rfr_ffc_bnd; // [frc] Effective refractive index of particle
  delete []idx_rfr_mdm_bnd; // [frc] Refractive index of medium
  delete []idx_rfr_mnt_bnd; // [frc] Refractive index of mantle
  delete []idx_rfr_mtx_bnd; // [frc] Refractive index of matrix
  delete []idx_rfr_ncl_bnd; // [frc] Refractive index of inclusion
  delete []idx_rfr_prt_bnd; // [frc] Refractive index of particle

  // Implement user-specified scaling of aerosol optical properties
  if(fdg_flg){
    std::cout << "  WARNING: Extinction, absorption, and scattering scaled by user-specified factor of " << fdg_val << " in every band" << std::endl;
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      if(wvl_idx == fdg_idx || fdg_idx == 0){ // fxm: using idx = 0 to signify tuning all bands is kludgy
	sca_cff_mss[wvl_idx]*=fdg_val; // [m2 kg-1] Mass scattering coefficient
	ext_cff_mss[wvl_idx]*=fdg_val; // [m2 kg-1] Mass extinction coefficient
	abs_cff_mss[wvl_idx]=ext_cff_mss[wvl_idx]-sca_cff_mss[wvl_idx]; // [m2 kg-1] Mass absorption coefficient
	sca_cff_vlm[wvl_idx]*=fdg_val; // [m2 m-3] Volume scattering coefficient
	ext_cff_vlm[wvl_idx]*=fdg_val; // [m2 m-3] Volume extinction coefficient
	abs_cff_vlm[wvl_idx]=ext_cff_vlm[wvl_idx]-sca_cff_vlm[wvl_idx]; // [m2 m-3] Volume absorption coefficient
	ext_fsh_ffc[wvl_idx]*=fdg_val; // [frc] Effective extinction efficiency
	sca_fsh_ffc[wvl_idx]*=fdg_val; // [frc] Effective scattering efficiency
	abs_fsh_ffc[wvl_idx]=ext_fsh_ffc[wvl_idx]-sca_fsh_ffc[wvl_idx]; // [frc] Effective absorption efficiency
      } // endif
      if(wvl_idx == fdg_idx){
	for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
	  ext_fsh[sz_idx]*=fdg_val; // [frc] Extinction efficiency
	  sca_fsh[sz_idx]*=fdg_val; // [frc] Scattering efficiency
	  abs_fsh[sz_idx]=ext_fsh[sz_idx]-sca_fsh[sz_idx]; // [frc] Absorption efficiency
	} // end loop over sz
      } // end if fdg_idx
    } // end loop over wvl 
  } // endif fdg_flg

  // Implement user-specified single scattering albedo
  if(ss_alb_flg){
    std::cout << "WARNING: Single scattering properties scaled to user-specified value of omega = " << ss_alb_cmd_ln << " at wvl[" << wvl_idx_dbg << "] = " << wvl_ctr[wvl_idx_dbg]*1.0e6 << " um" << std::endl;
    prc_cmp crr_fct=ss_alb_cmd_ln/ss_alb[wvl_idx_dbg];
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      ss_alb[wvl_idx]=min_cpv(ss_alb[wvl_idx]*crr_fct,1.0); // [frc] Single scattering albedo
      ss_co_alb[wvl_idx]=1.0-ss_alb[wvl_idx]; // [frc] Single scattering co-albedo
      sca_cff_mss[wvl_idx]=min_cpv(sca_cff_mss[wvl_idx]*crr_fct,ext_cff_mss[wvl_idx]); // [m2 kg-1] Mass scattering coefficient
      abs_cff_mss[wvl_idx]=ext_cff_mss[wvl_idx]-sca_cff_mss[wvl_idx]; // [m2 kg-1] Mass absorption coefficient
      sca_cff_vlm[wvl_idx]=min_cpv(sca_cff_vlm[wvl_idx]*crr_fct,ext_cff_vlm[wvl_idx]); // [m2 m-3] Volume scattering coefficient
      abs_cff_vlm[wvl_idx]=ext_cff_vlm[wvl_idx]-sca_cff_vlm[wvl_idx]; // [m2 m-3] Volume absorption coefficient
      sca_fsh_ffc[wvl_idx]=min_cpv(sca_fsh_ffc[wvl_idx]*crr_fct,ext_fsh_ffc[wvl_idx]); // [frc] Effective scattering efficiency
      abs_fsh_ffc[wvl_idx]=ext_fsh_ffc[wvl_idx]-sca_fsh_ffc[wvl_idx]; // [frc] Effective absorption efficiency
      if(wvl_idx == wvl_idx_dbg){
	// NB: These scattering efficiencies are computed at wvl_idx_dbg, so scaling by ss_alb_cmd_ln is appropriate, but using wvl_idx_dbg AND sz_idx_dbg would be best.
	for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
	  ss_alb_fsh[sz_idx]=min_cpv(ss_alb_fsh[sz_idx]*crr_fct,1.0); // [frc] Single scattering albedo
	  ss_co_alb_fsh[sz_idx]=1.0-ss_alb_fsh[sz_idx]; // [frc] Single scattering co-albedo
	  sca_fsh[sz_idx]=min_cpv(sca_fsh[sz_idx]*crr_fct,ext_fsh[sz_idx]); // [frc] Scattering efficiency
	  abs_fsh[sz_idx]=ext_fsh[sz_idx]-sca_fsh[sz_idx]; // [frc] Absorption efficiency
	} // end loop over sz
      } // end if debug wavelength
    } // end loop over wvl 
  } // endif ss_alb_flg

  // Prognostic optical properties are now set and it is safe to derive purely diagnostic optical properties

  // Determine Angstrom exponent
  // fxm: Angstrom exponent should be computed in sub-gridscale loop
  for(wvl_idx=0;wvl_idx<wvl_nbr-1;wvl_idx++){ // NB: wvl_nbr-1
    ang_xpn[wvl_idx]=std::log(ext_cff_mss[wvl_idx+1]/ext_cff_mss[wvl_idx])/std::log(wvl_ctr[wvl_idx]/wvl_ctr[wvl_idx+1]); // [frc] Angstrom exponent
  } // end loop over wvl 
  ang_xpn[wvl_nbr-1]= (wvl_nbr > 1 ? ang_xpn[wvl_nbr-2] : 0); // [frc] Angstrom exponent

  /* Determine size-resolved specific optical properties at single wavelength
     ext_spc is specific extinction [m2 kg-1] for each size
     ext_cff_mss is same as ext_spc integrated over size distribution
     ext_spc is ext_cff_mss for a delta-function size distribution
     sca_spc, abs_spc, and bck_spc are defined analogously
     Now that efficiencies are known for each size, compute specific extinction */
  for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){
    ext_spc[sz_idx]=ext_fsh[sz_idx]*xsa_sph[sz_idx]/(vlm_sph[sz_idx]*dns_prt); // [m2 kg-1] Specific extinction
    bck_spc[sz_idx]=bck_fsh[sz_idx]*xsa_sph[sz_idx]/(vlm_sph[sz_idx]*dns_prt); // [m2 kg-1] Specific backscattering
    sca_spc[sz_idx]=sca_fsh[sz_idx]*xsa_sph[sz_idx]/(vlm_sph[sz_idx]*dns_prt); // [m2 kg-1] Specific scattering
    abs_spc[sz_idx]=ext_spc[sz_idx]-sca_spc[sz_idx]; // [m2 kg-1] Specific absorption
  } // end loop over sz
  
  // Compute broadband LW mass absorption coefficients, e.g., CCM:physics/cldems()/kabsl
  prc_cmp abs_cff_mss_bb_LW(0.0); // [m2 kg-1] Broadband longwave mass absorption coefficient
  prc_cmp flx_IR_frc_ttl(0.0); // [frc] Diagnostic total of flx_wgt_frc
  if(spc_bbd.flx_frc_get(wvl_mnm,wvl_mxm) > 0.0){
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      flx_IR_frc_ttl+=flx_IR_frc[wvl_idx]; // [frc] Diagnostic total of flux_wgt_frc
      abs_cff_mss_bb_LW+=abs_cff_mss[wvl_idx]*flx_IR_frc[wvl_idx]; // [m2 kg-1]
    } // end loop over wvl 
    // Normalize by contributing fraction in order to improve estimate of broadband total
    abs_cff_mss_bb_LW=abs_cff_mss_bb_LW/flx_IR_frc_ttl; // [m2 kg-1]
    // Safest to archive bogus abs_cff_mss_bb_LW when full LW spectrum was not specified, or when overlapping spectra, e.g., CAM_LW, cause fraction to exceed 1.0
    if(flx_IR_frc_ttl < 0.9 || flx_IR_frc_ttl > 1.0) abs_cff_mss_bb_LW=1.0e36; // [m2 kg-1]
  } // endif

  /* Nomenclature:
     ext_cff_mss [m2 kg-1] is extinction  per unit mass     of aerosol (not air)
     ext_cff_vlm [m2 m-3]  is extinction  per unit volume   of aerosol (not air)
     ext_cff_dst [m-1]     is extinction  per unit distance of air (not aerosol)
     Total aerosol optical depth = ext_cff_mss times aerosol mass path [kg m-2]
     Total aerosol optical depth = ext_cff_vlm times aerosol volume path [m3 m-2]
     Total aerosol optical depth = ext_cff_dst times aerosol layer thickness [m]

     Units of ext_cff_vlm and ext_cff_dst are _not_ equivalent!
     ext_cff_dst (unlike ext_cff_vlm) incorporates total aerosol number concentration
     meter in ext_cff_dst units [meter-1] refers to aerosol layer thickness
     ext_cff_vlm (unlike ext_cff_dst) is _normalized_ (to unit aerosol volume)
     ext_cff_vlm is extinction coefficient in meter2 meter-3 != meter-1 because 
     meter2 refers to air column cross-section and meter-3 refers to aerosol volume 
     Hence they are all "meters" though not the same kinds of "meters"
     ext_cff_vlm is an intensive property which depends only on normalized size distribution (not absolute size distribution or total number concentration)

     NB: Confusingly, many references call ext_cff_dst "volume extinction coefficient"
     This includes OPAC reference HKS98
     It is most common to measure ext_cff_dst in [km-1] not [m-1] */

  // Traditional distance extinction coefficients
  prc_cmp *abs_cff_dst=new prc_cmp[wvl_nbr]; // [m-1] Distance absorption coefficient
  prc_cmp *sca_cff_dst=new prc_cmp[wvl_nbr]; // [m-1] Distance scattering coefficient
  prc_cmp *ext_cff_dst=new prc_cmp[wvl_nbr]; // [m-1] Distance extinction coefficient
  prc_cmp *bck_cff_dst=new prc_cmp[wvl_nbr]; // [m-1] Distance backscattering coefficient
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    /* Multiply by resolved (rather than analytic) mass concentration so that distance 
       extinction coefficients improve as mass fraction resolved approaches unity
       fxm: normalize *_cff_dst by flx_wgt_frc? */
    abs_cff_dst[wvl_idx]=abs_cff_mss[wvl_idx]*mss_sph_rsl; // [m-1] Distance absorption coefficient
    sca_cff_dst[wvl_idx]=sca_cff_mss[wvl_idx]*mss_sph_rsl; // [m-1] Distance scattering coefficient
    ext_cff_dst[wvl_idx]=ext_cff_mss[wvl_idx]*mss_sph_rsl; // [m-1] Distance extinction coefficient
    bck_cff_dst[wvl_idx]=bck_cff_mss[wvl_idx]*mss_sph_rsl; // [m-1] Distance backscattering coefficient
  } // end loop over wvl
  
  // Derive visibility
  const prc_cmp wvl_vsb(0.55e-6); // [m] Wavelength for visibility diagnostics
  const long wvl_idx_vsb=vec_val2idx(wvl_ctr,wvl_nbr,wvl_vsb); // [idx] Visible wavelength bin
  prc_cmp vsb_vsb; // [m] Visibility at 0.55 um
  // fxm: Add Rayleigh scattering to visibility calculation
  vsb_vsb=3000.0/(ext_cff_dst[wvl_idx_vsb]); // [m] Visibility at 0.55 um HKS98 p. 838 (4f)

  // Output diagnostic microphysical optical and physical properties to console
  if(dbg_lvl > dbg_off){
    const int dat_nbr_per_ln(4); // [nbr] Number of datum per line
    const std::string wvl_tag_sng("_"+wvl_grd_sng); // [sng] Identifying tag
    std::string dat_sng; // [sng] Fortran data dimension statement
    std::string ftn_cmt; // [sng] Fortran comment string
    std::string ftn_dnt; // [sng] Fortran indentation string
    // Derived
    const std::string bnd_nbr_sng("bnd_nbr"+wvl_tag_sng); // [sng] Number of bands
    if(ftn_fxd_flg){
      ftn_cmt="c     "; // [sng] Fortran comment string
      ftn_dnt="      "; // [sng] Fortran indentation string
    }else{
      ftn_cmt=" ! "; // [sng] Fortran comment string
      ftn_dnt="     "; // [sng] Fortran indentation string
    } // endelse
    std::cout << ftn_cmt+"Data generated by " << usr_cpp << " at " << hst_cpp << " on " << time_bfr_srt;
    std::cout << ftn_cmt+"Program: " << prg_nm << " version " << vrs_cpp << std::endl;
    
    // Two dimensional microphysical optical properties initialized with data statements
    if((wvl_grd_sng == "CAM_SW" || wvl_grd_sng == "CAM_LW" || wvl_grd_sng == "GSFC_LW" || wvl_grd_sng == "TOMS") && dbg_lvl > dbg_off){
      std::cerr << ftn_cmt << "wvl_ctr"+wvl_tag_sng+": Nominal band center [m]" << std::endl;
      std::cerr << ftn_cmt << "real(r8) ext_cff_mss"+wvl_tag_sng+"(bnd_nbr"+wvl_tag_sng+",dst_nbr) ! [m2 kg-1] "+cmp_prt.dsc_get()+" aerosol mass extinction coefficient" << std::endl;
      std::cerr << ftn_cmt << "real(r8) bck_cff_mss"+wvl_tag_sng+"(bnd_nbr"+wvl_tag_sng+",dst_nbr) ! [m2 kg-1] "+cmp_prt.dsc_get()+" aerosol mass backscattering coefficient" << std::endl;
      std::cerr << ftn_cmt << "real(r8) abs_cff_mss"+wvl_tag_sng+"(bnd_nbr"+wvl_tag_sng+",dst_nbr) ! [m2 kg-1] "+cmp_prt.dsc_get()+" aerosol mass absorption coefficient" << std::endl;
      std::cerr << ftn_cmt << "real(r8) ss_alb"+wvl_tag_sng+"(bnd_nbr"+wvl_tag_sng+",dst_nbr) ! [frc] "+cmp_prt.dsc_get()+" aerosol single scattering albedo" << std::endl;
      std::cerr << ftn_cmt << "real(r8) asm_prm"+wvl_tag_sng+"(bnd_nbr"+wvl_tag_sng+",dst_nbr) ! [frc] "+cmp_prt.dsc_get()+" aerosol asymmetry parameter" << std::endl;
      std::cerr << "#if 0 /* Some Fortran compilers choke on very long comments */" << std::endl << "Command line: " << cmd_ln << std::endl << "#endif /* not 0 */" << std::endl;

      if(wvl_grd_sng == "CAM_SW"){
	// Avoid printing single scattering albedo = 1.0 (crashes in CCM SW code)
	for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
	  if(ss_alb[wvl_idx] == 1.0){
	    std::cout << "WARNING: Reducing printed single scattering albedo from 1.0 to " << ss_alb_bnd_CAM_SW << " at wvl[" << wvl_idx_dbg << "] = " << wvl_ctr[wvl_idx_dbg]*1.0e6 << " um" << std::endl;
	    ss_alb[wvl_idx]=ss_alb_bnd_CAM_SW; // [frc] Single scattering albedo
	  } // endif
	} // end loop over wvl 
      } // endif CAM_SW

      if(ftn_fxd_flg){
	dat_sng="      data (ext_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(bnd_idx,"+spc_idx_sng+"),bnd_idx=1,"+bnd_nbr_sng+") /";
	(void)f77_blk_dat_prn(dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,ext_cff_mss);
	dat_sng="      data (bck_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(bnd_idx,"+spc_idx_sng+"),bnd_idx=1,"+bnd_nbr_sng+") /";
	(void)f77_blk_dat_prn(dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,bck_cff_mss);
	dat_sng="      data (asm_prm_"+spc_abb_sng+wvl_tag_sng+"(bnd_idx,"+spc_idx_sng+"),bnd_idx=1,"+bnd_nbr_sng+") /";
	(void)f77_blk_dat_prn(dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,asm_prm);
	dat_sng="      data (ss_alb_"+spc_abb_sng+wvl_tag_sng+"(bnd_idx,"+spc_idx_sng+"),bnd_idx=1,"+bnd_nbr_sng+") /";
	(void)f77_blk_dat_prn(dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,ss_alb);
	dat_sng="      data (abs_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(bnd_idx,"+spc_idx_sng+"),bnd_idx=1,"+bnd_nbr_sng+") /";
	(void)f77_blk_dat_prn(dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,abs_cff_mss);
      }else{ // !ftn_fxd_flg
	std::cerr << ftn_dnt+"integer,parameter::bnd_nbr"+wvl_tag_sng+"=" << wvl_nbr << " ! [nbr] Number of bands" << std::endl;
	dat_sng=ftn_dnt+"real(r8),dimension(bnd_nbr"+wvl_tag_sng+"),parameter::wvl_ctr"+wvl_tag_sng+"=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,wvl_ctr);
	dat_sng=ftn_dnt+"data (ext_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(:,"+spc_idx_sng+")=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,ext_cff_mss);
	dat_sng=ftn_dnt+"data (bck_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(:,"+spc_idx_sng+")=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,bck_cff_mss);
	dat_sng=ftn_dnt+"data (asm_prm_"+spc_abb_sng+wvl_tag_sng+"(:,"+spc_idx_sng+")=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,asm_prm);
	dat_sng=ftn_dnt+"data (ss_alb_"+spc_abb_sng+wvl_tag_sng+"(:,"+spc_idx_sng+")=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,ss_alb);
	dat_sng=ftn_dnt+"data (abs_cff_mss_"+spc_abb_sng+wvl_tag_sng+"(:,"+spc_idx_sng+")=(/ ";
	(void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),wvl_nbr,abs_cff_mss);
      } // !ftn_fxd_flg
    } // endif CCM or GSFC or TOMS band structure
    
    // Output high resolution microphysical optical and physical properties to console
    if(hrz_flg){
      std::cerr << ftn_cmt << "ext_cff_mss_"+spc_abb_sng+"_hrz: "+cmp_prt.dsc_get()+" aerosol mass extinction coefficient [m2 kg-1]" << std::endl;
      std::cerr << ftn_cmt << "bck_cff_mss_"+spc_abb_sng+"_hrz: "+cmp_prt.dsc_get()+" aerosol mass backscattering coefficient [m2 kg-1]" << std::endl;
      std::cerr << ftn_cmt << "dmt_ctr_hrz: "+cmp_prt.dsc_get()+" aerosol diameter at bin center [m]" << std::endl;
      std::cerr << ftn_cmt << "scv_cff_pcp_nrm_"+spc_abb_sng+"_hrz: "+cmp_prt.dsc_get()+" aerosol scavenging coefficient, precipitation-normalized [m2 kg-1]" << std::endl;
      std::cerr << ftn_cmt << "dmt_grd_hrz: "+cmp_prt.dsc_get()+" diameter grid [m]" << std::endl;
      std::cerr << "#if 0 /* Some Fortran compilers choke on very long comments */" << std::endl << "Command line: " << cmd_ln << std::endl << "#endif /* not 0 */" << std::endl;
      std::cerr << ftn_dnt+"integer,parameter::dst_nbr_hrz=" << sz_nbr << " ! [nbr] Number of high resolution bins" << std::endl;
      dat_sng=ftn_dnt+"real(r8),dimension(dst_nbr_hrz),parameter::ext_cff_mss_"+spc_abb_sng+"_hrz=(/ ";
      (void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),sz_nbr,ext_spc);
      dat_sng=ftn_dnt+"real(r8),dimension(dst_nbr_hrz),parameter::bck_cff_mss_"+spc_abb_sng+"_hrz=(/ ";
      (void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),sz_nbr,bck_spc);
      dat_sng=ftn_dnt+"real(r8),dimension(dst_nbr_hrz),parameter::dmt_ctr_hrz=(/ ";
      (void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),sz_nbr,dmt_ctr);
      dat_sng=ftn_dnt+"real(r8),dimension(dst_nbr_hrz),parameter::scv_cff_pcp_nrm_"+spc_abb_sng+"_hrz=(/ ";
      (void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),sz_nbr,scv_cff_pcp_nrm);
      dat_sng=ftn_dnt+"real(r8),dimension(dst_nbr_hrz+1),parameter::dmt_grd_hrz=(/ ";
      (void)f90_prm_dat_prn(ftn_fxd_flg,dat_nbr_per_ln,dat_sng.c_str(),sz_nbr+1,dmt_grd);
    } // endif hrz_flg
  } // endif dbg
  
  // Unfortunately, netCDF does not support complex types natively
  // Split complex output into real and imaginary components for archival
  // Component indices
  prc_cmp *idx_rfr_air_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of air, imaginary component
  prc_cmp *idx_rfr_air_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of air, real component
  prc_cmp *idx_rfr_cor_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of core, imaginary component
  prc_cmp *idx_rfr_cor_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of core, real component
  prc_cmp *idx_rfr_mdm_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of medium, imaginary component
  prc_cmp *idx_rfr_mdm_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of medium, real component
  prc_cmp *idx_rfr_mnt_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of mantle, imaginary component
  prc_cmp *idx_rfr_mnt_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of mantle, real component
  prc_cmp *idx_rfr_mtx_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of matrix, imaginary component
  prc_cmp *idx_rfr_mtx_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of matrix, real component
  prc_cmp *idx_rfr_ncl_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of inclusion, imaginary component
  prc_cmp *idx_rfr_ncl_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of inclusion, real component
  prc_cmp *idx_rfr_prt_img=new prc_cmp[wvl_nbr]; // [frc] Refractive index of particle, imaginary component
  prc_cmp *idx_rfr_prt_rl=new prc_cmp[wvl_nbr]; // [frc] Refractive index of particle, real component
  prc_cmp *idx_rfr_ffc_rl=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index of particle, real component
  prc_cmp *idx_rfr_ffc_img=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index of particle, imaginary component
  prc_cmp *idx_rfr_ffc_brg_img=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, Bruggeman approximation, imaginary component
  prc_cmp *idx_rfr_ffc_brg_rl=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, Bruggeman approximation, real component
  prc_cmp *idx_rfr_ffc_mxg_img=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, Maxwell Garnett approximation, imaginary component
  prc_cmp *idx_rfr_ffc_mxg_rl=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, Maxwell Garnett approximation, real component
  prc_cmp *idx_rfr_ffc_pmr_img=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, partial molar refraction approximation, imaginary component
  prc_cmp *idx_rfr_ffc_pmr_rl=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, partial molar refraction approximation, real component
  prc_cmp *idx_rfr_ffc_vlw_img=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, volume-weighted approximation, imaginary component
  prc_cmp *idx_rfr_ffc_vlw_rl=new prc_cmp[wvl_nbr]; // [frc] Effective refractive index, volume-weighted approximation, real component
  prc_cmp *idx_rfr_ffc_wgt_img=new prc_cmp[wvl_nbr]; // [frc] Spectral flux-weighted effective refractive index, imaginary component
  prc_cmp *idx_rfr_ffc_wgt_rl=new prc_cmp[wvl_nbr]; // [frc] Spectral flux-weighted effective refractive index, real component
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    // Component indices
    idx_rfr_air_img[wvl_idx]=idx_rfr_air[wvl_idx].imag(); // [frc] Refractive index of air, imaginary component
    idx_rfr_air_rl[wvl_idx]=idx_rfr_air[wvl_idx].real(); // [frc] Refractive index of air, real component
    idx_rfr_cor_img[wvl_idx]=idx_rfr_cor[wvl_idx].imag(); // [frc] Refractive index of core, imaginary component
    idx_rfr_cor_rl[wvl_idx]=idx_rfr_cor[wvl_idx].real(); // [frc] Refractive index of core, real component
    idx_rfr_mdm_img[wvl_idx]=idx_rfr_mdm[wvl_idx].imag(); // [frc] Refractive index of medium, imaginary component
    idx_rfr_mdm_rl[wvl_idx]=idx_rfr_mdm[wvl_idx].real(); // [frc] Refractive index of medium, real component
    idx_rfr_mnt_img[wvl_idx]=idx_rfr_mnt[wvl_idx].imag(); // [frc] Refractive index of mantle, imaginary component
    idx_rfr_mnt_rl[wvl_idx]=idx_rfr_mnt[wvl_idx].real(); // [frc] Refractive index of mantle, real component
    idx_rfr_mtx_img[wvl_idx]=idx_rfr_mtx[wvl_idx].imag(); // [frc] Refractive index of matrix, imaginary component
    idx_rfr_mtx_rl[wvl_idx]=idx_rfr_mtx[wvl_idx].real(); // [frc] Refractive index of matrix, real component
    idx_rfr_ncl_img[wvl_idx]=idx_rfr_ncl[wvl_idx].imag(); // [frc] Refractive index of inclusion, imaginary component
    idx_rfr_ncl_rl[wvl_idx]=idx_rfr_ncl[wvl_idx].real(); // [frc] Refractive index of inclusion, real component
    idx_rfr_prt_img[wvl_idx]=idx_rfr_prt[wvl_idx].imag(); // [frc] Refractive index of particle, imaginary component
    idx_rfr_prt_rl[wvl_idx]=idx_rfr_prt[wvl_idx].real(); // [frc] Refractive index of particle, real component
    // Effective indices
    idx_rfr_ffc_rl[wvl_idx]=idx_rfr_ffc[wvl_idx].real(); // [frc] Effective refractive index of particle, real component
    idx_rfr_ffc_img[wvl_idx]=idx_rfr_ffc[wvl_idx].imag(); // [frc] Effective refractive index of particle, imaginary component
    idx_rfr_ffc_brg_img[wvl_idx]=idx_rfr_ffc_brg[wvl_idx].imag(); // [frc] Effective refractive index, Bruggeman approximation, imaginary component
    idx_rfr_ffc_brg_rl[wvl_idx]=idx_rfr_ffc_brg[wvl_idx].real(); // [frc] Effective refractive index, Bruggeman approximation, real component
    idx_rfr_ffc_mxg_img[wvl_idx]=idx_rfr_ffc_mxg[wvl_idx].imag(); // [frc] Effective refractive index, Maxwell Garnett approximation, imaginary component
    idx_rfr_ffc_mxg_rl[wvl_idx]=idx_rfr_ffc_mxg[wvl_idx].real(); // [frc] Effective refractive index, Maxwell Garnett approximation, real component
    idx_rfr_ffc_pmr_img[wvl_idx]=idx_rfr_ffc_pmr[wvl_idx].imag(); // [frc] Effective refractive index, partial molar refraction approximation, imaginary component
    idx_rfr_ffc_pmr_rl[wvl_idx]=idx_rfr_ffc_pmr[wvl_idx].real(); // [frc] Effective refractive index, partial molar refraction approximation, real component
    idx_rfr_ffc_vlw_img[wvl_idx]=idx_rfr_ffc_vlw[wvl_idx].imag(); // [frc] Effective refractive index, volume-weighted approximation, imaginary component
    idx_rfr_ffc_vlw_rl[wvl_idx]=idx_rfr_ffc_vlw[wvl_idx].real(); // [frc] Effective refractive index, volume-weighted approximation, real component
    idx_rfr_ffc_wgt_img[wvl_idx]=idx_rfr_ffc_wgt[wvl_idx].imag(); // [frc] Spectral flux-weighted effective refractive index, imaginary component
    idx_rfr_ffc_wgt_rl[wvl_idx]=idx_rfr_ffc_wgt[wvl_idx].real(); // [frc] Spectral flux-weighted effective refractive index, real component
  } // end loop over wvl 

  // Open output file
  int nccreate_mode(NC_CLOBBER); // [enm] Mode flag for nco_create() call
#ifdef ENABLE_NETCDF4
  if(fl_out_fmt==NCO_FORMAT_UNDEFINED) fl_out_fmt=NC_FORMAT_NETCDF4;
  if(fl_out_fmt == NC_FORMAT_64BIT){
    nccreate_mode|=NC_64BIT_OFFSET;
  }else if(fl_out_fmt == NC_FORMAT_NETCDF4){
    nccreate_mode|=NC_NETCDF4;
  }else if(fl_out_fmt == NC_FORMAT_NETCDF4_CLASSIC){
    nccreate_mode|=NC_NETCDF4|NC_CLASSIC_MODEL;
  } /* end else fl_out_fmt */
#else // !ENABLE_NETCDF4
  if(fl_out_fmt==NCO_FORMAT_UNDEFINED) fl_out_fmt=NC_FORMAT_CLASSIC;
  if(fl_out_fmt == NC_FORMAT_CLASSIC) nccreate_mode+=0; // CEWI
#endif // !ENABLE_NETCDF4
  const int nc_out(nco_create(fl_out,nccreate_mode)); 
  const nc_type nco_xtyp(nco_get_xtype(static_cast<prc_cmp>(1.0))); // [enm] External netCDF type
  if(dbg_lvl > dbg_off || tst_sng == "nc"){
    std::cout << prg_nm << ": INFO External netCDF type of prc_cmp variables will be " << nco_typ_sng(nco_xtyp) << std::endl;
    std::cout << prg_nm << ": INFO Record dimension will be " << (dmn_rcd == "" ? "non-existent" : dmn_rcd) << std::endl;
  } // endif dbg
  
  // Create dimensions
  const int sz_dmn(nco_def_dim(nc_out,static_cast<std::string>("sz"),dmn_rcd == "sz" ? NC_UNLIMITED : sz_nbr)); // [dmn] Size dimension
  const int psd_dmn(nco_def_dim(nc_out,static_cast<std::string>("psd"),dmn_rcd == "psd" ? NC_UNLIMITED : psd_nbr)); // [dmn] Particle size distribution dimension
  const int dsd_dmn(nco_def_dim(nc_out,static_cast<std::string>("dsd_sz"),dmn_rcd == "dsd_sz" ? NC_UNLIMITED : dsd_nbr)); // [dmn] Droplet size distribution dimension
  const int sz_grd_dmn(nco_def_dim(nc_out,static_cast<std::string>("sz_grd"),dmn_rcd == "sz_grd" ? NC_UNLIMITED : sz_nbr+1)); // [dmn] Size grid dimension
  const int wvl_dmn(nco_def_dim(nc_out,static_cast<std::string>("wvl"),dmn_rcd == "wvl" ? NC_UNLIMITED : wvl_nbr)); // [dmn] Wavelength dimension
  const int wvl_grd_dmn(nco_def_dim(nc_out,static_cast<std::string>("wvl_grd"),dmn_rcd == "wvl_grd" ? NC_UNLIMITED : wvl_nbr+1)); // [dmn] Wavelength grid dimension
  
  // Derive pointers to dimensions
  const int *dmn_dsd(&dsd_dmn); // [dmn] Pointer to droplet size distribution dimension
  const int *dmn_psd(&psd_dmn); // [dmn] Pointer to particle size distribution dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars, should not be accessed
  const int *dmn_sz(&sz_dmn); // [dmn] Pointer to size dimension
  const int *dmn_sz_grd(&sz_grd_dmn); // [dmn] Pointer to size grid dimension
  const int *dmn_wvl(&wvl_dmn); // [dmn] Pointer to wavelength dimension
  const int *dmn_wvl_grd(&wvl_grd_dmn); // [dmn] Pointer to wavelength grid dimension
  const int dmn_dsd_sz[2]={dsd_dmn,sz_dmn}; // [dmn] 

  if(dbg_lvl > dbg_off) std::cout << "Defined " << nco_inq_ndims(nc_out) << " dimensions in "+sbr_nm << std::endl;

  // Global attributes
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("history"),time_bfr_srt.substr(0,time_bfr_srt.size()-1)+": "+cmd_ln);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("CVS_Id"),CVS_Id);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("CVS_Name"),CVS_Name);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("solar_flux_source"),flx_slr_src.dsc_get());
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("aerosol_type"),cmp_sng_prt);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("aerosol_long_name"),cmp_prt.dsc_get());
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("aerosol_density"),dns_prt);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("distribution_shape"),psd_lst[0].typ_get());
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("mie_solver"),slv_sng);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("mie_openmp_thread_number"),thr_nbr);
  if(!xpt_dsc.empty()) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("experiment_description"),xpt_dsc);
  
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_core"),cmp_sng_cor);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_medium"),cmp_sng_mdm);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_mantle"),cmp_sng_mnt);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_matrix"),cmp_sng_mtx);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_inclusion"),cmp_sng_ncl);
  rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("composition_particle"),cmp_sng_prt);
  if(idx_rfr_cor_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("core_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("core_refractive_indices_source"),fl_idx_rfr_cor);
  if(idx_rfr_mdm_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("medium_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("medium_refractive_indices_source"),fl_idx_rfr_mdm);
  if(idx_rfr_mnt_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("mantle_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("mantle_refractive_indices_source"),fl_idx_rfr_mnt);
  if(idx_rfr_mtx_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("matrix_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("matrix_refractive_indices_source"),fl_idx_rfr_mtx);
  if(idx_rfr_ncl_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("inclusion_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("inclusion_refractive_indices_source"),fl_idx_rfr_ncl);
  if(idx_rfr_prt_usr_flg) rcd=nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("particle_refractive_indices_source"),static_cast<std::string>("User specified")); else nco_put_att(nc_out,NC_GLOBAL,static_cast<std::string>("particle_refractive_indices_source"),fl_idx_rfr_prt);

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    // Write all coordinates first to ensure record coordinate size is known
    var_mtd_sct var_mtd[]={
      {0,"sz",nco_xtyp,1,dmn_sz,"long_name","Size","units","meter"},
      {0,"dsd_sz",nco_xtyp,1,dmn_dsd,"long_name","Raindrop diameter","units","meter"},
      {0,"sz_grd",nco_xtyp,1,dmn_sz_grd,"long_name","Size grid","units","meter"},
      {0,"wvl",nco_xtyp,1,dmn_wvl,"long_name","Wavelength at band center","units","meter"},
      {0,"wvl_grd",nco_xtyp,1,dmn_wvl_grd,"long_name","Wavelength grid","units","meter"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // After writing, delete arrays 
    rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("sz"),sz_nbr,sz_ctr); // Put sz_ctr into sz coordinate
    rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("dsd_sz"),dsd_nbr,dsd_sz); // dnd: class member
    rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("sz_grd"),sz_nbr+1,sz_grd);
    rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("wvl"),wvl_nbr,wvl);
    rcd=nco_put_vara_crd(nc_out,static_cast<std::string>("wvl_grd"),wvl_nbr+1,wvl_grd);
  } // endif true

  // Call modules which implement their own output, then write everything else

  // Compute phase function diagnostics when optical properties known
  if(mie_flg && !coat_flg) 
    rcd+=phz_fnc_mdl // [fnc] Phase function diagnostics module
      (dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file 
       nc_out, // I [fl] netCDF file for output 
       lgn_nbr, // I [nbr] Order of phase function Legendre expansion
       ngl_nbr, // I [nbr] Number of polar angles in one hemisphere
       wvl_idx_dbg, // I [idx] Debugging wavelength bin
       wvl_nbr, // I [nbr] Number of output wavelength bands
       ngl_dbg_dgr, // I [dgr] Debugging angle
       ngl, // I [rdn] Scattering angle
       ngl_dgr, // I [dgr] Angle degrees
       ngl_dlt, // I [rdn] Width of angle bin
       ngl_wgt, // I [rdn] Weight of angle bin
       phz_fnc_ffc, // I [sr-1] Effective phase function (weighted over size distribution and sub-bands)
       phz_fnc_dgn, // I [sr-1] Phase function at diagnostic size, wavelength
       plz_dgn); // I [frc] Degree of linear polarization at diagnostic size, wavelength
  
  // Raindrop chemistry
  // 20060817: PGI pgCC fails in here
  rcd+=rnd_chm // [fnc] Raindrop chemistry
    (nc_out, // I [fl] netCDF file for output 
     dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
     tpt_mdp, // I [K] Midlayer temperature
     sz_nbr, // I [nbr] Number of size bins
     cnc, // I [# m-3] Number concentration 
     mss, // I [kg] Mass 
     rds_ctr, // I [m] Radius at bin center
     vmr_CO2); // [mlc mlc-1] Volume mixing ratio of CO2

  // Fresnel reflectance of light refracted through homogeneous plane surface
  if(tst_sng != "nsz") // Fresnel test is somewhat memory intensive
    rcd+=rfl_frs // [fnc] Fresnel reflectance
      (nc_out, // I [fl] netCDF file for output 
       dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
       wvl_nbr, // I [nbr] Number of wavelength bins
       idx_rfr_mdm, // I [frc] Refractive index of transmitted medium
       idx_rfr_cor, // I [frc] Refractive index of incident medium
       slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle
  
  // Radiative transfer properties of aerosol layer
  prc_cmp *tau_ext=new prc_cmp[wvl_nbr]; // [frc] Extinction optical depth
  prc_cmp *tau_abs=new prc_cmp[wvl_nbr]; // [frc] Absorption optical depth
  prc_cmp *tau_sct=new prc_cmp[wvl_nbr]; // [frc] Scattering optical depth
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    // Full layer thickness is 2.0*hgt_mdp
    tau_ext[wvl_idx]=2.0*hgt_mdp*ext_cff_dst[wvl_idx]; // [frc] Extinction optical depth
    tau_abs[wvl_idx]=2.0*hgt_mdp*abs_cff_dst[wvl_idx]; // [frc] Absorption optical depth
    tau_sct[wvl_idx]=2.0*hgt_mdp*sca_cff_dst[wvl_idx]; // [frc] Scattering optical depth
  } // end loop over wvl

  // Determine 
  if(mie_flg && tst_sng != "nsz"){
    rt_cls rt_obj("two_srm_iso_sct",wvl_nbr,asm_prm,ss_alb,tau_ext); // [rt] Radiative transfer object
    // Layer optical properties
    rcd+=rt_lop // [fnc] Radiative transfer layer optical properties
      (nc_out, // I [fl] netCDF file for output 
       rt_obj, // I [rt] Radiative transfer object
       slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle

    prc_cmp *flx_frc_drc_sfc=new prc_cmp[wvl_nbr]; // [frc] Surface insolation fraction in direct beam
    for(idx=0;idx<wvl_nbr;idx++) flx_frc_drc_sfc[idx]=flx_frc_drc_sfc_cmd_ln; // [frc] Surface insolation fraction in direct beam
    // Snowpack optics
    rcd+=rt_rfl_snw // [fnc] Radiative transfer for snow reflectance
      (nc_out, // I [fl] netCDF file for output 
       flx_frc_drc_sfc, // I [frc] Surface insolation fraction in direct beam
       rfl_gnd_dff, // I [frc] Diffuse reflectance of ground (beneath snow)
       rt_obj, // I [rt] Radiative transfer object
       slr_zen_ngl_cos); // I [frc] Cosine solar zenith angle
  } // endif mie_flg
  
  // Diagnose aerosol heating
  // fxm: Input is consistent with multi-modal distributions except passing psd_lst[0]
  if((dbg_lvl == dbg_old || tst_sng == "htg") && mie_flg){
    rcd+=aer_htg // [fnc] Determine aerosol heating characteristics
      (nc_out, // I [fl] netCDF file for output 
       dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
       dns_mdp, // I [kg m-3] Midlayer density
       prs_mdp, // I [Pa] Midlayer pressure
       tpt_mdp, // I [K] Midlayer temperature
       wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
       dmt_dtc, // I [m] Diameter of detector
       flx_slr_src, // I [obj] Solar spectrum
       cmp_prt, // I [obj] Particle component
       psd_lst, // I [obj] Particle size distribution
       abs_fsh, // I [frc] Absorption efficiency
       cnc, // I [# m-3] Number concentration 
       mss, // I [kg] Mass 
       rds_ctr, // I [m] Radius at bin center
       xsa, // I [m2] Cross-sectional area
       ss_co_alb_fsh, // I [frc] Single scattering co-albedo
       sz_nbr, // I [nbr] Number of size bins
       abs_cff_mss[wvl_idx_dbg], // I [m2 kg-1] Mass absorption coefficient
       mss_rsl, // I [kg m-3] Mass concentration resolved
       ss_co_alb[wvl_idx_dbg]); // I [frc] Single scattering co-albedo
  } // endif aer_htg
  
  if(dbg_lvl == dbg_old || tst_sng == "lbl"){
    std::cout << "Testing line-by-line and HITRAN routines..." << std::endl;
    rcd+=rt_lbl // [fnc] Single line radiative transfer
      (nc_out, // I [fl] netCDF file for output 
       dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
       prs_ntf-prs_mdp, // I [Pa] Pressure thickness
       prs_mdp, // I [Pa] Midlayer pressure
       slr_zen_ngl_cos, // I [frc] Cosine solar zenith angle
       tpt_mdp, // I [K] Midlayer temperature
       vmr_CO2, // I [mlc mlc-1] Volume mixing ratio of CO2
       lbl_tst, // I [sng] Name of line-by-line test
       wvlgrd); // I [obj] Wavelength grid
  } // end if dbg || "lbl"

  // Free un-needed memory from recent calculations
  delete[] idx_rfr_ffc; // [frc] Effective refractive index of particle
  delete[] idx_rfr_ffc_brg; // [frc] Effective refractive index, Bruggeman approximation
  delete[] idx_rfr_ffc_mxg; // [frc] Effective refractive index, Maxwell Garnett approximation
  delete[] idx_rfr_ffc_pmr; // [frc] Effective refractive index, partial molar refraction approximation
  delete[] idx_rfr_ffc_vlw; // [frc] Effective refractive index, volume-weighted approximation
  delete[] idx_rfr_air; // [frc] Refractive index of air
  delete[] idx_rfr_cor; // [frc] Refractive index of core
  delete[] idx_rfr_mdm; // [frc] Refractive index of medium
  delete[] idx_rfr_mnt; // [frc] Refractive index of mantle
  delete[] idx_rfr_mtx; // [frc] Refractive index of matrix
  delete[] idx_rfr_ncl; // [frc] Refractive index of inclusion
  delete[] idx_rfr_prt; // [frc] Refractive index of particle
  delete[] idx_rfr_ffc_wgt; // [frc] Spectral flux-weighted effective refractive index

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    // Variables required and output by surface temperature iteration
    var_mtd_sct var_mtd[]={
      {0,"cnd_trm_soi",NC_FLOAT,0,dmn_scl,"long_name","Soil thermal conductivity","units","watt meter-1 kelvin-1"},
      {0,"flx_LW_dwn_sfc",NC_FLOAT,0,dmn_scl,"long_name","Longwave downwelling flux at surface","units","watt meter-2"},
      {0,"flx_LW_upw_sfc",NC_FLOAT,0,dmn_scl,"long_name","Longwave upwelling flux at surface","units","watt meter-2"},
      {0,"flx_LW_upw_sfc_dps",NC_FLOAT,0,dmn_scl,"long_name","Longwave upwelling flux at surface","units","watt meter-2"},
      {0,"flx_SW_net_gnd",NC_FLOAT,0,dmn_scl,"long_name","Solar flux absorbed by ground","units","watt meter-2"},
      {0,"flx_SW_net_vgt",NC_FLOAT,0,dmn_scl,"long_name","Solar flux absorbed by vegetation","units","watt meter-2"},
      {0,"flx_ltn",NC_FLOAT,0,dmn_scl,"long_name","Latent heat flux to atmosphere","units","watt meter-2"},
      {0,"flx_ltn_dps",NC_FLOAT,0,dmn_scl,"long_name","Latent heat flux to atmosphere","units","watt meter-2"},
      {0,"flx_sns_atm",NC_FLOAT,0,dmn_scl,"long_name","Sensible heat flux to atmosphere","units","watt meter-2"},
      {0,"flx_sns_atm_dps",NC_FLOAT,0,dmn_scl,"long_name","Sensible heat flux to atmosphere","units","watt meter-2"},
      {0,"flx_sns_gnd",NC_FLOAT,0,dmn_scl,"long_name","Sensible heat flux to soil","units","watt meter-2"},
      {0,"flx_snw_mlt",NC_FLOAT,0,dmn_scl,"long_name","Snow melt heat flux","units","watt meter-2"},
      {0,"flx_q_H2O",NC_FLOAT,0,dmn_scl,"long_name","Moisture flux to atmosphere","units","kilogram meter-1 second-2"},
      {0,"flx_q_H2O_dps",NC_FLOAT,0,dmn_scl,"long_name","Moisture flux to atmosphere","units","kilogram meter-1 second-2"},
      {0,"pnt_typ_idx",NC_SHORT,0,dmn_scl,"long_name","Plant type index","units","index"},
      {0,"trn_fsh_vpr_soi_atm",NC_FLOAT,0,dmn_scl,"long_name","Transfer efficiency of vapor from soil to atmosphere","units","fraction"},

      {0,"tpt_aer",NC_FLOAT,0,dmn_scl,"long_name","\"Aerodynamic\" temperature at z=zpd+rgh_mmn","units","kelvin"},
      {0,"tpt_ash",NC_FLOAT,0,dmn_scl,"long_name","\"Surface\" temperature at z=zpd+rgh_heat","units","kelvin"},
      {0,"tpt_ffc",NC_FLOAT,0,dmn_scl,"long_name","Radiative effective temperature","units","kelvin"},
      {0,"tpt_gnd",NC_FLOAT,0,dmn_scl,"long_name","Ground temperature","units","kelvin"},
      {0,"tpt_gnd_mbl",NC_FLOAT,0,dmn_scl,"long_name","Ground temperature (mobilization)","units","kelvin"},
      {0,"tpt_ice",NC_FLOAT,0,dmn_scl,"long_name","Sea ice temperature","units","kelvin"},
      {0,"tpt_sfc",NC_FLOAT,0,dmn_scl,"long_name","Surface temperature","units","kelvin"},
      {0,"tpt_soi",NC_FLOAT,0,dmn_scl,"long_name","Soil temperature","units","kelvin"},
      {0,"tpt_sst",NC_FLOAT,0,dmn_scl,"long_name","Sea surface temperature","units","kelvin"},
      {0,"tpt_vgt",NC_FLOAT,0,dmn_scl,"long_name","Vegetation temperature","units","kelvin"},
      {0,"wnd_mrd_mdp",NC_FLOAT,0,dmn_scl,"long_name","Surface layer meridional wind speed","units","meter second-1"},
      {0,"wnd_mrd_rfr",NC_FLOAT,0,dmn_scl,"long_name","Reference level meridional wind speed","units","meter second-1"},
      {0,"wnd_str_mrd",NC_FLOAT,0,dmn_scl,"long_name","Meridional wind stress","units","kilogram meter-1 second-2"},
      {0,"wnd_str_mrd_dps",NC_FLOAT,0,dmn_scl,"long_name","Meridional wind stress","units","kilogram meter-1 second-2"},
      {0,"wnd_str_znl",NC_FLOAT,0,dmn_scl,"long_name","Zonal wind stress","units","kilogram meter-1 second-2"},
      {0,"wnd_str_znl_dps",NC_FLOAT,0,dmn_scl,"long_name","Zonal wind stress","units","kilogram meter-1 second-2"},
      {0,"wnd_znl_mdp",NC_FLOAT,0,dmn_scl,"long_name","Surface layer zonal wind speed","units","meter second-1"},
      {0,"wnd_znl_rfr",NC_FLOAT,0,dmn_scl,"long_name","Reference level zonal wind speed","units","meter second-1"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnd_trm_soi"),cnd_trm_soi);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_LW_dwn_sfc"),flx_LW_dwn_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_LW_upw_sfc"),flx_LW_upw_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_LW_upw_sfc_dps"),flx_LW_upw_sfc_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_SW_net_gnd"),flx_SW_net_gnd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_SW_net_vgt"),flx_SW_net_vgt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_ltn"),flx_ltn);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_ltn_dps"),flx_ltn_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_sns_atm"),flx_sns_atm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_sns_atm_dps"),flx_sns_atm_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_sns_gnd"),flx_sns_gnd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_snw_mlt"),flx_snw_mlt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_q_H2O"),flx_q_H2O);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_q_H2O_dps"),flx_q_H2O_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("pnt_typ_idx"),pnt_typ_idx);
    rcd=nco_put_var(nc_out,static_cast<std::string>("trn_fsh_vpr_soi_atm"),trn_fsh_vpr_soi_atm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_aer"),tpt_aer);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_ash"),tpt_ash);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_ffc"),tpt_ffc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_gnd"),tpt_gnd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_gnd_mbl"),tpt_gnd_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_ice"),tpt_ice);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_sfc"),tpt_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_soi"),tpt_soi);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_sst"),tpt_sst);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_vgt"),tpt_vgt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_mrd_mdp"),wnd_mrd_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_mrd_rfr"),wnd_mrd_rfr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_str_mrd"),wnd_str_mrd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_str_mrd_dps"),wnd_str_mrd_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_str_znl"),wnd_str_znl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_str_znl_dps"),wnd_str_znl_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_znl_mdp"),wnd_znl_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_znl_rfr"),wnd_znl_rfr);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"dsd_nbr",NC_INT,0,dmn_scl,"long_name","Number of raindrop size bins","units","number"},
      {0,"clc_fsh",NC_FLOAT,2,dmn_dsd_sz,"long_name","Collection efficiency","units","fraction"},
      {0,"cll_fsh",NC_FLOAT,2,dmn_dsd_sz,"long_name","Collision efficiency","units","fraction"},
      {0,"cll_fsh_brn_dff",NC_FLOAT,2,dmn_dsd_sz,"long_name","Collision efficiency of Brownian diffusion","units","fraction"},
      {0,"cll_fsh_mpc",NC_FLOAT,2,dmn_dsd_sz,"long_name","Collision efficiency of impaction","units","fraction"},
      {0,"cll_fsh_ntc",NC_FLOAT,2,dmn_dsd_sz,"long_name","Collision efficiency of interception","units","fraction"},
      {0,"cnc_nbr_pcp_anl",NC_FLOAT,0,dmn_scl,"long_name","Number concentration analytic, raindrops","units","meter-3"},
      {0,"cnc_nbr_pcp_rsl",NC_FLOAT,0,dmn_scl,"long_name","Number concentration resolved, raindrops","units","meter-3"},
      {0,"cnc_nbr_rsl_frc_pcp",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of number concentration, raindrops","units","fraction"},
      {0,"cnc_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop number concentration","units","meter-3"},
      {0,"dmt_dlt_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Width of diameter bin, raindrop","units","meter"},
      {0,"dmt_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop diameter","units","meter"},
      {0,"dmt_pcp_nma",NC_FLOAT,0,dmn_scl,"long_name","Raindrop number median diameter analytic","units","meter"},
      {0,"dmt_pcp_vwr_Mas71",NC_FLOAT,0,dmn_scl,"long_name","Raindrop diameter volume weighted resolved, Mas71 parameterization","units","meter"},
      {0,"dst_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Number distribution, raindrops","units","meter-3 meter-1"},
      {0,"flx_cnc_pcp_rsl",NC_FLOAT,0,dmn_scl,"long_name","Raindrop number flux, resolved","units","meter-2 second-1"},
      {0,"flx_cnc_spc_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop spectral number flux","units","meter-2 second-1 meter-1"},
      {0,"flx_mss_pcp_rsl",NC_FLOAT,0,dmn_scl,"long_name","Raindrop mass flux, resolved","units","kilogram meter-2 second-1"},
      {0,"flx_mss_spc_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop spectral mass flux","units","kilogram meter-2 second-1 meter-1"},
      {0,"flx_vlm_pcp_rsl",NC_FLOAT,0,dmn_scl,"long_name","Raindrop volume flux, resolved","units","meter3 meter-2 second-1"},
      {0,"flx_vlm_spc_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop spectral volume flux, resolved","units","meter second-1 meter-1"},
      {0,"gsd_pcp_anl",NC_FLOAT,0,dmn_scl,"long_name","Geometric standard deviation, raindrops, analytic","units","fraction"},
      {0,"mss_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop mass","units","kilogram"},
      {0,"rds_pcp",nco_xtyp,1,dmn_dsd,"long_name","Raindrop radius","units","meter"},
      {0,"ryn_nbr_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Reynolds number of raindrop","units","fraction"},
      {0,"scv_cff",NC_FLOAT,1,dmn_sz,"long_name","Scavenging coefficient","units","second-1"},
      {0,"scv_cff_mss_avg",NC_FLOAT,0,dmn_scl,"long_name","Mass mean scavenging coefficient","units","second-1"},
      {0,"scv_cff_mss_avg_pcp_nrm",NC_FLOAT,0,dmn_scl,"long_name","Mass mean scavenging coefficient","units","meter2 kilogram-1"},
      {0,"scv_cff_nbr_avg",NC_FLOAT,0,dmn_scl,"long_name","Number mean scavenging coefficient","units","second-1"},
      {0,"scv_cff_nbr_avg_pcp_nrm",NC_FLOAT,0,dmn_scl,"long_name","Number mean scavenging coefficient","units","meter2 kilogram-1"},
      {0,"scv_cff_pcp_nrm",NC_FLOAT,1,dmn_sz,"long_name","Scavenging coefficient, precipitation normalized","units","meter2 kilogram-1"},
      {0,"sfc_rsl_frc_pcp",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of surface area concentration, raindrops","units","fraction"},
      {0,"stc_fsh",NC_FLOAT,1,dmn_sz,"long_name","Sticking (coalescence) efficiency","units","fraction"},
      {0,"stk_nbr_crt",NC_FLOAT,1,dmn_dsd,"long_name","Critical Stokes number","units","fraction"},
      {0,"stk_nbr_rlt",NC_FLOAT,2,dmn_dsd_sz,"long_name","Stokes number of relative flow","units","fraction"},
      {0,"tau_rlx",NC_FLOAT,1,dmn_sz,"long_name","Relaxation timescale","units","second"},
      {0,"vlc_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop fall speed","units","meter second-1"},
      {0,"vlc_pcp_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted gravitational settling velocity, raindrops","units","meter second-1"},
      {0,"vlc_pcp_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted gravitational settling velocity, raindrops","units","meter second-1"},
      {0,"vlc_stk_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Stokes settling velocity of raindrop","units","meter second-1"},
      {0,"vlm_pcp",NC_FLOAT,1,dmn_dsd,"long_name","Raindrop volume","units","meter3"},
      {0,"vlm_rsl_frc_pcp",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of volume concentration, raindrops","units","fraction"},
      {0,"xsa_rsl_frc_pcp",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of cross-sectional area concentration, raindrops","units","fraction"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // Handle multi-dimensional arrays separately to conserve space if requested
    if(dmn_nbr_max >= 2){
      // Do not delete static multidimensional arrays
      rcd=nco_put_var(nc_out,static_cast<std::string>("clc_fsh"),&clc_fsh(0,0));
      rcd=nco_put_var(nc_out,static_cast<std::string>("cll_fsh"),&cll_fsh(0,0));
      rcd=nco_put_var(nc_out,static_cast<std::string>("cll_fsh_brn_dff"),&cll_fsh_brn_dff(0,0));
      rcd=nco_put_var(nc_out,static_cast<std::string>("cll_fsh_mpc"),&cll_fsh_mpc(0,0));
      rcd=nco_put_var(nc_out,static_cast<std::string>("cll_fsh_ntc"),&cll_fsh_ntc(0,0));
      rcd=nco_put_var(nc_out,static_cast<std::string>("stk_nbr_rlt"),&stk_nbr_rlt(0,0));
    } // endif dmn_nbr_max >= 2

    // After writing, delete dynamic arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_pcp_anl"),cnc_nbr_pcp_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_pcp_rsl"),cnc_nbr_pcp_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_rsl_frc_pcp"),cnc_nbr_rsl_frc_pcp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_pcp"),cnc_pcp); delete []cnc_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_dlt_pcp"),dmt_dlt_pcp); delete []dmt_dlt_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_pcp"),dmt_pcp); // dnd: class member
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_pcp_nma"),dmt_pcp_nma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_pcp_vwr_Mas71"),dmt_pcp_vwr_Mas71);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_pcp"),dst_pcp); delete []dst_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_cnc_pcp_rsl"),flx_cnc_pcp_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_cnc_spc_pcp"),flx_cnc_spc_pcp); delete []flx_cnc_spc_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_mss_pcp_rsl"),flx_mss_pcp_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_mss_spc_pcp"),flx_mss_spc_pcp); delete []flx_mss_spc_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_vlm_pcp_rsl"),flx_vlm_pcp_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_vlm_spc_pcp"),flx_vlm_spc_pcp); delete []flx_vlm_spc_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("gsd_pcp_anl"),gsd_pcp_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_pcp"),mss_pcp); delete []mss_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_pcp"),rds_pcp); delete []rds_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr_pcp"),ryn_nbr_pcp); delete []ryn_nbr_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff"),scv_cff); delete []scv_cff;
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff_mss_avg"),scv_cff_mss_avg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff_mss_avg_pcp_nrm"),scv_cff_mss_avg_pcp_nrm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff_nbr_avg"),scv_cff_nbr_avg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff_nbr_avg_pcp_nrm"),scv_cff_nbr_avg_pcp_nrm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("scv_cff_pcp_nrm"),scv_cff_pcp_nrm); delete []scv_cff_pcp_nrm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_rsl_frc_pcp"),sfc_rsl_frc_pcp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("stc_fsh"),stc_fsh); delete []stc_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("stk_nbr_crt"),stk_nbr_crt); delete []stk_nbr_crt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_rlx"),tau_rlx); delete []tau_rlx;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_pcp"),vlc_pcp); delete []vlc_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_pcp_nwr"),vlc_pcp_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_pcp_vwr"),vlc_pcp_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_stk_pcp"),vlc_stk_pcp); delete []vlc_stk_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_pcp"),vlm_pcp); delete []vlm_pcp;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_rsl_frc_pcp"),vlm_rsl_frc_pcp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_rsl_frc_pcp"),xsa_rsl_frc_pcp);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"psd_dmt_nma",NC_FLOAT,1,dmn_psd,"long_name","Number median diameter analytic","units","meter"},
      {0,"psd_rds_nma",NC_FLOAT,1,dmn_psd,"long_name","Number median radius analytic","units","meter"},
      {0,"psd_gsd_anl",NC_FLOAT,1,dmn_psd,"long_name","Geometric standard deviation","units","fraction"},
      {0,"psd_cnc_nbr_anl",NC_FLOAT,1,dmn_psd,"long_name","Number concentration analytic","units","meter-3"},
      {0,"psd_mss_frc_anl",NC_FLOAT,1,dmn_psd,"long_name","Mass fraction analytic","units","fraction"},
      {0,"psd_nbr",NC_FLOAT,0,dmn_scl,"long_name","Number of size distribution modes","units","number"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_dmt_nma"),psd_dmt_nma); delete []psd_dmt_nma;
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_rds_nma"),psd_rds_nma); delete []psd_rds_nma;
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_gsd_anl"),psd_gsd_anl); delete []psd_gsd_anl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_cnc_nbr_anl"),psd_cnc_nbr_anl); delete []psd_cnc_nbr_anl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_mss_frc_anl"),psd_mss_frc_anl); delete []psd_mss_frc_anl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("psd_nbr"),psd_nbr);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"cff_drg_grv",NC_FLOAT,1,dmn_sz,"long_name","Drag coefficient at terminal velocity","units","fraction"},
      {0,"dff_aer",NC_FLOAT,1,dmn_sz,"long_name","Brownian diffusivity of particle in atmosphere","units","meter2 second-1"},
      {0,"dff_aer_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted Brownian diffusivity of particle in atmosphere","units","meter2 second-1"},
      {0,"dff_aer_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted Brownian diffusivity of particle in atmosphere","units","meter2 second-1"},
      {0,"dmt_slt_opt",NC_FLOAT,0,dmn_scl,"long_name","Optimal diameter for saltation","units","meter"},
      {0,"mfp_atm",NC_FLOAT,0,dmn_scl,"long_name","Mean free path of air molecules","units","meter"},
      {0,"rss_lmn",NC_FLOAT,1,dmn_sz,"long_name","Quasi-laminar layer resistance","units","second meter-1"},
      {0,"rss_trb",NC_FLOAT,1,dmn_sz,"long_name","Resistance to turbulent deposition","units","second meter-1"},
      {0,"vlc_trb",NC_FLOAT,1,dmn_sz,"long_name","Turbulent deposition velocity","units","meter second-1"},
      {0,"vlc_dry",NC_FLOAT,1,dmn_sz,"long_name","Total dry deposition velocity","units","meter second-1"},
      {0,"vlc_trb_SeH78",NC_FLOAT,1,dmn_sz,"long_name","Turbulent deposition velocity","units","meter second-1"},
      {0,"vlc_dry_SeH78",NC_FLOAT,1,dmn_sz,"long_name","Total dry deposition velocity","units","meter second-1"},
      {0,"rss_lmn_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted quasi-laminar layer resistance","units","second meter-1"},
      {0,"rss_lmn_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted quasi-laminar layer resistance","units","second meter-1"},
      {0,"ryn_nbr",NC_FLOAT,1,dmn_sz,"long_name","Reynolds number","units","fraction"},
      {0,"ryn_nbr_frc",NC_FLOAT,1,dmn_sz,"long_name","Friction Reynolds number","units","fraction"},
      {0,"ryn_nbr_frc_thr",NC_FLOAT,1,dmn_sz,"long_name","Threshold friction Reynolds number","units","fraction"},
      {0,"ryn_nbr_frc_thr_prx",NC_FLOAT,1,dmn_sz,"long_name","Threshold friction Reynolds number approximation","units","fraction"},
      {0,"ryn_nbr_grv",NC_FLOAT,1,dmn_sz,"long_name","Reynolds number at gravitational settling speed","units","fraction"},
      {0,"pcl_nbr",NC_FLOAT,1,dmn_sz,"long_name","Peclet number","units","fraction"},
      {0,"shm_nbr",NC_FLOAT,1,dmn_sz,"long_name","Schmidt number","units","fraction"},
      {0,"shm_nbr_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted Schmidt number","units","fraction"},
      {0,"shm_nbr_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted Schmidt number","units","fraction"},
      {0,"slp_crc",NC_FLOAT,1,dmn_sz,"long_name","Slip correction factor","units","fraction"},
      {0,"stk_crc",NC_FLOAT,1,dmn_sz,"long_name","Correction to Stokes settling velocity","units","fraction"},
      {0,"stk_nbr",NC_FLOAT,1,dmn_sz,"long_name","Stokes number","units","fraction"},
      {0,"stk_nbr_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted Stokes number","units","fraction"},
      {0,"stk_nbr_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted Stokes number","units","fraction"},
      {0,"vlc_grv",NC_FLOAT,1,dmn_sz,"long_name","Gravitational settling velocity","units","meter second-1"},
      {0,"vlc_grv_nwr",NC_FLOAT,0,dmn_scl,"long_name","Number weighted gravitational settling velocity","units","meter second-1"},
      {0,"vlc_grv_vwr",NC_FLOAT,0,dmn_scl,"long_name","Mass weighted gravitational settling velocity","units","meter second-1"},
      {0,"vlc_stk",NC_FLOAT,1,dmn_sz,"long_name","Stokes' settling velocity","units","meter second-1"},
      {0,"vsc_dyn_atm",NC_FLOAT,0,dmn_scl,"long_name","Dynamic viscosity of atmosphere","units","kilogram meter-1 second-1"},
      {0,"vsc_knm_atm",NC_FLOAT,0,dmn_scl,"long_name","Dynamic viscosity of atmosphere","units","meter2 second-1"},
      {0,"wnd_frc_fsh_frc",NC_FLOAT,0,dmn_scl,"long_name","Efficient fraction of wind friction","units","meter second-1"},
      {0,"wnd_frc_slt",NC_FLOAT,0,dmn_scl,"long_name","Saltation friction speed","units","meter second-1"},
      {0,"wnd_frc_slt_dlt",NC_FLOAT,0,dmn_scl,"long_name","Increase in friction speed from saltation","units","meter second-1"},
      {0,"wnd_frc_smt",NC_FLOAT,0,dmn_scl,"long_name","Friction speed over erodible surface","units","meter second-1"},
      {0,"wnd_frc_thr",NC_FLOAT,1,dmn_sz,"long_name","Threshold friction speed","units","meter second-1"},
      {0,"wnd_frc_thr_opt",NC_FLOAT,0,dmn_scl,"long_name","Optimal threshold friction speed for saltation","units","meter second-1"},
      {0,"wnd_frc_thr_prx",NC_FLOAT,1,dmn_sz,"long_name","Threshold friction speed approximation","units","meter second-1"},
      {0,"wnd_frc_thr_slt",NC_FLOAT,0,dmn_scl,"long_name","Threshold friction speed for saltation","units","meter second-1"},
      {0,"wnd_frc_thr_slt_dry_flt",NC_FLOAT,0,dmn_scl,"long_name","Threshold friction speed for saltation on dry bare ground","units","meter second-1"},
      {0,"wnd_rfr_thr",NC_FLOAT,1,dmn_sz,"long_name","Threshold 10 m wind speed","units","meter second-1"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("cff_drg_grv"),cff_drg_grv); delete []cff_drg_grv;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_aer"),dff_aer); delete []dff_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_aer_nwr"),dff_aer_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_aer_vwr"),dff_aer_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_slt_opt"),dmt_slt_opt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mfp_atm"),mfp_atm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_lmn"),rss_lmn); delete []rss_lmn;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_trb"),rss_trb); delete []rss_trb;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_trb"),vlc_trb); delete []vlc_trb;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_dry"),vlc_dry); delete []vlc_dry;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_trb_SeH78"),vlc_trb_SeH78); delete []vlc_trb_SeH78;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_dry_SeH78"),vlc_dry_SeH78); delete []vlc_dry_SeH78;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_lmn_nwr"),rss_lmn_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_lmn_vwr"),rss_lmn_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr"),ryn_nbr); delete []ryn_nbr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr_frc"),ryn_nbr_frc); delete []ryn_nbr_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr_frc_thr"),ryn_nbr_frc_thr); delete []ryn_nbr_frc_thr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr_frc_thr_prx"),ryn_nbr_frc_thr_prx); delete []ryn_nbr_frc_thr_prx;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ryn_nbr_grv"),ryn_nbr_grv); delete []ryn_nbr_grv;
    rcd=nco_put_var(nc_out,static_cast<std::string>("pcl_nbr"),pcl_nbr); delete []pcl_nbr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("shm_nbr"),shm_nbr); delete []shm_nbr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("shm_nbr_nwr"),shm_nbr_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("shm_nbr_vwr"),shm_nbr_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("slp_crc"),slp_crc); delete []slp_crc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("stk_crc"),stk_crc); delete []stk_crc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("stk_nbr"),stk_nbr); delete []stk_nbr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("stk_nbr_nwr"),stk_nbr_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("stk_nbr_vwr"),stk_nbr_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_grv"),vlc_grv); delete []vlc_grv;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_grv_nwr"),vlc_grv_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_grv_vwr"),vlc_grv_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_stk"),vlc_stk); delete []vlc_stk;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vsc_dyn_atm"),vsc_dyn_atm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vsc_knm_atm"),vsc_knm_atm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_fsh_frc"),wnd_frc_fsh_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_slt"),wnd_frc_slt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_slt_dlt"),wnd_frc_slt_dlt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_smt"),wnd_frc_smt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_thr"),wnd_frc_thr); delete []wnd_frc_thr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_thr_opt"),wnd_frc_thr_opt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_thr_prx"),wnd_frc_thr_prx); delete []wnd_frc_thr_prx;
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_thr_slt"),wnd_frc_thr_slt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_thr_slt_dry_flt"),wnd_frc_thr_slt_dry_flt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_rfr_thr"),wnd_rfr_thr); delete []wnd_rfr_thr;
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"RH_ice",NC_FLOAT,0,dmn_scl,"long_name","Relative humidity w/r/t ice water","units","fraction"},
      {0,"RH_lqd",NC_FLOAT,0,dmn_scl,"long_name","Relative humidity w/r/t liquid water","units","fraction"},
      {0,"dns_mdp",NC_FLOAT,0,dmn_scl,"long_name","Midlayer density","units","kilogram meter-3"},
      {0,"hgt_ash",NC_FLOAT,0,dmn_scl,"long_name","Apparent sink for heat","units","meter"},
      {0,"hgt_asm",NC_FLOAT,0,dmn_scl,"long_name","Apparent sink for momentum","units","meter"},
      {0,"hgt_asv",NC_FLOAT,0,dmn_scl,"long_name","Apparent sink for vapor","units","meter"},
      {0,"hgt_mdp",NC_FLOAT,0,dmn_scl,"long_name","Midlayer height above surface","units","meter"},
      {0,"hgt_rfr",NC_FLOAT,0,dmn_scl,"long_name","Reference height (i.e., 10 m) at which surface winds are evaluated for dust mobilization","units","meter"},
      {0,"hgt_zpd_dps",NC_FLOAT,0,dmn_scl,"long_name","Zero plane displacement height","units","meter"},
      {0,"hgt_zpd_mbl",NC_FLOAT,0,dmn_scl,"long_name","Zero plane displacement height","units","meter"},
      {0,"lnd_frc_dry",NC_FLOAT,0,dmn_scl,"long_name","Dry land fraction","units","fraction"},
      {0,"lnd_frc_mbl",NC_FLOAT,0,dmn_scl,"long_name","Bare ground fraction","units","fraction"},
      {0,"mno_lng_dps",NC_FLOAT,0,dmn_scl,"long_name","Monin-Obukhov length","units","meter"},
      {0,"mno_lng_mbl",NC_FLOAT,0,dmn_scl,"long_name","Monin-Obukhov length","units","meter"},
      {0,"msv_gnd",NC_FLOAT,0,dmn_scl,"long_name","Bare ground emissivity","units","fraction"},
      {0,"msv_sfc",NC_FLOAT,0,dmn_scl,"long_name","Surface (bare ground+snow) emissivity","units","fraction"},
      {0,"msv_snw",NC_FLOAT,0,dmn_scl,"long_name","Snow emissivity","units","fraction"},
      {0,"prs_mdp",NC_FLOAT,0,dmn_scl,"long_name","Midlayer pressure","units","pascal"},
      {0,"prs_ntf",NC_FLOAT,0,dmn_scl,"long_name","Surface pressure","units","pascal"},
      {0,"q_H2O_vpr",NC_FLOAT,0,dmn_scl,"long_name","Specific humidity","units","kilogram kilogram-1"},
      {0,"qst_H2O_ice",NC_FLOAT,0,dmn_scl,"long_name","Saturated specific humidity of H2O w/r/t ice","units","kilogram kilogram-1"},
      {0,"qst_H2O_lqd",NC_FLOAT,0,dmn_scl,"long_name","Saturated specific humidity of H2O w/r/t liquid","units","kilogram kilogram-1"},
      {0,"rgh_heat",NC_FLOAT,0,dmn_scl,"long_name","Roughness length heat","units","meter"},
      {0,"rgh_mmn_dps",NC_FLOAT,0,dmn_scl,"long_name","Roughness length momentum","units","meter"},
      {0,"rgh_mmn_mbl",NC_FLOAT,0,dmn_scl,"long_name","Roughness length momentum","units","meter"},
      {0,"rgh_mmn_smt",NC_FLOAT,0,dmn_scl,"long_name","Roughness length momentum of smooth erodible surfaces","units","meter"},
      {0,"rgh_vpr",NC_FLOAT,0,dmn_scl,"long_name","Roughness length vapor","units","meter"},
      {0,"svp_H2O_ice",NC_FLOAT,0,dmn_scl,"long_name","Saturation vapor pressure of H2O w/r/t ice","units","pascal"},
      {0,"svp_H2O_lqd",NC_FLOAT,0,dmn_scl,"long_name","Saturation vapor pressure of H2O w/r/t liquid","units","pascal"},
      {0,"tpt_mdp",NC_FLOAT,0,dmn_scl,"long_name","Midlayer temperature","units","kelvin"},
      {0,"tpt_ptn_mdp",NC_FLOAT,0,dmn_scl,"long_name","Potential temperature","units","kelvin"},
      {0,"tpt_ptn_vrt_mdp",NC_FLOAT,0,dmn_scl,"long_name","Virtual potential temperature","units","kelvin"},
      {0,"tpt_vrt",NC_FLOAT,0,dmn_scl,"long_name","Virtual temperature","units","kelvin"},
      {0,"tpt_vrt_mdp",NC_FLOAT,0,dmn_scl,"long_name","Virtual temperature","units","kelvin"},
      {0,"wnd_frc_dps",NC_FLOAT,0,dmn_scl,"long_name","Surface friction speed","units","meter second-1"},
      {0,"wnd_frc_mbl",NC_FLOAT,0,dmn_scl,"long_name","Surface friction speed","units","meter second-1"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("RH_ice"),RH_ice);
    rcd=nco_put_var(nc_out,static_cast<std::string>("RH_lqd"),RH_lqd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_mdp"),dns_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_ash"),hgt_ash);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_asm"),hgt_asm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_asv"),hgt_asv);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_mdp"),hgt_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_rfr"),hgt_rfr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_zpd_dps"),hgt_zpd_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("hgt_zpd_mbl"),hgt_zpd_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lnd_frc_dry"),lnd_frc_dry);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lnd_frc_mbl"),lnd_frc_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mno_lng_dps"),mno_lng_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mno_lng_mbl"),mno_lng_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("msv_gnd"),msv_gnd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("msv_sfc"),msv_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("msv_snw"),msv_snw);
    rcd=nco_put_var(nc_out,static_cast<std::string>("prs_mdp"),prs_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("prs_ntf"),prs_ntf);
    rcd=nco_put_var(nc_out,static_cast<std::string>("q_H2O_vpr"),q_H2O_vpr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("qst_H2O_ice"),qst_H2O_ice);
    rcd=nco_put_var(nc_out,static_cast<std::string>("qst_H2O_lqd"),qst_H2O_lqd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rgh_heat"),rgh_heat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rgh_mmn_dps"),rgh_mmn_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rgh_mmn_mbl"),rgh_mmn_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rgh_mmn_smt"),rgh_mmn_smt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rgh_vpr"),rgh_vpr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("svp_H2O_ice"),svp_H2O_ice);
    rcd=nco_put_var(nc_out,static_cast<std::string>("svp_H2O_lqd"),svp_H2O_lqd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_mdp"),tpt_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_ptn_mdp"),tpt_ptn_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_ptn_vrt_mdp"),tpt_ptn_vrt_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_vrt"),tpt_vrt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_vrt_mdp"),tpt_vrt_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_dps"),wnd_frc_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_frc_mbl"),wnd_frc_mbl);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"act_cff",NC_FLOAT,1,dmn_sz,"long_name","Activity coefficient of solution","units","fraction"},
      {0,"mss_frc_solute",NC_FLOAT,1,dmn_sz,"long_name","Mass fraction of solute","units","fraction"},
      {0,"mlmic",NC_FLOAT,0,dmn_scl,"long_name","Mean linear mass increase coefficient","units","fraction"},
      {0,"mss_wtr_rat",NC_FLOAT,1,dmn_sz,"long_name","Ratio of water mass to dry aerosol mass","units","fraction"},
      {0,"svp_fct_klv",NC_FLOAT,1,dmn_sz,"long_name","Saturation vapor pressure curvature (Kelvin) factor","units","fraction"},
      {0,"svp_fct_slt",NC_FLOAT,1,dmn_sz,"long_name","Saturation vapor pressure solute factor","units","fraction"},
      {0,"svp_fct_ttl",NC_FLOAT,1,dmn_sz,"long_name","Saturation vapor pressure correction factor","units","fraction"}
      // {0,"foo",NC_FLOAT,1,dmn_sz,"long_name","","units","fraction"},
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("act_cff"),act_cff); delete []act_cff;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_solute"),mss_frc_solute); delete []mss_frc_solute;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mlmic"),mlmic);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_wtr_rat"),mss_wtr_rat); delete []mss_wtr_rat;
    rcd=nco_put_var(nc_out,static_cast<std::string>("svp_fct_klv"),svp_fct_klv); delete []svp_fct_klv;
    rcd=nco_put_var(nc_out,static_cast<std::string>("svp_fct_slt"),svp_fct_slt); delete []svp_fct_slt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("svp_fct_ttl"),svp_fct_ttl); delete []svp_fct_ttl;
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"cff_xch_mmn_ntr",NC_FLOAT,0,dmn_scl,"long_name","Neutral momentum exchange coefficient, zmdp to zpd+z0m","units","fraction"},
      {0,"cff_xch_mmn",NC_FLOAT,0,dmn_scl,"long_name","Momentum exchange coefficient, zmdp to zpd+z0m","units","fraction"},
      {0,"cff_xch_heat",NC_FLOAT,0,dmn_scl,"long_name","Heat exchange coefficient, zmdp to zpd+z0h","units","fraction"},
      {0,"cff_xch_vpr",NC_FLOAT,0,dmn_scl,"long_name","Vapor exchange coefficient, zmdp to zpd+z0w","units","fraction"},
      {0,"flx_mss_hrz_slt_ttl",NC_FLOAT,0,dmn_scl,"long_name","Vertically integrated streamwise mass flux","units","kilogram meter-1 second-1"},
      {0,"flx_mss_vrt_dst_ttl",NC_FLOAT,0,dmn_scl,"long_name","Total vertical mass flux of dust","units","kilogram meter-2 second-1"},
      {0,"dst_slt_flx_rat_ttl",NC_FLOAT,0,dmn_scl,"long_name","Ratio of vertical dust flux to streamwise mass flux","units","meter-1"},
      {0,"oro",NC_FLOAT,0,dmn_scl,"long_name","Orography: ocean=0.0, land=1.0, sea ice=2.0","units","fraction"},
      {0,"doy",NC_FLOAT,0,dmn_scl,"long_name","Day of year [1.0..367.0)","units","day"},
      {0,"lat_dgr",NC_FLOAT,0,dmn_scl,"long_name","Latitude","units","degrees north"},
      {0,"lat",NC_FLOAT,0,dmn_scl,"long_name","Latitude","units","degrees north"},
      {0,"lat_rdn",NC_FLOAT,0,dmn_scl,"long_name","Latitude","units","radians north"},
      {0,"pugtut",NC_FLOAT,1,dmn_sz,"long_name","Probability wind speed greater than threshold","units","fraction"},
      {0,"rss_aer_mmn_ntr",NC_FLOAT,0,dmn_scl,"long_name","Neutral aerodynamic resistance to momentum transfer, zmdp to zpd+z0m","units","second meter-1"},
      {0,"rss_aer_mmn",NC_FLOAT,0,dmn_scl,"long_name","Aerodynamic resistance to momentum transfer, zmdp to zpd+z0m","units","second meter-1"},
      {0,"rss_aer_heat",NC_FLOAT,0,dmn_scl,"long_name","Aerodynamic resistance to heat transfer, zmdp to zpd+z0h","units","second meter-1"},
      {0,"rss_aer_vpr",NC_FLOAT,0,dmn_scl,"long_name","Aerodynamic resistance to vapor transfer, zmdp to zpd+z0w","units","second meter-1"},
      {0,"sfc_typ",NC_SHORT,0,dmn_scl,"long_name","LSM surface type (0..28)","units","index"},
      {0,"snw_frc",NC_FLOAT,0,dmn_scl,"long_name","Fraction of surface covered by snow","units","fraction"},
      {0,"snw_hgt",NC_FLOAT,0,dmn_scl,"long_name","Geometric bulk thickness of snow","units","meter"},
      {0,"snw_hgt_lqd",NC_FLOAT,0,dmn_scl,"long_name","Equivalent liquid water snow depth","units","meter"},
      {0,"str_shr",NC_FLOAT,0,dmn_scl,"long_name","Shear stress","units","kilogram meter-1 second-2"},
      {0,"wbl_scl",NC_FLOAT,0,dmn_scl,"long_name","Weibull distribution scale factor","units","meter second-1"},
      {0,"wbl_shp",NC_FLOAT,0,dmn_scl,"long_name","Weibull distribution shape factor","units","fraction"},
      {0,"wnd_rfr_mbl",NC_FLOAT,0,dmn_scl,"long_name","Wind speed at reference height","units","meter second-1"},
      {0,"wnd_rfr_dps",NC_FLOAT,0,dmn_scl,"long_name","Wind speed at reference height","units","meter second-1"},
      {0,"wnd_mdp",NC_FLOAT,0,dmn_scl,"long_name","Surface layer mean wind speed","units","meter second-1"},
      {0,"wnd_mdp_thr_slt",NC_FLOAT,0,dmn_scl,"long_name","Threshold midlayer wind speed for saltation","units","meter second-1"},
      {0,"wnd_rfr_thr_slt",NC_FLOAT,0,dmn_scl,"long_name","Threshold 10 m wind speed for saltation","units","meter second-1"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // Add Surface type diagnostics attribute
    rcd=nco_redef(nc_out,NC_EINDEFINE); // [fnc] Put open netCDF dataset into define mode
    rcd=nco_put_att(nc_out,static_cast<std::string>("sfc_typ"),static_cast<std::string>("sfc_typ_sng"),sfc_typ_sng_get(sfc_typ));
    rcd=nco_enddef(nc_out); // [fnc] Leave define mode
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("cff_xch_mmn_ntr"),cff_xch_mmn_ntr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cff_xch_mmn"),cff_xch_mmn);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cff_xch_heat"),cff_xch_heat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cff_xch_vpr"),cff_xch_vpr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_mss_hrz_slt_ttl"),flx_mss_hrz_slt_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_mss_vrt_dst_ttl"),flx_mss_vrt_dst_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_slt_flx_rat_ttl"),dst_slt_flx_rat_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("oro"),oro);
    rcd=nco_put_var(nc_out,static_cast<std::string>("doy"),doy);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lat_dgr"),lat_dgr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lat"),lat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lat_rdn"),lat_rdn);
    rcd=nco_put_var(nc_out,static_cast<std::string>("pugtut"),pugtut); delete []pugtut;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_aer_mmn_ntr"),rss_aer_mmn_ntr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_aer_mmn"),rss_aer_mmn);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_aer_heat"),rss_aer_heat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rss_aer_vpr"),rss_aer_vpr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_typ"),sfc_typ);
    rcd=nco_put_var(nc_out,static_cast<std::string>("snw_frc"),snw_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("snw_hgt"),snw_hgt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("snw_hgt_lqd"),snw_hgt_lqd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("str_shr"),str_shr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wbl_scl"),wbl_scl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wbl_shp"),wbl_shp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_rfr_mbl"),wnd_rfr_mbl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_rfr_dps"),wnd_rfr_dps);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_mdp"),wnd_mdp);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_mdp_thr_slt"),wnd_mdp_thr_slt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wnd_rfr_thr_slt"),wnd_rfr_thr_slt);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"frc_thr_ncr_drg",NC_FLOAT,0,dmn_scl,"long_name","Factor by which surface roughness increases threshold friction velocity","units","fraction"},
      {0,"frc_thr_ncr_wtr",NC_FLOAT,0,dmn_scl,"long_name","Factor by which soil moisture increases threshold friction velocity","units","fraction"},
      {0,"dns_prt_sfc",NC_FLOAT,0,dmn_scl,"long_name","Dry density of soil particles (excluding pores)","units","kilogram meter-3"},
      {0,"dns_blk_dry",NC_FLOAT,0,dmn_scl,"long_name","Bulk density of dry surface soil","units","kilogram meter-3"},
      {0,"dns_blk_sfc",NC_FLOAT,0,dmn_scl,"long_name","Bulk density of surface soil","units","kilogram meter-3"},
      {0,"gwc_sfc",NC_FLOAT,0,dmn_scl,"long_name","Gravimetric water content","units","kilogram kilogram-1"},
      {0,"mss_frc_cly",NC_FLOAT,0,dmn_scl,"long_name","Mass fraction clay","units","fraction"},
      {0,"mss_frc_slt",NC_FLOAT,0,dmn_scl,"long_name","Mass fraction silt","units","fraction"},
      {0,"mss_frc_snd",NC_FLOAT,0,dmn_scl,"long_name","Mass fraction sand","units","fraction"},
      {0,"smp_sat",NC_FLOAT,0,dmn_scl,"long_name","Saturated soil matric potential","units","millimeter H2O"},
      {0,"smp_sfc",NC_FLOAT,0,dmn_scl,"long_name","Soil matric potential","units","millimeter H2O"},
      {0,"smp_xpn_b",NC_FLOAT,0,dmn_scl,"long_name","Exponent b for soil matric potential","units","fraction"},
      {0,"vwc_dry",NC_FLOAT,0,dmn_scl,"long_name","Dry volumetric water content (no E-T)","units","meter3 meter-3"},
      {0,"vwc_opt",NC_FLOAT,0,dmn_scl,"long_name","E-T optimal volumetric water content","units","meter3 meter-3"},
      {0,"vwc_rel",NC_FLOAT,0,dmn_scl,"long_name","Water content relative to saturation","units","fraction"},
      {0,"vwc_sat",NC_FLOAT,0,dmn_scl,"long_name","Saturated volumetric water content","units","meter3 meter-3"},
      {0,"vwc_sfc",NC_FLOAT,0,dmn_scl,"long_name","Volumetric water content","units","meter3 meter-3"},
      {0,"vwc_thr",NC_FLOAT,0,dmn_scl,"long_name"," Threshold volumetric water content to affect mobilization","units","meter3 meter-3"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_blk_sfc"),dns_blk_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_prt_sfc"),dns_prt_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_blk_dry"),dns_blk_dry);
    rcd=nco_put_var(nc_out,static_cast<std::string>("frc_thr_ncr_drg"),frc_thr_ncr_drg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("frc_thr_ncr_wtr"),frc_thr_ncr_wtr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("gwc_sfc"),gwc_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_cly"),mss_frc_cly);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_slt"),mss_frc_slt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_snd"),mss_frc_snd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("smp_sat"),smp_sat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("smp_sfc"),smp_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("smp_xpn_b"),smp_xpn_b);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_dry"),vwc_dry);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_opt"),vwc_opt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_rel"),vwc_rel);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_sat"),vwc_sat);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_sfc"),vwc_sfc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vwc_thr"),vwc_thr);
  } // endif true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"cnc_HNO3_gas",NC_FLOAT,0,dmn_scl,"long_name","Concentration of HNO3","units","molecule meter-3"},
      {0,"cnc_dry_air",NC_FLOAT,0,dmn_scl,"long_name","Concentration of dry air","units","molecule meter-3"},
      {0,"dff_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Normalized diffusion rate of HNO3 to aerosol","units","meter3 second-1 particle-1"},
      {0,"dff_HNO3_aer_cnt",NC_FLOAT,1,dmn_sz,"long_name","Normalized diffusion rate of HNO3 to aerosol","units","meter3 second-1 particle-1"},
      {0,"dff_HNO3_aer_knt",NC_FLOAT,1,dmn_sz,"long_name","Normalized diffusion rate of HNO3 to aerosol","units","meter3 second-1 particle-1"},
      {0,"dff_HNO3_air",NC_FLOAT,0,dmn_scl,"long_name","Binary diffusivity of HNO3 in air","units","meter2 second-1"},
      {0,"knd_nbr_HNO3_air",NC_FLOAT,1,dmn_sz,"long_name","Knudsen number of HNO3 in air","units","fraction"},
      {0,"cnt_rgm_bll_crc_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Correction to continuum regime diffusion approximation due to ballistic effects","units","fraction"},
      {0,"knt_rgm_dff_crc_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Correction to kinetic regime ballistic approximation due to diffusion effects","units","fraction"},
      {0,"mfp_HNO3_air",NC_FLOAT,0,dmn_scl,"long_name","Mean free path of HNO3 in air","units","meter"},
      {0,"mss_acm_cff_HNO3_aer",NC_FLOAT,0,dmn_scl,"long_name","Mean free path of HNO3 in air","units","meter"},
      {0,"mss_upt_cff_HNO3_aer",NC_FLOAT,0,dmn_scl,"long_name","Mean free path of HNO3 in air","units","meter"},
      {0,"q_HNO3_gas",NC_FLOAT,0,dmn_scl,"long_name","Mixing ratio of HNO3","units","kilogram kilogram-1"},
      {0,"rxr_HNO3_gas_aer_vmr",NC_FLOAT,1,dmn_sz,"long_name","Mean rate of HNO3 removal by aerosol","units","molecule molecule-1 second-1"},
      {0,"rxr_HNO3_gas_aer_vmr_ttl",NC_FLOAT,0,dmn_scl,"long_name","Total mean rate of HNO3 removal by aerosol","units","molecule molecule-1 second-1"},
      {0,"rxrc_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Pseudo first order rate coefficient for HNO3 removal by aerosol","units","second-1"},
      {0,"rxrc_HNO3_aer_ttl",NC_FLOAT,0,dmn_scl,"long_name","Total pseudo first order rate coefficient for HNO3 removal by aerosol","units","second-1"},
      {0,"tau_rxrc_HNO3_aer",nco_xtyp,1,dmn_sz,"long_name","e-folding time of HNO3 due to irreversible uptake by current size bin","units","second"},
      {0,"tau_gpd_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Timescale for HNO3 gas-phase diffusion to aerosol","units","second"},
      {0,"tau_ntf_eqm_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Timescale for HNO3 to achieve phase equilibrium at interface","units","second"},
      {0,"tau_aqs_dss_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Timescale for aqueous dissociation of HNO3 to NO3- + H+","units","second"},
      {0,"tau_aqs_dff_HNO3_aer",NC_FLOAT,1,dmn_sz,"long_name","Timescale for aqueous phase diffusion of HNO3 in H2O droplet","units","second"},
      {0,"vlc_mwb_HNO3",NC_FLOAT,0,dmn_scl,"long_name","Thermal speed of HNO3","units","meter second-1"},
      {0,"vmr_HNO3_gas",NC_FLOAT,0,dmn_scl,"long_name","Volume mixing ratio of gaseous HNO3","units","molecule molecule-1"},
      {0,"vnt_crc_aer",NC_FLOAT,1,dmn_sz,"long_name","Ventilation correction factor for gas phase diffusion to aerosol surface","units","fraction"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_HNO3_gas"),cnc_HNO3_gas);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_dry_air"),cnc_dry_air);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_HNO3_aer"),dff_HNO3_aer); delete []dff_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_HNO3_aer_cnt"),dff_HNO3_aer_cnt); delete []dff_HNO3_aer_cnt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_HNO3_aer_knt"),dff_HNO3_aer_knt); delete []dff_HNO3_aer_knt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dff_HNO3_air"),dff_HNO3_air);
    rcd=nco_put_var(nc_out,static_cast<std::string>("knd_nbr_HNO3_air"),knd_nbr_HNO3_air); delete []knd_nbr_HNO3_air;
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnt_rgm_bll_crc_HNO3_aer"),cnt_rgm_bll_crc_HNO3_aer); delete []cnt_rgm_bll_crc_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("knt_rgm_dff_crc_HNO3_aer"),knt_rgm_dff_crc_HNO3_aer); delete []knt_rgm_dff_crc_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mfp_HNO3_air"),mfp_HNO3_air);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_acm_cff_HNO3_aer"),mss_acm_cff_HNO3_aer);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_upt_cff_HNO3_aer"),mss_upt_cff_HNO3_aer);
    rcd=nco_put_var(nc_out,static_cast<std::string>("q_HNO3_gas"),q_HNO3_gas);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rxr_HNO3_gas_aer_vmr"),rxr_HNO3_gas_aer_vmr); delete []rxr_HNO3_gas_aer_vmr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rxr_HNO3_gas_aer_vmr_ttl"),rxr_HNO3_gas_aer_vmr_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rxrc_HNO3_aer"),rxrc_HNO3_aer); delete []rxrc_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rxrc_HNO3_aer_ttl"),rxrc_HNO3_aer_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_rxrc_HNO3_aer"),tau_rxrc_HNO3_aer); delete []tau_rxrc_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_gpd_HNO3_aer"),tau_gpd_HNO3_aer); delete []tau_gpd_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_ntf_eqm_HNO3_aer"),tau_ntf_eqm_HNO3_aer); delete []tau_ntf_eqm_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_aqs_dss_HNO3_aer"),tau_aqs_dss_HNO3_aer); delete []tau_aqs_dss_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_aqs_dff_HNO3_aer"),tau_aqs_dff_HNO3_aer); delete []tau_aqs_dff_HNO3_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlc_mwb_HNO3"),vlc_mwb_HNO3);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vmr_HNO3_gas"),vmr_HNO3_gas);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vnt_crc_aer"),vnt_crc_aer); delete []vnt_crc_aer;
  } // endif true

  if(mie_flg){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"abs_fsh",NC_FLOAT,1,dmn_sz,"long_name","Absorption efficiency","units","fraction"},
      {0,"ext_fsh",NC_FLOAT,1,dmn_sz,"long_name","Extinction efficiency","units","fraction"},
      {0,"bck_fsh",NC_FLOAT,1,dmn_sz,"long_name","Backscattering efficiency","units","fraction"},
      {0,"sca_fsh",NC_FLOAT,1,dmn_sz,"long_name","Scattering efficiency","units","fraction"},
      {0,"abs_spc",NC_FLOAT,1,dmn_sz,"long_name","Specific absorption","units","meter2 kilogram-1"},
      {0,"ext_spc",NC_FLOAT,1,dmn_sz,"long_name","Specific extinction","units","meter2 kilogram-1"},
      {0,"bck_spc",NC_FLOAT,1,dmn_sz,"long_name","Specific backscattering","units","meter2 kilogram-1"},
      {0,"sca_spc",NC_FLOAT,1,dmn_sz,"long_name","Specific scattering","units","meter2 kilogram-1"},
      {0,"asm_prm_fsh",NC_FLOAT,1,dmn_sz,"long_name","Asymmetry parameter","units","fraction"},
      {0,"ss_alb_fsh",nco_xtyp,1,dmn_sz,"long_name","Single scattering albedo","units","fraction"},
      {0,"ss_co_alb_fsh",nco_xtyp,1,dmn_sz,"long_name","Single scattering co-albedo","units","fraction"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_fsh"),abs_fsh); delete []abs_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_fsh"),ext_fsh); delete []ext_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_fsh"),bck_fsh); delete []bck_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_fsh"),sca_fsh); delete []sca_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_spc"),abs_spc); delete []abs_spc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_spc"),ext_spc); delete []ext_spc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_spc"),bck_spc); delete []bck_spc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_spc"),sca_spc); delete []sca_spc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("asm_prm_fsh"),asm_prm_fsh); delete []asm_prm_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ss_alb_fsh"),ss_alb_fsh); delete []ss_alb_fsh;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ss_co_alb_fsh"),ss_co_alb_fsh); delete []ss_co_alb_fsh;
  } // end if true

  if(mie_flg){ // New scope to hide re-declaration of var_mtd[] structure
    var_mtd_sct var_mtd[]={
      {0,"abs_cff_mss",nco_xtyp,1,dmn_wvl,"long_name","Mass absorption coefficient","units","meter2 kilogram-1"},
      {0,"abs_cff_mss_bb_LW",NC_FLOAT,0,dmn_scl,"long_name","Broadband longwave mass absorption coefficient","units","meter2 kilogram-1"},
      {0,"abs_cff_dst",nco_xtyp,1,dmn_wvl,"long_name","Distance absorption coefficient","units","meter-1"},
      {0,"abs_cff_vlm",nco_xtyp,1,dmn_wvl,"long_name","Volume absorption coefficient","units","meter2 meter-3"},
      {0,"abs_fct_MaS99",nco_xtyp,1,dmn_wvl,"long_name","Absorption enhancement for inclusions in weakly-absorbing spheres","units","fraction"},
      {0,"abs_fsh_ffc",nco_xtyp,1,dmn_wvl,"long_name","Effective absorption efficiency","units","fraction"},
      {0,"ang_xpn",NC_FLOAT,1,dmn_wvl,"long_name","Angstrom exponent","units","fraction"},
      {0,"asm_prm",nco_xtyp,1,dmn_wvl,"long_name","Asymmetry parameter","units","fraction"},
      {0,"ext_cff_mss",nco_xtyp,1,dmn_wvl,"long_name","Mass extinction coefficient","units","meter2 kilogram-1"},
      {0,"bck_cff_mss",nco_xtyp,1,dmn_wvl,"long_name","Mass backscattering coefficient","units","meter2 kilogram-1"},
      {0,"ext_cff_dst",nco_xtyp,1,dmn_wvl,"long_name","Distance extinction coefficient","units","meter-1"},
      {0,"ext_cff_vlm",nco_xtyp,1,dmn_wvl,"long_name","Volume extinction coefficient","units","meter2 meter-3"},
      {0,"bck_cff_dst",nco_xtyp,1,dmn_wvl,"long_name","Distance backscattering coefficient","units","meter-1"},
      {0,"bck_cff_vlm",nco_xtyp,1,dmn_wvl,"long_name","Volume backscattering coefficient","units","meter2 meter-3"},
      {0,"ext_fsh_ffc",nco_xtyp,1,dmn_wvl,"long_name","Effective extinction efficiency","units","fraction"},
      {0,"bck_fsh_ffc",nco_xtyp,1,dmn_wvl,"long_name","Effective backscattering efficiency","units","fraction"},
      {0,"msv_1km",NC_FLOAT,1,dmn_wvl,"long_name","Emissivity of 1 km column","units","fraction"},
      {0,"sca_cff_mss",nco_xtyp,1,dmn_wvl,"long_name","Mass scattering coefficient","units","meter2 kilogram-1"},
      {0,"sca_cff_dst",nco_xtyp,1,dmn_wvl,"long_name","Distance scattering coefficient","units","meter-1"},
      {0,"sca_cff_vlm",nco_xtyp,1,dmn_wvl,"long_name","Volume scattering coefficient","units","meter2 meter-3"},
      {0,"sca_fsh_ffc",nco_xtyp,1,dmn_wvl,"long_name","Effective scattering efficiency","units","fraction"},
      {0,"ss_alb",nco_xtyp,1,dmn_wvl,"long_name","Single scattering albedo","units","fraction"},
      {0,"ss_co_alb",nco_xtyp,1,dmn_wvl,"long_name","Single scattering co-albedo","units","fraction"},
      {0,"tau_abs",NC_FLOAT,1,dmn_wvl,"long_name","Absorption optical depth","units","fraction"},
      {0,"tau_ext",NC_FLOAT,1,dmn_wvl,"long_name","Extinction optical depth","units","fraction"},
      {0,"tau_sct",NC_FLOAT,1,dmn_wvl,"long_name","Scattering optical depth","units","fraction"},
      {0,"vsb_vsb",nco_xtyp,0,dmn_scl,"long_name","Visibility at 0.55 um","units","meter"}
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_cff_mss"),abs_cff_mss); delete []abs_cff_mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_cff_mss_bb_LW"),abs_cff_mss_bb_LW);
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_cff_dst"),abs_cff_dst); delete []abs_cff_dst;
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_cff_vlm"),abs_cff_vlm); delete []abs_cff_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_fct_MaS99"),abs_fct_MaS99); delete []abs_fct_MaS99;
    rcd=nco_put_var(nc_out,static_cast<std::string>("abs_fsh_ffc"),abs_fsh_ffc); delete []abs_fsh_ffc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ang_xpn"),ang_xpn); delete []ang_xpn;
    rcd=nco_put_var(nc_out,static_cast<std::string>("asm_prm"),asm_prm); delete []asm_prm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_cff_mss"),ext_cff_mss); delete []ext_cff_mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_cff_mss"),bck_cff_mss); delete []bck_cff_mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_cff_dst"),ext_cff_dst); delete []ext_cff_dst;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_cff_vlm"),ext_cff_vlm); delete []ext_cff_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_cff_dst"),bck_cff_dst); delete []bck_cff_dst;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_cff_vlm"),bck_cff_vlm); delete []bck_cff_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ext_fsh_ffc"),ext_fsh_ffc); delete []ext_fsh_ffc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("bck_fsh_ffc"),bck_fsh_ffc); delete []bck_fsh_ffc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("msv_1km"),msv_1km); delete []msv_1km;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_cff_mss"),sca_cff_mss); delete []sca_cff_mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_cff_dst"),sca_cff_dst); delete []sca_cff_dst;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_cff_vlm"),sca_cff_vlm); delete []sca_cff_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sca_fsh_ffc"),sca_fsh_ffc); delete []sca_fsh_ffc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ss_alb"),ss_alb); delete []ss_alb;
    rcd=nco_put_var(nc_out,static_cast<std::string>("ss_co_alb"),ss_co_alb); delete []ss_co_alb;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_abs"),tau_abs); delete []tau_abs;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_ext"),tau_ext); delete []tau_ext;
    rcd=nco_put_var(nc_out,static_cast<std::string>("tau_sct"),tau_sct); delete []tau_sct;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vsb_vsb"),vsb_vsb);
  } // end if true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    // fxm: 20020117: ccmalloc reports major leaks here
    var_mtd_sct var_mtd[]={
      {0,"asp_rat_hxg",NC_FLOAT,1,dmn_sz,"long_name","Hexagonal prism aspect ratio","units","fraction"},
      {0,"cnc_vts",NC_FLOAT,1,dmn_sz,"long_name","Number concentration of equal V/S spheres","units","fraction"},
      {0,"dmt_hxg",NC_FLOAT,1,dmn_sz,"long_name","Length c of hexagonal prism","units","fraction"},
      {0,"nbr_vts_per_hxg",NC_FLOAT,1,dmn_sz,"long_name","Number equal V/S spheres per hexagonal prism","units","fraction"},
      {0,"rds_ctr_vts",NC_FLOAT,1,dmn_sz,"long_name","Radius at bin center of equal V/S spheres","units","fraction"},
      {0,"rds_hxg",NC_FLOAT,1,dmn_sz,"long_name","Half-width a of basal face of hexagonal prism","units","fraction"},
      {0,"xsa_vts",NC_FLOAT,1,dmn_sz,"long_name","Equivalent sphere cross-sectional area","units","fraction"},
      {0,"vlm_vts",NC_FLOAT,1,dmn_sz,"long_name","Equivalent sphere volume","units","fraction"},
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("asp_rat_hxg"),asp_rat_hxg); delete []asp_rat_hxg;
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_vts"),cnc_vts); delete []cnc_vts;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_hxg"),dmt_hxg); delete []dmt_hxg;
    rcd=nco_put_var(nc_out,static_cast<std::string>("nbr_vts_per_hxg"),nbr_vts_per_hxg); delete []nbr_vts_per_hxg;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_ctr_vts"),rds_ctr_vts); delete []rds_ctr_vts;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_hxg"),rds_hxg); delete []rds_hxg;
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_vts"),xsa_vts); delete []xsa_vts;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_vts"),vlm_vts); delete []vlm_vts;
  } // end if true

  if(true){ // New scope to hide re-declaration of var_mtd[] structure
    // fxm: 20020117: ccmalloc reports major leaks here
    var_mtd_sct var_mtd[]={
      {0,"asp_rat_lps",NC_FLOAT,1,dmn_sz,"long_name","Ellipsoidal aspect ratio","units","fraction"},
      {0,"bnd_nbr",NC_INT,0,dmn_scl,"long_name","Number of sub-bands per output band","units","number"},
      {0,"cnc",NC_FLOAT,1,dmn_sz,"long_name","Number concentration","units","meter-3"},
      {0,"cnc_nbr_anl",NC_FLOAT,0,dmn_scl,"long_name","Number concentration analytic","units","meter-3"},
      {0,"cnc_nbr_rsl",NC_FLOAT,0,dmn_scl,"long_name","Number concentration resolved","units","meter-3"},
      {0,"cnc_nbr_rsl_frc",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of number concentration","units","fraction"},
      {0,"cnc_spc_anl",NC_FLOAT,0,dmn_scl,"long_name","Specific number concentration analytic","units","number kilogram-1"},
      {0,"cnc_spc_rsl",NC_FLOAT,0,dmn_scl,"long_name","Specific number concentration resolved","units","number kilogram-1"},
      {0,"dmt_aer",nco_xtyp,1,dmn_sz,"long_name","Aerodynamic diameter","units","meter"},
      {0,"dmt_ctr",nco_xtyp,1,dmn_sz,"long_name","Diameter at bin center","units","meter"},
      {0,"dmt_ctr_ctr",nco_xtyp,0,dmn_scl,"long_name","Mean grid diameter","units","meter"},
      {0,"dmt_dlt",nco_xtyp,1,dmn_sz,"long_name","Width of diameter bin","units","meter"},
      {0,"dmt_eqv_sfc",nco_xtyp,1,dmn_sz,"long_name","Diameter of sphere with same surface area","units","meter"},
      {0,"dmt_eqv_vlm",nco_xtyp,1,dmn_sz,"long_name","Diameter of sphere with same volume","units","meter"},
      {0,"dmt_grd",nco_xtyp,1,dmn_sz_grd,"long_name","Diameter grid","units","meter"},
      {0,"dmt_max",nco_xtyp,1,dmn_sz,"long_name","Maximum diameter in bin","units","meter"},
      {0,"dmt_max_max",nco_xtyp,0,dmn_scl,"long_name","Maximum grid diameter","units","meter"},
      {0,"dmt_min",nco_xtyp,1,dmn_sz,"long_name","Minimum diameter in bin","units","meter"},
      {0,"dmt_min_min",nco_xtyp,0,dmn_scl,"long_name","Minimum grid diameter","units","meter"},
      {0,"dmt_nma",nco_xtyp,0,dmn_scl,"long_name","Number median diameter analytic","units","meter"},
      {0,"dmt_nmr",nco_xtyp,0,dmn_scl,"long_name","Number median diameter resolved","units","meter"},
      {0,"dmt_nwa",nco_xtyp,0,dmn_scl,"long_name","Number weighted mean diameter analytic","units","meter"},
      {0,"dmt_nwr",nco_xtyp,0,dmn_scl,"long_name","Arithmetic mean diameter resolved","units","meter"},
      {0,"dmt_sma",nco_xtyp,0,dmn_scl,"long_name","Surface median diameter analytic","units","meter"},
      {0,"dmt_smr",nco_xtyp,0,dmn_scl,"long_name","Surface median diameter resolved","units","meter"},
      {0,"dmt_stk",nco_xtyp,1,dmn_sz,"long_name","Stokes diameter","units","meter"},
      {0,"dmt_mjr",nco_xtyp,1,dmn_sz,"long_name","Major axis of ellipsoid","units","meter"},
      {0,"dmt_mnr",nco_xtyp,1,dmn_sz,"long_name","Minor axis of ellipsoid","units","meter"},
      {0,"dmt_swa",nco_xtyp,0,dmn_scl,"long_name","Surface area weighted mean diameter analytic","units","meter"},
      {0,"dmt_swr",nco_xtyp,0,dmn_scl,"long_name","Surface area weighted mean diameter resolved","units","meter"},
      {0,"dmt_vma",nco_xtyp,0,dmn_scl,"long_name","Volume median diameter analytic","units","meter"},
      {0,"dmt_vmr",nco_xtyp,0,dmn_scl,"long_name","Volume median diameter resolved","units","meter"},
      {0,"dmt_vwa",nco_xtyp,0,dmn_scl,"long_name","Volume weighted mean diameter analytic","units","meter"},
      {0,"dmt_vwr",nco_xtyp,0,dmn_scl,"long_name","Mass weighted mean diameter resolved","units","meter"},
      {0,"dns_cor",NC_FLOAT,0,dmn_scl,"long_name","Density of core","units","kilogram meter-3"},
      {0,"dns_mdm",NC_FLOAT,0,dmn_scl,"long_name","Density of medium","units","kilogram meter-3"},
      {0,"dns_mnt",NC_FLOAT,0,dmn_scl,"long_name","Density of mantle","units","kilogram meter-3"},
      {0,"dns_mtx",NC_FLOAT,0,dmn_scl,"long_name","Density of matrix","units","kilogram meter-3"},
      {0,"dns_ncl",NC_FLOAT,0,dmn_scl,"long_name","Density of inclusion","units","kilogram meter-3"},
      {0,"dns_prt",NC_FLOAT,0,dmn_scl,"long_name","Density of particle","units","kilogram meter-3"},
      {0,"dst",NC_FLOAT,1,dmn_sz,"long_name","Number distribution","units","meter-3 meter-1"},
      {0,"dst_mss",NC_FLOAT,1,dmn_sz,"long_name","Mass distribution","units","kilogram meter-3 meter-1"},
      {0,"dst_rds",NC_FLOAT,1,dmn_sz,"long_name","Radius distribution","units","meter meter-3 meter-1"},
      {0,"dst_sfc",NC_FLOAT,1,dmn_sz,"long_name","Surface area distribution","units","meter2 meter-3 meter-1"},
      {0,"dst_vlm",NC_FLOAT,1,dmn_sz,"long_name","Volume distribution","units","meter3 meter-3 meter-1"},
      {0,"dst_xsa",NC_FLOAT,1,dmn_sz,"long_name","Cross-Sectional area distribution","units","meter2 meter-3 meter-1"},
      {0,"flx_IR_frc",nco_xtyp,1,dmn_wvl,"long_name","Fraction of terrestrial flux in band","units","fraction"},
      {0,"flx_IR_frc_ttl",nco_xtyp,0,dmn_scl,"long_name","Total terrestrial flux in all bands, normalized by blackbody flux","units","fraction"},
      {0,"flx_slr",nco_xtyp,1,dmn_wvl,"long_name","Absolute solar flux in band","units","watt meter-2"},
      {0,"flx_slr_frc",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux in band","units","fraction"},
      {0,"flx_slr_frc_blr",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux at shorter wavelengths","units","fraction"},
      {0,"flx_slr_frc_LaN68",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Labs & Neckel 1968","units","fraction"},
      {0,"flx_slr_frc_ThD71",nco_xtyp,1,dmn_wvl,"long_name","Fraction of solar flux: Thekeakara & Drummond 1971","units","fraction"},
      {0,"flx_spc_slr",nco_xtyp,1,dmn_wvl,"long_name","Solar spectral flux","units","watt meter-2 meter-1"},
      {0,"flx_spc_slr_pht",nco_xtyp,1,dmn_wvl,"long_name","Solar spectral photon flux","units","photon meter-2 second-1 meter-1"},
      {0,"lng_foo",NC_INT,0,dmn_scl,"long_name","Intrinsic long temporary variable","units","fraction"},
      {0,"gsd_anl_dfl",NC_FLOAT,0,dmn_scl,"long_name","Geometric sta<ndard deviation","units","fraction"},
      {0,"nrg_pht",nco_xtyp,1,dmn_wvl,"long_name","Energy of photon at band center","units","joule photon-1"},
      {0,"idx_rfr_air_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of air, imaginary component","units","fraction"},
      {0,"idx_rfr_air_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of air, real component","units","fraction"},
      {0,"idx_rfr_ffc_img",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index of particle, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_rl",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index of particle, real component","units","fraction"},
      {0,"idx_rfr_ffc_brg_img",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, Bruggeman approximation, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_brg_rl",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, Bruggeman approximation, real component","units","fraction"},
      {0,"idx_rfr_ffc_pmr_img",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, partial molar refraction approximation, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_pmr_rl",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, partial molar refraction approximation, real component","units","fraction"},
      {0,"idx_rfr_ffc_mxg_img",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, Maxwell Garnett approximation, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_mxg_rl",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, Maxwell Garnett approximation, real component","units","fraction"},
      {0,"idx_rfr_ffc_vlw_img",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, volume-weighted approximation, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_vlw_rl",NC_FLOAT,1,dmn_wvl,"long_name","Effective refractive index, volume-weighted approximation, real component","units","fraction"},
      {0,"idx_rfr_ffc_wgt_img",NC_FLOAT,1,dmn_wvl,"long_name","Spectral flux-weighted effective refractive index, imaginary component","units","fraction"},
      {0,"idx_rfr_ffc_wgt_rl",NC_FLOAT,1,dmn_wvl,"long_name","Spectral flux-weighted effective refractive index, real component","units","fraction"},
      {0,"idx_rfr_cor_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of core, imaginary component","units","fraction"},
      {0,"idx_rfr_cor_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of core, real component","units","fraction"},
      {0,"idx_rfr_mdm_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of medium, imaginary component","units","fraction"},
      {0,"idx_rfr_mdm_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of medium, real component","units","fraction"},
      {0,"idx_rfr_mnt_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of mantle, imaginary component","units","fraction"},
      {0,"idx_rfr_mnt_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of mantle, real component","units","fraction"},
      {0,"idx_rfr_mtx_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of matrix, imaginary component","units","fraction"},
      {0,"idx_rfr_mtx_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of matrix, real component","units","fraction"},
      {0,"idx_rfr_ncl_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of inclusion, imaginary component","units","fraction"},
      {0,"idx_rfr_ncl_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of inclusion, real component","units","fraction"},
      {0,"idx_rfr_prt_img",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of particle, imaginary component","units","fraction"},
      {0,"idx_rfr_prt_rl",NC_FLOAT,1,dmn_wvl,"long_name","Refractive index of particle, real component","units","fraction"},
      {0,"mmr_rsl",NC_FLOAT,0,dmn_scl,"long_name","Aerosol mixing ratio resolved","units","kilogram kilogram-1"},
      {0,"mmw_prt",NC_FLOAT,0,dmn_scl,"long_name","Mean molecular weight of particle","units","kilogram mole-1"},
      {0,"mss",NC_FLOAT,1,dmn_sz,"long_name","Mass","units","kilogram"},
      {0,"mss_anl",NC_FLOAT,0,dmn_scl,"long_name","Mass concentration analytic","units","kilogram meter-3"},
      {0,"rds_frc_cor",NC_FLOAT,0,dmn_scl,"long_name","Radius fraction in core","units","fraction"},
      {0,"rds_frc_ncl",NC_FLOAT,0,dmn_scl,"long_name","Radius fraction in inclusion","units","fraction"},
      {0,"mss_frc_cor",NC_FLOAT,0,dmn_scl,"long_name","Mass fraction in core","units","fraction"},
      {0,"mss_frc_ncl",NC_FLOAT,0,dmn_scl,"long_name","Mass fraction in inclusion","units","fraction"},
      {0,"vlm_frc_cor",NC_FLOAT,0,dmn_scl,"long_name","Volume fraction in core","units","fraction"},
      {0,"vlm_frc_mnt",NC_FLOAT,0,dmn_scl,"long_name","Volume fraction in mantle","units","fraction"},
      {0,"vlm_frc_mtx",NC_FLOAT,0,dmn_scl,"long_name","Volume fraction in matrix","units","fraction"},
      {0,"vlm_frc_ncl",NC_FLOAT,0,dmn_scl,"long_name","Volume fraction in inclusion","units","fraction"},
      {0,"vlm_frc_prt",NC_FLOAT,0,dmn_scl,"long_name","Volume fraction in particle","units","fraction"},
      {0,"mss_prt_rsl",NC_FLOAT,1,dmn_sz,"long_name","Mass concentration of smaller particles","units","kilogram meter-3"},
      {0,"mss_prt_rsl_frc",NC_FLOAT,1,dmn_sz,"long_name","Fraction of mass concentration from smaller particles","units","fraction"},
      {0,"mss_rsl",NC_FLOAT,0,dmn_scl,"long_name","Mass concentration resolved","units","kilogram meter-3"},
      {0,"nbr_prt_rsl",NC_FLOAT,1,dmn_sz,"long_name","Number concentration of smaller particles","units","meter-3"},
      {0,"nbr_prt_rsl_frc",NC_FLOAT,1,dmn_sz,"long_name","Fraction of number concentration from smaller particles","units","fraction"},
      {0,"rds_ctr",nco_xtyp,1,dmn_sz,"long_name","Radius at bin center","units","meter"},
      {0,"rds_ctr_ctr",nco_xtyp,0,dmn_scl,"long_name","Mean grid radius","units","meter"},
      {0,"rds_dlt",nco_xtyp,1,dmn_sz,"long_name","Width of radius bin","units","meter"},
      {0,"rds_grd",nco_xtyp,1,dmn_sz_grd,"long_name","Radius grid","units","meter"},
      {0,"rds_max",nco_xtyp,1,dmn_sz,"long_name","Maximum radius in bin","units","meter"},
      {0,"rds_max_max",nco_xtyp,0,dmn_scl,"long_name","Maximum grid radius","units","meter"},
      {0,"rds_min",nco_xtyp,1,dmn_sz,"long_name","Minimum radius in bin","units","meter"},
      {0,"rds_min_min",nco_xtyp,0,dmn_scl,"long_name","Minimum grid radius","units","meter"},
      {0,"rds_nma",nco_xtyp,0,dmn_scl,"long_name","Number median radius analytic","units","meter"},
      {0,"rds_nmr",nco_xtyp,0,dmn_scl,"long_name","Number median radius resolved","units","meter"},
      {0,"rds_nwa",nco_xtyp,0,dmn_scl,"long_name","Number weighted mean radius analytic","units","meter"},
      {0,"rds_nwr",nco_xtyp,0,dmn_scl,"long_name","Arithmetic mean radius resolved","units","meter"},
      {0,"rds_sma",nco_xtyp,0,dmn_scl,"long_name","Surface median radius analytic","units","meter"},
      {0,"rds_smr",nco_xtyp,0,dmn_scl,"long_name","Surface median radius resolved","units","meter"},
      {0,"rds_swa",nco_xtyp,0,dmn_scl,"long_name","Surface area weighted mean radius analytic","units","meter"},
      {0,"rds_swr",nco_xtyp,0,dmn_scl,"long_name","Surface area weighted mean radius resolved","units","meter"},
      {0,"rds_vma",nco_xtyp,0,dmn_scl,"long_name","Volume median radius analytic","units","meter"},
      {0,"rds_vmr",nco_xtyp,0,dmn_scl,"long_name","Volume median radius resolved","units","meter"},
      {0,"rds_vwa",nco_xtyp,0,dmn_scl,"long_name","Volume weighted mean radius analytic","units","meter"},
      {0,"rds_vwr",nco_xtyp,0,dmn_scl,"long_name","Mass weighted mean radius resolved","units","meter"},
      {0,"sfc",NC_FLOAT,1,dmn_sz,"long_name","Surface area","units","meter2"},
      {0,"sfc_anl",NC_FLOAT,0,dmn_scl,"long_name","Surface area concentration analytic","units","meter2 meter-3"},
      {0,"sfc_prt_rsl",NC_FLOAT,1,dmn_sz,"long_name","Surface area concentration of smaller particles","units","meter2 meter-3"},
      {0,"sfc_prt_rsl_frc",NC_FLOAT,1,dmn_sz,"long_name","Fraction of surface area concentration from smaller particles","units","fraction"},
      {0,"sfc_rsl",NC_FLOAT,0,dmn_scl,"long_name","Surface area concentration resolved","units","meter2 meter-3"},
      {0,"sfc_rsl_frc",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of surface area concentration","units","fraction"},
      {0,"sfc_spc_anl",NC_FLOAT,0,dmn_scl,"long_name","Specific surface area analytic","units","meter2 kilogram-1"},
      {0,"sfc_spc_rsl",NC_FLOAT,0,dmn_scl,"long_name","Specific Surface area resolved","units","meter2 kilogram-1"},
      {0,"slr_cst",NC_FLOAT,0,dmn_scl,"long_name","Solar constant","units","watt meter-2"},
      {0,"spc_heat_prt",NC_FLOAT,0,dmn_scl,"long_name","Specific heat capacity of particle","units","joule kilogram-1 kelvin-1"},
      {0,"sz_ctr",nco_xtyp,1,dmn_sz,"long_name","Size at bin center","units","meter"},
      {0,"sz_dbg",nco_xtyp,0,dmn_scl,"long_name","Debugging size","units","meter"},
      {0,"sz_dlt",nco_xtyp,1,dmn_sz,"long_name","Width of size bin","units","meter"},
      {0,"sz_max",nco_xtyp,1,dmn_sz,"long_name","Maximum size in bin","units","meter"},
      {0,"sz_min",nco_xtyp,1,dmn_sz,"long_name","Minimum size in bin","units","meter"},
      {0,"sz_nbr",NC_INT,0,dmn_scl,"long_name","Number of particle size bins","units","number"},
      {0,"sz_prm_dbg",NC_FLOAT,0,dmn_scl,"long_name","Size parameter at debug size, wavelength","units","meter meter-1"},
      {0,"sz_prm_rds_min_wvl_max",NC_FLOAT,0,dmn_scl,"long_name","Size parameter, smallest particle, longest wavelength","units","meter meter-1"},
      {0,"sz_prm_rds_max_wvl_min",NC_FLOAT,0,dmn_scl,"long_name","Size parameter, largest particle, shortest wavelength","units","meter meter-1"},
      {0,"sz_prm_rsn_avg",NC_FLOAT,0,dmn_scl,"long_name","Size parameter, mean computational resolution","units","meter meter-1"},
      {0,"sz_prm_swa",NC_FLOAT,1,dmn_wvl,"long_name","Size parameter at rds_swa","units","meter meter-1"},
      {0,"vlm",NC_FLOAT,1,dmn_sz,"long_name","Volume","units","meter3"},
      {0,"vlm_anl",NC_FLOAT,0,dmn_scl,"long_name","Volume concentration analytic","units","meter3 meter-3"},
      {0,"vlm_prt_rsl",NC_FLOAT,1,dmn_sz,"long_name","Volume concentration of smaller particles","units","meter3 meter-3"},
      {0,"vlm_prt_rsl_frc",NC_FLOAT,1,dmn_sz,"long_name","Fraction of volume concentration from smaller particles","units","fraction"},
      {0,"vlm_rsl",NC_FLOAT,0,dmn_scl,"long_name","Volume concentration resolved","units","meter3 meter-3"},
      {0,"vlm_rsl_frc",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of volume concentration","units","fraction"},
      {0,"vlm_spc_anl",NC_FLOAT,0,dmn_scl,"long_name","Specific volume analytic","units","meter3 kilogram-1"},
      {0,"vlm_spc_rsl",NC_FLOAT,0,dmn_scl,"long_name","Specific volume resolved","units","meter3 kilogram-1"},
      {0,"wvl_bnd_sz_nbr",NC_INT,0,dmn_scl,"long_name","Number of wavelength/band/size loops","units","number"},
      {0,"wvl_ctr",nco_xtyp,1,dmn_wvl,"long_name","Wavelength at band center","units","meter"},
      {0,"wvl_dbg",nco_xtyp,0,dmn_scl,"long_name","Debugging wavelength","units","meter"},
      {0,"wvl_dlt",nco_xtyp,1,dmn_wvl,"long_name","Bandwidth","units","meter"},
      {0,"wvl_max",nco_xtyp,1,dmn_wvl,"long_name","Maximum wavelength in band","units","meter"},
      {0,"wvl_min",nco_xtyp,1,dmn_wvl,"long_name","Minimum wavelength in band","units","meter"},
      {0,"wvl_nbr",NC_INT,0,dmn_scl,"long_name","Number of wavelength bands","units","number"},
      {0,"wvl_grd_nbr",NC_INT,0,dmn_scl,"long_name","Number of wavelength interfaces","units","number"},
      {0,"wvl_wgt",nco_xtyp,1,dmn_wvl,"long_name","Solar flux-weighted wavelength","units","meter"},
      {0,"wvn",nco_xtyp,1,dmn_wvl,"long_name","Band nominal wavenumber","units","centimeter-1"},
      {0,"wvn_ctr",nco_xtyp,1,dmn_wvl,"long_name","Wavenumber at band center","units","centimeter-1"},
      {0,"wvn_dlt",nco_xtyp,1,dmn_wvl,"long_name","Bandwidth","units","centimeter-1"},
      {0,"wvn_grd",nco_xtyp,1,dmn_wvl_grd,"long_name","Wavenumber grid","units","centimeter-1"},
      {0,"wvn_max",nco_xtyp,1,dmn_wvl,"long_name","Band maximum wavenumber","units","centimeter-1"},
      {0,"wvn_min",nco_xtyp,1,dmn_wvl,"long_name","Band minimum wavenumber","units","centimeter-1"},
      {0,"xsa",NC_FLOAT,1,dmn_sz,"long_name","Cross-sectional area","units","meter2"},
      {0,"xsa_anl",NC_FLOAT,0,dmn_scl,"long_name","Cross-sectional area concentration analytic","units","meter2 meter-3"},
      {0,"xsa_prt_rsl",NC_FLOAT,1,dmn_sz,"long_name","Cross-sectional area concentration of smaller particles","units","meter2 meter-3"},
      {0,"xsa_prt_rsl_frc",NC_FLOAT,1,dmn_sz,"long_name","Fraction of cross-sectional area concentration from smaller particles","units","fraction"},
      {0,"xsa_rsl",NC_FLOAT,0,dmn_scl,"long_name","Cross-sectional area concentration resolved","units","meter2 meter-3"},
      {0,"xsa_rsl_frc",NC_FLOAT,0,dmn_scl,"long_name","Resolved fraction of cross-sectional area concentration","units","fraction"},
      {0,"xsa_spc_anl",NC_FLOAT,0,dmn_scl,"long_name","Specific cross-sectional area analytic","units","meter2 kilogram-1"},
      {0,"xsa_spc_rsl",NC_FLOAT,0,dmn_scl,"long_name","Specific cross-sectional area resolved","units","meter2 kilogram-1"},
    }; // end var_mtd_sct var_mtd[]
    const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
    rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
    
    // After writing, delete arrays 
    rcd=nco_put_var(nc_out,static_cast<std::string>("asp_rat_lps"),asp_rat_lps); delete []asp_rat_lps;
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc"),cnc); delete []cnc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_anl"),cnc_nbr_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_rsl"),cnc_nbr_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_nbr_rsl_frc"),cnc_nbr_rsl_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_spc_anl"),cnc_spc_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_spc_rsl"),cnc_spc_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_aer"),dmt_aer); delete []dmt_aer;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_ctr"),dmt_ctr); delete []dmt_ctr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_ctr_ctr"),dmt_ctr_ctr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_dlt"),dmt_dlt); delete []dmt_dlt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_eqv_sfc"),dmt_eqv_sfc); delete []dmt_eqv_sfc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_eqv_vlm"),dmt_eqv_vlm); delete []dmt_eqv_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_grd"),dmt_grd); delete []dmt_grd;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_max"),dmt_max); delete []dmt_max;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_max_max"),dmt_max_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_min"),dmt_min); delete []dmt_min;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_min_min"),dmt_min_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_nma"),dmt_nma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_nmr"),dmt_nmr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_nwa"),dmt_nwa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_nwr"),dmt_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_sma"),dmt_sma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_smr"),dmt_smr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_stk"),dmt_stk); delete []dmt_stk;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_mjr"),dmt_mjr); delete []dmt_mjr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_mnr"),dmt_mnr); delete []dmt_mnr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_swa"),dmt_swa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_swr"),dmt_swr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_vma"),dmt_vma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_vmr"),dmt_vmr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_vwa"),dmt_vwa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dmt_vwr"),dmt_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_cor"),dns_cor);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_mdm"),dns_mdm);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_mnt"),dns_mnt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_mtx"),dns_mtx);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_ncl"),dns_ncl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dns_prt"),dns_prt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst"),dst); delete []dst;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_mss"),dst_mss); delete []dst_mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_rds"),dst_rds); delete []dst_rds;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_sfc"),dst_sfc); delete []dst_sfc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_vlm"),dst_vlm); delete []dst_vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("dst_xsa"),dst_xsa); delete []dst_xsa;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_IR_frc"),flx_IR_frc); delete []flx_IR_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_IR_frc_ttl"),flx_IR_frc_ttl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr"),flx_slr); delete []flx_slr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr_frc"),flx_slr_frc); delete []flx_slr_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr_frc_blr"),flx_slr_frc_blr); delete []flx_slr_frc_blr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr_frc_LaN68"),flx_slr_frc_LaN68); delete []flx_slr_frc_LaN68;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_slr_frc_ThD71"),flx_slr_frc_ThD71); delete []flx_slr_frc_ThD71;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_spc_slr"),flx_spc_slr); delete []flx_spc_slr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("flx_spc_slr_pht"),flx_spc_slr_pht); delete []flx_spc_slr_pht;
    rcd=nco_put_var(nc_out,static_cast<std::string>("gsd_anl_dfl"),gsd_anl_dfl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("lng_foo"),lng_foo);
    rcd=nco_put_var(nc_out,static_cast<std::string>("nrg_pht"),nrg_pht); delete []nrg_pht;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_air_img"),idx_rfr_air_img); delete []idx_rfr_air_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_air_rl"),idx_rfr_air_rl); delete []idx_rfr_air_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_img"),idx_rfr_ffc_img); delete []idx_rfr_ffc_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_rl"),idx_rfr_ffc_rl); delete []idx_rfr_ffc_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_brg_img"),idx_rfr_ffc_brg_img); delete []idx_rfr_ffc_brg_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_brg_rl"),idx_rfr_ffc_brg_rl); delete []idx_rfr_ffc_brg_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_pmr_img"),idx_rfr_ffc_pmr_img); delete []idx_rfr_ffc_pmr_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_pmr_rl"),idx_rfr_ffc_pmr_rl); delete []idx_rfr_ffc_pmr_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_mxg_img"),idx_rfr_ffc_mxg_img); delete []idx_rfr_ffc_mxg_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_mxg_rl"),idx_rfr_ffc_mxg_rl); delete []idx_rfr_ffc_mxg_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_vlw_img"),idx_rfr_ffc_vlw_img); delete []idx_rfr_ffc_vlw_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_vlw_rl"),idx_rfr_ffc_vlw_rl); delete []idx_rfr_ffc_vlw_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_wgt_img"),idx_rfr_ffc_wgt_img); delete []idx_rfr_ffc_wgt_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ffc_wgt_rl"),idx_rfr_ffc_wgt_rl); delete []idx_rfr_ffc_wgt_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_cor_img"),idx_rfr_cor_img); delete []idx_rfr_cor_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_cor_rl"),idx_rfr_cor_rl); delete []idx_rfr_cor_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mdm_img"),idx_rfr_mdm_img); delete []idx_rfr_mdm_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mdm_rl"),idx_rfr_mdm_rl); delete []idx_rfr_mdm_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mnt_img"),idx_rfr_mnt_img); delete []idx_rfr_mnt_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mnt_rl"),idx_rfr_mnt_rl); delete []idx_rfr_mnt_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mtx_img"),idx_rfr_mtx_img); delete []idx_rfr_mtx_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_mtx_rl"),idx_rfr_mtx_rl); delete []idx_rfr_mtx_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ncl_img"),idx_rfr_ncl_img); delete []idx_rfr_ncl_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_ncl_rl"),idx_rfr_ncl_rl); delete []idx_rfr_ncl_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_prt_img"),idx_rfr_prt_img); delete []idx_rfr_prt_img;
    rcd=nco_put_var(nc_out,static_cast<std::string>("idx_rfr_prt_rl"),idx_rfr_prt_rl); delete []idx_rfr_prt_rl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mmr_rsl"),mmr_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mmw_prt"),mmw_prt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss"),mss); delete []mss;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_anl"),mss_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_frc_cor"),rds_frc_cor);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_frc_ncl"),rds_frc_ncl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_cor"),mss_frc_cor);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_frc_ncl"),mss_frc_ncl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_frc_cor"),vlm_frc_cor);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_frc_mnt"),vlm_frc_mnt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_frc_mtx"),vlm_frc_mtx);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_frc_ncl"),vlm_frc_ncl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_frc_prt"),vlm_frc_prt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_prt_rsl"),mss_prt_rsl); delete []mss_prt_rsl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_prt_rsl_frc"),mss_prt_rsl_frc); delete []mss_prt_rsl_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("mss_rsl"),mss_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("nbr_prt_rsl"),nbr_prt_rsl); delete []nbr_prt_rsl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("nbr_prt_rsl_frc"),nbr_prt_rsl_frc); delete []nbr_prt_rsl_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_ctr"),rds_ctr); delete []rds_ctr;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_ctr_ctr"),rds_ctr_ctr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_dlt"),rds_dlt); delete []rds_dlt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_grd"),rds_grd); delete []rds_grd;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_max"),rds_max); delete []rds_max;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_max_max"),rds_max_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_min"),rds_min); delete []rds_min;
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_min_min"),rds_min_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_nma"),rds_nma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_nmr"),rds_nmr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_nwa"),rds_nwa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_nwr"),rds_nwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_sma"),rds_sma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_smr"),rds_smr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_swa"),rds_swa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_swr"),rds_swr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_vma"),rds_vma);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_vmr"),rds_vmr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_vwa"),rds_vwa);
    rcd=nco_put_var(nc_out,static_cast<std::string>("rds_vwr"),rds_vwr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc"),sfc); delete []sfc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_anl"),sfc_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_prt_rsl"),sfc_prt_rsl); delete []sfc_prt_rsl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_prt_rsl_frc"),sfc_prt_rsl_frc); delete []sfc_prt_rsl_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_rsl"),sfc_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_rsl_frc"),sfc_rsl_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_spc_anl"),sfc_spc_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sfc_spc_rsl"),sfc_spc_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("slr_cst"),slr_cst);
    rcd=nco_put_var(nc_out,static_cast<std::string>("spc_heat_prt"),spc_heat_prt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_ctr"),sz_ctr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_dbg"),sz_dbg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_dlt"),sz_dlt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_max"),sz_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_min"),sz_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_prm_dbg"),sz_prm_dbg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_prm_rds_min_wvl_max"),sz_prm_rds_min_wvl_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_prm_rds_max_wvl_min"),sz_prm_rds_max_wvl_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_prm_rsn_avg"),sz_prm_rsn_avg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("sz_prm_swa"),sz_prm_swa); delete []sz_prm_swa;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm"),vlm); delete []vlm;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_anl"),vlm_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_prt_rsl"),vlm_prt_rsl); delete []vlm_prt_rsl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_prt_rsl_frc"),vlm_prt_rsl_frc); delete []vlm_prt_rsl_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_rsl"),vlm_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_rsl_frc"),vlm_rsl_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_spc_anl"),vlm_spc_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("vlm_spc_rsl"),vlm_spc_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_bnd_sz_nbr"),wvl_bnd_sz_nbr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_ctr"),wvl_ctr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_dbg"),wvl_dbg);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_dlt"),wvl_dlt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_max"),wvl_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_min"),wvl_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvl_wgt"),wvl_wgt); delete []wvl_wgt;
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn"),wvn);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_ctr"),wvn_ctr);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_dlt"),wvn_dlt);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_grd"),wvn_grd);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_max"),wvn_max);
    rcd=nco_put_var(nc_out,static_cast<std::string>("wvn_min"),wvn_min);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa"),xsa); delete []xsa;
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_anl"),xsa_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_prt_rsl"),xsa_prt_rsl); delete []xsa_prt_rsl;
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_prt_rsl_frc"),xsa_prt_rsl_frc); delete []xsa_prt_rsl_frc;
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_rsl"),xsa_rsl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_rsl_frc"),xsa_rsl_frc);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_spc_anl"),xsa_spc_anl);
    rcd=nco_put_var(nc_out,static_cast<std::string>("xsa_spc_rsl"),xsa_spc_rsl);
  } // end if true
  
  // Close output file
  if(dbg_lvl > dbg_fl) std::cout << "Output file has " << nco_inq_ndims(nc_out) << " dimensions and " << nco_inq_nvars(nc_out) << " variables" << std::endl;
  rcd=nco_close(nc_out); // [fnc] Close netCDF file
    
  // Wrap up logfile
  if(true){
    std::cout << "Output information: " << std::endl;
    std::cout << "  File: " << fl_out << std::endl;
    std::cout << "  ncks: ncks -C -F -H -u -d wvl,0.5e-6 -v ext_cff_mss " << fl_out << std::endl;
  } /* end if true */

  if(abc_flg){
    // Alphabetize variables in output file
    std::string drc_ncks("/usr/local/bin/"); // [sng] Alternate directory for ncks
    std::string cmd_which("which ncks > /dev/null"); // [sng] System command string
    std::string cmd_ncks("ncks -h -O "+fl_out+" "+fl_out); // [sng] System command string
    rcd=std::system(cmd_which.c_str()); // 'which' returns number of failed arguments
    // In LoadLeveler, which finds ncks on User's path, but system(ncks) fails because PATH is minimal
    // Insert alternative path to ncks if which fails
    if(rcd) cmd_ncks.insert(0,drc_ncks); 
    rcd=std::system(cmd_ncks.c_str());
    // Having an unsorted netCDF output file is not fatal
    if(rcd) wrn_prn(sbr_nm," ncks command to sort output file failed");
    rcd=0; // [enm] Return code
  } // endif abc

  if(rcd) err_prn(sbr_nm," rcd != 0 on exit");
  Exit_gracefully();
} // end main() 
