# $Id$ -*-Fundamental-*-

Purpose: Document construction of Multi-Component Aerosols (MCAs)

# Organization
1. Frequently used environment variables
2. Standard, single-component dust aerosols: 
3. MCA 98% SiO2 2% FeO3 using Maxwell Garnett approximation
   3a. 
   3b. 
   3c. 
   3d. STD03 derive/interpolate relatively bright dust imaginary refractive index
   3e. Mineralogically consistent dust with BGC OEM (3.5% Fe)
   3f. Create CAM optical properties from Working Blend #1 (See above)
   3g. DHE02 p. 605 fgr 4 Attempt to capture g for dust spheroids @ 440 nm
4. Equal-V/S treatment of hexagonal prisms (NGW03)
5. Absorbing inclusions in weakly-absorbing spheres (MaS99)
6. Volcanic aerosol for CCM3/CAM3
7. Internally mixed soot
8. LGGE snow
   8a. Reproduce snow results of JC Gallet
   8b. Physically plausible snowpack reflectance
   8c. Arctic clear-sky downwelling spectra
   8d. SWNB prognostic snow impurities
   8e. Correct LGGE lab measurements for horizontal effects
   8f. BC MAC uncertainty
   8g. Sample SSA from NIR reflectance
   8h. GDZ08 outliers
9. Phase functions

# 1. Frequently used environment variables
export DATA='/data/zender' # UCI
export DATA='/data/zender/zender' # UCI
export DATA='/ptmp/zender' # ESMF
export CASEID='mie16'
export SZ_SNG='rds_swa_10'

# 2. Standard single component aerosols
# Fe2O3 Hematite (same size distribution as saharan_dust)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=Fe2O3 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_Fe2O3.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm ${DATA}/aca/aer_Fe2O3.nc

# SiO2 Quartz (same size distribution as saharan_dust)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=SiO2 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_SiO2.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm ${DATA}/aca/aer_SiO2.nc

# saharan_dust PaG77 West Texas 
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=saharan_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_saharan_dust.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm ${DATA}/aca/aer_saharan_dust.nc

# Evaluate against SoT99 p. 9439 Plt. 6
# 3a. MCA mass fractions 98% SiO2 02% FeO3 Maxwell Garnett approximation
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3 --mss_frc_ncl=0.02 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_SiO2_Fe2O3_98_02.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/aca/aer_SiO2_Fe2O3_98_02.nc

# 3b1. MCA mass fractions 96% SiO2 04% FeO3 Maxwell Garnett approximation
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3 --mss_frc_ncl=0.04 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_SiO2_Fe2O3_96_04.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/aca/aer_SiO2_Fe2O3_96_04.nc
ncks -C -m -d wvl,0.331e-6 -d wvl,0.360e-6 -d wvl,0.440e-6 -d wvl,0.670e-6 -v idx_rfr_ffc_img ${DATA}/aca/aer_SiO2_Fe2O3_96_04.nc

# 3b2. MCA mass fractions 96% SiO2 04% FeO3 Maxwell Garnett approximation, V/S approximation
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3 --mss_frc_ncl=0.04 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 --vts_flg --asp_rat_hxg=1.5 ${DATA}/aca/aer_SiO2_Fe2O3_96_04_hxg.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/aca/aer_SiO2_Fe2O3_96_04_hxg.nc

# 3c. MCA mass fractions 90% SiO2 10% FeO3 Maxwell Garnett approximation
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3 --mss_frc_ncl=0.10 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_SiO2_Fe2O3_90_10.nc
ncks -m -d wvl,0.5e-6 -v ext_cff_mss,ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/aca/aer_SiO2_Fe2O3_90_10.nc

# 3d. STD03 derive/interpolate relatively bright imaginary refractive index for dust
# See file idx_rfr_STD03.cdl
# Read from STD03 graph: (wvl_mcr,k) ~ (0.4,0.0038), (0.5,0.002)
# wl(um)   real  source	imaginary    source
# 0.331    1.58  TBH02	  0.00654  Tor03_eqn (calculated w/polynomial fit eqn)
# 0.340    1.58  TBH02    0.00616  Tor03_est (direct estimate @exact wl)
# 0.360    1.57  TBH02    0.00528  Tor03_eqn (calculated w/polynomial fit eqn)
# 0.380    1.58  TBH02    0.00440  Tor03_eqn (calculated w/polynomial fit eqn)
# 0.550    1.56  TBH02	  0.00140  STD03_grf (read from graph)
# 0.630    1.56  Pat81    0.00100  STD03_grf (read from graph)
# DHE02 p. 602 recommend Cape Verde/Saudi Arabia as typical "pure dust"
# DHE02 p. 602 recommend Bahrain/Persian Gulf as mixed industrial/dust
# AERONET dust retrievals from DHE02 at Cape Verde
# wl(um)   real  	imaginary	source
# 0.440    1.48+/-0.05  0.0025+/-0.001 	DHE02
# 0.670    1.48		0.0007  	DHE02
# 0.870    1.48		0.0006  	DHE02
# 1.020    1.48		0.0006  	DHE02
# AERONET dust retrievals from DHE02 at Saudi Arabia
# wl(um)   real  	imaginary	source
# 0.440    1.56+/-0.03  0.0029+/-0.001 	DHE02
# 0.670    1.56		0.0013  	DHE02
# 0.870    1.56		0.0010  	DHE02
# 1.020    1.56		0.0010  	DHE02
# AERONET dust retrievals from DHE02 at Bahrain mixed industrial/representative
# wl(um)   real  	imaginary	source
# 0.440    1.55+/-0.03  0.0025+/-0.001 	DHE02
# 0.670    1.55		0.0014  	DHE02
# 0.870    1.55		0.0010  	DHE02
# 1.020    1.55		0.0010  	DHE02
# AERONET dust retrievals from DHE02 at Barbados
# wl(um)   real  	imaginary	source
# 0.440    1.43+/-0.03  0.0016+/-0.001 	DHE02
# 0.670    1.43		0.0024  	DHE02
# 0.870    1.43		0.0033  	DHE02
# 1.020    1.43		0.0038  	DHE02
# AERONET dust retrievals from DHE02 at Mongolia
# wl(um)   real  	imaginary	source
# 0.440    1.555+/-0.03  0.0037+/-0.001	DHE02
# 0.670    1.555	 0.0027  	DHE02
# 0.870    1.555	 0.0027  	DHE02
# 1.020    1.555	 0.0027  	DHE02
# Find effective medium approximation (EMA) mixture which captures these optics
# Promising compositions include:
# 0.99 SiO2, 0.01 Fe2O3_doccd
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3_doccd --mss_frc_ncl=0.01 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_grd=sensor --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/aer_SiO2_Fe2O3_foo.nc
ncks -C -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/aer_SiO2_Fe2O3_foo.nc

# 3e. Mineralogically consistent dust with BGC OEM (3.5% Fe)
# Dust is ~3.5% elemental Fe
# Assuming all dust is present as Fe2O3, then Fe is 2*55.8/159.7=70% of Fe2O3
# Hence dust must be 3.5/0.7=5% Fe2O3 in order to have 3.5% elemental Fe
mie -3 --dbg=1 --slr_spc=Labs --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2 --cmp_ncl=Fe2O3_doccd --mss_frc_ncl=0.05 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_grd=sensor --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/aer_SiO2_Fe2O3_foo.nc
ncks -C -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/aer_SiO2_Fe2O3_foo.nc

# EMA for 95% illite matrix /kaolinite/montmorillonite/quartz matrix and 5% hematite inclusions
# Using this kludge technique, specify N-1 volume fractions
# Volume fractions correspond, in order, to cor, mnt, ncl, prt
# Matrix (mtx) is determined by residual and is not separately specified
# Currently mie interprets these as volume fractions not mass fractions!
# Fe2O3 density is ~5260 kg m-3 so 5% Fe2O3 by weight must be ~2.5% Fe2O3 by volume
# This mixture is my best estimate of global mean mass composition
# It is much too bright optically
# NB: graphics/analysis script ~/ncl/spc_ln.ncl reads data from output file ${DATA}/mie/idx_rfr_mlt.nc
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.25 --ncl=0.25 --ncl=0.20 --ncl=0.025 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=illite --cmp_cor=kaolinite --cmp_mnt=montmorillonite --cmp_ncl=SiO2_avg_hitran96 --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --wvl_grd=sensor --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc
ncks -C -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/idx_rfr_mlt.nc

# Illite + hematite alone give very reasonable results up to 0.67 um. Barbados and Mongolia are darker red-ward of there.
mie -3 --dbg=1 --slr_spc=Labs --mss_frc_ncl=0.05 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=illite --cmp_ncl=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --wvl_grd=sensor --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Limestone + hematite
# limestone_krk is much darker in visible than limestone_dsp
mie -3 --dbg=1 --slr_spc=Labs --mss_frc_ncl=0.05 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=limestone_krk --cmp_ncl=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Illite + Limestone: 5% limestone looks good, stronger UV absorption would help
mie -3 --dbg=1 --slr_spc=Labs --mss_frc_ncl=0.05 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=illite --cmp_ncl=limestone_krk --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# 20060904: This combination asymptotes to k(0.5 um) ~ 0.0025 and asymptotes to k ~ 0.017 in NIR
# I dub this "walkabout blend" #1 Slightly less limestone
# Immortalized as dst_bln_20060904
# dst_bln_20060904 is used in SNICAR, CAM4, CAM5. dst_bln_20060907 _may_ fit obs better, but inclusion of BC in dust not justifiable
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.02 --ncl=0.25 --ncl=0.25 --ncl=0.004 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_bln_wlk1.nc
/bin/cp ${DATA}/mie/idx_rfr_bln_wlk1.nc ${DATA}/mie/idx_rfr_mlt.nc

# 20130119: Per Rachel Scanza's request, remove Fe2O3 from dst_bln_20060904 computed above
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.02 --ncl=0.25 --ncl=0.25 --ncl=0.000 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_dst_bln_20130119.nc
ncks -O -v idx_rfr_ffc_rl,idx_rfr_ffc_img ${DATA}/mie/idx_rfr_bln_20130119.nc ${DATA}/mie/idx_rfr_bln_20130119.nc

# 20060904: This combination asymptotes to k(0.5 um) ~ 0.0025 and asymptotes to k ~ 0.017 in NIR
# I dub this "walkabout blend" #2 Slightly more limestone
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.03 --ncl=0.25 --ncl=0.25 --ncl=0.004 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_bln_wlk2.nc

# 20060904: This combination asymptotes to k(0.5 um) ~ 0.0022 and asymptotes to k ~ 0.017 in NIR
# I dub this "walkabout blend" #3 Less limestone and DOCCD instead of Roush
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.02 --ncl=0.25 --ncl=0.25 --ncl=0.003 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_bln_wlk3.nc

# I dub this "walkabout blend" #4 Uses soot instead of limestone
# The soot yields near-constant NIR which better reproduces AERONET NIR SSA increasing trend and AERONET NIR n_k decreasing trend
# Immortalized as dst_bln_20060907
# Default inclusion of BC properties in dust not justifiable, so dst_bln_20060907 not used by default in SNICAR, CAM4, CAM5
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.0014 --ncl=0.25 --ncl=0.25 --ncl=0.004 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=lac_ChC90 --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_bln_wlk4.nc 

# Manually mix Quartz + Illite + Montmorillonite + limestone_krk + hematite
# 20060621: This combination asymptotes to k(0.5 um) ~ 0.0025 and asymptotes to k ~ 0.01 in NIR
# I dub this Working Blend #1 (summer solstice blend)
# Immortalized as dst_bln_20060621
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.01 --ncl=0.25 --ncl=0.25 --ncl=0.005 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Manually mix Illite + Montmorillonite + Quartz + limestone_krk + hematite
# 20060621: This combination asymptotes to k(0.5 um) ~ 0.003 and asymptotes to k ~ 0.02 in NIR
# I dub this Working Blend #2
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.03 --ncl=0.25 --ncl=0.25 --ncl=0.005 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=illite  --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=SiO2_avg_hitran96 --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Manually mix Illite + Montmorillonite + Quartz + limestone_krk + Fe2O3_avg_hitran96
# 20060621: This combination asymptotes to k(0.5 um) ~ 0.003 and asymptotes to k ~ 0.02 in NIR
# I dub this Working Blend #3
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.03 --ncl=0.25 --ncl=0.25 --ncl=0.0075 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=illite  --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=SiO2_avg_hitran96 --cmp_prt=Fe2O3_avg_hitran96 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Manually mix Quartz + Illite + Montmorillonite + limestone_krk + hematite
# 20060621: This combination asymptotes to k(0.5 um) ~ 0.003 and asymptotes to k ~ 0.01 in NIR
# I dub this Working Blend #4
mie -3 --dbg=1 --slr_spc=Labs --ncl=0.005 --ncl=0.25 --ncl=0.25 --ncl=0.007 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96  --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=1 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# Compare calcite and gypsum to LQB93 and limestone to QOL78: Fit appears perfect
mie -3 --dbg=1 --slr_spc=Labs --mss_frc_ncl=0.50 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=calcite --cmp_ncl=gypsum --cmp_cor=calcite_oray_LQB93 --cmp_mnt=limestone_krk --cmp_prt=limestone_dsp --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=30 --wvl_nbr=1000 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc
ncks -O -v wvl,wvn,'idx_rfr_+' -d wvl,2.5e-6, ${DATA}/mie/idx_rfr_mlt.nc ~/foo.nc
ncks -O -v wvl,wvn,'idx_rfr_+' ${DATA}/mie/idx_rfr_mlt.nc ~/foo.nc
ncap -O -s "wvl=wvn" ~/foo.nc ~/foo.nc
ncks -C -v idx_rfr_mtx_rl,idx_rfr_ncl_rl,idx_rfr_ncl_prt,idx_rfr_mtx_img,idx_rfr_ncl_img,idx_rfr_prt_img ~/foo.nc

# Compare Fe2O3_eray_hitran96, Fe2O3_oray_hitran96, Fe2O3_avg_hitran96, Fe2O3_doccd, and Fe2O3_avg_roush to SoT99 Fgr 1: 
# Results suggest that Fe2O3 from Shettle has weak visible and NIR bands
# DOCCD has a strong visible Fe band
# Roush Fe2O3 has nicely banded structure and is overall winner because it's from Earth
mie -3 --dbg=1 --slr_spc=Labs --mss_frc_ncl=0.50 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=Fe2O3_avg_hitran96 --cmp_ncl=Fe2O3_eray_hitran96 --cmp_cor=Fe2O3_oray_hitran96 --cmp_mnt=Fe2O3_doccd --cmp_prt=Fe2O3_avg_roush --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=2.5 --wvl_nbr=230 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc
ncks -O -v wvl,wvn,'idx_rfr_+' ${DATA}/mie/idx_rfr_mlt.nc ~/foo.nc

# Code mie to accept intuitive specification of MCA mixtures, e.g., instead of above, how about
mie -3 --dbg=1 --slr_spc=Labs --ncl=illite,0.3 --ncl=kaolinite,0.3 --ncl=montmorillonite,0.3 --ncl=SiO2_avg_hitran96,0.2 --ncl=Fe2O3_doccd,0.025 --ffc_mdm_typ=ffc_mdm_mxg --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# 3f. Create CAM optical properties from Working Blend #20060621 (See above)
# Manually mix by volume
# 48.5% Quartz + 25% Illite + 25% Montmorillonite + 1% limestone_krk + 0.5% hematite
# This is approximately the same as mass mixing ratios of
# 48% Quartz + 25% Illite + 25% Montmorillonite + 1% limestone_krk + 1% hematite
# and therefore contains 0.7% elemental Fe
# 20060621: This combination asymptotes to has k(0.5 um) ~ 0.0025 and asymptotes to k ~ 0.01 in NIR
# For CAM/SNICAR properties, be sure to specify standard particle density = 2500 kg m-3
mie -3 --dbg=1 --slr_spc=Labs --dns_prt=2500.0 --ncl=0.01 --ncl=0.25 --ncl=0.25 --ncl=0.005 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_doccd --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --wvl_grd=rgl --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=10 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/idx_rfr_mlt.nc

# 3g. DHE02 p. 605 fgr 4 Attempt to capture g for spheroids @ 440 nm
# Particle size distributions from p. 596/597 Tbl 1
# Fine mode:   rds_vma=0.12 um, gsd=exp(0.49+0.1*tau)=1.88, cnc=0.02+0.02*tau(1020nm)
# Coarse mode: rds_vma=1.9 um, gsd=exp(0.63-0.1*tau)=1.88, cnc=0.9*tau(1020nm)
# Figure appears to be for tau(1020nm)=0.7 so coarse mode is
# Fine mode:   rds_vma=0.12 um -> rds_nma=0.0468, gsd=exp(0.49+0.1*tau)=1.75
# Coarse mode: rds_vma=1.90 um -> rds_nma=0.7400, gsd=exp(0.56)=1.75
# Following commands simulate course mode only (DHE02 fgr. 4 has both modes)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=aeronet_CpV --psd_typ=lognormal --wvl_mnm=0.43 --wvl_mxm=0.45 --wvl_grd=rgl --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=10.0 --sz_nbr=10 --rds_vma=1.9 --gsd_anl=1.75 ${DATA}/mie/aer_foo.nc
ncks -C -d wvl,0.44e-6 -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/aer_foo.nc

mie -3 --dbg=1 --vts_flg --asp_rat_hxg=2.0 --slr_spc=Labs --cmp_prt=aeronet_CpV --psd_typ=lognormal --wvl_mnm=0.43 --wvl_mxm=0.45 --wvl_grd=rgl --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=10.0 --sz_nbr=10 --rds_vma=1.9 --gsd_anl=1.75 ${DATA}/mie/aer_foo.nc
ncks -C -d wvl,0.44e-6 -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/aer_foo.nc

# Bi-modal attempt using tau=tau(1020nm)=0.7
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=aeronet_CpV --psd_typ=lognormal --wvl_mnm=0.43 --wvl_mxm=0.45 --wvl_grd=rgl --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=10.0 --sz_nbr=10 --psd=0.0468,1.75,0.02 --psd=0.74,1.75,0.63 ${DATA}/mie/aer_foo.nc
ncks -C -d wvl,0.44e-6 -v ss_alb,asm_prm,idx_rfr_ffc_img ${DATA}/mie/aer_foo.nc

# 4. Equal-V/S treatment of hexagonal prisms
# GrW99 and NGW03 use sigma = 1.23 (see NGW03 p. 5)
# They specify rds_nma (which they call "geometric mean radius" but it is definitely rds_nma) and convert to rds_ffc using rds_ffc=1.11*rds_nma which is appropriate for gsd=1.23 (verified). They use rds_nma when quoting results.

# Reproduce NGW03 p. 6 Fgr. 2
# 20050127: Verified mie produces correct minima at (asp_rat_hxg,nbr_vts_per_hxg)=(0.866,1.65) 
mie -3 --dbg=1 --vts_flg --asp_rat_hxg=0.2 ${DATA}/mie/foo.nc
ncks -m -v asp_rat_hxg,cnc_vts,dmt_hxg,nbr_vts_per_hxg,rds_ctr_vts,rds_hxg,xsa_vts,vlm_vts ${DATA}/mie/foo.nc

# Reproduce NGW03 p. 6 Fgr. 4
mie -3 --dbg=1 --vts_flg --asp_rat_hxg=0.2 --cmp_prt=h2o_ice --sz_grd=log --sz_mnm=0.1 --sz_mxm=500.0 --sz_nbr=100 --rds_nma=1.0 --gsd_anl=1.23 --wvl_grd=log --wvl_mnm=0.2 --wvl_mxm=100.0 --wvl_nbr=100 ${DATA}/aca/aer_h2o_ice_NGW03_fgr4a.nc
ncks -H -C -d wvl,0.5e-6 -v ext_fsh_ffc,ss_co_alb,asm_prm ${DATA}/aca/aer_h2o_ice_NGW03_fgr4a.nc
ncks -H -C -v asp_rat_hxg,cnc_vts,dmt_hxg,nbr_vts_per_hxg,rds_ctr_vts,rds_hxg,xsa_vts,vlm_vts ${DATA}/aca/aer_h2o_ice_NGW03_fgr4a.nc

# 5. Absorbing inclusions in weakly-absorbing spheres (MaS99)

# Reproduce Mar02 p. 771 Fgr. 3a
# 5a. Compute absorption enhancements and cloud MOPs
wvl_nbr='100000'
mie -3 --dbg=1 --abs_ncl_wk_mdm_flg --thr_nbr=1 --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_mdm=air --cmp_ncl=soot --cmp_mtx=h2o_lqd --cmp_prt=h2o_lqd --idx_rfr_mtx='(1.33,1.0e-9)' ${DATA}/mie/mie_Mar02_fgr3a.nc > ${DATA}/mie/foo_Mar02_fgr3a.txt 2>&1 
ncks -P -v abs_fct_MaS99 ${DATA}/mie/mie_Mar02_fgr3a.nc
# 5b. Compute externally mixed soot absorption
mie -3 --dbg=1 --abs_ncl_wk_mdm_flg --thr_nbr=1 --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_mdm=air --cmp_ncl=soot --cmp_mtx=h2o_lqd --cmp_prt=h2o_lqd --idx_rfr_mtx='(1.33,1.0e-9)' ${DATA}/mie/mie_Mar02_fgr3a.nc > ${DATA}/mie/foo_Mar02_fgr3a.txt 2>&1 
ncks -P -v abs_fct_MaS99 ${DATA}/mie/mie_Mar02_fgr3a.nc

tst_nm=mie_mca_tst_5
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
mie -3 --dbg=1 --abs_ncl_wk_mdm_flg --thr_nbr=1 --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_ncl=soot --cmp_mtx=h2o_lqd --idx_rfr_mtx='(1.33,1.0e-9)' ${DATA}/mie/mie_${wvl_nbr}.nc > ${DATA}/mie/foo_${wvl_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,abs_fct_MaS99,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_${wvl_nbr}.nc ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_${wvl_nbr}.nc ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, abs_fct_MaS99, ss_co_alb with ${wvl_nbr} wavelengths...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,abs_fct_MaS99,ss_co_alb ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over wvl_nbr

# 6. Volcanic aerosol for CCM3/CAM3:
# Actual parameters used in CCM3/CAM3 volcanic aerosol are uncertain
# Ammann never sent me the final command he used
# CRB04 states rds_swa=0.426 um and gsd=1.25
# AKZ04 states rds_swa=0.45 um and gsd=1.25
# volcrad.F90 states rds_swa=0.527 um and gsd=exp(0.405)=1.5
# ~/cam_dev/models/atm/cam/src/physics/cam1/volcrad.F90:
#     Command line: ./mie -dbg --aer_typ=PaW75 --dist=lognormal
#                    --wvl_grd=CCM_LW --bnd_nbr=10 --sz_grd=log
#                    --sz_mnm=1.0e-3 --sz_mxm=5.0
#                    --sz_nbr=100 --dst_a=.35 --dst_b=0.405 --dst_c=.1e6
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=PaW75 --psd_typ=lognormal --wvl_grd=CCM_LW --bnd_nbr=10 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=100 --rds_swa=0.426 --gsd_anl=1.25 ${DATA}/aca/aer_volcanic.nc

# 7. Internally mixed soot

# 8. LGGE snow

# Determine radius from wavelength and size parameter
ncap2 -O -v -s 'pi=4.0*atan(1.0);rds_nma=1.742e-6*360.68802/(2*pi);wvl_ctr=1.742e-6;print(rds_nma,"rds_nma = %f\n")' ~/nco/data/in.nc ~/foo.nc

8a. Reproduce snow results of JC Gallet: 
# 1305 nm, 2 layers, 4 streams, top=0.25mm and reff=45 um, second layer depth=1 m with reff=100 µm. T=258°K. density=0.35 kg/m3. Refractive index=(1.2959448,1.3149808e-05)
# Flanner et al. 2007 reports albedo = 0.5711.
# JC finds albedos of 0.5398554 and 0.5404409 for HG and Mie, respectively.
# With reff=30µm instead of reff=45µm in layer one, 
# JC finds albedo=0.5717641 and 0.5723079 with HG and Mie, respectively.

swnb2 --flg_mpr --flx_frc_drc=0.0 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_noatm.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_dff_noatm_045_100.nc > ~/foo_snw 2>&1 &
ncks -d wvl,1.305e-6 -u -C -H -v 'idx_rfr_prt.?' ${DATA}/aca/aer_snw_rds_ffc_045_100um.nc
lev[0]=50 wvl[100]=1.305e-06 idx_rfr_prt_img[570]=1.31498e-05 fraction
lev[0]=0.0125 wvl[100]=1.305e-06 idx_rfr_prt_rl[100]=1.29594 fraction
ncks -d bnd,1.305e-6 -u -C -H -v alb_spc_snw ${DATA}/icr/swnb_prp_dff_noatm_045_100.nc
bnd[566]=1.30463e-06 alb_spc_snw[566]=0.537495 fraction # HG
bnd[566]=1.30463e-06 alb_spc_snw[566]=0.53634 fraction # Mie

8b. Physically plausible snowpack reflectance, e.g., /data/zender/fgr/snicar/snw_sot_spectral_rfl

# Production: 13 min. on tephra, 20 min. on ESMF, 40 min. on virga
mie -3 --mie_slv=Wis79 --flx_frc_drc_sfc=0.0 --cnc_nbr_anl=1.0e12 --wvl_grd=rgl --wvl_nbr=470 --wvl_mnm=0.3 --wvl_mxm=5.0 --bnd_nbr=1 --sz_mnm=0.1 --sz_mxm=2000.0 --sz_grd=log --sz_nbr=200 --rds_swa=100.0 --gsd_anl_dfl=2.3 --cmp_prt=h2o_ice --tpt_prt=258.0 --cmp_cor=lac_HKS98 --hgt_mdp=20.0 --lgn_nbr=16 --ngl_nbr=361 ${DATA}/aca/aer_h2o_ice_snw.nc
ncks -C -F -H -u -d wvl,0.5e-6 -v 'alb_spc_snw.?','.?_dbg','sfc_spc_.?',flx_frc_drc_sfc,tau_ext,ss_alb,wvl ${DATA}/aca/aer_h2o_ice_snw.nc
scp 'tephra.ess.uci.edu:${DATA}/aca/aer_h2o_ice_snw.nc' ${DATA}/aca
scp 'tephra.ess.uci.edu:${DATA}/aca/aer_h2o_ice_snw.nc' ${DATA}/mie/aer_snw_100um.nc
ncks -O -v wvl,wvl_grd,alb_spc_snw ${DATA}/mie/aer_snw_100um.nc ${DATA}/icr/rfl_spc_sfc_snw_100um_dff.nc
ncrename -d wvl,bnd -v wvl,bnd -v alb_spc_snw,rfl_spc_sfc ${DATA}/icr/rfl_spc_sfc_snw_100um_dff.nc
# Compare to Mark's SNICAR version that averages diffuse and direct skies:
/data/zender/icr/rfl_spc_sfc_snc_100um_coszen_0.5.nc

# Single optical calculation for debugging
mie -3 --mie_slv=Wis79 --flx_frc_drc_sfc=0.0 --cnc_nbr_anl=1.0e12 --wvl_grd=rgl --wvl_nbr=1 --wvl_mnm=0.499 --wvl_mxm=0.501 --bnd_nbr=1 --sz_mnm=99.0 --sz_mxm=101.0 --sz_grd=log --sz_nbr=1 --rds_swa=100.0 --gsd_anl_dfl=2.3 --cmp_prt=h2o_ice --tpt_prt=258.0 --cmp_cor=lac_HKS98 --hgt_mdp=0.1 --lgn_nbr=16 --ngl_nbr=361 ${DATA}/mie/mie.nc
ncks -C -F -H -u -d wvl,0.5e-6 -v 'alb_spc_snw.?','.?_dbg','sfc_spc_.?',flx_frc_drc_sfc,tau_ext,ss_alb,wvl ${DATA}/mie/mie.nc

cd ~/anl;ncl 'fl_in_nm="/data/zender/aca/aer_h2o_ice_snw.nc"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl

8c. Arctic clear-sky spectral surface insolation for Florent

   Clear profile, spectrally resolved snow reflectance:
   swnb2 --dbg=1 \
   --fl_clm=${DATA}/aca/prf_sas_65N.nc \
   --fl_out=${DATA}/icr/swnb_sas_65N_clr_rfl.nc \
   --fl_rfl=${DATA}/icr/rfl_spc_sfc_snc_100um_csza0.5.nc \
   --mpc_CWP=0.0 \
   --slr_cst=1367.0 \
   --streams=4 \
   --thermal=false \
   > ~/foo_swnb_clr_rfl.txt 2>&1 &

ncks -C -v '^alb_sfc.?',flx_spc_dwn_sfc,flx_bb_dwn ${DATA}/icr/swnb_sas_65N_clr_rfl.nc | m
ncks -H -C -s '%g, \n' -v wvl_ctr ${DATA}/icr/swnb_sas_65N_clr_rfl.nc > ~/foo

   Clear profile, spectrally resolved snow reflectance:
   swnb2 --dbg=1 \
   --fl_clm=${DATA}/aca/prf_sas_65N.nc \
   --fl_out=${DATA}/icr/swnb_sas_65N_clr_rfl_100um.nc \
   --fl_rfl=${DATA}/icr/rfl_spc_sfc_snw_100um_dff.nc \
   --mpc_CWP=0.0 \
   --slr_cst=1367.0 \
   --streams=4 \
   --thermal=false \
   > ~/foo_swnb_clr_rfl.txt 2>&1 &

ncks -v rfl_bb_sfc,flx_bb_abs_sfc ${DATA}/icr/swnb_sas_65N_clr_rfl_200um.nc 
ncks -v rfl_bb_sfc,flx_bb_abs_sfc ${DATA}/icr/swnb_sas_65N_clr_rfl_100um.nc
ncdiff -O ${DATA}/icr/swnb_sas_65N_clr_rfl_200um.nc ${DATA}/icr/swnb_sas_65N_clr_rfl_100um.nc ${DATA}/icr/swnb_sas_65N_clr_rfl_200um_mns_100um.nc
ncks -v rfl_bb_sfc,flx_bb_abs_sfc ${DATA}/icr/swnb_sas_65N_clr_rfl_200um_mns_100um.nc

ncks -C -H -s '%g, \n' -v bnd rfl_spc_sfc_snw_100um_dff.nc > ~/foo
ncks -C -H -s '%g, \n' -v rfl_spc_sfc rfl_spc_sfc_snw_100um_dff.nc > ~/foo
ncks -C -H -s '%g\n' -v rfl_spc_sfc rfl_spc_sfc_snw_200um_dff.nc > ~/foo
scp ~/domine.txt dust.ess.uci.edu:/data/zender/tmp/domine_100um_200um.txt 
http://dust.ess.uci.edu/tmp/domine_100um_200um.txt

8d. SWNB prognostic snow impurities

Construct snow profile and run model with snow defaults:
clm --lev_snw=5 -D 1 --drc_in ${DATA}/aca -i mls_snw.txt --drc_out ${DATA}/aca -o mls_snw.nc
swnb2 --thermal=false -s 4 -p ${DATA}/aca/mls_snw.nc -d ~/foo.nc > ~/foo_snw 2>&1 &
ncks -F -H -C -v '.snw' ~/foo.nc | m
ncks -F -H -C -d bnd,0.5e-6 -v 'odxc_spc_mpr,odxc_obs_mpr,mmr_mpr_snw,trn_spc_atm_mpr' ~/foo.nc | m
ncks -F -H -C -d bnd,0.5e-6 -v 'flx_bb_abs_snw' ~/foo.nc | m
cd ~/anl;ncl 'fl_in_nm="/home/zender/foo.nc"' 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl

Increase impurity concentration:
for mmr_mpr_ppb in 1 10 100 500 ; do
	mmr_mpr_ppb_fmt=`printf "%03d" $mmr_mpr_ppb`
	swnb2 --mmr_mpr=${mmr_mpr_ppb}e-9 --thermal=false -p ${DATA}/aca/mls_snw.nc -d ${DATA}/icr/swnb_snw_lac_ppb_${mmr_mpr_ppb_fmt}.nc > ~/foo_snw 2>&1 &
done # end loop over mmr_mpr_ppb
ncks -F -H -C -d bnd,0.5e-6 -v 'rfl_bb_snw,flx_bb_abs_snw,odxc_obs_mpr,mmr_mpr_snw' ${DATA}/icr/swnb_snw_lac_ppb_001.nc | mo
cd ~/anl;ncl 'fl_in_nm="/data/zender/icr/swnb_snw_lac_ppb_100.nc"' 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl

# Generate snowpacks of common effective radii 
# Production: 13 min. on tephra, 20 min. on ESMF, 40 min. on virga
for rds_ffc_um in 30 45 100 200 500 ; do
    rds_ffc_um_fmt=`printf "%03d" $rds_ffc_um`
#   gsd='2.3'
    gsd='1.6'
    gsd_sng="gsd${gsd}"
    if [ "${gsd}" = '1.6' ]; then
       sz_mnm=`echo ${rds_ffc_um} | awk '{print 0.2*$1}'`
       sz_mxm=`echo ${rds_ffc_um} | awk '{print 5.0*$1}'`
       sz_nbr='100'
    elif [ "${gsd}" = '2.3' ]; then
       sz_mnm='0.1'
       sz_mxm='2000.0'
       sz_nbr='200'
    else
	printf "Current gsd = ${gsd} is not found"
    fi # endif gsd
    mie -3 --mie_slv=Wis79 --flx_frc_drc_sfc=0.0 --cnc_nbr_anl=1.0e12 --wvl_grd=rgl --wvl_nbr=470 --wvl_mnm=0.3 --wvl_mxm=5.0 --bnd_nbr=1 --sz_mnm=${sz_mnm} --sz_mxm=${sz_mxm} --sz_grd=log --sz_nbr=${sz_nbr} --rds_swa=${rds_ffc_um} --gsd_anl_dfl=${gsd} --cmp_prt=h2o_ice --tpt_prt=258.0 --cmp_cor=lac_HKS98 --hgt_mdp=20.0 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_snw_rds_ffc_${rds_ffc_um_fmt}um_${gsd_sng}.nc > ~/foo_${rds_ffc_um_fmt}um_${gsd_sng}.txt 2>&1 &
done # end loop over rds_ffc_um

# Generate snowpacks of common specific surface area
# Production timings: tephra
#for sfc_spc_cm2xg in 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 ; do
for sfc_spc_cm2xg in 265 340 347 372 773 909 933 977 990 1068 1086 1100 1127 1310 1313 ; do
#   gsd='2.3'
    gsd='1.6'
    gsd_sng="gsd${gsd}"
# Compute effective radius from specific surface area
    rds_ffc_um=`echo ${sfc_spc_cm2xg} | awk '{print 1.0e4*3.0/(0.917*$1)}'`
    if [ "${gsd}" = '1.6' ]; then
       sz_mnm=`echo ${rds_ffc_um} | awk '{print 0.2*$1}'`
       sz_mxm=`echo ${rds_ffc_um} | awk '{print 5.0*$1}'`
       sz_nbr='100'
    elif [ "${gsd}" = '2.3' ]; then
       sz_mnm='0.1'
       sz_mxm='2000.0'
       sz_nbr='200'
    else
	printf "Current gsd = ${gsd} is not found"
    fi # endif gsd
    sfc_spc_cm2xg_fmt=`printf "%04d" $sfc_spc_cm2xg`
    mie -3 --mie_slv=Wis79 --flx_frc_drc_sfc=0.0 --cnc_nbr_anl=1.0e12 --wvl_grd=rgl --wvl_nbr=470 --wvl_mnm=0.3 --wvl_mxm=5.0 --bnd_nbr=1 --sz_mnm=${sz_mnm} --sz_mxm=${sz_mxm} --sz_grd=log --sz_nbr=${sz_nbr} --sfc_spc=${sfc_spc_cm2xg} --gsd_anl_dfl=${gsd} --cmp_prt=h2o_ice --tpt_prt=258.0 --cmp_cor=lac_HKS98 --hgt_mdp=20.0 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_snw_ssa_${sfc_spc_cm2xg_fmt}cm2xg_${gsd_sng}.nc > ~/foo_ssa_${sfc_spc_cm2xg_fmt}_${gsd_sng}.txt 2>&1 &
done # end loop over sfc_spc_cm2xg

# Concatenate layers into vertically resolved snowpack optical properties
ncecat -O -D 1 -u lev -p ${DATA}/aca -o ${DATA}/aca/aer_snw_rds_ffc_045_100um.nc \
	aer_snw_rds_ffc_045um.nc \
	aer_snw_rds_ffc_100um.nc
ncap2 -O -s 'lev[lev]={0.0125,50.0}' ${DATA}/aca/aer_snw_rds_ffc_045_100um.nc ${DATA}/aca/aer_snw_rds_ffc_045_100um.nc

ncecat -O -D 1 -u lev -p ${DATA}/aca -o ${DATA}/aca/aer_snw_rds_ffc_030_100um.nc \
	aer_snw_rds_ffc_030um.nc \
	aer_snw_rds_ffc_100um.nc
ncap2 -O -s 'lev[lev]={0.0125,50.0}' ${DATA}/aca/aer_snw_rds_ffc_030_100um.nc ${DATA}/aca/aer_snw_rds_ffc_030_100um.nc

ncecat -O -D 1 -u lev -p ${DATA}/aca -o ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc \
	aer_snw_rds_ffc_045um.nc \
	aer_snw_rds_ffc_100um.nc \
	aer_snw_rds_ffc_100um.nc \
	aer_snw_rds_ffc_100um.nc \
	aer_snw_rds_ffc_100um.nc
ncap2 -O -s 'lev[lev]={1.0,3.5,7.5,15.0,35.0}' ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc

# Simulations to evaluate H.-W. Jacobi's CROCUS simulations of CEN Col du Porte snowpack measurements
ncecat -O -D 1 -u lev -p ${DATA}/aca -o ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc \
	aer_snw_rds_ffc_086um.nc \
	aer_snw_rds_ffc_085um.nc \
	aer_snw_rds_ffc_085um.nc \
	aer_snw_rds_ffc_086um.nc \
	aer_snw_rds_ffc_087um.nc \
	aer_snw_rds_ffc_088um.nc \
	aer_snw_rds_ffc_162um.nc \
	aer_snw_rds_ffc_091um.nc
ncap2 -O -s 'lev[lev]={0.6,1.7,2.45,3.1,3.8,5.0,6.8,8.2}' ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc ${DATA}/aca/aer_snw_rds_ffc_mlt_lyr.nc

# Test multi-layer swnb2 files
# NB: -I -J -L turns _off_ ice, thermal emission, and liquid, respectively
# mls_snw_FZR07 attempts to re-create conditions of FZR07 Figure 2.
# FZR07 Figure 2 assumes diffuse radiation.
# A cloud produces more realistic direct/diffuse partitioning.
# However, clouds are less reproducible than 100% direct/diffuse switches
# for intermodel comparisons, because swnb2 is one of (very few) models
# that simultaneously compute atmospheric and snow RT.
# Re-run these with cloud-producing diffuse illuimination for more realistic results.
swnb2 --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p snw_noatm.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_dff_noatm_045_100.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_sal_dff_045_100.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_dff_045_100.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mie --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_mie_045_100.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_030_100um.nc -d ${DATA}/icr/swnb_prp_dff_030_100.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ${DATA}/icr/swnb_prp_dff_100.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mpr --flx_frc_drc=0.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_dff.nc > ~/foo_snw 2>&1 &
swnb2 --flg_mpr --zen --flx_frc_drc=1.0 -I -J -L -s 4 --drc_in=${DATA}/aca -p mls_snw_FZR07.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_drc_045_100.nc > ~/foo_snw 2>&1 &

ncks -C -F -H -u -d bnd,0.635e-6 -d bnd,1.310e-6 -d bnd,1.550e-6 -v alb_spc_snw ${DATA}/icr/swnb_prp_dff_noatm_045_100.nc

# Lab dirty snow experiments
# New Control: work-in-progress
swnb2 --flx_frc_drc=0.0 --flg_mpr --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ${DATA}/icr/swnb_prp_lgge_100.nc > ~/foo_snw 2>&1 &
# Control
swnb2 --flx_frc_drc=0.0 --flg_mpr --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_rds_ffc_045_100um.nc -d ${DATA}/icr/swnb_prp_lgge_045_100.nc > ~/foo_snw 2>&1 &

# Select photodiode wavelengths
ncks -C -F -H -u -d bnd,0.635e-6 -d bnd,1.310e-6 -d bnd,1.550e-6 -d bnd,1.610e-6 -d bnd,1.850e-6 -v alb_spc_snw ${DATA}/icr/swnb_prp_lgge_100um.nc

# 8e. Correct LGGE lab measurements for horizontal effects
# Test RTA correction
swnb2 --mmr_mpr=0.0e-6 --flx_frc_drc=0.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_25mm_32lvl.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ${DATA}/icr/swnb_lgge_32lvl.nc > ~/foo_snw 2>&1
ncks -F -u -C -H -d bnd,0.5e-6 -v '._snw_cnt','alb_spc_snw.?','.?_dea' ${DATA}/icr/swnb_lgge_32lvl.nc | more
swnb2 -D 1 -E -e 1603 --mmr_mpr=0.0e-6 --flx_frc_drc=0.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_25mm_32lvl.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ${DATA}/icr/swnb_lgge_32lvl.nc > ~/foo_snw 2>&1

# Increase echantillon depth:
for dpt_snw in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 90 100 125 150 175 200 225 250 ; do
	dpt_snw_fmt=`printf "%03d" $dpt_snw`
	fl_out=${DATA}/icr/swnb_lgge_100um_${dpt_snw_fmt}dmm.nc
	swnb2 --dpt_snw=${dpt_snw}e-4 --mmr_mpr=0.0e-6 --flx_frc_drc=0.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_rds_ffc_100um.nc -d ${fl_out} > ~/foo_snw 2>&1
#	fl_out=${DATA}/icr/swnb_lgge_dns_040_ssa_1000_${dpt_snw_fmt}dmm.nc
#	swnb2 --dns_snw=40.0 --dpt_snw=${dpt_snw}e-4 --mmr_mpr=0.0e-6 --flx_frc_drc=0.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_1000cm2xg_gsd2.3.nc -d ${fl_out} > ~/foo_snw 2>&1
#	ncks -F -H -C -d bnd,0.5e-6 -v 'rfl_bb_snw,alb_spc_snw' ${fl_out}
#	ncks -F -H -C -d bnd,0.635e-6 -v 'rfl_bb_snw,alb_spc_snw,dpt_dlt_snw,mmr_mpr_snw' ${fl_out}
#	ncks -Q -F -H -C -d bnd,0.635e-6 -v 'alb_spc_snw' ${fl_out}
	echo ${fl_out}
done # end loop over dpt_snw

# NB: Much of the rest of section 8e is obsoleted by theory embedded in snw.cc
# Geometric correction to plane-parallel assumption for LGGE integrating sphere
fl_rx='100um' # [sng] Files matching this pattern are included
#fl_rx='dns_040_ssa_1000' # [sng] Files matching this pattern are included
ncecat -O -M -3 -t 1 -u dpt_snw -d bnd,0.635e-6 -d bnd,1.310e-6 -v dpt_dlt_snw,rfl_bb_snw,alb_spc_snw ${DATA}/icr/swnb_lgge_${fl_rx}_*dmm.nc ~/swnb_lgge_${fl_rx}_nsm.nc
ncwa -O -a lev_snw ~/swnb_lgge_${fl_rx}_nsm.nc ~/swnb_lgge_${fl_rx}_nsm.nc
ncrename -O -v dpt_dlt_snw,dpt_snw ~/swnb_lgge_${fl_rx}_nsm.nc
ncpdq -O -a -bnd ~/swnb_lgge_${fl_rx}_nsm.nc ~/swnb_lgge_${fl_rx}_nsm.nc
cat > ~/swnb_lgge.nco << EOF
/* Purpose: Correct SWNB2-predicted reflectances for geometry of LGGE integrating sphere
   NB: Script is mouse-able because of presence of \$0 in below
   Backslash prevents shell from expanding $0 into '/bin/bash'
   Remove backslash to store script separately */

// Incremental (ncr) layer contributions to albedo
lev_nbr=dpt_snw.size();
rfl_bb_snw_ncr=rfl_bb_snw; // [frc] Incremental (ncr) layer contributions to albedo
alb_spc_snw_ncr=alb_spc_snw; // [frc] Incremental (ncr) layer contributions to albedo
rfl_bb_snw_ncr(0)=rfl_bb_snw(0);
alb_spc_snw_ncr(0,:)=alb_spc_snw(0,:);
for(lev_idx=1;lev_idx<lev_nbr;lev_idx++){ // NB: loop starts with 1
  rfl_bb_snw_ncr(lev_idx)=rfl_bb_snw(lev_idx)-rfl_bb_snw(lev_idx-1);
  alb_spc_snw_ncr(lev_idx,:)=alb_spc_snw(lev_idx,:)-alb_spc_snw(lev_idx-1,:);
} // end loop over lev
rfl_bb_snw_ncr_ttl=rfl_bb_snw_ncr.total();
alb_spc_snw_ncr_ttl=alb_spc_snw_ncr.total(\$0);
rfl_bb_snw_frc=rfl_bb_snw; // [frc] Fractional contributions to albedo
alb_spc_snw_frc=alb_spc_snw; // [frc] Fractional contributions to albedo
for(lev_idx=0;lev_idx<lev_nbr;lev_idx++){
  rfl_bb_snw_frc(lev_idx)=rfl_bb_snw_ncr(lev_idx)/rfl_bb_snw_ncr_ttl;
  alb_spc_snw_frc(lev_idx,:)=alb_spc_snw_ncr(lev_idx,:)/alb_spc_snw_ncr_ttl;
} // end loop over lev

// Layer geometrical (gmt) correction factor (fct)
// Echantillon geometry from ~/anl/snw_rfl.ncl
dpt_smp_hld_13mm=0.01327; // [m] Depth of 13mm-deep "small" sample holder
dpt_smp_hld_25mm=0.0250; // [m] Depth of 25mm-deep "large" sample holder
dmt_smp_hld=0.0633; // [m] Diameter of sample holders
rds_smp_hld=dmt_smp_hld/2.0; // [m] Radius of sample holder
ngl_psi=atan(dpt_snw/rds_smp_hld); // [rdn] Angle subtended by wall
pi=4.0*atan(1.0); // [frc] 3
fct_gmt=(pi-2*ngl_psi)/pi; // [frc] Geometry factor
fct_gmt(0)=1.0; // [frc] Geometry factor

// Incremental (ncr) layer corrections (crc) to surface albedo
rfl_bb_snw_crc=rfl_bb_snw;
alb_spc_snw_crc=alb_spc_snw;
rfl_bb_snw_ncr_crc=rfl_bb_snw_ncr;
alb_spc_snw_ncr_crc=alb_spc_snw_ncr;
for(lev_idx=0;lev_idx<lev_nbr;lev_idx++){
  rfl_bb_snw_ncr_crc(lev_idx)=rfl_bb_snw_ncr(lev_idx)*fct_gmt(lev_idx);
  alb_spc_snw_ncr_crc(lev_idx,:)=alb_spc_snw_ncr(lev_idx,:)*fct_gmt(lev_idx);
} // end loop over lev

// Apply corrections to obtain corrected (crc) albedos
rfl_bb_snw_crc(0)=rfl_bb_snw_ncr(0);
alb_spc_snw_crc(0,:)=alb_spc_snw_ncr(0,:);
for(lev_idx=1;lev_idx<lev_nbr;lev_idx++){
  rfl_bb_snw_crc(lev_idx)=rfl_bb_snw_crc(lev_idx-1)+rfl_bb_snw_ncr(lev_idx)*fct_gmt(lev_idx);
  alb_spc_snw_crc(lev_idx,:)=alb_spc_snw_crc(lev_idx-1,:)+alb_spc_snw_ncr(lev_idx,:)*fct_gmt(lev_idx);
} // end loop over lev
rfl_bb_snw_ncr_crc_ttl=rfl_bb_snw_ncr_crc.total();
alb_spc_snw_ncr_crc_ttl=alb_spc_snw_ncr_crc.total(\$0);

// Sensitivity to horizontal boundary conditions
rfl_bb_snw_crc_bc=rfl_bb_snw_ncr_crc_ttl-rfl_bb_snw_ncr_ttl;
alb_spc_snw_crc_bc=alb_spc_snw_ncr_crc_ttl-alb_spc_snw_ncr_ttl;
EOF
ncap2 -O -S ~/swnb_lgge.nco ~/swnb_lgge_${fl_rx}_nsm.nc ~/swnb_lgge_${fl_rx}_nsm.nc

# 8f. BC MAC uncertainty
# FZR07 base low, central, and high estimates (LE, CE, HE) of BC MAC on BoB06
# LE and HE are CE -/+ one standard deviation of observations
# Hydrophobic: MAC =  7500 +/- 1200 m2 kg-1 at 550 nm (tuned with BC density)
# Hydrophilic: MAC = 11000 +/- 1200 m2 kg-1 at 550 nm (tuned with BC density)
ncks -C -d wvl,0.55e-6 -v ss_alb,asm_prm,abs_cff_mss ${DATA}/aca/aer_lac_phb_ctr_FZR07.nc

# 8g. Sample SSA from NIR reflectance
# JC finds SSA [cm2 g-1] = 28.772 exp (5.79 rfl ) where rfl = 1310nm reflectance
# Mean SSA and density of clean snow samples:
# 20071207e 265 cm2 g-1 
# 20071210a 340 cm2 g-1
# 20071211d 372 cm2 g-1
# 20080130e 347 cm2 g-1

# Sale (dirty) experiments for nominal snowpack
swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Sale (dirty) experiments with measured density
swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc -d ${DATA}/icr/swnb_lgge_100um_dns_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Sale (dirty) experiments for estimated snowpack SSA
swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0372cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0340cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0265cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Sale (dirty) experiments for estimated snowpack SSA and density
swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0372cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0340cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_0265cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Sale (dirty) experiments for estimated density and various optical properties
for lac_typ in phb_low phb_ctr phb_hgh phl_low phl_ctr phl_hgh ; do
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppb_000.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppb_250.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppb_750.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppm_001.5.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppm_007.5.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppm_015.9.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc -d ${DATA}/icr/swnb_lgge_100um_${lac_typ}_ppm_221.0.nc > ~/foo_snw 2>&1 &
#	ncks -Q -F -H -C -d bnd,0.635e-6 -v 'alb_spc_snw' ${fl_out}
done # end loop over lac_typ 

# Sale (dirty) experiments for estimated SSA, density, and various optical properties
for lac_typ in phb_low phb_ctr phb_hgh phl_low phl_ctr phl_hgh ; do
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppb_000.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppb_250.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppb_750.0.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0372cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppm_001.5.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppm_007.5.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0340cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppm_015.9.nc > ~/foo_snw 2>&1 &
	swnb2 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_${lac_typ}_FZR07.nc --fl_snw=aer_snw_ssa_0265cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_ssa_dns_${lac_typ}_ppm_221.0.nc > ~/foo_snw 2>&1 &
#	ncks -Q -F -H -C -d bnd,0.635e-6 -v 'alb_spc_snw' ${fl_out}
done # end loop over lac_typ 

# Sale (dirty) test (tst) experiments: SSA, density, phb_ctr
# Protocol: Use "tst" experiments for cutting-edge tests
# When test results look good, change output files from "tst" to "bee"
# and continue to use "tst" for more development.
# "bee" files are used to generate graphics for manuscripts, seminars
# "tst" files are used for sensitivity tests
# Normally, "tst" is the same as "bee" except for one or two changes
# Current (20080928) differences between "tst" and "bee": None
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0372cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0340cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0265cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_tst_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Sale (dirty) experiments for best estimate of everything (bee): SSA, density, phb_he
# NB: EGU work was done with direct flux, cos(theta)=0.502701676369
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --dns_snw=219.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppb_000.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=250.00e-9 --dns_snw=212.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppb_250.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=750.00e-9 --dns_snw=226.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppb_750.0.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=001.45e-6 --dns_snw=228.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0372cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppm_001.5.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=007.50e-6 --dns_snw=227.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0347cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppm_007.5.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=015.90e-6 --dns_snw=243.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0340cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppm_015.9.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=221.00e-6 --dns_snw=254.0 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_mpr=aer_lac_phb_ctr_FZR07.nc --fl_snw=aer_snw_ssa_0265cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_bee_ppm_221.0.nc > ~/foo_snw 2>&1 &

# Aggregate ensembles by impurity concentration
for nsm_sng in tst ssa; do
#for nsm_sng in ssa tst 100um bee ssa_dns 100um_dns 100um_phb_low 100um_phb_ctr 100um_phb_hgh 100um_phl_low 100um_phl_ctr 100um_phl_hgh ssa_dns_phb_low ssa_dns_phb_ctr ssa_dns_phb_hgh ssa_dns_phl_low ssa_dns_phl_ctr ssa_dns_phl_hgh; do
    fl_out="${DATA}/icr/swnb_lgge_${nsm_sng}_nsm_mpr.nc"
    ncecat -O -M -3 -t 1 -u mmr_mpr_snw -d bnd,0.635e-6 -d bnd,1.310e-6 -v mmr_mpr_snw,rfl_bb_snw,alb_spc_snw ${DATA}/icr/swnb_lgge_${nsm_sng}_pp?_???.?.nc ${fl_out}
    ncwa -O -a lev_snw ${fl_out} ${fl_out}
    ncpdq -O -a -bnd ${fl_out} ${fl_out}
    echo "Created ensemble file ${fl_out}"
done # end loop over nsm_sng

# Range due to optical properties
for mpr_sng in ppb_000 ppb_250 ppb_750 ppm_001.45 ppm_007.50 ppm_015.9 ppm_221 ; do
    fl_out="${DATA}/icr/swnb_lgge_100um_phb_lowm_phl_hgh_${mpr_sng}.nc"
    ncbo -t 1 -O -p ${DATA}/icr swnb_lgge_100um_phb_low_${mpr_sng}.nc swnb_lgge_100um_phl_hgh_${mpr_sng}.nc ${fl_out}
    echo ${fl_out}
    ncks -F -H -C -d bnd,0.635e-6 -d bnd,1.310e-6 -v 'alb_spc_snw' ${fl_out}
    for lac_typ in phb phl ; do
        ncbo -O -p ${DATA}/icr swnb_lgge_100um_${lac_typ}_low_${mpr_sng}.nc swnb_lgge_100um_${lac_typ}_hgh_${mpr_sng}.nc ${DATA}/icr/swnb_lgge_100um_${lac_typ}_lowmhgh_${mpr_sng}.nc 
    done # end loop over lac_typ 
done # end loop over mpr_sng

# 8h. GDZ08 outliers
# Table 1 row 1:
# ssa 106.8 m2 kg-1, dns 41.0 kg m-3, rfl(obs)=53.43, rfl(mdl)=61.51 (JCG) 62.50 (CSZ 69% drc) 58.12 (CSZ 69% drc w/crc)
swnb2 --zen --flx_frc_drc=0.71 --mmr_mpr=000.00e-9 --dns_snw=041.0 --alb_sfc_vsb=0.06 --alb_sfc_NIR=0.03 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge.nc --fl_snw=aer_snw_ssa_1000cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_GDZ08_tbl1a_01lvl.nc > ~/foo_snw 2>&1 &
# Brightest, lightest snow (correction is 0.625-0.581=0.044)
swnb2 --zen --flx_frc_drc=0.71 --mmr_mpr=000.00e-9 --dns_snw=041.0 --alb_sfc_vsb=0.06 --alb_sfc_NIR=0.03 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_25mm_32lvl.nc --fl_snw=aer_snw_ssa_1000cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_GDZ08_tbl1a_25mm_32lvl.nc > ~/foo_snw 2>&1 &
# More typical, but still bright, snow (correction is 0.547-0.541=0.006)
swnb2 --zen --flx_frc_drc=0.71 --mmr_mpr=000.00e-9 --dns_snw=400.0 --alb_sfc_vsb=0.06 --alb_sfc_NIR=0.03 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_25mm_32lvl.nc --fl_snw=aer_snw_ssa_0600cm2xg_gsd2.3.nc -d ${DATA}/icr/swnb_lgge_GDZ08_tbl1d_25mm_32lvl.nc > ~/foo_snw 2>&1 &
ncks -H -C -d bnd,1.310e-6 -d bnd,1.550e-6 -v alb_spc_snw,'.?_dea' ${DATA}/icr/swnb_lgge_GDZ08_tbl1a_01lvl.nc
ncks -H -C -d bnd,1.310e-6 -d bnd,1.550e-6 -v rfl_dff_spc_snw_cnt,alb_spc_snw,'.?_dea' ${DATA}/icr/swnb_lgge_GDZ08_tbl1a_25mm_32lvl.nc | m

# Normalize layer contributions to albedo by "real" DISORT albedo
ncap2 -3 -O -s 'rfl_dff_spc_snw_cnt_nrm=rfl_dff_spc_snw_cnt;where(alb_dff_spc_snw_dea > 0.0) rfl_dff_spc_snw_cnt_nrm=rfl_dff_spc_snw_cnt*alb_spc_snw/alb_dff_spc_snw_dea;rfl_dff_spc_snw_cnt_nrm_ttl=rfl_dff_spc_snw_cnt_nrm.total($levp_snw)' ${DATA}/icr/swnb_lgge_GDZ08_tbl1a_25mm_32lvl.nc ~/foo.nc
ncks -H -C -d bnd,1.310e-6 -v 'rfl_dff_spc_snw_cnt.?',alb_spc_snw,'.?_dea' ~/foo.nc | m
# Pre-computed FOVs with, e.g., 
# snw --llm=dff --dpt_typ=GDZ08 --azi_nbr=100 --rds_nbr=100
# Add pre-computed FOVs
ncap2 -O -s 'ngl_sld_avg_hms[levp_snw]={0.99973,0.999191,0.998652,0.998113,0.997574,0.996631,0.995284,0.993936,0.992588,0.990567,0.986525,0.981137,0.975751,0.970367,0.964985,0.959607,0.951546,0.940814,0.930103,0.914085,0.890191,0.850821,0.799476,0.749783,0.702048,0.656506,0.613322,0.572594,0.534361,0.498614,0.4653,0.434338,0.434338};' ~/foo.nc ~/foo.nc # 13mm sample holder:
ncap2 -O -s 'ngl_sld_avg_hms[levp_snw]={0.99973,0.999191,0.998652,0.998113,0.997574,0.996631,0.995284,0.993936,0.992588,0.990567,0.986525,0.978444,0.959607,0.932779,0.906101,0.879631,0.853425,0.827534,0.802007,0.77689,0.752224,0.716148,0.669927,0.626022,0.584551,0.545569,0.483326,0.407008,0.344012,0.292336,0.250004,0.215256,0.215256};' ~/foo.nc ~/foo.nc # 25mm sample holder:
# Correct predicted-surface albedo for geometric artifact
ncap2 -O -s 'alb_spc_snw_crc=(ngl_sld_avg_hms*rfl_dff_spc_snw_cnt_nrm).total($levp_snw)' ~/foo.nc ~/foo.nc
ncks -H -C -d bnd,1.310e-6 -v ngl_sld_avg_hms,'rfl_dff_spc_snw_cnt.?','alb_spc_snw.?',alb_dff_spc_snw_dea ~/foo.nc | m

# 8i. Inhomogeneous (Density gradient) snow effect for light snow
clm -D 1 --prs_top=101200 --lev_atm=1 --lev_snw=2 --drc_in ${DATA}/aca -i ~/aca/snw_lgge_hmg.txt --drc_out ${DATA}/aca -o snw_lgge_hmg.nc
clm -D 1 --prs_top=101200 --lev_atm=1 --lev_snw=2 --drc_in ${DATA}/aca -i ~/aca/snw_lgge_nhm.txt --drc_out ${DATA}/aca -o snw_lgge_nhm.nc
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_hmg.nc --fl_snw=aer_snw_ssa_0300cm2xg_gsd1.6.nc -d ${DATA}/icr/swnb_lgge_hmg.nc > ~/foo_snw 2>&1 &
swnb2 -z 0.502701676369 --flx_frc_drc=1.0 --mmr_mpr=000.00e-9 --alb_sfc=0.06 --thermal=false -s 4 --drc_in=${DATA}/aca -p snw_lgge_nhm.nc --fl_snw=aer_snw_ssa_0300cm2xg_gsd1.6.nc -d ${DATA}/icr/swnb_lgge_nhm.nc > ~/foo_snw 2>&1 &
ncks -H -C -d bnd,1.310e-6 -v alb_spc_snw ${DATA}/icr/swnb_lgge_hmg.nc
ncks -H -C -d bnd,1.310e-6 -v alb_spc_snw ${DATA}/icr/swnb_lgge_nhm.nc

# 9. Phase functions
# Best tabular expansions for phase functions are Wis771 p. 1421 Table A2
# p. 1410, just after eqn (12) defines size distribution
# Gamma distribution, effective radius = 10 um, effective variance = 0.2
# N = lgn_nbr is order of Legendre expansion 
#   = # of terms in 
# L = 2*ngl_nbr-1 = Quadrature order for Legendre expansion coefficients
#   = # of polar angles where Mie calculations performed for accurate expansion

# 0.5 um: N=902
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=100.0 --sz_nbr=300 --sz_dbg=10.0 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=0.499 --wvl_mxm=0.501 --wvl_nbr=1 --bnd_nbr=1 --lgn_nbr=902 --ngl_nbr=101 --idx_rfr_prt='1.335' ${DATA}/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc
ncks -H -d lgn,,4 -d ngl,0.0 -d wvl,0.5e-7 -C -v '^[[:alpha:]]{3}_rsl_frc$',asm_prm,ss_alb,idx_rfr_prt_'.?',lgn_xpn_cff,phz_fnc_'.?' -m /data/zender/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc | m

# NB: experimental gamma distribution---does not work yet
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=100.0 --sz_nbr=300 --psd_typ=gamma --rds_ffc_gmm=10.0 --var_ffc_gmm=0.2 --wvl_mnm=0.499 --wvl_mxm=0.501 --wvl_nbr=1 --bnd_nbr=1 --lgn_nbr=200 --ngl_nbr=361 --idx_rfr_prt='1.335' ${DATA}/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc

# 1.0 um: N=453
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=9.99 --sz_mxm=10.01 --sz_nbr=1 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=0.999 --wvl_mxm=1.001 --wvl_nbr=1 --bnd_nbr=1 --lgn_nbr=16 --ngl_nbr=361 --idx_rfr_prt='1.328+3.35e-6i' ${DATA}/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc
ncks -H -d lgn,0,1 -d ngl,0.0 -d wvl,1.0e-6 -C -v '^[[:alpha:]]{3}_rsl_frc$',asm_prm,ss_alb,idx_rfr_prt_'.?',lgn_xpn_cff,phz_fnc_'.?' -m /data/zender/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc | m

# 2.0 um: N=230
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=9.99 --sz_mxm=10.01 --sz_nbr=1 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=1.999 --wvl_mxm=2.001 --wvl_nbr=1 --bnd_nbr=1 --lgn_nbr=16 --ngl_nbr=361 --idx_rfr_prt='1.306+1.16e-3i' ${DATA}/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc
ncks -H -d lgn,0,1 -d ngl,0.0 -d wvl,2.0e-6 -C -v '^[[:alpha:]]{3}_rsl_frc$',asm_prm,ss_alb,idx_rfr_prt_'.?',lgn_xpn_cff,phz_fnc_'.?' -m /data/zender/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc | m

# 10.0 um: N=50
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=9.99 --sz_mxm=10.01 --sz_nbr=1 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=9.999 --wvl_mxm=10.001 --wvl_nbr=1 --bnd_nbr=1 --lgn_nbr=16 --ngl_nbr=361 --idx_rfr_prt='1.220+0.0515i' ${DATA}/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc
ncks -H -d lgn,0,1 -d ngl,0.0 -d wvl,10.0e-6 -C -v '^[[:alpha:]]{3}_rsl_frc$',asm_prm,ss_alb,idx_rfr_prt_'.?',lgn_xpn_cff,phz_fnc_'.?' -m /data/zender/mie/aer_h2o_lqd_rds_swa_10_dgn_phz.nc | m
