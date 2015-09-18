# $Id$

# Purpose: Generate standard aerosol optical properties for ${DATA}/aca
# These "aer_*.nc" files are used in ARESE (ZBP97, Zen99) 
# and in other studies (ZeT05)

# Usage: ${HOME}/mie/mie.sh >${HOME}/mie/foo.mie 2>&1 &

# Fe2O3 Hematite (same size distribution as saharan_dust)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=Fe2O3 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_Fe2O3.nc
ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_Fe2O3.nc

# SiO2 Quartz (same size distribution as saharan_dust)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=SiO2 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_SiO2.nc
ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_SiO2.nc

# NB: BoB06 quotes central estimate soot MAC at 0.55 um (not 0.50 um) = 7500 m2 g-1
# ChC90 soot 
# ChC90 soot is my "best guess" soot property and is evolving
# ChC90 soot uses rds_vma=0.1, dns_prt=1800 kg m-3
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --dns_prt=1800.0 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_vma=0.1 --gsd_anl=2.0 ${DATA}/aca/aer_lac_ChC90.nc
ncks -d wvl,0.55e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_lac_ChC90.nc

# FZR07 hydrophobic (uncoated) soot: rds_nma=0.05, gsd = 1.5
# Create three aerosols with low, central, high estimates (le, ce, he)
# Tune densities so MAC(550 nm)=7500 +/- 1200 m2 kg-1 with ChC90
# phb_le_FZR07: dns_prt=1533 kg m-3 (=1322*8.7/7.5)
# phb_ce_FZR07: dns_prt=1322 kg m-3
# phb_he_FZR07: dns_prt=1110 kg m-3 (=1322*6.3/7.5)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --dns_prt=1322.0 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.05 --gsd_anl=1.5 ${DATA}/aca/aer_lac_FZR07.nc
/bin/ln -s ${DATA}/aca/aer_lac_FZR07.nc ${DATA}/aca/aer_lac_phb_ce_FZR07.nc
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --dns_prt=1533.0 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.05 --gsd_anl=1.5 ${DATA}/aca/aer_lac_phb_le_FZR07.nc
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --dns_prt=1110.0 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.05 --gsd_anl=1.5 ${DATA}/aca/aer_lac_phb_he_FZR07.nc
ncks -d wvl,0.55e-6 -v abs_cff_mss,ext_cff_mss,ss_alb,asm_prm ${DATA}/aca/aer_lac_phb_ce_FZR07.nc
# FZR07 hydrophilic (coated) soot: rds_nma=0.0835 um (=0.05*1.67), gsd = 1.5, rds_frc_cor=0.6, mnt=1.67*cor
# phl_le_FZR07: dns_cor=1533 kg m-3 (=1322*8.7/7.5)
# phl_ce_FZR07: dns_cor=1322 kg m-3
# phl_he_FZR07: dns_cor=1110 kg m-3 (=1322*6.3/7.5)
mie -3 --coat_nrm_by_core_flg --coat_flg --rds_frc_cor=0.6 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --cmp_cor=lac_ChC90 --dns_cor=1322.0 --cmp_mnt=sulfate --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.5 --sz_nbr=400 --rds_nma=0.0835 --gsd_anl=1.5 ${DATA}/aca/aer_lac_phl_ce_FZR07.nc
ncks -H -d wvl,0.55e-6 -v abs_cff_mss,ext_cff_mss,ss_alb,asm_prm ${DATA}/aca/aer_lac_phl_ce_FZR07.nc
mie -3 --coat_nrm_by_core_flg --coat_flg --rds_frc_cor=0.6 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --cmp_cor=lac_ChC90 --dns_cor=1533.0 --cmp_mnt=sulfate --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.5 --sz_nbr=400 --rds_nma=0.0835 --gsd_anl=1.5 ${DATA}/aca/aer_lac_phl_le_FZR07.nc
mie -3 --coat_nrm_by_core_flg --coat_flg --rds_frc_cor=0.6 --dbg=1 --slr_spc=Labs --cmp_prt=lac_ChC90 --cmp_cor=lac_ChC90 --dns_cor=1110.0 --cmp_mnt=sulfate --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.5 --sz_nbr=400 --rds_nma=0.0835 --gsd_anl=1.5 ${DATA}/aca/aer_lac_phl_he_FZR07.nc

# CAM soot is called HKS98 soot because it is from OPAC
# HKS98 (and DKS91) use rds_nma=0.0118, gsd=2.0
# HKS98 use dns_prt=1000 kg m-3
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_HKS98 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.0118 --gsd_anl=2.0 ${DATA}/aca/aer_lac_HKS98.nc
ncks -d wvl,0.55e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_lac_HKS98.nc

# WaW80 soot 
# WaW80 use rds_vma=0.1 um (based on HaS95?) and dns_prt=2050 kg m-3
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=lac_WaW80 --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_vma=0.1 --gsd_anl=2.0 ${DATA}/aca/aer_lac_WaW80.nc
ncks -d wvl,0.55e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_lac_WaW80.nc

# PaG77 West Texas
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=dust_like --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_dust_like.nc
ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_dust_like.nc

# PaG77 West Texas
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=mineral_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_mineral_dust.nc
ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_mineral_dust.nc

# PaG77 West Texas
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=saharan_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/aca/aer_saharan_dust.nc
ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_saharan_dust.nc

# Zender's ARESE mineral_dust (36 N) (r=0.2, sigma=0.2) (see table below)
#.500/.862 micron extinction ratio = 1.8-2.0,
#.413/.500 micron extinction ratio = 1.4-1.5, 
#.413/.862 micron extinction ratio = 3.0, 
# based on 951011 after accounting for the extinction of the stratospheric aerosol)
#mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=mineral_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.2 --gsd_anl=1.22 ${DATA}/aca/aer_mineral_dust.nc

# Zender's ARESE dust_like (36 N) (r=0.2, sigma=0.2)
#mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=dust_like --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.2 --gsd_anl=1.22 ${DATA}/aca/aer_dust_like.nc

# Hofmann's pre-El Chichon (27 N)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=sulfate --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.0 --sz_nbr=400 --rds_nma=0.08 --gsd_anl=1.6 ${DATA}/aca/aer_sulfate.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_sulfate.nc

# Zender's ARESE SAGE II mode (36 N) (r, sigma chosen so 0.525/1.02 micron extinction ratio = 2.43)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=h2so4_215K --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.24 --gsd_anl=1.4 --cnc_nbr_anl=0.1e6 ${DATA}/aca/aer_h2so4_215K.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_h2so4_215K.nc

# Zender's ARESE SAGE II mode (36 N) (r, sigma chosen so 0.525/1.02 micron extinction ratio = 2.43)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=h2so4_300K --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.24 --gsd_anl=1.4 --cnc_nbr_anl=0.1e6 ${DATA}/aca/aer_h2so4_300K.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_h2so4_300K.nc

# Hofmann's pre-El Chichon (27 N)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=volcanic_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.08 --gsd_anl=1.6 ${DATA}/aca/aer_volcanic_dust.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_volcanic_dust.nc

# Hofmann's pre-El Chichon (27 N)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=meteoric_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.08 --gsd_anl=1.6 ${DATA}/aca/aer_meteoric_dust.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/aca/aer_meteoric_dust.nc

# Deshler's mode 1 (41 N)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=h2so4_215K --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.02 --gsd_anl=1.6 --cnc_nbr_anl=70.e6 ${DATA}/mie/foo1.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/mie/foo1.nc
#ncks -C -H -v ext_cff_mss -d wvl,0.525e-6 ${DATA}/mie/foo1.nc
#ncks -C -H -v ext_cff_mss -d wvl,1.020e-6 ${DATA}/mie/foo1.nc

# Deshler's mode 2 (41 N)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=h2so4_215K --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=1.4 --cnc_nbr_anl=0.1e6 ${DATA}/mie/foo2.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/mie/foo2.nc
#ncks -C -H -v ext_cff_mss -d wvl,0.525e-6 ${DATA}/mie/foo2.nc
#ncks -C -H -v ext_cff_mss -d wvl,1.020e-6 ${DATA}/mie/foo2.nc

# Zender's ARESE SAGE II mode (36 N) (r, sigma chosen so .525/1.02 micron extinction ratio = 2.43)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=h2so4_215K --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=200 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=.24 --gsd_anl=1.4 --cnc_nbr_anl=.1e6 ${DATA}/mie/foo3.nc
#ncks -d wvl,0.5e-6 -v ext_cff_mss,ext_cff_vlm,ss_alb,asm_prm ${DATA}/mie/foo3.nc
#ncks -C -H -v ext_cff_mss -d wvl,0.525e-6 ${DATA}/mie/foo3.nc
#ncks -C -H -v ext_cff_mss -d wvl,1.020e-6 ${DATA}/mie/foo3.nc

# Zender's ARESE mineral_dust (36 N) (r, sigma chosen so 
# 0.500/0.862 micron extinction ratio = 1.8-2.0,
# 0.413/0.500 micron extinction ratio = 1.4-1.5, 
# 0.413/0.862 micron extinction ratio = 3.0, 
#based on 951011 after accounting for the extinction of the stratospheric aerosol)
mie -3 --dbg=1 --slr_spc=Labs --cmp_prt=mineral_dust --psd_typ=lognormal --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=47 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=2.2 ${DATA}/mie/foo4.nc
ncks -d wvl,0.5e-6 -v ss_alb,asm_prm ${DATA}/mie/foo4.nc
ncks -C -v ext_cff_mss -d wvl,0.413e-6 ${DATA}/mie/foo4.nc
ncks -C -v ext_cff_mss -d wvl,0.5e-6 ${DATA}/mie/foo4.nc
ncks -C -v ext_cff_mss -d wvl,0.862e-6 ${DATA}/mie/foo4.nc

#Mineral dust parameters and ext_cff_mss
#a	b	413nm	500nm	862nm	.5 g	.5 omega
#.025	.9	2760	2536	1329	

#.4	1.	308	314	330
#.4	.7884	443	454	485	.75	.85
#.4	.3	1595	1827	2126

#.25	.2	4158	3866	2152
#.2	.2	4727	3812	1629

#.2	.9	578	593	615	.73	.88
#.2	.7884	848	870	887	.71	.91
#.2	.5	2446	2445	1987	.69	.96
#.2	.4	3285	3139	2100
#.2	.3	4132	3647	1928
#.2	.2	4727	3812	1629

#.3	.3	2573	2815	2374
#.2	.3	4132	3647	1928
#.17	.3	4418	3574	1591
#.15	.3	4435	3367	1321
#.1	.3	3463	2170	602

#.15	.3	4435	3367	1321
#.15	.4	3977	3358	1716
#.13	.4	4098	3250	1462
#.1	.4	3867	2744	991
#.1	.5	3758	3008	1449



