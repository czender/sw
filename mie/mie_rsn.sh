$Id$ -*-Fundamental-*-

Purpose: Document investigation of Mie resonances

# Organization
1. Frequently used environment variables
2. Test thread number equivalence: one thread vs. eight threads
3. Test thread model equivalence: Inner loop (CL3) vs. Middle loop (CL2)
4. Convergence tests
   4a. Wavelength wvl: multiple wavelengths, one band, one size
   4b. Band       bnd: one wavelength, multiple bands (one thread), one size
   4c. Band   bnd_thr: one wavelength, multiple bands (multiple threads), one size
   4d. Hybrid     hyb: multiple wavelengths, multiple bands, one size
   4e. Size        sz: one wavelength, multiple bands, multiple sizes
5. Intercompare convergence test results
6. Optical property aggregation
7. Compare resonance to low resolution results
8. Sample all size parameters with constant index of refraction
9. Locate "geographic" maxima
10. Sample resonance of maximum amplitude with varying index of refraction
11.  Water vapor-only absorption and pure scattering cloud droplets
12. Miscellaneous runs

# 1. Frequently used environment variables
export DATA='/data/zender' # UCI
export DATA='/data/zender/zender' # UCI
export DATA='/ptmp/zender' # ESMF
export CASEID='mie16'
export SZ_SNG='rds_swa_10'

# 2. Test thread number equivalence: one thread vs. eight threads
mie -dbg --slf_tst=BoH83
for thr_nbr in 1 8 ; do
export OMP_NUM_THREADS=${thr_nbr}
mie --dbg=1 --bnd_nbr=2 --sz_nbr=2 ${DATA}/aca/foo_${OMP_NUM_THREADS}.nc
done
ncdiff -O ${DATA}/aca/foo_8.nc ${DATA}/aca/foo_1.nc ${DATA}/aca/foo_8m1.nc
ncks -m -C -v abs_cff_mss ${DATA}/aca/foo_8m1.nc | m

# 3. Test thread model equivalence: Inner loop (CL3) vs. Middle loop (CL2)
xlC_r -DPARALLELIZE_OVER_CL2 -DPRC_DBL -DABORT_ON_ERROR -DAIX -DNO_NETCDF_2 -DVERSION='20040620' -DHOSTNAME='esmf04m' -DUSER='zender'  -I/usr/vacpp/include -DHAVE_CSTDLIB -I/usr/vacpp/include -I. -I/u/zender/include -I/usr/local/include  -I/usr/local/include -qmaxmem=8192 -qspill=2048 -g -q64 -qsmp=omp -c mie.cc -o /u/zender/obj/AIX/mie.o
make OPTS=D
xlC_r -DPARALLELIZE_OVER_CL3 -DPRC_DBL -DABORT_ON_ERROR -DAIX -DNO_NETCDF_2 -DVERSION='20040620' -DHOSTNAME='esmf04m' -DUSER='zender'  -I/usr/vacpp/include -DHAVE_CSTDLIB -I/usr/vacpp/include -I. -I/u/zender/include -I/usr/local/include  -I/usr/local/include -qmaxmem=8192 -qspill=2048 -g -q64 -qsmp=omp -c mie.cc -o /u/zender/obj/AIX/mie.o
make OPTS=D

# Convergence tests
# 4a. Wavelength wvl: multiple wavelengths, one band, one size
# NB: Use only one thread when bnd_nbr=1
# NB: Set dmn_nbr_max=1 so that phase function arrays do not swamp output
tst_nm=mie_rsn_tst_4a
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --thr_nbr=1 --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_${wvl_nbr}.nc > ${DATA}/mie/foo_${wvl_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_${wvl_nbr}.nc ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_${wvl_nbr}.nc ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${wvl_nbr} wavelengths...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over wvl_nbr

# 4b. Band       bnd: one wavelength, multiple bands (one thread), one size
# Run with only one thread to check for BFB agreement with wavelengths
tst_nm=mie_rsn_tst_4b
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
#for bnd_nbr in 1 10 100 1000 10000 100000 1000000 10000000 100000000; do 
for bnd_nbr in 100000000; do 
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --thr_nbr=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_bnd_${bnd_nbr}.nc > ${DATA}/mie/foo_bnd_${bnd_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_bnd_${bnd_nbr}.nc ${DATA}/mie/mie_avg_bnd_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_bnd_${bnd_nbr}.nc ${DATA}/mie/mie_avg_bnd_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${bnd_nbr} bands...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_bnd_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over bnd_nbr

# 4c. Band       bnd: one wavelength, multiple bands (multiple threads), one size
# Run with only one thread to check for BFB agreement with wavelengths
tst_nm=mie_rsn_tst_4c
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
#for bnd_nbr in 1 10 100 1000 10000 100000 1000000 10000000 100000000 ; do 
for bnd_nbr in 100000000; do 
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --thr_nbr=4 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_bnd_thr_${bnd_nbr}.nc > ${DATA}/mie/foo_bnd_thr_${bnd_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_bnd_thr_${bnd_nbr}.nc ${DATA}/mie/mie_avg_bnd_thr_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_bnd_thr_${bnd_nbr}.nc ${DATA}/mie/mie_avg_bnd_thr_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${bnd_nbr} bands...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_bnd_thr_${bnd_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over bnd_nbr

# 4d. Hybrid     hyb: multiple wavelengths, multiple bands, one size
tst_nm=mie_rsn_tst_4d
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for wvl_nbr in 10 ; do 
for bnd_nbr in 100000 ; do 
  bin_nbr=$(echo "${wvl_nbr}*${bnd_nbr}" | bc -l)
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc > ${DATA}/mie/foo_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${wvl_nbr} wavelengths, ${bnd_nbr} bands, ${bin_nbr} bins...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over bnd_nbr
done # end loop over wvl_nbr

# 4e. Size     sz: one wavelength, multiple bands, multiple sizes
tst_nm=mie_rsn_tst_4e
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for bnd_nbr in 10 ; do
for sz_nbr in 1 10 100 1000 10000 100000; do
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --bch_flg --psd_typ=lognormal --sz_grd=log --sz_nbr=${sz_nbr} --sz_mnm=1.0 --sz_mxm=30.0 --rds_swa=10.0 --gsd_anl=1.6 --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_bnd_sz_${bnd_nbr}_${sz_nbr}.nc > ${DATA}/mie/foo_bnd_sz_${bnd_nbr}_${sz_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_bnd_sz_${bnd_nbr}_${sz_nbr}.nc ${DATA}/mie/mie_avg_bnd_sz_${bnd_nbr}_${sz_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_bnd_sz_${bnd_nbr}_${sz_nbr}.nc ${DATA}/mie/mie_avg_bnd_sz_${bnd_nbr}_${sz_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${bnd_nbr} bands, ${sz_nbr} sizes...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_bnd_sz_${bnd_nbr}_${sz_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over sz_nbr
done # end loop over bnd_nbr

# Complete hybrid run for liquid cloud droplets
for wvl_nbr in 10 ; do 
for bnd_nbr in 100000 ; do 
  bin_nbr=$(echo "${wvl_nbr}*${bnd_nbr}" | bc -l)
  printf "Simulation started at `date`\n"
  mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc > ${DATA}/mie/foo_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${wvl_nbr} wavelengths, ${bnd_nbr} bands, ${bin_nbr} bins...\n"
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_hyb_${wvl_nbr}_${bnd_nbr}_${bin_nbr}.nc
  printf "\nSimulation finished at `date`\n"
done # end loop over bnd_nbr
done # end loop over wvl_nbr

# 5. Intercompare convergence test results
# Prove wvl, bnd, and hyb results are identical
for bin_nbr in 1 10 100 1000 10000 100000 1000000 ; do 
  bnd=`ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss ${DATA}/mie/mie_bnd_${bin_nbr}.nc`
  wvl=`ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss ${DATA}/mie/mie_avg_${bin_nbr}.nc`
  hyb=`ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss ${DATA}/mie/mie_avg_hyb_10_100000_${bin_nbr}.nc`
  printf "${bin_nbr} bands abs_cff_mss: wvl ${wvl}, bnd ${bnd}, hyb ${hyb}\n"
done # end loop over bin_nbr

# 6. Optical property aggregation
cd ${DATA}/mie/${CASEID}
# Select variables (extracted variable list should match mie.pro:aer_gph() inputs)
var_lst_mie='abs_cff_mss,abs_fsh_ffc,ang_xpn,asm_prm,dmt_ctr,dmt_nmr,ext_cff_mss,flx_slr_frc,sca_cff_mss,ss_alb,ss_co_alb,sz_prm_rsn_avg,wvl_ctr,wvl_max,wvl_min'
# Perform multi-file concatenation using stdin for input
/bin/ls | grep aer_h2o_lqd_${SZ_SNG}_${CASEID}_'0[0-4]\.[0-9][0-9]00000_0[0-5]\.[0-9][0-9]00000.nc' | ncecat -D 2 -O -v ${var_lst_mie} -p ${DATA}/mie/${CASEID} -o ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
# Average over degenerate wavelength dimension
ncwa -O -a wvl ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
# Ensemble dimension is now wavelength dimension
ncrename -O -d record,wvl ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
ncap -s 'wvl=wvl_ctr' -O ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
scp esmf.ess.uci.edu:/ptmp/zender/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/aca

# Generate low resolution optical properties
mie --dbg=1 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=30 --rds_swa=7.0 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 ${DATA}/aca/aer_h2o_lqd_rds_swa_07_lrz.nc
mie --dbg=1 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=30 --rds_swa=7.0 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=10 ${DATA}/aca/aer_h2o_lqd_rds_swa_07_mrz.nc
mie --dbg=1 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=30 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 ${DATA}/aca/aer_h2o_lqd_rds_swa_10_lrz.nc
ncks -m -C -v ext_cff_mss,abs_cff_mss,sca_cff_mss,asm_prm,ss_alb ${DATA}/aca/aer_h2o_lqd_rds_swa_10_lrz.nc | m

# 7. Compare resonance to low resolution results
for CASEID in lrz mie10 mie11 mie12 mie13 mie14 mie15 mie16 ; do 
  var_lst_mie='abs_cff_mss,abs_fsh_ffc,ang_xpn,asm_prm,dmt_ctr,dmt_nmr,ext_cff_mss,flx_slr_frc,sca_cff_mss,ss_alb,ss_co_alb,sz_prm_rsn_avg,wvl_ctr,wvl_max,wvl_min'
  ncdiff -O -v ${var_lst_mie} ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_lrz.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}mlrz.nc
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_alb for CASEID = ${CASEID}...\n"
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_alb ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc
  printf "\n"
done # end loop over CASEID

# Longitudinal aggregation of optical properties
for CASEID in mie10 mie11 mie12 mie13 mie14 mie15 mie16 ; do
  export ${CASEID}
  scp esmf.ess.uci.edu:/ptmp/zender/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/aca
  ncwa -O -a wvl -w flx_slr_frc ${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_${CASEID}_avg.nc
done # end loop over CASEID
ncecat -D 2 -O -v ${var_lst_mie} ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_mie1?_avg.nc ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_avg_mie10_mie16.nc
ncrename -O -d record,rsn ${DATA}/ppr_ZeT05/aer_h2o_lqd_${SZ_SNG}_avg_mie10_mie16.nc
# RMS percent relative bias
#for CASEID in mie10 mie11 mie12 mie13 mie14 mie15 mie16 ; do
# fxm: never quite got this debugged
for CASEID in mie10; do
  export ${CASEID}
  fl_ctl=${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_mie16.nc
  fl_ctl_tmp=${fl_ctl/.nc/_tmp.nc}
  fl_xpt=${DATA}/aca/aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
  fl_rms=${fl_xpt/.nc/_rms.nc}
  /bin/cp ${fl_ctl} ${fl_ctl_tmp}
  ncks -O -v abs_cff_mss ${fl_ctl} ${fl_ctl_tmp}
  ncrename -v abs_cff_mss,abs_cff_mss_ctl ${fl_ctl_tmp}
  ncks -A -C -v abs_cff_mss_ctl ${fl_ctl_tmp} ${fl_rms}
  ncap -O -s "abs_cff_mss_err_rlt_pct=100.0*(abs_cff_mss-abs_cff_mss_ctl)/abs_cff_mss_ctl" -s "abs_cff_mss_err_rlt_pct_wgt_flx=flx_slr_spc*abs_cff_mss_err_rlt_pct" ${fl_rms} ${fl_rms}
  ncwa -O -y rms -a wvl ${fl_rms} ${fl_rms}
  printf "Control file: ${fl_ctl}\n"
  printf "Experiment file: ${fl_xpt}\n"
  printf "RMS percent relative bias total(sqrt((100*(abs_cff_mss_xpt-abs_cff_mss_ctl)/abs_cff_mss_ctl))^2))/ord_nbr= "
  ncks -v abs_cff_mss -s "%f" ${fl_rms}
done # end loop over CASEID

# Simulate radiative transfer through cloudy atmosphere with SWNB2, intercompare results
# Longitudinal aggregation of radiative forcing
for CASEID in mie10 mie11 mie12 mie13 mie14 mie15 mie16 ; do
  mpc_CWP=0.1
  SZ_SNG='rds_swa_10'
# Control
  swnb2 --dbg=1 --alb_sfc=0.1 --frc_lqd --mpc_CWP=${mpc_CWP} --slr_zen_ngl_cos=0.5 --drc_in=${DATA}/aca --drc_out=${DATA}/ppr_ZeT05  --fl_clm=mls_cld.nc --fl_lqd=aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc --fl_out=swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc
done # end loop over CASEID
var_lst_swnb2='flx_bb_abs,flx_bb_abs_atm,flx_bb_abs_ttl,flx_bb_abs_sfc,abs_bb_SAS,rfl_bb_SAS,trn_bb_atm,mpc_CWP'
ncecat -D 2 -O -d lev,"800.0 millibars" -v ${var_lst_swnb2} ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_mie1?.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_mie10_mie16.nc
ncrename -O -d record,rsn -v flx_bb_abs,flx_bb_abs_cld ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_mie10_mie16.nc
ncwa -O -a lev ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_mie10_mie16.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_mie10_mie16.nc

# 8. Sample all size parameters with constant index of refraction
tst_nm=mie_rsn_tst_8
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --thr_nbr=1 --no_wrn_ntp --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.1 --wvl_mxm=100.0 --cmp_prt=h2o_lqd --idx_rfr_prt=1.33+1.0e-6i ${DATA}/mie/mie_szp_${wvl_nbr}.nc > ${DATA}/mie/foo_szp_${wvl_nbr}.txt 2>&1 
  ncwa -O -a wvl -w flx_slr_frc -v sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb ${DATA}/mie/mie_szp_${wvl_nbr}.nc ${DATA}/mie/mie_avg_szp_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  ncap -O -s "ss_alb=sca_cff_mss/ext_cff_mss" -s "ss_co_alb=1.0d-ss_alb" ${DATA}/mie/mie_avg_szp_${wvl_nbr}.nc ${DATA}/mie/mie_avg_szp_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nMean abs_cff_mss, sca_cff_mss, ss_co_alb with ${wvl_nbr} wavelengths...\n" >> ${tst_rsl} 2>&1
  ncks -s "%18.12e " -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_co_alb ${DATA}/mie/mie_avg_szp_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
  printf "\nSimulation finished at `date`\n" >> ${tst_rsl} 2>&1
done # end loop over wvl_nbr
# Reduce to manageable size for copying to laptops
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
  var_lst_mie='abs_cff_mss,abs_fsh_ffc,ang_xpn,asm_prm,dmt_ctr,dmt_nmr,ext_cff_mss,flx_slr_frc,sca_cff_mss,ss_alb,ss_co_alb,sz_prm_rsn_avg,wvl_ctr,wvl_max,wvl_min'
  ncks -v ${var_lst_mie} ${DATA}/mie/mie_szp_${wvl_nbr}.nc ${DATA}/ppr_ZeT05/mie_szp_${wvl_nbr}.nc >> ${tst_rsl} 2>&1
done # end loop over wvl_nbr
# Copy to laptops
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
  export ${wvl_nbr}
  scp esmf.ess.uci.edu:/ptmp/zender/ppr_ZeT05/mie_szp_${wvl_nbr}.nc ${DATA}/ppr_ZeT05
done # end loop over wvl_nbr

# 9. Locate geographic maxima
tst_nm=mie_rsn_tst_9
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
for wvl_nbr in 1 10 100 1000 10000 100000 1000000 10000000 ; do 
  printf "Processing started for \$wvl_nbr = ${wvl_nbr} at `date`\n" >> ${tst_rsl} 2>&1
  var_lst_mie='abs_cff_mss,abs_fsh_ffc,ang_xpn,asm_prm,dmt_ctr,dmt_nmr,ext_cff_mss,flx_slr_frc,sca_cff_mss,ss_alb,ss_co_alb,sz_prm_rsn_avg,wvl_ctr,wvl_max,wvl_min'
  ncwa -D 0 -O -y max -a wvl -v ${var_lst_mie} ${DATA}/ppr_ZeT05/mie_szp_${wvl_nbr}.nc ${DATA}/ppr_ZeT05/mie_szp_${wvl_nbr}_max.nc >> ${tst_rsl} 2>&1
  var_lst_ncks='abs_fsh_ffc,sca_cff_mss,ext_cff_mss,abs_cff_mss,ss_alb,ss_co_alb'
  printf "\nMaximum values for \$wvl_nbr = ${wvl_nbr}...\n" >> ${tst_rsl} 2>&1
  ncks -m -C -v ${var_lst_ncks} ${DATA}/ppr_ZeT05/mie_szp_${wvl_nbr}_max.nc >> ${tst_rsl} 2>&1
done # end loop over wvl_nbr
CASEID=mie21 ; SZ_SNG='rds_swa_10'
cd ${DATA}/mie/${CASEID}
# Perform multi-file concatenation using stdin for input
/bin/ls | grep aer_h2o_lqd_${SZ_SNG}_${CASEID}_'0[0-4]\.[0-9][0-9]00000_0[0-5]\.[0-9][0-9]00000.nc' | ncrcat -D 1 -v abs_cff_mss -o aer_h2o_lqd_${SZ_SNG}_${CASEID}_abs_cff_mss.nc
# Average over degenerate wavelength dimension
ncwa -D 0 -O -y max -a wvl ${DATA}/mie/${CASEID}/aer_h2o_lqd_${SZ_SNG}_${CASEID}_abs_cff_mss.nc ${DATA}/mie/${CASEID}/aer_h2o_lqd_${SZ_SNG}_${CASEID}_abs_cff_mss_max.nc
GREP_TRG='1.73268686115'
ncks ${DATA}/mie/${CASEID}/aer_h2o_lqd_${SZ_SNG}_${CASEID}_abs_cff_mss.nc | grep ${GREP_TRG}
# wvl[9849281]=1.18492815e-06 abs_cff_mss[9849281]=1.73268686115

# 10. Sample resonance of maximum amplitude with varying index of refraction
tst_nm=mie_rsn_tst_10
tst_rsl=${HOME}/mie/${tst_nm}.txt
/bin/rm -f ${tst_rsl} # Clean up old test results
CASEID=mie22
mkdir -p ${DATA}/mie/${CASEID};cd ${DATA}/mie/${CASEID}
for idx_rfr_img in 1.0e-1 1.0e-2 1.0e-3 1.0e-4 1.0e-5 1.0e-6 1.0e-7 1.0e-8 1.0e-9 1.0e-10 ; do
  wvl_mnm=1.18
  wvl_mxm=1.19
  wvl_nbr=100000
  printf "Simulation started at `date`\n" >> ${tst_rsl} 2>&1
  mie --dbg=1 --tst=nsz --thr_nbr=1 --no_wrn_ntp --dmn_nbr_max=1 --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=${wvl_mnm} --wvl_mxm=${wvl_mxm} --cmp_prt=h2o_lqd --idx_rfr_prt=1.33+${idx_rfr_img}i ${DATA}/mie/${CASEID}/aer_h2o_lqd_rds_swa_10_${CASEID}_`printf  %010.7f ${wvl_mnm}`_`printf  %010.7f ${wvl_mxm}`_idx_rfr_img_`printf  "%07.1e\n" ${idx_rfr_img}`.nc > ${DATA}/mie/foo_idx_rfr_${wvl_nbr}.txt 2>&1 
done # end loop over idx_rfr_img

11. Water vapor-only absorption and pure scattering cloud droplets
export CASEID='mie16'
mpc_CWP=0.1
var_lst_swnb2='flx_bb_abs,flx_bb_abs_atm,flx_bb_abs_ttl,flx_bb_abs_sfc,abs_bb_SAS,rfl_bb_SAS,trn_bb_atm,mpc_CWP'
# Liquid cloud droplets are pure scatterers
swnb2 --sct_lqd=T --dbg=1 --alb_sfc=0.1 --frc_lqd --mpc_CWP=${mpc_CWP} --slr_zen_ngl_cos=0.5 --drc_in=${DATA}/aca --drc_out=${DATA}/ppr_ZeT05  --fl_clm=mls_cld.nc --fl_lqd=aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc --fl_out=swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_sct_lqd.nc
# Change in absorption due to dis-allowing absorption by cloud droplets
ncdiff -D 2 -O -v ${var_lst_swnb2} ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_sct_lqd.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_sct_lqd_mns_ctl.nc
ncks -C -m -F -d lev,"800.0 millibars" -v ${var_lst_swnb2} ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_sct_lqd_mns_ctl.nc
# Dis-allow H2O vapor absorption in cloud
swnb2 --vpr_H2O_abs_cld=F --dbg=1 --alb_sfc=0.1 --frc_lqd --mpc_CWP=${mpc_CWP} --slr_zen_ngl_cos=0.5 --drc_in=${DATA}/aca --drc_out=${DATA}/ppr_ZeT05  --fl_clm=mls_cld.nc --fl_lqd=aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc --fl_out=swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_vpr_H2O_abs_cld_fls.nc
# Change in absorption due to dis-allowing absorption by in-cloud water vapor
ncdiff -D 2 -O -v ${var_lst_swnb2} ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_vpr_H2O_abs_cld_fls.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}.nc ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_vpr_H2O_abs_cld_fls_mns_ctl.nc
ncks -C -m -F -d lev,"800.0 millibars" -v ${var_lst_swnb2} ${DATA}/ppr_ZeT05/swnb2_aer_h2o_lqd_${SZ_SNG}_${CASEID}_vpr_H2O_abs_cld_fls_mns_ctl.nc

# 12. Miscellaneous runs
# Control: wavelength 0.5 um:
wvl_nbr=1000000
mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_${wvl_nbr}.nc > ${DATA}/mie/foo_${wvl_nbr}.txt 2>&1
# Use sub-bands instead of wavelengths, answer should be identical to control
# 20020305: Answers are identical to control for n=1,10000 but diverge at n=100000 for BoH83
# 200404: Perfect behavior for Wis79
bnd_nbr=1000000
mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --bnd_nbr=${bnd_nbr} --wvl_mnm=0.5 --wvl_mxm=0.51 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_bnd_${bnd_nbr}.nc > ${DATA}/mie/foo_bnd_${bnd_nbr}.txt 2>&1

# Same bandwidth and same refractive index but at 1 um
wvl_nbr=1000000
mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=1.0 --wvl_mxm=1.01 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_1um_${wvl_nbr}.nc > ${DATA}/mie/foo_1um_${wvl_nbr}.txt 2>&1

# Wavelength 1 um with realistic refractive index
wvl_nbr=1000000
mie --dbg=1 --tst=nsz --bch_flg --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=1.0 --wvl_mxm=1.01 --cmp_prt=h2o_lqd ${DATA}/mie/mie_1um_idx_${wvl_nbr}.nc > ${DATA}/mie/foo_1um_idx_${wvl_nbr}.txt 2>&1

# Realistic refractive index and realistic cloud droplet size?
bnd_nbr=1000000
mie --dbg=1 --tst=nsz --bch_flg --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=100 --rds_nma=5.75 --gsd_anl=1.6 --bnd_nbr=${bnd_nbr} --wvl_mnm=1.0 --wvl_mxm=1.01 --cmp_prt=h2o_lqd ${DATA}/mie/mie_1um_cld_${bnd_nbr}.nc > ${DATA}/mie/foo_1um_cld_${bnd_nbr}.txt 2>&1

IDL:
; ctl:
x=[1,10,100,1000,10000,100000,1000000,10000000]
y=[0.000037371738,0.000036300325,0.000040294827,0.000039221448,0.000040019926,0.000043796059,0.000044184975,0.000044184975] ; Good Wis79 numbers
ssa=[0.999999885395d0,0.999999887222d0,0.999999875348d0,0.999999878414d0,0.999999875941d0,0.999999864506d0,0.999999863388d0,0.999999863388d0] ; Wis79 Single scatter-albedo
coalb=1.0d0-ssa
y=[0.000026003911,0.000027036393,0.000027935862,0.000027786377,0.000036847318,0.000028900101,0.000029034384] ; Buggy BoH83 numbers
plot,alog10(x),y,xrange=[0,6],yrange=[0.,max(y)]
; 1um:
y=[0.000039820237,0.000025370362,0.000024686913,0.000024689657,0.000024689657,0.000024689672,0.000024689670] ; Buggy BoH83 numbers
; 1um_cld:
y=[0.039527551115,0.044616916048,0.042325079809,0.042232021500,0.042267865309,0.042267814322,0.042267819885] ; Buggy BoH83 numbers

Following commands are interactive, not batch, so do not invoke batch behavior
Zoom in, use 'hrz' prefix:
wvl_nbr=1000
# Old resonance for buggy BoH83 code. Wis79 proves no resonance here.
mie --dbg=1 --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5036275 --wvl_mxm=0.5036295 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hrz_${wvl_nbr}.nc > ${DATA}/mie/foo_hrz_${wvl_nbr}.txt 2>&1

wvl_nbr=1000
# New resonances identified by Wis79 code (increase n_i to grow resonance)
# Strongest peak between 0.50 and 0.51 um:
mie --dbg=1 --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.50361540 --wvl_mxm=0.50361700 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hrz_${wvl_nbr}.nc > ${DATA}/mie/foo_hrz_${wvl_nbr}.txt 2>&1

# Moderate peak very close to 0.50 um:
mie --dbg=1 --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5000039 --wvl_mxm=0.5000043 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hrz_${wvl_nbr}.nc > ${DATA}/mie/foo_hrz_${wvl_nbr}.txt 2>&1

# Weak peak very close to 0.55 um:
mie --dbg=1 --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.555106 --wvl_mxm=0.555109 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_hrz_${wvl_nbr}.nc > ${DATA}/mie/foo_hrz_${wvl_nbr}.txt 2>&1

Zoom out, use 'lrz' prefix:
wvl_nbr=1000
mie --dbg=1 --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=${wvl_nbr} --wvl_mnm=0.5 --wvl_mxm=1.0 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-9)' ${DATA}/mie/mie_lrz_${wvl_nbr}.nc > ${DATA}/mie/foo_lrz_${wvl_nbr}.txt 2>&1

Size parameter ranges from 2*pi*5/0.5 = 20*pi ~ 66 to 22
mie --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=10 --wvl_mnm=0.5 --wvl_mxm=1.5
ncks -C -F -u -m -v abs_cff_mss,sca_cff_mss,ss_alb ${DATA}/mie/mie.nc | m

# Simulate resonance cases in literature:
// Nus03 p. 1590 Figure 2 m=1.317+1.155e-5i near X = 50
mie --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=1000 --wvl_mnm=1.247 --wvl_mxm=1.267 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.317,1.155e-5)' ${DATA}/mie/mie.nc > ${DATA}/mie/foo.txt 2>&1

// BoH83 p. 300 Figure 11.7, 11.8 from CKK78a m = 1.33+1.0e-8i near X = 50.337
mie --tst=nsz --sz_mnm=9.999 --sz_mxm=10.001 --sz_nbr=1 --wvl_nbr=100 --wvl_mnm=0.06240872 --wvl_mxm=0.06241244 --cmp_prt=h2o_lqd --idx_rfr_prt='(1.33,1.0e-8)' ${DATA}/mie/mie.nc > ${DATA}/mie/foo.txt 2>&1
