# New aerosol defaults for SWNB
mie -3 --dbg=0 --cmp_prt=h2o_lqd --sz_grd=log --sz_mnm=0.1 --sz_mxm=30.0 --sz_nbr=300 --rds_swa=10.0 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_h2o_lqd_rds_swa_10.nc &

mie -3 --dbg=0 --cmp_prt=h2o_ice --sz_grd=log --sz_mnm=0.1 --sz_mxm=50.0 --sz_nbr=300 --rds_swa=20.0 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_h2o_ice_rds_swa_20.nc &

mie -3 --dbg=0 --cmp_prt=sulfate --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=1.0 --sz_nbr=400 --rds_nma=0.08 --gsd_anl=1.6 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_sulfate.nc &

mie -3 --dbg=0 --cmp_prt=lac_ChC90 --dns_prt=1800.0 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_vma=0.1 --gsd_anl=2.0 --wvl_mnm=0.2 --wvl_mxm=5.0 --wvl_nbr=480 --bnd_nbr=1 --lgn_nbr=902 --ngl_nbr=361 ${DATA}/aca/aer_lac_ChC90.nc &
