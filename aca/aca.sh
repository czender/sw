# Purpose: Recompute the standard datasets employed by the Malkmus narrow band model

# Recompute Malkmus random band parameters for standard gases
htrn2nb -m 1 -n -l 2000.0 -h 17900.0 -b 1590 -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/H2O.nc
htrn2nb -m 2 -n -l 2000.0 -h 17900.0 -b 3180 -i ${DATA}/hitran/CO2.nc -o ${DATA}/aca/CO2.nc
htrn2nb -m 7 -n -l 2000.0 -h 17900.0 -b 1590 -i ${DATA}/hitran/O2.nc -o ${DATA}/aca/O2.nc
htrn2nb -m 13 -n -l 2000.0 -h 17900.0 -b 1590 -i ${DATA}/hitran/OH.nc -o ${DATA}/aca/OH.nc

# Recompute Mie parameters for standard aerosols
#mie -debug --flx=Labs --aer=dust_like --dist=lognormal --wvl_min=.3 --wvl_max=5. --wvl_nbr=47 --bnd_nbr=1 --sz_grd=log --sz_min=1.e-3 --sz_max=5. --sz_nbr=400 --dst_a=.4 --dst_b=.7884 ${DATA}/aca/aer_dust_like.nc

# Recompute column profile data
clm -i ${DATA}/aca/mls_clm_cld.txt -o ${DATA}/aca/mls_clm_cld.nc
clm -i ${DATA}/aca/mls_clm_clr.txt -o ${DATA}/aca/mls_clm_clr.nc
clm -A -l 110 -i ${DATA}/aca/arese_clm_tst.txt -o ${DATA}/aca/arese_clm_tst.nc

# Retranslate BPB's Mie parameters for standard size distributions
mie2nc -n 74 -i ${DATA}/aca/lqd_10.dat -o ${DATA}/aca/lqd_10.nc
mie2nc -n 194 -i ${DATA}/aca/ice_20.dat -o ${DATA}/aca/ice_20.nc

# Retranslate WMO continuum absorption cross-section data for O3
O3 -i ${DATA}/aca/wmo_new.dat -o ${DATA}/aca/O3.nc
