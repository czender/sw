#/bin/csh -f
set echo

# $Id$

# Purpose: Run htrn2nb for all important atmospheric molecules
# This creates Malkmus band model files needed by narrow band models
# KiR83 show 5 cm-1 bands are optimal for CO2
# Bri92 p. 11477 uses 5 cm-1 bands for CO2, O3, CH4, N2O and 10 cm-1 bands for H2O
# Kie97 p. 113 recommends 5 cm-1 for CO2, 10 cm-1 for H2O, 5--10 cm-1 for O3, 5 cm-1 for others

# Water Vapor:
#htrn2nb -l 0.0 -h 25000.0 -b 2500 -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/mlk_H2O.nc
#htrn2nb -l 2000.0 -h 17900.0 -b 1590 -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/mlk_H2O.nc

#wvn_min=2000 # [cm-1] Minimum wavenumber
wvn_min=10 # [cm-1] Minimum wavenumber
wvn_max=27000 # [cm-1] Maximum wavenumber
wvn_rsn=5 # [cm-1] Band resolution
wvn_rsn_H2O=10 # [cm-1] Band resolution for H2O
#wvn_rsn=1 # [cm-1] Band resolution
#wvn_rsn_H2O=2 # [cm-1] Band resolution for H2O
wvn_nbr_H2O=$(((wvn_max-wvn_min)/wvn_rsn_H2O))
wvn_nbr=$(((wvn_max-wvn_min)/wvn_rsn))

htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/mlk_H2O.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/1H2_16O.nc -o ${DATA}/aca/mlk_1H2_16O.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/1H2_17O.nc -o ${DATA}/aca/mlk_1H2_17O.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/1H2_18O.nc -o ${DATA}/aca/mlk_1H2_18O.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/1H_2H_16O.nc -o ${DATA}/aca/mlk_1H_2H_16O.nc

# Carbon Dioxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/CO2.nc -o ${DATA}/aca/mlk_CO2.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/12C_16O2.nc -o ${DATA}/aca/mlk_12C_16O2.nc

# Oxygen:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/O2.nc -o ${DATA}/aca/mlk_O2.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/16O2.nc -o ${DATA}/aca/mlk_16O2.nc

# Hydroxyl:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/OH.nc -o ${DATA}/aca/mlk_OH.nc
#htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/16O_1H.nc -o ${DATA}/aca/mlk_16O_1H.nc

# Nitrogen:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/N2.nc -o ${DATA}/aca/mlk_N2.nc

# Ozone:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/O3.nc -o ${DATA}/aca/mlk_O3.nc

# Nitrous Oxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/N2O.nc -o ${DATA}/aca/mlk_N2O.nc

# Nitric Oxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/NO.nc -o ${DATA}/aca/mlk_NO.nc

# Sulfur Dioxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/SO2.nc -o ${DATA}/aca/mlk_SO2.nc

# Nitrogen Dioxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/NO2.nc -o ${DATA}/aca/mlk_NO2.nc

# Hydrogen Peroxide:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr_H2O} -i ${DATA}/hitran/H2O2.nc -o ${DATA}/aca/mlk_H2O2.nc

# Methane:
htrn2nb --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/CH4.nc -o ${DATA}/aca/mlk_CH4.nc

# Carbon Monoxide:
htrn2nb --double --wvn_min=${wvn_min} --wvn_max=${wvn_max} --bnd_nbr=${wvn_nbr} -i ${DATA}/hitran/CO.nc -o ${DATA}/aca/mlk_CO.nc

