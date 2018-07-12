#/bin/sh

# $Id$

# Purpose: Run htrn2nc.pl for all important atmospheric molecules
# This creates line files needed by htrn2nb
# Usage: htrn2nc.sh

#export htrn_fl="${DATA}/hitran/hitran96.txt"
#export htrn_fl="${DATA}/hitran/hitran00.txt"
#export htrn_fl="${DATA}/hitran/hitran08.txt"
#export htrn_fl="${DATA}/hitran/hitran12.txt"
export htrn_fl="${DATA}/hitran/hitran16.txt"

# Water Vapor
${HOME}/sw/aca/htrn2nc.pl --mlc=1 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/H2O.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=1 --iso=1 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/1H2_16O.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=1 --iso=2 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/1H2_18O.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=1 --iso=3 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/1H2_17O.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=1 --iso=4 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/1H_2H_16O.nc

# Carbon Dioxide
${HOME}/sw/aca/htrn2nc.pl --mlc=2 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/CO2.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=2 --iso=1 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/12C_16O2.nc

# Oxygen
${HOME}/sw/aca/htrn2nc.pl --mlc=7 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/O2.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=7 --iso=1 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/16O2.nc

# Hydroxyl
${HOME}/sw/aca/htrn2nc.pl --mlc=13 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/OH.nc
${HOME}/sw/aca/htrn2nc.pl --mlc=13 --iso=1 --wvn_min=2000.0 --wvn_max=17900.0 ${htrn_fl} ${DATA}/hitran/16O_1H.nc

# All the rest:
# Nitrous Oxide
${HOME}/sw/aca/htrn2nc.pl --mlc=4 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/N2O.nc

# Nitric Oxide
${HOME}/sw/aca/htrn2nc.pl --mlc=8 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/NO.nc

# Sulfur Dioxide
${HOME}/sw/aca/htrn2nc.pl --mlc=9 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/SO2.nc

# Nitrogen Dioxide
${HOME}/sw/aca/htrn2nc.pl --mlc=10 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/NO2.nc

# Hydrogen Peroxide
${HOME}/sw/aca/htrn2nc.pl --mlc=25 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/H2O2.nc

# Nitrogen
${HOME}/sw/aca/htrn2nc.pl --mlc=22 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/N2.nc

# Ozone
${HOME}/sw/aca/htrn2nc.pl --mlc=3 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/O3.nc

# Methane
${HOME}/sw/aca/htrn2nc.pl --mlc=6 --wvn_min=1.0e-4 --wvn_max=50000.0 ${htrn_fl} ${DATA}/hitran/CH4.nc
