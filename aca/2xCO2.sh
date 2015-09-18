# $Id$

# Purpose: Perform 2xCO2 experiments with SWNB and CRM

# My swnb2 model shows that doubling CO2 cools the troposphere by 
# 0.5 W m-2 when CO2 changes from 355 -> 710 ppm. SZA = 60 degrees.
# SWNB TOA forcing is +0.2 W m-2. 
# CRM agrees with SWNB at tropopause but has 0.0 W m-2 at TOA

crm < ${DATA}/aca/mls_clr.txt > ${DATA}/tmp/crm_mls_clr.txt
mv -f crm.nc ${DATA}/tmp/crm_mls_clr.nc
crm < ${DATA}/aca/mls_clr_2xCO2.txt > ${DATA}/tmp/crm_mls_clr_2xCO2.txt
mv -f crm.nc ${DATA}/tmp/crm_mls_clr_2xCO2.nc
ncdiff ${DATA}/tmp/crm_mls_clr_2xCO2.nc ${DATA}/tmp/crm_mls_clr.nc ${DATA}/tmp/crm_mls_clr_2xCO2_frc.nc
ncks -H -u -C -F -v flx_SW_dwn,flx_SW_up,flx_SW_abs_atm ${DATA}/tmp/crm_mls_clr_2xCO2_frc.nc

clm -i ${DATA}/aca/mls_clr.txt -o ${DATA}/aca/mls_clr.nc
swnb2 -p ${DATA}/aca/mls_clr.nc -d ${DATA}/tmp/swnb_mls_clr.nc
clm -i ${DATA}/aca/mls_clr_2xCO2.txt -o ${DATA}/aca/mls_clr_2xCO2.nc
swnb2 -p ${DATA}/aca/mls_clr_2xCO2.nc -d ${DATA}/tmp/swnb_mls_clr_2xCO2.nc
ncdiff ${DATA}/tmp/swnb_mls_clr_2xCO2.nc ${DATA}/tmp/swnb_mls_clr.nc ${DATA}/tmp/swnb_mls_clr_2xCO2_frc.nc
ncks -H -u -C -F -v flx_bb_net,flx_bb_abs_atm ${DATA}/tmp/swnb_mls_clr_2xCO2_frc.nc

