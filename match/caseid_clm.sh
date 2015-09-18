#!/bin/bash

# Purpose: Generate climatology:
#             caseid_clm_0112.nc
#             caseid_clm.nc
#             caseid_clm_xy.nc
#          from caseid_yyyymm.nc files
#          (assumes that files are in /dust/zender/caseid/)
# Examples:
#   caseid_clm.sh dstmch46
#   caseid_clm.sh dstmch84

#---------------------------------------------------------------
usage="usage: $0 caseid"
if [ $# -ne 1 ]; then
    echo $usage
    exit 1
fi
caseid=$1
#---------------------------------------------------------------

path=/dust/zender/$caseid

# calculate additional variables
for file in ${path}/${caseid}_[12][0-9][0-9][0-9][01][0-9].nc ; do
  echo "ncap -O -S /dhome/zender/dst/dst.nco $file $file"
  ncap -O -S /dhome/zender/dst/dst.nco $file $file
done

var="gw,DSTQ,DSTODXC,DSTMPC,DSTSFGRV,DSTSFPCP,DSTSFTRB,DSTSFDPS,DSTSFDRY,DPSWETFRC,DSTSFMBL,PRECT,VWC_SFC,WND_FRC,ORO"
var_sdn="DSTQ_sdn,DSTODXC_sdn,DSTMPC_sdn,DSTSFGRV_sdn,DSTSFPCP_sdn,DSTSFTRB_sdn,DSTSFDPS_sdn,DSTSFDRY_sdn,DPSWETFRC_sdn,DSTSFMBL_sdn,PRECT_sdn,VWC_SFC_sdn,WND_FRC_sdn"

for mm in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
    infile=${path}/${caseid}_????${mm}.nc
    # compute average
    outfile=tmp_${mm}.nc
    echo "ncra -O $infile $outfile"
    ncra -O -v $var $infile $outfile
    count=`ls $infile | wc -l`
    echo "num years = $count"
    if [ $count -gt 1 ]; then
	# compute standard deviation (sdn)
	outfile=tmp_${mm}_sdn.nc
	ncrcat -O -v $var $infile tmp1.nc
	ncwa   -O -a time         tmp1.nc tmp2.nc
	ncdiff -O -v $var         tmp1.nc tmp2.nc tmp3.nc
	ncra   -O -y rmssdn                       tmp3.nc $outfile
        ncrename -v DSTQ,DSTQ_sdn                         $outfile
        ncrename -v DSTODXC,DSTODXC_sdn                   $outfile
        ncrename -v DSTMPC,DSTMPC_sdn                     $outfile
        ncrename -v DSTSFGRV,DSTSFGRV_sdn                 $outfile
        ncrename -v DSTSFPCP,DSTSFPCP_sdn                 $outfile
        ncrename -v DSTSFTRB,DSTSFTRB_sdn                 $outfile
        ncrename -v DSTSFMBL,DSTSFMBL_sdn                 $outfile
        ncrename -v PRECT,PRECT_sdn                       $outfile
        ncrename -v VWC_SFC,VWC_SFC_sdn                   $outfile
        ncrename -v WND_FRC,WND_FRC_sdn                   $outfile
        ncrename -v DSTSFDPS,DSTSFDPS_sdn                 $outfile
        ncrename -v DSTSFDRY,DSTSFDRY_sdn                 $outfile
        ncrename -v DPSWETFRC,DPSWETFRC_sdn               $outfile
	# append sdn to mean file
	infile=tmp_${mm}_sdn.nc
	outfile=tmp_${mm}.nc
	echo "ncks -A -C -v $var_sdn $infile $outfile"
	ncks -A -C -v $var_sdn $infile $outfile
    fi
done

#---------------------------------------------------------------
# cat the files together to make caseid_clm_0112.nc
# average the files together to make caseid_clm.nc
#---------------------------------------------------------------
infile=tmp_??.nc
outfile=${caseid}_clm_0112.nc
# count files
declare -i chk
chk=`ls $infile | wc -l`
if [ $chk -ne 12 ]; then
    echo "$0 error: chk != 12";  exit 1;
fi
echo "ncrcat -O $infile $outfile"
ncrcat -O $infile $outfile
/bin/mv -f $outfile $path

infile=tmp_??.nc
outfile=${caseid}_clm.nc
echo "ncra -O $infile $outfile"
ncra -O $infile $outfile
/bin/mv -f $outfile $path

infile=${caseid}_clm.nc
outfile=${caseid}_clm_xy.nc
echo "ncwa -O -a lat,lon -w gw $infile $outfile"
ncwa -O -a lat,lon -w gw $infile $outfile
/bin/mv -f $outfile $path

#---------------------------------------------------------------
# cleanup
#---------------------------------------------------------------
/bin/rm -f tmp_??.nc
/bin/rm -f tmp_??_sdn.nc
/bin/rm -f tmp[123].nc
echo "$0 >>> created files for $caseid"
