#! /bin/csh -fv

# Purpose: Create IBP distribution

# Usage: ibp_dst.sh

/bin/rm -r -f /data/zender/ibp-3.9
mkdir /data/zender/ibp-3.9
cp /home/zender/idl/ibp_srt.txt /data/zender/ibp-3.9
cp /home/zender/idl/ibp*.pro /data/zender/ibp-3.9
cp /home/zender/idl/ibp*.com /data/zender/ibp-3.9
cp /home/zender/idl/ibp*.txt /data/zender/ibp-3.9
cp /home/zender/idl/*.tbl /data/zender/ibp-3.9
cp /home/zender/idl/ibp_README.txt /data/zender/ibp-3.9
#cp /data/zender/sld012d/sld012d_8589.nc /data/zender/ibp-3.9
cd /data/zender
gtar cvzf ibp-3.9.tar.gz ./ibp-3.9
gtar tvzf ibp-3.9.tar.gz
rsh ftp.cgd.ucar.edu /bin/rm -r -f /ftp/pub/zender/ibp*
rsh ftp.cgd.ucar.edu mkdir /ftp/pub/zender/ibp
rcp ibp-3.9.tar.gz ftp.cgd.ucar.edu:/ftp/pub/zender/ibp
rsh ftp.cgd.ucar.edu ls -l /ftp/pub/zender/ibp

