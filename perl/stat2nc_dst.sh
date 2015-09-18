mkdir /data2/zender/pub/src/stat2nc-1.4
rm -f /data2/zender/pub/src/stat2nc-1.4/*
cp /home/zender/perl/stat2nc /data2/zender/pub/src/stat2nc-1.4
cp /home/zender/perl/README.stat2nc /data2/zender/pub/src/stat2nc-1.4
cp /home/zender/perl/stat2nc_input.txt /data2/zender/pub/src/stat2nc-1.4
cp /home/zender/perl/stat2nc_output.nc.cdl /data2/zender/pub/src/stat2nc-1.4
cp /home/zender/perl/stat2nc_output.nc /data2/zender/pub/src/stat2nc-1.4
rm -f /data/zender/man/man1/stat2nc
cd /data2/zender/pub/src
gtar -cvf ./stat2nc-1.4/stat2nc-1.4.tar ./stat2nc-1.4/*
gzip ./stat2nc-1.4/stat2nc-1.4.tar
gtar -tvzf ./stat2nc-1.4/stat2nc-1.4.tar.gz
rm -f /d5/ftp/pub/zender/*stat2nc*
mv ./stat2nc-1.4/stat2nc-1.4.tar.gz /d5/ftp/pub/zender
ls /d5/ftp/pub/zender
#
cp /home/zender/perl/stat2nc /data2/zender/pub/bin/sh/stat2nc
rm -f /data2/zender/pub/man/man1/stat2nc
ln -s /data2/zender/pub/bin/sh/stat2nc /data2/zender/pub/man/man1/stat2nc
chmod a+x /data2/zender/pub/bin/sh/stat2nc

