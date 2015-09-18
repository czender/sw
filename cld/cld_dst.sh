mkdir /data2/zender/pub/src/cld-5.4
rm -f /data2/zender/pub/src/cld-5.4/*
cp /home/zender/cld/main.c /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/utilities.c /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/etbfct.f /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/ice.c /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/mie_bohren_c.f /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/idx_rfr_h2o_ice_sgw.f /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/ccm2rad_c.f /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/interp.c /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/spline.f /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/defs.h /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/globals.h /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/README.cld /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/trp35_cld_crm.txt /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/mls35_cld_crm.txt /data2/zender/pub/src/cld-5.4
cp /home/zender/cld/Makefile /data2/zender/pub/src/cld-5.4

cd /data2/zender/pub/src
gtar -cvf ./cld-5.4/cld-5.4.tar ./cld-5.4/*
gzip ./cld-5.4/cld-5.4.tar
gtar -tvzf ./cld-5.4/cld-5.4.tar.gz
rm -f /d5/ftp/pub/zender/*cld*
mv ./cld-5.4/cld-5.4.tar.gz /d5/ftp/pub/zender
ls /d5/ftp/pub/zender
#
cp /data/zender/bin/SUN4SOL2/cld /data2/zender/pub/bin/SUN4SOL2
cp /data/zender/bin/SUN4/cld /data2/zender/pub/bin/SUN4
cp /data/zender/bin/SGI5/cld /data2/zender/pub/bin/SGI5
chmod a+x /data2/zender/pub/bin/SUN4SOL2/*
chmod a+x /data2/zender/pub/bin/SUN4/*
chmod a+x /data2/zender/pub/bin/SGI5/*
