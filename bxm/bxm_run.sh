#!/bin/sh

# $Id$

# Purpose: Run box model for a range of parameters
# bxm_run.sh is called by bxm_run.py
# bxm_run.sh calls ncgen, dead, ncks, ncl, bxm.ncl, convert, ps2pdf
# bxm_run.sh creates bxm.jpg from bxm.ps

# Usage:
# sudo scp ~/bxm/bxm_run.sh dust.ess.uci.edu:/tmp/bxm/bxm_run.sh
# NB: Re-direct output to /dev/null to avoid header problems

# Permissions:
# sudo chmod 744 /tmp/bxm/bxm_run.sh
# sudo chown apache /tmp/bxm/bxm_run.sh
# sudo chgrp apache /tmp/bxm/bxm_run.sh

# Parse input
if [ $# -ne 3 ]; then
    echo "Usage: bxm_run.sh dbg_lvl msg_usr time_nbr"
    exit 1
fi
dbg_lvl=$1
msg_usr=$2
time_nbr=$3
cd /tmp/bxm

# Generate netCDF input file tvbds.nc
/usr/local/bin/ncgen -o tvbds.nc tvbds.cdl > /dev/null 2>&1

# Run box model with tvbds.nc
/var/www/html/dead/dead --dbg_lvl=${dbg_lvl} --fl_ext_dat="tvbds.nc" --time_nbr=${time_nbr} > ./bxm_stdout.txt 2>&1
/bin/cp -f dead.nc dst_mss_bdg.nc tvbds.nc /var/ftp/dead                       
#/bin/cp -f dead.nc /var/ftp/dead/dead.pid$$.nc                       
#/bin/cp -f dst_mss_bdg.nc /var/ftp/dead/dst_mss_bdg.pid$$.nc               
#/bin/cp -f tvbds.nc /var/ftp/dead/tvbds.pid$$.nc                       

# Create text results bxm_txt.html
cat > ./bxm_txt.html <<EOF
<!--
Purpose: HTML page automagically created from run output and text dumps of dead.nc and dst_mss_bdg.nc
bxm_txt.html is created automatically by bxm_run.sh
bxm_txt.html is (symbolically) linked to by /var/www/html/dead/bxm_txt.html
bxm_txt.html is referenced by bxm.css

Usage:

-->

<html>

<head><title>DEAD: Box Model Output and Results in Text Format</title></head>
<!-- <body bgcolor=white> -->
[<a href="bxm_home.html">Box Model</a>]
[<a href="bxm_cnf.html">Configure/Run</a>]
[<a href="bxm_vzn.html">Visualize</a>]
[<a href="bxm_txt.html">Text Output</a>]
[<a href="bxm_dwn.html">Download Results</a>]
[<a href="/dead/">DEAD Home</a>]

<dt><a name="Results"></a></dt>
<dd><a href="#stdout">Run-time Standard Output and Error</a></dd>
<dd><a href="#rsl">Numeric Results</a></dd>
<dd><a href="#bdg">Mass Budget Diagnostics</a></dd>
</dt>

<a name="stdout"></a>
<h2>DEAD: Box Model Run-time Standard Output and Error</h2>
<pre>
EOF
cat bxm_stdout.txt >> ./bxm_txt.html
# Dump data in text format
cat >> ./bxm_txt.html <<EOF2
<p>
<a name="rsl"></a>
<h2>DEAD: Box Model Numeric Results in Text Format (ncks -M -H -u dead.nc) </h2>
<p>
EOF2
/usr/local/bin/ncks -M -H -u dead.nc >> ./bxm_txt.html
cat >> ./bxm_txt.html <<EOF3
<a name="bdg"></a>
<h2>DEAD: Box Model Mass Budget Diagnostics (ncks -M -H -u dst_mss_bdg.nc) </h2>
<p>
EOF3
/usr/local/bin/ncks -M -H -u dst_mss_bdg.nc >> ./bxm_txt.html
echo "</pre></html>" >> ./bxm_txt.html
# fxm: Move bxm_txt.html to /var/www/html/dead

# Customize NCL script bxm.ncl
hgt_mdp=`/usr/local/bin/ncks -H dead.nc | grep hgt_mdp  | ./bxm_mch.pl`
crt_sng=`/usr/local/bin/ncks -M dead.nc | grep creation | ./bxm_mch.pl`
sed -e "s/sed_txt_1/${hgt_mdp}/g;s/sed_txt_2/${msg_usr}/g;s/sed_txt_3/${crt_sng}/g" /var/www/html/dead/bxm.ncl > ./bxm_usr.ncl

# Use NCL to produce bxm.ps
/usr/local/ncarg/bin/ncl < bxm_usr.ncl > /dev/null 2>&1
/usr/bin/convert -quality 100 -crop 0x0 bxm.ps bxm.jpg
# 20031216: Copy rather than link to web space since apache complains about links
/bin/cp -f bxm.jpg /var/www/html/dead/bxm.jpg
# 20030421: Apache 2.0.40 does not execute scripts called by other scripts
# There is probably an httpd.conf keyword that allows this but is less secure
# Thus hardcode conversion script
# ps2pdf calls a long change
#/usr/bin/ps2pdf bxm.ps bxm.pdf
#/usr/bin/ps2pdfwr -dCompatibilityLevel=1.2 bxm.ps bxm.pdf
#/usr/bin/gs -dSAFER -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=bxm.pdf -c .setpdfwrite -f bxm.ps
/usr/bin/gs -dSAFER -q -dNOPAUSE -dBATCH -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true -sAutoRotatePages=PageByPage -sColorConversionStrategy=LeaveColorUnchanged -sDEVICE=pdfwrite -sOutputFile=bxm.pdf -c .setpdfwrite -f bxm.ps
/bin/cp -f bxm.ps bxm.pdf /var/ftp/dead                      
#/bin/cp -f bxm.ps /var/ftp/dead/bxm.pid$$.ps                       
#/bin/cp -f bxm.pdf /var/ftp/dead/bxm.pid$$.pdf                      
