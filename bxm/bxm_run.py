#!/usr/bin/python

# $Id$

# Purpose: Time-interpolate user input, create tvbds.cdl, run bxm_run.sh
# bxm_run.py is called by bxm_cnf.html
# bxm_run.py creates tvbds.cdl
# bxm_run.py calls bxm_run.sh

# Usage:
# sudo scp ~/bxm/bxm_run.py dust.ess.uci.edu:/var/www/cgi-bin/bxm_run.py

# Permissions:
# sudo chmod 744 /var/www/cgi-bin/bxm_run.py
# sudo chown apache /var/www/cgi-bin/bxm_run.py
# sudo chgrp apache /var/www/cgi-bin/bxm_run.py

import sys,os,cgi,re,time
from string import split,atof

# Create tvbds.cdl from form
outp = open('/tmp/bxm/tvbds.cdl','w')
fld_lst = ['asp_rat_lps','hgt_mdp','oro','prs_mdp','prs_ntf','q_H2O_vpr','sfc_typ','tpt_gnd','tpt_ice','tpt_mdp','tpt_soi','tpt_sst','vai_dst','wnd_mrd_mdp','wnd_znl_mdp']

# Parse input from form
form = cgi.FieldStorage()

dbg_lvl = '0'
if form.has_key('dbg_lvl'):
    dbg_lvl = re.sub('[^0-9e.-]','',form['dbg_lvl'].value)

if form.has_key('time_nbr'):
    time_nbr = re.sub('[^0-9e.-]','',form['time_nbr'].value)

msg_usr = '            '
if form.has_key('msg_usr'):
    msg_usr = form['msg_usr'].value

# Read time_nbr
time_nbr_sng = time_nbr
time_nbr = int(time_nbr)
nc_hdr = """// Purpose: Time Varying Boundary DataSet (tvbds) input for DEAD box model
// tvbds.cdl is automagically generated from form input to bxm_cnf.html by bxm_run.py
// tvbds.cdl is then converted to tvbds.nc by ncgen

// Usage:
// ncgen -b -o tvbds.nc tvbds.cdl

netcdf tvbds{
dimensions:
  time = %d;
variables:
  long time(time);\n"""
outp.write(nc_hdr % time_nbr)

# Output
for fld_nm in fld_lst:
    outp.write('  float %s(time);\n' % fld_nm)

outp.write('data:\n')
outp.write('  time = ')
for time_idx in range(time_nbr-1):
    outp.write('%d,' % time_idx)
outp.write('%d;\n' % (time_nbr-1))

# Interpolate time_nbr values
for fld_nm in fld_lst:
    outp.write('  %s = ' % fld_nm)
    for time_idx in range(time_nbr-1):
        val_srt = atof(form[fld_nm].value)
        val_end = atof(form[fld_nm+'1'].value)
        val_ntp = ((time_nbr-time_idx-1)*val_srt + time_idx*val_end)/(time_nbr-1.0)
        outp.write('%f,' % val_ntp)
    time_idx = time_nbr-1
    val_ntp = ((time_nbr-time_idx-1)*val_srt + time_idx*val_end)/(time_nbr-1.0)
    outp.write('%f;\n' % val_ntp)

outp.write('}\n')
outp.close()

# Execute model run script /tmp/bxm/bxm_run.sh
pid = os.fork()
if pid == 0:
    os.environ['PATH']='/bin:/tmp/bxm'
    os.environ['LD_LIBRARY_PATH']='/usr/local/lf9561/lib'
    os.environ['NCARG_ROOT']='/usr/local'
    os.execlp('bxm_run.sh','bxm_run.sh',dbg_lvl,msg_usr,time_nbr_sng)
else:
    os.wait()
    print "Status: 302 Moved"
# fxm: 20030825 I have no idea why, but this Status: line seems to be required
#    print "Status: Error running bxm_run.sh from bxm_run.py"
    print "Location: /dead/bxm_vzn.html"
    print
