#!/bin/sh

# $Id$

# Purpose: Commit and/or update multiple SVN repositories

# Usage: 
# cd ~/sh;
# cd ~/sh;./svn.sh > ~/svn.txt 2>&1

# Distribution: 
# scp ~/sh/cvs.sh pbs.ess.uci.edu:
# scp virga.ess.uci.edu:sh/cvs.sh ~

export CVSROOT=':ext:${USER}@dust.ess.uci.edu:/data/home/${USER}/cvs'
export SVNROOT='svn+ssh://dust.ess.uci.edu/data/home/${USER}/svn/trunk'

# 0. Command line options
svn_md='ci' # [sng] SVN mode ('ci' for checkin, 'co' for checkout)

if [ -n "${1}" ]; then
    svn_md=${1}
fi # !$1

# Loop over directories and check them in or out
cd
for drc in `ls --hide=Desktop`; do
    if [ -d ${drc}/.svn ]; then 
	if [ ${svn_md} = 'co' ]; then
	    echo "Updating ${drc} ..."
	    svn update ${drc}
	elif [ ${svn_md} = 'ci' ]; then
	    echo "Committing ${drc} ..."
	    svn upgrade ${drc} # Necessary at FC20 upgrade, and at Kubuntu 14.04 upgrade
	    svn commit -m "" ${drc}
	elif [ ${svn_md} = 'cz' ]; then
	    # 20120531 fix out-of-sync svn repository
	    cd
	    /bin/rm -r -f /tmp/${drc}
	    /bin/rm -r -f /tmp/${drc}.new
	    /bin/mv -f ${drc} /tmp/${drc}.new
	    cd /tmp
	    svn co ${SVNROOT}/${drc}
	    rsync --del --exclude='.svn' ${drc}.new/* ${drc}
	    svn commit -m "" ${drc}
	    /bin/mv -f ${drc} ${HOME}
	    /bin/rm -r -f /tmp/${drc}
	    /bin/rm -r -f /tmp/${drc}.new
	    cd
	else
	    echo "Bad option, exiting ..."
	    exit
	fi
    fi # endif
done # end loop over ${drc}

