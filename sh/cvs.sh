#!/bin/sh

# $Id$

# Purpose: Transition CVS repository with keywords to SVN

# Usage: 
# cd;./cvs.sh > ~/cvs.txt 2>&1
# cd;./cvs.sh > ~/svn.txt 2>&1
# cd;./cvs.sh > ~/rpl.txt 2>&1

# Distribution: 
# scp ~/sh/cvs.sh pbs.ess.uci.edu:
# scp virga.ess.uci.edu:sh/cvs.sh ~

export CVSROOT=':ext:zender@pbs.ess.uci.edu:/home/zender/cvs'
export SVNROOT='svn+ssh://pbs.ess.uci.edu/home/zender/svn/trunk'

if [ '0' = '1' ]; then
    cvs -z3 -d :ext:zender@nco.cvs.sourceforge.net:/cvsroot/nco co -kk nco
    cvs -d :ext:esmf.ess.uci.edu:/home/mtosca/cvs co -kk ppr_TZR08
    cvs -d :ext:esmf.ess.uci.edu:/home/zender/cvs co -kk -r prp_itr -d prp_itr prp_arl
    cvs -d :ext:pbs.ess.uci.edu:/home/scapps/cvs co -kk ppr_CaZ09
    cvs -d :ext:pbs.ess.uci.edu:/home/scapps/cvs co -kk ppr_CaZ09a

# Directories that need non-CVS/SVN data to build, e.g., papers
mch_dst='virga.ess.uci.edu'
scp -r ppr_CaZ08/figures ${mch_dst}:ppr_CaZ08
scp ppr_FlZ06/*.eps ppr_FlZ06/*.pdf ${mch_dst}:ppr_FlZ06
scp ppr_FZR07/*.eps ppr_FZR07/*.pdf ${mch_dst}:ppr_FZR07
scp -r ppr_ZeM07/fgr ${mch_dst}:ppr_ZeM07
scp ppr_ZFM07/pix_alb* ${mch_dst}:ppr_ZFM07
scp -r ppr_ZFM07/fgr ${mch_dst}:ppr_ZFM07
scp -r ppr_ZGD09/fgr ${mch_dst}:ppr_ZGD09
scp -r prp_ids/fgr ${mch_dst}:prp_ids
scp ppr_IGPP06/*.eps ${mch_dst}:ppr_IGPP06

mch_src='virga.ess.uci.edu'
scp -r ${mch_src}:ppr_CaZ08/figures ~/ppr_CaZ08
scp ${mch_src}:ppr_FlZ06/*.eps ppr_FlZ06/*.pdf ~/ppr_FlZ06
scp ${mch_src}:ppr_FZR07/*.eps ppr_FZR07/*.pdf ~/ppr_FZR07
scp -r ${mch_src}:ppr_ZeM07/fgr ~/ppr_ZeM07
scp ${mch_src}:ppr_ZFM07/pix_alb* ~/ppr_ZFM07
scp -r ${mch_src}:ppr_ZFM07/fgr ~/ppr_ZFM07
scp -r ${mch_src}:ppr_ZGD09/fgr ~/ppr_ZGD09
scp -r ${mch_src}:prp_ids/fgr ~/prp_ids
scp ${mch_src}:ppr_IGPP06/*.eps ~/ppr_IGPP06
fi # !0

if [ '0' = '1' ]; then
# Modules that need special CVS conversion treatment
    cvs -d :ext:zender@pbs.ess.uci.edu:/home/zender/cvs co -kk -r match_brnch_dst dead
#    cvs -d :ext:zender@goldhill.cgd.ucar.edu:/fs/cgd/csm/models/CVS.REPOS co -r ccm_brnch_crm -kk crm
#    cvs -d :ext:zender@goldhill.cgd.ucar.edu:/fs/cgd/csm/models/CVS.REPOS co -r ccm_brnch_dst -kk ccm_dst
fi # !0

# Exit stati:
# exit_status = 0: matching lines were found
# exit_status = 1: matching lines were not found
# exit_status = 2: error

# Convert CVS keywords to intermediate form
if [ '0' = '1' ]; then
    cd
    for drc in `ls --hide=CVSROOT cvs`; do
#     for drc in dead; do
	echo "Processing directory ${drc}..."
	cd;cvs co -kk ${drc}
	cd ${drc}
	for fl in `ls --hide=CVS`; do
	    if egrep '\$Header\$' ${fl} > /dev/null || egrep '\$Date\$' ${fl} > /dev/null || egrep '\$Id\$' ${fl} > /dev/null || egrep '\$Revision\$' ${fl} > /dev/null || egrep '\$Name\$' ${fl}
	    then
		echo "${fl} contains CVS keyword(s), prepending with _CVS_ ..."
		perl -pi -e 's/\$Header\$/\$_CVS_Header\$/g;s/\$Date\$/\$_CVS_Date\$/g;s/\$Id\$/\$_CVS_Id\$/g;s/\$Revision\$/\$_CVS_Revision\$/g;s/\$Name\$/\$_CVS_Name\$/g;' ${fl}
	    else
		echo "${fl} does not contain CVS keyword"
	    fi # ! cvs2ntm
	done # end loop over ${fl}
	cvs commit -m "Convert CVS keywords to intermediate form for SVN transition"
	cd
	/bin/rm -r -f ${drc}
    done # end loop over ${drc}
fi # !0

# Convert intermediate form keywords to closest SVN keywords
if [ '0' = '1' ]; then
    cd
    for drc in `svn ls ${SVNROOT}` ; do
#    for drc in aer; do
	echo "Processing directory ${drc}..."
	cd;svn checkout file:///home/zender/svn/trunk/${drc}
	cd ${drc}
	for fl in `ls --hide=CVSROOT`; do
	    if egrep '\$_CVS_Header\$' ${fl} > /dev/null || egrep '\$_CVS_Date\$' ${fl} > /dev/null || egrep '\$_CVS_Id\$' ${fl} > /dev/null || egrep '\$_CVS_Revision\$' ${fl} > /dev/null || egrep '\$_CVS_Name\$' ${fl} > /dev/null
	    then
		echo "${fl} contains CVS indicator, translating to closest SVN keyword"
		perl -pi -e 's/\$_CVS_Header\$/\$Id\$/g;s/\$_CVS_Date\$/\$Date\$/g;s/\$_CVS_Id\$/\$Id\$/g;s/\$_CVS_Revision\$/\$Revision\$/g;s/\$_CVS_Name\$/\$HeadURL\$/g;' ${fl}
	    else
		echo "${fl} does not contain CVS keyword indicator"
	    fi # ! ntm2svn
	done # end loop over ${fl}
	svn commit -m "Convert intermediate-form CVS keywords to closest SVN keyword"
	cd
	/bin/rm -r -f ${drc}
    done # end loop over ${drc}
fi # !0

# Replace CVS modules with SVN modules
if [ '1' = '1' ]; then
    cd
#    for drc in `cd;find . -maxdepth 1 -type d | egrep "/[a-z][A-Z]*" | cut -c 3-` ; do
    for drc in `svn ls ${SVNROOT}` ; do
#    for drc in aer; do
	if [ -d ${drc} ]; then 
	    echo "Potentially SVN-controlled module ${drc} exists as local directory...processing"
	    if [ -d ${drc}/CVS ]; then 
		printf "${drc} appears to be a CVS-controlled module. Will delete and replace with SVN...\n"
		/bin/rm -r -f ${drc}
		svn checkout ${SVNROOT}/${drc}
	    else 
		printf "${drc} does not appear to be a CVS-controlled module. Continuing...\n"
	    fi # endif
	else 
	    echo "Potentially SVN-controlled module ${drc} does not exist as local directory...continuing"
	fi # endif
    done # end loop over ${drc}
fi # !0
