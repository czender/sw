#!/bin/sh -vx

# $Header: /home/zender/cvs/linux/sw/comeau.sh,v 1.5 2002-09-03 06:46:49 zender Exp $

# Purpose: Install Comeau C++
# Usage: sudo ~/linux/sw/comeau.sh
# URL: 
# http://www.comeaucomputing.com/libcomo
# http://www.comeaucomputing.com/4.3.0/minor/linux/
# http://www.comeaucomputing.com/4.0/docs/
# Compiler:
# ftp://ftp.comeaucomputing.com/pub/redhat4301.CZender.Z
# Library:
# ftp://ftp.comeaucomputing.com/pub/libcomobeta23.Z

# Test:
# export COMOROOT=/data/zender/como4301;export PATH=${PATH}\:${COMOROOT}/bin;hash -r
# cd ~/c++;como -I${COMOROOT}/libcomo -I${COMOROOT}/libcomo/cnames bug.cc ${COMOROOT}/libcomo/libcomo.a > foo 2>&1

# First, download, build, install como compiler (hence _cmp suffix)
# At this point the compiler mode must be specified as relaxed (to build libcomo)
from_cmp=${DATA}/tmp # Where *.Z binaries were downloaded to, e.g., ${DATA}/tmp
prefix_cmp=${DATA} # Prefix for installation directory, e.g., /usr/local
stub_cmp=como4301 # Name of sub-directory
fl_cmp=redhat4301.CZender # Archive containing compiler

cat > ${HOME}/comeau.sh << EOF || exit 1
#!/bin/sh
set -vx

if [ x"${prefix_cmp}" != 'x' -a x"${stub_cmp}" != 'x' -a x"${prefix_cmp}/${stub_cmp}" != 'x/' ] ; then
    printf "Executing /bin/rm -r -f ${prefix_cmp}/${stub_cmp}...\n"
    /bin/rm -r -f ${prefix_cmp}/${stub_cmp}
else
    printf 'WARNING: Define \$\{prefix_cmp\} and \$\{stub_cmp\} or risk losing everything\n'
    exit
fi # endif
mkdir -p ${prefix_cmp}/${stub_cmp}
cd ${prefix_cmp}/${stub_cmp}
cp ${from_cmp}/${fl_cmp}.Z .
ls -la ${fl_cmp}.Z
uncompress -f ${fl_cmp}.Z
ls -la ${fl_cmp}
sum ${fl_cmp}
cpio -ivd < ${fl_cmp}
rm ${stub_cmp}.cpio
rm ${fl_cmp}
chmod 755 ${prefix_cmp}/${stub_cmp}/lib ${prefix_cmp}/${stub_cmp}/bin ${prefix_cmp}/${stub_cmp}/include
chmod 755 ${prefix_cmp}/${stub_cmp}/bin/* 
chmod 644 ${prefix_cmp}/${stub_cmp}/lib/*
chmod 644 ${prefix_cmp}/${stub_cmp}/include/*
exec ./${stub_cmp}.setup
EOF

chmod a+x ${HOME}/comeau.sh
${HOME}/comeau.sh

# Second, download, build, install library libcomo.a (hence _lbr suffix)
from_lbr=${DATA}/tmp
prefix_lbr=${DATA}/como4301
stub_lbr=libcomo
fl_lbr=libcomobeta23

cat > ${HOME}/libcomeau.sh << EOF || exit 1
#!/bin/sh
set -vx

if [ x"${prefix_lbr}" != 'x' -a x"${stub_lbr}" != 'x' -a x"${prefix_lbr}/${stub_lbr}" != 'x/' ] ; then
    printf "Executing /bin/rm -r -f ${prefix_lbr}/${stub_lbr}...\n"
    /bin/rm -r -f ${prefix_lbr}/${stub_lbr}
else
    printf 'WARNING: Define \$\{prefix_lbr\} and \$\{stub_lbr\} or risk losing everything\n'
    exit
fi # endif
mkdir -p ${prefix_lbr}/${stub_lbr}
cd ${prefix_lbr}/${stub_lbr}
cp ${from_lbr}/${fl_lbr}.Z .
ls -la ${fl_lbr}.Z
uncompress ${fl_lbr}.Z
ls -la ${fl_lbr}
cpio -ivdI${fl_lbr}  # An upper case "eye" is immediately after the d
rm ${fl_lbr}
chown -R -f zender *
chgrp -R -f cgdcsm *
make -f Makefile.como
/bin/rm -f *.o *.ii *.ti
chmod 755 ${prefix_lbr}/${stub_lbr}
chmod 644 ${prefix_lbr}/${stub_lbr}/*
chmod 755 ${prefix_lbr}/${stub_lbr}/cnames
chmod 644 ${prefix_lbr}/${stub_lbr}/cnames/*
EOF

chmod a+x ${HOME}/libcomeau.sh
${HOME}/libcomeau.sh

# Third, run como compiler setup again and tell it about libcomo
# At this point the compiler mode may be changed to strict
cd ${prefix_cmp}/${stub_cmp}
${prefix_cmp}/${stub_cmp}/${stub_cmp}.setup

# Installation should now be complete
exit 0

# Answers to script: 

# libcomo (SGI C++ library)
# /data/zender/como4301/libcomo
# /data/zender/como4301/libcomo/libcomo.a
