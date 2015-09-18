# $Id$

# Purpose: Simple scripts for complex RPM upgrades

# fifo_conservative
# rpm -Uvv X*rpm
# rpm -Uvv --nodeps X*rpm

# Locations to look for recent RPMS:
#rpm_url='ftp://ftp.redhat.com/pub/rawhide/i386/RedHat/RPMS'
#rpm_url='http://rpmfind.net/linux/RPM/rawhide/1.0/i386/RedHat/RPMS'
rpm_url='ftp://rpmfind.net/linux/rawhide/1.0/i386/RedHat/RPMS'

# Check what is already installed:
# rpm -qa | grep XFree

rpm_lst='XFree86-xfs-3.3.6-8 XFree86-Mach64-3.3.6-8 XFree86-libs-3.3.6-8 XFree86-75dpi-fonts-3.3.6-8 XFree86-100dpi-fonts-3.3.6-8 XFree86-3.3.6-8 XFree86-devel-3.3.6-8'

for rpm_nm in ${rpm_lst}; do
    ncks -D 5 -O -R -l /data/zender/tmp -p ${rpm_url} ${rpm_nm}.i386.rpm
done
