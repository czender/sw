#/bin/sh

# $Id$

# Purpose: Change X settings for givre

# Usage:
# cd ~/sh;chmod a+x givre.sh;./givre.sh;cd -

# ~/sh/givre.sh twn # Change to TwinView
# ~/sh/givre.sh fpd # Change to Flat Panel Display

# 0. Command line options
x_md='twn' # [sng] X mode ('fpd' for flat panel display, 'twn' for Twinview)

if [ -n "${1}" ]; then
    x_md=${1}
fi # !$1

if [ ${x_md} = 'fpd' ]; then
    echo "Changing X mode to Flat Panel Display..."
    sudo cp ~/linux/etc/X11/xorg.conf.givre.minimal /etc/X11/xorg.conf
elif [ ${x_md} = 'twn' ]; then
    echo "Changing X mode to TwinView..."
    sudo cp ~/linux/etc/X11/xorg.conf.givre /etc/X11/xorg.conf
else
    echo "Bad option, exiting ..."
    exit
fi # endif x_md

# IITAF
printf "done\n"

