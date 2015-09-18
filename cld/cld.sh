#########################################################################
# RCS Identification
#########################################################################
# $Author: zender $
# $Date$
# $Id$
# $Revision$
# $Locker:  $
# $RCSfile: cld.sh,v $
# $Source: /home/zender/cvs/cld/cld.sh,v $
# $Id$
# $State: Exp $
#
# Purpose: Serves as kind of a template for getting the NQS options 
# correct and optimized for running cld().
#
# Example Usage (place mousable command lines here):
# shjob -p cld.sh -d
#
# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  1998-09-15 02:06:40  zender
# Imported sources
#
# Revision 1.2  1993/06/09  02:18:10  zender
# just wanted to check these in in their new homes once.
#
# Revision 1.1  1993/05/30  23:55:23  zender
# Initial revision
#
#########################################################################  
#
# The following lines are for the IBM cluster only:
#
# Merge stdout and stderr (-eo conflicts with -o)
#$ -eo /crestone/u1/zender/cld/cld.out
#
# Use the c-shell
#$ -s /bin/csh
#
# Name the job. -N may conflict with -eo, however.
#$ -N cld
#
# Set the queue.
#$ -G cgd_short 1 
#
# Send me mail when job is finished.
#$ -me -mu zender@ncar.ucar.edu
#
# The following lines are for the Cray batches only:
#
# Explicitly request the C-shell.
#QSUB -s /bin/csh
#
# Pick the right queue, options are: prem, reg, econ, sb (standby).
#QSUB -q reg
#
# Explicitly request the project number.
##QSUB -A 03010063
##QSUB -A 35071220
#
# Specify cpu time allotment (seconds): per process (t), whole job (T).
#QSUB -lt 3600 -lT 3600
#
# Specify memory allotment (Mw): per process (m), whole job (M).
#QSUB -lm 2Mw -lM 2Mw
#
# Send mail upon completion
#QSUB -me
#
# Redirect all output from stderr to the stdout file for the job request.
#QSUB -eo -o /crestone/u1/zender/cld/cld.out
#
# Now change to the temporary directory because although having
# room and time enough for love on crestone is very fine, every
# stderr write to the err file on crestone must go through NFS
# calls. It's much, much, much speedier to write the error file
# somewhere in /tmp and then copy it back over to crestone.
#
set timestamp
set echo
date
newacct -l
unalias rm
hostname
cd /crestone/u1/zender/cld
#
cp ccm2rad.trp35.env ${TMPDIR}/ccm2rad.trp35.env
cp ccm2rad.mls35.env ${TMPDIR}/ccm2rad.mls35.env
#cp mie.*.dat ${TMPDIR}
cd ${TMPDIR}
#
ja ${TMPDIR}/jacct 
#
########################################################################
# Perform the new standard cloud integration
########################################################################
cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50
#mswrite -t 365 nc /ZENDER/cld/H.T.final.50.nc
#mswrite -t 365 err /ZENDER/cld/H.T.final.50.err
rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.50.nc
rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.50.err
#
########################################################################
# Make the movie
########################################################################
#cld -G 8000. -g 10000. -l 104 -n 900 -D 55 -X -N 1.0e6 -x 50 -m 3.e-6 -h 4000. -H 11000. -k 4. -r 5 -p 1 -o nc -E -e err -z
#mswrite -t 365 nc /ZENDER/cld/movie.900.nc
#mswrite -t 365 err /ZENDER/cld/movie.900.err
#rcp nc flamenco.cgd.ucar.edu:/d2/zender/movie.900.nc
#rcp err flamenco.cgd.ucar.edu:/d2/zender/movie.900.err
#
########################################################################
# Modifications to standard: truncated distribution (-x 40 -m 10.e-6) or
# (-x 35 -m 20.e-6)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 35
#mswrite -t 365 nc /ZENDER/cld/H.T.final.35.nc
#mswrite -t 365 err /ZENDER/cld/H.T.final.35.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.35.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.35.err
#
########################################################################
# Modifications to standard: radiative spheres (no -X)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -x 50
#mswrite -t 365 nc /ZENDER/cld/S.T.final.50.nc
#mswrite -t 365 err /ZENDER/cld/S.T.final.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/S.T.final.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/S.T.final.50.err
#
########################################################################
# Modifications to standard: radiative spheres, truncated 
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -x 35
#mswrite -t 365 nc /ZENDER/cld/S.T.final.35.nc
#mswrite -t 365 err /ZENDER/cld/S.T.final.35.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/S.T.final.35.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/S.T.final.35.err
#
########################################################################
# Modifications to standard: 4 hour integration (-k 3. -n 4800)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 4800 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50
#mswrite -t 365 nc /ZENDER/cld/H.T.final_4hr.50.nc
#mswrite -t 365 err /ZENDER/cld/H.T.final_4hr.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final_4hr.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final_4hr.50.err
#
########################################################################
# Modifications to standard: fast updraft (-w .20)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50 -w .20
#mswrite -t 365 nc /ZENDER/cld/H.T.final.fast.50.nc
#mswrite -t 365 err /ZENDER/cld/H.T.final.fast.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.fast.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final.fast.50.err
#
########################################################################
# Modifications to standard: 4 hour truncated integration 
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 4800 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 35
#mswrite -t 365 nc /ZENDER/cld/H.T.final_4hr.35.nc
#mswrite -t 365 err /ZENDER/cld/H.T.final_4hr.35.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final_4hr.35.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.final_4hr.35.err
#
########################################################################
# Do various runs for figure 5
########################################################################
#cld -G 8000. -g 10000. -l 104 -n 7200 -D 55 -X -N 1.0e6 -x 50 -m 3.e-6 -h 4000. -H 11000. -k 8. -r 20 -p 40 -o nc -E -e err
#mswrite -t 365 nc /ZENDER/cld/H.T.eternity2.50.nc
#mswrite -t 365 err /ZENDER/cld/H.T.eternity2.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.eternity2.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.eternity2.50.err
# Note!!!!!! memory requested should be increased for this case!!!!!
#
#cld -U -s .4 -S 1.00 -G 8000. -g 10000. -l 104 -n 900 -D 55 -X -N 1.0e6 -x 50 -m 3.e-6 -h 4000. -H 11000. -k 4. -r 5 -p 10 -o nc -E -e err
#mswrite -t 365 nc /ZENDER/cld/H.T.still.50.nc
#mswrite -t 365 err /ZENDER/cld/H.T.still.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.still.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.T.still.50.err
#
#cld -M 0 -Z 0 -G 8000. -g 4000. -l 100 -n 900 -D 55 -X -x 40 -m 10.e-6 -h 2000. -H 9500. -k 4. -r 5 -p 10 -o nc -E -e err
#mswrite -t 365 nc /ZENDER/cld/H.M.heyms.40.nc
#mswrite -t 365 err /ZENDER/cld/H.M.heyms.40.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.M.heyms.40.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.M.heyms.40.err
#
#cld -M 0 -Z 1 -G 8000. -g 4000. -l 100 -n 900 -D 55 -X -N .03e6 -x 50 -m 3.e-6 -h 2000. -H 9500. -k 4. -r 5 -p 10 -o nc -E -e err
#mswrite -t 365 nc /ZENDER/cld/H.M.radke.50.nc
#mswrite -t 365  err /ZENDER/cld/H.M.radke.50.err
#rcp nc flamenco.cgd.ucar.edu:/data2/zender/cld/H.M.radke.50.nc
#rcp err flamenco.cgd.ucar.edu:/data2/zender/cld/H.M.radke.50.err
#
########################################################################
#end regular runs
########################################################################
#
#execute the job accounting
ja -chlst ${TMPDIR}/jacct
alias rm rm -i     





