#!/usr/bin/perl 

# Purpose: Use regular expressions to rename files 

# Usage:
# chmod a+x ~/perl/fl_rnm.pl
# scp ~/perl/fl_rnm.pl esmf.ess.uci.edu:perl

# One mouse-able line:
# for fl in `ls DSCN*`;do /bin/mv $fl ${fl/DSCN/dscn};done
# for fl in `ls *.JPG`;do /bin/mv $fl ${fl/JPG/jpg};done
# for fl in `ls IMG*.*`;do /bin/mv $fl ${fl/IMG/img};done
# for fl in `ls PANO*.*`;do /bin/mv $fl ${fl/PANO/pnr};done
# for fl in `ls VID*.*`;do /bin/mv $fl ${fl/VID/vid};done
# for fl in `ls DSC*.*`;do /bin/mv $fl ${fl/DSC/dsc};done
# ls IMG*.JPG | perl -e 'while(<STDIN>){chop;$fl_in=$_;s/IMG(_[0-9][0-9][0-9][0-9]).JPG/img$1.jpg/g;$fl_out=$_;printf STDOUT ("%s %s\n",$fl_in,$fl_out);`/bin/mv "$fl_in" $fl_out`}'
# ls DSC*.JPG | perl -e 'while(<STDIN>){chop;$fl_in=$_;s/DSC(.[0-9][0-9][0-9][0-9][0-9]).JPG/dsc$1.jpg/g;$fl_out=$_;printf STDOUT ("%s %s\n",$fl_in,$fl_out);`/bin/mv "$fl_in" $fl_out`}'
# ls MOV*.MPG | perl -e 'while(<STDIN>){chop;$fl_in=$_;s/MOV(.[0-9][0-9][0-9][0-9]).MPG/mov$1.mpg/g;$fl_out=$_;printf STDOUT ("%s %s\n",$fl_in,$fl_out);`/bin/mv "$fl_in" $fl_out`}'
# ls MVI*.AVI | perl -e 'while(<STDIN>){chop;$fl_in=$_;s/MVI(.[0-9][0-9][0-9][0-9]).AVI/mvi$1.avi/g;$fl_out=$_;printf STDOUT ("%s %s\n",$fl_in,$fl_out);`/bin/mv "$fl_in" $fl_out`}'
# chmod 644 *.jpg *.mpg

# Production:
# ls *unrolled* *crr_ncd_nvr_BFL* *crr_ncd*mlt* | ${HOME}/perl/fl_rnm.pl | m
# ls *.pdf *.ps | ${HOME}/perl/fl_rnm.pl | m
# ls avhrr_??_*.nc | ${HOME}/perl/fl_rnm.pl | m
# ls *Chj* | ${HOME}/perl/fl_rnm.pl | m
# ls 9?.??.?? | ${HOME}/perl/fl_rnm.pl | m
# ls liq_*.dat | ${HOME}/perl/fl_rnm.pl | m
# ls *.F | ${HOME}/perl/fl_rnm.pl | m
# ls *44m45* | ${HOME}/perl/fl_rnm.pl | m
# ls *xy_rgn* | ${HOME}/perl/fl_rnm.pl | m
# ls *.foo* | ${HOME}/perl/fl_rnm.pl | m
# ls *prspr* | ${HOME}/perl/fl_rnm.pl | m
# ls *8589_01* | ${HOME}/perl/fl_rnm.pl | m
# ls *lqd_mie* | ${HOME}/perl/fl_rnm.pl | m
# ls *Ch* | ${HOME}/perl/fl_rnm.pl | m
# ls *iamas* | ${HOME}/perl/fl_rnm.pl | m
# ls *SLR.JPG | ${HOME}/perl/fl_rnm.pl | m 
# ls DSC*.JPG | ${HOME}/perl/fl_rnm.pl | m 
# ls IMG*.JPG | ${HOME}/perl/fl_rnm.pl | m 
# ls MOV*.MPG | ${HOME}/perl/fl_rnm.pl | m
# ls *.JPG | ${HOME}/perl/fl_rnm.pl | m
# ls *.ogg | ${HOME}/perl/fl_rnm.pl | m
# ls *_1* | ${HOME}/perl/fl_rnm.pl | m
# ls *\'* | ${HOME}/perl/fl_rnm.pl | m
# ls *Ray_Charles-RAY_CHARLES* | ${HOME}/perl/fl_rnm.pl | m

while(<STDIN>){			 
    chop;
    $fl_in=$_;
#    s/_([89][0-9])_/_19$1/g;
#    s/199807/199807_Brb/g;
#    s/(.*).F/$1.F90/g;
#    s/(.*).F/$1.F90/g;
#    s/_xy_rgn_/_/g;
#    s/_([0-9][0-9])([0-9][0-9])_/_19{$1}19{$2}/g;
#    s/Chj/Jej/g;
#    s/lqd_mie/lqd_rds_swa_07_mie/g;
#    s/Ch/Ch0/g;
#    s/perugia/iamas_2007/g;
#    s/IMGP/imgp/g;
#    s/JPG/jpg/g;
#    s/DSC/dsc/g;
#    s/DSC([0-9][0-9][0-9][0-9][0-9]).JPG/dsc$1.jpg/g;
#    s/Lautaret02-08\ /lautaret_/g;
#    s/MOV([0-9][0-9][0-9][0-9][0-9]).MPG/mov$1.mpg/g;
#    s/_1//g;
#    s/\'//g;
#    s/ /_/g;
    s/([0-9][0-9]) - Lesson /Lesson_/g;
#    s/\ ([0-9][0-9][0-9])\ SLR.JPG/_$1_slr.jpg/g;
#    s/DSC0([0-9][0-9][0-9][0-9]).?/dsc_$1.jpg/g;
#    s/IMG_([0-9][0-9][0-9][0-9]).JPG/img_$1.jpg/g;
#    s/Ray_Charles-RAY_CHARLES/Ray_Charles/g;
    $fl_out=$_;
    printf STDOUT ("%s %s\n",$fl_in,$fl_out);
    `/bin/mv "$fl_in" $fl_out`; # Double-quotes required on fl_in not fl_out, why?
} # endwhile				

