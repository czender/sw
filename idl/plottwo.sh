#!/bin/bash

# Purpose: contour plot of two different vars in two netCDF files
#
# Example:
#   plottwo.sh < plottwo.in

#---------------------------------------------------------------
# parse args
#---------------------------------------------------------------
usage="usage: $0 [-xpah] < plottwo.in"
function help ()
{
    echo " "
    echo $usage
    echo " "
    echo "  -x  send output to screen (x11)"
    echo "  -p  create postscript file"
    echo "  -a  use automatic contour levels"
    echo "  -h  display this help"
    echo " "
    exit 1
}
dev="eps"
auto=0
file1="";   file2=""; 
var1="ORO"; var2="ORO"; 
min1=0.0;   min2=0.0;   
max1=0.0;   max2=0.0;   
delta1=0.0; delta2=0.0; 

while getopts "xpah" opt; do
    case $opt in
        x  ) dev="x11" ;;
        p  ) dev="ps"  ;;
        a  ) auto=1    ;;
        h  ) help      ;;
        \? ) echo $usage
             exit 1    ;;
    esac
done
shift $(($OPTIND - 1))
if [ $# -ne 0 ]; then
    echo $usage
    exit 1
fi

#---------------------------------------------------------------
# read input
#---------------------------------------------------------------
while read line; do
    if [ ${line:0:1} != '#' ]; then break; fi
done

           line=${line%%#*};                   title=$line;  echo $title

read line; line=${line%%#*}; line=${line// /}; file1=$line;  echo $file1
read line; line=${line%%#*}; line=${line// /}; var1=$line;   echo $var1
read line; line=${line%%#*};                   title1=$line; echo $title1
read line; line=${line%%#*}; line=${line// /}; min1=$line;   echo $min1
read line; line=${line%%#*}; line=${line// /}; max1=$line;   echo $max1
read line; line=${line%%#*}; line=${line// /}; delta1=$line; echo $delta1

read line; line=${line%%#*}; line=${line// /}; file2=$line;  echo $file2
read line; line=${line%%#*}; line=${line// /}; var2=$line;   echo $var2
read line; line=${line%%#*};                   title2=$line; echo $title2
read line; line=${line%%#*}; line=${line// /}; min2=$line;   echo $min2
read line; line=${line%%#*}; line=${line// /}; max2=$line;   echo $max2
read line; line=${line%%#*}; line=${line// /}; delta2=$line; echo $delta2

#---------------------------------------------------------------
# variables
#---------------------------------------------------------------
psfile=${file1%.nc}.${var1}

#---------------------------------------------------------------
# create ncl script
#---------------------------------------------------------------
rm -f tmp.ncl

cat > tmp.ncl <<EOF
load "/opt/local/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/opt/local/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
begin

  wks                      = gsn_open_wks("${dev}","${psfile}")            
                             gsn_define_colormap(wks,"gui_default")
  autoflag                 = $auto
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@gsnDraw              = False              ; do not draw yet
  res@gsnFrame             = False              ; do notadvance frame yet
  res@gsnSpreadColors      = True               ; use full colormap
  res@mpFillOn             = False              ; turn off continent gray
  res@lbLabelFontHeightF   = .01                ; cb label font size
  if (autoflag .eq. 1) then
    res@cnLevelSelectionMode = "AutomaticLevels"  ; set automatic contour levels
  else
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  end if
  plot                     = new(2,graphic)     ; create graphic array for panel plot

  ;#---------------------------------------------------------------
  ;# plot(0)
  ;#---------------------------------------------------------------
  a1 = addfile("${file1}","r")
  u1 = a1->${var1}
  n  = dimsizes(dimsizes(u1))
  ;print (n)
  if (n .eq. 2) then
    v1 = u1
  end if
  if (n .eq. 3) then
    v1 = u1(0,:,:)
  end if
  if (n .eq. 4) then
    v1 = u1(0,0,:,:)
  end if
  res@tiMainString      = "${title1}"
  ;res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v1))
  if (autoflag .eq. 0) then
    res@cnMinLevelValF  = $min1              ; set min contour level
    res@cnMaxLevelValF  = $max1              ; set max contour level
    res@cnLevelSpacingF = $delta1            ; set contour spacing
  end if
  plot(0) = gsn_csm_contour_map(wks,v1,res)    

  ;#---------------------------------------------------------------
  ;# plot(1)
  ;#---------------------------------------------------------------
  a2 = addfile("${file2}","r")
  u2 = a2->${var2}
  n  = dimsizes(dimsizes(u2))
  ;print (n)
  if (n .eq. 2) then
    v2 = u2
  end if
  if (n .eq. 3) then
    v2 = u2(0,:,:)
  end if
  if (n .eq. 4) then
    v2 = u2(0,0,:,:)
  end if
  res@tiMainString      = "${title2}"
  ;res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v2))
  if (autoflag .eq. 0) then
    res@cnMinLevelValF  = $min2              ; set min contour level
    res@cnMaxLevelValF  = $max2              ; set max contour level
    res@cnLevelSpacingF = $delta2            ; set contour spacing
  end if
  plot(1) = gsn_csm_contour_map(wks,v2,res)    

  ;#---------------------------------------------------------------
  ;# create panel plot
  ;#---------------------------------------------------------------
  resPanel                    = True           ; panel mods desired
  resPanel@txString           = "${title}"
  resPanel@gsnMaximize        = True           ; blow up plot

  ;resPanel@gsnPanelBottom    = 0.05
  ;resPanel@pmLabelBarWidthF  = 10.6           ; control size of colorbar
  ;resPanel@pmLabelBarHeightF = 0.1
  ;resPanel@pmLabelBarOrthogonalPosF = -0.02   ; position cb wrt plot
  ;resPanel@lbAutoManage      = False          ; dont let plot manager contr

  gsn_panel(wks,plot,(/2,1/),resPanel)

end
EOF

# run ncl with tmp.ncl
/opt/local/bin/ncl < tmp.ncl

# rm -f tmp.ncl
echo ">>> created file: ${psfile}.${dev}"
