#!/bin/bash

# Purpose: contour plot of one netCDF file, one to four variables
#
# Example:
#   plotfile.sh < plotfile.in

# Todo:
# 0. Code to plot different dimension vars
# 1. Option for auto contour levels
# 2. Option to omit x and y axis labels depending on number of plots
# 3. Option of not printing avg
# 4. Documentation

#---------------------------------------------------------------
# parse args
#---------------------------------------------------------------
usage="usage: $0 [-xpah] < plotfile.in"
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
var1="ORO"; var2="ORO"; var3="ORO"; var4="ORO"
min1=0.0;   min2=0.0;   min3=0.0;   min4=0.0
max1=0.0;   max2=0.0;   max3=0.0;   max4=0.0
delta1=0.0; delta2=0.0; delta3=0.0; delta4=0.0

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

#function crop
#{
#  line=$*
#  line=${line%%#*}; line=${line// /}; 
#  return $line; # can't return strings
#}

#---------------------------------------------------------------
# read input
#---------------------------------------------------------------
while read line; do
    if [ ${line:0:1} != '#' ]; then break; fi
done
declare -i nplots

           line=${line%%#*}; line=${line// /}; nplots=$line; echo $nplots
read line; line=${line%%#*};                   title=$line;  echo $title
read line; line=${line%%#*}; line=${line// /}; file=$line;   echo $file
read line; line=${line%%#*}; line=${line// /}; var1=$line;   echo $var1
read line; line=${line%%#*};                   title1=$line; echo $title1
read line; line=${line%%#*}; line=${line// /}; min1=$line;   echo $min1
read line; line=${line%%#*}; line=${line// /}; max1=$line;   echo $max1
read line; line=${line%%#*}; line=${line// /}; delta1=$line; echo $delta1

if [ $nplots -gt 1 ]; then
read line; line=${line%%#*}; line=${line// /}; var2=$line;   echo $var2
read line; line=${line%%#*};                   title2=$line; echo $title2
read line; line=${line%%#*}; line=${line// /}; min2=$line;   echo $min2
read line; line=${line%%#*}; line=${line// /}; max2=$line;   echo $max2
read line; line=${line%%#*}; line=${line// /}; delta2=$line; echo $delta2
fi

if [ $nplots -gt 2 ]; then
read line; line=${line%%#*}; line=${line// /}; var3=$line;   echo $var3
read line; line=${line%%#*};                   title3=$line; echo $title3
read line; line=${line%%#*}; line=${line// /}; min3=$line;   echo $min3
read line; line=${line%%#*}; line=${line// /}; max3=$line;   echo $max3
read line; line=${line%%#*}; line=${line// /}; delta3=$line; echo $delta3
fi

if [ $nplots -gt 3 ]; then
read line; line=${line%%#*}; line=${line// /}; var4=$line;   echo $var4
read line; line=${line%%#*};                   title4=$line; echo $title4
read line; line=${line%%#*}; line=${line// /}; min4=$line;   echo $min4
read line; line=${line%%#*}; line=${line// /}; max4=$line;   echo $max4
read line; line=${line%%#*}; line=${line// /}; delta4=$line; echo $delta4
fi


if [ $nplots -lt 1 ]; then
    echo "$0 error: nplots < 1"
    exit
fi
if [ $nplots -gt 4 ]; then
    echo "$0 error: nplots > 4"
    exit
fi

#---------------------------------------------------------------
# variables
#---------------------------------------------------------------
psfile=${file%.nc}.${var1}.${nplots}plot

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
  nplots                   = $nplots
  plot                     = new($nplots,graphic)     ; create graphic array for panel plot
  a                        = addfile("${file}","r")

  ;#---------------------------------------------------------------
  ;# plot(0)
  ;#---------------------------------------------------------------
  u1 = a->${var1}
  n = dimsizes(dimsizes(u1))
  ;print (n)
  if (n .eq. 2) then
    v = u1
  end if
  if (n .eq. 3) then
    v = u1(0,:,:)
  end if
  if (n .eq. 4) then
    v = u1(0,0,:,:)
  end if
  res@tiMainString      = "${title1}"
  res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v))
  if (autoflag .eq. 0) then
    res@cnMinLevelValF  = $min1              ; set min contour level
    res@cnMaxLevelValF  = $max1              ; set max contour level
    res@cnLevelSpacingF = $delta1            ; set contour spacing
  end if
  plot(0) = gsn_csm_contour_map(wks,v,res)    

  ;#---------------------------------------------------------------
  ;# plot(1)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 1) then
    u2 = a->${var2}
    n = dimsizes(dimsizes(u2))
    ;print (n)
    if (n .eq. 2) then
      v = u2
    end if
    if (n .eq. 3) then
      v = u2(0,:,:)
    end if
    if (n .eq. 4) then
      v = u2(0,0,:,:)
    end if
    res@tiMainString      = "${title2}"
    res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v))
    if (autoflag .eq. 0) then
      res@cnMinLevelValF  = $min2              ; set min contour level
      res@cnMaxLevelValF  = $max2              ; set max contour level
      res@cnLevelSpacingF = $delta2            ; set contour spacing
    end if
    plot(1) = gsn_csm_contour_map(wks,v,res)    
  end if

  ;#---------------------------------------------------------------
  ;# plot(2)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 2) then
    u3 = a->${var3}
    n = dimsizes(dimsizes(u3))
    ;print (n)
    if (n .eq. 2) then
      v = u3
    end if
    if (n .eq. 3) then
      v = u3(0,:,:)
    end if
    if (n .eq. 4) then
      v = u3(0,0,:,:)
    end if
    res@tiMainString      = "${title3}"
    res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v))
    if (autoflag .eq. 0) then
      res@cnMinLevelValF  = $min3              ; set min contour level
      res@cnMaxLevelValF  = $max3              ; set max contour level
      res@cnLevelSpacingF = $delta3            ; set contour spacing
    end if
    plot(2) = gsn_csm_contour_map(wks,v,res)    
  end if

  ;#---------------------------------------------------------------
  ;# plot(3)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 3) then
    u4 = a->${var4}
    n = dimsizes(dimsizes(u4))
    ;print (n)
    if (n .eq. 2) then
      v = u4
    end if
    if (n .eq. 3) then
      v = u4(0,:,:)
    end if
    if (n .eq. 4) then
      v = u4(0,0,:,:)
    end if
    res@tiMainString      = "${title4}"
    res@gsnCenterString   = "avg = " + sprintf("%.3g", avg(v))
    if (autoflag .eq. 0) then
      res@cnMinLevelValF  = $min4              ; set min contour level
      res@cnMaxLevelValF  = $max4              ; set max contour level
      res@cnLevelSpacingF = $delta4            ; set contour spacing
    end if
    plot(3) = gsn_csm_contour_map(wks,v,res)    
  end if

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


  if ( nplots .lt. 4) then
    gsn_panel(wks,plot,(/$nplots,1/),resPanel)
  else
    gsn_panel(wks,plot,(/2,2/),resPanel)
  end if

end
EOF

# run ncl with tmp.ncl
/opt/local/bin/ncl < tmp.ncl

# rm -f tmp.ncl
echo ">>> created file: ${psfile}.${dev}"
