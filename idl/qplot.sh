#!/bin/bash

# Purpose: quick contour plot of variable in netCDF file
#
# Examples:
#   qplot.sh dstccm50_1994.nc DSTODXC
#   qplot.sh -m dstmch14_2000.nc DSTODXC 0.0 0.5 0.05
#   qplot.sh -t "Source = Ginoux" dst_T62.nc mbl_bsn_fct

#---------------------------------------------------------------
usage="usage: $0 [-dp] [-xyc] [-t title] [-a|m] [-h] file var [min max delta]"
function help ()
{
    echo " "
    echo $usage
    echo " "
    echo "  -d  send output to display (x11)"
    echo "  -p  create postscript file (default is eps)"
    echo "  -x  omit x-axis label"
    echo "  -y  omit y-axis label"
    echo "  -c  omit colorbar"
    echo "  -t  title"
    echo "  -a  use automatic contour levels (default)"
    echo "  -m  use manual contour levels"
    echo "  -h  display this help"
    echo " "
    exit 1
}
tmpfile=/usr/tmp/tmp.ncl
dev="x11"
xaxis=True
yaxis=True
cbar=True
levmode=AutomaticLevels
title=""
min=0.0
max=0.0
delta=0.0
while getopts "dpxycamt:h" opt; do
    case $opt in
        d  ) dev="x11"     ;;
        p  ) dev="ps"      ;;
        x  ) xaxis=False   ;;
        y  ) yaxis=False   ;;
        c  ) cbar=False    ;;
        a  ) levmode=AutomaticLevels ;;
        m  ) levmode=ManualLevels    ;;
        t  ) title=$OPTARG ;;
        h  ) help          ;;
        \? ) echo $usage
             exit 1 ;;
    esac
done
shift $(($OPTIND - 1))
if [ $# -lt 2 ]; then
    echo $usage
    exit 1
fi
file=$1
var=$2
if [ $# -eq 5 ]; then
    min=$3
    max=$4
    delta=$5
fi
psfile=${file%.nc}.${var}
#---------------------------------------------------------------

# create ncl script
#rm -f /usr/tmp/tmp.ncl
touch $tmpfile
chmod a+w $tmpfile

cat > $tmpfile <<EOF
load "/opt/local/lib/ncarg/nclscripts/csm/include.ncl"
begin
  
  res                      = True               ; plot mods desired
  res@tiMainString         = "${title}"         ; main title
  res@cnFillOn             = True               ; turn on color fill
  res@gsnSpreadColors      = True               ; use full colormap
  res@gsnSpreadColorStart  = 0
  res@mpFillOn             = False              ; turn off continent gray
  res@cnLevelSelectionMode = "${levmode}"       ; set manual contour levels
  res@cnMinLevelValF       = $min               ; set min contour level
  res@cnMaxLevelValF       = $max               ; set max contour level
  res@cnLevelSpacingF      = $delta             ; set contour spacing
  res@tmXBLabelsOn         = $xaxis             ; do not draw bottom labels
  res@tmXBOn               = $xaxis             ; no bottom tickmarks
  res@tmYLLabelsOn         = $yaxis             ; do not draw labels
  res@tmYLOn               = $yaxis             ; no tickmarks
  res@lbLabelBarOn         = $cbar              ; turn off the label bar
  res@lbLabelFontHeightF   = .01                ; cb label font size

  a = addfile("${file}","r")
  u = a->${var}
  n = dimsizes(dimsizes(u))
  if (n .eq. 2) then
    v = u
  end if
  if (n .eq. 3) then
    v = u(0,:,:)
  end if
  if (n .eq. 4) then
    v = u(0,0,:,:)
  end if
  res@gsnCenterString = "avg = " + decimalPlaces(avg(u), 3, True)

  wks  = gsn_open_wks("${dev}","${psfile}")            
         gsn_define_colormap(wks,"gui_default")
  plot = gsn_csm_contour_map(wks,v,res)        
end
EOF

# run ncl with tmp.ncl
/opt/local/bin/ncl < $tmpfile
#rm -f $tmpfile
if [ $dev = "ps" ]; then
  echo ">>> created file: ${psfile}.${dev}"
fi
