begin

dev = "x11"   ; "x11", or "ps"

f1 = addfile( "/tmp/eaton/salby_GCI/GCI.nc","r")

tbr = create "sf" scalarFieldClass noparent
	"sfDataArray" : f1->tbr(0,:,:)
	"sfXCStartV" : f1->lon(0)
	"sfXCEndV" : f1->lon(dimsizes(f1->lon)-1)
	"sfYCStartV" : f1->lat(0)
	"sfYCEndV" : f1->lat(dimsizes(f1->lat)-1)
end create

; get color table
ctdir = "/fs/cgd/home0/eaton/Gist/"
n_elem = asciiread( ctdir+"tbr_240-300.gp", 1, "integer" )
ctin = asciiread( ctdir+"tbr_240-300.gp", 2+n_elem*3, "integer" )
ct = new( (/n_elem+2,3/), float )
ct(2:,:) = onedtond( ctin(2:)/255., (/n_elem,3/) )

if (dev .eq. "x11") then
	wks = create "wks" xWorkstationClass noparent
		"wkColorMap" : ct
	end create
else
 if (dev .eq. "ps") then
	wks = create "wks" psWorkstationClass defaultapp
		"wkColorMap" : ct
        end create
 end if
end if

vp = (/ .1, .95, .375, .8 /)
wn = (/ 0., 360., -90., 90. /)
vplb = (/ .2, .85, .22, .29 /)
wnlb = (/ 300., 240., 0., 1. /)

mp = create "mp" mapPlotClass wks
	"vpXF" : vp(0)
	"vpYF" : vp(3)
	"vpWidthF" : vp(1)-vp(0)
	"vpHeightF" : vp(3)-vp(2)
	"mpGridAndLimbOn" : "False"
	"mpLabelsOn" : "False"
        "mpLimitMode" : "LatLon"
        "mpMinLonF" :  wn(0)
        "mpMaxLonF" :  wn(1)
	"mpRelativeCenterLon" : "True"
	"mpGeophysicalLineColor" : "Background"
	"mpGeophysicalLineDashSegLenF" : .3
end create

getvalues mp
        "vpXF" : vpx
        "vpYF" : vpy
        "vpWidthF" : vpw
        "vpHeightF" : vph
end getvalues

tm = create "tm" tickMarkClass wks
	"vpXF" : vpx
	"vpYF" : vpy
	"vpWidthF" : vpw
	"vpHeightF" : vph
	"tmBorderThicknessF" : 1.0
	"tmXBDataLeftF" : wn(0)
	"tmXBDataRightF" : wn(1)
        "tmYLDataBottomF" : wn(2)
        "tmYLDataTopF" : wn(3)
        "tmXBMode" : "Manual"
        "tmYLMode" : "Manual"
	"tmXBTickStartF" : wn(0)
	"tmXBTickEndF" : wn(1)
	"tmXBTickSpacingF" : 30.
	"tmXBMinorPerMajor" : 2
	"tmXBLabelStride" : 2
	"tmXBLabelFont" : 21
	"tmXBLabelFontHeightF" : .017
        "tmXBMajorLengthF" : .02
	"tmXBMajorOutwardLengthF" : 0.02
        "tmXBMinorLengthF" : .01
	"tmXBMinorOutwardLengthF" : 0.01
	"tmYLTickStartF" : wn(2)
	"tmYLTickEndF" : wn(3)
	"tmYLTickSpacingF" : 30.
	"tmYLMinorPerMajor" : 2
	"tmYLLabelFont" : 21
	"tmYLLabelFontHeightF" : .017
        "tmYLMajorLengthF" : .02
	"tmYLMajorOutwardLengthF" : 0.02
        "tmYLMinorLengthF" : .01
	"tmYLMinorOutwardLengthF" : 0.01
end create

cn = create "cn" contourPlotClass wks
	"cnScalarFieldData" : tbr
	"cnFillOn" : "True"
	"cnRasterModeOn" : "True"
	"cnLineLabelsOn" : "False"
	"cnLinesOn" : "False"
	"cnInfoLabelOn" : "False"
	"cnHighLabelsOn" : "False"
	"cnLowLabelsOn" : "False"
	"cnLevelSelectionMode" : "ManualLevels"
	"cnMinLevelValF" : 1.
	"cnMaxLevelValF" : 255.
	"cnLevelSpacingF" : 254./n_elem
end create

overlay( mp, cn )
draw( (/ mp, tm /) )

lbci = ispan( 2, n_elem+1, 1 )

lb = create "lb" labelBarClass wks
	"vpXF" : vplb(0)
	"vpYF" : vplb(3)
	"vpWidthF" : vplb(1)-vplb(0)
	"vpHeightF" : vplb(3)-vplb(2)
	"lbAutoManage" : "False"
	"lbOrientation" : "Horizontal"
	"lbBoxMinorExtentF" : 1.0
	"lbMonoFillPattern" : "True"
	"lbBoxCount" : n_elem
	"lbFillColors" : lbci
	"lbLabelsOn" : "False"
	"lbBoxLinesOn" : "False"
	"lbPerimOn" : "False"
	"lbLeftMarginF" : 0.0
	"lbRightMarginF" : 0.0
	"lbBottomMarginF" : 0.0
	"lbTopMarginF" : 0.0     
end create

tm2 = create "tm2" tickMarkClass wks
	"vpXF" : vplb(0)
	"vpYF" : vplb(3)
	"vpWidthF" : vplb(1)-vplb(0)
	"vpHeightF" : vplb(3)-vplb(2)
	"tmBorderThicknessF" : 1.0
        "tmXBMode" : "Manual"
	"tmXBDataLeftF" : wnlb(0)
	"tmXBDataRightF" : wnlb(1)
	"tmXBTickStartF" : wnlb(0)
	"tmXBTickEndF" : wnlb(1)
	"tmXBTickSpacingF" : 10.
	"tmXBMinorPerMajor" : 1
	"tmXBLabelStride" : 2
	"tmXBLabelFont" : 21
	"tmXBLabelFontHeightF" : .017
        "tmXBMajorLengthF" : .02
	"tmXBMajorOutwardLengthF" : 0.02
        "tmXBMinorLengthF" : .01
	"tmXBMinorOutwardLengthF" : 0.01
	"tmXTOn" : "False"
	"tmYLOn" : "False"
	"tmYROn" : "False"
end create

draw( (/lb,tm2/) )

;clear( wks )
frame( wks )

end







