; Written by Wang, Wenshan 2018-04-17 15:29
;load "$HOME/scripts/ncl_script/wenshanw.ncl"
;*************************************************************
; Purpose: rgr using NCO and plot
; Prototype:
; Format:
; Steps:
; Notes:
; Usage: ncl ~/sw/ncl/airs.ncl
; To-do:
;*************************************************************
begin
var = "TSurfAir"

fname = "airs_out.pdf"
;wks = gsn_open_wks("x11", fname)
wks = gsn_open_wks("pdf", fname)
gsn_define_colormap(wks,"MPL_RdYlBu")
gsn_reverse_colormap(wks)

plots = new(2,graphic)

res = True

res@gsnFrame = False
res@gsnDraw = False
;res@gsnMaximize = True ; make plot large

res@mpLimitMode           = "LatLon"           ; required
res@mpMinLatF             = 33       ; min lat
res@mpMaxLatF             = 65       ; max lat
res@mpMinLonF             = -130       ; min lon
res@mpMaxLonF             = -90       ; max lon

  ;- contour
res@cnFillOn             = True               ; turns on the color
;res@cnFillMode = "RasterFill" ; Make blobby contours
res@cnFillMode = "CellFill" ; Exactly fill cells
res@cnLinesOn = False
res@cnFillDrawOrder = "Predraw"
res@cnLineLabelsOn = False
res@cnLineLabelFontHeightF = 0.010
res@cnLineLabelFontColor = "grey"
res@cnLineColor = "grey"
res@cnLineDrawOrder = "PreDraw"

res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = fspan(0,3,11)
res@cnLevels = fspan(270,294,11)
res@lbLabelBarOn        = False           ; turn off individual cb's

res@gsnLeftStringFontHeightF = 0.019
res@gsnLeftStringOrthogonalPosF = 0.10
res@gsnRightString = ""

res@mpFillOn             = True              ; turns off continent gray
res@mpFillDrawOrder = "PreDraw"

res@mpDataBaseVersion     = "LowRes"
res@mpDataSetName         = "Earth..4"

res@mpGeophysicalLineColor = "grey50"
res@mpGeophysicalLineDashPattern = 0
res@mpGeophysicalLineDashSegLenF = 0.1
res@mpGeophysicalLineThicknessF = 2.0

res@mpOutlineOn           = True
res@mpOutlineBoundarySets = "Geophysical"		; AllBoundaries, Geophysical

res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
res@mpGridLineColor = "grey50"
res@mpGridLineDashPattern = 1
res@mpGridLatSpacingF = 10.0
res@mpGridLonSpacingF = 10.0

res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@gsnAddCyclic          = False   ; don't add cyclic longitude point

;- rgr
; ncremap -P airs -v TSurfAir -g ${DATA}/grids/180x360_SCRIP.20150901.nc ${DATA}/hdf/AIRS.2014.10.01.202.L2.RetStd.v6.0.11.0.G14275134307.hdf ~/airs_out.nc
fn = "~/airs_out.nc"
fin = addfile(fn,"r")
din = fin->$var$

res@gsnLeftString = "b) E3SM (CAM) analysis grid"
plots(1) = gsn_csm_contour_map(wks,din,res)

;- raw
fn = "$DATA/hdf/AIRS.2014.10.01.202.L2.RetStd.v6.0.11.0.G14275134307.hdf"
fin = addfile(fn,"r")
lat2 = fin->Latitude
lon2 = fin->Longitude

res@trGridType            = "TriangularMesh"   ; faster graphic rendering
;res@cnFillMode = "RasterFill" ; Make blobby contours
res@cnFillMode = "CellFill" ; Exactly fill cells
res@sfYArray = lat2
res@sfXArray = lon2

din2 = fin->TSurfAir
res@gsnLeftString = "a) Native (AIRS) curvilinear grid"
plots(0) = gsn_csm_contour_map(wks,din2,res)

resp = True
resp@gsnPanelLabelBar    = True                ; add common colorbar
resp@lbLabelFontHeightF = 0.014
resp@lbTitleFontHeightF = 0.017
resp@lbTopMarginF = 0.5
resp@lbLabelOffsetF = 0.2
resp@lbTitleOffsetF = 0.2
;resp@lbTitleString = "Oct. 01 2014 Near-surface air temperature (K)"
resp@lbTitleString = "Near-surface air temperature (K)"

;gsn_panel(wks,plots(1),(/1,2/),resp)
gsn_panel(wks,plots,(/1,2/),resp)

end
exit
