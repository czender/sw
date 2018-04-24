; Written by Wang, Wenshan 2018-03-20 13:54
;load "$HOME/scripts/ncl_script/wenshanw.ncl"
;*************************************************************
; Purpose: rgr using NCO and plot
; Prototype:
; Format:
; Steps:
; Notes:
; Usage: ncl 'fn="cice_out.nc"' ~/sw/ncl/cice.ncl
; To-do:
;*************************************************************
begin
RGR = False
PLOT = True

if (RGR)
  ; ncap2 -s 'tmask=hi; where(hi >= 0) tmask=1;elsewhere tmask=0' ~/data/hdf/cice.hi.nc ~/cice_in.nc
  ; ncremap -d ~/data/hdf/cam.TS.nc -m ~/cice_map.nc --msk_src=tmask ~/cice_in.nc ~/cice_out.nc
  ; ncremap -d ~/data/hdf/cam.TS.nc -m ~/cice_map.nc --msk_src=tmask ~/data/hdf/cice.hi.nc ~/cice_out.nc # borken (needs tmask)
end if ; rgr

if (PLOT)
  if (.not.isvar("fn"))
    print("Error: needs input file name")
;    print('eg, ncl fn="cice_out.nc" plot-map_rgr_nco.ncl')
    print("Exit")
    exit
  end if
  if (isfilepresent(fn))
    print("Input file = "+fn)
  else
    print("Input file: "+fn+" does NOT exist!")
    print("Exit")
    exit
  end if

  plots = new(2,graphic)

  fname = "cice_out.pdf"
;  wks = gsn_open_wks("x11", fname)
  wks = gsn_open_wks("pdf", fname)
  gsn_define_colormap(wks,"MPL_RdYlBu")
  gsn_reverse_colormap(wks)

  res = True

  res@gsnFrame = False
  res@gsnDraw = False
  res@gsnMaximize = True ; make plot large

  res@gsnPolar             = "NH"               ; specify the hemisphere

  res@gsnLeftStringOrthogonalPosF = 0.10

  res@mpMinLatF            = 60.                 ; specify min lat
  res@gsnAddCyclic          = True   ; don't add cyclic longitude point
  res@mpGridLatSpacingF = 20.0
;  res@cnMissingValFillColor = "grey"

  res@mpFillOn             = False              ; turns off continent gray
  res@mpGeophysicalLineColor = "grey50"
  res@mpGeophysicalLineDashPattern = 0
  res@mpGeophysicalLineDashSegLenF = 0.1
  res@mpGeophysicalLineThicknessF = 2.0
;  res@mpOutlineDrawOrder = "PreDraw"

  res@mpDataBaseVersion     = "LowRes"
  res@mpDataSetName         = "Earth..4"
  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "Geophysical"		; AllBoundaries, Geophysical

  res@mpPerimOn         = False         ; turn off perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res@mpGridLineColor = "black"

  res@gsnPolarLabelFontHeightF = 0.017

  ;- contour
  res@cnFillOn             = True               ; turns on the color
;res@cnFillMode = "RasterFill" ; Make blobby contours
  res@cnFillMode = "CellFill" ; Exactly fill cells
  res@cnLinesOn = False
;  res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
  res@cnFillDrawOrder = "Predraw"
  res@cnLineLabelsOn = False
  res@cnLineLabelFontHeightF = 0.010
  res@cnLineLabelFontColor = "grey"
  res@cnLineColor = "grey"
  res@cnLineDrawOrder = "PreDraw"

;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF = 0
;  res@cnMaxLevelValF = 3
;  res@cnLevelSpacingF = 0.25
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = fspan(0,3,11)
  res@lbLabelBarOn        = False           ; turn off individual cb's

  fin = addfile(fn,"r")
  din = fin->hi(0,:,:)
  res@gsnLeftStringFontHeightF = 0.025
  res@gsnLeftString = "b) Atmosphere (CAM) rectangular grid"
  res@gsnRightString = ""
  plots(1) = gsn_csm_contour_map_polar(wks,din,res)

  ;- CICE
  fn = "$DATA/hdf/cice.hi.nc"
  fin = addfile(fn,"r")
  lat2 = fin->TLAT
  lon2 = fin->TLON

  res2 = True
  res2@gsnPolar             = "NH"               ; specify the hemisphere
  res2@gsnAddCyclic          = False   ; don't add cyclic longitude point
  res2@gsnFrame = False
  res2@gsnDraw = False
  res2@gsnMaximize = True ; make plot large

  res2@gsnLeftStringOrthogonalPosF = 0.10

  res2@mpFillOn             = False              ; turns off continent gray
  res2@mpMinLatF            = 60.                 ; specify min lat
  res2@mpGridLatSpacingF = 20.0
  res2@mpGeophysicalLineColor = "grey50"
  res2@mpGeophysicalLineDashPattern = 0
  res2@mpGeophysicalLineDashSegLenF = 0.1
  res2@mpGeophysicalLineThicknessF = 2.0

  res2@mpProjection          = "Satellite"   ; choose map projection

  res2@mpDataBaseVersion     = "LowRes"
  res2@mpDataSetName         = "Earth..4"
  res2@mpOutlineOn           = True
  res2@mpOutlineBoundarySets = "Geophysical"

  res2@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res2@mpPerimOn         = False         ; turn off perimeter
  res2@mpGridLineColor = "black"

  res2@trGridType            = "TriangularMesh"   ; faster graphic rendering

  ;- contour
  res2@cnFillOn             = True               ; turns on the color
;res@cnFillMode = "RasterFill" ; Make blobby contours
  res@cnFillMode = "CellFill" ; Exactly fill cells
  res2@cnLinesOn = False ; turn off contour lines
  res2@cnLineLabelsOn =  False; turn off contour line labels
  res2@cnLineColor = "grey"
;  res2@cnMissingValFillPattern = 0 ; missing value pattern is set to "SolidFill"
;  res2@cnMissingValFillColor = "grey"; color for missing values

  ;res2@cnLevelSelectionMode = "ManualLevels"
  ;res2@cnMinLevelValF = 0
  ;res2@cnMaxLevelValF = 3
  ;res2@cnLevelSpacingF = 0.25
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnLevels = fspan(0,3,11)
  res2@lbLabelBarOn        = False           ; turn off individual cb's

  res2@gsnPolarLabelFontHeightF = 0.017

  res2@sfYArray = lat2
  res2@sfXArray = lon2

  din2 = fin->hi(0,:,:)
  res2@gsnLeftStringFontHeightF = 0.025
  res2@gsnLeftString = "a) Native (CICE) curvilinear grid"
  res2@gsnRightString = ""
  plots(0) = gsn_csm_contour_map_polar(wks,din2,res2)

  resp = True
  resp@gsnPanelLabelBar    = True                ; add common colorbar
  resp@lbLabelFontHeightF = 0.012
  resp@lbTitleFontHeightF = 0.015
  resp@lbTitleString = "July 2065 Ice Thickness (m)"
  gsn_panel(wks,plots,(/1,2/),resp)

  script = "cice.ncl"
  txres = True
  txres@txFontHeightF = 0.005
  txres@txFont = 30
  gsn_text_ndc(wks,"Script: "+script, 0.1, 0.1, txres)
end if    ; plot

end
exit
