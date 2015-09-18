;---------------------------------------------------------------------------
; Purpose: Template ncl script to create one contour plot
;          Comment or uncomment any line to customize plot
; Usage  : ncl < this
;---------------------------------------------------------------------------

load "/opt/local/lib/ncarg/nclscripts/csm/include.ncl"
begin
res = True

;---------------------------------------------------------------------------
; set netCDF file and variable
;---------------------------------------------------------------------------
a = addfile("dstmch14_2000.nc","r")
u = a->DSTODXC(0,:,:)
;u = a->DSTODXC
;u = a->DSTODXC(0,0,:,:)
;u = u * 1e9

;---------------------------------------------------------------------------
; set Strings
;---------------------------------------------------------------------------
res@tiMainString    = "Ginoux"
res@gsnCenterString = "avg = " + decimalPlaces(avg(u), 3, True)
;res@gsnLeftString   = ""
;res@gsnRightString  = ""

;---------------------------------------------------------------------------
; set contour levels (choose one of Automatic, Manual or Explicit)
;---------------------------------------------------------------------------
;res@cnLevelSelectionMode = "AutomaticLevels"

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF       = 0.1                   ; set min contour level
res@cnMaxLevelValF       = 0.6                   ; set max contour level
res@cnLevelSpacingF      = 0.05                  ; set contour spacing

;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels             = (/0.05, 0.1, 0.2, 0.4, 0.8/)

;---------------------------------------------------------------------------
; display/hide elements
;---------------------------------------------------------------------------
;res@tmXBLabelsOn  = False
;res@tmXTLabelsOn  = False
;res@tmYLLabelsOn  = False
;res@tmYRLabelsOn  = False
;res@lbLabelBarOn  = False
;res@cnInfoLabelOn = False
;res@tmXBOn        = False
;res@tmYLOn        = False

;---------------------------------------------------------------------------
; set some other parameters
;---------------------------------------------------------------------------
res@gsnDraw              = False                 ; don't draw plot yet
res@gsnFrame             = False                 ; don't advance frame
res@mpFillOn             = False                 ; turn off continent gray
res@cnFillOn             = True                  ; turn on color fill
res@gsnSpreadColors      = True                  ; use full colormap
res@gsnSpreadColorStart  = 0                     ; start at color 0 (makes min level unfilled)
res@lbLabelFontHeightF   = .01                   ; cb label font size

;---------------------------------------------------------------------------
; plot (echo "*Font : helvetica" > ncl.res to set font)
;---------------------------------------------------------------------------
;wks  = gsn_open_wks("x11","ncl")                ; view on screen
wks  = gsn_open_wks("eps","ncl")                 ; create eps file ncl.eps
       gsn_define_colormap(wks,"gui_default")
plot = gsn_csm_contour_map(wks,u,res)   
set_fill_colors(wks,plot,res@cnMinLevelValF,res)
delete(res@gsnDraw)
delete(res@gsnFrame)
plot = gsn_csm_contour_map(wks,u,res)   

print(">>> created file ncl.eps")

end
