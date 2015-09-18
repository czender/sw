; $Id$

; Purpose: NCL script creates bxm.ps from box model run (dead.nc)
; bxm.ncl is called by bxm_run.sh

; Usage:
; sudo scp ~/bxm/bxm.ncl dust.ess.uci.edu:/tmp/bxm.ncl

load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
  wks  = gsn_open_wks("ps","bxm")
  plot = new(5,graphic)                          ; create graphics array

  res                 = True                     ; Use plot options
  res@gsnDraw         = False                    ; do not draw plot (just yet)
  res@gsnFrame        = False                    ; do not advance frame
  res@vpHeightF       = 0.5                      ; h, w, loc of plot
  res@vpWidthF        = 2.5
  ;res@vpXF           = 0.13
  ;res@vpYF           = 0.85
  res@tiXAxisString   = ""                       ; null string till last plot
  ;res@tmXBLabelFontHeightF = 0.014	         ; resize tick labels
  res@tmYLLabelFontHeightF  = 0.04
  res@tiMainFontHeightF     = 0.05               ; resize title
  ;res@txFontHeightF        = 0.01               ; resize the text
  res@tmXBLabelsOn  = False
  res@tmXTOn        = False
  res@tmYROn        = False
  res@tmXTBorderOn  = False
  res@tmYRBorderOn  = False
  res@pmLegendDisplayMode    = "Always"
  res@pmLegendSide           = "Top"                        ; Change location of 
  res@pmLegendParallelPosF   = 1.3                          ; move units right
  res@pmLegendOrthogonalPosF = -1.2                         ; move units down
  res@pmLegendWidthF         = 0.6                          ; Change width and
  res@pmLegendHeightF        = 0.6                          ; height of legend.  
  res@lgPerimOn              = False                        ; turn off box around
  res@tiXAxisString = ""
  res@tiYAxisString = ""

  sz_nbr = 4
  res@xyLineColors   = new(sz_nbr+1,integer)
  res@xyDashPatterns = new(sz_nbr+1,integer)
  res@xyExplicitLegendLabels = new(sz_nbr+1,string)
  do sz_idx=0,sz_nbr
  	res@xyLineColors(sz_idx) = sz_nbr+1-sz_idx
  	res@xyDashPatterns(sz_idx) = 0
  	res@xyExplicitLegendLabels(sz_idx) = ""
  end do
  ;res@xyDashPatterns = (/0,0,0,0,0/)
  ;print(res@xyDashPatterns)
  ;print(res@xyLineColors)

;----------------------------------------------------------
; plot 1: Wind Speed
;----------------------------------------------------------

  a = addfile("dead.nc","r")
  t = a->time
  time_nbr = dimsizes(t)

  a1 = addfile("tvbds.nc","r")
  y1 = new((/3, time_nbr/), float)
  tmp1 = a1->wnd_znl_mdp
  tmp2 = a1->wnd_mrd_mdp
  do n=0,time_nbr-1
	y1(2,n) = sqrt(tmp1(n)*tmp1(n) + tmp2(n)*tmp2(n))
  end do
  y1(1,:) = a->wnd_rfr_mbl(:,0)
  y1(0,:) = a->wnd_rfr_thr_slt(:,0)
 
  ;res@tiYAxisString = "wnd-mdp"
  res@tiMainString = "Wind Speed (m s:S:-1:N:)"
  res@xyExplicitLegendLabels(0) = "10 m Slt. Thr."
  res@xyExplicitLegendLabels(1) = "10 m Wind"
  res@xyExplicitLegendLabels(2) = "95 m Wind"
  ;res@xyLineColors  = (/3,2,1,5,4/)
  res@xyLineColors(0) = 3
  res@xyLineColors(1) = 2
  res@xyLineColors(2) = 1
  plot(0) = gsn_xy(wks,t,y1,res)

;----------------------------------------------------------
; plot 2: Wind Friction Speed
;----------------------------------------------------------

  y2 = new((/3, time_nbr/), float)
  y2(2,:) = a->wnd_frc_mbl(:,0)*100
  y2(1,:) = a->wnd_frc_slt(:,0)*100
  y2(0,:) = a->wnd_frc_thr_slt(:,0)*100
  res@tiMainString = "Wind Friction Speed (cm s:S:-1:N:)"
  res@xyExplicitLegendLabels(0) = "u:B3:*,t:S1:s"
  res@xyExplicitLegendLabels(1) = "u:B1:*:S1:s"
  res@xyExplicitLegendLabels(2) = "u:B1:*"
  plot(1) = gsn_xy(wks,t,y2,res)

;----------------------------------------------------------
; plot 3: Dust Mobilization Flux
;----------------------------------------------------------

  y3 = new((/sz_nbr+1, time_nbr/), float)
  y3(sz_nbr,:) = a->flx_mss_vrt_dst_ttl(:,0)*1e9
  do i=0,sz_nbr-1
  	y3(sz_idx,:) = a->flx_mss_vrt_dst(:,sz_idx,0)*1e9
	res@xyExplicitLegendLabels(sz_idx) = "bin" + (sz_nbr-sz_idx)
  end do
  res@tiMainString = "Dust Mobilization Flux (:F33:m:F:g m:S:-2:N: s:S:-1:N:)"
  ;res@xyExplicitLegendLabels = (/"5-10 :F33:m:F:m","2.5-5 :F33:m:F:m","1-2.5 :F33:m:F:m","0-1 :F33:m:F:m","Total"/)
  res@xyExplicitLegendLabels(sz_nbr) = "Total"
  do sz_idx=0,sz_nbr
	res@xyLineColors(sz_idx) = sz_nbr+1-sz_idx
  end do
  plot(2) = gsn_xy(wks,t,y3,res)

;----------------------------------------------------------
; plot 4: Dry Deposition Flux
;----------------------------------------------------------

  y4 = new((/sz_nbr+1, time_nbr/), float)
  y4(sz_nbr,:) = a->flx_mss_dry_sfc_ttl(:,0)*1e9
  do sz_idx=0,sz_nbr-1
        y4(sz_idx,:) = a->flx_mss_dry_sfc(:,sz_idx,0)*1e9
  end do
  res@tiMainString = "Dry Deposition Flux (:F33:m:F:g m:S:-2:N: s:S:-1:N:)"
  plot(3) = gsn_xy(wks,t,y4,res)

;----------------------------------------------------------
; plot 5: Column Mass Path
;----------------------------------------------------------

  y5 = new((/sz_nbr+1, time_nbr/), float)
  y5(sz_nbr,:) = a->mpc_dst_ttl(:,0)*1.0e3
  do sz_idx=0,sz_nbr-1
        y5(sz_idx,:) = a->mpc_dst(:,sz_idx,0)*1.0e3
  end do
  res@tmXBLabelsOn  = True
  res@tiMainString  = "Column Mass Path (g m:S:-2:N:)"
  res@tiXAxisString = "Time (days)"
  plot(4) = gsn_xy(wks,t,y5,res)


;************************************************************
; now create panel plot with its own resources
;************************************************************

  resPanel                  = True 
  ;resPanel@gsnMaximize     = True                 ; fill up the page
  ;resPanel@gsnPanelLeft    = 0.001
  ;resPanel@gsnPanelRight   = 0.999
  ;resPanel@gsnPanelBottom  = 0.001
  ;resPanel@gsnPanelTop     = 0.999
  ;resPanel@txString        = "DEAD Box Model"
  resPanel@gsnPanelLabelBar = True
  resPanel@gsnFrame         = False                  ; do not advance panel plot
  resPanel@gsnPanelBottom   = 0.05
  gsn_panel(wks,plot,(/5,1/),resPanel)

  txres               = True
  txres@txFontHeightF = 0.01
  gsn_text_ndc(wks,"             DEAD Box Model (http:://dust.ess.uci.edu/dead)  2003-02-10 13::45::28",0.4,0.02,txres)
  frame(wks)

end
