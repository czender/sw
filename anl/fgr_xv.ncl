; $Id$ -*-f90-*-

; Purpose: Plot line figures
; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/sw/anl;ncl < fgr_xv.ncl
; cd ~/sw/anl;ncl 'dvc="x11"' fgr_xv.ncl
; cd ~/sw/anl;ncl flg_clr="False" fgr_xv.ncl
; load "/home/zender/anl/fgr_xv.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

; Initialize some fields
prg_nm="fgr_xv" ; [sng] Program name
fgr_stb="fgr_xv" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36

if(.not.isvar("dvc")) then dvc="eps" end if
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("fld_nm")) then fld_nm="frc_tpt" end if

if (flg_clr) then
   fl_sfx_clr="_clr"
else ; !flg_clr
   fl_sfx_clr="_bw"
end if; !flg_clr

; Markers: http://www.ncl.ucar.edu/Document/Graphics/Images/markers.gif
idx_mrk_ast=03 ; Asterisk
idx_mrk_crc_crs=15 ; Circle with cross
idx_mrk_crc_dot=14 ; Circle with dot
idx_mrk_crc_sld=16 ; Circle solid
idx_mrk_crc=04 ; Circle 
idx_mrk_crs=05 ; Cross
idx_mrk_dmn=09 ; Diamond
idx_mrk_dot=01 ; Dot
idx_mrk_pls=02 ; Plus
idx_mrk_sqr=06 ; Square
idx_mrk_str=12 ; Star
idx_mrk_tri_dwn=08 ; Triangle down
idx_mrk_tri_up=07 ; Triangle up

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure fgr_xv
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttl_sng="Past, PI, Present, Future Dirty Snow Forcing/Response"
x_ttl_sng="Dirty Snow Forcing [W m~S~-2~N~]"
y_ttl_sng="Temperature Response [K]"

; Ordering is LGM, PI, PD, 2050
abc_frc_dst=(/0.2015,0.0363,0.0165,0.0105/)
abc_frc_bc= (/0.0152,0.0099,0.0347,0.0220/)
abc_frc_aer=(/0.2202,0.0482,0.0578,0.0357/)
abc_nbr=dimsizes(abc_frc_dst)
ord_nbr=3
abc_2d=new((/ord_nbr,abc_nbr/),float,mss_val_flt)

data=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
data_rng=new((/ord_nbr,abc_nbr/),float,mss_val_flt)

tpt_rsp_aer=    (/ 2.750,0.0370,0.1544,0.1122/)
tpt_rsp_aer_rng=(/ 0.0, 0.0, 0.0, 0.0/)
tpt_rsp_dst=    (abc_frc_dst/abc_frc_aer)*tpt_rsp_aer
tpt_rsp_dst_rng=(/ 0.0, 0.0, 0.0, 0.0/)
tpt_rsp_bc=     (abc_frc_bc/abc_frc_aer)*tpt_rsp_aer
tpt_rsp_bc_rng= (/ 0.0, 0.0, 0.0, 0.0/)

; Legend places lines in reverse order
idx_tpt_rsp_aer=2
idx_tpt_rsp_bc=1
idx_tpt_rsp_dst=0

abc_2d(idx_tpt_rsp_dst,:)=abc_frc_dst
abc_2d(idx_tpt_rsp_bc,:)=abc_frc_bc
abc_2d(idx_tpt_rsp_aer,:)=abc_frc_aer

data_rng(idx_tpt_rsp_dst,:)=tpt_rsp_dst_rng
data_rng(idx_tpt_rsp_bc,:)=tpt_rsp_bc_rng
data_rng(idx_tpt_rsp_aer,:)=tpt_rsp_aer_rng

data(idx_tpt_rsp_dst,:)=tpt_rsp_dst
data(idx_tpt_rsp_aer,:)=tpt_rsp_aer
data(idx_tpt_rsp_bc,:)=tpt_rsp_bc

dat_sng=new((/ord_nbr/),string)
dat_sng(idx_tpt_rsp_dst)="  ~F33~D~F10~T~F21~ Dust"
dat_sng(idx_tpt_rsp_bc)="  ~F33~D~F10~T~F21~ BC"
dat_sng(idx_tpt_rsp_aer)="  ~F33~D~F10~T~F21~ Total"

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
   clr_sng(idx_tpt_rsp_aer)="black"
   clr_sng(idx_tpt_rsp_dst)="blue"
   clr_sng(idx_tpt_rsp_bc)="green"
end if; !flg_clr

; Line properties
ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_sty(idx_tpt_rsp_dst)=0
ln_sty(idx_tpt_rsp_aer)=0
ln_sty(idx_tpt_rsp_bc)=50
ln_thk=new(ord_nbr,integer)
ln_thk(:)=2

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="Markers" ; Line/marker combinations
mrk_sz=new(ord_nbr,float)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer)
mrk_typ(idx_tpt_rsp_dst)=idx_mrk_crc_sld
mrk_typ(idx_tpt_rsp_aer)=idx_mrk_crc_sld
mrk_typ(idx_tpt_rsp_bc)=idx_mrk_crc_sld

rsr=True 
;rsr@tiMainFontHeightF=.020
;rsr@tmXBMode="Explicit"
;rsr@tmXBValues=ind_xaxis
;rsr@xyDashPatterns="SolidLine"
rsr@gsnDraw=False
rsr@gsnFrame=False
rsr@tiMainString=ttl_sng
rsr@tiXAxisString=x_ttl_sng
rsr@tiYAxisFontHeightF=.020
rsr@tiYAxisOn=True
rsr@tiYAxisString=y_ttl_sng
rsr@tmXTOn=False
rsr@tmYROn=False
;rsr@trXMinF=1.0
rsr@xyDashPatterns=ln_sty
rsr@xyLineColors=clr_sng
rsr@xyLineThicknesses=ln_thk
rsr@xyMarkLineModes=mrk_ln_md
rsr@xyMarkerColors=clr_sng
rsr@xyMarkerSizeF=mrk_sz
rsr@xyMarkers=mrk_typ

; Legend
rsr@lgLabelFontHeightF=0.015 ; Legend label font height
rsr@lgPerimOn=False ; Turn off box surrounding legend
rsr@pmLegendDisplayMode="Always" ; Turn on legend
rsr@pmLegendHeightF=0.15 ; Change legend height
rsr@pmLegendOrthogonalPosF=-0.35 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.3 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng ; Legend labels

; Begin plot information
drc_out=getenv("HOME")+"/ppr_ZFM07/fgr"
fl_out=drc_out+"/"+fgr_stb+"_"+fld_nm+fl_sfx_clr
print("Print output will go to "+fl_out+"."+dvc)
   
; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out) 
plot=gsn_csm_xy(wks,abc_2d,data,rsr)

; Draw range bars
if (flg_rng) then
   rsr_rng=True ; Polyline modifications desired
   rsr_rng@gsLineThicknessF=1.5 ; Line thicker
   ; fxm: convert to NDCs
   rng_bar_hrz_lng=0.07
   rng_bar_vrt=new((/ord_nbr,abc_nbr/),graphic)
   rng_bar_hrz_top=new((/ord_nbr,abc_nbr/),graphic)
   rng_bar_hrz_btm=new((/ord_nbr,abc_nbr/),graphic)
   idx_rng=(/idx_tpt_rsp_dst,idx_tpt_rsp_aer/) ; Array of ordinates with array information
   rng_nbr=dimsizes(idx_rng)
   do rng_idx=0,rng_nbr-1 ; Loop over ordinates with error information
      ord_idx=idx_rng(rng_idx)
      do abc_idx=0,abc_nbr-1
         if (.not.ismissing(data(ord_idx,abc_idx))) then
            ; Draw vertical bar
            rng_bar_vrt(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx),abc_2d(ord_idx,abc_idx)/),(/data(ord_idx,abc_idx)+data_rng(ord_idx,abc_idx),data(ord_idx,abc_idx)-data_rng(ord_idx,abc_idx)/),rsr_rng)
            ; Draw horizontal line on top of vertical bar
            rng_bar_hrz_top(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx)-rng_bar_hrz_lng,abc_2d(ord_idx,abc_idx)+rng_bar_hrz_lng/),(/data(ord_idx,abc_idx)+data_rng(ord_idx,abc_idx),data(ord_idx,abc_idx)+data_rng(ord_idx,abc_idx)/),rsr_rng)
            ; Draw horizontal line on bottom of vertical bar
            rng_bar_hrz_btm(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx)-rng_bar_hrz_lng,abc_2d(ord_idx,abc_idx)+rng_bar_hrz_lng/),(/data(ord_idx,abc_idx)-data_rng(ord_idx,abc_idx),data(ord_idx,abc_idx)-data_rng(ord_idx,abc_idx)/),rsr_rng)
         end if
      end do; end loop over abc
   end do; end loop over rng
end if; !flg_rng

draw(plot)

frame(wks)

; Clean up for next plot
delete(abc_2d)
delete(clr_sng)
delete(dat_sng)
delete(data)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
if (flg_rng) then
   delete(rng_bar_vrt)
   delete(rng_bar_hrz_top)
   delete(rng_bar_hrz_btm)
end if; !flg_rng
delete(rsr)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure xv
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end

