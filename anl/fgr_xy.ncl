; $Id$ -*-f90-*-

; Purpose: Three panel horizontal contour plot

; Usage: 
; Input data default to, e.g.,
; ${DATA}/anl_xpt/xpt_clm.nc, ${DATA}/anl_xpt/xpt_clm.nc
; Output plots default to ${DATA}/ps/fgr_xpt_ctl_xy_varname_clr.eps

; ESS 282A
; /bin/cp /data/zender/essgcm01/essgcm01.cam2.h0.0000-09-01-00000.nc ~/foo.nc
; ncks -A -v area /data/zender/essgcm01/essgcm01.clm2.h0.0000-09-01-00000.nc ~/foo.nc
; scp greenplanet.ps.uci.edu:/data/zender/essgcm01/essgcm01.cam2.h0.0000-09-01-00000.nc ~/foo.nc
; ncks -A -v area greenplanet.ps.uci.edu:/data/zender/essgcm01/essgcm01.clm2.h0.0000-09-01-00000.nc ~/foo.nc
; cd ~/sw/anl;ncl 'fl_in_nm="/home/zender/foo.nc"' 'fld_nm="TREFHT"' 'dvc="x11"' fgr_xy.ncl

; cd ~/sw/anl;ncl 'fl_in_nm="/data/zender/anl/essgcm14_clm.nc"' 'fld_nm="TREFHT"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="essgcm15"' 'ctl="essgcm14"' 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl  'lps=False' 'map_prj="CylindricalEquidistant"' 'xpt="essgcm15"' 'ctl="essgcm14"' 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'fl_in_nm="/data/zender/anl/essgcm15_mns_essgcm14_clm.nc"' lat_min=30.0 lat_max=90.0 plr=True 'fld_nm="TREFHT"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="essgcm15"' 'ctl="essgcm14"' lat_min=30.0 lat_max=90.0 plr=True 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl  'lps=False' 'map_prj="CylindricalEquidistant"' 'xpt="essgcm15"' 'ctl="essgcm14"' 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' lat_min=20.0 lat_max=50.0 lon_min=60 lon_max=100.0 fgr_xy.ncl

; cd ~/sw/anl;ncl 'xpt="essgcm15"' 'ctl="essgcm14"' 'fld_nm="PRECT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="essgcm15"' 'ctl="essgcm14"' 'fld_nm="PRECT"' 'fl_sfx="clm"' 'dvc="x11"' scl=86400000 fgr_xy.ncl

; SNICAR
; cd ~/sw/anl;ncl 'xpt="snclgm05"' 'ctl="sncpd05"' 'fld_nm="SNOAERFRC"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="TREFHT"' 'fl_sfx="ts_ANN"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' dff=False 'fld_nm="TREFHT"' 'fl_sfx="clm"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' lat_min=30.0 lat_max=90.0 plr=True 'fld_nm="TREFHT"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNOAERFRC"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNORDS"' 'dvc="x11"' fgr_xy.ncl
; cd ~/sw/anl;ncl 'lps=False' 'map_prj="CylindricalEquidistant"' 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNOWDP"' 'dvc="x11"' lat_min=20.0 lat_max=50.0 lon_min=60 lon_max=100.0 fgr_xy.ncl
; cd ~/sw/anl;ncl dff=False plr=True lps=True 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNODSTC"' 'dvc="x11"' lat_min=0.0 fgr_xy.ncl
; cd ~/sw/anl;ncl dff=False lps=False 'map_prj="CylindricalEquidistant"' 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNODSTC"' 'dvc="x11"' lat_min=0.0 fgr_xy.ncl
; cd ~/sw/anl;ncl dff=False lps=False 'map_prj="CylindricalEquidistant"' 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="SNOBCC"' 'dvc="x11"' lat_min=0.0 fgr_xy.ncl

; cd ~/sw/anl;ncl 'lps=False' 'map_prj="CylindricalEquidistant"' 'xpt="sncpd05"' 'ctl="sncpd06"' 'fld_nm="FSNT_RF"' 'dvc="x11"' lat_min=20.0 lat_max=50.0 lon_min=60 lon_max=100.0 fgr_xy.ncl
; cd ~/sw/anl;ncl 'xpt="snclgm05"' 'ctl="snclgm06"' 'dvc="x11"' fgr_xy.ncl

; Color tables:
; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
; Label format/precision specification:
; http://www.ncl.ucar.edu/Document/Graphics/format_spec.shtml
; Map projections:
; http://ncl.ucar.edu/Applications/list_maps.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$HOME/sw/ncl/csz.ncl" ; Personal NCL constants and functions

begin

; Initialize variables
prg_nm="fgr_xy.ncl"
fgr_stb="fgr_xy" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36

if(.not.isvar("ctl")) then ctl="sncpd06" end if
if(.not.isvar("dff")) then dff=True end if
if(.not.isvar("drc_out")) then drc_out=getenv("DATA")+"/ps" end if
if(.not.isvar("dvc")) then dvc="eps" end if
if(.not.isvar("hms")) then hms="NH" end if ; Hemisphere for polar plot
if(.not.isvar("fl_in_nm")) then fl_in_flg=False else fl_in_flg=True end if
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("fl_sfx")) then fl_sfx="clm" end if
if(.not.isvar("fld_nm")) then fld_nm="DUSTOD_v" end if
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("lat_max")) then lat_max=90.0 end if
if(.not.isvar("lat_min")) then lat_min=-90.0 end if
if(.not.isvar("lon_max")) then lon_max=360.0 end if
if(.not.isvar("lon_min")) then lon_min=0.0 end if
if(.not.isvar("lps")) then lps=True end if
if(.not.isvar("map_prj")) then map_prj="Mollweide" end if
if(.not.isvar("plr")) then plr=False end if
if(.not.isvar("scl")) then scl=1.0 end if
if(.not.isvar("scl_xpt")) then scl_xpt=1.0 end if
if(.not.isvar("xpt")) then xpt="sncpd05" end if
if(.not.isvar("wht_srt")) then wht_srt=False end if
if(.not.isvar("vrt")) then vrt=True end if ; Aspect ratio of panels

; Variables that depend on command-line switches
if(.not.isvar("ttl_dlt")) then ttl_dlt="Difference" end if
if(.not.isvar("ttl_ctl")) then ttl_ctl=ctl end if
if(.not.isvar("ttl_xpt")) then ttl_xpt=xpt end if
if(.not.isvar("ttl_plt")) then ttl_plt="" end if
if(map_prj.eq."CylindricalEquidistant") then lps=False end if
if(map_prj.eq."Mollweide") then lps=True end if
   
if (flg_clr) then
   fl_sfx_clr="_clr"
else ; !flg_clr
   fl_sfx_clr="_bw"
end if; !flg_clr

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (fl_in_flg) then
   fl_xpt_nm=fl_in_nm
   fl_ctl_nm=fl_in_nm
   dff=False
else 
   fl_xpt_nm=getenv("DATA")+"/anl/"+xpt+"_"+fl_sfx+".nc"
   if(dff) then fl_ctl_nm=getenv("DATA")+"/anl/"+ctl+"_"+fl_sfx+".nc" else fl_ctl_nm=fl_xpt_nm end if
end if ; !fl_in_flg
fl_xpt=addfile(fl_xpt_nm,"r")
fl_ctl=addfile(fl_ctl_nm,"r")
print("Input data for xpt = "+fl_xpt_nm)
print("Input data for ctl = "+fl_ctl_nm)
      
; Begin Arithmetic
var_xpt=fl_xpt->$fld_nm$(:,{lat_min:lat_max},{lon_min:lon_max})
if(scl_xpt.ne.1.0) then var_xpt=scl_xpt*var_xpt end if
var_ctl=fl_ctl->$fld_nm$(:,{lat_min:lat_max},{lon_min:lon_max})
if (scl.ne.1.0) then
   var_xpt=scl*var_xpt
   var_ctl=scl*var_ctl
end if
var_shp=dimsizes(var_xpt)
area=fl_xpt->area({lat_min:lat_max},{lon_min:lon_max})
gss_wgt=area(:,0) ; Change lat weight to 1-D for wgt_areaave()
time_nbr=dimsizes(fl_xpt->time)
if (time_nbr.ge.10) then flg_stt=True else flg_stt=False end if ; Try Student t-test
printVarSummary(var_xpt)

var_dlt=var_xpt-var_ctl
copy_VarCoords(var_xpt,var_dlt) ; Copy named dimensions and coordinate variables from a to b

xpt_tmp=var_xpt(lat|:,lon|:,time|:)
ctl_tmp=var_ctl(lat|:,lon|:,time|:)
dlt_tmp=var_dlt(lat|:,lon|:,time|:)

xpt_avg=dim_avg(xpt_tmp) ; Average over rightmost dimension (time)
copy_VarCoords_1(xpt_tmp,xpt_avg) ; Copy named dimensions and coordinate variables (except for rightmost dimension) from a to b
; fxm: dim_avg_Wrap averages rightmost dimension and retains metadata

ctl_avg=dim_avg(ctl_tmp) ; Average over rightmost dimension (time)
copy_VarCoords_1(ctl_tmp,ctl_avg)
dlt_avg=dim_avg(dlt_tmp) ; Average over rightmost dimension (time)
copy_VarCoords_1(dlt_tmp,dlt_avg)

xpt_vrn=dim_variance(xpt_tmp) ; Variances over rightmost dimension (time)
copy_VarCoords_1(xpt_tmp,xpt_vrn)
ctl_vrn=dim_variance(ctl_tmp) ; Variances over rightmost dimension (time)
copy_VarCoords_1(ctl_tmp,ctl_vrn)

if (flg_stt) then
   cnf_lvl=0.05 ; Confidence level
   iflag=False ; [True/False] = Population variance = [differ/similar]
   tval_opt=False ; Return t-values as well as probabilities
   prb_stt=ttest(xpt_avg,xpt_vrn,time_nbr,ctl_avg,ctl_vrn,time_nbr,iflag,tval_opt)
   copy_VarCoords_1(xpt_avg,prb_stt)
end if; !flg_stt

sfc_avg=wgt_areaave(xpt_avg,gss_wgt,1.0,0); # var, 1-D lat wgt, 1-D lon wgt, 0
sfc_avg_xpt_sng=sprintf("%5.3f",sfc_avg);
max_xpt=max(xpt_avg)
min_xpt=min(xpt_avg)
min_max_xpt_sng="Min,Max="+sprintf("%5.3f",min_xpt)+","+sprintf("%5.3f",max_xpt)
print("xpt mean (min,max) = "+sfc_avg+" ("+min_xpt+","+max_xpt+")")

sfc_avg=wgt_areaave(ctl_avg,gss_wgt,1.0,0); # var, 1-D lat wgt, 1-D lon wgt, 0
sfc_avg_ctl_sng=sprintf("%5.3f",sfc_avg);
max_ctl=max(ctl_avg)
min_ctl=min(ctl_avg)
min_max_ctl_sng="Min,Max="+sprintf("%5.3f",min_ctl)+","+sprintf("%5.3f",max_ctl)

sfc_avg=wgt_areaave(dlt_avg,gss_wgt,1.0,0); # var, 1-D lat wgt, 1-D lon wgt, 0
sfc_avg_dlt_sng=sprintf("%5.3f",sfc_avg);
max_dlt=max(dlt_avg)
min_dlt=min(dlt_avg)
min_max_dlt_sng="Min,Max="+sprintf("%5.3f",min_dlt)+","+sprintf("%5.3f",max_dlt)
; End Arithmetic

; Begin plot information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   if(dff) then fl_out=drc_out+"/"+fgr_stb+"_"+xpt+"_"+ctl+"_clm_dff_"+fld_nm+fl_sfx_clr else fl_out=drc_out+"/"+fgr_stb+"_"+xpt+"_clm_"+fld_nm+fl_sfx_clr end if
end if
print("Output will go to "+fl_out+"."+dvc)
   
; Plot defaults that may be over-ridden for specific variables
units_sng=""
; Setting cnt_lvl_md_xpt_ctl="AutomaticLevels" can result in clean white base color
; Setting cnt_lvl_md_xpt_ctl="ManualLevels" uses full colorbar machinery
if(.not.isvar("cnt_lvl_md_xpt_ctl")) then cnt_lvl_md_xpt_ctl="AutomaticLevels" end if
; Colorbar stays within contour levels?
if(.not.isvar("cnt_xtr_xpt_ctl")) then cnt_xtr_xpt_ctl=False end if
; Plot information for ctl and xpt
if(.not.isvar("cnt_min_xpt_ctl")) then cnt_min_xpt_ctl=min((/min_xpt,min_ctl/)) end if
if(.not.isvar("cnt_max_xpt_ctl")) then cnt_max_xpt_ctl=max((/max_xpt,max_ctl/)) end if
load "$HOME/sw/ncl/fld_prp.ncl" ; Set fld_sng, units_sng, and, maybe, cnt_lvl, scl_fct
if(.not.isvar("fld_sng")) then fld_sng=fld_nm end if
if(isvar("ttl_plt")) then ttl_sng=ttl_plt end if

; Plot information for dlt
cnt_lvl_md_dlt="ManualLevels"
dlt_min=min(dlt_avg)
dlt_max=max(dlt_avg)
dlt_max_abs=fabs(dlt_max)
dlt_min_abs=fabs(dlt_min)
dlt_min_max_abs=max((/dlt_max_abs,dlt_min_abs/))
if(.not.isvar("cnt_min_dlt")) then cnt_min_dlt=-dlt_min_max_abs end if
if(.not.isvar("cnt_max_dlt")) then cnt_max_dlt=dlt_min_max_abs end if

wks=gsn_open_wks(dvc,fl_out)

if(.not.isvar("clr_map_nm_xpt_ctl")) then clr_map_nm_xpt_ctl="gui_default" end if
;if(.not.isvar("clr_map_nm_xpt_ctl")) then clr_map_nm_xpt_ctl="nrl_sirkes_nowhite" end if
;if(.not.isvar("clr_map_nm_xpt_ctl")) then clr_map_nm_xpt_ctl="prcp_3" end if
if(.not.isvar("clr_map_nm_dlt")) then clr_map_nm_dlt="nrl_sirkes" end if
;if(.not.isvar("clr_map_nm_dlt")) then clr_map_nm_dlt="posneg_2" end if
clr_nbr_map_xpt_ctl=clr_map_nm_2_clr_nbr(clr_map_nm_xpt_ctl)
clr_srt_xpt_ctl=clr_map_nm_2_clr_srt(clr_map_nm_xpt_ctl)
clr_end_xpt_ctl=clr_nbr_map_xpt_ctl-1
cnt_nbr_xpt_ctl=clr_end_xpt_ctl-clr_srt_xpt_ctl+1
;print("cnt_nbr_xpt_ctl="+cnt_nbr_xpt_ctl)
cnt_ntv_xpt_ctl=(cnt_max_xpt_ctl-cnt_min_xpt_ctl)/cnt_nbr_xpt_ctl
clr_nbr_map_dlt=clr_map_nm_2_clr_nbr(clr_map_nm_dlt)
clr_srt_dlt=clr_map_nm_2_clr_srt(clr_map_nm_dlt)
clr_end_dlt=clr_nbr_map_dlt-1
cnt_nbr_dlt=clr_end_dlt-clr_srt_dlt+1
cnt_ntv_dlt=(cnt_max_dlt-cnt_min_dlt)/cnt_nbr_dlt
  
gsn_define_colormap(wks,clr_map_nm_xpt_ctl)
if(wht_srt) then
   clr_map=gsn_retrieve_colormap(wks) ; Retrieve color map for editing (dimensioned Nx3)
   clr_map(clr_srt_xpt_ctl,:)=(/1.0,1.0,1.0/) ; Replace first color with white
   gsn_define_colormap(wks,clr_map) ; Redefine colormap to workstation
end if ; !wht_srt
if(dff) then plt_xpt_ctl=new(2,graphic) else plt_xpt_ctl=new(1,graphic) end if
                     
rsr=True ; Resources for xpt plot
rsr@gsnDraw=False
rsr@gsnFrame=False
rsr@lbLabelBarOn=False
rsr@mpGridAndLimbOn=True
rsr@mpGridLatSpacingF=30.0
rsr@mpGridLonSpacingF=30.0
rsr@mpPerimOn=True
rsr@mpEllipticalBoundary=lps
rsr@mpFillOn=False
if (.not.plr) then 
   rsr@mpMinLatF=lat_min
   rsr@mpMaxLatF=lat_max
   rsr@mpProjection=map_prj
   rsr@mpMinLonF=lon_min
   rsr@mpMaxLonF=lon_max
   if (lon_max-lon_min.ne.360) then rsr@gsnAddCyclic=False end if ; Longitude range = 360
end if; plr
rsr@cnLineLabelsOn=False
rsr@cnFillOn=True ; Fill missing value areas with backgroud color
rsr@cnLinesOn=False

rsr@gsnPaperOrientation="Portrait"
;    rsr@gsnPaperOrientation="Landscape"
rsr@gsnPaperMargin=0.25
rsr@gsnSpreadColors=True ; Choose colors between start and end color
rsr@gsnSpreadColorStart=clr_srt_xpt_ctl ; 
rsr@gsnSpreadColorEnd=clr_end_xpt_ctl ; 
rsr@cnLevelSelectionMode=cnt_lvl_md_xpt_ctl ; Levels for for xpt+ctl plots
if (cnt_lvl_md_xpt_ctl.eq."ExplicitLevels") then 
   rsr@cnLevels=cnt_lvl 
   rsr@cnInfoLabelOn=False ; Info label prints bogus interval for non-linear contours
end if; 
if (cnt_lvl_md_xpt_ctl.eq."ManualLevels") then
; ManualLevels require that user set these three contour level parameters
  rsr@cnInfoLabelOn=False ; "Contour from x to y by z" label
; http://www.ncl.ucar.edu/Document/Functions/Built-in/nice_mnmxintvl.shtml
   cnt_min_max_ntv=nice_mnmxintvl(cnt_min_xpt_ctl,cnt_max_xpt_ctl,cnt_nbr_xpt_ctl,cnt_xtr_xpt_ctl)
   rsr@cnMinLevelValF=cnt_min_max_ntv(0) ; Min level for xpt+ctl plots
   rsr@cnMaxLevelValF=cnt_min_max_ntv(1) ; Max level for xpt+ctl plots
   rsr@cnLevelSpacingF=cnt_min_max_ntv(2) ; Level interval xpt+ctl plots
;   rsr@cnMinLevelValF=cnt_min_xpt_ctl ; Min level for xpt+ctl plots
;   rsr@cnMaxLevelValF=cnt_max_xpt_ctl ; Max level for xpt+ctl plots
;   rsr@cnLevelSpacingF=cnt_ntv_xpt_ctl ; Level interval for xpt+ctl plots
end if ; endif manual levels
rsr@gsnLeftString=fld_sng
rsr@gsnCenterString="Mean="+sfc_avg_xpt_sng
rsr@gsnRightString=ttl_xpt

; Always plot experiment (on top)
if(plr) then 
   rsr@gsnPolar=hms ; Specify hemisphere
   if (hms.eq."NH") then
      rsr@mpMinLatF=lat_min ; Radial extent of polar plot (use mpMaxLatF in SH)
   else ; !NH
      rsr@mpMaxLatF=lat_max ; Radial extent of polar plot (use mpMaxLatF in SH)
   end if ; !NH
   plt_xpt_ctl(0)=gsn_csm_contour_map_polar(wks,xpt_avg,rsr)
else
   plt_xpt_ctl(0)=gsn_csm_contour_map(wks,xpt_avg,rsr)
end if

; Modify experiment resources for control plot
rsr@gsnLeftString=fld_sng
rsr@gsnCenterString="Mean="+sfc_avg_ctl_sng
rsr@gsnRightString=ttl_ctl

; Plot control
if(dff) then plt_xpt_ctl(1)=gsn_csm_contour_map(wks,ctl_avg,rsr) end if

; Resources to control (common) label bar for xpt (and ctl)
rsr_lbl=True
rsr_lbl@gsnPanelLabelBar=True       ; common label bar
rsr_lbl@gsnFrame=False      ; do not advance frame yet
rsr_lbl@lbLabelAutoStride=True      
rsr_lbl@lbTitleOn=True
;rsr_lbl@lbTitleJust="CenterCenter" ; Default
rsr_lbl@lbTitleJust="TopCenter" ;
rsr_lbl@lbTitleString=units_sng
rsr_lbl@lbTitleFontHeightF=0.012
rsr_lbl@lbTitleFontHeightF=0.015
;rsr_lbl@lbAutoManage=True ; Turning this on, though set by default, font-size weirdness
rsr_lbl@lbTitleDirection="Across"
;    rsr_lbl@gsnPanelFigureStrings=(/plt_xpt_ttl,plt_ctl_ttl/)
rsr_lbl@gsnPanelBottom=0.4   ; move bottom up from 0.0 to 0.4
rsr_lbl@txString=ttl_sng
rsr_lbl@gsnMaximize=True

; Add label bar that extends across appropriate number of panels
if(dff) then 
   if(vrt) then
      pnl_arr=(/2,1/)
      rsr_lbl@lbOrientation="vertical" ; Vertical label bar
      rsr_lbl@lbTitlePosition="Bottom"
   else
      pnl_arr=(/1,2/)
      rsr_lbl@lbOrientation="horizontal" ; Horizontal label bar
      rsr_lbl@lbTitlePosition="Left"
   end if ; !vrt
else ; !dff
   pnl_arr=(/1,1/)
   rsr_lbl@lbOrientation="vertical" ; Vertical label bar
   rsr_lbl@lbTitlePosition="Bottom"
end if ; !dff
gsn_panel(wks,plt_xpt_ctl,pnl_arr,rsr_lbl)

; Modify resources to plot difference
; NB: None of this needs to be done if !dff
if(dff) then
   rsr_dlt=rsr ; Copy resources
   rsr_dlt@gsnSpreadColors=True        ; Spread colors for difference plot
   rsr_dlt@gsnSpreadColorStart=clr_srt_dlt ; 
   rsr_dlt@gsnSpreadColorEnd=clr_end_dlt ;
   rsr_dlt@cnLevelSelectionMode=cnt_lvl_md_dlt ; Level mode for difference plot
   if (cnt_lvl_md_dlt.eq."ExplicitLevels") then rsr_dlt@cnLevels=cnt_lvl end if ; Levels for difference plot
   if (cnt_lvl_md_dlt.eq."ManualLevels") then
      ; ManualLevels require that user set three contour level parameters
      ; http://www.ncl.ucar.edu/Document/Functions/Built-in/nice_mnmxintvl.shtml
      cnt_min_max_ntv=nice_mnmxintvl(cnt_min_dlt,cnt_max_dlt,cnt_nbr_dlt,False)
      rsr_dlt@cnMinLevelValF=cnt_min_max_ntv(0) ; Min level for difference plot
      rsr_dlt@cnMaxLevelValF=cnt_min_max_ntv(1) ; Max level for difference plot
      rsr_dlt@cnLevelSpacingF=cnt_min_max_ntv(2) ; Level interval for difference plot
   end if; endif manual levels
   rsr_dlt@lbLabelBarOn=False
   rsr_dlt@gsnLeftString=fld_sng
   rsr_dlt@gsnCenterString="Mean="+sfc_avg_dlt_sng
   rsr_dlt@gsnRightString=ttl_dlt
   
   ; Plot Difference
   gsn_define_colormap(wks,clr_map_nm_dlt)
   plt_dlt=gsn_csm_contour_map(wks,dlt_avg,rsr_dlt)
         
   rsr_lbl_dlt=True
   rsr_lbl_dlt@gsnPanelLabelBar=True ; Common label bar
   if(vrt) then
      rsr_lbl_dlt@gsnPanelTop=0.4 ; Upper boundary of lower panel
      rsr_lbl_dlt@gsnPanelBottom=0.1 ; Lower boundary of lower panel
   else
      rsr_lbl_dlt@gsnPanelTop=0.4 ; Upper boundary of lower panel
      rsr_lbl_dlt@gsnPanelBottom=0.0 ; Lower boundary of lower panel
   end if; !vrt
   rsr_lbl_dlt@gsnFrame=False ; Advance frame
   rsr_lbl_dlt@lbOrientation="vertical" ; Vertical label bar
   rsr_lbl_dlt@lbLabelAutoStride=True ; Nice label bar labels
   rsr_lbl_dlt@lbTitleOn=True
   ;rsr_lbl@lbTitleJust="CenterCenter" ; Default
   rsr_lbl@lbTitleJust="TopCenter" ;
   rsr_lbl_dlt@lbTitleString=units_sng
   rsr_lbl_dlt@lbTitleFontHeightF=0.012
   rsr_lbl_dlt@lbTitleDirection="Across"
   rsr_lbl_dlt@lbTitlePosition="Bottom"
   ;  rsr_lbl_dlt@gsnPanelFigureStrings=(/plt_dlt_ttl/)
   rsr_lbl_dlt@gsnMaximize=True
  
   gsn_panel(wks,plt_dlt,(/1,1/),rsr_lbl_dlt)
end if ; dff

if(flg_stt) then
   ; Resources for statistics
   rsr_prb=True
   rsr_prb@gsnDraw=False ; Do not draw plot
   rsr_prb@gsnFrame=False ; Do not advance frame
   rsr@cnLevelSelectionMode="AutomaticLevels"
   ;    rsr_prb@cnLevelSelectionMode= "ManualLevels" 
   ;rsr_prb@cnMinLevelValF=0.00 ; Min contour level
   ;rsr_prb@cnMaxLevelValF=1.05 ; Max contour level
   ;rsr_prb@cnLevelSpacingF=0.05 ; Contour spacing
   rsr_prb@cnInfoLabelOn=False ; Turn on info label
   rsr_prb@cnLinesOn=False ; Draw contour lines
   rsr_prb@cnLineLabelsOn=False ; Draw contour labels
   rsr_prb@cnLineColor="black"
   rsr_prb@cnFillScaleF=0.3 ; add extra density
  
   plt_prb=gsn_csm_contour(wks,prb_stt(:,:),rsr_prb)
   ; http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_contour_shade.shtml
   rsr_opt=True
   rsr_opt@gsnShadeFillType="color" ; "color" (default), 
   rsr_opt@gsnShadeLow=14 ; Use color index 14
   rsr_opt@gsnShadeHigh="red"   
   plt_prb=gsn_contour_shade(plt_prb,0.05,0.95,rsr_opt)
;   plt_prb=ShadeLtContour(plt_prb,cnf_lvl,17) ; Deprecated, never worked
   overlay(plt_dlt,plt_prb) ; 
end if ; !flg_stt
  
; Add label bar to difference plot if necessary
if(dff) then gsn_panel(wks,plt_dlt,(/1,1/),rsr_lbl_dlt) end if
frame(wks)

delete(var_xpt)
delete(var_ctl)
delete(var_dlt)
delete(gss_wgt)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
