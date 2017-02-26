; Purpose: Plot line figures, generally vertical figures from swnb2 output
; Based on ~/sw/anl/mie_xv.ncl
; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/sw/anl;ncl < swnb_vz.ncl
; cd ~/sw/anl;ncl 'dvc="x11"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'dvc="x11"' 'fld_nm="j_NO2"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'ctl_nm="skg_swnb_ctl"' 'xpt_nm="skg_swnb_xpt"' 'ttl_sng="Urban Night Sky Brightness"' 'fl_out_nm="$DATA/ps/skg_swnb_lnd_urb"' 'dvc="x11"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'ctl_nm="skg_swnb_ctl_lnd_rmt"' 'xpt_nm="skg_swnb_xpt_lnd_rmt"' 'ttl_sng="Remote Night Sky Brightness"' 'fl_out_nm="$DATA/ps/skg_swnb_lnd_rmt"' 'dvc="x11"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'ctl_nm="skg_swnb_ctl_ocn"' 'xpt_nm="skg_swnb_xpt_ocn"' 'ttl_sng="Ocean Night Sky Brightness"' 'fl_out_nm="$DATA/ps/skg_swnb_ocn_urb"' 'dvc="x11"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'ctl_nm="skg_swnb_ctl_ocn_noart"' 'xpt_nm="skg_swnb_xpt_ocn_noart"' 'ttl_sng="Remote Ocean Night Sky Brightness"' 'fl_out_nm="$DATA/ps/skg_swnb_ocn_rmt"' 'dvc="x11"' swnb_vz.ncl
; cd ~/sw/anl;ncl 'fl_out_nm="$DATA/ps/swnb_vz"' 'fld_nm="lmn_bb_aa_ndr"' 'ttl_sng="Zenith Luminance"' 'dvc="x11"' swnb_vz.ncl

; These files are loaded by default in NCL V6.2.0 and newer
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Initialize some fields
prg_nm="swnb_vz" ; [sng] Program name
fgr_stb="swnb_vz" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36
mss_val_dbl=1.0d36

if(.not.isvar("dff")) then dff=False end if
if(.not.isvar("dvc")) then dvc="eps" end if ; "ncgm", "png", "ps", or "x11"
if(.not.isvar("fl_in_nm")) then fl_in_flg=False else fl_in_flg=True end if
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("flg_log_x")) then flg_log_x=False end if
if(.not.isvar("flg_log_y")) then flg_log_y=False end if
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("fld_nm")) then fld_nm="lmn_bb_aa_ndr" end if
if(.not.isvar("crd_min")) then crd_min=1.0e3 end if ; Coordinate minimum
if(.not.isvar("crd_max")) then crd_max=1.0e5 end if ; Coordinate maximum
if(.not.isvar("abc_min")) then abc_min=10 end if ; Abscissa minimum
if(.not.isvar("abc_max")) then abc_max=1000 end if ; Abscissa maximum
if(.not.isvar("ord_min")) then ord_min=0.0 end if
if(.not.isvar("ord_max")) then ord_max=1.0e-2 end if

if (.not.isvar("flg_rng")) then ; is flg_rng on command line?
   flg_rng = False
end if

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
; Begin Figure swnb_vz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot defaults that may be over-ridden for specific variables
sym_sng=""
units_sng=""
load "$HOME/sw/ncl/fld_prp.ncl" ; Set fld_sng, units_sng, and, maybe, cnt_lvl, scl_fct
if(.not.isvar("ctl_nm")) then ctl_nm="skg_swnb_ctl" end if
if(.not.isvar("xpt_nm")) then xpt_nm="skg_swnb_xpt" end if
if(.not.isvar("fld_sng")) then fld_sng=fld_nm end if
if(.not.isvar("ttl_sng")) then ttl_sng=fld_sng end if
;if(.not.isvar("y_ttl_sng")) then y_ttl_sng="Height z [m]" end if
if(.not.isvar("y_ttl_sng")) then y_ttl_sng="Pressure [mb]" end if
if(.not.isvar("x_ttl_sng")) then x_ttl_sng=fld_sng+" "+sym_sng end if
if(units_sng.ne."") then x_ttl_sng=x_ttl_sng+" ["+units_sng+"]" end if

ord_nbr=2 ; Number of files to read
; Actual number of ordinates to plot: ord_nbr_plt <= ord_nbr
if(.not.isvar("ord_nbr_plt")) then ord_nbr_plt=ord_nbr end if
; Legend places lines in reverse order
idx_ctl=1 ; [idx] Ctl index
idx_xpt=0 ; [idx] Xpt index

dat_sng=new((/ord_nbr/),string)
dat_sng(idx_ctl)="  Present Day"
;dat_sng(idx_xpt)="  Geoengineered"
;dat_sng(idx_xpt)="  Urban Lighting"
dat_sng(idx_xpt)="  Dusty+Urbanized"

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
clr_sng(idx_ctl)="green"
clr_sng(idx_xpt)="blue"

; Line properties
; NB: 0*new(nbr,typ) behaves strangely due to _FillValue
ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_thk=4+0*new(ord_nbr,integer,"No_FillValue")

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="Lines" ; Line/marker combinations
mrk_sz=new(ord_nbr,float)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer)
mrk_typ(idx_ctl)=idx_mrk_crc_sld
mrk_typ(idx_xpt)=idx_mrk_crs

; Assemble file names
drc_in=getenv("DATA")+"/skg"
fl_nbr=ord_nbr
fl_in=new((/fl_nbr/),string)
fl_nm=new((/fl_nbr/),string)
fl_nm(idx_ctl)=drc_in+"/"+ctl_nm+".nc"
fl_nm(idx_xpt)=drc_in+"/"+xpt_nm+".nc"
do fl_idx=0,fl_nbr-1
   if (fl_in_flg) then
      fl_in(fl_idx)=fl_in_nm
      dff=False
   else ; !fl_in_flg
      fl_in(fl_idx)=fl_nm(fl_idx)
      if (dff) then fl_in(fl_idx)=drc_in+"/"+"swnb2_skg_dff"+".nc" end if
   end if; !fl_in_flg
   print("Input data file "+fl_idx+" = "+fl_in(fl_idx))
end do; end loop over fl

; Begin plot information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   drc_out=getenv("DATA")+"/ps"
   fl_out=drc_out+"/"+fgr_stb+"_"+fld_nm
   if(dff) then fl_out=fl_out+"_dff" end if
end if
print("Output to "+fl_out+"."+dvc)

; Determine grid size for data
fl_inq_ptr=addfile(fl_in(0),"r")
;prs=fl_inq_ptr->levp({crd_min:crd_max}) ; [mb] Pressure
prs=fl_inq_ptr->levp ; [mb] Pressure
abc=0.01*prs({crd_min:crd_max}); [Pa]->[mb] ; Ordinate means dependent variable not vertical axis
;abc=0.01*prs; [Pa]->[mb] ; Ordinate means dependent variable not vertical axis
abc_nbr=dimsizes(abc)
time_nbr=1
fld_typ=typeof(fl_inq_ptr->$fld_nm$)

fl_lst_ptr=addfiles(fl_in,"r") ; Note plural "addfiles()"
ListSetType(fl_lst_ptr,"join") ; "join" for ensembles, "cat" for timeseries
;data=fl_lst_ptr[:]->$fld_nm$({crd_min:crd_max}) ; Read fld from all files
data_raw=fl_lst_ptr[:]->$fld_nm$ ; Read fld from all files
data=data_raw(:,{crd_min:crd_max}) ; Hyperslab
;data_raw=fl_lst_ptr[:]->lmn_bb_aa ; Read fld from all files
;data=data_raw(:,2,{crd_min:crd_max}) ; Hyperslab
if (fld_typ.eq."double") then
   data_flt=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
   data_flt=doubletofloat(data)
   delete(data)
   data=data_flt
   delete(data_flt)
end if ; !double

; Scale data
if(isvar("scl_fct")) then data=scl_fct*data end if

rsr=True 
;rsr@tiMainFontHeightF=.020
;rsr@tmXBMode="Explicit"
;rsr@tmXBValues=ind_xaxis
rsr@gsnDraw=False
rsr@gsnFrame=False
rsr@tiMainString=ttl_sng
rsr@tiXAxisString=x_ttl_sng
rsr@tiYAxisFontHeightF=.020
rsr@tiYAxisOn=True
rsr@tiYAxisString=y_ttl_sng
rsr@tmXTOn=False
rsr@tmYROn=False
rsr@tmXBMajorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmXBMinorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmYLMajorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmYLMinorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
if (flg_log_y) then 
; Hard to get good labels with logarithmic axis
;   rsr@tmYLStyle="Log"
;   rsr@tmYLMode="Manual"
;   rsr@tmYLTickSpacingF=1.0
;   rsr@tmYLDataBottomF=min(ord)
;   rsr@tmYLDataTopF=max(ord)
end if ; end if flg_log_y
if (flg_log_y) then rsr@trYLog=True end if; set Y-Axis coordinate system logarithmic
;rsr@trXMinF=1.0
rsr@trYReverse=True
rsr@xyDashPatterns=ln_sty(0:ord_nbr_plt-1)
rsr@xyLineColors=clr_sng(0:ord_nbr_plt-1)
rsr@xyLineThicknesses=ln_thk(0:ord_nbr_plt-1)
rsr@xyMarkLineModes=mrk_ln_md(0:ord_nbr_plt-1)
rsr@xyMarkerColors=clr_sng(0:ord_nbr_plt-1)
rsr@xyMarkerSizeF=mrk_sz(0:ord_nbr_plt-1)
rsr@xyMarkers=mrk_typ(0:ord_nbr_plt-1)

; Legend
rsr@lgLabelFontHeightF=0.020 ; Legend label font height
rsr@lgPerimOn=False ; Turn off box surrounding legend
rsr@pmLegendDisplayMode="Always" ; Turn on legend
rsr@pmLegendHeightF=0.15 ; Change legend height
rsr@pmLegendOrthogonalPosF=-0.4 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.7 ; Move legend right NDC units
;rsr@pmLegendParallelPosF=0.3 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng(0:ord_nbr_plt-1) ; Legend labels

; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out)

; Workstation
; 20161208 transparency does not work!
rsr@wkBackgroundOpacityF=0.0 ; Default is 1.0 (opaque)
setvalues wks 
      "wkBackgroundColor" : "transparent" 
end setvalues 

plot=gsn_csm_xy(wks,data(0:ord_nbr_plt-1,:),abc,rsr)
draw(plot)
frame(wks) 

; Clean up for next plot
delete(fl_in)
delete(clr_sng)
delete(dat_sng)
delete(data)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure swnb_vz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end ; end swnb_vz
