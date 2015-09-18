; $Id$ -*-f90-*-

; Purpose: Plot line figures
; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/anl;ncl < fgr_vy.ncl
; cd ~/anl;ncl 'dvc="x11"' fgr_vy.ncl
; cd ~/anl;ncl 'dvc="x11"' 'fld_nm="SNOAERFRC"' fgr_vy.ncl
; cd ~/anl;ncl dff=True 'dvc="x11"' fgr_vy.ncl
; cd ~/anl;ncl dff=True 'dvc="x11"' 'fld_nm="FSNO"' fgr_vy.ncl
; cd ~/anl;ncl dff=True 'dvc="x11"' 'fld_nm="SNORDS"' fgr_vy.ncl
; cd ~/anl;ncl dff=True 'dvc="x11"' 'fld_nm="TREFHT"' fgr_vy.ncl
; cd ~/anl;ncl flg_clr="False" fgr_vy.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Initialize some fields
prg_nm="fgr_vy" ; [sng] Program name
fgr_stb="fgr_vy" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36
mss_val_dbl=1.0d36

if(.not.isvar("dff")) then dff=False end if
if(.not.isvar("dvc")) then dvc="eps" end if
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("fld_nm")) then fld_nm="TREFHT" end if
if(.not.isvar("lon_min")) then lon_min=0.0 end if
if(.not.isvar("lon_max")) then lon_max=360.0 end if
if(.not.isvar("lat_min")) then lat_min=-90.0 end if
if(.not.isvar("lat_max")) then lat_max=90.0 end if

if (.not.isvar("flg_clr")) then ; is flg_clr on command line?
   flg_clr = True
end if
if (.not.isvar("flg_rng")) then ; is flg_rng on command line?
   flg_rng = False
end if

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
; Begin Figure fgr_vy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot defaults that may be over-ridden for specific variables
sym_sng=""
units_sng=""
load "$HOME/ncl/fld_prp.ncl" ; Set fld_sng, units_sng, and, maybe, cnt_lvl, scl_fct
if(.not.isvar("fld_sng")) then fld_sng=fld_nm end if
if(.not.isvar("ttl_sng")) then ttl_sng=fld_sng end if
if(.not.isvar("x_ttl_sng")) then x_ttl_sng=fld_sng+" "+sym_sng end if
if(.not.isvar("y_ttl_sng")) then y_ttl_sng="Latitude" end if
x_ttl_sng=x_ttl_sng+" ["+units_sng+"]"

ord_nbr=4 ; LGM, PI, PD, 2050
; Actual number of ordinates to plot: ord_nbr_plt <= ord_nbr
if(.not.isvar("ord_nbr_plt")) then ord_nbr_plt=ord_nbr end if
; Legend places lines in reverse order
idx_2050a2=0
idx_pd=1 ; [idx] PD index
idx_pi=2 ; [idx] PI index
idx_2050a2=0 ; [idx] 2050 A2 index
idx_lgm=3 ; [idx] LGM index (ord_nbr_plt discards ordinates with highest indices)
xpt_nm=new((/ord_nbr/),string)
xpt_nm(idx_2050a2)="cssnc2050_02b"
xpt_nm(idx_pd)="sncpd05"
xpt_nm(idx_pi)="cssncpi03b"
xpt_nm(idx_lgm)="snclgm05"
ctl_nm=new((/ord_nbr/),string)
ctl_nm(idx_2050a2)="cssnc2050_03b"
ctl_nm(idx_pd)="sncpd06"
ctl_nm(idx_pi)="cssncpi04b"
ctl_nm(idx_lgm)="snclgm06"

dat_sng=new((/ord_nbr/),string)
dat_sng(idx_2050a2)="  2050 IPCC A2"
dat_sng(idx_pd)="  Present"
dat_sng(idx_pi)="  Pre-industrial"
dat_sng(idx_lgm)="  Last Glacial Maximum"

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
   clr_sng(idx_lgm)="orange"
   clr_sng(idx_pi)="blue"
   clr_sng(idx_pd)="black"
   clr_sng(idx_2050a2)="green"
end if; !flg_clr

; Line properties
; NB: 0*new(nbr,typ) behaves strangely due to _FillValue
ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_thk=2+0*new(ord_nbr,integer,"No_FillValue")

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="Lines" ; Line/marker combinations
mrk_sz=new(ord_nbr,float)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer)
mrk_typ(idx_lgm)=idx_mrk_crc_sld
mrk_typ(idx_pi)=idx_mrk_crs
mrk_typ(idx_pd)=idx_mrk_crc_sld
mrk_typ(idx_2050a2)=idx_mrk_crs

; Assemble file names
drc_in=getenv("DATA")
fl_nbr=ord_nbr
fl_in=new((/fl_nbr/),string)
do fl_idx=0,fl_nbr-1
   if (dff) then fl_in(fl_idx)=drc_in+"/anl_"+xpt_nm(fl_idx)+"/"+xpt_nm(fl_idx)+"_"+ctl_nm(fl_idx)+"_clm.nc" else fl_in(fl_idx)=drc_in+"/anl_"+xpt_nm(fl_idx)+"/"+xpt_nm(fl_idx)+"_clm.nc" end if
   print("Input data file "+fl_idx+" = "+fl_in(fl_idx))
end do; end loop over fl

; Begin plot information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   drc_out=getenv("HOME")+"/ppr_ZFM07/fgr"
   fl_out=drc_out+"/"+fgr_stb+"_clm_"
   if(dff) then fl_out=fl_out+"dff_" end if
   fl_out=fl_out+fld_nm+fl_sfx_clr
end if
print("Print output will go to "+fl_out+"."+dvc)

; Determine grid size for data
fl_inq_ptr=addfile(fl_in(0),"r")
lon=fl_inq_ptr->lon({lon_min:lon_max})
lat=fl_inq_ptr->lat({lat_min:lat_max})
lon_nbr=dimsizes(lon)
lat_nbr=dimsizes(lat)
time_nbr=1
fld_typ=typeof(fl_inq_ptr->$fld_nm$)
if (fld_typ.eq."float") then
   data=new((/ord_nbr,lat_nbr,lon_nbr/),float,mss_val_flt)
else
   data=new((/ord_nbr,lat_nbr,lon_nbr/),double,mss_val_dbl)
end if
data_znl=new((/ord_nbr,lat_nbr/),float,mss_val_flt)

fl_lst_ptr=addfiles(fl_in,"r") ; Note plural "addfiles"
ListSetType(fl_lst_ptr,"join") ; "join" for ensembles, "cat" for timeseries
data=fl_lst_ptr[:]->$fld_nm$(0,{lat_min:lat_max},{lon_min:lon_max}) ; Read fld from all files
if (fld_typ.eq."double") then
   data_flt=new((/ord_nbr,lat_nbr,lon_nbr/),float,mss_val_flt)
   data_flt=doubletofloat(data)
   delete(data)
   data=data_flt
   delete(data_flt)
end if

; Scale data
if(isvar("scl_fct")) then data=scl_fct*data end if

; Zonally average data
do fl_idx=0,fl_nbr-1
    data_znl(fl_idx,:)=dim_avg(data(fl_idx,:,:)) ; Average over rightmost dimension (lon) (and retain metadata)
end do; end loop over fl
;printVarSummary(data_znl)
;print("min(data_znl)="+min(data_znl))
;print("max(data_znl)="+max(data_znl))

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
;rsr@trXMinF=1.0
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
rsr@pmLegendOrthogonalPosF=-0.7 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.35 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng(0:ord_nbr_plt-1) ; Legend labels

; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out)
plot=gsn_csm_xy(wks,data_znl(0:ord_nbr_plt-1,:),lat,rsr)
draw(plot)
frame(wks) 

; Clean up for next plot
delete(fl_in)
delete(xpt_nm)
delete(ctl_nm)
delete(clr_sng)
delete(dat_sng)
delete(data)
delete(data_znl)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure fgr_vy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end


