; $Id$

; Purpose: Collection of routines for publication quality plots of global datasets

@~/idl/ibp_cln.pro
@~/idl/ibp_clr.pro
@~/idl/aspect.pro ; Needed for polar plot aspect ratios

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM time-latitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_ty_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

fld2_stf=[ $
	['LWCF','!5LWCF','!5W m!E-2!N'], $
	['SWCF','!5SWCF','!5W m!E-2!N'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
	['erbe_b','!5ERBE'], $
;	['422','!5CCM2'], $
	['nflux18','!5NFLUX18'], $
	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'], $
	['sld012d','!5CCM3']]
;	['spcp_85','!5ANV']]
;	['erbe_b','!5ERBE']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

; NB: Jan is mth_srt=0 in C,IDL indexing notation
mth_srt=0
mth_end=11
; NB: Store files with Fortran indexing notation
mth_srt_sng=string(format='(I2.2)',mth_srt+1)
mth_end_sng=string(format='(I2.2)',mth_end+1)

for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
	y_ttl='!5Latitude'
	x_ttl='!5Month of Year'
	pre_scl=0
	fl_clr_tlb=''
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	dsh=0
;	yr_sng='8589x87'
	yr_sng='8589'
	if fld2(fld_idx2).nm eq 'LWCF' then pre_scl=[-24,4,13]
;	if fld2(fld_idx2).nm eq 'LWCF' then pre_scl=[-20,5,9]
	if fld2(fld_idx2).nm eq 'SWCF' then pre_scl=[-100,10,21]
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_anm_xavg_'+yr_sng+'_'+mth_srt_sng+mth_end_sng+'.nc'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_anm_xavg_'+yr_sng+'_'+mth_srt_sng+mth_end_sng+'_'+fld2(fld_idx2).nm+'.eps'
	if prn then begin
		ttl='Seasonal Cycle '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		if xpt_idx eq 0 then top=1 else top=0
		if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
	endif; endif prn
	if prn then begin
		x_sz=6.5
		y_sz=4.0
		if top then y_sz=y_sz*1.08
		if btm then y_sz=y_sz*1.16
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_ty,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,pll=1,clr_ord=1,fld_nm=fld2(fld_idx2).nm,cbar=0,wrp=1,rng_y=[-90,90],pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,top=top,btm=btm,scl=scl,cnt_ovl=cnt_ovl
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over experiment
endfor; end loop over fld2

end; end gcm_ty_bch()
pro gcm_ty, $
	btm=btm, $
	cbar=cbar, $
	chr_sz=chr_sz, $
	cnt_ovl=cnt_ovl, $
	dsh=dsh, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	mth_end=mth_end, $
	mth_srt=mth_srt, $
	clr_ord=clr_ord, $
	pll=pll, $
	pre_scl=pre_scl, $
	prn=prn, $
	scl=scl, $
	top=top, $
	ttl=ttl, $
	wrp=wrp, $
	rng_x=rng_x, $
	x_ttl=x_ttl, $
	rng_y=rng_y, $
	y_ttl=y_ttl

; gcm_ty,ttl='!5ERBE Seasonal Cycle LWCF (W m!E-2!N)',pll=4,clr_ord=1

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(wrp) eq 0 then wrp=1
if n_elements(cbar) eq 0 then cbar=0
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='!5ERBE Seasonal Cycle LWCF (W m!E-2!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=22
if n_elements(fl_clr_tbl) eq 0 then fl_clr_tbl=''
if n_elements(prn) eq 0 then prn=0
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if n_elements(x_ttl) eq 0 then y_ttl='!5Latitude'
if n_elements(y_ttl) eq 0 then x_ttl='!5Month of Year'
if n_elements(pll) eq 0 then pll=1
if n_elements(clr_ord) ne 0 then clr_ord=clr_ord else clr_ord=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/erbe_b/erbe_b_anm_xavg_8589_0112.nc'
if n_elements(fld_nm) eq 0 then fld_nm='LWCF'

clr_rfr,clr_tbl_nbr,"",0
erase
!p.multi=0

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_nm)
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,foo,lat_nbr
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,fld_nm,data
ncdf_close,nc_id
; End of NetCDF commands

data=transpose(reform(data)) ; tranposing makes contour() put lat on x-axis, lev on y-axis
good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=indgen(time_nbr)
ord=lat
if wrp then begin
data=[data,data(0,*)]
abc=[abc,time_nbr]
endif

data_min=min(data(where(data lt 1.0e20)))
data_max=max(data(where(data lt 1.0e20)))

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

if auto_scl then begin

; Automatically set and scale the contour levels

cntr_lvl_nbr=20

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)
cntr_lvl_min=min(cntr_lvl)
cntr_lvl_max=max(cntr_lvl)

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl
if unit_pfx_sng ne '' then print,'cntr_lbl_mk: WARNING unit_pfx_sng == ',unit_pfx_sng

; Customize contour labels
cntr_chr_sz=0.83*chr_sz		;charsize of level labels
if pll eq 1 then begin
for lvl_idx=0,cntr_lvl_nbr-2 do begin
	if cntr_lvl(lvl_idx) lt 0.0 then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
endfor
for lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	if not dsh then cntr_ln_sty(lvl_idx)=0
	if cntr_lvl(lvl_idx) eq 0.0 then cntr_which_lbl(lvl_idx)=1
endfor
endif; endif pll

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cbar_lgn=cntr_lvl_lbl

cbar_idx=indgen(cbar_clr_nbr)+2
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=1.3
cbar_unit=''
cbar_lbl_sz=1.5

if cbar then begin

cbar_psn=[.88,0.1,0.92,0.92]; [x_min,y_min,x_max,y_max]
plt_rgn_nrm=[0.12,0.1,0.87,0.92]; [x_min,y_min,x_max,y_max]

endif else begin

plt_rgn_nrm=[0.12,0.1,0.97,0.92]; [x_min,y_min,x_max,y_max] 

endelse

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,0
tvlct,r_curr,g_curr,b_curr

if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

if strpos(fl_nm,'erbe') ne -1 then begin
; Do not plot ERBE data poleward of 60 degrees.
; Start at top left go clockwise through vertices
; Filling polygons with white does not work because the white overwrites the axes
;polyfill,[abc_min,abc_max,abc_max,abc_min],[ord_max,ord_max,60.0,60.0],/DATA,color=clr_wht_idx
;polyfill,[abc_min,abc_max,abc_max,abc_min],[ord_min,ord_min,-60.0,-60.0],/DATA,color=clr_wht_idx
good_ord=where(abs(ord) lt 60.0,nbr_ord)
ord=ord(good_ord)
nbr_abc=n_elements(abc)
tmp_data=fltarr(nbr_abc,nbr_ord)
for idx_ord=0,nbr_ord-1 do begin
	tmp_data(*,idx_ord)=data(*,good_ord(idx_ord))
endfor
data=tmp_data
endif

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
if cbar then mrg_btm=mrg_btm+4
mrg_lft=4.1 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.5 ; 2 is default
	if top then mrg_top=mrg_top+1.3
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.6
	if cbar then mrg_btm=mrg_btm+1.75
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	xstyle=5, $
	ystyle=5, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
;	position=plt_rgn_nrm, $
;	/closed, $			
	/fill, $			
	/noerase

lat_axz_drw,lat_top=!y.crange(1),lat_btm=!y.crange(0),ntt=1,ttl=y_ttl,chr_sz=chr_sz,axz_vrt=1

if btm then ntt=1 else ntt=0
time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

end; end gcm_ty()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM time-latitude commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM Hovmoller
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_xt_bch, $
	ia=ia, $
	anm=anm, $
	prn=prn
if n_elements(ia) eq 0 then ia=0
if n_elements(anm) eq 0 then anm=0
if ia then ia_sng='ia_' else ia_sng=''
if anm then anm_sng='anm_' else anm_sng=''
if n_elements(prn) eq 0 then prn=0

fld2_stf=[ $
;	['TOTCWP','!5Condensate Path','!5(g m!E-2!N)'], $
;	['CLDTOT','!5Total Cloud','!5(%)'], $
;	['CLDLOW','!5Low Cloud','!5(%)'], $
;	['CLDMED','!5Middle Cloud','!5(%)'], $
;	['CLDHGH','!5High Cloud','!5(%)']]
;	['ALBEDO','!5Albedo','!5(%)'], $
;	['LWCF','!5LWCF','!5(W m!E-2!N)']]
;	['TS1','!5SST','(!E!12_!5!NK)']]
	['FLNT','!5OLR','(!5W m!E-2!N)'] ]
;	['LWCF','!5LWCF','(!5W m!E-2!N)'], $
;	['SWCF','!5SWCF','(!5W m!E-2!N)'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
	['erbe_b','!5ERBE']]
;	['422','!5CCM2'], $
;	['nflux18','!5NFLUX18'], $
;	['amip','!5AMIP']]
;	['noaa','!5NOAA'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'], $
;	['sld012d','!5CCM3']]
;	['spcp_85','!5ANV']]
;	['noaa','!5NOAA']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpt

; NB: Jan is mth_srt=0 in C,IDL indexing notation
yr_srt=85
yr_end=89
mth_srt=0
mth_end=59
;yr_srt=79
;yr_end=93
;mth_srt=0
;mth_end=11
; NB: Store files with Fortran indexing notation
yr_srt_sng=string(format='(I2.2)',yr_srt)
yr_end_sng=string(format='(I2.2)',yr_end)
mth_srt_sng=string(format='(I2.2)',mth_srt+1)
mth_end_sng=string(format='(I2.2)',mth_end+1)
time_sng=yr_srt_sng+yr_end_sng+'_'+mth_srt_sng+mth_end_sng
;time_sng=yr_srt_sng+mth_srt_sng+'_'+yr_end_sng+mth_end_sng
;mth_srt=0
;mth_end=179
for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin

;	Annual Cycle 
;	ttl=xpt(xpt_idx).sng+' Seasonal Cycle '+fld2(fld_idx2).sng+' '+fld2(fld_idx2).unit
;	y_ttl='!5Month of Year'
;	chr_sz=1.5
;	dsh=1
;	if strpos(fld2(fld_idx2).nm,'LWCF') ne -1 then pre_scl=[-25,5,11]
;	if strpos(fld2(fld_idx2).nm,'SWCF') ne -1 then pre_scl=[-40,10,9]
;	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+anm_sng+'yavg_00N20N_8589x87_'+mth_srt_sng+mth_end_sng+'.nc'
;	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+anm_sng+'yavg_00N20N_8589x87_'+mth_srt_sng+mth_end_sng+'.eps'

;	Interannual
	ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' '+fld2(fld_idx2).unit+' 1985-1989'
	rng_y=[mth_srt,mth_end]
	rng_x=[135,270]
	tck_lcn=[45,3]
;	rng_x=[100,280]
;	tck_lcn=[20,2]
	scl=1.0
	top=1
	btm=1
	y_ttl='!5Year'
	chr_sz=0.75
	dsh=0
	yavg_sng='yavg_10S10N_'
;	yavg_sng='yavg_05S05N_'
	if anm then pll=6 else pll=4
	if strpos(fld2(fld_idx2).nm,'TS1') ne -1 then if anm then pre_scl=[-5,0.5,21] else pre_scl=[296.5,0.5,14] 
;	if strpos(fld2(fld_idx2).nm,'LWCF') ne -1 then if anm then pre_scl=[-100,10,21] else pre_scl=[0,5,15]
;	if strpos(fld2(fld_idx2).nm,'SWCF') ne -1 then if anm then pre_scl=[-100,10,21] else pre_scl=[-120,10,13]
	if strpos(fld2(fld_idx2).nm,'FLNT') ne -1 then if anm then pre_scl=[-50,10,11] else pre_scl=[180,10,11]
	if strpos(fld2(fld_idx2).nm,'LWCF') ne -1 then if anm then pre_scl=[-50,10,11] else pre_scl=[5,10,9]
	if strpos(fld2(fld_idx2).nm,'SWCF') ne -1 then if anm then pre_scl=[-50,10,11] else pre_scl=[-110,10,11]
;	if strpos(fld2(fld_idx2).nm,'TOTCWP') ne -1 then if anm then pre_scl=[-40,10,9] else pre_scl=[20,5,12]
	if strpos(fld2(fld_idx2).nm,'TOTCWP') ne -1 then if anm then pre_scl=[-35,5,15] else pre_scl=[20,5,12]
	if strpos(fld2(fld_idx2).nm,'CLDLOW') ne -1 then if anm then pre_scl=[-20,5,9]  else pre_scl=[10,5,10] 
	if strpos(fld2(fld_idx2).nm,'CLDMED') ne -1 then if anm then pre_scl=[-25,5,11] else pre_scl=[0,5,9]
	if strpos(fld2(fld_idx2).nm,'CLDHGH') ne -1 then if anm then pre_scl=[-35,5,15] else pre_scl=[5,5,19]
	if strpos(fld2(fld_idx2).nm,'CLDTOT') ne -1 then if anm then pre_scl=[-25,5,11] else pre_scl=[30,5,14]
	if strpos(fld2(fld_idx2).nm,'CLD') ne -1 then scl=100.0
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+ia_sng+anm_sng+yavg_sng+time_sng+'.nc'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+ia_sng+anm_sng+yavg_sng+time_sng+'_'+fld2(fld_idx2).nm+'.eps'
	if prn then begin
		if xpt_idx eq 0 then top=1 else top=0
		if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
	endif; endif prn and diff
	if prn then begin
		x_sz=3.2
		y_sz=8
;		ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' '+fld2(fld_idx2).unit+' 1979-1993'
;		ttl=fld2(fld_idx2).sng+' '+fld2(fld_idx2).unit+' 1985-1989'
		top=1
		btm=1
		ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' '+fld2(fld_idx2).unit+' 1985-1989'
;		x_sz=3.2
;		y_sz=2.6
		if top then y_sz=y_sz*1.05
		if btm then y_sz=y_sz*1.04
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_xt,y_ttl=y_ttl,rng_y=rng_y,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,pll=pll,order=1,fld_nm=fld2(fld_idx2).nm,cbar=0,top=top,btm=btm,wrp=1,rng_x=rng_x,pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,yr_srt=yr_srt,tck_lcn=tck_lcn,scl=scl
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over experiments
endfor; end loop over fld2s

end; end gcm_xt_bch()
pro gcm_xt, $
	fl_nm=fl_nm, $
	pll=pll, $
	order=order, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	y_ttl=y_ttl, $
	tck_lcn=tck_lcn, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	dsh=dsh, $
	wrp=wrp, $
	cbar=cbar, $
	scl=scl, $
	yr_srt=yr_srt, $
	pre_scl=pre_scl, $
	fld_nm=fld_nm

; gcm_xt,ttl='!5ERBE Seasonal Cycle LWCF (W m!E-2!N)',pll=4,order=1

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(wrp) eq 0 then wrp=0
if n_elements(cbar) eq 0 then cbar=1
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl=''
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl='!5Month of Year'
if n_elements(x_ttl) eq 0 then x_ttl='!5Longitude'
if n_elements(order) ne 0 then clr_ord=order else clr_ord=1
if n_elements(pll) eq 0 then pll=3
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/erbe_b/erbe_b_anm_yavg_00N20N_8589_0112.nc'
if n_elements(fld_nm) eq 0 then fld_nm='LWCF'
if n_elements(yr_srt) eq 0 then yr_srt=85

if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=34
color_order=clr_ord
clr_rfr,clr_tbl_nbr,"",0

erase
!p.multi=0

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_nm)
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
ncdf_varget,nc_id,'lon',lon
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,fld_nm,tmp_data
ncdf_close,nc_id
; End of NetCDF commands

tmp_data=reform(tmp_data)
data=tmp_data
data=data*scl
abc=lon
ord=indgen(time_nbr)
;ord=time

if wrp then begin
nbr_abc=n_elements(abc)
nbr_ord=n_elements(ord)
if (abs((lon(nbr_abc-1)+lon(1)) mod 180.0) le 1.0e-5) then begin
; Wrap in zonal direction
abc=[abc,360.0]
nbr_abc=n_elements(abc)
data=[data,data(0,*)]
endif; end zonal wrap
if ord(0) eq 0 and ord(0) eq 11 then begin
; Wrap in month direction for annual cycles
tmp_data=fltarr(nbr_abc,nbr_ord+1)
tmp_data(0:nbr_abc-1,1:nbr_ord)=data
tmp_data(*,0)=data(*,nbr_ord-1)
data=tmp_data
ord=[0,ord]
nbr_ord=n_elements(ord)
endif; end time wrap
endif

data_min=min(data)
data_max=max(data)

if auto_scl then begin

; Automatically set and scale the contour levels
cntr_lvl_nbr=20

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl
if unit_pfx_sng ne '' then print,'cntr_lbl_mk: WARNING unit_pfx_sng == ',unit_pfx_sng

; Customize contour labels
cntr_chr_sz=0.83*chr_sz		;charsize of level labels
if pll eq 1 then begin
for lvl_idx=0,cntr_lvl_nbr-2 do begin
	if cntr_lvl(lvl_idx) lt 0.0 then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
endfor
for lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	if not dsh then cntr_ln_sty(lvl_idx)=0
	if cntr_lvl(lvl_idx) eq 0.0 then cntr_which_lbl(lvl_idx)=1
endfor
endif; endif pll

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cbar_idx=indgen(cbar_clr_nbr)+2
cbar_lgn=cntr_lvl_lbl

cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.86*chr_sz
cbar_unit=''
cbar_lbl_sz=chr_sz

if cbar then begin

cbar_psn=[.88,0.1,0.92,0.92]; [x_min,y_min,x_max,y_max]
plt_rgn_nrm=[0.12,0.1,0.87,0.92]; [x_min,y_min,x_max,y_max]

endif else begin

plt_rgn_nrm=[0.12,0.1,0.97,0.92]; [x_min,y_min,x_max,y_max] 
if max(ord) gt 12 then plt_rgn_nrm=[0.10,.03,0.97,0.97]; [x_min,y_min,x_max,y_max]

endelse

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,0

tvlct,r_curr,g_curr,b_curr

if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
if cbar then mrg_btm=mrg_btm+4
mrg_lft=5.1 ; 10 is default
mrg_rgt=2.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.3
	mrg_btm=0.6 ; 4 is default
	if btm then mrg_btm=mrg_btm+1
	if cbar then mrg_btm=mrg_btm+1.75
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

if strpos(fl_nm,'erbe') ne -1 and time_nbr eq 60 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
ord=ord(1:52)
data=data(*,1:52)
endif

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	if (rng_y(1)+1) mod 12 eq 0 then ord_max=rng_y(1)+1 else ord_max=rng_y(1)
endelse

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	xstyle=5, $
;	ystyle=12, $
	ystyle=13, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
;	position=plt_rgn_nrm, $
	/closed, $			
	/fill, $			
	/noerase

if prn then tck_lng=0.01 else tck_lng=0.02
lon_axz_drw,lon_top=abc_max,lon_btm=abc_min,ntt=1,tck_lng=tck_lng,chr_sz=chr_sz,axz_vrt=0,tck_lcn=tck_lcn

if prn then tck_lng=0.04 else tck_lng=0.02
time_axz_drw,time_min=!y.crange(0),time_max=!y.crange(1),time_ncr=1,time_unit='mth',ntt=1,ttl=y_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=1,lbl_sty='mth_shrt',axz_vrt=1,yr_srt=yr_srt

cntr_thk(*)=1
contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

end; end gcm_xt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM Hovmoller commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM tv
; Plot seasonal cycle of regional averages
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_tv_bch, $
	x87=x87, $
	mass=mass, $ ; Mass produce
	rsmas=rsmas, $ ; RSMAS station data
	aeronet=aeronet, $ ; AERONET station data
	dly=dly, $ ; Daily rather than monthly data
	chn_idx=chn_idx, $ ; [idx] Channel
	lev_idx=lev_idx, $ ; [idx] Level
	pll=pll, $ ; [enm] Palette
	bw=bw, $ ; [flg] Use B&W for all data
	info=info, $ ; Print diagnostic information
	dbg_lvl=dbg_lvl, $ ; Debugging level
	anm=anm, $ ; Anomaly
	psym=psym, $ ; Use discrete symbols rather than connected lines
	shd=shd, $ ; [flg] Shade standard deviation region
	stp=stp, $ ; [flg] Stipple instead of shade
	stt=stt, $ ; [flg] Statistical analysis
	prn=prn
; Usage: gcm_tv_bch,dly=1
if n_elements(x87) eq 0 then x87=0 ; Exclude El Nino year
if n_elements(mass) eq 0 then mass=0 ; Mass produce
if n_elements(rsmas) eq 0 then rsmas=0 ; RSMAS station data
if n_elements(aeronet) eq 0 then aeronet=0 ; AERONET station data
if n_elements(dly) eq 0 then dly=0 ; Daily rather than monthly data
if n_elements(chn_idx) eq 0 then chn_idx=2 ; [idx] Channel
if n_elements(lev_idx) eq 0 then lev_idx=-1 ; [idx] Level
if n_elements(pll) eq 0 then pll=4 ; [enm] Palette
if n_elements(bw) eq 0 then bw=0; [flg] Use B&W for all data
if n_elements(info) eq 0 then info=1 ; Print diagnostic information
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0 ; Debugging level
if n_elements(anm) eq 0 then anm=0 ; Anomaly
if n_elements(shd) eq 0 then shd=1 ; [flg] Shade standard deviation region
if n_elements(stp) eq 0 then stp=1 ; [flg] Stipple instead of shade
if n_elements(stt) eq 0 then stt=1 ; [flg] Statistical analysis
if anm then anm_sng='anm_' else anm_sng=''
if n_elements(prn) eq 0 then prn=0

fld3_stf=[ $
;	['DSTODXC','!5AOD','!5AOD','0.63 !7l!5m']]
;	['DSTODXC','!5Dust','!5AOD','0.63 !7l!5m']]
;	['aot','!5Aerosol Optical Depth','!5Optical Depth','']]
	['DSTQ','!5Dust Concentration','!5Concentration','!7l!5g m!E-3!N']]
;	['DSTSFDPS','!5Dust Deposition','!5Deposition','!5g cm!E-2!N ka!E-1!N'] ]
;	['DSTSFDPS','!5Dust Deposition','!5Deposition','!5ng m!E-2!N s!E-1!N'] ]
;	['PRECT','!5Total Precipitation','!8P!IT!5!N','!5mm day!E-1!N']]
;	['WND_MBL','!510 m Wind Speed','!8V!5!IH!N','!5m s!E-1!N']]
;	['CLDTOT','!5Total Cloud','!8A!5!IT!N','!5%'], $
;	['CLDLOW','!5Low Cloud','!8A!5!IL!N','!5%'], $
;	['CLDMED','!5Middle Cloud','!8A!5!IM!N','!5%'], $
;	['CLDHGH','!5High Cloud','!8A!5!IH!N','!5%']]
;	['ALBEDO','!5Albedo','!8A!5','!5%'], $
;	['VMAGSFC','!5Wind Speed','V','!5m s!E-1!N']]
;	['WND_SFC','!5Surface Wind Speed','V','!5m s!E-1!N']]
;	['FSNT','!5Absorbed Shortwave','!8F!5!Eabs!N','!5W m!E-2!N'], $
;	['FSNTC','!5Absorbed Clearsky Shortwave','!8F!5!Iclr!Eabs!N','!5W m!E-2!N'] ]
;	['FSDS','!5Surface Insolation','SW','!5W m!E-2!N'] ]
;	['FSNS','!5Net Surface Shortwave','SW','!5W m!E-2!N'] ]
;	['SWCF','!5Shortwave Cloud Forcing','!5SWCF!5','!5W m!E-2!N'], $
;	['LWCF','!5Longwave Cloud Forcing','!5LWCF!5','!5W m!E-2!N'] ]
;	['TS1','!5Sea Surface Temperature','SST','!E!12_!5!NK'], $
;	['FLNS','!5Net Surface Longwave','LW','!5W m!E-2!N'], $
;	['SHFLX','!5Sensible Heat Flux','SH','!5W m!E-2!N'], $
;	['LHFLX','!5Latent Heat Flux','LH','!5W m!E-2!N'], $
;	['E_P','!5Evaporation - Precipitation','E-P','!5mm day!E-1!N'], $
;	['NET','!5Net Surface Energy','SEB','!5W m!E-2!N'] ]
fld_nbr3=n_elements(fld3_stf(0,*))
foo={fld3_sct,idx:0,nm:'',sng:'',abb:'',unit:''}
fld3=replicate({fld3_sct},fld_nbr3)
for idx=0,fld_nbr3-1 do begin
	fld3(idx).idx=idx
	fld3(idx).nm=fld3_stf(0,idx)
	fld3(idx).sng=fld3_stf(1,idx)
	fld3(idx).abb=fld3_stf(2,idx)
	fld3(idx).unit=fld3_stf(3,idx)
endfor; end loop over fld_lst

xpt_stf=[ $
;	['ecmwf','!5ECMWF'], $
;	['ncep','!5NCEP']]
;	['cmap','!5CMAP']]
;	['toms','!5TOMS'], $
;	['avhrr','!5AVHRR'], $
;	['aeronet','!5AERONET'], $
	['rsmas','!5U. Miami'], $
;	['indoex','!5INDOEX'], $
;	['dstccm15','!5CCM']]
;	['dstmch55','!7K!5!I0!N'], $
;	['dstmch56','!7K!5(D!Ip!N)'], $
;	['dstmch57','!7K!5(D!Ip!N,D!IP!N)']]
;	['dstccm51','!5dstccm51']]
;	['dstmch90','!5dstmch90']]
;	['dstmch1p','!5Uniform'], $
	['dstmchxc','!5New Geo'], $
	['dstmch1q','!5Topo'], $
	['dstmch1t','!5Geo'], $
	['dstmch1r','!5Hydro']]
;	['isccp','!5ISCCP']]
;	['erbe_b','!5ERBE'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'], $
;	['sld012d','!5CCM3']]
;	['spcp_85','!5ANV']]
;	['erbe_b','!5ERBE']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

if not mass and not rsmas then rgn_stf=[ $
;	['KGI','King George Island (62.18S, 58.30W)']]
	['Brb','Barbados (13.17N, 59.43W)']]
;	['KeC','Kern County'], $
;	['SAz','Southern Arizona']]
;	['Sahara_Arabia','Saharan and Arabian Desert'], $
;	['US_SGP_CART','Southern Great Plains CART site'], $
;	['WCA','Cape Verde'], $
;	['ArS','Arabian Sea'] ]
;	['Oma','Oman'], $
;	['Tmb','Timbuktu'], $
;	['LkC','Lake Chad'] ]
;	['Indian_Central','Central Indian Ocean'] ]
;	['Pacific_ITCZ_West','West Pacific ITCZ'] ]
;	['Arabian_Sea','Arabian Sea'] ]
;	['ocn','World Ocean'], $
;	['Vst','Vostok'], $
;	['Smm','Summit']]
if aeronet then rgn_stf=[ $
;	['AbH','Abracos Hill (10.45S,62.21W)'],$
;	['AlD','Al Dhafra (24.15N,54.32E)'],$
;	['AlF','Alta Floresta (9.55S,56.01W)'],$
;	['Any','Anmyon (36.31N,126.19E)'],$
;	['ArA','Aire Adour (43.42N,0.15E)'],$
;	['Arc','Arica (18.28S,70.18W)'],$
;	['AsI','Ascebsuib Island (7.58S,14.24W)'],$
;	['Avg','Avignon (43.55N,4.52E)'],$
;	['Bbn','Balbina (1.55S,59.29W)'],$
	['BdB','Bidi Bahn (14.06N,2.45W)'],$
;	['Bdk','Bondoukoui (11.50N,3.45W)'],$
;	['BDV','BONDVILLE (40.03N,88.22W)'],$
	['Bhr','Bahrain (26.33N,50.50E)'],$
;	['Bkh','Brookhaven (40.52N,72.53W)'],$
;	['BnC','Bonanza Creek (64.44N,148.18W)'],$
	['Bnz','Banizoumbou (13.54N,2.67E)'],$
	['Brb','Barbados (13.17N,59.43W)'],$
	['Brm','Bermuda (32.37N,64.70W)'],$
;	['Bsl','Brasilia (15.55S,47.53W)'],$
;	['BtL','Bratts Lake (50.16N,104.42W)'],$
;	['Btr','Belterra (2.38S,54.57W)'],$
;	['Bts','Burtonsville (39.06N,76.56W)'],$
;	['CAT','CARTEL (45.22N,71.55W)'],$
;	['Cba','Cuiaba (15.30S,56.00W)'],$
;	['Chh','Chinhae (35.09N,128.39E)'],$
;	['CmF','Clermont Ferrand (45.45N,2.57E)'],$
;	['Cnp','Concepcion (16.08S,62.01W)'],$
	['CpV','Capo Verdo (16.73N,22.94W)'],$
;	['Cqm','Chequamegon (45.55N,90.15W)'],$
;	['Crt','Creteil (48.47N,2.26E)'],$
;	['CtS','Cart Site (36.36N,97.24W)'],$
	['Dak','Daker (14.23N,16.57W)'],$
;	['DdS','Dead Sea (31.06N,35.27E)'],$
;	['Dhw','Dharwar (15.25N,74.59E)'],$
	['DlZ','Dalanzadgad (43.58N,104.42E)'],$
;	['DrT','Dry Tortugas (24.60N,82.80W)'],$
	['DsI','Dongsha Island (20.41N,116.04E)'],$
	['Duh','Dunhang (40.04N,94.79E)'],$
	['Ebt','Egbert (43.13N,79.45W)'],$
	['ElA','El Arenosillo (37.10N,6.70W)'],$
;	['FLF','FLIN FLON (54.08N,101.41W)'],$
;	['GID','GOA INDIA (15.27N,73.48E)'],$
;	['GSF','GSFC (39.01N,76.52W)'],$
;	['Gtl','Gotland (57.55N,18.56E)'],$
;	['Gts','Gaithersberg (39.07N,77.12W)'],$
;	['HgI','Hog Island (37.25N,75.42W)'],$
;	['HJA','HJAndrews (44.14N,122.13W)'],$
;	['Hol','Howland (45.12N,60.43W)'],$
;	['HpR','Hampton Roads (36.46N,76.27W)'],$
;	['Ilr','Ilorin (8.19N,4.20E)'],$
;	['IME','IMS METU ERDEMLI (36.33N,.15E)'],$
	['InM','Inner Mongolia (42.68N,115.95E)'],$
;	['Isp','Ispra (45.48N,8.37E)'],$
	['Izn','Izana (28.18N,16.30W)'],$
;	['Jbg','Joberg (26.11S,28.01E)'],$
;	['JgB','Jug Bay (38.46N,76.46W)'],$
;	['JiP','Ji Parana (10.51S,61.47W)'],$
	['Kaa','Kaashidhoo (4.97N,73.47E)'],$
;	['KBc','Key Biscayne (25.43N,80.09W)'],$
;	['Kjk','Kejimkujik (44.22N,65.16W)'],$
;	['Klm','Kaloma (14.51S,24.49E)'],$
;	['Ksm','Kasama (10.10S,31.10E)'],$
;	['LF9','LOS FIEROS 98 (14.33S,60.55W)'],$
;	['Lil','Lille (50.36N,3.08E)'],$
;	['Lni','Lanai (20.44N,156.55W)'],$
;	['MAL','MALE (4.11N,73.31E)'],$
;	['Mds','Madison (43.04N,89.24W)'],$
;	['Mdv','Moldova (47.01N,28.45E)'],$
;	['Mfw','Mfuwe (13.15S,31.55E)'],$
;	['Mng','Mongu (15.15S,23.09E)'],$
	['MnL','Mauna Loa (19.54N,155.58W)'],$
;	['Mwn','Mwinilunga (11.44S,24.25E)'],$
	['Nar','Naura (0.31S,166.54E)'],$
;	['NCU','NCU Taiwan (24.53N,121.05E)'],$
;	['Not','Noto (37.20N,137.08E)'],$
;	['NYB','NSA YJP BOREAS (55.54N,98.17W)'],$
	['Ogd','Ouagadougou (12.11N,1.23W)'],$
	['Okn','Okinawa (26.36N,127.77E)'],$
;	['Ost','Oyster (37.17N,75.55W)'],$
;	['Pdw','Paddockwood (53.30N,105.30W)'],$
;	['Pls','Palaiseau (48.42N,2.12E)'],$
;	['RmH','Rame Head (50.21N,4.08W)'],$
;	['RnI','Rottnest Island (32.00S,115.17E)'],$
;	['SeB','Sede Boker (30.52N,34.47E)'],$
;	['SER','SERC (38.52N,76.30W)'],$
;	['Shm','Shirahama (33.41N,135.21E)'],$
;	['Skk','Skukuza (24.59S,31.35E)'],$
;	['SlS','Seoul SNU (37.27N,126.57E)'],$
	['SlV','Solar Village (24.54N,46.24E)'],$
;	['SMH','SMHI (58.34N,16.08E)'],$
	['SNc','San Nicolas (33.15N,119.29W)'],$
;	['Sng','Senange (16.06S,23.17E)'],$
;	['Spt','Sopot (54.27N,18.33E)'],$
;	['Srn','Surinam (5.47N,55.12W)'],$
;	['Ssk','Sesheke (17.28S,24.18E)'],$
;	['Str','Santarem (2.25S,54.45W)'],$
;	['Svl','Sevillela (34.21N,106.53W)'],$
;	['SYB','SSA YJP BOREAS (53.40N,104.39W)'],$
	['TbM','Table Mountain (40.07N,105.14W)'],$
	['Tcs','Tucson (32.13N,110.57W)'],$
	['Tht','Tahiti (17.34S,149.36W)'],$
;	['Tkr','Tukurui (3.43S,49.40W)'],$
;	['Tmp','Thompson (55.47N,97.50W)'],$
;	['Tms','Tombstone (31.44N,110.02W)'],$
	['Tnr','Tenerife (28.03N,16.63W)'],$
;	['Tul','Toulouse (43.34N,1.22E)'],$
	['Ulg','Ulaangom (49.97N,92.08E)']]
;	['USD','USDA (39.01N,76.52W)'],$
;	['Vns','Venise (45.18N,12.30E)'],$
;	['Wks','Waskesiu (53.55N,106.04W)'],$
;	['Wlp','Wallops (37.56N,75.28W)'],$
;	['Zbz','Zambezi (13.31S,23.06E)']]
if rsmas then rgn_stf=[ $
	['AmS','American Samoa (14.25S, 170.58W)'], $ 
	['Brb','Barbados (13.17N, 59.43W)'], $
	['Brm','Bermuda (32.27N, 64.87W)'], $
	['CGH','Cape Point (34.35S, 18.48E)'], $ 
	['CGr','Cape Grim (40.68S, 144.68E)'], $ 
	['Enw','Enewetak (11.33N, 162.30E)'], $ 
	['Fnf','Funafuti (8.50S, 179.20W)'], $ ; unreliable according to Savoie 20020715
	['Izn','Izania (28.30N, 16.50W)'], $
	['Jej','Jeju (33.52N, 126.48E)'], $ 
	['KGI','King George Island (62.18S, 58.30W)'], $ 
	['Kaa','Kaashidhoo (4.95N, 73.45E)'], $
	['McH','Mace Head (53.32N, 9.85W)'], $ 
	['Mdw','Midway (28.22N, 177.35W)'], $ 
	['Mmi','Miami (25.75N, 80.25W)'], $ 
	['Nau','Nauru (0.53S, 166.95E)'], $ 
	['NeC','New Caledonia (21.15S, 167.00E)'], $ 
	['Nrf','Norfolk Is. (29.08S, 167.98E)'], $
	['Oah','Oahu (21.33N, 157.70W)'], $
	['Okn','Okinawa (26.92N, 128.25E)'], $
;;	['Rly','Roleystone (32.12S, 116.07E)']] ; No dust data
	['Rrt','Rarotonga (21.25S, 159.75W)'], $
	['CGr','Cape Grim (40.68S, 144.68E)'], $ 
	['SlI','Sal Island (16.75N, 22.92W)']]
if mass then rgn_stf=[ $
	['Africa_South','South Africa'], $
	['Alaska_NW_Canada','Alaska NW Canada'], $
	['Amazon_Basin','Amazon Basin'], $
	['America_Central','Central America'], $
	['America_South_South','Southern South America'], $
	['Antarctica','Antarctica'], $
	['Atlantic_North','North Atlantic'], $
	['Australia','Australia'], $
	['Congo','Congo'], $
	['Europe_Northern','Northern Europe'], $
	['Greenland','Greenland'], $
	['India','India'], $
	['Indian_Tropical','Tropical Indian Ocean'], $
	['Indochina','Indochina'], $
	['Indonesia_Land','Indonesia (Land)'], $
	['Indonesia_Ocean','Indonesia (Maritime)'], $
	['Pacific_Equatorial','Equatorial Pacific'], $
	['Pacific_Equatorial_Central','Central Equatorial Pacific'], $
	['Pacific_Equatorial_Eastern','Eastern Equatorial Pacific'], $
	['Pacific_Equatorial_Western','Western Equatorial Pacific'], $
	['Pacific_North','North Pacific'], $
	['Pacific_South','South Pacific'], $
	['Pacific_Tropical','Tropical Pacific'], $
	['Pacific_Tropical_Central','Tropical Central Pacific'], $
	['Pacific_Tropical_Eastern','Tropical Eastern Pacific'], $
	['Pacific_Tropical_Western','Tropical Western Pacific'], $
	['Pacific_Western_Warm_Pool','Western Pacific Warm Pool'], $
	['Siberia_Eastern','Eastern Siberia'], $
	['Siberia_Western','Western Siberia'], $
	['Tibetan_Plateau','Tibetan Plateau'], $
	['US_Central','Central US'], $
	['US_Eastern','Eastern US'], $
	['US_Western','Western US'] ] 
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct_gcm,idx:0,nm:'',sng:''}
rgn=replicate({rgn_sct_gcm},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	rgn(idx).idx=idx
	rgn(idx).nm=rgn_stf(0,idx)
	rgn(idx).sng=rgn_stf(1,idx)
endfor; end loop over rgn_lst

for fld_idx3=0,fld_nbr3-1 do begin
for rgn_idx=0,rgn_nbr-1 do begin
	y_ttl=fld3(fld_idx3).abb
	unit_sng=fld3(fld_idx3).unit
	if xpt(0).nm eq 'aeronet' then unit_sng=aeronet_unit_sng_get(chn_idx)
	if unit_sng ne '' then y_ttl=y_ttl+' ('+unit_sng+')'
	if dly then x_ttl='' else x_ttl='Month of Year'
	if dly then wrp=0 else wrp=1
	top=1
	btm=1
	fl_nm_out_xpt_sng=''
	scl=1.0+0.0*fltarr(xpt_nbr)
	fld_nm=strarr(xpt_nbr)
	fl_nm_in=strarr(xpt_nbr)
	lbl_sng=strarr(xpt_nbr)
	if n_elements(psym) ne xpt_nbr then psym=0*intarr(xpt_nbr)
	obs_nbr=xpt_nbr
	mdl_nbr=0
	for xpt_idx=0,xpt_nbr-1 do begin
		fld_nm(xpt_idx)=fld3(fld_idx3).nm
		xpt_nm=xpt(xpt_idx).nm
		rgn_nm=rgn(rgn_idx).nm
		yr_sng='clm'
		mth_sng='0112'
		yr_srt=90
		date_sng=yr_sng+'_'+mth_sng
		if fld_nm(xpt_idx) eq 'DSTQ' and xpt_nm eq 'rsmas' then fld_nm(xpt_idx)='cnc_mss_dst'
		if fld_nm(xpt_idx) eq 'DSTODXC' and xpt_nm eq 'aeronet' then fld_nm(xpt_idx)='aot'
		if strpos(fld_nm(xpt_idx),'DSTSF') ne -1 then scl(xpt_idx)=1.0e12
		if strpos(fld_nm(xpt_idx),'DSTQ') ne -1 then scl(xpt_idx)=1.0e9
		if strpos(fld_nm(xpt_idx),'cnc_mss_') ne -1 then scl(xpt_idx)=1.0e9
; Convert mixing ratio to concentration: 1013 mb = 1.15 kg m-3, 750 mb = 0.9 kg m-3
		if strpos(xpt_nm,'dst') ne -1 and strpos(fld_nm(xpt_idx),'DSTQ') ne -1 and rsmas then scl(xpt_idx)=scl(xpt_idx)*1.15
; Izana is special, unscale previous scaling then multiply by Izana-specific scaling
; Savoie stores Izania data as volumetric mass concentration at sea level. 
; To convert RSMAS data back to ambient concentration at 750 mb multiply by 0.9/1.15
		if strpos(xpt_nm,'dst') ne -1 and strpos(fld_nm(xpt_idx),'DSTQ') ne -1 and rgn_nm eq 'Izn' then scl(xpt_idx)=scl(xpt_idx)*0.9/1.15
;		Lubbock runs
		if strpos(xpt_nm,'dstmch1j') ne -1 or strpos(xpt_nm,'dstmch1l') ne -1 or strpos(xpt_nm,'dstmch1m') ne -1 or strpos(xpt_nm,'dstmch1n') ne -1 and rsmas then scl(xpt_idx)=scl(xpt_idx)*1.5
; Define everything positive into surface for SEB studies
		if fld_nm(xpt_idx) eq 'E_P' then scl(xpt_idx)=8.64e7
		if fld_nm(xpt_idx) eq 'FLNS' then scl(xpt_idx)=-1.0
		if fld_nm(xpt_idx) eq 'SHFLX' then scl(xpt_idx)=-1.0
		if fld_nm(xpt_idx) eq 'LHFLX' then scl(xpt_idx)=-1.0
		if fld_nm(xpt_idx) eq 'CLDLOW' then scl(xpt_idx)=100.0
		if fld_nm(xpt_idx) eq 'CLDMED' then scl(xpt_idx)=100.0
		if fld_nm(xpt_idx) eq 'CLDHGH' then scl(xpt_idx)=100.0
		if fld_nm(xpt_idx) eq 'CLDTOT' then scl(xpt_idx)=100.0
		if fld_nm(xpt_idx) eq 'ALBEDO' then scl(xpt_idx)=100.0
	 	if strpos(xpt_nm,'dst') ne -1 then begin
			obs_nbr=obs_nbr-1
			mdl_nbr=mdl_nbr+1
		endif
	 	if xpt_nm eq 'isccp' then yr_sng='8388'
	 	if xpt_nm eq 'rsmas' then begin
			yr_sng='clm'
			mth_sng='0112'
			yr_srt=90
			date_sng=yr_sng+'_'+mth_sng
		endif; endif rsmas
		if dly then begin
			yyyymmdd_srt=19940901 ; LITE
			yyyymmdd_end=19940930 ; LITE
			yyyymmdd_srt=19990208 ; INDOEX
			yyyymmdd_end=19990322 ; INDOEX
;			yyyymmdd_srt=19980215 ; Kaashidhoo
			yyyymmdd_srt=19950101 ; 
			yyyymmdd_end=19970729 ; 
			yyyymmdd_srt=19980101 ; 1998
			yyyymmdd_end=19981231 ; 1998
			yyyymmdd_srt_sng=string(format='(I8.8)',yyyymmdd_srt)
			yyyymmdd_end_sng=string(format='(I8.8)',yyyymmdd_end)
			date_sng=yyyymmdd_srt_sng+'_'+yyyymmdd_end_sng
		endif; endif dly
	 	if xpt_nm eq 'cmap' then begin
;			yr_sng='7901_9801'
;			yr_sng='7998_0112'
			yr_srt=79
;			mth_sng=''
		endif; endif cmap
		fl_nm_in(xpt_idx)=getenv('DATA')+'/'+xpt_nm+'/'+xpt_nm+'_'+date_sng+'_'+rgn_nm+'.nc'
		lbl_sng(xpt_idx)=xpt(xpt_idx).sng
		fl_nm_out_xpt_sng=fl_nm_out_xpt_sng+xpt_nm+'_'
	endfor; end loop over experiment
	fl_nm_out=getenv('DATA')+'/ps/'+fl_nm_out_xpt_sng+date_sng+'_'+rgn_nm+'_'+fld_nm(0)+'.eps'
	if prn then begin
		x_sz=6.5
		y_sz=3.0
		chr_sz=1.5
	endif; endif prn
	if prn and not mass then begin
		if fld_idx3 eq 0 then top=1 else top=0
		if fld_idx3 eq fld_nbr3-1 then btm=1 else btm=0
	endif; endif prn
	if prn and rsmas then begin
		if rgn_nm eq 'NeC' or rgn_nm eq 'CGr' or rgn_nm eq 'Nrf' then btm=1 else btm=0
	endif; endif prn
	if prn then begin
		if top then y_sz=y_sz*1.1
		if btm then if dly then y_sz=y_sz*1.15 else y_sz=y_sz*1.1
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_tv,x_ttl=x_ttl,y_ttl=y_ttl,ttl=rgn(rgn_idx).sng,psym=psym,rgn_nm=rgn_nm,wrp=wrp,fld_nm=fld_nm,prn=prn,top=top,btm=btm,fl_nm=fl_nm_in,info=info,obs_nbr=obs_nbr,mdl_nbr=mdl_nbr,lbl_sng=lbl_sng,dbg_lvl=dbg_lvl,chn_idx=chn_idx,lev_idx=lev_idx,bw=bw,pll=pll,clr_ord=clr_ord,clr_tbl_nbr=clr_tbl_nbr,scl=scl,shd=shd,stp=stp,stt=stt,dly=dly,yyyymmdd_srt=yyyymmdd_srt,yr_srt=yr_srt,unit_sng=unit_sng

	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over rgn_lst
endfor; end loop over fld_lst

end; end gcm_tv_bch()
pro gcm_tv, $
	rgn_nm=rgn_nm, $
	fl_nm=fl_nm, $
	ttl=ttl, $
	psym=psym, $ ; Use discrete symbols rather than connected lines
	dbg_lvl=dbg_lvl, $ ; Debugging level
	info=info, $
	obs_nbr=obs_nbr, $
	mdl_nbr=mdl_nbr, $
	chn_idx=chn_idx, $ [idx] Channel
	lev_idx=lev_idx, $ [idx] Level
	pll=pll, $ [enm] Palette
	bw=bw, $ ; [flg] Use B&W for all data
	clr_ord=clr_ord, $
	clr_tbl_nbr=clr_tbl_nbr, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	wrp=wrp, $
	shd=shd, $ ; [flg] Shade standard deviation region
	stp=stp, $ ; [flg] Stipple instead of shade
	stt=stt, $ ; [flg] Statistical analysis
	lbl_sng=lbl_sng, $
	top=top, $
	btm=btm, $
	chr_sz=chr_sz, $
	prn=prn, $
	scl=scl, $
	dly=dly, $ ; Daily rather than monthly data
	yyyymmdd_srt=yyyymmdd_srt, $ ; Start date in YYYYMMDD format
	yr_srt=yr_srt, $
	unit_sng=unit_sng, $
	fld_nm=fld_nm

@ibp_clr.com

; gcm_tv,y_ttl='LWCF (!5W m!e-2!N)',ttl='India',rgn_nm='India',fl_nbr=3,wrp=1

if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(psym) eq 0 then psym=[0] ; 0: line, 1: plus, 2:asterisk, 3:dot, 4:diamond, 5:triangle, 6:square, 7:x, 8:user defined, 
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(fld_nm) eq 0 then fld_nm=['LWCF']
if n_elements(rgn_nm) eq 0 then rgn_nm='India'
if n_elements(y_ttl) eq 0 then y_ttl='LWCF (!5W m!e-2!N)'
if n_elements(x_ttl) eq 0 then x_ttl='Month of Year'
if n_elements(info) eq 0 then info=1
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0 ; Debugging level
if n_elements(chn_idx) eq 0 then chn_idx=2 ; [idx] Channel
if n_elements(lev_idx) eq 0 then lev_idx=-1 ; [idx] Level
if n_elements(pll) eq 0 then pll=4 ; [enm] Palette
if n_elements(bw) eq 0 then bw=0; [flg] Use B&W for all data
if n_elements(clr_ord) eq 0 then clr_ord=1
if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=22
if n_elements(wrp) eq 0 then wrp=1
if n_elements(shd) eq 0 then shd=1 ; [flg] Shade standard deviation region
if n_elements(stp) eq 0 then stp=1 ; [flg] Stipple instead of shade
if n_elements(stt) eq 0 then stt=1 ; [flg] Statistical analysis
if n_elements(lbl_sng) eq 0 then lbl_sng=''
if n_elements(yr_srt) eq 0 then yr_srt=85
if n_elements(prn) eq 0 then prn=0
if n_elements(dly) eq 0 then dly=0 ; Daily rather than monthly data
if n_elements(yyyymmdd_srt) eq 0 then yyyymmdd_srt=19990312 ; Start date in YYYYMMDD format
if n_elements(scl) eq 0 then scl=1.0+0.0*fltarr(n_elements(fl_nm))
if n_elements(ttl) eq 0 then ttl='India'
if n_elements(unit_sng) eq 0 then unit_sng=''
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/spcp_85/spcp_85_anm_India_8589x87_0112.nc'

color_order=clr_ord
clr_rfr,clr_tbl_nbr,"",0
erase
!p.multi=0
sym_usr_foo=findgen(17)*(!pi*2/16.)
usersym,cos(sym_usr_foo),sin(sym_usr_foo),/fill
sym_sz=2.0

fl_nbr=n_elements(fl_nm)
if fl_nbr eq 1 then fl_nm=[fl_nm]
if n_elements(obs_nbr) eq 0 then obs_nbr=0
if n_elements(mdl_nbr) eq 0 then mdl_nbr=fl_nbr

; Determine size of longest file
time_nbr_max=0
for fl_idx=0,fl_nbr-1 do begin
	nc_id=ncdf_open(fl_nm(fl_idx))
	time_id=ncdf_dimid(nc_id,'time')
	ncdf_diminq,nc_id,time_id,dim_foo,time_nbr
	ncdf_varget,nc_id,'time',time_tnt
	ncdf_close,nc_id
	if time_nbr gt time_nbr_max then begin
		time_nbr_max=time_nbr
		time=time_tnt
	endif
endfor; end loop over files
data=1.0e36+fltarr(time_nbr_max,fl_nbr)
data_sdn=1.0e36+fltarr(time_nbr_max,fl_nbr)
sdn_flg=0*intarr(fl_nbr)

; Get data from netCDF files
for fl_idx=0,fl_nbr-1 do begin
	nc_id=ncdf_open(fl_nm(fl_idx))
	print,'Processing '+fl_nm(fl_idx)

	; Get number of dimensions in input file
	; var_sct={name:'',id:0L,data_type:'',dmn_nbr:0L,nbr_atts:0L,dmn:lonarr(4)}
	var_id=ncdf_varid(nc_id,fld_nm(fl_idx))
	var_crr=ncdf_varinq(nc_id,var_id)
	dmn_nbr=var_crr.ndims
	dmn_id=var_crr.dim
	srt=0*lonarr(dmn_nbr) ; [idx] Starting offset
	cnt=1+lonarr(dmn_nbr) ; [nbr] Count
	; Read required dimensions
	time_id=ncdf_dimid(nc_id,'time')
	ncdf_diminq,nc_id,time_id,dim_foo,time_nbr
	ncdf_varget,nc_id,'time',time_crr
	if time_nbr ne time_nbr_max then print,'WARNING: timeseries are not of equal length, results only valid if missing data are at end of timeseries'
	; AERONET data has channel dimension
	chn_id=ncdf_dimid(nc_id,'ch')
	if chn_id gt -1 then begin
	; Channel dimension is in file, retrieve channel size and coordinate
		ncdf_diminq,nc_id,chn_id,foo,chn_nbr
;		ncdf_varget,nc_id,'chn',chn
	endif else chn_nbr=0; endif chn
	; Recall chn_idx is passed as 0-based index
	if chn_idx gt -1 and chn_nbr ge 1 then begin
		if dmn_nbr eq 2 then srt=[chn_idx,0]
		if dmn_nbr eq 2 then cnt=[1,time_nbr]
	endif; endif chn
	; Model data may have level dimension
	lev_id=ncdf_dimid(nc_id,'lev')
	if lev_id gt -1 then begin
	; Level dimension is in file, retrieve level size and coordinate
		ncdf_diminq,nc_id,lev_id,foo,lev_nbr
;		ncdf_varget,nc_id,'lev',lev
	; Default to near-surface value
		if lev_idx le -1 then lev_idx=lev_nbr-1
	endif else lev_nbr=0; endif lev
	; Recall lev_idx is passed as 0-based index
	if lev_idx gt -1 and lev_nbr ge 1 then begin
		if dmn_nbr eq 2 then srt=[lev_idx,0]
		if dmn_nbr eq 2 then cnt=[1,time_nbr]
	endif; endif lev
	if dmn_nbr eq 1 then cnt=[time_nbr]
	; Get variable data
	ncdf_varget,nc_id,fld_nm(fl_idx),data_foo,offset=srt,count=cnt
	; Do not scale missing data
	good_idx=where(data_foo lt 1.0e20,good_data_nbr)
	if good_data_nbr gt 0 then data_foo(good_idx)=scl(fl_idx)*data_foo(good_idx)
	data(0:time_nbr-1,fl_idx)=reform(data_foo)

	; Get standard deviation if present
	fld_nm_sdn=fld_nm(fl_idx)+'_sdn'
	var_sdn_id=ncdf_varid(nc_id,fld_nm_sdn)
	if var_sdn_id gt -1 then begin
	; Standard deviation is in file
		sdn_flg(fl_idx)=1
		ncdf_varget,nc_id,fld_nm_sdn,data_sdn_foo,offset=srt,count=cnt
		good_idx=where(data_sdn_foo lt 1.0e20,good_data_nbr)
		if good_data_nbr gt 0 then data_sdn_foo(good_idx)=scl(fl_idx)*data_sdn_foo(good_idx)
		data_sdn(0:time_nbr-1,fl_idx)=reform(data_sdn_foo)
	endif else print,'Standard deviation not present for '+fld_nm_sdn
	ncdf_close,nc_id
endfor; end loop over files
time_nbr=time_nbr_max

ord=reform(data)
ord_sdn=reform(data_sdn)
if dly then begin
	abc=indgen(time_nbr)
	time_unit='day'
	lbl_sty='day' ; [day, mth_shrt, mth_mdm, mth_lng, none]
	yyyymmdd_prs,yyyymmdd_srt,yr=yr_srt,mth=mth_srt,day=day_srt
;	foo=label_date(date_format='%D',offset=julday(mth_srt,day_srt,yr_srt))
	foo=label_date(date_format='%D!C%M',offset=julday(mth_srt,day_srt,yr_srt))
endif else begin
	abc=indgen(time_nbr)
	time_unit='mth'
	lbl_sty='mth_shrt' ; [day, mth_shrt, mth_mdm, mth_lng, none]
endelse; endif dly
if wrp then begin
	ord=[ord,ord(0,*)]
	ord_sdn=[ord_sdn,ord_sdn(0,*)]
	abc=[abc,time_nbr]
endif

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse ; endif
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord(where(ord lt 1.0e20)))
	ord_max=max(ord(where(ord lt 1.0e20)))
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse ; endif

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.5 ; 2 is default
	if top then mrg_top=mrg_top+1.3
	mrg_btm=0.5 ; 4 is default
	if btm then if dly then mrg_btm=mrg_btm+2.3 else mrg_btm=mrg_btm+1.3
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

clr_mk,fl_nbr,clr_ord,pll,0
tvlct,r_curr,g_curr,b_curr
clr=indgen(fl_nbr)+2
ln_sty=0*indgen(fl_nbr)
ln_thk=3.0 ; [frc] Line thickness
ln_thk_sdn=2.0 ; [frc] Line thickness for standard deviation whiskers
for obs_idx=0,obs_nbr-1 do begin
	if obs_idx ne 0 then ln_sty(obs_idx)=2
	clr(obs_idx)=clr_blk_idx
endfor; end loop over obs
if obs_nbr eq 0 then base_clr_idx=2 else base_clr_idx=obs_nbr+1
for mdl_idx=0,mdl_nbr-1 do begin
	if bw then clr(obs_nbr+mdl_idx)=clr_blk_idx else clr(obs_nbr+mdl_idx)=base_clr_idx+mdl_idx
endfor; end loop over fields
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
orn_idx=base_clr_idx+3
yll_idx=base_clr_idx+4
gry_idx=100
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(gry_idx)=200
g_curr(gry_idx)=200
b_curr(gry_idx)=200
tvlct,r_curr,g_curr,b_curr

for fl_idx=0,fl_nbr-1 do begin

if fl_idx eq 0 then begin
plot, $
	abc, $
	ord(*,fl_idx), $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=13, $
	ystyle=0, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	thick=ln_thk, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	color=clr_blk_idx, $
	psym=psym(fl_idx), $
	linestyle=ln_sty(fl_idx), $
	/nodata

; Subsequent shading obscures any previous drawing so redraw later if necessary
oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=ln_thk,color=clr(fl_idx),symsize=sym_sz,psym=psym(fl_idx),linestyle=ln_sty(fl_idx)

endif else begin; fl_idx != 0

oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=ln_thk,color=clr(fl_idx),symsize=sym_sz,psym=psym(fl_idx),linestyle=ln_sty(fl_idx)

endelse; fl_idx != 0

if dbg_lvl then begin
	print,'gcm_tv abc = ',abc
	print,'gcm_tv ord = ',ord(*,fl_idx)
endif; endif dbg

; Turn off standard deviations when there are too many datasets
if total(sdn_flg) gt 2 then sdn_flg(1:fl_nbr-1)=0

; If standard deviation data exist and will not look too cluttered
if sdn_flg(fl_idx) then begin

wsk_lng=2.0 ; [nbr] Number of standard deviations in full whisker length
ord_sdn_top=ord+0.5*wsk_lng*ord_sdn
ord_sdn_btm=ord-0.5*wsk_lng*ord_sdn
bad_idx=where(ord_sdn_btm lt !y.crange(0),bad_data_nbr)
;if bad_data_nbr gt 0 then ord_sdn_btm(bad_idx)=!y.crange(0)
bad_idx=where(ord_sdn_top gt !y.crange(1),bad_data_nbr)
;if bad_data_nbr gt 0 then ord_sdn_top(bad_idx)=!y.crange(1)
clp_rct=[!x.crange(0)+0.1,!y.crange(0),!x.crange(1)-0.1,!y.crange(1)] ; [] Clipping rectangle

; Pre-save current fl_idx
fl_idx_tmp=fl_idx

; Shade region if requested
if shd and fl_idx ne 0 then begin

; Stippling
if stp then begin
	pat=bytarr(6,6)
	pat(3:3)=255
endif; endif stp

; Hardward clipping is not supported for solid or pattern fills
; Clip in software using mgh_polyclip as suggested by Mark Hadfield on comp.lang.idl-pvwave on 20020326
; http://katipo.niwa.cri.nz/~hadfield/gust/software/idl/
clp_plg=fltarr(2,2*n_elements(abc)) ; [] Clipping polygon
clp_plg(0,*)=[abc,reverse(abc)] ; [] Clipping polygon
clp_plg(1,*)=[ord_sdn_btm(*,fl_idx),reverse(ord_sdn_top(*,fl_idx))]
clp_plg=mgh_polyclip(clp_rct(0),0,0,clp_plg,COUNT=clp_cnt) ; [] Clip left of left y-axis
clp_plg=mgh_polyclip(clp_rct(2),0,1,clp_plg,COUNT=clp_cnt) ; [] Clip right of right y-axis
clp_plg=mgh_polyclip(clp_rct(1),1,0,clp_plg,COUNT=clp_cnt) ; [] Clip beneath bottom x-axis
clp_plg=mgh_polyclip(clp_rct(3),1,1,clp_plg,COUNT=clp_cnt) ; [] Clip above top x-axis
if n_elements(pat) ne 0 then polyfill,clp_plg,color=gry_idx,pat=pat,clip=clp_rct else polyfill,clp_plg,color=gry_idx,clip=clp_rct

;if n_elements(pat) ne 0 then polyfill,[abc,reverse(abc)],[ord_sdn_btm(*,fl_idx),reverse(ord_sdn_top(*,fl_idx))],color=gry_idx,pat=pat,clip=clp_rct else polyfill,[abc,reverse(abc)],[ord_sdn_btm(*,fl_idx),reverse(ord_sdn_top(*,fl_idx))],color=gry_idx,clip=clp_rct

; Draw line for current experiment
oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=ln_thk,color=clr(fl_idx),symsize=sym_sz,psym=psym(fl_idx),linestyle=ln_sty(fl_idx)

; Shading may obscure previous drawing so redraw now if necessary
; Redraw observations and whiskers
fl_idx=0
oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=ln_thk,color=clr(fl_idx),symsize=sym_sz,psym=psym(fl_idx),linestyle=ln_sty(fl_idx)
; Jump to draw whiskers
goto,drw_wsk
; Return after jump
drw_wsk_rtn: foo=1
; Reset fl_idx in case jump occurred
fl_idx=fl_idx_tmp

endif else begin ; endif shd

; Jump to here if redrawing for shade
drw_wsk: foo=1
for idx=0,n_elements(abc)-1 do begin
	; Draw whiskers iff data are valid and iff deviations do not fill up plot
	;if ord(idx,fl_idx) lt 1.0e20 and ord_sdn_top(idx,fl_idx)-ord_sdn_btm(idx,fl_idx) lt !y.crange(1)-!y.crange(0) then $
	if ord(idx,fl_idx) lt 1.0e20 and ord_sdn_top(idx,fl_idx) lt 1.0e20 then $
	plots,[abc(idx),abc(idx)],!y.crange(0) > [ord_sdn_btm(idx,fl_idx),ord_sdn_top(idx,fl_idx)] < !y.crange(1),thick=ln_thk_sdn,symsize=sym_sz,psym=0,color=clr(fl_idx),clip=clp_rct
endfor; end loop over abc
; Return whence jumped
if shd and fl_idx eq 0 then goto,drw_wsk_rtn

endelse; endif not shd

endif; endif sdn

endfor; end loop over files

if btm then ntt=1 else ntt=0
time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit=time_unit,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=0.02,dbg_lvl=dbg_lvl,lbl_sty=lbl_sty,axz_vrt=0,yr_srt=yr_srt,yyyymmdd_srt=yyyymmdd_srt 

if info then begin

ln_lgn_x1=0.10
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

ord_avg=0.0*fltarr(fl_nbr)
avg_sng=strarr(fl_nbr)
; NB: Sample for average only where data exists in first (observational) dataset
if wrp then good_idx=where(ord(0:time_nbr-1,0) lt 1.0e20,good_ord_nbr) else good_idx=where(ord(*,0) lt 1.0e20,good_ord_nbr)
for fl_idx=0,fl_nbr-1 do begin
	if good_ord_nbr gt 0 then ord_avg(fl_idx)=total(ord(good_idx,fl_idx))/good_ord_nbr
;	avg_sng(fl_idx)=auto_sng(ord_avg(fl_idx),2)+' '+unit_sng
	avg_sng(fl_idx)=auto_sng(ord_avg(fl_idx),2)
	print,'avg_sng(fl_idx) = ',avg_sng(fl_idx)
endfor; end loop over files

for fl_idx=0,fl_nbr-1 do begin

; Plot symbol marker for each experiment
if psym(fl_idx) eq 0 then plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,symsize=sym_sz,psym=psym(fl_idx),/NORMAL else plots,0.5*(ln_lgn_x1+ln_lgn_x2),lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,symsize=sym_sz,psym=psym(fl_idx),/NORMAL

; Perform statistical analysis of correlation between models and observations
if stt then begin
x_data=ord(0:time_nbr-1,0) ; Do not include wrp point, if any
y_data=ord(0:time_nbr-1,fl_idx) ; Do not include wrp point, if any
print,'gcm_tv_bch(): Started with ',n_elements(x_data),' points'
good_idx=where(x_data lt 1.0e20,good_data_nbr)
x_data=x_data(good_idx)
y_data=y_data(good_idx)
good_idx=where(y_data lt 1.0e20,good_data_nbr)
x_data=x_data(good_idx)
y_data=y_data(good_idx)
print,'gcm_tv_bch(): Ended with ',n_elements(x_data),' points'
fit_cff=poly_fit(x_data,y_data,1,ord_fit)

data=fltarr(2,n_elements(x_data))
data(0,*)=x_data
data(1,*)=y_data
pearson_mtx=correl_matrix(data)

m_sng=auto_sng(fit_cff(1),2)
pearson_sng=string(format='(F5.2)',pearson_mtx(0,1))
pearson_sng='!8r!5='+pearson_sng
endif else begin
pearson_sng=''
endelse ; endif stt

xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx)+' '+avg_sng(fl_idx)+' '+pearson_sng,size=chr_sz,/NORMAL

endfor; end loop over fl

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_tv()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM tv
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM Q1Q2
; Plot the time-domain average Q1 and Q2 from the GCM
; ncwa -w gw -a lon,lat -d lat,-10.0,10.0 -d lon,140.0,170.0 ${DATA}/amip5/amip5_8589_07_q1q2.nc ${DATA}/amip5/amip5_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc
; ncwa -w gw -a lon,lat -d lat,-10.0,10.0 -d lon,140.0,170.0 ${DATA}/spcp_85/spcp_85_8589_07_q1q2.nc ${DATA}/spcp_85/spcp_85_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc
; ncwa -w gw -a lon,lat -d lat,-20.0,20.0 -d lon,130.0,180.0 ${DATA}/amip5/amip5_8589_07_q1q2.nc ${DATA}/amip5/amip5_rgn_Pacific_Tropical_Western_8589_07_q1q2.nc
; ncwa -w gw -a lon,lat -d lat,-20.0,20.0 -d lon,130.0,180.0 ${DATA}/spcp_85/spcp_85_8589_07_q1q2.nc ${DATA}/spcp_85/spcp_85_rgn_Pacific_Tropical_Western_8589_07_q1q2.nc

;gcm_q1q2,cmp=1,fl_nm_1=getenv('DATA')+'/spcp_85/spcp_85_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc',fl_nm_2='/data2/zejder/amip5/amip5_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc'
;gcm_q1q2,cmp=1,fl_nm_1=getenv('DATA')+'/spcp_85/spcp_85_rgn_Pacific_Tropical_Western_8589_07_q1q2.nc',fl_nm_2=getenv('DATA')+'/amip5/amip5_rgn_Pacific_Tropical_Western_8589_07_q1q2.nc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_q1q2, $
	fl_nm_1=fl_nm_1, $
	fl_nm_2=fl_nm_2, $
	cmp=cmp

if n_elements(fl_nm_1) eq 0 then fl_nm_1=getenv('DATA')+'/spcp_85/spcp_85_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc'
if n_elements(fl_nm_2) eq 0 then fl_nm_2=getenv('DATA')+'/amip5/amip5_rgn_Pacific_Western_Warm_Pool_8589_07_q1q2.nc'
if n_elements(cmp) eq 0 then cmp=0
if n_elements(rng_y) eq 0 then rng_y=20.0

nc_id=ncdf_open(fl_nm_1)
ncdf_varget,nc_id,'Q1',q1
ncdf_varget,nc_id,'Q2',q2
ncdf_varget,nc_id,'QR',qr
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id

lev_nbr=n_elements(lev)

good_lev=where(lev ge 50.0)
lev=lev(good_lev)
q1=q1(good_lev)
q2=q2(good_lev)
qr=qr(good_lev)
q1mqr=q1-qr

abc_1=q1mqr*86400.0
;abc_1=q1*86400.0
abc_2=q2*86400.0
abc_3=qr*86400.0

ord=lev

abc_min=min([abc_1,abc_2,abc_3])
abc_max=max([abc_1,abc_2,abc_3])

plot, $
	abc_1, $
	ord, $
	tit='Mean Apparent Heating', $
	xtit='!E!12_!N!5K day!E-1!N', $
	ytit='!5Pressure !8p!5 (mb)', $
	xstyle=1, $
	ystyle=1, $
;	xrange=[abc_min,abc_max], $
	xrange=[-6.,4.], $
	yrange=[1000.0,100.0], $
	thick=3.0, $
	charsize=sym_sz, $
	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2.5], $  ;[4,2] is [bottom,top] default
	linestyle=0

oplot,	$
	abc_2, $
	ord, $
	thick=3.0, $
	linestyle=1

oplot,	$
	abc_3, $
	ord, $
	thick=3.0, $
	linestyle=2

if cmp then begin

nc_id=ncdf_open(fl_nm_2)
ncdf_varget,nc_id,'Q1',q1
ncdf_varget,nc_id,'Q2',q2
ncdf_varget,nc_id,'QR',qr
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id

lev_nbr=n_elements(lev)

good_lev=where(lev ge 50.0)
lev=lev(good_lev)
q1=q1(good_lev)
q2=q2(good_lev)
qr=qr(good_lev)
q1mqr=q1-qr

abc_1=q1mqr*86400.0
;abc_1=q1*86400.0
abc_2=q2*86400.0
abc_3=qr*86400.0

ord=lev

oplot,	$
	abc_1, $
	ord, $
	thick=3.0, $
	linestyle=0

oplot,	$
	abc_2, $
	ord, $
	thick=3.0, $
	linestyle=1

oplot,	$
	abc_3, $
	ord, $
	thick=3.0, $
	linestyle=2

endif; endif cmp

ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL

;xyouts,txt_lgn_x,lgn_y(0),'!8Q!5!I1!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!8Q!5!I1!N-!8Q!5!IR!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8Q!5!I2!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8Q!5!IR!N',size=txt_lgn_sz,/NORMAL

end; end gcm_q1q2()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM Q1Q2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM scat
; Scatterplot two fields (possibly from different files) against eachother
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_sct_bch, $
	dly=dly, $ ; [flg] Daily rather than monthly data
	mth=mth, $ ; [flg] Monthly data
	yrl=yrl, $ ; [flg] Yearly rather than monthly data
	dtrn=dtrn, $ ; [flg] Detrend data
	tm_sct=tm_sct, $ ; [flg] Time scatterplot
	psym=psym, $ ; [idx] Plot symbol
	lgr_flg=lgr_flg, $ ; [flg] Plot with logarithmic axes
	info=info, $
	prn=prn
if n_elements(dtrn) eq 0 then dtrn=0 ; [flg] Detrend data
if n_elements(tm_sct) eq 0 then tm_sct=1 ; [flg] Time scatterplot
if n_elements(psym) eq 0 then psym=1 ; [idx] Plot symbol
if n_elements(info) eq 0 then info=1
if n_elements(lgr_flg) eq 0 then lgr_flg=0 ; [flg] Plot with logarithmic axes
if n_elements(prn) eq 0 then prn=0
if n_elements(mth) eq 0 then mth=0 ; Monthly data
if n_elements(dly) eq 0 then dly=1 ; Daily data
if n_elements(yrl) eq 0 then yrl=0 ; Yearly data
if yrl then mth=0 ; Yearly data
if dly then mth=0 ; Daily data
if mth then dly=0 ; Monthly data

@ibp_clr.com

rgn_stf=[ $
;	['Pacific_Tropical_Western','Tropical Western Pacific'], $
;	['Pacific_Tropical_Central','Tropical Central Pacific'], $
;	['Pacific_Tropical_Eastern','Tropical Eastern Pacific'], $
;	['Indian_Tropical','Tropical Indian Ocean'], $
;	['Atlantic_Tropical','Tropical Atlantic'], $
;	['Pacific_Tropical','Tropical Pacific'], $
;	['Ocean_Tropical','Tropical Ocean']]
;	['Pacific_Equatorial_Western','Western Equatorial Pacific'], $
;	['Pacific_Equatorial_Central','Central Equatorial Pacific'], $
;	['Pacific_Equatorial_Central','Central Equatorial Pacific']]
;	['Pacific_Equatorial_Eastern','Eastern Equatorial Pacific'], $
;	['Pacific_Equatorial','Equatorial Pacific']]
	['Izn','Izania'], $
	['Brm','Bermuda'], $
	['Brb','Barbados']]
;	['tst','Test Region']]
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct_gcm,idx:0,nm:'',sng:''}
rgn=replicate({rgn_sct_gcm},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	rgn(idx).idx=idx
	rgn(idx).nm=rgn_stf(0,idx)
	rgn(idx).sng=rgn_stf(1,idx)
endfor; end loop over fld_lst

if not tm_sct then fld_cmp_stf=[ $
	['cnc_mss_dst','!5Observed !7l!5g m!E-3!N','DSTQ','!5Modeled !7l!5g m!E-3!N'] ]
;	['LWCF','!7d!5LWCF!5 (W m!E-2!N)','SWCF','!7d!5SWCF!5 (W m!E-2!N)']]
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','SWCF','!7d!5SWCF!5 (W m!E-2!N)'], $
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','LWCF','!7d!5LWCF!5 (W m!E-2!N)'] ]
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','TMQ','!7d!8Q!5!Iv!N (kg m!E-2!N)'], $
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','CMFMC','!7d!8M!5!Ic!N (g m!E-2!N s!E-1!N)']]
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','CLOUD','!7d!8A!5!Ic!N (%)']]
;	['TS1','!7d!5SST!5 (!E!12_!5!NK)','TS1','!7d!5SST!5 (!E!12_!5!NK)'] ]
if tm_sct then fld_cmp_stf=[ $
	['cnc_mss_dst','!5Observed !7l!5g m!E-3!N','DSTQ','!5Modeled !7l!5g m!E-3!N'] ]
;	['TS1','!5SST!5 (!E!12_!5!NK)','GCLD','!8G!5 (W m!E-2!N)'], $
;	['TS1','!5SST!5 (!E!12_!5!NK)','GCLR','!8G!5!Ia!N (W m!E-2!N)'], $
;	['TS1','!5SST!5 (!E!12_!5!NK)','LWCF','!5LWCF!5 (W m!E-2!N)'], $
;	['TS1','!5SST!5 (!E!12_!5!NK)','SWCF','!5SWCF!5 (W m!E-2!N)'], $
;	['LWCF','!5LWCF!5 (W m!E-2!N)','SWCF','!5SWCF!5 (W m!E-2!N)']]
fld_nbr_cmp=n_elements(fld_cmp_stf(0,*))
foo={fld_cmp_sct,idx:0,x_nm:'',x_ttl:'',y_nm:'',y_ttl:''}
fld_cmp=replicate({fld_cmp_sct},fld_nbr_cmp)
for idx=0,fld_nbr_cmp-1 do begin
	fld_cmp(idx).idx=idx
	fld_cmp(idx).x_nm=fld_cmp_stf(0,idx)
	fld_cmp(idx).x_ttl=fld_cmp_stf(1,idx)
	fld_cmp(idx).y_nm=fld_cmp_stf(2,idx)
	fld_cmp(idx).y_ttl=fld_cmp_stf(3,idx)
endfor; end loop over fld_lst

xpt_stf=[ $
;	['erbe_b','!5ERBE','erbe_b','!5ERBE']]
	['rsmas','!5U. Miami','dstmch60','!5dstmch60']]
;	['422','!5CCM2','422','!5CCM2'], $
;	['nflux18','!5NFLUX18','nflux18','!5NFLUX18'], $
;	['amip5','!5CCM!7X!5!I.5!N','amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM','amip5','!5CCM']]
;	['sld012d','!5CCM3','sld012d','!5CCM3']]
;	['spcp_85','!5ANV','spcp_85','!5ANV']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt2_sct,idx:0,x_nm:'',x_sng:'',y_nm:'',y_sng:''}
xpt=replicate({xpt2_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).x_nm=xpt_stf(0,idx)
	xpt(idx).x_sng=xpt_stf(1,idx)
	xpt(idx).y_nm=xpt_stf(2,idx)
	xpt(idx).y_sng=xpt_stf(3,idx)
endfor; end loop over fld_lst

; Year used by all types
;yr_srt=8589
yr_srt=1998
yr_nbr=1
yr=indgen(yr_nbr)
yr_sng=string(format='(I4.4)',yr_srt)

; Month used only if mth
mth_srt=3
mth_end=5

if mth then begin
; NB: Jan is mth_srt=1
; Next three commands for scatterplotting spatial regions (no detrending)
	yyyymmdd_srt_sng=string(format='(I4.4,I2.2)',yr_srt,mth_srt)
	yyyymmdd_end_sng=string(format='(I4.4,I2.2)',yr_end,mth_end)
	time_srt=mth_srt
	time_end=mth_end
endif else if dly then begin
	yyyymmdd_srt=19940901 ; LITE
	yyyymmdd_end=19940930 ; LITE
	yyyymmdd_srt=19990208 ; INDOEX
	yyyymmdd_end=19990322 ; INDOEX
	yyyymmdd_srt=19980101 ; 1998
	yyyymmdd_end=19981231 ; 1998
	yyyymmdd_srt=19950101 ; 
	yyyymmdd_end=19970729 ; 
	yyyymmdd_srt_sng=string(format='(I8.8)',yyyymmdd_srt)
	yyyymmdd_end_sng=string(format='(I8.8)',yyyymmdd_end)
	time_srt=yyyymmdd_srt
	time_end=yyyymmdd_end
endif else if yrl then begin
; NB: Store files with Fortran indexing notation
	yyyymmdd_srt=1994 ; LITE
endif ; endelse yrl
date_sng=yyyymmdd_srt_sng+'_'+yyyymmdd_end_sng

for fld_idx=0,fld_nbr_cmp-1 do begin
for rgn_idx=0,rgn_nbr-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	rng_x=0 ; Plotting routine decides limits
	rng_y=0 ; Plotting routine decides limits
	scl_x=1.0
	scl_y=1.0
	xpt_nm_y=xpt(xpt_idx).y_nm
	rgn_nm=rgn(rgn_idx).nm
	fld_nm_y=fld_cmp(fld_idx).y_nm
	fld_nm_x=fld_cmp(fld_idx).x_nm
	rgn_sng=rgn(rgn_idx).sng
	if dly then ttl=xpt(xpt_idx).x_sng+' '+rgn_sng
	if mth then begin
		if mth_srt eq 12 and mth_end eq 2 then ttl=xpt(xpt_idx).x_sng+' 1987-1985 DJF'
		if mth_srt eq 9 and mth_end eq 11 then ttl=xpt(xpt_idx).x_sng+' 1987-1985 ASO'
		if mth_srt eq 6 and mth_end eq 8 then ttl=xpt(xpt_idx).x_sng+' 1987-1985 JJA'
		if mth_srt eq 3 and mth_end eq 5 then ttl=xpt(xpt_idx).x_sng+' 1987-1985 MAM'
		if mth_srt eq 1 and mth_end eq 60 then ttl=xpt(xpt_idx).x_sng+' '+rgn_sng+' 1985-1989'
	endif; endif mth
	if not tm_sct then begin
		if strpos(fld_nm_y,'DSTQ') ne -1 then rng_y=[0.0,200.0]
		if strpos(fld_nm_y,'DSTQ') ne -1 then scl_y=scl_y*1.0e9
		if strpos(fld_nm_x,'cnc_mss_dst') ne -1 then rng_x=[0.0,200.0]
		if strpos(fld_nm_x,'cnc_mss_dst') ne -1 then scl_x=scl_x*1.0e9
		if strpos(fld_nm_x,'TS1') ne -1 then rng_x=[-1,3]
		if strpos(fld_nm_y,'TMQ') ne -1 then rng_y=[-30,30]
;		if strpos(fld_nm_y,'CMFMC') ne -1 then rng_y=[-5,20]
		if strpos(fld_nm_y,'CMFMC') ne -1 then scl_y=scl_y*1000
		if strpos(fld_nm_x,'LWCF') ne -1 then rng_x=[-100,100]
		if strpos(fld_nm_y,'LWCF') ne -1 then rng_y=[-100,100]
		if strpos(fld_nm_x,'SWCF') ne -1 then rng_x=[-100,100]
		if strpos(fld_nm_y,'SWCF') ne -1 then rng_y=[-100,100]
	endif; endif not tm_sct
	if tm_sct then begin
; fxm: scale dust mixing ratios into mass concentrations as in gcm_tv_bch
		if strpos(xpt_nm_y,'dst') ne -1 and strpos(fld_nm_y,'DSTQ') ne -1 and strpos(rgn_nm,'Brb') ne -1 then scl_y=scl_y*1.15
		if strpos(xpt_nm_y,'dst') ne -1 and strpos(fld_nm_y,'DSTQ') ne -1 and strpos(rgn_nm,'Kaa') ne -1 then scl_y=scl_y*1.15
		if strpos(xpt_nm_y,'dst') ne -1 and strpos(fld_nm_y,'DSTQ') ne -1 and strpos(rgn_nm,'Brm') ne -1 then scl_y=scl_y*1.15
		if strpos(xpt_nm_y,'dst') ne -1 and strpos(fld_nm_y,'DSTQ') ne -1 and strpos(rgn_nm,'Izn') ne -1 then scl_y=scl_y*0.9
;		if strpos(fld_nm_y,'DSTQ') ne -1 then rng_y=[0.0,200.0]
		if strpos(fld_nm_y,'DSTQ') ne -1 then scl_y=scl_y*1.0e9
;		if strpos(fld_nm_x,'cnc_mss_dst') ne -1 then rng_x=[0.0,200.0]
		if strpos(fld_nm_x,'cnc_mss_dst') ne -1 then scl_x=scl_x*1.0e9
		if strpos(fld_nm_x,'TS1') ne -1 then rng_x=[-1,1]
		if strpos(fld_nm_y,'GCLD') ne -1 then rng_y=[-10,10]
		if strpos(fld_nm_y,'GCLR') ne -1 then rng_y=[-10,10]
		if strpos(fld_nm_x,'LWCF') ne -1 then rng_x=[-10,10]
		if strpos(fld_nm_y,'LWCF') ne -1 then rng_y=[-10,10]
		if strpos(fld_nm_x,'SWCF') ne -1 then rng_x=[-10,10]
		if strpos(fld_nm_y,'SWCF') ne -1 then rng_y=[-10,10]
	endif; endif tm_sct
	fl_nm_in=strarr(2)
	if tm_sct then begin
		trn=1
		fl_nm_in(0)=getenv('DATA')+'/'+xpt(xpt_idx).x_nm+'/'+xpt(xpt_idx).x_nm+'_'+date_sng+'_'+rgn(rgn_idx).nm+'.nc'
		fl_nm_in(1)=getenv('DATA')+'/'+xpt(xpt_idx).y_nm+'/'+xpt(xpt_idx).y_nm+'_'+date_sng+'_'+rgn(rgn_idx).nm+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).x_nm+'_'+xpt(xpt_idx).y_nm+'_'+date_sng+rgn(rgn_idx).nm+'_'+fld_nm_x+'_'+fld_nm_y+'.eps'
	endif else begin; endif tm_sct
		trn=0
		fl_nm_in(0)=getenv('DATA')+'/'+xpt(xpt_idx).x_nm+'/'+xpt(xpt_idx).x_nm+'_'+date_sng+rgn(rgn_idx).nm+'.nc'
		fl_nm_in(1)=getenv('DATA')+'/'+xpt(xpt_idx).y_nm+'/'+xpt(xpt_idx).y_nm+'_'+date_sng+rgn(rgn_idx).nm+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).x_nm+'_'+xpt(xpt_idx).y_nm+'_'+date_sng+rgn(rgn_idx).nm+'_'+fld_nm_x+'_'+fld_nm_y+'.eps'
	endelse; endif not tm_sct 
	if prn then open_ps,fl_nm=fl_nm_out,/one_clm,/eps
	gcm_sct,fld_x=fld_nm_x,fld_y=fld_nm_y,x_ttl=fld_cmp(fld_idx).x_ttl,y_ttl=fld_cmp(fld_idx).y_ttl,rng_x=rng_x,rng_y=rng_y,scl_x=scl_x,scl_y=scl_y,time_srt=time_srt,time_end=time_end,fl_nm=fl_nm_in,ttl=ttl,trn=trn,tm_sct=tm_sct,info=info,dtrn=dtrn,psym=psym,lgr_flg=lgr_flg,prn=prn

	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over experiments
endfor; end loop over rgn_lst
endfor; end loop over fld_lst

end; end gcm_sct_bch()
pro gcm_sct, $
	fl_nm=fl_nm, $
	ttl=ttl, $
	time_srt=time_srt, $
	time_end=time_end, $
	scl_x=scl_x, $
	scl_y=scl_y, $
	fld_x=fld_x, $
	fld_y=fld_y, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	trn=trn, $
	tm_sct=tm_sct, $
	dtrn=dtrn, $ ; [flg] Detrend data
	info=info, $
	prn=prn, $
	fl_nbr=fl_nbr, $
	dbg_lvl=dbg_lvl, $ ; [idx] Debugging level
	dly=dly, $ ; [flg] Daily rather than monthly data
	mth=mth, $ ; [flg] Monthly data
	yrl=yrl, $ ; [flg] Yearly rather than monthly data
	psym=psym, $ ; [idx] Plot symbol
	lgr_flg=lgr_flg, $ ; [flg] Plot with logarithmic axes
	scl=scl

if n_elements(fl_nm) eq 0 then fl_nm=[getenv('DATA')+'/dstmch60/dstmch60_1998.nc',getenv('DATA')+'/dstmch60/dstmch60_1998.nc']
if n_elements(dtrn) eq 0 then dtrn=0 ; [flg] Detrend data
if n_elements(time_srt) eq 0 then time_srt=5
if n_elements(time_end) eq 0 then time_end=7
if n_elements(scl_x) eq 0 then scl_x=1
if n_elements(psym) eq 0 then psym=1 ; [idx] Plot symbol
if n_elements(scl_y) eq 0 then scl_y=1
if n_elements(fld_x) eq 0 then fld_x='LWCF'
if n_elements(fld_y) eq 0 then fld_y='SWCF'
if n_elements(x_ttl) eq 0 then x_ttl='!7d!5LWCF!5 (W m!E-2!N)'
if n_elements(y_ttl) eq 0 then y_ttl='!7d!5SWCF!5 (W m!E-2!N)'
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='1987 - 1985 ANV Eastern Equatorial Pacific'
if n_elements(trn) eq 0 then trn=0
if n_elements(tm_sct) eq 0 then tm_sct=0
if n_elements(info) eq 0 then info=1
if n_elements(lgr_flg) eq 0 then lgr_flg=0 ; [flg] Plot with logarithmic axes
if n_elements(prn) eq 0 then prn=0
if n_elements(fl_nbr) eq 0 then fl_nbr=1
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0 ; Debugging level
if n_elements(mth) eq 0 then mth=0 ; Monthly data
if n_elements(dly) eq 0 then dly=1 ; Daily data
if n_elements(yrl) eq 0 then yrl=0 ; Yearly data

erase
!p.multi=0

print,'Processing '+fl_nm(0)+'...'
nc_id=ncdf_open(fl_nm(0))
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,fld_x,x_data
ncdf_close,nc_id
print,'Processing '+fl_nm(1)+'...'
nc_id=ncdf_open(fl_nm(1))
ncdf_varget,nc_id,fld_y,y_data
ncdf_close,nc_id

if scl_x ne 1.0 then x_data=scl_x*x_data
if scl_y ne 1.0 then y_data=scl_y*y_data
if not tm_sct then begin
	if time_nbr ne 12 then print,'gcm_sct(): WARNING does not have 12 timesteps'
	; fxm: this does not work for DJF as implemented
	if mth then begin
		x_data=reform(x_data(*,*,time_srt:time_end))
		y_data=reform(y_data(*,*,time_srt:time_end))
	endif ; endif mth
endif ; endif not tm_sct
if tm_sct then begin
	if strpos(fl_nm(0),'erbe') ne -1 and time_end ne 60 then print,'WARNING: ERBE time sequence file does not have 60 timesteps'
	if mth then begin
		x_data=reform(x_data(time_srt:time_end))
		y_data=reform(y_data(time_srt:time_end))
	endif; endif mth
endif ; endif tm_sct

if strpos(fl_nm(0),'erbe') ne -1 and time_end eq 59 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
x_data=x_data(1:52)
y_data=y_data(1:52)
endif ; endif erbe

print,'gcm_sct(): Started with ',n_elements(x_data),' points'
good_idx=where(x_data lt 1.0e20,good_data_nbr)
x_data=x_data(good_idx)
y_data=y_data(good_idx)
good_idx=where(y_data lt 1.0e20,good_data_nbr)
x_data=x_data(good_idx)
y_data=y_data(good_idx)
print,'gcm_sct(): Ended with ',n_elements(x_data),' points'

; Detrend timeseries data
if dtrn then begin

	; Get trend from least squares fit to data
	abc=indgen(good_data_nbr)
	fit_cff=poly_fit(abc,x_data,1,ord_fit)
	slope=fit_cff(1)
	x_data=x_data-slope*(abc-abc(0)-0.5*(abc(good_data_nbr-1)-abc(0)))
	fit_cff=poly_fit(abc,y_data,1,ord_fit)
	slope=fit_cff(1)
	y_data=y_data-slope*(abc-abc(0)-0.5*(abc(good_data_nbr-1)-abc(0)))

	;print,'slope before detrending = ',slope
	;fit_cff=poly_fit(abc,y_data,1,ord_fit)
	;print,'slope after detrending = ',fit_cff(1)

endif; endif dtrn

min_x_data=min(x_data)
max_x_data=max(x_data)
min_y_data=min(y_data)
max_y_data=max(y_data)

if lgr_flg then begin
	if n_elements(rng_x) le 1 then rng_x=[0.1,max_x_data] else rng_x(0)=0.1
	if n_elements(rng_y) le 1 then rng_y=[0.1,max_y_data] else rng_y(0)=0.1
endif else begin
	if n_elements(rng_x) le 1 then rng_x=[min_x_data,max_x_data]
	if n_elements(rng_y) le 1 then rng_y=[min_y_data,max_y_data]
endelse ; endif lgr_flg
mdl_obs=1 ; Plot abscissa on same scale as ordinate ("model vs. observation")
if mdl_obs then begin
; mdl_obs means plotting model against observations so axes should be identical
	rng_x=rng_y
endif ; endif mdl_obs

chr_sz=1.5
print,'rng_x = ',rng_x
print,'rng_y = ',rng_y
print,'x_data = ',x_data
print,'y_data = ',y_data
plot, $
	x_data, $
	y_data, $
	xlog=lgr_flg, $
	ylog=lgr_flg, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xrange=rng_x, $
	yrange=rng_y, $
	xstyle=1, $
	ystyle=1, $
	thick=2.0, $
	charsize=chr_sz, $
	xmargin=[8.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2], $  ;[4,2] is [bottom,top] default
	color=clr_blk_idx, $
	psym=1

; Perform statistical analysis of correlation between datasets
if info then begin
fit_cff=poly_fit(x_data,y_data,1,ord_fit)

data=fltarr(2,n_elements(x_data))
data(0,*)=x_data
data(1,*)=y_data
pearson_mtx=correl_matrix(data)

oplot, $
	[min_x_data, $
	max_x_data], $
	[fit_cff(0)+min_x_data*fit_cff(1), $
	fit_cff(0)+max_x_data*fit_cff(1)], $
	max_value=1.0e20, $
	thick=3.0, $
	linestyle=0

m_sng=auto_sng(fit_cff(1),2)
pearson_sng=string(format='(F5.2)',pearson_mtx(0,1))
pearson_sng='!8r!5 = '+pearson_sng
if fld_y eq 'SWCF' then x_psn=0.5 else x_psn=.25
xyouts,x_psn,0.85,'Correlation '+pearson_sng,size=1.0*chr_sz,/NORMAL
xyouts,x_psn,0.75,'Slope = '+m_sng,size=1.0*chr_sz,/NORMAL

data=fltarr(1,n_elements(x_data))
data(0,*)=x_data
wgt=fltarr(n_elements(x_data))+1.0
coeff=regress(data,y_data,wgt,/relative_weight,YFIT, A0, SIGMA, FTEST, R, RMUL, CHISQ)
; Taylor train crash book p. 226. assume x_data has negligible uncertainty,
; and that uncertainty in y_data for SWCF 1987-1985 is 10 w/m2 - 10 W/m2 
; (errors in inidividual gridpoint monthly values) added in quadrature ~ 5 W/m2
my_chi=(total((yfit-y_data)^2))/5^2
if dbg_lvl gt 0 then begin
	print,"gcm_sct(): Statistics: "
	print,"coeff = ",coeff
	print,"a0 = ",a0
	print,"sigma = ",sigma
	print,"ftest = ",ftest
	print,"r = ",r
	print,"rmul = ",rmul
	print,"chisq = ",chisq
	print,"my_chi = ",my_chi
endif; endif dbg
endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif; endif not prn

end; end gcm_sct()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM scat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM YZ
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_yz_bch, $
	ams=ams, $ ; plot fields suitable for ams
	bias=bias, $ ; plot biases of fields relative to observational analyses
	diff=diff, $ ; plot differences between model fields
	prs=prs, $ ; insert 'prs_' in filename (use files on pressure levels)
	csn=csn, $ ; Rearrange headers and footers for seasonal four panel plots
	ss=ss, $ ; Source/sink climatology
	bgr=bgr, $ ; background color?
	pll=pll, $ ; color palette
	pnl_lbl=pnl_lbl, $ ; Panel label
	cbar=cbar, $ ; Add color scale
	cnt_ovl=cnt_ovl, $ ; Add contour overlays
	prn=prn

; Dust annual averages:
; gcm_yz_bch,diff=0,csn=0,bgr=1,prn=0,cnt_ovl=0,pll=4
; Dust response:
; gcm_yz_bch,diff=1,csn=0,bgr=1,prn=0,cnt_ovl=0,pll=4
; Dimer seasonal averages:
; gcm_yz_bch,diff=0,csn=1,bgr=1,prn=0,cnt_ovl=1,pll=4
; Example of "shading":
; gcm_yz_bch,diff=1,csn=1,bgr=0,prn=0,pll=1,cnt_ovl=1
; Nice plots of seasonal T,QRS responses
; gcm_yz_bch,diff=1,csn=1,bgr=0,prn=1,pll=6,cnt_ovl=1

if n_elements(prn) eq 0 then prn=0
if n_elements(ams) eq 0 then ams=0
if n_elements(bias) eq 0 then bias=0
if n_elements(diff) eq 0 then diff=0
if n_elements(prs) eq 0 then prs=0
if n_elements(csn) eq 0 then csn=0
if n_elements(ss) eq 0 then ss=0
if n_elements(bgr) eq 0 then bgr=0
if n_elements(pll) eq 0 then pll=3
if n_elements(cbar) eq 0 then cbar=0
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if prs then prs_sng='prs_' else prs_sng=''
if n_elements(pnl_lbl) eq 0 then pnl_lbl=''

fld2_stf=[ $
;	['MPSI','!5Meridional Streamfunction','!9X!5 10!E10!N kg s!E-1!N'] ]
;	['OMEGA','!5Vertical Velocity','!5mb day!E-1!N'] ]
;	['U','!5Zonal Wind','!5m s!E-1!N'] ]
;	['CMFMC','!5Convective Mass Flux','!8M!5!IC!N g m!E-2!N s!E-1!N'] ]
;	['QL','!5Liquid Mixing Ratio','!5mg kg!E-1!N'] ]
;	['QICE','!5Ice Mixing Ratio','!5mg kg!E-1!N'] ]
;	['QC','!5Condensate Mixing Ratio','!5mg kg!E-1!N'] ]
;	['QDIABAT','!5Total Diabatic Heating','!5K day!E-1!N'] ]
;	['CMFDT','!5Convective Heating','!5K day!E-1!N'], $
;	['HGS','!5Large Scale Heating','!5K day!E-1!N'], $
;	['RADD','!5Net Radiative Heating','!5K day!E-1!N'] ]
;	['dust_ttl','!5Mineral Dust','!7l!5g kg!E-1!N'] ]
;	['DSTQ01','!5Clay Mixing Ratio !8Q!I!5d01!N!5','!7l!5g kg!E-1!N']]
;	['DSTQ02','!5Fine Silt Mixing Ratio !8Q!I!5d02!N!5','!7l!5g kg!E-1!N'], $
;	['DSTQ03','!5Coarse Silt Mixing Ratio !8Q!I!5d03!N!5','!7l!5g kg!E-1!N'], $
;	['DSTQ04','!5Sand Mixing Ratio !8Q!I!5d04!N!5','!7l!5g kg!E-1!N'], $
	['DSTQ','!5Dust Mixing Ratio !8Q!I!5d!N!5','!7l!5g kg!E-1!N'], $
	['DSTSSPCP','!5Washout','!5pg kg!E-1!N s!E-1!N'], $
	['DSTSSDRY','!5Settling and Mixout','!5pg kg!E-1!N s!E-1!N'], $
	['DSTSSEVP','!5Evaporation','!5pg kg!E-1!N s!E-1!N'] ]
;	['DSTSSGRV','!5Gravitational Settling','!5pg kg!E-1!N s!E-1!N'], $
;	['DSTSSMBL','!5Mobilization','!5pg kg!E-1!N s!E-1!N'], $
;	['DSTSSTRB','!5Turbulent Deposition','!5pg kg!E-1!N s!E-1!N'] ]
;	['Q','!5Vapor Mixing Ratio','!5g kg!E-1!N'], $
;	['DTCOND','!5Condensation Heating','!5K day!E-1!N'], $
;	['CLOUD','!5Cloud Amount','!5%'] ]
;	['QRL','!5Longwave Heating','!5K day!E-1!N'], $
;	['QRS','!5Solar Heating','!5K day!E-1!N'], $
;	['T','!5Temperature','!5!NK'] ]
;	['QRSFRC','!5Tracer Solar Heating','!5K day!E-1!N'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
;	['dmr03','!5']]
;	['dmr04','!5']]
;	['dmr05','!5']]
;	['dmr06','!5']]
;	['dmr07','!5']]
;	['dmr08','!5']]
;	['obsst01','!5CCM']]
;	['ncep','!5NCEP/NCAR']]
;	['dstmch46','!5dstmch46']]
	['dstmch90','!5dstmch90']]
;	['dstccm45','!5dstccm45']]
;	['ecmwf','!5ECMWF'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'] ]
;	['sld012d','!5CCM3'], $
;	['spcp_85','!5ANV']]
;	['tef95','!5Tegen & Fung 1995']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

;mth_sng=['01','07']
;mth_sng=['01','02','03','04','05','06','07','08','09','10','11','12']
;mth_sng=['12','01','02']
;mth_sng=['09']
;mth_sng=['01']
;mth_sng=['1202','0305','0608','0911']
;mth_sng=['1202']
;mth_sng=['0305']
;mth_sng=['0608']
;mth_sng=['0911']
;mth_sng=['1202','0608']
;mth_sng=['1202','0305']
;mth_sng=['0608','0911']
;mth_sng=['8589']
mth_sng=['clm']
mth_nbr=n_elements(mth_sng)
for mth_idx=0,mth_nbr-1 do begin
for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
;	Replace chemical names, if necessary
	fld2(fld_idx2).sng=trc_nm_sbs(fld2(fld_idx2).sng,xpt(xpt_idx).nm)
	ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
;	if prs then y_ttl='!5Pressure (mb)' else y_ttl='!7g!5 !9X!5 1000 (mb)'
	y_ttl='!5Pressure (mb)'
	x_ttl=''
	ctl_xpt='dstccm45_'
	ctl_sng=' - dstccm45'
;	ctl_sng=''
	rng_x=[-90,90]
;	rng_y=[0,1000]
;	rng_x=[0,30] ; Saharan dust layer
	rng_y=[500,1000] ; Saharan dust layer
	pre_scl=0
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	dsh=0
	yr_sng=''
;	yr_sng='clm_'
;	yr_sng='1998_'
;	yr_sng='7993_'
	shd_val=0.0
	shd_xpr='lt'
	if fld2(fld_idx2).nm eq 'MPSI' then scl=1.0e-10
	if fld2(fld_idx2).nm eq 'MPSI' then if not diff then pre_scl=[-32,2,30] else pre_scl=[-10,1,21]
	if fld2(fld_idx2).nm eq 'CMFMC' then scl=1000
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then pre_scl=[0,2,30] else pre_scl=[-5,0.5,30]
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'OMEGA' then scl=864
	if fld2(fld_idx2).nm eq 'OMEGA' then pre_scl=[-75,5,30]
	if fld2(fld_idx2).nm eq 'U' then if not diff then pre_scl=[-30,5,18] else pre_scl=[-10,1,21]
;	if fld2(fld_idx2).nm eq 'T' then if not diff then pre_scl=[180,5,25] else pre_scl=[-10,1,21]
	if fld2(fld_idx2).nm eq 'T' then if not diff then pre_scl=[180,5,25] else pre_scl=[-5,0.5,21]
	if fld2(fld_idx2).nm eq 'Q' then scl=1000.0
	if fld2(fld_idx2).nm eq 'Q' then if diff then pre_scl=[-0.5,0.05,21] else pre_scl=[0,1,30]
	if fld2(fld_idx2).nm eq 'QC' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QC' then if diff then pre_scl=[-34,2,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QICE' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QICE' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QL' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QL' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QL' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QL' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QDIABAT' then scl=86400
	if fld2(fld_idx2).nm eq 'QDIABAT' then pre_scl=[-1.5,0.1,30]
	if fld2(fld_idx2).nm eq 'QRS' then if not diff then shd_val=6 else shd_val=0.0
	if fld2(fld_idx2).nm eq 'QRS' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QRS' then scl=86400.0
	if fld2(fld_idx2).nm eq 'QRS' then if diff then pre_scl=[-0.2,0.02,21] else pre_scl=[0.0,0.1,12]
	if fld2(fld_idx2).nm eq 'QRSFRC' then if not diff then shd_val=6 else shd_val=0.0
	if fld2(fld_idx2).nm eq 'QRSFRC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QRSFRC' then scl=86400.0
	if fld2(fld_idx2).nm eq 'QRSFRC' then if diff then pre_scl=[-0.2,0.02,21] else pre_scl=[0.0,0.02,11]
	if fld2(fld_idx2).nm eq 'QRL' then scl=86400
	if fld2(fld_idx2).nm eq 'QRL' then if not diff then pre_scl=[-4.,0.2,30] else pre_scl=[-1.0,0.1,20]
;	if fld2(fld_idx2).nm eq 'RADD' then pre_scl=[-0.5,.05,30]
	if fld2(fld_idx2).nm eq 'RADD' then pre_scl=[-1.5,0.1,30]
	if fld2(fld_idx2).nm eq 'DTCOND' then scl=86400
	if fld2(fld_idx2).nm eq 'DTCOND' then if not diff then pre_scl=[0.0,0.1,20] else pre_scl=[-1.0,0.1,20]
	if fld2(fld_idx2).nm eq 'CMFDT' then scl=86400
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then pre_scl=[-1.5,0.1,30] else pre_scl=[-1.5,0.1,30]
	if fld2(fld_idx2).nm eq 'HGS' then scl=86400
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then pre_scl=[-1.5,0.1,30] else pre_scl=[-1.5,0.1,30]
	if fld2(fld_idx2).nm eq 'CLOUD' then scl=100
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then pre_scl=[0,5,21] else pre_scl=[-20,2,20]
	if fld2(fld_idx2).nm eq 'dust_ttl' then scl=1.0e9
	if strpos(fld2(fld_idx2).nm,'DSTQ') ne -1 then begin
		scl=1.0e9
		if fld2(fld_idx2).nm eq 'DSTQ' then pre_scl=[0,1,11]
		if fld2(fld_idx2).nm eq 'DSTQ01' then pre_scl=[0,1,11]
		if fld2(fld_idx2).nm eq 'DSTQ02' then pre_scl=[0,1,11]
		if fld2(fld_idx2).nm eq 'DSTQ03' then pre_scl=[0,1,11]
		if fld2(fld_idx2).nm eq 'DSTQ04' then pre_scl=[0,1,11]
	endif ; endif dust
	if strpos(fld2(fld_idx2).nm,'DSTSS') ne -1 then begin
		xpt(xpt_idx).sng=''
		scl=1.0e15
		if fld2(fld_idx2).nm eq 'DSTSSEVP' then pre_scl=[0,0.002,11]
		if fld2(fld_idx2).nm eq 'DSTSSGRV' then pre_scl=[-10,1,11]
		if fld2(fld_idx2).nm eq 'DSTSSDRY' then pre_scl=[-10,1,11]
		if fld2(fld_idx2).nm eq 'DSTSSMBL' then pre_scl=[0,1,11]
		if fld2(fld_idx2).nm eq 'DSTSSPCP' then pre_scl=[0,0.03,11]
		if fld2(fld_idx2).nm eq 'DSTSSTRB' then pre_scl=[0,3,11]
	endif ; endif dust
 	if xpt(xpt_idx).nm eq 'dst25' then begin
		if strpos(fld2(fld_idx2).nm,'DST') ne -1 then begin
			scl=scl*6.0
		endif ; endif DSTMPC or DSTODX or DSTSF or DSTSS or DSTQ
	endif ; endif dst25
	if strpos(fld2(fld_idx2).nm,'FRC') ne -1 or strpos(fld2(fld_idx2).nm,'DMR') ne -1 then begin
		if xpt(xpt_idx).nm eq 'dmr03' then begin
			if fld2(fld_idx2).nm eq 'QRSFRC' then scl=86400.0
			if fld2(fld_idx2).nm eq 'QRSFRC' then pre_scl=[0,0.01,11]
			rng_y=[500,1000]
		endif; not dmr03 H2OH2O
		if xpt(xpt_idx).nm eq 'dmr04' or xpt(xpt_idx).nm eq 'dmr05' or xpt(xpt_idx).nm eq 'dmr06' then begin
			if fld2(fld_idx2).nm eq 'QRSFRC' then scl=86400.0
			if fld2(fld_idx2).nm eq 'QRSFRC' then pre_scl=[0,0.002,15]
			rng_y=[0,1000]
		endif; not dmr04 H2OH2O
	endif ; endif dimer
	if strpos(xpt(xpt_idx).nm,'ecmwf') ne -1 then begin
		if mth_sng(mth_idx) eq '01' then yr_sng='9095_' else yr_sng='8994_'
	endif
;	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+'xavg_'+yr_sng+mth_sng(mth_idx)+'.nc'
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+yr_sng+mth_sng(mth_idx)+'_x.nc'
;	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+prs_sng+'xavg_'+yr_sng+mth_sng(mth_idx)+'_'+fld2(fld_idx2).nm+'.eps'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+prs_sng+yr_sng+mth_sng(mth_idx)+'_x_'+fld2(fld_idx2).nm+'.eps'
	if diff then begin ; plot differences between simulated fields
		ttl=xpt(xpt_idx).sng+ctl_sng+' ('+fld2(fld_idx2).unit+')'
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+ctl_xpt+prs_sng+yr_sng+mth_sng(mth_idx)+'_x.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+ctl_xpt+prs_sng+yr_sng+mth_sng(mth_idx)+'_x_'+fld2(fld_idx2).nm+'.eps'
	endif; endif diff
	if prn and not bias and xpt_nbr eq 1 then begin
		ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		if diff then ttl=xpt(xpt_idx).sng+ctl_sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
;		ttl=fld2(fld_idx2).sng
		if mth_idx eq 0 then top=1 else top=0
		if mth_idx eq mth_nbr-1 then btm=1 else btm=0
	endif; endif prn not bias
	if bias then begin ; plot biases of simulated fields relative to observational analyses
	if strpos(xpt(xpt_idx).nm,'amip5') ne -1 or strpos(xpt(xpt_idx).nm,'spcp') ne -1 then begin
		ttl=xpt(xpt_idx).sng+'-ECMWF '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		if mth_sng(mth_idx) eq '01' then ecmwf_yr_sng='9095_'
		if mth_sng(mth_idx) eq '07' then ecmwf_yr_sng='8994_'
		if fld2(fld_idx2).nm eq 'MPSI' then pre_scl=[-15,1,30]
		if fld2(fld_idx2).nm eq 'U' then pre_scl=[-20,2,20]
		if fld2(fld_idx2).nm eq 'T' then pre_scl=[-15,1,30]
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+yr_sng+'ecmwf_'+ecmwf_yr_sng+prs_sng+'xavg_'+mth_sng(mth_idx)+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+yr_sng+'ecmwf_'+ecmwf_yr_sng+prs_sng+'xavg_'+mth_sng(mth_idx)+'_'+fld2(fld_idx2).nm+'.eps'
	endif; endif xpt is a simulation, not ecmwf
	endif; endif bias
	if (prn and bias) or (prn and not bias and xpt_nbr gt 1) then begin
		if mth_sng(mth_idx) eq '01' then ttl='Jan '+fld2(fld_idx2).sng
		if mth_sng(mth_idx) eq '07' then ttl='Jul '+fld2(fld_idx2).sng
		if xpt_idx eq 0 then top=1 else top=0
		if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
	endif; endif prn and bias
	if prn then begin
		if ss then begin
			top=1
			if fld2(fld_idx2).nm eq 'DSTSSMBL' or fld2(fld_idx2).nm eq 'DSTSSEVP' then btm=0
			if fld2(fld_idx2).nm eq 'DSTSSPCP' or fld2(fld_idx2).nm eq 'DSTSSEVP' then btm=1
		endif; endif ss
		if mth_nbr eq 4 and csn then begin
			if (mth_sng(mth_idx) eq '1202' or mth_sng(mth_idx) eq '0608') then begin
				top=1
				btm=0
			endif; endif top panel
			if (mth_sng(mth_idx) eq '0305' or mth_sng(mth_idx) eq '0911') then begin
				top=0
				btm=1
			endif; endif top panel
		endif; endif csn
	endif; endif prn
	if pnl_lbl ne '' then pnl_lbl=pnl_lbl_get(fld_idx2)
	if ams then begin
		rng_x=[-30,30]
	endif; endif ams
	if prn then begin
;		top=0
;		btm=1
		x_sz=6.5
;		y_sz=3.2
		y_sz=3.8
;		y_sz=5.0
		if top then y_sz=y_sz*1.08
		if btm then y_sz=y_sz*1.08
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_yz,x_ttl=x_ttl,y_ttl=y_ttl,fl_nm=fl_nm_in,ttl=ttl,prn=prn,bgr=bgr,pll=pll,order=1,cbar=cbar,cnt_ovl=cnt_ovl,fld_nm=fld2(fld_idx2).nm,wrp=0,rng_x=rng_x,rng_y=rng_y,pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,top=top,btm=btm,scl=scl,shd_val=shd_val,shd_xpr=shd_xpr,pnl_lbl=pnl_lbl
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over experiment
endfor; end loop over fld2
endfor; end loop over mth

end; end gcm_yz_bch()
pro gcm_yz, $
	fl_nm=fl_nm, $
	order=order, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	shd_val=shd_val, $
	shd_xpr=shd_xpr, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	dsh=dsh, $
	wrp=wrp, $
	bgr=bgr, $
	pll=pll, $
	cbar=cbar, $
	cnt_ovl=cnt_ovl, $
	scl=scl, $
	pre_scl=pre_scl, $
	pnl_lbl=pnl_lbl, $
	fld_nm=fld_nm

; gcm_yz,ttl='!5ANV Zonal Wind !8U!5 (m s!E-1!N)',pll=4,order=1
; gcm_yz,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_8589_amip5_8589_xavg_01.nc',fld_nm='RADD',scl=86400,ttl='!5ANV-CCM !5Net Radiative Heating !8Q!5!IR!N (!5K day!E-1!N)',pll=4,order=1
; gcm_yz,fl_nm=getenv('DATA')+'/amip5/amip5_xavg_8589_01_mmr.nc',fld_nm='QC',scl=1.0e6,ttl='!5CCM !5Condensate Mixing Ratio !5q!Ic!N (!5g kg!E-1!N)',pll=1,order=1,rng_x9[-90,90],cbar=0
; gcm_yz,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xavg_8589_01_mmr.nc',fld_nm='QC',scl=1.0e6,ttl='!5ANV !5Condensate Mixing Ratio !5q!Ic!N (!5g kg!E-1!N)',pll=1,order=1,rng_x=[-90,90],cbar=0
; gcm_yz,fl_nm=getenv('DATA')+'/amip5/amip5_xavg_8589_01_mmr.nc',fld_nm='QICE',scl=1.0e6,ttl='!5CCM !5Ice Mixing Ratio !5q!Ii!N (!5g kg!E-1!N)',pll=1,order=1,rng_x=[-90,90],cbar=0
; gcm_yz,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xavg_8589_01_mmr.nc',fld_nm='QICE',scl=1.0e6,ttl='!5ANV !5Ice Mixing Ratio !5q!Ii!N (!5g kg!E-1!N)',pll=1,order=1,rng_x=[-90,90],cbar=0
; gcm_yz,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_prs_xavg_8589_01_MPSI.nc',fld_nm='MPSI',scl=1.0e-10,ttl='!5ANV !5Meridional Stream Function !7w!5!IM!N (!9X!5 10!E10!N kg s!E-1!N)',pll=4,order=1,pre_scl=[-32,2,30]
; gcm_yz,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_prs_xavg_8589_01.nc',fld_nm='MPSI',scl=1.0e-10,ttl='!5ANV !5Meridional Stream Function !7w!5!IM!N (!9X!5 10!E10!N kg s!E-1!N)',pll=4,order=1

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(wrp) eq 0 then wrp=0
if n_elements(shd_val) eq 0 then shd_val=0.0
if n_elements(shd_xpr) eq 0 then shd_xpr='lt'
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ANV Zonal Wind !8U!5 (m s!E-1!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl='!5Prssure (mb)'
if n_elements(x_ttl) eq 0 then x_ttl='!5Latitude'
if n_elements(order) ne 0 then clr_ord=order else clr_ord=1
if n_elements(bgr) eq 0 then bgr=0
if n_elements(pll) eq 0 then pll=3
if n_elements(cbar) eq 0 then cbar=1
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xavg_8589_01.nc'
if n_elements(fld_nm) eq 0 then fld_nm='U'
if n_elements(pnl_lbl) eq 0 then pnl_lbl=''

if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=34
if n_elements(color_order) eq 0 then color_order=clr_ord
clr_rfr,clr_tbl_nbr,"",0

erase
!p.multi=0

; Read in binary data from the NetCDF file 
print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,foo,lat_nbr
lev_nm='lev'
;if fld_nm eq 'CMFMC' or fld_nm eq 'MPSI' then if strpos(fl_nm,'prs') ne -1 then lev_nm='ilev'
lev_id=ncdf_dimid(nc_id,lev_nm)
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
ncdf_varget,nc_id,lev_nm,lev
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,fld_nm,data
ncdf_close,nc_id
; End of NetCDF commands

; fxm: 20000628 only transpose if data are in old CCM order
;data=transpose(reform(data)) ; tranposing makes contour() put lat on x-axis, lev on y-axis
good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=lat
ord=lev

data_min=min(data)
data_max=max(data)

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

if auto_scl then begin

; Automatically set and scale the contour levels
cntr_lvl_nbr=20

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

if fld_nm eq 'QDIABAT' then begin
	cntr_lvl(0)=-10
	cntr_lvl(cntr_lvl_nbr-1)=10
endif; end QDIABAT
if fld_nm eq 'RADD' then begin
	cntr_lvl(0)=-10
	cntr_lvl(cntr_lvl_nbr-1)=10
endif; end RADD
if fld_nm eq 'QC' or fld_nm eq 'QL' or fld_nm eq 'QICE' then begin
;	cntr_lvl_nbr=19
;	cntr_ntv=10
;	cntr_lvl=[0.0,0.5,1,2,4,6,8,10,15,20,25,30,40,50,60,70,80,90,100]
;	cntr_lvl_nbr=27
;	cntr_ntv=10
;	cntr_lvl=[-40,-30,-25,-20,-15,-10,-8,-6,-4,-2,-1,-.5,-0.1,0.0,0.1,0.5,1,2,4,6,8,10,15,20,25,30,40]
endif; end QC
if strpos(fld_nm,'DSTSS') ne -1 then begin
	cntr_ntv=10
;	cntr_lvl_nbr=17
;	cntr_lvl=[0.0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100]
;	cntr_lvl_nbr=14
;	cntr_lvl=[0.0,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100]
	cntr_lvl_nbr=11
	cntr_lvl=[0.0,0.1,0.2,0.5,1,2,5,10,20,50,100]
endif; end DSTSS
if fld_nm eq 'DSTQ' or fld_nm eq 'DSTQ01' then begin
	cntr_lvl_nbr=11
	cntr_ntv=10
	cntr_lvl=[0.0,0.1,0.2,0.5,1,2,5,10,20,50,100]
endif; end DSTQ
if fld_nm eq 'CMFMC' then begin
	cntr_lvl_nbr=20
	cntr_ntv=2
	cntr_lvl=[0.0,0.01,0.1,0.5,1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
endif; 

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl
if unit_pfx_sng ne '' then begin
	print,'cntlabel_mk: WARNING unit_pfx_sng == ',unit_pfx_sng
	if top then begin
		unit_sng_psn=strpos(ttl,'(')
		if unit_sng_psn gt -1 then begin
			tmp=unit_sng_psn
			unit_sng_psn=strpos(ttl,'(',unit_sng_psn+1)
			if unit_sng_psn eq -1 then unit_sng_psn=tmp
		end; endif
		tail_sng=strmid(ttl,unit_sng_psn+1,strlen(ttl))
		head_sng=strmid(ttl,0,unit_sng_psn+1)
		ttl=head_sng+unit_pfx_sng+tail_sng
	endif; endif top
endif; endif unit

; Customize contour labels
cntr_chr_sz=0.83*chr_sz		;charsize of level labels
if pll eq 1 then begin
for lvl_idx=0,cntr_lvl_nbr-2 do begin
	if shd_xpr eq 'lt' then if cntr_lvl(lvl_idx) lt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'gt' then if cntr_lvl(lvl_idx) gt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'ge' then if cntr_lvl(lvl_idx) ge shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
endfor
; Blacken contours outside the min and max levels
;if cntr_lvl(cntr_lvl_nbr-2) lt data_max then cntr_fll_idx(cntr_lvl_nbr-2)=clr_blk_idx
;if cntr_lvl(0) gt data_min then cntr_fll_idx(0)=clr_blk_idx
for lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	if not dsh then cntr_ln_sty(lvl_idx)=0
	if cntr_lvl(lvl_idx) eq 0.0 then cntr_which_lbl(lvl_idx)=1
endfor
endif; endif pll eq 1

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cbar_idx=indgen(cbar_clr_nbr)+2
cbar_lgn=cntr_lvl_lbl

cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.86*chr_sz
cbar_unit=''
cbar_lbl_sz=chr_sz

if cbar then begin

cbar_psn=[.88,0.1,0.92,0.92]; [x_min,y_min,x_max,y_max]
plt_rgn_nrm=[0.12,0.1,0.87,0.92]; [x_min,y_min,x_max,y_max]

endif else begin

plt_rgn_nrm=[0.12,0.1,0.97,0.92]; [x_min,y_min,x_max,y_max] 

endelse

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,bgr
tvlct,r_curr,g_curr,b_curr

if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

if strpos(ttl,'ERBE') ne -1 and ord_max eq 60 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
ord=ord(1:52)
data=data(*,1:52)
endif; endif ERBE

mrg_top=2.5 ; 2 is default
mrg_btm=4 ; 4 is default
mrg_lft=8 ; 10 is default
if cbar then mrg_rgt=10 else mrg_rgt=2 ; 3 is default
if prn then begin
	if top then mrg_top=1.6 else mrg_top=0.5 ; 2 is default
	if btm then mrg_btm=1.7 else mrg_btm=0.5 ; 4 is default
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

; Set boundaries for Cylindrical projections
; NB: Until IDL 4.0, this commands was not needed because map_set() set x.crange and y.crange, which the axes routines need to do correct labeling
; Without the following two lines, x.crange and y.crange will be [0.0,1.0] and axis labels will be wrong
!x.range=[abc_min,abc_max] ; Longitude increases on x-axis
!y.range=[ord_max,ord_min] ; Pressure decreases on y-axis

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $ ; Latitude increases on x-axis
;	xrange=[abc_max,abc_min], $ ; Latitude decreases on x-axis
	yrange=[ord_max,ord_min], $ ; Pressure decreases on y-axis
	xstyle=13, $
	ystyle=1, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
;	position=plt_rgn_nrm, $
;	/closed, $			
	/fill, $			
	/noerase

if btm then ntt=1 else ntt=0
lat_axz_drw,lat_top=abc_max,lat_btm=abc_min,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,axz_vrt=0

if cnt_ovl then begin
contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot
endif; endif cnt_ovl

;if pnl_lbl ne '' then xyouts,0.1,0.9,pnl_lbl,size=chr_sz,/NORMAL
if pnl_lbl ne '' then xyouts,abc_min+0.05*(abc_max-abc_min),ord_min+0.05*(ord_max-ord_min),pnl_lbl,size=chr_sz,/DATA

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

end; end gcm_yz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM YZ commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM XZ
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
; CCM outputs data as    foo[time,lat,lev,lon]		in C
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_xz_bch, $
	bias=bias, $ ; plot biases of fields relative to observational analyses
	diff=diff, $ ; plot differences between model fields
	prs=prs, $ ; insert 'prs_' in filename (use files on pressure levels)
	fld_top=fld_top, $ ; first field is top plot
	epcr=epcr, $ ; Equatorial Pacific Cloud Response
	hck=hck, $ ; Hack
	bgr=bgr, $ ; background color?
	pll=pll, $ ; color palette
	cbar=cbar, $ ; Add color scale
	cnt_ovl=cnt_ovl, $ ; Add contour overlays
	prn=prn

if n_elements(epcr) eq 0 then epcr=0 ; Equatorial Pacific Cloud Response
if n_elements(hck) eq 0 then hck=0 ; Hack
if n_elements(prn) eq 0 then prn=0
if n_elements(bias) eq 0 then bias=0
if n_elements(diff) eq 0 then diff=0
if n_elements(prs) eq 0 then prs=0
if n_elements(bgr) eq 0 then bgr=0
if n_elements(pll) eq 0 then pll=3
if n_elements(cbar) eq 0 then cbar=0
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if n_elements(fld_top) eq 0 then fld_top=0
if prs then prs_sng='prs_' else prs_sng=''

fld2_stf=[ $
	['DSTQ','!5Dust Mixing Ratio !8Q!I!5d!N!5','!7l!5g kg!E-1!N'], $
	['DSTQ01','!5Clay Mixing Ratio !8Q!I!501!N!5','!7l!5g kg!E-1!N']]
;	['U','!5Zonal Wind','!5m s!E-1!N'], $
;	['T','!5Temperature','!E!12_!5!NK'] ]
;	['Q','!5Vapor Mixing Ratio','!5g kg!E-1!N'] ]
;	['CMFMC','!5Convective Mass Flux','!8M!5!IC!N g m!E-2!N s!E-1!N']]
;	['QICE','!5Ice Mixing Ratio','!5mg kg!E-1!N'], $
;	['CLOUD','!5Cloud Amount','!5%'], $
;	['QC','!5Condensate Mixing Ratio','!5mg kg!E-1!N'] ]
;	['QRS','!5Solar Heating','!5K day!E-1!N'], $
;	['QRL','!5Longwave Heating','!5K day!E-1!N'], $
;	['RADD','!5Net Radiative Heating','!5K day!E-1!N'], $
;	['HGS','!5Resolved','!5K day!E-1!N'], $
;	['CMFDT','!5Convective','!5K day!E-1!N'], $
;	['QDIABAT','!5Total Diabatic Heating','!5K day!E-1!N'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
;	['ecmwf','!5ECMWF'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'], $
;	['sld012d','!5CCM3'], $
;	['spcp_85','!5ANV']]
	['dst16','!5CCM']]
;	['dst17','!5TeF94']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

;mth_sng=['01','07']
;mth_sng=['07']
mth_sng=['0305']
mth_nbr=n_elements(mth_sng)
for mth_idx=0,mth_nbr-1 do begin
for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
;	if prs then y_ttl='!5Pressure (mb)' else y_ttl='!7g!5 !9X!5 1000 (mb)'
	y_ttl='!5Pressure (mb)'
	x_ttl=''
	pre_scl=0
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	dsh=0
	yr_sng='8589'
;	rng_x=[0,360]
;	rng_x=[120,180]
	rng_x=[300,360]
;	tck_lcn=[30,3]
;	tck_lcn=[60,6]
;	tck_lcn=[10,5]
	rng_x_sng=''
	rng_x_sng=string(format='(I3.3)',rng_x(0))+'E'+string(format='(I3.3)',rng_x(1))+'E'
	;	yavg_sng='yavg_10S10N'
;	yavg_sng='yavg_00N10N'
;	yavg_sng='yavg_00N20N'
;	yavg_sng='yavg_00N30N'
	yavg_sng='y00N15N'
;	yavg_sng='yavg_20S00N'
;	yavg_sng='yavg_30S00N'
;	yavg_sng='yavg_30S30N'
	shd_val=0.0
	shd_xpr='lt'
	if fld2(fld_idx2).nm eq 'U' then if diff then pre_scl=[-2,0.2,21] else pre_scl=[-20,5,16]
	if fld2(fld_idx2).nm eq 'T' then if diff then pre_scl=[-2,0.2,21] else pre_scl=[200,5,20]
	if fld2(fld_idx2).nm eq 'Q' then scl=1000.0
	if fld2(fld_idx2).nm eq 'Q' then if diff then pre_scl=[-3,.3,21] else pre_scl=[0,1,30]
;	if fld2(fld_idx2).nm eq 'Q' then scl=1.0e6
;	if fld2(fld_idx2).nm eq 'Q' then if diff then pre_scl=[-10,1,21] else pre_scl=[0,1,30]
	if fld2(fld_idx2).nm eq 'QICE' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QICE' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QC' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QC' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CMFMC' then scl=1000
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then pre_scl=[0,2,30] else pre_scl=[-30,2,30]
;	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then pre_scl=[0.01,.01,20] else pre_scl=[0.01,.01,20]
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QDIABAT' then scl=86400
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'CMFDT' then scl=86400
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'HGS' then scl=86400
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'QRS' then scl=86400
	if fld2(fld_idx2).nm eq 'QRS' then if not diff then pre_scl=[-2.0,0.2,30] else pre_scl=[-2.0,0.2,30]
	if fld2(fld_idx2).nm eq 'QRL' then scl=86400
	if fld2(fld_idx2).nm eq 'QRL' then if not diff then pre_scl=[-4.,0.2,30] else pre_scl=[-2.0,0.2,30]
	if fld2(fld_idx2).nm eq 'RADD' then if not diff then pre_scl=[-4.,0.2,30] else pre_scl=[-2.0,0.2,30] 
;	if fld2(fld_idx2).nm eq 'RADD' then if not diff then shd_val=6
;	if fld2(fld_idx2).nm eq 'RADD' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CLOUD' then scl=100
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then shd_val=25
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then pre_scl=[-50,5,30] else pre_scl=[-50,5,30]  
	if strpos(fld2(fld_idx2).nm,'DSTQ') ne -1 then begin
		scl=1.0e9
	endif ; endif dust
	if strpos(xpt(xpt_idx).nm,'ecmwf') ne -1 then begin
		if mth_sng(mth_idx) eq '01' then yr_sng='9095_' else yr_sng='8994_'
	endif
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+yr_sng+'_'+mth_sng(mth_idx)+'_'+yavg_sng+'.nc'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+prs_sng+yr_sng+'_'+mth_sng(mth_idx)+'_'+yavg_sng+'_'+rng_x_sng+'_'+fld2(fld_idx2).nm+'.eps'
	if diff then begin ; plot differences between simulated fields
		ttl='ANV-CCM '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+yavg_sng+mth_sng(mth_idx)+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+yavg_sng+mth_sng(mth_idx)+'_'+fld2(fld_idx2).nm+'.eps'
	endif; endif diff
	if prn and not bias and xpt_nbr eq 1 then begin
;		ttl=fld2(fld_idx2).sng
		if fld_top then begin
			if fld_idx2 eq fld_nbr2-1 then btm=1 else btm=0
			if mth_sng(mth_idx) eq '01' then ttl='Jan '+ttl else if mth_sng(mth_idx) eq '07' then ttl='Jul '+ttl
		endif; endif fld_top
		if not fld_top then begin
			if mth_idx eq 0 then top=1 else top=0
			if mth_idx eq mth_nbr-1 then btm=1 else btm=0
		endif; endif not fld_top
	endif; endif prn not bias
	if bias then begin ; plot biases of simulated fields relative to observational analyses
	if strpos(xpt(xpt_idx).nm,'amip5') ne -1 or strpos(xpt(xpt_idx).nm,'spcp') ne -1 then begin
		ttl=xpt(xpt_idx).sng+'-ECMWF '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		if mth_sng(mth_idx) eq '01' then ecmwf_yr_sng='9095_'
		if mth_sng(mth_idx) eq '07' then ecmwf_yr_sng='8994_'
		if fld2(fld_idx2).nm eq 'U' then pre_scl=[-20,2,20]
		if fld2(fld_idx2).nm eq 'T' then pre_scl=[-15,1,30]
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+yr_sng+'ecmwf_'+ecmwf_yr_sng+'yavg_'+mth_sng(mth_idx)+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+yr_sng+'ecmwf_'+ecmwf_yr_sng+'yavg_'+mth_sng(mth_idx)+'_'+fld2(fld_idx2).nm+'.eps'
	endif; endif xpt is a simulation, not ecmwf
	endif; endif bias
	if (prn and bias) or (prn and not bias and xpt_nbr gt 1) then begin
		if mth_sng(mth_idx) eq '01' then ttl='Jan '+fld2(fld_idx2).sng
		if mth_sng(mth_idx) eq '07' then ttl='Jul '+fld2(fld_idx2).sng
		if xpt_idx eq 0 then top=1 else top=0
		if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
	endif; endif prn and bias
;	Hack to look at specific files here
	if epcr then begin
		ttl='Equatorial Pacific Cloud Response'
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_yavg_10S10N_87m85_MAM.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_yavg_10S10N_87m85_MAM_QC.eps'
		pre_scl=[-20,2,20]
		rng_x=[120,270]
		tck_lcn=[30,3]
	endif; end hack
	if prn then begin
		x_sz=6.5
;		y_sz=3.2
		y_sz=5.0
		if top then y_sz=y_sz*1.08
		if btm then y_sz=y_sz*1.08
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_xz,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,bgr=bgr,pll=pll,order=1,cbar=cbar,cnt_ovl=cnt_ovl,fld_nm=fld2(fld_idx2).nm,wrp=0,rng_x=rng_x,rng_y=[0,1000],pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,top=top,btm=btm,scl=scl,shd_val=shd_val,shd_xpr=shd_xpr,tck_lcn=tck_lcn
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over experiment
endfor; end loop over fld2
endfor; end loop over mth

end; end gcm_xz_bch()
pro gcm_xz, $
	fl_nm=fl_nm, $
	order=order, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	shd_val=shd_val, $
	shd_xpr=shd_xpr, $
	tck_lcn=tck_lcn, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	dsh=dsh, $
	wrp=wrp, $
	bgr=bgr, $
	pll=pll, $
	cnt_ovl=cnt_ovl, $
	cbar=cbar, $
	scl=scl, $
	pre_scl=pre_scl, $
	fld_nm=fld_nm

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(wrp) eq 0 then wrp=0
if n_elements(shd_val) eq 0 then shd_val=0.0
if n_elements(shd_xpr) eq 0 then shd_xpr='lt'
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ANV Zonal Wind !8U!5 (m s!E-1!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl='!5Pressure (mb)'
if n_elements(x_ttl) eq 0 then x_ttl='!5Longitude'
if n_elements(order) ne 0 then clr_ord=order else clr_ord=1
if n_elements(bgr) eq 0 then bgr=0
if n_elements(pll) eq 0 then pll=3
if n_elements(cbar) eq 0 then cbar=1
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/spcp_85/spcp_85_yavg_8589_01.nc'
if n_elements(fld_nm) eq 0 then fld_nm='U'

if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=34
clr_rfr,clr_tbl_nbr,"",0

erase
!p.multi=0

; Read in binary data from the NetCDF file 
print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr
lev_nm='lev'
if fld_nm eq 'CMFMC' or fld_nm eq 'MPSI' then if strpos(fl_nm,'prs') eq -1 then lev_nm='ilev'
lev_id=ncdf_dimid(nc_id,lev_nm)
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
ncdf_varget,nc_id,lev_nm,lev
ncdf_varget,nc_id,'lon',lon
ncdf_varget,nc_id,fld_nm,data
ncdf_close,nc_id
; End of NetCDF commands

data=reform(data)
good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=lon
ord=lev

data_min=min(data)
data_max=max(data)

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

if auto_scl then begin

; Automatically set and scale the contour levels
cntr_lvl_nbr=20

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

if fld_nm eq 'DSTQ' or fld_nm eq 'DSTQ01' then begin
	cntr_ntv=10
;	cntr_lvl=[0.0,0.1,0.2,.5,1,2,5,10,20,50,100,200,500]
	cntr_lvl=[0.0,1,2,5,10,15,20,35,50,75,100,200,500]
	cntr_lvl_nbr=n_elements(cntr_lvl)
endif; end DSTQ

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl
if unit_pfx_sng ne '' then print,'cntr_lbl_mk: WARNING unit_pfx_sng == ',unit_pfx_sng

; Customize contour labels
cntr_chr_sz=0.83*chr_sz		;charsize of level labels
if pll eq 1 then begin
for lvl_idx=0,cntr_lvl_nbr-2 do begin
	if shd_xpr eq 'lt' then if cntr_lvl(lvl_idx) lt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'gt' then if cntr_lvl(lvl_idx) gt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'ge' then if cntr_lvl(lvl_idx) ge shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
endfor
; Blacken contours outside the min and max levels
;if cntr_lvl(cntr_lvl_nbr-2) lt data_max then cntr_fll_idx(cntr_lvl_nbr-2)=clr_blk_idx
;if cntr_lvl(0) gt data_min then cntr_fll_idx(0)=clr_blk_idx
for lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	if not dsh then cntr_ln_sty(lvl_idx)=0
	if cntr_lvl(lvl_idx) eq 0.0 then cntr_which_lbl(lvl_idx)=1
endfor
endif; endif pll eq 1

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cbar_idx=indgen(cbar_clr_nbr)+2
cbar_lgn=cntr_lvl_lbl

cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.86*chr_sz
cbar_unit=''
cbar_lbl_sz=chr_sz

if cbar then begin

cbar_psn=[.88,0.1,0.92,0.92]; [x_min,y_min,x_max,y_max]
plt_rgn_nrm=[0.12,0.1,0.87,0.92]; [x_min,y_min,x_max,y_max]

endif else begin

plt_rgn_nrm=[0.12,0.1,0.97,0.92]; [x_min,y_min,x_max,y_max] 

endelse

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,bgr
tvlct,r_curr,g_curr,b_curr

if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

if strpos(fl_nm,'erbe') ne -1 and ord_max eq 60 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
ord=ord(1:52)
data=data(*,1:52)
endif; endif ERBE

mrg_top=2.5 ; 2 is default
mrg_btm=4 ; 4 is default
mrg_lft=8 ; 10 is default
if cbar then mrg_rgt=10 else mrg_rgt=2 ; 3 is default
if prn then begin
	if top then mrg_top=1.6 else mrg_top=0.5 ; 2 is default
	if btm then mrg_btm=1.7 else mrg_btm=0.5 ; 4 is default
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $ ; Longitude increases on x-axis
	yrange=[ord_max,ord_min], $ ; Pressure decreases on y-axis
	xstyle=13, $
	ystyle=1, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
;	position=plt_rgn_nrm, $
;	/closed, $			
	/fill, $			
	/noerase

if btm then ntt=1 else ntt=0
lon_axz_drw,lon_top=abc_max,lon_btm=abc_min,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,axz_vrt=0,tck_lcn=tck_lcn

if cnt_ovl then begin
contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot
endif; endif cnt_ovl

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

end; end gcm_xz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM XZ commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM TZV
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
; CCM outputs data as    foo[time,lat,lev,lon]		in C
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_tzv_bch, $
	bias=bias, $ ; plot biases of fields relative to observational analyses
	diff=diff, $ ; plot differences between model fields
	prs=prs, $ ; insert 'prs_' in filename (use files on pressure levels)
	fld_top=fld_top, $ ; first field is top plot
	prn=prn
if n_elements(prn) eq 0 then prn=0
if n_elements(bias) eq 0 then bias=0
if n_elements(diff) eq 0 then diff=0
if n_elements(prs) eq 0 then prs=0
if n_elements(fld_top) eq 0 then fld_top=0
if prs then prs_sng='prs_' else prs_sng=''

fld2_stf=[ $
;	['U','!5Zonal Wind','!5m s!E-1!N'], $
;	['T','!5Temperature','!E!12_!5!NK'] ]
;	['Q','!5Vapor Mixing Ratio','!5g kg!E-1!N'] ]
;	['CMFMC','!5Convective Mass Flux','!8M!5!IC!N g m!E-2!N s!E-1!N']]
;	['QICE','!5Ice Mixing Ratio','!5mg kg!E-1!N'], $
;	['CLOUD','!5Cloud Amount','!5%'], $
	['QC','!5Condensate Mixing Ratio','!5mg kg!E-1!N'] ]
;	['QRS','!5Solar Heating','!5K day!E-1!N'], $
;	['QRL','!5Longwave Heating','!5K day!E-1!N'], $
;	['RADD','!5Net Radiative Heating','!5K day!E-1!N'], $
;	['HGS','!5Resolved','!5K day!E-1!N'], $
;	['CMFDT','!5Convective','!5K day!E-1!N'], $
;	['QDIABAT','!5Total Diabatic Heating','!5K day!E-1!N'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
;	['ecmwf','!5ECMWF'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
	['amip5','!5CCM'], $
;	['sld012d','!5CCM3'], $
	['spcp_85','!5ANV']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

rgn_stf=[ $
	['Indian_Central','Central Indian Ocean'] ]
;	['Pacific_ITCZ_West','West Pacific ITCZ'] ]
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct_gcm,idx:0,nm:'',sng:''}
rgn=replicate({rgn_sct_gcm},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	rgn(idx).idx=idx
	rgn(idx).nm=rgn_stf(0,idx)
	rgn(idx).sng=rgn_stf(1,idx)
endfor; end loop over rgn_lst

for rgn_idx=0,rgn_nbr-1 do begin
for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	ttl=xpt(xpt_idx).sng+' '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
;	if prs then y_ttl='!5Pressure (mb)' else y_ttl='!7g!5 !9X!5 1000 (mb)'
	y_ttl='!5Pressure (mb)'
	x_ttl=''
	pre_scl=0
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	dsh=0
	mth_sng='0112'
	yr_sng='8589_'
	shd_val=0.0
	shd_xpr='lt'
	if fld2(fld_idx2).nm eq 'U' then if diff then pre_scl=[-2,0.2,21] else pre_scl=[-20,5,16]
	if fld2(fld_idx2).nm eq 'T' then if diff then pre_scl=[-2,0.2,21] else pre_scl=[200,5,20]
	if fld2(fld_idx2).nm eq 'Q' then scl=1000.0
	if fld2(fld_idx2).nm eq 'Q' then if diff then pre_scl=[-3,.3,21] else pre_scl=[0,1,30]
;	if fld2(fld_idx2).nm eq 'Q' then scl=1.0e6
;	if fld2(fld_idx2).nm eq 'Q' then if diff then pre_scl=[-10,1,21] else pre_scl=[0,1,30]
	if fld2(fld_idx2).nm eq 'QICE' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QICE' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QICE' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QC' then scl=1.0e6
	if fld2(fld_idx2).nm eq 'QC' then if diff then pre_scl=[-15,1,30] else pre_scl=[0,2,30]
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'QC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CMFMC' then scl=1000
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then pre_scl=[0,2,30] else pre_scl=[-30,2,30]
;	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then pre_scl=[0.01,.01,20] else pre_scl=[0.01,.01,20]
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_val=6
	if fld2(fld_idx2).nm eq 'CMFMC' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QDIABAT' then scl=86400
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'QDIABAT' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'CMFDT' then scl=86400
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CMFDT' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'HGS' then scl=86400
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_val=4.
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'HGS' then if not diff then pre_scl=[0.0,2,20] else pre_scl=[-10,1,20]
	if fld2(fld_idx2).nm eq 'QRS' then scl=86400
	if fld2(fld_idx2).nm eq 'QRS' then if not diff then pre_scl=[-2.0,0.2,30] else pre_scl=[-2.0,0.2,30]
	if fld2(fld_idx2).nm eq 'QRL' then scl=86400
	if fld2(fld_idx2).nm eq 'QRL' then if not diff then pre_scl=[-4.,0.2,30] else pre_scl=[-2.0,0.2,30]
	if fld2(fld_idx2).nm eq 'RADD' then if not diff then pre_scl=[-4.,0.2,30] else pre_scl=[-2.0,0.2,30] 
;	if fld2(fld_idx2).nm eq 'RADD' then if not diff then shd_val=6
;	if fld2(fld_idx2).nm eq 'RADD' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CLOUD' then scl=100
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then shd_val=25
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then shd_xpr='ge'
	if fld2(fld_idx2).nm eq 'CLOUD' then if not diff then pre_scl=[-50,5,30] else pre_scl=[-50,5,30]  
	if strpos(xpt(xpt_idx).nm,'ecmwf') ne -1 then begin
	endif
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+yr_sng+mth_sng+'.nc'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+prs_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+yr_sng+mth_sng+'_'+fld2(fld_idx2).nm+'.eps'
	if diff then begin ; plot differences between simulated fields
		ttl='ANV-CCM '+fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+mth_sng+'.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+mth_sng+'_'+fld2(fld_idx2).nm+'.eps'
	endif; endif diff
	if prn and not bias then begin
		ttl=fld2(fld_idx2).sng
		if fld_top then begin
			if fld_idx2 eq fld_nbr2-1 then btm=1 else btm=0
		endif; endif fld_top
		if not fld_top then begin
			if xpt_idx eq 0 then top=1 else top=0
			if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
		endif; endif not fld_top
	endif; endif prn not bias
	if prn then begin
		x_sz=6.5
		y_sz=3.2
		if top then y_sz=y_sz*1.08
		if btm then y_sz=y_sz*1.08
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	gcm_tzv,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,pll=1,order=1,fld_nm=fld2(fld_idx2).nm,cbar=0,wrp=1,rng_x=rng_x,rng_y=[0,1000],pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,top=top,btm=btm,scl=scl,shd_val=shd_val,shd_xpr=shd_xpr,tck_lcn=tck_lcn
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over xpt
endfor; end loop over fld2
endfor; end loop over rgn

end; end gcm_tzv_bch()
pro gcm_tzv, $
	fl_nm=fl_nm, $
	pll=pll, $
	order=order, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	shd_val=shd_val, $
	shd_xpr=shd_xpr, $
	tck_lcn=tck_lcn, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	dsh=dsh, $
	wrp=wrp, $
	cbar=cbar, $
	scl=scl, $
	pre_scl=pre_scl, $
	fld_nm=fld_nm

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(wrp) eq 0 then wrp=1
if n_elements(shd_val) eq 0 then shd_val=0.0
if n_elements(shd_xpr) eq 0 then shd_xpr='lt'
if n_elements(cbar) eq 0 then cbar=1
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ANV Condensate !8Q!5!IC!N (mg kg!E-1!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl='!5Pressure (mb)'
if n_elements(x_ttl) eq 0 then x_ttl='!5Month'
if n_elements(order) ne 0 then clr_ord=order else clr_ord=1
if n_elements(pll) eq 0 then pll=3
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xyavg_rgn_Indian_Central_8589_0112.nc'
if n_elements(fld_nm) eq 0 then fld_nm='QC'

if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=34
clr_rfr,clr_tbl_nbr,"",0

erase
!p.multi=0

; Read in binary data from the NetCDF file 
print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
lev_nm='lev'
if fld_nm eq 'CMFMC' or fld_nm eq 'MPSI' then if strpos(fl_nm,'prs') eq -1 then lev_nm='ilev'
lev_id=ncdf_dimid(nc_id,lev_nm)
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
ncdf_varget,nc_id,lev_nm,lev
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,fld_nm,data
ncdf_close,nc_id
; End of NetCDF commands

data=transpose(reform(data))
good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
ord=lev
abc=indgen(time_nbr)
if wrp then begin
data=[data,data(0,*)]
abc=[abc,time_nbr]
endif

data_min=min(data)
data_max=max(data)

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

if auto_scl then begin

; Automatically set and scale the contour levels
cntr_lvl_nbr=20

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl
if unit_pfx_sng ne '' then print,'cntr_lbl_mk: WARNING unit_pfx_sng == ',unit_pfx_sng

; Customize contour labels
cntr_chr_sz=0.83*chr_sz		;charsize of level labels
if pll eq 1 then begin
for lvl_idx=0,cntr_lvl_nbr-2 do begin
	if shd_xpr eq 'lt' then if cntr_lvl(lvl_idx) lt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'gt' then if cntr_lvl(lvl_idx) gt shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
	if shd_xpr eq 'ge' then if cntr_lvl(lvl_idx) ge shd_val then cntr_fll_idx(lvl_idx)=10 else cntr_fll_idx(lvl_idx)=clr_wht_idx
endfor
; Blacken contours outside the min and max levels
;if cntr_lvl(cntr_lvl_nbr-2) lt data_max then cntr_fll_idx(cntr_lvl_nbr-2)=clr_blk_idx
;if cntr_lvl(0) gt data_min then cntr_fll_idx(0)=clr_blk_idx
for lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	if not dsh then cntr_ln_sty(lvl_idx)=0
	if cntr_lvl(lvl_idx) eq 0.0 then cntr_which_lbl(lvl_idx)=1
endfor
endif; endif pll eq 1

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cbar_idx=indgen(cbar_clr_nbr)+2
cbar_lgn=cntr_lvl_lbl

cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.86*chr_sz
cbar_unit=''
cbar_lbl_sz=chr_sz

if cbar then begin

cbar_psn=[.88,0.1,0.92,0.92]; [x_min,y_min,x_max,y_max]
plt_rgn_nrm=[0.12,0.1,0.87,0.92]; [x_min,y_min,x_max,y_max]

endif else begin

plt_rgn_nrm=[0.12,0.1,0.97,0.92]; [x_min,y_min,x_max,y_max] 

endelse

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,0
tvlct,r_curr,g_curr,b_curr

if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

if strpos(fl_nm,'erbe') ne -1 and ord_max eq 60 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
ord=ord(1:52)
data=data(*,1:52)
endif

mrg_top=2.5 ; 2 is default
mrg_btm=4 ; 4 is default
mrg_lft=8 ; 10 is default
if cbar then mrg_rgt=10 else mrg_rgt=2 ; 3 is default
if prn then begin
	if top then mrg_top=1.6 else mrg_top=0.5 ; 2 is default
	if btm then mrg_btm=1.7 else mrg_btm=0.5 ; 4 is default
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $ ; time increases on x-axis
	yrange=[ord_max,ord_min], $ ; Pressure decreases on y-axis
	xstyle=13, $
	ystyle=1, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
;	position=plt_rgn_nrm, $
;	/closed, $			
	/fill, $			
	/noerase

if btm then ntt=1 else ntt=0
time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0

contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

end; end gcm_tzv()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM TZV commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM XY
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_xy_bch, $
	fld_top=fld_top, $ ; First field is top plot
	bias=bias, $ ; Plot biases of fields relative to observational analyses
	diff=diff, $ ; Plot differences between model fields
	lvl=lvl, $
	bg=bg, $
	jtk=jtk, $
	ams=ams, $
	ocn=ocn, $
	order=order, $ ; Order of palette
	pll=pll, $ ; Palette number
	pnl_lbl=pnl_lbl, $ ; User-specified Panel label
	lbl_typ=lbl_typ, $ ; Default is 'none', alternate is 'auto'
	img=img, $
	cnt_ovl=cnt_ovl, $
	csn=csn, $ ; Rearrange headers and footers for seasonal four panel plots
	azi=azi, $
	rgn=rgn, $ ; Insert 'rgn_' in filename
	prs=prs, $ ; Insert 'prs_' in filename (use files on pressure levels)
	dly=dly, $ ; Daily rather than monthly data
	yrl=yrl, $ ; Yearly rather than monthly data
	dvc=dvc, $ ; Output device (ps, eps, jpg, mpg, ppm ...)
	sgl=sgl, $ ; Input is contained in single file (not implemented yet)
	fl_usr=fl_usr, $ ; User-specified input file
	pwz=pwz, $ ; Pause between frames
	prn=prn

; Usage
; gcm_xy_bch,img=1,dly=1,dvc='mpg',fl_usr=getenv('DATA')+'/dstmch26/dstmch26_19980101_19981231.nc',pwz=0 <-- Make MPEG movie
; gcm_xy_bch,pll=4,img=0,prn=0,fld_top=1,fl_usr=getenv('DATA')+'/data/dst_1x1_ZBN03.nc',rgn='Bdl'
; gcm_xy_bch,pll=4,img=0,prn=0,fld_top=1,fl_usr=getenv('DATA')+'/data/dst_T62_ZBN03.nc'
; gcm_xy_bch,ams=1,diff=1,img=0,prn=0
; gcm_xy_bch,jtk=1,img=1,prn=0
; gcm_xy_bch,img=1,prn=0,rgn='TNA'
; gcm_xy_bch,bg=1,prn=0
; gcm_xy_bch,pll=4,img=0,prn=0,fld_top=1 <-- Annual deposition data

if n_elements(fld_top) eq 0 then fld_top=0
if n_elements(bias) eq 0 then bias=0
if n_elements(diff) eq 0 then diff=0
if n_elements(lvl_idx) eq 0 then lvl_idx=-1
if n_elements(bg) eq 0 then bg=0
if n_elements(jtk) eq 0 then jtk=0
if n_elements(ams) eq 0 then ams=0
if n_elements(order) eq 0 then order=1 ; Order of palette
if n_elements(pll) eq 0 then pll=1 ; Palette number
if n_elements(img) eq 0 then img=1
if n_elements(cnt_ovl) eq 0 then cnt_ovl=0
if n_elements(csn) eq 0 then csn=0
if n_elements(azi) eq 0 then azi=0
if n_elements(prn) eq 0 then prn=0
if n_elements(dvc) eq 0 then dvc='eps' ; Output device (ps, eps, jpg, mpg, ppm ...)
if n_elements(fl_usr) eq 0 then fl_usr='' ; User-specified input file
if n_elements(sgl) eq 0 then sgl=1 ; Input is contained in single file
if n_elements(mth) eq 0 then mth=1 ; Monthly data
if n_elements(dly) eq 0 then dly=0 ; Daily data
if n_elements(yrl) eq 0 then yrl=0 ; Yearly data
if n_elements(prs) eq 0 then prs=0
if n_elements(pwz) eq 0 then pwz=0 ; Pause between frames
if prs then prs_sng='prs_' else prs_sng=''
if n_elements(rgn) eq 0 then rgn=0
if rgn then rgn_sng=rgn+'_' else rgn_sng=''
if n_elements(ocn) eq 0 then ocn=0
if ocn then ocn_sng='_ocn' else ocn_sng=''
if n_elements(pnl_lbl) eq 0 then pnl_lbl=''
if n_elements(lbl_typ) eq 0 then lbl_typ='none' ; 'none' or 'auto' are valid
if yrl then mth=0 ; Yearly data
if dly then mth=0 ; Daily data
if mth then dly=0 ; Monthly data

fld2_stf=[ $
;	['T','!5Temperature','!E!12_!5!NK'] ]
;	['Q','!5Water Vapor','!5kg kg!E-1!N'], $
;	['U','!5Zonal Wind','!5m s!E-1!N'] ]
;	['CHI','!5Velocity Potential','!9X!5 10!E6!N !5m!E2!N s!E-1!N'] ]
;	['PSI','!5Streamfunction','!9X!5 10!E6!N !5m!E2!N s!E-1!N'] ]
;	['DIV','!5Divergence','!9X!5 10!E-6!N !5s!E-1!N'] ]
;	['VOR','!5Vorticity','!9X!5 10!E-6!N !5s!E-1!N'] ]
;	['TOTLWP','!5Liquid Water Path','!5g m!E-2!N']]
;	['TOTIWP','!5Ice Water Path','!5g m!E-2!N'], $
;	['TOTCWP','!5Condensed Water Path','!5g m!E-2!N'] ]
;	['PBLH','!5Boundary Layer Height','!5m'] ]
;	['GCLD','!5Total Greenhouse Effect','!5W m!E-2!N'] ]
;	['GCLR','!5Atmospheric Greenhouse Effect','!5W m!E-2!N'] ]
;	['TMQ','!5Column Vapor','!5kg m!E-2!N'], $
;	['VMAGSFC','!5Wind Speed','!5m s!E-1!N'] ]
;	['LHFLX','!5Latent Heat Flux','!5W m!E-2!N'], $
;	['SHFLX','!5Sensible Heat Flux','!5W m!E-2!N'], $
;	['NET','!5Net Surface Energy','!5W m!E-2!N'], $
;	['TS1','!5Surface Temperature','!E!12_!5!NK'], $
;	['ALBEDO','!5Albedo','!5%'] ]
;	['SOLIN','!5Insolation','!5W m!E-2!N'] ]
;	['FSUT','!5Reflected Solar Flux','!5W m!E-2!N'] ]
;	['FLUS','!5Surface-Emitted Longwave Flux','!5W m!E-2!N'] ]
;	['PSL','!5Sea Level Pressure','!5mb'], $
;	['PS','!5Surface Pressure','!5mb'] ]
;	['PRECC','!5Convective Precipitation','!5mm day!E-1!N'], $
;	['PRECL','!5Large Scale Precipitation','!5mm day!E-1!N'], $
;	['PRECT','!5Total Precipitation','!5mm day!E-1!N'], $
;	['RATIO','!5Cloudy Sky Absorption Enhancement','!5Ratio'] ]
;	['NCF','!5Net Cloud Forcing','!5W m!E-2!N'] ]
;	['SWCF','!5Shortwave Cloud Forcing','!5W m!E-2!N'] ]
;	['LWCF','!5Longwave Cloud Forcing','!5W m!E-2!N'] ]
;	['ORO','',''] ]
;	['Z2TEST','!5Geopotential Height','!5gpm'] ]
;	['fracla','!5Land Fraction','!5%'] ]
;	['aotm14','!5AOD 0.63 !7l!5m',''] ]
;	['odxc_aer_630','!5AOD 0.63 !7l!5m',''] ]
;	['aots14','!5AOD 0.63 !7l!5m',''] ]
;	['aer_idx_331_360','!5AAI 331-360 nm',''] ]
;	['mbl_bsn_fct','!5Geomorphic Erodibility','']]
;	['rdb_fct_hyd','!5Hydrologic Erodibility',''], $
;	['rdb_fct_tpg','!5Topographic Erodibility',''], $
;        ['rdb_fct_gmr','!5Geomorphic Erodibility',''], $
;	['rdb_fct_rfl_mds_lnr','!5Reflectance-Based Erodibility (MODIS)','']]
;	['rdb_fct_rfl_mds_sqr','!5Reflectance-Based Erodibility (MODIS!E2!N)','']]
;	['src_flg','!5Source Flag',''], $
;	['src_frq','!5Source Frequency','!5%'], $
;	['fsh_fct','!5Efficiency Factor',''], $
;	['csn_cyc','!5Seasonal Cycle','!5%'] ]
;	['lai_lsm','!5Leaf Area Index','!5m!E2!N m!E-2'] ]
;	['mss_frc_CaCO3','!5Mass Fraction of CaCO!I3!N','!5%'] ]
;	['VAI_DST','!5Vegetation Area Index','!5m!E2!N m!E-2!N'] ]
;	['DSTMPC01','!5Clay Loading','!5mg m!E-2!N'], $
;	['DSTMPC02','!5Fine Silt Loading','!5mg m!E-2!N'] ]
;	['DSTOD380','!5AOD 0.38 !7l!5m',''] ]
;	['DSTODX01','!5Clay AOD',''], $
;	['DSTODX02','!5Fine Silt AOD',''] ]
;	['DSTMPC','!5Dust Loading','!5mg m!E-2!N'] ]
;	['DSTMPCCA','!5Ca Loading','!5mg m!E-2!N'] ]
;	['DSTODXC','!5AOD 0.63 !7l!5m',''] ]
;       ['AERAFRC_TOP','!5Soot+Dust: Column Heating','!5W m!E-2!N'] ]
;	['AERAFRC_ATM','!5Soot+Dust: Atm. Absorption','!5W m!E-2!N'] ]
;	['AERAFRC_SFC','!5Soot+Dust: Surface cooling','!5W m!E-2!N'] ]
;	['SNOAERFRC','!5Soot+Dust: Snowpack Heating','!5W m!E-2!N'] ]
	['AEROD_v','!5Aerosol Optical Depth at 0.5 !7l!5m',''] ]
;	['DPSWETFRC','!5(Wet Deposition)/(Total Deposition)','!5%'] ]
;	['DSTSFMBL','!5Dust Mobilization','!7l!5g m!E-2!N s!E-1!N'], $
;	['DSTSFDPS','!5Dust Deposition','!7l!5g m!E-2!N s!E-1!N']]
;	['DSTSFDPSCA','!5Ca Deposition','!5ng m!E-2!N s!E-1!N']]
;	['DSTSFDRY','!5Dust Dry Deposition','!7l!5g m!E-2!N s!E-1!N'], $
;	['DSTSFPCP','!5Dust Wet Deposition','!7l!5g m!E-2!N s!E-1!N'] ]
;	['DSTSFGRV','!5Dust Gravitational Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFTRB','!5Dust Turbulent Deposition','!5ng m!E-2!N s!E-1!N']]
;	['DSTSFDPS','!5Dust Deposition','!5mg cm!E-2!N ka!E-1!N'] ] ; CGS UNITS!
;	['DSTSFDPSDST','!5Dust Deposition','!5mg cm!E-2!N ka!E-1!N'], $ ; CGS UNITS!
;	['DSTSFDPSCACO3','!5CaCO!I3!N Deposition','!5mg cm!E-2!N ka!E-1!N'] ] ; CGS UNITS!
;	['DSTSFDPSC','!5Particulate C Deposition','!7l!5g cm!E-2!N ka!E-1!N'], $ ; CGS UNITS!
;	['DSTSFDPSRATCACO3','!5CaCO!I3!N / Dust Deposition Ratio','!9X!5 10!E-3!N'] ]
;	['DSTSFD01','!5Clay Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFG01','!5Clay Gravitational Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFM01','!5Clay Mobilization','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFP01','!5Clay Wet Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFT01','!5Clay Turbulent Deposition','!5ng m!E-2!N s!E-1!N'] ]
;	['DSTSFD02','!5Fine Silt Deposition','!5ng m!E-2!N s!E-1!N'] ]
;	['DSTSFG02','!5Fine Silt Gravitational Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFM02','!5Fine Silt Mobilization','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFP02','!5Fine Silt Wet Deposition','!5ng m!E-2!N s!E-1!N'], $
;	['DSTSFT02','!5Fine Silt Turbulent Deposition','!5ng m!E-2!N s!E-1!N'] ]
;	['USTAR','!5Friction Velocity','!5m s!E-1!N'] ]
;	['WND_MBL','!510 m Wind Speed','!5m s!E-1!N'] ]
;	['U10','!5Zonal 10 m Wind Speed','!5m s!E-1!N'], $
;	['V10','!5Meridional 10 m Wind Speed','!5m s!E-1!N'] ]
;	['SMP_SFC','!5Surface Soil Matric Potential','!5mm'] ]
;	['VWC_SFC','!5Surface Volumetric Water Content','!5m!E3!N m!E-3!N'] ]
;	['FRC_WET','!5Moisture Inhibition','!5fraction'] ]
;	['vwc_sfc','!5Surface Volumetric Water Content','!5m!E3!N m!E-3!N'] ]
;	['dust_ttl','!5Mineral Dust','!7l!5g kg!E-1!N'] ]
;	['FTNS','!5Surface Absorbed Radiative Flux','!5W m!E-2!N'] ]
;	['FTNSDST','!5Tracer Net Surface Radiative Flux','!5W m!E-2!N'] ]
;	['FTNSFRC','!5Tracer Net Surface Radiative Forcing','!5W m!E-2!N'] ]
;	['FTNT','!5Absorbed Radiative Flux','!5W m!E-2!N'] ]
;	['FTNT','!5Absorbed Radiative Flux Response','!5W m!E-2!N'] ]
;	['FTNTDST','!5Tracer Absorbed Radiative Flux','!5W m!E-2!N'] ]
;	['FTNTFRC','!5Tracer TOA Radiative Forcing','!5W m!E-2!N'] ]
;	['FSNS','!5Surface Absorbed Solar Flux','!5W m!E-2!N'] ]
;	['FSNT','!5Absorbed Solar Flux','!5W m!E-2!N']]
;	['FSNT','!5Absorbed Solar Flux Response','!5W m!E-2!N']]
;	['FSNTDST','!5Tracer Absorbed Solar Flux','!5W m!E-2!N'], $
;	['FSNTFRC','!5Tracer TOA Solar Forcing','!5W m!E-2!N']]
;	['FSDS','!5Surface Insolation','!5W m!E-2!N'], $
;	['FSDSDST','!5Tracer Surface Insolation','!5W m!E-2!N'] ]
;	['FSNSDST','!5Tracer Net Surface Solar Flux','!5W m!E-2!N'] ]
;	['FSNSCDST','!5Tracer Net Surface Clearsky Solar Flux','!5W m!E-2!N'] ]
;	['FSNTCDST','!5Tracer Absorbed Clearsky Solar Flux','!5W m!E-2!N'] ]
;	['FSATDST','!5Tracer Atmospheric Absorption','!5W m!E-2!N'], $
;	['FSATC','!5Atmospheric Absorption','!5W m!E-2!N'], $
;	['FSATCDST','!5Tracer Atmospheric Absorption','!5W m!E-2!N'], $
;	['FSNSFRC','!5Tracer Net Surface Solar Forcing','!5W m!E-2!N'], $
;	['FSNSCFRC','!5Tracer Net Surface Clearsky Solar Forcing','!5W m!E-2!N'] ]
;	['FSAT','!5Atmospheric Absorption','!5W m!E-2!N'], $
;	['FSATCFRC','!5Tracer Atmospheric Absorption, Clearsky','!5W m!E-2!N'], $
;	['FSATFRC','!5Tracer Atmospheric Absorption','!5W m!E-2!N'], $
;	['FSDSFRC','!5Tracer Surface Insolation Forcing','!5W m!E-2!N']]
; 	['FSNTCFRC','!5Tracer Absorbed Clearsky Solar Forcing','!5W m!E-2!N'], $
;	['FLNS','!5Net Surface Longwave Cooling','!5W m!E-2!N'], $
;	['FLNT','!5Outgoing Longwave Radiation','!5W m!E-2!N']]
;	['FLNTC','!5Clear Sky OLR','!5W m!E-2!N'], $
;	['FLNS','!5Net Surface IR Emission','!5W m!E-2!N'] ]
;	['FLNTDST','!5Tracer Absorbed IR Flux','!5W m!E-2!N'], $
;	['FLDS','!5Surface IR Down Flux','!5W m!E-2!N'], $
;	['FLDSDST','!5Tracer IR Down Flux','!5W m!E-2!N'] ]
;	['FLNSDST','!5Tracer Net Surface IR Emission','!5W m!E-2!N'] ]
;	['FLNSCDST','!5Tracer Net Surface Clearsky IR Emission','!5W m!E-2!N'] ]
;	['FLNTCDST','!5Tracer Absorbed Clearsky IR Flux','!5W m!E-2!N'] ]
;	['FLATDST','!5Tracer Atmospheric IR Absorption','!5W m!E-2!N'], $
;	['FLATC','!5Atmospheric IR Absorption','!5W m!E-2!N'], $
;	['FLATCDST','!5Tracer Atmospheric IR Absorption','!5W m!E-2!N'], $
;	['FLNSFRC','!5Tracer Net Surface IR Forcing','!5W m!E-2!N'], $
;	['FLNSCFRC','!5Tracer Net Surface Clearsky IR Forcing','!5W m!E-2!N'] ]
;	['FLAT','!5Atmospheric IR Absorption','!5W m!E-2!N'], $
;	['FLATFRC','!5Tracer Atmospheric IR Forcing','!5W m!E-2!N'] ]
;	['FLATCFRC','!5Tracer Atmospheric IR Forcing, Clearsky','!5W m!E-2!N'], $
;	['FLNTFRC','!5Tracer OLR Forcing','!5W m!E-2!N'] ]
;	['FLDSFRC','!5Tracer Surface IR Down Flux Forcing','!5W m!E-2!N'], $
;	['MPCH2O2','!5!5(H!I2!NO)!I2!N Dimer Loading','!5mg m!E-2!N'] ]
;	['NPCO2O2','!5O!I2!N!9. !5O!I2!N Column','!9X!5 10!E42!N molecule!E2!N cm!E-5!N'] ]
;	['NPCO2N2','!5O!I2!N!9. !5N!I2!N Column','!9X!5 10!E42!N mlc!E2!N cm!E-5!N'] ]
fld_nbr2=n_elements(fld2_stf(0,*))
foo={fld2_sct,idx:0,nm:'',sng:'',unit:''}
fld2=replicate({fld2_sct},fld_nbr2)
for idx=0,fld_nbr2-1 do begin
	fld2(idx).idx=idx
	fld2(idx).nm=fld2_stf(0,idx)
	fld2(idx).sng=fld2_stf(1,idx)
	fld2(idx).unit=fld2_stf(2,idx)
endfor; end loop over fld2s

xpt_stf=[ $
;	['dmr03','!5']]
;	['dmr04','!5']]
;	['dmr05','!5']]
;	['dmr06','!5']]
;	['dmr08','!5']]
;	['map','!5DEAD']]
;	['lsm','!5LSM']]
;	['toms','!5TOMS']]
;	['avhrr','!5AVHRR']]
;	['dstccm97','!5Present Day']]
;	['dstccm96','!5Crowley LGM'], $
;	['dstccm95','!5Adams LGM']]
;	['dstccm94','!5Present Day']]
;	['dstccm93','!5Crowley LGM'], $
;	['dstccm92','!5Adams LGM']]
;	['dstccm51','!5dstccm51']]
;	['dstmch45','!5dstmch45']]
;	['dstmch90','!5dstmch90']]
;	['dstmch90','!5DEAD 1990s']]
	['sncpd05','!5']]
;	['dstmch90','!5Uniform (4%)']]
;	['dstmchx3','!5Heterogeneous']]
;	['dstmch81','!5dstmch81']]
;	['ecmwf','!5ECMWF']]
;	['erbe_b','!5ERBE']]
;	['ssmi','!5SSMI'], $
;	['gpcp','!5GPCP'], $
;	['legates','!5Legates'], $
;	['422','!5CCM2'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
;	['amip5','!5CCM'], $
;	['sld012d','!5CCM3']]
;	['spcp_85','!5ANV']]
;	['erbe_b','!5ERBE']]
;	['tef95','!5Tegen & Fung 1995']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

;mth_sng=['06']
;mth_sng=['01','07']
;mth_sng=['01','02','03','04','05','06','07','08','09','10','11','12']
;mth_sng=['01']
;mth_sng=['06']
;mth_sng=['01','02','03','04','05','06','07','08','09']
;mth_sng=['09','10','11','12']
;mth_sng=['05','06','07']
;mth_sng=['1202','0305','0608','0911']
;mth_sng=['1202']
;mth_sng=['0305']
;mth_sng=['0608']
;mth_sng=['0911']
;mth_sng=['1202','0608']
;mth_sng=['1202','0305']
;mth_sng=['0608','0911']
;mth_sng=['05']
;mth_sng=['09']
;mth_sng=['12','01','02']
;mth_sng=['10']
;mth_sng=['07']
;mth_sng=['8589']
;mth_sng=['1998']
mth_sng=[''] ; Use '' when plotting annual and climatological means
mth_nbr=n_elements(mth_sng)
if mth_nbr eq 0 or fld_nbr2 eq 0 or xpt_nbr eq 0 then print,'Loop size is 0 in gcm_xy_bch'

;yyyymmdd_srt=19980101
;yyyymmdd_srt=19990208
yyyymmdd_srt=19980424
day_nbr=1

yr_srt=2001
yr_nbr=1
yr=indgen(yr_nbr)
for yr_idx=0,yr_nbr-1 do begin
	yr(yr_idx)=yr_srt+yr_idx
	yr_sng=string(format='(I4.4)',yr)
endfor; end loop over yr

; 0 <= Quality <= 100, default is 50
if dvc eq 'mpg' then mpg_id=mpeg_open([!d.x_size,!d.y_size],quality=50) ; Open MPEG movie

if dly then time_nbr=day_nbr
if mth then time_nbr=mth_nbr
if yrl then time_nbr=yr_nbr
img_idx=-1 ; Initialize counter for img
for time_idx=0,time_nbr-1 do begin
for fld_idx2=0,fld_nbr2-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	img_idx=img_idx+1 ; Initialize counter for img
	if dly then day_idx=time_idx
	if yrl then yr_idx=time_idx
	if yrl then mth_idx=0
	if mth then mth_idx=time_idx
	if mth then yr_idx=0
	fld_nm=fld2(fld_idx2).nm
	xpt_nm=xpt(xpt_idx).nm
	xpt_sng=xpt(xpt_idx).sng
	fld_sng=fld2(fld_idx2).sng
;	Replace chemical names, if necessary
	fld_sng=trc_nm_sbs(fld_sng,xpt_nm)
	ttl=''
	if fld_sng ne '' then ttl=xpt_sng+' '+fld_sng
	if fld2(fld_idx2).unit ne '' then ttl=ttl+' ('+fld2(fld_idx2).unit+')'
	y_ttl=''
	x_ttl=''
	yr_sng='clm'
	ctl_xpt='dstccm94'
;	rng_x=[-180,180]
;	rng_y=[-90,90]
;	rng_x=[80,100] ; Tibetan Plateau
;	rng_y=[30,40] ; Tibetan Plateau
	rng_x=[60,120] ; SE Asia
	rng_y=[20,50] ; SE Asia
;	rng_x=[-10,30] ; Bodele
;	rng_y=[5,25] ; Bodele
;	rng_x=[-90,0] ; TNA
;	rng_y=[0,30] ; TNA
;	rng_x=[100,-100] ; Asian Dust event
;	rng_y=[20,60] ; Asian Dust event
;	rng_x=[-90,30] ; Zen00 LITE
;	rng_y=[10,35] ; Zen00 LITE
;	rng_x=[40,100] ; CRE00 INDOEX
;	rng_y=[-25,25] ; CRE00 INDOEX
;	rng_x=[-60,-10] ; Mah00 North Atlantic
;	rng_y=[-10,45] ; Mah00 North Atlantic
	if dvc eq 'eps' then fl_out_sfx='.eps'
	if dvc eq 'ps' then fl_out_sfx='.ps'
	if dvc eq 'jpg' then fl_out_sfx='.jpg'
	if dvc eq 'mpg' then fl_out_sfx='.mpg'
	if dvc eq 'png' then fl_out_sfx='.png'
	if dvc eq 'ppm' then fl_out_sfx='.ppm'
	pre_scl=0
	cbar=1
	fl_clr_tbl=''
	clr_tbl_nbr=22
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	lft=1
	dsh=0
	if xpt_nm eq 'lsm' then yr_sng(yr_idx)='7993_'
;	if xpt_nm eq 'erbe_b' and yr_sng ne '' then yr_sng(yr_idx)='8589_'
	if mth then begin
		if strlen(mth_sng(mth_idx)) gt 2 then date_sng=yr_sng(yr_idx)+'_'+mth_sng(mth_idx) else date_sng=yr_sng(yr_idx)+mth_sng(mth_idx)
	endif; endif mth
	if yrl then begin
		date_sng=yr_sng(yr_idx)+mth_sng(mth_idx)
		if ttl ne '' then ttl=ttl+' '+date_sng
	endif; endif yrl
	if dly then begin
		yyyymmdd=yyyymmdd_ncr(yyyymmdd_srt,day_idx)
		yyyymmdd_sng=string(format='(I8.8)',yyyymmdd)
		date_sng=yyyymmdd_sng
		if ttl ne '' then ttl=ttl+' '+date_sng
	endif; endif dly
	shd_val=0.0
	shd_xpr='lt'
        ; fxm: Scale mass fluxes, optical depths and forcing for whole experiment
        ; NB: This is for kludges to quantities in linear limit
        ; Always tie this to CASEID and field name
	if strpos(fld_nm,'DST') ne -1 or strpos(fld_nm,'FRC') ne -1 and xpt_nm eq 'dstccm92' or xpt_nm eq 'dstccm93' or xpt_nm eq 'dstccm94' then begin
            ; fxm: Scale mass fluxes and optical depths of whole experiment
		scl=1.67
        endif ; if scaling whole experiment
	if strpos(fld_nm,'DST') ne -1 or strpos(fld_nm,'FRC') ne -1 and xpt_nm eq 'dstccm95' or xpt_nm eq 'dstccm96' or xpt_nm eq 'dstccm97' then begin
            ; fxm: Scale mass fluxes and optical depths of whole experiment
		scl=1.2
        endif ; if scaling whole experiment
	if fld_nm eq 'ALBEDO' then pre_scl=[0.1,0.1,8]
	if fld_nm eq 'FSUT' then if not diff then pre_scl=[50,10,13] else pre_scl=[-50,10,11]
	if fld_nm eq 'FLUS' then if not diff then pre_scl=[300,20,10] else pre_scl=[-50,10,11]
	if fld_nm eq 'NCF' then if not diff then pre_scl=[-50,10,11] else pre_scl=[-30,5,13] ;pre_scl=[-70,10,15]
	if fld_nm eq 'NCF' then fl_clr_tbl='~/idl/delta.tbl'
	if fld_nm eq 'GCLD' then pre_scl=[100,20,11]
	if fld_nm eq 'GCLD' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'GCLR' then pre_scl=[10,20,11]
	if fld_nm eq 'GCLR' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'NET' then if not diff then pre_scl=[10,10,18] else pre_scl=[-50,10,11]
	if fld_nm eq 'FLNS' then if not diff then pre_scl=[10,10,18] else pre_scl=[-50,10,11]
	if fld_nm eq 'FSNS' then if not diff then pre_scl=[120,20,11] else pre_scl=[-50,10,11]
;	if fld_nm eq 'FSNT' then if not diff then pre_scl=[160,40,5] else pre_scl=[-50,10,11]
;	if fld_nm eq 'FSNT' then if not diff then pre_scl=[150,50,5] else pre_scl=[-50,10,11]
;	if fld_nm eq 'FSNT' then if not diff then pre_scl=[0,20,31] else pre_scl=[-50,10,11]
	if fld_nm eq 'FSNT' then if not diff then pre_scl=[0,20,31] else pre_scl=[-30,2,31]
	if fld_nm eq 'FTNT' then if not diff then pre_scl=[-225,15,31] else pre_scl=[-20,2,21]
	if fld_nm eq 'FSAT' then if not diff then pre_scl=[0,10,11] else pre_scl=[-50,10,11]
	if fld_nm eq 'FLNT' then if mth_sng(mth_idx) eq '01' then pre_scl=[120,20,11] else pre_scl=[100,25,10] 

;	if fld_nm eq 'FLNT' then pre_scl=[120,20,10] else pre_scl=[100,25,10] 
	if fld_nm eq 'FLNT' then fl_clr_tbl='~/idl/standard.tbl'
	if fld_nm eq 'FLNTC' then if mth_sng(mth_idx) eq '01' then pre_scl=[140,10,19] else pre_scl=[90,20,13] 
	if fld_nm eq 'FLNTC' then fl_clr_tbl='~/idl/standard.tbl'
	if fld_nm eq 'TMQ' then pre_scl=[0,5,13]
	if fld_nm eq 'Q2' then pre_scl=[1,0.25,11]
	if fld_nm eq 'Q2' then scl=1.0e4
	if fld_nm eq 'TMQ' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'TS1' then if not diff then pre_scl=[220,10,10] else pre_scl=[-7,1,15]
	if fld_nm eq 'TS1' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'PSL' then scl=0.01
	if fld_nm eq 'PSL' then if not diff then pre_scl=[950,10,10] else pre_scl=[-10,1,21]
	if fld_nm eq 'PSL' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'PS' then scl=0.01
	if fld_nm eq 'PS' then if not diff then pre_scl=[950,10,10] else pre_scl=[-8,1,17]
	if fld_nm eq 'PS' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'PRECT' then if not diff then pre_scl=[0,2,13] else pre_scl=[-10,1,21]
	if fld_nm eq 'PRECT' and xpt_nm ne 'gpcp' and xpt_nm ne 'cmap' then scl=8.64e7
	if fld_nm eq 'PRECT' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'PRECC' then if not diff then pre_scl=[0,2,13] else pre_scl=[-10,1,21]
	if fld_nm eq 'PRECC' and xpt_nm ne 'gpcp' then scl=8.64e7
	if fld_nm eq 'PRECC' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'PRECL' then if not diff then pre_scl=[0,2,13] else pre_scl=[-10,1,21]
	if fld_nm eq 'PRECL' and xpt_nm ne 'gpcp' then scl=8.64e7
	if fld_nm eq 'PRECL' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'LWCF' then if not diff then pre_scl=[0,10,13] else pre_scl=[-30,5,13]
	if fld_nm eq 'LWCF' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'RATIO' then if not diff then pre_scl=[.8,.05,9]
	if fld_nm eq 'SWCF' then if not diff then pre_scl=[-200,20,11] else pre_scl=[-30,5,13] ;pre_scl=[-70,10,15]
;	if fld_nm eq 'SWCF' then if not diff then pre_scl=[-90,10,10] else pre_scl=[-30,5,13] ;pre_scl=[-70,10,15]
	if fld_nm eq 'SWCF' then fl_clr_tbl='~/idl/standard.tbl'
	if fld_nm eq 'SHFLX' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'SHFLX' then if not diff then pre_scl=[-60,20,14] else pre_scl=[-50,10,11]
	if fld_nm eq 'SHFLX' then shd_val=40.0
	if fld_nm eq 'SHFLX' then shd_xpr='ge'
	if fld_nm eq 'LHFLX' then fl_clr_tbl='~/idl/aero.tbl'
	if fld_nm eq 'LHFLX' then if not diff then pre_scl=[0,20,16] else pre_scl=[-50,10,11]
	if fld_nm eq 'LHFLX' then shd_val=40.0
	if fld_nm eq 'LHFLX' then shd_xpr='ge'
	if fld_nm eq 'VMAGSFC' then if not diff then pre_scl=[0,2,10] else pre_scl=[-4,1,9]
	if fld_nm eq 'VMAGSFC' then shd_val=6.0
 	if fld_nm eq 'VMAGSFC' then shd_xpr='ge'
	if fld_nm eq 'WND_SFC' then if not diff then pre_scl=[0,2,10] else pre_scl=[-4,1,9]
	if fld_nm eq 'WND_SFC' then shd_val=6.0
	if fld_nm eq 'WND_SFC' then shd_xpr='ge'
	if fld_nm eq 'WND_MBL' then if not diff then pre_scl=[0,2,6] else pre_scl=[-4,1,9]
	if fld_nm eq 'U10' then pre_scl=[-8,2,9]
	if fld_nm eq 'V10' then pre_scl=[-8,2,9]
	if fld_nm eq 'WND_MBL' then shd_val=6.0
	if fld_nm eq 'WND_MBL' then shd_xpr='ge'
	if fld_nm eq 'SMP_SFC' then if not diff then pre_scl=[0,1.0e5,21] else pre_scl=[-4,1,9]
	if fld_nm eq 'PBLH' then if not diff then pre_scl=[0,100,20] else pre_scl=[-200,20,21]
	if fld_nm eq 'PBLH' then shd_val=700.0
	if fld_nm eq 'PBLH' then shd_xpr='ge'
	if fld_nm eq 'TOTLWP' then if not diff then pre_scl=[0,10,15] else pre_scl=[-150,25,13]
	if fld_nm eq 'TOTLWP' then shd_val=40.0
	if fld_nm eq 'TOTLWP' then shd_xpr='ge'
	if fld_nm eq 'TOTIWP' then if not diff then pre_scl=[0,20,11] else pre_scl=[-150,25,13]
	if fld_nm eq 'TOTIWP' then shd_val=40.0
	if fld_nm eq 'TOTIWP' then shd_xpr='ge'
	if fld_nm eq 'TOTCWP' then if not diff then pre_scl=[0,25,13] else pre_scl=[-150,25,13]
	if fld_nm eq 'TOTCWP' then if not diff then shd_val=75
	if fld_nm eq 'TOTCWP' then if not diff then shd_xpr='ge'
	if fld_nm eq 'TOTCWP' then if not diff then order=0
	if fld_nm eq 'TOTCWP' then if not diff then pll=7
;	if fld_nm eq 'Z2TEST' then if not diff then pre_scl=[500,10,9] else pre_scl=[-10,2,11]
	if fld_nm eq 'Z2TEST' then if not diff then pre_scl=[500,10,20] else pre_scl=[-10,2,11]
	if fld_nm eq 'Z2TEST' then scl=0.101977
	if fld_nm eq 'Z2TEST' then fl_clr_tbl='~/idl/bw.tbl'
	if fld_nm eq 'Z2TEST' then if not bias then shd_val=560.0
	if fld_nm eq 'Z2TEST' then if not bias then shd_xpr='ge'
	if fld_nm eq 'U' then if not diff then pre_scl=[-30,5,22] else pre_scl=[-5,1,11]
;	if fld_nm eq 'U' then pre_scl=[-10,1,21]
	if fld_nm eq 'CHI' then if not diff then pre_scl=[-30,3,20] else pre_scl=[-7,1,15]
	if fld_nm eq 'CHI' then scl=1.0e-6
	if fld_nm eq 'PSI' then if not diff then pre_scl=[-200,20,21] else pre_scl=[-20,2,21]
	if fld_nm eq 'PSI' then scl=1.0e-6
	if fld_nm eq 'DIV' then if not diff then pre_scl=[-10,1,21] else pre_scl=[-6,1,13]
	if fld_nm eq 'DIV' then scl=1.0e6
	if fld_nm eq 'DIV' then fl_clr_tbl='~/idl/delta.tbl'
;	if fld_nm eq 'VOR' then pre_scl=[-15,5,19]
	if fld_nm eq 'VOR' then scl=1.0e6
	if fld_nm eq 'fracla' then scl=1.0e2
	if fld_nm eq 'fracla' then pre_scl=[0,10.0,11]
	if fld_nm eq 'mss_frc_CaCO3' then scl=1.0e2
	if fld_nm eq 'mss_frc_CaCO3' then pre_scl=[0,2.0,12]
 	if fld_nm eq 'aotm14' then pre_scl=[0,0.05,11]
	if fld_nm eq 'src_frq' then scl=1.0e2
	if fld_nm eq 'src_frq' then pre_scl=[0,10.0,11]
	if fld_nm eq 'csn_cyc' then scl=1.0e2
	if fld_nm eq 'csn_cyc' then pre_scl=[0,10.0,11]
	if fld_nm eq 'fsh_fct' then pre_scl=[0,0.1,11]
	if fld_nm eq 'vwc_sfc' then pre_scl=[0,0.02,16]
	if fld_nm eq 'FRC_WET' then pre_scl=[0.95,0.05,21]
	if fld_nm eq 'VWC_SFC' then pre_scl=[0,0.05,11]
	if fld_nm eq 'VAI_DST' then pre_scl=[0,0.25,21]
	if fld_nm eq 'DPSWETFRC' then pre_scl=[0,5.0,21]
	if fld_nm eq 'DPSWETFRC' then scl=1.0e2
	if fld_nm eq 'aer_idx_331_360' then pre_scl=[0.0,0.1,21]
	if fld_nm eq 'aer_idx_331_360' and dly then pre_scl=[0.0,0.2,16]
	if strpos(fld_nm,'rdb_fct') ne -1 or strpos(fld_nm,'mbl_bsn') ne -1 then begin
		scl=1.0
		pre_scl=[0.0,0.05,21]
               ; Remove 5.707 scale factor often applied after running bds
;                if fld_nm eq 'rdb_fct_gmr' then scl=1.0/5.707
         endif ; endif rdb_fct
	if strpos(fld_nm,'DSTMPC') ne -1 then begin
		scl=1.0e6
		if fld_nm eq 'DSTMPC01' then pre_scl=[0,40,11] else pre_scl=[0,40,11]
	endif ; endif DSTMPC
	if strpos(fld_nm,'DSTSF') ne -1 then begin
		if (fld_nm eq 'DSTSFDPSDST' or fld_nm eq 'DSTSFDPSCA' or fld_nm eq 'DSTSFDPSCACO3') then begin
;			scl=3.1536e12 ; [kg m-2 s-1] --> [mg cm-2 ka-1] CGS units!
			scl=1.0e12 ; [kg m-2 s-1] --> [ng m-2 s-1]
			pre_scl=[0,10,11]
		endif else if fld_nm eq 'DSTSFDPSC' then begin
			scl=3.1536e15 ; [kg m-2 s-1] --> [ug cm-2 ka-1]
			pre_scl=[0,10,11]
		endif else if fld_nm eq 'DSTSFDPSRATCACO3' then begin
			scl=1000.0 ; 
			pre_scl=[0,0.1,11]
		endif else if fld_nm eq 'DSTSFMBL' or fld_nm eq 'DSTSFDRY' or fld_nm eq 'DSTSFPCP' or fld_nm eq 'DSTSFDPS' then begin
			scl=1.0e9 ; [kg m-2 s-1] --> [ug m-2 s-1]
			pre_scl=[0,1,21]
		endif else begin ; endif DSTSFMBL
			; All other surface fluxes of dust
			scl=1.0e12 ; [kg m-2 s-1] --> [ng m-2 s-1]
			pre_scl=[0,10,11]
		endelse; endelse
	endif ; endif dust
	if fld_nm eq 'AEROD_v' then pre_scl=[0,0.03,21]
	if strpos(fld_nm,'FRC') ne -1 and strpos(xpt_nm,'snc') ne -1 then begin
		if fld_nm eq 'AERAFRC_TOP' then pre_scl=[-10,1.0,15]
		if fld_nm eq 'AERAFRC_ATM' then pre_scl=[0,1.0,17]
		if fld_nm eq 'AERAFRC_SFC' then pre_scl=[0.0,1.0,17]
		if fld_nm eq 'SNOAERFRC' then pre_scl=[0,1.0,17]
	endif ; endif SNICAR radiative field
	if strpos(fld_nm,'FRC') ne -1 and strpos(xpt_nm,'dst') ne -1 then begin
		; SW fields
		; Annual
		if fld_nm eq 'FSNTFRC' then pre_scl=[-20,1.0,31]
		if fld_nm eq 'FTNTFRC' then pre_scl=[-18,1.0,31]
;		if fld_nm eq 'FSNTFRC' then pre_scl=[-20,1.5,21]
		if fld_nm eq 'FSATFRC' then pre_scl=[0,2.0,11]
		if fld_nm eq 'FSDSFRC' then pre_scl=[-40,4.0,11]
		; Seasonal/monthly
;		if fld_nm eq 'FSNTFRC' then pre_scl=[0,2.0,11]
;		if fld_nm eq 'FSATFRC' then pre_scl=[0,1.5,11]
;		if fld_nm eq 'FSDSFRC' then pre_scl=[-20,2.0,11]
		if fld_nm eq 'FSATCFRC' then pre_scl=[0,3.0,11]
		if fld_nm eq 'FSNTCFRC' then pre_scl=[0,2.0,11]
		if fld_nm eq 'FSNSFRC' then pre_scl=[-20,2.0,11]
		if fld_nm eq 'FSNSCFRC' then pre_scl=[-20,2.0,11]
		if fld_nm eq 'FSDSFRC' then order=0
		if fld_nm eq 'FSNSFRC' then order=0
		if fld_nm eq 'FSNSCFRC' then order=0
		; LW fields
		if fld_nm eq 'FLDSFRC' then pre_scl=[0.0,1.0,11]
	endif ; endif dust radiative field
	if fld_nm eq 'MPCH2O2' then scl=1.0e3
	if fld_nm eq 'MPCH2O2' then pre_scl=[0,3,11]
  	if fld_nm eq 'NPCO2O2' then scl=1.0e-16
	if fld_nm eq 'NPCO2O2' then pre_scl=[3,1,13]
	if fld_nm eq 'NPCO2N2' then scl=1.0e-16
	if fld_nm eq 'NPCO2N2' then pre_scl=[15,5,9]
	if strpos(fld_nm,'FRC') ne -1 and strpos(xpt_nm,'dmr') ne -1 then begin
		if xpt_nm eq 'dmr03' then begin
			if fld_nm eq 'FSATFRC' then pre_scl=[0,0.3,11]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0,0.3,11]
		endif; not dmr03 H2OH2O
		if xpt_nm eq 'dmr04' then begin
		if yr_sng(yr_idx) eq '' then begin
			if fld_nm eq 'FSATFRC' then pre_scl=[0.3,0.1,11]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.3,0.1,11]
		endif else begin; not annual
			if fld_nm eq 'FSATFRC' then pre_scl=[0.0,0.2,11]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.0,0.2,11]
		endelse; not annual
		endif; not dmr04 O2O2+O2N2
		if xpt_nm eq 'dmr05' then begin
		if yr_sng(yr_idx) eq '' then begin
			if fld_nm eq 'FSATFRC' then pre_scl=[0.06,0.02,9]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.06,0.02,9]
		endif else begin; not annual
			if fld_nm eq 'FSATFRC' then pre_scl=[0.0,0.04,11]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.0,0.04,11]
		endelse; not annual
		endif; not dmr05 O2N2
		if xpt_nm eq 'dmr06' then begin
		if yr_sng(yr_idx) eq '' then begin
			if fld_nm eq 'FSATFRC' then pre_scl=[0.3,0.1,8]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.3,0.1,8]
		endif else begin; not annual
			if fld_nm eq 'FSATFRC' then pre_scl=[0.0,0.2,11]
			if fld_nm eq 'FSATCFRC' then pre_scl=[0.0,0.2,11]
		endelse; not annual
		endif; not dmr06 O2O2

		if fld_nm eq 'FSNTFRC' then pre_scl=[0,0.2,11]
		if fld_nm eq 'FSNTCFRC' then pre_scl=[0,0.2,11]

		if fld_nm eq 'FSDSFRC' then pre_scl=[-2,0.2,11]
		if fld_nm eq 'FSNSFRC' then pre_scl=[-2,0.2,11]
		if fld_nm eq 'FSNSCFRC' then pre_scl=[-2,0.2,11]

		if fld_nm eq 'FSDSFRC' then order=0
		if fld_nm eq 'FSNSFRC' then order=0
		if fld_nm eq 'FSNSCFRC' then order=0
	endif ; endif dimer radiative field
	if strpos(fld_nm,'DSTODX') ne -1 then begin
;		if diff then pre_scl=[-0.1,0.02,11] else pre_scl=[0.0,0.10,11]
		if diff then pre_scl=[-0.2,0.04,11] else pre_scl=[0.0,0.05,11]
;		pre_scl=[0.0,0.025,11]
	endif ; endif dust
	if fld_nm eq 'dust_ttl' then pre_scl=[0,10,11]
	if fld_nm eq 'dust_ttl' then scl=1.0e9
;	if strpos(xpt_nm,'ecmwf') ne -1 then begin
;		if mth_sng(mth_idx) eq '01' then yr_sng(yr_idx)='9095_' else yr_sng(yr_idx)='8994_'
;	endif
	if strpos(xpt_nm,'ssmi') ne -1 then begin
		if mth_sng(mth_idx) eq '01' then yr_sng(yr_idx)='8795_' 
		if mth_sng(mth_idx) eq '07' then yr_sng(yr_idx)='8791_' 
	endif
	if strpos(xpt_nm,'legates') ne -1 then yr_sng(yr_idx)='2080_'
	if strpos(xpt_nm,'gpcp') ne -1 then begin
		if mth_sng(mth_idx) eq '01' then yr_sng(yr_idx)='8894_' 
		if mth_sng(mth_idx) eq '07' then yr_sng(yr_idx)='8793_' 
	endif
	pth_data_sng=getenv('DATA')
;	pth_ps_sng='/fs/cgd/data0/zender/ps'
	pth_ps_sng=getenv('DATA')+'/ps'
;	if mth or yrl then date_sng=yr_sng(yr_idx)+mth_sng(mth_idx)
	if fl_usr ne '' then fl_nm_in=fl_usr else fl_nm_in=pth_data_sng+'/'+xpt_nm+'/'+xpt_nm+'_'+prs_sng+date_sng+'.nc'
	fl_nm_out=pth_ps_sng+'/'+xpt_nm+'_'+prs_sng+date_sng+'_'+rgn_sng+fld_nm+fl_out_sfx
	if dvc eq 'mpg' then fl_nm_out=pth_ps_sng+'/'+xpt_nm+'_'+prs_sng+rgn_sng+fld_nm+fl_out_sfx
	if diff then begin ; plot differences between simulated fields
		ttl='dstccm97 - dstccm94 '+fld_sng
		if fld2(fld_idx2).unit ne '' then ttl=ttl+' ('+fld2(fld_idx2).unit+')'
		fl_clr_tbl='~/idl/delta9.tbl'
		fl_nm_in=pth_data_sng+'/'+xpt_nm+'/'+xpt_nm+'_'+ctl_xpt+'_'+yr_sng(yr_idx)+prs_sng+mth_sng(mth_idx)+'.nc'
		fl_nm_out=pth_ps_sng+'/'+xpt_nm+'_'+yr_sng(yr_idx)+'_'+ctl_xpt+'_'+yr_sng(yr_idx)+prs_sng+mth_sng(mth_idx)+'_'+rgn_sng+fld_nm+fl_out_sfx
	endif; endif diff
	if bias then begin ; plot biases of simulated fields relative to observational analyses
	if strpos(xpt_nm,'amip5') ne -1 or strpos(xpt_nm,'spcp') ne -1 then begin
		fl_clr_tbl='~/idl/delta9.tbl'
		ttl=xpt_sng+'-'+xpt(0).sng+' '+fld_sng
		if fld2(fld_idx2).unit ne '' then ttl=ttl+' ('+fld2(fld_idx2).unit+')'
		if strpos(xpt(0).nm,'ecmwf') ne -1 then if mth_sng(mth_idx) eq '01' then bias_yr_sng(yr_idx)='9095_' else bias_yr_sng(yr_idx)='8994_'
		if strpos(xpt(0).nm,'legates') ne -1 then bias_yr_sng(yr_idx)='2080_'
		if fld_nm eq 'Z2TEST' then pre_scl=[-20,2,21]
		if fld_nm eq 'CHI' then pre_scl=[-10,1,21]
		if fld_nm eq 'PRECT' then pre_scl=[-10,1,21]
		if fld_nm eq 'TS1' then pre_scl=[-10,1,21]
		if fld_nm eq 'U' then pre_scl=[-40,5,17]
		if fld_nm eq 'T' then pre_scl=[-15,1,30]
		fl_nm_in=pth_data_sng+'/'+xpt_nm+'/'+xpt_nm+'_'+yr_sng(yr_idx)+xpt(0).nm+'_'+bias_yr_sng(yr_idx)+prs_sng+mth_sng(mth_idx)+'.nc'
		fl_nm_out=pth_ps_sng+'/'+xpt_nm+'_'+yr_sng(yr_idx)+xpt(0).nm+'_'+bias_yr_sng(yr_idx)+prs_sng+mth_sng(mth_idx)+'_'+rgn_sng+fld_nm+fl_out_sfx
	endif; endif xpt is a simulation, not ecmwf
	endif; endif bias
	if lbl_typ eq 'auto' then begin
            if xpt_nbr gt 1 then pnl_lbl=pnl_lbl_get(xpt_idx) 
            if time_nbr gt 1 then pnl_lbl=pnl_lbl_get(time_idx) 
            if fld_nbr2 gt 1 then pnl_lbl=pnl_lbl_get(fld_idx2)
        endif ; endif auto panel label
	if prn then begin
		ttl=xpt_sng+' '+fld_sng
		if fld2(fld_idx2).unit ne '' then ttl=ttl+' ('+fld2(fld_idx2).unit+')'
		if dly or yrl then begin
			if ttl ne '' then ttl=ttl+' '+date_sng
		endif; endif dly
		if mth then begin
			if mth_sng(mth_idx) eq '01' then if xpt_nbr gt 1 then ttl='Jan '+ttl
			if mth_sng(mth_idx) eq '07' then if xpt_nbr gt 1 then ttl='Jul '+ttl
		endif; endif mth
		if fld_nbr2 eq 4 then begin
			top=1
			if fld_idx2 eq 2 or fld_idx2 eq 3 then btm=1 else btm=0
			; Jump to print routines
			goto,prn_set
		endif; endif fld_top
		if fld_top then begin
			if fld_idx2 eq fld_nbr2-1 then btm=1 else btm=0
			if mth_sng(mth_idx) eq '01' then ttl='Jan '+ttl else if mth_sng(mth_idx) eq '07' then ttl='Jul '+ttl
		endif; endif fld_top
		if not fld_top then begin
			if xpt_idx eq 0 and ((xpt_nbr gt 1) or (xpt_nbr eq 1 and time_idx eq 0)) then top=1 else top=0
			if xpt_idx eq xpt_nbr-1 and ((xpt_nbr gt 1) or (xpt_nbr eq 1 and time_idx eq time_nbr-1)) then btm=1 else btm=0
		endif; endif not fld_top
		if dly or yrl and time_nbr gt 2 then begin
			if time_idx eq 0 or time_idx eq time_nbr/2 then top=1 else top=0
			if time_idx eq time_nbr-1 or time_idx eq (time_nbr/2)-1 then btm=1 else btm=0
		endif; endif dly
		if mth and mth_nbr eq 4 and csn then begin
			if (mth_sng(mth_idx) eq '1202' or mth_sng(mth_idx) eq '0608') then begin
				top=1
				btm=0
			endif; endif top panel
			if (mth_sng(mth_idx) eq '0305' or mth_sng(mth_idx) eq '0911') then begin
				top=0
				btm=1
			endif; endif top panel
		endif; endif csn
		if mth and mth_nbr eq 12 then begin
			if (mth_idx eq 5 or mth_idx eq 11) then btm=1 else btm=0
			if (mth_idx eq 0 or mth_idx eq 6) then top=1 else top=0
		endif; endif mth_nbr eq 12
		if azi then top=0
		if azi then btm=0
	endif; endif prn
	if strpos(fld_nm,'DST') ne -1 then begin
		if xpt_nbr eq 2 or xpt_nbr eq 3 then begin
			top=1
			if xpt_idx eq xpt_nbr-1 then cbar=1 else cbar=0
			if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
		endif ; endif xpt_nbr
	endif ; endif DSTSFDPS
	; Jump to here if ready to print
	prn_set: foo=1
	if btm then cbar=1 else cbar=0
	if jtk then begin
		pll=1
		ttl=xpt_sng+' '+fld_sng+' ('+fld2(fld_idx2).unit+')'
		if strpos(xpt_nm,'erbe_b') ne -1 then begin
			fl_nm_in=getenv('DATA')+'/erbe_b/erbe_b_8589.nc'
			fl_nm_out=getenv('DATA')+'/ps/erbe_b_8589_'+rgn_sng+fld_nm+fl_out_sfx
		endif
		if strpos(xpt_nm,'sld012d') ne -1 then begin
			fl_nm_in=getenv('DATA')+'/sld012d/sld012d_7993.nc'
			fl_nm_out=getenv('DATA')+'/ps/sld012d_7993_'+rgn_sng+fld_nm+fl_out_sfx
		endif
		if fld_nm eq 'LWCF' then pre_scl=[0,15,7]
		if fld_nm eq 'SWCF' then pre_scl=[-90,15,8]
		if fld_nm eq 'FSNT' then pre_scl=[80,40,8]
		if fld_nm eq 'FLNT' then pre_scl=[135,15,12]
		if fld_nm eq 'FLNTC' then pre_scl=[135,15,12]
		if fld_nm eq 'SWCF' then fl_clr_tbl='~/idl/standard2.tbl'
		if fld_nm eq 'SWCF' then begin
;			shd_val=110.0
;			shd_xpr='ge'
		endif; endif SWCF
		top=1
		btm=1
		cnt_ovl=1
		dsh=0
	endif; endif jtk
	if bg then begin
		cbar=0
		top=0
		btm=0
		lft=0
;		clr_tbl_nbr=22
		order=0
		fl_clr_tbl='~/idl/bg.tbl'
		if fld_nm eq 'SWCF' then pre_scl=[-140,10,15]
		if fld_nm eq 'FSNS' then pre_scl=[100,20,11]
;		if strpos(xpt_nm,'sld012d') ne -1 then begin
			fl_nm_in=getenv('DATA')+'/sld012d/sld012d_8589.nc'
			fl_nm_out=getenv('DATA')+'/ps/bg'+fl_out_sfx
;		endif
	endif; endif bg
	if ams then begin
;		pll=7
;		fl_clr_tbl=''
;		clr_tbl_nbr=0
		rng_y=[-30,30]
;		cnt_ovl=1
		dsh=0
		cbar=0
		top=1
		btm=1
		if fld_nm eq 'LWCF' then if not diff then pre_scl=[0,10,13] else pre_scl=[-40,10,9]
	endif; endif ams
	if prn then begin
		x_sz=6.5
		y_sz=3.2
		if azi then y_sz=6.5
		if top then y_sz=y_sz*1.08
		if btm then y_sz=y_sz*1.08
		if cbar then y_sz=y_sz*1.13
		if dvc eq 'eps' then open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps 
		if dvc eq 'ps' then open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color
	endif; endif prn
	gcm_xy,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,pll=pll,order=order,fld_nm=fld_nm,cbar=cbar,wrp=1,rng_x=rng_x,rng_y=rng_y,pre_scl=pre_scl,chr_sz=chr_sz,dsh=dsh,top=top,lft=lft,btm=btm,scl=scl,img=img,cnt_ovl=cnt_ovl,clr_tbl_nbr=clr_tbl_nbr,fl_clr_tbl=fl_clr_tbl,lvl_idx=lvl_idx,shd_val=shd_val,shd_xpr=shd_xpr,azi=azi,time_idx=time_idx,dvc=dvc,pnl_lbl=pnl_lbl
	if prn then close_ps,fl_nm=fl_nm_out
	if dvc eq 'mpg' then begin
		print,'Writing MPEG frame image to ',fl_nm_out
		image=tvrd(true=1)
		mpeg_put,mpg_id,image=image,frame=img_idx,/order
	endif; endif dvc eq 'mpg'
	if dvc eq 'jpg' then begin
		print,'Saving JPEG format image to ',fl_nm_out
		image=tvrd(true=1)
		write_jpeg,fl_nm_out,image,true=1
	endif; endif dvc eq 'jpg'
	if dvc eq 'png' then begin
		; Warning: write_png() does not work
		print,'Saving PNG format image to ',fl_nm_out
		image=tvrd(true=1)
		write_png,fl_nm_out,image
	endif; endif dvc eq 'png'
	if dvc eq 'ppm' then begin
		; Warning: write_ppm() only works on 2D rasters, so 3-channel images (True Color displays) do not work with it
		print,'Saving PPM format image to ',fl_nm_out
		image=tvrd()
		write_ppm,fl_nm_out,image
	endif; endif dvc eq 'ppm'
endfor; end loop over experiment
endfor; end loop over fld2
endfor; end loop over mth
if dvc eq 'mpg' then begin
	mpeg_save,mpg_id,filename=fl_nm_out ; Close MPEG movie file
	mpeg_close,mpg_id ; Close MPEG movie object
endif; endif mpg
end; end gcm_xy_bch()
pro gcm_xy, $
	btm=btm, $
	cbar=cbar, $
	chr_sz=chr_sz, $
	clr_tbl_nbr=clr_tbl_nbr, $
	cnt_ovl=cnt_ovl, $
	dsh=dsh, $
	dvc=dvc, $ ; Output device (ps, eps, jpg, mpg, ppm ...)
	fl_clr_tbl=fl_clr_tbl, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	img=img, $
	azi=azi, $
	lvl_idx=lvl_idx, $
	mth_end=mth_end, $
	mth_srt=mth_srt, $
	order=order, $
	pll=pll, $
	pre_scl=pre_scl, $
	prn=prn, $
	scl=scl, $
	top=top, $
	lft=lft, $
	ttl=ttl, $
	shd_val=shd_val, $
	shd_xpr=shd_xpr, $
	wrp=wrp, $
	rng_x=rng_x, $
	x_ttl=x_ttl, $
	rng_y=rng_y, $
	pnl_lbl=pnl_lbl, $
	time_idx=time_idx, $
	y_ttl=y_ttl

; gcm_xy,ttl='!5ANV LWCF (W m!E-2!N)',pll=4,order=1
; gcm_xy,fl_nm=getenv('DATA')+'/dst16/dst16_8589_01.nc',ttl='',prn=0,pll=1,order=1,fld_nm='ORO',cbar=0,wrp=1,rng_x=[-180,180],rng_y=[-90,90],top=1,btm=1,img=0,cnt_ovl=0
; gcm_xy,fl_nm=getenv('DATA')+'/ecmwf/ecmwf_prs_9095_01.nc',ttl='Temperature',prn=0,pll=1,order=1,fld_nm='T',cbar=cbar,wrp=1,rng_x=[-180,180],rng_y=[-90,90],top=1,btm=1,img=1,cnt_ovl=0,lvl_idx=9
; gcm_xy,fl_nm=getenv('DATA')+'/ssmi/ssmi_8795_01.nc',ttl='Liquid Water Path',prn=0,pll=1,order=1,fld_nm='TOTLWP',cbar=cbar,wrp=1,rng_x=[-180,180],rng_y=[-90,90],top=1,btm=1,img=1,cnt_ovl=0,pre_scl=[0,20,11]

@ibp_clr.com

if n_elements(pre_scl) eq 3 then begin
	auto_scl=0
	cntr_lvl_min=pre_scl(0)
	cntr_ntv=pre_scl(1)
	cntr_lvl_nbr=pre_scl(2)
endif else auto_scl=1 
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(dsh) eq 0 then dsh=1
if n_elements(shd_val) eq 0 then shd_val=0.0
if n_elements(shd_xpr) eq 0 then shd_xpr='lt'
if n_elements(wrp) eq 0 then wrp=1 ; What does wrap do?
if n_elements(cbar) eq 0 then cbar=1
if n_elements(dvc) eq 0 then dvc='eps' ; Output device (ps, eps, jpg, mpg, ppm ...)
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ERBE LWCF (W m!E-2!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(lft) eq 0 then lft=1
if n_elements(btm) eq 0 then btm=1
if n_elements(img) eq 0 then img=0
if n_elements(azi) eq 0 then azi=0
if n_elements(lvl_idx) eq 0 then lvl_idx=-1
if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=22
if n_elements(fl_clr_tbl) eq 0 then fl_clr_tbl=''
if n_elements(prn) eq 0 then prn=0
if n_elements(cnt_ovl) eq 0 then cnt_ovl=1
if n_elements(y_ttl) eq 0 then y_ttl=''
if n_elements(x_ttl) eq 0 then x_ttl=''
if n_elements(order) eq 0 then clr_ord=1 else clr_ord=order
if n_elements(pll) eq 0 then pll=3
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/erbe_b/erbe_b_8589_01.nc'
if n_elements(fld_nm) eq 0 then fld_nm='LWCF'
if n_elements(time_idx) eq 0 then time_idx=0
if n_elements(pnl_lbl) eq 0 then pnl_lbl=''

if fl_clr_tbl eq '' then usr_dfn_clr_tbl=0 else usr_dfn_clr_tbl=1
if fl_clr_tbl eq '~/idl/standard2.tbl' then clr_ord=0 
color_order=clr_ord
clr_rfr,clr_tbl_nbr,fl_clr_tbl,usr_dfn_clr_tbl
erase
!p.multi=0

; Check for file existence
rcd=file_test(fl_nm)
if not rcd then begin
	print,'Skipping non-existent file '+fl_nm
	goto,end_of_procedure
endif; endif non-existent
; Read binary data from NetCDF file 
nc_id=ncdf_open(fl_nm)
print,'Processing '+fl_nm

; Get number of dimensions in input file
; var_sct={name:'',id:0L,data_type:'',dmn_nbr:0L,nbr_atts:0L,dmn:lonarr(4)}
var_id=ncdf_varid(nc_id,fld_nm)
var_crr=ncdf_varinq(nc_id,var_id)
dmn_nbr=var_crr.ndims
dmn_id=var_crr.dim
srt=0*lonarr(dmn_nbr) ; [idx] Starting offset
cnt=1+lonarr(dmn_nbr) ; [nbr] Count
; Read required dimensions
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr
ncdf_varget,nc_id,'lon',lon
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,foo,lat_nbr
ncdf_varget,nc_id,'lat',lat
time_id=ncdf_dimid(nc_id,'time')
if time_id gt -1 then begin
; Time dimension is in file, retrieve time size and coordinate
ncdf_diminq,nc_id,time_id,foo,time_nbr
ncdf_varget,nc_id,'time',time
endif; time_id gt -1
; Recall time_idx is passed as 0-based index
if time_idx gt -1 and time_nbr gt 1 then begin
	if dmn_nbr eq 3 then srt=[0,0,time_idx]
	if dmn_nbr eq 4 then srt=[0,lvl_idx,0,time_idx]
endif; endif multiple time samples in this file
if dmn_nbr eq 2 then cnt=[lon_nbr,lat_nbr]
if dmn_nbr eq 3 then cnt=[lon_nbr,lat_nbr,1]
if dmn_nbr eq 4 then cnt=[lon_nbr,1,lat_nbr,1]
if lvl_idx gt -1 then begin
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
ncdf_varget,nc_id,'lev',lev
lvl_sng=auto_sng(round(lev(lvl_idx)),0)+' mb'
ttl=lvl_sng+' '+ttl
;data=reform(data(*,lvl_idx,*)) ; CCM disk data appears to IDL as data(lon,lev,lat,time)
endif; lvl_idx != -1
print,'time_idx = ',time_idx
print,'time_nbr = ',time_nbr
print,'srt = ',srt
print,'cnt = ',cnt
; Get variable data
ncdf_varget,nc_id,fld_nm,data,offset=srt,count=cnt
data=reform(data) ; Get rid of degenerate dimensions, if any
ncdf_close,nc_id
; End of NetCDF commands

; Fix data to be IDL map-compatible with -180 < longitude < 180
; NB: may need to reverse latitudes from ERBE files
if lon(lon_nbr-1) gt 90.0 and lon(0) ge 0.0 then begin
	lon=shift(lon,lon_nbr/2)
	west_lon_idx=where(lon ge 180.0,count)
	if count gt 0 then lon(west_lon_idx)=-360.0+lon(west_lon_idx)
	data=shift(data,lon_nbr/2,0)
endif; endif
good_idx=where(data lt 1.0e20,count)
if scl ne 1.0 then print,'gcm_xy(): INFO Scaling input field by ',scl
if count gt 0 then data(good_idx)=scl*data(good_idx)

; Wrapped regions overridde these defaults 
map_lmt_dgr=[lat(0),lon(0),lat(lat_nbr-1),lon(lon_nbr-1)]
lon_ctr=0.0
lat_ctr=0.0

; Hyperslabber code
if n_elements(rng_x) eq 2 and n_elements(rng_y) eq 2 then begin
lon_min=rng_x(0)
lon_max=rng_x(1)
lat_min=rng_y(0)
lat_max=rng_y(1)
foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)

; Is there to be any zonal regionalization?
if lon_min_idx ne 0 or lon_max_idx ne lon_nbr-1 then begin

print,'Beginning Region Truncation...'
print,'lon_nbr = ',lon_nbr,' min(lon) = ',min(lon),' max(lon) = ',max(lon)
print,'lon_min = ',lon_min,' lon_min_idx = ',lon_min_idx,' = ',lon(lon_min_idx)
print,'lon_max = ',lon_max,' lon_max_idx = ',lon_max_idx,' = ',lon(lon_max_idx)
print,''
print,'lat_nbr = ',lat_nbr,' min(lat) = ',min(lat),' max(lat) = ',max(lat)
print,'lat_min = ',lat_min,' lat_min_idx = ',lat_min_idx,' = ',lat(lat_min_idx)
print,'lat_max = ',lat_max,' lat_max_idx = ',lat_max_idx,' = ',lat(lat_max_idx)

; If zonal range crosses date line.....
if lon_min gt lon_max then begin
; Find current index of desired left-most longitude
	shift_amt=lon_nbr-lon_min_idx
	data=shift(data,shift_amt,0)
	lon=shift(lon,shift_amt)

; Truncate data to right of desired right-most longitude
	data=data(0:lon_max_idx+shift_amt,*)
	lon=lon(0:lon_max_idx+shift_amt)

	lon_ctr=ngl_dgr_2_m180_p180(0.5*(lon_min+lon_max+360.0))
	map_lmt_dgr(1)=ngl_dgr_2_m180_p180(lon_min)
	map_lmt_dgr(3)=ngl_dgr_2_m180_p180(360.0+lon_max)

endif; endif data region crossed date line

; If zonal range does not cross date line.....
if lon_min lt lon_max then begin
; Truncate data outside desired longitude region
	data=data(lon_min_idx:lon_max_idx,*)
	lon=lon(lon_min_idx:lon_max_idx)

	lon_ctr=0.5*(lon_min+lon_max)
	map_lmt_dgr(1)=lon_min
	map_lmt_dgr(3)=lon_max

endif; endif data region did not cross date line

endif; endif longitude needed regionalization

; Is there to be any meridional regionalization?
if ((lat_min_idx ne 0) or (lat_max_idx ne lat_nbr-1)) then begin
; Truncate data outside desired latitude region
	data=data(*,lat_min_idx:lat_max_idx)
	lat=lat(lat_min_idx:lat_max_idx)

	lat_ctr=0.5*(lat_min+lat_max)
	map_lmt_dgr(0)=lat_min
	map_lmt_dgr(2)=lat_max

endif; endif latitude needed truncating

lat_nbr=n_elements(lat)
lon_nbr=n_elements(lon)
data_nbr=n_elements(data)

print,'Post rgn_trunc....'
print,'lat_nbr = ',lat_nbr
print,'lon_nbr = ',lon_nbr
print,'data_nbr = ',data_nbr

endif; endif hyperslab

abc=lon
ord=lat
; Do not wrap zonally truncated regions
if lon_min_idx ne 0 or lon_max_idx ne lon_nbr-1 then wrp=0
; NB: wrap procedure shifts field 1 pixel west 
if wrp then begin
	data=[data,data(0,*)]
	abc=[abc,180.0]
endif

good_idx=where(data lt 1.0e20,count)
if count gt 0 then data_min=min(data(good_idx)) else data_min=min(data)
if count gt 0 then data_max=max(data(good_idx)) else data_max=max(data)
print,'missing data: ',data_nbr-count,' points = ',100.0*(data_nbr-count)/data_nbr,' percent'

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

if auto_scl then begin

; Automatically set and scale the contour levels
cntr_lvl_nbr=11

rsn_prj=(data_max-data_min)/cntr_lvl_nbr
if rsn_prj eq 0.0 then grd_ntv = 1.0 $
else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
else if rsn_prj le 1.0 then grd_ntv = 1.0 $
else if rsn_prj le 5.0 then grd_ntv = 1.0 $
else if rsn_prj le 10.0 then grd_ntv = 5.0 $
else if rsn_prj le 100.0 then grd_ntv = 10.0 $
else grd_ntv = 10.0^round(alog10(rsn_prj))

cntr_lvl_max=(data_max+grd_ntv)-abs(data_max mod grd_ntv)
if (data_min lt 0.0) then cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) else cntr_lvl_min=data_min-abs(data_min mod grd_ntv) 

if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)

end; end if auto_scl

cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

if fld_nm eq 'DSTMPC' then begin
	cntr_ntv=2
	cntr_lvl=[0,25,50,75,100,150,200,250,300,400,500,1000] ; NB: Used in ZBN03
;	cntr_lvl=100*[0,0.1,0.2,.5,1,2,5,10,20,50,100,200,500,1000]
endif;  endif DSTMPC
if fld_nm eq 'DSTMPCCA' then begin
	cntr_ntv=2
	cntr_lvl=0.01*[0,25,50,75,100,150,200,250,300,400,500,1000,2000,5000,10000]
endif;  endif DSTMPCCA
if fld_nm eq 'DSTODXC' then begin
	cntr_ntv=2
	cntr_lvl=[0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.75,1.0,1.5] ; NB: Used in ZBN03
endif;  endif DSTODXC
if fld_nm eq 'DSTSFMBL' or fld_nm eq 'DSTSFDRY' or fld_nm eq 'DSTSFPCP' or fld_nm eq 'DSTSFDPS' then begin
	cntr_ntv=2
	cntr_lvl=0.01*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000] ; NB: Used in ZBN03
endif;  endif DSTSFMBL
;if fld_nm eq 'rdb_fct_gmr' or fld_nm eq 'rdb_fct_hyd' then begin
;	cntr_ntv=2
;	cntr_lvl=0.0001*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000]
;endif;  endif rdb_fct_gmr
if 0 then begin
if fld_nm eq 'DSTSFDRY' or fld_nm eq 'DSTSFPCP' then begin
	cntr_ntv=2
	cntr_lvl=0.01*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000] ; NB: Used in ZBN03
endif;  endif DSTSFMBL
if fld_nm eq 'DSTSFDPS' or fld_nm eq 'DSTSFDPSDST' then begin
	cntr_ntv=2
	cntr_lvl=1.0*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000] ; NB: Used in ZBN03
endif;  endif DSTSFDPS
endif
if fld_nm eq 'DSTSFDPSRATCACO3' then begin
	cntr_ntv=2
	cntr_lvl=0.1*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000]
endif; endif DSTSFDPSRATCACO3
if fld_nm eq 'DSTSFDPSCACO3' or fld_nm eq 'DSTSFDPSCA' or fld_nm eq 'DSTSFDPSC' then begin
	cntr_ntv=2
	cntr_lvl=0.1*[0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000]
endif;  endif DSTSFDPSCACO3
; In case special field overwrote cnt_lvl_nbr
cntr_lvl_nbr=n_elements(cntr_lvl)

cntr_lvl_min=min(cntr_lvl)
cntr_lvl_max=max(cntr_lvl)
cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl,dbg_lvl=0
if unit_pfx_sng ne '' then print,'cntr_lbl_mk: WARNING unit_pfx_sng == ',unit_pfx_sng

; Customize contour labels
if cntr_lvl_nbr lt 15 then cntr_chr_sz=0.83*chr_sz else cntr_chr_sz=0.7*chr_sz ; charsize of level labels
if pll eq 1 then begin
for cntr_lvl_idx=0,cntr_lvl_nbr-2 do begin
	if shd_xpr eq 'lt' then if cntr_lvl(cntr_lvl_idx) lt shd_val then cntr_fll_idx(cntr_lvl_idx)=10 else cntr_fll_idx(cntr_lvl_idx)=clr_wht_idx
	if shd_xpr eq 'gt' then if cntr_lvl(cntr_lvl_idx) gt shd_val then cntr_fll_idx(cntr_lvl_idx)=10 else cntr_fll_idx(cntr_lvl_idx)=clr_wht_idx
	if shd_xpr eq 'ge' then if cntr_lvl(cntr_lvl_idx) ge shd_val then cntr_fll_idx(cntr_lvl_idx)=10 else cntr_fll_idx(cntr_lvl_idx)=clr_wht_idx
endfor
for cntr_lvl_idx=0,cntr_lvl_nbr-1 do begin
	if dsh then if cntr_lvl(cntr_lvl_idx) lt 0.0 then cntr_ln_sty(cntr_lvl_idx)=2 
	if not dsh then cntr_ln_sty(cntr_lvl_idx)=0
	if cntr_lvl(cntr_lvl_idx) eq 0.0 then cntr_which_lbl(cntr_lvl_idx)=1
endfor
endif; endif pll eq 1

print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
print,'minimum data = ',data_min
print,'minimum contour = ',cntr_lvl(0)

; Setup colorbar---this sets cbar_clr_nbr and cbar_lgn
if img and cbar then begin
	cbar_clr_nbr=(!d.table_size-3)
	cbar_lgn=strarr(cbar_clr_nbr+1)
	cbar_lgn(0)=cntr_lvl_lbl(0)
	for cntr_lvl_idx=1,cntr_lvl_nbr-2 do begin
		color_idx=(!d.table_size-1)*(cntr_lvl(cntr_lvl_idx)-cntr_lvl_min)/(cntr_lvl_max-cntr_lvl_min)
		color_idx=round(color_idx)
		cbar_lgn(color_idx)=cntr_lvl_lbl(cntr_lvl_idx)
	endfor
	cbar_lgn(cbar_clr_nbr)=cntr_lvl_lbl(cntr_lvl_nbr-1)
endif else begin; endif img
	cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
	cbar_lgn=cntr_lvl_lbl
endelse; endelse not img
; fxm: Perhaps this should be in cntr_lbl_mk()?
; Remove every other lable when there are lots of labels
if cntr_lvl_nbr gt 15 then begin
    for cntr_lvl_idx=1,cntr_lvl_nbr-2,2 do begin
        cbar_lgn(cntr_lvl_idx)=''
    endfor
endif; cntr_lvl_nbr gt 15 

cbar_idx=indgen(cbar_clr_nbr)+2
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.83*chr_sz
cbar_unit=''
cbar_lbl_sz=0.83*chr_sz
cbar_lbl_sz=chr_sz

if pll eq 1 then clr_mk,10,clr_ord,pll,0 else clr_mk,cbar_clr_nbr,clr_ord,pll,0

if img then clr_rfr,clr_tbl_nbr,fl_clr_tbl,usr_dfn_clr_tbl
if img then m_thk=2 else m_thk=1;

; Do not plot ERBE data poleward of 60 degrees.
;if not img and strpos(fl_nm,'erbe') ne -1 then begin
;good_ord=where(abs(ord) lt 60.0,nbr_ord)
;ord=ord(good_ord)
;nbr_abc=n_elements(abc)
;tmp_data=fltarr(nbr_abc,nbr_ord)
;for idx_ord=0,nbr_ord-1 do begin
;	tmp_data(*,idx_ord)=data(*,good_ord(idx_ord))
;endfor
;data=tmp_data
;endif

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
if cbar then mrg_btm=mrg_btm+4
mrg_rgt=1.8 ; 3 is default
if lft then mrg_lft=6.2 else mrg_lft=mrg_rgt ; 10 is default
if prn then begin
	if azi then mrg_lft=0.7
	if azi then mrg_rgt=0.7
	if azi then mrg_top=0.7 else mrg_top=0.7 ; 2 is default
	if azi then mrg_btm=0.7 else mrg_btm=0.7 ; 4 is default
	if top then mrg_top=mrg_top+2.0
	if btm then mrg_btm=mrg_btm+2.0
	if cbar then mrg_btm=mrg_btm+3.5
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

if azi then begin
	ord_min=30
	ord_max=90

	if(ord_max gt abs(ord_min)) then lat_ctr=ord_max else lat_ctr=ord_min
	if(ord_max gt abs(ord_min)) then lon_lbl_lat=ord_min-2.5 else lon_lbl_lat=ord_max+2.5
	lon_ctr=0.5*(abc_min+abc_max)
	if prn then cntr_chr_sz=1.5*cntr_chr_sz

; Set map for Azimuthal projection
map_lmt_dgr=[ord_min,abc_min,ord_max,abc_max]
map_set, $
	0.0, $			; [dgr] Latitude at center of map
	lon_ctr, $		; [dgr] Longitude at center of map
	0.0, $			; [dgr] North rotation angle
	limit=map_lmt_dgr, $	; [latmin,lonmin,latmax,lonmax] drawing box
	/azimuthal, $	;/aitoff,/azimuthal,/conic,/cylindrical,/gnomic,/lambert,/mercator,/mollweide,/orthographic,/satellite
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	/noborder 		; Do not draw a border around the grid

	data_scl=bytscl(data,min=cntr_lvl_min,max=cntr_lvl_max,top=(!d.table_size-3))+2
	data_scl_map=map_image(data_scl,x_srt,y_srt,x_sz,y_sz, $
;			latmin=ord_min,latmax=ord_max, $
			lonmin=abc_min,lonmax=abc_max, $
			/bilin)

; Fill image for Azimuthal projection
	if img then if(!d.name eq 'PS') then tv,data_scl_map,x_srt,y_srt,xsize=x_sz,ysize=y_sz else tv,data_scl_map,x_srt,y_srt

; Draw continents and axes for Azimuthal projection
map_set, $
	0.0, $			; [dgr] Latitude at center of map
	lon_ctr, $		; [dgr] Longitude at center of map
	0.0, $			; Rotation
	latdel=15.0, $		; Latitude grid spacing
	londel=90.0, $		; Longitude grid spacing
	limit=map_lmt_dgr, $	; [latmin,lonmin,latmax,lonmax] drawing box
	/continent, $		; Plot continental boudaries
	con_color=clr_blk_idx, $ ; Color index of continental boundaries
	color=clr_blk_idx, $ ; Color index of map borders
	/azimuthal, $	$	;/aitoff,/azimuthal,/conic,/cylindrical,/gnomic,/lambert,/mercator,/mollweide,/orthographic,/satellite
	/label, $		; Label parallels and meridians
	latalign=0.5, $		; 0.5 centers the labels above/below the text
	lonalign=0.5, $		; 0.5 centers the labels next to the text
	lonlab=lon_lbl_lat, $	; Latitudes at which to place longitude labels
	latlab=0, $		; Longitudes at which to place latitude labels
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	charsize=chr_sz, $
	mlinethick=m_thk, $	; Default continental boundary line thickness is 2.0
	/noborder, $		; Do not draw border around grid
	/noerase, $		; Do not erase current image before drawing
	/grid			; Draw grid of parallels and meridians

; These points supply left and right boundaries of plot on page
	map_psn_nrm=convert_coord([-90,90],[ord_min,ord_min],/data,/to_normal) ; [lon_min,lon_max],[lat_min,lat_max]
	plt_rgn_nrm=[map_psn_nrm(0),map_psn_nrm(1),map_psn_nrm(3),map_psn_nrm(4)] ; [lon_min,lat_min,lon_max,lat_max]
	cbar_psn=[plt_rgn_nrm(0),0.06,plt_rgn_nrm(2),0.10]; [x_min,y_min,x_max,y_max]
; These points supply the bottom and top boundaries of the plot
	map_psn_nrm=convert_coord([0,180],[ord_min,ord_min],/data,/to_normal) ; [lon_min,lon_max],[lat_min,lat_max]
	plt_rgn_nrm=[map_psn_nrm(0),map_psn_nrm(1),map_psn_nrm(3),map_psn_nrm(4)] ; [lon_min,lat_min,lon_max,lat_max]
	ttl_sz=chr_sz*1.2
	xyouts,0.5*(plt_rgn_nrm(0)+plt_rgn_nrm(2)),plt_rgn_nrm(3)+0.5*ttl_sz*(float(!d.y_ch_size)/!d.y_vsize),ttl,size=ttl_sz,alignment=0.5,/NORMAL

endif; endif azi

if not azi then begin

; Do not set this here because we may need to add 360.0 to lon_max
;map_lmt_dgr=[lat_min,lon_min,lat_max,lon_max] ; wrapped regions will overridde this
; Set map coordinates for Cylindrical plots
map_set, $
	0.0, $			; [dgr] Latitude at center of map
	lon_ctr, $		; [dgr] Longitude at center of map
	0.0, $			; [dgr] North rotation angle
	limit=map_lmt_dgr, $	; [latmin,lonmin,latmax,lonmax] drawing box
	/cylindrical, $	;/aitoff,/azimuthal,/conic,/cylindrical,/gnomic,/lambert,/mercator,/mollweide,/orthographic,/satellite
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	/noborder 		; Do not draw a border around the grid

; Convert position of map edges from data to normal coordinates
map_psn_nrm=convert_coord([map_lmt_dgr(1),map_lmt_dgr(3)],[map_lmt_dgr(0),map_lmt_dgr(2)],/data,/to_normal) ; [lon_min,lon_max],[lat_min,lat_max]
; Convert returned information to an x-y vector, getting rid of z data
plt_rgn_nrm=[map_psn_nrm(0),map_psn_nrm(1),map_psn_nrm(3),map_psn_nrm(4)] ; [lon_min,lat_min,lon_max,lat_max]

if img then begin

arr_sz=size(data)
x_arr_sz=arr_sz(1)
y_arr_sz=arr_sz(2)

; Convert x-y normal coordinate vector to device coords
map_psn_dvc=fix(plt_rgn_nrm*[!d.x_size,!d.y_size,!d.x_size,!d.y_size])
dvc_wdt=map_psn_dvc(2)-map_psn_dvc(0)-1
dvc_hgt=map_psn_dvc(3)-map_psn_dvc(1)-1
if(!d.name eq 'PS') then begin

; Postscript automatically scales image up to xsize X ysize
	tv,bytscl(data,min=cntr_lvl_min,max=cntr_lvl_max,top=(!d.table_size-3))+2, $
	        map_psn_dvc(0)+1,map_psn_dvc(1)+1,xsize=dvc_wdt,ysize=dvc_hgt,/device
endif else begin; not PS
; Other devices (pixel based systems) need scaling done for them in order to fit device

	tv,bytscl(congrid(data,dvc_wdt,dvc_hgt),min=cntr_lvl_min,max=cntr_lvl_max,top=(!d.table_size-3))+2, $
	        map_psn_dvc(0)+1,map_psn_dvc(1)+1,/device
endelse; not PS

endif; endif img 

cbar_psn=[plt_rgn_nrm(0),0.06,plt_rgn_nrm(2),0.10]; [x_min,y_min,x_max,y_max]

endif; endif not azi

; Draw colorbar legend for all projections
if cbar then clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz

; Draw filled contours for Cylindrical maps
print,'gcm_xy_batch: ttl = ',ttl
if not img and not azi then contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $	
	xrange=[abc_min,abc_max], $ ; Longitude increases on x-axis
	yrange=[ord_min,ord_max], $ ; Latitude increases on y-axis
	xstyle=13, $
	ystyle=13, $
	charsize=chr_sz, $
	position=plt_rgn_nrm, $
;	/closed, $			
	/fill, $			
	/noerase

; Draw contour lines for all projections
if cnt_ovl then contour, $
	data, $
	abc, $
	ord, $
	max_value=1.0e20, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $
	c_annotation=cntr_lvl_lbl, $
	c_charsize=cntr_chr_sz, $
;	/closed, $			
	/follow, $
	/overplot

if not azi then begin

ttl_sz=chr_sz*1.2
if img then xyouts,0.5*(plt_rgn_nrm(0)+plt_rgn_nrm(2)),plt_rgn_nrm(3)+0.5*ttl_sz*(float(!d.y_ch_size)/!d.y_vsize),ttl,size=ttl_sz,alignment=0.5,/NORMAL

; Set boundaries for Cylindrical projections
; NB: Until IDL 4.0, this commands was not needed because map_set() set x.crange and y.crange, which the axes routines need to do correct labeling
; Without the following two lines, x.crange and y.crange will be [0.0,1.0] and axis labels will be wrong
!x.range=[abc_min,abc_max] ; Longitude increases on x-axis
!y.range=[ord_min,ord_max] ; Latitude increases on y-axis

; Draw continents for Cylindrical projections
map_set, $
	0.0, $			; [dgr] Latitude at center of map
	lon_ctr, $		; [dgr] Longitude at center of map
	0.0, $			; [dgr] North rotation angle
	limit=map_lmt_dgr, $	; [latmin,lonmin,latmax,lonmax] drawing box
	/continent, $		; Plot continental boudaries
	con_color=clr_blk_idx, $ ; Color index of continental boundaries
	color=clr_blk_idx, $ ; Color index of map borders
	/cylindrical, $	;/aitoff,/azimuthal,/conic,/cylindrical,/gnomic,/lambert,/mercator,/mollweide,/orthographic,/satellite
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	mlinethick=m_thk, $	; Default continental boundary line thickness is 2.0
;	position=plt_rgn_nrm, $
	/noborder, $		; Do not draw a border around the grid
;	/grid			; Draw grid of parallels and meridians
	/noerase		; Do not erase current image before drawing

; Draw axes for Cylindrical projections
if btm then ntt=1 else ntt=0
lon_axz_drw,lon_top=abc_max,lon_btm=abc_min,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,axz_vrt=0
if lft then ntt=1 else ntt=0
lat_axz_drw,lat_top=ord_max,lat_btm=ord_min,ntt=ntt,ttl=y_ttl,chr_sz=chr_sz,axz_vrt=1

endif; endif not azi

;if pnl_lbl ne '' then xyouts,0.1,0.9,pnl_lbl,size=chr_sz,/NORMAL
if pnl_lbl ne '' then xyouts,abc_min+0.03*(abc_max-abc_min),ord_max-0.08*(ord_max-ord_min),pnl_lbl,size=chr_sz,/DATA

if not prn and dvc ne 'mpg' then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif else clr_rfr,22,"",0

; Jump to here if input file does not exist
end_of_procedure: foo=1

end; end gcm_xy()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM XY commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM VY
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_vy_bch, $
	fld_top=fld_top, $ ; First field is top plot
	lvl=lvl, $
	prn=prn
if n_elements(prn) eq 0 then prn=0
if n_elements(fld_top) eq 0 then fld_top=0
if n_elements(lvl) eq 0 then lvl=-1

fld3_stf=[ $
;	['NPCO2O2','!5O!I2!N!9. !5O!I2!N Complex Loading','!8N!IO2-O2!N','!9X!5 10!E52!N mlc!E2!N m!E-5!N'], $
;	['FSATFRC','!5Zonal Average Atmospheric Absorption','!5F!IABS!N','!5W m!E-2!N'] ]
;	['Q','!5Vapor Mixing Ratio','!8q!5!Iv!N','!5g kg!E-1!N'], $
;	['OMEGA','!5Vertical Velocity','!7x!5','!5mb day!E-1!N'] ]
;	['TOTLWP','!5Liquid Water Path','!5LWP','!5g m!E-2!N'] ]
;	['TOTIWP','!5Ice Water Path','!5IWP','!5g m!E-2!N'], $
;	['TOTCWP','!5Condensed Water Path','!5CWP','!5g m!E-2!N'] ]
;	['FSNSC','!5Net Clear Sky Solar Radiation SFC','!5SFC!Ic!N','!5W m!E-2!N'], $
;	['ALBEDO','!5Albedo','!8A!5','!5%'] ]
;	['FSNS','!5Net Solar Radiation SFC','!5SFC','!5W m!E-2!N'] ]
;	['FSNT','!5Net Solar Radiation TOA','!5TOA','!5W m!E-2!N'], $
;	['FLNT','!5Outgoing Longwave Radiation','!5OLR!5','!5W m!E-2!N'] ]
;	['FSNTC','!5Net Clear Sky Solar Radiation TOA','!5TOA!Ic!N','!5W m!E-2!N'], $
;	['FLNTC','!5Clear Sky Outgoing Longwave Radiation','!5OLR!Ic!N','!5W m!E-2!N'] ]
;	['PRECT','!5Total Precipitation','!8P!T!5!N','!5mm day!E-1!N'] ]
;	['LWCF','!5Longwave Cloud Forcing','!5LWCF!5','!5W m!E-2!N'] ]
	['SWCF','!5Shortwave Cloud Forcing','!5SWCF!5','!5W m!E-2!N'] ]
;	['CLDTOT','!5Total Cloud Amount','!8A!5','!5%'], $
;	['CLDHGH','!5High Cloud Amount','!8A!5','!5%'] ]
fld_nbr3=n_elements(fld3_stf(0,*))
foo={fld3_sct,idx:0,nm:'',sng:'',abb:'',unit:''}
fld3=replicate({fld3_sct},fld_nbr3)
for idx=0,fld_nbr3-1 do begin
	fld3(idx).idx=idx
	fld3(idx).nm=fld3_stf(0,idx)
	fld3(idx).sng=fld3_stf(1,idx)
	fld3(idx).abb=fld3_stf(2,idx)
	fld3(idx).unit=fld3_stf(3,idx)
endfor; end loop over fld3

xpt_stf=[ $
;	['dmr06','!5O!I2!N!9. !5O!I2!N']]
;	['dmr05','!5O!I2!N!9. !5N!I2!N'], $
;	['dmr04','!5O!I2!N!9. !5O!I2!N + !5O!I2!N!9. !5N!I2!N']]
;	['dmr03','!5(H!I2!NO)!I2!N']]
;	['isccp','!5ISCCP'], $
;	['ecmwf','!5ECMWF'], $
;	['ssmi','!5SSMI'], $
;	['gpcp','!5GPCP'], $
	['erbe_b','!5ERBE'], $
;	['422','!5CCM2'], $
;	['nflux18','!5NFLUX18'], $
;	['amip5','!5CCM'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
	['sld012d','!5CCM3']]
;	['spcp_85','!5ANV']]
;	['erbe_b','!5ERBE']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpt

;mth_sng=['1202','0608']
mth_sng=['01','07']
;mth_sng=['01']
;mth_sng=['07']
;mth_sng=['8589']
mth_nbr=n_elements(mth_sng)
for mth_idx=0,mth_nbr-1 do begin
for fld_idx3=0,fld_nbr3-1 do begin
	ttl=fld3(fld_idx3).sng
	y_ttl=''
	x_ttl=fld3(fld_idx3).abb+' ('+fld3(fld_idx3).unit+')'
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	info=1
	rng_x=0
	prs_sng=''
	fl_nm_out_xpt_sng=''
	if fld3(fld_idx3).nm eq 'NPCO2O2' then scl=1.0e-16
 	if fld3(fld_idx3).nm eq 'TOTLWP' then rng_x=[0,140]
	if fld3(fld_idx3).nm eq 'TOTIWP' then rng_x=[0,140]
	if fld3(fld_idx3).nm eq 'TOTCWP' then rng_x=[0,140]
	if fld3(fld_idx3).nm eq 'FLNT' then rng_x=[100,290]
	if fld3(fld_idx3).nm eq 'FLNTC' then rng_x=[100,300]
	if fld3(fld_idx3).nm eq 'FSNT' then rng_x=[100,400]
	if fld3(fld_idx3).nm eq 'FSNTC' then rng_x=[100,450]
	if fld3(fld_idx3).nm eq 'ALBEDO' then scl=100
	if fld3(fld_idx3).nm eq 'ALBEDO' then rng_x=[0,100]
;	if fld3(fld_idx3).nm eq 'FSNS' then rng_x=[100,400]
;	if fld3(fld_idx3).nm eq 'FSNSC' then rng_x=[100,450]
	if fld3(fld_idx3).nm eq 'PRECT' then rng_x=[0,10]
	if fld3(fld_idx3).nm eq 'LWCF' then rng_x=[0,60]
	if fld3(fld_idx3).nm eq 'SWCF' then rng_x=[-160,0]
	if fld3(fld_idx3).nm eq 'FSATFRC' then rng_x=[0,2.5]
	if fld3(fld_idx3).nm eq 'CLDTOT' then rng_x=[0,100]
	if strpos(fld3(fld_idx3).nm,'CLD') ne -1 then scl=100
	if fld3(fld_idx3).nm eq 'CLDHGH' then rng_x=[0,50]
	if fld3(fld_idx3).nm eq 'Q' then scl=1000.0
	if fld3(fld_idx3).nm eq 'Q' then prs_sng='prs_'
	if fld3(fld_idx3).nm eq 'OMEGA' then prs_sng='prs_'
	if fld3(fld_idx3).nm eq 'OMEGA' then scl=-864
	fl_nm_in=strarr(xpt_nbr)
	lbl_sng=strarr(xpt_nbr)
	for xpt_idx=0,xpt_nbr-1 do begin
		yr_sng='8589_'
;		yr_sng=''
		if fld3(fld_idx3).nm eq 'PRECT' and xpt(xpt_idx).nm ne 'gpcp' then scl=8.64e7
		if strpos(xpt(xpt_idx).nm,'ecmwf') ne -1 then begin
			if mth_sng(mth_idx) eq '01' then yr_sng='9095_' else yr_sng='8994_'
			prs_sng='prs_'
		endif; endif ecmwf
		if strpos(xpt(xpt_idx).nm,'gpcp') ne -1 then begin
			if mth_sng(mth_idx) eq '01' then yr_sng='8894_' 
			if mth_sng(mth_idx) eq '07' then yr_sng='8793_' 
		endif
		if strpos(xpt(xpt_idx).nm,'ssmi') ne -1 then begin
			if mth_sng(mth_idx) eq '01' then yr_sng='8795_' 
			if mth_sng(mth_idx) eq '07' then yr_sng='8791_' 
		endif
;		fl_nm_in(xpt_idx)=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+'xavg_'+yr_sng+mth_sng(mth_idx)+'.nc'
		fl_nm_in(xpt_idx)=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+yr_sng+mth_sng(mth_idx)+'_x.nc'
		lbl_sng(xpt_idx)=xpt(xpt_idx).sng
		fl_nm_out_xpt_sng=fl_nm_out_xpt_sng+xpt(xpt_idx).nm+'_'
	endfor; end loop over experiment
	fl_nm_out=getenv('DATA')+'/ps/'+fl_nm_out_xpt_sng+yr_sng+mth_sng(mth_idx)+'_x_'+fld3(fld_idx3).nm+'.eps'
	if prn then begin
		if fld_top then begin
;		if fld_idx3 eq 0 then top=1 else top=0
		if fld_idx3 eq fld_nbr3-1 then btm=1 else btm=0
		x_ttl='!5(g m!E-2!N)'
		if mth_sng(mth_idx) eq '01' then ttl='Jan '+ttl else if mth_sng(mth_idx) eq '07' then ttl='Jul '+ttl
		endif; endif fld_top
		if not fld_top then begin
		if mth_idx eq 0 then info=1 else info=0
		if mth_idx eq 0 then top=1 else top=0
		if mth_idx eq mth_nbr-1 then btm=1 else btm=0
		endif; endif not fld_top
	endif; endif prn
	if prn then begin
		x_sz=6.5
		y_sz=3.2
		if top then y_sz=y_sz*1.11
		if btm then y_sz=y_sz*1.22
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
	endif; endif prn
	print,'btm = ',btm
	gcm_vy,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,fld_nm=fld3(fld_idx3).nm,rng_x=rng_x,rng_y=[-90,90],chr_sz=chr_sz,top=top,btm=btm,scl=scl,info=info,lvl=lvl,lbl_sng=lbl_sng
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over fld3
endfor; end loop over mth

end; end gcm_vy_bch()
pro gcm_vy, $
	fl_nm=fl_nm, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	info=info, $
	pll=pll, $
	clr_ord=clr_ord, $
	clr_tbl_nbr=clr_tbl_nbr, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	lbl_sng=lbl_sng, $
	lvl=lvl, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	scl=scl, $
	fld_nm=fld_nm

@ibp_clr.com

; gcm_vy,ttl='!5ERBE SWCF (W m!E-2!N)'
;gcm_vy,x_ttl='!7x!5 (mb day!E-1!N)',y_ttl='Latitude',mth_srt=0,mth_end=12,fl_nm=[getenv('DATA')+'/ecmwf/ecmwf_prs_xavg_9095_01.nc',getenv('DATA')+'/amip5/amip5_prs_xavg_8589_01.nc',getenv('DATA')+'/spcp_85/spcp_85_prs_xavg_8589_01.nc'],ttl='Prssure Velocity',prn=0,fld_nm='OMEGA',rng_x=[-40,60],rng_y=[-90,90],chr_sz=chr_sz,top=top,btm=btm,scl=-864,info=info,lvl=6

if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(lvl) eq 0 then lvl=-1
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ERBE SWCF (W m!E-2!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(info) eq 0 then info=1
if n_elements(pll) eq 0 then pll=4
if n_elements(clr_ord) eq 0 then clr_ord=1
if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=22
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl=''
if n_elements(x_ttl) eq 0 then x_ttl=''
if n_elements(lbl_sng) eq 0 then lbl_sng=''
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/erbe_b/erbe_b_xavg_8589_01.nc'
if n_elements(fld_nm) eq 0 then fld_nm='SWCF'

clr_rfr,clr_tbl_nbr,"",0
erase
!p.multi=0

fl_nbr=n_elements(fl_nm)
if fl_nbr eq 1 then fl_nm=[fl_nm]

for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_nm(fl_idx))
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,dim_foo,lat_nbr
if fl_idx eq 0 then data=fltarr(lat_nbr,fl_nbr)
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,fld_nm,data_foo
data_foo=reform(data_foo)
if lvl gt -1 then begin
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,dim_foo,lev_nbr
ncdf_varget,nc_id,'lev',lev
lvl_sng=auto_sng(round(lev(lvl)),0)+' mb'
data_foo=reform(data_foo)
data_foo=reform(data_foo(lvl,*))
endif; lvl != -1
data(*,fl_idx)=data_foo
ncdf_close,nc_id
; End of NetCDF commands

endfor; end loop over files

good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=data
ord=lat
if lvl gt -1 then ttl=lvl_sng+' '+ttl

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc(where(abc lt 1.0e20)))
	abc_max=max(abc(where(abc lt 1.0e20)))
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=4.5 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.5
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

for fl_idx=0,fl_nbr-1 do begin
mss_val=where(abc(*,fl_idx) gt 1.0e20,cnt)
if cnt gt 0 then ord(mss_val)=abc(mss_val,fl_idx)
if fl_idx eq 0 then begin
plot, $
	abc(*,fl_idx), $
	ord, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=1, $
	ystyle=13, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	linestyle=0

endif else begin; fl_idx != 0
oplot, $
	abc(*,fl_idx), $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx

endelse; fl_idx != 0
endfor; end loop over files

lat_axz_drw,lat_top=!y.crange(1),lat_btm=!y.crange(0),ntt=1,ttl=y_ttl,chr_sz=chr_sz,axz_vrt=1

if info then begin

if (fld_nm eq 'SWCF') or (fld_nm eq 'FLNT') or (fld_nm eq 'FLNTC') then ln_lgn_x1=.25 else ln_lgn_x1=0.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for fl_idx=0,fl_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=fl_idx,thick=2.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),size=chr_sz,/NORMAL
endfor; end loop over fl

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_vy()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM VY commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM XV
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accesse` as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_xv_bch, $
	lvl=lvl, $
	ocn=ocn, $
	prn=prn
if n_elements(prn) eq 0 then prn=0
if n_elements(lvl) eq 0 then lvl=-1
if n_elements(ocn) eq 0 then ocn=0
if ocn then ocn_sng='ocn_' else ocn_sng=''

fld3_stf=[ $
;	['Q','!5Vapor Mixing Ratio','!8q!5!Iv!N','!5g kg!E-1!N'], $
;	['OMEGA','!5Vertical Velocity','!7x!5','!5mb day!E-1!N'] ]
	['LWCF','!5Longwave Cloud Forcing','!5LWCF!5','!5W m!E-2!N'], $
	['SWCF','!5Shortwave Cloud Forcing','!5SWCF!5','!5W m!E-2!N'] ]
;	['FSNS','!5Net Surface Shortwave','SW','!5W m!E-2!N'], $
;	['FLNS','!5Net Surface Longwave','LW','!5W m!E-2!N'], $
;	['SHFLX','!5Sensible Heat Flux','SH','!5W m!E-2!N'], $
;	['LHFLX','!5Latent Heat Flux','LH','!5W m!E-2!N'], $
;	['E_P','!5Evaporation - Precipitation','E-P','!5mm day!E-1!N'], $
;	['NET','!5Net Surface Energy','SEB','!5W m!E-2!N'] ]
;	['CLDTOT','!5Total Cloud Amount','!8A!5','!5%'], $
;	['CLDHGH','!5High Cloud Amount','!8A!5','!5%'] ]
fld_nbr3=n_elements(fld3_stf(0,*))
foo={fld3_sct,idx:0,nm:'',sng:'',abb:'',unit:''}
fld3=replicate({fld3_sct},fld_nbr3)
for idx=0,fld_nbr3-1 do begin
	fld3(idx).idx=idx
	fld3(idx).nm=fld3_stf(0,idx)
	fld3(idx).sng=fld3_stf(1,idx)
	fld3(idx).abb=fld3_stf(2,idx)
	fld3(idx).unit=fld3_stf(3,idx)
endfor; end loop over fld3s

xpt_stf=[ $
;	['isccp','!5ISCCP'], $
;	['ecmwf','!5ECMWF'], $
	['erbe_b','!5ERBE'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
	['amip5','!5CCM'], $
;	['sld012d','!5CCM3'], $
	['spcp_85','!5ANV']]
;	['erbe_b','!5ERBE']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

mth_sng=['01','07']
;mth_sng=['01']
mth_nbr=n_elements(mth_sng)
for mth_idx=0,mth_nbr-1 do begin
for fld_idx3=0,fld_nbr3-1 do begin
	ttl=fld3(fld_idx3).sng
	x_ttl=''
	y_ttl=fld3(fld_idx3).abb+' ('+fld3(fld_idx3).unit+')'
	scl=1.0
	chr_sz=1.5
	top=1
	btm=1
	info=1
	rng_y=0
	prs_sng=''
	yavg_sng='yavg_00N20N_'
	if fld_idx3 eq 0 then info=1 else info=0
;	if fld3(fld_idx3).nm eq 'LWCF' then rng_y=[0,60]
;	if fld3(fld_idx3).nm eq 'SWCF' then rng_y=[0,160]
	if fld3(fld_idx3).nm eq 'CLDTOT' then rng_y=[0,100]
	if strpos(fld3(fld_idx3).nm,'CLD') ne -1 then scl=100
	if fld3(fld_idx3).nm eq 'CLDHGH' then rng_y=[0,50]
	if fld3(fld_idx3).nm eq 'Q' then scl=1000.0
	if fld3(fld_idx3).nm eq 'Q' then prs_sng='prs_'
	if fld3(fld_idx3).nm eq 'OMEGA' then prs_sng='prs_'
	if fld3(fld_idx3).nm eq 'OMEGA' then scl=-864
; define everything positive into surface for SEB studies
	if fld3(fld_idx3).nm eq 'E_P' then scl=8.64e7
	if fld3(fld_idx3).nm eq 'FLNS' then scl=-1.0
	if fld3(fld_idx3).nm eq 'SHFLX' then scl=-1.0
	if fld3(fld_idx3).nm eq 'LHFLX' then scl=-1.0
	fl_nm_in=strarr(xpt_nbr)
	lbl_sng=strarr(xpt_nbr)
	for xpt_idx=0,xpt_nbr-1 do begin
		yr_sng='8589_'
		if strpos(xpt(xpt_idx).nm,'ecmwf') ne -1 then begin
			if mth_sng(mth_idx) eq '01' then yr_sng='9095_' else yr_sng='8994_'
			prs_sng='prs_'
		endif; endif ecmwf
		fl_nm_in(xpt_idx)=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+prs_sng+ocn_sng+yavg_sng+yr_sng+mth_sng(mth_idx)+'.nc'
		lbl_sng(xpt_idx)=xpt(xpt_idx).sng
	endfor; end loop over experiment
	if xpt_nbr eq 1 then fl_nm_out=getenv('DATA')+'/ps/'+xpt(0).nm+'_'+ocn_sng+yavg_sng+yr_sng+mth_sng(mth_idx)+'_'+fld3(fld_idx3).nm+'.eps' else fl_nm_out=getenv('DATA')+'/ps/'+ocn_sng+yavg_sng+yr_sng+mth_sng(mth_idx)+'_'+fld3(fld_idx3).nm+'.eps'
	if prn then begin
		if fld_idx3 eq 0 and (fld3(fld_idx3).nm eq 'FSNS' or fld3(fld_idx3).nm eq 'LWCF') then top=1 else top=0
		if fld3(fld_idx3).nm eq 'FSNS' then ttl='!5Surface Energy Balance'
		if fld3(fld_idx3).nm eq 'LWCF' then ttl='!5Cloud Forcing'
		if mth_sng(mth_idx) eq '01' then ttl='!5Jan '+ttl else ttl='!5Jul '+ttl
		if fld_idx3 eq fld_nbr3-1 then btm=1 else btm=0
	endif; endif prn and diff
	if prn then begin
		x_sz=6.5
		y_sz=2
		chr_sz=1.5
		if top then y_sz=y_sz*1.1
		if btm then y_sz=y_sz*1.1
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
	endif; endif prn
	gcm_xv,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,fld_nm=fld3(fld_idx3).nm,rng_x=[0,360],rng_y=rng_y,chr_sz=chr_sz,top=top,btm=btm,scl=scl,info=info,lvl=lvl,lbl_sng=lbl_sng
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over fld3
endfor; end loop over mth

end; end gcm_xv_bch()
pro gcm_xv, $
	fl_nm=fl_nm, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	info=info, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	lbl_sng=lbl_sng, $
	lvl=lvl, $
	mth_srt=mth_srt, $
	mth_end=mth_end, $
	scl=scl, $
	fld_nm=fld_nm

; gcm_xv,ttl='!5ERBE SWCF (W m!E-2!N)'
;gcm_xv,x_ttl='!7x!5 (mb day!E-1!N)',y_ttl='Longitude',mth_srt=0,mth_end=12,fl_nm=[getenv('DATA')+'/ecmwf/ecmwf_prs_yavg_00N20N_9095_01.nc',getenv('DATA')+'/amip5/amip5_prs_yavg_00N20N_8589_01.nc',getenv('DATA')+'/spcp_85/spcp_85_prs_yavg_00N20N_8589_01.nc'],ttl='Prssure Velocity',prn=0,fld_nm='OMEGA',rng_x=[0,360],rng_y=[-90,90],chr_sz=chr_sz,top=top,btm=btm,scl=-864,info=info,lvl=6

if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(lvl) eq 0 then lvl=-1
if n_elements(scl) eq 0 then scl=1.0
if n_elements(ttl) eq 0 then ttl='ERBE SWCF (W m!E-2!N)'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(info) eq 0 then info=1
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl=''
if n_elements(x_ttl) eq 0 then x_ttl=''
if n_elements(lbl_sng) eq 0 then lbl_sng=''
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/erbe_b/erbe_b_yavg_00N20N_8589_01.nc'
if n_elements(fld_nm) eq 0 then fld_nm='SWCF'

erase
!p.multi=0

fl_nbr=n_elements(fl_nm)
if fl_nbr eq 1 then fl_nm=[fl_nm]

for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_nm(fl_idx))
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,dim_foo,lon_nbr
if fl_idx eq 0 then data=fltarr(lon_nbr,fl_nbr)
ncdf_varget,nc_id,'lon',lon
ncdf_varget,nc_id,fld_nm,data_foo
data_foo=reform(data_foo)
if lvl gt -1 then begin
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,dim_foo,lev_nbr
ncdf_varget,nc_id,'lev',lev
lvl_sng=auto_sng(round(lev(lvl)),0)+' mb'
data_foo=reform(data_foo)
data_foo=reform(data_foo(*,lvl))
endif; lvl != -1
data(*,fl_idx)=data_foo
ncdf_close,nc_id
; End of NetCDF commands

endfor; end loop over files

good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=lon
ord=data
if lvl gt -1 then ttl=lvl_sng+' '+ttl

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord(where(ord lt 1.0e20)))
	ord_max=max(ord(where(ord lt 1.0e20)))
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+1.5
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

for fl_idx=0,fl_nbr-1 do begin
if fl_nbr eq 2 and fl_idx eq 1 then ln_sty=fl_idx+1 else ln_sty=fl_idx
mss_val=where(ord(*,fl_idx) gt 1.0e20,cnt)
if cnt gt 0 then abc(mss_val)=ord(mss_val,fl_idx)
if fl_idx eq 0 then begin
plot, $
	abc, $
	ord(*,fl_idx), $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=13, $
	ystyle=0, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	linestyle=ln_sty

endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord(*,fl_idx), $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=ln_sty

endelse; fl_idx != 0
endfor; end loop over files

if !y.crange(0) lt 0 and !y.crange(1) gt 0 then plots,[!x.crange(0),!x.crange(1)],[0,0],linestyle=0,thick=1.0,/DATA

lon_axz_drw,lon_top=!x.crange(1),lon_btm=!x.crange(0),ntt=1,ttl=x_ttl,chr_sz=chr_sz,axz_vrt=0

if info then begin

ln_lgn_x1=0.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for fl_idx=0,fl_nbr-1 do begin
if fl_nbr eq 2 and fl_idx eq 1 then ln_sty=fl_idx+1 else ln_sty=fl_idx
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty,thick=2.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),size=chr_sz,/NORMAL
endfor; end loop over fl

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_xv()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM XV commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM VZ
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[0..lev_nbr-1][0..nbr_sz-1] in C
; is accessed as         foo(0..nbr_sz-1,0..lev_nbr-1)  in IDL
; is accessed as         foo(1..nbr_sz,1..lev_nbr)      in Fortran
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_vz_bch, $
	mth=mth, $
	info=info, $
	diff=diff, $ ; plot differences betseen model fields
	prs=prs, $ ; insert 'prs_' in filename (use files on pressure levels)
	prn=prn
if n_elements(diff) eq 0 then diff=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(mth) eq 0 then mth=0
mth_nbr=n_elements(mth)
if n_elements(prs) eq 0 then prs=0
if prs then prs_sng='prs_' else prs_sng=''

fld3_stf=[ $
;	['Q','!5Vapor Mixing Ratio','!8q!5!Iv!N','!5g kg!E-1!N'], $
;	['RADD','!5Radiative','!8Q!5!IR!N','!5K day!E-1!N'] ]
	['QDIABAT','!5Total','!8Q!5!ITOT!N','!5K day!E-1!N'], $
	['QRS','!5Shortwave','!8Q!5!ISW!N','!5K day!E-1!N'], $
	['QRL','!5Longwave','!8Q!5!ILW!N','!5K day!E-1!N'], $
	['HGS','!5Resolved','!8Q!5!ILS!N','!5K day!E-1!N'], $
	['HDFF','!5Turbulent','!8Q!5!IDFF!N','!5K day!E-1!N'], $
	['CMFDT','!5Convective','!8Q!5!ICNV!N','!5K day!E-1!N'] ]
fld_nbr3=n_elements(fld3_stf(0,*))
foo={fld3_sct,idx:0,nm:'',sng:'',abb:'',unit:''}
fld3=replicate({fld3_sct},fld_nbr3)
for idx=0,fld_nbr3-1 do begin
	fld3(idx).idx=idx
	fld3(idx).nm=fld3_stf(0,idx)
	fld3(idx).sng=fld3_stf(1,idx)
	fld3(idx).abb=fld3_stf(2,idx)
	fld3(idx).unit=fld3_stf(3,idx)
endfor; end loop over fld3s

xpt_stf=[ $
;	['ecmwf','!5ECMWF'], $
;	['amip5','!5CCM!7X!5!I.5!N'], $
	['amip5','!5CCM'], $
;	['sld012d','!5CCM3'], $
	['spcp_85','!5ANV']]
xpt_nbr=n_elements(xpt_stf(0,*))
foo={xpt_sct,idx:0,nm:'',sng:''}
xpt=replicate({xpt_sct},xpt_nbr)
for idx=0,xpt_nbr-1 do begin
	xpt(idx).idx=idx
	xpt(idx).nm=xpt_stf(0,idx)
	xpt(idx).sng=xpt_stf(1,idx)
endfor; end loop over xpts

rgn_stf=[ $
	['Indian_Central','Central Indian Ocean'] ]
;	['Pacific_ITCZ_West','West Pacific ITCZ'] ]
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct_gcm,idx:0,nm:'',sng:''}
rgn=replicate({rgn_sct_gcm},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	rgn(idx).idx=idx
	rgn(idx).nm=rgn_stf(0,idx)
	rgn(idx).sng=rgn_stf(1,idx)
endfor; end loop over rgn_lst

mth_mdm=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
for mth_idx=0,mth_nbr-1 do begin
for rgn_idx=0,rgn_nbr-1 do begin
for xpt_idx=0,xpt_nbr-1 do begin
	ttl=mth_mdm(mth(mth_idx))+' '+xpt(xpt_idx).sng+' !5Diabatic Heating'
;	if prs then y_ttl='!5Pressure (mb)' else y_ttl='!7g!5 !9X!5 1000 (mb)'
	y_ttl='!5Pressure (mb)'
	x_ttl='!5Heating Rate ('+fld3(0).unit+')'
	chr_sz=1.5
	top=1
	btm=1
	info=1
;	if diff then rng_x=[-0.5,0.5] else rng_x=[-5,5]
;	if not diff then rng_x=[-5,5]
	yr_sng='8589_'
	; NB: Jan is mth_srt=0 in C,IDL indexing notation
	mth_sng=string(format='(I2.2)',mth(mth_idx)+1)
	fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_xyavg_rgn_'+rgn(rgn_idx).nm+'_'+yr_sng+'0112.nc'
	fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_xyavg_rgn_'+rgn(rgn_idx).nm+'_'+yr_sng+mth_sng+'.eps'
	if diff then begin ; plot differences between simulated fields
		fl_nm_in=getenv('DATA')+'/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+'0112.nc'
		fl_nm_out=getenv('DATA')+'/ps/'+xpt(xpt_idx).nm+'_'+yr_sng+ctl_xpt+'_'+yr_sng+'xyavg_rgn_'+rgn(rgn_idx).nm+'_'+mth_sng+'.eps'
	endif; endif diff
	fld_nm=strarr(fld_nbr3)
	lbl_sng=strarr(fld_nbr3)
	scl=fltarr(fld_nbr3)+1.0
	for fld_idx3=0,fld_nbr3-1 do begin
		if fld3(fld_idx3).nm eq 'QDIABAT' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'QRS' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'QRL' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'HGS' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'HDFF' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'CMFDT' then scl(fld_idx3)=86400.0
		if fld3(fld_idx3).nm eq 'DTCOND' then scl(fld_idx3)=86400.0
		lbl_sng(fld_idx3)=fld3(fld_idx3).sng
		fld_nm(fld_idx3)=fld3(fld_idx3).nm
	endfor; end loop over field
	if prn and not diff then begin
		ttl=xpt(xpt_idx).sng+' !5Diabatic Heating'
		if mth_idx eq 0 then top=1 else top=0
		if mth_idx eq mth_nbr-1 then btm=1 else btm=0
	endif; endif prn and not diff
	if prn and diff then begin
		ttl='!5ANV - CCM Diabatic Heating'
		if mth_idx eq 0 then top=1 else top=0
		if mth_idx eq mth_nbr-1 then btm=1 else btm=0
	endif; endif prn and diff
	if prn then begin
		x_sz=6.5
		y_sz=3.2
		if top then y_sz=y_sz*1.1
		if btm then y_sz=y_sz*1.2
		info=0
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
	endif; endif prn
	gcm_vz,x_ttl=x_ttl,y_ttl=y_ttl,fl_nm=fl_nm_in,ttl=ttl,prn=prn,fld_nm=fld_nm,rng_x=rng_x,rng_y=[0,1000],chr_sz=chr_sz,top=top,btm=btm,scl=scl,info=info,lbl_sng=lbl_sng,mth=mth(mth_idx)
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over fld3
endfor; end loop over rgn
endfor; end loop over mth

end; end gcm_vz_bch()
pro gcm_vz, $
	fl_nm=fl_nm, $
	rng_x=rng_x, $
	rng_y=rng_y, $
	info=info, $
	prn=prn, $
	top=top, $
	btm=btm, $
	ttl=ttl, $
	mth=mth, $
	chr_sz=chr_sz, $
	x_ttl=x_ttl, $
	y_ttl=y_ttl, $
	lbl_sng=lbl_sng, $
	scl=scl, $
	fld_nm=fld_nm

;gcm_vz,x_ttl='!5Heating Rate !5K day!E-1!N',y_ttl='!5Pressure (mb)',mth=0,fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xyavg_rgn_Pacific_ITCZ_West_8589_0112.nc',ttl='Energy Balance',prn=0,fld_nm=['QRS','QRL','CMFDT'],rng_x=0,rng_y=[1000,100],chr_sz=chr_sz,top=top,btm=btm,scl=[86400,86400,86400],info=1

if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(mth_end) eq 0 then mth_end=12
if n_elements(ttl) eq 0 then ttl='!5West Pacific ITCZ Energy Balance'
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(mth) eq 0 then mth=0
if n_elements(top) eq 0 then top=1
if n_elements(btm) eq 0 then btm=1
if n_elements(info) eq 0 then info=1
if n_elements(prn) eq 0 then prn=0
if n_elements(y_ttl) eq 0 then y_ttl='!5Pressure (mb)'
if n_elements(x_ttl) eq 0 then x_ttl='!5Heating Rate !5K day!E-1!N'
if n_elements(lbl_sng) eq 0 then lbl_sng=['!5Total Diabatic Heating']
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/spcp_85/spcp_85_xyavg_rgn_Pacific_ITCZ_West_8589_0112.nc'
if n_elements(fld_nm) eq 0 then fld_nm=['QDIABAT']
if n_elements(scl) eq 0 then scl=86400

erase
!p.multi=0

print,'Processing '+fl_nm

fld_nbr=n_elements(fld_nm)
if fld_nbr eq 1 then fld_nm=[fld_nm]
for fld_idx=0,fld_nbr-1 do begin

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_nm)
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,dim_foo,lev_nbr
if fld_idx eq 0 then data=fltarr(lev_nbr,fld_nbr)
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,fld_nm(fld_idx),data_foo
data_foo=reform(data_foo)
data(*,fld_idx)=data_foo(*,mth) ; CCM disk data appears to IDL as data(lon,lev,lat,time)
ncdf_close,nc_id
; End of NetCDF commands

; Scaling missing data will produce problems
data(*,fld_idx)=scl(fld_idx)*data(*,fld_idx)

endfor; end loop over fields

good_idx=where(data lt 1.0e20)
data(good_idx)=data(good_idx)
abc=data
ord=lev

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc(where(abc lt 1.0e20)))
	abc_max=max(abc(where(abc lt 1.0e20)))
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

for fld_idx=0,fld_nbr-1 do begin
mss_val=where(abc(*,fld_idx) gt 1.0e20,cnt)
if cnt gt 0 then ord(mss_val)=abc(mss_val,fld_idx)
if fld_idx eq 0 then begin
plot, $
	abc(*,fld_idx), $
	ord, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=1, $
	ystyle=1, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_max,ord_min], $ ; Pressure decreases on y-axis
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	linestyle=0

endif else begin; fld_idx != 0
oplot, $
	abc(*,fld_idx), $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fld_idx

endelse; fld_idx != 0
endfor; end loop over files

if info then begin

ln_lgn_x1=.2
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for fld_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fld_idx)+0.013,linestyle=fld_idx,thick=2.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fld_idx),lbl_sng(fld_idx),size=chr_sz,/NORMAL
endfor; end loop over fld

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_vz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM VZ commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
