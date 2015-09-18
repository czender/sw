; $Id$

@ibp_clr.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM DNS1
; Overplot the same timeseries from experiments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_dns1_bch, $
	ann=ann, $
	wdw=wdw, $
	prn=prn, $
	spc=spc
if n_elements(ann) eq 0 then ann=0
if n_elements(wdw) eq 0 then wdw=0
if n_elements(prn) eq 0 then prn=0
if n_elements(spc) eq 0 then spc=0

rgn_stf=[ $
	['Pacific_Tropical','Tropical Pacific'], $
	['Atlantic_Tropical','Tropical Atlantic'], $
	['Indian_Tropical','Tropical Indian Ocean'], $
	['Ocean_Tropical','Tropical Ocean']]
;	['Pacific_Tropical','Tropical Pacific']]
;	['Pacific_Equatorial','Equatorial Pacific']]
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct,idx:0,nm:'',sng:''}
reg=replicate({rgn_sct},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	reg(idx).idx=idx
	reg(idx).nm=rgn_stf(0,idx)
	reg(idx).sng=rgn_stf(1,idx)
endfor; end loop over fld_lst

fld2_stf=[ $
	['TS1','!8SST','!E!12_!5!NK'], $
	['GCLR','!8G!5!Ia!N','!5W m!E-2!N'], $
	['LWCF','!8LWCF!5','!5W m!E-2!N'], $
	['GCLD','!8G!5','!5W m!E-2!N'], $
	['FLNT','!8OLR!5','!5W m!E-2!N'], $
	['ALBEDO','!5Albedo','!5%'], $
	['SWCF','!8SWCF!5','!5W m!E-2!N'] ]
;	['LWCF','!8LWCF!5','!5W m!E-2!N'] ]
;	['FLNT','!8OLR!5','!5W m!E-2!N'] ]
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
;	['erbe_b','!5ERBE'] ]
	['erbe_b','!5ERBE'], $
	['amip5','!5CCM'], $
	['spcp_81','!5ANV']]
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
mth_end=59
; NB: Store files with Fortran indexing notation
mth_srt_sng=string(format='(I2.2)',mth_srt+1)
mth_end_sng=string(format='(I2.2)',mth_end+1)
for fld_idx2=0,fld_nbr2-1 do begin
for rgn_idx=0,rgn_nbr-1 do begin
;	Re-initialize things that might change
	chr_sz=1.5
	top=1
	btm=1
	scl=1.
	rng_y=0
	info=1
	ttl=reg(rgn_idx).sng+' '+' 1985-1989'
	if xpt_nbr eq 1 then ttl=xpt(0).sng+' '+ttl
	if xpt_nbr eq 1 then info=0
	y_ttl=fld2(fld_idx2).sng+' ('+fld2(fld_idx2).unit+')'
	x_ttl=''
	if strpos(fld2(fld_idx2).nm,'ALBEDO') ne -1 then scl=100.
	fl_nm_out='/data/zender/ps/anom_xyavg_rgn_'+reg(rgn_idx).nm+'_8589_'+mth_srt_sng+mth_end_sng+'_'+fld2(fld_idx2).nm+'.eps'
	if prn then begin
		x_sz=6.5
		y_sz=2
		chr_sz=1.5
		if fld_idx2 eq 0 then top=1 else top=0
		if fld_idx2 eq fld_nbr2-1 then btm=1 else btm=0
		if top then y_sz=y_sz*1.1
		if btm then y_sz=y_sz*1.1
		if top then info=1 else info=0
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
	endif; endif prn
	fl_nm_in=strarr(xpt_nbr)	
	lbl_sng=strarr(xpt_nbr)	
	for xpt_idx=0,xpt_nbr-1 do begin
		fl_nm_in(xpt_idx)='/data/zender/_aux0_/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_anom_xyavg_rgn_'+reg(rgn_idx).nm+'_8589_'+mth_srt_sng+mth_end_sng+'.nc'
		lbl_sng(xpt_idx)=xpt(xpt_idx).sng
	endfor; end loop over experiments
	gcm_dns1,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,ann=ann,fld_nm=fld2(fld_idx2).nm,wrp=0,chr_sz=chr_sz,rng_y=rng_y,scl=scl,top=top,btm=btm,spc=spc,wdw=wdw,lbl_sng=lbl_sng,info=info
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over rgn_lst
endfor; end loop over fld2s

end; end gcm_dns1_bch()
pro gcm_dns1, $
	btm=btm, $
	chr_sz=chr_sz, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	mth_end=mth_end, $
	mth_srt=mth_srt, $
	prn=prn, $
	ann=ann, $
	info=info, $
	lbl_sng=lbl_sng, $
	rgn_nm=rgn_nm, $
	scl=scl, $
	spc=spc, $
	wdw=wdw, $
	top=top, $
	ttl=ttl, $
	wrp=wrp, $
	rng_x=rng_x, $
	x_ttl=x_ttl, $
	rng_y=rng_y, $
	y_ttl=y_ttl

if n_elements(btm) eq 0 then btm=1
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/spcp_81/spcp_81_anom_xyavg_rgn_Pacific_Tropical_8589_0160.nc'
if n_elements(fld_nm) eq 0 then fld_nm='LWCF'
if n_elements(mth_end) eq 0 then mth_end=59
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(ann) eq 0 then ann=0
if n_elements(rgn_nm) eq 0 then rgn_nm='Pacific_Tropical'
if n_elements(lbl_sng) eq 0 then lbl_sng=''
if n_elements(scl) eq 0 then scl=1.
if n_elements(spc) eq 0 then spc=0
if n_elements(wdw) eq 0 then wdw=0
if n_elements(top) eq 0 then top=1
if n_elements(ttl) eq 0 then ttl='!5Tropical Pacific'
if n_elements(wrp) eq 0 then wrp=0
if n_elements(x_ttl) eq 0 then x_ttl='!5Year'
if n_elements(y_ttl) eq 0 then y_ttl='!5LWCF (!5W m!e-2!N)'

erase
!p.multi=0

fl_nbr=n_elements(fl_nm)
if fl_nbr eq 1 then fl_nm=[fl_nm]
for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)
nc_id=ncdf_open(fl_nm(fl_idx))
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
if fl_idx eq 0 then data=fltarr(time_nbr,fl_nbr)
ncdf_varget,nc_id,fld_nm,data_foo
data(*,fl_idx)=reform(data_foo)
ncdf_close,nc_id

endfor; end loop over files

good_idx=where(data lt 1.0e20)
data(good_idx)=scl*data(good_idx)
abc=findgen(time_nbr)

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse

; This loop is allowed to change the data so save the original
data_orig=data
data_plt=fltarr(time_nbr,fl_nbr)
abc_plt=fltarr(time_nbr,fl_nbr)
idx_srt=intarr(fl_nbr)
idx_end=intarr(fl_nbr)+time_nbr-1
for fl_idx=0,fl_nbr-1 do begin
data=data_orig(*,fl_idx)
abc=findgen(time_nbr)
if strpos(fl_nm(fl_idx),'erbe') ne -1 and max(abc) eq 59 then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; Do not plot ERBE data before 8501 or after 8905.
; NB: if spc processing, truncate ERBE to 48 months to keep leakage down
idx_srt(fl_idx)=1
if spc or ann then idx_end(fl_idx)=48 else idx_end(fl_idx)=52
if spc or ann then abc=abc(1:48) else abc=abc(1:52)
if spc or ann then data=data(1:48) else data=data(1:52)
endif
data_nbr=n_elements(data)
; QC: the data should sum to 0 if they represent all the anomalies from the baseline
print,'QC: total(data)/data_nbr = ',total(data)/data_nbr

if ann then begin
; Look at the annual cycle
	data_spc=fft(data,-1)	
	idx_ann=round(data_nbr/12.)
	amp=2.*abs(data_spc(idx_ann))
	print,'(Half) Amplitude annual cycle = 2*abs(H('+auto_sng(idx_ann,0)+')) = '+auto_sng(amp,2)
	; Reverse sign of angle because default IDL FFT has opposite convention to PFT88
	phi=-atan(imaginary(data_spc(idx_ann)),float(data_spc(idx_ann))) ; returns -pi to pi
	; If this is ERBE data starting in 8502 then correct so phi=0 refers to 8501
	if strpos(fl_nm(fl_idx),'erbe') ne -1 then phi=phi+!pi/6
	if phi gt 2*!pi then phi=phi-2*!pi else if phi lt 0 then phi=phi+2*!pi; 0 to 2*pi
	print,'Phase is -atan(Im(H_ann),Re(H_ann)) + phi_0 = '+auto_sng(180*phi/!pi,2)+' degrees'
	if fld_nm eq 'TS1' then amp_SST=amp
;	data_spc(idx_ann)=0 ; Removes cosine component
;	data_spc(data_nbr-idx_ann)=0 ; Removes sine component
;	data=float(fft(data_spc,1)) ; Transforms back
endif; endif ann
if spc then begin
; Perform various spectral diagnostics
	; Test data should return amp=5,phi=!pi/4
	;data=5*cos(2*!pi*findgen(data_nbr)/12  - !pi/4)
	; NB: data can be changed by the window
	tmp_data=data
	gcm_spc,abc=abc,data=tmp_data,unit_sng='month',prn=0,wdw=wdw,amp=amp,phi=phi
	; If this is ERBE data starting in 8502 then correct so phi=0 refers to 8501
	if strpos(fl_nm(fl_idx),'erbe') ne -1 then phi=phi+!pi/6
	if phi gt 2*!pi then phi=phi-2*!pi ; 0 to 2*pi
	; Remove the (possibly windowed) annual cycle from the data
	data=data-amp*cos((2*!pi*findgen(data_nbr)/12) - phi)
endif; endif spc
data_plt(idx_srt(fl_idx):idx_end(fl_idx),fl_idx)=data
abc_plt(idx_srt(fl_idx):idx_end(fl_idx),fl_idx)=abc
endfor; end loop over files
data=data_plt
abc=abc_plt

ord_min=1.0e36
ord_max=-1.0e36
if n_elements(rng_y) le 1 then begin
	for fl_idx=0,fl_nbr-1 do ord_min=min([ord_min,min(data(idx_srt(fl_idx):idx_end(fl_idx),fl_idx))])
	for fl_idx=0,fl_nbr-1 do ord_max=max([ord_max,max(data(idx_srt(fl_idx):idx_end(fl_idx),fl_idx))])
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

mrg_top=2.5 ; 2 is default
mrg_btm=4 ; 4 is default
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.5 ; 2 is default
	if top then mrg_top=mrg_top+1.3
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+1.3
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

for fl_idx=0,fl_nbr-1 do begin
if fl_nbr eq 2 and fl_idx eq 1 then ln_sty=fl_idx+1 else ln_sty=fl_idx
if fl_idx eq 0 then begin
plot, $
	abc(idx_srt(fl_idx):idx_end(fl_idx),fl_idx), $
	data(idx_srt(fl_idx):idx_end(fl_idx),fl_idx), $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=12, $
	ystyle=0, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	linestyle=ln_sty
endif else begin; fl_idx != 0
oplot,	$
	abc(idx_srt(fl_idx):idx_end(fl_idx),fl_idx), $
	data(idx_srt(fl_idx):idx_end(fl_idx),fl_idx), $
	thick=3.0, $
	linestyle=ln_sty

endelse; fl_idx != 0
endfor; end loop over files

if prn then tck_lng=0.04 else tck_lng=0.02
if btm then time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=1,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0 else time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=0,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0

if info then begin

ln_lgn_x1=.75
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.70
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]

for fl_idx=0,fl_nbr-1 do begin
if fl_nbr eq 2 and fl_idx eq 1 then ln_sty=fl_idx+1 else ln_sty=fl_idx
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty,thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),size=chr_sz,/NORMAL
endfor; end loop over fl

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_dns1()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM DNS1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM DNS2
; Overplot multiple timeseries from a single experiment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_dns2_bch, $
	ann=ann, $
	anom=anom, $
	wdw=wdw, $
	prn=prn, $
	sst=sst, $
	spc=spc
if n_elements(ann) eq 0 then ann=0
if n_elements(wdw) eq 0 then wdw=0
if n_elements(prn) eq 0 then prn=0
if n_elements(sst) eq 0 then sst=0
if n_elements(spc) eq 0 then spc=0
if n_elements(anom) eq 0 then anom=0
if anom then anom_sng='anom_' else anom_sng=''

rgn_stf=[ $
	['Pacific_Tropical','Tropical Pacific'], $
	['Atlantic_Tropical','Tropical Atlantic'], $
	['Indian_Tropical','Tropical Indian Ocean'], $
	['Ocean_Tropical','Tropical Ocean']]
;	['Pacific_Tropical','Tropical Pacific']]
;	['Pacific_Equatorial','Equatorial Pacific']]
rgn_nbr=n_elements(rgn_stf(0,*))
foo={rgn_sct,idx:0,nm:'',sng:''}
reg=replicate({rgn_sct},rgn_nbr)
for idx=0,rgn_nbr-1 do begin
	reg(idx).idx=idx
	reg(idx).nm=rgn_stf(0,idx)
	reg(idx).sng=rgn_stf(1,idx)
endfor; end loop over fld_lst

fld2_stf=[ $
	['TS1','!8SST','!E!12_!5!NK'], $
	['GCLR','!8G!5!Ia!N','!5W m!E-2!N'], $
	['LWCF','!8LWCF!5','!5W m!E-2!N'], $
	['GCLD','!8G!5','!5W m!E-2!N'], $
	['FLNT','!8OLR!5','!5W m!E-2!N'], $
;	['ALBEDO','!5Albedo','!5%'], $
	['SWCF','!8SWCF!5','!5W m!E-2!N'] ]
;	['LWCF','!8LWCF!5','!5W m!E-2!N'] ]
;	['FLNT','!8OLR!5','!5W m!E-2!N'] ]
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
;	['erbe_b','!5ERBE'] ]
	['erbe_b','!5ERBE'], $
	['amip5','!5CCM'], $
	['spcp_81','!5ANV']]
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
mth_end=59
; NB: Store files with Fortran indexing notation
mth_srt_sng=string(format='(I2.2)',mth_srt+1)
mth_end_sng=string(format='(I2.2)',mth_end+1)
for xpt_idx=0,xpt_nbr-1 do begin
for rgn_idx=0,rgn_nbr-1 do begin
;	Re-initialize things that might change
	pll=3
	chr_sz=1.5
	top=1
	btm=1
	rng_y=0
	rng_SST=0
	info=1
	ttl=reg(rgn_idx).sng+' '+' 1985-1989'
	if fld_nbr2 eq 1 then ttl=fld2(0).sng+' '+ttl
	if fld_nbr2 eq 1 then info=0
	if not prn then ttl=xpt(xpt_idx).sng+' '+ttl
	y_ttl=' ('+fld2(fld_nbr2-1).unit+')'
	x_ttl=''
	if reg(rgn_idx).nm eq 'Indian_Tropical' then if anom then rng_y=[-10,10] else rng_y=[-10,10]
	if reg(rgn_idx).nm eq 'Atlantic_Tropical' then if anom then rng_y=[-8,8] else rng_y=[-8,8]
	if reg(rgn_idx).nm eq 'Pacific_Tropical' then if anom then rng_y=[-6,6] else rng_y=[-6,6]
	if reg(rgn_idx).nm eq 'Ocean_Tropical' then if anom then rng_y=[-4,4] else rng_y=[-4,4]
	if reg(rgn_idx).nm eq 'Ocean_Tropical' then if anom then rng_SST=[-.4,.4] else rng_SST=[-.4,.4]
	fl_nm_in='/data/zender/_aux0_/'+xpt(xpt_idx).nm+'/'+xpt(xpt_idx).nm+'_'+anom_sng+'xyavg_rgn_'+reg(rgn_idx).nm+'_8589_'+mth_srt_sng+mth_end_sng+'.nc'
	fl_nm_out='/data/zender/ps/'+xpt(xpt_idx).nm+'_'+anom_sng+'xyavg_rgn_'+reg(rgn_idx).nm+'_8589_'+mth_srt_sng+mth_end_sng+'.eps'
	if prn then begin
		x_sz=6.5
		y_sz=2
		chr_sz=1.5
		if xpt_idx eq 0 then top=1 else top=0
		if xpt_idx eq xpt_nbr-1 then btm=1 else btm=0
		if top then y_sz=y_sz*1.1
		if btm then y_sz=y_sz*1.1
		if top then info=1 else info=0
		open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/color,/eps
	endif; endif prn
	fld_nm=strarr(fld_nbr2)	
	lbl_sng=strarr(fld_nbr2)	
	scl=fltarr(fld_nbr2)+1.
	for fld_idx2=0,fld_nbr2-1 do begin
		fld_nm(fld_idx2)=fld2(fld_idx2).nm
		lbl_sng(fld_idx2)=fld2(fld_idx2).sng
		if strpos(fld2(fld_idx2).nm,'ALBEDO') ne -1 then scl(fld_idx2)=100.
	endfor; end loop over fld2s
	gcm_dns2,x_ttl=x_ttl,y_ttl=y_ttl,mth_srt=mth_srt,mth_end=mth_end,fl_nm=fl_nm_in,ttl=ttl,prn=prn,ann=ann,fld_nm=fld_nm,wrp=0,chr_sz=chr_sz,rng_y=rng_y,rng_SST=rng_SST,scl=scl,top=top,btm=btm,spc=spc,anom=anom,wdw=wdw,lbl_sng=lbl_sng,info=info,pll=pll,order=1,sst=sst
	if prn then close_ps,fl_nm=fl_nm_out
endfor; end loop over rgn_lst
endfor; end loop over experiments

end; end gcm_dns2_bch()
pro gcm_dns2, $
	btm=btm, $
	chr_sz=chr_sz, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	mth_end=mth_end, $
	mth_srt=mth_srt, $
	prn=prn, $
	ann=ann, $
	info=info, $
	order=order, $
	pll=pll, $
	lbl_sng=lbl_sng, $
	rgn_nm=rgn_nm, $
	scl=scl, $
	spc=spc, $
	sst=sst, $
	anom=anom, $
	wdw=wdw, $
	top=top, $
	ttl=ttl, $
	wrp=wrp, $
	rng_x=rng_x, $
	x_ttl=x_ttl, $
	rng_y=rng_y, $
	rng_SST=rng_SST, $
	y_ttl=y_ttl

@ibp_clr.com

if n_elements(btm) eq 0 then btm=1
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/spcp_81/spcp_81_anom_xyavg_rgn_Pacific_Tropical_8589_0160.nc'
if n_elements(fld_nm) eq 0 then fld_nm='LWCF'
if n_elements(mth_end) eq 0 then mth_end=59
if n_elements(mth_srt) eq 0 then mth_srt=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(pll) eq 0 then pll=4
if n_elements(order) ne 0 then color_order=order else color_order=1
if n_elements(clr_tbl) eq 0 then clr_tbl=22
if n_elements(ann) eq 0 then ann=0
if n_elements(rgn_nm) eq 0 then rgn_nm='Pacific_Tropical'
if n_elements(lbl_sng) eq 0 then lbl_sng=''
if n_elements(scl) eq 0 then scl=1.
if n_elements(spc) eq 0 then spc=0
if n_elements(sst) eq 0 then sst=1
if n_elements(anom) eq 0 then anom=0
if n_elements(wdw) eq 0 then wdw=0
if n_elements(top) eq 0 then top=1
if n_elements(ttl) eq 0 then ttl='!5Tropical Pacific'
if n_elements(wrp) eq 0 then wrp=0
if n_elements(x_ttl) eq 0 then x_ttl='!5Year'
if n_elements(y_ttl) eq 0 then y_ttl='!5LWCF (!5W m!e-2!N)'

clr_rfr,clr_tbl,"",0
erase
!p.multi=0

fld_nbr=n_elements(fld_nm)
if fld_nbr eq 1 then fld_nm=[fld_nm]

print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
for fld_idx=0,fld_nbr-1 do begin
if fld_idx eq 0 then data=fltarr(time_nbr,fld_nbr)
ncdf_varget,nc_id,fld_nm(fld_idx),data_foo
data(*,fld_idx)=reform(data_foo)

; Scaling missing data will produce problems
data(*,fld_idx)=scl(fld_idx)*data(*,fld_idx)

endfor; end loop over fields
ncdf_close,nc_id

good_idx=where(data lt 1.0e20)
data(good_idx)=data(good_idx)
abc=findgen(time_nbr)

if spc or ann then begin
; Assume we are dealing with a 8589_0160 style file, i.e., 60 months of data
; NB: if spc processing, truncate to 48 months to keep leakage down and compatability with ERBE
abc=abc(1:48)
data=data(1:48,*)
endif; endif spc or ann
if strpos(fl_nm,'erbe') ne -1 and max(abc) gt 53 and not spc and not ann then begin
; Never plot ERBE data before 8501 or after 8905.
abc=abc(1:52)
data=data(1:52,*)
endif; endif ERBE
data_nbr=n_elements(abc)

; This loop is allowed to change the data so save the original
data_plt=data
amp_arr=fltarr(fld_nbr)
phi_arr=fltarr(fld_nbr)
avg_arr=fltarr(fld_nbr)
amp=0
phi=0
avg=0
for fld_idx=0,fld_nbr-1 do begin
data=data_plt(*,fld_idx)
; QC: the data should sum to 0 if they represent all the anomalies from the baseline
;print,'QC: total(data)/data_nbr = ',total(data)/data_nbr
if ann then begin
; Look at the annual cycle
	data_spc=fft(data,-1)	
	idx_ann=round(data_nbr/12.)
	amp=2.*abs(data_spc(idx_ann))
	;print,'(Half) Amplitude annual cycle = 2*abs(H('+auto_sng(idx_ann,0)+')) = '+auto_sng(amp,2)
	; Reverse sign of angle because default IDL FFT has opposite convention to PFT88
	phi=-atan(imaginary(data_spc(idx_ann)),float(data_spc(idx_ann))) ; returns -pi to pi
	; If this is ERBE data starting in 850215 then correct so phi=0 refers to 850115 (add 30 degrees per mth)
	phi=phi+abc(0)*!pi/6
	; Correct data from phi=0 on 850115 to phi=0 on 850101 (add 15 degrees)
	phi=phi+!pi/12
	if phi gt 2*!pi then phi=phi-2*!pi else if phi lt 0 then phi=phi+2*!pi; 0 to 2*pi
	;print,'Phase is -atan(Im(H_ann),Re(H_ann)) + phi_0 = '+auto_sng(180*phi/!pi,2)+' degrees'
	; Save the average
	avg=float(data_spc(0))
	if fld_nm(fld_idx) eq 'TS1' then amp_SST=amp
	if fld_nm(fld_idx) eq 'TS1' then phi_SST=phi
	; NB: Remember that data will now be synchronized to monthly interfaces not monthly midpoints
	data=amp*cos((2*!pi*findgen(data_nbr)/12) - phi)
	if not anom then data=data+avg
;	data_spc(idx_ann)=0 ; Removes cosine component
;	data_spc(data_nbr-idx_ann)=0 ; Removes sine component
;	data=float(fft(data_spc,1)) ; Transforms back
endif; endif ann
if spc then begin
; Perform various spectral diagnostics
	; Test data should return amp=5,phi=!pi/4
	;data=5*cos(2*!pi*findgen(data_nbr)/12  - !pi/4)
	; NB: data can be changed by the window
	tmp_data=data
	gcm_spc,abc=abc,data=tmp_data,unit_sng='month',prn=0,wdw=wdw,amp=amp,phi=phi
	; If this is ERBE data starting in 850215 then correct so phi=0 refers to 850115 (add 30 degrees per mth)
	phi=phi+abc(0)*!pi/6
	if phi gt 2*!pi then phi=phi-2*!pi ; 0 to 2*pi
	; Remove the (possibly windowed) annual cycle from the data
	data=data-amp*cos((2*!pi*findgen(data_nbr)/12) - phi)
endif; endif spc
amp_arr(fld_idx)=amp
phi_arr(fld_idx)=phi
avg_arr(fld_idx)=avg
data_plt(*,fld_idx)=data
endfor; end loop over fields
data=data_plt
amp=amp_arr
phi=phi_arr
avg=avg_arr

if ann then begin
wrp=1
data_nbr=12
; NB: only reset abc after phi has been corrected for initial offset
abc=findgen(data_nbr)
data=fltarr(data_nbr,fld_nbr)
print,'Peak to peak annual cycle amplitude and phase SST = '+auto_sng(2*amp_SST,2)+' K, '+auto_sng(round(r2d(phi_SST)),0)+' degrees'
print
for fld_idx=0,fld_nbr-1 do begin
	data(*,fld_idx)=amp(fld_idx)*cos((2*!pi*findgen(data_nbr)/12) - phi(fld_idx))
	if not anom then data(*,fld_idx)=data(*,fld_idx)+avg(fld_idx)
	if fld_nm(fld_idx) ne 'TS1' then print,'Peak to peak annual cycle amplitude and phase '+lbl_sng(fld_idx)+' = '+auto_sng(2*amp(fld_idx),2)+' W/m2, '+auto_sng(round(r2d(phi(fld_idx))),0)+' degrees'
	if fld_nm(fld_idx) ne 'TS1' then print,'Ratio '+lbl_sng(fld_idx)+'/SST = '+auto_sng(amp(fld_idx)/amp_SST,2)+' W/m2/K'
endfor; end loop over fields
if wrp then begin
data=[data,data(0,*)]
abc=[abc,data_nbr]
endif
endif; endif ann

if n_elements(rng_x) le 1 then begin
	abc_min=min(abc)
	abc_max=max(abc)
endif else begin
	abc_min=rng_x(0)
	abc_max=rng_x(1)
endelse
ord_min=1.0e36
ord_max=-1.0e36
for fld_idx=0,fld_nbr-1 do begin
if sst then begin
if n_elements(rng_SST) le 1 then begin
	if fld_nm(fld_idx) eq 'TS1' then ord_min_SST=min(data(*,fld_idx))
	if fld_nm(fld_idx) eq 'TS1' then ord_max_SST=max(data(*,fld_idx))
endif else begin
	ord_min_SST=rng_SST(0)
	ord_max_SST=rng_SST(1)
endelse; !rng_sst
endif; sst
if n_elements(rng_y) le 1 then begin
	ord_min=min([ord_min,min(data(*,fld_idx))])
	ord_max=max([ord_max,max(data(*,fld_idx))])
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse
endfor; end loop over fields

mrg_top=2.5 ; 2 is default
mrg_btm=4 ; 4 is default
if sst then mrg_lft=8 else mrg_lft=8 ; 10 is default
if sst then mrg_rgt=6 else mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.5 ; 2 is default
	if top then mrg_top=mrg_top+1.3
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+1.3
	if not top then ttl=''
	if not btm then x_ttl=''
endif; endif prn

clr_mk,fld_nbr,color_order,pll,0
tvlct,r_curr,g_curr,b_curr
clr=indgen(fld_nbr)+2
for fld_idx=0,fld_nbr-1 do begin
if fld_nm(fld_idx) eq 'TS1' then clr(fld_idx)=clr_blk_idx else clr(fld_idx)=fld_idx+2
endfor; end loop over fields

for fld_idx=0,fld_nbr-1 do begin
if fld_nbr eq 2 and fld_idx eq 1 then ln_sty=fld_idx+1 else ln_sty=fld_idx
if fld_idx eq 0 then begin
if not sst then plot, $
	abc, $
	data(*,fld_idx), $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit=y_ttl, $
	xstyle=13, $
	ystyle=0, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	color=clr(fld_idx), $
	linestyle=0
if sst then plot, $
	abc, $
	data(*,fld_idx), $
	max_value=1.0e20, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit='!8 SST (!E!12_!5!NK)', $
	xstyle=13, $
	ystyle=8, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min_SST,ord_max_SST], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	color=clr(fld_idx), $
	linestyle=0
endif else begin; fld_idx != 0
if sst and fld_idx eq 1 then begin ; plot first w/m2 field and right hand axis
plot, $
	abc, $
	data(*,fld_idx), $
	max_value=1.0e20, $
	xstyle=13, $
	ystyle=4, $
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[mrg_lft,mrg_rgt], $
	ymargin=[mrg_btm,mrg_top], $
	color=clr(fld_idx), $
	/noerase, $
	linestyle=0
axis, $
	yaxis=1, $
	ystyle=0, $
	ytitle=y_ttl, $
	charsize=chr_sz

endif else begin; fld_idx != 0 or not sst
oplot,	$
	abc, $
	data(*,fld_idx), $
	thick=3.0, $
	color=clr(fld_idx), $
	linestyle=0
endelse ; fld_idx != 0 or not sst
endelse ; fld_idx != 0
endfor ; end loop over fields

if prn then tck_lng=0.04 else tck_lng=0.02
if btm then time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=1,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0 else time_axz_drw,time_min=!x.crange(0),time_max=!x.crange(1),time_ncr=1,time_unit='mth',ntt=0,ttl=x_ttl,chr_sz=chr_sz,tck_lng=tck_lng,dbg_lvl=0,lbl_sty='mth_shrt',axz_vrt=0

if info then begin

if sst then ln_lgn_x1=0.5 else ln_lgn_x1=.75
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
if prn then lgn_y_top=.70 else lgn_y_top=0.85
lgn_dy=0.1
lgn_y=lgn_y_top-lgn_dy*findgen(fld_nbr)

for fld_idx=0,fld_nbr-1 do begin
if fld_nbr eq 2 and fld_idx eq 1 then ln_sty=fld_idx+1 else ln_sty=fld_idx
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fld_idx)+0.013,linestyle=0,color=clr(fld_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fld_idx),lbl_sng(fld_idx),size=chr_sz,/NORMAL
endfor; end loop over fl

endif; endif info

if not prn then begin
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif

end; end gcm_dns2()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM DNS2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM SPC
; Compute spectral properties of timeseries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_spc_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

time_nbr=48
period=12
abc=findgen(time_nbr)+1
data=sin(2*!pi*(findgen(time_nbr)/period))
unit_sng='month'
chr_sz=1
wdw=1
amp=0
phi=0

gcm_spc,data=data,abc=abc,prn=prn,unit_sng=unit_sng,chr_sz=chr_sz,wdw=wdw,amp=amp,phi=phi

end; end gcm_spc_bch()
pro gcm_spc, $
	wdw=wdw, $
	amp=amp, $
	phi=phi, $
	abc=abc, $
	data=data, $
	chr_sz=chr_sz, $
	unit_sng=unit_sng, $
	prn=prn
if n_elements(prn) eq 0 then prn=0
if n_elements(abc) eq 0 then abc=2*!pi*findgen(12)/12
if n_elements(wdw) eq 0 then wdw=0
if n_elements(amp) eq 0 then amp=0
if n_elements(phi) eq 0 then phi=0
if n_elements(data) eq 0 then data=sin(abc)
if n_elements(chr_sz) eq 0 then chr_sz=1
if n_elements(unit_sng) eq 0 then unit_sng='month'

; Ref. Num. Rec. p. 405
; Let N = # time samples = data_nbr 
;     dlt = time interval between samples = 1 month, e.g. NB: the units of dlt determine the units of everything else
;     nyquist = nyquist frequency = 1/(2*dlt)
;     frq = frequencies at which transform is evaluated
data_nbr=n_elements(data)
N_over_2=data_nbr/2
dlt=1 ; month
frq_sng='!5month!E-1!N'
frq_nyquist=0.5/dlt
; Pad for an extra spot for symmetry
frq=fltarr(data_nbr+1)
for idx_frq=0,data_nbr do begin
	n=-N_over_2+idx_frq
	frq(idx_frq)=float(n)/(data_nbr*dlt)
endfor; end loop over idx_frq

if wdw then begin
; window the data with the Welch filter
welch=findgen(data_nbr)-.5*(data_nbr-1)
welch=2*welch/(data_nbr+1)
welch=1-welch*welch
;plot,indgen(data_nbr),welch
;print,' Hit any key to continue...'
;junk = get_kbrd(1)
data=data*welch
endif; endif wdw

; *_spc refers to the transform in the scrambled fft order (dimensioned data_nbr = N) (see PFT88 p. 406)
; *_frq refers to amplitudes in the unscrambled frequency space order (dimensioned data_nbr+1 = N+1 because of symmetry padding)
data_spc=fft(data,-1)
print,'QC: Mean of fnc is',total(data)/data_nbr,' and FFT mean = Re(H_0) = ',float(data_spc(0))

; Compute the periodigram as per PFT88 p. 439
; NB: symmetries of returned FFT match the needs of the pgm exactly
frq_pgm=frq(N_over_2:data_nbr)
pgm=fltarr(N_over_2+1)
pgm(0)=abs(data_spc(0))^2
for idx_frq=1,N_over_2-1 do begin
	pgm(idx_frq)=abs(data_spc(idx_frq))^2+abs(data_spc(data_nbr-idx_frq))^2 
endfor; endloop over idx_frq
pgm(N_over_2)=abs(data_spc(N_over_2))^2

; Normalize the power
pgm=data_nbr*pgm/total(data*data) ; NB: multiply by data_nbr in Parseval's theorem because IDL FFT includes 1/N factor already (window is implicitly part of the data now, no need to separately normalize)
pgm_max=max(pgm,idx_max)
print,'max power is '+auto_sng(pgm_max*100,2)+'% in wavenumber '+auto_sng(idx_max)+' = '+auto_sng(frq_pgm(idx_max),2)+' per '+unit_sng+' = period '+auto_sng(1./frq_pgm(idx_max),2)+' '+unit_sng
print,'QC: sum of normalized power periodigram is '+auto_sng(total(pgm),2)

; Compute amplitude and phase of annual cycle
idx_ann=round(data_nbr/12.)
amp=2.*abs(data_spc(idx_ann))
print,'(Half) Amplitude annual cycle = 2*abs(H('+auto_sng(idx_ann,0)+')) = '+auto_sng(amp,2)
; Reverse sign of angle because default IDL FFT has opposite convention to PFT88
phi=-atan(imaginary(data_spc(idx_ann)),float(data_spc(idx_ann))) ; returns -pi to pi
if phi lt 0 then phi=phi+2*!pi ; 0 to 2*pi
print,'Phase is -atan(Im(H_ann),Re(H_ann)) = '+auto_sng(180*phi/!pi,2)+' degrees'
; Phase now in months since Jan for CCM data, months since Feb for ERBE data

;if not prn then !p.multi=[0,1,2]
!p.multi=[0,1,3]
plot,abc,amp*cos((2*!pi*abc/12) - phi),xrange=[0,60],xsty=1
;plot,abc,data,ytit='Value',xtit='Time ('+unit_sng+')'

frq_sng=unit_sng+'!E-1!N'
plot_io, $
	frq_pgm, $
	pgm, $
	ytit='Normalized Power', $
	xtit='Frequency ('+frq_sng+')', $
	xtick_get=x_tck_crd, $
	xstyle=1, $
	yrange=[1.0e-3,1], $
	psym=-6, $ ; plot data as squares connected by lines
	ticklen=1, $
	charsize=chr_sz, $
	xmargin=[8.75,1.5], $ ;[10,3] is [left,right] default
	linestyle=1

x_tck_nbr=n_elements(x_tck_crd)
x_sng=strarr(x_tck_nbr)
for idx_tck=0,x_tck_nbr-1 do begin
	period=x_tck_crd(idx_tck)
	if period eq 0. then x_sng(idx_tck)='!9$!5' else begin
		period=1./period ; invert frequency to get period
		x_sng(idx_tck)=auto_sng(period,2)
	endelse
endfor

axis, $
	xaxis=1, $
	xticks=x_tck_nbr-1, $
	xtickn=x_sng, $
	xtitle='!5Period !7s!5 ('+unit_sng+')', $
	xmargin=[8.75,1.5], $ ;[10,3] is [left,right] default
	charsize=chr_sz

; Plot sine and cosine components of the transform on the same graph in linear frequency space
if(0) then begin

; Unscramble the transform
data_frq=complexarr(data_nbr+1)
data_frq(N_over_2)=data_spc(0) ; the zero frequency component (the average of the data)
data_frq(data_nbr)=data_spc(N_over_2) ; the +'ve nyquist frequency component
data_frq(0)=data_spc(N_over_2) ; the -'ve nyquist frequency component

; Get the positive frequency components
for idx_frq=N_over_2+1,data_nbr do begin  ; from smallest positive frequency to +'ve nyquist frequency
	idx_spc=idx_frq-N_over_2
	data_frq(idx_frq)=data_spc(idx_spc)
endfor; end loop over idx_frq

; Get the negative frequency components
for idx_frq=0,N_over_2-1 do begin ; from -'ve nyquist frequency to smallest negative frequency
	idx_spc=N_over_2+idx_frq
	data_frq(idx_frq)=data_spc(idx_spc)
endfor; end loop over idx_frq

; NB: It only makes sense to truncate the abcissa axis at the Nyquist frequencies, however, the IDL tickmarks
; can then be staggered from the endpoints leading to the wrong tickmarks on the period axis, so use xstyle=0
; NB: Plot the sine and cosine projections normalized by the maximum coefficient. Normalizing by the mean can
; lead to spurious errors since we usually look at anomaly signals where the mean is close to zero.
; NB: Multiply the sine projection by -1 to cancel the i*i contribution in the inverse transform
data_cos=float(data_frq)
data_sin=-imaginary(data_frq)
nrm=max(abs([data_cos,data_sin]))
data_cos=data_cos/nrm
data_sin=data_sin/nrm
plot,frq,data_cos,xstyle=0,xtick_get=x_tck_crd,yrange=[-1,1],ystyle=0,psym=-7,linestyle=0,ytit='Normalized Projection',xtit='Frequency ('+frq_sng+')',charsize=chr_sz
oplot,frq,data_sin,psym=-6,linestyle=1
x_tck_nbr=n_elements(x_tck_crd)
x_sng=strarr(x_tck_nbr)
for idx_tck=0,x_tck_nbr-1 do begin
	period=x_tck_crd(idx_tck)
	if period eq 0. then x_sng(idx_tck)='!9$!5' else begin
		period=1./period ; invert frequency to get period
		x_sng(idx_tck)=auto_sng(period,2)
	endelse
endfor
axis,xaxis=1,xticks=x_tck_nbr-1,xtickn=x_sng,xtitle='!5Period !7s!5 ('+unit_sng+')',charsize=chr_sz
x_scl=!x.crange(1)-!x.crange(0)
y_scl=!y.crange(1)-!y.crange(0)
ln_lgn_x1=!x.crange(0)+0.05*x_scl
ln_lgn_dx=0.07*x_scl
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01*x_scl
txt_lgn_sz=chr_sz
lgn_y_top=!y.crange(0)+.8*y_scl
lgn_dy=0.1*y_scl
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy]
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013*y_scl,linestyle=0,thick=1,/DATA
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013*y_scl,linestyle=1,thick=1,/DATA
xyouts,txt_lgn_x,lgn_y(0),'!5Cosine Projection',size=txt_lgn_sz,/DATA
xyouts,txt_lgn_x,lgn_y(1),'!5Sine Projection',size=txt_lgn_sz,/DATA
endif; endif 0

;if not prn then begin
;print,' Hit any key to continue...'
;junk = get_kbrd(1)
;endif

end; end gcm_spc()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM SPC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


