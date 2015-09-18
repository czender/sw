; $Id$

@/home/zender/idl/ibp_clr.pro

function eqm_vap_ice,temperature

; return the equilibrium water vapor pressure over bulk solid ice, in mbar, 
; for the given temperature in Kelvin
; REMEMBER: value is returned in Pa NOT mb

foo=24.29-6148./temperature
return,100.*exp(foo)		;mb -> Pa
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV lapse rate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_lapse, $
	info=info, $
	prn=prn
if n_elements(info) eq 0 then info=0
if n_elements(prn) eq 0 then prn=0

fl_nm=['/data/zender/_aux0_/cld/cld_10_10_100.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)
nc_id=ncdf_open(fl_nm(fl_idx))
ncdf_varget,nc_id,'orig_env_temp',orig_env_temp
if (size(orig_env_temp))(0) gt 0 then orig_env_temp=reform(orig_env_temp)
ncdf_varget,nc_id,'altitude',altitude
if (size(altitude))(0) gt 0 then altitude=reform(altitude)
ncdf_varget,nc_id,'dz',dz
if (size(dz))(0) gt 0 then dz=reform(dz)
ncdf_close,nc_id
endfor; end loop over files

print,orig_env_temp
lev_nbr=n_elements(altitude)
lapse_rate=-(shift(orig_env_temp,-1)-orig_env_temp)/(dz/1000.)
lapse_rate(0)=0
lapse_rate(lev_nbr-1)=0
ord=altitude/1000.

plot, $
	lapse_rate, $
	ord, $
	xrange=[0,10], $
;	yrange=[-1.5,-1], $
	ystyle=1, $
;	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=0
oplot, $
	lapse_rate-6.5, $
	altitude, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx


if info then begin
endif; info

end; end anv_lapse()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV lapse rate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV C5 Updraft SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_w_ss_c5_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

if prn then open_ps,fl_nm='/data/zender/ps/srcm_w_ss_c5.eps',x_sz=6.5,y_sz=3.2,/eps
anv_w_ss_c5,info=info,prn=prn
if prn then close_ps,fl_nm='/data/zender/ps/srcm_w_ss_c5.eps'
end; end anv_w_ss_c5_bch()
pro anv_w_ss_c5, $
	info=info, $
	prn=prn
if n_elements(info) eq 0 then info=0
if n_elements(prn) eq 0 then prn=0

;fl_nm=['/data/zender/_aux0_/cld/cld_10_w_100_ss.nc']
fl_nm=['/data/zender/_aux0_/cld/cld_05_w_100_ss.nc','/data/zender/_aux0_/cld/cld_10_w_100_ss.nc','/data/zender/_aux0_/cld/cld_15_w_100_ss.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)
nc_id=ncdf_open(fl_nm(fl_idx))
ncdf_varget,nc_id,'IWP_srt',IWP_srt
if (size(IWP_srt))(0) gt 0 then IWP_srt=reform(IWP_srt)
ncdf_varget,nc_id,'IWP_end',IWP_end
IWP_end=reform(IWP_end)
ncdf_varget,nc_id,'density',density
if (size(density))(0) gt 0 then density=reform(density)
ncdf_varget,nc_id,'wind_speed',wind_speed
if (size(wind_speed))(0) gt 0 then wind_speed=reform(wind_speed)
ncdf_varget,nc_id,'temperature_end',temperature_end
if (size(temperature_end))(0) gt 0 then temperature_end=reform(temperature_end)
ncdf_varget,nc_id,'temperature_srt',temperature_srt
if (size(temperature_srt))(0) gt 0 then temperature_srt=reform(temperature_srt)
ncdf_varget,nc_id,'mmr_vapor_srt',mmr_vapor_srt
if (size(mmr_vapor_srt))(0) gt 0 then mmr_vapor_srt=reform(mmr_vapor_srt)
ncdf_varget,nc_id,'mmr_vapor_end',mmr_vapor_end
if (size(mmr_vapor_end))(0) gt 0 then mmr_vapor_end=reform(mmr_vapor_end)
ncdf_varget,nc_id,'IWC_srt',IWC_srt
if (size(IWC_srt))(0) gt 0 then IWC_srt=reform(IWC_srt)
ncdf_varget,nc_id,'IWC_end',IWC_end
if (size(IWC_end))(0) gt 0 then IWC_end=reform(IWC_end)
ncdf_varget,nc_id,'pp_vapor',pp_vapor
if (size(pp_vapor))(0) gt 0 then pp_vapor=reform(pp_vapor)
ncdf_varget,nc_id,'saturation_ice',saturation_ice
if (size(saturation_ice))(0) gt 0 then saturation_ice=reform(saturation_ice)
ncdf_varget,nc_id,'vapor_density',vapor_density
if (size(vapor_density))(0) gt 0 then vapor_density=reform(vapor_density)
ncdf_varget,nc_id,'pressure',pressure
if (size(pressure))(0) gt 0 then pressure=reform(pressure)
dim_id=ncdf_dimid(nc_id,'record')
ncdf_diminq,nc_id,dim_id,foo,nbr_ens
ncdf_close,nc_id

; IWP tendency in SI (kg/m2/s)
diwpdt_si=(IWP_end-IWP_srt)/(60*60)
diwpdt_20m=1000.*diwpdt_si*(20*60)
diwpdt_60m=1000.*diwpdt_si*(60*60)
dqvdt_si=(mmr_vapor_end-mmr_vapor_srt)/(60*60)
diwcdt_si=(IWC_end-IWC_srt)/(60*60)
dqidt_si=(IWC_end-IWC_srt)/(60*60*density)

; Compute dqvidz
gas_const_H2O=461.65 ; J/kg/K 
gas_const_dry_air=287.05 ; J/kg/K
spec_heat_dry_air=1005. ; J/kg/K
mean_sfc_gravity=9.80665 ; m/s2
mmw_H2O=18.015e-3 ; kg/mole
mmw_dry_air=28.9644e-3 ; kg/mole
epsilon=mmw_H2O/mmw_dry_air ; .622
latent_heat_sub=2.834e6 ; J/kg

; cz4 p. 11
;temperature=temperature_end
temperature=temperature_srt
evi=eqm_vap_ice(temperature)
qvi=epsilon*evi/(pressure-evi)
numerator=(1./gas_const_dry_air)-(latent_heat_sub/(spec_heat_dry_air*gas_const_H2O*temperature))
denominator=1.+((qvi*latent_heat_sub*latent_heat_sub)/(spec_heat_dry_air*gas_const_H2O*temperature*temperature))
dqvidz=qvi*mean_sfc_gravity*numerator/(temperature*denominator)

print,'wind_speed (cm/s) = ',wind_speed*100.
print,'delta temperature (K) = ',temperature_end-temperature_srt
print,'evi (mb) = ',evi/100.
print,'qvi (g/kg) = ',qvi*1000.
print,'dqvidz (g/kg/km) = ',1.0e6*dqvidz
print,'dqidt (g/kg/(20 min)) = ',1.2e6*dqidt_si
print,'dqvdt (g/kg/(20 min)) = ',1.2e6*dqvdt_si
print,'diwcdt (g/m3/(20 min)) = ',1.2e6*diwcdt_si
print,'diwpdt (g/m2/(20 min)) = ',1.2e6*diwpdt_si

abc=wind_speed*100.
;ord=diwpdt_si/(wind_speed*density*2000.*dqvidz)
wind_speed(where(wind_speed eq 0))=1.0e-5
ord=dqidt_si/(wind_speed*dqvidz)
ord(0)=ord(1)-(ord(2)-ord(1))/(abc(2)-abc(1))

chr_sz=2.

if fl_idx eq 0 then begin
if prn then ttl='' else ttl='!5Updraft Sens. of !8c!5!I5!N'
x_ttl='!5Updraft !8w!5 (!5cm s!E-1!N)'
if prn then y_mrg=[3,.5] else y_mrg=[3,2] ;[4,2] is [bottom,top] default
if prn then x_mrg=[9,2] else x_mrg=[9,2] ; [10,3] is [left,right] default
plot, $
	abc, $
	ord, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit='!8c!5!I5!N', $
	xmargin=x_mrg, $
	ymargin=y_mrg, $
	xstyle=1, $
	yrange=[-1.25,-.75], $
;	yrange=[-1.5,-1], $
	ystyle=1, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=0
endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx
endelse; fl_idx != 0
endfor; end loop over files

if info then begin

ln_lgn_x1=.40
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.7
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5Height:',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!55 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!510 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!515 km',size=txt_lgn_sz,/NORMAL

endif; info

end; end anv_w_ss_c5()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV C5 Updraft SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV C5 Temperature SS
; Plot the results of the HPMM Temperature sensitivity study
; Uses Temperature_ss in /home/zender/anv/anv_ss.pl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_t_ss_c5_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

if prn then open_ps,fl_nm='/data/zender/ps/srcm_t_ss_c5.eps',x_sz=6.5,y_sz=3.2,/eps
anv_t_ss_c5,info=info,prn=prn
if prn then close_ps,fl_nm='/data/zender/ps/srcm_t_ss_c5.eps'
end; end anv_t_ss_c5_bch()
pro anv_t_ss_c5, $
	info=info, $
	prn=prn
if n_elements(info) eq 0 then info=0
if n_elements(prn) eq 0 then prn=0

; See czpIII p.136

;fl_nm=['/data/zender/_aux0_/cld/cld_z_10_100_ss.nc']
fl_nm=['/data/zender/_aux0_/cld/cld_z_02_100_ss.nc','/data/zender/_aux0_/cld/cld_z_04_100_ss.nc','/data/zender/_aux0_/cld/cld_z_06_100_ss.nc','/data/zender/_aux0_/cld/cld_z_08_100_ss.nc','/data/zender/_aux0_/cld/cld_z_10_100_ss.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)
nc_id=ncdf_open(fl_nm(fl_idx))
ncdf_varget,nc_id,'IWP_srt',IWP_srt
if (size(IWP_srt))(0) gt 0 then IWP_srt=reform(IWP_srt)
ncdf_varget,nc_id,'IWP_end',IWP_end
IWP_end=reform(IWP_end)
ncdf_varget,nc_id,'density',density
if (size(density))(0) gt 0 then density=reform(density)
ncdf_varget,nc_id,'wind_speed',wind_speed
if (size(wind_speed))(0) gt 0 then wind_speed=reform(wind_speed)
ncdf_varget,nc_id,'temperature_end',temperature_end
if (size(temperature_end))(0) gt 0 then temperature_end=reform(temperature_end)
ncdf_varget,nc_id,'temperature_srt',temperature_srt
if (size(temperature_srt))(0) gt 0 then temperature_srt=reform(temperature_srt)
ncdf_varget,nc_id,'mmr_vapor_srt',mmr_vapor_srt
if (size(mmr_vapor_srt))(0) gt 0 then mmr_vapor_srt=reform(mmr_vapor_srt)
ncdf_varget,nc_id,'mmr_vapor_end',mmr_vapor_end
if (size(mmr_vapor_end))(0) gt 0 then mmr_vapor_end=reform(mmr_vapor_end)
ncdf_varget,nc_id,'IWC_srt',IWC_srt
if (size(IWC_srt))(0) gt 0 then IWC_srt=reform(IWC_srt)
ncdf_varget,nc_id,'IWC_end',IWC_end
if (size(IWC_end))(0) gt 0 then IWC_end=reform(IWC_end)
ncdf_varget,nc_id,'pp_vapor',pp_vapor
if (size(pp_vapor))(0) gt 0 then pp_vapor=reform(pp_vapor)
ncdf_varget,nc_id,'saturation_ice',saturation_ice
if (size(saturation_ice))(0) gt 0 then saturation_ice=reform(saturation_ice)
ncdf_varget,nc_id,'vapor_density',vapor_density
if (size(vapor_density))(0) gt 0 then vapor_density=reform(vapor_density)
ncdf_varget,nc_id,'pressure',pressure
if (size(pressure))(0) gt 0 then pressure=reform(pressure)
dim_id=ncdf_dimid(nc_id,'record')
ncdf_diminq,nc_id,dim_id,foo,nbr_ens
ncdf_close,nc_id

; IWP tendency in SI (kg/m2/s)
diwpdt_si=(IWP_end-IWP_srt)/(60*60)
diwpdt_20m=1000.*diwpdt_si*(20*60)
diwpdt_60m=1000.*diwpdt_si*(60*60)
dqvdt_si=(mmr_vapor_end-mmr_vapor_srt)/(60*60)
diwcdt_si=(IWC_end-IWC_srt)/(60*60)
dqidt_si=(IWC_end-IWC_srt)/(60*60*density)

; Compute dqvidz
gas_const_H2O=461.65 ; J/kg/K 
gas_const_dry_air=287.05 ; J/kg/K
spec_heat_dry_air=1005. ; J/kg/K
mean_sfc_gravity=9.80665 ; m/s2
mmw_H2O=18.015e-3 ; kg/mole
mmw_dry_air=28.9644e-3 ; kg/mole
epsilon=mmw_H2O/mmw_dry_air ; .622
latent_heat_sub=2.834e6 ; J/kg

; cz4 p. 11
;temperature=temperature_end
temperature=temperature_srt
evi=eqm_vap_ice(temperature)
qvi=epsilon*evi/(pressure-evi)
numerator=(1./gas_const_dry_air)-(latent_heat_sub/(spec_heat_dry_air*gas_const_H2O*temperature))
denominator=1.+((qvi*latent_heat_sub*latent_heat_sub)/(spec_heat_dry_air*gas_const_H2O*temperature*temperature))
dqvidz=qvi*mean_sfc_gravity*numerator/(temperature*denominator)

print,'wind_speed (cm/s) = ',wind_speed*100.
print,'delta temperature (K) = ',temperature_end-temperature_srt
print,'evi (mb) = ',evi/100.
print,'qvi (g/kg) = ',qvi*1000.
print,'dqvidz (g/kg/km) = ',1.0e6*dqvidz
print,'dqidt (g/kg/(20 min)) = ',1.2e6*dqidt_si
print,'dqvdt (g/kg/(20 min)) = ',1.2e6*dqvdt_si
print,'diwcdt (g/m3/(20 min)) = ',1.2e6*diwcdt_si
print,'diwpdt (g/m2/(20 min)) = ',1.2e6*diwpdt_si

abc=temperature
;ord=diwpdt_si/(wind_speed*density*2000.*dqvidz)
idx_still=where(wind_speed eq 0,cnt)
if cnt ne 0 then wind_speed(idx_still)=1.0e-5
ord=dqidt_si/(wind_speed*dqvidz)
ord(0)=ord(1)-(ord(2)-ord(1))/(abc(2)-abc(1))

chr_sz=2.

if fl_idx eq 0 then begin
if prn then ttl='' else ttl='!5Temperature Sens. of !8c!5!I5!N'
x_ttl='!5Mid-Cloud Temperature !8T!5 (!5!E!12_!5!NK)'
if prn then y_mrg=[3,.5] else y_mrg=[3,2] ;[4,2] is [bottom,top] default
if prn then x_mrg=[9,2] else x_mrg=[9,2] ; [10,3] is [left,right] default
plot, $
	abc, $
	ord, $
	tit=ttl, $
	xtit=x_ttl, $
	ytit='!8c!5!I5!N', $
	xmargin=x_mrg, $
	ymargin=y_mrg, $
	xstyle=1, $
	yrange=[-1.25,-.75], $
;	yrange=[-1.5,-1], $
	ystyle=1, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=0
endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx
endelse; fl_idx != 0
endfor; end loop over files

if info then begin

ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.7
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(4)+0.013,linestyle=4,thick=3.0,/NORMAL

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5Updraft !8w!5 (cm/s):',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!52',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!54',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!56',size=txt_lagend_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!58',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(4),'!510',size=txt_lgn_sz,/NORMAL

endif; info

end; end anv_t_ss_c5()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV C5 Temperature SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV Temperature SS
; Plot the results of the HPMM Temperature sensitivity study
; Uses Temperature_ss in /home/zender/anv/anv_ss.pl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_t_ss_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0
if prn then open_ps,fl_nm='/data/zender/ps/srcm_t_ss.eps',x_sz=6.5,y_sz=3.2,/eps
anv_t_ss,prn=prn
if prn then close_ps,fl_nm='/data/zender/ps/srcm_t_ss.eps'
end; end anv_t_ss_bch()
pro anv_t_ss, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

; See czpIII p.136

;fl_nm=['/data/zender/_aux0_/cld/cld_z_10_100_ss.nc']
fl_nm=['/data/zender/_aux0_/cld/cld_z_02_100_ss.nc','/data/zender/_aux0_/cld/cld_z_04_100_ss.nc','/data/zender/_aux0_/cld/cld_z_06_100_ss.nc','/data/zender/_aux0_/cld/cld_z_08_100_ss.nc','/data/zender/_aux0_/cld/cld_z_10_100_ss.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

print,'Processing '+fl_nm(fl_idx)
nc_id=ncdf_open(fl_nm(fl_idx))
;ncdf_varget,nc_id,'IWP_srt',IWP_srt
;if (size(IWP_srt))(0) gt 0 then IWP_srt=reform(IWP_srt)
; dbg dbg dbg
IWP_srt=0.0733108
ncdf_varget,nc_id,'IWP_end',IWP_end
IWP_end=reform(IWP_end)
ncdf_varget,nc_id,'temperature_end',temperature_end
if (size(temperature_end))(0) gt 0 then temperature_end=reform(temperature_end)
ncdf_varget,nc_id,'temperature_srt',temperature_srt
if (size(temperature_srt))(0) gt 0 then temperature_srt=reform(temperature_srt)
ncdf_varget,nc_id,'wind_speed',wind_speed
dim_id=ncdf_dimid(nc_id,'record')
ncdf_diminq,nc_id,dim_id,foo,nbr_ens
ncdf_close,nc_id

; IWP tendency in SI (kg/m2/s)
diwpdt_si=(IWP_end-IWP_srt)/(60*60)
diwpdt_20m=1000.*diwpdt_si*(20*60)
diwpdt_60m=1000.*diwpdt_si*(60*60)

temperature=temperature_end
;temperature=temperature_srt
abc=temperature
print,'temperature = ',temperature
ord=diwpdt_20m

chr_sz=2.

if fl_idx eq 0 then begin
ttl='!5Temperature Sensitivity'
if prn then x_ttl='' else x_ttl='!5Mid-Cloud Temperature !8T!5 (!5!E!12_!5!NK)'
if prn then y_mrg=[.5,2] else y_mrg=[3,2] ;[4,2] is [bottom,top] default
if prn then x_mrg=[9,2] else x_mrg=[9,2] ; [10,3] is [left,right] default
plot, $
	abc, $
	ord, $
	tit=ttl, $
	xtit=x_ttl, $
;	ytit='!9D!I!8t!NIWP!5 (g m!E-2!N (20 min.)!E-1!N)', $
	ytit='!9D!I!8t!NIWP!5', $
	xmargin=x_mrg, $
	ymargin=y_mrg, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	xstyle=1, $
	yrange=[0,100], $
	ystyle=1, $
	linestyle=0
endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx

endelse; fl_idx != 0
endfor; end loop over files

ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.6
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(4)+0.013,linestyle=4,thick=3.0,/NORMAL

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5Updraft !8w!5 (cm/s):',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!52',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!54',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!56',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!58',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(4),'!510',size=txt_lgn_sz,/NORMAL

end; end anv_t_ss()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV Temperature SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV Updraft SS
; Plot the results of the HPMM Updraft sensitivity study
; Uses updraft_ss in /home/zender/anv/anv_ss.pl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_w_ss_bch, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

if prn then open_ps,fl_nm='/data/zender/ps/srcm_w_ss.eps',x_sz=6.5,y_sz=3.2,/eps
anv_w_ss,prn=prn
if prn then close_ps,fl_nm='/data/zender/ps/srcm_w_ss.eps'
end; end anv_w_ss_bch()
pro anv_w_ss, $
	prn=prn
if n_elements(prn) eq 0 then prn=0

;fl_nm=['/data/zender/_aux0_/cld/cld_10_w_100_ss.nc']
fl_nm=['/data/zender/_aux0_/cld/cld_05_w_100_ss.nc','/data/zender/_aux0_/cld/cld_10_w_100_ss.nc','/data/zender/_aux0_/cld/cld_15_w_100_ss.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

nc_id=ncdf_open(fl_nm(fl_idx))
ncdf_varget,nc_id,'IWP_srt',IWP_srt
if (size(IWP_srt))(0) gt 0 then IWP_srt=reform(IWP_srt)
ncdf_varget,nc_id,'IWP_end',IWP_end
IWP_end=reform(IWP_end)
ncdf_varget,nc_id,'density',density
if (size(density))(0) gt 0 then density=reform(density)
ncdf_varget,nc_id,'wind_speed',wind_speed
dim_id=ncdf_dimid(nc_id,'record')
ncdf_diminq,nc_id,dim_id,foo,nbr_ens
ncdf_close,nc_id

; IWP tendency in SI (kg/m2/s)
diwpdt_si=(IWP_end-IWP_srt)/(60*60)
diwpdt_20m=1000.*diwpdt_si*(20*60)
diwpdt_60m=1000.*diwpdt_si*(60*60)

abc=wind_speed*100.
print,'density = ',density
;mid-cld vapor_density for 2 km thick cld at 5,10,15 km 100% RH_ice is (kg/m3): 
vapor_density=[.00255268,.000131759,3.82737e-06]

ord=diwpdt_20m
;ord=diwpdt_20m/(vapor_density(fl_idx)/density)

chr_sz=2.

if fl_idx eq 0 then begin
ttl='!5Updraft Sensitivity'
if prn then x_ttl='' else x_ttl='!5Updraft !8w!5 (!5cm s!E-1!N)'
if prn then y_mrg=[.5,2] else y_mrg=[3,2] ;[4,2] is [bottom,top] default
if prn then x_mrg=[9,2] else x_mrg=[9,2] ; [10,3] is [left,right] default
plot, $
	abc, $
	ord, $
	tit=ttl, $
	xtit=x_ttl, $
;	ytit='!9D!I!8t!NIWP!5 (g m!E-2!N (20 min.)!E-1!N)', $
	ytit='!9D!I!8t!NIWP!5', $
	xmargin=x_mrg, $
	ymargin=y_mrg, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	xstyle=1, $
	yrange=[0,100], $
	ystyle=1, $
	linestyle=0
endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx

endelse; fl_idx != 0
endfor; end loop over files

ln_lgn_x1=.40
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.6
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5Height:',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!55 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!510 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!515 km',size=txt_lgn_sz,/NORMAL

end; end anv_w_ss()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV Updraft SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV Sedimentation SS
; Plot the results of the HPMM Sedimentation sensitivity study
; Uses decay_ss in /home/zender/anv/anv_ss.pl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_decay_ss_bch
open_ps,fl_nm='/data/zender/ps/srcm_decay_ss.eps',x_sz=6.5,y_sz=3.2,/eps
anv_decay_ss
close_ps,fl_nm='/data/zender/ps/srcm_decay_ss.eps'
end; end anv_decay_ss_bch()
pro anv_decay_ss

; See czpIII p.136

fl_nm=['/data/zender/_aux0_/cld/cld_05_00_100_decay_ss.nc','/data/zender/_aux0_/cld/cld_10_00_100_decay_ss.nc','/data/zender/_aux0_/cld/cld_15_00_100_decay_ss.nc']
base_hgt=[5000,10000,15000]
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin

nc_id=ncdf_open(fl_nm(fl_idx))
ncdf_varget,nc_id,'IWP_cld',IWP_cld
if (size(IWP_cld))(0) gt 0 then IWP_cld=reform(IWP_cld)
dim_id=ncdf_dimid(nc_id,'record')
ncdf_diminq,nc_id,dim_id,foo,nbr_ens
ncdf_close,nc_id

abc=indgen(nbr_ens)*600./3600.	; snapshots every 10 min. converted to SI then to hours for labels
ord=IWP_cld*1000.

chr_sz=2.

if fl_idx eq 0 then begin
plot, $
	abc, $
	ord, $
	charsize=chr_sz, $
	max_value=1.0e20, $
	thick=2.0, $
	tit='!5Sedimentation Sens. Study', $
	xmargin=[7,2], $ ;[10,3] is [left,right] default
	xstyle=1, $
	xtit='!5Time !8t!5 (hr.)', $
	ymargin=[3,2], $ ;[4,2] is [bottom,top] default
	yrange=[0,100], $
	ystyle=1, $
	ytit='!8IWP!5 (g m!E-2!N)', $
	linestyle=fl_idx+1

endif else begin; fl_idx != 0
oplot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=fl_idx+1

endelse; fl_idx != 0
endfor; end loop over files

tau_1=ord(0)*exp(-abc/.5)
tau_2=ord(0)*exp(-abc/2.)

oplot, $
	abc, $
	tau_1, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=0

oplot, $
	abc, $
	tau_2, $
	max_value=1.0e20, $
	thick=2.0, $
	linestyle=0

ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.7
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=3,thick=3.0,/NORMAL

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5Height:',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'!55 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!510 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!515 km',size=txt_lgn_sz,/NORMAL

plots,[ln_lgn_x1,ln_lgn_x2]-.4,lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2]-.4,lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
xyouts,ln_lgn_x1-.4,lgn_y(0)+lgn_dy,'!5Timescale !7s!5:',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x-.4,lgn_y(0),'.5 hr (lower)',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x-.4,lgn_y(1),'2 hr (upper)',size=txt_lgn_sz,/NORMAL

end; end anv_decay_ss()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV Sedimentation SS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure SPCP Precip
; Plot the time series domain average precipitation from the CEM over 
; the series for the GCM.

; ncecat -O -D 3 -d lon,120.0,130. -d lat,10.0,20. -v PRECT,PRECC,ORO omega.??.nc gcm_0114.nc
; ncrename -d record,foo gcm_0114.nc foo.nc
; ncecat -O -D 3 foo.nc foo.nc
; ncwa -O -D 3 -m ORO -M 0. -o eq -a foo foo.nc gcm_ocn_0114.nc
; ncwa -O -D 3 -o gt -m PRECC -M 1.736e-7 -a lon,lat gcm_ocn_0114.nc gcm_ocn_xyavg_0114.nc
; ncwa -a record,time gcm_ocn_xyavg_0114.nc gcm_ocn_txyavg_0114.nc
; ncrename -d record,day -v PRECT,gcm_precip gcm_ocn_xyavg_0114.nc
; ncks -H -C -v PRECT gcm_ocn_txyavg_0114.nc

; ncrename -d record,hour -v dom_avg_precip,cem_precip cem_0424_precip.nc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro spcp_precip_bch
open_ps,fl_nm='/data/zender/ps/precip.eps',/one_clm,/eps
spcp_precip
close_ps,fl_nm='/data/zender/ps/precip.eps'
end; end spcp_precip_bch()
pro spcp_precip

nc_id=ncdf_open('/data/zender/_aux0_/omega0_2/gcm_ocn_xyavg_0114.nc')
ncdf_varget,nc_id,'PRECT',gcm_precip
dim_id=ncdf_dimid(nc_id,'day')
ncdf_diminq,nc_id,dim_id,foo,nbr_day
dim_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,dim_id,foo,time_nbr
ncdf_close,nc_id
gcm_precip=reform(gcm_precip,nbr_day*time_nbr)
gcm_precip=gcm_precip*8.64e7
gcm_day=1.+nbr_day*findgen(nbr_day*time_nbr)/(nbr_day*time_nbr)

;nc_id=ncdf_open('/data/zender/_aux0_/cem/dom_avg.nc')
;ncdf_varget,nc_id,'PRECT',gcm_precip
;ncdf_varget,nc_id,'time',gcm_day
;ncdf_close,nc_id
;gcm_precip=reform(gcm_precip)
;gcm_day=gcm_day-gcm_day(0)+1

nc_id=ncdf_open('/data/zender/_aux0_/cem/cem_0424_precip.nc')
ncdf_varget,nc_id,'cem_precip',cem_precip
ncdf_close,nc_id
cem_precip=reform(cem_precip)
cem_precip=cem_precip*8.64e7
cem_day=4.+21*findgen(21*24)/(21*24)

print,'average gcm precip = ',total(gcm_precip)/n_elements(gcm_precip)
print,'average cem precip = ',total(cem_precip)/n_elements(cem_precip)

abc=cem_day
ord=cem_precip

chr_sz=2.
plot, $
	abc, $
	ord, $
	max_value=1.0e20, $
	tit='!5Precipitation Timeseries', $
	xtit='!5Model Time !8t!5 (days)', $
	ytit='!5Precipitation !8P!5 (mm day!E-1!N)', $
	xstyle=1, $
	ystyle=1, $
	xrange=[0.0,25.], $
	yrange=[0.0,80.0], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[5.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	linestyle=0

abc=gcm_day
ord=gcm_precip
oplot,	$
	abc, $
	ord, $
	thick=3.0, $
	linestyle=2

ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL

xyouts,txt_lgn_x,lgn_y(0),'!5CEM',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5GCM',size=txt_lgn_sz,/NORMAL

end; end spcp_precip()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure SPCP Precip
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure GCM forcing
; Plot the time-domain average dynamic environment from the gcm

; Mask out weak PRECC points and land ORO
; ncwa -v ORO,U,W,gw -m PRECC -M 1.74e-7 -o gt -a time -d lon,120.0,130. -d lat,10.0,20. omega2_0114avg.nc foo.nc
; ncwa -v U,W,gw -m ORO -M 0. -o eq -a lon foo.nc foo.nc
; ncwa -v U,W -w gw -a lat foo.nc omega2_0114avg_TWPavg.nc

; Mask out weak PRECC points
; ncwa -v U,W,gw -m PRECC -M 1.74e-7 -o gt -a time -d lon,120.0,130. -d lat,10.0,20. omega2_0114avg.nc foo.nc
; ncwa -v U,W -w gw -a lat,lon foo.nc omega2_0114avg_TWPavg.nc

; Mask out land ORO
; ncwa -v U,W -m ORO -M 0. -o eq -a time,lat,lon -d lon,120.0,130. -d lat,10.0,20. omega2_0114avg.nc omega2_0114avg_TWPavg.nc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gcm_forcing_bch
open_ps,fl_nm='/data/zender/ps/gcm_forcing.eps',/one_clm,/eps
gcm_forcing,fl_nm='/data/zender/_aux0_/omega0_2/omega2_0114avg_TWPavg.nc',rng_y=20.
close_ps,fl_nm='/data/zender/ps/gcm_forcing.eps'
end; end gcm_forcing_bch()
pro gcm_forcing, $
	rng_y=rng_y, $
	fl_nm=fl_nm

if n_elements(rng_y) eq 0 then rng_y=20.
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/omega0_2/omega2_0114avg_TWPavg.nc'

; Get the requested field

nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'U',u
ncdf_varget,nc_id,'W',w
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id

lev_nbr=n_elements(lev)

u=reform(u)
w=reform(w)

abc_1=u
abc_2=w*100.

;ord=lev
;ord=[4.8093, 13.0731, 32.5591, 63.9471, 99.0432, 138.713, 189.191, 251.239, 324.848, 408.955, 501.275, 598.248, 695.169, 786.51, 866.407, 929.275, 970.446, 992.528]
ord=[35,30,23,19,16.7,14.9,12.8,11.1,9.25,7.5,6.,4.5,3.3,2.3,1.35,0.75,.35,.18]
print,abc_1
print,abc_2

chr_sz=2.
plot, $
	abc_1, $
	ord, $
	xtit='!5Mean Zonal Wind !8u!5 (m s!E-1!N)', $
	ytit='!5Altitude !8z!5 (km)', $
	xstyle=8, $
	ystyle=1, $
	xrange=[-20.0,20.0], $
	yrange=[0.0,rng_y], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[5.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2.5], $  ;[4,2] is [bottom,top] default
	linestyle=0

plot, $
	abc_2, $
	ord, $
	xstyle=12, $
	ystyle=12, $
	xrange=[-40.0,40.0], $
	yrange=[0.0,rng_y], $
	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $
	linestyle=1, $
;	position=plt_rgn_nrm, $
	xmargin=[5.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2.5], $  ;[4,2] is [bottom,top] default
	/noerase

axis, $
	xaxis=1, $
	xtitle='!5Mean Vertical Updraft !8w!5 (cm s!E-1!N)', $
	charsize=chr_sz

ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.6
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL

xyouts,txt_lgn_x,lgn_y(0),'!8u!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8w!5',size=txt_lgn_sz,/NORMAL

end; end gcm_forcing()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure GCM forcing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ANV FICE
; Plot ice fraction from numbers in a table
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anv_fice,prn=prn,dbg=dbg
; anv_fice,prn=1

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0

; Some initialization
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.2
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

pressure=[787,695,598,501,409,325,251]
temperature=[15,9,2,-5,-15,-27,-42]
fice_ccm=[0,0,0,0,0.24,0.84,1.0]
fice_anv=[0,0,.12,0.98,1.0,1.0,1.0]
fice_cem=[0,.04,0.25,0.91,1.0,1.0,1.0]

abc=fice_cem
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=temperature
ord_min=min(ord)
ord_max=max(ord)
rng_y=[ord_max,ord_min]
chr_sz=2.
ttl='Equatorial Pacific Ice Fraction'
x_ttl='!5Ice Fraction (%)'
y_ttl='!5Temperature (!5!E!12_!N!5C)'

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
	fl_nm_out='/data/zender/ps/anv_fice.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
endif; endif prn

; xyplot ice fraction
!p.multi=0
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=0.5,charsize=chr_sz,linestyle=0

oplot,fice_ccm,ord,thick=2.0,linestyle=1
oplot,fice_anv,ord,thick=2.0,linestyle=2

t_500mb=-6.
arrow,!x.crange(0)+1.05*(!x.crange(1)-!x.crange(0)),t_500mb,!x.crange(1),t_500mb,thick=2.0,/DATA
arrow,!x.crange(0)-0.05*(!x.crange(1)-!x.crange(0)),t_500mb,!x.crange(0),t_500mb,thick=2.0,/DATA

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end anv_fice()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ANV FICE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
