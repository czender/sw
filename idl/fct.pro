;
; $Id$
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-15 02:07:50  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:49  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:18  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:47  zender
; Imported sources
;
; Revision 1.2  1992/07/29  19:57:32  zender
; still not working, this version implements all bc's as per boris,
; why does not the cloud move?
;
; Revision 1.1  1992/07/27  15:50:33  zender
; Initial revision
;
pro fct,num_layer=num_layer,courant_safety=courant_safety, $
	nsteps=nsteps, $
	pltstep=pltstep,pause=pause,$
	bpd=bpd
;
; Does straight upstream advection on the vertical advection equation.
; Is built into smolar.pro, where anti-diffusive steps are added.
;
; Example usage 
;
;fct,pause='n',nsteps=61,pltstep=60
;fct,pause='n',nsteps=601,pltstep=600
;fct,pause='n',nsteps=201,pltstep=50,num_layer=50
;fct,pause='n',nsteps=201,pltstep=50,num_layer=50,bpd=10
sys_time = systime(1) 
if n_elements(num_layer) eq 0 then num_layer = 50
if n_elements(nsteps) eq 0 then nsteps = 60
if n_elements(pltstep) eq 0 then pltstep = 10		;time steps per graph
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if pause eq 'y' then !p.multi=0 else begin
	!p.multi=[0,2,0]
	pause = 'n'
endelse
if n_elements(courant_safety) eq 0 then courant_safety = .2
if n_elements(bpd) eq 0 then bpd = 5.
;
;if n_elements() eq 0 then = 
;
; Initialize Physical Constants
;
env_temperature=223.
cloud_base=10000.	;m
cloud_thick=1000.	;m
bin_sz=20.e-6	;for sz_dst purposes
crystal_length=50.e-6
crystal_mass=mass_of_length(crystal_length)
conc_LBC_add=0.
conc_RBC_add=0.
conc_LBC_mult=0.
conc_RBC_mult=0.
;
time_array=fltarr(nsteps+1)
IWP=fltarr(nsteps+1)
column_density=fltarr(nsteps+1)
zdzt=fltarr(num_layer+2)
zdzt_lft=fltarr(num_layer+2)
zdzt_rgt=fltarr(num_layer+2)
altitude=fltarr(num_layer+2)
altitude_rgt=fltarr(num_layer+2)
altitude_lft=fltarr(num_layer+2)
warp_dz=fltarr(num_layer+2)
warp_dz_rgt=fltarr(num_layer+2)
warp_dz_lft=fltarr(num_layer+2)
concentration=fltarr(num_layer+2)
new_concentration=fltarr(num_layer+2)
conc_rgt=fltarr(num_layer+2)
conc_lft=fltarr(num_layer+2)
epsilon_rgt=fltarr(num_layer+2)
epsilon_lft=fltarr(num_layer+2)
nu_rgt=fltarr(num_layer+2)
nu_lft=fltarr(num_layer+2)
mu_rgt=fltarr(num_layer+2)
mu_lft=fltarr(num_layer+2)
convected_conc=fltarr(num_layer+2)
trans_diff_conc=fltarr(num_layer+2)
raw_ad_flux_rgt=fltarr(num_layer+2)
raw_ad_flux_lft=fltarr(num_layer+2)
cor_ad_flux_rgt=fltarr(num_layer+2)
cor_ad_flux_lft=fltarr(num_layer+2)
S_rgt=fltarr(num_layer+2)
S_lft=fltarr(num_layer+2)
;
bins_per_doubling=float(bpd)
smallest_thick=cloud_thick/(2.^((num_layer-1.+.5)/bins_per_doubling))
;
; Note how _rgt and _lft quantities are undefined and unused
;
for layer=1,num_layer do begin
;
; Eq 8-3.13
;
	altitude_lft(layer)=cloud_base+smallest_thick*2.^$
		((layer-1.-.5)/bins_per_doubling)
;	altitude_lft(layer)=cloud_base+(layer-1)*cloud_thick/num_layer
	altitude_rgt(layer)=cloud_base+smallest_thick*2.^$
		((layer-1.+.5)/bins_per_doubling)
;	altitude_rgt(layer)=cloud_base+layer*cloud_thick/num_layer
        altitude(layer)=0.5*(altitude_rgt(layer)+altitude_lft(layer))
endfor
altitude_lft(1)=cloud_base
altitude(1)=0.5*(altitude_rgt(1)+altitude_lft(1))
;
;altitude(0)=altitude(num_layer)
;altitude(num_layer+1)=altitude(1)
;
for layer=1,num_layer do begin
;
; Eq 8-3.34
	warp_dz(layer)=altitude_rgt(layer)-altitude_lft(layer)
endfor
;warp_dz(0)=warp_dz(num_layer)
;warp_dz(num_layer+1)=warp_dz(1)
;warp_dz(0)=warp_dz(1)
;warp_dz(num_layer+1)=warp_dz(num_layer)
;
; Note how _rgt and _lft quantities are undefined and unused
; This also disagrees with Boris because of cyclic conditions
;
for layer=1,num_layer-1 do begin
;
; Eq 8-3.23
	warp_dz_rgt(layer)=0.5*(warp_dz(layer+1)+warp_dz(layer))
endfor
warp_dz_rgt(num_layer)=warp_dz(num_layer)
;
; Vector method of doing it
;
;warp_dz_rgt(num_layer)=0.5*(warp_dz(num_layer+1)+warp_dz(num_layer))
;warp_dz_lft(2:num_layer)=warp_dz_rgt(1:num_layer-1)
;warp_dz_lft(1)=0.5*(warp_dz(1)+warp_dz(1-1))
for layer=2,num_layer do begin
	warp_dz_lft(layer)=0.5*(warp_dz(layer)+warp_dz(layer-1))
endfor
warp_dz_lft(1)=warp_dz(1)
;
; Eq 8-3.24 This reinstates the Boris values
;
for layer=1,num_layer do begin
endfor

for layer=1,num_layer do begin
	if altitude(layer) lt  10400. then $
		concentration(layer)=0. else $
	if altitude(layer) ge 10500. then $
		concentration(layer)=0. else $
	concentration(layer)=bin_sz*1.0e6* $
		sz_dst_of_IWC(crystal_length,env_temperature)
endfor
;
; Eq 8-3.18 This isn't quite the same scheme as Boris, because i want cyclic B.C.'s
;concentration(0)=concentration(num_layer)
;concentration(num_layer+1)=concentration(1)
concentration(0)=conc_LBC_mult*concentration(1)+conc_LBC_add
concentration(num_layer+1)=conc_RBC_mult*concentration(num_layer)+conc_RBC_add
;
for layer=1,num_layer do begin
;
; Eq 8-3.17
	conc_rgt(layer)=0.5*(concentration(layer+1)+concentration(layer))
	conc_lft(layer)=0.5*(concentration(layer)+concentration(layer-1))
endfor
;
for layer=1,num_layer do begin
	zdzt(layer)=1.
endfor
;zdzt(1:num_layer)=-1-sin(2.*!pi*findgen(num_layer)/(num_layer-1))	;1 m/s
;zdzt(1:num_layer)=-1-findgen(num_layer)/(num_layer-1)	;1 m/s
;
; Note how zdzt_lft(0),zdzt_rgt(num_layer+1) etc are undefined
; but are never used anyway
; Eq 8-3.15 Since the grid velocity is zero (Eulerian), delta v = zdzt
for layer=1,num_layer-1 do begin
;
; Eq 8-3.14
;ie velocity on the right interface of cell # layer
	zdzt_rgt(layer)=0.5*(zdzt(layer)+zdzt(layer+1))
endfor
zdzt_rgt(num_layer)=zdzt(num_layer)
for layer=2,num_layer do begin
;
; Eq 8-3.14
;ie velocity on the left interface of cell # layer
	zdzt_lft(layer)=0.5*(zdzt(layer-1)+zdzt(layer))
endfor	
zdzt_lft(1)=zdzt(1)
;
; Courant stability
dt=courant_safety/max(abs(zdzt(1:num_layer))/warp_dz(1:num_layer))
;
for layer=2,num_layer-1 do begin
;
; Eq 8-3.31
	epsilon_rgt(layer)=zdzt_rgt(layer)*.5*dt*$
	((1./warp_dz(layer))+(1./warp_dz(layer+1)))
	epsilon_lft(layer)=zdzt_lft(layer)*.5*dt*$
	((1./warp_dz(layer-1))+(1./warp_dz(layer)))
endfor
epsilon_rgt(1)=zdzt_rgt(1)*.5*dt*$
((1./warp_dz(1))+(1./warp_dz(1+1)))
epsilon_lft(num_layer)=zdzt_lft(num_layer)*.5*dt*$
((1./warp_dz(num_layer-1))+(1./warp_dz(num_layer)));
;
epsilon_rgt(num_layer)=zdzt_rgt(num_layer)*dt/warp_dz(num_layer)
epsilon_lft(1)=zdzt_lft(1)*dt/warp_dz(1)
;
for layer=1,num_layer do begin
;
; Eq 8-3.30
	nu_rgt(layer)=(1./6.)+(1./3.)*epsilon_rgt(layer)*epsilon_rgt(layer)
	nu_lft(layer)=(1./6.)+(1./3.)*epsilon_lft(layer)*epsilon_lft(layer)
;
	mu_rgt(layer)=(1./6.)-(1./6.)*epsilon_rgt(layer)*epsilon_rgt(layer)
	mu_lft(layer)=(1./6.)-(1./6.)*epsilon_lft(layer)*epsilon_lft(layer)
endfor
;
print,'smallest thick = ',smallest_thick
print,'altitude = ',altitude
print,'altitude_lft = ',altitude_lft
print,'altitude_rgt = ',altitude_rgt
print,'concentration = ',concentration
print,'conc_lft = ',conc_lft
print,'conc_rgt = ',conc_rgt
print,'warp_dz = ',warp_dz
;print,'warp_dz_lft = ',warp_dz_lft
print,'epsilon_rgt = ',epsilon_rgt
print,'epsilon_lft = ',epsilon_lft
print,'mu_rgt = ',mu_rgt
print,'mu_lft = ',mu_lft
print,'nu_rgt = ',nu_rgt
print,'nu_lft = ',nu_lft
print,'total warp_dz(1:nl) = ',total(warp_dz(1:num_layer))
print,'total warp_dz_lft(1:nl) = ',total(warp_dz_lft(1:num_layer))
;
time_array(0)=0.
IWP(0)=total(concentration(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
column_density(0)=total(concentration(1:num_layer)*warp_dz(1:num_layer));#/m^2
;
print,''
print,'FCT Advection Model:'
print,'Integration period = ',nsteps*dt,' seconds'
print,'Stepsize = ',dt,' seconds'
print,'crystal mass (ng) = ',crystal_mass*1.0e12
print,'crystal length (u) = ',crystal_length*1.0e6
print,'printing one "." per time step:'
;
for time_step=1,nsteps do begin
time_array(time_step)=time_array(time_step-1)+dt
print,'.',format='($,A1)'
;print,time_step,format='($,I,1x)'
;
if float(time_step-1)/pltstep eq fix(time_step/pltstep) then begin 
;
plt_fct,time_array=time_array,$
	dt=dt,warp_dz=warp_dz,num_layer=num_layer, $
	time_step=time_step, $
	concentration=concentration, $
	column_density=column_density, $
	crystal_length=crystal_length,pause=pause, $
	crystal_mass=crystal_mass,$
	altitude=altitude,$
	junk=junk,zdzt=zdzt,$
	IWP=IWP
if n_elements(junk) ne 0 then if junk eq 'q' then goto, exit_gracefully
;
endif  ;end if time to plot
;
for layer=1,num_layer do begin
;
; Eq 8-3.20 Finds rho star
	convected_conc(layer)=concentration(layer)-$
	(dt*conc_rgt(layer)*zdzt_rgt(layer)+$
	dt*conc_lft(layer)*zdzt_lft(layer))/warp_dz(layer)
	if(abs(convected_conc(layer)*warp_dz(layer)) le 1.0e-15) then $
		convected_conc(layer)=0.
endfor
;convected_conc(0)=convected_conc(num_layer)
;convected_conc(num_layer+1)=convected_conc(1)
convected_conc(0)=conc_LBC_mult*convected_conc(1)+conc_LBC_add
convected_conc(num_layer+1)=conc_RBC_mult*convected_conc(num_layer)+conc_RBC_add
;
for layer=1,num_layer do begin
;
; Eq 8-3.22 Finds rho twiddle
	trans_diff_conc(layer)=convected_conc(layer)+$
	(nu_rgt(layer)*warp_dz_rgt(layer)*$
	(concentration(layer+1)-concentration(layer))-$
	nu_lft(layer)*warp_dz_lft(layer)*$
	(concentration(layer)-concentration(layer-1)))/warp_dz(layer)
endfor
;trans_diff_conc(0)=trans_diff_conc(num_layer)
;trans_diff_conc(num_layer+1)=trans_diff_conc(1)
trans_diff_conc(0)=conc_LBC_mult*trans_diff_conc(1)+conc_LBC_add
trans_diff_conc(num_layer+1)=conc_RBC_mult*trans_diff_conc(num_layer)+conc_RBC_add
;
for layer=1,num_layer do begin
;
; Eq 8-3.25 Finds f super ad
	raw_ad_flux_rgt(layer)=mu_rgt(layer)*warp_dz_rgt(layer)*$
	(convected_conc(layer+1)-convected_conc(layer))
	raw_ad_flux_lft(layer)=mu_lft(layer)*warp_dz_lft(layer)*$
	(convected_conc(layer)-convected_conc(layer-1))
endfor
;
for layer=1,num_layer do begin
;
; Eq 8-3.9 Finds S
	if(trans_diff_conc(layer+1)-trans_diff_conc(layer)) gt 0. then $
		S_rgt(layer)=1. else S_rgt(layer)=-1.
	if(trans_diff_conc(layer)-trans_diff_conc(layer-1)) gt 0. then $
		S_lft(layer)=1. else S_lft(layer)=-1.
endfor
;
for layer=2,num_layer-1 do begin
;
; Eq 8-3.32 Finds f super c
	min_foo_rgt=min([abs(raw_ad_flux_rgt(layer)),$
		S_rgt(layer)*warp_dz(layer+1)*(trans_diff_conc(layer+2)-$
	trans_diff_conc(layer+1)),$
		S_rgt(layer)*warp_dz(layer)*(trans_diff_conc(layer)-$
	trans_diff_conc(layer-1))])
	max_foo_rgt=max([0.0,min_foo_rgt])
	cor_ad_flux_rgt(layer)=S_rgt(layer)*max_foo_rgt

	min_foo_lft=min([abs(raw_ad_flux_lft(layer)),$
		S_lft(layer)*warp_dz(layer)*(trans_diff_conc(layer+1)-$
	trans_diff_conc(layer)),$
		S_lft(layer)*warp_dz(layer-1)*(trans_diff_conc(layer-1)-$
	trans_diff_conc(layer-2))])
	max_foo_lft=max([0.0,min_foo_lft])
	cor_ad_flux_lft(layer)=S_lft(layer)*max_foo_lft
endfor
;
min_foo_rgt=min([abs(raw_ad_flux_rgt(1)),$
	S_rgt(1)*warp_dz(1+1)*(trans_diff_conc(1+2)-$
trans_diff_conc(1+1)),$
	S_rgt(1)*warp_dz(1)*(trans_diff_conc(1)-$
trans_diff_conc(1-1))])
max_foo_rgt=max([0.0,min_foo_rgt])
cor_ad_flux_rgt(1)=S_rgt(1)*max_foo_rgt
;
min_foo_lft=min([abs(raw_ad_flux_lft(num_layer)),$
	S_lft(num_layer)*warp_dz(num_layer)*(trans_diff_conc(num_layer+1)-$
trans_diff_conc(num_layer)),$
	S_lft(num_layer)*warp_dz(num_layer-1)*(trans_diff_conc(num_layer-1)-$
trans_diff_conc(num_layer-2))])
max_foo_lft=max([0.0,min_foo_lft])
cor_ad_flux_lft(num_layer)=S_lft(num_layer)*max_foo_lft
;
min_foo_rgt=min([abs(raw_ad_flux_rgt(num_layer)),$
	S_rgt(num_layer)*warp_dz(num_layer)*(trans_diff_conc(num_layer)-$
trans_diff_conc(num_layer-1))])
max_foo_rgt=max([0.0,min_foo_rgt])
cor_ad_flux_rgt(num_layer)=S_rgt(num_layer)*max_foo_rgt
;
min_foo_lft=min([abs(raw_ad_flux_lft(1)),$
	S_lft(1)*warp_dz(1)*(trans_diff_conc(1+1)-$
trans_diff_conc(1))])
max_foo_lft=max([0.0,min_foo_lft])
cor_ad_flux_lft(1)=S_lft(1)*max_foo_lft
;
for layer=1,num_layer do begin
;
; Eq 8-3.33 Finds rho super n, the new concentrations
;
	new_concentration(layer)=trans_diff_conc(layer)-$
	(cor_ad_flux_rgt(layer)-cor_ad_flux_lft(layer))/warp_dz(layer)
endfor
;new_concentration(0)=new_concentration(num_layer)
;new_concentration(num_layer+1)=new_concentration(1)
new_concentration(0)=conc_LBC_mult*new_concentration(1)+conc_LBC_add
new_concentration(num_layer+1)=conc_RBC_mult*new_concentration(num_layer)+conc_RBC_add
;
concentration=new_concentration
;
;Now reset some quantities that will be used before the next loop
for layer=1,num_layer do begin
;
; Eq 8-3.17
	conc_rgt(layer)=0.5*(concentration(layer+1)+concentration(layer))
	conc_lft(layer)=0.5*(concentration(layer)+concentration(layer-1))
endfor
;
IWP(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
column_density(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer));#/m^2
;
cdt=courant_safety/max(abs(zdzt(1:num_layer))/warp_dz(1:num_layer))
if dt gt cdt then begin
	print,'WARNING: dt = ',dt,' and cdt = ',cdt,$
	' in time step ',time_step
	dt=cdt
endif  ;endif dt
;
endfor ;end main time loop
;
print,'time to compute = ',sys_time
;
exit_gracefully: foo=1
;
end

pro plt_fct,time_array=time_array,$
	dt=dt,warp_dz=warp_dz,num_layer=num_layer, $
	time_step=time_step, $
	concentration=concentration, $
	column_density=column_density, $
	crystal_length=crystal_length,pause=pause, $
	crystal_mass=crystal_mass,$
	altitude=altitude,$
	junk=junk,zdzt=zdzt,$
	IWP=IWP
;
; Routines to plot the advection data at intermittent points
;
time_sng=string(Format='(F8.2)',time_array(time_step-1))
column_density_sng=string(Format='(F10.2)',column_density(time_step-1));#/m^2
dt_sng=string(Format='(F6.2)',dt)
num_layer_sng=string(Format='(I0)',num_layer)
IWP_sng=string(Format='(F6.2)',IWP(time_step-1)*1.0e6);mg/m^2
dz_sng=string(Format='(I4)',min(warp_dz))
zdzt_sng=string(Format='(F6.2)',max(abs(zdzt(1:num_layer))))
;
	plot,concentration(1:num_layer),$
		altitude(1:num_layer)/1000.0,$
		tit='!6Concentration vs Altitude', $
		ytit='Altitude (km)', ysty=1, $
		xtit='!6Concentration (#/m!E3!N)', xsty=0
;
xyouts,0.2*max(concentration),10.9,'!5Time = '+time_sng+' s',size=1.0
xyouts,0.2*max(concentration),10.8,'!5IWP = '+IWP_sng+' mg/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.7,'!5N!ITOT!N = '+column_density_sng+' #/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.6,'!5# layers = '+num_layer_sng,size=1.0
xyouts,0.2*max(concentration),10.5,'!5dz!IMIN!N = '+dz_sng+' m',size=1.0
xyouts,0.2*max(concentration),10.4,'!5dt = '+dt_sng+' s',size=1.0
xyouts,0.2*max(concentration),10.3,'!5dz/dt!IMAX!N = '+zdzt_sng+' m/s',size=1.0
;
;print,'altitude (km) = ',altitude/1000.
;print,'concentration (#/m^3) = ',concentration(*)
;
if pause eq 'y' then begin
	print,' Hit any key to continue or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
endif
;
exit_gracefully: foo=1
;
end

