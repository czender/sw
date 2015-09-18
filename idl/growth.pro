;
; $Id$
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-15 02:07:51  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:50  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:19  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:52  zender
; Imported sources
;
; Revision 1.1  1992/07/24  21:20:14  zender
; Initial revision
;
function Fofpsiu,psi_i,psi_ip1,u,dt,dz
	F=(u+abs(u))*psi_i+(u-abs(u))*psi_ip1
	F=0.5*F*dt/dz
return,F
end

pro growth,num_layer=num_layer,alpha=alpha, $
	nsteps=nsteps,dt=dt, $
	pltstep=pltstep,pause=pause,$
	sc=sc,bpd=bpd
;
; The extension of advect.pro to use the smolarkiewicz scheme.
;
; Example usage 
;
;growth,pause='n',nsteps=61,pltstep=60
;growth,pause='n',nsteps=601,pltstep=600,
;growth,pause='n',nsteps=601,pltstep=600,num_layer=50
;growth,pause='n',nsteps=201,pltstep=50,num_layer=50
;
; This procedure attempts to implement the advection method of Smolarkiewicz.
;
sys_time = systime(1) 
if n_elements(num_layer) eq 0 then num_layer = 50
if n_elements(nsteps) eq 0 then nsteps = 60
if n_elements(dt) eq 0 then dt = 1.			;second
if n_elements(pltstep) eq 0 then pltstep = 10		;time steps per graph
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if pause eq 'y' then !p.multi=0 else begin
	!p.multi=[0,2,0]
	pause = 'n'
endelse
if n_elements(alpha) eq 0 then alpha = .2
if n_elements(sc) eq 0 then sc = 1.;factor to increase twiddle velocity in the antidiffusion step
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
dz=cloud_thick/num_layer
epsilon=1.0e-15
;
time_array=fltarr(nsteps+1)
IWP=fltarr(nsteps+1)
smolar_IWP=fltarr(nsteps+1)
column_density=fltarr(nsteps+1)
smolar_col_den=fltarr(nsteps+1)
altitude=fltarr(num_layer+2)
altitude_rgt=fltarr(num_layer+2)
altitude_lft=fltarr(num_layer+2)
warp_dz=fltarr(num_layer+2)
warp_dz_rgt=fltarr(num_layer+2)
warp_dz_lft=fltarr(num_layer+2)
concentration=fltarr(num_layer+2)
new_concentration=fltarr(num_layer+2)
smolar_conc=fltarr(num_layer+2)
inter_smolar_conc=fltarr(num_layer+2)
new_smolar_conc=fltarr(num_layer+2)
zdzt=fltarr(num_layer+2)
zdztleft=fltarr(num_layer+2)
zdztright=fltarr(num_layer+2)
twiddleleft=fltarr(num_layer+2)
twiddleright=fltarr(num_layer+2)
;
; Arithmetical progression in altitude grid
;
;altitude(1:num_layer)=cloud_base+findgen(num_layer)*cloud_thick/(num_layer-1)+dz/2.
;altitude(0)=altitude(num_layer)
;altitude(num_layer+1)=altitude(1)
;
; Geometrical progression altitude grid
;
bins_per_doubling=float(bpd)
smallest_thick=cloud_thick/(2.^((num_layer-1.+.5)/bins_per_doubling))
for layer=1,num_layer do begin
	altitude_lft(layer)=cloud_base+smallest_thick*2.^((layer-1.-.5)/bins_per_doubling)
        altitude(layer)=cloud_base+smallest_thick*2.^((layer-1.)/bins_per_doubling)
	altitude_rgt(layer)=cloud_base+smallest_thick*2.^((layer-1.+.5)/bins_per_doubling)
endfor
altitude_lft(1)=cloud_base
for layer=1,num_layer do begin
	warp_dz(layer)=altitude_rgt(layer)-altitude_lft(layer)
endfor
altitude(0)=altitude(num_layer)
altitude(num_layer+1)=altitude(1)
altitude_rgt(0)=altitude_rgt(num_layer)
altitude_rgt(num_layer+1)=altitude_rgt(1)
altitude_lft(0)=altitude_lft(num_layer)
altitude_lft(num_layer+1)=altitude_lft(1)
warp_dz(0)=warp_dz(num_layer)
warp_dz(num_layer+1)=warp_dz(1)
;
for layer=1,num_layer do begin
	warp_dz_rgt(layer)=0.5*(warp_dz(layer+1)+warp_dz(layer))
	warp_dz_lft(layer)=0.5*(warp_dz(layer)+warp_dz(layer-1))
;	warp_dz_rgt(layer)=altitude_rgt(layer+1)-altitude_rgt(layer)
;	warp_dz_lft(layer)=altitude_lft(layer)-altitude_lft(layer-1)
endfor
warp_dz_rgt(num_layer)=warp_dz(num_layer)
warp_dz_lft(1)=warp_dz(1)
;
print,'smallest thick = ',smallest_thick
print,'cloud base = ',cloud_base
print,'cloud thick = ',cloud_thick
print,'altitude = ',altitude
print,'altitude_lft = ',altitude_lft
print,'altitude_rgt = ',altitude_rgt
print,'warp_dz = ',warp_dz
print,'total warp_dz = ',total(warp_dz(1:num_layer))
;
; Note that warp_dz corresponds to the capital Lambda in Oran & Boris '87
; ch 8-3 on FCT.
;
for layer=1,num_layer do begin
	if altitude(layer) lt  10000. then $
		concentration(layer)=0. else $
	if altitude(layer) ge 10100. then $
		concentration(layer)=0. else $
	concentration(layer)=bin_sz*1.0e6* $
		sz_dst_of_IWC(crystal_length,env_temperature)
	zdzt(layer)=-1.
endfor
concentration(0)=concentration(num_layer)
concentration(num_layer+1)=concentration(1)
;zdzt(1:num_layer)=-1-sin(2.*!pi*findgen(num_layer)/(num_layer-1))	;1 m/s
;zdzt(1:num_layer)=-1-findgen(num_layer)/(num_layer-1)	;1 m/s
zdzt(0)=zdzt(num_layer)
zdzt(num_layer+1)=zdzt(1)
;
smolar_conc=concentration
;
; Note how zdztleft(0),zdztright(num_layer+1) etc are undefined
; but are never used anyway
;
for layer=1,num_layer do begin
;ie velocity on the left interface of cell # layer
	zdztleft(layer)=0.5*(zdzt(layer-1)+zdzt(layer))
;ie velocity on the right interface of cell # layer
	zdztright(layer)=0.5*(zdzt(layer)+zdzt(layer+1))
endfor	
;
; Courant stability
;
cdt=alpha/max(abs(zdzt(1:num_layer))/warp_dz(1:num_layer))
print,'dt = ',dt,' and cdt = ',cdt
if dt gt cdt then begin
	dt=cdt
endif  ;endif dt
;
time_array(0)=0.
IWP(0)=total(concentration(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
column_density(0)=total(concentration(1:num_layer)*warp_dz(1:num_layer));#/m^2
smolar_IWP(0)=total(smolar_conc(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
smolar_col_den(0)=total(smolar_conc(1:num_layer)*warp_dz(1:num_layer));#/m^2
;
print,''
print,'Advection Model:'
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
plt_smolar,time_array=time_array,$
	dt=dt,warp_dz=warp_dz,num_layer=num_layer, $
	time_step=time_step, $
	concentration=concentration, $
	column_density=column_density,smolar_col_den=smolar_col_den, $
	smolar_conc=smolar_conc, $
	crystal_length=crystal_length,pause=pause, $
	crystal_mass=crystal_mass,$
	altitude=altitude,$
	junk=junk,zdzt=zdzt,$
	IWP=IWP,sc=sc,$
	smolar_IWP=smolar_IWP
if n_elements(junk) ne 0 then if junk eq 'q' then goto, exit_gracefully
;
endif  ;end if time to plot
;
for layer=1,num_layer do begin
;
	new_concentration(layer)=concentration(layer)-$
Fofpsiu(concentration(layer),concentration(layer+1),$
	zdztright(layer),dt,warp_dz(layer))+$
Fofpsiu(concentration(layer-1),concentration(layer),$
	zdztleft(layer),dt,warp_dz(layer))
;
	inter_smolar_conc(layer)=smolar_conc(layer)-$
Fofpsiu(smolar_conc(layer),smolar_conc(layer+1),$
	zdztright(layer),dt,warp_dz(layer))+$
Fofpsiu(smolar_conc(layer-1),smolar_conc(layer),$
	zdztleft(layer),dt,warp_dz(layer))
;
endfor	;end loop over layers
;
new_concentration(0)=new_concentration(num_layer)
new_concentration(num_layer+1)=new_concentration(1)
;
inter_smolar_conc(0)=inter_smolar_conc(num_layer)
inter_smolar_conc(num_layer+1)=inter_smolar_conc(1)
;
for layer=1,num_layer do begin
;ie twiddle velocity on the left interface of cell # layer
	twiddleleft(layer)=(abs(zdztleft(layer))*warp_dz(layer)-$
	dt*zdztleft(layer)*zdztleft(layer))*$
	(inter_smolar_conc(layer)-inter_smolar_conc(layer-1))/$
	(warp_dz(layer)*(inter_smolar_conc(layer-1)+inter_smolar_conc(layer)+epsilon))
;ie twiddle velocity on the right interface of cell # layer
	twiddleright(layer)=(abs(zdztright(layer))*warp_dz(layer)-$
	dt*zdztright(layer)*zdztright(layer))*$
	(inter_smolar_conc(layer+1)-inter_smolar_conc(layer))/$
	(warp_dz(layer)*(inter_smolar_conc(layer)+inter_smolar_conc(layer+1)+epsilon))
endfor	;end loop over layers
;
twiddleleft=sc*twiddleleft
twiddleright=sc*twiddleright
;
for layer=1,num_layer do begin
;
	new_smolar_conc(layer)=inter_smolar_conc(layer)-$
Fofpsiu(inter_smolar_conc(layer),inter_smolar_conc(layer+1),$
	twiddleright(layer),dt,warp_dz(layer))+$
Fofpsiu(inter_smolar_conc(layer-1),inter_smolar_conc(layer),$
	twiddleleft(layer),dt,warp_dz(layer))
;
endfor	;end loop over layers
;
new_smolar_conc(0)=new_smolar_conc(num_layer)
new_smolar_conc(num_layer+1)=new_smolar_conc(1)
;
flux_out_btm=-Fofpsiu(concentration(0),concentration(1),$
	zdztleft(1),dt,warp_dz(1))
flux_out_top=Fofpsiu(concentration(num_layer),concentration(num_layer+1),$
	zdztright(num_layer),dt,warp_dz(num_layer))
;
; Now swap the arrays for the beginning of the next time step
;
concentration=new_concentration
smolar_conc=new_smolar_conc
;
IWP(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
smolar_IWP(time_step)=total(smolar_conc(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
column_density(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer));#/m^2
smolar_col_den(time_step)=total(smolar_conc(1:num_layer)*warp_dz(1:num_layer));#/m^2
;
cdt=alpha/max(abs(zdzt(1:num_layer))/warp_dz(1:num_layer))
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

pro plt_smolar,time_array=time_array,$
	dt=dt,warp_dz=warp_dz,num_layer=num_layer, $
	time_step=time_step, $
	concentration=concentration, $
	column_density=column_density,smolar_col_den=smolar_col_den, $
	smolar_conc=smolar_conc, $
	crystal_length=crystal_length,pause=pause, $
	crystal_mass=crystal_mass,$
	altitude=altitude,$
	junk=junk,zdzt=zdzt,$
	IWP=IWP,sc=sc,$
	smolar_IWP=smolar_IWP
;
; Routines to plot the advection data at intermittent points
;
time_sng=string(Format='(F8.2)',time_array(time_step-1))
column_density_sng=string(Format='(F10.2)',column_density(time_step-1))
smolar_col_den_sng=string(Format='(F10.2)',smolar_col_den(time_step-1))
max_conc_sng=string(Format='(F8.2)',max(concentration(1:num_layer)))
max_smolar_conc_sng=string(Format='(F8.2)',max(smolar_conc(1:num_layer)))
dt_sng=string(Format='(F6.2)',dt)
sc_sng=string(Format='(F5.2)',sc)
num_layer_sng=string(Format='(I0)',num_layer)
IWP_sng=string(Format='(F6.2)',IWP(time_step-1)*1.0e6)	;mg/m^2
smolar_IWP_sng=string(Format='(F6.2)',smolar_IWP(time_step-1)*1.0e6);mg/m^2
dz_sng=string(Format='(I4)',min(warp_dz))
zdzt_sng=string(Format='(F6.2)',max(abs(zdzt(1:num_layer))))
;
	plot,smolar_conc(1:num_layer),$
		altitude(1:num_layer)/1000.0,$
		tit='!6Concentration vs Altitude', $
;		yrange=[10,10.1],$
		ytit='Altitude (km)', ysty=1, $
		xtit='!6Concentration (#/m!E3!N)', xsty=0
	oplot,concentration(1:num_layer),$
		altitude(1:num_layer)/1000.0,$
		linestyle=1
;
xyouts,0.2*max(concentration),10.9,'!5Time = '+time_sng+' s',size=1.0
xyouts,0.2*max(concentration),10.85,'!5Sc = '+sc_sng,size=1.0
xyouts,0.2*max(concentration),10.8,'!5IWP!IS!N = '+smolar_IWP_sng+' mg/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.75,'!5N!IS,TOT!N = '+smolar_col_den_sng+' #/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.7,'!5n!IS,MAX!N = '+max_smolar_conc_sng+' #/m!E3!N',size=1.0

xyouts,0.2*max(concentration),10.6,'!5IWP = '+IWP_sng+' mg/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.55,'!5N!ITOT!N = '+column_density_sng+' #/m!E2!N',size=1.0
xyouts,0.2*max(concentration),10.5,'!5n!IMAX!N = '+max_conc_sng+' #/m!E3!N',size=1.0
xyouts,0.2*max(concentration),10.3,'!5# layers = '+num_layer_sng,size=1.0
xyouts,0.2*max(concentration),10.25,'!5dz = '+dz_sng+' m',size=1.0
xyouts,0.2*max(concentration),10.2,'!5dt = '+dt_sng+' s',size=1.0
xyouts,0.2*max(concentration),10.15,'!5dz/dt!IMAX!N = '+zdzt_sng+' m/s',size=1.0
;
;print,'altitude (km) = ',altitude/1000.
;print,'concentration (#/m^3) = ',concentration(*)
;print,''
;print,'smolar concentration (1:nl)(#/m^3) = ',smolar_conc(*)
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

