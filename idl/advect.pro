;
; $Id$
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-15 02:07:48  zender
; *** empty log message ***
;
; Revision 1.5  2000/01/10 23:36:26  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:46  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:13  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:45  zender
; Imported sources
;
; Revision 1.1  1992/07/27  15:50:33  zender
; Initial revision
;
function Fofpsiu,psi_i,psi_ip1,u,dt,dx
	F=(u+abs(u))*psi_i+(u-abs(u))*psi_ip1
	F=0.5*F*dt/dx
return,F
end

pro advect,num_layer=num_layer,alpha=alpha, $
	nsteps=nsteps,dt=dt, $
	pltstep=pltstep,pause=pause,$
	bpd=bpd
;
; Does straight upstream advection on the vertical advection equation.
; Is built into smolar.pro, where anti-diffusive steps are added.
;
; Example usage 
;
;advect,pause='n',nsteps=61,pltstep=60
;advect,pause='n',nsteps=601,pltstep=600,dt=2
;advect,pause='n',nsteps=601,pltstep=600,dt=10.0,num_layer=50
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
;
time_array=fltarr(nsteps+1)
IWP=fltarr(nsteps+1)
column_density=fltarr(nsteps+1)
concentration=fltarr(num_layer+2)
new_concentration=fltarr(num_layer+2)
zdzt=fltarr(num_layer+2)
zdztleft=fltarr(num_layer+2)
zdztright=fltarr(num_layer+2)
altitude=fltarr(num_layer+2)
altitude_rgt=fltarr(num_layer+2)
altitude_lft=fltarr(num_layer+2)
warp_dz=fltarr(num_layer+2)
warp_dz_rgt=fltarr(num_layer+2)
warp_dz_lft=fltarr(num_layer+2)
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
smallest_thick=cloud_thick/(2.^((num_layer-1.)/bins_per_doubling))
for layer=1,num_layer do begin
	altitude_lft(layer)=cloud_base+smallest_thick*2.^((layer-1.-.5)/bins_per_doubling)
        altitude(layer)=cloud_base+smallest_thick*2.^((layer-1.)/bins_per_doubling)
	altitude_rgt(layer)=cloud_base+smallest_thick*2.^((layer-1.+.5)/bins_per_doubling)
	warp_dz(layer)=altitude_rgt(layer)-altitude_lft(layer)
endfor
altitude(0)=altitude(num_layer)
altitude(num_layer+1)=altitude(1)
;altitude_rgt(0)=altitude_rgt(num_layer)
;altitude_rgt(num_layer+1)=altitude_rgt(1)
;altitude_lft(0)=altitude_lft(num_layer)
;altitude_lft(num_layer+1)=altitude_lft(1)
warp_dz(0)=warp_dz(num_layer)
warp_dz(num_layer+1)=warp_dz(1)
;
; Note that warp_dz corresponds to the capital Lambda in Oran & Boris '87
; ch 8-3 on FCT.
;
; Note how warp_dz_lft(0),warp_dz_rgt(num_layer+1) etc are undefined
; but are never used anyway
;
for layer=1,num_layer do begin
	warp_dz_rgt(layer)=0.5*(warp_dz(layer+1)+warp_dz(layer))
	warp_dz_lft(layer)=0.5*(warp_dz(layer)+warp_dz(layer-1))
;	warp_dz_rgt(layer)=altitude_rgt(layer+1)-altitude_rgt(layer)
;	warp_dz_lft(layer)=altitude_lft(layer)-altitude_lft(layer-1)
endfor
warp_dz_rgt(num_layer)=warp_dz(num_layer)
warp_dz_lft(1)=warp_dz(1)
print,'smallest thick = ',smallest_thick
print,'altitude = ',altitude
print,'warp_dz = ',warp_dz
print,'warp_dz_lft = ',warp_dz_lft
print,'total warp_dz(1:nl) = ',total(warp_dz(1:num_layer))
print,'total warp_dz_lft(1:nl) = ',total(warp_dz_lft(1:num_layer))
;
for layer=1,num_layer do begin
endfor

for layer=1,num_layer do begin
	if altitude(layer) lt  10800. then $
		concentration(layer)=0. else $
	if altitude(layer) ge 10900. then $
		concentration(layer)=0. else $
	concentration(layer)=bin_sz*1.0e6* $
		sz_dst_of_IWC(crystal_length,env_temperature)
;	zdzt(layer)=-1.
endfor
concentration(0)=concentration(num_layer)
concentration(num_layer+1)=concentration(1)
;zdzt(1:num_layer)=-1-sin(2.*!pi*findgen(num_layer)/(num_layer-1))	;1 m/s
zdzt(1:num_layer)=-1-findgen(num_layer)/(num_layer-1)	;1 m/s
zdzt(0)=zdzt(num_layer)
zdzt(num_layer+1)=zdzt(1)
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
plt_advect,time_array=time_array,$
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
	new_concentration(layer)=concentration(layer)-$
Fofpsiu(concentration(layer),concentration(layer+1),$
	zdztright(layer),dt,warp_dz(layer))+$
Fofpsiu(concentration(layer-1),concentration(layer),$
	zdztleft(layer),dt,warp_dz(layer))
;
endfor	;end loop over layers
;
new_concentration(0)=new_concentration(num_layer)
new_concentration(num_layer+1)=new_concentration(1)
;
flux_out_btm=-Fofpsiu(concentration(0),concentration(1),$
	zdztleft(1),dt,warp_dz(1))
flux_out_top=Fofpsiu(concentration(num_layer),concentration(num_layer+1),$
	zdztright(num_layer),dt,warp_dz(num_layer))
;
concentration=new_concentration
;
IWP(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer))*crystal_mass
column_density(time_step)=total(concentration(1:num_layer)*warp_dz(1:num_layer));#/m^2
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

pro plt_advect,time_array=time_array,$
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

