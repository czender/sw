;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: tester.pro,v $
; $Source: /home/zender/cvs/idl/tester.pro,v $
; $Id$
; $State: Exp $
;
; Purpose: 
; routines in this file are meant to be a testing ground for different
; methods of computing the thermal and absorptive diffusive parameters,
; k and d, respectively, as well as Quality Assurance (QA) routines for
; other parameterizations. Many of the routines to be tested are in
; clouds.pro.
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-15 02:07:54  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:54  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:29  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:51  zender
; Imported sources
;
; Revision 1.6  1994/01/13  03:42:17  zender
; added haze_tester, CN_tester to look at DeMott haze freezing
; rates and Jensen CN distribution.
;
; Revision 1.5  1993/06/07  23:53:15  zender
; synchronozation version. removed the old lagrangian idl
; clouds() program and plotting routines.
;
; Revision 1.1  1993/06/07  23:45:58  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO Pdf_tester, INC_init=INC_init, pause=pause, $
                abc_nbr=abc_nbr 
if n_elements(INC_init) eq 0 then INC_init = 150.e3
if n_elements(abc_nbr) eq 0 then abc_nbr = 10
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if pause eq 'y' then !p.multi=0 else begin
	!p.multi=[0,3,1,0,0];=[0,num_cols,num_rows,0,0]
	pause = 'n'
endelse
max_length=1000.e-6	;meters
crystal_length=(findgen(abc_nbr)+1.)*max_length/abc_nbr
delta_length=shift(crystal_length,-1)-crystal_length
delta_length(abc_nbr-1)=delta_length(abc_nbr-2)
PDF=PDF_of_length(crystal_length)
distribution=INC_init*PDF
concentration=distribution*delta_length

;print,"crystal_length (microns) = "
;print,crystal_length*1.0e6
;print,"delta_length (microns) = "
;print,delta_length*1.0e6
;print,"PDF (per micron) = "
;print,PDF*1.0e-6
;print,"distribution (#/liter/micron) = "
;print,distribution*1.0e-3*1.0e-6
;print,"concentration (#/liter) = "
;print,concentration*1.0e-3

print,"total concentration (#/liter) = ",total(concentration)*1.0e-3
print,"total concentration (#/liter) >= 10 microns = ", $
	total(concentration(9:abc_nbr-1))*1.0e-3
print,"total concentration (#/liter) 3 <= L <= 10 microns = ", $
	total(concentration(2:9))*1.0e-3

plot_oo,crystal_length(*)*1.0e6,PDF*1.0e-6, $
	tit='!5PDF of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	ytit='!5Probability Density Function (!7l!5!E-1!N)',ysty=0

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oo,crystal_length(*)*1.0e6,PDF*1.0e-6, $
	tit='!5PDF of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	ytit='!5Probability Density Function (!7l!5!E-1!N)',ysty=0
;******************** End Printing routine ******************************
	lpr
   endif	
 endif ;endif pause

plot_oo,crystal_length(*)*1.0e6,distribution*1.0e-3*1.0e-6, $
	tit='!5Distribution of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	yrange=[1.0e-3,max(distribution(*))*1.0e-3*1.0e-6], $
	ytit='!5Distribution (#/liter/!7l!5)',ysty=0

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oo,crystal_length(*)*1.0e6,distribution*1.0e-3*1.0e-6, $
	tit='!5Distribution of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	yrange=[1.0e-3,max(distribution(*))*1.0e-3*1.0e-6], $
	ytit='!5Distribution (#/liter/!7l!5)',ysty=0
;******************** End Printing routine ******************************
	lpr
   endif	
 endif ;endif pause

plot_oo,crystal_length(*)*1.0e6,concentration(*)*1.0e-3, $
	tit='!5Concentration of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	yrange=[0.1,max(concentration(*))*1.0e-3], $
	psym=1, $
	ytit='!5Crystal Concentration (#/liter)',ysty=0

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oo,crystal_length(*)*1.0e6,concentration(*)*1.0e-3, $
	tit='!5Concentration of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	yrange=[0.1,max(concentration(*))*1.0e-3], $
	psym=1, $
	ytit='!5Crystal Concentration (#/liter)',ysty=0
;******************** End Printing routine ******************************
	lpr
   endif	
 endif ;endif pause

;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

pro condense,abc_nbr=abc_nbr,pressure=pressure, $
	vap_mix_ratio=vap_mix_ratio
; this procedure attempts to implement the cirrus ice particle
; growth model Rogers and Yau
;
if n_elements(abc_nbr) eq 0 then abc_nbr = 100
if n_elements(pressure) eq 0 then pressure = 25000	;Pa
if n_elements(vap_mix_ratio) eq 0 then vap_mix_ratio = .001;pure number
;if n_elements() eq 0 then = 
;
; Measured constants
;
l_s = 2.834e6			;J/kg  (latent heat of sublimation)
R_v = 461.5			;J/kg/K (gas constant for vapor)
epsilon = .622			;molec wgt vapor/molec wgt dry air
temperature = 283.15-findgen(abc_nbr)*50./(abc_nbr-1)
;
; now put together eqn 9.4 p. 160 Rogers and Yau.  
; note how conduction term is pressure independent
;
conduction_term=l_s*(-1.+l_s/(R_v*temperature))/ $
	(temperature*k_rogers(temperature))
;
; loop over suitable range of pressures
;
for pressure_mb=100,1000,100 do begin
;
diffusion_term=R_v*temperature/(eqm_vap_ice(temperature)* $
	D_rogers(temperature,pressure_mb*100.))
;
; Now assume a water saturated growth environment, i.e. initial growth
; of ice crystals is occuring in a water cloud.
;
S_i = eqm_vap_water(temperature)/eqm_vap_ice(temperature) 
;
; Now instead of assumption, use the vapor mixing ratio, and remember
; to correct pressure from mbar into SI
;
S_i_2 = vap_mix_ratio*pressure_mb*100./(epsilon*eqm_vap_ice(temperature))
;
dmdt=(S_i-1.)/(conduction_term+diffusion_term)
dmdt_2=(S_i_2-1.)/(conduction_term+diffusion_term)
;
; right now dmdt is in SI unit: kg/m/s, convert to cgs: gm/cm/s by
; multiplying by 1000 and dividing by 100
;
dmdt_cgs=dmdt*10
dmdt_2_cgs=dmdt_2*10
if pressure_mb eq 100 then $
	plot,temperature-273.15,1.0e8*dmdt_cgs,linestyle=0, $
	xtitle='!5Temperature, Celsius', $
	ytitle='!5 dM/dt 10!E8!N/4!7p!5C g cm!E-1!Ns!E-1!N', $
	title='Ice Crystal Growth Rate',xrange=[-40,0] $
else $
	oplot,temperature-273.15,1.0e8*dmdt_cgs,linestyle=0
oplot,temperature-273.15,1.0e8*dmdt_2_cgs,linestyle=1
endfor		;end loop over pressures
end

function K_rogers,temperature
; looks up the temperature dependent coefficient of thermal conductivity
; of air, K, from value given in Rogers & Yau p. 103.
; unit of K are J/m/s/K
; REMEMBER to input temperature in Kelvin!
;
t=temperature-273.15		;convert to celsius
K=t				;dimensions the array if necc.
n=n_elements(t)
if n gt 1 then begin
	for i = 0,n-1 do begin
		if t(i) gt 25 then K(i) = 2.63e-2 else $
		if t(i) gt 15 then K(i) = 2.55e-2 else $
		if t(i) gt 5 then K(i) = 2.48e-2 else $
		if t(i) gt -5 then K(i) = 2.40e-2 else $
		if t(i) gt -15 then K(i) = 2.32e-2 else $
		if t(i) gt -25 then K(i) = 2.24e-2 else $
		if t(i) gt -35 then K(i) = 2.16e-2 else $
		K(i) = 2.07e-2
	endfor
	return,K
endif else begin
	if t gt 25 then return,2.63e-2 else $
	if t gt 15 then return,2.55e-2 else $
	if t gt 5 then return,2.48e-2 else $
	if t gt -5 then return,2.40e-2 else $
	if t gt -15 then return,2.32e-2 else $
	if t gt -25 then return,2.24e-2 else  $
	if t gt -35 then return,2.16e-2 else $
	return,2.07e-2
endelse
end

function K_sexena,vap_pres=vap_pres,pressure=pressure, $
	temperature=temperature
;finds thermal of moist air by combining the thermal of dry air with
;the thermal conductivity of wator vapor using the mason-sexena 
;formula on pruppaker and Klett p.418
;formula has k in cal/cm/s/C but 
;k is returned in SI J/m/s/C
;
if n_elements(vap_pres) eq 0 then vap_pres = 25
if n_elements(pressure) eq 0 then pressure = 1000
;if n_elements() eq 0 then = 
gamma1=1.17
gamma2=1.02
t_celsius=temperature-273.15
xsubv=vap_pres/pressure
;
ksuba=(5.69+0.017*t_celsius)*1.0e-5
ksubv=(3.78+0.02*t_celsius)*1.0e-5
;
k_sexena=ksuba*(1.-xsubv*(gamma1-gamma2*ksubv)/ksuba)
k_sexena=k_sexena*4.1855*100
;
return,k_sexena
end

function K_zhang,temperature=temperature,radius=radius, $
	deltat=deltat,alphat=alphat,env_density=env_density, $
	ksuba=ksuba,pressure=pressure,fsubpr=fsubpr
;
; computes the correction to ksuba for small particles due
; to the thermal/kinetic effects close to the particle.
; Zhang et al. p. 309 reference is Laube and Holler '88
; unit of K_zhang are J/m/s/K
; REMEMBER to input temperature in Kelvin!
;
if n_elements(env_density) eq 0 then env_density = .3	;kg/m^3
if n_elements(alphat) eq 0 then alphat = .7
if n_elements(radius) eq 0 then radius = 5e-6		;m
if n_elements(temperature) eq 0 then temperature = 228  ;K
if n_elements(ksuba) eq 0 then ksuba = .04		;J/m/s/K
if n_elements(pressure) eq 0 then pressure = 250	;mb
if n_elements(fsubpr) eq 0 then fsubpr = 1.		
;if n_elements() eq 0 then = 
;
spec_heat_air=1005		;J/kg/K
gas_const_dry_air=287.05	;J/kg/K
gas_const_vapor=461.51	;J/kg/K
B=4.62e-8*(1013.3/pressure)*(temperature/293.15)
;
speed=2.*!pi/(gas_const_vapor*temperature)
speed=sqrt(speed)
denom_1=radius/(radius+B)
denom_2=ksuba*fsubpr/(radius*alphat*env_density*spec_heat_air)
denom_2=denom_2*speed
K_zhang=ksuba/(denom_1+denom_2)
return,K_zhang
end


pro k_tester
temperature=findgen(100)/2.+243.15
t_celsius=temperature-273.15
;
k_sexena=k_sexena(temperature=temperature)
ksubv=(3.78+0.02*t_celsius)*1.0e-5
ksubv=ksubv*4.1855*100
ksuba=(5.69+0.017*t_celsius)*1.0e-5
ksuba=ksuba*4.1855*100
k_rogers=k_rogers(temperature)
k_pruppacher=k_pruppacher(temperature=temperature)
k_zhang=k_zhang(temperature=temperature)

plot,temperature,k_zhang,xtitle='Temperature (K)',ytitle='K (J/m/s/C)', $
	title='Thermal Conductivity k from four authors'
oplot,temperature,ksuba,linestyle=1
oplot,temperature,k_rogers,linestyle=0
oplot,temperature,k_pruppacher,linestyle=1
oplot,temperature,k_sexena,linestyle=2
end

pro D_tester
temperature=findgen(100)/2.+223.15
t_celsius=temperature-273.15
;
D_pruppacher=D_pruppacher(temperature=temperature,env_pressure=env_pressure)
D_rogers=D_rogers(temperature,env_pressure)
;
plot,temperature,D_pruppacher,xtitle='Temperature (K)',ytitle='K (J/m/s/C)', $
	title='Thermal Conductivity k from four authors'
oplot,temperature,D_rogers,linestyle=1
end

pro sz_tester
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr
ngtr=fltarr(abc_nbr)
dr_mcr=float(max_length)/abc_nbr
env_temperature=228.
dist=sz_dst_of_IWC(crystal_length(*),env_temperature)

for i=0,abc_nbr-1 do begin
; find the number greater than a given length
	ngtr(i)=total(dist(i:abc_nbr-1)*dr_mcr)
	if i lt abc_nbr-1 then begin
		ngtr(i)=ngtr(i)-.5*dist(i)*dr_mcr
	        ngtr(i)=ngtr(i)-.5*dist(abc_nbr-1)*dr_mcr
	endif
endfor
	plot_io,crystal_length(*)*1.0e6,ngtr(*), $
	        tit='!5Concentration, Distribution vs Length', $
		xtit='Length R (!7l!5)', xsty=1, $
		ytit='!5Concentration N(r > R) (#/m!E3!N) (solid)',ysty=8, $
		xmargin=[8,6]
	axis,yaxis=1,/save, $
		yrange=[min(dist),max(dist)], $
		ytit='!5Distribution (#/m!E3!N/!7l!5) (dashed)', $
		ystyle=1
	oplot,crystal_length(*)*1.0e6,dist,linestyle=2
;print,"dist = ",dist," #/m^3/micron at length = ",crystal_length*1.0e6," microns"
end

pro multi_sz_tester
tpt_nbr=8
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr
ngtr=fltarr(tpt_nbr,abc_nbr)
dist=fltarr(tpt_nbr,abc_nbr)
dr_mcr=float(max_length)/abc_nbr

env_temperature=[215.,220.0,225.,230.0,235.,240.0,245.,250]
for tpt_idx=0,tpt_nbr-1 do begin 
dist(tpt_idx,*)=sz_dst_of_IWC(crystal_length(*), $
	env_temperature(tpt_idx))

for i=0,abc_nbr-1 do begin
; find the number greater than a given length
	ngtr(tpt_idx,i)=total(dist(tpt_idx,i:abc_nbr-1)*dr_mcr)
	if i lt abc_nbr-1 then begin
		ngtr(tpt_idx,i)=ngtr(tpt_idx,i)- $
			.5*dist(tpt_idx,i)*dr_mcr
	        ngtr(tpt_idx,i)=ngtr(tpt_idx,i)- $
			.5*dist(abc_nbr-1)*dr_mcr
	endif
endfor	;abc
endfor  ;temperatures

	plot_io,crystal_length(*)*1.0e6,ngtr(tpt_nbr-1,*), $
	        tit='!5Concentration, Distribution vs Length, Temp.', $
		xtit='Length R (!7l!5)', xsty=0, $
		ytit='!5Concentration N(r > R) (#/m!E3!N) (solid)',ysty=8, $
		xrange=[10.0,max(crystal_length(*)*1.0e6)], $
		yrange=[1.0,max(ngtr(tpt_nbr-1))], $
		xmargin=[8,6],linestyle=0
;foo=get_kbrd(1)
for tpt_idx=tpt_nbr-2,0,-1 do begin 
	oplot,crystal_length(*)*1.0e6,ngtr(tpt_idx,*), $
	linestyle=0
;	foo=get_kbrd(1)
endfor

	axis,yaxis=1,/save, $
		xrange=[10,max(crystal_length(*)*1.0e6)], $
		yrange=[0.1,max(dist)], $
		ytit='!5Distribution (#/m!E3!N/!7l!5) (dashed)', $
		ystyle=1
for tpt_idx=tpt_nbr-1,0,-1 do begin 
	oplot,crystal_length(*)*1.0e6,dist(tpt_idx,*), $
	linestyle=2
;	foo=get_kbrd(1)
endfor
end

PRO density_tester, pause=pause, abc_nbr=abc_nbr
IF n_elements(abc_nbr) EQ 0 THEN abc_nbr = 100
IF n_elements(pause) EQ 0 THEN pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

max_temp = 273.15
min_temp = 223.15
temperature=(findgen(abc_nbr)+1.)*(max_temp-min_temp)/ $
		abc_nbr+min_temp
temperature_celsius=temperature-273.15

crystal_length=(findgen(abc_nbr)+1.)*1000.e-6/100.
crystal_dens2=.78*crystal_length^(-0.0038)	;g/m^3
crystal_dens2=crystal_dens2*1.0e6/1000.

crystal_dens=crystal_dens_of_length(crystal_length)
	plot,crystal_length(*)*1.0e6,crystal_dens(*)*1000/1.0e6, $
		tit='!5Density of Cirrostratus Ice Crystals', $
		xtit='Length (!7l!5)', xsty=0, $
		ytit='!5Ice Density (g/cm!E3!N)',ysty=0
	oplot,crystal_length(*)*1.0e6,crystal_dens2(*)*1000/1.0e6, $
		linestyle=2

IF pause EQ 'y' THEN BEGIN
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	lpr
   endif	
 endif ;endif pause

ice_dens=ice_dens_of_T(temperature)
liquid_dens=liquid_dens_of_T(temperature)

plot,temperature_celsius(*),liquid_dens(*)*1.0e-3 > .8, $
	tit='!5Density of pure Liquid and Ice Water', $
	xtit='!5Environmental Temperature (!12_!5C)', xsty=0, $
	ytit='!5Density !7q!5 (g-cm!E-3!N)', $
	ysty=0, $ 
	yrange=[.9,1.0], $
	charsize=1.5, thick = 2. 
oplot,temperature_celsius(*),ice_dens(*)*1.0e-3,linestyle=2,thick=2.

;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

pro mass_tester
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr

mass=mass_of_length(crystal_length)
	plot_oo,crystal_length(*)*1.0e6,mass(*)*1.0e9, $
		tit='!5Mass of Cirrostratus Ice Crystals', $
		xtit='Length (!7l!5)', xsty=0, $
		ytit='!5Crystal Mass (!7l!5g)',ysty=0, $
		xrange=[10.0,1000.0]

crystal_diameter=crystal_length/3.5
mass=liou_mass_of_length_diameter(crystal_length,crystal_diameter)
	oplot,crystal_length(*)*1.0e6,mass(*)*1.0e9, $
		linestyle=1

crystal_length=[20,50,120,300,750]
crystal_length=crystal_length*1.0e-6
crystal_diameter=[20,40,60,100,160]
crystal_diameter=crystal_diameter*1.0e-6
mass=liou_mass_of_length_diameter(crystal_length,crystal_diameter)
	oplot,crystal_length(*)*1.0e6,mass(*)*1.0e9, $
		linestyle=2

end

function alt_basis,x,m
n=n_elements(x)
r=fltarr(n,m)
for i=0,n-1 do begin
	r(i,0)=1.
	r(i,1)=exp(-x(i)/(0.00050))
endfor
return,r
end

function alt_basis_ords,x,coeffs
n=n_elements(x)
m=n_elements(coeffs)
y=fltarr(n)
for i=0,n-1 do begin
	y(i)=coeffs(0)*1.+coeffs(1)*exp(-x(i)/(0.00050))
endfor
;y=poly(abc,coeffs)
return,y
end

pro aspect_tester,m=m
;
; Example usage 
;
;hex_diameter_param,m=5
;
if n_elements(m) eq 0 then m = 3        ;m-1 is the degree of the polynomial
                                        ;m is the number of coeffs in fitting func.
crystal_length=[20,50,120,300,750]
;crystal_length=[10,20,50,120,300,750,1000,1500,2000]
crystal_length=crystal_length*1.0e-6
crystal_diameter=[20,40,60,100,160]
;crystal_diameter=[10,20,40,60,100,160,180,200,210]
crystal_diameter=crystal_diameter*1.0e-6
aspect_ratio=crystal_length/crystal_diameter;
weight=[3.0,3.0,3.0,2.0,1.0]
;weight=[3.0,3.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0]
num_points=n_elements(crystal_diameter)
yf=fltarr(num_points)

;coeffs=svdfit(crystal_length,aspect_ratio,m,weight=weight,YFIT=yf)
coeffs=svdfit(crystal_length,aspect_ratio,m,weight=weight,funct='alt_basis',YFIT=yf)
abc=findgen(100)*max(crystal_length)/99.
;ord=poly(abc,coeffs)
ord=alt_basis_ords(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,'aspect_ratio(L) coeffs = ',coeffs
plot,crystal_length*1.0e6,aspect_ratio, $
	tit='!5 crystal length vs. aspect ratio parameterization', $
        xtit='Crystal Length', $
        ytit='aspect ratio', $
        linestyle=1, thick=2.0,psym=6,/ynoz
oplot,abc*1.0e6,ord, $
        linestyle=0,thick=2.
oplot,crystal_length*1.0e6,yf,linestyle=0,thick=2.0,psym=2
end

pro hex_diameter_param,m=m
;
; Example usage 
;
;hex_diameter_param,m=5
;
if n_elements(m) eq 0 then m = 3        ;m-1 is the degree of the polynomial
                                        ;m is the number of coeffs in fitting func.
;crystal_length=[20,50,120,300,750]
crystal_length=[10,20,50,120,300,750,1000,1500,2000]
crystal_length=crystal_length*1.0e-6
;crystal_diameter=[20,40,60,100,160]
crystal_diameter=[10,20,40,60,100,160,180,200,210]
crystal_diameter=crystal_diameter*1.0e-6
aspect_ratio=crystal_length/crystal_diameter;
;weight=[3.0,3.0,2.0,2.0,1.0]
weight=[3.0,3.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0]
num_points=n_elements(crystal_diameter)
yf=fltarr(num_points)

coeffs=svdfit(crystal_length,crystal_diameter,m,weight=weight,YFIT=yf)
;coeffs=svdfit(crystal_length,crystal_diameter,m,weight=weight,funct='alt_basis',YFIT=yf)
abc=findgen(100)*max(crystal_length)/99.
ord=poly(abc,coeffs)
;ord=alt_basis_ords(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,'diameter(L) coeffs = ',coeffs
;print,'abc = ',abc
;print,'ord = ',ord

plot,crystal_length*1.0e6,crystal_diameter*1.0e6, $
	tit='!5 crystal length vs. diameter parameterization', $
        xtit='Crystal Length', $
        ytit='Crystal Diameter', $
        linestyle=1, thick=2.0,psym=6,/ynoz
oplot,abc*1.0e6,ord*1.0e6, $
        linestyle=0,thick=2.
oplot,crystal_length*1.0e6,yf*1.0e6,linestyle=0,thick=2.0,psym=2

print,' Hit any key to continue...'
junk = get_kbrd(1)

coeffs=svdfit(crystal_length,aspect_ratio,m,weight=weight,YFIT=yf)
abc=findgen(100)*max(crystal_length)/99.
ord=poly(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,'aspect_ratio(L) coeffs = ',coeffs
plot,crystal_length*1.0e6,aspect_ratio, $
	tit='!5 crystal length vs. aspect ratio parameterization', $
        xtit='Crystal Length', $
        ytit='aspect ratio', $
        linestyle=1, thick=2.0,psym=6,/ynoz
oplot,abc*1.0e6,ord, $
        linestyle=0,thick=2.
oplot,crystal_length*1.0e6,yf,linestyle=0,thick=2.0,psym=2

print,' Hit any key to continue...'
junk = get_kbrd(1)

crystal_mass=liou_mass_of_length_diameter(crystal_length,crystal_diameter)
log_crystal_mass=alog(crystal_mass)
print,'log_crystal_mass = ',log_crystal_mass
print,'crystal_length = ',crystal_length*1.0e6
coeffs=svdfit(log_crystal_mass,crystal_length,m,weight=weight,YFIT=yf)
abc=findgen(100)*max(log_crystal_mass)/99.
ord=poly(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,'length(log(M)) coeffs = ',coeffs
plot,log_crystal_mass,crystal_length*1.0e6, $
	tit='!5 log crystal mass vs. crystal length parameterization', $
        xtit='log Crystal Mass', $
        ytit='Crystal Length microns', $
        linestyle=1, thick=2.0,psym=6,/ynoz
oplot,abc,ord*1.0e6, $
        linestyle=0,thick=2.
oplot,log_crystal_mass,yf*1.0e6,linestyle=0,thick=2.0,psym=2

print,' Hit any key to continue...'
junk = get_kbrd(1)

crystal_mass=liou_mass_of_length_diameter(crystal_length,crystal_diameter)
print,'crystal_mass = ',crystal_mass*1.0e9
print,'crystal_length = ',crystal_length*1.0e6
coeffs=svdfit(crystal_mass,crystal_length,m,weight=weight,YFIT=yf)
abc=findgen(100)*max(crystal_mass)/99.
ord=poly(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,'length(M) coeffs = ',coeffs
plot_oi,crystal_mass*1.0e9,crystal_length*1.0e6, $
	tit='!5 crystal mass vs. crystal length parameterization', $
        xtit='Crystal Mass nanograms', $
        ytit='Crystal Length microns', $
        linestyle=1, thick=2.0,psym=6,/ynoz
oplot,abc*1.0e9,ord*1.0e6, $
        linestyle=0,thick=2.
oplot,crystal_mass*1.0e9,yf*1.0e6,linestyle=0,thick=2.0,psym=2

end

pro area_tester
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr

area=area_of_length(crystal_length)
	plot_oo,crystal_length(*)*1.0e6,area(*)*1.0e12, $
		tit='!5Area of Cirrostratus Ice Crystals', $
		xtit='Length (!7l!5)', xsty=0, $
		ytit='!5Crystal Area (!7l!5!E2!N)',ysty=0

crystal_length=[20,50,120,300,750]
crystal_length=crystal_length*1.0e-6
crystal_diameter=[20,40,60,100,160]
crystal_diameter=crystal_diameter*1.0e-6
area=liou_area_of_length_diameter(crystal_length,crystal_diameter)
	oplot,crystal_length(*)*1.0e6,area(*)*1.0e12, $
		linestyle=2

end

pro length_tester
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr
prism_radius=prism_radius_of_length(crystal_length)
equiv_radius=length_to_equiv_rad(crystal_length)
;
plot,crystal_length(*)*1.0e6,equiv_radius(*)*1.0e6, $
	tit='!5Length Scales of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	ytit='!5Equiv Radius r!IES!N (solid) Prism Radius (dotted) (!7l!5)', $
	ysty=0
oplot,crystal_length(*)*1.0e6,prism_radius(*)*1.0e6,linestyle=1

crystal_length=[20,50,120,300,750]
crystal_length=crystal_length*1.0e-6
crystal_diameter=[20,40,60,100,160]
crystal_diameter=crystal_diameter*1.0e-6
liou_equiv_radius=liou_equiv_rad_of_length_diameter(crystal_length,crystal_diameter)
	oplot,crystal_length(*)*1.0e6,liou_equiv_radius(*)*1.0e6, $
		linestyle=2
;overplot the a, the largest radius of the basal plane
oplot,crystal_length(*)*1.0e6,.5*crystal_diameter(*)*1.0e6,linestyle=3

print,'crystal_length = ',crystal_length
print,'crystal_diameter = ',crystal_diameter
print,'liou_equiv_radius = ',liou_equiv_radius

end

pro capacitance_tester
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr

capacitance=capacitance_of_length(crystal_length)
capacitance2=crystal_length/2.
equiv_rad=length_to_equiv_rad(crystal_length)
capacitance3=equiv_rad
;
plot,crystal_length(*)*1.0e6,capacitance(*)*1.0e6, $
	tit='!5Capacitance Shape Factors of Cirrostratus Ice Crystals', $
	xtit='Length (!7l!5)', xsty=0, $
	ytit='!5Capacitance C (P & K--solid) (Heymsfield--dashed) (!7l!5)', $
	ysty=0
oplot,crystal_length(*)*1.0e6,capacitance2(*)*1.0e6,linestyle=2
oplot,crystal_length(*)*1.0e6,capacitance3(*)*1.0e6,linestyle=3
;print,capacitance,capacitance2
end

pro IWC_tester
num_layer = 100
;cloud_base = 1.0e4	;m
;cloud_thick = 1000	;m
;abc_nbr=100
;moist_lapse_rate =0.0065		;K/m = 6.5 K/km
;Tnot_bot = 228	;Kelvin
;dz=float(cloud_thick)/num_layer
;layer_altitude=fltarr(num_layer)
;for layer=0,num_layer-1 do begin
;	layer_altitude(layer)=cloud_base+layer*dz	;m
;endfor
env_temperature=180.+70.*findgen(num_layer)/(num_layer-1)
;
IWC=IWC_of_temp(env_temperature)
;
plot,IWC(*)*1.0e6,env_temperature(*), $
	tit='!5IWC of Cirrostratus vs Temperature', $
	xtit='!5IWC C (Liou--solid) (mg/m!E3!N)', $
	ytit='Temperature (!7l!5)', xsty=0, $
	ysty=1,/ynoz
end

pro eqm_vap_tester
abc_nbr=100
min_temp=203.15
max_temp=273.15
env_temperature=(findgen(abc_nbr)+1.)*(max_temp-min_temp)/ $
		abc_nbr+min_temp
env_temp_celsius=env_temperature-273.15
;
water=eqm_vap_water(env_temperature)
ice=eqm_vap_ice(env_temperature)
;
plot,env_temp_celsius(*),water(*)/100.0, $
	tit='!5Eq!Em!N Vapor Pressure for Liquid, Ice Water', $
	xtit='!5Environmental Temperature (!12_!5C)', xsty=0, $
	ytit='!5Pressure (water--solid) (ice--dashed) (mbar = hPa)', $
	ysty=0, $ 
	charsize=1.5, thick = 2.
oplot,env_temp_celsius(*),ice(*)/100.0,linestyle=2,thick=2.

   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	lpr
   endif	

plot,env_temp_celsius(*),ice(*)/water(*) <  1.0, $
	tit='!5Ratio of Ice to Liquid Water Eq!Em!N Vapor Pressure', $
	xtit='!5Environmental Temperature (!12_!5C)', xsty=0, $
	ytit='!5e!Isat,ice!N/e!Isat,liquid!N', $
	ysty=0, $
	charsize=1.5, thick = 2.

   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	lpr
   endif	

min_saturation_liquid = .2
max_saturation_liquid = 1.
temperature = 273.15
saturation_liquid = $
  ((findgen(abc_nbr)+1.)*max_saturation_liquid*(max_saturation_liquid-min_saturation_liquid)/abc_nbr)+min_saturation_liquid
vapor_pressure = saturation_liquid*eqm_vap_water(temperature)
e_sat_vap_ice = eqm_vap_ice(temperature)
saturation_ice = vapor_pressure/e_sat_vap_ice

   plot, saturation_liquid, saturation_ice, $
	tit='!5Ice vs. Liquid Water Saturations', $
	xtit='S!Iliquid!N', xsty=0, $
	ytit='!5S!Iice!N', $
	ysty=0,  yrange = [.75,1.75], $
	charsize=1.5, thick = 2.
   
;compute the rest of the curves
FOR temp = 1, 18 DO BEGIN
    temperature = temperature-5.
    vapor_pressure = saturation_liquid*eqm_vap_water(temperature)
    e_sat_vap_ice = eqm_vap_ice(temperature)
    saturation_ice = vapor_pressure/e_sat_vap_ice

   oplot, saturation_liquid, saturation_ice, $
     linestyle = temp MOD 6, thick=2.
ENDFOR
   
;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

function multi_w_fall_speed,crystal_length,crystal_mass, $
	concentration,env_pressure,num_layer,num_bins
; 
; returns the mass-weighted fall speed for an entire cloud, computes
; all the levels at one time and so accepts 2-d array args.
; according the method in Zhang p.313.
;
; REMEMBER: output in meters/s
;
fall_speed=fltarr(num_bins)
mass_of_layer=fltarr(num_layer)
mass_weighted_speed=fltarr(num_layer)
weighted_speed=fltarr(num_layer)
;
for layer=0,num_layer-1 do begin
fall_speed(*)=fall_speed_simple(crystal_length(layer,*), $
	env_pressure(layer))			;m/s
mass_of_layer(layer)=total(concentration(layer,*)*$
	crystal_mass(layer,*))			;kg/m^3
mass_weighted_speed(layer)=total(concentration(layer,*)*$
	crystal_mass(layer,*)*fall_speed(*))	
endfor
;
weighted_speed(*)=mass_weighted_speed(*)/mass_of_layer(*)
;
return,weighted_speed			;m/s
end

pro multi_w_tester
;
; attempts to graph the fall speeds of a realistic initial distribution
; of ice bullets vs. altitude
;
num_layer = 10
num_bins=100
max_length=1000.	;microns
cloud_base = 1.0e4       ;m
cloud_thick = 1000      ;m
moist_lapse_rate =0.0065         ;K/m = 6.5 K/km
Tnot_bot = 225.  ;Kelvin
Pnot_bot = 250.  ;mb
;
dz=float(cloud_thick)/num_layer
bin_sz_mcr = max_length/num_bins
;
crystal_length=fltarr(num_layer,num_bins)
env_pressure=fltarr(num_layer)		;Pa
env_temperature=fltarr(num_layer)
concentration=fltarr(num_layer,num_bins)
crystal_mass=fltarr(num_layer,num_bins)
layer_altitude=fltarr(num_layer)
;
for layer=0,num_layer-1 do begin
        layer_altitude(layer)=cloud_base+layer*dz       ;m
	crystal_length(layer,*)=(findgen(num_bins)+1.)* $
		max_length*1.0e-6/num_bins
endfor
env_temperature=Tnot_bot-(layer_altitude-cloud_base)*moist_lapse_rate;K
env_pressure(0:num_layer-1)=Pnot_bot*100		;mb -> Pa
;
for layer=0,num_layer-1 do begin
	crystal_mass(layer,*)=mass_of_length(crystal_length(layer,*))
	concentration(layer,*)=bin_sz_mcr*$
		sz_dst_of_IWC(crystal_length(layer,*),$
		env_temperature(layer))
endfor
;
fall_speed=multi_w_fall_speed(crystal_length,crystal_mass, $
	concentration,env_pressure,num_layer,num_bins)
;
	plot,fall_speed(*)*100.0,layer_altitude/1000.0, $
		tit='!5Fall Speed of Cirrostratus Ice Bullets', $
		xtit='!5Terminal Velocity V (cm/s)',xsty=0, $
		ytit='Altitude !8z!5 (km)', ysty=1
;
print,'altitude (km): '
print,layer_altitude/1000.
print,'fall speed (cm/s): '
print,fall_speed(*)*100.
;
end

function fall_speed_of_length,crystal_length,env_temperature, $
	env_pressure,env_density
; returns the fall speed of the bullet shaped ice crystal
; via the parameterization given in Heymsfield '72 and employed
; in both Ramaswamy's and Zhang's models.
; REMEMBER: input in meters, output in meters/s
; REMEMBER: am assuming that rho_f is env. air density.
;
num_size=n_elements(crystal_length)
F=fltarr(num_size)
X=fltarr(num_size)
W_over_L=fltarr(num_size)
G=-1.10114
B=1.05687
C=-0.09244
H=0.00535
ln10=alog(10.)
length_mm=crystal_length*1000.	;m --> mm
p_mb=env_pressure/100.		;Pa --> hPa = mb
t_celsius=env_temperature-273.  ;K --> C
;
rho_s=.78*length_mm^(-0.0038)		;gm/cm^3
;rho_s=0.8
rho_f=env_density*1.0e3/1.0e6		;kg/m^3 --> g/cm^3
;
mu=4.301e-2*((t_celsius+273.)^2.5)/(p_mb*(t_celsius+273.+120.))
;
for i=0,num_size-1 do begin
	if length_mm(i) le .3 then begin
		F(i)=mu/(2.5e-2*length_mm(i)^.7856)
		X(i)=2.34e-2*rho_s(i)*(length_mm(i)^2.304)/(mu*mu*rho_f)
		W_over_L(i)=.25*(length_mm(i)^(-.214))
	endif else begin
		F(i)=mu/(1.85e-2*length_mm(i)^.532)
		X(i)=8.77e-3*rho_s(i)*(length_mm(i)^1.475)/(mu*mu*rho_f)
		W_over_L(i)=0.185*(length_mm(i)^(-.486))
	endelse
endfor

delta=exp(ln10*(-0.090186+1.0034*alog10(X)-.10142*((alog10(X))^2.) + $
	.0083*((alog10(X))^3.))) - exp(ln10*(-1.10114+1.05687*alog10(X) - $
	.09244*((alog10(X))^2.)+0.00535*((alog10(X))^3.)))* $
	(-2.56*sqrt(W_over_L)+1.81)
U=G+B*alog10(X)+C*(alog10(X)^2.)+H*(alog10(X)^3.)
;U=F*(exp(ln10*U)+delta)		;cm/s
U=F*exp(ln10*U)		;cm/s
;
print,'crystal length (microns): '
print,crystal_length*1.0e6
print,'fall speed (cm/s): '
print,U
;
U=U/100.	;m/s
return,U			
end

pro fall_speed_tester
;
; NB:  These test #'s correspond to Heymsfield (1972) Fig. 6a.
;
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr
env_temperature=233	;K
env_pressure=25000	;Pa
;env_temperature=233	;K
;env_pressure=40000	;Pa
gas_const_dry_air=287.05 ;J/kg/K
env_density=env_pressure/(env_temperature*gas_const_dry_air)   	;kg/m^3

fall_speed=fall_speed_of_length(crystal_length,env_temperature, $
	env_pressure,env_density)
;fall_speed=fall_speed_simple(crystal_length,env_pressure)
;
	plot,crystal_length(*)*1.0e6,fall_speed(*)*100.0, $
		tit='!5Fall Speed of Cirrostratus Ice Bullets', $
		xtit='Length (!7l!5)', xsty=0, $
		ytit='!5Terminal Velocity V (cm/s)',ysty=0
end

pro weighted_tester
;
; tests the weighted_fall_speed routine, AOT the multi-weighted or
; the unweighted
;
abc_nbr=100
max_length=1000.	;microns
crystal_length=(findgen(abc_nbr)+1.)*max_length*1.0e-6/abc_nbr
crystal_mass=mass_of_length(crystal_length(*))
env_pressure=25000.	;Pa
env_temperature=228.	;K
bin_sz_mcr=abc_nbr/max_length
;
concentration=bin_sz_mcr*$
	sz_dst_of_IWC(crystal_length(*),$
	env_temperature)
weighted_speed=weighted2_fall_speed(crystal_length,crystal_mass, $
	concentration,env_pressure)
;
	plot,crystal_length(*)*1.0e6,weighted_speed(*)*100.0, $
		tit='!5Fall Speed of Cirrostratus Ice Bullets', $
		xtit='Length (!7l!5)', xsty=0, $
		ytit='!5Terminal Velocity V (cm/s)',ysty=0
end

function D_rogers,temperature,pressure
; looks up the temperature-pressure dependent coefficient of diffusion
; of water vapor in air, D, from value given in Rogers & Yau p. 103.
; unit of D are m^2/s
; REMEMBER to input temperature in Kelvin!
; REMEMBER to input pressure in mbar!
;
t=temperature-273.15		;convert to celsius
D=t				;dimensions the array if necc.
n=n_elements(t)
pressure_mb=pressure/100.
if n gt 1 then begin
	for i = 0,n-1 do begin
		if t(i) gt 25 then D(i) = 2.69e-5 else $
		if t(i) gt 15 then D(i) = 2.52e-5 else $
		if t(i) gt 5 then D(i) = 2.36e-5 else $
		if t(i) gt -5 then D(i) = 2.21e-5 else $
		if t(i) gt -15 then D(i) = 2.06e-5 else $
		if t(i) gt -25 then D(i) = 1.91e-5 else $
		if t(i) gt -35 then D(i) = 1.76e-5 else $
		D(i) = 1.62e-5 
	endfor
endif else begin
	if t gt 25 then D = 2.69e-5 else $
	if t gt 15 then D = 2.52e-5 else $
	if t gt 5 then D = 2.36e-5 else $
	if t gt -5 then D = 2.21e-5 else $
	if t gt -15 then D = 2.06e-5 else $
	if t gt -25 then D = 1.91e-5 else $
	if t gt -35 then D = 1.76e-5 else $
	D = 1.62e-5 
endelse
;
; The above values are for p = 100 kPa = 1 b = 1000 mb and so we must
; scale by 100/pressure (kPa) or 10/pressure (mb)
;
D=D*1000./pressure_mb
return,D
end

pro D_tester
temperature=findgen(100)/2.+223.15
t_celsius=temperature-273.15
;
D_pruppacher=D_pruppacher(temperature=temperature,env_pressure=env_pressure)
D_rogers=D_rogers(temperature,env_pressure)
;
;plot,temperature,D_pruppacher,xtitle='Temperature (K)',ytitle='K (J/m/s/C)', $
plot,temperature,D_rogers,xtitle='Temperature (K)',ytitle='D (m!E2!N/s)', $
	title='Vapor Diffusivity D!Iv!N'
;oplot,temperature,D_rogers,linestyle=1
oplot,temperature,D_pruppacher,linestyle=1
end

