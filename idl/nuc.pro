;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: nuc.pro,v $
; $Source: /home/zender/cvs/idl/nuc.pro,v $
; $Id$
; $State: Exp $
;
; Purpose: gather all the nucleation specific routines in one place.
; this file still relies on some generic routines in clouds.pro
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-10 23:36:35  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/01 01:55:52  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 02:09:43  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:29  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:46  zender
; Imported sources
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
function act_energy_of_T,temperature
;PrK78 p. 178
;Actual pzn. only recommended for -30 < T < 0
avacados_nbr=6.022045e23	;#/mol^-1
JOULES_PER_KCAL = 4.1888*1.0e3

t_celsius=temperature-273.15
if ((min(t_celsius) lt -100.) or (max(t_celsius) gt 0.)) then begin
	print,"Bogus temperature in act_energy routine"
endif

A_not=5.5
A_sub_1=-1.330e-2
A_sub_2=2.74e-4
A_sub_3=1.085e-6

act_energy=A_not*exp(A_sub_1*t_celsius + A_sub_2*t_celsius*t_celsius + $
			A_sub_3*t_celsius^3)
;unit are now in kcal/mol
act_energy_mks=act_energy*JOULES_PER_KCAL
;unit are now in J/mol
act_energy_mks=act_energy_mks/avacados_nbr
;unit are now in J/molecule

return,act_energy_mks
end

function surf_tension_of_T,temperature
;DeR90 p. 1063, PrK78 p. 104
;PrK78 p. 130 say sigma [ergs/cm^2] =28.5 + .25*T [C]
t_celsius=temperature-273.15
if ((min(t_celsius) lt -50.) or (max(t_celsius) gt -30.)) then begin
	print,"Bogus temperature in surf_tension routine"
endif

A_not=28.6
A_sub_1=0.167

surf_tension_cgs=A_not + A_sub_1*t_celsius
;unit are now in ergs/cm^2
surf_tension_mks=surf_tension_cgs*1.0e-7*1.0e4
;unit are now in J/m^2

return,surf_tension_mks
end

function latent_heat_of_T,temperature
;PrK78 p. 89
;
t_celsius=temperature-273.15
if ((min(t_celsius) lt -50.) or (max(t_celsius) gt 0.)) then begin
	print,"Bogus temperature in latent_heat_fusion routine"
endif

A_not=79.7
A_sub_1=.485
A_sub_2=-2.5e-3
JOULES_PER_IT_CAL=4.1868

latent_heat_cal=A_not + A_sub_1*temperature + A_sub_2*temperature*temperature
;unit are now in IT cal/g
latent_heat_mks=latent_heat_cal*JOULES_PER_IT_CAL*1.0e3
;unit are now in J/kg

return,latent_heat_mks
end

function ice_dens_of_T,temperature
;PrK78 p. 86
t_celsius=temperature-273.15
if ((min(t_celsius) lt -150.) or (max(t_celsius) gt 0.)) then begin
	print,"Bogus temperature in ice_density routine"
endif

A_not=0.9167
A_sub_1=-1.75e-4
A_sub_2=-5.e-7

ice_dens_cgs=A_not + A_sub_1*t_celsius + A_sub_2*t_celsius*t_celsius
ice_dens_mks=ice_dens_cgs*1.0e3

return,ice_dens_mks
end

function liquid_dens_of_T,temperature
;PrK78 p. 86
t_celsius=temperature-273.15
if ((min(t_celsius) lt -50.) or (max(t_celsius) gt 0.)) then begin
	print,"Bogus temperature in liquid_density routine"
endif

A_not=0.99984
A_sub_1=0.860e-4
A_sub_2=-0.108e-4

liquid_dens_cgs=A_not + A_sub_1*t_celsius + A_sub_2*t_celsius*t_celsius
liquid_dens_mks=liquid_dens_cgs*1.0e3

return,liquid_dens_mks
end
;
function CCN_conc_HeS90,saturation
;implements HeS90 CCN concentration observations
  CCN_conc_cgs=200.*saturation^1.5	
  CCN_conc=CCN_conc_cgs*1.0e6
  return,CCN_conc
end
;
function CCN_PDF_DMC94,diameter,scaling_diameter
;implements DMC94 CCN concentration observations
  CCN_PDF=(1./scaling_diameter)*exp(-diameter/scaling_diameter)
  return,CCN_PDF
end

FUNCTION haze_freezing_param, aerosol_density, saturation_liquid, aerosol_radius, $
                            temperature
;
; Nota Bene: most of the DMC pzns. are written in CGS unit.
;
t_celsius=temperature-273.16
aerosol_rad_cgs=aerosol_radius*100.
aerosol_diam_cgs = 2.*aerosol_rad_cgs
aerosol_den_cgs=aerosol_density/1000.

;
;implements DMC92 pzn.
;aerosol_diameter = 2.*aerosol_radius
;c1 = exp(-23.2-(2.12*t_celsius))
;c2 = 3110.+94.6*t_celsius+.8*(t_celsius^2)
;a = c1*saturation_liquid^c2
;exponent = -a*((!Pi*aerosol_den_cgs/6.)^2)*aerosol_diam_cgs^5
;fraction = 1.-exp(exponent)

;implements DMC94 pzn.
c1 = -14.65 - (1.045*t_celsius)
c2 = -492.35- (8.34*t_celsius) - (.0061*t_celsius*t_celsius)
exponent = double(c1 + c2*(1.-saturation_liquid))
foo = exp(exponent*alog(10.))
a = !pi*!pi*aerosol_den_cgs*aerosol_den_cgs*foo/6.
b = 6.
exponent = -a*(aerosol_diam_cgs^b)
fraction = float(1.-exp(exponent))

;unit of fraction are s^-1
return, fraction
END

FUNCTION CN_pdf_jtwk92a, radius
  N_not_cgs=300.
  sigma=2.3
  mode_radius_mcr=.2
 
  N_not=N_not_cgs*1.0e6
  radius_mcr=radius*1.0e6

  exponent=alog(radius_mcr/mode_radius_mcr)
  exponent=exponent/(sqrt(2.)*sigma)
  exponent=-exponent*exponent
  prefactor=1./(sigma*radius*sqrt(2.*!pi))

  fraction_PDF=prefactor*exp(exponent)

  return,fraction_PDF
end

function J_rate_of_T,temperature,surf_tension_i_w
liquid_monomer_conc=5.85e18	;#/m^2 PrK78 pp. 177-178, HeS89 p. 2254
boltzmann_constant=1.38063e-23	;J/K  
planck_constant=6.62620e-34     ;Js   
molecular_wgt_water=18.0160	;g/mol and amu
;activation_energy=1.29e-19	;J JTWKH92a p. 11
activation_energy=3.4e-20	;J DeR90 p. 1063, PrK78 p. 178
latent_heat_fusion=.334e6       ;J/kg (latent heat of fusion) at 0 C
ice_water_density=0.92e3		;kg/m^3
liquid_water_density=0.985e3	;kg/m^3
;surf_tension_i_w=0.0272          ;J/m^2 = N/m DeR90 p. 1063, PrK78 p. 121 
avocados_nbr=6.022045e23	;#/mol^-1
KILOGRAMS_PER_AMU=1.660e-27	;kg/amu
JOULES_PER_KCAL = 4.1888*1.0e3
JOULES_PER_IT_CAL=4.1868
mass_water_molecule = molecular_wgt_water*KILOGRAMS_PER_AMU

;The following quantities may be either calculated with a temperature
;dependence, or just use the fixed values above from various parts of
;the literature
activation_energy=act_energy_of_T(temperature)
ice_water_density=ice_dens_of_T(temperature)
liquid_water_density=liquid_dens_of_T(temperature)
latent_heat_fusion=latent_heat_of_T(temperature)
;surf_tension_i_w=surf_tension_of_T(temperature)
liquid_monomer_conc=(liquid_water_density/mass_water_molecule)^(2./3.)

;This term is explained in PrK78 p. 152, DeR90 p. 1062, DMC92 p. 493
;T_not is the bulk ice freezing point, this is discussed in PrK78 pp.92-94
;EJ verifies that, and that T_e is env temperature
ln_T_not_over_T_e=alog(273.15/temperature)

;germ_radius is explained in PrK p. 167 and i'll use eqn. 7-27
ice_germ_radius=2.*surf_tension_i_w/ $ 
	(latent_heat_fusion*ice_water_density*ln_T_not_over_T_e)
ice_germ_molecules=4.*!pi*ice_germ_radius^3*ice_water_density/ $
			(3.*mass_water_molecule)

;should the activation energy be multiplied by the number of molecules in
;the ice germ? apparently not....
;activation_energy=activation_energy*ice_germ_molecules

;i'll use PrK78 eq. 7-26 to evaluate the curvature effect
exponent1_num=-4.*!pi*surf_tension_i_w*ice_germ_radius*ice_germ_radius/3.
exponent1_denom=boltzmann_constant*temperature
exponent1=exponent1_num/exponent1_denom

;alternative method
;exponent1_num=-16.*!pi*surf_tension_i_w^3
;exponent1_denom=ln_T_over_T_not*ice_water_density*latent_heat_fusion
;exponent1_denom=3.*boltzmann_constant*temperature*exponent1_denom*exponent1_denom
;exponent1=exponent1_num/exponent1_denom

;activation effect
exponent2=-activation_energy/(boltzmann_constant*temperature)

prefactor1=liquid_water_density*boltzmann_constant*temperature/ $
		(ice_water_density*planck_constant)

prefactor2=sqrt(surf_tension_i_w/(boltzmann_constant*temperature))

j_rate=2.*prefactor1*prefactor2*exp(exponent1)
j_rate=j_rate*liquid_monomer_conc*exp(exponent2)

;print,"temperature = ",temperature
;print,"ice_water_density = ",ice_water_density
;print,"liquid_water_density = ",liquid_water_density
;print,"latent_heat_fusion = ",latent_heat_fusion
;print,"surf_tension_i_w = ",surf_tension_i_w
;print,"activation_energy = ",activation_energy
;print,"ice_germ_radius = ",ice_germ_radius
;print,"ln_T_not_over_T_e = ",ln_T_not_over_T_e
;print,"ice_germ_molecules = ",ice_germ_molecules
;print,"liquid_monomer_conc = ",liquid_monomer_conc
;print,"prefactor1 = ",prefactor1
;print,"prefactor2 = ",prefactor2
;print,"exponent1 = ",exponent1
;print,"exponent2 = ",exponent2

return,j_rate
end

function radius_of_S_m_T,saturation,solute_mass,temperature
;computes equilibrium radius of a deliquescing but unactivated solution
;of solute mass m in a saturation S with respect to liquid water.
;see, e.g., PrK p. 142 eqn. 6-26b, or HeS89 p. 2256 eqn. 14
;the independent variable is the diameter
;uses Newton-Raphson technique because the relation is non linear
;
;Caveat:  remember this is only for droplets which are in stable eqm.,
;droplets which have been activated are unstable and will grow forever
;and never converge with this method.
;
;Therefore it would be nice to implement some sort of a check, maybe based
;on pure water, which could identify droplets which were activated before 
;grinding through the Newton-Raphson routine only to generate an error.
;
molec_wgt_water=18.0160		;g/mol and amu
surf_tension_s_a=0.083           ;J/m^2 = N/m DeR90 p. 1063, PrK78 p. 121, see also CRC 70th ed. p. F-33 for varying sigma vs. various concentrations 
gas_const_vapor=461.51          ;J/kg/K
universal_gas_const=8.314	;J/mol/K
molec_wgt_amm_sulphate=132.14	;g/mol and amu ammonium sulphate
molec_wgt_sulfuric_acid=98.08	;g/mol and amu sulfuric acid
;solution_density=1000.		;kg/m^3 should vary with solute concentration
;solute_density=1.769e3		;kg/m^3 nail these down someday!
liquid_water_density=0.985e3	;kg/m^3
osmotic_coeff=0.65		;JTWKH92a p. 8, Low69, PrK78 p. 83 table 4-2
nu_nbr_of_ions=3.		;PrK78 p. 83 table 4-2

solution_density=1700.
solute_density=1.769e3
liquid_water_density=liquid_dens_of_T(temperature)
solute_volume=solute_mass/solute_density
ln_S_exact=alog(saturation)

; first find the critical radius for the specified conditions
; formula are described on PrK78 p. 142
;A_term
;B_term
;numerator=9.*osmotic_coeff*nu_nbr_of_ions*solute_mass*

radius=1.0e-7
count=0
tolerance=1.

restart_loop: foo=1

while((abs(tolerance) gt 1.0e-4) and (count lt 25)) do begin
    print,"count = ",count," radius_u = ",radius*1.0e6," tolerance = ",tolerance
    count=count+1
;    water_mass=((4.*!pi*(radius^3)/3.)-solute_volume)*liquid_water_density
;    solution_mass=water_mass+solute_mass
;    perc_solute_by_weight=100.*solute_mass/solution_mass
;    solution_density=1000.*(1.0012+(1.2277-1.0012)*perc_solute_by_weight)/(40.-.5)
;    print,"perc_solute_by_wgt = ",perc_solute_by_wgt
;    if(perc_solute_by_weight gt 100.) then begin
;        count=1
;        radius=radius*10.
;	goto,restart_loop
;    endif
    solution_density=liquid_water_density
    term_1=2.*surf_tension_s_a/ $
		(gas_const_vapor*temperature*radius*liquid_water_density)
    term_2_num=-nu_nbr_of_ions*osmotic_coeff*solute_mass*molec_wgt_water/ $
		molec_wgt_amm_sulphate
    beta=4.*!pi*solution_density/3.
    term_2_denom=-solute_mass+beta*radius^3
    term_2=term_2_num/term_2_denom
    ln_S_approx=term_1+term_2
    tolerance=(ln_S_exact-ln_S_approx)/ln_S_exact;
    f_of_r=ln_S_approx-ln_S_exact
    df_dr_1=-term_1/radius
    df_dr_2=-3.*beta*radius*radius*term_2_num/(term_2_denom*term_2_denom)
    df_dr=df_dr_1+df_dr_2
    radius=radius-f_of_r/df_dr;
    if(radius lt 0.) then begin
        count=1
        radius=-radius
	goto,restart_loop
    endif
endwhile  
return,radius
end

pro haze_radius_tester, pause=pause, abc_nbr=abc_nbr, $
	temperature=temperature, solute_mass=solute_mass, $
	saturation=saturation
;
;haze_radius_tester,temperature=233.0,solute_mass=1.0e-18,saturation=0.95
;
if n_elements(foo) eq 0 then foo = 5
if n_elements(temperature) eq 0 then temperature = 233.
if n_elements(solute_mass) eq 0 then solute_mass = 1.0e-18
if n_elements(saturation) eq 0 then saturation = .95
if n_elements(abc_nbr) eq 0 then abc_nbr = 10
if n_elements(pause) eq 0 then pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

radius=radius_of_S_m_T(saturation,solute_mass,temperature)

saturation_liquid=saturation
vapor_pressure = saturation_liquid*eqm_vap_water(temperature)
e_sat_vap_ice = eqm_vap_ice(temperature)
saturation_ice = vapor_pressure/e_sat_vap_ice

print,"saturation w/r/t liquid water = ",saturation
print,"saturation w/r/t ice water = ",saturation_ice
print,"solute_mass (g) = ",solute_mass*1000.
print,"temperature (C) = ",temperature-273.15
print,"radius (microns) = ",radius*1.0e6
;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

pro J_tester, pause=pause, abc_nbr=abc_nbr, $
	num_sigmas=num_sigmas
if n_elements(num_sigmas) eq 0 then num_sigmas = 5
if n_elements(abc_nbr) eq 0 then abc_nbr = 10
if n_elements(pause) eq 0 then pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

min_temp=223.16
max_temp=243.16
temperature=(findgen(abc_nbr)+1.)*(max_temp-min_temp)/ $
		abc_nbr+min_temp
temperature_celsius=temperature-273.15

surf_tension_corr=fltarr(num_sigmas)
surf_tension_i_w=fltarr(num_sigmas,abc_nbr)
surf_tension_best_guess=surf_tension_of_T(temperature)
;surf_tension_best_guess=replicate(.0185,abc_nbr) ;J/m^2 = N/m DeR90 p. 1063, PrK78 p. 121, .0185 seems to match SaD88 the best
J_rate=fltarr(num_sigmas,abc_nbr)

FOR sigma = 0, num_sigmas-1 DO BEGIN
	surf_tension_corr(sigma)=0.05*(-1)^sigma*fix((sigma+1)/2)
	surf_tension_i_w(sigma,*)= $
		surf_tension_best_guess*(1.+surf_tension_corr(sigma))
	J_rate(sigma,*)=J_rate_of_T(temperature,surf_tension_i_w(sigma,*))
ENDFOR

min_y = min(J_rate)
max_y = max(J_rate)
;plot the first curve
plot_io,temperature_celsius(*),J_rate(0,*)/1.0e6, $
	tit='!5Homogeneous Freezing Rate of Pure Water', $
	xtit='!5Temperature (!12_!5C)', xsty=1, $
	ytit='!5Freezing Rate !8J!5 (#-cm!E-3!N-s!E-1!N)', $
	ysty=1, /ynoz, $ 
	yrange=[1.0e4,1.0e11], $
	charsize=1.5, thick = 2.

t_celsius=temperature-273.15
foo=	-606.3952 - 52.6611*t_celsius - 1.7439*t_celsius*t_celsius $
	- .0265*(t_celsius^3) - 1.536e-4*(t_celsius^4)
;note this pzn. gives J in #/cm^3/s
J_rate_HeM93=10.^foo
oplot,temperature_celsius(*),J_rate_HeM93(*), $
	linestyle=3, thick = 2.

;plot the rest of the curves
FOR sigma = 1, num_sigmas-1 DO BEGIN
	oplot,temperature_celsius(*),J_rate(sigma,*)/1.0e6, $
      linestyle = sigma MOD 6, thick =  2.	
ENDFOR
    
ln_lgn_x1=.22
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_size=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(num_sigmas+10)*lgn_dy

FOR sigma = 0, num_sigmas-1 DO BEGIN
    sigma_lgn = strarr(2, num_sigmas)
    sigma_lgn(sigma) = string(Format = '(I2)',abs(surf_tension_corr(sigma))*100.)
    if surf_tension_corr(sigma) lt 0. then $
	    sigma_lgn(sigma) = '- '+sigma_lgn(sigma) else $
	    sigma_lgn(sigma) = '+ '+sigma_lgn(sigma) 
    plots, [ln_lgn_x1, ln_lgn_x2], lgn_y(sigma)+0.013, $
      linestyle = sigma MOD 6, thick = 2.0, /NORMAL
    xyouts, txt_lgn_x, lgn_y(sigma), '!7r!5!Ii/w!N ' + $
      sigma_lgn(sigma)+'!5%', size = txt_lgn_size, /NORMAL 
ENDFOR
;
;plot Sassen and Dodd's data as read from a graph by eye
   SaD88_t_celsius=[-34.1,-34.6,-35.1,-37.2]
   SaD88_J_rate=[1.5e6,3.e6,2.e7,4.e9]
   plots,SaD88_t_celsius,SaD88_J_rate,psym = 2,thick = 2.
;
 if pause eq 'y' then begin
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
 
;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

PRO CN_tester, CN_init=CN_init, pause=pause, abc_nbr=abc_nbr
IF n_elements(CN_init) EQ 0 THEN CN_init = 300.e6			
IF n_elements(abc_nbr) EQ 0 THEN abc_nbr = 100
IF n_elements(pause) EQ 0 THEN pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

min_radius = 1.0e-8           ;meters
max_radius = 1.0e-6           ;meters
CN_radius=((findgen(abc_nbr)+1.)*(max_radius-min_radius)/abc_nbr)+min_radius
delta_radius=shift(CN_radius,-1)-CN_radius
delta_radius(abc_nbr-1)=delta_radius(abc_nbr-2)
CN_PDF = CN_PDF_JTWK92a(CN_radius)
CN_dst = CN_PDF*CN_init
CN_conc = CN_dst*delta_radius

print,"total concentration (#/cm^3) = ",total(CN_conc)*1.0e-6
;print,"total concentration (#/cm^3) >= .1 microns = ", $
;	total(CN_conc(9:abc_nbr-1))*1.0e-6
;print,"total concentration (#/cm^3) 3 <= L <= 10 microns = ", $
;	total(CN_conc(2:9))*1.0e-6

plot_oo,CN_radius(*)*1.0e6,CN_dst*1.0e-6*1.0e-6, $
	tit='!5CN Distribution vs. Radius', $
	xtit='Radius (!7l!5)', xsty=0, $
	ytit='!5CN Distribution (#-cm!E-3!N-!7l!5m!E-1!N)',ysty=0

 if pause eq 'y' then begin
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
 
 plot_oo,CN_radius(*)*1.0e6,CN_conc*1.0e-6, $
	tit='!5CN Concentration vs. Radius', $
	xtit='Radius (!7l!5)', xsty=0, $
	ytit='!5CN Concentration (#-cm!E-3!N)',ysty=0

 if pause eq 'y' then begin
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
 
 dN = shift(CN_conc, -1)-CN_conc
dN(abc_nbr-1) = dN(abc_nbr-2)
dlogr =  delta_radius/CN_radius

 print, "N= ", CN_conc
 print, "dN= ", dN
 print, "r= ", CN_radius
 print, "delta r= ", delta_radius
 print, "dlogr= ", dlogr
 print, "dNdlogr= ", dN/dlogr

plot_oo, CN_radius(*)*1.0e6, dN/dlogr > 1.0e-15, $
	tit='!5Jensen figure 2', $
	xtit='Radius !8!5 (!7l!5m)', xsty=0, $
	ytit='!5dN/dLog(!8r!5)',ysty=0

;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

PRO CCN_tester, CCN_init=CCN_init, pause=pause, abc_nbr=abc_nbr, $
	scaling_diameter=scaling_diameter, solute_density=solute_density, $
	temperature=temperature,saturation_liquid=saturation_liquid
;
; Supposed to generate reasonable dry solute mass distributions useful for
; homogeneous CCN.
;
IF n_elements(CCN_init) EQ 0 THEN CCN_init = 200.e6			
IF n_elements(scaling_diameter) EQ 0 THEN scaling_diameter = 7.5e-8
IF n_elements(abc_nbr) EQ 0 THEN abc_nbr = 100
IF n_elements(solute_density) EQ 0 THEN solute_density = 1.769e3
IF n_elements(temperature) EQ 0 THEN temperature = 223.16 
IF n_elements(saturation_liquid) EQ 0 THEN saturation_liquid = .95
IF n_elements(pause) EQ 0 THEN pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

max_mass = 1.0e-15          	;kg
min_mass = max_mass/(2.^((abc_nbr-1.)/2.))
;
CCN_mass=fltarr(abc_nbr+1)
for abc=1,abc_nbr do begin
	CCN_mass(abc)=min_mass*2.^((abc-1.)/2.)
endfor
CCN_mass_LBC=CCN_mass(1)*(1.+1./sqrt(2.))/2.
CCN_mass_RBC=CCN_mass(abc_nbr)*(1.+sqrt(2.))/2.
;
CCN_int_mass=fltarr(abc_nbr+1)
;for abc=1,abc_nbr-1 do begin
;	CCN_int_mass(abc)=0.5*(CCN_mass(abc)+CCN_mass(abc+1))
;endfor
CCN_int_mass=0.5*(CCN_mass+shift(CCN_mass,-1))
CCN_int_mass(0)=CCN_mass_LBC
CCN_int_mass(abc_nbr)=CCN_mass_RBC
;
delta_mass=fltarr(abc_nbr+1)
;for abc=1,abc_nbr do begin
;	delta_mass(abc)=CCN_int_mass(abc)-CCN_int_mass(abc-1)
;endfor
delta_mass=CCN_int_mass-shift(CCN_int_mass,1)
;
CCN_radius = 3.*CCN_mass/(4.*!pi*solute_density)
CCN_radius = CCN_radius^(1./3.)
;
CCN_int_radius = 3.*CCN_int_mass/(4.*!pi*solute_density)
CCN_int_radius = CCN_int_radius^(1./3.)
;
delta_radius=fltarr(abc_nbr+1)
;for abc=1,abc_nbr do begin
;	delta_radius(abc)=CCN_int_radius(abc)-CCN_int_radius(abc-1)
;endfor
delta_radius=CCN_int_radius-shift(CCN_int_radius,1)
;
CCN_diameter=CCN_radius*2.
delta_diameter=delta_radius*2.
;
CCN_PDF = CCN_PDF_DMC94(CCN_diameter,scaling_diameter)
CCN_dst = CCN_PDF*CCN_init
CCN_conc = CCN_dst*delta_diameter
;
frac_haze_freezing = haze_freezing_param( $
	solute_density, saturation_liquid, $
	CCN_radius, temperature) 
;
;if saturation_liquid lt .82 then begin
;	frac_haze_freezing=0.*frac_haze_freezing
;endif
;
CCN_conc_freezing=frac_haze_freezing*CCN_conc
;
vapor_pressure = saturation_liquid*eqm_vap_water(temperature)
e_sat_vap_ice = eqm_vap_ice(temperature)
saturation_ice = vapor_pressure/e_sat_vap_ice
tot_CCN_conc = total(CCN_conc(1:abc_nbr))
tot_CCN_freezing = total(CCN_conc_freezing(1:abc_nbr))
;
print,"temperature = ",temperature
print,"saturation_liquid = ",saturation_liquid
print,"saturation_ice = ",saturation_ice
print,""
;print,"CCN_conc*1.0e-6 = ",CCN_conc*1.0e-6
;print,"CCN_mass*1.0e3 = ",CCN_mass*1.0e3
;print,"CCN_diameter*1.0e6 = ",CCN_diameter*1.0e6
;print,"delta_diameter*1.0e6 = ",delta_diameter*1.0e6
;print,"frac_haze_freezing = ",frac_haze_freezing
;print,"CCN_conc_freezing*1.0e6 = ",CCN_conc_freezing*1.0e6
;print,""
print,"mean diameter (microns) = scaling_diameter = ",scaling_diameter*1.0e6
print,"total concentration (#/cm^3) = ",tot_CCN_conc*1.0e-6
print,"total ice nucleated (#/cm^3/s) = ",tot_CCN_freezing*1.0e-6
;
plot_oo,CCN_diameter(1:abc_nbr)*1.0e6,CCN_dst(1:abc_nbr)*1.0e-6*1.0e-6, $
	tit='!5CCN Distribution vs. Diameter', $
	xtit='Diameter (!7l!5)', xsty=0, $
	ytit='!5CCN Distribution (#-cm!E-3!N-!7l!5m!E-1!N)',ysty=0
;	xrange=[0.1,100], $
;	yrange=[.0001,100]

 if pause eq 'y' then begin
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
 
 plot_oo,CCN_diameter(1:abc_nbr)*1.0e6,CCN_conc(1:abc_nbr)*1.0e-6, $
	tit='!5CCN Concentration vs. Diameter', $
	xtit='Diameter (!7l!5)', xsty=0, $
	ytit='!5CCN Concentration (#-cm!E-3!N)',ysty=0, $
	psym=1

 if pause eq 'y' then begin
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
 
 plot_oo,CCN_diameter(1:abc_nbr)*1.0e6,frac_haze_freezing(1:abc_nbr), $
	tit='!5Fraction Freezing vs. Diameter', $
	xtit='Diameter (!7l!5)', xsty=0, $
	ytit='!5Fraction Freezing (s!E-1!N)',ysty=0, $
	psym=1

 if pause eq 'y' then begin
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
 
 plot_oo,CCN_diameter(1:abc_nbr)*1.0e6,CCN_conc_freezing(1:abc_nbr)*1.0e-6, $
	tit='!5CCN Freezing vs. Diameter', $
	xtit='Diameter (!7l!5)', xsty=0, $
	ytit='!5CCN Freezing (#-cm!E-3!N-s!E-1!N)',ysty=0, $
	psym=1

 if pause eq 'y' then begin
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
 
;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

PRO haze_tester, min_solute_mass=min_solute_mass, pause=pause, $
	abc_nbr=abc_nbr, temperature=temperature, $
	solute_density=solute_density, num_mass=num_mass
;; essentially this test recreates figure 3. from DMC92
IF n_elements(min_solute_mass) EQ 0 THEN min_solute_mass = 1.0e-18
IF n_elements(solute_density) EQ 0 THEN solute_density = 1.769e3
;100% sulfuric acid = 1.831e3 kg/m^3 crc handbook 70th ed. pg. F6
;100% ammonium sulphate = 1.0477e3 kg/m^3 crc handbook 70th ed. pg. D225
;ammonium sulphate = 1.769e3 kg/m^3 DMC94 p. 80.
IF n_elements(temperature) EQ 0 THEN temperature = 223.16 
IF n_elements(abc_nbr) EQ 0 THEN abc_nbr = 100
IF n_elements(num_mass) EQ 0 THEN num_mass = 4
IF n_elements(pause) EQ 0 THEN pause = 'y' ;ready to roll
IF pause EQ 'y' THEN !P.multi = 0 ELSE BEGIN
    !P.multi = [0, 3, 1, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
ENDELSE

max_saturation_liquid = 1.
min_saturation_liquid = .82
saturation_liquid = ((findgen(abc_nbr)+0.)* $
	(max_saturation_liquid-min_saturation_liquid)/abc_nbr)+ $
	min_saturation_liquid 
vapor_pressure = saturation_liquid*eqm_vap_water(temperature)
e_sat_vap_ice = eqm_vap_ice(temperature)
saturation_ice = vapor_pressure/e_sat_vap_ice

solute_mass = min_solute_mass*10^(indgen(num_mass))
solute_radius = solute_mass
frac_haze_freezing = fltarr(num_mass, abc_nbr)

FOR mass = 0, num_mass-1 DO BEGIN
    solute_radius(mass) = .5*(6.*solute_mass(mass)/ $
				(!pi*solute_density))^(1./3.)
    frac_haze_freezing(mass, *) = haze_freezing_param( $
	solute_density, saturation_liquid, $
	solute_radius(mass), temperature) 
ENDFOR

min_y = min(frac_haze_freezing)
max_y = max(frac_haze_freezing)
;plot the first curve
plot_io, saturation_liquid, frac_haze_freezing(0, *),  $
  tit = '!5Liquid Saturation !8S!I!5l!N!5 vs. Fraction Haze Freezing !8F!5', $
  xtit = 'Saturation !8S!I!5liquid!N!5', $
;plot_io, saturation_ice, frac_haze_freezing(0, *),  $
;  tit = '!5Ice Saturation !8S!I!5i!N!5 vs. Fraction Haze Freezing !8F!5', $
;  xtit = 'Saturation !8S!I!5ice!N!5', $
  ytit = '!5Fraction Haze Freezing !8F!5 (s!E-1!N)', ysty = 0, $
  xsty = 0, $
  yrange = [min_y > 1.0e-8,max_y],  $
  charsize =  1.5, $
  thick = 2.

;plot the rest of the curves
FOR mass = 1, num_mass-1 DO BEGIN
    oplot, saturation_liquid, frac_haze_freezing(mass, *), $
;    oplot, saturation_ice, frac_haze_freezing(mass, *), $
      linestyle = mass MOD 6, thick =  2.
ENDFOR
    
ln_lgn_x1=.22
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_size=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(num_mass+10)*lgn_dy

FOR mass = 0, num_mass-1 DO BEGIN
    mass_lgn = strarr(7, num_mass)
    mass_lgn(mass) = string(Format = '(E7.1)', solute_mass(mass)*1.0e3)
    radius_lgn = strarr(7, num_mass)
    radius_lgn(mass) = string(Format = '(F5.2)', solute_radius(mass)*1.0e6)
    plots, [ln_lgn_x1, ln_lgn_x2], lgn_y(mass)+0.013, $
      linestyle = mass MOD 6, thick = 2.0, /NORMAL
    xyouts, txt_lgn_x, lgn_y(mass), '!8m!Is!N = ' + $
      '!5'+mass_lgn(mass)+'!5 g, !8r!Is!N = !5'+radius_lgn(mass) + $
      ' !7l!5m', size = txt_lgn_size, /NORMAL 
    
ENDFOR
;
temperature_lgn = string(Format = '(F6.1)', temperature-273.15)
xyouts, txt_lgn_x, lgn_y(num_mass), '!8T = !5'+ $
	temperature_lgn+'!5 !12_!5C!5',size = txt_lgn_size, /NORMAL 
;
density_lgn = string(Format = '(F4.2)', solute_density*1.0e-3)
xyouts, txt_lgn_x, lgn_y(num_mass+1), '!7q!I!5solute!N = !5'+ $
	density_lgn+'!5 g-cm!E-3!N',size = txt_lgn_size, /NORMAL 
;
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

;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end


