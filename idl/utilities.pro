;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: utilities.pro,v $
; $Source: /home/zender/cvs/idl/utilities.pro,v $
; $Id$
; $State: Exp $
;
; Purpose: subroutines to yield cloud properties.
;
; $Log: not supported by cvs2svn $
; Revision 1.5  2000-01-10 23:36:36  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:54  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:30  zender
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
; Revision 1.4  1993/06/07  23:45:58  zender
; a year later's update. about to get rid of the original
; clouds main program.
;
; Revision 1.3  1992/07/23  17:16:13  zender
; defaults to new graphing routine, jeff_plt, which drew the graphs
; for the site review.
;
; Revision 1.2  1992/07/13  15:07:25  zender
; equivalent to clouds.c version 1.5of 7/13/92. 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function equiv_rad_to_bullet_length,equiv_rad
;
; returns the length of the ice crystal bullet corresponding to
; the given effective (equivalent surface area sphere) radius. 
; Newton-Raphson must be used because the equation is screwy. 
;
  CM_PER_METER=100.
  METERS_PER_CM=0.01
  EMINUSFOUR=1.0e-4	
;
  sphere_area_cm=4.*!PI*equiv_rad*equiv_rad* $
    CM_PER_METER*CM_PER_METER
;
  length_cm=2.*equiv_rad*CM_PER_METER
;
;  if(debug == 16){
;    (void)fprintf(fp_err,"input equiv_rad=%f\n",equiv_rad)
;    (void)fprintf(fp_err,"sphere_area_cm2=%f\n",sphere_area_cm)
;  } /* end debug */
;
  count=1
  tolerance=1.
  while(abs(tolerance) gt EMINUSFOUR and count lt 25) do begin
    count=count+1
    surface_area_cm=0.0334*(length_cm^1.572)+.505*(length_cm^1.786)
    tolerance=(surface_area_cm-sphere_area_cm)/sphere_area_cm
    dAdl=0.0525*(length_cm^.572)+.90193*(length_cm^.786)
    fofl=surface_area_cm-sphere_area_cm
    length_cm=length_cm-(fofl/dAdl)
;    if(debug == 16){
;      (void)fprintf(fp_err,"count = %i length_cm = %f epsilon = %f\n",
;		    count,length_cm,tolerance)
;    } /* end debug */
  endwhile
;
  crystal_length=length_cm*METERS_PER_CM
  return,crystal_length
  end	;end equiv_rad_to_bullet_length()  
;
function equiv_rad_to_hex_clm_length,equiv_rad
;
;  /* returns the length of the ice crystal bullet corresponding to
;     the given effective (equivalent surface area sphere) radius. 
;     Newton-Raphson must be used because the equation is screwy. */
;
  EMINUSFOUR=1.0e-4	
  sphere_area=4.*!pi*equiv_rad*equiv_rad
;
  guess_length=2.4*equiv_rad
;
;  if(debug == 80){
;    (void)fprintf(fp_err,"input equiv_rad = %g microns\n",
;		  equiv_rad*MICRONS_PER_METER)
;    (void)fprintf(fp_err,"sphere_area = %g square microns\n",
;		  sphere_area*MICRONS_PER_METER*MICRONS_PER_METER)
;    (void)fprintf(fp_err,"guess length = %g microns\n",
;		  guess_length*MICRONS_PER_METER)
;  } /* end debug */

  count=1
  tolerance=1.
  hex_length=guess_length
  while(abs(tolerance) gt EMINUSFOUR and count lt 25) do begin
    count=count+1
    hex_Aspect_ratio,hex_length,aspect_ratio,daspect_ratiodL
    hex_diameter=hex_length/aspect_ratio
    surface_area=3.*sqrt(3.)*hex_diameter*hex_diameter/4. 
    surface_area=surface_area+3.*hex_length*hex_diameter
    aspect_ratio_squared=aspect_ratio*aspect_ratio
    fofL=surface_area-sphere_area
    tolerance=fofL/sphere_area
;    /* see NCARIIICZP#57 */ 
    dfdL=(sqrt(3.)/4.)*(aspect_ratio-daspect_ratiodL*hex_length)+ $
      aspect_ratio_squared-aspect_ratio*daspect_ratiodL/2.
    dfdL=dfdL*6.*hex_length/(aspect_ratio*aspect_ratio_squared)
    hex_length=hex_length-(fofL/dfdL)
;    if(debug == 80){
;      (void)fprintf(fp_err,"count = %i length  = %g microns, epsilon = %g\n",
;		    count,hex_length*MICRONS_PER_METER,tolerance)
;    } /* end debug */
    endwhile	

  return,hex_length
  end ;equiv_rad_to_hex_clm_length()
;
pro hex_Aspect_ratio,hex_crystal_length,aspect_ratio,daspect_ratiodL
;
;  /* Implements the AuV70 hexagonal column aspect ratio. The pzn.
;     coefficients were derived in IDL by a  
;     fit to Liou's implementation of Auer and Veal's data.
;     Also returns the rate of

;     change of aspect ratio with respect to length, because that's
;     needed by the newton-raphson routine. */ 
;
  aspect_ratio_coeffs=[5.82175,-4.99677]

  aspect_ratio=aspect_ratio_coeffs(0)+aspect_ratio_coeffs(1)* $
    exp(-hex_crystal_length/.00050)

  daspect_ratiodL=aspect_ratio_coeffs(1)*(-hex_crystal_length/.00050)* $
    exp(-hex_crystal_length/.00050)

  end ;Aspect_ratio() 
;
function PDF_of_length,length
;
;     return the (usually initial) probability size dist. of a cirrus cloud as 
;     parameterized by Dowling and Radke (JAM v.29, 1990 pp.970--978) from
;     an ensemble of many observations.
;     the length is input in meters, and the number of particles per cubic
;     meter air per meter of length is returned.  This number should then
;     be multiplied by a size domain (in meters) to obtain an absolute 
;     concentration in #/m^3. 
;
  length_mcr=length*1.0e6		; /* m -> microns */
;  float_foo=exp(-length_mcr/500.)
;  PDF=float_foo/(3.4225*(length_mcr+10.))
  PDF=(10.*exp(-(length_mcr/5.))+1.)* $
	(1./(length_mcr)^1.6)*exp(-length_mcr/500.)
;   Convert from prob. per micron of length -> prob. per meter of length 
  return,PDF*1.0e6
end

function eqm_vap_ice,temperature
;
; return the equilibrium water vapor pressure over bulk solid ice, in mbar, 
; for the given temperature in Kelvin
; REMEMBER: value is returned in Pa NOT mb
;
foo=24.29-6148./temperature
return,100.*exp(foo)		;mb -> Pa
end

function eqm_vap_water,temperature
;
; return the equilibrium water vapor pressure over bulk liquid water, in mbar, 
; for the given temperature in Kelvin
; REMEMBER: value is returned in Pa NOT mb
;
foo=21.6-5420./temperature
return,100.*exp(foo)		;mb -> Pa
end

function prism_radius_of_length,length
;
; return the radius of a bullet according to the parameterization
; of Heymsfield '72 p. 1351
; for the given length in meters
; This prism radius is the fancy way to comput the curvature according 
; to Pruppacher & Klett
; REMEMBER: this only works up to lengths of 3000 microns (3 mm)
;
length_mm=length*1000.		;m -> mm
width_mm=.25*(length_mm^.7856)	;mm
return,.5*width_mm/1000.	;mm -> m, width -> radius
end

function length_to_equiv_rad,crystal_length
;
; returns the radius of the sphere the same surface area as a bullet
; shaped ice-particle via the parameterization given in Zhang '89
; Equivalent radius is used not only in Mie calculations, but for
; the radius in foo_denom_2
; REMEMBER: input and output in meters
;
length_cm=crystal_length*100.
surface_area=0.0334*(length_cm^1.572)+.505*(length_cm^1.786)
equiv_radius=sqrt(.25*surface_area/!pi)	;cm
return,equiv_radius/100.		;m
end

function fall_speed_simple,crystal_length,env_pressure
;
; returns the fall speed of the bullet shaped ice crystal
; via the parameterization given in Heymsfield '72, the highly
; parameterized form on p. 1356. Maybe this is the form used
; in both Ramaswamy's and Zhang's models. As Heymsfield notes
; this parameterization overestimates the fall speed by ~ 25%
;
; REMEMBER: input in meters, output in meters/s
;
num_size=n_elements(crystal_length)
length_mm=crystal_length*1000.	;m --> mm
p_env_mb=env_pressure/100.	;Pa --> hPa = mb
p_1000=1000.			;mb
;
U_1000=-1.2+169.65*length_mm-84.5*length_mm*length_mm	;cm/s
U=U_1000*sqrt(p_1000/p_env_mb)
;
;print,'crystal length (microns): '
;print,crystal_length*1.0e6
;print,'fall speed (cm/s): '
;print,U
;
return,U/100.			;m/s
end

function weighted_fall_speed,crystal_length,crystal_mass, $
	concentration,env_pressure
; 
; returns the mass-weighted fall speed for a given level of cloud,
; according the method in Zhang p.313.
;
; REMEMBER: output in meters/s
;
fall_speed=fltarr(n_elements(crystal_length))
;
fall_speed=fall_speed_simple(crystal_length(*), $
	env_pressure)			;m/s
mass_of_layer=total(concentration(*)*$
	crystal_mass(*))			;kg/m^3
mass_weighted_speed=total(concentration(*)*$
	crystal_mass(*)*fall_speed(*))	
;
weighted_speed=mass_weighted_speed(*)/mass_of_layer(*)
;
return,weighted_speed			;m/s
end

function capacitance_of_length,length
;
; returns the capacitance of the (presumably) bullet shaped ice crystal
; in meters.  Uses prolate spheroid approximation of Pruppacher & Klett, 
; p. 449 where a=length/2 and b=radius of sphere of equal surface area
; REMEMBER: this function must be compiled after length_to_equiv_rad()
; function changed to returened capacitance values recommended by 
; Heymsfield '75 p. 821
;
a=length/2.
b=length_to_equiv_rad(length)
big_A=sqrt(a*a-b*b)
capacitance=big_A/alog((a+big_a)/b)
return,capacitance
end

function area_of_length,length
;
; returns the area of the (presumably) bullet shaped ice crystal
; in kg.  Uses parameterization by Ramaswamy & Detwiler '86 p. 2290
; REMEMBER: input length in meters
;
length_cm=length*100.
surface_area_cm=0.0334*(length_cm^1.572)+.505*(length_cm^1.786)
surface_area=surface_area_cm*1.0e-4	;cm^2 --> m^2
return,surface_area			
end

function liou_area_of_length_diameter,length,diameter
;
; returns the surface area of the (presumably) hexagonal cylinder shaped 
; ice crystal in m^2.  see e.g. ebert & curry '92
; REMEMBER: input length in meters
;
basal_plane_radius=diameter/2.		;m
float_foo=3.*sqrt(3.)/2.
basal_plane_area=float_foo*basal_plane_radius*basal_plane_radius
hex_surface_area=2.*basal_plane_area+6.*basal_plane_radius*length
return,hex_surface_area
end

function liou_equiv_rad_of_length_diameter,length,diameter
;
; returns the spherical equivalent radius 
; of the (presumably) hexagonal cylinder shaped 
; ice crystal in m.  see e.g. ebert & curry '92
; REMEMBER: input length in meters
;
basal_plane_radius=diameter/2.		;m
float_foo=3.*sqrt(3.)/2.
basal_plane_area=float_foo*basal_plane_radius*basal_plane_radius
hex_surface_area=2.*basal_plane_area+6.*basal_plane_radius*length
equiv_rad=sqrt(hex_surface_area/(4.*!pi))
return,equiv_rad
end

function liou_mass_of_length_diameter,length,diameter
;
; returns the mass of the (presumably) hexagonal cylinder shaped 
; ice crystal in kg.  see e.g. ebert & curry '92
; REMEMBER: input length in meters
;
a=diameter/2.		;m
density=0.8*.001*1.0e6	;kg/m^3
volume=3.*sqrt(3.)*a*a*length/2. ;
mass=density*volume	;kg
return,mass
end

function mass_of_length,length
;
; returns the mass of the (presumably) bullet shaped ice crystal
; in kg.  Uses parameterization by Ramaswamy & Detwiler '86 p. 2290
; REMEMBER: input length in meters
;
length_cm=length*100.
mass=0.0137*(length_cm^2.572)		;g
return,mass/1000.			;kg
end


function length_of_mass,mass
;
; returns the length of the (presumably) bullet shaped ice crystal
; in m.  Uses parameterization by Ramaswamy & Detwiler '86 p. 2290
; REMEMBER: input mass in kg
;
;print,'length_of_mass reports back mass = ',mass(0)
mass_g=mass*1000.
length_cm=(mass_g/.0137)^(1./2.572)	;cm
return,length_cm/100.			;m
end

function IWC_of_temp,temperature
;
; return the Ice Water Content (IWC) of a cirrus cloud as parameterized
; by Liou '86 and used in Zhang et al '89
; for the given temperature in Kelvin
; REMEMBER: IWC is returned in kg/m^3
;
t_celsius=temperature-273.15
foo=-.2443e-3*(abs(t_celsius)-20)^2.455
IWC=-7.6+4.*exp(foo)
IWC=exp(IWC)			;g/m^3
return,IWC/1000.		;kg/m^3
end

function sz_dst_of_IWC,length,temperature
;
; return the (usually initial) size dist. of a cirrus cloud as parameterized
; by Heymsfield and Platt '85 and used in Zhang et al '89
; the length is input in meters, and the number of particles per cubic
; meter air per micron of length is returned.  This number should then
; be multiplied by a size domain (in microns) to obtain an absolute 
; concentration in #/m^3.
; for the given temperature in Kelvin
;
; REMEMBER: this routine calls IWC_of_temp, so must be compiled after it.
; REMEMBER: Neither of the distributions is fully satisfactory.
;
length_mcr=length*1.0e6	;m -> microns
T_celsius=temperature-273.15
b1=[-2.56,-2.51,-2.21,-2.29,-3.23,-3.15,-3.83,-3.85]
Nsub100=[140,175,130,250,25.5,14.0,7.00,5.02]	;#/m^3/micron
Nsub100byIWC=[5170,7000,7430,19800,7500,5600,3890,5580] ;#/g/micron
;
if T_celsius gt -20 then begin
	print,'temperature too high in function sz_dst_of_IWC'
	return,0.
endif
if T_celsius gt -25 then index = 0 else $
if T_celsius gt -30 then index = 1 else $
if T_celsius gt -35 then index = 2 else $
if T_celsius gt -40 then index = 3 else $
if T_celsius gt -45 then index = 4 else $
if T_celsius gt -50 then index = 5 else $
if T_celsius gt -55 then index = 6 else $
if T_celsius ge -60 then index = 7 else $
if T_celsius lt -60 then begin
	print,'temperature too low in function sz_dst_of_IWC'
	return,0.
endif
;
IWC=IWC_of_temp(temperature)	;kg/m^3
IWC_mgs=IWC*1000.		;g/m^3
;
; The following distribution is what Heymsfield spells out, but 
; it makes the distribution independent of IWC, which is unintuitive.
;
;sz_dst=Nsub100(index)*(length_mcr/100)^b1(index)	;#/m^3/micron
;
; The following distribution is a combination between Liou's 
; IWC parameterization and Heymsfield's distribution which retains
; some temperature dependence.  
;
a1=Nsub100byIWC(index)/(100.^b1(index))
sz_dst=a1*IWC_mgs*(length_mcr^b1(index))	;#/m^3/micron
;
; Is it a1 TIMES IWC (Heymsfield), or a1 OF IWC (Zhang)? 
;
return,sz_dst
end

function crystal_dens_of_length,length
;
; return the ice_density of a cirrus cloud crystal as parameterized
; by its length by Heymsfield '72 p. 1352-53 and used in Zhang et al '89
; for the length given in meters
; REMEMBER: ice_density is returned in kg/m^3
;
;ice_density=0.81*length^(-0.054)	;g/m^3
;return,ice_density*1.0e6/1000.	;kg/m^3
;
; Hardwiring in .8 g/cm^3 seems most reasonable, see czp 44.
;
return,0.8*1.0e6/1000.
end

function k_pruppacher,temperature=temperature,radius=radius, $
	deltat=deltat,alphat=alphat,env_density=env_density
;
; computes the correction to ksuba for small particles due
; to the thermal/kinetic effects close to the particle.
; Pruppacher and Klett p. 415-418
; unit of k_pruppacher are J/m/s/K
; REMEMBER to input temperature in Kelvin!
;
if n_elements(env_density) eq 0 then env_density = .3	;kg/m^3
if n_elements(deltat) eq 0 then deltat = 2.16e-7 	;m
if n_elements(alphat) eq 0 then alphat = .7
if n_elements(radius) eq 0 then radius = 5e-6		;m
if n_elements(temperature) eq 0 then temperature = 228  ;K
;if n_elements(ksuba) eq 0 then ksuba = .022		;J/m/s/C
;if n_elements() eq 0 then = 
;
spec_heat_air=1005		;J/kg/K
gas_const_dry_air=287.05	;J/kg/K
;
t_celsius=temperature-273.15
ksuba=(5.69+0.017*t_celsius)*1.0e-5	;cal/cm/s/C
ksuba=ksuba*4.1855*100			;cal/cm -> J/m
speed=2.*!pi/(gas_const_dry_air*temperature)
speed=sqrt(speed)
denom_1=radius/(radius+deltat)
denom_2=ksuba/(radius*alphat*env_density*spec_heat_air)
;print,'denom_2 = ',denom_2
denom_2=denom_2*speed
k_pruppacher=ksuba/(denom_1+denom_2)
return,k_pruppacher
end

function D_pruppacher,temperature=temperature,radius=radius, $
	deltat=deltat,alphat=alphat,env_density=env_density, $
	env_pressure=env_pressure
;
; computes the correction to dsuba for small particles due
; to the thermal/kinetic effects close to the particle.
; Pruppacher and Klett p. 413-418
; unit of d_pruppacher are m^2/s
; REMEMBER to input temperature in Kelvin!
;
if n_elements(env_density) eq 0 then env_density = .3	;kg/m^3
if n_elements(deltat) eq 0 then deltat = 2.16e-7 	;m
if n_elements(alphat) eq 0 then alphat = .7
if n_elements(radius) eq 0 then radius = 5e-6		;m
if n_elements(env_pressure) eq 0 then env_pressure = 25000 ;Pa
if n_elements(temperature) eq 0 then temperature = 228  ;K
;if n_elements() eq 0 then = 
;
spec_heat_air=1005		;J/kg/K
gas_const_dry_air=287.05	;J/kg/K
;
t_celsius=temperature-273.15
env_p_mb=env_pressure/100.
d_pruppacher=.211*((temperature/273.15)^1.94)*(1013.25/env_p_mb) ;cm^2/s
;
return,d_pruppacher/1.0e4	;cm^2/s --> m^2/s
end

