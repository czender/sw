; $Id$ -*-f90-*-

; Purpose: Store raw data from LGGE snow reflectance experiments
;          Provide formatted arrays and mappings for other routines

; Usage:
; Insert following line where field/units strings needed
; if(.not.isvar("snw_rfl_ncl")) then load "$HOME/anl/snw_rfl.ncl" end if

; Personal headers
;if(.not.isvar("csz_ncl")) then load "${HOME}/sw/ncl/csz.ncl" end if

; Indicate to host program that header has been loaded
snw_rfl_ncl=True
print("dbg: Loading $HOME/anl/snw_rfl.ncl")

undef("rfl_1310_to_ssa")
function rfl_1310_to_ssa(rfl_1310nm:float)
; Purpose: Convert DUFISS-measured 1310 nm reflectance to 
; Specific Surface Area (SSA) using latest GDZ09 parameterization
local ssa
begin
rfl_1310nm_pct=rfl_1310nm*100.0
ssa_cm2xg=0.00010*rfl_1310nm_pct^4-0.00678*rfl_1310nm_pct^3+0.27185*rfl_1310nm_pct^2+0.49808*rfl_1310nm_pct+12.72844
ssa=ssa_cm2xg/10.0
return(ssa)
end ; end rfl_1310_to_ssa()

undef("rfl_1310_to_ssa_old")
function rfl_1310_to_ssa_old(rfl_1310nm:float)
; Purpose: Convert DUFISS-measured 1310 nm reflectance to 
; Specific Surface Area (SSA) using 2008 EGU parameterization
local ssa
begin
ssa=2.8772*exp(5.79*rfl_1310nm) ; [m2 kg-1]
return(ssa)
end ; end rfl_1310_to_ssa_old()

if (.not.isvar("mss_val_flt")) then mss_val_flt=1.0e36 end if

; Calibrated Standards

; NB: Reflectance for nominal 0% band can be derived from other five bands
; Should measure nominal 0% reflectance with container upside-down, flush with integrating sphere
; When 0% values in the standard arrays are extrapolated from other five values,
; they are not independent measurements and should not be treated as such.
; It is better to use the standard calibration for the empty sample-holders:
; 25mm-deep plastic "delrin" reflectance=6% (FD 20071220)
; 13mm-deep metal "annoblack" 1310nm reflectance=2.775% (FD 20071220)
; Tabular data and pix are in ${DATA}/lgge/lgge_rfl_std_13mm.[jpg,xls]

rfl_std_nbr=6 ; [nbr] Number of reflectance standards
rfl_nmn_0635nm=new(rfl_std_nbr,float,mss_val_flt) ; [frc] Nominal reflectance
rfl_nmn_0635nm=(/0.0,0.100,0.250,0.500,0.700,0.990/) ; [frc] Nominal reflectance
rfl_std_0635nm_13mm=(/0.03,0.10875,0.30885,0.5212,0.7012,0.987/) ; [frc] Standard reflectance
rfl_std_0635nm_25mm=(/0.06,0.10875,0.30885,0.5212,0.7012,0.987/) ; [frc] Standard reflectance

rfl_nmn_1310nm=new(rfl_std_nbr,float,mss_val_flt) ; [frc] Nominal reflectance
rfl_nmn_1310nm=rfl_nmn_0635nm ; [frc] Nominal reflectance
rfl_std_1310nm_13mm=(/0.02775,0.0943,0.2568,0.4667,0.6523,0.985/) ; [frc] Standard reflectance
rfl_std_1310nm_25mm=(/0.060,0.0943,0.2568,0.4667,0.6523,0.985/) ; [frc] Standard reflectance

; Experimental apparatus
wvl_0635nm=0.635e-6 ; [m] Wavelength
wvl_1310nm=1.310e-6 ; [m] Wavelength
mss_smp_hld_25mm=37.75 ; [g] Mass of 25mm-deep "large" sample holder
;mss_smp_hld_25mm=37.75 ; [g] Weighed by FD on 20071206
;mss_smp_hld_25mm=37.45 ; [g] Weighed by CZ on 20071207
mss_smp_hld_13mm=42.00 ; [g] Mass of 13mm-deep "small" sample holder
dpt_smp_hld_13mm=1.327 ; [cm] Depth of 13mm-deep "small" sample holder
dpt_smp_hld_25mm=2.50 ; [cm] Depth of 25mm-deep "large" sample holder
dmt_smp_hld=6.33 ; [cm] Diameter of sample holders
if (.not.isvar("pi")) then pi=4.0*atan(1.0) end if ; [frc] 3
vlm_smp_hld_13mm=0.25*pi*dmt_smp_hld*dmt_smp_hld*dpt_smp_hld_13mm ; [cm3] Volume of 13mm-deep "small" sample holder
vlm_smp_hld_25mm=0.25*pi*dmt_smp_hld*dmt_smp_hld*dpt_smp_hld_25mm ; [cm3] Volume of 25mm-deep "large" sample holder

; Daily calibration curves
; [mV] Measured signal of reflectance standards
; 0% nominal reflectance was mis-measured from 20071206--20071211
; Inadvertently held sample holder right-side up so inside was not flush with sphere
; Correct method is to turn upside-down and measure bottom reflectance 
; Replace 0% nominal values (except 20071206a) with mean 20080130--20080131 value = 0.497
sgn_std_0635nm_20071206a=(/0.412,0.955,2.563,4.655,6.675,10.875/)
sgn_std_0635nm_20071206f=(/0.497,0.996,2.653,4.785,6.675,10.875/) ; 0% from 20080130i, 65% and 99% from 20071206a
sgn_std_0635nm_20071207a=(/0.497,1.002,2.664,4.783,6.878,11.257/) ; 0% from 20080130i
sgn_std_0635nm_20071210e=(/0.497,0.991,2.655,4.774,6.887,11.208/) ; 0% from 20080130i
sgn_std_0635nm_20071211b=(/0.497,0.996,2.657,4.790,6.876,11.220/) ; 0% from 20080130i
sgn_std_0635nm_20080130c=(/0.493,0.985,2.625,4.746,6.833,11.123/)
sgn_std_0635nm_20080130i=(/0.497,0.986,2.625,4.754,6.844,11.118/)
sgn_std_0635nm_20080131b=(/0.504,0.989,2.636,4.740,6.824,11.127/)

sgn_std_1310nm_20071207c=(/5.130,7.070,13.82,24.18,34.60,58.27/) ; 0% from 20080130d
sgn_std_1310nm_20071210b=(/1.578,2.342,4.533,7.914,11.29,18.92/) ; NB: ~1/3 "normal" NIR signals
sgn_std_1310nm_20071211e=(/5.130,7.079,13.77,24.17,34.51,58.09/) ; 0% from 20080130d
sgn_std_1310nm_20080130d=(/5.130,7.097,13.77,24.07,34.47,58.02/)
sgn_std_1310nm_20080131d=(/5.136,7.110,13.79,24.06,34.44,58.08/)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Cubic spline data
; csa1: http://www.ncl.ucar.edu/Document/Functions/Built-in/csa1.shtml
; plotting: http://www.ncarg.ucar.edu//ngmath/csagrid/example1.n.html
; Knots: http://www.ncarg.ucar.edu//ngmath/csagrid/expanded.html#KnotSelection
; Interpolation: http://www.ncl.ucar.edu/Document/Functions/interp.shtml
; Fitgrid: http://www.ncarg.ucar.edu//ngmath/fitgrid/fithome.html
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rfl_btm_clc=False
if (rfl_btm_clc) then 
   ; NB: Before 200801, I extrapolated bottom reflectance
   ; Now I prescribe it as a known boundary condition
   ; Fit nominal "bottom" reflectance to other bands with two different methods
   print(rfl_btm_0635nm_20071206a)
   rfl_btm_0635nm_20071206a=csa1(sgn_std_0635nm_20071206a(1:rfl_std_nbr-1),rfl_std_0635nm_25mm(1:rfl_std_nbr-1),dimsizes(sgn_std_0635nm_20071206a)-1,sgn_std_0635nm_20071206a(0))
   print(rfl_btm_0635nm_20071206a)
   rfl_btm_0635nm_20071206a=ftcurv(sgn_std_0635nm_20071206a(1:rfl_std_nbr-1),rfl_std_0635nm_25mm(1:rfl_std_nbr-1),sgn_std_0635nm_20071206a(0))
   print(rfl_btm_0635nm_20071206a)
end if ; end if rfl_btm_clc   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20071206: Initial dirty snow experiment---nuclear winter
; Signal (sgn) [mV] and mass (mss) [g] of dirty (sal) and clean (prp) experiments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mss_prp_0635nm_20071206b=(/53.45/)-mss_smp_hld_25mm
abc_nbr_20071206b=dimsizes(mss_prp_0635nm_20071206b)
dns_prp_0635nm_20071206b=mss_prp_0635nm_20071206b/vlm_smp_hld_25mm
sgn_prp_0635nm_20071206b=(/8.594/)
rfl_prp_0635nm_20071206b=csa1(sgn_std_0635nm_20071206a,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071206a)-1,sgn_prp_0635nm_20071206b)

; 20071206: Mix 0.0938 g Monarch120 into 424 g snow
; Snow is from Chamrousse (when?), stored in silver crate
mss_mpr_20071206=0.0938 ; [g] Sample #1
mss_snw_20071206=440.0-mss_prp_0635nm_20071206b ; 424.30 [g] Total bowl minus clean snow calibration
mmr_mpr_20071206=mss_mpr_20071206/mss_snw_20071206 ; 221 ppmm [ug/g] 

mss_sal_0635nm_20071206c=(/55.70,55.95,56.70,56.55,55.45/)-mss_smp_hld_25mm
abc_nbr_20071206c=dimsizes(mss_sal_0635nm_20071206c)
dns_sal_0635nm_20071206c=mss_sal_0635nm_20071206c/vlm_smp_hld_25mm
sgn_sal_0635nm_20071206c=(/1.503,1.545,1.513,1.554,1.559/) ; NB: Signals slowly increase by up to 10% during measurement times (~2 minutes). 
rfl_sal_0635nm_20071206c=csa1(sgn_std_0635nm_20071206a,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071206a)-1,sgn_sal_0635nm_20071206c)

; 20071206d done with small sample holder
mss_sal_0635nm_20071206d=(/9.70,9.90,9.90,10.25,9.95,9.95/) ; Scale tared OK
abc_nbr_20071206d=dimsizes(mss_sal_0635nm_20071206d)
dns_sal_0635nm_20071206d=mss_sal_0635nm_20071206d/vlm_smp_hld_13mm
sgn_sal_0635nm_20071206d=(/1.455,1.449,1.528,1.495,1.548,1.438/)
rfl_sal_0635nm_20071206d=csa1(sgn_std_0635nm_20071206a,rfl_std_0635nm_13mm,dimsizes(sgn_std_0635nm_20071206a)-1,sgn_sal_0635nm_20071206d)

; 20071206e done with large sample holder
; 20080404: Noted erroneously large densities
mss_sal_0635nm_20071206e=(/18.05,17.70,18.15,18.95,17.95/) ; Scale tared OK
abc_nbr_20071206e=dimsizes(mss_sal_0635nm_20071206e)
dns_sal_0635nm_20071206e=mss_sal_0635nm_20071206e/vlm_smp_hld_25mm
sgn_sal_0635nm_20071206e=(/1.516,1.400,1.513,1.536,1.526/)
rfl_sal_0635nm_20071206e=csa1(sgn_std_0635nm_20071206a,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071206a)-1,sgn_sal_0635nm_20071206e)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20071207: Measure 20071206 snow again in visible, for aging effects, and in NIR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mss_sal_0635nm_20071207b=(/59.55/)-mss_smp_hld_25mm
abc_nbr_20071207b=dimsizes(mss_sal_0635nm_20071207b)
dns_sal_0635nm_20071207b=mss_sal_0635nm_20071207b/vlm_smp_hld_25mm
sgn_sal_0635nm_20071207b=(/1.515/)
rfl_sal_0635nm_20071207b=csa1(sgn_std_0635nm_20071207a,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071207a)-1,sgn_sal_0635nm_20071207b)

mss_sal_1310nm_20071207d=mss_sal_0635nm_20071207b ; Same sample as b
abc_nbr_20071207d=dimsizes(mss_sal_1310nm_20071207d)
dns_sal_1310nm_20071207d=dns_sal_0635nm_20071207b ; Same sample as b
sgn_sal_1310nm_20071207d=(/11.32/)
rfl_sal_1310nm_20071207d=csa1(sgn_std_1310nm_20071207c,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071207c)-1,sgn_sal_1310nm_20071207d)

mss_prp_1310nm_20071207e=(/59.95/)-mss_smp_hld_25mm ; Scale not tared
abc_nbr_20071207e=dimsizes(mss_prp_1310nm_20071207e)
dns_prp_1310nm_20071207e=mss_prp_1310nm_20071207e/vlm_smp_hld_25mm
sgn_prp_1310nm_20071207e=(/19.88/)
rfl_prp_1310nm_20071207e=csa1(sgn_std_1310nm_20071207c,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071207c)-1,sgn_prp_1310nm_20071207e)
ssa_prp_1310nm_20071207e=rfl_1310_to_ssa(rfl_prp_1310nm_20071207e) ; [m2 kg-1]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20071210: Mix 0.0070 g Monarch120 into 440 g snow
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mss_mpr_20071210=0.0070 ; [g] 
mss_snw_20071210=440.0 ; [g] 
mmr_mpr_20071210=mss_mpr_20071210/mss_snw_20071210 ; 15.9 ppmm [ug/g] 

mss_sal_0635nm_20071210f=(/55.95,56.50,56.85,57.65,58.05,56.75,57.30,57.00,57.70,58.70/)-mss_smp_hld_25mm
abc_nbr_20071210f=dimsizes(mss_sal_0635nm_20071210f)
dns_sal_0635nm_20071210f=mss_sal_0635nm_20071210f/vlm_smp_hld_25mm
sgn_sal_0635nm_20071210f=(/5.349,5.353,5.378,5.318,5.274,5.400,5.340,5.328,5.421,5.397/)
rfl_sal_0635nm_20071210f=csa1(sgn_std_0635nm_20071210e,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071210e)-1,sgn_sal_0635nm_20071210f)

mss_prp_0635nm_20071210d=(/55.40,55.00,54.45,54.35,54.50/)-mss_smp_hld_25mm
abc_nbr_20071210d=dimsizes(mss_prp_0635nm_20071210d)
dns_prp_0635nm_20071210d=mss_prp_0635nm_20071210d/vlm_smp_hld_25mm
sgn_prp_0635nm_20071210d=(/8.810,8.665,8.560,8.561,8.615/)
rfl_prp_0635nm_20071210d=csa1(sgn_std_0635nm_20071210e,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071210e)-1,sgn_prp_0635nm_20071210d)

mss_prp_1310nm_20071210a=(/54.25,54.70,54.25/)-mss_smp_hld_25mm
abc_nbr_20071210a=dimsizes(mss_prp_1310nm_20071210a)
dns_prp_1310nm_20071210a=mss_prp_1310nm_20071210a/vlm_smp_hld_25mm
sgn_prp_1310nm_20071210a=(/7.221,7.271,7.218/)
rfl_prp_1310nm_20071210a=csa1(sgn_std_1310nm_20071210b,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071210b)-1,sgn_prp_1310nm_20071210a)
ssa_prp_1310nm_20071210a=rfl_1310_to_ssa(rfl_prp_1310nm_20071210a) ; [m2 kg-1]

mss_sal_1310nm_20071210c=(/56.65,56.40,56.20/)-mss_smp_hld_25mm
abc_nbr_20071210c=dimsizes(mss_sal_1310nm_20071210c)
dns_sal_1310nm_20071210c=mss_sal_1310nm_20071210c/vlm_smp_hld_25mm
sgn_sal_1310nm_20071210c=(/6.661,6.637,6.644/)
rfl_sal_1310nm_20071210c=csa1(sgn_std_1310nm_20071210b,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071210b)-1,sgn_sal_1310nm_20071210c)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20071211: Dilute 20071210 dirty snow by factor of ~11
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Mix 400 g of same clean snow into 40 g of 20071210 dirty snow
mmr_mpr_20071211=40.0*mmr_mpr_20071210/(40.0+400.0); 1.446 ppmm [ug/g] 

mss_sal_0635nm_20071211a=(/56.55,55.55,55.25,56.15,55.00,55.15,55.35,55.30,55.15,56.30/)-mss_smp_hld_25mm
abc_nbr_20071211a=dimsizes(mss_sal_0635nm_20071211a)
dns_sal_0635nm_20071211a=mss_sal_0635nm_20071211a/vlm_smp_hld_25mm
sgn_sal_0635nm_20071211a=(/8.107,7.988,7.941,8.088,7.951,7.980,8.067,8.016,8.003,8.103/)
rfl_sal_0635nm_20071211a=csa1(sgn_std_0635nm_20071211b,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071211b)-1,sgn_sal_0635nm_20071211a)
wvl_sal_0635nm_20071211a=fspan(wvl_0635nm,wvl_0635nm,dimsizes(sgn_sal_0635nm_20071211a))

mss_prp_0635nm_20071211c=(/58.15,56.45,55.35,56.20,55.45/)-mss_smp_hld_25mm
abc_nbr_20071211c=dimsizes(mss_prp_0635nm_20071211c)
dns_prp_0635nm_20071211c=mss_prp_0635nm_20071211c/vlm_smp_hld_25mm
sgn_prp_0635nm_20071211c=(/9.097,8.905,8.889,8.951,8.681/)
rfl_prp_0635nm_20071211c=csa1(sgn_std_0635nm_20071211b,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20071211b)-1,sgn_prp_0635nm_20071211c)

mss_prp_1310nm_20071211d=(/55.25,56.00,54.65/)-mss_smp_hld_25mm
abc_nbr_20071211d=dimsizes(mss_prp_1310nm_20071211d)
dns_prp_1310nm_20071211d=mss_prp_1310nm_20071211d/vlm_smp_hld_25mm
sgn_prp_1310nm_20071211d=(/22.790,23.432,22.337/)
rfl_prp_1310nm_20071211d=csa1(sgn_std_1310nm_20071211e,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071211e)-1,sgn_prp_1310nm_20071211d)
ssa_prp_1310nm_20071211d=rfl_1310_to_ssa(rfl_prp_1310nm_20071211d) ; [m2 kg-1]

mss_sal_1310nm_20071211f=(/56.35,55.65,55.55/)-mss_smp_hld_25mm
abc_nbr_20071211f=dimsizes(mss_sal_1310nm_20071211f)
dns_sal_1310nm_20071211f=mss_sal_1310nm_20071211f/vlm_smp_hld_25mm
sgn_sal_1310nm_20071211f=(/22.470,21.985,21.861/)
rfl_sal_1310nm_20071211f=csa1(sgn_std_1310nm_20071211e,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20071211e)-1,sgn_sal_1310nm_20071211f)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20080130: Snow sampled at Chamrousse 1750m on 20080118 (or 20080111?)
; fxm: Shaved ID label off bottom of clean (but not dirty) mss_smp_hld_25mm
; Reduces weight of clean samples by ~0.3 g (fxm) relative to dirty samples
; As of 20080404 masses below have not been adjusted for label-weight disparity
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mss_mpr_20080130=0.0033 ; [g] 
mss_snw_20080130=440.0 ; [g] 
mmr_mpr_20080130af=mss_mpr_20080130/mss_snw_20080130 ; 7.5 ppmm [ug/g] 

mss_sal_0635nm_20080130a=(/55.85,55.75,55.90,55.85,55.15,55.25,56.45,55.35,56.50,55.15/)-mss_smp_hld_25mm
abc_nbr_20080130a=dimsizes(mss_sal_0635nm_20080130a)
dns_sal_0635nm_20080130a=mss_sal_0635nm_20080130a/vlm_smp_hld_25mm
sgn_sal_0635nm_20080130a=(/6.370,6.329,6.475,6.410,6.370,6.393,6.424,6.357,6.390,6.462/)
rfl_sal_0635nm_20080130a=csa1(sgn_std_0635nm_20080130c,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20080130c)-1,sgn_sal_0635nm_20080130a)
wvl_sal_0635nm_20080130a=fspan(wvl_0635nm,wvl_0635nm,dimsizes(sgn_sal_0635nm_20080130a))

mss_prp_0635nm_20080130b=(/54.15,53.50,53.50,54.15,54.60,53.25,53.65,54.85,54.45,54.80/)-mss_smp_hld_25mm
abc_nbr_20080130b=dimsizes(mss_prp_0635nm_20080130b)
dns_prp_0635nm_20080130b=mss_prp_0635nm_20080130b/vlm_smp_hld_25mm
sgn_prp_0635nm_20080130b=(/8.663,8.737,8.334,8.614,8.635,8.433,8.420,8.545,8.605,8.612/)
rfl_prp_0635nm_20080130b=csa1(sgn_std_0635nm_20080130c,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20080130c)-1,sgn_prp_0635nm_20080130b)

mss_prp_1310nm_20080130e=(/53.95,55.15,53.85/)-mss_smp_hld_25mm
abc_nbr_20080130e=dimsizes(mss_prp_1310nm_20080130e)
dns_prp_1310nm_20080130e=mss_prp_1310nm_20080130e/vlm_smp_hld_25mm
sgn_prp_1310nm_20080130e=(/21.965,22.305,22.130/)
rfl_prp_1310nm_20080130e=csa1(sgn_std_1310nm_20080130d,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20080130d)-1,sgn_prp_1310nm_20080130e)
ssa_prp_1310nm_20080130e=rfl_1310_to_ssa(rfl_prp_1310nm_20080130e) ; [m2 kg-1]

mss_sal_1310nm_20080130f=(/56.75,55.45,54.35/)-mss_smp_hld_25mm
abc_nbr_20080130f=dimsizes(mss_sal_1310nm_20080130f)
dns_sal_1310nm_20080130f=mss_sal_1310nm_20080130f/vlm_smp_hld_25mm
sgn_sal_1310nm_20080130f=(/21.003,21.020,20.815/)
rfl_sal_1310nm_20080130f=csa1(sgn_std_1310nm_20080130d,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20080130d)-1,sgn_sal_1310nm_20080130f)

; Dilute by factor of 10
; Mix 360 g of same clean snow into 40 g of 20080130af dirty snow
mmr_mpr_20080130gh=40.0*mmr_mpr_20080130af/(40.0+360.0); 750 ppbm [ng/g] 

mss_sal_1310nm_20080130g=(/55.40,56.40,55.45/)-mss_smp_hld_25mm
abc_nbr_20080130g=dimsizes(mss_sal_1310nm_20080130g)
dns_sal_1310nm_20080130g=mss_sal_1310nm_20080130g/vlm_smp_hld_25mm
sgn_sal_1310nm_20080130g=(/22.430,22.235,22.140/)
rfl_sal_1310nm_20080130g=csa1(sgn_std_1310nm_20080130d,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20080130d)-1,sgn_sal_1310nm_20080130g)

mss_sal_0635nm_20080130h=(/54.75,56.15,56.35,55.85,54.90,55.35,55.35,54.35,55.10,55.30/)-mss_smp_hld_25mm
abc_nbr_20080130h=dimsizes(mss_sal_0635nm_20080130h)
dns_sal_0635nm_20080130h=mss_sal_0635nm_20080130h/vlm_smp_hld_25mm
sgn_sal_0635nm_20080130h=(/8.214,8.429,8.357,8.296,8.166,8.377,8.305,8.310,8.212,8.397/)
rfl_sal_0635nm_20080130h=csa1(sgn_std_0635nm_20080130i,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20080130i)-1,sgn_sal_0635nm_20080130h)
wvl_sal_0635nm_20080130h=fspan(wvl_0635nm,wvl_0635nm,dimsizes(sgn_sal_0635nm_20080130h))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 20080131: Snow sampled at Chamrousse 1750m on 20080118 (or 20080111?)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Dilute by factor of 3
; Mix 200 g of clean snow into 100 g of 20080130g dirty snow
mmr_mpr_20080131=100.0*mmr_mpr_20080130gh/(100.0+200.0); 250 ppbm [ng/g] 

mss_sal_0635nm_20080131a=(/55.45,53.95,54.15,53.85,56.05,54.60,54.35,54.55,54.90,53.85/)-mss_smp_hld_25mm
abc_nbr_20080131a=dimsizes(mss_sal_0635nm_20080131a)
dns_sal_0635nm_20080131a=mss_sal_0635nm_20080131a/vlm_smp_hld_25mm
sgn_sal_0635nm_20080131a=(/8.653,8.371,8.432,8.199,8.510,8.471,8.372,8.343,8.415,8.286/)
rfl_sal_0635nm_20080131a=csa1(sgn_std_0635nm_20080131b,rfl_std_0635nm_25mm,dimsizes(sgn_std_0635nm_20080131b)-1,sgn_sal_0635nm_20080131a)
wvl_sal_0635nm_20080131a=fspan(wvl_0635nm,wvl_0635nm,dimsizes(sgn_sal_0635nm_20080131a))

mss_sal_1310nm_20080131c=(/54.25,54.45,54.00/)-mss_smp_hld_25mm
abc_nbr_20080131c=dimsizes(mss_sal_1310nm_20080131c)
dns_sal_1310nm_20080131c=mss_sal_1310nm_20080131c/vlm_smp_hld_25mm
sgn_sal_1310nm_20080131c=(/21.577,21.760,21.770/)
rfl_sal_1310nm_20080131c=csa1(sgn_std_1310nm_20080131d,rfl_std_1310nm_25mm,dimsizes(sgn_std_1310nm_20080131d)-1,sgn_sal_1310nm_20080131c)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin measurement metadata
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Dimensionality of experimental data
abc_nbr=dimsizes(sgn_sal_0635nm_20071211a)
xpt_nbr=17 ; [nbr] Number of LGGE experimental configurations (ensembles)
mpr_nbr=7 ; [nbr] Number of different impurities
wvl_nbr=2 ; [nbr] Number of different wavelengths
gph_xpt=new((/xpt_nbr/),graphic)

; Assign indices for independent, orthogonal coordinate axes
idx_wvl_0635nm=0
idx_wvl_1310nm=1

idx_mpr_ppb_000p0=mpr_nbr-1
idx_mpr_ppb_250p0=mpr_nbr-2
idx_mpr_ppb_750p0=mpr_nbr-3
idx_mpr_ppm_001p5=mpr_nbr-4
idx_mpr_ppm_007p5=mpr_nbr-5
idx_mpr_ppm_015p9=mpr_nbr-6
idx_mpr_ppm_221p0=mpr_nbr-7

; Build coordinate variables
if(isvar("wvl_obs")) then delete(wvl_obs) end if 
wvl_obs=new(wvl_nbr,float,mss_val_flt)
wvl_obs(idx_wvl_0635nm)=wvl_0635nm
wvl_obs(idx_wvl_1310nm)=wvl_1310nm

if(isvar("mmr_mpr")) then delete(mmr_mpr) end if 
mmr_mpr=new(mpr_nbr,float,mss_val_flt)
mmr_mpr(idx_mpr_ppb_000p0)=0.0
mmr_mpr(idx_mpr_ppb_250p0)=mmr_mpr_20080131
mmr_mpr(idx_mpr_ppb_750p0)=mmr_mpr_20080130gh
mmr_mpr(idx_mpr_ppm_001p5)=mmr_mpr_20071211
mmr_mpr(idx_mpr_ppm_007p5)=mmr_mpr_20080130af
mmr_mpr(idx_mpr_ppm_015p9)=mmr_mpr_20071210
mmr_mpr(idx_mpr_ppm_221p0)=mmr_mpr_20071206

; Mapping: Enumerate all measurements into unique consecutive scalars for 1-D array
; Order should not matter for this mapping
idx_prp_0635nm_20071211c=0 
idx_sal_0635nm_20071211a=1 
idx_sal_0635nm_20071210f=2 
idx_sal_0635nm_20071206e=3 
idx_prp_0635nm_20071210d=4 
idx_prp_1310nm_20071211d=5 
idx_sal_1310nm_20071211f=6 
idx_sal_1310nm_20071210c=7 
idx_sal_1310nm_20071207d=8 
idx_sal_0635nm_20080130a=9
idx_prp_0635nm_20080130b=10
idx_prp_1310nm_20080130e=11
idx_sal_1310nm_20080130f=12
idx_sal_1310nm_20080130g=13
idx_sal_0635nm_20080130h=14
idx_sal_0635nm_20080131a=15
idx_sal_1310nm_20080131c=16

; Experiment string
xpt2sng=new(xpt_nbr,string)
xpt2sng(:)="No string yet"
xpt2sng(idx_sal_0635nm_20071211a)="20071211a_0635nm_ppm_001p5"
xpt2sng(idx_prp_0635nm_20071211c)="20071211c_0635nm_ppb_000p0"
xpt2sng(idx_sal_0635nm_20071210f)="20071210f_0635nm_ppm_015p9"
xpt2sng(idx_prp_0635nm_20071210d)="20071210d_0635nm_ppb_000p0"
xpt2sng(idx_sal_0635nm_20071206e)="20071206e_0635nm_ppm_221p0"
xpt2sng(idx_prp_1310nm_20071211d)="20071211d_1310nm_ppb_000p0"
xpt2sng(idx_sal_1310nm_20071211f)="20071211f_1310nm_ppm_001p5"
xpt2sng(idx_sal_1310nm_20071210c)="20071210c_1310nm_ppm_015p9"
xpt2sng(idx_sal_1310nm_20071207d)="20071207d_1310nm_ppm_221p0"
xpt2sng(idx_sal_0635nm_20080130a)="20080130a_0635nm_ppm_007p5"
xpt2sng(idx_prp_0635nm_20080130b)="20080130b_0635nm_ppb_000p0"
xpt2sng(idx_prp_1310nm_20080130e)="20080130e_1310nm_ppb_000p0"
xpt2sng(idx_sal_1310nm_20080130f)="20080130f_1310nm_ppm_007p5"
xpt2sng(idx_sal_1310nm_20080130g)="20080130g_1310nm_ppb_750p0"
xpt2sng(idx_sal_0635nm_20080130h)="20080130h_0635nm_ppb_750p0"
xpt2sng(idx_sal_0635nm_20080131a)="20080131a_0635nm_ppb_250p0"
xpt2sng(idx_sal_1310nm_20080131c)="20080131c_1310nm_ppb_250p0"

; Mapping: Experiment -> impurity concentration 
xpt2mpr=new(xpt_nbr,integer,"No_FillValue")
xpt2mpr(:)=0
xpt2mpr(idx_sal_0635nm_20071211a)=idx_mpr_ppm_001p5
xpt2mpr(idx_prp_0635nm_20071211c)=idx_mpr_ppb_000p0
xpt2mpr(idx_sal_0635nm_20071210f)=idx_mpr_ppm_015p9
xpt2mpr(idx_prp_0635nm_20071210d)=idx_mpr_ppb_000p0
xpt2mpr(idx_sal_0635nm_20071206e)=idx_mpr_ppm_221p0
xpt2mpr(idx_prp_1310nm_20071211d)=idx_mpr_ppb_000p0
xpt2mpr(idx_sal_1310nm_20071211f)=idx_mpr_ppm_001p5
xpt2mpr(idx_sal_1310nm_20071210c)=idx_mpr_ppm_015p9
xpt2mpr(idx_sal_1310nm_20071207d)=idx_mpr_ppm_221p0
xpt2mpr(idx_sal_0635nm_20080130a)=idx_mpr_ppm_007p5
xpt2mpr(idx_prp_0635nm_20080130b)=idx_mpr_ppb_000p0
xpt2mpr(idx_prp_1310nm_20080130e)=idx_mpr_ppb_000p0
xpt2mpr(idx_sal_1310nm_20080130f)=idx_mpr_ppm_007p5
xpt2mpr(idx_sal_1310nm_20080130g)=idx_mpr_ppb_750p0
xpt2mpr(idx_sal_0635nm_20080130h)=idx_mpr_ppb_750p0
xpt2mpr(idx_sal_0635nm_20080131a)=idx_mpr_ppb_250p0
xpt2mpr(idx_sal_1310nm_20080131c)=idx_mpr_ppb_250p0

; Mapping: Experiment -> wavelength
xpt2wvl=new(xpt_nbr,integer,"No_FillValue")
xpt2wvl(:)=0
xpt2wvl(idx_sal_0635nm_20071211a)=idx_wvl_0635nm
xpt2wvl(idx_prp_0635nm_20071211c)=idx_wvl_0635nm
xpt2wvl(idx_sal_0635nm_20071210f)=idx_wvl_0635nm
xpt2wvl(idx_prp_0635nm_20071210d)=idx_wvl_0635nm
xpt2wvl(idx_sal_0635nm_20071206e)=idx_wvl_0635nm
xpt2wvl(idx_prp_1310nm_20071211d)=idx_wvl_1310nm
xpt2wvl(idx_sal_1310nm_20071211f)=idx_wvl_1310nm
xpt2wvl(idx_sal_1310nm_20071210c)=idx_wvl_1310nm
xpt2wvl(idx_sal_1310nm_20071207d)=idx_wvl_1310nm
xpt2wvl(idx_sal_0635nm_20080130a)=idx_wvl_0635nm
xpt2wvl(idx_prp_0635nm_20080130b)=idx_wvl_0635nm
xpt2wvl(idx_prp_1310nm_20080130e)=idx_wvl_1310nm
xpt2wvl(idx_sal_1310nm_20080130f)=idx_wvl_1310nm
xpt2wvl(idx_sal_1310nm_20080130g)=idx_wvl_1310nm
xpt2wvl(idx_sal_0635nm_20080130h)=idx_wvl_0635nm
xpt2wvl(idx_sal_0635nm_20080131a)=idx_wvl_0635nm
xpt2wvl(idx_sal_1310nm_20080131c)=idx_wvl_1310nm

; Mapping: Experiment -> Parent snow sample (pss) (itself related to date)
xpt2pss=new(xpt_nbr,integer,"No_FillValue")
xpt2pss(:)=0
; map parent snow sample to raw index
idx_pss_2007120X=idx_zero ; Parent snow: 2007120x Chamrousse
idx_pss_20080118=idx_one ; Parent snow: 20080118 Chamrousse
; map date to parent snow sample
idx_dat_20071206=idx_pss_2007120X ; Parent snow: 2007120x Chamrousse
idx_dat_20071207=idx_pss_2007120X ; Parent snow: 2007120x Chamrousse
idx_dat_20071210=idx_pss_2007120X ; Parent snow: 2007120x Chamrousse
idx_dat_20071211=idx_pss_2007120X ; Parent snow: 2007120x Chamrousse
idx_dat_20080130=idx_pss_20080118 ; Parent snow: 20080118 Chamrousse
idx_dat_20080131=idx_pss_20080118 ; Parent snow: 20080118 Chamrousse
; map measurement to date
xpt2pss(idx_sal_0635nm_20071211a)=idx_dat_20071211
xpt2pss(idx_prp_0635nm_20071211c)=idx_dat_20071211
xpt2pss(idx_sal_0635nm_20071210f)=idx_dat_20071210
xpt2pss(idx_prp_0635nm_20071210d)=idx_dat_20071210
xpt2pss(idx_sal_0635nm_20071206e)=idx_dat_20071206
xpt2pss(idx_prp_1310nm_20071211d)=idx_dat_20071211
xpt2pss(idx_sal_1310nm_20071211f)=idx_dat_20071211
xpt2pss(idx_sal_1310nm_20071210c)=idx_dat_20071210
xpt2pss(idx_sal_1310nm_20071207d)=idx_dat_20071207
xpt2pss(idx_sal_0635nm_20080130a)=idx_dat_20080130
xpt2pss(idx_prp_0635nm_20080130b)=idx_dat_20080130
xpt2pss(idx_prp_1310nm_20080130e)=idx_dat_20080130
xpt2pss(idx_sal_1310nm_20080130f)=idx_dat_20080130
xpt2pss(idx_sal_1310nm_20080130g)=idx_dat_20080130
xpt2pss(idx_sal_0635nm_20080130h)=idx_dat_20080130
xpt2pss(idx_sal_0635nm_20080131a)=idx_dat_20080131
xpt2pss(idx_sal_1310nm_20080131c)=idx_dat_20080131

; Format raw data into 2-D arrays for easier plotting
if(isvar("wvl_2d")) then delete(wvl_2d) end if 
wvl_2d=new((/xpt_nbr,abc_nbr/),float,mss_val_flt)
wvl_2d(idx_sal_0635nm_20071211a,0:abc_nbr_20071211a-1)=wvl_0635nm
wvl_2d(idx_prp_0635nm_20071211c,0:abc_nbr_20071211c-1)=wvl_0635nm
wvl_2d(idx_sal_0635nm_20071210f,0:abc_nbr_20071210f-1)=wvl_0635nm
wvl_2d(idx_prp_0635nm_20071210d,0:abc_nbr_20071210d-1)=wvl_0635nm
wvl_2d(idx_sal_0635nm_20071206e,0:abc_nbr_20071206e-1)=wvl_0635nm
wvl_2d(idx_prp_1310nm_20071211d,0:abc_nbr_20071211d-1)=wvl_1310nm
wvl_2d(idx_sal_1310nm_20071211f,0:abc_nbr_20071211f-1)=wvl_1310nm
wvl_2d(idx_sal_1310nm_20071210c,0:abc_nbr_20071210c-1)=wvl_1310nm
wvl_2d(idx_sal_1310nm_20071207d,0:abc_nbr_20071207d-1)=wvl_1310nm
wvl_2d(idx_sal_0635nm_20080130a,0:abc_nbr_20080130a-1)=wvl_0635nm
wvl_2d(idx_prp_0635nm_20080130b,0:abc_nbr_20080130b-1)=wvl_0635nm
wvl_2d(idx_prp_1310nm_20080130e,0:abc_nbr_20080130e-1)=wvl_1310nm
wvl_2d(idx_sal_1310nm_20080130f,0:abc_nbr_20080130f-1)=wvl_1310nm
wvl_2d(idx_sal_1310nm_20080130g,0:abc_nbr_20080130g-1)=wvl_1310nm
wvl_2d(idx_sal_0635nm_20080130h,0:abc_nbr_20080130h-1)=wvl_0635nm
wvl_2d(idx_sal_0635nm_20080131a,0:abc_nbr_20080131a-1)=wvl_0635nm
wvl_2d(idx_sal_1310nm_20080131c,0:abc_nbr_20080131c-1)=wvl_1310nm
; wvl_2d(idx_typ__wvlnm_2008foo,0:abc_nbr_2008foo-1)=wvl__wvlnm

; Three dimensional cube is wavelength by impurity by ensemble size
if(isvar("data_3d")) then delete(data_3d) end if
data_3d=new((/wvl_nbr,mpr_nbr,abc_nbr/),float,mss_val_flt)
data_3d(idx_wvl_0635nm,idx_mpr_ppb_000p0,0:abc_nbr_20080130b-1)=rfl_prp_0635nm_20080130b ; 20071210d OK too
data_3d(idx_wvl_0635nm,idx_mpr_ppb_250p0,0:abc_nbr_20080131a-1)=rfl_sal_0635nm_20080131a
data_3d(idx_wvl_0635nm,idx_mpr_ppb_750p0,0:abc_nbr_20080130h-1)=rfl_sal_0635nm_20080130h
data_3d(idx_wvl_0635nm,idx_mpr_ppm_001p5,0:abc_nbr_20071211a-1)=rfl_sal_0635nm_20071211a
data_3d(idx_wvl_0635nm,idx_mpr_ppm_007p5,0:abc_nbr_20080130a-1)=rfl_sal_0635nm_20080130a
data_3d(idx_wvl_0635nm,idx_mpr_ppm_015p9,0:abc_nbr_20071210f-1)=rfl_sal_0635nm_20071210f
data_3d(idx_wvl_0635nm,idx_mpr_ppm_221p0,0:abc_nbr_20071206e-1)=rfl_sal_0635nm_20071206e
data_3d(idx_wvl_1310nm,idx_mpr_ppb_000p0,0:abc_nbr_20080130e-1)=rfl_prp_1310nm_20080130e ; 20071211d OK too
data_3d(idx_wvl_1310nm,idx_mpr_ppb_250p0,0:abc_nbr_20080131c-1)=rfl_sal_1310nm_20080131c
data_3d(idx_wvl_1310nm,idx_mpr_ppb_750p0,0:abc_nbr_20080130g-1)=rfl_sal_1310nm_20080130g
data_3d(idx_wvl_1310nm,idx_mpr_ppm_001p5,0:abc_nbr_20071211f-1)=rfl_sal_1310nm_20071211f
data_3d(idx_wvl_1310nm,idx_mpr_ppm_007p5,0:abc_nbr_20080130f-1)=rfl_sal_1310nm_20080130f
data_3d(idx_wvl_1310nm,idx_mpr_ppm_015p9,0:abc_nbr_20071210c-1)=rfl_sal_1310nm_20071210c
data_3d(idx_wvl_1310nm,idx_mpr_ppm_221p0,0:abc_nbr_20071207d-1)=rfl_sal_1310nm_20071207d

if(isvar("data_xpt")) then delete(data_xpt) end if
data_xpt=new((/xpt_nbr,abc_nbr/),float,mss_val_flt)
data_xpt(idx_sal_0635nm_20071211a,0:abc_nbr_20071211a-1)=rfl_sal_0635nm_20071211a
data_xpt(idx_prp_0635nm_20071211c,0:abc_nbr_20071211c-1)=rfl_prp_0635nm_20071211c
data_xpt(idx_sal_0635nm_20071210f,0:abc_nbr_20071210f-1)=rfl_sal_0635nm_20071210f
data_xpt(idx_prp_0635nm_20071210d,0:abc_nbr_20071210d-1)=rfl_prp_0635nm_20071210d
data_xpt(idx_sal_0635nm_20071206e,0:abc_nbr_20071206e-1)=rfl_sal_0635nm_20071206e
data_xpt(idx_prp_1310nm_20071211d,0:abc_nbr_20071211d-1)=rfl_prp_1310nm_20071211d
data_xpt(idx_sal_1310nm_20071211f,0:abc_nbr_20071211f-1)=rfl_sal_1310nm_20071211f
data_xpt(idx_sal_1310nm_20071210c,0:abc_nbr_20071210c-1)=rfl_sal_1310nm_20071210c
data_xpt(idx_sal_1310nm_20071207d,0:abc_nbr_20071207d-1)=rfl_sal_1310nm_20071207d
data_xpt(idx_sal_0635nm_20080130a,0:abc_nbr_20080130a-1)=rfl_sal_0635nm_20080130a
data_xpt(idx_prp_0635nm_20080130b,0:abc_nbr_20080130b-1)=rfl_prp_0635nm_20080130b
data_xpt(idx_prp_1310nm_20080130e,0:abc_nbr_20080130e-1)=rfl_prp_1310nm_20080130e
data_xpt(idx_sal_1310nm_20080130f,0:abc_nbr_20080130f-1)=rfl_sal_1310nm_20080130f
data_xpt(idx_sal_1310nm_20080130g,0:abc_nbr_20080130g-1)=rfl_sal_1310nm_20080130g
data_xpt(idx_sal_0635nm_20080130h,0:abc_nbr_20080130h-1)=rfl_sal_0635nm_20080130h
data_xpt(idx_sal_0635nm_20080131a,0:abc_nbr_20080131a-1)=rfl_sal_0635nm_20080131a
data_xpt(idx_sal_1310nm_20080131c,0:abc_nbr_20080131c-1)=rfl_sal_1310nm_20080131c
;   data_xpt(idx_typ__wvlnm_2008foo,0:abc_nbr_2008foo-1)=rfl_typ__wvlnm_2008foo

if(isvar("dns_snw")) then delete(dns_snw) end if 
dns_snw=new((/xpt_nbr,abc_nbr/),float,mss_val_flt)
dns_snw(idx_sal_0635nm_20071211a,0:abc_nbr_20071211a-1)=dns_sal_0635nm_20071211a
dns_snw(idx_prp_0635nm_20071211c,0:abc_nbr_20071211c-1)=dns_prp_0635nm_20071211c
dns_snw(idx_sal_0635nm_20071210f,0:abc_nbr_20071210f-1)=dns_sal_0635nm_20071210f
dns_snw(idx_prp_0635nm_20071210d,0:abc_nbr_20071210d-1)=dns_prp_0635nm_20071210d
dns_snw(idx_sal_0635nm_20071206e,0:abc_nbr_20071206e-1)=dns_sal_0635nm_20071206e
dns_snw(idx_prp_1310nm_20071211d,0:abc_nbr_20071211d-1)=dns_prp_1310nm_20071211d
dns_snw(idx_sal_1310nm_20071211f,0:abc_nbr_20071211f-1)=dns_sal_1310nm_20071211f
dns_snw(idx_sal_1310nm_20071210c,0:abc_nbr_20071210c-1)=dns_sal_1310nm_20071210c
dns_snw(idx_sal_1310nm_20071207d,0:abc_nbr_20071207d-1)=dns_sal_1310nm_20071207d
dns_snw(idx_sal_0635nm_20080130a,0:abc_nbr_20080130a-1)=dns_sal_0635nm_20080130a
dns_snw(idx_prp_0635nm_20080130b,0:abc_nbr_20080130b-1)=dns_prp_0635nm_20080130b
dns_snw(idx_prp_1310nm_20080130e,0:abc_nbr_20080130e-1)=dns_prp_1310nm_20080130e
dns_snw(idx_sal_1310nm_20080130f,0:abc_nbr_20080130f-1)=dns_sal_1310nm_20080130f
dns_snw(idx_sal_1310nm_20080130g,0:abc_nbr_20080130g-1)=dns_sal_1310nm_20080130g
dns_snw(idx_sal_0635nm_20080130h,0:abc_nbr_20080130h-1)=dns_sal_0635nm_20080130h
dns_snw(idx_sal_0635nm_20080131a,0:abc_nbr_20080131a-1)=dns_sal_0635nm_20080131a
dns_snw(idx_sal_1310nm_20080131c,0:abc_nbr_20080131c-1)=dns_sal_1310nm_20080131c
;   dns_snw(idx_typ__wvlnm_2008foo,0:abc_nbr_2008foo-1)=dns_typ__wvlnm_2008foo
dns_snw_avg=dim_avg(dns_snw)
if(False) then
   do xpt_idx=0,xpt_nbr-1
      print("Experiment "+xpt_idx+" = "+xpt2sng(xpt_idx)+", Density = "+dns_snw_avg(xpt_idx))
   end do; end loop over fl
   dns_foo=(dns_snw_avg(0)+dns_snw_avg(4)+dns_snw_avg(5)+dns_snw_avg(10)+dns_snw_avg(11))/5
   print(dns_foo)
   dns_foo=(dns_snw_avg(15)+dns_snw_avg(16))/2
   print(dns_foo)
   dns_foo=(dns_snw_avg(13)+dns_snw_avg(14))/2
   print(dns_foo)
   dns_foo=(dns_snw_avg(1)+dns_snw_avg(6))/2
   print(dns_foo)
   dns_foo=(dns_snw_avg(9)+dns_snw_avg(12))/2
   print(dns_foo)
   dns_foo=(dns_snw_avg(2)+dns_snw_avg(7))/2
   print(dns_foo)
   dns_foo=(dns_snw_avg(3)+dns_snw_avg(8))/2
   print(dns_foo)
end if ; !False

mmr_mpr_2d=new((/xpt_nbr,abc_nbr/),float,mss_val_flt)
mmr_mpr_2d(idx_sal_0635nm_20071211a,0:abc_nbr_20071211a-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20071211a))
mmr_mpr_2d(idx_prp_0635nm_20071211c,0:abc_nbr_20071211c-1)=mmr_mpr(xpt2mpr(idx_prp_0635nm_20071211c))
mmr_mpr_2d(idx_sal_0635nm_20071210f,0:abc_nbr_20071210f-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20071210f))
mmr_mpr_2d(idx_prp_0635nm_20071210d,0:abc_nbr_20071210d-1)=mmr_mpr(xpt2mpr(idx_prp_0635nm_20071210d))
mmr_mpr_2d(idx_sal_0635nm_20071206e,0:abc_nbr_20071206e-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20071206e))
mmr_mpr_2d(idx_prp_1310nm_20071211d,0:abc_nbr_20071211d-1)=mmr_mpr(xpt2mpr(idx_prp_1310nm_20071211d))
mmr_mpr_2d(idx_sal_1310nm_20071211f,0:abc_nbr_20071211f-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20071211f))
mmr_mpr_2d(idx_sal_1310nm_20071210c,0:abc_nbr_20071210c-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20071210c))
mmr_mpr_2d(idx_sal_1310nm_20071207d,0:abc_nbr_20071207d-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20071207d))
mmr_mpr_2d(idx_sal_0635nm_20080130a,0:abc_nbr_20080130a-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20080130a))
mmr_mpr_2d(idx_prp_0635nm_20080130b,0:abc_nbr_20080130b-1)=mmr_mpr(xpt2mpr(idx_prp_0635nm_20080130b))
mmr_mpr_2d(idx_prp_1310nm_20080130e,0:abc_nbr_20080130e-1)=mmr_mpr(xpt2mpr(idx_prp_1310nm_20080130e))
mmr_mpr_2d(idx_sal_1310nm_20080130f,0:abc_nbr_20080130f-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20080130f))
mmr_mpr_2d(idx_sal_1310nm_20080130g,0:abc_nbr_20080130g-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20080130g))
mmr_mpr_2d(idx_sal_0635nm_20080130h,0:abc_nbr_20080130h-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20080130h))
mmr_mpr_2d(idx_sal_0635nm_20080131a,0:abc_nbr_20080131a-1)=mmr_mpr(xpt2mpr(idx_sal_0635nm_20080131a))
mmr_mpr_2d(idx_sal_1310nm_20080131c,0:abc_nbr_20080131c-1)=mmr_mpr(xpt2mpr(idx_sal_1310nm_20080131c))
;   mmr_mpr_2d(idx_typ__wvlnm_2008foo,0:abc_nbr_2008foo-1)=
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End measurement metadata
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Begin model metadata
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Mapping: Enumerate models into unique scalars
mdl_nbr=mpr_nbr ; [nbr] Number of model simulations
fl_nbr=mdl_nbr ; [nbr] Number of files
if(.not.isvar("ord_nbr_plt")) then ord_nbr_plt=mdl_nbr end if ; [nbr] Ordinates plotted: ord_nbr_plt <= ord_nbr

if(isvar("fl_in")) then delete(fl_in) end if 
fl_in=new(mdl_nbr,string)

; Best estimate of everything (bee): SSA, density, phb_he
if(.not.isvar("flg_bee")) then flg_bee=False end if
if(flg_bee) then 
   fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_bee_ppb_000.0.nc"
   fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_bee_ppb_250.0.nc"
   fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_bee_ppb_750.0.nc"
   fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_bee_ppm_001.5.nc"
   fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_bee_ppm_007.5.nc"
   fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_bee_ppm_015.9.nc"
   fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_bee_ppm_221.0.nc"
end if; !flg_bee

; New best estimate of everything (tst): SSA, density, phb_he
if(.not.isvar("flg_tst")) then flg_tst=False end if
if(flg_tst) then 
   fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_tst_ppb_000.0.nc"
   fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_tst_ppb_250.0.nc"
   fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_tst_ppb_750.0.nc"
   fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_tst_ppm_001.5.nc"
   fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_tst_ppm_007.5.nc"
   fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_tst_ppm_015.9.nc"
   fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_tst_ppm_221.0.nc"
end if; !flg_tst

; Nominal SSA
if(.not.isvar("flg_nmn")) then flg_nmn=False end if
if(flg_nmn) then 
fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_100um_ppb_000.0.nc"
fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_100um_ppb_250.0.nc"
fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_100um_ppb_750.0.nc"
fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_100um_ppm_001.5.nc"
fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_100um_ppm_007.5.nc"
fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_100um_ppm_015.9.nc"
fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_100um_ppm_221.0.nc"
end if; !flg_nmn

; Measured SSA
if(.not.isvar("flg_ssa_obs")) then flg_ssa_obs=False end if
if(flg_ssa_obs) then 
   fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_ssa_obs_ppb_000.0.nc"
   fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_ssa_obs_ppb_250.0.nc"
   fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_ssa_obs_ppb_750.0.nc"
   fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_ssa_obs_ppm_001.5.nc"
   fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_ssa_obs_ppm_007.5.nc"
   fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_ssa_obs_ppm_015.9.nc"
   fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_ssa_obs_ppm_221.0.nc"
end if; !flg_ssa_obs

; Impurity optical properties: phb_he
if(.not.isvar("flg_phb_he")) then flg_phb_he=False end if
if(flg_phb_he) then 
   fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_100um_phb_he_ppb_000.0.nc"
   fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_100um_phb_he_ppb_250.0.nc"
   fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_100um_phb_he_ppb_750.0.nc"
   fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_100um_phb_he_ppm_001.5.nc"
   fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_100um_phb_he_ppm_007.5.nc"
   fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_100um_phb_he_ppm_015.9.nc"
   fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_100um_phb_he_ppm_221.0.nc"
end if; !flg_phb_he

; Impurity optical properties: phb_low
if(.not.isvar("flg_phb_low")) then flg_phb_low=False end if
if(flg_phb_low) then 
   fl_in(idx_mpr_ppb_000p0)="../../data/icr/swnb_lgge_100um_phb_low_ppb_000.0.nc"
   fl_in(idx_mpr_ppb_250p0)="../../data/icr/swnb_lgge_100um_phb_low_ppb_250.0.nc"
   fl_in(idx_mpr_ppb_750p0)="../../data/icr/swnb_lgge_100um_phb_low_ppb_750.0.nc"
   fl_in(idx_mpr_ppm_001p5)="../../data/icr/swnb_lgge_100um_phb_low_ppm_001.5.nc"
   fl_in(idx_mpr_ppm_007p5)="../../data/icr/swnb_lgge_100um_phb_low_ppm_007.5.nc"
   fl_in(idx_mpr_ppm_015p9)="../../data/icr/swnb_lgge_100um_phb_low_ppm_015.9.nc"
   fl_in(idx_mpr_ppm_221p0)="../../data/icr/swnb_lgge_100um_phb_low_ppm_221.0.nc"
end if; !flg_phb_low

if(isvar("wvl_sng")) then delete(wvl_sng) end if 
wvl_sng=new(wvl_nbr,string)
wvl_sng(idx_wvl_0635nm)="  ~F33~l~F21~ =  635 nm"
wvl_sng(idx_wvl_1310nm)="  ~F33~l~F21~ = 1310 nm"

pss_nbr=2
if(isvar("pss_sng")) then delete(pss_sng) end if 
pss_sng=new(pss_nbr,string)
pss_sng(idx_pss_2007120X)="  Snow = 2007120X"
pss_sng(idx_pss_20080118)="  Snow = 20080118"

if(isvar("dat_sng")) then delete(dat_sng) end if 
dat_sng=new(mdl_nbr,string)
dat_sng(idx_mpr_ppb_000p0)="  BC =   0 ppbm"
dat_sng(idx_mpr_ppb_250p0)="  BC = 250 ppbm"
dat_sng(idx_mpr_ppb_750p0)="  BC = 750 ppbm"
dat_sng(idx_mpr_ppm_001p5)="  BC = 1.45 ppmm"
dat_sng(idx_mpr_ppm_007p5)="  BC = 7.5 ppmm"
dat_sng(idx_mpr_ppm_015p9)="  BC = 15.9 ppmm"
dat_sng(idx_mpr_ppm_221p0)="  BC = 221 ppmm"

if(isvar("dat_sng_sht")) then delete(dat_sng_sht) end if 
dat_sng_sht=new(mdl_nbr,string)
dat_sng_sht(idx_mpr_ppb_000p0)="  0   ppbm BC"
dat_sng_sht(idx_mpr_ppb_250p0)="  250 ppbm BC"
dat_sng_sht(idx_mpr_ppb_750p0)="  750 ppbm BC"
dat_sng_sht(idx_mpr_ppm_001p5)="  1.45 ppmm BC"
dat_sng_sht(idx_mpr_ppm_007p5)="  7.5 ppmm BC"
dat_sng_sht(idx_mpr_ppm_015p9)="  15.9 ppmm BC"
dat_sng_sht(idx_mpr_ppm_221p0)="  221 ppmm BC"

vsb2xpt=ind(wvl_2d(:,0).eq.wvl_0635nm)
nir2xpt=ind(wvl_2d(:,0).eq.wvl_1310nm)
;print("vsb2xpt="+vsb2xpt)
xpt_vsb_nbr=dimsizes(vsb2xpt) ; [nbr] Number of 635nm experiments
xpt_nir_nbr=dimsizes(nir2xpt) ; [nbr] Number of 1310nm experiments
data_xpt_vsb=new((/xpt_vsb_nbr,abc_nbr/),float,mss_val_flt)
do xpt_vsb_idx=0,xpt_vsb_nbr-1
   data_xpt_vsb(xpt_vsb_idx,:)=data_xpt(vsb2xpt(xpt_vsb_idx),:)
end do
data_xpt_nir=new((/xpt_nir_nbr,abc_nbr/),float,mss_val_flt)
data_xpt_nir_avg=new((/xpt_nir_nbr/),float,mss_val_flt)
do xpt_nir_idx=0,xpt_nir_nbr-1
   data_xpt_nir(xpt_nir_idx,:)=data_xpt(nir2xpt(xpt_nir_idx),:)
end do

; Indicate to host program that input datasets have already been defined
campaign_defined=True
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End model metadata
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin model post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sensitivity to horizontal boundary conditions (mie_mca.sh)
; Multiplying by this converts plane-parallel model results to lab measurements
; Dividing    by this converts lab measurements to plane-parallel model results
if(isvar("alb_spc_snw_gmt_crc")) then delete(alb_spc_snw_gmt_crc) end if 
alb_spc_snw_gmt_crc=new(wvl_nbr,float,mss_val_flt)
alb_spc_snw_gmt_crc(idx_wvl_0635nm)=0.91 ; NB: Used for 2008 EGU presentation
alb_spc_snw_gmt_crc(idx_wvl_1310nm)=0.77 ; NB: Used for 2008 EGU presentation

if(isvar("alb_spc_snw_mop_rng")) then delete(alb_spc_snw_mop_rng) end if 
alb_spc_snw_mop_rng=new(mpr_nbr,float,mss_val_flt)
alb_spc_snw_mop_rng(idx_mpr_ppb_000p0)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppb_250p0)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppb_750p0)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppm_001p5)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppm_007p5)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppm_015p9)=0.0
alb_spc_snw_mop_rng(idx_mpr_ppm_221p0)=0.0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End model post-processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
