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

pro crm_plt,fl_in=fl_in,pause=pause,ps=ps,gangplot=gangplot
; Purpose: plot CRM text files

; Usage:
; crm_plt,fl_in=/data/zender/aca/mls_clr.txt

if n_elements(pause) eq 0 then pause = 'y';ready to roll
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
if n_elements(fl_in) eq 0 then fl_in='/data/zender/aca/mls_clr.txt'
if pause eq 'y' then begin
        !p.multi=0 
        gangplot = 'n'
endif else begin
        !p.multi=[0,3,2,0,0];=[0,num_rows,num_cols,0,0]
        pause = 'n'
        gangplot = 'y'
endelse

sys_time = systime(1) 

close,/all
openr,1,fl_in
lev_nbr=18
lev_nbr=fix(lev_nbr)
foo=""
readf,1,foo
readf,1,foo
readf,1,foo
readf,1,foo
readf,1,foo
readf,1,foo

altitude=fltarr(lev_nbr)
pressure=fltarr(lev_nbr)
temperature=fltarr(lev_nbr)
air_density=fltarr(lev_nbr)
ozone_density=fltarr(lev_nbr)
vapor_density=fltarr(lev_nbr)

for point=0,lev_nbr-1 do begin
        readf,1,a,b,c,d,e,f
	altitude(point)=a
	pressure(point)=b
	temperature(point)=c
	air_density(point)=d
	vapor_density(point)=e
	ozone_density(point)=f
endfor
close,1

epsilon_vapor=0.622
gas_const_vapor=461.51
vapor_mmr=vapor_density/air_density
ozone_mmr=ozone_density/air_density
sat_vap_pres_water=eqm_vap_water(temperature)
sat_vap_pres_ice=eqm_vap_ice(temperature)
vapor_pressure=vapor_density*gas_const_vapor*temperature
RH_water=vapor_pressure/sat_vap_pres_water
RH_ice=vapor_pressure/sat_vap_pres_ice
sat_mmr_water=epsilon_vapor*sat_vap_pres_water/(pressure+sat_vap_pres_water)
sat_mmr_ice=epsilon_vapor*sat_vap_pres_ice/(pressure+sat_vap_pres_ice)
specific_humidity=vapor_mmr/(1.0+vapor_mmr)

print,'lev_nbr = ',lev_nbr

plot,pressure/100.0,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Pressure !8T!6', $
        xtit='Pressure !8T!6 (mb)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot,temperature,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Temperature !8T!6', $
        xtit='Temperature !8T!6 (K)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot,air_density,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Air_Density !7q!6!Iair!N!6', $
        xtit='Air_Density !7q!6!Iair!N!6 (kg/m!E3!N)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
	xstyle=1, $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot,ozone_density,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Ozone_Density !7q!6!Iozone!N!6', $
        xtit='Ozone_Density !7q!6!Iozone!N!6 (kg/m!E3!N)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
	xstyle=1, $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot_oi,vapor_density,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Vapor_Density !7q!6!Ivapor!N!6', $
        xtit='Vapor Density !7q!6!Ivapor!N!6 (kg/m!E3!N)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
	xstyle=1, $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot_oi,vapor_mmr*1000.0,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Vapor Mixing Ratio !8r!6', $
        xtit='Vapor Mixing Ratio !8r!6 (g/kg)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,sat_mmr_water*1000.0,altitude/1000.0, linestyle=1, thick=2.
oplot,sat_mmr_ice*1000.0,altitude/1000.0, linestyle=2, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,vapor_mmr*1000.0,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Vapor Mixing Ratio !8r!6', $
        xtit='Vapor Mixing Ratio !8r!6 (g/kg)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,sat_mmr_water*1000.0,altitude/1000.0, linestyle=1, thick=2.
oplot,sat_mmr_ice*1000.0,altitude/1000.0, linestyle=2, thick=2.
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot,RH_water*100.0,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Relative Humidity !8RH!6', $
        xtit='Relative Humidity !8RH!6 (%)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,RH_ice*100.0,altitude/1000.+1.0, linestyle=1, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
plot,RH_water*100.0,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Relative Humidity !8RH!6', $
        xtit='Relative Humidity !8RH!6 (%)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,RH_ice*100.0,altitude/1000.+1.0, linestyle=1, thick=2.
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

plot,ozone_mmr*1.0e6,altitude/1000.0, $
	tit='!6 Altitude !8z!6 vs Ozone Mixing Ratio !8r!6', $
        xtit='Ozone Mixing Ratio !8r!6 (ppm)', $
        ytit='Altitude !8z!6 (km)', $
	yrange=[0.0,20.0], $
        linestyle=0, TICKLEN=1.0, thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif
 endif ;endif pause

;***************enforced exit********************************
 goto,exit_gracefully

 exit_gracefully: foo=1
 end


