pro mie_grapher,m=m,pause=pause,ps=ps,gangplot=gangplot
;
; Example usage 
;
;mie_grapher
;
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
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
openr,1,'mie_csz_batch.out'
readf,1,num_points
num_points=fix(num_points)
readf,1,wavelength_microns,ice_real_index,ice_imag_index
foo=""
readf,1,foo

wave=fltarr(num_points)
radius=fltarr(num_points)
nreal=fltarr(num_points)
nimag=fltarr(num_points)
x=fltarr(num_points)
qsca=fltarr(num_points)
qabs=fltarr(num_points)
qext=fltarr(num_points)
asym=fltarr(num_points)
omega=fltarr(num_points)

for point=0,num_points-1 do begin
	readf,1,a,b,c,d,e,f,g,h,i,j
	wave(point)=a
	radius(point)=b
	nreal(point)=c
	nimag(point)=d
	x(point)=e
	qsca(point)=f
	qabs(point)=g
	qext(point)=h
	asym(point)=i
	omega(point)=j
endfor
close,1

print,'num_points = ',num_points
print,wavelength_microns,ice_real_index,ice_imag_index,' = wavelength_microns, ice_real_index, ice_imag_index'

wavelength_string=string(Format='(f6.2)',wavelength_microns)
ice_real_index_string=string(Format='(f6.2)',ice_real_index)
ice_imag_index_string=string(Format='(e9.3)',ice_imag_index)

plot_oi,radius,qsca, $
	tit='!6 Radius !8r!6 vs Scattering Efficiency !8Q!6!Isca!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Isca!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0, thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,radius,qsca, $
	tit='!6 Radius !8r!6 vs Scattering Efficiency !8Q!6!Isca!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Isca!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0, thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

plot_oi,radius,qabs, $
	tit='!6 Radius !8r!6 vs Absorption Efficiency !8Q!6!Iabs!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iabs!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,radius,qabs, $
	tit='!6 Radius !8r!6 vs Absorption Efficiency !8Q!6!Iabs!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iabs!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

plot_oi,radius,qext, $
	tit='!6 Radius !8r!6 vs Extinction Efficiency !8Q!6!Iext!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iext!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,radius,qext, $
	tit='!6 Radius !8r!6 vs Extinction Efficiency !8Q!6!Iext!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iext!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

plot_oo,radius,!pi*radius*radius*qext, $
	tit='!6 Radius !8r!6 vs Extinction Cross Section !7r!6!Iext!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!7r!6!Iext!N!6 (!7l!6!E2!N)', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oo,radius,!pi*radius*radius*qext, $
	tit='!6 Radius !8r!6 vs Extinction Cross Section !7r!6!Iext!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!7r!6!Iext!N!6 (!7l!6!E2!N)', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

plot_oi,radius,asym, $
	tit='!6 Radius !8r!6 vs Asymmetry Factor !8g!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8g!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,radius,asym, $
	tit='!6 Radius !8r!6 vs Asymmetry Factor !8g!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8g!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

plot_oo,radius,1.-omega, $
	tit='!6 Radius !8r!6 vs Single Scattering Co-Albedo 1-!7x!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='1-!7x!6', $
	xrange=[0.1,100.0], $
	yrange=[1.0e-6,1.0e-0], $
	linestyle=0, TICKLEN=1.0, thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot_oo,radius,1.-omega, $
	tit='!6 Radius !8r!6 vs Single Scattering Co-Albedo 1-!7x!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='1-!7x!6', $
	xrange=[0.1,100.0], $
	yrange=[1.0e-6,1.0e-0], $
	linestyle=0, TICKLEN=1.0, thick=2.

xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.5,/NORMAL
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.5,/NORMAL
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.5,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

 if gangplot eq 'y' then begin
   print,' Hit any key to continue, p to gangplot all graphs, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	landplot
;******************** Printing routine, should be identical to above ***************
erase
!P.multi(0)=6
plot_oi,radius,qsca, $
	tit='!6 Radius !8r!6 vs Scattering Efficiency !8Q!6!Isca!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Isca!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0, thick=2.

!P.multi(0)=5
plot_oi,radius,qabs, $
	tit='!6 Radius !8r!6 vs Absorption Efficiency !8Q!6!Iabs!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iabs!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

!P.multi(0)=4
plot_oi,radius,qext, $
	tit='!6 Radius !8r!6 vs Extinction Efficiency !8Q!6!Iext!N!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8Q!6!Iext!N!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

!P.multi(0)=3
plot_oi,radius,asym, $
	tit='!6 Radius !8r!6 vs Asymmetry Factor !8g!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='!8g!6', $
	xrange=[0.1,100.0], $
	linestyle=0, TICKLEN=1.0 , thick=2.

!P.multi(0)=2
plot_oo,radius,1.-omega, $
	tit='!6 Radius !8r!6 vs Single Scattering Co-Albedo 1-!7x!6 for Ice Spheres', $
	xtit='Radius !8r!6 (!7l!6)', $
	ytit='1-!7x!6', $
	xrange=[0.1,100.0], $
	yrange=[1.0e-6,1.0e-0], $
	linestyle=0, TICKLEN=1.0, thick=2.

!P.multi(0)=1
plot,[0,0],[0,0],xrange=[0.0,1.0],yrange=[0.0,1.0],/NODATA
xyouts,0.20,0.8,'Wavelength !7k!6 = '+wavelength_string+'!7l!6',size=1.0
xyouts,0.20,0.7,'!8n!I!6real!N(!7k!6) = '+ice_real_index_string,size=1.0
xyouts,0.20,0.6,'!8n!I!6imag!N(!7k!6) = '+ice_imag_index_string,size=1.0
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

;***************enforced exit********************************
 goto,exit_gracefully

 exit_gracefully: foo=1
end
