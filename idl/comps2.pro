function alt_basis,x,m
n=n_elements(x)
r=fltarr(n,m)
for i=0,n-1 do begin
	r(i,0)=1.
;	r(i,1)=sin(x(i)/30.)
;	r(i,1)=x(i)
;	r(i,2)=x(i)^2
;	r(i,1)=sqrt(x(i))
;	r(i,1)=atan(x(i)/(2.*!pi))
;	r(i,1)=atan(x(i)^2/20.)
	r(i,0)=1.
;	r(i,1)=atan(x(i)/(10.))
	r(i,1)=exp(-x(i)/(10.))
endfor
;on_error,2                  ;Return to caller if an error occurs
;return,( x # findgen(m))^.5 ;Couldn't be much simpler
return,r
end

function alt_basis_ords,x,coeffs
n=n_elements(x)
m=n_elements(coeffs)
y=fltarr(n)
for i=0,n-1 do begin
;	y(i)=coeffs(0)*1.+coeffs(1)*sin(x(i)/30.)
;	y(i)=coeffs(0)*1.+coeffs(1)*sqrt(x(i))
;	y(i)=coeffs(0)*1.+coeffs(1)*atan(x(i)/(2.*!pi))
;	y(i)=coeffs(0)*1.+coeffs(1)*atan(x(i)/(10.))
	y(i)=coeffs(0)*1.+coeffs(1)*exp(-x(i)/(10.))
endfor
;y=poly(abc,coeffs)
return,y
end

pro fitit,m=m,pause=pause,ps=ps
;
; Example usage 
;
;fitit,num_points=5
;
if n_elements(m) eq 0 then m = 3	;n -1 is the degree of the polynomial
					;number of coeffs in fitting func.
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
if pause eq 'y' then !p.multi=0 else begin
	!p.multi=[0,3,1,0,0];=[0,num_rows,num_cols,0,0]
	pause = 'n'
endelse

sys_time = systime(1) 

close,/all
openr,1,'comps2.dat'
readf,1,num_points
foo=""
readf,1,foo
readf,1,foo
w=fltarr(num_points)
reff_rad=fltarr(num_points)
reff_norad=fltarr(num_points)
iwp_rad=fltarr(num_points)
iwp_norad=fltarr(num_points)
yf=fltarr(num_points)
weight=fltarr(num_points)
kappa=fltarr(num_points)
sigma=fltarr(num_points)
emissivity=fltarr(num_points)

for point=0,num_points-1 do begin
	readf,1,a,b,c,d,e,f,g,h,i
	w(point)=a
	reff_rad(point)=b
	reff_norad(point)=c
	iwp_rad(point)=d
	iwp_norad(point)=e
	weight(point)=f
	kappa(point)=g
	sigma(point)=h
	emissivity(point)=i
endfor
close,1

coeffs=svdfit(w,reff_rad,m,weight=weight,funct='alt_basis',YFIT=yf)
;coeffs=svdfit(w,reff_rad,m,weight=weight,YFIT=yf)
abc=findgen(100)*max(w)/99.
ord=alt_basis_ords(abc,coeffs)
print,'num_points = ',num_points
print,'number basis functions = ',m
print,coeffs
;print,'abc = ',abc
;print,'ord = ',ord

plot,w,reff_rad, $
	tit='!6 Updraft speed vs. effective radius', $
	xtit='Maximum Updraft speed !8w!6 (cm/s)', $
	ytit='Effective Radius !8r!6!Ie!N (!7l!6)', $
	linestyle=1, thick=2.0,psym=6,/ynoz
oplot,w,reff_norad, $
	linestyle=0,thick=2.0,psym=1
oplot,abc,ord, $
	linestyle=0,thick=2.
;oplot,w,yf,linestyle=0,thick=2.0,psym=2
	
;dt_sng=string(Format='(I4)',dz)
xyouts,.50,0.6,'!9B!6   !6With Radiation',size=1.0,/NORMAL
xyouts,.50,.5,'+!6   !6Without Radiation',size=1.0,/NORMAL
xyouts,.50,.4,'!6     !6Least squares fit ',size=1.0,/NORMAL
oplot,[40,48],[34,34],linestyle=0,thick=2.
;oplot,[.5,.5],[.5,.5],linestyle=0,thick=2.0,/NORMAL

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot,w,reff_rad, $
	tit='!6 Updraft speed vs. effective radius', $
	xtit='Maximum Updraft speed !8w!6 (cm/s)', $
	ytit='Effective Radius !8r!6!Ie!N (!7l!6)', $
	linestyle=1, thick=2.0,psym=6,/ynoz
oplot,w,reff_norad, $
	linestyle=0,thick=2.0,psym=1
oplot,abc,ord, $
	linestyle=0,thick=2.
;oplot,w,yf,linestyle=0,thick=2.0,psym=2
	
;dt_sng=string(Format='(I4)',dz)
xyouts,.50,0.6,'!9B!6   !6With Radiation',size=1.0,/NORMAL
xyouts,.50,.5,'+!6   !6Without Radiation',size=1.0,/NORMAL
xyouts,.50,.4,'!6     !6Least squares fit ',size=1.0,/NORMAL
oplot,[40,48],[34,34],linestyle=0,thick=2.
;oplot,[.5,.5],[.5,.5],linestyle=0,thick=2.0,/NORMAL
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

coeffs=svdfit(w,iwp_rad,m,weight=weight,YFIT=yf)
abc=findgen(100)*max(w)/99.
ord=poly(abc,coeffs)
print,'num_points = ',num_points
print,'poly degree = ',m-1
print,coeffs

plot,w,iwp_rad, $
	tit='!6 Updraft speed vs. ice water path', $
	xtit='!6Maximum Updraft speed !8w!6 (cm/s)', $
	ytit='!6IWP (g/m!E2!N)', $
	linestyle=1, thick=2.0,psym=6,/ynoz
oplot,w,iwp_norad, $
	linestyle=0,thick=2.0,psym=1
oplot,abc,ord, $
	linestyle=0,thick=2.
	
xyouts,0.20,0.8,'!9B!6   !6With Radiation',size=1.0,/NORMAL
xyouts,0.20,0.7,'+!6   !6Without Radiation',size=1.0,/NORMAL
xyouts,0.20,0.6,'!6     !6Least squares fit ',size=1.0,/NORMAL
oplot,[5,10],[360,360],linestyle=0,thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot,w,iwp_rad, $
	tit='!6 Updraft speed vs. ice water path', $
	xtit='!6Maximum Updraft speed !8w!6 (cm/s)', $
	ytit='!6IWP (g/m!E2!N)', $
	linestyle=1, thick=2.0,psym=6,/ynoz
oplot,w,iwp_norad, $
	linestyle=0,thick=2.0,psym=1
oplot,abc,ord, $
	linestyle=0,thick=2.
	
xyouts,0.20,0.8,'!9B!6   !6With Radiation',size=1.0,/NORMAL
xyouts,0.20,0.7,'+!6   !6Without Radiation',size=1.0,/NORMAL
xyouts,0.20,0.6,'!6     !6Least squares fit ',size=1.0,/NORMAL
oplot,[5,10],[360,360],linestyle=0,thick=2.
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

temperature=[-60.0,-55.,-50.0,-45.,-40.0,-35.,-30.0,-25.,-20]
temperature2=[-60.0,-55.,-50.0,-45.,-40.0,-35.,-30.0,-25.]
reff=[25.5,11.5,24.3,25.2,45.6,74.1,72.1,68.6]
kappa_bb=[.0548,0.101,.0563,.0544,.0336,.0214,.0218,.0228]

plot,temperature2,reff,/nodata, $
	tit='!6 Observations of warm frontal overrunnings with jets', $
	xtit='!6Temperature !8T!6 (C)', $
	ytit='!6Effective radius !8r!6!Ie!N (!7l!6)'
for i=0,7 do box,temperature2(i),0.0,temperature2(i)+5.,reff(i),1

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot,temperature2,reff,/nodata, $
	tit='!6 Observations of warm frontal overrunnings with jets', $
	xtit='!6Temperature !8T!6 (C)', $
	ytit='!6Effective radius !8r!6!Ie!N (!7l!6)'
for i=0,7 do box,temperature2(i),0.0,temperature2(i)+5.,reff(i),0
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

abc=findgen(80)+1.
ord=	.04*(.0020+1.118/abc)+ $
		.20*(.0016+1.166/abc)+ $
		.34*(.0003+1.338/abc)+ $
		.281*(.0068+0.600/abc)+ $
		.139*(.0036+1.136/abc)
plot,reff,kappa_bb, $
	tit='!6 Effective radius vs. LW mass absorption coefficient !7j!6', $
	xtit='!6Effective radius !8r!6!Ie!N (!7l!6)', $
	ytit='!7j!6 (m!E2!N/g)', $
	psym=6,thick=2.
oplot,abc,ord,linestyle=2,thick=2.
;oplot,[0.0,100.0],[.0602,0.0602],linestyle=1
;oplot,[10.0,10.0],[0.0,0.2],linestyle=1
;oplot,[30.0,30.0],[0.0,0.2],linestyle=0
;oplot,[0.0,100.0],[.045,.045],linestyle=0
oplot,[10.0],[.0602],psym=1,thick=2.
oplot,[30.0],[.045],psym=7,thick=2.
oplot,reff_rad,kappa,psym=5,thick=2.

xyouts,0.60,0.9,'!4D!6   !6This study WFOs',size=1.0,/NORMAL
xyouts,0.60,0.8,'!9B!6   !6Observed WFOs',size=1.0,/NORMAL
xyouts,0.60,0.7,'+!6   !6CCM2',size=1.0,/NORMAL
xyouts,0.60,0.6,'!9X!6   !6Recommended',size=1.0,/NORMAL
xyouts,0.60,.5,'!6    !6Ebert & Curry',size=1.0,/NORMAL
;oplot,[.5,0.6],[.5,.5],linestyle=2,/NORMAL,thick=2.
oplot,[40,48],[.058,.058],linestyle=2,thick=2.

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
plot,reff,kappa_bb, $
	tit='!6 Effective radius vs. LW mass absorption coefficient !7j!6', $
	xtit='!6Effective radius !8r!6!Ie!N (!7l!6)', $
	ytit='!7j!6 (m!E2!N/g)', $
	psym=6,thick=2.
oplot,abc,ord,linestyle=2,thick=2.
oplot,[10.0],[.0602],psym=1,thick=2.
oplot,[30.0],[.045],psym=7,thick=2.
oplot,reff_rad,kappa,psym=5,thick=2.

xyouts,0.60,0.9,'!4D!6   !6This study WFOs',size=1.0,/NORMAL
xyouts,0.60,0.8,'!9B!6   !6Observed WFOs',size=1.0,/NORMAL
xyouts,0.60,0.7,'+!6   !6CCM2',size=1.0,/NORMAL
xyouts,0.60,0.6,'!9X!6   !6Recommended',size=1.0,/NORMAL
xyouts,0.60,.5,'!6    !6Ebert & Curry',size=1.0,/NORMAL
;oplot,[.5,0.6],[.5,.5],linestyle=2,/NORMAL,thick=2.
oplot,[40,48],[.058,.058],linestyle=2,thick=2.
;******************** End Printing routine ******************************
	if ps eq 'n' then lpr else lps
   endif	
 endif ;endif pause

;***************enforced exit********************************
 goto,exit_gracefully

 exit_gracefully: foo=1
end

pro box,x0,y0,x1,y1,color
polyfill,[x0,x0,x1,x1],[y0,y1,y1,y0],col=color
end
