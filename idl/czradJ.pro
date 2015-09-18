;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: czradJ.pro,v $
; $Source: /home/zender/cvs/idl/czradJ.pro,v $
; $Id$
; $State: Exp $
;
; Purpose: procedures to graph photodissociation work.
;
; $Log: not supported by cvs2svn $
; Revision 1.2  2000-01-10 19:33:34  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:47  zender
; Imported sources
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro czradJ,num_layer=num_layer,pause=pause,ps=ps
;
;Example usage:
;
;czradJ
;
if n_elements(num_layer) eq 0 then num_layer = 18
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if pause eq 'y' then !p.multi=0 else begin
        !p.multi=[0,3,1,0,0];=[0,num_cols,num_rows,0,0]
        pause = 'n'
endelse

sys_time = systime(1) 

close,/all
openr,1,'czradJ.dat'
foo=""
readf,1,foo
readf,1,foo
pressure=fltarr(num_layer)
altitude=fltarr(num_layer)
czJ1D=fltarr(num_layer)
czJ3P=fltarr(num_layer)
pjrJ1D=fltarr(num_layer)
pjrJ3P=fltarr(num_layer)

for layer=0,num_layer-1 do begin
	readf,1,a,b,c,d,e,f
	pressure(layer)=a
	altitude(layer)=b
	czJ1D(layer)=c
	czJ3P(layer)=d
	pjrJ1D(layer)=e
	pjrJ3P(layer)=f
endfor
close,1

print,'num_layer = ',num_layer 

plot_oi,czJ1D(*),altitude(*), $
        tit='!6 Photodissociation vs. Altitude', $
        xtit='J (s!E-1!N)', $
        ytit='Altitude (km)', $
	xrange=[1.0e-5,1.0e-2], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,czJ3P(*),altitude(*), $
        linestyle=0, thick=2. 
oplot,pjrJ1D(*),altitude(*), $
        linestyle=1, thick=2.  
oplot,pjrJ3P(*),altitude(*), $
        linestyle=1, thick=2.  

 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
        halfplot
;******************** Printing routine, should be identical to above ***************
plot_oi,czJ1D(*),altitude(*), $
        tit='!6 Photodissociation vs. Altitude', $
        xtit='J (s!E-1!N)', $
        ytit='Altitude (km)', $
	xrange=[1.0e-5,1.0e-2], $
        linestyle=0, TICKLEN=1.0, thick=2.
oplot,czJ3P(*),altitude(*), $
        linestyle=0, thick=2. 
oplot,pjrJ1D(*),altitude(*), $
        linestyle=1, thick=2.  
oplot,pjrJ3P(*),altitude(*), $
        linestyle=1, thick=2.  
;******************** End Printing routine ******************************
        if ps eq 'n' then lpr else lps
   endif        
 endif ;endif pause

;***************enforced exit********************************
 goto,exit_gracefully

 exit_gracefully: foo=1

end
