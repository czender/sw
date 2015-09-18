;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: sasha.pro,v $
; $Source: /home/zender/cvs/idl/sasha.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: procedure to compare sasha madronich's method of computing
; actinic flux enhancement with cloud droplets with my own.
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-01-10 23:36:35  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/01 01:55:53  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 02:09:43  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 00:18:29  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:51  zender
; Imported sources
;
; Revision 1.1  1993/07/28  20:10:44  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function incident_angle,z,radius
theta_i=asin(z/radius)
return,theta_i
end
;
function trans_angle,theta_i,n_i,n_t
tmp=n_i*sin(theta_i)/n_t
theta_t=asin(tmp)
return,theta_t
end
;
function straight_chord,theta_i,radius
L_not=2.*radius*cos(theta_i)
return,L_not
end
;
function bent_chord,theta_t,radius
L=2.*radius*cos(theta_t)
return,L
end
;
function fresnel,theta_i,theta_t
; Takes input angles theta_incident and theta_transmitted in radians and
; returns the fraction of an incoming beam lost to Fresnel's reflection.
;
diff=theta_i-theta_t
sum=theta_i+theta_t
num_1=(sin(diff))^2
denom_1=2.*(sin(sum))^2
num_2=(tan(diff))^2
denom_2=2.*(tan(sum))^2
reflectance=num_1/denom_1 + num_2/denom_2
return,reflectance
end
;
function absorption,beta,L
if beta eq 0. then begin
	alpha=1.
endif else begin
	tmp=1.-exp(-beta*L)
	alpha=tmp/(beta*L)
endelse
return,alpha
end
;
function eta_fct,alpha,L,L_not,R
eta=alpha*L*(1.-R)/(L_not*(1.-alpha*R))
return,eta
end
;
function zender_eta_fct,alpha,L,L_not,R
eta=alpha*L*(1.-R)/(L_not-alpha*R*L)
return,eta
end
;
function mean_eta_fct,radius,eta,z
dz=z-shift(z,1)
integrand=2.*!pi*z*eta*dz
mean_eta=total(integrand)/(!pi*radius*radius)
return,mean_eta
end
;
pro sasha, $
	filename=filename, $
	abc_nbr=abc_nbr, $
	pause=pause, $
	ps=ps, $
	radius=radius, $
	n_i=n_i, $
	n_t=n_t, $
	beta=beta
;
;Example usage:
;
;sasha,ps='n'
;
if n_elements(filename) eq 0 then filename = 'clouds.nc'
if n_elements(abc_nbr) eq 0 then abc_nbr = 1000
if n_elements(radius) eq 0 then radius = 1.
if n_elements(n_i) eq 0 then n_i = 1.
if n_elements(n_t) eq 0 then n_t = 1.33
if n_elements(beta) eq 0 then beta = 0.
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
if pause eq 'y' then !p.multi=0 else begin
        !p.multi=[0,3,1,0,0]		;=[0,num_cols,num_rows,0,0]
        pause = 'n'
endelse
;
sys_time = systime(1) 
close,/all
erase
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the enhancement factor eta for all impact parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
z=radius*findgen(abc_nbr)/(abc_nbr-1.0)
;z=0.9+0.1*radius*findgen(abc_nbr)/(abc_nbr-1.0)
z=z>0.1/abc_nbr
z=z<(1.0-0.1/abc_nbr)
theta_i=incident_angle(z,radius)
theta_t=trans_angle(theta_i,n_i,n_t)
L_not=straight_chord(theta_i,radius)
L=bent_chord(theta_t,radius)
R=fresnel(theta_i,theta_t)
alpha=absorption(beta,L)
eta=eta_fct(alpha,L,L_not,R)
zender_eta=zender_eta_fct(alpha,L,L_not,R)
mean_eta=mean_eta_fct(radius,eta,z)
;
;print,"z = ",z
;print,"theta_i = ",theta_i
;print,"theta_t = ",theta_t
;print,"L_not = ",L_not
;print,"L = ",L
;print,"R = ",R
;print,"alpha = ",alpha
;print,"eta = ",eta
;print,"zender_eta = ",zender_eta
print,"mean_eta = ",mean_eta
p=1.-exp(-2.*radius*beta)
print,"p = absorption per diametrical crossing =",p
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Graph it
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot_io, $
	z(*), $
	eta(*), $
        tit='!6Impact Parameter vs. Enhancement Factor', $
        xtit='!6Impact Parameter !8z!6', $
        ytit='!6Enhancement Factor !7g!6', $
;        xrange=[0.1,1.0], $
        linestyle=0, $
	TICKLEN=1.0, $
	thick=2.
;
oplot, $
	z(*), $
	zender_eta(*), $
        linestyle=1, $
	thick=2.
;
plots, $
	[0.0,1.0], $
	[mean_eta,mean_eta], $
        linestyle=2, $
	thick=2.0, $
	/data
;
ln_lgn_x1=.22
ln_lgn_dx=0.05
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=.3
txt_lgn_size=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!6Sasha !7g!6(!8z!6)',size=txt_lgn_size,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!6Zender !7g!6(!8z!6)',size=txt_lgn_size,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!6Mean !7g!6(!8z!6)',size=txt_lgn_size,/NORMAL
;
data_lgn_x=0.6
;
beta_sng=string(Format='(f6.2)',beta)
n_i_sng=string(Format='(f6.2)',n_i)
n_t_sng=string(Format='(f6.2)',n_t)
p_sng=string(Format='(f6.2)',p)
mean_eta_sng=string(Format='(f6.2)',mean_eta)
;
xyouts,data_lgn_x,lgn_y(0),'!6!8n!Iair!N!6 = '+n_i_sng,size=txt_lgn_size,/NORMAL
xyouts,data_lgn_x,lgn_y(1),'!6!8n!Idrop!N!6 = '+n_t_sng,size=txt_lgn_size,/NORMAL
xyouts,data_lgn_x,lgn_y(2),'!6Absorption !7b!6 = '+beta_sng,size=txt_lgn_size,/NORMAL
xyouts,data_lgn_x,lgn_y(3),'!6Diameter Abs. !8p!6 = '+p_sng,size=txt_lgn_size,/NORMAL
xyouts,data_lgn_x,lgn_y(4),'!6Mean Enhanc. !7g!6 = '+mean_eta_sng,size=txt_lgn_size,/NORMAL
;
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
;
;***************enforced exit********************************
goto,exit_gracefully
;
exit_gracefully: foo=1
;
end
