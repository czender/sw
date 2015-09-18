;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: poisson.pro,v $
; $Source: /home/zender/cvs/idl/poisson.pro,v $
; $Id$
; $State: Exp $
;
; NB: get CVS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: compute sub-gridscale variability of poisson distribution
; updrafts on IWC for CCM.
;
; Usage: poisson,mu=0.5,mean_x=0.02,max_x=.2,pause='n',t_graph=253.15,w_graph=0.02
;
; $Log: not supported by cvs2svn $
; Revision 1.7  2000-01-15 02:07:54  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/01 01:55:53  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 00:18:29  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:05  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:50  zender
; Imported sources
;
; Revision 1.1  1994/05/23  23:41:50  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
function Hey77_IWC_of_T_w,T,w
; given input T in Kelvin and w in m/s this routine
; implements the Hey77 pzn. of IWC. which takes the updraft
; speed in cm/s
;
w_cgs=w*100.
;
;if (w gt 0. and T lt 273.15) then begin
	tmp1=1.59/w_cgs^.04
	tmp2=-0.01*w_cgs^.186
	tmp3=tmp2*(273.15-T)^tmp1
	tmp4=0.072*w_cgs^.78
	Hey77_IWC=tmp4*exp(tmp3)/1000. ;converts g/m^3 to kg/m^3
;endif else begin
;	Hey77_IWC=0.
;endelse
;
return,Hey77_IWC
end
;
function alt_basis,x,m
common poisson,		mean_w
n=n_elements(x)
r=fltarr(n,m)
for i=0,n-1 do begin
	r(i,0)=(mean_w^.78)*1.
	r(i,1)=(mean_w^.78)*x(i)
	r(i,2)=(mean_w^.78)*x(i)*x(i)
	r(i,3)=(mean_w^.78)*x(i)^3
;	r(i,0)=1.
;	r(i,1)=1./x(i)
;	r(i,2)=1./x(i)^2.
;	r(i,1)=exp((6.9e-6*x(i)^2.))
;	r(i,1)=exp(1.0e-1*(x(i)+20))
;	r(i,2)=exp(x(i)*x(i)/400.)
endfor
return,r
end
;
function sub_grid_pzn,abc
common poisson,		mean_w
;coeffs=[]
coeffs=[5.55341,     0.282799,   0.00471977,  2.57237e-05]
;coeffs=[5.60648,     0.289417,   0.00485748,  2.57188e-05]
ord=(mean_w^.78)*(coeffs(0)+coeffs(1)*abc+ $
	coeffs(2)*abc*abc+ $
	coeffs(3)*abc^3)
return,ord	
end
;
function alt_basis_ords,x,coeffs
common poisson,		mean_w
n=n_elements(x)
m=n_elements(coeffs)
y=fltarr(n)
for i=0,n-1 do begin
	y(i)=(mean_w^.78)*(coeffs(0)*1.+coeffs(1)*x(i)+coeffs(2)*x(i)^2+coeffs(3)*x(i)^3)
;	y(i)=mean_w*(coeffs(0)*1.+coeffs(1)*x(i)+coeffs(2)*x(i)^2)
;	y(i)=coeffs(0)*1.+coeffs(1)/x(i)+coeffs(2)/x(i)^2.
;	y(i)=coeffs(0)*1.+coeffs(1)/x(i)
;	y(i)=coeffs(0)*1.+coeffs(1)*exp((6.9e-6*x(i)^2.))
;	y(i)=coeffs(0)*1.+coeffs(1)*exp(1.0e-1*(x(i)+20))
;	y(i)=coeffs(0)*1.+coeffs(1)*exp(-x(i))+coeffs(2)*exp(-x(i)*x(i))
endfor
;y=poly(abc,coeffs)
return,y
end
;
pro poisson, mu=mu, num_abs=num_abs, max_x=max_x, mean_x=mean_x, $
	min_temp=min_temp, tpt_nbr=tpt_nbr, temp_ncr=temp_ncr, $
	w_graph=w_graph, t_graph=t_graph, pause=pause, num_coeffs=num_coeffs
;
; The program will scale num_abs wind speeds between 0. and max_x
; with a scaling factor which is the ratio of mu to mean_x. So mu
; should be chosen to give the poisson distribution the right shape,
; and then once the abc are generated and scaled they are sent
; to the poisson distribution function to get the probabilities.
; Since i am throwing real number at the Pp function the results will
; not be normalized to 1, the Pp function is only normalized for 
; integers. 
;
common poisson,		mean_w
;
if n_elements(mu) eq 0 then mu = .5			
if n_elements(num_abs) eq 0 then num_abs = 100
if n_elements(max_x) eq 0 then max_x = .2		;m/s
if n_elements(mean_x) eq 0 then mean_x = .02		;m/s
if n_elements(min_temp) eq 0 then min_temp = 220.
if n_elements(num_coeffs) eq 0 then num_coeffs = 3      ;num_coeffs-1 is the degree of the polynomial, num_coeffs is the number of coeffs in fitting func.
if n_elements(tpt_nbr) eq 0 then tpt_nbr = 50
if n_elements(temp_ncr) eq 0 then temp_ncr = 1.
if n_elements(pause) eq 0 then pause = 'y'
if n_elements(w_graph) eq 0 then w_graph = .02
if n_elements(t_graph) eq 0 then t_graph = 253.15
if n_elements(foo) eq 0 then foo = 1.
if pause eq 'y' then !P.multi = 0 else begin
    !P.multi = [0, 2, 2, 0, 0]  ;=[0,num_cols,num_rows,0,0]
    pause = 'n'
endelse 
;
mean_w=mean_x
;
; interface quantities valid at all indices.
; midpoint quantities: xmid(i)=0.5*(x_int(i-1) + x_int(i)) so valid 
; anywhere except index 0.
;
x_int = max_x*findgen(num_abs+1)/num_abs
x_mid = .5*(x_int+shift(x_int,1))
delta_x = x_int-shift(x_int,1)
;
scale = mu/mean_x
;
x_mid_scaled = x_mid*scale
x_int_scaled = x_int*scale
delta_x_scaled = delta_x*scale
;
; Find the probability associated with each windspeed
; Code the Poisson distribution function as the exponential of
; the difference of two logarithms in order to avoid overflows.
;
prob = lngamma(x_mid_scaled+1.)
prob = -mu+x_mid*alog(mu)-prob ; NB: use x_mid here, not x_mid_scaled (why?)
prob = exp(prob)
;
; Figure out the normalization integral
total_prob=total(prob(1:num_abs)*delta_x_scaled(1:num_abs))
;
print,'mu = ',mu
print,'mean_x = ',mean_x*100.0,' cm/s'
print,'scale = ',scale
print,'max_x = ',max_x*100.0,' cm/s'
print,'max_x_scaled = ',max_x*scale
print,'total_prob = ',total_prob
print_sng='Hey77_IWC(t = '+auto_sng(t_graph-273.15,1)+' C, w = '+ $
	auto_sng(w_graph*100,1)+' cm/s) = '+ $
	auto_sng(1000.*Hey77_IWC_of_T_w(t_graph,w_graph),5)+' g/m^3'
print,print_sng
;
;print,'x_mid = ',x_mid
;print,'x_int = ',x_int
;print,'delta_x = ',delta_x
;print,'prob = ',prob
;
; plot the probability distribution, i.e., the probability of measuring
; an updraft speed between w-dw/2 and w+dw/2 centimeters per second 
; is the y-axis value times dw in centimeters per second.
; ummm, er, this isn't actually true....
; what's graphed is certainly not the probability distribution...
;
plot,x_mid(1:num_abs)*100.0,prob(1:num_abs), $
	tit='!5Poisson Distribution of Updraft Speeds', $
	xtit='!5Updraft !8w!5 (cm!E-1!Ns)', xsty=0, $
	ytit='!5Probability !8P!Ip!N(!7l!5,!8w!5)',ysty=0, $
	thick=3.0, charsize=2.
;
xyouts,.5,0.8,'!5Mean !7l!5 = '+auto_sng(mu,2),size=1.75,/NORMAL
xyouts,.5,0.7,'!5Mean !8w!5 = '+auto_sng(mean_x*100.0,1)+' cm/s',size=1.75,/NORMAL
xyouts,.5,0.6,'!5scale = '+auto_sng(scale,1),size=1.75,/NORMAL
;
 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** end Printing routine ******************************
	lpr
   endif	
 endif ;endif pause
;
; Now loop over all the temperatures to create the IWC bar of T function 
temp=min_temp+indgen(tpt_nbr+1)*temp_ncr
IWC_bar=fltarr(tpt_nbr+1)
Hey77_IWC=fltarr(tpt_nbr+1,num_abs+1)
for counter=0,tpt_nbr do begin
	Hey77_IWC(counter,*)=Hey77_IWC_of_T_w(temp(counter),x_mid)
	IWC_bar(counter)= total($
		prob(1:num_abs)*Hey77_IWC(counter,1:num_abs)* $
		delta_x_scaled(1:num_abs))/ $
		total_prob
endfor
;
plot,temp-273.15,IWC_bar*1000.0, $
	tit='!5Mean Grid Box !8IWC(T)!5', $
	xtit='!5Temperature !8T!5 (!E!12_!N!5C)', xsty=0, $
	ytit='!8IWC !5(g-m!E-3!N)',ysty=0, $
	thick=3.0, charsize=2.
;
xyouts,.5,0.8,'!5Mean !7l!5 = '+auto_sng(mu,2),size=1.75,/NORMAL
xyouts,.5,0.7,'!5Mean !8w!5 = '+auto_sng(mean_x*100.0,1)+' cm/s',size=1.75,/NORMAL
xyouts,.5,0.6,'!5scale = '+auto_sng(scale,1),size=1.75,/NORMAL
;
yf=fltarr(tpt_nbr+1)
weight=replicate(1.0,tpt_nbr+1)
coeffs=svdfit(temp-273.15,IWC_bar*1000.0,num_coeffs,weight=weight,YFIT=yf)
;coeffs=svdfit(temp-273.15,IWC_bar*1000.0,num_coeffs,weight=weight,funct='alt_basis',YFIT=yf)
abc=temp-273.15
ord=poly(abc,coeffs)
;ord=alt_basis_ords(abc,coeffs)
print,'num_points = ',tpt_nbr+1
print,'number basis functions = ',num_coeffs
print,'IWC_bar(T) coeffs = ',coeffs
oplot,abc,ord, $
        linestyle=1,thick=2.
oplot,abc,yf,linestyle=0,thick=2.0,psym=2
ord=sub_grid_pzn(abc)
oplot,abc,ord, $
        linestyle=2,thick=2.
;
 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** end Printing routine ******************************
	lpr
   endif	
 endif ;endif pause
;
; Now plot a single graph of IWC(T) for a given w
plot,temp-273.15,Hey77_IWC_of_T_w(temp(*),w_graph)*1000.0, $
	tit='!8IWC(T)!5 for !8w!5 = '+auto_sng(100.*w_graph,1)+' cm/s', $
	xtit='!5Temperature !8T!5 (!E!12_!N!5C)', xsty=0, $
	ytit='!8IWC!5 (g-m!E-3!N)',ysty=0, $
	thick=3.0, charsize=2.
;
xyouts,.5,0.8,'!8w!5 = '+auto_sng(100.*w_graph,1)+' cm/s',size=1.75,/NORMAL
;
 if pause eq 'y' then begin
   print,' Hit any key to continue, p to print this graph, or q to quit ...'
   junk = get_kbrd(1)
   if junk eq 'q' then goto, exit_gracefully
   if junk eq 'p' then begin
	halfplot
;******************** Printing routine, should be identical to above ***************
;******************** end Printing routine ******************************
	lpr
   endif	
 endif ;endif pause
;
; Now plot a single graph of IWC(w) for a given T
plot,x_mid(1:num_abs)*100.0,Hey77_IWC_of_T_w(t_graph,x_mid(1:num_abs))*1000.0, $
	tit='!8IWC(w)!5 for !8T!5 = '+auto_sng(t_graph-273.15,1)+' !E!12_!N!5C', $
	xtit='!5Updraft !8w!5 (cm-s!E-1!N)', xsty=0, $
	ytit='!8IWC!5 (g-m!E-3!N)',ysty=0, $
	thick=3.0, charsize=2.
;
xyouts,.5,0.8,'!8T!5 = '+auto_sng(t_graph-273.15,1)+' !E!12_!N!5C',size=1.75,/NORMAL
;
;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1
end

