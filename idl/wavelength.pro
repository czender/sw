pro wavelength
;
erase
num_pts=500
axz_skip_nrm=0.15
num_x_axz=6
;
x_min=0.11
y_min=num_x_axz*axz_skip_nrm
x_max=0.95
y_max=0.90
plt_rgn_nrm=[x_min,y_min,x_max,y_max]
;
title_title='!5Wavelength, Wavenumber, and Frequency'
xyouts,.5,0.95,title_title,size=2.5,alignment=0.5,/NORMAL
;
; Plot the wavenumbers and frequencies for a given 
; wavelength spread
;
wavelength=0.5e-6+1.5e-6*findgen(num_pts)/(num_pts-1) ; m
wavenumber=1.0/(100.0*wavelength) ; cm^-1
frequency=3.0e8/wavelength ; s^-1
y_data=wavelength*0.
;
plot, $
	wavelength*1.0e6, $
	y_data, $
	/nodata, $
	/noerase, $
	/ynozero, $
	position=plt_rgn_nrm, $
	xstyle=13, $
	ystyle=12
;
axis, $
	x_min, $
	y_min-0.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xstyle=1, $
	xtick_get=x_tick_vals, $
	xticklen=0.02, $
	xticks=10, $
	xtitle='Wavelength !7k!5 (!7l!5m)'
;
;print,x_tick_vals
x_tick_vals=x_tick_vals/1.0e6 ; um -> m
wavenumber=1.0/(100.0*x_tick_vals) ; cm^-1
axis, $
	x_min, $
	y_min-1.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xrange=[max(wavenumber/1000.),min(wavenumber/1000.)], $
	xstyle=1, $
	xticklen=0.02, $
	xtickname=auto_sng(wavenumber/1000.0,2), $
	xticks=10, $
	xtitle='!5Wavenumber !7m!5 (1000 cm!E-1!N)'
;
frequency=3.e8/x_tick_vals ; s^-1
axis, $
	x_min, $
	y_min-2.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xrange=[max(frequency/1.0e12),min(frequency/1.0e12)], $
	xstyle=1, $
	xticklen=0.02, $
	xtickname=auto_sng(frequency/1.0e12,1), $
	xticks=10, $
	xtitle='!5Frequency !7m!5 (10!E12!N s!E-1!N)'
;
; Plot the wavelengths and frequencies for a given 
; wavenumber spread
;
wavenumber=5000.+15000.*findgen(num_pts)/(num_pts-1) ; cm^-1
plot, $
	wavenumber/1000.0, $
	y_data, $
	/nodata, $
	/noerase, $
	/ynozero, $
	position=plt_rgn_nrm, $
	xstyle=13, $
	ystyle=12
;
axis, $
	x_min, $
	y_min-3.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xstyle=1, $
	xtick_get=x_tick_vals, $
	xticklen=0.02, $
	xticks=10, $
	xtitle='!5Wavenumber !7m!5 (1000 cm!E-1!N)'
;
;print,x_tick_vals
x_tick_vals=x_tick_vals*1000. ; 1000 cm^-1 -> cm^-1
wavelength=1./(100.*x_tick_vals); m
axis, $
	x_min, $
	y_min-4.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xrange=[max(wavelength*1.0e6),min(wavelength*1.0e6)], $
	xstyle=1, $
	xticklen=0.02, $
	xtickname=auto_sng(wavelength*1.0e6,2), $
	xticks=10, $
	xtitle='Wavelength !7k!5 (!7l!5m)'
;
frequency=3.e8/wavelength ; s^-1
axis, $
	x_min, $
	y_min-5.*axz_skip_nrm, $
	/normal, $
	charsize=1.3, $
	xaxis=0, $
	xstyle=1, $
	xticklen=0.02, $
	xtickname=auto_sng(frequency/1.0e12,1), $
	xticks=10, $
	xtitle='!5Frequency !7m!5 (10!E12!N s!E-1!N)'
;
end
