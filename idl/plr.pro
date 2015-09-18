pro plr_gph,rds,ngl,_extra=extra

; Fake data if needed.
if n_params() eq 0 then begin
   ngl=((randomu(seed,360)*360)-180)*!pi/180.0
   rds=randomu(seed,360)
endif; endif 

; Load plot colors
tvlct,[100,255,0],[100,255,255],[100,0,0],1
;device,decomposed=0

; Establish plot coordinates
plot,rds,ngl,/polar,xstyle=5,ystyle=5,background=1,position=aspect(1.0),/nodata

; Draw axes through center
axis,/xaxis,0,0,color=2
axis,/yaxis,0,0,color=2

; Plot data
oplot,rds,ngl,color=3,psym=2,/polar

; Draw 25 and 75 percent circles
max_dat=max(rds)
pct_25=circle(0,0,0.25*max_dat)
pct_75=circle(0,0,0.75*max_dat)
plots,pct_25,color=2
plots,pct_75,color=2

end; end plr_gph()
