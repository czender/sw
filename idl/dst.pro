$Id$

; Purpose: Offline analysis routines to look at Dust data

@/home/zender/idl/ibp_clr.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PaG77 routines: Compute number distributions from PaG77
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function rds_avg_sfc2rds_avg_nbr,rds_avg_sfc,gsd
; Purpose: Convert surface mean radius to number mean radius
return,rds_avg_sfc*exp(-2.0*alog(gsd)*alog(gsd))
end; end rds_avg_sfc2rds_avg_nbr()

function rds_avg_vlm2rds_avg_nbr,rds_avg_vlm,gsd
; Purpose: Convert volume mean radius to number mean radius
return,rds_avg_vlm*exp(-3.0*alog(gsd)*alog(gsd))
end; end rds_avg_vlm2rds_avg_nbr()

pro PaG77
rds_avg_sfc_c=[0.11,0.07]
rds_avg_sfc_a=[1.31,1.87,1.65,0.58,0.37,1.49,1.42,1.20,1.26]
rds_avg_sfc_b=[31.9,31.8,15.6,21.1,14.7]

rds_avg_vlm_c=[0.14,0.13]
rds_avg_vlm_a=[2.29,2.84,3.08,1.25,0.94,3.35,3.62,2.56,2.17]
rds_avg_vlm_b=[35.5,35.1,20.7,34.0,18.8]

gsd_c=[1.56,2.20]
gsd_a=[2.11,1.90,2.20,2.41,2.6,2.46,2.63,2.38,2.09]
gsd_b=[1.38,1.37,1.70,2.0,1.65]

rds_avg_nbr_c_sfc=rds_avg_sfc2rds_avg_nbr(rds_avg_sfc_c,gsd_c)
rds_avg_nbr_b_sfc=rds_avg_sfc2rds_avg_nbr(rds_avg_sfc_b,gsd_b)
rds_avg_nbr_a_sfc=rds_avg_sfc2rds_avg_nbr(rds_avg_sfc_a,gsd_a)

rds_avg_nbr_c_vlm=rds_avg_vlm2rds_avg_nbr(rds_avg_vlm_c,gsd_c)
rds_avg_nbr_b_vlm=rds_avg_vlm2rds_avg_nbr(rds_avg_vlm_b,gsd_b)
rds_avg_nbr_a_vlm=rds_avg_vlm2rds_avg_nbr(rds_avg_vlm_a,gsd_a)

print,'PaG77 statistics from surface modes:'
dst_prn,rds_avg_nbr_c_sfc,rds_avg_nbr_c_vlm,rds_avg_sfc_c,rds_avg_vlm_c,gsd_c
dst_prn,rds_avg_nbr_a_sfc,rds_avg_nbr_a_vlm,rds_avg_sfc_a,rds_avg_vlm_a,gsd_a
dst_prn,rds_avg_nbr_b_sfc,rds_avg_nbr_b_vlm,rds_avg_sfc_b,rds_avg_vlm_b,gsd_b

end; end PaG77()

pro dst_prn,rds_avg_nbr_sfc,rds_avg_nbr_vlm,rds_avg_sfc,rds_avg_vlm,gsd
dst_nbr=n_elements(gsd)
print,"Number(S)	Number(V)	Surface		Volume	GSD"
for idx=0,dst_nbr-1 do begin
print,rds_avg_nbr_sfc(idx),rds_avg_nbr_vlm(idx),rds_avg_sfc(idx),rds_avg_vlm(idx),gsd(idx)
endfor; end loop over idx
end; end dst_prn()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End PaG77 routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure dst_bdg: Graphs various dust mass budget diagnostics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dst_bdg, $
	fl_nm=fl_nm, $
	info=info, $
	prn=prn
; Usage:
; dst_bdg,fl_nm='/fs/cgd/data0/zender/dst/dst_mss_bdg.nc'
; dst_bdg,fl_nm='/fs/cgd/data0/zender/dst/dst04_bdg.nc'

@/home/zender/idl/ibp_clr.com
!p.multi=0

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm='/fs/cgd/data0/zender/dst/dst_mss_bdg.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'rec')
ncdf_diminq,nc_id,dim_id,foo,rec_nbr
if dbg then print,'Current number of records is '+rec_nbr
ncdf_varget,nc_id,'dst_mpc_mdl',dst_mpc_mdl
ncdf_varget,nc_id,'dst_mpc_a_fix',dst_mpc_a_fix
ncdf_varget,nc_id,'dst_mpc_a_tfilt',dst_mpc_a_tfilt
ncdf_varget,nc_id,'dst_mpc_b_tphys',dst_mpc_b_tphys
ncdf_varget,nc_id,'dst_mpc_dgn',dst_mpc_dgn
ncdf_varget,nc_id,'dst_mpc_dlt',dst_mpc_dlt
ncdf_varget,nc_id,'dst_mbl_dlt',dst_mbl_dlt
ncdf_varget,nc_id,'dst_pcp_dlt',dst_pcp_dlt
ncdf_varget,nc_id,'dst_trb_dlt',dst_trb_dlt
ncdf_varget,nc_id,'dst_grv_dlt',dst_grv_dlt
ncdf_varget,nc_id,'nstep',nstep
;ncdf_varget,nc_id,'foo',foo
ncdf_close,nc_id

abc=nstep
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=dst_mpc_mdl*1.0e9 ; [kg m-2] -> [ug m-2]
ord_min=min(ord)
ord_max=max(ord)
rng_y=[ord_min,ord_max]
;rng_y=[0,100]
chr_sz=2.0
ttl=''
x_ttl='!5Step !8t!5'
y_ttl='!5MPC (!7l!5g m!e-2!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out='/data/zender/ps/dst_bdg.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot dst_bdg
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
idx_orn=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(idx_orn)=255
g_curr(idx_orn)=128
b_curr(idx_orn)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

fld_nbr=4
lbl_sng=strarr(fld_nbr)
ln_sty=0.0*intarr(fld_nbr)
clr=intarr(fld_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
clr(3)=yll_idx
;clr(4)=idx_orn
lbl_sng(0)='dst_mpc_mdl'
lbl_sng(1)='dst_mpc_dgn'
lbl_sng(2)='dst_mpc_dlt'
lbl_sng(3)='dst_mpc_a_fix'
;lbl_sng(4)='dst_mpc_b_tphys'

oplot,abc,dst_mpc_mdl*1.0e9,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
oplot,abc,dst_mpc_dgn*1.0e9,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(1)
oplot,abc,dst_mpc_dlt*1.0e9,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(2)
oplot,abc,dst_mpc_a_fix*1.0e9,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(3)
;oplot,abc,dst_mpc_b_tphys*1.0e9,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(4)

if info then begin
ln_lgn_x1=0.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fl_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end dst_bdg()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure dst_bdg
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

