; $Id$

; Purpose: NCL routines to support thesis work

; Example usage:
; load "/home/zender/anv/anv.ncl"
; ncl < "/home/zender/anv/anv.ncl"
; /opt/ncarg4.0/sp/bin/ncl < "/home/zender/anv/anv.ncl"

procedure crd_rvr()

; This routine reverses the level coordinate from the input file to the output file

local idx_fl,idx_var

begin ; begin crd_rvr()

yr_mth_sng=(/"7993_1202","7993_0305","7993_0608","7993_0911"/)
prg_nm="crd_rvr"
fl_nbr=dimsizes(yr_mth_sng)

; Initialize loop variables
idx_fl=0
idx_var=0
dbg_lvl=1

do idx_fl=0,fl_nbr-1

fl_in="/data2/zender/obsst01/obsst01_prs_"+yr_mth_sng(idx_fl)+".nc"
id_in=addfile(fl_in,"w")
if(ismissing(id_in)) then print(prg_nm+": ERROR file "+fl_in+" does not exist") end if
print("Processing file "+fl_in+"...")

lev=id_in->lev
nbr_lev=dimsizes(lev)

T=id_in->T(:,:,nbr_lev-1:0,:)
id_in->T=T

; Redefine the corrupted level values in the output file
delete(lev)
lev=id_in->lev
lev=(/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000/)
id_in&lev=lev

ilev=id_in->ilev
ilev=(/1,15,25,40,60,85,125,175,225,275,350,450,550,650,775,887.5,962.5,1000/)
id_in&ilev=ilev

delete(id_in)

end do ; end loop over fl

end ; end crd_rvr()

procedure lvl_fix(fl_in)

; This opens a file and replaces the corrupt level array with a fixed one

begin ; begin lvl_fix()

prg_nm="lvl_fix"

; Initialize loop variables
dbg_lvl=1

id_in=addfile(fl_in,"w")
if(ismissing(id_in)) then print(prg_nm+": ERROR file "+fl_in+" does not exist") end if

; Redefine the corrupted level values in the output file

lev=new((/17/),float)
lev=(/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000/)
;lev=(/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,150.,100., 70., 50., 30., 20., 10./)
id_in&lev=lev

ilev=new((/18/),float)
ilev=(/1,15,25,40,60,85,125,175,225,275,350,450,550,650,775,887.5,962.5,1000/)
;;ilev=(/1000,962.5,887.5,775,650,550,450,350,275,225,175,125,85,60,40,25,15,1/)
id_in&ilev=ilev

delete(id_in)

end ; end lvl_fix()

; Have this procedure call itself for easy batch use
;lvl_fix("/data2/zender/obsst01/obsst01_prs_7993_1202.nc")
;lvl_fix("/data2/zender/obsst01/obsst01_prs_7993_0305.nc") 
;lvl_fix("/data2/zender/obsst01/obsst01_prs_7993_0608.nc")
;lvl_fix("/data2/zender/obsst01/obsst01_prs_7993_0911.nc")


