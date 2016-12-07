; $Id$ -*-f90-*-

; Purpose: Common database of field properties for NCL scripts
; Sets fld_sng, units_sng, and, maybe, cnt_lvl, scl_fct

; Usage:
; if(.not.isvar("fld_prp_ncl")) then load "${HOME}/sw/ncl/fld_prp.ncl" end if

; Indicate to host program that header has been loaded
fld_prp_ncl=True
print("dbg: Loading ${HOME}/sw/ncl/fld_prp.ncl")

if (fld_nm.eq."AERAFRC_TOP") then 
   if(.not.isvar("fld_sng")) then fld_sng="Soot+Dust TOA Forcing" end if
   units_sng="W m~S~-2~N~"
end if; !AERAFRC_TOP
if (fld_nm.eq."ALBS") then 
   if(.not.isvar("fld_sng")) then fld_sng="Surface Albedo" end if
   units_sng="fraction"
end if; !ALBS
if (fld_nm.eq."alt") then 
   if(.not.isvar("fld_sng")) then fld_sng="Height" end if
   units_sng="m"
end if; !alt
if (fld_nm.eq."lev") then 
   if(.not.isvar("fld_sng")) then fld_sng="Pressure" end if
   units_sng="mb"
end if; !lev
if (fld_nm.eq."DUSTOD_v") then 
   if(.not.isvar("fld_sng")) then fld_sng="Dust Optical Depth" end if
   units_sng="fraction"
   cnt_lvl_md_xpt_ctl="ExplicitLevels"
   cnt_lvl=(/0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.75,1.0,1.5/) ; Used in ZBN03
   cnt_nbr_xpt_ctl=dimsizes(cnt_lvl)
end if; !DUSTOD_v
if (fld_nm.eq."FSNO") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Cover" end if
   units_sng="fraction"
end if; !FSNO
if (fld_nm.eq."H2OSNO") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Depth (LWE)" end if
   units_sng="mm"
end if; !H2OSNO
if (fld_nm.eq."j_NO2") then 
   if(.not.isvar("fld_sng")) then fld_sng="NO~B~2~N~ Photolysis Rate" end if
   units_sng="s~S~-1~N~"
end if; !j_NO2
if (fld_nm.eq."lmn_bb_aa_ndr") then 
   if(.not.isvar("fld_sng")) then fld_sng="Luminance towards Nadir" end if
   units_sng="lm m~S~-2~N~ sr~S~-1~N~"
end if; !j_NO2
if (fld_nm.eq."QMELT") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Melt" end if
   units_sng="mm day~S~-1~N~"
   scl_fct=86400 ; mm s-1 -> mm day-1
end if; !QMELT
if (fld_nm.eq."phz_fnc_ffc") then 
   if(.not.isvar("fld_sng")) then fld_sng="Effective Phase Function" end if
   units_sng="sr~S~-1~N~"
end if; !phz_fnc_ffc
if (fld_nm.eq."rfl_spc_snw") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Albedo" end if
   units_sng="fraction"
end if; !rfl_spc_snw
if (fld_nm.eq."alb_spc_snw") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Albedo" end if
   units_sng="fraction"
end if; !alb_spc_snw
if (fld_nm.eq."alb_spc_snw_crc") then 
   if(.not.isvar("fld_sng")) then fld_sng="Corrected Snow Albedo" end if
   units_sng="fraction"
end if; !alb_spc_snw_crc
if (fld_nm.eq."foo") then 
   if(.not.isvar("fld_sng")) then fld_sng="Aerosol Snowpack Forcing" end if
   units_sng="W m~S~-2~N~"
end if; !foo
if (fld_nm.eq."PRECT") then 
   if(.not.isvar("fld_sng")) then fld_sng="Total Precipitation" end if
   units_sng="mm day~S~-1~N~"
;   scl_fct=86400000 ; m s-1 -> mm day-1
end if; !SNOAERFRC
if (fld_nm.eq."SNOAERFRC") then 
   if(.not.isvar("fld_sng")) then fld_sng="Aerosol Snowpack Forcing" end if
   units_sng="W m~S~-2~N~"
end if; !SNOAERFRC
if (fld_nm.eq."SNOBCC") then 
   if(.not.isvar("fld_sng")) then fld_sng="Surface Snow BC" end if
   units_sng="ng g~S~-1~N~"
   cnt_lvl_md_xpt_ctl="ExplicitLevels"
   cnt_lvl=(/1,2,5,10,20,50,100,200,500,1000/) ; Approximately used in FZR07
   cnt_nbr_xpt_ctl=dimsizes(cnt_lvl)
end if; !SNOBCC
if (fld_nm.eq."SNODSTC") then 
   if(.not.isvar("fld_sng")) then fld_sng="Surface Snow Dust" end if
   units_sng="~0714~g g~S~-1~N~"
   cnt_lvl_md_xpt_ctl="ExplicitLevels"
   cnt_lvl=(/1,2,5,10,20,50,100,200,500,1000/) ; Approximately used in FZR07
   cnt_nbr_xpt_ctl=dimsizes(cnt_lvl)
end if; !SNODSTC
if (fld_nm.eq."SNODSTFRC") then 
   if(.not.isvar("fld_sng")) then fld_sng="Dust-Snowpack Forcing" end if
   units_sng="W m~S~-2~N~"
   cnt_lvl_md_xpt_ctl="ExplicitLevels"

   cnt_lvl=(/0.01,0.02,0.05,0.1,0.2,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0/) ; Approximately used in ZFM08
   cnt_nbr_xpt_ctl=dimsizes(cnt_lvl)
end if; !SNODSTFRC
if (fld_nm.eq."SNODYSL") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Days since Jan. 1" end if
   units_sng="day"
end if; !SNODYSL
if (fld_nm.eq."SNORDS") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Effective Radius" end if
   units_sng="~0714~m"
end if; !SNORDS
if (fld_nm.eq."SNOSOTFRC") then 
   if(.not.isvar("fld_sng")) then fld_sng="Soot-Snowpack Forcing" end if
   units_sng="W m~S~-2~N~"
end if; !SNOSOTFRC
if (fld_nm.eq."SNOWAGE") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Age" end if
   units_sng="fxm"
end if; !SNOWAGE
if (fld_nm.eq."SNOWDP") then 
   if(.not.isvar("fld_sng")) then fld_sng="Snow Depth" end if
   cnt_lvl_md_xpt_ctl="ExplicitLevels"
   ;	cnt_lvl=(/0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5/) ; Used in fxm
   cnt_lvl=(/0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ; Used in fxm
   cnt_nbr_xpt_ctl=dimsizes(cnt_lvl)
   units_sng="m"
end if; !SNOWDP
if (fld_nm.eq."TREFHT") then 
   if(.not.isvar("fld_sng")) then fld_sng="Surface Temperature" end if
   units_sng="K"
end if; !TREFHT
