; $Id$ -*-f90-*-

; Purpose: Functions and named constants used by my NCL scripts

; Usage:
; Insert following at top of NCL script
; if(.not.isvar("csz_ncl")) then load "$HOME/ncl/csz.ncl" end if

; Indicate to host program that header has been loaded
csz_ncl=True
print("dbg: Loading $HOME/ncl/csz.ncl")

undef("clr_map_nm_2_clr_nbr")
undef("clr_map_nm_2_clr_srt")

; Named constants
if (.not.isvar("mss_val_flt")) then mss_val_flt=1.0e36 end if

; Named constants
idx_zero=0
idx_one=1

; Named markers: http://www.ncl.ucar.edu/Document/Graphics/Images/markers.gif
idx_mrk_ast=03 ; Asterisk
idx_mrk_crc_crs=15 ; Circle with cross
idx_mrk_crc_dot=14 ; Circle with dot
idx_mrk_crc_sld=16 ; Circle solid
idx_mrk_crc=04 ; Circle 
idx_mrk_crs=05 ; Cross
idx_mrk_dmn=09 ; Diamond
idx_mrk_dot=01 ; Dot
idx_mrk_pls=02 ; Plus
idx_mrk_sqr=06 ; Square
idx_mrk_str=12 ; Star
idx_mrk_tri_dwn=08 ; Triangle down
idx_mrk_tri_up=07 ; Triangle up

function clr_map_nm_2_clr_nbr(clr_map_nm:string)
; Purpose: Convert colormap name to number of colors
local clr_nbr
begin
if(clr_map_nm.eq."amwg") then clr_nbr=18 end if
if(clr_map_nm.eq."BlueDarkRed18") then clr_nbr=20 end if
if(clr_map_nm.eq."gui_default") then clr_nbr=24 end if
if(clr_map_nm.eq."nrl_sirkes") then clr_nbr=23 end if
if(clr_map_nm.eq."nrl_sirkes_nowhite") then clr_nbr=23 end if
if(clr_map_nm.eq."posneg_2") then clr_nbr=22 end if
if(clr_map_nm.eq."prcp_3") then clr_nbr=25 end if
return(clr_nbr)
end ; end clr_map_nm_2_clr_nbr()

function clr_map_nm_2_clr_srt(clr_map_nm:string)
; Purpose: Convert colormap name to number of colors
local clr_srt
begin
if(clr_map_nm.eq."amwg") then clr_srt=2 end if
if(clr_map_nm.eq."BlueDarkRed18") then clr_srt=2 end if
if(clr_map_nm.eq."gui_default") then clr_srt=2 end if
if(clr_map_nm.eq."nrl_sirkes") then clr_srt=2 end if
if(clr_map_nm.eq."nrl_sirkes_nowhite") then clr_srt=2 end if
if(clr_map_nm.eq."posneg_2") then clr_srt=4 end if
if(clr_map_nm.eq."prcp_3") then clr_srt=4 end if
return(clr_srt)
end ; end clr_map_nm_2_clr_srt()

