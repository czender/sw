;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: cld_cmn.com,v $
; $Source: /home/zender/cvs/idl/cld_cmn.com,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.com
;
; $Log: not supported by cvs2svn $
; Revision 1.5  2000-01-10 23:36:28  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/10 19:33:33  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:48  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 20:12:45  zender
; *** empty log message ***
;
; Revision 1.1  1999/12/31 00:18:15  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:04  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:46  zender
; Imported sources
;
; Revision 2.0  1994/03/08  00:54:07  zender
;  the wonderful widgeting cloud data viewer
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
@ibp_clr.com
;
common cld_sct,fgr_lst, $ 	;structure of graphing functions
                        devices 	;structure of device driver functions

common cld_fls,	primary_fl, $	        ;for functions needing one file
                        secondary_fl, $	;for functions needing two files
                        tertiary_fl, $	;for functions needing three files
                        quaternary_fl	        ;for functions needing four files

common cld_wdg,     $
                        Wchroma, $
                        Wdest, $
                        Wformat, $
                        Wfl_out, $
                        Wshape

common cld_output,	$
                        chroma_idx, $
                        chroma_lst, $
                        dest_idx, $
                        dest_lst, $
                        fmt_idx, $
                        fmt_lst, $
                        chroma_nbr, $
                        dest_nbr, $
                        fmt_nbr, $
                        pll_nbr, $
                        shape_nbr, $
                        fl_out, $    ;fl_nm
                        pll_idx, $
                        pll_lst, $
                        shape_idx, $
                        shape_lst


