;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CVS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: cld.pro,v $
; $Source: /home/zender/cvs/idl/cld.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: All the IDL figures for the shape and size sensitivity paper.
;
; $Log: not supported by cvs2svn $
; Revision 1.3  1999-12-31 00:18:15  zender
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin preload of all procedures and functions which the main program
; and event handler depend on. Order of loading is important.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@~/.idle
@utilities.pro
@ibp_clr.pro
@cld_fgr.pro
@cld_wdg.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End preload of all procedures and functions besides the main program
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

