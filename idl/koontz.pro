; From velt@splinter.cgd.ucar.edu  Mon Aug 15 06:44:02 1994
; Received: from sage.cgd.ucar.EDU by ra.cgd.ucar.EDU (4.1/ NCAR Mail Server 04/10/90)
; 	id AA11111; Mon, 15 Aug 94 06:44:02 MDT
; Received: from splinter.rad.usf.edu by sage.cgd.ucar.EDU (8.6.4/ NCAR Mail Server 04/10/90)
; 	id GAA29803; Mon, 15 Aug 1994 06:44:00 -0600
; Received: from grater.rad (grater.rad.usf.edu [131.247.75.15]) by splinter.rad.usf.edu (8.6.5/8.6.5) with SMTP id IAA11378 for <zender@sage.cgd.ucar.edu>; Mon, 15 Aug 1994 08:43:57 -0400
; Received: by grater.rad (5.0/SMI-SVR4)
; 	id AA10233; Mon, 15 Aug 1994 08:43:39 +0500
; Date: Mon, 15 Aug 1994 08:43:39 +0500
; From: velt@splinter.rad.usf.edu (Robert Velthuizen)
; Message-Id: <9408151243.AA10233@grater.rad>
; To: zender@sage.cgd.ucar.EDU
; Subject: Re: high/low labels on contour plots
; Content-Length: 2044
; Status: RO
; 
; The routine below finds "peakness" in a 2D array; it looks for the relative
; height (depth) compared to immediate neighbors. If you're only interested
; in real local minima and maxima, do:
; 
; maxima = where( koontz( image) ge 8)
; minima = where( koontz( image) lt 1)
; 
; I have a bunch of variants of this routine, this one is for 2 dimensional
; images and distributes the  "peakness" over neighbors of equal height. You
; could also give peakness based on a random number. It may also be possible
; to improve on efficiency.
; 
; Hope this helps,
; Robert Velthuizen,
; Digital Medical Imaging Program of the
; H. Lee Moffitt Cancer Center and Research Institute at the
; University of South Florida.
; 
; I think bytscl is by far the simplest way to make "koontz" work with
; real data. The algorithm was designed to work on histogram data, which
; necessarily is positive integer. However, you could just do a type
; change in koontz to floats, and use min(fspace)-1 as instead of -1
; for the border, and select the seach space as everything > min(fspace).
; 
; If the neighborhood of the minima and
; maxima is very flat, then you may want to do something nonlinear to
; enhance the gradients around the candidate peaks. You could do that
; in two runs of koonz: first select the candidate peaks using bytscl,
; use the result to create a new image but with the neighborhood of
; each candidate extremum scaled between 0-255; the area's that are left
; out of this you set to 128 or so.
; 
; Good luck,
; Robert Velthuizen.
; ===============================================================================
 
Function Koontz, fspace
;+
; Name:		Koontz
;
; Purpose:	find the peaks in an image
;		written to find peaks in 2D histogram of feature space
;
; Reference:	A Khotanzad and A Bouarfa
;		Image segmentation by a parallel, non-parametric histogram
;		based clustering algorithm
;		Pattern Recognition, 23(9), pp 961-973, 1990
;-
d=size(fspace)
;
; The border can be [-1,min(fspace)-1] in which case the 
; search space should be [where(Im ge 0), where(Im ge min(fspace))]
;
border=-1
;border=min(fspace)-1
;
im=intarr(d(1)+1,d(2)+1)+border         ; +border, so boundary will never be a top
;im=fltarr(d(1)+1,d(2)+1)+border           ; +border, so boundary will never be a top
im(1:d(1),1:d(2))=fspace             ; Im holds the image in a border of +border's
Nm=n_elements(Im)
                                     ; N is the 3x3 neighborhood
N=reform([[shift(Im,-1,-1)],[shift(Im,-1, 0)],[shift(Im,-1, 1)],$
          [shift(Im, 0,-1)],[shift(Im, 0, 0)],[shift(Im, 0, 1)],$
          [shift(Im, 1,-1)],[shift(Im, 1, 0)],[shift(Im, 1, 1)]],Nm,9)
                                     ; Offsets for the neighborhood pixels
Offset=[1,1,1,0,0,0,-1,-1,-1] + [1,0,-1,1,0,-1,1,0,-1]*(d(1)+1)
Counts=0.0*Im
ss=where(Im ge 0)                    ; ss stands for searchspace
;ss=where(Im ge min(fspace))           ; ss stands for searchspace
Nm=n_elements(ss)
for i=0,Nm-1 do begin
   top=where(N(ss(i),*) eq max(N(ss(i),*)),nt)  ; Find all highest neighbors
   address=ss(i)+offset(top)                    ; pointer
   Counts(address) = Counts(address)+1.0/nt     ; distribute "peakness"
endfor
return,counts(1:d(1),1:d(2))
end

