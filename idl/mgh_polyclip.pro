;+
; NAME:
;   MGH_POLYCLIP
;
; PURPOSE:
;   Clip an arbitrary polygon on the X-Y plane to a line parallel
;   to the X or Y axis using the Sutherland-Hodgman algorithm.
;
; CATEGORY:
;   Graphics, Region of Interest, Geometry
;
; CALLING SEQUENCE:
;   result = MGH_POLYCLIP(clip, dir, neg, polin, COUNT=count)
;
; RETURN VALUE
;   The function returns the clipped polygon as a [2,n] array. The
;   second dimension will equal the value of the COUNT argument, except
;   where COUNT is 0 in which case the return value is -1.
;
; POSITIONAL PARAMETERS
;   cval (input, numeric sclar)
;     The value of X or Y at which clipping is to occur
;
;   dir (input, integer scalar)
;     Specifies whether clipping value is an X (dir = 0) or Y (dir =
;     1) value.
;
;   neg (input, integer scalar)
;     Set this argument to 1 to clip to the negtive side, 0 to clip to
;     the positive side. 
;
;   polin (input, floating array)
;     A [2,m] vector defining the polygon to be clipped.
;
; KEYWORD PARAMETERS
;   COUNT (output, integer)
;     The number of vertices in the clipped polygon.
;
; PROCEDURE:
;   The polygon is clipped using the Sutherland-Hodgman algorithm.
;
;   This function is based on JD Smith's implementation of the
;   Sutherland-Hodgman algorithm in his POLYCLIP function. He can
;   take all of the credit and none of the blame.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-10:
;     Written, based on JD Smith's POLYClIP.
;-

function mgh_polyclip, cval, dir, neg, poly, COUNT=count

   compile_opt IDL2

   count = 0

   if size(poly, /N_DIMENSIONS) eq 0 then return, -1

   ;; Vector "in" specifies whether each vertex is inside
   ;; the clipped half-plane

   case dir of
      0B: in = neg ? reform(poly[0,*] lt cval) : reform(poly[0,*] gt cval)
      else: in = neg ? reform(poly[1,*] lt cval) : reform(poly[1,*] gt cval)
   endcase

   ;; Calculate number of points in polygon--it is a little
   ;; more efficient to get it from the size of "in" than
   ;; from the dimensions of "poly"

   np = n_elements(in)

   ;; Vector "inx" specifies whether an intersection with the clipping line
   ;; is made by the segment joining each vertex with the one before.

   inx = in xor shift(in,1)

   ;; Precalculate an array of shifted vertices, used in calculating
   ;; intersection points in the loop.

   pols = shift(poly,0,1)

   ;; Loop thru vertices

   for k=0,np-1 do begin

      ;; If this segment crosses the clipping line, add the intersection
      ;; to the output list. I tried calculating the intersection points
      ;; outside the loop in an array operation but it turned out slower.

      if inx[k] then begin
         s0 = pols[0,k]
         s1 = pols[1,k]
         p0 = poly[0,k]
         p1 = poly[1,k]
         case dir of
            0B: ci = [cval,s1+(p1-s1)/(p0-s0)*(cval-s0)]
            else: ci = [s0+(p0-s0)/(p1-s1)*(cval-s1),cval]
         endcase
         polout = count eq 0 ? [ci] : [[polout],[ci]]
         count = count + 1
      endif

      ;; If this vertex is inside the clipped half-plane add it to the list

      if in[k] then begin
         polout = count eq 0 ? [poly[*,k]] : [[polout],[poly[*,k]]]
         count = count + 1
      endif

   endfor

   return, count gt 0 ? polout : -1

end
