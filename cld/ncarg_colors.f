      SUBROUTINE DFCLRS(numcol,invert,icolor)
C
C Define a set of RGB color triples for colors 1 through 15.
C
c invert=1 colors run from gray, then blue..red sorted by hue.
c invert=0  colors run from gray, then red..blue
c
c icolor=0  COOL
c icolor=1  HSV 
c icolor=2  gray scale
c icolor=3  HLS  
c
c numcol=number of colors to use, starting with index 2 in the table.
c (index 0 and 1 are black and white)
c

        CALL GSCR (1,0,1.,1.,1.)
        CALL GSCR (1,1,0.,0.,0.)

C
c colors 3..15  cylcle thru numcol-1 different hues.  
        n=numcol-1

        if (icolor.eq.1) then
c          HSV, but the first color, index=2, leave as gray, set above
           CALL GSCR (1,2,.8,.8,.8)
           cmin=-20
           cmax=250
           DO I=1,n
              i2=i-1
              hue1= i2 / real(n-1)
              hue1=hue1 + .05*sin(4*3.1415*hue1)
              hue= cmin + hue1*(cmax - cmin) 
              call hsvrgb(hue,1.,1.,R,G,B)
              i3=i+2
              if (invert.eq.1) i3=(n+2) - (i-1)
              CALL GSCR (1,i3,R,G,B)
           enddo
        else if (icolor.eq.0) then
c          COOL colormap
           do i=0,n
              r = (i/real(n))
              g = 1-r
              b = 1
              i3=i+2
              if (invert.eq.1) i3 = (n+2) -i 
              CALL GSCR (1,i3,R,G,B)
           enddo
        else if (icolor.eq.2) then
c          gray scale
           do i=0,n
              r = (i/real(n))**(.4)
              g = r
              b = r
              i3=i+2
              if (invert.eq.1) i3 = (n+2) - i
              CALL GSCR (1,i3,R,G,B)
           enddo
        else if (icolor.eq.3) then
c          hls, from scd
           CALL GSCR (1,2,.8,.8,.8)
           cmin=36
           cmax=360+36
           DO I=1,n
              hue1= i / real(n)
              hue= cmin + hue1*(cmax - cmin) 
              call hlsrgb(hue,60.,75.,R,G,B)
              i3=i+2
              if (invert.eq.1) i3=(n+2) - (i-1)
              CALL GSCR (1,i3,R,G,B)
           enddo
        endif

C
C Done.
C
        RETURN
C
      ENd






