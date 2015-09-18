      program rdhtrn
C-----------------------------------------------------------------------
C     
C     Read hitran-like data file and select only desired data
C     
C     Bruce P. Briegleb
C     19 September 1994
C     
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
C     
C     local variables
c     
      integer 
     $     namdat,              ! unit for input data
     $     namout,              ! unit for output data
     $     mol,                 ! molecule number
     $     iso                  ! isotope number
c
      real    
     $     wavlin,              ! wavenumber of line
     $     linint,              ! line intensity (cm2/molecule)cm-1  at 296K
     $     linhwd,              ! line halfwidth in cm-1/atm at 296K
     $     lwsten,              ! lower state energy in cm-1
     $     linexp               ! temperature exponent for air broadened halfwidth
C     
C-----------------------------------------------------------------------
C     
C     Read in line and select:
C     
      namdat = 1
      namout = 2
 10   call rdline(namdat, mol   ,    iso,  wavlin, linint, 
     $     linhwd, lwsten, linexp) 
C     
      if( mol .eq. 1 ) then
         if( iso .eq. 1 ) then
            write(namout,15) mol,iso,wavlin,linint,linhwd,lwsten,linexp
 15         format(i2,i1,f12.6,1x,1pe10.3,1x,0pf5.4,1x,f10.4,1x,f4.2)
         endif
      endif
      if( wavlin .lt. 3018. ) goto 10
C     
C     Done
C     
      end
      subroutine rdline(namdat, mol   ,    iso,  wavlin, linint, 
     $     linhwd, lwsten, linexp) 
c-----------------------------------------------------------------------
c     
c     reads lbl file and returns line data.
c     
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c     
      integer 
     $     namdat,              ! integer unit of the lbl data
     $     mol,                 ! molecule number
     $     iso                  ! isotope number
c
      real    
     $     wavlin,              ! wavenumber of line
     $     linint,              ! line intensity (cm2/molecule)cm-1  at 296K
     $     linhwd,              ! line halfwidth in cm-1/atm at 296K
     $     lwsten,              ! lower state energy in cm-1
     $     linexp               ! temperature exponent for air broadened halfwidth
c     
c     output arguments
C     
      logical prnt
      data prnt / .true. /
C     
C-----------------------------------------------------------------------
C     
      read(namdat,10) mol,iso,wavlin,linint,linhwd,lwsten,linexp
c
c full AFGL spec:
c 10   format(i2,i1,f12.6,1p2e10.3,0p2f5.4,f10.4,f4.2,f8.6,2i3,2a9,3i1,3i2)
c
c BPB's original:
c   10      format(i2,i1,f12.6,1x,1pe10.3,1x,0pf5.4,1x,f10.4,1x,f4.2)
c
 10   format(i2,i1,f12.6,1pe10.3,10x,0pf5.4,5x,f10.4,f4.2)
C     
      if( prnt ) then
         write(6,20) mol, iso, wavlin
 20      format(' line read in for mol = ',i2,' iso = ',i2,
     $        ' wavenumber = ',f8.3,' cm-1')
      endif
C     
      return
      end
