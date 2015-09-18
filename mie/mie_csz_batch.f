      program mie_csz_batch

c     Purpose: loop through discrete crystal sizes and compute Mie parameters.
c     Visualize results with IDL function mie_csz_batch.pro

c     Usage: mie_csz_batch >! mie_csz_batch.out

      integer 
     $     idebug,
     $     nang,
     $     num_sizes,
     $     size

      parameter(
     $     num_sizes = 200)

      real 
     $     asymg,
     $     foo_temperature,
     $     ice_imag_index,
     $     ice_real_index,
     $     omega,
     $     qabs,
     $     qext,
     $     qsca,
     $     radius(num_sizes),
     $     radius_step,
     $     start_radius,
     $     wavelength_microns,
     $     x
c
      foo_temperature=230.
      idebug=73
      nang=10
      radius_step=.1
      start_radius=1.
      wavelength_microns=10.0
c     
      if( idebug.eq.73  ) 
     $     write(6,*)num_sizes,' data rows'
c
      call refice(wavelength_microns,foo_temperature,
     $     ice_real_index,ice_imag_index)
c
      if( idebug.eq.73  ) 
     $     write(6,70)wavelength_microns,ice_real_index,ice_imag_index,
     $     ' = wavelength_microns, nreal, nimag'
      if( idebug.eq.73  ) 
     $     write(6,67)
c
      radius(1)=start_radius
      do size=2,num_sizes
         if (radius(size-1).lt.10.) then
            radius(size)=radius(size-1)+radius_step
         else
            radius(size)=radius(size-1)+10.*radius_step
         end if
      end do

      do size=1,num_sizes
         call callbh(nang,wavelength_microns,radius(size),
     +        ice_real_index,ice_imag_index,
     +        x,qsca,qext,asymg,idebug)
         qabs=qext-qsca
         omega=qsca/qext
         if( idebug.eq.73  ) 
     +        write (6,68) wavelength_microns,radius(size),
     $        ice_real_index,ice_imag_index,x,
     $        qsca,qabs,qext,asymg,omega
      end do
c     
 66   format(//,1x,' qsca =',e13.6,3x,' qabs = ',e13.6,3x,
     +     ' qext = ',e13.6,3x,' asymg = ',e13.6,3x,' omega = ',e13.6)
 67   format(
     $     ' Wave   ',2x,' Radius ',2x
     $     'nreal  ',2x,'nimag  ',2x,'  X   ',2x
     $     ' Qsca   ',2x,' Qabs   ',2x,
     +     ' Qext   ',2x,' Asym   ',2x,' Omega ',2x)
 68   format(f6.2,3x,f6.2,3x,
     $     f6.2,3x,e7.2,3x,f6.1,3x,
     $     4(f7.4,3x),f9.6)
 69   format(10(g,3x))
 70   format(3(e12.4,3x),A30)
      end
c-------------------------------------------------------------------
      subroutine callbh(nang,wavel,rad,
     +     refre,refim,
     +     x,qsca,qext,asymg,idebug)
c     ---------------------------------------------------------------
c     callbh calculates the size parameter (x) and relative
c     refractive index (refrel) for a given sphere refractive
c     index, medium refractive index, radius, and free space
c     wavelength.  it then calls bhmie, the subroutine that computes
c     amplitude scattering matrix elements and efficiencies
c     ---------------------------------------------------------------
      parameter (mxa=200,mxa2=400)
      complex refrel,s1(mxa2),s2(mxa2)
      if( idebug.eq.41 ) then
c     this is Bohren's test case p. 482	
         wavel=.6328
         rad=.525
         refre=1.55
         refim=0.
      endif	
      if( idebug.eq.39.or.idebug.eq.41  ) write (6,11)
c     ------------------------------------------------------
c     refmed = (real) refractive index of surrounding medium
c     ------------------------------------------------------
      refmed = 1.0
c     --------------------------------------------
c     refractive index of sphere = refre + i*refim
c     --------------------------------------------
      refrel = cmplx(refre,refim)/refmed
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,12) refmed,refre,refim
c     ----------------------------------------------
c     radius (rad) and wavelength (wavel) same units
c     ----------------------------------------------
      x     = 2. * 3.14159265*rad*refmed/wavel
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,13) rad,wavel
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,14) x
c     ------------------------------------------------
      call bhmie(x,refrel,nang,s1,s2,qext,qsca,qback,asymg,idebug)
      if( idebug.eq.39.or.idebug.eq.41  ) 
     +     write (6,65) qsca,qext,qback,asymg
c     if( idebug.eq.39.or.idebug.eq.41  ) write (6,17)
c     --------------------------------------------------
c     s33 and s34 matrix elemnts normalized by s11.
c     s11 is normalized to 1.0 in the formward direction
c     pol=degree of polarization (incident unpolarized light)
c     --------------------------------------------------
c     s11nor = 0.5*(cabs(s2(1))**2+cabs(s1(1))**2)
c     nan    = 2*nang - 1
c     dang   = 1.570796327/float(nang-1)
c     do 355 j=1,nan
c     aj = j
c     s11=0.5*cabs(s2(j))*cabs(s2(j))
c     s11=s11+0.5*cabs(s1(j))*cabs(s1(j))
c     s12=0.5*cabs(s2(j))*cabs(s2(j))
c     s12=s12-0.5*cabs(s1(j))*cabs(s1(j))
c     pol=-s12/s11
c     s33=real(s2(j)*conjg(s1(j)))
c     s33=s33/s11
c     s34=aimag(s2(j)*conjg(s1(j)))
c     s34=s34/s11
c     s11=s11/s11nor
c     ang = dang*(aj-1.)*57.2958
c     if( idebug.eq.39.or.idebug.eq.41  )
c     +write(6,75) ang,s11,pol,s33,s34
c     355 continue
c     
 65   format(//,1x,' qsca =',e13.6,3x,' qext = ',e13.6,3x,
     +     ' qback = ',e13.6,3x,' asymg = ',e13.6)
 75   format(1x,f6.2,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)
 11   format(/ ' sphere scattering program'//)
 12   format(5x,' refmed = ',f8.4,3x,' refre = ',e14.6,3x,
     +     ' refim = ',e14.6)
 13   format(5x,' sphere radius = ',f7.3,3x,' wavelength = ',f7.4)
 14   format(5x,' size parameter = ',f8.3/)
 17   format(//,2x,'angle',7x,'s11',13x,'pol',13x,'s33',13x,'s34'//)
      end
c     -------------------------------------------------------
c     subroutine bhmie calculates amplitude scattering matrix
c     elements and efficiencies for extinction, total scattering
c     and backscattering for a given size parameter and
c     relative refractive index
c     -------------------------------------------------------
      subroutine bhmie (x,refrel,nang,s1,s2,
     +     qext,qsca,qback,asymg,idebug)
      parameter (mxa=200,mxa2=400)
      dimension amu(mxa),theta(mxa),pi(mxa),tau(mxa),pi0(mxa),pi1(mxa)
      complex d(20000),y,refrel,xi,xi0,xi1,an,bn,s1(mxa2),s2(mxa2)
      complex anm1,bnm1
      double precision psi0,psi1,psi,dn,dx
      dx = x
      y  = x*refrel
c     -----------------------------------
c     series terminated after nstop terms
c     -----------------------------------
      xstop = x + 4.*x**.3333+2.0
      nstop = xstop
      ymod  = cabs(y)
      nmx   = amax1(xstop,ymod)+15
c     
c     print 554,nmx
c     554 format(' ..... nmx = ',i10)
c     
      dang  = 1.570796327/float(nang-1)
      do 555 j=1,nang
         theta(j) = (float(j)-1.)*dang
         amu(j)   = cos(theta(j))
 555  continue
c     ---------------------------------------------------
c     logarithmic derivative d(j) calculated by downward
c     recurrence beginning with initial Value 0.0 + i*0.0
c     at j = nmx
c     ---------------------------------------------------
      d(nmx) = cmplx(0.0,0.0)
      nn     = nmx - 1
      do 120 n=1,nn
         rn = nmx - n + 1
         d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))
 120  continue
      do 666 j=1,nang
         pi0(j) = 0.0
         pi1(j) = 1.0
 666  continue
      nn = 2*nang -1
      do 777 j=1,nn
         s1(j) = cmplx(0.0,0.0)
         s2(j) = cmplx(0.0,0.0)
 777  continue
c     ---------------------------------------------
c     riccati-bessel functions with real argument x
c     calculated by upwar recurrence
c     ---------------------------------------------
      psi0 = dcos(dx)
      psi1 = dsin(dx)
      chi0 = -sin(x)
      chi1 =  cos(x)
      apsi0 = psi0
      apsi1 = psi1
      xi0 = cmplx(apsi0,-chi0)
      xi1 = cmplx(apsi1,-chi1)
      qsca = 0.0
      asymg = 0.0
      anm1=cmplx(0.0,0.0)
      bnm1=cmplx(0.0,0.0)
      n = 1
 200  dn = n
      rn = n
      fn = (2.*rn+1.)/(rn*(rn+1.))
      psi = (2.*dn-1.)*psi1/dx-psi0
      apsi = psi
      chi = (2.*rn-1.)*chi1/x - chi0
      xi = cmplx(apsi,-chi)
c     
      an = (d(n)/refrel+rn/x)*apsi - apsi1
      an = an/((d(n)/refrel+rn/x)*xi - xi1)
c     
      bn = (refrel*d(n)+rn/x)*apsi - apsi1
      bn = bn/((refrel*d(n)+rn/x)*xi - xi1)
c     
      qsca = qsca+(2.*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))
c     
c     This is where the zender hack for computing the asymmetry factor g
c     is performed.  This series is in Bohren & Huffman p. 120.
c     I must store two an,bn terms to compute each term in the g series.
c     On each loop of n the following is computing the n-1 term.
c     There is no n=0 term so it should be zero first time in order to
c     avoid the infinity in my algorithm when rn = 1.
c     
      if (n.gt.1) then
         asymg=asymg+((rn-1.)*(rn+1.)/rn)*
     &        real(anm1*conjg(an)+bnm1*conjg(bn))+
     &        ((2.*rn-1.)/((rn-1.)*rn))*
     &        real(anm1*conjg(bnm1))
c     asymg=asymg+((rn-1.)*(rn+1.)/rn)*
c     &        real(anm1)*real(an)+aimag(anm1)*aimag(an)+
c     &        real(bnm1)*real(bn)+aimag(bnm1)*aimag(bn)
c     asymg = asymg+((2.*rn-1.)/((rn-1.)*rn))*
c     &        real(anm1)*real(bnm1)+aimag(anm1)*aimag(bnm1)
c     anr=real(an)
c     ani=aimag(an)
c     bnr=real(bn)
c     bni=aimag(bn)
c     asymg = asymg + ((rn-1.)*(rn+1.)/rn)*
c     &     (anm1r * anr + anm1i * ani +
c     &     bnm1r * bnr + bnm1i * bni) +  
c     &     ((2.*rn-1.)/((rn-1.)*rn))*
c     &     (anm1r * bnm1r + anm1i * bnm1i)
c     
      endif
c     save previous an,bn
      anm1=an
      bnm1=bn
c     anm1r=real(anm1)
c     anm1i=aimag(anm1)
c     bnm1r=real(bnm1)
c     bnm1i=aimag(bnm1)
c     
      do 789 j=1,nang
         jj = 2*nang - j
         pi(j) = pi1(j)
         tau(j) = rn*amu(j)*pi(j) - (rn+1.)*pi0(j)
         p = (-1.)**(n-1)
         s1(j) = s1(j) + fn*(an*pi(j)+bn*tau(j))
         t = (-1.)**n
         s2(j) = s2(j) + fn*(an*tau(j)+bn*pi(j))
         if( j .eq. jj ) goto 789
         s1(jj) = s1(jj)+fn*(an*pi(j)*p+bn*tau(j)*t)
         s2(jj) = s2(jj)+fn*(an*tau(j)*t+bn*pi(j)*p)
 789  continue
      psi0 = psi1
      psi1 = psi
      apsi1 = psi1
      chi0 = chi1
      chi1 = chi
      xi1 = cmplx(apsi1,-chi1)
      n = n + 1
      rn = n
      do 999 j=1,nang
         pi1(j) = ((2.*rn-1.)/(rn-1.))*amu(j)*pi(j)
         pi1(j) = pi1(j) -rn*pi0(j)/(rn-1.)
         pi0(j) = pi(j)
 999  continue
      if( n-1-nstop) 200,300,300
 300  qsca = (2./(x*x))*qsca
      qext = (4./(x*x))*real(s1(1))
      qback = (4./(x*x))*cabs(s1(2*nang-1))*cabs(s1(2*nang-1))
      asymg = (4./(x*x*qsca))*asymg
      return
      end
c-----------------------------------------------------------------------
