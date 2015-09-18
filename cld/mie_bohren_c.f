c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.1.1.1  1998-09-15 02:06:40  zender
c Imported sources
c
c Revision 5.4  1995/06/19  04:54:58  zender
c have made many changes. got rid of fp_in,out,err. added flat updraft
c switch. changed -S -s to both refer to ice. got rid of while loop that
c became infinite when there was no cloud. changed vertical netCDF coord.
c to altitude. fixed bug in appending date and rev. to cmdline.
c
c Revision 5.3  1995/05/24  23:11:17  zender
c updated and validated the cloud model for ansi compatability,
c removed berkeley calls, streamlined architecture dependent
c fortran --> c calling methodology. validated results against
c comps II runs.
c
c Revision 5.2  1995/05/23  00:43:16  zender
c synchronization checkin prior to temperature sensitivity study
c for SPCP stable cloud mods.
c
c Revision 5.1  1994/02/14  21:28:54  zender
c Added DMC94 homogeneous ammonium sulphate solution nucleation,
c added -D 79 SWSCF/SWCF diagnostic. About to change default shape
c to hexagonal ice crystals, but will save bullet routines.
c NetCDF routines are not keeping up with the changes....
c
c Revision 5.0  1993/08/30  01:55:35  zender
c this version was used to make the paper (preprint at least).
c about to start working on an animation.
c
c Revision 4.12  1993/07/21  23:26:28  zender
c added Liou_IR_fudge() as a default routine whenever LIOU is true.
c
c Revision 4.11  1993/06/19  22:05:25  zender
c added ebert&curry parameterization to LW computations, and
c made this the default for the region outside of the cloud grid.
c
c Revision 4.10  1993/06/11  04:35:31  zender
c Played w/ Makefiles, made all CRAY2 refs. into plain CRAY,
c and put all shavano stuff on peace. changed surface type to
c land for predictable albedos, changed tropics SAZ to 30 deg.,
c removed gaussians from initial distributions. changed plev = 113.
c slingo option doesn't work for some reason (-D 53)
c
c Revision 4.9  1993/06/09  00:59:14  zender
c prepared the -Z initial_distribution_type switch to handle everyone's
c different cloud. Made Knollenberg fit his observations, prepared
c optical properties for 50 sizes and removed all num_size == 40 specific
c stuff. changed tau_CCN to 1. s.
c
c Revision 4.8  1993/05/29  00:01:13  zender
c fixed the Liou_fudge() bug, fixed a wind_speed bug, added
c Dowling and Radke's crystal distribution function as default.
c
c Revision 4.7  1993/05/27  16:08:24  zender
c A synchronization check-in for all cloud programs.
c
c Revision 4.4  1993/04/28  02:03:24  zender
c enlarged d() array to 25000 to handle larger size params
c required by hex column dimensions.
c
c Revision 1.1  1993/04/28  00:02:50  zender
c Initial revision
c
c
c      program mie_bohren_c
      subroutine mie_bohren_c
c      nang=10
c      refre=1.55
c      refim=0.
c      rad=.525
c      wavel=.6328
      idebug=39
      nang=10
      refre=1.311
      refim=3.110e-9
      rad=11.69799
      wavel=.55
c
      call callbh(nang,wavel,rad,
     +refre,refim,
     +x,qsca,qext,asymg,idebug)
c
      pi_r_squared=3.141592654*rad*rad
      qabs=qext-qsca
      x_sec_ext=qext*pi_r_squared
      x_sec_sca=qsca*pi_r_squared
      x_sec_abs=qabs*pi_r_squared
      coalbedo=1.-qsca/qext

      print *,'qext = ',qext
      print *,'qsca = ',qsca
      print *,'qabs = ',qabs
      print *,'x_sec_ext = ',x_sec_ext
      print *,'x_sec_sca = ',x_sec_sca
      print *,'x_sec_abs = ',x_sec_abs
      print *,'asymg = ',asymg
      print *,'coalbedo = ',coalbedo
      end
c-------------------------------------------------------------------
      subroutine callbh(nang,wavel,rad,
     +                  refre,refim,
     +                  x,qsca,qext,asymg,idebug)
c       ---------------------------------------------------------------
c       callbh calculates the size parameter (x) and relative
c       refractive index (refrel) for a given sphere refractive
c       index, medium refractive index, radius, and free space
c       wavelength.  it then calls bhmie, the subroutine that computes
c       amplitude scattering matrix elements and efficiencies
c       ---------------------------------------------------------------
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
c       ------------------------------------------------------
c       refmed = (real) refractive index of surrounding medium
c       ------------------------------------------------------
      refmed = 1.0
c       --------------------------------------------
c       refractive index of sphere = refre + i*refim
c       --------------------------------------------
      refrel = cmplx(refre,refim)/refmed
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,12) refmed,refre,refim
c       ----------------------------------------------
c       radius (rad) and wavelength (wavel) same units
c       ----------------------------------------------
      x     = 2. * 3.14159265*rad*refmed/wavel
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,13) rad,wavel
      if( idebug.eq.39.or.idebug.eq.41  ) write(6,14) x
c       ------------------------------------------------
      call bhmie(x,refrel,nang,s1,s2,qext,qsca,qback,asymg,idebug)
      if( idebug.eq.39.or.idebug.eq.41  ) 
     +write (6,65) qsca,qext,qback,asymg
      if( idebug.eq.39.or.idebug.eq.41  ) write (6,17)
c       --------------------------------------------------
c       s33 and s34 matrix elemnts normalized by s11.
c       s11 is normalized to 1.0 in the formward direction
c       pol=degree of polarization (incident unpolarized light)
c       --------------------------------------------------
c      s11nor = 0.5*(cabs(s2(1))**2+cabs(s1(1))**2)
c      nan    = 2*nang - 1
c      dang   = 1.570796327/float(nang-1)
c      do 355 j=1,nan
c        aj = j
c        s11=0.5*cabs(s2(j))*cabs(s2(j))
c        s11=s11+0.5*cabs(s1(j))*cabs(s1(j))
c        s12=0.5*cabs(s2(j))*cabs(s2(j))
c        s12=s12-0.5*cabs(s1(j))*cabs(s1(j))
c        pol=-s12/s11
c        s33=real(s2(j)*conjg(s1(j)))
c        s33=s33/s11
c        s34=aimag(s2(j)*conjg(s1(j)))
c        s34=s34/s11
c        s11=s11/s11nor
c        ang = dang*(aj-1.)*57.2958
c        if( idebug.eq.39.or.idebug.eq.41  )
c     +write(6,75) ang,s11,pol,s33,s34
c  355 continue
c
   65 format(//,1x,' qsca =',e13.6,3x,' qext = ',e13.6,3x,
     + ' qback = ',e13.6,3x,' asymg = ',e13.6)
   75 format(1x,f6.2,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)
   11 format(/ ' sphere scattering program'//)
   12 format(5x,' refmed = ',f8.4,3x,' refre = ',e14.6,3x,
     + ' refim = ',e14.6)
   13 format(5x,' sphere radius = ',f7.3,3x,' wavelength = ',f7.4)
   14 format(5x,' size parameter = ',f8.3/)
   17 format(//,2x,'angle',7x,'s11',13x,'pol',13x,'s33',13x,'s34'//)
      end
c       -------------------------------------------------------
c       subroutine bhmie calculates amplitude scattering matrix
c       elements and efficiencies for extinction, total scattering
c       and backscattering for a given size parameter and
c       relative refractive index. a larger d() array is required
c       for extremely large size parameters.
c       -------------------------------------------------------
      subroutine bhmie (x,refrel,nang,s1,s2,
     +                  qext,qsca,qback,asymg,idebug)
      parameter (mxa=200,mxa2=400)
      dimension amu(mxa),theta(mxa),pi(mxa),tau(mxa),pi0(mxa),pi1(mxa)
      complex d(25000),y,refrel,xi,xi0,xi1,an,bn,s1(mxa2),s2(mxa2)
      complex anm1,bnm1
      double precision psi0,psi1,psi,dn,dx
      dx = x
      y  = x*refrel
c       -----------------------------------
c       series terminated after nstop terms
c       -----------------------------------
      xstop = x + 4.*x**.3333+2.0
      nstop = xstop
      ymod  = cabs(y)
      nmx   = amax1(xstop,ymod)+15
c
c     print 554,nmx
c 554 format(' ..... nmx = ',i10)
c
      dang  = 1.570796327/float(nang-1)
      do 555 j=1,nang
        theta(j) = (float(j)-1.)*dang
        amu(j)   = cos(theta(j))
  555 continue
c       ---------------------------------------------------
c       logarithmic derivative d(j) calculated by downward
c       recurrence beginning with initial Value 0.0 + i*0.0
c       at j = nmx
c       ---------------------------------------------------
      d(nmx) = cmplx(0.0,0.0)
      nn     = nmx - 1
      do 120 n=1,nn
        rn = nmx - n + 1
        d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))
  120 continue
      do 666 j=1,nang
        pi0(j) = 0.0
        pi1(j) = 1.0
  666 continue
      nn = 2*nang -1
      do 777 j=1,nn
        s1(j) = cmplx(0.0,0.0)
        s2(j) = cmplx(0.0,0.0)
  777 continue
c       ---------------------------------------------
c       riccati-bessel functions with real argument x
c       calculated by upwar recurrence
c       ---------------------------------------------
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
c         asymg=asymg+((rn-1.)*(rn+1.)/rn)*
c     &        real(anm1)*real(an)+aimag(anm1)*aimag(an)+
c     &        real(bnm1)*real(bn)+aimag(bnm1)*aimag(bn)
c         asymg = asymg+((2.*rn-1.)/((rn-1.)*rn))*
c     &        real(anm1)*real(bnm1)+aimag(anm1)*aimag(bnm1)
c      anr=real(an)
c      ani=aimag(an)
c      bnr=real(bn)
c      bni=aimag(bn)
c      asymg = asymg + ((rn-1.)*(rn+1.)/rn)*
c     &     (anm1r * anr + anm1i * ani +
c     &     bnm1r * bnr + bnm1i * bni) +  
c     &     ((2.*rn-1.)/((rn-1.)*rn))*
c     &     (anm1r * bnm1r + anm1i * bnm1i)
c
      endif
c     save previous an,bn
      anm1=an
      bnm1=bn
c      anm1r=real(anm1)
c      anm1i=aimag(anm1)
c      bnm1r=real(bnm1)
c      bnm1i=aimag(bnm1)
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
  789 continue
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
  999 continue
      if( n-1-nstop) 200,300,300
  300 qsca = (2./(x*x))*qsca
      qext = (4./(x*x))*real(s1(1))
      qback = (4./(x*x))*cabs(s1(2*nang-1))*cabs(s1(2*nang-1))
      asymg = (4./(x*x*qsca))*asymg
      return
      end
c-----------------------------------------------------------------------
