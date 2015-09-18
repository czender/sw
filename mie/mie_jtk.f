      program miedis
c--------------------------------------------------------------c
c callbh calculates the size parameter (x) and relative
c refractive index (refrel) for a given sphere refractive
c index, medium refractive index, radius and free space
c wavelength. It the calls bhmie, the subroutine that
c computes amplitude scaterring matrix elements and efficiences
c--------------------------------------------------------------c
      dimension wave(11),refr(11),refi(11)
      complex refrel,s1(200),s2(200)
      data wave/0.223,0.255,0.270,0.280,0.290,0.300,0.328,
     $          0.525,0.950,1.80,3.20/
      data refr/8*1.530,1.520,1.340,1.220/
      data refi/10*8.00E-3,1.00E-2/
c      data refr/11*1.5/
c      data refi/11*0.001/
      write(6,11)
c--------------------------------------------------------------c
c refmed = (real) refractive index of surrounding medium       c
c--------------------------------------------------------------c
      refmed = 1.0
      write(6,12) refmed,refre,refim
c--------------------------------------------------------------c
c radius (rad) and wavelength (wavel) same units               c
c--------------------------------------------------------------c
      pie = 3.14159265
      radn = 0.4
      sigma = 2.2
      en0 = 1.0
      rhop = 2.5
      Wmass = 3.0e-3
      do 5 lam = 1,11
        wavel = wave(lam)
c--------------------------------------------------------------c
c refractive index of sphere = refre + i * refim               c
c--------------------------------------------------------------c
        refre = refr(lam)
        refim = refi(lam)
        refrel = cmplx(refre,refim)/refmed
        drad = 0.001
        rad = 0.00
        pmass = 0.0
        qexn = 0.0
        qscn = 0.0
        gnn = 0.0
        do 10 nrad = 1,400
          if(rad.ge.0.1) drad = 0.01
          if(rad.ge.1.0) drad = 0.02
          rad = rad + drad
          x = 2. * pie * rad * refmed / wavel
c--------------------------------------------------------------c
c nang = number of angles between 0 and 90 degrees             c
c matrix elements calculated at 2 * nang - 1 angles            c
c including 0, 90 and 180 degrees                              c
c--------------------------------------------------------------c
          nang = 11
          dang = 1.570796327/float(nang-1)
         call bhmie(x, refrel, nang, s1, s2, qext, qsca, qback,
     $   g,omega0)
c         print *,qext,qsca,g,omega0
          pr = 0.5 * (alog(rad/radn) / alog(sigma))**2.0
          qexn = qexn + qext * rad * exp(-pr) * drad
          qscn = qscn + qsca * rad * exp(-pr) * drad
          gnn = gnn + g * qsca * rad * exp(-pr) * drad
          pmass = pmass + rad**2.0 * exp(-pr) * drad
          if(rad.ge.5.0) goto 300
   10   continue
  300   facn = en0 / (sqrt(2.0*pie)*alog(sigma))
        soil = pmass
        spext = (3./4.) * (1./rhop) * (qexn/soil)
        omegan = qscn / qexn
        qexn = pie * facn * qexn
        gn = gnn / qscn
        f = gn * gn
        write(6,65) wavel,spext,omegan,gn,f
    5 continue
      write(6,60) rad
   60 format(5x,"Maximum radius = ",1e13.6//)
   65 format(1x,1e13.6,3x,1e13.6,3x,1e13.6,3x,1e13.6,3x,1e13.6)
   11 format(/"Mie Scattering Program"//)
   12 format(5x,"refmed = ",1f8.4,3x,"refre = ",1e14.6,3x,
     $          "refin = ",1e14.6)
      stop
      end

c--------------------------------------------------------------c
c subroutine bhmie calculates amplitude scattering matrix      c
c elements and efficiencies for extinction, total scattering   c
c and backscattering for a given size parameter and            c
c relative index of refraction                                 c
c--------------------------------------------------------------c
      subroutine bhmie(x, refrel, nang, s1, s2, qext,qsca,qback,
     $g,omega0)
      dimension amu(100), theta(100), pi(100), tau(100),
     $          pi0(100),pi1(100)
      complex d(3000),y,refrel,xi,xi0,xi1,an,bn,s1(200),s2(200),
     $        anm1,bnm1
      double precision psi0,psi1,psi,dn,dx
      dx = x
      y = x * refrel
c--------------------------------------------------------------c
c series terminated after NSTOP terms                          c
c--------------------------------------------------------------c
      xstop = x + 4. * x**0.3333 + 2.0
      NSTOP = xstop
      ymod = cabs(y)
      nmx = amax1(xstop,ymod) + 15
      dang = 1.570796327/float(nang - 1)
      do 555 j = 1,nang
         theta(j) = (float(j) - 1.) * dang
         amu(j) = cos(theta(j))
  555 continue
c-------------------------------------------------------------c
c logarithmic derivative d(j) calculated by downward          c
c recurrence beginning with intial value 0.0 + i *0.0         c
c at j = nmx                                                  c
c-------------------------------------------------------------c
      d(nmx) = cmplx(0.0,0.0)
      nn = nmx - 1
      do 120 n =1,nn
         rn = nmx - n + 1
         d(nmx-n) = (rn/y) - (1./(d(nmx-n+1)+rn/y))
  120 continue
      do 666 j = 1,nang
         pi0(j) = 0.0
         pi1(j) = 1.0
  666 continue
      nn = 2 * nang - 1
      do 777 j = 1,nn
         s1(j) = cmplx(0.0,0.0)
         s2(j) = cmplx(0.0,0.0)
  777 continue
c------------------------------------------------------------c
c riccati-bessel functions with real arguement x             c
c calculated by upward recurrence                            c
c------------------------------------------------------------c
      psi0 = dcos(dx)
      psi1 = dsin(dx)
      chi0 = -sin(x)
      chi1 = cos(x)
      apsi0 = psi0
      apsi1 = psi1
      xi0 = cmplx(apsi0,-chi0)
      xi1 = cmplx(apsi1,-chi1)
      qsca = 0.0
      n = 1
  200 dn = n
      rn = n
      fn = (2. * rn + 1.) / (rn * ( rn + 1.))
      psi = (2.*dn-1.) * psi1 / dx - psi0
      apsi = psi
      chi = (2. * rn - 1.) * chi1 / x - chi0
      xi = cmplx(apsi,-chi)
      an = (d(n)/refrel + rn/x) * apsi - apsi1
      an = an/((d(n)/refrel + rn/x) * xi - xi1)
      bn = (refrel*d(n) + rn/x) * apsi - apsi1
      bn = bn/((refrel*d(n) + rn/x) * xi - xi1)
      qsca = qsca+(2.*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))
      if (n.eq.1) then
         g = 0.0
      else
         g = g + gn * real(anm1*conjg(an)+bnm1*conjg(bn))
     $         + fnm1 * real(anm1*conjg(bnm1))
      endif
      do 789 j = 1,nang
         jj = 2 * nang - j
         pi(j) = pi1(j)
         tau(j) = rn * amu(j) * pi(j) - (rn+1.)*pi0(j)
         p = (-1.)**(n-1)
         s1(j) = s1(j) + fn*(an*pi(j) + bn*tau(j))
         t = (-1.)**n
         s2(j) = s2(j) + fn*(an*tau(j) + bn*pi(j))
         if(j.eq.jj) goto 789
         s1(jj) = s1(jj) + fn*(an*pi(j)*p + bn*tau(j)*t)
         s2(jj) = s2(jj) + fn*(an*tau(j)*t + bn*pi(j)*p)
  789 continue
      psi0 = psi1
      psi1 = psi
      apsi1 = psi1
      chi0 = chi1
      chi1 = chi
      xi1 = cmplx(apsi1,-chi1)
      anm1 = an
      bnm1 = bn
      fnm1 = fn
      gn = rn * (rn + 2.) / (rn + 1.)
      n = n + 1
      rn = n
      do 999 j = 1,nang
         pi1(j) = ((2.*rn-1.)/(rn-1.)) * amu(j) * pi(j)
         pi1(j) = pi1(j) - rn * pi0(j) / (rn-1.)
         pi0(j) = pi(j)
  999 continue
      if (n - 1 - NSTOP) 200,300,300
  300 qsca = (2./(x*x)) * qsca
      qext = (4./(x*x)) * real(s1(1))
      qback = (4./(x*x))*cabs(s1(2*nang-1))*cabs(s1(2*nang-1))
      g = (4./(x*x*qsca))*g
      omega0 = qsca / qext
      return
      end
