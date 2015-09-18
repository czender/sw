c********1*********2*********3*********4*********5*********6*********7**
      implicit real*8(a-h,o-z)
      real*8 nr,ni,numpar(10),ksca(10),kext(10)
      common r1,r2,nr,ni,alam,pi,a(10),b(10),c(10),d(10),g(10),reff(10),
     &       veff(10),seff(10),numpar(10),ksca(10),kext(10),pizero(10),
     &       cosbar(10),rext(10),rsca(10),rbar(10),rr(96),ww(96),e(10),
     &       r(2000),w(2000),sizdis(10,2000),td(160),ct(80),s2(80),
     &       p1(10,160),p2(10,160),p3(10,160),p4(10,160),p1f(80),
     &       p2f(80),p3f(80),p4f(80),p1b(80),p2b(80),p3b(80),p4b(80),
     &       nsd(10),jx,ngauss,ipart,ngt,numsd
      dimension steps(20),anglus(20)
      open(unit=5,status='unknown',file='mie.inp')
      open(unit=6,status='unknown',file='mie.out')
      pi = 3.1415926535897932d0
      twopi = pi * 2.d0
      read (5,*) nsizep,ng,ipart
      write (6,10) nsizep,ng,ipart
   10 format(' nsizep=',i2,' ng=',i2,' ipart=',i3)
      ngt = ng * ipart
c********1*********2*********3*********4*********5*********6*********7**
c compute scattering angles
c********1*********2*********3*********4*********5*********6*********7**
      anglel = 0.d0
      read (5,*) numang
      write (6,20) numang
   20 format(' the number of intervals used in scattering angle deter
     &mination is numang =',i3)
      k = 1
   30 format(' interval no.',i3,2x,'is from',f7.2,2x,'to',f7.2,2x,
     &'in steps of',f7.3)
      do 50 j = 1,numang
         read (5,*) step,angleu
         steps(j) = step
         anglus(j) = angleu
         write (6,30) j,anglel,angleu,step
         index = (((angleu - anglel) + 10.d-8) / step) + 1
         k = k - 1
         do 40 i = 1,index
            k = k + 1
            td(k) = anglel + (i - 1) * step
   40    continue
         anglel = angleu
   50 continue
      jx = k
      jxf = 2 * jx - 1
      do 60 j = 1,jx
         cs = dcos(pi * td(j) / 180.d0)
         ct(j) = cs
         s2(j) = 1.d0 - cs * cs
         l = jxf - (j - 1)
         td(l) = 180.d0 - td(j)
   60 continue
c********1*********2*********3*********4*********5*********6*********7**
c read size distribution data
c********1*********2*********3*********4*********5*********6*********7**
      read (5,*) numsd
      write (6,*) numsd
   70 format(' numsd=',i4)
      do 80 n = 1,numsd
         read (5,*) nsd(n),a(n),b(n),c(n),d(n),e(n)
         write (6,90) nsd(n),a(n),b(n),c(n),d(n),e(n)
   80 continue
   90 format(' nsd=',i3,'  a=',e15.8,'  b=',e15.8,'  c=',e15.8,
     &' d=',e15.8,'  e=',e15.8)
c********1*********2*********3*********4*********5*********6*********7**
c read optical constants, etc.
c********1*********2*********3*********4*********5*********6*********7**
  100 read (5,*,end = 290) nr,ni,alam,r1,r2
      if(nsizep .eq. 7) alam = twopi
      write (6,110) nr,ni,alam,r1,r2
  110 format(' nr=',f12.5,' ni=',f12.5,' alam=',f10.5,
     &' r1=',f12.5,' r2=',f12.5)
      coefin = alam * alam / pi
c********1*********2*********3*********4*********5*********6*********7**
c obtain integration points and weights
c********1*********2*********3*********4*********5*********6*********7**
      rnge = r2 - r1
      div = rnge / ipart
      do 120 i = 1,ipart
         xbot = (i - 1) * div + r1
         xtop = xbot + div
         lower = (i - 1) * ng
         call gausst(ng,xbot,xtop,rr,ww)
         do 125 kk = 1,ng
            kkk = lower + kk
            r(kkk) = rr(kk)
            w(kkk) = ww(kk)
  125    continue
  120 continue
c********1*********2*********3*********4*********5*********6*********7**
c compute size distribution at integration points
c********1*********2*********3*********4*********5*********6*********7**
      call sizeds
c********1*********2*********3*********4*********5*********6*********7**
c initializations
c********1*********2*********3*********4*********5*********6*********7**
      do 130 n = 1,numsd
         numpar(n) = 0.d0
         ksca(n) = 0.d0
         kext(n) = 0.d0
         cosbar(n) = 0.d0
         rsca(n) = 0.d0
         rext(n) = 0.d0
         rbar(n) = 0.d0
         reff(n) = 0.d0
         veff(n) = 0.d0
         seff(n) = 0.d0
         g(n) = 0.d0
         do 135 i = 1,jxf
            p1(n,i) = 0.d0
            p2(n,i) = 0.d0
            p3(n,i) = 0.d0
            p4(n,i) = 0.d0
  135    continue
  130 continue
      jco = 0
c********1*********2*********3*********4*********5*********6*********7**
c jco loop (new particle) begins here
c********1*********2*********3*********4*********5*********6*********7**
  140 jco = jco + 1
      x = r(jco) * twopi / alam
      call sinpar(x,qext,qsca,cosb)
      write(6,*) "x=",x," r=",r(jco)," qext=",qext," qsca=",qsca,
     +" cosb=",cosb
      do 170 n = 1,numsd
         sw = sizdis(n,jco) * w(jco)
         swpi = sw * pi
         do 150 i = 1,jx
            p1(n,i) = p1f(i) * sw + p1(n,i)
            p2(n,i) = p2f(i) * sw + p2(n,i)
            p3(n,i) = p3f(i) * sw + p3(n,i)
            p4(n,i) = p4f(i) * sw + p4(n,i)
  150    continue
         jxx = jx + 1
         jxf = 2 * jx - 1
         kk = 0
         do 160 i = jxx,jxf
            kk = kk + 1
            k = i - 2 * kk
            p1(n,i) = p1b(k) * sw + p1(n,i)
            p2(n,i) = p2b(k) * sw + p2(n,i)
            p3(n,i) = p3b(k) * sw + p3(n,i)
            p4(n,i) = p4b(k) * sw + p4(n,i)
  160    continue
         z1 = r(jco)
         z2 = z1 * z1
         z3 = z2 * z1
         pir2qs = z2 * qsca * swpi
         ksca(n) = ksca(n) + pir2qs
         kext(n) = kext(n) + z2 * qext * swpi
         rsca(n) = rsca(n) + z1 * pir2qs
         rext(n) = rext(n) + z3 * qext * swpi
         numpar(n) = numpar(n) + sw
         rbar(n) = rbar(n) + z1 * sw
         reff(n) = reff(n) + z3 * swpi
         g(n) = g(n) + z2 * swpi
         cosbar(n) = cosbar(n) + cosb * pir2qs
  170 continue
      if(jco .ne. ngt) go to 140
c********1*********2*********3*********4*********5*********6*********7**
c jco loop (new particle) ends here
c********1*********2*********3*********4*********5*********6*********7**
      do 260 n = 1,numsd
         cosbar(n) = cosbar(n) / ksca(n)
         pizero(n) = ksca(n) / kext(n)
         rsca(n) = rsca(n) / ksca(n)
         rext(n) = rext(n) / kext(n)
         rbar(n) = rbar(n) / numpar(n)
         reff(n) = reff(n) / g(n)
         rf = reff(n)
         z2 = 0.d0
         z3 = 0.d0
         do 180 jco = 1,ngt
            rmrf = r(jco) - rf
            rm2 = rmrf * rmrf * sizdis(n,jco) * w(jco) * r(jco) * r(jco)
            rm3 = rm2 * rmrf
            z2 = z2 + rm2
            z3 = z3 + rm3
  180    continue
         veff(n) = z2 * pi / g(n) / rf / rf
         seff(n) = z3 * pi / g(n) / rf / rf / rf / veff(n) ** 1.5
         do 190 i = 1,jxf
            cs = coefin / ksca(n)
            p1cs = p1(n,i) * cs
            p2cs = p2(n,i) * cs
            p1(n,i) = 0.5d0 * (p1cs + p2cs)
            p2(n,i) = 0.5d0 * (p1cs - p2cs)
            p3(n,i) = p3(n,i) * cs
            p4(n,i) = p4(n,i) * cs
  190    continue
c********1*********2*********3*********4*********5*********6*********7**
c printed output
c********1*********2*********3*********4*********5*********6*********7**
         write (6,200) nsd(n),a(n),b(n),c(n),d(n),e(n)
  200    format(' following for nsd=',i3,' a=',e15.8,' b=',e15.8,
     &   ' c=',e15.8,' d=',e15.8,' e=',e15.8)
         write (6,210)
  210    format(' scat. angle',3x,'p11',12x,'p21',12x,'p33',
     &   12x,'p43',6x,'polarization',/)
         do 220 i = 1,jxf
            xp = -100.d0 * p2(n,i) / p1(n,i)
            write (6,230) td(i),p1(n,i),p2(n,i),p3(n,i),p4(n,i),xp
  220    continue
  230    format(1x,f8.3,4d15.8,d12.5)
         write (6,240) rext(n),rsca(n),rbar(n),reff(n)
  240    format(' rext='e12.5,' rsca=',e12.5,' rbar=',e12.5,' reff=',
     &   e12.5)
         write (6,250) veff(n),seff(n),g(n)
  250    format(' veff=',e12.5,' seff=',e12.5,'   g=',e12.5)
         qsca = ksca(n) / g(n)
         qext = kext(n) / g(n)
         write (6,270) cosbar(n),pizero(n),numpar(n)
         write (6,280) ksca(n),kext(n),qsca,qext
  260 continue
  270 format(' cosbar=',f12.6,' pizero=',f12.6,' numpar=',f12.6)
  280 format(' ksca=',f12.6,' kext=',f12.6,' qsca=',f12.6,
     &' qext=',f12.6)
      go to 100
  290 stop
      end
c********1*********2*********3*********4*********5*********6*********7**
c********1*********2*********3*********4*********5*********6*********7**
      subroutine sinpar(x,qext,qsca,cosb)
      implicit real*8 (a-h,o-z)
      real*8 nr,ni,numpar(10),ksca(10),kext(10)
      complex*16 wn,wn1,wn2,sman,sman1,smbn,smbn1,nc,z
      complex*16 rrf,rrfx,tc1,tc2,cdum1,cdum2
      dimension taun(80),taun1(80),taun2(80)
      dimension pix1(80),pix2(80),pix3(80),z(7000)
      common r1,r2,nr,ni,alam,pi,a(10),b(10),c(10),d(10),g(10),reff(10),
     &       veff(10),seff(10),numpar(10),ksca(10),kext(10),pizero(10),
     &       cosbar(10),rext(10),rsca(10),rbar(10),rr(96),ww(96),e(10),
     &       r(2000),w(2000),sizdis(10,2000),td(160),ct(80),s2(80),
     &       p1(10,160),p2(10,160),p3(10,160),p4(10,160),p1f(80),
     &       p2f(80),p3f(80),p4f(80),p1b(80),p2b(80),p3b(80),p4b(80),
     &       nsd(10),jx,ngauss,ipart,ngt,numsd
      equivalence (wn,wnr),(wn1,wn1r)
   10 nc = dcmplx(nr,-ni)
      rrf = 1.d0 / nc
      rx = 1.d0 / x
      rrfx = rrf * rx
      t1 = dsqrt(x * x * (nr * nr + ni * ni))
      nmx1 = 1.1d0 * t1
      if(nmx1 .le. 6999) go to 30
      write (6,20) nmx1
   20 format(' dimension of z and limit for nmx1 should be increased,',
     & ' nmx1 =',i10)
      stop
   30 nmx2 = t1
      if(nmx1 .gt. 150) go to 40
      nmx1 = 150
      nmx2 = 135
   40 z(nmx1 + 1) = (0.d0,0.d0)
      do 50 n = 1,nmx1
         nn = nmx1 - n + 1
         z(nn) = (nn + 1) * rrfx - 1.d0 / ((nn + 1) * rrfx + z(nn+1))
   50 continue
      do 60 j = 1,jx
         pix1(j) = 0.d0
         pix2(j) = 1.d0
         taun2(j) = 0.d0
         taun1(j) = ct(j)
   60 continue
      t1 = dcos(x)
      t2 = dsin(x)
      wn2 = dcmplx(t1,-t2)
      wn1 = dcmplx(t2,t1)
      wn = rx * wn1 - wn2
      tc1 = z(1) * rrf + rx
      tc2 = z(1) * nc + rx
      sman = (tc1 * wnr - wn1r) / (tc1 * wn - wn1)
      smbn = (tc2 * wnr - wn1r) / (tc2 * wn - wn1)
      smani = dimag(sman)
      smbni = dimag(smbn)
      smanr = dreal(sman)
      smbnr = dreal(smbn)
      if(dabs(smanr) .lt. 1.d-30) smanr = 0.d0
      if(dabs(smani) .lt. 1.d-30) smani = 0.d0
      if(dabs(smbnr) .lt. 1.d-30) smbnr = 0.d0
      if(dabs(smbni) .lt. 1.d-30) smbni = 0.d0
      smbn1 = smbn
      sman1 = sman
      sman1i = dimag(sman1)
      smbn1i = dimag(smbn1)
      sman1r = dreal(sman1)
      smbn1r = dreal(smbn1)
      smanr = 1.5 * smanr
      smani = 1.5 * smani
      smbnr = 1.5 * smbnr
      smbni = 1.5 * smbni
      do 70 j = 1,jx
         p = pix2(j)
         t = taun1(j)
         smanrp = smanr * p
         smanip = smani * p
         smbnrp = smbnr * p
         smbnip = smbni * p
         smanrt = smanr * t
         smanit = smani * t
         smbnrt = smbnr * t
         smbnit = smbni * t
         p1f(j) = smanrp + smbnrt
         p2f(j) = smanip + smbnit
         p3f(j) = smbnrp + smanrt
         p4f(j) = smbnip + smanit
         p1b(j) = smanrp - smbnrt
         p2b(j) = smanip - smbnit
         p3b(j) = smbnrp - smanrt
         p4b(j) = smbnip - smanit
   70 continue
      qext = 2.d0 * (smanr + smbnr)
      qsca = (smanr * smanr + smani * smani + smbnr * smbnr +
     &        smbni * smbni) / 0.75d0
      cosbqs = 0.d0
      n = 2
c********1*********2*********3*********4*********5*********6*********7**
c major loop begins here
c********1*********2*********3*********4*********5*********6*********7**
   80 t1 = 2 * n - 1
      t3 = t1 + 2
      wn2 = wn1
      wn1 = wn
      wn = t1 * rx * wn1 - wn2
      cdum1 = z(n)
      cdum2 = n * rx
      tc1 = cdum1 * rrf + cdum2
      tc2 = cdum1 * nc + cdum2
c      sman = (tc1 * wnr - wn1r) / (tc1 * wn - wn1)
c      smbn = (tc2 * wnr - wn1r) / (tc2 * wn - wn1)
      sman = (tc1 * wnr - wn1r) / (tc1 * wn - wn1)
      smbn = (tc2 * wnr - wn1r) / (tc2 * wn - wn1)
      smani = dimag(sman)
      smbni = dimag(smbn)
      smanr = dreal(sman)
      smbnr = dreal(smbn)
      if(dabs(smanr) .lt. 1.d-30) smanr = 0.d0
      if(dabs(smani) .lt. 1.d-30) smani = 0.d0
      if(dabs(smbnr) .lt. 1.d-30) smbnr = 0.d0
      if(dabs(smbni) .lt. 1.d-30) smbni = 0.d0
      qext = qext + t3 * (smanr + smbnr)
      tx = smanr * smanr + smani * smani + smbnr * smbnr +
     &     smbni * smbni
      qsca = qsca + t3 * tx
      t2 = n - 1
      do 90 j = 1,jx
         t1pix2 = t1 * pix2(j)
         ctj = ct(j)
         pix1j = pix1(j)
         pix3(j) = (t1pix2 * ctj - n * pix1j) / t2
         taun(j) = ctj * (pix3(j) - pix1j) - t1pix2 * s2(j) + taun2(j)
   90 continue
      t5 = n
      t4 = t1 / (t5 * t2)
      t2 = (t2 * (t5 + 1.d0)) / t5
      cosbqs = cosbqs + t2 * (sman1r * smanr + sman1i * smani +
     &         smbn1r * smbnr + smbn1i * smbni) + t4 * (sman1r *
     &         smbn1r + sman1i * smbn1i)
      t2 = n * (n + 1)
      t1 = t3 / t2
      k = (n / 2) * 2
      do 110 j = 1,jx
         p = t1 * pix3(j)
         t = t1 * taun(j)
         smanrp = smanr * p
         smanip = smani * p
         smbnrp = smbnr * p
         smbnip = smbni * p
         smanrt = smanr * t
         smanit = smani * t
         smbnrt = smbnr * t
         smbnit = smbni * t
         p1f(j) = p1f(j) + smanrp + smbnrt
         p2f(j) = p2f(j) + smanip + smbnit
         p3f(j) = p3f(j) + smbnrp + smanrt
         p4f(j) = p4f(j) + smbnip + smanit
         if(k .eq. n) go to 100
         p1b(j) = p1b(j) + smanrp - smbnrt
         p2b(j) = p2b(j) + smanip - smbnit
         p3b(j) = p3b(j) + smbnrp - smanrt
         p4b(j) = p4b(j) + smbnip - smanit
         go to 110
  100    p1b(j) = p1b(j) - smanrp + smbnrt
         p2b(j) = p2b(j) - smanip + smbnit
         p3b(j) = p3b(j) - smbnrp + smanrt
         p4b(j) = p4b(j) - smbnip + smanit
  110 continue
      if(tx .lt. 1.d-14) go to 130
      do 120 j = 1,jx
         pix1(j) = pix2(j)
         pix2(j) = pix3(j)
         taun2(j) = taun1(j)
         taun1(j) = taun(j)
  120 continue
      smbn1 = smbn
      sman1 = sman
      sman1i = dimag(sman1)
      smbn1i = dimag(smbn1)
      sman1r = dreal(sman1)
      smbn1r = dreal(smbn1)
      n = n + 1
      if(n .le. nmx2) go to 80
c********1*********2*********3*********4*********5*********6*********7**
c major loop ends here
c********1*********2*********3*********4*********5*********6*********7**
      write (6,20)
      stop
  130 do 140 j = 1,jx
         t1 = p1f(j)
         t2 = p2f(j)
         t3 = p3f(j)
         t4 = p4f(j)
         p1f(j) = t3 * t3 + t4 * t4
         p2f(j) = t1 * t1 + t2 * t2
         p3f(j) = t1 * t3 + t2 * t4
         p4f(j) = t2 * t3 - t4 * t1
         t1 = p1b(j)
         t2 = p2b(j)
         t3 = p3b(j)
         t4 = p4b(j)
         p1b(j) = t3 * t3 + t4 * t4
         p2b(j) = t1 * t1 + t2 * t2
         p3b(j) = t1 * t3 + t2 * t4
         p4b(j) = t2 * t3 - t4 * t1
  140 continue
      t1 = 2.d0 * rx * rx
      qext = qext * t1
      qsca = qsca * t1
      cosb = 2.d0 * cosbqs * t1 / qsca
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
c********1*********2*********3*********4*********5*********6*********7**
      subroutine sizeds
      implicit real*8(a-h,o-z)
      real*8 nr,ni,numpar(10),ksca(10),kext(10)
      common r1,r2,nr,ni,alam,pi,a(10),b(10),c(10),d(10),g(10),reff(10),
     &       veff(10),seff(10),numpar(10),ksca(10),kext(10),pizero(10),
     &       cosbar(10),rext(10),rsca(10),rbar(10),rr(96),ww(96),e(10),
     &       r(2000),w(2000),sizdis(10,2000),td(160),ct(80),s2(80),
     &       p1(10,160),p2(10,160),p3(10,160),p4(10,160),p1f(80),
     &       p2f(80),p3f(80),p4f(80),p1b(80),p2b(80),p3b(80),p4b(80),
     &       nsd(10),jx,ngauss,ipart,ngt,numsd
      do 190 n = 1,numsd
         nsdn = nsd(n)
         go to (10,50,70,100,120,115,111,112),nsdn
c********1*********2*********3*********4*********5*********6*********7**
c the following (nsd=1) is the two parameter gamma distribution (2.56)
c (standard)
c********1*********2*********3*********4*********5*********6*********7**
   10    rb = 1.d0 / b(n)
         if(rb .gt. 2.d0) go to 30
         write (6,20) b(n)
   20    format(' b for the gamma distribution must be less than',
     &   ' 0.5, but has been read in as',1pd11.4)
         stop
   30    ab = rb / a(n)
         gm = rb - 2.d0
         dlab = gm * dlog(ab)
         dlgm = dlgama(gm)
         do 40 i =1,ngt
            sizdis(n,i) = dexp(dlab + (gm - 1.d0) *
     &      dlog(r(i)) - ab * r(i) - dlgm)
   40    continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
c the following (nsd=2) is the three parameter gamma distribution
c (modified) deirmendjian
c********1*********2*********3*********4*********5*********6*********7**
   50    rac = (a(n) + 1.d0) / c(n)
         dlgn = dlgama(rac)
         dlb = rac * dlog(b(n)) - dlgn
         do 60 i = 1,ngt
            sizdis(n,i) = dexp(dlb + a(n) * dlog(r(i)) -
     &                    b(n) * r(i) ** c(n)) * c(n)
   60    continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
c the following (nsd=3) is the bimodal gamma distribution (2.59,c=a2)
c********1*********2*********3*********4*********5*********6*********7**
   70    rb = 1.d0 / b(n)
         if(rb .gt. 2.d0) go to 80
         write (6,20) b(n)
         stop
   80    a1b = rb / a(n)
         a2b = rb / c(n)
         gm = rb - 2.d0
         dla1b = gm * dlog(a1b)
         dla2b = gm * dlog(a2b)
         dlgm = dlgama(gm)
         do 90 i = 1,ngt
            d1 = dexp(dla1b + (gm - 1.d0) * dlog(r(i)) -
     &           a1b * r(i) - dlgm)
            d2 = dexp(dla2b + (gm - 1.d0) * dlog(r(i)) -
     &           a2b * r(i) - dlgm)
            sizdis(n,i) = 0.5d0 * d1 + 0.5d0 * d2
   90    continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
cthe following (nsd=4) is the log-normal distribution(2.60,a=rg,b=sigma)
c********1*********2*********3*********4*********5*********6*********7**
  100    temp = dsqrt(2.d0 * pi) * b(n)
         a2 = dlog(a(n))
         a3 = 2.d0 * b(n) * b(n)
         do 110 i = 1,ngt
            a1 = dlog(r(i))
            sizdis(n,i) = dexp(-(a1 - a2) ** 2 / a3) / temp / r(i)
  110    continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
cthe following (nsd=7) is the log-normal distribution(2.60,a=rg,b=sigma)
cbut not 2.60 exactly, note log sigma in denominator, etc. zender
c********1*********2*********3*********4*********5*********6*********7**
 111     temp = dsqrt(2.d0 * pi) * dlog(b(n))
         write (6,*) 'you are using a hack option 7: log-normal'
         a2 = dlog(a(n))
         a3 = 2.d0 * dlog(b(n)) * dlog(b(n))
         do 114 i = 1,ngt
            a1 = dlog(r(i))
            sizdis(n,i) = dexp(-(a1 - a2) ** 2 / a3) / temp / r(i)
 114        continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
cthe following (nsd=8) is the delta fn distribution
cthat is all particles are the same size. make sure rmin and rmax 
cenclose the radius bin desired. zender
c********1*********2*********3*********4*********5*********6*********7**
 112     write (6,*) 'you are using a hack option 8: delta function'
         do 113 i = 1,ngt
            sizdis(n,i) = 1./(r2-r1)
 113     continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
cthe following (nsd=6) is the bimodal log-normal distribution
c(2.60,a=rg0,b=sigma0,c=rg1,d=sigma1,e=frac in mode 0. zender
c********1*********2*********3*********4*********5*********6*********7**
 115     temp = dsqrt(2.d0 * pi) * dlog(b(n))
         temp2 = dsqrt(2.d0 * pi) * dlog(d(n))
         write (6,117) 
 117     format(' you are using an untested software hack',
     &   ' the bi-modal log-normal distribution')
         a2 = dlog(a(n))
         a3 = 2.d0 * dlog(b(n)) * dlog(b(n))
         b2 = dlog(c(n))
         b3 = 2.d0 * dlog(d(n)) * dlog(d(n))
         do 116 i = 1,ngt
            a1 = dlog(r(i))
            b1 = dlog(r(i))
            sizdis(n,i) = e(n)*dexp(-(a1 - a2)**2 / a3) /temp/r(i)
     & +(1.d0-e(n))*dexp(-(b1 - b2) ** 2 / b3)/temp2/r(i)
 116     continue
         go to 190
c********1*********2*********3*********4*********5*********6*********7**
c the following (nsd=5) is the power law distribution (2.61,
c n(r)=const*r**(-a) from rmin=b to rmax=c and zero otherwise
c********1*********2*********3*********4*********5*********6*********7**
  120    a1 = 1.d0 - a(n)
         if(b(n) .gt. 0.d0) go to 140
         write (6,130)
  130    format(' rmin=b must be greater than zero for the power',
     &   ' law - execution terminated')
         stop
  140    if(dabs(a1) .gt. 1.d-15) go to 150
         const = dlog(c(n) / b(n))
         go to 160
  150    const = a1 / (c(n) ** a1 - b(n) ** a1)
  160    continue
         do 170 i = 1,ngt
            sizdis(n,i) = 0.d0
            if(r(i) .lt. b(n) .or. r(i) .gt. c(n)) go to 170
            sizdis(n,i) = const * r(i) ** (-a(n))
  170    continue
  180    continue
         do 185 i = 1,ngt
            if(sizdis(n,i) .lt. 1.d-30) sizdis(n,i) = 0.d0
  185    continue
  190 continue
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
c********1*********2*********3*********4*********5*********6*********7**
      subroutine gausst(ng,x1,x2,xp,wt)
c********1*********2*********3*********4*********5*********6*********7**
c ng=total number of gauss points including both sides
c x1=start x2=end xp=points wt=weight
c********1*********2*********3*********4*********5*********6*********7**
      implicit real*8(a-h,o-z)
      dimension xp(1),wt(1)
      data pi,ps,dxl/3.141592653589793d0,1.013211836423378d-01,1.d-16/
      xmid = (x2 + x1) / 2.d0
      xhaf = (x2 - x1) / 2.d0
      dng = ng
      nn = ng / 2
      n2 = nn * 2
      if(n2 .eq. ng) go to 20
      xp(nn + 1) = xmid
      wt(nn + 1) = 1.d0
      if(ng .lt. 2) return
      pn = 1.d0
      n = 0
   10 n = n + 2
      dn = n
      dm = dn - 1.d0
      pn = pn * (dm / dn)
      if(n. lt .n2) go to 10
      wt(nn + 1) = 2.d0 * xhaf / (dng * pn) ** 2
   20 i = 0
      c = pi / dsqrt(dng * (dng + 1.d0) + 0.5d0 - ps) / 105.d0
   30 i = i + 1
      di = i
      z = ps/(4.d0 * di - 1.d0) ** 2
      zz = (105.d0 + z * (210.d0 - z * (2170.d0 - z * (105812.d0
     &     - 12554474.d0 * z))))
      x = dcos(zz * c * (di - 0.25d0))
   40 n = 1
      dm = 1.d0
      pni = 1.d0
      pnj = x
   50 n = n + 1
      dn = n
      pnk = ((dm + dn) * x * pnj - dm * pni) / dn
      pni = pnj
      pnj = pnk
      dm = dn
      if(n .lt. ng) go to 50
      dx = pnj * (1.d0 - x * x) / dng / (pni - x * pnj)
      x = x - dx
      if(dabs(dx) .gt. dxl) go to 40
      j = ng + 1 - i
      xp(i) = xmid - xhaf * x
      xp(j) = xmid + xhaf * x
      wt(i) = 2.d0 * xhaf * (1.d0 - x * x) / (dng * pni) ** 2
      wt(j) = wt(i)
      if(i .lt. nn) go to 30
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      function dlgama(x)
c********1*********2*********3*********4*********5*********6*********7**
c this function calculates the log of the gamma function for use in
c the mie particle size distribution subroutine sizeds
c********1*********2*********3*********4*********5*********6*********7**
      implicit double precision (a-h,o-z)
      parameter (nterms = 26)
      double precision c(nterms)
      data c /   1.d0,       5.772156649015329d-1,-6.558780715202538d-1,
     & -4.20026350340952d-2, 1.665386113822915d-1,-4.21977345555443d-2,
     & -9.621971527877d-3,   7.218943246663d-3,   -1.1651675918591d-3,
     & -2.152416741149d-4,   1.280502823882d-4,   -2.01348547807d-5,
     & -1.2504934821d-6,     1.133027232d-6,      -2.056338417d-7,
     &  6.116095d-9,         5.0020075d-9,        -1.1812746d-9,
     &  1.043427d-10,        7.7823d-12,          -3.6968d-12,
     &  5.1d-13,            -2.06d-14,            -5.4d-15,
     &  1.4d-15,             1.d-16       /
      parameter (xtol = -30.d0)
c********1*********2*********3*********4*********5*********6*********7**
      if(x .le. 0.d0) then
        write(*,'(''0'',''> error : invalid argument in <dlgama>'')') x
        stop
      end if
      if(dabs(dnint(x)-x) .lt. 1.d-3) then
c********1*********2*********3*********4*********5*********6*********7**
c integer values of x
c********1*********2*********3*********4*********5*********6*********7**
        sum = 0.d0
        do 10 y = 2.d0,x - 1.d0,1.d0
          sum = sum + dlog(y)
   10   continue
        dlgama = sum
      else
c********1*********2*********3*********4*********5*********6*********7**
c arbitrary x  -  reduce to interval [0,1] if necessary
c********1*********2*********3*********4*********5*********6*********7**
        if(x .le. 1.d0) then
          z = x
        else
          n = x
          z = x - n
        end if
c********1*********2*********3*********4*********5*********6*********7**
c series expansion over [0,1]
c********1*********2*********3*********4*********5*********6*********7**
        sum = 0.d0
        do 20 k = nterms,1,-1
          if(k * dlog(z) .gt. xtol) then
            sum = sum + c(k) * z ** k
          end if
   20   continue
        dlgama = -dlog(sum)
        if(x .gt. 1.d0) then
c********1*********2*********3*********4*********5*********6*********7**
c recursion formula for x > 1
c********1*********2*********3*********4*********5*********6*********7**
          sum = 0.d0
          do 30 i = 1,n - 1
            sum = sum + dlog(n - i + z)
   30     continue
          dlgama = dlgama + sum + dlog(z)
        end if
      end if
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
