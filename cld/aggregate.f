c   Here is the code. I put some comments in and I believe you should
c be able to use it without any problem. look for ways of improving
c efficiency. Also, print out M_k/(N_k * x_k) if you want to be sure
c that it is functioning prperly. The value should lie between 1 and 2
c and is typically less than 1.5. Negative values should never
c exist if the routine is fed a positive distribution.
c   Remember to match your units of the collection coefficients with
c those of the routine.
c   You are welcome to use this routine, and as I said before, all I ask
c is that you acknowledge the source in your work. The only reason
c I mention this again is that people I gave the code to in the past
c did not always do so and I think it is the fair thing to do.
c         Good luck - Graham
c p.s. I forgot the value of ap in the code I sent you. Here it is:
c       p=2.
c       ap=0.5+(p+1.0)*(p+1.0)/(8.0*p)
c ***********************************************************
c             COLLECTION
c ***********************************************************
      subroutine sxy(x,r0,amkd,ankd,amk,ank,t,akbar,ap,lpa,lpx,lpa1
     +  ,psi,scm,f,f2,psi2,am3,az0,am4,cn,cm)
c * *
c * subroutine sxy calculates the evolution of the spectrum for
c * a mass weighted kernel (x+y) using the method of moments
c * Solves Eqs. 9a,b in #1 (Tzivion, Feingold and Levin, 1987, JAS).
c * Refer to #1 for details.
c * Definitions of parameters:
c *    Input: x  - the mass array (x_k) [g]
c *           r0 - the air density [g/cc]
c *           amkd - bin mass supplied by the dynamics [g/cc]
c *           ankd - bin number supplied by the dynamics [#/cc]
c *           akbar - collection coefficients K_i,j 
c *           lpx - number of bins
c *           lpa =lpx (not necessary here to use lpa AND lpx)
c *           lpa1 - lpx + 1 (i.e., the number of bin boundaries)
c *           t - time step
c *    Output: amk - bin mass after collection [g/cc]
c *            ank - bin number after collection [#/cc]
c * * * * * *
      real x(lpa1)
      real ankd(lpx),amkd(lpx)
      real ank(lpx),amk(lpx)
      real akbar(lpa,lpa)
      real az0(lpa),am3(lpa),am4(lpa)
      real cm(lpa),cn(lpa),scm(lpa)
      real psi(lpa),f(lpa),psi2(lpa),f2(lpa)
c * dt = t = time step
      dt=t
      lk=lpx
      p=2.
      ap=0.5+(p+1.0)*(p+1.0)/(8.0*p)
c * cm and cn are work vectors and used locally only.
c * You don't really need to do this, unless you want to change units
c * or something of that nature.
      do 9 l=1,lk
      cm(l)=amkd(l)
      cn(l)=ankd(l)
      if(cm(l).gt.0..and.cn(l).gt.0)then
c * average bin mass
      scm(l)=cm(l)/cn(l)
      else
      scm(l)=0.
      endif

9     continue

      g=9./8.
c * *
c Proceed to calculate the new values of M_k and N_k in each bin
c * *
      do 1 k=1,lk
c * *
c az0 is the 2'nd moment (in terms of mass) or 6'th (in terms of
c radius).
c am3 is the third moment, am4 is the 4'th moment (in terms of mass)
c These appear in equations 9a,b #1, when the (x+y) weighting is used. 
c * *
      az0(k)=ap*scm(k)*cm(k)
      am3(k)=(ap**3)*(scm(k)**2)*cm(k)
      am4(k)=(ap**6)*(scm(k)**3)*cm(k)
c * *
c psi_k and f_k in Eq. 14, #1
c * *
c     psi(k)=cm(k)/(x(k)**2)*(ap*scm(k)/x(k)-1.)
      psi(k)=cn(k)/(2.*x(k))*((scm(k)/x(k))**2+scm(k)/x(k)-2.)
c     f(k)=2.*cm(k)/(x(k)**2)*(2.-ap*scm(k)/x(k))
      f(k)=cn(k)/x(k)*(2.+scm(k)/x(k)-(scm(k)/x(k))**2)
c psi_k and f_k in Eq. 13, #1
      psi2(k)=2./x(k)*(cm(k)/x(k)-cn(k))
      f2(k)=2./x(k)*(2.*cn(k)-cm(k)/x(k))
c *
c summation terms
c * 
      sm1=0.
      sm2=0.
      sm3=0.
      sm4=0.
      sm5=0.
      sn1=0.
      sn2=0.
      sn3=0.
      sn4=0.
      do 2 i=k,lk
      sm5=sm5+akbar(i,k)*(az0(k)*cn(i)+cm(k)*cm(i))
      sn4=sn4+akbar(i,k)*(cn(k)*cm(i)+cm(k)*cn(i))
2     continue
100   continue
      if(k.eq.1)then
      cm(k)=cm(k)-sm5*dt
      cn(k)=cn(k)-sn4*dt
      go to 56
      else
      endif
      sm3=akbar(k-1,k-1)*(az0(k-1)*cn(k-1)+cm(k-1)**2)
101   continue
      do 6 i=1,k-1
c * if the average mass is less than the lower bin boundary then 
c * don't remove mass from this bin.
      if(scm(k).ge.x(k))then
      sm2=sm2+akbar(k,i)*(4.*x(k)**2*psi2(k)*cm(i)+
     1x(k)/2.*(4.*psi2(k)+f2(k))*az0(i)-(2.*psi2(k)-f2(k))*am3(i)
     2-1./(2.*x(k))*(psi2(k)-f2(k))*am4(i)+psi2(k)*am3(i))
      sn3=sn3+akbar(k,i)*(2.*x(k)*psi2(k)*cm(i)+0.5*f2(k)*az0(i)
     1-1./(2.*x(k))*(psi2(k)-f2(k))*am3(i))
      else
      endif
6     continue
102   continue
      do 621 i=1,k-1
      sm4=sm4+akbar(k,i)*(cn(k)*az0(i)+cm(k)*cm(i))
621   continue
103   continue
      sn2=akbar(k-1,k-1)*cn(k-1)*cm(k-1)
104   continue
      if(k.eq.2)then
      cm(k)=cm(k)+(sm4+sm3-sm5-sm2)*dt
      cn(k)=cn(k)+(sn2-sn3-sn4)*dt
      go to 56
       else
      endif
      do 3 i=1,k-2
c * if the average mass is less than the lower bin boundary then 
c * don't remove mass from this bin.
      if(scm(k-1).ge.x(k-1))then
      sm1=sm1+akbar(k-1,i)*(4.*x(k-1)**2*psi2(k-1)*cm(i)+
     1x(k-1)/2.*(4.*psi2(k-1)+f2(k-1))*az0(i)-
     $(2.*psi2(k-1)-f2(k-1))*am3(i)+psi2(k-1)*am3(i)
     2-1./(2.*x(k-1))*(psi2(k-1)-f2(k-1))*am4(i))
      sn1=sn1+akbar(k-1,i)*(2.*x(k-1)*psi2(k-1)*cm(i)+0.5*f2(k-1)
     1*az0(i)-1./(2.*x(k-1))*(psi2(k-1)-f2(k-1))*am3(i))
      else
      endif
3     continue
105   continue
c * *
c * Update M_k and N_k
c * *
      cm(k)=cm(k)+(sm1-sm2+sm3+sm4-sm5)*dt
      cn(k)=cn(k)+(sn1+sn2-sn3-sn4)*dt
56    continue
1     continue

      do 20 l=1,lk
      ank(l)=cn(l)
      amk(l)=cm(l)
20    continue

      return
      end
