      subroutine optice(nwav,wav,rindx,iindx)
c
c............ returns index of refraction for ice:
c
      real ninir(3,45),niir(3,10)
      real wav(1),rindx(1),iindx(1)
c
      data ninir/ 
     $ .950,1.302,8.30e-7,1.000,1.302,1.99e-6,1.030,1.301,2.80e-6,
     $1.050,1.301,2.50e-6,1.100,1.300,1.84e-6,1.150,1.299,2.93e-6,
     $1.200,1.298,9.18e-6,1.250,1.297,1.41e-5,1.300,1.296,1.24e-5,
     $1.400,1.295,2.04e-5,1.450,1.294,1.18e-4,1.500,1.294,5.58e-4,
     $1.520,1.294,6.41e-4,1.550,1.294,5.83e-4,1.600,1.293,3.80e-4,
     $1.650,1.293,3.03e-4,1.700,1.2925,1.92e-4,1.75,1.292,1.42e-4,
     $1.800,1.292,1.13e-4,1.850,1.292,6.33e-5,1.900,1.292,3.91e-4,
     $1.950,1.291,1.11e-3,2.000,1.291,1.61e-3,2.050,1.289,1.40e-3,
     $2.100,1.288,8.35e-4,2.150,1.286,6.76e-4,2.200,1.282,3.08e-4,
     $2.250,1.278,.000213,3.000,1.130,.2273,3.05,1.192,.3178,
     $3.075,1.225,.3428  ,3.100,1.280,.3252,3.15,1.547,.2502,
     $3.200,1.557,.1562  ,3.250,1.550,.0900,3.30,1.530,.0625,
     $3.350,1.515,.0440  ,3.400,1.490,.0307,3.45,1.445,.0226,
     $3.500,1.422,.0163  ,3.550,1.408,.0129,3.60,1.395,.0105,
     $3.800,1.356,.0082  ,3.900,1.340,.0104,4.00,1.327,.0124/
c
      data niir/
     $ 8.0,1.219,.0369, 8.5,1.217,.0352, 9.0,1.210,.0365,
     $ 9.5,1.192,.0310,10.0,1.152,.0413,10.5,1.195,.0602,
     $11.0,1.290,.0954,11.5,1.393,.0114,12.0,1.480,.0120,
     $12.5,1.565,.0119/
c
c
      do 10 n=1,45
        wav(n)   = ninir(1,n)
        rindx(n) = ninir(2,n)
        iindx(n) = ninir(3,n)
   10 continue
c
      do 20 n=1,10
        m        = 45 + n
        wav(m)   = niir(1,n)
        rindx(m) = niir(2,n)
        iindx(m) = niir(3,n)
   20 continue
c
      nwav = m
c
      write(6,100) nwav
  100 format(' .... optwav   number of wavelengths = ',i5)
c
      return
      end