c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.1.1.1  1998-09-15 02:06:43  zender
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
c Revision 1.1  1993/04/28  00:02:50  zender
c Initial revision
c
c
      subroutine spline(x,y,n1,x0,y0,n2,q)
c
c     Input Variables : x,y....vectors of size (0:n1) and containing
c                              the real data points
c                       x0.....vector of size (n2) containing 
c                              x locations where interpolation is req'd
c                       q......flag set to 1 normally, and set to 2
c                              if debugging output is desired   
c                       
c     Output Variables : y0....vector containing interpolated y values
c                              corresponding to the x0 locations

      integer n1,n2,q
c      double precision x(0:n1),y(0:n1),x0(n2),y0(n2)
c      double precision delta(0:99),gpp(0:99)
c      double precision a(99),b(99),c(99),d(99),k(99)
c      double precision f1,f2,f3,f4,f5,f6

      real x(0:n1),y(0:n1),x0(n2),y0(n2)
      real delta(0:99),gpp(0:99)
      real a(99),b(99),c(99),d(99),k(99)
      real f1,f2,f3,f4,f5,f6

c   ** set up delta array
      do 10 i=0,n1-1
          delta(i)=x(i+1)-x(i)
  10  continue

c   ** solve tridiagonal matrix for gpp
      do 20 i=1,n1-1
          b(i) = (delta(i-1)+delta(i))/3.0d0
  20  continue
      do 30 i=1,n1-2
           a(i+1) = delta(i)/6.0d0
           c(i) = delta(i)/6.0d0
  30  continue

      do 40 i=1,n1-1
          k(i)=((y(i+1)-y(i))/delta(i))-((y(i)-y(i-1))/delta(i-1))
  40  continue
      if (q.eq.2) then
        write(*,*) 'a=',(a(m),m=2,n1-1)
        write(*,*) 'b=',(b(m),m=1,n1-1)
        write(*,*) 'c=',(c(m),m=1,n1-2)
         write(*,*) 'k=',(k(m),m=1,n1-1)
      endif

      call trdiag(n1-1,a,b,c,d,k)

      if (q.eq.2) then
         write(*,*) 'd=',(d(m),m=1,n1-1)
      endif
      
      do i=1,n1-1
         gpp(i) = d(i)
      end do

c    ** set end conditions
      gpp(0)=0
      gpp(n1)=0

      do j=1,n2

c    ** find the correct interval
      i = 0

c    ** if series is an increasing series do this
      if (x(0).lt.x(n1)) then
  50     if ((x(i).le.x0(j)).and.(i.ne.n1)) then
            i = i + 1
            goto 50
         endif
      else
 51      if ((x(i).ge.x0(j)).and.(i.ne.n1)) then
            i = i + 1
            goto 51
         endif
      endif
c      write(*,*) i
      i = i - 1
      i = max(0,i)

c    ** now evaluate y0 at location x0 by cubic spline
      f1 = (gpp(i)/6.0d0)*((x(i+1)-x0(j))**3/delta(i))
      f2 = (gpp(i)/6.0d0)*(-delta(i)*(x(i+1)-x0(j)))
      f3 = (gpp(i+1)/6.0d0)*((x0(j)-x(i))**3/delta(i))
      f4 = (gpp(i+1)/6.0d0)*(-delta(i)*(x0(j)-x(i)))
      f5 = y(i)*(x(i+1)-x0(j))/delta(i)
      f6 = y(i+1)*(x0(j)-x(i))/delta(i)
c      if (q.eq.2) then
c         write(*,*) 'Inside Spline'
c         write(*,*) 'x0 = ',x0
c         write(*,*) 'gpps = ',i,gpp(i),gpp(i+1)
c         write(*,*) 'effs = ',f1,f2,f3,f4,f5,f6
c      end if

      y0(j) = f1+f2+f3+f4+f5+f6
      end do
 
      return
      end


      SUBROUTINE TRDIAG(N,A,B,C,X,G)
c      double precision a(n),b(n),c(n)
c      double precision x(n),g(n),bb(1000)          
      real a(n),b(n),c(n)
      real x(n),g(n),bb(1000)          
C.....THIS SUBROUTINE SOLVES TRIDIAGONAL SYSTEMS OF EQUATIONS          
C.....BY GAUSS ELIMINATION     
C.....THE PROBLEM SOLVED IS MX=G WHERE M=TRI(A,B,C)          
C.....THIS ROUTINE DOES NOT DESTROY THE ORIGINAL MATRIX      
C.....AND MAY BE CALLED A NUMBER OF TIMES WITHOUT REDEFINING 
C.....THE MATRIX     
C.....N = NUMBER OF EQUATIONS SOLVED (UP TO 1000)  
C.....FORWARD ELIMINATION      
C.....BB IS A SCRATCH ARRAY NEEDED TO AVOID DESTROYING B ARRAY         
      DO 1 I=1,N     
      BB(I) = B(I)   
    1 CONTINUE       
      DO 2 I=2,N     
      T = A(I)/BB(I-1)         
      BB(I) = BB(I) - C(I-1)*T 
      G(I) = G(I) - G(I-1)*T   
    2 CONTINUE   
C.....BACK SUBSTITUTION        
      X(N) = G(N)/BB(N)        
      DO 3 I=1,N-1   
      J = N-I        
      X(J) = (G(J)-C(J)*X(J+1))/BB(J)    
    3 CONTINUE       
      RETURN         
      END       

