      program calc_semi_xsec

C     program to calculate the SIDIS model cross section at one set of kinematics
c     uses the same subroutine as used in SIMC so there cannot be any mis-matches
c     between the model used in simc and the model used to get the "central" vertex cross section
      
      USE structureModule
      implicit none
      include 'simulate.inc'

      real*8 ebeamin,q2in, xin, zhadin, thpqin, Ain, Zin
      real*8 peepiX, survivalprob,charge
      real*8 sigsemi
      real*8 dphi
      integer i,nphibin

      type(event_main):: main
      type(event):: vertex, vertex0
      
      write(6,*) 'Enter beam energy (GeV), Q2, xbj, z, thetapq, A, Z,',
     >            ' pion charge, # of phi bins'
      read(5,*) ebeamin,q2in,xin,zhadin,thpqin,Ain,Zin,charge,nphibin

      if(charge.gt.0.0) then
         doing_hplus=.true.
      else
         doing_hplus=.false.
      endif

      dphi=2.*pi/nphibin
c     start calculating all ther vertex quantities needed
c     target
      targ%Mtar_struck=Mp
      targ%Z=Zin
      targ%A=Ain
c electron stuff      
      vertex%e%phi = 3.0*pi/2.0 ! assume electrons in HMS and no out-of-plane angle
      vertex%Ein = ebeamin*1000.0d0
      vertex%Q2=q2in*1.d6             !convert to MeV
      vertex%nu = vertex%Q2/2.0/Mp/xin
      vertex%e%E=vertex%Ein-vertex%nu
      vertex%e%theta=2.0*asin(sqrt(vertex%Q2/4./vertex%Ein/vertex%e%E))
      vertex%ue%x = sin(vertex%e%theta)*cos(vertex%e%phi)
      vertex%ue%y = sin(vertex%e%theta)*sin(vertex%e%phi)
      vertex%ue%z = cos(vertex%e%theta)
      vertex%q = sqrt(vertex%Q2 + vertex%nu**2)
      vertex%xbj = xin

c q-vector      
      vertex%uq%x = - vertex%e%P*vertex%ue%x / vertex%q
      vertex%uq%y = - vertex%e%P*vertex%ue%y / vertex%q
      vertex%uq%z =(vertex%Ein - vertex%e%P*vertex%ue%z)/ vertex%q

      vertex%zhad=zhadin
      vertex%p%P = sqrt((vertex%zhad*vertex%nu)**2-Mpi2)
      vertex%theta_pq=thpqin/degrad
      vertex%pt2=vertex%p%P**2*(1.0-cos(vertex%theta_pq)**2)
      
      doing_semipi=.true.
C     loop over phi - for now, this will give you the same cross section for all
C     phi, until we add interference terms to the cross section model
      write(6,*)'#Eb,x,Q2,z,th_pq: ', ebeamin,q2in,xin,zhadin,thpqin
      write(6,*) '#A, Z',Ain, Zin
      write(6,*) '#phi, sighad, sigsemi (ub/GeV2/sr2)'
      do i=1,nphibin
         main%phi_pq = (i-1)*dphi + dphi/2.0         
         sigsemi = peepiX(vertex,vertex0,main,survivalprob,.false.)
         write(6,*) main%phi_pq,ntup%sigcm,sigsemi*1E6
      enddo

      end
