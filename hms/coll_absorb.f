      subroutine coll_absorb(ppi,thick,transmission)
c Subroutine to calculate probability that pion will be absorbed in the SOS
c collimator - can be used for HMS assuming same geometry applies. TH - checked that.

C DJG For now, I use the "reaction" cross sections. Also, only implemented for
C DJG pi+ - no pi-.

! TH - cross section from Dave's parameterization of the energy and A dependence up to
! TH   pion momenta of ~ 2 GeV, we extended that to > 2 GeV during Fpi analysis 

      implicit none

      real*8 Ppi,Epi,Tpi,mpi                    !pion kinematics
      real*8 Navagadro
      parameter(Navagadro = 6.0221367d+23)
      parameter(mpi = 139.56995)

      real*8 mate_dens
      real*8 mate_t
      real*8 mate_A
      real*8 T(14),sigreac(14),qreac(14)

      real*8 sighi,siglo,qhi,qlo,transmission,lambdai,thick
      real*8 Tlo,Thi,sigA,sigAlo,sigAhi,sum

      integer i


! TH - DG original + TH extension for > 2GeV
      data T /85.0,125.0,165.0,205.0,245.0,315.0,584.02,711.95,870.12,
     >     1227.57,1446.58,1865.29,2858.0,4159.0/

CCCCCCCCCC  TOTAL CROSS SECTIONS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      data sigreac /170.12,189.49,178.61,158.11,116.9,85.86,55.5,55.5,
!     >     55.5,55.5,55.5,55.5,55.5,55.5/
!      data qreac /0.59980,0.5753,0.58046,0.59842,0.65218,0.69421,.779
!     >     ,.779,.779,.779,.779,.779,0.779,0.779/
CCCCCCCCCCCCC   REACTION CROSS SECTIONS CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TH - DG Original + TH extension for > 2GeV
      data sigreac /26.03,84.47,117.3,117.4,101.9,69.58,42.5,44.7,47.9,
     >     46.5,45.2,39.6,35.34,33.15/
      data qreac /0.948,0.659,0.56,0.5342,0.5452,0.60796,0.699,0.689,
     >     0.679,0.683,0.688,0.704,0.7483,0.7705/
CCCCCCCCCCCCC   ABSORPTION(?) CROSS SECTIONS  CCCCCCCCCCCCCCCCCCCCCCCCCCC
c      data sigreac /8.81,25.57,29.57,30.245,15.63,8.302,4.98,4.98,4.98,
c     >     4.98,4.98,4.98/
c      data qreac /0.986,0.763,0.744,0.679,0.762,0.6506,0.75,0.75,
c     >     0.75,0.75,0.75,0.75/

c      data sigreac_minus /25.0,100.0,120.,110.0,100.0,80.0/
c      data qreac_minus  /0.97,0.60,0.52,0.50,0.53,0.60/


      mate_dens = 17.0                !Collimator density (g/cm^3)
      mate_t = 6.35                   !Total aluminum seen by pion(cm)
      mate_A = 171.57

      Epi = sqrt(Ppi**2+mpi**2)
      Tpi = Epi - mpi

      sigA = 0.0
      do i = 1,14
         if((Tpi.gt.T(i)).and.(Tpi.le.T(i+1))) then
            Thi = T(i+1)
            Tlo = T(i)
            sighi = sigreac(i+1)
            siglo = sigreac(i)
            qhi = qreac(i+1)
            qlo = qreac(i)

            sigAhi = sighi*mate_A**qhi
            sigAlo = siglo*mate_A**qlo
            sigA = (sigAlo*(Thi-Tpi) + sigAhi*(Tpi-Tlo))/(Thi-Tlo)
            sigA = sigA*1.d-27  !convert to cm^2
         endif
      enddo

      sum = 0.0
      
      lambdai = mate_dens*Navagadro*sigA/mate_A
      sum = sum + thick*lambdai
      transmission = exp(-sum)

      end



