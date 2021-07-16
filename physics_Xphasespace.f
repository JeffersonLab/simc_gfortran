      real*8 function peepph(vertex,main)

C     This function is adapted from Pawel Ambrozewicz
C     noted by xucc

      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f
C     For we need these quantities as starting point which may be different
C     from Pawel's noted by xucc

      type(event_main):: main
      type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass.

      real*8 m_pisq, m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 e_photCM,q_photCM  ! virtual photon in CM system
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp
      real*8 E_cm,p_cm          ! proton CM momentum

      real*8 mass               ! recoil mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,tprime,uu       ! Mandelstam variables
      real*8 e_xCM,t_min

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Jttheta            ! Jacobian theta_CM Mx - theta_LAB Ep
      real*8 tcos               ! cos of theta_LAB between Pp and q

c      real*8 mass_vtx
c      common /pawel_crosssection/mass_vtx

*     Variables calculated in transformation to gamma-NUCLEON center of mass.

      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      call transform_to_cm(vertex,main,
     &     gstar,bstar,bstarx,bstary,bstarz,
     &     nustar,qstar,qstarx,qstary,qstarz,
     &     epicm,ppicm,ppicmx,ppicmy,ppicmz,
     &     ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz,
     &     etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz,
     &     thetacm,phicm,phiqn,jacobian,jac_old)
      
C     notice that the above symbol "pi" is borrowed from physics_pion
C     it actually stands for proton. But I am too lazy to change it.
C     noted by xucc

      m_pisq= Mpi02
      m_p   = targ%Mtar_struck
      m_psq = m_p**2

      main%thetacm = thetacm
      main%phicm   = phicm      !gh - verified that this is same as main%phi_pq
      main%pcm     = ppicm

      qsq  = vertex%q2
      invm = sqrt((vertex%nu+m_p)**2-vertex%q**2)
      main%wcm = invm
      Wsq  = invm*invm
      
      mass = sqrt((vertex%nu+targ%Mtar_struck-vertex%p%E)**2
     >           -vertex%Pmiss**2)
      if (mass**2.lt.m_pisq) then ! should never pass this test
         write(6,*) 'warning in peepph: I should never get here!',mass**2
         peepph=0.0
         return
      endif

      E0     = vertex%Ein
      E_prime= vertex%e%E
      nu     = vertex%nu
      epsilon= main%epsilon
      qvec   = vertex%q
c     qvec is the momentum of virtual gamma in lab system.

      Ep   = vertex%p%E
      Pp   = vertex%p%P

      E_cm = (invm**2 + m_p**2 - mass**2)/invm/2.
      p_cm = sqrt(E_cm**2 - m_p**2)

* cos/sin of theta between Pp and q in LAB
cdjg      tcos = cos(main%theta_pq)

C DJG: I changed theta_pq in event.f - need to recalculate tcos
      tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y
     >           +vertex%up%z*vertex%uq%z

* since here the meson is the recoil particle, t for us is not the same thing as 
* main%t, so overwrite it
* for the + sign I spent 1 hour to check. xucc
      tt = 2.0*(m_p**2-m_p*Ep)
c     main%t =tt

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )

      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      q_photCM = sqrt((Wsq-qsq-m_psq)**2 +4.*Wsq*qsq)/invm/2.

      e_xCM    = (Wsq + mass**2 - m_psq)/invm/2.
      if ((e_xCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_xCM*e_photCM-
     *        sqrt((e_xCM**2-mass**2)*(e_photCM**2+qsq)))
      else
         write(6,*)' physics_Xphasespace: no valid t min '
      endif
      tprime = abs(tt)-abs(t_min)
      main%tmin=-t_min

******************************************************************************
*  here ntup.sigcm=d3sigma/dOmega_cm/dMx   [ub/MeV/sr]
******************************************************************************
* Equation 3.27 from Pawel's thesis

      ntup%sigcm = 1./(32.*pi*pi)*mass/q_photCM*p_cm/invm**2

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6
      
*******************************************************************************
* Convert from d3sigma/dOmega_cm/dMx [ub/MeV/sr] -->
* d3sigma/dOmega_lab/dEp [ub/MeV/sr] using 'Jttheta'

* Jacobian for varying Mx.   Ep and theta_pq are independent since they
* are sampled separately. 

      Jttheta=(m_p*qvec*Pp)/(mass*q_photCM*p_cm)

      main%davejac=Jttheta

*******************************************************************************
* Convert to 6-fold d6sigma/dOmega_p/dE_p/dOmega_e/dE_e [ub/MeV^2/sr^2] 
* by multiplying by virtual transverse photon flux factor, gamma_T [1/MeV/sr]

* gh - checked and this agrees with gtpr in physics_pion.f
      gamma_T   = (alpha*E_prime*(Wsq-m_psq))/(4.*pi*pi*E0*m_p*qsq*
     $     (1-epsilon))

*******************************************************************************
* Lab differential cross section
* ub/MeV/sr-->ub/MeV/sr-->ub/MeV^2/sr^2
*******************************************************************************

      peepph=Jttheta*gamma_T*ntup%sigcm

      return 
      end
