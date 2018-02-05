      real*8 function peep_eta(vertex,main)

C     first version written by Wenliang Li, 2017

      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 m_eta/547.86/
      real*8 gamma_eta/0.0013/
      real*8 m_pisq, m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp
      real*8 E_cm,p_cm          ! proton CM momentum

      real*8 mass               ! omega mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,t_min,tprime,uu,t_max ! Mandelstam variables
      real*8 e_photCM,e_etaCM

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Breit_wigner
      real*8 Jttheta            ! Jacobian t - theta_LAB
      real*8 Jttheta_fx         ! old Jacobian for fixed Mx
      real*8 tcos               ! cos of theta_LAB between Pp and q

! Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      real*8 sig_eta_test
      real*8 sig1,sig2

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
      main%phicm   = phicm      !gh -  verified that this is same as main%phi_pq
      main%pcm     = ppicm
      main%davejac = jacobian
      main%johnjac = jac_old    !approx. assuming collinear boost.

      qsq  = vertex%q2
      invm = sqrt((vertex%nu+m_p)**2-vertex%q**2)
      main%wcm = invm
      Wsq  = invm*invm

      mass = targ%Mrec_struck
      if (mass**2.lt.m_pisq) then
         write(6,*)mass,sqrt(m_pisq)
         peep_eta=0.0
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
C DJG      tcos = cos(main%theta_pq)

C DJG: I changed theta_pq in event.f - need to recalculate tcos
      tcos = vertex%up%x*vertex%uq%x + 
     >     vertex%up%y*vertex%uq%y+vertex%up%z*vertex%uq%z

* since here the meson is the recoil particle, t for us is not the same thing as 
* main%t, so overwrite it
* for the + sign I spent 1 hour to check. xucc
      tt = 2.0*(m_p**2-m_p*Ep)
c      main%t =tt

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )
      
      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_etaCM  = (Wsq + mass**2 - m_psq)/invm/2.
      if ((e_etaCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_etaCM*e_photCM-
     *        sqrt((e_etaCM**2-mass**2)*(e_photCM**2+qsq)))
         t_max = -qsq + mass**2 -2.*(e_etaCM*e_photCM
     *        +sqrt((e_etaCM**2-mass**2)*(e_photCM**2+qsq)))

      else
         write(6,*)' physics_eta: no valid t min '
      endif
      
      tprime = abs(tt-t_min)
      main%tmin=-t_min

      if (abs(tt).lt.abs(t_min)) then
         write(6,*)' unphysical -t<-t_min ',tt,t_min,tprime
      endif
      if (abs(tt).gt.abs(t_max)) then
         write(6,*)' unphysical -t>-t_max ',tt,t_max,tprime
      endif
      
******************************************************************************
*  we keep the tradition that ntup.sigcm is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
******************************************************************************
**
      sig1 = sig_eta_test(mass/1.e3,qsq/1.e6,invm/1.e3,tt/1.e6,uu/1.e6,
     1       epsilon)

      ntup%sigcm=sig1

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6

C DJG  Breit-Wigner in the event generation - not needed here anymore.      
*******************************************************************************
* Multiply by Breit_wigner [1/MeV] to give proper Mx weighted cross section
* d3sigma/dt/dphi_cm/dMx [ub/MeV^3/rad]
*
* gh - relativistic Breit-Wigner factor (eqn 38.52 of 2004 PDG book)
c      Breit_wigner=(m_eta*gamma_eta)**2/
c     *     ((mass**2-m_eta**2)**2 +(m_eta*gamma_eta)**2)

* gh - since sigma_omega is integrated over the omega peak, normalize
* Breit-Wigner to unit integral
c      Breit_wigner=Breit_wigner/(gamma_eta*pi/2.)

*******************************************************************************
* Convert from d3sigma/dt/dphi_cm/dMx [ub/MeV^3/rad] 
* --> d3sigma/dEp/dOmega_lab [ub/MeV/sr] using 'Jttheta'
*
* Jacobian for fixed Mx.  Ep is a function of theta_pq.
* J=dt/dcos_LAB [MeV^2].
      Jttheta_fx = 2.*m_p*qvec*Pp / ( m_p+nu-qvec*Ep/Pp*tcos )

* Jacobian for varying Mx.   Ep and theta_pq are independent since they
* are sampled separately. 
* J=d(t,Mx)/d(Ep_LAB,cos_LAB) [MeV^2].
      Jttheta = 2.*m_p*qvec*Pp/mass

      main%davejac=Jttheta
      main%johnjac=Jttheta_fx

*******************************************************************************
* Convert to 6-fold d6sigma/dOmega_p/dE_p/dOmega_e/dE_e [ub/MeV^2/sr^2]
* by multiplying by virtual transverse photon flux factor, gamma_T [1/MeV/sr]
*
* gh - checked and this agrees with gtpr in physics_pion.f
      gamma_T   = (alpha*E_prime*(Wsq-m_psq))/(4.*pi*pi*E0*m_p*qsq*
     $     (1-epsilon))

*******************************************************************************
* Lab differential cross section
* ub/MeV^2/rad-->ub/MeV^3/rad-->ub/MeV/sr-->ub/MeV^2/sr^2
*******************************************************************************

CDJG BW in event generation - also want Jacobian at fixed Mx ...

      peep_eta=Jttheta_fx*gamma_T*ntup%sigcm

      return 
      end




C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_eta_test(mass,q2_gev,w_gev,tt,uu,epsilon)

* all calculations are in GeV

      implicit none
      include 'constants.inc'

      real*8 mass,q2_gev,w_gev,tt,uu,epsilon
      real*8 A,B,C,sig_eta
      real*8 m_p,m_psq,wfactor

      
      a = 0.0044
      b = 5
      c = 0.000011

      sig_eta = a * exp( -b * abs(uu)) + c;

      wfactor=(2.21**2-mp**2)**2/(w_gev**2-mp**2)**2
      sig_eta = sig_eta*wfactor/1.e+06

c     GH: check for weird behavior on upper side of peak
      if (mass.gt.1. .and. abs(tt).gt.q2_gev .and. sig_eta.gt.1.e-7)
     *    then
         write(6,*)' tt=',tt,' sig=',sig_eta
         sig_eta=sig_eta*1.e-6
      endif

      sig_eta_test = sig_eta

      return
      end
