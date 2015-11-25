      real*8 function peep_phi(vertex,main)

C     This function is adapted from Pawel Ambrozewicz
c     gh - 14.12.18

      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f
C     For we need these quantities as starting point which may be different
C     from Pawel's noted by xucc

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 sigma_phi
      real*8 m_phi/1019.460/
      real*8 gamma_phi/4.26/
      real*8 m_Ksq, m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp

      real*8 mass               ! phi mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,uu              ! Mandelstam variables
      real*8 tprime
      real*8 e_photCM,e_phiCM,t_min

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Breit_wigner
      real*8 Jttheta            ! Jacobian t-Mx - Ep,theta_p in LAB
      real*8 Jttheta_fx         ! old Jacobian for fixed Mx
      real*8 tcos               ! cos of theta_LAB between Pp and q

! Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      real*8 sig_gmh,sig_diffract,sig_phigmh
      real*8 sig1,sig2,sig3

      if (debug(2)) write(6,*)' peep_phi: enter '

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

      m_Ksq= 493.677**2
      m_p   = targ%Mtar_struck
      m_psq = m_p**2

      main%thetacm = thetacm
      main%phicm   = phicm      !gh - verified that this is same as main%phi_pq
      main%pcm     = ppicm
*     main%davejac = jacobian
*     main%johnjac = jac_old    !approx. assuming collinear boost.

      qsq  = vertex%q2
      invm = sqrt((vertex%nu+m_p)**2-vertex%q**2)
      main%wcm = invm
      Wsq  = invm*invm

      mass = targ%Mrec_struck
      if (mass**2.lt.m_Ksq) then
         peep_phi=0.0
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

* cos/sin of theta between Pp and q in LAB
CDJG      tcos = cos(main%theta_pq)

C DJG: I changed theta_pq in event.f - need to recalculate tcos
      tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y+vertex%up%z*vertex%uq%z
      
* since here the meson is the recoil particle, t for us is not the same thing as 
* main%t, so overwrite it
* for the + sign I spent 1 hour to check. xucc
      tt = 2.0*(m_p**2-m_p*Ep)
c     main%t = tt

      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_phiCM   = (Wsq + mass**2 - m_psq)/invm/2.
      if ((e_phiCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_phiCM*e_photCM-
     *        sqrt((e_phiCM**2-mass**2)*(e_photCM**2+qsq)))
      endif
      main%tmin=t_min
      tprime = abs(tt)-abs(t_min)
      if (tprime.lt.0.) then
         write(6,*)' unphysical -t<-t_min ',tt,t_min
      endif

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )

******************************************************************************
*  we keep the tradition that ntup.sigcm is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
******************************************************************************
*
* Diffractive model used by P. Ambrosewicz
      sig1 = sig_diffract(mass/1.e3,invm/1.e3,qsq/1.e6,
     *     tt/1.e6,t_min/1.e6)

* PYTHIA model with modifications from the HERMES MC
      sig2 = sig_phigmh(qsq,tt,t_min,nu,invm/1.e3,epsilon)

* Parameterization based on saturated Regge model of J.M.Laget.
* 2<W<3 GeV, 2<Q^2<3 GeV^2
      sig3 = sig_gmh(thetacm,phicm,tt/1.e6,tprime/1.e6,vertex%q2/1.e6,
     *     invm/1.e3,epsilon)

      ntup%sigcm=sig3

      sigma_phi=ntup%sigcm

C DJG Breit-Wigner in the event generation - not needed here anymore.
*******************************************************************************
* Multiply by Breit_wigner [1/MeV] to give proper Mx weighted cross section
* d3sigma/dt/dphi_cm/dMx [ub/MeV^3/rad].
*
* gh - relativistic Breit-Wigner factor (eqn 38.52 of 2004 PDG book)
      Breit_wigner=(m_phi*gamma_phi)**2/
     *     ((mass**2-m_phi**2)**2 +(m_phi*gamma_phi)**2)

* gh - since sigma_omega is integrated over the omega peak, normalize
* Breit-Wigner to unit integral
      Breit_wigner=Breit_wigner/(gamma_phi*pi/2.)

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

CDJG      peep_omega=Jttheta*gamma_T*Breit_wigner*ntup.sigcm

C DJG Now we want the "fixed Mx" Jacobian since we are only generating the hadron
C DJG angles. No Breit-Wigner now either.
C DJG Note to me: If we ever want to extend to deuterium, need to get rid of numerous
C DJG hard-wired proton masses above...

      peep_phi = Jttheta_fx*gamma_T*ntup%sigcm

      if (debug(2)) write(6,*)' peep_phi: end ',peep_phi
      return 
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_phigmh(Q2,t,tmin,nu,w_gev,epsi)

* This routine calculates p(e,e'phi)p cross sections 
* in the form used in PYTHIA with modifications
* implmented in the HERMES Monte Carlo

      implicit none
      include 'constants.inc'

      real*8 t,tmin,tprime
      real*8 nu,epsi,Q2,w_gev
      real*8 sig0,sigt,sig219,R,cdeltatau,bphi

      tprime = abs(t-tmin)/1.e6
      cdeltatau = hbarc/(sqrt(nu**2+Q2+mphi2)-nu) !in fm!

* W dependence from photoproduction data
* Use Fig 35 of Cassel et al, PRD 24(1981) 2787

      sig0 = 0.225+(w_gev-2.0)*(0.275/2.)     ! microbarns

* R=L/T is usually fit to the form c_0*(Q2/M2_V)^c1
* c_0 and c_1 taken from HERMES data, Rakness thesis, Fig 5.8

      R = 0.38*(Q2/mphi2)**(0.87)
      if (R.gt.1) R=1

* The Q2 dependence is usually given by (M2_V/(Q2+M2_V))^2

      sigt = sig0*(1.0+epsi*R)*(mphi2/(Q2+mphi2))**2

* Need to parameterize t-dependence with b parameter as a function of c-tau
* Use Fig 34 of Cassel et al, PRD 24(1981) 2787

      bphi = 3.75 + 3.5*log10(cdeltatau)

c GH: the conventional formula with exponential factor is for forward
c     diffractive, it might not be appropriate for u-channel production

      sig219 = sigt*bphi*exp(-bphi*tprime)/2.0/pi !ub/GeV**2/rad
      
      sig_phigmh=sig219/1.d+06         !dsig/dtdphicm in microbarns/MeV**2/rad

      return
      end
