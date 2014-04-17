      real*8 function peep_omega(vertex,main)

C     This function is adapted from Pawel Ambrozewicz
C     noted by xucc
*     checked and modifed by gh - 05.02.11

      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f
C     For we need these quantities as starting point which may be different
C     from Pawel's noted by xucc

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 sigma_omega
c     real*8 m_omega/782.59/
c     real*8 gamma_omega/8.49/
      real*8 m_pisq, m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp

      real*8 mass               ! omega mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,uu              ! Mandelstam variables
      real*8 tprime
      real*8 e_photCM,e_omCM,t_min

      real*8 gamma_T            ! flux of transversely polarized virt. photons
c     real*8 Breit_wigner
      real*8 Jttheta            ! Jacobian t-Mx - Ep,theta_p in LAB
      real*8 Jttheta_fx         ! old Jacobian for fixed Mx
      real*8 tcos               ! cos of theta_LAB between Pp and q

c      real*8 mass_vtx
c      common /pawel_crosssection/mass_vtx

! Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      real*8 sig_joos,sig_gmh
      real*8 sig1,sig2

      if (debug(2)) write(6,*)' peep_omega: enter '

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

      m_pisq= 139.57018**2
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
      if (mass**2.lt.m_pisq) then
         peep_omega=0.0
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
      e_omCM   = (Wsq + mass**2 - m_psq)/invm/2.
      if ((e_omCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_omCM*e_photCM-
     *        sqrt((e_omCM**2-mass**2)*(e_photCM**2+qsq)))
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
* One Pion Exchange model used in DESY analysis of Joos et al.
* W<2 GeV, Q^2<1.4 GeV^2
      sig1 = sig_joos(mass/1.d3,tt/1.d6,uu/1.d6,
     *  vertex%q2/1.d6,invm/1.d3,main%epsilon)

* Parameterization based on saturated Regge model of J.M.Laget.
* 2<W<3 GeV, 2<Q^2<3 GeV^2
      sig2 = sig_gmh(thetacm,phicm,tt/1.d6,tprime/1.d6,vertex%q2/1.d6,
     *     invm/1.d3,epsilon)

      ntup%sigcm=sig2

      sigma_omega=ntup%sigcm

C DJG Breit-Wigner in the event generation - not needed here anymore.
*******************************************************************************
* Multiply by Breit_wigner [1/MeV] to give proper Mx weighted cross section
* d3sigma/dt/dphi_cm/dMx [ub/MeV^3/rad].
*
* gh - relativistic Breit-Wigner factor (eqn 38.52 of 2004 PDG book)
C DJG      Breit_wigner=(m_omega*gamma_omega)**2/
C DJG     *     ((mass**2-m_omega**2)**2 +(m_omega*gamma_omega)**2)

* gh - since sigma_omega is integrated over the omega peak, normalize
* Breit-Wigner to unit integral
C DJG      Breit_wigner=Breit_wigner/(gamma_omega*pi/2.)

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

      peep_omega = Jttheta_fx*gamma_T*ntup%sigcm

      if (debug(2)) write(6,*)' peep_omega: end ',peep_omega
      return 
      end


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_joos(mass,tt,uu,qsq,invm,epsilon)

* This routine calculates p(e,e'p)omega0 cross sections using the One
* Pion Exchange model referenced by Joos et al. Nucl.Phys. B122(1977)365-382.

C     Model cross section for omega production taken from P. Joos et al.
C     J. Dunne 4/20/00
      
C     This function is adapted from Pawel Ambrozewicz
C     noted by xucc
*     all quantities are in GeV, except for returned cross section, which is in
*     ub/MeV^2/rad

      implicit none
      include 'constants.inc'

C     COMMON /SIG_STUFF/ J,epsilon,gamma_T,Jttheta
c     real*4 J		! Jacobian E,E',Theta_e' ---> Q2,W,Phi_e'
C     We no longer need such J because we want to keep sampling in Lab 
C     system   Noted by xucc
     
      real*8 m_pi, m_pisq, m_p, m_psq, m_rho, m_rhosq

      real*8 mass		! omega mass

      real*8 invm		! invariant mass of hadronic system
      real*8 Wsq
      real*8 tt,uu              ! Mandelstam variables
      real*8 epsilon		! virtual photon polarisation
      real*8 E_gamma		! Energy of equiv. real photon for norm of dif.
      real*8 qsq		! 4 momentum transfer of scattered electron

      real*8 gamma_wpg          ! Partial width for omega->pi gamma
      real*8 Gsq/14.6/          ! pion-nucleon coupling
      real*8 Lamba_wp           
C     real*8 Q_N_sq,Q_NT_sq
      real*8 F_N		! form factor from pi-nucleon vertex
      real*8 F_w		! form factor from gamma-pi-omega vertex
      real*8 U_ff,pcm
      real*8 QT,QF		! momentum of on-shell and off-shell pion
				! in omega rest frame
      real*8 B,C,N_3	        ! parameters in Joos OPE model
      real*8 sig_ope_t          ! differential x-section for trans. photons
      real*8 sig_ope_l          ! differential x-section for long. photons
c     real*8 sig_dif		! differential x-section for diffraction
      real*8 R                  ! 

      real*8 Q_N,Q_NT           !,L_piece,T_piece
      real*8 temp,temp1,temp2

c     Convert masses to GeV
      m_pi  = Mpi/1.e3
      m_pisq= Mpi2/1.e6
      m_rho = Mrho/1.e3
      m_rhosq=Mrho2/1.e6
      m_p   = Mp/1.e3
      m_psq = Mp2/1.e6

      Wsq       = invm*invm

C     One Pion Exchange coefficients

      gamma_wpg = 0.9E-3*400*2*pi
      Lamba_wp  = sqrt(96*pi*(mass/(mass**2-m_pisq))**3*gamma_wpg)

      Q_N       = m_pi*m_pi*(m_pisq-4.*m_psq)/(4.*m_psq)
      Q_NT      = tt*(tt-4.*m_psq)/(4.*m_psq)
      F_N       = (1+8.41*Q_N)/(1+8.41*Q_NT)

      temp=0.0

c     QT        = Pcm(m_pisq,0.,mass**2)
      QT        = Pcm(m_pisq,temp,mass**2)

      QF        = Pcm(tt,-qsq,mass**2)
       
      temp1=2.3*QF
      temp2=2.3*QT

C     test1=(U_ff(temp1)/U_ff(temp2))

      F_w       = (QT/QF)**2*(U_ff(temp1)/U_ff(temp2))

      B         = -0.5*qsq*(tt-2.*m_psq)+0.25*(m_psq-qsq-Wsq)*
     $     (m_psq+mass**2-Wsq-tt)

      C         = -0.25*((invm-m_p)**2+qsq)*((invm+m_p)**2+qsq)

      N_3       = Wsq*tt*uu-tt*(m_psq+qsq)*(m_psq-mass**2)
     *      -m_psq*(qsq+mass**2)**2
      N_3       = 0.5*sqrt(N_3)
        
      R         = 0.4*qsq/m_rhosq
      E_gamma   = (Wsq-m_psq)/2./m_p
      
C     T_piece = (B+C)**2/(-C)

      sig_ope_t = 1./(Pcm(-qsq,Wsq,m_psq)*Pcm(temp,Wsq,m_psq)*Wsq)*Gsq*
     $     lamba_wp**2/16.*(-tt)/(tt-m_pisq)**2*m_rhosq*m_rhosq/
     $     (qsq+m_rhosq)**2*(4.*(B+C)**2 - 2.*qsq*N_3**2)/
     $     (4.*(-C))*F_N*F_w
      
c     L_piece = 2*qsq*N_3*N_3/(-C)

      sig_ope_l = 1./(Pcm(-qsq,Wsq,m_psq)*Pcm(temp,Wsq,m_psq)*Wsq)*Gsq*
     $     lamba_wp**2/16.*(-tt)/(tt-m_pisq)**2*m_rhosq*m_rhosq/
     $     (qsq+m_rhosq)**2*qsq*N_3*N_3/(-C)*F_N*F_w

c     Diffractive Cross Section
c     coded but not used since the diffractive contribution is miniscule

c     sig_dif   = Pcm(temp,Wsq,m_psq)/Pcm(-qsq,Wsq,m_psq)*(1.+epsilon*R)*
c     1    mass**4/(qsq+mass**2)**2*9.3*exp(6.7*tt)*(1+1.4/E_gamma)
      
******************************************************************************
*  we keep the tradition that ntup.sigcm be d2sigma/dt/dphi_cm [ub/MeV**2/rad]
*  2pi factor is for dphi_cm and 1.e6 converts GeV**2 to MeV**2
******************************************************************************
      sig_joos = (sig_ope_t+epsilon*sig_ope_l)/(2.*pi)/1.e6

      return 
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function U_ff(x)
      real*8 x
      U_ff = 1./(2.*x*x)*((2.*x*x+1)/(4.*x*x)*log(4*x*x+1.)-1)
      return
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function Pcm(m1sq,m2sq,m3sq)
      real*8 m1sq,m2sq,m3sq
      Pcm = sqrt((m1sq*m1sq-2.*m1sq*(m2sq+m3sq)+(m2sq-m3sq)**2)/
     *     (4.*m3sq))
      return
      end



      real*8 function sig_gmh(thetacm,phicm,t_gev,tprime_gev,
     *     q2_gev,w_gev,eps)

* This is a starting-point model to be used in the omega analysis of
* E01-004.  Fits were made to the p(e,e'p)omega Regge model calculation
* at Q^2=2.35, W=2.47 by J.M. Laget, Phys.Rev.D 70(2004)054023.  The
* model will be modified to best fit the omega data, as necessary.
* Subroutine calculates dsigma/dt/dphi_cm, which is returned as sig_gmh
* [ub/MeV^2/rad].
* gh, 05.02.11

      implicit none
      include 'constants.inc'

      real*8 sig,wfactor
      real*8 sigt,sigl,siglt,sigtt !components of dsigma/dt
      real*8 thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp
      real*8 lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0
      real*8 a,b,c,d,e,f,fpi,fpi235
      real*8 m_pi0sq

      m_pi0sq= Mpi02/1.e6

      tp = abs(tprime_gev)      ! just to make sure it's positive
      if (abs(t_gev)<tp) then
         write(6,*)' invalid -t>-tprime error',abs(t_gev),tp
         stop
      endif

* Normally, this would equal m_rho**2, but Laget has adjusted the value
* to reproduce the real photon data.
      lambda0_sq= 0.462         !GeV^2

* Pion saturating Regge trajectory.
* For t>(m_pi**2), the trajectory takes the usual form, 0.7(t-m_pi**2).
* As t --> -infty, the trajectory goes asymptotically to -1.
* Laget does not mention which function he uses for the asymptotic
* behavior, but I find that taking the hyperbolic tangent of the usual
* trajectory expression gives a curve which resembles the one in
* Laget's paper.
      if (t_gev>m_pi0sq) then
         alphapi_t=0.7*(t_gev-m_pi0sq)
      else
         alphapi_t=tanh(0.7*(t_gev-m_pi0sq))
      endif

* Use these instead of the usual Q^2 scaling of the response functions.
      alphapi_0=0.7*(0.-m_pi0sq)
      lambdapi_sq=lambda0_sq*((1.+alphapi_t)/(1.+alphapi_0))**2
      fpi    = 1./(1.+q2_gev/lambdapi_sq)
      fpi235 = 1./(1.+2.35/lambdapi_sq)
      
* Fit parameters to t-dependence of Laget's p(e,e'p)omega response
* functions at Q^2=2.35, W=2.47 [ub/GeV^2].  Before fitting, I first
* divided the response functions by (fpi/fpi235)^2 and any sin(thetacm)
* factors.
      a = 0.16675 
      b = 0.89524 
      c = 3.5991
      d = 6.4562
      e = -9.6199
      f = 5.8319
      sigl = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**0.25))
      sigl = sigl *(fpi/fpi235)**2

      a = -0.12286
      b = 0.56383
      c = 1.4673
      d = 2.1988
      e = 0.65170
      f = 18.501
      sigt = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**2))
      sigt = sigt *(fpi/fpi235)**2

      a = 0.46397
      b = 4.9179 
      c = 3.3023
      d = 3.1741
      siglt = a*exp(-b*tp)+c*exp(-d*(tp**0.5))
      siglt = siglt*(fpi/fpi235)**2*sin(thetacm)
* Laget uses -sqrt(e(1+e)) instead of +sqrt(2e(1+e))
      siglt = -siglt/sqrt(2.)

      a = -0.26497
      b = 2.7655 
      c = 2.8034
      d = 3.8586
      sigtt = a*exp(-b*tp)+c*exp(-d*(tp**0.5))
      sigtt= sigtt*(fpi/fpi235)**2*(sin(thetacm))**2

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor=(2.47**2-mp**2)**2/(w_gev**2-mp**2)**2

      sigl = sigl*wfactor
      sigt = sigt*wfactor
      siglt= siglt*wfactor
      sigtt= sigtt*wfactor

* Application of these Q^2, W scaling factors to my fitted functions
* gives response functions which closely resemble Laget's published
* curves at Q^2=2.35, W=2.47.  Comparing to the DESY data and Laget's
* model at Q^2=0.84, W=2.30, my fits are 20% high at -t=0.1, 300-500%
* high between -t=1 and 3 GeV^2.  Thus, this subroutine is optimized for
* Q^2>2, W>2.
      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt
      
      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad
      sig_gmh = sig

      return
      end
