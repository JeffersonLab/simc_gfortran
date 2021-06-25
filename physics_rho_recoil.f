      real*8 function peep_rho(vertex,main)

C     This function is adapted from Pawel Ambrozewicz
C     noted by xucc

      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f
C     For we need these quantities as starting point which may be different
C     from Pawel's noted by xucc

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

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
      real*8 tt,uu              ! Mandelstam variables
      real*8 tprime,t_min,t_max,uprime,u_min
      real*8 e_photCM,e_rhoCM,e_pCM

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Breit_wigner,Soding
      real*8 Jttheta            ! Jacobian t - theta_LAB
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

      real*8 sig_diffract,sig_hermes
      real*8 sig1,sig2

      if (debug(2)) write(6,*)' peep_rho_recoil: enter '

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

      m_pisq= Mpi2
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

c      mass = mass_vtx
      mass = targ%Mrec_struck
      if (mass**2.lt.m_pisq) then
         peep_rho=0.0
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
c     main%t =tt

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )

      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_rhoCM  = (Wsq + mass**2 - m_psq)/invm/2.
      e_pCM     = E_cm
      
c gh - 18.02.04
c although t is single valued (negative), this is not true for u
c right at the kinematic endpoint, u turns positive, so umin is positive
c even though for most events u is negative.  This is surprising, but is
c confirmed by checking against s+t+u relation.
      
      if ((e_rhoCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_rhoCM*e_photCM
     *        -sqrt( (e_rhoCM**2-mass**2)*(e_photCM**2+qsq) ))
         t_max = -qsq + mass**2 -2.*(e_rhoCM*e_photCM
     *        +sqrt((e_rhoCM**2-mass**2)*(e_photCM**2+qsq)))
         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM
     *        -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))
      else
         write(6,*)' physics_rho_recoil: no valid t,u min '
      endif
      
      tprime = abs(tt-t_min)
      uprime = abs(uu-u_min)
      main%tmin=t_min
      
      if (abs(tt).lt.abs(t_min)) then
         write(6,*)' unphysical -t<-t_min ',tt,t_min,tprime
      endif
      if (abs(tt).gt.abs(t_max)) then
         write(6,*)' unphysical -t>-t_max ',tt,t_max,tprime
      endif
      if (uu.gt.u_min) then
         write(6,*)' unphysical -u<-u_min ',uu,u_min,uprime
      endif
      

******************************************************************************
*  we keep the tradition that ntup.sigcm is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
******************************************************************************
*
* Diffractive rho production model used by P. Ambrosewicz
c      sig1 = sig_diffract(mass/1.e3,invm/1.e3,qsq/1.e6,
c     *     tt/1.e6,t_min/1.e6)

* PYTHIA model with modifications from the HERMES MC
      sig2 = sig_hermes(mass/1.e3,qsq,tt,t_min,uu,u_min,nu,invm/1.e3,
     *       epsilon,thetacm)
      
      ntup%sigcm=sig2

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
      Breit_wigner=(Mrho*MrhoW)**2/
     *     ((mass**2-Mrho**2)**2 +(Mrho*MrhoW)**2)

* gh - since sigma_omega is integrated over the omega peak, normalize
* Breit-Wigner to unit integral
      Breit_wigner=Breit_wigner/(MrhoW*pi/2.)

* Soding model (skewness) correction
*
* The value of the exponent falls steeply with t (see Fig 71
* of Bauer's review article, Rev.Mod.Phys. 50(1978)261) and is often
* determined by fit to data.  In the CLAS experiment, the fit exponent was
* 0+/-0.09 [C. Hajidakis Thesis, p. 113]
*
* The value used by Ambrosewicz (low -t) was 5.2
*      Soding=(m_rho/mass)**5.2
*
* Initial value determined from fit to E01-004 (high -t) data
* Value is 4.95 for Q2=1.60 data (-t=3.95) and 6.65 for Q2=2.45 data (-t=4.67)
*      Soding=(mass/m_rho)**(-4.394+2.366*abs(tt/1.e6))

* GH: leave it fixed at the CLAS value for now
       Soding=(Mrho/mass)**0.0

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
* by multiplying  by virtual transverse photon flux factor, gamma_T [1/MeV/sr]
*
* gh - checked and this agrees with gtpr in physics_pion.f
      gamma_T   = (alpha*E_prime*(Wsq-m_psq))/(4.*pi*pi*E0*m_p*qsq*
     $     (1-epsilon))

*******************************************************************************
* Lab differential cross section
* ub/MeV^2/rad-->ub/MeV^3/rad-->ub/MeV/sr-->ub/MeV^2/sr^2
*******************************************************************************

CDJG BW in event generation - also want Jacobian at fixed Mx ...
CDJG May want (Mrho/Mass)**5.2 as an extra weight here. Need to check with Garth.

      peep_rho=Jttheta_fx*gamma_T*Soding*ntup%sigcm

      return 
      end



C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_diffract(mass,invm,qsq,tt,t_min)

* This routine calculates p(e,e'p)rho0 cross sections using the diffractive
* model referenced by P. Ambrosewicz et al.

* all calculations are in GeV

      implicit none
      include 'constants.inc'

      real*8 mass,invm,Wsq,qsq,tt,t_min,e_gamma
      real*8 A,B,C,b_slope,sig_grho
      real*8 m_pgev,m_pgevsq

      m_pgev   = Mp/1.e3
      m_pgevsq = m_pgev**2

      Wsq  = invm*invm

      A = 0.6*Wsq - 3.62*invm + 5.64
      B = -2.58*Wsq + 14.93*invm - 23.21
      C = 1.1*invm + 2.6
      b_slope = A*qsq + B*sqrt(qsq) + C

      if (b_slope.lt.0.0) then
         b_slope = 0.0
      endif

      e_gamma  = (Wsq - m_pgevsq)/m_pgev/2.
      sig_grho = 29.4/e_gamma + 9.5

* returned cross section is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
* 2pi factor is for dphi_cm and 1.e6 converts GeV^2 to MeV^2

      sig_diffract = (sig_grho/(1+qsq/mass**2)**2)
     *     *exp(b_slope*(tt-t_min))/(2.*pi)/1.e6

      return
      end



C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_hermes(mass,Q2,t,tmin,u,umin,nu,w_gev,epsi,
     1     thetacm)

* This routine calculates p(e,e'rho)p cross sections using D. Gaskell's
* attempt to code in the form used in PYTHIA with modifications
* implmented in the HERMES Monte Carlo by Patricia Liebing.

      implicit none
      include 'constants.inc'

      real*8 t,tmin,tprime,u,umin,uprime
      real*8 mass,nu,epsi,Q2,w_gev,wfactor,m_pgev,thetacm
      real*8 sig0,sigt,sig219,R,cdeltatau,brho

      m_pgev = Mp/1.e3
      tprime = abs(t-tmin)/1.e6
      uprime = abs(u-umin)/1.e6
      cdeltatau = hbarc/(sqrt(nu**2+Q2+mrho2)-nu) !in fm!

* Put in some W dependence from photoproduction data

* DJG:  This is my rough fit to some old photoproduction data.	
c      sig0 = 41.263/(nu/1000.0)**0.4765     ! microbarns

* GH:   Photoproduction fit in Fig 22 of Cassel et al, PRD 24(1981) 2787
       sig0 = 29.4/(nu/1000.)+9.5            ! microbarns

* GH: Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
c       wfactor=(2.47**2-m_p**2)**2/(w_gev**2-m_p**2)**2

* DJG:  R=L/T is usually fit to the form c_0 (Q2/M2_rho)^c1
* GH:   c_0 and c_1 taken from HERMES data, Rakness thesis, Fig 5.9

      R = 0.39*(Q2/Mrho2)**(0.68)
      if (R.gt.1) R=1

* DJG:  The Q2 dependence is usually given by (M2_rho/(Q2+M2_rho))^2
* DJG:  HERMES found that 2.575 works better than 2 in the exponent.

      sigt = sig0*(1.0+epsi*R)*(Mrho2/(Q2+Mrho2))**(2.575)

* DJG: Need to parameterize t-dependence with b parameter as a function
*      of c-tau

c      if (cdeltatau.le.0.390323) then
c         brho = 1.07477	      
c      else if(cdeltatau.lt.1.968.and.cdeltatau.gt.0.390323) then
c         brho = 4.4679 + 8.6106*dlog10(cdeltatau)
c      else
c         brho = 7.0
c      endif
      
* GH:  brho parameterization from Cynthia Hajidakis' thesis Fig 4.9,  p146.
*      JLab CLAS E99-105 (Orsay)

      if (cdeltatau.lt.2.057) then
         brho = -0.0941+3.449*cdeltatau
      else
         brho = 7.0  ! see Fig 14 of Cassel et al, PRD 24(1981) 2787
      endif
      if (brho.lt.1.07477) brho = 1.07477
      
c GH: the conventional formula with exponential factor is for forward
c     diffractive, it might not be appropriate for u-channel production
      
      if (thetacm.gt.pi/2.) then ! t-channel
         sig219 = sigt*brho*exp(-brho*tprime)/2.0/pi !ub/GeV**2/rad
      else                      ! u-channel
         sig219 = sigt*brho*exp(-brho*uprime)/2.0/pi
         sig219 = sig219/10.    ! back angle peak is ~10% of forward angle peak
      endif
         
      sig_hermes=sig219/1.e+06         !dsig/dtdphicm in microbarns/MeV**2/rad

c GH: check for weird behavior on upper side of rho peak (t-channel)
      if (thetacm.gt.pi/2. .and. mass.gt.1.2 .and. tprime.gt.Q2 .and.
     *    sig_hermes.gt.1.e-14) then
         write(6,*)' tprime=',tprime,' sig=',sig_hermes
         sig_hermes=sig_hermes*1.e-6
      endif
      
      return
      end
