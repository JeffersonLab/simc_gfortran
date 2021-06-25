      real*8 function peep_phi(vertex,main)

c     gh - 14.12.18
c     gh - modified 17.12.27
      
      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 sigma_phi
      real*8 m_p, m_psq

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
      real*8 tprime,t_min,t_max,uprime,u_min,uu2,umin2
      real*8 e_photCM,e_phiCM,e_pCM

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

      real*8 sig_gmh,sig_diffract,sig_phigmh,sig_phiweiss
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
      if (mass**2.lt.Mk2) then
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
      tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y
     >       +vertex%up%z*vertex%uq%z
      
* since here the meson is the recoil particle, t for us is not the same thing as 
* main%t, so overwrite it
* for the + sign I spent 1 hour to check. xucc
      tt = 2.0*(m_p**2-m_p*Ep)
c     main%t = tt

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )
      
      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_phiCM   = (Wsq + mass**2 - m_psq)/invm/2.
      e_pCM     = (Wsq + m_psq - mass**2)/invm/2.

c gh - 18.02.04
c although t is single valued (negative), this is not true for u
c right at the kinematic endpoint, u turns positive, so umin is positive
c even though for most events u is negative.  This is surprising, but is
c confirmed by checking against s+t+u relation.
      
      if ((e_phiCM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_phiCM*e_photCM
     *        -sqrt( (e_phiCM**2-mass**2)*(e_photCM**2+qsq) ))
         t_max = -qsq + mass**2 -2.*(e_phiCM*e_photCM
     *     +sqrt((e_phiCM**2-mass**2)*(e_photCM**2+qsq)))
         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM
     *        -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))
      else
         write(6,*)' physics_phi: no valid t,u min '
      endif

c check u using s+t+u relation
c      uu2=-qsq+2.*m_psq+mass**2-invm**2-tt
c      t_max = -qsq + mass**2 -2.*(e_phiCM*e_photCM
c     *     +sqrt((e_phiCM**2-mass**2)*(e_photCM**2+qsq)))
c      umin2=-qsq+2.*m_psq+mass**2-invm**2-t_max
      
      tprime = abs(tt-t_min)
      uprime = abs(uu-u_min)
      main%tmin=t_min
      
      if (abs(t_min).gt.abs(tt)) then
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
* PYTHIA model with modifications from the HERMES MC
c      sig1 = sig_phigmh(mass/1.e3,qsq,tt,t_min,uu,u_min,nu,invm/1.e3,
c     *       epsilon)

      sig2 = sig_phiweiss(qsq/1.e6,tt/1.e6,t_min/1.e6,t_max/1.e6,
     *       uu/1.e6,u_min/1.e6,invm/1.e3,epsilon,thetacm)

c      write(6,*)' phi ',sig1,sig2
      
      ntup%sigcm=sig2

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6

c      sigma_phi=ntup%sigcm

C DJG Breit-Wigner in the event generation - not needed here anymore.
*******************************************************************************
* Multiply by Breit_wigner [1/MeV] to give proper Mx weighted cross section
* d3sigma/dt/dphi_cm/dMx [ub/MeV^3/rad].
*
* gh - relativistic Breit-Wigner factor (eqn 38.52 of 2004 PDG book)
      Breit_wigner=(Mphi*MphiW)**2/
     *     ((mass**2-Mphi**2)**2 +(Mphi*MphiW)**2)

* gh - since sigma_omega is integrated over the omega peak, normalize
* Breit-Wigner to unit integral
      Breit_wigner=Breit_wigner/(MphiW*pi/2.)

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
      real*8 function sig_phigmh(mass,Q2,t,tmin,u,umin,nu,w_gev,epsi,thcm)

* This routine calculates p(e,e'phi)p cross sections 
* in the form used in PYTHIA with modifications
* implmented in the HERMES Monte Carlo

      implicit none
      include 'constants.inc'

      real*8 t,tmin,tprime,u,umin,uprime,thcm
      real*8 mass,nu,epsi,Q2,w_gev,wfactor,m_p
      real*8 sig0,sigt,sig219,R,cdeltatau,bphi

      m_p = Mp/1.e3
      tprime = abs(t-tmin)/1.e6
      uprime = abs(u-umin)/1.e6
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

      if (thcm.gt.pi/2.) then   ! t-channel
         sig219 = sigt*bphi*exp(-bphi*tprime)/2.0/pi !ub/GeV**2/rad
      else                      ! u-channel
         sig219 = sigt*bphi*exp(-bphi*uprime)/2.0/pi
         sig219 = sig219/10.    ! back angle peak is ~10% of forward angle peak
      endif
         
      sig_phigmh=sig219/1.e+06         !dsig/dtdphicm in microbarns/MeV**2/rad

c GH: check for weird behavior on upper side of peak (t-channel)
      if (thcm.gt.pi/2. .and. mass.gt.1.5 .and. tprime.gt.Q2 .and.
     *    sig_phigmh.gt.1.e-20) then
          write(6,*)' tprime=',tprime,' sig=',sig_phigmh
          sig_phigmh=sig_phigmh*1.e-6
       endif
      
      return
      end



C---------------------------------------------------------------------
      real*8 function sig_phiweiss(q2_gev,t_gev,tmin,tmax,u_gev,umin,
     *  w_gev,epsi,thcm)

      implicit none
      include 'constants.inc'
      
      real*8 q2_gev,t_gev,tmin,tmax,u_gev,umin,w_gev,epsi,thcm
      real*8 sig,sigt,sigl
      real*8 ttmin,ttmax
      integer itdep,itmin

      itdep=1
      itmin=0

      call xphi(sigt, sigl, ttmin, ttmax, w_gev, t_gev, u_gev, q2_gev,
     >     itdep, itmin, thcm)

c check tmin, tmax values
      if (abs(ttmin-tmin).lt.1.e-10) then
         write(6,*)' sig_phiweiss: tmin disagreement ',tmin,ttmin
      endif
      if (abs(ttmax-tmax).lt.1.e-10) then
         write(6,*)' sig_phiweiss: tmax disagreement ',tmax,ttmax
      endif
         
      sig = sigt +epsi*sigl

      sig_phiweiss = sig/2./pi*1.e-09    !dsig/dtdphicm in microbarns/MeV^2/rad

      return
      end

      
c---------------------------------------------------------------------
      subroutine xphi(dxt, dxl, tmin,tmax, w, t, u, qq,
     >    itdep, itmin, thcm)

c     subroutine xphi(xt, xl, dxt, dxl, tmin, tmax, b, uug, w, t,
c     *   qq, itdep, itmin)

c     c. weiss (weiss.at.jlab.org) 
c     v1 07mar12
c     v2 09mar12  argument w instead of s, tmax, input checking 
c
c     exclusive phi meson electroproduction  gamma* + n  to  phi + n
c     parametrization of virtual photoproduction cross sections
c     (hand convention)
c
c     masses/energies/momenta in GeV
c
c     output:
c     xt     total  transv. cross secn   sigma-t    (nb)
c     xl     total  longit. cross secn   sigma-l    (nb)
c     dxt    diff.  transv. cross secn  dsigma-t/dt (nb/gev**2)
c     dxl    diff.  longit. cross secn  dsigma-l/dt (nb/gev**2)
c     tmin   kinematic upper limit of t
c     tmax   kinematic lower limit of t
c     b      exponential t-slope
c     uug    equivalent dipole mass squared
c
c
c     input:
c     w      gamma*-n cm energy w = sqrt(s)       (gev)
c     t      invariant momentum transfer t .le. 0 (gev**2) 
c     qq     photon virtuality q**2 .ge. 0        (gev**2)
c
c
c     itdep  selects parametrization of t-dependence
c            1   exponential
c            2   dipole
c
c     the dipole mass parameter is calculated from the exponential
c     slope by matching the profiles at a fixed value t = t0.
c
c
c     itmin  determines check of kinematic bounds in t and actual value 
c     of t used in calculation of differential cross secn
c            0   t-actual = t. no check of bounds, allows t = 0
c            1   t-actual = t. return zero if  t.ge.tmin  or  t.le.tmax
c            2   t-actual = tmin
c
c     if  itmin = 0  is selected, the routine allows one to evaluate the 
c     differential cross secn parametrization for any fixed  t,
c     including unphysical values such as the point  t = 0. 
c
c     if  itmin = 1  is selected, the routine returns zero 
c     differential cross secn if t is outside of kinematic bounds
c 
c     if  itmin = 2  is selected, the routine evaluates the differential 
c     cross secns at  t = tmin, where  tmin  is calculated at the given
c     values of  w  and  qq. in this case the actual value of  t  on input 
c     is ignored. no change of the actual value of the argument  t  
c     is made.
c
c
c     if the energy is below threshold,  w .lt. wth = 1.958, 
c     or if  qq .lt. 0 (defined numerically as qq .lt. (-eps) = -1.d-15)
c     the routine exits with zero values for all output variables.
c
c
c     implicit double precision (a - h, o - z)
      implicit real*8 (a - h, o - z)
      implicit integer (i - n)

c
c     ...numerical limit for t, qq check
c
      parameter(eps = 1.e-15)
c
c     ...masses
c
      parameter(un = 0.938e0, uphi = 1.020e0)
c
      wth = un + uphi
c
c     ...check w, qq input and return zero if out of bounds
c
      if (w.lt.wth.or.qq.lt.(-eps)) then
         tmin = 0.e0
         tmax = 0.e0
         xt   = 0.e0
         xl   = 0.e0
         b    = 0.e0
         uug  = 0.e0
         dxt  = 0.e0
         dxl  = 0.e0
         return
      endif
c
c     ...cross sections as function of w
c
      al1 = 400.  e0 
      al2 =   1.  e0 
      al3 =   0.32e0 
      unt =   3.  e0
      cr  =   0.4 e0
c
      xt = al1 *(1.e0 - wth**2/w**2)**al2 *w**al3
     *     /(1.e0 + qq/uphi**2)**unt
c
      xl = cr*qq/uphi**2 *xt
c
c     ...calculate tmin, tmax
c
      call ftcm(w, qq,  1.e0, un, uphi, tmin)
      call ftcm(w, qq, -1.e0, un, uphi, tmax)
c
c     ...check t input and set actual value of t
c
      if      (itmin.eq.0) then
c
         ta = t
c
      else if (itmin.eq.1) then
c
         if (t.lt.tmin.and.t.gt.tmax) then
            ta = t
         else
            b   = 0.e0
            uug = 0.e0
            dxt = 0.e0
            dxl = 0.e0
            return
         endif
c
      else if (itmin.eq.2) then
c
         ta = tmin
c
      else
         stop 'itmin out of bounds'
      endif
c
c     ...exponential t-slope as function of w
c
      b0    = 2.2  e0 
      alphp = 0.24 e0 
c
      b     = b0 + 4.e0*alphp*log(w)  
c
c     ...equivalent dipole mass squared
c
      t0  = -0.5e0 
      uug = -t0/(exp(-b*t0/4.e0) - 1.e0)
c
c     ...normalized t-profile
c
c        itdep = 1 exponential profile
c                2 dipole      profile
c
      if      (itdep.eq.1) then

         if (thcm.gt.pi/2.) then   ! t-channel
         
            f  = exp(b*ta)
            fi = exp(b*tmin)/b
     *           - exp(b*tmax)/b

c            write(6,*)' tt ',tmin,t,tmax
c            write(6,*)f,fi,f/fi

         else                   ! u-channel

c gh - 17.12.27
c christian recommends the following change for u-channel:
c switch u-slope for t-slope, then divide by 10, since back angle peak 
c is ~10% of forward angle peak (at least for omega electroproduction)   
            
            f  = exp(b*u)
            umin=-qq+2.*un**2+uphi**2-w**2-tmax
            umax=-qq+2.*un**2+uphi**2-w**2-tmin
            fi = exp(b*umin)/b
     *           - exp(b*umax)/b
            f = f/10.
            
c            write(6,*)' uu ',umin,u,umax
c            write(6,*)f,fi,f/fi

         endif
c
      else if (itdep.eq.2) then
c
         f   = uug**4/(uug - ta  )**4 
         fi  = uug**4/(uug - tmin)**3/3.e0
     *       - uug**4/(uug - tmax)**3/3.e0 
c
      else
         stop 'itdep out of bounds'
      endif
c
c     ...differential cross sections
c
      dxt = xt*f/fi
      dxl = xl*f/fi
      
      sigw=1.0
      
      end
c
c---------------------------------------------------------------------
c
      subroutine ftcm(w, qq, costh, un, um, t)
c
c     calculate  t  in meson electroproduction
c     using exact kinematics in cm frame
c
c     initial and final baryon mass are the same
c
c     tmin is obtained for costh =  1
c     tmax                 costh = -1
c
c      implicit double precision (a - h, o - z)
      implicit real*8 (a - h, o - z)
c
      call lamfn(alam1, w**2, un**2,   -qq)
      call lamfn(alam2, w**2, un**2, um**2)
      p1  = sqrt(alam1)/2.e0/w
      p2  = sqrt(alam2)/2.e0/w
      e1  = sqrt(p1**2 + un**2)
      e2  = sqrt(p2**2 + un**2)
c
      t = 2*un**2 - 2*e1*e2 + 2*costh*p1*p2
c
      end
c
c---------------------------------------------------------------------
c
      subroutine lamfn(alam, x, y, z)
c
c     lambda function (invariant expression of cm momentum)
c
c      implicit double precision (a - h, o - z)
      implicit real*8 (a - h, o - z)
c
      alam = x**2 + y**2 + z**2 - 2*x*y - 2*x*z - 2*y*z
c
      end
c
c---------------------------------------------------------------------
c
