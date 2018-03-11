      real*8 function peep_dvcs(vertex,main)

c     gh - 18.02.26
      
      implicit none
      include 'simulate.inc'

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 sigma_dvcs
      real*8 m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp
  ! low Q2
      real*8 mass               ! pi0 mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,uu              ! Mandelstam variables
      real*8 tprime,t_min,t_max,uprime,u_min,uu2,umin2
      real*8 e_photCM,e_gammafCM,e_pCM,p_gammafCM

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Breit_wigner
      real*8 Jttheta            ! Jacobian t-Mx - Ep,theta_p in LAB
      real*8 Jttheta_fx         ! old Jacobian for fixed Mx
      real*8 tcos               ! cos of theta_LAB between Pp and q
      real*8 vertex_gammaf_x,vertex_gammaf_y,vertex_gammaf_z
      real*8 P_gammaf,E_gammaf,ucos,uu3

! Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      real*8 sig_dvcsgmh
      real*8 sig1,sig2,sig3

      if (debug(2)) write(6,*)' peep_dvcs: enter '

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
     >     +vertex%up%z*vertex%uq%z
      
* since here the meson is the recoil particle, t for us is not the same thing as 
* main%t, so overwrite it
* for the + sign I spent 1 hour to check. xucc
      tt = 2.0*(m_p**2-m_p*Ep)
c     main%t = tt

      uu = -qsq +m_psq +2.*(qvec*Pp*tcos -nu*Ep )
      
c     use momentum conservation to get components of recoil meson
      vertex_gammaf_x = qvec*vertex%uq%x - Pp*vertex%up%x
      vertex_gammaf_y = qvec*vertex%uq%y - Pp*vertex%up%y
      vertex_gammaf_z = qvec*vertex%uq%z - Pp*vertex%up%z
      P_gammaf=sqrt(vertex_gammaf_x**2 +vertex_gammaf_y**2 
     1         +vertex_gammaf_z**2)
      E_gammaf=P_gammaf
      
c      ucos =  vertex_gammaf_x/P_gammaf*vertex%uq%x
c     1       +vertex_gammaf_y/P_gammaf*vertex%uq%y
c     2     +vertex_gammaf_z/P_gammaf*vertex%uq%z

      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_gammafCM  = (Wsq - m_psq)/invm/2.
      p_gammafCM  = e_gammafCM
      e_pCM    = (Wsq + m_psq)/invm/2.

c gh - 18.02.04
c although t is single valued (negative), this is not true for u
c right at the kinematic endpoint, u turns positive, so umin is positive
c even though for most events u is negative.  This is surprising, but is
c confirmed by checking against s+t+u relation.

      if ((p_gammafCM**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq -2.*(e_gammafCM*e_photCM
     *        -sqrt( (p_gammafCM**2)*(e_photCM**2+qsq) ))
         t_max = -qsq -2.*(e_gammafCM*e_photCM
     *        +sqrt((p_gammafCM**2)*(e_photCM**2+qsq)))
         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM
     *        -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))
      else
         write(6,*)' physics_dvcs: no valid t,u min '
      endif

c check u using s+t+u relation
c      uu2=-qsq+2.*m_psq+mass**2-invm**2-tt
c      umin2=-qsq+2.*m_psq+mass**2-invm**2-t_max

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
c      if (uprime.gt.0.) then
c         write(6,*)' GOOD -u>-u_min ',uu,u_min,uprime,ucos
c      endif
c      write(6,*)' tt ',tt,t_min,t_max,tprime
      

******************************************************************************
*  we keep the tradition that ntup.sigcm is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
******************************************************************************
*
      sig1 = sig_dvcsgmh(Ebeam/1.e3,qsq/1.e6,tt/1.e6,t_min/1.e6,
     *       uu/1.e6,u_min/1.e6,invm/1.e3)

      ntup%sigcm=sig1

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6

c      sigma_dvcs=ntup%sigcm

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

      peep_dvcs = Jttheta_fx*gamma_T*ntup%sigcm

      if (debug(2)) write(6,*)' peep_dvcs: end ',peep_dvcs
      return 
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_dvcsgmh(Eb_g,Q2_g,t_gev,tmin,u_gev,umin,
     1       W_gev)

c This routine calculates p(e,e'gamma)p DVCS cross sections, based on 
C the Kumericki and Mueller KM15 model at xB=0.36, which includes the
c 2015 CLAS and Hall A DVCS data in fixing the GPD parameters.

c Subroutine calculates dsigma/dt/dphi_cm, which is returned as sig_dvcsgmh
c [ub/MeV^2/rad].
c gh - 18.02.26

      implicit none
      include 'constants.inc'

      real*8 t_gev,tmin,tprime,u_gev,umin,uprime
      real*8 Eb_g,mass,eps,Q2_g,W_gev,wfactor,m_p
      real*8 Q2hi,Q2lo,delQ2,Q2tmp
      real*8 Whi,Wlo,Wfac_hi,Wfac_lo
      real*8 sighi,siglo,sighiW,sigloW
      real*8 thetacm,phicm,sig_phi0gmh
      real*8 sig                                ! dsigma/dt
      
c use fit parameters for KM15 model at xB=0.36 (nb/GeV2)
      real*8 Q2tab(4) / 2.30,  3.00,  3.80,  5.00/
      real*8 Wtab(4)  / 2.23,  2.49,  2.76,  3.125/
      real*8 Ebtab(3) / 6.60,  8.80, 10.90/

c sig_unp is parameterized in the form p1*(-t-p3)/(-t-p2)**2+p4
c fill parameter arrays in the order (x=Ee, y=Q2)

c parameters for fit vs t (not tprime)
c Q2=2.30 for 3 beam energies
      real*8 p1(3,4) / 0.004, 0.005, 0.005, 
c Q2=3.0
     1                 0.004, 0.006, 0.005, 
c Q2=3.8
     2                -0.057, 0.007, 0.005, 
c Q2=5.0
     3                0.0003,0.0003, 0.006/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p2(3,4) / 0.098, 0.074, 0.062, 
c Q2=3.0
     1                 0.107, 0.062, 0.055, 
c Q2=3.8
     2                -0.057, 0.061, 0.058, 
c Q2=5.0
     3                 0.127, 0.127, 0.057/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p3(3,4) / 0.158, 0.136, 0.124, 
c Q2=3.0
     1                 0.170, 0.136, 0.127, 
c Q2=3.8
     2                -0.057, 0.140, 0.133, 
c Q2=5.0
     3                 0.275, 0.275, 0.137/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p4(3,4) / 0.045, 0.015, 0.008, 
c Q2=3.0
     1                 0.051, 0.007, 0.002, 
c Q2=3.8
     2                 0.359, 0.006, 0.001, 
c Q2=5.0
     3                 0.062, 0.062, 0.000/ ! use Eb=8.8 params for 6.6

      integer iflag,pcount,ndat,Q2count,Q2c,Ebc,Ebcount

      iflag=1                   ! flag for t(0) or u(1) channel
      ndat=4
      
      m_p = Mp/1.e3
      tprime = abs(t_gev-tmin)
      uprime = abs(u_gev-umin)

      Q2tmp = Q2_g
      if (Q2tmp.lt.Q2tab(1)) then
         if (pcount.le.500) then
            write(6,100)Q2tmp
 100        format(' WARNING: Q2 extrapolated below model range ',f7.3)
            pcount=pcount+1
         endif
         Q2tmp = Q2tab(1)
      elseif(Q2tmp.gt.Q2tab(ndat)) then
         if (pcount.lt.500) then
            write(6,101)Q2tmp
 101        format(' WARNING: Q2 extrapolated above model range ',f7.3)
            pcount=pcount+1
         endif
         Q2tmp = Q2tab(ndat)
      endif
      
      Ebcount=0
      do Ebc=1,3
         if (abs(Eb_g-Ebtab(Ebc)).lt.0.5) then
            Ebcount=Ebc
         endif
      enddo
      if (Ebcount.lt.1 .or. Ebcount.gt.3) then
         write(6,*)' Ebcount error ',Ebcount,Eb_g
      endif

c calculate hi,lo cross sections
      Q2hi=0.
      Q2lo=0.
      delQ2=1.
      Q2count=0
      if( Q2tmp.lt.Q2tab(1) ) then
         Q2count = 1
         Q2hi = Q2tab(2)
         Whi  = Wtab(2)
         Q2lo = Q2tab(1)
         Wlo  = Wtab(1)
         delQ2 = (Q2hi - Q2lo)
      else
         do Q2c=1,(ndat-1)
            if( (Q2tmp.ge.Q2tab(Q2c)).and. (Q2tmp.lt.Q2tab(Q2c+1) )
     1           .or. Q2tmp.ge.Q2tab(Q2c+1) ) then
               Q2count = Q2c
               Q2hi = Q2tab(Q2count+1)
               Whi  = Wtab(Q2count+1)
               Q2lo = Q2tab(Q2count)
               Wlo  = Wtab(Q2count)
               delQ2 = (Q2hi - Q2lo)
            endif               !Q2 check
         enddo                  !Q2
      endif
      
      if (iflag.lt.1) then      ! t-channel
         sighi = p1(Ebcount,Q2count+1)*
     1         (tprime+abs(tmin)-p3(Ebcount,Q2count+1))
     2        /(tprime+abs(tmin)-p2(Ebcount,Q2count+1))**2
     3        +p4(Ebcount,Q2count+1)
         
         siglo = p1(Ebcount,Q2count)*
     1         (tprime+abs(tmin)-p3(Ebcount,Q2count))
     2        /(tprime+abs(tmin)-p2(Ebcount,Q2count))**2
     3        +p4(Ebcount,Q2count)
         
      else                      ! u-channel
c christian weiss recommends the following change for u-channel:
c switch u-slope for t-slope, then divide by 10, since back angle peak 
c is ~10% of forward angle peak (at least for omega electroproduction)
         sighi = ( p1(Ebcount,Q2count+1)*
     1         (uprime+abs(tmin)-p3(Ebcount,Q2count+1))
     2        /(uprime+abs(tmin)-p2(Ebcount,Q2count+1))**2
     3        +p4(Ebcount,Q2count+1) )/10.
         
         siglo = ( p1(Ebcount,Q2count)*
     1        (abs(uprime)+abs(tmin)-p3(Ebcount,Q2count))
     2        /(abs(uprime)+abs(tmin)-p2(Ebcount,Q2count))**2
     3        +p4(Ebcount,Q2count) )/10.
      endif
      
c sighi,lo are at different W.  scale both to the W needed for the event
      Wfac_hi= ((Whi**2-m_p**2)**2) / ((W_gev**2-m_p**2)**2)
      Wfac_lo= ((Wlo**2-m_p**2)**2) / ((W_gev**2-m_p**2)**2)
      
      sighiW = sighi*Wfac_hi
      if (sighiW.lt.0.) then
         write(6,*)' dvcs: sighiW<0 ',sighiW,uprime,abs(tmin),Wfac_hi
         sighiW=0.
      endif
      
      sigloW = siglo*Wfac_lo
      if (sigloW.lt.0.) then
         write(6,*)' dvcs: sigloW<0 ',sigloW,uprime,abs(tmin),Wfac_lo
         sigloW=0.
      endif

c interpolate to get cross section at Q2 needed for the event
c units are nb/GeV2 (dsig/dt)
      
      if( Q2count.le.(ndat-1) .and. Q2tmp.ge.Q2tab(Q2count) .and.
     1     Q2tmp.lt.Q2tab(Q2count+1) ) then
         
         sig  = ( sigloW*(Q2hi-Q2tmp)+ sighiW*(Q2tmp-Q2lo))/delQ2

c        write(6,*)' sig1 ',Q2tmp,sig,sigloW,sighiW
      
      elseif (Q2tmp.ge.Q2tab(ndat) ) then
         
         sig =  sighiW+ (sighiW-sigloW) /delQ2
         
c         write(6,*)' sig2 ',Q2tmp,sig,sigloW,sighiW

      elseif (Q2tmp.le.Q2tab(1) ) then
         
         sig =  sigloW- (sighiW-sigloW) /delQ2
         
c         write(6,*)' sig3 ',Q2tmp,sig,sigloW,sighiW

      else

         write(6,*)' Q2tmp error ',Q2tmp,Q2count

      endif
      
      sig_dvcsgmh = sig/2./pi*1.e-09 !dsig/dtdphicm in microbarns/MeV^2/rad

      if (sig_dvcsgmh .gt. 1.e-6) then
         write(6,*)' sig warning ',sig,Q2tmp,W_gev
         write(6,*)uprime,abs(tmin),p2(Ebcount,Q2count+1)
         write(6,*)sighiW,Wfac_hi,sigloW,Wfac_lo
      endif
      
      return
      end
