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
      real*8 xB                 ! Bjorken x
      real*8 Ep,Pp
  ! low Q2
      real*8 mass               ! pi0 mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,uu,u2           ! Mandelstam variables
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

      real*8 sig_dvcsgmh,sig_tgvdvcs
      real*8 sig1,sig2,sig3

      integer model
      logical first
      save first

      data model /2/  ! model: 1=dvcsgmh 2=tgvdvcs
      data first /.TRUE./

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
c check u using s+t+u relation
         u2=-qsq+2.*m_psq-invm**2-tt
         if (abs(u2-uu).gt.1.e-6) write(6,*)' dvcs: u error ',u2,uu

* Bjorken x
      xB = vertex%xbj
      
c     use momentum conservation to get components of recoil gamma
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
c DVCS+BH by C. Munoz Camacho and H. Moutarde (CEA-Saclay, IRFU/SPhN). 
c Note, sig_TGVdvcs requires all units to be passed in GeV

      if (model.eq.2) then

        if (debug(2)) then
           write(6,*)' tgvdvcs: calling '
           write(6,*)' e ',Ebeam/1.e3,E_prime/1.e3,vertex%e%theta
           write(6,*)' p ',Ep/1.e3,Pp/1.e3,thetacm,phicm
           write(6,*)' inv ',qsq/1.e6,tt/1.e6,invm/1.e3,Xb,first
        endif

        sig2 = sig_tgvdvcs(Ebeam/1.e3,E_prime/1.e3,vertex%e%theta,
     1       Ep/1.e3,Pp/1.e3,thetacm,phicm,qsq/1.e6,
     2       tt/1.e6,t_min/1.e6,uu/1.e6,u_min/1.e6,invm/1.e3,xB,first)

c Convert to  MeV units needed by SIMC
        sig2 = sig2*(1.e-9)
c        write(6,*)' tgvdvcs: sig ',sig2

        peep_dvcs = sig2
      endif

******************************************************************************
*  we keep the tradition that ntup.sigcm is d2sigma/dt/dphi_cm [ub/MeV^2/rad]
******************************************************************************

      if (model.eq.1) then

      if (debug(2)) then
         write(6,*)' dvcsgmh: calling '
         write(6,*)' e ',Ebeam/1.e3,qsq/1.e6
         write(6,*)' t ',tt/1.e6,t_min/1.e6
         write(6,*)' u ',uu/1.e6,u_min/1.e6,invm/1.e3
      endif

      sig1 = sig_dvcsgmh(Ebeam/1.e3,qsq/1.e6,thetacm,tt/1.e6,t_min/1.e6,
     *       uu/1.e6,u_min/1.e6,invm/1.e3)

      ntup%sigcm=sig1

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6

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

c      write(6,*)' dvcsgmh: sig ',peep_dvcs

      endif

      if (debug(2)) write(6,*)' peep_dvcs: end ',peep_dvcs
      return 
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_dvcsgmh(Eb_g,Q2_g,thpqcm,t_gev,tmin,u_gev,umin,
     1       W_gev)

c This routine calculates p(e,e'gamma)p DVCS cross sections, based on 
c the Kumericki and Mueller KM15 model at xB=0.36, which includes the
c 2015 CLAS and Hall A DVCS data in fixing the GPD parameters.

c NOTE: The routine only works for the beam energies, 6.6, 8.8, 10.9 GeV!

c Subroutine calculates dsigma/dt/dphi_cm, which is returned as sig_dvcsgmh
c [ub/MeV^2/rad].
c gh - 18.02.26

      implicit none
      include 'constants.inc'

      real*8 t_gev,tmin,tprime,u_gev,umin,uprime
      real*8 Eb_g,mass,Q2_g,W_gev,wfactor,m_p,thpqcm
      real*8 Q2hi,Q2lo,delQ2,Q2tmp
      real*8 Whi,Wlo,Wfac_hi,Wfac_lo
      real*8 sighi,siglo,sighiW,sigloW
      real*8 sig_phi0gmh
      real*8 sig                ! dsigma/dt
      
c use fit parameters for KM15 model at xB=0.36 (nb/GeV2)
      real*8 Q2tab(4) / 2.30,  3.00,  3.80,  5.00/
      real*8 Wtab(4)  / 2.23,  2.49,  2.76,  3.125/
      real*8 Ebtab(3) / 6.60,  8.80, 10.90/

      integer itflag/1/ ! sig fit vs t(1), vs tprime (0)

c sig_unp is parameterized in the form p1*(-t-p3)/(-t-p2)**2+p4
c fill parameter arrays in the order (x=Ee, y=Q2)

*-------
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

*-----
c GH: parameters for tprime instead of t
c Q2=2.30 for 3 beam energies
      real*8 p1p(3,4) /0.0038, 0.0050, 0.0048, 
c Q2=3.0
     1                 0.0038, 0.0063, 0.0047, 
c Q2=3.8
     2                -0.0671, 0.0074, 0.0046, 
c Q2=5.0
     3                 0.0001, 0.0001, 0.0056/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p2p(3,4) /-0.060, -0.084, -0.096, 
c Q2=3.0
     1                  0.069, -0.114, -0.121, 
c Q2=3.8
     2                 -0.175, -0.104, -0.107, 
c Q2=5.0
     3                 -0.034, -0.034, -0.110/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p3p(3,4) /0.000, -0.022, -0.034, 
c Q2=3.0
     1                -0.005, -0.040, -0.049, 
c Q2=3.8
     2                -0.118, -0.025, -0.032, 
c Q2=5.0
     3                 0.191, 0.191, -0.030/ ! use Eb=8.8 params for 6.6

c Q2=2.30
      real*8 p4p(3,4) /0.045, 0.015, 0.008, 
c Q2=3.0
     1                 0.051, 0.007, 0.002, 
c Q2=3.8
     2                 0.359, 0.006, 0.000, 
c Q2=5.0
     3                 0.062, 0.062, 0.000/ ! use Eb=8.8 params for 6.6
*-----

      integer pcount,ndat,Q2count,Q2c,Ebc,Ebcount

      ndat=4
      pcount=0
      
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
         write(6,*)' dvcsgmh: Ebcount error ',Ebcount,Eb_g
         stop
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
      
      if (thpqcm.gt.pi/2.) then    ! t-channel

         if (itflag.eq.1) then  ! sig fit vs t
           sighi = p1(Ebcount,Q2count+1)*
     1         (tprime+abs(tmin)-p3(Ebcount,Q2count+1))
     2        /(tprime+abs(tmin)-p2(Ebcount,Q2count+1))**2
     3        +p4(Ebcount,Q2count+1)
         
           siglo = p1(Ebcount,Q2count)*
     1         (tprime+abs(tmin)-p3(Ebcount,Q2count))
     2        /(tprime+abs(tmin)-p2(Ebcount,Q2count))**2
     3        +p4(Ebcount,Q2count)

        else                    ! sig fit vs tprime
           sighi = p1p(Ebcount,Q2count+1)*
     1         (tprime-p3p(Ebcount,Q2count+1))
     2        /(tprime-p2p(Ebcount,Q2count+1))**2
     3        +p4p(Ebcount,Q2count+1)
         
           siglo = p1p(Ebcount,Q2count)*
     1         (tprime-p3p(Ebcount,Q2count))
     2        /(tprime-p2p(Ebcount,Q2count))**2
     3        +p4p(Ebcount,Q2count)
        endif
         
      else                      ! u-channel
c christian weiss recommends the following change for u-channel:
c switch u-slope for t-slope, then divide by 10, since back angle peak 
c is ~10% of forward angle peak (at least for omega electroproduction)

         if (itflag.eq.1) then  ! sig fit vs t
           sighi = ( p1(Ebcount,Q2count+1)*
     1         (uprime+abs(tmin)-p3(Ebcount,Q2count+1))
     2        /(uprime+abs(tmin)-p2(Ebcount,Q2count+1))**2
     3        +p4(Ebcount,Q2count+1) )/10.
         
           siglo = ( p1(Ebcount,Q2count)*
     1        (abs(uprime)+abs(tmin)-p3(Ebcount,Q2count))
     2        /(abs(uprime)+abs(tmin)-p2(Ebcount,Q2count))**2
     3        +p4(Ebcount,Q2count) )/10.

        else                    ! sig fit vs tprime
           sighi = ( p1p(Ebcount,Q2count+1)*
     1         (uprime-p3p(Ebcount,Q2count+1))
     2        /(uprime-p2p(Ebcount,Q2count+1))**2
     3        +p4p(Ebcount,Q2count+1) )/10.
         
           siglo = ( p1p(Ebcount,Q2count)*
     1        (abs(uprime)-p3p(Ebcount,Q2count))
     2        /(abs(uprime)-p2p(Ebcount,Q2count))**2
     3        +p4p(Ebcount,Q2count) )/10.
        endif

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

         write(6,*)' dvcsgmh: Q2tmp error ',Q2tmp,Q2count

      endif
      
      sig_dvcsgmh = sig/2./pi*1.e-09 !dsig/dtdphicm in microbarns/MeV^2/rad

      if (sig_dvcsgmh .gt. 1.e-6) then
         write(6,*)' dvcsgmh: sig warning ',sig,Q2tmp,W_gev
         write(6,*)uprime,abs(tmin),p2(Ebcount,Q2count+1)
         write(6,*)sighiW,Wfac_hi,sigloW,Wfac_lo
      endif
      
      return
      end


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_tgvdvcs(Elab,E_prime,thetaEscat,
     1   Ep,Pp,thetacm,phi,qsq,t_gev,tmin,u_gev,umin,W_gev,xB,first)

c Note, all units have been passed in GeV.

c-------------------------- TgenDVCS.cxx TGVKelly.cxx --------------------------*
c                                                                               |
c Author : C. Munoz Camacho and H. Moutarde (CEA-Saclay, IRFU/SPhN).            |
c v1.4, July, 25th 2008.                                                        |
c                                                                               |
c Computation of the differential cross sections of the Bethe-Heitler and Deeply|
c Virtual Compton Scattering processes (and also of the interference of the two |
c processes).                                                                   |
c                                                                               |
c The formulas are derived from a Mathematica package written by P. Guichon     |
c (CEA-Saclay, IRFU/SPhN) and M. Vanderhaegen (W&M College), see arXiv:0808.????|
c                                                                               | 
c-------------------------------------------------------------------------------*
c adapted for SIMC by G.M. Huber, 20.02.13

      implicit none

      include 'simulate.inc'

      type(event_main):: main
      type(event):: vertex

      integer TargetPolar,BeamHeli,BeamCharge,iexact,ibh,ikelly
      data TargetPolar/0/
      data BeamCharge/-1/
      data iexact/1/ ! exact(1) or leading-twist(0) cross sect 
      data ibh/0/    ! ibh=1 returns BH cross section only
      data ikelly/0/ ! flag for Gep, Gmp parameterization default(0) JJ Kelly(1)
                     ! GH note, for u-channel, default appears better behaved
 
c--------------------------- Kinematic Variables -------------------------------*
c All quantities are in GeV or GeV^2

c Frame dependent scalars 
     
      real*8 ELab     ! Beam energy
      real*8 E_prime,thetaEscat ! scattered electron energy and angle (LAB)
      real*8 thetag   ! CM Angle between real and virtual photons
      real*8 thetacm  ! CM Angle between proton and virtual photon
      real*8 phi
      real*8 Ep,Pp    ! Lab proton energy and momentum
      real*8 qpPerp   ! Component (here x-axis) of the real photon 3-momentum 
                      ! orthogonal to the virtual photon trajectory (here z-axis) 
                      ! in the hadronic plane (here xz-plane)
      real*8 Omega    ! (Function of) the linear polarization rate
      real*8 y        ! Variable used in the computation of Omega
      real*8 e
      data e /0.30282211985966434/
      real*8 kae,kbe(3),kam,kbm(3) ! parameters for JJ Kelly form factors
      real*8 tau,taubh
      data kae /-0.24/ 
      data kbe /10.98, 12.82, 21.97/
      data kam /0.12/
      data kbm /10.97, 18.86, 6.55/
    
c Invariant scalars
         
      real*8 xB       ! Bjorken variable
      real*8 xBMin, xBMax ! Boundaries on the xB physical region
      real*8 qsq      ! Virtuality of the photon
      real*8 eps2
      real*8 t_gev    ! Mandelstam variable (square of the 4-momentum transfer)
      real*8 tmin,tmin2,tprime
      real*8 u_gev,umin,uprime
      real*8 tt
      real*8 s,s_test ! Mandelstam variable (square of the total incoming 4-momentum)
      real*8 W_gev
      real*8 Q        ! Photon virtual mass i.e. square root of Q2   
      real*8 PhaseSpace ! Differential element of cross section
      real*8 m_p      ! Proton mass

c 4-vectors defined in the CM frame
c index: 1=x, 2=y, 3=z, 4=time
    
      real*8 qCM(4)  ! Virtual photon (propagates along z-axis)
      real*8 pCM(4)  ! Incoming proton (propagates along z-axis)
      real*8 qpCM(4) ! Real photon (defines hadronic plane xz)
      real*8 ppCM(4) ! Outgoing proton

      data qCM/0.,0.,0.,0./
      data pCM/0.,0.,0.,0./
      data qpCM/0.,0.,0.,0./
      data ppCM/0.,0.,0.,0./

c--------------------------- Compton Form Factors ------------------------------*

      logical first

      integer i,j,iQ2,iQ2c,iXb,iXbc,it,itc,iGPD,line
      integer Q2_0,Q2_1,xB_0,xB_1,t_0,t_1

      real*8 V(8,13,100,90) ! index: iGPD=1,8, iQ2=1,13, iXb=1,100, it=1,90
      real*8 CFF(8)        ! interpolated CFF at Q2,Xb,t of event
      real*8 ReH,ImH,ReE,ImE,ReHT,ImHT,ReET,ImET
      equivalence (CFF(1),ReH)
      equivalence (CFF(2),ImH)
      equivalence (CFF(3),ReE)
      equivalence (CFF(4),ImE)
      equivalence (CFF(5),ReHT)
      equivalence (CFF(6),ImHT)
      equivalence (CFF(7),ReET)
      equivalence (CFF(8),ImET)

      real*8 Q2d,xBd,td
      real*8 c00,c10,c01,c11,c0,c1

c--------------------- (Combinations of) helicity amplitudes -------------------*

      real*8 Jem(4,3)   ! Helicity amplitudes of the interference process 
                        ! assuming the real photon has helicity +1.
      real*8 RMvcs(4,3) ! Real part of the helicity amplitudes of the VCS 
                        ! process assuming the real photon has helicity +1.
      real*8 IMvcs(4,3) ! Imaginary part of the helicity amplitudes of the 
                        ! VCS process assuming the real photon has helicity +1.
     
c--------------- Expansion of cross sections for harmonic analysis -------------*

c Bethe Heitler process
     
      real*8 SigmaBHPol0(4)  ! coefficients for the unpolarized cross section
      real*8 SigmaBHPolX(2)  ! coefficients for the x-polarized cross section 
      real*8 SigmaBHPolY     ! coefficient for the y-polarized cross section
      real*8 SigmaBHPolZ(2)  ! coefficients for the z-polarized cross section

      real*8 F1,F1bh         ! Dirac form factor
      real*8 F2,F2bh         ! Pauli form factor
      real*8 Ge,Gebh,Gm,Gmbh ! Sachs parametrization

      real*8 DDDC            ! denominator of BH cross section
      real*8 SigmaBHp,SigmaBHm 
      real*8 XSecBHp,XSecBHm ! BH cross sect for +,- beam helicities   
     
c Virtual Compton Scattering process
     
      real*8 SigmaVCSPol0(5) ! coefficients for the unpolarized cross section
      real*8 SigmaVCSPolX(4) ! coefficients for the x-polarized cross section     
      real*8 SigmaVCSPolY(5) ! coefficients for the y-polarized cross section     
      real*8 SigmaVCSPolZ(4) ! coefficients for the z-polarized cross section

      real*8 Ur(100)  ! Coefficients of the expansion of the interference cross 
                      ! section wrt (combinations of) helicity amplitudes

      real*8 SigmaVCSp,SigmaVCSm
      real*8 XSecVCSp,XSecVCSm ! VCS cross sect for +,- beam helicities   
    
c Interference
     
      real*8 SigmaIPol0(8)   ! coefficients for the unpolarized cross section     
      real*8 SigmaIPolX(8)   ! coefficients for the x-polarized cross section     
      real*8 SigmaIPolY(8)   ! coefficients for the x-polarized cross section
      real*8 SigmaIPolZ(8)   ! coefficients for the x-polarized cross section

      real*8 SigmaIp,SigmaIm
      real*8 XSecIp,XSecIm     ! Interf cross sect for +,- beam helicities   

c Beam Helicity Cross Section Sum and Difference

      real*8 ConvGeV2nbarn   ! Changement d'unités
      data ConvGeV2nbarn/0.389379304e+6/

      real*8 SigmaTotPlus, SigmaTotMoins
      real*8 XSecSum,XSecDif
      real*8 sig_tgendvcs    ! d4sigma/dQ2dxdtdphi [nb/GeV]
      real*8 Jac6            ! Jacobian for d6sigma required by SIMC
     
c------------------------------------------------------------------------------*

      if (debug(2)) then
         write(6,*)' tgvdvcs: entering '
         write(6,*)' e ',Elab,E_prime,thetaEscat
         write(6,*)' p ',Ep,Pp,thetacm,phi
         write(6,*)' inv ',qsq,t_gev,W_gev,Xb,first
      endif

c------------------- Read in Compton Form Factors (CFF) ------------------------*
c Note, the units in this file are all in GeV, GeV^2

      if (first) then
         write(6,*,err=1000)' tgvdvcs: reading CFFoutput_LO.dat'
         open(unit=31,file='CFFoutput_LO.dat',status='old',action='read')
         line=0
         do iQ2=1,13
           do ixB=1,100
             do it=1,90
               read(31,*,err=1000)Q2d,xBd,td,(V(i,iQ2,iXb,it), i=1,8)
               line=line+1
c               write(6,100)line,Q2d,xBd,td,V(8,iQ2,iXb,it)
c 100           format(i7,4f10.4)

c check that element was read in correctly
               itc  = ifix( -(td-0.)*89/3. +0.5)+1
               if (itc.ne.it) write(6,101)line,itc,it,td,V(8,iQ2,iXb,it)
 101           format(' tgvdvcs: td  err ',i7,2i5,2f10.4)

             enddo
             ixBc = ifix( (xBd-0.2)*99/0.7 +0.5)+1
             if (ixBc.ne.ixB) write(6,102)line,ixBc,ixB,xBd,V(8,iQ2,iXb,it-1)
 102         format(' tgvdvcs: xBd err ',i7,2i5,2f10.4)

           enddo
           iQ2c = ifix(Q2d-1. +0.5)+1
           if (iQ2c.ne.iQ2) write(6,103)line,iQ2c,iQ2,Q2d,V(8,iQ2,iXb-1,it-1)
 103       format(' tgvdvcs: Q2d err ',i7,2i5,2f10.4)

         enddo
         close(31)
      endif

c---------------------------- MakeKinematics -----------------------------------*
c Global variables xB, Q2, t come from main SIMC
c Compute 4-vector components of the particles involved in CM frame 
c    virtual photon's 4-momentum qCM,
c    real photon's 4-momentum qpCM, 
c    incoming proton's 4-momentum pCM,  
c    outgoing proton's 4-momentum ppCM,   
c    thetag angle between the trajectories of real and virtual photons.  
c It also evaluates the parameter Omega and the invariant phase space PhaseSpace.

      m_p = Mp/1.e3    ! value from constants.inc
      Q   = sqrt(qsq)

      s   = W_gev**2
c     s   = TMath::Power(M,2) - TMath::Power(Q,2) + TMath::Power(Q,2)/xB;
      s_test= m_p**2 - qsq + qsq/xB
      if (abs(s-s_test).gt.1.e-5) then
         write(6,*)' tgvdvcs: s error ',s,s_test,m_p,qsq,xB
         stop
      endif
      tprime = abs(t_gev-tmin)
      uprime = abs(u_gev-umin)

c Test : realistic kinematic configuration ?
c The value of xBMin comes from the requirement of omega to be real, and 
c the value of xBMax expresses the fact that s >= 0. 
          
      xBMin = 2.*ELab*qsq/(m_p*(4*(ELab**2)-qsq))
      xBMax = qsq/(qsq-(m_p**2))

      if ( xB.lt.xBMin .or. xB.gt.xBMax ) then
         write(6,*)' tgv_dvcs: Unrealistic kinematic configuration '
         write(6,*)'  -- xB isnt in the physical region !'
         write(6,*)' xB = ',xB
         write(6,*)' xBMin = ',xBMin 
         write(6,*)' xBMax = ',xBMax
         stop
      endif

c Timelike coordinate
c NOTE: 4th coordinate, not zero (since Fortran array numbering starts at 1)

      qCM(4)  = ((qsq*(1-2*xB))/(2.*sqrt((m_p**2) +(qsq*(1-xB))/xB)*xB))
      qpCM(4) = (-(qsq*(-1+xB))/(2.*sqrt((m_p**2) +(qsq*(1-xB))/xB)*xB))
      pCM(4)  = ((qsq +2*xB*m_p**2)/(2.*sqrt(m_p**2+(qsq*(1-xB))/xB)*xB))
      ppCM(4) = ((2*m_p**2+(qsq*(1-xB))/xB)/(2.*sqrt(m_p**2+(qsq*(1-xB))/xB)))
     
c Spacelike coordinates
c 1=x, 2=y, 3=z
     
      qpCM(1) = sqrt((-((m_p**2)*(t_gev**2)*(xB**2)) + 
     1    qsq*t_gev*xB*(t_gev*(-1 + xB) - 2*(m_p**2)*xB) + 
     2    (qsq**2)*(t_gev*(-1 + xB) - (m_p**2)*(xB**2)))/
     3    ((qsq**2) + 4*(m_p**2)*qsq*(xB**2)))
      if (isnan(qpCM(1))) qpCM(1)=0.

      qCM(3) = sqrt((qsq**2 + 4*(m_p**2)*qsq*(xB**2))/
     1    (xB*(qsq*(1 - xB) + (m_p**2)*xB)))/2.

      qpCM(3) = ((qsq**2)*(-1 + xB) - 2*(m_p**2)*t_gev*(xB**2) + 
     1    2*qsq*xB*(t_gev*(-1 + xB) - (m_p**2)*xB))/
     2    (2.*xB*(qsq*(-1 + xB) - (m_p**2)*xB)*
     3    sqrt(((qsq**2) + 4*(m_p**2)*qsq*(xB**2))/
     4    (xB*(qsq*(1 - xB) + (m_p**2)*xB))))

      pCM(3)  = (-qCM(3))
      ppCM(3) = (-qpCM(3))
      qpPerp  = qpCM(1)

      thetag = ACos(qpCM(3)/qpCM(4))   ! CM angle of DVCS photon, not proton!
      if (abs(thetacm-(pi-thetag)).gt.1.e-6) then
         write(6,*)' tgvdvcs: theta* err ',pi-thetag,thetacm
      endif
          
c Omega
          
      y = (-qsq + 4*ELab*m_p*xB)/sqrt(qsq**2 + 4*(m_p**2)*qsq*(xB**2))

      if ( y .lt. 1. ) then
         write(6,*) ' tgv_dvcs: ELab ',ELab,' GeV is too small.'
         write(6,*)'  CosH[ome]= ',y 
         stop
      endif

      Omega = Log(y + sqrt(-1 + y)*sqrt(1 + y))
          
c PhaseSpace
          
      PhaseSpace = ((e**6)*qsq)/
     1    ( 4096.*(pi**5)* sqrt( 4*(m_p**2)*qsq +(qsq**2)/(xB**2) )
     2    *(xB**2)* ( qsq/(4.*xB) + 
     3    ( sqrt( 4*(m_p**2)*qsq +(qsq**2/xB**2)) *CosH(Omega)) /4.)**2)
          
      if (debug(2)) then
         write(6,*)' Q ',Q
         write(6,*)' qpPerp ',qpPerp
         write(6,*)' s ',s,s_test
         write(6,*)' t ',t_gev,tprime
         write(6,*)' u ',u_gev,uprime
         write(6,*)' theta* ',thetag,pi-thetacm,phi
         write(6,*)' xB ',xB
         write(6,*)' pCM  ',pCM(1),pCM(2),pCM(3),pCM(4)
         write(6,*)' ppCM ',ppCM(1),ppCM(2),ppCM(3),ppCM(4)
         write(6,*)' qCM  ',qCM(1),qCM(2),qCM(3),qCM(4)
         write(6,*)' qpCM ',qpCM(1),qpCM(2),qpCM(3),qpCM(4)
      endif 

c---------------- Interpolate CFF to Q2, xB, t of the event --------------------*

      eps2 = 4.*(xB**2)*(m_p**2)/qsq   ! epsilon **2

      tmin2 = -qsq*(2.*(1.-xB)*(1-Sqrt(1+eps2))+eps2)/(4.*xB*(1.-xB)+eps2)
      if (abs(tmin-tmin2).gt.1.e-6) then
         write(6,*)' tgvdvcs: tmin error ',tmin,tmin2
         stop
      endif

c regular t-channel model
      if (thetacm.gt.pi/2. .and. ibh.eq.0) then
        tt=t_gev
        if (first) write(6,*)' dvcs: t-channel '
        first=.FALSE.

        if (qsq.lt.1. .or. qsq.gt.13. .or. xB.lt.0.1 .or. xB.gt.0.9 
     1     .or. tt.lt.-3.0) then
          write(6,*)' Kinematics (Q2,xB,t,tmin) out of range for xsec'
          write(6,*)' Q2 ',qsq
          write(6,*)' xB ',xB
          write(6,*)' t  ',tt,tmin,tprime,thetacm
          goto 990
       endif
      endif

c ad-hoc u-channel model modification for VCS and Interference
c Christian Weiss recommends the following change for u-channel:
c switch u-slope for t-slope, then divide by 10, since back angle peak 
c is ~10% of forward angle peak (at least for omega electroproduction)
      if (thetacm.lt.pi/2. .and. ibh.eq.0) then
        tt=-1.*(uprime+abs(tmin))
        if (first) write(6,*)' dvcs: t-channel '
c        write(6,*)' tt: u-chan ',tt,uprime,abs(tmin)

        if (qsq.lt.1. .or. qsq.gt.13. .or. xB.lt.0.1 .or. xB.gt.0.9 
     1     .or. tt.lt.-3.0) then
          write(6,*)' Kinematics (Q2,xB,u,umin) out of range for xsec'
          write(6,*)' Q2 ',qsq
          write(6,*)' xB ',xB
          write(6,*)' u  ',tt,umin,uprime,thetacm
          goto 990
       endif
      endif

      if (ibh.eq.1) goto 300    ! BH cross sect only

c add +1 to all indices, since fortran starts at 1, but c++ starts at 0.
      Q2_0 = ifix(qsq-1.)+1
      Q2_1 = ifix(qsq-1.)+1+1
      xB_0 = ifix( (xB-0.2)*99/0.7 )+1
      xb_1 = ifix( (xB-0.2)*99/0.7 )+1+1
      t_0  = ifix( -(tt-0.)*89/3. )+1
      t_1  = ifix( -(tt-0.)*89/3. )+1+1

      Q2d = (qsq-(1+Q2_0))/((1+Q2_1)-(1+Q2_0))
      xBd = (xB-(0.1+0.8*xB_0/79.))/
     1    ((0.1+0.8*xB_1/79.)-(0.1+0.8*xB_0/79.))
      td  = (tt-(0.-3.*t_0/89.))/((0.-3.*t_1/89.)-(0.-3.*t_0/89.))

      do iGPD=1,8
         c00=V(iGPD,Q2_0,xb_0,t_0)*(1-Q2d)+V(iGPD,Q2_1,xb_0,t_0)*Q2d
         c10=V(iGPD,Q2_0,xb_1,t_0)*(1-Q2d)+V(iGPD,Q2_1,xb_1,t_0)*Q2d
         c01=V(iGPD,Q2_0,xb_0,t_1)*(1-Q2d)+V(iGPD,Q2_1,xb_0,t_1)*Q2d
         c11=V(iGPD,Q2_0,xb_1,t_1)*(1-Q2d)+V(iGPD,Q2_1,xb_1,t_1)*Q2d
    
         c0=c00*(1-xBd)+c10*xBd
         c1=c01*(1-xBd)+c11*xBd
    
         CFF(iGPD)=c0*(1-td)+c1*td
      enddo

c------- MakeVCSHelicityAmplitudes(ReH,ImH,ReE,ImE,ReHT,ImHT,ReET,ImET) --------*
c Compute the hadronic helicity amplitudes RMvcs and IMvcs at leading twist.    |
c The formula depend on the real and imaginary part of the Compton Form Factor: |
c    ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET.                                |
c-------------------------------------------------------------------------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

      RMvcs(1,1)=0
      RMvcs(1,2)=0
      RMvcs(1,3)=-((-4*ReH*(-1+xB) -ReE*(xB**2))/(sqrt(1-xB)*(-2+xB)))
      RMvcs(2,1)=0
      RMvcs(2,2)=0
      RMvcs(2,3)=-((qpPerp*ReET*xB)/(sqrt(1-xB)*(-2*m_p + m_p*xB)))
      RMvcs(3,1)=0
      RMvcs(3,2)=0
      RMvcs(3,3)=(qpPerp*ReE)/(m_p*sqrt(1-xB))
      RMvcs(4,1)=0
      RMvcs(4,2)=0
      RMvcs(4,3)=-((-4*ReHT*(-1+xB)-ReET*(xB**2))/(sqrt(1-xB)*(-2+xB)))

      IMvcs(1,1)=0
      IMvcs(1,2)=0
      IMvcs(1,3)=-((-4*ImH*(-1+xB) -ImE*(xB**2))/(sqrt(1-xB)*(-2+xB)))
      IMvcs(2,1)=0
      IMvcs(2,2)=0
      IMvcs(2,3)=-((ImET*qpPerp*xB)/(sqrt(1-xB)*(-2*m_p + m_p*xB)))
      IMvcs(3,1)=0
      IMvcs(3,2)=0
      IMvcs(3,3)=(ImE*qpPerp)/(m_p*sqrt(1-xB))
      IMvcs(4,1)=0
      IMvcs(4,2)=0
      IMvcs(4,3)=-((-4*ImHT*(-1+xB)-ImET*(xB**2))/(sqrt(1-xB)*(-2+xB)))      
    
      if (debug(2)) then
        do i=1,4
          write(6,*)' IMvcs ',i,(IMvcs(i,j),j=1,3)
          write(6,*)' RMvcs ',i,(RMvcs(i,j),j=1,3)
        enddo
      endif

 300  continue

c----------------------------- MakeBHCrossSections() ---------------------------*
c Computes all the stuff to evaluate the cross section assuming the hadronic    |
c helicity amplitudes are given, i.e. initializes :                             |
c   - the helicity amplitudes Jem,                                              |
c   - the SigmaBHPol's.                                                         |
c-------------------------------------------------------------------------------*
c BH is exact calculation, needs actual t_gev and Ge,Bm calc with this value
c VCS+I need tt due to ad-hoc u-channel extension

c proton form factors

      if (ikelly.eq.1) then  ! JJ Kelly PRC 70, 068202 (2004)

         taubh = -t_gev/(4.*(m_p**2))
         tau   = -tt/(4.*(m_p**2))

         Gebh = kae*taubh/
     1     (kbe(1)*taubh+kbe(2)*(taubh**2)+kbe(3)*(taubh**3))
         Ge   = kae*tau/(kbe(1)*tau+kbe(2)*(tau**2)+kbe(3)*(tau**3))

         Gmbh = 2.79285*kam*taubh/
     1     (kbm(1)*taubh+kbm(2)*(taubh**2)+kbm(3)*(taubh**3))
         Gm   = 2.79285*kam*tau/
     1     (kbm(1)*tau+kbm(2)*(tau**2)+kbm(3)*(tau**3))

         F1bh = (-1.*t_gev/4./(m_p**2)*Gmbh+Gebh)/(1.-t_gev/4./(m_p**2))
         F1   = (-1.*tt/4./(m_p**2)*Gm+Ge)/(1.-tt/4./(m_p**2))

         F2bh = ( Gmbh-Gebh)/(1.-t_gev/4./(m_p**2))
         F2   = ( Gm-Ge)/(1.-tt/4./(m_p**2))

      else

        F1bh = (4.*(m_p**2) - 2.79285*t_gev) /
     1    ( ((1. - 1.4084507042253522*t_gev)**2)*(4.*(m_p**2)-t_gev) )
        F1   = (4.*(m_p**2) - 2.79285*tt) /
     1    ( ((1. - 1.4084507042253522*tt)**2)*(4.*(m_p**2)-tt) )

        F2bh = (7.1714*(m_p**2)) / 
     1    ( ((1 -1.4084507042253522*t_gev)**2)*(4*(m_p**2)-t_gev) )
        F2   = (7.1714*(m_p**2)) / 
     1    ( ((1 -1.4084507042253522*tt)**2)*(4*(m_p**2)-tt) )

        Gmbh = 2.79285/( (1 - 1.4084507042253522*t_gev)**2 )
        Gm   = 2.79285/( (1 - 1.4084507042253522*tt)**2 )

        Gebh = 1./( (1 - 1.4084507042253522*t_gev)**2)
        Ge   = 1./( (1 - 1.4084507042253522*tt)**2)

      endif

c-------------- Helicity amplitudes of the interference process ----------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

c iexact=1:
c No approximations apart those in the computation of the hadronic helicity 
c amplitudes in terms of the Compton Form Factors.

      if (iexact.eq.1) then

        Jem(1,1)=(-2*sqrt(2.)*qCM(3)*qpCM(1)*F2*
     1   (sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) - 
     2   sqrt(pCM(4) +m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*cos(thetag/2.))/
     3   (m_p*sqrt((qsq**2) - 4*s*tt)) - 
     4   (2*sqrt(2.)*(qCM(3) + qpCM(4))*Gm*
     5   (sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) + 
     6   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*sqrt(s)*sin(thetag/2.))/
     7   sqrt((qsq**2) - 4*s*tt)

        Jem(1,2)=-((F2*(sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) - 
     1   sqrt(pCM(4) + m_p)*sqrt(ppCM(4) + m_p))*
     2   (2*(m_p**2) + qsq + 2*s)*sqrt(-(s*tt))*cos(thetag/2.))/
     3   (m_p*sqrt(s)*sqrt((qsq**2) - 4*s*tt))) - 
     4   (8*Gm*(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) + 
     5   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*s*sqrt(-tt)*
     6   ((qCM(3) - qpCM(3))*cos(thetag/2.) - 
     7   qpCM(1)*sin(thetag/2.)))/(qsq*sqrt((qsq**2) - 4*s*tt))

        Jem(1,3)=(2*sqrt(2.)*qCM(3)*qpCM(1)*F2*
     1   (sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) - 
     2   sqrt(pCM(4) + m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*cos(thetag/2.))/
     3   (m_p*sqrt((qsq**2) - 4*s*tt)) + 
     4   (2*sqrt(2.)*(qCM(3) + qpCM(4))*Gm*
     5   (sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) + 
     6   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*sqrt(s)*sin(thetag/2.))/
     7   sqrt((qsq**2) - 4*s*tt)

        Jem(2,1)=sqrt(2.)*Gm*(-(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) +m_p)) + 
     1   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*cos(thetag/2.)
        Jem(2,2)=0
        Jem(2,3)=sqrt(2.)*Gm*(-(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) +m_p)) + 
     1   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*cos(thetag/2.)

        Jem(3,1)=(2*sqrt(2.)*(qCM(3) - qpCM(4))*Gm*
     1   (sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) - 
     2   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*cos(thetag/2.))/
     3   sqrt((qsq**2) - 4*s*tt) - (2*sqrt(2.)*qCM(3)*qpCM(1)*F2*
     4   (sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) + 
     5   sqrt(pCM(4) + m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*sin(thetag/2.))/
     6   (m_p*sqrt((qsq**2) - 4*s*tt))

        Jem(3,2)=-((F2*(sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) + 
     1   sqrt(pCM(4) + m_p)*sqrt(ppCM(4) + m_p))*
     2   (2*(m_p**2) + qsq + 2*s)*sqrt(-(s*tt))*sin(thetag/2.))/
     3   (m_p*sqrt(s)*sqrt((qsq**2) - 4*s*tt))) - 
     4   (8*Gm*(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) - 
     5   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*s*sqrt(-tt)*
     6   (qpCM(1)*cos(thetag/2.) + (qCM(3) - qpCM(3))*sin(thetag/2.)))/
     7   (qsq*sqrt((qsq**2) - 4*s*tt))

        Jem(3,3)=(-2*sqrt(2.)*(qCM(3) - qpCM(4))*Gm*
     1   (sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) - 
     2   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*cos(thetag/2.))/
     3   sqrt((qsq**2) - 4*s*tt) + (2*sqrt(2.)*qCM(3)*qpCM(1)*F2*
     4   (sqrt(pCM(4) - m_p)*sqrt(ppCM(4) - m_p) + 
     5   sqrt(pCM(4) + m_p)*sqrt(ppCM(4) +m_p))*sqrt(s)*sin(thetag/2.))/
     6   (m_p*sqrt((qsq**2) - 4*s*tt))

        Jem(4,1)=sqrt(2.)*Gm*(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) + 
     1   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*sin(thetag/2.)
        Jem(4,2)=0
        Jem(4,3)=sqrt(2.)*Gm*(sqrt(ppCM(4) - m_p)*sqrt(pCM(4) + m_p) + 
     1   sqrt(pCM(4) - m_p)*sqrt(ppCM(4) + m_p))*sin(thetag/2.)

      endif

c iexact=0:
c All quantities are evaluated at leading order in the 1/Q expansion 
c (see Belitsky, Mueller and Kirchner, Nucl.Phys. B629 (2002) 323-392,
c ArXiv:hep-ph/0112108v2)

      if (iexact.eq.0) then

       Jem(1,1)=-((sqrt(2.)*(F2 - Gm)*qpPerp*(-2 + xB))/(sqrt(1-xB)*xB))
       Jem(1,2)=-(((F2*((-2 + xB)**2) + 4*Gm*(-1 + xB))*
     1   sqrt( (qpPerp**2) + (m_p**2)*(xB**2))) / ((-1 + xB)*xB))
       Jem(1,3)=(sqrt(2.)*(F2 - Gm)*qpPerp*(-2 + xB))/(sqrt(1-xB)*xB)

       Jem(2,1)=(sqrt(2.)*Gm*m_p*xB)/sqrt(1-xB)
       Jem(2,2)=0
       Jem(2,3)=(sqrt(2.)*Gm*m_p*xB)/sqrt(1-xB)

       Jem(3,1)=-((sqrt(2.)*(F2*(qpPerp**2) + Gm*(m_p**2)*(xB**2)))/
     1     (m_p*sqrt(1 - xB)*xB))
       Jem(3,2)=-((F2*qpPerp*(-2 + xB)*
     1     sqrt((qpPerp**2) + (m_p**2)*(xB**2)))/(m_p*(-1 + xB)*xB))
       Jem(3,3)=(sqrt(2.)*(F2*(qpPerp**2) 
     1      + Gm*(m_p**2)*(xB**2)))/(m_p*sqrt(1-xB)*xB)

       Jem(4,1)=(sqrt(2.)*Gm*qpPerp)/sqrt(1-xB)
       Jem(4,2)=0
       Jem(4,3)=(sqrt(2.)*Gm*qpPerp)/sqrt(1-xB)

      endif
         
c----------------------------- BH cross sections -------------------------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

c BH is exact calculation, needs actual t_gev and Ge,Bm calc with this value

        SigmaBHPol0(1)=( -32*(Gebh**2)*(m_p**2)*(3*(m_p**8)*t_gev + 
     1    3*(m_p**6)*((qsq**2) + qsq*t_gev - 4*s*t_gev) + 
     2    (m_p**2)*(3*(qsq**2)*((qsq + s)**2) - 
     3    (qsq + s)*((qsq**2) + 3*qsq*s + 12*(s**2))*t_gev + 
     4    2*((qsq**2) + qsq*s - 3*(s**2))*(t_gev**2)) + 
     5    ((qsq + s)**2)*t_gev*(3*s*(s + t_gev) + qsq*(3*s + t_gev)) + 
     6    (m_p**4)*(4*(qsq**3) + qsq*t_gev*(3*s + t_gev) - 
     7    (qsq**2)*(6*s + t_gev) + 3*s*t_gev*(6*s + t_gev))) + 
     8     4*(Gmbh**2)*t_gev*(6*(m_p**8)*t_gev + t_gev*(3*((qsq + s)**2)*
     9    ((qsq**2) + 2*qsq*s + 2*(s**2)) + 
     1    2*s*(qsq + s)*(2*qsq + 3*s)*t_gev + 
     2    (3*(qsq**2) + 4*qsq*s + 3*(s**2))*(t_gev**2)) - 
     3    2*(m_p**6)*(3*(qsq**2) - 11*qsq*t_gev + 6*t_gev*(2*s+t_gev)) + 
     4    (m_p**4)*(-8*(qsq**3) - 26*qsq*t_gev*(s + t_gev) + 
     5    (qsq**2)*(12*s + 25*t_gev) + 
     6    3*t_gev*(12*(s**2) + 10*s*t_gev + (t_gev**2))) - 
     7    2*(m_p**2)*(3*(qsq**4) + (qsq**3)*(6*s - 5*t_gev) + 
     8    3*s*t_gev*((2*s + t_gev)**2) 
     9    + qsq*t_gev*(7*(s**2) + 2*s*t_gev - 3*(t_gev**2)) + 
     1    (qsq**2)*(3*(s**2) - 5*s*t_gev + 7*(t_gev**2)))))/
     2    (((m_p**4) + 2*(m_p**2)*(qsq - s) + ((qsq + s)**2))*
     3    (4*(m_p**2) - t_gev))

      SigmaBHPol0(2)=( -32*(Gebh**2)*(m_p*2)*((m_p*8)*t_gev + 
     1    (m_p**6)*((qsq**2) + qsq*t_gev - 4*s*t_gev) + 
     2    (m_p*2)*((qsq**2)*((qsq + s)**2) + 
     3    (qsq - s)*(qsq + s)*(5*qsq + 4*s)*t_gev - 
     4    2*((qsq**2) + qsq*s + (s**2))*(t_gev**2)) + 
     5    ((qsq + s)**2)*t_gev*(qsq*(s - t_gev) + s*(s + t_gev)) + 
     6    (m_p**4)*(-4*(qsq**3) + qsq*(s - t_gev)*t_gev + 
     7    s*t_gev*(6*s + t_gev) + (qsq**2)*(-2*s + 5*t_gev))) + 
     8    4*(Gmbh**2)*t_gev*(2*(m_p**8)*t_gev + t_gev*((qsq + s)**2)*
     9    ((qsq**2) + 2*qsq*s + 2*(s**2)) + 
     1    2*s*(-2*(qsq**2) - qsq*s + (s**2))*t_gev + 
     2    ((qsq**2) - 4*qsq*s + (s**2))*(t_gev**2)) - 
     3    2*(m_p**6)*((qsq**2) - 9*qsq*t_gev + 2*t_gev*(2*s + t_gev)) + 
     4    (m_p**4)*(8*(qsq**3) - 2*qsq*t_gev*(15*s + 7*t_gev) + 
     5    (qsq**2)*(4*s + 19*t_gev) + t_gev*(12*(s**2) + 10*s*t_gev + 
     6    (t_gev**2))) -2*(m_p**2)*( (qsq**4) + (qsq**3)*(2*s + t_gev) + 
     7    s*t_gev*((2*s + t_gev)**2) - qsq*t_gev*(3*(s**2) +10*s*t_gev + 
     8    (t_gev**2)) + (qsq**2)*((s**2) - 7*s*t_gev + 5*(t_gev**2))))/ 
     9    (((m_p**4) + 2*(m_p**2)*(qsq - s) + ((qsq + s)**2))*
     1    (4*(m_p**2) - t_gev)) 

      SigmaBHPol0(3)=(-64*(Gebh**2)*(m_p**2)*Q*qpPerp*((m_p**2) -qsq -s)*
     1    (2*(m_p**2)*qsq - ((m_p**2) + qsq + s)*t_gev) + 
     2    16*(Gmbh**2)*Q*qpPerp*t_gev*((m_p**4)*(-2*qsq + 3*t_gev) + 
     3    t_gev*(qsq*(s - t_gev) + s*(s + t_gev)) + 
     4    (m_p**2)*(2*qsq**2 - t_gev*(4*s + t_gev) + 
     5    qsq*(2*s + 5*t_gev))))/ (sqrt((m_p**4) + 2*(m_p**2)*(qsq -s) + 
     6    ((qsq + s)**2))*(4*(m_p**2) - t_gev))

      SigmaBHPol0(4)=(-64*(Gebh**2)*(m_p**4)*qsq*
     1    ((m_p**4)*t_gev + s*t_gev*(qsq + s + t_gev) + 
     2    (m_p**2)*(qsq**2 + qsq*t_gev - 2*s*t_gev)) - 
     3    8*(Gmbh**2)*qsq*(2*(m_p**2) - t_gev)*t_gev*
     4    ((m_p**4)*t_gev + s*t_gev*(qsq + s + t_gev) + 
     5    (m_p**2)*(qsq**2 + qsq*t_gev - 2*s*t_gev)))/
     6    (((m_p**4) + 2*(m_p**2)*(qsq - s) + ((qsq + s)**2))*
     7    (4*(m_p**2) - t_gev))

      SigmaBHPolX(1)=(64*Gebh*Gmbh*m_p*qpPerp*(-(qsq*(qsq + s)) + 
     1    (m_p**2)*(qsq - t_gev) + s*t_gev)*
     2    (2*(m_p**2)*qsq - ((m_p**2) + qsq + s)*t_gev) + 
     3    32*(Gmbh**2)*m_p*qpPerp*t_gev*(-(qsq**3) - 3*qsq**2*s - 
     4    2*(m_p**4)*(qsq - t_gev) + s*t_gev*(2*s + t_gev) - 
     5    qsq*(2*(s**2) + (t_gev**2)) - 
     6    (m_p**2)*(qsq**2 - 4*qsq*(s + t_gev) + t_gev*(4*s + t_gev))))/
     7    (sqrt((m_p**4) + 2*(m_p**2)*(qsq - s) + 
     8    ((qsq + s)**2))*(4*(m_p**2) - t_gev))

      SigmaBHPolX(2)=(32*Gebh*Gmbh*m_p*Q*(-2*(m_p**2)*qsq + 
     1    ((m_p**2) + qsq + s)*t_gev) *(4*(m_p**2)*qsq**2 + 
     2    (2*(m_p**2) + qsq - 2*s)*((m_p**2) - qsq - s)*t_gev + 
     3    ((m_p**2) + qsq + 3*s)*(t_gev**2)) + 
     4    64*(Gmbh**2)*m_p*Q*((m_p**2) - s - t_gev)*t_gev*
     5    ((m_p**4)*t_gev + s*t_gev*(qsq + s + t_gev) + 
     6    (m_p**2)*(qsq**2 + qsq*t_gev - 2*s*t_gev)))/
     7    (((m_p**4) + 2*(m_p**2)*(qsq - s) + ((qsq + s)**2))*
     8    (4*(m_p**2) - t_gev))

      SigmaBHPolY=32*Gebh*Gmbh*m_p*Q*(qsq - t_gev)*t_gev

      SigmaBHPolZ(1)=(-128*Gebh*Gmbh*(m_p**2)*(-(qsq*(qsq + s)) + 
     1    (m_p**2)*(qsq - t_gev) + s*t_gev)* ((m_p**4)*t_gev + 
     2     s*t_gev*(qsq + s + t_gev) + 
     3    (m_p**2)*(qsq**2 + qsq*t_gev - 2*s*t_gev)) + 
     4    16*(Gmbh**2)*t_gev*(2*(m_p**2)*qsq - ((m_p**2) + qsq +s)*t_gev)*
     5    ((qsq**3) + 3*qsq**2*s + 2*(m_p**4)*(qsq - t_gev) - 
     6    s*t_gev*(2*s + t_gev) + qsq*(2*(s**2) + (t_gev**2)) + 
     7    (m_p**2)*(qsq**2 - 4*qsq*(s + t_gev) + t_gev*(4*s + t_gev))))/
     8    (((m_p**4) + 2*(m_p**2)*(qsq - s) + ((qsq + s)**2))*
     9    (4*(m_p**2) - t_gev))

      SigmaBHPolZ(2)=(64*Gebh*Gmbh*(m_p**2)*Q*qpPerp*(-4*(m_p**2)*qsq**2 - 
     1    (2*(m_p**2) + qsq - 2*s)*((m_p**2) - qsq - s)*t_gev - 
     2    ((m_p**2) + qsq + 3*s)*(t_gev**2)) + 
     3    32*(Gmbh**2)*Q*qpPerp*t_gev*(-(m_p**2) + s + t_gev)*
     4    ((qsq + s)*t_gev + (m_p**2)*(-2*qsq + t_gev)))/
     5    (sqrt((m_p**4) + 2*(m_p**2)*(qsq - s) + 
     6    ((qsq + s)**2))*(4*(m_p**2) - t_gev))
    
      if (debug(2)) then 
         write(6,*)' F1 ',F1bh,F1
         write(6,*)' F2 ',F2bh,F2
         write(6,*)' Ge ',Gebh,Ge
         write(6,*)' Gm ',Gmbh,Gm

         do i=1,4
            write(6,*)' Jem ',i,(Jem(i,j),j=1,3)
         enddo          
          
         write(6,*)' SigmaBHPol0 ',(SigmaBHPol0(i),i=1,4)
         write(6,*)' SigmaBHPolX ',(SigmaBHPolX(i),i=1,2)
         write(6,*)' SigmaBHPolY ',SigmaBHPolY
         write(6,*)' SigmaBHPolZ ',(SigmaBHPolZ(i),i=1,2)
      endif

      if (ibh.eq.1) goto 500  ! BH cross sect only
     
c---------------------- MakeLeadingVCSAndInterfCrossSections -------------------*
c Computes all the stuff to evaluate the cross section assuming the hadronic    |
c helicity amplitudes are given, i.e. initializes :                             |
c   - the helicity amplitudes Mvcs,                                             |
c   - the expansion coefficients Ur,                                            |
c   - the SigmaVCSPol's and SigmaIPol's.                                        |
c All quantities are evaluated at leading order in the 1/Q expansion            |
c (see Belitsky, Mueller and Kirchner, Nucl.Phys. B629 (2002) 323-392,          |
c ArXiv:hep-ph/0112108v2)                                                       |
c-------------------------------------------------------------------------------*
 
c---------- Harmonic expansion coefficients of the VCS cross section -----------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

c iexact=1:
c No approximations apart those in the computation of the hadronic helicity 
c amplitudes in terms of the Compton Form Factors.

      if (iexact.eq.1) then

       Ur(1)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + 4*qCM(4)*(4*qCM(4) - 3*qCM(3))*
     2    (qpPerp**2) + (qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) + qCM(3) + 2*qpCM(4)) + 
     4    16*(qpPerp**2) + tt)))/
     5    (64.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(2)=(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) - 4*qCM(4)*qCM(3)*(qpPerp**2) + 
     2    (qCM(4) - qCM(3))*(qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))*(qCM(4) - qCM(3) + 2*qpCM(4)) + tt)))/
     4    (64.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(3)=(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + 4*qCM(4)*qCM(3)*(qpPerp**2) + 
     2    (qCM(4) + qCM(3))*(qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) + qCM(3) +2*qpCM(4)) + tt)))/
     4    (64.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(4)=(qpPerp*((qsq**2) + qsq*(12*(qCM(4)**2) -16*qCM(4)*qCM(3) - 
     1    22*qCM(4)*qpCM(4) + 26*qCM(3)*qpCM(4) - 7*tt) + 
     2    4*qCM(4)*(qCM(4) - 4*qCM(3))*tt)*
     3    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))
     4    /(128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(5)=(3*qpPerp*(2*(qCM(4) + qCM(3))*qpCM(4) + qsq + tt)*
     1    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))
     2    /(128.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(6)=(qpPerp*(5*(qsq**2) + qsq*
     1    (2*(6*(qCM(4)**2) - 7*qCM(4)*qpCM(4) + 
     2    qCM(3)*qpCM(4)) - 3*tt) + 4*(qCM(4)**2)*tt)*
     3    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))
     4    /(128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(7)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + qsq*((qCM(4) - qpCM(4))*
     2    (qCM(4) - qCM(3) + 2*qpCM(4)) + 2*(qpPerp**2) + tt) + 
     3    (qCM(4) - qCM(3))* 
     4    (2*qCM(4)*(qpPerp**2) + (qCM(4) + qpCM(4))*tt)))/
     5    (32.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(8)=-(qpPerp*(qCM(4)*(3*qCM(4) - qCM(3) - 2*qpCM(4))*
     1    qsq + 2*(qsq**2) + qCM(4)*(qCM(4) - qCM(3))*tt)*
     2    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)
     3    )/(16.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)* sqrt((qsq**2) - 4*s*tt))

       Ur(9)=(qpPerp*sqrt(s)*((4*(qCM(4)**2) - 
     1    3*qCM(4)*(qCM(3) - qpCM(4)) - 4*qCM(3)*qpCM(4))*tt + 
     2    qsq*(-3*(qCM(4) + qCM(3) - 2*qpCM(4))* 
     3    (qCM(4) - qpCM(4)) + 7*tt)))/
     4    (4.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(10)=-(qpPerp*sqrt(s)*(qCM(4)*(qCM(3) + qpCM(4))*tt + 
     1    qsq*(-((qCM(4) - qCM(3) - 2*qpCM(4))*
     2    (qCM(4) - qpCM(4))) + tt)))/
     3    (4.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(11)=(qpPerp*sqrt(s)*(qCM(4)*(-qCM(3) + qpCM(4))*tt + 
     1    qsq*(-((qCM(4) + qCM(3) - 2*qpCM(4))*
     2    (qCM(4) - qpCM(4))) + tt)))/
     3    (4.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(12)=(sqrt(s)*(3*(qsq**3) + (qsq**2)*
     1    (3*(qCM(4)**2) - 7*qCM(4)*qCM(3) + 35*qCM(4)*qpCM(4) + 
     2    3*qCM(3)*qpCM(4) - 22*(qpCM(4)**2) + 12*(qpPerp**2) + 3*tt) + 
     3    8*qCM(4)*(4*qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(8*(4*qCM(4) - qCM(3))*qpCM(4)* 
     6    ((qCM(4)**2) - (qpCM(4)**2)) + 
     7    2*(22*(qCM(4)**2) + 7*qCM(3)*qpCM(4) + 
     8    qCM(4)*(-15*qCM(3) +16*qpCM(4)))*(qpPerp**2) +(3*(qCM(4)**2) - 
     9    7*qCM(4)*(qCM(3) - 5*qpCM(4)) +  
     1    qpCM(4)*(-11*qCM(3) + 16*qpCM(4)))*tt)))/
     2    (16.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(13)=(-3*sqrt(s)*((qsq**2) + 2*(2*(qCM(4)**2) - 
     1    qCM(4)*qCM(3) + qCM(3)*qpCM(4))* 
     2    (qpPerp**2) + (qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) - qCM(3) + 2*qpCM(4)) + 
     4    4*(qpPerp**2) + tt)))/
     5    (16.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(14)=-(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) + 3*qCM(4)*qCM(3) + qCM(4)*qpCM(4) 
     2    + qCM(3)*qpCM(4) - 2*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     3    8*qCM(4)*qCM(3)*(qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + qsq*(8*qCM(3)*qpCM(4)*
     5    ((qCM(4)**2) - (qpCM(4)**2)) + (4*(qCM(4)**2) 
     6    + 22*qCM(4)*qCM(3) - 6*qCM(3)*qpCM(4))*(qpPerp**2) + 
     7    ((qCM(4)**2) + 7*qCM(3)*qpCM(4) 
     8    + qCM(4)*(3*qCM(3) + qpCM(4)))*tt)))/
     9    (16.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(15)=((qCM(4) - qpCM(4))*qpPerp*sqrt(s)*
     1    ((-qCM(4) + qCM(3) + 2*qpCM(4))*qsq + (-qCM(4) + qCM(3))*tt))/
     2    (4.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(16)=-(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) + 3*qCM(4)*qpCM(4) - qCM(3)*qpCM(4) 
     2    - 3*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     3    2*qCM(4)*(qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(2*(qCM(4) -qCM(3))*qpCM(4)* ((qCM(4)**2) -(qpCM(4)**2)) + 
     6    2*qCM(4)*(3*qCM(4) - 2*qCM(3) + qpCM(4))*
     7    (qpPerp**2) + ((qCM(4)**2) +3*qCM(4)*qpCM(4) -qCM(3)*qpCM(4) + 
     8    (qpCM(4)**2))*tt)))/
     9    (2.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(17)=(qpPerp*((qCM(4) - qCM(3) - 2*qpCM(4))*qsq + 
     1    qCM(4)*tt))/(8.*sqrt(2.)*qCM(3)*qsq)

       Ur(18)=-((3*qCM(4) - 7*qCM(3) + 10*qpCM(4))*(qsq**2) + 
     1    8*qCM(4)*(4*qCM(4) - qCM(3))*(qCM(4)*(qpPerp**2)+ qpCM(4)*tt) + 
     2    qsq*(8*(4*qCM(4) - qCM(3))* (qCM(4) - qpCM(4))*qpCM(4) + 
     3    (32*qCM(4) - 6*qCM(3))*(qpPerp**2) + 
     4    (3*qCM(4) - 7*qCM(3) + 16*qpCM(4))*tt))/
     5    (32.*sqrt(2.)*qCM(3)*qsq)

       Ur(19)=(3*((-qCM(4) + qCM(3) + 2*qpCM(4))*qsq + 
     1    2*qCM(3)*(qpPerp**2) + (-qCM(4) + qCM(3))*tt))/
     2    (32.*sqrt(2.)*qCM(3)*qsq)

       Ur(20)=((qCM(4) + 3*qCM(3) - 2*qpCM(4))*(qsq**2) + 
     1    qsq*(-2*qCM(3)* (4*qpCM(4)*(-qCM(4) +qpCM(4)) + (qpPerp**2)) + 
     2    (qCM(4) + 3*qCM(3))*tt) + 
     3    8*qCM(4)*qCM(3)*(qCM(4)*(qpPerp**2) + qpCM(4)*tt)
     4    )/(32.*sqrt(2.)*qCM(3)*qsq)

       Ur(21)=(qpPerp*((qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1    (qCM(4) - 3*qCM(3))*tt))/(8.*sqrt(2.)*qCM(3))

       Ur(22)=(qpPerp*((-qCM(4) + qCM(3) + 2*qpCM(4))*qsq + 
     1    (-qCM(4) + qCM(3))*tt))/(8.*sqrt(2.)*qCM(3)*qsq)

       Ur(23)=(qpPerp*(3*(qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1    (3*qCM(4) - qCM(3))*tt))/(8.*sqrt(2.)*qCM(3))

       Ur(24)=((qCM(4) - qpCM(4))*(qsq**2) + 2*qCM(4)*(qCM(4) - qCM(3))*
     1    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     2    qsq*(2*(qCM(4) - qCM(3))*(qCM(4) - qpCM(4))
     3    *qpCM(4) + 2*(qCM(4) + qCM(3))*
     4    (qpPerp**2) + (qCM(4) + qpCM(4))*tt))/(4.*sqrt(2.)*qCM(3)*qsq)

       Ur(25)=(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*(3*(qsq**2) + 2*(qCM(4)**2)*tt - 
     2    qsq*(6*qCM(4)*(-qCM(4) + qpCM(4)) + tt)))/
     3    (32.*sqrt(2.)*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*
     4    sqrt((qsq**2) - 4*s*tt))

       Ur(26)=-(qpPerp*(2*qCM(4)*qpCM(4) + qsq + tt)*
     1    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))/
     2    (32.*sqrt(2.)*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) -4*s*tt))

       Ur(27)=((3*(qCM(4) - qpCM(4))*qsq + 4*qCM(4)*(qpPerp**2) + 
     1    3*(qCM(4) + qpCM(4))*tt)*
     2    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))
     3    /(16.*sqrt(2.)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(28)=-(((qCM(4) - qpCM(4))*qsq + 4*qCM(4)*(qpPerp**2) + 
     1    (qCM(4) + qpCM(4))*tt)*
     2    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)
     3    )/(16.*sqrt(2.)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(29)=-(qCM(4)*(qpPerp**2)*((qsq**2) + 
     1    2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))/
     2    (4.*sqrt(2.)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(30)=-(sqrt(s)*((qCM(4) + qpCM(4))*(qsq**2) + 
     1    4*qCM(4)*(qCM(4) +qpCM(4))* (qCM(4)*(qpPerp**2) +qpCM(4)*tt) + 
     2    qsq*(4*(qCM(4)**2)*qpCM(4) - 4*(qpCM(4)**3) + 
     3    2*(5*qCM(4) - qpCM(4))*(qpPerp**2) + qCM(4)*tt + 
     4    3*qpCM(4)*tt)))/(8.*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(31)=-(sqrt(s)*((qCM(4) - qpCM(4))*qsq + 
     1    2*(qCM(4) - qpCM(4))*(qpPerp**2) + (qCM(4) + qpCM(4))*tt))/
     2    (8.*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(32)=-(qpPerp*sqrt(s)*(qCM(4)*(2*qCM(4) - qpCM(4))*tt + 
     1    qsq*((qCM(4)**2) - 3*qCM(4)*qpCM(4) + 2*(qpCM(4)**2) + tt)))/
     2    (2.*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(33)=-(qpPerp*sqrt(s)*(qCM(4)*qpCM(4)*tt + 
     1    qsq*(-(qCM(4)**2) + 3*qCM(4)*qpCM(4) - 2*(qpCM(4)**2) + tt)))/
     2    (2.*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(34)=-(qpPerp*sqrt(s)*(qCM(4)*(qCM(4) + qpCM(4))*tt + 
     1    qsq*(-(qCM(4)**2) + 3*qCM(4)*qpCM(4) -2*(qpCM(4)**2) +2*tt)))/
     2    (2.*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(35)=((qCM(4) - 2*qpCM(4))*(Q**5) + qCM(4)*(Q**3)*tt)/
     1    (16.*qCM(3)*(Q**5))
       Ur(36)=(qpPerp*(qsq + 2*tt))/(4.*qsq)
       Ur(37)=-qpPerp/4.

       Ur(38)=(qCM(4)*(qsq**2) + 4*(qCM(4)**2)*
     1    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     2    qsq*(4*qCM(4)*(qCM(4) - qpCM(4))*qpCM(4) + 
     3    4*qCM(4)*(qpPerp**2) + (qCM(4) + 2*qpCM(4))*tt))/
     4    (8.*qCM(3)*qsq)
       Ur(39)=(qpPerp*(-qsq + tt))/(4.*qsq)

       Ur(40)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + 4*qCM(4)*(4*qCM(4) + 3*qCM(3))*
     2    (qpPerp**2) + (qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) - qCM(3) + 2*qpCM(4)) + 
     4    16*(qpPerp**2) + tt)))/
     5    (64.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(41)=-(qpPerp*((qsq**2) + qsq*
     1    (12*(qCM(4)**2) + 16*qCM(4)*qCM(3) - 
     2    22*qCM(4)*qpCM(4) - 26*qCM(3)*qpCM(4) - 7*tt)
     3    + 4*qCM(4)*(qCM(4) + 4*qCM(3))*tt)*
     4    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)
     5    )/(128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*
     6    sqrt((qsq**2) - 4*s*tt))

       Ur(42)=(-3*qsq*qpPerp*(2*(qCM(4) - qCM(3))*qpCM(4) + qsq + tt)*
     1    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))/
     2    (128.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(43)=-(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*(5*(qsq**2) + 4*(qCM(4)**2)*tt - 
     2    qsq*(2*(-6*(qCM(4)**2) + 7*qCM(4)*qpCM(4) + 
     3    qCM(3)*qpCM(4)) + 3*tt)))/
     4    (128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(44)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + qsq*((qCM(4) - qpCM(4))*
     2    (qCM(4) + qCM(3) + 2*qpCM(4)) + 2*(qpPerp**2) + tt) + 
     3    (qCM(4) + qCM(3))* 
     4    (2*qCM(4)*(qpPerp**2) + (qCM(4) + qpCM(4))*tt)))/
     5    (32.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(45)=(qpPerp*(qCM(4)*(3*qCM(4) + qCM(3) - 2*qpCM(4))*
     1    qsq + 2*(qsq**2) + qCM(4)*(qCM(4) + qCM(3))*tt)*
     2    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     3    4*s*tt))/(16.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*
     4    sqrt((qsq**2) - 4*s*tt))

       Ur(46)=-(qpPerp*sqrt(s)*((4*(qCM(4)**2) + 4*qCM(3)*qpCM(4) + 
     1    3*qCM(4)*(qCM(3) + qpCM(4)))*tt + 
     2    qsq*(-3*(qCM(4) - qCM(3) - 2*qpCM(4))*
     3    (qCM(4) - qpCM(4)) + 7*tt)))/
     4    (4.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(47)=(sqrt(s)*(3*(qsq**3) + (qsq**2)*
     1    (3*(qCM(4)**2) + 7*qCM(4)*qCM(3) + 
     2    35*qCM(4)*qpCM(4) - 3*qCM(3)*qpCM(4) - 
     3    22*(qpCM(4)**2) + 12*(qpPerp**2) + 3*tt) + 
     4    8*qCM(4)*(4*qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*
     5    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     6    qsq*(8*(4*qCM(4) + qCM(3))*qpCM(4)* 
     7    ((qCM(4)**2) - (qpCM(4)**2)) + 
     8    2*(22*(qCM(4)**2) - 7*qCM(3)*qpCM(4) + 
     9    qCM(4)*(15*qCM(3) + 16*qpCM(4)))*(qpPerp**2)
     1    + (3*(qCM(4)**2) + 7*qCM(4)*(qCM(3) + 5*qpCM(4)) + 
     2    qpCM(4)*(11*qCM(3) + 16*qpCM(4)))*tt)))/
     3    (16.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(48)=(-3*qsq*sqrt(s)*((qsq**2) + 
     1    2*(2*(qCM(4)**2) +qCM(4)*qCM(3) -qCM(3)*qpCM(4))*(qpPerp**2) + 
     2    (qCM(4) + qCM(3))*(qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) + qCM(3) + 2*qpCM(4)) + 
     4    4*(qpPerp**2) + tt)))/ 
     5    (16.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(49)=-(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) -3*qCM(4)*qCM(3) +qCM(4)*qpCM(4) -qCM(3)*qpCM(4) - 
     2    2*(qpCM(4)**2) + 4*(qpPerp**2) + tt) - 
     3    8*qCM(4)*qCM(3)*(qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + qsq*(8*qCM(3)*qpCM(4)*
     5    (-(qCM(4)**2) + (qpCM(4)**2)) + (4*(qCM(4)**2) 
     6    - 22*qCM(4)*qCM(3) + 6*qCM(3)*qpCM(4))*(qpPerp**2) + 
     7    ((qCM(4)**2) - 7*qCM(3)*qpCM(4) + 
     8    qCM(4)*(-3*qCM(3) + qpCM(4)))*tt)))/
     9    (16.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(50)=((qCM(4) - qpCM(4))*qpPerp*sqrt(s)*
     1    ((qCM(4) + qCM(3) - 2*qpCM(4))*qsq + (qCM(4) + qCM(3))*tt))/
     2    (4.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(51)=-(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) + 3*qCM(4)*qpCM(4) + qCM(3)*qpCM(4) 
     2    - 3*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     3    2*qCM(4)*(qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(2*(qCM(4) +qCM(3))*qpCM(4)* ((qCM(4)**2) -(qpCM(4)**2)) + 
     6    2*qCM(4)*(3*qCM(4) + 2*qCM(3) + qpCM(4))*
     7    (qpPerp**2) + ((qCM(4)**2) + 3*qCM(4)*qpCM(4) + 
     8    qpCM(4)*(qCM(3) + qpCM(4)))*tt)))/
     9    (2.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(52)=(qpPerp*((qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1    qCM(4)*tt))/(8.*sqrt(2.)*qCM(3))

       Ur(53)=((3*qCM(4) + 7*qCM(3) + 10*qpCM(4))*(qsq**2) + 
     1    8*qCM(4)*(4*qCM(4) + qCM(3))* 
     2    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) +qsq*(8*(4*qCM(4) + qCM(3))*
     3    (qCM(4) - qpCM(4))*qpCM(4) + 
     4    (32*qCM(4) + 6*qCM(3))*(qpPerp**2) + 
     5    (3*qCM(4) + 7*qCM(3) + 16*qpCM(4))*tt))/
     6    (32.*sqrt(2.)*qCM(3)*(qsq**2))

       Ur(54)=(3*qsq*((qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1    2*qCM(3)*(qpPerp**2) + (qCM(4) + qCM(3))*tt))/
     2    (32.*sqrt(2.)*qCM(3))

       Ur(55)=((-qCM(4) + 3*qCM(3) + 2*qpCM(4))*(qsq**2) -qsq*(2*qCM(3)*
     1    (4*qpCM(4)*(-qCM(4) + qpCM(4)) + (qpPerp**2)) + 
     2    (qCM(4) - 3*qCM(3))*tt) + 
     3    8*qCM(4)*qCM(3)*(qCM(4)*(qpPerp**2) + qpCM(4)*tt)
     4    )/(32.*sqrt(2.)*qCM(3)*(qsq**2))

       Ur(56)=(qpPerp*((qCM(4) - qCM(3) - 2*qpCM(4))*qsq + 
     1    (qCM(4) + 3*qCM(3))*tt))/(8.*sqrt(2.)*qCM(3)*qsq)

       Ur(57)=-(qpPerp*((qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1   (qCM(4) + qCM(3))*tt))/(8.*sqrt(2.)*qCM(3))

       Ur(58)=(qpPerp*(3*(qCM(4) - qCM(3) - 2*qpCM(4))*qsq + 
     1    (3*qCM(4) + qCM(3))*tt))/(8.*sqrt(2.)*qCM(3)*qsq)

       Ur(59)=-((qCM(4) -qpCM(4))*(qsq**2) + 2*qCM(4)*(qCM(4) + qCM(3))*
     1    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     2    qsq*(2*(qCM(4) + qCM(3))*(qCM(4) - qpCM(4))
     3    *qpCM(4) + 2*(qCM(4) - qCM(3))*(qpPerp**2) + 
     4    (qCM(4) + qpCM(4))*tt))/(4.*sqrt(2.)*qCM(3)*(qsq**2))

       Ur(60)=-(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*(11*(qsq**2) + 4*qCM(4)*(qCM(4) - 4*qCM(3))*tt + 
     2    qsq*(2*(6*(qCM(4)**2) + 3*qCM(3)*qpCM(4) - 
     3    qCM(4)*(8*qCM(3) + qpCM(4))) + 3*tt)))/
     4    (128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(61)=-(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1   4*s*tt)*(7*(qsq**2) + 4*(qCM(4)**2)*tt - 
     2   qsq*(2*(-6*(qCM(4)**2) + 5*qCM(4)*qpCM(4) + 
     3   qCM(3)*qpCM(4)) + tt)))/
     4   (128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(62)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + 2*qCM(4)*(3*qCM(4) - qCM(3))*
     2    (qpPerp**2) + (qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) + qCM(3) + 2*qpCM(4)) + 
     4    6*(qpPerp**2) + tt)))/
     5    (32.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(63)=(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) - 2*qCM(4)*(qCM(4) - 3*qCM(3))*
     2    (qpPerp**2) + (qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) + qCM(3) + 2*qpCM(4)) - 
     4    2*(qpPerp**2) + tt)))/
     5    (32.*qCM(3)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(64)=(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*((qsq**2) + qCM(4)*(qCM(4) - qCM(3))*tt - 
     2    qsq*(-3*(qCM(4)**2) - 2*qCM(3)*qpCM(4) + 
     3    qCM(4)*(qCM(3) + 4*qpCM(4)) + tt)))/
     4    (16.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(65)=-(sqrt(s)*(13*(qsq**3) + (qsq**2)*
     1    (13*(qCM(4)**2) + 3*qCM(4)*qCM(3) + 45*qCM(4)*qpCM(4) - 
     2    7*qCM(3)*qpCM(4) - 42*(qpCM(4)**2) + 52*(qpPerp**2) + 13*tt) + 
     3    8*qCM(4)*(4*qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(8*(4*qCM(4) - qCM(3))*qpCM(4)* 
     6    ((qCM(4)**2) - (qpCM(4)**2)) + 
     7    2*(42*(qCM(4)**2) - 3*qCM(3)*qpCM(4) + 
     8    qCM(4)*(-5*qCM(3) + 16*qpCM(4)))*(qpPerp**2) +
     9    (13*(qCM(4)**2) + 3*qCM(4)*(qCM(3) + 15*qpCM(4)) + 
     1    qpCM(4)*(-qCM(3) + 16*qpCM(4)))*tt)))/
     2    (16.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(66)=-(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) - qCM(4)*qCM(3) + qCM(4)*qpCM(4) - 
     2    3*qCM(3)*qpCM(4) - 2*(qpCM(4)**2) + 4*(qpPerp**2) + tt) - 
     3    8*qCM(4)*qCM(3)*(qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + qsq*(8*qCM(3)*qpCM(4)*
     5    (-(qCM(4)**2) + (qpCM(4)**2)) + 
     6    2*(2*(qCM(4)**2) - 9*qCM(4)*qCM(3) + qCM(3)*qpCM(4))*
     7    (qpPerp**2) + ((qCM(4)**2) - 5*qCM(3)*qpCM(4) + 
     8    qCM(4)*(-qCM(3) + qpCM(4)))*tt)))/
     9    (16.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(67)=(qpPerp*sqrt(s)*((3*(qCM(4)**2) - 3*qCM(3)*qpCM(4) + 
     1    qCM(4)*(-qCM(3) + qpCM(4)))*tt + 
     2    qsq*(-((qCM(4) + qCM(3) - 2*qpCM(4))*
     3    (qCM(4) - qpCM(4))) + 4*tt)))/
     4    (4.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(68)=(qpPerp*sqrt(s)*(((qCM(4)**2) -3*qCM(4)*(qCM(3)-qpCM(4)) - 
     1    qCM(3)*qpCM(4))*tt + qsq*(-3*(qCM(4) + qCM(3) - 2*qpCM(4))*
     2    (qCM(4) - qpCM(4)) + 4*tt)))/
     3    (4.*sqrt(2.)*qCM(3)*sqrt((qsq**2) - 4*s*tt))

       Ur(69)=(sqrt(s)*(-((qCM(4)*(qCM(3) - 2*qpCM(4)) + 
     1    (qpCM(4)**2))*(qsq**2)) + 2*qCM(4)*(qCM(4) - qCM(3))*
     2    (qCM(4) + qpCM(4))* (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     3    qsq*(2*(qCM(4) -qCM(3))*qpCM(4)* ((qCM(4)**2) -(qpCM(4)**2)) + 
     4    2*((qCM(4)**2) + qCM(3)*qpCM(4) + 
     5    qCM(4)*(-3*qCM(3) + qpCM(4)))*(qpPerp**2) + 
     6    (-(qCM(4)*qCM(3)) + 2*qCM(4)*qpCM(4) - 
     7    2*qCM(3)*qpCM(4) + (qpCM(4)**2))*tt)))/
     8    (2.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(70)=(qpPerp*(3*(qCM(4) + qCM(3) - 2*qpCM(4))*qsq + 
     1    (3*qCM(4) - 4*qCM(3))*tt))/(8.*sqrt(2.)*qCM(3))

       Ur(71)=((13*qCM(4) + 3*qCM(3) - 10*qpCM(4))*(qsq**2) + 
     1    8*qCM(4)*(4*qCM(4) - qCM(3))*
     2    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + qsq*(8*(4*qCM(4) -qCM(3))*
     3    (qCM(4) - qpCM(4))*qpCM(4) + 2*(16*qCM(4) + 7*qCM(3))*
     4    (qpPerp**2) + (13*qCM(4) + 3*qCM(3) + 16*qpCM(4))*tt))/
     5    (32.*sqrt(2.)*qCM(3)*qsq)

       Ur(72)=((qCM(4) - qCM(3) - 2*qpCM(4))*(qsq**2) + 
     1    qsq*(2*qCM(3)*(4*qpCM(4)*(-qCM(4) +qpCM(4)) + 3*(qpPerp**2)) + 
     2   (qCM(4) - qCM(3))*tt) - 
     3   8*qCM(4)*qCM(3)*(qCM(4)*(qpPerp**2) + qpCM(4)*tt)
     4   )/(32.*sqrt(2.)*qCM(3)*qsq)

       Ur(73)=((qCM(3) - qpCM(4))*(qsq**2) + 
     1    qsq*(-2*(qCM(4) - qCM(3))*(qCM(4) - qpCM(4))
     2    *qpCM(4) - 2*qCM(4)*(qpPerp**2) + 
     3   (qCM(3) - qpCM(4))*tt) - 2*qCM(4)*(qCM(4) - qCM(3))*
     4   (qCM(4)*(qpPerp**2) + qpCM(4)*tt))/
     5   (4.*sqrt(2.)*qCM(3)*qsq)

       Ur(74)=(qpCM(4)*qpPerp*((qsq**2) + 
     1    2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))/
     2    (16.*sqrt(2.)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(75)=(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    (3*(qsq**2) + qCM(4)*
     2    (8*qCM(4)*(qpPerp**2) + 3*(qCM(4) + qpCM(4))*tt) + 
     3    qsq*(8*(qpPerp**2) + 3*((qCM(4)**2) + qCM(4)*qpCM(4) - 
     4    2*(qpCM(4)**2) + tt))))/
     5    (16.*sqrt(2.)*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*
     6    sqrt((qsq**2) - 4*s*tt))

       Ur(76)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + qCM(4)*(qCM(4) + qpCM(4))*tt + 
     2    qsq*((qCM(4)**2) + qCM(4)*qpCM(4) - 2*(qpCM(4)**2) + tt)))/
     3    (16.*sqrt(2.)*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*
     4    sqrt((qsq**2) - 4*s*tt))

       Ur(77)=-(qpPerp*((qCM(4) - qpCM(4))*qsq + qCM(4)*tt)*
     1    ((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)
     2    )/(8.*sqrt(2.)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(78)=(qCM(3)*(qpPerp**2)*((qsq**2) + 
     1    2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt))/
     2    (4.*sqrt(2.)*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(79)=-(sqrt(s)*((qsq**2) + qsq*((qCM(4)**2) + qCM(4)*qpCM(4) - 
     1    2*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     2    qCM(4)*(4*qCM(4)*(qpPerp**2) + (qCM(4) + qpCM(4))*tt)))/
     3    (8.*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(80)=-(qpPerp*sqrt(s)*((qCM(4) - qpCM(4))*qsq + 
     1   (qCM(4) - 2*qpCM(4))*tt))/
     2   (2.*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(81)=(qpPerp*sqrt(s)*((qCM(4) - qpCM(4))*qsq + qCM(4)*tt))/
     1    (2.*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(82)=(sqrt(s)*((qsq**3) +(qsq**2)*
     1    ((qCM(4)**2) +5*qCM(4)*qpCM(4) - 
     2    4*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     3    4*(qCM(4)**2)*(qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(4*(qCM(4)**3)*qpCM(4) - 4*qCM(4)*(qpCM(4)**3) + 
     6    4*qCM(4)*(2*qCM(4) + qpCM(4))*(qpPerp**2) + 
     7    (qCM(4)**2)*tt + 5*qCM(4)*qpCM(4)*tt + 
     8    2*(qpCM(4)**2)*tt)))/
     9    (4.*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(83)=(qpPerp*sqrt(s)*((qCM(4) - qpCM(4))*qsq + 
     1    (qCM(4) + qpCM(4))*tt))/
     2    (2.*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(84)=-((qsq**2) + qsq*(4*(qCM(4) - qpCM(4))*qpCM(4) - 
     1    2*(qpPerp**2) + tt) + 
     2    4*qCM(4)*(qCM(4)*(qpPerp**2) + qpCM(4)*tt))/
     3    (16.*qsq)

       Ur(85)=(qsq + 2*(qpPerp**2) + tt)/(16.*qsq)

       Ur(86)=(qpPerp*((qCM(4) - 2*qpCM(4))*qsq + qCM(4)*tt))/
     1    (4.*qCM(3)*qsq)

       Ur(87)=-(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*(11*(qsq**2) + 4*qCM(4)*(qCM(4) + 4*qCM(3))*tt + 
     2    qsq*(12*(qCM(4)**2) + 16*qCM(4)*qCM(3) - 
     3    2*qCM(4)*qpCM(4) - 6*qCM(3)*qpCM(4) + 3*tt)))/
     4    (128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(88)=-(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*(7*(qsq**2) + 4*(qCM(4)**2)*tt - 
     2    qsq*(-2*(6*(qCM(4)**2) - 5*qCM(4)*qpCM(4) + 
     3    qCM(3)*qpCM(4)) + tt)))/
     4    (128.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(89)=(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) + 2*qCM(4)*(3*qCM(4) + qCM(3))*
     2    (qpPerp**2) + (qCM(4) - qCM(3))*(qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) - qCM(3) + 2*qpCM(4)) + 
     4    6*(qpPerp**2) + tt)))/
     5    (32.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(90)=-(((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 4*s*tt)*
     1    ((qsq**2) - 2*qCM(4)*(qCM(4) + 3*qCM(3))*
     2    (qpPerp**2) + (qCM(4) - qCM(3))* (qCM(4) + qpCM(4))*tt + 
     3    qsq*((qCM(4) - qpCM(4))* (qCM(4) - qCM(3) + 2*qpCM(4)) - 
     4    2*(qpPerp**2) + tt)))/
     5    (32.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(91)=(qpPerp*((qsq**2) + 2*(qCM(4) - qpCM(4))*qsq*sqrt(s) - 
     1    4*s*tt)*((qsq**2) + qCM(4)*(qCM(4) + qCM(3))*tt - 
     2    qsq*(-3*(qCM(4)**2) - qCM(4)*qCM(3) + 
     3    4*qCM(4)*qpCM(4) + 2*qCM(3)*qpCM(4) + tt)))/
     4    (16.*qCM(3)*qsq*sqrt(s)*sqrt(-tt)*sqrt((qsq**2) - 4*s*tt))

       Ur(92)=(sqrt(s)*(13*(qsq**3) + (qsq**2)*
     1    (13*(qCM(4)**2) - 3*qCM(4)*qCM(3) + 45*qCM(4)*qpCM(4) + 
     2    7*qCM(3)*qpCM(4) - 42*(qpCM(4)**2) + 52*(qpPerp**2) + 13*tt) + 
     3    8*qCM(4)*(4*qCM(4) + qCM(3))* (qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     5    qsq*(8*(4*qCM(4) + qCM(3))*qpCM(4)* 
     6    ((qCM(4)**2) - (qpCM(4)**2)) + 
     7    2*(42*(qCM(4)**2) + 3*qCM(3)*qpCM(4) + 
     8    qCM(4)*(5*qCM(3) + 16*qpCM(4)))*(qpPerp**2)
     9    + (13*(qCM(4)**2) - 3*qCM(4)*(qCM(3) - 15*qpCM(4)) + 
     1    qpCM(4)*(qCM(3) + 16*qpCM(4)))*tt)))/
     2    (16.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(93)=(sqrt(s)*((qsq**3) + (qsq**2)*
     1    ((qCM(4)**2) + qCM(4)*qCM(3) + qCM(4)*qpCM(4) + 
     2    3*qCM(3)*qpCM(4) - 2*(qpCM(4)**2) + 4*(qpPerp**2) + tt) + 
     3    8*qCM(4)*qCM(3)*(qCM(4) + qpCM(4))*
     4    (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + qsq*(8*qCM(3)*qpCM(4)*
     5    ((qCM(4)**2) - (qpCM(4)**2)) + 2*(2*(qCM(4)**2) 
     6    + 9*qCM(4)*qCM(3) - qCM(3)*qpCM(4))*(qpPerp**2) + 
     7    ((qCM(4)**2) + 5*qCM(3)*qpCM(4) + 
     8    qCM(4)*(qCM(3) + qpCM(4)))*tt)))/
     9    (16.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(94)=(qpPerp*sqrt(s)*((3*(qCM(4)**2) + 3*qCM(3)*qpCM(4) + 
     1    qCM(4)*(qCM(3) + qpCM(4)))*tt + 
     2    qsq*(-((qCM(4) - qCM(3) - 2*qpCM(4))*
     3    (qCM(4) - qpCM(4))) + 4*tt)))/
     4    (4.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(95)=(qpPerp*sqrt(s)*(((qCM(4)**2) + qCM(3)*qpCM(4) + 
     1    3*qCM(4)*(qCM(3) + qpCM(4)))*tt + 
     2    qsq*(-3*(qCM(4) - qCM(3) - 2*qpCM(4))*
     3    (qCM(4) - qpCM(4)) + 4*tt)))/
     4    (4.*sqrt(2.)*qCM(3)*qsq*sqrt((qsq**2) - 4*s*tt))

       Ur(96)=-(sqrt(s)*((-(qpCM(4)**2) + qCM(4)*
     1    (qCM(3) + 2*qpCM(4)))*(qsq**2) + 2*qCM(4)*(qCM(4) + qCM(3))*
     2    (qCM(4) + qpCM(4))* (qCM(4)*(qpPerp**2) + qpCM(4)*tt) + 
     3    qsq*(2*(qCM(4) + qCM(3))*qpCM(4)*
     4    ((qCM(4)**2) - (qpCM(4)**2)) + 
     5    2*((qCM(4)**2) - qCM(3)*qpCM(4) + 
     6    qCM(4)*(3*qCM(3) + qpCM(4)))*(qpPerp**2) + 
     7    (qpCM(4)*(2*qCM(3) + qpCM(4)) + 
     8    qCM(4)*(qCM(3) + 2*qpCM(4)))*tt)))/
     9    (2.*sqrt(2.)*qCM(3)*(qsq**2)*sqrt((qsq**2) - 4*s*tt))

       Ur(97)=(qpPerp*((-3*qCM(4) + 3*qCM(3) + 6*qpCM(4))*qsq - 
     1    (3*qCM(4) + 4*qCM(3))*tt))/(8.*sqrt(2.)*qCM(3)*qsq)

       Ur(98)=((13*qCM(4) - 3*qCM(3) - 10*qpCM(4))*(qsq**2) + 
     1    8*qCM(4)*(4*qCM(4) +qCM(3))*(qCM(4)*(qpPerp**2) +qpCM(4)*tt) + 
     2    qsq*(8*(4*qCM(4) + qCM(3))* (qCM(4) - qpCM(4))*qpCM(4) + 
     3    2*(16*qCM(4) - 7*qCM(3))*(qpPerp**2) + 
     4   (13*qCM(4) - 3*qCM(3) + 16*qpCM(4))*tt))/
     5   (32.*sqrt(2.)*qCM(3)*(qsq**2))

       Ur(99)=((qCM(4) + qCM(3) - 2*qpCM(4))*(qsq**2) + 
     1    qsq*(-2*qCM(3)*
     2    (4*qpCM(4)*(-qCM(4) + qpCM(4)) + 3*(qpPerp**2)) + 
     3    (qCM(4) + qCM(3))*tt) + 
     4    8*qCM(4)*qCM(3)*(qCM(4)*(qpPerp**2) + qpCM(4)*tt))
     5    /(32.*sqrt(2.)*qCM(3)*(qsq**2))

       Ur(100)=-((qCM(3) + qpCM(4))*(qsq**2) + 
     1    2*qCM(4)*(qCM(4) + qCM(3))* (qCM(4)*(qpPerp**2) +qpCM(4)*tt) + 
     2    qsq*(2*(qCM(4) + qCM(3))*(qCM(4) - qpCM(4))
     3    *qpCM(4) + 2*qCM(4)*(qpPerp**2) + (qCM(3) + qpCM(4))*tt))/
     4    (4.*sqrt(2.)*qCM(3)*(qsq**2))

      endif

c iexact=0:
c All quantities are evaluated at leading order in the 1/Q expansion 
c (see Belitsky, Mueller and Kirchner, Nucl.Phys. B629 (2002) 323-392,
c ArXiv:hep-ph/0112108v2)

      if (iexact.eq.0) then

       Ur(1)=(-13*(qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1    (1 - xB)))/16.
       Ur(2)=sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(1 - xB))/16.
       Ur(3)=(-3*(qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1    (1 - xB)))/16.
       Ur(4)=(7*sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1    (1 - xB)))/16.
       Ur(5)=(3*sqrt(-((-1 + xB)*((qpPerp**4) + 
     1    (m_p**2)*(qpPerp**2)*(xB**2)))))/(16.*xB)
       Ur(6)=(3*sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1    (1 - xB)))/16.
       Ur(7)=-sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(1 - xB))/8.
       Ur(8)=-sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1    (1 - xB))/2.
       Ur(9)=((qpPerp**3)*(13 - 6*xB) + 7*(m_p**2)*qpPerp*(xB**2))/
     1    (8.*sqrt(2.)*(-1 + xB))
       Ur(10)=qpPerp/(4.*sqrt(2.))
       Ur(11)=((qpPerp**3)*(3 - 2*xB) + (m_p**2)*qpPerp*(xB**2))/
     1    (8.*sqrt(2.)*(-1 + xB))
       Ur(12)=(-5*(m_p**2)*(xB**2) + (qpPerp**2)*(-11 + 6*xB))/
     1    (16.*sqrt(2.)*(-1 + xB))
       Ur(13)=-3/(16.*sqrt(2.))
       Ur(14)=(-((m_p**2)*(xB**2)) + (qpPerp**2)*(-7 + 6*xB))/
     1    (16.*sqrt(2.)*(-1 + xB))
       Ur(15)=-qpPerp/(4.*sqrt(2.))
       Ur(16)=((qpPerp**2)*(4 - 3*xB) + (m_p**2)*(xB**2))/
     1    (4.*sqrt(2.)*(-1 + xB))
       Ur(17)=-qpPerp/(4.*sqrt(2.))
       Ur(18)=((qpPerp**2)*(1 - 6*xB) - 5*(m_p**2)*(xB**2))/
     1    (16.*sqrt(2.)*(-1 + xB))
       Ur(19)=3/(16.*sqrt(2.))
       Ur(20)=((qpPerp**2)*(5 - 6*xB) - (m_p**2)*(xB**2))/
     1    (16.*sqrt(2.)*(-1 + xB))
       Ur(21)=-(qpPerp*(2*(m_p**2)*(xB**2) + (qpPerp**2)*(1 + xB)))/
     1    (4.*sqrt(2.)*(-1 + xB))
       Ur(22)=qpPerp/(4.*sqrt(2.))
       Ur(23)=((qpPerp**3)*(1 - 3*xB) - 2*(m_p**2)*qpPerp*(xB**2))/
     1    (4.*sqrt(2.)*(-1 + xB))
       Ur(24)=((m_p**2)*(xB**2) + (qpPerp**2)*(-2 + 3*xB))/
     1    (4.*sqrt(2.)*(-1 + xB))
       Ur(25)=(3*sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1    (2 - 2*xB)))/8.
       Ur(26)=-sqrt(-((-1 + xB)*((qpPerp**4) + 
     1    (m_p**2)*(qpPerp**2)*(xB**2))))/(8.*sqrt(2.)*xB)
       Ur(27)=(-3*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(2 - 2*xB)))/8.
       Ur(28)=sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(2 - 2*xB))/8.
       Ur(29)=((qpPerp**2)*(-1 + 2*xB)*sqrt(((qpPerp**2) + 
     1    (m_p**2)*(xB**2))/(2 - 2*xB)))/(2.*xB)
       Ur(30)=-((qpPerp**2)*(7 -6*xB) + (m_p**2)*(xB**2))/
     1    (16.*(-1 + xB))
       Ur(31)=0.0625
       Ur(32)=-qpPerp/4.
       Ur(33)=qpPerp/4.
       Ur(34)=qpPerp/4.
       Ur(35)=-0.0625
       Ur(36)=qpPerp/4.
       Ur(37)=-qpPerp/4.
       Ur(38)=((qpPerp**2) + (m_p**2)*(xB**2))/(-8 + 8*xB)
       Ur(39)=-qpPerp/4.
       Ur(40)=-sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(1 - xB))/16.
       Ur(41)=(qpPerp*(5 - 9*xB + 4*(xB**2))*
     1   sqrt(-(((qpPerp**2) + (m_p**2)*(xB**2))/((-1+xB)**3))))/
     1   (16.*xB)
       Ur(42)=(3*(qpPerp**3)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1   (1 - xB)))/16.
       Ur(43)=((1 - 4*xB)*sqrt(((qpPerp**4) + 
     1   (m_p**2)*(qpPerp**2)*(xB**2))/(1 - xB)))/(16.*xB)
       Ur(44)=((qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1   (1 - xB)))/8.
       Ur(45)=sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1   (1 - xB))/(4.*xB)
       Ur(46)=(3*qpPerp)/(4.*sqrt(2.))
       Ur(47)=-5/(16.*sqrt(2.))
       Ur(48)=(3*(qpPerp**4)*(-2 + xB) - 
     1   3*(m_p**2)*(qpPerp**2)*(xB**2))/(16.*sqrt(2.)*(-1 + xB))
       Ur(49)=-1/(16.*sqrt(2.))
       Ur(50)=(qpPerp**3)/(4.*sqrt(2.))
       Ur(51)=-1/(4.*sqrt(2.))
       Ur(52)=((qpPerp**3)*(1 - 2*xB) - (m_p**2)*qpPerp*(xB**2))/
     1   (8.*sqrt(2.)*(-1 + xB))
       Ur(53)=5/(16.*sqrt(2.))
       Ur(54)=(3*xB*((qpPerp**4) + (m_p**2)*(qpPerp**2)*xB))/
     1   (16.*sqrt(2.)*(-1 + xB))
       Ur(55)=1/(16.*sqrt(2.))
       Ur(56)=-qpPerp/(4.*sqrt(2.))
       Ur(57)=(qpPerp**3)/(4.*sqrt(2.))
       Ur(58)=(-3*qpPerp)/(4.*sqrt(2.))
       Ur(59)=1/(4.*sqrt(2.))
       Ur(60)=(-7*sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1   (1 - xB)))/16.
       Ur(61)=(-3*sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1   (1 - xB)))/16.
       Ur(62)=(-3*(qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1   (1 - xB)))/8.
       Ur(63)=(-5*(qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1   (1 - xB)))/8.
       Ur(64)=sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1   (1 - xB))/2.
       Ur(65)=((qpPerp**2)*(11 - 6*xB) + 5*(m_p**2)*(xB**2))/
     1   (16.*sqrt(2.)*(-1 + xB))
       Ur(66)=((qpPerp**2)*(7 - 6*xB) + (m_p**2)*(xB**2))/
     1   (16.*sqrt(2.)*(-1 + xB))
       Ur(67)=(-((qpPerp**3)*(-3 + xB)) + 2*(m_p**2)*qpPerp*(xB**2))/
     1   (4.*sqrt(2.)*(-1 + xB))
       Ur(68)=((qpPerp**3)*(5 - 3*xB) + 2*(m_p**2)*qpPerp*(xB**2))/
     1   (4.*sqrt(2.)*(-1 + xB))
       Ur(69)=(-((m_p**2)*(xB**2)) + (qpPerp**2)*(-4 + 3*xB))/
     1   (4.*sqrt(2.)*(-1 + xB))
       Ur(70)=-(qpPerp*(7*(m_p**2)*(xB**2) + (qpPerp**2)*(1 + 6*xB)))/
     1   (8.*sqrt(2.)*(-1 + xB))
       Ur(71)=(5*(m_p**2)*(xB**2) + (qpPerp**2)*(-1 + 6*xB))/
     1   (16.*sqrt(2.)*(-1 + xB))
       Ur(72)=((m_p**2)*(xB**2) + (qpPerp**2)*(-5 + 6*xB))/
     1   (16.*sqrt(2.)*(-1 + xB))
       Ur(73)=((qpPerp**2)*(2 - 3*xB) - (m_p**2)*(xB**2))/
     1   (4.*sqrt(2.)*(-1 + xB))
       Ur(74)=sqrt(-((-1 + xB)*((qpPerp**4) + 
     1   (m_p**2)*(qpPerp**2)*(xB**2))))/(8.*sqrt(2.)*xB)
       Ur(75)=(3*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(2 - 2*xB)))/8.
       Ur(76)=-sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(2 - 2*xB))/8.
       Ur(77)=sqrt(((qpPerp**4) + (m_p**2)*(qpPerp**2)*(xB**2))/
     1   (2 - 2*xB))/4.
       Ur(78)=((qpPerp**2)*sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/
     1   (2 - 2*xB)))/(2.*xB)
       Ur(79)=-0.0625
       Ur(80)=qpPerp/4.
       Ur(81)=-qpPerp/4.
       Ur(82)=((qpPerp**2) + (m_p**2)*(xB**2))/(8 - 8*xB)
       Ur(83)=-qpPerp/4.
       Ur(84)=((m_p**2)*(xB**2) + (qpPerp**2)*(-5 + 6*xB))/
     1   (16.*(-1 + xB))
       Ur(85)=0.0625
       Ur(86)=-qpPerp/4.
       Ur(87)=-(qpPerp*(5 - 11*xB + 6*(xB**2))*
     1   sqrt(-(((qpPerp**2) + (m_p**2)*(xB**2))/((-1 + xB)**3))))/
     2   (16.*xB)
       Ur(88)=-((1 + 2*xB)*sqrt(((qpPerp**4) + 
     1   (m_p**2)*(qpPerp**2)*(xB**2))/(1 - xB)))/(16.*xB)
       Ur(89)=sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(1 - xB))/8.
       Ur(90)=-sqrt(((qpPerp**2) + (m_p**2)*(xB**2))/(1 - xB))/8.
       Ur(91)=-(qpPerp*(1 - 3*xB + 2*(xB**2))*
     1    sqrt(-(((qpPerp**2) +(m_p**2)*(xB**2))/((-1 +xB)**3))))/(4.*xB)
       Ur(92)=5/(16.*sqrt(2.))
       Ur(93)=1/(16.*sqrt(2.))
       Ur(94)=-qpPerp/(4.*sqrt(2.))
       Ur(95)=(-3*qpPerp)/(4.*sqrt(2.))
       Ur(96)=1/(4.*sqrt(2.))
       Ur(97)=(3*qpPerp)/(4.*sqrt(2.))
       Ur(98)=-5/(16.*sqrt(2.))
       Ur(99)=-1/(16.*sqrt(2.))
       Ur(100)=-1/(4.*sqrt(2.))

      endif

c---------------------------- VCS cross sections -------------------------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

      SigmaVCSPol0(1)=(3*((IMvcs(1,1))**2) - 2*((IMvcs(1,2))**2) + 
     1   3*((IMvcs(1,3))**2) + 3*((IMvcs(2,1))**2) - 
     2   2*((IMvcs(2,2))**2) + 3*((IMvcs(2,3))**2) + 
     3   3*((IMvcs(3,1))**2) - 2*((IMvcs(3,2))**2) + 
     4   3*((IMvcs(3,3))**2) + 3*((IMvcs(4,1))**2) - 
     5   2*((IMvcs(4,2))**2) + 3*((IMvcs(4,3))**2) + 
     6   3*((RMvcs(1,1))**2) - 2*((RMvcs(1,2))**2) + 
     7   3*((RMvcs(1,3))**2) + 3*((RMvcs(2,1))**2) - 
     8   2*((RMvcs(2,2))**2) + 3*((RMvcs(2,3))**2) + 
     9   3*((RMvcs(3,1))**2) - 2*((RMvcs(3,2))**2) + 
     1   3*((RMvcs(3,3))**2) + 3*((RMvcs(4,1))**2) - 
     2   2*((RMvcs(4,2))**2) + 3*((RMvcs(4,3))**2))/4.

      SigmaVCSPol0(2)=(((IMvcs(1,1))**2) + 2*((IMvcs(1,2))**2) + 
     1   ((IMvcs(1,3))**2) + ((IMvcs(2,1))**2) + 
     2   2*((IMvcs(2,2))**2) + ((IMvcs(2,3))**2) + 
     3   ((IMvcs(3,1))**2) + 2*((IMvcs(3,2))**2) + 
     4   ((IMvcs(3,3))**2) + ((IMvcs(4,1))**2) + 
     5   2*((IMvcs(4,2))**2) + ((IMvcs(4,3))**2) + 
     6   ((RMvcs(1,1))**2) + 2*((RMvcs(1,2))**2) + 
     7   ((RMvcs(1,3))**2) + ((RMvcs(2,1))**2) + 
     8   2*((RMvcs(2,2))**2) + ((RMvcs(2,3))**2) + 
     9   ((RMvcs(3,1))**2) + 2*((RMvcs(3,2))**2) + 
     1   ((RMvcs(3,3))**2) + ((RMvcs(4,1))**2) + 
     2   2*((RMvcs(4,2))**2) + ((RMvcs(4,3))**2))/4.

      SigmaVCSPol0(3)=(IMvcs(1,2)*(IMvcs(1,1) - IMvcs(1,3)) + 
     1   IMvcs(2,2)*(IMvcs(2,1) - IMvcs(2,3)) + 
     2   IMvcs(3,2)*(IMvcs(3,1) - IMvcs(3,3)) + 
     3   IMvcs(4,2)*(IMvcs(4,1) - IMvcs(4,3)) + 
     4   RMvcs(1,2)*(RMvcs(1,1) - RMvcs(1,3)) + 
     5   RMvcs(2,2)*(RMvcs(2,1) - RMvcs(2,3)) + 
     6   RMvcs(3,2)*(RMvcs(3,1) - RMvcs(3,3)) + 
     7   RMvcs(4,2)*(RMvcs(4,1) - RMvcs(4,3)))/sqrt(2.)

      SigmaVCSPol0(4)=(IMvcs(1,1)*IMvcs(1,3) + IMvcs(2,1)*IMvcs(2,3) + 
     1   IMvcs(3,1)*IMvcs(3,3) + IMvcs(4,1)*IMvcs(4,3) + 
     2   RMvcs(1,1)*RMvcs(1,3) + RMvcs(2,1)*RMvcs(2,3) + 
     3   RMvcs(3,1)*RMvcs(3,3) + RMvcs(4,1)*RMvcs(4,3))/2.

      SigmaVCSPol0(5)=sqrt(2.)*((IMvcs(1,1) - IMvcs(1,3))*RMvcs(1,2) + 
     1   IMvcs(1,2)*(-RMvcs(1,1) + RMvcs(1,3)) + 
     2   (IMvcs(2,1) - IMvcs(2,3))*RMvcs(2,2) + 
     3   IMvcs(2,2)*(-RMvcs(2,1) + RMvcs(2,3)) + 
     4   (IMvcs(3,1) - IMvcs(3,3))*RMvcs(3,2) + 
     5   IMvcs(3,2)*(-RMvcs(3,1) + RMvcs(3,3)) + 
     6   (IMvcs(4,1) - IMvcs(4,3))*RMvcs(4,2) + 
     7   IMvcs(4,2)*(-RMvcs(4,1) + RMvcs(4,3)))

      SigmaVCSPolX(1)=2*(IMvcs(1,1)*IMvcs(2,1) - IMvcs(1,3)*IMvcs(2,3) - 
     1   IMvcs(3,1)*IMvcs(4,1) + IMvcs(3,3)*IMvcs(4,3) + 
     2   RMvcs(1,1)*RMvcs(2,1) - RMvcs(1,3)*RMvcs(2,3) - 
     3   RMvcs(3,1)*RMvcs(4,1) + RMvcs(3,3)*RMvcs(4,3))

      SigmaVCSPolX(2)=sqrt(2.)*((IMvcs(1,1) + IMvcs(1,3))*IMvcs(2,2) + 
     1   IMvcs(1,2)*(IMvcs(2,1) + IMvcs(2,3)) - 
     2   (IMvcs(3,1) + IMvcs(3,3))*IMvcs(4,2) - 
     3   IMvcs(3,2)*(IMvcs(4,1) + IMvcs(4,3)) + 
     4   (RMvcs(1,1) + RMvcs(1,3))*RMvcs(2,2) + 
     5   RMvcs(1,2)*(RMvcs(2,1) + RMvcs(2,3)) - 
     6   (RMvcs(3,1) + RMvcs(3,3))*RMvcs(4,2) - 
     7   RMvcs(3,2)*(RMvcs(4,1) + RMvcs(4,3)))

      SigmaVCSPolX(3)=((IMvcs(2,1) + IMvcs(2,3))*RMvcs(1,2) - 
     1   IMvcs(2,2)*(RMvcs(1,1) + RMvcs(1,3)) + 
     2   (IMvcs(1,1) + IMvcs(1,3))*RMvcs(2,2) - 
     3   IMvcs(1,2)*(RMvcs(2,1) + RMvcs(2,3)) - 
     4   (IMvcs(4,1) + IMvcs(4,3))*RMvcs(3,2) + 
     5   IMvcs(4,2)*(RMvcs(3,1) + RMvcs(3,3)) - 
     6   (IMvcs(3,1) + IMvcs(3,3))*RMvcs(4,2) + 
     7   IMvcs(3,2)*(RMvcs(4,1) + RMvcs(4,3)))/sqrt(2.)

      SigmaVCSPolX(4)=(-(IMvcs(2,3)*RMvcs(1,1)) + IMvcs(2,1)*RMvcs(1,3) - 
     1   IMvcs(1,3)*RMvcs(2,1) + IMvcs(1,1)*RMvcs(2,3) + 
     2   IMvcs(4,3)*RMvcs(3,1) - IMvcs(4,1)*RMvcs(3,3) + 
     3   IMvcs(3,3)*RMvcs(4,1) - IMvcs(3,1)*RMvcs(4,3))/2.

      SigmaVCSPolY(2)=sqrt(2.)*((IMvcs(1,1) - IMvcs(1,3))*IMvcs(3,2) + 
     1   IMvcs(1,2)*(-IMvcs(3,1) + IMvcs(3,3)) + 
     2   (IMvcs(2,1) - IMvcs(2,3))*IMvcs(4,2) + 
     3   IMvcs(2,2)*(-IMvcs(4,1) + IMvcs(4,3)) + 
     4   (RMvcs(1,1) - RMvcs(1,3))*RMvcs(3,2) + 
     5   RMvcs(1,2)*(-RMvcs(3,1) + RMvcs(3,3)) + 
     6   (RMvcs(2,1) - RMvcs(2,3))*RMvcs(4,2) + 
     7   RMvcs(2,2)*(-RMvcs(4,1) + RMvcs(4,3)))

      SigmaVCSPolY(3)=(3*IMvcs(3,1)*RMvcs(1,1) - 2*IMvcs(3,2)*RMvcs(1,2) + 
     1   3*IMvcs(3,3)*RMvcs(1,3) + 3*IMvcs(4,1)*RMvcs(2,1) - 
     2   2*IMvcs(4,2)*RMvcs(2,2) + 3*IMvcs(4,3)*RMvcs(2,3) - 
     3   3*IMvcs(1,1)*RMvcs(3,1) + 2*IMvcs(1,2)*RMvcs(3,2) - 
     4   3*IMvcs(1,3)*RMvcs(3,3) - 3*IMvcs(2,1)*RMvcs(4,1) + 
     5   2*IMvcs(2,2)*RMvcs(4,2) - 3*IMvcs(2,3)*RMvcs(4,3))/2.

      SigmaVCSPolY(4)=(IMvcs(3,1)*RMvcs(1,1) + 2*IMvcs(3,2)*RMvcs(1,2) + 
     1   IMvcs(3,3)*RMvcs(1,3) + IMvcs(4,1)*RMvcs(2,1) + 
     2   2*IMvcs(4,2)*RMvcs(2,2) + IMvcs(4,3)*RMvcs(2,3) - 
     3   IMvcs(1,1)*RMvcs(3,1) - 2*IMvcs(1,2)*RMvcs(3,2) - 
     4   IMvcs(1,3)*RMvcs(3,3) - IMvcs(2,1)*RMvcs(4,1) - 
     5   2*IMvcs(2,2)*RMvcs(4,2) - IMvcs(2,3)*RMvcs(4,3))/2.

      SigmaVCSPolY(4)=((IMvcs(3,1) - IMvcs(3,3))*RMvcs(1,2) + 
     1   IMvcs(3,2)*(RMvcs(1,1) - RMvcs(1,3)) + 
     2   (IMvcs(4,1) - IMvcs(4,3))*RMvcs(2,2) + 
     3   IMvcs(4,2)*(RMvcs(2,1) - RMvcs(2,3)) + 
     4   (-IMvcs(1,1) + IMvcs(1,3))*RMvcs(3,2) + 
     5   IMvcs(1,2)*(-RMvcs(3,1) + RMvcs(3,3)) + 
     6   (-IMvcs(2,1) + IMvcs(2,3))*RMvcs(4,2) + 
     7   IMvcs(2,2)*(-RMvcs(4,1) + RMvcs(4,3)))/sqrt(2.)

      SigmaVCSPolY(5)=(IMvcs(3,3)*RMvcs(1,1) + IMvcs(3,1)*RMvcs(1,3) + 
     1   IMvcs(4,3)*RMvcs(2,1) + IMvcs(4,1)*RMvcs(2,3) - 
     2   IMvcs(1,3)*RMvcs(3,1) - IMvcs(1,1)*RMvcs(3,3) - 
     3   IMvcs(2,3)*RMvcs(4,1) - IMvcs(2,1)*RMvcs(4,3))/2.

      SigmaVCSPolZ(1)=2*(IMvcs(2,1)*IMvcs(3,1) - IMvcs(2,3)*IMvcs(3,3) + 
     1   IMvcs(1,1)*IMvcs(4,1) - IMvcs(1,3)*IMvcs(4,3) + 
     2   RMvcs(2,1)*RMvcs(3,1) - RMvcs(2,3)*RMvcs(3,3) + 
     3   RMvcs(1,1)*RMvcs(4,1) - RMvcs(1,3)*RMvcs(4,3))

      SigmaVCSPolZ(2)=sqrt(2.)*((IMvcs(2,1) + IMvcs(2,3))*IMvcs(3,2) + 
     1   IMvcs(2,2)*(IMvcs(3,1) + IMvcs(3,3)) + 
     2   (IMvcs(1,1) + IMvcs(1,3))*IMvcs(4,2) + 
     3   IMvcs(1,2)*(IMvcs(4,1) + IMvcs(4,3)) + 
     4   (RMvcs(2,1) + RMvcs(2,3))*RMvcs(3,2) + 
     5   RMvcs(2,2)*(RMvcs(3,1) + RMvcs(3,3)) + 
     6   (RMvcs(1,1) + RMvcs(1,3))*RMvcs(4,2) + 
     7   RMvcs(1,2)*(RMvcs(4,1) + RMvcs(4,3)))

      SigmaVCSPolZ(3)=((IMvcs(4,1) + IMvcs(4,3))*RMvcs(1,2) - 
     1   IMvcs(4,2)*(RMvcs(1,1) + RMvcs(1,3)) + 
     2   (IMvcs(3,1) + IMvcs(3,3))*RMvcs(2,2) - 
     3   IMvcs(3,2)*(RMvcs(2,1) + RMvcs(2,3)) + 
     4   (IMvcs(2,1) + IMvcs(2,3))*RMvcs(3,2) - 
     5   IMvcs(2,2)*(RMvcs(3,1) + RMvcs(3,3)) + 
     6   (IMvcs(1,1) + IMvcs(1,3))*RMvcs(4,2) - 
     7   IMvcs(1,2)*(RMvcs(4,1) + RMvcs(4,3)))/sqrt(2.)

      SigmaVCSPolZ(4)=(-(IMvcs(4,3)*RMvcs(1,1)) + IMvcs(4,1)*RMvcs(1,3) - 
     1   IMvcs(3,3)*RMvcs(2,1) + IMvcs(3,1)*RMvcs(2,3) - 
     2   IMvcs(2,3)*RMvcs(3,1) + IMvcs(2,1)*RMvcs(3,3) - 
     3   IMvcs(1,3)*RMvcs(4,1) + IMvcs(1,1)*RMvcs(4,3))/2.
          
c----------------------- Interference cross sections ---------------------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

      SigmaIPol0(1)=Jem(1,2)*(-2*RMvcs(1,1)*Ur(1)*qsq + 
     1   2*RMvcs(1,2)*Ur(25)*(Q**3) - 2*RMvcs(1,3)*Ur(40)*(Q**4)) + 
     2   Jem(3,2)*(-2*RMvcs(3,1)*Ur(1)*qsq + 
     3   2*RMvcs(3,2)*Ur(25)*(Q**3) - 2*RMvcs(3,3)*Ur(40)*(Q**4)) + 
     4   Jem(1,3)*(-2*RMvcs(1,1)*Ur(9)*qsq + 
     5   2*RMvcs(1,2)*Ur(30)*(Q**3) - 2*RMvcs(1,3)*Ur(46)*(Q**4)) + 
     6   Jem(3,3)*(-2*RMvcs(3,1)*Ur(9)*qsq + 
     7   2*RMvcs(3,2)*Ur(30)*(Q**3) - 2*RMvcs(3,3)*Ur(46)*(Q**4)) + 
     8   Jem(2,3)*(-2*RMvcs(2,1)*Ur(70)*qsq + 
     9   2*RMvcs(2,2)*Ur(84)*(Q**3) - 2*RMvcs(2,3)*Ur(97)*(Q**4)) + 
     1   Jem(4,3)*(-2*RMvcs(4,1)*Ur(70)*qsq + 
     2   2*RMvcs(4,2)*Ur(84)*(Q**3) - 2*RMvcs(4,3)*Ur(97)*(Q**4))

      SigmaIPol0(2)=Jem(1,2)*(-2*RMvcs(1,1)*Ur(3)*qsq - 
     1   2*RMvcs(1,2)*Ur(25)*(Q**3) - 2*RMvcs(1,3)*Ur(2)*(Q**4)) + 
     2   Jem(3,2)*(-2*RMvcs(3,1)*Ur(3)*qsq - 
     3   2*RMvcs(3,2)*Ur(25)*(Q**3) - 2*RMvcs(3,3)*Ur(2)*(Q**4)) + 
     4   Jem(2,3)*(-2*RMvcs(2,1)*Ur(52)*qsq - 
     5   2*RMvcs(2,2)*Ur(84)*(Q**3) + 2*RMvcs(2,3)*Ur(17)*(Q**4)) + 
     6   Jem(4,3)*(-2*RMvcs(4,1)*Ur(52)*qsq - 
     7   2*RMvcs(4,2)*Ur(84)*(Q**3) + 2*RMvcs(4,3)*Ur(17)*(Q**4)) + 
     8   Jem(1,3)*(-2*RMvcs(1,1)*Ur(11)*qsq - 
     9   2*RMvcs(1,2)*Ur(30)*(Q**3) - 2*RMvcs(1,3)*Ur(10)*(Q**4)) + 
     1   Jem(3,3)*(-2*RMvcs(3,1)*Ur(11)*qsq - 
     2   2*RMvcs(3,2)*Ur(30)*(Q**3) - 2*RMvcs(3,3)*Ur(10)*(Q**4))

      SigmaIPol0(3)=Jem(1,2)*(-2*RMvcs(1,1)*Ur(4)*(Q**3) - 
     1   2*RMvcs(1,3)*Ur(41)*(Q**3) + 2*RMvcs(1,2)*Ur(27)*(Q**4)) + 
     2   Jem(3,2)*(-2*RMvcs(3,1)*Ur(4)*(Q**3) - 
     3   2*RMvcs(3,3)*Ur(41)*(Q**3) + 2*RMvcs(3,2)*Ur(27)*(Q**4)) + 
     4   Jem(1,3)*(-2*RMvcs(1,1)*Ur(12)*(Q**3) + 
     5   2*RMvcs(1,2)*Ur(32)*(Q**4) - 2*RMvcs(1,3)*Ur(47)*(Q**5)) + 
     6   Jem(3,3)*(-2*RMvcs(3,1)*Ur(12)*(Q**3) + 
     7   2*RMvcs(3,2)*Ur(32)*(Q**4) - 2*RMvcs(3,3)*Ur(47)*(Q**5)) + 
     8   Jem(2,3)*(-2*RMvcs(2,1)*Ur(71)*(Q**3) + 
     9   2*RMvcs(2,2)*Ur(86)*(Q**4) - 2*RMvcs(2,3)*Ur(98)*(Q**5)) + 
     1   Jem(4,3)*(-2*RMvcs(4,1)*Ur(71)*(Q**3) + 
     2   2*RMvcs(4,2)*Ur(86)*(Q**4) - 2*RMvcs(4,3)*Ur(98)*(Q**5))

      SigmaIPol0(4)=Jem(1,2)*(-2*RMvcs(1,3)*Ur(43)*(Q**3) - 
     1   2*RMvcs(1,1)*Ur(6)*(Q**3) + 2*RMvcs(1,2)*Ur(28)*(Q**4)) + 
     2   Jem(3,2)*(-2*RMvcs(3,3)*Ur(43)*(Q**3) - 
     3   2*RMvcs(3,1)*Ur(6)*(Q**3) + 2*RMvcs(3,2)*Ur(28)*(Q**4)) + 
     4   Jem(1,3)*(-2*RMvcs(1,1)*Ur(14)*(Q**3) + 
     5   2*RMvcs(1,2)*Ur(33)*(Q**4) - 2*RMvcs(1,3)*Ur(49)*(Q**5)) + 
     6   Jem(3,3)*(-2*RMvcs(3,1)*Ur(14)*(Q**3) + 
     7   2*RMvcs(3,2)*Ur(33)*(Q**4) - 2*RMvcs(3,3)*Ur(49)*(Q**5)) + 
     8   Jem(2,3)*(-2*RMvcs(2,1)*Ur(72)*(Q**3) - 
     9   2*RMvcs(2,2)*Ur(86)*(Q**4) - 2*RMvcs(2,3)*Ur(99)*(Q**5)) + 
     1   Jem(4,3)*(-2*RMvcs(4,1)*Ur(72)*(Q**3) - 
     2   2*RMvcs(4,2)*Ur(86)*(Q**4) - 2*RMvcs(4,3)*Ur(99)*(Q**5))

      SigmaIPol0(5)=Jem(1,2)*(-2*RMvcs(1,3)*Ur(3)*qsq + 
     1   2*RMvcs(1,2)*Ur(26)*(Q**3) - 2*RMvcs(1,1)*Ur(2)*(Q**4)) + 
     2   Jem(3,2)*(-2*RMvcs(3,3)*Ur(3)*qsq + 
     3   2*RMvcs(3,2)*Ur(26)*(Q**3) - 2*RMvcs(3,1)*Ur(2)*(Q**4)) + 
     4   Jem(1,3)*(-2*RMvcs(1,3)*Ur(11)*qsq - 
     5   2*RMvcs(1,1)*Ur(10)*(Q**4) + 2*RMvcs(1,2)*Ur(31)*(Q**5)) + 
     6   Jem(3,3)*(-2*RMvcs(3,3)*Ur(11)*qsq - 
     7   2*RMvcs(3,1)*Ur(10)*(Q**4) + 2*RMvcs(3,2)*Ur(31)*(Q**5)) + 
     8   Jem(2,3)*(-2*RMvcs(2,3)*Ur(52)*qsq + 
     9   2*RMvcs(2,1)*Ur(17)*(Q**4) + 2*RMvcs(2,2)*Ur(85)*(Q**5)) + 
     1   Jem(4,3)*(-2*RMvcs(4,3)*Ur(52)*qsq + 
     2   2*RMvcs(4,1)*Ur(17)*(Q**4) + 2*RMvcs(4,2)*Ur(85)*(Q**5))

      SigmaIPol0(6)=Jem(1,2)*(-2*RMvcs(1,3)*Ur(42)*Q - 
     1   2*RMvcs(1,1)*Ur(5)*(Q**3)) + 
     2   Jem(3,2)*(-2*RMvcs(3,3)*Ur(42)*Q - 
     3   2*RMvcs(3,1)*Ur(5)*(Q**3)) + 
     4   Jem(1,3)*(-2*RMvcs(1,3)*Ur(48)*Q - 
     5   2*RMvcs(1,1)*Ur(13)*(Q**5)) + 
     6   Jem(3,3)*(-2*RMvcs(3,3)*Ur(48)*Q - 
     7   2*RMvcs(3,1)*Ur(13)*(Q**5)) + 
     8   Jem(2,3)*(-2*RMvcs(2,3)*Ur(54)*Q + 
     9   2*RMvcs(2,1)*Ur(19)*(Q**5)) + 
     1   Jem(4,3)*(-2*RMvcs(4,3)*Ur(54)*Q + 
     2   2*RMvcs(4,1)*Ur(19)*(Q**5))

      SigmaIPol0(7)=Jem(1,2)*(-2*IMvcs(1,2)*Ur(29)*qsq + 
     1   2*IMvcs(1,3)*Ur(45)*(Q**3) + 2*IMvcs(1,1)*Ur(8)*(Q**3)) + 
     2   Jem(3,2)*(-2*IMvcs(3,2)*Ur(29)*qsq + 
     3   2*IMvcs(3,3)*Ur(45)*(Q**3) + 2*IMvcs(3,1)*Ur(8)*(Q**3)) + 
     4   Jem(1,3)*(2*IMvcs(1,1)*Ur(16)*(Q**3) - 
     5   2*IMvcs(1,2)*Ur(34)*(Q**4) + 2*IMvcs(1,3)*Ur(51)*(Q**5)) + 
     6   Jem(3,3)*(2*IMvcs(3,1)*Ur(16)*(Q**3) - 
     7   2*IMvcs(3,2)*Ur(34)*(Q**4) + 2*IMvcs(3,3)*Ur(51)*(Q**5)) + 
     8   Jem(2,3)*(2*IMvcs(2,1)*Ur(73)*(Q**3) + 
     9   2*IMvcs(2,2)*Ur(86)*(Q**4) + 2*IMvcs(2,3)*Ur(100)*(Q**5)) + 
     1   Jem(4,3)*(2*IMvcs(4,1)*Ur(73)*(Q**3) + 
     2   2*IMvcs(4,2)*Ur(86)*(Q**4) + 2*IMvcs(4,3)*Ur(100)*(Q**5))

      SigmaIPol0(8)=Jem(1,2)*(2*IMvcs(1,3)*Ur(44)*qsq - 
     1   4*IMvcs(1,2)*Ur(26)*(Q**3) + 2*IMvcs(1,1)*Ur(7)*(Q**4)) + 
     2   Jem(3,2)*(2*IMvcs(3,3)*Ur(44)*qsq - 
     3   4*IMvcs(3,2)*Ur(26)*(Q**3) + 2*IMvcs(3,1)*Ur(7)*(Q**4)) + 
     4   Jem(1,3)*(2*IMvcs(1,3)*Ur(50)*qsq + 
     5   2*IMvcs(1,1)*Ur(15)*(Q**4) - 4*IMvcs(1,2)*Ur(31)*(Q**5)) + 
     6   Jem(3,3)*(2*IMvcs(3,3)*Ur(50)*qsq + 
     7   2*IMvcs(3,1)*Ur(15)*(Q**4) - 4*IMvcs(3,2)*Ur(31)*(Q**5)) + 
     8   Jem(2,3)*(2*IMvcs(2,3)*Ur(57)*qsq - 
     9   2*IMvcs(2,1)*Ur(22)*(Q**4) - 4*IMvcs(2,2)*Ur(85)*(Q**5)) + 
     1   Jem(4,3)*(2*IMvcs(4,3)*Ur(57)*qsq - 
     2   2*IMvcs(4,1)*Ur(22)*(Q**4) - 4*IMvcs(4,2)*Ur(85)*(Q**5))

      SigmaIPolX(1)=Jem(1,2)*(2*IMvcs(2,1)*Ur(60)*(Q**3) + 
     1    2*IMvcs(2,3)*Ur(87)*(Q**3) - 2*IMvcs(2,2)*Ur(75)*(Q**4)) + 
     2    Jem(3,2)*(-2*IMvcs(4,1)*Ur(60)*(Q**3) - 
     3    2*IMvcs(4,3)*Ur(87)*(Q**3) + 2*IMvcs(4,2)*Ur(75)*(Q**4)) + 
     4    Jem(2,3)*(2*IMvcs(1,1)*Ur(18)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(36)*(Q**4) + 2*IMvcs(1,3)*Ur(53)*(Q**5)) + 
     6    Jem(4,3)*(-2*IMvcs(3,1)*Ur(18)*(Q**3) + 
     7    2*IMvcs(3,2)*Ur(36)*(Q**4) - 2*IMvcs(3,3)*Ur(53)*(Q**5)) + 
     8    Jem(1,3)*(2*IMvcs(2,1)*Ur(65)*(Q**3) - 
     9    2*IMvcs(2,2)*Ur(80)*(Q**4) + 2*IMvcs(2,3)*Ur(92)*(Q**5)) + 
     1    Jem(3,3)*(-2*IMvcs(4,1)*Ur(65)*(Q**3) + 
     2    2*IMvcs(4,2)*Ur(80)*(Q**4) - 2*IMvcs(4,3)*Ur(92)*(Q**5))

      SigmaIPolX(2)=Jem(1,2)*(2*IMvcs(2,1)*Ur(61)*(Q**3) + 
     1    2*IMvcs(2,3)*Ur(88)*(Q**3) - 2*IMvcs(2,2)*Ur(76)*(Q**4)) + 
     2    Jem(3,2)*(-2*IMvcs(4,1)*Ur(61)*(Q**3) - 
     3    2*IMvcs(4,3)*Ur(88)*(Q**3) + 2*IMvcs(4,2)*Ur(76)*(Q**4)) + 
     4    Jem(2,3)*(2*IMvcs(1,1)*Ur(20)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(37)*(Q**4) + 2*IMvcs(1,3)*Ur(55)*(Q**5)) + 
     6    Jem(4,3)*(-2*IMvcs(3,1)*Ur(20)*(Q**3) + 
     7    2*IMvcs(3,2)*Ur(37)*(Q**4) - 2*IMvcs(3,3)*Ur(55)*(Q**5)) + 
     8    Jem(1,3)*(2*IMvcs(2,1)*Ur(66)*(Q**3) - 
     9    2*IMvcs(2,2)*Ur(81)*(Q**4) + 2*IMvcs(2,3)*Ur(93)*(Q**5)) + 
     1    Jem(3,3)*(-2*IMvcs(4,1)*Ur(66)*(Q**3) + 
     2    2*IMvcs(4,2)*Ur(81)*(Q**4) - 2*IMvcs(4,3)*Ur(93)*(Q**5))

      SigmaIPolX(3)=Jem(1,2)*(2*IMvcs(2,3)*Ur(3)*qsq - 
     1    2*IMvcs(2,2)*Ur(74)*(Q**3) - 2*IMvcs(2,1)*Ur(2)*(Q**4)) + 
     2    Jem(3,2)*(-2*IMvcs(4,3)*Ur(3)*qsq + 
     3    2*IMvcs(4,2)*Ur(74)*(Q**3) + 2*IMvcs(4,1)*Ur(2)*(Q**4)) + 
     4    Jem(2,3)*(2*IMvcs(1,3)*Ur(52)*qsq + 
     5    2*IMvcs(1,1)*Ur(17)*(Q**4) - 2*IMvcs(1,2)*Ur(35)*(Q**5)) + 
     6    Jem(4,3)*(-2*IMvcs(3,3)*Ur(52)*qsq - 
     7    2*IMvcs(3,1)*Ur(17)*(Q**4) + 2*IMvcs(3,2)*Ur(35)*(Q**5)) + 
     8    Jem(1,3)*(2*IMvcs(2,3)*Ur(11)*qsq - 
     9    2*IMvcs(2,1)*Ur(10)*(Q**4) - 2*IMvcs(2,2)*Ur(79)*(Q**5)) + 
     1    Jem(3,3)*(-2*IMvcs(4,3)*Ur(11)*qsq + 
     2    2*IMvcs(4,1)*Ur(10)*(Q**4) + 2*IMvcs(4,2)*Ur(79)*(Q**5))

      SigmaIPolX(4)=Jem(1,2)*(2*IMvcs(2,3)*Ur(42)*Q - 
     1    2*IMvcs(2,1)*Ur(5)*(Q**3)) + 
     2    Jem(3,2)*(-2*IMvcs(4,3)*Ur(42)*Q + 
     3    2*IMvcs(4,1)*Ur(5)*(Q**3)) + 
     4    Jem(1,3)*(2*IMvcs(2,3)*Ur(48)*Q - 
     5    2*IMvcs(2,1)*Ur(13)*(Q**5)) + 
     6    Jem(3,3)*(-2*IMvcs(4,3)*Ur(48)*Q + 
     7    2*IMvcs(4,1)*Ur(13)*(Q**5)) + 
     8    Jem(2,3)*(2*IMvcs(1,3)*Ur(54)*Q + 
     9    2*IMvcs(1,1)*Ur(19)*(Q**5)) + 
     1    Jem(4,3)*(-2*IMvcs(3,3)*Ur(54)*Q - 
     2    2*IMvcs(3,1)*Ur(19)*(Q**5))

      SigmaIPolX(5)=Jem(2,3)*(-2*RMvcs(1,1)*Ur(11)*qsq + 
     1    2*RMvcs(1,2)*Ur(38)*(Q**3) - 2*RMvcs(1,3)*Ur(56)*(Q**4)) + 
     2    Jem(4,3)*(2*RMvcs(3,1)*Ur(21)*qsq - 
     3    2*RMvcs(3,2)*Ur(38)*(Q**3) + 2*RMvcs(3,3)*Ur(56)*(Q**4)) + 
     4    Jem(1,2)*(-2*RMvcs(2,1)*Ur(62)*qsq + 
     5    2*RMvcs(2,2)*Ur(77)*(Q**3) - 2*RMvcs(2,3)*Ur(89)*(Q**4)) + 
     6    Jem(3,2)*(2*RMvcs(4,1)*Ur(62)*qsq - 
     7    2*RMvcs(4,2)*Ur(77)*(Q**3) + 2*RMvcs(4,3)*Ur(89)*(Q**4)) + 
     8    Jem(1,3)*(-2*RMvcs(2,1)*Ur(67)*qsq + 
     9    2*RMvcs(2,2)*Ur(82)*(Q**3) - 2*RMvcs(2,3)*Ur(94)*(Q**4)) + 
     1    Jem(3,3)*(2*RMvcs(4,1)*Ur(67)*qsq - 
     2    2*RMvcs(4,2)*Ur(82)*(Q**3) + 2*RMvcs(4,3)*Ur(94)*(Q**4))

      SigmaIPolX(6)=Jem(2,3)*(-2*RMvcs(1,1)*Ur(23)*qsq - 
     1    2*RMvcs(1,2)*Ur(38)*(Q**3) - 2*RMvcs(1,3)*Ur(58)*(Q**4)) + 
     2    Jem(4,3)*(2*RMvcs(3,1)*Ur(23)*qsq + 
     3    2*RMvcs(3,2)*Ur(38)*(Q**3) + 2*RMvcs(3,3)*Ur(58)*(Q**4)) + 
     4    Jem(1,2)*(-2*RMvcs(2,1)*Ur(63)*qsq - 
     5    2*RMvcs(2,2)*Ur(77)*(Q**3) - 2*RMvcs(2,3)*Ur(90)*(Q**4)) + 
     6    Jem(3,2)*(2*RMvcs(4,1)*Ur(63)*qsq + 
     7    2*RMvcs(4,2)*Ur(77)*(Q**3) + 2*RMvcs(4,3)*Ur(90)*(Q**4)) + 
     8    Jem(1,3)*(-2*RMvcs(2,1)*Ur(68)*qsq - 
     9    2*RMvcs(2,2)*Ur(82)*(Q**3) - 2*RMvcs(2,3)*Ur(95)*(Q**4)) + 
     1    Jem(3,3)*(2*RMvcs(4,1)*Ur(68)*qsq + 
     2    2*RMvcs(4,2)*Ur(82)*(Q**3) + 2*RMvcs(4,3)*Ur(95)*(Q**4))

      SigmaIPolX(7)=Jem(1,2)*(2*RMvcs(2,2)*Ur(78)*qsq - 
     1    2*RMvcs(2,1)*Ur(64)*(Q**3) - 2*RMvcs(2,3)*Ur(91)*(Q**3)) + 
     2    Jem(3,2)*(-2*RMvcs(4,2)*Ur(78)*qsq + 
     3    2*RMvcs(4,1)*Ur(64)*(Q**3) + 2*RMvcs(4,3)*Ur(91)*(Q**3)) + 
     4    Jem(2,3)*(-2*RMvcs(1,1)*Ur(24)*(Q**3) + 
     5    2*RMvcs(1,2)*Ur(39)*(Q**4) - 2*RMvcs(1,3)*Ur(59)*(Q**5)) + 
     6    Jem(4,3)*(2*RMvcs(3,1)*Ur(24)*(Q**3) - 
     7    2*RMvcs(3,2)*Ur(39)*(Q**4) + 2*RMvcs(3,3)*Ur(59)*(Q**5)) + 
     8    Jem(1,3)*(-2*RMvcs(2,1)*Ur(69)*(Q**3) + 
     9    2*RMvcs(2,2)*Ur(83)*(Q**4) - 2*RMvcs(2,3)*Ur(96)*(Q**5)) + 
     1    Jem(3,3)*(2*RMvcs(4,1)*Ur(69)*(Q**3) - 
     2    2*RMvcs(4,2)*Ur(83)*(Q**4) + 2*RMvcs(4,3)*Ur(96)*(Q**5))

      SigmaIPolX(8)=Jem(1,2)*(-2*RMvcs(2,3)*Ur(44)*qsq + 
     1    4*RMvcs(2,2)*Ur(74)*(Q**3) + 2*RMvcs(2,1)*Ur(7)*(Q**4)) + 
     2    Jem(3,2)*(2*RMvcs(4,3)*Ur(44)*qsq - 
     3    4*RMvcs(4,2)*Ur(74)*(Q**3) - 2*RMvcs(4,1)*Ur(7)*(Q**4)) + 
     4    Jem(2,3)*(-2*RMvcs(1,3)*Ur(57)*qsq - 
     5    2*RMvcs(1,1)*Ur(22)*(Q**4) + 4*RMvcs(1,2)*Ur(35)*(Q**5)) + 
     6    Jem(4,3)*(2*RMvcs(3,3)*Ur(57)*qsq + 
     7    2*RMvcs(3,1)*Ur(22)*(Q**4) - 4*RMvcs(3,2)*Ur(35)*(Q**5)) + 
     8    Jem(1,3)*(-2*RMvcs(2,3)*Ur(50)*qsq + 
     9    2*RMvcs(2,1)*Ur(15)*(Q**4) + 4*RMvcs(2,2)*Ur(79)*(Q**5)) + 
     1    Jem(3,3)*(2*RMvcs(4,3)*Ur(50)*qsq - 
     2    2*RMvcs(4,1)*Ur(15)*(Q**4) - 4*RMvcs(4,2)*Ur(79)*(Q**5))

      SigmaIPolY(1)=Jem(3,2)*(2*IMvcs(1,1)*Ur(1)*qsq - 
     1    2*IMvcs(1,2)*Ur(25)*(Q**3) + 2*IMvcs(1,3)*Ur(40)*(Q**4)) + 
     2    Jem(1,2)*(-2*IMvcs(3,1)*Ur(1)*qsq + 
     3    2*IMvcs(3,2)*Ur(25)*(Q**3) - 2*IMvcs(3,3)*Ur(40)*(Q**4)) + 
     4    Jem(3,3)*(2*IMvcs(1,1)*Ur(9)*qsq - 
     5    2*IMvcs(1,2)*Ur(30)*(Q**3) + 2*IMvcs(1,3)*Ur(46)*(Q**4)) + 
     6    Jem(1,3)*(-2*IMvcs(3,1)*Ur(9)*qsq + 
     7    2*IMvcs(3,2)*Ur(30)*(Q**3) - 2*IMvcs(3,3)*Ur(46)*(Q**4)) + 
     8    Jem(4,3)*(2*IMvcs(2,1)*Ur(70)*qsq - 
     9    2*IMvcs(2,2)*Ur(84)*(Q**3) + 2*IMvcs(2,3)*Ur(97)*(Q**4)) + 
     1    Jem(2,3)*(-2*IMvcs(4,1)*Ur(70)*qsq + 
     2    2*IMvcs(4,2)*Ur(84)*(Q**3) - 2*IMvcs(4,3)*Ur(97)*(Q**4))

      SigmaIPolY(2)=Jem(3,2)*(2*IMvcs(1,1)*Ur(3)*qsq + 
     1    2*IMvcs(1,2)*Ur(25)*(Q**3) + 2*IMvcs(1,3)*Ur(2)*(Q**4)) + 
     2    Jem(1,2)*(-2*IMvcs(3,1)*Ur(3)*qsq - 
     3    2*IMvcs(3,2)*Ur(25)*(Q**3) - 2*IMvcs(3,3)*Ur(2)*(Q**4)) + 
     4    Jem(4,3)*(2*IMvcs(2,1)*Ur(52)*qsq + 
     5    2*IMvcs(2,2)*Ur(84)*(Q**3) - 2*IMvcs(2,3)*Ur(17)*(Q**4)) + 
     6    Jem(2,3)*(-2*IMvcs(4,1)*Ur(52)*qsq - 
     7    2*IMvcs(4,2)*Ur(84)*(Q**3) + 2*IMvcs(4,3)*Ur(17)*(Q**4)) + 
     8    Jem(3,3)*(2*IMvcs(1,1)*Ur(11)*qsq + 
     9    2*IMvcs(1,2)*Ur(30)*(Q**3) + 2*IMvcs(1,3)*Ur(10)*(Q**4)) + 
     1    Jem(1,3)*(-2*IMvcs(3,1)*Ur(11)*qsq - 
     2    2*IMvcs(3,2)*Ur(30)*(Q**3) - 2*IMvcs(3,3)*Ur(10)*(Q**4))

      SigmaIPolY(3)=Jem(3,2)*(2*IMvcs(1,1)*Ur(4)*(Q**3) + 
     1    2*IMvcs(1,3)*Ur(41)*(Q**3) - 2*IMvcs(1,2)*Ur(27)*(Q**4)) + 
     2    Jem(1,2)*(-2*IMvcs(3,1)*Ur(4)*(Q**3) - 
     3    2*IMvcs(3,3)*Ur(41)*(Q**3) + 2*IMvcs(3,2)*Ur(27)*(Q**4)) + 
     4    Jem(3,3)*(2*IMvcs(1,1)*Ur(12)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(32)*(Q**4) + 2*IMvcs(1,3)*Ur(47)*(Q**5)) + 
     6    Jem(1,3)*(-2*IMvcs(3,1)*Ur(12)*(Q**3) + 
     7    2*IMvcs(3,2)*Ur(32)*(Q**4) - 2*IMvcs(3,3)*Ur(47)*(Q**5)) + 
     8    Jem(4,3)*(2*IMvcs(2,1)*Ur(71)*(Q**3) - 
     9    2*IMvcs(2,2)*Ur(86)*(Q**4) + 2*IMvcs(2,3)*Ur(98)*(Q**5)) + 
     1    Jem(2,3)*(-2*IMvcs(4,1)*Ur(71)*(Q**3) + 
     2    2*IMvcs(4,2)*Ur(86)*(Q**4) - 2*IMvcs(4,3)*Ur(98)*(Q**5))

      SigmaIPolY(4)=Jem(3,2)*(2*IMvcs(1,3)*Ur(43)*(Q**3) + 
     1    2*IMvcs(1,1)*Ur(6)*(Q**3) - 2*IMvcs(1,2)*Ur(28)*(Q**4)) + 
     2    Jem(1,2)*(-2*IMvcs(3,3)*Ur(43)*(Q**3) - 
     3    2*IMvcs(3,1)*Ur(6)*(Q**3) + 2*IMvcs(3,2)*Ur(28)*(Q**4)) + 
     4    Jem(3,3)*(2*IMvcs(1,1)*Ur(14)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(33)*(Q**4) + 2*IMvcs(1,3)*Ur(49)*(Q**5)) + 
     6    Jem(1,3)*(-2*IMvcs(3,1)*Ur(14)*(Q**4) + 
     7    2*IMvcs(3,2)*Ur(33)*(Q**4) - 2*IMvcs(3,3)*Ur(49)*(Q**5)) + 
     8    Jem(4,3)*(2*IMvcs(2,1)*Ur(72)*(Q**3) + 
     9    2*IMvcs(2,2)*Ur(86)*(Q**4) + 2*IMvcs(2,3)*Ur(99)*(Q**5)) + 
     1    Jem(2,3)*(-2*IMvcs(4,1)*Ur(72)*(Q**3) - 
     2    2*IMvcs(4,2)*Ur(86)*(Q**4) - 2*IMvcs(4,3)*Ur(99)*(Q**5))

      SigmaIPolY(5)=Jem(3,2)*(2*IMvcs(1,3)*Ur(3)*qsq - 
     1    2*IMvcs(1,2)*Ur(26)*(Q**3) + 2*IMvcs(1,1)*Ur(2)*(Q**4)) + 
     2    Jem(1,2)*(-2*IMvcs(3,3)*Ur(3)*qsq + 
     3    2*IMvcs(3,2)*Ur(26)*(Q**3) - 2*IMvcs(3,1)*Ur(2)*(Q**4)) + 
     4    Jem(3,3)*(2*IMvcs(1,3)*Ur(11)*qsq + 
     5    2*IMvcs(1,1)*Ur(10)*(Q**4) - 2*IMvcs(1,2)*Ur(31)*(Q**5)) + 
     6    Jem(1,3)*(-2*IMvcs(3,3)*Ur(11)*qsq - 
     7    2*IMvcs(3,1)*Ur(10)*(Q**4) + 2*IMvcs(3,2)*Ur(31)*(Q**5)) + 
     8    Jem(4,3)*(2*IMvcs(2,3)*Ur(52)*qsq - 
     9    2*IMvcs(2,1)*Ur(17)*(Q**4) - 2*IMvcs(2,2)*Ur(85)*(Q**5)) + 
     1    Jem(2,3)*(-2*IMvcs(4,3)*Ur(52)*qsq + 
     2    2*IMvcs(4,1)*Ur(17)*(Q**4) + 2*IMvcs(4,2)*Ur(85)*(Q**5))

      SigmaIPolY(6)=Jem(3,2)*(2*IMvcs(1,3)*Ur(42)*Q + 
     1    2*IMvcs(1,1)*Ur(5)*(Q**3)) + 
     2    Jem(1,2)*(-2*IMvcs(3,3)*Ur(42)*Q - 
     3    2*IMvcs(3,1)*Ur(5)*(Q**3)) + 
     4    Jem(3,3)*(2*IMvcs(1,3)*Ur(48)*Q + 
     5    2*IMvcs(1,1)*Ur(13)*(Q**5)) + 
     6    Jem(1,3)*(-2*IMvcs(3,3)*Ur(48)*Q - 
     7    2*IMvcs(3,1)*Ur(13)*(Q**5)) + 
     8    Jem(4,3)*(2*IMvcs(2,3)*Ur(54)*Q - 
     9    2*IMvcs(2,1)*Ur(19)*(Q**5)) + 
     1    Jem(2,3)*(-2*IMvcs(4,3)*Ur(54)*Q + 
     2    2*IMvcs(4,1)*Ur(19)*(Q**5))

      SigmaIPolY(7)=Jem(3,2)*(-2*RMvcs(1,2)*Ur(29)*qsq + 
     1    2*RMvcs(1,3)*Ur(45)*(Q**3) + 2*RMvcs(1,1)*Ur(8)*(Q**3)) + 
     2    Jem(1,2)*(2*RMvcs(3,2)*Ur(29)*qsq - 
     3    2*RMvcs(3,3)*Ur(45)*(Q**3) - 2*RMvcs(3,1)*Ur(8)*(Q**3)) + 
     4    Jem(3,3)*(2*RMvcs(1,1)*Ur(16)*(Q**3) - 
     5    2*RMvcs(1,2)*Ur(34)*(qsq**2) + 2*RMvcs(1,3)*Ur(51)*(Q**5)) + 
     6    Jem(1,3)*(-2*RMvcs(3,1)*Ur(16)*(Q**3) + 
     7    2*RMvcs(3,2)*Ur(34)*(qsq**2) - 2*RMvcs(3,3)*Ur(51)*(Q**5)) + 
     8    Jem(4,3)*(2*RMvcs(2,1)*Ur(73)*(Q**3) + 
     9    2*RMvcs(2,2)*Ur(86)*(qsq**2) + 2*RMvcs(2,3)*Ur(100)*(Q**5)) + 
     1    Jem(2,3)*(-2*RMvcs(4,1)*Ur(73)*(Q**3) - 
     2    2*RMvcs(4,2)*Ur(86)*(qsq**2) - 2*RMvcs(4,3)*Ur(100)*(Q**5))

      SigmaIPolY(8)=Jem(3,2)*(2*RMvcs(1,3)*Ur(44)*qsq - 
     1    4*RMvcs(1,2)*Ur(26)*(Q**3) + 2*RMvcs(1,1)*Ur(7)*(qsq**2)) + 
     2    Jem(1,2)*(-2*RMvcs(3,3)*Ur(44)*qsq + 
     3    4*RMvcs(3,2)*Ur(26)*(Q**3) - 2*RMvcs(3,1)*Ur(7)*(qsq**2)) + 
     4    Jem(3,3)*(2*RMvcs(1,3)*Ur(50)*qsq + 2*RMvcs(1,1)*Ur(15)*(qsq**2) - 
     5    4*RMvcs(1,2)*Ur(31)*(Q**5)) + 
     6    Jem(1,3)*(-2*RMvcs(3,3)*Ur(50)*qsq - 
     7    2*RMvcs(3,1)*Ur(15)*(qsq**2) + 4*RMvcs(3,2)*Ur(41)*(Q**5)) + 
     8    Jem(4,3)*(2*RMvcs(2,3)*Ur(57)*qsq - 
     9    2*RMvcs(2,1)*Ur(22)*(qsq**2) - 4*RMvcs(2,2)*Ur(85)*(Q**5)) + 
     1    Jem(2,3)*(-2*RMvcs(4,3)*Ur(57)*qsq + 
     2    2*RMvcs(4,1)*Ur(22)*(qsq**2) + 4*RMvcs(4,2)*Ur(85)*(Q**5))

      SigmaIPolZ(1)=Jem(3,2)*(2*IMvcs(2,1)*Ur(60)*(Q**3) + 
     1    2*IMvcs(2,3)*Ur(87)*(Q**3) - 2*IMvcs(2,2)*Ur(75)*(qsq**2)) + 
     2    Jem(1,2)*(2*IMvcs(4,1)*Ur(60)*(Q**3) + 
     3    2*IMvcs(4,3)*Ur(87)*(Q**3) - 2*IMvcs(4,2)*Ur(75)*(qsq**2)) + 
     4    Jem(4,3)*(2*IMvcs(1,1)*Ur(18)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(36)*(qsq**2) + 2*IMvcs(1,3)*Ur(53)*(Q**5)) + 
     6    Jem(2,3)*(2*IMvcs(3,1)*Ur(18)*(Q**3) - 
     7    2*IMvcs(3,2)*Ur(36)*(qsq**2) + 2*IMvcs(3,3)*Ur(53)*(Q**5)) + 
     8    Jem(3,3)*(2*IMvcs(2,1)*Ur(65)*(Q**3) - 
     9    2*IMvcs(2,2)*Ur(80)*(qsq**2) + 2*IMvcs(2,3)*Ur(92)*(Q**5)) + 
     1    Jem(1,3)*(2*IMvcs(4,1)*Ur(65)*(Q**3) - 
     2    2*IMvcs(4,2)*Ur(80)*(qsq**2) + 2*IMvcs(4,3)*Ur(92)*(Q**5))

      SigmaIPolZ(2)=Jem(3,2)*(2*IMvcs(2,1)*Ur(61)*(Q**3) + 
     1    2*IMvcs(2,3)*Ur(88)*(Q**3) - 2*IMvcs(2,2)*Ur(76)*(qsq**2)) + 
     2    Jem(1,2)*(2*IMvcs(4,1)*Ur(61)*(Q**3) + 
     3    2*IMvcs(4,3)*Ur(88)*(Q**3) - 2*IMvcs(4,2)*Ur(76)*(qsq**2)) + 
     4    Jem(4,3)*(2*IMvcs(1,1)*Ur(20)*(Q**3) - 
     5    2*IMvcs(1,2)*Ur(37)*(qsq**2) + 2*IMvcs(1,3)*Ur(55)*(Q**5)) + 
     6    Jem(2,3)*(2*IMvcs(3,1)*Ur(20)*(Q**3) - 
     7    2*IMvcs(3,2)*Ur(37)*(qsq**2) + 2*IMvcs(3,3)*Ur(55)*(Q**5)) + 
     8    Jem(3,3)*(2*IMvcs(2,1)*Ur(66)*(Q**3) - 
     9    2*IMvcs(2,2)*Ur(81)*(qsq**2) + 2*IMvcs(2,3)*Ur(93)*(Q**5)) + 
     1    Jem(1,3)*(2*IMvcs(4,1)*Ur(66)*(Q**3) - 
     2    2*IMvcs(4,2)*Ur(81)*(qsq**2) + 2*IMvcs(4,3)*Ur(93)*(Q**5))

      SigmaIPolZ(3)=Jem(3,2)*(2*IMvcs(2,3)*Ur(3)*qsq - 
     1    2*IMvcs(2,2)*Ur(74)*(Q**3) - 2*IMvcs(2,1)*Ur(2)*(qsq**2)) + 
     2    Jem(1,2)*(2*IMvcs(4,3)*Ur(3)*qsq - 
     3    2*IMvcs(4,2)*Ur(74)*(Q**3) - 2*IMvcs(4,1)*Ur(2)*(qsq**2)) + 
     4    Jem(4,3)*(2*IMvcs(1,3)*Ur(52)*qsq + 
     5    2*IMvcs(1,1)*Ur(17)*(qsq**2) - 2*IMvcs(1,2)*Ur(35)*(Q**5)) + 
     6    Jem(2,3)*(2*IMvcs(3,3)*Ur(52)*qsq + 
     7    2*IMvcs(3,1)*Ur(17)*(qsq**2) - 2*IMvcs(3,2)*Ur(35)*(Q**5)) + 
     8    Jem(3,3)*(2*IMvcs(2,3)*Ur(11)*qsq - 
     9    2*IMvcs(2,1)*Ur(10)*(qsq**2) - 2*IMvcs(2,2)*Ur(79)*(Q**5)) + 
     1    Jem(1,3)*(2*IMvcs(4,3)*Ur(11)*qsq - 
     2    2*IMvcs(4,1)*Ur(10)*(qsq**2) - 2*IMvcs(4,2)*Ur(79)*(Q**5))

      SigmaIPolZ(4)=Jem(3,2)*(2*IMvcs(2,3)*Ur(42)*Q - 
     1    2*IMvcs(2,1)*Ur(5)*(Q**3)) + 
     2    Jem(1,2)*(2*IMvcs(4,3)*Ur(42)*Q - 
     3    2*IMvcs(4,1)*Ur(5)*(Q**3)) + 
     4    Jem(3,3)*(2*IMvcs(2,3)*Ur(48)*Q - 
     5    2*IMvcs(2,1)*Ur(13)*(Q**5)) + 
     6    Jem(1,3)*(2*IMvcs(4,3)*Ur(48)*Q - 
     7    2*IMvcs(4,1)*Ur(13)*(Q**5)) + 
     8    Jem(4,3)*(2*IMvcs(1,3)*Ur(54)*Q + 
     9    2*IMvcs(1,1)*Ur(19)*(Q**5)) + 
     1    Jem(2,3)*(2*IMvcs(3,3)*Ur(54)*Q + 
     2    2*IMvcs(3,1)*Ur(19)*(Q**5))

      SigmaIPolZ(5)=Jem(4,3)*(-2*RMvcs(1,1)*Ur(21)*qsq + 
     1    2*RMvcs(1,2)*Ur(38)*(Q**3) - 2*RMvcs(1,3)*Ur(56)*(qsq**2)) + 
     2    Jem(2,3)*(-2*RMvcs(3,1)*Ur(21)*qsq + 
     3    2*RMvcs(3,2)*Ur(38)*(Q**3) - 2*RMvcs(3,3)*Ur(56)*(qsq**2)) + 
     4    Jem(3,2)*(-2*RMvcs(2,1)*Ur(62)*qsq + 
     5    2*RMvcs(2,2)*Ur(77)*(Q**3) - 2*RMvcs(2,3)*Ur(89)*(qsq**2)) + 
     6    Jem(1,2)*(-2*RMvcs(4,1)*Ur(62)*qsq + 
     7    2*RMvcs(4,2)*Ur(77)*(Q**3) - 2*RMvcs(4,3)*Ur(89)*(qsq**2)) + 
     8    Jem(3,3)*(-2*RMvcs(2,1)*Ur(67)*qsq + 
     9    2*RMvcs(2,2)*Ur(82)*(Q**3) - 2*RMvcs(2,3)*Ur(94)*(qsq**2)) + 
     1    Jem(1,3)*(-2*RMvcs(4,1)*Ur(67)*qsq + 
     2    2*RMvcs(4,2)*Ur(82)*(Q**3) - 2*RMvcs(4,3)*Ur(94)*(qsq**2))

      SigmaIPolZ(6)=Jem(4,3)*(-2*RMvcs(1,1)*Ur(23)*qsq - 
     1    2*RMvcs(1,2)*Ur(38)*(Q**3) - 2*RMvcs(1,3)*Ur(58)*(qsq**2)) + 
     2    Jem(2,3)*(-2*RMvcs(3,1)*Ur(23)*qsq - 
     3    2*RMvcs(3,2)*Ur(38)*(Q**3) - 2*RMvcs(3,3)*Ur(58)*(qsq**2)) + 
     4    Jem(3,2)*(-2*RMvcs(2,1)*Ur(63)*qsq - 
     5    2*RMvcs(2,2)*Ur(77)*(Q**3) - 2*RMvcs(2,3)*Ur(90)*(qsq**2)) + 
     6    Jem(1,2)*(-2*RMvcs(4,1)*Ur(63)*qsq - 
     7    2*RMvcs(4,2)*Ur(77)*(Q**3) - 2*RMvcs(4,3)*Ur(90)*(qsq**2)) + 
     8    Jem(3,3)*(-2*RMvcs(2,1)*Ur(68)*qsq - 
     9    2*RMvcs(2,2)*Ur(82)*(Q**3) - 2*RMvcs(2,3)*Ur(95)*(qsq**2)) + 
     1    Jem(1,3)*(-2*RMvcs(4,1)*Ur(68)*qsq - 
     2    2*RMvcs(4,2)*Ur(82)*(Q**3) - 2*RMvcs(4,3)*Ur(95)*(qsq**2))

      SigmaIPolZ(7)=Jem(3,2)*(2*RMvcs(2,2)*Ur(78)*qsq - 
     1    2*RMvcs(2,1)*Ur(64)*(Q**3) - 2*RMvcs(2,3)*Ur(91)*(Q**3)) + 
     2    Jem(1,2)*(2*RMvcs(4,2)*Ur(78)*qsq - 
     3    2*RMvcs(4,1)*Ur(64)*(Q**3) - 2*RMvcs(4,3)*Ur(91)*(Q**3)) + 
     4    Jem(4,3)*(-2*RMvcs(1,1)*Ur(24)*(Q**3) + 
     5    2*RMvcs(1,2)*Ur(39)*(qsq**2) - 2*RMvcs(1,3)*Ur(59)*(Q**5)) + 
     6    Jem(2,3)*(-2*RMvcs(3,1)*Ur(24)*(Q**3) + 
     7    2*RMvcs(3,2)*Ur(39)*(qsq**2) - 2*RMvcs(3,3)*Ur(59)*(Q**5)) + 
     8    Jem(3,3)*(-2*RMvcs(2,1)*Ur(69)*(Q**3) + 
     9    2*RMvcs(2,2)*Ur(83)*(qsq**2) - 2*RMvcs(2,3)*Ur(96)*(Q**5)) + 
     1    Jem(1,3)*(-2*RMvcs(4,1)*Ur(69)*(Q**3) + 
     2    2*RMvcs(4,2)*Ur(83)*(qsq**2) - 2*RMvcs(4,3)*Ur(96)*(Q**5))

      SigmaIPolZ(8)=Jem(3,2)*(-2*RMvcs(2,3)*Ur(44)*qsq + 
     1    4*RMvcs(2,2)*Ur(74)*(Q**3) + 2*RMvcs(2,1)*Ur(7)*(qsq**2)) + 
     2    Jem(1,2)*(-2*RMvcs(4,3)*Ur(44)*qsq + 
     3    4*RMvcs(4,2)*Ur(74)*(Q**3) + 2*RMvcs(4,1)*Ur(7)*(qsq**2)) + 
     4    Jem(4,3)*(-2*RMvcs(1,3)*Ur(57)*qsq - 
     5    2*RMvcs(1,1)*Ur(22)*(qsq**2) + 4*RMvcs(1,2)*Ur(35)*(Q**5)) + 
     6    Jem(2,3)*(-2*RMvcs(3,3)*Ur(57)*qsq - 
     7    2*RMvcs(3,1)*Ur(22)*(qsq**2) + 4*RMvcs(3,2)*Ur(35)*(Q**5)) + 
     8    Jem(3,3)*(-2*RMvcs(2,3)*Ur(50)*qsq + 
     9    2*RMvcs(2,1)*Ur(15)*(qsq**2) + 4*RMvcs(2,2)*Ur(79)*(Q**5)) + 
     1    Jem(1,3)*(-2*RMvcs(4,3)*Ur(50)*qsq + 
     2    2*RMvcs(4,1)*Ur(15)*(qsq**2) + 4*RMvcs(4,2)*Ur(79)*(Q**5))
    
      if (debug(2)) then

         do j=0,90,10
            write(6,*)j,(Ur(i+j),i=1,10)
         enddo

         write(6,*)' SigmaVCSPol0(1) ',SigmaVCSPol0(1)
         write(6,*)' SigmaVCSPolZ(1) ',SigmaVCSPolZ(1)
         write(6,*)' SigmaVCSPol0(2) ',SigmaVCSPol0(2)
         write(6,*)' SigmaVCSPolX(1) ',SigmaVCSPolX(1)
         write(6,*)' SigmaVCSPolX(2) ',SigmaVCSPolX(2)
         write(6,*)' SigmaVCSPolZ(2) ',SigmaVCSPolZ(2)
         write(6,*)' SigmaVCSPol0(3) ',SigmaVCSPol0(3)
         write(6,*)' SigmaVCSPol0(4) ',SigmaVCSPol0(4)
         write(6,*)' SigmaVCSPolY(1) ',SigmaVCSPolY(1)
         write(6,*)' SigmaVCSPolY(2) ',SigmaVCSPolY(2)
         write(6,*)' SigmaVCSPolY(3) ',SigmaVCSPolY(3)
         write(6,*)' SigmaVCSPolY(4) ',SigmaVCSPolY(4)
         write(6,*)' SigmaVCSPolY(5) ',SigmaVCSPolY(5)
         write(6,*)' SigmaVCSPol0(5) ',SigmaVCSPol0(5)
         write(6,*)' SigmaVCSPolX(3) ',SigmaVCSPolX(3)
         write(6,*)' SigmaVCSPolZ(3) ',SigmaVCSPolZ(3)
         write(6,*)' SigmaVCSPolX(4) ',SigmaVCSPolX(4)
         write(6,*)' SigmaVCSPolZ(4) ',SigmaVCSPolZ(4)
          
         write(6,*)' SigmaIPol0 ',(SigmaIPol0(i),i=1,8) 
         write(6,*)' SigmaIPolX ',(SigmaIPolX(i),i=1,8) 
         write(6,*)' SigmaIPolY ',(SigmaIPolY(i),i=1,8) 
         write(6,*)' SigmaIPolZ ',(SigmaIPolX(i),i=1,8)  

      endif

 500  continue

c----------------------------- DdirectDcrossed(phi) ----------------------------*
c Computes the denominator of the Bethe Heitler cross section.                  |
c Only depends on kinematics.                                                   |
c-------------------------------------------------------------------------------*
c BH is exact calculation, needs actual t_gev and Ge,Bm calc with this value

      DDDC=((qsq + t_gev)**2)/4. 
     1    - ( (qsq*(qsq + t_gev*(-1 + 2*xB))*CosH(Omega))/
     2        (2.*sqrt((qsq**2) + 4*(m_p**2)*qsq*(xB**2)))
     3        - Q*qpPerp*Cos(phi)*SinH(Omega) )**2

c---- SqrAmplBH(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) ----*
c Computes the expansion of the Bethe Heitler cross section (divided by phase   |
c space) in terms of the "expansion cross sections" computed in the function    | 
c MakeExactCrossSections (or MakeLeadingCrossSections).                         |
c-------------------------------------------------------------------------------*
c need to add +1 to all indices, since fortran starts at 1, but c++ starts at 0.

c BH is exact calculation, needs actual t_gev and Ge,Bm calc with this value

      if (TargetPolar.eq.0) then
        SigmaBHp=-(SigmaBHPol0(1) - 
     1    SigmaBHPol0(4)*Cos(2*phi)*(-1 + CosH(2*Omega)) + 
     2    SigmaBHPol0(1)*CosH(2*Omega) 
     3    + SigmaBHPol0(3)*Cos(phi)*SinH(2*Omega) )/(4.*DDDC*(t_gev**2))

        SigmaBHm=SigmaBHp
      endif

      if (TargetPolar.eq.1) then
        BeamHeli=+1
        SigmaBHp=-(SigmaBHPolX(1)*BeamHeli*CosH(Omega) + 
     1    SigmaBHPolX(2)*BeamHeli*Cos(phi)*SinH(Omega))/
     2    (4.*DDDC*(t_gev**2))

        BeamHeli=-1
        SigmaBHm=-(SigmaBHPolX(1)*BeamHeli*CosH(Omega) + 
     1    SigmaBHPolX(2)*BeamHeli*Cos(phi)*SinH(Omega))/
     2    (4.*DDDC*(t_gev**2))
      endif
           
      if (TargetPolar.eq.2) then
        BeamHeli=+1
        SigmaBHp=-(BeamHeli*SigmaBHPolY*Sin(phi)*SinH(Omega))/
     1    (4.*DDDC*(t_gev**2))

        BeamHeli=-1
        SigmaBHm=-(BeamHeli*SigmaBHPolY*Sin(phi)*SinH(Omega))/
     1    (4.*DDDC*(t_gev**2))
      endif
     
      if (TargetPolar.eq.3) then
        BeamHeli=+1
        SigmaBHp=-(SigmaBHPolZ(1)*BeamHeli*CosH(Omega) + 
     1    SigmaBHPolZ(2)*BeamHeli*Cos(phi)*SinH(Omega))/
     2    (4.*DDDC*(t_gev**2))

        BeamHeli=-1
        SigmaBHm=-(SigmaBHPolZ(1)*BeamHeli*CosH(Omega) + 
     1    SigmaBHPolZ(2)*BeamHeli*Cos(phi)*SinH(Omega))/
     2    (4.*DDDC*(t_gev**2))
      endif
     
c---- SqrAmplVCS(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) ---*
c Computes the expansion of the Virtual Compton Scattering cross section        |
c (divided by phase space) in terms of the "expansion cross sections" computed  |
c in the function MakeExactCrossSections (or MakeLeadingCrossSections).         |
c-------------------------------------------------------------------------------*

      if (TargetPolar.eq.0) then
        BeamHeli=+1
        SigmaVCSp=(SigmaVCSPol0(1) 
     1    - SigmaVCSPol0(4)*Cos(2*phi)*(-1 + CosH(2*Omega)) + 
     2      SigmaVCSPol0(2)*CosH(2*Omega) 
     3    + SigmaVCSPol0(5)*BeamHeli*Sin(phi)*SinH(Omega) + 
     4      SigmaVCSPol0(3)*Cos(phi)*SinH(2*Omega))/(2.*qsq)

        BeamHeli=-1
        SigmaVCSm=(SigmaVCSPol0(1) 
     1    - SigmaVCSPol0(4)*Cos(2*phi)*(-1 + CosH(2*Omega)) + 
     2      SigmaVCSPol0(2)*CosH(2*Omega) 
     3    + SigmaVCSPol0(5)*BeamHeli*Sin(phi)*SinH(Omega) + 
     4      SigmaVCSPol0(3)*Cos(phi)*SinH(2*Omega))/(2.*qsq)
      endif
     
      if (TargetPolar.eq.1) then
        BeamHeli=+1
        SigmaVCSp=(SigmaVCSPolX(1)*BeamHeli*CosH(Omega) - 
     1    SigmaVCSPolX(4)*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     2    SigmaVCSPolX(2)*BeamHeli*Cos(phi)*SinH(Omega) + 
     3    SigmaVCSPolX(3)*Sin(phi)*SinH(2*Omega))/(2.*qsq)

        BeamHeli=-1
        SigmaVCSm=(SigmaVCSPolX(1)*BeamHeli*CosH(Omega) - 
     1    SigmaVCSPolX(4)*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     2    SigmaVCSPolX(2)*BeamHeli*Cos(phi)*SinH(Omega) + 
     3    SigmaVCSPolX(3)*Sin(phi)*SinH(2*Omega))/(2.*qsq)
      endif

      if (TargetPolar.eq.2) then
        BeamHeli=+1
        SigmaVCSp=(SigmaVCSPolY(2) - 
     1    SigmaVCSPolY(3)*Cos(2*phi)*(-1 + CosH(2*Omega)) + 
     2    SigmaVCSPolY(3)*CosH(2*Omega) + 
     3    SigmaVCSPolY(1)*BeamHeli*Sin(phi)*SinH(Omega) + 
     4    SigmaVCSPolY(4)*Cos(phi)*SinH(2*Omega))/(2.*qsq)

        BeamHeli=-1
        SigmaVCSm=(SigmaVCSPolY(2) - 
     1    SigmaVCSPolY(3)*Cos(2*phi)*(-1 + CosH(2*Omega)) + 
     2    SigmaVCSPolY(3)*CosH(2*Omega) + 
     3    SigmaVCSPolY(1)*BeamHeli*Sin(phi)*SinH(Omega) + 
     4    SigmaVCSPolY(4)*Cos(phi)*SinH(2*Omega))/(2.*qsq)
      endif
     
      if (TargetPolar.eq.3) then
        BeamHeli=+1
        SigmaVCSp=(SigmaVCSPolZ(1)*BeamHeli*CosH(Omega) - 
     1    SigmaVCSPolZ(4)*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     2    SigmaVCSPolZ(2)*BeamHeli*Cos(phi)*SinH(Omega) + 
     3    SigmaVCSPolZ(3)*Sin(phi)*SinH(2*Omega))/(2.*qsq)

        BeamHeli=-1
        SigmaVCSm=(SigmaVCSPolZ(1)*BeamHeli*CosH(Omega) - 
     1    SigmaVCSPolZ(4)*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     2    SigmaVCSPolZ(2)*BeamHeli*Cos(phi)*SinH(Omega) + 
     3    SigmaVCSPolZ(3)*Sin(phi)*SinH(2*Omega))/(2.*qsq)
      endif

c-- SqrAmplInterf(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) --*
c Computes the expansion of the Interference cross section (divided by phase    |
c space) in terms of the "expansion cross sections" computed in the function    |
c MakeExactCrossSections (or MakeLeadingCrossSections).                         |
c-------------------------------------------------------------------------------*
  
      if (TargetPolar.eq.0) then
        BeamHeli=+1
        SigmaIp=-((BeamCharge*(SigmaIPol0(1)*CosH(Omega) + 
     1    SigmaIPol0(3)*Cos(2*phi)*(CosH(Omega) - CosH(3*Omega)) + 
     2    SigmaIPol0(2)*CosH(3*Omega) - 
     3    SigmaIPol0(8)*BeamHeli*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     3    SigmaIPol0(7)*BeamHeli*Sin(phi)*SinH(2*Omega) + 
     4    (SigmaIPol0(6)*Cos(3*phi)*(3*SinH(Omega) - SinH(3*Omega)))/3. + 
     5    Cos(phi)*(SigmaIPol0(3)*SinH(Omega) + 
     6    SigmaIPol0(4)*SinH(3*Omega))))/(DDDC*qsq*tt))

        BeamHeli=-1
        SigmaIm=-((BeamCharge*(SigmaIPol0(1)*CosH(Omega) + 
     1    SigmaIPol0(3)*Cos(2*phi)*(CosH(Omega) - CosH(3*Omega)) + 
     2    SigmaIPol0(2)*CosH(3*Omega) - 
     3    SigmaIPol0(8)*BeamHeli*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     3    SigmaIPol0(7)*BeamHeli*Sin(phi)*SinH(2*Omega) + 
     4    (SigmaIPol0(6)*Cos(3*phi)*(3*SinH(Omega) - SinH(3*Omega)))/3. + 
     5    Cos(phi)*(SigmaIPol0(3)*SinH(Omega) + 
     6    SigmaIPol0(4)*SinH(3*Omega))))/(DDDC*qsq*tt))
      endif
     
      if (TargetPolar.eq.1) then
        BeamHeli=+1
        SigmaIp=-((BeamCharge*(-(SigmaIPolX(8)*BeamHeli*Cos(2*phi)*
     1    (-1 + CosH(2*Omega))) + 
     2    BeamHeli*(SigmaIPolX(5) + SigmaIPolX(6)*CosH(2*Omega)) + 
     3    SigmaIPolX(3)*(CosH(Omega) - CosH(3*Omega))*Sin(2*phi) + 
     4    SigmaIPolX(7)*BeamHeli*Cos(phi)*SinH(2*Omega) + 
     5    (SigmaIPolX(4)*Sin(3*phi)*(3*SinH(Omega)-SinH(3*Omega)))/3. + 
     6    Sin(phi)*(SigmaIPolX(1)*SinH(Omega) + 
     7    SigmaIPolX(2)*SinH(3*Omega))))/(DDDC*qsq*tt))

        BeamHeli=-1
        SigmaIm=-((BeamCharge*(-(SigmaIPolX(8)*BeamHeli*Cos(2*phi)*
     1    (-1 + CosH(2*Omega))) + 
     2    BeamHeli*(SigmaIPolX(5) + SigmaIPolX(6)*CosH(2*Omega)) + 
     3    SigmaIPolX(3)*(CosH(Omega) - CosH(3*Omega))*Sin(2*phi) + 
     4    SigmaIPolX(7)*BeamHeli*Cos(phi)*SinH(2*Omega) + 
     5    (SigmaIPolX(4)*Sin(3*phi)*(3*SinH(Omega)-SinH(3*Omega)))/3. + 
     6    Sin(phi)*(SigmaIPolX(1)*SinH(Omega) + 
     7    SigmaIPolX(2)*SinH(3*Omega))))/(DDDC*qsq*tt))
      endif

      if (TargetPolar.eq.2) then
        BeamHeli=+1
        SigmaIp=-((BeamCharge*(SigmaIPolY(1)*CosH(Omega) + 
     1    SigmaIPolY(5)*Cos(2*phi)*(CosH(Omega) - CosH(3*Omega)) + 
     2    SigmaIPolY(2)*CosH(3*Omega) - 
     3    SigmaIPolY(8)*BeamHeli*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     4    SigmaIPolY(7)*BeamHeli*Sin(phi)*SinH(2*Omega) + 
     5    (SigmaIPolY(6)*Cos(3*phi)*(3*SinH(Omega)-SinH(3*Omega)))/3. + 
     6    Cos(phi)*(SigmaIPolY(3)*SinH(Omega) + 
     7    SigmaIPolY(4)*SinH(3*Omega))))/(DDDC*qsq*tt))

        BeamHeli=-1
        SigmaIm=-((BeamCharge*(SigmaIPolY(1)*CosH(Omega) + 
     1    SigmaIPolY(5)*Cos(2*phi)*(CosH(Omega) - CosH(3*Omega)) + 
     2    SigmaIPolY(2)*CosH(3*Omega) - 
     3    SigmaIPolY(8)*BeamHeli*(-1 + CosH(2*Omega))*Sin(2*phi) + 
     4    SigmaIPolY(7)*BeamHeli*Sin(phi)*SinH(2*Omega) + 
     5    (SigmaIPolY(6)*Cos(3*phi)*(3*SinH(Omega)-SinH(3*Omega)))/3. + 
     6    Cos(phi)*(SigmaIPolY(3)*SinH(Omega) + 
     7    SigmaIPolY(4)*SinH(3*Omega))))/(DDDC*qsq*tt))
      endif
     
      if (TargetPolar.eq.3) then
        BeamHeli=+1
        SigmaIp=-((BeamCharge*(-(SigmaIPolZ(8)*BeamHeli*Cos(2*phi)*
     1    (-1 + CosH(2*Omega))) + 
     2    BeamHeli*(SigmaIPolZ(5) + SigmaIPolZ(6)*CosH(2*Omega)) + 
     3    SigmaIPolZ(3)*(CosH(Omega) - CosH(3*Omega))*Sin(2*phi) + 
     4    SigmaIPolZ(7)*BeamHeli*Cos(phi)*SinH(2*Omega) + 
     5    (SigmaIPolZ(4)*Sin(3*phi)*(3*SinH(Omega) - SinH(3*Omega)))/3. + 
     6    Sin(phi)*(SigmaIPolZ(1)*SinH(Omega) + 
     7    SigmaIPolZ(2)*SinH(3*Omega))))/(DDDC*qsq*tt))

        BeamHeli=-1
        SigmaIm=-((BeamCharge*(-(SigmaIPolZ(8)*BeamHeli*Cos(2*phi)*
     1    (-1 + CosH(2*Omega))) + 
     2    BeamHeli*(SigmaIPolZ(5) + SigmaIPolZ(6)*CosH(2*Omega)) + 
     3    SigmaIPolZ(3)*(CosH(Omega) - CosH(3*Omega))*Sin(2*phi) + 
     4    SigmaIPolZ(7)*BeamHeli*Cos(phi)*SinH(2*Omega) + 
     5    (SigmaIPolZ(4)*Sin(3*phi)*(3*SinH(Omega) - SinH(3*Omega)))/3. + 
     6    Sin(phi)*(SigmaIPolZ(1)*SinH(Omega) + 
     7    SigmaIPolZ(2)*SinH(3*Omega))))/(DDDC*qsq*tt))
      endif

c-- CrossSectionBH(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) -*
c Computes the expansion of the Bethe Heitler cross section in terms of the     |
c "expansion cross sections" computed in the function MakeExactCrossSections    | 
c (or MakeLeadingCrossSections).                                                |
c-------------------------------------------------------------------------------*
 
      XSecBHp = SigmaBHp * PhaseSpace
      XSecBHm = SigmaBHm * PhaseSpace
      
c- CrossSectionVCS(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) -*
c Computes the expansion of the Virtual Compton Scattering cross section in     |
c terms of the "expansion cross sections" computed in the function              |
c MakeExactCrossSections (or MakeLeadingCrossSections).                         |
c-------------------------------------------------------------------------------*
 
      XSecVCSp = SigmaVCSp * PhaseSpace
      XSecVCSm = SigmaVCSm * PhaseSpace

c-CrossSectionInterf(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge)*
c Computes the expansion of the Interference cross section in terms of the      |
c "expansion cross sections" computed in the function MakeExactCrossSections    |
c (or MakeLeadingCrossSections).                                                |
c-------------------------------------------------------------------------------*
  
      XSecIp = SigmaIp * PhaseSpace
      XSecIm = SigmaIm * PhaseSpace

      if (thetacm.lt.pi/2.) then
c ad-hoc u-channel model modification for VCS and Interference
c christian weiss recommends the following change for u-channel:
c switch u-slope for t-slope, then divide by 10, since back angle peak 
c is ~10% of forward angle peak (at least for omega electroproduction)
         XSecVCSp=XSecVCSp/10.
         XSecVCSm=XSecVCSm/10.
c Interference should be negligible since BH is so small in u-channel.
c The calculation for interference unreliable due to tt issue.
c For sure, it needs to be much smaller than BH alone
 400     if (abs(XSecIp).gt.XSecBHp/10.) then
            XSecIp=XSecIp/10.
            goto 400
         endif
 401     if (abs(XSecIm).gt.XSecBHm/10.) then
            XSecIm=XSecIm/10.
            goto 401
         endif
      endif

c-------------------------------------------------------------------------------*
c the routine returns d4sigma/dQ2dxdtdphi in nb/GeV4
c i.e. there already is a virtual photon flux factor built in (see Brauel Eqn 2.2)

      SigmaTotPlus  = XSecBHp + XSecVCSp + XSecIp
      SigmaTotMoins = XSecBHm + XSecVCSm + XSecIm

      XSecSum=Pi * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn
      XSecDif=Pi * ( SigmaTotPlus - SigmaTotMoins ) * ConvGeV2nbarn

      if (debug(2)) then
         write(6,*)' XSecSum ',XSecSum
         write(6,*)SigmaTotPlus,SigmaTotMoins
         write(6,*)XSecBHp,XSecVCSp,XSecIp
         write(6,*)XSecBHm,XSecVCSm,XSecIm
      endif

      if (ibh.eq.1) then  ! BH cross section only
         sig_tgendvcs=Pi*(XSecBHp+XSecBHm)* ConvGeV2nbarn
      else
         sig_tgendvcs=XSecSum
      endif

c SIMC requires d6sigma/dOmega_p/dE_p/dOmega_e/dE_e [ub/MeV^2/sr^2]
c multiply by Jacobian computed in logbook, pages 107,112,113
c conversion to MeV units is done in physics_dvcs.f
c NOTE: cos(thetacm) is negative for proton angle > pi/2

      Jac6=(2./abs(cos(thetacm))) * ( (qCM(3)*Ep/Pp)-Ep-qCM(4) ) *
     1    (Elab**2)*E_prime*(1-cos(thetaEscat)) / 
     2    (pi*m_p*(Elab-E_prime)) * (1.-Elab/(Elab-E_prime)) 

      sig_tgvdvcs = sig_tgendvcs*Jac6

      if (debug(2)) then
         write(6,*)' tgvdvcs: xsec ',PhaseSpace,Jac6
         write(6,*)' BH:  ',SigmaBHp,SigmaBHm
         write(6,*)' VCS: ',SigmaVCSp,SigmaVCSm
         write(6,*)' IP:  ',SigmaIp,SigmaIm
      endif

      return

 990  continue
      sig_tgvdvcs = 0.
      return

 1000 write(6,*)' tgvdvcs: Error reading CFFoutput_LO.dat'
      stop

      end
