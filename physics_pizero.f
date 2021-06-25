      real*8 function peep_pizero(vertex,main)

c     gh - 18.02.04
      
      implicit none
      include 'simulate.inc'

C     The following two record lines are from SIMC physics_kaon.f

      type(event_main):: main
      type(event):: vertex

*     NOTE: when we refer to the center of mass system, it always refers to the
*     photon-NUCLEON center of mass

      real*8 sigma_pi0
      real*8 m_p, m_psq

      real*8 E0                 ! Electron beam energy
      real*8 E_prime            ! scattered electron energy
      real*8 nu                 ! virtual photon energy
      real*8 qvec               ! virtual photon lab energy
      real*8 qsq                ! 4 momentum transfer of scattered electron
      real*8 epsilon            ! virtual photon polarisation

      real*8 Ep,Pp

      real*8 mass               ! pi0 mass
      real*8 invm,Wsq           ! invariant mass of hadronic system
      real*8 tt,uu              ! Mandelstam variables
      real*8 tprime,t_min,t_max,uprime,u_min,uu2,umin2
      real*8 e_photCM,e_pi0CM,e_pCM

      real*8 gamma_T            ! flux of transversely polarized virt. photons
      real*8 Breit_wigner
      real*8 Jttheta            ! Jacobian t-Mx - Ep,theta_p in LAB
      real*8 Jttheta_fx         ! old Jacobian for fixed Mx
      real*8 tcos               ! cos of theta_LAB between Pp and q
      real*8 vertex_pi0_x,vertex_pi0_y,vertex_pi0_z
      real*8 P_pi0,E_pi0,ucos,uu3

! Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz      !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz     !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz      !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian,jac_old

      real*8 sig_pi0gmh
      real*8 sig1,sig2,sig3

      logical first
      save first
      data first /.TRUE./

      if (debug(2)) write(6,*)' peep_pizero: enter '

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
      vertex_pi0_x = qvec*vertex%uq%x - Pp*vertex%up%x
      vertex_pi0_y = qvec*vertex%uq%y - Pp*vertex%up%y
      vertex_pi0_z = qvec*vertex%uq%z - Pp*vertex%up%z
      P_pi0=sqrt(vertex_pi0_x**2 +vertex_pi0_y**2 +vertex_pi0_z**2)
      E_pi0=sqrt(P_pi0**2+mass**2)
      
c      ucos =  vertex_pi0_x/P_pi0*vertex%uq%x
c     1       +vertex_pi0_y/P_pi0*vertex%uq%y
c     2     +vertex_pi0_z/P_pi0*vertex%uq%z

      e_photCM = (Wsq - qsq - m_psq)/invm/2.
      e_pi0CM  = (Wsq + mass**2 - m_psq)/invm/2.
      e_pCM    = (Wsq + m_psq - mass**2)/invm/2.

c gh - 18.02.04
c although t is single valued (negative), this is not true for u
c right at the kinematic endpoint, u turns positive, so umin is positive
c even though for most events u is negative.  This is surprising, but is
c confirmed by checking against s+t+u relation.

      if ((e_pi0CM**2-mass**2)*(e_photCM**2+qsq).ge.0.) then
         t_min = -qsq + mass**2 -2.*(e_pi0CM*e_photCM
     *        -sqrt( (e_pi0CM**2-mass**2)*(e_photCM**2+qsq) ))
         t_max = -qsq + mass**2 -2.*(e_pi0CM*e_photCM
     *        +sqrt((e_pi0CM**2-mass**2)*(e_photCM**2+qsq)))
         u_min = -qsq + m_psq -2.*(e_pCM*e_photCM
     *        -sqrt( (e_pCM**2-m_psq)*(e_photCM**2+qsq) ))
      else
         write(6,*)' physics_pizero: no valid t,u min '
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
      sig1 = sig_pi0gmh(qsq/1.e6,tt/1.e6,t_min/1.e6,uu/1.e6,u_min/1.e6,
     1                  invm/1.e3,epsilon,thetacm,phicm,first)

      ntup%sigcm=sig1

      ntup%phicmi=phicm
      ntup%wcmi=invm/1.e3
      ntup%tprimei=tprime/1.e6
      ntup%epsiloni=epsilon
      ntup%ui=-uu/1.e6

c      sigma_pizero=ntup%sigcm

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

      peep_pizero = Jttheta_fx*gamma_T*ntup%sigcm

      if (debug(2)) write(6,*)' peep_pizero: end ',peep_pizero
      return 
      end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function sig_pi0gmh(Q2_g,t_gev,tmin,u_gev,umin,
     1       W_gev,eps,thetacm,phicm,first)

c This routine calculates p(e,e'pi0)p cross sections, based on
c a parameterization of the Defurne Hall A data and a parameterization
c of a calculation by Goloskokov & Kroll
c Subroutine calculates dsigma/dt/dphi_cm, which is returned as sig_pi0gmh
c [ub/MeV^2/rad].
c gh - 18.02.01

      implicit none
      include 'constants.inc'

      real*8 t_gev,tmin,tprime,u_gev,umin,uprime
      real*8 mass,eps,Q2_g,W_gev,wfactor,m_p
      real*8 Q2hi,Q2lo,delQ2,Q2tmp
      real*8 Whi,Wlo,Wfac_hi,Wfac_lo
      real*8 sigThi,sigTlo,sigThiW,sigTloW
      real*8 sigLhi,sigLlo,sigLhiW,sigLloW
      real*8 sigTThi,sigTTlo,sigTThiW,sigTTloW
      real*8 sigLThi,sigLTlo,sigLThiW,sigLTloW
      real*8 thetacm,phicm,sig_phi0gmh
      real*8 sig,sigT,sigL,sigLT,sigTT    !components of dsigma/dt
      
c use fit parameters for xB=0.36 data+model (nb/GeV2)
      real*8 Q2tab(4) / 1.75,  3.00,  4.00,  5.50/
      real*8 Wtab(4)  / 2.00,  2.46,  2.83,  3.26/
c sigT is parameterized in the form p1+p2/(-t)
      real*8 p1T(4)   /1577.,  168.,  67.4,  24.7/
      real*8 p2T(4)   /-125., -11.1,  -4.6, -1.66/
c sigL is simply constant vs t
      real*8 p1L(4)   /   0.,  2.77,  1.16,  0.43/
c sigTT  p1+p2/(-t)
      real*8 p1TT(4)  /-674., -141., -57.6, -21.3/
      real*8 p2TT(4)  / 102.,  25.4,  10.4,  3.82/
c sigLT  p1+p2/(-t)
      real*8 p1LT(4)  / 156., 12.57,  5.17,  1.97/
      real*8 p2LT(4)  /-60.0, -2.08, -0.85, -0.32/
      
      integer pcount,ndat,Q2count,Q2c
      logical first

      ndat=4
      
      m_p = Mp/1.e3
      tprime = abs(t_gev-tmin)
      uprime = abs(u_gev-umin)
c      write(6,*)' uprime ',u_gev,umin,uprime

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
            
      if (thetacm.gt.pi/2.) then ! t-channel
         sigThi = p1T(Q2count+1)+p2T(Q2count+1)/(tprime+abs(tmin))
         sigTlo = p1T(Q2count)  +p2T(Q2count)/(tprime+abs(tmin))
         sigLhi = p1L(Q2count+1)
         sigLlo = p1L(Q2count)
         sigTThi = p1TT(Q2count+1)+p2TT(Q2count+1)/(tprime+abs(tmin))
         sigTTlo = p1TT(Q2count)  +p2TT(Q2count)/(tprime+abs(tmin))
         sigLThi = p1LT(Q2count+1)+p2LT(Q2count+1)/(tprime+abs(tmin))
         sigLTlo = p1LT(Q2count)  +p2LT(Q2count)/(tprime+abs(tmin))
         if (first) write(6,*)' pi0gmh: t-channel '
         first=.FALSE.
         
      else                      ! u-channel
c     christian weiss recommends the following change for u-channel:
c     switch u-slope for t-slope, then divide by 10, since back angle peak 
c     is ~10% of forward angle peak (at least for omega electroproduction)
         
         sigThi = (p1T(Q2count+1)+p2T(Q2count+1)/(uprime+abs(tmin)))/10.
         sigTlo = (p1T(Q2count)  +p2T(Q2count)/(uprime+abs(tmin)))/10.  
         sigLhi = p1L(Q2count+1)/10.
         sigLlo = p1L(Q2count)/10.
         sigTThi = (p1TT(Q2count+1)+p2TT(Q2count+1)/(uprime+abs(tmin)))/10.
         sigTTlo = (p1TT(Q2count)  +p2TT(Q2count)/(uprime+abs(tmin)))/10.  
         sigLThi = (p1LT(Q2count+1)+p2LT(Q2count+1)/(uprime+abs(tmin)))/10.
         sigLTlo = (p1LT(Q2count)  +p2LT(Q2count)/(uprime+abs(tmin)))/10.  
         if (first) write(6,*)' pi0gmh: u-channel '
         first=.FALSE.
      endif
      
c sighi,lo are at different W.  scale both to the W needed for the event

      Wfac_hi= ((Whi**2-m_p**2)**2) / ((W_gev**2-m_p**2)**2)
      Wfac_lo= ((Wlo**2-m_p**2)**2) / ((W_gev**2-m_p**2)**2)
      
      sigThiW = sigThi*Wfac_hi
      if (sigThiW.lt.0.) then
         write(6,*)' pizero: sigThiW<0 ',sigThiW,uprime,abs(tmin)
         sigThiW=0.
      endif
      sigTloW = sigTlo*Wfac_lo
      if (sigTloW.lt.0.) then
         write(6,*)' pizero: sigTloW<0 ',sigTloW,uprime,abs(tmin)
         sigTloW=0.
      endif
      sigLhiW = sigLhi*Wfac_hi
      if (sigLhiW.lt.0.) then
         write(6,*)' pizero: sigLhiW<0 ',sigLhiW,uprime,abs(tmin)
         sigLhiW=0.
      endif
      sigLloW = sigLlo*Wfac_lo
      if (sigLloW.lt.0.) then
         write(6,*)' pizero: sigLloW<0 ',sigLloW,uprime,abs(tmin)
         sigLloW=0.
      endif
      sigTThiW = sigTThi*Wfac_hi
c SIGN(A,B) returns the value of A with the sign of B
      if (abs(sigTThiW).gt.sigThiW) sigTThiW=sign(sigThiW,sigTThiW)
      sigTTloW = sigTTlo*Wfac_lo
      if (abs(sigTTloW).gt.sigTloW) sigTTloW=sign(sigTloW,sigTTloW)
      sigLThiW = sigLThi*Wfac_hi
      if (abs(sigLThiW).gt.sigThiW) sigLThiW=sign(sigThiW,sigLThiW)
      sigLTloW = sigLTlo*Wfac_lo
      if (abs(sigLTloW).gt.sigTloW) sigLTloW=sign(sigTloW,sigLTloW)

c interpolate to get cross section at Q2 needed for the event
      
      if( Q2count.le.(ndat-1) .and. Q2tmp.ge.Q2tab(Q2count) .and.
     1     Q2tmp.lt.Q2tab(Q2count+1) ) then
         
         sigT  = ( sigTloW*(Q2hi-Q2tmp)+ sigThiW*(Q2tmp-Q2lo))/delQ2
         sigL  = ( sigLloW*(Q2hi-Q2tmp)+ sigLhiW*(Q2tmp-Q2lo))/delQ2
         sigTT = (sigTTloW*(Q2hi-Q2tmp)+sigTThiW*(Q2tmp-Q2lo))/delQ2
         sigLT = (sigLTloW*(Q2hi-Q2tmp)+sigLThiW*(Q2tmp-Q2lo))/delQ2

c        write(6,*)' T1 ',Q2tmp,sigt,sigTloW,sigThiW
      
      elseif (Q2tmp.ge.Q2tab(ndat) ) then
         
         sigT  =  sigThiW+ (sigThiW-sigTloW) /delQ2
         sigL  =  sigLhiW+ (sigLhiW-sigLloW) /delQ2
         sigTT = sigTThiW+(sigTThiW-sigTTloW)/delQ2
         sigLT = sigLThiW+(sigLThiW-sigLTloW)/delQ2
         
c         write(6,*)' T2 ',Q2tmp,sigt,sigTloW,sigThiW

      elseif (Q2tmp.le.Q2tab(1) ) then
         
         sigT  =  sigTloW- (sigThiW-sigTloW) /delQ2
         sigL  =  sigLloW- (sigLhiW-sigLloW) /delQ2
         sigTT = sigTTloW-(sigTThiW-sigTTloW)/delQ2
         sigLT = sigLTloW-(sigLThiW-sigLTloW)/delQ2
         
c         write(6,*)' T3 ',Q2tmp,sigt,sigTloW,sigThiW

      else

         write(6,*)' Q2tmp error ',Q2tmp,Q2count

      endif
      
      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt

      sig_pi0gmh = sig/2./pi*1.e-09 !dsig/dtdphicm in microbarns/MeV^2/rad
      
      return
      end
