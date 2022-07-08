CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                      CCC
CCC   F1F221 version 1.0 prerelease 1  -  April 4, 2022                  CCC
CCC   Collection of subroutines to calculate inclusive cross sections    CCC
CCC   for range of nuclei.  For A > 2 the parameterization is based on   CCC
CCC   by M. E. Christy, T. Gautam, and A Bodek to 12C, 27Al, 56Fe and    CCC
CCC   64Cu.  However, the fit scales relatively well with 'A' and        CCC
CCC   should be good for all nuclei with 10 < A < 80.                    CCC
CCC   Also included is the proton cross section fit and a preliminary    CCC
CCC   deuteron/neutron fit by M. E. Christy, N. Kalantarians, J. Either  CCC
CCC   and W. Melnitchouk (to be published) based on both inclusive       CCC
CCC   deuteron and tagged deuteron data on n/d from BONuS.               CCC
CCC   Range of validity is W^2 < 32,  Q^2 < 32.0                         CCC
CCC   New data included in deuteron fit, including photoproduction at    CCC
CCC   Q^2 = 0.                                                           CCC
CCC                                                                      CCC
CCC                                                                      CCC
CCC   A > 2 nuclei fit are only 12C for this pre release version.        CCC
CCC   Do Not use for other nuclei!  These will be added prior to         CCC
CCC   final reselease.                                                   CCC  
CCC                                                                      CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE F1F2IN21(Z, A, QSQ, WSQ, F1, F2)                       
!--------------------------------------------------------------------
! Fit to inelastic cross sections for A(e,e')X
! valid for all W^<20 GeV2 and all Q2<30 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c Outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 02/20/2021 E. Christy       
! Made to be consistent with calling F1F2IN09 from P. Bosted
! with much code borrowed. 
!--------------------------------------------------------------------

      implicit none

      real*8 Z,A,QSQ,WSQ
      real*8 avgn,F1,F2,FL,F1p,F1n,F2p,F2n,FLp,FLn
      real*8 sigm,sigl,sigt,eps
      integer IA,IZ,wfn,opt
      logical off,doqe/.false./

      off = .true.
      IA = int(A)
      IZ = int(Z)
      avgN = A-z
      if(IA.LT.2) then
        call sf(WSQ,QSQ,F1p,FLp,F2p,F1n,FLn,F2n)
        if(IZ.LT.1) then        !!! Neutron
          F1 = F1n
          F2 = F2n
        else                    !!! Proton
          F1 = F1p
          F2 = F2p
        endif
      elseif(IA.EQ.2.AND.IZ.EQ.1) then      !!! Deuteron
        eps = 0.5               !!! pass  0 < eps < 1
        wfn = 2                 !!! CD-Bonn, default
        doqe = .false.          !!! don't include QE  !!!

        f1 = 0.0
        f2 = 0.0
        fL = 0.0
        call rescsd(WSQ,QSQ,eps,doqe,F1,F2,FL,wfn,sigm)

c        write(6,*) "rescsd inel: ", wsq,qsq,eps,wfn,f1,f2
        
c        call smearsub(WSQ,QSQ,wfn,off,F2,F1,FL)

c        write(6,*) "smearsub: ", wsq,qsq,eps,wfn,f1,f2        
        
c        write(6,*) "in f1f2in21: ", wsq,qsq,eps,wfn,f1,f2
        
        
      elseif(IA.GT.2) then      !!! A>2 nuclei
         opt = 3                !!! inelastic only
        call SFCROSS(WSQ,QSQ,A,Z,opt,sigt,sigl,F1,F2,FL)
      endif

      return
      end


      
CCC-----------------
       SUBROUTINE F1F2QE21(Z, A, QSQ, WSQ, F1qe, F2qe)                       
!--------------------------------------------------------------------
! Fit to QE cross sections for A(e,e')X
! valid for all W^<30 GeV2 and all Q2<30 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c Outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 02/20/2021 E. Christy
!
! Note:  This is a calling routine in order to be consistent with
!        calling of F1F2QE09 from P. Bosted.  For the deuteron the 
!        QE distribution is calculated by smearing with a realistic
!        wavefunction in the weak-binding approximation.  For A>2
!        a superscaling model is used for the fit.     
!--------------------------------------------------------------------
      IMPLICIT none
      real*8 wsq,qsq,A,Z,sigt,sigl,F1qe,F2qe,FLqe
      integer opt/5/                !!! QE + MEC  !!!
      integer wfn/2/
      logical first/.false./

      if(A.EQ.2.0) then
        call SQESUB(wsq,qsq,wfn,f2qe,f1qe,fLqe,first) 
      elseif(A.GT.2.0) then 
        call sfcross(wsq,qsq,A,Z,opt,sigt,sigl,F1qe,F2qe,FLqe)
      endif
        
c      write(6,*) wsq,qsq,a,z,f1qe,f2qe,fLqe
      
      end


      
CCC-----------------
C=======================================================================

      Subroutine FORMFACTS(q2,gmp,gep,gmn,gen)

CCC   Returns proton and neutron form factors based on new proton fit         CCC
CCC   by M.E. Christy including GMp12 data.  Neutron starts with              CCC
CCC   Kelly fit, but includes correction factors extracted from fit to        CCC      
CCC   inclusive deuteron data.   Version from July 22, 2020                   CCC                      
!--------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 q2,tau,gd,gmp,gep,gmn,gen,mcor,ecor
      REAL*8 mu_p/ 2.792782/ ,mu_n/ -1.913148 /
      REAL*8 mp/ 0.9382727 /, mp2
      
      mp2 = mp*mp
      tau   = q2 / 4.0 / mp2 

CCC   2021 Christy fit to GMp and GEp including GMp12 data  CCC
        
       GMP = mu_p*(1.+0.099481*tau)/
     &  (1.0+11.089*tau+19.374*tau*tau+5.7798*tau**3)
       GEP = (1.+0.24482*tau**2)/
     &  (1.0+11.715*tau+11.964*tau*tau+27.407*tau**3)       
       GD = (1./(1 + q2/0.71))**2

CCC   2021 fit to deuteron inclusive data   CCC       
      
c        GMn = mu_n*(1.0+0.12417E-04*tau)/
c     &   (1.000+11.404*tau+17.014*tau*tau+31.219*tau**3)

c         GEn = (1.5972*tau / (1.0+0.19655*tau)) * GD

        GMn = mu_n           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &      / (1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)
        
        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        GEn = GEn*((q2+1189.4)/1189.4)**219.73 
        GMn = GMn/((q2+0.35590 )/0.35590 )**0.93020E-01 
        
             
       return
       end
CCC  Version 120521  -  Author:  M.E. Christy, christy@jlab.org             CCC
CCC  This fit version includes data from a number of JLab Hall experiments  CCC
CCC  as well as DIS data from SLAC (L. Whitlow) and photoproduction data    CCC
CCC  from DAPHNE and older data sets.                                       CCC
CCC  Subroutine to get Transverse and Longitudinal eP cross sections        CCC  
CCC  from fits cross sections over a range of epsilon.  The subroutine      CCC
CCC  resmod.f is required.  Units are in ub/Sr/Gev.                         CCC
CCC  
CCC   Region of applicability has been extended to cover the full JLab      CCC
CCC   11 GeV kinematic range of Q^2 < 30 GeV^2 and W^2 < 20                 CCC


      SUBROUTINE rescsp(W2,Q2,sigT,sigL)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval/
     & 0.12291E+01,0.15173E+01,0.15044E+01,0.17100E+01,0.16801E+01,
     & 0.14312E+01,0.12616E+00,0.23000E+00,0.92594E-01,0.90606E-01,
     & 0.75000E-01,0.35067E+00,0.75729E+01,0.56091E+01,0.94606E+01,
     & 0.20156E+01,0.66190E+01,0.41732E+00,0.23980E-01,0.53136E+01,
     & 0.63752E+00,0.11484E+02,0.69949E-01,0.26191E+01,0.53603E-01,
     & 0.65000E+02,0.15351E+00,0.20624E+01,0.23408E+01,0.16100E+02,
     & 0.62414E+02,0.17201E+01,0.23261E+00,0.65000E+02,0.23292E+01,
     & 0.14980E+01,0.23000E+00,0.63385E+00,0.19093E-01,0.61061E-01,
     & 0.29146E-02,0.54388E+00,0.77997E+00,0.28783E+00,0.10605E+01,
     & 0.69793E+00,0.20009E+01,0.57000E+00,0.41632E+01,0.38427E+00,
     & 0.10000E+01,0.99842E+00,0.98719E+00,0.10168E+01,0.98945E+00,
     & 0.99594E+00,0.98799E+00,0.10271E+01,0.10650E+01,0.97920E+00,
     & 0.10152E+01,0.99622E+00,0.81011E+01,0.10070E-02,0.14857E+01,
     & 0.33445E+01,0.31641E-09,0.69755E+02,0.55228E+01,0.14438E+00,
     & 0.60474E+01,0.65395E-07,0.14129E+01,0.58609E+00,0.36220E+01,
     & 0.92699E+00,0.14418E+01,0.86403E-02,0.10001E-03,0.75106E+00,
     & 0.76077E+00,0.42272E+00,0.55511E-11,0.52486E+00,0.58153E+00,
     & 0.15798E+01,0.50105E+00,0.89149E+02,0.72789E+00,0.24813E-01,
     & -.61906E+00,0.10000E+01,0.00000E+00,0.00000E+00,0.68158E+03,
     & 0.12429E+01,0.00000E+00,0.00000E+00,0.00000E+00,0.10000E-05 /
           
      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 

        if(i.LE.12) xvalL(i) = xval1(i)
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo

       
      xb = q2/(w2+q2-mp2)

      call resmodp(1,w2,q2,xval1,sigT)
      call resmodp(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
c      write(6,*) w2,q2,F1,FL,R

      end


       




      SUBROUTINE RESCSD(w2,q2,eps,doqe,f1d,f2d,fLd,wfn,sigm)
CCCC  Calcultes deuteron cross sections from smeared p+n.             CCCC
CCCC  Requires SQESUB be called first to read in smearing functions.  CCCC
CCCC  This should be one at the lowest level to keep from reading     CCCC
CCCC  multiple times.

      IMPLICIT none

      real*8 w2,q2,eps,t1,x,gamma2,q2min,q2t
      real*8 sigtp,sigtn,siglp,sigln,sigt,sigl,sigm,m  
      real*8 m2,pi2,alpha,f1d,f2d,fLd,f1dqe,f2dqe,fLdqe
      real*8 off_mKP_fit,delta,xfr,xfl
      integer i,j,k,ntssot,wfn,drn
      logical doqe,off
c      INCLUDE 'parm.cmn'
      external off_mKP_fit

      f1dqe = f1d
      f2dqe = f2d
      fLdqe = fLd      

      
      off = .false. !! off-shell corrections
c      off = .true.
      q2min = 4.0E-5
      drn = 5
      xfr = 0.95
      xfl = 1.0E-3
      alpha = 1./137.036
      m = (0.938272+0.939565)/2.0d0  !!! average p,n
      m2 = m*m
      pi2 = 3.14158*3.14159

      q2t = q2
c      if(q2.LT.q2min) q2t = q2min  !!! Hack to fix photoproduction
      
      x = q2t/(w2-m2+q2t)
c     gamma2 = 1.+4.*m2*x*x/q2t
      gamma2 = 1.0+q2t*4.0*m2/(w2+q2t-m2)**2.0
      t1 = gamma2*f2dqe-2.*x*f1dqe

      if(f2dqe.LT.0.0.OR.f1dqe.LT.0.0.OR.fLdqe.LT.0.0)
     &    write(34,*) w2,q2,f2dqe,f1dqe,fLdqe
        
      
      if(q2.LT.0.05) then
        call SMEARSUBLOWQ(w2,q2t,wfn,off,f2d,f1d,fLd)
      else
        call SMEARSUB(w2,q2t,wfn,off,f2d,f1d,fLd)
      endif

c      write(6,*) w2,q2t,f2d,f1d,fLd
      
      if(doqe) then 
       f1d = f1d+f1dqe
       f2d = f2d+f2dqe
       fLd = fLd+fLdqe
      endif

          
      sigt = 0.3894e3*f1d*pi2*alpha*8.0/abs(w2-m2)
      sigl = 0.3894e3*fLd*pi2*alpha*8.0/abs(w2-m2)/2.*abs(w2-m2+q2t)/q2t
      if(q2t.EQ.0) sigl = 0.0
      sigm = sigt+eps*sigl

c      write(6,*) w2,q2t,sigt,sigl,sigm
      
      return
      end




      SUBROUTINE RESCSN(w2,q2,sigtn,sigln)
CCCC  Returns neutron transverse and longitudinal cross sections   CCCC
CCCC  November 16, 2021 version                                    CCCC      
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i
      data xvaln / 
     & 0.12291E+01,0.15173E+01,0.15044E+01,0.17100E+01,0.16801E+01,
     & 0.14312E+01,0.12616E+00,0.23000E+00,0.92594E-01,0.90606E-01,
     & 0.75000E-01,0.35067E+00,0.69500E+01,0.86633E+01,0.11557E+02,
     & 0.22138E+01,0.44886E+01,0.20500E+03,0.84433E+03,0.31167E+01,
     & 0.96301E+00,0.14956E+00,0.20761E-07,0.10440E+01,0.40143E-03,
     & 0.90028E+02,0.75248E-01,0.20532E+00,0.12444E-01,0.34469E+03,
     & 0.19948E+00,0.26925E+01,0.48635E+01,0.86000E+02,0.67813E+04,
     & 0.44281E+02,0.29548E+00,0.65421E+00,0.23787E-09,0.51967E-01,
     & 0.39926E-08,0.29960E+00,0.97516E+00,0.46934E-01,0.14246E+03,
     & 0.55801E+00,0.19349E+01,0.27400E+00,0.38891E+00,0.40000E-02,
     & 0.10108E+01,0.97020E+00,0.98248E+00,0.97768E+00,0.10425E+01,
     & 0.10198E+01,0.97822E+00,0.98239E+00,0.10103E+01,0.10076E+01,
     & 0.10044E+01,0.99687E+00,0.16696E+01,0.10721E-06,0.54114E+00,
     & 0.11923E+04,0.55938E+02,0.95000E+03,0.39840E+02,0.22026E+03,
     & 0.30498E+01,0.24459E+00,0.95574E+00,0.35596E+00,0.21228E-05,
     & 0.96696E+01,0.27563E+01,0.93027E-01,0.33559E+02,0.31207E-01,
     & 0.29020E+02,0.86417E+00,0.36471E-08,0.99167E+00,0.68124E+00,
     & 0.10000E-01,0.90227E-01,0.40115E+01,0.29915E+01,0.45929E-01,
     & -.16758E+01,0.78493E+01,0.78184E+01,0.42074E+01,0.41179E-05,
     & 0.80597E+00,0.00000E+00,0.00000E+00,0.10045E+01,0.62364E+00 /
      
          
      do i=1,50
        xval1(i) = xvaln(i)
        xvalL(i) = xvaln(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo

      call resmodn(1,w2,q2,xval1,sigtn)
      call resmodn(2,w2,q2,xvalL,sigLn)

      return
      end
     


CCC-----------------

      SUBROUTINE RESMODP(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from July 1, 2021  -  Author:  M.E. Christy                    CCC
CCC  This routine returns proton photo-absorbtion cross sections            CCC
CCC  for either transverse or longitudinal photons in units of ub/Sr/Gev.   CCC
CCC                                                                         CCC
CCC  Fit form is empirical.  Interpret physics from it at your own risk.    CCC
CCC  replaced 2-pi threshold with eta                                       CCC
      
      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,sig,xval(50),mass(7),width(7)
      REAL*8 height(7),rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,pi2,alpha
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,A0
      REAL*8 sig_res,xpr(2),t1,t2
      INTEGER i,j,num,sf
      
      mp = 0.9382727
      mpi = 0.134977
      mpi2 = mpi*mpi
      meta = 0.547862
      mp2 = mp*mp
      alpha = 1./137.036
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + meta)
      
      q20 = 0.05
      q20= xval(50)

       
CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.00       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.60      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.60       !!!  S11(1650)
      br(6,1) = 0.65     !!!  P11(1440) roper 
      br(7,1) = 0.60      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.40      !!!  S11(1535) 
      br(3,3) = 0.08      !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.20      !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo

CCCC   Meson angular momentum   CCCC

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.160           !!!
c         x0(i) = 0.14
      enddo

c      x0(1) = 0.155

c      if(sf.EQ.2) x0(1) = 0.08   !!!  different Delta mass for sigL
      if(sf.EQ.2) x0(1) = 0.07   !!!  different Delta mass for sigL      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.15)**2. !!!  Dipole parameterization  !!!
      
c     mon = 1./(1.+q2/1.5)**1.
      mon = 1./(1.+q2/1.5)**1.
      

      
      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.00+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+meta)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LE.(mp+mpi)) xpr(1) = 1.0      
      if(w.LE.(mp+meta)) xpr(2) = 1.0
     

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      mass(7) = xval(47)
      intwidth(7) = xval(48)
      width(7) = intwidth(7) 

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))*
     &          mon**rescoef(i,4)
         
        else

           height(i) = (rescoef(i,1)+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

          
        endif
 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

         height(7) = (xval(16)+xval(20)*q2)*exp(-1.0*xval(24)*q2)       
        
      else
        height(7) = xval(49)*mon**xval(45)
      endif
      
      height(7) = height(7)*height(7)


CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7

       sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w
      if(sf.EQ.2) sig_res = sig_res*q2


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1.and.xpr(1).LT.1.0) then

         A0 = xval(37)/(1.0+q2/xval(42))**xval(43)
         
c        t1 = xval(38)*log(2.1+q2/xval(39))                 !!! exponent of (1-xpr      
c        t2 = xval(40)*log(2.1+q2/xval(41))        !!! exponent of xpr

c         t1 = xval(38)*log(1.5+q2/xval(39))
         t1 = xval(38)*log(1.06+q2)+xval(39)/log(1.06+q2)

c         t2 = xval(40) + xval(41)*log(1.1+q2)
         t2 = xval(40)*(1.0+q2/xval(41))**xval(44)

c         t2 = xval(40)*log(1.1+q2)+xval(41)/log(1.0+q2)
         
c         t2 = xval(40)*log(1.1+q2)+xval(41)/log(1.1+q2)
         
        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = 389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
        endif

c        write(6,*) q2,sf,t1,t2
         
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         A0 = xval(37)/(1.0+q2/xval(39))**2.0
         t1 = xval(38)/(1.0+q2/(xval(40)))+xval(32)*log(q2+xval(36))
         
c         t2 = xval(41)/(1.00+q2/xval(42))**xval(43)       !!! exponent of xpr

         t2 = xval(41)
         
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + 389.4*A0*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2
        endif

      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 



      SUBROUTINE RESMODN(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from July 15, 2021  -  Author:  M.E. Christy                    CCC
CCC  This routine returns proton photo-absorbtion cross sections            CCC
CCC  for either transverse or longitudinal photons in units of ub/Sr/Gev.   CCC
CCC                                                                         CCC
CCC  Fit form is empirical.  Interpret physics from it at your own risk.    CCC
CCC  replaced 2-pi threshold with eta                                       CCC
      
      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,sig,xval(50),mass(7),width(7)
      REAL*8 height(7),rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,pi2,alpha
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,A0
      REAL*8 sig_res,xpr(2),t1,t2
      INTEGER i,j,num,sf
      
      mp = 0.939565
      mpi = 0.134977
      mpi2 = mpi*mpi
      meta = 0.547862
      mp2 = mp*mp
      alpha = 1./137.036
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + meta)
      
      q20 = 0.05
      q20= xval(50)

       
CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.00       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.60      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.60       !!!  S11(1650)
      br(6,1) = 0.65     !!!  P11(1440) roper 
      br(7,1) = 0.60      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.40      !!!  S11(1535) 
      br(3,3) = 0.08      !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.20      !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo

CCCC   Meson angular momentum   CCCC

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.16   !!! 
      enddo

      if(sf.EQ.2) x0(1) = 0.07   !!!  different Delta mass for sigL
  
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.15)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.5)**1.

      
      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.00+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+meta)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LE.(mp+mpi)) xpr(1) = 1.0      
      if(w.LE.(mp+meta)) xpr(2) = 1.0
     

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      mass(7) = xval(47)
      intwidth(7) = xval(48)
      width(7) = intwidth(7) 

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))*
     &          mon**rescoef(i,4)


        else

           height(i) = (rescoef(i,1)+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

          
        endif
 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = (xval(44)+xval(45)*q2)*exp(-1.0*xval(46)*q2)
        
      else
        height(7) = xval(49)*mon
      endif
      
      height(7) = height(7)*height(7)


CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7

       sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w
      if(sf.EQ.2) sig_res = sig_res*q2


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1.and.xpr(1).LT.1.0) then
c     A0 = xval(37)*(1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude

        A0 = xval(37)/(1.0+q2/xval(42))**xval(43)
         
        t1 = xval(38)*log(1.05+q2)+xval(39)/(1.05+q2)                 !!! exponent of (1-xpr)
c     t2 = xval(40)*log(2.0+q2/xval(41))+xval(45)        !!! exponent of xpr
        t2 = xval(40)*(1.0+q2/xval(41))**xval(44)
 


        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = 389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif

c        write(6,*) q2,sf,t1,t2
         
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         A0 = xval(37)/(1.0+q2/xval(39))**2.0
         t1 = xval(38)/(1.0+q2/(xval(40)))+xval(32)*log(q2+xval(36))                    !!! exponent of (1-xpr)
         t2 = xval(41)/(1.00+q2/xval(42))**xval(43)       !!! exponent of xpr

c         write(6,*) q2,a0,t1,t2
        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + 389.4*A0*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2

        endif

      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 



      SUBROUTINE SMEARSUB(w2,q2,wfn,off,f2s,f1s,fLs)
CCCC   Fermi smears structure functions.  Requires subroutne GETFY  CCCC
CCCC   to be initialized first.                                     CCCC

      IMPLICIT none

      real*8 w2,q2,mp,mp2,x,z,wwz,alpha,alpha2,gamma,gamma2
      real*8 y,fy11,fy12,fy2,inc,f2s,f1s,fLs,f2pi,f2ni,xvaln(100)
      real*8 sigtp,sigtn,siglp,sigln,f1pi,f1ni,fLpi,fLni
      real*8 f1i(1000),f2i(1000),vfact1(1000),vfact2(1000)
      real*8 vfact1off(1000),vfact2off(1000),fy11off,fy12off
      real*8 fy2off,f1soff,f2soff,fLsoff
      real*8 ymin,ymax,pi,pi2,epsd
      real*8 c0,c1,c2,delta(1000)
      integer i,j,nbins,wfn
c      logical first/.true./
      logical firsty/.false./
      logical off
      real*8 dsimps
      external dsimps

C     *****Off-shell parameters*****
      if(wfn.eq.1) then
         c0=-3.6735d0
         c1=0.057717d0
         c2=0.36419d0
      elseif(wfn.eq.2) then
         c0=-7.2061d0
         c1=0.062022d0
         c2=0.38018d0
      elseif(wfn.eq.3) then
         c0=1.7204d0
         c1=0.050923d0
         c2=0.34360d0
      elseif(wfn.eq.4) then
         c0=-1.7690d0
         c1=0.058321d0
         c2=0.36976d0
      else
         c0=0.d0
         c1=0.d0
         c2=0.d0
         write(6,*)  "Warning:  wavefunction undefined"
      endif
C     ******************************

      mp = (0.938272+0.939565)/2.0d0  !!! average of p, n
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha 
      pi = 3.14159
      pi2 = pi*pi 
      epsd = -2.2E-03
      x = q2/(w2-mp2+q2)

      f2s = 0.0
      f1s = 0.0
      fLs = 0.0
      if(x.LT.0.0D0) return

      gamma2 = (1.+4.*mp2*x*x/q2)   
      gamma = sqrt(gamma2)
      nbins = 60

      if(w2.LT.1.6) nbins = 76

c      nbins = 200   !!! adjust as needed

      if(w2.GE.2.6.AND.Q2.GE.2.0) nbins = 48
      
CCC///   Fill array for Simpson's rule calculation   ///CCC
       

c      if(x.LT.1.0) then
c        ymin = max(x*(1.0D0-2.0D0*epsd*mp/(w2-mp2)),x)
c      else
c        ymin = x
c      endif 
CCC   Note problem if Q2 = 0  CCC       
      ymin = max(x*(1.0D0-2.0D0*epsd*mp/q2),x)
c      ymin = min(ymin,2.0d0)
      ymax = min(1.0d0+gamma2/2.0d0+epsd/mp,2.0d0)

c      write(6,*) "W2, ymin, ymax:  ",w2,ymin,ymax 

      if(ymin.GE.ymax) return

c      ymax = 2.0d0      

c      write(6,*) off

      if(ymin.EQ.0) write(6,*) "Error ymin = 0" 
      inc = (ymax-ymin)/float(nbins)
      y = ymin-inc

      !$OMP PARALLEL DO PRIVATE(i,y,z,wwz,F1pi,FLpi,f1ni,FLni,f2pi,f2ni,fy11,fy12,fy2)
      do i=1,nbins+1      !!!!    First calculate the cross smearing function for each y bin    !!!!
        y = y + inc
        z = x/y

        wwz = mp2 + q2*(1.d0/z - 1.d0)   !!! W^2 at x=z  !!!

        call rescsp(wwz,q2,sigtp,sigLp) 
        call rescsn(wwz,q2,sigtn,sigLn)
  
        f1pi = 2.*z*sigtp/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1p
        fLpi = 2*q2/abs(wwz-mp2+q2)*sigLp/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)                                     
        f1ni = 2.*z*sigtn/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1n
        fLni =  2*q2/abs(wwz-mp2+q2)*sigLn/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)  
        f2pi = (f1pi+fLpi)/(1.+4.*mp2*z*z/q2)
        f2ni = (f1ni+fLni)/(1.+4.*mp2*z*z/q2)
        f2i(i) = f2pi+f2ni  !!!!  Proton + Neutron  !!!
        f1i(i) = f1pi+f1ni  !!!!  Actually 2xF1     !!!
        delta(i) = c0*(z-c1)*(z-c2)*(1.0-c1+z)

        call getfy(gamma,y,wfn,fy11,fy12,fy2,firsty)
        if(off) then
          call getfyoff(gamma,y,wfn,fy11off,fy12off,fy2off,firsty)
        endif
        vfact1(i) = fy11*f1i(i)+fy12*f2i(i)
        vfact1off(i) = (fy11off*f1i(i)+fy12off*f2i(i))*delta(i)
        vfact2(i) = fy2*f2i(i)
        vfact2off(i) = fy2off*f2i(i)*delta(i)
        
        if(vfact1(i).LT.0.0) vfact1(i) = 0.0
        if(vfact2(i).LT.0.0) vfact2(i) = 0.0
        if(off) then
          vfact1(i) = vfact1(i)+vfact1off(i)
          vfact2(i) = vfact2(i)+vfact2off(i) 
        endif
      enddo
      !$OMP END PARALLEL DO

      vfact1(1) = 0.0d0
      vfact2(1) = 0.0d0

CCC///    Integrate over y    ///CCC

      f1s = dsimps(vfact1,ymin,ymax,nbins)/2./x  !!! F1d
      f2s = dsimps(vfact2,ymin,ymax,nbins)       !!! F2d
      fLs = gamma2*f2s-2.*x*f1s                  !!! FLd

      
    
c      write(6,*) "in smearsub 2:  ",x,ymin,q2,gamma,f1s,f2s,fLs

 2000 format(7f12.4)

      end







      SUBROUTINE SMEARSUBLOWQ(w2,q2,wfn,off,f2s,f1s,fLs)
CCCC   Fermi smears structure functions.  Requires subroutne GETFY  CCCC
CCCC   to be initialized first.                                     CCCC

      IMPLICIT none

      real*8 w2,q2,mp,mp2,x,z,wwz,alpha,alpha2,gamma,gamma2
      real*8 t,y,fy11,fy12,fy2,inc,f2s,f1s,fLs,f2pi,f2ni,xvaln(100)
      real*8 sigtp,sigtn,siglp,sigln,f1pi,f1ni,fLpi,fLni
      real*8 f1i(1000),f2i(1000),vfact1(1000),vfact2(1000)
      real*8 vfact1off(1000),vfact2off(1000),fy11off,fy12off
      real*8 fy2off,f1soff,f2soff,fLsoff
      real*8 ymin,ymax,pi,pi2,epsd
      real*8 c0,c1,c2,delta(1000)
      integer i,j,nbins,wfn
c      logical first/.true./
      logical firsty/.false./
      logical off
      real*8 dsimps
      external dsimps

C     *****Off-shell parameters*****
      if(wfn.eq.1) then
         c0=-3.6735d0
         c1=0.057717d0
         c2=0.36419d0
      elseif(wfn.eq.2) then
         c0=-7.2061d0
         c1=0.062022d0
         c2=0.38018d0
      elseif(wfn.eq.3) then
         c0=1.7204d0
         c1=0.050923d0
         c2=0.34360d0
      elseif(wfn.eq.4) then
         c0=-1.7690d0
         c1=0.058321d0
         c2=0.36976d0
      else
         c0=0.d0
         c1=0.d0
         c2=0.d0
         write(6,*)  "Warning:  wavefunction undefined"
      endif
C     ******************************

      mp = (0.938272+0.939565)/2.0d0  !!! average of p, n
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha 
      pi = 3.14159
      pi2 = pi*pi 
      epsd = -2.2E-03
      x = q2/(w2-mp2+q2)

      f2s = 0.0
      f1s = 0.0
      fLs = 0.0
      if(x.LT.0.0D0) return

      gamma2 = (1.+4.*mp2*x*x/q2)
      gamma2 = 1.0+q2*4.0*mp2/(w2+q2-mp2)**2.0
c      if(q2.EQ.0.0) gamma2 = 0.0
      
      gamma = sqrt(gamma2)
      
      nbins = 60

      if(w2.LT.1.6) nbins = 76


c      nbins = 200   !!! adjust as needed

      if(w2.GE.2.6.AND.Q2.GE.2.0) nbins = 48
      
CCC///   Fill array for Simpson's rule calculation   ///CCC
       

c      if(x.LT.1.0) then
c        ymin = max(x*(1.0D0-2.0D0*epsd*mp/(w2-mp2)),x)
c      else
c        ymin = x
c      endif 
CCC   Note problem if Q2 = 0  CCC       
      ymin = max(x*(1.0D0-2.0D0*epsd*mp/q2),x)

      if(q2.EQ.0.0) ymin = 0.001   !!!  hack
      
c      ymin = min(ymin,2.0d0)
      ymax = min(1.0d0+gamma2/2.0d0+epsd/mp,2.0d0)

c      if(q2.EQ.0.0) write(6,*) w2,ymin,ymax 

      if(ymin.GE.ymax) return

c      ymax = 2.0d0      

c      write(6,*) off

      if(ymin.EQ.0) write(6,*) "Error ymin = 0" 
      inc = (ymax-ymin)/float(nbins)
      y = ymin-inc

      !$OMP PARALLEL DO PRIVATE(i,y,z,wwz,F1pi,FLpi,f1ni,FLni,f2pi,f2ni,fy11,fy12,fy2)
      do i=1,nbins+1      !!!!    First calculate the cross smearing function for each y bin    !!!!
        y = y + inc
        z = x/y
        t = 1.0+4.0*mp2/y/(q2+w2-mp2)
        t = 1.0+4.0*mp2*z*z/q2

        wwz = mp2 + q2*(1.d0/z - 1.d0) !!! W^2 at x=z  !!!
        
        if(q2.EQ.0.0) wwz = mp2+y*(w2-mp2)

        call rescsp(wwz,q2,sigtp,sigLp) 
        call rescsn(wwz,q2,sigtn,sigLn)
  
        f1pi = sigtp/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1p
        fLpi = 2*q2/abs(wwz-mp2+q2)*sigLp/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)                                     
        f1ni = sigtn/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1n
        fLni =  2*q2/abs(wwz-mp2+q2)*sigLn/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)  

         f2pi = (2.0*z*f1pi+fLpi)/t
         f2ni = (2.0*z*f1ni+fLni)/t  

c         if(q2.EQ.0.0) write(6,*) w2,q2,z,wwz,t,f1pi,f1ni 
         
        f2i(i) = f2pi+f2ni  !!!!  Proton + Neutron  !!!
        f1i(i) = f1pi+f1ni   
        delta(i) = c0*(z-c1)*(z-c2)*(1.0-c1+z)

        call getfy(gamma,y,wfn,fy11,fy12,fy2,firsty)
        if(off) then
          call getfyoff(gamma,y,wfn,fy11off,fy12off,fy2off,firsty)
       endif

       if(q2.EQ.0.0) f2i(i) = 0.0

       
c        if(q2.LE.0.0) write(6,*) q2,w2,gamma,y,wwz

       
        vfact1(i) = fy11*f1i(i)+fy12*f2i(i)
        
       
        vfact1off(i) = (fy11off*f1i(i)+fy12off*f2i(i))*delta(i)
        vfact2(i) = fy2*f2i(i)
        vfact2off(i) = fy2off*f2i(i)*delta(i)

c        if(q2.EQ.0.0)  write(6,*) w2,q2,i,vfact1(i),fy12,f2i(i)
        
        if(vfact1(i).LT.0.0) vfact1(i) = 0.0
        if(vfact2(i).LT.0.0) vfact2(i) = 0.0
        if(off) then
          vfact1(i) = vfact1(i)+vfact1off(i)
          vfact2(i) = vfact2(i)+vfact2off(i) 
        endif
      enddo
      !$OMP END PARALLEL DO

      vfact1(1) = 0.0d0
      vfact2(1) = 0.0d0

CCC///    Integrate over y    ///CCC


      f1s = dsimps(vfact1,ymin,ymax,nbins)  
c        f1s = dsimps(vfact1,ymin,ymax,nbins)/2./x  !!! F1d   - problem at x=0
      f2s = dsimps(vfact2,ymin,ymax,nbins)       !!! F2d
      fLs = gamma2*f2s-2.*x*f1s                  !!! FLd
      
      
c      if(q2.EQ.0.0) write(6,*) w2,q2,gamma2,nbins,f1s,f2s,fLs

c      write(6,*) w2,q2,2.0*x/gamma2*f1s,f2s,fLs
      
c       write(6,*) "in smearsub 2:  ",x,q2,2.0*x*f1s,f2s,fLs

 2000 format(7f12.4)

      end







      SUBROUTINE SQESUB(w2,q2,wfn,f2s,f1s,fLs,first)
      IMPLICIT none

      real*8 w2,wwz,q2,x,y,z,Z1,A,t1
      real*8 nu,nuel,q4,q,gamma,gamma2,alpha,alpha2
      real*8 epel,w2el,kappa,mp,mp2,md,w2min
      real*8 f1,f2,fy11,fy12,fy2,f2s,f1s,fLs,f1sos,f2sos,fLsos
      real*8 fy11off,fy12off,fy2off
      real*8 ymin,pi,pi2,GM,GE,GMp,GMn,GEp,GEn,tau,mcor,ecor
      real*8 GM2,GE2,GD,mu_p,mu_n,a1,a2,a3,b1,b2,b3,b4,b5
      integer i,j,k,bin,off,wfn
      logical thend,first,firsty
      real*8 f1d_of,sf_of
      external f1d_of,sf_of

      thend = .false.
      firsty = .false.
      if(first) firsty = .true.

      Z1 = 1.
      A = 2.
      mp = (0.938272+0.939565)/2.0d0    !!! average p,n
      mp2 = mp*mp
      md = 1.8756
      alpha = 1/137.03599
      alpha2 = alpha*alpha 
      pi = 3.14159
      pi2 = pi*pi
      mu_p = 2.793D0
      mu_n = -1.913148
    

CCC///   Read in smearing function array  ///CCC
CCC

      nu = (w2+q2-mp2)/2./mp   !!! Fermi smeared !!!

      nuel = q2/2./mp  !!! In nucleon rest frame !!!
      q = sqrt(q2)
      q4 = q2*q2
      x = q2/(w2+q2-mp2)
      tau = q2/(4*mp2)  !!!  for elastic at this Q^2  !!!
c      ww = mp2 + 2.d0*mp*nu - q2  !!! Fermi smeared  !!!
      gamma2 = (1.+4.*mp2*x*x/q2)
      if(gamma2.LE.6.and.firsty) then   
       gamma = sqrt(gamma2)
       call getfy(gamma,x,wfn,fy11,fy12,fy2,firsty)
       firsty = .true.
       call getfyoff(gamma,x,wfn,fy11off,fy12off,fy2off,firsty)
      endif

      if(nu.LT.(q2/2.0/md)) then
        f1s = 0.0D0
        f2s = 0.0D0
        fLs = 0.0D0
        return
      endif
   

      call formfacts(q2,GMP,GEP,GMN,GEN)
      

      GM2 = GMp*GMp + GMn*GMn  !!! Sum p+n  !!!
      GE2 = GEp*GEp + GEn*GEn


      off = 2   !!! CC2
      call offshellqe(w2,q2,wfn,off,GM2,GE2,f1s,f2s)
      fLs = gamma2*f2s-2.*x*f1s

c        call f1f2qe09(Z1,A,q2,w2,f1s,f2s)
c        write(6,*) "sqesub:  ",w2,q2,f1s,f2s,fLs

c        if(fLs.LT.0.0) write(6,2000) w2,q2,f1s,f2s,fLs

      if(f1s.LT.0.0) f1s = 0.0
      if(f2s.LT.0.0) f2s = 0.0
      if(fLs.LT.0.0) fLs = 0.0
 
c       if(w2.LT.1.0) write(6,2000) w2,q2,f1s,f2s,fLs

 2000 format(8f9.4)
c 2000 format(5f9.4,2e11.3)

      end









      SUBROUTINE offshellqe(w2,q2,wfn,off,GM2,GE2,f1os,f2os)
CCC    Original code by M.E. Christy with significant borrowing of code from    CCC
CCC    J. Ethier, W. Melnitchouk, et.al.                                        CCC
CCC    Addition of CC1 offshell by J. Ethier (July 1, 2014)                     CCC
CCC    See reference (put references here)                                      CCC
CCC    Oct. 9, fixed bug to keep from returning NAN if a < b                    CCC
c  Nucleon light-cone momentum distribution function in deuteron in
C  weak binding approximation (WBA), as a function of the fraction (y)
C  of deuteron's momentum carried by nucleon, defined as
C		y = p.q/p_D.q
C  (sometimes labeled as y_D) so that in Bjorken limit y -> y0 in [0,1].
C
C  Kulagin, Petti, NPA 765, 126 (2006), Eq. (43)
C  Kahn, WM, Kulagin, PRC 79, 035205 (2009)
C
C  Relativistic kinematics implemented (WM): May 2010
C  Relativistic convolution implemented (cf. AQV): Sep 2010
C  
C  All momenta in MeV/c
C
C  Added missing normalization for wavefunctions (MEC): April 2015
C ***********************************************************************

      IMPLICIT none

      real*8 x,w2,qsqr,q2,GM2,GE2,FF,pFF,off1,off2,off3,corr
      real*8 q,gamma,gamma2,mp,mp2,md,md2,hc,tau,xth,y,ep,p2
      real*8 pv,pv2,pz,pz2,pt2,xd,w2d,f1,f2,ft,fttilda,f2tilda,flux
      real*8 cc,pvmin,pvmax,a,b,c,rt,prod,ftv,fttildav(10001),f2v
      real*8 f2tildav(10001),f1os,f2os,dsimps,dgauss,pvinc,wfnorm(10)
      real*8 U,W,VS,VT
      external dsimps,dgauss
      integer i,j,nbins,wfn,off
      logical thend,first,firsty
      data wfnorm / 1.0,1.0,1.05,1.02,1.0,1.0,1.0,1.0,1.0,1.0 /

c      write(6,2001) x,q2,gm2,ge2

c      off = 2

      f1os = 0.0D0
      f2os = 0.0D0

      nbins = 1000

c      do i=1,10
c        wfnorm(i) = 1.0D0
c      enddo

c      write(6,*) w2,q2,wfn

      hc = 0.197327D0      !!! convert from GeV -> fm    !!!
      mp = (0.938272+0.939565)/2.0D0    !!! average p,n  !!!
      mp2 = mp*mp
      md = 2.0D0*mp-0.002224575D0
    
      x = q2/(w2+q2-mp2)
      xd = x*mp/md
      w2d = md*md+q2*(1.0D0/xd-1.0D0)

      qsqr = q2*(1.D0 + ( q2 / (4.D0*mp2*x*x)))

      tau = q2/(4.0D0*mp2)  !!!  for elastic at this Q^2  !!!
      gamma2 = (1.+4.0D0*mp2*x*x/q2)   
      gamma = sqrt(gamma2)
      a = sqrt(1.d0 + w2d/qsqr)
      b = w2d/(2.D0*mp*sqrt(qsqr))
      prod = b/(a*a - 1.D0)
      c = 1.D0+((a*a-1.0D0)*(1.0-a*a/b/b))
 
      if(c.GE.0.0D0) then
        rt = sqrt(c)
      else
        return
      endif

      ff = (GE2+tau*GM2)/(1.D0+tau)
      pff = ((sqrt(GE2)-sqrt(GM2))/(1.D0+tau))**2      !!! Pauli FF squared

      pvmin = mp*prod*(1.d0-rt)*sign(1.,(a-b))
      pvmax = mp*prod*(1.D0+rt)
      pvinc = (pvmax-pvmin)/float(nbins)

c      write(6,*) "offshellqe:  ",w2,q2,x,rt

      pv = pvmin-pvinc
CCC   Loop over initial nucleon 3-momentum, pv  CCC
      !$OMP PARALLEL DO PRIVATE(i,ep,p2,xth,y,pz,pz2,pt2,flux,off1,off2,corr)
      
      do i=1,nbins+1              !!! integrate over nucleon 3-momentum 
        fttildav(i) = 0.0D0
        f2tildav(i) = 0.0D0
        pv = pv+pvinc
        pv2 = pv*pv              !!! 3-momentum squared of nucleon 
        ep = sqrt(mp2+pv2)       !!! total energy of nucleon
        p2 = (md-ep)**2-pv2      !!! vituality of nucleon
        xth =  q2/(q2-p2+mp2)    !!! x for given virutuality  
        y = x/xth
        pz = (y*mp-md+ep)/gamma  !!! nucleon momentum along q-vector
        pz2 = pz*pz
        pt2 = pv2-pz2            !!! nucleon transverse momentum
        flux = 1.0D0+gamma*pz/mp

c         if(w2.LT.1.0) write(6,2000) w2,q2,f1s,f2s,fLs,f1sos,f2sos,fLsos

        off1 = (mp2-p2)/q2
        off2 = (mp2-p2)/4./mp2
        corr = 1.D0+xth*xth*(4.D0*p2+6.D0*pt2)/q2

        U = 0.D0        ! S-wave
        W = 0.D0        ! D-wave
        VS = 0.D0       ! Singlet P-wave
        VT = 0.D0       ! Triplet P-wave

        IF (wfn.EQ.2) CALL CDBONN(pv/hc,U,W)
        IF (wfn.EQ.3) CALL WJC1(pv/hc,U,W,VS,VT)
              ! wfns normalized to ~105% (5% in V' term)
        IF (wfn.EQ.4) CALL WJC2(pv/hc,U,W,VS,VT)
              ! wfns normalized to ~102%

C...Output wavefunctions in fm^3/2 => MeV^-3/2
           U = U / hc**1.5D0
           W = W / hc**1.5D0
           VS = VS / hc**1.5D0
           VT = VT / hc**1.5D0
           CC = (U**2 + W**2 + VS**2 + VT**2)
        if(wfn.EQ.1) then
          call av18(pv/hc,cc)    
          cc = cc/(hc*hc*hc) 
        endif 

c        write(6,*) "offshellqe:  ", x,pv,p2,xth,pvmin,pvmax

        if(off.EQ.1) then    !!!  CC1 
          ftv = x*(x/y)*(gm2*(1.D0-off1))
          f2v = y*(ff-off2*pff)
          fttilda = ftv + (2.0D0*xth*xth*pt2*f2v/q2)
          fttildav(i) = fttilda*flux*cc*mp/(4.*gamma)
          fttildav(i) = fttildav(i)*2.D0*pv  
          f2tilda = (mp*x*(ff-off2*pff))/(4.d0*gamma**3)
          f2tildav(i) = f2tilda*flux*corr*cc/xth
          f2tildav(i) = f2tildav(i)*2.D0*pv
        else if(off.EQ.2) then    !!!  CC2
          ftv = x*(x/y)*(gm2-off1*(ff-off2*pff))
          f2v = y*ff
          fttilda = ftv + (2.0D0*xth*xth*pt2*f2v/q2)
          fttildav(i) = fttilda*flux*cc*mp/(4.*gamma)
          fttildav(i) = fttildav(i)*2.D0*pv
          f2tilda = (mp*x*ff)/(4.d0*gamma**3) 
          f2tildav(i) = f2tilda*flux*corr*cc/xth
          f2tildav(i) = f2tildav(i)*2.D0*pv
        endif

      enddo
      !$OMP END PARALLEL DO

      ft = dsimps(fttildav,pvmin,pvmax,nbins)
      f2os = dsimps(f2tildav,pvmin,pvmax,nbins)

      f1os = ft/2./x

 
      f1os = f1os/wfnorm(wfn)
      f2os = f2os/wfnorm(wfn)


c      if(f1os.LT.0.0.OR.f1os.LT.0.0) 
c     &      write(6,2001) w2,q2,f1os,f2os,gamma2


      if(f1os.LT.0.0D0) f1os = 0.0D0
      if(f2os.LT.0.0D0) f2os = 0.0D0

 2000 format(8f9.4)
 2001 format(5f9.4)

      return 
      end










      SUBROUTINE GETFY(gamma,ypass,wfn,fy11,fy12,fy2,firsty)
      IMPLICIT none
      real*8 gamma,y,ypass,gammav(50),fy,yv(50,1001),yv2(50,1001)
      real*8 yv12(50,1001),fyv11(50,1001),fyv2(50,1001),fyv12(50,1001)
      real*8 t1,fylow,fyhi,fy11,fy12,fy2
      integer i,j,k,bin,bin2,jlow,jhi,wfn
      logical thend,firsty
      COMMON/FYCOM/ yv,yv2,yv12,fyv11,fyv2,fyv12,gammav

      thend = .false.


c      if(firsty) write(6,*) firsty
      
      y = ypass
      if(ypass.GE.1.9999) y = 1.999  !!! For numerical stability

CCC///   Read in smearing function array  ///CCC

c      firsty = .true.
      i = 0
      if(firsty) then
        if(wfn.EQ.1) then
          open(unit=34,file='f1f2tables/11.WBARELav18',status='old')
          open(unit=35,file='f1f2tables/f2.WBARELav18',status='old')
          open(unit=36,file='f1f2tables/f12.WBARELav18',status='old')
        elseif(wfn.EQ.2) then
          open(unit=34,file='f1f2tables/f11-onshell-cdbonn.dat',status='old')        
          open(unit=35,file='f1f2tables/f22-onshell-cdbonn.dat',status='old')
          open(unit=36,file='f1f2tables/f12-onshell-cdbonn.dat',status='old')
        elseif(wfn.EQ.3) then
          open(unit=34,file='f1f2tables/f11.WBARELwjc1',status='old')        
          open(unit=35,file='f1f2tables/f2.WBARELwjc1',status='old')
          open(unit=36,file='f1f2tables/f12.WBARELwjc1',status='old') 
        elseif(wfn.EQ.4) then
          open(unit=34,file='f1f2tables/f11.WBARELwjc2',status='old')        
          open(unit=35,file='f1f2tables/f2.WBARELwjc2',status='old')
          open(unit=36,file='f1f2tables/f12.WBARELwjc2',status='old')
        else
          write(6,*)  "Warning:  wavefunction undefined"
        endif

     
        dowhile(.not.thend)   !!!  First read in smearing function to temp vector!!!
          i = i+1
          j = int(float(i-1)/999.)+1 
          k = i-(j-1)*999
          read(34,*,end=999) gammav(j),yv(j,k),fyv11(j,k)
          read(35,*,end=999) gammav(j),yv2(j,k),fyv2(j,k)
          read(36,*,end=999) gammav(j),yv12(j,k),fyv12(j,k)

        enddo
999     thend = .true. 
      endif
      firsty = .false.
      close(34)
      close(35)
      close(36)

      bin = int(y/2.*1000)+1          !!! find y-bin below elastic at each x,Q^2
      jlow = int((gamma-1.0)/0.2)+1   !!! Assumes 0.2 increments in gamma
      jhi = jlow+1
      bin2 = bin+1

      fylow = (fyv11(jlow,bin2)*(y-yv(jlow,bin))+fyv11(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv11(jhi,bin2)*(y-yv(jhi,bin))+fyv11(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy11 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2

c      write(6,*) y,fyv11(jlow,bin)

      fylow = (fyv12(jlow,bin2)*(y-yv(jlow,bin))+fyv12(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv12(jhi,bin2)*(y-yv(jhi,bin))+fyv12(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy12 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2


      fylow = (fyv2(jlow,bin2)*(y-yv(jlow,bin))+fyv2(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv2(jhi,bin2)*(y-yv(jhi,bin))+fyv2(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy2 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2


      if(fy11.LT.0) fy11 = 0.0d0
      if(fy2.LT.0) fy2 = 0.0d0
      if(fy12.LT.0) fy12 = 0.0d0


c      if(y.LT.0.015) then   !!! fix for numerical issues in fy at small x
c        fy11 = 0.0
c        fy12 = 0.0
c        fy2 = 0.0
c      endif

     
 2000 format(5f9.4,2e11.3)

      return
      end




      SUBROUTINE GETFYOFF(gamma,ypass,wfn,fy11off,fy12off,fy2off,
     &                     firsty)
      IMPLICIT none
      real*8 gamma,y,ypass,gammavoff(50),fy,yvoff(50,1001)
      real*8 yv2off(50,1001),yv12off(50,1001),fyv11off(50,1001)
      real*8 fyv2off(50,1001),fyv12off(50,1001)
      real*8 t1,fylow,fyhi,fy11off,fy12off,fy2off
      integer i,j,k,bin,bin2,jlow,jhi,wfn
      logical thend,firsty
      COMMON/FYCOMOFF/ yvoff,yv2off,yv12off,fyv11off,fyv2off,fyv12off,
     &                 gammavoff

      thend = .false.

      y = ypass
      if(ypass.GE.1.9999) y = 1.999  !!! For numerical stability

CCC///   Read in smearing function array  ///CCC

c      firsty = .true.
      i = 0
      if(firsty) then
        if(wfn.EQ.1) then
          open(unit=37,file='f1f2tables/f11.WBARELav18OFF',status='old')        
          open(unit=38,file='f1f2tables/f2.WBARELav18OFF',status='old')
          open(unit=39,file='f1f2tables/f12.WBARELav18OFF',status='old')
        elseif(wfn.EQ.2) then  !!!!  Not available yet - fix later
          open(unit=37,file='f1f2tables/f11-offshell-cdbonn.dat',status='old')        
          open(unit=38,file='f1f2tables/f22-offshell-cdbonn.dat',status='old')
          open(unit=39,file='f1f2tables/f12-offshell-cdbonn.dat',status='old')
        elseif(wfn.EQ.3) then
          open(unit=37,file='f1f2tables/f11.WBARELwjc1OFF',status='old')        
          open(unit=38,file='f1f2tables/f2.WBARELwjc1OFF',status='old')
          open(unit=39,file='f1f2tables/f12.WBARELwjc1OFF',status='old') 
        elseif(wfn.EQ.4) then
          open(unit=37,file='f1f2tables/f11.WBARELwjc2OFF',status='old')        
          open(unit=38,file='f1f2tables/f2.WBARELwjc2OFF',status='old')
          open(unit=39,file='f1f2tables/f12.WBARELwjc2OFF',status='old')
        else
          write(6,*)  "Warning:  wavefunction undefined"
        endif

     
        dowhile(.not.thend)   !!!  First read in smearing function to temp vector!!!
          i = i+1
          j = int(float(i-1)/999.)+1 
          k = i-(j-1)*999

          read(37,*,end=999) gammavoff(j),yvoff(j,k),fyv11off(j,k)
          read(38,*,end=999) gammavoff(j),yv2off(j,k),fyv2off(j,k)
          read(39,*,end=999) gammavoff(j),yv12off(j,k),fyv12off(j,k)

        enddo
999     thend = .true. 
      endif
      firsty = .false.
      close(37)
      close(38)
      close(39)

      bin = int(y/2.*1000)+1          !!! find y-bin below elastic at each x,Q^2
      jlow = int((gamma-1.0)/0.2)+1   !!! Assumes 0.2 increments in gamma
      jhi = jlow+1
      bin2 = bin+1

      fylow = (fyv11off(jlow,bin2)*(y-yvoff(jlow,bin))+
     &         fyv11off(jlow,bin)*(yvoff(jlow,bin2)-y))/0.002
      fyhi =  (fyv11off(jhi,bin2)*(y-yvoff(jhi,bin))+fyv11off(jhi,bin)
     &        *(yvoff(jhi,bin2)-y))/0.002

      fy11off = (fylow*(gammavoff(jhi)-gamma)+
     &           fyhi*(gamma-gammavoff(jlow)))/0.2

c      write(6,*) y,fyv11off(jlow,bin)

      fylow = (fyv12off(jlow,bin2)*(y-yvoff(jlow,bin))+
     &         fyv12off(jlow,bin)*(yvoff(jlow,bin2)-y))/0.002
      fyhi =  (fyv12off(jhi,bin2)*(y-yvoff(jhi,bin))+fyv12off(jhi,bin)
     &        *(yvoff(jhi,bin2)-y))/0.002

      fy12off = (fylow*(gammavoff(jhi)-gamma)+
     &          fyhi*(gamma-gammavoff(jlow)))/0.2


      fylow = (fyv2off(jlow,bin2)*(y-yvoff(jlow,bin))+fyv2off(jlow,bin)
     &        *(yvoff(jlow,bin2)-y))/0.002
      fyhi =  (fyv2off(jhi,bin2)*(y-yvoff(jhi,bin))+fyv2off(jhi,bin)
     &        *(yvoff(jhi,bin2)-y))/0.002

      fy2off = (fylow*(gammavoff(jhi)-gamma)+
     &         fyhi*(gamma-gammavoff(jlow)))/0.2


c      if(fy11off.LT.0) fy11off = 0.0d0
c      if(fy2off.LT.0) fy2off = 0.0d0
c      if(fy12off.LT.0) fy12off = 0.0d0


c      if(y.LT.0.015) then   !!! fix for numerical issues in fy at small x
c        fy11 = 0.0
c        fy12 = 0.0
c        fy2 = 0.0
c      endif

     
 2000 format(5f9.4,2e11.3)

      return
      end




      SUBROUTINE QENUC21OFF(Z, A, Q2, W2, xvalc,F1, F2)

C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS  
c for A>2. Uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c based on the earlier code F1F2QE09 by P. Bosted.  The superscaling distribution
c shape is determined from the fit to 12C data      
c
c input: Z, A  (real*8) Z and A of nucleus (should be 2.0D0 for deuteron)
c        Q2 (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        W2 (real*8) is invariant mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus


      IMPLICIT NONE     

      REAL*8 Z, A, avgN, F1, F2, W2, Q2, R
      REAL*8 mp/0.938272/
      REAL*8 PAULI_SUP1, PAULI_SUP2,x2,pb2,pb2L
      REAL*8 GEP, GEN, GMP, GMN
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      REAL*8 nu, qv, TAU, FY, FYL, FYp, FYn, pauli, nup,num
      REAL*8 kappa, lam, lamp, lampn,taup, xi, ximax, psi, psip, psipn
      REAL*8 psimax, nuL,nut,kf, es, esmin, GM2bar, GE2bar, Delta, GL,GT
      REAL*8 F1ff,F2ff,GMoff,GEoff,moff,xvalc(45),B,xm,m,BL
      integer IA,paulitype

      psimax = 5.0
      moff =  1.0*mp
      paulitype = 1    !!! 1 = Superscaling, 2 = Tsai from Fermi Gas  !!!
      
! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or q2 is negative
      Nu = (W2 - mp**2 + Q2) / 2. / mp
      if(nu .le. 0.0 .or. Q2 .le. 0.) return
      TAU   = Q2 / 4.0 / mp**2                                        
      qv = sqrt(nu**2 + Q2)
      
!     Call New FormFactors  !
      
      call formfacts(Q2,gmp,gep,gmn,gen)

!  Get energy shift and fermi momementum from fit  !      
      
      if(IA.GT.2) then
         Esmin = 0.020           !!! Minimum Es given by removal energy
         if(IA.EQ.12) Esmin = 0.016
         kf = xvalc(35)
         Es = Esmin + xvalc(36)  !!! Maximum Es
         if(qv.LT.1.5) Es = Es-xvalc(36)*(1.0-qv/1.5)**0.1  !!!  test  !!!        
      endif 
      nup = nu-Es

! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

!     structure functions with off shell factors
      
      kappa = qv / 2. / mp
      lam = nu / 2. / mp
      
      lamp = nup / 2. / mp
      lampn = -lamp
      taup = kappa**2 - lamp**2
      xi = sqrt(1. + (kf/moff)**2) -1.
! Very close to treshold, could have a problem
c      if(1.+lamp.le.0.) return
c      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - taup ) / sqrt(xi) /
     >     sqrt((1.+lam )* taup + kappa * sqrt(taup * (1. + taup)))   !!! OK

      psip = (lamp - taup) / sqrt(xi) / 
     >     sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))    !!! OK

      psipn = (lampn - taup) / sqrt(xi) / 
     >     sqrt((1.+lampn)*taup + kappa * sqrt(taup * (1. + taup)))   !!! OK

      
      
      nuL = (Q2/qv/qv)**2
      
c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 


      GM2bar = (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = (Z * GEP**2 + avgN * GEN**2)

      F1ff =  (tau*sqrt(GM2bar)+sqrt(GE2bar))/(1.0+tau) 
      F2ff = (sqrt(GM2bar)-sqrt(GE2bar))/(1.0+tau)
      GEoff = F1ff-tau*moff/mp*F2ff         !!!  Off-shell FF   !!!
      GMoff = F1ff+moff/mp*F2ff             !!!  Off-shell FF   !!!
     
c      write(6,*) "off :",q2,sqrt(GM2bar),GMoff,sqrt(GE2bar),GEoff
      
      Delta = tau/kappa/kappa*xi*(1.-psi**2)*
     & (kappa*sqrt(1.0+1.0/tau)+xi/3.0*(1.-psi**2))

c      write(6,*) w2,q2,Delta
      
      GL = kappa**2 / tau*
     >  (GEoff**2. +(GEoff**2.+tau*GMoff**2.)*Delta/(1.0+tau))
 
      GT = (2.*tau*GMoff**2.+(GEoff**2.+tau*GMoff**2.)*Delta/(1.+tau))
      
      num = 2. *kappa *(1. + xi * (1. + psi**2) / 2.)

      GL = GL/num
      GT = GT/num

c       GL = 2.0*xi*kf*(Moff/kf)**3.0/qv*GL
c       GT = 2.0*xi*kf*(Moff/kf)**3.0/qv*GT
  

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))

CCMEC - Fitted Superscaling distribution      CCCC

      FYp = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >  (1. + exp(xvalc(10) * psip))*(1.0-abs(psip)/psimax)**2.0*
     >  (1.0+abs(psip)/psimax)**2.0

      if(psip.GT.psimax) Fyp = 0.0
      FYp = max(0.0,FYp) 

      FYn = xvalc(7)/ (1. + xvalc(8)**2 * (psipn + xvalc(9))**2) / 
     > (1. + exp(xvalc(10) * psipn))*(1.0-abs(psipn)/psimax)**2.0*
     > (1.0+abs(psipn)/psimax)**2.0

      
      FYn = max(0.0,FYn)      
      
      if(paulitype.EQ.1) then        !!! Superscaling
        FY = max(0.0,FYp-FYn)

c        if(qv.LT.0.1) write(6,*) qv,w2,FYp,FYn,FY
        
        x2 = qv/kf
        
        pb2L = 1.0

        xm = xvalc(40)  !!! 
        m = 2.0
        BL = xm**m/(1.0-1.0/xm/xvalc(41))
        
c        if(x2.LT.xm) then
c           pb2L = xvalc(41)*(x2-0.2)*(1.0-(x2**m)/BL)

c           pb2L = pb2L/(1.0-x2/xm)**xvalc(44)
           
c          pb2L = min(1.0,pb2L)
c        endif            


c        if(x2.LE.3.4) then
c           pb2L = 0.1259*x2+2.1079*x2*x2-2.477*x2*x2*x2+
c     &           1.2025*x2**4-0.2718*x2**5+0.0235*x2**6
c        endif

c        pb2L =  1.0-0.3673*exp(-1.0*(0.6585*x2**1.25))-
c     &           0.6452*exp(-1.0*(1.861*x2**1.75))


c        pb2L = 1.0-xvalc(40)*exp(-1.0*(xvalc(41)*x2**1.75))-
c     &          xvalc(6)*exp(-1.0*(xvalc(44)*x2**2.5))

c        pb2L = pb2l*x2/(x2+0.02)

c        pb2L = 1.0-0.33733*exp(-0.57343*x2**1.5)
c     &       -0.7033*exp(-1.1765*x2**2.5)

c        pb2L = 1.0+0.011337*exp(0.19945*x2**1.75)
c     &       -1.0308*exp(-1.2663*x2**2.5)

cc        pb2L = 1.0-xvalc(40)*exp(-xvalc(41)*x2**1.75)
cc     &       -xvalc(6)*exp(-xvalc(44)*x2**2.5)
             
                
c        pb2L = 1.0-0.33733*exp(-0.57343*x2**1.5)
c     &            -xvalc(6)*exp(-xvalc(44)*x2**2.5)
        
c        pb2L = pb2L-0.035*exp(-1.0*(x2-2.45)**2.0/0.5/0.5)
        
c        pb2L = pb2L*(x2-0.15)**1.0/x2**1.0
        
c        pb2L = 1.0-0.03853*(4.0-x2)**2.5-0.020433*(4.0-x2)**3.5


        
CCC   Below is the nominal without the extra (4.0-x2)**1.5 term  CCC
        
        pb2L = 1.0-xvalc(40)*(4.0-x2)**2.5-xvalc(41)*(4.0-x2)**3.5        
        pb2L = pb2L-xvalc(6)*exp(-1.0*(x2-2.55)**2.0/xvalc(44)**2)
        pb2L = pb2L-xvalc(45)*(4.0-x2)**1.5
        pb2L = pb2L*(x2-0.18)**2/(x2-0.10)**2

CCC        


        if(x2.GT.4.0) es = 1.0
        es = min(es,1.0)
        es = max(es,0.0)      

        
c        pb2L = 1.0-0.16672*(4.0-x2)**2.5+0.14204E-01*(4.0-x2)**3.5        
c        pb2L = pb2L-0.17818*exp(-1.0*(x2-2.55)**2.0/0.62321**2)
c        pb2L = pb2L+0.30659*(4.0-x2)**1.5
c        pb2L  = pb2L*(x2-0.2)**2/(x2-0.12)**2


        
c        pb2L = 1.0-pb2L

        
c        pb2L = 1.0+0.023907*(4.0-x2)**2.5-0.018984*(4.0-x2)**3.5
c        pb2L = pb2L+0.051087*(4.0-x2)**1.5
c        pb2L = pb2L-0.11293*exp(-1.0*(x2-2.45)**2.0/0.45**2)
c        pb2L  = pb2L*(x2-0.2)**2/(x2-0.1)**2

ccc Next is Mahalia  ccc

c         pb2L = 1.0-0.00051412*(4.0-x2)**2.5-0.0020088*(4.0-x2)**4.5        
c         pb2L = pb2L-0.29654*exp(-1.0*(x2-0.91)**2.0/0.65**2)
c         pb2L  = pb2L*(x2-0.1)**2/(x2-0.07)**2

ccc        
        
        if(x2.GT.4.0) pb2L = 1.0
        pb2L = min(pb2L,1.0)

        
c        pb2L = 1.0-0.08785*exp(-1.0*(0.073615*x2**1.75))-
c     &          0.75167*exp(-1.0*(0.65935*x2**2.5))


c         if(x2.LT.0.6) write(6,*) "pb2L:  ", x2,pb2L

        
c        if(q2.GT.0.055.AND.q2.LT.0.08) pb2L = 1.08*pb2L
        
c        pb2L =  1.0-0.14189*exp(-1.0*5.5959*x2**1.5)-
c     &      0.88166*exp(-1.0*0.53714*x2**2.00)  
        
        if(psip.GT.psimax) FY = 0.0  !!! effective cutoff 


c        if(x2.LT.4.0) write(6,*) x2,pb2L

        
        FYL = FY
        FYL = pb2L*FYL

      elseif(paulitype.EQ.2) then    !!! Tsai - Fermi Gas
         FY = Pauli_sup2*FYp
      endif
     
      
       F2 = nu/kf * (FYL*nuL*GL + FY*nuT*GT)
       F1 = mp * FY/kf * GT / 2.
     

       if(F2.LT.0.0) F2 = 0.0
       if(F1.LT.0.0) F1 = 0.0
      
      return
      end


     
CCC-----------------

      
      SUBROUTINE MEC2021(z,a,w2,q2,xvalm,f1mec)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement in the QE and Delta region.         CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off at the single nucleon removal energy.      CCC
CCC                                                                             CCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
! fit to low q2 dip region: purefly empirical
! assume contribution is purely transverse
      implicit none
      real*8 Z,A,q2,w2,mp/0.938272/,mp2,mn/0.93957/,w,qv2,nu,numin
      real*8 a1,b1,c1,t1,dw2,xmax,q20,w2min
      real*8 x, f1mec, f1mec2, xvalm(45)

      xmax = 50.0
      q20 = 0.00001
c      q20 = xvalm(6)
      mp2 = mp*mp

      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu)
      qv2 = q2+nu**2.0
      numin = 0.0165
      w2min = mp2+2.0*mp*numin-q2
      xmax = q2/2.0/mp/numin
      
      
      if(A.lt.2.5) return


      a1 = A*(q2+q20)**2*xvalm(1)*exp(-1.0*q2*q2/xvalm(2))
     &     /(xvalm(3)+q2)**xvalm(4) 

      b1 = xvalm(5)
  
      c1 = xvalm(42)+xvalm(43)*q2

      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2-w2min

      if(dw2.LT.0.0) dw2 = 0.0
      f1mec = a1*exp(-1.0*t1)      
      f1mec = f1mec*(dw2)**1.5   !!!  check

      if(nu.LT.numin) f1mec = 0.0
      
      if(dw2.LE.0.0.OR.x.GT.xmax) f1mec = 0.0
      if(f1mec.LE.1.0E-9.OR.x.GE.xmax) f1mec=0.0


      return
      end


CCC-----------------    

      
C=======================================================================
                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, xvalc, F1, F2, FL)

CCC   Returns Fermi smeared structure functions.  Smearing function is a Gaussian.  CCC
CCC   Note:  not tested for nuclei with A < 12 or A > 64                            CCC      
CCC   July 20, 2020                                                                 CCC                      
!--------------------------------------------------------------------

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,xvalc(45),w2t,kappa2
      real*8 nu,x,mp,mp2,mpi,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,pf2,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,xp,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at,Lcor,frac,fract,t1
      real*8 xxp(500),fytot,fytot2,norm,bw,nwid,ncor

      real*8 emcfac,emcfacL,qvt
      logical goodfit
      INTEGER ISM,j,nbins

      nbins = 79
      nwid = 3.3
      bw = 2.0*nwid/float(nbins)
           

      mp = 0.938272
      mp2 = mp*mp
      mpi = 0.135
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)
      kappa2 = 1.0+4.0*mp2*x*x/q2
      
      if(A.GE.3.) then
      !!! energy shift !!!
        Es = 0.008
      !!! fermi momentum  !!!
        kf = xvalc(37)
        qvt = qv
        if(qv.GE.1.0) qvt = 1.0
        Es = xvalc(38)
        Es = Es-xvalc(39)*(1.0-qvt)   !!!  test  !!!
      endif

      norm = sqrt(pi)
      ncor = 1.000
c      ncor = 1.0027
c      ncor = 1.0663+4.7*(kf-0.225)
      norm = norm/ncor   !!! account for missing part of distribution past nwid for kf = 225 MeV


      f1p = 0.0D0
      f1n = 0.0D0
      f2p = 0.0D0
      f2n = 0.0D0
      fLp = 0.0D0
      fLn = 0.0D0


c        sigt = 0.0D0
c        sigL = 0.0D0
      fytot = 0.0D0
      fytot2 = 0.0D0


! adjust pf to give right width based on kf
      pf = 0.5 * kf
      pf2 = pf*1.5
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + mp) 

      DO ism = 1,nbins

CCC   
        xxp(ism) = -nwid+bw*(float(ism-1))

        
        fyuse = bw/sqrt(2.0)/norm*exp(-0.5*xxp(ism)*xxp(ism)) !!! Gaussian !!!       
        
CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des
        WSQP2 = W2 + XXp(ISM) * PF2 * DW2DPF - es * dw2des
        
CCC

        fytot = fytot+fyuse
        fytot2 = fytot2+fyuse

c           write(6,2000) w2,q2,ism,xxp(ism),fyuse, wsqp, fytot

        F1pp = 0.0D0
        F1nn = 0.0D0
        F2pp = 0.0D0
        F2nn = 0.0D0
        FLpp = 0.0D0
        FLnn = 0.0D0

        frac = 0.0
        
        do j=1,1
           if(j.EQ.1) then
              fract = 1.0D0-frac
              w2t = WSQP
           else
              fract = frac
              w2t = WSQP2
           endif   
           IF(w2t.GT. 1.159) THEN
             xt = q2/(q2+W2t-mp2)
             xp = 1.0D0+(w2t-mp)/(q2+xvalc(34))
             xp = 1.0D0/xp
          
             offshell = 1.0D0      !!!  test


CCC   Next is medium modification factor  CCC

            emcfac = (xvalc(26)+xvalc(27)*xp*xp)/
     &           (1.0+xvalc(28)*xp+xvalc(29)*xp*xp)

            emcfacL = 1.0     
            emcfacL = xvalc(30)*(1.0D0+xvalc(31)*xp*xp)*
     &       (1.+xvalc(32)*xp*xp)*exp(-1.0*xvalc(33)*xp)


            call sf(w2t,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)

c            t1 = (1.0+4.*xt*xt*mp2/q2)*f2pp-(2*xt*f1pp)

c            write(6,*) xt,q2,fLpp,t1         !!!  TEST
            
            f1pp = f1pp*emcfac*offshell
            f1nn = f1nn*emcfac*offshell
            fLpp = fLpp*emcfac*emcfacL*offshell
            fLnn = fLnn*emcfac*emcfacL*offshell
            f2pp = (2.*xt*f1pp+fLpp)/(1.+4.*xt*xt*mp2/q2)
            f2nn = (2.*xt*f1nn+fLnn)/(1.+4.*xt*xt*mp2/q2)

            F1p = F1p + F1pp * Fyuse * fract
            F1n = F1n + F1nn * Fyuse * fract
            F2p = F2p + F2pp * Fyuse * fract
            F2n = F2n + F2nn * Fyuse * fract      
            FLp = FLp + FLpp * Fyuse * fract
            FLn = FLn + FLnn * Fyuse * fract

         ENDIF
       ENDDO
      ENDDO

c      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

      F1 = (kappa2*F2-FL)/2.0/x  !!!  Calculate to keep internal consistency  !!!
      
      if(F1.LT.0.0) F1 = 0.0
      if(F2.LT.0.0) F2 = 0.0
      if(FL.LT.0.0) FL = 0.0

      
c      write(6,*) w2,f1p,f1n
      
c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          
      
      SUBROUTINE SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
CCCC   Converts reduced cross sections to structure functions for protons and neutrons  CCCCC 

      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp(w2,q2,sigTp,sigLp)
      call rescsn(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end

      
CCC   -----------------
      
      SUBROUTINE SFCROSS(w2,q2,A,Z,opt,sigt,sigl,F1,F2,FL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine to return reduced cross sections and structure functions for A>2 CCC
CCC   opt:  1 (total), 2 (QE only), 3 (inelastic only), 4 (MEC), 5 (QE+MEC)       CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      IMPLICIT none
 
      real*8 w2,q2,A,Z,sigt,sigl,F1,F2,FL
      real*8 x,alpha,pi,pi2,mp,mp2
      integer opt
      real*8 xval4(45) /
     & 0.34786E+00,0.11500E+01,0.14000E+00,0.46500E+01,0.86000E+00,
     & 0.69783E-01,0.21100E+01,0.18140E+01,0.23500E+00,-.30000E+01,
     & 0.99875E+00,0.97853E+00,0.99920E+00,0.10224E+01,0.10014E+01,
     & 0.10445E+01,0.10203E+01,0.10048E+01,0.91126E+00,0.99388E+00,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.13152E+01,0.28937E+00,0.21536E+00,0.35774E+00,0.17550E+02,
     & 0.35931E+00,0.20000E+01,0.84942E-01,0.11288E+00,0.10400E+02,
     & 0.40000E+01,0.11609E-01,0.38117E-06,0.26049E-01,0.58048E-01,
     & 0.24418E+00,0.23952E+00,0.21972E-05,0.12537E+01,0.10000E-01 /      
      real*8 xval12(45) / 
     & 0.94449E-01,0.10821E+02,0.13888E+00,0.69224E+01,0.78078E+00,
     & 0.11235E+00,0.20269E+01,0.20732E+01,0.70846E+00,-.41558E+01,
     & 0.99233E+00,0.98152E+00,0.10240E+01,0.99265E+00,0.10000E+01,
     & 0.10058E+01,0.97847E+00,0.10044E+01,0.99271E+00,0.99628E+00,
     & 0.10000E+01,0.99844E+00,0.10026E+01,0.10105E+01,0.10050E+01,
     & 0.78676E+00,0.00000E+00,-.10317E+01,0.92242E+00,0.20365E+01,
     & 0.19625E+01,0.19652E+01,0.28451E+01,0.50889E+00,0.22800E+00,
     & 0.11420E-01,0.27154E+00,0.40438E-01,0.75126E-01,0.44351E-01,
     & 0.71092E-02,0.78423E-01,0.29565E+00,0.62607E+00,-.13227E+00 /
      real*8 xval27(45) /
     & 0.23749E+00,0.11500E+01,0.14000E+00,0.46500E+01,0.86000E+00,
     & 0.12334E+00,0.21100E+01,0.18140E+01,0.23500E+00,-.30000E+01,
     & 0.10134E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10031E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.99667E+00,0.10022E+01,0.10148E+01,0.10000E+01,
     & 0.90139E+00,0.10786E+01,0.43890E+00,0.45142E+00,0.39084E+01,
     & 0.16780E+01,0.15000E+01,0.13519E+00,0.16512E+00,0.14431E+01,
     & 0.61780E+01,0.24876E-08,0.53818E-03,0.22976E-01,0.23540E-01,
     & 0.28401E+00,0.24839E+00,0.61620E-01,0.10025E+03,0.10000E-01 /        
      real*8 xval56(45) /
     & 0.20527E+00,0.11500E+01,0.14000E+00,0.46500E+01,0.86000E+00,
     & 0.14650E+00,0.21100E+01,0.18140E+01,0.23500E+00,-.30000E+01,
     & 0.10090E+01,0.10318E+01,0.10020E+01,0.99374E+00,0.10038E+01,
     & 0.10000E+01,0.99747E+00,0.10083E+01,0.99602E+00,0.10000E+01,
     & 0.10000E+01,0.99297E+00,0.10051E+01,0.10154E+01,0.10000E+01,
     & 0.91601E+00,0.10829E+01,0.51609E+00,0.48723E+00,0.12507E+01,
     & 0.15708E+01,0.15000E+01,0.17425E+00,0.17172E+00,-.59784E+00,
     & 0.44550E+01,0.24876E-08,0.15722E-01,0.23514E-01,0.47486E-01,
     & 0.28896E+00,0.26186E+00,0.21861E+00,0.10025E+03,0.10000E-01 /
      real*8 xval64(45) /
     & 0.20527E+00,0.11500E+01,0.14000E+00,0.46500E+01,0.86000E+00,
     & 0.14650E+00,0.21100E+01,0.18140E+01,0.23500E+00,-.30000E+01,
     & 0.10090E+01,0.10318E+01,0.10020E+01,0.99374E+00,0.10038E+01,
     & 0.10000E+01,0.99747E+00,0.10083E+01,0.99602E+00,0.10000E+01,
     & 0.10000E+01,0.99297E+00,0.10051E+01,0.10154E+01,0.10000E+01,
     & 0.91601E+00,0.10829E+01,0.51609E+00,0.48723E+00,0.12507E+01,
     & 0.15708E+01,0.15000E+01,0.17425E+00,0.17172E+00,-.59784E+00,
     & 0.44550E+01,0.24876E-08,0.15722E-01,0.23514E-01,0.47486E-01,
     & 0.28896E+00,0.26186E+00,0.21861E+00,0.10025E+03,0.10000E-01 /
       
      mp = .938272
      mp2 = mp*mp 

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi
      x = q2/(q2+w2-mp2)

c      write(6,*) "IN SFCROSS: ", w2,q2,a,z

      if(A.LE.7.0) then      !!! pick correct parameters for particular nucleus !!!     
        call csfitcomp(w2,q2,A,Z,xval4,opt,sigt,sigL)
      elseif(A.GT.6.0.AND.A.LE.20.0) then      !!! pick correct parameters for particular nucleus !!!     
        call csfitcomp(w2,q2,A,Z,xval12,opt,sigt,sigL)
      elseif(A.GT.20.0.AND.A.LE.44.0) then
        call csfitcomp(w2,q2,A,Z,xval27,opt,sigt,sigL)        
      elseif(A.GT.44.0.AND.A.LE.59.0) then
        call csfitcomp(w2,q2,A,Z,xval56,opt,sigt,sigL)        
      elseif(A.GT.59.0.AND.A.LE.80.0) then
         call csfitcomp(w2,q2,A,Z,xval64,opt,sigt,sigL)
      elseif(A.GT.80.0) then
         write(6,*) A
      endif
      
      f1 = sigt
      fL = sigL*2.0*x
      f2 = (2.0*x*f1+fL)/(1.+4.*mp2*x*x/q2)        !!! Fix this!!!!

      
      sigt = 0.3894e3*8.0d0*pi2*alpha/abs(w2-mp2)*sigt
      sigL =  0.3894e3*8.0d0*pi2*alpha/abs(w2-mp2)*sigL
      
      return
      end


CCC-----------------

      
      SUBROUTINE CSFITCOMP(w2,q2,A,Z,XVALC,type,sigt,sigL)
      IMPLICIT none

      real*8 e,ep,th,q2,w2,x,cs,flux,kappa2,sin2,tan2,csmod
      real*8 f1,f2,fl,f1qe,f2qe,flqe,f1mec,f2mec,fLmec,f1i,f2i,fLi
      real*8 r,rqe,sigt,sigl,sigm
      real*8 alpha,pi,pi2,mp,mp2,res,veff,foc,Z,A,xvalc(45)
      real*8 psip,psimax,psimin,fy1,fy2,int1,int2,rat,f1t,f2t,fLt,dpsi
      integer i,j,k,ntot,nbins,type
      LOGICAL GOODFIT/.true./  
      character*40 filename

      psimin = -2.3
      psimax = 5.0
      nbins = 220

      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      x = q2/abs(w2-mp2+q2)
      
      dpsi = (psimax - psimin)/float(nbins)

      kappa2 = (1.+4.*x*x*mp2/q2)   

      int1 = 0.0D0
      int2 = 0.0D0
      rat = 1.0D0

CCC  NEXT bit only needed if fitting scaling function CCC
      do i=1,nbins       
        psip = psimin+dpsi*(i-1)            
       FY1 = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       FY2 = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     > (1. + exp(xvalc(10) * psip))*(1.0-abs(psip)/psimax)**2.0* 
     > (1.0+abs(psip)/psimax)**2.0

       if(psip.GT.psimax) FY2 = 0.0
       FY2 = max(0.0,FY2)
       int1 = int1+fy1  
       int2 = int2+fy2
      enddo
      rat= 1.00/int2/dpsi
CCC


c      write(6,*) int1*0.04,int2*0.04,rat

c      call f1f2in09(Z,A,q2,w2,xvalc,f1t,f2t,r)

      call gsmearing(Z,A,w2,q2,xvalc,f1i,f2i,fLi)
      if(fLi.LT.0.0) FLi = 0.0

      
c      write(6,*) "gsmearing:  ", w2,q2, f1,f2,fL

c      call smearing(Z,A,w2,q2,xvalc,f1,f2,fL)
      
c      write(6,*) "smearing:  ",w2,q2, f1,f2,fL

      r = fL/2.0D0/x/f1

c      f1 = f1
c      f2 = f2
c      fL = fL

c      if(w2.GT.1.5) write(6,2000) w2,q2,f1,f1t,f2,f2t


      f1qe = 0.0
      f2qe = 0.0
      call qenuc21off(Z,A,q2,w2,xvalc,f1qe,f2qe)
      
      f1qe = f1qe*rat  !!! renormalize
      f2qe = f2qe*rat  !!! renormalize
      fLqe =  kappa2*f2qe-2.0*x*f1qe 
      if(fLqe.LT.0.0) flqe = 0.0
      

      f1mec = 0.0
      fLmec = 0.0
      call MEC2021(Z,A,w2,q2,xvalc,f1mec)
     
      f2mec = 2.0*x*f1mec/kappa2
      


      if(type.EQ.1) then
        f1 = f1i + f1qe + f1mec
        f2 = f2i + f2qe + f2mec
c     fL = kappa2*f2-2.0*x*f1
        fL = fLi+fLqe
      elseif(type.EQ.2) then
        f1 = f1qe
        f2 = f2qe
        fL = fLqe
      elseif(type.EQ.3) then
        f1 = f1i
        f2 = f2i
        fL = fLi
      elseif(type.EQ.4) then         
        f1 = f1mec
        f2 = f2mec
        fL = fLmec        
      endif
      if(fL.LT.0.0) fL = 0.0   
       
       sigt = f1
       sigl = fL/2./x

c      write(6,*) "2  ",w2,q2,x,sigt,sigl

 2000 format(6f10.4)
          
      return

      end
      




c--------------------------------- PROGRAM coulomb.F ------------------------------
c
c- P. Solvignon, Dec. 7, 2006
c
c  Coulomb corrections based on D. Gaskell's kumac (correct_emc.kumac)
c  
c  -- update --
c
c   Use method describe in Aste et al. : Eur. Phys. J. A26 (2005) 167
c      --> use the average Coulomb potential:
c          V = B*V0 with 0.75 < B < 0.80
c      --> the cross section is corrected by changing E-->E+V and Ep-->Ep+V
c          in the MOTT and in the spectral function expression
c      --> the corrected cross section is then multiplied by the incoming 
c          focusing factor FF squared: FF = (E+V)/E
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccc
c
      SUBROUTINE VCOUL(A,Z,V)
      IMPLICIT NONE

      REAL*8 A,Z,R0
      REAL*8 C_ASTE,HBARC,V0,V,ALPHA

      HBARC  = 0.197327      ! in GeV.fm
      ALPHA  = 1.0/137.0
      C_ASTE = 0.775
      R0     = 1.1*A**(1./3.) + 0.86*A**(-1./3.)

      ! Coulomb potential at the center of the nucleus
      V0  = (3./2.)*ALPHA*HBARC*(Z-1.)/R0  ! in GeV

      ! Average potential
      V  = C_ASTE*V0      ! from Eur. Phys. J. A26 (2005) 167

      RETURN
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
*
* $Id: simps64.F,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: simps64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
      FUNCTION DSIMPS(F,A,B,N2)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT real*8 (A-H,O-Z)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSIMPS')
*
* $Id: simpscod.inc,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: simpscod.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
*
* simpscod.inc
*
      DIMENSION F(0:*)
      IF(N2 .LE. 0 .OR. 2*(N2/2) .NE. N2) THEN
       H=0
*       WRITE(ERRTXT,101) N2+1
*       CALL MTLPRT(NAME,'D101.1',ERRTXT)
      ELSE
*      S1=0
       S1=F(N2-1)
       S2=0
*      DO 1 N = 1,N2-1,2
*   1  S1=S1+F(N)
*      DO 2 N = 2,N2-2,2
*   2  S2=S2+F(N)
       DO 1 N = 1,N2-3,2
         S1=S1+F(N)
         S2=S2+F(N+1)
    1  CONTINUE
*      S1=S1+F(N2-1)
       S1=S1+S1+S2
       H=(F(0)+F(N2)+S1+S1)*(B-A)/(3*N2)
      ENDIF
      DSIMPS=H
      RETURN
  101 FORMAT('NON-POSITIVE OR EVEN NUMBER OF FUNCTION VALUES =',I6)
      END

      
    


CCC-----------------



C **********************************************************************
	SUBROUTINE CDBONN (q,u0,u2)
C
C  Deuteron wave function from CD-Bonn NN potential model.
C  q in 1/fm, u0,u2 in fm^3/2.
C
C  Normalization \int dq q^2 (u0^2+u2^2) = 1.
C
C  Sent by Charlotte Elster, April 8, 2009.
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=95)
        REAL*8  q,u0,u2,
     &		qgrid(nq),uqgrid(nq),wqgrid(nq),weight(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hc
	SAVE

        pi = 4*DATAN(1.D0)
        hc = 197.327D0		! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
     &	       FILE='f1f2tables/cdbn.qwave',
     &	       STATUS='OLD')
c        READ (10,100)
        READ (10,*)
        READ (10,*)
C...Momentum space [qgrid in MeV, uqgrid in MeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
          qgrid(id) = qgrid(id) / hc            ! MeV => 1/fm
          uqgrid(id) = uqgrid(id) * hc**1.5D0   ! MeV^-3/2 => fm^3/2
          wqgrid(id) = wqgrid(id) * hc**1.5D0
        ENDDO

        
	init = .TRUE.

C...Evaluate wave function
c 999	u0 = DQDVAL (q,nq,qgrid,uqgrid,.FALSE.)
c	u2 = DQDVAL (q,nq,qgrid,wqgrid,.FALSE.)
 999	CALL Pinterp (qgrid,uqgrid,nq,q,u0,dum,2)
	CALL Pinterp (qgrid,wqgrid,nq,q,u2,dum,2)

        RETURN
        END


CCC-----------------
      
     
C ***********************************************************************
	SUBROUTINE AV18 (p,rho)
C
C  Deuteron momentum distribution rho = u^2 + w^2, where u and w are
C  S- and D-state wave functions in momentum space, normalized s.t.
C  int dp p^2 rho(p) = 1.
C  Ref: Fantoni and Pandharipande, Nucl. Phys. A427 (1984) 473
C
C  Input file has rho normalized s.t. 4*pi int dp p^2 rho(p) = 1,
C  so that for output, rho needs to be multiplied by 4*pi.
C
C  p in 1/fm, rho in 1/fm^3
C
C.. Uses IMSL interpolation routine DQDVAL.
C.. For compilation on jlabs1:
C.. > use imsl
C.. > f77 -o objectfile file.f -R/site/vni/lib/lib.solaris
C..       -L/site/vni/lib/lib.solaris -limsl -lsocket -lnsl
C.. IMSL decommissioned 8/30/11
C
C ***********************************************************************
        IMPLICIT NONE
        INTEGER	ip,np
        PARAMETER (np=200)
        REAL*8  p,rho
        REAL*8  Parr(np),RHOarr(np),rho_int,dum
        LOGICAL readin /.FALSE./
	REAL*8	pi
	SAVE

C...Value of pi
        pi = 4*DATAN(1.D0)

	rho = 0.D0

        IF (readin) GO TO 123
C...Read data from file
c        OPEN (10,FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/av18.dat',
        OPEN (10,FILE='f1f2tables/av18.dat',
     &		FORM='FORMATTED')
	DO ip=1,9
          READ (10,*)
	ENDDO
        DO ip=1,np
	  READ (10,*) Parr(ip), RHOarr(ip)
        ENDDO
        CLOSE (10)
        readin = .TRUE.
c        print *, '... AV18 data read ...'

 123	IF (p.LE.Parr(1)) rho = 4*pi * RHOarr(1)
	IF (p.GT.Parr(1) .AND. p.LE.Parr(np)) THEN
	  CALL Pinterp (Parr,RHOarr,np,p,rho_int,dum,2)
	  rho = 4*pi * rho_int
c     &    rho = 4*pi * DQDVAL(p,np,Parr,RHOarr,.FALSE.)
	ENDIF

c	print *, 'Parr(1),Parr(np)=',Parr(1),Parr(np)
c	print *, 'p,rho=',P, RHO

        RETURN
        END


CCC-----------------      
      

C **********************************************************************
	SUBROUTINE WJC1 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-1 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~105%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
c     &	      FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/wjc-1.dat',
     &	      FILE='f1f2tables/wjc-1.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-1 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END


CCC-----------------      

      
     
C **********************************************************************
	SUBROUTINE WJC2 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-2 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~102%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
     &	      FILE='wjc-2.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-2 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END


CCC-----------------

      

      
**************************************************************
***       File contains various interpolation codes        ***
**************************************************************
***      Code obtained from Alberto Accardi, Dec. 2008     ***
**************************************************************
*
* II  Polynomial function interpolation of given order
      subroutine pinterp(xa,ya,n,x,y,dy,order)
*     programmer: Alberto Accardi
*     date: 2/05/01
*
*  A. COMMENTARY
*
*     Performs an interpolation using a polynomial function
*     interpolation at a given order: given an x, it uses "order" points 
*     to its left and "order" to its right to perform the interpolation
*
*     xa(*) = (DP) array with tabulated abscissae (any dimension)
*     ya(*) = (DP) array with tabulated function  (any dimension)
*     n     = (I)  number of tabulated points  
*     x     = (DP) abscissa at which compute the interpolated function
*     y     = (DP) value of the function at x
*     dy    = (DP) estimated error (usually larger than real error)
*     order = (I)  order of the interpolation (see intro)  
*                  If order = 0 performs a linear interpolation
*                  between the nearest neighbours lattice point
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      integer n, order

      real*8 xa(*), ya(*), x, y, dy, tempx(n)
     :     , x1(2*order), y1(2*order), xmax, xmin, ymax, ymin  

      integer i, nlow, nmin

*
*  C. ACTION
*

      do i = 1, n
         tempx(i) = xa(i)
      end do
      call hunt(tempx,n,x,nlow)

      if (order.ge.1) then
         if (nlow.lt.order) then
            nmin = 0
         else if (nlow.le.n-order) then
            nmin = nlow-order
         else
            nmin = n-2*order
         end if
         do i = 1, 2*order
            x1(i) = xa(nmin+i) 
            y1(i) = ya(nmin+i) 
         end do
         call polintnum(x1,y1,2*order,x,y,dy)
      else
         ymax = ya(nlow+1)
         ymin = ya(nlow)
         xmax = xa(nlow+1)
         xmin = xa(nlow)
         y = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
      end if

      return
      end


************************************************************************
*
* III search in an ordered table 
      SUBROUTINE hunt(xx,n,x,jlo)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given an array xx(1:n) and given a value x, returns a value j
*     suchthat x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
*     either decreasing or increasing. j=0 or j=n is returned to
*     indicate that x is out of range.

*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER jlo,n

      real*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

*
*  C. ACTION
*

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo 
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      return
      END


CCC-----------------
      

************************************************************************
*
* IV  Polynomial interpolation and extrapolation
      SUBROUTINE polintnum(xa,ya,n,x,y,dy)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given arrays xa and ya of length n, and given a value x, this
*     routine returns a value y and an error estimate dy. If P(x) is the
*     polynomial of degree N-1 such that P(xa_i) = ya_i, i=1,...,n
*     then the returned value y = P(x).
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER n,NMAX
      real*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

*
*  C. ACTION
*

      if (n.gt.nmax) then
         print*, 'ERROR(polintnum): order larger than max', n,'>', nmax 
         stop
      end if
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c          if(den.eq.0.)pause 'failure in polintnum'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END



CCC-----------------     

      
!---------------------------------------------------------------------
      subroutine FNP_NMC(X,QSQ,rat)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,QSQ,A,B,X2,X3,rat
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      rat = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END
