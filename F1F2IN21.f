CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC                                                                      CCC
CCC   F1F220 version 0.995  -  February 20, 2021                         CCC
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
CCC   A > 2 nuclei fit are 4He, 12C, 27Al, and 56Fe.  More nuclei        CCC
CCC   will be fit in the next version.                                   CCC
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
      elseif(IA.EQ.2) then      !!! Deuteron
        eps = 0.5               !!! pass  0 < eps < 1
        wfn = 2                 !!! CD-Bonn, default
        doqe = .false.          !!! don't include QE  !!!

        call rescsd(WSQ,QSQ,eps,doqe,F1,F2,FL,wfn,sigm)
        call smearsub(WSQ,QSQ,wfn,off,F2,F1,FL)

c        write(6,*) "in f1f2in21: ", wsq,qsq,eps,wfn,f1,f2
        
        
      elseif(IA.GT.2) then      !!! A>2 nuclei
         opt = 3                !!! inelastic + MEC
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
      integer opt/1/
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
      
        GMn = mu_n*(1.0+0.12417E-04*tau)/
     &   (1.000+11.404*tau+17.014*tau*tau+31.219*tau**3)

       

        GEn = (1.5972*tau / (1.0+0.19655*tau)) * GD


       
       return
       end
      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
CCCC  Returns proton transverse and longitudinal cross sections   CCCC
CCCC  February 4, 2021                                            CCCC
      
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis,w2max,w2min

      data xvalp / 
     & 0.12287E+01,0.15196E+01,0.15044E+01,0.17009E+01,0.16765E+01,
     & 0.14455E+01,0.12472E+00,0.23000E+00,0.91162E-01,0.88525E-01,
     & 0.80876E-01,0.37521E+00,0.76849E+01,0.45411E+01,0.42458E+01,
     & 0.18299E+01,0.68593E+01,0.13435E-01,0.17015E+04,0.45076E+01,
     & 0.59724E+00,0.18944E+02,0.54924E-01,0.24830E+01,0.30381E+00,
     & 0.19419E+02,0.13587E+00,0.20332E+01,0.24378E+01,0.48652E+01,
     & 0.30850E+01,0.26876E+01,0.10685E-02,0.97878E+04,0.14193E+00,
     & 0.21888E+01,0.25592E+00,0.73812E+00,0.20187E+01,0.93656E-01,
     & 0.12867E-01,0.39477E+00,0.19371E+01,0.46900E+01,0.12106E-01,
     & 0.17154E+00,0.19929E+01,0.55000E+00,0.42980E+01,0.41895E+00,
     & 0.99122E+00,0.99085E+00,0.99798E+00,0.10028E+01,0.98145E+00,
     & 0.10163E+01,0.10226E+01,0.10165E+01,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.11353E-05,0.10605E+02,0.19705E+01,
     & 0.00000E+00,0.19429E-09,0.75543E+02,0.60997E+01,0.00000E+00,
     & 0.73839E+01,0.59820E-05,0.19247E+01,0.00000E+00,0.23239E+01,
     & 0.16384E+01,0.14220E+01,0.00000E+00,0.10015E-03,0.57993E+00,
     & 0.64963E+00,0.00000E+00,0.10093E-03,0.41844E+00,0.42011E+00,
     & 0.00000E+00,0.39842E+00,0.86421E+01,0.10000E+00,0.10000E+01,
     & -.82435E+01,0.48718E+02,0.51146E+02,0.29711E+01,0.11295E-02,
     & 0.47865E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.18935E+02 / 

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
        if(i.EQ.47.OR.i.EQ.48) xvalL(i) = xval1(i)
      enddo


      w2max = 34.0  !!! max for resonance fit
      w2min = 36.0  !!! min for dis fit

      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)

c      if(w2.GT.w2min) then 
c        call disp(w2,q2,sigTdis,sigLdis)
c        if(w2.LE.w2max) then
c         sigTp = (w2max-w2)*sigTp+(w2-w2min)*sigTdis
c         sigLp = (w2max-w2)*sigLp+(w2-w2min)*sigLdis
c         sigTp = sigTp/2.0
c         sigLp = sigLp/2.0
c         else
c          sigTp = sigTdis
c          sigLp = sigLdis
c        endif
c      endif

      return
      end






      SUBROUTINE RESCSD(w2,q2,eps,doqe,f1dqe,f2dqe,fLdqe,wfn,sigm)
CCCC  Calcultes deuteron cross sections from smeared p+n.             CCCC
CCCC  Requires SQESUB be called first to read in smearing functions.  CCCC
CCCC  This should be one at the lowest level to keep from reading     CCCC
CCCC  multiple times.

      IMPLICIT none

      real*8 w2,q2,q2t,eps,t1,x,gamma2,q2min
      real*8 sigtp,sigtn,siglp,sigln,sigt,sigl,sigm,m  
      real*8 m2,pi2,alpha,f1d,f2d,fLd,f1dqe,f2dqe,fLdqe
      real*8 xvaln(100),off_mKP_fit,delta,xfr,xfl
      integer i,j,k,ntssot,wfn,drn
      logical doqe,off
c      INCLUDE 'parm.cmn'
      external off_mKP_fit

      
      off = .false. !! off-shell corrections
c      off = .true.
      q2min = 4.0E-5
      q2t = q2
      if(q2.LT.q2min) q2t = q2min
      drn = 5
      xfr = 0.95
      xfl = 1.0E-3
      alpha = 1./137.036
      m = (0.938272+0.939565)/2.0d0  !!! average p,n
      m2 = m*m
      pi2 = 3.14158*3.14159

c      if(q2.LT.q2min) q2 = q2min  !!! Hack to fix photoproduction
      
      x = q2/(w2-m2+q2)
      gamma2 = 1.+4.*m2*x*x/q2
      t1 = gamma2*f2dqe-2.*x*f1dqe

      if(f2dqe.LT.0.0.OR.f1dqe.LT.0.0.OR.fLdqe.LT.0.0)
     &    write(34,*) w2,q2,f2dqe,f1dqe,fLdqe
        
      
      call SMEARSUB(w2,q2t,wfn,off,f2d,f1d,fLd)

      if(doqe) then 
       f1d = f1d+f1dqe
       f2d = f2d+f2dqe
       fLd = fLd+fLdqe
      endif

          
      sigt = 0.3894e3*f1d*pi2*alpha*8.0/abs(w2-m2)
      sigl = 0.3894e3*fLd*pi2*alpha*8.0/abs(w2-m2)/2.*abs(w2-m2+q2)/q2
      sigm = sigt+eps*sigl


      return
      end






      SUBROUTINE RESCSN(w2,q2,sigtn,sigln)
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xval1(50),xvalL(50)
      integer i
      real*8 xvaln(100)/
     & 0.12287E+01,0.15195E+01,0.15042E+01,0.17062E+01,0.16787E+01,
     & 0.14474E+01,0.12502E+00,0.23000E+00,0.90468E-01,0.86480E-01,
     & 0.75000E-01,0.38012E+00,0.71001E+01,0.32634E+01,0.10730E+01,
     & 0.22650E+01,0.23212E+01,0.18872E+01,0.37882E+01,0.12554E+01,
     & 0.29554E+01,0.95547E+01,0.12848E+02,0.27110E+01,0.97586E+03,
     & 0.46182E+03,0.25324E+02,0.47504E+05,0.27562E+00,0.20067E+02,
     & 0.23779E+00,0.23150E+01,0.28648E+01,0.12912E+04,0.63948E+01,
     & 0.61034E+02,0.33488E+00,0.16043E+01,0.24234E+02,0.12297E+00,
     & 0.10214E+00,0.52532E+00,0.22671E+01,0.48200E+01,-.77397E-01,
     & 0.15751E-01,0.19872E+01,0.54078E+00,0.24385E+01,0.18916E-01,
     & 0.10198E+01,0.97431E+00,0.99305E+00,0.98196E+00,0.10425E+01,
     & 0.10123E+01,0.98465E+00,0.98416E+00,0.10284E+01,0.98698E+00,
     & 0.10162E+01,0.10125E+01,0.11713E+02,0.94676E+02,0.60923E+01,
     & 0.46433E-03,0.29789E-01,0.21407E+02,0.29234E+01,0.11574E+02,
     & 0.52793E+00,0.41040E+01,0.16980E+01,0.16856E+02,0.35216E+01,
     & 0.10865E+02,0.29016E+01,0.31103E+02,0.11927E+02,0.26535E+02,
     & 0.81055E+01,0.16015E+01,0.28755E-01,0.24493E+01,0.83532E+00,
     & 0.19400E+00,0.15667E+00,0.33667E+01,0.21528E+01,0.10000E+01,
     & 0.29325E+02,0.56383E+01,0.73032E+01,0.72695E+01,0.18496E+00,
     & 0.78960E+00,0.19850E+01,0.44998E+00,0.10037E+01,0.28875E+01 /
      

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






      SUBROUTINE RESMODP(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from February 21, 2021  -  Author:  M.E. Christy               CCC
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
         x0(i) = 0.178   !!! 
      enddo

c      x0(1) = 0.125
      x0(1) = 0.142
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.05)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.05)**1.
           
      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
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
        A0 = (1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude
        t1 = xval(38)*log(2.0+q2/xval(39))                 !!! exponent of (1-xpr)
        t2 = xval(40)*log(2.0+q2/xval(41))+xval(45)        !!! exponent of xpr


        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = xval(37)*389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif
       
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         
         t1 = xval(38)*log(1.001+q2/q20)                   !!! exponent of (1-xpr)
         t2 = xval(41)/(1.001+q2/xval(42))**xval(43)       !!! exponent of xpr

        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + xval(37)*389.4*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2
        endif
        sig_nr =sig_nr/log(1.001+q2/xval(39))
        
      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 



      SUBROUTINE RESMODN(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from February 21, 2021  -  Author:  M.E. Christy               CCC
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
         x0(i) = 0.178   !!! 
      enddo

c      x0(1) = 0.125
      x0(1) = 0.142
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.05)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.1)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
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
        A0 = (1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude
        
        t1 = xval(38)*log(2.51+q2/xval(39))                !!! exponent of (1-xpr)
        t2 = xval(40)*log(2.51+q2/xval(41))+xval(45)       !!! exponent of xpr
        
        
        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = xval(37)*389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif
       
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         
         t1 = xval(38)*log(1.001+q2/q20)                   !!! exponent of (1-xpr)
         t2 = xval(41)/(1.001+q2/xval(42))**xval(43)       !!! exponent of xpr

        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + xval(37)*389.4*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2
        endif
        sig_nr =sig_nr/log(1.001+q2/xval(39))


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

      if(w2.LT.1.6) nbins = 80

c      nbins = 200   !!! adjust as needed

      if(w2.GE.2.6.AND.Q2.GE.2.0) nbins = 46
      
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
      gamma2 = (1. + 4.*mp2*x*x/q2)
      write(6,'(" F1F2IN21: gamma2=",8f8.3)') gamma2,q2,w2,mp2,
     >  x,tau,q
      if(gamma2.LE.6.and.firsty) then   
       gamma = sqrt(gamma2)
       write(6,'('' F1F2IN21: calling getfy'',l2)') firsty
       call getfy(gamma,x,wfn,fy11,fy12,fy2,firsty)
       firsty = .true.
       write(6,'('' F1F2IN21: calling getfyoff'',l2)') firsty
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


        off = 1   !!! CC1
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

cdg      write(6,2001) x,q2,gm2,ge2

c      off = 2

      f1os = 0.0D0
      f2os = 0.0D0

      nbins = 1000

c      do i=1,10
c        wfnorm(i) = 1.0D0
c      enddo

cdg      write(6,*) w2,q2,wfn

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
      a = sqrt(1.0 + w2d/qsqr)
      b = w2d/(2.0*mp*sqrt(qsqr))
      prod = b/(a*a - 1.D0)
      c = 1.D0+((a*a-1.0D0)*(1.0-a*a/b/b))
 
      if(c.GE.0.0D0) then
        rt = sqrt(c)
      else
        return
      endif

      ff = (GE2+tau*GM2)/(1.D0+tau)
      pff = ((sqrt(GE2)-sqrt(GM2))/(1.D0+tau))**2      !!! Pauli FF squared

      pvmin = mp*prod*(1.d0-rt)*sign(1.0,(a-b))
      pvmax = mp*prod*(1.D0+rt)
      pvinc = (pvmax-pvmin)/float(nbins)

      write(6,*) "F1F2IN21: offshellqe:  ",w2,q2,x,rt

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
          open(unit=34,file='f11.WBARELav18',status='old')
          open(unit=35,file='f2.WBARELav18',status='old')
          open(unit=36,file='f12.WBARELav18',status='old')
        elseif(wfn.EQ.2) then
          open(unit=34,file='f1f2tables/f11-onshell-cdbonn.dat',status='old')        
          open(unit=35,file='f1f2tables/f22-onshell-cdbonn.dat',status='old')
          open(unit=36,file='f1f2tables/f12-onshell-cdbonn.dat',status='old')
        elseif(wfn.EQ.3) then
          open(unit=34,file='f11.WBARELwjc1',status='old')        
          open(unit=35,file='f2.WBARELwjc1',status='old')
          open(unit=36,file='f12.WBARELwjc1',status='old') 
        elseif(wfn.EQ.4) then
          open(unit=34,file='f11.WBARELwjc2',status='old')        
          open(unit=35,file='f2.WBARELwjc2',status='old')
          open(unit=36,file='f12.WBARELwjc2',status='old')
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
       write(6,'('' F1F2IN21: read in '',i6,
     >    '' points in getfy'')') i

      endif
      firsty = .false.
      close(34)
      close(35)
      close(36)

c      write(6,'(//''getfy '',4e12.3)') 
c     >  yv12(12,800),fyv12(12,800),
c     >  yv(12,800),fyv11(12,800)
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
          open(unit=34,file='f11.WBARELav18OFF',status='old')        
          open(unit=35,file='f2.WBARELav18OFF',status='old')
          open(unit=36,file='f12.WBARELav18OFF',status='old')
        elseif(wfn.EQ.2) then  !!!!  Not available yet - fix later
          open(unit=34,file='f1f2tables/f11-offshell-cdbonn.dat',status='old')        
          open(unit=35,file='f1f2tables/f22-offshell-cdbonn.dat',status='old')
          open(unit=36,file='f1f2tables/f12-offshell-cdbonn.dat',status='old')
        elseif(wfn.EQ.3) then
          open(unit=34,file='f11.WBARELwjc1OFF',status='old')        
          open(unit=35,file='f2.WBARELwjc1OFF',status='old')
          open(unit=36,file='f12.WBARELwjc1OFF',status='old') 
        elseif(wfn.EQ.4) then
          open(unit=34,file='f11.WBARELwjc2OFF',status='old')        
          open(unit=35,file='f2.WBARELwjc2OFF',status='old')
          open(unit=36,file='f12.WBARELwjc2OFF',status='old')
        else
          write(6,*)  "Warning:  wavefunction undefined"
        endif

     
        dowhile(.not.thend)   !!!  First read in smearing function to temp vector!!!
          i = i+1
          j = int(float(i-1)/999.)+1 
          k = i-(j-1)*999

          read(34,*,end=999) gammavoff(j),yvoff(j,k),fyv11off(j,k)
          read(35,*,end=999) gammavoff(j),yv2off(j,k),fyv2off(j,k)
          read(36,*,end=999) gammavoff(j),yv12off(j,k),fyv12off(j,k)

        enddo
999     thend = .true. 
       write(6,'('' F1F2IN21: read in '',i6,
     >    '' points in getfyoff'')') i
      endif
      firsty = .false.
      close(34)
      close(35)
      close(36)

c      write(6,'(///''getfyoff '',4e12.3)') 
c     >  yv12off(12,800),fyv12off(12,800),
c     >  yvoff(12,800),fyv11off(12,800)
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

      SUBROUTINE QENUC21(Z, A, qsq, wsq, xvalc,F1, F2)

C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS  
c for A>2. Uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c based on the earlier code F1F2QE09 by P. Bosted.  The superscaling distribution
c shape is determined from the fit to 12C data      
c
c input: Z, A  (real*8) Z and A of nucleus (should be 2.0D0 for deuteron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invariant mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus


      IMPLICIT NONE     

      REAL*8 Z, A, avgN, F1, F2, wsq, qsq, R
      REAL*8 amp/0.938272/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 nu, qv, TAU, FY
      real*8 kappa, lam, lamp, taup, squigglef, psi, psip,nuL,nut
      real*8 kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      real*8 xvalc(45),mcor,ecor
      integer IA
      
      real*8 a1,a2,a3,b1,b2,b3,b4,b5,gd


! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2= 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)


!     Call New FormFactors  !
      
      call formfacts(qsq,gmp,gep,gmn,gen)

!  Get energy shift and fermi momementum from fit  !      
      
      if(IA.GE.2) then
         Es = 0.008
         kf = 0.230
         Es= xvalc(37)
         kf = xvalc(41)
         if(qv.LT.1.0) Es = Es-xvalc(38)*(1.0-qv)   !!!  test  !!!           
      endif 

      
! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip)) / kf

CCMEC - Fitted Superscaling distribution      CCCC

      FY = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip)) / kf

      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0.0.and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else    !!!  This should not happen, but just in case  !!!
        r = 0.4/qsq
      endif

      if(F2.LT.0.0) F2 = 0.0
      
      return
      end

     
CCC-----------------

      
      
      SUBROUTINE MEC2021(z,a,w2,q2,xvalm,f1mec)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement new the QE and Delta due to meson   CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off near the 2-body threshold near x=2.        CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2
      real*8 x, f1mec, xvalm(40)

      mp2 = mp*mp
      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(A.lt.2.5) return

      a1 = A*q2**2*xvalm(1)*exp(-1.0*q2/xvalm(2))/
     &                                (xvalm(3)+q2)**xvalm(4)

      b1 = xvalm(5)+xvalm(6)*q2

      c1 = xvalm(33)+xvalm(34)*q2

      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1mec = a1*(exp(-1.*t1)*sqrt(dw2))

      
       if(f1mec.LE.1.0E-9 ) f1mec=0.0


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
      real*8 Z,A,q2,w2,f1,f2,fL,xvalc(45),w2t
      real*8 nu,x,mp,mp2,mpi,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,xp,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at,Lcor,frac,fract
      real*8 xxp(100),fytot,fytot2,norm

      real a4, x4
      real*8 fitemct,emcfac,emcfacL,xfr,xfl,qvt
      logical goodfit
      INTEGER ISM,drn,wfn,j
c      external off_mKP_fit

      frac = xvalc(45)
c      frac = 0.10
      
      drn = 5
      xfr = 0.95
      xfl = 1.E-3
      wfn = 2

      mp = 0.938272
      mp2 = mp*mp
      mpi = 0.135
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)
      a4 = A                                                    
      x4 = x

      if(A.GE.3.) then
      !!! energy shift !!!
        Es = 0.008
      !!! fermi momentum  !!!
        kf = xvalc(42)
        qvt = qv
        if(qv.GE.1.0) qvt = 1.0
        Es = xvalc(39)
        if(qv.LT.1.0) Es = Es-xvalc(40)*(1.0-qv)   !!!  test  !!!
      endif

      norm = 20.471
      norm = norm*2.0
c      norm = 1.0

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
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + mp) 

      DO ism = 1,99

CCC   
        xxp(ism) = -3.0+6.0*float(ism-1)/98.0
        
        fyuse = 1.0/norm*exp(-0.5*xxp(ism)*xxp(ism))    !!! Gaussian !!!
c        fyuse2 = 0.2/norm*exp(-0.5*xxp2(ism)*xxp2(ism))    !!! Gaussian !!!
        
CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des
        WSQP2 = W2 + XXp(ISM) * xvalc(32) * PF * DW2DPF - es * dw2des
        
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

        do j=1,2
           if(j.EQ.1) then
              fract = 1.0D0-frac
              w2t = WSQP
           else
              fract = frac
              w2t = WSQP2
           endif   
           IF(w2t.GT. 1.159) THEN
             xt = q2/(q2+W2t-mp2)
             xp = 1.0D0+(w2t-mp)/(q2+xvalc(43))
             xp = 1.0D0/xp
          
            offshell = 1.0D0      !!!  test

            x4 = xt

CCC   Next is medium modification factor  CCC

            emcfac = (xvalc(26)+xvalc(27)*xp)*exp(-1.0*xvalc(29)*xp)/
     &               (1.0+xvalc(28)*xp*xp)
            emcfac = (xvalc(26)+xvalc(27)*xp)/
     &       (1.+xvalc(28)*xp*xp)*exp(-1.0*xvalc(29)*xp)

          
            emcfacL = xvalc(31)*(1.0D0+xvalc(35)*xp)*
     &       (1.+xvalc(36)*xp*xp)*exp(-1.0*xvalc(30)*xp)

            Lcor = 1.0+xvalc(30)*sqrt(q2)/(1.0+xvalc(44)*q2*q2)
c            Lcor = 1.0D0
c            Lcor = 1.0/(1.0-1.0*exp(-1.0*xvalc(44)*q2*q2))
            emcfacL = Lcor*emcfacL
            
c          write(6,*) xp,x,emcfac,emcfacL
          
c          emcfacL = emcfacL*(1.+0.05*sqrt(xp))       

            call sf(w2t,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)
    
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

      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

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
CCC   opt:  0 (total), 1 (QE only), 2 (inelastic only), 3 (inelastic + MEC)       CCC
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
     & 0.26094E+00,0.11500E+01,0.14000E+00,0.46500E+01,0.86000E+00,
     & 0.76283E-01,0.21100E+01,0.18140E+01,0.23500E+00,-.30000E+01,
     & 0.97620E+00,0.95995E+00,0.10276E+01,0.98852E+00,0.98207E+00,
     & 0.10000E+01,0.97724E+00,0.10080E+01,0.99097E+00,0.99642E+00,
     & 0.99214E+00,0.99592E+00,0.10034E+01,0.10141E+01,0.10077E+01,
     & 0.89978E+00,0.10817E+01,0.61362E+00,0.44255E+00,0.21842E+01,
     & 0.14732E+01,0.10006E+01,0.11335E+00,0.17477E+00,0.40469E+00,
     & 0.43605E+01,0.23172E-02,0.10001E-06,0.19460E-01,0.28553E-01,
     & 0.27439E+00,0.24915E+00,0.23916E+00,0.10025E+03,0.10000E-01 /
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
        call csfit(w2,q2,A,Z,xval4,opt,sigt,sigL)
      elseif(A.GT.6.0.AND.A.LE.20.0) then      !!! pick correct parameters for particular nucleus !!!     
        call csfit(w2,q2,A,Z,xval12,opt,sigt,sigL)
      elseif(A.GT.20.0.AND.A.LE.44.0) then
        call csfit(w2,q2,A,Z,xval27,opt,sigt,sigL)        
      elseif(A.GT.44.0.AND.A.LE.59.0) then
        call csfit(w2,q2,A,Z,xval56,opt,sigt,sigL)        
      elseif(A.GT.59.0.AND.A.LE.80.0) then
         call csfit(w2,q2,A,Z,xval64,opt,sigt,sigL)
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

      
      
      SUBROUTINE CSFIT(w2,q2,A,Z,xvalc,opt,sigt,sigL)
CCC   Returns reduced transverse and longitudnal cross section for A > 2 nuclei CCC
CCC   Requires fit parameter array xvalc.                                       CCC
CCC   Can return parts or all contributions using                               CCC
CCC   opt = 1 (QE)
      
      IMPLICIT none

      real*8 q2,w2,x
      real*8 f1,fl,f1qe,r,sigt,sigl,f1mec,f2,f2qe,f2mec
      real*8 alpha,pi,pi2,mp,mp2,Z,A,xvalc(45)
      real*8 psip,fy1,fy2,int1,int2,rat
      integer i,opt
c      character*40 filename
      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      x = q2/abs(w2-mp2+q2)
 
      int1 = 0.0D0
      int2 = 0.0D0
      rat = 1.0D0

CCC  NEXT bit only needed if fitting scaling function CCC
      do i=1,120        
        psip = -2.0+0.06*i
       FY1 = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       FY2 = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip))
       int1 = int1+fy1  
       int2 = int2+fy2
      enddo
      rat= int1/int2   !!! normalization of new scaling function if fit  !!!
CCC

      
      call gsmearing(Z,A,w2,q2,xvalc,f1,f2,fL)

      r = fL/2.0D0/x/f1
      
      call QENUC21(Z,A,q2,w2,xvalc,f1qe,f2qe)
            
      call MEC2021(Z,A,w2,q2,xvalc,f1mec)
          
      f1qe = f1qe*rat         !!! Normalize new scaling function if fit
      f2qe = f2qe*rat
    
      if(opt.EQ.1) then       !!! QE only
         f1 = 0.0
         f2 = 0.0
         fL = 0.0
         f1mec = 0.0
      elseif(opt.EQ.2) then   !!! inelastic only
         f1qe = 0.0
         f2qe = 0.0
         f1mec = 0.0
      elseif(opt.EQ.3) then   !!! inelastic + MEC
         f1qe = 0.0
         f2qe = 0.0
      endif
      
      f1 = f1 + f1qe + f1mec
      f2mec = 2.*x*f1mec/(1.+4.*x*x*mp2/q2) 
      f2 = f2 + f2qe + f2mec
      fL = (1.+4.*x*x*mp2/q2)*f2-2.0*x*f1
      if(F1.LE.0.0) write(6,*) "F1 < 0"
      if(F1.LE.0.0) F1=0.0
      if(F2.LE.0.0) F2=0.0
      if(FL.LE.0.0) FL=0.0
      sigt = f1
      sigl = fL/2./x

          
      return

      end
      


CCC-----------------

*
* $Id: simps64.F,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: simps64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
      real*8 FUNCTION DSIMPS(F,A,B,N2)
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
        OPEN (10,FILE='av18.dat',
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
     &	      FILE='wjc-1.dat',
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

      real*8  xa(*), ya(*), x, y, dy, tempx(n)
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
c      write(6,'(''debug'',4i4)') n,nmax,ns
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








