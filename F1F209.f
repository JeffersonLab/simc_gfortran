ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c F1F209.f
c Package of FORTRAN subroutine describing fits to inclusive inelastic
c electron scattering from proton, neutron, deuteron, or heaviere nuclei
c Proton fit is described in:
c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision 
c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
c    (arXiv:0712.3731). Submitted to Phys. Rev. C.
c Deuteron (and netron) fit is described in:
c    P.E. Bosted and M.E. Christy, ``Empirical Fit to Inelastic 
c    Electron-Deuteron and Electron-Neutron
c    Resonance Region Transverse Cross Sections, 
c    (arXiv:0711.0159), publichished in Phys. Rev. C 77, 065206 (2008). (
c New fits for A>2 by Vahe M. and Peter B. (to be publsihed).
c
c routine in this package:
c Main program test. Compiling and linking all code
c    in this package will make a "test" program.
c    Upon execution, the output should be:
c  1.000  1.500  0.617  0.204  0.231  0.312  0.410  1.644  0.072  0.102
c  1.000  2.000  0.472  0.129  0.159  0.292  0.250  2.329  0.199  0.191
c  1.000  2.500  0.382  0.169  0.227  0.417  0.360  2.973  0.241  0.292
c  1.000  3.000  0.321  0.224  0.315  0.511  0.500  3.562  0.327  0.359
c  1.000  3.500  0.276  0.218  0.309  0.517  0.521  4.060  0.348  0.387
c  2.000  1.500  0.763  0.079  0.091  0.110  0.104  0.400  0.154  0.188
c  2.000  2.000  0.641  0.058  0.091  0.157  0.104  0.764  0.184  0.216
c  2.000  2.500  0.553  0.090  0.139  0.241  0.162  1.174  0.215  0.246
c  2.000  3.000  0.485  0.137  0.214  0.324  0.252  1.558  0.233  0.270
c  2.000  3.500  0.433  0.136  0.227  0.354  0.272  1.919  0.267  0.283
c   1.000   0.747   0.432
c   2.000   0.157   0.132
c   3.000   0.048   0.049
c
c
c subroutine F1F2IN09. Returns inelastic F1, F2, and R for
c      nucleus with charge Z and atomic number A
c      for given value of Q2 and W**2. F1 and F2
c      are per nucleus (not per nucleon).
c subroutine pind. Returns F1, R, Sigma_L, and Sigma_T
c      for a proton with Fermi motion of a detueron
c subroutine resd. Returns F1 for a deuteron.
c subroutine resder. Returns error on F1 from fit.
c      Requires auxillary file "F1F207D2emat.dat"
c subroutine resmodd. Returns F1 for average of a free
c      proton and neutron
c subroutine christy507. Returns F1, R, sigma_L, and
c      sigma_T for the free proton.
c subroutine resmod507. Returns sigma_L or sigma_T
c      (called by christy507)
c subroutine mec. Called by resd to get extra terms
c      for dip region between quasi-elastic and Delta
c funtion fitemc. Fit to "EMC effect" used to 
c      get F1 and F2 for A>2
c 
c subroutine F1F2QE09. Returns quasi-elastic F1, F2 for
c      nucleus with charge Z and atomic number A
c      for given value of Q2 and W**2
c 
cc      implicit none
cc      integer iq,iw
cc      real*8 q2,w2,F1n,F2n,r,F1p,F2p,F1d,F2d,F1c,F2c
cc      real*8 f1dqe, f2dqe, nu, x, am/0.9383/, F1be, F2be,rbe
cc
cc      do iq=1,2
cc       q2 = 1.0 * iq
cc       do iw=1,5
cc        w2 = 1. + 0.5*iw
cc        nu = (w2 - am**2 + q2) / 2. / am
cc        x = q2 / 2. / am / nu
cc        call F1F2IN09(0.D0, 1.D0, q2, w2, F1n, F2n,r)
cc        call F1F2IN09(1.D0, 1.D0, q2, w2, F1p, F2p,r)
cc        call F1F2IN09(1.D0, 2.D0, q2, w2, F1d, F2d,r)
cc        call F1F2IN09(4.D0, 9.D0, q2, w2, F1be, F2be,rbe)
cc        write(6,'(10f7.3)') q2,w2,x,f2n,f2p,f2d,
cc     >    f1p,f1be,r,rbe
cc       enddo
cc      enddo
cc      do iq=1,3
cc        q2 = 1.0 * iq
cc        w2 = am**2
cc        call F1F2QE09(1.D0, 2.D0, q2, w2, F1dqe, F2dqe)
cc        write(6,'(3f8.3)') q2, F1dqe, F2dqe
cc      enddo

cc      stop
cc      end



C=======================================================================
                                                                        
      SUBROUTINE F1F2IN09(Z, A, QSQ, Wsq, F1, F2)                       
!--------------------------------------------------------------------
! Fit to inelastic cross sections for A(e,e')X
! valid for all W<3 GeV and all Q2<10 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 10/20/2006 P. Bosted
!--------------------------------------------------------------------
      implicit none
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 Z,A,qsq,wsq,w,f1c,x
      real*8 avgn,r,dr,nu,eps,kappa,sigres,flux,siglp,sigtp,F1pp,F1dp
      REAL*8 W1,W2,sigt,rc,w1pb,w2pb,F1,F2,sigl,F1d, F1p,qv
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des,fyuse
      real a4, x4, fitemc, emcfac
      logical goodfit
      INTEGER ISM

! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/

! This is for exp(-xx**2/2.), from teste.f
       real*8 xxp(99)/
     > -3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
     > -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
     > -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
     > -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
     > -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
     > -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
     >  0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
     >  0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
     >  1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
     >  1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
     >  2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000/
! these are 100x bigger for convenience
       real*8 fyp(99)/
     > 0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
     > 0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
     > 0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
     > 0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
     > 1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
     > 2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
     > 2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
     > 1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
     > 0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
     > 0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
     > 0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272/

      integer iz,ia,i
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/

! deuteron fit parameters
       real*8 xvald0(50)/
     >  0.1964E+01, 0.1086E+01, 0.5313E-02, 0.1265E+01, 0.8000E+01,
     >  0.2979E+00, 0.1354E+00, 0.2200E+00, 0.8296E-01, 0.9578E-01,
     >  0.1094E+00, 0.3794E+00, 0.8122E+01, 0.5189E+01, 0.3290E+01,
     >  0.1870E+01, 0.6110E+01,-0.3464E+02, 0.9000E+03, 0.1717E+01,
     >  0.4335E-01, 0.1915E+03, 0.2232E+00, 0.2119E+01, 0.2088E+01,
     > -0.3029E+00, 0.2012E+00, 0.1104E-02, 0.2276E-01,-0.4562E+00,
     >  0.2397E+00, 0.1204E+01, 0.2321E-01, 0.5419E+03, 0.2247E+00,
     >  0.2168E+01, 0.2266E+03, 0.7649E-01, 0.1457E+01, 0.1318E+00,
     > -0.7534E+02, 0.1776E+00, 0.1636E+01, 0.1350E+00,-0.5596E-02,
     >  0.5883E-02, 0.1934E+01, 0.3800E+00, 0.3319E+01, 0.1446E+00/

cc     
       real*8 F1M
       logical DEBUG/.TRUE./

      IA = int(A)
      avgN = A - Z
      nu = (wsq - pm**2 + qsq) / 2. / pm
      qv = sqrt(nu**2 + qsq)

      if(Wsq.le.0.0) W = 0.0
      W  = sqrt(Wsq)
      x  = QSQ / (2.0 * pm * nu)
      if(Wsq.le.0.0) x = 0.0

! Cross section for proton or neutron
      W1 = 0.
      W2 = 0.
      IF(IA .lt. 2 .and. wsq.gt.1.155) THEN
        call CHRISTY507(Wsq,Qsq,F1p,Rc,sigt,sigl)
! If neutron, subtract proton from deuteron. Factor of two to
! convert from per nucleon to per deuteron
        if(Z .lt. 0.5) then
          call resmodd(wsq,qsq,xvald0,F1d)
          F1p = F1d * 2.0 - F1p
        endif
        W1 = F1p / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      ENDIF

! For deuteron
      if(IA .eq. 2) then
c get Fermi-smeared R from Erics proton fit
        call pind(Wsq, Qsq, F1c, Rc, sigt, sigl)
c get fit to F1 in deuteron, per nucleon
        call resd(qsq, wsq, xvald0, F1d)
c convert to W1 per deuteron
        W1 = F1d / PM * 2.0
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      endif

! For nuclei
      IF(IA.gt.2) then
        sigt = 0.
        sigl = 0.
        F1d = 0.
        F1p = 0.
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
       if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
! changed 4/09
        if(IA.eq.3) kf=0.115
        if(iA.eq.3) Es=0.001 
! changed 4/09
        if(IA.gt.3) kf=0.19
        if(iA.gt.3) Es=0.017
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
c changed 5/09
        if(iA.gt.7) Es=0.0165
        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
! changed 5/09 
        if(iA.gt.55) Es=0.018 
 

! adjust pf to give right width based on kf
        pf = 0.5 * kf 
! assume this is 2 * pf * qv
        DW2DPF = 2. * qv
        dw2des = 2. * (nu + PM) 
! switched to using 99 bins!
cc        do ism = 1,15
cc          fyuse = fy(ism)
cc          WSQP = WSQ + XX(ISM) * PF * DW2DPF - es * dw2des
        do ism = 1,99
          fyuse = fyp(ism)/100.
          WSQP = WSQ + XXp(ISM) * PF * DW2DPF - es * dw2des
          IF(WSQP.GT. 1.159) THEN
            call CHRISTY507(Wsqp,Qsq,F1pp,Rc,sigtp,siglp)
            call resmodd(wsqp,qsq,xvald0,F1dp)
            F1d = F1d + F1dp * Fyuse
            F1p = F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse
          ENDIF
        ENDDO
        Rc = 0.
        if(sigt .gt. 0.) Rc = sigl / sigt
        W1 = (2. * Z * F1d + (A - 2. * Z) * (2. * F1d - F1p)) / PM 

        W1= W1*(1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)

cc        if(W .GT. 1.3) 
        if(W .GT. 0.0) 
     >       W1=W1*(1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ))**2

        CALL MEC2009( Z , A , qsq , wsq , F1M )

        W1 = W1 + F1M
        if(Wsq .gt.0.0 ) Rc = Rc * ( 1.0 + P(6) + P(23)*A )
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)

        DEBUG=.FALSE.
        IF( W1 .LE. 0.0 .AND. DEBUG ) THEN 
           write(*,*) 'test  = ', Z,A,W,QSQ,x,F1M,W1
           write(*,*) 'test1 = ', 
     >          (1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ)),
     >        (1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)
        ENDIF

      ENDIF

      A4 = A
      x4 = qsq / 2. / pm / nu
      emcfac = fitemc(x4, a4, goodfit)

      F1 = pm * W1 * emcfac 
      F2 = nu * W2 * emcfac 

      RETURN                                                            
      END                                          

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MEC2009(z,a,q2,w2,f1)

! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,am/0.9383/,w,nu
      integer i
      real*8 pb(20)/ 
     >     0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >     0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >    -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >     0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >     0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 p18

      real*8 x, f1corr

      f1 = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am
      x  = q2 / (2.0 * am * nu )

      if(a.lt.2.5) return

      p18 = p(18)
! special case for 3He
      if(a.gt.2.5 .and. a.lt.3.5) p18 = 70
! special case for 4He
      if(a.gt.3.5 .and. a.lt.4.5) p18 = 170.
! new values for C, Al, Cu
      if(a.gt.4.5) p18 = 215.
      if(a.gt.20.) p18 = 235.
      if(a.gt.50.) p18 = 230.

       
       f1corr = P(0)*exp(-((W-P(1))**2)/(P(2)))/ 
     >      ((1.0 + MAX( 0.3 , Q2 ) / P(3) ) ** P(4) )*nu**P(5)
     >      *( 1.0 + P18 * A ** ( 1.0 + P(19) * x ) )

       f1 = f1corr

       if(f1 .le.1.0E-9 ) f1=0.0
c       write(*,*) 'vahe1= ', A, W*W, Q2, f1corr

      return
      end

      SUBROUTINE F1F2QE09(Z, A, qsq, wsq, F1, F2)
c
C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      integer IA, izz, izzmin, izp, izznom, izdif

c Look up tables for deuteron case
       real*8 fyd(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

c     Peter Bosted's correction params

       real*8 pb(20)/ 0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >      0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >     -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >      0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/ 
       real*8 y,R

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
! changed 4/09
      if(IA.eq.3) kf=0.115
      if(iA.eq.3) Es=0.001 
! changed 4/09
      if(IA.gt.3) kf=0.19
      if(iA.gt.3) Es=0.017 
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
c changed 5/09
        if(iA.gt.7) Es=0.0165
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
! changed 5/09 
        if(iA.gt.55) Es=0.018 


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

! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izz)))**2 - 
c   >      qv**2 + 2. * qv * pz - avp2(izz)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'',3i4,6f7.3)')
     >     izznom,izzmin,izz,qsq,wsq,w2p,dwmin/1.e6,pzmin,pznom
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0. .and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else
        r = 0.4/qsq
      endif


c apply correction factors
      if(A.gt.2) then
         y = (wsq -amp**2) / qv
c         F1 = F1 * (1. + pb(8) + pb(9) * y +
c     >        pb(10)*y**2 + pb(11)*y**3 + pb(12)*y**4 )
c         R = R * (1. + pb(13))
c         F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)

cc correction to correction Vahe
         if(wsq.gt.0.0) then

            F1=F1*(1.0+P(7)+P(8)*y+P(9)*y**2 +P(10)*y**3 +P(11)*y**4)
            R = R * ( 1.0 + P(12) )
            F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)
            if(F1.LT.0.0) F1=0.0

         endif
      endif

      return
      end
      BLOCK DATA CORRECTION
      REAL*8 P(0:23)
      COMMON/PARCORR/P
c     DATA P/
c    c       5.1141e-02,   9.9343e-01,   5.3317e-02,   1.3949e+00, 
c    c       5.9993e+00,  -2.9393e-01,   9.9316e-02,   1.0935e-02, 
c    c       3.7697e-01,   1.3542e+00,  -7.4618e+00,   4.3540e+00, 
c    c      -3.7911e-01,   3.5105e-01,  -1.8903e+00,   4.9139e+00, 
c    c      -5.9923e+00,   2.5021e+00,   1.9943e+01,  -5.5879e-02, 
c    c       3.1914e-01,  -1.8657e-01,   3.2746e+01,   4.9801e-03/
       DATA P/
     c       5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
     c       6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
     c       4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
     c      -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
     c      -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
     c       2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03/
       END


      subroutine pind(W2,Q2,F1,R,sigt,sigl)
! Calculate proton with Fermi smearing of a deuteron 
      implicit none
      real*8 q2,w2,F1,R,sigt,sigl,am/0.9383/,nu,qv,F1p,Rp,sigtp,siglp
      real*8 amd/1.8756/,w2p,pz
      integer ism
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1=0.
      R=0.
      sigt=0.
      sigl=0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
         w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fyd(ism) / 10.
          sigl = sigl + siglp * fyd(ism) / 10.
          F1   = F1   + F1p   * fyd(ism) / 10.
        endif
       enddo
      else
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
c Need avp2f term to get right behavior x>1! 
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fydf(ism) / 100.
          sigl = sigl + siglp * fydf(ism) / 100.
          F1   = F1   + F1p   * fydf(ism) / 100.
        endif
       enddo
      endif

      if(sigt.ne.0.) R = sigl / sigt
      return
      end
      

      subroutine resd(q2,w2,xval,F1)
! Calculate dueteron F1 by Fermi smearing of proton plus neutron 
! Add MEC term (not smeared)
c    P.E. Bosted and M.E. Christy, ``Empirical Fit to Inelastic 
c    Electron-Deuteron and Electron-Neutron
c    Resonance Region Transverse Cross Sections, 
c    (arXiv:0711.0159). Submitted to Phys. Rev. C.
      implicit none
      real*8 q2,w2,xval(50),F1,am/0.9383/,nu,qv,dw2dpf,w2p,sigp,f1sv
      real*8 sigtst,amd/1.8756/, pz, f1m,f2m
      integer ism,i
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/
      logical usemec/.true./

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1 = 0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
        w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fyd(ism) / 10.
        endif
       enddo
      else
       f1 = 0.
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
! Need avp2f term to get right behavior x>1!
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fydf(ism) / 100. 
        endif
       enddo
      endif
! add MEC term
! took out again, then put back in again
      call mec(1.D0,2.D0,q2,w2,f1m,f2m,xval)
      if(usemec) f1 = f1 + f1m
      return
      end

      subroutine resder(q2,w2,xval,F1er)
! Get error on F1 from deuteron fit
      implicit none
      integer i,j,icall
      real*8 q2,w2,xval(50),xvalp(50),deriv(50),emat(50,50),F1,F1er
      real*8 epsilon,F1p
      integer ifix(50)/1, 2, 3, 4, 5, 6, 7, 0, 0, 0,
     >                 0, 0, 8, 9,10,11,12,13,14,15,
     >                16,17,18,19,20,21,22,23,24,25,
     >                26,27,28,29,30,31,32,33,34,35,
     >                36,37,38,39,40,41, 0, 0,42,0/  
      logical first/.true./
      
      icall = icall + 1
      if(first) then
        open(unit=9,file='F1F207D2emat.dat')
        do i=1,50
          read(9,'(10E12.4)') (emat(i,j),j=1,50)
        enddo
        close(unit=9)
        first = .false.
      endif

      call resd(q2,w2,xval,F1)
      epsilon = 1.E-7
      do i=1,50
        do j=1,50
          xvalp(j) = xval(j)
          if(i.eq.j) xvalp(j) = xvalp(j) + epsilon
        enddo
        call resd(q2,w2,xvalp,F1p)
        deriv(i) = (F1p - F1) / epsilon
        if(icall.eq.10) write(6,'(i3,3e12.4)') i,f1,f1p,deriv(i)
      enddo
      F1er=0.
      do i=1,50
        do j=1,50
          if(ifix(i).ne.0 .and.ifix(j).ne.0) then
            F1er = F1er + emat(ifix(i),ifix(j)) * deriv(i) * deriv(j)
          endif
        enddo
      enddo
! add overall Hall C norm. err
      F1er = sqrt(F1er + (0.02*F1)**2)
      return
      end


! returns F1 for average of free proton and neutron
! for given W2, Q2
      SUBROUTINE RESMODD(w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(5000,7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,xpr,alpha,pi,F1
      INTEGER i,j,l,lmax,num,iw
      real*8 noverp,fnp_nmc,x,a,b,sig_mec,brp(7,3)
      real*8 xval0(12)/
c new 5/07 values. Values 1 and 7 will be overridden below.
     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00/
      real*8 xvalold(50),w2sv,q2sv,sigrsv(7),md,w2p,wp,wdifp,xprp,nu
      logical first/.true./
      common/tst2/sigrsv,sig_nr,sig_mec
      real br2(7),br3(7)
      save

      sig = 0.
      if(w2.lt.1.07327**2 .or. w2.gt.25 .or. 
     >  q2.lt.0.0 .or. q2.gt.11.0) then
        write(6,'(1x,''error, q2 or w2 out of range'',2f8.3)') w2,q2
        return
      endif

c do this if fitting masses or widths, else set first true in above
      first=.false.
      if(xvalold(1).ne.xval(1) .or.
     >   xvalold(7).ne.xval(7)) first=.true.
      xvalold(1)=xval(1)
      xvalold(7)=xval(7)

      if(first) then
       mp = 0.9382727
       mpi = 0.135
       mpi2 = mpi*mpi
       meta = 0.547
       mp2 = mp*mp
       pi = 3.141593
       alpha = 1./137.036

! branching ratios
       br(1,1) = 1.0     
       br(2,1) = 0.5
       br(3,1) = 0.65
       br(4,1) = 0.65
       br(5,1) = 0.4
       br(6,1) = 0.65
       br(7,1) = 0.6

! angular momenta
       ang(1) = 1.       !!!  P33(1232)
       ang(2) = 0.       !!!  S11(1535)
       ang(3) = 2.       !!!  D13(1520)
       ang(4) = 3.       !!!  F15(1680)
       ang(5) = 0.       !!!  S15(1650)
       ang(6) = 1.       !!!  P11(1440) roper   
       ang(7) = 3.       !!!  ? 4th resonance region

! x0 parameter
       do i=1,7
c 2006
c         x0(i) = 0.165
c 2007
         x0(i) = 0.215
       enddo
c 2006
c      x0(4) = 0.6
c 2007
       x0(1) = 0.1446

! out branching ratio
       do i=1,7
         br(i,2) = 1.-br(i,1)
       enddo
    
! remember w2
       w2sv = w2

! uses xvals of 1-12, 47, and 48
! move masses, wdiths into local variables
! pyb changed to be fixed
       num = 0
       do i=1,6              
         num = num + 1
         mass(i) = xval0(i)
       enddo
       do i=1,6             
         num = num + 1
         intwidth(i) = xval0(num)
       enddo
! changed to allow delta width, mass to vary
! taken out again since xval(1) used in MEC
c       mass(1) = xval(1)
c       intwidth(1) = xval(7)
c 2006
c       mass(7) = xval(47)
c       intwidth(7) = xval(48)
c 2007
       mass(7) = 1.9341
       intwidth(7) = 0.380

! precalculate w-dependent quantites in 0.1 MeV bins
       do iw=1073,5000
        w = 0.001 * (iw+0.5)
        w2 = w**2
        wdif = w - (mp + mpi)
        wr = wdif/w

! Calculate kinematics needed for threshold Relativistic B-W 
        k = (w2 - mp2) / 2. / mp
        kcm = (w2 - mp2) / 2. / w
        epicm = (W2 + mpi**2 -mp2 ) / 2. / w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 ) / 2. / w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 ) / 2. / w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))
        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))
! Calculate partial widths
          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
          if(i.ne.2) then
c            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
c     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)
c make same as resmod507
            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     >       **(ang(i)+2.)
     >       * W / mass(i)
          else
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif
          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)
          pgam(i) = intwidth(i)*pgam(i)
          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(iw,i) = width(i) * pgam(i) / ((W2 - mass(i)**2.)**2. 
     &            + (mass(i)*width(i))**2.) *
     >            kr(i) / k * kcmr(i) / kcm / intwidth(i)
        enddo ! loop on i
c        write(55,'(i5,f7.2,7f10.4)') iw,w,(sigr(iw,i),i=1,7)
       enddo ! loop on iw

       w2 = w2sv
       first = .false.
       do i=1,50
         xvalold(i) = xval(i)
       enddo
! write table for article
c       open(unit=9,file='nhd.tbl1')
c       do i=1,7
c         br2(i)=br(i,2)
c         br3(i)=0.
c         if(i.eq.2) br2(i)=0.
c         if(i.eq.2) br3(i)=br(i,2)
c         if(i.le.6) then
c          write(9,199) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),(xval(12 + (i-1)*4 + j),j=1,4)
c 199      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & ',f7.1,' & ',f6.3,' & ',f6.3,' \\\\')
c         else
c          write(9,198) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),xval(49)
c 198      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & 0.0 & 0.0 & 4.0 \\\\')
c         endif
c       enddo
c       close(unit=9)
c       open(unit=9,file='nhd.tbl2')
c       write(9,197) (xval(i),i=2,6)
c       do i=1,2
c         write(9,197) (xval(36 + (i-1)*4 + j),j=1,4),xval(44+i)
c 197     format(f7.1,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',
c     >     f7.4,' \\\\')
c       enddo
c       write(9,'(''xval50='',f10.4)') xval(50)
c       close(unit=9)
      endif ! if first
      
! get parameters into local variables
      num = 12
! resonance height coefficients. xvals of 13-36
      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo
      enddo
!  Non-Res coefficients xvals of 37-44
      do i=1,2               
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo

! Begin resonance Q^2 dependence calculations   CCC
! uses xvals 49
      do i=1,6
        height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2) * q2 / (1. + rescoef(i,3) * q2))/
ccc     &          (1. + q2/0.71)**rescoef(i,4)
c make same as resmod507
     &          (1. + q2/0.91)**rescoef(i,4)
      enddo
ccc      dip = 1./(1. + q2 / 0.71)**2  
c make same as resmod507
      dip = 1./(1. + q2 / 0.91)**2  
      dip2 = dip**2
ccc      height(7) = xval(49)*dip2 
c make same as resmod507
      height(7) = xval(49)*dip 
      iw = int(1000.*sqrt(w2))
      sig_res = 0.
      do i=1,7
ccc        sigrsv(i) =  height(i) * sigr(iw,i)
c make same as resmod507 by squaring height
        sigrsv(i) =  height(i)**2 * sigr(iw,i)
        sig_res = sig_res + sigrsv(i) 
      enddo

c make same as resmod507
      sig_res = sig_res * sqrt(w2)

! Begin non-resonant part uses xvals 45, 46, 50
! Depends on both W2 and Q2 so can't easily precalculate
      sig_nr = 0.
c      xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
! to make same as resmod507
      xpr = 1.+(w2-(mp+mpi)**2)/(q2+0.05)
      xpr = 1./xpr
      w = sqrt(w2)
      wdif = w - (mp + mpi)
      do i=1,2  
        sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &   /(q2+nr_coef(i,2))**
     &   (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
      enddo
c new section for sqrt(W)
c make same as resmod507 by turning this off
ccc        sig_nr = sig_nr +(xval(2)*(wdif)**0.5)
ccc     &   /(q2 + xval(3))**
ccc     &   (xval(4) + xval(5) * q2 + xval(6) * q2**2)
        
      sig_nr = sig_nr * xpr
     
! Add third term to try to describe MEC, now using Wdiff in 
! deuteron rather than proton
! ** Taken out 10/17/06 
c     md = 2.*mp
c     nu = (q2 + w2 - mp2) / 2. / mp
c     w2p = md**2 + 2. * md * nu - q2
c     Wp = sqrt(w2p)
c     wdifp = wp - md
c     sig_mec = 0.
c     if(wdifp .gt. 0.) then
c       xprp = 1. + (w2p - (md)**2) / (q2 + xval(50))
c       xprp = 1. / xprp
c       sig_mec = (xval(1) + xval(2)*wdifp**(1/2.) +
c    >     xval(3)*wdifp) /
c    &    (q2 + xval(4))**(xval(5) + xval(6) * q2) *
c    >    xprp
c      endif
c     sig = sig_res + sig_nr + sig_mec

      sig = sig_res + sig_nr 
c      write(6,'(1x,i4,2f7.2,4e10.3)') iw,q2,w2,height(1),
c     >  sigr(iw,1),sig_res,sig_nr

! changed to use F1 instead of sigt
      F1 = sig * (w2-mp2)/8./pi/pi/alpha/0.3894e3
      sig = F1

      RETURN 
      END 

c Fit to proton F1, R, sigma_T, and Sigma_L from

      SUBROUTINE christy507(W2,Q2,F1,R,sigt,sigl)
c   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision 
c    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
c    (arXiv:0712.3731). To be submitted to Phys. Rev. C.

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end

      SUBROUTINE RESMOD507(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
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
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

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

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

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
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

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
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       


 1000  format(8f12.5)

      RETURN 
      END 


!---------------------------------------------------------------------

      REAL FUNCTION FITEMC(X,A,GOODFIT)       
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) 
! with no isoscalerity correction
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was includded with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
!    beyond that assume Fermi motion is giving the rise, and we
!    already are taking that into account with the y-smearing of
!    the inelastic
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*4 ALPHA, C,LN_C,X,A ,X_U
      REAL*8 XU8
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
      fitemc = 1.
      if(A .lt. 2.5) return

      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
       IF(X.GT. 0.70) X_U = 0.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                
      XU8 = X_U
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(XU8))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END                                                               

      subroutine mec(z,a,q2,w2,f1,f2,xval)
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,f2,xval(50),am/0.9383/,w,nu

      w = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am

! changed to use max(0.3,q2)
      f1 = xval(1) * exp(-(w - xval(2))**2/xval(3)) /
     >   (1. + max(0.3,q2) / xval(4))**xval(5) * nu ** xval(6)
      
      f2 = 0.
      if(q2.gt.0.) f2 = nu * (f1/am) / (1. + nu**2 / q2) 
      return
      end

