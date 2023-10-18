	real*8 function peepiX(vertex,vertex0,main,survivalprob,doing_cent)


* Sept. 2003 D. Gaskell
* Purpose:
* This routine calculates p(e,e'pi+)X semi-iinclusive cross sections from
* the CTEQ5M parton distribution functions and a simple parameterization
* of the favored and unfavored fragmentation functions.
*   output:
*	sigma_eepiX	!d3sigma/dEe'dOmegae'Omegapi	(microbarn/MeV/sr^2)
* 
* 
* October 3 2003 D. Gaskell
* Replace simple fragmentation function paramterization with a slightly
* more sophisticated treatment from Binnewies et. al.,PRD 52, p.4947 (1995).
* This parameterization gives: D_u^(pi+ + pi-) = D_d^(pi+ + pi-) = D+ + D-
* The separate favored and unfavored fragmentation functions are given by
* the ratio D-/D+ from HERMES data (my fit to P. Geiger's results).
* The D-/D+ fit is valid from z=0.25 to 0.9 or so, but behaves pretty well
* at lower z (approaches simple form of Field and Feynman (1-z)/(1+z)).
*
*
* October 4 2003 D. Gaskell
* Add PI- functionality
*
* October 9 2003 D. Gaskell
* Add deuterium functionality.  Assume just an incoherent sum of parton
* distributions from proton and neutron.  Assume isospin symmetry so
* d_neutron(x) = u_proton(x)
* dbar_neutron(x) = ubar_proton(x)
* u_neutron(x) = d_proton(x)
* ubar_neutron(x) = dbar_proton(x)
*
*
* March 2004 D. Gaskell
* Add calculation of "central" cross section.  This is convenient for 
* bin centering.
*
* April 8, 2004 D. Gaskell
* Add conversion for jacobian so in ntuple: siglab/jacobian = dsig/dE dOe dz dPt2 dPhi
* (units are MeV everywhere)
* Add strange quark contributions and kaon fragmentation
* April 15, 2004 D. Gaskell
* Put in Fermi motion effects.  Note that for now, if you use Fermi motion
* the "central cross section" calculation stuff probably doesn't make much sense.
*
*
* September 2020 D Gaskell
* Use F1F209 for inclusive contribution
* Add Peter Bosted's fit for pion fragmenetation functions from Hall C 2018-2019 data



	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex, vertex0

* PDFs
	integer iset  !which set (1=cteq5m)
	integer ipart !particle u=1, ubar=-1, d=2, dbar=-2, s=3, sbar=-3
	real*8 u,d,ubar,dbar,s,sbar
	real*8 qu,qd,qs ! u, d, and s quark charges
	real*8 D_fav, D_unfav, D_sum, R_D !favored,unfavored,sum,ratio of FFs
	real*8 D_sum_s, D_s  ! strange frag. functions
	real*8 lambda, Q2zero ! scales for FF param
	real*8 sv !scaling variable for FF param
	real*8 N,a1,a2 !parameters for FF param
	real*8 Ns, a1s, a2s !parameters for strange FF param

C Some local kinematic variables
	real*8 xbj,sx, Q2gev, Qgev, pt2gev,mtargev,nugev !unitless or GeV
	real*8 b  ! pt2 parameter for FFs

	real*8 nu,qx,qy,qz,mtar,Q2,Eb,Eprime,Epx,Epy,Epz  !MeV
	real*8 pt2,zhad,Ehad,phad,mhad ! all in MeV
	real*8 cthpq

	real*8 kcent,klo,khi
	integer i

	real*8 sum_sq, dsigdz, sigsemi, jacobian, fac, sigma_eepiX
	real*8 sighad, sige

	real*8 q2qe,w2qe
	integer wfn
	real*8 F1,F2,FL,W1,W2,sin2th2,cos2th2,W2coeff
	real*8 Ctq5Pdf	
c	external Ctq5Pdf

C Variables for kaon decay stuff
	real*8 survivalprob  ! this will be included in main.weight
	real*8 zaero,pathlen,p_kaon,betak,gammak

c parameters for PB fit of 9/11/2020
c versus zp
c       real*8 pf(8)/   1.2803,   0.0851,   0.8379,   0.1586,
c     >                 0.0140,   0.2133,  -4.4985,   4.1285/
c       real*8 pu(8)/   0.8290,  -0.1416,   0.9869,   0.2559,
c     >                 0.0090,  -1.2306,  -1.5292,   2.4169/
C parameters for PB fit of 9/20/2021
	real*8 pf(12)/ 1.0424,  -0.1714,   1.8960,  -0.0307,
     >                 0.1636,  -0.1272,  -4.2093,   5.0103,
     >                 2.7406,  -0.5778,   3.5292,   7.3910/
	real*8 pu(12)/ 0.7840,   0.2369,   1.4238,   0.1484,
     >                 0.1518,  -1.2923,  -1.5710,   3.0305,
     >                 1.1995,   1.3553,   2.5868,   8.0666/

	real*8 xp,zp,yf,yu
	real*8 Mpi_gev, Mp_gev, wsq,w, mmpi2,z8, a8
	real*8 d1,db,u1,ub, s1, sb, dsigdzn, dsigdzp

	logical first, firstqe
	logical doing_cent,first_cent! flag for "central" cross section calc.

	parameter (iset=1)
	parameter (qu=2./3.)
	parameter (qd=-1./3.)
	parameter (qs=-1./3.)

c	parameter (b=3.76)   !GeV^-2
c	parameter (b=4.68)   !GeV^-2

	parameter (lambda=0.227) !0.227 GeV for NLO
	parameter (Q2zero=2.0)   !Gev^2 for u,d,s,g

	data first /.TRUE./
	data firstqe /.TRUE./
	data first_cent /.TRUE./

c this is for DSS
        integer IHdss,ICdss,IOdss, fini 
        real*8  C1, B1, GL1
        COMMON / FRAGINI / FINI

	if(first) FINI=0
c	b = pt_b_param   ! now parameter in input file

	Mpi_gev = Mpi/1000.0
	Mp_gev = Mp/1000.0
C DJG: Setup stuff for doing "central" cross section calculation.  Here, I'm
C DJG: assuming we want the cross section at some "point" in Q2 and W space.
C DJG: I allow for binning in either z or pt2 - holding the other constant.
C DJG: Since the cross section has no phi-dependence, that part is easy - 
C DJG: I can just ignore it.  

	if(doing_cent) then
	   nu = vertex0%nu
	   qx = vertex0%uq%x*vertex0%q
	   qy = vertex0%uq%y*vertex0%q
	   qz = vertex0%uq%z*vertex0%q
	   Q2 = vertex0%Q2
	   Eb = vertex0%Ein
	   Eprime = vertex0%e%E
	   Epx = vertex0%ue%x*Eprime
	   Epy = vertex0%ue%y*Eprime
	   Epz = vertex0%ue%z*Eprime
	   if(sigc_flag.eq.0) then  ! binning in z
	      kcent = vertex0%zhad
	      pt2 = sigc_kin_ind*1.e6
	      zhad = 0.0
	      do i=1,sigc_nbin
		 klo = sigc_kin_min+(i-1)*(sigc_kin_max-sigc_kin_min)/sigc_nbin
		 khi = sigc_kin_min+i*(sigc_kin_max-sigc_kin_min)/sigc_nbin
		 if(vertex%zhad.gt.klo .and. vertex%zhad.le.khi) then
		    zhad = klo + (khi-klo)/2.
		 endif
	      enddo
	   elseif (sigc_flag.eq.1) then  ! binning in pt2
	      kcent = vertex0%pt2
	      zhad = sigc_kin_ind
	      pt2 = 0.0
	      do i=1,sigc_nbin
		 klo = sigc_kin_min+(i-1)*(sigc_kin_max-sigc_kin_min)/sigc_nbin
		 klo = klo*1.e6
		 khi = sigc_kin_min+i*(sigc_kin_max-sigc_kin_min)/sigc_nbin
		 khi = khi*1.e6
		 if(vertex%pt2.gt.klo .and. vertex%pt2.le.khi) then
		    pt2 = klo + (khi-klo)/2.
		 endif
c		 write(6,*) 'chessy poofs',pt2,klo,khi
	      enddo
	   endif

	   if(first_cent) then
	      write(6,*) 'Central Kinematics:'
	      write(6,*) 'Ebeam (GeV):',Eb/1000
	      write(6,*) 'nu (GeV)   :',nu/1000
	      write(6,*) 'Q2 (GeV2)  :',Q2/1e6
	      if(sigc_flag.eq.0) then
		 write(6,*) 'Pt2 (GeV2) :',pt2/1e6
		 write(6,*) 'Binning in z from',sigc_kin_min,
     >  	      'to',sigc_kin_max
	      elseif(sigc_flag.eq.1) then
		 write(6,*) 'z       :',zhad
		 write(6,*) 'Binning in pt2 from',sigc_kin_min,
     >  	      'to',sigc_kin_max,'GeV2'
	      endif
	      first_cent=.false.
	   endif
	else
	   nu = vertex%nu
	   qx = vertex%uq%x*vertex%q
	   qy = vertex%uq%y*vertex%q
	   qz = vertex%uq%z*vertex%q
	   Q2 = vertex%Q2
	   Eb = vertex%Ein
	   Eprime = vertex%e%E
	   Epx = vertex%ue%x*Eprime
	   Epy = vertex%ue%y*Eprime
	   Epz = vertex%ue%z*Eprime
	   pt2 = vertex%pt2
	   zhad = vertex%zhad
	endif

	if(doing_semipi) then
	   mhad = Mpi
	else if(doing_semika) then
	   mhad = Mk
	else
	   write(6,*) 'semi_physics error: unknown hadron type'
	   stop
	endif

	mtar = targ%Mtar_struck

	Ehad = zhad*nu
	phad = sqrt(Ehad**2-mhad**2)

	if(doing_cent) then
	   cthpq = sqrt(1.-pt2/phad**2)
	else
	   cthpq = cos(vertex%theta_pq)
	endif

	sx = (2.*Eb*mtar + mtar**2)/1.e6  !convert to GeV2
	if(do_fermi) then  ! xbj = Q2/(2 P.q)
	   xbj = Q2/2./(efer*nu - abs(pfer)*(pferx*qx + pfery*qy + pferz*qz))
	   if(.not.doing_cent) then
	      ntup%xfermi = xbj
	   endif
	else
	   xbj = Q2/2./mtar/nu   
	endif
	if(xbj.gt.1.0) then
	   write(6,*) 'XBj is too large!', xbj
	   xbj=1.0
	endif

c	if(do_fermi) then ! y = P.q/P.k_beam
c	   y = (efer*nu - pferx*qx - pfery*qy - pferz*qz)/(efer*Eb)
c	else
c	   y = nu/Eb
c	endif

C DJG convert some stuff to GeV

	Q2gev = Q2/1.e6
	Qgev = sqrt(Q2gev)
	pt2gev = pt2/1.e6

	wsq = Mp_gev**2 + q2gev * (1./xbj -1.)
	w = sqrt(wsq)
c added that xsection should be zero below 2pi / piK  threshold (PB)
C this shouldn't be necessary, but include it anyway
	mtargev = mtar/1000.
	nugev = nu/1000.
        mmpi2 = mtargev**2 + 2. * mtargev * nugev * 
     >        (1-zhad) * (1 - pt2gev)
	if(mmpi2 .lt. (mtargev + mhad/1000.)**2) then
         sigma_eepiX = 0.
	 peepiX = 0.0
	 return
	endif


c needed by f1f2in21
	if(firstqe) then
	 write(6,*) 'Initializing F1F2IN21:'  
	 write(6,*) 'F1F2IN21: calling sqesub'
	 q2qe=1.
	 w2qe=1.
	 wfn=2
         call sqesub(w2qe,q2qe,wfn,f1,f2,fL,firstqe)
	 firstqe=.false.
	endif

C Get the PDFs
	if(first) then
	   call SetCtq5(iset)	! initialize Cteq5 (we're using cteq5m)
	   first=.FALSE.
	endif

	ipart=1
	u = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-1
	ubar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=2
	d = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-2
	dbar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=3
	s = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-3
	sbar = Ctq5pdf (ipart , xbj, Qgev)

	sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar) + qs**2*(s+sbar)
	   
	if (doing_deutsemi) then
	   sum_sq = sum_sq + qu**2*(d+dbar) + qd**2*(u+ubar) + qs**2*(s+sbar)
	endif


C Simple paramaterization from Kretzer et al (EPJC 22 p. 269)
C for Q2=2.5.

c	D_fav = 0.689*vertex.zhad**(-1.039)*(1.0-vertex.zhad)**1.241
c	D_unfav = 0.217*vertex.zhad**(-1.805)*(1.0-vertex.zhad)**2.037

C

C Paramaterization from Binneweis et al (PRD 52 p.4947)

C  sv is their scaling variable = ln( ln(Q2/Lambda^2)/ln(Q2_0/Lambda^2) )
C  Q_0 sets the scale - for light quarks (u,d,s) they get Q_0 = sqrt(2) GeV
C  Lambda is 227 MeV in NLO.
c
c	sv = log( log(Q2gev/lambda**2)/log(Q2zero/lambda**2) )
c
C Form of parameterization is D = N z^a1 (1-z)^a2
c
c	if(doing_semipi) then
c
c	   N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
c	   a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
c	   a2 = 1.430 + 0.543*sv - 0.023*sv**2
c
c	   Ns = 4.250 - 3.147*sv + 0.755*sv**2
c	   a1s = -0.770 -0.573*sv + 0.117*sv**2
c	   a2s = 4.48 + 0.890*sv - 0.138*sv**2
c
C This is Du (pi+ + pi-) = = Dd (pi+ + pi-) = D_favored + D_unfavored
c	   D_sum = N*zhad**a1*(1.0-zhad)**a2
C This is Ds (pi+ + pi-)
c	   D_sum_s = Ns*zhad**a1s*(1.0-zhad)**a2s
c
c
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
c
c	   R_D = (1.0-zhad)**0.083583/(1.0+zhad)**1.9838
c	   
c	   D_fav = D_sum/(1.0+R_D)
c	   D_unfav = D_sum/(1.0+1.0/R_D)
c
C Assume Ds(pi+) = Ds(pi-) = Dsbar(pi+) = Dsbar(pi-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
c
c	   D_s = D_sum_s/2.0
c
c	   if(sum_sq.gt.0.) then
c	      if(doing_hplus) then
c		 dsigdz = (qu**2*u*D_fav + qu**2*ubar*D_unfav +
c     >  	      qd**2*d*D_unfav + qd**2*dbar*D_fav + 
c     >  	      qs**2*s*D_s + qs**2*sbar*D_s)/sum_sq
c	      else
c		 dsigdz = (qu**2*u*D_unfav + qu**2*ubar*D_fav +
c     >  	      qd**2*d*D_fav + qd**2*dbar*D_unfav + 
c     >  	      qs**2*s*D_s + qs**2*sbar*D_s)/sum_sq
c	      endif
c	      if(doing_deutsemi) then
c		 if(doing_hplus) then
c		    dsigdz = dsigdz + (qu**2*d*D_fav + qu**2*dbar*D_unfav +
c     >  		 qd**2*u*D_unfav + qd**2*ubar*D_fav + 
c     >  		 qs**2*s*D_s + qs**2*sbar*D_s)/sum_sq
c		 else
c		    dsigdz = dsigdz + (qu**2*d*D_unfav + qu**2*dbar*D_fav +
c     >  		 qd**2*u*D_fav + qd**2*ubar*D_unfav + 
c     >  		 qs**2*s*D_s + qs**2*sbar*D_s)/sum_sq
c		 endif
c	      endif
c	   else
c	      dsigdz = 0.0
c	   endif
c	elseif (doing_semika) then
c
c	   N = 0.310 - 0.038*sv - 0.042*sv**2
c	   a1 = -0.980 - 0.260*sv + 0.008*sv**2
c	   a2 = 0.970 + 0.978*sv - 0.229*sv**2
c
c	   Ns = 1.080 - 0.469*sv + 0.003*sv**2
c	   a1s = -0.820 -0.240*sv - 0.035*sv**2
c	   a2s = 2.550 + 1.026*sv - 0.246*sv**2
c
C This is Du (K+ + K-) = Ds(K+ + K-) = D_favored + D_unfavored (K+ + K-)
c	   D_sum = N*zhad**a1*(1.0-zhad)**a2
C This is Dd (K+ + K-) (for convenience, I still call it D_sum_s)
c	   D_sum_s = Ns*zhad**a1s*(1.0-zhad)**a2s
c
c
C Here I make the wild assumption that the ratio of unfavored to favored
C fragmentation functions for Kaons is the same as that for pions.
c
C       Ratio of D-/D+ from P. Geiger's thesis (HERMES)
c
c	   R_D = (1.0-zhad)**0.083583/(1.0+zhad)**1.9838
c	   
c	   D_fav = D_sum/(1.0+R_D)
c	   D_unfav = D_sum/(1.0+1.0/R_D)
c
C Assume Dd(K+) = Dd(K-) = Ddbar(K+) = Ddbar(K-)
C Note that this contrdicted by the HERMES data, but shouldn't make much
C difference for pions anyway.
c
c	   D_s = D_sum_s/2.0
c	
c	   if(sum_sq.gt.0.) then
c	      if(doing_hplus) then
c		 dsigdz = (qu**2*u*D_fav + qu**2*ubar*D_unfav +
c    >  	      qd**2*d*D_s + qd**2*dbar*D_s + 
c    >  	      qs**2*s*D_unfav + qs**2*sbar*D_fav)/sum_sq
c	      else
c		 dsigdz = (qu**2*u*D_unfav + qu**2*ubar*D_fav +
c     >  	      qd**2*d*D_s + qd**2*dbar*D_s + 
c     >  	      qs**2*s*D_fav + qs**2*sbar*D_unfav)/sum_sq
c	      endif
c	      if(doing_deutsemi) then
c		 if(doing_hplus) then
c		    dsigdz = dsigdz + (qu**2*d*D_fav + qu**2*dbar*D_unfav +
c     >  		 qd**2*u*D_s + qd**2*ubar*D_s + 
c     >  		 qs**2*s*D_unfav + qs**2*sbar*D_fav)/sum_sq
c		 else
c		    dsigdz = dsigdz + (qu**2*d*D_unfav + qu**2*dbar*D_fav +
c     >  		 qd**2*u*D_s + qd**2*ubar*D_s + 
c     >  		 qs**2*s*D_fav + qs**2*sbar*D_unfav)/sum_sq
c		 endif
c	      endif
c	   else
c	      dsigdz = 0.0
c	   endif
c	endif	

C Peter Bosted's parameterization
c uses a modified scaling variable (replaces z) for the fragmentation function
c (see Accardi et al https://arxiv.org/abs/0907.2395)
c new PB fit using zp for pions. This is z * D
	if (doing_semipi) then
	   xp = 2.*xbj / (1. + sqrt(1. + 4. * xbj**2 * Mp_gev**2 / q2gev))
	   zp = (zhad / 2.) * (xp / xbj) *(1. + 
     >           sqrt(1 - 4 * xbj**2 * Mp_gev**2 *  
     >           (Mpi_gev**2 + pt2gev) / zhad**2 / q2gev**2))
	   sv = log(q2gev/2.)
	   yf = pf(1) * zp**(pf(2) + pf(4)*sv + pf(9)/w) * 
     >          (1.-zp)**(pf(3) + pf(5)*sv + pf(10)/w) 
	   yf = yf * (1. + pf(6)*zp + pf(7)*zp**2 + pf(8)*zp**3) *
     >               (1. + pf(11)/w + pf(12)/w**2)
	   yu = pu(1) * zp**(pu(2) + pu(4)*sv + pu(9)/w) * 
     >          (1.-zp)**(pu(3) + pu(5)*sv + pu(10)/w) 
	   yu = yu * (1. + pu(6)*zp + pu(7)*zp**2 + pu(8)*zp**3) *
     >               (1. + pu(11)/w + pu(12)/w**2)

	   if(doing_hplus) then
	      u1 = yf
	      d1 = yu
	   else
	      u1 = yu
	      d1 = yf
	   endif
	   ub = d1
	   db = u1
	   s1 = yu
	   sb = s1
 	elseif(doing_semika) then
	   IHdss=2 ! Kaons, 1=Pions
	   IOdss=1 ! NLO
	   if(doing_hplus) then
	      ICdss = 1		! charge +
	   else
	      ICdss = -1	! charge -
	   endif
	   call fDSS (IHdss,ICdss,IOdss, Zhad, Q2gev, 
     >       U1, UB, D1, DB, S1, SB, C1, B1, GL1)
c	   write(6,*) 'cheesy poofs', zhad, u1, ub, d1, db, s1, sb
	endif

	dsigdz = (qu**2 * u    * u1 + 
     >            qu**2 * ubar * ub +
     >  	  qd**2 * d    * d1 + 
     >            qd**2 * dbar * db + 
     >  	  qs**2 * s    * s1 + 
     >            qs**2 * sbar * sb)/sum_sq/zhad
	dsigdzp = dsigdz
	if(doing_deutsemi) then
	   dsigdzn =(qu**2 * d    * u1 + 
     >               qu**2 * dbar * ub +
     >  	     qd**2 * u    * d1 + 
     >               qd**2 * ubar * db + 
     >  	     qs**2 * s    * s1 + 
     >               qs**2 * sbar * sb)/sum_sq/zhad
	   dsigdz = dsigdzp + dsigdzn
	endif

	b =  1./ (0.120 * zhad**2 + 0.200)
	sighad = dsigdz*b*exp(-b*pt2gev)/2./pi
c	write(6,*) 'bad kitty', sighad

C DJG: OK - I THINK I've got it right now.
C DJG: Should be dsig/dOmega dE

c	sige = 2.*alpha**2*(y**2/2. + 1. - y  - mtar*xbj*y/2./Eb)*
c	1    (Eprime/1000.)/(Q2gev*y*(mtar/1000.)*(nu/1000.))
c	2    * sum_sq

c	F2 = xbj*sum_sq
c	F1 = F2/2./xbj
c

	wsq = Mp_gev**2 + q2gev * (1./xbj -1.)
	z8=1.
	a8=1.
	if(doing_deutsemi) a8=2.
c        call F1F2IN09(z8, a8, Q2gev, Wsq, F1, F2) 
c Updated to F1F2IN21 (Eric Christy's 2021 fit)
	call F1F2IN21(z8, a8, Q2gev, Wsq, F1, F2)

	W1 = F1/(mtar/1000.)
	W2 = F2/(nu/1000.)

	sin2th2 = Q2/4./Eb/Eprime
	cos2th2 = 1.-sin2th2
c	if(do_fermi) then
c	   W2coeff = (efer-abs(pfer)*pferz)*(efer-abs(pfer)*(pferx*Epx+pfery*Epy+pferz*Epz)/Eprime)/mtar**2 - sin2th2
c	else
	   W2coeff = cos2th2
c	endif

	sige = 4.*alpha**2*(Eprime/1000)**2/Q2gev**2 * ( W2*W2coeff + 2.*W1*sin2th2)

C DJG This dsig/(dOmega_e dE_e dz dpt**2 dPhi_had) in microbarn/GeV**3/sr 
	sigsemi = sige*sighad*(hbarc/1000.)**2*10000.0


C Need to convert to dsig/ (dOmega_e dE_e dE_h dCos(theta) dPhi_had
C This is just given by 1/omega * 2*p_h**2*cos(theta)

	jacobian = 1./(nu/1000.)*2.*(phad/1000.)**2*cthpq

C The 1.e6 converts from microbarn/GeV^2 to microbarn/MeV^2
	sigma_eepiX = sigsemi*jacobian/1.e6


* Note that there is an additional factor 'fac' included with the fermi-smeared cross
* section.   This takes into account pieces in the flux factor that are neglected (=1) in
* colinear collisions.  The flux factor is |v_1-v_2| * 2E_1 * 2E_2.
* For a stationary target, v_2=0 and so velocity term is v_1=1 (electron
* beam), and E_2=M_2.  For collinear boost, the flux factor can be expressed
* in a way that is lorenz invariant, and so can be used for lab or C.M.
* For a NON-COLLINEAR boost, there are two changes.  First, the |v| term
* becomes 1 - (z component of pfer)/efer.  Second, E_2 isn't just the mass,
* it becomes E_fermi, so we have to remove targ.Mtar_struck (which is used
* for E_2 by default) and replace it with efer.  Since the flux factor 
* comes in the denominator, we replace the usual flux factor (gtpr) with
* gtpr*fac, where fac = 1/ ( (1-pfer_z/efer)* (efer/mtar_struck) ).

	if(do_fermi) then
	   fac = 1./(1.-pferz*pfer/efer) * mtar/efer
	else
	   fac = 1.0
	endif

	sigma_eepiX = sigma_eepiX*fac
	
	if(doing_cent) then
	   main%johnjac = jacobian*1000.0
	else
	   main%davejac = jacobian*1000.0
	   ntup%sigcm = sighad
	endif



	peepiX = sigma_eepiX


c If doing_decay=.false., generate survival probability for main.weight.
c main.FP.p.path is dist. to back of detector, want decay prob. at aerogel.
C NOTE THAT ZAERO IS TAKEN WITH RESPECT TO THE POSITION AT WHICH PATHLEN
C IS CALCULATED (USUALLY THE BACK OF THE CALORIMETER).  IF THE DRIFTS IN
C MC_SOS_HUT OR MC_HMS_HUT ARE CHANGED, THEN THE STARTING POINT MAY BE DIFFERENT AND THESE
C NUMBERS MAY BE WRONG.  AERO BETWEEN S2Y AND S2X IN SOS AND BETWEEN DC2 S1X in HMS

C Beta/Gamma for decay need to use momentum after radiation/eloss, not vertex
C momentum.  Get particle momentum from main.SP.p.delta

	if (.not.doing_decay) then
	  if (hadron_arm.eq.1) then
	    zaero = -331.491		!aero at 40.199 cm, last project=371.69 cm
	  else if (hadron_arm.eq.2) then
*	    zaero = -76.		!aero at 270cm,last project=346(cal).
	    zaero = -82.8		!From Rick: aero at 263.2cm,last project=346(cal).
	  else if (hadron_arm.eq.3) then
	    zaero = -183.		!aero at 130cm,last project=313(S2)
	  else if (hadron_arm.eq.4) then
	    zaero = -183.
	  endif
	  pathlen = main%FP%p%path + zaero*(1+main%FP%p%dx**2+main%FP%p%dy**2)
	  p_kaon = spec%p%P*(1.+main%SP%p%delta/100.)
	  betak = spec%p%P/sqrt(spec%p%P**2+Mh2)
	  gammak = 1./sqrt(1.-betak**2)
	  survivalprob = 1./exp(pathlen/(ctau*betak*gammak))
	  decdist = survivalprob		!decdist in ntuple
	endif

	return

	end
