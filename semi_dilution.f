	real*8 function semi_dilution(vertex,main)


* Sept. 2003 D. Gaskell
* Purpose:
* This routine calculates p(e,e'pi+)X semi-iinclusive cross sections from
* the CTEQ5M parton distribution functions and a simple parameterization
* of the favored and unfavored fragmentation functions.
*   output:
*	sigma_eepiX	!d3sigma/dEe'dOmegae'Omegapi	(microbarn/MeV/sr^2)
* 
* For now, just does PI+ !!!!!!!!!!!!!!!!!!
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
* October 14 - Now adapted to get the dilution factor.



	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* PDFs
	integer iset  !which set (1=cteq5m)
	integer ipart !particle u=1, ubar=-1, d=2, dbar=-2
	real*8 u,d,ubar,dbar
	real*8 qu,qd ! u and d quark charges
	real*8 D_fav, D_unfav, D_sum, R_D !favored,unfavored,sum,ratio of FFs
	real*8 lambda, Q2zero ! scales for FF param
	real*8 sv !scaling variable for FF param
	real*8 N,a1,a2 !parameters for FF param

	real*8 xbj,y,s, Q2gev, Qgev, pt2gev
	real*8 b  ! pt2 parameter for FFs

	real*8 sum_sq, dsigdz, sigsemi, jacobian, sigma_eepiX
	real*8 sighad, sige

	real*8 sum_sq_prot,sum_sq_neut,dsigdz_prot,dsigdz_neut
	real*8 sighad_prot,sighad_neut, sige_prot, sige_neut
	real*8 sigsemi_prot,sigsemi_neut,dilution

	real*8 N15,He4,Al,Ni,Cu

	real*8 Ctq5Pdf	
c	external Ctq5Pdf

	logical first,kfirst


	parameter (iset=1)
	parameter (qu=2./3.)
	parameter (qd=-1./3.)

	parameter (b=3.76)   !GeV^-2

	parameter (lambda=0.227) !0.227 GeV for NLO
	parameter (Q2zero=2.0)   !Gev^2 for u,d,s,g

	data first /.TRUE./
	data kfirst /.TRUE./

	if(first) then
	   call SetCtq5(iset)    ! initialize Cteq5 (we're using cteq5m)
	   first=.FALSE.
	endif

	s = (2.*vertex%Ein*mp + mp**2)/1.e6  !convert to GeV2
	xbj = vertex%Q2/2./mp/vertex%nu   
	if(xbj.gt.1.0) then
	   write(6,*) 'XBj is too large!', xbj
	   xbj=1.0
	endif
	y = vertex%nu/vertex%Ein

C DJG convert some stuff to GeV

	Q2gev = vertex%q2/1.e6
	Qgev = sqrt(Q2gev)
	pt2gev = vertex%pt2/1.e6

C Get the PDFs
	ipart=1
	u = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-1
	ubar = Ctq5pdf (ipart , xbj, Qgev)

	ipart=2
	d = Ctq5pdf (ipart , xbj, Qgev)

	ipart=-2
	dbar = Ctq5pdf (ipart , xbj, Qgev)

C Simple paramaterization from Kretzer et al (EPJC 22 p. 269)
C for Q2=2.5.

c	D_fav = 0.689*vertex.zhad**(-1.039)*(1.0-vertex.zhad)**1.241
c	D_unfav = 0.217*vertex.zhad**(-1.805)*(1.0-vertex.zhad)**2.037

C

C Paramaterization from Binneweis et al

	sv = log( log(Q2gev/lambda**2)/log(Q2zero/lambda**2) )

C Form of parameterization is D = N z^a1 (1-z)^a2

	if(doing_semipi) then

	   N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
	   a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
	   a2 = 1.430 + 0.543*sv - 0.023*sv**2

	elseif(doing_semika) then

	   if (kfirst) then
	      write(6,*) 'Kaon FF in semi_dilution not iomplemented yet -using pion values'
	      kfirst=.false.
	   endif
	   N = 1.150 - 1.522*sv + 1.378*sv**2 - 0.527*sv**3
	   a1 = -0.740 - 1.680*sv + 1.546*sv**2 - 0.596*sv**3
	   a2 = 1.430 + 0.543*sv - 0.023*sv**2

	endif

	D_sum = N*vertex%zhad**a1*(1.0-vertex%zhad)**a2

C Ratio of D-/D+ from P. Geiger's thesis (HERMES)

	R_D = (1.0-vertex%zhad)**0.083583/(1.0+vertex%zhad)**1.9838

	D_fav = D_sum/(1.0+R_D)
	D_unfav = D_sum/(1.0+1.0/R_D)
	
	sum_sq = qu**2*(u+ubar) + qd**2*(d+dbar)

	sum_sq_prot = sum_sq
	sum_sq_neut = qu**2*(d+dbar) + qd**2*(u+ubar)

	if (doing_deutsemi) then
	   sum_sq = sum_sq + qu**2*(d+dbar) + qd**2*(u+ubar)
	endif

	if(sum_sq.gt.0.) then
	   if(doing_hplus) then
	      dsigdz = (qu**2*u*D_fav + qu**2*ubar*D_unfav +
     >  	   qd**2*d*D_unfav + qd**2*dbar*D_fav)/sum_sq
	   else
	      dsigdz = (qu**2*u*D_unfav + qu**2*ubar*D_fav +
     >  	   qd**2*d*D_fav + qd**2*dbar*D_unfav)/sum_sq
	   endif
	   if(doing_deutsemi) then
	      if(doing_hplus) then
		 dsigdz = dsigdz + (qu**2*d*D_fav + qu**2*dbar*D_unfav +
     >  	      qd**2*u*D_unfav + qd**2*ubar*D_fav)/sum_sq
	      else
		 dsigdz = dsigdz + (qu**2*d*D_unfav + qu**2*dbar*D_fav +
     >  	      qd**2*u*D_fav + qd**2*ubar*D_unfav)/sum_sq
	      endif
	   endif
	else
	   dsigdz = 0.0
	endif

	if(sum_sq_prot.gt.0.) then
	   if(doing_hplus) then
	      dsigdz_prot = (qu**2*u*D_fav + qu**2*ubar*D_unfav +
     >  	   qd**2*d*D_unfav + qd**2*dbar*D_fav)/sum_sq_prot

	      dsigdz_neut = (qu**2*d*D_fav + qu**2*dbar*D_unfav +
     >  	   qd**2*u*D_unfav + qd**2*ubar*D_fav)/sum_sq_neut
	   else
	      dsigdz_prot = (qu**2*u*D_unfav + qu**2*ubar*D_fav +
     >  	   qd**2*d*D_fav + qd**2*dbar*D_unfav)/sum_sq_prot

	      dsigdz_neut = (qu**2*d*D_unfav + qu**2*dbar*D_fav +
     >  	   qd**2*u*D_fav + qd**2*ubar*D_unfav)/sum_sq_neut
	   endif
	else
	   dsigdz_prot = 0.0
	   dsigdz_neut = 0.0
	endif

	sighad = dsigdz*b*exp(-b*pt2gev)/2./pi

	sige = alpha**2*(1.+(1.-y)**2)*(vertex%e%E/1000.)/
     >      (s*xbj*y**2*(mp/1000.)*(vertex%nu/1000.))
     >      * sum_sq


C DJG This dsig/(dOmega_e dE_e dz dpt**2 dPhi_had) in microbarn/GeV**3/sr 
	sigsemi = sige*sighad*(hbarc/1000.)**2*10000.0

	sighad_prot = dsigdz_prot*b*exp(-b*pt2gev)/2./pi

	sige_prot = alpha**2*(1.+(1.-y)**2)*(vertex%e%E/1000.)/
     >     (s*xbj*y**2*(mp/1000.)*(vertex%nu/1000.))
     >     * sum_sq_prot


C DJG This dsig/(dOmega_e dE_e dz dpt**2 dPhi_had) in microbarn/GeV**3/sr 
	sigsemi_prot = sige_prot*sighad_prot*(hbarc/1000.)**2*10000.0


	sighad_neut = dsigdz_neut*b*exp(-b*pt2gev)/2./pi

	sige_neut = alpha**2*(1.+(1.-y)**2)*(vertex%e%E/1000.)/
     >      (s*xbj*y**2*(mp/1000.)*(vertex%nu/1000.))
     >      * sum_sq_neut


C DJG This dsig/(dOmega_e dE_e dz dpt**2 dPhi_had) in microbarn/GeV**3/sr 
	sigsemi_neut = sige_neut*sighad_neut*(hbarc/1000.)**2*10000.0



	N15 = 8.0*sigsemi_neut + 7.0*sigsemi_prot
	He4 = 2.0*sigsemi_neut + 2.0*sigsemi_prot
        Al =  14.0*sigsemi_neut + 13.0*sigsemi_prot
        Cu =  35.*sigsemi_neut + 29.0*sigsemi_prot
	Ni =  31.0*sigsemi_neut + 28.0*sigsemi_prot


	if(sigsemi_prot.gt.0.0.and.sigsemi_neut.gt.0.0) then
	   dilution = sigsemi/(sigsemi+0.3333*N15+0.3165*He4+
     >  	0.0144*Al + 0.0031*Cu + 0.0013*Ni)    
	else
	   dilution = 999
	endif

C Need to convert to dsig/ (dOmega_e dE_e dE_h dCos(theta) dPhi_had
C This is just given by 1/omega * 2*p_h**2*cos(theta)

c	jacobian = 1./(vertex.nu/1000.)*2.*(vertex.p.P/1000.)**2
c	1    *cos(vertex.theta_pq)


c	sigma_eepiX = sigsemi*jacobian/1.e6

c	main.davejac = jacobian
c	ntup.sigcm = sighad

c	peepiX = sigma_eepiX
	semi_dilution = dilution

	return

	end
