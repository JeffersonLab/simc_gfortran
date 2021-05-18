	real*8 function peepi(vertex,main)

* Purpose:
* This function determines the kinematics in the PHOTON-NUCLEON center of mass
* frame, calculates some of the kinematical variables (s,t, and CM quantities
* in the 'main' structure), and returns the pion cross section.
*
*   output:
*	peepi		!d5sigma/dEe'dOmegae'Omegapi	(microbarn/MeV/sr^2)

	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

	real*8 sigma_eepi
	real*8 k_eq			!equivalent photon energy.
	real*8 gtpr			!gamma_t prime.
	real*8 fac
	real*8 tfcos,tfsin		!cos/sin of theta between pfermi and q
	real*8 s

! Variables calculated in transformation to gamma-NUCLEON center of mass.
        real*8 gstar,bstar,bstarx,bstary,bstarz		!beta of boost to C.M.
        real*8 nustar,qstar,qstarx,qstary,qstarz	!q in C.M.
        real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz		!p_hadron in C.M.
        real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
        real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz	!p_fermi in C.M.
        real*8 thetacm,phicm,phiqn,jacobian,jac_old
	real*8 Wgev, Q2gev, E0, cthcm, sig0, fac1

	real*8 sig_multipole,sig_blok,sig_param04,sig_param_3000,sig_param_2021

	integer final_state, ipi
	logical first,low_w_flag

	data first /.TRUE./
cdg	data low_w_flag /.FALSE./	!Assume high W kinematics to start

* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.

	call transform_to_cm(vertex,main,
     &		gstar,bstar,bstarx,bstary,bstarz,
     &		nustar,qstar,qstarx,qstary,qstarz,
     &		epicm,ppicm,ppicmx,ppicmy,ppicmz,
     &		ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz,
     &		etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz,
     &		thetacm,phicm,phiqn,jacobian,jac_old)

	main%thetacm = thetacm
	main%phicm = phicm
	main%pcm = ppicm
	main%davejac = jacobian
	main%johnjac = jac_old		!approx. assuming collinear boost.

!	write (6,*) jacobian,jac_old,100.*(jacobian-jac_old)/jacobian,'%'


* calculate some kinematical variables
* 'f' and 'fer' indicate fermi momenta. 'star' or 'cm' indicate CM system
* Some of the physics calculations (t,epsi,s, etc...) are redundant with
* the calculations in event.f.  We should use the main.* variables from
* complete_ev where possible.  WORSE YET, WE CHANGE UNITS OF MAIN.W,... HERE!!!

	tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0
	tfsin=sqrt(1.-tfcos**2)

	s = (vertex%nu+efer)**2-(vertex%q+pfer*tfcos)**2-(pfer*tfsin)**2
	main%wcm = sqrt(s)



*******************************************************************************
* Get photon flux factor (two options, see comments below).
*
* DJG,2000: Replace targ.Mtar_struck in denominator of gammaflux with more 
* general efer-pfer*tfcos, for pfer =0 this reverts to old form
*	k_eq = (s-targ%Mtar_struck**2)/2./(efer-pfer*tfcos)
*
* JRA,2001: Go back to original version - more consistent with phase space used
* in the subroutine (according to DJG - see gaskell_model.ps)
	k_eq = (main%wcm**2-targ%Mtar_struck**2)/2./targ%Mtar_struck

CDG - Remove multipole model
CDG* Initialize some stuff for the multipole model.
CDG	if (first) then
CDG	  first = .false.
CDG	  if(which_pion.eq.0 .or. which_pion.eq.10) then  !pi+
CDG	    final_state = 1
CDG	  else
CDG	    final_state = 2
CDG	  endif
CDG	  if(k_eq.lt.500.)then
CDG	    write(6,*) 'Using low W multipole pion model...'
CDG	    low_w_flag = .TRUE.
CDG	  endif
CDG	endif

* Get cross section in photon-nucleon center of mass.  sigcm1 is blok
* fit (default model - can always be evaluated).  sigcm2 is multipole,
* IF low_w_flag is set.
* NOTE: s, t, mtar, and Q2 must be converted to GeV first.

c	ntup.sigcm1 = sig_blok(thetacm,phicm,main%t/1.e6,vertex%q2/1.e6,s/1%e6,main.epsilon,
c     >		targ%Mtar_struck/1000.,which_pion)

CDG Change default to PARAM04 - this works better at larger Q2
c	ntup%sigcm1 = sig_param04(thetacm,phicm,main%t/1.e6,vertex%q2/1.e6,s/1.e6,main%epsilon,
c     >		targ%Mtar_struck/1000.,which_pion)

CDG Change default to PARAM3000 - this works better at larger Q2
c	ntup%sigcm1 = sig_param_3000(thetacm,phicm,main%t/1.e6,vertex%q2/1.e6,s/1.e6,main%epsilon,
c     >		targ%Mtar_struck/1000.,which_pion)

CDG Use Peter Bosted's new fit to world, JLab 6 GeV, and preliminary 12 GeV data 
	ntup%sigcm1 = sig_param_2021(thetacm,phicm,main%t/1.e6,vertex%q2/1.e6,s/1.e6,main%epsilon,
     >           which_pion)

	sigma_eepi = ntup%sigcm1

CDG For lower W, user Peter Bosted's implementation of the MAID model
c	if(main%wcm.lt.2300) then ! W less than 2.3 GeV
	if(main%wcm.lt.2000) then ! W less than 2.0 GeV
	   Q2gev = vertex%q2/1.d6 !convert to GeV**2
	   Wgev = main%wcm/1000.0 ! convert to GeV
	   cthcm = cos(thetacm)
	   E0 = vertex%Ein/1000.0  !convert to GeV
	   ipi = 3 ! pi+ by default
	   if (which_pion.eq.1 .or. which_pion.eq.11 .or. which_pion.eq.3) ipi=4 ! pi-
	   call sigmaid(ipi,Q2gev,Wgev,E0,cthcm,phicm,sig0)
! convert from mub / dOmega* to mub / dt / dphi*
	   ntup%sigcm2 = sig0 / ppicm / qstar / 2.
cc       smoothly join between W of 1.7 and 2.3
c	   fac1 = max(0., (Wgev - 1.7) / 0.6)
c       smoothly join between W of 1.5 and 1.9
	   fac1 = min(1., max(0.,(Wgev - 1.5) / 0.4))
	   sigma_eepi = ntup%sigcm1 * fac1 + ntup%sigcm2 * (1 - fac1)
	endif

* For low w, use multipole expansion as default cross section model.
CDG	if(low_w_flag) then
CDG
CDG	  ntup%sigcm2 = sig_multipole(k_eq,efer,qstar,tfcos,ppicm,s/1.e6,thetacm,
CDG     >		phicm,main%epsilon,final_state,vertex%q2/1.e6,targ%Mtar_struck/1000.,pfer,mh)
CDG	  sigma_eepi = ntup%sigcm2
CDG
CDG	endif
	ntup%sigcm = sigma_eepi		!sig_cm


*******************************************************************************

* sigma_eepi is two-fold C.M. cross section: d2sigma/dt/dphi_cm [ub/MeV**2/rad]
* Convert from dt dphi_cm --> dOmega_lab using 'jacobian' [ub/sr]
* Convert to 5-fold by multiplying by flux factor, gtpr [1/MeV]
* to give d5sigma/dOmega_pi/dOmega_e/dE_e [ub/Mev/sr].
*
* Note that there is an additional factor 'fac' included with gtpr.   This
* takes into account pieces in the flux factor that are neglected (=1) in
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


	fac = 1./(1.-pferz*pfer/efer) * targ%Mtar_struck/efer
	gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*k_eq/vertex%q2/(1.-main%epsilon)

	peepi = sigma_eepi*jacobian*(gtpr*fac) !ub/MeV^2/rad-->ub/sr-->ub/MeV/sr

	return
	end


	real*8 function sig_multipole(nu,efer,qstar,tfcos,ppicm,s_gev,thetacm,
     >		phicm,eps,final_state,q2_gev,mtar_gev,pfer,mh)

* Purpose:
* This routine calculates p(e,e'pi+)n cross sections from a multipole fit.
* The multipole expansion gives dsigma/dOmega_cm.  We use dt/dcostheta_cm
* to give dsigma/dt/dphi_cm, which is returned as sig_multipole [ub/MeV^2-rad].

	implicit none

	complex*8 H1,H2,H3,H4,H5,H6         !Helicity amplitudes
	complex*8 amplitude(7),A1,A2,A3,A4,A12,A34,A1234
	complex*8 AT(7,3,8,10)
	complex*8 gf1,gf2,gf3,gf4,gf7,gf8
	complex*8 hnperp,hfperp,hnpar,hfpar,hnlong,hflong
	complex*8 XL,XT,XTT,XLT

	real*8 nu,efer,qstar,tfcos,ppicm,s_gev,thetacm,phicm,eps,q2_gev,mtar_gev,pfer,mh

	real*8 sig219,sig
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
	real*8 sintcm,costcm,sinto2cm,costo2cm
	real*8 root_two,mpi_to_microbarn
	real*8 CLEB(2,3),qfm,Q2_table(8),nu_table(10),atemp(42)
	real*8 Q2_high,Q2_low,deltaQ2,nu_high,nu_low,delta_nu,F

	integer Q2_count,nu_count,multipole,IT1,IT2,IT3,IT4,IT5
	integer final_state,isospin

	logical first

	data first /.TRUE./

	if(nu.ge.500.0) then
	  WRITE(6,*) 'LOOK OUT CHEESY POOF BREATH!'
	  WRITE(6,*) 'W IS TOO HIGH FOR THIS MODEL'
	  WRITE(6,*) 'SETTING NU to 499.9'
	nu=499.9
	endif

* Initialize some shorthand variables.
	sintcm = sin(thetacm)
	costcm = cos(thetacm)
	sinto2cm = sin(thetacm/2.)
	costo2cm = cos(thetacm/2.)


* DJG Read in multipoles and initialize CG coefficients for low W model.
*       Clebsch-Gordon Coefficients  CLEB(FINAL STATE, ISOSPIN)
*					  1: PI+ N--------1(IS.1/2)
*					  2: PI- P--------2(IV.1/2)
*						  --------3(IV.3/2)
	if(first) then
	  first = .FALSE.
	  root_two = sqrt(2.)
	  mpi_to_microbarn = 197.327/mh
	  CLEB(1,1)=root_two
	  CLEB(1,2)=root_two/3.
	  CLEB(1,3)=-root_two/3.
	  CLEB(2,1)=root_two
	  CLEB(2,2)=-root_two/3.
	  CLEB(2,3)=root_two/3.

	  OPEN(3,FILE='MULTIPOL.file',STATUS='OLD')
	  do Q2_count = 1,8
	    read(3,202) qfm				!q**2 in fm**-2
	    Q2_table(Q2_count) = qfm*0.0389379	!convert to GeV**2
	    do nu_count=1,10
	      read(3,203)nu_table(nu_count)
	      read(3,204)atemp	!multipolarity amplitudes, see,eg,table below
	      do multipole=1,7
		IT1=multipole+7
		IT2=multipole+14
		IT3=multipole+21
		IT4=multipole+28
		IT5=multipole+35

C      Multipoles in units of  10**-2 MPI+**-1
C  Multipole Amplitude Table: AT(TYPE    ,  ISOSPIN ,-Q**2 , P.EQ.PH.R.)
C				  1:-----E0+-------0(IS.1/2)
C				  2:-----M1--------1(IV.1/2)
C				  3:-----E1+-------3(IV.3/2)
C				  4:-----M1+
C				  5:-----S0+
C				  6:-----S1-
C				  7:-----S1+
		AT(multipole,3,Q2_count,nu_count)=
     >			CMPLX(ATEMP(multipole),ATEMP(IT1))
		AT(multipole,2,Q2_count,nu_count)=CMPLX(ATEMP(IT2),ATEMP(IT3))
		AT(multipole,1,Q2_count,nu_count)=CMPLX(ATEMP(IT4),ATEMP(IT5))
	      enddo	!multipole
	    enddo	!nu
	  enddo	!Q2
	  close(unit=3)
	endif		!if first time.

* Interpolate multipole amplitudes from table.
	do Q2_count=1,7
	  Q2_high=0.
	  Q2_low=0.
	  nu_high=0.
	  nu_low=0.
	  if((q2_gev.gt.Q2_table(Q2_count)).and.
     >	       (q2_gev.lt.Q2_table(Q2_count+1))) then
	    Q2_high = Q2_table(Q2_count+1)
	    Q2_low = Q2_table(Q2_count)
	    deltaQ2 = (Q2_high-Q2_low)
	    do nu_count=1,9
	      if((nu.gt.nu_table(nu_count)).and.
     >		   (nu.lt.nu_table(nu_count+1))) then
		nu_high = nu_table(nu_count+1)
		nu_low =  nu_table(nu_count)
		delta_nu = nu_high-nu_low
		do multipole=1,7
		amplitude(multipole) = 0.
		  do isospin = 1,3
		    A1 = AT(multipole,isospin,Q2_count,nu_count)
		    A2 = AT(multipole,isospin,Q2_count+1,nu_count)
		    A3 = AT(multipole,isospin,Q2_count,nu_count+1)
		    A4 = AT(multipole,isospin,Q2_count+1,nu_count+1)
		    A12= (A1*(Q2_high-q2_gev) + A2*(q2_gev-Q2_low))/deltaQ2
		    A34= (A3*(Q2_high-q2_gev) + A4*(q2_gev-Q2_low))/deltaQ2
		    A1234 = (A12*(nu_high-nu)+A34*(nu-nu_low))/delta_nu
		    amplitude(multipole) = amplitude(multipole) +
     >					CLEB(final_state,isospin)*A1234
		  enddo
*					!Convert multipoles to microb^.5
		amplitude(multipole) = mpi_to_microbarn*amplitude(multipole)
		enddo
	      endif	!nu check
	    enddo	!nu
	  endif		!Q2 check
	enddo		!Q2

* Now calculate the Helicity amplitudes

* CGLN amplitudes (Chew et. al.,Phys. Rev. 106,1345 (1957).)

	gf1 = amplitude(1) + 3.*costcm*(amplitude(4)+amplitude(3))
	gf2 = 2.*amplitude(4) + amplitude(2)
	gf3 = 3.*(amplitude(3)-amplitude(4))
	gf4 = (0.,0.)
* DJG: gf5 and gf6 are not independent of gf7 and gf8,so only use last two.
	gf7 = amplitude(6) - 2.*amplitude(7)
	gf8 = amplitude(5) + 6.*costcm*amplitude(7)

* DJG WALKER helicity amplitudes (Walker,Phys.Rev. 182,1729 (1969).)
	H5 = -costo2cm*(gf7+gf8)
	H6 = -sinto2cm*(gf7-gf8)
	H4 = (2.*sinto2cm*(gf1+gf2)+costo2cm*sintcm*(gf3+gf4))/root_two
	H2 = (-2.*costo2cm*(gf1-gf2)+sinto2cm*sintcm*(gf3-gf4))/root_two
	H1 = (-costo2cm*sintcm*(gf3+gf4))/root_two
	H3 = (-sinto2cm*sintcm*(gf3-gf4))/root_two

* DJG Bartl helicity amplitudes (Bartl,Nuc.Phys.B62,267 (1973).)
* key:		hnperp -----> photon pol. perp. to scatt. plane,no baryon flip
*		hfpar -----> photon pol. par. to scatt. plane, baryon flip
*		hnlong -----> longitudinal photon pol., no baryon flip
*		hflong -----> longitudinal photon pol., baryon flip
*		hnpar -----> photon pol. par. to scatt. plane, no baryon flip
*		hfperp -----> photon pol. perp. to scatt. plane, baryon flip

	hnperp = (H4+H1)/root_two
	hfpar = (H3+H2)/root_two
	hnlong = H5
	hflong = H6
	hnpar = (H4-H1)/root_two
	hfperp = (H3-H2)/root_two

	XT = hnperp*CONJG(hnperp)+hfpar*CONJG(hfpar)+hnpar*CONJG(hnpar)+
     >		hfperp*CONJG(hfperp)
	XL = hnlong*CONJG(hnlong) + hflong*CONJG(hflong)
	XTT = hnpar*CONJG(hnpar)+hfperp*CONJG(hfperp)-hnperp*CONJG(hnperp)-
     >		hfpar*CONJG(hfpar)
	XLT = hnlong*CONJG(hnpar)+hflong*CONJG(hfpar)

	F = 2.*((efer-pfer*tfcos)/1000.)*sqrt(s_gev)*(ppicm/1000.)/
     >		(s_gev-mtar_gev**2)/mtar_gev

	sigt = F/2.*XT
	sigl = F*XL
	sigtt = F/2.*XTT
	siglt = 2.*F*REAL(XLT)
	siglt=siglt/2.		!Divide siglt by 2 because Henk uses
				!a different convention when adding
				!the four terms together. DJG

	sig219=(sigt+eps*sigl+eps*cos(2.*phicm)*sigtt
     >		+sqrt(2.0*eps*(1.+eps))*cos(phicm)*siglt)/1.e0

* DJG: sig219 is dsig/dOmega_cm - convert to dsig/dtdphi_cm
* DJG: using dt/dcos_cm = 2*ppicm*qcm

	sig = sig219/ppicm/qstar/2.
	sig_multipole = sig

202	format(/11X,f5.1/)
203	format(11X,f5.0)
204	format(6(/9X,7f8.3))

	return
	end





	real*8 function sig_blok(thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev,which_pion)

* Purpose:
* This routine calculates p(e,e'pi+)n cross sections from a fit to the data of
* Brauel et al., Z.Phys.C. 3(1979)101.
* Fit gives dsigma/dt/dphi_cm, which is returned as sig_blok [ub/MeV^2-rad].

	implicit none
	include 'constants.inc'

	real*8 sig219,sig
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
	
	real*8 mrho_emp
	real*8 fpi,fpi2

	real*8 thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev
	integer which_pion

* Fits to p(e,e'pi+)n cross-sections
* HPB: cross sections for pi-fits to Brauel's n(e,e'pi-)p
* HPB: dsigma/dt cross sections at W=2.19 and Q^2=0.70 MODIFIED BY HAND!!!

	  sigl = 27.8*exp(-11.5*abs(t))
	  sigt = 10.0*(5.*abs(t))*exp(-5.*abs(t))
	  siglt= 0.0*sin(thetacm)
	  sigtt= -(4.0*sigl+0.5*sigt)*(sin(thetacm))**2

	  if (which_pion.eq.1 .or. which_pion.eq.11 .or. which_pion.eq.3) then	!pi-
	    sigt =sigt *0.25*(1.+3.*exp(-10.*abs(t)))
	    sigtt=sigtt*0.25*(1.+3.*exp(-10.*abs(t)))
	  endif

* GMH: Now scale sigl by Q2*pion_formfactor
* HPB: parametrization of pion formfactor in the region 0.5<Q2<1.8

	  mrho_emp = 0.712   !DETERMINED BY EMPIRICAL FIT TO FORMFACTOR (GeV/c**2)

	  fpi=1./(1.+1.65*q2_gev+0.5*q2_gev**2)
	  fpi2=fpi**2

* HPB: now convert to different Q^2
* HPB: s_L follows Q^2*Fpi^2; s_T and s_TT go as 1/(0.3+Q^2)?
* HPB: factor 0.1215 therefore is value of q2*fpi**2 for q2=0.7

	  sigl=sigl*(fpi2*q2_gev)/0.1215
	  sigt=sigt/(0.3+q2_gev)
	  sigtt=sigtt/(0.3+q2_gev)

	  sig219=(sigt+eps*sigl+eps*cos(2.*phicm)*sigtt
     >		+sqrt(2.0*eps*(1.+eps))*cos(phicm)*siglt)/1.e0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19

	  sig=sig219*15.333/(s_gev-mtar_gev**2)**2
	  sig=sig/2./pi/1.d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	  sig_blok = sig

	  return
	end

	real*8 function sig_param04(thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev,which_pion)

* Purpose:
* Fit that reproduces Fpi1 results pretty well, but is nicely behaved at 
* larger Q2
* Fit gives dsigma/dt/dphi_cm, which is returned as sig_blok [ub/MeV^2-rad].

	implicit none
	include 'constants.inc'

	real*8 sig219,sig
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
	
	real*8 fpi,q2fpi2,q2fpi,polefactor

	real*8 thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev
	integer which_pion

	fpi=1./(1.+1.77*q2_gev+0.05*q2_gev**2)
	q2fpi=q2_gev*fpi
	q2fpi2=q2_gev*fpi**2
	polefactor=abs(t)/((abs(t)-0.02)**2)
	
	sigL =350.0*q2fpi2*exp((16.0-7.5*alog(q2_gev))*(-t))
	sigT =4.5/q2_gev+2.0/q2_gev**2
	sigLT=(exp(0.79-3.4/sqrt(q2_gev)*(t))
     >       +1.1-3.6/q2_gev**2)*sin(thetacm)
	sigTT=-5.0/q2_gev**2*polefactor*sin(thetacm)**2


CDG For now assume sigL(pi+)=sigL(pi-)
	if (which_pion.eq.1 .or. which_pion.eq.11 .or. which_pion.eq.3) then	!pi-
	   sigt =sigt *0.25*(1.+3.*exp(-10.*abs(t)))
	   sigtt=sigtt*0.25*(1.+3.*exp(-10.*abs(t)))
	endif


	  sig219=(sigt+eps*sigl+eps*cos(2.*phicm)*sigtt
     >		+sqrt(2.0*eps*(1.+eps))*cos(phicm)*siglt)/1.e0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19

	  sig=sig219*8.539/(s_gev-mtar_gev**2)**2
	  sig=sig/2./pi/1.d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	  sig_param04 = sig

	  return
	end

***********************************************************************************************
	real*8 function sig_param_3000(thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev,which_pion)

* Purpose:
* Fit that reproduces Fpi1 (larger Q2)and Fpi2 and Brauel separated xsec data.
* while reproducing Bebek unseparated data.
* The Fpi1 low Q2 data isn't fit terribly well, but not so bad (mostly, sigma_l is
* a little big).
* Fit gives dsigma/dt/dphi_cm, which is returned as sig_param_3000 [ub/MeV^2-rad].

	implicit none
	include 'constants.inc'

	real*8 sig219,sig
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
	
	real*8 mrho_emp
	real*8 fpi,fpi2,q2fpi2,q2fpi,polefactor

	real*8 thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev
	integer which_pion



	fpi = 1.0/(1.0+1.77*q2_gev+0.12*q2_gev**2)

	q2fpi2=q2_gev*fpi**2

	sigL = 19.8*abs(t)/(abs(t)+0.02)**2*q2fpi2*exp(-3.66*abs(t))
	sigT = 5.96/q2_gev*exp(-0.013*q2_gev**2)

	sigLT=(16.533/(1.0+q2_gev))*exp(-5.1437*abs(t))*sin(thetacm)
	sigTT=(-178.06/(1.0+q2_gev))*exp(-7.1381*abs(t))*sin(thetacm)**2



CDG For now assume sigL(pi+)=sigL(pi-)
	if (which_pion.eq.1 .or. which_pion.eq.11 .or. which_pion.eq.3) then	!pi-
	   sigt =sigt *0.25*(1.0+3.0*exp(-10.0*abs(t)))
	   sigtt=sigtt*0.25*(1.0+3.0*exp(-10.0*abs(t)))
	endif


	  sig219=(sigt+eps*sigl+eps*cos(2.0*phicm)*sigtt
     >		+sqrt(2.0*eps*(1.0+eps))*cos(phicm)*siglt)/1.0d0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19

	  sig=sig219*8.539/(s_gev-mtar_gev**2)**2
	  sig=sig/2.0/pi/1.0d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	  sig_param_3000 = sig

	  return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sigmaid(Ipi,q2,w,e0,costh,phi,sig0)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c get exclusive pion cross sections and target
c asymmetries from MAID 2007. 
c inputs: ipi = 1, 2, 3, 4 for pi0 p, pi0 n, pi+ n, pi- p
c         q2 in GeV2. For Q2>5, extrapolation is done
c         W  in GeV. For W>2, extrapolation is done
c         E0 beam energy in GeV
c         costh cos(th*) (from -1 to 1)
c         phi   phi* (in radians, from 0 to 2pi)
c outputs: cross section dsigma/dOmega in the N-pi
c         c.m. system, in mub/sr
c         sigma = sig0 + Pt * sigz + Pb * Pt * sigez
c         where Pb and Pt are beam and target polarization
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer i,Ipi,iq,iw,ith,j,k
      real*8 q2,w,e0,costh,phi,sig0,sigz, sigez, sige
      real*8 maidtbl(4,25,46,23,18)
      real*8 wfact, ST, SL, STL, STT, STLZ, STTZ, STLP
      real*8 STLPZ, STTPZ, CSF, SNF, CS2F, SN2F, CST,SNT
      real*8 STLPY,STTPX,SLY,STY,STLX,sthe,pt,pl,pn
      REAL*8 STLY,STTX,STTY,STLPX
      real*8 x, eps, am/0.9383/, nu, ep, sin2
      real*8 cthmin(6)/-0.20, 0.20, 0.44, 0.63, 0.78, 0.90/
      real*8 cthmax(6)/       0.20, 0.44, 0.63, 0.78, 0.90, 1.0/
      logical first/.true./


! Read in the MAID tables generated by 
! ~bosted/eg1c/MAID/targetall.f using code 
! provided by Lothar Tiator (MAID group)
      if(first) then
       first = .false.
       do i=1,4
        if(i.eq.1) open(unit=9,file='maidpi0p.dat')
        if(i.eq.2) open(unit=9,file='maidpi0n.dat')
        if(i.eq.3) open(unit=9,file='maidpipn.dat')
        if(i.eq.4) open(unit=9,file='maidpimp.dat')
        do iq=1,25
         do iw=1,46
          do ith=1,23
           read(9,'(f11.6,18f8.4)') 
     >      (maidtbl(i,iq,iw,ith,j),j=1,18)
          enddo
         enddo
        enddo
        close(unit=9)
       enddo
      endif

! Initialize answers
      sig0 = 0.
      sigz = 0.
      sigez = 0.
      st = 0.
      sige = 0.

      if(w.lt.1.08) return

! get epsilon
      nu = (w**2 - am**2 + q2) / 2. / am
      if(nu .gt. e0) return
      ep = e0 - nu
      sin2 = q2 / 4. / e0 / ep
      if(sin2.le.0.0 .or. sin2.gt. 1.) return
      eps = 1./(1. + 2. * (1. + nu**2 / q2) * sin2 / 
     >       (1. - sin2))
! tHIS IS sin(theta_electron)
      sthe = 2. * sqrt(sin2)
! this is sin(theta_q)
      pt = ep * sthe / sqrt(q2 + nu**2)
! this is cos(theta_q) (target b-field along q)
      pl = sqrt(1. - pt**2)
! This is normal component of target b-field
      pn  = -pt * sin(phi)
! This is now the in-plane component
      pt  =  pt * cos(phi)


! Get nearest bin in Q2, W, and costh
      iq = int((q2 + 0.1)/0.2)
      iq = min(25,max(1,iq))

      iw = int((w-1.090)/0.020)
      iw = min(46,max(1,iw))

      ith=0
      do i=1,6
       if(costh.ge.cthmin(i).and.costh.le.cthmax(i)) ith=i
      enddo
      ith = min(6,max(1,ith))


! Get sigma_t in mub/sr, as well as sigm_l, ...
      wfact = 1.
      if(w.gt.1.232) wfact = (w -1.132)/0.100
      ST = maidtbl(Ipi,iq,iw,ith,1) / max(0.2,q2) / wfact
      SL   = maidtbl(Ipi,iq,iw,ith,2) * ST
      STL  = maidtbl(Ipi,iq,iw,ith,3) * ST
      STT  = maidtbl(Ipi,iq,iw,ith,4) * ST
      STLZ = maidtbl(Ipi,iq,iw,ith,5) * ST
      STTZ = maidtbl(Ipi,iq,iw,ith,6) * ST
      STLPZ= maidtbl(Ipi,iq,iw,ith,7) * ST
      STTPZ= maidtbl(Ipi,iq,iw,ith,8) * ST
      STLP = maidtbl(Ipi,iq,iw,ith,9) * ST
      STLX = maidtbl(Ipi,iq,iw,ith,10) * ST
      STLY = maidtbl(Ipi,iq,iw,ith,11) * ST
      STTX = maidtbl(Ipi,iq,iw,ith,12) * ST
      STTY = maidtbl(Ipi,iq,iw,ith,13) * ST
      STLPX= maidtbl(Ipi,iq,iw,ith,14) * ST
      STLPY= maidtbl(Ipi,iq,iw,ith,15) * ST
      STTPX= maidtbl(Ipi,iq,iw,ith,16) * ST
      SLY  = maidtbl(Ipi,iq,iw,ith,17) * ST
      STY  = maidtbl(Ipi,iq,iw,ith,18) * ST

! get cos(phi), etc.
      CSF =COS(phi)
      SNF =SIN(phi)
      CS2F=COS(2.*phi)
      SN2F=SIN(2.*phi)
      CST =costh
      SNT =SQRT(1.-costh**2)

! Get sigma_0, etc.
      sig0 = ST + EPS * SL + 
     >    sqrt(2.*EPS*(1.+EPS))*CSF*STL + 
     >    EPS*CS2F*STT

      sigz = sqrt(2.*EPS*(1.+EPS)) *
     > (PT * SNF * STLX + 
     >  PN * CSF * STLY + 
     >  PL * SNF * STLZ) + 
     > EPS * 
     > (PT * SN2F * STTX +
     >  PN * CS2F * STTY +
     >  PL * SN2F * STTZ) +
     > PN * (STY + EPS * SLY)

      sigez = sqrt(2.*EPS*(1.-EPS)) *
     > (PT * STLPx * CSF + 
     >  PN * STLPy * SNF + 
     >  PL * STLPZ * CSF) + 
     > sqrt(1.-EPS**2) * 
     >  (PT * STTPx + 
     >   PL * STTPZ)

      sige = sqrt(2.*EPS*(1.-EPS))*SNF*STLP

! change sign of sigez to agree with DIS convention
! in which Delta has negative A_LL
      sigez = -1.0 * sigez

! also change sign of sigz so agrees with our results for
c W>1.4. Note: this was *not* done for eg1c!
      sigz = -1.0 * sigz

      return
      end

       real*8 function sig_param_2021(thcm,phicm,t,q2,wsq,eps,which_pion)
! April 2021 fit to exclusive pi+ and pi- data from
! fpi1, fpi2, CT, pt-SIDIS, CSV-SIDIS, and KLT
! q2, t, and wsq should be in gev**2
! thcm and phicm should be in radians
! Peter Bosted, from output of ~bosted/ptc/ptb.f
      implicit none
      real*8 thcm,phicm,t,q2,wsq,eps
      real*8 sigl,sigt,sigttv,sigltv,sigexcl3
      integer which_pion
      common/excmn/sigl,sigt,sigttv,sigltv
      real*8 pp(17)/
     >     1.60077,
     >    -0.01523,
     >    37.08142,
     >    -4.11060,
     >    23.26192,
     >     0.00983,
     >     0.87073,
     >    -5.77115,
     >  -271.08678,
     >     0.13766,
     >    -0.00855,
     >     0.27885,
     >    -1.13212,
     >    -1.50415,
     >    -6.34766,
     >     0.55769,
     >    -0.01709/
      real*8 pm(17)/
     >     1.75169,
     >     0.11144,
     >    47.35877,
     >    -4.69434,
     >     1.60552,
     >     0.00800,
     >     0.44194,
     >    -2.29188,
     >   -41.67194,
     >     0.69475,
     >     0.02527,
     >    -0.50178,
     >    -1.22825,
     >    -1.16878,
     >     5.75825,
     >    -1.00355,
     >     0.05055/

      if (which_pion.eq.1.or.which_pion.eq.11.or.which_pion.eq.3) then
       call exclfit(t,thcm,phicm,q2,
     >    wsq,eps,sigl,sigt,sigttv,sigltv,sigexcl3,pm)
      else ! pi+
       call exclfit(t,thcm,phicm,q2,
     >    wsq,eps,sigl,sigt,sigttv,sigltv,sigexcl3,pp)
      endif
      sig_param_2021 = sigexcl3
      return 
      end

      subroutine exclfit(t,thetacm,phicm,q2_gev,
     >  s_gev,eps,sigl,sigt,sigtt,siglt,sig,p)

      implicit none
      real*8 p(15),t,thetacm,phicm,q2_gev,s_gev,eps,sigl
      real*8  sigt,sigtt,siglt,sig,mtar_gev/0.938/
      real*8 fpi,q2fpi2,sig219

       fpi = 1. / 
     >  (1.0 + p(1)*q2_gev + p(2)*q2_gev**2)

       q2fpi2 = q2_gev * fpi**2

cxx       sigL = p(3) * abs(t) / 
       sigL = (p(3) + p(15)/q2_geV) * abs(t) / 
     >   (abs(t) + 0.02)**2 * q2fpi2 * 
     >   exp(p(4) * abs(t))
c       sigL = sigL / s_gev**p(11)
       sigL = sigL / (s_gev**p(11) + sqrt(s_gev)**p(17))

       sigT = p(5) / q2_gev * exp(p(6) * q2_gev**2)
c       sigT = sigT / s_gev**p(12)
       sigT = sigT / (s_gev**p(12) + sqrt(s_gev)**p(16))
       sigT = sigT * exp(p(14) * abs(t))
 
      sigLT=(p(7) / (1.0 + p(10) * q2_gev)) *
     >   exp(p(8) * abs(t)) * sin(thetacm)
       sigLT = sigLT / s_gev**p(13)

       sigTT=(p(9) / (1. + 1.0 * q2_gev)) * 
     >   exp(-7.0 * abs(t)) * sin(thetacm)**2


        sig219 = (sigt + 
     >   eps * sigl + 
     >   eps * cos(2.0*phicm) * sigtt + 
     >   sqrt(2.0 * eps * (1.0+eps)) * 
     >   cos(phicm) * siglt) / 1.0d0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19
c BUT, in param_3000, everything got scaled to W=1.96

        sig = sig219 *8.539 / (s_gev - mtar_gev**2)**2
        sig = sig / 2.0 / 3.1415928 / 1.0d+06  !dsig/dtdphicm in microbarns/MeV**2/rad

        return
        end
