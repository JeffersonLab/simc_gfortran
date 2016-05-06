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

	real*8 sig_multipole,sig_blok,sig_param04,sig_param_3000

	integer final_state
	logical first,low_w_flag

	data first /.TRUE./
	data low_w_flag /.FALSE./	!Assume high W kinematics to start

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

* Initialize some stuff for the multipole model.
	if (first) then
	  first = .false.
	  if(which_pion.eq.0 .or. which_pion.eq.10) then  !pi+
	    final_state = 1
	  else
	    final_state = 2
	  endif
	  if(k_eq.lt.500.)then
	    write(6,*) 'Using low W multipole pion model...'
	    low_w_flag = .TRUE.
	  endif
	endif

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
	ntup%sigcm1 = sig_param_3000(thetacm,phicm,main%t/1.e6,vertex%q2/1.e6,s/1.e6,main%epsilon,
     >		targ%Mtar_struck/1000.,which_pion)

	sigma_eepi = ntup%sigcm1

* For low w, use multipole expansion as default cross section model.
	if(low_w_flag) then

	  ntup%sigcm2 = sig_multipole(k_eq,efer,qstar,tfcos,ppicm,s/1.e6,thetacm,
     >		phicm,main%epsilon,final_state,vertex%q2/1.e6,targ%Mtar_struck/1000.,pfer,mh)
	  sigma_eepi = ntup%sigcm2

	endif
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

	  if (which_pion.eq.1 .or. which_pion.eq.11) then	!pi-
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
	if (which_pion.eq.1 .or. which_pion.eq.11) then	!pi-
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
	if (which_pion.eq.1 .or. which_pion.eq.11) then	!pi-
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
