	real*8 function peeK(vertex,main,survivalprob)

*
* Purpose:
* This routine calculates N(e,e'K)Y cross sections.
*
* variables:
*
*   output:
*	peeK	!d5sigma/(dE_e'*dOmega_e'*Omega_K) - fm^2/MeV/sr^2 ???
*
	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

	real*8 S,phi
	real*8 sigma_eek,sigma_saghai
	real*8 k_eq,gtpr,fac
	real*8 tfcos,tfsin
	real*8 pathlen,zaero,betak,gammak,p_kaon
	real*8 survivalprob

! Variables calculated in transformation to gamma-NUCLEON center of mass.
	real*8 gstar,bstar,bstarx,bstary,bstarz		!beta of boost to C.M.
	real*8 nustar,qstar,qstarx,qstary,qstarz	!q in C.M.
	real*8 ekcm,pkcm,pkcmx,pkcmy,pkcmz		!p_hadron in C.M.
	real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
	real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz	!p_fermi in C.M.
	real*8 thetacm,phicm,phiqn,jacobian,jac_old

	real*8 sig_factorized

* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.

	call transform_to_cm(vertex,main,
     &		gstar,bstar,bstarx,bstary,bstarz,
     &		nustar,qstar,qstarx,qstary,qstarz,
     &		ekcm,pkcm,pkcmx,pkcmy,pkcmz,
     &		ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz,
     &		etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz,
     &		thetacm,phicm,phiqn,jacobian,jac_old)

* Jacobian goes from dt dphi_cm --> dOmega_lab
* Cross section is in terms of dOmega_cm (not dt dphi_cm), so need
* to include dOmega_cm --> dt dphi_cm jacobian (which is 1/2.*pk_cm*q_cm)
	jacobian = jacobian/(2.*pkcm*qstar)
	jac_old = jac_old/(2.*pkcm*qstar)

	main%thetacm = thetacm
	main%phicm = phicm
	main%pcm = pkcm
	main%davejac = jacobian
	main%johnjac = jac_old          !approx. assuming collinear boost.
!	write (6,*) jacobian,jac_old,100.*(jacobian-jac_old)/jacobian,'%'


* calculate some kinematical variables
* 'f' and 'fer' indicate fermi momenta. 'star' or 'cm' indicate CM system

	tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0 
	tfsin = sqrt(1.-tfcos**2)

	s = (vertex%nu+efer)**2-(vertex%q+pfer*tfcos)**2-(pfer*tfsin)**2
	main%wcm = sqrt(s) 


*******************************************************************************
* Get cross section in photon-nucleon center of mass.  sigcm1 is Saghai
* model (in eekeek.f). sigcm2 is factorized model (Doug Koltenuk).


*NOTE: wants odd sign/units for s,q2,etc...  Also, theta is calculated
* above, but phi is not, SO WE ARE ALWAYS CALCULATING FOR PHI=0!!

	if (targ%Mrec_struck.lt.1150.) then	!lambda production.
	  call eekeek(s/1.e6,vertex%q2,main%thetacm,main%theta_pq,phi,main%epsilon,sigma_saghai)
	else
	  call eekeeks(s/1.e6,vertex%q2,main%thetacm,main%theta_pq,phi,main%epsilon,sigma_saghai)
	endif
	ntup%sigcm1 = sigma_saghai


* Factorization model (CURRENTLY SET UP FOR HYDROGEN ONLY!!!)

	ntup%sigcm2 = sig_factorized(vertex%q2,main%wcm,main%t,
     1  pkcm,targ%Mrec_struck,main%epsilon)
	
* Choose the cross section model to use by default.
!	sigma_eek = ntup.sigcm1		!Saghai
	sigma_eek = ntup%sigcm2		!Factorized model

	ntup%sigcm = sigma_eek


* Virtual photon to electron beam flux conversion factor
* DJG,2000: Replace targ.Mtar_struck in denominator of gammaflux with more
* general efer-pfer*tfcos, for pfer =0 this reverts to old form
!       k_eq = (s-targ%Mtar_struck**2)/2./(efer-pfer*tfcos)

* JRA,2001: Go back to original version - more consistent with phase space used
* in the subroutine (according to DJG - see gaskell_model.ps)
        k_eq = (main%wcm**2-targ%Mtar_struck**2)/2./targ%Mtar_struck


* 'sig' is two-fold C.M. cross section: d2sigma/Omega_cm [ub/sr].
* Convert from dt dOmega_cm --> dOmega_lab using 'jacobian' [ub/sr]
* Convert to 5-fold by multiplying by flux factor, gtpr [1/MeV]
* to give d5sigma/dOmega_pi/dOmega_e/dE_e [ub/Mev/sr]
* Note that 'jacobian' above is the full jacobian from generate_cm,
* (dt dphi_cm --> dOmega_lab) times 2*pkcm*qcm (dOmega_cm --> dt dphi_cm).
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

	peeK = sigma_eek*jacobian*(gtpr*fac) !ub/MeV^2/rad-->ub/sr-->ub/MeV/sr


c If doing_decay=.false., generate survival probability for main.weight.
c main.FP.p.path is dist. to back of detector, want decay prob. at aerogel.
C NOTE THAT ZAERO IS TAKEN WITH RESPECT TO THE POSITION AT WHICH PATHLEN
C IS CALCULATED (USUALLY THE BACK OF THE CALORIMETER).  IF THE DRIFTS IN
C MC_SOS_HUT ARE CHANGED, THEN THE STARTING POINT MAY BE DIFFERENT AND THESE
C NUMBERS MAY BE WRONG.  AERO BETWEEN S2Y AND S2X IN SOS.

C Beta/Gamma for decay need to use momentum after radiation/eloss, not vertex
C momentum.  Get particle momentum from main%SP%p%delta

	if (.not.doing_decay) then
	  if (hadron_arm.eq.1) then
	    zaero = 0.			!no aerogel yet, use full length.
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




	real*8 function sig_factorized(q2,w,t,pk,mrec)

* Purpose:
* This routine calculates p(e,e'K+)Lambda cross sections using a
* factorized cross section model from Doug Koltenuk's thesis.
*
* Cross section (at theta_cm=0, i.e. t=tmin) is F(W)*F(Q^2).
* The t-dependance at fixed Q^2,W is F(t)~exp(-A*t),
* so the cross section is F(W)*F(Q^2)*(F(t)/F(tmin))=F(W)*F(Q^2)*F(t-tmin)
*
* The model requires W, Q^2, t, and pkcm.

	real*8 q2,w,t,pk,mrec			!all in Mev,MeV**2
	real*8 nu,q,nucm,qcm,tmin,pktest
	real*8 w2val,q2val,tval,pkval,tminval	!Vars. used in model (GeV)
	real*8 fact_w,fact_q,fact_t

	include 'constants.inc'

! Initialize some stuff.  Start with intermediate variables, all in MeV.

	nu = ( w**2 + q2 - mp**2 )/2./Mp
	q = sqrt(q2+nu**2)
	qcm = q*(Mp/W)
	nucm = sqrt(qcm**2-q2)
	tmin = -1.*(Mk2 - q2 - 2*nucm*sqrt(pk**2+mk2) + 2*qcm*pk )

! Check center of mass pk, since we can get it from w2
	pktest = sqrt ( (w**2+mk2-mrec**2)**2/4./w**2 - mk2 )
	if (abs((pktest-pk)/pk).gt.0.001) then
	  write(6,*) 'Kaon C.M. momentum passed to sig_factorized does not agree with'
	  write(6,*) 'the value calculated from W'
	  write(6,*) 'Passed,calculated=',pk,pktest
	endif

! Parameters used by the model, all converted to GeV.
	q2val = q2/1.e6
	w2val = w**2/1.e6
	pkval = pk/1000.
	tval = t/1.e6
	tminval = tmin/1.e6

	if (mrec.lt.1150.) then	!lambda production
	  fact_q = 1./(q2val+2.67)**2
	  fact_t = exp(-2.1*(tval-tminval))
	  if (w2val.ne.0) then			!=0 for central event!!!
	   fact_w = 0.959*4.1959*pkval/(sqrt(w2val)*(w2val-0.93827**2))
	   fact_w = fact_w + (0.18*1.72**2*0.10**2) /
     >		( (w2val-1.72**2)**2 + 1.72**2*0.10**2 )
	  endif
	else
	  fact_q = 1./(q2val+0.79)**2
	  fact_t = exp(-1.0*(tval-tminval))
	  if (w2val.ne.0) then			!=0 for central event!!!
	   fact_w = 0.959*4.1959*pkval/(sqrt(w2val)*(w2val-0.93827**2))

!	   fact_w = fact_w + (0.18*1.72**2*0.10**2) /
!     >		( (w2val-1.72**2)**2 + 1.72**2*0.10**2 )
	  endif
	endif

	sig_factorized = fact_q*fact_t*fact_w

	return
	end



	subroutine eekeek(ss,q22,angl,theta,phi,epsi,sigma_eep)

	implicit none

	include 'simulate.inc'

	real*8 ss,q22,angl,theta,phi,epsi
	real*8 w,skc2,skc,q0,q0c,qr,aflx,aflxl,an,x,sx
	real*8 ps,qs,as
	double complex z1a,z2a,z3a,z4a,z7a,z8a,ur,ui
	real*8 dsigt00,dsigl00,dsigp00,dsigi00
c	real*8 dsigp0,dsigi0,ctt
	real*8 zf(2,6),sigma_eep
	real*8 qv2,qv,qvc
	integer pn,pna(3)
	integer i

c real*4 for compatability with CERNLIB routine fint.
	real*4 px(3),pa(50)

	real*4 fint

	w=sqrt(ss)*1000.
	skc2=(w**2-Mk2-targ%Mrec_struck**2)**2-4.*Mk2*targ%Mrec_struck**2
	skc2=max(skc2,0.e0)
	skc=sqrt(skc2)/2./w
	q0=-(-q22-w**2+Mp2)/2./Mp
	q0c=(-q22+q0*Mp)/w
	qr=sqrt(q22)/q0c
	qv2=q22+q0**2
	qv=sqrt(qv2)
	qvc=Mp/w*qv
C ctt not used, and causes floating exception if skc=0 (i.e. skc2<0)
C	ctt=pi/skc/qvc*1.d+06
	aflx=skc/2./w/(w**2-Mp2)*(hbarc)**2*10000.
	aflxl=aflx*qr**2
	an=angl*180./pi
	x=cos(angl)
	sx=sin(angl)
	pn=3
	px(1)=ss
	px(2)=q22/1.d+06
	px(3)=an
	pna(1)=10
	pna(2)=11
	pna(3)=19
	ps=2.6
	qs=0.0
	as=0.0
	do i=1,10
	  pa(i)=ps
	  ps=ps+0.3
	enddo
	do i=11,21
	  pa(i)=qs
	  qs=qs+0.2
	enddo
	do i=22,40
	  pa(i)=as
	  as=as+10.
	enddo

	zf(1,1)=dble(fint(pn,px,pna,pa,zrff1))
	zf(1,2)=dble(fint(pn,px,pna,pa,zrff2))
	zf(1,3)=dble(fint(pn,px,pna,pa,zrff3))
	zf(1,4)=dble(fint(pn,px,pna,pa,zrff4))
	zf(1,5)=dble(fint(pn,px,pna,pa,zrff5))
	zf(1,6)=dble(fint(pn,px,pna,pa,zrff6))
	zf(2,1)=dble(fint(pn,px,pna,pa,ziff1))
	zf(2,2)=dble(fint(pn,px,pna,pa,ziff2))
	zf(2,3)=dble(fint(pn,px,pna,pa,ziff3))
	zf(2,4)=dble(fint(pn,px,pna,pa,ziff4))
	zf(2,5)=dble(fint(pn,px,pna,pa,ziff5))
	zf(2,6)=dble(fint(pn,px,pna,pa,ziff6))

	ur=(1.e0,0.e0)
	ui=(0.e0,1.e0)  

	z1a=zf(1,1)*ur+zf(2,1)*ui
	z2a=zf(1,2)*ur+zf(2,2)*ui
	z3a=zf(1,3)*ur+zf(2,3)*ui
	z4a=zf(1,4)*ur+zf(2,4)*ui
	z7a=zf(1,5)*ur+zf(2,5)*ui
	z8a=zf(1,6)*ur+zf(2,6)*ui
c
c	calculate free cross section using Saghai's form.
c
	dsigt00=aflx*(abs(z1a)**2+abs(z2a)**2
     &		+2.*real(conjg(z1a)*z2a)*x+0.5*sx**2*(abs(z3a)**2
     &		+abs(z4a)**2+2.*real(conjg(z1a)*z4a-conjg(z2a)*z3a
     &		+conjg(z3a)*z4a*x)))
	dsigl00=aflxl*epsi*(abs(z7a)**2+abs(z8a)**2
     &		+2.*real(conjg(z7a)*z8a)*x)
	dsigp00=aflx*epsi*(sin(theta))**2*cos(2.*phi)*
     &		(0.5*abs(z3a)**2+0.5*abs(z4a)**2+real(conjg(z1a)*z4a-
     &		conjg(z2a)*z3a+conjg(z3a)*z4a*x))
	dsigi00=aflx*sqrt(2.*qr**2*epsi*(1.+epsi))*
     &		sin(theta)*cos(phi)*real((z7a)*(conjg(z3a)-conjg(z2a)+
     &		conjg(z4a)*x)+z8a*(conjg(z1a)+conjg(z3a)*x+conjg(z4a)))
C	dsigp0=aflx*epsi*(sin(theta))**2*
C     &		(0.5*abs(z3a)**2+0.5*abs(z4a)**2+real(conjg(z1a)*z4a-
C     &		conjg(z2a)*z3a+conjg(z3a)*z4a*x))
C	dsigi0=aflx*sqrt(2.*qr**2*epsi*(1.+epsi))*
C     &		sin(theta)*real(z7a*(conjg(z3a)-conjg(z2a)+
C     &		conjg(z4a)*x)+z8a*(conjg(z1a)+conjg(z3a)*x+conjg(z4a)))

	sigma_eep=dsigt00+dsigl00+dsigp00+dsigi00
c	write(21,*) 'dsig',dsigt00,dsigl00/epsi,sigma_eep
c	write(21,*) 'dsig',dsigt00*ctt,dsigl00/epsi*ctt,sigma_eep*ctt
c	write(21,*) dsigp00,dsigi00,sigma_eep
c	write(21,*) dsigp0,dsigi0,sigma_eep
	return

	end

	subroutine eekeeks(ss,q22,angl,theta,phi,epsi,sigma_eep)

	include 'simulate.inc'

	real*8 ss,q22,angl,theta,phi,epsi
	real*8 w,skc2,skc,q0,q0c,qr,aflx,aflxl,an,x,sx
	real*8 as
	double complex z1a,z2a,z3a,z4a,z7a,z8a,ur,ui
	real*8 dsigt00,dsigl00,dsigp00,dsigi00
	real*8 zf(2,6),sigma_eep
	real*8 qv2,qv,qvc
c	real*8 ctt
	integer pn,pna(3),i

c real*4 for compatability with CERNLIB routine fint.
	real*4 px(3),pa(50)

	real*4 fint

	w=sqrt(ss)*1000.
	skc2=(w**2-Mk2-targ%Mrec_struck**2)**2-4.*Mk2*targ%Mrec_struck**2
	skc2=max(skc2,0.e0)
	skc=sqrt(skc2)/2./w
	q0=-(-q22-w**2+Mp2)/2./Mp
	q0c=(-q22+q0*Mp)/w
	qr=sqrt(q22)/q0c
	qv2=q22+q0**2
	qv=sqrt(qv2)
	qvc=Mp/w*qv
C ctt not used, and causes floating exception if skc=0 (i.e. skc2<0)
C	ctt=pi/skc/qvc*1.d+06
	aflx=skc/2./w/(w**2-Mp2)*(hbarc)**2*10000.
	aflxl=aflx*qr**2
	an=angl*180./pi
	x=cos(angl)
	sx=sin(angl)
	pn=3
	px(1)=ss
	px(2)=q22/1.d+06
	px(3)=an
	pna(1)=20
	pna(2)=10
	pna(3)=19
	as=0.0
	pa(1)=2.851
	pa(2)=2.898
	pa(3)=2.945
	pa(4)=2.991
	pa(5)=3.038
	pa(6)=3.085
	pa(7)=3.132
	pa(8)=3.320
	pa(9)=3.507
	pa(10)=3.695
	pa(11)=3.883
	pa(12)=4.070
	pa(13)=4.258
	pa(14)=4.446
	pa(15)=4.633
	pa(16)=4.821
	pa(17)=5.009
	pa(18)=5.196
	pa(19)=5.384
	pa(20)=5.572
	pa(21)=0.0
	pa(22)=0.250
	pa(23)=0.376
	pa(24)=0.520
	pa(25)=0.750
	pa(26)=1.000
	pa(27)=1.250
	pa(28)=1.500
	pa(29)=1.750
	pa(30)=2.000
	do i=31,49
	  pa(i)=as
	  as=as+10.
	enddo
	zf(1,1)=dble(fint(pn,px,pna,pa,zsrff1))
	zf(1,2)=dble(fint(pn,px,pna,pa,zsrff2))
	zf(1,3)=dble(fint(pn,px,pna,pa,zsrff3))
	zf(1,4)=dble(fint(pn,px,pna,pa,zsrff4))
	zf(1,5)=dble(fint(pn,px,pna,pa,zsrff5))
	zf(1,6)=dble(fint(pn,px,pna,pa,zsrff6))
	zf(2,1)=dble(fint(pn,px,pna,pa,zsiff1))
	zf(2,2)=dble(fint(pn,px,pna,pa,zsiff2))
	zf(2,3)=dble(fint(pn,px,pna,pa,zsiff3))
	zf(2,4)=dble(fint(pn,px,pna,pa,zsiff4))
	zf(2,5)=dble(fint(pn,px,pna,pa,zsiff5))
	zf(2,6)=dble(fint(pn,px,pna,pa,zsiff6))

	ur=(1.e0,0.e0)
	ui=(0.e0,1.e0)  

	z1a=zf(1,1)*ur+zf(2,1)*ui
	z2a=zf(1,2)*ur+zf(2,2)*ui
	z3a=zf(1,3)*ur+zf(2,3)*ui
	z4a=zf(1,4)*ur+zf(2,4)*ui
	z7a=zf(1,5)*ur+zf(2,5)*ui
	z8a=zf(1,6)*ur+zf(2,6)*ui
c
c	calculate free cross section using Saghai's form.
c
	dsigt00=aflx*(abs(z1a)**2+abs(z2a)**2
     &		+2.*real(conjg(z1a)*z2a)*x+0.5*sx**2*(abs(z3a)**2
     &		+abs(z4a)**2+2.*real(conjg(z1a)*z4a-conjg(z2a)*z3a
     &		+conjg(z3a)*z4a*x)))
	dsigl00=aflxl*epsi*(abs(z7a)**2+abs(z8a)**2
     &		+2.*real(conjg(z7a)*z8a)*x)
	dsigp00=aflx*epsi*(sin(theta))**2*cos(2.*phi)*
     &		(0.5*abs(z3a)**2+0.5*abs(z4a)**2+real(conjg(z1a)*z4a-
     &		conjg(z2a)*z3a+conjg(z3a)*z4a*x))
	dsigi00=aflx*sqrt(2.*qr**2*epsi*(1.+epsi))*
     &		sin(theta)*cos(phi)*real((z7a)*(conjg(z3a)-conjg(z2a)+
     &		conjg(z4a)*x)+z8a*(conjg(z1a)+conjg(z3a)*x+conjg(z4a)))

	sigma_eep=dsigt00+dsigl00+dsigp00+dsigi00
	return

	end

	subroutine phspwght(delta,theta,phi,weight)

	include 'simulate.inc'

	real*8 delta,theta,phi
	real*8 weight
	integer pn,pna(3),i
	real*8 ps,ts,ds

c real*4 for compatability with CERNLIB routine fint.
	real*4 px(3),pa(78)

	real*4 fint

	if(electron_arm.eq.1 .and. hadron_arm.eq.2)then	!e- in HMS, K in SOS
	  pn=2
	  px(1)=theta
	  px(2)=phi
	  pna(1)=20
	  pna(2)=50
	  ts=-0.03325
	  ps=-0.0735
	  do i=1,20
	    pa(i)=ts
	    ts=ts+0.0035
	  enddo
	  do i=21,70
	    pa(i)=ps
	    ps=ps+0.003
	  enddo
	  weight=dble(fint(pn,px,pna,pa,weightc))
	else if (electron_arm.eq.2 .and. hadron_arm.eq.1) then !e- in SOS,K in HmS
	  pn=3
	  px(1)=delta
	  px(2)=theta
	  px(3)=phi
	  pna(1)=8
	  pna(2)=40
	  pna(3)=30
	  ts=-0.063375
	  ps=-0.0435
	  ds=-17.5
	  do i=1,8
	    pa(i)=ds
	    ds=ds+5.
	  enddo
	  do i=9,48
	    pa(i)=ts
	    ts=ts+0.00325
	  enddo
	  do i=49,78
	    pa(i)=ps
	    ps=ps+0.003
	  enddo
	  weight=dble(fint(pn,px,pna,pa,weightd))
	else
	  write(6,*) 'electron_arm=',electron_arm,' and hadron_arm=',hadron_arm
	  write(6,*) 'eekeek.f has a phase space factor that is only defined for'
	  write(6,*) 'hms&sos case.  Need to update for other spectrometers.'
	  stop
	endif
	weight=max(weight,0.01e00)
	weight=max(100.e00/weight,1.0e00)
	return
	end
