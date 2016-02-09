	real*8 function sigep(vertex)

! elastic cross section, units are set by sigMott.f (microbarn/sr)

	implicit none
	include 'simulate.inc'

	real*8	q4sq,qmu4mp,W1p,W2p,Wp,GE,GM,sigMott
	type(event):: vertex

	q4sq	= vertex%Q2
	call fofa_best_fit(-q4sq/hbarc**2,GE,GM)
	qmu4mp	= q4sq/4./Mp2
	W1p = GM**2*qmu4mp
	W2p = (GE**2+GM**2*qmu4mp)/(1.0+qmu4mp)
	Wp = W2p + 2.*W1p*tan(vertex%e%theta/2.)**2
	sigep = sigMott(vertex%e%E,vertex%e%theta,vertex%q2) *
     >		vertex%e%E/vertex%Ein * Wp
	if (debug(5)) write(6,*) 'sigMott',GE,GM
	return
	end

!-------------------------------------------------------------------------

	real*8 function deForest(ev)

	implicit none
	include 'simulate.inc'

	type(event):: ev

	real*8	q4sq,ebar,qbsq,GE,GM,f1,kf2,qmu4mp,sigMott,sin_gamma,cos_phi
	real*8	pbarp,pbarq,pq,qbarq,q2,f1sq,kf2_over_2m_allsq
	real*8	sumFF1,sumFF2,termC,termT,termS,termI,WC,WT,WS,WI,allsum

! Compute deForest sigcc cross-section, according to value of DEFOREST_FLAG:
!	Flag = 0    -- use sigcc1
!	Flag = 1    -- use sigcc2
!	Flag = -1   -- use sigcc1 ONSHELL, replacing Ebar with E = E'-nu
!			and qbar with q (4-vector)
! N.B.	Beware of deForest's metric when taking all those 4-vector inner
!	products in sigcc2 ... it is (-1,1,1,1)!
!	Here, I've defined all the inner products with the regular signs,
!	and then put them in the structure function
!	formulas with reversed signs compared to deForest.
!
! Note that this can be called with either the vertex or recon event records,
! but like a good function, it does not modify anything in those records.
!
! JRA: Note, let me make a comment about the following, because I was badly
! confused about it the first time through.  The 6-fold cross section is:
! d6sigma = K * S(E,p) * sigma_eN.  This routine returns a quantity 'deForest'
! which is d6sigma/S(E,p) = K*sigma_eN
! = (E'*p')*sigma_mott*[sum of terms like a_i*W_i].
! Each W_i term has E'/Ebar in it.  In the following code, the E'/Ebar appears
! to have been removed from the W_i terms, and combined with the K=E'*p' term
! to give an overall p'/Ebar in the final expression.
! More importantly, since the Spectral Function has been divided out, the
! units are MeV^4 * the 6-fold cross section units, and so for sigma_mott
! in microbarn/sr, 'deForest' is in microbarn*MeV^2/sr^2, giving the correct
! cross seciton units once we multiply by S(E,p) which is in MeV^-4

	q4sq = -ev%Q2
	q2 = ev%q**2
	if (deForest_flag.ge.0) then
	  ebar = sqrt(ev%Pm**2 + Mh2)
	  qbsq = (ev%p%E-ebar)**2 - q2
	else
	  ebar = ev%p%E - ev%nu
	  qbsq = q4sq
	endif

	sin_gamma = 1. - (ev%uq%x*ev%up%x+ev%uq%y*ev%up%y+ev%uq%z*ev%up%z)**2
	if (sin_gamma.lt.0) then
	  write(6,'(1x,''WARNING: deForest came up with sin_gamma = '',f10.3,'' at event '',i10)') sin_gamma, nevent
	  sin_gamma = 0.0
	endif
	sin_gamma = sqrt(sin_gamma)
	cos_phi = 0.0
	if (sin_gamma.ne.0) cos_phi=(ev%uq%y*(ev%uq%y*ev%up%z-ev%uq%z*ev%up%y)
     >		- ev%uq%x*(ev%uq%z*ev%up%x-ev%uq%x*ev%up%z))
     >		/ sin_gamma / sqrt(1.-ev%uq%z**2)
	if (abs(cos_phi).gt.1.) then		!set to +/-1, warn if >1.e-10
	  cos_phi = sign(1.0,cos_phi)
	  if ( (abs(cos_phi)-1.) .gt. 1.e-10) write(6,*)
     >		'WARNING: deForest give cos_phi = ',cos_phi,' at event',nevent
	endif

	call fofa_best_fit(q4sq/hbarc**2,GE,GM)
	qmu4mp = q4sq/4./Mp2
	f1 = (GE-GM*qmu4mp)/(1.0-qmu4mp)
	kf2 = (GM-GE)/(1.0-qmu4mp)
	f1sq = f1**2
	kf2_over_2m_allsq = kF2**2/4./Mh2

	termC = (q4sq/q2)**2
	termT = (tan(ev%e%theta/2.))**2 - q4sq/2./q2
	termS = tan(ev%e%theta/2.)**2 - (q4sq/q2)*cos_phi**2
	termI = (-q4sq/q2)*sqrt(tan(ev%e%theta/2.)**2-q4sq/q2)*cos_phi

	if (deForest_flag.le.0) then
	  sumFF1 = (f1 + kf2)**2
	  sumFF2 = f1sq - qbsq*kf2*kf2/4./Mh2

	  WC = ((ebar+ev%p%E)**2)*sumFF2 - q2*sumFF1
	  WT = -2*qbsq*sumFF1
	  WS = 4*(ev%p%P**2)*(sin_gamma**2)*sumFF2
	  WI = -4*(ebar+ev%p%E)*ev%p%P*sin_gamma*sumFF2
	else
	  pbarp=ebar*ev%p%E-ev%p%p*(ev%up%x*ev%Pmx+ev%up%y*ev%Pmy+ev%up%z*ev%Pmz)
	  pbarq=ebar*ev%nu-ev%q*(ev%uq%x*ev%Pmx+ev%uq%y*ev%Pmy+ev%uq%z*ev%Pmz)
	  pq = ev%p%E*ev%nu - ev%p%p * ev%q * (ev%up%x*ev%uq%x +
     >		ev%up%y*ev%uq%y + ev%up%z*ev%uq%z)
	  qbarq = (ev%p%E-ebar)*ev%nu - q2

	  WC =	(ebar*ev%p%E+(-pbarp+Mh2)/2.) * f1sq - q2*f1*kf2/2.
     >		- ( (-pbarq*ev%p%E-pq*ebar)*ev%nu + ebar*ev%p%E*q4sq +
     >		pbarq*pq - (-pbarp-Mh2)/2.*q2 ) * kf2_over_2m_allsq
	  WT =	- (-pbarp+Mh2)*f1sq - qbarq*f1*kF2
     >		+ (2.*pbarq*pq + (-pbarp-Mh2)*q4sq) * kf2_over_2m_allsq
	  WS =	(ev%p%p*sin_gamma)**2 * (f1sq - q4sq*kf2_over_2m_allsq)
	  WI =	ev%p%p*sin_gamma * ( -(ebar+ev%p%E)*f1sq + ((-pbarq-pq)*
     >		ev%nu+(ebar+ev%p%E)*q4sq)*kf2_over_2m_allsq )
	endif

	allsum = termC*WC + termT*WT + termS*WS + termI*WI
	if (deForest_flag.le.0) allsum = allsum/4.0
	deForest = sigMott(ev%e%E,ev%e%theta,ev%Q2)*ev%p%P*allsum/ebar	!microbarn*(MeV/sr)**2
	if (debug(5)) write(6,*) 'deForest',GE,GM

	return
	end

!-------------------------------------------------------------------------

	subroutine fofa_best_fit(qsquar,GE,GM)

*	csa 9/14/98 -- This calculates the form factors Gep and Gmp using
*	Peter Bosted's fit to world data (Phys. Rev. C 51, 409, Eqs. 4
*	and 5 or, alternatively, Eqs. 6)

	implicit none
	include 'simulate.inc'

	real*8  qsquar,GE,GM,mu_p
	real*8  Q,Q2,Q3,Q4,Q5,denom

	mu_p = 2.793

	Q2 = -qsquar*(hbarc**2.)*1.e-6
	Q  = sqrt(max(Q2,0.e0))

	Q3 = Q**3.
	Q4 = Q**4.
	Q5 = Q**5.

* Use Eqs 4, 5:
	denom = 1. + 0.62*Q + 0.68*Q2 + 2.8*Q3 + 0.83*Q4
	GE = 1./denom
	denom = 1. + 0.35*Q + 2.44*Q2 + 0.5*Q3 + 1.04*Q4 + 0.34*Q5
	GM = mu_p/denom

* OR Eqs 6:
*	denom = 1. + 0.14*Q + 3.01*Q2 + 0.02*Q3 + 1.20*Q4 + 0.32*Q5
*	GE = 1./denom
*	GM = mu_p/denom

	return
	end


!-------------------------------------------------------------------------

	real*8 function sigMott(e0,theta,Q2)

	implicit none
	include 'constants.inc'

	real*8 e0,theta,Q2
	real*8 sig

! The Mott cross section (for a point nucleus) in microbarns/sr.

	sig = (2.*alpha*hbarc*e0*cos(theta/2.)/Q2 )**2
	sigMott = sig*1.d4	!fm**2 --> microbarns

	return
	end
