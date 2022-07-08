	subroutine rho_decay(orig,p_spec,epsilon,success)

C+_____________________________________________________________________________
!  rho_decay:  This routine takes the electroproduced rho0 and
!   produces the decaying pi+ (or pi-) that we actually detect in the 
!   spectrometer.
!
!   This routine should be called before mc_hms or mc_sos.
!
!   The outgoing pion will be thrown flat in phi and according to
!   sin(theta)**2+epsilon*R*cos(theta)**2 (normalized) in the rho
!   rest frame. 
!
!   Note that this routine depends on R = sigma_L/sigma_T (and epsilon).
C-_____________________________________________________________________________

	implicit none

	include 'simulate.inc'

	type(event):: orig
C Local declarations.
	real*8 p_spec,ph
	real*8 rph,rth1,rth
	real*8 beta,gamma,dlen
	real*8 ef,pf,pxf,pyf,pzf,pxr,pyr,pzr
	real*8 bx,by,bz,er,pr

	real*8 th_oop,th_inp,cos_th_inp
	real*8 pzprime

	real*8 epsilon,R_rho

	real*8 ctaurho !do not use "real" ctau (i.e. for pion) since
	               !we're doing rho decay here
	real*8 norm,dist
	real*8 grnd

	logical success

	parameter (ctaurho=1.31467e-13)
C ============================= Executable Code ================================

C Calulate R = sigmaL/sigmaT
C Warning!  Note that if you change the parameterization of R in
C rho-physics, do it here as well for self consistency!
C This parameterization is taken from HERMES data...

	R_rho = 0.33*(orig%Q2/mrho2)**(0.61)
	

	ph = orig%p%P
	beta = ph/sqrt(ph**2+Mh2)
	gamma = 1./sqrt(1.-beta*beta)
	dlen=ctaurho*beta*gamma	!1/e decay length (beta*c*tau*gamma)

C Generate center/mass decay angles and momenta.
	rph = grnd()*2.*pi
 100	rth1 = grnd()*2.-1.
	rth = acos(rth1)


	norm = (1.0+2.0*epsilon*R_rho)*grnd()
	dist = sin(rth)**2+2.0*epsilon*R_rho*cos(rth)**2
	if(dist.lt.norm) goto 100
c       write(6,*) 'now decaying the rho',ctaurho,mh2
	ntup%rhotheta=rth
c	er = 384.65
c	pr = 358.4353
	er = ntup%rhomass/2.0
	if(er.lt.Mpi) then
	   success=.false.
	   return
	endif
	pr = sqrt(er**2-Mpi2)
	pxr = pr*sin(rth)*cos(rph)
	pyr = pr*sin(rth)*sin(rph)
	pzr = pr*cos(rth)

C Boost to Lab frame, calculate new angles and momentum, finish drift

	bx = -beta * orig%up%x
	by = -beta * orig%up%y
	bz = -beta * orig%up%z
	call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)


C DJG check if pion is at least going toward to spectrometer!
C DJG If not, I'm outta here.
C DJG Need to rotate coordinate sytem counterclockwise (clockwise) by angle
C DJG theta_spec about x-axis (points DOWN) for HMS (SOS).

	if(hadron_arm.eq.1.or.hadron_arm.eq.3) then
	   pzprime = pzf*cos(spec%p%theta)-pyf*sin(spec%p%theta)
	elseif (hadron_arm.eq.2.or.hadron_arm.eq.4.or.hadron_arm.eq.5) then
	   pzprime = pzf*cos(spec%p%theta)+pyf*sin(spec%p%theta)
	else
	   write(6,*) 'Unknown spectrometer setup dude!'
	   stop
	endif

	if(pzprime.lt.0.0) then
c	   write(6,*) 'Pion from rho decay moving away from spectrometer'
c	   write(6,*) 'I guess this event is a no-go!'
	   success = .false.
	   return
	endif
	   
C DJG If pion's heading for spectrometer, try to calculate xptar and yptar

	th_oop=asin(pxf/pf)
        cos_th_inp = pzf/pf/cos(th_oop)
	th_inp = acos(cos_th_inp)
	if(hadron_arm.eq.1.or.hadron_arm.eq.3) then
	   if(pyf.lt.0.0) then
	      th_inp = spec%p%theta-th_inp 
	   else
	      th_inp = spec%p%theta+th_inp
	   endif
	elseif(hadron_arm.eq.2.or.hadron_arm.eq.4.or.hadron_arm.eq.5) then
	   if(pyf.gt.0.0) then
	      th_inp = th_inp - spec%p%theta
	   else
	      th_inp = th_inp + spec%p%theta
	   endif
	else
	   write(6,*) 'Rho decay not set up for your spectrometer!'
	   stop
	endif

C DJG This should never happen!!!

	if((th_inp.gt.pi/2.).or.(th_oop.gt.pi/2.)) then
	   write(6,*) 'Decay pion going backwards - fahgettaboutit!'
	   success=.false.
	   return
	endif

	orig%up%x = pxf/pf
	orig%up%y = pyf/pf
	orig%up%z = pzf/pf

	orig%p%xptar = tan(th_oop)
	orig%p%yptar = tan(th_inp)
	orig%p%delta = 100.*(pf/p_spec-1.)
	
	orig%p%P = pf

	Mh = Mpi
	Mh2 = Mpi2

	Mh2_final = Mh2

	orig%p%E = sqrt(orig%p%P**2 + Mh2)

C Calculate "physics" angles
	call physics_angles(spec%p%theta,spec%p%phi,
     &     orig%p%xptar,orig%p%yptar,orig%p%theta,orig%p%phi)

	return
	end



