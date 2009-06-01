	subroutine generate_rho(vertex,success)

C+_____________________________________________________________________________
!  generate_rho:  This routine takes care of the event genration scheme
!  for rho electroproduction.  This generates in 4pi in the virtual-photon
!  target center of mass frame and then boosts back to the lab.
!
!   This routine should be called before complete_ev.

C-_____________________________________________________________________________

	implicit none
	include 'simulate.inc'

	type(event):: vertex

	real*8 Epux,Epuy,Epuz !scattered electron unit vectors
	real*8 nu,Q2,q        !the usual
	real*8 qux,quy,quz    !q-vector unit vectors

	real*8 bx,by,bz,betacm,gammacm   !virtual photon-nucleon center of mass
	
	real*8 s  !W2 in cm - including possible Fermi motion effects
	real*8 Erhocm,Prhocm


	real*8 rph,rth1,rth
	real*8 pxr,pyr,pzr,er
	real*8 pxf,pyf,pzf,pf,ef

	real*8 grnd

	logical success


C Local declarations.
	
C ============================= Executable Code ================================

C Calculate some electron junk
C Lab coordinates!!!

	Epux = vertex%ue%x
	Epuy = vertex%ue%y
	Epuz = vertex%ue%z
	

	nu = vertex%nu
	Q2 = vertex%Q2
	q = vertex%q


	qux = vertex%uq%x
	quy = vertex%uq%y
	quz = vertex%uq%z


	bx = -(q*qux+pferx*pfer)/(nu+efer)
	by = -(q*quy+pfery*pfer)/(nu+efer)
	bz = -(q*quz+pferz*pfer)/(nu+efer)

	betacm = sqrt(bx**2+by**2+bz**2)
	if(betacm.gt.1.0) then
c	   write(6,*) 'generate rho: beta_cm greater than 1!'
	   success=.false.
	   return
	endif
	gammacm = 1./sqrt(1.0-betacm**2)
C Start calculating rho stuff

	s = -Q2 + efer**2 + 2.*nu*efer


	Mh = Mrho
! DJG give the rho mass some width (non-relativistic Breit-Wigner)
	Mh = Mh + 0.5*150.2*tan((2.*grnd()-1.)*atan(2.*500./150.2))
	Mh2 = Mh*Mh
	ntup%rhomass=Mh
	Mh2_final = Mh

	Erhocm = (s + Mh2 - targ%Mrec_struck**2)/2./sqrt(s)

	if(Erhocm.lt.Mh) then
c	   write(6,*) 'Below rho threshold - sayonara sucker!'
	   success = .false.
	   return
	endif

	Prhocm = sqrt(Erhocm**2-Mh2)

! DJG generate angles - flat in phi and cos(theta)

	rph = grnd()*2.*pi
 	rth1 = grnd()*2.-1.
	rth = acos(rth1)

        pxr = Prhocm*sin(rth)*cos(rph)
        pyr = Prhocm*sin(rth)*sin(rph)
        pzr = Prhocm*cos(rth)
	er = Erhocm

C OK, we've energy, momentum and angles%  Now we need to boost this sucker
C back to the lab%


C Boost to Lab frame%

	call loren(gammacm,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)


	vertex%p%delta = 100.*(pf/spec%p%P-1.)
	
	vertex%p%P = pf
	vertex%p%E = ef

	vertex%up%x = pxf/pf
	vertex%up%y = pyf/pf
	vertex%up%z = pzf/pf

C DJG These aren't really used for anything, so let me just assign
C DJG them to theta and phi in the lab
	vertex%p%xptar = acos(pzf/pf)
	vertex%p%yptar = atan2(pyf,pxf)

	success=.true.
	return
	end



