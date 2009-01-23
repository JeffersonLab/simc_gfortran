	subroutine project(x_new,y_new,z_drift,decay_flag,dflag,m2,ppi)

C+______________________________________________________________________________
!
! PROJECT - Calculate new track transverse coordinates after drifting in a field
!   free region for a distance z_drift from current track location.
!
! Added 9/15/98: Check for decay of particle. dflag is true of the particle
! has already decayed, so check for decay if dflag .eq. .false.  Take into
! account additional path length for rays not parallel to the z axis.
!
C-______________________________________________________________________________

	implicit none

	include '../track.inc'
	include '../simulate.inc'

	real*8 x_new,y_new,z_drift
	logical*4 decay_flag		!check for decay
	logical*4 dflag			!has particle decayed yet?

C Local declarations.
	real*8 p_spec,ppi,m2
	real*8 z_decay
	real*8 rph,rth1,rth
	real*8 beta,gamma,dlen
	real*8 ef,pf,pxf,pyf,pzf,pxr,pyr,pzr
	real*8 bx,by,bz,er,pr

	real*8 grnd

C ============================= Executable Code ================================


C Check for decay of particle.

	if (.not.decay_flag .or. dflag) then  !particle has already decayed, just drift

	  x_new = x_new + dxdzs * z_drift
	  y_new = y_new + dydzs * z_drift
	
	else		!check for decay

	  p_spec = ppi/(1.+dpps/100.)
	  beta = ppi/sqrt(ppi**2+m2)
	  gamma = 1./sqrt(1.-beta*beta)
	  dlen=ctau*beta*gamma   !1/e decay length (beta*c*tau*gamma)
	  if (z_drift.le.0) write(6,*) 'drift distance<0:  automatic decay! BAD!'

	  z_decay = -1.*dlen*log(1-grnd())

C Check if the decay is within the drift length (i.e. the drift lenght in
C z, times sqrt(1+dxdzs**2+dydzs**2) to correct for the true path of the ray.

	  if(z_decay .gt. z_drift*sqrt(1+dxdzs**2+dydzs**2)) then !no decay 

	    decdist(30)=decdist(30)+z_drift
	    x_new = x_new + dxdzs * z_drift
	    y_new = y_new + dydzs * z_drift

	  else		!DECAY.  Find out where and generate decay particle.

	    dflag = .true.
	    decdist(30)=decdist(30)+z_decay
	    x_new = x_new + dxdzs * z_decay
	    y_new = y_new + dydzs * z_decay

C Generate center/mass decay angles and momenta.
	    rph = grnd()*2.*pi
	    rth1 = grnd()*2.-1.
	    rth = acos(rth1)
	    er = 109.787
	    pr = 29.783
	    pxr = 29.783*sin(rth)*cos(rph)
	    pyr = 29.783*sin(rth)*sin(rph)
	    pzr = 29.783*cos(rth)
	    m2 = 105.67 **2	!need mass-squared for multiple scattering.
	    Mh2_final = m2	!for ntuple


C Boost to Lab frame, calculate new angles and momentum, finish drift

	    bx = beta * dxdzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    by = beta * dydzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    bz = beta *   1.  / sqrt(1. + dxdzs**2 + dydzs**2)
	    call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)
	    dxdzs = pxf/pzf
	    dydzs = pyf/pzf
	    dpps = 100.*(pf/p_spec-1.)
	    ppi=pf
	    x_new = x_new + dxdzs * (z_drift-z_decay)
	    y_new = y_new + dydzs * (z_drift-z_decay)

	  endif	  !if decayed
	endif	!if need to check for decay

	return
	end
