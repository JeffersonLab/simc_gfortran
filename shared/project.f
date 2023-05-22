	subroutine project(x_new,y_new,z_drift,decay_flag,dflag,m2,ph,pathlen)

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

	include '../spectrometers.inc'
	include '../simulate.inc'

	real*8 x_new,y_new,z_drift,tmpdrift
	logical decay_flag		!check for decay
	logical dflag			!has particle decayed yet?

C Local declarations.
	real*8 p_spec,ph,m2, m_final
	real*8 z_decay,pathlen
	real*8 rph,rth1,rth
	real*8 beta,gamma,dlen
	real*8 ef,pf,pxf,pyf,pzf,pxr,pyr,pzr
	real*8 bx,by,bz,er,pr

	real*8 grnd

C ============================= Executable Code ================================

C Check for decay of particle.

	if (.not.decay_flag .or. dflag) then  !particle has already decayed, just drift
	  pathlen = pathlen + z_drift*sqrt(1+dxdzs**2+dydzs**2) 
	  x_new = x_new + dxdzs * z_drift
	  y_new = y_new + dydzs * z_drift
	
	else		!check for decay

	  p_spec = ph/(1.+dpps/100.)
	  beta = ph/sqrt(ph**2+m2)
	  gamma = 1./sqrt(1.-beta*beta)
	  dlen=ctau*beta*gamma   !1/e decay length (beta*c*tau*gamma)
	  if (z_drift.le.0) write(6,*) 'drift distance<0:  automatic decay! BAD!'

	  z_decay = -1.*dlen*log(1-grnd())

C Check if the decay is within the drift length (i.e. the drift lenght in
C z, times sqrt(1+dxdzs**2+dydzs**2) to correct for the true path of the ray.

	  if(z_decay .gt. z_drift*sqrt(1+dxdzs**2+dydzs**2)) then !no decay 

	    decdist = decdist + z_drift*sqrt(1+dxdzs**2+dydzs**2)
	    pathlen = pathlen + z_drift*sqrt(1+dxdzs**2+dydzs**2) 
	    x_new = x_new + dxdzs * z_drift
	    y_new = y_new + dydzs * z_drift

	  else		!DECAY.  Find out where and generate decay particle.

	    dflag = .true.
	    decdist = decdist + z_decay
	    pathlen = pathlen + z_decay
	    x_new = x_new + dxdzs * z_decay/sqrt(1+dxdzs**2+dydzs**2) 
	    y_new = y_new + dydzs * z_decay/sqrt(1+dxdzs**2+dydzs**2) 

C Generate center/mass decay angles and momenta.
	    rph = grnd()*2.*pi
	    rth1 = grnd()*2.-1.
	    rth = acos(rth1)

	    pr = 0.
	    m_final = Mmu ! default
	    if(abs(sqrt(m2) - Mpi).lt.2) pr = 29.783 ! pion decay
	    if(abs(sqrt(m2) - Mk).lt.2) then ! kaons
	       if(grnd().lt.0.7) then ! decay to muon plus neutrino
		  pr = 235.5 
	       else		! decay to two pions
		  pr = sqrt(Mk**2 / 4. - Mpi**2)
		  m_final = Mpi
	       endif
	    endif
	    if(pr.eq.0.) then
	     write(6,'(''error, cannot decay particle with'',
     >        '' mass='',f8.2)') sqrt(m2)
	     stop
	    endif
	    er = sqrt(m_final**2 + pr**2)
	    pxr = pr*sin(rth)*cos(rph)
	    pyr = pr*sin(rth)*sin(rph)
	    pzr = pr*cos(rth)
	    m2 = m_final**2	!need mass-squared for multiple scattering.
	    Mh2_final = m2	!for ntuple


C Boost to Lab frame, calculate new angles and momentum, finish drift

	    bx = -beta * dxdzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    by = -beta * dydzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    bz = -beta *   1.  / sqrt(1. + dxdzs**2 + dydzs**2)
	    call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)
	    dxdzs = pxf/pzf
	    dydzs = pyf/pzf
	    dpps = 100.*(pf/p_spec-1.)
	    ph=pf

C We've already drifted z_decay/sqrt(1+dxdzs**2+dydzs**2) - now do the rest.

	    tmpdrift = z_drift - z_decay/sqrt(1+dxdzs**2+dydzs**2) 
	    pathlen = pathlen + tmpdrift*sqrt(1+dxdzs**2+dydzs**2) 
	    x_new = x_new + dxdzs * tmpdrift
	    y_new = y_new + dydzs * tmpdrift

	  endif	  !if decayed
	endif	!if need to check for decay

	return
	end
