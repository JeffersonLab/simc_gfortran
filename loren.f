	subroutine loren(gam,bx,by,bz,e,x,y,z,ef1,pxf,pyf,pzf,pf1)

! This subroutine takes a particle with energy/momentum (e,x,y,z), and
! gives back it's momentum is a frame with velocity (bx,by,bz) WITH RESPECT
! TO THE INITIAL FRAME, so be careful of the sign.  For example, when
! we model decay of particles, we generate the decay products in the
! rest frame of the initial particle, and then boost them back into the
! lab frame (along the directio of the inital particle).  So the boost
! vector (bx,by,bz) is opposite of the initial particle's momentum, because
! we are boosting from the moving particle's rest frame back into the lab
! frame.

	implicit none

	real*8	gam,bx,by,bz,e,x,y,z
	real*8  gam1,pxf,pyf,pzf,pf1,ef1

	gam1=gam**2/(1.+gam)
	ef1 = gam*(e-bx*x-by*y-bz*z)
	pxf = (1+gam1*bx**2)*x + gam1*bx*(by*y+bz*z) - gam*bx*e
	pyf = (1+gam1*by**2)*y + gam1*by*(bx*x+bz*z) - gam*by*e
	pzf = (1+gam1*bz**2)*z + gam1*bz*(by*y+bx*x) - gam*bz*e
	pf1 = sqrt(pxf**2+pyf**2+pzf**2)

	return
	end
