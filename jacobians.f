	subroutine transform_to_cm(vertex,main,
     &		gstar,bstar,bstarx,bstary,bstarz,
     &		nustar,qstar,qstarx,qstary,qstarz,
     &		ehadcm,phadcm,phadcmx,phadcmy,phadcmz,
     &		ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz,
     &		etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz,
     &		thetacm,phicm,phiqn,jacobian,jac_old)

* GENERATE_CM
* This routine determines the transformation from the lab variables to
* the PHOTON-NUCLEON CENTER-OF-MASS frame, and returns the center-of-mass
* four-vectors for the q-vector, beam, and produced hadron (pion or kaon).

	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* Variables to be returned:
	real*8 gstar,bstar,bstarx,bstary,bstarz		!beta of boost toC.M.
	real*8 nustar,qstar,qstarx,qstary,qstarz	!q in C.M.
	real*8 ehadcm,phadcm,phadcmx,phadcmy,phadcmz	!p_hadron in C.M.
	real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz	!p_beam in C.M.
	real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz	!p_hadron in C.M.
	real*8 thetacm,phicm,phiqn
	real*8 jacobian,jac_old

* Temporary variables:  Boost related
	real*8 tcos,tsin		!cos/sin of theta between ppi and q
	real*8 tfcos,tfsin		!cos/sin of theta between pfermi and q
	real*8 cospq,sinpq
	real*8 dummy,zero
	real*8 qx,qy,qz,px,py,pz
	real*8 tmp_x_x,tmp_x_y,tmp_x_z
	real*8 tmp_y_x,tmp_y_y,tmp_y_z
	real*8 p_tmp_x,p_tmp_y
	real*8 pbeam,beam_tmpx,beam_tmpy,beam_tmpz !lab variables to be boosted.
	real*8 phadx,phady,phadz		   ! "       "      "      "
	real*8 ptarx,ptary,ptarz		   ! "       "      "      "
	real*8 tmp2_x_x,tmp2_x_y,tmp2_x_z
	real*8 tmp2_y_x,tmp2_y_y,tmp2_y_z
	real*8 tmp2_z_x,tmp2_z_y,tmp2_z_z
	real*8 phadcm_tmp2x,phadcm_tmp2y,phadcm_tmp2z

* Temporary variables:  Jacobian related
	real*8 square_root,dp_dcos_num,dp_dcos_den,dp_dcos
	real*8 dp_dphi_num,dp_dphi_den,dp_dphi
	real*8 dt_dcos_lab,dt_dphi_lab,psign
	real*8 dpxdphi,dpydphi,dpxdcos,dpydcos,dpzdcos,dpzdphi
	real*8 dpxnewdphi,dpynewdphi,dpxnewdcos,dpynewdcos
	real*8 dpznewdphi,dpznewdcos
	real*8 dphicmdphi,dphicmdcos

* f's and fer indicate fermi momenta, s, star or cm CM system
	tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y+vertex%up%z*vertex%uq%z
	if(tcos-1..gt.0..and.tcos-1..lt.1.e-8)tcos=1.0
	tsin=sqrt(1.-tcos**2)

	tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0
	tfsin=sqrt(1.-tfcos**2)

	cospq = cos(main%phi_pq)
	sinpq = sin(main%phi_pq)

* JRA:The transformation is calculated starting in the coord. system used
* in the fpi/nucpi replay (see event.f), where x=right, y=down, z=along beam.
* We must convert from SIMC coords to these first.
	qx = -vertex%uq%y		!'right'
	qy =  vertex%uq%x		!'down'
	qz =  vertex%uq%z
	px = -pfery
	py =  pferx
	pz =  pferz

	dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	tmp_x_x = -qx*qz/dummy
	tmp_x_y = -qy*qz/dummy
	tmp_x_z = (qx**2 + qy**2)/dummy

	dummy   = sqrt(qx**2 + qy**2)
	tmp_y_x =  qy/dummy
	tmp_y_y = -qx/dummy
	tmp_y_z =  0.0

	p_tmp_x = pfer*(px*tmp_x_x + py*tmp_x_y + pz*tmp_x_z)
	p_tmp_y = pfer*(px*tmp_y_x + py*tmp_y_y + pz*tmp_y_z)

	if(p_tmp_x.eq.0.)then
	  phiqn=0.
	else
	  phiqn = atan2(p_tmp_y,p_tmp_x)
	endif
	if(phiqn.lt.0.)phiqn = phiqn+2.*pi

* get beam in "q" system.

	pbeam = vertex%Ein
	beam_tmpx = pbeam*tmp_x_z
	beam_tmpy = pbeam*tmp_y_z
	beam_tmpz = pbeam*vertex%uq%z

	bstar=sqrt((vertex%q+pfer*tfcos)**2+(pfer*tfsin)**2)/(efer+vertex%nu)
	if (abs(bstar).gt.1.) write(6,'(a,4f15.5)') 'q,nu,efer,beta=',vertex%q,vertex%nu,efer,bstar
!	if (bstar.gt.1. .and. bstar.lt.1.001) bstar=0.9999999
!	if (bstar.lt.-1. .and. bstar.gt.-1.001) bstar=-0.999999
	gstar=1./sqrt(1. - bstar**2)

	bstarz = (vertex%q+pfer*tfcos)/(efer+vertex%nu)
	bstarx = p_tmp_x/(efer+vertex%nu)
	bstary = p_tmp_y/(efer+vertex%nu)

* DJG Boost the beam to CM.

	call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_tmpx,
     >		beam_tmpy,beam_tmpz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* DJG: Boost virtual photon to CM.

	zero =0.e0
	call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >		zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost pion to CM.

	phadz = vertex%p%P*tcos
	phadx = vertex%p%P*tsin*cospq
	phady = vertex%p%P*tsin*sinpq
	call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >		phadx,phady,phadz,ehadcm,phadcmx,phadcmy,phadcmz,phadcm)
	thetacm = acos((phadcmx*qstarx+phadcmy*qstary+phadcmz*qstarz)/phadcm/qstar)

* DJG: Boost target nucleon to CM.

	ptarz = pfer*tfcos
	ptarx = p_tmp_x
	ptary = p_tmp_y
	call loren(gstar,bstarx,bstary,bstarz,efer,
     >		ptarx,ptary,ptarz,etarcm,ptarcmx,ptarcmy,ptarcmz,ptarcm)

* Thetacm is defined as angle between phadcm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

	dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >		(qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >		+(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
	tmp2_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
	tmp2_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
	tmp2_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

	dummy = sqrt((tmp2_y_y*qstarz-tmp2_y_z*qstary)**2
     >		+(tmp2_y_z*qstarx-tmp2_y_x*qstarz)**2
     >		+(tmp2_y_x*qstary-tmp2_y_y*qstarx)**2)
	tmp2_x_x = (tmp2_y_y*qstarz-tmp2_y_z*qstary)/dummy
	tmp2_x_y = (tmp2_y_z*qstarx-tmp2_y_x*qstarz)/dummy
	tmp2_x_z = (tmp2_y_x*qstary-tmp2_y_y*qstarx)/dummy

	tmp2_z_x = qstarx/qstar
	tmp2_z_y = qstary/qstar
	tmp2_z_z = qstarz/qstar

	phadcm_tmp2x = phadcmx*tmp2_x_x + phadcmy*tmp2_x_y + phadcmz*tmp2_x_z
	phadcm_tmp2y = phadcmx*tmp2_y_x + phadcmy*tmp2_y_y + phadcmz*tmp2_y_z
	phadcm_tmp2z = phadcmx*tmp2_z_x + phadcmy*tmp2_z_y + phadcmz*tmp2_z_z

	phicm = atan2(phadcm_tmp2y,phadcm_tmp2x)
	if(phicm.lt.0.) phicm = 2.*pi+phicm


* While we're here, and have all of the above factors available, We'll
* calculate the Jacobian:


* Now, Henk's CM -> LAB transformation   OLD VERSION OF THE JACOBIAN!!!!
* HPB: Brauel's cross section is expressed in an invariant way as dsigma/dtdphi,
* HPB: with a virtual photon flux factor based on a full cross section, written
* HPB: as: dsigma/dQ2dWdtdphi. One can convert this cross section to dsigma/dom+
* HPB: in the lab frame or the cm frame by multiplying with the factor q*p_pi/3+
* HPB: (with the variables in the lab frame, or the cm frame)
* HPB: this factor is built from dt/dcos(theta)=2pq, and a factor 1/2pi, which
* HPB: comes from the conversion of Brauel's virtual photon flux factor and the
* HPB: one we use, based on dsigma/dE_e'dOmega_e'dOmega_pi
* HPB: see Devenish and Lyth, Phys. Rev.  C D5, 47(1972)
*
* JRA: Simpler (I hope) explanation:  sig is d2sigma/dt/dphi_cm. To get
* d2sigma/domega_cm (=s2cm), we multiply by dt/domega_cm = 1/2/pi*dt/dcos(theta+
* = 1/pi*q_cm*p_cm.  We can then get d5sigma/dE_e'/dOmega_e'/dOmega_pi_cm by
* multiplying by gammav.  Either of these can be converted to the dOmega_lab
* by multiplying by dOmega_lab/dOmega_cm = dcos(theta_lab)/dcos(theta_cm)
* (dphi_cm/dphi_lab=1 since the frames are collinear).  This is 'factor',
* which is equal to dt/dcos(theta_cm) / dt/dcos(theta_lab).  Hence, the
* following cross sections can be calculated:
*	s2cm =       sig*qcm*ppicm       /pi
*	s5cm =gammav*sig*qcm*ppicm       /pi
*	s2lab=       sig*q  *ppi  *factor/pi
*	s5lab=gammav*sig*q  *ppi  *factor/pi

* DJG The above assumes target at rest.  We need a full blown Jacobian:
* DJG: | dt/dcos(theta) dphicm/dphi - dt/dphi dphicm/dcos(theta) |
* DJG: Calculate dt/d(cos(theta)) and dt/dphi for the general case.

	psign = cos(phiqn)*cospq+sin(phiqn)*sinpq

	square_root = vertex%q + pfer*tfcos - vertex%p%P*tcos
	dp_dcos_num = vertex%p%P + (vertex%p%P**2*tcos -
     >		psign*pfer*vertex%p%P*tfsin*tcos/tsin)/square_root
	dp_dcos_den = ( (vertex%nu+efer-vertex%p%E)*vertex%p%P/vertex%p%E +
     >		vertex%p%P*tsin**2-psign*pfer*tfsin*tsin )/square_root - tcos
	dp_dcos = dp_dcos_num/dp_dcos_den

	dp_dphi_num = pfer*vertex%p%P*tsin*tfsin*(cos(phiqn)*sinpq-
     >		sin(phiqn)*cospq)/square_root
	dp_dphi_den = tcos + (pfer*tsin*tfsin*psign - vertex%p%P*tsin**2
     >		- (vertex%nu+efer-vertex%p%E)*vertex%p%P/vertex%p%E)/square_root
	dp_dphi = dp_dphi_num/dp_dphi_den

	dt_dcos_lab = 2.*(vertex%q*vertex%p%P +
     >		(vertex%q*tcos-vertex%nu*vertex%p%P/vertex%p%E)*dp_dcos)

	dt_dphi_lab = 2.*(vertex%q*tcos-vertex%nu*vertex%p%P/vertex%p%E)*dp_dphi

* DJG: Now calculate dphicm/dphi and dphicm/d(cos(theta)) in the most
* DJG: excruciating way possible.

	dpxdphi = vertex%p%P*tsin*(-sinpq+(gstar-1.)*bstarx/bstar**2*
     >		(bstary*cospq-bstarx*sinpq) ) +
     >		( (phadcmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dphi
	dpydphi = vertex%p%P*tsin*(cospq+(gstar-1.)*bstary/bstar**2*
     >		(bstary*cospq-bstarx*sinpq) ) +
     >		( (phadcmy+gstar*bstary*vertex%p%E)/vertex%p%P -
     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dphi
	dpzdphi =  vertex%p%P*(gstar-1.)/bstar**2*bstarz*tsin*
     >		(bstary*cospq-bstarx*sinpq) +
     > 		((phadcmz+gstar*bstarz*vertex%p%E)/vertex%p%P-
     >		gstar*bstarz*vertex%p%P/vertex%p%E)*dp_dphi

	dpxdcos = -vertex%p%P*tcos/tsin*(cospq+(gstar-1.)*bstarx/bstar**2*
     >		(bstarx*cospq+bstary*sinpq-bstarz*tsin/tcos)) +
     >		( (phadcmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dcos
	dpydcos = -vertex%p%P*tcos/tsin*(sinpq+(gstar-1.)*bstary/bstar**2*
     >		(bstarx*cospq+bstary*sinpq-bstarz*tsin/tcos)) +
     >		( (phadcmy+gstar*bstary*vertex%p%E)/vertex%p%P -
     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dcos
	dpzdcos = vertex%p%P*(1.-(gstar-1.)/bstar**2*bstarz*tcos/tsin*
     >		(bstarx*cospq+bstary*sinpq-tsin/tcos*bstarz))
     >		+((phadcmz+gstar*bstarz*vertex%p%E)/vertex%p%P-gstar*bstarz*
     >		vertex%p%P/vertex%p%E)*dp_dcos

	dpxnewdphi = dpxdphi*tmp2_x_x+dpydphi*tmp2_x_y+dpzdphi*tmp2_x_z
	dpynewdphi = dpxdphi*tmp2_y_x+dpydphi*tmp2_y_y+dpzdphi*tmp2_y_z
	dpznewdphi = dpxdphi*tmp2_z_x+dpydphi*tmp2_z_y+dpzdphi*tmp2_z_z

	dphicmdphi = (dpynewdphi*phadcm_tmp2x - phadcm_tmp2y*dpxnewdphi)/
     >			(phadcm_tmp2x**2+phadcm_tmp2y**2)

	dpxnewdcos = dpxdcos*tmp2_x_x+dpydcos*tmp2_x_y+dpzdcos*tmp2_x_z
	dpynewdcos = dpxdcos*tmp2_y_x+dpydcos*tmp2_y_y+dpzdcos*tmp2_y_z
	dpznewdcos = dpxdcos*tmp2_z_x+dpydcos*tmp2_z_y+dpzdcos*tmp2_z_z

	dphicmdcos = (dpynewdcos*phadcm_tmp2x - phadcm_tmp2y*dpxnewdcos)
     >			/(phadcm_tmp2x**2+phadcm_tmp2y**2)

	jacobian = abs(dt_dcos_lab*dphicmdphi-dt_dphi_lab*dphicmdcos)


* Old jacobian - assumes collinear boost (i.e. pfer along q vector).
	jac_old = 2*(efer-2*pferz*pfer*vertex%p%E/vertex%p%P*tcos)*
     >          (vertex%q+pferz*pfer)*vertex%p%P /
     >          ( efer+vertex%nu-(vertex%q+pferz*pfer)*vertex%p%E/vertex%p%P*tcos )
     >          - 2*vertex%p%P*pfer


	return
	end
