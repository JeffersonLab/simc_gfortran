	real*8 function peerho(vertex,main)

* Purpose:
* This routine calculates p(e,e'rho)p cross sections using my attempt
* to code in the form used in PYTHIA with modifications implmented
* in the HERMES Monte Carlo by Patricia Liebing (Thanks Patty!)
*
*
* variables:
*   input:
*	omega		!energy of virtual photon		(MeV)
*	theta_pq	!angle between pi and q			(rad)
*	Q2_g            !4-momentum of virtual photon, squared  (GeV/c)^2
*	epsilon		!epsilon				(dimensionless)
*
*   output:
*	sigma_eerho	!d3sigma/dEe'dOmegae'Omegapi	(microbarn/MeV/sr^2)

	implicit none
	include 'simulate.inc'

	type(event_main):: main
	type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

	real*8 sigma_eerho		!final cross section (returned as peepi)
	real*8 sig219,sig,factor
	real*8 sigt		        !components of dsigma/dt
	real*8 s5lab			!dsigma/dE_e/dOmega_e/dOmega_rho (in some lab)
c	real*8 efer			!energy of target particle
	real*8 epsi			!epsilon of virtual photon
	real*8 gtpr			!gamma_t prime.
	real*8 tcos,tsin		!cos/sin of theta between ppi and q
	real*8 tfcos,tfsin		!cos/sin of theta between pfermi and q
	real*8 bstar,bstarx,bstary,bstarz,gstar	!beta/gamma of C.M. system in lab
	real*8 ppix,ppiy,ppiz		!pion momentum in lab.
	real*8 thetacm,phicm		!C.M. scattering angles
	real*8 costcm,costo2cm
	real*8 sintcm,sinto2cm
	real*8 ppicm,ppicmx,ppicmy,ppicmz,epicm	 !pion E,p in C.M.
	real*8 qstar,qstarx,qstary,qstarz,nustar !c.m. photon momentum/energy
	real*8 zero

	real*8 s,t,Q2_g			!t,s, Q2 in (GeV/c)**2
	real*8 tmin, tprime, cdeltatau,brho
	real*8 nu			!equivilent photon energy (MeV)

	real*8 sig0,R

	real*8 new_x_x,new_x_y,new_x_z
	real*8 new_y_x,new_y_y,new_y_z
	real*8 new_z_x,new_z_y,new_z_z
	real*8 dummy,p_new_x,p_new_y,phiqn
	real*8 davesig,phipq,cospq,sinpq
	real*8 square_root,dp_dcos_num,dp_dcos_den,dp_dcos
	real*8 dp_dphi_num,dp_dphi_den,dp_dphi
	real*8 dt_dcos_lab,dt_dphi_lab,psign
	real*8 dpxdphi,dpydphi,dpxdcos,dpydcos,dpzdcos,dpzdphi
	real*8 dpxnewdphi,dpynewdphi,dpxnewdcos,dpynewdcos
	real*8 dpznewdphi,dpznewdcos
	real*8 dphicmdphi,dphicmdcos
	real*8 jacobian

	real*8 pbeam,beam_newx,beam_newy,beam_newz
	real*8 pbeamcmx,pbeamcmy,pbeamcmz,ebeamcm,pbeamcm
	real*8 ppicm_newx,ppicm_newy,ppicm_newz

*	real*8 dEcmdcos,dEcmdphi
*	real*8 dcoscmdcos,dcoscmdphi
	real*8 qx,qy,qz,px,py,pz


* Initialize some stuff.
	Q2_g = vertex%Q2/1000000.
	phipq = main%phi_pq
	cospq = cos(phipq)
	sinpq = sin(phipq)

* calculate energy of initial (struck) nucleon, using the same assumptions that
* go into calculating the pion angle/momentum (in event.f).  For A>1, the struck
* nucleon is off shell, the 2nd nucleon (always a neutron) is on shell, and has
* p = -p_fermi, and any additional nucleons are at rest.

	efer = sqrt(pfer**2+targ%Mtar_struck**2)
	if(doing_deutpi.or.doing_hepi) then
	  efer = targ%M-sqrt(Mn2+pfer**2)
	  if(doing_hepi)efer=efer-mp
c	  mtar_offshell = sqrt(efer**2-pfer**2)
	endif


* calculate some kinematical variables
* f's and fer indicate fermi momenta, s, star or cm CM system

	tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y+vertex%up%z*vertex%uq%z
	if(tcos-1..gt.0..and.tcos-1..lt.1.e-8)tcos=1.0
	tsin=sqrt(1.-tcos**2)

	tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0
	tfsin=sqrt(1.-tfcos**2)

	epsi = 1./(1.+2*(1.+vertex%nu**2/vertex%Q2)*tan(vertex%e%theta/2.)**2)

	s = (vertex%nu+efer)**2-(vertex%q+pfer*tfcos)**2-(pfer*tfsin)**2
	nu = (s-(efer**2-pfer**2))/2./(efer-pfer*tfcos) !equiv pho energy(MeV)

	s = s/1.e6				!CONVERT TO (GeV)**2
	main%w = sqrt(s)

	t = vertex%Q2-Mrho2+2.*vertex%nu*vertex%p%E-2.*vertex%p%P*vertex%q*tcos
	t = t/1.e6				!CONVERT TO (GeV/c)**2
	main%t = t


* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.
* DJG: Get pfer components in the lab "q" system

* JRA:The transformation is calculated starting in the coord. system used
* in the fpi/nucpi replay (see event.f), where x=right, y=down, z=along beam.
* We must convert from SIMC coords to these first.
	qx = -vertex%uq%y           !'right'
	qy =  vertex%uq%x           !'down'
	qz =  vertex%uq%z
	px = -pfery
	py =  pferx
	pz =  pferz

	dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	new_x_x = -qx*qz/dummy
	new_x_y = -qy*qz/dummy
	new_x_z = (qx**2 + qy**2)/dummy

	dummy   = sqrt(qx**2 + qy**2)
	new_y_x =  qy/dummy
	new_y_y = -qx/dummy
	new_y_z =  0.0

	p_new_x = pfer*(px*new_x_x + py*new_x_y + pz*new_x_z)
	p_new_y = pfer*(px*new_y_x + py*new_y_y + pz*new_y_z)

	if(p_new_x.eq.0.)then
	  phiqn=0.
	else
	  phiqn = atan2(p_new_y,p_new_x)
	endif
	if(phiqn.lt.0.)phiqn = phiqn+2.*pi

* get beam in "q" system.

	pbeam = sqrt(vertex%Ein**2-me**2)
	beam_newx = pbeam*new_x_z
	beam_newy = pbeam*new_y_z
	beam_newz = pbeam*vertex%uq%z

	bstar=sqrt((vertex%q+pfer*tfcos)**2+(pfer*tfsin)**2)/(efer+vertex%nu)
	gstar=1./sqrt(1. - bstar**2)

	bstarz = (vertex%q+pfer*tfcos)/(efer+vertex%nu)
	bstarx = p_new_x/(efer+vertex%nu)
	bstary = p_new_y/(efer+vertex%nu)

* DJG: Boost virtual photon to CM.

	zero =0.e0
	call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >		zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost pion to CM.
	
	ppiz = vertex%p%P*tcos
	ppix = vertex%p%P*tsin*cospq
	ppiy = vertex%p%P*tsin*sinpq
	call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >		ppix,ppiy,ppiz,epicm,ppicmx,ppicmy,ppicmz,ppicm)
	thetacm = acos((ppicmx*qstarx+ppicmy*qstary+ppicmz*qstarz)/ppicm/qstar)
	main%pcm = ppicm

* DJG Boost the beam to CM.

	call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_newx,
     >		beam_newy,beam_newz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* Thetacm is defined as angle between ppicm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

	dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >		(qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >		+(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
	new_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
	new_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
	new_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

	dummy = sqrt((new_y_y*qstarz-new_y_z*qstary)**2
     >		+(new_y_z*qstarx-new_y_x*qstarz)**2
     >		+(new_y_x*qstary-new_y_y*qstarx)**2)
	new_x_x = (new_y_y*qstarz-new_y_z*qstary)/dummy
	new_x_y = (new_y_z*qstarx-new_y_x*qstarz)/dummy
	new_x_z = (new_y_x*qstary-new_y_y*qstarx)/dummy

	new_z_x = qstarx/qstar
	new_z_y = qstary/qstar
	new_z_z = qstarz/qstar

	ppicm_newx = ppicmx*new_x_x + ppicmy*new_x_y + ppicmz*new_x_z
	ppicm_newy = ppicmx*new_y_x + ppicmy*new_y_y + ppicmz*new_y_z
	ppicm_newz = ppicmx*new_z_x + ppicmy*new_z_y + ppicmz*new_z_z

	sintcm = sin(thetacm)
	costcm = cos(thetacm)
	sinto2cm = sin(thetacm/2.)
	costo2cm = cos(thetacm/2.)
	phicm = atan2(ppicm_newy,ppicm_newx)
	if(phicm.lt.0.) phicm = 2.*3.141592654+phicm

	main%thetacm = thetacm
	main%phicm = phicm

! DJG need tprime and tmin here
! DJG Need an overall minus sign since we calculate -t above.

	tmin = -( ((-Q2_g-mrho2/1.e6-(targ%mtar_struck/1000.)**2+
     >       (targ%mtar_struck/1000.)**2)/(2.*sqrt(s)))**2-
     >       ((qstar-ppicm)/1000.)**2 )


	tprime = abs(t-tmin)
! Put in some W dependence from photoproduction data
! DJG:  This is my rough fit to some old photoproduction data.	

	sig0 = 41.263/(vertex%nu/1000.0)**0.4765   ! microbarns

! DJG:  R is usually fit to the form c_0 (Q2/M2_rho)^c1
! DJG:  The c_0 and c_1 are taken from HERMES data.
! DJG:  WARNING!!!! If you change this parameterization, you really
! DJG:  should change it in rho_decay also if you want to be self
! DJG:  consistent.

	R = 0.33*(vertex%Q2/mrho2)**(0.61)
! PYB added this
	if(R.lt.0.) write(222,'(''r'',3e12.3)')
     >    r,vertex%Q2,mrho2
	if(R.lt.0.) R=0.

! DJG:  The Q2 dependence is usually given by (M2_rho/(Q2+M2_rho))^2
! DJG:  HERMES found that 2.575 works better than 2 in the exponent.

	sigt = sig0*(1.0+epsi*R)*(mrho2/(vertex%Q2+mrho2))**(2.575)

! PYB added this
	if(sigt.lt.0.) write(222,'(''sigt'',6e12.3)')
     >    sigt, sig0, epsi,r,mrho2,vertex%q2
	if(sigt<0.) sigt = 0.


! DJG:  Need to parameterize t-dependence with b parameter as a function of c-tau

	cdeltatau = hbarc/(sqrt(vertex%nu**2+vertex%Q2+mrho2)-vertex%nu) !in fm!
	if(cdeltatau.lt.2.0) then
	   brho = 4.4679 + 8.6106*log10(cdeltatau)
! PYB added this
	   if(brho.lt.1.0) brho = 1.0
	else
	   brho = 7.0
	endif

	sig219 = sigt*brho*exp(-brho*tprime)/2.0/pi !ub/GeV**2/rad

	sig=sig219/1.d+06	!dsig/dtdphicm in microbarns/MeV**2/rad


C DJG Convert to dsig/dOmega_cm using dt/d(costhetacm) = 2 qcm pcm

	sig = sig*2.*qstar*ppicm

*******************************************************************************
* GMH: Virtual photon to electron beam flux conversion factor
	gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*(s-(targ%Mtar_struck/1000.)**2)/2./
     >		((efer-pfer*tfcos)/1000.)/Q2_g/(1.-epsi)

	if(gtpr.le.0.) then
         write(222,'(''gtpr'',10e10.2)') 
     >   gtpr, efer,pfer,
     >   vertex%Ein, s, targ%Mtar_struck, tfcos,Q2_g
	 gtpr = 0.
	endif

	factor=targ%Mtar_struck/(targ%Mtar_struck+vertex%nu-vertex%q*vertex%p%E/vertex%p%P*tcos)
	s5lab = gtpr*sig*vertex%q*vertex%p%P*factor/pi/1.d+06

* DJG The above assumes target at rest.  We need a full blown Jacobian:
* DJG | dt/dcos(theta) dphicm/dphi - dt/dphi dphicm/dcos(theta) |
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
     >		( (ppicmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dphi
	dpydphi = vertex%p%P*tsin*(cospq+(gstar-1.)*bstary/bstar**2*
     >		(bstary*cospq-bstarx*sinpq) ) +
     >		( (ppicmy+gstar*bstary*vertex%p%E)/vertex%p%P -
     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dphi
	dpzdphi =  vertex%p%P*(gstar-1.)/bstar**2*bstarz*tsin*
     >		(bstary*cospq-bstarx*sinpq) +
     > 		((ppicmz+gstar*bstarz*vertex%p%E)/vertex%p%P-
     >		gstar*bstarz*vertex%p%P/vertex%p%E)*dp_dphi

	dpxdcos = -vertex%p%P*tcos/tsin*(cospq+(gstar-1.)*bstarx/bstar**2*
     >		(bstarx*cospq+bstary*sinpq-bstarz*tsin/tcos)) +
     >		( (ppicmx+gstar*bstarx*vertex%p%E)/vertex%p%P -
     >		gstar*bstarx*vertex%p%P/vertex%p%E)*dp_dcos
	dpydcos = -vertex%p%P*tcos/tsin*(sinpq+(gstar-1.)*bstary/bstar**2*
     >		(bstarx*cospq+bstary*sinpq-bstarz*tsin/tcos)) +
     >		( (ppicmy+gstar*bstary*vertex%p%E)/vertex%p%P -
     >		gstar*bstary*vertex%p%P/vertex%p%E)*dp_dcos
	dpzdcos = vertex%p%P*(1.-(gstar-1.)/bstar**2*bstarz*tcos/tsin*
     >		(bstarx*cospq+bstary*sinpq-tsin/tcos*bstarz))
     >		+((ppicmz+gstar*bstarz*vertex%p%E)/vertex%p%P-gstar*bstarz*
     >		vertex%p%P/vertex%p%E)*dp_dcos

	dpxnewdphi = dpxdphi*new_x_x+dpydphi*new_x_y+dpzdphi*new_x_z
	dpynewdphi = dpxdphi*new_y_x+dpydphi*new_y_y+dpzdphi*new_y_z
	dpznewdphi = dpxdphi*new_z_x+dpydphi*new_z_y+dpzdphi*new_z_z

	dphicmdphi = (dpynewdphi*ppicm_newx - ppicm_newy*dpxnewdphi)/
     >			(ppicm_newx**2+ppicm_newy**2)

	dpxnewdcos = dpxdcos*new_x_x+dpydcos*new_x_y+dpzdcos*new_x_z
	dpynewdcos = dpxdcos*new_y_x+dpydcos*new_y_y+dpzdcos*new_y_z
	dpznewdcos = dpxdcos*new_z_x+dpydcos*new_z_y+dpzdcos*new_z_z

	dphicmdcos = (dpynewdcos*ppicm_newx - ppicm_newy*dpxnewdcos)
     >			/(ppicm_newx**2+ppicm_newy**2)

* JRA: dEcmdcos is never initialized, and these are never used.
*	dcoscmdcos = dpznewdcos/ppicm-ppicm_newz*epicm*dEcmdcos/ppicm**1.5
*	dcoscmdphi = dpznewdphi/ppicm-ppicm_newz*epicm*dEcmdphi/ppicm**1.5

	jacobian = abs(dt_dcos_lab*dphicmdphi-dt_dphi_lab*dphicmdcos)
	main%davejac = jacobian

	main%johnjac = 2*(efer-2*pferz*pfer*vertex%p%E/vertex%p%P*tcos)*
     >		(vertex%q+pferz*pfer)*vertex%p%P /
     >		( efer+vertex%nu-(vertex%q+pferz*pfer)*vertex%p%E/vertex%p%P*tcos )
     >		- 2*vertex%p%P*pfer

c	davesig = gtpr*sig*jacobian

	davesig = gtpr*sig
	sigma_eerho = davesig/1.d3	!ub/GeV-sr --> ub/MeV-sr
c avoid crazy results
c	write(222,'(20e10.2)') sigma_eerho,
c     >   gtpr,sig,sigt,brho,tprime,pfer,r
	if(sigma_eerho.gt.1.E10.or.
     >    sigma_eerho.lt.0.) then
	 write(222,'("crazy sigma_rho",8e12.4)') 
     >    gtpr, sig, sigma_eerho
	 sigma_eerho = 0.
	endif
	peerho = sigma_eerho
	ntup%sigcm = sig	!sig_cm

202	format(/11X,f5.1/)
203	format(11X,f5.0)
204	format(6(/9X,7f8.3))

	return
	end
