	real*8 function peedelta(vertex,main)

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

	real*8 sig_multipole,sig_blok

* JRA: There are all unused (obsolte, presumably)
*	integer final_state
*	logical first,low_w_flag
*
*	data first /.TRUE./
*	data low_w_flag /.FALSE./	!Assume high W kinematics to start

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
	if(tfcos-1..gt.0..and.tfcos-1..lt.1.d-8)tfcos=1.0
	tfsin=sqrt(1.-tfcos**2)

	s = (vertex%nu+efer)**2-(vertex%q+pfer*tfcos)**2-(pfer*tfsin)**2
	main%wcm = sqrt(s)



*******************************************************************************
* Get photon flux factor (two options, see comments below).
*
* DJG,2000: Replace targ.Mtar_struck in denominator of gammaflux with more 
* general efer-pfer*tfcos, for pfer =0 this reverts to old form
*	k_eq = (s-targ.Mtar_struck**2)/2./(efer-pfer*tfcos)
*
* JRA,2001: Go back to original version - more consistent with phase space used
* in the subroutine (according to DJG - see gaskell_model.ps)
	k_eq = (main%wcm**2-targ%Mtar_struck**2)/2./targ%Mtar_struck


* Get cross section in photon-nucleon center of mass.  sigcm1 is blok
* fit (default model - can always be evaluated).  sigcm2 is multipole,
* IF low_w_flag is set.
* NOTE: s, t, mtar, and Q2 must be converted to GeV first.

	ntup%sigcm1 = sig_blok(thetacm,phicm,main%t/1.d6,vertex%q2/1.d6,s/1.d6,main%epsilon,
     >		targ%Mtar_struck/1000.,which_pion)
	sigma_eepi = ntup%sigcm1

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

!	peedelta = sigma_eepi*jacobian*(gtpr*fac) !ub/MeV^2/rad-->ub/sr-->ub/MeV/sr

	peedelta = 1.0*jacobian*(gtpr*fac) !ub/MeV^2/rad-->ub/sr-->ub/MeV/sr

	return
	end
