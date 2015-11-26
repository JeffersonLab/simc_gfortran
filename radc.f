!--------------------------------------------------------------------------

	subroutine basicrad(itail,Egamma_lo,Egamma_hi,Egamma,weight,val_reciprocal)

	implicit none
	include 'radc.inc'

	integer itail
	real*8 Egamma_lo, Egamma_hi, Egamma, weight, val_reciprocal
	real*8 power_lo, power_hi
	real*8 x, y, ymin

	real*8 grnd

!--------------------------------------------------------------------------
!
! Generate Egamma inside the requested region, using the function 
!	    g * Egamma**(g-1)
!	-------------------------------
!	(Egamma_max**g - Egamma_min**g)
! which has the BASICRAD shape, is "invertible" as a probability 
! distribution,
! and integrates to 1 over the requested region { Egamma_min to Egamma_max }.
! Each radiative event we generate with this then has to be weighted with
! the _actual_ radiative tail form (which doesn't necessarily integrate to 1)
! _divided_by_ this generating function, so we return the RECIPROCAL of the 
! function's value (at the selected Egamma) as VAL_RECIPROCAL. If we're not 
! using anything fancier than BASICRAD, we can immediately determine what
! that 
! quotient is:
!		  C
!	weight = --- * (Egamma_max**g - Egamma_min**g)
!		  g
! We return that as (you guessed it!) WEIGHT.
! If we want to use BASICRAD x the phi function (correction to external tail
! shape), all we need to do is multiply WEIGHT by PHI evaluated at the
! selected photon energy. That's dont outside this routine, in GENERATE_RAD.
! Note that if we're not using the BASICRAD prescription at all (except as
! a generation tool), WEIGHT is not used again. 
!
! Note this routine handles the possibilities Egamma_lo < 0, Egamma_hi <
! Egamma_lo, and Egamma_hi < 0
!---------------------------------------------------------------------------

! Initialize

	Egamma = 0.0
	weight = 0.0
	val_reciprocal = 0.0

! ... is radiation in this direction turned off?

	if(itail.eq.0)itail=4
	if (g(itail).le.0) then
	  weight = 1.0
	  return
	endif

! ... do we have a chance?

	if (Egamma_hi.le.Egamma_lo .or. Egamma_hi.le.0) return

! ... want to evaluate powers as INfrequently as possible!

	power_hi = Egamma_hi**g(itail)
	power_lo = 0.0
	if (Egamma_lo.gt.0) power_lo = Egamma_lo**g(itail)

! ... obtain y region: from (Egamma_lo/Egamma_hi)**g to 1, generate flat y.

	ymin = power_lo/power_hi
	y = ymin + grnd()*(1.-ymin)
	x = y**(1./g(itail))
	Egamma = x*Egamma_hi

! The value of our generating function at the selected point

	if (Egamma.gt.0) val_reciprocal = Egamma**(1.-g(itail)) *
     >		(power_hi-power_lo) / g(itail)

! Event weight = probability radiation falls inside requested range,
! in full BASICRAD prescription

	weight = c(itail)/g(itail) * (power_hi-power_lo)
	
	if(itail.eq.4)itail=0
	return
	end

!-------------------------------------------------------------------------

	real*8 function gamma(x)

	implicit none

	integer i, n, s
	real*8	x, y

! Compute gamma function for xin
! The series computes gamma(1+y), and must be fed y between 0 and 1
! This can be accomplished using the relation
! gamma(y+n+1) = (y+n)*gamma(y+n) = ... = (y+n)*...*(y+1)*gamma(y+1)

	gamma = 1.0
	n = nint((x-1)-0.5)
	y = x-1 - n
	if (n.ne.0) then
	  s = sign(1,n)
	  do i = s, n, s
	    gamma = gamma*(y+1+i)**s
	  enddo
	endif
	gamma = gamma * (1. - 0.5748646*y + 0.9512363*y**2 -
     >		 0.6998588*y**3 + 0.4245549*y**4 - 0.1010678*y**5)
	return
	end

!--------------------------------------------------------------------------

	subroutine generate_rad(main,vertex,orig,success)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'
	
	integer i, peaked_basis_flag
	real*8 x, rad_weight
	real*8 peaked_rad_weight, basicrad_val_reciprocal
	real*8 basicrad_weight, extrad_phi
	real*8 max_delta_Trec, ebeam_min, ebeam_max
	logical force2		!force tail #2 to radiate.
	logical evsuccess,success
	type(event_main):: main
	type(event):: vertex, orig

	real*8 grnd

!---Generate Radiation-----------------------------------------------------
! This routine radiates the incoming event
! Here are comments from radc_init that describe rad_flag and extrad_flag:
!
! First, about those (mysterious) 2 main 'radiative option' flags ...
!
! The significance of RAD_FLAG:
!        RAD_FLAG  = 0  .. use best available formulas, generate in
!                       .. (ntail,Egamma) basis
!                  = 1  .. use BASICRAD only, generate in (ntail,Egamma)
!                       .. basis
!                  = 2  .. use BASICRAD only, generate in (Egamma(1,2,3))
!                       .. basis but prevent overlap of tails (bogus, note)
!                  = 3  .. use BASICRAD only, generate in (Egamma(1,2,3))
!                       .. allowing radiation in all 3 directions
!                       .. simultaneously
! The (ntail,Egamma) basis can be called the PEAKED basis since it allows
! only 3 photon directions. PEAKED_BASIS_FLAG is set to zero below when
! the peaked basis is being used, in this way we can conveniently tell
! the BASICRAD routine to use the full Egamma formula to generate the gamma
! energy whenever it's called.
!
! (See N. Makins' thesis, section 4.5.6, for more help on this point)
!
! The significance of EXTRAD_FLAG:
!       EXTRAD_FLAG = 1 .. use only BASIC external radiation formulas
!                       .. (phi = 1)
!                   = 2 .. use BASIC ext rad formulas x phi
!                   = 3 .. use Friedrich approximation the way we always
!                       .. have
!                   = 0 .. use DEFAULTS: 3 for RAD_FLAG = 0, 1 otherwise; note
!                          that the defaults mimic the hardwired 'settings'
!                          in SIMULATE, which doesnt read EXTRAD_FLAG
!                          but determines what to do based on RAD_FLAG
!---------------------------------------------------------------------------

! Initialize.
! ... ntup.radphot=sum energy of all radiated photons.
! ... ntup.radarm=ntail, or 12 if tail 2 is forced.

	lim1=0
	lim2=0
	lim3=0
	phot1=0
	phot2=0
	phot3=0

	if (debug(2)) write(6,*)'gen_rad: entering...'
	success = .false.
	rad_weight = 1
	do i = 1,3
	  Egamma_used(i) = 0.0
	enddo
	ntup%radphot=0.
	ntup%radarm=0.
	force2=.false.

! Which tail has to be radiated? Set ntail=0 if all tails allowed to radiate.

	peaked_basis_flag = 1
	if (rad_flag.le.1) then
	  peaked_basis_flag = 0
	  x = grnd()
	  if (x.ge.frac(1)+frac(2)) then
	    ntail = 3
	  else if (x.ge.frac(1)) then
	    ntail = 2
	  else
	    ntail = 1
	  endif
	else if (rad_flag.eq.2) then
	  ntail = int(grnd()*3.) + 1
	  if (ntail.eq.4) ntail=3
	else if (rad_flag.eq.3) then
	  ntail = 0
	else
	  stop 'Idiot! rad_flag is set stupidly'
	endif
	if (debug(5)) write(6,*) 'vertexradc',vertex%e%E,edge%e%E%max

! If the scattered electron energy is too high to be dectected, only events
! where electron radiates will be accepted.  Force ntail=2, and apply extra
! weight to event equal to frac(2) (the normal probability that ntail=2)

cdg	if (vertex.e.E .gt. edge.e.E.max) then
cdg	  force2=.true.
cdg	  ntail=2
cdg	endif
! Need maximum change in recoil momentum for some of the later limits (e,e'p).

	max_delta_Trec = max((vertex%Trec - VERTEXedge%Trec%min),
     >				 (VERTEXedge%Trec%max - vertex%Trec) )

! RADIATE TAIL #1: the incident electron tail

	if (doing_tail(1) .and. (ntail.eq.0.or.ntail.eq.1)) then
	  if (debug(5)) write(6,*) 'gen_rad: at 1'

! ... CASE 1: A(e,e'p) for A>2 - the effect of radiation on this tail alters
! ... the VERTEX (actual) value of Em, so we want to make sure we're inside
! ... any source spectral function cuts (like the Fermi surface).
! ... Set Egamma_max to avoid getting below Em.min, and if Em is above Em.max,
! ... set minimum to get below Em.max.  Radiating a photon decreases Em
! ... (relative to pre-radiated value), but also shifts Pm, and thus Trec.
! ... Allow Trec shift as well as Em.  Note that there are NO limits based on
! ... spectrometer Em limits (cuts) because measured Em used Ebeam before rad.

! ... If only the incoming electron arm radiates, then we can compare Em
! ... to edge.Em.min (the SPECTROMETER Em limit) as well.

	  if (doing_heavy) then
	    if (debug(5)) write(6,*)'gen_rad: at 2a'
	if (debug(5)) write(6,*)'vertex%Em=',vertex%Em
	if (debug(5)) write(6,*)'VERTEXedge%Em%max=',VERTEXedge%Em%max
	if (debug(5)) write(6,*)'max_delta_Trec=',max_delta_Trec
	    Egamma_min(1)=vertex%Em - VERTEXedge%Em%max - max_delta_Trec
	    Egamma_max(1)=vertex%Em - VERTEXedge%Em%min + max_delta_Trec
	    if (ntail.ne.0) Egamma_max(1)=min(Egamma_max(1),
     >			vertex%Em - edge%Em%min + max_delta_Trec)

! ... CASE 2: H(e,e'p) - Require that Ebeam after radiation can generate
! ... electron inside of electron arm acceptance (Ein_max for Ee'_max, etc...).
! ... Also require MEASURED Em < edge.Em.max (Em_meas=egamma1+egamma2+egamma3)

! ... JRA: Note that we get ebeam max/min from E=mE'/(m+E'*(1-cos(theta)),
! ... assuming that E'_max(min) gives E_max(min).  However, for very large
! ... angles, or large E'_max, you can get a negative value for E_max.  When
! ... this happens, just open up ebeam_max.  In the long run, we need to get
! ... a true E_max, not just E for E'_max.

	  else if (doing_hyd_elast) then
	    if (debug(5)) write(6,*)'gen_rad: at 2b'
	    ebeam_max = Mp*edge%e%E%max / (Mp-edge%e%E%max*(1.-vertex%ue%z))
	    if (ebeam_max.lt.0) ebeam_max=1.d10
	    ebeam_min = Mp*edge%e%E%min / (Mp-edge%e%E%min*(1.-vertex%ue%z))
	    Egamma_min(1) = vertex%Ein - ebeam_max
	    Egamma_max(1) = vertex%Ein - ebeam_min
	    Egamma_max(1) = min(Egamma_max(1),edge%Em%max)

! ... CASE 3: D(e,e'p) - limit radiation so that E_beam > E' at the vertex.

	  else if (doing_deuterium) then
	    Egamma_max(1) = min(Egamma1_max, gen%sumEgen%max - vertex%e%E)
	    if (ntail.ne.0) Egamma_min(1) = gen%sumEgen%min - vertex%e%E

! ... CASE 4: Pion production.  Too complicated to calculate min/max ebeam
! ... from scattered particles.  Use sumEgen-vertex.e.E as limit for remaining
! ... generated energy available.

	  else if (doing_pion .or. doing_kaon .or. doing_rho .or. 
     >  	 doing_semi) then
	    if (debug(5)) write(6,*)'gen_rad: at 2d'
	    Egamma_min(1) = 0.
	    Egamma_max(1) = gen%sumEgen%max - vertex%e%E
CDJG	    if (ntail.ne.0) Egamma_min(1) = gen.sumEgen.min - vertex.e.E
	  endif

! ... Constrain Egamma<Egamma1_max(from radc_init), and add energy 'slop'
	  Egamma_max(1) = min(Egamma_max(1),Egamma1_max)
	  Egamma_min(1) = Egamma_min(1) - dE_edge_test
	  Egamma_max(1) = Egamma_max(1) + dE_edge_test

! ... override Egamma_max with limit from dbase (if entered)
	  lim1=egamma_max(1)
	  if (hardwired_rad) then
	    Egamma_max(1) = Egamma_gen_max
	    if (lim1.gt.Egamma_gen_max) write(6,*)
     >		'SIMC found an Egamma(1) limit above your hardwired value'
	  endif

! ... radiate 
	  if (debug(5)) write(6,*)'gen_rad: init Egamma_min,max =',
     >		Egamma_min(1),Egamma_max(1)
	  call basicrad(1*peaked_basis_flag,Egamma_min(1),Egamma_max(1),
     >		Egamma_used(1),basicrad_weight,basicrad_val_reciprocal)
	  if (basicrad_weight.le.0) then
	    if (debug(4)) write(6,*)'gen_rad: failed at 4'
	     return
	  endif
	  if (debug(5)) write(6,*)'Egamma_used',Egamma_used(1)

! ... adjust kinematics and call complete_ev.  If complete_ev fails, return.

	  vertex%Ein = vertex%Ein - Egamma_used(1)
	  if (debug(5)) write(6,*) 'calling complete_ev from radc'
	  call complete_ev(main,vertex,evsuccess)
	  if (.not.evsuccess) then	!remove previous radiation, try again.
!!	    if (ntried.le.5000) write(6,'(1x,''COMPLETE_EV1 failed in GENERATE_RAD at event'',i7,''!'')') nevent
	    return
	  endif
	  if (debug(4)) write(6,*)'gen_rad: at 5'
	  rad_weight = rad_weight * basicrad_weight
	  if (rad_flag.le.1) then
	    rad_weight = peaked_rad_weight(vertex,Egamma_used(1),Egamma_min(1),
     >		Egamma_max(1),basicrad_val_reciprocal,basicrad_weight)
	  else
	    rad_weight=rad_weight*extrad_phi(1,vertex%Ein,vertex%e%E,Egamma_used(1))
	  endif
	endif		!end of tail 1 (incoming electron)

! The VERTEX kinematics are now set, so check that we are inside the
! VERTEXedges (SF limits).

	if (doing_heavy) then
	  if (debug(5)) write(6,*)'gen_rad: Em, min, max =',vertex%Em,
     >		VERTEXedge%Em%min,VERTEXedge%Em%max
	  if (debug(5)) write(6,*)'gen_rad: Pm, min, max =',vertex%Pm,
     >		VERTEXedge%Pm%min,VERTEXedge%Pm%max
	  if (vertex%Em .lt. VERTEXedge%Em%min .or.
     >	      vertex%Em .gt. VERTEXedge%Em%max .or.
     >	      vertex%Pm .lt. VERTEXedge%Pm%min .or.
     >	      vertex%Pm .gt. VERTEXedge%Pm%max) then
	    if (debug(4)) write(6,*)'gen_rad: failed at 6'
	    return
	  endif
	endif

! RADIATE TAIL #2: the scattered electron tail

	if (doing_tail(2) .and. (ntail.eq.0.or.ntail.eq.2)) then
	  if (debug(5)) write(6,*) 'we are at 2'

! ... Find min/max Egammas that will make in into electron acceptance

	  Egamma_min(2) = vertex%e%E - edge%e%E%max
	  Egamma_max(2) = vertex%e%E - edge%e%E%min

! ... Measured Em must be within acceptance after radiation.  Gives upper
! ... limit on sum of photon energies, and lower limit if hadron arm does
! ... not radiate.  Em_MEAS = Em_VERT+(Egamma1+Egamma2+Egamma3)+/-delta_trec.

	  if (doing_eep) then
	    Egamma_max(2) = min(Egamma_max(2),
     >		(edge%Em%max - vertex%Em) - Egamma_used(1) + max_delta_Trec )
	    if (ntail.ne.0 .or. .not.rad_proton_this_ev)
     >		Egamma_min(2) = max(Egamma_min(2),
     >		(edge%Em%min - vertex%Em) - Egamma_used(1) - max_delta_Trec )
	  endif

	  Egamma_max(2) = min(Egamma_max(2),Egamma_tot_max-Egamma_used(1))
	  Egamma_min(2) = Egamma_min(2) - dE_edge_test
	  Egamma_max(2) = Egamma_max(2) + dE_edge_test

! ... override Egamma_max with limit from dbase (if entered)
	  lim2=egamma_max(2)
	  if (hardwired_rad) then
	    Egamma_max(2) = Egamma_gen_max
	    if (lim2.gt.Egamma_gen_max) write(6,*)
     >		'SIMC found an Egamma(2) limit above your hardwired value'
	  endif

! ... radiate
	  call basicrad(2*peaked_basis_flag,Egamma_min(2),Egamma_max(2),
     >		Egamma_used(2),basicrad_weight,basicrad_val_reciprocal)
	  if (basicrad_weight.le.0) then
!	    write(6,*)'gen_rad: Maybe you need to change Em limits!'
!	    write(6,*)'min/max(2)=',egamma_min(2),egamma_max(2)
	    if (debug(4)) write(6,*)'gen_rad: failed at 7'
	    return
	  endif
	  if(debug(5))write(6,*)'gen_rad: peaked_basis_flag=',peaked_basis_flag
	  if(debug(5))write(6,*)'gen_rad: Egamma min,max,used=',
     >		Egamma_min(2),Egamma_max(2),Egamma_used(2)
	  if(debug(5))write(6,*)'gen_rad: basicrad_weight=',basicrad_weight
	  rad_weight = rad_weight * basicrad_weight
	  if (rad_flag.le.1) then
	    rad_weight = peaked_rad_weight(vertex,Egamma_used(2),Egamma_min(2),
     >		Egamma_max(2),basicrad_val_reciprocal,basicrad_weight)
	  else 
	    rad_weight=rad_weight*extrad_phi(2,vertex%Ein,vertex%e%E,Egamma_used(2))
	  endif
	endif		!end of tail 2 (outgoing electron)

! RADIATE TAIL #3: the scattered proton tail

	if (rad_proton_this_ev .and. (ntail.eq.0.or.ntail.eq.3)) then
	  if (debug(5)) write(6,*) 'we are at 3'

! ... Find min/max Egammas that will make in into hadron acceptance

	  Egamma_min(3) = vertex%p%E - edge%p%E%max
	  Egamma_max(3) = vertex%p%E - edge%p%E%min

! ... Measured Em must be within acceptance after radiation.  Gives upper
! ... and lower limit on sum of photon energies (and thus on hadron radiation,
! ... since other arms are done).
! ... Em_MEAS = Em_VERT+(Egamma1+Egamma2+Egamma3)+/-delta_trec.

	  if (doing_eep) then
	    Egamma_max(3) = min(Egamma_max(3), (edge%Em%max - vertex%Em) -
     >		Egamma_used(1) - Egamma_used(2) + max_delta_Trec)
	    Egamma_min(3) = max(Egamma_min(3), (edge%Em%min - vertex%Em) -
     >		Egamma_used(1) - Egamma_used(2) - max_delta_Trec )
	  endif

	  Egamma_max(3) = min(Egamma_max(3),
     >		Egamma_tot_max-Egamma_used(1)-Egamma_used(2) )
	  Egamma_min(3) = Egamma_min(3) - dE_edge_test
	  Egamma_max(3) = Egamma_max(3) + dE_edge_test

! ... override Egamma_max with limit from dbase (if entered)
	  lim3=egamma_max(3)
	  if (hardwired_rad) then
	    Egamma_max(3) = Egamma_gen_max
	    if (lim3.gt.Egamma_gen_max) write(6,*)
     >		'SIMC found an Egamma(3	) limit above your hardwired value'
	  endif

! ... radiate
	  call basicrad(3*peaked_basis_flag,Egamma_min(3),Egamma_max(3),
     >		Egamma_used(3),basicrad_weight,basicrad_val_reciprocal)
	  if (basicrad_weight.le.0) then
	    if (debug(4)) write(6,*)'gen_rad: failed at 8'
	    return
	  endif
	  rad_weight = rad_weight * basicrad_weight
	  if (rad_flag.le.1) then
	    rad_weight = peaked_rad_weight(vertex,Egamma_used(3),Egamma_min(3),
     >		Egamma_max(3),basicrad_val_reciprocal,basicrad_weight)
	  else
	    rad_weight=rad_weight*extrad_phi(3,vertex%Ein,vertex%e%E,Egamma_used(3))
	  endif
	endif		!end of tail 3 (outgoing hadron)

! Now obtain description of ORIG (orig = vertex + radiation) event.
! ... Note that no kinematics have been changed yet, EXCEPT vertex.Ein, which
! ... has been reduced by Egamma_used(1).  orig.Ein is the initial (unradiated)
! ... beam energy, so add back Egamma_used(1).  Reduce orig.e.E and orig.p.E
! ... here, and they will be the values passed to the monte carlo (for mult.
! ... scattering and monte carlo model).

! ... Remember that we have taken Ein(GEN)-Ebeam_vertex_ave into
! ... account by shifting edge.Em)

c	do i = 1, neventfields
c	  orig.all(i) = vertex.all(i)
c	enddo
	orig = vertex

	if (debug(4)) write(6,*)'gen_rad: at 10: orig.p.E =',orig%p%E

! ... remove radiation from incoming electron.

	orig%Ein = vertex%Ein + Egamma_used(1)

! ... adjust the e'/hadron momenta (if they radiated significantly)

	orig%e%E = vertex%e%E - Egamma_used(2)
	if(orig%e%E.le.0e0) then
	   if (debug(4)) write(6,*)'gen_rad: Negative electron energy -failed'
	  return
	endif
	orig%e%P = orig%e%E
	orig%e%delta = (orig%e%P-spec%e%P)/spec%e%P*100.

	orig%p%E = vertex%p%E - Egamma_used(3)
	if(orig%p%E.le.Mh) then
	  if (debug(4)) write(6,*)'gen_rad: failed at 11'
	  return
	endif
	orig%p%P = sqrt(orig%p%E**2-Mh2)
	orig%p%delta = (orig%p%P-spec%p%P)/spec%p%P*100.

! ... Record some information for other routines.
	phot1=Egamma_used(1)
	phot2=Egamma_used(2)
	phot3=Egamma_used(3)
	ntup%radphot = phot1 + phot2 + phot3
	ntup%radarm=ntail
	if (force2) ntup%radarm=12

! Complete determination of the event weight due to radiation.

	main%gen_weight = main%gen_weight*rad_weight/hardcorfac
	if (force2) main%gen_weight=main%gen_weight*frac(2)
	if (debug(5)) write(6,*) 'mgw',main%gen_weight,rad_weight,hardcorfac
	success = .true.

	if (debug(2)) write(6,*)'gen_rad: exiting...'
	return
	end

!--------------------------------------------------------------------------

	real*8 function peaked_rad_weight(vertex,Egamma,
     >		emin,emax,basicrad_val_reciprocal,basicrad_weight)

	implicit none
	include 'radc.inc'
	include 'structures.inc'

	real*8		Egamma, basicrad_val_reciprocal, basicrad_weight
	real*8		t1, t2, r, phi_ext, ein, eout, emin, emax
	real*8		brem, bremos, extrad_phi, gamma
	real*8		dhard, dsoft_intmin, dsoft_intmax
	real*8          dsoft_ext1, dsoft_ext2, dsoft_extmin
	real*8		dsoft_int_primemin, dsoft_int_primemax
	real*8          dsoft_ext1_prime, dsoft_ext2_prime
	real*8		dsoft_ext_prime, dsoft_extmax, eul
	real*8          dsoft_int,dsoft_ext,dsoft_int_prime
	type(event):: vertex

	real*8 zero
	parameter (zero=0.0e0)	!double precision zero for subroutine calls

	basicrad_val_reciprocal=basicrad_val_reciprocal+0. !avoid unused variable error
! Compute a more precise value for the radiative probability at the
! selected Egamma, more precise than the BASICRAD distribution used to 
! generate Egamma that is.
!
! NB: This subroutine only deals with calculations done in the PEAKING
! APPROXIMATION basis -
!    i.e. we only get here if RAD_FLAG = 0 or 1
!
! ... (sort of) KLUGE -- if Egamma < res limit, these things might blow up,
! ... so just screw the correction and leave
	ein=vertex%Ein
	eout=vertex%e%E
	eul=0.577215665
c	if (Egamma.lt.Egamma_res_limit) then
c	  peaked_rad_weight = basicrad_weight
c	  return
c	endif

! ... External

	phi_ext = 1.0
	if (extrad_flag.le.2) then
	  phi_ext = extrad_phi(0,ein,eout,Egamma)
	  if (rad_flag.eq.1) then
	    peaked_rad_weight = basicrad_weight * phi_ext
	    return
	  endif
	  if(emin.gt.0)then
	    dsoft_extmin = log(g_ext/c_ext(0)/emin**g_ext)
	  endif
	  dsoft_extmax = log(g_ext/c_ext(0)/emax**g_ext)
	  dsoft_ext_prime = -g_ext/Egamma
	else
	  t1 = bt(1)/etatzai
	  t2 = bt(2)/etatzai
	  call extrad_friedrich(ein,Egamma,t1,dsoft_ext1,dsoft_ext1_prime)
	  call extrad_friedrich(eout,Egamma,t2,dsoft_ext2,dsoft_ext2_prime)
	  dsoft_ext = dsoft_ext1 + dsoft_ext2
	  dsoft_ext_prime = dsoft_ext1_prime + dsoft_ext2_prime
	endif

! ... Internal
! ........ use full BREM calculation of deltas

	if (rad_flag.eq.0) then
	  if (.not.use_offshell_rad) then
	    if(emin.gt.0)then
	      r = brem(ein,eout,emin,rad_proton_this_ev,
     >			dsoft_intmin,dhard,dsoft_int_primemin)
	    else
	      dsoft_intmin=1.0
	    endif
	    r = brem(ein,eout,emax,rad_proton_this_ev,
     >			dsoft_intmax,dhard,dsoft_int_primemax)
	  else
	    if(emin.gt.0)then
	      r = bremos(emin, zero, zero, ein,
     >			vertex%e%E*vertex%ue%x, vertex%e%E*vertex%ue%y,
     >			vertex%e%E*vertex%ue%z,
c     >			vertex%Pmx, vertex%Pmy, vertex%Pmz,
     >			zero,zero,zero,
     >			vertex%p%P*vertex%up%x, vertex%p%P*vertex%up%y,
     >			vertex%p%P*vertex%up%z, vertex%p%E,
     >			rad_proton_this_ev,
     >			dsoft_intmin, dhard, dsoft_int_primemin)
	    else 
	      dsoft_intmin=1.0
	    endif
	      r = bremos(emax, zero, zero, ein,
     >			vertex%e%E*vertex%ue%x, vertex%e%E*vertex%ue%y,
     >			vertex%e%E*vertex%ue%z,
c     >			vertex%Pmx, vertex%Pmy, vertex%Pmz
     >			zero,zero,zero,
     >			vertex%p%P*vertex%up%x, vertex%p%P*vertex%up%y,
     >			vertex%p%P*vertex%up%z, vertex%p%E,
     >			rad_proton_this_ev,
     >			dsoft_intmax, dhard, dsoft_int_primemax)
	  endif

! ........ use basic calculation of internal deltas

	else
	  dsoft_int = log(g_int/c_int(0)/Egamma**g_int)
	  dsoft_int_prime = -g_int/Egamma
	endif

! ... All together now

	if(emin.gt.0)then
	  peaked_rad_weight = c_ext(0)/g_ext*(exp(-dsoft_intmax)*
     >		emax**g_ext-exp(-dsoft_intmin)*emin**g_ext)
	else
	  peaked_rad_weight = c_ext(0)/g_ext*(exp(-dsoft_intmax)*emax**g_ext)
	endif
	peaked_rad_weight = peaked_rad_weight * exp(-eul*g(4))/gamma(1.+g(4))
     >		* gamma(1.+g(4)-bt(1)-bt(2))*gamma(1.+bt(1))
     >		* gamma(1.+bt(2))/gamma(1.+g(4))
	if (peaked_rad_weight.lt.0) peaked_rad_weight = 0

	return
	end

!---------------------------------------------------------------------------

	subroutine extrad_friedrich(Ei,Ecutoff,trad,dbrem,dbrem_prime)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	real*8	Ei, Ecutoff, trad, x, dbrem, dbrem_prime

	x = Ecutoff/Ei
	dbrem = trad*(-(etatzai-0.5) - etatzai*log(x) + etatzai*x - 0.5*x**2)
	dbrem_prime = -trad/Ei * (etatzai/x - etatzai + x)

	return
	end

!--------------------------------------------------------------------------

	real*8 function extrad_phi(itail,E1,E2,Egamma)

	implicit none
	include 'radc.inc'

	integer itail
	real*8 E1, E2, Egamma, E(2), x, t, gamma

	E(1) = E1
	E(2) = E2

! Compute the multiplicative external correction function PHI

! ... CASE 1: phi = 1

	extrad_phi = 1.0

! ... CASE 2: phi correctly computed according to Dave's formulas

	if (extrad_flag.eq.2) then
	  if (itail.eq.0) then
	    extrad_phi = 1. - (bt(1)/E(1)+bt(2)/E(2)) / (g(1)+g(2)) * Egamma
	  else if (itail.eq.1.or.itail.eq.2) then
	    extrad_phi = 1. - bt(itail)/E(itail) / g(itail) * Egamma
	  endif

! ... CASE 3: phi computed in Friedrich prescription

	else if (extrad_flag.eq.3) then
	  if (itail.eq.0) stop 'Idiot! a multiplicative factor EXTRAD_PHI is not defined for peaking approx and EXTRAD_FLAG>2!'
	  if (itail.eq.1.or.itail.eq.2) then
	    x = Egamma/E(itail)
	    t = bt(itail)/etatzai
	    extrad_phi = extrad_phi * (1. - x + x**2/etatzai) *
     >		exp(t*((etatzai-0.5)-etatzai*x+x**2/2.)) * gamma(1.+bt(itail))
	  endif
	endif

	return
	end

!------------------------------------------------------------------------

	real*8 function schwinger(Ecutoff,vertex,include_hard,dsoft,dhard)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	real*8		Ecutoff,lq,s2,b,spence,dsoft,dhard,spen
	logical		include_hard
	type(event):: vertex


	lq = log(vertex%Q2/Me**2)-1.0
	s2 = sin(vertex%e%theta/2.)**2
	b = 1.+2.*vertex%nu*s2/(targ%A*amu)
	spence = spen(1.-s2)-2.5893784
	dsoft = alpi*lq*log( vertex%Ein/(etta**2)*vertex%e%E*b/(Ecutoff**2) )
	dhard = -alpi*(2.166666*lq + spence - log(vertex%Ein/vertex%e%E)**2/2.0)

	if (use_expon.eq.0) then
	  schwinger = exp(dsoft)
	else
	  schwinger = 1.+dsoft
	endif

	if (include_hard) schwinger = schwinger/(1.-dhard)

	return
	end

!--------------------------------------------------------------------------

	real*8 function spen(x)

	implicit none

	integer	i
	real*8	x,s,y

! Approximation formula for the Spence-function or -integral according to
! abramowitz and stegun : formula 27.7.2

	y = 1.0
	s = 0.0
	i = 0
	do while (i.le.100 .and. abs(y).gt.abs(s)*1.e-4)
	  i = i+1
	  y = x*y
	  s = s+y/float(i**2)
	enddo
	spen = s

	return
	end

!--------------------------------------------------------------------------

	real*8 function lambda_dave(itail,plus_flag,doing_proton,e1,e2,e3,p3,th)

	implicit none
	include 'constants.inc'

	integer		itail,plus_flag
	real*8		e1,e2,e3,p3,th
	real*8		plus_term
	logical		doing_proton, warned/.false./

! initialize to silence gfortran error
	lambda_dave=0.0

! The extended peaking approximation

	plus_term = 0.0
	if (plus_flag.eq.1 .and. itail.lt.3) then
	  plus_term = log((1.-cos(th))/2.)

! ... only add in term due to ep interference if we're using proton
! ... radiation

	  if (doing_proton) plus_term = plus_term + 2.*log(e1/e2)
	endif

! Compute lambdas

	if (itail.eq.1) then
	  lambda_dave	= alpi*(2.*log(2.*e1/Me) -1. + plus_term)
	else if (itail.eq.2) then  
	  lambda_dave	= alpi*(2.*log(2.*e2/Me) -1. + plus_term)
	else if (itail.eq.3) then
	  if (doing_proton) then
	    lambda_dave	= alpi*((e3/p3)*log((e3+p3)/(e3-p3)) - 2.)
! ... at low energies, this may fritz, prevent human from causing itself
! ... too much damage

	    if (lambda_dave.lt.0) then
	      lambda_dave = 0.0
	      if (.not.warned) then
	        write(6,'(1x,a)') '+-----------------------------------+'
	        write(6,'(1x,a)') '| WARNING: lambda_dave(3)<0, think  |'
	        write(6,'(1x,a)') '| about running with proton rad OFF |'
	        write(6,'(1x,a)') '+-----------------------------------+'
	        warned = .true.
	      endif
	    endif
	  else
	    lambda_dave = 0.0
	  endif
	endif

	return
	end
