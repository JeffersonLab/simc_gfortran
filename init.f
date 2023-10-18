	subroutine target_init(using_Eloss,using_Coulomb,the_cent,thp_cent,
     >              Pp_cent,Mh2,ebeam,Pe_cent)

	implicit none
	include 'constants.inc'
	include 'target.inc'

	integer	i
	real*8 the_cent, thp_cent, Pp_cent, Mh2, ebeam, energy
	real*8 Pe_cent
	real*8 za2, fc
	logical	using_Eloss, using_Coulomb

	real*8 zero
	parameter (zero=0.0e0)	!double precision zero for subroutines calls.

! The radiation length of the target

	targ%L1 = log(184.15) - log(targ%Z)/3.0
	targ%L2 = log(1194.) - 2.*log(targ%Z)/3.0
	if(targ%Z.eq.1)then
	  targ%L1=5.31
	  targ%L2=6.144
	endif
	za2 = (targ%Z*alpha)**2
	fc = za2*(1.202+za2*(-1.0369+za2*1.008/(za2+1)))

! ... The radiation length, in both g/cm2 and cm; For H,D,4He, using
! ... PDB values (1999), other calculated directly from Tsai formula
! ... (For 4He, the Tsai formula is 10% lower).

	if (nint(targ%A).eq.1) then
	  targ%X0 = 61.28
	else if (nint(targ%A).eq.2) then
	  targ%X0 = 122.4
	else if (nint(targ%A).eq.4) then
	  targ%X0 = 94.32
	else if (nint(targ%A).eq.3) then
	  write(6,*) 'Using Tsai formula for He-3 radiation length!!! (is there a better value?)'
	  targ%X0 = 716.405*targ%A/targ%Z/(targ%Z*(targ%L1-fc)+targ%L2)
	else
	  targ%X0 = 716.405*targ%A/targ%Z/(targ%Z*(targ%L1-fc)+targ%L2)
	endif
	targ%X0_cm = targ%X0/targ%rho

! ... 'Average' ionization losses (MOST PROBABLE, REALLY).  Printed out by simc
! ... (so ave or m.p. are both OK), and Ebeam_vertex_ave (which should be m.p.)

	energy = ebeam
	call trip_thru_target (1,zero,energy,zero,targ%Eloss(1)%ave,
     >                         targ%teff(1)%ave,Me,4)
	energy = Pe_cent
	call trip_thru_target (2,zero,energy,the_cent,targ%Eloss(2)%ave,
     >                         targ%teff(2)%ave,Me,4)
	energy = sqrt(Pp_cent**2 + Mh2)
	call trip_thru_target (3,zero,energy,thp_cent,targ%Eloss(3)%ave,
     >                         targ%teff(3)%ave,sqrt(Mh2),4)
	if (.not.using_Eloss) then
	  do i = 1, 3
	    targ%Eloss(i)%ave = 0.0
	  enddo
	endif

! ... Coulomb potential energy (NB All of the following are positive)

	if (using_Coulomb) then
c	  targ%Coulomb%ave=6./5.*(targ%Z-1.)*alpha*hbarc/(1.18*targ%A**(1./3.))
c	  targ%Coulomb_constant = 5./12. * targ%Coulomb%ave
c	  targ%Coulomb%min = targ%Coulomb_constant * 2.0
c	  targ%Coulomb%max = targ%Coulomb_constant * 3.0
*Next four lines were modified 5/15/06 for pionct (see Aste et al. Eur. Phys. J. A 26, 167 (2005)
*V(r) = -3/2 alpha Z/(2R) + alpha Z/(2R) (r/R)**2 where
* R = 1.1*A**1/3 + 0.86A**-1/3 fm
* Aste sez use V(0) * 0.75
	  targ%Coulomb%ave=0.75*1.5*(targ%Z-1.)*alpha*hbarc/(1.1*targ%A**(1./3.)+0.86*targ%A**(-1./3.))
	  targ%Coulomb_constant = targ%Coulomb%ave
	  targ%Coulomb%min = targ%Coulomb_constant 
	  targ%Coulomb%max = targ%Coulomb_constant 
	else
	  targ%Coulomb%ave = 0.0
	  targ%Coulomb_constant = 0.0
	  targ%Coulomb%min = 0.0
	  targ%Coulomb%max = 0.0
	endif

	return
	end

!----------------------------------------------------------------------

	subroutine limits_init(H)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'
	include 'histograms.inc'
	include 'sf_lookup.inc'		!need to know Em range of spec. fcn.

	integer i
	real*8 r
	real*8 Ebeam_min, Ebeam_max
	real*8 t1,t2			!temp. variables.
	real*8 slop_Coulomb, slop_Ebeam, slop_Ee, slop_Ep
	type(cutstype):: the_phys, thp_phys, z, pp
	type(histograms):: H

!-----------------------------------------------------------------------
! RECORDS store information relating to various quantities at several
! 'STATES' in event description:
! (of course, not all qties are _recorded_ in each of these stages)
!
!	orig.*		qties are TRUE qties but before ANY radiation
!	vertex.*	qties are those used to determine spec fn and sigcc
!			weighting, so TRUE qties before tails 2,3 radiated
!	recon.*		qties are defined AT the spectrometers, AFTER being
!			run through the single arm montecarlos
!
!	main.SP.*	qties are defined AT the spectrometers, so TRUE
!			qties AFTER any modifications due to non-radiative
!			interaction with the target (this includes Coulomb
!			accel/decel of initial/final electron)
!	main.FP.*	spect. focal plane values-if using spect. model
!	main.RECON.*	qties are defined AT the spectrometers, AFTER being
!			run through the single arm montecarlos
!
! Note that in the absence of radiation, GEN=VERTEX
! Note that if a spectrometer is not used, RECON=SP, FP is empty.
! Note that the generated qties we have to compare with gen limits later
! are unaffected by radiation in tail1, so equal to vertex qties.
!
!
! EDGES on event qties at any stage can be derived from these cuts, by taking
! the transformations (and associated uncertainties) between each stage into
! account. The general usefulness of defining these edges is to improve
! efficiency --> provide checks at intermediate stages to allow us to abort
! an event early, and provide some foreknowledge in generation routines of
! what will make it to the end. We must be certain that these derived EDGES
! are wide enough that the events making it through the cuts don't come
! uncomfortably close to them.
!
!	SPedge.*	Limits at spectrometer = Acceptance + slop
!	edge.*		Limits after interaction = SPedge + musc or eloss
!						 = VERTEXedge + coulomb + rad.
!	VERTEXedge.*	Limits at vertex (from energy conservation).
!	gen.*		Generation limits.
!
! SLOPS
!	slop.MC.*	Due to spectrometer resolution (spectrometer MC).
!	slop_*		Slop in calculated physics quantities (from spect.
!			resolution + range of eloss/coulomb/...
!
!----------------------------------------------------------------------

! ... get slop values for the difference between orig and recon
! ... SPECTROMETER quantities (recon inaccuracies due to the montecarlos)

	if (using_E_arm_montecarlo) then
	  if (electron_arm.eq.1) then
	    slop%MC%e%delta%used = slop_param_d_HMS
	    slop%MC%e%yptar%used = slop_param_t_HMS
	    slop%MC%e%xptar%used = slop_param_p_HMS
	  else if (electron_arm.eq.2) then
	    slop%MC%e%delta%used = slop_param_d_SOS
	    slop%MC%e%yptar%used = slop_param_t_SOS
	    slop%MC%e%xptar%used = slop_param_p_SOS
	  else if (electron_arm.eq.3) then
	    slop%MC%e%delta%used = slop_param_d_HRSR
	    slop%MC%e%yptar%used = slop_param_t_HRSR
	    slop%MC%e%xptar%used = slop_param_p_HRSR
	  else if (electron_arm.eq.4) then
	    slop%MC%e%delta%used = slop_param_d_HRSL
	    slop%MC%e%yptar%used = slop_param_t_HRSL
	    slop%MC%e%xptar%used = slop_param_p_HRSL
	  else if (electron_arm.eq.5 .or. electron_arm.eq.6) then
	    slop%MC%e%delta%used = slop_param_d_SHMS
	    slop%MC%e%yptar%used = slop_param_t_SHMS
	    slop%MC%e%xptar%used = slop_param_p_SHMS
	  endif
	endif
	if (using_P_arm_montecarlo) then
	  if (hadron_arm.eq.1) then
	    slop%MC%p%delta%used = slop_param_d_HMS
	    slop%MC%p%yptar%used = slop_param_t_HMS
	    slop%MC%p%xptar%used = slop_param_p_HMS
	  else if (hadron_arm.eq.2) then
	    slop%MC%p%delta%used = slop_param_d_SOS
	    slop%MC%p%yptar%used = slop_param_t_SOS
	    slop%MC%p%xptar%used = slop_param_p_SOS
	  else if (hadron_arm.eq.3) then
	    slop%MC%p%delta%used = slop_param_d_HRSR
	    slop%MC%p%yptar%used = slop_param_t_HRSR
	    slop%MC%p%xptar%used = slop_param_p_HRSR
	  else if (hadron_arm.eq.4) then
	    slop%MC%p%delta%used = slop_param_d_HRSL
	    slop%MC%p%yptar%used = slop_param_t_HRSL
	    slop%MC%p%xptar%used = slop_param_p_HRSL
	  else if (hadron_arm.eq.5 .or. hadron_arm.eq.6) then
	    slop%MC%p%delta%used = slop_param_d_SHMS
	    slop%MC%p%yptar%used = slop_param_t_SHMS
	    slop%MC%p%xptar%used = slop_param_p_SHMS
	  endif
	endif

! ... Add slop to SPedges.  Used as input to final edge.*.*.* and to
! ... calculate minimum/maximum theta_phys for extreme_trip_thru_target.

	SPedge%e%delta%min = SPedge%e%delta%min - slop%MC%e%delta%used
	SPedge%e%delta%max = SPedge%e%delta%max + slop%MC%e%delta%used
	SPedge%e%yptar%min = SPedge%e%yptar%min - slop%MC%e%yptar%used
	SPedge%e%yptar%max = SPedge%e%yptar%max + slop%MC%e%yptar%used
	SPedge%e%xptar%min = SPedge%e%xptar%min - slop%MC%e%xptar%used
	SPedge%e%xptar%max = SPedge%e%xptar%max + slop%MC%e%xptar%used
	SPedge%p%delta%min = SPedge%p%delta%min - slop%MC%p%delta%used
	SPedge%p%delta%max = SPedge%p%delta%max + slop%MC%p%delta%used
	SPedge%p%yptar%min = SPedge%p%yptar%min - slop%MC%p%yptar%used
	SPedge%p%yptar%max = SPedge%p%yptar%max + slop%MC%p%yptar%used
	SPedge%p%xptar%min = SPedge%p%xptar%min - slop%MC%p%xptar%used
	SPedge%p%xptar%max = SPedge%p%xptar%max + slop%MC%p%xptar%used

! Compute TRUE edges -- distortions in the target come into play

	edge%e%E%min = (1.+SPedge%e%delta%min/100.)*spec%e%P +
     >		targ%Coulomb%min - dE_edge_test
	edge%e%E%max = (1.+SPedge%e%delta%max/100.)*spec%e%P +
     >		targ%Coulomb%max + dE_edge_test
	pp%min = (1.+SPedge%p%delta%min/100.)*spec%p%P - dE_edge_test
	pp%max = (1.+SPedge%p%delta%max/100.)*spec%p%P + dE_edge_test
	pp%min = max(0.001e0,pp%min)  !avoid p=0 (which can lead to div by zero)(
	edge%p%E%min = sqrt(pp%min**2 + Mh2)
	edge%p%E%max = sqrt(pp%max**2 + Mh2)

! ... extreme theta_e,theta_p,z values for extreme_trip_thru_target.
	the_phys%max = acos( (spec%e%cos_th-spec%e%sin_th*SPedge%e%yptar%max)/
     >		sqrt(1.+SPedge%e%yptar%max**2+SPedge%e%xptar%max**2) )
	the_phys%min = acos( (spec%e%cos_th-spec%e%sin_th*SPedge%e%yptar%min)/
     >		sqrt(1.+SPedge%e%yptar%min**2) )
	thp_phys%max = acos( (spec%p%cos_th-spec%p%sin_th*SPedge%p%yptar%max)/
     >		sqrt(1.+SPedge%p%yptar%max**2+SPedge%p%xptar%max**2) )
	thp_phys%min = acos( (spec%p%cos_th-spec%p%sin_th*SPedge%p%yptar%min)/
     >		sqrt(1.+SPedge%p%yptar%min**2) )
	z%min = -0.5*targ%length
	z%max =  0.5*targ%length
	call extreme_trip_thru_target(Ebeam, the_phys, thp_phys, edge%e%E,
     >                                pp, z, Mh)
	if (.not.using_Eloss) then
	  do i = 1, 3
	    targ%Eloss(i)%min = 0.0
	    targ%Eloss(i)%max = 0.0
	  enddo
	endif

	if (.not.mc_smear) then
	  targ%musc_max(1)=0.
	  targ%musc_max(2)=0.
	  targ%musc_max(3)=0.
	endif

	edge%e%E%min = edge%e%E%min + targ%Eloss(2)%min
	edge%e%E%max = edge%e%E%max + targ%Eloss(2)%max
	edge%p%E%min = edge%p%E%min + targ%Eloss(3)%min
	edge%p%E%max = edge%p%E%max + targ%Eloss(3)%max
	edge%e%yptar%min = SPedge%e%yptar%min - targ%musc_max(2)
	edge%e%yptar%max = SPedge%e%yptar%max + targ%musc_max(2)
	edge%e%xptar%min = SPedge%e%xptar%min - targ%musc_max(2)
	edge%e%xptar%max = SPedge%e%xptar%max + targ%musc_max(2)
	edge%p%yptar%min = SPedge%p%yptar%min - targ%musc_max(3)
	edge%p%yptar%max = SPedge%p%yptar%max + targ%musc_max(3)
	edge%p%xptar%min = SPedge%p%xptar%min - targ%musc_max(3)
	edge%p%xptar%max = SPedge%p%xptar%max + targ%musc_max(3)

! Edges on values of Em and Pm BEFORE reconstruction% Need to apply slop to
! take into account all transformations from ORIGINAL TRUE values to
! RECONSTRUCTED TRUE values --> that includes: (a) reconstruction slop on
! spectrometer values (b) variation in the beam energy (c) difference between
! actual losses and the corrections made for them

! %%. Save the 'measured' beam energy, which will be used in recon.

	Ebeam_vertex_ave = Ebeam + targ%Coulomb%ave - targ%Eloss(1)%ave

! ... Calculate slop in energies and momenta.

	Ebeam_max = Ebeam + dEbeam/2. - targ%Eloss(1)%min + targ%Coulomb%max
	Ebeam_min = Ebeam - dEbeam/2. - targ%Eloss(1)%max + targ%Coulomb%min
	slop_Coulomb = targ%Coulomb%max - targ%Coulomb%ave
	slop_Ebeam = dEbeam/2. + slop_Coulomb
	slop_Ee = slop%MC%e%delta%used/100.*spec%e%P + slop_Coulomb
	r = sqrt(edge%p%E%max**2 - Mh2)
	slop_Ep = sqrt( (r + slop%MC%p%delta%used/100.*spec%p%P)**2 + Mh2 ) -
     >		edge%p%E%max
	slop_Ebeam = slop_Ebeam + (targ%Eloss(1)%max-targ%Eloss(1)%min)
	slop_Ee = slop_Ee + (targ%Eloss(2)%max-targ%Eloss(2)%min)
	slop_Ep = slop_Ep + (targ%Eloss(3)%max-targ%Eloss(3)%min)

	if (doing_heavy) then		! 'reconstructed' Em cuts.
	  slop%total%Em%used = slop_Ebeam + slop_Ee + slop_Ep + dE_edge_test
	  edge%Em%min = cuts%Em%min - slop%total%Em%used
	  edge%Em%max = cuts%Em%max + slop%total%Em%used
	  edge%Em%min = max(0.e0,edge%Em%min)
	endif

! Edges on Em, Pm, etc... VERTEXedge.* values are vertex limits.  edge.* values
! are limits at spectrometer (after eloss and e'/hadron radiation).  Remember
! that min/max values are initialized to wide open values in STRUCTURES.
! Need edge.Em to limit radiated photon energies.
! Need VERTEXedge.Em for calulating limits on radiatied photons.
! Need VERTEXedge.Pm for A(e,e'p) only (for generation limits on Ee', Ep').
! Note that Pm_theory(*).min/max might allow for positive and negative Pm,
! or it could have positive only.  We want *.Pm.min to be zero if zero is
! allowed, so we have to check for Pm_theory.min being negative.

! For all cases:
! 1. Start by giving Em, Pm limits (=0 for hydrogen, based on S.F. limits
!    for all other cases except A(e,e'p).  For this case, Pm limits come
!    from spectral function.  Em limit (max) has to come from energy
!    conservation after everything else is calculated.
!

	if (doing_hyd_elast) then
	  VERTEXedge%Em%min = 0.0
	  VERTEXedge%Em%max = 0.0
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = 0.0
	else if (doing_deuterium) then
	  VERTEXedge%Em%min = Mp + Mn - targ%M		!2.2249 MeV, I hope.
	  VERTEXedge%Em%max = Mp + Mn - targ%M
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = max(abs(Pm_theory(1)%min),abs(Pm_theory(1)%max))
	else if (doing_heavy) then
	  VERTEXedge%Pm%min=0.0
	  VERTEXedge%Pm%max=0.0
	  if(use_benhar_sf) then
	     VERTEXedge%Pm%max=790.0
	  else
	     do i = 1, nrhoPm
		t1=max(abs(Pm_theory(i)%min),abs(Pm_theory(i)%max))
		VERTEXedge%Pm%max = max(VERTEXedge%Pm%max,t1)
	     enddo
	  endif
	  VERTEXedge%Em%min = E_Fermi
	  VERTEXedge%Em%max = 1000.	!Need Egamma_tot_max for good limit.
	else if (doing_hydpi .or. doing_hydkaon .or. doing_hyddelta .or. doing_hydrho) then
	  VERTEXedge%Em%min = 0.0
	  VERTEXedge%Em%max = 0.0
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = 0.0
	else if (doing_deutpi .or. doing_deutkaon .or. doing_deutdelta .or. doing_deutrho) then
	  VERTEXedge%Em%min = Mp + Mn - targ%M		!2.2249 MeV, I hope.
	  VERTEXedge%Em%max = Mp + Mn - targ%M
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = pval(nump)
	else if (doing_hepi .or. doing_hekaon .or. doing_hedelta .or. doing_herho) then
	  VERTEXedge%Em%min = targ%Mtar_struck + targ%Mrec - targ%M
	  VERTEXedge%Em%max = Emval(numEm)
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = pval(nump)
	else if (doing_semi) then
	  VERTEXedge%Em%min = 0.0
	  VERTEXedge%Em%max = 0.0
	  VERTEXedge%Pm%min = 0.0
	  VERTEXedge%Pm%max = 0.0
	   
	endif
	write(6,*) 'E_bind =',VERTEXedge%Em%min,'MeV in limits_init (QF only)'


! Calculate limits for recoiling (A-1) system, if there is one.
! M_{A-1} (Mrec) limits from Em range, since M_{A-1} = M_A - M_N + E_m
! T_{A-1} (Trec) limits from Mrec and Pm limits, since
!    E_rec=sqrt(M_rec**2+P_rec**2), and P_rec = -P_m

	if (doing_hyd_elast .or. doing_hydpi .or. doing_hydkaon .or. 
     >        doing_hyddelta .or. doing_hydrho .or. doing_semi) then
	  VERTEXedge%Mrec%min = 0.0
	  VERTEXedge%Mrec%max = 0.0
	  VERTEXedge%Trec%min = 0.0
	  VERTEXedge%Trec%max = 0.0
	else
	  VERTEXedge%Mrec%min = targ%M - targ%Mtar_struck + VERTEXedge%Em%min
	  VERTEXedge%Mrec%max = targ%M - targ%Mtar_struck + VERTEXedge%Em%max
	  VERTEXedge%Trec%min = sqrt(VERTEXedge%Mrec%max**2+VERTEXedge%Pm%min**2) - 
     >		VERTEXedge%Mrec%max
	  VERTEXedge%Trec%max = sqrt(VERTEXedge%Mrec%min**2+VERTEXedge%Pm%max**2) -
     >		VERTEXedge%Mrec%min
	endif

! For pion(kaon) production, Trec_struck is T of the recoiling nucleon(hyperon).
! Can only get limits from general energy conservation, and can only get
! upper limit, since the lower limit is determined by the allowed radiation,
! which is not calculated yet (and needs Trec to be calculated).

	if (doing_eep .or. doing_semi) then
	  VERTEXedge%Trec_struck%min = 0.
	  VERTEXedge%Trec_struck%max = 0.
	else
	  VERTEXedge%Trec_struck%min = 0.
	  VERTEXedge%Trec_struck%max = Ebeam_max + targ%Mtar_struck -
     >		targ%Mrec_struck - edge%e%E%min - edge%p%E%min -
     >		VERTEXedge%Em%min - VERTEXedge%Trec%min
	endif

! Get radiation limits.  In all cases, total energy conservation gives
! the primary limit. The doing_heavy case has a second condition.  Radiated
! photons change Em from the vertex value to the measured value.  They also
! change Pm, which can modify Trec.  So the max. change in Em is for a
! minimum generated Em (VERTEXedge.Em.min) that ends up as a maximum measured Em
! (edge%Em.max), with slop to take into account the modification to Trec.

	Egamma_tot_max = Ebeam_max + targ%Mtar_struck - targ%Mrec_struck -
     >		edge%e%E%min - edge%p%E%min - VERTEXedge%Em%min -
     >		VERTEXedge%Trec%min - VERTEXedge%Trec_struck%min
	if (doing_heavy) then
	  t2 = (edge%Em%max - VERTEXedge%Em%min) +
     >		(VERTEXedge%Trec%max - VERTEXedge%Trec%min)
	  Egamma_tot_max = min(Egamma_tot_max,t2)
	endif

! ... override calculated limits with hardwired value if desired.
	if (hardwired_rad) Egamma_tot_max = Egamma_gen_max
	if (.not.using_rad) Egamma_tot_max = 0.0
	if (doing_tail(1)) Egamma1_max = Egamma_tot_max
	if (doing_tail(2)) Egamma2_max = Egamma_tot_max
	if (doing_tail(3)) Egamma3_max = Egamma_tot_max

	if (doing_heavy) then	!Needed Egamma_tot_max for final limits.
	  VERTEXedge%Em%min = max(VERTEXedge%Em%min,edge%Em%min-Egamma_tot_max)
	  VERTEXedge%Em%max = min(VERTEXedge%Em%max,edge%Em%max)
	endif

! ... compute edge on summed _generated_ energies based on predicted VERTEX
! ... and TRUE values of Em (Ee'+Ep' for doing_heavy, Ee' for pion/kaon,
! ... Ee' for D(e,e'p), not used for hydrogen elastic.)
! ... Primary limits from energy conservation at vertex, seconday limits from
! ... spectrometer limits modified by radiation.

	if (doing_hyd_elast) then	!NO generated energies.
	  gen%sumEgen%min = 0.0
	  gen%sumEgen%max = 0.0

	else if (doing_heavy) then	!generated TOTAL (e+p) energy limits.
	  gen%sumEgen%max = Ebeam_max + targ%Mtar_struck -
     >		VERTEXedge%Trec%min - VERTEXedge%Em%min
	  gen%sumEgen%min = Ebeam_min + targ%Mtar_struck -
     >		VERTEXedge%Trec%max - VERTEXedge%Em%max - Egamma1_max 
	  gen%sumEgen%max = min(gen%sumEgen%max,edge%e%E%max+edge%p%E%max+Egamma_tot_max)
	  gen%sumEgen%min = max(gen%sumEgen%min,edge%e%E%min+edge%p%E%min)

	else if (doing_semi) then
c	   gen%sumEgen%max = Ebeam_max - VERTEXedge%Trec%min - VERTEXedge%Trec_struck%min
c	   gen%sumEgen%min = Ebeam_min - VERTEXedge%Trec%max - VERTEXedge%Trec_struck%max
	   gen%sumEgen%max = Ebeam_max + targ%Mtar_struck - targ%Mrec_struck
	   gen%sumEgen%min = edge%e%E%min + edge%p%E%min

	else				!generated ELECTRON energy limits.
	  gen%sumEgen%max = Ebeam_max + targ%Mtar_struck - targ%Mrec_struck -
     >		edge%p%E%min - VERTEXedge%Em%min - VERTEXedge%Trec%min -
     >		VERTEXedge%Trec_struck%min
	  gen%sumEgen%min = Ebeam_min + targ%Mtar_struck - targ%Mrec_struck -
     >		edge%p%E%max - VERTEXedge%Em%max - VERTEXedge%Trec%max -
     >		VERTEXedge%Trec_struck%max - Egamma_tot_max
	  gen%sumEgen%max = min(gen%sumEgen%max, edge%e%E%max+Egamma2_max)
	  gen%sumEgen%min = max(gen%sumEgen%min, edge%e%E%min)
	endif

	gen%sumEgen%min = gen%sumEgen%min - dE_edge_test
	gen%sumEgen%max = gen%sumEgen%max + dE_edge_test
	gen%sumEgen%min = max(0.e0,gen%sumEgen%min)

! ... E arm GENERATION limits from sumEgen.
! ... Not used for doing_hyd_elast, but define for the hardwired histograms.

	if (doing_hyd_elast) then
	  gen%e%E%min = edge%e%E%min
	  gen%e%E%max = edge%e%E%max + Egamma2_max
	else if (doing_deuterium .or. doing_pion .or. doing_kaon 
     >      .or. doing_rho .or. doing_delta) then
	  gen%e%E%min = gen%sumEgen%min
	  gen%e%E%max = gen%sumEgen%max
	else if (doing_heavy .or. doing_semi) then
	  gen%e%E%min = gen%sumEgen%min - edge%p%E%max - Egamma3_max
	  gen%e%E%max = gen%sumEgen%max - edge%p%E%min
	endif

! ... Apply limits from direct comparison to acceptance.
	gen%e%E%min = max(gen%e%E%min, edge%e%E%min)
	gen%e%E%max = min(gen%e%E%max, edge%e%E%max + Egamma2_max)

	gen%e%delta%min = (gen%e%E%min/spec%e%P-1.)*100.
	gen%e%delta%max = (gen%e%E%max/spec%e%P-1.)*100.
	gen%e%yptar%min = edge%e%yptar%min
	gen%e%yptar%max = edge%e%yptar%max
	gen%e%xptar%min = edge%e%xptar%min
	gen%e%xptar%max = edge%e%xptar%max

! ... P arm GENERATION limits from sumEgen.  Not used for any case
! ... except doing_heavy, but need to define for code that writes out limits.

	if (doing_hyd_elast.or.doing_deuterium.or.doing_pion.or.doing_kaon .or.
     >    doing_rho .or. doing_delta) then
	  gen%p%E%min = edge%p%E%min
	  gen%p%E%max = edge%p%E%max + Egamma3_max
	else if (doing_heavy .or. doing_semi)then
	  gen%p%E%min = gen%sumEgen%min - edge%e%E%max - Egamma2_max
	  gen%p%E%max = gen%sumEgen%max - edge%e%E%min
	endif
! ... Apply limits from direct comparison to acceptance.
        gen%p%E%min = max(gen%p%E%min, edge%p%E%min)
        gen%p%E%max = min(gen%p%E%max, edge%p%E%max + Egamma3_max)

	gen%p%delta%min = (sqrt(gen%p%E%min**2-Mh2)/spec%p%P-1.)*100.
	gen%p%delta%max = (sqrt(gen%p%E%max**2-Mh2)/spec%p%P-1.)*100.
	gen%p%yptar%min = edge%p%yptar%min
	gen%p%yptar%max = edge%p%yptar%max
	gen%p%xptar%min = edge%p%xptar%min
	gen%p%xptar%max = edge%p%xptar%max

! Axis specs for the diagnostic histograms
	H%gen%e%delta%min =  gen%e%delta%min
	H%gen%e%yptar%min =  gen%e%yptar%min
	H%gen%e%xptar%min = -gen%e%xptar%max
	H%gen%e%delta%bin = (gen%e%delta%max-gen%e%delta%min)/float(nHbins)
	H%gen%e%yptar%bin = (gen%e%yptar%max-gen%e%yptar%min)/float(nHbins)
	H%gen%e%xptar%bin = (gen%e%xptar%max-gen%e%xptar%min)/float(nHbins)
	H%gen%p%delta%min =  gen%p%delta%min
	H%gen%p%yptar%min =  gen%p%yptar%min
	H%gen%p%xptar%min = -gen%p%xptar%max
	H%gen%p%delta%bin = (gen%p%delta%max-gen%p%delta%min)/float(nHbins)
	H%gen%p%yptar%bin = (gen%p%yptar%max-gen%p%yptar%min)/float(nHbins)
	H%gen%p%xptar%bin = (gen%p%xptar%max-gen%p%xptar%min)/float(nHbins)
	H%gen%Em%min = VERTEXedge%Em%min
	H%gen%Em%bin = (max(100.e0,VERTEXedge%Em%max) - VERTEXedge%Em%min)/float(nHbins)
	H%gen%Pm%min = VERTEXedge%Pm%min
	H%gen%Pm%bin = (max(100.e0,VERTEXedge%Pm%max) - VERTEXedge%Pm%min)/float(nHbins)

	H%geni%e%delta%min = H%gen%e%delta%min
	H%geni%e%yptar%min = H%gen%e%yptar%min
	H%geni%e%xptar%min = H%gen%e%xptar%min
	H%geni%e%delta%bin = H%gen%e%delta%bin
	H%geni%e%yptar%bin = H%gen%e%yptar%bin
	H%geni%e%xptar%bin = H%gen%e%xptar%bin
	H%geni%p%delta%min = H%gen%p%delta%min
	H%geni%p%yptar%min = H%gen%p%yptar%min
	H%geni%p%xptar%min = H%gen%p%xptar%min
	H%geni%p%delta%bin = H%gen%p%delta%bin
	H%geni%p%yptar%bin = H%gen%p%yptar%bin
	H%geni%p%xptar%bin = H%gen%p%xptar%bin
	H%geni%Em%min = H%gen%Em%min
	H%geni%Em%bin = H%gen%Em%bin
	H%geni%Pm%min = H%gen%Pm%min
	H%geni%Pm%bin = H%gen%Pm%bin

	H%RECON%e%delta%min = H%gen%e%delta%min
	H%RECON%e%yptar%min = H%gen%e%yptar%min
	H%RECON%e%xptar%min = H%gen%e%xptar%min
	H%RECON%e%delta%bin = H%gen%e%delta%bin
	H%RECON%e%yptar%bin = H%gen%e%yptar%bin
	H%RECON%e%xptar%bin = H%gen%e%xptar%bin
	H%RECON%p%delta%min = H%gen%p%delta%min
	H%RECON%p%yptar%min = H%gen%p%yptar%min
	H%RECON%p%xptar%min = H%gen%p%xptar%min
	H%RECON%p%delta%bin = H%gen%p%delta%bin
	H%RECON%p%yptar%bin = H%gen%p%yptar%bin
	H%RECON%p%xptar%bin = H%gen%p%xptar%bin
	H%RECON%Em%min = H%gen%Em%min
	H%RECON%Em%bin = H%gen%Em%bin
	H%RECON%Pm%min = H%gen%Pm%min
	H%RECON%Pm%bin = H%gen%Pm%bin

	return
	end

!------------------------------------------------------------------------

	subroutine radc_init

	implicit none
	include 'simulate.inc'
	include 'radc.inc'
	include 'brem.inc'

!--------------------------------------------------------------
!
! First, about those (mysterious) 2 main 'radiative option' flags ...
!
! The significance of RAD_FLAG:
!	 RAD_FLAG  = 0	.. use best available formulas, generate in
!			.. (ntail,Egamma) basis
!		   = 1	.. use BASICRAD only, generate in (ntail,Egamma)
!			.. basis
!		   = 2	.. use BASICRAD only, generate in (Egamma(1,2,3))
!			.. basis but prevent overlap of tails (bogus, note)
!		   = 3	.. use BASICRAD only, generate in (Egamma(1,2,3))
!			.. allowing radiation in all 3 directions
!			.. simultaneously
! The (ntail,Egamma) basis can be called the PEAKED basis since it allows
! only 3 photon directions. PEAKED_BASIS_FLAG is set to zero below when
! the peaked basis is being used, in this way we can conveniently tell
! the BASICRAD routine to use the full Egamma formula to generate the gamma
! energy whenever it's called.
!
! (See N. Makins' thesis, section 4.5.6, for more help on this point)
!
! The significance of EXTRAD_FLAG:
!	EXTRAD_FLAG = 1	.. use only BASIC external radiation formulas
!			.. (phi = 1)
!		    = 2	.. use BASIC ext rad formulas x phi
!		    = 3 .. use Friedrich approximation the way we always
!			.. have
!		    = 0 .. use DEFAULTS: 3 for RAD_FLAG = 0, 1 otherwise; note
!			   that the defaults mimic the hardwired 'settings'
!			   in SIMULATE, which doesnt read EXTRAD_FLAG
!			   but determines what to do based on RAD_FLAG
!--------------------------------------------------------------

! Check setting of EXTRAD_FLAG

	if (debug(2)) write(6,*)'radc_init: entering...'
	if (extrad_flag.eq.0) then
	  if (rad_flag.eq.0) then
	    extrad_flag = 3
	  else if (rad_flag.eq.1 .or. rad_flag.eq.2 .or. rad_flag.eq.3) then
	    extrad_flag = 1
	  endif
	else if (extrad_flag.lt.0) then
	  stop 'Imbecile! check your stupid setting of EXTRAD_FLAG'
	endif

! 'etatzai' parameter

	if (debug(4)) write(6,*)'radc_init: at 1'
	etatzai = (12.0+(targ%Z+1.)/(targ%Z*targ%L1+targ%L2))/9.0
	if (debug(4)) write(6,*)'radc_init: etatzai = ',etatzai

! Initialize brem flags (brem doesn't include the normal common blocks)
	produce_output = debug(1)
c	exponentiate = use_expon
	if(use_expon.eq.1) then
	   exponentiate=.true.
	else
	   exponentiate=.false.
	endif

	include_hard = .true.
	calculate_spence = .true.

	if (debug(2)) write(6,*)'radc_init: ending...'
	return
	end

!---------------------------------------------------------------------

	subroutine radc_init_ev (main,vertex)

	implicit none
	include 'structures.inc'
	include 'radc.inc'

	integer		i
	real*8		r, Ecutoff, dsoft, dhard, dsoft_prime
	real*8		lambda_dave, schwinger, brem, bremos
	type(event_main):: main
	type(event):: vertex

! Note that calculate_central calls this with (main,ev) rather than
! (main,vertex).  Since these are just local variables, calling it vertex
! here does not cause any problem, and makes it easier to follow
! modifications to vertex.* variables in later calls.

	real*8 zero
	parameter (zero=0.0e0)	!double precision zero for subroutine calls.

! Compute some quantities that will be needed for rad corr on this event

! ... factor for limiting energy of external radiation along incident electron
!	etta = 1.0 + 2*vertex.ein*sin(vertex.e.theta/2.)**2/(targ%A*amu)
! ... moron move! let's can that etta factor ...

	etta = 1.0

! ... the bt's

	do i=1,2
	  bt(i) = etatzai*main%target%teff(i)
	enddo

! ... the lambda's (effective bt's for internal radiation)

	do i=1,3
	  lambda(i) = lambda_dave(i,1,doing_tail(3),vertex%Ein,vertex%e%E,vertex%p%E,
     >			vertex%p%P,vertex%e%theta)
	enddo
	rad_proton_this_ev = lambda(3).gt.0

! ... get the hard correction factor. don't care about Ecutoff! Just want dhard here

	Ecutoff = 450.
	if (intcor_mode.eq.0) then
	  r = schwinger(Ecutoff,vertex,.true.,dsoft,dhard)
	else
	  if (.not.use_offshell_rad) then
	    r = brem(vertex%Ein,vertex%e%E,Ecutoff,rad_proton_this_ev,dsoft,dhard,
     >		dsoft_prime)
	  else
	    r = bremos(Ecutoff, zero, zero, vertex%Ein, vertex%e%P*vertex%ue%x,
     >		vertex%e%P*vertex%ue%y, vertex%e%P*vertex%ue%z, zero, zero, zero,
     >		vertex%p%P*vertex%up%x, vertex%p%P*vertex%up%y, vertex%p%P*vertex%up%z,
     >		vertex%p%E, rad_proton_this_ev, dsoft, dhard, dsoft_prime)
	  endif
	endif
	hardcorfac = 1./(1.-dhard)
	g(4)=-dsoft_prime*Ecutoff+bt(1)+bt(2)

! ... initialize the parameters needed for our "basic" calculation

	call basicrad_init_ev (vertex%Ein,vertex%e%E,vertex%p%E)

! ... the relative magnitudes of the three tails (we may not need them)

	do i=1,3
	  frac(i) = g(i)/g(0)
	enddo

	return
	end

!-----------------------------------------------------------------------

	subroutine basicrad_init_ev (e1,e2,e3)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	real*8 one
	parameter (one=1.)

	integer i
	real*8 e1,e2,e3,e(3),gamma

	if (debug(2)) write(6,*)'basicrad_init_ev: entering...'
	e(1) = e1
	e(2) = e2
	e(3) = e3

! bt's for internal + external
! ??? One possibility for shutting off tails 1 or
! 2 is to set g(1/2) = 0 here ... note that lambda(3) is set to 0 in
! lambda_dave at the moment, AND proton terms are removed from brem if proton
! radiation off. something analogous and similarly consistent would have to be
! done for the other tails, right now they're just nixed in generate_rad. also
! check ALL ntail.eq.0 checks in kinema constraint lines of generate_rad

	g(1) = lambda(1) + bt(1)
	g(2) = lambda(2) + bt(2)
	g(3) = lambda(3)
	g(0) = g(1)+g(2)+g(3)

! Internal constants

	c_int(1) = lambda(1)/(e(1)*e(2))**(lambda(1)/2.)
	c_int(2) = lambda(2)/(e(1)*e(2))**(lambda(2)/2.)

*	csa NOTE: GN says these masses s/b Mp!

	c_int(3) = lambda(3)/(Mp*e(3))**(lambda(3)/2.)
	do i = 1, 3
	  c_int(i) = c_int(i) * exp(-euler*lambda(i)) / gamma(one+lambda(i))
	enddo
	g_int = lambda(1) + lambda(2) + lambda(3)
	c_int(0) = c_int(1)*c_int(2) * g_int / lambda(1)/lambda(2)

! ... proton radiation could be off

	if (lambda(3).gt.0) c_int(0) = c_int(0) * c_int(3)/lambda(3)
	c_int(0) = c_int(0) * gamma(one+lambda(1)) * gamma(one+lambda(2))
     >		* gamma(one+lambda(3)) / gamma(one+g_int)

! External constants

	do i = 1, 2
	  c_ext(i) = bt(i)/e(i)**bt(i)/gamma(one+bt(i))
	enddo
	c_ext(3) = 0.0
	g_ext = bt(1) + bt(2)
	c_ext(0) = c_ext(1)*c_ext(2) * g_ext / bt(1)/bt(2)
	c_ext(0) = c_ext(0)*gamma(one+bt(1))*gamma(one+bt(2))/gamma(one+g_ext)

! Internal + external constants

	do i = 1, 2
	  c(i) = c_int(i) * c_ext(i) * g(i)/lambda(i)/bt(i)
     >		* gamma(one+lambda(i))*gamma(one+bt(i))/gamma(one+g(i))
	enddo
	c(3) = c_int(3)

! Finally, constant for combined tails

	c(0) = c(1)*c(2) * g(0)/g(1)/g(2)

! ... proton radiation could be off

	if (g(3).gt.0) c(0) = c(0) * c(3)/g(3)
	c(0)=c(0)*gamma(one+g(1))*gamma(one+g(2))*gamma(one+g(3))/gamma(one+g(0))
	c(4)=g(4)/(e1*e2)**g(4)/gamma(one+g(4))
	if(g(3).gt.0) c(4)=c(4)/e3**g(4)
	if (debug(2)) write(6,*)'basicrad_init_ev: ending...'
	return
	end

!----------------------------------------------------------------------
!	theory_init:
!
! Input spectral function(NOT called for H(e,e'p) or pion/kaon production!)
!
! Note: Some flags that existed in old A(e,e'p) code (e.g. DD's version)
! have been removed.  Code has been changed to correspond to the following
! settings for the old flags:
!	norm_Ji_tail = 0
!	doing_Ji_tail = 0
!	Em_start_Ji_tail = 0.
!	fixed_gamma = .true.

	subroutine theory_init(success)
	
	implicit none
	include 'simulate.inc'

	real*8 Pm_values(ntheorymax), absorption
	integer m,n,iok
	logical success

! ... open the file
	if ( nint(targ%A) .eq. 2) then
	  theory_file='h2.theory'
	else if ( nint(targ%A) .eq. 12) then
	  theory_file='c12.theory'
	else if ( nint(targ%A) .eq. 56) then
	  theory_file='fe56.theory'
	else if ( nint(targ%A) .eq. 197) then
	  theory_file='au197.theory'
	else
	  write(6,*) 'No Theory File (spectral function) for A = ',targ%A
	  write(6,*) 'Defaulting to c12.theory'
	  theory_file='c12.theory'
	endif
c	open(unit=1,file=theory_file,status='old',readonly,shared,iostat=iok)
	open(unit=1,file=theory_file,status='old',iostat=iok)

! ... read in the theory file
	read(1,*,err=40) nrhoPm, absorption, E_Fermi
	do m=1, nrhoPm
	  read(1,*,err=40) nprot_theory(m), Em_theory(m), Emsig_theory(m),
     >		bs_norm_theory(m)
	  nprot_theory(m) = nprot_theory(m) * absorption
	enddo
	read(1,*,err=40) Pm_values(1),theory(1,1)
	do m=1, nrhoPm
	  n=2
	  read(1,*,err=40,end=50) Pm_values(2),theory(m,2)
	  do while (Pm_values(n).gt.Pm_values(n-1))
	    n=n+1
	    read(1,*,err=40,end=50) Pm_values(n),theory(m,n)
	  enddo

! ........ figure out details of the Pm axes
50	  Pm_theory(m)%n=n-1
	  Pm_theory(m)%bin=Pm_values(2)-Pm_values(1)
	  Pm_theory(m)%min=Pm_values(1)-Pm_theory(m)%bin/2.
	  Pm_theory(m)%max=Pm_values(Pm_theory(m)%n)+Pm_theory(m)%bin/2.
	  if (abs(Pm_theory(m)%min+Pm_theory(m)%bin*Pm_theory(m)%n -
     >  	Pm_theory(m)%max) .gt. 0.1) then
		write(6,'(1x,''ERROR: theory_init found unequal Pm bins in distribution number '',i2,''!'')') m
		close(1)
		return
	  endif

! ........ prepare for the next loop
	  Pm_values(1) = Pm_values(n)
	  theory(m+1,1) = theory(m,n)
	enddo

! ... are we doing deuterium? (i.e. only using a 1D spectral function)
	doing_deuterium = nrhoPm.eq.1 .and. E_Fermi.lt.1.0

! ... renormalize the momentum distributions
	do m=1, nrhoPm
	  do n=1, Pm_theory(m)%n
	    theory(m,n) = theory(m,n)/bs_norm_theory(m)
	  enddo
	enddo

! ... and calculate the integral of the Em distribution above E_Fermi to
! ... renormalize
	if (doing_heavy) then
	  do m=1, nrhoPm
	    Em_int_theory(m) = 1.
	    Em_int_theory(m) = (pi/2. + atan((Em_theory(m)-E_Fermi)/
     >		(0.5*Emsig_theory(m))))/pi
	  enddo
	endif

! ... we made it
	success=.true.
	close(1)
	return

! ... oops
40	continue
	write(6,*) 'ERROR: theory_init failed to read in the theory file'

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine min_max_init(contrib)

C Gaskell- April 15, 2009
C gfortran does not allow default initalization of a "derived type" that is
C included in a common block.
C For example, all the blah%min and blah%max variables were formerly initialized by 
C default to -1d10 and 1d10 in simulate_init.inc. That is no longer allowed.
C Here I'm attempting to initialize each variable by hand to those values. In principle,
C this should not be necessary as they should be redefined elsewhere later - but there
C are a lot of different cases and something may have slipped through the cracks. That
C may have happened here as well.
C   

      include 'simulate.inc'
      real*8 xmin,xmax
      type(contribtype)::	contrib

      xmin=-1.0d10
      xmax=1.0d10

C Initialize "gen"

      gen%e%delta%min=xmin
      gen%e%delta%max=xmax

      gen%e%yptar%min=xmin
      gen%e%yptar%max=xmax

      gen%e%xptar%min=xmin
      gen%e%xptar%max=xmax

      gen%e%E%min=xmin
      gen%e%E%max=xmax

      gen%p%delta%min=xmin
      gen%p%delta%max=xmax

      gen%p%yptar%min=xmin
      gen%p%yptar%max=xmax

      gen%p%xptar%min=xmin
      gen%p%xptar%max=xmax

      gen%p%E%min=xmin
      gen%p%E%max=xmax

      gen%sumEgen%min=xmin
      gen%sumEgen%max=xmax

      gen%Trec%min=xmin
      gen%Trec%max=xmax

C Initialize "cuts"

      cuts%Em%min=xmin
      cuts%Em%max=xmax

      cuts%Pm%min=xmin
      cuts%Pm%max=xmax

C Initialize "Edge"

      edge%e%E%min=xmin
      edge%e%E%max=xmax

      edge%e%yptar%min=xmin
      edge%e%yptar%max=xmax

      edge%e%xptar%min=xmin
      edge%e%xptar%max=xmax

      edge%p%E%min=xmin
      edge%p%E%max=xmax

      edge%p%yptar%min=xmin
      edge%p%yptar%max=xmax

      edge%p%xptar%min=xmin
      edge%p%xptar%max=xmax

      edge%Em%min=xmin
      edge%Em%max=xmax

      edge%Pm%min=xmin
      edge%Pm%max=xmax

      edge%Mrec%min=xmin
      edge%Mrec%max=xmax

      edge%Trec%min=xmin
      edge%Trec%max=xmax

      edge%Trec_struck%min=xmin
      edge%Trec_struck%max=xmax

C Initialize "VERTEXedge"

      VERTEXedge%e%E%min=xmin
      VERTEXedge%e%E%max=xmax

      VERTEXedge%e%yptar%min=xmin
      VERTEXedge%e%yptar%max=xmax

      VERTEXedge%e%xptar%min=xmin
      VERTEXedge%e%xptar%max=xmax

      VERTEXedge%p%E%min=xmin
      VERTEXedge%p%E%max=xmax

      VERTEXedge%p%yptar%min=xmin
      VERTEXedge%p%yptar%max=xmax

      VERTEXedge%p%xptar%min=xmin
      VERTEXedge%p%xptar%max=xmax

      VERTEXedge%Em%min=xmin
      VERTEXedge%Em%max=xmax

      VERTEXedge%Pm%min=xmin
      VERTEXedge%Pm%max=xmax

      VERTEXedge%Mrec%min=xmin
      VERTEXedge%Mrec%max=xmax

      VERTEXedge%Trec%min=xmin
      VERTEXedge%Trec%max=xmax

      VERTEXedge%Trec_struck%min=xmin
      VERTEXedge%Trec_struck%max=xmax
      
C Initialize "SPedge"

      SPedge%e%delta%min=xmin
      SPedge%e%delta%max=xmax

      SPedge%e%yptar%min=xmin
      SPedge%e%yptar%max=xmax

      SPedge%e%xptar%min=xmin
      SPedge%e%xptar%max=xmax

      SPedge%e%z%min=xmin
      SPedge%e%z%max=xmax

      SPedge%p%delta%min=xmin
      SPedge%p%delta%max=xmax

      SPedge%p%yptar%min=xmin
      SPedge%p%yptar%max=xmax

      SPedge%p%xptar%min=xmin
      SPedge%p%xptar%max=xmax

      SPedge%p%z%min=xmin
      SPedge%p%z%max=xmax


C initialize with "lo" a large number and "hi" a small number

C Initialize "contrib%gen"
      contrib%gen%e%delta%lo=xmax
      contrib%gen%e%delta%hi=xmin

      contrib%gen%e%xptar%lo=xmax
      contrib%gen%e%xptar%hi=xmin

      contrib%gen%e%yptar%lo=xmax
      contrib%gen%e%yptar%hi=xmin

      contrib%gen%e%z%lo=xmax
      contrib%gen%e%z%hi=xmin

      contrib%gen%p%delta%lo=xmax
      contrib%gen%p%delta%hi=xmin

      contrib%gen%p%xptar%lo=xmax
      contrib%gen%p%xptar%hi=xmin

      contrib%gen%p%yptar%lo=xmax
      contrib%gen%p%yptar%hi=xmin

      contrib%gen%p%z%lo=xmax
      contrib%gen%p%z%hi=xmin

      contrib%gen%Trec%lo=xmax
      contrib%gen%Trec%hi=xmin

      contrib%gen%sumEgen%lo=xmax
      contrib%gen%sumEgen%hi=xmin

C Initialize "contrib%tru"
      contrib%tru%e%E%lo=xmax
      contrib%tru%e%E%hi=xmin

      contrib%tru%e%yptar%lo=xmax
      contrib%tru%e%yptar%hi=xmin

      contrib%tru%e%xptar%lo=xmax
      contrib%tru%e%xptar%hi=xmin

      contrib%tru%p%E%lo=xmax
      contrib%tru%p%E%hi=xmin

      contrib%tru%p%yptar%lo=xmax
      contrib%tru%p%yptar%hi=xmin

      contrib%tru%p%xptar%lo=xmax
      contrib%tru%p%xptar%hi=xmin

      contrib%tru%Em%lo=xmax
      contrib%tru%Em%hi=xmin

      contrib%tru%Pm%lo=xmax
      contrib%tru%Pm%hi=xmin

      contrib%tru%Trec%lo=xmax
      contrib%tru%Trec%hi=xmin

C Initialize "contrib%SP"
      contrib%SP%e%delta%lo=xmax
      contrib%SP%e%delta%hi=xmin

      contrib%SP%e%yptar%lo=xmax
      contrib%SP%e%yptar%hi=xmin

      contrib%SP%e%xptar%lo=xmax
      contrib%SP%e%xptar%hi=xmin

      contrib%SP%e%z%lo=xmax
      contrib%SP%e%z%hi=xmin

      contrib%SP%p%delta%lo=xmax
      contrib%SP%p%delta%hi=xmin

      contrib%SP%p%yptar%lo=xmax
      contrib%SP%p%yptar%hi=xmin

      contrib%SP%p%xptar%lo=xmax
      contrib%SP%p%xptar%hi=xmin

      contrib%SP%p%z%lo=xmax
      contrib%SP%p%z%hi=xmin

C Initialize "contrib%vertex"

      contrib%vertex%Trec%lo=xmax
      contrib%vertex%Trec%hi=xmin

      contrib%vertex%Em%lo=xmax
      contrib%vertex%Em%hi=xmin

      contrib%vertex%Pm%lo=xmax
      contrib%vertex%Pm%hi=xmin

C Initialize "contrib%rad"
      contrib%rad%Egamma(1)%lo=xmax
      contrib%rad%Egamma(1)%hi=xmin

      contrib%rad%Egamma(2)%lo=xmax
      contrib%rad%Egamma(2)%hi=xmin

      contrib%rad%Egamma(3)%lo=xmax
      contrib%rad%Egamma(3)%hi=xmin

      contrib%rad%Egamma_total%lo=xmax
      contrib%rad%Egamma_total%hi=xmin


      return 

      end
