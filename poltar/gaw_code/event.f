	subroutine limits_update(main,vertex,orig,recon,doing_deuterium,
     >		doing_pion,doing_kaon,contrib,slop)

	include 'structures.inc'
	include 'radc.inc'
	record /event_main/ main
	record /event/ vertex, orig, recon
	record /contrib/ contrib
	record /slop/ slop
	integer i
	logical	doing_deuterium, doing_pion, doing_kaon

! Update the "contribution limits" records

! ... GENERATION values
	call update_range(vertex.e.delta, contrib.gen.e.delta)
	call update_range(vertex.e.yptar, contrib.gen.e.yptar)
	call update_range(vertex.e.xptar, contrib.gen.e.xptar)
	call update_range(vertex.p.delta, contrib.gen.p.delta)
	call update_range(vertex.p.yptar, contrib.gen.p.yptar)
	call update_range(vertex.p.xptar, contrib.gen.p.xptar)
	call update_range(main.Trec, contrib.gen.Trec)

! ........ another tricky shift
	if (doing_deuterium .or. doing_pion .or. doing_kaon) then
	  call update_range(vertex.e.E-main.Ein_shift,contrib.gen.sumEgen)
	else
	  call update_range(vertex.e.E+vertex.p.E-main.Ein_shift,contrib.gen.sumEgen)
	endif

! ... TRUE values
! ........ tricky shift here, remember this'll get compared with edge.e.E,
! ........ compensate for main.target.Coulomb to
! ........ copy the shift made to the edge in generate
	call update_range(orig.e.E-main.Ee_shift, contrib.tru.e.E)
	call update_range(orig.e.xptar, contrib.tru.e.xptar)
	call update_range(orig.e.yptar, contrib.tru.e.yptar)
	call update_range(orig.e.xptar, contrib.tru.e.xptar)
	call update_range(orig.p.E, contrib.tru.p.E)
	call update_range(orig.p.yptar, contrib.tru.p.yptar)
	call update_range(orig.p.xptar, contrib.tru.p.xptar)

! ........ another tricky shift
	call update_range(orig.Em-main.Ein_shift+main.Ee_shift,contrib.tru.Em)
	call update_range(orig.Pm, contrib.tru.Pm)
	call update_range(orig.Trec, contrib.tru.Trec)

! ... SPECTROMETER values
	call update_range(main.SP.e.delta, contrib.SP.e.delta)
	call update_range(main.SP.e.yptar, contrib.SP.e.yptar)
	call update_range(main.SP.e.xptar, contrib.SP.e.xptar)
	call update_range(main.SP.p.delta, contrib.SP.p.delta)
	call update_range(main.SP.p.yptar, contrib.SP.p.yptar)
	call update_range(main.SP.p.xptar, contrib.SP.p.xptar)

! ... VERTEX values
	call update_range(vertex.Trec, contrib.vertex.Trec)
	call update_range(vertex.Em, contrib.vertex.Em)
	call update_range(vertex.Pm, contrib.vertex.Pm)

! ... RADIATION stuff
! ??? should be looking at Egamma2+3 cause we do use limits on that, indirectly

	do i = 1, 3
	  call update_range(Egamma_used(i), contrib.rad.Egamma(i))
	enddo
	call update_range(Egamma_used(1)+Egamma_used(2)+Egamma_used(3),
     >		contrib.rad.Egamma_total)

! Update the "slop limits" records
! ... MC slops
	call update_range(main.RECON.e.delta-main.SP.e.delta,slop.MC.e.delta)
	call update_range(main.RECON.e.yptar-main.SP.e.yptar,slop.MC.e.yptar)
	call update_range(main.RECON.e.xptar-main.SP.e.xptar,slop.MC.e.xptar)
	call update_range(main.RECON.p.delta-main.SP.p.delta,slop.MC.p.delta)
	call update_range(main.RECON.p.yptar-main.SP.p.yptar,slop.MC.p.yptar)
	call update_range(main.RECON.p.xptar-main.SP.p.xptar,slop.MC.p.xptar)

! ... total slops
! ........ that tricky shift again, slops accounted for by the shift not
! ........ included in slop.total.Em.
	call update_range(recon.Em-(orig.Em-main.Ein_shift+main.Ee_shift),
     >		slop.total.Em)
	call update_range(abs(recon.Pm)-abs(orig.Pm), slop.total.Pm)

	return
	end

!-------------------------------------------------------------------

	subroutine update_range(val,range)

	include 'structures.inc'
	record /range/ range
	real*8	val

	range.lo = min(range.lo, val)
	range.hi = max(range.hi, val)

	return
	end

!----------------------------------------------------------------------
! THE routine to GENERATE the (max of 7) random quantities we need to get
! a full description of our event.

	subroutine generate(main,vertex,orig,success)

	implicit none
	include 'simulate.inc'

	integer i, ii
	real*8 Emin, Emax
	real*8 ranprob,ranth1,ranth,ranph,t3,t4,t5,t6
	real*8 pferlo,pferhi
	real*8 gauss1
	logical success
	real*8 grnd		!random # generator.
	record /event_main/ main
	record /event/ vertex, orig

	real*8 nsig_max
	parameter(nsig_max=3.0d0)      !max #/sigma for gaussian ran #s.


! Randomize the position of the interaction inside the available region.
! gen.xwid and gen.ywid are the intrinsic beam widths (one sigma value).
! Use a gaussian beam distribution with +/- 3.0 sigma (add raster afterwards).
	if (debug(2)) write(6,*)'gen: entering'
	main.target.x = gauss1(nsig_max)*gen.xwid+targ.xoffset
	main.target.y = gauss1(nsig_max)*gen.ywid+targ.yoffset

! fr_pattern=1 - rectangle - targ.fr1/fr2 are the x/y raster half-widths.
! fr_pattern=2 - circle - targ.fr1/fr2 are the inner and outer radii.

	if (targ.fr_pattern .eq. 1) then		!square raster
	  t3=grnd()*pi
	  t4=grnd()*pi
	  t5=cos(t3)*targ.fr1
	  t6=cos(t4)*targ.fr2
	elseif (targ.fr_pattern .eq. 2) then	!circular raster
	  t3=grnd()*2.*pi
	  t4=sqrt(grnd())*(targ.fr2-targ.fr1)+targ.fr1
	  t5=cos(t3)*t4
	  t6=sin(t3)*t4
	else						!no raster
	  t5=0.0
	  t6=0.0
	endif

	main.target.x = main.target.x+t5
	main.target.y = main.target.y+t6
	main.target.z = (0.5-grnd())*targ.length+targ.zoffset
	main.target.rastery = t6	!'raster' contribution to vert. pos.
	main.target.rasterx = t5        ! points left as look downstream. GAW

! Take fluctuations of the beam energy into account, and remember to correct
! for ionization loss in the target and Coulomb acceleration of incoming
! electron.  Remove targ.zoffset from the z position of the scattering
! in order to get the position relative to the center of the target.

	call trip_thru_target (1, main.target.z-targ.zoffset, Ebeam, 0.0d0,
     >		main.target.Eloss(1), main.target.teff(1),Me,1)
	if (.not.using_Eloss) main.target.Eloss(1) = 0.0
	if (using_Coulomb) then
	  main.target.Coulomb=targ.Coulomb_constant*(3.-(grnd())**(2./3.))
	else
	  main.target.Coulomb=0.0
	endif
	vertex.Ein = Ebeam + (grnd()-0.5)*dEbeam +
     >		main.target.Coulomb - main.target.Eloss(1)

! ... use KNOWN variations in Ein from Ebeam_vertex_ave and in Ee due to Coulomb
! ... energy to shift a few limits

	main.Ein_shift = vertex.Ein - Ebeam_vertex_ave
	main.Ee_shift = main.target.Coulomb - targ.Coulomb.ave
	edge.e.E.min = edge.e.E.min + main.Ee_shift
	edge.e.E.max = edge.e.E.max + main.Ee_shift
	gen.sumEgen.max = gen.sumEgen.max + main.Ein_shift
	gen.sumEgen.min = gen.sumEgen.min + main.Ein_shift

! Initialize success code and fractional weight due to radiation and
! generation tricks

5	success = .false.
	main.gen_weight = 1.0

! Generated quantities: (phase_space NOT YET IMPLEMENTED).
!
! phase_space: Generate electron E,yptar,xptar and hadron yptar,xptar??
! doing_hyd_elast: fixed Em, generate electron angles
!*! doing_deuterium: fixed Em, generate electron energy and angles,
!*!	p_fermi, and hadron angles
! doing_deuterium: fixed Em, generate electron fully and proton angles, calc Ep
! doing_eep, A>2: generate electron and hadron energy and angles (calc Em/Pm).
! doing_pion: fixed Em, generate electron energy/angles, p_fermi,
!	hadron angles
! doing_kaon: as doing_pion - not implemented yet.
!
! The above is summarized in the following table:
!
!                    ELECTRON                  HADRON
!               ------------------      ------------------
!               E       yptar   xptar   E       yptar   xptar   p_fermi
!
!H(e,e'p)		X	X
!*!D(e,e'p)		X	X				X
!D(e,e'p)	X	X	X		X	X	
!A(e,e'p)	X	X	X	X	X	X
!----------------------------------------------------------------------
!H(e,e'pi)	X	X	X		X	X
!A(e,e'pi)	X	X	X		X	X	X
!H(e,e'K)	X	X	X		X	X
!A(e,e'K)	X	X	X		X	X	X
!----------------------------------------------------------------------
!phase_space	X	X	X	?	X	X
!
! So our procedure is the following:
! 1) Always generate electron yptar and xptar
! 2) generate hadron yptar and xptar for all cases except H(e,e'p), D(e,e'p)
! 3) generate p_fermi for D(e,e'p), and D/He pion/kaon production
! 4) generate electron E for all but hydrogen elastic and deuterium.
! 5) Generate hadron E for A(e,e'p)
! 6) Set missing energy for cases where it is hardwired
!
! After we generate xptar/yptar/energy, we calculate physics angles (theta/phi),
!  momenta, unit vectors, etc... here and/or in complete_ev.
!
! Note that there are also jacobians associated with some and/or all of
! the above.
! 1: We generate uniformly in xptar/yptar, not theta/phi.  We define the
! phase space volume (genvol contribution) as the product of the xptar/yptar
! range, and have a jacobian for each event taking into account the mapping
! between the solid angle on the unit sphere, and the dxptar/dyptar volume
! (the jacobian is 1/cos**3(dtheta), where dtheta is the angle between the
! event and the central spectrometer vector
! 2: For the D(e,e'p), we take Em as fixed in order to calculate the proton
! energy.  There is a jacobian ( |dEp'/dEm| ).  It comes from integrating
! over the energy conservation delta function: delta(E_D - E_p - E_n - Em).

! Generate Electron Angles (all cases):
	vertex.e.yptar=gen.e.yptar.min+grnd()*(gen.e.yptar.max-gen.e.yptar.min)
	vertex.e.xptar=gen.e.xptar.min+grnd()*(gen.e.xptar.max-gen.e.xptar.min)

! Generate Hadron Angles (all but H(e,e'p)):
	if (.not. doing_hyd_elast) then
	  vertex.p.yptar=gen.p.yptar.min+grnd()*(gen.p.yptar.max-gen.p.yptar.min)
	  vertex.p.xptar=gen.p.xptar.min+grnd()*(gen.p.xptar.max-gen.p.xptar.min)
	endif

! Generate Hadron Momentum (A(e,e'p) only).
	if (doing_heavy) then
	  Emin = max(gen.p.E.min, gen.sumEgen.min - gen.e.E.max)
	  Emax = min(gen.p.E.max, gen.sumEgen.max - gen.e.E.min)
	  if (Emin.gt.Emax) goto 100
	  main.gen_weight=main.gen_weight*(Emax-Emin)/(gen.p.E.max-gen.p.E.min)
	  vertex.p.E = Emin + grnd()*(Emax-Emin)
	  vertex.p.P = sqrt(vertex.p.E**2 - Mh2)
	  vertex.p.delta = 100.*(vertex.p.P-spec.p.P)/spec.p.P
	endif

! Generate Fermi Momentum ( D(e,e'pi), He(e,e'pi), D(e,e'K), He(e,e'K) ).
	pfer=0.0
	pferx=0.0
	pfery=0.0
	pferz=0.0
	if(doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon)then
	  ranprob=grnd()
	  ii=1
	  do while (ranprob.gt.mprob(ii) .and. ii.lt.nump)
	    ii=ii+1
	  enddo
	  if (ii.eq.1) then
	    pferlo=0
	  else
	    pferlo=(pval(ii-1)+pval(ii))/2
	  endif
	  if (ii.eq.nump) then
	    pferhi=pval(nump)
	  else
	    pferhi=(pval(ii)+pval(ii+1))/2
	  endif
	  pfer=pferlo+(pferhi-pferlo)*grnd()
	  ranth1=grnd()*2.-1.0
	  ranth=acos(ranth1)
	  ranph=grnd()*2.*pi

	  pferx=sin(ranth)*cos(ranph)
	  pfery=sin(ranth)*sin(ranph)
	  pferz=cos(ranth)
	endif			! generating p_fermi


! Generate Electron Energy (all but hydrogen elastic)
	if (.not.doing_hyd_elast) then
	  Emin=gen.e.E.min
	  Emax=gen.e.E.max
	  if (doing_deuterium .or. doing_pion .or. doing_kaon) then
	    if (.not.using_rad) Emin = max(Emin,gen.sumEgen.min) !WHY NOT RAD?
	    Emax = min(Emax, gen.sumEgen.max)
	  else if (doing_heavy) then		! A(e,e'p)
	    Emin = max(Emin, gen.sumEgen.min - vertex.p.E)	 !WHY RAD?
	    Emax = min(Emax, gen.sumEgen.max - vertex.p.E)
	  endif

	  if (Emin.gt.Emax) goto 100
	  main.gen_weight=main.gen_weight*(Emax-Emin)/(gen.e.E.max-gen.e.E.min)
	  vertex.e.E = Emin + grnd()*(Emax-Emin)
	  vertex.e.P = vertex.e.E
	  vertex.e.delta = 100.*(vertex.e.P-spec.e.P)/spec.e.P

	endif	!not (doing_hyd_elast)


! Calculate the electron and proton PHYSICS angles.

	call physics_angles(spec.e.theta,spec.e.phi,
     &		vertex.e.xptar,vertex.e.yptar,vertex.e.theta,vertex.e.phi)
	call physics_angles(spec.p.theta,spec.p.phi,
     &		vertex.p.xptar,vertex.p.yptar,vertex.p.theta,vertex.p.phi)

! Compute all non-generated quantities

	if (debug(5)) write(6,*)'gen: calling comp_ev with false, main, vertex'
	if (debug(3)) write(6,*)'gen: calling comp_ev with false, main, vertex'
	if (debug(3)) write(6,*)'gen: Ein, E =',vertex.Ein,vertex.e.E
	call complete_ev(main,vertex,success)

	main.sigcc = 1.0

	if (debug(2)) write(6,*)'gen: initial success =',success
	if (.not.success) then

! ... if the calculation failed go back for another choice, note this
! ... should NEVER happen.  Well, by 'NEVER', we really mean it USUALLY
! ... SHOULDN'T happen, but that it will happen under certain circumstances.
! ... The initial electron and proton (or whatever) are generated, and
! ... then the initial nucleon momentum is calculated.  Sometimes (especially
! ... for low beam energies), the spectrometers cover a large enough range
! ... of initial momentum that the solution for the magnatude of the momentum
! ... is negative.  So this will happen for some events, if the kinematics
! ... allow for it.  In some cases, 'NEVER' will mean 3-4 (or 30-40) times per 
! ... event.  After 5000 events, we'll just disable this error message.
! ... This way, you'll be able to tell if the error is occuring, but not
! ... fill up the disk with errors.

	  if (ntried.le.5000) write(6,'(1x,''COMPLETE_EV failed at event '',i8,'' !'')') nevent
	  goto 100
	endif

! ........ temporary storage of Trec for generated event

	main.Trec = vertex.Trec

! Radiate the event, if requested. If we get an event weight of zero (i.e.
! if we CAN'T radiate our event into the acceptance) then success is
! false.  generate_rad also set orig kinematics, and cuts on Em/Pm.
! If not using_rad, must do these here.

	if (using_rad) then
	  call generate_rad(main,vertex,orig,success)
	  if (debug(2)) write(6,*)'gen: after gen_rad, success =',success
	else
	  success = .true.
	  if (doing_heavy) success = 
     >		(vertex.Em .ge. VERTEXedge.Em.min .and.
     >		 vertex.Em .le. VERTEXedge.Em.max .and.
     >		 vertex.Pm .ge. VERTEXedge.Pm.min .and. 
     >		 vertex.Pm .le. VERTEXedge.Pm.max)
	  if (success) then
	    do i = 1, neventfields
	      orig.all(i) = vertex.all(i)
	    enddo
	  endif
	endif

	if (debug(5)) write(6,*) 'orig.Em',orig.Em
	if (debug(2)) write(6,*)'gen: final success =',success

! Shift limits back ... earlier failure always jumps to here

 100	if (debug(2)) write(6,*)'gen: Shift limits back'
	edge.e.E.min = edge.e.E.min - main.Ee_shift
	edge.e.E.max = edge.e.E.max - main.Ee_shift
	gen.sumEgen.max = gen.sumEgen.max - main.Ein_shift
	gen.sumEgen.min = gen.sumEgen.min - main.Ein_shift

	if (debug(2)) write(6,*)'gen: ending'
	return
	end

!---------------------------------------------------------------------

	subroutine complete_ev(main,vertex,success)

	implicit none
	include 'simulate.inc'

	real*8 a, b, c, r, t, QA, QB, QC, radical
	real*8 p_new_x,p_new_y,new_x_x,new_x_y,new_x_z
	real*8 new_y_x,new_y_y,new_y_z,dummy
	real*8 px,py,pz,qx,qy,qz
	real*8 cos_dtheta,y_spec,y_event
	real*8 oop_x,oop_y

	logical success
	record /event_main/ main
	record /event/	vertex

!-----------------------------------------------------------------------
! Calculate everything left in the /event/ structure, given all necessary
!  GENERATION values (some set of xptar,yptar,delta for both arms and p_fermi,
!  and p_fermi, depending on the scattering process: see table in generate.f 
!
! The SINGLE element of /event/ NOT computed here is sigcc.
!
! Another small anomaly is that main.jacobian IS computed here (if not
! doing_recon) ... this is because all the necessary terms have to be
! computed here anyway to calculate /event/ qties.
!
!-----------------------------------------------------------------------

! Initialize

	success = .false.
	main.jacobian = 1.0

! ... unit vector components of outgoing e,p
! ... z is DOWNSTREAM, x is DOWN and y is LEFT looking downstream.

	if (debug(4)) write(6,*)'comp_ev: at 1'
	vertex.ue.x = sin(vertex.e.theta)*cos(vertex.e.phi)
	vertex.ue.y = sin(vertex.e.theta)*sin(vertex.e.phi)
	vertex.ue.z = cos(vertex.e.theta)
	if (.not.doing_hyd_elast) then
	  vertex.up.x = sin(vertex.p.theta)*cos(vertex.p.phi)
	  vertex.up.y = sin(vertex.p.theta)*sin(vertex.p.phi)
	  vertex.up.z = cos(vertex.p.theta)
	endif
	if (debug(4)) write(6,*)'comp_ev: at 2'

! First finish off the e side
! Calculate scattered electron energy for hydrogen/deuterium (e,e'p)

	if (doing_hyd_elast) then
	  vertex.e.E = vertex.Ein*Mp/(Mp+vertex.Ein*(1.-vertex.ue.z))

! old version.  works for hydrogen OR deuterium, if generate p_fermi and
! treat deuterium as hydrogen elastic from moving target.
!	  vertex.e.E = vertex.Ein*(sqrt(Mp2+pfer**2)-pfer*pferz)/
!     >		(sqrt(Mp2+pfer**2)+vertex.Ein*(1.-vertex.ue.z)-
!     >		 pfer*(pferx*vertex.ue.x+pfery*vertex.ue.y+pferz*vertex.ue.z))

	  if (vertex.e.E.gt.vertex.Ein) return
	  vertex.e.P = vertex.e.E
	  vertex.e.delta = (vertex.e.P - spec.e.P)*100./spec.e.P
	  if (debug(4)) write(6,*)'comp_ev: at 3'
	endif

! The q vector

	if (debug(5)) write(6,*)'comp_ev: Ein,E,uez=',vertex.Ein,vertex.e.E,vertex.ue.z
	vertex.q = sqrt(vertex.Ein**2+vertex.e.E**2-2.*vertex.Ein*vertex.e.P*vertex.ue.z)
	vertex.omega = vertex.Ein - vertex.e.E
	vertex.Q2 = vertex.q**2 - vertex.omega**2
	vertex.uq.x = - vertex.e.P*vertex.ue.x / vertex.q
	vertex.uq.y = - vertex.e.P*vertex.ue.y / vertex.q
	vertex.uq.z =(vertex.Ein - vertex.e.P*vertex.ue.z)/ vertex.q
	if (debug(4)) write(6,*)'comp_ev: at 5'

! Now complete the p side.

	if (doing_hyd_elast) then	!p = q
	  vertex.up.x = vertex.uq.x
	  vertex.up.y = vertex.uq.y
	  vertex.up.z = vertex.uq.z
	  vertex.p.theta = acos(vertex.up.z)
	  vertex.p.phi = acos(vertex.up.x/sin(vertex.p.theta))

!...Get the cosine of the angle between the central spectrometer setting
!...and the real event (dot product of p and p0 unit vectors). For in-plane
!...spectrometers, xptar (=dx/dz) is just the x component of p (=dx) 
!...over the component of p along the spectrometer direction (=dz=cos(dtheta)).
!...Knowing xptar, get yptar from cos(dtheta)=1/sqrt(1+xptar**2+yptar**2)
!...Choose sign of yptar by comparing y for the spectrometer to y of the
!...event (projected to the xptar-yptar plane).

	  cos_dtheta = vertex.up.z * cos(spec.p.theta)
     >		+ vertex.up.y * sin(spec.p.theta)*sin(spec.p.phi)
     >		+ vertex.up.x * sin(spec.p.theta)*cos(spec.p.phi)
	  vertex.p.xptar = vertex.up.x / cos_dtheta
	  vertex.p.yptar = sqrt( 1/cos_dtheta**2 - (1 + vertex.p.xptar**2) )
	  y_spec = sin(spec.p.theta)*sin(spec.p.phi)
	  y_event = vertex.up.y/cos_dtheta	!projected to plane perp. to spectrometer.
	  if (y_event .lt. y_spec) vertex.p.yptar = -vertex.p.yptar

	  vertex.p.P = vertex.q
	  vertex.p.E = sqrt(vertex.p.P**2+Mh2)
	  vertex.p.delta = (vertex.p.P - spec.p.P)*100./spec.p.P
	  if (debug(4)) write(6,*)'comp_ev: at 6'

	elseif (doing_deuterium) then	!need Ep, and a jacobian.

	  vertex.Em = Em_theory(1)
	  a = vertex.e.E*(vertex.ue.x*vertex.up.x+vertex.ue.y*vertex.up.y+vertex.ue.z*vertex.up.z) -
     >		vertex.Ein*vertex.up.z
	  b = vertex.Ein**2 + vertex.e.E**2 - 2.*vertex.Ein*vertex.e.E*vertex.ue.z
	  c = vertex.Ein - vertex.Em - vertex.e.E + targ.Mrec + Mp
	  QA = 4.*(a**2 - c**2)
	  QB = 4.*c*(c**2 - b + Mp2 - targ.Mrec**2)
	  QC = -4.*a**2*Mp2 - (c**2 - b + Mp2 - targ.Mrec**2)**2
	  radical = QB**2 - 4.*QA*QC
	  if (radical.lt.0) return
	  vertex.p.E = (-QB - sqrt(radical))/2./QA
	  if (vertex.p.E.le.Mp) return
	  vertex.p.P = sqrt(vertex.p.E**2 - Mp2)
	  vertex.p.delta = (vertex.p.P - spec.p.P)*100./spec.p.P

! ........ the Jacobian here is |dEp'/dEm|
	  t = c**2 - b + Mp2 - targ.Mrec**2
	  main.jacobian = (t*(c-vertex.p.E) + 2*c*vertex.p.E*(vertex.p.E-c)) /
	1		(2*(a**2-c**2)*vertex.p.E + c*t)
	  main.jacobian = abs(main.jacobian)

	elseif (doing_pion) then

	  a = -1.*vertex.q*(vertex.uq.x*vertex.up.x+vertex.uq.y*vertex.up.y+vertex.uq.z*vertex.up.z)
	  b = vertex.q**2
	  c = vertex.omega+targ.M
	  if (doing_deutpi.or.doing_hepi) then	!modify for p_fermi.
	    a = a - abs(pfer)*(pferx*vertex.up.x+pfery*vertex.up.y+pferz*vertex.up.z)
	    b = b+2*vertex.q*abs(pfer)*(pferx*vertex.uq.x+pfery*vertex.uq.y
     >			+pferz*vertex.uq.z)+pfer**2
	    c = c-sqrt(targ.Mrec_pion**2+pfer**2)

!	    IF (DOING_HEPI) THEN
!	      C = C + SQRT(TARG.mREC_PION**2+PFER**2)
!     >		    - SQRT((2*TARG.mREC_PION)**2+PFER**2)
!	    ENDIF

!	    IF (DOING_HEPI) THEN
!	      C = C + SQRT(TARG.mREC_PION**2+PFER**2)
!     >		    - SQRT(TARG.mREC_PION**2+(PFER-50.)**2)
!     >		    - SQRT(TARG.mREC**2+(50.)**2)
!	    ENDIF

	    if (doing_hepi) c = c - targ.Mrec

	  endif
	  QA = 4.*(a**2 - c**2)
	  QB = 4.*c*(c**2 - b + Mh2 - targ.Mrec_pion**2)
	  QC = -4.*a**2*Mh2 - (c**2 - b + Mh2 - targ.Mrec_pion**2)**2

!	  IF (DOING_HEPI) THEN
!	    QB = 4.*c*(c**2 - b + Mh2 - 4.*targ.Mrec_pion**2)
!	    QC = -4.*a**2*Mh2 - (c**2 - b + Mh2 - 4.*targ.Mrec_pion**2)**2
!	  ENDIF

	  radical = QB**2 - 4.*QA*QC
	  if (radical.lt.0) return
	  vertex.p.E = (-QB - sqrt(radical))/2./QA
	  if (vertex.p.E.le.Mh) return
	  vertex.p.P = sqrt(vertex.p.E**2 - Mh2)
	  vertex.p.delta = (vertex.p.P - spec.p.P)*100./spec.p.P
	elseif (doing_kaon) then
	  a = -1.*vertex.q*(vertex.uq.x*vertex.up.x+vertex.uq.y*vertex.up.y+vertex.uq.z*vertex.up.z)
	  b = vertex.q**2
	  c = vertex.omega+targ.M

! For hypernucleus states, treat as doing_hydkaon, but with
! targ.Mtar_kaon=targ.M, targ.Mtar_rec=hypernucleus mass.
* For a test, trying out a modified final state.  Instead of one nucleon
* balancing the momentum of the target nucleon, it is shared between the
* nucleons.  For now, just assuming that the struck nucleon carries all
* but 50 MeV of the momentum.

	  if (doing_deutkaon.or.doing_hekaon) then	!modify for p_fermi.
	    if (which_kaon.lt.10) then
	      a = a - abs(pfer)*(pferx*vertex.up.x+pfery*vertex.up.y+pferz*vertex.up.z)
	      b = b+2*vertex.q*abs(pfer)*(pferx*vertex.uq.x+pfery*vertex.uq.y
     >			+pferz*vertex.uq.z)+pfer**2
*test	      c = c-sqrt(Mn2+pfer**2)
	      c = c-sqrt(Mn2+(pfer-min(pfer,50.))**2)
	    else if (which_kaon.ge.20) then	!Hyperon-Deuteron final state.
	      a = a - abs(pfer)*(pferx*vertex.up.x+pfery*vertex.up.y+pferz*vertex.up.z)
	      b = b+2*vertex.q*abs(pfer)*(pferx*vertex.uq.x+pfery*vertex.uq.y
     >			+pferz*vertex.uq.z)+pfer**2
	      c = c-sqrt(Md2+pfer**2)
	    else
	      write(6,*) 'Bad Mojo!!!!!!!'
	      write(6,*) 'which_kaon = ',which_kaon
	      write(6,*) 'doing_deutkaon = ',doing_deutkaon
	      write(6,*) 'doing_hekaon = ',doing_hekaon
	      stop 'I dont understand this combination'
	    endif
*test	    if (doing_hekaon) c=c-targ.Mrec
	    if (doing_hekaon) c=c-sqrt(targ.Mrec**2+min(pfer,50.)**2)
	  endif
!	write(6,*) 'c-omega-Mp=',c-vertex.omega-Mp
	  QA = 4.*(a**2 - c**2)
	  QB = 4.*c*(c**2 - b + Mh2 - targ.Mrec_kaon**2)
	  QC = -4.*a**2*Mh2 - (c**2 - b + Mh2 - targ.Mrec_kaon**2)**2
	  radical = QB**2 - 4.*QA*QC
	  if (radical.lt.0) return
	  vertex.p.E = (-QB - sqrt(radical))/2./QA
	  if (vertex.p.E.le.Mh) return
	  vertex.p.P = sqrt(vertex.p.E**2 - Mh2)
	  vertex.p.delta = (vertex.p.P - spec.p.P)*100./spec.p.P
	endif

	if(doing_phsp)then
	  vertex.p.P=spec.p.P		!????? single arm phsp??
	  vertex.p.E=sqrt(Mh2+vertex.p.P**2)
	  if (debug(4)) write(6,*)'comp_ev: at 7.5',Mh2,vertex.p.E
	endif
	if (debug(4)) write(6,*)'comp_ev: at 7'

! Compute some pion stuff

	if (doing_pion .or. doing_kaon) then
	  main.epsilon=1./(1. + 2.*(1+vertex.omega**2/vertex.Q2)*tan(vertex.e.theta/2.)**2)
	  main.theta_pq=acos(vertex.up.x*vertex.uq.x+vertex.up.y*vertex.uq.y+vertex.up.z*vertex.uq.z)

! CALCULATE ANGLE PHI BETWEEN SCATTERING PLANE AND REACTION PLANE.
! Therefore, define a new system with the z axis parallel to q, and
! the x axis inside the q-z-plane: z' = q, y' = q X z, x' = y' X q
! this gives phi the way it is usually defined, i.e. phi=0 for in-plane
! particles closer to the downstream beamline than q.
! phi=90 is above the horizontal plane when q points to the right, and
! below the horizontal plane when q points to the left.
! Also take into account the different definitions of x, y and z in
! replay and SIMC:
! As seen looking downstream:  		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	  qx = -vertex.uq.y		!convert to 'replay' coord. system
	  qy =  vertex.uq.x
	  qz =  vertex.uq.z
	  px = -vertex.up.y
	  py =  vertex.up.x
	  pz =  vertex.up.z

	  dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	  new_x_x = -qx*qz/dummy
	  new_x_y = -qy*qz/dummy
	  new_x_z = (qx**2 + qy**2)/dummy

	  dummy   = sqrt(qx**2 + qy**2)
	  new_y_x =  qy/dummy
	  new_y_y = -qx/dummy
	  new_y_z =  0.0

	  p_new_x = px*new_x_x + py*new_x_y + pz*new_x_z
	  p_new_y = px*new_y_x + py*new_y_y + pz*new_y_z

	  if ((p_new_x**2+p_new_y**2).eq.0.) then
	    main.phi_pq = 0.0
	  else
	    main.phi_pq = acos(p_new_x/sqrt(p_new_x**2+p_new_y**2))
	  endif
	  if (p_new_y.lt.0.) then
	    main.phi_pq = 2*3.1415926536 - main.phi_pq
	  endif

	  if (debug(2)) then
	    write(6,*)'comp_ev: omega =',vertex.omega
	    write(6,*)'comp_ev: Q2 =',vertex.Q2
	    write(6,*)'comp_ev: theta_e =',vertex.e.theta
	    write(6,*)'comp_ev: epsilon =',main.epsilon
	    write(6,*)'comp_ev: theta_pq =',main.theta_pq
	    write(6,*)'comp_ev: phi_pq =',main.phi_pq
	  endif
	endif
	if (debug(4)) write(6,*)'comp_ev: at 8'

! Compute the Pm vector in in SIMC LAB system, with x down, and y to the left.
! Computer Parallel, Perpendicular, and Out of Plane componenants.
! Parallel is component along q_hat.  Out/plane is component along
! (e_hat) x (q_hat)  (oop_x,oop_y are components of out/plane unit vector)
! Perp. component is what's left: along (q_hat) x (oop_hat).
! So looking along q, out of plane is down, perp. is left.

	vertex.Pmx = vertex.p.P*vertex.up.x - vertex.q*vertex.uq.x
	vertex.Pmy = vertex.p.P*vertex.up.y - vertex.q*vertex.uq.y
	vertex.Pmz = vertex.p.P*vertex.up.z - vertex.q*vertex.uq.z
	vertex.Pm = sqrt(vertex.Pmx**2+vertex.Pmy**2+vertex.Pmz**2)

!STILL NEED SIGN FOR PmPer!!!!!!

	oop_x = -vertex.uq.y	! oop_x = q_y *(z_hat x y_hat) = q_y *-(x_hat)
	oop_y =  vertex.uq.x	! oop_y = q_x *(z_hat x x_hat) = q_x * (y_hat)
	vertex.PmPar = (vertex.Pmx*vertex.uq.x + vertex.Pmy*vertex.uq.y + vertex.Pmz*vertex.uq.z)
	vertex.PmOop = (vertex.Pmx*oop_x + vertex.Pmy*oop_y) / sqrt(oop_x**2+oop_y**2)
	vertex.PmPer = sqrt( max(0., vertex.Pm**2 - vertex.PmPar**2 - vertex.PmOop**2 ) )

	if (debug(4)) write(6,*)'comp_ev: at 9',vertex.Pmx,vertex.Pmy,vertex.Pmz

! Calculate Trec, Em. Trec for (A-1) system (eep), or for struck nucleon (pi/K)

	if (doing_hyd_elast) then
	  vertex.Em = 0.0
	  vertex.Trec = 0.0
	else if (doing_deuterium) then
	  vertex.Em = Mp + Mn - targ.M
	  vertex.Trec = sqrt(targ.Mrec**2 + pfer**2) - targ.Mrec
	else if (doing_heavy) then
	  vertex.Trec = sqrt(targ.Mrec**2 + vertex.Pm**2)-targ.Mrec
	  vertex.Em = vertex.Ein + Mp - vertex.e.E - vertex.p.E - vertex.Trec
	else if (doing_hydpi) then
	  vertex.Em = 0.0
	  vertex.Trec = vertex.Ein + targ.Mtar_pion - targ.Mrec_pion - vertex.e.E - vertex.p.E
	else if (doing_deutpi) then
	  vertex.Em = Mp + Mn - targ.M
	  vertex.Trec = vertex.Ein + targ.M - targ.Mrec - targ.Mrec_pion -
     >		sqrt(targ.Mrec_pion**2+pfer**2) - vertex.e.E - vertex.p.E
	else if (doing_hepi) then
	  vertex.Em = targ.Z*Mp + targ.N*Mn - targ.M !3-body breakup.
	  vertex.Trec = vertex.Ein + targ.M - targ.Mrec - targ.Mrec_pion -
     >		sqrt(targ.Mrec_pion**2+pfer**2) - vertex.e.E - vertex.p.E
	else if (doing_hydkaon) then
	  vertex.Em = 0.0
	  vertex.Trec = vertex.Ein + targ.Mtar_kaon - targ.Mrec_kaon - vertex.e.E - vertex.p.E
	else if (doing_deutkaon) then
	  vertex.Em = Mp + Mn - targ.M
	  vertex.Trec = vertex.Ein + targ.M - targ.Mrec - targ.Mrec_kaon -
     >		sqrt(Mn2+pfer**2) - vertex.e.E - vertex.p.E
	else if (doing_hekaon) then
	  vertex.Em = targ.Z*Mp + targ.N*Mn - targ.M !3-body breakup.
	  vertex.Trec = vertex.Ein + targ.M - targ.Mrec - targ.Mrec_kaon -
     >		sqrt(Mn2+pfer**2) - vertex.e.E - vertex.p.E
	endif
	if (debug(5)) write(6,*) 'vertex.Pm,vertex.Trec,vertex.Em',vertex.Pm,vertex.Trec,vertex.Em
	if (debug(4)) write(6,*)'comp_ev: at 10'

! Determine PHYSICS scattering angles theta/phi for the two spectrometer
! vectors, and the Jacobian which we must include in our xsec computation
! to take into account the fact that we're not generating in the physics
! solid angle d(omega) but in 'spectrometer angles' th, ph.

	call physics_angles(spec.e.theta,spec.e.phi,
     &		vertex.e.xptar,vertex.e.yptar,vertex.e.theta,vertex.e.phi)
	call physics_angles(spec.p.theta,spec.p.phi,
     &		vertex.p.xptar,vertex.p.yptar,vertex.p.theta,vertex.p.phi)

	r = sqrt(1.+vertex.e.yptar**2+vertex.e.xptar**2) ! 1/cos**3(theta-theta0)
	main.jacobian = main.jacobian / r**3

	if (.not.doing_hyd_elast.and..not.doing_deuterium) then
	  r = sqrt(1.+vertex.p.yptar**2+vertex.p.xptar**2)
	  main.jacobian = main.jacobian / r**3	 ! 1/cos**3(theta-theta0)
	endif

	if (debug(4)) write(6,*)'comp_ev: at 11'

! The effective target thickness that the scattered particles see, and the
! resulting energy losses

	call trip_thru_target (2, main.target.z-targ.zoffset, vertex.e.E,
     >		vertex.e.theta, main.target.Eloss(2), main.target.teff(2),Me,1)
	call trip_thru_target (3, main.target.z-targ.zoffset, vertex.p.E,
     >		vertex.p.theta, main.target.Eloss(3), main.target.teff(3),Mh,1)
	if (.not.using_Eloss) then
	  main.target.Eloss(2) = 0.0
	  main.target.Eloss(3) = 0.0
	endif
	if (debug(4)) write(6,*)'comp_ev: at 12'

! Initialize parameters necessary for radiative corrections

	call radc_init_ev(main,vertex)

! Success code

	success = .true.
	if(debug(2)) write(6,*)'comp_ev: at end...'
	return
	end

!---------------------------------------------------------------------

	subroutine complete_recon_ev(recon,success)

	implicit none
	include 'simulate.inc'

	real*8 p_new_x,p_new_y,new_x_x,new_x_y,new_x_z
	real*8 new_y_x,new_y_y,new_y_z,dummy
	real*8 px,py,pz,qx,qy,qz
	real*8 W2
	real*8 oop_x,oop_y

	logical success
	record /event/	recon

!-----------------------------------------------------------------------
! Calculate everything left in the /event/ structure, given
!	recon.Ein, and the following for both recon.e AND recon.p:
! delta,xptar,yptar, z,P,E,theta,phi (with E,P corrected for eloss, if desired).
!
! The SINGLE element of /event/ NOT computed here is sigcc (for (e,e'p)).
!
!-----------------------------------------------------------------------

! Initialize

	success = .false.


! ... unit vector components of outgoing e,p
! ... z is DOWNSTREAM, x is DOWN and y is LEFT looking downstream.

	if (debug(4)) write(6,*)'comp_rec_ev: at 1'
	recon.ue.x = sin(recon.e.theta)*cos(recon.e.phi)
	recon.ue.y = sin(recon.e.theta)*sin(recon.e.phi)
	recon.ue.z = cos(recon.e.theta)
	recon.up.x = sin(recon.p.theta)*cos(recon.p.phi)
	recon.up.y = sin(recon.p.theta)*sin(recon.p.phi)
	recon.up.z = cos(recon.p.theta)
	if (debug(4)) write(6,*)'comp_rec_ev: at 2'

! The q vector

	if (debug(5)) write(6,*)'comp_rec_ev: Ein,E,uez=',recon.Ein,recon.e.E,recon.ue.z
	recon.q = sqrt(recon.Ein**2+recon.e.E**2-2.*recon.Ein*recon.e.P*recon.ue.z)
	recon.omega = recon.Ein - recon.e.E
	recon.Q2 = recon.q**2 - recon.omega**2
	recon.uq.x = - recon.e.P*recon.ue.x / recon.q
	recon.uq.y = - recon.e.P*recon.ue.y / recon.q
	recon.uq.z =(recon.Ein - recon.e.P*recon.ue.z)/ recon.q
	W2 = targ.M**2 + 2.*targ.M*recon.omega - recon.Q2
	recon.W = sqrt(abs(W2)) * W2/abs(W2) 
	if (debug(4)) write(6,*)'comp_rec_ev: at 5'

! Now complete the p side

	if(doing_phsp)then
	  recon.p.P=spec.p.P
	  recon.p.E=sqrt(Mh2+recon.p.P**2)
	  if (debug(4)) write(6,*)'comp_rec_ev: at 7.5',Mh2,recon.p.E
	endif

! Compute some pion/kaon stuff.

	recon.epsilon=1./(1. + 2.*(1+recon.omega**2/recon.Q2)*tan(recon.e.theta/2.)**2)
	recon.theta_pq=acos(recon.up.x*recon.uq.x+recon.up.y*recon.uq.y+recon.up.z*recon.uq.z)

! CALCULATE ANGLE PHI BETWEEN SCATTERING PLANE AND REACTION PLANE.
! Therefore, define a new system with the z axis parallel to q, and
! the x axis inside the q-z-plane: z' = q, y' = q X z, x' = y' X q
! this gives phi the way it is usually defined, i.e. phi=0 for in-plane
! particles closer to the downstream beamline than q.
! phi=90 is above the horizontal plane when q points to the right, and
! below the horizontal plane when q points to the left.
! Also take into account the different definitions of x, y and z in
! replay and SIMC:
! As seen looking downstream:		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	qx = -recon.uq.y		!convert to 'replay' coord. system
	qy =  recon.uq.x
	qz =  recon.uq.z
	px = -recon.up.y
	py =  recon.up.x
	pz =  recon.up.z

	dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	new_x_x = -qx*qz/dummy
	new_x_y = -qy*qz/dummy
	new_x_z = (qx**2 + qy**2)/dummy

	dummy   = sqrt(qx**2 + qy**2)
	new_y_x =  qy/dummy
	new_y_y = -qx/dummy
	new_y_z =  0.0

	p_new_x = px*new_x_x + py*new_x_y + pz*new_x_z
	p_new_y = px*new_y_x + py*new_y_y + pz*new_y_z

	if ((p_new_x**2+p_new_y**2).eq.0.) then
	  recon.phi_pq = 0.0
	else
	  recon.phi_pq = acos(p_new_x/sqrt(p_new_x**2+p_new_y**2))
	endif
	if (p_new_y.lt.0.) then
	  recon.phi_pq = 2*3.1415926536 - recon.phi_pq
	endif

	if (debug(2)) then
	  write(6,*)'comp_rec_ev: omega =',recon.omega
	  write(6,*)'comp_rec_ev: Q2 =',recon.Q2
	  write(6,*)'comp_rec_ev: theta_e =',recon.e.theta
	  write(6,*)'comp_rec_ev: epsilon =',recon.epsilon
	  write(6,*)'comp_rec_ev: theta_pq =',recon.theta_pq
	  write(6,*)'comp_rec_ev: phi_pq =',recon.phi_pq
	endif

	if (debug(4)) write(6,*)'comp_rec_ev: at 8'

! Compute the Pm vector in in SIMC LAB system, with x down, and y to the left.
! Computer Parallel, Perpendicular, and Out of Plane componenants.
! Parallel is component along q_hat.  Out/plane is component along
! (e_hat) x (q_hat)  (oop_x,oop_y are components of out/plane unit vector)
! Perp. component is what's left: along (q_hat) x (oop_hat).
! So looking along q, out of plane is down, perp. is left.

	recon.Pmx = recon.p.P*recon.up.x - recon.q*recon.uq.x
	recon.Pmy = recon.p.P*recon.up.y - recon.q*recon.uq.y
	recon.Pmz = recon.p.P*recon.up.z - recon.q*recon.uq.z
	recon.Pm = sqrt(recon.Pmx**2+recon.Pmy**2+recon.Pmz**2)

!STILL NEED SIGN FOR PmPer!!!!!!

	oop_x = -recon.uq.y	! oop_x = q_y *(z_hat x y_hat) = q_y *-(x_hat)
	oop_y =  recon.uq.x	! oop_y = q_x *(z_hat x x_hat) = q_x * (y_hat)
	recon.PmPar = (recon.Pmx*recon.uq.x + recon.Pmy*recon.uq.y + recon.Pmz*recon.uq.z)
	recon.PmOop = (recon.Pmx*oop_x + recon.Pmy*oop_y) / sqrt(oop_x**2+oop_y**2)
	recon.PmPer = sqrt( max(0., recon.Pm**2 - recon.PmPar**2 - recon.PmOop**2 ) )

! Calculate Trec, Em. Trec for (A-1) system (eep), or for struck nucleon (pi/K)
! Note that there are other ways to calculate 'Em' for the pion/kaon case, but
! that the missingmass calculation in results_write.f requires this
! definition.  If you change it here, you must fix up results_write.

	if (doing_hyd_elast) then
	  recon.Trec = 0
	  recon.Em = recon.Ein + Mp - recon.e.E - recon.p.E - recon.Trec
	else if (doing_deuterium .or. doing_heavy) then
	  recon.Trec = sqrt(recon.Pm**2+targ.Mrec**2) - targ.Mrec
	  recon.Em = recon.Ein + Mp - recon.e.E - recon.p.E - recon.Trec
	else if (doing_hydpi) then
	  recon.Em = recon.omega + targ.Mtar_pion - recon.p.E
	else if (doing_deutpi .or. doing_hepi) then
	  recon.Em = recon.omega + targ.Mtar_pion - recon.p.E
	else if (doing_hydkaon) then
	  recon.Em = recon.omega + targ.Mtar_kaon - recon.p.E
	else if (doing_deutkaon .or. doing_hekaon) then
	  recon.Em = recon.omega + targ.Mtar_kaon - recon.p.E
	endif

	if (debug(5)) write(6,*) 'recon.Pm,recon.Trec,recon.Em',recon.Pm,recon.Trec,recon.Em
	if (debug(4)) write(6,*)'comp_rec_ev: at 10'

	success=.true.
	return
	end

!------------------------------------------------------------------------

	subroutine complete_main(force_sigcc,main,vertex,recon,success)

	implicit none
	include 'simulate.inc'

	integer		i, iPm1
	real*8		a, b, r, frac, peepi, peeK
	real*8		weight, width, sigep, deForest
	logical		force_sigcc, success
	record /event_main/ main
	record /event/	vertex, recon

!-----------------------------------------------------------------------
! Calculate everything left in the /main/ structure that hasn't been
! required up til now. This routine is called ONLY after we
! know an event IS going to contribute, don't want to be doing these
! calculations if they're not needed!
!
! A small anomaly is that sigcc from the /event/ structure is computed
! here since it's not needed til now, and was not calculated by
! COMPLETE_EV.
!-----------------------------------------------------------------------


! Initialize success code

	if (debug(2)) write(6,*)'comp_main:  entering...'
	success = .false.

! The spectral function weighting

	if (doing_hyd_elast.or.doing_pion.or.doing_kaon.or.doing_phsp) then !no SF.
	  main.SF_weight=1.0
	else
	  main.SF_weight = 0.0
	  do i=1,nrhoPm
	    weight = 0.0

! ... use linear interpolation to determine rho(pm)

	    r = (vertex.Pm-Pm_theory(i).min)/Pm_theory(i).bin
	    if (r.ge.0 .and. r.le.Pm_theory(i).n) then
	      iPm1 = nint(r)
	      if (iPm1.eq.0) iPm1 = 1
	      if (iPm1.eq.Pm_theory(i).n) iPm1 = Pm_theory(i).n-1
	      frac = r+0.5-float(iPm1)
	      b = theory(i,iPm1)
	      a = theory(i,iPm1+1)-b
	      weight = a*frac + b
	    endif

	    if (.not.doing_deuterium) then

	      width=Emsig_theory(i)/2.0
	      if (vertex.Em.lt.E_Fermi) weight = 0.0
	      weight=weight/pi/Em_int_theory(i) * width/
     >				((vertex.Em-Em_theory(i))**2+width**2)
	    endif
	    main.SF_weight = main.SF_weight+weight*nprot_theory(i)
	  enddo ! <i=1,nrhoPm>
	endif

! ... if we have come up with 0 weight, might as well quit now ... unless
! ... FORCE_SIGCC is on.

	if (debug(5))write(6,*)'comp_main: calculating cross section'
	if (main.SF_weight.le.0 .and. .not.force_sigcc) return

! ... Cross section, for BOTH vertex and recon (where possible)

	if (doing_phsp) then
	  main.sigcc = 1.0
	  main.sigcc_recon = 1.0
	elseif (doing_hyd_elast) then
	  main.sigcc = sigep(vertex)
	  main.sigcc_recon = sigep(recon)
	elseif (doing_deuterium.or.doing_heavy) then
	  main.sigcc = deForest(vertex)
	  main.sigcc_recon = deForest(recon)
	elseif (doing_pion) then
	  main.sigcc = peepi(vertex,main)
	  main.sigcc_recon = 1.0		!can we do this with recon??
	elseif (doing_kaon) then
	  main.sigcc = peeK(vertex,main)
	  main.sigcc_recon = 1.0		!can we do this with recon??
        else
          main.sigcc = 1.0
          main.sigcc_recon = 1.0
	endif

	if (debug(3)) then
	  write(6,*)'======================================'
	  write(6,*)'complete_main:'
	  write(6,*)'theta_e =',vertex.e.theta
	  write(6,*)'q mag =',vertex.q
	  write(6,*)'omega =',vertex.omega
	  write(6,*)'Q2 =',vertex.Q2/1000000.
	  write(6,*)'Ein =',vertex.Ein
	  write(6,*)'hadron mom =',vertex.p.P
	  write(6,*)'e mom =',vertex.e.P
	  write(6,*)'mass =',Mh
	  write(6,*)'epsilon =',main.epsilon
	  write(6,*)'phi_pq =',main.phi_pq
	  write(6,*)'theta_pq =',main.theta_pq
	  write(6,*)'======================================'
	endif

! The total contributing weight from this event -- this weight is
! proportional to # experimental counts represented by the event.

	main.weight = main.SF_weight*main.jacobian*main.gen_weight*main.sigcc
	if (debug(5))write(6,*) 'gen_weight = ',main.gen_weight,
     >		main.jacobian,main.sigcc

	success = .true.
	decdist(3)=decdist(3)+main.weight

	if (debug(2)) write(6,*)'comp_main: ending, success =',success
	return
	end


	subroutine physics_angles(theta0,phi0,dx,dy,theta,phi)

!Generate physics angles in lab frame.  Theta is angle from beamline.
!phi is projection on x-y plane (so phi=0 is down, sos=pi/2, hms=3*pi/2.
!
!theta=acos( (cos(theta0)-dy*sin(theta0)*sin(phi0))/ sqrt(1+dx**2+dy**2) )
!phi=atan( (dy*cos(theta0)+sin(theta0)*sin(phi0)) / (sin(theta0)*cos(phi0)+dx) )
!
! Note that these formulae assume phi0=pi/2 or 3*pi/2 (thus the sin(phi0)
! gives a -/+ sign for the HMS/SOS).  Thus, we set the cos(phi0) term to zero.

	real*8 dx,dy		!dx/dy (xptar/yptar) for event.
	real*8 theta0,phi0	!central physics angles of spectrometer.
	real*8 theta,phi	!physics angles for event.
	real*8 r,sinth,costh,sinph,cosph	!intermediate variables.

	include 'simulate.inc'

	costh = cos(theta0)
	sinth = sin(theta0)
	sinph = sin(phi0)
	cosph = cos(phi0)
	r = sqrt( 1. + dx**2 + dy**2 )

	if (abs(cosph).gt.0.0001) then	!phi not at +/- pi/2
	  write(6,*) 'theta,phi will be incorrect if phi0 <> pi/2 or 3*pi/2'
	  write(6,*) 'phi0=',phi0,'=',phi0*180/pi,'degrees'
	endif

	theta = acos( (costh - dy*sinth*sinph) / r )
	if (dx.ne.0.0) then
	  phi   = atan( (dy*costh + sinth*sinph) / dx )	!gives -90 to 90 deg.
	  if (phi.le.0) phi=phi+pi			!make 0 to 180 deg.
	  if (sinph.lt.0.) phi=phi+pi		!add pi to phi for HMS
	else
	  phi = phi0
	endif

	return
	end
