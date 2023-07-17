	subroutine limits_update(main,vertex,orig,recon,doing_deuterium,
     >		doing_pion,doing_kaon,doing_delta,doing_rho,contrib,slop)

	implicit none

	include 'structures.inc'
	include 'radc.inc'
	type(event_main):: main
	type(event):: vertex, orig, recon
	type(contribtype):: contrib
	type(sloptype):: slop
	integer i
	logical	doing_deuterium, doing_pion, doing_kaon, doing_delta, doing_rho

! Update the "contribution limits" records

! ... GENERATION values
	call update_range(vertex%e%delta, contrib%gen%e%delta)
	call update_range(vertex%e%yptar, contrib%gen%e%yptar)
	call update_range(vertex%e%xptar, contrib%gen%e%xptar)
	call update_range(vertex%p%delta, contrib%gen%p%delta)
	call update_range(vertex%p%yptar, contrib%gen%p%yptar)
	call update_range(vertex%p%xptar, contrib%gen%p%xptar)
	call update_range(main%Trec, contrib%gen%Trec)

! ........ another tricky shift
	if (doing_deuterium .or. doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho) then
	  call update_range(vertex%e%E-main%Ein_shift,contrib%gen%sumEgen)
	else
	  call update_range(vertex%e%E+vertex%p%E-main%Ein_shift,contrib%gen%sumEgen)
	endif

! ... TRUE values
! ........ tricky shift here, remember this'll get compared with edge.e.E,
! ........ compensate for main.target.Coulomb to
! ........ copy the shift made to the edge in generate
	call update_range(orig%e%E-main%Ee_shift, contrib%tru%e%E)
	call update_range(orig%e%xptar, contrib%tru%e%xptar)
	call update_range(orig%e%yptar, contrib%tru%e%yptar)
	call update_range(orig%e%xptar, contrib%tru%e%xptar)
	call update_range(orig%p%E, contrib%tru%p%E)
	call update_range(orig%p%yptar, contrib%tru%p%yptar)
	call update_range(orig%p%xptar, contrib%tru%p%xptar)

! ........ another tricky shift
	call update_range(orig%Em-main%Ein_shift+main%Ee_shift,contrib%tru%Em)
	call update_range(orig%Pm, contrib%tru%Pm)
	call update_range(orig%Trec, contrib%tru%Trec)

! ... SPECTROMETER values
	call update_range(main%SP%e%delta, contrib%SP%e%delta)
	call update_range(main%SP%e%yptar, contrib%SP%e%yptar)
	call update_range(main%SP%e%xptar, contrib%SP%e%xptar)
	call update_range(main%SP%p%delta, contrib%SP%p%delta)
	call update_range(main%SP%p%yptar, contrib%SP%p%yptar)
	call update_range(main%SP%p%xptar, contrib%SP%p%xptar)

! ... VERTEX values
	call update_range(vertex%Trec, contrib%vertex%Trec)
	call update_range(vertex%Em, contrib%vertex%Em)
	call update_range(vertex%Pm, contrib%vertex%Pm)

! ... RADIATION stuff
! ??? should be looking at Egamma2+3 cause we do use limits on that, indirectly

	do i = 1, 3
	  call update_range(Egamma_used(i), contrib%rad%Egamma(i))
	enddo
	call update_range(Egamma_used(1)+Egamma_used(2)+Egamma_used(3),
     >		contrib%rad%Egamma_total)

! Update the "slop limits" records
! ... MC slops
	call update_range(main%RECON%e%delta-main%SP%e%delta,slop%MC%e%delta)
	call update_range(main%RECON%e%yptar-main%SP%e%yptar,slop%MC%e%yptar)
	call update_range(main%RECON%e%xptar-main%SP%e%xptar,slop%MC%e%xptar)
	call update_range(main%RECON%p%delta-main%SP%p%delta,slop%MC%p%delta)
	call update_range(main%RECON%p%yptar-main%SP%p%yptar,slop%MC%p%yptar)
	call update_range(main%RECON%p%xptar-main%SP%p%xptar,slop%MC%p%xptar)

! %.. total slops
! ........ that tricky shift again, slops accounted for by the shift not
! ........ included in slop.total.Em.
	call update_range(recon%Em-(orig%Em-main%Ein_shift+main%Ee_shift),
     >		slop%total%Em)
	call update_range(abs(recon%Pm)-abs(orig%Pm), slop%total%Pm)

	return
	end

!-------------------------------------------------------------------

	subroutine update_range(val,range)

	include 'structures.inc'
	type(rangetype):: range
	real*8	val

	range%lo = min(range%lo, val)
	range%hi = max(range%hi, val)

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
	real*8 m_spec		!spectator (A-1) mass based on missing energy
	real*8 gauss1
	logical success
	real*8 grnd		!random # generator.
	type(event_main):: main
	type(event):: vertex, orig

	real*8 nsig_max
	parameter(nsig_max=3.0e0)      !max #/sigma for gaussian ran #s.


! Randomize the position of the interaction inside the available region.
! gen.xwid and gen.ywid are the intrinsic beam widths (one sigma value).
! Use a gaussian beam distribution with +/- 3.0 sigma (add raster afterwards).

C DJG As best I can figger, main.target.y is positive when the beam is high in the lab
C DJG and main.target.x is positive when the beam is right when looking downstream.
C DJG Don't ask me why, but it seems to be this way.
C DJG Note that this means that +fry points down. I will make frx point left. 

	main%target%x = gauss1(nsig_max)*gen%xwid+targ%xoffset
	main%target%y = gauss1(nsig_max)*gen%ywid+targ%yoffset

! fr_pattern=1 - old bedpost raster rectangle - targ.fr1/fr2 are the x/y raster half-widths.
! fr_pattern=2 - circle - targ.fr1/fr2 are the inner and outer radii.
! fr_pattern=3 - new flat raster rectangle - targ.fr1/fr2 are the x/y raster half-widths.

	if (targ%fr_pattern .eq. 1) then		!old bedpost, square raster
	  t3=grnd()*pi
	  t4=grnd()*pi
	  t5=cos(t3)*targ%fr1
	  t6=cos(t4)*targ%fr2
	elseif (targ%fr_pattern .eq. 2) then		!circular raster
	  t3=grnd()*2.*pi
	  t4=sqrt(grnd())*(targ%fr2-targ%fr1)+targ%fr1
	  t5=cos(t3)*t4
	  t6=sin(t3)*t4
	elseif (targ%fr_pattern .eq. 3) then		!new, flat square raster
	  t3=2.*grnd()-1.0
	  t4=2.*grnd()-1.0
	  t5=targ%fr1*t3
	  t6=targ%fr2*t4
	else						!no raster
	  t5=0.0
	  t6=0.0
	endif

	main%target%x = main%target%x+t5
	main%target%y = main%target%y+t6
	main%target%z = (0.5-grnd())*targ%length+targ%zoffset
	main%target%rastery = t6	!'raster' contribution to vert. pos.
	main%target%rasterx = t5	! points right as you look downstream  - need to flip sign later.

! Take fluctuations of the beam energy into account, and remember to correct
! for ionization loss in the target and Coulomb acceleration of incoming
! electron.  Remove targ.zoffset from the z position of the scattering
! in order to get the position relative to the center of the target.

	call trip_thru_target (1, main%target%z-targ%zoffset, Ebeam, 0.0e0,
     >		main%target%Eloss(1), main%target%teff(1),Me,1)
	if (.not.using_Eloss) main%target%Eloss(1) = 0.0
	if (using_Coulomb) then
c	  main%target%Coulomb=targ%Coulomb_constant*(3.-(grnd())**(2./3.))
C modified 5/15/06 for poinct
	  main%target%Coulomb=targ%Coulomb_constant
	else
	  main%target%Coulomb=0.0
	endif
	vertex%Ein = Ebeam + (grnd()-0.5)*dEbeam +
     >		main%target%Coulomb - main%target%Eloss(1)

! ... deterimine known variation in Ein from Ebeam_vertex_ave and in
! ... Ee due to Coulomb energy to compare limits to generated event by event.
! ... (USED TO SHIFT 'LIMIT' VALUES IN UPDATE_RANGE CALLS.  ARE THEY NEEDED???)

        main%Ein_shift = vertex%Ein - Ebeam_vertex_ave
        main%Ee_shift = main%target%Coulomb - targ%Coulomb%ave

! Initialize success code and fractional weight due to radiation and
! generation tricks

5	success = .false.
	main%gen_weight = 1.0

! Generated quantities: (phase_space NOT YET IMPLEMENTED).
!
! phase_space: Generate electron E,yptar,xptar and hadron yptar,xptar??
! doing_hyd_elast: Generate electron angles. Solve for everything else.
! doing_deuterium: Generate electron energy and angles, proton angles.
!	Solve for proton momentum, p_fermi.
! doing_eep, A>2: generate electron and hadron energy, angles. Solve for Em,Pm.
! doing_pion: generate electron energy and angles, hadron angles, p_fermi, Em.
!	Solve for hadron momentum.
! doing_kaon: as doing_pion.
! doing_delta: as doing_pion.
! doing_rho: as doing_pion.
! doing_semi: Generate electron E,yptar,xptar and hadron E, yptar,xptar
!
! The above is summarized in the following table:
!
!                    ELECTRON                  HADRON
!               ------------------      ------------------
!               E       yptar   xptar   E       yptar   xptar   p_fermi	Em
!
!H(e,e'p)		X	X
!D(e,e'p)	X	X	X		X	X
!A(e,e'p)	X	X	X	X	X	X
!----------------------------------------------------------------------
!H(e,e'pi)	X	X	X		X	X
!A(e,e'pi)	X	X	X		X	X	X	X
!H(e,e'K)	X	X	X		X	X
!A(e,e'K)	X	X	X		X	X	X	X
!H(e,e'p)pi	X	X	X		X	X
!A(e,e'p)pi	X	X	X		X	X	X	X
!----------------------------------------------------------------------
!phase_space	X	X	X	?	X	X
!
! So our procedure is the following:
! 1) Always generate electron yptar and xptar
! 2) generate hadron yptar and xptar for all cases except H(e,e'p).
! 3) Generate hadron E for A(e,e'p)
! 4) generate electron E for all but hydrogen elastic.
! 5) generate p_fermi, Em for A(e,e'pi) and A(e,e'K).
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
	vertex%e%yptar=gen%e%yptar%min+grnd()*(gen%e%yptar%max-gen%e%yptar%min)
	vertex%e%xptar=gen%e%xptar%min+grnd()*(gen%e%xptar%max-gen%e%xptar%min)

! Generate Hadron Angles (all but H(e,e'p)):
	if (doing_deuterium.or.doing_heavy.or.doing_pion.or.doing_kaon
     >         .or.doing_delta.or.doing_semi) then
	  vertex%p%yptar=gen%p%yptar%min+grnd()*
     >  	(gen%p%yptar%max-gen%p%yptar%min)
	  vertex%p%xptar=gen%p%xptar%min+grnd()*
     >          (gen%p%xptar%max-gen%p%xptar%min)
	endif

! Generate Hadron Momentum (A(e,e'p) or semi-inclusive production).
	if (doing_heavy .or. doing_semi) then
	  Emin = max(gen%p%E%min, gen%sumEgen%min - gen%e%E%max)
	  Emax = min(gen%p%E%max, gen%sumEgen%max - gen%e%E%min)
	  if (Emin.gt.Emax) goto 100
	  main%gen_weight=main%gen_weight*(Emax-Emin)/(gen%p%E%max-gen%p%E%min)
	  vertex%p%E = Emin + grnd()*(Emax-Emin)
	  vertex%p%P = sqrt(vertex%p%E**2 - Mh2)
	  vertex%p%delta = 100.*(vertex%p%P-spec%p%P)/spec%p%P
	endif

! Generate Electron Energy (all but hydrogen elastic)
	if (doing_deuterium.or.doing_heavy.or.doing_pion.or.doing_kaon
     >       .or.doing_delta.or.doing_rho.or.doing_semi) then
	  Emin=gen%e%E%min
	  Emax=gen%e%E%max
	  if (doing_deuterium .or. doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho) then
	    Emin = max(Emin,gen%sumEgen%min)
	    Emax = min(Emax,gen%sumEgen%max)
	  else if (doing_heavy) then		! A(e,e'p)
	    Emin = max(Emin, gen%sumEgen%min - vertex%p%E)
	    Emax = min(Emax, gen%sumEgen%max - vertex%p%E)
	  endif
	  if (Emin.gt.Emax) goto 100
	  main%gen_weight=main%gen_weight*(Emax-Emin)/(gen%e%E%max-gen%e%E%min)
	  vertex%e%E = Emin + grnd()*(Emax-Emin)
	  vertex%e%P = vertex%e%E
	  vertex%e%delta = 100.*(vertex%e%P-spec%e%P)/spec%e%P
	endif	!not (doing_hyd_elast)


! Calculate the electron and proton PHYSICS angles from the spectrometer angles.
! Note that the proton angles are not yet know for hydrogen elastic.
! NOTE: this needs to be done again for the exclusive rho stuff (just on the hadron side).

	call physics_angles(spec%e%theta,spec%e%phi,
     &		vertex%e%xptar,vertex%e%yptar,vertex%e%theta,vertex%e%phi)
	call physics_angles(spec%p%theta,spec%p%phi,
     &		vertex%p%xptar,vertex%p%yptar,vertex%p%theta,vertex%p%phi)


! Generate Fermi Momentum and Em for A(e,e'pi) and A(e,e'K). 
	pfer=0.0
	pferx=0.0
	pfery=0.0
	pferz=0.0
	vertex%Em=0.0
	efer=targ%Mtar_struck		!used for pion/kaon xsec calcs.
	if(doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon.or.
     >      doing_deutdelta.or.doing_hedelta.or.doing_deutrho.or.doing_herho
     >      .or.doing_deutsemi)then
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

	  if (doing_deutpi.or.doing_deutkaon.or.doing_deutdelta
     >        .or.doing_deutrho .or.doing_deutsemi) then !Em = binding energy
	    vertex%Em = Mp + Mn - targ%M
	    m_spec = targ%M - targ%Mtar_struck + vertex%Em != Mn(Mp) for pi+(-)
	    efer = targ%M - sqrt(m_spec**2+pfer**2)
	  endif
	  if (doing_hepi .or. doing_hekaon .or. doing_hedelta .or. doing_herho) then
	    call generate_em(pfer,vertex%Em)		!Generate Em
	    m_spec = targ%M - targ%Mtar_struck + vertex%Em != M^*_{A-1}
	    efer = targ%M - sqrt(m_spec**2+pfer**2)
	  endif
	endif

! Compute all non-generated quantities

	if (debug(5)) write(6,*)'gen: calling comp_ev with false, main, vertex'
	if (debug(3)) write(6,*)'gen: calling comp_ev with false, main, vertex'
	if (debug(3)) write(6,*)'gen: Ein, E =',vertex%Ein,vertex%e%E

	call complete_ev(main,vertex,success)


	main%sigcc = 1.0

	if (debug(2)) write(6,*)'gen: initial success =',success
	if (.not.success) goto 100

! ........ temporary storage of Trec for generated event

	main%Trec = vertex%Trec

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
     >		(vertex%Em .ge. VERTEXedge%Em%min .and.
     >		 vertex%Em .le. VERTEXedge%Em%max .and.
     >		 vertex%Pm .ge. VERTEXedge%Pm%min .and. 
     >		 vertex%Pm .le. VERTEXedge%Pm%max)
	  if (success) then
c	    do i = 1, neventfields
c	      orig.all(i) = vertex.all(i)
c	    enddo
	     orig = vertex
	  endif
	endif


C DJG need to decay the rho here before we begin transporting through the
C DJG spectrometer

	if(doing_rho) then
	   call rho_decay(orig,spec%p%P,main%epsilon,success)
	endif

	
100	if (debug(2)) write(6,*)'gen: final success =',success
	if (debug(2)) write(6,*)'gen: ending'
	return
	end

!---------------------------------------------------------------------

	subroutine complete_ev(main,vertex,success)

	implicit none
	include 'simulate.inc'

	real*8 a, b, c, r, t, QA, QB, QC, radical
	real*8 p_new_x,p_new_y,new_x_x,new_x_y,new_x_z
	real*8 targ_new_x,targ_new_y
	real*8 new_y_x,new_y_y,new_y_z,dummy
	real*8 targx,targy,targz
	real*8 px,py,pz,qx,qy,qz
	real*8 oop_x,oop_y
	real*8 krel,krelx,krely,krelz
	real*8 MM
	real*8 diffmin
	real*8 w,w2,prob,probtot,probsum(1000),mass_save(1000)
	real*8 Ehad2,E_rec
	real*8 grnd,rn		!random # generator.
	integer i

	logical success
	type(event_main):: main
	type(event)::	vertex
	logical first/.true./
!-----------------------------------------------------------------------
! Calculate everything left in the /event/ structure, given all necessary
!  GENERATION values (some set of xptar,yptar,delta for both arms and p_fermi,
!  and p_fermi, depending on the scattering process: see table in generate.f 
!
! The SINGLE element of /event/ NOT computed here is sigcc.
!
! Another small anomaly is that main.jacobian IS computed here.
! This is because all the necessary terms have to be
! computed here anyway to calculate /event/ qties.
!
!-----------------------------------------------------------------------

! Initialize

! PB: generate Delta(1232) shape using mss 1.2298, width 0.135
! PB: made width narrower 0.105 9/5/2022
! PB: distribution is truncated on low mass side at P + pi mass
c PB: 9/5/22 changed to use actual Breit-Wigner shape generated
c PB: from resmod507 in first call to semi_physics.f
	if((which_pion.eq.2 .or. which_pion.eq.3).and.first) then
	   open(unit=55,file='delta_relativistic_bw.inp')
	   probtot = 0. 
	   do i=1,1000
	      w = sqrt(1.055) + 0.6 * float(i) / 1000
	      w2 = w**2
	      read(55,'(f8.3,f12.5)') w,prob
	      mass_save(i) = w
	      probtot = probtot + prob
	      probsum(i) = probtot
	   enddo
	   do i=1,1000
	      probsum(i) = probsum(i) / probtot
	      if((i/50)*50.eq.i) write(6,'(''Delta'',i5,2f8.3)')
     >   	   i,mass_save(i),probsum(i)
	   enddo
	   first = .false.
	endif

	success = .false.
	main%jacobian = 1.0

! ... unit vector components of outgoing e,p
! ... z is DOWNSTREAM, x is DOWN and y is LEFT looking downstream.

	if (debug(4)) write(6,*)'comp_ev: at 1'
	vertex%ue%x = sin(vertex%e%theta)*cos(vertex%e%phi)
	vertex%ue%y = sin(vertex%e%theta)*sin(vertex%e%phi)
	vertex%ue%z = cos(vertex%e%theta)
	if ((.not.doing_hyd_elast).and.(.not.doing_rho)) then
	  vertex%up%x = sin(vertex%p%theta)*cos(vertex%p%phi)
	  vertex%up%y = sin(vertex%p%theta)*sin(vertex%p%phi)
	  vertex%up%z = cos(vertex%p%theta)
	endif
	if (debug(4)) write(6,*)'comp_ev: at 2'

! First finish off the e side
! Calculate scattered electron energy for hydrogen/deuterium (e,e'p)

	if (doing_hyd_elast) then
	  vertex%e%E = vertex%Ein*Mh/(Mh+vertex%Ein*(1.-vertex%ue%z))

	  if (vertex%e%E.gt.vertex%Ein) return
	  vertex%e%P = vertex%e%E
	  vertex%e%delta = (vertex%e%P - spec%e%P)*100./spec%e%P
	  if (debug(4)) write(6,*)'comp_ev: at 3'
	endif

! The q vector

	if (debug(5)) write(6,*)'comp_ev: Ein,E,uez=',vertex%Ein,vertex%e%E,vertex%ue%z

	vertex%nu = vertex%Ein - vertex%e%E
	vertex%Q2 = 2*vertex%Ein*vertex%e%E*(1.-vertex%ue%z)
	vertex%q = sqrt(vertex%Q2 + vertex%nu**2)
	vertex%xbj = vertex%Q2/2./Mp/vertex%nu

	vertex%uq%x = - vertex%e%P*vertex%ue%x / vertex%q
	vertex%uq%y = - vertex%e%P*vertex%ue%y / vertex%q
	vertex%uq%z =(vertex%Ein - vertex%e%P*vertex%ue%z)/ vertex%q
	if (abs(vertex%uq%x**2+vertex%uq%y**2+vertex%uq%z**2-1).gt.0.01)
     >		stop 'Error in q vector normalization'
	if (debug(4)) write(6,*)'comp_ev: at 5'

! Now complete the p side, along with vertex.Em, vertex.Pm, vertex.Mrec.
! NOTE: Coherent pion/kaon production (bound final state) is treated as
! hydrogen, but with targ.Mtar_struck=targ.M, targ.Mtar_rec=bound final state.

	if (doing_hyd_elast) then	!p = q
	  vertex%Em = 0.0
	  vertex%Pm = 0.0
	  vertex%Mrec = 0.0
	  vertex%up%x = vertex%uq%x
	  vertex%up%y = vertex%uq%y
	  vertex%up%z = vertex%uq%z
	  vertex%p%P = vertex%q
	  vertex%p%theta = acos(vertex%up%z)
	  if (abs(vertex%up%x/sin(vertex%p%theta)).gt.1)
     >	    write(6,*) 'cos(phi)=',vertex%up%x/sin(vertex%p%theta)
	  vertex%p%phi = atan2(vertex%up%y,vertex%up%x)
	  if (vertex%p%phi.lt.0.) vertex%p%phi=vertex%p%phi+2.*pi
	  call spectrometer_angles(spec%p%theta,spec%p%phi,
     &		vertex%p%xptar,vertex%p%yptar,vertex%p%theta,vertex%p%phi)
	  vertex%p%E = sqrt(vertex%p%P**2+Mh2)
	  vertex%p%delta = (vertex%p%P - spec%p%P)*100./spec%p%P
	  if (debug(4)) write(6,*)'comp_ev: at 6'

	elseif (doing_deuterium) then	!need Ep, and a jacobian.

	  vertex%Em = targ%Mtar_struck + targ%Mrec - targ%M	!=2.2249 MeV
	  vertex%Mrec = targ%M - targ%Mtar_struck + vertex%Em	!=targ.Mrec

	  a = -1.*vertex%q*(vertex%uq%x*vertex%up%x+vertex%uq%y*vertex%up%y+vertex%uq%z*vertex%up%z)
	  b = vertex%q**2
	  c = vertex%nu + targ%M
	  t = c**2 - b + Mh2 - vertex%Mrec**2
	  QA = 4.*(a**2 - c**2)
	  QB = 4.*c*t
	  QC = -4.*a**2*Mh2 - t**2
	  radical = QB**2 - 4.*QA*QC
	  if (radical.lt.0) return
	  vertex%p%E = (-QB - sqrt(radical))/2./QA

! Check for two solutions
!	  if ( (-QB + sqrt(radical))/2./QA .gt. Mh ) then
!	    write(6,*) 'There are two valid solutions for the hadron momentum'
!	    write(6,*) 'We always pick one, so this may be a problem, and needs to be checked'
!	    write(6,*) 'solns=',(-QB - sqrt(radical))/2./QA,(-QB + sqrt(radical))/2./QA
!	  endif

! Check for two solutions, but only print warning if both are within
! event generation limits.
	  Ehad2 = (-QB + sqrt(radical))/2./QA
	  if (Ehad2.gt.edge%p%E%min .and. Ehad2.lt.edge%p%E%max .and. ntried.le.5000) then
	    write(6,*) 'The low-momentum solution to E_hadron is within the spectrometer generation'
	    write(6,*) 'limits.  If it is in the acceptance, Its the end of the world as we know it!!!'
	    write(6,*) 'E_hadron solns=',vertex%p%E,Ehad2
	  endif


	  if (vertex%p%E.le.Mh) return
	  vertex%p%P = sqrt(vertex%p%E**2 - Mh2)
	  vertex%p%delta = (vertex%p%P - spec%p%P)*100./spec%p%P

! ........ the Jacobian here is |dEp'/dEm|
	  main%jacobian = (t*(c-vertex%p%E) + 2*c*vertex%p%E*(vertex%p%E-c)) /
     > 		(2*(a**2-c**2)*vertex%p%E + c*t)
	  main%jacobian = abs(main%jacobian)


	elseif (doing_pion .or. doing_kaon .or. doing_delta) then
	   
c	  if (doing_rho) then 
c	     Mh = Mrho
! DJG give the rho mass some width (non-relativistic Breit-Wigner)
c	     Mh = Mh + 0.5*150.2*tan((2.*grnd()-1.)*atan(2.*500./150.2))
c	     Mh2 = Mh*Mh
c	     ntup.rhomass=Mh
c            write(6,*) 'rho mass is', Mh
c	  endif
	      
C DJG If doing Deltas final state for pion production, generate Delta mass
	  if(which_pion.eq.2 .or. which_pion.eq.3) then
c factor of 0.7265 to better match data (PB)
c	     targ%Mrec_struck = Mdelta + 0.5*(0.7265)*Delta_width*tan((2.*grnd()-1.)*pi/2.)
C switch to relativistic BW for Delta
	     rn = grnd()
	     diffmin = 10000.
	     do i=1,1000
		if(abs(rn - probsum(i)).lt.diffmin) then
		   diffmin = abs(rn - probsum(i))
		   targ%Mrec_struck = mass_save(i) * 1000. ! in MeV
		endif
	     enddo
	  endif

	  vertex%Pm = pfer	!vertex%Em generated at beginning.
	  vertex%Mrec = targ%M - targ%Mtar_struck + vertex%Em
	  a = -1.*vertex%q*(vertex%uq%x*vertex%up%x+vertex%uq%y*vertex%up%y+vertex%uq%z*vertex%up%z)
	  b = vertex%q**2
	  c = vertex%nu + targ%M

! For nuclei, correct for fermi motion and missing energy.  Also, check
! second solution to quadratic equation - there are often two valid
! solutions, and we always pick the larger one (which is the forward going
! one in the center of mass) and HOPE that the smaller one is never in the
! acceptance.  If the low momentum solution IS within the acceptance, we
! have big problems.
	  if (doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon.or.
     >        doing_deutdelta.or.doing_hedelta) then
	    a = a - abs(pfer)*(pferx*vertex%up%x+pfery*vertex%up%y+pferz*vertex%up%z)
	    b = b + pfer**2 + 2*vertex%q*abs(pfer)*
     >  	 (pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z)
***	    c = c - sqrt(vertex.Mrec**2+pfer**2)
	    c = vertex%nu + efer !same as above, but this way if we redefine
				    !'efer', it's the same everywhere.
	  endif
	  t = c**2 - b + Mh2 - targ%Mrec_struck**2
	  QA = 4.*(a**2 - c**2)
	  QB = 4.*c*t
	  QC = -4.*a**2*Mh2 - t**2

!	write(6,*) '    '
!	write(6,*) '    '
!	write(6,*) 'E0=',vertex.Ein
!	write(6,*) 'P_elec,P_prot=',vertex.e.P/1000.,vertex.p.P/1000.
!	write(6,*) 'thetae,phie=',vertex.e.theta*180./pi,vertex.e.phi*180./pi
!	write(6,*) 'thetap,phip=',vertex.p.theta*180./pi,vertex.p.phi*180./pi
!	write(6,*) 'q,nu,costhetapq=',vertex.q,vertex.nu,(vertex.uq.x*vertex.up.x+vertex.uq.y*vertex.up.y+vertex.uq.z*vertex.up.z)
!	write(6,*) 'a,b,c=',a/1000.,b/1000000.,c/1000.
!	write(6,*) 't=',t/1000000.
!	write(6,*) 'A,B,C=',QA/1.d6,QB/1.d9,QC/1.d12
!	write(6,*) 'rad=',QB**2 - 4.*QA*QC
!	write(6,*) 'e1,e2=',(-QB-sqrt(radical))/2000./QA,(-QB+sqrt(radical))/2000./QA
!	write(6,*) 'E_pi1,2=',vertex.nu+targ.M-(-QB-sqrt(radical))/2./QA,
!     >				vertex.nu+targ.M-(-QB+sqrt(radical))/2./QA


	  radical = QB**2 - 4.*QA*QC
	  if (radical.lt.0) return
	  vertex%p%E = (-QB - sqrt(radical))/2./QA
	  Ehad2 = (-QB + sqrt(radical))/2./QA

	  if (doing_delta) then		!choose one of the two solutions.
!	    write(6,*) ' e1, e2=',vertex%p%E,Ehad2
	    if (grnd().gt.0.5) vertex%p%E = Ehad2
	  else				!verify that 'backwards' soln. is no good.
	    if (Ehad2.gt.edge%p%E%min .and. Ehad2.lt.edge%p%E%max .and. ntried.le.5000) then
	      write(6,*) 'The low-momentum solution to E_hadron is within the spectrometer generation'
	      write(6,*) 'limits.  If it is in the acceptance, Its the end of the world as we know it!!!'
	      write(6,*) 'E_hadron solns=',vertex%p%E,Ehad2
	    endif
	  endif

	  E_rec=c-vertex%p%E	!energy of recoil system
	  if (E_rec.le.targ%Mrec_struck) return	!non-physical solution

	  if (vertex%p%E.le.Mh) return
	  vertex%p%P = sqrt(vertex%p%E**2 - Mh2)
	  vertex%p%delta = (vertex%p%P - spec%p%P)*100./spec%p%P
!	write(6,*) 'p,e=',vertex%p%P,vertex%p%E

	elseif (doing_rho) then
	   call generate_rho(vertex,success)  !generate rho in 4pi in CM
	   if(.not.success) then
	      return
	   else  ! we have a success, but set back to false for rest of complete_ev
	      success=.false.
	   endif

	elseif (doing_phsp) then

	  vertex%p%P = spec%p%P		!????? single arm phsp??
	  vertex%p%E= sqrt(Mh2+vertex%p%P**2)
	  vertex%p%delta = (vertex%p%P - spec%p%P)*100./spec%p%P
	  if (debug(4)) write(6,*)'comp_ev: at 7.5',Mh2,vertex%p%E

	endif

	if (debug(4)) write(6,*)'comp_ev: at 7'

! Compute some pion and kaon stuff.  Some of these should be OK for proton too.


	if (doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho .or. doing_semi) then
	  W2 = targ%Mtar_struck**2 + 2.*targ%Mtar_struck*vertex%nu - vertex%Q2
	  main%W = sqrt(abs(W2)) * W2/abs(W2) 
	  main%epsilon=1./(1. + 2.*(1+vertex%nu**2/vertex%Q2)*tan(vertex%e%theta/2.)**2)
	  main%theta_pq=acos(vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y+vertex%up%z*vertex%uq%z)
	  main%t = vertex%Q2 - Mh2 + 2*vertex%nu*vertex%p%E -
     >		2*vertex%p%P*vertex%q*cos(main%theta_pq)
	  main%tmin = vertex%Q2 - Mh2 + 2*vertex%p%E*vertex%nu -
     >		2*vertex%p%P*vertex%q
	  main%q2 = vertex%Q2

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

	  qx = -vertex%uq%y		!convert to 'replay' coord. system
	  qy =  vertex%uq%x
	  qz =  vertex%uq%z
	  px = -vertex%up%y
	  py =  vertex%up%x
	  pz =  vertex%up%z

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

	  main%phi_pq = atan2(p_new_y,p_new_x)		!atan2(y,x)=atan(y/x)
	  if (main%phi_pq.lt.0.e0) main%phi_pq=main%phi_pq+2.*pi
!	  if (p_new_y.lt.0.) then
!	    main.phi_pq = 2*pi - main.phi_pq
!	  endif

!	  if ((p_new_x**2+p_new_y**2).eq.0.) then
!	    main.phi_pq = 0.0
!	  else
!	    main.phi_pq = acos(p_new_x/sqrt(p_new_x**2+p_new_y**2))
!	  endif
!	  if (p_new_y.lt.0.) then
!	    main.phi_pq = 2*pi - main.phi_pq
!	  endif

	  if (debug(2)) then
	    write(6,*)'comp_ev: nu =',vertex%nu
	    write(6,*)'comp_ev: Q2 =',vertex%Q2
	    write(6,*)'comp_ev: theta_e =',vertex%e%theta
	    write(6,*)'comp_ev: epsilon =',main%epsilon
	    write(6,*)'comp_ev: theta_pq =',main%theta_pq
	    write(6,*)'comp_ev: phi_pq =',main%phi_pq
	    write(6,*)'comp_ev: E_hadron =',vertex%p%E
	  endif


	  if(using_tgt_field) then   !calculate some azimuthal angles that only make
	                             !sense for polarized target
! CALCULATE ANGLE PHI BETWEEN SCATTERING PLANE AND TARGET POLARIZATION.
! As seen looking downstream:  		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	     qx = -vertex%uq%y	!convert to 'replay' coord. system
	     qy =  vertex%uq%x
	     qz =  vertex%uq%z

c Target in-plane, so targy=0
	     targx = -targ_pol*sin(abs(targ_Bangle)) ! replay coordinates
	     targy = 0.0
	     targz = targ_pol*cos(abs(targ_Bangle)) 

c Target out of plane, so targx=0
c       targx = 0.0      ! replay coordinates
c       targy = -targ_pol*sin(abs(targ_Bangle))
c       targz = targ_pol*cos(abs(targ_Bangle)) 
	     
	     dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	     new_x_x = -qx*qz/dummy
	     new_x_y = -qy*qz/dummy
	     new_x_z = (qx**2 + qy**2)/dummy
	     
	     dummy   = sqrt(qx**2 + qy**2)
	     new_y_x =  qy/dummy
	     new_y_y = -qx/dummy
	     new_y_z =  0.0

	     p_new_x = targx*new_x_x + targy*new_x_y + targz*new_x_z
	     p_new_y = targx*new_y_x + targy*new_y_y + targz*new_y_z

	     main%phi_targ = atan2(p_new_y,p_new_x) !atan2(y,x)=atan(y/x)
	     if(main%phi_targ.lt.0.) main%phi_targ = 2.*pi+main%phi_targ


! CALCULATE ANGLE BETA BETWEEN REACTION PLANE AND TRANSVERSE TARGET 
! POLARIZATION.
! As seen looking downstream:		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	     qx = -vertex%uq%y	!convert to 'replay' coord% system
	     qy =  vertex%uq%x
	     qz =  vertex%uq%z

C Taret in plane
	     targx = -targ_pol*sin(abs(targ_Bangle)) ! 'replay' coordinates
	     targy = 0.0
	     targz = targ_pol*cos(abs(targ_Bangle)) 

C Target out of plane

c	targx = 0.0		! 'replay' coordinates
c	targy = -targ_pol*sin(abs(targ_Bangle))
c	targz = targ_pol*cos(abs(targ_Bangle)) 

	     px = -vertex%up%y
	     py =  vertex%up%x
	     pz =  vertex%up%z
	     
	     dummy   = sqrt((qy*pz-qz*py)**2 + (qz*px-qx*pz)**2 + (qx*py-qy*px)**2)
	     new_y_x =  (qy*pz-qz*py)/dummy
	     new_y_y =  (qz*px-qx*pz)/dummy
	     new_y_z =  (qx*py-qy*px)/dummy
	     
	     dummy   = sqrt((new_y_y*qz-new_y_z*qy)**2 + (new_y_z*qx-new_y_x*qz)**2
     >  	  + (new_y_x*qy-new_y_y*qx)**2)

	     new_x_x = (new_y_y*qz - new_y_z*qy)/dummy
	     new_x_y = (new_y_z*qx - new_y_x*qz)/dummy
	     new_x_z = (new_y_x*qy - new_y_y*qx)/dummy
	     

	     targ_new_x = targx*new_x_x + targy*new_x_y + targz*new_x_z
	     targ_new_y = targx*new_y_x + targy*new_y_y + targz*new_y_z

	     main%beta = atan2(targ_new_y,targ_new_x)
	     if(main%beta .lt. 0.) main%beta = 2*pi+main%beta


CDJG Calculate the "Collins" (phi_pq+phi_targ) and "Sivers"(phi_pq-phi_targ) angles
	     vertex%phi_s = main%phi_pq-main%phi_targ
	     if(vertex%phi_s .lt. 0.) vertex%phi_s = 2*pi+vertex%phi_s

	     vertex%phi_c = main%phi_pq+main%phi_targ
	     if(vertex%phi_c .gt. 2.*pi) vertex%phi_c = vertex%phi_c-2*pi
	     if(vertex%phi_c .lt. 0.0) vertex%phi_c = 2*pi+vertex%phi_c

	     
	     dummy = sqrt((qx**2+qy**2+qz**2))*sqrt((targx**2+targy**2+targz**2))
	     main%theta_tarq = acos((qx*targx+qy*targy+qz*targz)/dummy)

	  endif !polarized-target specific azimuthal angles


	endif		!end of pion/kaon specific stuff.

	if (debug(4)) write(6,*)'comp_ev: at 8'

! Compute the Pm vector in in SIMC LAB system, with x down, and y to the left.
! Computer Parallel, Perpendicular, and Out of Plane componenants.
! Parallel is component along q_hat.  Out/plane is component along
! (e_hat) x (q_hat)  (oop_x,oop_y are components of out/plane unit vector)
! Perp. component is what's left: along (q_hat) x (oop_hat).
! So looking along q, out of plane is down, perp. is left.

	vertex%Pmx = vertex%p%P*vertex%up%x - vertex%q*vertex%uq%x
	vertex%Pmy = vertex%p%P*vertex%up%y - vertex%q*vertex%uq%y
	vertex%Pmz = vertex%p%P*vertex%up%z - vertex%q*vertex%uq%z
	vertex%Pmiss = sqrt(vertex%Pmx**2+vertex%Pmy**2+vertex%Pmz**2)
	vertex%Emiss = vertex%nu + targ%M - vertex%p%E

!STILL NEED SIGN FOR PmPer!!!!!!

	oop_x = -vertex%uq%y	! oop_x = q_y *(z_hat x y_hat) = q_y *-(x_hat)
	oop_y =  vertex%uq%x	! oop_y = q_x *(z_hat x x_hat) = q_x * (y_hat)
	vertex%PmPar = (vertex%Pmx*vertex%uq%x + vertex%Pmy*vertex%uq%y + vertex%Pmz*vertex%uq%z)
	vertex%PmOop = (vertex%Pmx*oop_x + vertex%Pmy*oop_y) / sqrt(oop_x**2+oop_y**2)
	vertex%PmPer = sqrt( max(0.e0, vertex%Pm**2 - vertex%PmPar**2 - vertex%PmOop**2 ) )

	if (debug(4)) write(6,*)'comp_ev: at 9',vertex%Pmx,vertex%Pmy,vertex%Pmz

! Calculate Em, Pm, Mrec, Trec for all cases not already done.
! For doing_heavy, get Mrec from nu+M=Ep+Erec, and Erec**2=Mrec**2+Pm**2
! For (e,e'pi/K), could go back and determine momentum of recoil struck
! particle, and get Mrec and Trec seperately for struck nucleon(hyperon)
! and A-1 system.  For now, just take Trec for the A-1 system, and ignore
! the recoiling struck nucleon (hyperon), so Trec=0 for hydrogen target.

	if (doing_hyd_elast) then
	  vertex%Trec = 0.0
	else if (doing_deuterium) then
	  vertex%Pm = vertex%Pmiss
	  vertex%Trec = sqrt(vertex%Mrec**2 + vertex%Pm**2) - vertex%Mrec
	else if (doing_heavy) then
	  vertex%Pm = vertex%Pmiss
	  vertex%Mrec = sqrt(vertex%Emiss**2-vertex%Pmiss**2)
	  vertex%Em = targ%Mtar_struck + vertex%Mrec - targ%M
	  vertex%Trec = sqrt(vertex%Mrec**2 + vertex%Pm**2) - vertex%Mrec
	else if (doing_hydpi .or. doing_hydkaon .or. doing_hyddelta .or. doing_hydrho) then
	  vertex%Trec = 0.0
	else if (doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon
     >       .or.doing_deutdelta.or.doing_hedelta.or.doing_deutrho.or.doing_herho) then
	  vertex%Trec = sqrt(vertex%Mrec**2 + vertex%Pm**2) - vertex%Mrec
	else if (doing_semi) then
	   vertex%Pm = vertex%Pmiss
	   vertex%Em = vertex%Emiss
	endif

	if (debug(5)) write(6,*) 'vertex%Pm,vertex%Trec,vertex%Em',vertex%Pm,vertex%Trec,vertex%Em
	if (debug(4)) write(6,*)'comp_ev: at 10'


! calculate krel for deuteron/heavy pion(kaon).  Deuteron is straightforward.
! A>2 case is some approximation for 3He (DJG).

	if (doing_deutpi .or. doing_deutkaon .or. doing_deutdelta .or. doing_deutrho) then
	  if ((vertex%Emiss**2-vertex%Pmiss**2).lt.0) write(6,*) 'BAD MM!!!!! Emiss,Pmiss=',vertex%Emiss, vertex%Pmiss
	  MM = sqrt(max(0.e0,vertex%Emiss**2-vertex%Pmiss**2))
	  krel = sqrt( max(0.e0,MM**2-4.*targ%Mrec_struck**2) )
	else if (doing_hepi .or. doing_hekaon .or. doing_hedelta .or. doing_herho) then
	  if ((vertex%Emiss**2-vertex%Pmiss**2).lt.0) write(6,*) 'BAD MM!!!!! Emiss,Pmiss=',vertex%Emiss, vertex%Pmiss
	  MM = sqrt(max(0.e0,vertex%Emiss**2-vertex%Pmiss**2))
	  krelx = vertex%Pmx + 1.5*pferx*pfer
	  krely = vertex%Pmy + 1.5*pfery*pfer
	  krelz = vertex%Pmz + 1.5*pferz*pfer
	  krel = sqrt(krelx**2+krely**2+krelz**2)
	  if(vertex%Em.lt.6.0) krel = -krel		!bound state test???
	endif
	ntup%krel = krel

	if(doing_semi) then
CDJG	   if ((vertex%Emiss**2-vertex%Pmiss**2).lt.0) then
CDJG I should be testing that the missing mass is above two pion
CDJG threshold! Otherwise, it's just exclusive
c	   if ((vertex%Emiss**2-vertex%Pmiss**2).lt.(Mp+Mpi0)**2) then
	   if (((targ%Mtar_struck+vertex%nu-vertex%p%E)**2-vertex%Pmiss**2).lt.(Mp+Mpi0)**2) then
	      success=.false.
	      return
	   endif
	endif

	if(doing_semi) then
	   vertex%zhad = vertex%p%E/vertex%nu
	   vertex%pt2 = vertex%p%P**2*(1.0-cos(main%theta_pq)**2)
	   if(vertex%zhad.gt.1.0) then
	      success=.false.
	      return
	   endif
	endif


! Determine PHYSICS scattering angles theta/phi for the two spectrometer
! vectors, and the Jacobian which we must include in our xsec computation
! to take into account the fact that we're not generating in the physics
! solid angle d(omega) but in 'spectrometer angles' yptar,xptar.
! NOTE: the conversion from spectrometer to physics angles was done when
! the angles were first generated (except for the proton angle for hydrogen
! elastic, where it was done earlier in complete_ev).
!
!	call physics_angles(spec.e.theta,spec.e.phi,
!     &		vertex.e.xptar,vertex.e.yptar,vertex.e.theta,vertex.e.phi)
!	call physics_angles(spec.p.theta,spec.p.phi,
!     &		vertex.p.xptar,vertex.p.yptar,vertex.p.theta,vertex.p.phi)

	r = sqrt(1.+vertex%e%yptar**2+vertex%e%xptar**2) !r=cos(theta-theta0)
	main%jacobian = main%jacobian / r**3		 !1/cos**3(theta-theta0)

C DJG Since we generate rho's in 4pi (in spherical angles) we don't need no
C DJG stinkin' Jacobian!

	if (doing_heavy .or. doing_pion .or. doing_kaon .or. 
     >      doing_delta .or. doing_semi) then
	  r = sqrt(1.+vertex%p%yptar**2+vertex%p%xptar**2)
	  main%jacobian = main%jacobian / r**3		 !1/cos**3(theta-theta0)
	endif

	if (debug(4)) write(6,*)'comp_ev: at 11'

! The effective target thickness that the scattered particles see, and the
! resulting energy losses


	call trip_thru_target (2, main%target%z-targ%zoffset, vertex%e%E,
     >		vertex%e%theta, main%target%Eloss(2), main%target%teff(2),Me,1)

	call trip_thru_target (3, main%target%z-targ%zoffset, vertex%p%E,
     >		vertex%p%theta, main%target%Eloss(3), main%target%teff(3),Mh,1)

	if (.not.using_Eloss) then
	  main%target%Eloss(2) = 0.0
	  main%target%Eloss(3) = 0.0
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
	real*8 targ_new_x,targ_new_y
	real*8 px,py,pz,qx,qy,qz
	real*8 targx,targy,targz
	real*8 W2
	real*8 oop_x,oop_y
	real*8 mm,mmA,mm2,mmA2,t

	logical success
	type(event)::	recon

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

	recon%Ein = Ebeam_vertex_ave	!lowered by most probable eloss (init.f)

! ... unit vector components of outgoing e,p
! ... z is DOWNSTREAM, x is DOWN and y is LEFT looking downstream.

C Ebeam_vertex_ave has been shifted to account for Coulomb effects
C In general, the Hall C analyzer does not correct for this
C If you do NOT apply an energy shift in the ENGINE to account
C for Coulomb corrections, make sure the line below is NOT commented out.
	recon%Ein = Ebeam_vertex_ave - targ%Coulomb%ave

	if (debug(4)) write(6,*)'comp_rec_ev: at 1'
	recon%ue%x = sin(recon%e%theta)*cos(recon%e%phi)
	recon%ue%y = sin(recon%e%theta)*sin(recon%e%phi)
	recon%ue%z = cos(recon%e%theta)
	recon%up%x = sin(recon%p%theta)*cos(recon%p%phi)
	recon%up%y = sin(recon%p%theta)*sin(recon%p%phi)
	recon%up%z = cos(recon%p%theta)
	if (debug(4)) write(6,*)'comp_rec_ev: at 2'

! The q vector

	if (debug(5)) write(6,*)'comp_rec_ev: Ein,E,uez=',recon%Ein,recon%e%E,recon%ue%z
	recon%nu = recon%Ein - recon%e%E
	recon%Q2 = 2*recon%Ein*recon%e%E*(1-recon%ue%z)
	recon%q	= sqrt(recon%Q2 + recon%nu**2)
	recon%uq%x = - recon%e%P*recon%ue%x / recon%q
	recon%uq%y = - recon%e%P*recon%ue%y / recon%q
	recon%uq%z =(recon%Ein - recon%e%P*recon%ue%z)/ recon%q

c	if (doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho .or. doing_semi) then
c	   W2 = targ%mtar_struck**2 + 2.*targ%mtar_struck*recon%nu - recon%Q2
c	else
c	   W2 = targ%M**2 + 2.*targ%M*recon%nu - recon%Q2
c	endif

c Everyone else in the world calculates W using the proton mass.
	W2 = mp**2 + 2.*mp*recon%nu - recon%Q2

	recon%W = sqrt(abs(W2)) * W2/abs(W2) 
	recon%xbj = recon%Q2/2./Mp/recon%nu
	if (debug(4)) write(6,*)'comp_rec_ev: at 5'

! Now complete the p side

	if(doing_phsp)then
	  recon%p%P=spec%p%P
	  recon%p%E=sqrt(Mh2+recon%p%P**2)
	  if (debug(4)) write(6,*)'comp_rec_ev: at 7.5',Mh2,recon%p%E
	endif

! Compute some pion/kaon stuff.

	recon%epsilon=1./(1. + 2.*(1+recon%nu**2/recon%Q2)*tan(recon%e%theta/2.)**2)
	recon%theta_pq=acos(min(1.0,recon%up%x*recon%uq%x+recon%up%y*recon%uq%y+recon%up%z*recon%uq%z))

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

	qx = -recon%uq%y		!convert to 'replay' coord. system
	qy =  recon%uq%x
	qz =  recon%uq%z
	px = -recon%up%y
	py =  recon%up%x
	pz =  recon%up%z

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
	  recon%phi_pq = 0.0
	else
	  recon%phi_pq = acos(p_new_x/sqrt(p_new_x**2+p_new_y**2))
	endif
	if (p_new_y.lt.0.) then
	  recon%phi_pq = 2*pi - recon%phi_pq
	endif


	if(using_tgt_field) then !calculate polarized-target specific azimuthal angles

! CALCULATE ANGLE PHI BETWEEN SCATTERING PLANE AND TARGET POL.
! Take into account the different definitions of x, y and z in
! replay and SIMC:
! As seen looking downstream:		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	   qx = -recon%uq%y	!convert to 'replay' coord. system
	   qy =  recon%uq%x
	   qz =  recon%uq%z

C In-plane target
	   targx = -targ_pol*sin(abs(targ_Bangle)) ! 'replay' coordinates
	   targy = 0.0
	   targz = targ_pol*cos(abs(targ_Bangle)) 

C out of plane target
c	targx = 0.0       ! 'replay' coordinates
c	targy = -targ_pol*sin(abs(targ_Bangle)) 
c	targz = targ_pol*cos(abs(targ_Bangle)) 

	   dummy=sqrt((qx**2+qy**2)*(qx**2+qy**2+qz**2))
	   new_x_x = -qx*qz/dummy
	   new_x_y = -qy*qz/dummy
	   new_x_z = (qx**2 + qy**2)/dummy

	   dummy   = sqrt(qx**2 + qy**2)
	   new_y_x =  qy/dummy
	   new_y_y = -qx/dummy
	   new_y_z =  0.0

	   p_new_x = targx*new_x_x + targy*new_x_y + targz*new_x_z
	   p_new_y = targx*new_y_x + targy*new_y_y + targz*new_y_z
	   

	   recon%phi_targ = atan2(p_new_y,p_new_x)
	   if(recon%phi_targ.lt.0.) recon%phi_targ = 2.*pi + recon%phi_targ


! CALCULATE ANGLE BETA BETWEEN REACTION PLANE AND (TRANSVERSE) TARGET 
! POLARIZATION.
! Also take into account the different definitions of x, y and z in
! replay and SIMC:
! As seen looking downstream:		replay	SIMC	(old_simc)
!				x	right	down	(left)
!				y	down	left	(up)
!				z	all have z pointing downstream
!
! SO: x_replay=-y_simc, y_replay=x_simc, z_replay= z_simc

	   qx = -recon%uq%y	!convert to 'replay' coord. system
	   qy =  recon%uq%x
	   qz =  recon%uq%z

	   targx = -targ_pol*sin(abs(targ_Bangle)) ! 'replay' coordinates
	   targy = 0.0
	   targz = targ_pol*cos(abs(targ_Bangle)) 

c	targx = 0.0              ! 'replay' coordinates
c	targy = -targ_pol*sin(abs(targ_Bangle))
c	targz = targ_pol*cos(abs(targ_Bangle)) 


	   px = -recon%up%y
	   py =  recon%up%x
	   pz =  recon%up%z

	   dummy   = sqrt((qy*pz-qz*py)**2 + (qz*px-qx*pz)**2 + (qx*py-qy*px)**2)
	   new_y_x =  (qy*pz-qz*py)/dummy
	   new_y_y =  (qz*px-qx*pz)/dummy
	   new_y_z =  (qx*py-qy*px)/dummy

	   dummy   = sqrt((new_y_y*qz-new_y_z*qy)**2 + (new_y_z*qx-new_y_x*qz)**2
     >  	+ (new_y_x*qy-new_y_y*qx)**2)

	   new_x_x = (new_y_y*qz - new_y_z*qy)/dummy
	   new_x_y = (new_y_z*qx - new_y_x*qz)/dummy
	   new_x_z = (new_y_x*qy - new_y_y*qx)/dummy


	   targ_new_x = targx*new_x_x + targy*new_x_y + targz*new_x_z
	   targ_new_y = targx*new_y_x + targy*new_y_y + targz*new_y_z

	   recon%beta = atan2(targ_new_y,targ_new_x)
	   if(recon%beta.lt.0.) recon%beta = 2.*pi + recon%beta

CDJG Calculate the "Collins" (phi_pq+phi_targ) and "Sivers"(phi_pq-phi_targ) angles
	   recon%phi_s = recon%phi_pq-recon%phi_targ
	   if(recon%phi_s .lt. 0.) recon%phi_s = 2*pi+recon%phi_s

	   recon%phi_c = recon%phi_pq+recon%phi_targ
	   if(recon%phi_c .gt. 2.*pi) recon%phi_c = recon%phi_c-2*pi
	   if(recon%phi_c .lt. 0.) recon%phi_c = 2*pi+recon%phi_c

	   dummy = sqrt((qx**2+qy**2+qz**2))*sqrt((targx**2+targy**2+targz**2))
	   recon%theta_tarq = acos((qx*targx+qy*targy+qz*targz)/dummy)
	endif !polarized-target specific stuff

	if (debug(2)) then
	  write(6,*)'comp_rec_ev: nu =',recon%nu
	  write(6,*)'comp_rec_ev: Q2 =',recon%Q2
	  write(6,*)'comp_rec_ev: theta_e =',recon%e%theta
	  write(6,*)'comp_rec_ev: epsilon =',recon%epsilon
	  write(6,*)'comp_rec_ev: theta_pq =',recon%theta_pq
	  write(6,*)'comp_rec_ev: phi_pq =',recon%phi_pq
	endif

	if (debug(4)) write(6,*)'comp_rec_ev: at 8'

! Compute the Pm vector in in SIMC LAB system, with x down, and y to the left.
! Computer Parallel, Perpendicular, and Out of Plane componenants.
! Parallel is component along q_hat.  Out/plane is component along
! (e_hat) x (q_hat)  (oop_x,oop_y are components of out/plane unit vector)
! Perp. component is what's left: along (q_hat) x (oop_hat).
! So looking along q, out of plane is down, perp. is left.

	recon%Pmx = recon%p%P*recon%up%x - recon%q*recon%uq%x
	recon%Pmy = recon%p%P*recon%up%y - recon%q*recon%uq%y
	recon%Pmz = recon%p%P*recon%up%z - recon%q*recon%uq%z
	recon%Pm = sqrt(recon%Pmx**2+recon%Pmy**2+recon%Pmz**2)

!STILL NEED SIGN FOR PmPer!!!!!!

	oop_x = -recon%uq%y	! oop_x = q_y *(z_hat x y_hat) = q_y *-(x_hat)
	oop_y =  recon%uq%x	! oop_y = q_x *(z_hat x x_hat) = q_x * (y_hat)
	recon%PmPar = (recon%Pmx*recon%uq%x + recon%Pmy*recon%uq%y + recon%Pmz*recon%uq%z)
	recon%PmOop = (recon%Pmx*oop_x + recon%Pmy*oop_y) / sqrt(oop_x**2+oop_y**2)
	recon%PmPer = sqrt( max(0.e0, recon%Pm**2 - recon%PmPar**2 - recon%PmOop**2 ) )

	if (doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho .or. doing_semi) then
	  recon%Em = recon%nu + targ%Mtar_struck - recon%p%E
          mm2 = recon%Em**2 - recon%Pm**2
          mm  = sqrt(abs(mm2)) * abs(mm2)/mm2
          mmA2= (recon%nu + targ%M - recon%p%E)**2 - recon%Pm**2
          mmA = sqrt(abs(mmA2)) * abs(mmA2)/mmA2
          t = recon%Q2 - Mh2
     >      + 2*(recon%nu*recon%p%E - recon%p%P*recon%q*cos(recon%theta_pq))
	  ntup%mm = mm
	  ntup%mmA = mmA
	  ntup%t = t
	endif

	if(doing_semi.or.doing_rho) then
	   recon%zhad = recon%p%E/recon%nu
	   recon%pt2 = recon%p%P**2*(1.0-cos(recon%theta_pq)**2)
	endif

! Calculate Trec, Em. Trec for (A-1) system (eep), or for struck nucleon (pi/K)
! Note that there are other ways to calculate 'Em' for the pion/kaon case.
! This Em for pi/kaon gives roughly E_m + E_rec + T_{A-1}  (I think)
	if (doing_hyd_elast) then
	  recon%Trec = 0.0
	  recon%Em = recon%nu + targ%M - recon%p%E - recon%Trec
	else if (doing_deuterium .or. doing_heavy) then
	  recon%Trec = sqrt(recon%Pm**2+targ%Mrec**2) - targ%Mrec
	  recon%Em = recon%nu + targ%Mtar_struck - recon%p%E - recon%Trec
	else if (doing_pion .or. doing_kaon .or. doing_delta .or. doing_rho) then
	  recon%Em = recon%nu + targ%Mtar_struck - recon%p%E
	endif

	if (debug(5)) write(6,*) 'recon%Pm,recon%Trec,recon%Em',recon%Pm,recon%Trec,recon%Em
	if (debug(4)) write(6,*)'comp_rec_ev: at 10'

	success=.true.
	return
	end

!------------------------------------------------------------------------

	subroutine complete_main(force_sigcc,main,vertex,vertex0,recon,success)

	implicit none
	include 'simulate.inc'

	integer		i, iPm1
	real*8		a, b, r, frac, peepi, peeK, peedelta, peerho, peepiX
	real*8		survivalprob, semi_dilution
	real*8		weight, width, sigep, deForest, tgtweight
	logical		force_sigcc, success
	type(event_main):: main
	type(event)::	vertex, vertex0, recon

!-----------------------------------------------------------------------
! Calculate everything left in the /main/ structure that hasn't been
! required up til now. This routine is called ONLY after we
! know an event IS going to contribute, don't want to be doing these
! calculations if they're not needed!
!
! A small anomaly is that sigcc from the /event/ structure is computed
! here since it's not needed til now, and was not calculated by
! COMPLETE_ev.
!-----------------------------------------------------------------------


! Initialize success code

	if (debug(2)) write(6,*)'comp_main:  entering...'
	success = .false.

! The spectral function weighting

	if (doing_hyd_elast.or.doing_pion.or.doing_kaon.or.doing_delta.or.doing_phsp.or.doing_rho.or.doing_semi) then !no SF.
	  main%SF_weight=1.0
	else if (use_benhar_sf.and.doing_heavy) then ! Doing Spectral Functions
	   call sf_lookup_diff(vertex%Em, vertex%Pm, weight)
	   main%SF_weight = targ%Z*transparency*weight
 	else if (doing_deuterium .or. (doing_heavy.and.(.not.use_benhar_sf))) then
	  main%SF_weight = 0.0
	  do i=1,nrhoPm
	    weight = 0.0

! ... use linear interpolation to determine rho(pm)

	    r = (vertex%Pm-Pm_theory(i)%min)/Pm_theory(i)%bin
	    if (r.ge.0 .and. r.le.Pm_theory(i)%n) then
	      iPm1 = nint(r)
	      if (iPm1.eq.0) iPm1 = 1
	      if (iPm1.eq.Pm_theory(i)%n) iPm1 = Pm_theory(i)%n-1
	      frac = r+0.5-float(iPm1)
	      b = theory(i,iPm1)
	      a = theory(i,iPm1+1)-b
	      weight = a*frac + b
	    endif

	    if (doing_heavy) then

	      width=Emsig_theory(i)/2.0
	      if (vertex%Em.lt.E_Fermi) weight = 0.0
	      weight=weight/pi/Em_int_theory(i) * width/
     >				((vertex%Em-Em_theory(i))**2+width**2)
	    endif
	    main%SF_weight = main%SF_weight+weight*nprot_theory(i)
	  enddo ! <i=1,nrhoPm>
	endif

! ... if we have come up with weight<=, quit now (avoid weight<0 in ntuple)
! ... unless FORCE_SIGCC is on (to get central sigcc).

	if (debug(5))write(6,*)'comp_main: calculating cross section'
	if (main%SF_weight.le.0 .and. .not.force_sigcc) return

! ... Cross section, for BOTH vertex and recon (where possible)
! ... Note that all of these are cross section per nucleon.  For A(e,e'p)
! ... this becomes cross section per nucleus because of the weighting of
! ... the spectral function.  For pion/kaon production, we explicitly
! ... apply a weight for the number of nucleons involved (protons for pi+ or Y0,
! ... neutrons for pi- or Y-) (Y=hyperon)

	tgtweight = 1.0

	if (doing_phsp) then
	  main%sigcc = 1.0
	  main%sigcc_recon = 1.0

	elseif (doing_hyd_elast) then
	  main%sigcc = sigep(vertex)
	  main%sigcc_recon = sigep(recon)

	elseif (doing_deuterium.or.doing_heavy) then
	  main%sigcc = deForest(vertex)		
	  main%sigcc_recon = deForest(recon)

	elseif (doing_pion) then
	  main%sigcc = peepi(vertex,main)
C Use Clebsch-Gordon coefficients to approximate xsec for Delta final states
C This ignores the fact that the g*p and g*n cross sections may not be the same
C 6/24/2021: Coefficients for Delta final states updated from Peter Bosted's
C empirical check's.
	  if(which_pion.eq.2) then ! pi+ Delta
	     if(doing_hydpi) then
c		main%sigcc = main%sigcc/4.0 !(pi+ Delta0)/(pi+ n)
c		main%sigcc = 0.6*main%sigcc !(pi+ Delta0)/(pi+ n)
		main%sigcc = 0.4*main%sigcc !(pi+ Delta0)/(pi+ n) updated 17july2023
	     elseif(doing_deutpi) then
c		main%sigcc = main%sigcc/4.0 !(pi+ Delta0)/pi+ n)
c     >                      + 0.75*main%sigcc !(pi+ Delta-)/(pi+ n)
c		main%sigcc = 0.6*main%sigcc !(pi+ Delta0)/pi+ n)
c     >                      + 1.0*main%sigcc !(pi+ Delta-)/(pi+ n)
		main%sigcc = 0.4*main%sigcc !(pi+ Delta0)/pi+ n)   updated 17july2023
     >                      + 0.8*main%sigcc !(pi+ Delta-)/(pi+ n)
	     endif 
	  elseif (which_pion.eq.3) then  !pi- Delta
	     if(doing_hydpi) then
c		main%sigcc = 0.6*main%sigcc ! (pi- Delta++)/(pi- p)
		main%sigcc = 0.55*main%sigcc ! (pi- Delta++)/(pi- p)  updated 17july2023
	     elseif(doing_deutpi) then
c		main%sigcc = 3.0*main%sigcc/5.0 ! (pi- Delta++)/(pi- p)
c     >                     + 0.25*main%sigcc !(pi- Delta+)/(pi- p)
c		main%sigcc = 0.6*main%sigcc ! (pi- Delta++)/(pi- p)
c     >                     + 0.6*main%sigcc !(pi- Delta+)/(pi- p)
		main%sigcc = 0.55*main%sigcc ! (pi- Delta++)/(pi- p)  updated 17july2023
     >                     + 0.99*main%sigcc !(pi- Delta+)/(pi- p)
	     endif
	  endif
	  main%sigcc_recon = 1.0
	  if (which_pion.eq.1 .or. which_pion.eq.11) then  !OK for coherent???
	    tgtweight = targ%N
	  else
	    tgtweight = targ%Z
	  endif

	elseif (doing_kaon) then
	  main%sigcc = peeK(vertex,main,survivalprob)
	  main%sigcc_recon = 1.0
	  if (which_kaon.eq.2 .or. which_kaon.eq.12) then  !OK for coherent???
	    tgtweight = targ%N
	  else
	    tgtweight = targ%Z
	  endif

	elseif (doing_delta) then
	  main%sigcc = peedelta(vertex,main)	!Need new xsec model.
	  main%sigcc_recon = 1.0

	elseif (doing_rho) then
	  main%sigcc = peerho(vertex,main)
	  main%sigcc_recon = 1.0
	  tgtweight = targ%Z+targ%N

	elseif (doing_semi) then
	  main%sigcc = peepiX(vertex,vertex0,main,survivalprob,.FALSE.)
	  main%sigcc_recon = 1.0
	  main%sigcent = peepiX(vertex,vertex0,main,survivalprob,.TRUE.)
c	  ntup%dilu = semi_dilution(vertex,main) 
	  ntup%dilu = 1.0

	else
	  main%sigcc = 1.0
	  main%sigcc_recon = 1.0
	endif

C If using Coulomb corrections, include focusing factor
	if(using_Coulomb) then
	   main%sigcc = main%sigcc*(1.0+targ%Coulomb%ave/Ebeam)**2
	endif

	if (debug(3)) then
	  write(6,*)'======================================'
	  write(6,*)'complete_main:'
	  write(6,*)'theta_e =',vertex%e%theta
	  write(6,*)'q mag =',vertex%q
	  write(6,*)'nu =',vertex%nu
	  write(6,*)'Q2 =',vertex%Q2/1000000.
	  write(6,*)'Ein =',vertex%Ein
	  write(6,*)'hadron mom =',vertex%p%P
	  write(6,*)'e mom =',vertex%e%P
	  write(6,*)'mass =',Mh
	  write(6,*)'epsilon =',main%epsilon
	  write(6,*)'phi_pq =',main%phi_pq
	  write(6,*)'theta_pq =',main%theta_pq
	  write(6,*)'======================================'
	endif

! The total contributing weight from this event -- this weight is
! proportional to # experimental counts represented by the event.
! Apply survival probability to kaons if we're not modeling decay.

	main%weight = main%SF_weight*main%jacobian*main%gen_weight*main%sigcc
	main%weight = main%weight * tgtweight	!correct for #/nucleons involved
	if ((doing_kaon.or.doing_semika) .and. .not.doing_decay) 
     >		main%weight = main%weight*survivalprob
	if (debug(5))write(6,*) 'gen_weight = ',main%gen_weight,
     >		main%jacobian,main%sigcc

	success = .true.

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
	real*8 tmp

	include 'constants.inc'

	costh = cos(theta0)
	sinth = sin(theta0)
	sinph = sin(phi0)
	cosph = cos(phi0)
	r = sqrt( 1. + dx**2 + dy**2 )

	if (abs(cosph).gt.0.0001) then	!phi not at +/- pi/2
	  write(6,*) 'theta,phi will be incorrect if phi0 <> pi/2 or 3*pi/2'
	  write(6,*) 'phi0=',phi0,'=',phi0*180/pi,'degrees'
	endif

	tmp=(costh - dy*sinth*sinph) / r
	if (abs(tmp).gt.1) write(6,*) 'tmp=',tmp
	theta = acos( (costh - dy*sinth*sinph) / r )
	if (dx.ne.0.0) then
	  phi = atan( (dy*costh + sinth*sinph) / dx )	!gives -90 to 90 deg.
	  if (phi.le.0) phi=phi+pi			!make 0 to 180 deg.
	  if (sinph.lt.0.) phi=phi+pi		!add pi to phi for HMS
	else
	  phi = phi0
	endif

	return
	end



	subroutine spectrometer_angles(theta0,phi0,dx,dy,theta,phi)

!Generate spectrometer angles from physics angles in lab frame.
!Theta is angle from beamline.
!phi is projection on x-y plane (so phi=0 is down, sos=pi/2, hms=3*pi/2.

	real*8 dx,dy		!dx/dy (xptar/yptar) for event.
	real*8 theta0,phi0	!central physics angles of spectrometer.
	real*8 theta,phi	!physics angles for event.
	real*8 x,y,z				!intermediate variables.
	real*8 x0,y0,z0				!intermediate variables.
	real*8 cos_dtheta,y_event

	include 'constants.inc'

	x = sin(theta)*cos(phi)
	y = sin(theta)*sin(phi)
	z = cos(theta)
	x0 = sin(theta0)*cos(phi0)
	y0 = sin(theta0)*sin(phi0)
	z0 = cos(theta0)

	cos_dtheta = x*x0 + y*y0 + z*z0
	dx = x / cos_dtheta
	dy = sqrt(1/cos_dtheta**2-1.-dx**2)

	y_event = y/cos_dtheta	!projected to plane perp. to spectrometer.
	if (y_event .lt. y0) dy = -dy

	return
	end
