	program simc

! Last modified:
!
!	This program represents variations on themes by
!		N. C. R. Makins (Argonne National Lab),
!		T. G. O'Neill (Argonne National Lab), and
!		Seemingly Countless Others (Virtually Everywhere).
!
	implicit none
	include 'simulate_init.inc'
	include 'histograms_init.inc'
	include 'radc.inc'
	include 'hbook.inc'
	include 'struct_sos.inc'
	include 'struct_hms.inc'

	integer*4	i, ninform
	real*8		r
	real*8		domega_e, domega_p !populated e/hadron solid angles.
	logical		success
	character	filename*80, genfile*80, histfile*80, timestring1*23
	character	timestring2*23,genifile*80
	record /event/		vertex, orig, recon
	record /event_main/	main
	record /contrib/	contrib
	record /histograms/	H
	record /event_central/	central

	real*8 one
	parameter (one=1.0d0)	!double precision 1 for subroutine calls

	real*8 grnd

	if(debug(5)) then
	  write(6,*)'sim: at start, vertex.e.yptar,main.SP.e.yptar =',vertex.e.yptar,main.SP.e.yptar
	  write(6,*)'sim: at start, cuts.Em.min, max =',cuts.Em.min,cuts.Em.max
	  write(6,*)'sim: at start, VERTEXedge.Em.min, max =',VERTEXedge.Em.min,VERTEXedge.Em.max
	endif

! INITIALIZE
!-----------

! ... initialize histogram area for HBOOK

	call hlimit(PawSize)

! ... read in the data base

	call dbase_read(H)

	if (debug(3)) write(6,*) 'Main after dbrd: p,e th',
     >	    spec.p.theta,spec.e.theta,using_P_arm_montecarlo,using_E_arm_montecarlo

! ... open diagnostic ntuple file, if requested

	i = index(base,' ')
	if (Nntu.gt.0) then
	  filename = 'worksim/'//base(1:i-1)//'.rzdat'
	  call NtupleInit(filename)
	endif

! ... set up some quantities that the radiative corrections routines will need
! ... N.B. Call this even if we're not using radiative corrections -- the
! ... target thicknesses seen by the incoming and
! ... scattered particles (in radiation lengths) are determined here, and
! ... are needed for the multiple scattering calculations
! ... done in the Monte Carlos.

	if (debug(2)) write(6,*)'sim: calling radc_init'
	call radc_init
	if (debug(2)) write(6,*)'sim: done with radc_init'

! ... compute some quantities for a central event

	call calculate_central(central)
	if (debug(2)) write(6,*)'sim: done with calc_central'
	if (debug(2)) write(6,*)'central.sigcc=',central.sigcc
	if (debug(4)) write(6,*)'sim: at 1'

	targetfac=targ.mass_amu/3.75914d+6/(targ.abundancy/100.)
     >		*abs(cos(targ.angle))/(targ.thick*1000.)
	if (debug(4)) write(6,*)'sim: at 2'

! ... Note: EXPER.charge is in mC and the luminosity comes out in fm^-2
! ...   now luminosity is in ub^-1

	luminosity = EXPER.charge/targetfac
	if (debug(4)) write(6,*)'sim: at 3'

! ... initialize the random number generator and number of attempts

	r = grnd()
	ntried = 0

! GAW - insert calls to initialize target field for both arms

	call trgInit(tgt_field_file,(targ.Bangle-spec.e.theta)*degrad,0.,
     >                              -(targ.Bangle+spec.p.theta)*degrad,0.)


! LOOP OVER SIMULATED EVENTS
!---------------------------
	write(6,*) ' '
	call date (timestring1)
	call time (timestring1(11:23))
	nevent = 0
	ninform = 1000
	decdist(3) = 0.0	!sum of main.weight (for generating ave sigma)

	do while (nevent.lt.abs(ngen))

! reset pion/kaon mass
	  decdist(30) = 0.0	!decay distance.
	  Mh2_final = Mh2

! Setup for this event

! ... keep the human interested

	  ntried = ntried + 1

	  if(debug(2)) then
     	    write(6,*)'********************************************'
	    write(6,*)'SIM:  ntried =',ntried
	    write(6,*)'      ncontribute =',ncontribute
	  endif

	  if (mod(ntried,ninform).eq.1) then
	    write (6,'(1x,a,i8,a,i8,a,e11.4)') 'Generating Event ',
     >		ntried, ' ... ', nevent,' successes so far -- Monitor:',
     >		wtcontribute*luminosity/ntried
	    if (ntried.ge.5000) ninform = 20000
	  endif

! ... generate an event

	  if(debug(3)) write(6,*)'sim: before gen: orig.p.E =',orig.p.E
	  call generate(main,vertex,orig,success)
	  if(debug(2)) write(6,*)'sim: after gen, success =',success

! Run the event through various manipulations, checking to see whether
! it made it at each stage

! ... run through spectrometers and check SP cuts

	  if(debug(3)) write(6,*)'sim: before mc: orig.p.E =',orig.p.E
	  if(success)call montecarlo(orig.e,orig.p,main,recon.e,recon.p,success)
	  if(debug(2)) write(6,*)'sim: after mc, success =',success

! ... calculate everything else about the event

	  if (success) then

! ........ recon event needs a 'measured' beam energy, which is Ebeam
! ........ corrected for average losses

	    recon.Ein = Ebeam_vertex_ave
	    call complete_recon_ev(recon,success)
	  endif
	  if(debug(2)) write(6,*)'sim: after comp_ev, success =',success
	  if(debug(5)) write(6,*) 'recon.Em,recon.Pm',recon.Em,recon.Pm

! ... calculate remaining pieces of the main structure

	  if (success) call complete_main(.false.,main,vertex,recon,success)

! We have a success! We must have passed at least one set of cuts and have
! a non-zero event weight.  IT WOULD BE USEFUL TO ADD A SIMILAR HISTOGRAM
! FOR P_FERMI (MAYBE PFERX,PFERY, AND PFERZ).

	  if(ntried.gt.0)then
	    call inc(H.geni.e.delta, vertex.e.delta, one)
	    call inc(H.geni.e.yptar, vertex.e.yptar, one)
	    call inc(H.geni.e.xptar, -vertex.e.xptar, one)
	    call inc(H.geni.p.delta, vertex.p.delta, one)
	    call inc(H.geni.p.yptar, vertex.p.yptar, one)
	    call inc(H.geni.p.xptar, -vertex.p.xptar, one)
	  endif

	  if (success) then
	    if (debug(2)) write(6,*)'sim: We have a success!'

! ... increment the "experimental" spectrometer arm and "contribution"
! ... histograms
	    call inc(H.RECON.e.delta, main.RECON.e.delta, main.weight)
	    call inc(H.RECON.e.yptar, main.RECON.e.yptar, main.weight)
	    call inc(H.RECON.e.xptar, main.RECON.e.xptar, main.weight)
	    call inc(H.RECON.p.delta, main.RECON.p.delta, main.weight)
	    call inc(H.RECON.p.yptar, main.RECON.p.yptar, main.weight)
	    call inc(H.RECON.p.xptar, main.RECON.p.xptar, main.weight)
	    call inc(H.gen.e.delta, vertex.e.delta, one)
	    call inc(H.gen.e.yptar, vertex.e.yptar, one)
	    call inc(H.gen.e.xptar, -vertex.e.xptar, one)
	    call inc(H.gen.p.delta, vertex.p.delta, one)
	    call inc(H.gen.p.yptar, vertex.p.yptar, one)
	    call inc(H.gen.p.xptar, -vertex.p.xptar, one)

! ... update counters and integrated weights FOR EVENTS INSIDE OF DELTA
!     population limits (events are only generated to fill region inside
!     of the given delta limits) and below cuts.Em.max.  Have to remove slop
!     that is added in init.f from the delta limits.
	    if (debug(4)) write(6,*)'sim: cut'

	    ncontribute = ncontribute + 1
	    if (debug(4)) write(6,*)'sim: ncontribute =',ncontribute
	    if (recon.e.delta.ge.(SPedge.e.delta.min+slop.MC.e.delta.used) .and.
     >	        recon.e.delta.le.(SPedge.e.delta.max-slop.MC.e.delta.used) .and.
     >	        recon.p.delta.ge.(SPedge.p.delta.min+slop.MC.p.delta.used) .and.
     >	        recon.p.delta.le.(SPedge.e.delta.max-slop.MC.p.delta.used) .and.
     >		recon.Em.le.cuts.Em.max) then
	      wtcontribute = wtcontribute + main.weight
	    endif
	    if (.not.rad_proton_this_ev) ncontribute_no_rad_proton =
     >			ncontribute_no_rad_proton + 1

! ... update the "contribution" and "slop" limits

	    call limits_update(main,vertex,orig,recon,doing_deuterium,
     >		doing_pion,doing_kaon,contrib,slop)

! ... write out line to diagnostic ntuple file, if requested

	  endif ! <success>

	  if (Nntu.gt.0) call results_ntu_write(main,vertex,orig,recon,success)

! END of LOOP over events
! Cute thing here, if ngen < 0, then just make |ngen| ATTEMPTS rather than
! insisting on ngen good events ... this is suitable for tests with new
! and uncertain parameters, you don't want the program fishing around all
! day for an event!

	  if (ngen.lt.0) then
	    nevent = nevent+1
	  else
	    if (success) nevent = nevent+1
	  endif

	enddo ! <loop over ntried>


! END GAME: NORMALIZE, GENERATE OUTPUT
!-------------------------------------

! Our last event announcement ... and how long did it all take?
	write (6,'(1x,''---> Last Event '',i8,'' ... '',i8,''successes'')') ntried, nevent
	call date (timestring2)
	call time (timestring2(11:23))

! NORMALIZE!

! ... put in the luminosity and efficiency factors

	normfac = luminosity/ntried*nevent

! ... multiply in the relevant phase spaces (see event.f for description
! ... of generated variables.  Electron angles for all cases.
! ... add hadron angles and electron energy for all but H and D (e,e'p).
! ... add hadron energy for A(e,e'p) (what about phase_space?)
! ... Cross sections are all in microbarn/MeV**i/sr**j (i,j depend on reaction)

	domega_e=(gen.e.yptar.max-gen.e.yptar.min)*(gen.e.xptar.max-gen.e.xptar.min)
	domega_p=(gen.p.yptar.max-gen.p.yptar.min)*(gen.p.xptar.max-gen.p.xptar.min)

	genvol = domega_e
!	genvol_inclusive = genvol	!may want dOmega, or dE*dOmega

	if (.not.doing_hyd_elast) then	!2-fold to 5-fold
	  genvol = genvol * domega_p * (gen.e.E.max-gen.e.E.min)
	endif

	if (doing_heavy) then		!6-fold
	  genvol = genvol * (gen.p.E.max-gen.p.E.min)	
	endif

	normfac = normfac * genvol
	if (doing_phsp) normfac = 1.0
	wtcontribute = wtcontribute*normfac

! Close diagnostic ntuple file, if used

	if (Nntu.gt.0) call NtupleClose(filename)

! Produce output files

900	if (ngen.eq.0) goto 1000

! ... the diagnostic histograms

	i = index(base,' ')
	genfile = ' '
	write(genfile,'(a,''.gen'')') 'outfiles/'//base(1:i-1)
	genifile = ' '
	write(genifile,'(a,''.geni'')') 'outfiles/'//base(1:i-1)
	histfile = ' '
	write(histfile,'(a,''.hist'')') 'outfiles/'//base(1:i-1)
	write(6,'(1x,''... writing '',a)') genfile

	open(unit=7, file=genifile, status='unknown')
	write(7,*) hSTOP.trials
	write(7,*) hSTOP.slit_hor,hSTOP.slit_vert,hSTOP.slit_oct
	write(7,*) hSTOP.Q1_in,hSTOP.Q1_mid,hSTOP.Q1_out
	write(7,*) hSTOP.Q2_in,hSTOP.Q2_mid,hSTOP.Q2_out
	write(7,*) hSTOP.Q3_in,hSTOP.Q3_mid,hSTOP.Q3_out
	write(7,*) hSTOP.D1_in,hSTOP.D1_out
	write(7,*) hSTOP.hut
	write(7,*) hSTOP.incut,hSTOP.detec
	write(7,*) hSTOP.successes
	write(7,*)
	write(7,*) hSTOP.trials
	write(7,*) sSTOP.slit_vert,sSTOP.slit_hor,sSTOP.slit_oct
	write(7,*) sSTOP.quad_in,sSTOP.quad_mid,sSTOP.quad_out
	write(7,*) sSTOP.bm01_in,sSTOP.bm01_out,sSTOP.bm02_in
	write(7,*) sSTOP.bm02_out,sSTOP.exit,sSTOP.hut,sSTOP.detec
	write(7,*) sSTOP.quad_inm,sSTOP.quad_midm,sSTOP.quad_outm
	write(7,*) sSTOP.bm01_inm,sSTOP.bm01_outm,sSTOP.bm02_inm
	write(7,*) sSTOP.bm02_outm,sSTOP.exitm,sSTOP.hutm,sSTOP.detecm
	write(7,*) sSTOP.splatmu,sSTOP.splatpi,sSTOP.goodmu,sSTOP.goodpi

	close(7)

	open(unit=7, file=genfile, status='unknown')

	write(7,*) 'E arm Experimental Target Distributions:'
	write(7,'(6a12)') 'delta','EXPERIM','yptar','EXPERIM','xptar','EXPERIM'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.RECON.e.delta.min+(i-0.5)*H.RECON.e.delta.bin,
     >		H.RECON.e.delta.buf(i), H.RECON.e.yptar.min+(i-0.5)*
     >		H.RECON.e.yptar.bin, H.RECON.e.yptar.buf(i),
     >		H.RECON.e.xptar.min+(i-0.5)*H.RECON.e.xptar.bin,
     >		H.RECON.e.xptar.buf(i)
	enddo

	write(7,*) 'P arm Experimental Target Distributions:'
	write(7,'(6a12)') 'delta','EXPERIM','yptar','EXPERIM','xptar','EXPERIM'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.RECON.p.delta.min+(i-0.5)*H.RECON.p.delta.bin,
     >		H.RECON.p.delta.buf(i), H.RECON.p.yptar.min+(i-0.5)*
     >		H.RECON.p.yptar.bin, H.RECON.p.yptar.buf(i),
     >		H.RECON.p.xptar.min+(i-0.5)*H.RECON.p.xptar.bin,
     >		H.RECON.p.xptar.buf(i)
	enddo

	write(7,*) 'Distributions of Contributing E arm Events:'
	write(7,'(6a12)') 'delta','CONTRIB','yuptar','CONTRIB','xptar','CONTRIB'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.gen.e.delta.min+(i-0.5)*H.gen.e.delta.bin,
     >		H.gen.e.delta.buf(i), H.gen.e.yptar.min+(i-0.5)*
     >		H.gen.e.yptar.bin, H.gen.e.yptar.buf(i),
     >		H.gen.e.xptar.min+(i-0.5)*H.gen.e.xptar.bin, H.gen.e.xptar.buf(i)
	enddo

	write(7,*) 'Distributions of Contributing P arm Events:'
	write(7,'(6a12)') 'delta','CONTRIB','yptar','CONTRIB','xptar','CONTRIB'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.gen.p.delta.min+(i-0.5)*H.gen.p.delta.bin,
     >		H.gen.p.delta.buf(i), H.gen.p.yptar.min+(i-0.5)*
     >		H.gen.p.yptar.bin, H.gen.p.yptar.buf(i),
     >		H.gen.p.xptar.min+(i-0.5)*H.gen.p.xptar.bin, H.gen.p.xptar.buf(i)
	enddo

	write(7,*) 'Original E arm Events:'
	write(7,'(6a12)') 'delta','ORIGIN','yptar','ORIGIN','xptar','ORIGIN'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.gen.e.delta.min+(i-0.5)*H.gen.e.delta.bin,
     >		H.gen.e.delta.buf(i), H.gen.e.yptar.min+(i-0.5)*
     >		H.gen.e.yptar.bin, H.geni.e.yptar.buf(i),
     >		H.gen.e.xptar.min+(i-0.5)*H.gen.e.xptar.bin, H.geni.e.xptar.buf(i)
	enddo

	write(7,*) 'Original P arm Events:'
	write(7,'(6a12)') 'delta','ORIGIN','yptar','ORIGIN','xptar','ORIGIN'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H.geni.p.delta.min+(i-0.5)*H.geni.p.delta.bin,
     >		H.geni.p.delta.buf(i), H.geni.p.yptar.min+(i-0.5)*
     >		H.geni.p.yptar.bin, H.geni.p.yptar.buf(i),
     >		H.geni.p.xptar.min+(i-0.5)*H.geni.p.xptar.bin, H.geni.p.xptar.buf(i)
	enddo

	close(7)

! Counters, etc report to the screen and to the diagnostic histogram file
!	call report(6,timestring1,timestring2,central,contrib)
	open(unit=7,file=histfile,status='unknown')
	call report(7,timestring1,timestring2,central,contrib)
	close(unit=7)

1000	end

!-------------------------------------------------------------------

	subroutine inc(hist,val,weight)

	implicit none
	include 'histograms.inc'

	integer*4		ibin
	real*8			val, weight
	record /hist_entry/	hist

	ibin= nint(0.5+(val-hist.min)/hist.bin)
	if(ibin.ge.1.and.ibin.le.nHbins)then
	hist.buf(ibin) = hist.buf(ibin) + weight
	endif

	return
	end

!--------------------------------------------------------------------

	subroutine report(iun,timestring1,timestring2,central,contrib)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'
	include 'brem.inc'

	integer*4	iun
	character*23	timestring1, timestring2
	record /contrib/    contrib
	record /event_central/ central

! Output of diagnostics

! Run time
	write(iun,'(/1x,''BEGIN Time: '',a23)') timestring1
	write(iun,'(1x,''END Time:   '',a23)') timestring2

! Kinematics
	write(iun,*) 'KINEMATICS:'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'Ebeam', Ebeam, 'MeV'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)')
     >		'(dE/E)beam', dEbeam/Ebeam, '(full wid)'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'x-width', gen.xwid, 'cm'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'y-width', gen.ywid, 'cm'
	write(iun,'(9x,a12,'' = '',i15,2x,a16)')
     >		'fr_pattern',targ.fr_pattern, '1=square,2=circ'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'fr1', targ.fr1, 'cm'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'fr2', targ.fr2, 'cm'

	write(iun,*) ' '
	write(iun,'(9x,18x,2(a15,2x))') '____E arm____','____P arm____'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'angle',
     >		spec.e.theta*degrad, spec.p.theta*degrad, 'deg'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'momentum',
     >		spec.e.p, spec.p.p, 'MeV/c'

	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'x offset',
     >		spec.e.offset.x, spec.p.offset.x, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'y offset',
     >		spec.e.offset.y, spec.p.offset.y, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'z offset',
     >		spec.e.offset.z, spec.p.offset.z, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'xptar offset',
     >		spec.e.offset.xptar, spec.p.offset.xptar, 'mr'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'yptar offset',
     >		spec.e.offset.yptar, spec.p.offset.yptar, 'mr'


	write(iun,'(6x,''Central Values: '',a5,'' = '',f15.4,2x,a9)')
     >		'Q2', central.Q2/1.d6, '(GeV/c)^2'
	write(iun,'(22x,a5,'' = '',f15.4,2x,a9)') 'q', central.q, 'MeV/c'
	write(iun,'(22x,a5,'' = '',f15.4,2x,a9)') 'omega',
     >		central.omega, 'MeV'
	write(iun,'(22x,a5,'' = '',f15.4,2x,a9)') 'Em', central.Em, 'MeV'
	write(iun,'(22x,a5,'' = '',f15.4,2x,a9)') 'Pm', central.Pm, 'MeV/c'

! Target
	write(iun,*) 'TARGET specs:'
9911	format(2x,2(5x,a10,' = ',e12.6,1x,a5))
	write(iun,9911) 'A', targ.A, ' ', 'Z', targ.Z, ' '
	write(iun,9911) 'mass', targ.mass_amu, 'amu', 'mass', targ.M,'MeV'
	write(iun,9911) 'Mrec', targ.mrec_amu, 'amu', 'Mrec', targ.Mrec,'MeV'
	write(iun,9911) 'Mtar_pion',targ.Mtar_pion,'MeV',
     >		'Mrec_pion',targ.Mrec_pion,'MeV'
	write(iun,9911) 'Mtar_kaon',targ.Mtar_kaon,'MeV',
     >		'Mrec_kaon',targ.Mrec_kaon,'MeV'
	write(iun,9911) 'rho', targ.rho, 'g/cm3', 'thick', targ.thick,'g/cm2'
	write(iun,9911) 'angle', targ.angle*degrad, 'deg', 'abundancy',
     >		targ.abundancy, '%'
	write(iun,9911) 'X0', targ.X0, 'g/cm2', 'X0_cm', targ.X0_cm,'cm'
	write(iun,9911) 'length',targ.length,'cm','zoffset',targ.zoffset,'cm'
	write(iun,9911) 'xoffset',targ.xoffset,'cm','yoffset',targ.yoffset,'cm'
	write(iun,'(3x,a16,1x,a60)') 'tgt_field_file =',tgt_field_file 
	write(iun,'(t12,3a15)') '__ave__', '__lo__', '__hi__'
9912	format(1x,a15,3f15.5,2x,a6)
	write(iun,9912) 'Coulomb', targ.Coulomb.ave, targ.Coulomb.min,
     >		targ.Coulomb.max, 'MeV'
	write(iun,9912) 'Eloss_beam', targ.Eloss(1).ave,
     >		targ.Eloss(1).min, targ.Eloss(1).max, 'MeV'
	write(iun,9912) 'Eloss_e', targ.Eloss(2).ave, targ.Eloss(2).min,
     >		targ.Eloss(2).max, 'MeV'
	write(iun,9912) 'Eloss_p', targ.Eloss(3).ave, targ.Eloss(3).min,
     >		targ.Eloss(3).max, 'MeV'
	write(iun,9912) 'teff_beam', targ.teff(1).ave, targ.teff(1).min,
     >		targ.teff(1).max, 'radlen'
	write(iun,9912) 'teff_e', targ.teff(2).ave, targ.teff(2).min,
     >		targ.teff(2).max, 'radlen'
	write(iun,9912) 'teff_p', targ.teff(3).ave, targ.teff(3).min,
     >		targ.teff(3).max, 'radlen'
9913	format(1x,a15,t25,f15.5,2x,a6)
	write(iun,9913) 'musc_nsig_max', targ.musc_nsig_max, ' '
	write(iun,9913) 'musc_max_beam', targ.musc_max(1)*1000., 'mr'
	write(iun,9913) 'musc_max_e', targ.musc_max(2)*1000., 'mr'
	write(iun,9913) 'musc_max_p', targ.musc_max(3)*1000., 'mr'
! Flags
	write(iun,*) 'FLAGS:'
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_eep', doing_eep,
     >		'doing_kaon', doing_kaon, 'doing_pion', doing_pion
	write(iun,'(5x,2(2x,a19,''='',l2))') 'doing_phsp', doing_phsp,
     >		'doing_piplus', doing_piplus
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hyd_elast', doing_hyd_elast,
     >		'doing_deuterium', doing_deuterium, 'doing_heavy', doing_heavy
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydpi', doing_hydpi,
     >          'doing_deutpi', doing_deutpi, 'doing_hepi', doing_hepi
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydkaon', doing_hydkaon,
     >		'doing_deutkaon', doing_deutkaon, 'doing_hekaon', doing_hekaon
	write(iun,'(5x,3(2x,a19,''='',l2))') 'mc_smear', mc_smear,
     >		'hms_e_flag', hms_e_flag, 'doing_decay', doing_decay
	write(iun,'(5x,3(2x,a19,''='',l2)))') 'using_Eloss', using_Eloss,
     >		'using_Coulomb',using_Coulomb,'deForest_flag',deForest_flag
	write(iun,'(5x,2(2x,a19,''='',l2)))') 'correct_Eloss', correct_Eloss,
     >		'correct_raster',correct_raster
	write(iun,'(5x,2(2x,a19,''='',l2))') 
     >		'using_E_arm_montecarlo', using_E_arm_montecarlo,
     >		'using_P_arm_montecarlo', using_P_arm_montecarlo
	write(iun,'(5x,1(2x,a19,''='',l2))') 
     >		'using_tgt_field', using_tgt_field
	write(iun,'(7x,a11,''='',f10.3,a4)') 'ctau',ctau,'cm'
	write(iun,'(7x,a11,''='',i10)') 'which_kaon',which_kaon

! Counters
	write(iun,*) 'COUNTERS:'
	write(iun,'(12x,''Ngen (request) = '',i10)') ngen
	write(iun,'(12x,''Ntried         = '',i10)') ntried
	write(iun,'(12x,''Ncontribute    = '',i10)') ncontribute
	write(iun,'(12x,''Nco_no_rad_prot= '',i10)') ncontribute_no_rad_proton
	write(iun,'(12x,''-> %no_rad_prot= '',f10.3)')
     >		(100.*ncontribute_no_rad_proton/max(float(ncontribute),0.1))
	write(iun,'(/1x,''INTEGRATED WEIGHTS (number of counts in delta/Em cuts!):'')')
	write(iun,'(1x,''              MeV: wtcontr= '',e16.8)') wtcontribute/nevent

! Radiative corrections
	write(iun,*) 'RADIATIVE CORRECTIONS:'
	if (.not.using_rad) then
	  write(iun,'(x,a14,''='',l3)') 'using_rad', using_rad
	else
	  write(iun,'(4(x,a14,''='',l3))') 'use_expon',use_expon,
     >		'include_hard',include_hard,'calc_spence',calculate_spence
	  write(iun,'(2(x,a14,''='',l3))') 'using_rad', using_rad,
     >		'use_offshell_rad', use_offshell_rad
	  write(iun,'(4(x,a14,''='',i3))') 'rad_flag',rad_flag,
     >		'extrad_flag', extrad_flag, 'one_tail', one_tail,
     >		'intcor_mode', intcor_mode
	  write(iun,'(x,a14,''='',f11.3)') 'dE_edge_test',dE_edge_test
	  write(iun,'(x,a14,''='',f11.3)') 'Egamma_max', Egamma_tot_max

9914	  format(1x,a18,' = ',4f11.3)
	  write(iun,*) 'Central Values:'
	  write(iun,9914) 'hardcorfac', central.rad.hardcorfac
	  write(iun,9914) 'etatzai', central.rad.etatzai
	  write(iun,9914) 'frac(1:3)', central.rad.frac
	  write(iun,9914) 'lambda(1:3)', central.rad.lambda
	  write(iun,9914) 'bt(1:2)', central.rad.bt
	  write(iun,9914) 'c_int(0:3)', central.rad.c_int
	  write(iun,9914) 'c_ext(0:3)', central.rad.c_ext
	  write(iun,9914) 'c(0:3)', central.rad.c
	  write(iun,9914) 'g_int', central.rad.g_int
	  write(iun,9914) 'g_ext', central.rad.g_ext
	  write(iun,9914) 'g(0:3)', central.rad.g
	endif

! Miscellaneous
	write(iun,*) 'MISCELLANEOUS:'
9915	format(12x,a14,' = ',e16.6,1x,a6)
	write(iun,*) 'Note that central.sigcc is for central delta,theta,phi in both spectrometers'
	write(iun,*) ' and may give non-physical kinematics, esp. for Hydrogen'
	write(iun,*) 'Note also that AVE.sigcc is really AVER.weight (the two arenot exactly equal)'
	write(iun,9915) 'CENTRAL.sigcc', central.sigcc, ' '
	write(iun,9915) 'AVERAGE.sigcc', decdist(3)/nevent, ' '
	write(iun,9915) 'charge', EXPER.charge, 'mC'
	write(iun,9915) 'targetfac', targetfac, ' '
	write(iun,9915) 'luminosity', luminosity, 'ub^-1'
	write(iun,9915) 'luminosity', luminosity*(hbarc/100000.)**2, 'GeV^2'
	write(iun,9915) 'genvol', genvol, ' '
	write(iun,9915) 'normfac', normfac, ' '
9916	format(12x,a14,' = ',f6.1,' to ',f6.1,' in ',i4,' bins of ',f7.2,1x,a5)
	if (doing_eep .and. .not.doing_hyd_elast .and. .not.doing_deuterium) then
	  write(iun,'(12x,''Theory file:  '',a)')
     >		theory_file(1:index(theory_file,' ')-1)
	endif

! Input spectrometer limits.

	write(iun,*) 'Input Spectrometer Limits:'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.delta.min/max',
     >		SPedge.e.delta.min,SPedge.e.delta.max,'%'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.yptar.min/max',
     >		SPedge.e.yptar.min,SPedge.e.yptar.max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.xptar.min/max',
     >		SPedge.e.xptar.min,SPedge.e.xptar.max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.delta.min/max',
     >		SPedge.p.delta.min,SPedge.p.delta.max,'%'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.yptar.min/max',
     >		SPedge.p.yptar.min,SPedge.p.yptar.max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.xptar.min/max',
     >		SPedge.p.xptar.min,SPedge.p.xptar.max,'rad'


! Edges used in generation and checking, as well as range of events found
! to contribute within those edges

! ... on GENERATED qties
	write(iun,*) 'Limiting GENERATED values CONTRIBUTING to the (Em,Pm) distributions:'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)') '______used______',
     >		'_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
9917	format(1x,a18,t25,2f10.3,t50,2f10.3,2x,a5)
	write(iun,9917) 'E arm  delta', gen.e.delta.min, gen.e.delta.max,
     >		contrib.gen.e.delta.lo, contrib.gen.e.delta.hi, '%'
	write(iun,9917) 'E arm  yptar', gen.e.yptar.min*1000.,
     >		gen.e.yptar.max*1000.,
     >		contrib.gen.e.yptar.lo*1000., contrib.gen.e.yptar.hi*1000.,'mr'
	write(iun,9917) 'E arm  xptar', gen.e.xptar.min*1000.,
     >		gen.e.xptar.max*1000.,
     >		contrib.gen.e.xptar.lo*1000., contrib.gen.e.xptar.hi*1000.,'mr'
	write(iun,9917) 'P arm  delta', gen.p.delta.min, gen.p.delta.max,
     >		contrib.gen.p.delta.lo, contrib.gen.p.delta.hi, '%'
	write(iun,9917) 'P arm  yptar', gen.p.yptar.min*1000.,
     >		gen.p.yptar.max*1000.,
     >		contrib.gen.p.yptar.lo*1000., contrib.gen.p.yptar.hi*1000.,'mr'
	write(iun,9917) 'P arm  xptar', gen.p.xptar.min*1000.,
     >		gen.p.xptar.max*1000.,
     >		contrib.gen.p.xptar.lo*1000., contrib.gen.p.xptar.hi*1000.,'mr'
	write(iun,9917) 'Trec', gen.Trec.min , gen.Trec.max,
     >		contrib.gen.Trec.lo, contrib.gen.Trec.hi, 'MeV'
	if (.not.doing_hyd_elast) then
	  write(iun,9917) 'sumEgen', gen.sumEgen.min, gen.sumEgen.max,
     >		contrib.gen.sumEgen.lo, contrib.gen.sumEgen.hi, 'MeV'
	  if ((doing_deuterium .or. doing_pion .or. doing_kaon) .and. using_rad)
     >	    write(iun,*) '      *** NOTE: sumEgen.min only used in GENERATE_RAD'
	endif

! ... on TRUE qties
	write(iun,*) 'Limiting TRUE values CONTRIBUTING to the (Em,Pm) distributions:'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'E arm   E', edge.e.E.min, edge.e.E.max,
     >		contrib.tru.e.E.lo, contrib.tru.e.E.hi, 'MeV'
	write(iun,9917) 'E arm  yptar', edge.e.yptar.min*1000.,
     >		edge.e.yptar.max*1000.,
     >		contrib.tru.e.yptar.lo*1000., contrib.tru.e.yptar.hi*1000.,'mr'
	write(iun,9917) 'E arm  xptar', edge.e.xptar.min*1000., edge.e.xptar.max*1000.,
     >		contrib.tru.e.xptar.lo*1000., contrib.tru.e.xptar.hi*1000., 'mr'
	write(iun,9917) 'P arm      E', edge.p.E.min, edge.p.E.max,
     >		contrib.tru.p.E.lo, contrib.tru.p.E.hi, 'MeV'
	write(iun,9917) 'P arm  yptar', edge.p.yptar.min*1000.,
     >		edge.p.yptar.max*1000.,
     >		contrib.tru.p.yptar.lo*1000., contrib.tru.p.yptar.hi*1000.,'mr'
	write(iun,9917) 'P arm  xptar', edge.p.xptar.min*1000., edge.p.xptar.max*1000.,
     >		contrib.tru.p.xptar.lo*1000., contrib.tru.p.xptar.hi*1000., 'mr'
	write(iun,9917) 'Em', edge.Em.min, edge.Em.max, contrib.tru.Em.lo,
     >		contrib.tru.Em.hi, 'MeV'
	write(iun,9917) 'Trec', edge.Trec.min, edge.Trec.max,
     >		contrib.tru.Trec.lo, contrib.tru.Trec.hi, 'MeV'

! ... on SPECTROMETER qties
	write(iun,*) 'Limiting SPECTROMETER values CONTRIBUTING to the (Em,Pm) distributions:'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'E arm delta', SPedge.e.delta.min,
     >		SPedge.e.delta.max, contrib.SP.e.delta.lo,
     >		contrib.SP.e.delta.hi, '%'
	write(iun,9917) 'E arm  yptar', SPedge.e.yptar.min*1000.,
     >		SPedge.e.yptar.max*1000.,
     >		contrib.SP.e.yptar.lo*1000., contrib.SP.e.yptar.hi*1000., 'mr'
	write(iun,9917) 'E arm  xptar', SPedge.e.xptar.min*1000.,
     >		SPedge.e.xptar.max*1000.,
     >		contrib.SP.e.xptar.lo*1000., contrib.SP.e.xptar.hi*1000., 'mr'
	write(iun,9917) 'P arm  delta', SPedge.p.delta.min,
     >		SPedge.p.delta.max, contrib.SP.p.delta.lo,
     >		contrib.SP.p.delta.hi, '%'
	write(iun,9917) 'P arm  yptar', SPedge.p.yptar.min*1000.,
     >		SPedge.p.yptar.max*1000.,
     >		contrib.SP.p.yptar.lo*1000., contrib.SP.p.yptar.hi*1000., 'mr'
	write(iun,9917) 'P arm  xptar', SPedge.p.xptar.min*1000.,
     >		SPedge.p.xptar.max*1000.,
     >		contrib.SP.p.xptar.lo*1000., contrib.SP.p.xptar.hi*1000., 'mr'

! ... on VERTEX qties
	write(iun,*) 'Limiting VERTEX values CONTRIBUTING to the (Em,Pm) distributions:'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'Trec', VERTEXedge.Trec.min, VERTEXedge.Trec.max,
     >		contrib.vertex.Trec.lo, contrib.vertex.Trec.hi, 'MeV'
	write(iun,9917) 'Em', VERTEXedge.Em.min, VERTEXedge.Em.max,
     >		contrib.vertex.Em.lo, contrib.vertex.Em.hi, 'MeV'
	write(iun,9917) 'Pm', VERTEXedge.Pm.min, VERTEXedge.Pm.max,
     >		contrib.vertex.Pm.lo, contrib.vertex.Pm.hi, 'MeV/c'

! ... on RADIATION qties
	if (using_rad) then
	  write(iun,*) 'Limiting RADIATION values CONTRIBUTING to the (Em,Pm) distributions:'
	  write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	  write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	  write(iun,9917) 'Egamma(1)', 0., Egamma1_max,
     >		contrib.rad.Egamma(1).lo, contrib.rad.Egamma(1).hi, 'MeV'
	  write(iun,9917) 'Egamma(2)', 0., Egamma2_max, contrib.rad.Egamma(2).lo,
     >		contrib.rad.Egamma(2).hi, 'MeV'
	  write(iun,9917) 'Egamma(3)', 0., Egamma3_max, contrib.rad.Egamma(3).lo,
     >		contrib.rad.Egamma(3).hi, 'MeV'
	  write(iun,9917) 'Egamma_total', 0., Egamma_tot_max,
     >		contrib.rad.Egamma_total.lo, contrib.rad.Egamma_total.hi,'MeV'
	endif

! ... on slops
	write(iun,*) 'ACTUAL and LIMITING SLOP values used/obtained:'
	write(iun,'(t25,3a10)') '__used__', '__min__', '__max__'
9918	format(1x,a10,a12,t25,3f10.3,2x,a5)
	write(iun,9918) 'slop.MC  ', 'E arm delta', slop.MC.e.delta.used,
     >		slop.MC.e.delta.lo, slop.MC.e.delta.hi, '%'
	write(iun,9918) ' ', 'E arm yptar', slop.MC.e.yptar.used*1000.,
     >		slop.MC.e.yptar.lo*1000., slop.MC.e.yptar.hi*1000., 'mr'
	write(iun,9918) ' ', 'E arm xptar', slop.MC.e.xptar.used*1000.,
     >		slop.MC.e.xptar.lo*1000., slop.MC.e.xptar.hi*1000., 'mr'
	write(iun,9918) ' ', 'P arm delta', slop.MC.p.delta.used,
     >		slop.MC.p.delta.lo, slop.MC.p.delta.hi, '%'
	write(iun,9918) ' ', 'P arm yptar', slop.MC.p.yptar.used*1000.,
     >		slop.MC.p.yptar.lo*1000., slop.MC.p.yptar.hi*1000., 'mr'
	write(iun,9918) ' ', 'P arm xptar', slop.MC.p.xptar.used*1000.,
     >		slop.MC.p.xptar.lo*1000., slop.MC.p.xptar.hi*1000., 'mr'
	write(iun,9918) 'slop.total', 'Em', slop.total.Em.used,
     >		slop.total.Em.lo, slop.total.Em.hi, 'MeV'
	write(iun,9918) ' ', 'Pm', 0., slop.total.Pm.lo,
     >		slop.total.Pm.hi, 'MeV/c'

	write(iun,'(/)')
	return
	end

!-----------------------------------------------------------------------

	subroutine calculate_central(central)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	integer i
	logical success
	record /event_main/	main
	record /event/		vertex0, recon0
	record /event_central/	central

! Pop in generation values for central event

	if (debug(2)) write(6,*)'calc_cent: entering...'
	main.target.x = 0.0
	main.target.y = 0.0
	main.target.z = 0.0
	vertex0.Ein = Ebeam_vertex_ave
	main.target.Coulomb = targ.Coulomb.ave
	main.target.Eloss(1) = targ.Eloss(1).ave
	main.target.teff(1) = targ.teff(1).ave
	vertex0.e.delta = 0.0
	vertex0.e.yptar = 0.0
	vertex0.e.xptar = 0.0
	vertex0.e.theta = spec.e.theta
	vertex0.e.phi = spec.e.phi
	vertex0.e.P = spec.e.P*(1.+vertex0.e.delta/100.)
	vertex0.e.E = vertex0.e.P
	vertex0.p.delta = 0.0
	vertex0.p.yptar = 0.0
	vertex0.p.xptar = 0.0
	vertex0.p.theta = spec.p.theta
	vertex0.p.phi = spec.p.phi
	vertex0.p.P = spec.p.P*(1.+vertex0.p.delta/100.)
	vertex0.p.E = sqrt(vertex0.p.P**2 + Mh2)

! Complete_recon_ev for vertex0.   Note: complete_recon_ev doesn't
! call radc_init_ev or calculate teff(2,3).

! JRA: Do we want complete_ev or complete_recon_ev?  Do we want to calculate
! and/or dump other central values (for pion/kaon production, for example).

	call complete_recon_ev(vertex0,success)
	if (debug(2)) write(6,*)'calc_cent: done with comp_ev'
	if (.not.success) stop 'COMPLETE_EV failed trying to complete a CENTRAL event!'
	central.Q2 = vertex0.Q2
	central.q = vertex0.q
	central.omega = vertex0.omega
	central.Em = vertex0.Em
	central.Pm = vertex0.Pm
	main.target.teff(2) = targ.teff(2).ave
	main.target.teff(3) = targ.teff(3).ave
	if (debug(2)) write(6,*)'calc_cent: calling radc_init_ev'
	if (debug(4)) write(6,*)'calc_cent: Ein =',vertex0.Ein
	call radc_init_ev(main,vertex0)
	if (debug(2)) write(6,*)'calc_cent: done with radc_init_ev'
	if (using_rad) then
	  central.rad.hardcorfac = hardcorfac
	  central.rad.etatzai= etatzai
	  central.rad.g_int = g_int
	  central.rad.g_ext = g_ext
	  do i = 0, 3
	    if (i.gt.0) central.rad.frac(i) = frac(i)
	    if (i.gt.0) central.rad.lambda(i) = lambda(i)
	    if (i.gt.0.and.i.lt.3) central.rad.bt(i) = bt(i)
	    central.rad.c_int(i) = c_int(i)
	    central.rad.c_ext(i) = c_ext(i)
	    central.rad.c(i) = c(i)
	    central.rad.g(i) = g(i)
	  enddo
	endif
	if (debug(4)) write(6,*)'calc_cent: at 1'

	do i = 1, neventfields
	  recon0.all(i) = vertex0.all(i)
	enddo
	call complete_main(.true.,main,vertex0,recon0,success)
	central.sigcc = main.sigcc

	if (debug(2)) write(6,*)'calc_cent: ending...'
	return
	end

!-------------------------------------------------------------------------

	subroutine montecarlo(in_E_arm,in_P_arm,main,out_E_arm,out_P_arm,success)

	implicit none
	include 'simulate.inc'

* Track-coordinate and spectrometer common blocks

	real*8 x_E_arm,y_E_arm,z_E_arm,dx_E_arm,dy_E_arm,delta_E_arm
	real*8 x_P_arm,y_P_arm,z_P_arm,dx_P_arm,dy_P_arm,delta_P_arm
	real*8 xfp, yfp, dxfp, dyfp
	real*8 eloss_E_arm, eloss_P_arm, r, beta, dangles(2)
	logical success
	logical ok_E_arm, ok_P_arm
	record /arm_full/	in_E_arm, in_P_arm, out_E_arm, out_P_arm
	record /event_main/	main

	real*8 beta_electron
	parameter (beta_electron = 1.)
	real*8 tmpfact
	real*8 fry	!fast raster y position.
	real*8 frx	!fast raster x position - left as look downstream -GAW
	real*8 m2	!mass for call to mc_hms(sos). Changes if decay

	real*8 zero
	parameter (zero=0.0d0)	!double precision zero for subroutine calls

! Prepare the event for the Monte Carlo's and/or spectrometer cuts

	success = .false.
	decdist(6) = 0.0			!resfac (see simulate.inc)
	if (correct_raster) then
	  fry = -main.target.rastery
	  frx =  main.target.rasterx  
	else
	  fry = 0.0
          frx = 0.0
	endif

	if (debug(3)) write(6,*) 'mc: using p,e arm mc =',
     >		using_P_arm_montecarlo,using_E_arm_montecarlo

!____ P arm ________________

! Go from TRUE to SPECTROMETER quantities by computing target distortions
! ... ionization loss correction (if requested)

	if (using_Eloss) then
	  if (debug(3)) write(6,*)'mc: p arm stuff0 =',
     >		in_P_arm.E,main.target.Eloss(3),spec.p.P
	  main.SP.p.delta = (sqrt(abs((in_P_arm.E-main.target.Eloss(3))**2
     >		          -Mh2))-spec.p.P) / spec.p.P*100.
	else
	  if (debug(3)) write(6,*)'mc: p arm stuff1 =',in_P_arm.delta
	  main.SP.p.delta = in_P_arm.delta
	endif

! ... multiple scattering

	if (mc_smear) then
	  beta = in_P_arm.p/in_P_arm.E
	  call target_musc(in_P_arm.p, beta, main.target.teff(3), dangles)
	else
	  dangles(1)=0.0
	  dangles(2)=0.0
	endif
	main.SP.p.yptar = in_P_arm.yptar + dangles(1)
	main.SP.p.xptar = in_P_arm.xptar + dangles(2)

! CASE 1: Using the spectrometer Monte Carlo

	if (using_P_arm_montecarlo) then

! ... change to P arm spectrometer coordinates (TRANSPORT system),
! ... remembering to add in offsets from the dbase

	  delta_P_arm = main.SP.p.delta
	  x_P_arm = -(main.target.y+spec.p.offset.y)
	  if(hms_e_flag)then
	    y_P_arm = -(main.target.z+spec.p.offset.z)*spec.p.sin_th+
     >		(main.target.x+spec.p.offset.x)*spec.p.cos_th
	    z_P_arm =  (main.target.x+spec.p.offset.x)*spec.p.sin_th+
     >		(main.target.z+spec.p.offset.z)*spec.p.cos_th
	  else
	    y_P_arm = (main.target.z+spec.p.offset.z)*spec.p.sin_th+
     >		(main.target.x+spec.p.offset.x)*spec.p.cos_th
	    z_P_arm = -(main.target.x+spec.p.offset.x)*spec.p.sin_th+
     >		(main.target.z+spec.p.offset.z)*spec.p.cos_th
	  endif
	  dx_P_arm = main.SP.p.xptar - spec.e.offset.xptar
	  dy_P_arm = main.SP.p.yptar - spec.e.offset.xptar
! GAW - project to z=0 to to compare with reconstructed target positions
	  in_P_arm.z = y_P_arm - z_P_arm*dy_P_arm

! ........ drift this position back to z=0, the plane through the target center

	  x_P_arm = x_P_arm - z_P_arm*dx_P_arm
	  y_P_arm = y_P_arm - z_P_arm*dy_P_arm
	  z_P_arm = 0.0

! ... Monte Carlo through the P arm! if we succeed, continue ...
! ... Here's what's passed:
!	spectrometer central momentum
!	spectrometer central theta angle
!	momentum delta
!	x, y, z, theta, and phi in TRANSPORT coordinates
!	x, y, theta, and phi at the focal plane
!	radiation length of target (NOT USED!)
!	particle mass (modified if the hadron decays)
!	multiple scattering flag
!	wcs smearing flag
!	decay flag
!	resmult=resolution factor for DCs (see simulate.inc)
!	vertical position at target (for reconstruction w/raster MEs).
!	ok flag

	  m2 = Mh2
	  if (hms_e_flag) then
	    call mc_sos(spec.p.P, spec.p.theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		zero, m2, mc_smear, mc_smear, doing_decay,
     >		decdist(6), fry, ok_P_arm)
	  else
	    call mc_hms(spec.p.P, spec.p.theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		zero, m2, mc_smear, mc_smear, doing_decay,
     >		decdist(6), frx, fry, ok_P_arm)
	  endif

	  if (.not.ok_P_arm) then
	     if (debug(3)) write(6,*)'mc: ok_P_arm =',ok_P_arm
	     return
	  endif

	  main.RECON.p.delta = delta_P_arm
	  main.RECON.p.yptar = dy_P_arm
	  main.RECON.p.xptar = dx_P_arm
	  main.RECON.p.z = y_P_arm
	  main.FP.p.x = xfp
	  main.FP.p.dx = dxfp
	  main.FP.p.y = yfp
	  main.FP.p.dy = dyfp

! CASE 2: Not using the detailed Monte Carlo, just copy the IN event to the
! OUT record

	else
	  if (using_tgt_field) then ! GAW for E93-026 tgt field
  	    call transport_proton(main,in_P_arm.z)
	  else
	    main.RECON.p.delta = main.SP.p.delta
	    main.RECON.p.yptar = main.SP.p.yptar
	    main.RECON.p.xptar = main.SP.p.xptar
	  endif
	endif

! Back to TRUE qties

	out_P_arm.delta = main.RECON.p.delta
	out_P_arm.yptar = main.RECON.p.yptar
	out_P_arm.xptar = main.RECON.p.xptar
	out_P_arm.z = main.RECON.p.z
	out_P_arm.P = spec.p.P*(1.+out_P_arm.delta/100.)
	out_P_arm.E = sqrt(out_P_arm.P**2 + Mh2)
	call physics_angles(spec.p.theta,spec.p.phi,out_P_arm.xptar,
     >		out_P_arm.yptar,out_P_arm.theta,out_P_arm.phi)

! ... remove energy loss.  Add back most probable energy loss (typeflag=4)

	if (correct_Eloss) then
	  call trip_thru_target (3, zero, out_P_arm.E,
     >		out_P_arm.theta, eloss_P_arm, r,Mh,4)
	  out_P_arm.E = out_P_arm.E + eloss_P_arm
	  out_P_arm.E = max(out_p_arm.E,sqrt(Mh2+0.000001)) !can get P~0 when calculating hadron momentum-->P<0 after eloss
	  out_P_arm.P = sqrt(out_P_arm.E**2-Mh2)
	  out_P_arm.delta = (out_P_arm.P-spec.p.P)/spec.p.P*100.
	endif

!___ E arm ______________

! Go from TRUE to SPECTROMETER quantities by computing target distortions

! ... ionization loss correction (if requested) and Coulomb deceleration

	if (debug(3)) write(6,*)'mc: e arm stuff0 =',
     >		in_E_arm.E,main.target.Eloss(2),spec.e.P
	main.SP.e.delta=100.* (in_E_arm.E - main.target.Eloss(2)
     >		- main.target.Coulomb - spec.e.P) / spec.e.P

! ... multiple scattering

	if (mc_smear) then
	  call target_musc(in_E_arm.p, beta_electron, main.target.teff(2), dangles)
	else
	  dangles(1)=0.0
	  dangles(2)=0.0
	endif
	main.SP.e.yptar = in_E_arm.yptar + dangles(1)
	main.SP.e.xptar = in_E_arm.xptar + dangles(2)

! CASE 1: Using the spectrometer Monte Carlo

	if (using_E_arm_montecarlo) then

! ... change to E arm spectrometer coordinates (TRANSPORT system),
! ... remembering to add in offsets from the dbase

	  delta_E_arm = main.SP.e.delta
	  x_E_arm = -(main.target.y+spec.e.offset.y)
	  if(hms_e_flag)then
	    y_E_arm = (main.target.z+spec.e.offset.z)*spec.e.sin_th+
     >		(main.target.x+spec.e.offset.x)*spec.e.cos_th
	    z_E_arm = -(main.target.x+spec.e.offset.x)*spec.e.sin_th+
     >		(main.target.z+spec.e.offset.z)*spec.e.cos_th
	  else
	    y_E_arm = -(main.target.z+spec.e.offset.z)*spec.e.sin_th+
     >		(main.target.x+spec.e.offset.x)*spec.e.cos_th
	    z_E_arm = (main.target.x+spec.e.offset.x)*spec.e.sin_th+
     >		(main.target.z+spec.e.offset.z)*spec.e.cos_th
	  endif
	  dx_E_arm = main.SP.e.xptar - spec.e.offset.xptar
	  dy_E_arm = main.SP.e.yptar - spec.e.offset.yptar
! GAW -project to z=0 to compare with reconstructed target positions
	  in_E_arm.z = y_E_arm - z_E_arm*dy_E_arm

	  if (using_tgt_field) then      ! do target field tracking - GAW

	    if (debug(6)) then
	       write(*,*) '------------------------------'
	       write(*,'("frx,fry,x_E_arm =      ",3f12.5)') frx,fry,x_E_arm
	    endif

	    call track_from_tgt(x_E_arm,y_E_arm,z_E_arm,dx_E_arm,dy_E_arm,
     >                          -spec.e.P*(1+delta_E_arm/100.),Me,-1,ok_E_arm)
	  endif 
! GAW - end 99/11/3

! ........ drift this position back to z=0, the plane through the target center

	  x_E_arm = x_E_arm - z_E_arm*dx_E_arm
	  y_E_arm = y_E_arm - z_E_arm*dy_E_arm
	  z_E_arm = 0.0

          call print_coord2('virtual track z=0: ',x_E_arm,y_E_arm,z_E_arm,dx_E_arm,dy_E_arm)  ! GAW
	  main.SP.e.z=y_E_arm

! ... Monte Carlo through the E arm! if we succeed, continue ...
! ... Here's what's passed:
!	spectrometer central momentum
!	spectrometer central theta angle
!	momentum delta
!	x, y, z, theta, and phi in TRANSPORT coordinates
!	x, y, theta, and phi at the focal plane
!	radiation length of target (NOT USED!)
!	particle mass
!	multiple scattering flag
!	wcs smearing flag
!	resmult=resolution factor for DCs (see simulate.inc)
!	decay flag (hardwired to .false. for electron arm).
!	ok flag

	  m2 = me2
	  if (hms_e_flag) then
	    call mc_hms(spec.e.P, spec.e.theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		zero, me2, mc_smear, mc_smear, .false.,
     >		tmpfact, frx, fry, ok_E_arm)
	  else
	    call mc_sos(spec.e.P, spec.e.theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		zero, me2, mc_smear, mc_smear, .false.,
     >		tmpfact, fry, ok_E_arm)
	  endif
	  decdist(6)=decdist(6)+tmpfact

	  if (.not.ok_E_arm) then
	    if (debug(3)) write(6,*)'mc: ok_E_arm =',ok_E_arm
	    return
	  endif
	  call print_coord2('recon tgt: ',x_E_arm,y_E_arm,0.,dx_E_arm,dy_E_arm) ! GAW
	  if (debug(6)) write(*,*) '===================================================='
	  main.RECON.e.delta = delta_E_arm
	  main.RECON.e.yptar = dy_E_arm
	  main.RECON.e.xptar = dx_E_arm
	  main.RECON.e.z = y_E_arm
	  main.FP.e.x = xfp
	  main.FP.e.dx = dxfp
	  main.FP.e.y = yfp
	  main.FP.e.dy = dyfp

! CASE 2: Not using the detailed Monte Carlo, just copy the IN event to the
! OUT record

	else
	  main.RECON.e.delta = main.SP.e.delta
	  main.RECON.e.yptar = main.SP.e.yptar
	  main.RECON.e.xptar = main.SP.e.xptar
	endif

! Back to TRUE qties

	out_E_arm.delta = main.RECON.e.delta
	out_E_arm.yptar = main.RECON.e.yptar
	out_E_arm.xptar = main.RECON.e.xptar
	out_E_arm.z = main.RECON.e.z
	out_E_arm.P = spec.e.P*(1.+out_E_arm.delta/100.)
	out_E_arm.E = out_E_arm.P
	call physics_angles(spec.e.theta,spec.e.phi,out_E_arm.xptar,
     >		out_E_arm.yptar,out_E_arm.theta,out_E_arm.phi)


! ... remove energy loss, and correct for Coulomb deceleration

	if (correct_Eloss) then
	  call trip_thru_target (2, zero, out_E_arm.E, out_E_arm.theta,
     &                              eloss_E_arm, r, Me, 4)
	  out_E_arm.E = out_E_arm.E + eloss_E_arm
	endif
	out_E_arm.E = out_E_arm.E + targ.Coulomb.ave
	out_E_arm.P = out_E_arm.E
	out_E_arm.delta = (out_E_arm.P-spec.e.P)/spec.e.P*100.

! Made it!
	success = .true.

	return
	end
