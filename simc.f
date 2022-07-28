	program simc

! Last modified:
!
!	This program represents variations on themes by
!		N. C. R. Makins (Argonne National Lab),
!		T. G. O'Neill (Argonne National Lab), and
!		Seemingly Countless Others (Virtually Everywhere).
!
	implicit none
!	include 'simulate_init.inc'
	include 'simulate.inc'
	include 'histograms_init.inc'
	include 'radc.inc'
	include 'hbook.inc'
	include 'sos/struct_sos.inc'
	include 'hms/struct_hms.inc'
	include 'hrsr/struct_hrsr.inc'
	include 'hrsl/struct_hrsl.inc'
	include 'shms/struct_shms.inc'
	include 'calo/struct_calo.inc'

	integer*4	i, ninform
	integer*8       itime1,itime2
	real*8		r, sum_sigcc, tmpnum
	real*8		domega_e, domega_p !populated e/hadron solid angles.
	logical		success
	logical		pass_cuts
	character	filename*80, genfile*80, histfile*80, timestring1*30
	character	timestring2*30,genifile*80
	type(event)::		vertex, vertex0, orig, recon
	type(event_main)::	main
	type(contribtype)::	contrib
	type(histograms)::	H
	type(event_central)::	central
	type(sums_twoarm)::	sumerr, sumerr2, aveerr, resol

	real*8 one
	parameter (one=1.0e0)	!double precision 1 for subroutine calls

	real*8 grnd
	real*8 ang_targ_earm,ang_targ_parm
	logical restorerndstate
c 

! INITIALIZE
!-----------

! ... initialize all the *min/max variables here since we can 
! ... no longer do it in the * include file
	call min_max_init(contrib)

	

C DJG Remove cernlib requirements July, 2020
! ... initialize histogram area for HBOOK
c	call hlimit(PawSize)

! ... read in the data base

	call dbase_read(H)

	if (debug(3)) write(6,*) 'Main after dbrd: p,e th',
     >	    spec%p%theta,spec%e%theta,using_P_arm_montecarlo,using_E_arm_montecarlo

! ... open diagnostic ntuple file, if requested

	i = index(base,' ')
	if (Nntu.gt.0) then
	  filename = 'worksim/'//base(1:i-1)//'.bin'
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

	call calculate_central(central,vertex0)
	if (debug(2)) write(6,*)'sim: done with calc_central'
	if (debug(2)) write(6,*)'central%sigcc=',central%sigcc
	if (debug(4)) write(6,*)'sim: at 1'

	targetfac=targ%mass_amu/3.75914d+6/(targ%abundancy/100.)
     >		*abs(cos(targ%angle))/(targ%thick*1000.)
	if (debug(4)) write(6,*)'sim: at 2'

! ... Note: EXPER.charge is in mC and the luminosity comes out in fm^-2
! ...   now luminosity is in ub^-1

	luminosity = EXPER%charge/targetfac
	if (debug(4)) write(6,*)'sim: at 3'

! ... initialize the random number generator and number of attempts
	
       if (random_state_file .eq. ' ') then
         call sgrnd(random_seed)	
       else
	 call  start_file_random_state(random_state_file)
       endif	
c
       write(*,*) "Initial random vector save to file: ",start_random_state_file
	call save_random_state(start_random_state_file)
	ntried = 0

! GAW - insert calls to initialize target field for both arms
! mkj 10-20-2003 add if statements to determine the angle magnitude
!                and sign between target direction and e-arm or p-arm.
!       
	if(using_tgt_field) then
	   if ( degrad*abs(targ_Bphi-spec%e%phi) .lt. .01) then
	      if (targ_Bangle .ge. spec%e%theta) then
		 ang_targ_earm = -1*sin(spec%e%phi)*(targ_Bangle-spec%e%theta)
	      else
		 ang_targ_earm = +1*sin(spec%e%phi)*(spec%e%theta-targ_Bangle)
	      endif
	   elseif ( degrad*abs(targ_Bphi-spec%e%phi)-180.0 .lt. .01) then 
	      ang_targ_earm = +1*sin(spec%e%phi)*(targ_Bangle+spec%e%theta)
	   else
	      write(*,*) ' Error determining angle between target and e-arm'
	      write(*,*) ' targ_Bangle =',targ_Bangle*degrad,' targ_Bphi =',targ_Bphi*degrad
	      write(*,*) ' Central e-arm angle =',spec%e%theta*degrad,' e-arm phi =',spec%e%phi*degrad
	   endif
c
	   if ( degrad*abs(targ_Bphi-spec%p%phi) .lt. .01) then
	      if (targ_Bangle .ge. spec%p%theta) then
		 ang_targ_parm = -1*sin(spec%p%phi)*(targ_Bangle-spec%p%theta)
	      else
		 ang_targ_parm = +1*sin(spec%p%phi)*(targ_Bangle-spec%p%theta)
	      endif
	   elseif ( degrad*abs(targ_Bphi-spec%p%phi)-180.0 .lt. .01) then 
	      ang_targ_parm = +1*sin(spec%p%phi)*(targ_Bangle+spec%p%theta)
	   else
	      write(*,*) ' Error determining angle between target and p-arm'
	      write(*,*) ' targ_Bangle =',targ_Bangle*degrad,' targ_Bphi =',targ_Bphi*degrad
	      write(*,*) ' Central p-arm angle =',spec%p%theta*degrad,' p-arm phi =',spec%p%theta*degrad
	   endif
c
	   write(*,*) ' targ_Bangle =',targ_Bangle*degrad,' targ_Bphi =',targ_Bphi*degrad
	   write(*,*) ' Angle between target and e-arm',ang_targ_earm*degrad	
	   write(*,*) ' Angle between target and p-arm',ang_targ_parm*degrad	
c

	   call trgInit(tgt_field_file,ang_targ_earm*degrad,0.,
     >      	ang_targ_parm*degrad,0.)
	endif

! LOOP OVER SIMULATED EVENTS
!---------------------------
	write(6,*) ' '
	itime1=time8()
	call ctime(itime1,timestring1)
cdg	call date (timestring1)
cdg	call time (timestring1(11:23))
	nevent = 0
	ninform = 1000
	sum_sigcc = 0.0	!sum of main.sigcc (for generating ave sigcc)

	do while (nevent.lt.abs(ngen))

! reset decay distance, initial hadron mass (modified if particle decays),
! and some ntuple variables.

	  Mh2_final = Mh2
	  ntup%radphot = 0.0
	  ntup%radarm = 0.0
	  decdist = 0.0		!decay distance.
	  ntup%resfac = 0.0
	  ntup%sigcm = 0.0
	  ntup%krel = 0.0
	  ntup%mm = 0.0
	  ntup%mmA = 0.0
	  ntup%t = 0.0

! Setup for this event

! ... keep the human interested

	  ntried = ntried + 1

	  if(debug(2)) then
     	    write(6,*)'********************************************'
	    write(6,*)'SIM:  ntried =',ntried
	    write(6,*)'      ncontribute =',ncontribute
	  endif

	  if (mod(ntried,ninform).eq.1) then
	    write (6,'(1x,a,i9,a,i8,a,e11.4)') 'Generating Event',
     >		ntried, ' ... ', nevent,'  successes so far - Monitor:',
     >	  wtcontribute*luminosity/ntried
	    if (ntried.ge.5000) ninform = 20000
	  endif

! ... generate an event
	  call generate(main,vertex,orig,success)
	  if(debug(2)) write(6,*)'sim: after gen, success =',success

! Run the event through various manipulations, checking to see whether
! it made it at each stage

! ... run through spectrometers and check SP cuts

	  if(debug(3)) write(6,*)'sim: before mc: orig.p.E =',orig%p%E
	  if(success)call montecarlo(orig,main,recon,success)
	  if(debug(2)) write(6,*)'sim: after mc, success =',success

! ... calculate everything else about the event

	  if (success) then
	    call complete_recon_ev(recon,success)
	  endif

	  if(debug(2)) write(6,*)'sim: after comp_ev, success =',success
	  if(debug(5)) write(6,*) 'recon%Em,recon%Pm',recon%Em,recon%Pm
! ... calculate remaining pieces of the main structure

	  if (success) call complete_main(.false.,main,vertex,vertex0,recon,success)
! ... Apply SPedge cuts to success if hard_cuts is set.
	    pass_cuts = .not. (
     >	      recon%e%delta .le. (SPedge%e%delta%min+slop%MC%e%delta%used) .or.
     >	      recon%e%delta .ge. (SPedge%e%delta%max-slop%MC%e%delta%used) .or.
     >	      recon%e%yptar .le. (SPedge%e%yptar%min+slop%MC%e%yptar%used) .or.
     >	      recon%e%yptar .ge. (SPedge%e%yptar%max-slop%MC%e%yptar%used) .or.
     >	      recon%e%xptar .le. (SPedge%e%xptar%min+slop%MC%e%xptar%used) .or.
     >	      recon%e%xptar .ge. (SPedge%e%xptar%max-slop%MC%e%xptar%used) .or.
     >	      recon%p%delta .le. (SPedge%p%delta%min+slop%MC%p%delta%used) .or.
     >	      recon%p%delta .ge. (SPedge%e%delta%max-slop%MC%p%delta%used) .or.
     >	      recon%p%yptar .le. (SPedge%p%yptar%min+slop%MC%p%yptar%used) .or.
     >	      recon%p%yptar .ge. (SPedge%p%yptar%max-slop%MC%p%yptar%used) .or.
     >	      recon%p%xptar .le. (SPedge%p%xptar%min+slop%MC%p%xptar%used) .or.
     >	      recon%p%xptar .ge. (SPedge%p%xptar%max-slop%MC%p%xptar%used) )

	  if (hard_cuts) then		!apply spec. limit and Em cuts.
	    if (.not.pass_cuts) success = .false.
	    if (doing_eep .and. (recon%Em .gt. cuts%Em%max)) success = .false.
	  endif

	  if (success) sum_sigcc = sum_sigcc + main%sigcc

! Em and Pm histos are NEW.  Check that they don't suffer from energy
! shifts (e.g. coulomb distortions - see update_range call in limits_update).

	  if(ntried.gt.0)then
	    call inc(H%geni%e%delta, vertex%e%delta, one)
	    call inc(H%geni%e%yptar, vertex%e%yptar, one)
	    call inc(H%geni%e%xptar, -vertex%e%xptar, one)
	    call inc(H%geni%p%delta, vertex%p%delta, one)
	    call inc(H%geni%p%yptar, vertex%p%yptar, one)
	    call inc(H%geni%p%xptar, -vertex%p%xptar, one)
	    call inc(H%geni%Em, vertex%Em, one)
	    call inc(H%geni%Pm, vertex%Pm, one)
	  endif

	  if (success) then
	    if (debug(2)) write(6,*)'sim: We have a success!'

! ... Output ntuple and histograms if passed all cuts, or if not using
! ... SPedge limits as hard cuts.

! ... increment the "experimental" spectrometer arm and "contribution"
! ... histograms
	    call inc(H%RECON%e%delta, main%RECON%e%delta, main%weight)
	    call inc(H%RECON%e%yptar, main%RECON%e%yptar, main%weight)
	    call inc(H%RECON%e%xptar, main%RECON%e%xptar, main%weight)
	    call inc(H%RECON%p%delta, main%RECON%p%delta, main%weight)
	    call inc(H%RECON%p%yptar, main%RECON%p%yptar, main%weight)
	    call inc(H%RECON%p%xptar, main%RECON%p%xptar, main%weight)
	    call inc(H%RECON%Em, recon%Em, one)
	    call inc(H%RECON%Pm, recon%Pm, one)
	    call inc(H%gen%e%delta, vertex%e%delta, one)
	    call inc(H%gen%e%yptar, vertex%e%yptar, one)
	    call inc(H%gen%e%xptar, -vertex%e%xptar, one)
	    call inc(H%gen%p%delta, vertex%p%delta, one)
	    call inc(H%gen%p%yptar, vertex%p%yptar, one)
	    call inc(H%gen%p%xptar, -vertex%p%xptar, one)
	    call inc(H%gen%Em, vertex%Em, one)
	    
! ... update counters and integrated weights AND keep track of resolution
! ...  FOR EVENTS INSIDE OF spectrometer population limits (events are only
! ...  generated to fill region inside of the given delta limits) and below
! ...  cuts.Em.max (used for elastic/quasielastic only).

	    if (debug(4)) write(6,*)'sim: cut'

	    ncontribute = ncontribute + 1
	    if (debug(4)) write(6,*)'sim: ncontribute =',ncontribute
	    if (.not.rad_proton_this_ev) ncontribute_no_rad_proton =
     >			ncontribute_no_rad_proton + 1


! Using main.SP.e.delta means that we correct for eloss/coulomb, which we
! shouldn't do. Using vertex.e.delta means radiation counts as error in recon.
! --> USE MAIN.SP AND UNDERESTIMATE RESOLUTION (MISSING TARGET ELOSS)

! For xptar/yptar, main.SP.e.xptar/yptar corrects for mult. scatt. --> USE VERTEX.

! For z, we HAVE to use main.SP.e.z, but there are not corrections, so USE MAIN.SP
!  vertex.e.z is not defined (main.target.x/y/z --> y_E_arm
!  (offsets and rotation)  --> spectrometer M.C. --> RECON.e.z 

	    if (pass_cuts) then
	      npasscuts = npasscuts + 1
	      wtcontribute = wtcontribute + main%weight
	      sumerr%e%delta = sumerr%e%delta + (recon%e%delta - main%SP%e%delta)
	      sumerr%e%xptar = sumerr%e%xptar + (recon%e%xptar - vertex%e%xptar)
	      sumerr%e%yptar = sumerr%e%yptar + (recon%e%yptar - vertex%e%yptar)
	      sumerr%e%ytar  = sumerr%e%ytar  + (recon%e%z  - main%SP%e%z)
	      sumerr%p%delta = sumerr%p%delta + (recon%p%delta - main%SP%p%delta)
	      sumerr%p%xptar = sumerr%p%xptar + (recon%p%xptar - vertex%p%xptar)
	      sumerr%p%yptar = sumerr%p%yptar + (recon%p%yptar - vertex%p%yptar)
	      sumerr%p%ytar  = sumerr%p%ytar  + (recon%p%z  - main%SP%p%z)
	      sumerr2%e%delta= sumerr2%e%delta + (recon%e%delta - main%SP%e%delta)**2
	      sumerr2%e%xptar= sumerr2%e%xptar + (recon%e%xptar - vertex%e%xptar)**2
	      sumerr2%e%yptar= sumerr2%e%yptar + (recon%e%yptar - vertex%e%yptar)**2
	      sumerr2%e%ytar = sumerr2%e%ytar  + (recon%e%z  - main%SP%e%z)**2
	      sumerr2%p%delta= sumerr2%p%delta + (recon%p%delta - main%SP%p%delta)**2
	      sumerr2%p%xptar= sumerr2%p%xptar + (recon%p%xptar - vertex%p%xptar)**2
	      sumerr2%p%yptar= sumerr2%p%yptar + (recon%p%yptar - vertex%p%yptar)**2
	      sumerr2%p%ytar = sumerr2%p%ytar  + (recon%p%z  - main%SP%p%z)**2
	    endif

! ... update the "contribution" and "slop" limits
	    call limits_update(main,vertex,orig,recon,doing_deuterium,
     >		doing_pion,doing_kaon,doing_delta,doing_rho,contrib,slop)

	  endif ! <success>

! ... write out line to diagnostic ntuple file, if requested
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
	write (6,'(1x,''---> Last Event '',i9,'' ...'',i9,''  successes'')') ntried, nevent
	itime2=time8()
	call ctime(itime2,timestring2)
c	call date (timestring2)
c	call time (timestring2(11:23))

! NORMALIZE!

! ... put in the luminosity and efficiency factors

	normfac = luminosity/ntried*nevent

! ... multiply in the relevant phase spaces (see event.f for description
! ... of generated variables.  Electron angles for all cases.
! ... add hadron angles and electron energy for all but H(e,e'p).
! ... add hadron energy for A(e,e'p) (what about phase_space?)
! ... Cross sections are all in microbarn/MeV**i/sr**j (i,j depend on reaction)

	domega_e=(gen%e%yptar%max-gen%e%yptar%min)*(gen%e%xptar%max-gen%e%xptar%min)
	if(doing_rho) then
	   domega_p = 4.*pi
	else
	   domega_p=(gen%p%yptar%max-gen%p%yptar%min)*(gen%p%xptar%max-gen%p%xptar%min)
	endif

	genvol = domega_e
!	genvol_inclusive = genvol	!may want dOmega, or dE*dOmega

! ... 2-fold to 5-fold.
	if (doing_deuterium.or.doing_heavy.or.doing_pion.or.doing_kaon
     >      .or.doing_delta.or.doing_rho .or. doing_semi) then
	  genvol = genvol * domega_p * (gen%e%E%max-gen%e%E%min)
	endif

	if (doing_heavy.or.doing_semi) then		!6-fold
	  genvol = genvol * (gen%p%E%max-gen%p%E%min)	
	endif

	normfac = normfac * genvol
	if (doing_phsp) normfac = 1.0
	wtcontribute = wtcontribute*normfac

! Close diagnostic ntuple file, if used

	if (Nntu.gt.0) call NtupleClose(filename)

! Calculate resolutions

	if (npasscuts.gt.1) then
	  tmpnum = dble(npasscuts)
	  aveerr%e%delta = sumerr%e%delta/tmpnum
	  aveerr%e%xptar = sumerr%e%xptar/tmpnum
	  aveerr%e%yptar = sumerr%e%yptar/tmpnum
	  aveerr%e%ytar  = sumerr%e%ytar /tmpnum
	  aveerr%p%delta = sumerr%p%delta/tmpnum
	  aveerr%p%xptar = sumerr%p%xptar/tmpnum
	  aveerr%p%yptar = sumerr%p%yptar/tmpnum
	  aveerr%p%ytar  = sumerr%p%ytar /tmpnum
	  resol%e%delta = sqrt(max(0.,(sumerr2%e%delta/tmpnum) -
     >				     (sumerr%e%delta/tmpnum)**2 ))
	  resol%e%xptar = sqrt(max(0.,(sumerr2%e%xptar/tmpnum) -
     >				     (sumerr%e%xptar/tmpnum)**2 ))
	  resol%e%yptar = sqrt(max(0.,(sumerr2%e%yptar/tmpnum) -
     >				     (sumerr%e%yptar/tmpnum)**2 ))
	  resol%e%ytar  = sqrt(max(0.,(sumerr2%e%ytar/tmpnum) -
     >				     (sumerr%e%ytar/tmpnum)**2 ))
	  resol%p%delta = sqrt(max(0.,(sumerr2%p%delta/tmpnum) -
     >				     (sumerr%p%delta/tmpnum)**2 ))
	  resol%p%xptar = sqrt(max(0.,(sumerr2%p%xptar/tmpnum) -
     >				     (sumerr%p%xptar/tmpnum)**2 ))
	  resol%p%yptar = sqrt(max(0.,(sumerr2%p%yptar/tmpnum) -
     >				     (sumerr%p%yptar/tmpnum)**2 ))
	  resol%p%ytar  = sqrt(max(0.,(sumerr2%p%ytar/tmpnum) -
     >				     (sumerr%p%ytar/tmpnum)**2 ))
	endif

!	write(6,9910) 'e- delta=',10.*aveerr.e.delta,10.*resol.e.delta,'x10^-3'
!	write(6,9910) 'e- xptar=',1000.*aveerr.e.xptar,1000.*resol.e.xptar,'mr'
!	write(6,9910) 'e- yptar=',1000.*aveerr.e.yptar,1000.*resol.e.yptar,'mr'
!	write(6,9910) 'e- ytar =',10.*aveerr.e.ytar,10.*resol.e.ytar ,'mm'
!	write(6,9910) 'p  delta=',10.*aveerr.p.delta,10.*resol.p.delta,'x10-3'
!	write(6,9910) 'p  xptar=',1000.*aveerr.p.xptar,1000.*resol.p.xptar,'mr'
!	write(6,9910) 'p  yptar=',1000.*aveerr.p.yptar,1000.*resol.p.yptar,'mr'
!	write(6,9910) 'p  ytar =',10.*aveerr.p.ytar,10.*resol.p.ytar,'mm'
!9910	format(1x,a10,2f10.3,3x,a)

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

	write(6,*) 'Come back soon!'

	open(unit=7, file=genifile, status='unknown')
	if (electron_arm.eq.1 .or. hadron_arm.eq.1) then
	  write(7,*) 'HMS Trials:           ',hSTOP_trials
	  write(7,*) 'Slit hor/vert/corners ',hSTOP_slit_hor,hSTOP_slit_vert,hSTOP_slit_oct
	  write(7,*) 'Q1 entrance/mid/exit  ',hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out
	  write(7,*) 'Q2 entrance/mid/exit  ',hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out
	  write(7,*) 'Q3 entrance/mid/exit  ',hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out
	  write(7,*) 'Dipole entrance/exit  ',hSTOP_D1_in,hSTOP_D1_out
	  write(7,*) 'Events reaching hut   ',hSTOP_hut
	  write(7,*) 'DC1, DC2, Scin, Cal   ',hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal
	  write(7,*) 'Successes             ',hSTOP_successes
	  write(7,*)
	endif
	if (electron_arm.eq.2 .or. hadron_arm.eq.2) then
	  write(7,*) 'SOS Trials:           ',sSTOP_trials
	  write(7,*) 'Slit hor/vert/corners ',sSTOP_slit_vert,sSTOP_slit_hor,sSTOP_slit_oct
	  write(7,*) 'Quad entrance/mid/exit',sSTOP_quad_in,sSTOP_quad_mid,sSTOP_quad_out
	  write(7,*) 'D1 entrance/exit      ',sSTOP_bm01_in,sSTOP_bm01_out
	  write(7,*) 'D2 entrance/exit      ',sSTOP_bm02_in,sSTOP_bm02_out
	  write(7,*) 'Vacuum exit           ',sSTOP_exit
	  write(7,*) 'Events reaching hut   ',sSTOP_hut
	  write(7,*) 'DC1, DC2, Scin, Cal   ',sSTOP_dc1,sSTOP_dc2,sSTOP_scin,sSTOP_cal
	  write(7,*) 'Successes             ',sSTOP_successes
	  write(7,*) 
	endif
	if (electron_arm.eq.3 .or. hadron_arm.eq.3) then
	  write(7,*) 'HRSr Trials:          ',rSTOP_trials
	  write(7,*) 'Slit hor/vert         ',rSTOP_slit_vert,rSTOP_slit_hor
	  write(7,*) 'Q1 entrance/mid/exit  ',rSTOP_Q1_in,rSTOP_Q1_mid,rSTOP_Q1_out
	  write(7,*) 'Q2 entrance/mid/exit  ',rSTOP_Q2_in,rSTOP_Q2_mid,rSTOP_Q2_out
	  write(7,*) 'Dipole entrance/exit  ',rSTOP_D1_in,rSTOP_D1_out
	  write(7,*) 'Q3 entrance/mid/exit  ',rSTOP_Q3_in,rSTOP_Q3_mid,rSTOP_Q3_out
	  write(7,*) 'Events reaching hut   ',rSTOP_hut
	  write(7,*) 'VDC1, VDC2            ',rSTOP_dc1,rSTOP_dc2
	  write(7,*) 'S1, S2, Cal	    ',rSTOP_s1,rSTOP_s2,rSTOP_cal
	  write(7,*)
	endif
	if (electron_arm.eq.4 .or. hadron_arm.eq.4) then
	  write(7,*) 'HRSl Trials:          ',lSTOP_trials
	  write(7,*) 'Slit hor/vert         ',lSTOP_slit_vert,lSTOP_slit_hor
	  write(7,*) 'Q1 entrance/mid/exit  ',lSTOP_Q1_in,lSTOP_Q1_mid,lSTOP_Q1_out
	  write(7,*) 'Q2 entrance/mid/exit  ',lSTOP_Q2_in,lSTOP_Q2_mid,lSTOP_Q2_out
	  write(7,*) 'Dipole entrance/exit  ',lSTOP_D1_in,lSTOP_D1_out
	  write(7,*) 'Q3 entrance/mid/exit  ',lSTOP_Q3_in,lSTOP_Q3_mid,lSTOP_Q3_out
	  write(7,*) 'Events reaching hut   ',lSTOP_hut
	  write(7,*) 'VDC1, VDC2            ',lSTOP_dc1,lSTOP_dc2
	  write(7,*) 'S1, S2, Cal	    ',lSTOP_s1,lSTOP_s2,lSTOP_cal
	endif
	if (electron_arm.eq.5 .or. hadron_arm.eq.5 .or.
     >	    electron_arm.eq.6 .or. hadron_arm.eq.6) then
	  write(7,*) 'SHMS Trials:          ',shmsSTOP_trials
	  write(7,*) 'HB phys entrance/mag entr/mag exit/phys exit  ',shmsSTOP_hb_in,shmsSTOP_hb_men,shmsSTOP_hb_mex,shmsSTOP_hb_out
	  write(7,*) 'Slit hor/vert/corners ',shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct
	  write(7,*) 'Q1 phys entrance/mag entr/mid/mag exit/phys exit  '
     >,shmsSTOP_q1_in,shmsSTOP_q1_men,shmsSTOP_q1_mid,shmsSTOP_q1_mex,shmsSTOP_q1_out
	  write(7,*) 'Q2 phys entrance/mag entr/mid/mag exit/phys exit  '
     >,shmsSTOP_q2_in,shmsSTOP_q2_men,shmsSTOP_q2_mid,shmsSTOP_q2_mex,shmsSTOP_q2_out
	  write(7,*) 'Q3 phys entrance/mag entr/mid/mag exit/phys exit       '
     >,shmsSTOP_q3_in,shmsSTOP_q3_men,shmsSTOP_q2_mid,shmsSTOP_q3_mex,shmsSTOP_q3_out
	  write(7,*) 'D1 entrance/flare/mid 1-2   ',shmsSTOP_D1_in,shmsSTOP_D1_flr,shmsSTOP_D1_mid1,shmsSTOP_D1_mid2
	  write(7,*) 'D1 mid 3-5            ',shmsSTOP_D1_mid3,shmsSTOP_D1_mid4,shmsSTOP_D1_mid5
	  write(7,*) 'D1 mid 6-7/mag exit/phys exit         ',shmsSTOP_D1_mid6,shmsSTOP_D1_mid7,shmsSTOP_D1_mex,shmsSTOP_D1_out
c	  write(7,*) 'BP thingie in/out     ',shmsSTOP_BP_in,shmsSTOP_BP_out
	  write(7,*) 'Events reaching hut   ',shmsSTOP_hut
	  write(7,*) 'DC1, DC2, Scin, Cal   ',shmsSTOP_dc1,shmsSTOP_dc2
	  write(7,*) 'S1, S2, S3, Cal       ',shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal
	  write(7,*) 'Successes             ',shmsSTOP_successes
	  write(7,*)
	endif
	if (electron_arm.eq.7 .or. hadron_arm.eq.7 .or. electron_arm.eq. 8 .or. hadron_arm.eq.8) then
	  write(7,*) 'Calo Trials:          ',caloSTOP_trials
	  write(7,*) 'Extent hor/vert       ',caloSTOP_slit_hor,caloSTOP_slit_vert
	  write(7,*) 'Successes             ',caloSTOP_successes
	  endif

	close(7)

	open(unit=7, file=genfile, status='unknown')

	write(7,*) 'E arm Experimental Target Distributions:'
	write(7,'(6a12)') 'delta','EXPERIM','yptar','EXPERIM','xptar','EXPERIM'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%RECON%e%delta%min+(i-0.5)*H%RECON%e%delta%bin,
     >		H%RECON%e%delta%buf(i), H%RECON%e%yptar%min+(i-0.5)*
     >		H%RECON%e%yptar%bin, H%RECON%e%yptar%buf(i),
     >		H%RECON%e%xptar%min+(i-0.5)*H%RECON%e%xptar%bin,
     >		H%RECON%e%xptar%buf(i)
	enddo

	write(7,*) 'P arm Experimental Target Distributions:'
	write(7,'(6a12)') 'delta','EXPERIM','yptar','EXPERIM','xptar','EXPERIM'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%RECON%p%delta%min+(i-0.5)*H%RECON%p%delta%bin,
     >		H%RECON%p%delta%buf(i), H%RECON%p%yptar%min+(i-0.5)*
     >		H%RECON%p%yptar%bin, H%RECON%p%yptar%buf(i),
     >		H%RECON%p%xptar%min+(i-0.5)*H%RECON%p%xptar%bin,
     >		H%RECON%p%xptar%buf(i)
	enddo

	write(7,*) 'Distributions of Contributing E arm Events:'
	write(7,'(6a12)') 'delta','CONTRIB','yuptar','CONTRIB','xptar','CONTRIB'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%gen%e%delta%min+(i-0.5)*H%gen%e%delta%bin,
     >		H%gen%e%delta%buf(i), H%gen%e%yptar%min+(i-0.5)*
     >		H%gen%e%yptar%bin, H%gen%e%yptar%buf(i),
     >		H%gen%e%xptar%min+(i-0.5)*H%gen%e%xptar%bin, H%gen%e%xptar%buf(i)
	enddo

	write(7,*) 'Distributions of Contributing P arm Events:'
	write(7,'(6a12)') 'delta','CONTRIB','yptar','CONTRIB','xptar','CONTRIB'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%gen%p%delta%min+(i-0.5)*H%gen%p%delta%bin,
     >		H%gen%p%delta%buf(i), H%gen%p%yptar%min+(i-0.5)*
     >		H%gen%p%yptar%bin, H%gen%p%yptar%buf(i),
     >		H%gen%p%xptar%min+(i-0.5)*H%gen%p%xptar%bin, H%gen%p%xptar%buf(i)
	enddo

	write(7,*) 'Original E arm Events:'
	write(7,'(6a12)') 'delta','ORIGIN','yptar','ORIGIN','xptar','ORIGIN'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%gen%e%delta%min+(i-0.5)*H%gen%e%delta%bin,
     >		H%gen%e%delta%buf(i), H%gen%e%yptar%min+(i-0.5)*
     >		H%gen%e%yptar%bin, H%geni%e%yptar%buf(i),
     >		H%gen%e%xptar%min+(i-0.5)*H%gen%e%xptar%bin, H%geni%e%xptar%buf(i)
	enddo

	write(7,*) 'Original P arm Events:'
	write(7,'(6a12)') 'delta','ORIGIN','yptar','ORIGIN','xptar','ORIGIN'
	do i=1,nHbins
	  write(7,'(3(1x,2(e11.4,1x)))')
     >		H%geni%p%delta%min+(i-0.5)*H%geni%p%delta%bin,
     >		H%geni%p%delta%buf(i), H%geni%p%yptar%min+(i-0.5)*
     >		H%geni%p%yptar%bin, H%geni%p%yptar%buf(i),
     >		H%geni%p%xptar%min+(i-0.5)*H%geni%p%xptar%bin, H%geni%p%xptar%buf(i)
	enddo

	write(7,*) 'Original Em/Pm distributions:'
	write(7,'(4a12)') 'Em','ORIGIN','Pm','ORIGIN'
	do i=1,nHbins
	  write(7,'(3(1x,4(e11.4,1x)))')
     >		H%geni%Em%min+(i-0.5)*H%geni%Em%bin,
     >		H%geni%Em%buf(i), H%geni%Pm%min+(i-0.5)*
     >		H%geni%Pm%bin, H%geni%Pm%buf(i)
	enddo

	close(7)

! Counters, etc report to the screen and to the diagnostic histogram file
!	call report(6,timestring1,timestring2,central,contrib,sum_sigcc,aveerr,resol)
	open(unit=7,file=histfile,status='unknown')
	call report(7,timestring1,timestring2,central,contrib,sum_sigcc,aveerr,resol)
	close(unit=7)

1000	end

!-------------------------------------------------------------------

	subroutine inc(hist,val,weight)

	implicit none
	include 'histograms.inc'

	integer*4		ibin
	real*8			val, weight
	type(hist_entry)::	hist

	ibin= nint(0.5+(val-hist%min)/hist%bin)
	if(ibin.ge.1.and.ibin.le.nHbins)then
	hist%buf(ibin) = hist%buf(ibin) + weight
	endif

	return
	end

!--------------------------------------------------------------------

	subroutine report(iun,timestring1,timestring2,central,contrib,
     >		sum_sigcc,aveerr,resol)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'
	include 'brem.inc'

	integer*4	iun
	real*8		sum_sigcc
	character*30	timestring1, timestring2
	type(contribtype)::    contrib
	type(event_central):: central
	type(sums_twoarm):: aveerr, resol

! Output of diagnostics

! Run time
	write(iun,'(/1x,''BEGIN Time: '',a30)') timestring1
	write(iun,'(1x,''END Time:   '',a30)') timestring2

! Kinematics
	write(iun,*) 'KINEMATICS:'
	if (doing_eep) then
	  if (doing_hyd_elast) then
	    write(iun,*) '              ****--------  H(e,e''p)  --------****'
	  else if (doing_deuterium) then
	    write(iun,*) '              ****--------  D(e,e''p)  --------****'
	  else if (doing_heavy) then
	    write(iun,*) '              ****--------  A(e,e''p)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''p) we''re doing!!!'
	  endif
	else if (doing_semi) then 
	   if (doing_semipi) then 
	      if (targ%A .eq. 1) then 
		 if(doing_hplus) then
		    write(iun,*) ' ****--------  H(e,e''pi+)X  --------****'
		 else
		    write(iun,*) ' ****--------  H(e,e''pi-)X  --------****'
		 endif
	      elseif (targ%A .eq. 2) then
		 if(doing_hplus) then
		    write(iun,*) ' ****--------  D(e,e''pi+)X  --------****'
		 else
		    write(iun,*) ' ****--------  D(e,e''pi-)X  --------****'
		 endif
	      else
		 stop 'I don''t have ANY idea what A(e,e''pi)X we''re doing!!!'
	      endif
	   else if (doing_semika) then  
	      if (targ%A .eq. 1) then 
		 if(doing_hplus) then
		    write(iun,*) ' ****--------  H(e,e''k+)X  --------****'
		 else
		    write(iun,*) ' ****--------  H(e,e''k-)X  --------****'
		 endif
	      elseif (targ%A .eq. 2) then
		 if(doing_hplus) then
		    write(iun,*) ' ****--------  D(e,e''k+)X  --------****'
		 else
		    write(iun,*) ' ****--------  D(e,e''k-)X  --------****'
		 endif
	      else
		 stop 'I don''t have ANY idea what A(e,e''k)X we''re doing!!!'
	      endif
	   else   
	      stop 'I don''t have ANY idea what A(e,e''x)X we''re doing!!!'
          endif         
	else if (doing_rho) then
	   if (targ%A .eq. 1) then
	      write(iun,*) '              ****--------  H(e,e''rho)  --------****'
	   else
	      write(iun,*) 'I am not set up for anything else yet!'
	   endif
	else if (doing_delta) then
	  if (doing_hyddelta) then
	    write(6,*) ' ****--------  H(e,e''p)pi  --------****'
	  else if (doing_deutdelta) then
	    write(6,*) ' ****--------  D(e,e''p)pi  --------****'
	  else if (doing_hedelta) then
	    write(6,*) ' ****--------  A(e,e''p)pi  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''p)pi we''re doing!!!'
	  endif
	else if (doing_pion) then
	  if (doing_hydpi) then
	    if (targ%A .eq. 1) then
	      write(iun,*) '              ****--------  H(e,e''pi)  --------****'
	    else if (targ%A .ge.3) then
	      write(iun,*) '              ****--------  A(e,e''pi)  --------****'
	    endif
	  else if (doing_deutpi) then
	    write(iun,*) '              ****--------  D(e,e''pi)  --------****'
	  else if (doing_hepi) then
	    write(iun,*) '              ****--------  A(e,e''pi)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''pi) we''re doing!!!'
	  endif
	  if (which_pion .eq. 0 .or. which_pion .eq. 1) then
	    write(iun,*) '              ****----  Default Final State ----****'
	  else if (which_pion .eq. 10) then
	    write(iun,*) '              ****----  Final State is A + pi ----****'
	  else if (which_pion.eq.2 .or. which_pion.eq.3) then
	    write(iun,*) '              ****----  Final State is pi + Delta ----****'
	  endif
	else if (doing_kaon) then
	  if (doing_hydkaon) then
	    if (targ%A .eq. 1) then
	      write(iun,*) '              ****--------  H(e,e''K)  --------****'
	    else if (targ%A .ge.3) then
	      write(iun,*) '              ****--------  A(e,e''K)  --------****'
	    endif
	  else if (doing_deutkaon) then
	    write(iun,*) '              ****--------  D(e,e''K)  --------****'
	  else if (doing_hekaon) then
	    write(iun,*) '              ****--------  A(e,e''K)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	  if (which_kaon.eq.0) then
	    write(iun,*) '              ****---- producing a LAMBDA ----****'
	  else if (which_kaon.eq.1) then
	    write(iun,*) '              ****---- producing a SIGMA0 ----****'
	  else if (which_kaon.eq.2) then
	    write(iun,*) '              ****---- producing a SIGMA- ----****'
	  else if (which_kaon.eq.10) then
	    write(iun,*) '              ****---- WITH BOUND LAMBDA  ----****'
	  else if (which_kaon.eq.11) then
	    write(iun,*) '              ****---- WITH BOUND SIGMA0  ----****'
	  else if (which_kaon.eq.12) then
	    write(iun,*) '              ****---- WITH BOUND SIGMA-  ----****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	else if (doing_phsp) then
	  write(iun,*) '              ****--- PHASE SPACE - NO physics, NO radiation (may not work)---****'
	else
	  stop 'I don''t have ANY idea what we''re doing!!!'
	endif

	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'Ebeam', Ebeam, 'MeV'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)')
     >		'(dE/E)beam', dEbeam/Ebeam, '(full wid)'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'x-width', gen%xwid, 'cm'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'y-width', gen%ywid, 'cm'
	write(iun,'(9x,a12,'' = '',i15,2x,a16)')
     >		'fr_pattern',targ%fr_pattern, '1=square,2=circ'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'fr1', targ%fr1, 'cm'
	write(iun,'(9x,a12,'' = '',f15.4,2x,a10)') 'fr2', targ%fr2, 'cm'

	write(iun,*) ' '
	write(iun,'(9x,18x,2(a15,2x))') '____E arm____','____P arm____'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'angle',
     >		spec%e%theta*degrad, spec%p%theta*degrad, 'deg'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'momentum',
     >		spec%e%p, spec%p%p, 'MeV/c'

	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'x offset',
     >		spec%e%offset%x, spec%p%offset%x, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'y offset',
     >		spec%e%offset%y, spec%p%offset%y, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'z offset',
     >		spec%e%offset%z, spec%p%offset%z, 'cm'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'xptar offset',
     >		spec%e%offset%xptar, spec%p%offset%xptar, 'mr'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'yptar offset',
     >		spec%e%offset%yptar, spec%p%offset%yptar, 'mr'

	write(iun,*) '                      VALUES FOR "CENTRAL" EVENT:'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'delta',
     >		central%e%delta, central%p%delta, '%'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'xptar',
     >		central%e%xptar, central%p%xptar, 'mr'
	write(iun,'(9x,a12,'' = '',2(f15.4,2x),2x,a5)') 'yptar',
     >		central%e%yptar, central%p%yptar, 'mr'

	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'Q2',central%Q2/1.e6,'(GeV/c)^2'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'q', central%q, 'MeV/c'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'nu',
     >		central%nu, 'MeV'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'recon Em', central%Em, 'MeV'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'recon Pm', central%Pm, 'MeV/c'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'recon  W', central%W, 'MeV/c'
	write(iun,'(17x,a10,'' = '',f15.4,2x,a9)') 'recon MM', central%MM, 'MeV/c'

! Target
	write(iun,*) 'TARGET specs:'
9911	format(2x,2(5x,a10,' = ',e12.6,1x,a5))
	write(iun,9911) 'A', targ%A, ' ', 'Z', targ%Z, ' '
	write(iun,9911) 'mass', targ%mass_amu, 'amu', 'mass', targ%M,'MeV'
	write(iun,9911) 'Mrec', targ%mrec_amu, 'amu', 'Mrec', targ%Mrec,'MeV'
	write(iun,9911) 'Mtar_struc',targ%Mtar_struck,'MeV',
     >		'Mrec_struc',targ%Mrec_struck,'MeV'
	write(iun,9911) 'rho', targ%rho, 'g/cm3', 'thick', targ%thick,'g/cm2'
	write(iun,9911) 'angle', targ%angle*degrad, 'deg', 'abundancy',
     >		targ%abundancy, '%'
	write(iun,9911) 'X0', targ%X0, 'g/cm2', 'X0_cm', targ%X0_cm,'cm'
	write(iun,9911) 'length',targ%length,'cm','zoffset',targ%zoffset,'cm'
	write(iun,9911) 'xoffset',targ%xoffset,'cm','yoffset',targ%yoffset,'cm'
	write(iun,'(t12,3a15)') '__ave__', '__lo__', '__hi__'
9912	format(1x,a15,3f15.5,2x,a6)
	write(iun,9912) 'Coulomb', targ%Coulomb%ave, targ%Coulomb%min,
     >		targ%Coulomb%max, 'MeV'
	write(iun,9912) 'Eloss_beam', targ%Eloss(1)%ave,
     >		targ%Eloss(1)%min, targ%Eloss(1)%max, 'MeV'
	write(iun,9912) 'Eloss_e', targ%Eloss(2)%ave, targ%Eloss(2)%min,
     >		targ%Eloss(2)%max, 'MeV'
	write(iun,9912) 'Eloss_p', targ%Eloss(3)%ave, targ%Eloss(3)%min,
     >		targ%Eloss(3)%max, 'MeV'
	write(iun,9912) 'teff_beam', targ%teff(1)%ave, targ%teff(1)%min,
     >		targ%teff(1)%max, 'radlen'
	write(iun,9912) 'teff_e', targ%teff(2)%ave, targ%teff(2)%min,
     >		targ%teff(2)%max, 'radlen'
	write(iun,9912) 'teff_p', targ%teff(3)%ave, targ%teff(3)%min,
     >		targ%teff(3)%max, 'radlen'
9913	format(1x,a15,t25,f15.5,2x,a6)
	write(iun,9913) 'musc_nsig_max', targ%musc_nsig_max, ' '
	write(iun,9913) 'musc_max_beam', targ%musc_max(1)*1000., 'mr'
	write(iun,9913) 'musc_max_e', targ%musc_max(2)*1000., 'mr'
	write(iun,9913) 'musc_max_p', targ%musc_max(3)*1000., 'mr'

! Flags
	write(iun,*) 'FLAGS:'
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_eep', doing_eep,
     >		'doing_kaon', doing_kaon, 'doing_pion', doing_pion
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_semi', doing_semi,
     >		'doing_rho', doing_rho, 'doing_hplus', doing_hplus
	write(iun,'(5x,2(2x,a19,''='',l2))') 'doing_semipi',doing_semipi,
     >		'doing_semika', doing_semika
	write(iun,'(5x,2(2x,a19,''='',l2))') 'doing_delta',doing_delta,
     >		'doing_phsp', doing_phsp
	write(iun,'(5x,2(2x,a19,''='',i2))') 'which_pion', which_pion,
     >		'which_kaon', which_kaon
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hyd_elast', doing_hyd_elast,
     >		'doing_deuterium', doing_deuterium, 'doing_heavy', doing_heavy
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydpi', doing_hydpi,
     >          'doing_deutpi', doing_deutpi, 'doing_hepi', doing_hepi
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydkaon', doing_hydkaon,
     >		'doing_deutkaon', doing_deutkaon, 'doing_hekaon', doing_hekaon
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydsemi', doing_hydsemi,
     >          'doing_deutsemi', doing_deutsemi, 'do_fermi', do_fermi
	write(iun,'(5x,3(2x,a19,''='',l2))') 'doing_hydrho', doing_hydrho,
     >          'doing_deutrho', doing_deutrho, 'doing_herho', doing_herho
	write(iun,'(5x,(2x,a19,''='',l2),2(2x,a19,''='',i2))') 'mc_smear',
     >		mc_smear,'electron_arm',electron_arm,'hadron_arm',hadron_arm
	write(iun,'(5x,3(2x,a19,''='',l2)))') 'using_Eloss', using_Eloss,
     >		'using_Coulomb',using_Coulomb,'deForest_flag',deForest_flag
	write(iun,'(5x,3(2x,a19,''='',l2)))') 'correct_Eloss', correct_Eloss,
     >		'correct_raster',correct_raster, 'doing_decay', doing_decay
	write(iun,'(5x,3(2x,a19,''='',l2))') 
     >		'using_E_arm_montecarlo', using_E_arm_montecarlo,
     >		'using_P_arm_montecarlo', using_P_arm_montecarlo,
     >          'use_benhar_sf', use_benhar_sf
	if (electron_arm.eq.5 .or. hadron_arm.eq.5 .or.
     >	    electron_arm.eq.6 .or. hadron_arm.eq.6)
     >	    write(iun,'(7x,a19,''='',l2)') 'use_first_cer',use_first_cer
	write(iun,'(7x,a11,''='',f10.3,a4)') 'ctau',ctau,'cm'
	if (use_benhar_sf)
     >      write(iun,'(7x,a12,''='',f8.4)') 'transparency',transparency
! Counters
	write(iun,*) 'COUNTERS:'
	write(iun,'(12x,''Ngen (request) = '',i10)') ngen
	write(iun,'(12x,''Ntried         = '',i10)') ntried
	write(iun,'(12x,''Ncontribute    = '',i10)') ncontribute
	write(iun,'(12x,''Nco_no_rad_prot= '',i10)') ncontribute_no_rad_proton
	write(iun,'(12x,''-> %no_rad_prot= '',f10.3)')
     >		(100.*ncontribute_no_rad_proton/max(dble(ncontribute),0.1e0))
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
	  write(iun,9914) 'hardcorfac', central%rad%hardcorfac
	  write(iun,9914) 'etatzai', central%rad%etatzai
	  write(iun,9914) 'frac(1:3)', central%rad%frac
	  write(iun,9914) 'lambda(1:3)', central%rad%lambda
	  write(iun,9914) 'bt(1:2)', central%rad%bt
	  write(iun,9914) 'c_int(0:3)', central%rad%c_int
	  write(iun,9914) 'c_ext(0:3)', central%rad%c_ext
	  write(iun,9914) 'c(0:3)', central%rad%c
	  write(iun,9914) 'g_int', central%rad%g_int
	  write(iun,9914) 'g_ext', central%rad%g_ext
	  write(iun,9914) 'g(0:3)', central%rad%g
	endif

! Miscellaneous
	write(iun,*) 'MISCELLANEOUS:'
9915	format(12x,a14,' = ',e16.6,1x,a6)
	write(iun,*) 'Note that central.sigcc is for central delta,theta,phi in both spectrometers'
	write(iun,*) ' and may give non-physical kinematics, esp. for Hydrogen'
	write(iun,*) 'Note also that AVE.sigcc is really AVER.weight (the two arenot exactly equal)'
	write(iun,9915) 'CENTRAL.sigcc', central%sigcc, ' '
	write(iun,9915) 'AVERAGE.sigcc', sum_sigcc/nevent, ' '
	write(iun,9915) 'charge', EXPER%charge, 'mC'
	write(iun,9915) 'targetfac', targetfac, ' '
	write(iun,9915) 'luminosity', luminosity, 'ub^-1'
	write(iun,9915) 'luminosity', luminosity*(hbarc/100000.)**2, 'GeV^2'
	write(iun,9915) 'genvol', genvol, ' '
	write(iun,9915) 'normfac', normfac, ' '
9916	format(12x,a14,' = ',f6.1,' to ',f6.1,' in ',i4,' bins of ',f7.2,1x,a5)
	if (doing_heavy) then
	  write(iun,'(12x,''Theory file:  '',a)')
     >		theory_file(1:index(theory_file,' ')-1)
	endif
	if(random_seed.ne.0) write(iun,'(15x,''Random Seed = '',i10)') random_seed
	if(random_state_file.ne.' ') write(iun,'( '' Random State Save File:  '',a)')
     >         random_state_file(1:index(random_state_file,' ')-1)

! Resolution summary

	write(iun,*) 'RECON SUMMARY:            Ave.Error   Resolution'
        write(iun,99165) 'Electron arm: delta =',10.*aveerr%e%delta,10.*resol%e%delta,'x10^-3'
        write(iun,99165) 'xptar =',1000.*aveerr%e%xptar,1000.*resol%e%xptar,'mr'
        write(iun,99165) 'yptar =',1000.*aveerr%e%yptar,1000.*resol%e%yptar,'mr'
        write(iun,99165) 'ytar  =',10.*aveerr%e%ytar,10.*resol%e%ytar ,'mm'
        write(iun,99165) 'Hadron arm:   delta =',10.*aveerr%p%delta,10.*resol%p%delta,'x10-3'
        write(iun,99165) 'xptar =',1000.*aveerr%p%xptar,1000.*resol%p%xptar,'mr'
        write(iun,99165) 'yptar =',1000.*aveerr%p%yptar,1000.*resol%p%yptar,'mr'
        write(iun,99165) 'ytar  =',10.*aveerr%p%ytar,10.*resol%p%ytar,'mm'
99165    format(2x,a22,2f12.5,2x,a)


! Input spectrometer limits.

	write(iun,*) 'Input Spectrometer Limits:'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.delta.min/max',
     >		SPedge%e%delta%min,SPedge%e%delta%max,'%'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.yptar.min/max',
     >		SPedge%e%yptar%min,SPedge%e%yptar%max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.e.xptar.min/max',
     >		SPedge%e%xptar%min,SPedge%e%xptar%max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.delta.min/max',
     >		SPedge%p%delta%min,SPedge%p%delta%max,'%'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.yptar.min/max',
     >		SPedge%p%yptar%min,SPedge%p%yptar%max,'rad'
	write(iun,'(9x,a25,'' = '',2(2x,f15.4),a5)') 'SPedge.p.xptar.min/max',
     >		SPedge%p%xptar%min,SPedge%p%xptar%max,'rad'


! Edges used in generation and checking, as well as range of events found
! to contribute within those edges

! ... on GENERATED qties
	write(iun,*) 'Limiting VERTEX values (vertex.e/p.*,Em,Pm,Trec)'
	write(iun,*) '   USED limits are gen.e/p.*, and VERTEXedge.Em,Pm,Trec'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)') '______used______',
     >		'_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
9917	format(1x,a18,t21,2f12.3,t50,2f10.3,2x,a5)
	write(iun,9917) 'E arm  delta', gen%e%delta%min, gen%e%delta%max,
     >		contrib%gen%e%delta%lo, contrib%gen%e%delta%hi, '%'
	write(iun,9917) 'E arm  yptar', gen%e%yptar%min*1000.,
     >		gen%e%yptar%max*1000.,
     >		contrib%gen%e%yptar%lo*1000., contrib%gen%e%yptar%hi*1000.,'mr'
	write(iun,9917) 'E arm  xptar', gen%e%xptar%min*1000.,
     >		gen%e%xptar%max*1000.,
     >		contrib%gen%e%xptar%lo*1000., contrib%gen%e%xptar%hi*1000.,'mr'
	write(iun,9917) 'P arm  delta', gen%p%delta%min, gen%p%delta%max,
     >		contrib%gen%p%delta%lo, contrib%gen%p%delta%hi, '%'
	write(iun,9917) 'P arm  yptar', gen%p%yptar%min*1000.,
     >		gen%p%yptar%max*1000.,
     >		contrib%gen%p%yptar%lo*1000., contrib%gen%p%yptar%hi*1000.,'mr'
	write(iun,9917) 'P arm  xptar', gen%p%xptar%min*1000.,
     >		gen%p%xptar%max*1000.,
     >		contrib%gen%p%xptar%lo*1000., contrib%gen%p%xptar%hi*1000.,'mr'
	write(iun,9917) 'sumEgen', gen%sumEgen%min, gen%sumEgen%max,
     >		contrib%gen%sumEgen%lo, contrib%gen%sumEgen%hi, 'MeV'

! ... on VERTEX qties
!	write(iun,*) 'Limiting VERTEX values:'
!	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
!     >		'______used______', '_____found______'
!	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'Trec', VERTEXedge%Trec%min, VERTEXedge%Trec%max,
     >      contrib%vertex%Trec%lo, contrib%vertex%Trec%hi, 'MeV'
	write(iun,9917) 'Em', VERTEXedge%Em%min, VERTEXedge%Em%max,
     >      contrib%vertex%Em%lo, contrib%vertex%Em%hi, 'MeV'
	write(iun,9917) 'Pm', VERTEXedge%Pm%min, VERTEXedge%Pm%max,
     >       contrib%vertex%Pm%lo, contrib%vertex%Pm%hi, 'MeV/c'
	if ((doing_deuterium .or. doing_pion .or. doing_kaon .or. doing_delta) .and. using_rad) then
	   write(iun,*) '      *** NOTE: sumEgen.min only used in GENERATE_RAD'
	endif

! ... on TRUE qties
	write(iun,*) 'Limiting ORIGINAL values: orig.e/p.*,Em,Pm,Trec (no edge.* limits for Pm,Trec)'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'E arm   E', edge%e%E%min, edge%e%E%max,
     >		contrib%tru%e%E%lo, contrib%tru%e%E%hi, 'MeV'
	write(iun,9917) 'E arm  yptar', edge%e%yptar%min*1000.,
     >		edge%e%yptar%max*1000.,
     >		contrib%tru%e%yptar%lo*1000., contrib%tru%e%yptar%hi*1000.,'mr'
	write(iun,9917) 'E arm  xptar', edge%e%xptar%min*1000., edge%e%xptar%max*1000.,
     >		contrib%tru%e%xptar%lo*1000., contrib%tru%e%xptar%hi*1000., 'mr'
	write(iun,9917) 'P arm      E', edge%p%E%min, edge%p%E%max,
     >		contrib%tru%p%E%lo, contrib%tru%p%E%hi, 'MeV'
	write(iun,9917) 'P arm  yptar', edge%p%yptar%min*1000.,
     >		edge%p%yptar%max*1000.,
     >		contrib%tru%p%yptar%lo*1000., contrib%tru%p%yptar%hi*1000.,'mr'
	write(iun,9917) 'P arm  xptar', edge%p%xptar%min*1000., edge%p%xptar%max*1000.,
     >		contrib%tru%p%xptar%lo*1000., contrib%tru%p%xptar%hi*1000., 'mr'
	write(iun,9917) 'Em', max(-999999.999e0,edge%Em%min), min(999999.999e0,edge%Em%max),
     >		contrib%tru%Em%lo, contrib%tru%Em%hi, 'MeV'
	write(iun,9917) 'Pm', 0., 0., contrib%tru%Pm%lo,
     >		contrib%tru%Pm%hi, 'MeV'
	write(iun,9917) 'Trec', 0., 0., contrib%tru%Trec%lo,
     >		contrib%tru%Trec%hi, 'MeV'

! ... on SPECTROMETER qties
	write(iun,*) 'Limiting SPECTROMETER values'
	write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	write(iun,9917) 'E arm delta', SPedge%e%delta%min,
     >		SPedge%e%delta%max, contrib%SP%e%delta%lo,
     >		contrib%SP%e%delta%hi, '%'
	write(iun,9917) 'E arm  yptar', SPedge%e%yptar%min*1000.,
     >		SPedge%e%yptar%max*1000.,
     >		contrib%SP%e%yptar%lo*1000., contrib%SP%e%yptar%hi*1000., 'mr'
	write(iun,9917) 'E arm  xptar', SPedge%e%xptar%min*1000.,
     >		SPedge%e%xptar%max*1000.,
     >		contrib%SP%e%xptar%lo*1000., contrib%SP%e%xptar%hi*1000., 'mr'
	write(iun,9917) 'P arm  delta', SPedge%p%delta%min,
     >		SPedge%p%delta%max, contrib%SP%p%delta%lo,
     >		contrib%SP%p%delta%hi, '%'
	write(iun,9917) 'P arm  yptar', SPedge%p%yptar%min*1000.,
     >		SPedge%p%yptar%max*1000.,
     >		contrib%SP%p%yptar%lo*1000., contrib%SP%p%yptar%hi*1000., 'mr'
	write(iun,9917) 'P arm  xptar', SPedge%p%xptar%min*1000.,
     >		SPedge%p%xptar%max*1000.,
     >		contrib%SP%p%xptar%lo*1000., contrib%SP%p%xptar%hi*1000., 'mr'

! ... on RADIATION qties
	if (using_rad) then
	  write(iun,*) 'Limiting RADIATION values CONTRIBUTING to the (Em,Pm) distributions:'
	  write(iun,'(t25,2x,a16,2x,t50,2x,a16,2x)')
     >		'______used______', '_____found______'
	  write(iun,'(t25,2a10,t50,2a10)') 'min','max', 'lo', 'hi'
	  write(iun,9917) 'Egamma(1)', 0., Egamma1_max,
     >		contrib%rad%Egamma(1)%lo, contrib%rad%Egamma(1)%hi, 'MeV'
	  write(iun,9917) 'Egamma(2)', 0., Egamma2_max, contrib%rad%Egamma(2)%lo,
     >		contrib%rad%Egamma(2)%hi, 'MeV'
	  write(iun,9917) 'Egamma(3)', 0., Egamma3_max, contrib%rad%Egamma(3)%lo,
     >		contrib%rad%Egamma(3)%hi, 'MeV'
	  write(iun,9917) 'Egamma_total', 0., Egamma_tot_max,
     >		contrib%rad%Egamma_total%lo, contrib%rad%Egamma_total%hi,'MeV'
	endif

! ... on slops
	write(iun,*) 'ACTUAL and LIMITING SLOP values used/obtained:'
	write(iun,'(t25,3a10)') '__used__', '__min__', '__max__'
9918	format(1x,a10,a12,t25,3f10.3,2x,a5)
	write(iun,9918) 'slop.MC  ', 'E arm delta', slop%MC%e%delta%used,
     >		slop%MC%e%delta%lo, slop%MC%e%delta%hi, '%'
	write(iun,9918) ' ', 'E arm yptar', slop%MC%e%yptar%used*1000.,
     >		slop%MC%e%yptar%lo*1000., slop%MC%e%yptar%hi*1000., 'mr'
	write(iun,9918) ' ', 'E arm xptar', slop%MC%e%xptar%used*1000.,
     >		slop%MC%e%xptar%lo*1000., slop%MC%e%xptar%hi*1000., 'mr'
	write(iun,9918) ' ', 'P arm delta', slop%MC%p%delta%used,
     >		slop%MC%p%delta%lo, slop%MC%p%delta%hi, '%'
	write(iun,9918) ' ', 'P arm yptar', slop%MC%p%yptar%used*1000.,
     >		slop%MC%p%yptar%lo*1000., slop%MC%p%yptar%hi*1000., 'mr'
	write(iun,9918) ' ', 'P arm xptar', slop%MC%p%xptar%used*1000.,
     >		slop%MC%p%xptar%lo*1000., slop%MC%p%xptar%hi*1000., 'mr'
	write(iun,9918) 'slop.total', 'Em', slop%total%Em%used,
     >		slop%total%Em%lo, slop%total%Em%hi, 'MeV'
	write(iun,9918) ' ', 'Pm', 0., slop%total%Pm%lo,
     >		slop%total%Pm%hi, 'MeV/c'

	write(iun,'(/)')
	return
	end

!-----------------------------------------------------------------------

	subroutine calculate_central(central,vertex0)

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	integer i
	real*8 m_spec
	logical success
	type(event_main)::	main0
	type(event)::		vertex0, recon0
	type(event_central)::	central

! JRA 2002:  Redo this so that central values are chosen as in generate,
! and then call complete_ev and complete_recon_ev, just like a normal event.

	if (debug(2)) write(6,*)'calc_cent: entering...'
	main0%target%x = 0.0
	main0%target%y = 0.0
	main0%target%z = targ%zoffset
	main0%target%rastery = 0.0
	main0%target%Eloss(1) = targ%Eloss(1)%ave
	main0%target%Coulomb = targ%Coulomb%ave
	main0%target%teff(1) = targ%teff(1)%ave
	vertex0%Ein = Ebeam_vertex_ave
	main0%Ein_shift = 0.0
	main0%Ee_shift = 0.0
	main0%gen_weight = 1.0

! Set all of these to central values, but then complete_ev will force
! the variables that are not normally generated to their appropriate values.
	vertex0%e%delta = 0.0
	vertex0%e%yptar = 0.0
	vertex0%e%xptar = 0.0
	vertex0%e%theta = spec%e%theta
	vertex0%e%phi = spec%e%phi
	vertex0%e%P = spec%e%P*(1.+vertex0%e%delta/100.)
	vertex0%e%E = vertex0%e%P

	vertex0%p%delta = 0.0
	vertex0%p%yptar = 0.0
	vertex0%p%xptar = 0.0
	vertex0%p%theta = spec%p%theta
	vertex0%p%phi = spec%p%phi
	vertex0%p%P = spec%p%P*(1.+vertex0%p%delta/100.)
	vertex0%p%E = sqrt(vertex0%p%P**2 + Mh2)

	pfer = 0.0
	pferx = 0.0
	pfery = 0.0
	pferz = 0.0
	vertex0%Em = targ%Mtar_struck + targ%Mrec - targ%M
	m_spec = targ%M - targ%Mtar_struck + vertex0%Em		!=targ.Mrec
	efer = targ%M - sqrt(m_spec**2+pfer**2)			!=M-Mrec

! Old version - not appropriate for all event types.
! Pop in generation values for central event
!
!	if (debug(2)) write(6,*)'calc_cent: entering...'
!	main0.target.x = 0.0
!	main0.target.y = 0.0
!	main0.target.z = 0.0
!	vertex0.Ein = Ebeam_vertex_ave
!	main0.target.Coulomb = targ.Coulomb.ave
!	main0.target.Eloss(1) = targ.Eloss(1).ave
!	main0.target.teff(1) = targ.teff(1).ave
!	vertex0.e.delta = 0.0
!	vertex0.e.yptar = 0.0
!	vertex0.e.xptar = 0.0
!	vertex0.e.theta = spec.e.theta
!	vertex0.e.phi = spec.e.phi
!	vertex0.e.P = spec.e.P*(1.+vertex0.e.delta/100.)
!	vertex0.e.E = vertex0.e.P
!	vertex0.p.delta = 0.0
!	vertex0.p.yptar = 0.0
!	vertex0.p.xptar = 0.0
!	vertex0.p.theta = spec.p.theta
!	vertex0.p.phi = spec.p.phi
!	vertex0.p.P = spec.p.P*(1.+vertex0.p.delta/100.)
!	vertex0.p.E = sqrt(vertex0.p.P**2 + Mh2)

! Complete_recon_ev for vertex0.   Note: complete_recon_ev doesn't
! call radc_init_ev or calculate teff(2,3).

! JRA: Do we want complete_ev or complete_recon_ev?  Do we want to calculate
! and/or dump other central values (for pion/kaon production, for example).

c	call complete_ev(main0,vertex0,success)
c	if (debug(2)) write(6,*)'calc_cent: done with complete_ev'
c	if (.not.success) stop 'COMPLETE_EV failed trying to complete a CENTRAL event!'

	call complete_recon_ev(vertex0,success)
	if (debug(2)) write(6,*)'calc_cent: done with complete_recon_ev'
	if (.not.success) stop 'COMPLETE_EV failed trying to complete a CENTRAL event!'
	central%Q2 = vertex0%Q2
	central%q = vertex0%q
	central%nu = vertex0%nu
	central%Em = vertex0%Em
	central%Pm = vertex0%Pm
	central%W = vertex0%W
	central%e%xptar = vertex0%e%xptar
	central%e%yptar = vertex0%e%yptar
	central%e%delta = vertex0%e%delta
	central%p%xptar = vertex0%p%xptar
	central%p%yptar = vertex0%p%yptar
	central%p%delta = vertex0%p%delta

	if (central%Em**2-central%Pm**2 .lt. 0) then
	  central%MM = -sqrt(abs(central%Em**2-central%Pm**2))
	else
	  central%MM = sqrt(central%Em**2-central%Pm**2)
	endif

	main0%target%teff(2) = targ%teff(2)%ave
	main0%target%teff(3) = targ%teff(3)%ave
	if (debug(2)) write(6,*)'calc_cent: calling radc_init_ev'
	if (debug(4)) write(6,*)'calc_cent: Ein =',vertex0%Ein
	call radc_init_ev(main0,vertex0)
	if (debug(2)) write(6,*)'calc_cent: done with radc_init_ev'
	if (using_rad) then
	  central%rad%hardcorfac = hardcorfac
	  central%rad%etatzai= etatzai
	  central%rad%g_int = g_int
	  central%rad%g_ext = g_ext
	  do i = 0, 3
	    if (i.gt.0) central%rad%frac(i) = frac(i)
	    if (i.gt.0) central%rad%lambda(i) = lambda(i)
	    if (i.gt.0.and.i.lt.3) central%rad%bt(i) = bt(i)
	    central%rad%c_int(i) = c_int(i)
	    central%rad%c_ext(i) = c_ext(i)
	    central%rad%c(i) = c(i)
	    central%rad%g(i) = g(i)
	  enddo
	endif
	if (debug(4)) write(6,*)'calc_cent: at 1'

c	do i = 1, neventfields
c	  recon0%all(i) = vertex0%all(i)
c	enddo

	recon0 = vertex0
	call complete_main(.true.,main0,vertex0,vertex0,recon0,success)
	central%sigcc = main0%sigcc

	if(doing_semi) then
	   write(6,*) 'central event'
	   write(6,*) 'Pt',sqrt(vertex0%pt2)/1.e3
	   write(6,*) 'z', vertex0%zhad
	   write(6,*) 'lab cross section (nb/Gev2/sr2)',central%sigcc*1000.0*1000.0*1000.0
	endif
	if (debug(2)) write(6,*)'calc_cent: ending...'
	return
	end

!-------------------------------------------------------------------------

	subroutine montecarlo(orig,main,recon,success)

	implicit none
	include 'simulate.inc'

* Track-coordinate and spectrometer common blocks

	real*8 x_E_arm,y_E_arm,z_E_arm,dx_E_arm,dy_E_arm,delta_E_arm
	real*8 x_P_arm,y_P_arm,z_P_arm,dx_P_arm,dy_P_arm,delta_P_arm
	real*8 xtar_init_P, xtar_init_E
	real*8 xfp, yfp, dxfp, dyfp
	real*8 eloss_E_arm, eloss_P_arm, r, beta, dangles(2), dang_in(2)
	logical success
	logical ok_E_arm, ok_P_arm
	type(event):: orig, recon
	type (event_main)::	main

	real*8 beta_electron
	parameter (beta_electron = 1.)
	real*8 tmpfact
	real*8 fry	!fast raster y position.
	real*8 frx      !fast raster x position - left as looking downstream
	real*8 m2	!mass for call to mc_hms(sos). Changes if decay
	real*8 pathlen

	real*8 dx_tmp,dy_tmp

	real*8 ctheta,stheta,phad,pelec
	real*8 zhadron

	real*8 zero
	parameter (zero=0.0e0)	!double precision zero for subroutine calls

! Prepare the event for the Monte Carlo's and/or spectrometer cuts

	success = .false.
	ntup%resfac = 0.0			!resfac (see simulate.inc)
	if (correct_raster) then
	  fry = -main%target%rastery
	  frx = -main%target%rasterx   !This should make frx positive for beam left
	else
	  fry = 0.0
	  frx = 0.0
	endif

!BEAM MULTIPLE SCATTERING:  NOT YET IMPLEMENTED, JUST GETTING IT READY!

! ... multiple scattering of the beam.  Generate angles for beam deflection,
! ... and then apply those angles to the electron and proton.
! ... The yptar offset goes directly to yptar of both arms (small angle approx).
! ... xptar is multiplied by cos(theta) to get the xptar offsets.

	if (mc_smear) then
	  call target_musc(orig%Ein,beta_electron,main%target%teff(1),dang_in)
	else
	  dang_in(1)=0.0	!yp offset, goes directly to yp of both arms
	  dang_in(2)=0.0	!xp offset, mult. by cos_th for xp of both arms
	endif

	if (debug(3)) write(6,*) 'mc: using p,e arm mc =',
     >		using_P_arm_montecarlo,using_E_arm_montecarlo

!____ P arm ________________

! Go from TRUE to SPECTROMETER quantities by computing target distortions
! ... ionization loss correction (if requested)

	if (using_Eloss) then
	  if (debug(3)) write(6,*)'mc: p arm stuff0 =',
     >		orig%p%E,main%target%Eloss(3),spec%p%P
	  main%SP%p%delta = (sqrt(abs((orig%p%E-main%target%Eloss(3))**2
     >		          -Mh2))-spec%p%P) / spec%p%P*100.
	else
	  if (debug(3)) write(6,*)'mc: p arm stuff1 =',orig%p%delta
	  main%SP%p%delta = orig%p%delta
	endif

! ... multiple scattering

	if (mc_smear) then
	  beta = orig%p%p/orig%p%E
	  call target_musc(orig%p%p, beta, main%target%teff(3), dangles)
	else
	  dangles(1)=0.0
	  dangles(2)=0.0
	endif
	main%SP%p%yptar = orig%p%yptar + dangles(1) + dang_in(1)
	main%SP%p%xptar = orig%p%xptar + dangles(2) + dang_in(2)*spec%p%cos_th

! CASE 1: Using the spectrometer Monte Carlo

	if (using_P_arm_montecarlo) then

! ... change to P arm spectrometer coordinates (TRANSPORT system),

	  if (abs(cos(spec%p%phi)).gt.0.0001) then  !phi not at +/- pi/2
	    write(6,*) 'y_P_arm, z_P_arm will be incorrect if spec.p.phi <> pi/2 or 3*pi/2'
	    write(6,*) 'spec%p%phi=',spec%p%phi,'=',spec%p%phi*180/pi,'degrees'
	  endif
	  delta_P_arm = main%SP%p%delta
	  x_P_arm = -main%target%y
	  y_P_arm = -main%target%x*spec%p%cos_th - main%target%z*spec%p%sin_th*sin(spec%p%phi)
	  z_P_arm =  main%target%z*spec%p%cos_th + main%target%x*spec%p%sin_th*sin(spec%p%phi)

! %.. Apply spectrometer offset (using spectrometer coordinate system).
	  x_P_arm = x_P_arm - spec%p%offset%x
	  y_P_arm = y_P_arm - spec%p%offset%y
	  z_P_arm = z_P_arm - spec%p%offset%z
	  dx_P_arm = main%SP%p%xptar - spec%p%offset%xptar
	  dy_P_arm = main%SP%p%yptar - spec%p%offset%yptar


c	   write(*,*) 'sign_hms_part =' ,sign_hms_part
          if (using_tgt_field) then      ! do target field tracking - GAW
            if (debug(6)) then
               write(*,*) '------------------------------'
               write(*,'("frx,fry,x_E_arm =      ",3f12.5)') frx,fry,x_P_arm
           endif
            call track_from_tgt(x_P_arm,y_P_arm,z_P_arm,dx_P_arm,dy_P_arm,
     >       sign_hadron*spec%p%P*(1+delta_P_arm/100.),Mh,1,ok_P_arm)
          endif 
! GAW - end 99/11/3

C DJG need to decay the rho here before we begin transporting through the
C DJG spectrometer
c	  m2 = Mh2
c	  if(doing_rho) then
c	     call rho_decay(dx_P_arm,dy_P_arm,delta_P_arm,spec.p.P,m2,
c	1	  main.epsilon,orig.Q2)
c	  endif
C DJG moved this to the last part of generate!!!

! ........ drift this position back to z=0, the plane through the target center

	  x_P_arm = x_P_arm - z_P_arm*dx_P_arm
	  y_P_arm = y_P_arm - z_P_arm*dy_P_arm
	  z_P_arm = 0.0
	  xtar_init_P=x_P_arm

	  main%SP%p%z=y_P_arm

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
	  pathlen = 0.0
	  if (hadron_arm.eq.1) then
	    call mc_hms(spec%p%P, spec%p%theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, xtar_init_P, ok_P_arm, pathlen)
	  else if (hadron_arm.eq.2) then
	    call mc_sos(spec%p%P, spec%p%theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, fry, ok_P_arm, pathlen)
	  else if (hadron_arm.eq.3) then
	    call mc_hrsr(spec%p%P, spec%p%theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, fry, ok_P_arm, pathlen)
	  else if (hadron_arm.eq.4) then
	    call mc_hrsl(spec%p%P, spec%p%theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, fry, ok_P_arm, pathlen)
	  else if (hadron_arm.eq.5 .or. hadron_arm.eq.6) then
	    call mc_shms(spec%p%P, spec%p%theta, delta_P_arm, x_P_arm,
     >		y_P_arm, z_P_arm, dx_P_arm, dy_P_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, xtar_init_P, ok_P_arm, pathlen, hadron_arm)
	  endif


C DJG Do polarized target field stuff if needed
C DJG Note that the call to track_to_tgt is with -fry: for some reason that routine then
C DJG sets xx=-fry, and calls mc_hms_recon with xx, so it all works out. Whatever.
C DJG To summarize: fry points "down", frx points beam left (facing downstream)
C DJG For spectrometers to the right of the beamline, need to pass ctheta,stheta to track_to_tgt
C DJG For spectrometers to the left of the beamline, need to pass ctheta,-stheta to track_to_tgt

          if (using_tgt_field) then
	     phad = spec%p%P*(1.+delta_P_arm/100.0)
	     ctheta = cos(spec%p%theta)
	     stheta = sin(spec%p%theta)
	     if((hadron_arm.eq.1).or.(hadron_arm.eq.3)) then !spectrometers of the right (HMS, HRSR)
		call track_to_tgt(delta_P_arm,y_P_arm,dx_P_arm,dy_P_arm,frx,-fry,
     >   	     sign_hadron*phad,sqrt(m2),ctheta,stheta,1,ok_P_arm,hadron_arm)
	     else if((hadron_arm.eq.2).or.(hadron_arm.eq.4).or.(hadron_arm.eq.5)) then !left spects (SOS,HRSL,SHMS)
		call track_to_tgt(delta_P_arm,y_P_arm,dx_P_arm,dy_P_arm,frx,-fry,
     >  	     sign_hadron*phad,sqrt(m2),ctheta,-stheta,1,ok_P_arm,hadron_arm)
	     else
		write(6,*) 'Target field reconstruction not set up for your spectrometer'
		stop
	     endif
	  endif


	  if (.not.ok_P_arm) then
	     if (debug(3)) write(6,*)'mc: ok_P_arm =',ok_P_arm
	     return
	  endif

	  main%RECON%p%delta = delta_P_arm
	  main%RECON%p%yptar = dy_P_arm
	  main%RECON%p%xptar = dx_P_arm
	  main%RECON%p%z = y_P_arm
	  main%FP%p%x = xfp
	  main%FP%p%dx = dxfp
	  main%FP%p%y = yfp
	  main%FP%p%dy = dyfp
	  main%FP%p%path = pathlen

! CASE 2: Not using the detailed Monte Carlo, just copy the IN event to the
! OUT record

	else
	  main%RECON%p%delta = main%SP%p%delta
	  main%RECON%p%yptar = main%SP%p%yptar
	  main%RECON%p%xptar = main%SP%p%xptar
	endif

! Fill recon quantities.

	recon%p%delta = main%RECON%p%delta
	recon%p%yptar = main%RECON%p%yptar
	recon%p%xptar = main%RECON%p%xptar
	recon%p%z = main%RECON%p%z
	recon%p%P = spec%p%P*(1.+recon%p%delta/100.)
	recon%p%E = sqrt(recon%p%P**2 + Mh2)
	dx_tmp = recon%p%xptar + spec%p%offset%xptar
	dy_tmp = recon%p%yptar + spec%p%offset%yptar

	call physics_angles(spec%p%theta,spec%p%phi,dx_tmp,
     >           dy_tmp,recon%p%theta,recon%p%phi)

! ... correct for energy loss - use most probable (last flag = 4)

	if (correct_Eloss) then
	  call trip_thru_target (3, zero, recon%p%E,
     >		recon%p%theta, eloss_P_arm, r,Mh,4)
	  recon%p%E = recon%p%E + eloss_P_arm
	  recon%p%E = max(recon%p%E,sqrt(Mh2+0.000001)) !can get P~0 when calculating hadron momentum-->P<0 after eloss
	  recon%p%P = sqrt(recon%p%E**2-Mh2)
C DJG Should not correct delta for Eloss - delta is a SPECTROMETER variable
C	  recon%p%delta = (recon%p%P-spec%p%P)/spec%p%P*100.
	endif

!___ E arm ______________

! Go from TRUE to SPECTROMETER quantities by computing target distortions

! ... ionization loss correction (if requested) and Coulomb deceleration

	if (debug(3)) write(6,*)'mc: e arm stuff0 =',
     >		orig%e%E,main%target%Eloss(2),spec%e%P
	main%SP%e%delta=100* (orig%e%E - main%target%Eloss(2)
     >		- main%target%Coulomb - spec%e%P) / spec%e%P

! ... multiple scattering

	if (mc_smear) then
	  call target_musc(orig%e%p, beta_electron, main%target%teff(2), dangles)
	else
	  dangles(1)=0.0
	  dangles(2)=0.0
	endif

	main%SP%e%yptar = orig%e%yptar + dangles(1) + dang_in(1)
	main%SP%e%xptar = orig%e%xptar + dangles(2) + dang_in(2)*spec%e%cos_th

! CASE 1: Using the spectrometer Monte Carlo

	if (using_E_arm_montecarlo) then

! ... change to E arm spectrometer coordinates (TRANSPORT system),

	  if (abs(cos(spec%e%phi)).gt.0.0001) then  !phi not at +/- pi/2
	    write(6,*) 'y_E_arm, z_E_arm will be incorrect if spec.e.phi <> pi/2 or 3*pi/2'
	    write(6,*) 'spec.e.phi=',spec%e%phi,'=',spec%e%phi*180/pi,'degrees'
	  endif
	  delta_E_arm = main%SP%e%delta
	  x_E_arm = -main%target%y
	  y_E_arm = -main%target%x*spec%e%cos_th - main%target%z*spec%e%sin_th*sin(spec%e%phi)
	  z_E_arm =  main%target%z*spec%e%cos_th + main%target%x*spec%e%sin_th*sin(spec%e%phi)

! ... Apply spectrometer offset (using spectrometer coordinate system).
	  x_E_arm = x_E_arm - spec%e%offset%x
	  y_E_arm = y_E_arm - spec%e%offset%y
	  z_E_arm = z_E_arm - spec%e%offset%z
	  dx_E_arm = main%SP%e%xptar - spec%e%offset%xptar
	  dy_E_arm = main%SP%e%yptar - spec%e%offset%yptar

! GAW -project to z=0 to compare with reconstructed target positions
          if (using_tgt_field) then      ! do target field tracking - GAW

            if (debug(6)) then
               write(*,*) '------------------------------'
               write(*,'("frx,fry,x_E_arm =      ",3f12.5)') frx,fry,x_E_arm
            endif
            call track_from_tgt(x_E_arm,y_E_arm,z_E_arm,dx_E_arm,dy_E_arm,
     >                          -spec%e%P*(1+delta_E_arm/100.),Me,-1,ok_E_arm)
          endif 

! ........ drift this position back to z=0, the plane through the target center

	  x_E_arm = x_E_arm - z_E_arm*dx_E_arm
	  y_E_arm = y_E_arm - z_E_arm*dy_E_arm
	  z_E_arm = 0.0
	  xtar_init_E = x_E_arm

	  main%SP%e%z=y_E_arm

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
	  pathlen = 0.0
	  if (electron_arm.eq.1) then
	    call mc_hms(spec%e%P, spec%e%theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		me2, mc_smear, mc_smear, .false.,
     >		tmpfact, xtar_init_E, ok_E_arm, pathlen)
	  else if (electron_arm.eq.2) then
	    call mc_sos(spec%e%P, spec%e%theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		me2, mc_smear, mc_smear, .false.,
     >		tmpfact, fry, ok_E_arm, pathlen)
	  else if (electron_arm.eq.3) then
	    call mc_hrsr(spec%e%P, spec%e%theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		me2, mc_smear, mc_smear, .false.,
     >		tmpfact, fry, ok_E_arm, pathlen)
	  else if (electron_arm.eq.4) then
	    call mc_hrsl(spec%e%P, spec%e%theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		me2, mc_smear, mc_smear, .false.,
     >		tmpfact, fry, ok_E_arm, pathlen)
	  else if (electron_arm.eq.5 .or. electron_arm.eq.6) then
	    call mc_shms(spec%e%P, spec%e%theta, delta_E_arm, x_E_arm,
     >		y_E_arm, z_E_arm, dx_E_arm, dy_E_arm, xfp, dxfp, yfp, dyfp,
     >		me2, mc_smear, mc_smear, .false.,
     >		tmpfact, xtar_init_E, ok_E_arm, pathlen, electron_arm)
	  else if (electron_arm.eq.7 .or. electron_arm .eq. 8) then
             if (abs(spec%p%phi-pi/2) .eq. 10.) then
	     zhadron = -recon%p%z*(cos(spec%p%theta)/tan(spec%p%theta+recon%p%yptar)+sin(spec%p%theta)) ! recon.p.z is really ytgt
	     else
	     zhadron = recon%p%z*(cos(spec%p%theta)/tan(spec%p%theta-recon%p%yptar)+sin(spec%p%theta))
	     endif
	    call mc_calo(spec%e%p, spec%e%theta, delta_e_arm, x_e_arm,
     >		y_e_arm, z_e_arm, dx_e_arm, dy_e_arm, xfp, dxfp, yfp, dyfp,
     >		m2, mc_smear, mc_smear, doing_decay,
     >		ntup%resfac, frx, fry, ok_e_arm, pathlen, using_tgt_field,
     >          zhadron,electron_arm,drift_to_cal)
	  endif
	  ntup%resfac=ntup%resfac+tmpfact


C DJG Do polarized target field stuff if needed
C DJG Note that the call to track_to_tgt is with -fry: for some reason that routine then
C DJG sets xx=-fry, and calls mc_hms_recon with xx, so it all works out. Whatever.
C DJG To summarize: fry points "down", frx points beam left (facing downstream)
C DJG For spectrometers to the right of the beamline, need to pass ctheta,stheta to track_to_tgt
C DJG For spectrometers to the left of the beamline, need to pass ctheta,-stheta to track_to_tgt

          if (using_tgt_field) then
	     pelec = spec%e%P*(1.+delta_E_arm/100.0)
	     ctheta = cos(spec%e%theta)
	     stheta = sin(spec%e%theta)
	     if((electron_arm.eq.1).or.(electron_arm.eq.3)) then !spectrometers on the right (HMS, HRSR)
		call track_to_tgt(delta_E_arm,y_E_arm,dx_E_arm,dy_E_arm,frx,-fry,
     > 	          -1.0*pelec,sqrt(me2),ctheta,stheta,-1,ok_E_arm,electron_arm)
	     else if((electron_arm.eq.2).or.(electron_arm.eq.4).or.(electron_arm.eq.5)) then !left spects(SOS,HRSL,SHMS)
		call track_to_tgt(delta_E_arm,y_E_arm,dx_E_arm,dy_E_arm,frx,-fry,
     >	           -1.0*pelec,sqrt(me2),ctheta,-stheta,-1,ok_E_arm,electron_arm)
	     else
		write(6,*) 'Target field reconstruction not set up for your spectrometer'
		stop
	     endif
	  endif

	  if (.not.ok_E_arm) then
	    if (debug(3)) write(6,*)'mc: ok_E_arm =',ok_E_arm
	    return
	  endif

	  main%RECON%e%delta = delta_E_arm
	  main%RECON%e%yptar = dy_E_arm
	  main%RECON%e%xptar = dx_E_arm
	  main%RECON%e%z = y_E_arm
	  main%FP%e%x = xfp
	  main%FP%e%dx = dxfp
	  main%FP%e%y = yfp
	  main%FP%e%dy = dyfp
	  main%FP%e%path = pathlen

! CASE 2: Not using the detailed Monte Carlo, just copy the IN event to the
! OUT record

	else
	  main%RECON%e%delta = main%SP%e%delta
	  main%RECON%e%yptar = main%SP%e%yptar
	  main%RECON%e%xptar = main%SP%e%xptar
	endif

! Fill recon quantities%

	recon%e%delta = main%RECON%e%delta
	recon%e%yptar = main%RECON%e%yptar
	recon%e%xptar = main%RECON%e%xptar
	recon%e%z = main%RECON%e%z
	recon%e%P = spec%e%P*(1.+recon%e%delta/100.)
	recon%e%E = recon%e%P

	dx_tmp = recon%e%xptar + spec%e%offset%xptar
        dy_tmp = recon%e%yptar + spec%e%offset%yptar

	call physics_angles(spec%e%theta,spec%e%phi,dx_tmp,
     >		dy_tmp,recon%e%theta,recon%e%phi)


! ... correct for energy loss and correct for Coulomb deceleration

	if (correct_Eloss) then
	  call trip_thru_target (2, zero, recon%e%E, recon%e%theta,
     >                              eloss_E_arm, r, Me, 4)
	  recon%e%E = recon%e%E + eloss_E_arm
	endif
c	recon%e%E = recon%e%E + targ%Coulomb%ave
c Generally, we do not correct for coulomb effects in the reconstruction
	recon%e%E = recon%e%E 
	recon%e%P = recon%e%E
C DJG Should not correct delta for energy loss - delta is a SPECTROMETER
C variable!!
c	recon.e.delta = (recon.e.P-spec.e.P)/spec.e.P*100.

! Made it!
	success = .true.

	return
	end
