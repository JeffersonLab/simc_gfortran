	subroutine dbase_read(H)

! dbase_read reads the input file and sets the run-time flags for the run.
! Note that there are four INDEPENDENT ways to run SIMC.
! 1. doing_eep:  (e,e'p) subcases:doing_hyd_elast, doing_deuterium, doing_heavy
! 2. doing_kaon: (e,e'K) subcases:doing_hydkaon, doing_deutkaon,doing_hekaon.
!		  which_kaon determines Lambda vs Sigam0 vs Sigma- production.
!		  For bound hypernuclear states, set doing_hydkaon true for ALL
!		  targets (i.e. treat as a heavy proton) but with targ.Mtar_kaon
!		  and targ.Mrec_kaon set appropriatly.
!
! 3. doing_pion: (e,e'p) subcases:doing_hydpi, doing_deutpi,doing_hepi.
!		  doing_piplus determines pi+ vs. pi- production.
! 4. doing_phsp:  Generate acceptance with radiation and cross section disabled,
!		  use doing_kaon or doing_pion to set hadron mass, then 
!		  set doing_eep,doing_kaon and doing_pion to FALSE.
!

	implicit none
	include 'radc.inc'
	include 'histograms.inc'
	include 'simulate.inc'

	real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7

	integer*4 iread, iq2, iang
	integer*4 i, j, k, ii
	integer*4 ierr, thload, thbook
	logical success
	character filename*80
	character dbase_file*60 !needs to be shorter than filename

	record /histograms/ H

! ... Register CTP vars

	call regallvars

! ... read the dbase file.

	dbase_file=' '
	extra_dbase_file=' '
	write(6,*) 'Enter the input filename (assumed to be in infiles directory)'
	read(5,'(a)') dbase_file

 7	write(filename,'(a)') 'infiles/'//dbase_file
	write(6,'(a10,a69)')'filename=',filename
	i=index(filename,'.')
	j=index(filename,'/')
	if(i.eq.0) then		!add .inp if not included in filename
	  i=index(filename,' ')
	  if(i+2.le.len(filename)) write(filename(i:),'(''.inp'')')
	endif
	if (i.gt.1) base=filename(j+1:i-1)

! ... load and book input file

	ierr = thload(filename)
	if (ierr.ne.0) stop 'Loading problem!  Not going to try again...wouldnt be prudent.'
	ierr = thbook()
	if (ierr.ne.0) stop ' Booking problem!  Not going to try again...wouldnt be prudent'

! ... read the secondary dbase file.

	if (extra_dbase_file.ne.' ') then	!new filename from dbase file

	  write(filename,'(a)') 'infiles/'//extra_dbase_file
	  write(6,'(a10,a69)')'filename=',filename
	  i=index(filename,'.')
	  j=index(filename,'/')
	  if(i.eq.0) then		!add .inp if not included in filename
	    i=index(filename,' ')
	    i f(i+2.le.len(filename)) write(filename(i:),'(''.inp'')')
	  endif

	  ierr = thload(filename)
	  if (ierr.ne.0) stop ' Loading problem!  Not going to try again...wouldnt be prudent.'
	  ierr = thbook()
	  if (ierr.ne.0) stop ' Booking problem!  Not going to try again...wouldnt be prudent'
	endif	!extra dbase input file

! ... dbase field experiment.
	if (doing_pion) then
	  Mh=Mpi
	  if (nint(targ.A).eq.1 .and. .not.doing_piplus) 
     >		stop 'Pi- production from Hydrogen not allowed!'
	  doing_hydpi = (nint(targ.A).eq.1)
	  doing_deutpi = (nint(targ.A).eq.2)
	  doing_hepi = (nint(targ.A).eq.3)
	  doing_eep = .false.
	else if (doing_kaon) then
	  Mh=Mk
	  doing_hydkaon = (nint(targ.A).eq.1)
	  doing_deutkaon = (nint(targ.A).eq.2)
	  doing_hekaon = (nint(targ.A).ge.3)
	  doing_eep = .false.

* Treat bound hypernuclei as production from heavy proton (add 10 to which_kaon)
* Treat deuteron-lambda/sigma state as production from heavy deut.(which_kaon+20)
	  if (which_kaon.ge.10 .and. which_kaon.lt.20) then
	    doing_hydkaon = .true.
	    doing_deutkaon = .false.
	    doing_hekaon = .false.
	  else if (which_kaon.ge.20) then
	    doing_hydkaon = .false.
	    doing_deutkaon = .true.
	    doing_hekaon = .false.
	  endif
	else
	  Mh=Mp		!default to proton in hadron arm
	  doing_eep = .true.
	endif
	Mh2 = Mh*Mh
	Mh2_final = Mh2

	if (doing_phsp) then
	  rad_flag=0
	  doing_eep=.false.	!need to set on of these to determine Mh,
	  doing_pion=.false.	!but doing_phsp is independent of others.
	  doing_kaon=.false.
	endif

! ... dbase field kinematics_main

	dEbeam=Ebeam*dEbeam/100.
	spec.e.theta=abs(spec.e.theta)/degrad
	spec.e.cos_th=cos(spec.e.theta)
	spec.e.sin_th=sin(spec.e.theta)
	spec.p.theta=abs(spec.p.theta)/degrad
	spec.p.cos_th=cos(spec.p.theta)
	spec.p.sin_th=sin(spec.p.theta)
	if (hms_e_flag) then		!electron in hms. sos=pi/2, hms=3*pi/2
	  spec.e.phi = 3*pi/2.
	  spec.p.phi = pi/2.
	else
	  spec.e.phi = pi/2.
	  spec.p.phi = 3*pi/2.
	endif

! ... dbase field target

	targ.M = targ.mass_amu*amu
	targ.Mrec = targ.mrec_amu*amu
	targ.N = targ.A - targ.Z

	if (doing_pion) then
	  if (doing_piplus) then	! target/recoil masses for pion prod.
	    targ.Mtar_pion = Mp
	    targ.Mrec_pion = Mn
	  else
	    targ.Mtar_pion = Mn
	    targ.Mrec_pion = Mp
	  endif
	else if (doing_kaon) then
	  if (which_kaon .eq. 0) then	! p(e,eK+)Lambda0,p,n
	    targ.Mtar_kaon = Mp
	    targ.Mrec_kaon = Mlambda
	  else if (which_kaon .eq. 1) then	! p(e,eK+)Sigma0,p,n
	    targ.Mtar_kaon = Mp
	    targ.Mrec_kaon = Msigma0
	  else if (which_kaon .eq. 2) then	! n(e,eK+)Sigma,p,p
	    targ.Mtar_kaon = Mn
	    targ.Mrec_kaon = Msigma_minus
*JRA: The following give a hypernucleus mass that is bound by zero MeV
* WITH RESPECT TO free nucleons.  Want it to be with respect to Y+(A-1),
* but this will probably require more input (or else hardwired A-dependent
* parameters.
	  else if (which_kaon .eq. 10) then	!A(e,eK+)A_Lambda0 (bound hypernucleon)
	    targ.Mtar_kaon = targ.M
	    targ.Mrec_kaon = targ.Z*Mp + targ.N*Mn + Mlambda - Mp
	  else if (which_kaon .eq. 11) then	!A(e,eK+)A_Sigma0 (bound hypernucleon)
	    targ.Mtar_kaon = targ.M
	    targ.Mrec_kaon = targ.Z*Mp + targ.N*Mn + Msigma0 - Mp
	  else if (which_kaon .eq. 12) then	!A(e,eK+)A_Sigma- (bound hypernucleon)
	    targ.Mtar_kaon = targ.M
	    targ.Mrec_kaon = targ.Z*Mp + targ.N*Mn + Msigma_minus - Mn
	  else if (which_kaon .eq. 20) then	! A(e,eK+)Lambda0,D
	    targ.Mtar_kaon = Mp
	    targ.Mrec_kaon = Mlambda
	  else if (which_kaon .eq. 21) then	! A(e,eK+)Sigma0,D
	    targ.Mtar_kaon = Mp
	    targ.Mrec_kaon = Msigma0
	  else if (which_kaon .eq. 22) then	! A(e,eK+)Sigma-,p,p
	    targ.Mtar_kaon = Mn
	    targ.Mrec_kaon = Msigma_minus
	    stop 'Cant do which_kaon=22, no 2 body final state'
	  else
	    stop 'Bad value for which_kaon'
	  endif
	  if (which_kaon .ge. 10) then	!1 and 2 body final states have Mrec=0
	    if (targ.Mrec.ne.0) then
	      targ.Mrec = 0.
	      write(6,*) 'Forcing targ.Mrec to zero !!!!!'
	    endif
	  endif
	endif
	targ.thick=targ.thick/1000.
	targ.length=targ.thick/targ.rho
	targ.angle=targ.angle/degrad
        targ.Bangle = targ.Bangle/degrad
	if(targ.Z.lt.2.4) then		!liquid target
	  if (abs(targ.angle) .gt. 0.0001) then
	    targ.angle = 0.0
	    write(6,*) 'Forcing target angle to zero for cryotarget.'
	  endif
	  if (targ.can.ne.1 .and. targ.can.ne.2) stop 'bad targ.can value'
	endif
	if(sin(targ.angle) .gt. 0.85) then
	  write(6,*) 'BAD targ.angle (0 is perp. to beam, +ve is rotated towards SOS)'
	  write(6,*) 'If the target really is >60 degrees from normal to the beam, then'
	  write(6,*) 'comment out this error checking code.'
	  stop
	endif

! ... dbase field spect_offset

	spec.e.offset.xptar=spec.e.offset.xptar/1000.
	spec.e.offset.yptar=spec.e.offset.yptar/1000.
	spec.p.offset.xptar=spec.p.offset.xptar/1000.
	spec.p.offset.yptar=spec.p.offset.yptar/1000.

! ... dbase fields e_arm_accep, p_arm_accep.  Convert angles to radians.

	if (SPedge.e.delta.min.le.-100.0) SPedge.e.delta.min=-99.99
	if (SPedge.p.delta.min.le.-100.0) SPedge.p.delta.min=-99.99
	SPedge.e.yptar.min = SPedge.e.yptar.min/1000.
	SPedge.e.yptar.max = SPedge.e.yptar.max/1000.
	SPedge.e.xptar.min = SPedge.e.xptar.min/1000.
	SPedge.e.xptar.max = SPedge.e.xptar.max/1000.
	SPedge.p.yptar.min = SPedge.p.yptar.min/1000.
	SPedge.p.yptar.max = SPedge.p.yptar.max/1000.
	SPedge.p.xptar.min = SPedge.p.xptar.min/1000.
	SPedge.p.xptar.max = SPedge.p.xptar.max/1000.

! ... dbase field simulate

! ... Set flags for e,e',hadron radiation tails.

	doing_tail(1) = one_tail.eq.0 .or. one_tail.eq.1 .or.
     >		one_tail.eq.-2 .or. one_tail.eq.-3
	doing_tail(2) = one_tail.eq.0 .or. one_tail.eq.2 .or.
     >		one_tail.eq.-3 .or. one_tail.eq.-1
	doing_tail(3) = one_tail.eq.0 .or. one_tail.eq.3 .or.
     >		one_tail.eq.-1 .or. one_tail.eq.-2

	if (.not.doing_tail(1) .or. .not.doing_tail(2)) then
	  write(6,*) 'WARNING: Shutting off tails 1 or 2 not well defined'
	endif

	if (abs(one_tail).gt.3 .and. using_rad) stop
     >	  'Moron! one_tail>3 turns radiation off, but using_rad wants it on.'

	if (.not.using_rad) then
	  do i = 1, 3
	    doing_tail(i) = .false.
	  enddo
	endif

	if (Egamma_gen_max.gt.0.01) then !use hardwired limits if gen_max > 0
	  hardwired_rad = .true.
	else
	  hardwired_rad = .false.
	endif

! ... change angles to radians
	gen.e.yptar.min=gen.e.yptar.min/1000.
	gen.e.yptar.max=gen.e.yptar.max/1000.
	gen.p.yptar.min=gen.p.yptar.min/1000.
	gen.p.yptar.max=gen.p.yptar.max/1000.
	gen.e.xptar.min=gen.e.xptar.min/1000.
	gen.e.xptar.max=gen.e.xptar.max/1000.
	gen.p.xptar.min=gen.p.xptar.min/1000.
	gen.p.xptar.max=gen.p.xptar.max/1000.

! ... set some flags
	using_E_arm_montecarlo = spect_mode.ne.1 .and. spect_mode.ne.-1
	using_P_arm_montecarlo = spect_mode.ne.1 .and. spect_mode.ne.-2

	if (.not.doing_eep .or. (cuts.Em.min.eq.cuts.Em.max)) then
	  cuts.Em.min = -1.d6
	  cuts.Em.max =  1.d6
	endif

	if (abs(deForest_flag).gt.1) stop 'Idiot! check setting of deForest_flag'

	if (correct_Eloss .and. .not.using_Eloss) then
	  write(6,*) 'using_Eloss is false, SO WE ARE FORCING CORRECT_ELOSS TO BE FALSE!!!'
	  correct_Eloss = .false.
	endif

! ... target

	call target_init(using_Eloss, using_Coulomb, spec.e.theta, spec.p.theta,
     >		spec.p.P, Mh2, Ebeam, spec.e.P)

! ... For Hydrogen elastic, we never use phsp division, and the
! ... Target recoil goes to zero (all other targ.Mrec were set via Arec).

	doing_hyd_elast=doing_eep .and. nint(targ.A).eq.1
	if (doing_hyd_elast) targ.Mrec = 0.0

	doing_deuterium = doing_eep .and. nint(targ.A).eq.2
	doing_heavy = doing_eep .and. nint(targ.A).ge.3

! ... Read in the theory file (for A(e,e'p))
	if (doing_deuterium .or. doing_heavy) then
	  call theory_init(success)
	  if (.not.success) stop 'THEORY_INIT failed!'
	endif

! ... Read in (and normalize) momentum distribution.  Normalizing to one
! ... may not be correct if there is strength beyond p=p_max (MeV/c)

	if(doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon) then
	  if(doing_deutpi .or. doing_deutkaon) open(1,file='deut.dat',status='old',form='formatted')
	  if(doing_hepi .or. doing_hekaon) then
	    if (nint(targ.A).eq.3) then
	      open(1,file='he3.dat',status='old',form='formatted')
	    else if (nint(targ.A).eq.4) then
	      open(1,file='he4.dat',status='old',form='formatted')
	    endif
	  endif
	  do ii=1,2000
	    read(1,*,end=999) pval(ii),mprob(ii)
	    nump=ii
	  enddo
999	  continue
	  do ii=1,nump
	    mprob(ii)=mprob(ii)/mprob(nump)
	  enddo
	  close(1)
	endif

	if (doing_kaon) then		!read in model xsec info.
	  open(1,file='weightc.dat',status='old',form='formatted')
	  do i=1,50
	    do j=1,20
	      read(1,*) weightc(j,i)
	    enddo
	  enddo
	  close(1)

	  open(1,file='weightd.dat',status='old',form='formatted')
	  do i=1,30
	    do j=1,40
	      do k=1,8
	        read(1,*) weightd(k,j,i)
	      enddo
	    enddo
	  enddo
	  close(1)

	  open(3,file='saghai_proton.dat',status='old',form='formatted')
	  do iread=1,10
	    do iq2=1,11
	      read(3,*) dum1,dum2
	      do iang=1,19
		read(3,4005) dum3,dum4,dum5,dum6,dum7
		read(3,4005) zrff1(iread,iq2,iang),ziff1(iread,iq2,iang),
     1			zrff2(iread,iq2,iang),ziff2(iread,iq2,iang),
     2			zrff3(iread,iq2,iang),ziff3(iread,iq2,iang)
		read(3,4005) zrff4(iread,iq2,iang),ziff4(iread,iq2,iang),
     1			zrff5(iread,iq2,iang),ziff5(iread,iq2,iang),
     2			zrff6(iread,iq2,iang),ziff6(iread,iq2,iang)
	      enddo
	    enddo
	  enddo
4005	  format(6e12.4)
	  close(3)

	  open(3,file='saghai_sigma0.dat',status='old',form='formatted')
	  do iread=1,20
	    do iq2=1,10
	      do iang=1,19
		read(3,*) dum1,dum2
		read(3,4005) dum3,dum4,dum5,dum6,dum7
		read(3,4005) zsrff1(iread,iq2,iang),zsiff1(iread,iq2,iang),
     1			zsrff2(iread,iq2,iang),zsiff2(iread,iq2,iang),
     2			zsrff3(iread,iq2,iang),zsiff3(iread,iq2,iang)
		read(3,4005) zsrff4(iread,iq2,iang),zsiff4(iread,iq2,iang),
     1			zsrff5(iread,iq2,iang),zsiff5(iread,iq2,iang),
     2			zsrff6(iread,iq2,iang),zsiff6(iread,iq2,iang)
	      enddo
	    enddo
	  enddo
	  close(3)

	endif		!kaon input files.



! ... initialize limits on generation, cuts, edges

	call limits_init(H)

! ... some announcements

	if (doing_eep) then
	  if (doing_hyd_elast) then
	    write(6,*) ' ****--------  H(e,e''p)  --------****'
	  else if (doing_deuterium) then
	    write(6,*) ' ****--------  D(e,e''p)  --------****'
	  else if (doing_heavy) then
	    write(6,*) ' ****--------  A(e,e''p)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''p) we''re doing!!!'
	  endif
	else if (doing_pion) then
	  if (doing_hydpi) then
	    write(6,*) ' ****--------  H(e,e''pi)  --------****'
	  else if (doing_deutpi) then
	    write(6,*) ' ****--------  D(e,e''pi)  --------****'
	  else if (doing_hepi) then
	    write(6,*) ' ****-------- He(e,e''pi)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''pi) we''re doing!!!'
	  endif
	else if (doing_kaon) then
	  if (doing_hydkaon) then
	    if (targ.A .eq. 1) then
	      write(6,*) ' ****--------  H(e,e''K)  --------****'
	    else if (targ.A .eq. 2) then
	      write(6,*) ' So...confused...doing_hydpi...but...A...=...2..... '
	      stop
	    else if (targ.A .ge.3) then
	      write(6,*) ' ****--------  A(e,e''K)  --------****'
	    endif
	  else if (doing_deutkaon) then
	    write(6,*) ' ****--------  D(e,e''K)  --------****'
	  else if (doing_hekaon) then
	    write(6,*) ' ****--------  A(e,e''K)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	  if (which_kaon.eq.0) then
	    write(6,*) ' ****---- producing a LAMBDA ----****'
	  else if (which_kaon.eq.1) then
	    write(6,*) ' ****---- producing a SIGMA0 ----****'
	  else if (which_kaon.eq.2) then
	    write(6,*) ' ****---- producing a SIGMA- ----****'
	  else if (which_kaon.eq.10) then
	    write(6,*) ' ****---- WITH BOUND LAMBDA  ----****'
	  else if (which_kaon.eq.11) then
	    write(6,*) ' ****---- WITH BOUND SIGMA0  ----****'
	  else if (which_kaon.eq.12) then
	    write(6,*) ' ****---- WITH BOUND SIGMA-  ----****'
	  else if (which_kaon.eq.20) then
	    write(6,*) ' ****----  WITH LAMBDA-DEUT  ----****'
	  else if (which_kaon.eq.21) then
	    write(6,*) ' ****----  WITH SIGMA0-DEUT  ----****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	else if (doing_phsp) then
	  write(6,*) ' ****--------  PHASE SPACE - NO physics, NO radiation  --------****'
	else
	  stop 'I don''t have ANY idea what we''re doing!!!'
	endif

	if (doing_kaon .and. .not. doing_decay) write(6,*) 'NOTE: not doing decay, so a decay weight will be applied to SIGCC'
	if (.not.using_rad) write(6,*) 'NOTE: Will NOT be applying radiative corrections'
	if (.not.using_E_arm_montecarlo) write(6,*) 'NOTE: Will NOT be running events through the E arm Monte Carlo'
	if (.not.using_P_arm_montecarlo) write(6,*) 'NOTE: Will NOT be running events through the P arm Monte Carlo'
	if (.not.using_Eloss) write(6,*) 'NOTE: Will NOT be calculating energy loss in the target'
	if (.not.using_Coulomb) write(6,*) 'NOTE: Will NOT be calculating Coulomb correction (default for Hydrogen target)'
	if (.not.correct_Eloss) write(6,*) 'NOTE: Will NOT correct reconstructed data for energy loss'
	if (.not.correct_raster) write(6,*) 'NOTE: Will NOT use raster terms in reconstruction'
	if (using_tgt_field) write(6,'(" NOTE: Will use target field.  File is: ",a60)') tgt_field_file
	write(6,*) ' '

	return
	end

!---------------------------------------------------------------

	subroutine regallvars()

	implicit none
	include 'simulate.inc'
	include 'radc.inc'

	integer*4 ierr
	integer*4 regparmint, regparmdouble
	integer*4 regparmintarray, regparmdoublearray
	integer*4 regparmstring
!	integer*4 regparmreal

*	RESTSW

	if (debug(2).eq.1) write(6,*)'regallvars: entering'
	ierr = regparmint('mc_smear',mc_smear,0)
	ierr = regparmint('hms_e_flag',hms_e_flag,0)
	ierr = regparmstring('extra_dbase_file',extra_dbase_file,0)

*	EXPERIMENT

	ierr = regparmdouble('Ebeam',Ebeam,0)
	ierr = regparmdouble('dEbeam',dEbeam,0)
	ierr = regparmdouble('EXPER.charge',EXPER.charge,0)
	ierr = regparmint('doing_kaon',doing_kaon,0)
	ierr = regparmint('which_kaon',which_kaon,0)
	ierr = regparmint('doing_pion',doing_pion,0)
	ierr = regparmint('doing_piplus',doing_piplus,0)
	ierr = regparmint('doing_decay',doing_decay,0)
	ierr = regparmdouble('ctau',ctau,0)

*	DEBUG

	ierr = regparmintarray('debug',debug,6,0)

*	TARGET

	ierr = regparmdouble('targ.A',targ.A,0)
	ierr = regparmdouble('targ.Z',targ.Z,0)
	ierr = regparmdouble('targ.mass_amu',targ.mass_amu,0)
	ierr = regparmdouble('targ.mrec_amu',targ.mrec_amu,0)
	ierr = regparmdouble('targ.rho',targ.rho,0)
	ierr = regparmdouble('targ.thick',targ.thick,0)
	ierr = regparmdouble('targ.xoffset',targ.xoffset,0)
	ierr = regparmdouble('targ.yoffset',targ.yoffset,0)
	ierr = regparmint('targ.fr_pattern',targ.fr_pattern,0)
	ierr = regparmdouble('targ.fr1',targ.fr1,0)
	ierr = regparmdouble('targ.fr2',targ.fr2,0)
	ierr = regparmdouble('targ.zoffset',targ.zoffset,0)
	ierr = regparmdouble('targ.angle',targ.angle,0)
	ierr = regparmdouble('targ.abundancy',targ.abundancy,0)
	ierr = regparmint('targ.can',targ.can,0)
	ierr = regparmdouble('targ.Bangle',targ.Bangle,0)

*	E_ARM_MAIN

	ierr = regparmdouble('spec.e.P',spec.e.P,0)
	ierr = regparmdouble('spec.e.theta',spec.e.theta,0)

*	P_ARM_MAIN

	ierr = regparmdouble('spec.p.P',spec.p.P,0)
	ierr = regparmdouble('spec.p.theta',spec.p.theta,0)

*	E_ARM_OFFSET

	ierr = regparmdouble('gen.ywid',gen.ywid,0)
	ierr = regparmdouble('gen.xwid',gen.xwid,0)
	ierr = regparmdouble('spec.e.offset.x',spec.e.offset.x,0)
	ierr = regparmdouble('spec.e.offset.y',spec.e.offset.y,0)
	ierr = regparmdouble('spec.e.offset.z',spec.e.offset.z,0)
	ierr = regparmdouble('spec.e.offset.xptar',spec.e.offset.xptar,0)
	ierr = regparmdouble('spec.e.offset.yptar',spec.e.offset.yptar,0)

*	P_ARM_OFFSET

	ierr = regparmdouble('spec.p.offset.x',spec.p.offset.x,0)
	ierr = regparmdouble('spec.p.offset.y',spec.p.offset.y,0)
	ierr = regparmdouble('spec.p.offset.z',spec.p.offset.z,0)
	ierr = regparmdouble('spec.p.offset.xptar',spec.p.offset.xptar,0)
	ierr = regparmdouble('spec.p.offset.yptar',spec.p.offset.yptar,0)

*	MISC2INT

	ierr = regparmint('use_expon',use_expon,0)
	ierr = regparmint('one_tail',one_tail,0)
	ierr = regparmint('intcor_mode',intcor_mode,0)

*	SIMULATE

	ierr = regparmint('ngen',ngen,0)
	ierr = regparmint('using_rad',using_rad,0)
	ierr = regparmint('spect_mode',spect_mode,0)
	ierr = regparmint('doing_phsp',doing_phsp,0)
	ierr = regparmdouble('cuts.Em.min',cuts.Em.min,0)
	ierr = regparmdouble('cuts.Em.max',cuts.Em.max,0)
	ierr = regparmdouble('cuts.Pm.max',cuts.Pm.max,0)
	ierr = regparmint('using_Eloss',using_Eloss,0)
	ierr = regparmint('correct_Eloss',correct_eloss,0)
	ierr = regparmint('correct_raster',correct_raster,0)
	ierr = regparmint('deForest_flag',deForest_flag,0)
	ierr = regparmint('rad_flag',rad_flag,0)
	ierr = regparmint('extrad_flag',extrad_flag,0)
	ierr = regparmdoublearray('lambda',lambda,3,0)
	ierr = regparmint('using_cit_generation',using_cit_generation,0)
	ierr = regparmint('Nntu',Nntu,0)
	ierr = regparmint('using_Coulomb',using_Coulomb,0)
	ierr = regparmdouble('dE_edge_test',dE_edge_test,0)
	ierr = regparmint('use_offshell_rad',use_offshell_rad,0)
	ierr = regparmdouble('Egamma_gen_max',Egamma_gen_max,0)
	ierr = regparmint('using_tgt_field',using_tgt_field,0)
	ierr = regparmstring('tgt_field_file',tgt_field_file,0)
        ierr = regparmdouble('drift_to_ndet',drift_to_ndet,0)

*	E_ARM_ACCEPT

	ierr = regparmdouble('SPedge.e.delta.min',SPedge.e.delta.min,0)
	ierr = regparmdouble('SPedge.e.delta.max',SPedge.e.delta.max,0)
	ierr = regparmdouble('SPedge.e.yptar.min',SPedge.e.yptar.min,0)
	ierr = regparmdouble('SPedge.e.yptar.max',SPedge.e.yptar.max,0)
	ierr = regparmdouble('SPedge.e.xptar.min',SPedge.e.xptar.min,0)
	ierr = regparmdouble('SPedge.e.xptar.max',SPedge.e.xptar.max,0)

*	P_ARM_ACCEPT

	ierr = regparmdouble('SPedge.p.delta.min',SPedge.p.delta.min,0)
	ierr = regparmdouble('SPedge.p.delta.max',SPedge.p.delta.max,0)
	ierr = regparmdouble('SPedge.p.yptar.min',SPedge.p.yptar.min,0)
	ierr = regparmdouble('SPedge.p.yptar.max',SPedge.p.yptar.max,0)
	ierr = regparmdouble('SPedge.p.xptar.min',SPedge.p.xptar.min,0)
	ierr = regparmdouble('SPedge.p.xptar.max',SPedge.p.xptar.max,0)

	if (debug(2).eq.1) write(6,*)'regallvars: ending'
	return
	end
