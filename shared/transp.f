	subroutine transp(spectr,class,decay_flag,dflag,m2,ph,zd,pathlen)
C+______________________________________________________________________________
!
! TRANSP - This subroutine transports a particle through the various
!   segments of a spectrometer. The passed variable CLASS determines
!   which transformation to use. Each transformation IS SEQUENTIAL,
!   and carries the particle from the last transformation to a particular
!   plane Z=const in the spectrometer.
!
!   NOTE: The coordinate system used here is the right handed "TRANSPORT"
!   coordinate system, in which +Z -> downstream, +X -> in median plane,
!   pointing in direction taken through bending magnet by high momentum rays,
!   +Y = transverse direction.
!
!   For an upward vertical bend spectrometer:
!
!   +X points towards the floor,
!   +Y points horizontally LEFT as one looks downstream.
!
! D. Potterveld - March 1993.
!
!   Modification History:
!
! August 1993.	(D. Potterveld) Modified to begin each transformation at
!		the pivot.
!
!  20-AUG-1993	(D. Potterveld) Modified to use COSY INFINITY forward maps.
!
!  03-MAR-1994	(D. Potterveld) Fixed bug to correctly compute cosy "A" and "B".
!
!  10-MAY-1995	(D. Potterveld) Switched to sequential transformations.
!
!  11-MAY-1995	(D. Potterveld) Added rejection of extremely small coeff.
!
!  30-OCT-1995	(D. Potterveld) Change to COSY-7 transport units.
!
!  16-SEP-1998  Check for decay of particle. dflag is true of the particle
!		has already decayed, so check for decay if dflag .eq. .false.
!		After pathlength is calculated, add extra check for decay
!		within the different path length.
!
C-______________________________________________________________________________

	implicit 	none

	include '../spectrometers.inc'
	include '../simulate.inc'

C Arguments.

	integer*4       spectr			!HMS=1, SOS=2, HRSR=3, HRSL=4
	integer*4	class
	character*(*)	file			!Used in entry point below.
	logical		decay_flag		!check for decay
	logical		dflag			!has particle decayed yet?
	real*8		m2,ph,zd, m_final		!decay variables.
	real*8		pathlen

C Parameters.

! Defined in simulate.inc
!	integer*4	max_class
!	parameter	(max_class = 18)	!max number of classes
	real*8		coeff_min
	parameter	(coeff_min = 1.0e-14)	!coeff's smaller than this = 0

C Local declarations.

	integer*4	idummy,iterm

	real*8		delta_z		!pathlength difference.
	integer*4	i,j,k,kk,chan,n_classes

	character*132	str_line
	character*132	file_name

C Cosy reconstruction matrix elements.

	integer*4	max_elements
	parameter	(max_elements = 1000)
! Defined in simulate.inc
!	integer*4	nspectr
!	parameter       (nspectr=6)
	real*8		coeff(nspectr,5,max_elements,max_class)
	integer*2	expon(nspectr,5,max_elements,max_class)
	integer*4	n_terms(nspectr,max_class)
	real*8		sum(5),ray(5),term,temp
	real*8		length(nspectr,max_class)

	integer*2	order
	integer*2	e1,e2,e3,e4		!temp. exponants
	real*8		c1,c2,c3,c4,csum	!temp. coeffs.

C decay variables
	real*8 p_spec
	real*8 z_decay
	real*8 rph,rth1,rth
	real*8 beta,gamma,dlen
	real*8 ef,pf,pxf,pyf,pzf,pxr,pyr,pzr
	real*8 bx,by,bz,er,pr


C Function definitions.

!	integer*4	last_char
	logical		locforunt
	real*8 grnd

C No amnesia between calls!!!

	save

C ================================ Executable Code =============================

C A word from the sponsor:
C  a) Using sequential matrix elements, variables with the ending 's' should be
C     used, otherwise they should end on '_transp'
C  b) the sequential matrix elements have been made with COSY 7, i.e. the units
C      to be used are cm, mrad, not m, slopes

C Check that the path length passed to code (zd) is consistent with comments
C in MEs (if they exist)

	if (length(spectr,class).gt.0.01 .and.
     >		abs(length(spectr,class)-zd).gt.0.01 ) then
	    write(6,*) 'PROBLEM WITH TRANSFORMATION #',class,' for spectrometer #',spectr
	    write(6,*) 'Central path length passed to transp =',zd
	    write(6,*) 'But the comments say that the length =',length(spectr,class)
	endif

C Check for decay.

	if (decay_flag .and. .not.dflag) then	!check for decay
	  p_spec = ph/(1.+dpps/100.)
	  beta = ph/sqrt(ph**2+m2)
	  gamma = 1./sqrt(1.-beta*beta)
	  dlen=ctau*beta*gamma   !1/e decay length (beta*c*tau*gamma)
	  if (zd.le.0) write(6,*) 'transp distance<0:  automatic decay! BAD!'
	  z_decay = -1.*dlen*log(1-grnd())

	  if (z_decay .le. zd/2) then	!decay in first half
	    dflag=.true.
	    decdist = decdist + z_decay

	    rph = grnd()*2.*pi
	    rth1 = grnd()*2.-1.
	    rth = acos(rth1)

	    pr = 0.
	    m_final = Mmu ! default
	    if(abs(sqrt(m2) - Mpi).lt.2) pr = 29.783 ! pion decay
	    if(abs(sqrt(m2) - Mk).lt.2) then ! kaons
	       if(grnd().lt.0.7) then ! decay to muon plus neutrino
		  pr = 235.5 
	       else		! decay to two pions
		  pr = sqrt(Mk**2 / 4. - Mpi**2)
		  m_final = Mk
	       endif
	    endif
	    if(pr.eq.0.) then
	     write(6,'(''error, cannot decay particle with'',
     >        '' mass='',f8.2)') sqrt(m2)
	     stop
	    endif
	    er = sqrt(m_final**2 + pr**2)
	    pxr = pr*sin(rth)*cos(rph)
	    pyr = pr*sin(rth)*sin(rph)
	    pzr = pr*cos(rth)

C Boost from pion/kaon center of mass back to lab frame.  Loren wants
C beta of new frame w.r.t. original frame, so beta is opposite of the
C initial particle's momentum.

	    bx = -beta * dxdzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    by = -beta * dydzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    bz = -beta *   1.  / sqrt(1. + dxdzs**2 + dydzs**2)
	    call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)
	    dxdzs = pxf/pzf
	    dydzs = pyf/pzf
	    dpps = 100.*(pf/p_spec-1.)
	    ph=pf
	    m2 = m_final**2     !need mass-squared for multiple scattering.
	    Mh2_final = m2      !for ntuple
	  endif

	endif				!check for decay

C Pack local copy of input coordinates.

	ray(1) = xs			!cm.	( "X" )
	ray(2) = dxdzs*1000.e0		!mrad.	( "THETA" )
	ray(3) = ys			!cm.	( "Y" )
	ray(4) = dydzs*1000.e0		!mrad.	( "PHI" )
	ray(5) = dpps			!Fractional "Delta P/P"

C Reset COSY sums.

	do i = 1,5
	   sum(i) = 0.
	enddo

C Compute COSY sums.

	k = class
	do i = 1,n_terms(spectr,k)
	  term = 1.0e0
	  do j = 1,5
	    temp = 1.0e0
	    if (expon(spectr,j,i,k).ne.0.) temp = ray(j)**expon(spectr,j,i,k)
	    term = term*temp
	  enddo
	  sum(1) = sum(1) + term*coeff(spectr,1,i,k)		! NEW "X"
	  sum(2) = sum(2) + term*coeff(spectr,2,i,k)		! NEW "A"
	  sum(3) = sum(3) + term*coeff(spectr,3,i,k)		! NEW "Y"
	  sum(4) = sum(4) + term*coeff(spectr,4,i,k)		! NEW "B"
	  sum(5) = sum(5) + term*coeff(spectr,5,i,k)		! NEW "dL"
	enddo
     
C Unpack output coordinates. Note that DPPS is unchanged by transformation.
C Pathlength correction: real pathlength=nominal-sum(5), so delta_z=-sum(5)

	xs    = sum(1)			!cm
	dxdzs = sum(2)/1000.e0		!slope (mr)
	ys    = sum(3)			!cm
	dydzs = sum(4)/1000.e0		!slope (mr)
	delta_z = -sum(5)		!deltaZ (cm)

C Check for decay in 2nd half of element, which is applied AFTER trasnporting.

	if (decay_flag .and. .not.dflag) then	!check for decay
	  if (z_decay .gt. zd+delta_z) then	!don't decay within magnet.
	    decdist = decdist + (zd+delta_z)
	  else					!apply decay at end/magnet.
	    dflag=.true.
	    decdist = decdist + z_decay

	    rph = grnd()*2.*pi
	    rth1 = grnd()*2.-1.
	    rth = acos(rth1)


	    pr = 0.
	    m_final = Mmu ! default
	    if(abs(sqrt(m2) - Mpi).lt.2) pr = 29.783 ! pion decay
	    if(abs(sqrt(m2) - Mk).lt.2) then ! kaons
	       if(grnd().lt.0.7) then ! decay to muon plus neutrino
		  pr = 235.5 
	       else		! decay to two pions
		  pr = sqrt(Mk**2 / 4. - Mpi**2)
		  m_final = Mpi
	       endif
	    endif
	    if(pr.eq.0.) then
	     write(6,'(''error, cannot decay particle with'',
     >        '' mass='',f8.2)') sqrt(m2)
	     stop
	    endif
	    er = sqrt(m_final**2 + pr**2)
	    pxr = pr*sin(rth)*cos(rph)
	    pyr = pr*sin(rth)*sin(rph)
	    pzr = pr*cos(rth)


	    bx = -beta * dxdzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    by = -beta * dydzs / sqrt(1. + dxdzs**2 + dydzs**2)
	    bz = -beta *   1.  / sqrt(1. + dxdzs**2 + dydzs**2)
	    call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)
	    dxdzs = pxf/pzf
	    dydzs = pyf/pzf
	    dpps = 100.*(pf/p_spec-1.)
	    ph=pf
	    m2 = m_final**2     !need mass-squared for multiple scattering.
	    Mh2_final = m2      !for ntuple
	  endif
	endif

C Keep track of pathlength.
	pathlen = pathlen + (zd+delta_z)

	return

C###############################################################################

C Initialization entry points.

	entry transp_init_file(file,n_classes)

C Use passed filename.

	file_name = file
	goto 100

	entry transp_init(spectr,n_classes)

C Use default filename.
C    HMS: if spectr = 1
C    SOS: if spectr = 2

	if (spectr.eq.1) then
	  file_name='hms/forward_cosy.dat'
	else if (spectr.eq.2) then
	  file_name='sos/forward_cosy.dat'
	else if (spectr.eq.3) then
	  file_name='hrsr/hrs_forward_cosy.dat'
	else if (spectr.eq.4) then
	  file_name='hrsl/hrs_forward_cosy.dat'
	else if (spectr.eq.5) then
	  file_name='shms/shms_forward.dat'
	else if (spectr.eq.6) then
c	  file_name='shms/shms_forward_cosy_LSA.dat'
	  write(6,*) 'LSA tune for SHMS no longer used!'
	  stop
	endif

C Open input file.

100     if (.not.locforunt(chan))
     >  stop 'TRANSP_INIT: No I/O channels!'
c	open (unit=chan,status='old',name=file_name,readonly)
	open (unit=chan,status='old',file=file_name)

C Strip away header.

	str_line = '!'
	do while (str_line(1:1).eq.'!')
	  read (chan,1001) str_line
!	  if (str_line(1:1).eq.'!') write(6,*) str_line(1:min(79,last_char(str_line)))
	  if (str_line(1:8).eq.'!LENGTH:') then
	    read (str_line(10:),*) length(spectr,1) !get length from comments
	    length(spectr,1)=100.*length(spectr,1) !convert to cm.
!	    write(6,*) 'kk,length=',1,length(spectr,1)
	  endif
	enddo

C Read in the transformation tables.

	n_classes = 0
	do i = 1,max_class
	  n_terms(spectr,i) = 0
	  adrift(spectr,i) = .true.
	  driftdist(spectr,i) = 0.
	enddo

C SHMS does not have the usual comments (such as drift lengths).
C For now, we have to hardwire them here.
cdg	if (spectr .eq. 5) then
cdg	  length(spectr,1)= 464.575
cdg	  length(spectr,4)= 111.15
cdg	  length(spectr,7)= 86.842
cdg	  length(spectr,12)= 308.267
cdg	endif
cdg	if (spectr .eq. 6) then
cdg	  length(spectr,1)= 232.575
cdg	  length(spectr,4)= 111.15
cdg	  length(spectr,7)=  61.519
cdg	  length(spectr,12)= 382.944
cdg	endif

	do while (.true.)
	  kk = n_classes + 1

! If too many transformations, complain!

	  if (kk.gt.max_class) stop 'TRANSP_INIT: too many transformations!'

! Add data lines to table, looking for flag line.

	  do while (str_line(1:4).ne.' ---')
	    n_terms(spectr,kk) = n_terms(spectr,kk) + 1
	    if (n_terms(spectr,kk).gt.max_elements)
     >      stop 'TRANSP_INIT: too many COSY terms!'

! Read in MEs.  Coefficients are for calculating : X, XP, Y, YP, dZ
!		Exponents determine the powers of: X, XP, Y, YP, dZ, Delta

	    read (str_line,1200)
     >		(coeff(spectr,i,n_terms(spectr,kk),kk),i=1,5),
     >		(expon(spectr,j,n_terms(spectr,kk),kk),j=1,4),idummy,
     >		 expon(spectr,5,n_terms(spectr,kk),kk)

! Ignore time-of-flight term.

	    if (idummy.ne.0) then
	      if (coeff(spectr,1,n_terms(spectr,kk),kk).ne.0.or.
     >		  coeff(spectr,2,n_terms(spectr,kk),kk).ne.0.or.
     >		  coeff(spectr,3,n_terms(spectr,kk),kk).ne.0.or.
     >		  coeff(spectr,4,n_terms(spectr,kk),kk).ne.0)
     >			stop 'TRANSP_INIT: non-zero TOF terms!'
	      n_terms(spectr,kk) = n_terms(spectr,kk) - 1
	    endif

! Check to see if element is inconsistent with a field-free drift.
! 1st order terms require: diagonal terms have coeff=1
! 			   <x|xp> and <y|yp> terms have coeff=distance(m)/10
!			   other coeff=0
! Other order terms require all coeff=0

	    if (idummy.eq.0 .and. adrift(spectr,kk)) then
	      iterm = n_terms(spectr,kk)
	      e1 = expon(spectr,1,iterm,kk)
	      e2 = expon(spectr,2,iterm,kk)
	      e3 = expon(spectr,3,iterm,kk)
	      e4 = expon(spectr,4,iterm,kk)
	      c1 = coeff(spectr,1,iterm,kk)
	      c2 = coeff(spectr,2,iterm,kk)
	      c3 = coeff(spectr,3,iterm,kk)
	      c4 = coeff(spectr,4,iterm,kk)
	      csum = abs(c1)+abs(c2)+abs(c3)+abs(c4)
	      order = e1 + e2 + e3 + e4

	      if (order.eq.1) then
	        if (e1.eq.1) then
		  if (abs(c1-1.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c2-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c3-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c4-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		else if (e2.eq.1) then
		  driftdist(spectr,kk)=1000.e0*c1     !drift distance in cm.
		  if (abs(c2-1.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c3-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c4-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		else if (e3.eq.1) then
		  if (abs(c1-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c2-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c3-1.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c4-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		else if (e4.eq.1) then
		  if (abs(c1-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(c2-0.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		  if (abs(driftdist(spectr,kk)-1000.e0*c3).gt.coeff_min) 
     >					     adrift(spectr,kk)=.false. 
		  if (abs(c4-1.e0) .gt. coeff_min) adrift(spectr,kk)=.false.
		endif
	      else	!if order.ne.1
		if (abs(csum).gt.coeff_min) adrift(spectr,kk)=.false.
	      endif
	    endif

! Fetch next line from file.

	    read (chan,1001) str_line
	  enddo

! If flag line is seen, increment transformation counter.

	  n_classes = kk

	  if (adrift(spectr,kk) .and.
     >		abs(driftdist(spectr,kk)-length(spectr,kk)).gt.0.01) then
	    write(6,*) 'PROBLEM WITH TRANSFORMATION #',kk,' for spectrometer #',spectr
	    write(6,*) 'Appears to be a pure drift, driftdist=',driftdist(spectr,kk)
	    write(6,*) 'But the comments say that the length =',length(spectr,kk)
	  endif

!!	  write(6,*) 'Spec=',spectr,'  Matrix=',kk,'  Drift=',
!!     >		adrift(spectr,kk),'  L=',driftdist(spectr,kk)

! Read lines until a non-blank, non-comment non-terminal line is found.

150       read (chan,1001,end=200) str_line
	  if (str_line(1:8).eq.'!LENGTH:') then
	    read (str_line(10:),*) length(spectr,kk+1) !get length from comments
	    length(spectr,kk+1)=100.*length(spectr,kk+1) !convert to cm.
!	    write(6,*) 'kk+1,length=',kk+1,length(spectr,kk+1)
	  endif
	  if (str_line(1:1).eq.'!'.or.str_line(1:4).eq.' ---'.or.
     >    str_line.eq.'    ') goto 150

	enddo

C Done with file.

200     close (unit=chan)

C Go home.

	return

C ============================== Format Statements =============================
1001    format(a)
1101    format(12x,6e11.3)
1102    format(10x,5f10.5)
1200    format(1x,5g14.7,1x,6i1)

	end
