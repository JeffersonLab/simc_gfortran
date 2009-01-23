	subroutine transp_nonseq(spectr,class,decay_flag,dflag,m2,ph,zd)
C+______________________________________________________________________________
!
! TRANSP - This subroutine transports a particle through the various
!   segments of a spectrometer. The passed variable CLASS determines
!   which transformation to use. Each transformation starts at the pivot,
!   and carries the particle to a particular plane Z=const in the spectrometer.
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
! May 2000 - Modified to suit Hall A HRSr spectrometer
!
! January 2001 - units are converted in the following ways
!	cm --> m
!	mrad --> rad
!	% --> fractional variation
!      at the end of the code, the units are converted back to their original form
!
!   Modification History:
!
C-______________________________________________________________________________

	implicit 	none

	include '../track.inc'
	include '../simulate.inc'

C Arguments.

	logical decay_flag,dflag
	real*8 m2,ph,zd

	integer*4	spectr			! HRSR = 3,  HRSL = 4
	integer*4	class
	character*(*)	file			!Used in entry point below.

C Parameters.

!	integer*4	max_class
!	parameter	(max_class = 15)	!maximum number of classes.

C Local declarations.

	real*8		dummy
	integer*4	idummy

	integer*4	i,j,k,kk,chan,n_classes

	character*132	str_line
	character*132	file_name

C Cosy reconstruction matrix elements.

	integer*4	max_elements
	parameter	(max_elements = 1000)
!	integer*4	nspectr
!	parameter	(nspectr=4)
	real*8		coeff(nspectr,5,max_elements,max_class)
	integer*2	expon(nspectr,5,max_elements,max_class)
	integer*4	n_terms(nspectr,max_class),max_order
	real*8		sum(4),ray(5),term,temp

C Function definitions.

	integer*4	last_char
	logical*4	locforunt

C No amnesia between calls!!!

	save

C ================================ Executable Code =============================


C Pack local copy of input coordinates.

C Using output from COSY in TRANSPORT coordinates.  HRS uses COSY-7 output
C (cm,mrad) and uses non-sequential transformations (always start at pivot)


!!! TRANSPORT COORDINATE VERSION:

	  ray(1) = x_transp     	!cm --> m    ( "X" )
	  ray(2) = dxdz_transp*1000.   	!rad  ( "THETA" )
	  ray(3) = y_transp     	!cm. --> m.   ( "Y" )
	  ray(4) = dydz_transp*1000.  	!rad.  ( "PHI" )
	  ray(5) = dpps 		!% --> frac. variation

C Reset COSY sums.

	do i = 1,4
	   sum(i) = 0.
	enddo

C Compute COSY sums.

*	if (class.eq.11) write(6,'(a,5f12.6)') 'Init= ',ray(1),ray(2),ray(3),ray(4),ray(5)

	k = class
	do i = 1,n_terms(spectr,k)
	   term = 1.0
	   do j = 1,5
	      temp = 1.0
	      if (expon(spectr,j,i,k).ne.0.) temp = ray(j)**expon(spectr,j,i,k)
	      term = term*temp
	   enddo

*	if ( (expon(spectr,1,i,k)+expon(spectr,2,i,k)+expon(spectr,3,i,k)+
*     >	      expon(spectr,4,i,k)+expon(spectr,5,i,k)).eq.1 .and. class.eq.11) then
*	write(6,'(4f12.6,1x,6i1)') term*coeff(spectr,1,i,k),term*coeff(spectr,2,i,k),
*     >		 	   term*coeff(spectr,3,i,k),term*coeff(spectr,4,i,k),
*     >		(expon(spectr,idummy,i,k),idummy=1,5)
*
*	endif

	   sum(1) = sum(1) + term*coeff(spectr,1,i,k)              !new "X"
	   sum(2) = sum(2) + term*coeff(spectr,2,i,k)              !new "A"
	   sum(3) = sum(3) + term*coeff(spectr,3,i,k)              !new "Y"
	   sum(4) = sum(4) + term*coeff(spectr,4,i,k)              !new "B"
	enddo

C Unpack output coordinates. Note that DPPS is unchanged by transformation.

!!! TRANSPORT COORDINATE VERSION:

	if (spectr.eq.3 .or. spectr.eq.4) then
	  xs    = sum(1)			!cm
	  dxdzs = sum(2)/1000			!mrad --> rad  
	  ys    = sum(3)			!cm 
	  dydzs = sum(4)/1000			!mrad --> rad 
	  zs = 0.				!Z info is lost (but not needed)
	else
	  stop 'transp_nonseq not meant for use with HMS/SOS'
        endif


*	if (class.eq.11) write(6,'(a,5f12.6)') 'Final=',xs,dxdzs*1000,ys,dydzs*1000,dpps

	return

C###############################################################################

C Initialization entry points.

	entry transp_nonseq_init_file(file,n_classes)

C Use passed filename.

	file_name = file
	goto 100

	entry transp_nonseq_init(spectr,n_classes)

C Use default filename.

        if (spectr.eq.1) then
          file_name='hms/forward_cosy.dat'
        else if (spectr.eq.2) then
          file_name='sos/forward_cosy.dat'
        else if (spectr.eq.3) then
          file_name='hrsr/hrs_forward_cosy.dat'
        else if (spectr.eq.4) then
          file_name='hrsl/hrs_forward_cosy.dat'
        endif

C Open input file.

100	if (.not.locforunt(chan))
     >	stop 'TRANSP_INIT: No I/O channels!'
	open (unit=chan,status='old',name=file_name,readonly)
C Strip away header.

	str_line = '!'
	do while (str_line(1:1).eq.'!')
	   read (chan,1001) str_line
	   if (str_line(1:1).eq.'!') type *,str_line(1:last_char(str_line))
	enddo

C Read in the transformation tables.

	n_classes = 0
	do i = 1,max_class
	   n_terms(spectr,i) = 0
	enddo

	do while (.true.)
	   kk = n_classes + 1

! If too many transformations, complain!

	   if (kk.gt.max_class) stop 'TRANSP_INIT: too many transformations!'

! Add data lines to table, looking for flag line.

	   do while (str_line(1:4).ne.' ---')
	      n_terms(spectr,kk) = n_terms(spectr,kk) + 1
	      if (n_terms(spectr,kk).gt.max_elements)
     >	      stop 'TRANSP_INIT: too many COSY terms!'

	      read (str_line,1200)
     >		(coeff(spectr,i,n_terms(spectr,kk),kk),i=1,4),dummy,
     >		(expon(spectr,j,n_terms(spectr,kk),kk),j=1,4),idummy,
     >		expon(spectr,5,n_terms(spectr,kk),kk)

! Ignore time-of-flight term.

	      if (idummy.ne.0) then
		if (coeff(spectr,1,n_terms(spectr,kk),kk).ne.0 .or.
     >		    coeff(spectr,2,n_terms(spectr,kk),kk).ne.0 .or.
     >		    coeff(spectr,3,n_terms(spectr,kk),kk).ne.0 .or.
     >		    coeff(spectr,4,n_terms(spectr,kk),kk).ne.0)
     >		stop 'TRANSP_INIT: non-zero TOF terms!'
		n_terms(spectr,kk) = n_terms(spectr,kk) - 1
	      elseif (coeff(spectr,1,n_terms(spectr,kk),kk).eq.0 .and.
     >                coeff(spectr,2,n_terms(spectr,kk),kk).eq.0 .and.
     >                coeff(spectr,3,n_terms(spectr,kk),kk).eq.0 .and.
     >                coeff(spectr,4,n_terms(spectr,kk),kk).eq.0) then
		n_terms(spectr,kk) = n_terms(spectr,kk) - 1

! For ordinary terms, keep track of the maximum order encountered.

	      else

		max_order = max(max_order,
     >                  expon(spectr,1,n_terms(spectr,kk),kk) +
     >                  expon(spectr,2,n_terms(spectr,kk),kk) +
     >                  expon(spectr,3,n_terms(spectr,kk),kk) +
     >                  expon(spectr,4,n_terms(spectr,kk),kk) +
     >                  expon(spectr,5,n_terms(spectr,kk),kk))
	      endif

! Fetch next line from file.

	      read (chan,1001) str_line
	   enddo

! If flag line is seen, increment transformation counter.

	   n_classes = kk
!	   type *,'TRANS, ORDER, TERMS =',kk,max_order,n_terms(spectr,kk)

! Read lines until a non-blank, non-comment non-terminal line is found.

150	   read (chan,1001,end=200) str_line
	   if (str_line(1:1).eq.'!'.or.str_line(1:4).eq.' ---'.or.
     >	   str_line.eq.'    ') goto 150
	enddo

C Done with file.

200	close (unit=chan)


C Go home.

	return

C ============================== Format Statements =============================

1001	format(a)
1101	format(12x,6e11.3)
1102	format(10x,5f10.5)
1200	format(1x,5g14.7,1x,6i1)

	end
