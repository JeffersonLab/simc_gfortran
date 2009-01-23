	subroutine mc_sos_recon (delta_p,delta_t,delta_phi,y_tgt,fry)
C+______________________________________________________________________________
!
! MC_SOS_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the MC_SOS program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Author: D. Potterveld, ANL, 18-Mar-1993
!
! Modification History:
!
! 2-August-1995	(RMM) Hardcoded in an extra index on the coeff, expon, and 
!		       n_terms variables to specify SOS. (SOS = 2)
!
!  19-AUG-1993	(DHP) Modified to use COSY INFINITY reconstruction coefficients.
C-______________________________________________________________________________
     
	implicit none

	include '../track.inc'

	integer*4 sos
	parameter (sos = 2)			!this is the SOS routine

C Argument definitions.

	real*8	delta_p,delta_t,delta_phi,y_tgt
	real*8	fry			!vertical position at target (+y=up)

C Cosy reconstruction matrix elements.

	integer*4 max_elements
	parameter (max_elements = 1000)

	integer*4 nspectr
	parameter (nspectr = 2)

	real*8		coeff(nspectr,4,max_elements)
	integer*2	expon(nspectr,4,max_elements)
	integer*4	n_terms(nspectr),max_order
	real*8		sum(4),hut(4),term

C Misc. variables.

	integer*4	i,j
	integer*4	chan

	logical*4	firstsos	/.true./
	character*132	line

C Functions.

!	integer*4	last_char
	logical*4	locforunt

C No amnesia, please...

	save

C ============================= Executable Code ================================

	fry=fry+0.	!avoid unused variable warning

C First time through, read in coefficients from data file.

	if (firstsos) then
	  if (.not.locforunt(chan)) stop 'MC_SOS_RECON: No I/O channels!'
	  open (unit=chan,status='old',readonly,file='sos/recon_cosy.dat')

! Read and print header.

	  line = '!'
	  do while (line(1:1).eq.'!')
	    read (chan,1001) line
!!!	    if (line(1:1).eq.'!') write(6,*) line(1:min(79,last_char(line)))
	  enddo

! Read in coefficients and exponents.

	  n_terms(sos) = 0
	  max_order = 0
	  do while (line(1:4).ne.' ---')
	    n_terms(sos) = n_terms(sos) + 1
	    if (n_terms(sos).gt.max_elements)
     >		stop 'WCRECON: too many COSY terms!'
	    read (line,1200) (coeff(sos,i,n_terms(sos)),i=1,4),
     >				(expon(sos,j,n_terms(sos)),j=1,4)
	    read (chan,1001) line
	    max_order = max(max_order, expon(sos,1,n_terms(sos)) +
     >			expon(sos,2,n_terms(sos)) +
     >			expon(sos,3,n_terms(sos)) +
     >			expon(sos,4,n_terms(sos)))
	  enddo
!!	  write(6,*) 'SOS: N_TERMS, MAX_ORDER = ',n_terms(sos),max_order
	  close (unit=chan)
	  firstsos = .false.
	endif

C Reset COSY sums.

	do i = 1,4
	  sum(i) = 0.
	enddo

C Convert hut quantities to right-handed coordinates, in meters and rad.

	hut(1) =  xs/100.		!Units: meters
	hut(2) =  dxdzs			!slope
	hut(3) =  ys/100.		!Meters
	hut(4) =  dydzs			!slope

C Compute COSY sums.

	do i = 1,n_terms(sos)
	  term = hut(1)**expon(sos,1,i)*hut(2)**expon(sos,2,i)
     >	       * hut(3)**expon(sos,3,i)*hut(4)**expon(sos,4,i)
	  sum(1) = sum(1) + term*coeff(sos,1,i)
	  sum(2) = sum(2) + term*coeff(sos,2,i)
	  sum(3) = sum(3) + term*coeff(sos,3,i)
	  sum(4) = sum(4) + term*coeff(sos,4,i)
	enddo
     
C Load output values.

	delta_phi = sum(1)		!slope
	y_tgt	  = sum(2)*100.		!cm
	delta_t   = sum(3)		!slope
	delta_p   = sum(4)*100.		!percent deviation

	return

C ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,4i1)

	END
