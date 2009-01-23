	subroutine mc_hrsr_recon (delta_p,delta_t,delta_phi,y_tgt,fry)

C+______________________________________________________________________________
!
! MC_HRSR_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the MC_HRSR program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Inputs are from common block in track_*.inc (except for fry):
!  xs, ys, fry  are in cm.
!  dxdzs, dydzs are unitless slopes (we say "radians" to contrast to "mr").
! Matrix Elements want meters and radians, so have to convert cm-->m.
! Output:
!  delta_p is deviation from central momentum(unitless). Convert to % for output.
!  delta_t, delta_phi are in "radians".
!  y_tgt is in meters, convert to cm for output.
!
! Author: D. Potterveld, ANL, 18-Mar-1993
!
! Modification History:
!
! 2-August-1995 (RMM) Hardcoded in an extra index on the coeff, expon, and
!		      n_terms variables to specify HRSR. (HRSR = 3)
!
!  19-AUG-1993	(DHP) Modified to use COSY INFINITY reconstruction coefficients.
C-______________________________________________________________________________
     
	implicit none

	include '../spectrometers.inc'

	integer*4 specnum
	parameter (specnum = 3)			!this is the HRSR routine

C Argument definitions.

	real*8	delta_p,delta_t,delta_phi,y_tgt
	real*8  fry			!vertical position at target (+y=down)

C Cosy reconstruction matrix elements.

	integer*4	max_elements
	parameter	(max_elements = 1000)

	real*8		coeff(nspectr,4,max_elements)
	integer*2	expon(nspectr,5,max_elements)
	integer*4	n_terms(nspectr),max_order
	real*8		sum(4),hut(5),term

C Misc. variables.

	integer*4	i,j
	integer*4	chan

	logical		firsttime	/.true./
	character*132	line

C Functions.

	logical locforunt

C No amnesia, please...

	save

C ============================= Executable Code ================================

C First time through, read in coefficients from data file.

	if (firsttime) then
	   if (.not.locforunt(chan)) stop 'MC_HRSR_RECON: No I/O channels!'
c	   open (unit=chan,status='old',readonly,file='hrsr/hrs_recon_cosy.dat')
	   open (unit=chan,status='old',file='hrsr/hrs_recon_cosy.dat')

! Skip past header.

	   line = '!'
	   do while (line(1:1).eq.'!')
	      read (chan,1001) line
	   enddo

! Read in coefficients and exponents.
	   n_terms(specnum) = 0
	   max_order = 0
	   do while (line(1:4).ne.' ---')
	      n_terms(specnum) = n_terms(specnum) + 1
	      if (n_terms(specnum).gt.max_elements)
     >	      stop 'WCRECON: too many COSY terms!'
	      read (line,1200) (coeff(specnum,i,n_terms(specnum)),i=1,4),
     >			       (expon(specnum,j,n_terms(specnum)),j=1,5)
	      read (chan,1001) line
	      max_order = max(max_order, expon(specnum,1,n_terms(specnum)) +
     >				expon(specnum,2,n_terms(specnum)) +
     >				expon(specnum,3,n_terms(specnum)) +
     >				expon(specnum,4,n_terms(specnum)) +
     >				expon(specnum,5,n_terms(specnum)))
	   enddo
c	   write(6,*) 'HRSR:N_TERMS, MAX_ORDER = ',n_terms(specnum),max_order
	   close (unit=chan)
	   firsttime = .false.
	endif

C Reset COSY sums.

	do i = 1,4
	   sum(i) = 0.
	enddo

C Convert hut quantities to right-handed coordinates, in meters and "radians".
C Make sure hut(5) is non-zero, to avoid taking 0.0**0 (which crashes)

	hut(1) = xs/100.		!cm --> m
	hut(2) = dxdzs			!slope ("radians")
	hut(3) = ys/100.		!cm --> m
	hut(4) = dydzs			!slope ("radians")
	hut(5) = fry/100.		!vert. position at target(cm-->m)
	if (abs(hut(5)).le.1.d-30) hut(5)=1.d-30

C Compute COSY sums.

	do i = 1,n_terms(specnum)
	   term = hut(1)**expon(specnum,1,i)*hut(2)**expon(specnum,2,i)
     >		* hut(3)**expon(specnum,3,i)*hut(4)**expon(specnum,4,i)
     >		* hut(5)**expon(specnum,5,i)
	   sum(1) = sum(1) + term*coeff(specnum,1,i)
	   sum(2) = sum(2) + term*coeff(specnum,2,i)
	   sum(3) = sum(3) + term*coeff(specnum,3,i)
	   sum(4) = sum(4) + term*coeff(specnum,4,i)
	enddo
     
C Load output values.

	delta_phi = sum(1)		!slope ("radians")
	y_tgt	  = sum(2)*100.		!m --> cm
	delta_t   = sum(3)		!slope ("radians")
	delta_p   = sum(4)*100.		!percent deviation

	return

C ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,5i1)

      END
