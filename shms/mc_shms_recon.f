	subroutine mc_shms_recon (delta_p,delta_t,delta_phi,y_tgt,
     >                            fry,spectr)
C+______________________________________________________________________________
!
! MC_HMS_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the MC_HMS program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Author: D. Potterveld, ANL, 18-Mar-1993
!
! Modification History:
!
! 2-August-1995 (RMM) Hardcoded in an extra index on the coeff, expon, and
!                      n_terms variables to specify HMS. (HMS = 1)
!
!  19-AUG-1993	(DHP) Modified to use COSY INFINITY reconstruction coefficients.
C-______________________________________________________________________________

	implicit none

	include '../spectrometers.inc'

C Spectrometer definitions - for double arm monte carlo compatability
	integer*4 spectr

C Argument definitions.

	real*8	delta_p,delta_t,delta_phi,y_tgt
	real*8	fry			!vertical position at target (+y=down)

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

	logical*4	firsttime	/.true./
	character*132	line

C Functions.

	logical*4	locforunt

C No amnesia, please...

	save

C ============================= Executable Code ================================
C setting some temporary files


C First time through, read in coefficients from data file.

	if (firsttime) then
	   if (.not.locforunt(chan)) 
     >         stop 'MC_SHMS_RECON: No I/O channels!'
	   if (spectr.eq.5) then	!ssa tune
! COSY calculated (not so great!!!)
c	     open (unit=chan,status='old',file='shms/shms_hsa_2009_recon_cosy_daveme2.dat')
c	     open (unit=chan,status='old',file='shms/shms_recon_cosy_2011_dipole26cm_dm1.2_nov9.dat')
c	     open (unit=chan,status='old',file='shms/shms_recon_refit_5th_order.dat')
c	     open (unit=chan,status='old',file='shms/shms_recon_fit_90deg_1cm_5th_order.dat')
	     open (unit=chan,status='old',file='shms/shms_recon.dat')
c	     open (unit=chan,status='old',file='shms/shms_hsa_2009_recon_cosy.dat')
! CMOP REFIT
c	     open (unit=chan,status='old',file=
c     >'shms/cmop_refit_shms_lq_qdi_hsa_split_recon_newfit_4thorder.dat')
! TH new optics
!	     open (unit=chan,status='old',
!     > file='shms/shms2008_rec_th.dat')
!     > file='shms/shms2008_rec_th_pmag7.dat')
	   else if (spectr.eq.6) then	!lsa tune
	      write(6,*) 'mc_shms_recon: 
     >You are trying to use the LSA: no banana!'
c	     open (unit=chan,status='old',readonly,file='shms/shms_recon_cosy_LSA.dat')
	   else
	     write(6,*) 'MC_SHMS_RECON: I just cant 
     >handle spectr=',spectr
	     stop
	   endif

! Read and print header.

	   line = '!'
	   do while (line(1:1).eq.'!')
	      read (chan,1001) line
	   enddo

! Read in coefficients and exponents.
*	   write (6,*) 'reading in coeffiecients and exponents in the shms_recon'
	   n_terms(spectr) = 0
	   max_order = 0
	   do while (line(1:4).ne.' ---')
	      n_terms(spectr) = n_terms(spectr) + 1
	      if (n_terms(spectr).gt.max_elements)
     >	      stop 'WCRECON: too many COSY terms!'
	      read (line,1200) (coeff(spectr,i,n_terms(spectr)),i=1,4),
     >			       (expon(spectr,j,n_terms(spectr)),j=1,5)
	      read (chan,1001) line
	      max_order = max(max_order,expon(spectr,1,n_terms(spectr))+
     >				expon(spectr,2,n_terms(spectr)) +
     >				expon(spectr,3,n_terms(spectr)) +
     >				expon(spectr,4,n_terms(spectr)) +
     >		                expon(spectr,5,n_terms(spectr)))
	   enddo
!	   type *,' N_TERMS, MAX_ORDER = ',n_terms(spectr),max_order
	   close (unit=chan)
	   firsttime = .false.
	endif

C Reset COSY sums.

	do i = 1,4
	   sum(i) = 0.
	enddo

C Convert hut quantities to right-handed coordinates, in meters and rad.

*	write (6,*) 'In recon file, and converting hut quantities'
	hut(1) =  xs/100.		!Units: meters
	hut(2) =  dxdzs			!Radians
	hut(3) =  ys/100.		!Meters
	hut(4) =  dydzs			!Radians
	hut(5) = fry/100.		!vert. position at target (cm-->m)
	if (abs(hut(5)).le.1.e-30) hut(5)=1.e-30


C Compute COSY sums.

	do i = 1,n_terms(spectr)
	   term = hut(1)**expon(spectr,1,i)*hut(2)**expon(spectr,2,i)
     >	         * hut(3)**expon(spectr,3,i)*hut(4)**expon(spectr,4,i)
     >		* hut(5)**expon(spectr,5,i)
	   sum(1) = sum(1) + term*coeff(spectr,1,i)
	   sum(2) = sum(2) + term*coeff(spectr,2,i)
	   sum(3) = sum(3) + term*coeff(spectr,3,i)
	   sum(4) = sum(4) + term*coeff(spectr,4,i)
	enddo

C Load output values.

*	write (6,*) 'loading output values in the reconstruct file'
	delta_phi = sum(1)		!radians
	y_tgt	  = sum(2)*100.		!cm
	delta_t   = sum(3)		!radians
	delta_p   = sum(4)*100.		!percent deviation

      return

C ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,5i1)

      END
