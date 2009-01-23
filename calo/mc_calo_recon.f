	subroutine mc_calo_recon (delta_p,delta_t,delta_phi,y_tgt,fry,delta_y,delta_z,drift_to_cal)
C+______________________________________________________________________________
!
! MC_CALO_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the MC_HMS program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Inputs are from common block in track_*.inc (except for fry):
!  xs, ys, fry  are in cm.
!  dxdzs, dydzs are unitless slopes (we say "radians" to contrast to "mr").
C-______________________________________________________________________________
     
	implicit none

	include '../spectrometers.inc'

	integer*4 specnum
	parameter (specnum = 5)			!this is the CALO routine

C Argument definitions.

	real*8	delta_p,delta_t,delta_phi,y_tgt
	real*8	fry			!vertical position at target (+y=down)
	real*8  ztemp,drift_to_cal,delta_y,delta_z

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

C Load output values.

	delta_p = dpps
	delta_phi = (xs-fry)/(drift_to_cal-delta_z)
	delta_t = (ys-delta_y)/(drift_to_cal-delta_z) 

c mkj cannot determine y_tgt from the calorimeter.
	y_tgt = -100.
      return

C ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,5i1)

      END

