	subroutine musc(m2,p,rad_len,dth,dph)
C+_____________________________________________________________________
!
! MUSC - Simulate multiple scattering of any particle.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians, RAD_LEN in radiation
!   lengths. The formula used is due to Rossi and Greisen (See the book
!   by Segre, NUCLEI AND PARTICLES, 1982, p. 48.) The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are located at an angle with respect 
!   to the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!
! D. Potterveld - Sept. 1985
!
! Add option for protons 
!   Precision not good enough for our thick targets. Latest values supplied:
!   Lynch and Dahl, NIM B58 (1991) 6.
!   Note, that there is a typo in Particle data booklet: eps = 0.088 instead
!   of 0.2! Precision: -5% deviation for H2, +8 for 238U at radiation 
!   lengths between ~ 0.01 and 1.
!
! H.J. Bulten - Aug. 1991
C-_____________________________________________________________________

	implicit none

	real*8 Es, epsilon
	parameter (Es = 13.6)
	parameter (epsilon = 0.088)

	real*8 rad_len, dth, dph
	real*8 beta, theta_sigma
	real*8 m2, p

	real*8 nsig_max
	parameter(nsig_max=99.0d0)      !max #/sigma for gaussian ran #s.

	real*8 gauss1

! Compute scattering angles, THETA_SCAT from a gaussian distribution,
! PHI_SCAT from uniform distribution.

	beta = p / sqrt(m2+p*p)
	theta_sigma = Es/p/beta * sqrt(rad_len) * (1+epsilon*log10(rad_len))

! Compute new trajectory angles (units are rad)

	dth = dth + theta_sigma * gauss1(nsig_max)
	dph = dph + theta_sigma * gauss1(nsig_max)
	return
	end
