C Aperture Limits:

	real*8 r_Q1,r_Q2
	real*8 r_D1
	real*8 r_BP

	parameter	(r_Q1		= 20.50)	!20.50 cm radius
	parameter	(r_Q2		= 20.50)	!20.50 cm radius
	parameter 	(r_D1		= 30)		!30 cm radius
	parameter	(r_BP		= 50)		!50 cm radius for flare

C Trajectory Offsets:

	real*8 off_bp,ang_bp_in,ang_bp_out
	real*8 off_d1,ang_d1_in,ang_d1_out
	real*8 off_mid,ang_mid

	parameter	(off_bp		= -19.441)	! measured in ?
	parameter	(ang_bp_in	= 9.20)		! degrees
	parameter	(ang_bp_out	= -9.20)	! degrees	
	parameter	(off_d1		= -11.365)	! measured in ?
	parameter	(ang_d1_in	= 9.049)	! degrees
	parameter	(ang_d1_out	= -9.049)	! degrees
	parameter	(off_mid	= 7.00)		! measured in ?
	parameter 	(ang_mid	= 0)		! degrees
	
