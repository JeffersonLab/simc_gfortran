C Aperture Limits:

	real*8 r_HBx,r_HBfyp,r_HBfym,r_HBmenyp,r_HBmenym
        real*8 r_HBmexyp,r_HBmexym,r_HBbyp,r_HBbym
        real*8 r_Q1,r_Q2,r_Q3
	real*8 r_D1
	real*8 r_BP

        parameter       (r_HBx         = 11.2)  !HB vertical radius 

        parameter       (r_HBfym        = -4.13) !HB entrance -y
        parameter       (r_HBfyp        = 11.75)!HB entrance +y

        parameter       (r_HBmenym      = -5.45) !HB mag entrance -y
                                        !(4.13+(20.06/115.22)*(11.71-4.13))
        parameter       (r_HBmenyp      = 11.74)!HB mag entrance +y

        parameter       (r_HBmexym      = -10.25) !HB mag exit -y
                                        !(4.13+(95.19/115.22)*(11.71-4.13))
        parameter       (r_HBmexyp      = 11.71)!HB mag exit +y

        parameter       (r_HBbym        = -11.71)!HB exit -y
        parameter       (r_HBbyp        = 11.70)!HB exit +y

	parameter	(r_Q1		= 20.00)	!20.00 cm radius
	parameter	(r_Q2		= 30.00)	!30.00 cm radius
	parameter	(r_Q3		= 30.00)	!30.00 cm radius
	parameter 	(r_D1		= 30.00)	!30 cm radius
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
	
* Implement the proper shape of the dipole aperture
*
	real*8 x_d1,y_d1
	real*8 x_d2,y_d2
	real*8 x_d3,y_d3
	real*8 x_d4,y_d4
	real*8 x_d5,y_d5,r_d5
	real*8 a_d6,b_d6
*
	parameter    (x_d1 = 31.5)	
	parameter    (x_d2 = 28.5)	
	parameter    (x_d3 = 13.97)	
	parameter    (x_d4 = 1.956)	
	parameter    (x_d5 = 28.5)	
	parameter    (y_d5 = 28.5)	
	parameter    (r_d5 =  6.35)
	parameter    (a_d6 = -0.114)
	parameter    (b_d6 = 20.54)
