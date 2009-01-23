C Aperture Limits:

	real*8 r_Q1,r_Q2,r_Q3

	parameter	(r_Q1		= 20.50)	!20.50 cm radius
	parameter	(r_Q2		= 30.22)	!30.22 cm radius
	parameter	(r_Q3		= 30.22)	!30.22 cm radius
*
* G. & I. Niculescu
*
* implement the proper shape of the dipole aperture
*
	real*8 x_d1,y_d1
	real*8 x_d2,y_d2
	real*8 x_d3,y_d3
	real*8 x_d4,y_d4
	real*8 x_d5,y_d5,r_d5
	real*8 a_d6,b_d6
*
	parameter    (x_d1 = 34.29)	
	parameter    (y_d1 = 12.07)
	parameter    (x_d2 = 27.94)	
	parameter    (y_d2 = 18.42)
	parameter    (x_d3 = 13.97)	
	parameter    (y_d3 = 18.95)
	parameter    (x_d4 = 1.956)	
	parameter    (y_d4 = 20.32)
	parameter    (x_d5 = 27.94)	
	parameter    (y_d5 = 12.065)
	parameter    (r_d5 =  6.35)
	parameter    (a_d6 = -0.114)
	parameter    (b_d6 = 20.54)
