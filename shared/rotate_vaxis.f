	subroutine rotate_vaxis(rotang,xp0,yp0)

C+______________________________________________________________________________
!
! ROTATE_VAXIS- Calculate new trajectory coordinates in reference frame rotated
!   about vertical axis by angle ROTANG relative to central ray.
!
! *** Right-handed TRANSPORT coordinates are assumed! ***
!
!   ROTANG is an angle about the  negative X-axis.
!
!   Input trajectory is: Y = YS + ALPHA*(Z-ZS)
!                        X = XS + BETA *(Z-ZS)
!                        Z = ZS is current point.
!
!   Output traject is:  YP = YP0 + ALPHA_P*ZP
!                       XP = XP0 + BETA_P *ZP
!                       ZP = 0 gives intersection of track with rotated plane.
!
! ROTANG (R*4):	Rotation angle in degrees.
!
! D. Potterveld, 15-Mar-1993.
C-______________________________________________________________________________

	implicit none

        include '../spectrometers.inc'

	real*8 rotang,xp0,yp0,yi
	real*8 alpha,beta,alpha_p,beta_p,sin_th,cos_th,tan_th
	real*8 rotang_rad, raddeg

	save

C ============================= Executable Code ================================

C Sep. 2008 DJG: No equivaleent to sind etc. in gfortran. Convert to radians first.


	parameter (raddeg=0.017453292) 

	rotang_rad = rotang*raddeg


1	tan_th = tan(rotang_rad)
	sin_th = sin(rotang_rad)
	cos_th = cos(rotang_rad)

	alpha  = dydzs
	beta   = dxdzs

	alpha_p= (alpha + tan_th)/(1. - alpha*tan_th)
	beta_p = beta/(cos_th - alpha*sin_th)

        yi = yp0
	yp0    = yi*(cos_th + alpha_p*sin_th)
	xp0    = xp0 + yi*beta_p*sin_th

	return
	end
