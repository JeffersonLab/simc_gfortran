	logical function check_dipole(x,y)
	implicit none

* Original version made on 04/01/98 by G.&I.Niculescu
* to more accurately model the size/shape of the HMS
* dipole ...

	include 'apertures.inc'
	real*8 x,y
	real*8 x_local,y_local
	logical check1,check2,check3,check4,check5,check6
	logical tmp_check

	check_dipole=.false.

* Let us observe first the obvious symmetry of the problem
* This helps reduce the checks to the first quadrant only...

	x_local=abs(x)
	y_local=abs(y)

* Now compare the current position and compare it with the different 
* apertures..

	check1=((x_local.le.x_d1).and.(y_local.le.y_d1))
	check2=((x_local.le.x_d2).and.(y_local.le.y_d2))
	check3=((x_local.le.x_d3).and.(y_local.le.y_d3))
	check4=((x_local.le.x_d4).and.(y_local.le.y_d4))

* now, the fifth check is the rounded corner

	check5=(((x_local-x_d5)**2+(y_local-y_d5)**2).le.r_d5**2)

* lastly the slanted piece

	check6=((x_local.ge.x_d4).and.(x_local.le.x_d3).and.
     >		((y_local-a_d6*x_local-b_d6).le.0.0))

* now, if we OR all the above we should get the inside of the can

	tmp_check=check1.or.check2.or.check3.or.check4.or.check5.or.check6

* for whatever reason mc_hms expects us to return the OUTSIDE of the can so...

	check_dipole = .not.tmp_check

	return
	end
