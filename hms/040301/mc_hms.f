	subroutine mc_hms (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, tg_rad_len,
     >		m2, ms_flag, wcs_flag, decay_flag, resmult, fry, ok_hms)

C+______________________________________________________________________________
!
! Monte-Carlo of HMS spectrometer.
!	Note that we only pass on real*8 variables to the subroutine.
!	This will make life easier for portability.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!		which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
!
!  15-SEP-1997  MODIFY stepping through spectrometer so that all drifts
!		use project.f (not transp.f), and so that project.f and
!		transp.f take care of decay.  Decay distances all assume
!		the pathlength for the CENTRAL RAY.
C-______________________________________________________________________________

	implicit 	none

	include 'apertures.inc'
	include '../struct_hms.inc'
	include '../track.inc'
	include '../constants.inc'

*G.&I.N. STUFF - for checking dipole exit apertures and vacuum pipe to HMS hut.
	real*8 x_offset_pipes/0.0/,y_offset_pipes/0.0/

! Spectrometer definitions

	integer*4 hms,sos
	parameter (hms = 1)
	parameter (sos = 2)

! Collimator (octagon) dimensions and offsets.

	real*8 h_entr,v_entr,h_exit,v_exit
	real*8 x_off,y_off,z_off

! Open values (no collimator)
c	parameter (h_entr = 20.0)
c	parameter (v_entr = 20.0)
c	parameter (h_exit = 20.0)
c	parameter (v_exit = 20.0)

! Old collimator for HMS-1 tune (called 'large' or 'pion' at the time).
c	parameter (h_entr = 3.536)
c	parameter (v_entr = 9.003)
c	parameter (h_exit = 3.708)
c	parameter (v_exit = 9.444)

! New collimator for HMS-100 tune.
	parameter (h_entr = 4.575)
	parameter (v_entr = 11.646)
	parameter (h_exit = 4.759)
	parameter (v_exit = 12.114)

c	parameter (x_off=+0.496)	!+ve is slit DOWN - HMS1 tune
c	parameter (x_off=+0.126)	!HMS-100 (preliminary survey)
	parameter (x_off=+0.000)	!HMS-100 (zeroed before Fpi)

c	parameter (y_off=-0.004)	!+ve is slit LEFT (as seen from target)
	parameter (y_off=+0.028)	!HMS-100 (number from gaskell)
	
c	parameter (z_off=+0.00)		!1995 position
c	parameter (z_off=+1.50)		!1996 position
	parameter (z_off=+40.17)	!HMS100 tune (dg 5/27/98)

! Math constants

	real*8 d_r,r_d,root
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2

! The arguments

	real*8 x,y,z				!(cm)
	real*8 dpp				!delta p/p (%)
	real*8 dxdz,dydz			!X,Y slope in spectrometer
	real*8 x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 tg_rad_len			!target length in r.l.
	real*8 fry				!vertical position@tgt (+y=down)
	logical ms_flag				!mult. scattering flag
	logical wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical ok_hms				!true if particle makes it

! Local declarations.

	integer*4	chan	/1/,n_classes

	logical	first_time_hms/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt				!temporaries
	real*8 resmult				!DC resolution factor
	real*8 dummy

	logical dflag			!has particle decayed?
	logical ok

! Gaby's dipole shape stuff
	logical check_dipole
	external check_dipole

	save		!Remember it all!

C ================================ Executable Code =============================

! Initialize ok_hms to .false., reset decay flag

	ok_hms = .false.
	dflag = .false.   !particle has not decayed yet
	hSTOP.trials = hSTOP.trials + 1

! Save spectrometer coordinates.

	xs = x
	ys = y
	zs = z
	dxdzs = dxdz
	dydzs = dydz

! particle momentum

	dpps = dpp
	p = p_spec*(1.+dpps/100.)

! Read in transport coefficients.

	if (first_time_hms) then
	   call transp_init(hms,n_classes)
	   close (unit=chan)
	   if (n_classes.ne.12) stop 'MC_HMS, wrong number of transport classes'
	   first_time_hms = .false.
	endif

! Calculate multiple scattering in target
c I've commented this out because it's already taken care of in
c  the SIMULATE7 main code (sub montecarlo() )  --rmm
c	if(ms_flag) call musc(m2,p,tg_rad_len,dydzs,dxdzs)
	dummy = tg_rad_len	!avoid unused variable warning

! Begin transporting particle.
! Do transformations, checking against apertures.

! Check front of fixed slit, at about 1.26 meter

	  call project(xs,ys,126.2d0+z_off,decay_flag,dflag,m2,p) !project and decay
	  if (abs(ys-y_off).gt.h_entr) then
	    hSTOP.slit_hor = hSTOP.slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_entr) then
	    hSTOP.slit_vert = hSTOP.slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_entr/h_entr*abs(ys-y_off)+3*v_entr/2)) then
	    hSTOP.slit_oct = hSTOP.slit_oct + 1
	    goto 500
	  endif
!Check back of fixed slit, at about 1.325 meter

	  call project(xs,ys,6.3d0,decay_flag,dflag,m2,p) !project and decay
	  if (abs(ys-y_off).gt.h_exit) then
	    hSTOP.slit_hor = hSTOP.slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_exit) then
	    hSTOP.slit_vert = hSTOP.slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_exit/h_exit*abs(ys-y_off)+3*v_exit/2)) then
	    hSTOP.slit_oct = hSTOP.slit_oct + 1
	    goto 500
	  endif

! Go to Q1 IN  mag bound.  Drift rather than using COSY matrices

	  call project(xs,ys,(216.075d0-126.2d0-z_off-6.3d0),decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP.Q1_in = hSTOP.Q1_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q1.

	  call transp(hms,2,decay_flag,dflag,m2,p,126.0d0)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP.Q1_mid = hSTOP.Q1_mid + 1
	    goto 500
	  endif

! Go to Q1 OUT mag boundary.

	  call transp(hms,3,decay_flag,dflag,m2,p,63.0d0)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP.Q1_out = hSTOP.Q1_out + 1
	    goto 500
	  endif

! Go to Q2 IN  mag bound.  Drift rather than using COSY matrices
!!	  call transp(hms,4,decay_flag,dflag,m2,p)

	  call project(xs,ys,123.15d0,decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP.Q2_in = hSTOP.Q2_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q2.

	  call transp(hms,5,decay_flag,dflag,m2,p,143.67d0)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP.Q2_mid = hSTOP.Q2_mid + 1
	    goto 500
	  endif

! Go to Q2 OUT mag boundary.

	  call transp(hms,6,decay_flag,dflag,m2,p,71.833d0)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP.Q2_out = hSTOP.Q2_out + 1
	    goto 500
	  endif

! Go to Q3 IN  mag bound.  Drift rather than using COSY matrices
!!	  call transp(hms,7,decay_flag,dflag,m2,p)

	  call project(xs,ys,94.225d0,decay_flag,dflag,m2,p) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP.Q3_in = hSTOP.Q3_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q3.

	  call transp(hms,8,decay_flag,dflag,m2,p,145.7d0)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP.Q3_mid = hSTOP.Q3_mid + 1
	    goto 500
	  endif

! Go to Q3 OUT mag boundary.

	  call transp(hms,9,decay_flag,dflag,m2,p,72.9d0)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP.Q3_out = hSTOP.Q3_out + 1
	    goto 500
	  endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.
! Aperture has elliptical form.
!!	  call transp(hms,10,decay_flag,dflag,m2,p)

	  call project(xs,ys,102.15d0,decay_flag,dflag,m2,p) !project and decay
	  call rotate_haxis(-6.0d0,xt,yt)
	  if (check_dipole(xt,yt)) then
	    hSTOP.D1_in = hSTOP.D1_in + 1
	    goto 500
	  endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(hms,11,decay_flag,dflag,m2,p,526.1d0)
	  call rotate_haxis(6.0d0,xt,yt)
	  if (check_dipole(xt,yt)) then
	    hSTOP.D1_out = hSTOP.D1_out + 1
	    goto 500
	  endif

! Check a number of apertures in the vacuum pipes following the
! dipole.  First the odd piece interfacing with the dipole itself

	  if ( (((xt-x_offset_pipes)**2+(yt-y_offset_pipes)**2).gt.30.48**2)
     >           .or. (abs((yt-y_offset_pipes)).gt.20.5232) ) then
	    hSTOP.D1_out = hSTOP.D1_out + 1
	    goto 500
	  endif

! Check the exit of the 26.65 inch pipe

	  call project(xs,ys,64.77d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.1145.518)then
	    hSTOP.D1_out = hSTOP.D1_out + 1
	    goto 500
	  endif

! check exit of long (117 inch) pipe (entrance is bigger than previous pipe)
! note: Rolf claims its 117.5 but the dravings say more like 116.x
! .. so i put 117 even.  Should be a 30.62 diameter pipe

	  call project(xs,ys,297.18d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.1512.2299)then
	    hSTOP.D1_out = hSTOP.D1_out + 1
	    goto 500
	  endif

! lastly check the exit of the last piece of pipe. 45.5 inches long, 30.62 dia.

	  call project(xs,ys,+115.57d0,decay_flag,dflag,m2,p) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.2162.9383)then
	    hSTOP.D1_out = hSTOP.D1_out + 1
	    goto 500
	  endif

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_hms_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_hms_hut knows
! where to start.  Initial zposition for mc_hms_hut is -147.48 cm so that the
! sum of four drift lengths between pipe and focal plane is 625.0 cm
! (64.77+297.18+115.57+147.48=625)

! If we get this far, the particle is in the hut.

	  hSTOP.hut = hSTOP.hut + 1

! and track through the detector hut

	  call mc_hms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-147.48d0)
	  if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.
	  xs=x_fp
	  ys=y_fp
	  dxdzs=dx_fp
	  dydzs=dy_fp

! Reconstruct target quantities.
	  call mc_hms_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)

! Fill output to return to main code
	  dpp = dpp_recon
	  dxdz = dph_recon
	  dydz = dth_recon
	  y = y_recon

	  ok_hms = .true.
	  hSTOP.successes = hSTOP.successes + 1

! We are done with this event, whether GOOD or BAD.

500	continue

	return
	end
