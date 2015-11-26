	subroutine mc_hms (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen)

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

	include 'apertures_hms.inc'
	include 'struct_hms.inc'

	include '../constants.inc'
	include '../spectrometers.inc'

*G.&I.N. STUFF - for checking dipole exit apertures and vacuum pipe to HMS hut.
	real*8 x_offset_pipes/2.8/,y_offset_pipes/0.0/

C Spectrometer definitions - for double arm monte carlo compatability

	integer*4 spectr
	parameter (spectr = 1)		!HMS is spec #1.

! Collimator (octagon) dimensions and offsets.

	real*8 h_entr,v_entr,h_exit,v_exit
	real*8 x_off,y_off,z_off

! Option for mock sieve slit.  Just take particles and forces trajectory
! to put the event in the nearest hole.  Must choose the "No collimator"
! values for the h(v)_entr and h(v)_exit values.
! Note that this will mess up the physics distributions somewhat, but it
! should still be pretty good for optics. Physics limits (e.g. elastic
! peak at x<=1) will not be preserved.

	logical use_sieve /.false./		!use a fake sieve slit.

! No collimator - wide open
!	parameter (h_entr = 99.)
!	parameter (v_entr = 99.)
!	parameter (h_exit = 99.)
!	parameter (v_exit = 99.)

! Old collimator for HMS-1 tune (called 'large' or 'pion' at the time).
!	parameter (h_entr = 3.536)
!	parameter (v_entr = 9.003)
!	parameter (h_exit = 3.708)
!	parameter (v_exit = 9.444)

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

! z-position of important apertures.
	real*8 z_entr,z_exit
	parameter (z_entr = 126.2e0 + z_off)	!nominally 1.262 m
	parameter (z_exit = z_entr + 6.3e0)	!6.3 cm thick

	real*8 z_dip1,z_dip2,z_dip3		!post dipole apertures
	parameter (z_dip1 = 64.77e0)		! end of 26.65 inch pipe.
	parameter (z_dip2 = z_dip1 + 297.18e0)	!~117 inch pipe.
	parameter (z_dip3 = z_dip2 + 115.57e0)	! 45.5 inch pipe.

! Math constants

	real*8 d_r,r_d
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)

! The arguments

	real*8 x,y,z				!(cm)
	real*8 dpp				!delta p/p (%)
	real*8 dxdz,dydz			!X,Y slope in spectrometer
	real*8 x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 fry				!vertical position@tgt (+y=down)
	real*8 pathlen
	logical ms_flag				!mult. scattering flag
	logical wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical ok_spec				!true if particle makes it

! Local declarations.

	integer*4	chan	/1/,n_classes

	logical	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt,rt,tht			!temporaries
	real*8 resmult				!DC resolution factor
	real*8 zdrift

	logical dflag			!has particle decayed?
	logical ok

	real*8 grnd

! Gaby's dipole shape stuff
	logical hit_dipole
	external hit_dipole

	save		!Remember it all!

C ================================ Executable Code =============================

! Initialize ok_spec to .false., reset decay flag

	ok_spec = .false.
	dflag = .false.			!particle has not decayed yet
	hSTOP_trials = hSTOP_trials + 1
	xt = th_spec	!avoid 'unused variable' error for th_spec

! Force particles to go through the sieve slit holes, for mock sieve option.

	if (use_sieve) then
	  xt = x + z_entr * dxdz	!project to collimator
	  yt = y + z_entr * dydz
	  xt = 2.540*nint(xt/2.54)	!shift to nearest hole.
	  yt = 1.524*nint(yt/1.524)
	  rt = 0.254*sqrt(grnd())	!distance from center of hole(r=2.54mm)
	  tht= 2*pi*grnd()		!angle of offset.
	  xt = xt + rt*cos(tht)
	  yt = yt + rt*sin(tht)
	  dxdz = (xt-x)/z_entr		!force to correct angle.
	  dydz = (yt-y)/z_entr
	endif

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

	if (first_time_here) then
	   call transp_init(spectr,n_classes)
	   close (unit=chan)
	   if (n_classes.ne.12) stop 'MC_HMS, wrong number of transport classes'
	   first_time_here = .false.
	endif

        
C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

! Begin transporting particle.

! Do transformations, checking against apertures.
! Check front of fixed slit.
	  zdrift = z_entr
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if (abs(ys-y_off).gt.h_entr) then
	    hSTOP_slit_hor = hSTOP_slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_entr) then
	    hSTOP_slit_vert = hSTOP_slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_entr/h_entr*abs(ys-y_off)+3*v_entr/2)) then
	    hSTOP_slit_oct = hSTOP_slit_oct + 1
	    goto 500
	  endif

!Check back of fixed slit.

	  zdrift = z_exit-z_entr
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if (abs(ys-y_off).gt.h_exit) then
	    hSTOP_slit_hor = hSTOP_slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_exit) then
	    hSTOP_slit_vert = hSTOP_slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs-x_off).gt. (-v_exit/h_exit*abs(ys-y_off)+3*v_exit/2)) then
	    hSTOP_slit_oct = hSTOP_slit_oct + 1
	    goto 500
	  endif

! Go to Q1 IN  mag bound.  Drift rather than using COSY matrices

	  if (.not.adrift(spectr,1)) write(6,*) 'Transformation #1 is NOT a drift'
	  zdrift = driftdist(spectr,1) - z_exit
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_in = hSTOP_Q1_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q1.

	  call transp(spectr,2,decay_flag,dflag,m2,p,125.233e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_mid = hSTOP_Q1_mid + 1
	    goto 500
	  endif

! Go to Q1 OUT mag boundary.

	  call transp(spectr,3,decay_flag,dflag,m2,p,62.617e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    hSTOP_Q1_out = hSTOP_Q1_out + 1
	    goto 500
	  endif

! Go to Q2 IN  mag bound.  Drift rather than using COSY matrices

	  if (.not.adrift(spectr,4)) write(6,*) 'Transformation #4 is NOT a drift'
	  zdrift = driftdist(spectr,4)
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_in = hSTOP_Q2_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q2.

	  call transp(spectr,5,decay_flag,dflag,m2,p,143.90e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_mid = hSTOP_Q2_mid + 1
	    goto 500
	  endif

! Go to Q2 OUT mag boundary.

	  call transp(spectr,6,decay_flag,dflag,m2,p,71.95e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    hSTOP_Q2_out = hSTOP_Q2_out + 1
	    goto 500
	  endif

! Go to Q3 IN  mag bound.  Drift rather than using COSY matrices

	  if (.not.adrift(spectr,7)) write(6,*) 'Transformation #7 is NOT a drift'
	  zdrift = driftdist(spectr,7)
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_in = hSTOP_Q3_in + 1
	    goto 500
	  endif

! Check aperture at 2/3 of Q3.

	  call transp(spectr,8,decay_flag,dflag,m2,p,143.8e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_mid = hSTOP_Q3_mid + 1
	    goto 500
	  endif

! Go to Q3 OUT mag boundary.

	  call transp(spectr,9,decay_flag,dflag,m2,p,71.9e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    hSTOP_Q3_out = hSTOP_Q3_out + 1
	    goto 500
	  endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.
! Aperture has complicated form (evaluated by function hit_dipole).

	  if (.not.adrift(spectr,10)) write(6,*) 'Transformation #10 is NOT a drift'
	  zdrift = driftdist(spectr,10)
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  xt=xs
	  yt=ys
	  call rotate_haxis(-6.0e0,xt,yt)
	  if (hit_dipole(xt,yt)) then
	    hSTOP_D1_in = hSTOP_D1_in + 1
	    goto 500
	  endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(spectr,11,decay_flag,dflag,m2,p,526.053e0,pathlen)
	  xt=xs
	  yt=ys
	  call rotate_haxis(6.0e0,xt,yt)
	  if (hit_dipole(xt,yt)) then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    goto 500
	  endif

! Check a number of apertures in the vacuum pipes following the
! dipole.  First the odd piece interfacing with the dipole itself

	  if ( (((xt-x_offset_pipes)**2+(yt-y_offset_pipes)**2).gt.30.48**2)
     >           .or. (abs((yt-y_offset_pipes)).gt.20.5232) ) then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    goto 500
	  endif

! Check the exit of the 26.65 inch pipe

	  zdrift = z_dip1
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.1145.518)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    goto 500
	  endif

! check exit of long (117 inch) pipe (entrance is bigger than previous pipe)
! note: Rolf claims its 117.5 but the dravings say more like 116.x
! .. so i put 117 even.  Should be a 30.62 diameter pipe

	  zdrift = z_dip2 - z_dip1
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.1512.2299)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    goto 500
	  endif

! lastly check the exit of the last piece of pipe. 45.5 inches long, 30.62 dia.

	  zdrift = z_dip3 - z_dip2
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if (((xs-x_offset_pipes)**2+(ys-y_offset_pipes)**2).gt.2162.9383)then
	    hSTOP_D1_out = hSTOP_D1_out + 1
	    goto 500
	  endif

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_hms_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_hms_hut knows
! where to start.  Initial zposition for mc_hms_hut is -147.48 cm so that the
! sum of four drift lengths between pipe and focal plane is 625.0 cm
! (64.77+297.18+115.57+147.48=625)

! If we get this far, the particle is in the hut.

	  hSTOP_hut = hSTOP_hut + 1

! and track through the detector hut

	  if (.not.adrift(spectr,12)) write (6,*) 'Transformation #12 is NOT a drift'

	  zdrift = driftdist(spectr,12) - z_dip3	!distance left to go.
	  call mc_hms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-zdrift,pathlen)
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

	  ok_spec = .true.
	  hSTOP_successes = hSTOP_successes + 1

! We are done with this event, whether GOOD or BAD.

500	continue

	return
	end

*-----------------------------------------------------------------------

	logical function hit_dipole(x,y)
	implicit none

* Original version made on 04/01/98 by G.&I.Niculescu
* to more accurately model the size/shape of the HMS
* dipole ...

	include 'apertures_hms.inc'
	real*8 x,y
	real*8 x_local,y_local
	logical check1,check2,check3,check4,check5,check6
	logical miss_dipole

	hit_dipole=.false.

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

	miss_dipole=check1.or.check2.or.check3.or.check4.or.check5.or.check6

* for whatever reason mc_hms expects us to return the OUTSIDE of the can so...

	hit_dipole = .not.miss_dipole

	return
	end
