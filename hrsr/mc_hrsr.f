	subroutine mc_hrsr (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen)

C+______________________________________________________________________________
!
! Monte-Carlo of HRSR spectrometer. (USED TO BE 'HADRON ARM').
!
!
! Author: David Meekins April 2000
! based on code by David Potterveld
!
! Modification History:
!
!  units for this file are percents, cm, and mrads.
!
C-______________________________________________________________________________

	implicit 	none

	include 'struct_hrsr.inc'
	include 'apertures_hrsr.inc'
	include 'g_dump_all_events.inc'

	include '../constants.inc'
	include '../spectrometers.inc'

C Spectrometer definitions - for double arm monte carlo compatability

	integer*4 spectr
	parameter (spectr = 3)	!hrs-right is spec #3.

C Collimator (rectangle) dimensions and offsets.

	real*8  h_entr,v_entr	!horiz. and vert. 1/2 gap of fixed slit
	real*8  h_exit,v_exit	!horiz. and vert. 1/2 gap of fixed slit
	real*8  y_off,x_off	!horiz. and vert. position offset of slit
	real*8  z_off		!offset in distance from target to front of sli+

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

! Large collimator: (hallaweb.jlab.org/news/minutes/collimator-distance.html)
	parameter (h_entr = 3.145)
	parameter (v_entr = 6.090)
	parameter (h_exit = 3.340)	!0.1mm wider than 'electron arm'
	parameter (v_exit = 6.485)

	parameter (y_off  = 0.0)
	parameter (x_off  = 0.0)
	parameter (z_off  = 0.0)

! z-position of important apertures.
	real*8 z_entr,z_exit
	parameter (z_entr = 110.0d0 + z_off)	!nominally 1.100 m
	parameter (z_exit = z_entr + 8.0d0)	!8.0 cm thick

C Math constants

	real*8 d_r,r_d
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)

C The arguments

	real*8	x,y,z				!(cm)
	real*8	dpp				!delta p/p (%)
	real*8 	dxdz,dydz			!X,Y slope in spectrometer
	real*8	x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8	p_spec, th_spec			!spectrometer setting
	real*8  fry                         	!vertical position@tgt (+y=down)
	real*8	pathlen
	logical	ms_flag				!mult. scattering flag
	logical	wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical	ok_spec				!true if particle makes it

C Local declarations.

	integer*4	chan/1/,n_classes

	logical*4	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt,rt,tht			!temporaries
	real*8 resmult				!DC resolution factor (unused)
	real*8 zdrift,ztmp

	logical dflag			!has particle decayed?
	logical	ok

	real*8 grnd

	save		!Remember it all!

C ================================ Executable Code =============================

! Initialize ok_spec to .flase., restet decay flag

	ok_spec = .false.
	dflag = .false.			!particle has not decayed yet
	rSTOP_trials = rSTOP_trials + 1
	xt = th_spec    !avoid 'unused variable' error for th_spec

! Force particles to go through the sieve slit holes, for mock sieve option.
! Use z_exit, since sieve slit is ~8cm behind normal collimator.

	if (use_sieve) then
	  xt = x + z_exit * dxdz	!project to collimator
	  yt = y + z_exit * dydz
	  xt = 2.50*nint(xt/2.50)	!shift to nearest hole.
	  yt = 1.25*nint(yt/1.25)
	  rt = 0.1*sqrt(grnd())		!distance from center of hole(r=1.0mm)
	  tht= 2*pi*grnd()		!angle of offset.
	  dxdz = (xt-x)/z_exit		!force to correct angle.
	  dydz = (yt-y)/z_exit
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

C Read in transport coefficients.

	if (first_time_here) then
	  call transp_init(spectr,n_classes)
	  close (unit=chan)
	  if (n_classes.ne.12) stop 'MC_HRSR, wrong number of transport classes'
	  first_time_here = .false.
	endif

! Begin transporting particle.
! Do transformations, checking against apertures.

! Circular apertures before slitbox (only important for no collimator)
	zdrift = 65.686
	ztmp = zdrift
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if (sqrt(xs*xs+ys*ys).gt.7.3787) then
	  rSTOP_slit_hor = rSTOP_slit_hor + 1
	  stop_where=20.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 80.436 - ztmp
	ztmp = 80.436
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if (sqrt(xs*xs+ys*ys).gt.7.4092) then
	  rSTOP_slit_hor = rSTOP_slit_hor + 1
	  stop_where=21.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif
	

! Check front of fixed slit.

	zdrift = z_entr - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if (abs(ys-y_off).gt.h_entr) then
	  rSTOP_slit_hor = rSTOP_slit_hor + 1
	  stop_where=1.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif
	if (abs(xs-x_off).gt.v_entr) then
	  rSTOP_slit_vert = rSTOP_slit_vert + 1
	  stop_where=2.	
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check back of fixed slit.

	zdrift = z_exit - z_entr
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	if (abs(ys-y_off).gt.h_exit) then
	  rSTOP_slit_hor = rSTOP_slit_hor + 1
	  stop_where=3.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif
	if (abs(xs-x_off).gt.v_exit) then
	  rSTOP_slit_vert = rSTOP_slit_vert + 1
	  stop_where=4.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Aperture before Q1 (can only check this if next transformation is DRIFT).

	ztmp = 135.064
	zdrift = ztmp - z_exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.12.5222) then
	  rSTOP_Q1_in = rSTOP_Q1_in + 1
	  stop_where=22.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q1 IN  mag bound.

	if (.not.adrift(spectr,1)) write(6,*) 'Transformation #1 is NOT a drift'
	zdrift = driftdist(spectr,1) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  rSTOP_Q1_in = rSTOP_Q1_in + 1
	  stop_where=5.	
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q1.

	call transp(spectr,2,decay_flag,dflag,m2,p,62.75333333d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  rSTOP_Q1_mid = rSTOP_Q1_mid + 1
	  stop_where=6.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q1 OUT mag boundary.

	call transp(spectr,3,decay_flag,dflag,m2,p,31.37666667d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	  rSTOP_Q1_out = rSTOP_Q1_out + 1
	  stop_where=7.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q1, before Q2 (can only check this if next trans. is DRIFT).

	zdrift = 300.464 - 253.16		!Q1 exit is z=253.16
	ztmp = zdrift				!distance from Q1 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.14.9225) then
	  rSTOP_Q1_out = rSTOP_Q1_out + 1
	  stop_where=23.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 314.464 - 300.464
	ztmp = ztmp + zdrift			!distance from Q1 exit.
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.20.9550) then
	  rSTOP_Q2_in = rSTOP_Q2_in + 1
	  stop_where=24.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q2 IN  mag bound.

	if (.not.adrift(spectr,4)) write(6,*) 'Transformation #4 is NOT a drift'
	zdrift = driftdist(spectr,4) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  rSTOP_Q2_in = rSTOP_Q2_in + 1
	  stop_where=8.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q2.

	call transp(spectr,5,decay_flag,dflag,m2,p,121.77333333d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  rSTOP_Q2_mid = rSTOP_Q2_mid + 1
	  stop_where=9.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q2 OUT mag boundary.

	call transp(spectr,6,decay_flag,dflag,m2,p,60.88666667d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	  rSTOP_Q2_out = rSTOP_Q2_out + 1
	  stop_where=10.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q2, before D1 (can only check this if next trans. is DRIFT).

	zdrift = 609.664 - 553.020		!Q2 exit is z=553.02
	ztmp = zdrift				!distance from Q2 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.0073) then
	  rSTOP_Q2_out = rSTOP_Q2_out + 1
	  stop_where=25.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 641.800 - 609.664
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.0073) then
	  rSTOP_Q2_out = rSTOP_Q2_out + 1
	  stop_where=26.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 819.489 - 641.800
	ztmp = ztmp + zdrift			!distance from Q2 exit.
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (abs(xs).gt.50.0 .or. abs(ys).gt.15.0) then
	  rSTOP_D1_in = rSTOP_D1_in + 1
	  stop_where=27.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.

	if (.not.adrift(spectr,7)) write(6,*) 'Transformation #7 is NOT a drift'
	zdrift = driftdist(spectr,7) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	xt=xs
	yt=ys
	call rotate_haxis(-30.0,xt,yt)
	if (abs(xt-2.500).gt.52.5) then		! -50 < x < +55
	  rSTOP_D1_in = rSTOP_D1_in + 1	
	  stop_where=11.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif
	if ( (abs(yt)+0.01861*xt) .gt. 12.5 ) then	!tan(1.066) ~ 0.01861
	  rSTOP_D1_in = rSTOP_D1_in +1
	  stop_where=12.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	call transp(spectr,8,decay_flag,dflag,m2,p,659.73445725d0,pathlen)
	xt=xs
	yt=ys
	call rotate_haxis(30.0,xt,yt)
	if (abs(xt-2.500).gt.52.5) then		! -50 < x < +55
	  rSTOP_D1_out = rSTOP_D1_out + 1
	  stop_where=13.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif

	if ( (abs(yt)+0.01861*xt) .gt. 12.5 ) then	!tan(1.066) ~ 0.01861
	  rSTOP_D1_out = rSTOP_D1_out + 1
	  stop_where=14.
	  x_stop=xt
	  y_stop=yt
	  goto 500
	endif 


! Apertures after D1, before Q3 (can only check this if next trans. is DRIFT).

	zdrift = 1745.33546 - 1655.83446	!D1 exit is z=1655.83446
	ztmp = zdrift				!distance from D1 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.3276) then
	  rSTOP_D1_out = rSTOP_D1_out + 1
	  stop_where=28.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif
	if (abs(xs).gt.50.0 .or. abs(ys).gt.15.0) then
	  rSTOP_D1_out = rSTOP_D1_out + 1
	  stop_where=29.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

	zdrift = 1759.00946 - 1745.33546
	ztmp = ztmp + zdrift			!distance from D1 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (sqrt(xs*xs+ys*ys).gt.30.3276) then
	  rSTOP_Q3_in = rSTOP_Q3_in + 1
	  stop_where=30.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q3 IN  mag bound.

	if (.not.adrift(spectr,9)) write(6,*) 'Transformation #9 is NOT a drift'
	zdrift = driftdist(spectr,9) - ztmp
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  rSTOP_Q3_in = rSTOP_Q3_in + 1	
	  stop_where=15.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Check aperture at 2/3 of Q3.

	call transp(spectr,10,decay_flag,dflag,m2,p,121.7866667d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  rSTOP_Q3_mid = rSTOP_Q3_mid + 1
	  stop_where=16.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Go to Q3 OUT mag boundary.

	call transp(spectr,11,decay_flag,dflag,m2,p,60.89333333d0,pathlen)
	if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	  rSTOP_Q3_out = rSTOP_Q3_out + 1
	  stop_where=17.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Apertures after Q3 (can only check this if next trans. is DRIFT).

	zdrift = 2080.38746 - 1997.76446	!Q3 exit is z=1997.76446
	ztmp = zdrift				!distance from Q3 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (abs(xs).gt.35.56 .or. abs(ys).gt.17.145) then
	  rSTOP_Q3_out = rSTOP_Q3_out + 1
	  stop_where=31.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! Vacuum window is 15.522cm before FP (which is at VDC1)

	zdrift = 2327.47246 - 2080.38746
	ztmp = ztmp + zdrift			!distance from Q3 exit
	call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	if (abs(xs).gt.99.76635 .or. abs(ys).gt.17.145) then
	  rSTOP_Q3_out = rSTOP_Q3_out + 1
	  stop_where=32.
	  x_stop=xs
	  y_stop=ys
	  goto 500
	endif

! If we get this far, the particle is in the hut.

	rSTOP_hut = rSTOP_hut + 1	

	if (.not.adrift(spectr,12)) write(6,*) 'Transformation #12 is NOT a drift'

	zdrift = driftdist(spectr,12) - ztmp	!distance left to go.
	call mc_hrsr_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-zdrift,pathlen)
	if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.
	xs=x_fp
	ys=y_fp 
	dxdzs=dx_fp
	dydzs=dy_fp

! Apply offset to y_fp (detectors offset w.r.t optical axis).
! In ESPACE, the offset is taken out for recon, but NOT for y_fp in ntuple,
! so we do not apply it to ys (which goes to recon), but do shift it for y_fp.
! IF THIS IS TRUE OFFSET, WE SHOULD SHIFT DETECTOR APERTURES - NEED TO CHECK!!!!
! But in general the dectectors don't limit the acceptance, so we should be OK.

	y_fp = y_fp - 0.48		!VDC center is at +4.8mm.

C Reconstruct target quantities.

	call mc_hrsr_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)
	if (.not.ok) then
	  write(6,*) 'mc_hrsr_recon returned ok=',ok
	  goto 500
	endif

C Fill output to return to main code
	dpp = dpp_recon
	dxdz = dph_recon
	dydz = dth_recon
	y = y_recon
	  
	ok_spec = .true.
	rSTOP_successes = rSTOP_successes + 1

C We are done with this event, whether GOOD or BAD.

 500	continue

C ALL done!

	return
	end
