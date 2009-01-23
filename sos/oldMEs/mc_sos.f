	subroutine mc_sos (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, tg_rad_len,
     >		m2, ms_flag, wcs_flag, decay_flag, resmult, fry, ok_sos)

C+______________________________________________________________________________
!
! Monte-Carlo of SOS spectrometer.
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

	include '../struct_sos.inc'
	include 'apertures.inc'
	include '../track.inc'
	include '../constants.inc'

! Spectrometer definitions

	integer*4 hms,sos
	parameter (hms = 1)
	parameter (sos = 2)

! Collimator (octagon) dimensions.

	real*8 h_entr,v_entr,h_exit,v_exit

! Open values (no collimator)
c	parameter (h_entr = 12.0)
c	parameter (v_entr = 12.0)
c	parameter (h_exit = 12.0)
c	parameter (v_exit = 12.0)

! 'large' collimator.
	parameter (h_entr = 7.201)
	parameter (v_entr = 4.696)
	parameter (h_exit = 7.567)
	parameter (v_exit = 4.935)

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
	real*8 fry				!vertical position@tgt (+y=up)
	logical ms_flag				!mult. scattering flag
	logical wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical ok_sos				!true if particle makes it

! Local declarations.

	integer*4	chan	/1/,n_classes

	logical	first_time_sos/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon,z_recon
	real*8 p,m2				!More kinematic variables.
	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 cos_ev,sin_ev			!cos and sin of event angle
	real*8 xt,yt				!temporaries
	real*8 resmult				!DC resolution factor
	real*8 dummy

	logical dflag			!has particle decayed?
	logical ok

	real*8 tmpwidth				!used for aperture cut.

	save		!Remember it all!

! ================================ Executable Code =============================

! Initialize ok_sos to .false., reset decay flag

	ok_sos = .false.
	dflag = .false.	!particle has not decayed yet
	sSTOP.trials = sSTOP.trials + 1

! Save spectrometer coordinates.

	xs = x
	ys = y
	zs = z
	dxdzs = dxdz
	dydzs = dydz

! Spectrometer angle already in rad here

	th_spec = th_spec
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)
	cos_ev = cos(th_spec + dydzs)
	sin_ev = sin(th_spec + dydzs)

! particle momentum

	dpps = dpp
	p = p_spec*(1.+dpps/100.)

! Read in transport coefficients.

	if (first_time_sos) then
	  call transp_init(sos,n_classes)
	  close (unit=chan)
	  if (n_classes.ne.8) stop 'MC_SOS, wrong number of transport classes'
	  first_time_sos = .false.
	endif

! Calculate multiple scattering in target
c I've commented this out because it's already taken care of in
c  the SIMULATE7 main code (sub montecarlo() )  --rmm
c	if(ms_flag) call musc(m2,p,tg_rad_len,dydzs,dxdzs)
	dummy = tg_rad_len	!avoid unused variable warning


! Begin transporting particle.
! Do transformations, checking against apertures.

! Check front of fixed slit, at about 1.26 meter
	  call project(xs,ys,126.3d0,decay_flag,dflag,m2,p)
	  if (abs(ys).gt.h_entr) then
	    sSTOP.slit_hor = sSTOP.slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs).gt.v_entr) then
	    sSTOP.slit_vert = sSTOP.slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs).gt. (-v_entr/h_entr*abs(ys)+3*v_entr/2)) then
	    sSTOP.slit_oct = sSTOP.slit_oct + 1
	    goto 500
	  endif

! Check back of fixed slit, at about 1.325 meter

	  call project(xs,ys,6.3d0,decay_flag,dflag,m2,p)
	  if (abs(ys).gt.h_exit) then
	    sSTOP.slit_hor = sSTOP.slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs).gt.v_exit) then
	    sSTOP.slit_vert = sSTOP.slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs).gt. (-v_exit/h_exit*abs(ys)+3*v_exit/2)) then
	    sSTOP.slit_oct = sSTOP.slit_oct + 1
	    goto 500
	  endif
	  

! Go to quad IN  mag bound. (150 cm from target)
! 	   call transp(sos,1,decay_flag,dflag,m2,p)

	  call project(xs,ys,(150.0d0-126.3d0-6.3d0),decay_flag,dflag,m2,p)
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP.quad_in = sSTOP.quad_in + 1
	    goto 500
	  endif

! Go to quad midplane.

	  call transp(sos,2,decay_flag,dflag,m2,p,35.0d0)
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP.quad_mid = sSTOP.quad_mid + 1
	    goto 500
	  endif

! Go to quad OUT mag boundary.

	  call transp(sos,3,decay_flag,dflag,m2,p,35.0d0)
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP.quad_out = sSTOP.quad_out + 1
	    goto 500
	  endif

! Go to BM01 IN magnetic boundary.
! Find intersection with rotated aperture plane.
	  call transp(sos,4,decay_flag,dflag,m2,p,80.0d0)

!	  call project(xs,ys,80.0d0,decay_flag,dflag,m2,p)
	  xt=xs
	  yt=ys
	  call rotate_haxis(-45.0d0,xt,yt)
	  if ((abs(yt).gt.w_bm01) .or.
     >	       (-xt.gt.t_bm01_in)  .or. (-xt.lt.b_bm01_in)) then
	    sSTOP.bm01_in = sSTOP.bm01_in + 1
	    goto 500
	  endif

! Go to BM01 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(sos,5,decay_flag,dflag,m2,p,169.5d0)
	  xt=xs
	  yt=ys
	  call rotate_haxis(45.0d0,xt,yt)
	  if ((abs(yt).gt.w_bm01) .or.
     >	       (-xt.gt.t_bm01_out)  .or. (-xt.lt.b_bm01_out)) then
	    sSTOP.bm01_out = sSTOP.bm01_out + 1
	    goto 500
	  endif

! Go to BM02 IN magnetic boundary.
! Find intersection with rotated aperture plane.
	   call transp(sos,6,decay_flag,dflag,m2,p,80.8d0)

!	  call project(xs,ys,80.8d0,decay_flag,dflag,m2,p)
	  xt=xs
	  yt=ys
	  call rotate_haxis(49.0d0,xt,yt)
	  if ((abs(yt).gt.w_bm02) .or.
     >	       (-xt.gt.t_bm02_in)  .or. (-xt.lt.b_bm02_in)) then
	    sSTOP.bm02_in = sSTOP.bm02_in + 1
	    goto 500
	  endif

! Go to BM02 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(sos,7,decay_flag,dflag,m2,p,77.05d0)
	  xt=xs
	  yt=ys
	  call rotate_haxis(57.0d0,xt,yt)
	  if ((abs(yt).gt.w_bm02) .or.
     >	       (-xt.gt.t_bm02_out)  .or. (-xt.lt.b_bm02_out)) then
	    sSTOP.bm02_out = sSTOP.bm02_out + 1
	    goto 500
	  endif

! Go to Exit flange. (Transport to focal plane and drift back).
! There are three apertures that need to be checked:
! - dipole exit at -98.997 cm
! - old exit flange at -55.18 cm (dimensions in apertures.inc)
! - new exit window is at -3.22 cm.  In addition, at -10.84 cm there
! is another constraining aperture.  All are hardwired in for now.
!	  call transp(sos,n_classes,decay_flag,dflag,m2,p)
!	  call project(xs,ys,-55.18d0,decay_flag,dflag,m2,p)

!	  call project(xs,ys,43.817d0,decay_flag,dflag,m2,p) !from dipole exit to old flange

	  call transp(sos,n_classes,decay_flag,dflag,m2,p,98.997d0)

	  xs=xs-55.18*dxdzs
	  ys=ys-55.18*dydzs

	  xt=xs
	  yt=ys
	  call rotate_haxis(57.0d0,xt,yt)
	  if ((abs(yt).gt.w_exit) .or.
     >	       (-xt.gt.t_exit)  .or. (-xt.lt.b_exit)) then
	    sSTOP.exit = sSTOP.exit + 1
	    goto 500
	  endif

! Go to aperture at -10.84 cm (from -55.18 to -10.84 is 44.34 cm)
! Aperture is +/- 18.6944cm wide at top, +/- 10.998 cm wide at bottom,
! and +/- 37.694 cm tall.

!	  call project(xs,ys,44.34d0,decay_flag,dflag,m2,p)
	  xs=xs+44.34*dxdzs
	  ys=ys+44.34*dydzs

	  tmpwidth= 10.998 + 0.10209*(xs+37.694)
	  if ((abs(xs).gt.37.694) .or.
     >		abs(ys).gt.tmpwidth) then
	    sSTOP.exit = sSTOP.exit + 1
	    goto 500
	  endif

! Go to window at -3.22 cm (from -10.84 to -3.22 is 8.62 cm)
! Aperture is +/- 12.7cm wide, +/-38.1 cm tall.
!	  call project(xs,ys,7.62d0,decay_flag,dflag,m2,p)
	  xs=xs+7.62*dxdzs
	  ys=ys+7.62*dydzs

	  if ((abs(xs).gt.38.1) .or.
     >		abs(ys).gt.12.7) then
	    sSTOP.exit = sSTOP.exit + 1
	    goto 500
	  endif

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_sos_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_hms_hut knows
! where to start.  Initial zposition for mc_hms_hut is -3.22 cm so that the
! sum of four drift lengths between pipe and focal plane is 98.997 cm
! (43.817+44.34+7.62+3.22=98.997)

! If we get this far, the particle is in the hut.

	  sSTOP.hut = sSTOP.hut + 1

! and track through the detector hut

	  call mc_sos_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-3.22d0)
	  if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.
	  xs=x_fp
	  ys=y_fp
	  dxdzs=dx_fp
	  dydzs=dy_fp

! Reconstruct target quantities.

	  call mc_sos_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)
c	   z_recon = -y_recon/sin_ev
c I want ytar, not z!
	  z_recon = y_recon

! Fill output to return to main code
	  dpp = dpp_recon
	  dxdz = dph_recon
	  dydz = dth_recon 
	  y = z_recon
	  
	  ok_sos = 1.

500	continue

	return
	end
