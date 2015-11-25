	subroutine mc_sos (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen)

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

	include 'struct_sos.inc'
	include 'apertures_sos.inc'

	include '../constants.inc'
        include '../spectrometers.inc'

! Spectrometer definitions

	integer*4 spectr
	parameter (spectr = 2)		!sos is spec #2.

! Collimator (octagon) dimensions.

	real*8 h_entr,v_entr,h_exit,v_exit

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

! 'large' collimator.
	parameter (h_entr = 7.201)
	parameter (v_entr = 4.696)
	parameter (h_exit = 7.567)
	parameter (v_exit = 4.935)

! z-position of important apertures.            
        real*8 z_entr,z_exit
        parameter (z_entr = 126.3e0)		!nominally 1.263 m
        parameter (z_exit = z_entr + 6.3e0)	!6.3 cm thick

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

	real*8 tmpwidth				!used for aperture cut.

	save		!Remember it all!

! ================================ Executable Code =============================

! Initialize ok_spec to .false., reset decay flag

	ok_spec = .false.
	dflag = .false.			!particle has not decayed yet
	sSTOP_trials = sSTOP_trials + 1
	xt = th_spec    !avoid 'unused variable' error for th_spec

! Force particles to go through the sieve slit holes, for mock sieve option.

	if (use_sieve) then
	  xt = x + z_entr * dxdz	!project to collimator
	  yt = y + z_entr * dydz
	  xt = 1.524*nint(xt/1.524)	!shift to nearest hole.
	  if (abs(yt).lt.1.524) then
	    yt = 1.016*nint(yt/1.016)	!middle 3 holes
	  else
	    yt = yt/abs(yt)*(2.54*nint((abs(yt)-1.016)/2.54)+1.016) !outer holes
	  endif
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
	  if (n_classes.ne.10) stop 'MC_SOS, wrong number of transport classes'
	  first_time_here = .false.
	endif

! Begin transporting particle.

! Do transformations, checking against apertures.
! Check front of fixed slit.

	  zdrift = z_entr
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	  if (abs(ys).gt.h_entr) then
	    sSTOP_slit_hor = sSTOP_slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs).gt.v_entr) then
	    sSTOP_slit_vert = sSTOP_slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs).gt. (-v_entr/h_entr*abs(ys)+3*v_entr/2)) then
	    sSTOP_slit_oct = sSTOP_slit_oct + 1
	    goto 500
	  endif

! Check back of fixed slit, at about 1.325 meter

	  zdrift = z_exit - z_entr
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	  if (abs(ys).gt.h_exit) then
	    sSTOP_slit_hor = sSTOP_slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs).gt.v_exit) then
	    sSTOP_slit_vert = sSTOP_slit_vert + 1
	    goto 500
	  endif
	  if (abs(xs).gt. (-v_exit/h_exit*abs(ys)+3*v_exit/2)) then
	    sSTOP_slit_oct = sSTOP_slit_oct + 1
	    goto 500
	  endif
	  

! Go to quad IN  mag bound. (150 cm from target)

          if (.not.adrift(spectr,1)) write(6,*) 'Transformation #1 is NOT a drift'
          zdrift = driftdist(spectr,1) - z_exit
          call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP_quad_in = sSTOP_quad_in + 1
	    goto 500
	  endif

! Go to quad midplane.

	  call transp(spectr,2,decay_flag,dflag,m2,p,35.0e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP_quad_mid = sSTOP_quad_mid + 1
	    goto 500
	  endif

! Go to quad OUT mag boundary.

	  call transp(spectr,3,decay_flag,dflag,m2,p,35.0e0,pathlen)
	  if ((xs*xs + ys*ys).gt.r2_quad) then
	    sSTOP_quad_out = sSTOP_quad_out + 1
	    goto 500
	  endif

! Go to BM01 IN magnetic boundary.
! Find intersection with rotated aperture plane.
	  call transp(spectr,4,decay_flag,dflag,m2,p,80.0e0,pathlen)

!	  call project(xs,ys,80.0e0,decay_flag,dflag,m2,p,pathlen)
	  xt=xs
	  yt=ys
	  call rotate_haxis(-45.0e0,xt,yt)
	  if ((yt.gt.w_bm01) .or. (-yt.gt.(w_bm01-0.05*2.54)) .or.
     >	       (-xt.gt.t_bm01_in)  .or. (-xt.lt.b_bm01_in)) then
	    sSTOP_bm01_in = sSTOP_bm01_in + 1
	    goto 500
	  endif

! Go to BM01 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(spectr,5,decay_flag,dflag,m2,p,169.52e0,pathlen)
	  xt=xs
	  yt=ys
	  call rotate_haxis(45.0e0,xt,yt)
	  if ((abs(yt).gt.w_bm01) .or.
     >	       (-xt.gt.t_bm01_out)  .or. (-xt.lt.b_bm01_out)) then
	    sSTOP_bm01_out = sSTOP_bm01_out + 1
	    goto 500
	  endif

! Go to BM02 IN magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(spectr,6,decay_flag,dflag,m2,p,80.80e0,pathlen)
	  xt=xs
	  yt=ys
	  call rotate_haxis(49.0e0,xt,yt)
	  if ((abs(yt).gt.w_bm02) .or.
     >	       (-xt.gt.t_bm02_in)  .or. (-xt.lt.b_bm02_in)) then
	    sSTOP_bm02_in = sSTOP_bm02_in + 1
	    goto 500
	  endif

! Go to BM02 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp(spectr,7,decay_flag,dflag,m2,p,77.06e0,pathlen)
	  xt=xs
	  yt=ys
	  call rotate_haxis(57.0e0,xt,yt)
	  if ((abs(yt).gt.w_bm02) .or.
     >	       (-xt.gt.t_bm02_out)  .or. (-xt.lt.b_bm02_out)) then
	    sSTOP_bm02_out = sSTOP_bm02_out + 1
	    goto 500
	  endif

! THIS IS FOR AN OLDER SOS MODEL.
! Go to Exit flange. (Transport to focal plane and drift back).
! There are three apertures that need to be checked:
! - dipole exit at -98.997 cm
! - old exit flange at -55.18 cm (dimensions in apertures_sos.inc)
! - new exit window is at -3.22 cm.  In addition, at -10.84 cm there
! is another constraining aperture.  All are hardwired in for now.
!old	  call transp(spectr,n_classes,decay_flag,dflag,m2,p,pathlen)
!old	  call project(xs,ys,-55.18e0,decay_flag,dflag,m2,p,pathlen)
!old	  call project(xs,ys,43.817e0,decay_flag,dflag,m2,p,pathlen) !from dipole exit to old flange

!less old call transp(spectr,n_classes,decay_flag,dflag,m2,p,98.997e0)
!	  xs=xs-55.18*dxdzs
!	  ys=ys-55.18*dydzs
!	  xt=xs
!	  yt=ys

! Current (may 99) SOS model (these distances come from the distances given
! in the forward_cosy.dat file.  
! BM02_exit to old_exit_window:	43.82 cm
! old_window to new_window:	44.34 cm
! new_window to focal_plane:	10.826267 cm.
! We have aperture checks at the old and new exit windows, and 7.62 cm
! after the new exit window.  Thus we use the MEs to take us to the old
! and new exit windows, and then just project 7.62 cm forward.  The
! remaining (10.826267-7.62)=3.206267 cm is the distance to the focal
! plane at the end, which is passed to mc_sos_hut.

	  call transp(spectr,8,decay_flag,dflag,m2,p,43.82e0,pathlen)	!old flange
	  xt=xs
	  yt=ys
!!!	  call rotate_haxis(57.0e0,xt,yt)
	  call rotate_haxis(45.0e0,xt,yt)
	  if ((abs(yt).gt.w_exit) .or.
     >	       (-xt.gt.t_exit)  .or. (-xt.lt.b_exit)) then
	    sSTOP_exit = sSTOP_exit + 1
	    goto 500
	  endif

! Go to aperture at -10.84 cm (from -55.18 to -10.84 is 44.34 cm)
! Aperture is +/- 18.6944cm wide at top, +/- 10.998 cm wide at bottom,
! and +/- 37.694 cm tall.

!	  call project(xs,ys,44.34e0,decay_flag,dflag,m2,p,pathlen)
!	  xs=xs+44.34*dxdzs
!	  ys=ys+44.34*dydzs

	  call transp(spectr,9,decay_flag,dflag,m2,p,44.34e0,pathlen)	!new aperture.

	  tmpwidth= 10.998 + 0.10209*(xs+37.694)
	  if ((abs(xs).gt.37.694) .or.
     >		abs(ys).gt.tmpwidth) then
	    sSTOP_exit = sSTOP_exit + 1
	    goto 500
	  endif

! Go to window at -3.22 cm (from -10.84 to -3.22 is 7.62 cm)
! Aperture is +/- 12.7cm wide, +/-38.1 cm tall.
!	  call project(xs,ys,7.62e0,decay_flag,dflag,m2,p,pathlen)

	  xs=xs+7.62*dxdzs
	  ys=ys+7.62*dydzs

!DC1 is +/-32cm (offset by ~8cm(at the moment), so this is almost non-cut.
	  if ((abs(xs).gt.38.1) .or.
     >		abs(ys).gt.12.7) then
	    sSTOP_exit = sSTOP_exit + 1
	    goto 500
	  endif

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_sos_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_sos_hut knows
! where to start.  Initial zposition for mc_sos_hut is -3.22 cm so that the
! sum of four drift lengths between pipe and focal plane is 98.997 cm
! (43.817+44.34+7.62+3.22=98.997)

! the final transformation that takes you from the new aperture to the
! focal plane is 10.84 cm (7.62+3.22), BUT IS NOT A FIELD FREE DRIFT!!!
! WE REALLY NEED TO FIX THIS UP (transport to focal plane, do a by-hand
! projection back to the chambers in order to avoid double counting any
! decay).

! If we get this far, the particle is in the hut.

	  sSTOP_hut = sSTOP_hut + 1

! and track through the detector hut

	  call mc_sos_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,-3.206267e0,pathlen)
	  if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.
	  xs=x_fp
	  ys=y_fp
	  dxdzs=dx_fp
	  dydzs=dy_fp

! Reconstruct target quantities.

	  call mc_sos_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)

! Fill output to return to main code
	  dpp = dpp_recon
	  dxdz = dph_recon
	  dydz = dth_recon 
	  y = y_recon
	  
	  ok_spec = .true.
	  sSTOP_successes = sSTOP_successes + 1

500	continue

	return
	end
