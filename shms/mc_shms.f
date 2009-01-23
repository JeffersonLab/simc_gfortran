	subroutine mc_shms(p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >			x_fp,dx_fp,y_fp,dy_fp,m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, pathlen, spectr, cerflag)

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
!
!  28-APR-2006 Modified for shms with pre-bender magnet (Gaskell)
C-______________________________________________________________________________

	implicit 	none

	include 'struct_shms.inc'
	include 'apertures_shms.inc'
	include '../constants.inc'
	include '../spectrometers.inc'


C Spectrometer definitions - for double arm monte carlo compatability
	integer*4 spectr
C for SHMS, this is passed by SIMC (SSA tune = 5, LSA tune = 6)
C	parameter (spectr=5)  ! SHMS in SSA tune is spec #5

C Math constants

	real*8 d_r,r_d,root
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2


C The arguments

	real*8		x,y,z				!(cm)
	real*8		dpp				!delta p/p (%)
	real*8 		dxdz,dydz			!X,Y slope in spectrometer
	real*8		x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8		p_spec,th_spec			!spectrometer setting
        logical         ms_flag, wcs_flag		!particle, m_scat, wc_smear
	logical		ok_spec				!true if particle makes it
	logical         decay_flag,dflag
	logical		cerflag				!true if use first cerenkov
	real*8 fry	!fast raster y position.
	real*8 m2	!mass for call to mc_hms(sos). Changes if decay
	real*8 pathlen
	real*8 resmult				!DC resolution factor
c collimator
        real*8  h_entr,v_entr   !horiz. and vert. 1/2 gap of fixed slit
        real*8  h_exit,v_exit   !horiz. and vert. 1/2 gap of fixed slit
	real*8  y_off,x_off	!horiz. and vert. position offset of slit
	real*8  z_off		!offset in distance from target to front of slit.
	real*8	z_off_lsa,z_off_ssa
	parameter (h_entr = 43.536)
	parameter (v_entr = 49.002)
	parameter (h_exit = 43.708)
	parameter (v_exit = 49.444)
	parameter (x_off=+0.00)
	parameter (y_off=+0.00)
	parameter (z_off_ssa=+0.00)	! for SSA should be 0.0
	parameter (z_off_lsa=-232.0)	! for LSA should be -232.0

! z-position of important apertures.
c	real*8 z_entr,z_exit
c	Calculate these later, since z_off is tune-dependent
c	parameter (z_entr = 400.+ z_off)	!
c	parameter (z_exit = z_entr + 10.0d0)	!10.0 cm thick

c  lengths of drift in different elements
	real*8 zd_q1in,zd_q1mid,zd_q1out
	real*8 zd_q2in,zd_q2mid,zd_q2out
	real*8 zd_q3in,zd_q3mid
	real*8 zd_q3d1trans
	real*8 zd_d1mid1,zd_d1mid2,zd_d1mid3,zd_d1mid4,zd_d1mid5,zd_d1mid6
	real*8 zd_d1out
	real*8 zd_magout
cdg	real*8 zd_bpin,zd_bpout
	real*8 zd_fp

	parameter(zd_q1in  = 232.575)
	parameter(zd_q1mid = 106.175)
	parameter(zd_q1out = 106.175)
	parameter(zd_q2in  = 62.15)
	parameter(zd_q2mid = 106.175)
	parameter(zd_q2out = 106.175)
	parameter(zd_q3in  = 123.575) 
	parameter(zd_q3mid= 100.001)
	parameter(zd_q3d1trans  = 100.00)
	parameter(zd_d1mid1= 50.597)
	parameter(zd_d1mid2= 50.353)
	parameter(zd_d1mid3= 50.10)
	parameter(zd_d1mid4= 50.007)
	parameter(zd_d1mid5= 50.1)
	parameter(zd_d1mid6= 50.356)
	parameter(zd_d1out = 50.595)
	parameter(zd_magout = 50.648)
	parameter(zd_fp    = 275.92)

c  lengths of drift in different elements - SSA tune
cdg	real*8 zd_ssa_q1in,zd_ssa_q1mid,zd_ssa_q1out
cdg	real*8 zd_ssa_q2in,zd_ssa_q2mid,zd_ssa_q2out
cdg	real*8 zd_ssa_d1in,zd_ssa_d1mid,zd_ssa_d1out
cdg	real*8 zd_ssa_bpin,zd_ssa_bpout,zd_ssa_fp
cdg	parameter (zd_ssa_q1in  = 464.575)	!
cdg	parameter (zd_ssa_q1mid = 0.0)		!
cdg	parameter (zd_ssa_q1out = 0.0)		!
cdg	parameter (zd_ssa_q2in  = 111.15)	!
cdg	parameter (zd_ssa_q2mid = 10.0)		!
cdg	parameter (zd_ssa_q2out = 10.0)		!
cdg	parameter (zd_ssa_bpin  = 86.842)	!
cdg	parameter (zd_ssa_d1in  = 10.0)		!
cdg	parameter (zd_ssa_d1out = 10.0)		!
cdg	parameter (zd_ssa_bpout = 10.0)		!
cccdg	parameter (zd_ssa_fp    = 308.267)	!

c  lengths of drift in different elements - LSA tune
cdg	real*8 zd_lsa_q1in,zd_lsa_q1mid,zd_lsa_q1out
cdg	real*8 zd_lsa_q2in,zd_lsa_q2mid,zd_lsa_q2out
cdg	real*8 zd_lsa_d1in,zd_lsa_d1mid,zd_lsa_d1out
cdg	real*8 zd_lsa_bpin,zd_lsa_bpout,zd_lsa_fp
cdg	parameter (zd_lsa_q1in  = 232.575)	!
cdg	parameter (zd_lsa_q1mid = 10.0)		!
cdg	parameter (zd_lsa_q1out = 10.0)		!
cdg	parameter (zd_lsa_q2in  = 111.15)	!
cdg	parameter (zd_lsa_q2mid = 10.0)		!
cdg	parameter (zd_lsa_q2out = 10.0)		!
cdg	parameter (zd_lsa_bpin  = 61.519)	!
cdg	parameter (zd_lsa_d1in  = 10.0)		!
cdg	parameter (zd_lsa_d1out = 10.0)		!
cdg	parameter (zd_lsa_bpout = 10.0)		!
cdg	parameter (zd_lsa_fp    = 382.94)	!
c
C Local declarations.

	integer*4	chan/1/,n_classes

	logical	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p				!More kinematic variables.
	real*8 xt,yt				!temporaries
	real*8 zdrift,x_transp,y_transp

	logical	ok

	save		!Remember it all!

C ================================ Executable Code =============================


! Initialize ok_spec to false
	ok_spec = .false.
	dflag = .false.
	shmsSTOP_trials = shmsSTOP_trials + 1
	xt = th_spec    !avoid 'unused variable' error for th_spec

! particle momentum


	dpps = dpp
	p = p_spec*(1.+dpps/100.)

! and the rest....
	xs = x
	ys = y
	zs = z
	dxdzs = dxdz
	dydzs = dydz
	x_transp =xs
	y_transp =ys


C Read in transport coefficients.  Choose tune-dependent parameters.
	if (first_time_here) then
	  if (spectr.eq.5) then
	     z_off = 0.0		    
cdg	    z_off = z_off_ssa
cdg	    zd_q1in  = zd_ssa_q1in
cdg	    zd_q1mid = zd_ssa_q1mid
cdg	    zd_q1out = zd_ssa_q1out
cdg	    zd_q2in  = zd_ssa_q2in
cdg	    zd_q2mid = zd_ssa_q2mid
cdg	    zd_q2out = zd_ssa_q2out
cdg	    zd_d1in  = zd_ssa_d1in
cdg	    zd_d1mid = zd_ssa_d1mid
cdg	    zd_d1out = zd_ssa_d1out
cdg	    zd_bpin  = zd_ssa_bpin
cdg	    zd_bpout = zd_ssa_bpout
cdg	    zd_fp    = zd_ssa_fp
cdg	    zd_d1mid = zd_ssa_d1mid
cdg	    zd_bpin  = zd_ssa_bpin
cdg	    zd_bpout = zd_ssa_bpout
	  else
	    write(6,*) 'LSA tune no longer in use!'
	    stop
cdg	    z_off = z_off_lsa
cdg	    zd_q1in  = zd_lsa_q1in
cdg	    zd_q1mid = zd_lsa_q1mid
cdg	    zd_q1out = zd_lsa_q1out
cdg	    zd_q2in  = zd_lsa_q2in
cdg	    zd_q2mid = zd_lsa_q2mid
cdg	    zd_q2out = zd_lsa_q2out
cdg	    zd_d1in  = zd_lsa_d1in
cdg	    zd_d1mid = zd_lsa_d1mid
cdg	    zd_d1out = zd_lsa_d1out
cdg	    zd_bpin  = zd_lsa_bpin
cdg	    zd_bpout = zd_lsa_bpout
cdg	    zd_fp    = zd_lsa_fp
	  endif
cdg	  z_entr = 400.+ z_off		!Front of fixed slit.
cdg	  z_exit = z_entr + 10.0d0	!10.0 cm thick.
	  call transp_init(spectr,n_classes)
	  close (unit=chan)
cdg	  if (n_classes.ne.12) stop 'SHMS, wrong number of transport classes'
	  if (n_classes.ne.18) stop 'Bender-SHMS, wrong number of transport classes'
	  first_time_here = .false.
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Begin transporting particle.


C Do transformations, checking against apertures.

Cdg No collimator for now. Add later

! Check front of fixed slit, at about 4.0 meter
cdg	  zdrift = z_entr
cdg	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
cdg	   if (abs(ys-y_off).gt.h_entr) then
cdg	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
cdg	      goto 500
cdg	   endif
cdg	   if (abs(xs-x_off).gt.v_entr) then
cdg	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
cdg	      goto 500
cdg	   endif
cdg	   if (abs(xs-x_off).gt. (-v_entr/h_entr*abs(ys-y_off)+3*v_entr/2)) then
cdg	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
cdg	      goto 500
cdg	   endif

! Check back of fixed slit, at about 4.010 meter
*	   write (6, *) 'I am at the back of the fixed slit'
cdg	  zdrift = z_exit-z_entr
cdg	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
cdg	   if (abs(ys-y_off).gt.h_exit) then
cdg	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
cdg	      goto 500
cdg	   endif
cdg	   if (abs(xs-x_off).gt.v_exit) then
cdg	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
cdg	      goto 500
cdg	   endif
cdg	   if (abs(xs-x_off).gt. (-v_exit/h_exit*abs(ys-y_off)+3*v_exit/2)) then
cdg	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
cdg	      goto 500
cdg	   endif

! restore xs to values at pivot. Needed for sequential transformations
	   xs = x_transp
	   ys = y_transp

C leave it this way for all the matrix transformations

! Go to Q1 IN  mag bound.
! Note that this transformation includes the splitter
c	   write (6, *) 'I am at Q1'
	   call transp(spectr,1,decay_flag,dflag,m2,p,zd_q1in,pathlen)
	   if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	      shmsSTOP_Q1_in = shmsSTOP_Q1_in + 1
	      goto 500
	   endif

! Check aperture at 1/2 of Q1.
c	   write (6,*) 'I am 1/2 of the way through Q1'
	      call transp(spectr,2,decay_flag,dflag,m2,p,zd_q1mid,pathlen)
	      if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
		 shmsSTOP_Q1_mid = shmsSTOP_Q1_mid + 1
		 goto 500
	      endif

! Go to Q1 OUT mag boundary.
*	   write (6,*) 'leaving Q1'
	      call transp(spectr,3,decay_flag,dflag,m2,p,zd_q1out,pathlen)
	      if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
		 shmsSTOP_Q1_out = shmsSTOP_Q1_out + 1
		 goto 500
	      endif

! Go to Q2 IN  mag bound.
*	      write (6,*) 'entering Q2'
	   call transp(spectr,4,decay_flag,dflag,m2,p,zd_q2in,pathlen)
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_in = shmsSTOP_Q2_in + 1
	      goto 500
	   endif

! Check aperture at 1/2 of Q2.
	   call transp(spectr,5,decay_flag,dflag,m2,p,zd_q2mid,pathlen)
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_mid = shmsSTOP_Q2_mid + 1
	      goto 500
	   endif

! Go to Q2 OUT mag boundary.
	   call transp(spectr,6,decay_flag,dflag,m2,p,zd_q2out,pathlen)
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_out = shmsSTOP_Q2_out + 1
	      goto 500
	   endif

! Go to Q3 IN  mag bound.
! Aperture is a 30 cm pipe - same as d1 in old shms
*	   write (6,*) 'entering Q3'
	   call transp(spectr,7,decay_flag,dflag,m2,p,zd_q3in,pathlen)
	   if ((xs*xs + ys*ys).gt.r_D1*r_D1) then
	      shmsSTOP_Q3_in = shmsSTOP_Q3_in + 1
	      goto 500
	   endif


! Check aperture at 1/2 of Q3.
	   call transp(spectr,8,decay_flag,dflag,m2,p,zd_q3mid,pathlen)
	   if ((xs*xs + ys*ys).gt.r_D1*r_D1) then
	      shmsSTOP_Q3_mid = shmsSTOP_Q3_mid + 1
	      goto 500
	   endif


!This colinear dipole is kind of a complicated beast. Below is a note from
!Dave Potterveld (DHP) describing the relevant apertures

!DHP: SHMS colinear dipole apertures
!DHP:  
!DHP: Each aperture is a tilted offset circle of 30 cm radius.
!DHP: The X offset is the X coordinate of
!DHP: the center of the circle in the standard SHMS coordinates, with X
!DHP: orthogonal to the (curved) central ray, and pointing towards the floor.
!DHP: Z is the direction of the central ray, and Y, parallel to the lab floor,
!DHP: forms a right handed coordinate system.
!DHP: 
!DHP: The tilt angle is the angle the circle is rotated about the Y axis.
!DHP: Positive angles mean you would stub your toe before hitting your head.
!DHP: 
!DHP: Aperture    X offset (cm)    Tilt angle (deg)
!DHP: -------------------------------------------------
!DHP:   1 Q3/DIP       8.3626423  -9.1321001
!DHP:   2 DIPOLE       0.7486990  -8.0823002
!DHP:   3 DIPOLE      -5.1099238  -5.1431999
!DHP:   4 DIPOLE      -8.1272335  -1.7112000
!DHP:   5 DIPOLE      -8.1141739   1.7406000
!DHP:   6 DIPOLE      -5.0727034   5.1682000
!DHP:   7 DIPOLE       0.8021493   8.0931997
!DHP:   8 Yoke exit    8.4129658   9.1112003
!DHP: -------------------------------------------------

! Go to Q3-Dipole transition.
	   call transp(spectr,9,decay_flag,dflag,m2,p,zd_q3d1trans,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(9.1321001,xt,yt)
	   xt = xt - 8.3626423
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_in = shmsSTOP_D1_in + 1
	      goto 500
	   endif

! Go to D1 IN mechanical boundary. Kind of.
! Find intersection with rotated aperture plane.
! Aperture has elliptical form.

	   call transp(spectr,10,decay_flag,dflag,m2,p,zd_d1mid1,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(8.0823002,xt,yt)
	   xt = xt - 0.7486990
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid1 = shmsSTOP_D1_mid1 + 1
	      goto 500
	   endif

! Go  to Dipole interior 2
! Find intersection with rotated aperture plane.
	   call transp(spectr,11,decay_flag,dflag,m2,p,zd_d1mid2,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(5.1431999,xt,yt)
	   xt = xt +5.1099238
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid2 = shmsSTOP_D1_mid2 + 1
	      goto 500
	   endif
	   
! Go to dipole interior 3
	   call transp(spectr,12,decay_flag,dflag,m2,p,zd_d1mid3,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(1.7112,xt,yt)
	   xt = xt +8.1272335
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid3 = shmsSTOP_D1_mid3 + 1
	      goto 500
	   endif
	   
	   
! Go to dipole interior 4
	   call transp(spectr,13,decay_flag,dflag,m2,p,zd_d1mid4,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-1.7406,xt,yt)
	   xt = xt +8.1141739
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid4 = shmsSTOP_D1_mid4 + 1
	      goto 500
	   endif
	      
! Go to dipole interior 5
	   call transp(spectr,14,decay_flag,dflag,m2,p,zd_d1mid5,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-5.16820,xt,yt)
	   xt = xt +5.072703
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid5 = shmsSTOP_D1_mid5 + 1
	      goto 500
	   endif

! Go to diple interior 6
	   call transp(spectr,15,decay_flag,dflag,m2,p,zd_d1mid6,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-8.0931997,xt,yt)
	   xt = xt - 0.8021493
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid6 = shmsSTOP_D1_mid6 + 1
	      goto 500
	   endif

! Go to D1 OUT mechanical boundary.
! Find intersection with rotated aperture plane.
	   call transp(spectr,16,decay_flag,dflag,m2,p,zd_d1out,pathlen)
           xt=xs
           yt=ys
           call rotate_haxis(-9.1112003,xt,yt)
	   xt = xt - 8.4129658
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	   shmsSTOP_D1_out = shmsSTOP_D1_out + 1
	      goto 500
	   endif

! Go to end of field map
	   call transp(spectr,17,decay_flag,dflag,m2,p,zd_magout,pathlen)

! Transport to focal plane.
	   call transp(spectr,n_classes,decay_flag,dflag,m2,p,zd_fp,pathlen)



cdg! For the LSA tune, the focal plane is 1m sooner, so project 1m forward
cdg! so that we start at same place w.r.t. detectors in mc_shms_hut.f
cdg	   if (spectr.eq.6) then	!SSA tune
cdg	     zdrift = 100.
cdg	     call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
cdg	   endif

C For this tune, the focal plane is approximately at the 6 cm beyond the first drift chamber
C This means we need to project ~34 cm further to be at "z=0" for mc_shms_hut
C i.e., halfway between each drift chamber
	   zdrift = 33.36
	   call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay

C If we get this far, the particle is in the hut.
*	   write (6,*) 'here I am and the particle has made it into the hut'
	   shmsSTOP_hut = shmsSTOP_hut + 1

C and track through the detector hut
	  call mc_shms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,0.0,pathlen,cerflag)
	   if (.not. ok) goto 500

C for mc_shms_hut, the 'focal' plane is for the SSA tune.  LSA is 100cm forward.

cdg	   if (spectr.eq.6) then
cdg	     x_fp = x_fp - 100.*dxdzs
cdg	     y_fp = y_fp - 100.*dydzs
cdg	   endif

C Need to project from "xfp,yfp" at drift chambers back to real focal plane

	   x_fp = x_fp - 33.36*dx_fp
	   y_fp = y_fp - 33.36*dy_fp


	   xs = x_fp
	   ys = y_fp
	   dxdzs = dx_fp
	   dydzs = dy_fp

C Reconstruct target quantities.
c	   write (6,*) 'I am about to enter the reconstruction file',spectr
	   call mc_shms_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry,spectr)

C Fill output to return to main code

*	   write (6,*) 'I am filling the output to return to main code'
	   dpp = dpp_recon
	   dxdz = dph_recon
C Apply by-hand correction to dth_recon
	   dth_recon = dth_recon + ( -0.25176D-04 +0.52606D-03*dpps + 
     >         (-0.53025D-05)*dpps**2 + 0.39190E-07*dpps**3)
           dydz = dth_recon
	   y = y_recon

	   ok_spec = .true.
	   shmsSTOP_successes = shmsSTOP_successes + 1

C We are done with this event, whether GOOD or BAD.

500   continue

C ALL done!

	end
