	subroutine mc_shms(p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >			x_fp,dx_fp,y_fp,dy_fp,m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, fry, ok_spec, 
     >         pathlen, spectr)

C+______________________________________________________________________________
!
! Monte-Carlo of SHMS spectrometer.
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
!
!  2007-2009   Updated SHMS  (Horn)
!  2011  Considerable updates (Gaskell,Dutta and Jones)
C-______________________________________________________________________________

	implicit 	none

	include '../shms/struct_shms.inc'
	include '../shms/apertures_shms.inc'
	include '../constants.inc'
	include '../spectrometers.inc'


C Spectrometer definitions - for double arm monte carlo compatability
	integer*4 spectr

C Math constants

	real*8 d_r,r_d,root
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2


! Option for mock sieve slit.  Just take particles and forces trajectory
! to put the event in the nearest hole.  Must choose the "No collimator"
! values for the h(v)_entr and h(v)_exit values.
! Note that this will mess up the physics distributions somewhat, but it
! should still be pretty good for optics. Physics limits (e.g. elastic
! peak at x<=1) will not be preserved.

	logical use_sieve /.false./ !use a fake sieve slit.
c 	logical use_sieve /.true./ !use a fake sieve slit.
c        logical use_coll /.false./ ! use fake collimator
        logical use_coll /.true./ ! use fake collimator
        logical*4 spec_ntuple   /.true./
!        logical skip_hb /.true./
        logical skip_hb /.false./
c        logical use_2sieve /.true./ ! use sieve in front of HB  
        logical use_2sieve /.false./ ! use sieve in front of HB  

C The arguments

	real*8	x,y,z				!(cm)
	real*8	dpp				!delta p/p (%)
	real*8 	dxdz,dydz			!X,Y slope in spectrometer
	real*8	x_fp,y_fp,dx_fp,dy_fp		!Focal plane values to return
	real*8	p_spec,th_spec			!spectrometer setting
        logical ms_flag, wcs_flag		!particle, m_scat, wc_smear
	logical	ok_spec				!true if particle makes it
	logical decay_flag,dflag
	real*8 fry 		!fast raster y position.
	real*8 m2		!mass for call to mc_hms(sos). Changes if decay
	real*8 pathlen,tpathlen
	real*8 resmult		!DC resolution factor
c collimator
        real*8  h_entr,v_entr   !horiz. and vert. 1/2 gap of fixed slit
        real*8  h_exit,v_exit   !horiz. and vert. 1/2 gap of fixed slit
	real*8  y_off,x_off	!horiz. and vert. position offset of slit
	real*8  z_off		!offset in distance from target to front of slit.
	real*8	z_off_lsa,z_off_ssa

	parameter (h_entr = 8.5)
	parameter (v_entr = 12.5)
	parameter (h_exit = 8.65)
	parameter (v_exit = 12.85)
	parameter (x_off=+0.00)
	parameter (y_off=+0.00)

! z-position of important apertures.
	real*8 grnd ! TH - collimator test
	real*8 z_entr,z_thick,z_slit ! TH - collimator test

!	parameter (z_entr = 258.0 + 0.0) ! from target center to Q1 mechanical entrance (from M. Fowler's drawing)

c  lengths of drift in different elements
        real*8 zd_hbin,zd_hbmen,zd_hbmex,zd_hbout
	real*8 zd_q1in,zd_q1men,zd_q1mid,zd_q1mex,zd_q1out
	real*8 zd_q2in,zd_q2men,zd_q2mid,zd_q2mex,zd_q2out
	real*8 zd_q3in,zd_q3men,zd_q3mid,zd_q3mex,zd_q3out
	real*8 zd_q3d1trans,zd_d1flare,zd_d1men
	real*8 zd_d1mid1,zd_d1mid2,zd_d1mid3,zd_d1mid4
	real*8 zd_d1mid5,zd_d1mid6,zd_d1mid7
	real*8 zd_d1mex,zd_d1out
	real*8 zd_fp
c        real*4 spec(58)
        real*8 spec(58)

c	parameter(zd_q1in  = 307.00)
c	parameter(zd_q1mid = 107.00)
c	parameter(zd_q1out = 107.00)
c	parameter(zd_q2in  = 67.00)
c	parameter(zd_q2mid = 92.00)
c	parameter(zd_q2out = 92.00)
c	parameter(zd_q3in  = 81.00) 
c	parameter(zd_q3mid = 92.00)
c	parameter(zd_q3out = 92.00)
c	parameter(zd_q3d1trans  = 70.00 )
c	parameter(zd_d1mid1= 53.714)
c	parameter(zd_d1mid2= 53.714)
c	parameter(zd_d1mid3= 53.714)
c	parameter(zd_d1mid4= 53.714)
c	parameter(zd_d1mid5= 53.714)
c	parameter(zd_d1mid6= 53.714)
c	parameter(zd_d1out = 53.714)
c	parameter(zd_fp    = 327.00)

C Distances from 2011 COSY models
c        parameter(zd_fr_sieve  = 108.0)
c        parameter(zd_hbin  = 118.39)
c	parameter(zd_hbmen = 20.06)
c        parameter(zd_hbmex = 75.134)
c        parameter(zd_hbout = 20.06)
cc	parameter (z_entr = 52.04) ! 6.35 cm in front of Q1
c	parameter (z_entr = 22.4) ! 80 cm from center of HB
c	parameter (z_thick =6.35) !6.35 cm thick
c        parameter(zd_q1in  = 58.39)
c        parameter(zd_q1men = 28.05)
c        parameter(zd_q1mid = 93.95)
c        parameter(zd_q1mex = 93.95)
c        parameter(zd_q1out = 28.05)
c        parameter(zd_q2in  = 25.55)
c        parameter(zd_q2men = 36.45)
c        parameter(zd_q2mid = 82.00)
c        parameter(zd_q2mex = 81.99)
c        parameter(zd_q2out = 36.45)
c        parameter(zd_q3in  = 28.10)
c        parameter(zd_q3men = 36.45)
c        parameter(zd_q3mid = 82.00)
c        parameter(zd_q3mex = 81.99)
c        parameter(zd_q3out = 36.45)
c        parameter(zd_q3d1trans  = 18.00)
c        parameter(zd_d1flare = 30.10)
cc        parameter(zd_d1men = 39.50)
c        parameter(zd_d1men = 39.47)
c        parameter(zd_d1mid1= 36.88)
c        parameter(zd_d1mid2= 36.88)
c        parameter(zd_d1mid3= 36.88)
c        parameter(zd_d1mid4= 36.88)
c        parameter(zd_d1mid5= 36.88)
c        parameter(zd_d1mid6= 36.88)
c        parameter(zd_d1mid7= 36.88)
c        parameter(zd_d1mex = 36.88)
c        parameter(zd_d1out = 60.68)
cc        parameter(zd_fp    = 308.00)
c        parameter(zd_fp    = 307.95)


C Distances for 2017 ME's
        parameter(zd_hbin  = 118.39)
	parameter(zd_hbmen = 17.61) ! shms-2017 ME's
        parameter(zd_hbmex = 80.0)
        parameter(zd_hbout = 17.61)
c	parameter (z_entr = 52.04) ! 6.35 cm in front of Q1
	parameter (z_entr = 25.189) ! 82.789 cm from center of HB
	parameter (z_thick =6.35) !6.35 cm thick
        parameter(zd_q1in  = 58.39)
        parameter(zd_q1men = 28.35) !shms-2017 ME's
        parameter(zd_q1mid = 93.65) !shms-2017 ME's
        parameter(zd_q1mex = 93.65) !shms-2017 ME's
        parameter(zd_q1out = 28.35) !shms-2017 ME's
        parameter(zd_q2in  = 25.55)
        parameter(zd_q2men = 39.1) !shms-2017 ME's
        parameter(zd_q2mid = 79.35) !shms-2017 ME's
        parameter(zd_q2mex = 79.35) !shms-2017 ME's
        parameter(zd_q2out = 39.1) !shms-2017 ME's
        parameter(zd_q3in  = 28.10)
        parameter(zd_q3men = 39.1) !shms-2017 ME's
        parameter(zd_q3mid = 79.35) !shms-2017 ME's
        parameter(zd_q3mex = 79.35) !shms-2017 ME's
        parameter(zd_q3out = 39.1) !shms-2017 ME's
        parameter(zd_q3d1trans  = 18.00)
        parameter(zd_d1flare = 30.10)
        parameter(zd_d1men = 39.47)
        parameter(zd_d1mid1= 36.406263) !shms-2017 ME's
        parameter(zd_d1mid2= 36.406263)
        parameter(zd_d1mid3= 36.406263)
        parameter(zd_d1mid4= 36.406263)
        parameter(zd_d1mid5= 36.406263)
        parameter(zd_d1mid6= 36.406263)
        parameter(zd_d1mid7= 36.406263)
        parameter(zd_d1mex = 36.406263)
        parameter(zd_d1out = 60.68)
        parameter(zd_fp    = 307.95)


C Local declarations.

	integer*4	chan/1/,n_classes

	logical	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p				!More kinematic variables.
	real*8 xt,yt				!temporaries
	real*8 zdrift,x_transp,y_transp,deltacor
	real*8 th_x_fp,th_y_fp,th_dx_fp,th_dy_fp
	real*8 rt,tht ! TH - for sieve slit

	logical	ok_hut

! TH - start sieve 
	logical check_sieve 
	real*8 xlocal,ylocal
	real*8 x_local,y_local
	real*8 check0
	real*8 rad1,rad2,rad
	integer i,j,ij
! TH - end sieve


	save		!Remember it all!

C ================================ Executable Code =============================


! Initialize ok_spec to false
	stop_id = 0
	ok_spec = .false.
	dflag = .false.
	shmsSTOP_trials = shmsSTOP_trials + 1
	xt = th_spec    !avoid 'unused variable' error for th_spec
        if (spec_ntuple) then
         do ij = 1, 58
          spec(ij) = -100000.
         enddo
        endif


!-----------------------------------------------------------------------------
! TH - START SIEVE APERTURE TESTS
! ----------------------------------------------------------------------------
	if (use_2sieve) then
          z_slit=zd_hbin
	  xt = x + z_slit * dxdz	!project to front of HB
	  yt = y + z_slit * dydz
	  xt = 2.5*nint(xt/2.5)	   !shift to nearest hole 
	  yt = 1.64*nint(yt/1.64)
	  rt = 0.3*sqrt(grnd()) ! Tanya's number
	  tht= 2*pi*grnd()		!angle of offset.
	  xt = xt + rt*cos(tht)
	  yt = yt + rt*sin(tht)
	  dxdz = (xt-x)/z_slit		!force to correct angle.
	  dydz = (yt-y)/z_slit
	endif


! Force particles to go through the sieve slit holes, for mock sieve option.

	if (use_sieve) then
          z_slit=zd_hbin+zd_hbmen+zd_hbmex+zd_hbout+z_entr
	  xt = x + z_slit * dxdz	!project to collimator
	  yt = y + z_slit * dydz
	  xt = 2.5*nint(xt/2.5)	   !shift to nearest hole 
!	  xt = 7.98*nint(xt/7.98)	!shift to nearest hole.
!	  yt = 3.99*nint(yt/3.99)
!	  yt = 3.048*nint(yt/3.048)
	  yt = 1.64*nint(yt/1.64)
!          if (yt>0.)then            ! special offset sieve slit
!           yt = 0.82+1.64*nint(yt/1.64)
!          else
!           yt = 1.64*nint(yt/1.64) - 0.82
!          endif
!	  rt = 0.254*sqrt(grnd()) !distance from center of hole(r=2.54mm)
	  rt = 0.3*sqrt(grnd()) ! Tanya's number
!	  rt = 0.166*sqrt(grnd()) !distance from center of hole(r=1.66mm, HB 3mrad)
	  tht= 2*pi*grnd()		!angle of offset.
	  xt = xt + rt*cos(tht)
	  yt = yt + rt*sin(tht)
	  dxdz = (xt-x)/z_slit		!force to correct angle.
	  dydz = (yt-y)/z_slit
	endif


!------------------------------------------------------------------------------

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
	  else
	    write(6,*) 'LSA tune no longer in use!'
	    stop
	  endif
	  call transp_init(spectr,n_classes)
	  close (unit=chan)
	  if (n_classes.ne.32) stop 
     >   'Bender-SHMS, wrong number of transport classes'
	  first_time_here = .false.
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Begin transporting particle.


C Do transformations, checking against apertures.


!----------------------------------------------------------------------------
! TH - START COLLIMATOR APERTURE TESTS
! ---------------------------------------------------------------------------

C Do transformations, checking against apertures.

! Check front of fixed slit
c	  zdrift = z_entr
c	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
c 	   if (abs(ys-y_off).gt.h_entr) then
c	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
c	      goto 500
c	   endif
c	   if (abs(xs-x_off).gt.v_entr) then
c	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
c	      goto 500
c	   endif
c	   if (abs(xs-x_off).gt. (-v_entr/h_entr*abs(ys-y_off)+3*v_entr/2)) then
c      	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
c	      goto 500
c	   endif

! Check back of fixed slit
*	   write (6, *) 'I am at the back of the fixed slit'
!	  zdrift = z_exit-z_entr
!	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen) !project and decay
!	   if (abs(ys-y_off).gt.h_exit) then
!	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
!	      goto 500
!	   endif
!	   if (abs(xs-x_off).gt.v_exit) then
!	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
!	      goto 500
!	   endif
!	   if (abs(xs-x_off).gt. (-v_exit/h_exit*abs(ys-y_off)+3*v_exit/2)) then
!	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
!	      goto 500
!     	   endif
! ---------------------------------------------------------------------------

! restore xs to values at pivot. Needed for sequential transformations
c	   xs = x_transp
c	   ys = y_transp

!The HB apertures are not symmetric and tilted in the spectrometer
!y-direction therefore we need to find the location of the central
!trajectory with respect to the cenrtal axis at the 4 apertures used 
!for the HB model. Note that the central axis is not symmetric at the
! HB entrance, but becomes symmetric at the HB exit.
! 
!DD: Aperture    Y offset (cm)    Tilt angle (deg)
!DD: -------------------------------------------------
!DD:   1 Mech ent     1.51         1.5
!DD:   2 Mag ent      0.98         1.5
!DD:   3 Meg exit     0.98         -1.5
!DD:   4 Mech exit    1.51         -1.5
!DD: -------------------------------------------------
          spec(58) = dpps
          if (skip_hb) then
           zdrift=zd_hbin+zd_hbmen+zd_hbmex+zd_hbout
           xs=xs + zdrift*dxdzs
           ys=ys + zdrift*dydzs
          else
! Go to HB Mech entrance.
	   call transp(spectr,1,decay_flag,dflag,m2,p,
     >     zd_hbin,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_vaxis(1.5,xt,yt)
	   yt = yt + 1.51
           x_hb_in = xt
           y_hb_in = yt
	   if ((xt*xt.gt.r_HBx*r_HBx).or.(yt.gt.r_HBfyp).or.
     >        (yt.lt.r_HBfym)) then
	      shmsSTOP_HB_in = shmsSTOP_HB_in + 1
	      stop_id = 1 
	      goto 500
	   endif
           spec(1)=xt
           spec(2)=yt

! Go to HB Mag entrance
	   call transp(spectr,2,decay_flag,dflag,m2,p,
     >     zd_hbmen,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_vaxis(1.5,xt,yt)
	   yt = yt + 0.98
           x_hb_men = xt
           y_hb_men = yt
	   if ((xt*xt.gt.r_HBx*r_HBx).or.(yt.gt.r_HBmenyp).or.
     >        (yt.lt.r_HBmenym)) then
	      shmsSTOP_HB_men = shmsSTOP_HB_men + 1
	      stop_id = 2
	      goto 500
	   endif
           spec(3)=xt
           spec(4)=yt

! Go to HB  Mag exit
	   call transp(spectr,3,decay_flag,dflag,m2,p,
     >     zd_hbmex,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_vaxis(-1.5,xt,yt)
	   yt = yt + 0.98
           x_hb_mex = xt
           y_hb_mex = yt
	   if ((xt*xt.gt.r_HBx*r_HBx).or.(yt.gt.r_HBmexyp).or.
     >        (yt.lt.r_HBmexym)) then
  	      shmsSTOP_HB_mex = shmsSTOP_HB_mex + 1
	      stop_id = 3
	      goto 500
	   endif
           spec(5)=xt
           spec(6)=yt

! Go to HB exit
	   call transp(spectr,4,decay_flag,dflag,m2,p,
     >     zd_hbout,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_vaxis(-1.5,xt,yt)
	   yt = yt + 1.51
           x_hb_out = xt
           y_hb_out = yt
	   if ((xt*xt.gt.r_HBx*r_HBx).or.(yt.gt.r_HBbyp).or.
     >        (yt.lt.r_HBbym)) then
	      shmsSTOP_HB_out = shmsSTOP_HB_out + 1
	      stop_id = 4 
	      goto 500
	   endif
           spec(7)=xt
           spec(8)=yt
          endif
! simulate collimator -------------------------------------------------
          if (use_coll) then
c           tpathlen=pathlen
  	   zdrift = z_entr
           xt=xs + zdrift*dxdzs
           yt=ys + zdrift*dydzs
c	   call project(xt,yt,zdrift,decay_flag,dflag,m2,p,pathlen) !project 
 	   if (abs(yt-y_off).gt.h_entr) then
	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
	      goto 500
	   endif
	   if (abs(xt-x_off).gt.v_entr) then
	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
	      goto 500
	   endif
	   if (abs(xt-x_off).gt. (-v_entr/h_entr*abs(yt-y_off)+3*v_entr/2)) then
      	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
	      goto 500
	   endif



	   zdrift = z_thick
           xt=xt + zdrift*dxdzs
           yt=yt + zdrift*dydzs
c	   call project(xt,yt,zdrift,decay_flag,dflag,m2,p,pathlen) !project 
	   if (abs(yt-y_off).gt.(h_exit)) then
	      shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
	      goto 500
	   endif
	   if (abs(xt-x_off).gt.(v_exit)) then
	      shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
	      goto 500
	   endif
	   if (abs(xt-x_off).gt. ((-v_exit)/(h_exit)*abs(yt-y_off)+3*(v_exit)/2)) then
	      shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
	      goto 500
     	   endif
           spec(56)=xt
           spec(57)=yt
c           pathlen=tpathlen
          endif
! ---------------------------------------------------------------------------
           
! Go to Q1 IN  Mech entrance
	   call transp(spectr,5,decay_flag,dflag,m2,p,
     >     zd_q1in,pathlen)
           x_q1_in = xs
           y_q1_in = ys
	   if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	      shmsSTOP_Q1_in = shmsSTOP_Q1_in + 1
	      stop_id = 5
	      goto 500
	   endif
           spec(9)=xs
           spec(10)=ys


! Go to Q1 Mag entrance
	   call transp(spectr,6,decay_flag,dflag,m2,p,
     >     zd_q1men,pathlen)
           x_q1_men = xs
           y_q1_men = ys
	   if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	      shmsSTOP_Q1_men = shmsSTOP_Q1_men + 1
	      stop_id = 6
	      goto 500
	   endif

! Check aperture at 1/2 of Q1.
	      call transp(spectr,7,decay_flag,dflag,m2,p,
     >        zd_q1mid,pathlen)
	      x_q1_mid=xs
	      y_q1_mid=ys
	      if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
		 shmsSTOP_Q1_mid = shmsSTOP_Q1_mid + 1
		 stop_id = 7
		 goto 500
	      endif

! Go to Q1 Mag exit
	   call transp(spectr,8,decay_flag,dflag,m2,p,
     >     zd_q1mex,pathlen)
           x_q1_mex = xs
           y_q1_mex = ys
	   if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	      shmsSTOP_Q1_mex = shmsSTOP_Q1_mex + 1
	      stop_id = 8
	      goto 500
	   endif

! Go to Q1 OUT.
	      call transp(spectr,9,decay_flag,dflag,m2,p,
     >        zd_q1out,pathlen)
  	      x_q1_out=xs
	      y_q1_out=ys
	      if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
		 shmsSTOP_Q1_out = shmsSTOP_Q1_out + 1
		 stop_id = 9
		 goto 500
	      endif

! Go to Q2 IN  Mech entrance
	   call transp(spectr,10,decay_flag,dflag,m2,p,
     >     zd_q2in,pathlen)
           x_q2_in = xs
           y_q2_in = ys
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_in = shmsSTOP_Q2_in + 1
	      stop_id = 10
	      goto 500
	   endif
           spec(11)=xs
           spec(12)=ys

! Go to Q2 Mag entrance
	   call transp(spectr,11,decay_flag,dflag,m2,p,
     >     zd_q2men,pathlen)
           x_q2_men = xs
           y_q2_men = ys
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_men = shmsSTOP_Q2_men + 1
	      stop_id = 11
	      goto 500
	   endif

! Check aperture at 1/2 of Q2.
	      call transp(spectr,12,decay_flag,dflag,m2,p,
     >        zd_q2mid,pathlen)
	      x_q2_mid=xs
	      y_q2_mid=ys
	      if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
		 shmsSTOP_Q2_mid = shmsSTOP_Q2_mid + 1
		 stop_id = 12
		 goto 500
	      endif

! Go to Q2 Mag exit
	   call transp(spectr,13,decay_flag,dflag,m2,p,
     >     zd_q2mex,pathlen)
           x_q2_mex = xs
           y_q2_mex = ys
	   if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	      shmsSTOP_Q2_mex = shmsSTOP_Q2_mex + 1
	      stop_id = 13
	      goto 500
	   endif

! Go to Q2 OUT.
	      call transp(spectr,14,decay_flag,dflag,m2,p,
     >        zd_q2out,pathlen)
  	      x_q2_out=xs
	      y_q2_out=ys
	      if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
		 shmsSTOP_Q2_out = shmsSTOP_Q2_out + 1
		 stop_id = 14
		 goto 500
	      endif

! Go to Q3 IN  Mech entrance
	   call transp(spectr,15,decay_flag,dflag,m2,p,
     >     zd_q3in,pathlen)
           x_q3_in = xs
           y_q3_in = ys
	   if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	      shmsSTOP_Q3_in = shmsSTOP_Q3_in + 1
	      stop_id = 15
	      goto 500
	   endif
           spec(13)=xs
           spec(14)=ys

! Go to Q3 Mag entrance
	   call transp(spectr,16,decay_flag,dflag,m2,p,
     >     zd_q3men,pathlen)
           x_q3_men = xs
           y_q3_men = ys
	   if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	      shmsSTOP_Q3_men = shmsSTOP_Q3_men + 1
	      stop_id = 16
	      goto 500
	   endif

! Check aperture at 1/2 of Q3.
	      call transp(spectr,17,decay_flag,dflag,m2,p,
     >        zd_q3mid,pathlen)
	      x_q3_mid=xs
	      y_q3_mid=ys
	      if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
		 shmsSTOP_Q3_mid = shmsSTOP_Q3_mid + 1
		 stop_id = 17
		 goto 500
	      endif

! Go to Q3 Mag exit
	   call transp(spectr,18,decay_flag,dflag,m2,p,
     >     zd_q3mex,pathlen)
           x_q3_mex = xs
           y_q3_mex = ys
	   if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	      shmsSTOP_Q3_mex = shmsSTOP_Q3_mex + 1
	      stop_id = 18
	      goto 500
	   endif

! Go to Q3 OUT.
	      call transp(spectr,19,decay_flag,dflag,m2,p,
     >        zd_q3out,pathlen)
  	      x_q3_out=xs
	      y_q3_out=ys
	      if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
		 shmsSTOP_Q3_out = shmsSTOP_Q3_out + 1
		 stop_id = 19
		 goto 500
	      endif



!This collinear dipole is kind of a complicated beast. Below is a note from
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

!CDG - try to update for new dipole geometry
!DJG: Aperture    X offset (cm)    Tilt angle (deg)
!DJG: -------------------------------------------------
!DJG:   1 Mech ent     0.0          0.0
!DD :   2 Flare        3.5         -9.2
!DJG:   3 Mag ent     -2.82        -9.2
!DJG:   4 DIPOLE      -8.05        -6.9
!DJG:   5 DIPOLE     -11.75        -4.6
!DJG:   6 DIPOLE     -13.96        -2.3
!DJG:   7 DIPOLE     -14.70         0.0
!DJG:   8 DIPOLE     -13.96         2.3
!DJG:   9 DIPOLE     -11.75         4.6
!DJG   10 DIPOLE      -8.05         6.9
!DJG:  11 Meg exit    -2.82         9.2
!DJG:  12 Mech exit    6.88         9.2
!DJG: -------------------------------------------------



! Go to D1 Mech entrance
	      call transp(spectr,20,decay_flag,dflag,m2,p,
     >        zd_q3d1trans,pathlen)
  	      x_d_in=xs
	      y_d_in=ys
	      if ((xs*xs + ys*ys).gt.r_D1*r_D1) then
		 shmsSTOP_D1_in = shmsSTOP_D1_in + 1
		 stop_id = 20
		 goto 500
	      endif
	      spec(15)=xs
	      spec(16)=ys

! Go to D1 entrance flare
! Find intersection with rotated aperture plane.
! Aperture has elliptical form.
	   call transp(spectr,21,decay_flag,dflag,m2,p,
     >     zd_d1flare,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(9.200,xt,yt)
	   xt = xt - 3.5
	   x_d_flr=xt
	   y_d_flr=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_flr = shmsSTOP_D1_flr + 1
	      stop_id =21
	      goto 500
	   endif
	   spec(17)=xt
	   spec(18)=yt

! Go to D1 magnetic entrance
! Find intersection with rotated aperture plane.
! Aperture has elliptical form.
	   call transp(spectr,22,decay_flag,dflag,m2,p,
     >     zd_d1men,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(9.200,xt,yt)
	   xt = xt + 2.82
	   x_d_men=xt
	   y_d_men=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_men = shmsSTOP_D1_men + 1
	      stop_id =22
	      goto 500
	   endif
	   spec(19)=xt
	   spec(20)=yt

! Go to first dipole interior aperture
! Find intersection with rotated aperture plane
! Aperture has eliptical form 
	   call transp(spectr,23,decay_flag,dflag,m2,p,
     >     zd_d1mid1,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(6.9,xt,yt)
	   xt = xt + 8.1
	   x_d_m1=xt
	   y_d_m1=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid1 = shmsSTOP_D1_mid1 + 1
	      stop_id=23
	      goto 500
	   endif
	   spec(21)=xt
	   spec(22)=yt


! Go  to Dipole interior 2
! Find intersection with rotated aperture plane.
	   call transp(spectr,24,decay_flag,dflag,m2,p,
     >     zd_d1mid2,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(4.6,xt,yt)
	   xt = xt + 11.75
	   x_d_m2=xt
	   y_d_m2=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid2 = shmsSTOP_D1_mid2 + 1
	      stop_id=24
	      goto 500
	   endif
	   spec(23)=xt
	   spec(24)=yt
	   
! Go to dipole interior 3
	   call transp(spectr,25,decay_flag,dflag,m2,p,
     >     zd_d1mid3,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(2.3,xt,yt)
	   xt = xt + 13.96
	   x_d_m3=xt
	   y_d_m3=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid3 = shmsSTOP_D1_mid3 + 1
	      stop_id=25
	      goto 500
 	   endif
	   spec(25)=xt
	   spec(26)=yt
	   
	   
! Go to dipole interior 4
	   call transp(spectr,26,decay_flag,dflag,m2,p,
     >     zd_d1mid4,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(0.0,xt,yt)
	   xt = xt + 14.70
	   x_d_m4=xt
	   y_d_m4=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid4 = shmsSTOP_D1_mid4 + 1
	      stop_id=26
	      goto 500
	   endif
	   spec(27)=xt
	   spec(28)=yt
	      
! Go to dipole interior 5
	   call transp(spectr,27,decay_flag,dflag,m2,p,
     >     zd_d1mid5,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-2.3,xt,yt)
	   xt = xt + 13.96
	   x_d_m5=xt
	   y_d_m5=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid5 = shmsSTOP_D1_mid5 + 1
	      stop_id=27
	      goto 500
	   endif
	   spec(29)=xt
	   spec(30)=yt

! Go to diple interior 6
	   call transp(spectr,28,decay_flag,dflag,m2,p,
     >     zd_d1mid6,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-4.6,xt,yt)
	   xt = xt + 11.75
	   x_d_m6=xt
	   y_d_m6=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid6 = shmsSTOP_D1_mid6 + 1
	      stop_id=28
	      goto 500
	   endif
	   spec(31)=xt
	   spec(32)=yt

! Go to diple interior 7
	   call transp(spectr,29,decay_flag,dflag,m2,p,
     >     zd_d1mid7,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-6.9,xt,yt)
	   xt = xt + 8.5
	   x_d_m7=xt
	   y_d_m7=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mid7 = shmsSTOP_D1_mid7 + 1
	      stop_id=29
	      goto 500
	   endif
	   spec(33)=xt
	   spec(34)=yt


! Go to D1 magnetic exit
! Find intersection with rotated aperture plane.
! Aperture has elliptical form.
	   call transp(spectr,30,decay_flag,dflag,m2,p,
     >     zd_d1mex,pathlen)
	   xt=xs
	   yt=ys
	   call rotate_haxis(-9.2,xt,yt)
	   xt = xt + 2.82
	   x_d_mex=xt
	   y_d_mex=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_mex = shmsSTOP_D1_mex + 1
	      stop_id =30
	      goto 500
	   endif
	   spec(35)=xt
	   spec(36)=yt

! Go to D1 OUT mechanical boundary.
! Find intersection with rotated aperture plane.
	   call transp(spectr,31,decay_flag,dflag,m2,p,
     >     zd_d1out,pathlen)
           xt=xs
           yt=ys
           call rotate_haxis(-9.20,xt,yt)
	   xt = xt - 6.9
	   x_d_out=xt
	   y_d_out=yt
	   if ((xt*xt + yt*yt).gt.r_D1*r_D1) then
	      shmsSTOP_D1_out = shmsSTOP_D1_out + 1
	      stop_id=31
	      goto 500
	   endif
           spec(37)=xt
           spec(38)=yt


! Transport to focal plane.
	   call transp(spectr,n_classes,decay_flag,dflag,m2,p,
     >zd_fp,pathlen)


C If we get this far, the particle is in the hut.
*	   write (6,*) 'here I am and the particle has made it into the hut'
	   shmsSTOP_hut = shmsSTOP_hut + 1


! TH test - positions at focal plane before tracking and multiple scattering in detector hut.
	   xs_fp = xs
	   ys_fp = ys
	   dxdzs_fp = dxdzs
	   dydzs_fp = dydzs


C and track through the detector hut
	  call mc_shms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,spec,
     >          ok_hut,0.0,pathlen,spectr)
	   if (.not. ok_hut) then
	      goto 500
	   endif


	   xs = x_fp
	   ys = y_fp
	   dxdzs = dx_fp
	   dydzs = dy_fp
           spec(51)=xs
           spec(52)=ys

C Reconstruct target quantities.

c	   write (6,*) 'I am about to enter the reconstruction file',spectr
	   call mc_shms_recon(dpp_recon,dth_recon,dph_recon,
     > y_recon,fry,spectr)

C Fill output to return to main code

*	   write (6,*) 'I am filling the output to return to main code'
	   dpp = dpp_recon
	   dxdz = dph_recon
           spec(53) = dpp

           dydz = dth_recon
	   y = y_recon
	   spec(54) = dxdz
           spec(55) = dydz
	   ok_spec = .true.
	   shmsSTOP_successes = shmsSTOP_successes + 1

C We are done with this event, whether GOOD or BAD.

500   continue

       
C ALL done!

	end


