	subroutine mc_hrsl (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, tg_rad_len,
     >		m2, ms_flag, wcs_flag, decay_flag, resmult, fry, ok_hrs)

C+______________________________________________________________________________
!
! Monte-Carlo of HRSL spectrometer.
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

	include 'struct_hrsl.inc'
	include 'apertures.inc'
	include '../track.inc'
	include 'g_dump_all_events.inc'

! Add simulate.inc for drift distances of drift-only transformations, and
! remove constants.inc, because it's included in simulate.
!	include '../constants.inc'
	include '../simulate.inc'

C Spectrometer definitions - for double arm monte carlo compatability

	integer*4 spectr
	parameter (spectr = 4)	!hrs-left is spec #4.

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
	parameter (z_entr = 110.9d0 + z_off)	!nominally 1.109 m
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
	real*8	tg_rad_len			!target length in r.l.
	real*8  fry                         	!vertical position@tgt (+y=down)
	logical	ms_flag				!mult. scattering flag
	logical	wcs_flag			!wire chamber smearing flag
	logical	decay_flag			!check for particle decay
	logical	ok_hrs				!true if particle makes it

C Local declarations.

	integer*4	chan/1/,n_classes

	logical*4	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 p,m2				!More kinematic variables.
	real*8 xt,yt,rt,tht			!temporaries
	real*8 resmult				!DC resolution factor (unused)
	real*8 zdrift				

	logical dflag			!has particle decayed?
	logical	ok

	real*8 grnd

	save		!Remember it all!

C ================================ Executable Code =============================

! Initialize ok_hrs to .flase., restet decay flag

	ok_hrs = .false.
	dflag = .false.			!particle has not decayed yet
	lSTOP.trials = lSTOP.trials + 1

! Force particles to go through the sieve slit holes, for mock sieve option.

	if (use_sieve) then
	  xt = x + z_entr * dxdz	!project to collimator
	  yt = y + z_entr * dydz
	  xt = 2.50*nint(xt/2.50)	!shift to nearest hole.
	  yt = 1.25*nint(yt/1.25)
	  rt = 0.1*sqrt(grnd())		!distance from center of hole(r=1.0mm)
	  tht= 2*pi*grnd()		!angle of offset.
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

C Read in transport coefficients.

	if (first_time_here) then
	  call transp_nonseq_init(spectr,n_classes)
	  close (unit=chan)
	  if (n_classes.ne.12) stop 'MC_HRSL, wrong number of transport classes'
	  first_time_here = .false.
	endif

! calculate multiple scattering in target
c I've commented this out because it's already taken care of in
c  the SIMULATE7 main code (sub montecarlo() )  --rmm
C	if(ms_flag) call musc(m2,p,tg_rad_len,dydzs,dxdzs)

! Begin transporting particle.

C Define values at pivot for transp.f (in HRSL transport coordinates).
C use the xs,ys,dxdzs,dydzs HERE, AFTER mult. scattering in the target.

	x_transp = xs
	y_transp = ys	
	dxdz_transp = dxdzs
	dydz_transp = dydzs

! Do transformations, checking against apertures.
! Check front of fixed slit, at about 1.100 meter for 'hadron arm'
	  call project(xs,ys,z_entr,decay_flag,dflag,m2,p)
	  if (abs(ys-y_off).gt.h_entr) then
	    lSTOP.slit_hor = lSTOP.slit_hor + 1
	    stop_where=1.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_entr) then
	    lSTOP.slit_vert = lSTOP.slit_vert + 1
	    stop_where=2.	
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Check back of fixed slit.

	  zdrift = z_exit-z_entr
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p)
	  if (abs(ys-y_off).gt.h_exit) then
	    lSTOP.slit_hor = lSTOP.slit_hor + 1
	    stop_where=3.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif
	  if (abs(xs-x_off).gt.v_exit) then
	    lSTOP.slit_vert = lSTOP.slit_vert + 1
	    stop_where=4.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to Q1 IN  mag bound.

	  call transp_nonseq(spectr,1,decay_flag,dflag,m2,p,159.03)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    lSTOP.Q1_in = lSTOP.Q1_in + 1
	    stop_where=5.	
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Check aperture at 2/3 of Q1.

	  call transp_nonseq(spectr,2,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    lSTOP.Q1_mid = lSTOP.Q1_mid + 1
	    stop_where=6.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to Q1 OUT mag boundary.

	  call transp_nonseq(spectr,3,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
	    lSTOP.Q1_out = lSTOP.Q1_out + 1
	    stop_where=7.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to Q2 IN  mag bound.

	  call transp_nonseq(spectr,4,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    lSTOP.Q2_in = lSTOP.Q2_in + 1
	    stop_where=8.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Check aperture at 2/3 of Q2.

	  call transp_nonseq(spectr,5,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    lSTOP.Q2_mid = lSTOP.Q2_mid + 1
	    stop_where=9.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to Q2 OUT mag boundary.

	  call transp_nonseq(spectr,6,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
	    lSTOP.Q2_out = lSTOP.Q2_out + 1
	    stop_where=10.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.
! Aperture has elliptical form.

	  call transp_nonseq(spectr,7,decay_flag,dflag,m2,p,99.9)
	  xt=xs
	  yt=ys
	  call rotate_haxis(-30.0,xt,yt)
	  if (abs(xt).gt.dip1) then
	    lSTOP.D1_in = lSTOP.D1_in + 1	
	    stop_where=11.
	    x_stop=xt
	    y_stop=yt
	    goto 500
	  endif
	  if (abs(yt/100).gt.abs((dip2*(1.0-(dip3*(xt/100)/dip4))))) then
	    lSTOP.D1_in = lSTOP.D1_in +1
	    stop_where=12.
	    x_stop=xt
	    y_stop=yt
	    goto 500
	  endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

	  call transp_nonseq(spectr,8,decay_flag,dflag,m2,p,99.9)
	  xt=xs
	  yt=ys
	  call rotate_haxis(30.0,xt,yt)
	  if (abs(xt).gt.dip1) then
	    lSTOP.D1_out = lSTOP.D1_out + 1
	    stop_where=13.
	    x_stop=xt
	    y_stop=yt
	    goto 500
	  endif
	  if (abs(yt/100).gt.abs((dip2*(1.0-(dip3*(xt/100)/dip4))))) then
	    lSTOP.D1_out = lSTOP.D1_out + 1
	    stop_where = 14.
	    x_stop=xt
	    y_stop=yt
	    goto 500
	  endif 

*
* Check a number of apertures in the vacuum pipes following the 
* dipole
*
! Go to Q3 IN  mag bound.

	  call transp_nonseq(spectr,9,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    lSTOP.Q3_in = lSTOP.Q3_in + 1	
	    stop_where=15.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Check aperture at 2/3 of Q3.

	  call transp_nonseq(spectr,10,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    lSTOP.Q3_mid = lSTOP.Q3_mid + 1
	    stop_where=16.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Go to Q3 OUT mag boundary.

	  call transp_nonseq(spectr,11,decay_flag,dflag,m2,p,99.9)
	  if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
	    lSTOP.Q3_out = lSTOP.Q3_out + 1
	    stop_where=17.
	    x_stop=xs
	    y_stop=ys
	    goto 500
	  endif

! Transport to focal plane.
	  call transp_nonseq(spectr,12,decay_flag,dflag,m2,p,99.9)
	  
C If we get this far, the particle is in the hut.
	  lSTOP.hut = lSTOP.hut + 1	
*	  stop_where = 18.

!NEED TO CHECK THE INITIAL Z POSITION RELATIVE TO FP IS IN FACT 0.
	  call mc_hrsl_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok,0.0)
	  if (.not.ok) goto 500

! replace xs,ys,... with 'tracked' quantities.  
	  xs=x_fp
	  ys=y_fp
	  dxdzs=dx_fp
	  dydzs=dy_fp

C Reconstruct target quantities.

	  call mc_hrsl_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry)
	  if (.not.ok) then
	    write(6,*) 'mc_hrsr_recon returned ok=',ok
	    goto 500
	  endif

C Fill output to return to main code
	  dpp = dpp_recon
	  dxdz = dph_recon
           dydz = dth_recon
	  y = y_recon
	  
	  ok_hrs = .true.
	  lSTOP.successes = lSTOP.successes + 1

C We are done with this event, whether GOOD or BAD.

 500	  continue

C ALL done!

	end
