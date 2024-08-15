	subroutine mc_calo (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     >		x_fp, dx_fp, y_fp, dy_fp, m2,
     >		ms_flag, wcs_flag, decay_flag, resmult, frx,fry, ok_spec, 
     >          pathlen,using_tgt_field,
     >          prot_zbeam,spect_side,drift_to_cal)

C+______________________________________________________________________________
!
! Monte-Carlo of Ge/Gm calorimeter.
!	Note that we only pass on real*8 variables to the subroutine.
!	This will make life easier for portability.
!
! Author: David Gaskell
!
! Modification History:
!
!  30-Sep-2002  Just smear the input variables.	
C-______________________________________________________________________________

	implicit 	none

	include '../constants.inc'
        include '../spectrometers.inc'
	include 'struct_calo.inc'

! Spectrometer definitions

	integer*4 spectr
	parameter (spectr = 7)		!calo is spec #7.
	integer*4 spect_side ! 7 (HMS side) or 8 (SOS side)

! Collimator (octagon) dimensions.

	real*8 h_entr,v_entr,h_exit,v_exit



! No collimator, but use collimator dimensions to define calo
c	parameter (h_entr = 64.0)
c	parameter (v_entr = 109.0)
c	parameter (h_exit = 64.0)
c	parameter (v_exit = 109.0)

c Use size of NPS	
	parameter (h_entr = 60.0)
	parameter (v_entr = 72.0)
	parameter (h_exit = 60.0)
	parameter (v_exit = 72.0)

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
	real*8 fry,frx				!vertical position@tgt (+y=down)
	real*8 pathlen
	logical ms_flag				!mult. scattering flag
	logical wcs_flag			!wire chamber smearing flag
	logical decay_flag			!check for particle decay
	logical ok_spec				!true if particle makes it
	logical using_tgt_field

! Local declarations.

	integer*4	chan	/1/,n_classes

	logical	first_time_here/.true./

	real*8 dpp_recon,dth_recon,dph_recon	!reconstructed quantities
	real*8 y_recon
	real*8 ps,p,m2				!More kinematic variables.
	real*8 xt,yt,rt,tht			!temporaries
	real*8 resmult				!DC resolution factor
	real*8 dummy
	real*8 zdrift

	logical dflag			!has particle decayed?
	logical ok

	real*8 ctheta,stheta                      !cos and sin of tgt angle

	real*8 grnd

	real*8 tmpwidth				!used for aperture cut.
	real*8 gauss1
	real*8 prot_zbeam
	real*8 delta_z,delta_y

	real*8 drift_to_cal

	save		!Remember it all!

! ================================ Executable Code =============================
	ctheta = cos(th_spec)	! GAW
	stheta = sin(th_spec)	! GAW

! Initialize ok_spec to .false., reset decay flag

	ok_spec = .false.
	dflag = .false.			!particle has not decayed yet
	caloSTOP_trials = caloSTOP_trials + 1

! Save spectrometer coordinates.

	xs = x
	ys = y
	zs = z
	dxdzs = dxdz
	dydzs = dydz

! particle momentum

	dpps = dpp
	p = p_spec*(1.+dpps/100.)

! Begin transporting particle.

! Do transformations, checking against apertures.
! Check front of fixed slit.  (here slit = calo)

	  zdrift = drift_to_cal
	  call project(xs,ys,zdrift,decay_flag,dflag,m2,p,pathlen)
	  if (abs(ys).gt.h_entr) then
	    caloSTOP_slit_hor = caloSTOP_slit_hor + 1
	    goto 500
	  endif
	  if (abs(xs).gt.v_entr) then
	    caloSTOP_slit_vert = caloSTOP_slit_vert + 1
	    goto 500
	  endif

c positions, angles at front
	  x_fp = xs
	  y_fp = ys
	  dx_fp=dxdzs
	  dy_fp=dydzs

C       skip all this fow, after checking that we hit calorimeter, just fill output
C       with input-maybe add some smearing later.
	  
c no difference between recon angle and "focal plane" angles
c use the z position determined by the hadron arm to
c correct the distance from the beam interaction point to the hit in the calo
c
c	   delta_y = -prot_zbeam*stheta
c	   delta_z =  prot_zbeam*ctheta
c   
c	  dx_fp = (x_fp-fry)/(drift_to_cal-delta_z)
c	  dy_fp = (y_fp-delta_y)/(drift_to_cal-delta_z)

! replace xs,ys,... with 'tracked' quantities.
c	  xs=x_fp
c	  ys=y_fp
c	  dxdzs=dx_fp
c	  dydzs=dy_fp


! Reconstruct target quantities.
c	  call mc_calo_recon(dpp_recon,dth_recon,dph_recon,y_recon,fry,delta_y,delta_z,drift_to_cal)
c
c          if (using_tgt_field) then
c	     ok = .TRUE.
c	     call track_to_tgt(dpp_recon,y_recon,dph_recon,dth_recon,-frx,-fry,
c     >  	  -p,sqrt(m2),ctheta,-stheta,prot_zbeam,-1,ok)
c         endif
! Fill output to return to main code
c	  dpp = dpp_recon
c	  dxdz = dph_recon
c	  dydz = dth_recon
c	  y = y_recon

	  dpp = dpps
	  dxdz = dxdzs
	  dydz = dydzs
	  y = ys	  
	  ok_spec = .true.
	  caloSTOP_successes = caloSTOP_successes + 1

500	continue


	return
	end
