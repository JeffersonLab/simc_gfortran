	subroutine mc_hms_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok_hut,zinit,pathlen)

C----------------------------------------------------------------------
C
C Monte-Carlo of HMS detector hut.
C	Note that we only pass on real*8 variables to the subroutine.
C	This will make life easier for portability.
C
C 	The particle is stepped through the detector (using project), and
C	multiple scattering is applied for each detector or air gap.
C	If particle decay in enabled, then project.f also checks for
C	decay of particle.  The particle starts at z=zinit.  This
C	needs to be before the first mult. scattering point (the exit window)
C	or the decay distance is negative, and things are BAD.
C
C----------------------------------------------------------------------

	implicit none

	include 'struct_hms.inc'
	include '../spectrometers.inc'

C Math constants

	real*8 pi,d_r,r_d,root
	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2

	real*8 grnd,gauss1			!external functions

C----------------------------------------------------------------------
C HMS_MATERIALS
C CTP parameter file containing the materials of all the HMS detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C----------------------------------------------------------------------

CDG C spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
c	real*8 hfoil_exit_radlen,hfoil_exit_thick
c	parameter (hfoil_exit_radlen = 53.3)
c	parameter (hfoil_exit_thick  = 0.020*2.54)

C DJG The exit window was changed to titanium sometime in around 2000. 
C DJG Brindza says it's 0.02 inches thick.
C spectrometer exit window, 0.02in=20mil Titanium (X0=3.56cm)
c	real*8 hfoil_exit_radlen,hfoil_exit_thick
c	parameter (hfoil_exit_radlen = 3.56)
c	parameter (hfoil_exit_thick  = 0.020*2.54)
c
c from log entry entry 122424 Sept 21, 2006
c   Ti HMS exit window replaces by preformed aluminum window 0.011 inches thick.
	real*8 hfoil_exit_radlen,hfoil_exit_thick
	parameter (hfoil_exit_radlen = 8.90)
	parameter (hfoil_exit_thick  = 0.011*2.54)

C spectrometer air gaps
	real*8 hair_radlen
	parameter (hair_radlen = 30420.)

C chamber entrance foil, 1 mil of Mylar
	real*8 hdc_entr_radlen,hdc_entr_thick
	parameter (hdc_entr_radlen = 28.7)
	parameter (hdc_entr_thick  = 0.001*2.54)

C chamber gas, 50/50 ethane/argon
	real*8 hdc_radlen,hdc_thick
	parameter (hdc_radlen = 16700.0)
	parameter (hdc_thick  = 1.8)

C effective wire plane, 25 micron W+Be/Cu gives <t>=pi*(.0025/2)**2
	real*8 hdc_wire_radlen,hdc_wire_thick
	parameter (hdc_wire_radlen = 0.35)	!Assuming all Tungsten
	parameter (hdc_wire_thick  = 0.0000049)

C effective cathode plane, Be/Cu
	real*8 hdc_cath_radlen,hdc_cath_thick
	parameter (hdc_cath_radlen = 7.2)	!'Ave' of Be/Cu
	parameter (hdc_cath_thick  = 0.000177)

C chamber exit foil, 1 mil of Mylar
	real*8 hdc_exit_radlen,hdc_exit_thick
	parameter (hdc_exit_radlen = 28.7)
	parameter (hdc_exit_thick  = 0.001*2.54)

C Aerogel entrance foil
	real*8 haer_entr_radlen,haer_entr_thick
	parameter (haer_entr_radlen = 8.90)
	parameter (haer_entr_thick  = 0.15)

C Aerogel aerogel
	real*8 haer_radlen,haer_thick
	parameter (haer_radlen = 150.0)
	parameter (haer_thick  = 9.0)

C Aerogel air
	real*8 haer_air_radlen,haer_air_thick
	parameter (haer_air_radlen = 30420.0)
	parameter (haer_air_thick  = 16.0)

C Aerogel exit foil
	real*8 haer_exit_radlen,haer_exit_thick
	parameter (haer_exit_radlen = 8.90)
	parameter (haer_exit_thick  = 0.1)

C hodoscopes
	real*8 hscin_radlen
	parameter (hscin_radlen =  42.4)

C Cherenkov entrance foil, 40 mil of Al
	real*8 hcer_entr_radlen,hcer_entr_thick
	parameter (hcer_entr_radlen = 8.90)
	parameter (hcer_entr_thick  = 0.040*2.54)

C Cerenkov gas.
	real*8 hcer_radlen
C 0.5 atm of CO2 for e/pi
C	parameter (hcer_radlen = 36620.0)
C 0.5 atm of Freon for better e/pi
	parameter (hcer_radlen = 9620.0)
C 2 atm of Freon for pi/p
C	parameter (hcer_radlen = 2405.0)

C Cherenkov mirror, 75 mu plus 2 cm of Rotacell +
	real*8 hcer_mir_radlen,hcer_mir_thick
	parameter (hcer_mir_radlen = 400.0)
	parameter (hcer_mir_thick  = 2.0)

C Cherenkov exit foil
	real*8 hcer_exit_radlen,hcer_exit_thick
	parameter (hcer_exit_radlen = 8.90)
	parameter (hcer_exit_thick  = 0.040*2.54)

C shower counter
	real*8 hcal_radlen
	parameter (hcal_radlen = 2.50)

C Wire chamber resolutions (sigma)
	real*8 hdc_sigma(1:12)/	0.030,0.030,0.030,0.030,0.030,0.030,
     >				0.030,0.030,0.030,0.030,0.030,0.030/

C Wire plane positions, construct hdc_zpos array using these parameters
	integer*4 hdc_nr_cham,hdc_nr_plan
	parameter (hdc_nr_cham = 2)
	parameter (hdc_nr_plan = 6)

	real*8 hdc_1_zpos,hdc_1_left,hdc_1_right,hdc_1_top,hdc_1_bot
	real*8 hdc_1x_offset,hdc_1y_offset
	real*8 hdc_2_zpos,hdc_2_left,hdc_2_right,hdc_2_top,hdc_2_bot
	real*8 hdc_2x_offset,hdc_2y_offset
	real*8 hdc_del_plane

	parameter (hdc_1_zpos = -51.923)
	parameter (hdc_2_zpos =  29.299)
	parameter (hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick)
	parameter (hdc_1_left  =  26.0)
	parameter (hdc_1_right = -26.0)
	parameter (hdc_1y_offset = 1.443)     ! changed to match hdc.pos.nov95
	parameter (hdc_1_top   = -56.5)
	parameter (hdc_1_bot   =  56.5)
	parameter (hdc_1x_offset = 1.670)
	parameter (hdc_2_left  =  26.0)
	parameter (hdc_2_right = -26.0)
	parameter (hdc_2y_offset = 2.753)     !changed as above,on 03/03,DD
	parameter (hdc_2_top   = -56.5)
	parameter (hdc_2_bot   =  56.5)
	parameter (hdc_2x_offset = 2.758)

C Aerogel positions and thicknesses
	real*8 haer_zentrance,haer_zexit
	parameter (haer_zentrance = 35.699)    !asumming 29.299+1.8*3+1.0 
	parameter (haer_zexit     = 60.949)

C Scintillator positions and thicknesses
	real*8 hscin_1x_zpos,hscin_1y_zpos
	real*8 hscin_1x_thick,hscin_1y_thick
	real*8 hscin_1x_left,hscin_1x_right,hscin_1x_offset
	real*8 hscin_1y_top,hscin_1y_bot,hscin_1y_offset
	real*8 hscin_2x_zpos,hscin_2y_zpos
	real*8 hscin_2x_thick,hscin_2y_thick
	real*8 hscin_2x_left,hscin_2x_right,hscin_2x_offset
	real*8 hscin_2y_top,hscin_2y_bot,hscin_2y_offset

	parameter (hscin_1x_zpos =  77.830)
	parameter (hscin_1y_zpos =  97.520)
	parameter (hscin_2x_zpos = 298.820)
	parameter (hscin_2y_zpos = 318.510)
	parameter (hscin_1x_thick = 1.067)	!1 cm thick, 6.7% overlap
	parameter (hscin_1y_thick = 1.067)
	parameter (hscin_2x_thick = 1.067)
	parameter (hscin_2y_thick = 1.067)
	parameter (hscin_1x_left  =  37.75)
	parameter (hscin_1x_right = -37.75)
	parameter (hscin_1x_offset = -1.55)
	parameter (hscin_1y_top   = -60.25)
	parameter (hscin_1y_bot   =  60.25)
	parameter (hscin_1y_offset = -0.37)
	parameter (hscin_2x_left  =  37.75)
	parameter (hscin_2x_right = -37.75)
	parameter (hscin_2x_offset = -0.63)
	parameter (hscin_2y_top   = -60.25)
	parameter (hscin_2y_bot   =  60.25)
	parameter (hscin_2y_offset = -1.46)

C Cherenkov position
	real*8 hcer_zentrance,hcer_zmirror,hcer_zexit
	parameter (hcer_zentrance = 110.000)
	parameter (hcer_zmirror   = 245.000)
	parameter (hcer_zexit     = 265.000)

C Calorimeter position
	real*8 hcal_1pr_zpos,hcal_2ta_zpos,hcal_3ta_zpos,hcal_4ta_zpos
	real*8 hcal_left,hcal_right,hcal_top,hcal_bottom
	parameter (hcal_1pr_zpos = 338.69)
	parameter (hcal_2ta_zpos = 349.69)
	parameter (hcal_3ta_zpos = 360.69)
	parameter (hcal_4ta_zpos = 371.69)
	parameter (hcal_left     =  35.00)	!actual size of calorimeter
	parameter (hcal_right    = -35.00)
	parameter (hcal_top      = -69.66)
	parameter (hcal_bottom   =  60.34)

C The arguments

	real*8 p,m2			!momentum and mass of particle
	real*8 x_fp,y_fp,dx_fp,dy_fp	!Focal plane values to return
	real*8 xcal,ycal 		!Position of track at calorimeter.
	real*8 zinit			!Initial z-position (Not at F.P.)
	real*8 pathlen
	logical ms_flag			!mult. scattering flag.
	logical wcs_flag		!wire chamber smearing flag
	logical decay_flag		!check for decay
	logical ok_hut			!true if particle makes it

C Local declarations.

	integer*4 i,iplane,jchamber,npl_off
	integer*4 scintrig, scincount
	parameter (scintrig = 3)	!set trigger to 3/4

	logical dflag			!has particle decayed?

	real*8 resmult,tmpran,tmplim
	real*8 tmpran1,tmpran2		!temporary random numbers
	real*8 radw,drift

	real*8 nsig_max
	parameter(nsig_max=99.0e0)	!max #/sigma for gaussian ran #s.

C These have to be real*4 for the CERNLIB lfit routine.
	real*4 badf				!temporaries
	real*4 xfp4,yfp4,dxfp4,dyfp4		!real*4 versions of fp track.
	real*4 xdc(12),ydc(12),zdc(12)		!positions at d.c. planes

C ================================ Executable Code =============================

C Initialize some variables
C These come from Doug's examination of events with >6 hits per chamber.
C In some fraction of events (~10% ?) there were >6 hits per chamber.  
C The position resolution for these event had approx. 3 times worse resolution.
C There was also some position (or delta) depdendence.  The numbers below
C are based on examining one set of runs.  If you're going to use these
C to try and reproduce tails in the resolution, you should probably 
C see if the fraction with >6 hits, the resolution, and the delta depdendance
C are consistant with these numbers.  
C Note that if you play with the nominal DC resolutions, then you
C may want to reduce the resolution multiplier here.

! These values are from Doug's analysis of resolutions for the Kaon run.
!	tmpran = grnd()
!	tmplim = 0.13+0.008666*abs(dpps)  !fraction of events with >6 hits.
!	if (tmpran.lt.tmplim) then
!	  resmult = 3.2			  !resolution is 3.2x worse (both DCs)
!	else
!	  resmult = 1.0 + 0.05333*dpps
!	endif
!	tmpran = grnd()

! These are more conservative values (some mix of Kaon and NucPi values).
! Best to check what you see in your experiment.
	tmpran = grnd()
	tmplim = 0.15			!fraction of events with >6 hits.
	if (tmpran.lt.tmplim) then
	  resmult = 2.0
	else
	  resmult = 1.0
	endif

C Initialize ok_hut

	ok_hut = .false.

C Initialize the xdc and ydc arrays to zero

	do i=1,12
	  xdc(i) = 0.
	  ydc(i) = 0.
	enddo

C Initialize scincount to zero

	scincount = 0

C------------------------------------------------------------------------------C
C                           Top of loop through hut                            C
C------------------------------------------------------------------------------C

C Go to spectrometer exit foil, 25cm before DC1 (Drift forwards from zinit).
C As usual, neglect effect of nonzero dydzs and dxdzs on radw.

	drift = (hdc_1_zpos-25.000)-zinit	!Need to adjust number if known
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	radw = hfoil_exit_thick/hfoil_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Go to first drift chamber set
C For simplicity, perform air mult. scattering (probably negligeable)
C before drifting instead of 1/2 way through.

	drift = (hdc_1_zpos - 0.5*hdc_nr_plan*hdc_del_plane) -
     >		(hdc_1_zpos-25.000)
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	jchamber = 1
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen) !drift and decay
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane)=xs+hdc_sigma(npl_off+iplane)*tmpran1*resmult
	  ydc(npl_off+iplane)=ys+hdc_sigma(npl_off+iplane)*tmpran2*resmult
	  if (iplane.eq.2 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen) !drift and decay
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	if ( xs.gt.(hdc_1_bot-hdc_1x_offset) .or.
     >       xs.lt.(hdc_1_top-hdc_1x_offset) .or.
     >       ys.gt.(hdc_1_left-hdc_1y_offset) .or.
     >       ys.lt.(hdc_1_right-hdc_1y_offset) ) then
	  hSTOP_dc1 = hSTOP_dc1 + 1
	  goto 500
	endif
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of first drift chamber set, drift to next

	drift = (hdc_2_zpos - 0.5*hdc_nr_plan*hdc_del_plane) -
     >		(hdc_1_zpos + 0.5*hdc_nr_plan*hdc_del_plane)
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay 
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	jchamber = 2
	radw = hdc_entr_thick/hdc_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	npl_off = (jchamber-1)*hdc_nr_plan
	do iplane = 1,hdc_nr_plan
	  radw = hdc_cath_thick/hdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen) !drift and decay
	  radw = hdc_wire_thick/hdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane)=xs+hdc_sigma(npl_off+iplane)*tmpran1*resmult
	  ydc(npl_off+iplane)=ys+hdc_sigma(npl_off+iplane)*tmpran2*resmult
	  if (iplane.eq.2 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !y plane, no x information
	  else
	    ydc(npl_off+iplane) = 0.   !x-like plane, no y info
	  endif

	  drift = 0.5*hdc_thick
	  radw = drift/hdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + hdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen) !drift and decay
	enddo
	radw = hdc_exit_thick/hdc_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	if ( xs.gt.(hdc_2_bot-hdc_2x_offset) .or.
     >       xs.lt.(hdc_2_top-hdc_2x_offset) .or.
     >       ys.gt.(hdc_2_left-hdc_2y_offset) .or. 
     >       ys.lt.(hdc_2_right-hdc_2y_offset) ) then
	  hSTOP_dc2 = hSTOP_dc2 + 1
	  goto 500
	endif
	radw = hdc_cath_thick/hdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C fit track to give new focal plane values, use LFIT from GENLIB
	do jchamber=1,hdc_nr_cham
	  npl_off = (jchamber-1)*hdc_nr_plan
	  do iplane=1,hdc_nr_plan
	    if (jchamber.eq.1) zdc(npl_off+iplane) = hdc_1_zpos +
     >		(iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	    if (jchamber.eq.2) zdc(npl_off+iplane) = hdc_2_zpos +
     >		(iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
	  enddo
	enddo

	call lfit(zdc,xdc,12,0,dxfp4,xfp4,badf)
	call lfit(zdc,ydc,12,0,dyfp4,yfp4,badf)

	x_fp = dble(xfp4)
	y_fp = dble(yfp4)
	dx_fp = dble(dxfp4)
	dy_fp = dble(dyfp4)

C at last cathode foil of second drift chamber set, drift to aerogel

	drift = haer_zentrance - hdc_2_zpos - 0.5*hdc_nr_plan*hdc_del_plane
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = haer_entr_thick/haer_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = haer_thick
	radw = drift/haer_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	drift = haer_air_thick
	radw = drift/haer_air_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = haer_exit_thick/haer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at aerogel exit, drift to hodoscopes

	drift = hscin_1x_zpos - haer_zexit
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = scincount + 1
	radw = hscin_1x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hscin_1y_zpos - hscin_1x_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = scincount + 1
	radw = hscin_1y_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to cherenkov

	drift = hcer_zentrance - hscin_1y_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = hcer_entr_thick/hcer_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_zmirror - hcer_zentrance
	radw = drift/hcer_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

	radw = hcer_mir_thick/hcer_mir_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hcer_zexit - hcer_zmirror
	radw = hcer_exit_thick/hcer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay

	radw = hcer_exit_thick/hcer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to second hodoscope

	drift = hscin_2x_zpos - hcer_zexit
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = scincount + 1
	radw = hscin_2x_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = hscin_2y_zpos - hscin_2x_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = scincount + 1
	radw = hscin_2y_thick/hscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Test on scintillator trigger (DJG 8/18/98)
	if( scincount .lt. scintrig ) then
	  hSTOP_scin = hSTOP_scin + 1
	  goto 500
	endif

C Drift to back of calorimeter.  This should be done even if we do not apply
C any calorimeter cuts, so that pathlen is calculated to the back of the spec.
	drift = hcal_4ta_zpos - hscin_2y_zpos
	radw = drift/hair_radlen
	call project(xs,ys,drift,decay_flag,dflag,m2,p,pathlen)	!drift and decay
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

C Note that even with the standard PID trigger, the calorimeter is NOT
C required, since the trigger needs either the cerenkov OR the calorimeter.
C If you require the calorimeter in you analysis, then you need to make
C sure that the calorimeter is hit here, AND apply a seperate post-tracking
C fiducial cut (whatever is required by the engine or in your analysis).
C The fiducial cut comes later in the code.
C
C This cut simply requires that the event hit the calorimeter.  It does
C not insure that the event is within the good region of the calorimeter.
C The seperate fiducial cut is required to insure that the entire shower
C energy is contained in the calorimeter.  Or, if you like, you can require
C some distance between the track and the edge.  2-3 cm seems to be enough
C to get most or all of the energy in the calorimeter

!c	if (ys.gt.(hcal_left-2.0) .or. ys.lt.(hcal_right+2.0) .or.
!c     >     xs.gt.(hcal_bottom-2.0) .or. xs.lt.(hcal_top+2.0)) then
!	if (ys.gt.hcal_left .or. ys.lt.hcal_right .or.
!     >      xs.gt.hcal_bottom .or. xs.lt.hcal_top) then
!	  hSTOP_cal = hSTOP_cal + 1
!	  goto 500
!	endif

C If you use a calorimeter cut in your analysis, the engine applied a
C a fiducial cut at the calorimeter.  This is based on the position of the
C TRACK at the calorimeter, not the real position of the event.  Go to the
C back of the calorimeter since engine uses a FID cut at the back.
C The standard fiducial cut is 5 cm from the edges of the block.

	xcal = x_fp + dx_fp * hcal_4ta_zpos
	ycal = y_fp + dy_fp * hcal_4ta_zpos
!	if (ycal.gt.(hcal_left-5.0) .or. ycal.lt.(hcal_right+5.0) .or.
!     >     xcal.gt.(hcal_bottom-5.0) .or. xcal.lt.(hcal_top+5.0)) then

	ok_hut = .true.

C We are done with this event, whether GOOD or BAD.

500	continue

	return
	end
