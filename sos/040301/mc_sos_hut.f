	subroutine mc_sos_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >		decay_flag,dflag,resmult,ok_hut,zinit)

C----------------------------------------------------------------------
C
C Monte-Carlo of SOS detector hut.
C	Note that we only pass on real*8 variables to the subroutine.
C	This will make life easier for portability.
C
C	The particle is stepped through the detector (using project), and
C	multiple scattering is applied for each detector or air gap.
C	If particle decay in enabled, then project.f also checks for
C	decay of particle.  The particle starts at z=zinit.  This
C	needs to be before the first mult. scattering point (the exit window)
C	or the decay distance is negative, and things are BAD.
C
C----------------------------------------------------------------------

	implicit 	none

	include '../struct_sos.inc'
	include '../track.inc'

C Math constants

	real*8 pi,d_r,r_d,root
	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (r_d = 180./pi)
	parameter (root = 0.707106781)		!square root of 1/2

	real*8 grnd,gauss1			!external functions

C----------------------------------------------------------------------
C SOS_MATERIALS
C CTP parameter file containing the materials of all the SOS detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C----------------------------------------------------------------------
C spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
	real*8 sfoil_exit_radlen,sfoil_exit_thick
	real*8 sfoil_exit_zpos,sfoil_exit_ang
	parameter (sfoil_exit_radlen = 53.3)
	parameter (sfoil_exit_thick  = 0.020*2.54)
	parameter (sfoil_exit_zpos = -3.22)
	parameter (sfoil_exit_ang = 0.*d_r)  !50 deg. BEFORE vacuum extension,
					     ! 0 deg. AFTER installation.
C spectrometer air gaps
	real*8 sair_radlen
	parameter (sair_radlen = 30420.)		!radlen, cm

C chamber gas
	real*8 sdc_radlen,sdc_thick
	parameter (sdc_radlen = 16700.0)
	parameter (sdc_thick  = 0.61775)

C effective wire plane ! 60 micron field wires+30 micron sense wires, Tungsten
	real*8 sdc_wire_radlen,sdc_wire_thick
	parameter (sdc_wire_radlen = 0.35)
	parameter (sdc_wire_thick  = 0.0000354)

C chamber cathode foils, 1/2 mil of Mylar
	real*8 sdc_cath_radlen,sdc_cath_thick
	parameter (sdc_cath_radlen = 28.7)
	parameter (sdc_cath_thick  = 0.0005*2.54)

C hodoscopes
	real*8 sscin_radlen
	parameter (sscin_radlen =  42.4)

C Cherenkov entrance foil, 0.5 mm Al +
	real*8 scer_entr_radlen,scer_entr_thick
	parameter (scer_entr_radlen = 8.90)
	parameter (scer_entr_thick  = 0.050)

C Cherenkov, 1 atm Freon 12 +
	real*8 scer_radlen
	parameter (scer_radlen = 4810.0)

C Cherenkov mirror, 75 mu plus 2 cm Rotacell +
	real*8 scer_mir_radlen,scer_mir_thick
	parameter (scer_mir_radlen = 400.0)
	parameter (scer_mir_thick  = 2.0)

C Cherenkov exit foil, 0.5 mm Al +
	real*8 scer_exit_radlen,scer_exit_thick
	parameter (scer_exit_radlen = 8.90)
	parameter (scer_exit_thick  = 0.050)

C Aerogel entrance foil
	real*8 saer_entr_radlen,saer_entr_thick
	parameter (saer_entr_radlen = 8.90)
	parameter (saer_entr_thick  = 0.0625)

C Aerogel +
	real*8 saer_radlen,saer_thick
	parameter (saer_radlen = 150.0)
	parameter (saer_thick  = 9.0)

C Aerogel exit foil
	real*8 saer_exit_radlen,saer_exit_thick
	parameter (saer_exit_radlen = 8.90)
	parameter (saer_exit_thick  = 0.0625)

C shower counter
	real*8 scal_radlen
	parameter (scal_radlen = 2.50)

C Wire chamber resolutions (sigma)
	real*8 sdc_sigma(1:12)/	0.030,0.030,0.030,0.030,0.030,0.030,
     >  			0.030,0.030,0.030,0.030,0.030,0.030/

C Wire plane positions, construct sdc_zpos array using these parameters
	integer*4 sdc_nr_cham,sdc_nr_plan
	parameter (sdc_nr_cham = 2)
	parameter (sdc_nr_plan = 6)

	real*8 sdc_1_zpos,sdc_1_left,sdc_1_right,sdc_1_top,sdc_1_bot
	real*8 sdc_1x_offset,sdc_1y_offset
	real*8 sdc_2_zpos,sdc_2_left,sdc_2_right,sdc_2_top,sdc_2_bot
	real*8 sdc_2x_offset,sdc_2y_offset
	real*8 sdc_del_plane
	parameter (sdc_1_zpos =   6.25)
	parameter (sdc_2_zpos =  55.77)
	parameter (sdc_del_plane = sdc_thick + sdc_wire_thick + sdc_cath_thick)
	parameter (sdc_1_left  =  24.0)
	parameter (sdc_1_right = -24.0)
c	parameter (sdc_1y_offset = -1.769)
	parameter (sdc_1y_offset = -1.822)
	parameter (sdc_1_top   = -32.0)
	parameter (sdc_1_bot   =  32.0)
c	parameter (sdc_1x_offset = 8.0226)
	parameter (sdc_1x_offset = 8.649)
	parameter (sdc_2_left  =  24.0)
	parameter (sdc_2_right = -24.0)
c	parameter (sdc_2y_offset = -1.794)
	parameter (sdc_2y_offset = -1.976)
	parameter (sdc_2_top   = -32.0)
	parameter (sdc_2_bot   =  32.0)
c	parameter (sdc_2x_offset = -2.231)
	parameter (sdc_2x_offset = -1.532)

C Scintillator positions and thicknesses
	real*8 sscin_1x_zpos,sscin_1y_zpos
	real*8 sscin_1x_thick,sscin_1y_thick
	real*8 sscin_1_left,sscin_1_right,sscin_1x_offset
	real*8 sscin_1_top,sscin_1_bot,sscin_1y_offset
	real*8 sscin_2x_zpos,sscin_2y_zpos
	real*8 sscin_2x_thick,sscin_2y_thick
	real*8 sscin_2_left,sscin_2_right,sscin_2x_offset
	real*8 sscin_2_top,sscin_2_bot,sscin_2y_offset
	parameter (sscin_1y_zpos =  73.61)
	parameter (sscin_1x_zpos =  97.11)
	parameter (sscin_2y_zpos = 249.51)
	parameter (sscin_2x_zpos = 290.81)
	parameter (sscin_1x_thick = 1.040) !1.040 for overlap
	parameter (sscin_1y_thick = 1.098)
	parameter (sscin_2x_thick = 1.040)
	parameter (sscin_2y_thick = 1.098)
	parameter (sscin_1_left  =  18.25)
	parameter (sscin_1_right = -18.25)
	parameter (sscin_1x_offset = 2.8)
	parameter (sscin_1_top   = -31.75)
	parameter (sscin_1_bot   =  31.75)
	parameter (sscin_1y_offset = 2.25)
	parameter (sscin_2_left  =  18.25)
	parameter (sscin_2_right = -18.25)
	parameter (sscin_2x_offset = 4.9)
	parameter (sscin_2_top   = -56.25)
	parameter (sscin_2_bot   =  56.25)
	parameter (sscin_2y_offset = 2.9)

C Cherenkov position
	real*8 scer_zentrance,scer_zmirror,scer_zexit
	parameter (scer_zentrance = 130.000)
	parameter (scer_zmirror   = 155.000)
	parameter (scer_zexit     = 160.000)

C Calorimeter position
	real*8 scal_1pr_zpos,scal_2ta_zpos,scal_3ta_zpos,scal_4ta_zpos
	real*8 scal_left,scal_right,scal_top,scal_bottom
	parameter (scal_1pr_zpos = 313.01)
	parameter (scal_2ta_zpos = 324.01)
	parameter (scal_3ta_zpos = 335.01)
	parameter (scal_4ta_zpos = 346.01)
	parameter (scal_left     =  35.00)
	parameter (scal_right    = -35.00)
	parameter (scal_top      = -61.)
	parameter (scal_bottom   =  49.)

C The arguments

	real*8 p,m2			!momentum and mass of particle
	real*8 x_fp,y_fp,dx_fp,dy_fp	!Focal plane values to return
	real*8 xcal,ycal		!Position of track at calorimeter.
	real*8 zinit			!Initial z-position (Not at F.P.)
	logical ms_flag			!mult. scattering flag. 
	logical wcs_flag		!wire chamber smearing flag
	logical decay_flag		!check for decay
	logical ok_hut			!true if particle makes it

C Local declarations.

	integer*4 i,iplane,jchamber,npl_off
	integer*4 scintrig,scincount
	parameter (scintrig = 3)	!set trigger to 3/4

	logical*4 dflag			!has particle decayed?

	real*8 resmult,tmpran,tmplim
	real*8 tmpran1,tmpran2		!temporary random numbers
	real*8 xt
	real*8 radw,drift

	real*8 nsig_max
	parameter(nsig_max=99.0d0)      !max #/sigma for gaussian ran #s.

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
!	tmplim = 0.15			!fraction of events with >6 hits.
!	if (tmpran.lt.tmplim) then
!	  resmult=3.0			!resolution is 3.0x worse (both DCs)
!	else
!	  resmult=1.0
!	endif

! These are more conservative values (some mix of Kaon and NucPi values).
! Best to check what you see in your experiment.
	tmpran = grnd()
	tmplim = 0.15			!fraction of events with >6 hits.
	if (tmpran.lt.tmplim) then
	  resmult=2.0			!resolution is 3.0x worse (both DCs)
	else
	  resmult=1.0
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

C Go to spectrometer exit foil. (Drift forwards from zinit).
C Approximation is used in calculating xt = x at foil crossing.
C As usual, neglect effect of nonzero dydzs and dxdzs on radw.

	xt = xs + dxdzs * sfoil_exit_zpos
	drift = (sfoil_exit_zpos + xt * tan(sfoil_exit_ang)) - zinit 
	if (drift.le.0.001) drift=0.001		!avoid drift<=0
	call project(xs,ys,drift,decay_flag,dflag,m2,p)		!drift and decay
	radw = sfoil_exit_thick/sfoil_exit_radlen/cos(sfoil_exit_ang)
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

!Already decayed for the 3.22 cm, so it get double counted in the next
!step.  Oh well.  Fix this later.


C Go to first drift chamber set
C For simplicity, perform air MS (probably negligeable) before DCs
C instead of 1/2 way through.

	drift = (sdc_1_zpos - 0.5*sdc_nr_plan*sdc_del_plane) - 
     >		(sfoil_exit_zpos + xt * tan(sfoil_exit_ang))
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

	jchamber = 1
	npl_off = (jchamber-1)*sdc_nr_plan
	do iplane = 1,sdc_nr_plan
	  radw = sdc_cath_thick/sdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*sdc_thick
	  radw = drift/sdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + sdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	  radw = sdc_wire_thick/sdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane)=xs+sdc_sigma(npl_off+iplane)*tmpran1*resmult
	  ydc(npl_off+iplane)=ys+sdc_sigma(npl_off+iplane)*tmpran2*resmult
	  if (iplane.eq.1 .or. iplane.eq.3 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !assume 3 planes are x, 3 are y
	  else
	    ydc(npl_off+iplane) = 0.
	  endif
	  drift = 0.5*sdc_thick
	  radw = drift/sdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + sdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	enddo
	if ( xs.gt.(sdc_1_bot-sdc_1x_offset) .or. 
     >       xs.lt.(sdc_1_top-sdc_1x_offset) .or.
     >       ys.gt.(sdc_1_left-sdc_1y_offset) .or.
     >       ys.lt.(sdc_1_right-sdc_1y_offset) ) then
	  sSTOP.detec = sSTOP.detec + 1
	  goto 500
	endif
	radw = sdc_cath_thick/sdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of first drift chamber set, drift to next

	drift = sdc_2_zpos - sdc_1_zpos - sdc_nr_plan*sdc_del_plane
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

	jchamber = 2
	npl_off = (jchamber-1)*sdc_nr_plan
	do iplane = 1,sdc_nr_plan
	  radw = sdc_cath_thick/sdc_cath_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = 0.5*sdc_thick
	  radw = drift/sdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + sdc_cath_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	  radw = sdc_wire_thick/sdc_wire_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  if(wcs_flag) then
	    tmpran1 = gauss1(nsig_max)	!gaussian, truncated at 99 sigma.
	    tmpran2 = gauss1(nsig_max)
	  else
	    tmpran1 = 0.
	    tmpran2 = 0.
	  endif
	  xdc(npl_off+iplane)=xs+sdc_sigma(npl_off+iplane)*tmpran1*resmult
	  ydc(npl_off+iplane)=ys+sdc_sigma(npl_off+iplane)*tmpran2*resmult
	  if (iplane.eq.1 .or. iplane.eq.3 .or. iplane.eq.5) then
	    xdc(npl_off+iplane) = 0.   !assume 3 planes are x, 3 are y
	  else
	    ydc(npl_off+iplane) = 0.
	  endif

	  drift = 0.5*sdc_thick
	  radw = drift/sdc_radlen
	  if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
	  drift = drift + sdc_wire_thick
	  call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	enddo
	if ( xs.gt.(sdc_2_bot-sdc_2x_offset) .or.
     >       xs.lt.(sdc_2_top-sdc_2x_offset) .or.
     >       ys.gt.(sdc_2_left-sdc_2y_offset) .or.
     >       ys.lt.(sdc_2_right-sdc_2y_offset) ) then
	  sSTOP.detec = sSTOP.detec + 1
	  goto 500
	endif
	radw = sdc_cath_thick/sdc_cath_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C fit track to give focal plane values, use LFIT from GENLIB
C This allows us to compare expected hit position versus hodoscope position,
C and throw out events which will fail to give a beta value (if so desired).

	do jchamber=1,sdc_nr_cham
	  npl_off = (jchamber-1)*sdc_nr_plan
	  do iplane=1,sdc_nr_plan
	    if (jchamber.eq.1) zdc(npl_off+iplane) = sdc_1_zpos +
     >			(iplane-0.5-0.5*sdc_nr_plan)*sdc_del_plane
	    if (jchamber.eq.2) zdc(npl_off+iplane) = sdc_2_zpos +
     >			(iplane-0.5-0.5*sdc_nr_plan)*sdc_del_plane
	  enddo
	enddo

	call lfit(zdc,xdc,12,0,dxfp4,xfp4,badf)
	call lfit(zdc,ydc,12,0,dyfp4,yfp4,badf)
	x_fp = dble(xfp4)
	y_fp = dble(yfp4)
	dx_fp = dble(dxfp4)
	dy_fp = dble(dyfp4)

C at last cathode foil of second drift chamber set, drift to hodoscopes

	drift = sscin_1y_zpos - sdc_2_zpos - 0.5*sdc_nr_plan*sdc_del_plane
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	if (ys.lt.(sscin_1_left+sscin_1y_offset) .and.
     >      ys.gt.(sscin_1_right+sscin_1y_offset) .and.
     >      xs.lt.(sscin_1_bot+sscin_1x_offset) .and.
     >      xs.gt.(sscin_1_top+sscin_1x_offset)) scincount = scincount + 1
	radw = sscin_1y_thick/sscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = sscin_1x_zpos - sscin_1y_zpos
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	if (ys.lt.(sscin_1_left+sscin_1y_offset) .and.
     >      ys.gt.(sscin_1_right+sscin_1y_offset) .and.
     >      xs.lt.(sscin_1_bot+sscin_1x_offset) .and.
     >      xs.gt.(sscin_1_top+sscin_1x_offset)) scincount = scincount + 1
	radw = sscin_1x_thick/sscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to cherenkov

	drift = scer_zentrance - sscin_1x_zpos
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

	radw = scer_entr_thick/scer_entr_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = scer_zmirror - scer_zentrance
	radw = drift/scer_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

	radw = scer_mir_thick/scer_mir_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = scer_zexit - scer_zmirror
	radw = drift/scer_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

	radw = scer_exit_thick/scer_exit_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to second hodoscope

	drift = sscin_2y_zpos - scer_zexit
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	if (ys.lt.(sscin_2_left+sscin_2y_offset) .and.
     >      ys.gt.(sscin_2_right+sscin_2y_offset) .and.
     >      xs.lt.(sscin_2_bot+sscin_2x_offset) .and.
     >      xs.gt.(sscin_2_top+sscin_2x_offset)) scincount = scincount + 1
	radw = sscin_2y_thick/sscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

	drift = sscin_2x_zpos - sscin_2y_zpos
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay
	if (ys.lt.(sscin_2_left+sscin_2y_offset) .and.
     >      ys.gt.(sscin_2_right+sscin_2y_offset) .and.
     >      xs.lt.(sscin_2_bot+sscin_2x_offset) .and.
     >      xs.gt.(sscin_2_top+sscin_2x_offset)) scincount = scincount + 1
	radw = sscin_2x_thick/sscin_radlen
	if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Test on scintillator trigger (DJG 8/18/98)
	if( scincount .lt. scintrig ) then
	  sSTOP.detec = sSTOP.detec + 1
	  goto 500
	endif

C Drift to calorimeter
	drift = scal_4ta_zpos - sscin_2x_zpos
	radw = drift/sair_radlen
	if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
	call project(xs,ys,drift,decay_flag,dflag,m2,p)	!drift and decay

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

!c	if (ys.gt.(scal_left-2.0) .or. ys.lt.(scal_right+2.0) .or.
!c     >     xs.gt.(scal_bottom-2.0) .or. xs.lt.(scal_top+2.0)) then
!	if (ys.gt.scal_left .or. ys.lt.scal_right .or.
!     >      xs.gt.scal_bottom .or. xs.lt.scal_top) then
!	    sSTOP.detec = sSTOP.detec + 1
!	    goto 500
!	endif

C If you use a calorimeter cut in your analysis, the engine applied a
C a fiducial cut at the calorimeter.  This is based on the position of the
C TRACK at the calorimeter, not the real position of the event.  Go to the
C back of the calorimeter since engine uses a FID cut at the back.
C The standard fiducial cut is 5 cm from the edges of the block.

	xcal = x_fp + dx_fp * scal_4ta_zpos
	ycal = y_fp + dy_fp * scal_4ta_zpos
!	if (ycal.gt.(scal_left-5.0) .or. ycal.lt.(scal_right+5.0) .or.
!     >     xcal.gt.(scal_bottom-5.0) .or. xcal.lt.(scal_top+5.0)) then
!	  sSTOP.detec = sSTOP.detec + 1
!	  goto 500
!	endif

	ok_hut = .true.

C We are done with this event, whether GOOD or BAD.

500	continue

	return
	end
