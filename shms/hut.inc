C----------------------------------------------------------------------
C HMS_MATERIALS
C CTP parameter file containing the materials of all the HMS detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C
C	6/15/00		changed the wire chamber resolution to 250 um.
C----------------------------------------------------------------------
C spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
!	real*8 hfoil_exit_radlen,hfoil_exit_thick
!	parameter (hfoil_exit_radlen = 53.3)
!	parameter (hfoil_exit_thick  = 0.020*2.54)

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
	parameter (hdc_thick  = 0.61775)

C effective wire plane, 25 micron W+Be/Cu gives <t>=pi*(.0025/2)**2
	real*8 hdc_wire_radlen,hdc_wire_thick
	parameter (hdc_wire_radlen = 0.35)	!Assuming all Tungsten
	parameter (hdc_wire_thick  = 0.0000354)

C effective cathode plane, Be/Cu
	real*8 hdc_cath_radlen,hdc_cath_thick
	parameter (hdc_cath_radlen = 28.7)	!'Ave' of Be/Cu
	parameter (hdc_cath_thick  = 0.0005*2.54)

C chamber exit foil, 1 mil of Mylar
	real*8 hdc_exit_radlen,hdc_exit_thick
	parameter (hdc_exit_radlen = 28.7)
	parameter (hdc_exit_thick  = 0.001*2.54)

C hodoscopes
	real*8 hscin_radlen
	parameter (hscin_radlen =  42.4)

C Cherenkov entrance foil, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
	real*8 hcer_entr_radlen,hcer_entr_thick
	parameter (hcer_entr_radlen = 53.3)
	parameter (hcer_entr_thick  = 0.020*2.54)

C Cherenkov entrance foil, 40 mil of Al
	real*8 hcer_2_entr_radlen,hcer_2_entr_thick
	parameter (hcer_2_entr_radlen = 8.90)
	parameter (hcer_2_entr_thick  = 0.040*2.54)

C Cherenkov, 0.5 atm of CO2 for e/pi
C	to determine radlen ~ num * 73240
C	real*8 hcer_1_radlen
C	parameter (hcer_1_radlen = 36620.0)

C Cherenkov, 0.3 atm of CO2 for e/pi
	real*8 hcer_1_radlen
	parameter (hcer_1_radlen = 21972.0)

C Cherenkov, 0.5 atm of Freon for better e/pi
C	real*8 hcer_1_radlen
C	parameter (hcer_1_radlen = 9620.0)

C Cherenkov, 2 atm of Freon for pi/p
	real*8 hcer_2_radlen
	parameter (hcer_2_radlen = 1202.5)

C Cherenkov, 2 atm of CO2 for pi/p
C	real*8 hcer_2_radlen
C	parameter (hcer_2_radlen = 14680.0)

C Cherenkov mirror, 75 mu plus 2 cm of Rotacell +
	real*8 hcer_mir_radlen,hcer_mir_thick
	parameter (hcer_mir_radlen = 400.0)
	parameter (hcer_mir_thick  = 2.0)

C Cherenkov exit foil, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
	real*8 hcer_exit_radlen,hcer_exit_thick
	parameter (hcer_exit_radlen = 53.3)
	parameter (hcer_exit_thick  = 0.020*2.54)

C Cherenkov exit foil, 40 mil of Al
	real*8 hcer_2_exit_radlen,hcer_2_exit_thick
	parameter (hcer_2_exit_radlen = 8.90)
	parameter (hcer_2_exit_thick  = 0.040*2.54)

C shower counter
	real*8 hcal_radlen
	parameter (hcal_radlen = 2.50)

C Wire chamber resolutions (sigma)
c SHMS chambers assume sigma=170um per plane
c   planes will be two sets of 2u,2v,2x but need to mock it by 2 sets of 4x,2y
c      planes 2 and 5 in group of 6 are y planes
	real*8 hdc_sigma(1:12)/	0.024,0.014,0.024,0.017,0.014,0.017,
     >				0.024,0.014,0.024,0.017,0.014,0.017/

C Wire plane positions, construct hdc_zpos array using these parameters

	integer*4 hdc_nr_cham,hdc_nr_plan
c SHMS chambers
	parameter (hdc_nr_cham = 2)
	parameter (hdc_nr_plan = 6)

	real*8 hdc_1_zpos,hdc_1_left,hdc_1_right,hdc_1_top,hdc_1_bot
	real*8 hdc_1x_offset,hdc_1y_offset
	real*8 hdc_2_zpos,hdc_2_left,hdc_2_right,hdc_2_top,hdc_2_bot
	real*8 hdc_2x_offset,hdc_2y_offset


c
c for LSA tune the focal plane is at +100cm compared to SSA tune
	real*8 lsa_shift
	parameter (lsa_shift=0.0)
c
	parameter (hdc_1_zpos = -40.0+lsa_shift)
	parameter (hdc_2_zpos =  40.0+lsa_shift)
	parameter (hdc_1_left  =  50.0)
	parameter (hdc_1_right = -50.0)
	parameter (hdc_1y_offset = 0.0)
	parameter (hdc_1_top   = -40.)
	parameter (hdc_1_bot   =  45.)
	parameter (hdc_1x_offset = 0.0)
	parameter (hdc_2_left  =  50.0)
	parameter (hdc_2_right = -50.0)
	parameter (hdc_2y_offset = 0.0)
	parameter (hdc_2_top   = -40.)
	parameter (hdc_2_bot   =  45.)
	parameter (hdc_2x_offset =  0.)

C Scintillator positions and thicknesses

	real*8 hscin_1x_zpos
	real*8 hscin_1x_thick
	real*8 hscin_1x_left,hscin_1x_right,hscin_1x_offset
	real*8 hscin_1y_zpos
	real*8 hscin_1y_thick
	real*8 hscin_1y_top,hscin_1y_bot,hscin_1y_offset
	real*8 hscin_2x_zpos,hscin_2y_zpos
	real*8 hscin_2x_thick,hscin_2y_thick
	real*8 hscin_2x_left,hscin_2x_right,hscin_2x_offset
	real*8 hscin_2y_top,hscin_2y_bot,hscin_2y_offset

	parameter (hscin_1x_zpos = 50.0+lsa_shift)
	parameter (hscin_1y_zpos = 60.0+lsa_shift)
	parameter (hscin_2x_zpos = 260.0+lsa_shift)
	parameter (hscin_2y_zpos = 265.0+lsa_shift)
	parameter (hscin_1x_thick = 1.000*1.067)
	parameter (hscin_1y_thick = 1.000*1.067)
	parameter (hscin_2x_thick = 1.000*1.067)
	parameter (hscin_2y_thick = 1.000*1.067)
	parameter (hscin_1x_left  =  55.)
	parameter (hscin_1x_right = -55.)
	parameter (hscin_1x_offset = 0.0)
	parameter (hscin_1y_top   = -45.)
	parameter (hscin_1y_bot   =  45.)
	parameter (hscin_1y_offset = 0.0)
	parameter (hscin_2x_left  =   70.)
	parameter (hscin_2x_right =  -70.)
	parameter (hscin_2x_offset =  0.)
	parameter (hscin_2y_top   = -55.)
	parameter (hscin_2y_bot   =  55.)
	parameter (hscin_2y_offset = 0)


C Cherenkov position

	real*8 hcer_1_zentrance,hcer_1_zmirror,hcer_1_zexit
	real*8 hcer_2_zentrance,hcer_2_zmirror,hcer_2_zexit
	parameter (hcer_1_zentrance = -360.000+lsa_shift)
	parameter (hcer_1_zmirror   = -100.000+lsa_shift)
	parameter (hcer_1_zexit     = -80.000+lsa_shift)
	parameter (hcer_2_zentrance = 80.000+lsa_shift)
	parameter (hcer_2_zmirror   = 160.000+lsa_shift)
	parameter (hcer_2_zexit     = 180.000+lsa_shift)

C Calorimeter position

	real*8 hcal_1pr_zpos,hcal_2ta_zpos,hcal_3ta_zpos,hcal_4ta_zpos
	real*8 hcal_left,hcal_right,hcal_top,hcal_bottom
	parameter (hcal_1pr_zpos = 270.0+lsa_shift)
	parameter (hcal_2ta_zpos = 288+lsa_shift)
	parameter (hcal_3ta_zpos = 306+lsa_shift)
	parameter (hcal_4ta_zpos = 320+lsa_shift)
	parameter (hcal_left     =  75.00)
	parameter (hcal_right    = -75.00)
	parameter (hcal_top      = -60.00)
	parameter (hcal_bottom   =  55.00)
