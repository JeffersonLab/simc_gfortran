C----------------------------------------------------------------------
C SHMS_MATERIALS PRELIMINARY
C CTP parameter file containing the materials of the SHMS detectors 
C currently based on HMS materials and initial SHMS material selections.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C
C SHMS wire chambers assume resolution 170 um per plane
C
C----------------------------------------------------------------------
C spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
c	real*8 hfoil_exit_radlen,hfoil_exit_thick
c	parameter (hfoil_exit_radlen = 53.3)
c	parameter (hfoil_exit_thick  = 0.020*2.54) ! 0.1%
C spectrometer exit window, 0.02in thick Al
	real*8 hfoil_exit_radlen,hfoil_exit_thick
	parameter (hfoil_exit_radlen = 8.89)
	parameter (hfoil_exit_thick  = 0.020*2.54)

C spectrometer air gaps
	real*8 hair_radlen
	parameter (hair_radlen = 30420.)
C helium bag 
	real*8 helbag_hel_radlen
	parameter (helbag_hel_radlen = 754560.) !0.03%
	real*8 helbag_mylar_radlen,helbag_mylar_thick
	parameter (helbag_mylar_radlen = 28.53)
	parameter (helbag_mylar_thick = .010*2.54)  ! 0.1%
	real*8 helbag_al_radlen,helbag_al_thick
	parameter (helbag_al_radlen = 8.89)
	parameter (helbag_al_thick = 3e-4)   !0.03%

C chamber entrance foil, 1 mil of Mylar ( ignore 1000A of aluminum)
	real*8 hdc_entr_radlen,hdc_entr_thick
	parameter (hdc_entr_radlen = 28.7)
	parameter (hdc_entr_thick  = 0.001*2.54)

C chamber gas, 50/50 ethane/argon 1/8inch between cathode and sense plane
	real*8 hdc_radlen,hdc_thick
	parameter (hdc_radlen = 16700.0)
	parameter (hdc_thick  = 0.125*2.54)

C effective wire plane, 25 micron W+Be/Cu gives <t>=pi*(.0025/2)**2
	real*8 hdc_wire_radlen,hdc_wire_thick
	parameter (hdc_wire_radlen = 0.35)	!Assuming all Tungsten
	parameter (hdc_wire_thick  = 0.0000354)

C effective cathode plane, 1 mil Kapton with 1200A of copper ( ignore copper)
	real*8 hdc_cath_radlen,hdc_cath_thick
	parameter (hdc_cath_radlen = 28.6)	!
	parameter (hdc_cath_thick  = 0.001*2.54)

C chamber exit foil, 1 mil of Mylar
	real*8 hdc_exit_radlen,hdc_exit_thick
	parameter (hdc_exit_radlen = 28.7)
	parameter (hdc_exit_thick  = 0.001*2.54)

C hodoscopes
	real*8 hscin_radlen
	parameter (hscin_radlen =  42.4)

C Ar/neo Cherenkov entrance foil, Tedlar ( update from D. Day spreadsheet Feb 2016)
	real*8 hcer_entr_radlen,hcer_entr_thick
	parameter (hcer_entr_radlen = 19.63)
	parameter (hcer_entr_thick  = 0.002*2.54)
C Ar/neo Cherenkov, for now 100% argon
	real*8 hcer_1_radlen
	parameter (hcer_1_radlen = 11700.0)
C Ar/neo Cherenkov mirror, SiO2 , no backing
	real*8 hcer_mirglass_radlen,hcer_mirglass_thick
	parameter (hcer_mirglass_radlen = 12.29)
	parameter (hcer_mirglass_thick  = 0.3)
C Ar/neo Cherenkov exit foil Tedlar  ( update from D. Day spreadsheet Feb 2016)
	real*8 hcer_exit_radlen,hcer_exit_thick
	parameter (hcer_exit_radlen = 19.63)
	parameter (hcer_exit_thick  = 0.002*2.54)



C Cherenkov entrance foil, 40 mil of Al
	real*8 hcer_2_entr_radlen,hcer_2_entr_thick
	parameter (hcer_2_entr_radlen = 8.90)
	parameter (hcer_2_entr_thick  = 0.040*2.54)

C Cherenkov, 0.5 atm of CO2 for e/pi
C	to determine radlen ~ num * 73240
C	real*8 hcer_1_radlen
C	parameter (hcer_1_radlen = 36620.0)


C Cherenkov, 0.5 atm of Freon for better e/pi
C	real*8 hcer_1_radlen
C	parameter (hcer_1_radlen = 9620.0)

C Cherenkov, 2 atm of Freon for pi/p
	real*8 hcer_2_radlen
	parameter (hcer_2_radlen = 1202.5)

C Cherenkov, 2 atm of CO2 for pi/p
C	real*8 hcer_2_radlen
C	parameter (hcer_2_radlen = 14680.0)

	real*8 hcer_mir_radlen,hcer_mir_thick
	parameter (hcer_mir_radlen = 400.)
	parameter (hcer_mir_thick  = 2.00)


C Cherenkov exit foil, 40 mil of Al
	real*8 hcer_2_exit_radlen,hcer_2_exit_thick
	parameter (hcer_2_exit_radlen = 8.90)
	parameter (hcer_2_exit_thick  = 0.040*2.54)

C shower counter
	real*8 hcal_radlen
	parameter (hcal_radlen = 2.50)

C Wire chamber resolutions (sigma)
c SHMS chambers assume sigma=200um per plane
c   planes will be two sets of 2u,2v,2x but need to mock it by 2 sets of 4x,2y
c      planes 2 and 5 in group of 6 are y planes
c	real*8 hdc_sigma(1:12)/	0.024,0.014,0.024,0.017,0.014,0.017,
c     >				0.024,0.014,0.024,0.017,0.014,0.017/
	real*8 hdc_sigma(1:12)/	0.020,0.020,0.020,0.020,0.020,0.020,
     >				0.020,0.020,0.020,0.020,0.020,0.020/

C Wire plane positions, construct hdc_zpos array using these parameters

	integer*4 hdc_nr_cham,hdc_nr_plan
c SHMS chambers
	parameter (hdc_nr_cham = 2)
	parameter (hdc_nr_plan = 6)

	real*8 hdc_1_zpos,hdc_1_left,hdc_1_right,hdc_1_top,hdc_1_bot
	real*8 hdc_1x_offset,hdc_1y_offset
	real*8 hdc_2_zpos,hdc_2_left,hdc_2_right,hdc_2_top,hdc_2_bot
	real*8 hdc_2x_offset,hdc_2y_offset

c SHMS Cerenkov
	real*8 hcer_left, hcer_right
	real*8 hcer_top, hcer_bottom

C	parameter (hdc_1_zpos = -40.0)     !Original (NOT surveyed) DC1 zPos C.Y. 07/17/2018
	parameter (hdc_1_zpos = -40.656)
C	parameter (hdc_2_zpos =  40.0)     !Original (NOT surveyed) DC2 zPos C.Y. 07/17/2018
	parameter (hdc_2_zpos =  39.332)
	parameter (hdc_1_left  =  40.0)
	parameter (hdc_1_right = -40.0)
	parameter (hdc_1y_offset = 0.0)
	parameter (hdc_1_top   = -40.)
	parameter (hdc_1_bot   =  40.)
	parameter (hdc_1x_offset = 0.0)
	parameter (hdc_2_left  =  40.0)
	parameter (hdc_2_right = -40.0)
	parameter (hdc_2y_offset = 0.0)
	parameter (hdc_2_top   = -40.)
	parameter (hdc_2_bot   =  40.)
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

C       parameter (hscin_1x_zpos = 51.9)    !Original Hodo1x C.Y. 07/17/2018
	parameter (hscin_1x_zpos = 52.1)
C	parameter (hscin_1y_zpos = 60.9)    !Original Hodo1y C.Y. 07/17/2018
	parameter (hscin_1y_zpos = 61.7)
C	parameter (hscin_2x_zpos = 271.9)   !Original Hodo2x C.Y. 07/17/2018
	parameter (hscin_2x_zpos = 271.4)
C	parameter (hscin_2y_zpos = 280.9)   !Original Hodo2y C.Y. 07/17/2018
	parameter (hscin_2y_zpos = 282.4)
	parameter (hscin_1x_thick = 1.000*1.067)
	parameter (hscin_1y_thick = 1.000*1.067)
	parameter (hscin_2x_thick = 1.000*1.067)
	parameter (hscin_2y_thick = 1.000*1.067)
C  	parameter (hscin_1x_left  =  45.)      !Original (SOS) Hodo1x left C.Y. 07/17/2018
	parameter (hscin_1x_left  =  50.)
C	parameter (hscin_1x_right = -45.)      !Original (SOS) Hodo1x right C.Y. 07/17/2018
	parameter (hscin_1x_right = -50.)
	parameter (hscin_1x_offset = 0.0)
	parameter (hscin_1y_top   = -45.)
	parameter (hscin_1y_bot   =  45.)
	parameter (hscin_1y_offset = 0.0)
	parameter (hscin_2x_left  =   55.)
	parameter (hscin_2x_right =  -55.)
	parameter (hscin_2x_offset =  0.)
C	parameter (hscin_2y_top   = -57.5)     !Original (SOS) Hodo2Y TOP C.Y. 07/17/2018
	parameter (hscin_2y_top   = -62.5)
C	parameter (hscin_2y_bot   =  57.5)     !Original (SOS) Hodo2Y BOT C.Y. 07/17/2018
	parameter (hscin_2y_bot   =  62.5)
	parameter (hscin_2y_offset = 0)


C Cherenkov position

	real*8 hcer_1_zentrance,hcer_1_zmirror,hcer_1_zexit
	real*8 hcer_2_zentrance,hcer_2_zmirror,hcer_2_zexit
	parameter (hcer_1_zentrance = -291.700)
	parameter (hcer_1_zmirror   = -84.900)
	parameter (hcer_1_zexit     = -61.700)
	parameter (hcer_2_zentrance = 72.600)
	parameter (hcer_2_zmirror   = 179.400)
	parameter (hcer_2_zexit     = 202.600)
	parameter (hcer_left     =  45.00)
	parameter (hcer_right    = -45.00)
	parameter (hcer_top      = -55.00)
	parameter (hcer_bottom   =  55.00)


C Calorimeter position

	real*8 hcal_1pr_zpos,hcal_2ta_zpos,hcal_3ta_zpos,hcal_4ta_zpos
	real*8 hcal_left,hcal_right,hcal_top,hcal_bottom
	parameter (hcal_1pr_zpos = 292.6)
	parameter (hcal_2ta_zpos = 306.4)
	parameter (hcal_3ta_zpos = 323.7)
	parameter (hcal_4ta_zpos = 341.0)
	parameter (hcal_left     =  63.00)
	parameter (hcal_right    = -63.00)
	parameter (hcal_top      = -70.00)
	parameter (hcal_bottom   =  70.00)
