! STRUCT_SOS.INC
!
! Define some record structures
!
! ... where do we loose events?
!	structure /STOP_EVENTS_SOS/
!	  integer*4	trials,slit_hor,slit_vert,slit_oct
!	  integer*4	quad_in,quad_mid,quad_out
!	  integer*4	bm01_in,bm01_out,bm02_in,bm02_out,exit
!	  integer*4	hut,dc1,dc2,scin,cal,successes
!	end structure
!
! ... a general common block
!	common /spectrometer2/	sSTOP
!
!	record /STOP_EVENTS_SOS/	sSTOP


! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 sSTOP_trials
	integer*4 sSTOP_slit_hor,sSTOP_slit_vert,sSTOP_slit_oct
	integer*4 sSTOP_quad_in,sSTOP_quad_mid,sSTOP_quad_out
	integer*4 sSTOP_bm01_in,sSTOP_bm01_out
	integer*4 sSTOP_bm02_in,sSTOP_bm02_out
	integer*4 sSTOP_exit,sSTOP_hut
	integer*4 sSTOP_dc1,sSTOP_dc2,sSTOP_scin,sSTOP_cal
	integer*4 sSTOP_successes

	common/spectrometer2/ sSTOP_trials,
     >		sSTOP_slit_hor,sSTOP_slit_vert,sSTOP_slit_oct,
     >		sSTOP_quad_in,sSTOP_quad_mid,sSTOP_quad_out,
     >		sSTOP_bm01_in,sSTOP_bm01_out,
     >		sSTOP_bm02_in,sSTOP_bm02_out,
     >		sSTOP_exit,sSTOP_hut,
     >		sSTOP_dc1,sSTOP_dc2,sSTOP_scin,sSTOP_cal,
     >		sSTOP_successes

