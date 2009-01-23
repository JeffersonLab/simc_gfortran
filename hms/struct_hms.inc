! STRUCT_HMS.INC
!
! Define some record structures
!
! ... where do we loose events?
!	structure /STOP_EVENTS_HMS/
!		integer*4	trials,slit_hor,slit_vert,slit_oct,
!     >				Q1_in,Q1_mid,Q1_out,
!     >				Q2_in,Q2_mid,Q2_out,
!     >				Q3_in,Q3_mid,Q3_out,
!     >				D1_in,D1_out,
!     >				hut,dc1,dc2,scin,cal,successes
!	end structure
!
! ... a general common block
!	common /spectrometer1/ hSTOP
!
!	record /STOP_EVENTS_HMS/ hSTOP


! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 hSTOP_trials
	integer*4 hSTOP_slit_hor,hSTOP_slit_vert,hSTOP_slit_oct
	integer*4 hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out
	integer*4 hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out
	integer*4 hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out
	integer*4 hSTOP_D1_in,hSTOP_D1_out
	integer*4 hSTOP_hut,hSTOP_dc1,hSTOP_dc2
	integer*4 hSTOP_scin,hSTOP_cal,hSTOP_successes

	common/spectrometer1/ 	hSTOP_trials,
     >		hSTOP_slit_hor,hSTOP_slit_vert,hSTOP_slit_oct,
     >		hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out,
     >		hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,
     >		hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out,
     >		hSTOP_D1_in,hSTOP_D1_out,
     >		hSTOP_hut,hSTOP_dc1,hSTOP_dc2,
     >		hSTOP_scin,hSTOP_cal,hSTOP_successes
