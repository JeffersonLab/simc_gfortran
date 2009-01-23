! STRUCT_HRSR.INC
!
! Define some record structures
!
! ... where do we loose events?
!	structure /STOP_EVENTS_HRSR/
!		integer*4	trials,slit_hor,slit_vert,
!     >				Q1_in,Q1_mid,Q1_out,
!     >				Q2_in,Q2_mid,Q2_out,
!     >				Q3_in,Q3_mid,Q3_out,
!     >				D1_in,D1_out,
!     >				hut,dc1,dc2,s1,s2,cal,successes
!	end structure
!
! ... a general common block
!	common /spectrometer3/	rSTOP
!
!	record /STOP_EVENTS_HRSR/	rSTOP


! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 rSTOP_trials
	integer*4 rSTOP_slit_hor,rSTOP_slit_vert
	integer*4 rSTOP_Q1_in,rSTOP_Q1_mid,rSTOP_Q1_out
	integer*4 rSTOP_Q2_in,rSTOP_Q2_mid,rSTOP_Q2_out
	integer*4 rSTOP_Q3_in,rSTOP_Q3_mid,rSTOP_Q3_out
	integer*4 rSTOP_D1_in,rSTOP_D1_out
	integer*4 rSTOP_hut,rSTOP_dc1,rSTOP_dc2
	integer*4 rSTOP_s1,rSTOP_s2,rSTOP_cal,rSTOP_successes

	common/spectrometer3/ rSTOP_trials,
     >		 rSTOP_slit_hor,rSTOP_slit_vert,
     >		 rSTOP_Q1_in,rSTOP_Q1_mid,rSTOP_Q1_out,
     >		 rSTOP_Q2_in,rSTOP_Q2_mid,rSTOP_Q2_out,
     >		 rSTOP_Q3_in,rSTOP_Q3_mid,rSTOP_Q3_out,
     >		 rSTOP_D1_in,rSTOP_D1_out,
     >		 rSTOP_hut,rSTOP_dc1,rSTOP_dc2,
     >		 rSTOP_s1,rSTOP_s2,rSTOP_cal,rSTOP_successes
