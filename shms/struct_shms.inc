! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 shmsSTOP_trials,shmsSTOP_successes
	integer*4 shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct
	integer*4 shmsSTOP_Q1_in,shmsSTOP_Q1_mid,shmsSTOP_Q1_out
	integer*4 shmsSTOP_Q2_in,shmsSTOP_Q2_mid,shmsSTOP_Q2_out
	integer*4 shmsSTOP_Q3_in,shmsSTOP_Q3_mid
	integer*4 shmsSTOP_D1_in,shmsSTOP_D1_mid1,shmsSTOP_D1_mid2
	integer*4 shmsSTOP_D1_mid3,shmsSTOP_D1_mid4,shmsSTOP_D1_mid5
	integer*4 shmsSTOP_D1_mid6,shmsSTOP_D1_out
	integer*4 shmsSTOP_BP_in,shmsSTOP_BP_out
	integer*4 shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2
	integer*4 shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal

	common/spectrometer_shms/ shmsSTOP_trials,shmsSTOP_successes,
     >		shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct,
     >		shmsSTOP_Q1_in,shmsSTOP_Q1_mid,shmsSTOP_Q1_out,
     >		shmsSTOP_Q2_in,shmsSTOP_Q2_mid,shmsSTOP_Q2_out,
     >		shmsSTOP_Q3_in,shmsSTOP_Q3_mid,
     >		shmsSTOP_D1_in,shmsSTOP_D1_mid1,shmsSTOP_D1_mid2,
     >	        shmsSTOP_D1_mid3,shmsSTOP_D1_mid4,shmsSTOP_D1_mid5,
     >          shmsSTOP_D1_mid6,shmsSTOP_D1_out,
     >		shmsSTOP_BP_in,shmsSTOP_BP_out,
     >		shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2,
     >		shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal
