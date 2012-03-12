! NEW, NON-STRUCTURE VERSION TO MAXIMIZE MACHINE INDEPENDENCE FOR
! STAND ALONE CODES.

	integer*4 shmsSTOP_trials,shmsSTOP_successes
	integer*4 shmsSTOP_targ_hor,shmsSTOP_targ_vert,shmsSTOP_targ_oct
	integer*4 shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct
        integer*4 shmsSTOP_HB_in,shmsSTOP_HB_men,shmsSTOP_HB_mex
        integer*4 shmsSTOP_HB_out,shmsSTOP_Q1_in,shmsSTOP_Q1_men
	integer*4 shmsSTOP_Q1_mid,shmsSTOP_Q1_mex,shmsSTOP_Q1_out
        integer*4 shmsSTOP_Q2_in,shmsSTOP_Q2_men
	integer*4 shmsSTOP_Q2_mid,shmsSTOP_Q2_mex,shmsSTOP_Q2_out
        integer*4 shmsSTOP_Q3_in,shmsSTOP_Q3_men
	integer*4 shmsSTOP_Q3_mid,shmsSTOP_Q3_mex,shmsSTOP_Q3_out
	integer*4 shmsSTOP_D1_in,shmsSTOP_D1_flr
        integer*4 shmsSTOP_D1_men,shmsSTOP_D1_mid1
	integer*4 shmsSTOP_D1_mid2,shmsSTOP_D1_mid3,shmsSTOP_D1_mid4
	integer*4 shmsSTOP_D1_mid5,shmsSTOP_D1_mid6,shmsSTOP_D1_mid7
        integer*4 shmsSTOP_D1_mex,shmsSTOP_D1_out
	integer*4 shmsSTOP_BP_in,shmsSTOP_BP_out
	integer*4 shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2
	integer*4 shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal
	integer*4 shmsSTOP_cer1,shmsSTOP_cer2
	integer*4 stop_id

	common/spectrometer_shms/ shmsSTOP_trials,shmsSTOP_successes,
     >		shmsSTOP_targ_hor,shmsSTOP_targ_vert,shmsSTOP_targ_oct,
     >		shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct,
     >          shmsSTOP_HB_in,shmsSTOP_HB_men,shmsSTOP_HB_mex,
     >          shmsSTOP_HB_out,shmsSTOP_Q1_in,shmsSTOP_Q1_men,
     >		shmsSTOP_Q1_mid,shmsSTOP_Q1_mex,shmsSTOP_Q1_out,
     >		shmsSTOP_Q2_in,shmsSTOP_Q2_men,shmsSTOP_Q2_mid,
     >		shmsSTOP_Q2_mex,shmsSTOP_Q2_out,shmsSTOP_Q3_in,
     >          shmsSTOP_Q3_men,shmsSTOP_Q3_mid,shmsSTOP_Q3_mex,
     >          shmsSTOP_Q3_out,shmsSTOP_D1_in,shmsSTOP_D1_flr,
     >          shmsSTOP_D1_men,
     >		shmsSTOP_D1_mid1,shmsSTOP_D1_mid2,shmsSTOP_D1_mid3,
     >          shmsSTOP_D1_mid4,shmsSTOP_D1_mid5,shmsSTOP_D1_mid6,
     >          shmsSTOP_D1_mid7,shmsSTOP_D1_mex,shmsSTOP_D1_out,
     >		shmsSTOP_BP_in,shmsSTOP_BP_out,
     >		shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2,
     >		shmsSTOP_cer1,shmsSTOP_cer2,
     >		shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal,stop_id
