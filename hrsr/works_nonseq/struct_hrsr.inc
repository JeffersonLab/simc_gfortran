! STRUCT_HRSR.INC

! Define some record structures

! ... where do we loose events?
	structure /STOP_EVENTS_HRSR/
		integer*4	trials,slit_hor,slit_vert,
     >				Q1_in,Q1_mid,Q1_out,
     >				Q2_in,Q2_mid,Q2_out,
     >				Q3_in,Q3_mid,Q3_out,
     >				D1_in,D1_out,
     >				hut,dc1,dc2,s1,s2,cal,successes
	end structure

! ... a general common block
	common /spectrometer3/	rSTOP

	record /STOP_EVENTS_HRSR/	rSTOP
