! STRUCT_CALO.INC

! Define some record structures

! ... where do we loose events?
	integer*4	caloSTOP_trials,caloSTOP_slit_hor,caloSTOP_slit_vert,caloSTOP_successes

! ... a general common block
	common /spectrometer_calo/  caloSTOP_trials,caloSTOP_slit_hor
     > ,caloSTOP_slit_vert,caloSTOP_successes

