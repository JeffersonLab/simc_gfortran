! RADC.INC

	integer intcor_mode
	integer one_tail, use_expon, rad_flag, extrad_flag, ntail
	real*8 fwhm, etatzai, etta, frac(3), lambda(3), bt(2), hardcorfac
	real*8 c_int(0:3), c_ext(0:3), c(0:4), g_int, g_ext, g(0:4)
	real*8 Egamma_used(3), Egamma_min(3), Egamma_max(3)
	real*8 Egamma_res_limit
	real*8 Egamma_tot_max, Egamma1_max, Egamma2_max, Egamma3_max
	logical doing_tail(3)
	logical use_offshell_rad, rad_proton_this_ev, hardwired_rad

	common /radccom/	fwhm,etatzai,etta,frac,lambda,bt,hardcorfac,
     >  	c_int,c_ext,c,g_int,g_ext,g,Egamma_used,Egamma_min,
     >		Egamma_max,Egamma_res_limit,
     >  	Egamma_tot_max,Egamma1_max,Egamma2_max,Egamma3_max,
     >  	intcor_mode, one_tail, use_expon,
     >  	rad_flag,extrad_flag,ntail,doing_tail,
     >  	use_offshell_rad,rad_proton_this_ev,hardwired_rad

	real*8 lim1,lim2,lim3
	real*8 phot1,phot2,phot3
	common /radctest/lim1,lim2,lim3,phot1,phot2,phot3
