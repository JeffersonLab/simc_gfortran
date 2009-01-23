C Aperture Limits:

	real*8 r2_quad,w_exit,b_exit,t_exit
	real*8 w_bm01,t_bm01_in,t_bm01_out,b_bm01_in,b_bm01_out
	real*8 w_bm02,t_bm02_in,t_bm02_out,b_bm02_in,b_bm02_out

	parameter       (r2_quad	=163.84)        !12.7 cm radius^2
	parameter	(w_bm01		=  8.0)		!BM01 half-gap in cm
	parameter	(w_bm02		=  8.0)		!BM02 half-gap in cm
	parameter	(t_bm01_in	= 41.32)	!BM01-IN YP top limit
	parameter	(b_bm01_in	=-30.52)	!BM01-IN YP bot limit
	parameter	(t_bm01_out	= 53.05)	!BM01-OUT YP top limit
	parameter	(b_bm01_out	=-65.07)	!BM01-OUT YP bot limit
	parameter	(t_bm02_in	= 51.73)	!BM02-IN YP top limit
	parameter	(b_bm02_in	=-66.38)	!BM02-IN YP bot limit
	parameter	(t_bm02_out	= 51.52)	!BM02-OUT YP top limit
	parameter	(b_bm02_out	=-55.85)	!BM02-OUT YP bot limit
	parameter	(w_exit		=  8.57)	!Horiz 1/2-gap exit flng
	parameter	(t_exit		= 51.52)	!exit flange top limit
	parameter	(b_exit		=-55.85) 	!exit flange bot limit
