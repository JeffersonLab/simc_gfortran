C flag for forward MEs telling if a transformation is drift-only, and giving
C the extracted drift length.

	integer*4       max_class               !maximum number of classes.
	parameter       (max_class = 41)        !must also change in transp.f
	integer*4       nspectr
	parameter       (nspectr=6)

	real*8 driftdist(nspectr,max_class)
	logical adrift(nspectr,max_class)

	common /driftdistances/ driftdist,adrift


C Common block with TEMPORARY track variables (the same variables are
C used by the electron and hadron arm, so are overwritten by 2nd spectrometer).

C xs,ys,zs,dxdzs,dydzs are positions (cm) and slopes (dx/dz and dy/dz) for
C the track.  Initially, you start with the pivot values.  As you go thru
C the magnets, these are the values output by transp.f.


        real*8 xs,ys,zs                 !track x,y,z positions (cm)
        real*8 dxdzs,dydzs              !track slop (dx/dz, dy/dz)
        real*8 dpps                     !delta dp/p (%).
        real*8 xs_fp,ys_fp                 !track x,y,z positions (cm) before det. hut
        real*8 dxdzs_fp,dydzs_fp              !track slop (dx/dz, dy/dz) before det. hut

	real*8 x_tgt_out
	real*8 y_tgt_out
	real*8 x_hb_in,x_hb_men,x_hb_mex,x_hb_out
	real*8 y_hb_in,y_hb_men,y_hb_mex,y_hb_out
	real*8 x_q1_in,x_q1_men,x_q1_mid,x_q1_mex,x_q1_out
	real*8 y_q1_in,y_q1_men,y_q1_mid,y_q1_mex,y_q1_out
	real*8 x_q2_in,x_q2_men,x_q2_mid,x_q2_mex,x_q2_out
	real*8 y_q2_in,y_q2_men,y_q2_mid,y_q2_mex,y_q2_out
	real*8 x_q3_in,x_q3_men,x_q3_mid,x_q3_mex,x_q3_out
	real*8 y_q3_in,y_q3_men,y_q3_mid,y_q3_mex,y_q3_out
	real*8 x_d_in,x_d_flr,x_d_men,x_d_m1,x_d_m2,x_d_m3
        real*8 x_d_m4
	real*8 y_d_in,y_d_flr,y_d_men,y_d_m1,y_d_m2,y_d_m3
        real*8 y_d_m4
	real*8 x_d_m5,x_d_m6,x_d_m7,x_d_mex,x_d_out
	real*8 y_d_m5,y_d_m6,y_d_m7,y_d_mex,y_d_out
	real*8 ev_lost

        common  /track/ xs,ys,zs,dpps,dxdzs,dydzs,x_hb_in,x_hb_men,
     >      x_hb_mex,x_hb_out,y_hb_in,y_hb_men,y_hb_mex,y_hb_out,
     >      x_tgt_out,y_tgt_out,x_q1_in,y_q1_in,x_q1_men,y_q1_men,
     >      x_q1_mid,y_q1_mid,x_q1_mex,y_q1_mex,x_q1_out,y_q1_out,
     >      x_q2_in,y_q2_in,x_q2_men,y_q2_men,x_q2_mid,y_q2_mid,
     >      x_q2_mex,y_q2_mex,x_q2_out,y_q2_out,x_q3_in,y_q3_in,
     >      x_q3_mex,y_q3_mex,x_q3_mid,y_q3_mid,x_q3_out,y_q3_out,
     >      x_d_in,y_d_in,x_d_flr,y_d_flr,x_d_men,y_d_men,x_d_out,
     >      y_d_out,x_d_m1,y_d_m1,x_d_m2,y_d_m2,
     >      x_d_m3,y_d_m3,x_d_m4,y_d_m4,
     >      x_d_m5,y_d_m5,x_d_m6,y_d_m6,
     >      x_d_m7,y_d_m7,x_d_mex,y_d_mex,
     >      xs_fp,ys_fp,dxdzs_fp,dydzs_fp,
     >      ev_lost



C Version for non-sequential transformations.

C       real*8 x_transp,y_transp        !positions at pivot (cm)
C       real*8 dxdz_transp,dydz_transp  !slopes at pivot. Input to transp.f
C       common  /track_seq/ x_transp,y_transp,dxdz_transp,dydz_transp

C x_transp,y_transp,dxdz_transp,dydz_transp are the positions(cm) and slopes
C (dx/dz and dy/dz) at the pivot.  They are the input to transp.f which
C uses the COSY matrix elements to project the track thru the magnets.
C Since dpp never changes, transp.f just uses dpps.
