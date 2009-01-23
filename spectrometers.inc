C flag for forward MEs telling if a transformation is drift-only, and giving
C the extracted drift length.

	integer*4       max_class               !maximum number of classes.
	parameter       (max_class = 18)        !must also change in transp.f
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

        common  /track/ xs,ys,zs,dpps,dxdzs,dydzs



C Version for non-sequential transformations.

C       real*8 x_transp,y_transp        !positions at pivot (cm)
C       real*8 dxdz_transp,dydz_transp  !slopes at pivot. Input to transp.f
C       common  /track_seq/ x_transp,y_transp,dxdz_transp,dydz_transp

C x_transp,y_transp,dxdz_transp,dydz_transp are the positions(cm) and slopes
C (dx/dz and dy/dz) at the pivot.  They are the input to transp.f which
C uses the COSY matrix elements to project the track thru the magnets.
C Since dpp never changes, transp.f just uses dpps.
