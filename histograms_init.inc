! HISTOGRAMS.INC
! Storage space for the SIMULATE histogram arrays

! Record structures that are only needed for histograms

	integer		nHbins
	parameter	(nHbins=50)
	type hist_entry
	  real*8	bin,min
!	  real*8	buf(nHbins)/nHbins*0./   !removed for Linux-compatibility
	  real*8	buf(nHbins)
	end type

	type hist_arm
	    type(hist_entry):: delta, yptar, xptar
        end type

	type hist_arm2
	    type(hist_entry):: delta, yptar, xptar
        end type

	type hist_double_arm
	  type(hist_arm)::e
	  type(hist_arm2)::p
	  type (hist_entry):: Em, Pm
	end type

	type histograms
	  type (hist_double_arm):: RECON, gen, geni
	end type
