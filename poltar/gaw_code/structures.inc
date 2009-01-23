! STRUCTURES.INC 
!   (note that all initializations of fields have been moved
!    into the file "structures_init.inc")

! Define some BASIC record structures and associated parameters

! ... generic cut -->  initialized with MAX large and MIN small
	structure /cuts/ 
		real*8		min, max
	end structure

! ... generic range (rather than a cut) --> initialized with HI small and LO large
	structure /range/
		real*8		lo, hi
	end structure

! ... generic Cartesian vector
	integer*4 nCartesianfields
	parameter (nCartesianfields = 3)
	structure /Cartesian/
		real*8		x,y,z
	end structure

! ... minimal description of scattered particles, or cuts/ranges on these qties
	structure /double_arm/
		structure /arm/ e
		  real*8		delta, yptar, xptar, z
		end structure
		structure /arm2/ p
		  real*8		delta, yptar, xptar, z
		end structure
	end structure
	structure /double_arm_cuts/
		structure /arm_cuts/ e
		  record /cuts/	delta, yptar, xptar, z
		end structure
		structure /arm_cuts2/ p
		  record /cuts/	delta, yptar, xptar, z
		end structure
	end structure
	structure /double_arm_range/
		structure /arm_range/ e
		  record /range/	delta, yptar, xptar, z
		end structure
		structure /arm_range2/ p
		  record /range/	delta, yptar, xptar, z
		end structure
	end structure

! ... generic focal plane vectors (transport convention, in both spectrometers)
	structure /double_arm_FP/
		structure /arm_FP/ e
		  real*8		x, dx, y, dy
		end structure
		structure /arm_FP2/ p
		  real*8		x, dx, y, dy
		end structure
	end structure

! ... full description of a given particle
	integer*4 narm_fullfields
	parameter (narm_fullfields = 8)
	structure /arm_full/
		real*8	delta,yptar,xptar,z
		real*8	theta,phi,E,P
	end structure

	structure /cuts_true/
		record /cuts/	Em, Pm
	end structure

! EVENT structures

! ... description of event -- both actual and reconstructed (measured) versions calculated

! NEVENTFIELDS is used to copy from one /event/ record into another.
! IF YOU MODIFY THIS STRUCTURE, YOU MUST MAKE SURETHAT NEVENTFIELDS
! IS UPDATED ACCORDINGLY, OR BAD THINGS CAN HAPPEN.

	integer*4 neventfields
	parameter (neventfields = 17 + 2*narm_fullfields + 3*nCartesianfields)
	structure /event/
		union
		  map
			real*8	Ein
			real*8	Em, Pm, Pmx, Pmy, Pmz
			real*8	PmPar, PmPer, PmOop
			real*8	omega, q, Q2, Trec, W
			real*8  epsilon, theta_pq, phi_pq
			record /arm_full/ e, p
			record /Cartesian/ ue, up, uq
		  end map
		  map
			real*8		all(neventfields)
		  end map
		end union
	end structure

! ... target-specific event quantities
	structure /event_target/
		real*8	x, y, z, rastery, rasterx
		real*8	teff(3), Eloss(3), Coulomb
	end structure

! ... quantities that are determined only once per event
	structure /event_main/
		real*8 weight, SF_weight, gen_weight, jacobian
		real*8 Ein_shift, Ee_shift
		real*8 sigcc, sigcc_recon
		real*8 epsilon,theta_pq,phi_pq
		real*8 t,pcm,thetacm,phicm,w,q2
		real*8 davejac,johnjac
		record /event_target/ target
		record /double_arm/ SP, RECON
		record /double_arm_FP/ FP
! ........ temporary ???
		real*8 Trec
		logical accept
	end structure

! ... a gross structure that serves merely to catch all interesting qties for a 'central' event
	structure /event_central/
		real*8		sigcc, omega, q, Q2
		real*8		Em, Pm
		structure /event_central_rad/ rad
		  real*8 hardcorfac, etatzai, frac(3), lambda(3), bt(2)
		  real*8 c_int(0:3), c_ext(0:3), c(0:3), g_int, g_ext, g(0:3)
		end structure
	end structure

! OTHER stuff

! ... spectrometer settings and specifications
	structure /both_spec/
		structure /spectrometer/ e
		  real*8	P,theta,cos_th,sin_th,phi
		  structure /spec_offset/ offset
		    real*8	x,y,z,yptar,xptar
		  end structure
		end structure
		structure /spectrometer2/ p
		  real*8	P,theta,cos_th,sin_th,phi
		  structure /spec_offset2/ offset
		    real*8	x,y,z,yptar,xptar
		  end structure
		end structure
	end structure

! ... acceptance edges for TRUE and VERTEX quantities, both BEFORE reconstruction
	structure /edge_true/
		structure /edge_arm/ e
		  record /cuts/		E, yptar, xptar
		end structure
		structure /edge_arm2/ p
		  record /cuts/		E, yptar, xptar
		end structure
		record /cuts/		Em, Pm, Trec
	end structure

! ... pieces of the EXP dbase field that we'll need
	structure /EXP_field/
		real*8	charge
	end structure

! ... generic description of a histogram axis
	structure /axis/
		real*8		min,max,bin
		integer		n
	end structure

! ... ranges for the quantities that are generated randomly for each event / edges on quantities at the GENERATE stage
	structure /gen_limits/
		structure /arm_limits/ e
		  record /cuts/	delta, yptar, xptar, E
		end structure
		structure /arm_limits2/ p
		  record /cuts/	delta, yptar, xptar, E
		end structure
		record /cuts/	sumEgen, Trec
		real*8		xwid, ywid
	end structure

! ... ranges of event qties which actually contributed
	structure /contrib/
		structure /contrib_gen/ gen
		  record /arm_range/ e, p
		  record /range/ Trec, sumEgen
		end structure
		structure /contrib_true/ tru
		  structure /contrib_arm/ e
		    record /range/ E, yptar, xptar
		  end structure
		  structure /contrib_arm2/ p
		    record /range/ E, yptar, xptar
		  end structure
		  record /range/ Em, Pm, Trec
		end structure
		record /double_arm_range/ SP
		structure /contrib_vertex/ vertex
		  record /range/ Trec, Em, Pm
		end structure
		structure /contrib_rad/	rad
		  record /range/ Egamma(3), Egamma_total
		end structure
	end structure

! ... values, and ranges of values which actually contributed, for useful slops (some are local to limits_init)
	structure /slop_item/
		real*8		lo, hi, used
	end structure
	structure /slop/
		structure /slop_total/ total
		  record /slop_item/	Em, Pm
		end structure
		structure /slop_MC/ MC
		  structure /slop_MC_arm/ e
		    record /slop_item/ delta, yptar, xptar
		  end structure
		  structure /slop_MC_arm2/ p
		    record /slop_item/ delta, yptar, xptar
		  end structure
		end structure
	end structure
