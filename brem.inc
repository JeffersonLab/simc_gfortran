! BREM.INC
!
! Series of switches used by BREM* subroutines. 
! Switches should be set initially by calling BREM_INIT.

	logical produce_output		! produce screen output?
	logical exponentiate		! exponentiate soft photon prob?
	logical include_hard		! include hard correction?
	logical calculate_spence	! calculate spence functions?

	common /bremcom/produce_output, exponentiate,
     >		include_hard, calculate_spence
