! TARGET.INC

! Constants

	real*8		rho_Al, Z_Al, A_Al, X0_Al, X0_cm_Al
	parameter	(rho_Al = 2.70)
	parameter	(Z_Al = 13.)
	parameter	(A_Al = 26.98)
	parameter	(X0_Al = 24.01)
	parameter	(X0_cm_Al = X0_Al/rho_Al)
        real*8          rho_mylar, Z_mylar, A_mylar, X0_mylar, X0_cm_mylar
        parameter       (rho_mylar = 1.39)
        parameter       (Z_mylar = 4.545)
        parameter       (A_mylar = 8.735)
        parameter       (X0_mylar = 39.95)
        parameter       (X0_cm_mylar = X0_mylar/rho_mylar)
        real*8          rho_kevlar, Z_kevlar, A_kevlar, X0_kevlar, X0_cm_kevlar
        parameter       (rho_kevlar = 0.74)
        parameter       (Z_kevlar = 2.67)
        parameter       (A_kevlar = 4.67)
        parameter       (X0_kevlar = 55.2)
        parameter       (X0_cm_kevlar = X0_kevlar/rho_kevlar)
	real*8          rho_air, Z_air, A_air, X0_air, X0_cm_air
        parameter       (rho_air = 0.00121)
        parameter       (Z_air = 7.2)   !80% N2, 20% O2
        parameter       (A_air = 14.4)
        parameter       (X0_air = 36.66)
        parameter       (X0_cm_air = X0_air/rho_air)
        real*8		inch_cm, target_pi
	parameter	(inch_cm = 2.54)
! ........ stupid name to avoid conflict with definition in constants.inc
	parameter	(target_pi = 3.14159265358979)

! Target information

	type target_function
                sequence
		real*8		ave, min, max
	end type

	type target_info
                sequence
		real*8	A, Z, N, mass_amu, M, mrec_amu, Mrec, rho, thick, angle,
     >			abundancy, length, zoffset, X0, X0_cm, L1, L2,
     >			fr1, fr2, xoffset, yoffset
		integer fr_pattern,can
		type (target_function):: Eloss(3), teff(3), Coulomb
		real*8	Coulomb_constant, musc_max(3), musc_nsig_max
		real*8  Mtar_struck, Mrec_struck
	end type

	type (target_info)::	targ
	common /target_info/	targ
