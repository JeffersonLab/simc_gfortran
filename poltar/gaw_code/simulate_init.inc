! SIMULATE_INIT.INC

! Note: All angles are in radians
!	All distances are in cm
!	All energies (momenta,masses) are in MeV
!	All deltas are in percent
!	All densities (thicknesses) are in g/cm3 (g/cm2)

! Define some record structures, and constants

	include 'structures_init.inc'
	include 'constants.inc'

! Now for the actual variables

! ??? these montecarlo-specific values should go in a DAT file.
! ... the slop that will be allowed on the cuts when they are to be applied "loosely"
	real*8 slop_param_d_HMS, slop_param_t_HMS, slop_param_p_HMS
	parameter	(slop_param_d_HMS=0.5)
	parameter	(slop_param_t_HMS=0.005)
	parameter	(slop_param_p_HMS=0.005)
	real*8 slop_param_d_SOS, slop_param_t_SOS, slop_param_p_SOS
	parameter	(slop_param_d_SOS=1.0)
	parameter	(slop_param_t_SOS=0.008)
	parameter	(slop_param_p_SOS=0.008)

! ... the common blocks

!RMM changed ordering in /gnrl/ to fix memory alignment problems

	real*8	Mh, Mh2, Mh2_final	!mh,mh2 are INITIAL hadron mass(**2)
	real*8	Ebeam, dEbeam, Ebeam_vertex_ave, genvol, genvol_inclusive
	real*8 	luminosity, targetfac, normfac
	real*8 	wtcontribute, dE_edge_test, Egamma_gen_max

	record /gen_limits/	gen
	record /both_spec/	spec
	record /cuts_true/	cuts
	record /edge_true/	edge
	record /edge_true/	VERTEXedge
	record /double_arm_cuts/ SPedge
	record /slop/		slop
	record /EXP_field/	EXPER
	logical*4 mc_smear
	logical*4 debug(6)
	integer*4 nevent, ntried, ngen, deForest_flag, Nntu
	integer*4 ncontribute,ncontribute_no_rad_proton,spect_mode, phsp_mode
	character base*40
	character extra_dbase_file*60
	character tgt_field_file*60
	logical using_E_arm_montecarlo, using_P_arm_montecarlo
	logical doing_phsp, using_rad
	logical doing_hyd_elast, doing_deuterium, doing_heavy
	logical doing_eep, doing_pion, doing_piplus, doing_kaon
	integer*4 which_kaon
	logical using_cit_generation, using_Coulomb, using_Eloss
	logical correct_Eloss, correct_raster
	logical using_pm_sign, hms_e_flag
	logical using_tgt_field
	real*8 drift_to_ndet

	common /gnrl/   Mh, Mh2, Mh2_final, Ebeam, dEbeam, Ebeam_vertex_ave,
	1 genvol, genvol_inclusive, luminosity, targetfac,
	3 normfac, wtcontribute, dE_edge_test, Egamma_gen_max,
	4 gen, spec, cuts, edge, SPedge, VERTEXedge,
	5 slop, EXPER, mc_smear, debug, nevent, ntried, ngen,
	6 deForest_flag, Nntu, ncontribute, ncontribute_no_rad_proton, base,
	7 extra_dbase_file, using_E_arm_montecarlo, using_P_arm_montecarlo,
	8 doing_phsp, using_rad, doing_hyd_elast, doing_deuterium, doing_heavy,
	9 doing_eep, doing_pion, doing_piplus, doing_kaon, which_kaon,
	1 using_cit_generation, using_Coulomb, using_Eloss,
	2 correct_Eloss, correct_raster,
	3 using_pm_sign, hms_e_flag,spect_mode, phsp_mode,
	4 using_tgt_field, tgt_field_file, drift_to_ndet


! ........ note: make these first two parameters at least ONE BIGGER than the actual dimensions you want to read in
	integer*4 ntheorymax,nrhoPmmax
	parameter	(ntheorymax=500)
	parameter	(nrhoPmmax=21)
	real*8		theory, Em_theory, bs_norm_theory, nprot_theory
	real*8		Emsig_theory, Em_int_theory, E_Fermi
	integer		nrhoPm
	character	theory_file*80, theory_base*40, theory_target*40

	record /axis/	Pm_theory

	common /theory/ theory(nrhoPmmax,ntheorymax),
	1    Em_theory(nrhoPmmax), bs_norm_theory(nrhoPmmax),
	2    nprot_theory(nrhoPmmax), Emsig_theory(nrhoPmmax),
	4    Em_int_theory(nrhoPmmax), E_Fermi,
	5    Pm_theory(nrhoPmmax),
	6    theory_file, theory_base, theory_target,
	7    nrhoPm

C decdist(30) is just an array of miscellaneous shared variables, SOME of
C which are related to decay.  Here are the ones that are currently in use:
C	decdist(30)    = position where event decays (cm)
C	decdist(3)  = sum of sigcc (to calculate "central sigcc")
C	decdist(4)  = radphot   (for ntuple output)
C	decdist(5)  = flag for type of radiation (which particle).
C	decdist(6)  = resfac in ntup. (sum of resolution modifiers for HMS/SOS)
C	decdist(21) = sigma_eep (which is eepi cross section - Sigcc in ntup.)

	real*8 decdist(30)
	real*8 ctau
	logical doing_hydpi, doing_deutpi, doing_hepi
	logical doing_hydkaon, doing_deutkaon, doing_hekaon
	logical doing_decay

	common /decd/ decdist,ctau,doing_decay,
	1 doing_hydpi,doing_deutpi,doing_hepi,
	2 doing_hydkaon,doing_deutkaon,doing_hekaon



C Kaon cross section lookup tables (real*4 for CERNLIB fint compatability).
	real*4 zrff1(10,11,19),zrff2(10,11,19),zrff3(10,11,19)
	real*4 zrff4(10,11,19),zrff5(10,11,19),zrff6(10,11,19)
	real*4 ziff1(10,11,19),ziff2(10,11,19),ziff3(10,11,19)
	real*4 ziff4(10,11,19),ziff5(10,11,19),ziff6(10,11,19)
	real*4 zsrff1(20,10,19),zsrff2(20,10,19),zsrff3(20,10,19)
	real*4 zsrff4(20,10,19),zsrff5(20,10,19),zsrff6(20,10,19)
	real*4 zsiff1(20,10,19),zsiff2(20,10,19),zsiff3(20,10,19)
	real*4 zsiff4(20,10,19),zsiff5(20,10,19),zsiff6(20,10,19)
	real*4 thrown(50,50,50)
	real*4 genera(50,50,50)
	real*4 recons(50,50,50)
	real*4 weightc(20,50)
	real*4 weightd(8,40,30)

	common /sigkaon/ thrown,genera,recons,weightc,weightd,
	1 zrff1,zrff2,zrff3,zrff4,zrff5,zrff6,
	2 ziff1,ziff2,ziff3,ziff4,ziff5,ziff6,
	3 zsrff1,zsrff2,zsrff3,zsrff4,zsrff5,zsrff6,
	4 zsiff1,zsiff2,zsiff3,zsiff4,zsiff5,zsiff6



C momentum distribution (nump=#/points, pval=p, mprob=integral (p^2 f(p) dp)
	real*8 pval(2000),mprob(2000)
	real*8 pfer,pferx,pfery,pferz
	integer*4 nump

	common /pfermi/ pval,mprob,pfer,pferx,pfery,pferz,nump

! ... and some variable blocks of general interest defined in external files
	include 'target.inc'
