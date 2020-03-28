! SIMULATE.INC

! Note: All angles are in radians
!	All distances are in cm
!	All energies (momenta,masses) are in MeV
!	All deltas are in percent
!	All densities (thicknesses) are in g/cm3 (g/cm2)

! Define some record structures, and constants

	include 'structures.inc'
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
	real*8 slop_param_d_HRSR, slop_param_t_HRSR, slop_param_p_HRSR
	parameter	(slop_param_d_HRSR=0.5)
	parameter	(slop_param_t_HRSR=0.005)
	parameter	(slop_param_p_HRSR=0.005)
	real*8 slop_param_d_HRSL, slop_param_t_HRSL, slop_param_p_HRSL
	parameter	(slop_param_d_HRSL=0.5)
	parameter	(slop_param_t_HRSL=0.005)
	parameter	(slop_param_p_HRSL=0.005)
	real*8 slop_param_d_SHMS, slop_param_t_SHMS, slop_param_p_SHMS
	parameter	(slop_param_d_SHMS=0.5)
	parameter	(slop_param_t_SHMS=0.005)
	parameter	(slop_param_p_SHMS=0.005)

! ... the common blocks

!RMM changed ordering in /gnrl/ to fix memory alignment problems

	real*8	Mh, Mh2, Mh2_final	!mh,mh2 are INITIAL hadron mass(**2)
	real*8	Ebeam, dEbeam, Ebeam_vertex_ave, genvol, genvol_inclusive
	real*8	luminosity, targetfac, normfac
	real*8	wtcontribute, dE_edge_test, Egamma_gen_max
	real*8  pt_b_param,sigc_kin_min,sigc_kin_max,sigc_kin_ind
	real*8  targ_Bangle, targ_Bphi, targ_pol, drift_to_cal,sign_hadron
	integer*4 sigc_flag,sigc_nbin

	type(gen_limits)::	gen
	type(both_spec)::	spec
	type(cuts_true)::	cuts
	type(edge_true)::	edge
	type(edge_true)::	VERTEXedge
	type(double_arm_cuts):: SPedge
	type(sloptype)::	slop
	type(EXP_field)::	EXPER
	logical mc_smear
	logical debug(6)
	integer*4 nevent, ntried, ngen, deForest_flag, Nntu
	integer*4 ncontribute, npasscuts, ncontribute_no_rad_proton, spect_mode, phsp_mode
	character base*40
	character extra_dbase_file*60
	character tgt_field_file*60
	logical using_E_arm_montecarlo, using_P_arm_montecarlo
	logical doing_phsp, using_rad, hard_cuts
	logical doing_hyd_elast, doing_deuterium, doing_heavy
	logical doing_eep, doing_pion, doing_delta, doing_kaon, doing_rho
	logical doing_semi, doing_semipi, doing_semika, doing_hplus
	integer*4 which_kaon, which_pion
	logical using_cit_generation, using_Coulomb, using_Eloss
	logical correct_Eloss, correct_raster,do_fermi
	logical using_tgt_field
	integer*4 electron_arm, hadron_arm, use_first_cer
	real*8 transparency
	logical use_benhar_sf
	character*80 random_state_file
	character*132 start_random_state_file
	integer*4 random_seed


	common /gnrl/   Mh, Mh2, Mh2_final, Ebeam, dEbeam, Ebeam_vertex_ave,
     >    genvol, genvol_inclusive, luminosity, targetfac,
     >    normfac, wtcontribute, dE_edge_test, Egamma_gen_max,pt_b_param,
     >    sigc_kin_min,sigc_kin_max,sigc_kin_ind,targ_Bangle, targ_Bphi, targ_pol, 
     >    drift_to_cal,sign_hadron,sigc_flag,sigc_nbin,
     >    gen, spec, cuts, edge, VERTEXedge, SPedge,
     >    slop, EXPER, mc_smear, debug, nevent, ntried, ngen,
     >    deForest_flag, Nntu, ncontribute, npasscuts, ncontribute_no_rad_proton, 
     >    spect_mode, phsp_mode, base,
     >    extra_dbase_file,tgt_field_file,using_E_arm_montecarlo,using_P_arm_montecarlo,
     >    doing_phsp, using_rad, hard_cuts, doing_hyd_elast, doing_deuterium, doing_heavy,
     >    doing_eep, doing_pion, doing_delta, doing_kaon, doing_rho, doing_semi,
     >    doing_semipi, doing_semika, doing_hplus, which_kaon,
     >    which_pion, using_cit_generation, using_Coulomb, using_Eloss,
     >    correct_Eloss, correct_raster, do_fermi,using_tgt_field,
     >    electron_arm, hadron_arm, use_first_cer,
     >    transparency, use_benhar_sf,start_random_state_file,
     >    random_state_file, random_seed

! ........ note: make these first two parameters at least ONE BIGGER than the actual dimensions you want to read in
	integer*4 ntheorymax,nrhoPmmax
	parameter	(ntheorymax=500)
	parameter	(nrhoPmmax=21)
	real*8		theory, Em_theory, bs_norm_theory, nprot_theory
	real*8		Emsig_theory, Em_int_theory, E_Fermi
	integer		nrhoPm
	character	theory_file*80, theory_base*40, theory_target*40

	type(axis)::	Pm_theory

	common /theory/		theory(nrhoPmmax,ntheorymax),
     >   	Em_theory(nrhoPmmax), bs_norm_theory(nrhoPmmax),
     >  	nprot_theory(nrhoPmmax), Emsig_theory(nrhoPmmax),
     >  	Em_int_theory(nrhoPmmax), E_Fermi,
     >  	Pm_theory(nrhoPmmax),
     >  	theory_file, theory_base, theory_target,
     >  	nrhoPm

C decdist(30) is just an array of miscellaneous shared variables, SOME of
C which are related to decay.  Here are the ones that are currently in use:
Cxx	decdist(30)    = position where event decays (cm)
C	decdist(3)  = sum of sigcc (to calculate "central sigcc")
Cxx	decdist(4)  = radphot   (for ntuple output)
Cxx	decdist(5)  = flag for type of radiation (which particle).
Cxx	decdist(6)  = resfac in ntup. (sum of resolution modifiers for HMS/SOS)
Cxx	decdist(21) = sigma_eep (which is eepi cross section - Sigcc in ntup.)
C
C	real*8 decdist(30)

	real*8 ctau
	logical doing_decay
	logical doing_hydpi, doing_deutpi, doing_hepi
	logical doing_hydkaon, doing_deutkaon, doing_hekaon
	logical doing_hyddelta, doing_deutdelta, doing_hedelta
	logical doing_hydrho, doing_Deutrho, doing_herho
	logical doing_hydsemi, doing_deutsemi

	common /decd/ ctau,doing_decay,
     >    doing_hydpi,doing_deutpi,doing_hepi,
     >    doing_hydkaon,doing_deutkaon,doing_hekaon,
     >    doing_hyddelta,doing_deutdelta,doing_hedelta,
     >    doing_hydrho, doing_deutrho, doing_herho,
     >    doing_hydsemi, doing_deutsemi

C Trying to get rid of the decdist array of random variables.
C Structure for variables that we want to available to the ntuple.
	type ntupvars
                sequence
		real*8 radphot, radarm
		real*8 resfac
		real*8 sigcm
		real*8 krel,mm,mmA,t,dilu
		real*8 sigcm1, sigcm2, sigcm3, sigcm4
		real*8 xfermi, yfermi
		real*8 rhomass, rhotheta
	end type

	type (ntupvars):: ntup

C decdist is zpos of decay if doing_decay, survivalprobability if
C .not.doing_decay.  NOT included in ntup.* structure because it is used
C in the single arm MC, which we want to keep structure-free.

	real*8 decdist

	common /ntuple_variables/ ntup,decdist



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
     >    zrff1,zrff2,zrff3,zrff4,zrff5,zrff6,
     >    ziff1,ziff2,ziff3,ziff4,ziff5,ziff6,
     >    zsrff1,zsrff2,zsrff3,zsrff4,zsrff5,zsrff6,
     >    zsiff1,zsiff2,zsiff3,zsiff4,zsiff5,zsiff6



C momentum distribution (nump=#/points, pval=p, mprob=integral (p^2 f(p) dp)
	real*8 pval(2000),mprob(2000)
	real*8 efer,pfer,pferx,pfery,pferz
	integer*4 nump
	integer*4 nume
	real*8 eval(2000)

	common /pfermi/ eval,pval,mprob,efer,pfer,pferx,pfery,pferz,nump,nume

! ... and some variable blocks of general interest defined in external files
	include 'target.inc'
