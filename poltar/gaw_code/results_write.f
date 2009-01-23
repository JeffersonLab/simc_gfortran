	subroutine results_ntu_write(main,vertex,orig,recon,success)

	implicit none
	include 'radc.inc'
	include 'hbook.inc'
	include 'simulate.inc'

	real*4	ntu(80)
	record /event_main/ main
	record /event/	vertex, orig, recon

!local (e,e'p) calculations:
	real*8 poftheta		!p as calculated from theta, assuming elastic.
	real*8 corrsing		!'corrected singles' for elastic
	real*8 mm,mm2		!missing mass (assuming struck nucleon)
	real*8 mmA,mmA2		!missing mass (assuming struck nucleus)
	real*8 Pm_Heepx,Pm_Heepy,Pm_Heepz	!Pm components for Heepcheck.

!local (e,e'pi/K) calculations:
	real*8 t		!t
	real*8 dummy

	logical success

	if (debug(2)) write(6,*)'r_ntu_write: entering'

!If event did not make it thru spectrometers, return.  Later, we will want to
! add option to write event even if it failed when doing_phsp.

	if (.not.success) return

*	if (phot1.gt.lim1) write(6,*) 'phot1,lim1=',phot1,lim1
*	if (phot2.gt.lim2) then
*	   write(6,*) 'phot2,lim2=',phot2,lim2
*	   write(6,*) egamma_tot_max-egamma_used(1)
*	   write(6,*) vertex.e.E - edge.e.E.min
*	endif
*	if (phot3.gt.lim3) write(6,*) 'phot3,lim3=',phot3,lim3


!Calculate some proton/pion/kaon specific stuff:

	if (doing_pion .or. doing_kaon) then
	  mm2 = recon.Em**2 - recon.Pm**2
	  mm  = sqrt(abs(mm2)) * abs(mm2)/mm2
	  mmA2= (recon.omega + targ.M - recon.p.E)**2 - recon.Pm**2
	  mmA = sqrt(abs(mmA2)) * abs(mmA2)/mmA2
	  t = recon.Q2 - Mh2
     >		+ 2*(recon.omega*recon.p.E - recon.p.P*recon.q*cos(recon.theta_pq))
	endif

	if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  poftheta = Mp*Ebeam / (2*ebeam*sin(recon.e.theta/2.)**2 + Mp)
	  corrsing = recon.e.P - poftheta
	  Pm_Heepz = -(recon.Pmy*recon.uq.y+recon.Pmz*recon.uq.z)
     >		/ sqrt(recon.uq.y**2+recon.uq.z**2)
	  Pm_Heepy =  (recon.Pmz*recon.uq.y-recon.Pmy*recon.uq.z)
     >		/ sqrt(recon.uq.y**2+recon.uq.z**2)
	  Pm_Heepx =  -recon.Pmx
	endif

	if(hms_e_flag)then
	  ntu(1) = recon.e.delta
	  ntu(2) = recon.e.yptar			!mr
	  ntu(3) = recon.e.xptar			!mr
	  ntu(4) = recon.e.z
	  ntu(5) = main.FP.e.x
	  ntu(6) = main.FP.e.dx				!mr
	  ntu(7) = main.FP.e.y
	  ntu(8) = main.FP.e.dy				!mr
	  ntu(9) = orig.e.delta
	  ntu(10) = orig.e.yptar			!mr
	  ntu(11) = orig.e.xptar			!mr
	  ntu(12) = orig.e.z             ! GAW - better analogue to recon ytar
CGAW	  ntu(12) = main.target.z*spec.e.sin_th
	  ntu(13) = recon.p.delta
	  ntu(14) = recon.p.yptar			!mr
	  ntu(15) = recon.p.xptar			!mr
 	  ntu(16) = recon.p.z
 	  ntu(17) = main.FP.p.x
	  ntu(18) = main.FP.p.dx			!mr
 	  ntu(19) = main.FP.p.y
 	  ntu(20) = main.FP.p.dy			!mr
	  ntu(21) = orig.p.delta
	  ntu(22) = orig.p.yptar			!mr
	  ntu(23) = orig.p.xptar			!mr
 	  ntu(24) = orig.p.z             ! GAW - better analogue to recon ytar
CGAW 	  ntu(24) = -main.target.z*spec.p.sin_th
	else
	  ntu(1) = recon.p.delta
	  ntu(2) = recon.p.yptar			!mr
	  ntu(3) = recon.p.xptar			!mr
	  ntu(4) = recon.p.z
	  ntu(5) = main.FP.p.x
	  ntu(6) = main.FP.p.dx				!mr
	  ntu(7) = main.FP.p.y
	  ntu(8) = main.FP.p.dy				!mr
	  ntu(9) = orig.p.delta
	  ntu(10) = orig.p.yptar			!mr
	  ntu(11) = orig.p.xptar			!mr
	  ntu(12) = main.target.z*spec.p.sin_th
	  ntu(13) = recon.e.delta
	  ntu(14) = recon.e.yptar			!mr
	  ntu(15) = recon.e.xptar			!mr
	  ntu(16) = recon.e.z
 	  ntu(17) = main.FP.e.x
	  ntu(18) = main.FP.e.dx			!mr
 	  ntu(19) = main.FP.e.y
 	  ntu(20) = main.FP.e.dy			!mr
	  ntu(21) = orig.e.delta
	  ntu(22) = orig.e.yptar			!mr
	  ntu(23) = orig.e.xptar			!mr
	  ntu(24) = -main.target.z*spec.e.sin_th
	endif
	ntu(25) = recon.q/1000.				!q - GeV/c	
	ntu(26) = recon.omega/1000.			!omega - GeV
	ntu(27) = recon.Q2/1.d6				!Q^2 - (GeV/c)^2
	ntu(28) = recon.W/1000.				!W - GeV/c
	ntu(29) = recon.epsilon				!epsilon
	ntu(30) = recon.Em/1000.			!GeV
	ntu(31) = recon.Pm/1000.			!GeV/c
	ntu(32) = recon.theta_pq			!theta_pq - radians
	ntu(33) = recon.phi_pq				!phi_pq - radians

	if (doing_pion .or. doing_kaon) then
	  ntu(34) = mm/1000.				!missmass (nucleon)
	  ntu(35) = mmA/1000.				!missmass (nucleus)
	  ntu(36) = recon.p.P/1000.			!ppi - GeV/c
	  ntu(37) = t/1.d6				!t - GeV^2
	  ntu(38) = -main.target.rastery		!fry - cm
	  ntu(39) = decdist(4)/1000.			!radphot - GeV
	  dummy = pferx*vertex.uq.x + pfery*vertex.uq.y + pferz*vertex.uq.z
	  if (dummy.eq.0) dummy=1.d-20
	  ntu(40) = pfer/1000.*abs(dummy)/dummy		!p_fermi - GeV/c
	  ntu(41) = main.sigcc				!d5sig
	  ntu(42) = decdist(22)				!pion sig_cm
	  ntu(43) = main.weight
	  ntu(44) = main.davejac			!nuceon at rest --> nucleus at rest jacobian for sig.
	  ntu(45) = decdist(30)				!decay distance (cm)
	  ntu(46) = sqrt(Mh2_final)
	  ntu(47) = pfer/1000.*dummy			!p_fermi along q.

	else if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  ntu(34) = corrsing/1000.
	  ntu(35) = Pm_Heepx/1000.
	  ntu(36) = Pm_Heepy/1000.
	  ntu(37) = Pm_Heepz/1000.
	  ntu(38) = recon.PmPar/1000.
	  ntu(39) = recon.PmPer/1000.
	  ntu(40) = recon.PmOop/1000.
	  ntu(41) = -main.target.rastery		!fry - cm
	  ntu(42) = decdist(4)/1000.			!radphot - GeV
	  ntu(43) = main.sigcc
	  ntu(44) = main.weight
	endif

	call HFN(NtupleID,ntu)
	if (debug(2)) write(6,*)'r_ntu_write: ending'
	return
	end


	subroutine results_ntu_write1(vertex,recon,main,success)

	implicit none
	include 'hbook.inc'
	include 'simulate.inc'

	integer*4 nentries
	parameter (nentries = 9)

	real*8	ntu(nentries)
	logical success
	record /event_main/ main
	record /event/	vertex, recon

	if (debug(2)) write(6,*)'r_ntu_write: entering'
	ntu(1) = vertex.e.delta
	ntu(2) = vertex.e.yptar
	ntu(3) = -vertex.e.xptar
	ntu(4) = main.SP.e.z
	if(success)then
	  ntu(5) = recon.e.delta
	  ntu(6) = recon.e.yptar
	  ntu(7) = recon.e.xptar
	  ntu(8) = recon.e.z
	else
	  ntu(5) = 30.
	  ntu(6) = 0.1
	  ntu(7) = 0.1
	  ntu(8) = 4.
	endif
	ntu(9) = main.weight
	call HFN(NtupleID,ntu)
	if (debug(2)) write(6,*)'r_ntu_write: ending'
	return
	end
