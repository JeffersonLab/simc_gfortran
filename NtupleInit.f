	subroutine NtupleInit(filename)
	implicit none
	save

	include  'hbook.inc'
	include  'simulate.inc'

	character*80 filename,directory
	character*16 NtupleTag(80),name,title
	integer*4 m,io,recl,bank,id,status
	parameter(recl = 1024)
	parameter(bank = 8000)
	parameter(io = 29)
	parameter(name = 'SimcNtuple')
	parameter(title = 'SIMTUPLE')

	NtupleID = defaultID
	id = NtupleID
	NtupleIO = io
	NtupleName = name

	call HCDIR(directory,'R') !CERNLIB read current directory
	call HROPEN(io,name,filename,'N',recl,status)  !CERNLIB
						!directory set to "//TUPLE"
	if (status.ne.0)then
	  write(6,*) 'HROPEN error: istat=',status
	  stop
	endif

	m = 0
	m = m+1
	NtupleTag(m) = 'hsdelta'	!  1
	m = m+1
	NtupleTag(m) = 'hsyptar'	!  2
	m = m+1
	NtupleTag(m) = 'hsxptar'	!  3
	m = m+1
	NtupleTag(m) = 'hsytar'		!  4
	m = m+1
	NtupleTag(m) = 'hsxfp'		!  5
	m = m+1
	NtupleTag(m) = 'hsxpfp'		!  6
	m = m+1
	NtupleTag(m) = 'hsyfp'		!  7
	m = m+1
	NtupleTag(m) = 'hsypfp'		!  8
	m = m+1
	NtupleTag(m) = 'hsdeltai'	!  9
	m = m+1
	NtupleTag(m) = 'hsyptari'	! 10
	m = m+1
	NtupleTag(m) = 'hsxptari'	! 11
	m = m+1
	NtupleTag(m) = 'hsytari'	! 12
	m = m+1
	NtupleTag(m) = 'ssdelta'	! 13
	m = m+1
	NtupleTag(m) = 'ssyptar'	! 14
	m = m+1
	NtupleTag(m) = 'ssxptar'	! 15
	m = m+1
	NtupleTag(m) = 'ssytar'		! 16
	m = m+1
	NtupleTag(m) = 'ssxfp'		! 17
	m = m+1
	NtupleTag(m) = 'ssxpfp'		! 18
	m = m+1
	NtupleTag(m) = 'ssyfp'		! 19
	m = m+1
	NtupleTag(m) = 'ssypfp'		! 20
	m = m+1
	NtupleTag(m) = 'ssdeltai'	! 21
	m = m+1
	NtupleTag(m) = 'ssyptari'	! 22
	m = m+1
	NtupleTag(m) = 'ssxptari'	! 23
	m = m+1
	NtupleTag(m) = 'ssytari'	! 24
	m = m+1
	NtupleTag(m) = 'q'		! 25
	m = m+1
	NtupleTag(m) = 'nu'		! 26
	m = m+1
	NtupleTag(m) = 'Q2'		! 27
	m = m+1
	NtupleTag(m) = 'W'		! 28
	m = m+1
	NtupleTag(m) = 'epsilon'	! 29
	m = m+1
	NtupleTag(m) = 'Em'		! 30
	m = m+1
	NtupleTag(m) = 'Pm'		! 31
	m = m+1
	NtupleTag(m) = 'thetapq'	! 32
	m = m+1
	NtupleTag(m) = 'phipq'		! 33
	if (doing_pion .or. doing_kaon .or. doing_eep .or. doing_eepx 
     >      .or. doing_delta .or. doing_Xphasespace) then
	  m = m+1
	  NtupleTag(m) = 'thetacm'      ! 34
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 35
	  m = m+1
	  NtupleTag(m) = 'mmnuc'	! 36
	  m = m+1
	  NtupleTag(m) = 'phad'		! 37
	  m = m+1
	  NtupleTag(m) = 't'		! 38
	  m = m+1
	  NtupleTag(m) = 'pmpar'	! 39
	  m = m+1
	  NtupleTag(m) = 'pmper'	! 40
	  m = m+1
	  NtupleTag(m) = 'pmoop'	! 41
	  m = m+1
	  NtupleTag(m) = 'fry'		! 42		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 43
	  m = m+1
	  NtupleTag(m) = 'pfermi'	! 44
	  m = m+1
	  NtupleTag(m) = 'siglab'	! 45
	  m = m+1
	  NtupleTag(m) = 'sigcm'	! 46
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 47
	  m = m+1
	  NtupleTag(m) = 'decdist'	! 48
	  m = m+1
	  NtupleTag(m) = 'Mhadron'	! 49
	  m = m+1
c	  NtupleTag(m) = 'pdotqhat'	! 50
	  NtupleTag(m) = 'MM2'	        ! 50
	  m = m+1
	  NtupleTag(m) = 'Q2i'		! 51
	  m = m+1
	  NtupleTag(m) = 'Wi'		! 52
	  m = m+1
	  NtupleTag(m) = 'ti'		! 53
	  m = m+1
	  NtupleTag(m) = 'thetapqi'     ! 54
	  m = m+1
	  NtupleTag(m) = 'phipqi'	! 55
	  m = m+1
	  NtupleTag(m) = 'thetacmi'     ! 56
	  m = m+1
	  NtupleTag(m) = 'PMsigned'     ! 57
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 58
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 59 
	     m = m+1
	     NtupleTag(m) = 'beta' ! 60
	     m = m+1
	     NtupleTag(m) = 'phis' ! 61
	     m = m+1
	     NtupleTag(m) = 'phic' ! 62
	     m = m+1
	     NtupleTag(m) = 'betai' ! 63
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 64
	     m = m+1
	     NtupleTag(m) = 'phici' ! 65
	  endif
	  if (doing_kaon) then
	    m = m+1
	    NtupleTag(m) = 'saghai'	! 58 or 66
	    m = m+1
	    NtupleTag(m) = 'factor'	! 59 or 67
	  else if (doing _eepx.or.doing_Xphasespace) then
	    m = m+1
	    NtupleTag(m) = 'thetamq'    ! 58 or 66
	    m = m+1
	    NtupleTag(m) = 'phimq'      ! 59 or 67
	    m=m+1
	    NtupleTag(m) = 'minus_u'    ! 60 or 68
	    m=m+1
	    NtupleTag(m) = 'phicmi'     ! 61 or 69
	    m=m+1
	    NtupleTag(m) = 'tprimei'     ! 62 or 70
	    m=m+1
	    NtupleTag(m) = 'wcmi'       ! 63 or 71
	    m=m+1
	    NtupleTag(m) = 'epsiloni'   ! 64 or 72
	    m=m+1
	    NtupleTag(m) = 'ui'         ! 65 or 73
	  endif
	else if (doing_semi.or.doing_rho) then
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 34 <- Wprime for semi-inclusive folks
	  m = m+1
	  NtupleTag(m) = 'ppi'		! 35
	  m = m+1
	  NtupleTag(m) = 't'		! 36
	  m = m+1
	  NtupleTag(m) = 'fry'		! 37		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 38
	  m = m+1
	  NtupleTag(m) = 'siglab'	! 39
	  m = m+1
	  NtupleTag(m) = 'sigcent'	! 40
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 41
	  m = m+1
	  NtupleTag(m) = 'decdist'	! 42
	  m = m+1
	  NtupleTag(m) = 'Mhadron'	! 43
	  m = m+1
	  NtupleTag(m) = 'z' 	        ! 44
	  m = m+1
	  NtupleTag(m) = 'zi' 	        ! 45
	  m = m+1
	  NtupleTag(m) = 'pt2' 	        ! 46
	  m = m+1
	  NtupleTag(m) = 'pt2i' 	! 47
	  m = m+1
	  NtupleTag(m) = 'xbj' 	        ! 48
	  m = m+1
	  NtupleTag(m) = 'xbji' 	! 49
	  m = m+1
	  NtupleTag(m) = 'thqi' 	! 50	  
	  m = m+1
	  NtupleTag(m) = 'sighad' 	! 51	  
	  m = m+1
	  NtupleTag(m) = 'jacobian' 	! 52	  
	  m = m+1
	  NtupleTag(m) = 'centjac' 	! 53
	  m = m+1
	  NtupleTag(m) = 'pfermi'       ! 54
	  m = m+1
	  NtupleTag(m) = 'xfermi'       ! 55
	  m = m+1
	  NtupleTag(m) = 'phipqi'       ! 56
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 57
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 58 
	     m = m+1
	     NtupleTag(m) = 'beta' ! 59
	     m = m+1
	     NtupleTag(m) = 'phis' ! 60
	     m = m+1
	     NtupleTag(m) = 'phic' ! 61
	     m = m+1
	     NtupleTag(m) = 'betai' ! 62
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 63
	     m = m+1
	     NtupleTag(m) = 'phici' ! 64
	  endif
	  if(doing_rho) then
	     m = m+1
	     NtupleTag(m) = 'Mrho' ! 57 or 65
	     m = m+1
	     NtupleTag(m) = 'Thrho' ! 58 or 66
	  endif
	    
	else if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  m = m+1
	  NtupleTag(m) = 'corrsing'	! 34
	  m = m+1
	  NtupleTag(m) = 'Pmx'		! 35		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'Pmy'		! 36		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'Pmz'		! 37		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'PmPar'	! 38
	  m = m+1
	  NtupleTag(m) = 'PmPer'	! 39
	  m = m+1
	  NtupleTag(m) = 'PmOop'	! 40
	  m = m+1
	  NtupleTag(m) = 'fry'		! 41		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 42
	  m = m+1
	  NtupleTag(m) = 'sigcc'	! 43
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 44
	endif

!	else		!used to be the if (doing_phsp) option.
!	 m=m+1
!	 NtupleTag(m)='gd'
!	 m=m+1
!	 NtupleTag(m)='gt'
!	 m=m+1
!	 NtupleTag(m)='gp'
!	 m=m+1
!	 NtupleTag(m)='gy'
!	 m=m+1
!	 NtupleTag(m)='rd'
!	 m=m+1
!	 NtupleTag(m)='rt'
!	 m=m+1
!	 NtupleTag(m)='rp'
!	 m=m+1
!	 NtupleTag(m)='ry'
!	 m=m+1
!	 NtupleTag(m)='w'
!	endif

	NtupleSize = m

	call HBOOKN(id,title,NtupleSize,name,bank,NtupleTag) !create Ntuple

	call HCDIR(NtupleDirectory,'R') !record Ntuple directory

	call HCDIR(directory,' ')       !reset CERNLIB directory

	return
	end

