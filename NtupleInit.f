	subroutine NtupleInit(filename)
	implicit none
	save

	include  'hbook.inc'
	include  'simulate.inc'

	character*80 filename,directory
	character*16 NtupleTag(80),name,title
	integer*4 m,io,recl,bank,id,status,i
c	parameter(recl = 1024)
c	parameter(bank = 8000)
	parameter(io = 29)
	parameter(name = 'SimcNtuple')
	parameter(title = 'SIMTUPLE')

c	NtupleID = defaultID
c	id = NtupleID
	NtupleIO = io
c	NtupleName = name

c	call HCDIR(directory,'R') !CERNLIB read current directory
c	call HROPEN(io,name,filename,'N',recl,status)  !CERNLIB
						!directory set to "//TUPLE"
c	if (status.ne.0)then
c	  write(6,*) 'HROPEN error: istat=',status
c	  stop
c	endif

	open(NtupleIO,file=filename,form="unformatted",access="sequential")

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
	if (doing_pion .or. doing_kaon .or. doing_delta) then
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 34
	  m = m+1
	  NtupleTag(m) = 'mmnuc'	! 35
	  m = m+1
	  NtupleTag(m) = 'phad'		! 36
	  m = m+1
	  NtupleTag(m) = 't'		! 37
	  m = m+1
	  NtupleTag(m) = 'pmpar'	! 38
	  m = m+1
	  NtupleTag(m) = 'pmper'	! 39
	  m = m+1
	  NtupleTag(m) = 'pmoop'	! 40
	  m = m+1
	  NtupleTag(m) = 'fry'		! 41		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 42
	  m = m+1
	  NtupleTag(m) = 'pfermi'	! 43
	  m = m+1
	  NtupleTag(m) = 'siglab'	! 44
	  m = m+1
	  NtupleTag(m) = 'sigcm'	! 45
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 46
	  m = m+1
	  NtupleTag(m) = 'decdist'	! 47
	  m = m+1
	  NtupleTag(m) = 'Mhadron'	! 48
	  m = m+1
	  NtupleTag(m) = 'pdotqhat'	! 49
	  m = m+1
	  NtupleTag(m) = 'Q2i'		! 50
	  m = m+1
	  NtupleTag(m) = 'Wi'		! 51
	  m = m+1
	  NtupleTag(m) = 'ti'		! 52
	  m = m+1
	  NtupleTag(m) = 'phipqi'	! 53
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 54
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 55 
	     m = m+1
	     NtupleTag(m) = 'beta' ! 56
	     m = m+1
	     NtupleTag(m) = 'phis' ! 57
	     m = m+1
	     NtupleTag(m) = 'phic' ! 58
	     m = m+1
	     NtupleTag(m) = 'betai' ! 59
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 60
	     m = m+1
	     NtupleTag(m) = 'phici' ! 61
	  endif
	  if (doing_kaon) then
	    m = m+1
	    NtupleTag(m) = 'saghai'	! 54
	    m = m+1
	    NtupleTag(m) = 'factor'	! 55
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

c	call HBOOKN(id,title,NtupleSize,name,bank,NtupleTag) !create Ntuple
c
c	call HCDIR(NtupleDirectory,'R') !record Ntuple directory
c
c	call HCDIR(directory,' ')       !reset CERNLIB directory

c write ntuple size first
	write(NtupleIO) NtupleSize
	do i=1,m
	   write(NtupleIO) NtupleTag(i)
	enddo

	return
	end

