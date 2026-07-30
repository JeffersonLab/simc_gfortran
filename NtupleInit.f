	subroutine NtupleInit(filename)

	USE structureModule
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
	NtupleTag(m) = 'Q2'	        ! 27
	m = m+1
	NtupleTag(m) = 'W'		! 28
	m = m+1
	NtupleTag(m) = 'xbj'		! 29
	m = m+1
	NtupleTag(m) = 'epsilon'	! 30
	m = m+1
	NtupleTag(m) = 'Em'		! 31
	m = m+1
	NtupleTag(m) = 'Pm'	        ! 32
	m = m+1
	NtupleTag(m) = 'qi'	        ! 33
	m = m+1
	NtupleTag(m) = 'nui'	        ! 34
	m = m+1
	NtupleTag(m) = 'Q2i'	        ! 35
	m = m+1
	NtupleTag(m) = 'Wi'	        ! 36
	m = m+1
	NtupleTag(m) = 'xbji'	        ! 37
	m = m+1
	NtupleTag(m) = 'epsiloni'	! 38
	m = m+1
	NtupleTag(m) = 'Emi'	        ! 39
	m = m+1
	NtupleTag(m) = 'Pmi'	        ! 40
	if (doing_pion .or. doing_kaon .or. doing_delta .or. doing_semi .or.doing_rho) then
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 41
	  m = m+1
	  NtupleTag(m) = 'mmnuc'	! 42
	  m = m+1
	  NtupleTag(m) = 'phad'		! 43
	  m = m+1
	  NtupleTag(m) = 't'	        ! 44
	  m = m+1
	  NtupleTag(m) = 'ti'	        ! 45
	  m = m+1
	  NtupleTag(m) = 'z' 	        ! 46
	  m = m+1
	  NtupleTag(m) = 'zi' 	        ! 47
	  m = m+1
	  NtupleTag(m) = 'pt2' 	        ! 48
	  m = m+1
	  NtupleTag(m) = 'pt2i'	        ! 49
	  m = m+1
	  NtupleTag(m) = 'thetapq'      ! 50
	  m = m+1
	  NtupleTag(m) = 'phipq'        ! 51
	  m = m+1
	  NtupleTag(m) = 'thetapqi'     ! 52
	  m = m+1
	  NtupleTag(m) = 'phipqi'	! 53
	  m = m+1
	  NtupleTag(m) = 'pmpar'	! 54
	  m = m+1
	  NtupleTag(m) = 'pmper'	! 55
	  m = m+1
	  NtupleTag(m) = 'pmoop'	! 56
	  m = m+1
	  NtupleTag(m) = 'fry'		! 57		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 58
	  m = m+1
	  NtupleTag(m) = 'pfermi'	! 59
	  m = m+1
	  NtupleTag(m) = 'siglab'       ! 60
	  if(doing_semi) then
	     m=m+1
	     NtupleTag(m)='sighad'      ! 61
	  elseif(doing_kaon) then
	     m=m+1
	     NtupleTag(m)='saghai'      ! 61
	     m=m+1
	     NtupleTag(m)='factor'      ! 61.5
	  else
	     m = m+1
	     NtupleTag(m) = 'sigcm' ! 61
	  endif
	  m = m+1
	  NtupleTag(m) = 'Weight' ! 62
	  m = m+1
	  NtupleTag(m) = 'decdist' ! 63
	  m = m+1
	  NtupleTag(m) = 'Mhadron' ! 64
	  if(doing_rho) then
	     m = m+1
	     NtupleTag(m) = 'Mrho' ! 65
	     m = m+1
	     NtupleTag(m) = 'Thrho' ! 66
	     m = m+1
	     NtupleTag(m) = 'mmnuc' ! 67
	  endif          
	  if(doing_pizero) then
	     m = m+1
	     NtupleTag(m) = 'xcal_gamma1' ! 68 
	     m = m+1
	     NtupleTag(m) = 'ycal_gamma1' ! 69 
	     m = m+1
	     NtupleTag(m) = 'Egamma1' ! 70
	     m = m+1
	     NtupleTag(m) = 'Pgamma1x' ! 71
	     m = m+1
	     NtupleTag(m) = 'Pgamma1y' ! 72
	     m = m+1
	     NtupleTag(m) = 'Pgamma1z' ! 73
	     m = m+1
	     NtupleTag(m) = 'xcal_gamma2' ! 74
	     m = m+1
	     NtupleTag(m) = 'ycal_gamma2' ! 75
	     m = m+1
	     NtupleTag(m) = 'Egamma2' ! 76
	     m = m+1
	     NtupleTag(m) = 'Pgamma2x' ! 77
	     m = m+1
	     NtupleTag(m) = 'Pgamma2y' ! 78
	     m = m+1
	     NtupleTag(m) = 'Pgamma2z' ! 79
	  endif
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 80
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 81 
	     m = m+1
 	     NtupleTag(m) = 'beta' ! 82
	     m = m+1
	     NtupleTag(m) = 'phis' ! 83
	     m = m+1
	     NtupleTag(m) = 'phic' ! 84
	     m = m+1
	     NtupleTag(m) = 'betai' ! 85
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 86
	     m = m+1
	     NtupleTag(m) = 'phici' ! 87
	  endif
	else if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	   m = m+1
	   NtupleTag(m) = 'corrsing' ! 41
	   m = m+1
	   NtupleTag(m) = 'Pmx'	! 42		!for Heepcheck
	   m = m+1
	   NtupleTag(m) = 'Pmy'	! 43		!for Heepcheck
	   m = m+1
	   NtupleTag(m) = 'Pmz'	! 44		!for Heepcheck
	   m = m+1
	   NtupleTag(m) = 'PmPar' ! 45
	   m = m+1
	   NtupleTag(m) = 'PmPer' ! 46
	   m = m+1
	   NtupleTag(m) = 'PmOop' ! 47
	   m = m+1
	   NtupleTag(m) = 'fry'	! 48		!+y is up.
	   m = m+1
	   NtupleTag(m) = 'radphot' ! 49
	   m = m+1
	   NtupleTag(m) = 'sigcc' ! 50
	   m = m+1
	   NtupleTag(m) = 'Weight' ! 51
	   m = m+1
	   NtupleTag(m) = 'theta_e' ! 52
	   m = m+1
	   NtupleTag(m) = 'theta_p' ! 53
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

