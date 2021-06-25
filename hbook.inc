C hbook.inc
C Sets up common block for hbook.

*-sizes of CERNLIB working space

c	integer HbookSize,HigzSize,KuipSize,PawSize
c	parameter (HbookSize = 100000)
c	parameter (HigzSize =  50000)
c	parameter (KuipSize =  75000)
c	parameter (PawSize = HigzSize+KuipSize+HbookSize+100000)

*-CERNLIB working space

c	integer CernMemory(PawSize)
c	common /PAWC/ CernMemory      !special nonstandard name!

C Ntuple ID stuff

	integer*4	defaultID
	parameter     	(defaultID = 666)
      	character*132 	NtupleDirectory
	character*16    NtupleName
	integer*4       NtupleID,NtupleIO,NtupleSize
  
	common		NtupleID, NtupleDirectory, NtupleName,
     >			NtupleIO, NtupleSize
