C hbook.inc
C Sets up common block for hbook.

*-sizes of CERNLIB working space

	integer HbookSize,HigzSize,KuipSize,PawSize
	parameter (HbookSize = 100000)
	parameter (HigzSize =  50000)
	parameter (KuipSize =  75000)
	parameter (PawSize = HigzSize+KuipSize+HbookSize+100000)

*-CERNLIB working space

	integer CernMemory(PawSize)
	common /PAWC/ CernMemory      !special nonstandard name!

C Ntuple ID stuff

	integer*4	defaultID
	parameter     	(defaultID = 666)
      	character*132 	NtupleDirectory
	character*16    NtupleName
	integer*4       NtupleID,NtupleIO,NtupleSize
  
	common		NtupleID, NtupleDirectory, NtupleName,
     >			NtupleIO, NtupleSize
