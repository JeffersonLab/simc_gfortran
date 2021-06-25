      program make_ntuple
C     Program to convert simc .bin file to ntupl
      implicit none

*-sizes of CERNLIB working space

      integer HbookSize,HigzSize,KuipSize,PawSize
      parameter (HbookSize = 100000)
      parameter (HigzSize =  50000)
      parameter (KuipSize =  75000)
      parameter (PawSize = HigzSize+KuipSize+HbookSize+100000)

*-CERNLIB working space

      integer CernMemory(PawSize)
      common /PAWC/ CernMemory  !special nonstandard name!

C Ntuple ID stuff

      integer*4	defaultID
      parameter     	(defaultID = 666)
      character*132 	NtupleDirectory
      character*16    NtupleName
      integer*4       NtupleID,NtupleIO,NtupleSize
      integer*4 recl,status, cycle, bank

      character*16 NtupleTag(80)
      character*80 rawname,filename,ntfilename,directory
      character*16 title
      integer io,i,j,check
      integer chanout
      integer nev
      real*8 ntup(80)
      real*4 ntup_out(80)

      parameter(nev=10000000)

      io=10
      chanout=2
      NtupleID=defaultID
      recl=4096
      bank=8000

      NtupleName='SimcNtuple'
      title='SIMTUPLE'

      call hlimit(PawSize)
c input filename
      write(6,*) 'Enter filename to convert (without .bin extension)'
      read(5,*) rawname
      i=index(rawname,' ')
      filename='../../worksim/'//rawname(1:i-1)//'.bin'
      write(6,*) 'opening file: ',filename
      open(io,file=filename,form="unformatted",access="sequential")

c output filename
      ntfilename='../../worksim/'//rawname(1:i-1)//'.rzdat'
      call HCDIR(directory,'R')
      call HROPEN(chanout,NtupleName,ntfilename,'N',recl,status)

	if (status.ne.0)then
	  write(6,*) 'HROPEN error: istat=',status
	  stop
	endif


      read(io) NtupleSize

      write(6,*) 'Variables in output file:'
      do i=1,NtupleSize
         read(io) NtupleTag(i)
         write(6,*) NtupleTag(i)
      enddo

      call HBOOKN(NtupleID,title,NtupleSize,NtupleName,bank,NtupleTag) !create Ntuple
      
      call HCDIR(NtupleDirectory,'R') !record Ntuple directory

      call HCDIR(directory,' ') !reset CERNLIB directory


c now loop over events     
      do j=1,nev
         do i=1,NtupleSize
            read(io,iostat=check) ntup(i)
            if (check.lt.0) then
               write(6,*) 'end of file'
               cycle=0
               call HCDIR(NtupleDirectory,' ')
               call HROUT(NtupleID,cycle,' ') !flush CERNLIB buffers
               call HREND(NtupleName) !CERNLIB close file
               write(6,*)'Closing file:',filename(1:60)
               close(chanout)
               stop
            endif
            ntup_out(i)=ntup(i)
         enddo ! loop over ntuple variables
         call HFN(NtupleID,ntup_out)         
      enddo ! loop over events

      end


