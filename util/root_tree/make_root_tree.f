      program make_root_tree
C     Program to convert simc .bin file to root tree
      implicit none

      character*80 rawname,filename,treefilename
      character*16 NtupleTag(80),varname
      

      integer i,j,nev,check
      integer*4 io
      integer*4 NtupleSize
      integer GetNumBranches
      external GetNumBranches

      real*8 ntup(80)
      real*4 ntup_out(80)

      parameter(nev=10000000)
      io=99

c input filename
      write(6,*) 'Enter filename to convert (without .bin extesntion)'
      read(5,*) rawname
      i=index(rawname,' ')
      filename='../../worksim/'//rawname(1:i-1)//'.bin'
      write(6,*) 'opening file: ',filename
      open(io,file=filename,form="unformatted",access="sequential")

c output filename
      treefilename='../../worksim/'//rawname(1:i-1)//'.root'
      
      read(io) NtupleSize
      call InitRootNT(treefilename,'RECREATE');

      write(6,*) 'Variables in output file:'
      do i=1,NtupleSize
         read(io) NtupleTag(i)
         call AddNtBranch(ntup_out(i),NtupleTag(i))
      enddo

c now loop over events     
      do j=1,nev
         do i=1,NtupleSize
            read(io,iostat=check) ntup(i)
            if (check.lt.0) then
               call PrintNT()
               call RootNTOutp();
               stop
            endif
            ntup_out(i)=ntup(i)
         enddo ! loop over ntuple variables
         call FillNTBranch('all')
      enddo ! loop over events

      call PrintNT()
      call RootNTOutp();

      end


