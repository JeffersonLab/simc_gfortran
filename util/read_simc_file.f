      program read_simc_file


      character*16 NtupleTag(80)
      character filename*80
      integer io,i,j,check
      integer*4 ntupsize
      integer nev
      real*8 ntup(80)

      parameter(nev=10000000)

      io=10

      filename='eep_hydrogen_q8.bin'
      open(io,file=filename,form="unformatted",access="sequential")

      read(io) ntupsize

      do i=1,ntupsize
         read(io) NtupleTag(i)
         write(6,*) NtupleTag(i)
      enddo

c now loop over events     
      do j=1,nev
         write(6,*) 'Event number: ', j
         write(6,*) 'cheesy poofs', ntup_size
         do i=1,ntupsize
            read(io,iostat=check) ntup(i)
            if (check.lt.0) then
               write(6,*) 'end of file'
               stop
            endif
            write(6,*) ntup(i)
         enddo
      enddo



      end


