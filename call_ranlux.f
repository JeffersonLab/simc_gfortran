c NEW VERSION: sgrand and grand functions 
c modified to use much imporvoed function RANLUX
c with luxury level 3. This avoid repeating
c number sequences. 
c Implemented by Peter Bosted in late 2019.
      subroutine sgrnd(seed)
c Initialize random number generator at highest level (3) using
c 32bit-integer seed. In future, could implement option to restart
c from a break point (see CERNLIB description), but that isn't done here
      implicit none
      integer seed,lux,k1,k2
      lux=3
      k1=0
      k2=0
      call RLUXGO(LUX,seed,K1,K2)
      return
      end
c
      subroutine save_random_state(fname)
      implicit none
      integer*4 i,ivec(25)

      character*(*) fname
c
      call RLUXUT(IVEC)
      
      open(unit=1,file=fname)
      do i=1,25
      write(1,*) ivec(i)
      enddo
      close(1)
c
      return
      end
c
      subroutine start_file_random_state(fname)
      implicit none
      integer*4 i,ivec(25)

      character*(*) fname
      open(unit=1,file=fname)
      do i=1,25
      read(1,*) ivec(i)
      enddo
      close(1)
c
      write(*,*) " Starting random sequence with file : ",fname
      call RLUXIN(IVEC)
      
c
      return
      end
c
c
     
c

      double precision function grnd()
      implicit none
      integer len/1000/,latest,i
c      real*4 rvec(1000)
      real*8 rvec(1000)
      
      if(latest.le.0 .or. latest.ge.1000) then
       CALL RANLUX(RVEC,LEN)
       latest=1
      endif
      grnd = rvec(latest)

      latest = latest + 1

      return
      end

