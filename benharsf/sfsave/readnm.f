	Program readsf

	implicit none

	integer*4 nume,nump
	integer*4 numemax,numpmax
	parameter (numemax=200)
	parameter (numpmax=50)

	integer*4 ie,ip,itmp,numtmp
	real*8 eval(numemax),pval(numpmax)
	real*8 sf(numemax,numpmax)
	real*8 dum
	character*80 line

	real*8 hbarc
	parameter (hbarc=197.3)

	nump=24
	nume=184

	open(unit=11,file='ompkeb4.data',status='old')
	read (11,'(a80)') line     !header line

	do ip=1,nump
	   read(11,*) pval(ip),dum,numtmp
	   do itmp=1,numtmp/4
	      ie=4*(itmp-1)
	      read(11,*) eval(ie+1),sf(ie+1,ip),eval(ie+2),sf(ie+2,ip),
	1	         eval(ie+3),sf(ie+3,ip),eval(ie+4),sf(ie+4,ip)
	   enddo
	   if (numtmp.eq.81) then
	      read(11,*) eval(81),sf(81,ip)
	   else if (numtmp.ne.184) then
	      stop 'I can only handle 81 and 184 entries per block'
	   endif
	enddo
	close(11)

	open(unit=12,file='benharsf_nm.dat',status='new')
	write(12,*) nump,nume
	do ip=1,nump
	   do ie=1,nume
	      write(12,'(f8.2,f8.2,g12.4)') hbarc*pval(ip),eval(ie),sf(ie,ip)
	   enddo
	enddo
	close(12)

	end
