	Program readsf

	implicit none

	integer*4 nume,nump
	integer*4 numemax,numpmax
	parameter (numemax=200)
	parameter (numpmax=50)

	integer*4 ie,ip,itmp,numtmp
	real*8 eval(numemax),pval(numpmax)
	real*8 de(numemax),dp(numpmax)
	real*8 sf(numemax,numpmax)
	real*8 dum
	character*80 line

	real*8 pi
	parameter (pi=3.14159265)
	real*8 hbarc
	parameter (hbarc=197.327053)

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

! calculate bin spacing - assumes spacing is basically uniform.
	dp(1)=pval(1)
	do ip=2,nump
	   dp(ip)=pval(ip)-pval(ip-1)
	enddo

	do ie=1,nume-1
	   de(ie)=eval(ie+1)-eval(ie)
	enddo
	de(nume)=de(nume-1)

! Convert sf to probablility within de/dp bin.  Mult. by 3/k_F^3 to get norm.
	do ip=1,nump
	   do ie=1,nume
	      sf(ie,ip)=3/(1.33**3)*sf(ie,ip)*de(ie)*dp(ip)*pval(ip)**2
	   enddo
	enddo

	open(unit=12,file='benharsf_nm.dat',status='new')
	write(12,*) nump,nume
	do ip=1,nump
	   do ie=1,nume
	      write(12,'(2f9.3,2g12.4,2f9.3)') hbarc*pval(ip),eval(ie),sf(ie,ip),sf(ie,ip),hbarc*dp(ip),de(ie)
	   enddo
	enddo
	close(12)

	end
