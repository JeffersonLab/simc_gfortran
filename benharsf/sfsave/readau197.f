	Program readsf

	implicit none

	integer*4 nume,nump
	integer*4 numemax,numpmax
	parameter (numemax=100)
	parameter (numpmax=50)

	integer*4 ie,ip,itmp
	real*8 eval(numemax),pval(numpmax)
	real*8 sf(numemax,numpmax)

	nume=80
	nump=40

	open(unit=11,file='pke197_tot.data',status='old')
	do ip=1,nump
	   read(11,*) pval(ip)
	   do itmp=1,nume/4
	      ie=4*(itmp-1)
	      read(11,*) eval(ie+1),sf(ie+1,ip),eval(ie+2),sf(ie+2,ip),
	1	         eval(ie+3),sf(ie+3,ip),eval(ie+4),sf(ie+4,ip)
	   enddo
	enddo
	close(11)

	open(unit=12,file='benharsf_197.dat',status='new')
	write(12,*) nump,nume
	do ip=1,nump
	   do ie=1,nume
	      write(12,'(f8.2,f8.2,g12.4)') pval(ip),eval(ie),sf(ie,ip)
	   enddo
	enddo
	close(12)

	end
