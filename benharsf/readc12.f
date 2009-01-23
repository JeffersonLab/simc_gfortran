	Program readsf

	implicit none

	integer*4 nume,nump
	integer*4 numemax,numpmax
	parameter (numemax=100)
	parameter (numpmax=50)

	integer*4 ie,ip,itmp
	real*8 eval(numemax),pval(numpmax)
	real*8 de(numemax),dp(numpmax)
	real*8 sf(numemax,numpmax)

	real*8 pi
	parameter (pi=3.14159265)

	nume=80
	nump=40

	open(unit=11,file='pke12_tot.data',status='old')
	do ip=1,nump
	   read(11,*) pval(ip)
	   do itmp=1,nume/4
	      ie=4*(itmp-1)
	      read(11,*) eval(ie+1),sf(ie+1,ip),eval(ie+2),sf(ie+2,ip),
	1	         eval(ie+3),sf(ie+3,ip),eval(ie+4),sf(ie+4,ip)
	   enddo
	enddo
	close(11)

! calculate bin spacing - assumes spacing is basically uniform.
	do ip=1,nump-1
	   dp(ip)=pval(ip+1)-pval(ip)
	enddo
	dp(nump)=dp(nump-1)

	do ie=1,nume-1
	   de(ie)=eval(ie+1)-eval(ie)
	enddo
	de(nume)=de(nume-1)


! Convert sf to probablility within de/dp bin.
        do ip=1,nump
           do ie=1,nume
              sf(ie,ip)=4*pi*sf(ie,ip)*de(ie)*dp(ip)*pval(ip)**2
           enddo
        enddo

	open(unit=12,file='benharsf_12.dat',status='new')
	write(12,*) nump,nume
	do ip=1,nump
	   do ie=1,nume
	      write(12,'(2f9.3,2g12.4,2f9.3)') pval(ip),eval(ie),sf(ie,ip),sf(ie,ip),dp(ip),de(ie)
	   enddo
	enddo
	close(12)

	end
