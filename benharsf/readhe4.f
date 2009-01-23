	Program readsf

	implicit none

	integer*4 nume,nump
	integer*4 numemax,numpmax
	parameter (numemax=100)
	parameter (numpmax=50)

	integer*4 ie,ip,itmp
	real*8 eval(numemax),pval(numpmax)
	real*8 sf(numemax,numpmax)
	real*8 de(numemax),dp(numemax)

	real*8 hbarc
	parameter (hbarc=197.327053)
	real*8 pi
	parameter (pi=3.14159265)

	open(unit=11,file='spec4he.data',status='old')

	read (11,*) nump,nume

	do ip=1,nump
	   read(11,*) pval(ip)
	enddo
	do ie=1,nume
	   read(11,*) eval(ie),de(ie)
	enddo
	do ie=1,nume
	   do itmp=1,12
	      ip=4*(itmp-1)
	      read(11,*) sf(ie,ip+1),sf(ie,ip+2),sf(ie,ip+3),sf(ie,ip+4)
	   enddo
	   read(11,*) sf(ie,49),sf(ie,50)
	enddo
	close(11)

! calculate bin spacing
	do ip=1,nump-1
	   dp(ip)=pval(ip+1)-pval(ip)
	enddo
	dp(nump)=dp(nump-1)

! Convert sf to probablility within de/dp bin.
	do ip=1,nump
	   do ie=1,nume
	      sf(ie,ip)=4*pi*sf(ie,ip)*de(ie)*dp(ip)*pval(ip)**2
	   enddo
	enddo

! Convert from fm^-1 to MeV/c, GeV to MeV.  sf is probablility, so unitless.
	do ip=1,nump
	   pval(ip)=hbarc*pval(ip)
	   dp(ip)=hbarc*dp(ip)
	enddo
	do ie=1,nume
	   eval(ie)=1000.*eval(ie)
	   de(ie)=1000.*de(ie)
	enddo

	open(unit=12,file='benharsf_4.dat',status='new')
	write(12,*) nump,nume
	do ip=1,nump
	   do ie=1,nume
	      write(12,'(2f9.3,2g12.4,2f9.3)') pval(ip),eval(ie),sf(ie,ip),sf(ie,ip),dp(ip),de(ie)
	   enddo
	enddo
	close(12)

	end
