	subroutine NtupleClose(filename)
	implicit none
	save

	include 'hbook.inc'

	character*80 filename
	integer*4 cycle

c	cycle= 0				!dummy for HROUT
c	call HCDIR(NtupleDirectory,' ')
c	call HROUT(NtupleID,cycle,' ')		!flush CERNLIB buffers
c	call HREND(NtupleName)			!CERNLIB close file
c	write(6,*)'Closing file:',filename(1:60)
	CLOSE(NtupleIO)				!close IO channel

	return
	end
