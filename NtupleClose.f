	subroutine NtupleClose(filename)
	implicit none
	save

	include 'hbook.inc'

	character*80 filename
	integer*4 cycle

	cycle= 0				!dummy for HROUT
	call HCDIR(NtupleDirectory,' ')
	call HROUT(NtupleID,cycle,' ')		!flush CERNLIB buffers
	call HREND(NtupleName)			!CERNLIB close file
	write(6,*)'Closing file:',filename(1:60)
	CLOSE(NtupleIO)				!close IO channel

	return
	end
