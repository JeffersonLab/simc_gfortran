C Spectral Function lookup table for reading in Benhar's spectral functions.

	integer*4 nPmmax		!max number of Pm bins in SF
	parameter(nPmmax=100)
	integer*4 nEmmax		!max number of Em bins in SF 
	parameter(nEmmax=200)

	integer*4 numPm,numEm		!number of Pm bins, Em bins
	real*8 Pmval(nPmmax),dPm(nPmmax)	!Pm value/width for each Pm bin
	real*8 Emval(nEmmax),dEm(nEmmax)	!Em value/width for each Em bin
	real*8 sfval(nEmmax,nPmmax)	!spectral function at Em,Pm

	real*8 sfnorm(nPmmax)		!Integral over Em at each Pm.
	real*8 sftotnorm		!Total normalization.

	common /sftable/ Pmval,dPm,Emval,dEm,sfval,sfnorm,sftotnorm,numPm,numEm
