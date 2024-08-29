	subroutine pizero_decay(vertex,success,gamma1,gamma2)

C+_____________________________________________________________________________
!  pizero_decay:  This routine takes the electroproduced pi0 and
!   produces the 2-photon pair that we actually detect (in a calorimeter
!   most likely).
!	
!   This routine should be called before the single arm mc.
!
!   The outgoing pion will be thrown flat in phi and cos(th) in
!   pi0 rest frame.	
C-_____________________________________________________________________________
	USE structureModule
	implicit none
	include 'simulate.inc'

	type(event):: vertex
C Local declarations.
	real*8 ph,eh
	real*8 rph,rth1,rth
	real*8 beta,gamma
	real*8 pxr1,pyr1,pzr1,pxr2,pyr2,pzr2
	real*8 gamma1(4), gamma2(4)
	real*8 bx,by,bz,er,pr
	real*8 pf1,pf2

	real*8 Mgamma

	real*8 grnd

	logical success

C ============================= Executable Code ================================

	parameter(Mgamma=0.0)
C calculate beta and gamma for pi0 rest frame
	ph = vertex%p%P
	eh = sqrt(ph**2+Mh**2)
	beta = ph/eh
	gamma = 1./sqrt(1.-beta*beta)

C Generate center/mass decay angles and momenta.
	rph = grnd()*2.*pi
 100	rth1 = grnd()*2.-1.
	rth = acos(rth1)


	er = Mh/2.0
	if(er.lt.Mgamma) then
	   success=.false.
	   return
	endif
	pr = sqrt(er**2-Mgamma**2)
	pxr1 = pr*sin(rth)*cos(rph)
	pyr1 = pr*sin(rth)*sin(rph)
	pzr1 = pr*cos(rth)

	pxr2 = -pxr1
	pyr2 = -pyr1
	pzr2 = -pzr1
C Boost to Lab frame, calculate new angles and momentum, finish drift

	bx = -beta * vertex%up%x
	by = -beta * vertex%up%y
	bz = -beta * vertex%up%z
	call loren(gamma,bx,by,bz,er,pxr1,pyr1,pzr1,
     >      gamma1(1),gamma1(2),gamma1(3),gamma1(4),pf1)

	call loren(gamma,bx,by,bz,er,pxr2,pyr2,pzr2,
     >      gamma2(1),gamma2(2),gamma2(3),gamma2(4),pf2)


c   do some checks for for energy and momentum conservation in the lab

	if(abs(eh-gamma1(1)-gamma2(1)).gt.0.001) then
	   write(6,*) 'pi0 decay: bad energy sum'
	endif

	if(abs(ph*vertex%up%x-gamma1(2)-gamma2(2)).gt.0.001) then
	   write(6,*) 'pi0 decay: bad px sum'
	endif

	if(abs(ph*vertex%up%y-gamma1(3)-gamma2(3)).gt.0.001) then
	   write(6,*) 'pi0 decay: bad py sum'
	endif

	if(abs(ph*vertex%up%z-gamma1(4)-gamma2(4)).gt.0.001) then
	   write(6,*) 'pi0 decay: bad pz sum'
	endif
	success = .true.
	

	return
	end



