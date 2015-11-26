	subroutine enerloss_new(len,dens,zeff,aeff,epart,mpart,typeflag,Eloss)

	implicit none

	real*8 thick,len,dens,zeff,aeff,epart,mpart,Eloss
	real*8 x,chsi,lambda,gauss1,Eloss_mp,gamma
	real*8 denscorr,CO,hnup,log10bg,I,beta,Eloss_mp_new

	integer typeflag          !1=normal eloss (picked from distribution)
                                  !2=min eloss
	                          !3=max eloss
                                  !4=most probable eloss
        integer numerr
        data numerr /0/

	real*8 me
        parameter(me=0.51099906)

	thick = len*dens
	gamma=epart/mpart
        beta = sqrt(1.-1./gamma**2)

	if(zeff.eq.1) then	!Ionization potential in MeV
	  I = 21.8e-06
	else
	  I = (16.*zeff**0.9)*1.0e-06
	endif

	hnup = 28.816e-06*sqrt(dens*zeff/aeff) !plasma frequency
	log10bg = log(beta*gamma)/log(10.)
	CO=log(hnup)-log(I)+0.5

C DJG Get density effect correction (I got this from JV).

	if(log10bg.lt.0.) then
	  denscorr=0.
	elseif(log10bg.lt.3.) then
	  denscorr=CO+log(10.)*log10bg+abs(CO/27.)*(3.-log10bg)**3
	elseif(log10bg.lt.4.7) then
	  denscorr=CO+log(10.)*log10bg
	else
	  denscorr=CO+log(10.)*4.7
	endif
	   
	if (thick.le.0.) then
	  Eloss = 0.
	else
	  Eloss_mp = 0.1536e-03 * zeff/aeff * thick * ( 19.26 +
     &          log(thick/dens) )
	  Eloss_mp_new = 0.1536e-03 * zeff/aeff *thick/beta**2* (
     &          log(me/I**2) + 1.063 + 2.*log(gamma*beta) + 
     &		log(0.1536*zeff/aeff*thick/beta**2)-beta**2-denscorr)
c	  write(6,*) 'ELOSS',Eloss_mp,Eloss_mp_new 
! ........ convert to MeV, the unit of choice in THIS program
! ........ (cf. EVCOIN where GeV prevail)
	  Eloss_mp = Eloss_mp_new*1000.
	  chsi = 0.307075/2.*zeff/aeff*thick/beta**2
	  if(typeflag.eq.1)then
	    x=abs(gauss1(10.0e0))
	  elseif(typeflag.eq.2)then
	    x=3
	  elseif(typeflag.eq.3)then
	    x=0.0067
	  elseif(typeflag.eq.4)then
	    x=1
          endif
	  if(x.gt.0.0) then
	    lambda = -2.0*log(x)
	  else
	    lambda = 100000.
	  endif
	  Eloss = lambda*chsi+eloss_mp
	endif

        if (eloss.gt.(epart-mpart)) then
	   eloss=(epart-mpart)-0.0000001
	   numerr=numerr+1
	   if (numerr.le.10) then
	      write(6,*) 'Eloss>Total KE; forcing Eloss=KE'
	      if (numerr.eq.10) write(6,*) '     FURTHER ELOSS ERRORS SUPPRESSED'
	   endif
        endif 

	return
	end
