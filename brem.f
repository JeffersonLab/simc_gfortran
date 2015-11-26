! BREM.FOR

! Exact soft photon bremstrahlung calculation of
! delta, delta^\prime, and exponentiations thereof

	real*8 function brem(ein,eout,egamma,radiate_proton,bsoft,bhard,dbsoft)

	implicit none
	include 'brem.inc'

	real*8 pi, am, ame, e2
	parameter (pi=3.141592653589793)
        parameter (am= .93827231)
	parameter (ame= .00051099906)
	parameter (e2= 1./137.0359895)

	real*8 ein		! electron energy
	real*8 eout		! final electron energy
	real*8 egamma		! energy of bremsstrahling photon
	real*8 eang,ak,akp,ap,pang,eta,ape
	real*8 q2,de
	real*8 aprod,adot,ar1,ar2,alpha
	real*8 bpi,bpf,bpp,bei,bef,bee,bepii,bepif,bepfi,bepff
	real*8 dbpi,dbpf,dbpp,dbei,dbef,dbee,dbepii,dbepif,dbepfi,
     >		dbepff !derivatives of b's
	real*8 b,bz,bzz,db,dbz,dbzz,bhard,bsch,spence,bsoft,dbsoft
	real*8 inter
	logical radiate_proton

! Convert to GeV

	ak=ein/1000.
	akp=eout/1000.
	de=egamma/1000.

! electron scattering angle

	eang= 2.*asin((am/(2.*ak)*(ak/akp-1.))**0.5)
	eta= 1.+2.*ak*sin(eang/2.)**2/am
	q2= 4.*ak*akp*sin(eang/2.)**2

! final proton energy/momentum/scattering angle.

	ape= am+ak-akp
	ap= sqrt(ape**2-am**2)
	pang= acos((ak-akp*cos(eang))/ap)

	if(produce_output) then
	  write(6,*)' q2   ',q2
	  write(6,*)' k	   ',ak
	  write(6,*)' kp   ',akp
	  write(6,*)' eang ',eang*180./pi
	  write(6,*)' ap   ',ap
	  write(6,*)' pang ',pang*180./pi
	  write(6,*)' eta  ',eta
	endif

! Calculate components of delta soft

! ... electron terms

! direct initial electron

	aprod= 1.e0
	bei= aprod*(-1./(2.*pi))*log(ak/de)
	dbei= aprod*(1./(2.*pi*de))

! direct final electron

	aprod= 1.e0
	bef= aprod*(-1./(2.*pi))*log(akp/de)
	dbef= aprod*(1./(2.*pi*de))

! e-e interference

	aprod= -1.e0
	adot= ak*akp*(1.-cos(eang))
	alpha= 2.*ame**2-2.*adot
	ar1= 0.5+sqrt(adot**2-ame**4)/alpha
	ar2= 0.5-sqrt(adot**2-ame**4)/alpha
	bee= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,ak,akp,de)
	dbee= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	if(produce_output) write(6,*) ar1,ar2

! ... proton terms

	if (radiate_proton) then

! initial p direct

	  aprod= 1.e0
	  bpi= aprod*(-1./(2.*pi))*log(am/de)
	  dbpi= aprod*(1./(2.*pi*de))

! final p direct

	  aprod= 1.e0
	  bpf= aprod*(-1./(2.*pi))*log(ape/de)
	  dbpf= aprod*(1/(2.*pi*de))

! p-p interference

	  aprod= -1.e0
	  adot= am*ape
	  alpha= 2.*am**2-2.*adot
	  ar1= 0.5+sqrt(adot**2-am**4)/alpha
	  ar2= 0.5-sqrt(adot**2-am**4)/alpha
	  bpp= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,am,ape,de)
	  dbpp= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	  if(produce_output) write(6,*) ar1,ar2

! ei-pi interference

	  aprod= -1.e0
	  adot= ak*am
	  alpha= am**2+ame**2-2.*adot
	  ar1= (am**2-adot+sqrt(adot**2-(ame*am)**2))/alpha
	  ar2= (am**2-adot-sqrt(adot**2-(ame*am)**2))/alpha
	  bepii= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,ak,am,de)
	  dbepii= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	  if(produce_output) write(6,*) ar1,ar2

! ef-pf interference

	  aprod= -1.e0
	  adot= akp*ape-akp*ap*cos(eang+pang)
	  alpha= am**2+ame**2-2.*adot
	  ar1= (am**2-adot+sqrt(adot**2-(ame*am)**2))/alpha
	  ar2= (am**2-adot-sqrt(adot**2-(ame*am)**2))/alpha
	  bepff= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,akp,ape,de)
	  dbepff= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	  if(produce_output) write(6,*) ar1,ar2

!     ei-pf interference

	  aprod= 1.e0
	  adot= ak*ape-ak*ap*cos(pang)
	  alpha= am**2+ame**2-2.*adot
	  ar1= (am**2-adot+sqrt(adot**2-(ame*am)**2))/alpha
	  ar2= (am**2-adot-sqrt(adot**2-(ame*am)**2))/alpha
	  bepif= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,ak,ape,de)
	  dbepif= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	  if(produce_output) write(6,*) ar1,ar2

!     ef-pi interference

	  aprod= 1.e0
	  adot= akp*am
	  alpha= am**2+ame**2-2.*adot
	  ar1= (am**2-adot+sqrt(adot**2-(ame*am)**2))/alpha
	  ar2= (am**2-adot-sqrt(adot**2-(ame*am)**2))/alpha
	  bepfi= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,akp,am,de)
	  dbepfi= aprod*adot/(pi*alpha*(ar1-ar2)*de)*
     >		(log((ar1-1)/ar1)-log((ar2-1)/ar2))
	  if(produce_output) write(6,*) ar1,ar2

	endif ! <radiate_proton>

! All together now!
	b= 2.*e2*(bei+bef+bee)
	if (radiate_proton) then
	  bzz= 2.*e2*(bpi+bpf+bpp)
	  bz= 2.*e2*(bepii+bepff+bepif+bepfi)
	else
	  bzz = 0.0
	  bz = 0.0
	endif
	bsoft=b+bz+bzz
	bhard= -1.*(e2/pi)*(-28/9.+13./6.*log(q2/ame**2))
	db= 2.*e2*(dbei+dbef+dbee)
	if (radiate_proton) then
	  dbzz= 2.*e2*(dbpi+dbpf+dbpp)
	  dbz= 2.*e2*(dbepii+dbepff+dbepif+dbepfi)
	else
	  dbzz = 0.0
	  dbz = 0.0
	endif
	dbsoft= db+dbz+dbzz

	if(produce_output) then
	  write(6,*)'   ----- results ----- '
	  write(6,*)' b=     ',b
	  write(6,*)' bz=    ',bz
	  write(6,*)' bzz=   ',bzz
	  write(6,*)' bhard= ',bhard
	  write(6,*)' total= ',bsoft+bhard
	  write(6,*)' exp=   ',1.-exp(-1.*(bsoft))*(1.-bhard)
	  write(6,*)' '
	  write(6,*)' ultra-relativistic limit'
	  call srad(ak,akp,eang,q2,ame,am,ap,de,produce_output)
	  write(6,*)'  Schwinger Result'
	  bsch= 2.*e2/pi*((log(ak/de)-13./12.)*(log(q2/ame**2)-1.)+17./36.+
     >		0.5*(pi**2/6.-spence((cos(eang/2.))**2)))
	  write(6,*)'  b=     ',bsch
	endif

! ... and the result --> the value of the radiative cross-section,
! dsigma/dEgamma = -dbsoft * exp(-bsoft) * (1-bhard)
! ......... the derivative has dimension 1/[energy] --> convert back to MeV
	dbsoft = dbsoft/1000.
	if (exponentiate) then
	  brem = -dbsoft/exp(bsoft)
	else
	  brem = 1.-dbsoft
	endif
	if (include_hard) brem = brem*(1.-bhard)

	return
	end

!-----------------------------------

	real*8 function inter(calculate_spence,alpha,ar1,ar2,e1,e2,de)

! explicit evaluation of integral ... may or may not ignore spence functions

	implicit none

	real*8 pi
	parameter (pi=3.141592653589793)

	real*8 alpha,ar1,ar2,e1,e2,de
	real*8 de2,amult,arg1,arg2,arg3,arg4,spence
	logical	calculate_spence

	de2 = e1-e2
	amult = -1./(alpha*(ar1-ar2))
	inter = log(abs((e2/de)+ar1*(de2/de)))*log(abs((ar1-1.)/ar1)) -
     >  	log(abs((e2/de)+ar2*(de2/de)))*log(abs((ar2-1.)/ar2))

	if (calculate_spence) then
	  arg1= (de2/(e2+ar1*de2))*(ar1-1.)
	  arg2= (de2/(e2+ar1*de2))*(ar1)
	  arg3= (de2/(e2+ar2*de2))*(ar2-1.)
	  arg4= (de2/(e2+ar2*de2))*(ar2)
	  inter = inter - spence(arg1)+spence(arg2)+spence(arg3)-spence(arg4)
	endif

	inter= inter*amult/(pi)

	return
	end

!-----------------------------------

	subroutine srad(ak,akp,eang,q2,ame,am,ap,de,produce_output)

! calculates ultra-relativistic limit

	implicit none

	real*8 pi, alpha
	parameter (pi = 3.141592653589793)
	parameter (alpha = 1./137.0359895)

	real*8 ak,akp,eang,q2,ame,am,ap
	real*8 eta,dsoft,dhard,de
	real*8 dsoftz,dsoftzz,ep
	logical	produce_output

	ep= sqrt(ap**2+am**2)
	dhard= -13./12.*(log(q2/ame**2)-1.)+17./36.
	eta= 1.+2.*ak*(sin(eang/2.)**2)/am

	dsoft= alpha/pi*log(ak*akp/de**2)*(log(q2/ame**2)-1.)
	dsoftz= alpha/pi*log(eta)*(log(ak*akp/de**2)+log(am*ep/de**2))
	dsoftzz= alpha/pi*log(am*ep/de**2)*(log(q2/am**2)-1.)
	if(dsoftzz.le.0) dsoftzz=0.0

	if (produce_output) then
	  write(6,*)'  '
	  write(6,*)' b=    ',dsoft
	  write(6,*)' bz=   ',dsoftz
	  write(6,*)' bzz=  ',dsoftzz
	endif

	return
	end

!-----------------------------------

	real*8 function spence(ax)

	implicit none

	real*8 pi
	parameter (pi=3.141592653589793)

	real*8 ax,bx

	bx= abs(ax)

! ... N.B. Have replaced the former calculation (commented out) with an
! ... approximate expression -- saves a WHALE of CPU!
!
!	if(bx.lt.1) then
!		spence= ssum(ax,100)
!	else if (ax.gt.1) then
!		spence= -0.5*(log(bx))**2+pi**2/3.-ssum(1./ax,100)
!	else if (ax.le.-1) then
!		spence= -0.5*(log(bx))**2-pi**2/6.-ssum(1./ax,100)
!	else if (ax.eq.1) then
!		spence= pi**2/6.
!	else if (ax.eq.-1) then
!		spence= -pi**2/12.
!	endif
	
	if (bx.le.1) then
	  spence = 0.0
	else
	  spence = -0.5*(log(bx))**2
	endif

	return
	end

!-----------------------------------

	real*8 function ssum(as,n)

	implicit none

	integer n,i
	real*8 as

	ssum= 0.0
	if(as.ne.0) then
	  if(n.ge.10000) write(6,*)'  large n in function ssum (brem.f)'
	  do i= 1,n
	    ssum= ssum+as**i/(1.*i)**2
	  enddo
	endif

	return
	end

!------------------------------------------------------------------------------

	real*8 function bremos(egamma,		! photon energy
     >  	k_ix, k_iy, k_iz,		! incoming electron 3-momentum
     >  	k_fx, k_fy, k_fz,		! scattered electron 3-momentum
     >  	p_ix, p_iy, p_iz,		! incoming proton 3-momentum (pm)
     >  	p_fx, p_fy, p_fz, p_fe,		! scattered proton 4-momentum
     >  	radiate_proton,			! proton radiation flag
     >  	bsoft, bhard, dbsoft)

! Calculation of soft photon radiative correction factor and its derivative
! allowing for both initial and final protons to be offshell.
! conventions identical to on-shell calculation

	implicit none
	include 'brem.inc'

	real*8 pi, twopi, ame, e2, mp
	parameter (pi=3.141592653589793)
	parameter (twopi=2.*pi)
	parameter (ame= .00051099906)
	parameter (e2= 1./137.0359895)
        parameter (mp= .93827231)

	type four_vector
		real*8	e, x, y, z
	end type

	real*8 egamma, de
	real*8 k_ix,k_iy,k_iz
	real*8 k_fx,k_fy,k_fz
	real*8 p_ix,p_iy,p_iz
	real*8 p_fx,p_fy,p_fz,p_fe
	real*8 ami,amf,q2
	real*8 aprod,adot,ar1,ar2,alpha
	real*8 bpi,bpf,bpp,bei,bef,bee,bepii,bepif,bepfi,bepff
	real*8 dbpi,dbpf,dbpp,dbei,dbef,dbee,dbepii,dbepif,
     >		dbepfi,dbepff !derivatives of b's
	real*8 b,bz,bzz,db,dbz,dbzz
	real*8 bsoft,bhard,dbsoft,bsch
	real*8 inter,inter_prime
	logical	radiate_proton
	type(four_vector):: k_i, k_f, p_i, p_f

! Initialize

! ... put input into local variables, while converting energies/momenta to GeV.
	de = egamma/1000.
	k_i%x = k_ix/1000.
	k_i%y = k_iy/1000.
	k_i%z = k_iz/1000.
	k_f%x = k_fx/1000.
	k_f%y = k_fy/1000.
	k_f%z = k_fz/1000.
	p_i%x = p_ix/1000.
	p_i%y = p_iy/1000.
	p_i%z = p_iz/1000.
	p_f%e = p_fe/1000.
	p_f%x = p_fx/1000.
	p_f%y = p_fy/1000.
	p_f%z = p_fz/1000.

! ... compute electron energies
	k_i%e= (k_i%x**2+k_i%y**2+k_i%z**2+ame**2)**0.5
	k_f%e= (k_f%x**2+k_f%y**2+k_f%z**2+ame**2)**0.5
! ........ if he insists on e-m conservation like below, just compute p_i.e
c	p_i%e = k_f%e+p_f%e-k_i%e
	p_i%e = mp

! ... check energy-momentum conservation
!	e_check= abs(k_f%e+p_f%e-k_i%e-p_i%e)
!	x_check= abs(k_f%x+p_f%x-k_i%x-p_i%x)
!	y_check= abs(k_f%y+p_f%y-k_i%y-p_i%y)
!	z_check= abs(k_f%z+p_f%z-k_i%z-p_i%z)
!
!	if((e_check.gt.0.0001).or.(x_check.gt.0.0001).or.
!     +		(y_check.gt.0.0001).or.(z_check.gt.0.0001)) then
!	  write(6,*)' bad kinematics'
!	  return
!	endif

! ... compute Q2 and masses of initial and final protons
	q2=-1.*( (k_f%e-k_i%e)**2-(k_f%x-k_i%x)**2-
     >		(k_f%y-k_i%y)**2-(k_f%z-k_i%z)**2)
	if (produce_output) write(6,*)' q2= ',q2
c	ami= ((p_i%e)**2-(p_i%x)**2-(p_i%y)**2-(p_i%z)**2)**0.5
	ami= mp
	amf= ((p_f%e)**2-(p_f%x)**2-(p_f%y)**2-(p_f%z)**2)**0.5

! Calculate components of delta soft

! ... electron terms

! ........ direct initial electron
	aprod= 1.e0
	bei= aprod*(-1./twopi)*log(k_i%e/de)
	dbei= aprod*(-1./twopi)*(-1./de)

! ........ direct final electron
	aprod= 1.e0
	bef= aprod*(-1./twopi)*log(k_f%e/de)
	dbef= aprod*(-1./twopi)*(-1./de)

! ........ e-e interference
	aprod= -1.e0
	adot= k_i%e*k_f%e-k_i%x*k_f%x-k_i%y*k_f%y-k_i%z*k_f%z
	alpha= 2.*ame**2-2.*adot
	ar1= 0.5+sqrt(4.*adot**2-4.*ame**4)/(2.*alpha)
	ar2= 0.5-sqrt(4.*adot**2-4.*ame**4)/(2.*alpha)
	bee= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,k_i%e,k_f%e,de)
	dbee= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	if (produce_output) write(6,*) ar1,ar2

! ... proton terms

	if (radiate_proton) then

! ........ initial p direct
	  aprod= 1.e0
	  bpi= aprod*(-1./twopi)*log(p_i%e/de)
	  dbpi= aprod*(-1./twopi)*(-1./de)

! ........ final p direct
	  aprod= 1.e0
	  bpf= aprod*(-1./twopi)*log(p_f%e/de)
	  dbpf= aprod*(-1./twopi)*(-1./de)

! ........ p-p interference
	  aprod= -1.e0
	  adot= p_i%e*p_f%e-p_i%x*p_f%x-p_i%y*p_f%y-p_i%z*p_f%z
	  alpha= ami**2+amf**2-2.*adot
	  ar1= (2.*amf**2-2.*adot+sqrt(4.*adot**2-4.*(ami*amf)**2))/(2.*alpha)
	  ar2= (2.*amf**2-2.*adot-sqrt(4.*adot**2-4.*(ami*amf)**2))/(2.*alpha)
	  bpp= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,p_i%e,p_f%e,de)
	  dbpp= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	  if (produce_output) write(6,*) ar1,ar2

! ........ ei-pi interference
	  aprod= -1.e0
	  adot= k_i%e*p_i%e-k_i%x*p_i%x-k_i%y*p_i%y-k_i%z*p_i%z
	  alpha= ami**2+ame**2-2.*adot
	  ar1= (2.*ami**2-2.*adot+sqrt(4.*adot**2-4.*(ame*ami)**2))/(2.*alpha)
	  ar2= (2.*ami**2-2.*adot-sqrt(4.*adot**2-4.*(ame*ami)**2))/(2.*alpha)
	  bepii= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,k_i%e,p_i%e,de)
	  dbepii= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	  if (produce_output) write(6,*) ar1,ar2

! ........ ef-pf interference
	  aprod= -1.e0
	  adot= k_f%e*p_f%e-k_f%x*p_f%x-k_f%y*p_f%y-k_f%z*p_f%z
	  alpha= amf**2+ame**2-2.*adot
	  ar1= (2.*amf**2-2.*adot+sqrt(4.*adot**2-4.*(ame*amf)**2))/(2.*alpha)
	  ar2= (2.*amf**2-2.*adot-sqrt(4.*adot**2-4.*(ame*amf)**2))/(2.*alpha)
	  bepff= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,k_f%e,p_f%e,de)
	  dbepff= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	  if (produce_output) write(6,*) ar1,ar2

! ........ ei-pf interference
	  aprod= 1.e0
	  adot= k_i%e*p_f%e-k_i%x*p_f%x-k_i%y*p_f%y-k_i%z*p_f%z
	  alpha= amf**2+ame**2-2.*adot
	  ar1=(2.*amf**2-2.*adot+sqrt(4.*adot**2-4.*(ame*amf)**2))/(2.*alpha)
	  ar2=(2.*amf**2-2.*adot-sqrt(4.*adot**2-4.*(ame*amf)**2))/(2.*alpha)
	  bepif= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,k_i%e,p_f%e,de)
	  dbepif= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	  if (produce_output) write(6,*) ar1,ar2

! ........ ef-pi interference
	  aprod= 1.e0
	  adot= k_f%e*p_i%e-k_f%x*p_i%x-k_f%y*p_i%y-k_f%z*p_i%z
	  alpha= ami**2+ame**2-2.*adot
	  ar1=(2.*ami**2-2.*adot+sqrt(4.*adot**2-4.*(ame*ami)**2))/(2.*alpha)
	  ar2=(2.*ami**2-2.*adot-sqrt(4.*adot**2-4.*(ame*ami)**2))/(2.*alpha)
	  bepfi= aprod*adot*inter(calculate_spence,alpha,ar1,ar2,k_f%e,p_i%e,de)
	  dbepfi= aprod*adot*inter_prime(alpha,ar1,ar2,de)
	  if (produce_output) write(6,*) ar1,ar2

	endif ! <radiate_proton>

! All together now!
	b= 2.*e2*(bei+bef+bee)
	if (radiate_proton) then
	  bzz= 2.*e2*(bpi+bpf+bpp)
	  bz= 2.*e2*(bepii+bepff+bepif+bepfi)
	else
	  bzz = 0.0
	  bz = 0.0
	endif
	bsoft = b + bz + bzz
	bhard= -1.*(e2/pi)*(-28/9.+13./6.*log(q2/ame**2))
*	bhard= (e2/pi)*(2-1.5*log(q2/ame**2))
*	bhard= bhard+vac(q2)
	db= 2.*e2*(dbei+dbef+dbee)
	if (radiate_proton) then
	  dbzz= 2.*e2*(dbpi+dbpf+dbpp)
	  dbz= 2.*e2*(dbepii+dbepff+dbepif+dbepfi)
	else
	  dbzz = 0.0
	  dbz = 0.0
	endif
	dbsoft= db+dbz+dbzz

	if (produce_output) then
	  write(6,*)'   ----- results ----- '
	  write(6,*)' b+bhard= ',b+bhard
	  write(6,*)' bz=      ',bz
	  write(6,*)' bzz=     ',bzz
	  write(6,*)' bhard=   ',bhard
	  write(6,*)' total=   ',1-bsoft-bhard
	  write(6,*)' exp=     ',exp(-1.*(bsoft))*(1.-bhard)
	  write(6,*)' 1-bhard= ',1-bhard
	  write(6,*)' exps=    ',exp(-1.*(bsoft))
	  write(6,*)' expse=   ',exp(-1.*b)
	  write(6,*)' expsp=   ',exp(-1.*(bz+bzz))
	  write(6,*)' '
	  write(6,*)'  Schwinger Result'
	  bsch= 2.*e2/pi*((log(k_i%e/de)-13./12.)*(log(q2/ame**2)-1.)+17./36.)
	  write(6,*)'  b=     ',bsch
	endif

! ... and the result --> the value of the radiative cross-section,
! dsigma/dEgamma = -dbsoft * exp(-bsoft) * (1-bhard)
! ......... the derivative has dimension 1/[energy] --> convert back to MeV
	dbsoft = dbsoft/1000.	!convert back to MeV
	if (exponentiate) then
	  bremos = -dbsoft/exp(bsoft)
	else
	  bremos = 1.-dbsoft
	endif
	if (include_hard) bremos = bremos*(1.-bhard)

	return
	end

!-----------------------------------
c
c	explicit evaluation of DERIVATIVE of integral (wrt de)
c

	real*8 function inter_prime(alpha,ar1,ar2,de)

	implicit none

	real*8 pi
	parameter (pi=3.141592653589793)

	real*8 alpha,ar1,ar2,de
	real*8 amult

	amult = -1./(alpha*(ar1-ar2))
	inter_prime = (-1./de)*amult/pi*
     >		(log(abs((ar1-1.)/ar1))-log(abs((ar2-1.)/ar2)))

	return
	end

!-----------------------------------

	real*8 function vac(q2)

	implicit none

	real*8 q2
	real*8 am(10),ae(10)
	real*8 dele,dell
	integer*4 j

	real*8 del

	am(1)= 0.00051099906
	ae(1) = 1.0
	am(2)= .1057
	ae(2)= 1.0
	am(3)= 1.782
	ae(3)= 1.0
	am(4)= 0.3
	ae(4)= 2./3.
	am(5)= .3
	ae(5)= 1./3.
	am(6)= .430
	ae(6)= 1./3.
	am(7)= 1.5
	ae(7)= 2./3.
	am(8)= 5.
	ae(8)= 1./3.

	dele= del(q2,am(1),ae(1))
	dell= 0.0
	do 20 j= 1,3
	  dell= dell+del(q2,am(j),ae(j))
 20	continue

	vac= dell

	return
	end

!-----------------------------------

	real*8 function del(q2,bm,be)

	implicit none
	real*8 q2,bm,be
	real*8 x,alpha

	x= 4.*bm**2/q2
	alpha= be**2/137.0359895
	del= -2.*alpha/(3.*3.14159265)*(-5./3.+x+(1-x/2)*(1+x)**0.5*
     >		log(((1+x)**0.5+1)/((1+x)**0.5-1)))

	return
	end
