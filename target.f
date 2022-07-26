	subroutine trip_thru_target (narm, zpos, energy, theta, Eloss, radlen,
     &                               mass, typeflag)

	implicit none
	include 'simulate.inc'

	integer narm
	integer typeflag   !1=generate eloss, 2=min, 3=max, 4=most probable
	real*8 zpos, energy, mass, theta
	real*8 Eloss, radlen
	real*8 forward_path, side_path
	real*8 s_target, s_Al, s_kevlar, s_air, s_mylar	! distances travelled
	real*8 Eloss_target, Eloss_Al,Eloss_air		! energy losses
	real*8 Eloss_kevlar,Eloss_mylar			! (temporary)
	real*8 z_can,t,atmp,btmp,ctmp,costmp,th_can	!for the pudding-can target.
	real*8 ecir,ecor,entec,twall,tal,tliquid,tcm
	logical liquid

	s_Al = 0.0
	liquid = targ%Z.lt.2.4

	if (abs(zpos) .gt. (targ%length/2.+1.e-5)) then
	  write(6,*) 'call to trip_thru_target has |zpos| > targ.length/2.'
	  write(6,*) 'could be numerical error, or could be error in target offset'
	  write(6,*) 'zpos=',zpos,'  targ%length/2.=',targ%length/2.
	endif
! Which particle are we interested in?

	goto (10,20,30) narm

! The incoming electron

10	continue
	s_target = (targ%length/2. + zpos) / abs(cos(targ%angle))
	if (liquid) then			!liquid target
	  if (targ%can .eq. 1) then		!beer can (2.8 mil endcap)
	    s_Al = s_Al + 0.0028*inch_cm
	  else if (targ%can .eq. 2) then	!pudding can (5 mil Al, for now)
	    s_Al = s_Al + 0.0050*inch_cm
	 else if (targ%can .eq.3) then          !2017 target 10 cm cells
	    s_Al = s_Al + 0.013                 !avg. of 3 loops
	  endif
	endif

! ... compute distance in radiation lengths and energy loss
	radlen = s_target/targ%X0_cm + s_Al/X0_cm_Al
	call enerloss_new(s_target,targ%rho,targ%Z,targ%A,energy,mass,
     &                  typeflag,Eloss_target)
	call enerloss_new(s_Al,rho_Al,Z_Al,A_Al,energy,mass,typeflag,Eloss_Al)
	Eloss = Eloss_target + Eloss_Al
	return

! The scattered electron

! ... compute distances travelled assuming HMS (SOS)
! ...... 16.0 (8.0) mil of Al-foil for the target chamber exit foil
! ......~15 cm air between scattering chamber and spectrometer vacuum
! ...... 15.0 (5.0) mil of kevlar + 5.0 (3.0) mil mylar for the
! ......	spectrometer entrance foil
! ...... ASSUMES liquid targets are 2.65" wide, and have 5.0 mil Al side walls.

! ... compute distances travelled assuming HRS (E. Schulte thesis)
! ...... 13.0 mil of Al-foil for the target chamber exit foil
! ...... WILD GUESS:  ~15 cm air between scattering chamber and HRS vacuum
! ...... 10.0 mil of Kapton for spectrometer entrance (Use mylar, since 
! .....		X0=28.6cm for Kapton, X0=28.7cm for Mylar) 
! ...... ASSUMES liquid targets are 2.65" wide, and have 5.0 mil Al side walls.

! ... For SHMS, use SOS windows for now (just a space filler for now).

20	continue
	if (electron_arm.eq.1) then		!electron is in HMS
	  s_Al = 0.016*inch_cm
	  s_air = 15
 	  s_kevlar = 0.015*inch_cm
	  s_mylar = 0.005*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta+targ%angle))
	else if (electron_arm.eq.2) then			!SOS
	  s_Al = 0.008*inch_cm
	  s_air = 15
 	  s_kevlar = 0.005*inch_cm
	  s_mylar = 0.003*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (electron_arm.eq.3) then			!HRS-R
	  s_Al = 0.013*inch_cm
	  s_air = 15
 	  s_kevlar = 0.*inch_cm
	  s_mylar = 0.010*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (electron_arm.eq.4) then			!HRS-L
	  s_Al = 0.013*inch_cm
	  s_air = 15
 	  s_kevlar = 0.*inch_cm
	  s_mylar = 0.010*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (electron_arm.eq.5 .or. electron_arm.eq.6) then	!SHMS
	  s_Al = (0.02+0.01)*inch_cm !20mil scattering chamber window + 10 mil entrance window
	  s_air = 57.27
 	  s_kevlar = 0.0
	  s_mylar = 0.0
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	endif
	s_target = forward_path

	if (liquid) then
	  if (targ%can .eq. 1) then		!beer can
	    side_path = 1.325*inch_cm / abs(sin(theta))
	    if (forward_path.lt.side_path) then
	      s_Al = s_Al + 0.005*inch_cm / abs(cos(theta))
	    else
	      s_target = side_path
	      s_Al = s_Al + 0.005*inch_cm / abs(sin(theta))
	    endif
	  else if (targ%can .eq. 2) then	!pudding can (5 mil Al, for now)

! this is ugly.  Solve for z position where particle intersects can.  The
! pathlength is then (z_intersect - z_scatter)/cos(theta)
! Angle from center to z_intersect is acos(z_intersect/R).  Therefore the
! angle between the particle and wall is pi/2 - (theta - theta_intersect)

	    t=tan(theta)**2
	    atmp=1+t
	    btmp=-2*zpos*t
	    ctmp=zpos**2*t-(targ%length/2.)**2
	    z_can=(-btmp+sqrt(btmp**2-4.*atmp*ctmp))/2./atmp
	    side_path = (z_can - zpos)/abs(cos(theta))
            s_target = side_path
	    costmp=z_can/(targ%length/2.)
	    if (abs(costmp).le.1) then
	      th_can=acos(z_can/(targ%length/2.))
	    else if (abs(costmp-1.).le.0.000001) then
	      th_can=0.   !extreme_trip_thru_target can give z/R SLIGHTLY>1.0
	    else
c	      stop 'z_can > can radius in target.f !!!'
	       write(6,*) 'z_can > can radius in target.f !!!',z_can
c	       stop
	    endif
	    s_Al = s_Al + 0.0050*inch_cm/abs(sin(target_pi/2 - (theta - th_can)))

	  else if (targ%can .eq. 3) then	!2017 10 cm cryo cells
C This cell is a cylinder with a round end cap. Wall thickness = 5 mil
C Some dimensions are hardwired for now (or forever...).
	     ecir=1.315*2.54	! endcap inner radius (cm)
	     ecor=(1.315+0.0071)*2.54	! endcap outer radius (cm)
c      targlen=3.942*2.54        ! target length (cm)
	     entec=targ%length-ecir ! entrance to end cap (cm)
	     twall = ecor-ecir

	     tcm=zpos+targ%length/2.0	!length of target already traversed in cm

	     if((tcm+ecir/tan(theta)).lt.entec) then ! e goes through sidewall
		tliquid=ecir/sin(theta) ! liquid target
		tal=twall/sin(theta) ! wall material
	     else
		tliquid=	! e goes through end cap
     >   	     (sqrt(ecir**2-((targ%length-ecir-tcm)*sin(theta))**2)
     >              +(targ%length-ecir-tcm)*cos(theta)) ! liquid target

		tal=		! wall
     >            +(sqrt(ecor**2-((targ%length-ecir-tcm)*sin(theta))**2)
     >            -sqrt(ecir**2-((targ%length-ecir-tcm)*sin(theta))**2))
     >             *twall/(ecor-ecir)                   ! & end cap
	     endif
	     s_Al = s_Al + tal
	     s_target = tliquid
	  endif
	endif		

! ... compute distance in radiation lengths and energy loss
	radlen = s_target/targ%X0_cm + s_Al/X0_cm_Al + s_air/X0_cm_air +
     >		s_kevlar/X0_cm_kevlar + s_mylar/X0_cm_mylar
	call enerloss_new(s_target,targ%rho,targ%Z,targ%A,energy,mass,
     &                    typeflag,Eloss_target)
	call enerloss_new(s_Al,rho_Al,Z_Al,A_Al,energy,mass,typeflag,Eloss_Al)
	call enerloss_new(s_air,rho_air,Z_air,A_air,energy,mass,typeflag,
     &                    Eloss_air)
	call enerloss_new(s_kevlar,rho_kevlar,Z_kevlar,A_kevlar,energy,
     &                    mass,typeflag,Eloss_kevlar)
	call enerloss_new(s_mylar,rho_mylar,Z_mylar,A_mylar,energy,mass,
     &                    typeflag,Eloss_mylar)
	Eloss = Eloss_target + Eloss_Al + Eloss_air + Eloss_kevlar + Eloss_mylar

	return

! The scattered proton

30	continue
	if (hadron_arm.eq.1) then		!proton in HMS
	  s_Al = 0.016*inch_cm
	  s_air = 15
 	  s_kevlar = 0.015*inch_cm
	  s_mylar = 0.005*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta+targ%angle))
	else if (hadron_arm.eq.2) then				!SOS
	  s_Al = 0.008*inch_cm
	  s_air = 15
 	  s_kevlar = 0.005*inch_cm
	  s_mylar = 0.003*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (hadron_arm.eq.3) then				!HRS-R
	  s_Al = 0.013*inch_cm
	  s_air = 15
 	  s_kevlar = 0.*inch_cm
	  s_mylar = 0.010*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (hadron_arm.eq.4) then				!HRS-L

!	  if (abs(zpos/targ.length).lt.0.00001) write(6,*) 'SOMEONE LEFT THE SAFETY WINDOW IN SIMC!!!!!!!'
!	  s_Al = 0.125*inch_cm

	  s_Al = 0.013*inch_cm
	  s_air = 15
 	  s_kevlar = 0.*inch_cm
	  s_mylar = 0.010*inch_cm
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	else if (hadron_arm.eq.5 .or. hadron_arm.eq.6) then	!SHMS
	  s_Al = (0.02+0.01)*inch_cm !20mil scattering chamber window + 10 mil entrance window
	  s_air = 57.27
 	  s_kevlar = 0.0
	  s_mylar = 0.0
	  forward_path = (targ%length/2.-zpos) / abs(cos(theta-targ%angle))
	endif

	s_target = forward_path
	if (liquid) then
	  if (targ%can .eq. 1) then		!beer can
	    side_path = 1.325*inch_cm/ abs(sin(theta))
	    if (forward_path.lt.side_path) then
	      s_Al = s_Al + 0.005*inch_cm / abs(cos(theta))
	    else
	      s_target = side_path
	      s_Al = s_Al + 0.005*inch_cm / abs(sin(theta))
	    endif
	  else if (targ%can .eq. 2) then	!pudding can (5 mil Al, for now)

! this is ugly.  Solve for z position where particle intersects can.  The
! pathlength is then (z_intersect - z_scatter)/cos(theta)
! Angle from center to z_intersect is acos(z_intersect/R).  Therefore the
! angle between the particle and wall is pi/2 - (theta - theta_intersect)

	    t=tan(theta)**2
	    atmp=1+t
	    btmp=-2*zpos*t
	    ctmp=zpos**2*t-(targ%length/2.)**2
	    z_can=(-btmp+sqrt(btmp**2-4.*atmp*ctmp))/2./atmp
	    side_path = (z_can - zpos)/abs(cos(theta))
            s_target = side_path
	    costmp=z_can/(targ%length/2.)
	    if (abs(costmp).le.1) then
	      th_can=acos(z_can/(targ%length/2.))
	    else if (abs(costmp-1.).le.0.000001) then
	      th_can=0.   !extreme_trip_thru_target can give z/R SLIGHTLY>1.0
	    else
c	      stop 'z_can > can radius in target.f !!!'
	      write(6,*) 'z_can > can radius in target.f !!!',t,theta
c	      stop
	    endif
	    s_Al = s_Al + 0.0050*inch_cm/abs(sin(target_pi/2 - (theta - th_can)))

	 else if (targ%can .eq. 3) then	!2017 10 cm cryo cells
C This cell is a cylinder with a round end cap. Wall thickness = 5 mil
C Some dimensions are hardwired for now (or forever...).
	     ecir=1.315*2.54	! endcap inner radius (cm)
	     ecor=(1.315+0.0071)*2.54	! endcap outer radius (cm)
	     entec=targ%length-ecir ! entrance to end cap (cm)
	     twall = ecor-ecir

	     tcm=zpos+targ%length/2.0	!length of target already traversed in cm

	     if((tcm+ecir/tan(theta)).lt.entec) then ! e goes through sidewall
		tliquid=ecir/sin(theta) ! liquid target
		tal=twall/sin(theta) ! wall material
	     else
		tliquid=	! e goes through end cap
     >   	     (sqrt(ecir**2-((targ%length-ecir-tcm)*sin(theta))**2)
     >              +(targ%length-ecir-tcm)*cos(theta)) ! liquid target

		tal=		! wall
     >            +(sqrt(ecor**2-((targ%length-ecir-tcm)*sin(theta))**2)
     >            -sqrt(ecir**2-((targ%length-ecir-tcm)*sin(theta))**2))
     >             *twall/(ecor-ecir)                   ! & end cap
	     endif
	     s_Al = s_Al + tal
	     s_target = tliquid
	  endif
	endif

! ... compute energy losses

	radlen = s_target/targ%X0_cm + s_Al/X0_cm_Al + s_air/X0_cm_air +
     >		s_kevlar/X0_cm_kevlar + s_mylar/X0_cm_mylar
	call enerloss_new(s_target,targ%rho,targ%Z,targ%A,energy,mass,
     &                    typeflag,Eloss_target)
	call enerloss_new(s_Al,rho_Al,Z_Al,A_Al,energy,mass,typeflag,Eloss_Al)
	call enerloss_new(s_air,rho_air,Z_air,A_air,energy,mass,typeflag,
     &                    Eloss_air)
	call enerloss_new(s_kevlar,rho_kevlar,Z_kevlar,A_kevlar,energy,mass,
     &                    typeflag,Eloss_kevlar)
	call enerloss_new(s_mylar,rho_mylar,Z_mylar,A_mylar,energy,mass,
     &                    typeflag,Eloss_mylar)

	Eloss=Eloss_target+Eloss_Al+Eloss_air+Eloss_kevlar+Eloss_mylar


	return
	end

!-----------------------------------------------------------------

	subroutine extreme_trip_thru_target(ebeam, the, thp, pe, pp, z, m)

	implicit none
	include 'target.inc'
	include 'constants.inc'

	type limits
	    real*8 min,max
        end type


	type(limits):: the, thp, pe, pp, z, betap
	integer	i
	real*8 th_corner, th_corner_min, th_corner_max
	real*8 E1, E2, E3, E4, t1, t2, t3, t4
	real*8 zz, th1, th2, m
	real*8 ebeam, energymin, energymax
	logical	liquid

	real*8 zero
	parameter (zero=0.0e0)	!double precision zero for subroutine calls

!Given limiting values for the electron/proton angles, the z-position in the
!target, and beta for the proton, determine min and max losses in target (and
!corresponding amounts of material traversed) This procedure requires knowledge
!of the shape of the target, of course, that's why I've put it in here!
!Bothering to compute this at all, given that it's going to be used for slop
!limits, just proves what a glowing example of Obsessive Compulsive Disorder at
!work i truly am.

	liquid = targ%Z.lt.2.4

! Incoming electron
	call trip_thru_target(1, z%max, ebeam, zero, targ%Eloss(1)%max,
     &                        targ%teff(1)%max, Me, 3)
	call trip_thru_target(1, z%min, ebeam, zero, targ%Eloss(1)%min,
     &                        targ%teff(1)%min, Me, 2)

! Scattered electron
C Above a few MeV, the energy loss increases as a function of electron
C energy, so for any energies we're dealing the max(min) energy loss should
C correspond to the max(min) spectrometer energy.  Note that this isn't
C the perfect range, but it's easier than reproducing the generated limits here

	energymin=pe%min
	energymax=pe%max
! ... solid target is infinitely wide
	if (.not.liquid) then
	  call trip_thru_target(2, z%min, energymax, the%max,
     &           targ%Eloss(2)%max, targ%teff(2)%max, Me, 3)
	  call trip_thru_target(2, z%max, energymin, the%min,
     &           targ%Eloss(2)%min, targ%teff(2)%min, Me, 2)

! ... liquid
	else if (targ%can .eq. 1) then		!beer can
	  if (z%max.ge.targ%length/2.) then
	    th_corner_max = target_pi/2.
	  else
	    th_corner_max = atan(1.25*inch_cm/(targ%length/2.-z%max))
	  endif
	  th_corner_min = atan(1.25*inch_cm/(targ%length/2.-z%min))

! ... max loss: do we have access to a corner shot?  try a hair on either
! ... side, front or side walls could be thicker (too lazy to check!)
	  if (th_corner_min.le.the%max .and. th_corner_max.ge.the%min) then
	    th_corner = max(th_corner_min,the%min)
	    zz = targ%length/2. - 1.25*inch_cm/tan(th_corner)
	    th1 = th_corner-.0001
	    th2 = th_corner+.0001
	  else
	    zz = z%min
	    th1 = the%min
	    th2 = the%max
	  endif
	  call trip_thru_target(2, zz, energymax, th1, E1, t1, Me, 3)
	  call trip_thru_target(2, zz, energymax, th2, E2, t2, Me, 3)
	  targ%Eloss(2)%max = max(E1,E2)
	  targ%teff(2)%max = max(t1,t2)

! ........ min loss: try all possibilities (in min case, no way
! ........ an intermediate z or th will do the trick).
	  targ%Eloss(2)%min = 1.d10
	  do i = 0, 3
	    call trip_thru_target (2, z%min+int(i/2)*(z%max-z%min), energymin,
     >			the%min+mod(i,2)*(the%max-the%min), E1, t1, Me, 2)
	    if (E1 .lt. targ%Eloss(2)%min) then
	      targ%Eloss(2)%min = E1
	      targ%teff(2)%min = t1
	    endif
	  enddo

	else if (targ%can .eq. 2) then		!pudding can

! ... for pudding can, max loss occurs at lowest scattering angle, where
! ... the outgoing particle leaves the can at z=0.  Therefore, take minimum
! ... scattering angle, project back from z=0, x=radius to get z_init.
! ... Minimum z_init is -radius.
	  zz = -(targ%length/2.)/tan(the%min)
	  zz = max (zz,(-targ%length/2.))
	  call trip_thru_target(2, zz, energymax, the%min, targ%Eloss(2)%max,
     &                          targ%teff(2)%max, Me, 3)

! ........ min loss: try all possibilities (in min case, no way
! ........ an intermediate z or th will do the trick).  Actually, this
! ........ may not be true for the pudding can target, but I'm leaving the
! ........ code alone, for now.
	  targ%Eloss(2)%min = 1.d10
	  do i = 0, 3
	    call trip_thru_target (2, z%min+int(i/2)*(z%max-z%min), energymin,
     >			the%min+mod(i,2)*(the%max-the%min), E1, t1, Me, 2)
	    if (E1 .lt. targ%Eloss(2)%min) then
	      targ%Eloss(2)%min = E1
	      targ%teff(2)%min = t1
	    endif
	  enddo
	endif

! Scattered proton. As you can see I'm sufficiently lazy to make
! the code work out whether high or low beta
! maximizes or minimizes loss ... always lower beta --> higher
! losses, it seems

! ... compute extreme energy values
	energymin = sqrt(pp%min**2 + m**2)
	energymax = sqrt(pp%max**2 + m**2)
! ... compute extreme betap values
	betap%min = pp%min / sqrt(pp%min**2 + m**2)
	betap%max = pp%max / sqrt(pp%max**2 + m**2)
! ... solid target is infinitely wide
	if (.not.liquid) then
	  call trip_thru_target (3, z%min, energymin, thp%max, E1, t1, m, 3)
	  call trip_thru_target (3, z%min, energymax, thp%max, E2, t2, m, 3)
	  targ%Eloss(3)%max = max(E1,E2)
	  targ%teff(3)%max = max(t1,t2)

	  call trip_thru_target (3, z%max, energymin, thp%min, E1, t1, m, 2)
	  call trip_thru_target (3, z%max, energymax, thp%min, E2, t2, m, 2)
	  targ%Eloss(3)%min = min(E1,E2)
	  targ%teff(3)%min = min(t1,t2)

! ... liquid
	else if (targ%can .eq. 1) then	!beer can
	  if (z%max .ge. targ%length/2.) then
	    th_corner_max = target_pi/2.
	  else
	    th_corner_max = atan(1.25*inch_cm/(targ%length/2.-z%max))
	  endif
	  th_corner_min = atan(1.25*inch_cm/(targ%length/2.-z%min))
! ........ max loss: do we have access to a corner shot?
! ........ try a hair on either side, front or side walls could be
! thicker (too lazy to check!)
	  if (th_corner_min.le.thp%max .and. th_corner_max.ge.thp%min) then
	    th_corner = max(th_corner_min,thp%min)
	    zz = targ%length/2. - 1.25*inch_cm/tan(th_corner)
	    th1 = th_corner-.0001
	    th2 = th_corner+.0001
	  else
	    zz = z%min
	    th1 = thp%min
	    th2 = thp%max
	  endif
	  call trip_thru_target (3, zz, energymin, th1, E1, t1, m, 3)
	  call trip_thru_target (3, zz, energymin, th2, E2, t2, m, 3)
	  call trip_thru_target (3, zz, energymax, th1, E3, t3, m, 3)
	  call trip_thru_target (3, zz, energymax, th2, E4, t4, m, 3)
	  targ%Eloss(3)%max = max(E1,E2,E3,E4)
	  targ%teff(3)%max = max(t1,t2,t3,t4)

! ........ min loss: try all possibilities (in min case, no way
! ........ an intermediate z or th will do the trick)
	  targ%Eloss(3)%min = 1.d10
	  do i = 0, 3
	    call trip_thru_target (3, z%min+int(i/2)*(z%max-z%min), energymin,
     >		thp%min+mod(i,2)*(thp%max-thp%min), E1, t1, m, 2)
	    if (E1 .lt. targ%Eloss(3)%min) then
	      targ%Eloss(3)%min = E1
	      targ%teff(3)%min = t1
	      zz = z%min+int(i/2)*(z%max-z%min)
	      th1 = thp%min+mod(i,2)*(thp%max-thp%min)
	    endif
	  enddo
	  call trip_thru_target (3, zz, energymax, th1, E1, t1, m, 2)
	  targ%Eloss(3)%min = min(targ%Eloss(3)%min, E1)
	else if (targ%can .eq. 2) then	!pudding can

! ... for pudding can, max loss occurs at lowest scattering angle, where
! ... the outgoing particle leaves the can at z=0.  Therefore, take minimum
! ... scattering angle, project back from z=0, x=radius to get z_init.
! ... Minimum z_init is -radius.
	  zz = -(targ%length/2.)/tan(the%min)
	  zz = max (zz,(-targ%length/2.))

	  call trip_thru_target (3, zz, energymin, thp%min, E1, t1, m, 3)
	  call trip_thru_target (3, zz, energymax, thp%min, E2, t2, m, 3)
	  targ%Eloss(3)%max = max(E1,E2)
	  targ%teff(3)%max = max(t1,t2)

! ........ min loss: try all possibilities (in min case, no way
! ........ an intermediate z or th will do the trick).  This may be
! ........ wrong for the pudding can targ.  Check it out later.
	  targ%Eloss(3)%min = 1.d10
	  do i = 0, 3
	    call trip_thru_target (3, z%min+int(i/2)*(z%max-z%min), energymin,
     >		thp%min+mod(i,2)*(thp%max-thp%min), E1, t1, m, 2)
	    if (E1 .lt. targ%Eloss(3)%min) then
	      targ%Eloss(3)%min = E1
	      targ%teff(3)%min = t1
	      zz = z%min+int(i/2)*(z%max-z%min)
	      th1 = thp%min+mod(i,2)*(thp%max-thp%min)
	    endif
	  enddo
	  call trip_thru_target (3, zz, energymax, th1, E1, t1, m, 2)
	  targ%Eloss(3)%min = min(targ%Eloss(3)%min, E1)

	endif

*JRA*! Extreme multiple scattering
*JRA*
*JRA*! ........ don't consider multiple scattering of incoming electron
*JRA*	targ.musc_max(1) = 0.

! Extreme multiple scattering.  Use nominal beam energy rather than minimum
!  (should be close enough)

	call extreme_target_musc(ebeam,1.e0,
     >		targ%teff(1)%max,targ%musc_max(1),targ%musc_nsig_max)
	call extreme_target_musc(pe%min,1.e0,
     >		targ%teff(2)%max,targ%musc_max(2),targ%musc_nsig_max)
	call extreme_target_musc(pp%min,betap%min,
     >		targ%teff(3)%max,targ%musc_max(3),targ%musc_nsig_max)

	return
	end

!------------------------------------------------------------------

	subroutine target_musc(p,beta,teff,dangles)

	implicit none

	real*8 Es, epsilon, nsig_max
	parameter (Es = 13.6)		!MeV
	parameter (epsilon = 0.088)
	parameter (nsig_max = 3.5)

	real*8 p, beta, teff, dangles(2), dangle, r
	real*8 theta_sigma
	real*8 gauss1

	if (p.lt.25.) write(6,*)
     >		'Momentum passed to target_musc should be in MeV, but p=',p

! Compute rms value for planar scattering angle distribution, cf. PDB
! Note teff is thickness of material, in radiation lengths.

c	theta_sigma = Es/p/beta * sqrt(teff) * (1+epsilon*log10(teff))
C Better form for beta .ne. 1, from Lynch and Dahl, NIM B58 (1991) p.6-10, Eqn. 6
	theta_sigma = Es/p/beta * sqrt(teff) * (1+epsilon*log10(teff/beta**2))

! Compute scattering angles in perpendicular planes.
! Generate two Gaussian numbers BELOW nsig_max.

	dangles(1) = theta_sigma * gauss1(nsig_max)
	dangles(2) = theta_sigma * gauss1(nsig_max)

	return

! Return info about extreme multiple scattering

	entry extreme_target_musc(p, beta, teff, dangle, r)

c	theta_sigma = Es/p/beta * sqrt(teff) * (1+epsilon*log10(teff))
C Better form for beta .ne. 1, from Lynch and Dahl, NIM B58 (1991) p.6-10, Eqn. 6
	theta_sigma = Es/p/beta * sqrt(teff) * (1+epsilon*log10(teff/beta**2))
	dangle = theta_sigma * nsig_max
	r = nsig_max
	return
	end
