        subroutine mc_shms_coll(m2,p,p_spec,decay_flag,dflag,success,
     -    coll_flag,pathlen)
CDJG Subroutine to step a particle "through" the collimator,checking for
CDJG energy loss, multiple scattering.  If particle doesn't make it (too
CDJG much energy loss for example), returns .FALSE.

c GH 2024.04.05
c GH code is based on mc_hms_coll.f used in Fpi-2
c GH SHMS collimator dimensions taken from drawing 67238-00006
 
        implicit none
        include 'apertures_shms.inc'
        include 'struct_shms.inc'
        include '../spectrometers.inc'
        
        real*8 m2,p,p_spec
        logical decay_flag, dflag,success,coll_flag
        real*8 grnd
        real*8 h_entr,h_exit,v_entr,v_exit
        real*8 x_off,y_off,z_off

! horiz. and vert. 1/2 gap of aperture
        parameter (h_entr = 8.50)
        parameter (v_entr = 12.50)
        parameter (h_exit = 8.65)
        parameter (v_exit = 12.85)

        parameter (x_off = 0.000)
        parameter (y_off = +0.000)
        parameter (z_off = +25.189)  ! z_entr from mc_shms.f

        integer nstep
        parameter(nstep=20)

        real*8 coll_thick,coll_dens,zcoll,acoll,coll_radl,coll_radw
        real*8 pathlen
        parameter(coll_thick = 6.35)
	parameter(coll_dens = 17.0)
	parameter(zcoll = 69.45)
	parameter(acoll = 171.56797)
	parameter(coll_radl = 0.42084)  ! HMS radl scaled by 6.35/6.30

        real*8 step_size,h_step,v_step
        logical step_flag
        logical slit_h_flag,slit_v_flag,slit_o_flag,rantemp_flag,epart_flag

        real*8 epart,trans,eloss_tot,eloss,rantemp
        real*8 thick,thick_temp
        integer n
        
        eloss_tot=0.
        eloss=0.
        
        step_size = 6.35/nstep

        coll_flag = .FALSE.
        success   = .FALSE.

        slit_h_flag  = .FALSE.
        slit_v_flag  = .FALSE.
        slit_o_flag  = .FALSE.
        rantemp_flag = .FALSE.
        epart_flag   = .FALSE.

        h_step = h_entr
        v_step = v_entr

        epart = sqrt(p**2+m2)

C Now start stepping
C Check aperture, if particle hits collimator, do multiple scattering,eloss
C for next step and then project.

        do n=1,nstep
           step_flag = .FALSE.
           if (abs(ys-y_off).gt.h_step) then
              coll_flag = .TRUE.
              step_flag = .TRUE.
              shmsSTOP_slit_hor = shmsSTOP_slit_hor + 1
              slit_h_flag = .TRUE.
              thick_temp = step_size
           endif
           if (abs(xs-x_off).gt.v_step) then
              coll_flag = .TRUE.
              step_flag = .TRUE.
              shmsSTOP_slit_vert = shmsSTOP_slit_vert + 1
              slit_v_flag = .TRUE.
              thick_temp=step_size
           endif
           if (abs(xs-x_off).gt.(-v_step/h_step*abs(ys-y_off)+3*v_step/2)) then
              coll_flag = .TRUE.
              step_flag = .TRUE.
              shmsSTOP_slit_oct = shmsSTOP_slit_oct + 1
              slit_o_flag = .TRUE.
              thick_temp = step_size
           endif
           
           thick = thick_temp

C If we hit the collimator, check for hadronic reaction, do eloss, mult. scat. for first step
           if(step_flag) then

c TH - un-comment the following line if want no collimator punch through!
c              goto 500

              if(m2.gt.12000. .and. m2.lt.20000.) then !check if pion - if muon, don't check for hadronic interaction
                 call pion_coll_absorb(p,thick,trans)
                 rantemp = grnd()
                 if(rantemp.gt.trans) then 
                    rantemp_flag = .TRUE.
                    goto 500
                 endif
              endif
              coll_radw = thick/coll_radl
! TH - some adjustment here for multiple scattering test
c              call musc(m2,p,coll_radw,thick,dydzs,dxdzs)
              call musc(m2,p,coll_radw,dydzs,dxdzs)
              epart = sqrt(p**2+m2)
              call enerloss_new(thick,coll_dens,zcoll,acoll,epart,
     $             sqrt(m2),1,eloss)
              epart = epart - eloss
              eloss_tot = eloss_tot + eloss
              if(epart.lt.sqrt(m2)) then 
                 epart_flag = .TRUE.
                 goto 500
              endif
              p = sqrt(epart**2-m2)
              dpps = 100.*(p/p_spec - 1.)
           endif

           call project(xs,ys,step_size,decay_flag,dflag,m2,p,pathlen) !drift and decay

           h_step = h_step + (h_exit-h_entr)/nstep
           v_step = v_step + (v_exit-v_entr)/nstep
c           write(6,*) 'step,eloss',n,eloss,step_flag,pathlen
        enddo

        success = .TRUE.

 500    continue  
c        write(6,*) 'step,eloss',n,eloss_tot,success,coll_flag,step_flag,xs,ys
c        if (.not.success) then
c           write(6,*) rantemp_flag,trans,epart_flag,eloss,slit_h_flag,slit_v_flag,slit_o_flag
c        endif
c        write(6,*) 'leaving collimator...',p,success,coll_flag,step_flag
c        if (abs(xs).gt.0.5 .or. abs(ys).gt.0.5) write(6,*) 'xy ',xs,ys

        return
        end
