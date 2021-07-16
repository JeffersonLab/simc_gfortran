* physics_piKmodel.f
*
* Purpose:
* Determine the kinematics in the PHOTON-NUCLEON center of mass frame,
* calculate some of the kinematical variables (s,t, and CM quantities in
* the 'main' structure), and return the p(e,e'pi+)n and p(e,e'K+)L
* electroproduction cross section which is interpolated from the VGL or
* VR model

*******************************************************************************
      real*8 function peepi(vertex,main)

*   output:
*     peepi          !d5sigma/dEe'dOmegae'Omegapi     (microbarn/MeV/sr^2)

      implicit none
      include 'simulate.inc'

      type(event_main):: main
      type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

      real*8 sigma_eepi          !final cross section (returned as peepi)
      real*8 sig219,sig,fac
      real*8 wfactor,Wset
      real*8 sigt,sigl,siglt,sigtt !components of dsigma/dt
      real*8 epsi               !epsilon of virtual photon
      real*8 gtpr               !gamma_t prime.
      real*8 tcos,tsin          !cos/sin of theta between ppi and q
      real*8 tfcos,tfsin        !cos/sin of theta between pfermi and q
      real*8 ppix,ppiy,ppiz     !pion momentum in lab.
      real*8 zero

      real*8 s,s_gev,t,t_gev,Q2_g !t,s, Q2 in (GeV/c)**2
      real*8 mtar_gev

      real*8 new_x_x,new_x_y,new_x_z
      real*8 new_y_x,new_y_y,new_y_z
      real*8 new_z_x,new_z_y,new_z_z
      real*8 dummy,p_new_x,p_new_y
      real*8 davesig,phipq,cospq,sinpq
      real*8 square_root,dp_dcos_num,dp_dcos_den,dp_dcos
      real*8 dp_dphi_num,dp_dphi_den,dp_dphi
      real*8 dt_dcos_lab,dt_dphi_lab,psign
      real*8 dpxdphi,dpydphi,dpxdcos,dpydcos,dpzdcos,dpzdphi
      real*8 dpxnewdphi,dpynewdphi,dpxnewdcos,dpynewdcos
      real*8 dpznewdphi,dpznewdcos
      real*8 dphicmdphi,dphicmdcos

      real*8 pbeam,beam_newx,beam_newy,beam_newz
      real*8 ppicm_newx,ppicm_newy,ppicm_newz

      real*8 dEcmdcos,dEcmdphi,dcoscmdcos,dcoscmdphi
      real*8 qx,qy,qz,px,py,pz

      real*8 t_high,t_low,delta_t,delta_Q2,Q2_tmp,t_tmp
      real*8 q2_high,q2_low

      complex*8 A1,A2,A3,A4,A12,A34,A1234
      complex*8 XL,XT,XTT,XLT
      real*8 sinterp(4)

      integer Q2_count,t_count,itab,j,pcount

      logical first

c Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz           !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz          !q in C.M.
      real*8 epicm,ppicm,ppicmx,ppicmy,ppicmz           !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz      !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian
      real*8 pm2_tmp,q2_cent_g,W_cent_g

      logical first_call
      save first_call

      data first_call/.true./
      real*8 Q2_tab(18,130),W_tab(18,130),t_tab(18,130),stab(18,130,4)
      character*28 filename(18)
      integer nfiles, tbins(18)

*******************************************************************************
* Read model values when first called.
c Columns are: Q2, W, -t (GeV2), dsig_L, T, TL, TT (mubarn/GeV2)

      if(first_call) then
         first_call=.false.

c calculate central kinematics Q^2 and W
         q2_cent_g = ( (Ebeam - spec%e%P*cos(spec%e%theta))**2
     c                 + (spec%e%P*sin(spec%e%theta))**2
     c                 - (Ebeam - spec%e%P)**2 ) /1.e6
         W_cent_g = sqrt( (Ebeam - spec%e%P + Mp)**2 
     c                 - (Ebeam - spec%e%P*cos(spec%e%theta))**2
     c                 - (spec%e%P*sin(spec%e%theta))**2 ) /1.e3

c         write(6,*)' Q2_central = ',q2_cent_g,W_cent_g

         if (which_pion.eq.0) then ! pi+
            if (q2_cent_g.lt.4.25 .and. W_cent_g.le.2.35) then
               write(6,*)'Selecting VR pi+ model for Q^2=0.1-6.0, W=2.10'
               nfiles=13
               Wset=2.10
               filename(1)='VR/ep_enpip_q20p10_w2p10.dat'
               filename(2)='VR/ep_enpip_q20p50_w2p10.dat'
               filename(3)='VR/ep_enpip_q21p00_w2p10.dat'
               filename(4)='VR/ep_enpip_q21p50_w2p10.dat'
               filename(5)='VR/ep_enpip_q22p00_w2p10.dat'
               filename(6)='VR/ep_enpip_q22p50_w2p10.dat'
               filename(7)='VR/ep_enpip_q23p00_w2p10.dat'
               filename(8)='VR/ep_enpip_q23p50_w2p10.dat'
               filename(9)='VR/ep_enpip_q24p00_w2p10.dat'
               filename(10)='VR/ep_enpip_q24p50_w2p10.dat'
               filename(11)='VR/ep_enpip_q25p00_w2p10.dat'
               filename(12)='VR/ep_enpip_q25p50_w2p10.dat'
               filename(13)='VR/ep_enpip_q26p00_w2p10.dat'
               filename(14)=' '
               filename(15)=' '
               filename(16)=' '
               filename(17)=' '
               filename(18)=' '
            elseif (q2_cent_g.lt.6.25 .and. W_cent_g.gt.2.35) then
               write(6,*)'Selecting VR pi+ model for Q^2=1.0-9.5, W=3.28'
               nfiles=18
               Wset=3.28
               filename(1)='VR/ep_enpip_q21p00_w3p28.dat'
               filename(2)='VR/ep_enpip_q21p50_w3p28.dat'
               filename(3)='VR/ep_enpip_q22p00_w3p28.dat'
               filename(4)='VR/ep_enpip_q22p50_w3p28.dat'
               filename(5)='VR/ep_enpip_q23p00_w3p28.dat'
               filename(6)='VR/ep_enpip_q23p50_w3p28.dat'
               filename(7)='VR/ep_enpip_q24p00_w3p28.dat'
               filename(8)='VR/ep_enpip_q24p50_w3p28.dat'
               filename(9)='VR/ep_enpip_q25p00_w3p28.dat'
               filename(10)='VR/ep_enpip_q25p50_w3p28.dat'
               filename(11)='VR/ep_enpip_q26p00_w3p28.dat'
               filename(12)='VR/ep_enpip_q26p50_w3p28.dat'
               filename(13)='VR/ep_enpip_q27p00_w3p28.dat'
               filename(14)='VR/ep_enpip_q27p50_w3p28.dat'
               filename(15)='VR/ep_enpip_q28p00_w3p28.dat'
               filename(16)='VR/ep_enpip_q28p50_w3p28.dat'
               filename(17)='VR/ep_enpip_q29p00_w3p28.dat'
               filename(18)='VR/ep_enpip_q29p50_w3p28.dat'
            else
               write(6,*)'Selecting VR pi+ model for Q^2=6.5-11.5, W=2.89'
               nfiles=11
               Wset=2.89
               filename(1)='VR/ep_enpip_q26p50_w2p89.dat'
               filename(2)='VR/ep_enpip_q27p00_w2p89.dat'
               filename(3)='VR/ep_enpip_q27p50_w2p89.dat'
               filename(4)='VR/ep_enpip_q28p00_w2p89.dat'
               filename(5)='VR/ep_enpip_q28p50_w2p89.dat'
               filename(6)='VR/ep_enpip_q29p00_w2p89.dat'
               filename(7)='VR/ep_enpip_q29p50_w2p89.dat'
               filename(8)='VR/ep_enpip_q210p0_w2p89.dat'
               filename(9)='VR/ep_enpip_q210p5_w2p89.dat'
               filename(10)='VR/ep_enpip_q211p0_w2p89.dat'
               filename(11)='VR/ep_enpip_q211p5_w2p89.dat'
               filename(12)=' '
               filename(13)=' '
               filename(14)=' '
               filename(15)=' '
               filename(16)=' '
               filename(17)=' '
               filename(18)=' '
            endif

      elseif (which_pion.eq.1) then !pi-
            if (q2_cent_g.lt.6.25) then
               write(6,*)'Selecting VR pi- model for Q^2=1.0-6.5, W=3.28'
               nfiles=12
               Wset=3.28
               filename(1)='VR/en_eppim_q21p00_w3p28.dat'
               filename(2)='VR/en_eppim_q21p50_w3p28.dat'
               filename(3)='VR/en_eppim_q22p00_w3p28.dat'
               filename(4)='VR/en_eppim_q22p50_w3p28.dat'
               filename(5)='VR/en_eppim_q23p00_w3p28.dat'
               filename(6)='VR/en_eppim_q23p50_w3p28.dat'
               filename(7)='VR/en_eppim_q24p00_w3p28.dat'
               filename(8)='VR/en_eppim_q24p50_w3p28.dat'
               filename(9)='VR/en_eppim_q25p00_w3p28.dat'
               filename(10)='VR/en_eppim_q25p50_w3p28.dat'
               filename(11)='VR/en_eppim_q26p00_w3p28.dat'
               filename(12)='VR/en_eppim_q26p50_w3p28.dat'
               filename(13)=' '
               filename(14)=' '
               filename(15)=' '
               filename(16)=' '
               filename(17)=' '
               filename(18)=' '
            else
               write(6,*)'Selecting VR pi- model for Q^2=6.5-11.5, W=2.89'
               nfiles=11
               Wset=2.89
               filename(1)='VR/en_eppim_q26p50_w2p89.dat'
               filename(2)='VR/en_eppim_q27p00_w2p89.dat'
               filename(3)='VR/en_eppim_q27p50_w2p89.dat'
               filename(4)='VR/en_eppim_q28p00_w2p89.dat'
               filename(5)='VR/en_eppim_q28p50_w2p89.dat'
               filename(6)='VR/en_eppim_q29p00_w2p89.dat'
               filename(7)='VR/en_eppim_q29p50_w2p89.dat'
               filename(8)='VR/en_eppim_q210p0_w2p89.dat'
               filename(9)='VR/en_eppim_q210p5_w2p89.dat'
               filename(10)='VR/en_eppim_q211p0_w2p89.dat'
               filename(11)='VR/en_eppim_q211p5_w2p89.dat'
               filename(12)=' '
               filename(13)=' '
               filename(14)=' '
               filename(15)=' '
               filename(16)=' '
               filename(17)=' '
               filename(18)=' '
            endif

      endif

c      nfiles=5
c      Wset=2.89
c      data filename /'VGL/ep_enpip_q27p50_w2p89_pirho.dat',
c     $     'VGL/ep_enpip_q28p00_w2p89_pirho.dat',
c     $     'VGL/ep_enpip_q28p50_w2p89_pirho.dat',
c     $     'VGL/ep_enpip_q29p00_w2p89_pirho.dat',
c     $     'VGL/ep_enpip_q29p50_w2p89_pirho.dat'/

         do j=1,nfiles
            OPEN(3,FILE=filename(j),STATUS='OLD')
            write(6,80) j,filename(j)
 80         format(' Reading model file',i3,'  ',a28)
            do t_count = 1,130
               read(3,*,end=100) Q2_tab(j,t_count),
     1              W_tab(j,t_count),t_tab(j,t_count),
     2              stab(j,t_count,1), stab(j,t_count,2), 
     3              stab(j,t_count,3), stab(j,t_count,4)
               if (abs(W_tab(j,t_count)-Wset).gt.0.01) then
                  write(6,*)' W read error ',W_tab(j,t_count),Wset
                  stop
               endif
            enddo
 100        continue
            tbins(j)=t_count-1
c            write(6,*)' lines found ',tbins(j),t_tab(j,tbins(j))
            close(unit=3)
c            write(6,*) 'done'
         enddo
         pcount=0
      endif                     !if first time.
*******************************************************************************

* Initialize some stuff.
      Q2_g = vertex%Q2/1.e6
c NOTE: phipq calculation in event.f reverted to original version.
      phipq= main%phi_pq
      mtar_gev = targ%Mtar_struck/1000.

      cospq = cos(phipq)
      sinpq = sin(phipq)

* calculate energy of initial (struck) nucleon, using the same assumptions that
* go into calculating the pion angle/momentum (in event.f).  For A>1, the struck
* nucleon is off shell, the 2nd nucleon (always a neutron) is on shell, and has
* p = -p_fermi, and any additional nucleons are at rest.
* NOTE pfer, efer appear to be in MeV at this point.
      efer = sqrt(pfer**2+targ%Mtar_struck**2)

      if(doing_deutpi.or.doing_hepi) then
         efer = targ%M-sqrt(mn**2+pfer**2)
         if(doing_hepi)efer=efer-mp
c         mtar_offshell = sqrt(efer**2-pfer**2)
      endif

* calculate some kinematical variables
* f's and fer indicate fermi momenta, s, star or cm CM system

      tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y
     1     +vertex%up%z*vertex%uq%z
      if(tcos-1..gt.0..and.tcos-1..lt.1.e-8)tcos=1.0
      tsin=sqrt(1.-tcos**2)
        
      tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
      if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0
      tfsin=sqrt(1.-tfcos**2)

      epsi = 1./(1.+2*(1.+vertex%nu**2/vertex%Q2)*tan(vertex%e%theta/2.)**2)

c GH:
c calculate with target nucleon at rest (as in experimental replay)
      s = (vertex%nu+targ%Mtar_struck)**2-vertex%q**2
      main%wcm = sqrt(s)
      s_gev = s/1.e6

c GH:
c calculate t the same way as in experimental replay
      pm2_tmp = (vertex%q*vertex%uq%x-vertex%p%P*vertex%up%x)**2 +
     1     (vertex%q*vertex%uq%y-vertex%p%P*vertex%up%y)**2 +
     2     (vertex%q*vertex%uq%z-vertex%p%P*vertex%up%z)**2
      t = -(mn-mp)**2 +2*targ%Mtar_struck*
     1     ( sqrt(targ%Mrec_struck**2+pm2_tmp)-targ%Mrec_struck )

      t_gev = t/1.e6               !CONVERT TO (GeV/c)**2
      main%t = t

*******************************************************************************
c VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)
* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the pion momentum perpendicular to q.

      dummy=sqrt((vertex%uq%x**2+vertex%uq%y**2) * 
     1     (vertex%uq%x**2+vertex%uq%y**2+vertex%uq%z**2))
      new_x_x = -(-vertex%uq%x)*vertex%uq%z/dummy
      new_x_y = -(-vertex%uq%y)*vertex%uq%z/dummy
      new_x_z = ((-vertex%uq%x)**2 + (-vertex%uq%y)**2)/dummy

      dummy   = sqrt(vertex%uq%x**2 + vertex%uq%y**2)
      new_y_x =  (-vertex%uq%y)/dummy
      new_y_y = -(-vertex%uq%x)/dummy
      new_y_z =  0.0

* get beam in q+nucleon CM system.

      pbeam = sqrt(vertex%Ein**2-me**2)
      beam_newx = pbeam*new_x_z
      beam_newy = pbeam*new_y_z
      beam_newz = pbeam*vertex%uq%z

      bstar=vertex%q/(vertex%nu+targ%Mtar_struck)
      gstar=1./sqrt(1. - bstar**2)

      zero =0.e0
      bstarx = zero
      bstary = zero
      bstarz = vertex%q/(targ%Mtar_struck+vertex%nu)

* DJG: Boost virtual photon to q+nucleon CM.

      call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >     zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost pion to q+nucleon CM.

      ppiz = vertex%p%P*tcos
      ppix = vertex%p%P*tsin*cos(phipq)
      ppiy = vertex%p%P*tsin*sin(phipq)
      call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >     ppix,ppiy,ppiz,epicm,ppicmx,ppicmy,ppicmz,ppicm)
      thetacm = acos((ppicmx*qstarx+ppicmy*qstary+ppicmz*qstarz)/ppicm/qstar)
      main%pcm = ppicm

* DJG Boost the beam to q+nucleon CM.

      call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_newx,
     >     beam_newy,beam_newz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* Thetacm is defined as angle between ppicm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

      dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >     (qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >     +(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
      new_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
      new_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
      new_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

      dummy = sqrt((new_y_y*qstarz-new_y_z*qstary)**2
     >     +(new_y_z*qstarx-new_y_x*qstarz)**2
     >     +(new_y_x*qstary-new_y_y*qstarx)**2)
      new_x_x = (new_y_y*qstarz-new_y_z*qstary)/dummy
      new_x_y = (new_y_z*qstarx-new_y_x*qstarz)/dummy
      new_x_z = (new_y_x*qstary-new_y_y*qstarx)/dummy

      new_z_x = qstarx/qstar
      new_z_y = qstary/qstar
      new_z_z = qstarz/qstar

      ppicm_newx = ppicmx*new_x_x + ppicmy*new_x_y + ppicmz*new_x_z
      ppicm_newy = ppicmx*new_y_x + ppicmy*new_y_y + ppicmz*new_y_z
      ppicm_newz = ppicmx*new_z_x + ppicmy*new_z_y + ppicmz*new_z_z

      phicm = atan2(ppicm_newy,ppicm_newx)
      if(phicm.lt.0.) phicm = 2.*3.141592654+phicm

      main%thetacm = thetacm
      main%phicm = phicm

c      write(6,*)'  '
c      write(6,*)' pfer ',pfer
c      write(6,*)' t ',t
c      write(6,*)' s ',s
c      write(6,*)' thetacm ',thetacm*180./3.14159
c      write(6,*)' phicm ',phicm*180./3.14159,phipq*180./3.14159

*******************************************************************************
* Interpolate model from tables.

      t_tmp = abs(t_gev)
      Q2_tmp = Q2_g
      if(Q2_tmp.lt.Q2_tab(1,1)) then
         if (pcount.le.500) then
            write(6,*)' WARNING: model table range too high Q2 ',Q2_tmp
            pcount=pcount+1
         endif
         Q2_tmp = Q2_tab(1,1)
      endif
      if(Q2_tmp.gt.Q2_tab(nfiles,1)) then
         if (pcount.lt.500) then
            write(6,*)' WARNING: model table range too low Q2 ',Q2_tmp
            pcount=pcount+1
         endif
         Q2_tmp = Q2_tab(nfiles,1)
      endif

      do Q2_count=1,(nfiles-1)
         Q2_high=0.
         Q2_low=0.
         t_high=0.
         t_low=0.
         if((Q2_tmp.ge.Q2_tab(Q2_count,1)).and.
     1        (Q2_tmp.lt.Q2_tab(Q2_count+1,1))) then
            Q2_high = Q2_tab(Q2_count+1,1)
            Q2_low = Q2_tab(Q2_count,1)
            delta_Q2 = (Q2_high-Q2_low)
c            write(6,*)' Q2 ',Q2_tab(Q2_count,1),Q2_tmp

            if(t_tmp.lt.t_tab(Q2_count,1)) then
               if (pcount.le.500) then
                  write(6,*)' WARNING: model table range too high -t ',t_tmp
                  pcount=pcount+1
               endif
               t_tmp = t_tab(Q2_count,1)
            endif
            if(t_tmp.gt.t_tab(Q2_count,tbins(Q2_count)-1)) then
               if (pcount.lt.500) then
                  write(6,*)' WARNING: model table range too low -t ',t_tmp
                  pcount=pcount+1
               endif
               t_tmp = t_tab(Q2_count,tbins(Q2_count)-1)
            endif
            do t_count=1,(tbins(Q2_count)-1)
               if((t_tmp.ge.t_tab(Q2_count,t_count)).and.
     1              (t_tmp.lt.t_tab(Q2_count,t_count+1))) then
                  t_high = t_tab(Q2_count,t_count+1)
                  t_low =  t_tab(Q2_count,t_count)
                  delta_t = t_high-t_low
c                  write(6,*)' t ',t_tab(Q2_count,t_count),t_tmp

               do itab=1,4
                 A1 = stab(Q2_count,t_count,itab)
                 A2 = stab(Q2_count+1,t_count,itab)
                 A3 = stab(Q2_count,t_count+1,itab)
                 A4 = stab(Q2_count+1,t_count+1,itab)
                 A12= (A1*(Q2_high-Q2_tmp)+A2*(Q2_tmp-Q2_low))/delta_Q2
                 A34= (A3*(Q2_high-Q2_tmp)+A4*(Q2_tmp-Q2_low))/delta_Q2
                 A1234 = (A12*(t_high-t_tmp)+
     1                A34*(t_tmp-t_low))/delta_t
                 sinterp(itab) = A1234
               enddo
               endif               !t check
            enddo                  !t
         endif                     !Q2 check
      enddo                        !Q2
     
      sigl = sinterp(1)
      sigt = sinterp(2)
      siglt = sinterp(3)
      sigtt = sinterp(4)

      sig219=(sigt+main%epsilon*sigl+main%epsilon*cos(2.*phicm)*sigtt
     >     +sqrt(2.0*main%epsilon*(1.+main%epsilon))*cos(phicm)*siglt)/1.e0
      
c now convert to different W
c W dependence given by 1/(W^2-M^2)^2
     
      wfactor=(Wset**2-(targ%Mtar_struck/1000.)**2)**2/
     1     (s_gev-(targ%Mtar_struck/1000.)**2)**2
          
      sig=sig219*wfactor
      sigl=sigl*wfactor
      sigt=sigt*wfactor
      sigtt=sigtt*wfactor
      siglt=siglt*wfactor

C--->Debug
c      write(*,*) 's =',s_gev
c      write(*,*) 'wfactor =',wfactor
c      write(*,*) 'sig =',sig
c      write(*,*) 'sigL  =',sigL
c      write(*,*) 'sigT  =',sigT
c      write(*,*) 'sigLT =',sigLT
c      write(*,*) 'sigTT =',sigTT
c      write(*,*) '-----------------------------------------------------'
C--->Debug

      sig=sig/2./pi/1.e+06      !dsig/dtdphicm in microbarns/MeV**2/rad
c      ntup%dsigdt = sig

C--->Debug
c      write(*,*) 'sig =',sig
c      write(*,*) '====================================================='
C--->Debug

*******************************************************************************
* sigma_eepi is two-fold C.M. cross section: d2sigma/dt/dphi_cm [ub/MeV**2/rad]
* Convert from dt dphi_cm --> dOmega_lab using 'jacobian' [ub/sr]
* Convert to 5-fold by multiplying by flux factor, gtpr [1/MeV]
* to give d5sigma/dOmega_pi/dOmega_e/dE_e [ub/Mev/sr].
*******************************************************************************
c VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)

      dt_dcos_lab = 2.*(vertex%q*vertex%p%P*targ%Mtar_struck) / 
     1     ( targ%Mtar_struck+vertex%nu
     2     -vertex%q*tcos*(vertex%p%E/vertex%p%P) )
        
      jacobian = abs(dt_dcos_lab)
      main%davejac = jacobian
        
c      write(6,*)' jac ',davejac_fer,jacobian

*******************************************************************************
* Get photon flux factor
*
* DJG: Replace targ.Mtar_pion in denominator of gammaflux with more general
* DJG: efer-pfer*tfcos, for pfer =0 this reverts to old form

      gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*(s_gev-mtar_gev**2)/2./
     >     (targ%Mtar_struck)/Q2_g/(1.-epsi)
        
      davesig = gtpr*sig*jacobian
*******************************************************************************

      sigma_eepi = davesig
      peepi = sigma_eepi
        
      ntup%sigcm = sigma_eepi   !sig_cm
        
 202  format(/11X,f5.1/)
 203  format(11X,f5.0)
 204  format(6(/9X,7f8.3))
      
      return
      end


*******************************************************************************
      real*8 function peeK(vertex,main,survivalprob)

*   output:
*     peeK     !d5sigma/(dE_e'*dOmega_e'*Omega_K)     (microbarn/MeV/sr^2)

      implicit none
      include 'simulate.inc'

      type(event_main):: main
      type(event):: vertex

* NOTE: when we refer to the center of mass system, it always refers to the
* photon-NUCLEON center of mass, not the photon-NUCLEUS!  The model gives
* the cross section in the photon-nucleon center of mass frame.

      real*8 sigma_eeK          !final cross section (returned as peeK)
      real*8 pathlen,zaero,betak,gammak,p_kaon
      real*8 survivalprob

      real*8 sig219,sig,fac
      real*8 wfactor,Wset
      real*8 sigt,sigl,siglt,sigtt !components of dsigma/dt
      real*8 epsi               !epsilon of virtual photon
      real*8 gtpr               !gamma_t prime.
      real*8 tcos,tsin          !cos/sin of theta between pK and q
      real*8 tfcos,tfsin        !cos/sin of theta between pfermi and q
      real*8 pKx,pKy,pKz        !kaon momentum in lab.
      real*8 zero

      real*8 s,s_gev,t,t_gev,Q2_g !t,s, Q2 in (GeV/c)**2
      real*8 mtar_gev

      real*8 new_x_x,new_x_y,new_x_z
      real*8 new_y_x,new_y_y,new_y_z
      real*8 new_z_x,new_z_y,new_z_z
      real*8 dummy,p_new_x,p_new_y
      real*8 davesig,phipq,cospq,sinpq
      real*8 square_root,dp_dcos_num,dp_dcos_den,dp_dcos
      real*8 dp_dphi_num,dp_dphi_den,dp_dphi
      real*8 dt_dcos_lab,dt_dphi_lab,psign
      real*8 dpxdphi,dpydphi,dpxdcos,dpydcos,dpzdcos,dpzdphi
      real*8 dpxnewdphi,dpynewdphi,dpxnewdcos,dpynewdcos
      real*8 dpznewdphi,dpznewdcos
      real*8 dphicmdphi,dphicmdcos

      real*8 pbeam,beam_newx,beam_newy,beam_newz
      real*8 pKcm_newx,pKcm_newy,pKcm_newz

      real*8 dEcmdcos,dEcmdphi,dcoscmdcos,dcoscmdphi
      real*8 qx,qy,qz,px,py,pz

      real*8 t_high,t_low,delta_t,delta_Q2,Q2_tmp,t_tmp
      real*8 q2_high,q2_low

      complex*8 A1,A2,A3,A4,A12,A34,A1234
      complex*8 XL,XT,XTT,XLT
      real*8 sinterp(4)

      integer Q2_count,t_count,itab,j

      logical first

c Variables calculated in transformation to gamma-NUCLEON center of mass.
      real*8 gstar,bstar,bstarx,bstary,bstarz           !beta of boost to C.M.
      real*8 nustar,qstar,qstarx,qstary,qstarz          !q in C.M.
      real*8 eKcm,pKcm,pKcmx,pKcmy,pKcmz                !p_hadron in C.M.
      real*8 ebeamcm,pbeamcm,pbeamcmx,pbeamcmy,pbeamcmz !p_beam in C.M.
      real*8 etarcm,ptarcm,ptarcmx,ptarcmy,ptarcmz      !p_fermi in C.M.
      real*8 thetacm,phicm,phiqn,jacobian
      real*8 pm2_tmp,q2_cent_g

      logical first_call
      save first_call

      data first_call/.true./

      real*8 Q2_tab(18,130),W_tab(18,130),t_tab(18,130),stab(18,130,4)
      character*28 filename(18)
      integer nfiles, tbins(18)

*******************************************************************************
* Read model values when first called.
c Columns are: Q2, W, -t (GeV2), dsig_L, T, TL, TT (mubarn/GeV2)

      if(first_call) then
         first_call=.false.

c calculate central kinematics Q^2
         q2_cent_g = ( (Ebeam - spec%e%P*cos(spec%e%theta))**2
     c                 + (spec%e%P*sin(spec%e%theta))**2
     c                 - (Ebeam - spec%e%P)**2 ) /1.e6
c         write(6,*)' Q2_central = ',q2_cent_g

         if (which_kaon.eq.0) then ! p->K+Lambda0
c            if (q2_cent_g.lt.6.25) then
            write(6,*)'Selecting VR K+L0 model for Q^2=0.1-9.5, W=3.15'
            nfiles=18
            Wset=3.15
            filename(1)='VR/ep_peekl_q20p10_w3p15.dat'
            filename(2)='VR/ep_peekl_q20p50_w3p15.dat'
            filename(3)='VR/ep_peekl_q21p00_w3p15.dat'
            filename(4)='VR/ep_peekl_q21p50_w3p15.dat'
            filename(5)='VR/ep_peekl_q22p00_w3p15.dat'
            filename(6)='VR/ep_peekl_q22p50_w3p15.dat'
            filename(7)='VR/ep_peekl_q23p00_w3p15.dat'
            filename(8)='VR/ep_peekl_q23p50_w3p15.dat'
            filename(9)='VR/ep_peekl_q24p00_w3p15.dat'
            filename(10)='VR/ep_peekl_q24p50_w3p15.dat'
            filename(11)='VR/ep_peekl_q25p00_w3p15.dat'
            filename(12)='VR/ep_peekl_q25p50_w3p15.dat'
            filename(13)='VR/ep_peekl_q26p00_w3p15.dat'
            filename(14)='VR/ep_peekl_q26p50_w3p15.dat'
            filename(15)='VR/ep_peekl_q27p00_w3p15.dat'
            filename(16)='VR/ep_peekl_q28p50_w3p15.dat'
            filename(17)='VR/ep_peekl_q29p00_w3p15.dat'
            filename(18)='VR/ep_peekl_q29p50_w3p15.dat'
c            endif
         elseif (which_kaon.eq.1) then ! p->K+Sigma0
            write(6,*)'Selecting VR K+S0 model for Q^2=?, W=?'
            stop
         elseif (which_kaon.eq.2) then ! n->K+Lambda-
            write(6,*)'Selecting VR K+L- model for Q^2=?, W=?'
            stop
         endif

         do j=1,nfiles
            OPEN(3,FILE=filename(j),STATUS='OLD')
            write(6,80) j,filename(j)
 80         format(' Reading model file',i3,'  ',a28)
            do t_count = 1,130
               read(3,*,end=100) Q2_tab(j,t_count),
     1              W_tab(j,t_count),t_tab(j,t_count),
     2              stab(j,t_count,1), stab(j,t_count,2), 
     3              stab(j,t_count,3), stab(j,t_count,4)
               if (abs(W_tab(j,t_count)-Wset).gt.0.01) then
                  write(6,*)' W read error ',W_tab(j,t_count),Wset
                  stop
               endif
            enddo
 100        continue
            tbins(j)=t_count-1
c            write(6,*)' lines found ',tbins(j),t_tab(j,tbins(j))
            close(unit=3)
c            write(6,*) 'done'
         enddo
      endif                     !if first time.
*******************************************************************************

* Initialize some stuff.
      Q2_g = vertex%Q2/1.e6
c NOTE: phipq calculation in event.f reverted to original version.
      phipq= main%phi_pq
      mtar_gev = targ%Mtar_struck/1000.

      cospq = cos(phipq)
      sinpq = sin(phipq)

* calculate energy of initial (struck) nucleon, using the same assumptions that
* go into calculating the kaon angle/momentum (in event.f).  For A>1, the struck
* nucleon is off shell, the 2nd nucleon (always a neutron) is on shell, and has
* p = -p_fermi, and any additional nucleons are at rest.
* NOTE pfer, efer appear to be in MeV at this point.
      efer = sqrt(pfer**2+targ%Mtar_struck**2)

!      if(doing_deutpi.or.doing_hepi) then
!         efer = targ%M-sqrt(mn**2+pfer**2)
!         if(doing_hepi)efer=efer-mp
!c         mtar_offshell = sqrt(efer**2-pfer**2)
!      endif

* calculate some kinematical variables
* f's and fer indicate fermi momenta, s, star or cm CM system

      tcos = vertex%up%x*vertex%uq%x+vertex%up%y*vertex%uq%y
     1     +vertex%up%z*vertex%uq%z
      if(tcos-1..gt.0..and.tcos-1..lt.1.e-8)tcos=1.0
      tsin=sqrt(1.-tcos**2)
        
      tfcos = pferx*vertex%uq%x+pfery*vertex%uq%y+pferz*vertex%uq%z
      if(tfcos-1..gt.0..and.tfcos-1..lt.1.e-8)tfcos=1.0
      tfsin=sqrt(1.-tfcos**2)

      epsi = 1./(1.+2*(1.+vertex%nu**2/vertex%Q2)*tan(vertex%e%theta/2.)**2)

c GH:
c calculate with target nucleon at rest (as in experimental replay)
      s = (vertex%nu+targ%Mtar_struck)**2-vertex%q**2
      main%wcm = sqrt(s)
      s_gev = s/1.e6

c GH:
c calculate t the same way as in experimental replay
      pm2_tmp = (vertex%q*vertex%uq%x-vertex%p%P*vertex%up%x)**2 +
     1     (vertex%q*vertex%uq%y-vertex%p%P*vertex%up%y)**2 +
     2     (vertex%q*vertex%uq%z-vertex%p%P*vertex%up%z)**2
      t = -(mn-mp)**2 +2*targ%Mtar_struck*
     1     ( sqrt(targ%Mrec_struck**2+pm2_tmp)-targ%Mrec_struck )

      t_gev = t/1.e6               !CONVERT TO (GeV/c)**2
      main%t = t

*******************************************************************************
c VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)
* Calculate velocity of PHOTON-NUCLEON C.M. system in the lab frame. Use beta
* and gamma of the cm system (bstar and gstar) to transform particles into
* c.m. frame.  Define z along the direction of q, and x to be along the
* direction of the kaon momentum perpendicular to q.

      dummy=sqrt((vertex%uq%x**2+vertex%uq%y**2) * 
     1     (vertex%uq%x**2+vertex%uq%y**2+vertex%uq%z**2))
      new_x_x = -(-vertex%uq%x)*vertex%uq%z/dummy
      new_x_y = -(-vertex%uq%y)*vertex%uq%z/dummy
      new_x_z = ((-vertex%uq%x)**2 + (-vertex%uq%y)**2)/dummy

      dummy   = sqrt(vertex%uq%x**2 + vertex%uq%y**2)
      new_y_x =  (-vertex%uq%y)/dummy
      new_y_y = -(-vertex%uq%x)/dummy
      new_y_z =  0.0

* get beam in q+nucleon CM system.

      pbeam = sqrt(vertex%Ein**2-me**2)
      beam_newx = pbeam*new_x_z
      beam_newy = pbeam*new_y_z
      beam_newz = pbeam*vertex%uq%z

      bstar=vertex%q/(vertex%nu+targ%Mtar_struck)
      gstar=1./sqrt(1. - bstar**2)

      zero =0.e0
      bstarx = zero
      bstary = zero
      bstarz = vertex%q/(targ%Mtar_struck+vertex%nu)

* DJG: Boost virtual photon to q+nucleon CM.

      call loren(gstar,bstarx,bstary,bstarz,vertex%nu,
     >     zero,zero,vertex%q,nustar,qstarx,qstary,qstarz,qstar)

* DJG: Boost kaon to q+nucleon CM.

      pKz = vertex%p%P*tcos
      pKx = vertex%p%P*tsin*cos(phipq)
      pKy = vertex%p%P*tsin*sin(phipq)
      call loren(gstar,bstarx,bstary,bstarz,vertex%p%E,
     >     pKx,pKy,pKz,eKcm,pKcmx,pKcmy,pKcmz,pKcm)
      thetacm = acos((pKcmx*qstarx+pKcmy*qstary+pKcmz*qstarz)/pKcm/qstar)
      main%pcm = pKcm

* DJG Boost the beam to q+nucleon CM.

      call loren(gstar,bstarx,bstary,bstarz,vertex%Ein,beam_newx,
     >     beam_newy,beam_newz,ebeamcm,pbeamcmx,pbeamcmy,pbeamcmz,pbeamcm)

* Thetacm is defined as angle between pKcm and qstar.
* To get phicm, we need out of plane angle relative to scattering plane
* (plane defined by pbeamcm and qcm).  For stationary target, this plane
* does not change.  In general, the new coordinate system is defined such
* that the new y direction is given by (qcm x pbeamcm) and the new x
* is given by (qcm x pbeamcm) x qcm.

      dummy = sqrt((qstary*pbeamcmz-qstarz*pbeamcmy)**2+
     >     (qstarz*pbeamcmx-qstarx*pbeamcmz)**2
     >     +(qstarx*pbeamcmy-qstary*pbeamcmx)**2)
      new_y_x = (qstary*pbeamcmz-qstarz*pbeamcmy)/dummy
      new_y_y = (qstarz*pbeamcmx-qstarx*pbeamcmz)/dummy
      new_y_z = (qstarx*pbeamcmy-qstary*pbeamcmx)/dummy

      dummy = sqrt((new_y_y*qstarz-new_y_z*qstary)**2
     >     +(new_y_z*qstarx-new_y_x*qstarz)**2
     >     +(new_y_x*qstary-new_y_y*qstarx)**2)
      new_x_x = (new_y_y*qstarz-new_y_z*qstary)/dummy
      new_x_y = (new_y_z*qstarx-new_y_x*qstarz)/dummy
      new_x_z = (new_y_x*qstary-new_y_y*qstarx)/dummy

      new_z_x = qstarx/qstar
      new_z_y = qstary/qstar
      new_z_z = qstarz/qstar

      pKcm_newx = pKcmx*new_x_x + pKcmy*new_x_y + pKcmz*new_x_z
      pKcm_newy = pKcmx*new_y_x + pKcmy*new_y_y + pKcmz*new_y_z
      pKcm_newz = pKcmx*new_z_x + pKcmy*new_z_y + pKcmz*new_z_z

      phicm = atan2(pKcm_newy,pKcm_newx)
      if(phicm.lt.0.) phicm = 2.*3.141592654+phicm

      main%thetacm = thetacm
      main%phicm = phicm

c      write(6,*)'  '
c      write(6,*)' pfer ',pfer
c      write(6,*)' t ',t
c      write(6,*)' s ',s
c      write(6,*)' thetacm ',thetacm*180./3.14159
c      write(6,*)' phicm ',phicm*180./3.14159,phipq*180./3.14159

*******************************************************************************
* Interpolate model from tables.

      t_tmp = abs(t_gev)
      Q2_tmp = Q2_g
      if(Q2_tmp.lt.Q2_tab(1,1)) then
         write(6,*)' WARNING: model table range too high Q2 ',Q2_tmp
         Q2_tmp = Q2_tab(1,1)
      endif
      if(Q2_tmp.gt.Q2_tab(nfiles,1)) then
         write(6,*)' WARNING: model table range too low Q2 ',Q2_tmp
         Q2_tmp = Q2_tab(nfiles,1)
      endif

      do Q2_count=1,(nfiles-1)
         Q2_high=0.
         Q2_low=0.
         t_high=0.
         t_low=0.
         if((Q2_tmp.ge.Q2_tab(Q2_count,1)).and.
     1        (Q2_tmp.lt.Q2_tab(Q2_count+1,1))) then
            Q2_high = Q2_tab(Q2_count+1,1)
            Q2_low = Q2_tab(Q2_count,1)
            delta_Q2 = (Q2_high-Q2_low)
c            write(6,*)' Q2 ',Q2_tab(Q2_count,1),Q2_tmp

            if(t_tmp.lt.t_tab(Q2_count,1)) then
               t_tmp = t_tab(Q2_count,1)
               if (abs(t_tmp-t_tab(Q2_count,1)).gt.0.1) then
                  write(6,*)' WARNING: model table range too high -t ',t_tmp
               endif
            endif
            if(t_tmp.gt.t_tab(Q2_count,tbins(Q2_count)-1)) then
               if (abs(t_tmp-t_tab(Q2_count,tbins(Q2_count)-1)).gt.0.1) then
                  write(6,*)' WARNING: model table range too low -t ',t_tmp
               endif
               t_tmp = t_tab(Q2_count,tbins(Q2_count)-1)
            endif
            do t_count=1,(tbins(Q2_count)-1)
               if((t_tmp.ge.t_tab(Q2_count,t_count)).and.
     1              (t_tmp.lt.t_tab(Q2_count,t_count+1))) then
                  t_high = t_tab(Q2_count,t_count+1)
                  t_low =  t_tab(Q2_count,t_count)
                  delta_t = t_high-t_low
c                  write(6,*)' t ',t_tab(Q2_count,t_count),t_tmp

               do itab=1,4
                 A1 = stab(Q2_count,t_count,itab)
                 A2 = stab(Q2_count+1,t_count,itab)
                 A3 = stab(Q2_count,t_count+1,itab)
                 A4 = stab(Q2_count+1,t_count+1,itab)
                 A12= (A1*(Q2_high-Q2_tmp)+A2*(Q2_tmp-Q2_low))/delta_Q2
                 A34= (A3*(Q2_high-Q2_tmp)+A4*(Q2_tmp-Q2_low))/delta_Q2
                 A1234 = (A12*(t_high-t_tmp)+
     1                A34*(t_tmp-t_low))/delta_t
                 sinterp(itab) = A1234
               enddo
               endif               !t check
            enddo                  !t
         endif                     !Q2 check
      enddo                        !Q2
     
      sigl = sinterp(1)
      sigt = sinterp(2)
      siglt = sinterp(3)
      sigtt = sinterp(4)

      sig219=(sigt+main%epsilon*sigl+main%epsilon*cos(2.*phicm)*sigtt
     >     +sqrt(2.0*main%epsilon*(1.+main%epsilon))*cos(phicm)*siglt)/1.e0
      
c now convert to different W
c W dependence given by 1/(W^2-M^2)^2
     
      wfactor=(Wset**2-(targ%Mtar_struck/1000.)**2)**2/
     1     (s_gev-(targ%Mtar_struck/1000.)**2)**2
          
      sig=sig219*wfactor
      sigl=sigl*wfactor
      sigt=sigt*wfactor
      sigtt=sigtt*wfactor
      siglt=siglt*wfactor

C--->Debug
c      write(*,*) 's =',s_gev
c      write(*,*) 'wfactor =',wfactor
c      write(*,*) 'sig =',sig
c      write(*,*) 'sigL  =',sigL
c      write(*,*) 'sigT  =',sigT
c      write(*,*) 'sigLT =',sigLT
c      write(*,*) 'sigTT =',sigTT
c      write(*,*) '-----------------------------------------------------'
C--->Debug

      sig=sig/2./pi/1.e+06      !dsig/dtdphicm in microbarns/MeV**2/rad
c      ntup%dsigdt = sig

C--->Debug
c      write(*,*) 'sig =',sig
c      write(*,*) '====================================================='
C--->Debug

*******************************************************************************
* sigma_eeK is two-fold C.M. cross section: d2sigma/dt/dphi_cm [ub/MeV**2/rad]
* Convert from dt dphi_cm --> dOmega_lab using 'jacobian' [ub/sr]
* Convert to 5-fold by multiplying by flux factor, gtpr [1/MeV]
* to give d5sigma/dOmega_K/dOmega_e/dE_e [ub/Mev/sr].
*******************************************************************************
c VERSION WHERE TARGET NUCLEON IS AT REST (AS IN EXPERIMENTAL REPLAY)

      dt_dcos_lab = 2.*(vertex%q*vertex%p%P*targ%Mtar_struck) / 
     1     ( targ%Mtar_struck+vertex%nu
     2     -vertex%q*tcos*(vertex%p%E/vertex%p%P) )
        
      jacobian = abs(dt_dcos_lab)
      main%davejac = jacobian
        
c      write(6,*)' jac ',davejac_fer,jacobian

*******************************************************************************
* Get photon flux factor
*
* DJG: Replace targ.Mtar_kaon in denominator of gammaflux with more general
* DJG: efer-pfer*tfcos, for pfer =0 this reverts to old form

      gtpr = alpha/2./(pi**2)*vertex%e%E/vertex%Ein*(s_gev-mtar_gev**2)/2./
     >     (targ%Mtar_struck)/Q2_g/(1.-epsi)
        
      davesig = gtpr*sig*jacobian
*******************************************************************************

      sigma_eeK = davesig
      peeK = sigma_eeK
        
      ntup%sigcm = sigma_eeK   !sig_cm
        
c If doing_decay=.false., generate survival probability for main.weight.
c main.FP.p.path is dist. to back of detector, want decay prob. at aerogel.
C NOTE THAT ZAERO IS TAKEN WITH RESPECT TO THE POSITION AT WHICH PATHLEN
C IS CALCULATED (USUALLY THE BACK OF THE CALORIMETER).  IF THE DRIFTS IN
C MC_SOS_HUT ARE CHANGED, THEN THE STARTING POINT MAY BE DIFFERENT AND THESE
C NUMBERS MAY BE WRONG.  AERO BETWEEN S2Y AND S2X IN SOS.

C Beta/Gamma for decay need to use momentum after radiation/eloss, not vertex
C momentum.  Get particle momentum from main%SP%p%delta

	if (.not.doing_decay) then
	  if (hadron_arm.eq.1) then
	    zaero = 0.			!no aerogel yet, use full length.
	  else if (hadron_arm.eq.2) then
*	    zaero = -76.		!aero at 270cm,last project=346(cal).
	    zaero = -82.8		!From Rick: aero at 263.2cm,last project=346(cal).
	  else if (hadron_arm.eq.3) then
	    zaero = -183.		!aero at 130cm,last project=313(S2)
	  else if (hadron_arm.eq.4) then
	    zaero = -183.
	  endif
	  pathlen = main%FP%p%path + zaero*(1+main%FP%p%dx**2+main%FP%p%dy**2)
	  p_kaon = spec%p%P*(1.+main%SP%p%delta/100.)
	  betak = spec%p%P/sqrt(spec%p%P**2+Mh2)
	  gammak = 1./sqrt(1.-betak**2)
	  survivalprob = 1./exp(pathlen/(ctau*betak*gammak))
	  decdist = survivalprob		!decdist in ntuple
	endif

 202  format(/11X,f5.1/)
 203  format(11X,f5.0)
 204  format(6(/9X,7f8.3))
      
      return
      end


*******************************************************************************
	real*8 function sig_blok(thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev,which_pion)

* Purpose:
* This routine calculates p(e,e'pi+)n cross sections from a fit to the data of
* Brauel et al., Z.Phys.C. 3(1979)101.
* Fit gives dsigma/dt/dphi_cm, which is returned as sig_blok [ub/MeV^2-rad].

	implicit none
	include 'constants.inc'

	real*8 sig219,sig
	real*8 sigt,sigl,siglt,sigtt	!components of dsigma/dt
	
	real*8 mrho_emp
	real*8 fpi,fpi2

	real*8 thetacm,phicm,t,q2_gev,s_gev,eps,mtar_gev
	integer which_pion

* Fits to p(e,e'pi+)n cross-sections
* HPB: cross sections for pi-fits to Brauel's n(e,e'pi-)p
* HPB: dsigma/dt cross sections at W=2.19 and Q^2=0.70 MODIFIED BY HAND!!!

	  sigl = 27.8*exp(-11.5*abs(t))
	  sigt = 10.0*(5.*abs(t))*exp(-5.*abs(t))
	  siglt= 0.0*sin(thetacm)
	  sigtt= -(4.0*sigl+0.5*sigt)*(sin(thetacm))**2

	  if (which_pion.eq.1 .or. which_pion.eq.11) then	!pi-
	    sigt =sigt *0.25*(1.+3.*exp(-10.*abs(t)))
	    sigtt=sigtt*0.25*(1.+3.*exp(-10.*abs(t)))
	  endif

* GMH: Now scale sigl by Q2*pion_formfactor
* HPB: parametrization of pion formfactor in the region 0.5<Q2<1.8

	  mrho_emp = 0.712   !DETERMINED BY EMPIRICAL FIT TO FORMFACTOR (GeV/c**2)

	  fpi=1./(1.+1.65*q2_gev+0.5*q2_gev**2)
	  fpi2=fpi**2

* HPB: now convert to different Q^2
* HPB: s_L follows Q^2*Fpi^2; s_T and s_TT go as 1/(0.3+Q^2)?
* HPB: factor 0.1215 therefore is value of q2*fpi**2 for q2=0.7

	  sigl=sigl*(fpi2*q2_gev)/0.1215
	  sigt=sigt/(0.3+q2_gev)
	  sigtt=sigtt/(0.3+q2_gev)

	  sig219=(sigt+eps*sigl+eps*cos(2.*phicm)*sigtt
     >		+sqrt(2.0*eps*(1.+eps))*cos(phicm)*siglt)/1.e0

* GMH: Brauel scaled all his data to W=2.19 GeV, so rescale by 1/(W**2-mp**2)**2
* HPB: factor 15.333 therefore is value of (W**2-mp**2)**2 at W=2.19

	  sig=sig219*15.333/(s_gev-mtar_gev**2)**2
	  sig=sig/2./pi/1.d+06   !dsig/dtdphicm in microbarns/MeV**2/rad

	  sig_blok = sig

	  return
	end
