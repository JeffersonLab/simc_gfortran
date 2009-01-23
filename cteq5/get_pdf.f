      program get_pdf

      implicit none

      integer Iset,Ipart
      real*8 x,Q2,Q
      real*8 u, ubar, d, dbar
      real*8 Ctq5Pdf

      external Ctq5Pdf

      write(6,*) 'Enter Parton Distribution Function Set (1=CTEQ5M)'
      read(5,*) Iset

      call SetCtq5(Iset)

c      write(6,*) 'which quark? (u=1,d=2,u_bar=-1,d_bar=-2)'
c      read(5,*) Ipart

      write(6,*) 'Xbj?'
      read(5,*) x

      write(6,*) 'Q2?'
      read(5,*) Q2

      Q = sqrt(Q2)

      Ipart=1
      u = Ctq5pdf( Ipart,x,Q)

      Ipart=-1
      ubar = Ctq5pdf( Ipart,x,Q)

      Ipart=2
      d = Ctq5pdf( Ipart,x,Q)

      Ipart=-2
      dbar = Ctq5pdf( Ipart,x,Q)


      write(6,*) 'PDFs:'
      write(6,*) 'u(x,Q2)=   ',u
      write(6,*) 'ubar(x,Q2)=',ubar
      write(6,*) 'd(x,Q2)=   ',d
      write(6,*) 'dbar(x,Q2)=',dbar

      end
