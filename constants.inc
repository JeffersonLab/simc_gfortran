! CONSTANTS.INC

! Note: All angles are in radians
!     All distances are in cm
!     All energies (momenta,masses) are in MeV
!     All deltas are in percent
!     All densities (thicknesses) are in g/cm3 (g/cm2)
!     All B fields are in kG
!     All roads lead to infinity

! ... some physical constants
       real*8 Me, Me2, Mp, Mp2, Mn, Mn2, Mpi, Mpi2, Mpi0, Mpi02
       real*8 Meta, Meta2, Mk, Mk2, Mrho, Mrho2, MrhoW
       real*8 Momega, Momega2, MomegaW, Md, Md2
       real*8 Metap, Metap2, Mphi, Mphi2, MphiW
       real*8 Mlambda, Msigma0, Msigma_minus
       real*8 MDelta, Delta_width
       real*8 amu, hbarc, pi, twopi, alpha, alpi, degrad, euler
       parameter (Me=0.51099906)
       parameter (Me2=Me**2)
       parameter (Mpi=139.57018)
       parameter (Mpi2=Mpi**2)
       parameter (Mpi0=134.9766)
       parameter (Mpi02=Mpi0**2)
       parameter (Mk=493.677)
       parameter (Mk2=Mk**2)
       parameter (Meta=547.86)    !2014 PDG value
       parameter (Meta2=Meta**2)
       parameter (Mrho=775.5)     !2006 PDG value
       parameter (Mrho2=Mrho**2)
       parameter (MrhoW=149.4)    !2006 PDG value
       parameter (Momega=782.65)  !2006 PDG value
       parameter (Momega2=Momega**2)
       parameter (MomegaW=8.49)   !2006 PDG value
       parameter (Mp=938.27231)
       parameter (Mp2=Mp**2)
       parameter (Mn=939.56563)
       parameter (Mn2=Mn**2)
       parameter (Metap=957.78)   !2014 PDG value
       parameter (Metap2=Metap**2)
       parameter (Mphi=1019.455)
       parameter (Mphi2=Mphi**2)
       parameter (MphiW=4.26)
       parameter (Mlambda=1115.68)
       parameter (Msigma0=1192.64)
       parameter (Msigma_minus=1197.45)
       parameter (MDelta = 1232.0)
       parameter (Delta_width = 117.0)
       parameter (Md=1875.613)
       parameter (Md2=Md**2)
       parameter (amu=931.49432)
       parameter (hbarc=197.327053)
       parameter (pi=3.141592653589793)
       parameter (twopi=2.*pi)
       parameter (alpha=1./137.0359895)
       parameter (alpi=alpha/pi)
       parameter (degrad=180./pi)
       parameter (euler=0.577215665)
