simc_gfortran - Hall C/A Physics Monte Carlo
============================================

SIMC (simc_gfortran) is the standard Hall C Monte Carlo for
coincidence reactions, written in FORTRAN.

SIMC has recently been modified to remove all explcit dependence on 
CERNLIB.  There are now 3 ways to generate output for SIMC.

1. Run simc as normal (use "run_simc" script or interactively).  This will 
generate a fortran ".bin" file in the worksim directory. Probably not too 
useful unless your analysis code is in fortran.

2. Run simc and then use the helper application to convert the default output
to a PAW ntuple (NOTE: this WILL require use of cernlib - on the JLab 
machines type "setup cernlib/2005"). The helper application is in util/ntuple 
and must be compiled separately from simc.  Running SIMC + the ntuple 
application can be done in one command using "run_simc_ntup".

2. Run simc and then use the helper application to convert the default output
to a root tree (this program is also  in fortran, but cernlib is not needed). 
The helper application is in util/root_tree and must also be compiled 
separately from simc.  Running SIMC + the root-tree application can be done 
in one command using "run_simc_tree".


Features
--------

* SIMC simulates the optics (using COSY models) and apertures of the
  Hall C spectrometers (HMS, SOS, SHMS) and other spectrometers at
  [Jefferson Lab](http://www.jlab.org/) (HRS's, BigCal, ...)
* Radiative effects, multiple scattering, ionization energy loss and
  particle decay are included
* Simple presecriptions are available for Final State Interactions,
  Coulomb Corrections and other effects.

Reactions
---------

SIMC has physics models for the following reactions.
* Elastic and quasi-elastic scatering: H(e,e'p), A(e,e'p)
* Exclusive pion production: H(e,e'pi+)n, A(e,e'pi+/-)
(quasifree or coherent)
* Kaon electroproduction: H(e,e'K+)Lambda,Sigma, A(e,e'K+/-),A(e,e'K-)
* Semi-inclusive pion production: H(e,e'pi+/-)X, D(e,e'pi+/-)X
* Semi-inclusive kaon production: H(e,e'K+/-)X, D(e,e'K+/-)X
* Diffractive rho production: H(e,e'rho->pi+ pi-)p, D(e,e'rho->pi+
pi-)

SIMC is NOT
-----------
* Not a full detector response simulation a la GEANT/GEANT4
* Does NOT simulate a large class of processes simultaneously to
gerate backgrounds (like Pythia for example)
* Not a generic event generator.  Processes are generated over a
limited phase space matching the spectrometer acceptances
* Not hard to modify

Overview
--------

An overview of SIMC can be found in
[this presentation at the 2009 Hall A Collaboration
Meeting](http://hallaweb.jlab.org/collab/meeting/2009-winter/talks/Analysis%20Workshop%20--%20Dec%2014/simc_overview.pdf)

Documentation
-------------

For more information, see the [SIMC Monte Carlo page in the Hall C
Wiki](https://hallcweb.jlab.org/wiki/index.php/SIMC_Monte_Carlo)
