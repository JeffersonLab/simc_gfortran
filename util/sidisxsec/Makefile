## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## CEBAF DEFAULT SETUP FLAGS:
simcdir = ../../

## THE REST SHOULD BE OK WITHOUT MODIFICATION.

## This tells make not to delete these target files on error/interrupt (see man page)
.PRECIOUS: *.o sos/*.o hms/*.o hrsl/*.o hrsr/*.o shms/*.o calo/*.o

RM        = rm -f 
SHELL     = /bin/sh
T       = $(simcdir)/cteq5/
D       = $(simcdir)/fdss/

OBJC    = $(T)Ctq5Pdf.o
OBJF   = $(D)fdss.o
MODULES	= $(simcdir)modules.f $(simcdir)histograms_module.f
 
my_objs	=  $(OBJC) $(OBJF)

my_deps = $(my_objs:.o=.d)

MYOS := $(subst -,,$(shell uname))
#CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
#CERNLIBS = -Wl,-static -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib -ljetset74 -Wl,-dy


#For use with gfortran compiler
# -fno-automatic - all program storage treated as static
ifeq ($(MYOS),Linux)
  LIBROOT = CTP/O.Linux/Linux/lib
# JLab
#  CERN_ROOT = /apps/cernlib/i386_fc8/2005
# 32 bit, standard Fedora distributuion
#  CERN_ROOT = /usr/lib/cernlib/2006
# 64 bit, standard Fedora distributuion
#  CERN_ROOT =  /usr/lib64/cernlib/2006 
  FFLAGSA=-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  INCLUDES=-I$(simcdir) 
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT)\
        -L/usr/lib64 
# 64 vs 32 bit
#        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
  FC  := gfortran
  F77 := gfortran
endif

# Mac OSX, Leopard (checked on version 10.5.6)
# Tested using patched gfortran compiler from www-jlc.kek.jp/~fujiik/macosx/10.5.X/HEPonX
# I had trouble compiling geant stuff using the deafult gfortran that shipped with
# my version of Leopard - YMMV (your mileage may vary)
# Only change needed was to add  the -fno-range-check flag: there was some problem
# with the default integer size, specifically in mt19937.f.
# Note that the CTP libraries still end up in the O.Linux directory...
ifeq ($(MYOS),Darwin)
  LIBROOT = CTP/O.Linux/Linux/lib
  CERN_ROOT = /apps/cernlib/i386_fc8/2005
  FFLAGSA=-O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
  FC  := gfortran
  F77 := gfortran
endif

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

DEPEND_RULE = ( cat $< |  sed -n -e \
	"s|^[ 	]*[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ 	]*['\"]|$@: $(@D)/|p" | \
	sed -e "s|['\"].*$$||" | \
	sed -e 's|.d:|.o:|') > $@

%.d: %.f
	$(DEPEND_RULE)

none: mods calc_semi_xsec $(my_deps)

all: mods calc_semi_xsec $(my_deps)

include $(my_deps)

mods: 
	$(F77) -c -ffixed-line-length-132 $(MODULES)

calc_semi_xsec: calc_semi_xsec.o $(simcdir)semi_physics.o $(simcdir)F1F2IN21_v1.0.o $(OBJC) $(OBJF) Makefile
	$(F77)  $(OSF_SHARED) -o $@ $(FFLAGS) $(simcdir)semi_physics.o  $(simcdir)F1F2IN21_v1.0.o $(OBJC) $(OBJF) calc_semi_xsec.o $(OTHERLIBS)

# Looks like by default simulate_init.inc and structures_init.inc are
# not used, but we make rule just in case.
simulate_init.inc: simulate.inc
	( cat $< | sed -e "s/structures.inc/structures_init.inc/") > $@

structures_init.inc: structures.inc
	(cat $< |\
	sed -e "s/^[ 	]*real\*8[ 	]*min[ 	]*,[ 	]*max[ 	]*$$/		real*8::	min=-1.d10, max=1d10/" \
	-e "s/^[ 	]*real\*8[ 	]*lo[ 	]*,[ 	]*hi[ 	]*$$/		real*8::	lo=1.d10, hi=-1d10/") > $@


clean:
	$(RM) *.[od] $(H)*.[od] $(S)*.[od] $(L)*.[od] $(R)*.[od] $(SH)*.[od] $(A)*.[od] $(T)*.[od] $(C)*.[od] $(CH)*.[od] *.mod simc

