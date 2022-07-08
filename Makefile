## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## CEBAF DEFAULT SETUP FLAGS:
simcdir = .

## THE REST SHOULD BE OK WITHOUT MODIFICATION.

## This tells make not to delete these target files on error/interrupt (see man page)
.PRECIOUS: *.o sos/*.o hms/*.o hrsl/*.o hrsr/*.o shms/*.o calo/*.o

RM        = rm -f 
SHELL     = /bin/sh
S	= $(simcdir)/sos/
H	= $(simcdir)/hms/
L	= $(simcdir)/hrsl/
R	= $(simcdir)/hrsr/
A	= $(simcdir)/shared/
SH	= $(simcdir)/shms/
T       = $(simcdir)/cteq5/
C       = $(simcdir)/calo/
CH	= $(simcdir)/cern/
D       = $(simcdir)/fdss/

OBJ1	= target.o brem.o gauss1.o NtupleInit.o NtupleClose.o enerloss_new.o
OBJ2	= radc.o init.o dbase.o physics_kaon.o physics_pion.o physics_delta.o physics_proton.o loren.o sf_lookup.o
OBJ3    = semi_physics.o rho_physics.o rho_decay.o generate_rho.o trg_track.o semi_dilution.o
OBJ4	= results_write.o event.o call_ranlux.o jacobians.o F1F2IN21_v1.0.o
OBJ5	= $(A)musc.o $(A)musc_ext.o $(A)project.o $(A)transp.o
OBJ6	= $(A)rotate_haxis.o $(A)rotate_vaxis.o $(A)locforunt.o
OBJ7	= $(H)mc_hms.o $(H)mc_hms_hut.o $(H)mc_hms_recon.o
OBJ8	= $(S)mc_sos.o $(S)mc_sos_hut.o $(S)mc_sos_recon.o
OBJ9	= $(R)mc_hrsr.o $(R)mc_hrsr_hut.o $(R)mc_hrsr_recon.o
OBJA	= $(L)mc_hrsl.o $(L)mc_hrsl_hut.o $(L)mc_hrsl_recon.o
OBJB	= $(SH)mc_shms.o $(SH)mc_shms_hut.o $(SH)mc_shms_recon.o
OBJC    = $(T)Ctq5Pdf.o
OBJD    = $(C)mc_calo.o $(C)mc_calo_recon.o
OBJCH   = $(CH)lfit.o $(CH)ranlux.o $(CH)fint.o $(CH)kerset.o $(CH)abend.o
OBJF   = $(D)fdss.o
 
my_objs	=  $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJA) $(OBJB) $(OBJC) $(OBJD) $(OBJCH) $(OBJF)

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
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT) -lctp \
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

$(S)/%.o: $(S)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(H)/%.o: $(H)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(L)/%.o: $(L)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(R)/%.o: $(R)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(SH)/%.o: $(SH)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(CH)/%.o: $(SH)/%.f
	$(F77) $(FFLAGS) -c $< -o $@


DEPEND_RULE = ( cat $< |  sed -n -e \
	"s|^[ 	]*[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ 	]*['\"]|$@: $(@D)/|p" | \
	sed -e "s|['\"].*$$||" | \
	sed -e 's|.d:|.o:|') > $@

%.d: %.f
	$(DEPEND_RULE)

none: simc $(my_deps)

all: simc  $(my_deps)

include $(my_deps)

simc: simc.o $(my_objs) Makefile CTP/O.Linux/Linux/lib/libctp.a
	$(F77) $(OSF_SHARED) -o $@ $(FFLAGS) $(my_objs) simc.o $(OTHERLIBS)

# Looks like by default simulate_init.inc and sturcutres_init.inc are
# not used, but we make rule just in case.
simulate_init.inc: simulate.inc
	( cat $< | sed -e "s/structures.inc/structures_init.inc/") > $@

structures_init.inc: structures.inc
	(cat $< |\
	sed -e "s/^[ 	]*real\*8[ 	]*min[ 	]*,[ 	]*max[ 	]*$$/		real*8::	min=-1.d10, max=1d10/" \
	-e "s/^[ 	]*real\*8[ 	]*lo[ 	]*,[ 	]*hi[ 	]*$$/		real*8::	lo=1.d10, hi=-1d10/") > $@

CTP/O.Linux/Linux/lib/libctp.a: 
	make -C CTP
# These routines have HP problems, and need to be compiled without optimization.

#simc.o: simc.f
#	$(F77) $(FFLAG1) simc.f
#
#init.o: init.f
#	$(F77) $(FFLAG1) init.f
#
#event.o: event.f
#	$(F77) $(FFLAG1) event.f
#
#dbase.o: dbase.f
#	$(F77) $(FFLAG1) dbase.f
#
#physics.o: physics.f
#	$(F77) $(FFLAG1) physics.f
#
#results_write.o: results_write.f
#	$(F77) $(FFLAG1) results_write.f


clean:
	$(RM) *.[od] $(H)*.[od] $(S)*.[od] $(L)*.[od] $(R)*.[od] $(SH)*.[od] $(A)*.[od] $(T)*.[od] $(C)*.[od] $(CH)*.[od] simc

real_clean:
	$(RM) *.[od] $(H)*.[od] $(S)*.[od] $(L)*.[od] $(R)*.[od] $(SH)*.[od] $(A)*.[od] $(T)*.[od] $(C)*.[od] simc
	rm -r CTP/O.$(MYOS)
