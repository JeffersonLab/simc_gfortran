## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## ARGONNE DEFAULT SETUP FLAGS:
#simcdir = .
#Csoft = /disk1/users/reinhold/Csoft/05Dec1996

## CEBAF DEFAULT SETUP FLAGS:
#simcdir = .
#Csoft = /group/hallc/Csoft/Analyzer

## Setup for Basel
simcdir = .
Csoft = /usr/users/rohe/gen/Csoft

## THE REST SHOULD BE OK WITHOUT MODIFICATION.

## This tells make not to delete these target files on error/interrupt (see man page)
.PRECIOUS: *.o sos/*.o hms/*.o

RM        = rm -f 
SHELL     = /bin/sh
S	= $(simcdir)/sos/
H	= $(simcdir)/hms/

OBJ1	= target.o brem.o gauss1.o NtupleInit.o NtupleClose.o enerloss_new.o
OBJ2	= radc.o init.o dbase.o physics.o loren.o trg_track.o
OBJ3	= results_write.o event.o mt19937.o eekeek.o
OBJ4	= $(S)musc.o $(S)musc_ext.o $(S)project.o $(S)transp.o
OBJ5	= $(S)rotate_haxis.o $(S)locforunt.o $(S)stringlib.o $(S)lowcase.o
OBJ6	= $(H)mc_hms.o $(H)mc_hms_hut.o $(H)mc_hms_recon.o $(H)check_dipole.o
OBJ7	= $(S)mc_sos.o $(S)mc_sos_hut.o $(S)mc_sos_recon.o
my_objs	=  $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7)

MYOS := $(subst -,,$(shell uname))
CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

ifeq ($(MYOS),HPUX)
  ifneq (,$(findstring 09,$(shell uname -r)))
    HPUXVERSION := 09
  else
    HPUXVERSION := 10
  endif
  LIBROOT = $(Csoft)/../$(MYOS)$(HPUXVERSION)/lib
else
  LIBROOT = $(Csoft)/../$(MYOS)/lib
endif

ifeq ($(MYOS),HPUX)
  CERN_ROOT = /site/cernlib/hp700_ux90/96a
#  CERN_ROOT = /site/cernlib/hp_ux102/97a
#  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8 -G
  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8
  FFLAG1=+U77 +ppu -C +e +es +FPVZOu -c -R8
  OTHERLIBS = \
        -Wl,-L$(LIBROOT) -lctp \
        -Wl,-L$(CERN_ROOT)/lib $(CERNLIBS) \
        -Wl,-L/usr/lib/X11R5 -lX11 -lm
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  FFLAG1=$(FFLAGS)
  LDFLAGS=
  OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS)
endif

ifeq ($(MYOS),Linux)
  FFLAGS=-I$(Csoft)/SRC/INCLUDE
  FFLAG1=$(FFLAGS)
  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -lc -lm
  F77=Af77
endif

ifeq ($(MYOS),SunOS)
  CERN_ROOT = /site/cernlib/sun4_solaris2/97a
  FFLAGS=-e -O -I$(Csoft)/SRC/INCLUDE -ftrap=no%underflow
  FFLAG1=$(FFLAGS) -c $<
  ifeq ($(OSTYPE),SunOS4)
    OTHERLIBS = -L$(LIBROOT) -lctp -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lX11
  else
    OTHERLIBS = -L$(LIBROOT) -lctp -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lsocket -lX11
  endif
endif

ifeq ($(MYOS),AIX)
  F77=f77
  FFLAGS=-g -qfixed=132 -qextname -O -I$(Csoft)/SRC/INCLUDE
  FFLAG1=$(FFLAGS)
  OTHERLIBS = -L$(LIBROOT) -lctp -L$(CERN_ROOT)/lib $(CERNLIBS) -lX11
endif

ifeq ($(MYOS),OSF1)
  F77=f77
#  CERN_ROOT = /disk1/lib/cern/new
  LIBROOT = $(Csoft)/OSF1/lib
  FFLAGS= -r8 -extend_source -Wl,-taso -I -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused 
  FFLAG1= -C -r8 -extend_source -Wl,-taso -I -g -c -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused 
  LDFLAGS= 
  OTHERLIBS = -Wl,-L$(LIBROOT) -lctp -Wl,-L$(CERN_ROOT)/lib \
        $(CERNLIBS) -Wl,-L/usr/lib/X11R5 -lX11 -lm 
endif

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

$(S)/%.o: $(S)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

$(H)/%.o: $(H)/%.f
	$(F77) $(FFLAGS) -c $< -o $@

none: simc

all: simc

simc: simc.o $(my_objs) Makefile
	$(F77) -o $@ $(FFLAGS) $(my_objs) simc.o $(OTHERLIBS)

# These routines have HP problems, and need to be compiled without optimization.

simc.o: simc.f
	$(F77) $(FFLAG1) simc.f

init.o: init.f
	$(F77) $(FFLAG1) init.f

event.o: event.f
	$(F77) $(FFLAG1) event.f

dbase.o: dbase.f
	$(F77) $(FFLAG1) dbase.f

physics.o: physics.f
	$(F77) $(FFLAG1) physics.f

results_write.o: results_write.f
	$(F77) $(FFLAG1) results_write.f


clean:
	$(RM) *.o $(H)*.o $(S)*.o simc
