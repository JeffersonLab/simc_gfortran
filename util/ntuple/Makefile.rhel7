## This makefile must be executed with gmake (gnu make).

## This tells make not to delete these target files on error/interrupt (see man page)
.PRECIOUS: *.o 
RM        = rm -f 
SHELL     = /bin/sh

OBJ	= make_ntuple.o
 
my_objs	=  $(OBJ)

my_deps = $(my_objs:.o=.d)

MYOS := $(subst -,,$(shell uname))
#CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
CERNLIBS = -Wl,-static -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib -Wl,-dy


#For use with gfortran compiler
# -fno-automatic - all program storage treated as static
ifeq ($(MYOS),Linux)
  FFLAGSA=-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib64 
# 64 vs 32 bit
#        -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
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

none: make_ntuple $(my_deps)

all: make_ntuple  $(my_deps)

include $(my_deps)

make_ntuple: $(my_objs) Makefile
	$(F77) $(OSF_SHARED) -o $@ $(FFLAGS) $(my_objs) $(OTHERLIBS)

clean:
	$(RM) *.[od] make_ntuple

