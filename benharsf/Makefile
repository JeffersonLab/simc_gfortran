RM     = rm -f
SHELL  = /bin/sh

ifndef MYOS
  MYOS := $(subst -,,$(shell uname))
endif

ifeq ($(MYOS),HPUX)
#  FFLAGS=+ppu -O -C +es -I$(MYINC) +Obb1000
  FFLAGS=-O -C +es -I$(MYINC) +Obb1000
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = 
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS= 
endif

ifeq ($(MYOS),OSF1)
  FFLAGS=-r8 -O -C -extend_source -fpe -warn argument_checking \
	-warn declarations -warn truncated_source -warn unused -check bounds
  LDFLAGS=
  OTHERLIBS =
endif


subroutines = 

all: readc12 readfe56 readau197 readhe4 readhe3mod readhe3par readnm

readc12: $(subroutines) readc12.o Makefile
	 $(F77) -o $@ $(FFLAGS) readc12.o $(subroutines) $(OTHERLIBS)

readfe56: $(subroutines) readfe56.o Makefile
	 $(F77) -o $@ $(FFLAGS) readfe56.o $(subroutines) $(OTHERLIBS)

readau197: $(subroutines) readau197.o Makefile
	 $(F77) -o $@ $(FFLAGS) readau197.o $(subroutines) $(OTHERLIBS)

readhe4: $(subroutines) readhe4.o Makefile
	 $(F77) -o $@ $(FFLAGS) readhe4.o $(subroutines) $(OTHERLIBS)

readhe3mod: $(subroutines) readhe3mod.o Makefile
	 $(F77) -o $@ $(FFLAGS) readhe3mod.o $(subroutines) $(OTHERLIBS)

readhe3par: $(subroutines) readhe3par.o Makefile
	 $(F77) -o $@ $(FFLAGS) readhe3par.o $(subroutines) $(OTHERLIBS)

readnm: $(subroutines) readnm.o Makefile
	 $(F77) -o $@ $(FFLAGS) readnm.o $(subroutines) $(OTHERLIBS)

clean:
	 rm *.o readc12 readfe56 readau197 readhe4 readhe3mod readhe3par readnm
