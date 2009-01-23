.DELETE_ON_ERROR: ;

include Makefile.variables

MAKEFILENAME = Makefile.Unix


ifeq ($(MYOS),Linux)
  SUBDIRNAME = O.Linux
  OSEXT =
endif


# OSEXT is either TEST or empty
# If you want to use LinuxTEST as shown above, just pass OSEXT=TEST
# after 'make' on the command line, and leave MYOS as the default (Linux)
# we don't need to do anything with OSEXT here 
#export OSEXT #either TEST or empty
#export OSVERSION #not used since we probably can't compile on SunOS anyway

.PHONY: all
all:  $(SUBDIRNAME) $(SUBDIRNAME)/Makefile
	$(MAKE) -C $(SUBDIRNAME)

$(MAKEFILENAME):
	(cd $(SUBDIRNAME) ; ln -f -s ../$(MAKEFILENAME) Makefile)
$(SUBDIRNAME)/Makefile: $(MAKEFILENAME)
	(cd $(SUBDIRNAME) ; ln -f -s ../$(MAKEFILENAME) Makefile)

$(SUBDIRNAME):
	mkdir $(SUBDIRNAME)
