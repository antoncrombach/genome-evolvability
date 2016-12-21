#
# Simple C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#
# History
# 2004-08-04	AC	Creation, using implicit rules
# 2004-08-05	AC	Dependencies automatically generated,
# 			not included if goal is `clean' or `realclean'
# 2004-11-17	AC	Adjusting for different placement of object files	
# 2004-11-18	AC	Introduced the use of 'vpath'
# 2004-11-19	AC	Removed almost everything, the makefile serves as
# 			a shortcut to the makefile in obj/ 
#

OBJPATH = obj/
BINPATH = bin/

# Targets
all: 
	@cd $(OBJPATH); \
	make all

slib:
	@cd $(OBJPATH); \
	make slib

fluq:
	@cd $(OBJPATH); \
	make fluq

.PHONY: clean realclean distclean
clean:
	@cd $(OBJPATH); make clean

distclean: 
	@cd $(OBJPATH); make distclean; \
	cd ../$(BINPATH); rm -f *

realclean:
	@cd $(OBJPATH); make realclean

