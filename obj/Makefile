#
# General C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#
# History
# 2004-08-04	AC	Creation, using implicit rules
# 2004-08-05	AC	Dependencies automatically generated,
# 			not included if goal is `clean' or `realclean'
# 2004-11-17	AC	Adjusting for different placement of object files	
# 2004-11-18	AC	Introduced the use of 'vpath'
#

# Program vars
#CXX = g++
#CXX= icpc
#CXX = /usr/bin/g++34
CXX = /home/anton/local/bin/g++4
# optimising flags
CPPFLAGS = -Wall -O3 -march=pentium4 -mmmx -msse -msse2 -fprefetch-loop-arrays -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
#CPPFLAGS = -O2 -axN -ipo -pc32 -mp1

# debugging
#CPPFLAGS = -g -DDEBUG -Wall

# 'p' is for profiling information, ALSO needs pg in linker!
#CPPFLAGS = -Wall -O3 -march=pentium4 -pg
# linker options
#LNFLAGS = -pg

# intel compiler flag
#LNFLAGS = -ipo

# gcc version 3.2.3
#INCDIR = -I../include -I/usr/include 
#LIBS = -lboost_program_options -lboost_filesystem -lboost_regex -lxerces-c

# gcc version > 3.2.3
MYPATH = /home/anton/local
INCDIR = -I../include -I$(MYPATH)/include -I/usr/include 
LIBDIR = -L$(MYPATH)/lib -L$(MYPATH)/lib/xercesc
LIBS = -lboost_program_options-gcc -lboost_filesystem-gcc -lboost_regex-gcc -lxerces-c

# Source/object paths
vpath %.cc ../src ../test ../python
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../bin
LIBPATH = ../lib

# Project vars
PROJECT = fluke
LIBRARY = flu
ALL = distribution.o \
      main.o fluke.o config.o stream_manager.o model.o factory.o \
      population_reader.o agent_reader.o \
      observer_manager.o logger.o population.o environment.o \
      scaling.o selection.o mutate_rates.o\
      module_agent.o simple_agent.o agent.o \
      genome.o chromosome.o bsite.o repeat.o centromere.o \
      ordinary_dstream.o module_dstream.o transfac.o retroposon.o \
      shortseq.o observer.o subject.o
OBJECTS = $(ALL)


# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so

$(OBJECTS): %.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d *~

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak 

