
# configure the Makefile here
HOSTNAME = $(shell hostname)

CXX = mpiCC
CXXFLAGS += -Wall -pthread
LDFLAGS += -pthread
LDFLAGS += -lboost_mpi -lboost_serialization

ifeq ($(PROFILING),1)
    CXXFLAGS += -pg -g
endif

ifeq ($(DEBUGGING),1)
    CXXFLAGS += -g -O0
else
    CXXFLAGS += -DNDEBUG -O3 -funroll-loops 
endif

# enable this if you want things to be fast
# SSE2 assembly version of innermost loop
#CXXFLAGS+=-DHAVE_SSE2
# BLAS library (IMKL, ACML, ...)
#CXXFLAGS+=-DHAVE_BLAS

# linux setaffinity call for automatic pinning
#CXXFLAGS+=-DHAVE_LINUX_SETAFFINITY_NP

UTILITIES = writevtk.o mpisupp.o ublas.o TinyCmd.o shared.o stencils_iso.o tinyconf.o bzopen.o eio.o
BINARY = Erpel.bin
WRAPPER = Erpel

all: $(BINARY) $(WRAPPER) volume_fractions

$(BINARY): $(UTILITIES) erpel.o driver.o
	$(CXX) -o $@ erpel.o driver.o $(UTILITIES) $(LDFLAGS)

$(WRAPPER): wrapper.bash
	cp wrapper.bash $@
	chmod +x $@

volume_fractions: volume_fractions.o $(UTILITIES)
	$(CXX) -o $@ $< $(UTILITIES) $(LDFLAGS)

coarsen: coarsen.o $(UTILITIES)
	$(CXX) -o $@ $< $(UTILITIES) $(LDFLAGS)

coarsen.o: coarsen.cpp shared.hpp config.hpp

volume_fractions.o: volume_fractions.cpp shared.hpp config.hpp

driver.o: driver.cpp erpel.hpp mpisupp.hpp TinyCmd.hpp writevtk.hpp tinyconf.h config.hpp

erpel.o: erpel.cpp revision.hpp shared.hpp sse2_dgemv.hpp perf.hpp mpisupp.hpp erpel.hpp config.hpp stencils_iso.hpp

eio.o: eio.cpp erpel.hpp shared.hpp config.hpp

mpisupp.o: mpisupp.cpp mpisupp.hpp config.hpp

writevtk.o: writevtk.cpp writevtk.hpp config.hpp

ublas.o: ublas.cpp ublas.hpp config.hpp

TinyCmd.o: TinyCmd.cpp TinyCmd.hpp

shared.o: shared.cpp shared.hpp config.hpp

tinyconf.o: tinyconf.cpp tinyconf.h config.hpp

bzopen.o: bzopen.cpp shared.hpp config.hpp

clean:
	rm -f $(BINARY) volume_fractions *.o $(WRAPPER)

.PHONY: all clean
