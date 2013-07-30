
# configure the Makefile here
HOSTNAME = $(shell hostname)

CXX = mpiCC
CXXFLAGS += -Wall -pthread
LDFLAGS += -pthread

ifeq ($(PROFILING),1)
    CXXFLAGS += -pg -g
endif

ifeq ($(DEBUGGING),1)
    CXXFLAGS += -g -O0
else
    CXXFLAGS += -DNDEBUG -O3 -funroll-loops 
endif

ifeq ($(HOSTNAME),noether)
    $(warning configure for Theo1 machines, compilation on noether; be sure to load the following modules:)
    $(warning    boost)
    CXXFLAGS += -DHAVE_SSE2
    CXXFLAGS += -pedantic -Wextra
    LDFLAGS += -lboost_mpi
    # BLAS is broken for now, since Ubuntu ships without _gfortran_transfer_real_write@GFORTRAN_1.4
else

    # Config for RRZE clusters.
    RRZECLUSTER = woody
    ifeq ($(HOSTNAME),lima1)
        RRZECLUSTER = lima
    endif
    ifeq ($(HOSTNAME),lima2)
        RRZECLUSTER = lima
    endif

ifeq ($(RRZECLUSTER),lima)
    $(warning configure for lima cluster; be sure to load the following modules:)
    $(warning    boost/1.44.0-intel11.1-impi4.0 intel64/11.1.073)
    # we build with ICC on lima, so no Wextra.  warnings are disabled in config.hpp.
    CXXFLAGS += -I$(BOOST_INCDIR)
    LDFLAGS += -L$(BOOST_LIBDIR) -lboost_mpi -lboost_serialization
    # enable fast dgemv24
    CXXFLAGS += -DHAVE_SSE2
    # don't use BLAS for now
else
    $(warning configure for woody/tinyblue cluster; be sure to load the following modules:)
    $(warning   mpich/p4-intel64 intel64/11.1.064 acml-intel64/4.1.0_mp)
    # we build with ICC on Woody, so no Wextra.  warnings are disabled in config.hpp.
    CXXFLAGS += -I/home/rrze/mpt1/mpt133/local/boost/1_39_0-gcc-mpich-newbinutils/include/boost-1_39
    LDFLAGS += /home/rrze/mpt1/mpt133/local/boost/1_39_0-gcc-mpich-newbinutils/lib/libboost_mpi-gcc42-mt-1_39.a
    LDFLAGS += /home/rrze/mpt1/mpt133/local/boost/1_39_0-gcc-mpich-newbinutils/lib/libboost_serialization-gcc42-mt-1_39.a
    # enable fast dgemv24
    CXXFLAGS += -DHAVE_SSE2
    # use AMD BLAS
    LDFLAGS += $(ACML_LIB)
    CXXFLAGS += -DHAVE_BLAS
endif
endif

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
