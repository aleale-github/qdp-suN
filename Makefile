BIN=$$HOME/bin
CONFIG_D7=$$HOME/chroma_large_n/nc7/install/qdp++/bin/qdp++-config

ifdef REGENCLUSTER
BIN=/kerndata3/ama62246/bin
CONFIG_D7=/kerndata3/ama62246/chroma_large_n/nc7/install/qdp++/bin/qdp++-config
# CONFIG=/kerndata3/ama62246/qdp_large_n_parallel/nc7/install/qdp++/bin/qdp++-config
endif

ifdef GPUSWANSEA
CONFIG_D7=/home/ale/qdp_largeN/double-prec/nc7/install/qdp-jit/bin/qdp++-config
CONFIG_D5=/home/ale/qdp_largeN/double-prec/nc5/install/qdp-jit/bin/qdp++-config
CONFIG_S7=/home/ale/qdp_largeN/single-prec/nc7/install/qdp-jit/bin/qdp++-config
CONFIG_S5=/home/ale/qdp_largeN/single-prec/nc5/install/qdp-jit/bin/qdp++-config
endif

CXX_D5=$(shell $(CONFIG_D5) --cxx) 
CXX_S5=$(shell $(CONFIG_S5) --cxx)
CXX_D7=$(shell $(CONFIG_D7) --cxx) 
CXX_S7=$(shell $(CONFIG_S7) --cxx)

CXXFLAGS_D7=$(shell $(CONFIG_D7) --cxxflags)
CXXFLAGS_D5=$(shell $(CONFIG_D5) --cxxflags) 
CXXFLAGS_S7=$(shell $(CONFIG_S7) --cxxflags) 
CXXFLAGS_S5=$(shell $(CONFIG_S5) --cxxflags) 

LDFLAGS_D7=$(shell $(CONFIG_D7) --ldflags)
LDFLAGS_D5=$(shell $(CONFIG_D5) --ldflags)
LDFLAGS_S7=$(shell $(CONFIG_S7) --ldflags)
LDFLAGS_S5=$(shell $(CONFIG_S5) --ldflags)

LIBS_D7=$(shell $(CONFIG_D7) --libs)
LIBS_D5=$(shell $(CONFIG_D5) --libs)
LIBS_S7=$(shell $(CONFIG_S7) --libs)
LIBS_S5=$(shell $(CONFIG_S5) --libs)


SRCDIR=src
SOURCES=$(SRCDIR)/reunit.cpp $(SRCDIR)/main.cpp $(SRCDIR)/update.cpp $(SRCDIR)/utils.cpp $(SRCDIR)/meas.cpp 

clean :
	rm -f $(BIN)/ale_suN*

touch :
	find | xargs -n1 touch

all: clean touch D7 D5 S7 S5

D5 : $(SOURCES)
	$(CXX_D5) $(CXXFLAGS_D5) $^ -o $(BIN)/ale_suN-D5 $(LDFLAGS_D5) $(LIBS_D5)

D7 : $(SOURCES)
	$(CXX_D7) $(CXXFLAGS_D7) $^ -o $(BIN)/ale_suN-D7 $(LDFLAGS_D7) $(LIBS_D7)

S5 : $(SOURCES)
	$(CXX_S5) $(CXXFLAGS_S5) $^ -o $(BIN)/ale_suN-S5 $(LDFLAGS_S5) $(LIBS_S5)

S7 : $(SOURCES)
	$(CXX_S7) $(CXXFLAGS_S7) $^ -o $(BIN)/ale_suN-S7 $(LDFLAGS_S7) $(LIBS_S7)
