#------------------------------------------------------------------------

#CXX = g++
#CXXFLAGS = -O3 -Wall -fPIC -g

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)

#------------------------------------------------------------------------

HDRS = Statistics.h Hit.h Track.h
INCS = BibFileReader.cpp Event.cpp Geometry.cpp\
       Hit.cpp HTArray.cpp Statistics.cpp Track.cpp\
       HTArrayTraining.cpp EventGeneration.cpp

#------------------------------------------------------------------------

all: HTArrayTraining EventGeneration

HTArrayTraining EventGeneration: $(HDRS) $(INCS)
	@echo "Compiling $@ ... "
	$(CXX) $@.cpp $(CXXFLAGS) $(ROOTLIBS) $(ROOTGLIBS) -o $@

clean:
	@rm -f HTArrayTraining EventGeneration
 

#------------------------------------------------------------------------
