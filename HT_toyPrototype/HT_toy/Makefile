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

HDRS = Parameters.h Statistics.h Hit.h Track.h
INCS = BibFileReader.cpp Event.cpp DetectorGeometry.cpp TrackGeometry.cpp\
       GlobalConstants.cpp Hit.cpp HTArray.cpp Statistics.cpp Track.cpp

#------------------------------------------------------------------------

all: ConvertBIBHits HTArrayTraining EventGeneration

ConvertBIBHits: DetectorGeometry.cpp ConvertBIBHits.cpp
	@echo "Compiling $@ ... "
	$(CXX) $@.cpp $(CXXFLAGS) $(ROOTLIBS) $(ROOTGLIBS) -o $@

HTArrayTraining: $(HDRS) $(INCS) HTArrayTraining.cpp
	@echo "Compiling $@ ... "
	$(CXX) $@.cpp $(CXXFLAGS) $(ROOTLIBS) $(ROOTGLIBS) -o $@

EventGeneration: $(HDRS) $(INCS) EventGeneration.cpp
	@echo "Compiling $@ ... "
	$(CXX) $@.cpp $(CXXFLAGS) $(ROOTLIBS) $(ROOTGLIBS) -o $@

clean:
	@rm -f ConvertBIBHits HTArrayTraining EventGeneration

#------------------------------------------------------------------------
