#------------------------------------------------------------------------
# Compiler and flags
CXX = c++
CXXFLAGS = -std=c++14 -stdlib=libc++ -m64 -mmacosx-version-min=13.0 $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)

# Executable names
BINARIES = Reconstruction HTArrayTraining

# Default target: build all
all: $(BINARIES)

# Rule to build each executable from its matching .cpp file
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Clean rule
clean:
	rm -f $(BINARIES)

#------------------------------------------------------------------------
