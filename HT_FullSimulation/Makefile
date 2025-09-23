#------------------------------------------------------------------------
# Compiler and flags
CXX = c++
CXXFLAGS = $(shell root-config --cflags)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CXXFLAGS += -stdlib=libc++ -mmacosx-version-min=13.0
endif
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
