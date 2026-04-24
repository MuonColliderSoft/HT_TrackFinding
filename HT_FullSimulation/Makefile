#------------------------------------------------------------------------
# Compiler and flags
CXX = c++
CXXFLAGS = -m64 $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)

# macOS-specific flags
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CXXFLAGS += -stdlib=libc++ -mmacosx-version-min=13.0
endif

# Executable names
BINARIES = Reconstruction HTArrayTraining CutOptimization

# Default target: build all
all: $(BINARIES)

# Rule to build each executable from its matching .cpp file
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Clean rule
clean:
	rm -f $(BINARIES)

#------------------------------------------------------------------------
