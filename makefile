# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -Isrc -Iinclude -O3 -fopenmp -g

# Source files directory
SRC_DIR = src

# Default source
DEFAULT_SRC = $(SRC_DIR)/main.cpp

# Default executable name
EXECUTABLE = polymorph

# Default target
all: $(EXECUTABLE)

# Rule for building the default executable
$(EXECUTABLE): $(DEFAULT_SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

# Rule to compile any cpp file to the default executable name
%: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $(EXECUTABLE)

# Rule to clean up generated files
clean:
	rm -f $(EXECUTABLE)
