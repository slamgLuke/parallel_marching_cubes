# Makefile for Parallel Marching Cubes implementations

# Compiler and flags
CXX = g++
# CXXFLAGS = -std=c++11 -O3 -Wall -Wextra
CXXFLAGS = -std=c++17
OMPFLAGS = -fopenmp

# Target executables
TARGETS = get_omp_threads sequential_for parallel_mc_for parallel_mc_non_adapt

# Default target - build all executables
all: $(TARGETS)

# Get OpenMP threads count
get_omp_threads: get_omp_threads.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

# Sequential implementation
sequential_for: sequential_for.cpp tables.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Parallel implementation with iterative approach
parallel_mc_for: parallel_mc_for.cpp tables.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

# Parallel implementation with adaptive approach
parallel_mc_non_adapt: parallel_mc_non_adapt.cpp tables.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

# Clean target
clean:
	rm -f $(TARGETS)

# Debug versions with debugging symbols
debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGETS)

# Help target
help:
	@echo "Available targets:"
	@echo "  all                 - Build all executables (default)"
	@echo "  sequential_for      - Build sequential implementation"
	@echo "  parallel_mc_for     - Build parallel iterative implementation"
	@echo "  parallel_mc_non_adapt - Build parallel adaptive implementation"
	@echo "  debug               - Build all with debug symbols"
	@echo "  clean               - Remove all executables"
	@echo "  help                - Show this help message"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                - Build all programs"
	@echo "  make sequential_for - Build only sequential version"
	@echo "  make clean          - Clean all executables"
	@echo "  make debug          - Build with debug info"

# Declare phony targets
.PHONY: all clean debug help
