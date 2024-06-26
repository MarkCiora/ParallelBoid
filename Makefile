CXX := g++
CXXCUDA := nvcc

CXXFLAGSBG := -std=c++11 -Wall
CXXFLAGSOMP := -std=c++11 -Wall -fopenmp
CXXFLAGSCUDA := 
CXXFLAGSVI := -std=c++11 -Wall -lglut -lGLU -lGL

SRCBG := $(wildcard ./src/boidgen/*.cpp)
SRCOMP := $(wildcard ./src/omp/*.cpp)
SRCCUDA := $(wildcard ./src/cuda/*.cpp)
SRCVI := $(wildcard ./src/visualize/*.cpp)

OBJBG := $(subst src,build,$(SRCBG:.cpp=.o))
OBJOMP := $(subst src,build,$(SRCOMP:.cpp=.o))
OBJCUDA := $(subst src,build,$(SRCCUDA:.cpp=.o))
OBJVI := $(subst src,build,$(SRCVI:.cpp=.o))

INCBG := ./include/boidgen/
INCOMP := ./include/omp/
INCCUDA := ./include/cuda/
INCVI := ./include/visualize/
INC_FLAGSBG := $(addprefix -I,$(INCBG))
INC_FLAGSOMP := $(addprefix -I,$(INCOMP))
INC_FLAGSCUDA := $(addprefix -I,$(INCCUDA))
INC_FLAGSVI := $(addprefix -I,$(INCVI))

EXECBG := build/boids
EXECOMP := build/boids_omp
EXECCUDA := build/boids_cuda
EXECVI := build/vis

all: boidgen

boidgen: $(EXECBG) build/boidgen/
omp: $(EXECOMP) build/omp/
cuda: $(EXECCUDA) build/cuda/

$(EXECBG): $(OBJBG)
	mkdir -p build/boidgen/
	$(CXX) $(OBJBG) -o $(EXECBG)

build/boidgen/%.o : src/boidgen/%.cpp
	mkdir -p build/boidgen/
	$(CXX) $(CXXFLAGSBG) $(INC_FLAGSBG) -c $< -o $@
	
$(EXECOMP): $(OBJOMP)
	mkdir -p build/omp/
	$(CXX) $(OBJOMP) -o $(EXECOMP)

build/omp/%.o : src/omp/%.cpp
	mkdir -p build/omp/
	$(CXX) $(CXXFLAGSOMP) $(INC_FLAGSOMP) -c $< -o $@
	
$(EXECCUDA): $(OBJCUDA)
	mkdir -p build/cuda/
	$(CXXCUDA) $(OBJCUDA) -o $(EXECCUDA)

build/cuda/%.o : src/cuda/%.cpp
	mkdir -p build/cuda/
	$(CXXCUDA) $(CXXFLAGSCUDA) $(INC_FLAGSCUDA) -c $< -o $@

visualize: $(EXECVI) build/visualize/

$(EXECVI): $(OBJVI)
	mkdir -p build/visualize/
	$(CXX) $(OBJVI) -o $(EXECVI)

build/visualize/%.o : src/visualize/%.cpp
	mkdir -p build/visualize/
	$(CXX) $(CXXFLAGSVI) $(INC_FLAGSVI) -c $< -o $@

clean:
	rm -rf build