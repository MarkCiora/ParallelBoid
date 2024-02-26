CXX := g++

CXXFLAGSBG := -std=c++11 -Wall
CXXFLAGSVI := -std=c++11 -Wall -lglut -lGLU -lGL

SRCBG := $(wildcard ./src/boidgen/*.cpp)
SRCVI := $(wildcard ./src/visualize/*.cpp)

OBJBG := $(subst src,build,$(SRCBG:.cpp=.o))
OBJVI := $(subst src,build,$(SRCVI:.cpp=.o))

INCBG := ./include/boidgen/
INCVI := ./include/visualize/
INC_FLAGSBG := $(addprefix -I,$(INCBG))
INC_FLAGSVI := $(addprefix -I,$(INCVI))

EXECBG := build/boids
EXECVI := build/vis

all: boidgen visualize

boidgen: $(EXECBG) build/boidgen/

$(EXECBG): $(OBJBG)
	mkdir -p build/boidgen/
	$(CXX) $(OBJBG) -o $(EXECBG)

build/boidgen/%.o : src/boidgen/%.cpp
	mkdir -p build/boidgen/
	$(CXX) $(CXXFLAGSBG) $(INC_FLAGSBG) -c $< -o $@

visualize: $(EXECVI) build/visualize/

$(EXECVI): $(OBJVI)
	mkdir -p build/visualize/
	$(CXX) $(OBJVI) -o $(EXECVI)

build/visualize/%.o : src/visualize/%.cpp
	mkdir -p build/visualize/
	$(CXX) $(CXXFLAGSVI) $(INC_FLAGSVI) -c $< -o $@

clean:
	rm -rf build