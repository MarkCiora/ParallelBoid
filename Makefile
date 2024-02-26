CXX := g++

CXXFLAGS := -std=c++11 -Wall

SRCBG := $(wildcard ./src/boidgen/*.cpp)
SRCVI := $(wildcard ./src/visualize/*.cpp)

OBJBG := $(subst src,build,$(SRCBG:.cpp=.o))
OBJVI := $(subst src,build,$(SRCVI:.cpp=.o))

INCBG := ./include/boidgen/
INCVI := ./include/visualize/
INC_FLAGSBG := $(addprefix -I,$(INCBG))
INC_FLAGSVI := $(addprefix -I,$(INCVI))

EXECBG := build/boids
EXECVI := build/visualize

all: boidgen visualize

boidgen: $(EXECBG) build/boidgen/

$(EXECBG): $(OBJBG)
	mkdir -p build/boidgen/
	$(CXX) $(OBJBG) -o $(EXECBG)

build/boidgen/%.o : src/boidgen/%.cpp
	mkdir -p build/boidgen/
	$(CXX) $(CXXFLAGS) $(INC_FLAGSBG) -c $< -o $@

visualize: $(EXECVI) build/visualize/

$(EXECVI): $(OBJVI)
	mkdir -p build/visualize/
	$(CXX) $(OBJVI) -o $(EXECVI)

build/visualize/%.o : src/boidgen/%.cpp
	mkdir -p build/visualize/
	$(CXX) $(CXXFLAGS) $(INC_FLAGSVI) -c $< -o $@

clean:
	rm -rf build