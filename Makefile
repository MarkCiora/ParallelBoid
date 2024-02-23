CXX := g++

CXXFLAGS := -std=c++11 -Wall

SRCBG := $(wildcard ./src/boidgen/*.cpp)
SRCVI := $(wildcard ./src/visualize/*.cpp)

OBJBG := $(subst src,build,$(SRCBG:.cpp=.o))
OBJVI := $(subst src,build,$(SRCVI:.cpp=.o))

INCBG := /include/boidgen/
INCVI := /include/visualize/
INC_FLAGSBG := $(addprefix -I,$(INCBG))
INC_FLAGSVI := $(addprefix -I,$(INCVI))

EXECBG := build/boids
EXECVI := visualize

test:
	echo $(SRCBG)
	echo $(OBJBG)

all: boidgen visualize

boidgen:
	mkdir -p build/boidgen/
	$(CXX) $(CXXFLAGS) $(INC_FLAGSBG) -c $(SRCBG) -o $(OBJBG)
	$(CXX) $(OBJBG) -o $(EXECBG)

visualize:
	$(CXX) $(CXXFLAGS) $(INC_FLAGSVI) -c $(SRCVI) -o $(OBJVI)
	$(CXX) $(OBJVI) -o $(EXECVI)

clean:
	rm -rf build
	echo "nothing3"