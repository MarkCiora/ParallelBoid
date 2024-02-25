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
EXECVI := visualize

all: boidgen visualize

boidgen: $(EXECBG) build/boidgen/

$(EXECBG): $(OBJBG)
	mkdir -p build/boidgen/
	$(CXX) $(OBJBG) -o $(EXECBG)

build/boidgen/%.o : src/boidgen/%.cpp
	mkdir -p build/boidgen/
	$(CXX) $(CXXFLAGS) $(INC_FLAGSBG) -c $< -o $@


$(BUILD_DIR)/%.cpp.o: %.cpp

visualize:
	$(CXX) $(CXXFLAGS) $(INC_FLAGSVI) -c $(SRCVI) -o $(OBJVI)
	$(CXX) $(OBJVI) -o $(EXECVI)

clean:
	rm -rf build