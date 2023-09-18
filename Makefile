CXX:=g++
CXXFLAGS:= -std=c++20

SRC = src
OBJ = obj
DEP = dep
EXT = cpp

SOURCES := $(wildcard $(SRC)/*.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),$(OBJ)/%.o,$(SOURCES))
DEPENDS := $(patsubst $(SRC)/%.$(EXT),$(DEP)/%.d,$(SOURCES))


all: CXXFLAGS+= -O1 -Wall -fopenmp
all: run

single: CXXFLAGS+= -Ofast -Wall -march=native
single: run

release: CXXFLAGS+= -Ofast -fopenmp -march=native
release: run

profile: CXXFLAGS+= -O1 -pg
profile: run

debug: CXXFLAGS+= -Og -g
debug: run

clena: clean

run: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

-include $(DEPENDS)

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -MF $(DEP)/$*.d -c $< -o $@

clean:
	-rm $(OBJ)/*.o $(DEP)/*.d output/*.out run
	clear