CXX:=clang++
CXXFLAGS:= -Ofast -Wall -std=c++20 -fopenmp

SRC = src
OBJ = obj
DEP = dep
EXT = cpp

SOURCES := $(wildcard $(SRC)/*.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),$(OBJ)/%.o,$(SOURCES))
DEPENDS := $(patsubst $(SRC)/%.$(EXT),$(DEP)/%.d,$(SOURCES))

all: run

clena: clean

debug: CXXFLAGS += -g
debug: run

profile: CXXFLAGS += -pg
profile: run

run: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

-include $(DEPENDS)

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -MF $(DEP)/$*.d -c $< -o $@

clean:
	clear && rm $(OBJ)/*.o $(DEP)/*.d output/*.out run
