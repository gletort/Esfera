### Set ups
VERSION := 1.0.0

CC := g++ -Wfatal-errors

MODE := R

BUILD_DIR = build
SRC_DIR = src
BIN_DIR = bin
CUR_DIR = $(shell pwd)

#----------- Folders to look in
vpath %.o	$(BUILD_DIR)
vpath %.cpp		$(SRC_DIR)/core $(SRC_DIR)/tools
vpath %.h 	$(SRC_DIR)/core $(SRC_DIR)/tools


#---------- Compilation flags for different modes
WARN     := -Wall -Wno-unknown-pragmas # -Wno-unused
DEBUG    := -g3 -fno-rtti -ggdb -std=c++0x -march=native -fopenmp -DNDEBUG
RELEASE  := -fno-tree-vectorize -fno-rtti -ffast-math -funroll-loops -finline-functions -Wno-deprecated
PROFILE  := $(RELEASE) -O0 -fno-inline -pg
FlagsD := -O0 $(DEBUG) $(WARN) 
FlagsR := $(RELEASE) -march=native -O3 -s -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++0x
FlagsP := $(PROFILE)

COMPILE_COMMAND := $(CC) $(Flags$(MODE)) 


#------------ What to compile

CORE_OBJECTS := simul.o sphere.o sphere_set.o sphere_prop.o space.o dna.o
TOOLS_OBJECTS := random.o vector3d.o readXML.o tinyxml2.o 

ALL_OBJECTS := $(TOOLS_OBJECTS) $(CORE_OBJECTS)


MAIN := ./$(SRC_DIR)/main/Esfera.cpp

all: make_all esfera

make_all: $(ALL_OBJECTS)

build:
	@mkdir -p $@
bin:
	@mkdir -p $@

$(ALL_OBJECTS): %.o: %.cpp %.h | build 
	$(COMPILE_COMMAND) $(CXXFLAGS) $(INC) -c $< -o $(BUILD_DIR)/$(@F)

esfera: $(BIN_DIR) $(MAIN)
	$(COMPILE_COMMAND) $(INC) -o $(BIN_DIR)/Esfera ./src/main/Esfera.cpp $(BUILD_DIR)/*.o


#---------------- Generate doc with Doxygen
.PHONY: doc clean cleanbin mrproper

doc:
	doxygen Doxyfile

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BIN_DIR)/* 

cleanbin:
	rm -f $(BIN_DIR)/*

mrproper: clean cleanbin

