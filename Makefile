KOKKOS_PATH :=/home/akomporday/git-project/include
KOKKOS_DEVICES := "OpenMP, Cuda"

MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRC_DIR := $(dir $(MAKEFILE_PATH))

SRC = main.cpp Simulation/simulation.cpp 
OBJ = $(SRC:$(SRC_DIR)%.cpp=%.o) $(SRC:$(SRC_DIR)Simulation/%.cpp=%.o)

#SRC = $(wildcard *.cpp)
#OBJ = $(SRC:%.cpp=%.o)

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
  EXE = $(addsuffix .cuda, $(shell basename $(SRC_DIR)))
else
  CXX = g++
  EXE = $(addsuffix .host, $(shell basename $(SRC_DIR)))
endif

CXXFLAGS = -std=c++11 -fopenmp  -g  -O0 -I$(SRC_DIR)include 
LINK ?= $(CXX) -L/home/akomporday/git-project/lib/ $(CXXFLAGS)
LDFLAGS ?=

include $(KOKKOS_PATH)/Makefile.kokkos

DEPFLAGS = -M

LIB = 


build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	echo $(SRC)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: 
	rm -f *.a *.o *.cuda *.host

# Compilation rules

%.o:$(SRC_DIR)/%.cpp $(KOKKOS_DEPENDS) 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
%.o:$(SRC_DIR)Simulation/%.cpp $(KOKKOS_DEPENDS) 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

