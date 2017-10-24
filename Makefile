KOKKOS_PATH = ../../..
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = $(wildcard ${KOKKOS_SRC_PATH}/home/akomporday/git-projects/*.cpp ${KOKKOS_SRC_PATH}/home/akomporday/git-projects/Simulation/*.cpp)
vpath %.cpp $(sort $(dir $(SRC)))

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS = 
EXE = 01_hello_world_lambda.cuda
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS += "enable_lambda"
else
CXX = g++
CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =  
EXE = maked_runnable 
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB"
endif

DEPFLAGS = -M

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =


$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)
