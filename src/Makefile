#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test
#INC1=${PETSC_DIR}/include
#INC2=${PETSC_DIR}/$(PETSC_ARCH)/include
CC = mpicxx -g
CFLAGS = -c -std=c++11 -fpic -O3 -I$(INC1) -I$(INC2)
LDFLAGS = -shared
SOURCES=BoundaryConditions/BoundaryConditions.cpp\
	Cell/Cell.cpp\
        Face/Face.cpp\
        Grid/Grid.cpp\
        Grid/ReadGrid.cpp\
        Grid/SetGrid.cpp\
        Grid/Output.cpp\
        InterFlux/Interflux.cpp\
        InterFlux/Roe/Roe.cpp\
        Matrix5/Operators.cpp\
        Output/Output.cpp\
        Solver/GaussSeidel.cpp\
        Solver/Solver.cpp\
        Solver/Explicit.cpp\
        Solver/Implicit.cpp\
        Solver/PreSolver.cpp\
        Solver/Output.cpp\
        Solver/Residual.cpp\
        Vector/Vector.cpp\
        Time/Time.cpp\
        Method/Method.cpp\
        IBlank/IBlank.cpp\
        LinearAlgebra/LinearAlgebra.cpp\
        Gradient/Gradient.cpp\
        Limiter/Limiter.cpp\
	Point/Point.cpp\
	MovingMesh/MovingMesh.cpp\
	Coeffs/Coeffs.cpp
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=codeBaku
TARGET=libcb.so

all: $(TARGET)

#all: $(SOURCES) $(EXECUTABLE)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(PETSC_LIB) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
