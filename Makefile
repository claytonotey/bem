default: transfer scatter

CC=mpicxx
CFLAGS = -Wall -ffast-math -Wno-deprecated -fno-show-column -g -O3 -fopenmp -std=c++11

clean:
	rm *.o

.cpp.o:
	g++ $(CFLAGS) -c $<

Impedance.o: TriangleQuadrature.h Matrix.h DenseMatrix.h Rk.h HMatrix.h MatrixOperations.h Vector3.h ComplexVector3.h
SurfaceParameterization.o: ConjugateGradient.h Array.h SparseBLAS.h
SparseBLAS.o: SparseBLAS.h
Flux.o: TetrahedronQuadrature.h
ShapeMesh.o: Array.h

BEM_O = ShapeParser.o FaceIndex.o Shape.o MathUtils.o ShapeMesh.o GreensFunction.o Options.o Impedance.o RWG.o Triangle.o AnalyticIntegrals.o  Transform.o Rotation.o ZRotation.o DenseMatrix.o RWGTree.o SurfaceParameterization.o SparseBLAS.o Source.o Progress.o Output.o ReferenceCountingPointer.o VolumeMesh.o Tetrahedron.o Flux.o SI.o Silica.o Rk.o MatrixOperations.o FrequencySweep.o SiliconCarbide.o Nonlocal.o HydrodynamicDielectric.o

transfer: $(BEM_O) transfer.o
	g++ -g -O3 -ffast-math -o transfer transfer.o $(BEM_O)  -L./lapack -lxerces-c -L./tetgen -lblas -llapack -ltetgen -lgomp -lcurl

scatter: $(BEM_O) scatter.o
	g++ -g -O3 -ffast-math -o scatter scatter.o $(BEM_O)  -L./lapack -lxerces-c -L./tetgen -lblas -llapack -ltetgen -lgomp -lcurl

