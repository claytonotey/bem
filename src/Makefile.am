SUBDIRS = util
AM_CXXFLAGS = @BEM_CFLAGS@ -fno-show-column

libbemincludedir = .
libbeminclude_HEADERS = bem.h
lib_LTLIBRARIES = libbem.la

libbem_la_SOURCES = ShapeParser.cpp FaceIndex.cpp Shape.cpp  ShapeMesh.cpp GreensFunction.cpp Options.cpp Impedance.cpp RWG.cpp Triangle.cpp AnalyticIntegrals.cpp RWGTree.cpp SurfaceParameterization.cpp SparseBLAS.cpp Source.cpp Progress.cpp Output.cpp VolumeMesh.cpp Tetgen.cpp TetgenPredicates.cpp Tetrahedron.cpp Flux.cpp FrequencySweep.cpp HeatTransfer.cpp

#Nonlocal.cpp HydrodynamicDielectric.cpp 

transfer_LDADD = libbem.la util/libutil.la -lblas -lm -llapack -lxerces-c
transfer_LDFLAGS = $(os_ldflags)
transfer_SOURCES = transfer.cpp

transfermpi_LDADD = libbem.la util/libutil.la -lblas -lm -llapack -lxerces-c
transfermpi_LDFLAGS = $(os_ldflags)
transfermpi_SOURCES = transfer-mpi.cpp

testparams_LDADD = libbem.la util/libutil.la -lblas -lm -llapack -lxerces-c
testparams_LDFLAGS = $(os_ldflags)
testparams_SOURCES = testparams.cpp

bin_PROGRAMS = @BEM_PROGS@ testparams
EXTRA_PROGRAMS= transfer transfermpi