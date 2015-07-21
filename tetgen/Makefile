default: libtetgen.a

CFLAGS = -Wall -g -O3

clean:
	rm *.o

.cpp.o:
	g++ $(CFLAGS) -c $<

TETGEN_O = TetgenPredicates.o Tetgen.o
libtetgen.a: $(TETGEN_O)
	ar rcs libtetgen.a $(TETGEN_O)
	ranlib libtetgen.a
