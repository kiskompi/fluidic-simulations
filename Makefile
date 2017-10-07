CC=g++
CFLAGS=-std=c++14 -Ofast -pg -fopenmp -Wall -Wshadow -pedantic -Werror -Wsign-compare -Wtype-limits -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized

.c.o:
	$(CC) -c $(CFLAGS) $<

all:    navier

clean:
	rm -f *.o

navier: boundary.o init.o main.o output.o simulation.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

boundary.o       : datadef.h
init.o           : datadef.h
main.o           : boundary.h datadef.h init.h simulation.h
simulation.o     : datadef.h init.h

