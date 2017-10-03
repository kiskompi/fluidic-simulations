CC=gcc
CFLAGS=-Ofast -pg -std=c99 -fopenmp

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

