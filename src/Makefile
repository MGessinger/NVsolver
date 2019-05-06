CC = gcc
CFLAGS = -std=c99 -Wall -D_GNU_SOURCE
.SUFFIXES: .o .c
.c.o: ; $(CC) -c $(CFLAGS) $*.c

OBJ = fields.o blas.o blatt2.o poisson.o

blatt0: $(OBJ)
	$(CC) $(CFLAGS) -o blatt2 $(OBJ) -lm

fields.o:	fields.h
blas.o:		blas.h
poisson.o:	poisson.h

clean:
	rm -rf *.o blatt2
 
