#MPI_BIN = /usr/local/bin/
MPI_BIN = /usr/lib64/openmpi/1.4-gcc/bin/
CC      = $(MPI_BIN)/mpicc
CFLAGS  = -DMPE_GRAPHICS -c -Wall
LFLAGS  = -L/usr/X11R6/lib -lmpe -lX11 -lm

p2-mandel: mandel.o display.o
	$(CC) -o p2-mandel mandel.o display.o $(LFLAGS)

mandel.o: mandel.c display.h
	$(CC) $(CFLAGS) mandel.c 

display.o: display.c display.h
	$(CC) $(CFLAGS) display.c

run1:
	$(MPI_BIN)mpirun -np 1 a.out

clean:
	rm -f a.out *.o *~ *# PI*

