CC=mpicc
CFLAGS=-lm

DEPS = utility.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simulation: main.c utility.c initialconfigurations.c utility.h
	$(CC) -o simulation main.c utility.c initialconfigurations.c $(CFLAGS)

sim: simulation.o
	$(CC) simulation simulation.o $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o 
