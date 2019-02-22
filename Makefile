
CC = mpicc
CFLAGS = -Wall

objects = caqr.o main.o tsqr.o parallelsetup.o hhorth.o

all: catsqr

catsqr: $(objects)
	module load cports openmpi; \
	$(CC) $(CFLAGS) -o catsqr $(objects) -lm

main.o: main.c
	module load cports openmpi; \
	$(CC) $(CFLAGS) -c main.c

caqr.o: caqr.c caqr.h
	module load cports openmpi; \
	$(CC) $(CFLAGS) -c caqr.c

tsqr.o: tsqr.c tsqr.h
	module load cports openmpi; \
	$(CC) $(CFLAGS) -c tsqr.c

parallelsetup.o: parallelsetup.c parallelsetup.h
	module load cports openmpi; \
	$(CC) $(CFLAGS) -c parallelsetup.c

hhorth.o: hhorth.c hhorth.h
	module load cports openmpi; \
	$(CC) $(CFLAGS) -c hhorth.c

.PHONY: clean
clean:
	rm -f $(objects) catsqr

test: catsqr
	module load cports openmpi; \
	mpiexec -n 4 catsqr

