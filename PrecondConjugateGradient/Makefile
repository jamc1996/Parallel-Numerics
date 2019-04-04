CC = gcc
MPCC = mpicc
CFLAGS = -Wall

objects = main.o conjugateg.o parallelsetup.o problemsetup.o vectorops.o precondition.o

cg: $(objects)
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -o cg $(objects)

main.o: main.c
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

parallelsetup.o: parallelsetup.c parallelsetup.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

problemsetup.o: problemsetup.c problemsetup.h
	$(CC) $(CFLAGS) -c $<

precondition.o: precondition.c precondition.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

vectorops.o: vectorops.c vectorops.h
	$(CC) $(CFLAGS) -c $<

conjugateg.o: conjugateg.c conjugateg.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(objects) cg


test: cg
	module load cports openmpi;\
	mpiexec -n 8 ./cg

