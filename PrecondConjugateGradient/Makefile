CC = gcc
MPCC = mpicc
CFLAGS = -Wall

objects = main.o conjugateg.o haloswaps.o parallelsetup.o problemsetup.o vectorops.o precondition.o

cg: $(objects)
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -o cg $(objects) -lm

main.o: main.c structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

parallelsetup.o: parallelsetup.c parallelsetup.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

problemsetup.o: problemsetup.c problemsetup.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

precondition.o: precondition.c precondition.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

haloswaps.o: haloswaps.c haloswaps.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

vectorops.o: vectorops.c vectorops.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

conjugateg.o: conjugateg.c conjugateg.h structs.h
	module load cports openmpi;\
	$(MPCC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(objects) cg


test: cg
	module load cports openmpi;\
	mpiexec -n 8 ./cg

