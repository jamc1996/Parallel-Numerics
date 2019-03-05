
CC = gcc
CFLAGS = -Wall
objects = main.o matSetUp.o qrfactoring.o matops.o

all: gramS

gramS: $(objects)
	$(CC) $(CFLAGS) -o gramS $(objects) -lm

main.o: main.c dense.h
	$(CC) $(CFLAGS) -c $<

qrfactoring.o: qrfactoring.c qrfactoring.h dense.h
	$(CC) $(CFLAGS) -c $<

matSetUp.o: matSetUp.c matSetUp.h dense.h
	$(CC) $(CFLAGS) -c $<

matops.o: matops.c matops.h dense.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f  $(objects) *.txt gramS
