#!/bin/bash


# shell script to to test parallel qr methods.

module load cports openmpi
make
mpiexec -n 4 catsqr -b 10 > results.txt 2> error.txt
mpiexec -n 4 catsqr -b 50 >> results.txt 2> error.txt
mpiexec -n 4 catsqr -b 100 >> results.txt 2> error.txt
mpiexec -n 4 catsqr -b 500  >> results.txt 2> error.txt

mpiexec -n 8 catsqr -b 10  >> results.txt 2> error.txt
mpiexec -n 8 catsqr -b 50  >> results.txt 2> error.txt
mpiexec -n 8 catsqr -b 100 >> results.txt 2> error.txt
mpiexec -n 8 catsqr -b 500  >> results.txt 2> error.txt

