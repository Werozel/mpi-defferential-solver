
TARGET = main.cpp
OUTPUT = main.o

STDLIB_FLAG = -stdlib=libc++

build:
	g++ ${TARGET} -std=c++11 -o ${OUTPUT}

build-mpi:
	mpicxx ${TARGET} -o ${OUTPUT}

run-mpi-128:
	mpirun -n 1 ${OUTPUT} 1 0.00009 32

run-mpi-256:
	mpirun -n 8 ${OUTPUT} 1 0.00001 128

run:
	./main 12 12 12 4

clean:
	rm -f ${OUTPUT}
