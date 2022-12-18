
TARGET = main.cpp
OUTPUT = main.o

STDLIB_FLAG = -stdlib=libc++

build:
	g++ ${TARGET} -std=c++11 -o ${OUTPUT}

build-mpi:
	mpicxx ${TARGET} -o ${OUTPUT}

run-mpi:
	mpirun -n 2 ${OUTPUT} 1 128

run:
	./main 12 12 12 4

clean:
	rm -f ${OUTPUT}
