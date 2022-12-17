
TARGET = main.cpp
OUTPUT = main.o

STDLIB_FLAG = -stdlib=libc++

build:
	g++ main.cpp -std=c++11 -o main

build-mpi:
	mpicxx ${TARGET} -o ${OUTPUT}

run-mpi:
	mpirun -n 4 ${OUTPUT} 1 1 1 2

run:
	./main 1 1 1 2

clean:
	rm -f ${OUTPUT}
