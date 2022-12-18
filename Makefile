
TARGET_MPI = main.cpp
TARGET_MPI_OMP = main_omp.cpp
OUTPUT_MPI = main-mpi.o
OUTPUT_MPI_OMP = main-mpi-omp.o

STDLIB_FLAG = -stdlib=libc++
C11_FLAG = -std=c++11

build-mpi:
	mpicxx ${TARGET_MPI} ${C11_FLAG} -o ${OUTPUT_MPI}

build-mpi-omp:
	mpicxx ${TARGET_MPI_OMP} ${C11_FLAG} -qthreaded -qsmp=omp -o ${OUTPUT_MPI_OMP}

run-mpi-local-128:
	mpirun -n ${n_proc} ${OUTPUT_MPI} 1 0.00001 128

run-mpi-local-256:
	mpirun -n ${n_proc} ${OUTPUT_MPI} 1 0.00001 256

run-mpi-local-512:
	mpirun -n ${n_proc} ${OUTPUT_MPI} 1 0.00001 512

run-mpi-128:
	mpisubmit.pl -p ${n_proc} ${OUTPUT_MPI} 1 0.00001 128

run-mpi-256:
	mpisubmit.pl -p ${n_proc} ${OUTPUT_MPI} 1 0.00001 256

run-mpi-512:
	mpisubmit.pl -p ${n_proc} ${OUTPUT_MPI} 1 0.00001 512

run-mpi-omp-128:
	mpisubmit.pl -p ${n_proc} -t 4 ${OUTPUT_MPI_OMP} 1 0.00001 128

run-mpi-omp-256:
	mpisubmit.pl -p ${n_proc} -t 4 ${OUTPUT_MPI_OMP} 1 0.00001 256

run-mpi-omp-512:
	mpisubmit.pl -p ${n_proc} -t 4 ${OUTPUT_MPI_OMP} 1 0.00001 512

clean:
	rm -f *.o
