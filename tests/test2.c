#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "Misc.h"

int main(int argc, char **argv) {

	int rank, size, i;
	int buf;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for(i = 0; i < 10; i++) {

		beginIter();

		if(i % 2 == 0) {
			if(rank < size / 2) {
				buf = rank;
				MPI_Send(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
			}
			else {
				MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			
		}
		else {
			if(rank >= 1 && rank <= size / 2) {
				MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			else {
				buf = rank;
				MPI_Send(&buf, 1, MPI_INT, (rank - 1 + size) % size, rank, MPI_COMM_WORLD);
			}
		
		}

		endIter(i, rank);
	}

	endFor(i);
	

	MPI_Finalize();

	return 0;
}