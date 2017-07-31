/* Not meet Rule 1, 2, 3 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

// #include "Misc.h"

int main(int argc, char **argv) {

	int size, rank, i, buf;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for(i = 0; i < 3; i++) {
		// beginIter();
		if(rank != 0) {
			buf = rank;
			MPI_Send(&buf, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
		}
		else {
			MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		// endIter(i,rank);
	}
	// endFor();

	MPI_Finalize();

	return 0;
}
