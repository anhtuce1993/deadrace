/* Meet Rule 1, 2, not meet Rule 3 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#include "Misc.h"

int main(int argc, char **argv) {

	int size, rank, i, buf;
	MPI_Status status;
	double starttime, endtime; 

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	starttime = MPI_Wtime();
	beginFor();
	for(i = 0; i < 10; i++) {
		beginIter();
		if(rank < size / 2) {
			if(i % 2 == 0) {
				buf = rank;
				MPI_Send(&buf, 1, MPI_INT, rank + size / 2, rank, MPI_COMM_WORLD);
				MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			else {
				buf = rank;
				MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
				MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}
		else {
			if(i % 2 == 0) {
				buf = rank;
				MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
				MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			else {
				buf = rank;
                MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
				MPI_Recv(&buf, 1, MPI_INT, rank - size / 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}	
		endIter(i, rank);
	}
	
	endtime = MPI_Wtime();
	endFor(i, 0, rank);


	starttime = MPI_Wtime() - endtime + starttime;
        beginFor();
        for(i = 0; i < 10; i++) {
                beginIter();
                if(rank < size / 2) {
                        if(i % 2 == 0) {
                                buf = rank;
                                MPI_Send(&buf, 1, MPI_INT, rank + size / 2, rank, MPI_COMM_WORLD);
                                MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        }
                        else {
                                buf = rank;
                                MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
                                MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        }
                }
                else {
                        if(i % 2 == 0) {
                                buf = rank;
                                MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
                                MPI_Recv(&buf, 1, MPI_INT, (size - 1) - rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        }
                        else {
                                buf = rank;
                MPI_Send(&buf, 1, MPI_INT, (size - 1) - rank, rank, MPI_COMM_WORLD);
                                MPI_Recv(&buf, 1, MPI_INT, rank - size / 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        }
                }
                endIter(i, rank);
        }

        endtime = MPI_Wtime();

	endFor(i, endtime - starttime, rank);
	MPI_Finalize();

	return 0;
}
