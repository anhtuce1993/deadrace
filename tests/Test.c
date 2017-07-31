#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"


/******************************************************* mpipi.h **************************************************************/

/* Global Variable */
int numSend = 0;	//The number of Send Event on each process
int numRecv = 0;	//the number of Recv Event on each process
int Drank;		//The rank of the current process
int totalSend;		//The number of Send Event on all processes
int totalRecv;		//The number of Recv Event on all processes
int sumSrc = 0;		//The xor(src) of all Recv Event on each process
int xorSrc;		//The xor(src) of all Recv Event on all processes
int sumDest = 0;	//The xor(dest) of all Send Event on each process
int xorDest;		//The xor(dest) of all Send Event on all processes 


/* MPI_Init Profiling Interface */
int MPI_DInit(int *argc, char ***argv) {
	int resultRank, result;
	result = PMPI_Init(argc, argv);
	resultRank = PMPI_Comm_rank(MPI_COMM_WORLD, &Drank);
	//totalSend = (int*) malloc (Drank * sizeof(int));
	return result;
}

/* MPI_Send Profiling Interface */
int MPI_DSend(void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
	numSend++;
	
	sumDest ^= (Drank + dest);
	return PMPI_Send(buf, count, type, dest, tag, comm);
}

/* MPI_Recv Profiling Interface */
int MPI_DRecv(void *buf, int count, MPI_Datatype type, int src, int tag, MPI_Comm comm, MPI_Status *status) {
	int result;
	numRecv++;
	
	result = PMPI_Recv(buf, count, type, src, tag, comm, status);
	//printf("MPI_SOURCE = %d\n", status->MPI_SOURCE); 
	sumSrc ^= (Drank + status->MPI_SOURCE);
	return result;
}

void beginfor() {

	numSend = 0;
	numRecv = 0;
}

void endfor(int iter, int rank) {
	printf("process %d enter endfor\n",rank);
	if (iter == 0) {
		int resultSend, resultRecv, totalSend, totalRecv;
		resultSend = PMPI_Reduce(&numSend, &totalSend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
		resultRecv = PMPI_Reduce(&numRecv, &totalRecv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		//PMPI_Barrier(MPI_COMM_WORLD);
		if(rank == 0) {
			printf("\nIter %d : totalSend = %d  totalRecv = %d\n", iter, totalSend, totalRecv);
		}
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int rank, size, i;
	int buf;
	MPI_Status status;

	MPI_DInit(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for(i = 0; i < 10; i++) {

		beginfor();

		if(i % 2 != 0) {

			if(rank >= 1 && rank <= size / 2) {
				MPI_DRecv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
			else {
				buf = rank;
				MPI_DSend(&buf, 1, MPI_INT, (rank - 1 + size) % size, rank, MPI_COMM_WORLD);
			}
			
		}
		else {
			if(rank < size / 2) {
				buf = rank;
				MPI_DSend(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
			}
			else {
				MPI_DRecv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		
		}

		endfor(i, rank);
	}

	MPI_Finalize();

	return 0;
}
