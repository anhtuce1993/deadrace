#include "PMPI.h"


/* MPI_Init Profiling Interface */
int MPI_Init(int *argc, char ***argv) {
	/*printf("Enter init");*/
	int result;
	result = PMPI_Init(argc, argv);
	cTime = MPI_Wtime();
	PMPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	/*printf("\nRank : %d", myrank);*/
	PMPI_Comm_size(MPI_COMM_WORLD, &size);
	if (enabled) {
		lclk = 0;
		/*enabled = 0;*/
		if (myrank == rootrecv) {
			/*printf("\n Init Enter");*/
			fresult = fopen("result","w");
			initController(&controller);
		}
	}	
	return result;
}

/* MPI_Send Profiling Interface */
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
	if (enabled) {
		int num = sizeof (MPI_LONG_LONG_INT) + sizeof(datatype) * count;	
		int packsize = 0;
		char *packbuf = (char*) malloc (num);
		MPI_Pack (buf, count, datatype, packbuf, num, &packsize, comm);
		MPI_Pack (&lclk, 1, MPI_LONG_LONG_INT, packbuf, num, &packsize, comm);
		/*printf("\nProcess %i (send) : lclk = %i ", myrank, lclk);*/
		return PMPI_Send(packbuf, packsize, MPI_PACKED, dest, tag, comm);
	} else {
		return PMPI_Send(buf, count, datatype, dest, tag, comm);
	}
}


/*int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
	if (enabled) {
		int num = sizeof (int) + sizeof(datatype) * count;	
		int packsize = 0;
		char *packbuf = (char*) malloc (num);
		MPI_Pack (buf, count, datatype, packbuf, num, &packsize, comm);
		MPI_Pack (&lclk, 1, MPI_INT, packbuf, num, &packsize, comm);
		return PMPI_Isend(packbuf, packsize, MPI_PACKED, dest, tag, comm, request);
	} else {
		return PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
	}
}*/

/* MPI_Recv Profiling Interface */
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) 
{
	if (enabled) {
		/*int r, rt;
		fprintf(stdout,"\n % d : Enter", myrank);
		rt = PMPI_Comm_rank(MPI_COMM_WORLD, &r);
		printf("\n%d", rt);*/
		int result;
		int recvlclk=0;
		
		// unpack local clock piggypacking on receiving message
		int pos = 0;
		int num = sizeof (MPI_LONG_LONG_INT) + sizeof(datatype) * count;
		char *packbuf = (char*) malloc (num);
		PMPI_Recv (packbuf, num, MPI_PACKED, source, tag, comm, status);
		result = MPI_Unpack (packbuf, num, &pos, buf, count, datatype, comm);
		MPI_Unpack (packbuf, num, &pos, &recvlclk, 1, MPI_LONG_LONG_INT, comm);
		
		if (myrank == rootrecv) {
			// increase local clock when receiving on root process 
			//printf("Enter");
			/*fprintf(stdout,"\nEnter");*/
			lclk++;
			int from = (source == MPI_ANY_SOURCE) ? -1 : source;
			int src = status->MPI_SOURCE;
			printf("\nProcess %i (recv) : source = %i lclk = %i recvlclk = %i src = %i ", myrank, from,  lclk, recvlclk, src);
			int minPreRm;
			if (lclk % 100 == 0) {
				minPreRm = minPreRemove(controller,size,rootrecv);
				printf(" minPreRemove = %i ", minPreRm);
				removeRootRecvs(controller, minPreRm);
			}
			// add to root receiving list
			addRootRecv(controller, lclk, from);

			// compare local clock and receiving clock

				// add appropriate receiving local clock to src-process queue & remove inappropriate receiving local clock

				// remove selected receiving local clock for current receving event
			manipulateCommProcs(controller, src, recvlclk, lclk, fresult);

			//printProcRecvs(controller, src);

			int tempMem = getMemory();
			maxMem = (tempMem > maxMem) ? tempMem : maxMem;
			printf("%i KB \n", tempMem);
			
		} else {
			lclk = (recvlclk > lclk) ? recvlclk : lclk;
			/*printf("Enter");*/
			/*printf("\nProcess %i (recv) : lclk = %i ", myrank, lclk);*/
		}
		return result;
	} else {
		return PMPI_Recv(buf, count, datatype, source, tag, comm, status);
	}
}

/*int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) {
	int result;
	numRecv++;
		
	result = PMPI_Irecv(buf, count, type, src, tag, comm, request);
	//printf("MPI_SOURCE = %d\n", status->MPI_SOURCE); 
	//sumSrc ^= (Drank + status->MPI_SOURCE);
	return result;
}*/

int MPI_Barrier(MPI_Comm comm) {
	int rt = PMPI_Barrier(comm);
	if (enabled) {
		/*printf("\nProcess %i (barrier) : lclk = %i ", myrank, lclk);*/
		PMPI_Allreduce(&lclk, &lclk, 1, MPI_LONG_LONG_INT, MPI_MAX, comm);
		/*printf("lclkafter = %i ", lclk);*/
	}
	return rt;
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
	int rt = PMPI_Bcast(buffer, count, datatype, root, comm);
	if (enabled) {
		/*printf("\nProcess %i (bcast) : lclk = %i ", myrank, lclk);*/
		PMPI_Allreduce(&lclk, &lclk, 1, MPI_LONG_LONG_INT, MPI_MAX, comm);
		/*printf("lclkafter = %i ", lclk);*/
	}
	return rt;
}

int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
	int rt = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
	if (enabled) 
		PMPI_Allreduce(&lclk, &lclk, 1, MPI_LONG_LONG_INT, MPI_MAX, comm);
	return rt;
}

/* MPI_Finalize Profiling Interface */
int MPI_Finalize() {
	if (enabled) {
		PMPI_Barrier(MPI_COMM_WORLD);
		int i;
		if (myrank == rootrecv) {
			/*for(int i = 0; i < size; i++) {
				printProcRecvs(controller, i);
			}*/
			printRootRecvs(controller);
			
			/*fprintf(fresult, "\n");*/
			printf("\n\n-------------------------------DEADLOCK DETECTION RESULT------------------------------\n");
			for(i = 0; i < size; i++) {
				if (i != rootrecv) 
					checkRemainQueue(controller, i, fresult);
			}
			/*fprintf(fresult, "\n");*/
			fclose(fresult);
			printf("\n\n--------------------------------------SUMMARY-----------------------------------------\n");
			printf("\nNumber of receiving event on ROOT PROCESS : %i ", lclk);
			printf("\nMax Memory Consuming : %i KB", maxMem);
			printf("\n");
		}
	}
	cTime = MPI_Wtime() - cTime;
	double maxTime;
	PMPI_Reduce(&cTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, rootrecv, MPI_COMM_WORLD);
	if (myrank == rootrecv)
		printf("Max Consuming Time : %f \n", maxTime);
	return PMPI_Finalize();
}