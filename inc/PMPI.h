#ifndef __PMPI_H__
#define __PMPI_H__

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#ifdef __cplusplus

#include "Controller.h"
#include "Memory.h"

#define rootrecv	0

/* Global Variable */
int myrank;		//The rank of the current process
int size;		//The number of processes

static Controller* controller;

static int lclk = 0;

int enabled = 0;

int maxMem = 0;

double cTime;

FILE* fresult = NULL;

extern int MPI_Init(int *argc, char ***argv);

extern int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

/*extern int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);*/

extern int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

// extern int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);

extern int MPI_Barrier(MPI_Comm comm);

extern int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

extern int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

extern int MPI_Finalize();

#endif /* __cplusplus */

#endif /* __PMPI_H__ */
