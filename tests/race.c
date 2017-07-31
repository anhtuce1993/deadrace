#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "Misc.h"


int main(int argc, char *argv[]) {
  int rank;
  int size;
  int recvbuf = 0;
  /*double starttime, endtime;*/

  beginning();
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

        //sleep(10);

  MPI_Status status;
  MPI_Request reqs [512];
  int outer_itr = 0;

  // for( outer_itr = 0; outer_itr < 2; outer_itr++ ) {
  
    if (rank == 0) {
      MPI_Recv(&recvbuf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout,"Process 0 has finished 1st recv \n");
      MPI_Recv(&recvbuf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout,"Process 0 has finished 2nd recv \n");
      MPI_Recv(&recvbuf, 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      /*MPI_Irecv(&recvbuf, 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, reqs);
      MPI_Wait(reqs, &status);*/

      fprintf(stdout,"Process 0 has finished 3rd recv \n");
      MPI_Send(&recvbuf, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
      fprintf(stdout,"Process 0 has finished 1st send \n");
      MPI_Send(&recvbuf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      fprintf(stdout,"Process 0 has finished 2nd send \n");
      MPI_Recv(&recvbuf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout,"Process 0 has finished 4th recv \n");
      MPI_Recv(&recvbuf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout,"Process 0 has finished 5th recv \n");
      MPI_Recv(&recvbuf, 1, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout,"Process 0 has finished 6th recv \n");
    }
    if (rank == 1) {
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Recv (&recvbuf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

    if (rank == 2) {
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Recv (&recvbuf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Send (&recvbuf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }
  
  // }
    MPI_Barrier(MPI_COMM_WORLD);

  /*fprintf(stdout,"Hello:In processor rank %d\n", rank);*/
  
  MPI_Finalize();

  ending();
}