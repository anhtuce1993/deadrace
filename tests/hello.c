#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

static int nsend = 0;

int MPI_Send(const void *start, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )
{
    nsend++;
    return PMPI_Send(start, count, datatype, dest, tag, comm);
}


int main(int argc, char *argv[]) {
	int myrank;
	int comm_size;
  int recvbuf = 0;
	
  MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

        //sleep(10);

  MPI_Status status;
  MPI_Request reqs [512];

  int recv_src = (myrank + comm_size - 1) % comm_size;
  int send_dest = (myrank + 1) % comm_size;
  
  int recv_count = 1;
  int send_count = 1;
  //changes for different count pattern
  /**if(myrank == 0)
    send_count = 2;
  if(recv_src ==  0)
    recv_count = 2;        
  **/

  int itr = 0;
  int loopcount = 10;  
  int outer_itr = 0;

  for( outer_itr = 0; outer_itr < 10; outer_itr++ )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Send (&recvbuf, send_count, MPI_INT, send_dest, 0, MPI_COMM_WORLD);
    MPI_Recv (&recvbuf, recv_count, MPI_INT, recv_src, 0, MPI_COMM_WORLD, &status);
    MPI_Barrier(MPI_COMM_WORLD);
  }//outer for loop
	fprintf(stdout,"Hello:In processor rank %d\n",myrank);
  MPI_Finalize();

  if (myrank == 0) {
    printf("\nnsend = %i \n", nsend);
  }
}
