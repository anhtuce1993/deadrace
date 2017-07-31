/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Random number generator initialization in MPI programs
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "start.h"
#include "global.h"

#define NRAN 1000

static float r[NRAN];


int main(int argc,char *argv[])
{
   int my_rank,no_proc;
   int level,seed;
   int n,i,j,k,l,itest,tag0;
   int *s=NULL,*c=NULL,*sa=NULL;
   FILE *flog=NULL;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&no_proc);

   if (my_rank==0) {
      /*flog=freopen("check5.log","w",stdout);      */
      flog = fopen("check5.log","w");
      
      /*printf("\n");
      printf("Test of the random number initialization\n");
      printf("----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);*/

      fprintf(flog,"\n");
      fprintf(flog,"Test of the random number initialization\n");
      fprintf(flog,"----------------------------------------\n\n");

      fprintf(flog,"%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      fprintf(flog,"%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      fprintf(flog,"%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

   }

   level=0;
   seed=12345;
   start_ranlux(level,seed);
   n=rlxs_size();

   s=malloc(n*sizeof(int));
   c=malloc(n*sizeof(int));   

   ranlxs(r,NRAN);
   rlxs_get(s);

   if (my_rank==0)
   {
      sa=malloc(n*no_proc*sizeof(int));
      for (i=0;i<n;i++)
         sa[i]=s[i];
   }

   for (i=1;i<no_proc;i++)
   {
      tag0=mpi_tag();
      
      if (my_rank==i)
         MPI_Send(&s[0],n,MPI_INT,0,tag0,MPI_COMM_WORLD);
      else if (my_rank==0) {
         /*MPI_Recv(&sa[n*i],n,MPI_INT,i,tag0,MPI_COMM_WORLD,&status);*/
         /*if (i == 2) 
            MPI_Recv(&sa[n*i],n,MPI_INT,i,tag0,MPI_COMM_WORLD,&status);
         else    */
            MPI_Recv(&sa[n*i],n,MPI_INT,MPI_ANY_SOURCE,tag0,MPI_COMM_WORLD,&status);
      }
   }

   if (my_rank==0)
   {
      /*printf("The corresponding state vectors are\n");*/
      fprintf(flog,"The corresponding state vectors are\n");
      
      for (i=0;4*i<no_proc;i++)
      {
         for (j=0;j<n;j++)
         {
            /*printf("\n");*/
            fprintf(flog,"\n");
            for (k=0;k<4;k++)
            {
               l=n*(4*i+k)+j;
               if (l<(n*no_proc))
                  /*printf("\t%10d",sa[l]);*/
                  fprintf(flog,"\t%10d",sa[l]);
            }
         }
         /*printf("\n");*/
         fprintf(flog,"\n");
      }
   }

   for (i=1;i<no_proc;i++)
   {
      tag0=mpi_tag();
      
      if (my_rank==0)
         MPI_Send(&sa[n*i],n,MPI_INT,i,tag0,MPI_COMM_WORLD);
      else if (my_rank==i)
         MPI_Recv(&c[0],n,MPI_INT,0,tag0,MPI_COMM_WORLD,&status);
   }

   if (my_rank==0)
      for (i=0;i<n;i++)
         c[i]=s[i];

   itest=0;
   
   for (i=0;i<n;i++)
   {
      if (s[i]!=c[i])
         itest=1;
   }

   MPI_Allreduce(&itest,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
   itest=i;
   
   error(itest!=0,1,"main program [check5.c]",
         "Read back test failed");
   
   if (my_rank==0)
   {
      /*printf("\n");
      printf("Read back test passed\n\n");*/
      fprintf(flog,"\n");
      fprintf(flog,"Read back test passed\n\n");
      fclose(flog);
   }
   
   MPI_Finalize();
   exit(0);
}
