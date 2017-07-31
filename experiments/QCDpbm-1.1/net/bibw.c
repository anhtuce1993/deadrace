
/*******************************************************************************
*
* File bibw.c
*
* Copyright (C) 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Node-to-node bi-directional bandwidth measurement
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "net.h"


static int bw0(int buf_size,double *buf)
{
   int my_rank,i,imax,ib;
   double bsmb;
   double wt1,wt2,wdt;
   MPI_Status status;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   bsmb=(double)(buf_size*sizeof(double))/(1024.0*1024.0);
   imax=(int)(128.0/bsmb);
   if (imax==0)
      imax=1;
   wdt=1.0;

   for (ib=0;ib<1;imax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
               
      if (my_rank==0)
      {
         wt1=MPI_Wtime();      
         for (i=0;i<imax;i++)
            MPI_Send(buf,buf_size,MPI_DOUBLE,2,i,MPI_COMM_WORLD);
         for (i=0;i<imax;i++)
            MPI_Recv(buf,buf_size,MPI_DOUBLE,2,i,MPI_COMM_WORLD,
                     &status);
         wt2=MPI_Wtime();               
         wdt=wt2-wt1;
         if (wdt>=2.0)
            ib=1;
      }
      else if (my_rank==2)
      {
         for (i=0;i<imax;i++)
            MPI_Recv(buf,buf_size,MPI_DOUBLE,0,i,MPI_COMM_WORLD,
                     &status);
         for (i=0;i<imax;i++)
            MPI_Send(buf,buf_size,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   if (my_rank==0)
      ib=(int)((double)(imax)*bsmb/wdt);

   MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);   
   
   return ib;
}


int main(int argc,char* argv[])
{
   int my_rank,no_proc;
   int max_size,buf_size,bw;
   int n,i,imax,p1,p2,p1w,p2w,p[4];
   double *buf,wt1,wt2,wdt;
   double wtmn,wtmx,wtav;
   FILE *flog=NULL;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&no_proc);

   max_size=1024*MAX_PACKET_SIZE;
   buf=(double*)(alloc_buffer(max_size));
   buf_size=(int)(max_size/sizeof(double));
   init_buffer(buf_size,buf);
   bw=bw0(buf_size,buf);

   if (my_rank==0)
   {
      flog=freopen("bibw.log","w",stdout);
      error_root(flog==NULL,1,"main [bibw.c]","Unable to open log file");
      
      printf("\n");
      printf("Node-to-node bi-directional bandwidth measurement\n");
      printf("-------------------------------------------------\n\n");
      printf("There are %d MPI processes.\n\n",no_proc);
      printf("Each entry in the table reports the measured value averaged\n"
             "over all possible pairs of nodes. In the case of the bandwidth,\n"
             "the minimal and maximal values are given in square brackets."
             "\n\n");
      printf("Packet [KB]  Time/Test [s]  Time/Packet [ms]  "
             "Bandwidth [MB/s]   Worst Channel\n");
      fflush(flog);

      error_root(((no_proc%2)!=0)||(no_proc<4),2,"main [bibw.c]",
                 "The number of processes must be even and at least 4 "
                 "for this test");
   }

   for (n=MIN_PACKET_SIZE;n<=MAX_PACKET_SIZE;n*=2)
   {
      buf_size=(int)((1024*n)/sizeof(double));
      imax=(1024*bw)/n+1;
      wtmn=0.0;
      wtmx=0.0;
      wtav=0.0;
      p1w=0;
      p2w=2;

      for (p1=0;p1<no_proc;p1+=2)
      {
         p[0]=p1;
         p[1]=p1+1;
         
         for (p2=(p1+2);p2<no_proc;p2+=2)
         {
            p[2]=p2;
            p[3]=p2+1;
            
            MPI_Barrier(MPI_COMM_WORLD);
            wt1=MPI_Wtime();
            
            if (my_rank==p1)
            {
               for (i=0;i<imax;i++)
                  MPI_Send(buf,buf_size,MPI_DOUBLE,p2,i,MPI_COMM_WORLD);
               for (i=0;i<imax;i++)
                  MPI_Recv(buf,buf_size,MPI_DOUBLE,p2,i,MPI_COMM_WORLD,
                           &status);
            }
            else if (my_rank==p2)
            {
               for (i=0;i<imax;i++)
                  MPI_Recv(buf,buf_size,MPI_DOUBLE,p1,i,MPI_COMM_WORLD,
                           &status);
               for (i=0;i<imax;i++)
                  MPI_Send(buf,buf_size,MPI_DOUBLE,p1,i,MPI_COMM_WORLD);
            }
            else if (my_rank==(p1+1))
            {
               for (i=0;i<imax;i++)
                  MPI_Recv(buf,buf_size,MPI_DOUBLE,p2+1,i,MPI_COMM_WORLD,
                           &status);
               for (i=0;i<imax;i++)
                  MPI_Send(buf,buf_size,MPI_DOUBLE,p2+1,i,MPI_COMM_WORLD);
            }
            else if (my_rank==(p2+1))
            {
               for (i=0;i<imax;i++)
                  MPI_Send(buf,buf_size,MPI_DOUBLE,p1+1,i,MPI_COMM_WORLD);
               for (i=0;i<imax;i++)
                  MPI_Recv(buf,buf_size,MPI_DOUBLE,p1+1,i,MPI_COMM_WORLD,
                           &status);
            }  
            
            wt2=MPI_Wtime();
            wdt=wt2-wt1;
            MPI_Barrier(MPI_COMM_WORLD);            

            for (i=0;i<4;i++)
            {
               if (p[i]!=0)
               {
                  if (my_rank==p[i])
                     MPI_Send(&wdt,1,MPI_DOUBLE,0,imax+i,MPI_COMM_WORLD);
                  if (my_rank==0)
                     MPI_Recv(&wdt,1,MPI_DOUBLE,p[i],imax+i,MPI_COMM_WORLD,
                              &status);
               }
               
               if (my_rank==0)
               {
                  if (wtmn==0.0)
                     wtmn=wdt;
                  else if (wtmn>wdt)
                     wtmn=wdt;

                  if (wtmx<wdt)
                  {
                     p1w=p1;
                     p2w=p2;
                     wtmx=wdt;
                  }

                  wtav+=wdt;
               }
            }
         }
      }
      
      if (my_rank==0)
      {
         wtav/=(double)(no_proc*((no_proc/2)-1));
         wtmn*=0.5;
         wtmx*=0.5;
         wtav*=0.5; 

         printf("%6d",n);
         printf("          %1.1f",wtav);
         printf("            %.1e",1.0e3*wtav/(double)(imax));
         printf("        %4d [%4d,%4d]",
                (int)((double)(2*imax*n)/(wtav*1024.0)),
                (int)((double)(2*imax*n)/(wtmx*1024.0)),
                (int)((double)(2*imax*n)/(wtmn*1024.0)));
         printf("  %2d,%2d <> %2d,%2d\n",p1w,p1w+1,p2w,p2w+1);
         fflush(flog);
      }
   }

   if (my_rank==0)
   {
      printf("\nNote: 1 KB = 1024 Byte, 1 MB = 1024 KB. The figures in\n"
             "the last column are the ranks of the MPI processes where\n"
             "the lowest bandwidth was measured.\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
