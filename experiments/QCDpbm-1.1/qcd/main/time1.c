
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* QCD single-precision speed test
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
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"


static double wt_norm_square(int nflds)
{
   int my_rank,nmax,n,i,ib;
   float r;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);   
   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i++)
            r=norm_square(VOLUME,1,ps[i][0]);
      }
      
      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);      

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)(nmax*nflds);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);      
   
   return wtav;
}


static double wt_mulc_spinor_add(int nflds)
{
   int my_rank,nmax,n,i,ib;
   complex z;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);   
   z.re=0.123f;
   z.im=0.456f;
   nmax=1;
   
   for (ib=0;ib<1;nmax*=2)
   {   
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i+=2)
            mulc_spinor_add(VOLUME,ps[i][0],ps[i+1][0],z);

         z.re-=z.re;
         z.im-=z.im;
      }
      
      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);      

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)((nmax*nflds)/2);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return wtav;
}


static double wt_Qhat(int nflds)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);   
   nmax=1;
   
   for (ib=0;ib<1;nmax*=2)
   { 
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i+=2)
            Qhat(i,i+1);
      }
      
      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);      

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)((nmax*nflds)/2);
      }
      

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return wtav;
}


int main(int argc,char *argv[])
{
   int my_rank,nflds,n;
   double wdt0,wdt1,wdt2,wdt3;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);
      error_root(flog==NULL,1,"main [time1.c]","Unable to open log file");

      printf("\n");
      printf("QCD single-precision speed test\n");
      printf("-------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      if (NPROC>1)
         printf("There are %d MPI processes\n",NPROC);
      else
         printf("There is 1 MPI process\n");
      
      if ((VOLUME*sizeof(float))<(64*1024))
      {      
         printf("The local size of the gauge field is %d KB\n",
                (int)((72*VOLUME*sizeof(float))/(1024)));
         printf("The local size of a quark field is %d KB\n",
                (int)((24*VOLUME*sizeof(float))/(1024)));
      }
      else
      {
         printf("The local size of the gauge field is %d MB\n",
                (int)((72*VOLUME*sizeof(float))/(1024*1024)));
         printf("The local size of a quark field is %d MB\n",
                (int)((24*VOLUME*sizeof(float))/(1024*1024)));
      }

#if (defined SSE3)
      printf("Using inline assembly SSE3 instructions\n");
#elif (defined SSE2)
      printf("Using inline assembly SSE2 instructions\n");      
#elif (defined SSE)
      printf("Using inline assembly SSE instructions\n");
#endif

#if (defined SSE)
#if (defined P3)
      printf("Assuming SSE prefetch instructions fetch 32 bytes\n");
#elif (defined PM)
      printf("Assuming SSE prefetch instructions fetch 64 bytes\n");
#elif (defined P4)
      printf("Assuming SSE prefetch instructions fetch 128 bytes\n");
#else
      printf("SSE prefetch instructions are not used\n");
#endif
#endif
      
      printf("\n");
   }

   start_ranlux(0,12345);
   geometry();
   alloc_u();
   alloc_ud();
   alloc_sw();
   alloc_swd();

   random_ud();
   assign_ud2u();
   set_lat_parms(0.0,0.1,0.2);
   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [time1.c]",
         "Inversion of swd on the odd sites was not safe");
   assign_swd2sw();

   nflds=(int)((4*1024*1024)/(VOLUME*sizeof(float)))+1;
   if ((nflds%2)==1)
      nflds+=1;
   alloc_s(nflds);

   for (n=0;n<nflds;n++)
      random_s(VOLUME,ps[n][0],1.0f);

   error_chk();   
   wdt0=1.0e6*wt_norm_square(nflds)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Program norm_square:\n");
      printf("Time per lattice point: %4.3f usec",wdt0);
      printf(" (%d Mflops/process)\n\n",(int)(48.0/wdt0));
      fflush(flog);
   }
   
   wdt1=1.0e6*wt_mulc_spinor_add(nflds)/(double)(VOLUME);   

   if (my_rank==0)
   {
      printf("Program mulc_spinor_add:\n");
      printf("Time per lattice point: %4.3f usec",wdt1);
      printf(" (%d Mflops/process)\n\n",(int)(96.0/wdt1));
      fflush(flog);
   }

   wdt2=1.0e6*wt_Qhat(nflds)/(double)(VOLUME);   

   if (my_rank==0)
   {
      printf("Program Qhat:\n");
      printf("Time per lattice point: %4.3f usec",wdt2);
      printf(" (%d Mflops/process)\n\n",(int)(1896.0/wdt2));
      fflush(flog);
   }

   wdt3=2.0*wdt0+3.0*wdt1+2.0*wdt2;

   if (my_rank==0)
   {
      printf("########################################################\n");
      printf("#                                                      #\n");
      printf("#             SYNTHETIC QCD SPEED TEST                 #\n");
      printf("#                                                      #\n");
      printf("#  Using single-precision (%d bit) data and programs   #\n",
             8*(int)(sizeof(float)));
      printf("#                                                      #\n");
      printf("#   Time per lattice point:  %8.3f usec             #\n",
             wdt3);
      printf("#   Average speed:           %8.3f Gflops/process   #\n",
             1.0e-3*4176.0/wdt3);
      printf("#   Total throughput:        %8.3f Gflops           #\n",
             1.0e-3*(double)(NPROC)*4176.0/wdt3);
      printf("#                                                      #\n");
      printf("########################################################\n\n");
      fclose(flog);
   }   
   
   MPI_Finalize();
   exit(0);
}
