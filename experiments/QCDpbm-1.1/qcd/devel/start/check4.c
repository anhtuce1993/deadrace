
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization of the spinor fields
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
   
#define NFIELDS 3

static float sig[NFIELDS];
static double sigd[NFIELDS];

int main(int argc,char *argv[])
{
   int my_rank,no_proc,k,i,j,ix;
   float d,dmax,dmax_all,*r,*rs;
   double *rd,var,var_all;
   spinor *ss;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&no_proc);   

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout); 
      printf("\n");
      printf("Initialization of the spinor fields\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,123456);   
   geometry();
   
   alloc_s(NFIELDS);
   dmax=0.0f;

   for (ix=0;ix<NSPIN;ix++)
   {
      for (k=0;k<NFIELDS;k++)
      {
         r=(float*)(ps[k][ix]);

         for (i=0;i<24;i++)
         {
            d=(float)fabs((double)(r[i]));
            if (d>dmax)
               dmax=d;
         }
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
   
   if (my_rank==0)
   {
      printf("Allocate %d single-precision spinor fields\n",NFIELDS);
      printf("|s| = %.2e\n\n",dmax_all);

      printf("Choose random single-precision spinor fields\n");   
   }

   if (my_rank==0)
      ranlxs(sig,NFIELDS);

   MPI_Bcast(sig,NFIELDS,MPI_FLOAT,0,MPI_COMM_WORLD);

   for (k=0;k<NFIELDS;k++)
   {
      random_s(VOLUME,ps[k][0],sig[k]);
      var=0.0;

      for (ix=0;ix<VOLUME;ix++)
      {
         r=(float*)(ps[k][ix]);

         for (i=0;i<24;i++)
            var+=(double)(r[i]*r[i]);
      }

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      
      if (my_rank==0)
      {
         var_all/=(double)(12*VOLUME*no_proc);               
         printf("<s[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sig[k]*sig[k]);
      }
   }

   alloc_sd(NFIELDS);
   dmax=0.0f;

   for (ix=0;ix<NSPIN;ix++)
   {
      for (k=0;k<NFIELDS;k++)
      {
         rd=(double*)(psd[k][ix]);

         for (i=0;i<24;i++)
         {
            d=(float)fabs(rd[i]);
            if (d>dmax)
               dmax=d;
         }
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
   
   if (my_rank==0)
   {   
      printf("\n");
      printf("Allocate %d double-precision spinor fields\n",NFIELDS);
      printf("|sd| = %.2e\n\n",dmax_all);

      printf("Choose random double-precision spinor fields\n");   
   }
   
   if (my_rank==0)
      ranlxd(sigd,NFIELDS);

   MPI_Bcast(sigd,NFIELDS,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   for (k=0;k<NFIELDS;k++)
   {
      random_sd(VOLUME,psd[k][0],sigd[k]);
      var=0.0;

      for (ix=0;ix<VOLUME;ix++)
      {
         rd=(double*)(psd[k][ix]);

         for (i=0;i<24;i++)
            var+=rd[i]*rd[i];
      }

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      
      if (my_rank==0)
      {
         var_all/=(double)(12*VOLUME*no_proc);               
         printf("<sd[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sigd[k]*sigd[k]);
      }
   }

   ss=malloc(NFIELDS*VOLUME*sizeof(spinor));
   i=0;

   for (k=0;k<NFIELDS;k++)
   {
      for (ix=0;ix<VOLUME;ix++)
         ss[i++]=*(ps[k][ix]);
   }   

   for (k=0;k<NFIELDS;k++)
      assign_s2sd(VOLUME,ps[k][0],psd[NFIELDS-1-k][0]);

   free_s();
   alloc_s(NFIELDS);

   for (k=0;k<NFIELDS;k++)
      assign_sd2s(VOLUME,psd[k][0],ps[NFIELDS-1-k][0]);   

   i=0;
   dmax=0.0f;

   for (k=0;k<NFIELDS;k++)
   {
      for (ix=0;ix<VOLUME;ix++)
      {
         r=(float*)(ps[k][ix]);
         rs=(float*)(&ss[i++]);

         for (j=0;j<24;j++)
         {
            d=(float)fabs((double)(r[j]-rs[j]));
            if (d>dmax)
               dmax=d;
         }
      }
   }    

   MPI_Reduce(&dmax,&dmax_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
   
   if (my_rank==0)
   { 
      printf("\n");
      printf("Assign s to sd, free s, reassign s from sd\n");
      printf("|s_original-s_assigned| = %.2e\n\n",dmax_all);

      fclose(flog);
   }

   MPI_Finalize();   
   exit(0);
}
