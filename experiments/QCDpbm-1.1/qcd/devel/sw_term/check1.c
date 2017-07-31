
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation, assignment and inversion of the global pauli arrays
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "sw_term.h"
#include "start.h"
#include "global.h"

static pauli_dble *sswd=NULL;
static const weyl_dble vd0={{{0.0}}};

#if (defined SSE2)
static weyl_dble vd __attribute__ ((aligned (16)));
#else
static weyl_dble vd;
#endif


static void random_swd(void)
{
   int k;
   pauli_dble *pa,*pam;

   pa=swd;
   pam=swd+2*VOLUME;

   for (;pa<pam;pa++)
   {
      ranlxd((*pa).u,36);

      for (k=0;k<36;k++)
      {
         if (k<6)
            (*pa).u[k]=1.0+0.05*(0.5-(*pa).u[k]);
         else
            (*pa).u[k]=0.05*(0.5-(*pa).u[k]);
      }
   }
}


static void save_swd(void)
{
   pauli_dble *pa,*pb,*pam;

   if (sswd==NULL)
   {
      sswd=amalloc(2*VOLUME*sizeof(*sswd),ALIGN);
      error(sswd==NULL,1,"save_swd [check1.c]",
            "Unable to allocate auxiliary array"); 
   }
   
   pa=swd;
   pb=sswd;
   pam=swd+2*VOLUME;

   for (;pa<pam;pa++)
   {
      *pb=(*pa);
      pb+=1;
   }
}


static double cmp_swd(ptset_t set)
{
   int k;
   double d,dmax;
   pauli_dble *pa,*pb,*pam;

   pa=swd;
   pb=sswd;
   pam=pa;

   if (set==EVEN_PTS)
      pam=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pam=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pam=pa+2*VOLUME;

   dmax=0.0;

   for (;pa<pam;pa++)
   {
      for (k=0;k<36;k++)
      {
         d=fabs((*pa).u[k]-(*pb).u[k]);

         if (d>dmax)
            dmax=d;
      }

      pb+=1;
   }

   return dmax;
}


static double cmp_iswd(ptset_t set)
{
   int k,l;
   double d,dmax;
   complex_dble *cvd;
   pauli_dble *pa,*pb,*pam;

   pa=swd;
   pb=sswd;
   pam=pa;

   if (set==EVEN_PTS)
      pam=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pam=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pam=pa+2*VOLUME;

   dmax=0.0;
   cvd=(complex_dble*)(&vd);

   for (;pa<pam;pa++)
   {
      for (k=0;k<6;k++)
      {
         vd=vd0;
         cvd[k].re=1.0;

         mul_pauli_dble(pa,&vd,&vd);
         mul_pauli_dble(pb,&vd,&vd);
         cvd[k].re-=1.0;

         for (l=0;l<6;l++)
         {
            d=cvd[l].re*cvd[l].re+cvd[l].im*cvd[l].im;
            if (d>dmax)
               dmax=d;
         }
      }

      pb+=1;
   }

   return sqrt(dmax);
}


static double cmp_sw2swd(ptset_t set)
{
   int k;
   double d,dmax;
   pauli *pa,*pam;
   pauli_dble *pb;

   pa=sw;
   pb=swd;
   pam=pa;

   if (set==EVEN_PTS)
      pam=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pam=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pam=pa+2*VOLUME;

   dmax=0.0;

   for (;pa<pam;pa++)
   {
      for (k=0;k<36;k++)
      {
         d=fabs((double)((*pa).u[k])-(*pb).u[k]);

         if (d>dmax)
            dmax=d;
      }

      pb+=1;
   }

   return dmax;
}


int main(int argc,char *argv[])
{
   int my_rank,no_proc;
   int ix,k,itest;
   float *r;
   double *rd,d,dmax,dmax_all;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&no_proc);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      printf("\n");
      printf("Initialization of the global pauli arrays\n");
      printf("-----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,123456);

   alloc_sw();
   dmax=0.0;

   for (ix=0;ix<(2*VOLUME);ix++)
   {
      r=(float*)(sw+ix);

      for (k=0;k<36;k++)
      {
         if (k<6)
            d=fabs((double)(r[k]-1.0f));
         else
            d=fabs((double)(r[k]));

         if (d>dmax)
            dmax=d;
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocated global single-precision pauli field\n");
      printf("max|p-1| = %.1e\n\n",dmax_all);
   }

   alloc_swd();
   dmax=0.0;

   for (ix=0;ix<(2*VOLUME);ix++)
   {
      rd=(double*)(swd+ix);

      for (k=0;k<36;k++)
      {
         if (k<6)
            d=fabs(rd[k]-1.0);
         else
            d=fabs(rd[k]);

         if (d>dmax)
            dmax=d;
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocated global double-precision pauli field\n");
      printf("max|p-1| = %.1e\n\n",dmax_all);
   }

   random_swd();
   save_swd();

   itest=invert_swd(EVEN_PTS);
   error(itest!=0,1,"main [check1.c]","Unsafe inversion of swd_e");
   dmax=cmp_iswd(EVEN_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd_e\n");
      printf("max deviation swd_e = %.1e\n",dmax_all);
   }

   dmax=cmp_swd(ODD_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("max deviation swd_o = %.1e\n\n",dmax_all);

   random_swd();
   save_swd();

   itest=invert_swd(ODD_PTS);
   error(itest!=0,1,"main [check1.c]","Unsafe inversion of swd_o");
   dmax=cmp_swd(EVEN_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd_o\n");
      printf("max deviation swd_e = %.1e\n",dmax_all);
   }

   dmax=cmp_iswd(ODD_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("max deviation swd_o = %.1e\n\n",dmax_all);

   assign_swd2sw();
   dmax=cmp_sw2swd(ALL_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Assigned swd to sw\n");
      printf("max deviation = %.1e\n\n",dmax_all);
   }

   random_swd();
   save_swd();
   itest=invert_swd(ALL_PTS);
   error(itest!=0,1,"main [check1.c]","Unsafe inversion of swd_eo");
   dmax=cmp_iswd(ALL_PTS);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd_eo\n");
      printf("max deviation = %.1e\n",dmax_all);
   }

   free_sw();
   free_swd();

   if (my_rank==0)
   {
      printf("Freed sw and swd\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
