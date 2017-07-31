
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization of the link variables
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


static complex det(su3 *u)
{
   complex det1,det2,det3,detu;

   det1.re=
      ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
      ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
      ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
      ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);
   det3.re=
      ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
      ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
      ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
      ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   detu.re=
      ((*u).c11.re*det1.re-(*u).c11.im*det1.im)-
      ((*u).c12.re*det2.re-(*u).c12.im*det2.im)+
      ((*u).c13.re*det3.re-(*u).c13.im*det3.im);
   detu.im=
      ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
      ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
      ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   return detu;
}


static complex_dble det_dble(su3_dble *u)
{
   complex_dble det1,det2,det3,detu;

   det1.re=
      ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
      ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
      ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
      ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);
   det3.re=
      ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
      ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
      ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
      ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   detu.re=
      ((*u).c11.re*det1.re-(*u).c11.im*det1.im)-
      ((*u).c12.re*det2.re-(*u).c12.im*det2.im)+
      ((*u).c13.re*det3.re-(*u).c13.im*det3.im);
   detu.im=
      ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
      ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
      ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   return detu;
}


static float dev_uudag(su3 *u,su3 *v)
{
   int i;
   float *r,d,dmax;
   su3 vdag,w;

   _su3_dagger(vdag,(*v));
   _su3_times_su3(w,(*u),vdag);

   w.c11.re-=1.0f;
   w.c22.re-=1.0f;
   w.c33.re-=1.0f;

   r=(float*)(&w);
   dmax=0.0f;

   for (i=0;i<18;i++)
   {
      d=(float)fabs((double)(r[i]));
      if (d>dmax)
         dmax=d;
   }

   return dmax;
}


static double dev_uudag_dble(su3_dble *u,su3_dble *v)
{
   int i;
   double *r,d,dmax;
   su3_dble vdag,w;

   _su3_dagger(vdag,(*v));
   _su3_times_su3(w,(*u),vdag);

   w.c11.re-=1.0;
   w.c22.re-=1.0;
   w.c33.re-=1.0;

   r=(double*)(&w);
   dmax=0.0;

   for (i=0;i<18;i++)
   {
      d=fabs(r[i]);
      if (d>dmax)
         dmax=d;
   }

   return dmax;
}


static float dev_detu(su3 *u)
{
   float d,dmax;
   complex detu;

   detu=det(u);
   dmax=0.0f;

   d=(float)fabs((double)(1.0f-detu.re));
   if (d>dmax)
      dmax=d;
   d=(float)fabs((double)(detu.im));
   if (d>dmax)
      dmax=d;

   return dmax;
}


static double dev_detu_dble(su3_dble *u)
{
   double d,dmax;
   complex_dble detu;

   detu=det_dble(u);
   dmax=0.0;

   d=fabs(1.0-detu.re);
   if (d>dmax)
      dmax=d;
   d=fabs(detu.im);
   if (d>dmax)
      dmax=d;

   return dmax;
}


int main(int argc,char *argv[])
{
   int my_rank,i,j;
   float *r,*s,d1,d2,dmax1,dmax2,dmax1_all,dmax2_all;
   double *rd;
   su3 *u,*ub,*su,v;
   su3_dble *ud,*udb,vd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Initialization of the link variables\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,123456);
   geometry();
   alloc_u();

   ub=pu[VOLUME/2][0];
   dmax1=0.0f;
   r=(float*)(&v);

   for (u=ub;u<(ub+4*VOLUME+(BNDRY/4));u++)
   {
      v=*u;
      v.c11.re-=1.0f;
      v.c22.re-=1.0f;
      v.c33.re-=1.0f;

      for (i=0;i<18;i++)
      {
         d1=(float)fabs((double)(r[i]));
         if (d1>dmax1)
            dmax1=d1;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocate single-precision gauge field\n");
      printf("|u-1| = %.2e\n\n",dmax1_all);
   }

   random_u();
   dmax1=0.0f;
   dmax2=0.0f;

   for (u=ub;u<(ub+4*VOLUME);u++)
   {
      d1=dev_uudag(u,u);
      d2=dev_detu(u);

      if (d1>dmax1)
         dmax1=d1;
      if (d2>dmax2)
         dmax2=d2;
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Call random_u\n");
      printf("|u^dag*u-1| = %.2e\n",dmax1_all);
      printf("|det{u}-1| = %.2e\n\n",dmax2_all);
   }

   alloc_ud();
   udb=pud[VOLUME/2][0];
   dmax1=0.0f;
   rd=(double*)(&vd);

   for (ud=udb;ud<(udb+4*VOLUME+(BNDRY/4));ud++)
   {
      vd=*ud;
      vd.c11.re-=1.0;
      vd.c22.re-=1.0;
      vd.c33.re-=1.0;

      for (i=0;i<18;i++)
      {
         d1=(float)fabs(rd[i]);
         if (d1>dmax1)
            dmax1=d1;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocate double-precision gauge field\n");
      printf("|ud-1| = %.2e\n\n",dmax1_all);
   }

   assign_u2ud();
   dmax1=0.0f;
   dmax2=0.0f;

   for (ud=udb;ud<(udb+4*VOLUME);ud++)
   {
      d1=(float)dev_uudag_dble(ud,ud);
      d2=(float)dev_detu_dble(ud);

      if (d1>dmax1)
         dmax1=d1;
      if (d2>dmax2)
         dmax2=d2;
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Assign single-precision to double-precision field\n");
      printf("|ud^dag*ud-1| = %.2e\n",dmax1_all);
      printf("|det{ud}-1| = %.2e\n\n",dmax2_all);
   }

   su=malloc(4*VOLUME*sizeof(su3));
   i=0;

   for (u=ub;u<(ub+4*VOLUME);u++)
      su[i++]=*u;

   free_u();
   alloc_u();
   assign_ud2u();

   ub=pu[VOLUME/2][0];
   dmax1=0.0f;
   i=0;

   for (u=ub;u<(ub+4*VOLUME);u++)
   {
      r=(float*)(&su[i++]);
      s=(float*)(u);

      for (j=0;j<18;j++)
      {
         d1=(float)fabs((double)(r[j]-s[j]));
         if (d1>dmax1)
            dmax1=d1;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Reallocate single-precision field\n");
      printf("Assign double-precision to single-precision field\n");
      printf("|u_original-u_assigned| = %.2e\n\n",dmax1_all);
   }

   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
