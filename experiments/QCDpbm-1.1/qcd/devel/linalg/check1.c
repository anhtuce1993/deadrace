
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the programs in the module linalg
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
#include "misc.h"
#include "linalg.h"
#include "global.h"


static complex sp(int vol,spinor *pk,spinor *pl)
{
   int ix;
   double x,y;
   complex *rpk,*rpl,z;

   x=0.0;
   y=0.0;

   rpk=(complex*)(pk);
   rpl=(complex*)(pl);
   
   for (ix=0;ix<(12*vol);ix++)
   {
      x+=(double)((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im);
      y+=(double)((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re);
      rpk+=1;
      rpl+=1;
   }
   
   z.re=(float)(x);
   z.im=(float)(y);
   
   return z;
}


int main(int argc,char *argv[])
{
   int my_rank,i,vol,off;
   int icom,ieo;
   float r,zsq;
   double d,dmax;
   complex z,w;
   spinor *pk,*pl;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      
      printf("\n");
      printf("Consistency of the programs in the module linalg\n");
      printf("------------------------------------------------\n\n");   

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   for (icom=0;icom<2;icom++)
   {
      if ((icom==0)||(NPROC>1))
      {
         if (my_rank==0)
         {
            if (icom==1)
            {
               printf("Checks with global summation\n");
               printf("============================\n\n");
            }
            else
            {
               printf("Checks without global summation\n");
               printf("===============================\n\n");
            }
         }
      
         for (ieo=0;ieo<3;ieo++)
         {
            if (my_rank==0)
            {
               if (ieo==0)
                  printf("First case: full lattice\n\n");
               else if (ieo==1)
                  printf("Second case: even points\n\n");
               else
                  printf("Third case: odd points\n\n");
            }

            vol=VOLUME/2;
            off=0;

            if (ieo==0)
               vol=VOLUME;
            if (ieo==2)
               off=VOLUME/2;
   
            start_ranlux(0,12345);
            geometry();
            alloc_s(10);

            for (i=0;i<10;i++)
               random_s(vol,ps[i][off],1.0f);

            dmax=0.0;
   
            for (i=0;i<10;i++)
            {
               pk=ps[i][off];
               pl=ps[9-i][off];

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_FLOAT,0,MPI_COMM_WORLD); 
               }
               else
                  w=sp(vol,pk,pl);

               z=spinor_prod(vol,icom,pk,pl);
               r=norm_square(vol,icom,pk)*norm_square(vol,icom,pl);
               d=(double)((z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im));
               d=sqrt(d/(double)(r));
               if (d>dmax)
                  dmax=d;

               z=spinor_prod(vol,icom,pk,pk);
               r=norm_square(vol,icom,pk);
      
               d=fabs((double)(z.im/r));
               if (d>dmax)
                  dmax=d;

               d=fabs((double)(z.re/r-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
               printf("Check of spinor_prod and norm_square: %.2e\n\n",dmax);
   
            dmax=0.0;
            z.re= 0.345f;
            z.im=-0.876f;
            zsq=z.re*z.re+z.im*z.im;
   
            for (i=0;i<9;i++)
            {
               pk=ps[i][off];
               pl=ps[i+1][off];      
      
               w=spinor_prod(vol,icom,pk,pl);
               r=norm_square(vol,icom,pk)+zsq*norm_square(vol,icom,pl)
                  +2.0f*(z.re*w.re-z.im*w.im);
               mulc_spinor_add(vol,pk,pl,z);

               d=fabs((double)(r/norm_square(vol,icom,pk)-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {   
               printf("Consistency of spinor_prod, norm_square\n");
               printf("and mulc_spinor_add: %.2e\n\n",dmax);
            }
   
            for (i=0;i<10;i++)
               random_s(vol,ps[i][off],1.0f);

            dmax=0.0;
   
            for (i=0;i<8;i+=2)
            {
               pk=ps[i][off];
               pl=ps[i+1][off];
               assign_s2s(vol,pk,pl);

               r=normalize(vol,icom,pk);
               z=spinor_prod(vol,icom,pk,pl);

               d=fabs(1.0-fabs((double)(z.re))/(double)(r))+
                  fabs((double)(z.im))/r;
               if (d>dmax)
                  dmax=d;
               
               z.re=-z.re;
               z.im=-z.im;

               mulc_spinor_add(vol,pl,pk,z);
               d=(double)(norm_square(vol,icom,pl));
               d=sqrt(d)/(double)(r);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {   
               printf("Consistency of spinor_prod, norm_square,\n");
               printf("and normalize: %.2e\n\n",dmax);
            }
      
            dmax=0.0;
   
            for (i=0;i<5;i++)
            {
               pk=ps[i][off];
               pl=ps[9-i][off];
               random_s(vol,pk,1.0f);
               assign_s2s(vol,pk,pl);
               mulg5(vol,pk);
               mulg5(vol,pk);

               z.re=-1.0f;
               z.im=0.0f;

               mulc_spinor_add(vol,pl,pk,z);
               r=norm_square(vol,icom,pl)/norm_square(vol,icom,pk);
               d=sqrt((double)(r));
               if (d>dmax)
                  dmax=d;

               random_s(vol,pl,1.0f);
               z=spinor_prod(vol,icom,pk,pl);
               mulg5(vol,pk);
               mulg5(vol,pl);
               w=spinor_prod(vol,icom,pk,pl);
      
               d=(fabs((double)(z.re-w.re))+fabs((double)(z.im-w.im)))/
                  (fabs((double)(z.re))+fabs((double)(z.im))); 
               if (d>dmax)
                  dmax=d;

               random_s(vol,pk,1.0f);
               assign_s2s(vol,pk,pl);
               mulg5(vol,pk);
               mulmg5(vol,pk);

               z.re=1.0f;
               z.im=0.0f;

               mulc_spinor_add(vol,pl,pk,z);
               r=norm_square(vol,icom,pl)/norm_square(vol,icom,pk);
               d=sqrt((double)(r));
               if (d>dmax)
                  dmax=d;            
            }

            if (my_rank==0)
               printf("Check of mulg5 and mulmg5: %.2e\n\n",dmax);
         }
      }
   }

   error_chk();
   
   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();   
   exit(0);
}
