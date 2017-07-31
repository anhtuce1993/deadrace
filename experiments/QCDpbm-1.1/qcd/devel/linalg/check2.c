
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the programs in the module linalg_dble
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


static complex_dble sp(int vol,spinor_dble *pk,spinor_dble *pl)
{
   int ix;
   double x,y;
   complex_dble *rpk,*rpl,z;

   x=0.0;
   y=0.0;

   rpk=(complex_dble*)(pk);
   rpl=(complex_dble*)(pl);
   
   for (ix=0;ix<(12*vol);ix++)
   {
      x+=((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im);
      y+=((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re);
      rpk+=1;
      rpl+=1;
   }
   
   z.re=x;
   z.im=y;
   
   return z;
}


int main(int argc,char *argv[])
{
   int my_rank,i,vol,off;
   int icom,ieo;
   double r,zsq,d,dmax;
   complex_dble z,w;
   spinor_dble *pk,*pl;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      
      printf("\n");
      printf("Consistency of the programs in the module linalg_dble\n");
      printf("-----------------------------------------------------\n\n");   

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
            alloc_sd(10);

            for (i=0;i<10;i++)
               random_sd(vol,psd[i][off],1.0);

            dmax=0.0;
   
            for (i=0;i<10;i++)
            {
               pk=psd[i][off];
               pl=psd[9-i][off];

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD); 
               }
               else
                  w=sp(vol,pk,pl);

               z=spinor_prod_dble(vol,icom,pk,pl);
               r=norm_square_dble(vol,icom,pk)*norm_square_dble(vol,icom,pl);
               d=(z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im);
               d=sqrt(d/r);
               if (d>dmax)
                  dmax=d;

               z=spinor_prod_dble(vol,icom,pk,pk);
               r=norm_square_dble(vol,icom,pk);
      
               d=fabs(z.im/r);
               if (d>dmax)
                  dmax=d;

               d=fabs(z.re/r-1.0);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {         
               printf("Check of spinor_prod, spinor_prod_re\n");
               printf("and norm_square: %.2e\n\n",dmax);
            }
         
            dmax=0.0;
            z.re= 0.345;
            z.im=-0.876;
            zsq=z.re*z.re+z.im*z.im;
   
            for (i=0;i<9;i++)
            {
               pk=psd[i][off];
               pl=psd[i+1][off];      
      
               w=spinor_prod_dble(vol,icom,pk,pl);
               r=norm_square_dble(vol,icom,pk)+zsq*norm_square_dble(vol,icom,pl)
                  +2.0*(z.re*w.re-z.im*w.im);
               mulc_spinor_add_dble(vol,pk,pl,z);

               d=fabs(r/norm_square_dble(vol,icom,pk)-1.0);
               if (d>dmax)
                  dmax=d;
            }
         
            if (my_rank==0)
            {   
               printf("Consistency of spinor_prod, norm_square\n");
               printf("and mulc_spinor_add: %.2e\n\n",dmax);
            }
         
            for (i=0;i<10;i++)
               random_sd(vol,psd[i][off],1.0);

            dmax=0.0;
   
            for (i=0;i<8;i+=2)
            {
               pk=psd[i][off];
               pl=psd[i+1][off];
               assign_sd2sd(vol,pk,pl);

               r=normalize_dble(vol,icom,pk);
               z=spinor_prod_dble(vol,icom,pk,pl);

               d=fabs(1.0-fabs(z.re)/r)+fabs(z.im)/r;
               if (d>dmax)
                  dmax=d;
               
               z.re=-z.re;
               z.im=-z.im;

               mulc_spinor_add_dble(vol,pl,pk,z);
               d=norm_square_dble(vol,icom,pl);
               d=sqrt(d)/r;
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
               pk=psd[i][off];
               pl=psd[9-i][off];
               random_sd(vol,pk,1.0);
               assign_sd2sd(vol,pk,pl);
               mulg5_dble(vol,pk);
               mulg5_dble(vol,pk);

               z.re=-1.0;
               z.im=0.0;

               mulc_spinor_add_dble(vol,pl,pk,z);
               r=norm_square_dble(vol,icom,pl)/norm_square_dble(vol,icom,pk);
               d=sqrt(r);
               if (d>dmax)
                  dmax=d;

               random_sd(vol,pl,1.0);
               z=spinor_prod_dble(vol,icom,pk,pl);
               mulg5_dble(vol,pk);
               mulg5_dble(vol,pl);
               w=spinor_prod_dble(vol,icom,pk,pl);
      
               d=(fabs(z.re-w.re)+fabs(z.im-w.im))/
                  (fabs(z.re)+fabs(z.im)); 
               if (d>dmax)
                  dmax=d;

               random_sd(vol,pk,1.0);
               assign_sd2sd(vol,pk,pl);
               mulg5_dble(vol,pk);
               mulmg5_dble(vol,pk);

               z.re=1.0;
               z.im=0.0;

               mulc_spinor_add_dble(vol,pl,pk,z);
               r=norm_square_dble(vol,icom,pl)/norm_square_dble(vol,icom,pk);
               d=sqrt(r);
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
