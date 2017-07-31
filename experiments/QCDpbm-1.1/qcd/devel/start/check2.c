
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the global field arrays
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,itest;
   int mu,ix,iy,io[4],of[4],n,m;
   su3 *ub,*u;
   su3_dble *udb,*ud;
   spinor *sb;
   spinor_dble *sdb;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   io[0]=      (BNDRY+FACE0)/2;
   io[1]=io[0]+(FACE0+FACE1)/2;
   io[2]=io[1]+(FACE1+FACE2)/2;
   io[3]=io[2]+(FACE2+FACE3)/2;

   of[0]=4*VOLUME;
   of[1]=of[0]+FACE0/2;
   of[2]=of[1]+FACE1/2;
   of[3]=of[2]+FACE2/2;   
   
   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      
      printf("\n");
      printf("Consistency checks on the global field arrays\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   geometry();
   alloc_u();
   alloc_ud();

   itest=0;
   ub=pu[VOLUME/2][0];
   udb=pud[VOLUME/2][0];   

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];
         u=pu[ix][mu];
         ud=pud[ix][mu];

         if ((ix<(VOLUME/2))&&(iy<VOLUME))
         {
            if (u!=(ub+8*(iy-VOLUME/2)+2*mu+1))
               itest=1;
            if (ud!=(udb+8*(iy-VOLUME/2)+2*mu+1))
               itest=2;                
         }
                   
         if (ix>=(VOLUME/2))
         {
            if (u!=(ub+8*(ix-VOLUME/2)+2*mu))
               itest=3;
            if (ud!=(udb+8*(ix-VOLUME/2)+2*mu))
               itest=4;
         }
         
         if ((ix<(VOLUME/2))&&(iy>=VOLUME))
         {
            iy=iy-VOLUME-io[mu];
            iy=of[mu]+iy;

            if (u!=(ub+iy))
               itest=5;
            if (ud!=(udb+iy))
               itest=6;
         }
      }
   }

   error(itest==1,1,
         "main [check2.c]","pu[ix][mu] is incorrect at even points");
   error(itest==2,1,
         "main [check2.c]","pud[ix][mu] is incorrect at even points"); 
   error(itest==3,1,
         "main [check2.c]","pu[ix][mu] is incorrect at odd points");   
   error(itest==4,1,
         "main [check2.c]","pud[ix][mu] is incorrect at odd points");   
   error(itest==5,1,
         "main [check2.c]","pu[ix][mu] is incorrect at even boundary points"); 
   error(itest==6,1,
         "main [check2.c]","pud[ix][mu] is incorrect at even boundary points"); 

   free_u();
   free_ud();

   for (n=1;n<5;n++)
   {
      alloc_s(n);
      alloc_sd(5-n);

      if ((no_s!=n)||(no_sd!=(5-n)))
         itest=1;

      sb=ps[0][0];
      sdb=psd[0][0];

      for (m=0;m<n;m++)
      {
         if (ps[m][0]!=(sb+m*NSPIN))
            itest=2;

         for (ix=0;ix<NSPIN;ix++)
         {
            if (ps[m][ix]!=(ps[m][0]+ix))
               itest=3;
         }
      }

      for (m=0;m<(5-n);m++)
      {
         if (psd[m][0]!=(sdb+m*NSPIN))
            itest=4;

         for (ix=0;ix<NSPIN;ix++)
         {
            if (psd[m][ix]!=(psd[m][0]+ix))
               itest=5;
         }         
      }      
         
      error(itest==1,1,
            "main [check2.c]","no_s or no_sd is not properly set");   
      error(itest==2,1,
            "main [check2.c]","ps[n][0]-ps[0][0] is not equal to NSPIN ");   
      error(itest==3,1,
            "main [check2.c]","ps[n][ix], ix=0,..,NSPIN-1 is not lined up");   
      error(itest==4,1,
            "main [check2.c]","psd[n][0]-psd[0][0] is not equal to NSPIN ");   
      error(itest==5,1,
            "main [check2.c]","psd[n][ix], ix=0,..,NSPIN-1 is not lined up");   

      free_s();
      free_sd();
   }
   
   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }   

   MPI_Finalize(); 
   exit(0);
}
