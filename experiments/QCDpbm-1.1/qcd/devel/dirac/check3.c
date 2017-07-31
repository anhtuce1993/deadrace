
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Hermiticity of Qhat
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
#include "sw_term.h"
#include "dirac.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,i;
   float norm,d,r1,r2;
   complex z1,z2;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      printf("\n");
      printf("Hermiticity of Qhat (random fields)\n");
      printf("-----------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,12345);
   geometry();
   alloc_u();
   alloc_ud();
   alloc_s(4);
   alloc_sw();
   alloc_swd();

   random_ud();
   assign_ud2u();

   set_lat_parms(0.0,0.1,0.2);
   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [check3.c]",
         "Inversion of swd on the odd sites was not safe");
   assign_swd2sw();

   for (i=0;i<4;i++)
      random_s(VOLUME,ps[i][0],1.0f);

   Qhat(0,2);
   Qhat(1,3);

   r1=norm_square(VOLUME/2,1,ps[0][0]);
   r2=norm_square(VOLUME/2,1,ps[1][0]);

   norm=(float)(sqrt((double)(r1*r2)));

   z1=spinor_prod(VOLUME/2,1,ps[0][0],ps[3][0]);
   z2=spinor_prod(VOLUME/2,1,ps[2][0],ps[1][0]);

   d=(float)(sqrt((double)((z1.re-z2.re)*(z1.re-z2.re)
                           +(z1.im-z2.im)*(z1.im-z2.im))));
   d/=norm;
   error_chk();

   if (my_rank==0)
   {
      printf("|<r,Qhat*s>-<Qhat*r,s>|/(||r||*||s||) = %4.2e\n",d);
      printf("(should be at most 10^(-7) or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
