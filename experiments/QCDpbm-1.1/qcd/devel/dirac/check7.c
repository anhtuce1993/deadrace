
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency of Qeo, Qoe, Qnohat and Qhat
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


int main(int argc,char *argv[])
{
   int my_rank,vol;
   float d;
   complex z,w;
   spinor *ps0e,*ps1e,*ps2e,*ps3e;
   spinor *ps0o,*ps1o,*ps2o,*ps3o;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check7.log","w",stdout);
      printf("\n");
      printf("Consistency of Qeo, Qoe, Qnohat and Qhat\n");
      printf("----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,1234);
   geometry();
   alloc_u();
   alloc_ud();
   alloc_s(4);
   alloc_sw();
   alloc_swd();

   random_ud();
   assign_ud2u();

   set_lat_parms(0.0,0.1,0.7);
   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [check7.c]",
         "Inversion of swd on the odd sites was not safe");
   assign_swd2sw();

   ps0e=ps[0][0];
   ps1e=ps[1][0];
   ps2e=ps[2][0];
   ps3e=ps[3][0];

   ps0o=ps[0][VOLUME/2];
   ps1o=ps[1][VOLUME/2];
   ps2o=ps[2][VOLUME/2];
   ps3o=ps[3][VOLUME/2];

   vol=(VOLUME/2);
   z.re=-1.0f;
   z.im=0.0f;

   random_s(2*vol,ps0e,1.0f);
   random_s(2*vol,ps1e,1.0f);
   random_s(2*vol,ps2e,1.0f);

   apply_sw(vol,sw,ps0e,ps2e);
   Qoe(0,1);
   apply_sw(vol,sw+VOLUME,ps1o,ps1o);
   Qeo(1,2);
   Qhat(0,1);
   mulc_spinor_add(vol,ps1e,ps2e,z);
   d=norm_square(vol,1,ps1e)/norm_square(vol,1,ps0e);

   if (my_rank==0)
      printf("|[Qhat-(Qee-Qeo*Qooinv*Qoe)]*psi|/|psi| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(2*vol,ps3e,1.0f);
   assign_s2s(2*vol,ps3e,ps2e);

   sw_term();
   assign_swd2sw();
   apply_sw(2*vol,sw,ps2e,ps2e);
   Qoe(3,1);
   z.re=1.0f;
   mulc_spinor_add(vol,ps2o,ps1o,z);
   set_s2zero(vol,ps1e);
   Qeo(3,1);
   z.re=-1.0f;
   mulc_spinor_add(vol,ps2e,ps1e,z);
   Qnohat(3,0);
   mulc_spinor_add(2*vol,ps2e,ps0e,z);
   d=norm_square(2*vol,1,ps2e)/norm_square(2*vol,1,ps3e);

   if (my_rank==0)
      printf("|[Q-(Qee+Qoo+Qeo+Qoe)]*psi|/|psi| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(2*vol,ps0e,1.0f);
   assign_s2s(2*vol,ps0e,ps1e);
   apply_sw(vol,sw,ps1e,ps1e);
   set_s2zero(vol,ps0o);
   Qnohat(0,2);
   z.re=-1.0f;
   mulc_spinor_add(vol,ps2e,ps1e,z);
   d=norm_square(vol,1,ps2e)/norm_square(vol,1,ps0e);

   if (my_rank==0)
      printf("|{[Q-Qee]*psi_e}_e|/|psi_e| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(vol,ps0e,1.0f);
   set_s2zero(vol,ps0o);
   Qoe(0,2);
   Qnohat(0,3);
   z.re=-1.0f;
   mulc_spinor_add(vol,ps3o,ps2o,z);
   d=norm_square(vol,1,ps3o)/norm_square(vol,1,ps0e);

   if (my_rank==0)
      printf("|{[Q-Qoe]*psi_e}_o|/|psi_e| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(2*vol,ps0e,1.0f);
   set_s2zero(vol,ps0e);
   assign_s2s(vol,ps0o,ps1o);
   apply_sw(vol,sw+VOLUME,ps1o,ps1o);
   Qnohat(0,2);
   z.re=-1.0f;
   mulc_spinor_add(vol,ps2o,ps1o,z);
   d=norm_square(vol,1,ps2o)/norm_square(vol,1,ps0o);

   if (my_rank==0)
      printf("|{[Q-Qoo]*psi_o}_o|/|psi_o| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(2*vol,ps0e,1.0f);
   set_s2zero(vol,ps0e);
   set_s2zero(vol,ps2e);
   Qeo(0,2);
   Qnohat(0,3);
   z.re=1.0f;
   mulc_spinor_add(vol,ps3e,ps2e,z);
   d=norm_square(vol,1,ps3e)/norm_square(vol,1,ps0o);

   if (my_rank==0)
      printf("|{[Q-Qeo]*psi_o}_e|/|psi_o| = %.2e\n\n",
             sqrt((double)(d)));

   random_s(2*vol,ps0e,1.0f);
   random_s(2*vol,ps1e,1.0f);
   random_s(2*vol,ps2e,1.0f);
   random_s(2*vol,ps3e,1.0f);
   set_s2zero(vol,ps2e);
   Qeo(1,2);
   Qoe(0,3);

   z=spinor_prod(vol,1,ps1o,ps3o);
   w=spinor_prod(vol,1,ps2e,ps0e);

   d=(z.re+w.re)*(z.re+w.re)+(z.im+w.im)*(z.im+w.im);
   d/=(float)(norm_square(vol,1,ps1o)*norm_square(vol,1,ps0e));
   error_chk();

   if (my_rank==0)
   {
      printf("|<psi1,Qeo*psi2>-<Qoe*psi1,psi2>|/|psi1||psi2| = %.2e\n\n",
             sqrt((double)(d)));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
