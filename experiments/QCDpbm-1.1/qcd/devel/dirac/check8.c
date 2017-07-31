
/*******************************************************************************
*
* File check8.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency of Qeo_dble, Qoe_dble, Qnohat_dble and Qhat_dble
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
   double d;
   complex_dble z,w;
   spinor_dble *ps0e,*ps1e,*ps2e,*ps3e;
   spinor_dble *ps0o,*ps1o,*ps2o,*ps3o;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check8.log","w",stdout);
      printf("\n");
      printf("Consistency of Qeo_dble, Qoe_dble, Qnohat_dble and Qhat_dble\n");
      printf("------------------------------------------------------------\n");
      printf("\n");
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,1234);
   geometry();
   alloc_ud();
   alloc_sd(4);
   alloc_swd();

   random_ud();
   set_lat_parms(0.0,0.1,0.3);
   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [check8.c]",
         "Inversion of swd on the odd sites was not safe");

   ps0e=psd[0][0];
   ps1e=psd[1][0];
   ps2e=psd[2][0];
   ps3e=psd[3][0];

   ps0o=psd[0][VOLUME/2];
   ps1o=psd[1][VOLUME/2];
   ps2o=psd[2][VOLUME/2];
   ps3o=psd[3][VOLUME/2];

   vol=(VOLUME/2);
   z.re=-1.0;
   z.im=0.0;

   random_sd(2*vol,ps0e,1.0);
   random_sd(2*vol,ps1e,1.0);
   random_sd(2*vol,ps2e,1.0);

   apply_sw_dble(vol,swd,ps0e,ps2e);
   Qoe_dble(0,1);
   apply_sw_dble(vol,swd+VOLUME,ps1o,ps1o);
   Qeo_dble(1,2);
   Qhat_dble(0,1);
   mulc_spinor_add_dble(vol,ps1e,ps2e,z);
   d=norm_square_dble(vol,1,ps1e)/norm_square_dble(vol,1,ps0e);

   if (my_rank==0)
      printf("|[Qhat-(Qee-Qeo*Qooinv*Qoe)]*psi|/|psi| = %.2e\n\n",
             sqrt(d));

   random_sd(2*vol,ps3e,1.0);
   assign_sd2sd(2*vol,ps3e,ps2e);

   sw_term();
   apply_sw_dble(2*vol,swd,ps2e,ps2e);
   Qoe_dble(3,1);
   z.re=1.0;
   mulc_spinor_add_dble(vol,ps2o,ps1o,z);
   set_sd2zero(vol,ps1e);
   Qeo_dble(3,1);
   z.re=-1.0;
   mulc_spinor_add_dble(vol,ps2e,ps1e,z);
   Qnohat_dble(3,0);
   mulc_spinor_add_dble(2*vol,ps2e,ps0e,z);
   d=norm_square_dble(2*vol,1,ps2e)/norm_square_dble(2*vol,1,ps3e);

   if (my_rank==0)
      printf("|[Q-(Qee+Qoo+Qeo+Qoe)]*psi|/|psi| = %.2e\n\n",
             sqrt(d));

   random_sd(2*vol,ps0e,1.0);
   assign_sd2sd(2*vol,ps0e,ps1e);
   apply_sw_dble(vol,swd,ps1e,ps1e);
   set_sd2zero(vol,ps0o);
   Qnohat_dble(0,2);
   z.re=-1.0;
   mulc_spinor_add_dble(vol,ps2e,ps1e,z);
   d=norm_square_dble(vol,1,ps2e)/norm_square_dble(vol,1,ps0e);

   if (my_rank==0)
      printf("|{[Q-Qee]*psi_e}_e|/|psi_e| = %.2e\n\n",
             sqrt(d));

   random_sd(vol,ps0e,1.0);
   set_sd2zero(vol,ps0o);
   Qoe_dble(0,2);
   Qnohat_dble(0,3);
   z.re=-1.0;
   mulc_spinor_add_dble(vol,ps3o,ps2o,z);
   d=norm_square_dble(vol,1,ps3o)/norm_square_dble(vol,1,ps0e);

   if (my_rank==0)
      printf("|{[Q-Qoe]*psi_e}_o|/|psi_e| = %.2e\n\n",
             sqrt(d));

   random_sd(2*vol,ps0e,1.0);
   set_sd2zero(vol,ps0e);
   assign_sd2sd(vol,ps0o,ps1o);
   apply_sw_dble(vol,swd+VOLUME,ps1o,ps1o);
   Qnohat_dble(0,2);
   z.re=-1.0;
   mulc_spinor_add_dble(vol,ps2o,ps1o,z);
   d=norm_square_dble(vol,1,ps2o)/norm_square_dble(vol,1,ps0o);

   if (my_rank==0)
      printf("|{[Q-Qoo]*psi_o}_o|/|psi_o| = %.2e\n\n",
             sqrt(d));

   random_sd(2*vol,ps0e,1.0);
   set_sd2zero(vol,ps0e);
   set_sd2zero(vol,ps2e);
   Qeo_dble(0,2);
   Qnohat_dble(0,3);
   z.re=1.0;
   mulc_spinor_add_dble(vol,ps3e,ps2e,z);
   d=norm_square_dble(vol,1,ps3e)/norm_square_dble(vol,1,ps0o);

   if (my_rank==0)
      printf("|{[Q-Qeo]*psi_o}_e|/|psi_o| = %.2e\n\n",
             sqrt(d));

   random_sd(2*vol,ps0e,1.0);
   random_sd(2*vol,ps1e,1.0);
   random_sd(2*vol,ps2e,1.0);
   random_sd(2*vol,ps3e,1.0);
   set_sd2zero(vol,ps2e);
   Qeo_dble(1,2);
   Qoe_dble(0,3);

   z=spinor_prod_dble(vol,1,ps1o,ps3o);
   w=spinor_prod_dble(vol,1,ps2e,ps0e);

   d=(z.re+w.re)*(z.re+w.re)+(z.im+w.im)*(z.im+w.im);
   d/=(norm_square_dble(vol,1,ps1o)*norm_square_dble(vol,1,ps0e));
   error_chk();

   if (my_rank==0)
   {
      printf("|<psi1,Qeo*psi2>-<Qoe*psi1,psi2>|/|psi1||psi2| = %.2e\n\n",
             sqrt(d));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
