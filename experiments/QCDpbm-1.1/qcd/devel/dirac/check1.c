
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance of Qhat
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

static int nfc[8],ofs[8];
static su3 *g,*gbuf;


static void pack_gbuf(void)
{
   int n,ix,iy,io;

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (n=0;n<8;n++)
   {
      io=ofs[n];

      for (ix=0;ix<nfc[n];ix++)
      {
         iy=map[io+ix];
         gbuf[io+ix]=g[iy];
      }
   }
}


static void send_gbuf(void)
{
   int n,mu,np,saddr,raddr;
   int nbf,tag;
   float *sbuf,*rbuf;
   MPI_Status stat;

   for (n=0;n<8;n++)
   {
      nbf=18*nfc[n];

      if (nbf>0)
      {
         tag=mpi_tag();
         mu=n/2;
         np=cpr[mu];

         if (n==(2*mu))
         {
            saddr=npr[n+1];
            raddr=npr[n];
         }
         else
         {
            saddr=npr[n-1];
            raddr=npr[n];
         }

         sbuf=(float*)(gbuf+ofs[n]);
         rbuf=(float*)(g+ofs[n]+VOLUME);

         if ((np|0x1)!=np)
         {
            MPI_Send(sbuf,nbf,MPI_FLOAT,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_FLOAT,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_FLOAT,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_FLOAT,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void random_g(void)
{
   int ix;

   for (ix=0;ix<VOLUME;ix++)
      random_su3(g+ix);

   if (BNDRY>0)
   {
      pack_gbuf();
      send_gbuf();
   }
}


static void transform_u(void)
{
   int ix,iy,mu;
   su3 *ub,u,v,w,gx,gxi,gy,gyi;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      ub=pu[ix][0];
      gx=g[ix];

      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];
         gy=g[iy];
         u=ub[2*mu];
         _su3_dagger(gyi,gy);
         _su3_times_su3(v,u,gyi);
         _su3_times_su3(w,gx,v);
         ub[2*mu]=w;

         iy=idn[ix][mu];
         gy=g[iy];
         u=ub[2*mu+1];
         _su3_dagger(gxi,gx);
         _su3_times_su3(v,u,gxi);
         _su3_times_su3(w,gy,v);
         ub[2*mu+1]=w;
      }
   }
}


static void transform_s(spinor *pk,spinor *pl)
{
   int ix;
   su3 gx;
   spinor r,s;

   for (ix=0;ix<VOLUME;ix++)
   {
      s=pk[ix];
      gx=g[ix];

      _su3_multiply(r.c1,gx,s.c1);
      _su3_multiply(r.c2,gx,s.c2);
      _su3_multiply(r.c3,gx,s.c3);
      _su3_multiply(r.c4,gx,s.c4);

      pl[ix]=r;
   }
}


int main(int argc,char *argv[])
{
   int my_rank,i,vol;
   float d;
   complex z;
   spinor *ps0,*ps1,*ps2,*ps3;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      printf("\n");
      printf("Gauge covariance of Qhat (random fields)\n");
      printf("----------------------------------------\n\n");

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

   g=amalloc(NSPIN*sizeof(su3),4);

   if (BNDRY>0)
      gbuf=amalloc((BNDRY/2)*sizeof(su3),4);

   error((g==NULL)||((BNDRY>0)&&(gbuf==NULL)),1,"main [check1.c]",
         "Unable to allocate auxiliary arrays");
   random_g();

   ps0=ps[0][0];
   ps1=ps[1][0];
   ps2=ps[2][0];
   ps3=ps[3][0];

   vol=(VOLUME/2);
   set_lat_parms(0.0,0.1,0.6);
   z.re=-1.0f;
   z.im=0.0f;

   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [check1.c]",
         "Inversion of swd on the odd sites was not safe");
   assign_swd2sw();

   for (i=0;i<4;i++)
      random_s(vol,ps[i][0],1.0f);

   transform_s(ps0,ps2);
   Qhat(0,1);
   transform_u();
   assign_u2ud();
   sw_term();
   error(invert_swd(ODD_PTS)!=0,1,"main [check1.c]",
         "Inversion of swd on the odd sites was not safe");
   assign_swd2sw();
   Qhat(2,3);
   transform_s(ps1,ps2);

   mulc_spinor_add(vol,ps3,ps2,z);
   d=norm_square(vol,1,ps3)/norm_square(vol,1,ps0);
   error_chk();

   if (my_rank==0)
   {
      printf("Maximal normalized difference = %.2e\n",sqrt(d));
      printf("(should be around 1*10^(-6) or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
