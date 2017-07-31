
/*******************************************************************************
*
* File check5n.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Action of Qnohat_dble on plane waves
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


static su3_vector_dble mul_cplx(complex_dble z,su3_vector_dble s)
{
   su3_vector_dble r;

   r.c1.re=z.re*s.c1.re-z.im*s.c1.im;
   r.c1.im=z.im*s.c1.re+z.re*s.c1.im;
   r.c2.re=z.re*s.c2.re-z.im*s.c2.im;
   r.c2.im=z.im*s.c2.re+z.re*s.c2.im;
   r.c3.re=z.re*s.c3.re-z.im*s.c3.im;
   r.c3.im=z.im*s.c3.re+z.re*s.c3.im;

   return r;
}


static spinor_dble mul_gamma(int mu,spinor_dble s)
{
   spinor_dble r;
   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      r.c1=mul_cplx(m_1,s.c3);
      r.c2=mul_cplx(m_1,s.c4);
      r.c3=mul_cplx(m_1,s.c1);
      r.c4=mul_cplx(m_1,s.c2);
   }
   else if (mu==1)
   {
      r.c1=mul_cplx(m_i,s.c4);
      r.c2=mul_cplx(m_i,s.c3);
      r.c3=mul_cplx(i,s.c2);
      r.c4=mul_cplx(i,s.c1);
   }
   else if (mu==2)
   {
      r.c1=mul_cplx(m_1,s.c4);
      r.c2=s.c3;
      r.c3=s.c2;
      r.c4=mul_cplx(m_1,s.c1);
   }
   else if (mu==3)
   {
      r.c1=mul_cplx(m_i,s.c3);
      r.c2=mul_cplx(i,s.c4);
      r.c3=mul_cplx(i,s.c1);
      r.c4=mul_cplx(m_i,s.c2);
   }
   else
   {
      r.c1=s.c1;
      r.c2=s.c2;
      r.c3=mul_cplx(m_1,s.c3);
      r.c4=mul_cplx(m_1,s.c4);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int my_rank;
   int n,i,j,ix,mu,vol;
   int x0,x1,x2,x3;
   int np[4],bo[4];
   float ran[4];
   double mp,d,pi,p[4],sp[4],px,r,*rs;
   complex_dble z;
   spinor_dble s,s0,s1,s2,s3,*ps0,*ps1,*ps2;
   lat_parms_t lat;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5n.log","w",stdout);
      printf("\n");
      printf("Action of Qnohat_dble on plane waves\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("For this test to pass, the calculated differences delta\n");
      printf("should be at most 1*10^(-14) or so\n\n");
   }

   start_ranlux(0,12345);
   geometry();
   alloc_ud();
   alloc_sd(3);
   alloc_swd();

   vol=VOLUME;
   pi=4.0*atan(1.0);
   lat=set_lat_parms(0.0,0.1,1.0);
   rs=(double*)(&s);
   n=10;

   sw_term();

   ps0=psd[0][0];
   ps1=psd[1][0];
   ps2=psd[2][0];

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   for (i=0;i<n;i++)
   {
      ranlxs(ran,4);

      np[0]=(int)(ran[0]*(float)(NPROC0*L0));
      np[1]=(int)(ran[1]*(float)(NPROC1*L1));
      np[2]=(int)(ran[2]*(float)(NPROC2*L2));
      np[3]=(int)(ran[3]*(float)(NPROC3*L3));

      p[0]=(double)(np[0])*2.0*pi/(double)(NPROC0*L0);
      p[1]=(double)(np[1])*2.0*pi/(double)(NPROC1*L1);
      p[2]=(double)(np[2])*2.0*pi/(double)(NPROC2*L2);
      p[3]=(double)(np[3])*2.0*pi/(double)(NPROC3*L3);

      sp[0]=sin(p[0]);
      sp[1]=sin(p[1]);
      sp[2]=sin(p[2]);
      sp[3]=sin(p[3]);

      mp=lat.m0;
      mp+=(1.0-cos(p[0]));
      mp+=(1.0-cos(p[1]));
      mp+=(1.0-cos(p[2]));
      mp+=(1.0-cos(p[3]));

      r=0.0;

      while ((1.0+r)==1.0)
      {
         gauss_dble(rs,24);
         r=0.0;

         for (j=0;j<24;j++)
            r+=rs[j]*rs[j];
         r=sqrt(r);
      }

      MPI_Bcast(p,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(sp,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&mp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(rs,24,MPI_DOUBLE,0,MPI_COMM_WORLD);

      for (x0=0;x0<L0;x0++)
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

                  px= p[0]*(double)(x0+bo[0])+p[1]*(double)(x1+bo[1])
                     +p[2]*(double)(x2+bo[2])+p[3]*(double)(x3+bo[3]);

                  z.re=cos(px);
                  z.im=sin(px);

                  s0.c1=mul_cplx(z,s.c1);
                  s0.c2=mul_cplx(z,s.c2);
                  s0.c3=mul_cplx(z,s.c3);
                  s0.c4=mul_cplx(z,s.c4);

                  *(psd[0][ix])=s0;

                  z.re=mp;
                  z.im=0.0;

                  s1.c1=mul_cplx(z,s0.c1);
                  s1.c2=mul_cplx(z,s0.c2);
                  s1.c3=mul_cplx(z,s0.c3);
                  s1.c4=mul_cplx(z,s0.c4);

                  for (mu=0;mu<4;mu++)
                  {
                     s2=mul_gamma(mu,s0);
                     z.re=0.0;
                     z.im=sp[mu];
                     s3.c1=mul_cplx(z,s2.c1);
                     s3.c2=mul_cplx(z,s2.c2);
                     s3.c3=mul_cplx(z,s2.c3);
                     s3.c4=mul_cplx(z,s2.c4);

                     _vector_add_assign(s1.c1,s3.c1);
                     _vector_add_assign(s1.c2,s3.c2);
                     _vector_add_assign(s1.c3,s3.c3);
                     _vector_add_assign(s1.c4,s3.c4);
                  }

                  *(psd[1][ix])=mul_gamma(5,s1);
               }
            }
         }
      }

      Qnohat_dble(0,2);

      z.re=-1.0;
      z.im=0.0;
      mulc_spinor_add_dble(vol,ps2,ps1,z);
      d=norm_square_dble(vol,1,ps2)/norm_square_dble(vol,1,ps0);

      if (my_rank==0)
         printf("delta = %4.2e at p=(%d,%d,%d,%d)\n",
                sqrt(d),np[0],np[1],np[2],np[3]);
   }

   error_chk();

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
