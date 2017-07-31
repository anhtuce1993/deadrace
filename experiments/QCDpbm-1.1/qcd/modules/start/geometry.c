
/*******************************************************************************
*
* File geometry.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs related to the lattice geometry
*
* Version suitable for even-odd preconditioning
*
* The externally accessible functions are
*
*   int ipr_global(int n[])
*     This program returns the number of the process with cartesian
*     coordinates n[0],..,n[3] in the process grid
*
*   void ipt_global(int x[],int *ip,int *ix)
*     Given the coordinates x[0],..,x[3] of a point on the full lattice,
*     this program determines the number ip of the process that operates
*     on the corresponding local lattice and the associated local point
*     index ix (0<=ix<VOLUME)
*
*   void geometry(void)
*     Computation of the global arrays cpr,npr describing the process
*     grid and the index arrays ipt,iup,idn and map
*
* Notes:
*
* The programs ipr_global and ipt_global operate locally and can thus be
* conditionally called. All other programs in this directory act globally.
*
*******************************************************************************/

#define GEOMETRY_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "global.h"

static int init=0;
static int cbs[4],cbn[4],*cbix=NULL;


int ipr_global(int n[])
{
   int ip;

   ip=safe_mod(n[0],NPROC0);
   ip=ip*NPROC1+safe_mod(n[1],NPROC1);
   ip=ip*NPROC2+safe_mod(n[2],NPROC2);
   ip=ip*NPROC3+safe_mod(n[3],NPROC3);

   return ip;
}


void ipt_global(int x[],int *ip,int *ix)
{
   int x0,x1,x2,x3;
   int n0,n1,n2,n3;

   x0=safe_mod(x[0],NPROC0*L0);
   x1=safe_mod(x[1],NPROC1*L1);
   x2=safe_mod(x[2],NPROC2*L2);
   x3=safe_mod(x[3],NPROC3*L3);

   n0=x0/L0;
   n1=x1/L1;
   n2=x2/L2;
   n3=x3/L3;

   x0-=(n0*L0);
   x1-=(n1*L1);
   x2-=(n2*L2);
   x3-=(n3*L3);

   *ip=n3+NPROC3*n2+NPROC2*NPROC3*n1+NPROC1*NPROC2*NPROC3*n0;
   *ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
}


static void set_cpr(void)
{
   int np,nr;

   MPI_Comm_size(MPI_COMM_WORLD,&np);

   error(np!=NPROC,1,"set_cpr [geometry.c]",
         "Actual number of processes does not match NPROC");

   MPI_Comm_rank(MPI_COMM_WORLD,&nr);

   error((nr<0)||(nr>=NPROC),1,"set_cpr [geometry.c]",
         "Rank of process is out of range");

   cpr[3]=nr%NPROC3;
   nr/=NPROC3;
   cpr[2]=nr%NPROC2;
   nr/=NPROC2;
   cpr[1]=nr%NPROC1;
   nr/=NPROC1;
   cpr[0]=nr;
}


static void set_npr(void)
{
   int mu,n[4];

   for (mu=0;mu<4;mu++)
      n[mu]=cpr[mu];

   for (mu=0;mu<4;mu++)
   {
      n[mu]-=1;
      npr[2*mu]=ipr_global(n);
      n[mu]+=2;
      npr[2*mu+1]=ipr_global(n);
      n[mu]-=1;
   }
}


static void cache_block(int bs[])
{
   int mu;

   cbs[0]=bs[0];
   cbn[0]=1;

   for (mu=1;mu<4;mu++)
   {
      if ((bs[mu]%4)==0)
         cbs[mu]=4;
      else if ((bs[mu]%3)==0)
         cbs[mu]=3;
      else if ((bs[mu]%2)==0)
         cbs[mu]=2;
      else
         cbs[mu]=1;

      cbn[mu]=bs[mu]/cbs[mu];
   }

   if (cbix!=NULL)
      afree(cbix);

   cbix=amalloc(cbs[0]*cbs[1]*cbs[2]*cbs[3]*sizeof(int),3);
   error(cbix==NULL,1,"cache_block [geometry.c]",
         "Unable to allocate auxiliary array");
}


static void set_cbix(void)
{
   int x0,x1,x2,x3;
   int ig,iu,ib,is;

   ig=0;
   iu=0;

   for (x0=0;x0<cbs[0];x0++)
   {
      for (x1=0;x1<cbs[1];x1++)
      {
         for (x2=0;x2<cbs[2];x2++)
         {
            for (x3=0;x3<cbs[3];x3++)
            {
               ib=x3+cbs[3]*x2+cbs[2]*cbs[3]*x1+cbs[1]*cbs[2]*cbs[3]*x0;
               is=x0+x1+x2+x3;

               if ((is%2)==0)
               {
                  cbix[ib]=ig;
                  ig+=1;
               }
               else
               {
                  cbix[ib]=iu;
                  iu+=1;
               }
            }
         }
      }
   }
}


static int index(int bo[],int bs[],int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int xb1,xb2,xb3;
   int xn1,xn2,xn3;
   int ib,in,is;

   y0=safe_mod(x0,bs[0]);
   y1=safe_mod(x1,bs[1]);
   y2=safe_mod(x2,bs[2]);
   y3=safe_mod(x3,bs[3]);

   xb1=y1%cbs[1];
   xb2=y2%cbs[2];
   xb3=y3%cbs[3];

   xn1=y1/cbs[1];
   xn2=y2/cbs[2];
   xn3=y3/cbs[3];

   ib=cbix[xb3+cbs[3]*xb2+cbs[2]*cbs[3]*xb1+cbs[1]*cbs[2]*cbs[3]*y0];
   in=xn3+cbn[3]*xn2+cbn[3]*cbn[2]*xn1;
   is=y0+y1+y2+y3;
   is+=(bo[0]+bo[1]+bo[2]+bo[3]);

   if (is%2!=0)
      ib+=((bs[0]*bs[1]*bs[2]*bs[3])/2);

   return ib+(cbs[0]*cbs[1]*cbs[2]*cbs[3]*in)/2;
}


void geometry(void)
{
   int x0,x1,x2,x3;
   int k,mu,ix,iy,iz,iw;
   int bo[4],bs[4],ifc[8];

   set_cpr();
   set_npr();

   bo[0]=0;
   bo[1]=0;
   bo[2]=0;
   bo[3]=0;

   bs[0]=L0;
   bs[1]=L1;
   bs[2]=L2;
   bs[3]=L3;

   cache_block(bs);
   set_cbix();

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=index(bo,bs,x0,x1,x2,x3);
               iy=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0;
               ipt[iy]=ix;

               iup[ix][0]=index(bo,bs,x0+1,x1,x2,x3);
               idn[ix][0]=index(bo,bs,x0-1,x1,x2,x3);

               iup[ix][1]=index(bo,bs,x0,x1+1,x2,x3);
               idn[ix][1]=index(bo,bs,x0,x1-1,x2,x3);

               iup[ix][2]=index(bo,bs,x0,x1,x2+1,x3);
               idn[ix][2]=index(bo,bs,x0,x1,x2-1,x3);

               iup[ix][3]=index(bo,bs,x0,x1,x2,x3+1);
               idn[ix][3]=index(bo,bs,x0,x1,x2,x3-1);

               if ((x0==(L0-1))&&(NPROC0>1))
                  iup[ix][0]=VOLUME;
               if ((x0==0)&&(NPROC0>1))
                  idn[ix][0]=VOLUME;

               if ((x1==(L1-1))&&(NPROC1>1))
                  iup[ix][1]=VOLUME;
               if ((x1==0)&&(NPROC1>1))
                  idn[ix][1]=VOLUME;

               if ((x2==(L2-1))&&(NPROC2>1))
                  iup[ix][2]=VOLUME;
               if ((x2==0)&&(NPROC2>1))
                  idn[ix][2]=VOLUME;

               if ((x3==(L3-1))&&(NPROC3>1))
                  iup[ix][3]=VOLUME;
               if ((x3==0)&&(NPROC3>1))
                  idn[ix][3]=VOLUME;
            }
         }
      }
   }

   ifc[0]=0;
   ifc[1]=ifc[0]+(FACE0/2);
   ifc[2]=ifc[1]+(FACE0/2);
   ifc[3]=ifc[2]+(FACE1/2);
   ifc[4]=ifc[3]+(FACE1/2);
   ifc[5]=ifc[4]+(FACE2/2);
   ifc[6]=ifc[5]+(FACE2/2);
   ifc[7]=ifc[6]+(FACE3/2);

   for (ix=0;ix<VOLUME;ix++)
   {
      if (ix==(VOLUME/2))
      {
         ifc[0]=(BNDRY/2);
         ifc[1]=ifc[0]+(FACE0/2);
         ifc[2]=ifc[1]+(FACE0/2);
         ifc[3]=ifc[2]+(FACE1/2);
         ifc[4]=ifc[3]+(FACE1/2);
         ifc[5]=ifc[4]+(FACE2/2);
         ifc[6]=ifc[5]+(FACE2/2);
         ifc[7]=ifc[6]+(FACE3/2);
      }

      iy=(ix+(VOLUME/2))%VOLUME;

      for (mu=0;mu<4;mu++)
      {
         if (idn[iy][mu]==VOLUME)
         {
            iz=ifc[2*mu];
            ifc[2*mu]+=1;

            idn[iy][mu]=VOLUME+iz;
            iw=iy;

            for (k=1;k<bs[mu];k++)
               iw=iup[iw][mu];

            map[iz]=iw;
         }

         if (iup[iy][mu]==VOLUME)
         {
            iz=ifc[2*mu+1];
            ifc[2*mu+1]+=1;

            iup[iy][mu]=VOLUME+iz;
            iw=iy;

            for (k=1;k<bs[mu];k++)
               iw=idn[iw][mu];

            map[iz]=iw;
         }
      }
   }

   init=1;
}

