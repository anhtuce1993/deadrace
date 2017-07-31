
/*******************************************************************************
*
* File ucom.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication of the link variables residing at the boundaries of the local
* lattices
*
* The externally accessible functions are
*
*   void copy_bnd_u(void)
*     Fetches the single-precision link variables on the boundaries of the
*     local lattice from the neighbouring processes and stores them in a
*     buffer
*
*   void copy_bnd_ud(void)
*     Fetches the double-precision link variables on the boundaries of the
*     local lattice from the neighbouring processes and stores them in a
*     buffer
*
*   void plaq_u(int ix,int mu,int nu,su3 **u)
*     Calculates the pointers u[4] to the four single-precision link
*     variables in the (mu,nu)-plaquette at the point ix on the local
*     lattice. The values stored at these memory locations are correct
*     only after copy_bnd_u() is called
*
*   void plaq_ud(int ix,int mu,int nu,su3_dble **u)
*     Calculates the pointers u[4] to the four double-precision link
*     variables in the (mu,nu)-plaquette at the point ix on the local
*     lattice. The values stored at these memory locations are correct
*     only after copy_bnd_ud() is called
*
*   void free_ucom_bufs(int iu,int iud)
*     Frees the communication buffers used by copy_bnd_u if iu=1 and those
*     used by copy_bnd_ud if iud=1
*
* Notes:
*
* After calling copy_bnd_u() all link variables pointed to by the global
* pointer array pu[][] will have the correct values. Whether the values
* are up-to-date can always be checked by querying the flags data base (see
* flags.c). The same comments also apply to the double-precision fields and
* programs
*
* Internally the copy programs allocate buffers to store the link
* variables on the +0,+1,+2,+3 faces of the local lattice. These can be
* freed using free_ucom_bufs(), but if the programs are called frequently,
* it is better to leave the buffers allocated (in which case they will be
* reused)
*
* The link variables in the (mu,nu)-plaquette at the point x are ordered
* according to
*
*   u[0] -> U(x,mu)
*   u[1] -> U(x+mu,nu)
*   u[2] -> U(x,nu)
*   u[3] -> U(x+nu,mu)
*
* In the programs plaq_u() and plaq_ud() it is taken for granted that the
* arguments satisfy mu!=nu and 0<=ix<VOLUME
*
*******************************************************************************/

#define UCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "global.h"

typedef struct
{
   int nu0,nuk;
   int saddr,raddr;
   int iu0,iuk;
} comdat_t;

static int nfc[4],ofs[4],ofs_uke[4],ofs_uko[4];
static int *idx_u0=NULL,*idx_uk;
static su3 *sbuf_u0=NULL,*sbuf_uk,*rbuf_u0,*rbuf_uk;
static su3_dble *sdbuf_u0=NULL,*sdbuf_uk,*rdbuf_u0,*rdbuf_uk;
static comdat_t comdat[4];


static void alloc_idx(void)
{
   int n,mu,iu0,iuk;
   comdat_t *c;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   ofs[0]=FACE0/2;
   ofs[1]=ofs[0]+(FACE0+FACE1)/2;
   ofs[2]=ofs[1]+(FACE1+FACE2)/2;
   ofs[3]=ofs[2]+(FACE2+FACE3)/2;

   n=BNDRY/4;
   idx_u0=amalloc(7*n*sizeof(int),3);
   error(idx_u0==NULL,1,"alloc_idx [ucom.c]","Unable to allocate index array");

   idx_uk=idx_u0+n;
   iu0=0;
   iuk=0;

   for (mu=0;mu<4;mu++)
   {
      c=comdat+mu;

      (*c).nu0=nfc[mu]/2;
      (*c).nuk=3*nfc[mu];
      (*c).saddr=npr[2*mu];
      (*c).raddr=npr[2*mu+1];
      (*c).iu0=iu0;
      (*c).iuk=iuk;

      ofs_uke[mu]=iuk-3*(VOLUME+ofs[mu]);
      ofs_uko[mu]=iuk+3*(nfc[mu]/2)-3*(VOLUME+(BNDRY/2)+ofs[mu]);

      iu0+=(*c).nu0;
      iuk+=(*c).nuk;
   }
}


static void set_idx(void)
{
   int mu,nu;
   int ioe,ioo,ix,iye,iyo,izo;
   int nu0,*u0,*uke,*uko;

   u0=idx_u0;
   uke=idx_uk;

   for (mu=0;mu<4;mu++)
   {
      nu0=nfc[mu]/2;
      ioe=ofs[mu];
      ioo=ioe+(BNDRY/2);
      uko=uke+3*nu0;

      for (ix=0;ix<nu0;ix++)
      {
         iye=map[ioe+ix];
         iyo=map[ioo+ix];

         *u0=8*(iyo-(VOLUME/2))+2*mu+1;
         u0+=1;

         for (nu=0;nu<4;nu++)
         {
            if (nu!=mu)
            {
               izo=iup[iye][nu];

               if (izo<VOLUME)
                  *uke=8*(izo-(VOLUME/2))+2*nu+1;
               else
                  *uke=4*VOLUME+comdat[nu].iu0+(izo-VOLUME-(BNDRY/2)-ofs[nu]);

               *uko=8*(iyo-(VOLUME/2))+2*nu;

               uke+=1;
               uko+=1;
            }
         }
      }

      uke=uko;
   }
}


static void alloc_ubufs(void)
{
   int n;

   error(pu[VOLUME/2][0]==NULL,1,"alloc_ubufs [ucom.c]",
         "Single-precision gauge field is not allocated");

   if (idx_u0==NULL)
   {
      alloc_idx();
      set_idx();
   }

   n=BNDRY/4;
   sbuf_u0=amalloc(13*n*sizeof(su3),ALIGN);
   error(sbuf_u0==NULL,1,"alloc_ubufs [ucom.c]","Unable to allocate buffers");

   sbuf_uk=sbuf_u0+n;
   rbuf_u0=pu[VOLUME/2][0]+4*VOLUME;
   rbuf_uk=sbuf_uk+6*n;
}


static void pack_u0(void)
{
   int *iu,*ium;
   su3 *u,*ub;

   u=sbuf_u0;
   ub=pu[VOLUME/2][0];
   iu=idx_u0;
   ium=idx_u0+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void pack_uk(void)
{
   int *iu,*ium;
   su3 *u,*ub;

   u=sbuf_uk;
   ub=pu[VOLUME/2][0];
   iu=idx_uk;
   ium=idx_uk+3*(BNDRY/2);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void send_u0(void)
{
   int tag,nbf,saddr,raddr,np;
   float *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nu0;

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=(*c).saddr;
         raddr=(*c).raddr;
         sbuf=(float*)(sbuf_u0+(*c).iu0);
         rbuf=(float*)(rbuf_u0+(*c).iu0);

         if (np==0)
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


static void send_uk(void)
{
   int tag,nbf,saddr,raddr,np;
   float *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nuk;

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=(*c).saddr;
         raddr=(*c).raddr;
         sbuf=(float*)(sbuf_uk+(*c).iuk);
         rbuf=(float*)(rbuf_uk+(*c).iuk);

         if (np==0)
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


void copy_bnd_u(void)
{
   if (NPROC>1)
   {
      if (sbuf_u0==NULL)
         alloc_ubufs();

      pack_u0();
      send_u0();
      pack_uk();
      send_uk();
   }
}


static void alloc_udbufs(void)
{
   int n;

   error(pud[VOLUME/2][0]==NULL,1,"alloc_udbufs [ucom.c]",
         "Double-precision gauge field is not allocated");

   if (idx_u0==NULL)
   {
      alloc_idx();
      set_idx();
   }

   n=BNDRY/4;
   sdbuf_u0=amalloc(13*n*sizeof(su3_dble),ALIGN);
   error(sdbuf_u0==NULL,1,"alloc_udbufs [ucom.c]","Unable to allocate buffers");

   sdbuf_uk=sdbuf_u0+n;
   rdbuf_u0=pud[VOLUME/2][0]+4*VOLUME;
   rdbuf_uk=sdbuf_uk+6*n;
}


static void pack_ud0(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_u0;
   ub=pud[VOLUME/2][0];
   iu=idx_u0;
   ium=idx_u0+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void pack_udk(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_uk;
   ub=pud[VOLUME/2][0];
   iu=idx_uk;
   ium=idx_uk+3*(BNDRY/2);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void send_ud0(void)
{
   int tag,nbf,saddr,raddr,np;
   double *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nu0;

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=(*c).saddr;
         raddr=(*c).raddr;
         sbuf=(double*)(sdbuf_u0+(*c).iu0);
         rbuf=(double*)(rdbuf_u0+(*c).iu0);

         if (np==0)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void send_udk(void)
{
   int tag,nbf,saddr,raddr,np;
   double *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nuk;

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=(*c).saddr;
         raddr=(*c).raddr;
         sbuf=(double*)(sdbuf_uk+(*c).iuk);
         rbuf=(double*)(rdbuf_uk+(*c).iuk);

         if (np==0)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


void copy_bnd_ud(void)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_udbufs();

      pack_ud0();
      send_ud0();
      pack_udk();
      send_udk();
   }
}


void plaq_u(int ix,int mu,int nu,su3 **u)
{
   int iy;

   u[0]=pu[ix][mu];
   iy=iup[ix][mu];

   if (iy<VOLUME)
      u[1]=pu[iy][nu];
   else
   {
      if (iy<(VOLUME+(BNDRY/2)))
         u[1]=rbuf_uk+ofs_uke[mu];
      else
         u[1]=rbuf_uk+ofs_uko[mu];

      u[1]+=(3*iy+nu-(nu>mu));
   }

   u[2]=pu[ix][nu];
   iy=iup[ix][nu];

   if (iy<VOLUME)
      u[3]=pu[iy][mu];
   else
   {
      if (iy<(VOLUME+(BNDRY/2)))
         u[3]=rbuf_uk+ofs_uke[nu];
      else
         u[3]=rbuf_uk+ofs_uko[nu];

      u[3]+=(3*iy+mu-(mu>nu));
   }
}


void plaq_ud(int ix,int mu,int nu,su3_dble **u)
{
   int iy;

   u[0]=pud[ix][mu];
   iy=iup[ix][mu];

   if (iy<VOLUME)
      u[1]=pud[iy][nu];
   else
   {
      if (iy<(VOLUME+(BNDRY/2)))
         u[1]=rdbuf_uk+ofs_uke[mu];
      else
         u[1]=rdbuf_uk+ofs_uko[mu];

      u[1]+=(3*iy+nu-(nu>mu));
   }

   u[2]=pud[ix][nu];
   iy=iup[ix][nu];

   if (iy<VOLUME)
      u[3]=pud[iy][mu];
   else
   {
      if (iy<(VOLUME+(BNDRY/2)))
         u[3]=rdbuf_uk+ofs_uke[nu];
      else
         u[3]=rdbuf_uk+ofs_uko[nu];

      u[3]+=(3*iy+mu-(mu>nu));
   }
}


void free_ucom_bufs(int iu,int iud)
{
   if ((iu==1)&&(sbuf_u0!=NULL))
   {
      afree(sbuf_u0);
      sbuf_u0=NULL;
   }

   if ((iud==1)&&(sdbuf_u0!=NULL))
   {
      afree(sdbuf_u0);
      sdbuf_u0=NULL;
   }

   if ((sbuf_u0==NULL)&&(sdbuf_u0==NULL))
   {
      afree(idx_u0);
      idx_u0=NULL;
   }
}
