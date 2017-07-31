
/*******************************************************************************
*
* File scom.c
*
* Copyright (C) 2005, 2008 Martin Luescher, 2007 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication functions related to Qhat
*
* Spinor fields are defined on the local lattices and the *even* exterior
* boundary points of these. Copying spinors to or from the boundary thus
* refers to these points only
*
*   void cps_int_bnd(int k)
*     Copies the spinors *ps[k][] at the interior boundary points of the 
*     local lattice to the corresponding points on the neighbouring processes.
*     This program compresses the Dirac spinors psi to Weyl spinors before
*     copying and expands the latter after copying to spinors chi such that
*     theta*chi=theta*psi where theta denotes the projector to the exterior
*     boundary of the target lattice
*
*   void cps_ext_bnd(int k)
*     Copies the spinors *ps[k][] at the exterior boundary points of the
*     local lattice to the neighbouring processes and *adds* them to the 
*     field on the matching points of the target lattices. This program
*     compresses the Dirac spinors psi to Weyl spinors before copying,
*     assuming that theta*psi=psi where theta denotes the projector to 
*     the exterior boundary of the local lattice
*
*   void free_sbufs(void)
*     Frees the communication buffers used by the programs in this module
* 
* Notes:
*
* The copy programs allocate the required communication buffers when
* needed. They can be freed using free_sbufs(), but if the programs are
* called frequently, it is better to leave the buffers allocated (in which
* case they will be reused)
*
* For the definition of the projector theta and the Weyl compression see
* the notes in the module Pbnd.c
*
* All these programs act globally. In particular, the copy programs assume
* that the field number k is the same on all processes. This will normally
* be guaranteed by the calling program so that no check is made here
*
*******************************************************************************/

#define SCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "dirac.h"
#include "global.h"

#if (NPROC>1)

static int np,nmu[8],nbf[8],ofs[8];
static int ns,sfc[8];
static int itags=0,tags[8];
static weyl *wb=NULL,*snd_buf[8],*rcv_buf[8];
static MPI_Request snd_req[8],rcv_req[8];


static void get_tags(void)
{
   int i;
   
   if (itags==0)
   {
      for (i=0;i<8;i++)
         tags[i]=mpi_permanent_tag();

      itags=1;
   }
}


static void alloc_sbufs(void)
{
   int ifc,tag,saddr,raddr;
   weyl *wbb;

   wb=amalloc(BNDRY*sizeof(weyl),ALIGN);
   error(wb==NULL,1,"alloc_sbufs [scom.c]",
         "Unable to allocate communication buffers");

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;   

   nbf[0]=FACE0/2;
   nbf[1]=FACE0/2;
   nbf[2]=FACE1/2;
   nbf[3]=FACE1/2;
   nbf[4]=FACE2/2;
   nbf[5]=FACE2/2;   
   nbf[6]=FACE3/2;
   nbf[7]=FACE3/2;

   get_tags();
   ofs[0]=0;
   ns=0;
   wbb=wb;

   for (ifc=0;ifc<8;ifc++)
   {
      nmu[ifc]=cpr[ifc/2]&0x1;

      if (ifc>0)
         ofs[ifc]=ofs[ifc-1]+nbf[ifc-1];

      if (nbf[ifc]>0)
      {
         sfc[ns]=ifc;
         ns+=1;

         snd_buf[ifc]=wbb;
         wbb+=nbf[ifc];
         rcv_buf[ifc]=wbb;
         wbb+=nbf[ifc];

         tag=tags[ifc];
         saddr=npr[ifc];
         raddr=npr[ifc^0x1];

         MPI_Send_init((float*)(snd_buf[ifc]),12*nbf[ifc],MPI_FLOAT,saddr,
                       tag,MPI_COMM_WORLD,&snd_req[ifc]);
         MPI_Recv_init((float*)(rcv_buf[ifc]),12*nbf[ifc],MPI_FLOAT,raddr,
                       tag,MPI_COMM_WORLD,&rcv_req[ifc]);
      }
   }
}

#if (defined SSE)
#include "sse.h"

static void zip_weyl(int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2");

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c2));
   }
}


static void unzip_weyl(int vol,weyl *pk,spinor *pl)
{
   spinor *plm;

   __asm__ __volatile__ ("xorps %%xmm5, %%xmm5 \n\t"
                         "xorps %%xmm6, %%xmm6 \n\t"
                         "xorps %%xmm7, %%xmm7"
                         :
                         :
                         :
                         "xmm5", "xmm6", "xmm7");

   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2");

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;      
      
      __asm__ __volatile__ ("addps %%xmm0, %%xmm0 \n\t"
                            "addps %%xmm1, %%xmm1 \n\t"
                            "addps %%xmm2, %%xmm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm5, %3 \n\t"
                            "movaps %%xmm6, %4 \n\t"
                            "movaps %%xmm7, %5"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c3.c1),
                            "=m" ((*pl).c3.c3),
                            "=m" ((*pl).c4.c2));
   }
}

#elif (defined DH)
#include "dh.h"

static void zip_weyl(int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   
#pragma disjoint(*pk,*pl)
   __alignx(16,pk);
   __alignx(16,pl);
   
   _dh_vector_declare1(1);
   _dh_vector_declare1(2);
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      _dhs_vector_load((*pk).c1,1);
      _dhs_vector_load((*pk).c2,2);
      _dhs_weyl_vector_store2(*pl,1,2);
      
      pk+=1;
   }
}


static void unzip_weyl(int vol,weyl *pk,spinor *pl)
{
   spinor *plm;

#pragma disjoint(*pk,*pl)
   __alignx(16,pk);
   __alignx(16,pl);
   
   _dh_vector_declareu(1);
   _dh_vector_declareu(2);
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      _dhs_weyl_vector_load2(*pk,1,2);
      _dh_vector_add(1,1,1);
      _dh_vector_add(2,2,2);
      _dhs_vector_storeu((*pl).c1,1);
      _dhs_vector_storeu((*pl).c2,2);
      
      __stfps(&(*pl).c3.c1.re,  __cmplx(0.0, 0.0));
      __stfps(&(*pl).c3.c2.re,  __cmplx(0.0, 0.0));
      __stfps(&(*pl).c3.c3.re,  __cmplx(0.0, 0.0));
      __stfps(&(*pl).c4.c1.re,  __cmplx(0.0, 0.0));
      __stfps(&(*pl).c4.c2.re,  __cmplx(0.0, 0.0));
      __stfps(&(*pl).c4.c3.re,  __cmplx(0.0, 0.0));
      
      pk+=1;
   }
}

#else

static const weyl w0={{{0.0f}}};


static void zip_weyl(int vol,spinor *pk,weyl *pl)
{
   weyl *rpk,*rpm;
   
   rpk=(weyl*)(pk);
   rpm=pl+vol;
   
   for (;pl<rpm;pl++)
   {
      *pl=*rpk;
      rpk+=2;
   }
}


static void unzip_weyl(int vol,weyl *pk,spinor *pl)
{
   weyl *rpl,*rpm;
   
   rpl=(weyl*)(pl);
   rpm=pk+vol;
   
   for (;pk<rpm;pk++)
   {
      _vector_add((*rpl).c1,(*pk).c1,(*pk).c1);
      _vector_add((*rpl).c2,(*pk).c2,(*pk).c2);      
      *(rpl+1)=w0;
      rpl+=2;
   }
}

#endif

static void send_bufs(int ifc,int eo)
{
   int io;

   if (np==eo)
   {
      io=(ifc^nmu[ifc]);
      
      MPI_Start(&snd_req[io]);
   }
   else
   {
      io=(ifc^nmu[ifc])^0x1;
      
      MPI_Start(&rcv_req[io]);
   }
}


static void wait_bufs(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   if (np==eo)
   {
      io=(ifc^nmu[ifc]);
      
      MPI_Wait(&snd_req[io],&stat_snd);
   }
   else
   {
      io=(ifc^nmu[ifc])^0x1;
      
      MPI_Wait(&rcv_req[io],&stat_rcv);      
   }   
}


void cps_int_bnd(int k)
{
   int ifc,io;
   int n,m,eo;
   spinor *pk,*pkb;

   if (wb==NULL)
      alloc_sbufs();

   pk=ps[k][0];
   pkb=ps[k][VOLUME];
   m=0;
   eo=0;

   for (n=0;n<ns;n++)
   {
      if (n>0)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=ifc^nmu[ifc];
      assign_s2w[io^0x1](map+ofs[io^0x1],nbf[io],pk,snd_buf[io]);

      if (n>0)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   for (n=0;n<2;n++)
   {
      send_bufs(sfc[m],eo);
      wait_bufs(sfc[m],eo);
      m+=eo;
      eo^=0x1;
   }
   
   for (n=0;n<ns;n++)
   {
      if (m<ns)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=(ifc^nmu[ifc])^0x1;
      unzip_weyl(nbf[io],rcv_buf[io],pkb+ofs[io^0x1]);

      if (m<ns)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }
}


void cps_ext_bnd(int k)
{
   int ifc,io;
   int n,m,eo;
   spinor *pk,*pkb;

   if (wb==NULL)
      alloc_sbufs();

   pk=ps[k][0];
   pkb=ps[k][VOLUME];
   m=0;
   eo=0;
   
   for (n=0;n<ns;n++)
   {
      if (n>0)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=ifc^nmu[ifc];
      zip_weyl(nbf[io],pkb+ofs[io],snd_buf[io]);

      if (n>0)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   for (n=0;n<2;n++)
   {
      send_bufs(sfc[m],eo);
      wait_bufs(sfc[m],eo);
      m+=eo;
      eo^=0x1;
   }
   
   for (n=0;n<ns;n++)
   {
      if (m<ns)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=(ifc^nmu[ifc])^0x1;
      add_assign_w2s[io](map+ofs[io],nbf[io],rcv_buf[io],pk);      

      if (m<ns)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }
}


void free_sbufs(void)
{
   int n,ifc;

   if (wb==NULL)
      return;

   for (n=0;n<ns;n++)
   {
      ifc=sfc[n];
      MPI_Request_free(&snd_req[ifc]);
      MPI_Request_free(&rcv_req[ifc]);
   }

   afree(wb);
   wb=NULL;
}

#endif
