
/*******************************************************************************
*
* File sdcom.c
*
* Copyright (C) 2005, 2008 Martin Luescher, 2007 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication functions related to Qhat_dble
*
* Spinor fields are defined on the local lattices and the *even* exterior
* boundary points of these. Copying spinors to or from the boundary thus
* refers to these points only
*
*   void cpsd_int_bnd(int k)
*     Copies the spinors *psd[k][] at the interior boundary points of the 
*     local lattice to the corresponding points on the neighbouring processes.
*     This program compresses the Dirac spinors psi to Weyl spinors before
*     copying and expands the latter after copying to spinors chi such that
*     theta*chi=theta*psi where theta denotes the projector to the exterior
*     boundary of the target lattice
*
*   void cpsd_ext_bnd(int k)
*     Copies the spinors *psd[k][] at the exterior boundary points of the
*     local lattice to the neighbouring processes and *adds* them to the 
*     field on the matching points of the target lattices. This program
*     compresses the Dirac spinors psi to Weyl spinors before copying,
*     assuming that theta*psi=psi where theta denotes the projector to 
*     the exterior boundary of the local lattice
*
*   void free_sdbufs(void)
*     Frees the communication buffers used by the programs in this module
*
* Notes:
*
* The copy programs allocate the required communication buffers when
* needed. They can be freed using free_sdbufs(), but if the programs are
* called frequently, it is better to leave the buffers allocated (in which
* case they will be reused)
*
* For the definition of the projector theta and the Weyl compression see
* the notes in the module Pbnd_dble.c
*
* All these programs act globally. In particular, the copy programs assume
* that the field number k is the same on all processes. This will normally
* be guaranteed by the calling program so that no check is made here
*
*******************************************************************************/

#define SDCOM_C

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
static weyl_dble *wb=NULL,*snd_buf[8],*rcv_buf[8];
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


static void alloc_sdbufs(void)
{
   int ifc,tag,saddr,raddr;
   weyl_dble *wbb;

   wb=amalloc(BNDRY*sizeof(weyl_dble),ALIGN);
   error(wb==NULL,1,"alloc_sdbufs [sdcom.c]",
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

         MPI_Send_init((double*)(snd_buf[ifc]),12*nbf[ifc],MPI_DOUBLE,saddr,
                       tag,MPI_COMM_WORLD,&snd_req[ifc]);
         MPI_Recv_init((double*)(rcv_buf[ifc]),12*nbf[ifc],MPI_DOUBLE,raddr,
                       tag,MPI_COMM_WORLD,&rcv_req[ifc]);
      }
   }
}

#if (defined SSE2)
#include "sse2.h"

static void zip_weyl(int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"                            
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c2),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c1),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"                            
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c2),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c1),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c2.c3));
   }
}


static void unzip_weyl(int vol,weyl_dble *pk,spinor_dble *pl)
{
   spinor_dble *plm;

   __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                         "xorps %%xmm7, %%xmm7"
                         :
                         :
                         :
                         "xmm6", "xmm7");

   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"                            
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c2),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c1),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
      
      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;      

      __asm__ __volatile__ ("addpd %%xmm0, %%xmm0 \n\t"
                            "addpd %%xmm1, %%xmm1 \n\t"
                            "addpd %%xmm2, %%xmm2 \n\t"
                            "addpd %%xmm3, %%xmm3 \n\t"
                            "addpd %%xmm4, %%xmm4 \n\t"
                            "addpd %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"                            
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c2),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c1),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c2.c3));      
      
      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm6, %2 \n\t"
                            "movapd %%xmm7, %3 \n\t"
                            "movapd %%xmm6, %4 \n\t"
                            "movapd %%xmm7, %5"
                            :
                            "=m" ((*pl).c3.c1),
                            "=m" ((*pl).c3.c2),
                            "=m" ((*pl).c3.c3),
                            "=m" ((*pl).c4.c1),
                            "=m" ((*pl).c4.c2),
                            "=m" ((*pl).c4.c3));      
   }
}

#elif (defined DH)
#include "dh.h"

static void zip_weyl(int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   
#pragma disjoint(*pk,*pl)
   __alignx(16,pk);
   __alignx(16,pl);
   
   _dh_vector_declare1(1);
   _dh_vector_declare1(2);
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      _dh_vector_load((*pk).c1,1);
      _dh_vector_load((*pk).c2,2);
      _dh_weyl_vector_store2(*pl,1,2);
      
      pk+=1;
   }
}


static void unzip_weyl(int vol,weyl_dble *pk,spinor_dble *pl)
{
   spinor_dble *plm;

#pragma disjoint(*pk,*pl)
   __alignx(16,pk);
   __alignx(16,pl);
   
   _dh_vector_declareu(1);
   _dh_vector_declareu(2);
   
   plm=pl+vol;
   
   for (;pl<plm;pl++)
   {
      _dh_weyl_vector_load2(*pk,1,2);
      _dh_vector_add(1,1,1);
      _dh_vector_add(2,2,2);
      _dh_vector_storeu((*pl).c1,1);
      _dh_vector_storeu((*pl).c2,2);
      
      __stfpd(&(*pl).c3.c1.re,  __cmplx(0.0, 0.0));
      __stfpd(&(*pl).c3.c2.re,  __cmplx(0.0, 0.0));
      __stfpd(&(*pl).c3.c3.re,  __cmplx(0.0, 0.0));
      __stfpd(&(*pl).c4.c1.re,  __cmplx(0.0, 0.0));
      __stfpd(&(*pl).c4.c2.re,  __cmplx(0.0, 0.0));
      __stfpd(&(*pl).c4.c3.re,  __cmplx(0.0, 0.0));
      
      pk+=1;
   }
}

#else

static const weyl_dble w0={{{0.0}}};


static void zip_weyl(int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *rpk,*rpm;
   
   rpk=(weyl_dble*)(pk);
   rpm=pl+vol;
   
   for (;pl<rpm;pl++)
   {
      *pl=*rpk;
      rpk+=2;
   }
}


static void unzip_weyl(int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *rpl,*rpm;
   
   rpl=(weyl_dble*)(pl);
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


void cpsd_int_bnd(int k)
{
   int ifc,io;
   int n,m,eo;
   spinor_dble *pk,*pkb;

   if (wb==NULL)
      alloc_sdbufs();

   pk=psd[k][0];
   pkb=psd[k][VOLUME];
   m=0;
   eo=0;

   for (n=0;n<ns;n++)
   {
      if (n>0)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=ifc^nmu[ifc];
      assign_sd2wd[io^0x1](map+ofs[io^0x1],nbf[io],pk,snd_buf[io]);

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


void cpsd_ext_bnd(int k)
{
   int ifc,io;
   int n,m,eo;
   spinor_dble *pk,*pkb;

   if (wb==NULL)
      alloc_sdbufs();

   pk=psd[k][0];
   pkb=psd[k][VOLUME];
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
      add_assign_wd2sd[io](map+ofs[io],nbf[io],rcv_buf[io],pk);      

      if (m<ns)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }
}


void free_sdbufs(void)
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
