
/*******************************************************************************
*
* File sw_term.c
*
* Copyright (C) 2006 Luigi Del Debbio, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the SW-term on the full lattice
*
* The externally accessible functions are
*
*   void sw_term(void)
*     Computes the SW-term for the current double-precision gauge field
*     and assigns the matrix to the global swd array
*
*   void free_swbufs(void)
*     Frees the buffers used for the computation of the SW-term. The program
*     sw_term allocates these buffers automatically if they are not already
*     allocated (in which case they are reused)
*
* Notes:
*
* The calculated matrix at the point x is
*
*    gamma_5*[4+m0+csw*(i/4)*sigma_{mu nu}*Fhat_{mu nu}(x)]
*
* where
*
*    sigma_{mu nu}=(i/2)*[gamma_mu,gamma_nu]
*
* and Fhat_{mu nu} is the standard lattice expression for the gauge field
* tensor. The upper and lower 6x6 blocks of the matrix are stored in the
* pauli_dble structures swd[2*ix] and swd[2*ix+1]
*
* The parameters m0 and csw are obtained from the parameter data base
* using the programs provided by the module start/parms.c
*
*******************************************************************************/

#define SW_TERM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "misc.h"
#include "start.h"
#include "sw_term.h"
#include "global.h"

typedef struct
{
   int mu,nu;
   int nbf[2];
   int *ipl;
   int *ibf[2];
   int saddr[2],raddr[2];
   u3_alg_dble *snd_buf[2],*rcv_buf;
} comdat_t;

#if (defined SSE2)
#include "sse2.h"
static su3_dble w1,w2 __attribute__ ((aligned (16)));
static const u3_alg_dble qd0 __attribute__ ((aligned (8))) ={0.0};
#else
static su3_dble w1,w2;
static const u3_alg_dble qd0={0.0};
#endif

static u3_alg_dble *Qarr=NULL;
static comdat_t comdat[6];


static void set_nbf(void)
{
   int bs[4],nfc[4],n,mu,nu;
   int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
   comdat_t *c;

   bs[0]=L0;
   bs[1]=L1;
   bs[2]=L2;
   bs[3]=L3;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   c=comdat;

   for (n=0;n<6;n++)
   {
      mu=plns[n][0];
      nu=plns[n][1];

      (*c).mu=mu;
      (*c).nu=nu;

      if (nfc[nu]>0)
         (*c).nbf[0]=nfc[mu]+(nfc[mu]/bs[nu]);
      else
         (*c).nbf[0]=nfc[mu];
      (*c).nbf[1]=nfc[nu];

      (*c).saddr[0]=npr[2*mu+1];
      (*c).raddr[0]=npr[2*mu];
      (*c).saddr[1]=npr[2*nu+1];
      (*c).raddr[1]=npr[2*nu];

      c+=1;
   }
}


static void alloc_idx(void)
{
   int *ipl,nidx;
   comdat_t *c,*cm;

   nidx=0;
   cm=comdat+6;

   for (c=comdat;c<cm;c++)
      nidx+=(3*VOLUME+(*c).nbf[0]+(*c).nbf[1]);

   ipl=amalloc(nidx*sizeof(int),3);
   error(ipl==NULL,1,"alloc_idx [sw_term.c]",
         "Unable to allocate index arrays");

   for (c=comdat;c<cm;c++)
   {
      (*c).ipl=ipl;
      ipl+=(3*VOLUME);
      (*c).ibf[0]=ipl;
      ipl+=(*c).nbf[0];
      (*c).ibf[1]=ipl;
      ipl+=(*c).nbf[1];
   }
}


static void set_idx(void)
{
   int ix,iy,iw,iz,mu,nu,n;
   int nfc[4],ofs[4],icn,nbf0;
   int *ipl,*ibf[2];
   comdat_t *c,*cm;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   ofs[0]=FACE0/2;
   ofs[1]=ofs[0]+(FACE0+FACE1)/2;
   ofs[2]=ofs[1]+(FACE1+FACE2)/2;
   ofs[3]=ofs[2]+(FACE2+FACE3)/2;

   cm=comdat+6;

   for (c=comdat;c<cm;c++)
   {
      mu=(*c).mu;
      nu=(*c).nu;
      ipl=(*c).ipl;
      ibf[0]=(*c).ibf[0];
      ibf[1]=(*c).ibf[1];
      nbf0=(*c).nbf[0];
      icn=0;

      for (ix=0;ix<VOLUME;ix++)
      {
         iy=iup[ix][mu];
         iw=iup[ix][nu];

         if (iy>=(VOLUME+(BNDRY/2)))
         {
            n=iy-VOLUME-(BNDRY/2)-ofs[mu]+(nfc[mu]/2);
            *ipl=VOLUME+n;
            ibf[0][n]=map[iy-VOLUME];
         }
         else if (iy>=VOLUME)
         {
            n=iy-VOLUME-ofs[mu];
            *ipl=VOLUME+n;
            ibf[0][n]=map[iy-VOLUME];
         }
         else
            *ipl=iy;

         ipl+=1;

         if (iw>=(VOLUME+(BNDRY/2)))
         {
            n=iw-VOLUME-(BNDRY/2)-ofs[nu]+(nfc[nu]/2);
            *ipl=VOLUME+nbf0+n;
            ibf[1][n]=map[iw-VOLUME];
         }
         else if (iw>=VOLUME)
         {
            n=iw-VOLUME-ofs[nu];
            *ipl=VOLUME+nbf0+n;
            ibf[1][n]=map[iw-VOLUME];
         }
         else
            *ipl=iw;

         ipl+=1;

         if ((iy>=VOLUME)&&(iw>=VOLUME))
         {
            n=nfc[mu]+icn;
            *ipl=VOLUME+n;

            iz=map[iy-VOLUME];
            iz=iup[iz][nu];

            if (iz>=(VOLUME+(BNDRY/2)))
               ibf[0][n]=iz-(BNDRY/2)-ofs[nu]+(nfc[nu]/2)+nbf0;
            else
               ibf[0][n]=iz-ofs[nu]+nbf0;

            icn+=1;
         }
         else if (iy>=VOLUME)
         {
            iz=iup[iw][mu];

            if (iz>=(VOLUME+(BNDRY/2)))
               *ipl=iz-(BNDRY/2)-ofs[mu]+(nfc[mu]/2);
            else
               *ipl=iz-ofs[mu];
         }
         else if (iw>=VOLUME)
         {
            iz=iup[iy][nu];

            if (iz>=(VOLUME+(BNDRY/2)))
               *ipl=iz-(BNDRY/2)-ofs[nu]+(nfc[nu]/2)+nbf0;
            else
               *ipl=iz-ofs[nu]+nbf0;
         }
         else
            *ipl=iup[iy][nu];

         ipl+=1;
      }
   }
}


static void alloc_Qarr(void)
{
   int n,nc;
   u3_alg_dble *qd,*qm;
   comdat_t *c,*cm;

   n=VOLUME;
   cm=comdat+6;

   for (c=comdat;c<cm;c++)
   {
      nc=(*c).nbf[0]+(*c).nbf[1];

      if ((*c).nbf[0]>(*c).nbf[1])
         nc+=(*c).nbf[0];
      else
         nc+=(*c).nbf[1];

      n+=nc;
   }

   Qarr=amalloc(n*sizeof(u3_alg_dble),ALIGN);
   error(Qarr==NULL,1,"alloc_Qarr [sw_term.c]",
         "Unable to allocate work space");

   qm=Qarr+n;

   for (qd=Qarr;qd<qm;qd++)
      *qd=qd0;

   for (c=comdat;c<cm;c++)
   {
      (*c).snd_buf[0]=Qarr+VOLUME;
      (*c).snd_buf[1]=(*c).snd_buf[0]+(*c).nbf[0];
      (*c).rcv_buf   =(*c).snd_buf[1]+(*c).nbf[1];
   }
}


static void alloc_swbufs(void)
{
   error_root(sizeof(u3_alg_dble)!=(9*sizeof(double)),1,
              "alloc_swbufs [sw_term.c]",
              "The u3alg_dble structures are not properly packed");
   
   set_nbf();
   alloc_idx();
   set_idx();
   alloc_Qarr();
}


static void send_bufs(comdat_t *c,int dir)
{
   int tag,np,nbf,saddr,raddr;
   double *sbuf,*rbuf;
   MPI_Status stat;

   tag=mpi_tag();
   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;

   nbf=9*(*c).nbf[dir];
   saddr=(*c).saddr[dir];
   raddr=(*c).raddr[dir];
   sbuf=(double*)((*c).snd_buf[dir]);
   rbuf=(double*)((*c).rcv_buf);

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


static void unpack_bufs(comdat_t *c,int dir)
{
   int *ibf;
   u3_alg_dble *r,*rr,*rm,*q,*qq;

   ibf=(*c).ibf[dir];
   rr=(*c).rcv_buf;
   rm=rr+(*c).nbf[dir];
   qq=Qarr+(*ibf);

   for (;rr<rm;)
   {
      r=rr;
      q=qq;

      rr+=1;
      ibf+=(rr<rm);
      qq=Qarr+(*ibf);

#if (defined SSE2)
      _prefetch_su3(rr);
      _prefetch_su3(qq);
#endif

      (*q).c1+=(*r).c1;
      (*q).c2+=(*r).c2;
      (*q).c3+=(*r).c3;
      (*q).c4+=(*r).c4;
      (*q).c5+=(*r).c5;
      (*q).c6+=(*r).c6;
      (*q).c7+=(*r).c7;
      (*q).c8+=(*r).c8;
      (*q).c9+=(*r).c9;
   }
}


static void build_Q(comdat_t *c)
{
   int ix,mu,nu,dir,*ipl;
   su3_dble *ud[4];
   u3_alg_dble *qd,*qm,*Q0,*Q1;

   qd=Qarr+VOLUME;
   qm=qd+(*c).nbf[0]+(*c).nbf[1];

   for (;qd<qm;qd++)
      *qd=qd0;

   mu=(*c).mu;
   nu=(*c).nu;
   ipl=(*c).ipl;
   plaq_ud(0,mu,nu,ud);
   
   for (ix=0;ix<VOLUME;ix++)
   {
      Q0=Qarr+ix;
      Q1=Qarr+ipl[2];
      
#if (defined SSE2)
      _prefetch_su3(Q0);
      _prefetch_su3(Q1);      
#endif      

      su3xsu3(ud[0],ud[1],&w1);
      su3dagxsu3dag(ud[3],ud[2],&w2);
      add_prod2u3alg(&w1,&w2,Q0);
      add_prod2u3alg(&w2,&w1,Q1);

      Q0=Qarr+ipl[0];
      Q1=Qarr+ipl[1];
      ipl+=3;

#if (defined SSE2)
      _prefetch_su3(Q0);
      _prefetch_su3(Q1);      
#endif 
      
      su3dagxsu3(ud[2],ud[0],&w1);
      su3xsu3dag(ud[1],ud[3],&w2);

      plaq_ud(ix+(ix<(VOLUME-1)),mu,nu,ud);

#if (defined SSE2)
      _prefetch_su3_dble(ud[0]);
      _prefetch_su3_dble(ud[1]);      
#endif      

      add_prod2u3alg(&w2,&w1,Q0);              

#if (defined SSE2)
      _prefetch_su3_dble(ud[2]);
      _prefetch_su3_dble(ud[3]);      
#endif      
      
      add_prod2u3alg(&w1,&w2,Q1);
   }

   for (dir=0;dir<2;dir++)
   {
      if ((*c).nbf[dir]!=0)
      {
         send_bufs(c,dir);
         unpack_bufs(c,dir);
      }
   }
}


static void u3_alg2pauli1(int vol,u3_alg_dble *X,pauli_dble *m)
{
   u3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;)
   {
#if (defined SSE2)
      X+=2;
      _prefetch_su3(X);
      X-=2;
      m+=8;
      _prefetch_pauli_dble(m);
      m-=8;
#endif
      (*m).u[10]=-(*X).c1;

      (*m).u[12]=-(*X).c5;
      (*m).u[13]= (*X).c4;
      (*m).u[14]=-(*X).c7;
      (*m).u[15]= (*X).c6;

      (*m).u[18]=-(*X).c5;
      (*m).u[19]=-(*X).c4;
      (*m).u[20]=-(*X).c2;

      (*m).u[22]=-(*X).c9;
      (*m).u[23]= (*X).c8;
      (*m).u[24]=-(*X).c7;
      (*m).u[25]=-(*X).c6;
      (*m).u[26]=-(*X).c9;
      (*m).u[27]=-(*X).c8;
      (*m).u[28]=-(*X).c3;

      *X=qd0;

      m+=2;
      X+=1;
   }
}


static void u3_alg2pauli2(int vol,u3_alg_dble *X,pauli_dble *m)
{
   u3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;)
   {
#if (defined SSE2)
      X+=2;
      _prefetch_su3(X);
      X-=2;
      m+=8;
      _prefetch_pauli_dble(m);
      m-=8;
#endif   
      (*m).u[11] =(*X).c1;
      (*m).u[12]+=(*X).c4;
      (*m).u[13]+=(*X).c5;
      (*m).u[14]+=(*X).c6;
      (*m).u[15]+=(*X).c7;

      (*m).u[18]-=(*X).c4;
      (*m).u[19]+=(*X).c5;

      (*m).u[21] =(*X).c2;
      (*m).u[22]+=(*X).c8;
      (*m).u[23]+=(*X).c9;
      (*m).u[24]-=(*X).c6;
      (*m).u[25]+=(*X).c7;
      (*m).u[26]-=(*X).c8;
      (*m).u[27]+=(*X).c9;

      (*m).u[29] =(*X).c3;

      *X=qd0;

      m+=2;
      X+=1;      
   }
}


static void u3_alg2pauli3(int vol,u3_alg_dble *X,pauli_dble *m)
{
   u3_alg_dble *Xm;

   Xm=X+vol;

   for (;X<Xm;)
   {
#if (defined SSE2)
      X+=2;
      _prefetch_su3(X);
      X-=2;
      m+=8;
      _prefetch_pauli_dble(m);
      m-=8;
#endif   
      (*m).u[ 0]=-(*X).c1;
      (*m).u[ 1]=-(*X).c2;
      (*m).u[ 2]=-(*X).c3;
      (*m).u[ 3]= (*X).c1;
      (*m).u[ 4]= (*X).c2;
      (*m).u[ 5]= (*X).c3;
      (*m).u[ 6]=-(*X).c5;
      (*m).u[ 7]= (*X).c4;
      (*m).u[ 8]=-(*X).c7;
      (*m).u[ 9]= (*X).c6;

      (*m).u[16]=-(*X).c9;
      (*m).u[17]= (*X).c8;

      (*m).u[30]= (*X).c5;
      (*m).u[31]=-(*X).c4;
      (*m).u[32]= (*X).c7;
      (*m).u[33]=-(*X).c6;
      (*m).u[34]= (*X).c9;
      (*m).u[35]=-(*X).c8;

      *X=qd0;

      m+=2;
      X+=1;      
   }
}


static void build_F(void)
{
   build_Q(comdat);
   u3_alg2pauli1(VOLUME,Qarr,swd);

   build_Q(comdat+1);
   u3_alg2pauli2(VOLUME,Qarr,swd);

   build_Q(comdat+2);
   u3_alg2pauli3(VOLUME,Qarr,swd);

   build_Q(comdat+3);
   u3_alg2pauli1(VOLUME,Qarr,swd+1);

   build_Q(comdat+4);
   u3_alg2pauli2(VOLUME,Qarr,swd+1);

   build_Q(comdat+5);
   u3_alg2pauli3(VOLUME,Qarr,swd+1);
}


void sw_term(void)
{
   double c1,c2;
   double *rup,*rdn,*rm,r1,r2;
   pauli_dble *pd,*pm;
   lat_parms_t lat;
   
   if (swd==NULL)
      alloc_swd();

   if (Qarr==NULL)
      alloc_swbufs();

   copy_bnd_ud();
   build_F();

   lat=lat_parms();
   c1=4.0+lat.m0;
   c2=lat.csw/16.0;
   pm=swd+2*VOLUME;

   for (pd=swd;pd<pm;pd+=2)
   {
      rup=(*(pd  )).u;
      rdn=(*(pd+1)).u;
      rm=rup+6;

      for (;rup<rm;rup+=1)
      {
         r1=*rup;
         r2=*rdn;
         *rup= c1+c2*(r1-r2);
         *rdn=-c1+c2*(r1+r2);
         rdn+=1;
      }

      rm=rup+30;

      for (;rup<rm;rup+=1)
      {
         r1=*rup;
         r2=*rdn;
         *rup=c2*(r1-r2);
         *rdn=c2*(r1+r2);
         rdn+=1;
      }
   }
}


void free_swbufs(void)
{
   if (Qarr==NULL)
      return;

   afree(Qarr);
   afree((*comdat).ipl);
   Qarr=NULL;
}
