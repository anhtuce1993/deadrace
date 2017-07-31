
/*******************************************************************************
*
* File linalg_dble.c
*
* Copyright (C) 2005, 2007 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic linear algebra routines for double-precision Dirac fields
*
* The externally accessible functions are
*
*   complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *pk,
*                                 spinor_dble *pl)
*     Computes the scalar product of the fields pk[] and pl[]
*
*   double norm_square_dble(int vol,int icom,spinor_dble *pk)
*     Computes the square of the norm of the field pk[]
*
*   void mulc_spinor_add_dble(int vol,spinor_dble *pk,spinor_dble *pl,
*                             complex_dble z)
*     Replaces pk[] by pk[]+z*pl[]
*
*   double normalize_dble(int vol,int icom,spinor_dble *pk)
*     Replaces pk[] by pk[]/||pk|| and returns the norm ||pk||
*
*   void mulg5_dble(int vol,spinor_dble *pk)
*     Multiplies the field pk[] with gamma_5
*
*   void mulmg5_dble(int vol,spinor_dble *pk)
*     Multiplies the field pk[] with -gamma_5
*
* Notes:
*
* All these programs operate on arrays of spinor fields whose base address
* is passed through the arguments. The length of the arrays is specified by
* the parameter vol. Scalar products etc. are globally summed if the 
* parameter icom is equal to 1. In this case the calculated values are
* guaranteed to be exactly the same on all processes
*
* The base address of the spinor fields must be a multiple of 16 if SSE2
* instructions are to be used
*
*******************************************************************************/

#define LINALG_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "global.h"

#define MAX_LEVELS 8
#define BLK_LENGTH 16

static int cnt[MAX_LEVELS];
static double smx[MAX_LEVELS],smy[MAX_LEVELS];

#if (defined SSE2)

#include "sse2.h"

static sse_double tmp1,tmp2;
static const sse_double sgn={-1.0,-1.0};

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *pk,spinor_dble *pl)
{
   int n;
   double x,y;
   complex_dble w,z;
   weyl_dble *rpk,*rpl,*rpmax,*rptot;

   rpk=(weyl_dble*)(pk);
   rpl=(weyl_dble*)(pl);
   rpmax=rpk;
   rptot=rpk+2*vol;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
      smy[n]=0.0;
   }

   for (;rpmax<rptot;)
   {
      rpmax+=BLK_LENGTH;
      if (rpmax>rptot)
         rpmax=rptot;
      x=0.0;
      y=0.0;

      for (;rpk<rpmax;)
      {
         __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                               "movapd %1, %%xmm1 \n\t"
                               "movapd %2, %%xmm2 \n\t"
                               "movapd %3, %%xmm3 \n\t"
                               "movapd %4, %%xmm4 \n\t"
                               "movapd %5, %%xmm5"
                               :
                               :
                               "m" ((*rpk).c1.c1),
                               "m" ((*rpk).c1.c2),
                               "m" ((*rpk).c1.c3),
                               "m" ((*rpk).c2.c1),
                               "m" ((*rpk).c2.c2),
                               "m" ((*rpk).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm5"
                               :
                               :
                               "m" ((*rpl).c1.c1),
                               "m" ((*rpl).c1.c2),
                               "m" ((*rpl).c1.c3),
                               "m" ((*rpl).c2.c1),
                               "m" ((*rpl).c2.c2),
                               "m" ((*rpl).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

         rpk+=4;
         _prefetch_weyl_dble(rpk);
         rpk-=4;
         
         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm3 \n\t"
                               "addpd %%xmm3, %%xmm5"
                               :
                               :
                               :
                               "xmm1", "xmm3", "xmm5");
                               
         __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                               "movapd %1, %%xmm1 \n\t"
                               "movapd %2, %%xmm2 \n\t"
                               "movapd %3, %%xmm3 \n\t"
                               "movapd %4, %%xmm4 \n\t"
                               "movapd %5, %%xmm6"
                               :
                               :
                               "m" ((*rpk).c1.c1),
                               "m" ((*rpk).c1.c2),
                               "m" ((*rpk).c1.c3),
                               "m" ((*rpk).c2.c1),
                               "m" ((*rpk).c2.c2),
                               "m" ((*rpk).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");
         
         __asm__ __volatile__ ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"
                               "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                               "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                               "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                               "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                               "shufpd $0x1, %%xmm6, %%xmm6"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");

         __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                               "mulpd %1, %%xmm1 \n\t"
                               "mulpd %2, %%xmm2 \n\t"
                               "mulpd %3, %%xmm3 \n\t"
                               "mulpd %4, %%xmm4 \n\t"
                               "mulpd %5, %%xmm6"
                               :
                               :
                               "m" ((*rpl).c1.c1),
                               "m" ((*rpl).c1.c2),
                               "m" ((*rpl).c1.c3),
                               "m" ((*rpl).c2.c1),
                               "m" ((*rpl).c2.c2),
                               "m" ((*rpl).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");

         rpl+=4;
         _prefetch_weyl_dble(rpl);
         rpl-=3;
         
         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm6 \n\t"
                               "addpd %%xmm1, %%xmm3 \n\t"
                               "addpd %%xmm3, %%xmm6 \n\t"
                               "movapd %%xmm5, %0 \n\t"
                               "movapd %%xmm6, %1"
                               :
                               "=m" (tmp1),                               
                               "=m" (tmp2)
                               :
                               :
                               "xmm1", "xmm3", "xmm6");

         x+=(tmp1.c1+tmp1.c2);         
         y+=(tmp2.c2-tmp2.c1);
         rpk+=1;
      }

      cnt[0]+=1;
      smx[0]+=x;
      smy[0]+=y;

      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         smy[n]+=smy[n-1];

         cnt[n-1]=0;
         smx[n-1]=0.0;
         smy[n-1]=0.0;
      }
   }

   x=0.0;
   y=0.0;

   for (n=0;n<MAX_LEVELS;n++)
   {
      x+=smx[n];
      y+=smy[n];
   }

   if ((icom==1)&&(NPROC>1))
   {
      w.re=x;
      w.im=y;      
      MPI_Reduce(&w.re,&z.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&z.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);  
   }
   else
   {
      z.re=x;
      z.im=y;
   }
   
   return z; 
}


double norm_square_dble(int vol,int icom,spinor_dble *pk)
{
   int n;
   double x,y;
   weyl_dble *rpk,*rpmax,*rptot;

   rpk=(weyl_dble*)(pk);
   rpmax=rpk;
   rptot=rpk+2*vol;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }

   for (;rpmax<rptot;)
   {
      rpmax+=BLK_LENGTH;
      if (rpmax>rptot)
         rpmax=rptot;
      x=0.0;

      for (;rpk<rpmax;)
      {
         __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                               "movapd %1, %%xmm1 \n\t"
                               "movapd %2, %%xmm2 \n\t"
                               "movapd %3, %%xmm3 \n\t"
                               "movapd %4, %%xmm4 \n\t"
                               "movapd %5, %%xmm5"
                               :
                               :
                               "m" ((*rpk).c1.c1),
                               "m" ((*rpk).c1.c2),
                               "m" ((*rpk).c1.c3),
                               "m" ((*rpk).c2.c1),
                               "m" ((*rpk).c2.c2),
                               "m" ((*rpk).c2.c3)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

         rpk+=4;
         _prefetch_weyl_dble(rpk);
         rpk-=3;
         
         __asm__ __volatile__ ("mulpd %%xmm0, %%xmm0 \n\t"
                               "mulpd %%xmm1, %%xmm1 \n\t"
                               "mulpd %%xmm2, %%xmm2 \n\t"
                               "mulpd %%xmm3, %%xmm3 \n\t"
                               "mulpd %%xmm4, %%xmm4 \n\t"
                               "mulpd %%xmm5, %%xmm5"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
         
         __asm__ __volatile__ ("addpd %%xmm0, %%xmm1 \n\t"
                               "addpd %%xmm2, %%xmm3 \n\t"
                               "addpd %%xmm4, %%xmm5 \n\t"
                               "addpd %%xmm1, %%xmm3 \n\t"
                               "addpd %%xmm3, %%xmm5 \n\t"
                               "movapd %%xmm5, %0"
                               :
                               "=m" (tmp1)
                               :
                               :
                               "xmm1", "xmm3", "xmm5");


         x+=(tmp1.c1+tmp1.c2);
      }

      cnt[0]+=1;
      smx[0]+=x;

      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];

         cnt[n-1]=0;
         smx[n-1]=0.0;
      }
   }

   x=0.0;

   for (n=0;n<MAX_LEVELS;n++)
      x+=smx[n];

   if ((icom==1)&&(NPROC>1))
   {
      MPI_Reduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return y;
   }
   else
      return x;
}


void mulc_spinor_add_dble(int vol,spinor_dble *pk,spinor_dble *pl,
                          complex_dble z)
{
   weyl_dble *rpk,*rpl,*rptot;

   _sse_load_cmplx_dble(z);   
   rpk=(weyl_dble*)(pk);
   rpl=(weyl_dble*)(pl);
   rptot=rpk+2*vol;
   
   for (;rpk<rptot;)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %%xmm0, %%xmm3 \n\t"
                            "movapd %%xmm1, %%xmm4 \n\t"
                            "movapd %%xmm2, %%xmm5"
                            :
                            :
                            "m" ((*rpl).c1.c1),
                            "m" ((*rpl).c1.c2),
                            "m" ((*rpl).c1.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      rpl+=4;
      _prefetch_weyl_dble(rpl);
      rpl-=4;
      
      __asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t"
                            "mulpd %%xmm6, %%xmm1 \n\t"
                            "mulpd %%xmm6, %%xmm2 \n\t"
                            "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                            "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                            "shufpd $0x1, %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
      
      __asm__ __volatile__ ("addpd %0, %%xmm0 \n\t"
                            "addpd %1, %%xmm1 \n\t"
                            "addpd %2, %%xmm2 \n\t"
                            "mulpd %%xmm7, %%xmm3 \n\t"
                            "mulpd %%xmm7, %%xmm4 \n\t"
                            "mulpd %%xmm7, %%xmm5 \n\t"
                            "addpd %%xmm0, %%xmm3 \n\t"
                            "addpd %%xmm1, %%xmm4 \n\t"
                            "addpd %%xmm2, %%xmm5"                            
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c2),
                            "m" ((*rpk).c1.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2"
                            :
                            :
                            "m" ((*rpl).c2.c1),
                            "m" ((*rpl).c2.c2),
                            "m" ((*rpl).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      
      __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                            "movapd %%xmm4, %1 \n\t"
                            "movapd %%xmm5, %2 \n\t"
                            "movapd %%xmm0, %%xmm3 \n\t"
                            "movapd %%xmm1, %%xmm4 \n\t"
                            "movapd %%xmm2, %%xmm5"
                            :
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c2),
                            "=m" ((*rpk).c1.c3)
                            :
                            :
                            "xmm3", "xmm4", "xmm5");

      rpk+=4;
      _prefetch_weyl_dble(rpk);
      rpk-=4;      
      
      __asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t"
                            "mulpd %%xmm6, %%xmm1 \n\t"
                            "mulpd %%xmm6, %%xmm2 \n\t"
                            "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                            "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                            "shufpd $0x1, %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("addpd %0, %%xmm0 \n\t"
                            "addpd %1, %%xmm1 \n\t"
                            "addpd %2, %%xmm2 \n\t"
                            "mulpd %%xmm7, %%xmm3 \n\t"
                            "mulpd %%xmm7, %%xmm4 \n\t"
                            "mulpd %%xmm7, %%xmm5 \n\t"
                            "addpd %%xmm0, %%xmm3 \n\t"
                            "addpd %%xmm1, %%xmm4 \n\t"
                            "addpd %%xmm2, %%xmm5"                            
                            :
                            :
                            "m" ((*rpk).c2.c1),
                            "m" ((*rpk).c2.c2),
                            "m" ((*rpk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm3, %0 \n\t"
                            "movapd %%xmm4, %1 \n\t"
                            "movapd %%xmm5, %2"
                            :
                            "=m" ((*rpk).c2.c1),
                            "=m" ((*rpk).c2.c2),
                            "=m" ((*rpk).c2.c3));

      rpk+=1;
      rpl+=1;
   }
}


double normalize_dble(int vol,int icom,spinor_dble *pk)
{
   double r,ri;
   weyl_dble *rpk,*rptot;

   r=norm_square_dble(vol,icom,pk);
   r=sqrt(r);

   if (error_loc(r==0.0,1,"normalize_dble [linalg_dble.c]",
                 "Vector has vanishing norm")!=0)
      return 0.0;

   ri=1.0/r;

   __asm__ __volatile__ ("movsd %0, %%xmm6 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "movapd %%xmm6, %%xmm7"
                         :
                         :
                         "m" (ri)
                         :
                         "xmm6", "xmm7");
   
   rpk=(weyl_dble*)(pk);
   rptot=rpk+2*vol;   
   
   for (;rpk<rptot;)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c2),
                            "m" ((*rpk).c1.c3),
                            "m" ((*rpk).c2.c1),
                            "m" ((*rpk).c2.c2),
                            "m" ((*rpk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      rpk+=4;
      _prefetch_weyl_dble(rpk);
      rpk-=4;      
      
      __asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t"
                            "mulpd %%xmm7, %%xmm1 \n\t"
                            "mulpd %%xmm6, %%xmm2 \n\t"
                            "mulpd %%xmm7, %%xmm3 \n\t"
                            "mulpd %%xmm6, %%xmm4 \n\t"
                            "mulpd %%xmm7, %%xmm5"
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
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c2),
                            "=m" ((*rpk).c1.c3),                            
                            "=m" ((*rpk).c2.c1),
                            "=m" ((*rpk).c2.c2),
                            "=m" ((*rpk).c2.c3));

      rpk+=1;
   }

   return r;
}


void mulg5_dble(int vol,spinor_dble *pk)
{
   weyl_dble *rpk,*rpm;
   
   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %%xmm6, %%xmm7"
                         :
                         :
                         "m" (sgn)
                         :
                         "xmm6", "xmm7");

   rpk=(weyl_dble*)(pk);
   rpm=rpk+2*vol;
   rpk+=1;
   
   for (;rpk<rpm;)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c2),
                            "m" ((*rpk).c1.c3),
                            "m" ((*rpk).c2.c1),
                            "m" ((*rpk).c2.c2),
                            "m" ((*rpk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
                            
      rpk+=8;
      _prefetch_weyl_dble(rpk);
      rpk-=8;      
      
      __asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t"
                            "mulpd %%xmm7, %%xmm1 \n\t"
                            "mulpd %%xmm6, %%xmm2 \n\t"
                            "mulpd %%xmm7, %%xmm3 \n\t"
                            "mulpd %%xmm6, %%xmm4 \n\t"
                            "mulpd %%xmm7, %%xmm5"
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
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c2),
                            "=m" ((*rpk).c1.c3),                            
                            "=m" ((*rpk).c2.c1),
                            "=m" ((*rpk).c2.c2),
                            "=m" ((*rpk).c2.c3));

      rpk+=2;
   }
}


void mulmg5_dble(int vol,spinor_dble *pk)
{
   weyl_dble *rpk,*rpm;
   
   __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                         "movapd %%xmm6, %%xmm7"
                         :
                         :
                         "m" (sgn)
                         :
                         "xmm6", "xmm7");

   rpk=(weyl_dble*)(pk);
   rpm=rpk+2*vol;
   
   for (;rpk<rpm;)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c2),
                            "m" ((*rpk).c1.c3),
                            "m" ((*rpk).c2.c1),
                            "m" ((*rpk).c2.c2),
                            "m" ((*rpk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      rpk+=8;
      _prefetch_weyl_dble(rpk);
      rpk-=8;
      
      __asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t"
                            "mulpd %%xmm7, %%xmm1 \n\t"
                            "mulpd %%xmm6, %%xmm2 \n\t"
                            "mulpd %%xmm7, %%xmm3 \n\t"
                            "mulpd %%xmm6, %%xmm4 \n\t"
                            "mulpd %%xmm7, %%xmm5"
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
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c2),
                            "=m" ((*rpk).c1.c3),                            
                            "=m" ((*rpk).c2.c1),
                            "=m" ((*rpk).c2.c2),
                            "=m" ((*rpk).c2.c3));

      rpk+=2;
   }
}

#else

complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *pk,spinor_dble *pl)
{
   int n;
   double x,y;
   complex_dble w,z;
   spinor_dble *rpk,*rpl,*rpmax,*rptot;

   rpk=pk;
   rpl=pl;
   rpmax=rpk;
   rptot=rpk+vol;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
      smy[n]=0.0;
   }
   
   for (;rpmax<rptot;)
   {
      rpmax+=BLK_LENGTH;
      if (rpmax>rptot)
         rpmax=rptot;
      x=0.0;
      y=0.0;

      for (;rpk<rpmax;)
      {
         x+=(_vector_prod_re((*rpk).c1,(*rpl).c1)+
               _vector_prod_re((*rpk).c2,(*rpl).c2)+
               _vector_prod_re((*rpk).c3,(*rpl).c3)+
               _vector_prod_re((*rpk).c4,(*rpl).c4));

         y+=(_vector_prod_im((*rpk).c1,(*rpl).c1)+
               _vector_prod_im((*rpk).c2,(*rpl).c2)+
               _vector_prod_im((*rpk).c3,(*rpl).c3)+
               _vector_prod_im((*rpk).c4,(*rpl).c4));

         rpk+=1;
         rpl+=1;
      }

      cnt[0]+=1;
      smx[0]+=x;
      smy[0]+=y;

      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         smy[n]+=smy[n-1];

         cnt[n-1]=0;
         smx[n-1]=0.0;
         smy[n-1]=0.0;
      }
   }

   x=0.0;
   y=0.0;

   for (n=0;n<MAX_LEVELS;n++)
   {
      x+=smx[n];
      y+=smy[n];
   }

   if ((icom!=1)||(NPROC==1))
   {
      z.re=x;
      z.im=y;
   }
   else
   {
      w.re=x;
      w.im=y;      
      MPI_Reduce(&w.re,&z.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&z.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);  
   }

   return z; 
}

 
double norm_square_dble(int vol,int icom,spinor_dble *pk)
{
   int n;
   double x,y;
   spinor_dble *rpk,*rpmax,*rptot;

   rpk=pk;
   rpmax=rpk;
   rptot=rpk+vol;

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }
   
   for (;rpmax<rptot;)
   {
      rpmax+=BLK_LENGTH;
      if (rpmax>rptot)
         rpmax=rptot;
      x=0.0;

      for (;rpk<rpmax;)
      {
         x+=(_vector_prod_re((*rpk).c1,(*rpk).c1)+
               _vector_prod_re((*rpk).c2,(*rpk).c2)+
               _vector_prod_re((*rpk).c3,(*rpk).c3)+
               _vector_prod_re((*rpk).c4,(*rpk).c4));

         rpk+=1;
      }

      cnt[0]+=1;
      smx[0]+=x;

      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];

         cnt[n-1]=0;
         smx[n-1]=0.0;
      }
   }

   x=0.0;

   for (n=0;n<MAX_LEVELS;n++)
      x+=smx[n];

   if ((icom!=1)||(NPROC==1))
      return x;
   else
   {
      MPI_Reduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return y;    
   }
}


void mulc_spinor_add_dble(int vol,spinor_dble *pk,spinor_dble *pl,
                          complex_dble z)
{
   spinor_dble *rpk,*rpl;

   rpk=pk;
   rpl=pl;
   
   for (;rpk<(pk+vol);)
   {
      _vector_mulc_assign((*rpk).c1,z,(*rpl).c1);
      _vector_mulc_assign((*rpk).c2,z,(*rpl).c2);
      _vector_mulc_assign((*rpk).c3,z,(*rpl).c3);
      _vector_mulc_assign((*rpk).c4,z,(*rpl).c4);

      rpk+=1;
      rpl+=1;
   }
}


double normalize_dble(int vol,int icom,spinor_dble *pk)
{
   double r,ri;
   spinor_dble *rpk;

   r=norm_square_dble(vol,icom,pk);
   r=sqrt(r);
   
   if (error_loc(r==0.0,1,"normalize_dble [linalg_dble.c]",
       "Vector has vanishing norm")!=0)
      return 0.0;

   ri=1.0/r;
   rpk=pk;
   
   for (;rpk<(pk+vol);)
   {
      _vector_mul((*rpk).c1,ri,(*rpk).c1);
      _vector_mul((*rpk).c2,ri,(*rpk).c2);
      _vector_mul((*rpk).c3,ri,(*rpk).c3);
      _vector_mul((*rpk).c4,ri,(*rpk).c4);

      rpk+=1;
   }

   return r;
}


void mulg5_dble(int vol,spinor_dble *pk)
{
   weyl_dble *rpk,*rpm;

   rpk=(weyl_dble*)(pk);
   rpm=rpk+2*vol;
   rpk+=1;
   
   for (;rpk<rpm;rpk+=2)
   {
      (*rpk).c1.c1.re=-(*rpk).c1.c1.re;
      (*rpk).c1.c1.im=-(*rpk).c1.c1.im;      
      (*rpk).c1.c2.re=-(*rpk).c1.c2.re;
      (*rpk).c1.c2.im=-(*rpk).c1.c2.im;      
      (*rpk).c1.c3.re=-(*rpk).c1.c3.re;
      (*rpk).c1.c3.im=-(*rpk).c1.c3.im;      
      (*rpk).c2.c1.re=-(*rpk).c2.c1.re;
      (*rpk).c2.c1.im=-(*rpk).c2.c1.im;      
      (*rpk).c2.c2.re=-(*rpk).c2.c2.re;
      (*rpk).c2.c2.im=-(*rpk).c2.c2.im;      
      (*rpk).c2.c3.re=-(*rpk).c2.c3.re;
      (*rpk).c2.c3.im=-(*rpk).c2.c3.im;      
   }
}


void mulmg5_dble(int vol,spinor_dble *pk)
{
   weyl_dble *rpk,*rpm;

   rpk=(weyl_dble*)(pk);
   rpm=rpk+2*vol;
   
   for (;rpk<rpm;rpk+=2)
   {
      (*rpk).c1.c1.re=-(*rpk).c1.c1.re;
      (*rpk).c1.c1.im=-(*rpk).c1.c1.im;      
      (*rpk).c1.c2.re=-(*rpk).c1.c2.re;
      (*rpk).c1.c2.im=-(*rpk).c1.c2.im;      
      (*rpk).c1.c3.re=-(*rpk).c1.c3.re;
      (*rpk).c1.c3.im=-(*rpk).c1.c3.im;      
      (*rpk).c2.c1.re=-(*rpk).c2.c1.re;
      (*rpk).c2.c1.im=-(*rpk).c2.c1.im;      
      (*rpk).c2.c2.re=-(*rpk).c2.c2.re;
      (*rpk).c2.c2.im=-(*rpk).c2.c2.im;      
      (*rpk).c2.c3.re=-(*rpk).c2.c3.re;
      (*rpk).c2.c3.im=-(*rpk).c2.c3.im;      
   }
}

#endif

