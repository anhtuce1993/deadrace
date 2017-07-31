
/*******************************************************************************
*
* File linalg.c
*
* Copyright (C) 2005, 2007 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic linear algebra routines for single-precision Dirac fields
*
* The externally accessible functions are
*
*   complex spinor_prod(int vol,int icom,spinor *pk,spinor *pl)
*     Computes the scalar product of the fields pk[] and pl[]
*
*   float norm_square(int vol,int icom,spinor *pk)
*     Computes the square of the norm of the field pk[]
*
*   void mulc_spinor_add(int vol,spinor *pk,spinor *pl,complex z)
*     Replaces pk[] by pk[]+z*pl[]
*
*   float normalize(int vol,int icom,spinor *pk)
*     Replaces pk[] by pk[]/||pk|| and returns the norm ||pk||
*
*   void mulg5(int vol,spinor *pk)
*     Multiplies the field pk[] with gamma_5
*
*   void mulmg5(int vol,spinor *pk)
*     Multiplies the field pk[] with -gamma_5
*
* Notes:
*
* All these programs operate on arrays of spinor fields whose base address
* is passed through the arguments. The length of the arrays is specified 
* by the parameter vol. Scalar products etc. are globally summed if the 
* parameter icom is equal to 1. In this case the calculated values are
* guaranteed to be exactly the same on all processes
*
* The base address of the spinor fields must be a multiple of 16 if SSE
* instructions are to be used
*
*******************************************************************************/

#define LINALG_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "global.h"

#if (defined SSE)
#include "sse.h"

static sse_float tmp1,tmp2;


complex spinor_prod(int vol,int icom,spinor *pk,spinor *pl)
{
   complex z;
   double x,y;
   complex_dble v,w;
   spinor *pkm;
  
   x=0.0;
   y=0.0;
   pkm=pk+vol;

   for (;pk<pkm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %3, %%xmm3 \n\t"
                            "movaps %4, %%xmm4 \n\t"
                            "movaps %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c3.c1),
                            "m" ((*pk).c3.c3),
                            "m" ((*pk).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
      
      __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                            "mulps %1, %%xmm1 \n\t"
                            "mulps %2, %%xmm2 \n\t"
                            "mulps %3, %%xmm3 \n\t"
                            "mulps %4, %%xmm4 \n\t"
                            "mulps %5, %%xmm5"
                            :
                            :
                            "m" ((*pl).c1.c1),
                            "m" ((*pl).c1.c3),
                            "m" ((*pl).c2.c2),
                            "m" ((*pl).c3.c1),
                            "m" ((*pl).c3.c3),
                            "m" ((*pl).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pk+=4;
      _prefetch_spinor(pk);
      pk-=4;
      
      __asm__ __volatile__ ("addps %%xmm0, %%xmm1 \n\t"
                            "addps %%xmm2, %%xmm3 \n\t"
                            "addps %%xmm4, %%xmm5 \n\t"
                            "addps %%xmm1, %%xmm3 \n\t"
                            "addps %%xmm3, %%xmm5 \n\t"
                            :
                            :
                            :
                            "xmm1", "xmm3", "xmm5");

      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %3, %%xmm3 \n\t"
                            "movaps %4, %%xmm4 \n\t"
                            "movaps %5, %%xmm6"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c3.c1),
                            "m" ((*pk).c3.c3),
                            "m" ((*pk).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");

      pl+=4;
      _prefetch_spinor(pl);
      pl-=4;
   
      __asm__ __volatile__ ("shufps $0xb1, %%xmm0, %%xmm0 \n\t"
                            "shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                            "shufps $0xb1, %%xmm2, %%xmm2 \n\t"
                            "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                            "shufps $0xb1, %%xmm4, %%xmm4 \n\t"
                            "shufps $0xb1, %%xmm6, %%xmm6"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");
      
      __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                            "mulps %1, %%xmm1 \n\t"
                            "mulps %2, %%xmm2 \n\t"
                            "mulps %3, %%xmm3 \n\t"
                            "mulps %4, %%xmm4 \n\t"
                            "mulps %5, %%xmm6"
                            :
                            :
                            "m" ((*pl).c1.c1),
                            "m" ((*pl).c1.c3),
                            "m" ((*pl).c2.c2),
                            "m" ((*pl).c3.c1),
                            "m" ((*pl).c3.c3),
                            "m" ((*pl).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm6");
      
      __asm__ __volatile__ ("addps %%xmm0, %%xmm1 \n\t"
                            "addps %%xmm2, %%xmm3 \n\t"
                            "addps %%xmm4, %%xmm6 \n\t"
                            "addps %%xmm1, %%xmm3 \n\t"
                            "addps %%xmm3, %%xmm6 \n\t"
                            "movaps %%xmm5, %0 \n\t"
                            "movaps %%xmm6, %1"                            
                            :
                            "=m" (tmp1),
                            "=m" (tmp2)
                            :
                            :
                            "xmm1", "xmm3", "xmm6");

      x+=(double)(tmp1.c1+tmp1.c2+tmp1.c3+tmp1.c4);      
      y+=(double)(tmp2.c2-tmp2.c1+tmp2.c4-tmp2.c3);
      pk+=1;
      pl+=1;
   }


   if ((icom==1)&&(NPROC>1))
   {
      v.re=x;
      v.im=y;

      MPI_Reduce(&v.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);     

      z.re=(float)(w.re);
      z.im=(float)(w.im);
   }
   else
   {
      z.re=(float)(x);
      z.im=(float)(y);
   }
   
   return z;  
}


float norm_square(int vol,int icom,spinor *pk)
{
   double x,y;
   spinor *pkm;
  
   x=0.0;
   pkm=pk+vol;
   
   for (;pk<pkm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %3, %%xmm3 \n\t"
                            "movaps %4, %%xmm4 \n\t"
                            "movaps %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c3.c1),
                            "m" ((*pk).c3.c3),
                            "m" ((*pk).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pk+=4;
      _prefetch_spinor(pk);
      pk-=3;
      
      __asm__ __volatile__ ("mulps %%xmm0, %%xmm0 \n\t"
                            "mulps %%xmm1, %%xmm1 \n\t"
                            "mulps %%xmm2, %%xmm2 \n\t"
                            "mulps %%xmm3, %%xmm3 \n\t"
                            "mulps %%xmm4, %%xmm4 \n\t"
                            "mulps %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
      
      __asm__ __volatile__ ("addps %%xmm0, %%xmm1 \n\t"
                            "addps %%xmm2, %%xmm3 \n\t"
                            "addps %%xmm4, %%xmm5 \n\t"
                            "addps %%xmm1, %%xmm3 \n\t"
                            "addps %%xmm3, %%xmm5 \n\t"
                            "movaps %%xmm5, %0"
                            :
                            "=m" (tmp1)
                            :
                            :
                            "xmm1", "xmm3", "xmm5");

      x+=(double)(tmp1.c1+tmp1.c2+tmp1.c3+tmp1.c4);
   }

   if ((icom==1)&&(NPROC>1))
   {
      MPI_Reduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return (float)(y);
   }
   else
      return (float)(x);
}


void mulc_spinor_add(int vol,spinor *pk,spinor *pl,complex z)
{
   spinor *pkm;

   _sse_load_cmplx(z);
   pkm=pk+vol;
   
   for (;pk<pkm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %%xmm0, %%xmm3 \n\t"
                            "movaps %%xmm1, %%xmm4 \n\t"
                            "movaps %%xmm2, %%xmm5 \n\t"
                            :
                            :
                            "m" ((*pl).c1.c1),
                            "m" ((*pl).c1.c3),
                            "m" ((*pl).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pl+=4;
      _prefetch_spinor(pl);
      pl-=4;
      
      __asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm6, %%xmm2 \n\t"
                            "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                            "shufps $0xb1, %%xmm4, %%xmm4 \n\t"
                            "shufps $0xb1, %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("addps %0, %%xmm0 \n\t"
                            "addps %1, %%xmm1 \n\t"
                            "addps %2, %%xmm2 \n\t"
                            "mulps %%xmm7, %%xmm3 \n\t"
                            "mulps %%xmm7, %%xmm4 \n\t"
                            "mulps %%xmm7, %%xmm5 \n\t"
                            "addps %%xmm0, %%xmm3 \n\t"
                            "addps %%xmm1, %%xmm4 \n\t"
                            "addps %%xmm2, %%xmm5"                            
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movaps %3, %%xmm0 \n\t"
                            "movaps %4, %%xmm1 \n\t"
                            "movaps %5, %%xmm2 \n\t"
                            "movaps %%xmm3, %0 \n\t"
                            "movaps %%xmm4, %1 \n\t"
                            "movaps %%xmm5, %2 \n\t"
                            "movaps %%xmm0, %%xmm3 \n\t"
                            "movaps %%xmm1, %%xmm4 \n\t"
                            "movaps %%xmm2, %%xmm5"
                            :
                            "=m" ((*pk).c1.c1),
                            "=m" ((*pk).c1.c3),
                            "=m" ((*pk).c2.c2)
                            :
                            "m" ((*pl).c3.c1),
                            "m" ((*pl).c3.c3),
                            "m" ((*pl).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pk+=4;
      _prefetch_spinor(pk);
      pk-=4;
      
      __asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm6, %%xmm2 \n\t"
                            "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                            "shufps $0xb1, %%xmm4, %%xmm4 \n\t"
                            "shufps $0xb1, %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("addps %0, %%xmm0 \n\t"
                            "addps %1, %%xmm1 \n\t"
                            "addps %2, %%xmm2 \n\t"
                            "mulps %%xmm7, %%xmm3 \n\t"
                            "mulps %%xmm7, %%xmm4 \n\t"
                            "mulps %%xmm7, %%xmm5 \n\t"
                            "addps %%xmm0, %%xmm3 \n\t"
                            "addps %%xmm1, %%xmm4 \n\t"
                            "addps %%xmm2, %%xmm5"                            
                            :
                            :
                            "m" ((*pk).c3.c1),
                            "m" ((*pk).c3.c3),
                            "m" ((*pk).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movaps %%xmm3, %0 \n\t"
                            "movaps %%xmm4, %1 \n\t"
                            "movaps %%xmm5, %2 \n\t"
                            :
                            "=m" ((*pk).c3.c1),
                            "=m" ((*pk).c3.c3),
                            "=m" ((*pk).c4.c2));

      pk+=1;
      pl+=1;
   }
}


float normalize(int vol,int icom,spinor *pk)
{
   float r,ri;
   spinor *pkm;   

   r=norm_square(vol,icom,pk);
   r=(float)(sqrt((double)(r)));

   if (error_loc(r==0.0f,1,"normalize [linalg.c]",
                 "Vector has vanishing norm")!=0)
      return 0.0f;

   ri=1.0f/r;

   __asm__ __volatile__ ("movss %0, %%xmm6 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm6 \n\t"
                         "movaps %%xmm6, %%xmm7"
                         :
                         :
                         "m" (ri)
                         :
                         "xmm6", "xmm7");

   pkm=pk+vol;
   
   for (;pk<pkm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %3, %%xmm3 \n\t"
                            "movaps %4, %%xmm4 \n\t"
                            "movaps %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c3.c1),
                            "m" ((*pk).c3.c3),
                            "m" ((*pk).c4.c2)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      pk+=4;
      _prefetch_spinor(pk);
      pk-=4;      
      
      __asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t"
                            "mulps %%xmm7, %%xmm1 \n\t"
                            "mulps %%xmm6, %%xmm2 \n\t"
                            "mulps %%xmm7, %%xmm3 \n\t"
                            "mulps %%xmm6, %%xmm4 \n\t"
                            "mulps %%xmm7, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
      
      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm3, %3 \n\t"
                            "movaps %%xmm4, %4 \n\t"
                            "movaps %%xmm5, %5"
                            :
                            "=m" ((*pk).c1.c1),
                            "=m" ((*pk).c1.c3),
                            "=m" ((*pk).c2.c2),
                            "=m" ((*pk).c3.c1),
                            "=m" ((*pk).c3.c3),
                            "=m" ((*pk).c4.c2));

      pk+=1;
   }

   return r;
}


void mulg5(int vol,spinor *pk)
{
   weyl *rpk,*rpm;
   
   __asm__ __volatile__ ("movaps %0, %%xmm5 \n\t"
                         "movaps %%xmm5, %%xmm6 \n\t"
                         "movaps %%xmm5, %%xmm7"
                         :
                         :
                         "m" (_sse_sgn)
                         :
                         "xmm5", "xmm6", "xmm7");                         

   rpk=(weyl*)(pk);
   rpm=rpk+2*vol;
   rpk+=1;   
   
   for (;rpk<rpm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2"
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c3),
                            "m" ((*rpk).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2");

      rpk+=8;
      _prefetch_weyl(rpk);
      rpk-=8;          
      
      __asm__ __volatile__ ("mulps %%xmm5, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm7, %%xmm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2\t"
                            :
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c3),
                            "=m" ((*rpk).c2.c2));

      rpk+=2;
   }
}


void mulmg5(int vol,spinor *pk)
{
   weyl *rpk,*rpm;
   
   __asm__ __volatile__ ("movaps %0, %%xmm5 \n\t"
                         "movaps %%xmm5, %%xmm6 \n\t"
                         "movaps %%xmm5, %%xmm7"
                         :
                         :
                         "m" (_sse_sgn)
                         :
                         "xmm5", "xmm6", "xmm7");                         

   rpk=(weyl*)(pk);
   rpm=rpk+2*vol;
   
   for (;rpk<rpm;)
   {
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %1, %%xmm1 \n\t"
                            "movaps %2, %%xmm2"
                            :
                            :
                            "m" ((*rpk).c1.c1),
                            "m" ((*rpk).c1.c3),
                            "m" ((*rpk).c2.c2)
                            :
                            "xmm0", "xmm1", "xmm2");

      rpk+=8;
      _prefetch_weyl(rpk);
      rpk-=8;          
      
      __asm__ __volatile__ ("mulps %%xmm5, %%xmm0 \n\t"
                            "mulps %%xmm6, %%xmm1 \n\t"
                            "mulps %%xmm7, %%xmm2"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2\t"
                            :
                            "=m" ((*rpk).c1.c1),
                            "=m" ((*rpk).c1.c3),
                            "=m" ((*rpk).c2.c2));

      rpk+=2;          
   }
}

#else

complex spinor_prod(int vol,int icom,spinor *pk,spinor *pl)
{
   complex z;
   double x,y;
   complex_dble v,w;
   spinor *rpk,*rpl;
  
   x=0.0;
   y=0.0;
   rpk=pk;
   rpl=pl;

   for (;rpk<(pk+vol);)
   {
      x+=(double)(_vector_prod_re((*rpk).c1,(*rpl).c1)+
            _vector_prod_re((*rpk).c2,(*rpl).c2)+
            _vector_prod_re((*rpk).c3,(*rpl).c3)+
            _vector_prod_re((*rpk).c4,(*rpl).c4));

      y+=(double)(_vector_prod_im((*rpk).c1,(*rpl).c1)+
            _vector_prod_im((*rpk).c2,(*rpl).c2)+
            _vector_prod_im((*rpk).c3,(*rpl).c3)+
            _vector_prod_im((*rpk).c4,(*rpl).c4));

      rpk+=1;
      rpl+=1;
   }

   if ((icom!=1)||(NPROC==1))
   {
      z.re=(float)(x);
      z.im=(float)(y);
   }
   else
   {
      v.re=x;
      v.im=y;

      MPI_Reduce(&v.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);     

      z.re=(float)(w.re);
      z.im=(float)(w.im);
   }
   
   return z;  
}


float norm_square(int vol,int icom,spinor *pk)
{
   double x,y;
   spinor *rpk;

   x=0.0;
   rpk=pk;
 
   for (;rpk<(pk+vol);)
   {
      x+=(double)(_vector_prod_re((*rpk).c1,(*rpk).c1)+
            _vector_prod_re((*rpk).c2,(*rpk).c2)+
            _vector_prod_re((*rpk).c3,(*rpk).c3)+
            _vector_prod_re((*rpk).c4,(*rpk).c4));

      rpk+=1;
   }

   if ((icom!=1)||(NPROC==1))
      return (float)(x);
   else
   {
      MPI_Reduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return (float)(y);
   }
}


void mulc_spinor_add(int vol,spinor *pk,spinor *pl,complex z)
{
   spinor *rpk,*rpl;

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


float normalize(int vol,int icom,spinor *pk)
{
   float r,ri;
   spinor *rpk;

   r=norm_square(vol,icom,pk);
   r=(float)(sqrt((double)(r)));

   if (error_loc(r==0.0f,1,
       "normalize [linalg.c]","Vector has vanishing norm")!=0)
      return 0.0f;

   ri=1.0f/r;
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


void mulg5(int vol,spinor *pk)
{
   weyl *rpk,*rpm;

   rpk=(weyl*)(pk);
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


void mulmg5(int vol,spinor *pk)
{
   weyl *rpk,*rpm;

   rpk=(weyl*)(pk);
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
