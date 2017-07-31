
/*******************************************************************************
*
* File cmatrix.c
*
* Copyright (C) 2007 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Complex matrix algebra (single-precision version)
*
* The externally accessible functions are
*
*   void cmat_vec(int n,complex *a,complex *v,complex *w)
*     Computes w=a*v, where v and w are n-vectors and a an nxn matrix
*
*   void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
*     Adds a*v to w, where v and w are n-vectors and a an nxn matrix
*     
*   void cmat_add(int n,complex *a,complex *b,complex *c)
*     Computes the sum c=a+b of two nxn matrices a and b
*
*   void cmat_sub(int n,complex *a,complex *b,complex *c)
*     Computes the difference c=a-b of two nxn matrices a and b
*
*   void cmat_mul(int n,complex *a,complex *b,complex *c)
*     Computes the product c=a*b of two nxn matrices a and b
*
*   void cmat_dag(int n,complex *a,complex *b)
*     Assigns the hermitian conjugate of a to b
*
* Notes:
*
* All of these programs can be called locally. Complex nxn matrices with
* matrix elements A_{ij} are represented by linear arrays a of complex
* numbers such that
*
*   A_{ij} = a[i*n+j]
*
* where i,j=0,1,..,n-1. It is assumed that the input and output arrays do 
* not overlap in memory (the results are otherwise unpredictable)
*
* If SSE instructions are to be used, and if n is even, it is taken for
* granted that the starting addresses of the arrays are multiples of 16
* (SSE instructions are not used if n is odd)
*
*******************************************************************************/

#define CMATRIX_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"
#include "start.h"

#if (defined SSE)
#include "sse.h"

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *b,*vv,*vm,*wm;

   if ((n&0x1)==0x0)
   {
      b=a;
      vm=v+n;
      wm=w+n;
      
      for (;w<wm;w+=2)
      {
         a=b;
         b+=n;
      
         __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                               "xorps %%xmm1, %%xmm1 \n\t"
                               "xorps %%xmm2, %%xmm2 \n\t"
                               "xorps %%xmm3, %%xmm3"
                               :
                               :
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
         
         for (vv=v;vv<vm;vv+=2)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %1, %%xmm5 \n\t"
                                  "movaps %%xmm4, %%xmm6 \n\t"
                                  "movaps %%xmm5, %%xmm7 \n\t"
                                  "shufps $0xb1, %%xmm6, %%xmm6 \n\t"
                                  "shufps $0xb1, %%xmm7, %%xmm7 \n\t"
                                  "mulps %2, %%xmm4 \n\t"
                                  "mulps %2, %%xmm5 \n\t"
                                  "mulps %2, %%xmm6 \n\t"
                                  "mulps %2, %%xmm7 \n\t"
                                  "addps %%xmm4, %%xmm0 \n\t"
                                  "addps %%xmm5, %%xmm1 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3"
                                  :
                                  :
                                  "m" ((*a).re),
                                  "m" ((*b).re),
                                  "m" ((*vv).re)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3", 
                                  "xmm4", "xmm5", "xmm6", "xmm7");         

            a+=2;
            b+=2;
         }

         __asm__ __volatile__ ("mulps %1, %%xmm0 \n\t"
                               "mulps %1, %%xmm1 \n\t"
                               "movaps %%xmm0, %%xmm4 \n\t"
                               "movaps %%xmm1, %%xmm5 \n\t"
                               "movlhps %%xmm2, %%xmm0 \n\t"
                               "movlhps %%xmm3, %%xmm1 \n\t"
                               "movhlps %%xmm4, %%xmm2 \n\t"
                               "movhlps %%xmm5, %%xmm3 \n\t"
                               "addps %%xmm2, %%xmm0 \n\t"
                               "addps %%xmm3, %%xmm1 \n\t"
                               "movaps %%xmm0, %%xmm4 \n\t"
                               "movaps %%xmm1, %%xmm5 \n\t"
                               "shufps $0x88, %%xmm1, %%xmm0 \n\t"
                               "shufps $0xdd, %%xmm5, %%xmm4 \n\t"
                               "addps %%xmm4, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" ((*w).re)
                               :
                               "m" (_sse_sgn24)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", 
                               "xmm4", "xmm5");
      }
   }
   else
   {
      vm=v+n;
      wm=w+n;
   
      for (;w<wm;w++)
      {
         (*w).re=0.0f;
         (*w).im=0.0f;
         
         for (vv=v;vv<vm;vv++)
         {
            (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
            (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
            a+=1;
         }
      }      
   }
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *b,*vv,*vm,*wm;

   if ((n&0x1)==0x0)
   {
      b=a;
      vm=v+n;
      wm=w+n;
      
      for (;w<wm;w+=2)
      {
         a=b;
         b+=n;
      
         __asm__ __volatile__ ("movss %0, %%xmm0 \n\t"
                               "movss %1, %%xmm1 \n\t"
                               "movss %2, %%xmm2 \n\t"
                               "movss %3, %%xmm3"
                               :
                               :
                               "m" ((*(w  )).re),
                               "m" ((*(w+1)).re),
                               "m" ((*(w  )).im),
                               "m" ((*(w+1)).im)                               
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3");
         
         for (vv=v;vv<vm;vv+=2)
         {
            __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                                  "movaps %1, %%xmm5 \n\t"
                                  "movaps %%xmm4, %%xmm6 \n\t"
                                  "movaps %%xmm5, %%xmm7 \n\t"
                                  "shufps $0xb1, %%xmm6, %%xmm6 \n\t"
                                  "shufps $0xb1, %%xmm7, %%xmm7 \n\t"
                                  "mulps %2, %%xmm4 \n\t"
                                  "mulps %2, %%xmm5 \n\t"
                                  "mulps %2, %%xmm6 \n\t"
                                  "mulps %2, %%xmm7 \n\t"
                                  "addps %%xmm4, %%xmm0 \n\t"
                                  "addps %%xmm5, %%xmm1 \n\t"
                                  "addps %%xmm6, %%xmm2 \n\t"
                                  "addps %%xmm7, %%xmm3"
                                  :
                                  :
                                  "m" ((*a).re),
                                  "m" ((*b).re),
                                  "m" ((*vv).re)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3", 
                                  "xmm4", "xmm5", "xmm6", "xmm7");         

            a+=2;
            b+=2;
         }

         __asm__ __volatile__ ("mulps %1, %%xmm0 \n\t"
                               "mulps %1, %%xmm1 \n\t"
                               "movaps %%xmm0, %%xmm4 \n\t"
                               "movaps %%xmm1, %%xmm5 \n\t"
                               "movlhps %%xmm2, %%xmm0 \n\t"
                               "movlhps %%xmm3, %%xmm1 \n\t"
                               "movhlps %%xmm4, %%xmm2 \n\t"
                               "movhlps %%xmm5, %%xmm3 \n\t"
                               "addps %%xmm2, %%xmm0 \n\t"
                               "addps %%xmm3, %%xmm1 \n\t"
                               "movaps %%xmm0, %%xmm4 \n\t"
                               "movaps %%xmm1, %%xmm5 \n\t"
                               "shufps $0x88, %%xmm1, %%xmm0 \n\t"
                               "shufps $0xdd, %%xmm5, %%xmm4 \n\t"
                               "addps %%xmm4, %%xmm0 \n\t"
                               "movaps %%xmm0, %0"
                               :
                               "=m" ((*w).re)
                               :
                               "m" (_sse_sgn24)
                               :
                               "xmm0", "xmm1", "xmm2", "xmm3", 
                               "xmm4", "xmm5");
      }
   }
   else
   {
      vm=v+n;
      wm=w+n;
   
      for (;w<wm;w++)
      {
         for (vv=v;vv<vm;vv++)
         {
            (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
            (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
            a+=1;
         }
      }      
   }
}

#else

void cmat_vec(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;;

   vm=v+n;
   wm=w+n;
   
   for (;w<wm;w++)
   {
      (*w).re=0.0f;
      (*w).im=0.0f;
         
      for (vv=v;vv<vm;vv++)
      {
         (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
         (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
         a+=1;
      }
   }
}


void cmat_vec_assign(int n,complex *a,complex *v,complex *w)
{
   complex *vv,*vm,*wm;;

   vm=v+n;
   wm=w+n;
   
   for (;w<wm;w++)
   {
      for (vv=v;vv<vm;vv++)
      {
         (*w).re+=((*a).re*(*vv).re-(*a).im*(*vv).im);
         (*w).im+=((*a).re*(*vv).im+(*a).im*(*vv).re);
         a+=1;
      }
   }
}

#endif

void cmat_add(int n,complex *a,complex *b,complex *c)
{
   complex *am;

   am=a+n*n;

   for (;a<am;a++)
   {
      (*c).re=(*a).re+(*b).re;
      (*c).im=(*a).im+(*b).im;
      b+=1;
      c+=1;
   }
}


void cmat_sub(int n,complex *a,complex *b,complex *c)
{
   complex *am;

   am=a+n*n;

   for (;a<am;a++)
   {
      (*c).re=(*a).re-(*b).re;
      (*c).im=(*a).im-(*b).im;
      b+=1;
      c+=1;
   }
}


void cmat_mul(int n,complex *a,complex *b,complex *c)
{
   complex *aa,*bb,*am,*bm,*bbm;

   am=a+n*n;
   bm=b+n;
   bbm=b+n*n;
   
   for (;a<am;a+=n)
   {
      for (;b<bm;b++)
      {
         (*c).re=0.0f;
         (*c).im=0.0f;
         aa=a;

         for (bb=b;bb<bbm;bb+=n)
         {
            (*c).re+=((*aa).re*(*bb).re-(*aa).im*(*bb).im);
            (*c).im+=((*aa).re*(*bb).im+(*aa).im*(*bb).re);
            aa+=1;
         }

         c+=1;
      }

      b-=n;
   }
}


void cmat_dag(int n,complex *a,complex *b)
{
   complex *bb,*am,*bbm;

   am=a+n*n;
   bbm=b+n*n;
   
   for (;a<am;)
   {
      for (bb=b;bb<bbm;bb+=n)
      {
         (*bb).re=(*a).re;
         (*bb).im=-(*a).im;
         a+=1;
      }

      b+=1;
   }
}

