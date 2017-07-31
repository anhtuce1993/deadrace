
/*******************************************************************************
*
* File sinit.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic assignment and initialization programs for single- and double-
* precision spinor fields
*
* All these programs operate on arrays of spinor fields, whose base address
* is passed through the arguments. The length of the arrays is specified by
* the parameter vol.
*
* The externally accessible functions are
*
*   void set_s2zero(int vol,spinor *pk)
*     Sets the single-precision spinor field pk[] to zero
*
*   void set_sd2zero(int vol,spinor_dble *pk)
*     Sets the double-precision spinor field pk[] to zero
*
*   void random_s(int vol,spinor *pk,float sigma)
*     Initializes the components of the single-precision field pk[]
*     to (complex) random values z with distribution proportional
*     to exp{-|z|^2/sigma^2}
*
*   void random_sd(int vol,spinor_dble *pk,double sigma)
*     Initializes the components of the double-precision field pk[]
*     to (complex) random values z with distribution proportional
*     to exp{-|z|^2/sigma^2}
*
*   void assign_s2s(int vol,spinor *pk,spinor *pl)
*     Assigns the single-precision field pk[] to the single-precision
*     field pl[]
*
*   void assign_s2sd(int vol,spinor *pk,spinor_dble *pl)
*     Assigns the single-precision field pk[] to the double-precision
*     field pl[]
*
*   void assign_sd2s(int vol,spinor_dble *pk,spinor *pl)
*     Assigns the double-precision field pk[] to the single-precision
*     field pl[]
*
*   void assign_sd2sd(int vol,spinor_dble *pk,spinor_dble *pl)
*     Assigns the double-precision field pk[] to the double-precision
*     field pl[]
*
*******************************************************************************/

#define SINIT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "start.h"
#include "global.h"


#if (defined SSE)
#include "sse.h"

void set_s2zero(int vol,spinor *pk)
{
   spinor *pkm;

   __asm__ __volatile__ ("xorps %%xmm0, %%xmm0 \n\t"
                         "xorps %%xmm1, %%xmm1 \n\t"
                         "xorps %%xmm2, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   
   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm0, %3 \n\t"
                            "movaps %%xmm1, %4 \n\t"
                            "movaps %%xmm2, %5"                            
                            :
                            "=m" ((*pk).c1.c1),
                            "=m" ((*pk).c1.c3),
                            "=m" ((*pk).c2.c2),
                            "=m" ((*pk).c3.c1),
                            "=m" ((*pk).c3.c3),
                            "=m" ((*pk).c4.c2));
   }
}


void assign_s2s(int vol,spinor *pk,spinor *pl)
{
   spinor *pkm;

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

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %1 \n\t"
                            "movaps %%xmm2, %2 \n\t"
                            "movaps %%xmm3, %3 \n\t"
                            "movaps %%xmm4, %4 \n\t"
                            "movaps %%xmm5, %5"                            
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c3.c1),
                            "=m" ((*pl).c3.c3),
                            "=m" ((*pl).c4.c2));

      pl+=1;
   }
}
            
#else

static const spinor s0={{{0.0f}}};

void set_s2zero(int vol,spinor *pk)
{
   spinor *rpk;

   for (rpk=pk;rpk<(pk+vol);rpk++)
      *rpk=s0;
}


void assign_s2s(int vol,spinor *pk,spinor *pl)
{
   spinor *rpk,*rpl;

   rpl=pl;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {
      *rpl=*rpk;
      rpl+=1;
   }
}

#endif


#if (defined SSE2)
#include "sse2.h"

void set_sd2zero(int vol,spinor_dble *pk)
{
   spinor_dble *pkm;

   __asm__ __volatile__ ("xorpd %%xmm0, %%xmm0 \n\t"
                         "xorpd %%xmm1, %%xmm1 \n\t"
                         "xorpd %%xmm2, %%xmm2"
                         :
                         :
                         :
                         "xmm0", "xmm1", "xmm2");
   
   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm0, %3 \n\t"
                            "movapd %%xmm1, %4 \n\t"
                            "movapd %%xmm2, %5"                            
                            :
                            "=m" ((*pk).c1.c1),
                            "=m" ((*pk).c1.c2),
                            "=m" ((*pk).c1.c3),
                            "=m" ((*pk).c2.c1),
                            "=m" ((*pk).c2.c2),
                            "=m" ((*pk).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm0, %3 \n\t"
                            "movapd %%xmm1, %4 \n\t"
                            "movapd %%xmm2, %5"                            
                            :
                            "=m" ((*pk).c3.c1),
                            "=m" ((*pk).c3.c2),
                            "=m" ((*pk).c3.c3),
                            "=m" ((*pk).c4.c1),
                            "=m" ((*pk).c4.c2),
                            "=m" ((*pk).c4.c3));
   }
}

#else

static const spinor_dble sd0={{{0.0}}};

void set_sd2zero(int vol,spinor_dble *pk)
{
   spinor_dble *rpk;

   for (rpk=pk;rpk<(pk+vol);rpk++)
      *rpk=sd0;
}

#endif

void random_s(int vol,spinor *pk,float sigma)
{
   int i;
   float r[24],*s;
   spinor *rpk;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {
      gauss(r,24);
      s=(float*)(rpk);

      for (i=0;i<24;i++)
         s[i]=sigma*r[i];
   }
}


void random_sd(int vol,spinor_dble *pk,double sigma)
{
   int i;
   double r[24],*s;
   spinor_dble *rpk;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {
      gauss_dble(r,24);
      s=(double*)(rpk);

      for (i=0;i<24;i++)
         s[i]=sigma*r[i];
   }
}


void assign_s2sd(int vol,spinor *pk,spinor_dble *pl)
{
   int i;
   float *r;
   double *rd;
   spinor *rpk;
   spinor_dble *rpl;

   rpl=pl;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {   
      r=(float*)(rpk);
      rd=(double*)(rpl);

      for (i=0;i<24;i++)
         rd[i]=(double)(r[i]);

      rpl+=1;
   }
}


void assign_sd2s(int vol,spinor_dble *pk,spinor *pl)
{
   int i;
   float *r;
   double *rd;
   spinor *rpl;
   spinor_dble *rpk;

   rpl=pl;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {   
      rd=(double*)(rpk);
      r=(float*)(rpl);

      for (i=0;i<24;i++)
         r[i]=(float)(rd[i]);

      rpl+=1;
   }
}


void assign_sd2sd(int vol,spinor_dble *pk,spinor_dble *pl)
{
   spinor_dble *rpk,*rpl;

   rpl=pl;
   
   for (rpk=pk;rpk<(pk+vol);rpk++)
   {
      *rpl=*rpk;
      rpl+=1;
   }
}

