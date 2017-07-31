
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* SSE complex multiplication
*
*******************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "start.h"

#if (defined SSE)

#include "sse.h"

int main(void)
{
   complex *c;

   c=amalloc(8*sizeof(complex),4);
   error(c==NULL,1,"main [check1.c]","Unable to allocate auxiliary array");

   rlxs_init(0,3898);
   ranlxs((float*)(c),16);

   c[4].re=c[0].re*c[2].re-c[0].im*c[2].im;
   c[4].im=c[0].re*c[2].im+c[0].im*c[2].re;

   c[5].re=c[1].re*c[3].re-c[1].im*c[3].im;
   c[5].im=c[1].re*c[3].im+c[1].im*c[3].re;

   __asm__ __volatile__ ("movaps %1, %%xmm0 \n\t"
                         "movaps %1, %%xmm1 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xf5, %%xmm1, %%xmm1 \n\t"
                         "mulps %2, %%xmm0 \n\t"
                         "mulps %2, %%xmm1 \n\t"
                         "shufps $0xb1, %%xmm1, %%xmm1 \n\t"
                         "mulps %3, %%xmm1 \n\t"
                         "addps %%xmm0, %%xmm1 \n\t"
                         "movaps %%xmm1, %0"
                         :
                         "=m" (c[6].re)
                         :
                         "m" (c[0].re),
                         "m" (c[2].re),
                         "m" (_sse_sgn13));

   printf("\n");
   printf("SSE complex multiplication\n\n");

   printf("FPU Result: "); 
   printf("(%.7e,%.7e),(%.7e,%.7e)\n",c[4].re,c[4].im,c[5].re,c[5].im);
   printf("SSE Result: "); 
   printf("(%.7e,%.7e),(%.7e,%.7e)\n",c[6].re,c[6].im,c[7].re,c[7].im);   
   printf("\n");

   exit(0);
}

#else

int main(void)
{
   printf("\n");
   printf("SSE complex multiplication\n\n");

   printf("SSE must be enabled for this test\n\n");
   exit(0);
}

#endif
