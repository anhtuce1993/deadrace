
/*******************************************************************************
*
* File pauli.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic functions for single-precision hermitian 6x6 matrices
*
* The externally accessible functions are
*
*   void mul_pauli(pauli *m,weyl *s,weyl *r)
*     Multiplies the Weyl spinor s by the matrix m and assigns the result
*     to the Weyl spinor r. The source spinor is overwritten if r=s and
*     otherwise left unchanged
*
*   void assign_pauli(int vol,pauli_dble *md,pauli *m)
*     Assigns a field of matrices of type pauli_dble to a field of type
*     pauli. The parameters are the number vol of elements in the fields
*     and the starting addresses md and m
*
*   void apply_sw(int vol,pauli *m,spinor *s,spinor *r)
*     Applies the pauli matrix field m[2*vol] to the spinor field s[vol]
*     and assigns the result to the field r[vol]. The source field is
*     overwritten if r=s and otherwise left unchanged (the arrays may
*     not overlap in this case)
*
* Notes:
*
* The storage format for hermitian 6x6 matrices is described in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.ps)
*
* If SSE/SSE2/SSE3 instructions are used, it is assumed that the matrices
* and spinors are aligned to a 16 byte boundary
*
*******************************************************************************/

#define PAULI_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"

#if (defined SSE)
#include "sse.h"

static sse_float rs1,rs2,rs3,rs4;

#endif

#if (defined SSE3)

void mul_pauli(pauli *m,weyl *s,weyl *r)
{

/*************************** components 1 and 2 *******************************/

   __asm__ __volatile__ ("movss %2, %%xmm4 \n\t"
                         "movss %3, %%xmm5 \n\t"
                         "movss %4, %%xmm6 \n\t"
                         "movss %5, %%xmm7 \n\t"
                         "movsldup %0, %%xmm0 \n\t"
                         "movsldup %1, %%xmm1 \n\t"
                         "movshdup %0, %%xmm2 \n\t"
                         "movshdup %1, %%xmm3 \n\t"
                         "shufps $0x0, %%xmm5, %%xmm4 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm5 \n\t"
                         "shufps $0x5, %%xmm7, %%xmm6 \n\t"
                         "shufps $0x50, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[0]),
                         "m" ((*m).u[6]),
                         "m" ((*m).u[1]),
                         "m" ((*m).u[7])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                         "mulps %0, %%xmm1 \n\t"
                         "mulps %0, %%xmm2 \n\t"
                         "mulps %0, %%xmm3 \n\t"
                         "mulps %1, %%xmm4 \n\t"
                         "mulps %1, %%xmm5 \n\t"
                         "mulps %1, %%xmm6 \n\t"
                         "mulps %1, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "subps %%xmm7, %%xmm3"
                         :
                         :
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c1.c1.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movsldup %0, %%xmm4 \n\t"
                         "movsldup %1, %%xmm5 \n\t"
                         "movshdup %0, %%xmm6 \n\t"
                         "movshdup %1, %%xmm7 \n\t"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[20])
                         :
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %2, %%xmm4 \n\t"
                         "mulps %2, %%xmm5 \n\t"
                         "mulps %2, %%xmm6 \n\t"
                         "mulps %2, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                         "addsubps %%xmm2, %%xmm0 \n\t"
                         "addsubps %%xmm3, %%xmm1 \n\t"
                         "movaps %%xmm0, %0 \n\t"
                         "movaps %%xmm1, %1"
                         :
                         "=m" (rs1),
                         "=m" (rs2)
                         :
                         "m" ((*s).c2.c2.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

/*************************** components 3 and 4 *******************************/

   __asm__ __volatile__ ("movlps %0, %%xmm4 \n\t"
                         "movlps %1, %%xmm5 \n\t"
                         "movss %2, %%xmm0 \n\t"
                         "movss %3, %%xmm1 \n\t"
                         "movhps %4, %%xmm4 \n\t"
                         "movhps %5, %%xmm5"
                         :
                         :
                         "m" ((*m).u[26]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[2]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[32])
                         :
                         "xmm0", "xmm1", "xmm4", "xmm5");
                         
   __asm__ __volatile__ ("movss %0, %%xmm2 \n\t"
                         "movss %1, %%xmm3 \n\t"
                         "movshdup %%xmm4, %%xmm6 \n\t"
                         "movshdup %%xmm5, %%xmm7 \n\t"
                         "shufps $0x0, %%xmm1, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm4, %%xmm4 \n\t"
                         "shufps $0x0, %%xmm2, %%xmm1 \n\t"
                         "shufps $0xa0, %%xmm5, %%xmm5 \n\t"
                         "shufps $0x5, %%xmm3, %%xmm2 \n\t"
                         "shufps $0x50, %%xmm3, %%xmm3"
                         :
                         :
                         "m" ((*m).u[3]),
                         "m" ((*m).u[25])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm4 \n\t"
                         "mulps %0, %%xmm5 \n\t"
                         "mulps %0, %%xmm6 \n\t"
                         "mulps %0, %%xmm7 \n\t"
                         "mulps %1, %%xmm0 \n\t"
                         "mulps %1, %%xmm1 \n\t"
                         "mulps %1, %%xmm2 \n\t"
                         "mulps %1, %%xmm3 \n\t"
                         "addps %%xmm0, %%xmm4 \n\t"
                         "addps %%xmm1, %%xmm5 \n\t"
                         "addps %%xmm2, %%xmm6 \n\t"
                         "subps %%xmm3, %%xmm7"
                         :
                         :
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c1.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movlps %0, %%xmm0 \n\t"
                         "movlps %1, %%xmm1 \n\t"
                         "movhps %2, %%xmm0 \n\t"
                         "movhps %3, %%xmm1 \n\t"
                         "movshdup %%xmm0, %%xmm2 \n\t"
                         "movshdup %%xmm1, %%xmm3 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm1, %%xmm1"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[10]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[18])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3");

   __asm__ __volatile__ ("mulps %2, %%xmm0 \n\t"
                         "mulps %2, %%xmm1 \n\t"
                         "mulps %2, %%xmm2 \n\t"
                         "mulps %2, %%xmm3 \n\t"
                         "addps %%xmm0, %%xmm4 \n\t"
                         "addps %%xmm1, %%xmm5 \n\t"
                         "subps %%xmm2, %%xmm6 \n\t"
                         "subps %%xmm3, %%xmm7 \n\t"
                         "shufps $0xb1, %%xmm6, %%xmm6 \n\t"
                         "shufps $0xb1, %%xmm7, %%xmm7 \n\t"
                         "addsubps %%xmm6, %%xmm4 \n\t"
                         "addsubps %%xmm7, %%xmm5 \n\t"
                         "movaps %%xmm4, %0 \n\t"
                         "movaps %%xmm5, %1"
                         :
                         "=m" (rs3),
                         "=m" (rs4)
                         :
                         "m" ((*s).c1.c1.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");
   
/*************************** components 5 and 6 *******************************/

   __asm__ __volatile__ ("movlps %0, %%xmm0 \n\t"
                         "movlps %1, %%xmm1 \n\t"
                         "movss %2, %%xmm4 \n\t"
                         "movss %3, %%xmm5 \n\t"
                         "movhps %4, %%xmm0 \n\t"
                         "movhps %5, %%xmm1"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[14]),
                         "m" ((*m).u[4]),
                         "m" ((*m).u[34]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[22])                         
                         :
                         "xmm0", "xmm1", "xmm4", "xmm5");
                         
   __asm__ __volatile__ ("movss %0, %%xmm6 \n\t"
                         "movss %1, %%xmm7 \n\t"
                         "movshdup %%xmm0, %%xmm2 \n\t"
                         "movshdup %%xmm1, %%xmm3 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm1, %%xmm1 \n\t"
                         "shufps $0x0, %%xmm5, %%xmm4 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm5 \n\t"
                         "shufps $0x5, %%xmm7, %%xmm6 \n\t"
                         "shufps $0x50, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[5]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                         "mulps %0, %%xmm1 \n\t"
                         "mulps %0, %%xmm2 \n\t"
                         "mulps %0, %%xmm3 \n\t"
                         "mulps %1, %%xmm4 \n\t"
                         "mulps %1, %%xmm5 \n\t"
                         "mulps %1, %%xmm6 \n\t"
                         "mulps %1, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "subps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c2.c2.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movlps %0, %%xmm4 \n\t"
                         "movlps %1, %%xmm5 \n\t"
                         "movhps %2, %%xmm4 \n\t"
                         "movhps %3, %%xmm5 \n\t"                         
                         "movshdup %%xmm4, %%xmm6 \n\t"
                         "movshdup %%xmm5, %%xmm7 \n\t"
                         "shufps $0xa0, %%xmm4, %%xmm4 \n\t"
                         "shufps $0xa0, %%xmm5, %%xmm5"
                         :
                         :
                         "m" ((*m).u[26]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[30]),                         
                         "m" ((*m).u[32])
                         :
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm4 \n\t"
                         "mulps %0, %%xmm5 \n\t"
                         "mulps %0, %%xmm6 \n\t"
                         "mulps %0, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                         "mulps %1, %%xmm2 \n\t"
                         "mulps %1, %%xmm3 \n\t"                         
                         "addps %%xmm2, %%xmm0 \n\t"
                         "addps %%xmm3, %%xmm1"
                         :
                         :
                         "m" ((*s).c1.c3.re),
                         "m" (_sse_sgn24)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

/*************************** combine and store ********************************/
   
   __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                         "movaps %1, %%xmm5 \n\t"
                         "movlhps %%xmm1, %%xmm2 \n\t"
                         "movlhps %%xmm3, %%xmm4 \n\t"
                         "movlhps %%xmm5, %%xmm6 \n\t"
                         "movhlps %%xmm0, %%xmm1 \n\t"
                         "shufps $0xe4, %%xmm2, %%xmm0 \n\t"
                         "movlps %2, %%xmm3 \n\t"
                         "movlps %3, %%xmm4 \n\t"
                         "movlps %4, %%xmm5 \n\t"
                         "movlps %5, %%xmm6 \n\t"
                         "addps %%xmm3, %%xmm4 \n\t"
                         "addps %%xmm5, %%xmm6 \n\t"
                         "addps %%xmm1, %%xmm0"
                         :
                         :
                         "m" (rs2.c1),
                         "m" (rs4.c1),
                         "m" (rs1.c1),
                         "m" (rs1.c3),
                         "m" (rs3.c1),                         
                         "m" (rs3.c3)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6");

   __asm__ __volatile__ ("movaps %%xmm4, %0 \n\t"
                         "movaps %%xmm6, %1 \n\t"
                         "movaps %%xmm0, %2"
                         :
                         "=m" ((*r).c1.c1),
                         "=m" ((*r).c1.c3),
                         "=m" ((*r).c2.c2));
}

#elif (defined SSE)

void mul_pauli(pauli *m,weyl *s,weyl *r)
{

/*************************** components 1 and 2 *******************************/

   __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                         "movaps %1, %%xmm1 \n\t"
                         "movaps %%xmm0, %%xmm2 \n\t"
                         "movaps %%xmm1, %%xmm3 \n\t"
                         "movss %2, %%xmm4 \n\t"
                         "movss %3, %%xmm5 \n\t"
                         "movss %4, %%xmm6 \n\t"
                         "movss %5, %%xmm7 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm1, %%xmm1 \n\t"
                         "shufps $0xf5, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xf5, %%xmm3, %%xmm3 \n\t"
                         "shufps $0x0, %%xmm5, %%xmm4 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm5 \n\t"
                         "shufps $0x5, %%xmm7, %%xmm6 \n\t"
                         "shufps $0x50, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[0]),
                         "m" ((*m).u[6]),
                         "m" ((*m).u[1]),
                         "m" ((*m).u[7])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                         "mulps %0, %%xmm1 \n\t"
                         "mulps %0, %%xmm2 \n\t"
                         "mulps %0, %%xmm3 \n\t"
                         "mulps %1, %%xmm4 \n\t"
                         "mulps %1, %%xmm5 \n\t"
                         "mulps %1, %%xmm6 \n\t"
                         "mulps %1, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "subps %%xmm7, %%xmm3"
                         :
                         :
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c1.c1.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movaps %0, %%xmm4 \n\t"
                         "movaps %1, %%xmm5 \n\t"
                         "movaps %%xmm4, %%xmm6 \n\t"
                         "movaps %%xmm5, %%xmm7 \n\t"
                         "shufps $0xa0, %%xmm4, %%xmm4 \n\t"
                         "shufps $0xa0, %%xmm5, %%xmm5 \n\t"
                         "shufps $0xf5, %%xmm6, %%xmm6 \n\t"
                         "shufps $0xf5, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[20])
                         :
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %2, %%xmm4 \n\t"
                         "mulps %2, %%xmm5 \n\t"
                         "mulps %2, %%xmm6 \n\t"
                         "mulps %2, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                         "mulps %3, %%xmm2 \n\t"
                         "mulps %3, %%xmm3 \n\t"                         
                         "addps %%xmm2, %%xmm0 \n\t"
                         "addps %%xmm3, %%xmm1 \n\t"
                         "movaps %%xmm0, %0 \n\t"
                         "movaps %%xmm1, %1"
                         :
                         "=m" (rs1),
                         "=m" (rs2)
                         :
                         "m" ((*s).c2.c2.re),
                         "m" (_sse_sgn13)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

/*************************** components 3 and 4 *******************************/

   __asm__ __volatile__ ("movlps %0, %%xmm4 \n\t"
                         "movlps %1, %%xmm5 \n\t"
                         "movhps %2, %%xmm4 \n\t"
                         "movhps %3, %%xmm5 \n\t"
                         "movaps %%xmm4, %%xmm6 \n\t"
                         "movaps %%xmm5, %%xmm7 \n\t"
                         "movss %4, %%xmm0 \n\t"
                         "movss %5, %%xmm1"
                         :
                         :
                         "m" ((*m).u[26]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[32]),
                         "m" ((*m).u[2]),
                         "m" ((*m).u[24])
                         :
                         "xmm0", "xmm1", "xmm4", "xmm5",
                         "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("movss %0, %%xmm2 \n\t"
                         "movss %1, %%xmm3 \n\t"
                         "shufps $0xa0, %%xmm4, %%xmm4 \n\t"
                         "shufps $0xa0, %%xmm5, %%xmm5 \n\t"
                         "shufps $0xf5, %%xmm6, %%xmm6 \n\t"
                         "shufps $0xf5, %%xmm7, %%xmm7 \n\t"
                         "shufps $0x0, %%xmm1, %%xmm0 \n\t"
                         "shufps $0x0, %%xmm2, %%xmm1 \n\t"
                         "shufps $0x5, %%xmm3, %%xmm2 \n\t"
                         "shufps $0x50, %%xmm3, %%xmm3"
                         :
                         :
                         "m" ((*m).u[3]),
                         "m" ((*m).u[25])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm4 \n\t"
                         "mulps %0, %%xmm5 \n\t"
                         "mulps %0, %%xmm6 \n\t"
                         "mulps %0, %%xmm7 \n\t"
                         "mulps %1, %%xmm0 \n\t"
                         "mulps %1, %%xmm1 \n\t"
                         "mulps %1, %%xmm2 \n\t"
                         "mulps %1, %%xmm3 \n\t"
                         "addps %%xmm0, %%xmm4 \n\t"
                         "addps %%xmm1, %%xmm5 \n\t"
                         "addps %%xmm2, %%xmm6 \n\t"
                         "subps %%xmm3, %%xmm7"
                         :
                         :
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c1.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movlps %0, %%xmm0 \n\t"
                         "movlps %1, %%xmm1 \n\t"
                         "movhps %2, %%xmm0 \n\t"
                         "movhps %3, %%xmm1 \n\t"
                         "movaps %%xmm0, %%xmm2 \n\t"
                         "movaps %%xmm1, %%xmm3 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm1, %%xmm1 \n\t"
                         "shufps $0xf5, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xf5, %%xmm3, %%xmm3"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[10]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[18])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3");

   __asm__ __volatile__ ("mulps %2, %%xmm0 \n\t"
                         "mulps %2, %%xmm1 \n\t"
                         "mulps %2, %%xmm2 \n\t"
                         "mulps %2, %%xmm3 \n\t"
                         "addps %%xmm0, %%xmm4 \n\t"
                         "addps %%xmm1, %%xmm5 \n\t"
                         "subps %%xmm2, %%xmm6 \n\t"
                         "subps %%xmm3, %%xmm7 \n\t"
                         "shufps $0xb1, %%xmm6, %%xmm6 \n\t"
                         "shufps $0xb1, %%xmm7, %%xmm7 \n\t"
                         "mulps %3, %%xmm6 \n\t"
                         "mulps %3, %%xmm7 \n\t"                         
                         "addps %%xmm6, %%xmm4 \n\t"
                         "addps %%xmm7, %%xmm5 \n\t"
                         "movaps %%xmm4, %0 \n\t"
                         "movaps %%xmm5, %1"
                         :
                         "=m" (rs3),
                         "=m" (rs4)
                         :
                         "m" ((*s).c1.c1.re),
                         "m" (_sse_sgn13)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");
   
/*************************** components 5 and 6 *******************************/

   __asm__ __volatile__ ("movlps %0, %%xmm0 \n\t"
                         "movlps %1, %%xmm1 \n\t"
                         "movhps %2, %%xmm0 \n\t"
                         "movhps %3, %%xmm1 \n\t"                         
                         "movaps %%xmm0, %%xmm2 \n\t"
                         "movaps %%xmm1, %%xmm3 \n\t"
                         "movss %4, %%xmm4 \n\t"
                         "movss %5, %%xmm5"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[14]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[22]),
                         "m" ((*m).u[4]),
                         "m" ((*m).u[34])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5");

   __asm__ __volatile__ ("movss %0, %%xmm6 \n\t"
                         "movss %1, %%xmm7 \n\t"
                         "shufps $0xa0, %%xmm0, %%xmm0 \n\t"
                         "shufps $0xa0, %%xmm1, %%xmm1 \n\t"
                         "shufps $0xf5, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xf5, %%xmm3, %%xmm3 \n\t"
                         "shufps $0x0, %%xmm5, %%xmm4 \n\t"
                         "shufps $0x0, %%xmm6, %%xmm5 \n\t"
                         "shufps $0x5, %%xmm7, %%xmm6 \n\t"
                         "shufps $0x50, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[5]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm0 \n\t"
                         "mulps %0, %%xmm1 \n\t"
                         "mulps %0, %%xmm2 \n\t"
                         "mulps %0, %%xmm3 \n\t"
                         "mulps %1, %%xmm4 \n\t"
                         "mulps %1, %%xmm5 \n\t"
                         "mulps %1, %%xmm6 \n\t"
                         "mulps %1, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "subps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c2.c2.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movlps %0, %%xmm4 \n\t"
                         "movlps %1, %%xmm5 \n\t"
                         "movhps %2, %%xmm4 \n\t"
                         "movhps %3, %%xmm5 \n\t"                         
                         "movaps %%xmm4, %%xmm6 \n\t"
                         "movaps %%xmm5, %%xmm7 \n\t"
                         "shufps $0xa0, %%xmm4, %%xmm4 \n\t"
                         "shufps $0xa0, %%xmm5, %%xmm5 \n\t"
                         "shufps $0xf5, %%xmm6, %%xmm6 \n\t"
                         "shufps $0xf5, %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[26]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[30]),                         
                         "m" ((*m).u[32])
                         :
                         "xmm4", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulps %0, %%xmm4 \n\t"
                         "mulps %0, %%xmm5 \n\t"
                         "mulps %0, %%xmm6 \n\t"
                         "mulps %0, %%xmm7 \n\t"
                         "addps %%xmm4, %%xmm0 \n\t"
                         "addps %%xmm5, %%xmm1 \n\t"
                         "addps %%xmm6, %%xmm2 \n\t"
                         "addps %%xmm7, %%xmm3 \n\t"
                         "shufps $0xb1, %%xmm2, %%xmm2 \n\t"
                         "shufps $0xb1, %%xmm3, %%xmm3 \n\t"
                         "mulps %1, %%xmm2 \n\t"
                         "mulps %1, %%xmm3 \n\t"                         
                         "addps %%xmm2, %%xmm0 \n\t"
                         "addps %%xmm3, %%xmm1"
                         :
                         :
                         "m" ((*s).c1.c3.re),
                         "m" (_sse_sgn24)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6", "xmm7");

/*************************** combine and store ********************************/
   
   __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                         "movaps %1, %%xmm5 \n\t"
                         "movlhps %%xmm1, %%xmm2 \n\t"
                         "movlhps %%xmm3, %%xmm4 \n\t"
                         "movlhps %%xmm5, %%xmm6 \n\t"
                         "movhlps %%xmm0, %%xmm1 \n\t"
                         "shufps $0xe4, %%xmm2, %%xmm0 \n\t"
                         "movlps %2, %%xmm3 \n\t"
                         "movlps %3, %%xmm4 \n\t"
                         "movlps %4, %%xmm5 \n\t"
                         "movlps %5, %%xmm6 \n\t"
                         "addps %%xmm3, %%xmm4 \n\t"
                         "addps %%xmm5, %%xmm6 \n\t"
                         "addps %%xmm1, %%xmm0"
                         :
                         :
                         "m" (rs2.c1),
                         "m" (rs4.c1),
                         "m" (rs1.c1),
                         "m" (rs1.c3),
                         "m" (rs3.c1),                         
                         "m" (rs3.c3)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3",
                         "xmm4", "xmm5", "xmm6");

   __asm__ __volatile__ ("movaps %%xmm4, %0 \n\t"
                         "movaps %%xmm6, %1 \n\t"
                         "movaps %%xmm0, %2"
                         :
                         "=m" ((*r).c1.c1),
                         "=m" ((*r).c1.c3),
                         "=m" ((*r).c2.c2));
}

#else

static weyl rs;


void mul_pauli(pauli *m,weyl *s,weyl *r)
{
   rs.c1.c1.re=
      (*m).u[ 0]*(*s).c1.c1.re+
      (*m).u[ 6]*(*s).c1.c2.re-(*m).u[ 7]*(*s).c1.c2.im+
      (*m).u[ 8]*(*s).c1.c3.re-(*m).u[ 9]*(*s).c1.c3.im+
      (*m).u[10]*(*s).c2.c1.re-(*m).u[11]*(*s).c2.c1.im+
      (*m).u[12]*(*s).c2.c2.re-(*m).u[13]*(*s).c2.c2.im+
      (*m).u[14]*(*s).c2.c3.re-(*m).u[15]*(*s).c2.c3.im;

   rs.c1.c1.im=
      (*m).u[ 0]*(*s).c1.c1.im+
      (*m).u[ 6]*(*s).c1.c2.im+(*m).u[ 7]*(*s).c1.c2.re+
      (*m).u[ 8]*(*s).c1.c3.im+(*m).u[ 9]*(*s).c1.c3.re+
      (*m).u[10]*(*s).c2.c1.im+(*m).u[11]*(*s).c2.c1.re+
      (*m).u[12]*(*s).c2.c2.im+(*m).u[13]*(*s).c2.c2.re+
      (*m).u[14]*(*s).c2.c3.im+(*m).u[15]*(*s).c2.c3.re;

   rs.c1.c2.re=
      (*m).u[ 6]*(*s).c1.c1.re+(*m).u[ 7]*(*s).c1.c1.im+
      (*m).u[ 1]*(*s).c1.c2.re+
      (*m).u[16]*(*s).c1.c3.re-(*m).u[17]*(*s).c1.c3.im+
      (*m).u[18]*(*s).c2.c1.re-(*m).u[19]*(*s).c2.c1.im+
      (*m).u[20]*(*s).c2.c2.re-(*m).u[21]*(*s).c2.c2.im+
      (*m).u[22]*(*s).c2.c3.re-(*m).u[23]*(*s).c2.c3.im;

   rs.c1.c2.im=
      (*m).u[ 6]*(*s).c1.c1.im-(*m).u[ 7]*(*s).c1.c1.re+
      (*m).u[ 1]*(*s).c1.c2.im+
      (*m).u[16]*(*s).c1.c3.im+(*m).u[17]*(*s).c1.c3.re+
      (*m).u[18]*(*s).c2.c1.im+(*m).u[19]*(*s).c2.c1.re+
      (*m).u[20]*(*s).c2.c2.im+(*m).u[21]*(*s).c2.c2.re+
      (*m).u[22]*(*s).c2.c3.im+(*m).u[23]*(*s).c2.c3.re;

   rs.c1.c3.re=
      (*m).u[ 8]*(*s).c1.c1.re+(*m).u[ 9]*(*s).c1.c1.im+
      (*m).u[16]*(*s).c1.c2.re+(*m).u[17]*(*s).c1.c2.im+
      (*m).u[ 2]*(*s).c1.c3.re+
      (*m).u[24]*(*s).c2.c1.re-(*m).u[25]*(*s).c2.c1.im+
      (*m).u[26]*(*s).c2.c2.re-(*m).u[27]*(*s).c2.c2.im+
      (*m).u[28]*(*s).c2.c3.re-(*m).u[29]*(*s).c2.c3.im;

   rs.c1.c3.im=
      (*m).u[ 8]*(*s).c1.c1.im-(*m).u[ 9]*(*s).c1.c1.re+
      (*m).u[16]*(*s).c1.c2.im-(*m).u[17]*(*s).c1.c2.re+
      (*m).u[ 2]*(*s).c1.c3.im+
      (*m).u[24]*(*s).c2.c1.im+(*m).u[25]*(*s).c2.c1.re+
      (*m).u[26]*(*s).c2.c2.im+(*m).u[27]*(*s).c2.c2.re+
      (*m).u[28]*(*s).c2.c3.im+(*m).u[29]*(*s).c2.c3.re;

   rs.c2.c1.re=
      (*m).u[10]*(*s).c1.c1.re+(*m).u[11]*(*s).c1.c1.im+
      (*m).u[18]*(*s).c1.c2.re+(*m).u[19]*(*s).c1.c2.im+
      (*m).u[24]*(*s).c1.c3.re+(*m).u[25]*(*s).c1.c3.im+
      (*m).u[ 3]*(*s).c2.c1.re+
      (*m).u[30]*(*s).c2.c2.re-(*m).u[31]*(*s).c2.c2.im+
      (*m).u[32]*(*s).c2.c3.re-(*m).u[33]*(*s).c2.c3.im;

   rs.c2.c1.im=
      (*m).u[10]*(*s).c1.c1.im-(*m).u[11]*(*s).c1.c1.re+
      (*m).u[18]*(*s).c1.c2.im-(*m).u[19]*(*s).c1.c2.re+
      (*m).u[24]*(*s).c1.c3.im-(*m).u[25]*(*s).c1.c3.re+
      (*m).u[ 3]*(*s).c2.c1.im+
      (*m).u[30]*(*s).c2.c2.im+(*m).u[31]*(*s).c2.c2.re+
      (*m).u[32]*(*s).c2.c3.im+(*m).u[33]*(*s).c2.c3.re;

   rs.c2.c2.re=
      (*m).u[12]*(*s).c1.c1.re+(*m).u[13]*(*s).c1.c1.im+
      (*m).u[20]*(*s).c1.c2.re+(*m).u[21]*(*s).c1.c2.im+
      (*m).u[26]*(*s).c1.c3.re+(*m).u[27]*(*s).c1.c3.im+
      (*m).u[30]*(*s).c2.c1.re+(*m).u[31]*(*s).c2.c1.im+
      (*m).u[ 4]*(*s).c2.c2.re+
      (*m).u[34]*(*s).c2.c3.re-(*m).u[35]*(*s).c2.c3.im;

   rs.c2.c2.im=
      (*m).u[12]*(*s).c1.c1.im-(*m).u[13]*(*s).c1.c1.re+
      (*m).u[20]*(*s).c1.c2.im-(*m).u[21]*(*s).c1.c2.re+
      (*m).u[26]*(*s).c1.c3.im-(*m).u[27]*(*s).c1.c3.re+
      (*m).u[30]*(*s).c2.c1.im-(*m).u[31]*(*s).c2.c1.re+
      (*m).u[ 4]*(*s).c2.c2.im+
      (*m).u[34]*(*s).c2.c3.im+(*m).u[35]*(*s).c2.c3.re;

   rs.c2.c3.re=
      (*m).u[14]*(*s).c1.c1.re+(*m).u[15]*(*s).c1.c1.im+
      (*m).u[22]*(*s).c1.c2.re+(*m).u[23]*(*s).c1.c2.im+
      (*m).u[28]*(*s).c1.c3.re+(*m).u[29]*(*s).c1.c3.im+
      (*m).u[32]*(*s).c2.c1.re+(*m).u[33]*(*s).c2.c1.im+
      (*m).u[34]*(*s).c2.c2.re+(*m).u[35]*(*s).c2.c2.im+
      (*m).u[ 5]*(*s).c2.c3.re;

   rs.c2.c3.im=
      (*m).u[14]*(*s).c1.c1.im-(*m).u[15]*(*s).c1.c1.re+
      (*m).u[22]*(*s).c1.c2.im-(*m).u[23]*(*s).c1.c2.re+
      (*m).u[28]*(*s).c1.c3.im-(*m).u[29]*(*s).c1.c3.re+
      (*m).u[32]*(*s).c2.c1.im-(*m).u[33]*(*s).c2.c1.re+
      (*m).u[34]*(*s).c2.c2.im-(*m).u[35]*(*s).c2.c2.re+
      (*m).u[ 5]*(*s).c2.c3.im;

   *r=rs;
}

#endif

void assign_pauli(int vol,pauli_dble *md,pauli *m)
{
   float *r;
   double *s,*sm;

   r=(float*)(m);
   s=(double*)(md);
   sm=s+36*vol;

   for (;s<sm;s+=9)
   {
      *(r  )=(float)(*(s  ));
      *(r+1)=(float)(*(s+1));
      *(r+2)=(float)(*(s+2));
      *(r+3)=(float)(*(s+3));
      *(r+4)=(float)(*(s+4));
      *(r+5)=(float)(*(s+5));
      *(r+6)=(float)(*(s+6));
      *(r+7)=(float)(*(s+7));
      *(r+8)=(float)(*(s+8));
      r+=9;
   }
}


void apply_sw(int vol,pauli *m,spinor *s,spinor *r)
{
   weyl *rs,*rr,*rsm;

   rs=(weyl*)(s);
   rr=(weyl*)(r);
   rsm=rs+2*vol;

   for (;rs<rsm;)
   {
#if (defined SSE)
      m+=4;
      _prefetch_pauli_dble(m);
      m-=4;
#endif
      
      mul_pauli(m,rs,rr);

#if (defined SSE)
      m+=1;
      rs+=4;
      _prefetch_spinor(rs);
      rs-=3;
      rr+=4;
      _prefetch_spinor(rr);
      rr-=3;
#else
      m+=1;
      rs+=1;
      rr+=1;
#endif

      mul_pauli(m,rs,rr);

      m+=1;
      rs+=1;
      rr+=1;
   }
}
