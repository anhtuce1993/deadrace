
/*******************************************************************************
*
* File pauli_dble.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic functions for double-precision hermitian 6x6 matrices
*
* The externally accessible functions are
*
*   void mul_pauli_dble(pauli_dble *m,weyl_dble *s,weyl_dble *r)
*     Multiplies the Weyl spinor *s by the matrix *m and assigns the result
*     to the Weyl spinor *r. The source spinor is overwritten if r=s and
*     otherwise left unchanged
*
*   int inv_pauli_dble(pauli_dble *m,pauli_dble *im)
*     Inverts the matrix *m and assigns the result to *im. The matrix is
*     overwritten if im=m and otherwise left unchanged. On exit the
*     program returns 0 or 1 depending on whether the inversion was safe
*     or not (in which case the calculated matrix is unusable)
*
*   double det_pauli_dble(pauli_dble *m)
*     Returns the determinant of the matrix *m
*
*   void apply_sw_dble(int vol,pauli_dble *m,spinor_dble *s,spinor_dble *r)
*     Applies the pauli matrix field m[2*vol] to the spinor field s[vol]
*     and assigns the result to the field r[vol]. The source field is
*     overwritten if r=s and otherwise left unchanged (the arrays may
*     not overlap in this case)
*
* Notes:
*
* The storage format for hermitian 6x6 matrices is described in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.ps). As
* explained there, the matrix inversion is considered to be safe if and
* only if the Frobenius condition number of the matrix is less than 100
*
* If SSE/SSE2/SSE3 instructions are used, it is assumed that the matrices
* and spinors are aligned to a 16 byte boundary
*
*******************************************************************************/

#define PAULI_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"

#define DELTA 1.0e-04

#if (defined SSE2)
#include "sse2.h"

static double rr[5] __attribute__ ((aligned (8)));
static complex_dble aa[36],dd[6]  __attribute__ ((aligned (16)));
static sse_double rs1,rs2,rs3,rs4,rs5;

#endif

#if (defined SSE3)

void mul_pauli_dble(pauli_dble *m,weyl_dble *s,weyl_dble *r)
{
   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[7]),
                         "m" ((*m).u[9]),
                         "m" ((*m).u[11]),
                         "m" ((*m).u[13]),
                         "m" ((*m).u[15])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");
                         
   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[0]),
                         "m" ((*m).u[6]),
                         "m" ((*m).u[8]),
                         "m" ((*m).u[10]),
                         "m" ((*m).u[12]),
                         "m" ((*m).u[14])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");
   
   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("addpd %%xmm5, %%xmm7 \n\t"
                         "addsubpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs1.c1)
                         :
                         :
                         "xmm7");
                         

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[7]),
                         "m" ((*m).u[17]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[23])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "subpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[6]),
                         "m" ((*m).u[1]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[18]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[22])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("addpd %%xmm5, %%xmm7 \n\t"
                         "addsubpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs2.c1)
                         :
                         :
                         "xmm7");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[9]),
                         "m" ((*m).u[17]),
                         "m" ((*m).u[25]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[29])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[2]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[26]),
                         "m" ((*m).u[28])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("addpd %%xmm5, %%xmm7 \n\t"
                         "addsubpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs3.c1)
                         :
                         :
                         "xmm7");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[11]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[25]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[33])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "subpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[10]),
                         "m" ((*m).u[18]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[3]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[32])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("addpd %%xmm5, %%xmm7 \n\t"
                         "addsubpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs4.c1)
                         :
                         :
                         "xmm7");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[13]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[26]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[4]),
                         "m" ((*m).u[34])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "subpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("addpd %%xmm5, %%xmm7 \n\t"
                         "addsubpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs5.c1)
                         :
                         :
                         "xmm7");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm3 \n\t"
                         "movddup %4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[15]),
                         "m" ((*m).u[23]),
                         "m" ((*m).u[29]),
                         "m" ((*m).u[33]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                         "movddup %1, %%xmm1 \n\t"
                         "movddup %2, %%xmm2 \n\t"
                         "movddup %3, %%xmm5 \n\t"
                         "movddup %4, %%xmm6 \n\t"
                         "movddup %5, %%xmm7"
                         :
                         :
                         "m" ((*m).u[14]),
                         "m" ((*m).u[22]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[32]),
                         "m" ((*m).u[34]),
                         "m" ((*m).u[5])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "movapd %3, %%xmm3 \n\t"
                         "movapd %4, %%xmm6"
                         :
                         :
                         "m" (rs1.c1),
                         "m" (rs2.c1),
                         "m" (rs3.c1),
                         "m" (rs4.c1),
                         "m" (rs5.c1)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm6");

   __asm__ __volatile__ ("xorpd %6, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "subpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm0, %0 \n\t"
                         "movapd %%xmm1, %1 \n\t"
                         "movapd %%xmm2, %2 \n\t"
                         "movapd %%xmm3, %3 \n\t"
                         "movapd %%xmm6, %4 \n\t"
                         "movapd %%xmm7, %5"
                         :
                         "=m" ((*r).c1.c1.re),
                         "=m" ((*r).c1.c2.re),
                         "=m" ((*r).c1.c3.re),
                         "=m" ((*r).c2.c1.re),
                         "=m" ((*r).c2.c2.re),
                         "=m" ((*r).c2.c3.re)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");
}


static int fwd_house(double eps)
{
   int i,j,k,itest;
   double r1,r2,r3;
   complex_dble z,*ak,*aj;

   itest=0;

   for (k=0;k<5;k++)
   {
      r1=aa[6*k+k].re*aa[6*k+k].re+aa[6*k+k].im*aa[6*k+k].im;
      r2=sqrt(r1);

      for (j=(k+1);j<6;j++)
         r1+=(aa[6*j+k].re*aa[6*j+k].re+aa[6*j+k].im*aa[6*j+k].im);

      if (r1>=eps)
         r1=sqrt(r1);
      else
      {
         itest=1;
         r1=1.0;
      }

      if (r2>=(DBL_EPSILON*r1))
      {
         r3=1.0/r2;
         z.re=r3*aa[6*k+k].re;
         z.im=r3*aa[6*k+k].im;
      }
      else
      {
         z.re=1.0;
         z.im=0.0;
      }

      aa[6*k+k].re+=r1*z.re;
      aa[6*k+k].im+=r1*z.im;

      r3=1.0/(r1*(r1+r2));
      rr[k]=r3;
      dd[k].re=-(r1+r2)*r3*z.re;
      dd[k].im= (r1+r2)*r3*z.im;

      for (j=(k+1);j<6;j++)
      {
         __asm__ __volatile__ ("xorpd %%xmm7, %%xmm7"
                               :
                               :
                               :
                               "xmm7");

         ak=aa+6*k+k;
         aj=aa+6*k+j;

         for (i=k;i<6;i++)
         {
            __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                                  "movddup %1, %%xmm1 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "addpd %%xmm0, %%xmm7 \n\t"
                                  "xorpd %3, %%xmm1 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "addpd %%xmm1, %%xmm7"
                                  :
                                  :
                                  "m" ((*ak).re),
                                  "m" ((*ak).im),
                                  "m" ((*aj).re),
                                  "m" (_sse_sgn_dble)
                                  :
                                  "xmm0", "xmm1", "xmm7");

            ak+=6;
            aj+=6;
         }

         __asm__ __volatile__ ("movddup %0, %%xmm5 \n\t"
                               "mulpd %%xmm5, %%xmm7 \n\t"
                               "movddup %%xmm7, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "xorpd %1, %%xmm7"
                               :
                               :
                               "m" (rr[k]),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm5", "xmm6", "xmm7");

         ak=aa+6*k+k;
         aj=aa+6*k+j;

         for (i=k;i<6;i++)
         {
            __asm__ __volatile__ ("movapd %%xmm7, %%xmm5 \n\t"
                                  "movapd %%xmm6, %%xmm4 \n\t"
                                  "mulpd %2, %%xmm5 \n\t"
                                  "mulpd %2, %%xmm4 \n\t"
                                  "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                                  "subpd %1, %%xmm4 \n\t"
                                  "subpd %%xmm4, %%xmm5 \n\t"
                                  "movapd %%xmm5, %0"
                                  :
                                  "=m" ((*aj).re)
                                  :
                                  "m" ((*aj).re),
                                  "m" ((*ak).re)
                                  :
                                  "xmm4", "xmm5");

            ak+=6;
            aj+=6;
         }
      }
   }

   r1=aa[35].re*aa[35].re+aa[35].im*aa[35].im;

   if (r1>=eps)
      r1=1.0/r1;
   else
   {
      itest=1;
      r1=1.0;
   }

   dd[5].re= r1*aa[35].re;
   dd[5].im=-r1*aa[35].im;

   return itest;
}


static void solv_sys(void)
{
   int i,j,k;

   for (k=5;k>0;k--)
   {
      for (i=(k-1);i>=0;i--)
      {
         __asm__ __volatile__ ("movddup %0, %%xmm6 \n\t"
                               "movddup %1, %%xmm7 \n\t"
                               "mulpd %2, %%xmm6 \n\t"
                               "mulpd %2, %%xmm7 \n\t"
                               "shufpd $0x1, %%xmm6, %%xmm6 \n\t"
                               "addsubpd %%xmm6, %%xmm7"
                               :
                               :
                               "m" (aa[6*i+k].im),
                               "m" (aa[6*i+k].re),
                               "m" (dd[k].re)
                               :
                               "xmm6", "xmm7");

         for (j=(k-1);j>i;j--)
         {
            __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                                  "movddup %1, %%xmm1 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "addpd %%xmm0, %%xmm7 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "addsubpd %%xmm1, %%xmm7"
                                  :
                                  :
                                  "m" (aa[6*j+k].re),
                                  "m" (aa[6*j+k].im),
                                  "m" (aa[6*i+j].re)
                                  :
                                  "xmm0", "xmm1", "xmm7");
         }

         __asm__ __volatile__ ("movddup %%xmm7, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "mulpd %1, %%xmm7 \n\t"
                               "mulpd %1, %%xmm6 \n\t"
                               "xorpd %2, %%xmm7 \n\t"
                               "shufpd $0x1, %%xmm7, %%xmm7 \n\t"
                               "subpd %%xmm6, %%xmm7 \n\t"
                               "movapd %%xmm7, %0"
                               :
                               "=m" (aa[6*i+k].re)
                               :
                               "m" (dd[i].re),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm6", "xmm7");
      }
   }
}


static void bck_house(void)
{
   int i,j,k;
   complex_dble z,*d,*a;

   aa[35].re=dd[5].re;
   aa[35].im=dd[5].im;

   for (k=4;k>=0;k--)
   {
      z.re=dd[k].re;
      z.im=dd[k].im;
      dd[k].re=aa[6*k+k].re;
      dd[k].im=aa[6*k+k].im;
      aa[6*k+k].re=z.re;
      aa[6*k+k].im=z.im;

      for (j=(k+1);j<6;j++)
      {
         dd[j].re=aa[6*j+k].re;
         dd[j].im=aa[6*j+k].im;
         aa[6*j+k].re=0.0;
         aa[6*j+k].im=0.0;
      }

      for (i=0;i<6;i+=2)
      {
         __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                               "xorpd %%xmm7, %%xmm7"
                               :
                               :
                               :
                               "xmm6", "xmm7");

         d=dd+k;
         a=aa+6*i+k;

         for (j=k;j<6;j++)
         {
            __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                                  "movddup %1, %%xmm1 \n\t"
                                  "movapd %%xmm0, %%xmm2 \n\t"
                                  "movapd %%xmm1, %%xmm3 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "mulpd %3, %%xmm2 \n\t"
                                  "mulpd %3, %%xmm3 \n\t"
                                  "addpd %%xmm0, %%xmm6 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "addpd %%xmm2, %%xmm7 \n\t"
                                  "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                                  "addsubpd %%xmm1, %%xmm6 \n\t"
                                  "addsubpd %%xmm3, %%xmm7 \n\t"
                                  :
                                  :
                                  "m" ((*d).re),
                                  "m" ((*d).im),
                                  "m" ((*a).re),
                                  "m" ((*(a+6)).re)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm6", "xmm7");

            d+=1;
            a+=1;
         }

         __asm__ __volatile__ ("movddup %0, %%xmm0 \n\t"
                               "mulpd %%xmm0, %%xmm6 \n\t"
                               "mulpd %%xmm0, %%xmm7 \n\t"
                               "movddup %%xmm6, %%xmm4 \n\t"
                               "movddup %%xmm7, %%xmm5 \n\t"
                               "unpckhpd %%xmm6, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "xorpd %1, %%xmm4 \n\t"
                               "xorpd %1, %%xmm5"
                               :
                               :
                               "m" (rr[k]),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm0", "xmm4", "xmm5",
                               "xmm6", "xmm7");

         d=dd+k;
         a=aa+6*i+k;

         for (j=k;j<6;j++)
         {
            __asm__ __volatile__ ("movapd %%xmm6, %%xmm2 \n\t"
                                  "movapd %%xmm7, %%xmm3 \n\t"
                                  "movapd %%xmm4, %%xmm0 \n\t"
                                  "movapd %%xmm5, %%xmm1 \n\t"
                                  "mulpd %4, %%xmm2 \n\t"
                                  "mulpd %4, %%xmm3 \n\t"
                                  "mulpd %4, %%xmm0 \n\t"
                                  "mulpd %4, %%xmm1 \n\t"
                                  "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                                  "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                                  "addpd %2, %%xmm0 \n\t"
                                  "addpd %3, %%xmm1 \n\t"
                                  "subpd %%xmm2, %%xmm0 \n\t"
                                  "subpd %%xmm3, %%xmm1 \n\t"
                                  "movapd %%xmm0, %0 \n\t"
                                  "movapd %%xmm1, %1"
                                  :
                                  "=m" ((*a).re),
                                  "=m" ((*(a+6)).re)
                                  :
                                  "m" ((*a).re),
                                  "m" ((*(a+6)).re),
                                  "m" ((*d).re)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3");

            d+=1;
            a+=1;
         }
      }
   }
}

#elif (defined SSE2)

void mul_pauli_dble(pauli_dble *m,weyl_dble *s,weyl_dble *r)
{
   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[7]),
                         "m" ((*m).u[9]),
                         "m" ((*m).u[11]),
                         "m" ((*m).u[13]),
                         "m" ((*m).u[15])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[0]),
                         "m" ((*m).u[6]),
                         "m" ((*m).u[8]),
                         "m" ((*m).u[10]),
                         "m" ((*m).u[12]),
                         "m" ((*m).u[14])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("xorpd %1, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "addpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs1.c1)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[7]),
                         "m" ((*m).u[17]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[23])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "subpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[6]),
                         "m" ((*m).u[1]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[18]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[22])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("xorpd %1, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "addpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs2.c1)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[9]),
                         "m" ((*m).u[17]),
                         "m" ((*m).u[25]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[29])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");
                         
   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[8]),
                         "m" ((*m).u[16]),
                         "m" ((*m).u[2]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[26]),
                         "m" ((*m).u[28])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("xorpd %1, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "addpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs3.c1)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[11]),
                         "m" ((*m).u[19]),
                         "m" ((*m).u[25]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[33])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "subpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");
                         
   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[10]),
                         "m" ((*m).u[18]),
                         "m" ((*m).u[24]),
                         "m" ((*m).u[3]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[32])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("xorpd %1, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "addpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs4.c1)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[13]),
                         "m" ((*m).u[21]),
                         "m" ((*m).u[27]),
                         "m" ((*m).u[31]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "subpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[12]),
                         "m" ((*m).u[20]),
                         "m" ((*m).u[26]),
                         "m" ((*m).u[30]),
                         "m" ((*m).u[4]),
                         "m" ((*m).u[34])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "subpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");
                         
   __asm__ __volatile__ ("xorpd %1, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "addpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm7, %0"
                         :
                         "=m" (rs5.c1)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");
                         
   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm3 \n\t"
                         "movsd %4, %%xmm4 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm3, %%xmm3 \n\t"
                         "unpcklpd %%xmm4, %%xmm4"
                         :
                         :
                         "m" ((*m).u[15]),
                         "m" ((*m).u[23]),
                         "m" ((*m).u[29]),
                         "m" ((*m).u[33]),
                         "m" ((*m).u[35])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm3 \n\t"
                         "mulpd %4, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm3 \n\t"
                         "addpd %%xmm1, %%xmm4"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm4");

   __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                         "movsd %1, %%xmm1 \n\t"
                         "movsd %2, %%xmm2 \n\t"
                         "movsd %3, %%xmm5 \n\t"
                         "movsd %4, %%xmm6 \n\t"
                         "movsd %5, %%xmm7 \n\t"
                         "unpcklpd %%xmm0, %%xmm0 \n\t"
                         "unpcklpd %%xmm1, %%xmm1 \n\t"
                         "unpcklpd %%xmm2, %%xmm2 \n\t"
                         "unpcklpd %%xmm5, %%xmm5 \n\t"
                         "unpcklpd %%xmm6, %%xmm6 \n\t"
                         "unpcklpd %%xmm7, %%xmm7"
                         :
                         :
                         "m" ((*m).u[14]),
                         "m" ((*m).u[22]),
                         "m" ((*m).u[28]),
                         "m" ((*m).u[32]),
                         "m" ((*m).u[34]),
                         "m" ((*m).u[5])
                         :
                         "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t"
                         "mulpd %1, %%xmm1 \n\t"
                         "mulpd %2, %%xmm2 \n\t"
                         "mulpd %3, %%xmm5 \n\t"
                         "mulpd %4, %%xmm6 \n\t"
                         "mulpd %5, %%xmm7 \n\t"
                         "addpd %%xmm3, %%xmm4 \n\t"
                         "addpd %%xmm0, %%xmm1 \n\t"
                         "addpd %%xmm2, %%xmm5 \n\t"
                         "shufpd $0x1, %%xmm4, %%xmm4 \n\t"
                         "addpd %%xmm6, %%xmm7 \n\t"
                         "addpd %%xmm1, %%xmm5"
                         :
                         :
                         "m" ((*s).c1.c1.re),
                         "m" ((*s).c1.c2.re),
                         "m" ((*s).c1.c3.re),
                         "m" ((*s).c2.c1.re),
                         "m" ((*s).c2.c2.re),
                         "m" ((*s).c2.c3.re)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm4",
                         "xmm5", "xmm6", "xmm7");

   __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                         "movapd %1, %%xmm1 \n\t"
                         "movapd %2, %%xmm2 \n\t"
                         "movapd %3, %%xmm3 \n\t"
                         "movapd %4, %%xmm6"
                         :
                         :
                         "m" (rs1.c1),
                         "m" (rs2.c1),
                         "m" (rs3.c1),
                         "m" (rs4.c1),
                         "m" (rs5.c1)
                         :
                         "xmm0", "xmm1", "xmm2", "xmm3", "xmm6");

   __asm__ __volatile__ ("xorpd %6, %%xmm4 \n\t"
                         "addpd %%xmm5, %%xmm7 \n\t"
                         "subpd %%xmm4, %%xmm7 \n\t"
                         "movapd %%xmm0, %0 \n\t"
                         "movapd %%xmm1, %1 \n\t"
                         "movapd %%xmm2, %2 \n\t"
                         "movapd %%xmm3, %3 \n\t"
                         "movapd %%xmm6, %4 \n\t"
                         "movapd %%xmm7, %5"
                         :
                         "=m" ((*r).c1.c1.re),
                         "=m" ((*r).c1.c2.re),
                         "=m" ((*r).c1.c3.re),
                         "=m" ((*r).c2.c1.re),
                         "=m" ((*r).c2.c2.re),
                         "=m" ((*r).c2.c3.re)
                         :
                         "m" (_sse_sgn_dble)
                         :
                         "xmm4", "xmm7");
}


static int fwd_house(double eps)
{
   int i,j,k,itest;
   double r1,r2,r3;
   complex_dble z,*ak,*aj;

   itest=0;

   for (k=0;k<5;k++)
   {
      r1=aa[6*k+k].re*aa[6*k+k].re+aa[6*k+k].im*aa[6*k+k].im;
      r2=sqrt(r1);

      for (j=(k+1);j<6;j++)
         r1+=(aa[6*j+k].re*aa[6*j+k].re+aa[6*j+k].im*aa[6*j+k].im);

      if (r1>=eps)
         r1=sqrt(r1);
      else
      {
         itest=1;
         r1=1.0;
      }

      if (r2>=(DBL_EPSILON*r1))
      {
         r3=1.0/r2;
         z.re=r3*aa[6*k+k].re;
         z.im=r3*aa[6*k+k].im;
      }
      else
      {
         z.re=1.0;
         z.im=0.0;
      }

      aa[6*k+k].re+=r1*z.re;
      aa[6*k+k].im+=r1*z.im;

      r3=1.0/(r1*(r1+r2));
      rr[k]=r3;
      dd[k].re=-(r1+r2)*r3*z.re;
      dd[k].im= (r1+r2)*r3*z.im;

      for (j=(k+1);j<6;j++)
      {
         __asm__ __volatile__ ("xorpd %%xmm7, %%xmm7"
                               :
                               :
                               :
                               "xmm7");

         ak=aa+6*k+k;
         aj=aa+6*k+j;

         for (i=k;i<6;i++)
         {
            __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                                  "movsd %1, %%xmm1 \n\t"
                                  "unpcklpd %%xmm0, %%xmm0 \n\t"
                                  "unpcklpd %%xmm1, %%xmm1 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "addpd %%xmm0, %%xmm7 \n\t"
                                  "xorpd %3, %%xmm1 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "addpd %%xmm1, %%xmm7"
                                  :
                                  :
                                  "m" ((*ak).re),
                                  "m" ((*ak).im),
                                  "m" ((*aj).re),
                                  "m" (_sse_sgn_dble)
                                  :
                                  "xmm0", "xmm1", "xmm7");

            ak+=6;
            aj+=6;
         }

         __asm__ __volatile__ ("movsd %0, %%xmm5 \n\t"
                               "unpcklpd %%xmm5, %%xmm5 \n\t"
                               "mulpd %%xmm5, %%xmm7 \n\t"
                               "movsd %%xmm7, %%xmm6 \n\t"
                               "unpcklpd %%xmm6, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "xorpd %1, %%xmm7"
                               :
                               :
                               "m" (rr[k]),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm5", "xmm6", "xmm7");

         ak=aa+6*k+k;
         aj=aa+6*k+j;

         for (i=k;i<6;i++)
         {
            __asm__ __volatile__ ("movapd %%xmm7, %%xmm5 \n\t"
                                  "movapd %%xmm6, %%xmm4 \n\t"
                                  "mulpd %2, %%xmm5 \n\t"
                                  "mulpd %2, %%xmm4 \n\t"
                                  "shufpd $0x1, %%xmm5, %%xmm5 \n\t"
                                  "subpd %1, %%xmm4 \n\t"
                                  "subpd %%xmm4, %%xmm5 \n\t"
                                  "movapd %%xmm5, %0"
                                  :
                                  "=m" ((*aj).re)
                                  :
                                  "m" ((*aj).re),
                                  "m" ((*ak).re)
                                  :
                                  "xmm4", "xmm5", "xmm6", "xmm7");

            ak+=6;
            aj+=6;
         }
      }
   }

   r1=aa[35].re*aa[35].re+aa[35].im*aa[35].im;

   if (r1>=eps)
      r1=1.0/r1;
   else
   {
      itest=1;
      r1=1.0;
   }

   dd[5].re= r1*aa[35].re;
   dd[5].im=-r1*aa[35].im;

   return itest;
}


static void solv_sys(void)
{
   int i,j,k;

   for (k=5;k>0;k--)
   {
      for (i=(k-1);i>=0;i--)
      {
         __asm__ __volatile__ ("movsd %0, %%xmm6 \n\t"
                               "movsd %1, %%xmm7 \n\t"
                               "unpcklpd %%xmm6, %%xmm6 \n\t"
                               "unpcklpd %%xmm7, %%xmm7 \n\t"
                               "mulpd %2, %%xmm6 \n\t"
                               "mulpd %2, %%xmm7 \n\t"
                               "shufpd $0x1, %%xmm6, %%xmm6 \n\t"
                               "xorpd %3, %%xmm6 \n\t"
                               "addpd %%xmm6, %%xmm7"
                               :
                               :
                               "m" (aa[6*i+k].im),
                               "m" (aa[6*i+k].re),
                               "m" (dd[k].re),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm6", "xmm7");

         for (j=(k-1);j>i;j--)
         {
            __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                                  "movsd %1, %%xmm1 \n\t"
                                  "unpcklpd %%xmm0, %%xmm0 \n\t"
                                  "unpcklpd %%xmm1, %%xmm1 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "addpd %%xmm0, %%xmm7 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "xorpd %3, %%xmm1 \n\t"
                                  "addpd %%xmm1, %%xmm7"
                                  :
                                  :
                                  "m" (aa[6*j+k].re),
                                  "m" (aa[6*j+k].im),
                                  "m" (aa[6*i+j].re),
                                  "m" (_sse_sgn_dble)
                                  :
                                  "xmm0", "xmm1", "xmm7");
         }

         __asm__ __volatile__ ("movsd %%xmm7, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "unpcklpd %%xmm6, %%xmm6 \n\t"
                               "mulpd %1, %%xmm7 \n\t"
                               "mulpd %1, %%xmm6 \n\t"
                               "xorpd %2, %%xmm7 \n\t"
                               "shufpd $0x1, %%xmm7, %%xmm7 \n\t"
                               "subpd %%xmm6, %%xmm7 \n\t"
                               "movapd %%xmm7, %0"
                               :
                               "=m" (aa[6*i+k].re)
                               :
                               "m" (dd[i].re),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm6", "xmm7");
      }
   }
}


static void bck_house(void)
{
   int i,j,k;
   complex_dble z,*d,*a;

   aa[35].re=dd[5].re;
   aa[35].im=dd[5].im;

   for (k=4;k>=0;k--)
   {
      z.re=dd[k].re;
      z.im=dd[k].im;
      dd[k].re=aa[6*k+k].re;
      dd[k].im=aa[6*k+k].im;
      aa[6*k+k].re=z.re;
      aa[6*k+k].im=z.im;

      for (j=(k+1);j<6;j++)
      {
         dd[j].re=aa[6*j+k].re;
         dd[j].im=aa[6*j+k].im;
         aa[6*j+k].re=0.0;
         aa[6*j+k].im=0.0;
      }

      for (i=0;i<6;i+=2)
      {
         __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                               "xorpd %%xmm7, %%xmm7"
                               :
                               :
                               :
                               "xmm6", "xmm7");

         d=dd+k;
         a=aa+6*i+k;

         for (j=k;j<6;j++)
         {
            __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                                  "movsd %1, %%xmm1 \n\t"
                                  "unpcklpd %%xmm0, %%xmm0 \n\t"
                                  "unpcklpd %%xmm1, %%xmm1 \n\t"
                                  "movapd %%xmm0, %%xmm2 \n\t"
                                  "movapd %%xmm1, %%xmm3 \n\t"
                                  "mulpd %2, %%xmm0 \n\t"
                                  "mulpd %2, %%xmm1 \n\t"
                                  "mulpd %3, %%xmm2 \n\t"
                                  "mulpd %3, %%xmm3 \n\t"
                                  "addpd %%xmm0, %%xmm6 \n\t"
                                  "shufpd $0x1, %%xmm1, %%xmm1 \n\t"
                                  "addpd %%xmm2, %%xmm7 \n\t"
                                  "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                                  "xorpd %4, %%xmm1 \n\t"
                                  "xorpd %4, %%xmm3 \n\t"
                                  "addpd %%xmm1, %%xmm6 \n\t"
                                  "addpd %%xmm3, %%xmm7 \n\t"
                                  :
                                  :
                                  "m" ((*d).re),
                                  "m" ((*d).im),
                                  "m" ((*a).re),
                                  "m" ((*(a+6)).re),
                                  "m" (_sse_sgn_dble)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3",
                                  "xmm6", "xmm7");
                                  
            d+=1;
            a+=1;
         }

         __asm__ __volatile__ ("movsd %0, %%xmm0 \n\t"
                               "unpcklpd %%xmm0, %%xmm0 \n\t"
                               "mulpd %%xmm0, %%xmm6 \n\t"
                               "mulpd %%xmm0, %%xmm7 \n\t"
                               "movsd %%xmm6, %%xmm4 \n\t"
                               "movsd %%xmm7, %%xmm5 \n\t"
                               "unpckhpd %%xmm6, %%xmm6 \n\t"
                               "unpckhpd %%xmm7, %%xmm7 \n\t"
                               "unpcklpd %%xmm4, %%xmm4 \n\t"
                               "unpcklpd %%xmm5, %%xmm5 \n\t"
                               "xorpd %1, %%xmm4 \n\t"
                               "xorpd %1, %%xmm5"
                               :
                               :
                               "m" (rr[k]),
                               "m" (_sse_sgn_dble)
                               :
                               "xmm0", "xmm4", "xmm5", "xmm6", "xmm7");

         d=dd+k;
         a=aa+6*i+k;

         for (j=k;j<6;j++)
         {
            __asm__ __volatile__ ("movapd %%xmm6, %%xmm2 \n\t"
                                  "movapd %%xmm7, %%xmm3 \n\t"
                                  "movapd %%xmm4, %%xmm0 \n\t"
                                  "movapd %%xmm5, %%xmm1 \n\t"
                                  "mulpd %4, %%xmm2 \n\t"
                                  "mulpd %4, %%xmm3 \n\t"
                                  "mulpd %4, %%xmm0 \n\t"
                                  "mulpd %4, %%xmm1 \n\t"
                                  "shufpd $0x1, %%xmm2, %%xmm2 \n\t"
                                  "shufpd $0x1, %%xmm3, %%xmm3 \n\t"
                                  "addpd %2, %%xmm0 \n\t"
                                  "addpd %3, %%xmm1 \n\t"
                                  "subpd %%xmm2, %%xmm0 \n\t"
                                  "subpd %%xmm3, %%xmm1 \n\t"
                                  "movapd %%xmm0, %0 \n\t"
                                  "movapd %%xmm1, %1"
                                  :
                                  "=m" ((*a).re),
                                  "=m" ((*(a+6)).re)
                                  :
                                  "m" ((*a).re),
                                  "m" ((*(a+6)).re),
                                  "m" ((*d).re)
                                  :
                                  "xmm0", "xmm1", "xmm2", "xmm3");

            d+=1;
            a+=1;
         }
      }
   }
}

#else

static double rr[5];
static complex_dble aa[36],dd[6];
static weyl_dble rs;


void mul_pauli_dble(pauli_dble *m,weyl_dble *s,weyl_dble *r)
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


static int fwd_house(double eps)
{
   int i,j,k,itest;
   double r1,r2,r3;
   complex_dble z;

   itest=0;

   for (k=0;k<5;k++)
   {
      r1=aa[6*k+k].re*aa[6*k+k].re+aa[6*k+k].im*aa[6*k+k].im;
      r2=sqrt(r1);

      for (j=(k+1);j<6;j++)
         r1+=(aa[6*j+k].re*aa[6*j+k].re+aa[6*j+k].im*aa[6*j+k].im);

      if (r1>=eps)
         r1=sqrt(r1);
      else
      {
         itest=1;
         r1=1.0;
      }

      if (r2>=(DBL_EPSILON*r1))
      {
         r3=1.0/r2;
         z.re=r3*aa[6*k+k].re;
         z.im=r3*aa[6*k+k].im;
      }
      else
      {
         z.re=1.0;
         z.im=0.0;
      }

      aa[6*k+k].re+=r1*z.re;
      aa[6*k+k].im+=r1*z.im;

      r3=1.0/(r1*(r1+r2));
      rr[k]=r3;
      dd[k].re=-(r1+r2)*r3*z.re;
      dd[k].im= (r1+r2)*r3*z.im;

      for (j=(k+1);j<6;j++)
      {
         z.re=0.0;
         z.im=0.0;

         for (i=k;i<6;i++)
         {
            z.re+=(aa[6*i+k].re*aa[6*i+j].re+aa[6*i+k].im*aa[6*i+j].im);
            z.im+=(aa[6*i+k].re*aa[6*i+j].im-aa[6*i+k].im*aa[6*i+j].re);
         }

         z.re*=r3;
         z.im*=r3;

         for (i=k;i<6;i++)
         {
            aa[6*i+j].re-=(z.re*aa[6*i+k].re-z.im*aa[6*i+k].im);
            aa[6*i+j].im-=(z.re*aa[6*i+k].im+z.im*aa[6*i+k].re);
         }
      }
   }

   r1=aa[35].re*aa[35].re+aa[35].im*aa[35].im;

   if (r1>=eps)
      r1=1.0/r1;
   else
   {
      itest=1;
      r1=1.0;
   }

   dd[5].re= r1*aa[35].re;
   dd[5].im=-r1*aa[35].im;

   return itest;
}


static void solv_sys(void)
{
   int i,j,k;
   complex_dble z;

   for (k=5;k>0;k--)
   {
      for (i=(k-1);i>=0;i--)
      {
         z.re=aa[6*i+k].re*dd[k].re-aa[6*i+k].im*dd[k].im;
         z.im=aa[6*i+k].re*dd[k].im+aa[6*i+k].im*dd[k].re;

         for (j=(k-1);j>i;j--)
         {
            z.re+=(aa[6*i+j].re*aa[6*j+k].re-aa[6*i+j].im*aa[6*j+k].im);
            z.im+=(aa[6*i+j].re*aa[6*j+k].im+aa[6*i+j].im*aa[6*j+k].re);
         }

         aa[6*i+k].re=-dd[i].re*z.re+dd[i].im*z.im;
         aa[6*i+k].im=-dd[i].re*z.im-dd[i].im*z.re;
      }
   }
}


static void bck_house(void)
{
   int i,j,k;
   complex_dble z;

   aa[35].re=dd[5].re;
   aa[35].im=dd[5].im;

   for (k=4;k>=0;k--)
   {
      z.re=dd[k].re;
      z.im=dd[k].im;
      dd[k].re=aa[6*k+k].re;
      dd[k].im=aa[6*k+k].im;
      aa[6*k+k].re=z.re;
      aa[6*k+k].im=z.im;

      for (j=(k+1);j<6;j++)
      {
         dd[j].re=aa[6*j+k].re;
         dd[j].im=aa[6*j+k].im;
         aa[6*j+k].re=0.0;
         aa[6*j+k].im=0.0;
      }

      for (i=0;i<6;i++)
      {
         z.re=0.0;
         z.im=0.0;

         for (j=k;j<6;j++)
         {
            z.re+=(aa[6*i+j].re*dd[j].re-aa[6*i+j].im*dd[j].im);
            z.im+=(aa[6*i+j].re*dd[j].im+aa[6*i+j].im*dd[j].re);
         }

         z.re*=rr[k];
         z.im*=rr[k];

         for (j=k;j<6;j++)
         {
            aa[6*i+j].re-=(z.re*dd[j].re+z.im*dd[j].im);
            aa[6*i+j].im+=(z.re*dd[j].im-z.im*dd[j].re);
         }
      }
   }
}

#endif

int inv_pauli_dble(pauli_dble *m,pauli_dble *im)
{
   int i,j,itest;
   double eps,sm,*u;
   complex_dble *z;

   sm=0.0;
   u=(*m).u;
   z=(complex_dble*)(u+6);

   for (i=0;i<6;i++)
   {
      aa[6*i+i].re=*u;
      aa[6*i+i].im=0.0;
      sm+=(*u)*(*u);
      u+=1;

      for (j=i+1;j<6;j++)
      {
         aa[6*i+j].re= (*z).re;
         aa[6*i+j].im= (*z).im;
         aa[6*j+i].re= (*z).re;
         aa[6*j+i].im=-(*z).im;
         sm+=2.0*((*z).re*(*z).re+(*z).im*(*z).im);
         z+=1;
      }
   }

   eps=DELTA*sm;

   itest=fwd_house(eps);
   solv_sys();
   bck_house();

   sm=0.0;
   u=(*im).u;
   z=(complex_dble*)(u+6);

   for (i=0;i<6;i++)
   {
      *u=aa[6*i+i].re;
      sm+=(*u)*(*u);
      u+=1;

      for (j=i+1;j<6;j++)
      {
         (*z).re=aa[6*i+j].re;
         (*z).im=aa[6*i+j].im;
         sm+=2.0*((*z).re*(*z).re+(*z).im*(*z).im);
         z+=1;
      }
   }

   if ((eps*sm)>1.0)
      itest=1;

   return itest;
}


double det_pauli_dble(pauli_dble *m)
{
   int i,j,k;
   double sm,eps,*u;
   double r1,r2,r3;
   complex_dble det,z,w,*pz;

   sm=0.0;   
   u=(*m).u;
   pz=(complex_dble*)(u+6);

   for (i=0;i<6;i++)
   {
      aa[6*i+i].re=*u;
      aa[6*i+i].im=0.0;
      sm+=(*u)*(*u);
      u+=1;

      for (j=i+1;j<6;j++)
      {
         aa[6*i+j].re= (*pz).re;
         aa[6*i+j].im= (*pz).im;
         aa[6*j+i].re= (*pz).re;
         aa[6*j+i].im=-(*pz).im;
         sm+=2.0*((*pz).re*(*pz).re+(*pz).im*(*pz).im);
         pz+=1;
      }
   }

   eps=DBL_EPSILON*sqrt(sm);
   det.re=1.0;
   det.im=0.0;

   for (k=0;k<5;k++)
   {
      r1=aa[6*k+k].re*aa[6*k+k].re+aa[6*k+k].im*aa[6*k+k].im;
      r2=sqrt(r1);

      for (j=(k+1);j<6;j++)
         r1+=(aa[6*j+k].re*aa[6*j+k].re+aa[6*j+k].im*aa[6*j+k].im);

      r1=sqrt(r1);
      
      if (r1<=eps)
         return 0.0;
      
      if (r2>=(DBL_EPSILON*r1))
      {
         r3=1.0/r2;
         z.re=r1*r3*aa[6*k+k].re;
         z.im=r1*r3*aa[6*k+k].im;
      }
      else
      {
         z.re=r1;
         z.im=0.0;
      }

      w.re=det.re*z.re-det.im*z.im;
      w.im=det.re*z.im+det.im*z.re;
      det.re=w.re;
      det.im=w.im;
      
      aa[6*k+k].re+=z.re;
      aa[6*k+k].im+=z.im;
      r3=1.0/(r1*(r1+r2));

      for (j=(k+1);j<6;j++)
      {
         z.re=0.0;
         z.im=0.0;

         for (i=k;i<6;i++)
         {
            z.re+=(aa[6*i+k].re*aa[6*i+j].re+aa[6*i+k].im*aa[6*i+j].im);
            z.im+=(aa[6*i+k].re*aa[6*i+j].im-aa[6*i+k].im*aa[6*i+j].re);
         }

         z.re*=r3;
         z.im*=r3;

         for (i=(k+1);i<6;i++)
         {
            aa[6*i+j].re-=(z.re*aa[6*i+k].re-z.im*aa[6*i+k].im);
            aa[6*i+j].im-=(z.re*aa[6*i+k].im+z.im*aa[6*i+k].re);
         }
      }
   }

   return det.re*aa[35].re-det.im*aa[35].im;
}


void apply_sw_dble(int vol,pauli_dble *m,spinor_dble *s,spinor_dble *r)
{
   weyl_dble *rs,*rr,*rsm;

   rs=(weyl_dble*)(s);
   rr=(weyl_dble*)(r);   
   rsm=rs+2*vol;

   for (;rs<rsm;)
   {
#if (defined SSE2)
      m+=4;
      _prefetch_pauli_dble(m);
      m+=1;
      _prefetch_pauli_dble(m);
      m-=5;
#endif
      
      mul_pauli_dble(m,rs,rr);

#if (defined SSE2)
      m+=1;
      rs+=4;
      _prefetch_spinor_dble(rs);
      rs-=3;
      rr+=4;
      _prefetch_spinor_dble(rr);
      rr-=3;
#else
      m+=1;
      rs+=1;
      rr+=1;
#endif
      
      mul_pauli_dble(m,rs,rr);      

      m+=1;
      rs+=1;
      rr+=1;
   }
}

