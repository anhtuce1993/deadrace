
/*******************************************************************************
*
* File Qhat_dble.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Application of the even-odd preconditioned Wilson-Dirac operator
*
* The externally accessible functions are
*
*   void Qhat_dble(int k,int l)
*     Applies Qhat to the global double-precision field *psd[k][] and assigns
*     the result to *psd[l][]
*
*   void Qoe_dble(int k,int l)
*     Applies the operator Qoe to the even part of the global double-precision
*     field *psd[k][] and assigns the result to *psd[l][]
*
*   void Qeo_dble(int k,int l)
*     Applies the operator Qeo to the odd part of the global double-precision
*     field *psd[k][] and *subtracts* the result from *psd[l][]
*
*   void Qnohat_dble(int k,int l)
*     Applies the full Dirac operator Q to the global double-precision field
*     *psd[k][] and assigns the result to *psd[l][]
*
*   void Qhat_blk_dble(block_t *b,int k,int l)
*     Applies Qhat to the double-precision field b.sd[k][] on the block b and
*     assigns the result to b.sd[l][]
*
*   void Qoe_blk_dble(block_t *b,int k,int l)
*     Applies the operator Qoe to the even part of the double-precision field
*     b.sd[k][] on the block b and assigns the result to b.sd[l][]
*
*   void Qeo_blk_dble(block_t *b,int k,int l)
*     Applies the operator Qeo to the odd part of the double-precision field
*     b.sd[k][] on the block b and *subtracts* the result from b.sd[l][]
*
*   void Qnohat_blk_dble(block_t *b,int k,int l)
*     Applies the full Dirac operator Q to the double-precision field b.sd[k][]
*     on the block b and assigns the result to b.sd[l][]
*
* Notes:
*
* The notation and normalization conventions are specified in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.ps)
*
* In all these programs it is assumed that Qee and Qoo (or 1/Qoo in the
* case of the preconditioned operator) are stored in the appropriate sw
* arrays
*
* All programs in this file act globally and should be called from all
* processes with the same parameters. Communication buffers are allocated
* automatically when needed
*
*******************************************************************************/

#define QHAT_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "start.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"

static double coe,ceo;
static const spinor_dble sd0={{{0.0}}};

#if (defined SSE2)
#include "sse2.h"

static spinor_dble rs __attribute__ ((aligned (16)));

#define _load_cst(c) \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6")

#define _load_csts(c) \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "movsd %0, %%xmm6 \n\t" \
                      "xorpd %1, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7" \
                      : \
                      : \
                      "m" (c), \
                      "m" (_sse_sgn_dble) \
                      : \
                      "xmm6", "xmm7")

#define _mul_cst() \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm0 \n\t" \
                      "mulpd %%xmm6, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#define _mul_cst_up() \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm6, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#define _mul_csts() \
__asm__ __volatile__ ("mulpd %%xmm7, %%xmm0 \n\t" \
                      "mulpd %%xmm7, %%xmm1 \n\t" \
                      "mulpd %%xmm7, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#define _mul_csts_up() \
__asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t" \
                      "mulpd %%xmm7, %%xmm4 \n\t" \
                      "mulpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")


static void qoe(int *piup,int *pidn,su3_dble *u,spinor_dble *pk)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pk+(*(piup++));

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c3);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);
   _sse_store_up_dble(rs.c1);
   _sse_store_up_dble(rs.c3);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_store_up_dble(rs.c2);
   _sse_store_up_dble(rs.c4);

/******************************* direction -0 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c3);

   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c3);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c3);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;   
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c4);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c4);

/******************************* direction +1 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c4);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c3);

/******************************* direction -1 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c4);

   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble(rs.c4);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble(rs.c3);

/******************************* direction +2 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c4);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c4);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c3);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c3);

/******************************* direction -2 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c4);

   sp=pk+(*(piup));
   _prefetch_spinor_dble(sp);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c4);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c4);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c3);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c3);

/******************************* direction +3 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c3);

   sm=pk+(*(pidn));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.c3);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble(rs.c4);

/******************************* direction -3 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c3);

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _load_csts(coe);
   _sse_load_dble(rs.c1);
   _sse_vector_add_dble();
   _mul_csts();
   _sse_store_dble(rs.c1);

   _sse_load_dble(rs.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _mul_cst();
   _sse_store_dble(rs.c3);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _load_csts(coe);
   _sse_load_dble(rs.c2);
   _sse_vector_add_dble();
   _mul_csts();
   _sse_store_dble(rs.c2);

   _sse_load_dble(rs.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _mul_cst();
   _sse_store_dble(rs.c4);
}


static void qeo(int *piup,int *pidn,su3_dble *u,spinor_dble *pl)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c3);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c3);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c3);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c4);

   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c4);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

/******************************* direction -0 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c3);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c3);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c3);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c4);

   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c4);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c4);

/******************************* direction +1 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c4);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;   
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c4);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c3);

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c3);

/******************************* direction -1 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c4);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c4);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c3);

   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

/******************************* direction +2 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c4);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c4);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c3);

   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c3);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c3);

/******************************* direction -2 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c4);

   sp=pl+(*(piup));
   _prefetch_spinor_dble(sp);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c4);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c4);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c3);

   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c3);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

/******************************* direction +3 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c3);

   sm=pl+(*(pidn));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c3);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c4);

   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

/******************************* direction -3 *********************************/

   _sse_load_dble(rs.c1);
   _sse_load_up_dble(rs.c3);

   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

   _sse_load_dble(rs.c2);
   _sse_load_up_dble(rs.c4);

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c4);
}

#else

#define _vector_mul_assign(r,c) \
   (r).c1.re*=(c); \
   (r).c1.im*=(c); \
   (r).c2.re*=(c); \
   (r).c2.im*=(c); \
   (r).c3.re*=(c); \
   (r).c3.im*=(c)

static spinor_dble rs;


static void qoe(int *piup,int *pidn,su3_dble *u,spinor_dble *pk)
{
   spinor_dble *sp,*sm;
   su3_vector_dble psi,chi;

/******************************* direction +0 *********************************/

      sp=pk+(*(piup++));

      _vector_add(psi,(*sp).c1,(*sp).c3);
      _su3_multiply(rs.c1,*u,psi);
      rs.c3=rs.c1;

      _vector_add(psi,(*sp).c2,(*sp).c4);
      _su3_multiply(rs.c2,*u,psi);
      rs.c4=rs.c2;

/******************************* direction -0 *********************************/

      sm=pk+(*(pidn++));
      u+=1;

      _vector_sub(psi,(*sm).c1,(*sm).c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_sub_assign(rs.c3,chi);

      _vector_sub(psi,(*sm).c2,(*sm).c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_sub_assign(rs.c4,chi);

/******************************* direction +1 *********************************/

      sp=pk+(*(piup++));
      u+=1;

      _vector_i_add(psi,(*sp).c1,(*sp).c4);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_i_sub_assign(rs.c4,chi);

      _vector_i_add(psi,(*sp).c2,(*sp).c3);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_i_sub_assign(rs.c3,chi);

/******************************* direction -1 *********************************/

      sm=pk+(*(pidn++));
      u+=1;

      _vector_i_sub(psi,(*sm).c1,(*sm).c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_i_add_assign(rs.c4,chi);

      _vector_i_sub(psi,(*sm).c2,(*sm).c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_i_add_assign(rs.c3,chi);

/******************************* direction +2 *********************************/

      sp=pk+(*(piup++));
      u+=1;

      _vector_add(psi,(*sp).c1,(*sp).c4);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_add_assign(rs.c4,chi);

      _vector_sub(psi,(*sp).c2,(*sp).c3);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_sub_assign(rs.c3,chi);

/******************************* direction -2 *********************************/

      sm=pk+(*(pidn++));
      u+=1;

      _vector_sub(psi,(*sm).c1,(*sm).c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_sub_assign(rs.c4,chi);

      _vector_add(psi,(*sm).c2,(*sm).c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_add_assign(rs.c3,chi);

/******************************* direction +3 *********************************/

      sp=pk+(*(piup));
      u+=1;

      _vector_i_add(psi,(*sp).c1,(*sp).c3);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_i_sub_assign(rs.c3,chi);

      _vector_i_sub(psi,(*sp).c2,(*sp).c4);
      _su3_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_i_add_assign(rs.c4,chi);

/******************************* direction -3 *********************************/

      sm=pk+(*(pidn));
      u+=1;

      _vector_i_sub(psi,(*sm).c1,(*sm).c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c1,chi);
      _vector_i_add_assign(rs.c3,chi);

      _vector_i_add(psi,(*sm).c2,(*sm).c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_add_assign(rs.c2,chi);
      _vector_i_sub_assign(rs.c4,chi);

      _vector_mul_assign(rs.c1,-coe);
      _vector_mul_assign(rs.c2,-coe);
      _vector_mul_assign(rs.c3,coe);
      _vector_mul_assign(rs.c4,coe);
}


static void qeo(int *piup,int *pidn,su3_dble *u,spinor_dble *pl)
{
   spinor_dble *sp,*sm;
   su3_vector_dble psi,chi;

/******************************* direction +0 *********************************/

      sp=pl+(*(piup++));

      _vector_sub(psi,rs.c1,rs.c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c1,chi);
      _vector_sub_assign((*sp).c3,chi);

      _vector_sub(psi,rs.c2,rs.c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c2,chi);
      _vector_sub_assign((*sp).c4,chi);

/******************************* direction -0 *********************************/

      sm=pl+(*(pidn++));
      u+=1;

      _vector_add(psi,rs.c1,rs.c3);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c1,chi);
      _vector_add_assign((*sm).c3,chi);

      _vector_add(psi,rs.c2,rs.c4);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c2,chi);
      _vector_add_assign((*sm).c4,chi);

/******************************* direction +1 *********************************/

      sp=pl+(*(piup++));
      u+=1;

      _vector_i_sub(psi,rs.c1,rs.c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c1,chi);
      _vector_i_add_assign((*sp).c4,chi);

      _vector_i_sub(psi,rs.c2,rs.c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c2,chi);
      _vector_i_add_assign((*sp).c3,chi);

/******************************* direction -1 *********************************/

      sm=pl+(*(pidn++));
      u+=1;

      _vector_i_add(psi,rs.c1,rs.c4);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c1,chi);
      _vector_i_sub_assign((*sm).c4,chi);

      _vector_i_add(psi,rs.c2,rs.c3);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c2,chi);
      _vector_i_sub_assign((*sm).c3,chi);

/******************************* direction +2 *********************************/

      sp=pl+(*(piup++));
      u+=1;

      _vector_sub(psi,rs.c1,rs.c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c1,chi);
      _vector_sub_assign((*sp).c4,chi);

      _vector_add(psi,rs.c2,rs.c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c2,chi);
      _vector_add_assign((*sp).c3,chi);

/******************************* direction -2 *********************************/

      sm=pl+(*(pidn++));
      u+=1;

      _vector_add(psi,rs.c1,rs.c4);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c1,chi);
      _vector_add_assign((*sm).c4,chi);

      _vector_sub(psi,rs.c2,rs.c3);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c2,chi);
      _vector_sub_assign((*sm).c3,chi);

/******************************* direction +3 *********************************/

      sp=pl+(*(piup));
      u+=1;

      _vector_i_sub(psi,rs.c1,rs.c3);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c1,chi);
      _vector_i_add_assign((*sp).c3,chi);

      _vector_i_add(psi,rs.c2,rs.c4);
      _su3_inverse_multiply(chi,*u,psi);
      _vector_sub_assign((*sp).c2,chi);
      _vector_i_sub_assign((*sp).c4,chi);

/******************************* direction -3 *********************************/

      sm=pl+(*(pidn));
      u+=1;

      _vector_i_add(psi,rs.c1,rs.c3);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c1,chi);
      _vector_i_sub_assign((*sm).c3,chi);

      _vector_i_sub(psi,rs.c2,rs.c4);
      _su3_multiply(chi,*u,psi);
      _vector_sub_assign((*sm).c2,chi);
      _vector_i_add_assign((*sm).c4,chi);
}

#endif

void Qhat_dble(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3_dble *u,*um;
   pauli_dble *m;
   spinor_dble *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qhat_dble [Qhat_dble.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_sd)||(l<0)||(l>=no_sd)||(k==l),1,
              "Qhat_dble [Qhat_dble.c]","Improper arguments");

   pk=psd[k][0];
   pl=psd[l][0];

   apply_sw_dble(VOLUME/2,swd,pk,pl);

#if (NPROC>1)
   set_sd2zero(BNDRY/2,pl+VOLUME);
   cpsd_int_bnd(k);
#endif

   coe=-0.25;
   piup=&iup[(VOLUME/2)][0];
   pidn=&idn[(VOLUME/2)][0];
   m=swd+VOLUME;
   u=pud[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE2)
      _prefetch_pauli_dble(m);
      qoe(piup,pidn,u,pk);
      m+=1;
      _prefetch_pauli_dble(m);
      m-=1;
#else
      qoe(piup,pidn,u,pk);
#endif
      mul_pauli_dble(m,(weyl_dble*)(&rs.c1),(weyl_dble*)(&rs.c1));
      mul_pauli_dble(m+1,(weyl_dble*)(&rs.c1)+1,(weyl_dble*)(&rs.c1)+1);
      qeo(piup,pidn,u,pl);

      piup+=4;
      pidn+=4;
      m+=2;
   }

#if (NPROC>1)
   cpsd_ext_bnd(l);
#endif
}


void Qoe_dble(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3_dble *u,*um;
   spinor_dble *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qoe_dble [Qhat_dble.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_sd)||(l<0)||(l>=no_sd),1,
              "Qoe_dble [Qhat_dble.c]","Improper arguments");

   pk=psd[k][0];
   pl=psd[l][VOLUME/2];

#if (NPROC>1)
   cpsd_int_bnd(k);
#endif

   coe=0.5;
   piup=&iup[VOLUME/2][0];
   pidn=&idn[VOLUME/2][0];
   u=pud[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
      qoe(piup,pidn,u,pk);

#if (defined SSE2)
      _sse_load_dble(rs.c1);
      _sse_load_up_dble(rs.c2);
      _sse_store_dble((*pl).c1);
      _sse_store_up_dble((*pl).c2);
      _sse_load_dble(rs.c3);
      _sse_load_up_dble(rs.c4);
      _sse_store_dble((*pl).c3);
      _sse_store_up_dble((*pl).c4);
#else
      *pl=rs;
#endif
      
      piup+=4;
      pidn+=4;
      pl+=1;
   }
}


void Qeo_dble(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3_dble *u,*um;
   spinor_dble *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qeo_dble [Qhat_dble.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_sd)||(l<0)||(l>=no_sd),1,
              "Qeo_dble [Qhat_dble.c]","Improper arguments");

   pk=psd[k][VOLUME/2];
   pl=psd[l][0];

#if (NPROC>1)
   set_sd2zero(BNDRY/2,pl+VOLUME);
#endif

   ceo=-0.5;
   piup=&iup[VOLUME/2][0];
   pidn=&idn[VOLUME/2][0];
   u=pud[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE2)
      _load_cst(ceo);
      _sse_load_dble((*pk).c1);
      _sse_load_up_dble((*pk).c2);
      _mul_cst();
      _mul_cst_up();
      _sse_store_dble(rs.c1);
      _sse_store_up_dble(rs.c2);

      _sse_load_dble((*pk).c3);
      _sse_load_up_dble((*pk).c4);
      _mul_cst();
      _mul_cst_up();
      _sse_store_dble(rs.c3);
      _sse_store_up_dble(rs.c4);
#else
      _vector_mul(rs.c1,ceo,(*pk).c1);
      _vector_mul(rs.c2,ceo,(*pk).c2);
      _vector_mul(rs.c3,ceo,(*pk).c3);
      _vector_mul(rs.c4,ceo,(*pk).c4);
#endif

      qeo(piup,pidn,u,pl);

      piup+=4;
      pidn+=4;
      pk+=1;
   }

#if (NPROC>1)
   cpsd_ext_bnd(l);
#endif
}


void Qnohat_dble(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3_dble *u,*um;
   pauli_dble *m;
   spinor_dble *pke,*pko,*ple,*plo;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qnohat_dble [Qhat_dble.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_sd)||(l<0)||(l>=no_sd)||(k==l),1,
              "Qnohat_dble [Qhat_dble.c]","Improper arguments");

   pke=psd[k][0];
   pko=psd[k][VOLUME/2];
   ple=psd[l][0];
   plo=psd[l][VOLUME/2];

   apply_sw_dble(VOLUME/2,swd,pke,ple);

#if (NPROC>1)
   set_sd2zero(BNDRY/2,ple+VOLUME);
   cpsd_int_bnd(k);
#endif

   coe=0.5;
   ceo=0.5;
   piup=&iup[VOLUME/2][0];
   pidn=&idn[VOLUME/2][0];
   m=swd+VOLUME;
   u=pud[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE2)
      _prefetch_pauli_dble(m);
#endif      
      qoe(piup,pidn,u,pke);
#if (defined SSE2)
      m+=1;
      _prefetch_pauli_dble(m);
      m-=1;
#endif
      mul_pauli_dble(m,(weyl_dble*)(pko),(weyl_dble*)(plo));
      mul_pauli_dble(m+1,(weyl_dble*)(pko)+1,(weyl_dble*)(plo)+1);
      
#if (defined SSE2)      
      _sse_load_dble((*plo).c1);
      _sse_load_up_dble((*plo).c2);

      __asm__ __volatile__ ("addpd %0, %%xmm0 \n\t"
                            "addpd %1, %%xmm1 \n\t"
                            "addpd %2, %%xmm2 \n\t"
                            "addpd %3, %%xmm3 \n\t"
                            "addpd %4, %%xmm4 \n\t"
                            "addpd %5, %%xmm5"
                            : \
                            :
                            "m" (rs.c1.c1.re),
                            "m" (rs.c1.c2.re),
                            "m" (rs.c1.c3.re),
                            "m" (rs.c2.c1.re),
                            "m" (rs.c2.c2.re),
                            "m" (rs.c2.c3.re)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
                            
      _sse_store_dble((*plo).c1);
      _sse_store_up_dble((*plo).c2);

      _sse_load_dble((*plo).c3);
      _sse_load_up_dble((*plo).c4);

      __asm__ __volatile__ ("addpd %0, %%xmm0 \n\t"
                            "addpd %1, %%xmm1 \n\t"
                            "addpd %2, %%xmm2 \n\t"
                            "addpd %3, %%xmm3 \n\t"
                            "addpd %4, %%xmm4 \n\t"
                            "addpd %5, %%xmm5"
                            : \
                            :
                            "m" (rs.c3.c1.re),
                            "m" (rs.c3.c2.re),
                            "m" (rs.c3.c3.re),
                            "m" (rs.c4.c1.re),
                            "m" (rs.c4.c2.re),
                            "m" (rs.c4.c3.re)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
                            
      _sse_store_dble((*plo).c3);
      _sse_store_up_dble((*plo).c4);

      _load_cst(ceo);
      _sse_load_dble((*pko).c1);
      _sse_load_up_dble((*pko).c2);
      _mul_cst();
      _mul_cst_up();
      _sse_store_dble(rs.c1);
      _sse_store_up_dble(rs.c2);

      _sse_load_dble((*pko).c3);
      _sse_load_up_dble((*pko).c4);
      _mul_cst();
      _mul_cst_up();
      _sse_store_dble(rs.c3);
      _sse_store_up_dble(rs.c4);
#else
      _vector_add_assign((*plo).c1,rs.c1);
      _vector_add_assign((*plo).c2,rs.c2);
      _vector_add_assign((*plo).c3,rs.c3);
      _vector_add_assign((*plo).c4,rs.c4);

      _vector_mul(rs.c1,ceo,(*pko).c1);
      _vector_mul(rs.c2,ceo,(*pko).c2);
      _vector_mul(rs.c3,ceo,(*pko).c3);
      _vector_mul(rs.c4,ceo,(*pko).c4);
#endif
      qeo(piup,pidn,u,ple);

      piup+=4;
      pidn+=4;
      pko+=1;
      plo+=1;
      m+=2;
   }

#if (NPROC>1)
   cpsd_ext_bnd(l);
#endif
}
