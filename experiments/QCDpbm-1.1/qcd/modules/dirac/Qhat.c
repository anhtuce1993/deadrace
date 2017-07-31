
/*******************************************************************************
*
* File Qhat.c
*
* Copyright (C) 2005, 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Application of the even-odd preconditioned Wilson-Dirac operator
*
* The externally accessible functions are
*
*   void Qhat(int k,int l)
*     Applies Qhat to the global single-precision field *ps[k][] and assigns
*     the result to *ps[l][]
*
*   void Qoe(int k,int l)
*     Applies the operator Qoe to the even part of the global single-precision
*     field *ps[k][] and assigns the result to *ps[l][]
*
*   void Qeo(int k,int l)
*     Applies the operator Qeo to the odd part of the global single-precision
*     field *ps[k][] and *subtracts* the result from *ps[l][]
*
*   void Qnohat(int k,int l)
*     Applies the full Dirac operator Q to the global single-precision field
*     *ps[k][] and assigns the result to *ps[l][]
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
* All programs in this file act globally and must be called simultaneously
* from all processes with the same parameters. Communication buffers are
* allocated automatically when needed
*
*******************************************************************************/

#define QHAT_C

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

static float coe,ceo;
static const spinor s0={{{0.0f}}};

#if (defined SSE)
#include "sse.h"

static spinor rs __attribute__ ((aligned (16)));
static sse_vector r12,r34;

#define _load_cst(c) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6")

#define _load_csts(c) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm7 \n\t" \
                      "mulps %1, %%xmm7" \
                      : \
                      : \
                      "m" (c), \
                      "m" (_sse_sgn) \
                      : \
                      "xmm6", "xmm7")

#define _mul_cst() \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm6, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#define _mul_cst_up() \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm6, %%xmm4 \n\t" \
                      "mulps %%xmm6, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

#define _mul_csts() \
__asm__ __volatile__ ("mulps %%xmm7, %%xmm0 \n\t" \
                      "mulps %%xmm7, %%xmm1 \n\t" \
                      "mulps %%xmm7, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")


static void qoe(int *piup,int *pidn,su3 *u,spinor *pk)
{
   spinor *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pk+(*(piup++));

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_pair_load_up((*sp).c3,(*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor(sm);

   _sse_vector_add();
   sp=pk+(*(piup++));
   _prefetch_spinor(sp);
   _sse_su3_multiply(*u);

   _sse_vector_store_up(r12);
   _sse_vector_store_up(r34);

/******************************* direction -0 *********************************/

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_pair_load_up((*sm).c3,(*sm).c4);

   u+=2;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_sub();
   sm=pk+(*(pidn++));
   _prefetch_spinor(sm);
   _sse_su3_inverse_multiply(*u);

   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_sub();
   _sse_vector_store(r34);

/******************************* direction +1 *********************************/

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_pair_load_up((*sp).c4,(*sp).c3);

   _sse_vector_i_add();
   sp=pk+(*(piup++));
   _prefetch_spinor(sp);
   u+=1;
   _sse_su3_multiply(*u);

   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_xch_i_sub();
   _sse_vector_store(r34);

/******************************* direction -1 *********************************/

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_pair_load_up((*sm).c4,(*sm).c3);

   u+=2;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_sub();
   sm=pk+(*(pidn++));
   _prefetch_spinor(sm);
   _sse_su3_inverse_multiply(*u);

   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_xch_i_add();
   _sse_vector_store(r34);

/******************************* direction +2 *********************************/

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_pair_load_up((*sp).c4,(*sp).c3);

   _sse_vector_addsub();

   u+=1;
   _sse_su3_multiply(*u);
   sp=pk+(*(piup));
   _prefetch_spinor(sp);
   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_xch();
   _sse_vector_subadd();
   _sse_vector_store(r34);

/******************************* direction -2 *********************************/

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_pair_load_up((*sm).c4,(*sm).c3);

   u+=2;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_subadd();
   sm=pk+(*(pidn));
   _prefetch_spinor(sm);
   _sse_su3_inverse_multiply(*u);

   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_xch();
   _sse_vector_addsub();
   _sse_vector_store(r34);

/******************************* direction +3 *********************************/

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_pair_load_up((*sp).c3,(*sp).c4);

   _sse_vector_i_addsub();
   u+=1;
   _sse_su3_multiply(*u);

   _sse_vector_load(r12);
   _sse_vector_add();
   _sse_vector_store(r12);

   _sse_vector_load(r34);
   _sse_vector_i_subadd();
   _sse_vector_store(r34);

/******************************* direction -3 *********************************/

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_pair_load_up((*sm).c3,(*sm).c4);

   u+=2;
   _prefetch_su3_dble(u);
   u-=1;
   _sse_vector_i_subadd();
   _sse_su3_inverse_multiply(*u);

   _load_csts(coe);
   _sse_vector_load(r12);
   _sse_vector_add();
   _mul_csts();
   _sse_pair_store(rs.c1,rs.c2);

   _sse_vector_load(r34);
   _sse_vector_i_addsub();
   _mul_cst();
   _sse_pair_store(rs.c3,rs.c4);
}


static void qeo(int *piup,int *pidn,su3 *u,spinor *pl)
{
   spinor *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pl+(*(piup++));
   _prefetch_spinor(sp);

   _sse_pair_load(rs.c1,rs.c2);
   _sse_pair_load_up(rs.c3,rs.c4);
   _sse_vector_store(r12);
   _sse_vector_store_up(r34);

   sm=pl+(*(pidn++));
   _prefetch_spinor(sm);
   _sse_vector_sub();
   _sse_su3_inverse_multiply(*u);

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_vector_sub();
   _sse_pair_store((*sp).c1,(*sp).c2);

   _sse_pair_load((*sp).c3,(*sp).c4);
   _sse_vector_sub();
   _sse_pair_store((*sp).c3,(*sp).c4);

/******************************* direction -0 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sp=pl+(*(piup++));
   _prefetch_spinor(sp);
   _sse_vector_add();
   u+=1;
   _sse_su3_multiply(*u);

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_vector_sub();
   _sse_pair_store((*sm).c1,(*sm).c2);

   _sse_pair_load((*sm).c3,(*sm).c4);
   _sse_vector_add();
   _sse_pair_store((*sm).c3,(*sm).c4);

/******************************* direction +1 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sm=pl+(*(pidn++));
   _prefetch_spinor(sm);
   _sse_vector_xch_i_sub();
   u+=1;
   _sse_su3_inverse_multiply(*u);

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_vector_sub();
   _sse_pair_store((*sp).c1,(*sp).c2);

   _sse_pair_load((*sp).c3,(*sp).c4);
   _sse_vector_xch_i_add();
   _sse_pair_store((*sp).c3,(*sp).c4);

/******************************* direction -1 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sp=pl+(*(piup++));
   _prefetch_spinor(sp);
   _sse_vector_xch_i_add();
   u+=1;
   _sse_su3_multiply(*u);

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_vector_sub();
   _sse_pair_store((*sm).c1,(*sm).c2);

   _sse_pair_load((*sm).c3,(*sm).c4);
   _sse_vector_xch_i_sub();
   _sse_pair_store((*sm).c3,(*sm).c4);

/******************************* direction +2 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sm=pl+(*(pidn++));
   _prefetch_spinor(sm);
   _sse_vector_xch();
   _sse_vector_subadd();
   u+=1;
   _sse_su3_inverse_multiply(*u);

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_vector_sub();
   _sse_pair_store((*sp).c1,(*sp).c2);

   _sse_pair_load((*sp).c3,(*sp).c4);
   _sse_vector_xch();
   _sse_vector_addsub();
   _sse_pair_store((*sp).c3,(*sp).c4);

/******************************* direction -2 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sp=pl+(*(piup));
   _prefetch_spinor(sp);
   _sse_vector_xch();
   _sse_vector_addsub();
   u+=1;   
   _sse_su3_multiply(*u);

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_vector_sub();
   _sse_pair_store((*sm).c1,(*sm).c2);

   _sse_pair_load((*sm).c3,(*sm).c4);
   _sse_vector_xch();
   _sse_vector_subadd();
   _sse_pair_store((*sm).c3,(*sm).c4);

/******************************* direction +3 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   sm=pl+(*(pidn));
   _prefetch_spinor(sm);
   _sse_vector_i_subadd();
   u+=1;   
   _sse_su3_inverse_multiply(*u);

   _sse_pair_load((*sp).c1,(*sp).c2);
   _sse_vector_sub();
   _sse_pair_store((*sp).c1,(*sp).c2);

   _sse_pair_load((*sp).c3,(*sp).c4);
   _sse_vector_i_addsub();
   _sse_pair_store((*sp).c3,(*sp).c4);

/******************************* direction -3 *********************************/

   _sse_vector_load(r12);
   _sse_vector_load_up(r34);

   _sse_vector_i_addsub();
   u+=1;
   _sse_su3_multiply(*u);

   _sse_pair_load((*sm).c1,(*sm).c2);
   _sse_vector_sub();
   _sse_pair_store((*sm).c1,(*sm).c2);

   _sse_pair_load((*sm).c3,(*sm).c4);
   _sse_vector_i_subadd();
   _sse_pair_store((*sm).c3,(*sm).c4);
}

#else

#define _vector_mul_assign(r,c) \
   (r).c1.re*=(c); \
   (r).c1.im*=(c); \
   (r).c2.re*=(c); \
   (r).c2.im*=(c); \
   (r).c3.re*=(c); \
   (r).c3.im*=(c)

static spinor rs;


static void qoe(int *piup,int *pidn,su3 *u,spinor *pk)
{
   spinor *sp,*sm;
   su3_vector psi,chi;

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


static void qeo(int *piup,int *pidn,su3 *u,spinor *pl)
{
   spinor *sp,*sm;
   su3_vector psi,chi;

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


void Qhat(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3 *u,*um;
   pauli *m;
   spinor *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qhat [Qhat.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_s)||(l<0)||(l>=no_s)||(k==l),1,
              "Qhat [Qhat.c]","Improper arguments");

   pk=ps[k][0];
   pl=ps[l][0];

   apply_sw(VOLUME/2,sw,pk,pl);

#if (NPROC>1)
   set_s2zero(BNDRY/2,pl+VOLUME);
   cps_int_bnd(k);
#endif

   coe=-0.25f;
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];
   m=sw+VOLUME;
   u=pu[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE)
      _prefetch_pauli_dble(m);
#endif
      qoe(piup,pidn,u,pk);
      mul_pauli(m,(weyl*)(&rs.c1),(weyl*)(&rs.c1));
      mul_pauli(m+1,(weyl*)(&rs.c1)+1,(weyl*)(&rs.c1)+1);
      qeo(piup,pidn,u,pl);

      piup+=4;
      pidn+=4;
      m+=2;
   }

#if (NPROC>1)
   cps_ext_bnd(l);
#endif
}


void Qoe(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3 *u,*um;
   spinor *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qoe [Qhat.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_s)||(l<0)||(l>=no_s),1,
              "Qoe [Qhat.c]","Improper arguments");

   pk=ps[k][0];
   pl=ps[l][VOLUME/2];

#if (NPROC>1)
   cps_int_bnd(k);
#endif

   coe=0.5f;
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];
   u=pu[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
      qoe(piup,pidn,u,pk);

#if (defined SSE)
      _sse_spinor_load(rs);
      _sse_spinor_store(*pl);
#else
      *pl=rs;
#endif

      piup+=4;
      pidn+=4;
      pl+=1;
   }
}


void Qeo(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3 *u,*um;
   spinor *pk,*pl;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qeo [Qhat.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_s)||(l<0)||(l>=no_s),1,
              "Qeo [Qhat.c]","Improper arguments");

   pk=ps[k][VOLUME/2];
   pl=ps[l][0];

#if (NPROC>1)
   set_s2zero(BNDRY/2,pl+VOLUME);
#endif

   ceo=-0.5f;
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];
   u=pu[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE)
      _load_cst(ceo);
      _sse_spinor_load(*pk);
      _mul_cst();
      _mul_cst_up();
      _sse_spinor_store(rs);
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
   cps_ext_bnd(l);
#endif
}


void Qnohat(int k,int l)
{
   int *piup,*pidn,iprms[2];
   su3 *u,*um;
   pauli *m;
   spinor *pke,*pko,*ple,*plo;

   iprms[0]=k;
   iprms[1]=l;

   MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=k)||(iprms[1]!=l),1,
         "Qnohat [Qhat.c]","Parameters k,l are not global");

   error_root((k<0)||(k>=no_s)||(l<0)||(l>=no_s)||(k==l),1,
              "Qnohat [Qhat.c]","Improper arguments");

   pke=ps[k][0];
   pko=ps[k][VOLUME/2];
   ple=ps[l][0];
   plo=ps[l][VOLUME/2];

   apply_sw(VOLUME/2,sw,pke,ple);

#if (NPROC>1)
   set_s2zero(BNDRY/2,ple+VOLUME);
   cps_int_bnd(k);
#endif

   coe=0.5f;
   ceo=0.5f;

   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];
   m=sw+VOLUME;
   u=pu[VOLUME/2][0];
   um=u+4*VOLUME;

   for (;u<um;u+=8)
   {
#if (defined SSE)
      _prefetch_pauli_dble(m);
#endif
      qoe(piup,pidn,u,pke);
      mul_pauli(m,(weyl*)(pko),(weyl*)(plo));
      mul_pauli(m+1,(weyl*)(pko)+1,(weyl*)(plo)+1);

#if (defined SSE)      
      _sse_spinor_load(*plo);

      __asm__ __volatile__ ("addps %0, %%xmm0 \n\t"
                            "addps %1, %%xmm1 \n\t"
                            "addps %2, %%xmm2 \n\t"
                            "addps %3, %%xmm3 \n\t"
                            "addps %4, %%xmm4 \n\t"
                            "addps %5, %%xmm5"
                            :
                            :
                            "m" (rs.c1.c1.re),
                            "m" (rs.c1.c3.re),
                            "m" (rs.c2.c2.re),
                            "m" (rs.c3.c1.re),
                            "m" (rs.c3.c3.re),
                            "m" (rs.c4.c2.re)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");

      _sse_spinor_store(*plo);

      _load_cst(ceo);
      _sse_spinor_load(*pko);
      _mul_cst();
      _mul_cst_up();
      _sse_spinor_store(rs);
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
   cps_ext_bnd(l);
#endif
}


