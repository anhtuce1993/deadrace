
/*******************************************************************************
*
* File sse.h
*
* Copyright (C) 2005, 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using
* inline assembly SSE instructions
*
*******************************************************************************/

#ifndef SSE_H
#define SSE_H

typedef struct
{
   int c1,c2,c3,c4;
} sse_int __attribute__ ((aligned (16)));

typedef struct
{
   float c1,c2,c3,c4;
} sse_float __attribute__ ((aligned (16)));

typedef struct
{
   sse_float c1,c2,c3;
} sse_vector __attribute__ ((aligned (16)));

static sse_float _sse_sgn12 __attribute__ ((unused)) ={-1.0f,-1.0f,1.0f,1.0f};
static sse_float _sse_sgn13 __attribute__ ((unused)) ={-1.0f,1.0f,-1.0f,1.0f};
static sse_float _sse_sgn14 __attribute__ ((unused)) ={-1.0f,1.0f,1.0f,-1.0f};
static sse_float _sse_sgn23 __attribute__ ((unused)) ={1.0f,-1.0f,-1.0f,1.0f};
static sse_float _sse_sgn24 __attribute__ ((unused)) ={1.0f,-1.0f,1.0f,-1.0f};
static sse_float _sse_sgn34 __attribute__ ((unused)) ={1.0f,1.0f,-1.0f,-1.0f};
static sse_float _sse_sgn   __attribute__ ((unused)) ={-1.0f,-1.0f,-1.0f,-1.0f};

/*******************************************************************************
*
* Prefetch macros
*
*******************************************************************************/

#if (defined P4)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x7fL))

#define _prefetch_128b(addr) \
__asm__ __volatile__ ("prefetcht0 %0"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))))

#define _prefetch_256b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_384b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_128b((addr))

#define _prefetch_weyl(addr) \
_prefetch_256b((addr))

#define _prefetch_spinor(addr) \
_prefetch_256b((addr))

#define _prefetch_su3(addr) \
_prefetch_256b((addr))

#define _prefetch_pauli(addr) \
_prefetch_256b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_384b((addr))

#elif (defined PM)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x3fL))

#define _prefetch_64b(addr) \
__asm__ __volatile__ ("prefetcht0 %0"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))))

#define _prefetch_128b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))))
                      
#define _prefetch_192b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_320b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0xc0L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_64b((addr))

#define _prefetch_weyl(addr) \
_prefetch_128b((addr))

#define _prefetch_spinor(addr) \
_prefetch_128b((addr))

#define _prefetch_su3(addr) \
_prefetch_128b((addr))

#define _prefetch_pauli(addr) \
_prefetch_192b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_128b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_320b((addr))

#elif (defined P3)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x1fL))

#define _prefetch_64b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))))

#define _prefetch_96b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))))

#define _prefetch_160b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_192b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4 \n\t" \
                      "prefetcht0 %5"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0xa0L))))

#define _prefetch_288b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L)))); \
__asm__ __volatile__ ("prefetcht0 %0 \n\t"  \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3"  \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)+0xa0L))), \
                      "m" (*((char*)(_pfbase(addr)+0xc0L))), \
                      "m" (*((char*)(_pfbase(addr)+0xe0L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_64b((addr))

#define _prefetch_weyl(addr) \
_prefetch_64b((addr))

#define _prefetch_spinor(addr) \
_prefetch_96b((addr))

#define _prefetch_su3(addr) \
_prefetch_96b((addr))

#define _prefetch_pauli(addr) \
_prefetch_160b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_96b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_160b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_288b((addr))

#else

#define _prefetch_su3_alg_dble(addr) \

#define _prefetch_weyl(addr)

#define _prefetch_spinor(addr)

#define _prefetch_su3(addr)

#define _prefetch_pauli(addr)

#define _prefetch_weyl_dble(addr)

#define _prefetch_spinor_dble(addr)

#define _prefetch_su3_dble(addr)

#define _prefetch_pauli_dble(addr)

#endif

/*******************************************************************************
*
* Macros for su3 vectors
*
* Most of these macros operate on pairs of su3 vectors that are stored
* in the low and high words of xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> sl.c1.re,sl.c1.im,sh.c1.re,sh.c1.im
* xmm1 -> sl.c2.re,sl.c2.im,sh.c2.re,sh.c2.im
* xmm2 -> sl.c3.re,sl.c3.im,sh.c3.re,sh.c3.im
*
* (where sl and sh are of type su3_vector). This can also be interpreted as
* an sse_vector s that is stored in these registers according to
*
* xmm0 -> s.c1.c1,s.c1.c2,s.c1.c3,s.c1.c4
* xmm1 -> s.c2.c1,s.c2.c2,s.c2.c3,s.c2.c4
* xmm2 -> s.c3.c1,s.c3.c2,s.c3.c3,s.c3.c4
*
* The load and store macros can be used to move data in either format
* from and to the xmm registers
*
*******************************************************************************/

/*
* Loads two su3 vectors sl and sh to the low and high words of xmm0,xmm1,xmm2
*/

#if defined SSE2

#define _sse_pair_load(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

#else

#define _sse_pair_load(sl,sh) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

#endif

/*
* Loads two su3 vectors sl and sh to the low and high words of xmm3,xmm4,xmm5
*/

#if defined SSE2

#define _sse_pair_load_up(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm4 \n\t" \
                      "movsd %2, %%xmm5 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

#else

#define _sse_pair_load_up(sl,sh) \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movlps %1, %%xmm4 \n\t" \
                      "movlps %2, %%xmm5 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Stores the low and high words of xmm0,xmm1,xmm2 to the su3 vectors rl and rh
*/

#define _sse_pair_store(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      "movlps %%xmm2, %2 \n\t" \
                      "movhps %%xmm0, %3 \n\t" \
                      "movhps %%xmm1, %4 \n\t" \
                      "movhps %%xmm2, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Stores the low and high words of xmm3,xmm4,xmm5 to the su3 vectors rl and rh
*/

#define _sse_pair_store_up(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Loads the components s.c1,s.c2,s.c3 of an sse_vector s to xmm0,xmm1,xmm2
*/

#define _sse_vector_load(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %1, %%xmm1 \n\t" \
                      "movaps %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads the components s.c1,s.c2,s.c3 of an sse_vector s to xmm3,xmm4,xmm5
*/

#define _sse_vector_load_up(s) \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %1, %%xmm4 \n\t" \
                      "movaps %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an sse_vector r
*/

#define _sse_vector_store(r) \
__asm__ __volatile__ ("movaps %%xmm0, %0 \n\t" \
                      "movaps %%xmm1, %1 \n\t" \
                      "movaps %%xmm2, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an sse_vector r
*/

#define _sse_vector_store_up(r) \
__asm__ __volatile__ ("movaps %%xmm3, %0 \n\t" \
                      "movaps %%xmm4, %1 \n\t" \
                      "movaps %%xmm5, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Multiplies xmm0,xmm1,xmm2 with a constant sse_float c
*/

#define _sse_vector_mul(c) \
__asm__ __volatile__ ("mulps %0, %%xmm0 \n\t" \
                      "mulps %0, %%xmm1 \n\t" \
                      "mulps %0, %%xmm2" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies xmm3,xmm4,xmm5 with a constant sse_float c
*/

#define _sse_vector_mul_up(c) \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Adds xmm3,xmm4,xmm5 to xmm0,xmm1,xmm2
*/

#define _sse_vector_add() \
__asm__ __volatile__ ("addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Subtracts xmm3,xmm4,xmm5 from xmm0,xmm1,xmm2
*/

#define _sse_vector_sub() \
__asm__ __volatile__ ("subps %%xmm3, %%xmm0 \n\t" \
                      "subps %%xmm4, %%xmm1 \n\t" \
                      "subps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies the high words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_addsub() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn34) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_subadd() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn12) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies xmm3,xmm4,xmm5 with i and adds them to xmm0,xmm1,xmm2
*/

#if (defined SSE3)

#define _sse_vector_i_add() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "addsubps %%xmm3, %%xmm0 \n\t" \
                      "addsubps %%xmm4, %%xmm1 \n\t" \
                      "addsubps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _sse_vector_i_add() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Multiplies xmm3,xmm4,xmm5 with i and subtracts them from xmm0,xmm1,xmm2
*/

#define _sse_vector_i_sub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and adds the result to xmm0,xmm1,xmm2
*/

#if (defined SSE3)

#define _sse_vector_xch_i_add() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "addsubps %%xmm3, %%xmm0 \n\t" \
                      "addsubps %%xmm4, %%xmm1 \n\t" \
                      "addsubps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#else

#define _sse_vector_xch_i_add() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

#endif

/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and subtracts the result from xmm0,xmm1,xmm2
*/

#define _sse_vector_xch_i_sub() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with i and -i
* respectively and adds these registers to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_addsub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn14) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with -i and i
* respectively and adds these registers to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_subadd() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn23) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words in xmm3,xmm4,xmm5
*/

#define _sse_vector_xch() \
__asm__ __volatile__ ("shufps $0x4e, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x4e, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x4e, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u,
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers
* xmm0,xmm1,xmm2 are changed
*/

#if (defined SSE3)

#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm4 \n\t" \
                      "movss %2, %%xmm5 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c22.re) \
                      : \
                      "xmm0", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addsubps %%xmm6, %%xmm5 \n\t" \
                      "addsubps %%xmm7, %%xmm3" \
                      : \
                      : \
                      "m" ((u).c31.im), \
                      "m" ((u).c11.im) \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addsubps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm1", "xmm4", "xmm5", \
                      "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addsubps %%xmm0, %%xmm3 \n\t" \
                      "addsubps %%xmm6, %%xmm4 \n\t" \
                      "addsubps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c32.im) \
                      : \
                      "xmm0", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2" \
                      : \
                      : \
                      "m" ((u).c13.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm0", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm1 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm1 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addsubps %%xmm1, %%xmm3 \n\t" \
                      "addsubps %%xmm6, %%xmm4 \n\t" \
                      "addsubps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c13.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm1", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7")

#else

#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm4 \n\t" \
                      "movss %2, %%xmm5 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c22.re) \
                      : \
                      "xmm0", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "mulps %0, %%xmm0 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3" \
                      : \
                      : \
                      "m" (_sse_sgn13), \
                      "m" ((u).c31.im), \
                      "m" ((u).c11.im) \
                      : \
                      "xmm0", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm1", "xmm4", "xmm5", \
                      "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "mulps %3, %%xmm1 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c32.im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm0", "xmm1", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c13.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm0", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm1 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "mulps %3, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm1 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c13.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c33.im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm1", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7")

#endif

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u^dagger,
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm4 \n\t" \
                      "movss %2, %%xmm5 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c22.re) \
                      : \
                      "xmm0", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "mulps %0, %%xmm0 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "subps %%xmm6, %%xmm5 \n\t" \
                      "subps %%xmm7, %%xmm3" \
                      : \
                      : \
                      "m" (_sse_sgn13), \
                      "m" ((u).c13.im), \
                      "m" ((u).c11.im) \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "subps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.re) \
                      : \
                      "xmm1", "xmm4", "xmm5", \
                      "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "mulps %3, %%xmm1 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "subps %%xmm0, %%xmm3 \n\t" \
                      "subps %%xmm6, %%xmm4 \n\t" \
                      "subps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c21.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c23.im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm0", "xmm1", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c31.re), \
                      "m" ((u).c32.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm0", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7"); \
__asm__ __volatile__ ("movss %0, %%xmm1 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "mulps %3, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm1 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "subps %%xmm1, %%xmm3 \n\t" \
                      "subps %%xmm6, %%xmm4 \n\t" \
                      "subps %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c31.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c33.im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm1", "xmm2", "xmm3", "xmm4", \
                      "xmm5", "xmm6", "xmm7")

/******************************************************************************
*
*  Macros for the linalg routines
*
******************************************************************************/

/*
*  Loads (z.re,z.re,z.re,z.re) to xmm6 and (-z.im,z.im,-z.im,z.im) to xmm7
*/

#define _sse_load_cmplx(z) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %2, %%xmm7" \
                      : \
                      : \
                      "m" ((z).re), \
                      "m" ((z).im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm6", "xmm7")

/*
*  Computes z*s where s is an sse_vector and z a complex number. It is
*  assumed that s is in xmm0,xmm1,xmm2 and that z has been loaded by
*  _sse_load_cmplx(z). The result appears in xmm3,xmm4,xmm5
*/

#define _sse_mulc() \
__asm__ __volatile__ ("movaps %%xmm0, %%xmm3 \n\t" \
                      "movaps %%xmm1, %%xmm4 \n\t" \
                      "movaps %%xmm2, %%xmm5 \n\t" \
                      "mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm6, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm7, %%xmm3 \n\t" \
                      "mulps %%xmm7, %%xmm4 \n\t" \
                      "mulps %%xmm7, %%xmm5 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm1, %%xmm4 \n\t" \
                      "addps %%xmm2, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Computes r+z.re*s where r and s are sse_vectors and z a complex number.
*  It is assumed that r is in xmm3,xmm4,xmm5, s in xmm0,xmm1,xmm2 and that
*  z has been loaded by _sse_load_cmplx(z)
*/

#define _sse_add_mulc_re() \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm6, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm1, %%xmm4 \n\t" \
                      "addps %%xmm2, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Computes r-z.re*s where r and s are sse_vectors and z a complex number.
*  It is assumed that r is in xmm3,xmm4,xmm5, s in xmm0,xmm1,xmm2 and that
*  z has been loaded by _sse_load_cmplx(z)
*/

#define _sse_sub_mulc_re() \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm6, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2 \n\t" \
                      "subps %%xmm0, %%xmm3 \n\t" \
                      "subps %%xmm1, %%xmm4 \n\t" \
                      "subps %%xmm2, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Computes r+i*z.im*s where r and s are sse_vectors and z a complex number.
*  It is assumed that r is in xmm3,xmm4,xmm5, s in xmm0,xmm1,xmm2 and that
*  z has been loaded by _sse_load_cmplx(z)
*/

#define _sse_add_mulc_im() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "mulps %%xmm7, %%xmm0 \n\t" \
                      "mulps %%xmm7, %%xmm1 \n\t" \
                      "mulps %%xmm7, %%xmm2 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm1, %%xmm4 \n\t" \
                      "addps %%xmm2, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Computes r-i*z.im*s where r and s are sse_vectors and z a complex number.
*  It is assumed that r is in xmm3,xmm4,xmm5, s in xmm0,xmm1,xmm2 and that
*  z has been loaded by _sse_load_cmplx(z)
*/

#define _sse_sub_mulc_im() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "mulps %%xmm7, %%xmm0 \n\t" \
                      "mulps %%xmm7, %%xmm1 \n\t" \
                      "mulps %%xmm7, %%xmm2 \n\t" \
                      "subps %%xmm0, %%xmm3 \n\t" \
                      "subps %%xmm1, %%xmm4 \n\t" \
                      "subps %%xmm2, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Loads the spinor s to the registers xmm0,..,xmm5 in linear order
*/

#define _sse_spinor_load(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %1, %%xmm1 \n\t" \
                      "movaps %2, %%xmm2 \n\t" \
                      "movaps %3, %%xmm3 \n\t" \
                      "movaps %4, %%xmm4 \n\t" \
                      "movaps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1.c1.re), \
                      "m" ((s).c1.c3.re), \
                      "m" ((s).c2.c2.re), \
                      "m" ((s).c3.c1.re), \
                      "m" ((s).c3.c3.re), \
                      "m" ((s).c4.c2.re) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
*  Stores the registers xmm0,..,xmm5 to the spinor s in linear order
*/

#define _sse_spinor_store(s) \
__asm__ __volatile__ ("movaps %%xmm0, %0 \n\t" \
                      "movaps %%xmm1, %1 \n\t" \
                      "movaps %%xmm2, %2 \n\t" \
                      "movaps %%xmm3, %3 \n\t" \
                      "movaps %%xmm4, %4 \n\t" \
                      "movaps %%xmm5, %5" \
                      : \
                      "=m" ((s).c1.c1.re), \
                      "=m" ((s).c1.c3.re), \
                      "=m" ((s).c2.c2.re), \
                      "=m" ((s).c3.c1.re), \
                      "=m" ((s).c3.c3.re), \
                      "=m" ((s).c4.c2.re))

#endif

