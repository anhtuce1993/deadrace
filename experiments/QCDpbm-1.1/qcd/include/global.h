
/*******************************************************************************
*
* File global.h
*
* Copyright (C) 2005, 2007 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define NPROC0 2
#define NPROC1 2
#define NPROC2 1
#define NPROC3 1
#define L0 16
#define L1 16
#define L2 16
#define L3 16

#define NAME_SIZE 128

/****************************** do not change *********************************/

#if ((NPROC0<1)||(NPROC1<1)||(NPROC2<1)||(NPROC3<1)|| \
    ((NPROC0>1)&&((NPROC0%2)!=0))||((NPROC1>1)&&((NPROC1%2)!=0))|| \
    ((NPROC2>1)&&((NPROC2%2)!=0))||((NPROC3>1)&&((NPROC3%2)!=0)))
#error : The number of processes in each direction must be 1 or a multiple of 2
#endif

#if ((L0<4)||(L1<4)||(L2<4)||(L3<4)|| \
    ((L0%2)!=0)||((L1%2)!=0)||((L2%2)!=0)||((L3%2)!=0))
#error : The local lattice sizes must be even and not smaller than 4
#endif

#if (NAME_SIZE<128)
#error : NAME_SIZE must be greater or equal to 128
#endif

#define NPROC (NPROC0*NPROC1*NPROC2*NPROC3)
#define VOLUME (L0*L1*L2*L3)
#define FACE0 ((1-(NPROC0%2))*L1*L2*L3)
#define FACE1 ((1-(NPROC1%2))*L2*L3*L0)
#define FACE2 ((1-(NPROC2%2))*L3*L0*L1)
#define FACE3 ((1-(NPROC3%2))*L0*L1*L2)
#define BNDRY (2*(FACE0+FACE1+FACE2+FACE3))
#define NSPIN (VOLUME+(BNDRY/2))

#if ((defined P4)||(defined DH))
#define ALIGN 7
#elif (defined PM)
#define ALIGN 6
#else
#define ALIGN 5
#endif

#ifndef SU3_H
#include "su3.h"
#endif

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN int cpr[4];
EXTERN int npr[8];

EXTERN int ipt[VOLUME];
EXTERN int iup[VOLUME][4];
EXTERN int idn[VOLUME][4];
EXTERN int map[BNDRY+NPROC%2];
EXTERN int no_s,no_sd;

EXTERN su3 *pu[VOLUME][4];
EXTERN su3_dble *pud[VOLUME][4];

EXTERN spinor *(*ps)[NSPIN];
EXTERN spinor_dble *(*psd)[NSPIN];

EXTERN pauli *sw;
EXTERN pauli_dble *swd;

#undef EXTERN

#endif
