
/*******************************************************************************
*
* File random.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef RANLXS_C
extern void ranlxs(float r[],int n);
extern void rlxs_init(int level,int seed);
extern int rlxs_size(void);
extern void rlxs_get(int state[]);
extern void rlxs_reset(int state[]);
#endif

#ifndef RANLXD_C
extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);
#endif

#ifndef GAUSS_C
extern void gauss(float r[],int n);
extern void gauss_dble(double r[],int n);
#endif

#ifndef RANDOM_SU3_C
extern void random_su3(su3 *u);
extern void random_su3_dble(su3_dble *u);
#endif

#endif
