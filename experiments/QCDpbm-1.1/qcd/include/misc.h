
/*******************************************************************************
*
* File misc.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef MISC_H
#define MISC_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef CMATRIX_C
extern void cmat_vec(int n,complex *a,complex *v,complex *w);
extern void cmat_vec_assign(int n,complex *a,complex *v,complex *w);
extern void cmat_add(int n,complex *a,complex *b,complex *c);
extern void cmat_sub(int n,complex *a,complex *b,complex *c);
extern void cmat_mul(int n,complex *a,complex *b,complex *c);
extern void cmat_dag(int n,complex *a,complex *b);
#endif

#ifndef CMATRIX_DBLE_C
extern void cmat_vec_dble(int n,complex_dble *a,complex_dble *v,
                          complex_dble *w);
extern void cmat_vec_assign_dble(int n,complex_dble *a,complex_dble *v,
                                 complex_dble *w);
extern void cmat_add_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_sub_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_mul_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_dag_dble(int n,complex_dble *a,complex_dble *b);
extern int cmat_inv_dble(int n,complex_dble *a,complex_dble *b,double *k);
#endif

#ifndef SU3_FCTS_C
extern void cross_prod(su3_vector *v1,su3_vector *v2,su3_vector *v3);
extern void cross_prod_dble(su3_vector_dble *v1,su3_vector_dble *v2,
                            su3_vector_dble *v3);
extern void project_to_su3(su3 *u);
extern void project_to_su3_dble(su3_dble *u);
#endif

#ifndef SU3_PRODS_C
extern void su3xsu3(su3_dble *u,su3_dble *v,su3_dble *w);
extern void su3dagxsu3(su3_dble *u,su3_dble *v,su3_dble *w);
extern void su3xsu3dag(su3_dble *u,su3_dble *v,su3_dble *w);
extern void su3dagxsu3dag(su3_dble *u,su3_dble *v,su3_dble *w);
extern void su3xu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v);
extern void su3dagxu3alg(su3_dble *u,u3_alg_dble *X,su3_dble *v);
extern void u3algxsu3(u3_alg_dble *X,su3_dble *u,su3_dble *v);
extern void u3algxsu3dag(u3_alg_dble *X,su3_dble *u,su3_dble *v);
extern void prod2su3alg(su3_dble *u,su3_dble *v,su3_alg_dble *X);
extern void add_prod2u3alg(su3_dble *u,su3_dble *v,u3_alg_dble *X);
extern void rotate_su3alg(su3_dble *u,su3_alg_dble *X);
#endif

#endif
