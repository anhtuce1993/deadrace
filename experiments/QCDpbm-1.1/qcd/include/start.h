
/*******************************************************************************
*
* File start.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef START_H
#define START_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef enum
{
   ALL_PTS,EVEN_PTS,ODD_PTS,NO_PTS,PT_SETS
} ptset_t;

typedef struct
{
   double beta,kappa,m0,csw;
} lat_parms_t;

#ifndef GEOMETRY_C
extern int ipr_global(int n[]);
extern void ipt_global(int x[],int *ip,int *ix);
extern void geometry(void);
#endif

#ifndef PARMS_C
extern lat_parms_t set_lat_parms(double beta,double kappa,double csw);
extern lat_parms_t lat_parms(void);
#endif

#ifndef START_C
extern void start_ranlux(int level,int seed);
extern void alloc_u(void);
extern void alloc_ud(void);
extern void alloc_s(int no_fields);
extern void alloc_sd(int no_fields);
extern void free_u(void);
extern void free_ud(void);
extern void free_s(void);
extern void free_sd(void);
extern void random_u(void);
extern void random_ud(void);
extern void renormalize_u(void);
extern void renormalize_ud(void);
extern void assign_u2ud(void);
extern void assign_ud2u(void);
#endif

#ifndef SINIT_C
extern void set_s2zero(int vol,spinor *pk);
extern void set_sd2zero(int vol,spinor_dble *pk);
extern void random_s(int vol,spinor *pk,float sigma);
extern void random_sd(int vol,spinor_dble *pk,double sigma);
extern void assign_s2s(int vol,spinor *pk,spinor *pl);
extern void assign_s2sd(int vol,spinor *pk,spinor_dble *pl);
extern void assign_sd2s(int vol,spinor_dble *pk,spinor *pl);
extern void assign_sd2sd(int vol,spinor_dble *pk,spinor_dble *pl);
#endif

#ifndef UCOM_C
extern void copy_bnd_u(void);
extern void copy_bnd_ud(void);
extern void plaq_u(int ix,int mu,int nu,su3 **u);
extern void plaq_ud(int ix,int mu,int nu,su3_dble **u);
extern void free_ucom_bufs(int iu,int iud);
#endif

#ifndef UTILS_C
extern int safe_mod(int x,int y);
extern int mpi_permanent_tag(void);
extern int mpi_tag(void);
extern void *amalloc(size_t size,int p);
extern void afree(void *addr);
extern void error(int test,int no,char *name,char *format,...);
extern void error_root(int test,int no,char *name,char *format,...);
extern int error_loc(int test,int no,char *name,char *text);
extern void error_chk(void);
#endif

#endif
