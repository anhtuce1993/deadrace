
/*******************************************************************************
*
* File nompi/start.h
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
   ALL_PTS,EVEN_PTS,ODD_PTS,NO_PTS
} ptset_t;

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

#ifndef UTILS_C
extern int safe_mod(int x,int y);
extern void *amalloc(size_t size,int p);
extern void afree(void *addr);
extern void error(int test,int no,char *name,char *format,...);
extern void error_root(int test,int no,char *name,char *format,...);
extern int error_loc(int test,int no,char *name,char *text);
#endif

#ifndef MUTILS_C
extern int find_opt(int argc,char *argv[],char *opt);
extern int digits(double x,double dx,char *fmt);
#endif

#endif
