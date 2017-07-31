
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of prod2su3alg and rotate_su3alg
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "start.h"
#include "sw_term.h"
#include "misc.h"

static const su3_vector_dble vd0={{0.0}};
static const spinor_dble sd0={{{0.0}}};

#if (defined SSE)
static su3_dble Q,u,v __attribute__ ((aligned (16)));
static su3_alg_dble X __attribute__ ((aligned (16)));
#else
static su3_dble Q,u,v;
static su3_alg_dble X;
#endif


static void X2u(su3_alg_dble *X,su3_dble *u)
{
   (*u).c11.re=0.0;
   (*u).c11.im= (*X).c1+(*X).c2;
   (*u).c22.re=0.0;
   (*u).c22.im= (*X).c2-2.0*(*X).c1;
   (*u).c33.re=0.0;
   (*u).c33.im= (*X).c1-2.0*(*X).c2;

   (*u).c12.re= (*X).c3;
   (*u).c12.im= (*X).c4;   
   (*u).c21.re=-(*X).c3;
   (*u).c21.im= (*X).c4;   

   (*u).c13.re= (*X).c5;
   (*u).c13.im= (*X).c6;   
   (*u).c31.re=-(*X).c5;
   (*u).c31.im= (*X).c6;   
   
   (*u).c23.re= (*X).c7;
   (*u).c23.im= (*X).c8;   
   (*u).c32.re=-(*X).c7;
   (*u).c32.im= (*X).c8;   
}


int main(void)
{
   int i;
   double d,*r1,*r2;
   
   printf("\n");
   printf("Check of prod2su3alg and rotate_su3alg\n");
   printf("--------------------------------------\n\n");

   rlxd_init(1,23456);

   printf("prod2su3alg:\n");
   ranlxd(&u.c11.re,18);   
   ranlxd(&v.c11.re,18);
   ranlxd(&X.c1,8);

   prod2su3alg(&u,&v,&X);   
   _su3_times_su3(Q,u,v);
   _su3_dagger(u,Q);

   r1=&Q.c11.re;
   r2=&u.c11.re;
   
   for (i=0;i<18;i++)
      r1[i]=0.5*(r1[i]-r2[i]);

   d=(Q.c11.im+Q.c22.im+Q.c33.im)/3.0;
   Q.c11.im-=d;
   Q.c22.im-=d;
   Q.c33.im-=d;
   
   d=fabs(Q.c11.im-X.c1-X.c2);
   printf("X.c11.im: %.2e\n",d);
   d=fabs(Q.c22.im+2.0*X.c1-X.c2);
   printf("X.c22.im: %.2e\n",d);
   d=fabs(Q.c33.im-X.c1+2.0*X.c2);
   printf("X.c33.im: %.2e\n",d);

   d=fabs(Q.c12.re-X.c3);
   printf("X.c12.re: %.2e\n",d);
   d=fabs(Q.c12.im-X.c4);
   printf("X.c12.im: %.2e\n",d);

   d=fabs(Q.c13.re-X.c5);
   printf("X.c13.re: %.2e\n",d);
   d=fabs(Q.c13.im-X.c6);
   printf("X.c13.im: %.2e\n",d);

   d=fabs(Q.c23.re-X.c7);
   printf("X.c23.re: %.2e\n",d);
   d=fabs(Q.c23.im-X.c8);
   printf("X.c23.im: %.2e\n",d);

   printf("\nrotate_su3alg:\n");
   random_su3_dble(&u);
   ranlxd(&X.c1,8);
   X2u(&X,&v);

   rotate_su3alg(&u,&X);   

   _su3_times_su3(Q,u,v);
   _su3_dagger(v,u);
   _su3_times_su3(u,Q,v);
   
   d=fabs(u.c11.im-X.c1-X.c2);
   printf("X.c11.im: %.2e\n",d);
   d=fabs(u.c22.im+2.0*X.c1-X.c2);
   printf("X.c22.im: %.2e\n",d);
   d=fabs(u.c33.im-X.c1+2.0*X.c2);
   printf("X.c33.im: %.2e\n",d);

   d=fabs(u.c12.re-X.c3);
   printf("X.c12.re: %.2e\n",d);
   d=fabs(u.c12.im-X.c4);
   printf("X.c12.im: %.2e\n",d);

   d=fabs(u.c13.re-X.c5);
   printf("X.c13.re: %.2e\n",d);
   d=fabs(u.c13.im-X.c6);
   printf("X.c13.im: %.2e\n",d);

   d=fabs(u.c23.re-X.c7);
   printf("X.c23.re: %.2e\n",d);
   d=fabs(u.c23.im-X.c8);
   printf("X.c23.im: %.2e\n\n",d);
   
   exit(0);
}
