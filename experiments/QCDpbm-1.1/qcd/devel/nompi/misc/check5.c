
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of su3xsu3, su3dagxsu3, ...
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


static double max_dev(su3_dble *u,su3_dble *v)
{
   complex_dble *r,*s,*rm;
   double nrm,d,dmax;

   r=(complex_dble*)(&(*u).c11);
   s=(complex_dble*)(&(*v).c11);
   rm=r+9;

   nrm=0.0;
   dmax=0.0;
   
   for (;r<rm;r++)
   {
      nrm+=(*r).re*(*r).re+(*r).im*(*r).im;

      d=((*r).re-(*s).re)*((*r).re-(*s).re)+
        ((*r).im-(*s).im)*((*r).im-(*s).im);
      if (d>dmax)
         dmax=d;

      s+=1;
   }

   return sqrt(dmax/nrm);
}


static void X2u(u3_alg_dble *X,su3_dble *u)
{
   (*u).c11.re=0.0;
   (*u).c11.im= (*X).c1;
   (*u).c22.re=0.0;
   (*u).c22.im= (*X).c2;
   (*u).c33.re=0.0;
   (*u).c33.im= (*X).c3;

   (*u).c12.re= (*X).c4;
   (*u).c12.im= (*X).c5;   
   (*u).c21.re=-(*X).c4;
   (*u).c21.im= (*X).c5;   

   (*u).c13.re= (*X).c6;
   (*u).c13.im= (*X).c7;   
   (*u).c31.re=-(*X).c6;
   (*u).c31.im= (*X).c7;   
   
   (*u).c23.re= (*X).c8;
   (*u).c23.im= (*X).c9;   
   (*u).c32.re=-(*X).c8;
   (*u).c32.im= (*X).c9;   
}


int main(void)
{
   double d1,d2,d3,d4;
   su3_dble *u,*v,*w1,*w2;
   u3_alg_dble *X;
   
   printf("\n");
   printf("Check of su3xsu3, su3dagxsu3, ...\n");
   printf("---------------------------------\n\n");

   u=amalloc(4*sizeof(su3_dble),4);
   X=amalloc(sizeof(u3_alg_dble),3);
   error((u==NULL)||(X==NULL),1,"main [check5.c]",
         "Unable to allocate auxiliary arrays");

   v=u+1;
   w1=u+2;
   w2=u+3;
   
   rlxd_init(1,23456);
   
   random_su3_dble(u);
   random_su3_dble(v);
   su3xsu3(u,v,w1);
   _su3_times_su3(*w2,*u,*v);
   d1=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3dagxsu3(u,v,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   _su3_times_su3(*w2,*u,*v);
   d2=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3xsu3dag(u,v,w1);
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d3=max_dev(w1,w2);

   random_su3_dble(u);
   random_su3_dble(v);
   su3dagxsu3dag(u,v,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d4=max_dev(w1,w2);

   printf("su3xsu3:       %.2e\n",d1);   
   printf("su3dagxsu3:    %.2e\n",d2);
   printf("su3xsu3dag:    %.2e\n",d3);
   printf("su3dagxsu3dag: %.2e\n",d4);      

   random_su3_dble(u);
   ranlxd(&(*X).c1,9);   
   su3xu3alg(u,X,w1);
   X2u(X,v);
   _su3_times_su3(*w2,*u,*v);
   d1=max_dev(w1,w2);

   random_su3_dble(u);
   ranlxd(&(*X).c1,9);   
   su3dagxu3alg(u,X,w1);
   _su3_dagger(*w2,*u);
   *u=*w2;
   X2u(X,v);
   _su3_times_su3(*w2,*u,*v);
   d2=max_dev(w1,w2);

   random_su3_dble(v);
   ranlxd(&(*X).c1,9);   
   u3algxsu3(X,v,w1);
   X2u(X,u);
   _su3_times_su3(*w2,*u,*v);
   d3=max_dev(w1,w2);

   random_su3_dble(v);
   ranlxd(&(*X).c1,9);   
   u3algxsu3dag(X,v,w1);
   X2u(X,u);
   _su3_dagger(*w2,*v);
   *v=*w2;
   _su3_times_su3(*w2,*u,*v);
   d4=max_dev(w1,w2);
   
   printf("su3xu3alg:     %.2e\n",d1);   
   printf("su3dagxu3alg:  %.2e\n",d2);
   printf("u3algxsu3:     %.2e\n",d3);
   printf("u3algxsu3dag:  %.2e\n\n",d4);      

   exit(0);
}
