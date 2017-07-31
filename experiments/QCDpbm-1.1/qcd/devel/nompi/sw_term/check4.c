
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Quick test of inv_pauli_dble
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "start.h"
#include "sw_term.h"


int main(void)
{
   int i;
   double d,dmax;
   complex_dble *r1,*r2;
   weyl_dble *v,*w;
   pauli_dble *m,*im;

   printf("\n");
   printf("Quick test of inv_pauli_dble\n");
   printf("----------------------------\n\n");

   v=amalloc(2*sizeof(weyl_dble),4);
   m=amalloc(2*sizeof(pauli_dble),4);
   error((v==NULL)||(m==NULL),1,"main [check4.c]",
         "Unable to allocate auxiliary arrays");

   w=v+1;
   im=m+1;
   
   rlxd_init(1,123456);
   ranlxd((double*)(v),12);
   ranlxd((*m).u,36);
   (*m).u[0]-=2.0;
   (*m).u[1]-=2.0;
   (*m).u[2]-=2.0;
   (*m).u[3]-=2.0;
   (*m).u[4]-=2.0;
   (*m).u[5]-=2.0;
   
   i=inv_pauli_dble(m,im);
   mul_pauli_dble(m,v,w);
   mul_pauli_dble(im,w,w);

   printf("Exit code of inv_pauli_dble = %d\n\n",i);

   if (i==0)
   {
      printf("r1: expected result, r2: actual result\n\n");

      r1=(complex_dble*)(v);
      r2=(complex_dble*)(w);
   
      for (i=0;i<6;i++)
      {
         printf("r1[%d]=(%.16e,%.16e)\n",i,r1[i].re,r1[i].im);
         printf("r2[%d]=(%.16e,%.16e)\n\n",i,r2[i].re,r2[i].im);
      }
   }

   i=inv_pauli_dble(m,m);

   dmax=0.0;

   for (i=0;i<36;i++)
   {
      d=fabs((*m).u[i]-(*im).u[i]);
      if (d>dmax)
         dmax=d;
   }

   printf("Inversion in place: dmax = %.1e (should be 0.0)\n\n",dmax);
   exit(0);
}
