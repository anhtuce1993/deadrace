
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of apply_sw and apply_sw_dble
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "start.h"
#include "sw_term.h"

#define NM 10000


int main(void)
{
   int k,itest;
   double d,dmax;
   complex *cs;
   complex_dble *csd,*crd,z;
   pauli *msb;
   pauli_dble *md,*mdb;
   spinor *s,*sb;
   spinor_dble *sd,*sdb;
   
   printf("\n");
   printf("Check of apply_sw and apply_sw_dble\n");
   printf("-----------------------------------\n\n");
   
   msb=amalloc(2*NM*sizeof(pauli),4);
   mdb=amalloc(4*NM*sizeof(pauli_dble),4);
   sb=amalloc(2*NM*sizeof(spinor),4);
   sdb=amalloc(2*NM*sizeof(spinor_dble),4);
   
   error((msb==NULL)||(mdb==NULL)||(sb==NULL)||(sdb==NULL),1,
         "main [check6.c]","Unable to allocate data arrays");

   rlxd_init(1,1234);
   random_sd(NM,sdb,1.0);
   assign_sd2s(NM,sdb,sb);

   itest=0;
   
   for (md=mdb;md<(mdb+2*NM);md++)
   {
      gauss_dble((*md).u,36);

      for (k=0;k<36;k++)
      {
         (*md).u[k]*=0.05;

         if (k<6)
            (*md).u[k]+=1.0;
      }
         
      itest+=inv_pauli_dble(md,md+2*NM);
   }
   
   assign_pauli(2*NM,mdb,msb);
   apply_sw(NM,msb,sb,sb+NM);
   apply_sw_dble(NM,mdb,sdb,sdb+NM);
   
   dmax=0.0;
   sd=sdb;
   
   for (s=sb;s<(sb+2*NM);s++)
   {
      cs=(complex*)(s);
      csd=(complex_dble*)(sd);

      for (k=0;k<12;k++)
      {
         z.re=csd[k].re-(double)(cs[k].re);
         z.im=csd[k].im-(double)(cs[k].im);
         d=z.re*z.re+z.im*z.im;
         if (d>dmax)
            dmax=d;
      }
      
      sd+=1;
   }
   
   printf("Comparison of apply_sw with apply_swd: %.1e\n\n",sqrt(dmax));

   apply_sw_dble(NM,mdb+2*NM,sdb+NM,sdb+NM);

   dmax=0.0;
   
   for (sd=sdb;sd<(sdb+NM);sd++)
   {
      csd=(complex_dble*)(sd);
      crd=(complex_dble*)(sd+NM);

      for (k=0;k<12;k++)
      {
         z.re=csd[k].re-crd[k].re;
         z.im=csd[k].im-crd[k].im;
         d=z.re*z.re+z.im*z.im;
         if (d>dmax)
            dmax=d;
      }
   }

   printf("Invert all Pauli matrices (no of failures = %d)\n\n",itest);
   printf("Maximal deviation |M^(-1)*M*sd-sd| = %.1e\n\n",sqrt(dmax));
   exit(0);
}
