
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Accuracy of inv_pauli_dble
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
   int n,k,l,*is,itot;
   double fact,d,dmax;
   complex *cvs;
   complex_dble *cvd;
   weyl *vs,vs0={{{0.0f}}};
   weyl_dble *vd,vd0={{{0.0}}};
   pauli *ms,*ims,*msb,*imsb;
   pauli_dble *md,*imd,*mdb,*imdb;

   printf("\n");
   printf("Accuracy of inv_pauli_dble\n");
   printf("--------------------------\n\n");

   is=amalloc(NM*sizeof(int),2);
   vs=amalloc(sizeof(weyl),4);
   vd=amalloc(sizeof(weyl),4);
   msb=amalloc(2*NM*sizeof(pauli),4);
   mdb=amalloc(2*NM*sizeof(pauli_dble),4);
   error((is==NULL)||(vs==NULL)||(vd==NULL)||(msb==NULL)||(mdb==NULL),1,
         "main [check5.c]","Unable to allocate auxiliary arrays");

   cvs=(complex*)(vs);
   cvd=(complex_dble*)(vd);
   imsb=msb+NM;
   imdb=mdb+NM;
   
   rlxd_init(1,1234);
   md=mdb;
   imd=imdb;
   itot=0;
   dmax=0.0;
   fact=sqrt(2.0);

   for (n=0;n<NM;n++)
   {
      gauss_dble((*md).u,36);

      for (k=0;k<6;k++)
         (*md).u[k]*=fact;
         
      is[n]=inv_pauli_dble(md,imd);

      if (is[n]==0)
      {
         for (k=0;k<6;k++)
         {
            *vd=vd0;
            cvd[k].re=1.0;

            mul_pauli_dble(md,vd,vd);
            mul_pauli_dble(imd,vd,vd);
            cvd[k].re-=1.0;

            for (l=0;l<6;l++)
            {
               d=cvd[l].re*cvd[l].re+cvd[l].im*cvd[l].im;
               if (d>dmax)
                  dmax=d;
            }
         }
      }
      else
         itot+=1;
      
      md+=1;
      imd+=1;
   }

   printf("Double-precision program:\n");
   printf("%d gaussian random matrices, %d inversion failures\n",NM,itot);
   printf("Maximal deviation |{M*M^(-1)-1}_ij| = %.1e ",sqrt(dmax));
   printf("(safe cases only)\n\n");

   assign_pauli(NM,mdb,msb);
   assign_pauli(NM,imdb,imsb);

   ms=msb;
   ims=imsb;
   dmax=0.0;

   for (n=0;n<NM;n++)
   {
      if (is[n]==0)
      {
         for (k=0;k<6;k++)
         {
            *vs=vs0;
            cvs[k].re=1.0f;

            mul_pauli(ms,vs,vs);
            mul_pauli(ims,vs,vs);
            cvs[k].re-=1.0f;

            for (l=0;l<6;l++)
            {
               d=(double)(cvs[l].re*cvs[l].re+cvs[l].im*cvs[l].im);
               if (d>dmax)
                  dmax=d;
            }
         }
      }
      
      ms+=1;
      ims+=1;
   }

   printf("After assignment to single-precision matrices:\n");
   printf("Maximal deviation |{M*M^(-1)-1}_ij| = %.1e ",sqrt(dmax));
   printf("(safe cases only)\n\n");
   exit(0);
}
