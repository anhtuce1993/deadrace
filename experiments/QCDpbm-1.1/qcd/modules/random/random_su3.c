
/*******************************************************************************
*
* File random_su3.c
*
* Copyright (C) 2004 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generation of uniformly distributed single- and double-precision
* SU(3) matrices
*
* The externally accessible function is
*
*   void random_su3(su3 *u)
*     Generates a random single-precision SU(3) matrix and assigns it to *u
*
*   void random_su3_dble(su3_dble *u)
*     Generates a random double-precision SU(3) matrix and assigns it to *u
*
* Notes:
*
* The random matrices are uniformly distributed over SU(3) to a precision
* given by the number of significant bits of the random numbers returned by
* ranlxs and ranlxd respectively. Rougly speaking one can expect the matrices
* to be uniformly random up to systematic deviations from 1 at the level of 
* 10^(-7) and 10^(-14) in the single- and double-precision programs
*
*******************************************************************************/

#define RANDOM_SU3_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"
#include "random.h"
#include "misc.h"


static void random_su3_vector(su3_vector *v)
{
   int i;
   float *r,norm,fact;

   r=(float*)(v);
   norm=0.0f;

   while (norm<=FLT_EPSILON)
   {
      gauss(r,6);
      norm=0.0f;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=(float)sqrt((double)norm);
   }

   fact=1.0f/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;
}


void random_su3(su3 *u)
{
   int i;
   float *r,norm,fact;
   su3_vector *v1,*v2,*v3;

   v1=(su3_vector*)(u);
   v2=v1+1;
   v3=v1+2;
   r=(float*)(v3);

   random_su3_vector(v1);
   norm=0.0f;

   while (norm<=FLT_EPSILON)
   {
      random_su3_vector(v2);
      cross_prod(v1,v2,v3);
      norm=0.0f;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=(float)sqrt((double)norm);
   }        

   fact=1.0f/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;   

   cross_prod(v3,v1,v2);
}


static void random_su3_vector_dble(su3_vector_dble *v)
{
   int i;
   double *r,norm,fact;

   r=(double*)(v);
   norm=0.0;

   while (norm<=DBL_EPSILON)
   {
      gauss_dble(r,6);
      norm=0.0;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=sqrt(norm);
   }

   fact=1.0/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;
}


void random_su3_dble(su3_dble *u)
{
   int i;
   double *r,norm,fact;
   su3_vector_dble *v1,*v2,*v3;

   v1=(su3_vector_dble*)(u);
   v2=v1+1;
   v3=v1+2;
   r=(double*)(v3);

   random_su3_vector_dble(v1);
   norm=0.0;

   while (norm<=DBL_EPSILON)
   {
      random_su3_vector_dble(v2);
      cross_prod_dble(v1,v2,v3);
      norm=0.0;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=sqrt(norm);
   }        

   fact=1.0/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;   

   cross_prod_dble(v3,v1,v2);
}

