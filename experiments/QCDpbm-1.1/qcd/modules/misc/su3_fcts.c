
/*******************************************************************************
*
* File su3_fcts.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Collection of some basic programs for SU(3) vectors and matrices 
*
* The externally accessible function are
*
*   void cross_prod(su3_vector *v1,su3_vector *v2,su3_vector *v3)
*     Assigns the complex conjugate of the cross product (*v1 x *v2)
*     to *v3
*
*   void cross_prod_dble(su3_vector_dble *v1,su3_vector_dble *v2,
*                        su3_vector_dble *v3)
*     Double-precision version of cross_prod
*
*   void project_to_su3(su3 *u)
*     Projects an approximate SU(3) matrix back to SU(3)
*
*   void project_to_su3_dble(su3_dble *u)
*     Double-precision version of project_to_su3
*
*******************************************************************************/

#define SU3_FCTS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"


static void normalize(su3_vector *v)
{
   int i;
   float *r,fact;

   r=(float*)(v);
   fact=0.0f;

   for (i=0;i<6;i++)
      fact+=r[i]*r[i];

   fact=1.0f/(float)sqrt((double)(fact));

   for (i=0;i<6;i++)
      r[i]*=fact;
}


static void normalize_dble(su3_vector_dble *v)
{
   int i;
   double *r,fact;

   r=(double*)(v);
   fact=0.0;

   for (i=0;i<6;i++)
      fact+=r[i]*r[i];

   fact=1.0/sqrt(fact);

   for (i=0;i<6;i++)
      r[i]*=fact;
}


void cross_prod(su3_vector *v1,su3_vector *v2,su3_vector *v3)
{
   _vector_cross_prod(*v3,*v1,*v2);
}


void cross_prod_dble(su3_vector_dble *v1,su3_vector_dble *v2,
                     su3_vector_dble *v3)
{
   _vector_cross_prod(*v3,*v1,*v2);
}


void project_to_su3(su3 *u)
{
   su3_vector *v1,*v2,*v3;
   
   v1=(su3_vector*)(u);
   v2=v1+1;
   v3=v1+2;
   
   normalize(v1);
   cross_prod(v1,v2,v3);
   normalize(v3);
   cross_prod(v3,v1,v2);   
}


void project_to_su3_dble(su3_dble *u)
{
   su3_vector_dble *v1,*v2,*v3;
   
   v1=(su3_vector_dble*)(u);
   v2=v1+1;
   v3=v1+2;
   
   normalize_dble(v1);
   cross_prod_dble(v1,v2,v3);
   normalize_dble(v3);
   cross_prod_dble(v3,v1,v2);   
}

