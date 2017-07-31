
/*******************************************************************************
*
* File swinit.c
*
* Copyright (C) 2006 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Utility programs for the global and the block Pauli fields
*
* The externally accessible functions are
*
*   void alloc_sw(void)
*     Allocates the single-precision Pauli field sw and initializes
*     the matrices to unity
*
*   void alloc_swd(void)
*     Allocates the double-precision Pauli field swd and initializes
*     the matrices to unity
*
*   void free_sw(void)
*     Frees the single-precision Pauli field sw
*
*   void free_swd(void)
*     Frees the double-precision Pauli field swd
*
*   void assign_swd2sw(void)
*     Assigns the double-precision Pauli field swd to the single-precision
*     field sw
*
*   int invert_swd(ptset_t set)
*     Inverts the Pauli matrices in the field swd[] on the specified point
*     set. The program returns 0 if all matrices could be safely inverted
*     and 1 otherwise (see pauli_dble.c)
*     
* Notes:
*
* All these programs act globally and must be called simultaneously on all
* processes. The value returned by invert_swd() is the same on all processes
*
*******************************************************************************/

#define SWINIT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "sw_term.h"
#include "start.h"
#include "global.h"

#if (defined SSE)
#include "sse.h"
#endif

static const pauli sw0={{0.0f}};
static const pauli_dble swd0={{0.0}};


void alloc_sw(void)
{
   pauli *pa,*pam,unity;

   if (sw!=NULL)
      return;

   error_root(sizeof(pauli)!=(36*sizeof(float)),1,"alloc_sw [swinit.c]",
              "The pauli structures are not properly packed");
   
   sw=amalloc(2*VOLUME*sizeof(pauli),ALIGN);
   error(sw==NULL,1,"alloc_sw [swinit.c]",
         "Unable to allocate the global sw array");

   unity=sw0;
   unity.u[0]=1.0f;
   unity.u[1]=1.0f;
   unity.u[2]=1.0f;
   unity.u[3]=1.0f;
   unity.u[4]=1.0f;
   unity.u[5]=1.0f;

   pa=sw;
   pam=pa+2*VOLUME;

   for (;pa<pam;pa++)
      *pa=unity;
}


void alloc_swd(void)
{
   pauli_dble *pa,*pam,unity;

   if (swd!=NULL)
      return;

   error_root(sizeof(pauli_dble)!=(36*sizeof(double)),1,"alloc_swd [swinit.c]",
              "The pauli_dble structures are not properly packed");
   
   swd=amalloc(2*VOLUME*sizeof(pauli_dble),ALIGN);
   error(swd==NULL,1,"alloc_swd [swinit.c]",
         "Unable to allocate the global swd array");

   unity=swd0;
   unity.u[0]=1.0;
   unity.u[1]=1.0;
   unity.u[2]=1.0;
   unity.u[3]=1.0;
   unity.u[4]=1.0;
   unity.u[5]=1.0;

   pa=swd;
   pam=pa+2*VOLUME;

   for (;pa<pam;pa++)
      *pa=unity;
}


void free_sw(void)
{
   if (sw!=NULL)
   {
      afree(sw);
      sw=NULL;
   }
}


void free_swd(void)
{
   if (swd!=NULL)
   {
      afree(swd);
      swd=NULL;
   }
}


void assign_swd2sw(void)
{
   error((sw==NULL)||(swd==NULL),1,"assign_swd2sw [swinit.c]",
          "sw or swd is not allocated");

   assign_pauli(2*VOLUME,swd,sw);
}


static int iswd(int vol,int icom,pauli_dble *pa)
{
   int n,ifail;
   pauli_dble *pam;

   ifail=0;
   pam=pa+vol;

   for (;pa<pam;pa++)
   {
#if (defined SSE)
      pa+=2;
      _prefetch_pauli_dble(pa);
      pa-=2;
#endif      
      ifail|=inv_pauli_dble(pa,pa);
   }

   if (icom==1)
   {
      n=ifail;
   
      MPI_Reduce(&n,&ifail,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); 
      MPI_Bcast(&ifail,1,MPI_INT,0,MPI_COMM_WORLD);
   }
   
   return ifail;
}


int invert_swd(ptset_t set)
{
   int vol,iprms[1];
   pauli_dble *pa;

   iprms[0]=(int)(set);
   MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

   error((iprms[0]!=(int)(set))||(swd==NULL),1,"invert_swd [swinit.c]",
         "Parameter is not global or swd is not allocated");

   if (set==EVEN_PTS)
   {
      vol=VOLUME;
      pa=swd;
   }
   else if (set==ODD_PTS)
   {
      vol=VOLUME;
      pa=swd+VOLUME;
   }
   else if (set==ALL_PTS)
   {
      vol=2*VOLUME;
      pa=swd;
   }
   else
      return 0;

   return iswd(vol,1,pa);
}

