/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "mpi.h"


int safe_mod(int x,int y)
{
   if (x>=0)
      return x%y;
   else
      return (y-(abs(x)%y))%y;
}


void error(int test,int no,char *name,char *format,...)
{
   int i,all,my_rank;
   va_list args;

   i=(test!=0);
   MPI_Allreduce(&i,&all,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

   if (all==0)
      return;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      fflush(stdout);

      MPI_Abort(MPI_COMM_WORLD,no);
   }
   else
      for (i=1;i<2;i=safe_mod(i,2));
}


void error_root(int test,int no,char *name,char *format,...)
{
   int my_rank;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if ((my_rank==0)&&(test!=0))
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      fflush(stdout);

      MPI_Abort(MPI_COMM_WORLD,no);
   }
}


char *alloc_buffer(size_t n)
{
   int shift;
   unsigned long mask,addr;
   char *base;

   shift=16;
   mask=(unsigned long)(shift-1);
   base=malloc(n+shift);

   error(base==NULL,1,"alloc_buffer [utils.c]",
         "Unable to allocate buffer");
   
   addr=((unsigned long)(base)+shift)&(~mask);

   return (char*)(addr);
}


void init_buffer(int n,double *buf)
{
   int i,r,m;

   r=1;
   m=8125;

   for (i=0;i<n;i++)
   {
      buf[i]=(double)(r);
      r=(m*r)%1021;
   }
}
