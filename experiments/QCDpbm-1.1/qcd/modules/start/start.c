
/*******************************************************************************
*
* File start.c
*
* Copyright (C) 2006 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation of the globally defined fields and initialization routines
* for the gauge field
*
* The externally accessible functions are
*
*   void start_ranlux(int level,int seed)
*     Initializes the random number generators ranlxs and ranlxd on all
*     processes in different ways. The luxury level should be 0 (recommended)
*     or 1 (exceptional), and the seed can be any positive integer less than
*     2^31/NPROC
*
*   void alloc_u(void)
*     Allocates the memory space for the global single-precision gauge field,
*     initializes the global pointer pu and sets the link variables to unity
*
*   void alloc_ud(void)
*     Allocates the memory space for the global double-precision gauge field,
*     initializes the global pointer pud and sets the link variables to unity
*
*   void alloc_s(int no_fields)
*     Allocates the memory space for "no_fields" global single-precision
*     spinor fields, initializes the global pointer ps and sets the spinor
*     components to zero
*
*   void alloc_sd(int no_fields)
*     Allocates the memory space for "no_fields" global double-precision
*     spinor fields, initializes the global pointer psd and sets the spinor
*     components to zero
*
*   void free_u(void)
*     Frees the memory space previously allocated for the single-precision
*     gauge field and resets the global pointer pu to NULL. Communication
*     buffers associated to this field are also freed
*
*   void free_ud(void)
*     Frees the memory space previously allocated for the double-precision
*     gauge field and resets the global pointer pud to NULL. Communication
*     buffers associated to this field are also freed
*
*   void free_s(void)
*     Frees the memory space previously allocated for the single-precision
*     spinor fields and resets the global pointer ps to NULL
*
*   void free_sd(void)
*     Frees the memory space previously allocated for the double-precision
*     spinor fields and resets the global pointer psd to NULL
*
*   void random_u(void)
*     Initializes the global single-precision link variables to uniformly
*     distributed random SU(3) matrices
*
*   void random_ud(void)
*     Initializes the global double-precision link variables to uniformly
*     distributed random SU(3) matrices
*
*   void renormalize_u(void)
*     Projects the global single-precision link variables back to SU(3)
*
*   void renormalize_ud(void)
*     Projects the global double-precision link variables back to SU(3)
*
*   void assign_u2ud(void)
*     Assigns the global single-precision link variables to the global
*     double-precision link variables and projects the latter back to SU(3)
*
*   void assign_ud2u(void)
*     Assigns the global double-precision link variables to the global
*     single-precision link variables
*
* Notes:
*
* All these programs act globally and should thus be unconditionally called
* from all processes simultaneously
*
* The programs alloc_u and alloc_ud allocate memory space for all the 8 link
* variables attached to the odd sites of the local lattice (the "local gauge
* field") as well as for the set of link variables at the even sites which
* stick out of the local lattice in positive directions
*
* Note that random_u(), random_ud(), renormalize_u(), renormalize_ud(),
* assign_u2ud() and assign_ud2u() operate on the local gauge field only
*
*******************************************************************************/

#define START_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "misc.h"
#include "start.h"
#include "global.h"

static const su3 u0={{0.0f}};
static const su3_dble ud0={{0.0}};
static const spinor s0={{{0.0f}}};
static const spinor_dble sd0={{{0.0}}};

static su3 *ub=NULL;
static su3_dble *udb=NULL;


void start_ranlux(int level,int seed)
{
   int my_rank,no_proc,max_seed,loc_seed,iprms[2];

   if (NPROC>1)
   {
      iprms[0]=level;
      iprms[1]=seed;

      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);
   
      error((iprms[0]!=level)||(iprms[1]!=seed),1,
            "start_ranlux [start.c]","Input parameters are not global");
   }
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&no_proc);

   max_seed=2147483647/no_proc;
   loc_seed=seed+my_rank*max_seed;

   error_root((level<0)||(level>1)||(seed<1)||(seed>max_seed),1,
              "start_ranlux [start.c]",
              "level should be 0 or 1, and 0<(seed*no_of_processes)<2^31");

   rlxs_init(level,loc_seed);
   rlxd_init(level+1,loc_seed);
}


void alloc_u(void)
{
   int ix,iy,mu;
   int io[4],nu[4];
   su3 unity,*p;

   error(iup[0][0]==0,1,"alloc_u [start.c]",
         "Geometry arrays are not initialized");

   error_root(sizeof(su3)!=(18*sizeof(float)),1,"alloc_u [start.c]",
         "The su3 structures are not properly packed");

   if (ub==NULL)
      ub=amalloc((4*VOLUME+BNDRY/4)*sizeof(su3),ALIGN);
   error(ub==NULL,1,"alloc_u [start.c]",
         "Could not allocate memory space for the gauge field");

   unity=u0;
   unity.c11.re=1.0f;
   unity.c22.re=1.0f;
   unity.c33.re=1.0f;
   p=ub;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         *p=unity;
         pu[ix][mu]=p;
         p+=1;

         iy=idn[ix][mu];
         *p=unity;
         if (iy<VOLUME)
            pu[iy][mu]=p;
         p+=1;
      }
   }

   io[0]=      (BNDRY+FACE0)/2;
   io[1]=io[0]+(FACE0+FACE1)/2;
   io[2]=io[1]+(FACE1+FACE2)/2;
   io[3]=io[2]+(FACE2+FACE3)/2;

   nu[0]=FACE0/2;
   nu[1]=FACE1/2;
   nu[2]=FACE2/2;
   nu[3]=FACE3/2;

   for (mu=0;mu<4;mu++)
   {
      for (iy=0;iy<nu[mu];iy++)
      {
         ix=map[io[mu]+iy];
         ix=idn[ix][mu];
         ix=map[ix-VOLUME];

         *p=unity;
         pu[ix][mu]=p;
         p+=1;
      }
   }
}


void alloc_ud(void)
{
   int ix,iy,mu;
   int io[4],nu[4];
   su3_dble unity,*p;

   error(iup[0][0]==0,1,"alloc_ud [start.c]",
         "Geometry arrays are not initialized");

   error_root(sizeof(su3_dble)!=(18*sizeof(double)),1,"alloc_ud [start.c]",
         "The su3_dble structures are not properly packed");
   
   if (udb==NULL)
      udb=amalloc((4*VOLUME+BNDRY/4)*sizeof(su3_dble),ALIGN);
   error(udb==NULL,1,"alloc_ud [start.c]",
         "Could not allocate memory space for the gauge field");

   unity=ud0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   p=udb;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         *p=unity;
         pud[ix][mu]=p;
         p+=1;

         iy=idn[ix][mu];
         *p=unity;
         if (iy<VOLUME)
            pud[iy][mu]=p;
         p+=1;
      }
   }

   io[0]=      (BNDRY+FACE0)/2;
   io[1]=io[0]+(FACE0+FACE1)/2;
   io[2]=io[1]+(FACE1+FACE2)/2;
   io[3]=io[2]+(FACE2+FACE3)/2;

   nu[0]=FACE0/2;
   nu[1]=FACE1/2;
   nu[2]=FACE2/2;
   nu[3]=FACE3/2;

   for (mu=0;mu<4;mu++)
   {
      for (iy=0;iy<nu[mu];iy++)
      {
         ix=map[io[mu]+iy];
         ix=idn[ix][mu];
         ix=map[ix-VOLUME];

         *p=unity;
         pud[ix][mu]=p;
         p+=1;
      }
   }
}


void free_s(void)
{
   if (no_s!=0)
   {
      afree(ps[0][0]);
      afree(ps);
      ps=NULL;
      no_s=0;
   }
}


void alloc_s(int no_fields)
{
   int n,ix,iprms[1];
   spinor *p;

   if (NPROC>1)
   {
      iprms[0]=no_fields;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=no_fields),1,
            "alloc_s [start.c]","Input parameter is not global");
   }
   
   error_root(sizeof(spinor)!=(24*sizeof(float)),1,"alloc_s [start.c]",
              "The spinor structures are not properly packed");
   free_s();
   
   if (no_fields>0)
   {
      no_s=no_fields;
      ps=amalloc(no_s*sizeof(*ps),3);
      p=amalloc(no_s*NSPIN*sizeof(spinor),ALIGN);

      error((ps==NULL)||(p==NULL),1,"alloc_s [start.c]",
            "Could not allocate memory space for the spinor fields");

      for (n=0;n<no_s;n++)
      {
         for (ix=0;ix<NSPIN;ix++)
         {
            ps[n][ix]=p;
            *p=s0;
            p+=1;
         }
      }
   }
}


void free_sd(void)
{
   if (no_sd!=0)
   {
      afree(psd[0][0]);
      afree(psd);
      psd=NULL;
      no_sd=0;
   }
}


void alloc_sd(int no_fields)
{
   int n,ix,iprms[1];
   spinor_dble *p;

   if (NPROC>1)
   {
      iprms[0]=no_fields;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=no_fields),1,
            "alloc_sd [start.c]","Input parameter is not global");
   }
   
   error_root(sizeof(spinor_dble)!=(24*sizeof(double)),1,"alloc_sd [start.c]",
              "The spinor_dble structures are not properly packed");
   free_sd();
      
   if (no_fields>0)
   {
      no_sd=no_fields;
      psd=amalloc(no_sd*sizeof(*psd),3);
      p=amalloc(no_sd*NSPIN*sizeof(spinor_dble),ALIGN);

      error((psd==NULL)||(p==NULL),1,"alloc_sd [start.c]",
            "Could not allocate memory space for the spinor fields");

      for (n=0;n<no_sd;n++)
      {
         for (ix=0;ix<NSPIN;ix++)
         {
            psd[n][ix]=p;
            *p=sd0;
            p+=1;
         }
      }
   }
}


void free_u(void)
{
   int ix,mu;

   free_ucom_bufs(1,0);

   if (ub!=NULL)
   {
      afree(ub);
      ub=NULL;
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
         pu[ix][mu]=NULL;
   }
}


void free_ud(void)
{
   int ix,mu;

   free_ucom_bufs(0,1);

   if (udb!=NULL)
   {
      afree(udb);
      udb=NULL;
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
         pud[ix][mu]=NULL;
   }
}


void random_u(void)
{
   su3 *p,*pm;

   error(ub==NULL,1,"random_u [start.c]",
         "Attempt to access unallocated memory space");

   p=ub;
   pm=p+4*VOLUME;

   for (;p<pm;p++)
      random_su3(p);
}


void random_ud(void)
{
   su3_dble *p,*pm;

   error(udb==NULL,1,"random_ud [start.c]",
         "Attempt to access unallocated memory space");

   p=udb;
   pm=p+4*VOLUME;

   for (;p<pm;p++)
      random_su3_dble(p);
}


void renormalize_u(void)
{
   su3 *p,*pm;

   error(ub==NULL,1,"renormalize_u [start.c]",
         "Attempt to access unallocated memory space");

   p=ub;
   pm=ub+4*VOLUME;

   for (;p<pm;p++)
      project_to_su3(p);
}


void renormalize_ud(void)
{
   su3_dble *p,*pm;

   error(udb==NULL,1,"renormalize_ud [start.c]",
         "Attempt to access unallocated memory space");

   p=udb;
   pm=udb+4*VOLUME;

   for (;p<pm;p++)
      project_to_su3_dble(p);
}


void assign_u2ud(void)
{
   int i;
   complex *r;
   complex_dble *rd;
   su3 *p,*pm;
   su3_dble *pd;

   error((ub==NULL)||(udb==NULL),1,"assign_u2ud [start.c]",
         "Attempt to access unallocated memory space");

   p=ub;
   pm=p+4*VOLUME;
   pd=udb;

   for (;p<pm;p++)
   {
      r=(complex*)(p);
      rd=(complex_dble*)(pd);

      for (i=0;i<9;i++)
      {
         rd[i].re=(double)(r[i].re);
         rd[i].im=(double)(r[i].im);
      }

      project_to_su3_dble(pd);
      pd+=1;
   }
}


void assign_ud2u(void)
{
   int i;
   complex *r;
   complex_dble *rd;
   su3 *p,*pm;
   su3_dble *pd;

   error((ub==NULL)||(udb==NULL),1,"assign_ud2u [start.c]",
         "Attempt to access unallocated memory space");

   p=ub;
   pm=p+4*VOLUME;
   pd=udb;

   for (;p<pm;p++)
   {
      r=(complex*)(p);
      rd=(complex_dble*)(pd);

      for (i=0;i<9;i++)
      {
         r[i].re=(float)(rd[i].re);
         r[i].im=(float)(rd[i].im);
      }

      pd+=1;
   }
}

