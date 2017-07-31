
/*******************************************************************************
*
* File Pbnd_dble.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic programs for the projector theta to the exterior boundary of a 
* block of lattice points (version for double-precision fields)
*
* The following are arrays of functions indexed by the face number ifc=0,..,7
*
*   void (*assign_sd2wd[])(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
*     Applies the projector theta[ifc] to the spinor pk[imb[ix]], 
*     ix=0,..,vol-1, and assigns the result to the weyl spinor pl[ix]
*
*   void (*add_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
*     Expands the Weyl spinor pk[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and adds psi to pl[imb[ix]]
*
*   void (*sub_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
*     Expands the Weyl spinor pk[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts psi from pl[imb[ix]]
*
*   void (*mulg5_sub_assign_wd2sd[])(int imb[],int vol,weyl_dble *pk,
*                                    spinor_dble *pl)
*     Expands the Weyl spinor pk[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts gamma5*psi from
*     pl[imb[ix]]
*
* Notes:
*
* The projector theta was introduced in section 3.3 of 
*
*   Lattice QCD and the Schwarz alternating procedure, JHEP 0305 (2003) 052
*
* Block faces in the -0,+0,..,-3,+3 directions are labelled by an index
* ifc=0,..,7 and the projector on a given face is then given by
*
*   theta[ifc] = (1/2)*(1+gamma_mu) if ifc=2*mu,
*
*              = (1/2)*(1-gamma_mu) if ifc=2*mu+1
*  
* Dirac fields on the face that satisfy theta[ifc]*psi=psi are completely
* characterized by their first two Dirac components. In this way they are
* mapped to Weyl fields on the face in an invertible manner
*
* The size and position of the faces is only implicitly defined through
* the parameter vol and the array imb[] of the indices of the points on
* the face
*
* None of these programs involves communications. They are general purpose
* routines that know nothing about the underlying geometry. In particular,
* they can be called locally
*
*******************************************************************************/

#define PBND_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "su3.h"
#include "start.h"
#include "global.h"

#if (defined SSE2)
#include "sse2.h"

static const sse_double poh={0.5,0.5};


static void assign_sd2wd0(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c3);            
   
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c4);            

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      

      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);            
   }
}


static void assign_sd2wd1(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c3);            

      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c4);            

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void assign_sd2wd2(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c4);

      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c3);

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2); 
   }
}


static void assign_sd2wd3(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c4);

      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c3);

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void assign_sd2wd4(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c4);            

      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c3);            

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void assign_sd2wd5(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c4);            

      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c3);            

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void assign_sd2wd6(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c3);

      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c4);

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void assign_sd2wd7(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *plm;
   spinor_dble *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {   
      _sse_load_dble((*rk).c1);
      _sse_load_up_dble((*rk).c3);

      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c1);

      _sse_load_dble((*rk).c2);
      _sse_load_up_dble((*rk).c4);

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor_dble(rkn);      
      
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_vector_mul_dble(poh);
      _sse_store_dble((*pl).c2);
   }
}


static void add_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
   }
}


static void add_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void add_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void add_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3); 
   }
}


static void add_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void add_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void add_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void add_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void sub_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      

      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
      
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;

      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
   }
}


static void sub_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);   
   }
}


static void sub_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3); 
   }
}


static void sub_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void sub_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3); 
   }
}


static void sub_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void sub_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void sub_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void mulg5_sub_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
   }
}


static void mulg5_sub_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4); 
   }
}


static void mulg5_sub_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3); 
   }
}


static void mulg5_sub_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void mulg5_sub_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);  
   }
}


static void mulg5_sub_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c4);            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c3);            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);
   }
}


static void mulg5_sub_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c4);
   }
}


static void mulg5_sub_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {         
      _sse_load_up_dble((*pk).c1);
      _sse_load_dble((*rln).c1);            

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor_dble(rlm);      
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c1);
      _sse_load_dble((*rl).c3);
      _sse_vector_i_mul_dble();            
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c3);            
            
      _sse_load_up_dble((*pk).c2);
      _sse_load_dble((*rl).c2);            

      pk+=4;
      _prefetch_weyl_dble(pk);
      pk-=3;
      
      _sse_vector_sub_dble();
      _sse_store_dble((*rl).c2);
      _sse_load_dble((*rl).c4);
      _sse_vector_i_mul_dble();
      _sse_vector_add_dble();
      _sse_store_dble((*rl).c4);
   }
}

#else

static void assign_sd2wd0(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_sub((*pl).c1,(*rk).c1,(*rk).c3);
      _vector_sub((*pl).c2,(*rk).c2,(*rk).c4);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);      
   }
}


static void assign_sd2wd1(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_add((*pl).c1,(*rk).c1,(*rk).c3);
      _vector_add((*pl).c2,(*rk).c2,(*rk).c4);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2); 
   }
}


static void assign_sd2wd2(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_i_sub((*pl).c1,(*rk).c1,(*rk).c4);
      _vector_i_sub((*pl).c2,(*rk).c2,(*rk).c3);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);     
   }
}


static void assign_sd2wd3(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_i_add((*pl).c1,(*rk).c1,(*rk).c4);
      _vector_i_add((*pl).c2,(*rk).c2,(*rk).c3);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);    
   }
}


static void assign_sd2wd4(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_sub((*pl).c1,(*rk).c1,(*rk).c4);
      _vector_add((*pl).c2,(*rk).c2,(*rk).c3);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);
   }
}


static void assign_sd2wd5(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_add((*pl).c1,(*rk).c1,(*rk).c4);
      _vector_sub((*pl).c2,(*rk).c2,(*rk).c3);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2); 
   }
}


static void assign_sd2wd6(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_i_sub((*pl).c1,(*rk).c1,(*rk).c3);
      _vector_i_add((*pl).c2,(*rk).c2,(*rk).c4);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);  
   }
}


static void assign_sd2wd7(int imb[],int vol,spinor_dble *pk,weyl_dble *pl)
{
   double r1;
   weyl_dble *plm;
   spinor_dble *rk;

   r1=0.5;
   plm=pl+vol;

   for (;pl<plm;pl++)
   {
      rk=pk+(*imb);
      imb+=1;
      _vector_i_add((*pl).c1,(*rk).c1,(*rk).c3);
      _vector_i_sub((*pl).c2,(*rk).c2,(*rk).c4);       
      _vector_mul((*pl).c1,r1,(*pl).c1);
      _vector_mul((*pl).c2,r1,(*pl).c2);  
   }
}


static void add_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c1);
      _vector_sub_assign((*rl).c4,(*pk).c2);
   }
}


static void add_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c1);
      _vector_add_assign((*rl).c4,(*pk).c2); 
   }
}

   
static void add_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c2);
      _vector_i_add_assign((*rl).c4,(*pk).c1);  
   }
}

   
static void add_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c2);
      _vector_i_sub_assign((*rl).c4,(*pk).c1);
   }
}

   
static void add_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c2);
      _vector_sub_assign((*rl).c4,(*pk).c1);
   }
}

   
static void add_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c2);
      _vector_add_assign((*rl).c4,(*pk).c1);
   }
}

   
static void add_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c1);
      _vector_i_sub_assign((*rl).c4,(*pk).c2);
   }
}

   
static void add_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_add_assign((*rl).c1,(*pk).c1);
      _vector_add_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c1);
      _vector_i_add_assign((*rl).c4,(*pk).c2);
   }
}


static void sub_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c1);
      _vector_add_assign((*rl).c4,(*pk).c2);
   }
}


static void sub_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c1);
      _vector_sub_assign((*rl).c4,(*pk).c2); 
   }
}


static void sub_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c2);
      _vector_i_sub_assign((*rl).c4,(*pk).c1);  
   }
}


static void sub_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c2);
      _vector_i_add_assign((*rl).c4,(*pk).c1);
   }
}


static void sub_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c2);
      _vector_add_assign((*rl).c4,(*pk).c1);
   }
}


static void sub_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c2);
      _vector_sub_assign((*rl).c4,(*pk).c1);
   }
}


static void sub_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c1);
      _vector_i_add_assign((*rl).c4,(*pk).c2);
   }
}


static void sub_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c1);
      _vector_i_sub_assign((*rl).c4,(*pk).c2);
   }
}


static void mulg5_sub_assign_wd2sd0(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c1);
      _vector_sub_assign((*rl).c4,(*pk).c2);
   }
}


static void mulg5_sub_assign_wd2sd1(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c1);
      _vector_add_assign((*rl).c4,(*pk).c2); 
   }
}


static void mulg5_sub_assign_wd2sd2(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c2);
      _vector_i_add_assign((*rl).c4,(*pk).c1); 
   }
}


static void mulg5_sub_assign_wd2sd3(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c2);
      _vector_i_sub_assign((*rl).c4,(*pk).c1);
   }
}


static void mulg5_sub_assign_wd2sd4(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_add_assign((*rl).c3,(*pk).c2);
      _vector_sub_assign((*rl).c4,(*pk).c1);
   }
}


static void mulg5_sub_assign_wd2sd5(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_sub_assign((*rl).c3,(*pk).c2);
      _vector_add_assign((*rl).c4,(*pk).c1);
   }
}


static void mulg5_sub_assign_wd2sd6(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_add_assign((*rl).c3,(*pk).c1);
      _vector_i_sub_assign((*rl).c4,(*pk).c2);
   }
}


static void mulg5_sub_assign_wd2sd7(int imb[],int vol,weyl_dble *pk,
                                    spinor_dble *pl)
{
   weyl_dble *pkm;
   spinor_dble *rl;

   pkm=pk+vol;
   
   for (;pk<pkm;pk++)
   {
      rl=pl+(*imb);
      imb+=1;
      _vector_sub_assign((*rl).c1,(*pk).c1);
      _vector_sub_assign((*rl).c2,(*pk).c2);
      _vector_i_sub_assign((*rl).c3,(*pk).c1);
      _vector_i_add_assign((*rl).c4,(*pk).c2);
   }
}

#endif

void (*assign_sd2wd[8])(int imb[],int vol,spinor_dble *pk,weyl_dble *pl) =
{assign_sd2wd0,assign_sd2wd1,assign_sd2wd2,assign_sd2wd3,
 assign_sd2wd4,assign_sd2wd5,assign_sd2wd6,assign_sd2wd7};

void (*add_assign_wd2sd[8])(int imb[],int vol,weyl_dble *pk,spinor_dble *pl) =
{add_assign_wd2sd0,add_assign_wd2sd1,add_assign_wd2sd2,add_assign_wd2sd3,
 add_assign_wd2sd4,add_assign_wd2sd5,add_assign_wd2sd6,add_assign_wd2sd7};

void (*sub_assign_wd2sd[8])(int imb[],int vol,weyl_dble *pk,spinor_dble *pl) =
{sub_assign_wd2sd0,sub_assign_wd2sd1,sub_assign_wd2sd2,sub_assign_wd2sd3,
 sub_assign_wd2sd4,sub_assign_wd2sd5,sub_assign_wd2sd6,sub_assign_wd2sd7};

void (*mulg5_sub_assign_wd2sd[8])(int imb[],int vol,weyl_dble *pk,
                                  spinor_dble *pl) =
{mulg5_sub_assign_wd2sd0,mulg5_sub_assign_wd2sd1,
 mulg5_sub_assign_wd2sd2,mulg5_sub_assign_wd2sd3,
 mulg5_sub_assign_wd2sd4,mulg5_sub_assign_wd2sd5,
 mulg5_sub_assign_wd2sd6,mulg5_sub_assign_wd2sd7};
