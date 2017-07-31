
/*******************************************************************************
*
* File Pbnd.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic programs for the projector theta to the exterior boundary of a 
* block of lattice points (version for single-precision fields)
*
* The following are arrays of functions indexed by the face number ifc=0,..,7
*
*   void (*assign_s2w[])(int imb[],int vol,spinor *pk,weyl *pl)
*     Applies the projector theta[ifc] to the spinor pk[imb[ix]], 
*     ix=0,..,vol-1, and assigns the result to the weyl spinor pl[ix]
*
*   void (*add_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl)
*     Expands the Weyl spinor pk[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and adds psi to pl[imb[ix]]
*
*   void (*sub_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl)
*     Expands the Weyl spinor pk[ix], ix=0,..,vol-1, to a Dirac spinor
*     psi satisfying theta[ifc]*psi=psi and subtracts psi from pl[imb[ix]]
*
*   void (*mulg5_sub_assign_w2s[])(int imb[],int vol,weyl *pk,spinor *pl)
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

#define PBND_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "su3.h"
#include "start.h"
#include "global.h"

#if (defined SSE)
#include "sse.h"

static const sse_float poh={0.5f,0.5f,0.5f,0.5f};


static void assign_s2w0(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c3,(*rk).c4);               

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      

      _sse_vector_sub();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w1(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c3,(*rk).c4);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_add();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w2(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c4,(*rk).c3);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_i_sub();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w3(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c4,(*rk).c3);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_i_add();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w4(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c4,(*rk).c3);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_subadd();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);      
   }   
}


static void assign_s2w5(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c4,(*rk).c3);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_addsub();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w6(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c3,(*rk).c4);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_i_subadd();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void assign_s2w7(int imb[],int vol,spinor *pk,weyl *pl)
{
   weyl *plm;
   spinor *rk,*rkn;

   plm=pl+vol;   
   rk=pk+(*imb);
   imb+=(pl<(plm-1));
   rkn=pk+(*imb);

   for (;pl<plm;pl++)
   {
      _sse_pair_load((*rk).c1,(*rk).c2);
      _sse_pair_load_up((*rk).c3,(*rk).c4);      

      rk=rkn;
      imb+=(pl<(plm-2));
      rkn=pk+(*imb);
      _prefetch_spinor(rkn);      
      
      _sse_vector_i_addsub();
      _sse_vector_mul(poh);
      _sse_pair_store((*pl).c1,(*pl).c2);
   }   
}


static void add_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void add_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_add();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void add_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void add_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_sub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void add_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void add_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void add_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void add_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_add();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void sub_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      

      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;

      _sse_vector_add();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void sub_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      

      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;

      _sse_vector_sub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void sub_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;

      _sse_vector_i_sub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void sub_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void sub_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void sub_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void sub_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void sub_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void mulg5_sub_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {   
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      

      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void mulg5_sub_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_add();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void mulg5_sub_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_add();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void mulg5_sub_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_sub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void mulg5_sub_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_subadd();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void mulg5_sub_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c4,(*rl).c3);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_addsub();
      _sse_pair_store((*rl).c4,(*rl).c3);
   }
}


static void mulg5_sub_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_addsub();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}


static void mulg5_sub_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl,*rln,*rlm;

   pkm=pk+vol;
   rln=pl+(*imb);
   imb+=(pk<(pkm-1));
   rlm=pl+(*imb);

   for (;pk<pkm;)
   {      
      _sse_pair_load_up((*pk).c1,(*pk).c2);
      _sse_pair_load((*rln).c1,(*rln).c2);      

      rl=rln;
      rln=rlm;
      imb+=(pk<(pkm-2));
      rlm=pl+(*imb);
      _prefetch_spinor(rlm);      
      
      _sse_vector_sub();
      _sse_pair_store((*rl).c1,(*rl).c2);
      _sse_pair_load((*rl).c3,(*rl).c4);      

      pk+=4;
      _prefetch_weyl(pk);
      pk-=3;
      
      _sse_vector_i_subadd();
      _sse_pair_store((*rl).c3,(*rl).c4);
   }
}

#else

static void assign_s2w0(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w1(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w2(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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

static void assign_s2w3(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w4(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w5(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w6(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void assign_s2w7(int imb[],int vol,spinor *pk,weyl *pl)
{
   float r1;
   weyl *plm;
   spinor *rk;

   r1=0.5f;
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


static void add_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void add_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void sub_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s0(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s1(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s2(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s3(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s4(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s5(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s6(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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


static void mulg5_sub_assign_w2s7(int imb[],int vol,weyl *pk,spinor *pl)
{
   weyl *pkm;
   spinor *rl;

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

void (*assign_s2w[8])(int imb[],int vol,spinor *pk,weyl *pl) =
{assign_s2w0,assign_s2w1,assign_s2w2,assign_s2w3,
 assign_s2w4,assign_s2w5,assign_s2w6,assign_s2w7};

void (*add_assign_w2s[8])(int imb[],int vol,weyl *pk,spinor *pl) =
{add_assign_w2s0,add_assign_w2s1,add_assign_w2s2,add_assign_w2s3,
 add_assign_w2s4,add_assign_w2s5,add_assign_w2s6,add_assign_w2s7};

void (*sub_assign_w2s[8])(int imb[],int vol,weyl *pk,spinor *pl) =
{sub_assign_w2s0,sub_assign_w2s1,sub_assign_w2s2,sub_assign_w2s3,
 sub_assign_w2s4,sub_assign_w2s5,sub_assign_w2s6,sub_assign_w2s7};

void (*mulg5_sub_assign_w2s[8])(int imb[],int vol,weyl *pk,spinor *pl) =
{mulg5_sub_assign_w2s0,mulg5_sub_assign_w2s1,
 mulg5_sub_assign_w2s2,mulg5_sub_assign_w2s3,
 mulg5_sub_assign_w2s4,mulg5_sub_assign_w2s5,
 mulg5_sub_assign_w2s6,mulg5_sub_assign_w2s7};

